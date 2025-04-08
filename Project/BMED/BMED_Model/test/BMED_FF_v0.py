# Load Modules
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from sklearn.metrics import r2_score
from sklearn.preprocessing import StandardScaler
import os
import matplotlib.pyplot as plt
import pickle

class BMEDDataset(Dataset):
    '''BMED 데이터로부터 PyTorch 데이터셋 생성 - 피드포워드 네트워크용'''
    def __init__(self, data_path, fold_idx=0, n_folds=5, train=True, scalers=None):
        self.train = train
        self.fold_idx = fold_idx
        self.n_folds = n_folds

        # Spline 데이터 로드
        self.df = pd.read_excel(data_path, sheet_name='spline_data')
        
        # 스케일러 초기화 또는 전달받은 스케일러 사용
        if scalers is None:
            self.feature_scaler = StandardScaler()
            self.mol_change_scaler = StandardScaler()
            self.state_scaler = StandardScaler()
        else:
            self.feature_scaler = scalers['feature']
            self.mol_change_scaler = scalers['mol_change']
            self.state_scaler = scalers['state']

        # 실험 번호별로 데이터 분리
        self.dic_spline = {}
        for exp in self.df['exp'].unique():
            self.dic_spline[exp] = self.df[self.df['exp'] == exp].sort_values('t')
        
        # k-fold를 위한 실험 분할
        all_exps = list(self.dic_spline.keys())
        np.random.seed(42)
        np.random.shuffle(all_exps)
        np.random.seed(None)
        
        # 실험을 n_folds 개의 그룹으로 나누기
        fold_size = len(all_exps) // n_folds
        test_start = fold_idx * fold_size
        test_end = (fold_idx + 1) * fold_size if fold_idx < n_folds - 1 else len(all_exps)
        
        if train:
            # 현재 fold를 제외한 모든 실험을 훈련 데이터로 사용
            self.exps_to_use = all_exps[:test_start] + all_exps[test_end:]
        else:
            # 현재 fold의 실험을 테스트 데이터로 사용
            self.exps_to_use = all_exps[test_start:test_end]
            
        # 데이터 준비
        self.prepare_data()

        # 실험별 인덱스 저장 (배치 처리용)
        self.exp_indices = []
        cur_idx = 0
        for exp_id in self.exps_to_use:
            exp_data = self.dic_spline[exp_id]
            # 각 데이터 포인트에 대해 하나의 샘플 생성 (시퀀스 길이 필요 없음)
            n_samples = len(exp_data) - 1  # 각 데이터 포인트에 대해 다음 상태가 필요하므로 -1
            indices = list(range(cur_idx, cur_idx + n_samples))
            self.exp_indices.append((exp_id, indices))
            cur_idx += n_samples

    def get_scalers(self):
        """스케일러 반환"""
        return {
            'feature': self.feature_scaler,
            'mol_change': self.mol_change_scaler,
            'state': self.state_scaler
        }

    def prepare_data(self):
        features_list = []
        mol_changes_targets_list = []
        states_targets_list = []

        for exp in self.exps_to_use:
            exp_data = self.dic_spline[exp]

            # 각 시간 단계에 대해 하나의 샘플 생성
            for i in range(len(exp_data) - 1):  # 마지막 데이터 포인트는 다음 상태가 없으므로 제외
                cur = exp_data.iloc[i]
                next = exp_data.iloc[i+1]
                
                # 특성 벡터 (현재 상태)
                feature = [
                    cur['T'], cur['V'], cur['E'],
                    cur['CF_LA'], cur['CF_K'], cur['CA_LA'], cur['CB_K'],
                    cur['VF'], cur['VA'], cur['VB']
                ]
                
                # 시간 간격
                dt = next['t'] - cur['t']

                # 몰수 계산 (mol)
                cur_NA_LA = cur['CA_LA'] * cur['VA']
                cur_NB_K = cur['CB_K'] * cur['VB']
                next_NA_LA = next['CA_LA'] * next['VA']
                next_NB_K = next['CB_K'] * next['VB']

                # 몰 변화량 (mol/hr)
                dNLA = (next_NA_LA - cur_NA_LA) / dt
                dNK = (next_NB_K - cur_NB_K) / dt

                # 물 변화량 (L/hr)
                dVA = (next['VA'] - cur['VA']) / dt
                dVB = (next['VB'] - cur['VB']) / dt

                # 변화량 벡터
                mol_change = [dNLA, dNK, dVA, dVB]
                
                # 다음 상태 벡터
                state = [
                    next['T'], next['V'], next['E'],
                    next['CF_LA'], next['CF_K'], next['CA_LA'], next['CB_K'],
                    next['VF'], next['VA'], next['VB']
                ]
                
                features_list.append(feature)
                mol_changes_targets_list.append(mol_change)
                states_targets_list.append(state)
                
        # 데이터 배열로 변환
        self.features_raw = np.array(features_list, dtype=np.float32)
        self.mol_changes_targets_raw = np.array(mol_changes_targets_list, dtype=np.float32)
        self.states_targets_raw = np.array(states_targets_list, dtype=np.float32)
        
        # 데이터 정규화
        if self.train:  # 학습 데이터일 경우에만 스케일러 학습
            self.feature_scaler.fit(self.features_raw)
            self.mol_change_scaler.fit(self.mol_changes_targets_raw)
            self.state_scaler.fit(self.states_targets_raw)
        
        # 정규화 적용
        self.features = self.feature_scaler.transform(self.features_raw)
        self.mol_changes_targets = self.mol_change_scaler.transform(self.mol_changes_targets_raw)
        self.states_targets = self.state_scaler.transform(self.states_targets_raw)

        print(f"Data loaded: {len(self.features)} samples, feature shape: {self.features.shape}")
        print(f"Mol change targets shape: {self.mol_changes_targets.shape}")
        print(f"State targets shape: {self.states_targets.shape}")

    def __len__(self):
        return len(self.features)
                        
    def __getitem__(self, idx):
        # 3D 텐서가 아니라 2D 텐서를 반환 (배치 차원, 특성 차원만 유지)
        return (
            torch.FloatTensor(self.features[idx]).unsqueeze(0),  # [1, 10] 형태로 변환 (유연성을 위해)
            torch.FloatTensor(self.mol_changes_targets[idx]),
            torch.FloatTensor(self.states_targets[idx]).unsqueeze(0)  # [1, 10] 형태로 변환
        )
    
    def inverse_transform_predictions(self, mol_pred, state_pred=None):
        """예측값을 원래 스케일로 변환"""
        if isinstance(mol_pred, torch.Tensor):
            mol_pred = mol_pred.detach().cpu().numpy()
        
        # 예측된 몰 변화량을 원래 스케일로 변환
        if len(mol_pred.shape) == 2:  # [batch_size, 4] 형태
            mol_orig = self.mol_change_scaler.inverse_transform(mol_pred)
        else:  # 다른 형태일 경우
            # 2D 형태로 변환
            orig_shape = mol_pred.shape
            mol_pred_2d = mol_pred.reshape(-1, mol_pred.shape[-1])
            
            # 원래 스케일로 변환
            mol_orig_2d = self.mol_change_scaler.inverse_transform(mol_pred_2d)
            
            # 원래 형태로 복원
            mol_orig = mol_orig_2d.reshape(orig_shape)
        
        # state 예측값이 있는 경우에도 변환
        if state_pred is not None:
            if isinstance(state_pred, torch.Tensor):
                state_pred = state_pred.detach().cpu().numpy()
            
            # 예측된 상태를 원래 스케일로 변환
            if len(state_pred.shape) == 2:  # [batch_size, 10] 형태
                state_orig = self.state_scaler.inverse_transform(state_pred)
            else:  # 다른 형태일 경우
                orig_shape = state_pred.shape
                state_pred_2d = state_pred.reshape(-1, state_pred.shape[-1])
                state_orig_2d = self.state_scaler.inverse_transform(state_pred_2d)
                state_orig = state_orig_2d.reshape(orig_shape)
            
            return mol_orig, state_orig
        
        return mol_orig

class PhysicsLayer(nn.Module):
    '''
    물리적 법칙을 반영하는 레이어
    몰 변화량 예측값을 받아 물질수지와 부피 변화를 계산
    '''
    def __init__(self, time_step=0.1):
        super(PhysicsLayer, self).__init__()
        self.time_step = time_step # unit: hour

    def forward(self, mol_changes, features):
        '''현재 시간 단계 상태 추출'''
        T = features[..., 0:1] # 현재 상태 온도 [C]
        V = features[..., 1:2] # 현재 상태 전압 [V]
        E = features[..., 2:3] # 현재 상태 전해질 농도 [mol/L]
        CFLA = features[..., 3:4] # 현재 상태 Feed LA 농도 [mol/L]
        CFK = features[..., 4:5] # 현재 상태 Feed K+ 농도 [mol/L]
        CALA = features[..., 5:6] # 현재 상태 Acid LA 농도 [mol/L]
        CBK = features[..., 6:7] # 현재 상태 Base K+ 농도 [mol/L]
        VF = features[..., 7:8] # 현재 상태 Feed 부피 [L]
        VA = features[..., 8:9] # 현재 상태 Acid 부피 [L]
        VB = features[..., 9:10] # 현재 상태 Base 부피 [L]
        
        '''현재 시간 단계 몰 계산'''
        NFLA = CFLA * VF # 현재 상태 Feed LA 몰수 [mol]
        NFK = CFK * VF # 현재 상태 Feed K+ 몰수 [mol]
        NALA = CALA * VA # 현재 상태 Acid LA 몰수 [mol]
        NBK = CBK * VB # 현재 상태 Base K+ 몰수 [mol]

        '''변화량 상태 추출'''
        JLA = mol_changes[..., 0:1] # 시간당 LA 몰 변화량 [mol/h]
        JK = mol_changes[..., 1:2] # 시간당 K+ 몰 변화량 [mol/h]
        JVA = mol_changes[..., 2:3] # 시간당 물 부피 변화량 (Acid) [L/h]
        JVB = mol_changes[..., 3:4] # 시간당 물 부피 변화량 (Base) [L/h]

        '''time step에서의 변화량 계산'''
        dLA = JLA * self.time_step # LA 몰 변화량 [mol]
        dK = JK * self.time_step # K+ 몰 변화량[mol]
        dVA = JVA * self.time_step # 물 부피 변화량 (Acid) [L]
        dVB = JVB * self.time_step # 물 부피 변화량 (Base) [L]

        '''부피 업데이트'''
        nVF = VF - dVA - dVB # 다음 상태 Feed 부피 [L]
        nVA = VA + dVA # 다음 상태 Acid 부피 [L]
        nVB = VB + dVB # 다음 상태 Base 부피 [L]

        '''몰수 업데이트'''
        nNFLA = NFLA - dLA # 다음 상태 Feed LA 몰수 [mol]
        nNFK = NFK - dK # 다음 상태 Feed K+ 몰수 [mol]
        nNALA = NALA + dLA # 다음 상태 Acid LA 몰수 [mol]
        nNBK = NBK + dK # 다음 상태 Base K+ 몰수 [mol]
        
        '''농도 업데이트'''
        eps = 1e-6 # 0으로 나누기 방지
        nCFLA = nNFLA / (nVF + eps) # 다음 상태 Feed LA 농도 [mol/L]
        nCFK = nNFK / (nVF + eps) # 다음 상태 Feed K+ 농도 [mol/L]
        nCALA = nNALA / (nVA + eps) # 다음 상태 Acid LA 농도 [mol/L]
        nCBK = nNBK / (nVB + eps) # 다음 상태 Base K+ 농도 [mol/L]

        '''결과 출력'''
        new_states = torch.cat([
            T, V, E, nCFLA, nCFK, nCALA, nCBK, nVF, nVA, nVB
        ], dim=-1)

        return new_states

class MembraneSystemModel(nn.Module):
    '''멤브레인 시스템 모델링을 위한 Physics-Informed 피드포워드 모델'''
    def __init__(self, hidden_units=64, hidden_layers=3, time_step=0.1):
        super(MembraneSystemModel, self).__init__()
        self.hidden_units = hidden_units
        self.hidden_layers = hidden_layers
        self.time_step = time_step

        # 피드포워드 네트워크 구성
        layers = []
        # 입력 레이어
        layers.append(nn.Linear(10, hidden_units))  # [T, V, E, CF_LA, CF_K, CA_LA, CB_K, VF, VA, VB]
        layers.append(nn.ReLU())
        
        # 은닉 레이어
        for _ in range(hidden_layers - 1):
            layers.append(nn.Linear(hidden_units, hidden_units))
            layers.append(nn.ReLU())
        
        # 출력 레이어
        layers.append(nn.Linear(hidden_units, 4))  # 4개의 flux 값 출력
        
        # 시퀀셜 모듈로 레이어 구성
        self.flux_predictor = nn.Sequential(*layers)
        
        # Physics layer 유지
        self.physics_layer = PhysicsLayer(time_step=self.time_step)

    def forward(self, x, _=None):  # hidden 상태는 이제 사용하지 않음
        '''특성 추출'''
        # 시퀀스 차원 제거 (배치 처리를 위해)
        x_flat = x.squeeze(1)  # [batch_size, 10]
        
        '''몰 변화량 예측'''
        mol_changes = self.flux_predictor(x_flat)  # [batch_size, 4]
        
        '''물리 법칙 적용'''
        new_states = self.physics_layer(mol_changes, x_flat)  # [batch_size, 10]
        
        # API 호환성을 위해 None 반환 (기존 코드에서 hidden 상태가 사용되던 자리)
        return mol_changes, new_states, None

class MembraneSystemTrainer:
    """
    멤브레인 시스템 모델 훈련 및 시뮬레이션을 위한 클래스
    """
    def __init__(self, model, device='cuda' if torch.cuda.is_available() else 'cpu', epochs=100):
        self.model = model.to(device)
        self.device = device
        self.epochs = epochs
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=0.001)
        self.criterion = nn.MSELoss()
    
    def train(self, train_loader, val_loader=None, patience=100):
        """
        모델 훈련 - 피드포워드 네트워크를 위해 수정
        """
        self.model.train()
        train_losses = []
        val_losses = []
        
        # Early stopping 변수 초기화
        best_val_loss = float('inf')
        best_model_state = None
        counter = 0
        
        for epoch in range(self.epochs):
            epoch_loss = 0.0
            
            # 전체 배치에 대해 훈련 (실험별 그룹화 제거, 단순화)
            for batch_idx, (features, mol_change_targets, state_targets) in enumerate(train_loader):
                # 데이터를 디바이스로 이동
                features = features.to(self.device)
                state_targets = state_targets.to(self.device)
                
                # 그래디언트 초기화
                self.optimizer.zero_grad()
                
                # 순전파 (hidden 상태 미사용)
                mol_changes, state_predictions, _ = self.model(features)
                
                # 손실 계산
                loss = self.criterion(state_predictions, state_targets.squeeze(1))
                
                # 역전파
                loss.backward()
                
                # 그래디언트 클리핑
                torch.nn.utils.clip_grad_norm_(self.model.parameters(), max_norm=1.0)
                
                # 최적화 스텝
                self.optimizer.step()
                
                # 손실 누적
                epoch_loss += loss.item()
            
            # 에폭 평균 손실
            avg_loss = epoch_loss / len(train_loader)
            train_losses.append(avg_loss)

            # 검증 부분
            if val_loader is not None:
                self.model.eval()
                val_loss = 0.0
                
                with torch.no_grad():
                    for features, _, state_targets in val_loader:
                        features = features.to(self.device)
                        state_targets = state_targets.to(self.device)
                        
                        # 순전파 (hidden 상태 미사용)
                        _, state_predictions, _ = self.model(features)
                        
                        # 손실 계산
                        loss = self.criterion(state_predictions, state_targets.squeeze(1))
                        
                        # 손실 누적
                        val_loss += loss.item()
                
                # 검증 평균 손실
                avg_val_loss = val_loss / len(val_loader)
                val_losses.append(avg_val_loss)
                
                # Epoch, Train Loss, Val Loss 출력
                print(f'Epoch {epoch+1}/{self.epochs}, Train Loss: {avg_loss:.6f}, Val Loss: {avg_val_loss:.6f}')
                
                # Early stopping 체크
                if avg_val_loss < best_val_loss:
                    best_val_loss = avg_val_loss
                    best_model_state = self.model.state_dict().copy()
                    counter = 0
                else:
                    counter += 1
                    print(f'EarlyStopping counter: {counter} out of {patience}')
                    if counter >= patience:
                        print(f'Early stopping: 검증 손실이 {patience}번 연속으로 향상되지 않았습니다.')
                        # 최적의 모델 가중치로 복원
                        self.model.load_state_dict(best_model_state)
                        return train_losses, val_losses
                
                self.model.train()
            else:
                # 검증을 수행하지 않는 경우 출력
                print(f'Epoch {epoch+1}/{self.epochs}, Train Loss: {avg_loss:.6f}')
        
        # 학습이 완료되었을 때 최적의 모델 가중치로 복원
        if best_model_state is not None:
            self.model.load_state_dict(best_model_state)
            
        return train_losses, val_losses
    
    def evaluate_model(self, test_loader, test_dataset):
        self.model.eval()
        test_loss = 0
        test_predictions = []
        test_actuals = []

        with torch.no_grad():
            for features, mol_change_targets, state_targets in test_loader:
                features = features.to(self.device)
                mol_change_targets = mol_change_targets.to(self.device)
                state_targets = state_targets.to(self.device)

                # 모델 예측 (hidden 상태 추가)
                mol_change_predictions, state_predictions, _ = self.model(features)

                # 손실 계산
                loss = self.criterion(state_predictions, state_targets.squeeze(1))
                test_loss += loss.item()

                # 예측값과 실제값 저장
                test_predictions.append(state_predictions.cpu().numpy())
                test_actuals.append(state_targets.cpu().numpy())
         
        # 전체 테스트 세트에 대한 평균 손실
        avg_test_loss = test_loss / len(test_loader)                
                
        # 예측값과 실제값을 하나의 배열로 결합
        test_predictions = np.concatenate(test_predictions, axis=0)
        test_actuals = np.concatenate(test_actuals, axis=0)

        # R² 점수 계산 및 Parity Plot 생성
        r2_scores = calculate_test_r2_and_create_parity_plots(test_dataset, test_predictions, test_actuals)
        
        # 각 변수별 R² 점수 계산 - 차원 체크 및 수정
        r2_per_var = {}
        for i, var_name in enumerate(['T', 'V', 'E', 'CF_LA', 'CF_K', 'CA_LA', 'CB_K', 'VF', 'VA', 'VB']):
            # 차원 체크
            if test_actuals.ndim == 3:
                actuals_flatten = test_actuals[:, :, i].flatten()
            else:  # 2차원인 경우
                actuals_flatten = test_actuals[:, i].flatten()
            
            if test_predictions.ndim == 3:
                preds_flatten = test_predictions[:, :, i].flatten()
            else:  # 2차원인 경우
                preds_flatten = test_predictions[:, i].flatten()
            
            r2 = r2_score(actuals_flatten, preds_flatten)
            r2_per_var[var_name] = r2
        
        return {
            'test_loss': avg_test_loss,
            'predictions': test_predictions,
            'actuals': test_actuals,
            'r2_scores': r2_scores,
            'r2_per_variable': r2_per_var
        }

class ExperimentBatchSampler(torch.utils.data.Sampler):
    """
    실험 단위로 배치를 구성하는 샘플러
    각 실험은 연속적으로 처리되며 시계열 순서가 유지됨
    """
    def __init__(self, dataset, batch_size):
        self.dataset = dataset
        self.batch_size = batch_size
        self.exp_indices = dataset.exp_indices
        self.current_batch = []
        
    def __iter__(self):
        # 실험 ID 리스트 생성 및 섞기 (실험 순서는 랜덤)
        exp_ids = [exp_id for exp_id, _ in self.exp_indices]
        np.random.shuffle(exp_ids)
        
        # 각 실험마다
        for exp_id in exp_ids:
            # 해당 실험의 인덱스 가져오기
            exp_indices = None
            for e_id, indices in self.exp_indices:
                if e_id == exp_id:
                    exp_indices = indices
                    break
            
            if exp_indices is None or len(exp_indices) == 0:
                continue
                
            # 이 실험 내 인덱스들을 배치 크기로 분할
            num_batches = (len(exp_indices) + self.batch_size - 1) // self.batch_size
            
            # 실험 내 배치들을 순서대로 생성
            for i in range(num_batches):
                start_idx = i * self.batch_size
                end_idx = min((i + 1) * self.batch_size, len(exp_indices))
                self.current_batch = exp_indices[start_idx:end_idx]
                
                # 배치가 비어있지 않은 경우에만 반환
                if self.current_batch:
                    yield self.current_batch
    
    def __len__(self):
        # 대략적인 배치 수 계산
        total_batches = 0
        for _, indices in self.exp_indices:
            total_batches += (len(indices) + self.batch_size - 1) // self.batch_size
        return total_batches

    def get_current_batch_indices(self):
        """현재 배치의 인덱스 목록 반환"""
        return self.current_batch

def data_loaders_from_csv(data_path, batch_size=32, train_ratio=0.8):
    """
    CSV 파일에서 훈련 및 검증 데이터 로더 생성 (피드포워드 네트워크용)
    """
    # 먼저 학습 데이터셋 생성하여 스케일러 학습
    train_dataset = BMEDDataset(data_path, train=True, train_ratio=train_ratio)
    
    # 학습된 스케일러 가져오기
    scalers = train_dataset.get_scalers()
    
    # 테스트 데이터셋에 학습된 스케일러 전달
    test_dataset = BMEDDataset(
        data_path, 
        train=False, 
        train_ratio=train_ratio,
        scalers=scalers
    )
    
    # 배치 샘플러 및 데이터 로더 생성
    train_batch_sampler = ExperimentBatchSampler(train_dataset, batch_size)
    
    train_loader = DataLoader(
        train_dataset, 
        batch_sampler=train_batch_sampler
    )
    
    test_loader = DataLoader(
        test_dataset, 
        batch_size=batch_size, 
        shuffle=False
    )
    
    return train_loader, test_loader

def calculate_test_r2_and_create_parity_plots(test_dataset, predictions, actuals, save_dir='./result'):
    """
    테스트 세트의 각 변수별 R² 점수를 계산하고 Parity Plot을 생성하여 저장하는 함수
    
    Parameters:
    -----------
    test_dataset : BMEDDataset
        테스트 데이터셋 (원래 스케일로 변환 시 필요)
    predictions : numpy.ndarray
        테스트 세트에 대한 모델 예측값
    actuals : numpy.ndarray
        테스트 세트의 실제값
    save_dir : str
        결과를 저장할 디렉토리 경로
    """
    # 디렉토리 생성
    os.makedirs(save_dir, exist_ok=True)
    
    # 변수 이름 정의
    var_names = ['T', 'V', 'E', 'CF_LA', 'CF_K', 'CA_LA', 'CB_K', 'VF', 'VA', 'VB']
    
    # 배열 형태 체크
    if predictions.ndim == 3:
        pred_2d = predictions.reshape(-1, predictions.shape[-1])
    else:  # 이미 2차원이면 그대로 사용
        pred_2d = predictions
        
    if actuals.ndim == 3:
        actual_2d = actuals.reshape(-1, actuals.shape[-1])
    else:  # 이미 2차원이면 그대로 사용
        actual_2d = actuals
    
    # 원래 스케일로 변환
    pred_orig = test_dataset.state_scaler.inverse_transform(pred_2d)
    actual_orig = test_dataset.state_scaler.inverse_transform(actual_2d)
    
    # 테스트 세트 R² 점수 저장할 딕셔너리
    r2_scores = {}
    
    # 각 변수별 R² 점수 계산 및 Parity Plot 생성
    for i, var_name in enumerate(var_names):
        # 테스트 세트 R² 계산
        r2 = r2_score(actual_orig[:, i], pred_orig[:, i])
        r2_scores[var_name] = r2
        
        # Parity Plot 생성
        plt.figure(figsize=(8, 8))
        
        # 테스트 데이터 플롯
        plt.scatter(actual_orig[:, i], pred_orig[:, i], alpha=0.5)
        
        # 대각선 (perfect prediction) 추가
        min_val = min(actual_orig[:, i].min(), pred_orig[:, i].min())
        max_val = max(actual_orig[:, i].max(), pred_orig[:, i].max())
        plt.plot([min_val, max_val], [min_val, max_val], 'r--')
        
        # 그래프 제목 및 레이블 설정
        plt.title(f'{var_name} Test Set Parity Plot (R² = {r2:.4f})')
        plt.xlabel(f'Actual {var_name}')
        plt.ylabel(f'Predicted {var_name}')
        plt.grid(True, alpha=0.3)
        
        # 그래프 저장
        plt.tight_layout()
        plt.savefig(os.path.join(save_dir, f'parity_plot_{var_name}.png'), dpi=300)
        plt.close()
    
    # 모든 R² 점수를 텍스트 파일로 저장
    with open(os.path.join(save_dir, 'test_r2_scores.txt'), 'w') as f:
        f.write('Variable,Test R2 Score\n')
        for var_name, r2 in r2_scores.items():
            f.write(f'{var_name},{r2:.6f}\n')
    
    return r2_scores

def save_model_with_scalers(model, scalers, save_dir='result', model_name='bmed_model'):
    """
    모델 구조, 가중치 및 스케일러를 저장하는 함수
    
    Parameters:
    -----------
    model : nn.Module
        저장할 PyTorch 모델
    scalers : dict
        스케일러 딕셔너리 {'feature': feature_scaler, 'mol_change': mol_change_scaler, 'state': state_scaler}
    save_dir : str
        저장할 디렉토리 경로
    model_name : str
        저장할 모델 이름
    """
    # 디렉토리 생성
    os.makedirs(save_dir, exist_ok=True)
    
    # 모델 구조 및 가중치 저장 (전체 모델)
    model_path = os.path.join(save_dir, f'{model_name}.pt')
    torch.save(model, model_path)
    
    # 모델 가중치만 저장 (더 작은 파일 크기)
    weights_path = os.path.join(save_dir, f'{model_name}_weights.pt')
    torch.save(model.state_dict(), weights_path)
    
    # 스케일러 저장
    scalers_path = os.path.join(save_dir, f'{model_name}_scalers.pkl')
    with open(scalers_path, 'wb') as f:
        pickle.dump(scalers, f)
    
    print(f"모델 저장 완료:")
    print(f"- 전체 모델: {model_path}")
    print(f"- 모델 가중치: {weights_path}")
    print(f"- 스케일러: {scalers_path}")

def load_model_with_scalers(save_dir='result', model_name='bmed_model', device='cuda' if torch.cuda.is_available() else 'cpu'):
    """
    저장된 모델과 스케일러를 불러오는 함수
    
    Parameters:
    -----------
    save_dir : str
        저장된 디렉토리 경로
    model_name : str
        저장된 모델 이름
    device : str
        모델을 불러올 디바이스 ('cuda' 또는 'cpu')
        
    Returns:
    --------
    model : nn.Module
        불러온 모델
    scalers : dict
        불러온 스케일러 딕셔너리
    """
    # 모델 경로
    model_path = os.path.join(save_dir, f'{model_name}.pt')
    
    # 스케일러 경로
    scalers_path = os.path.join(save_dir, f'{model_name}_scalers.pkl')
    
    # 모델 불러오기
    model = torch.load(model_path, map_location=device)
    
    # 스케일러 불러오기
    with open(scalers_path, 'rb') as f:
        scalers = pickle.load(f)
    
    return model, scalers

def predict_with_saved_model(input_data, save_dir='result', model_name='bmed_model'):
    """
    저장된 모델을 사용하여 예측을 수행하는 함수
    
    Parameters:
    -----------
    input_data : numpy.ndarray 또는 torch.Tensor
        입력 데이터 (이미 정규화된 형태 또는 원본 형태)
    save_dir : str
        모델이 저장된 디렉토리 경로
    model_name : str
        저장된 모델 이름
        
    Returns:
    --------
    mol_changes : numpy.ndarray
        예측된 몰 변화량
    new_states : numpy.ndarray
        예측된 새로운 상태
    """
    # 모델 및 스케일러 불러오기
    model, scalers = load_model_with_scalers(save_dir, model_name)
    model.eval()
    
    # 입력 데이터가 텐서가 아니면 변환
    if not isinstance(input_data, torch.Tensor):
        # 데이터가 정규화되지 않았다면 정규화
        if hasattr(input_data, 'shape') and len(input_data.shape) == 3:
            # 3D 데이터를 2D로 변환
            orig_shape = input_data.shape
            input_data_2d = input_data.reshape(-1, input_data.shape[-1])
            
            # 정규화
            input_data_scaled = scalers['feature'].transform(input_data_2d)
            
            # 다시 3D로 변환
            input_data = input_data_scaled.reshape(orig_shape)
        
        # numpy 배열을 텐서로 변환
        input_data = torch.FloatTensor(input_data)
    
    # 모델 디바이스에 데이터 넣기
    device = next(model.parameters()).device
    input_data = input_data.to(device)
    
    # 예측 (hidden 상태 추가)
    with torch.no_grad():
        mol_changes, new_states, _ = model(input_data)
    
    # CPU로 이동 및 numpy 변환
    mol_changes = mol_changes.cpu().numpy()
    new_states = new_states.cpu().numpy()
    
    return mol_changes, new_states

def train_model_with_kfold(data_path, n_folds=5, hidden_units=64, hidden_layers=3, 
                          time_step=0.1, batch_size=32, epochs=100, save_dir='result'):
    """
    k-fold cross validation을 사용한 모델 훈련
    """
    fold_results = []
    
    for fold in range(n_folds):
        print(f"\n=== Fold {fold+1}/{n_folds} ===")
        
        # 현재 fold의 훈련 데이터셋 생성
        train_dataset = BMEDDataset(
            data_path=data_path,
            fold_idx=fold,
            n_folds=n_folds,
            train=True
        )
        
        # 스케일러 가져오기
        scalers = train_dataset.get_scalers()
        
        # 현재 fold의 테스트 데이터셋 생성
        test_dataset = BMEDDataset(
            data_path=data_path,
            fold_idx=fold,
            n_folds=n_folds,
            train=False,
            scalers=scalers
        )
        
        # 데이터 로더 생성
        train_batch_sampler = ExperimentBatchSampler(train_dataset, batch_size)
        train_loader = DataLoader(train_dataset, batch_sampler=train_batch_sampler)
        test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)
        
        # 새로운 모델 초기화
        model = MembraneSystemModel(
            hidden_units=hidden_units,
            hidden_layers=hidden_layers,
            time_step=time_step
        )
        
        trainer = MembraneSystemTrainer(model, epochs=epochs)
        
        # 모델 훈련
        train_losses, val_losses = trainer.train(
            train_loader, 
            test_loader,
            patience=100
        )
        
        # 모델 평가
        eval_results = trainer.evaluate_model(test_loader, test_dataset)
        
        # 결과 저장
        fold_results.append({
            'fold': fold,
            'test_loss': eval_results['test_loss'],
            'r2_scores': eval_results['r2_per_variable']
        })
        
        # 현재 fold의 모델 저장
        save_model_with_scalers(
            model, 
            scalers, 
            save_dir=save_dir,
            model_name=f'bmed_model_fold_{fold}'
        )
    
    # k-fold 결과 출력
    print("\n=== K-Fold Cross Validation 결과 ===")
    avg_test_loss = np.mean([r['test_loss'] for r in fold_results])
    print(f"평균 테스트 손실: {avg_test_loss:.6f}")
    
    # 변수별 평균 R² 점수 계산
    var_names = ['T', 'V', 'E', 'CF_LA', 'CF_K', 'CA_LA', 'CB_K', 'VF', 'VA', 'VB']
    for var in var_names:
        avg_r2 = np.mean([r['r2_scores'][var] for r in fold_results])
        print(f"{var} 평균 R² 점수: {avg_r2:.6f}")
    
    return fold_results

# 새로운 재귀적 시뮬레이션 함수 추가
def simulate_with_model(model, initial_state, steps=50, scalers=None, time_step=0.1):
    """
    초기 상태에서 시작하여 모델로 시스템 시뮬레이션을 수행
    
    Parameters:
    -----------
    model : MembraneSystemModel
        학습된 피드포워드 모델
    initial_state : numpy.ndarray
        초기 상태 벡터 [T, V, E, CF_LA, CF_K, CA_LA, CB_K, VF, VA, VB]
    steps : int
        시뮬레이션 단계 수
    scalers : dict
        정규화/역정규화를 위한 스케일러
    time_step : float
        시뮬레이션 시간 간격 (시간)
    
    Returns:
    --------
    states_history : numpy.ndarray
        시뮬레이션된 상태 기록 (steps+1, 10)
    mol_changes_history : numpy.ndarray
        예측된 몰 변화량 기록 (steps, 4)
    t_history : numpy.ndarray
        시뮬레이션 시간 기록 (steps+1)
    """
    model.eval()
    device = next(model.parameters()).device
    
    # 초기 상태를 NumPy 배열로 변환
    if isinstance(initial_state, torch.Tensor):
        initial_state = initial_state.cpu().numpy()
    initial_state = np.asarray(initial_state).reshape(-1)
    
    # 초기 상태 정규화
    if scalers is not None:
        initial_state_norm = scalers['feature'].transform(initial_state.reshape(1, -1))
        initial_state_norm = initial_state_norm.reshape(-1)
    else:
        initial_state_norm = initial_state
    
    # 현재 상태: [1, 10] 형태 (2D)
    current_state = torch.FloatTensor(initial_state_norm).reshape(1, -1)

    # 모델 입력을 위해 3D로 변환 (배치 차원, 시퀀스 차원, 특성 차원)
    current_state = current_state.unsqueeze(1)  # [1, 1, 10] 형태 (3D)
    
    # 결과 저장 리스트
    if scalers is not None:
        states_history = [initial_state.copy()]  # 원래 스케일로 저장
    else:
        states_history = [initial_state.copy()]
    
    mol_changes_history = []
    t_history = [0.0]  # 시간 기록 시작
    current_time = 0.0
    
    # 시뮬레이션 실행
    with torch.no_grad():
        for step in range(steps):
            # 예측 (hidden 상태 없음)
            mol_changes, new_state, _ = model(current_state)
            
            # 결과 저장 (역정규화)
            if scalers is not None:
                mol_changes_np = mol_changes.cpu().squeeze().numpy().reshape(1, -1)
                new_state_np = new_state.reshape(1, -1).cpu().numpy()
                
                mol_changes_orig = scalers['mol_change'].inverse_transform(mol_changes_np)[0]
                next_state_orig = scalers['state'].inverse_transform(new_state_np)[0]
                
                mol_changes_history.append(mol_changes_orig)
                states_history.append(next_state_orig)
            else:
                mol_changes_history.append(mol_changes.cpu().squeeze().numpy())
                states_history.append(new_state.cpu().squeeze().numpy())
            
            # 다음 상태를 입력으로 설정
            current_state = torch.FloatTensor(new_state.cpu().numpy()).reshape(1, 1, -1).to(device)
            
            # 시간 업데이트
            current_time += time_step
            t_history.append(current_time)
    
    return np.array(states_history), np.array(mol_changes_history), np.array(t_history)

# 시뮬레이션 결과 시각화 함수 추가
def plot_simulation_results(states_history, t_history, experiment_data=None, save_path=None):
    """
    시뮬레이션 결과를 시각화하는 함수
    
    Parameters:
    -----------
    states_history : numpy.ndarray
        시뮬레이션 상태 히스토리
    t_history : numpy.ndarray
        시뮬레이션 시간 히스토리
    experiment_data : pandas.DataFrame
        비교를 위한 실험 데이터 (선택 사항)
    save_path : str
        결과를 저장할 경로 (선택 사항)
    """
    # 변수 이름
    var_names = ['T', 'V', 'E', 'CF_LA', 'CF_K', 'CA_LA', 'CB_K', 'VF', 'VA', 'VB']
    var_labels = {
        'T': 'Temperature [°C]',
        'V': 'Voltage [V]',
        'E': 'Electrolyte Conc. [mol/L]',
        'CF_LA': 'Feed LA Conc. [mol/L]',
        'CF_K': 'Feed K+ Conc. [mol/L]',
        'CA_LA': 'Acid LA Conc. [mol/L]',
        'CB_K': 'Base K+ Conc. [mol/L]',
        'VF': 'Feed Volume [L]',
        'VA': 'Acid Volume [L]',
        'VB': 'Base Volume [L]'
    }
    
    # 그래프 생성
    fig, axes = plt.subplots(5, 2, figsize=(15, 20))
    axes = axes.flatten()
    
    for i, var in enumerate(var_names):
        ax = axes[i]
        
        # 시뮬레이션 결과 그리기
        ax.plot(t_history, states_history[:, i], 'b-', linewidth=2, label='Simulation')
        
        # 실험 데이터가 있으면 함께 그리기
        if experiment_data is not None:
            ax.scatter(experiment_data['t'], experiment_data[var], 
                       color='red', marker='o', label='Experiment')
        
        ax.set_xlabel('Time [h]')
        ax.set_ylabel(var_labels.get(var, var))
        ax.set_title(f'{var} vs Time')
        ax.grid(True, alpha=0.3)
        ax.legend()
    
    plt.tight_layout()
    
    # 저장 경로가 있으면 저장
    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"그래프가 {save_path}에 저장되었습니다.")
    
    plt.show()

# 시뮬레이션 실행 함수 추가 (편의 함수)
def run_simulation_from_saved_model(initial_state, steps=50, time_step=0.1, 
                                   save_dir='result', model_name='bmed_model',
                                   experiment_data=None, plot_results=True, save_plot_path=None):
    """
    저장된 모델을 불러와 시뮬레이션 실행 및 결과 시각화
    
    Parameters:
    -----------
    initial_state : numpy.ndarray
        초기 상태 벡터
    steps : int
        시뮬레이션 단계 수
    time_step : float
        시뮬레이션 시간 간격
    save_dir : str
        모델이 저장된 디렉토리
    model_name : str
        모델 이름
    experiment_data : pandas.DataFrame
        비교를 위한 실험 데이터 (선택 사항)
    plot_results : bool
        결과를 그래프로 표시할지 여부
    save_plot_path : str
        그래프 저장 경로 (선택 사항)
    
    Returns:
    --------
    states_history : numpy.ndarray
        시뮬레이션된 상태 기록
    mol_changes_history : numpy.ndarray
        예측된 몰 변화량 기록
    t_history : numpy.ndarray
        시간 기록
    """
    # 모델 및 스케일러 불러오기
    model, scalers = load_model_with_scalers(save_dir, model_name)
    
    # 시뮬레이션 실행
    states_history, mol_changes_history, t_history = simulate_with_model(
        model, initial_state, steps, scalers, time_step
    )
    
    # 결과 시각화
    if plot_results:
        plot_simulation_results(states_history, t_history, experiment_data, save_plot_path)
    
    return states_history, mol_changes_history, t_history

# 사용 예시 함수
def example_usage(data_path):
    """
    모델 학습, 저장 및 불러오기 예시
    """
    # 1. 모델 학습 및 저장
    model, trainer, _, _, _ = train_model_from_csv(
        data_path=data_path,
        hidden_units=128,        # 은닉층 유닛 수
        hidden_layers=4,         # 은닉층 수
        time_step=0.1,
        batch_size=32,
        epochs=100,
        train_ratio=0.8,
    )
    
    # 2. 저장된 모델 불러오기
    loaded_model, scalers = load_model_with_scalers()
    
    # 3. 테스트 데이터 생성
    test_dataset = BMEDDataset(data_path, train=False, scalers=scalers)
    test_loader = DataLoader(test_dataset, batch_size=1)
    
    # 4. 모델 예측
    for features, _, _ in test_loader:
        mol_changes, new_states = predict_with_saved_model(features.numpy())
        print("Predicted mol changes shape:", mol_changes.shape)
        print("Predicted new states shape:", new_states.shape)
        break

if __name__ == "__main__":
    data_path = './data/BMED_data_v6+spline.xlsx'
    
    # 5-fold cross validation 실행
    fold_results = train_model_with_kfold(
        data_path=data_path,
        n_folds=5,
        hidden_units=256,
        hidden_layers=20,
        time_step=0.1,
        batch_size=16,
        epochs=1000
    )