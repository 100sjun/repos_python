# Load Modules
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader
from sklearn.metrics import r2_score

class BMEDDataset(Dataset):
    '''BMED 데이터로부터 PyTorch 데이터셋 생성'''
    def __init__(self, data_path, sequence_length=10, train=True, train_ratio=0.8):
        self.sequence_length = sequence_length
        self.train = train

        # Spline 데이터 로드
        self.df = pd.read_excel(data_path, sheet_name='spline_data')

        # 실험 번호별로 데이터 분리
        self.dic_spline = {}
        for exp in self.df['exp'].unique():
            self.dic_spline[exp] = self.df[self.df['exp'] == exp].sort_values('t')
        
        # 훈련/테스트 셋 분리 (실험 번호 기준)
        all_exps = list(self.dic_spline.keys())
        np.random.seed(42)
        np.random.shuffle(all_exps)
        np.random.seed(None)
        
        split_idx = int(len(all_exps) * train_ratio)
        
        if train:
            self.exps_to_use = all_exps[:split_idx]
        else:
            self.exps_to_use = all_exps[split_idx:]
            
        # Data preparation
        self.prepare_data()

        # Save the indices of each experiments for batch processing
        self.exp_indices = []
        cur_idx = 0
        for exp_id in self.exps_to_use:
            exp_data = self.dic_spline[exp_id]
            n_samples = len(exp_data) - sequence_length
            indices = list(range(cur_idx, cur_idx + n_samples))
            self.exp_indices.append((exp_id, indices))
            cur_idx += n_samples

    def prepare_data(self):
        features_list = []
        mol_changes_targets_list= []
        states_targets_list = []

        for exp in self.exps_to_use:
            exp_data = self.dic_spline[exp]

            # Input data as LSTM sequence
            for i in range(len(exp_data) - self.sequence_length):
                seq_data = exp_data.iloc[i:i+self.sequence_length]

                # Feature vector
                features = []
                for _, row in seq_data.iterrows():
                    feature = [
                        row['T'], row['V'], row['E'],
                        row['CF_LA'], row['CF_K'], row['CA_LA'], row['CB_K'],
                        row['VF'], row['VA'], row['VB']
                    ]
                    features.append(feature)
                
                # Target vectors
                mol_change_targets = []
                state_targets = []
                for j in range(self.sequence_length-1):
                    cur = seq_data.iloc[j]
                    next = seq_data.iloc[j+1]
                    
                    # Time interval
                    dt = next['t'] - cur['t']

                    # Mole of LA and K+ (mol)
                    cur_NA_LA = cur['CA_LA'] * cur['VA']
                    cur_NB_K = cur['CB_K'] * cur['VB']
                    next_NA_LA = next['CA_LA'] * next['VA']
                    next_NB_K = next['CB_K'] * next['VB']

                    # Mole change (mol/hr)
                    dNLA = (next_NA_LA - cur_NA_LA) / dt
                    dNK = (next_NB_K - cur_NB_K) / dt

                    # Water change (L/hr)
                    dVA = (next['VA'] - cur['VA']) / dt
                    dVB = (next['VB'] - cur['VB']) / dt

                    # Change vector
                    mol_change = [dNLA, dNK, dVA, dVB]
                    mol_change_targets.append(mol_change)

                    # State vector
                    state = [
                        next['T'], next['V'], next['E'],
                        next['CF_LA'], next['CF_K'], next['CA_LA'], next['CB_K'],
                        next['VF'], next['VA'], next['VB']
                    ]
                    state_targets.append(state)

                # End point of the sequence
                cur = seq_data.iloc[-1]
                if i + self.sequence_length < len(exp_data):
                    next = exp_data.iloc[i + self.sequence_length]

                    # Time interval 
                    dt = next['t'] - cur['t']

                    # Mole of LA and K+ (mol)
                    cur_NA_LA = cur['CA_LA'] * cur['VA']
                    cur_NB_K = cur['CB_K'] * cur['VB']
                    next_NA_LA = next['CA_LA'] * next['VA']
                    next_NB_K = next['CB_K'] * next['VB']

                    # Mole change (mol/hr)
                    dNLA = (next_NA_LA - cur_NA_LA) / dt
                    dNK = (next_NB_K - cur_NB_K) / dt

                    # Water change (L/hr)
                    dVA = (next['VA'] - cur['VA']) / dt
                    dVB = (next['VB'] - cur['VB']) / dt

                    # Change vector
                    mol_change = [dNLA, dNK, dVA, dVB]
                    mol_change_targets.append(mol_change)

                    # State vector
                    state = [
                        next['T'], next['V'], next['E'],
                        next['CF_LA'], next['CF_K'], next['CA_LA'], next['CB_K'],
                        next['VF'], next['VA'], next['VB']
                    ]
                    state_targets.append(state)
                else:
                    mol_change_targets.append(mol_change_targets[-1])
                    state_targets.append(state_targets[-1])

                features_list.append(features)
                mol_changes_targets_list.append(mol_change_targets)
                states_targets_list.append(state_targets)
        self.features = np.array(features_list, dtype=np.float32)
        self.mol_changes_targets = np.array(mol_changes_targets_list, dtype=np.float32)
        self.states_targets = np.array(states_targets_list, dtype=np.float32)

        print(f"Data loaded: {len(self.features)} sequences, feature shape: {self.features.shape}")
        print(f"Mol change targets shape: {self.mol_changes_targets.shape}")
        print(f"State targets shape: {self.states_targets.shape}")

    def __len__(self):
        return len(self.features)
                        
    def __getitem__(self, idx):
        return (
            torch.FloatTensor(self.features[idx]),
            torch.FloatTensor(self.mol_changes_targets[idx]),
            torch.FloatTensor(self.states_targets[idx])
        )
    
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
    '''멤브레인 시스템 모델링을 위한 Physics-Informed LSTM 모델'''
    def __init__(self, lstm_units=64, lstm_layer=3, fc_units=64, fc_layer=5,time_step=0.1, sequence_length=10):
        super(MembraneSystemModel, self).__init__()
        self.lstm_units = lstm_units
        self.lstm_layer = lstm_layer
        self.fc_units = fc_units
        self.fc_layer = fc_layer
        self.time_step = time_step
        self.sequence_length = sequence_length

        '''LSMT layer'''
        self.lstm = nn.LSTM(
            input_size=10, # [T, V, E, CF_LA, CF_K, CA_LA, CB_K, VF, VA, VB]
            hidden_size=self.lstm_units,
            num_layers = self.lstm_layer,
            batch_first=True, #(batch, seq, feature) 형태 입력
        )

        '''FC layer'''
        self.fc = nn.Sequential(
            nn.Linear(self.lstm_units, self.fc_units),
            nn.ReLU(),
            *[nn.Sequential(
                nn.Linear(self.fc_units, self.fc_units),
                nn.ReLU()
            ) for _ in range(self.fc_layer - 2)],
            nn.Linear(self.fc_units, 4)  # 4개의 flux 값 출력 (LA 몰 변화량, K+ 몰 변화량, Acid 방향 물 변화량, Base 방향 물 변화량)
        )

        '''Physics layer'''
        self.physics_layer = PhysicsLayer(time_step=self.time_step)

    def forward(self, x):
        '''LSTM 처리'''
        lstm_out, _ = self.lstm(x)

        '''몰 변화량 예측 - 각 시간 단계별로 동일한 FC 적용'''
        batch_size, seq_len, _ = lstm_out.size()
        mol_changes_list = []

        for i in range(seq_len): # 시간별 LSTM out으로 flux 예측
            step_lstm_out = lstm_out[:, i, :]
            step_mol_changes = self.fc(step_lstm_out)
            mol_changes_list.append(step_mol_changes.unsqueeze(1))
        
        # 모든 시간 단계의 몰 변화량 예측 결과 결합
        mol_changes = torch.cat(mol_changes_list, dim=1)

        '''물리 법칙 적용: 질량 보존'''
        new_states_list = []

        for i in range(seq_len):
            new_state = self.physics_layer(mol_changes[:, i, :], x[:, i, :]) # x는 현재 상태의 state를 의미함
            new_states_list.append(new_state.unsqueeze(1))
        
        # 모든 시간 단계 결과 결합
        new_states = torch.cat(new_states_list, dim=1)

        return mol_changes, new_states

class MembraneSystemTrainer:
    """
    멤브레인 시스템 모델 훈련 및 시뮬레이션을 위한 클래스
    """
    def __init__(self, model, device='cuda' if torch.cuda.is_available() else 'cpu',epochs=100):
        self.model = model.to(device)
        self.device = device
        self.epochs = epochs
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=0.001)
        # 손실 함수는 상태함수를 기준으로 실시
        self.criterion = nn.MSELoss()
    
    def train(self, train_loader, val_loader=None):
        """
        모델 훈련
        
        입력:
        - train_loader: 훈련 데이터 로더
        - val_loader: 검증 데이터 로더 (옵션)
        - epochs: 훈련 에폭 수
        - mol_change_weight: 몰 변화량 손실 가중치
        - state_weight: 상태 손실 가중치
        """
        self.model.train()
        train_losses = []
        val_losses = []

        for epoch in range(self.epochs):
            epoch_loss = 0.0
            
            for batch_idx, (features, mol_change_targets, state_targets) in enumerate(train_loader):
                features = features.to(self.device)
                mol_change_targets = mol_change_targets.to(self.device)
                state_targets = state_targets.to(self.device)

                # initialize gradients
                self.optimizer.zero_grad()

                # forward pass
                _, state_predictions = self.model(features)

                # loss calculation
                loss = self.criterion(state_predictions, state_targets)

                # back propagation and optimization
                loss.backward()
                self.optimizer.step()

                epoch_loss += loss.item()

            # epoch average loss
            avg_loss = epoch_loss / len(train_loader)
            train_losses.append(avg_loss)

            if val_loader is not None:
                self.model.eval()
                val_loss = 0.0

                with torch.no_grad():
                    for features, mol_change_targets, state_targets in val_loader:
                        features = features.to(self.device)
                        mol_change_targets = mol_change_targets.to(self.device)
                        state_targets = state_targets.to(self.device)
                        
                        _, state_predictions = self.model(features)
                        
                        loss = self.criterion(state_predictions, state_targets)
                        val_loss += loss.item()
                
                avg_val_loss = val_loss / len(val_loader)
                val_losses.append(avg_val_loss)

                print(f'Epoch {epoch+1}/{self.epochs}, Train Loss: {avg_loss:.6f}, Val Loss: {avg_val_loss:.6f}')
                
                self.model.train()
            else:
                print(f'Epoch {epoch+1}/{self.epochs}, Train Loss: {avg_loss:.6f}')

        return train_losses, val_losses
    
    def evaluate_model(self, test_loader):
        self.model.eval()
        test_loss = 0
        predictions = []
        actuals = []

        with torch.no_grad():
            for features, mol_change_targets, state_targets in test_loader:
                features = features.to(self.device)
                mol_change_targets = mol_change_targets.to(self.device)
                state_targets = state_targets.to(self.device)

                # 모델 예측
                mol_change_predictions, state_predictions = self.model(features)

                # 손실 계산
                loss = self.criterion(state_predictions, state_targets)
                test_loss += loss.item()

                predictions.append(state_predictions.cpu().numpy())
                actuals.append(state_targets.cpu().numpy())
                
                # 예측값과 실제값 저장
                predictions.append(state_predictions.cpu().numpy())
                actuals.append(state_targets.cpu().numpy())
         
        # 전체 테스트 세트에 대한 평균 손실
        avg_test_loss = test_loss / len(test_loader)                
                
        # 예측값과 실제값을 하나의 배열로 결합
        predictions = np.concatenate(predictions, axis=0)
        actuals = np.concatenate(actuals, axis=0)

        # R2 점수 계산
        r2_scores = {}
        for i, var_name in enumerate(['CF_LA', 'CF_K', 'CA_LA', 'acid_k', 'base_la', 'CB_K', 'VF', 'VA', 'VB']):
            r2 = r2_score(actuals[:, :, i], predictions[:, :, i])
            r2_scores[var_name] = r2
        
        return {
            'test_loss': avg_test_loss,
            'predictions': predictions,
            'actuals': actuals,
            'r2_scores': r2_scores
        }

class ExperimentBatchSampler(torch.utils.data.Sampler):
    """
    실험 단위로 배치를 구성하는 샘플러
    """
    def __init__(self, dataset, batch_size):
        self.dataset = dataset
        self.batch_size = batch_size
        
    def __iter__(self):
        # 실험별로 인덱스를 섞어 배치 구성
        all_indices = []
        for exp_id, indices in self.dataset.exp_indices:
            # 각 실험 내 인덱스는 순서대로 유지
            all_indices.extend(indices)
            
        # 배치 크기에 맞게 인덱스 그룹화
        batches = [all_indices[i:i + self.batch_size] 
                   for i in range(0, len(all_indices), self.batch_size)]
        
        # 배치 순서는 섞음
        np.random.shuffle(batches)
        
        for batch in batches:
            yield batch
    
    def __len__(self):
        return (len(self.dataset) + self.batch_size - 1) // self.batch_size
    
def data_loaders_from_csv(data_path, sequence_length=10, batch_size=2, train_ratio=0.8):
    """
    CSV 파일에서 훈련 및 검증 데이터 로더 생성
    """
    train_dataset = BMEDDataset(data_path, sequence_length=sequence_length, train=True, train_ratio=train_ratio)
    test_dataset = BMEDDataset(data_path, sequence_length=sequence_length, train=False, train_ratio=train_ratio)
    
    # 실험 단위 배치 샘플러 사용
    train_batch_sampler = ExperimentBatchSampler(train_dataset, batch_size)
    
    train_loader = DataLoader(
        train_dataset, 
        batch_sampler=train_batch_sampler
    )
    
    # 테스트 데이터는 순서대로 (shuffle=False)
    test_loader = DataLoader(
        test_dataset, 
        batch_size=batch_size, 
        shuffle=False
    )
    
    return train_loader, test_loader

def train_model_from_csv(data_path, lstm_units=64, lstm_layer = 3, fc_units=64, fc_layer=5, time_step=0.1, sequence_length=10, 
                        batch_size=32, epochs=100, train_ratio=0.8):
    """
    CSV 파일에서 모델 훈련을 위한 편의 함수
    """
    # 데이터 로더 생성
    train_loader, test_loader = data_loaders_from_csv(
        data_path, sequence_length, batch_size, train_ratio
    )
    # 모델 생성
    model = MembraneSystemModel(lstm_units=lstm_units, lstm_layer=lstm_layer, fc_units=fc_units, fc_layer=fc_layer, time_step=time_step, sequence_length=sequence_length)
    trainer = MembraneSystemTrainer(model,epochs=epochs)
    
    # 모델 훈련
    train_losses, val_losses = trainer.train(
        train_loader, 
        test_loader
    )
    
    return model, trainer, train_losses, val_losses

if __name__ == "__main__":
    # CSV 파일 경로
    data_path = '../data/BMED_data_v6+spline.xlsx'
    
    # 모델 훈련
    model, trainer, train_losses, val_losses = train_model_from_csv(
        data_path=data_path,
        lstm_units=64,
        lstm_layer=3,
        fc_units=64,
        fc_layer=5,
        sequence_length=3,  # 이전 3개 시점 데이터 사용
        epochs = 10,
        time_step=0.1,
        batch_size=2,
        train_ratio=0.8,
    )