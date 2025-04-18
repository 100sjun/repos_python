import torch
import torch.nn as nn
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from torch.utils.data import Dataset, DataLoader
from sklearn.metrics import r2_score
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
import os

class PhysicsLayer(nn.Module):
    """
    물리적 법칙을 반영하는 레이어
    몰 변화량 예측값을 받아 물질수지와 부피 변화를 계산
    """
    def __init__(self, time_step=0.1):
        super(PhysicsLayer, self).__init__()
        self.time_step = time_step  # 시간 단위 [hour]

    def forward(self, mol_changes, features):
        # 현재 시간 단계 상태 추출
        current_temp = features[..., 0:1]
        current_voltage = features[..., 1:2]
        current_electrolyte = features[..., 2:3]
        
        # 현재 부피 추출
        current_feed_volume = features[..., 5:6]  # L
        current_acid_volume = features[..., 6:7]  # L
        current_base_volume = features[..., 7:8]  # L
        
        # 현재 농도 추출 (mol/L)
        current_feed_la = features[..., 8:9]  
        current_feed_k = features[..., 9:10]
        current_acid_la = features[..., 10:11]
        # current_acid_k는 항상 0 (Acid에는 K+가 들어갈 수 없음)
        # current_base_la는 항상 0 (Base에는 LA가 들어갈 수 없음)
        current_base_k = features[..., 11:12]
        
        # 몰 변화량 추출 (mol/h)
        la_mol_change = mol_changes[..., 0:1]  # Feed->Acid 방향 LA 몰 변화량
        k_mol_change = mol_changes[..., 1:2]   # Feed->Base 방향 K+ 몰 변화량
        water_acid_change = mol_changes[..., 2:3]  # Feed->Acid 방향 물 부피 변화량 (L/h)
        water_base_change = mol_changes[..., 3:4]  # Feed->Base 방향 물 부피 변화량 (L/h)
        
        # 한 시간 단계 동안의 물질 전달량 계산 (mol)
        la_transfer = la_mol_change * self.time_step
        k_transfer = k_mol_change * self.time_step
        
        # 한 시간 단계 동안의 물 전달량 계산 (L)
        water_to_acid = water_acid_change * self.time_step
        water_to_base = water_base_change * self.time_step
        
        # 부피 업데이트
        new_feed_volume = current_feed_volume - water_to_acid - water_to_base
        new_acid_volume = current_acid_volume + water_to_acid
        new_base_volume = current_base_volume + water_to_base
        
        # 몰수 계산 및 업데이트
        feed_la_mol = current_feed_la * current_feed_volume
        feed_k_mol = current_feed_k * current_feed_volume
        acid_la_mol = current_acid_la * current_acid_volume
        acid_k_mol = torch.zeros_like(acid_la_mol)  # Acid에는 K+가 없음
        base_la_mol = torch.zeros_like(current_base_k)  # Base에는 LA가 없음
        base_k_mol = current_base_k * current_base_volume
        
        # 몰수 전달 적용
        new_feed_la_mol = feed_la_mol - la_transfer
        new_feed_k_mol = feed_k_mol - k_transfer
        new_acid_la_mol = acid_la_mol + la_transfer
        new_acid_k_mol = torch.zeros_like(acid_la_mol)  # 항상 0
        new_base_la_mol = torch.zeros_like(current_base_k)  # 항상 0
        new_base_k_mol = base_k_mol + k_transfer
        
        # 새로운 농도 계산 (0으로 나누기 방지)
        eps = 1e-6  # 작은 값 추가하여 0으로 나누기 방지
        new_feed_la_conc = new_feed_la_mol / (new_feed_volume + eps)
        new_feed_k_conc = new_feed_k_mol / (new_feed_volume + eps)
        new_acid_la_conc = new_acid_la_mol / (new_acid_volume + eps)
        new_acid_k_conc = torch.zeros_like(new_acid_la_conc)  # 항상 0
        new_base_la_conc = torch.zeros_like(new_base_k_mol)  # 항상 0
        new_base_k_conc = new_base_k_mol / (new_base_volume + eps)
        
        # 결과 출력: 새로운 농도, 새로운 부피
        new_states = torch.cat([
            new_feed_la_conc, new_feed_k_conc,
            new_acid_la_conc, new_acid_k_conc,
            new_base_la_conc, new_base_k_conc,
            new_feed_volume, new_acid_volume, new_base_volume
        ], dim=-1)
        
        return new_states
    
class MembraneSystemModel(nn.Module):
    """
    멤브레인 시스템 모델링을 위한 Physics-Informed LSTM 모델
    """
    def __init__(self, lstm_units=64, mlp_hidden_units=[128, 64], time_step=0.1, sequence_length=3):
        super(MembraneSystemModel, self).__init__()
        self.lstm_units = lstm_units
        self.mlp_hidden_units = mlp_hidden_units
        self.time_step = time_step
        self.sequence_length = sequence_length
        
        # LSTM 층
        self.lstm = nn.LSTM(
            input_size=12,  # [온도, 전압, 전해질, Ci, Ki, VF, VA, VB, CF_LA, CF_K, CA_LA, CB_K]
            hidden_size=self.lstm_units,
            batch_first=True,  # (batch, seq, feature) 형태 입력
            num_layers=1
        )
        
        # 몰 변화량 예측을 위한 다층 퍼셉트론 정의
        mlp_layers = []
        input_size = self.lstm_units
        
        # MLP 히든 레이어 구성
        for hidden_size in self.mlp_hidden_units:
            mlp_layers.append(nn.Linear(input_size, hidden_size))
            mlp_layers.append(nn.ReLU())
            input_size = hidden_size
        
        # 출력 레이어 (4개 값: LA몰 변화량, K+몰 변화량, 물->Acid 부피 변화량, 물->Base 부피 변화량)
        mlp_layers.append(nn.Linear(input_size, 4))
        
        # MLP 레이어를 Sequential 컨테이너로 구성
        self.mol_change_mlp = nn.Sequential(*mlp_layers)
        
        # 물리 법칙 층
        self.physics_layer = PhysicsLayer(time_step=self.time_step)
        
    def forward(self, x):
        # LSTM 처리
        lstm_out, _ = self.lstm(x)
        
        # 몰 변화량 예측 - 각 시간 단계별로 동일한 MLP 적용
        batch_size, seq_len, _ = lstm_out.size()
        mol_changes_list = []
        
        for i in range(seq_len):
            # 각 시간 단계의 LSTM 출력을 MLP에 입력
            step_lstm_out = lstm_out[:, i, :]
            step_mol_changes = self.mol_change_mlp(step_lstm_out)
            mol_changes_list.append(step_mol_changes.unsqueeze(1))
        
        # 모든 시간 단계의 몰 변화량 예측 결과 결합
        mol_changes = torch.cat(mol_changes_list, dim=1)
        
        # 물리 법칙 적용 - 모든 시간 단계를 처리하기 위한 루프
        new_states_list = []
        
        for i in range(seq_len):
            new_state = self.physics_layer(mol_changes[:, i, :], x[:, i, :])
            new_states_list.append(new_state.unsqueeze(1))
        
        # 모든 시간 단계의 결과를 결합
        new_states = torch.cat(new_states_list, dim=1)
        
        return mol_changes, new_states
    
class BMEDDataset(Dataset):
    """
    BMED CSV 데이터를 위한 PyTorch 데이터셋
    """
    def __init__(self, csv_path, sequence_length=1, train=True, train_ratio=0.8):
        # CSV 파일 로드
        self.df = pd.read_csv(csv_path)
        
        # 실험 번호별로 데이터 분리
        self.experiment_data = {}
        for exp in self.df['exp'].unique():
            self.experiment_data[exp] = self.df[self.df['exp'] == exp].sort_values('t')
        
        # 훈련/테스트 셋 분리 (실험 번호 기준)
        all_exps = list(self.experiment_data.keys())
        # 랜덤 시드 고정으로 동일한 분할 보장
        np.random.seed(42)  
        np.random.shuffle(all_exps)
        np.random.seed(None)  # 랜덤 시드 초기화
        
        split_idx = int(len(all_exps) * train_ratio)
        
        if train:
            self.exps_to_use = all_exps[:split_idx]
        else:
            self.exps_to_use = all_exps[split_idx:]
        
        self.sequence_length = sequence_length
        self.train = train
        
        # 데이터 준비
        self.prepare_data()
        
        # 실험별 시작 인덱스 저장 - 배치 구성 시 사용
        self.exp_indices = []
        current_idx = 0
        self.exp_id_to_indices = {}
        
        for exp_id in self.exps_to_use:
            exp_data = self.experiment_data[exp_id]
            
            if len(exp_data) <= self.sequence_length:
                continue
                
            n_samples = len(exp_data) - self.sequence_length
            indices = list(range(current_idx, current_idx + n_samples))
            self.exp_id_to_indices[exp_id] = indices
            self.exp_indices.append((exp_id, indices))
            current_idx += n_samples
    
    def prepare_data(self):
        # 특성 및 목표 변수 데이터
        features_list = []
        mol_change_targets_list = []
        state_targets_list = []
        
        for exp in self.exps_to_use:
            exp_data = self.experiment_data[exp]
            
            # 시간 간격이 일정하지 않을 수 있으므로 시간 순으로 정렬된 데이터 사용
            if len(exp_data) <= self.sequence_length:
                continue
                
            for i in range(len(exp_data) - self.sequence_length):
                # 시퀀스 길이만큼의 입력 데이터
                seq_data = exp_data.iloc[i:i+self.sequence_length]
                
                # 특성 벡터 구성
                features = []
                for _, row in seq_data.iterrows():
                    # [온도, 전압, 전해질, 초기LA농도, 초기K+농도, 
                    #  Feed부피, Acid부피, Base부피, 
                    #  Feed LA농도, Feed K+농도, Acid LA농도, Base K+농도]
                    feature = [
                        row['T'], row['V'], row['E'], 
                        row['Ci'], row['Ki'],
                        row['VF'], row['VA'], row['VB'],
                        row['CF_LA'], row['CF_K'], row['CA_LA'], row['CB_K']
                    ]
                    features.append(feature)
                
                # 시퀀스 내 각 시간 단계에 대한 목표 변수 준비
                # (각 시간 단계의 t+1 데이터가 필요함)
                mol_change_targets = []
                state_targets = []
                
                for j in range(self.sequence_length-1):
                    current = seq_data.iloc[j]  # 현재 시간 단계
                    next_row = seq_data.iloc[j+1]  # 다음 시간 단계
                    
                    # 시간 간격
                    dt = next_row['t'] - current['t']
                    
                    # LA 및 K+ 몰수 계산 (농도 * 부피)
                    current_la_feed = current['CF_LA'] * current['VF']
                    current_la_acid = current['CA_LA'] * current['VA']
                    current_k_feed = current['CF_K'] * current['VF']
                    current_k_base = current['CB_K'] * current['VB']
                    
                    next_la_feed = next_row['CF_LA'] * next_row['VF']
                    next_la_acid = next_row['CA_LA'] * next_row['VA']
                    next_k_feed = next_row['CF_K'] * next_row['VF']
                    next_k_base = next_row['CB_K'] * next_row['VB']
                    
                    # 몰 변화량 계산 (mol/h)
                    la_mol_change = (next_la_acid - current_la_acid) / dt
                    k_mol_change = (next_k_base - current_k_base) / dt
                    
                    # 부피 변화량 계산 (L/h)
                    water_acid_change = (next_row['VA'] - current['VA']) / dt
                    water_base_change = (next_row['VB'] - current['VB']) / dt
                    
                    mol_change = [la_mol_change, k_mol_change, water_acid_change, water_base_change]
                    mol_change_targets.append(mol_change)
                    
                    # 상태 변수: 농도와 부피
                    state = [
                        next_row['CF_LA'], next_row['CF_K'],
                        next_row['CA_LA'], 0,  # Acid에는 K+가 없음
                        0, next_row['CB_K'],  # Base에는 LA가 없음
                        next_row['VF'], next_row['VA'], next_row['VB']
                    ]
                    state_targets.append(state)
                
                # 마지막 시간 단계와 그 다음 데이터
                current = seq_data.iloc[-1]  # 시퀀스의 마지막 시간 단계
                if i + self.sequence_length < len(exp_data):
                    next_row = exp_data.iloc[i + self.sequence_length]  # 다음 시간 단계
                    
                    # 시간 간격
                    dt = next_row['t'] - current['t']
                    
                    # LA 및 K+ 몰수 계산 (농도 * 부피)
                    current_la_feed = current['CF_LA'] * current['VF']
                    current_la_acid = current['CA_LA'] * current['VA']
                    current_k_feed = current['CF_K'] * current['VF']
                    current_k_base = current['CB_K'] * current['VB']
                    
                    next_la_feed = next_row['CF_LA'] * next_row['VF']
                    next_la_acid = next_row['CA_LA'] * next_row['VA']
                    next_k_feed = next_row['CF_K'] * next_row['VF']
                    next_k_base = next_row['CB_K'] * next_row['VB']
                    
                    # 몰 변화량 계산 (mol/h)
                    la_mol_change = (next_la_acid - current_la_acid) / dt
                    k_mol_change = (next_k_base - current_k_base) / dt
                    
                    # 부피 변화량 계산 (L/h)
                    water_acid_change = (next_row['VA'] - current['VA']) / dt
                    water_base_change = (next_row['VB'] - current['VB']) / dt
                    
                    mol_change = [la_mol_change, k_mol_change, water_acid_change, water_base_change]
                    mol_change_targets.append(mol_change)
                    
                    # 상태 변수: 농도와 부피
                    state = [
                        next_row['CF_LA'], next_row['CF_K'],
                        next_row['CA_LA'], 0,  # Acid에는 K+가 없음
                        0, next_row['CB_K'],  # Base에는 LA가 없음
                        next_row['VF'], next_row['VA'], next_row['VB']
                    ]
                    state_targets.append(state)
                else:
                    # 시퀀스의 마지막 데이터가 데이터셋의 마지막이면, 이전 단계 변화량을 그대로 사용
                    if len(mol_change_targets) > 0:
                        mol_change_targets.append(mol_change_targets[-1])
                        state_targets.append(state_targets[-1])
                    else:
                        # 시퀀스가 길이 1인 경우 처리
                        mol_change_targets.append([0, 0, 0, 0])
                        state = [
                            current['CF_LA'], current['CF_K'],
                            current['CA_LA'], 0,
                            0, current['CB_K'],
                            current['VF'], current['VA'], current['VB']
                        ]
                        state_targets.append(state)
                
                features_list.append(features)
                mol_change_targets_list.append(mol_change_targets)
                state_targets_list.append(state_targets)
        
        self.features = np.array(features_list, dtype=np.float32)
        self.mol_change_targets = np.array(mol_change_targets_list, dtype=np.float32)
        self.state_targets = np.array(state_targets_list, dtype=np.float32)
        
        print(f"Data loaded: {len(self.features)} sequences, feature shape: {self.features.shape}")
        print(f"Mol change targets shape: {self.mol_change_targets.shape}")
        print(f"State targets shape: {self.state_targets.shape}")
    
    def __len__(self):
        return len(self.features)
    
    def __getitem__(self, idx):
        return (
            torch.FloatTensor(self.features[idx]),
            torch.FloatTensor(self.mol_change_targets[idx]),
            torch.FloatTensor(self.state_targets[idx])
        )

class MembraneSystemTrainer:
    """
    멤브레인 시스템 모델 훈련 및 시뮬레이션을 위한 클래스
    """
    def __init__(self, model, device='cuda' if torch.cuda.is_available() else 'cpu'):
        self.model = model.to(device)
        self.device = device
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=0.001)
        # 몰 변화량과 상태 예측을 위한 손실 함수
        self.mol_change_criterion = nn.MSELoss()
        self.state_criterion = nn.MSELoss()
    
    def train(self, train_loader, val_loader=None, epochs=100, mol_change_weight=1.0, state_weight=1.0):
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
        
        for epoch in range(epochs):
            epoch_loss = 0.0
            
            for batch_idx, (features, mol_change_targets, state_targets) in enumerate(train_loader):
                features = features.to(self.device)
                mol_change_targets = mol_change_targets.to(self.device)
                state_targets = state_targets.to(self.device)
                
                # 그래디언트 초기화
                self.optimizer.zero_grad()
                
                # 순전파
                mol_change_predictions, state_predictions = self.model(features)
                
                # 손실 계산 - 모든 시간 단계에 대한 손실
                mol_change_loss = self.mol_change_criterion(mol_change_predictions, mol_change_targets)
                state_loss = self.state_criterion(state_predictions, state_targets)
                total_loss = mol_change_weight * mol_change_loss + state_weight * state_loss
                
                # 역전파 및 최적화
                total_loss.backward()
                self.optimizer.step()
                
                epoch_loss += total_loss.item()
            
            # 에폭 평균 손실
            avg_train_loss = epoch_loss / len(train_loader)
            train_losses.append(avg_train_loss)
            
            # 검증 손실 계산 (있는 경우)
            if val_loader is not None:
                self.model.eval()
                val_loss = 0.0
                
                with torch.no_grad():
                    for features, mol_change_targets, state_targets in val_loader:
                        features = features.to(self.device)
                        mol_change_targets = mol_change_targets.to(self.device)
                        state_targets = state_targets.to(self.device)
                        
                        mol_change_predictions, state_predictions = self.model(features)
                        
                        mol_change_loss = self.mol_change_criterion(mol_change_predictions, mol_change_targets)
                        state_loss = self.state_criterion(state_predictions, state_targets)
                        total_loss = mol_change_weight * mol_change_loss + state_weight * state_loss
                        
                        val_loss += total_loss.item()
                
                avg_val_loss = val_loss / len(val_loader)
                val_losses.append(avg_val_loss)
                
                print(f'Epoch {epoch+1}/{epochs}, Train Loss: {avg_train_loss:.6f}, Val Loss: {avg_val_loss:.6f}')
                
                self.model.train()
            else:
                print(f'Epoch {epoch+1}/{epochs}, Train Loss: {avg_train_loss:.6f}')
        
        return train_losses, val_losses
    
    def evaluate_experiment(self, exp_data, initial_conditions, operation_params):
        """
        특정 실험 데이터와 모델 예측 결과를 비교 평가
        
        입력:
        - exp_data: 단일 실험 데이터프레임
        - initial_conditions: 초기 조건
        - operation_params: 운전 조건
        
        출력:
        - 실제 데이터와 모델 예측의 비교 결과
        """
        self.model.eval()
        
        # 실험 시간 포인트
        time_points = exp_data['t'].values
        max_time = time_points[-1]
        
        # 0.1시간 간격으로 시뮬레이션
        sim_time_points = np.arange(0, max_time + self.model.time_step, self.model.time_step)
        num_steps = len(sim_time_points) - 1
        
        # 시뮬레이션 실행
        sim_results = self.simulate(initial_conditions, operation_params, num_steps)
        
        # 실제 데이터에 가장 가까운 시뮬레이션 시간점 찾기
        comparison = {}
        for i, t in enumerate(time_points):
            # 가장 가까운 시뮬레이션 시간 인덱스 찾기
            closest_idx = np.argmin(np.abs(sim_time_points - t))
            
            # 실제 값과 예측 값 비교
            comparison[t] = {
                'actual': {
                    'CF_LA': exp_data.iloc[i]['CF_LA'],
                    'CA_LA': exp_data.iloc[i]['CA_LA'],
                    'CF_K': exp_data.iloc[i]['CF_K'],
                    'CB_K': exp_data.iloc[i]['CB_K'],
                    'VF': exp_data.iloc[i]['VF'],
                    'VA': exp_data.iloc[i]['VA'],
                    'VB': exp_data.iloc[i]['VB']
                },
                'predicted': {
                    'CF_LA': sim_results['feed_la'][closest_idx],
                    'CA_LA': sim_results['acid_la'][closest_idx],
                    'CF_K': sim_results['feed_k'][closest_idx],
                    'CB_K': sim_results['base_k'][closest_idx],
                    'VF': sim_results['feed_volume'][closest_idx],
                    'VA': sim_results['acid_volume'][closest_idx],
                    'VB': sim_results['base_volume'][closest_idx]
                }
            }
            
            # 백분율 오차 계산
            error = {}
            for key in comparison[t]['actual'].keys():
                actual = comparison[t]['actual'][key]
                predicted = comparison[t]['predicted'][key]
                if actual != 0:
                    error[key] = abs((predicted - actual) / actual) * 100
                else:
                    error[key] = float('inf') if predicted != 0 else 0
            
            comparison[t]['error'] = error
        
        return {
            'comparison': comparison,
            'simulation': sim_results,
            'exp_time_points': time_points,
            'sim_time_points': sim_time_points
        }
    
    def simulate(self, initial_conditions, operation_params, num_steps):
        """
        시스템 시뮬레이션 실행
        
        입력:
        - initial_conditions: 초기 농도 및 부피 조건 (딕셔너리)
        - operation_params: 운전 조건 (온도, 전압, 전해질 농도) (딕셔너리)
        - num_steps: 시뮬레이션 단계 수
        
        출력:
        - 시뮬레이션 결과 (딕셔너리)
        """
        self.model.eval()
        
        # 초기 조건 설정
        feed_la_init = initial_conditions['feed_la']
        feed_k_init = initial_conditions['feed_k']
        acid_la_init = initial_conditions.get('acid_la', 0)
        acid_k_init = 0  # Acid에는 K+가 없음
        base_la_init = 0  # Base에는 LA가 없음
        base_k_init = initial_conditions.get('base_k', 0)
        
        feed_vol_init = initial_conditions['feed_volume']
        acid_vol_init = initial_conditions['acid_volume']
        base_vol_init = initial_conditions['base_volume']
        
        # 운전 조건
        temperature = operation_params['temperature']
        voltage = operation_params['voltage']
        electrolyte = operation_params['electrolyte']
        ci = operation_params['ci']  # 초기 LA 농도
        ki = operation_params['ki']  # 초기 K+ 농도
        
        # 결과 저장 배열
        results = {
            'time': np.arange(num_steps + 1) * self.model.time_step,
            'feed_la': np.zeros(num_steps + 1),
            'feed_k': np.zeros(num_steps + 1),
            'acid_la': np.zeros(num_steps + 1),
            'base_k': np.zeros(num_steps + 1),
            'feed_volume': np.zeros(num_steps + 1),
            'acid_volume': np.zeros(num_steps + 1),
            'base_volume': np.zeros(num_steps + 1),
            'la_mol_change': np.zeros(num_steps),
            'k_mol_change': np.zeros(num_steps),
            'water_acid_change': np.zeros(num_steps),
            'water_base_change': np.zeros(num_steps)
        }
        
        # 초기 값 설정
        results['feed_la'][0] = feed_la_init
        results['feed_k'][0] = feed_k_init
        results['acid_la'][0] = acid_la_init
        results['base_k'][0] = base_k_init
        results['feed_volume'][0] = feed_vol_init
        results['acid_volume'][0] = acid_vol_init
        results['base_volume'][0] = base_vol_init
        
        # 시계열 입력을 위한 버퍼 초기화 (sequence_length 크기의 슬라이딩 윈도우)
        sequence_buffer = []
        for _ in range(self.model.sequence_length):
            sequence_buffer.append([
                temperature, voltage, electrolyte, ci, ki,
                feed_vol_init, acid_vol_init, base_vol_init,
                feed_la_init, feed_k_init, acid_la_init, base_k_init
            ])
        
        with torch.no_grad():
            # 시뮬레이션 실행
            for i in range(num_steps):
                # NumPy 배열을 PyTorch 텐서로 변환
                input_tensor = torch.FloatTensor([sequence_buffer]).to(self.device)  # [batch_size=1, seq_len, features]
                
                # 예측 실행
                mol_change_pred, state_pred = self.model(input_tensor)
                
                # 마지막 시간 단계에 대한 예측만 사용
                mol_change_np = mol_change_pred.squeeze()[-1].cpu().numpy() 
                state_np = state_pred.squeeze()[-1].cpu().numpy()
                
                # 결과 저장
                results['la_mol_change'][i] = mol_change_np[0]
                results['k_mol_change'][i] = mol_change_np[1]
                results['water_acid_change'][i] = mol_change_np[2]
                results['water_base_change'][i] = mol_change_np[3]
                
                # 다음 시간 단계의 상태 저장
                results['feed_la'][i+1] = state_np[0]
                results['feed_k'][i+1] = state_np[1]
                results['acid_la'][i+1] = state_np[2]
                # Acid K+ 및 Base LA는 항상 0
                results['base_k'][i+1] = state_np[5]
                results['feed_volume'][i+1] = state_np[6]
                results['acid_volume'][i+1] = state_np[7]
                results['base_volume'][i+1] = state_np[8]
                
                # 시퀀스 버퍼 업데이트 (가장 오래된 항목 제거, 새 항목 추가 - 슬라이딩 윈도우)
                sequence_buffer.pop(0)
                sequence_buffer.append([
                    temperature, voltage, electrolyte, ci, ki,
                    results['feed_volume'][i+1], results['acid_volume'][i+1], results['base_volume'][i+1],
                    results['feed_la'][i+1], results['feed_k'][i+1],
                    results['acid_la'][i+1], results['base_k'][i+1]
                ])
        
        return results
    
    def plot_parity_and_calculate_r2(self, eval_results):
        """
        실험값과 계산값의 R2 값을 계산
        """
        from sklearn.metrics import r2_score
        
        # 모든 실험의 데이터 수집
        actual_values = {
            'CF_LA': [], 'CA_LA': [], 'CF_K': [], 'CB_K': [],
            'VF': [], 'VA': [], 'VB': []
        }
        predicted_values = {
            'CF_LA': [], 'CA_LA': [], 'CF_K': [], 'CB_K': [],
            'VF': [], 'VA': [], 'VB': []
        }
        
        # 데이터 수집
        for exp_results in eval_results.values():
            comparison = exp_results['comparison']
            for t in comparison:
                for key in actual_values:
                    actual_values[key].append(comparison[t]['actual'][key])
                    predicted_values[key].append(comparison[t]['predicted'][key])
        
        # 농도 변수
        conc_vars = ['CF_LA', 'CA_LA', 'CF_K', 'CB_K']
        # 부피 변수
        vol_vars = ['VF', 'VA', 'VB']
        
        # 전체 R2 계산 및 출력
        print("\n=== Overall R² Scores ===")
        total_r2 = 0
        n_vars = 0
        
        print("\nConcentration R² Scores:")
        conc_r2 = 0
        for var in conc_vars:
            r2 = r2_score(actual_values[var], predicted_values[var])
            print(f"{var}: {r2:.4f}")
            conc_r2 += r2
            total_r2 += r2
            n_vars += 1
        print(f"Average Concentration R²: {conc_r2/len(conc_vars):.4f}")
        
        print("\nVolume R² Scores:")
        vol_r2 = 0
        for var in vol_vars:
            r2 = r2_score(actual_values[var], predicted_values[var])
            print(f"{var}: {r2:.4f}")
            vol_r2 += r2
            total_r2 += r2
            n_vars += 1
        print(f"Average Volume R²: {vol_r2/len(vol_vars):.4f}")
        
        print(f"\nOverall Average R²: {total_r2/n_vars:.4f}")
        
        return {
            'concentration_r2': conc_r2/len(conc_vars),
            'volume_r2': vol_r2/len(vol_vars),
            'overall_r2': total_r2/n_vars
        }

def evaluate_model_on_experiments(csv_path, trainer, output_dir=None):
    """
    모든 실험에 대해 모델 평가 및 결과 시각화
    """
    df = pd.read_csv(csv_path)
    results = {}
    
    # 각 실험별로 평가
    for exp in df['exp'].unique():
        exp_data = df[df['exp'] == exp].sort_values('t')
        
        # 첫 번째 행 데이터로 초기 조건 설정
        first_row = exp_data.iloc[0]
        
        initial_conditions = {
            'feed_la': first_row['CF_LA'],
            'feed_k': first_row['CF_K'],
            'acid_la': first_row['CA_LA'],
            'base_k': first_row['CB_K'],
            'feed_volume': first_row['VF'],
            'acid_volume': first_row['VA'],
            'base_volume': first_row['VB']
        }
        
        operation_params = {
            'temperature': first_row['T'],
            'voltage': first_row['V'],
            'electrolyte': first_row['E'],
            'ci': first_row['Ci'],
            'ki': first_row['Ki']
        }
        
        # 실험 평가
        eval_result = trainer.evaluate_experiment(exp_data, initial_conditions, operation_params)
        results[exp] = eval_result
    
    # 결과 시각화 및 R2 계산
    r2_scores = trainer.plot_parity_and_calculate_r2(results)
    
    return results, r2_scores

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
    

def data_loaders_from_csv(csv_path, sequence_length=1, batch_size=32, train_ratio=0.8):
    """
    CSV 파일에서 훈련 및 검증 데이터 로더 생성
    """
    train_dataset = BMEDDataset(csv_path, sequence_length=sequence_length, train=True, train_ratio=train_ratio)
    test_dataset = BMEDDataset(csv_path, sequence_length=sequence_length, train=False, train_ratio=train_ratio)
    
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

def train_model_from_csv(csv_path, lstm_units=64, sequence_length=1, time_step=0.1, 
                        batch_size=32, epochs=100, train_ratio=0.8,
                        mol_change_weight=1.0, state_weight=1.0):
    """
    CSV 파일에서 모델 훈련을 위한 편의 함수
    """
    # 데이터 로더 생성
    train_loader, test_loader = data_loaders_from_csv(
        csv_path, sequence_length, batch_size, train_ratio
    )
    
    # 모델 생성
    model = MembraneSystemModel(lstm_units=lstm_units, time_step=time_step, sequence_length=sequence_length)
    trainer = MembraneSystemTrainer(model)
    
    # 모델 훈련
    train_losses, val_losses = trainer.train(
        train_loader, 
        test_loader, 
        epochs=epochs,
        mol_change_weight=mol_change_weight,
        state_weight=state_weight
    )
    
    return model, trainer, train_losses, val_losses

if __name__ == "__main__":
    # CSV 파일 경로
    csv_path = 'BMED_data_v5.csv'
    
    # 모델 훈련
    model, trainer, train_losses, val_losses = train_model_from_csv(
        csv_path=csv_path,
        lstm_units=64,
        sequence_length=3,  # 이전 3개 시점 데이터 사용
        time_step=0.01,
        batch_size=16,
        epochs=200,
        train_ratio=0.8,
        mol_change_weight=1.0,
        state_weight=1.0
    )
    
    # 학습 곡선 시각화
    plt.figure(figsize=(10, 5))
    plt.plot(train_losses, label='Train Loss')
    plt.plot(val_losses, label='Validation Loss')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.title('Training and Validation Loss')
    plt.legend()
    plt.grid(True)
    plt.savefig('training_curve.png')
    plt.close()
    
    # 모델 저장
    torch.save(model.state_dict(), 'membrane_model.pth')
    
    # 모든 실험에 대한 평가 및 결과 시각화
    results, r2_scores = evaluate_model_on_experiments(csv_path, trainer)
    
    # 결과 출력
    print("\n=== Final Results ===")
    print(f"Average Concentration R²: {r2_scores['concentration_r2']:.4f}")
    print(f"Average Volume R²: {r2_scores['volume_r2']:.4f}")
    print(f"Overall Average R²: {r2_scores['overall_r2']:.4f}")
