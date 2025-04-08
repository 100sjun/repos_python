# Load Modules
import torch
import torch.nn as nn
import pandas as pd
from torch.utils.data import Dataset, DataLoader
import numpy as np

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
        VF = features[..., 3:4] # 현재 상태 Feed 부피 [L]
        VA = features[..., 4:5] # 현재 상태 Acid 부피 [L]
        VB = features[..., 5:6] # 현재 상태 Base 부피 [L]
        CFLA = features[..., 6:7] # 현재 상태 Feed LA 농도 [mol/L]
        CFK = features[..., 7:8] # 현재 상태 Feed K+ 농도 [mol/L]
        CALA = features[..., 8:9] # 현재 상태 Acid LA 농도 [mol/L]
        CBK = features[..., 9:10] # 현재 상태 Base K+ 농도 [mol/L]
        
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
            input_size=10, # [T, V, E, VF, VA, VB, CF_LA, CF_K, CA_LA, CB_K]
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
    
class BMEDDataset(Dataset):
    '''BMED CSV 데이터를 위한 PyTorch 데이터셋'''
    def __init__(self, csv_path, sequence_length=10, dt=0.1, train=True, train_ratio=0.8):
        '''CSV 파일 로드'''
        self.df = pd.read_csv(csv_path)
        self.sequence_length = sequence_length
        self.dt = dt
        
        '''실험별 데이터 분리'''
        self.experiment_data = {}
        for exp in self.df['exp'].unique():
            self.experiment_data[exp] = self.df[self.df['exp'] == exp].sort_values('t')

        '''훈련/테스트 셋 분리'''
        all_exps = list(self.experiment_data.keys())
        np.random.seed(42)
        np.random.shuffle(all_exps)
        np.random.seed(None)

        split_idx = int(len(all_exps) * train_ratio)
        self.exps_to_use = all_exps[:split_idx] if train else all_exps[split_idx:]

        '''보간된 데이터 준비 및 인덱스 구성'''
        self.prepare_interpolated_data()

    def linear_interpolation(self, exp_data):
        t_exp = exp_data['t'].values
        t_interp = np.arange(t_exp[0], t_exp[-1], self.dt)
        interpolated_data = {}

        '''운전조건'''
        for col in ['T','V','E']:
            interpolated_data[col] = np.interp(t_interp, t_exp, exp_data[col].values)
        
        '''몰 농도'''
        for col in ['CF_LA','CF_K','CA_LA','CB_K']:
            interpolated_data[col] = np.interp(t_interp, t_exp, exp_data[col].values)
        
        '''부피'''
        for col in ['VF','VA','VB']:
            interpolated_data[col] = np.interp(t_interp, t_exp, exp_data[col].values)

        return t_interp, interpolated_data

    def prepare_interpolated_data(self):
        '''보간된 데이터 준비 및 인덱스 구성'''
        self.interpolated_data = {}
        self.exp_indices = []  # 이것만 유지

        for exp_id in self.exps_to_use:
            exp_data = self.experiment_data[exp_id]

            if len(exp_data) < 2:
                continue

            '''보간 실시'''
            t_interp, interpolated_data = self.linear_interpolation(exp_data)
            self.interpolated_data[exp_id] = {
                't': t_interp,
                **interpolated_data
            }
            
            '''시퀀스 생성을 위한 인덱스 계산: 슬라이딩 윈도우'''
            sequences = []
            for i in range(len(t_interp) - self.sequence_length + 1):
                sequence_indices = list(range(i, i + self.sequence_length))
                sequences.append(sequence_indices)

            if sequences:  # 시퀀스가 하나 이상 있는 경우
                self.exp_indices.append((exp_id, sequences))

    def __len__(self):
        return sum(len(sequences) for _, sequences in self.exp_indices)

    def __getitem__(self, idx):
        '''idx에 해당하는 실험과 sequence 찾기'''
        current_pos = 0
        for exp_id, sequence_list in self.exp_indices:
            if current_pos + len(sequence_list) > idx:
                '''현재 실험에서 해당 sequence 찾기'''
                sequence = sequence_list[idx - current_pos]
                
                '''현재 상태 특성 벡터 구성'''
                features = torch.zeros((self.sequence_length, 10))  # 10 = 특성 수
                for i, t in enumerate(sequence):
                    features[i] = torch.tensor([
                        self.interpolated_data[exp_id]['T'][t],
                        self.interpolated_data[exp_id]['V'][t],
                        self.interpolated_data[exp_id]['E'][t],
                        self.interpolated_data[exp_id]['VF'][t],
                        self.interpolated_data[exp_id]['VA'][t],
                        self.interpolated_data[exp_id]['VB'][t],
                        self.interpolated_data[exp_id]['CF_LA'][t],
                        self.interpolated_data[exp_id]['CF_K'][t],
                        self.interpolated_data[exp_id]['CA_LA'][t],
                        self.interpolated_data[exp_id]['CB_K'][t]
                    ])

                '''몰 변화량 계산'''
                mol_changes = torch.zeros((self.sequence_length, 4))  # 4 = 변화량 수
                for i in range(self.sequence_length - 1):
                    current_t = sequence[i]
                    next_t = sequence[i + 1]
                    
                    # LA 및 K+ 몰수 계산 (현재)
                    current_la_acid = (self.interpolated_data[exp_id]['CA_LA'][current_t] * 
                                     self.interpolated_data[exp_id]['VA'][current_t])
                    current_k_base = (self.interpolated_data[exp_id]['CB_K'][current_t] * 
                                    self.interpolated_data[exp_id]['VB'][current_t])
                    
                    # LA 및 K+ 몰수 계산 (다음)
                    next_la_acid = (self.interpolated_data[exp_id]['CA_LA'][next_t] * 
                                  self.interpolated_data[exp_id]['VA'][next_t])
                    next_k_base = (self.interpolated_data[exp_id]['CB_K'][next_t] * 
                                 self.interpolated_data[exp_id]['VB'][next_t])
                    
                    # 몰 변화량 계산
                    mol_changes[i] = torch.tensor([
                        (next_la_acid - current_la_acid) / self.dt,  # LA 몰 변화율
                        (next_k_base - current_k_base) / self.dt,    # K+ 몰 변화율
                        (self.interpolated_data[exp_id]['VA'][next_t] - 
                         self.interpolated_data[exp_id]['VA'][current_t]) / self.dt,  # Acid 부피 변화율
                        (self.interpolated_data[exp_id]['VB'][next_t] - 
                         self.interpolated_data[exp_id]['VB'][current_t]) / self.dt   # Base 부피 변화율
                    ])
                
                # 마지막 시점의 변화량은 이전 변화량과 동일하게 설정
                mol_changes[-1] = mol_changes[-2]

                '''다음 상태 구성'''
                next_states = torch.zeros((self.sequence_length, 10))  # 10 = 상태 변수 수
                for i, t in enumerate(sequence):
                    if t + 1 < len(self.interpolated_data[exp_id]['t']):
                        next_t = t + 1
                        next_states[i] = torch.tensor([
                            self.interpolated_data[exp_id]['T'][next_t],
                            self.interpolated_data[exp_id]['V'][next_t],
                            self.interpolated_data[exp_id]['E'][next_t],
                            self.interpolated_data[exp_id]['CF_LA'][next_t],
                            self.interpolated_data[exp_id]['CF_K'][next_t],
                            self.interpolated_data[exp_id]['CA_LA'][next_t],
                            self.interpolated_data[exp_id]['CB_K'][next_t],
                            self.interpolated_data[exp_id]['VF'][next_t],
                            self.interpolated_data[exp_id]['VA'][next_t],
                            self.interpolated_data[exp_id]['VB'][next_t]
                        ])
                    else:
                        # 마지막 시점은 현재 상태 유지
                        next_states[i] = next_states[i-1]

                return features, mol_changes, next_states
            
            current_pos += len(sequence_list)
        
        raise IndexError("Index out of range")
                
                
                
            
