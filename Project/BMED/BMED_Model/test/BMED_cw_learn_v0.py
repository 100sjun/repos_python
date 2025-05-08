# Load modules
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader, Subset
from sklearn.preprocessing import StandardScaler
import pandas as pd
import json

# Fix the random seed
np.random.seed(42)
torch.manual_seed(42)

# Feedforward network for migration prediction
class MigrationPredictor(nn.Module):
    def __init__(self, hidden_nodes = 64, hidden_layers = 3, dropout = 0.2):
        super().__init__()
        
        n_features = 7
        n_outputs = 4

        # Layer configuration
        layers = []
        # input layer
        layers.append(nn.Linear(n_features, hidden_nodes))
        layers.append(nn.ReLU())
        layers.append(nn.Dropout(dropout))

        # hidden layers
        for _ in range(hidden_layers - 1):
            layers.append(nn.Linear(hidden_nodes, hidden_nodes))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(dropout))

        # output layer
        layers.append(nn.Linear(hidden_nodes, n_outputs))

        self.model = nn.Sequential(*layers)

    def forward(self, x):
        return self.model(x)
    
# Physical layers for state update
class PhysicalLayer:
    def __init__(self, dt = 0.1):
        self.dt = dt
    
    def update_state(self, cur_states, migrations):
        # Current States
        T = cur_states[0]
        V = cur_states[1]
        E = cur_states[2]
        CF_LA = cur_states[3]
        CA_LA = cur_states[4]
        CF_K = cur_states[5]
        CB_K = cur_states[6]
        VF = cur_states[7]
        VA = cur_states[8]
        VB = cur_states[9]
        
        # Migration
        dNLA = migrations[0] * self.dt
        dNK = migrations[1] * self.dt
        dVA = migrations[2] * self.dt
        dVB = migrations[3] * self.dt

        # Fixed variables
        nT = T
        nV = V
        nE = E     
        
        # New Volumes
        nVF = VF - dVA - dVB
        nVA = VA + dVA
        nVB = VB + dVB

        # New Concentrations
        nCF_LA = (CF_LA * VF - dNLA) / nVF
        nCA_LA = (CA_LA * VA + dNLA) / nVA
        nCF_K = (CF_K * VF - dNK) / nVF
        nCB_K = (CB_K * VB + dNK) / nVB

        # New States
        new_states = cur_states.clone()
        new_states[0] = nT
        new_states[1] = nV
        new_states[2] = nE
        new_states[3] = nCF_LA
        new_states[4] = nCA_LA
        new_states[5] = nCF_K
        new_states[6] = nCB_K
        new_states[7] = nVF
        new_states[8] = nVA
        new_states[9] = nVB
        
        return new_states
    
# Generate dataset
class BMEDDataset(Dataset):
    def __init__(self, dict_spline):
        self.states = ['T', 'V', 'E', 'CF_LA', 'CA_LA', 'CF_K', 'CB_K', 'VF', 'VA', 'VB']
        self.experiments = []

        for exp_id, exp_data in dict_spline.items():
            # Save whole data of each experiment in one sample
            exp_array = exp_data[self.states].values
            times = exp_data['t'].values
            self.experiments.append({
                'init_state': torch.tensor(exp_array[0], dtype = torch.float32), # initial state
                'measured_state': torch.tensor(exp_array, dtype = torch.float32), # whole measurements
                'times': torch.tensor(times, dtype = torch.float32) # time points
            })
    
    def __len__(self):
        return len(self.experiments)
    
    def __getitem__(self, idx):
        return self.experiments[idx]
    
# Physics-informed neural network
class BMEDModel(nn.Module):
    def __init__(self, hidden_nodes = 192, hidden_layers = 3, dt = 0.1, scaler = None, dropout = 0.2):
        super().__init__()
        self.migration_predictor = MigrationPredictor(hidden_nodes, hidden_layers, dropout = dropout)
        self.physical_layer = PhysicalLayer(dt)
        self.scaler = scaler
        self.dt = dt

    def forward(self, init_state, times):

        cur_state = init_state # batch size 1
        cur_time = 0.0
        pred_states = []
        measured_indices = []

        times = times
        times_np = times[0].numpy()
        max_time = times_np[-1]
        measured_indices.append(0)

        # 초기 상태 저장
        pred_states.append(cur_state)

        while cur_time < max_time:
            # input_feature에 해당하는 변수만 정규화
            input_state = cur_state[:, :7]  # input feature 추출, 2차원 유지지
            
            norm_input = self.scaler.transform(input_state.detach().numpy())
            norm_input = torch.tensor(norm_input)
            
            # 상태 예측
            migration = self.migration_predictor(norm_input)  # (1, 6) -> (1, 3)
            cur_state = self.physical_layer.update_state(cur_state[0], migration[0]).unsqueeze(0)  # (1,8)
            pred_states.append(cur_state)  # (1, 8)
            cur_time += self.dt

            # 측정 시간과 매칭
            for t in times_np:
                if abs(cur_time - t) < self.dt/2:
                    measured_indices.append(len(pred_states) - 1)

        # 현재 배치의 예측 상태들을 스택
        pred_states = torch.cat(pred_states, dim=0)  # (n_steps, 8)

        return pred_states, measured_indices
    
# Declare custom loss
def custom_loss(pred_states, measured_indices, measured_states):

    # default weight
    default_wt = {
        'T': 0,
        'V': 0,
        'E': 0,
        'CF_LA': 1,
        'CA_LA': 1,
        'CF_K': 0.1,
        'CB_K': 0.1,
        'VF': 0.5,
        'VA': 1,
        'VB': 0.1
    }

    wt_tensor = torch.tensor([
        default_wt['T'], default_wt['V'], default_wt['E'], 
        default_wt['CF_LA'], default_wt['CA_LA'], default_wt['CF_K'], default_wt['CB_K'],
        default_wt['VF'], default_wt['VA'], default_wt['VB']
    ])

    total_loss = 0
    for idx, measured_state in zip(measured_indices, measured_states[0]):
        predicted_state = pred_states[idx]
        sq_errors = (predicted_state - measured_state) ** 2
        wt_errors = sq_errors * wt_tensor

        total_loss += torch.mean(wt_errors)

    return total_loss

# Evaluate the model
def evaluate_model(model, val_loader):
    model.eval()
    total_loss = 0
    num_exp = 0

    with torch.no_grad():
        for exp in val_loader:
            init_state = exp['init_state']
            measured_states = exp['measured_state']
            times = exp['times']

            pred_states, measured_indices = model(init_state, times)
            

            loss = custom_loss(pred_states, measured_indices, measured_states)
            total_loss += loss.item()
            num_exp += 1

    avg_loss = total_loss / num_exp

    return avg_loss

# Train the model
def train_model(model, train_loader, val_loader,epochs = 100, learning_rate = 0.001, weight_decay = 1e-5, clip_norm = 0.1, schd_factor = 0.5):
    optimizer = torch.optim.Adam(model.parameters(), lr = learning_rate, weight_decay= weight_decay)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, patience = 5, factor = schd_factor, min_lr = 1e-6
    )

    patience = 50
    min_delta = 0.0001

    best_val_loss = float('inf')
    patience_counter = 0
    best_model_state = None
    stopped_epoch = epochs

    for epoch in range(epochs):
        model.train()
        train_loss = 0
        num_batches = 0

        for exp in train_loader:
            optimizer.zero_grad()

            init_state = exp['init_state']
            measured_state = exp['measured_state']
            times = exp['times']

            # Simulation
            pred_state, measured_indices = model(init_state, times)

            # Loss
            loss = custom_loss(pred_state, measured_indices, measured_state)
            train_loss += loss.item()
            num_batches += 1

            loss.backward()

            # gradient clipping
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm = clip_norm)
            optimizer.step()

        avg_train_loss = train_loss / num_batches
        val_loss = evaluate_model(model, val_loader)
        
        # 학습률 스케줄러 업데이트
        scheduler.step(val_loss)
        current_lr = optimizer.param_groups[0]['lr']

        # 매 에포크마다 검증 손실이 개선될 때마다 모델 저장
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            patience_counter = 0
            best_model_state = model.state_dict().copy()
            # 모델 저장
            torch.save(model.state_dict(), f'best_model_epoch.pt')
            print(f'에포크 {epoch+1}/{epochs} - 훈련 손실: {avg_train_loss:.6f}, 검증 손실: {val_loss:.6f}, 학습률: {current_lr:.8f}, 인내심 카운터: {patience_counter}')
            print(f'검증 손실 개선! 최고 검증 손실: {best_val_loss:.6f}, 모델 저장됨')
        else:
            patience_counter += 1
            print(f'에포크 {epoch+1}/{epochs} - 훈련 손실: {avg_train_loss:.6f}, 검증 손실: {val_loss:.6f}, 학습률: {current_lr:.8f}, 인내심 카운터: {patience_counter}')

        # Early stopping 조건 추가
        if patience_counter >= patience:
            print(f'조기 종료! {patience}번의 에포크 동안 검증 손실이 개선되지 않았습니다.')
            stopped_epoch = epoch + 1
            # 최적의 모델 상태로 복원
            model.load_state_dict(best_model_state)
            break

    # 학습이 완료된 후 최적의 모델 상태로 복원
    if best_model_state is not None:
        model.load_state_dict(best_model_state)
    
    print(f'훈련 완료! 최종 에포크: {stopped_epoch}, 최고 검증 손실: {best_val_loss:.6f}')
    return stopped_epoch


# Load raw data
df = pd.read_csv('BMED_data_v8.csv')

# Split the data by experiment number
dict_spline = {}
for exp in df['exp'].unique():
    dict_spline[exp] = df[df['exp'] == exp].sort_values('t')

# Declare scaler
scaler = StandardScaler()
col_to_scale = ['T', 'V', 'E', 'CF_LA', 'CA_LA', 'CF_K', 'CB_K']
scaler.fit(df[col_to_scale].values)

# Load best hyperparameter set from json
with open('hpOpt_checkpoint_v1.json', 'r') as f:
    opt_hp_para = json.load(f)

hidden_nodes = opt_hp_para['best_params']['hidden_nodes']
hidden_layers = opt_hp_para['best_params']['hidden_layers']
lr = opt_hp_para['best_params']['lr']
wd = opt_hp_para['best_params']['weight_decay']
epochs = 10000
dropout = opt_hp_para['best_params']['dropout']/100
clip_norm = opt_hp_para['best_params']['clip_norm']/100
schd_factor = opt_hp_para['best_params']['schd_factor']/100
dt = 0.1

# Generate dataset
dataset = BMEDDataset(dict_spline)
dataset[0]

# Split dataset into train and test
dataset_size = len(dataset)

train_size = int(dataset_size * 0.8)
val_size = int(dataset_size * 0.15)
test_size = dataset_size - train_size - val_size

idx = list(range(dataset_size))
np.random.shuffle(idx)

train_idx = idx[:train_size]
val_idx = idx[train_size:train_size + val_size]
test_idx = idx[train_size + val_size:]

# Generate subsets
train_dataset = Subset(dataset, train_idx)
val_dataset = Subset(dataset, val_idx)
test_dataset = Subset(dataset, test_idx)

# Generate dataloader
train_loader = DataLoader(train_dataset, batch_size = 1, shuffle = True)
val_loader = DataLoader(val_dataset, batch_size = 1, shuffle = False)
test_loader = DataLoader(test_dataset, batch_size = 1, shuffle = False)

# Generate model
model = BMEDModel(hidden_nodes = hidden_nodes, hidden_layers = hidden_layers, dt = dt, scaler = scaler, dropout = dropout)

# 모델 훈련
train_model(model, train_loader, val_loader, epochs = epochs, learning_rate = lr, weight_decay = wd, clip_norm = clip_norm, schd_factor = schd_factor)

