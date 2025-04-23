# Load Modules
import pandas as pd
from sklearn.preprocessing import StandardScaler
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader, Subset
import matplotlib.pyplot as plt
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score
import optuna
from optuna.samplers import GPSampler
import pickle
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

# Generate PINN
class BMEDModel(nn.Module):
    def __init__(self, hidden_nodes = 32, hidden_layers = 3, dt = 0.1, scaler = None, dropout = 0.2):
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

        for exp in train_loader:
            optimizer.zero_grad()

            init_state = exp['init_state']
            measured_state = exp['measured_state']
            times = exp['times']

            # Simulation
            pred_state, measured_indices = model(init_state, times)

            # Loss
            loss = custom_loss(pred_state, measured_indices, measured_state)

            loss.backward()

            # gradient clipping
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm = clip_norm)
            optimizer.step()

        val_loss = evaluate_model(model, val_loader)

        if val_loss < best_val_loss - min_delta:
            best_val_loss = val_loss
            patience_counter = 0
            best_model_state = model.state_dict().copy()
        else:
            patience_counter += 1

        # Early stopping 조건 추가
        if patience_counter >= patience:
            stopped_epoch = epoch + 1
            # 최적의 모델 상태로 복원
            model.load_state_dict(best_model_state)
            break

                # 학습이 완료된 후 최적의 모델 상태로 복원
    if best_model_state is not None:
        model.load_state_dict(best_model_state)

    return stopped_epoch

# Data preparation
def prepare_data(df_path, exp_idx = None):
    # Load raw data
    df = pd.read_csv(df_path)

    if exp_idx is not None:
        df = df[df['exp'].isin(exp_idx)]

    # split the data by experiment number
    dict_spline = {}
    for exp in df['exp'].unique():
        dict_spline[exp] = df[df['exp'] == exp].sort_values('t')

    # scaler
    scaler = StandardScaler()
    col_to_scale = ['T', 'V', 'E', 'CF_LA', 'CA_LA', 'CF_K', 'CB_K']
    scaler.fit(df[col_to_scale].values)
    return dict_spline, scaler

# Save data during calculation
def save_callback(study, trial):
    with open('hpOpt_study.pkl', 'wb') as f:
        pickle.dump(study, f)

    # 2. best value만 json으로 저장
    best_trial_info = {
        'best_params': study.best_trial.params,
        'best_value': study.best_trial.value,
        'best_avg_epoch': study.best_trial.user_attrs['avg_epoch']
    }

    with open('hpOpt_checkpoint.json', 'w') as f:
        json.dump(best_trial_info, f, indent=4)

# Optuna objective function
def objective(trial):
    # load data
    df_path = './data/BMED_data_v8.csv'
    dict_spline, scaler = prepare_data(df_path)

    # Hyperparameter 설정
    hidden_nodes = trial.suggest_int('hidden_nodes', 16, 256, step = 16)
    hidden_layers = trial.suggest_int('hidden_layers', 1, 10, step = 1)
    lr = trial.suggest_float('lr', 1e-6, 1e-2, log = True)
    wd = trial.suggest_float('weight_decay', 1e-6, 1e-3, log = True)
    epochs = trial.suggest_int('epochs', 100, 1000, step = 100)
    dropout = trial.suggest_int('dropout', 20, 50, step = 10) / 100
    clip_norm = trial.suggest_int('clip_norm', 10, 100, step = 10) / 100
    schd_factor = trial.suggest_int('schd_factor', 10, 50, step = 10) / 100
    dt = 0.1

    # 5-fold 설정
    kfolds = 5
    kf = KFold(n_splits=kfolds, shuffle=True, random_state=42)

    # dataset 생성
    dataset = BMEDDataset(dict_spline)

    # k-fold 교차 검증 스코어 저장
    fold_scores = []
    fold_epochs = []

    for fold, (train_idx, val_idx) in enumerate(kf.split(dataset)):
        # kfold로 전체 dataset에 대해서 train과 validation subset을 생성
        train_subset = Subset(dataset, train_idx)
        val_subset = Subset(dataset, val_idx)

        # Dataloader 생성
        train_loader = DataLoader(train_subset, batch_size = 1, shuffle = True)
        val_loader = DataLoader(val_subset, batch_size = 1, shuffle = False)

        # 모델 생성
        model = BMEDModel(hidden_nodes=hidden_nodes, hidden_layers=hidden_layers, dt = dt, scaler=scaler, dropout=dropout)

        # 모델 학습
        stopped_epoch = train_model(
            model = model, train_loader = train_loader, val_loader=val_loader, 
            epochs = epochs, learning_rate = lr, weight_decay = wd, clip_norm = clip_norm, schd_factor = schd_factor)   
        fold_epochs.append(stopped_epoch)

        # 모델 평가
        avg_loss = evaluate_model(model, val_loader)  
        fold_scores.append(avg_loss)
    
    fold_epoch_str = ', '.join([f'fold {i+1}: {epoch}' for i, epoch in enumerate(fold_epochs)])
    avg_epoch = sum(fold_epochs) / len(fold_epochs)
    print(f'[{fold_epoch_str}, avg_epoch: {avg_epoch:.1f}]')

    trial.set_user_attr('avg_epoch', avg_epoch)
    
    return np.mean(fold_scores)

# Main funciton
sampler = GPSampler(n_startup_trials=10, seed=42)
n_trials = 100
study = optuna.create_study(
    study_name = 'hpOpt',
    direction='minimize',
    sampler=sampler,
    load_if_exists = True)

# optimize the hyperparameters
study.optimize(objective, n_trials=n_trials, n_jobs=5, callbacks=[save_callback])

# return the best hyperparameters
best_params = study.best_trial.params
best_value = study.best_value

# optimize the model with the best hyperparameters
results = {
    'best_params': best_params,
    'best_value': best_value,
    'study': study
}