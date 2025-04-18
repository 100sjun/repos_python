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
import json
import pickle

# Random seed
np.random.seed(42)
torch.manual_seed(42)

# Feedforward network for migration prediction
class MigrationPredictor(nn.Module):
    def __init__(self, hidden_nodes = 64, hidden_layers = 3):
        super().__init__()
        
        n_features = 7
        n_outputs = 4

        # Layer configuration
        layers = []
        # input layer
        layers.append(nn.Linear(n_features, hidden_nodes))
        layers.append(nn.ReLU())

        # hidden layers
        for _ in range(hidden_layers - 1):
            layers.append(nn.Linear(hidden_nodes, hidden_nodes))
            layers.append(nn.ReLU())

        # output layer
        layers.append(nn.Linear(hidden_nodes, n_outputs))

        self.model = nn.Sequential(*layers)

    def forward(self, x):
        return self.model(x)

# Physical Layers for State update
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
    
class BMEDModel(nn.Module):
    def __init__(self, hidden_nodes = 32, hidden_layers = 3, dt = 0.1, scaler = None):
        super().__init__()
        self.migration_predictor = MigrationPredictor(hidden_nodes, hidden_layers)
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
    
def custom_loss(pred_states, measured_indices, measured_states):
    total_loss = 0
    for idx, measured_state in zip(measured_indices, measured_states[0]):
        predicted_state = pred_states[idx]
        total_loss += torch.mean((predicted_state - measured_state) ** 2)

    return total_loss

def train_model(model, train_loader, epochs = 100, learning_rate = 0.001, weight_decay = 1e-5):
    optimizer = torch.optim.Adam(model.parameters(), lr = learning_rate, weight_decay= weight_decay)

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
            optimizer.step()

def r2_calculator(pred_states, measured_indices, measured_states):

    pred_value = pred_states[measured_indices][:,3:].numpy()
    target_value = measured_states[0][:,3:].numpy()
    r2_scores = []
    for col in range(pred_value.shape[1]):
        col_r2 = r2_score(pred_value[:,col], target_value[:,col])
        r2_scores.append(col_r2)

    return np.mean(r2_scores)

def evaluate_model(model, val_loader):
    model.eval()
    total_r2 = 0
    total_samples = 0

    with torch.no_grad():
        for exp in val_loader:
            init_state = exp['init_state']
            measured_states = exp['measured_state']
            times = exp['times']

            pred_states, measured_indices = model(init_state, times)
            r2 = r2_calculator(pred_states, measured_indices, measured_states)

            total_r2 += r2
            total_samples += 1

    avg_r2 = total_r2 / total_samples
    return avg_r2

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

def plot_results(model, dataset, exp_idx=0):
    # 데이터 준비
    exp = dataset.experiments[exp_idx]
    init_state = torch.tensor(exp['init_state'], dtype=torch.float32).unsqueeze(0)
    times = torch.tensor(exp['times'], dtype=torch.float32).unsqueeze(0)
    measured_state = torch.tensor(exp['measured_state'], dtype=torch.float32)
    
    # 예측
    with torch.no_grad():
        pred_states, measured_indices = model(init_state, times)
    
    # 예측을 위한 시간 포인트 생성 (dt = 0.1 간격)
    t_pred = torch.arange(0, times[0][-1].item() + 0.1, 0.1)
    
    # 변수 이름과 단위
    var_names = {
        'T': 'Temperature (°C)',
        'V': 'Voltage (V)',
        'E': 'Electric Field (V/cm)',
        'CF_LA': 'Feed LA Conc. (M)',
        'CA_LA': 'Acid LA Conc. (M)',
        'CF_K': 'Feed K Conc. (M)',
        'CB_K': 'Base K Conc. (M)',
        'VF': 'Feed Volume (L)',
        'VA': 'Acid Volume (L)',
        'VB': 'Base Volume (L)'
    }
    
    # 그래프 그리기
    fig, axes = plt.subplots(4, 3, figsize=(15, 20))
    axes = axes.ravel()
    
    for i, (var, label) in enumerate(var_names.items()):
        ax = axes[i]
        # 실제 측정값 (점으로 표시)
        ax.plot(times[0].numpy(), measured_state[:, i].numpy(), 
                'bo', label='Measured', markersize=6)
        # 예측값 (연속 선으로 표시)
        ax.plot(t_pred.numpy(), pred_states[:, i].numpy(), 
                'r-', label='Predicted', linewidth=2)
        
        ax.set_title(label, fontsize=12, pad=10)
        ax.set_xlabel('Time (hr)', fontsize=10)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()

def save_callback(study, trial):
    with open('hpOpt_study.pkl', 'wb') as f:
        pickle.dump(study, f)

    # 2. best value만 json으로 저장
    best_trial_info = {
        'best_params': study.best_trial.params,
        'best_value': study.best_trial.value,
    }

    with open('hpOpt_checkpoint.json', 'w') as f:
        json.dump(best_trial_info, f, indent=4)

    print(f"[Trial {trial.number}] Checkpoint saved: best value = {study.best_trial.value}")

def objective(trial):
    # load data
    df_path = './data/BMED_data_v8.csv'
    dict_spline, scaler = prepare_data(df_path,exp_idx=[1,2,3,4,5])

    # Hyperparameter 설정
    hidden_nodes = trial.suggest_int('hidden_nodes', 16, 128, step = 16)
    hidden_layers = trial.suggest_int('hidden_layers', 1, 10)
    lr = trial.suggest_float('lr', 1e-5, 1e-2, log = True)
    wd = trial.suggest_float('weight_decay', 1e-6, 1e-3, log = True)
    epochs = trial.suggest_int('epochs', 5, 10)
    dt = 0.1

    # 5-fold 설정
    kfolds = 5
    kf = KFold(n_splits=kfolds, shuffle=True, random_state=42)

    # dataset 생성
    dataset = BMEDDataset(dict_spline)

    # k-fold 교차 검증 스코어 저장
    fold_scores = []

    for fold, (train_idx, val_idx) in enumerate(kf.split(dataset)):
        # kfold로 전체 dataset에 대해서 train과 validation subset을 생성
        train_subset = Subset(dataset, train_idx)
        val_subset = Subset(dataset, val_idx)

        # Dataloader 생성
        train_loader = DataLoader(train_subset, batch_size = 1, shuffle = True)
        val_loader = DataLoader(val_subset, batch_size = 1, shuffle = False)

        # 모델 생성
        model = BMEDModel(hidden_nodes=hidden_nodes, hidden_layers=hidden_layers, dt = dt, scaler=scaler)

        # 모델 학습
        train_model(model = model, train_loader = train_loader, epochs = epochs, learning_rate = lr, weight_decay = wd)   

        # 모델 평가
        r2_score = evaluate_model(model, val_loader)  
        fold_scores.append(r2_score)
    
    return np.mean(fold_scores)

# generate study with gaussian process sampler
sampler = GPSampler(n_startup_trials=23, seed=42)
n_trials = 50
study = optuna.create_study(
    study_name = 'hpOpt',
    direction='maximize',
    sampler=sampler,
    load_if_exists = True)

# optimize the hyperparameters
study.optimize(objective, n_trials=n_trials, n_jobs=23, callbacks=[save_callback])

# return the best hyperparameters
best_params = study.best_trial.params
best_value = study.best_value

# optimize the model with the best hyperparameters
results = {
    'best_params': best_params,
    'best_value': best_value,
    'study': study
}

save_data = {
            'best_params': results['best_params'],
            'best_r2_score': float(results['best_value']),
        }

with open('hpOpt.json', 'w') as f:
    json.dump(save_data, f, indent=4)

with open('hpOpt.pkl', 'wb') as f:
    pickle.dump(results, f)

