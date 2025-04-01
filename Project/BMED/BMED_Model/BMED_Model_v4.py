import pandas as pd
import torch
import torch.nn as nn
import numpy as np
import matplotlib.pyplot as plt
import os
import optuna
import json
from sklearn.model_selection import KFold
from torch.utils.data import DataLoader, TensorDataset
import torch.optim.lr_scheduler as lr_scheduler
from pathlib import Path
import seaborn as sns

# Fix OpenMP runtime duplicate initialization error
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

# Load data
df = pd.read_csv('BMED_data_v5.csv')

# Prepare data for all experiments
X_batches = []
S_batches = []
Y_batches = []
times_list = []

max_points = max(len(df[df['exp'] == exp_id]) for exp_id in df['exp'].unique())

for exp_id in df['exp'].unique():
    df_exp = df[df['exp'] == exp_id]
    
    X = df_exp[['T','V','E','Ci','Ki']].iloc[0]
    S = df_exp[['NF_LA','NA_LA','NF_K','NB_K','VF','VA','VB']].iloc[0]
    Y = df_exp[['NF_LA','NA_LA','NF_K','NB_K','VF','VA','VB']].values
    times = df_exp['t'].values
    
    # Pad Y with the last values to match max_points
    if len(Y) < max_points:
        pad_length = max_points - len(Y)
        Y = np.pad(Y, ((0, pad_length), (0, 0)), mode='edge')
        times = np.pad(times, (0, pad_length), mode='edge')
    
    X_batches.append(torch.FloatTensor(X.values))
    S_batches.append(torch.FloatTensor(S.values))
    Y_batches.append(torch.FloatTensor(Y))
    times_list.append(times)

# Convert to batches
X_batch = torch.stack(X_batches)
S_batch = torch.stack(S_batches)
Y_batch = torch.stack(Y_batches)

# Model parameter settings
x_size = X_batch[0].shape[0]  # Update x_size based on first experiment
s_size = S_batch[0].shape[0]  # Update s_size based on first experiment
dt = 0.1
max_time = times_list[0][-1]  # First experiment's end time

# CUDA 사용 가능 여부 확인 및 device 설정
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"Using device: {device}")

# 데이터를 GPU로 이동
X_batch = X_batch.to(device)
S_batch = S_batch.to(device)
Y_batch = Y_batch.to(device)

class BMEDLSTM(nn.Module):
    def __init__(self, x_size, s_size, hidden_size, num_layers=2, dropout_rate=0.2):
        super(BMEDLSTM, self).__init__()
        self.hidden_size = hidden_size
        self.num_layers = num_layers
        self.s_size = s_size
        
        self.lstm = nn.LSTM(x_size + s_size, hidden_size, num_layers, 
                           batch_first=True, dropout=dropout_rate)
        self.dropout = nn.Dropout(dropout_rate)
        self.fc = nn.Linear(hidden_size, s_size)
        
    def forward(self, x, s0, dt, max_time):
        batch_size = x.size(0)
        h0 = torch.zeros(self.num_layers, batch_size, self.hidden_size).to(x.device)
        c0 = torch.zeros(self.num_layers, batch_size, self.hidden_size).to(x.device)
        
        # 초기 상태 설정
        current_s = s0.unsqueeze(1)  # [batch_size, 1, s_size]
        outputs = []
        res_step = []
        cal_times = []
        cal_t = 0.0
        
        # dt 간격으로 계산 진행
        while cal_t <= max_time:
            # 현재 입력 x와 상태 s 결합
            x_expanded = x.unsqueeze(1)  # [batch_size, 1, x_size]
            combined_input = torch.cat([x_expanded, current_s], dim=2)  # [batch_size, 1, x_size + s_size]
            
            # LSTM 순전파
            lstm_out, (h0, c0) = self.lstm(combined_input, (h0, c0))
            lstm_out = self.dropout(lstm_out)  # Dropout 적용
            step_output = self.fc(lstm_out)
            
            # 상태 업데이트 (S(n+1) = S(n) + output*dt)
            current_s = current_s + step_output * dt
            
            res_step.append(step_output)
            outputs.append(current_s)
            cal_times.append(cal_t)
            cal_t += dt
        
        # 모든 출력을 시퀀스로 결합
        outputs = torch.cat(outputs, dim=1)  # [batch_size, num_steps, s_size]
        cal_times = torch.tensor(cal_times).to(x.device)
        
        return outputs, cal_times, res_step

def train_model(model, X_train, S_train, Y_train, X_val, S_val, Y_val, val_idx, params, dt, criterion, optimizer, run_idx, total_runs, fold_idx):
    model.train()
    patience = int(params['epochs'] * 0.1)
    best_val_loss = float('inf')
    patience_counter = 0
    
    for epoch in range(params['epochs']):
        # Training
        model.train()
        train_loss = 0
        for i in range(0, len(X_train), params['batch_size']):
            batch_X = X_train[i:i+params['batch_size']]
            batch_S = S_train[i:i+params['batch_size']]
            batch_Y = Y_train[i:i+params['batch_size']]
            batch_times = times_list[i:i+params['batch_size']]
            
            optimizer.zero_grad()
            # 각 실험의 실제 종결 시간 사용
            batch_max_time = max(times[-1] for times in batch_times)
            outputs, cal_times, _ = model(batch_X, batch_S, dt, batch_max_time)
            
            loss = 0
            for j in range(len(batch_X)):
                for t_val in batch_times[j]:
                    idx = torch.abs(cal_times - t_val).argmin()
                    t_idx = np.where(batch_times[j] == t_val)[0][0]
                    loss += criterion(outputs[j:j+1, idx], batch_Y[j:j+1, t_idx])
            
            loss.backward()
            optimizer.step()
            train_loss += loss.item()
        
        
        # Validation
        model.eval()
        with torch.no_grad():
            # validation set의 실제 종결 시간 사용
            val_max_time = max(times[-1] for times in [times_list[i] for i in val_idx])
            outputs, cal_times, _ = model(X_val, S_val, dt, val_max_time)
            val_loss = 0
            for j in range(len(X_val)):
                for t_idx, t_val in enumerate(times_list[val_idx[j]]):
                    idx = torch.abs(cal_times - t_val).argmin()
                    val_loss += criterion(outputs[j:j+1, idx], Y_val[j, t_idx].unsqueeze(0))
        
        # Early stopping with validation loss
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            patience_counter = 0
            best_model_state = model.state_dict()
        else:
            patience_counter += 1
            if patience_counter >= patience:
                model.load_state_dict(best_model_state)
                break
        
        # 10 에포크마다 진행 상황 업데이트
        if (epoch + 1) % 1 == 0:
            print(f"\rRun [{run_idx}/{total_runs}] Fold [{fold_idx}/5] Epoch: {epoch + 1}/{params['epochs']}, "
                  f"Train Loss: {train_loss:.4f}, Val Loss: {val_loss:.4f}", end="", flush=True)
    
    print()  # 줄바꿈
    return model, best_val_loss

def objective(trial):
    # 하이퍼파라미터 정의 - 더 작은 범위로 조정
    params = {
        'hidden_size': trial.suggest_int('hidden_size', 8, 16),  # 범위 축소
        'num_layers': trial.suggest_int('num_layers', 2, 4),      # 범위 축소
        'dropout_rate': trial.suggest_float('dropout_rate', 0.1, 0.3),
        'learning_rate': trial.suggest_float('learning_rate', 1e-4, 1e-3, log=True),
        'batch_size': trial.suggest_int('batch_size', 4, 8),      # 범위 축소
        'epochs': trial.suggest_int('epochs', 10, 20)            # 범위 축소
    }
    
    # 5-fold cross validation
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    cv_scores = []
    
    for fold, (train_idx, val_idx) in enumerate(kf.split(X_batch)):
        # 데이터 분할
        X_train = X_batch[train_idx]
        S_train = S_batch[train_idx]
        Y_train = Y_batch[train_idx]
        X_val = X_batch[val_idx]
        S_val = S_batch[val_idx]
        Y_val = Y_batch[val_idx]
        
        # 모델 초기화 및 GPU로 이동
        model = BMEDLSTM(x_size, s_size, params['hidden_size'], 
                        params['num_layers'], params['dropout_rate']).to(device)
        
        # 학습
        criterion = nn.MSELoss()
        optimizer = torch.optim.Adam(model.parameters(), lr=params['learning_rate'])
        
        model, val_loss = train_model(model, X_train, S_train, Y_train, X_val, S_val, Y_val, val_idx,
                                    params, dt, criterion, optimizer, trial.number + 1, 100, fold + 1)
        cv_scores.append(val_loss.cpu())
    
    # numpy로 변환하여 평균 계산
    return np.mean([score.numpy() for score in cv_scores])

# Define y_names for plotting
y_names = ['NF_LA', 'NA_LA', 'NF_K', 'NB_K', 'VF', 'VA', 'VB']

def plot_results(model, X, S, Y, times, exp_id, save_path):
    model.eval()
    with torch.no_grad():
        # GPU에서 CPU로 데이터 이동
        X = X.to(device)
        S = S.to(device)
        exp_max_time = times[-1]
        predictions, cal_times, fluxes = model(X, S, dt, exp_max_time)
        
        # Convert to numpy (GPU -> CPU -> numpy)
        cal_times_np = cal_times.cpu().numpy()
        predictions_np = predictions[0].cpu().numpy()
        fluxes_np = torch.cat(fluxes, dim=1)[0].cpu().numpy()
        
        # Plot
        plt.figure(figsize=(15, 10))
        for i, name in enumerate(y_names):
            plt.subplot(3, 3, i+1)
            plt.scatter(times, Y[0, :, i].numpy(), c='b', label='Experimental', alpha=0.6)
            plt.plot(cal_times_np, predictions_np[:, i], 'r-', label='Model Prediction', linewidth=1.5)
            plt.title(f'{name}')
            plt.xlabel('Time (t)')
            plt.ylabel('Value')
            plt.legend()
            plt.grid(True)
        
        plt.tight_layout()
        plt.savefig(f'{save_path}/exp_{exp_id}.png')
        plt.close()
        
        # Save results to CSV
        results_df = pd.DataFrame()
        results_df['time'] = cal_times_np
        for i, name in enumerate(y_names):
            results_df[name] = predictions_np[:, i]
            results_df[f'd{name}/dt'] = fluxes_np[:, i]
        results_df.to_csv(f'{save_path}/exp_{exp_id}.csv', index=False)

# Create directories for results
results_dir = Path('results')
results_dir.mkdir(exist_ok=True)

# Hyperparameter optimization
study = optuna.create_study(direction='minimize')
study.optimize(objective, n_trials=10)

# Save optimization results
with open('LSTM_opt.json', 'w') as f:
    json.dump(study.best_params, f, indent=4)

# Print optimization results
print("\nBest hyperparameters:")
for key, value in study.best_params.items():
    print(f"{key}: {value}")

# Final training with best hyperparameters
best_params = study.best_params
best_params['epochs'] *= 10  # Increase epochs by 10

# Split data for final training (8:2)
train_size = int(0.8 * len(X_batch))
X_train = X_batch[:train_size]
S_train = S_batch[:train_size]
Y_train = Y_batch[:train_size]
X_test = X_batch[train_size:]
S_test = S_batch[train_size:]
Y_test = Y_batch[train_size:]

# 최종 모델 초기화 및 GPU로 이동
final_model = BMEDLSTM(x_size, s_size, best_params['hidden_size'], 
                      best_params['num_layers'], best_params['dropout_rate']).to(device)

# Train final model
criterion = nn.MSELoss()
optimizer = torch.optim.Adam(final_model.parameters(), lr=best_params['learning_rate'])

# Create test indices for validation
test_indices = list(range(train_size, len(X_batch)))

final_model, _ = train_model(final_model, X_train, S_train, Y_train, X_test, S_test, Y_test, 
                           test_indices, best_params, dt, criterion, optimizer, 1, 1, 0)

# Save final model
torch.save({
    'model_state_dict': final_model.state_dict(),
    'optimizer_state_dict': optimizer.state_dict(),
    'hyperparameters': best_params,
}, 'final_model.pth')

# Plot and save results for each experiment
for i in range(len(X_batch)):
    plot_results(final_model, X_batch[i:i+1], S_batch[i:i+1], 
                Y_batch[i:i+1], times_list[i], i, results_dir)

print("\nTraining completed. Results saved in 'results' directory.")






