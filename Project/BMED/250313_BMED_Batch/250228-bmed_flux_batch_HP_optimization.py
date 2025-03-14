# import modules
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import optuna
from optuna.samplers import GPSampler
import warnings
import matplotlib.pyplot as plt
import numpy as np
import json


# 경고 무시
warnings.filterwarnings("ignore", category=UserWarning)

# Hyperparameter Trial
t_trial = 100

# load data
class RawDataLoader():
    def __init__(self, path='250314_BMED_train data_v4.xlsx'):
        self.path = path
        self.X_data, self.Y_data = self.RawData()

    def RawData(self):
        df = pd.read_excel(self.path, sheet_name='Sheet2')
        X_data = df[['T','V','E','CF_LA','CF_K','CA_LA','CB_K']].values
        Y_data = df[['dNLA','dNK','dVF','dVA','dVB']].values
        return X_data, Y_data
    
    def FoldData(self):
        kf = KFold(n_splits=5, shuffle=True, random_state=42)
        folds = []
        for train_index, test_index in kf.split(self.X_data):
            # Split the data into training and test sets
            X_train, X_test = self.X_data[train_index], self.X_data[test_index]
            Y_train, Y_test = self.Y_data[train_index], self.Y_data[test_index]
            
            # Normalize the data
            scaler_X = StandardScaler()
            scaler_Y = StandardScaler()
            
            X_train_scaled = scaler_X.fit_transform(X_train)
            X_test_scaled = scaler_X.transform(X_test)
            
            Y_train_scaled = scaler_Y.fit_transform(Y_train)
            Y_test_scaled = scaler_Y.transform(Y_test)

            X_train_tensor = torch.FloatTensor(X_train_scaled).to('cuda')
            X_test_tensor = torch.FloatTensor(X_test_scaled).to('cuda')
            Y_train_tensor = torch.FloatTensor(Y_train_scaled).to('cuda')
            Y_test_tensor = torch.FloatTensor(Y_test_scaled).to('cuda')

            folds.append((X_train_tensor, X_test_tensor, Y_train_tensor, Y_test_tensor, scaler_X, scaler_Y))
        return folds
    
# Customize the NN architecture
class CustomModel(nn.Module):
    def __init__(self, hidden_layers=2, hidden_nodes = 8):
        super().__init__()
        layers = []
        nodes = 7
        for _ in range(hidden_layers):
            layers.append(nn.Linear(nodes, hidden_nodes))
            layers.append(nn.ReLU())
            nodes = hidden_nodes
        layers.append(nn.Linear(hidden_nodes, 5))
        self.hidden = nn.Sequential(*layers)

    def forward(self, x):
        return self.hidden(x)

# Hyperparameter optimization
class NNOpt():
    def __init__(self, hidden_layers=2, hidden_nodes = 8, learning_rate=0.001, num_epochs=1000, batch_size=256, weight_decay=1e-5):
        self.hidden_layers= hidden_layers
        self.hidden_nodes = hidden_nodes
        self.learning_rate = learning_rate
        self.num_epochs = num_epochs
        self.batch_size = batch_size
        self.weight_decay = weight_decay
        self.model = CustomModel(hidden_layers=self.hidden_layers, hidden_nodes=self.hidden_nodes).to('cuda')
        self.criterion = nn.MSELoss()
        self.optimizer = optim.Adam(self.model.parameters(), lr=self.learning_rate, weight_decay=self.weight_decay)
        self.train_losses = []
        self.test_losses = []

    def train(self, X_train, Y_train, X_test, Y_test, fold, trial_num, scaler_Y=None):
        X_train_gpu = X_train
        Y_train_gpu = Y_train
        X_test_gpu = X_test
        Y_test_gpu = Y_test

        dataset = TensorDataset(X_train_gpu, Y_train_gpu)
        dataloader = DataLoader(dataset, batch_size=self.batch_size, shuffle=True)

        for epoch in range(self.num_epochs):
            self.model.train()
            epoch_loss = 0
            for X_batch, Y_batch in dataloader:
                X_batch = X_batch.to('cuda')
                Y_batch = Y_batch.to('cuda')
                self.optimizer.zero_grad()
                train_outputs = self.model(X_batch)
                train_loss = self.criterion(train_outputs, Y_batch)
                train_loss.backward()
                self.optimizer.step()
                epoch_loss += train_loss.item()
            
            # 학습 곡선을 위한 loss 기록
            with torch.no_grad():
                train_pred = self.model(X_train_gpu)
                test_pred = self.model(X_test_gpu)
                self.train_losses.append(self.criterion(train_pred, Y_train_gpu).item())
                self.test_losses.append(self.criterion(test_pred, Y_test_gpu).item())
                
            if (epoch + 1) % 10 == 0:
                print(f'\rTrial [{trial_num}/{t_trial-1}], fold: {fold+1}/5, Epoch [{epoch + 1}/{self.num_epochs}]', end='',flush=True)
        
        with torch.no_grad():
            test_outputs = self.model(X_test_gpu)
            
            # R2 점수 계산
            y_mean = torch.mean(Y_test_gpu, dim=0)
            ss_tot = torch.sum((Y_test_gpu - y_mean) ** 2, dim=0)
            ss_res = torch.sum((Y_test_gpu - test_outputs) ** 2, dim=0)
            r2 = torch.mean(1 - ss_res / ss_tot)  # 모든 출력의 평균 R2 점수
            
        return r2.item()  # R2 점수 반환

# Objective function
def objective(trial):
    hidden_layers = trial.suggest_int('hidden_layers', 1, 20)
    hidden_nodes = trial.suggest_int('hidden_nodes', 8, 128)
    learning_rate = trial.suggest_float('learning_rate', 0.00001, 0.1, log=True)
    num_epochs = trial.suggest_int('num_epochs', 100, 1000)
    batch_size = trial.suggest_int('batch_size', 16, 256)
    weight_decay = trial.suggest_float('weight_decay', 1e-7, 1e-2, log=True)

    data = RawDataLoader()
    folds = data.FoldData()

    r2_scores = []

    for fold, (X_train, X_test, Y_train, Y_test, _, _) in enumerate(folds):
        model = NNOpt(
            hidden_layers=hidden_layers, 
            hidden_nodes=hidden_nodes, 
            learning_rate=learning_rate, 
            num_epochs=num_epochs, 
            batch_size=batch_size, 
            weight_decay=weight_decay)
        r2_score = model.train(X_train, Y_train, X_test, Y_test, fold, trial.number)
        r2_scores.append(r2_score)
    print('\n')
        
    return sum(r2_scores) / len(r2_scores)

# Optuna Study 설정 부분 수정
def create_study():
    sampler = GPSampler()
    study_name = "optimization_study"
    storage_name = "sqlite:///study.db"  # SQLite DB에 결과 저장
    return optuna.create_study(
        study_name=study_name,
        storage=storage_name,
        direction="maximize",
        sampler=sampler,
        load_if_exists=True  # 기존 study가 있으면 불러오기
    )

# 최적화 진행 상황을 저장하는 콜백 함수 수정
def save_study_callback(study, trial):
    # 현재 trial이 지금까지의 최고 성능을 달성한 경우에만 저장
    if study.best_trial.number == trial.number:
        print(f"\n새로운 최적 성능 발견: {study.best_value:.4f}")
        print(f"최적 파라미터: {study.best_params}")
        
        # 최적 결과를 JSON 파일로 저장
        best_result = {
            "best_value": study.best_value,
            "best_params": study.best_params,
            "n_trials": len(study.trials),
            "trial_number": trial.number
        }
        with open('optimization_checkpoint.json', 'w') as f:
            json.dump(best_result, f, indent=4)

# 최적화 실행 부분 수정
def run_optimization(study):
    study.optimize(objective, n_trials=t_trial, callbacks=[save_study_callback])

def plot_optimization_history(study):
    plt.figure(figsize=(10, 6))
    plt.plot(study.trials_dataframe()['number'], study.trials_dataframe()['value'])
    plt.xlabel('Trial number')
    plt.ylabel('Objective value (R2 score)')
    plt.title('Optimization History')
    plt.savefig('optimization_history.png')
    plt.close()

def train_final_model(best_params, X_data, Y_data, scaler_X, scaler_Y):
    # 데이터 분할 (80:20)
    train_size = int(0.8 * len(X_data))
    indices = np.random.permutation(len(X_data))
    train_idx, test_idx = indices[:train_size], indices[train_size:]
    
    X_train, X_test = X_data[train_idx], X_data[test_idx]
    Y_train, Y_test = Y_data[train_idx], Y_data[test_idx]
    
    # 데이터 정규화
    X_train_scaled = scaler_X.fit_transform(X_train)
    X_test_scaled = scaler_X.transform(X_test)
    Y_train_scaled = scaler_Y.fit_transform(Y_train)
    Y_test_scaled = scaler_Y.transform(Y_test)
    
    # 텐서 변환
    X_train_tensor = torch.FloatTensor(X_train_scaled).to('cuda')
    X_test_tensor = torch.FloatTensor(X_test_scaled).to('cuda')
    Y_train_tensor = torch.FloatTensor(Y_train_scaled).to('cuda')
    Y_test_tensor = torch.FloatTensor(Y_test_scaled).to('cuda')
    
    # 모델 학습
    model = NNOpt(**best_params)
    model.train(X_train_tensor, Y_train_tensor, X_test_tensor, Y_test_tensor, 0, 0, scaler_Y)
    
    return model, (X_train_tensor, Y_train_tensor, X_test_tensor, Y_test_tensor)

def plot_learning_curves(model):
    plt.figure(figsize=(10, 6))
    plt.plot(model.train_losses, label='Training Loss')
    plt.plot(model.test_losses, label='Validation Loss')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.title('Learning Curves')
    plt.legend()
    plt.savefig('learning_curves.png')
    plt.close()

def plot_predictions(model, data_tensors, scaler_Y, title_suffix=''):
    X_train_tensor, Y_train_tensor, X_test_tensor, Y_test_tensor = data_tensors
    
    with torch.no_grad():
        train_pred = model.model(X_train_tensor).cpu().numpy()
        test_pred = model.model(X_test_tensor).cpu().numpy()
        
    # 스케일 복원
    Y_train = scaler_Y.inverse_transform(Y_train_tensor.cpu().numpy())
    Y_test = scaler_Y.inverse_transform(Y_test_tensor.cpu().numpy())
    train_pred = scaler_Y.inverse_transform(train_pred)
    test_pred = scaler_Y.inverse_transform(test_pred)
    
    output_names = ['dNLA', 'dNK', 'dVF', 'dVA', 'dVB']
    
    for i, name in enumerate(output_names):
        plt.figure(figsize=(10, 6))
        
        # R2 점수 계산
        def calculate_r2(y_true, y_pred):
            ss_res = np.sum((y_true - y_pred) ** 2)
            ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
            r2 = 1 - (ss_res / ss_tot)
            return r2
        
        r2_train = calculate_r2(Y_train[:, i], train_pred[:, i])
        r2_test = calculate_r2(Y_test[:, i], test_pred[:, i])
        
        # Training data
        plt.scatter(Y_train[:, i], train_pred[:, i], alpha=0.5, label=f'Training (R² = {r2_train:.3f})')
        # Test data
        plt.scatter(Y_test[:, i], test_pred[:, i], alpha=0.5, label=f'Test (R² = {r2_test:.3f})')
        
        # Add diagonal line
        min_val = min(min(Y_train[:, i]), min(train_pred[:, i]))
        max_val = max(max(Y_train[:, i]), max(train_pred[:, i]))
        plt.plot([min_val, max_val], [min_val, max_val], 'r--')
        
        plt.xlabel(f'Actual {name}')
        plt.ylabel(f'Predicted {name}')
        plt.title(f'Prediction vs Actual for {name} {title_suffix}')
        plt.legend()
        plt.savefig(f'prediction_plot_{name}{title_suffix}.png')
        plt.close()



# 모델 저장
def save_model(model, scaler_X, scaler_Y, best_params, filepath='final_model.pth'):
    model_state = {
        'model_state_dict': model.model.state_dict(),
        'optimizer_state_dict': model.optimizer.state_dict(),
        'hyperparameters': {
            'hidden_layers': best_params['hidden_layers'],
            'hidden_nodes': best_params['hidden_nodes'],
            'learning_rate': best_params['learning_rate'],
            'num_epochs': best_params['num_epochs'],
            'batch_size': best_params['batch_size'],
            'weight_decay': best_params['weight_decay']
        },
        'train_losses': model.train_losses,
        'test_losses': model.test_losses,
        'scalers': {
            'scaler_X': scaler_X,
            'scaler_Y': scaler_Y
        }
    }
    torch.save(model_state, filepath)
    print(f"Model saved to {filepath}")

def load_model(filepath='final_model.pth'):
    model_state = torch.load(filepath)
    
    # 하이퍼파라미터로 모델 초기화
    model = NNOpt(**model_state['hyperparameters'])
    
    # 모델 가중치 로드
    model.model.load_state_dict(model_state['model_state_dict'])
    model.optimizer.load_state_dict(model_state['optimizer_state_dict'])
    
    # 학습 이력 로드
    model.train_losses = model_state['train_losses']
    model.test_losses = model_state['test_losses']
    
    # 스케일러 로드
    scaler_X = model_state['scalers']['scaler_X']
    scaler_Y = model_state['scalers']['scaler_Y']
    
    print(f"Model loaded from {filepath}")
    return model, scaler_X, scaler_Y

if __name__ == "__main__":
    # study 객체 생성
    study = create_study()
    
    # 최적화 실행
    run_optimization(study)
    
    # 최적화 히스토리 플롯
    plot_optimization_history(study)
    
    # 최적 파라미터로 최종 모델 학습
    best_params = study.best_params
    best_params['num_epochs'] = best_params['num_epochs'] * 10  # 최적화된 epoch의 10배로 설정
    
    data_loader = RawDataLoader()
    X_data, Y_data = data_loader.X_data, data_loader.Y_data
    scaler_X, scaler_Y = StandardScaler(), StandardScaler()
    final_model, data_tensors = train_final_model(best_params, X_data, Y_data, scaler_X, scaler_Y)
    
    # 학습 곡선 플롯
    plot_learning_curves(final_model)
    
    # 예측 결과 플롯
    plot_predictions(final_model, data_tensors, scaler_Y)
    
    # 모델 저장
    save_model(final_model, scaler_X, scaler_Y, best_params)
