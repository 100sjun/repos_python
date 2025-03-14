import torch
import torch.nn as nn
import numpy as np

class CustomModel(nn.Module):
    def __init__(self, hidden_layers=4, hidden_nodes=128):
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

class PredictionModel:
    def __init__(self, model_path='bmed_flux_batch_NN_model_v0.pth'):
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.model_state = torch.load(model_path,weights_only=False)
        
        # 하이퍼파라미터 로드
        hyperparameters = self.model_state['hyperparameters']
        self.model = CustomModel(
            hidden_layers=hyperparameters['hidden_layers'],
            hidden_nodes=hyperparameters['hidden_nodes']
        ).to(self.device)
        
        # 모델 가중치 로드
        self.model.load_state_dict(self.model_state['model_state_dict'])
        self.model.eval()  # 예측 모드로 설정
        
        # 스케일러 로드
        self.scaler_X = self.model_state['scalers']['scaler_X']
        self.scaler_Y = self.model_state['scalers']['scaler_Y']

    def predict(self, X):
        """
        입력:
        X: shape (n_samples, 7) - [T, V, E, CF_LA, CF_K, CA_LA, CB_K]
        
        출력:
        예측값: shape (n_samples, 5) - [dNLA, dNK, dVF, dVA, dVB]
        """
        # 입력 데이터 스케일링
        X_scaled = self.scaler_X.transform(X)
        
        # 텐서 변환 및 예측
        with torch.no_grad():
            X_tensor = torch.FloatTensor(X_scaled).to(self.device)
            y_pred_scaled = self.model(X_tensor)
            
        # 예측값 역스케일링
        y_pred = self.scaler_Y.inverse_transform(y_pred_scaled.cpu().numpy())
        
        return y_pred