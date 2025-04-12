import torch
import torch.nn as nn
import numpy as np
from sklearn.metrics import r2_score

'''
멤브레인 시스템 모델 훈련
'''

class MembraneSystemTrainer:
    def __init__(self, model, device = 'cuda' if torch.cuda.is_available() else 'cpu', epochs=100, lr = 0.001, rstop = 0.1, weight_decay = 1e-2):
        self.model = model.to(device)
        self.device = device
        self.epochs = epochs
        self.lr = lr
        self.rstop = rstop
        self.weight_decay = weight_decay
        self.optimizer = torch.optim.Adam(self.model.parameters(), lr=self.lr, weight_decay=self.weight_decay)
        self.criterion = self.custom_max_mse_loss

    def train(self, train_loader, val_loader=None):
        # model 훈련
        self.model.train()
        train_losses = []
        val_losses = []

        # Ealry stopping 변수 초기화
        best_val_loss = float('inf')
        best_model_state = None
        counter = 0
        patience = self.epochs * self.rstop

        # R2
        r2_final = None
        r2_min = None

        for epoch in range(self.epochs):
            epoch_loss = 0.0

            # training part
            for features, migrations, states in train_loader:
                # move data to device
                features = features.to(self.device)
                states = states.to(self.device)

                # initialize gradient
                self.optimizer.zero_grad()

                # forward pass
                migration, state = self.model(features)

                # calculate loss
                loss = self.criterion(state, states)

                # backpropagation
                loss.backward()

                # optimization
                self.optimizer.step()

                # accumulate loss
                epoch_loss += loss.item()

            # epoch average loss
            avg_loss = epoch_loss / len(train_loader)
            train_losses.append(avg_loss)

            # validation
            if val_loader is not None:
                self.model.eval()
                val_loss = 0.0
                val_preds = []
                val_actuals = []

                with torch.no_grad():
                    for features, migrations, states in val_loader:
                        features = features.to(self.device)
                        states = states.to(self.device)

                        # forward pass
                        migration, state = self.model(features)

                        # calculate loss
                        loss = self.criterion(state, states)
                        
                        # accumulate loss
                        val_loss += loss.item()

                        # save actual and predicted values
                        val_preds.append(state.cpu().numpy())
                        val_actuals.append(states.cpu().numpy())

                # validation average loss
                avg_val_loss = val_loss / len(val_loader)
                val_losses.append(avg_val_loss)

                val_preds = np.concatenate(val_preds, axis=0)
                val_actuals = np.concatenate(val_actuals, axis=0)

                # calculate R2 at this epoch
                r2_values = {}
                var_names = ['T', 'V', 'E', 'CF_LA', 'CF_K', 'CA_LA', 'CB_K', 'VF', 'VA', 'VB']

                for i, var in enumerate(var_names):
                    r2 = r2_score(val_actuals[:, i], val_preds[:, i])
                    r2_values[var] = r2
                
                min_r2 = min(r2_values.values())
                worst_var = min(r2_values.keys(), key= lambda k: r2_values[k])

                # Early stopping 체크
                if avg_val_loss < best_val_loss:
                    best_val_loss = avg_val_loss
                    best_model_state = self.model.state_dict().copy()
                    # 최고 성능 모델의 R2 값 저장
                    r2_final = r2_values.copy()
                    r2_min= min_r2
                    counter = 0
                else:
                    counter += 1
                    print(f'', end='')
                    if counter >= patience:
                        # restore best model parameters
                        self.model.load_state_dict(best_model_state)
                        return train_losses, val_losses, r2_final, r2_min, worst_var
                # Print epoch, train loss, val loss
                # print(f'\rEpoch {epoch+1}/{self.epochs}, Train Loss: {avg_loss:.6f}, Val Loss: {avg_val_loss:.6f}, Best Val Loss: {best_val_loss:.6f}, Min R2: {r2_min:.6f}, Worst Var: {worst_var}, EarlyStopping counter: {counter} out of {patience}                   ', end='')

                self.model.train() # set the train mode after validation
            # else:
            #     print(f'\rEpoch {epoch+1}/{self.epochs}, Train Loss: {avg_loss:.6f}', end='')

        # All epochs completion
        if best_model_state is not None:
            self.model.load_state_dict(best_model_state)

        return train_losses, val_losses, r2_final, r2_min, worst_var
    
    def custom_max_mse_loss(self,predictions, targets):

        squared_diff  = (predictions - targets) ** 2

        if squared_diff.dim() > 1:  # 배치 차원이 있는 경우
            # 각 변수(특성)별 MSE 계산
            mse_per_var = torch.mean(squared_diff, dim=0)  # 배치 차원을 따라 평균
        else:  # 배치 차원이 없는 경우 (이미 평균화된 경우)
            mse_per_var = squared_diff
        
        # 변수 중 가장 큰 MSE 반환
        return torch.max(mse_per_var)
