# import modules
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt

# load data
class RawDataLoader():
    def __init__(self, path='BMED_train_data_v2.xlsx'):
        self.path = path
        self.X_data, self.Y_data = self.RawData()
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    def RawData(self):
        df = pd.read_excel(self.path, sheet_name='Sheet2')
        X_data = df[['T','V','E','CF_LA','CF_K','CA_LA','CB_K']].values
        Y_data = df[['dNLA','dNK','dVF','dVA','dVB']].values
        return X_data, Y_data
    
    def PrepareData(self, test_size=0.2, random_state=42):
        # Split the data into training and test sets
        X_train, X_test, Y_train, Y_test = train_test_split(
            self.X_data, self.Y_data, 
            test_size=test_size, 
            random_state=random_state
            )
        
        # Normalize the data
        scaler_X = StandardScaler()
        scaler_Y = StandardScaler()

        X_train_scaled = scaler_X.fit_transform(X_train)
        X_test_scaled = scaler_X.transform(X_test)

        Y_train_scaled = scaler_Y.fit_transform(Y_train)
        Y_test_scaled = scaler_Y.transform(Y_test)

        # Convert to PyTorch tensors and move to appropriate device
        X_train_tensor = torch.FloatTensor(X_train_scaled).to(self.device)
        X_test_tensor = torch.FloatTensor(X_test_scaled).to(self.device)
        Y_train_tensor = torch.FloatTensor(Y_train_scaled).to(self.device)
        Y_test_tensor = torch.FloatTensor(Y_test_scaled).to(self.device)

        return X_train_tensor, X_test_tensor, Y_train_tensor, Y_test_tensor, scaler_X, scaler_Y
    
# Customize the NN architecture
class CustomModel(nn.Module):
    def __init__(self, hidden_layers=10, hidden_nodes = 85):
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

# Set the hyperparameters
class NNmodel():
    def __init__(self, hidden_layers=10, hidden_nodes = 85, learning_rate=0.000266709606875957, num_epochs=9271, batch_size=16, weight_decay=5.2171808135365869e-05, title='bmed_flux_batch_NN_model.pth'):
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.hidden_layers = hidden_layers
        self.hidden_nodes = hidden_nodes
        self.learning_rate = learning_rate
        self.num_epochs = num_epochs
        self.batch_size = batch_size
        self.weight_decay = weight_decay
        self.model = CustomModel(hidden_layers=self.hidden_layers, hidden_nodes=self.hidden_nodes).to(self.device)
        self.criterion = nn.MSELoss()
        self.optimizer = optim.Adam(self.model.parameters(), lr=self.learning_rate, weight_decay=self.weight_decay)
        self.train_losses = []
        self.title = title
        self.scaler_X = None
        self.scaler_Y = None

        if torch.cuda.is_available():
            print(f'Using GPU: {torch.cuda.get_device_name(0)}')
            print(f'Memory Usage: {torch.cuda.memory_allocated(0)/1024**2:.2f}MB')
        

    def train(self, X_train, Y_train, X_test, Y_test):
        try:
            X_train_gpu = X_train
            Y_train_gpu = Y_train
            X_test_gpu = X_test
            Y_test_gpu = Y_test

            dataset = TensorDataset(X_train_gpu, Y_train_gpu)
            dataloader = DataLoader(dataset, batch_size=self.batch_size, shuffle=True)
            
            for epoch in range(self.num_epochs):
                self.model.train()
                epoch_loss = 0
                batch_count = 0

                for X_batch, Y_batch in dataloader:
                    X_batch = X_batch.to(self.device)
                    Y_batch = Y_batch.to(self.device)
                    self.optimizer.zero_grad()
                    train_outputs = self.model(X_batch)
                    train_loss = self.criterion(train_outputs, Y_batch)
                    train_loss.backward()
                    self.optimizer.step()

                    epoch_loss += train_loss.item()
                    batch_count += 1
                # Calculate the average loss for the epoch
                avg_loss = epoch_loss / batch_count
                self.train_losses.append(avg_loss)

                if (epoch + 1) % 10 == 0:
                    print(f'\rEpoch [{epoch}/{self.num_epochs}] Train Loss: {avg_loss:.4f}', end='',flush=True)
                
                # 주기적으로 메모리 정리
                if torch.cuda.is_available():
                    torch.cuda.empty_cache()
            
            # Model evaluation with test set
            self.model.eval()
            with torch.no_grad():
                test_outputs = self.model(X_test_gpu)
                test_loss = self.criterion(test_outputs, Y_test_gpu)
            
            # Visualize the learning curves
            self.plot_learning_curves()

            # Save the model
            self.save_model(self.title)

            return test_loss.item()
        
        except RuntimeError as e:
            print(f"GPU 메모리 에러: {e}")
            # 메모리 해제
            del X_train_gpu, Y_train_gpu
            torch.cuda.empty_cache()
            raise
            
        finally:
            # 학습 완료 후 메모리 정리
            if torch.cuda.is_available():
                torch.cuda.empty_cache()

    def save_model(self, title):
        model_state = {
            'model_state_dict': self.model.state_dict(),
            'optimizer_state_dict': self.optimizer.state_dict(),
            'hyperparameters': {
                'train_losses': self.train_losses,
                'hidden_layers': self.hidden_layers,
                'hidden_nodes': self.hidden_nodes,
                'learning_rate': self.learning_rate,
                'num_epochs': self.num_epochs,
                'batch_size': self.batch_size,
                'weight_decay': self.weight_decay
            },
            'scalers': {
                'scaler_X': self.scaler_X,
                'scaler_Y': self.scaler_Y
            },
            'train_losses': self.train_losses
        }
        torch.save(model_state, title)
        print(f"Model saved to {title}")            
        
    def load_model(self, filepath):
        model_state = torch.load(filepath)

        # load hyperparameters
        hyperparameters = model_state['hyperparameters']
        self.hidden_layers = hyperparameters['hidden_layers']
        self.hidden_nodes = hyperparameters['hidden_nodes']
        self.learning_rate = hyperparameters['learning_rate']
        self.num_epochs = hyperparameters['num_epochs']
        self.batch_size = hyperparameters['batch_size']
        self.weight_decay = hyperparameters['weight_decay']

        # load connection weights of the model
        self.model.load_state_dict(model_state['model_state_dict'])
        self.optimizer.load_state_dict(model_state['optimizer_state_dict'])
        self.train_losses = model_state['train_losses']

        # load scalers
        self.scaler_X = model_state['scalers']['scaler_X']
        self.scaler_Y = model_state['scalers']['scaler_Y']

        print(f"Model loaded from {filepath}")
        
    
    def plot_learning_curves(self):
        plt.figure(figsize=(10,6))
        plt.plot(self.train_losses, label='Train Loss')
        plt.xlabel('Epoch')
        plt.ylabel('Loss')
        plt.title('Learning Curves')
        plt.legend()
        plt.grid(True)
        plt.show()

    def __del__(self):
        # 객체가 삭제될 때 메모리 정리
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

if __name__ == "__main__":
    # 데이터 로드
    data_loader = RawDataLoader()
    X_train, X_test, Y_train, Y_test, scaler_X, scaler_Y = data_loader.PrepareData()

    # 모델 생성 및 학습
    model = NNmodel(title='bmed_flux_batch_NN_model_v0.pth')
    
    # 스케일러 저장
    model.scaler_X = scaler_X
    model.scaler_Y = scaler_Y
    
    # 학습 시작 (인자 순서 수정)
    print("학습을 시작합니다...")
    test_loss = model.train(X_train, Y_train, X_test, Y_test)  # 순서 수정됨
    print(f"\n최종 테스트 손실: {test_loss:.4f}")
