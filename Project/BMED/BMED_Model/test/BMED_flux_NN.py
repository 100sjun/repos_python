import torch
import torch.nn as nn
import torch.optim as optim
import pandas as pd
import numpy as np
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
import optuna
import matplotlib.pyplot as plt
from torch.utils.data import Dataset, DataLoader
import warnings
warnings.filterwarnings('ignore')

# Set random seed for reproducibility
torch.manual_seed(42)
np.random.seed(42)

# Check if CUDA is available
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
print(f"Using device: {device}")

# Constants
DT = 0.01  # Fixed time step

# Custom Dataset
class BMEDDataset(Dataset):
    def __init__(self, X, y, t, exp_ids):
        self.X = torch.FloatTensor(X)
        self.y = torch.FloatTensor(y)
        self.t = torch.FloatTensor(t)
        self.exp_ids = torch.LongTensor(exp_ids)
        
    def __len__(self):
        return len(self.X)
    
    def __getitem__(self, idx):
        return self.X[idx], self.y[idx], self.t[idx], self.exp_ids[idx]

def calculate_cumulative_output(model, X, t, exp_ids):
    """Calculate cumulative output for each experiment"""
    unique_exps = torch.unique(exp_ids)
    outputs = []
    
    for exp in unique_exps:
        exp_mask = exp_ids == exp
        exp_X = X[exp_mask]
        exp_t = t[exp_mask]
        
        # 시간 순으로 정렬
        t_order = torch.argsort(exp_t)
        exp_X = exp_X[t_order]
        exp_t = exp_t[t_order]
        
        # flux 계산
        flux = model(exp_X)
        
        # 누적값 계산 (0부터 현재 시점까지의 flux 적분)
        n_steps = torch.round(exp_t / DT).long()  # 각 시점까지 필요한 step 수
        cumulative = torch.zeros_like(flux)
        
        for i in range(len(exp_t)):
            if i == 0:
                steps = n_steps[i]
                cumulative[i] = flux[i] * steps * DT
            else:
                steps = n_steps[i] - n_steps[i-1]
                cumulative[i] = cumulative[i-1] + flux[i] * steps * DT
        
        outputs.append(cumulative)
    
    return torch.cat(outputs, dim=0)

def load_data():
    # Load data
    df = pd.read_csv('BMED_data_v5.csv')
    
    # Input features
    X = df[['T', 'V', 'E', 'Ci', 'Ki']].values
    
    # Output features
    y = df[['NF_LA', 'NA_LA', 'NF_K', 'NB_K', 'VF', 'VA', 'VB']].values
    
    # Time and experiment IDs
    t = np.round(df['t'].values * 10) / 10  # Round to 0.1
    exp_ids = df['exp'].values
    
    # Scale input features
    scaler = StandardScaler()
    X = scaler.fit_transform(X)
    
    return X, y, t, exp_ids, scaler

# Custom loss function with R2 score
class R2Loss(nn.Module):
    def __init__(self):
        super(R2Loss, self).__init__()
        
    def forward(self, pred, target):
        # Calculate mean of target values
        target_mean = torch.mean(target, dim=0)
        
        # Calculate total sum of squares (TSS)
        tss = torch.sum((target - target_mean) ** 2, dim=0)
        
        # Calculate residual sum of squares (RSS)
        rss = torch.sum((target - pred) ** 2, dim=0)
        
        # Calculate R2 score for each output
        r2 = 1 - (rss / tss)
        
        # Return negative mean R2 (to minimize)
        return -torch.mean(r2)

# Neural Network Model
class BMEDModel(nn.Module):
    def __init__(self, input_size, hidden_sizes, output_size, dropout_rate, activation_type='ReLU'):
        super(BMEDModel, self).__init__()
        layers = []
        
        # Activation function selection
        if activation_type == 'ReLU':
            activation = nn.ReLU()
        elif activation_type == 'LeakyReLU':
            activation = nn.LeakyReLU(0.1)
        else:  # ELU
            activation = nn.ELU()
        
        # Input layer
        layers.append(nn.Linear(input_size, hidden_sizes[0]))
        layers.append(activation)
        layers.append(nn.BatchNorm1d(hidden_sizes[0]))
        layers.append(nn.Dropout(dropout_rate))
        
        # Hidden layers
        for i in range(len(hidden_sizes)-1):
            layers.append(nn.Linear(hidden_sizes[i], hidden_sizes[i+1]))
            layers.append(activation)
            layers.append(nn.BatchNorm1d(hidden_sizes[i+1]))
            layers.append(nn.Dropout(dropout_rate))
        
        # Output layer
        layers.append(nn.Linear(hidden_sizes[-1], output_size))
        
        self.model = nn.Sequential(*layers)
        
    def forward(self, x):
        return self.model(x)

def objective(trial):
    # Hyperparameters to optimize
    n_layers = trial.suggest_int('n_layers', 2, 8)
    hidden_sizes = []
    for i in range(n_layers):
        hidden_sizes.append(trial.suggest_int(f'hidden_size_{i}', 32, 128))
    
    learning_rate = trial.suggest_loguniform('learning_rate', 1e-5, 1e-2)
    batch_size = trial.suggest_categorical('batch_size', [16, 32, 64, 128])
    epochs = trial.suggest_int('epochs', 50, 200)
    dropout_rate = trial.suggest_float('dropout_rate', 0.1, 0.5)
    patience = trial.suggest_int('patience', 5, 20)
    activation = trial.suggest_categorical('activation', ['ReLU', 'LeakyReLU', 'ELU'])
    
    # Load data
    X, y, t, exp_ids, _ = load_data()
    
    # 5-fold cross validation
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    cv_scores = []
    
    for fold, (train_idx, val_idx) in enumerate(kf.split(X)):
        # Split data
        X_train, X_val = X[train_idx], X[val_idx]
        y_train, y_val = y[train_idx], y[val_idx]
        t_train, t_val = t[train_idx], t[val_idx]
        exp_ids_train, exp_ids_val = exp_ids[train_idx], exp_ids[val_idx]
        
        # Create datasets
        train_dataset = BMEDDataset(X_train, y_train, t_train, exp_ids_train)
        val_dataset = BMEDDataset(X_val, y_val, t_val, exp_ids_val)
        
        train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
        val_loader = DataLoader(val_dataset, batch_size=batch_size)
        
        # Initialize model
        model = BMEDModel(
            input_size=5, 
            hidden_sizes=hidden_sizes, 
            output_size=7,
            dropout_rate=dropout_rate,
            activation_type=activation
        ).to(device)
        
        optimizer = optim.Adam(model.parameters(), lr=learning_rate)
        scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=5, verbose=True)
        criterion = R2Loss()
        
        # Training loop
        best_val_loss = float('inf')
        patience_counter = 0
        best_model_state = None
        
        for epoch in range(epochs):
            model.train()
            train_loss = 0
            for batch_X, batch_y, batch_t, batch_exp_ids in train_loader:
                batch_X, batch_y, batch_t, batch_exp_ids = batch_X.to(device), batch_y.to(device), batch_t.to(device), batch_exp_ids.to(device)
                
                optimizer.zero_grad()
                # Calculate flux
                flux = model(batch_X)
                # Calculate cumulative output
                output = calculate_cumulative_output(model, batch_X, batch_t, batch_exp_ids)
                
                loss = criterion(output, batch_y)
                loss.backward()
                optimizer.step()
                
                train_loss += loss.item()
            
            # Validation
            model.eval()
            val_loss = 0
            with torch.no_grad():
                for batch_X, batch_y, batch_t, batch_exp_ids in val_loader:
                    batch_X, batch_y, batch_t, batch_exp_ids = batch_X.to(device), batch_y.to(device), batch_t.to(device), batch_exp_ids.to(device)
                    output = calculate_cumulative_output(model, batch_X, batch_t, batch_exp_ids)
                    val_loss += criterion(output, batch_y).item()
            
            # Learning rate scheduling
            scheduler.step(val_loss)
            
            if val_loss < best_val_loss:
                best_val_loss = val_loss
                patience_counter = 0
                best_model_state = model.state_dict().copy()
            else:
                patience_counter += 1
                if patience_counter >= patience:
                    break
        
        # Load best model state
        model.load_state_dict(best_model_state)
        cv_scores.append(best_val_loss)
    
    return np.mean(cv_scores)

def train_model(hidden_sizes, learning_rate, batch_size, epochs, dropout_rate, activation):
    # Load data
    X, y, t, exp_ids, scaler = load_data()
    
    # Create datasets
    dataset = BMEDDataset(X, y, t, exp_ids)
    train_size = int(0.8 * len(dataset))
    test_size = len(dataset) - train_size
    train_dataset, test_dataset = torch.utils.data.random_split(dataset, [train_size, test_size])
    
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    test_loader = DataLoader(test_dataset, batch_size=batch_size)
    
    # Initialize model
    model = BMEDModel(
        input_size=5, 
        hidden_sizes=hidden_sizes, 
        output_size=7,
        dropout_rate=dropout_rate,
        activation_type=activation
    ).to(device)
    
    optimizer = optim.Adam(model.parameters(), lr=learning_rate)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min', factor=0.5, patience=5, verbose=True)
    criterion = R2Loss()
    
    # Training history
    train_losses = []
    test_losses = []
    best_test_loss = float('inf')
    best_model_state = None
    
    # Training loop
    for epoch in range(epochs):
        model.train()
        train_loss = 0
        for batch_X, batch_y, batch_t, batch_exp_ids in train_loader:
            batch_X, batch_y, batch_t, batch_exp_ids = batch_X.to(device), batch_y.to(device), batch_t.to(device), batch_exp_ids.to(device)
            
            optimizer.zero_grad()
            output = calculate_cumulative_output(model, batch_X, batch_t, batch_exp_ids)
            
            loss = criterion(output, batch_y)
            loss.backward()
            optimizer.step()
            
            train_loss += loss.item()
        
        # Test evaluation
        model.eval()
        test_loss = 0
        all_preds = []
        all_targets = []
        
        with torch.no_grad():
            for batch_X, batch_y, batch_t, batch_exp_ids in test_loader:
                batch_X, batch_y, batch_t, batch_exp_ids = batch_X.to(device), batch_y.to(device), batch_t.to(device), batch_exp_ids.to(device)
                output = calculate_cumulative_output(model, batch_X, batch_t, batch_exp_ids)
                test_loss += criterion(output, batch_y).item()
                
                all_preds.extend(output.cpu().numpy())
                all_targets.extend(batch_y.cpu().numpy())
        
        # Learning rate scheduling
        scheduler.step(test_loss)
        
        train_losses.append(train_loss / len(train_loader))
        test_losses.append(test_loss / len(test_loader))
        
        if test_loss < best_test_loss:
            best_test_loss = test_loss
            best_model_state = model.state_dict().copy()
        
        if (epoch + 1) % 10 == 0:
            print(f'Epoch [{epoch+1}/{epochs}], Train Loss: {train_loss/len(train_loader):.4f}, Test Loss: {test_loss/len(test_loader):.4f}')
    
    # Load best model state
    model.load_state_dict(best_model_state)
    return model, train_losses, test_losses, all_preds, all_targets

def plot_results(train_losses, test_losses, all_preds, all_targets):
    # Plot training curves
    plt.figure(figsize=(15, 5))
    
    # Loss curves
    plt.subplot(1, 2, 1)
    plt.plot(train_losses, label='Train Loss')
    plt.plot(test_losses, label='Test Loss')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.title('Training and Test Loss Over Time')
    plt.legend()
    plt.grid(True)
    
    # Scatter plot of predictions vs targets
    plt.subplot(1, 2, 2)
    all_preds = np.array(all_preds)
    all_targets = np.array(all_targets)
    plt.scatter(all_targets, all_preds, alpha=0.5)
    plt.plot([all_targets.min(), all_targets.max()], [all_targets.min(), all_targets.max()], 'r--')
    plt.xlabel('True Values')
    plt.ylabel('Predicted Values')
    plt.title('Predicted vs True Values')
    plt.grid(True)
    
    plt.tight_layout()
    plt.savefig('training_results.png')
    plt.close()

def main():
    # Optuna optimization
    study = optuna.create_study(direction='minimize')
    study.optimize(objective, n_trials=50)
    
    print("Best trial:")
    trial = study.best_trial
    print("  Value: ", trial.value)
    print("  Params: ")
    for key, value in trial.params.items():
        print(f"    {key}: {value}")
    
    # Extract best hyperparameters
    n_layers = trial.params['n_layers']
    hidden_sizes = [trial.params[f'hidden_size_{i}'] for i in range(n_layers)]
    learning_rate = trial.params['learning_rate']
    batch_size = trial.params['batch_size']
    epochs = trial.params['epochs']
    dropout_rate = trial.params['dropout_rate']
    activation = trial.params['activation']
    
    # Train final model with 10x epochs
    model, train_losses, test_losses, all_preds, all_targets = train_model(
        hidden_sizes=hidden_sizes,
        learning_rate=learning_rate,
        batch_size=batch_size,
        epochs=epochs * 10,
        dropout_rate=dropout_rate,
        activation=activation
    )
    
    # Plot results
    plot_results(train_losses, test_losses, all_preds, all_targets)
    
    # Save model
    torch.save(model.state_dict(), 'bmed_model.pth')

if __name__ == "__main__":
    main()
