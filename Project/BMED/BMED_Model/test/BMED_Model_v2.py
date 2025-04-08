import pandas as pd
import torch
import torch.nn as nn
import numpy as np
import matplotlib.pyplot as plt
import os
# Fix OpenMP runtime duplicate initialization error
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

df = pd.read_csv('BMED_data_v5.csv')
df1 = df[df['exp'] == df['exp'].unique()[0]]

X = df1[['T','V','E','Ci','Ki']].iloc[0]
S = df1[['NF_LA','NA_LA','NF_K','NB_K','VF','VA','VB']].iloc[0]
Y = df1[['NF_LA','NA_LA','NF_K','NB_K','VF','VA','VB']]
times = df1['t'].values

# Convert data to tensors
X = torch.FloatTensor(X.values)
S = torch.FloatTensor(S.values)
Y = torch.FloatTensor(Y.values)

# Set sequence length to Y data length
seq_length = Y.shape[0]

# Prepare training data
X_batch = X.unsqueeze(0).repeat(1, seq_length, 1)  # [1, seq_length, x_size]
S_batch = S.unsqueeze(0)  # [1, s_size]
Y_batch = Y.unsqueeze(0)  # [1, seq_length, s_size]

class BMEDLSTM(nn.Module):
    def __init__(self, x_size, s_size, hidden_size, num_layers=2):
        super(BMEDLSTM, self).__init__()
        self.hidden_size = hidden_size
        self.num_layers = num_layers
        self.s_size = s_size
        
        self.lstm = nn.LSTM(x_size + s_size, hidden_size, num_layers, batch_first=True)
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
            combined_input = torch.cat([x[:, 0:1], current_s], dim=2)  # 모든 시간에서 같은 x 사용
            
            # LSTM 순전파
            lstm_out, (h0, c0) = self.lstm(combined_input, (h0, c0))
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

# Model parameter settings
x_size = X.shape[0]
s_size = S.shape[0]
hidden_size = 64
dt = 0.01
max_time = times[-1]

# Initialize model
model = BMEDLSTM(x_size, s_size, hidden_size)

# Define loss function and optimizer
criterion = nn.MSELoss()
optimizer = torch.optim.Adam(model.parameters(), lr=0.0001)

# Training
num_epochs = 100
for epoch in range(num_epochs):
    model.train()
    optimizer.zero_grad()
    
    # Forward pass
    outputs, cal_times = model(X_batch, S_batch, dt, max_time)
    
    # Find closest calculation results to actual data times
    loss = 0
    for i, t in enumerate(times):
        # Find index of closest time
        idx = torch.abs(cal_times - t).argmin()
        loss += criterion(outputs[:, idx:idx+1], Y_batch[:, i:i+1])
    
    # Backward pass
    loss.backward()
    optimizer.step()
    
    if (epoch + 1) % 10 == 0:
        print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.4f}')

# Prediction
model.eval()
with torch.no_grad():
    predictions, cal_times = model(X_batch, S_batch, dt, max_time)
    print("\nPrediction Results:")
    print(predictions[0])

# Visualization
plt.figure(figsize=(15, 10))
y_names = ['NF_LA', 'NA_LA', 'NF_K', 'NB_K', 'VF', 'VA', 'VB']  # Column names of Y data
for i in range(len(y_names)):
    plt.subplot(3, 3, i+1)
    # Plot experimental data as scatter points
    plt.scatter(times, Y_batch[0, :, i].numpy(), c='b', label='Experimental', alpha=0.6)
    # Plot model predictions as continuous line
    plt.plot(cal_times.numpy(), predictions[0, :, i].numpy(), 'r-', label='Model Prediction', linewidth=1.5)
    plt.title(f'{y_names[i]}')
    plt.xlabel('Time (t)')
    plt.ylabel('Value')
    plt.legend()
    plt.grid(True)

plt.tight_layout()
plt.show()






