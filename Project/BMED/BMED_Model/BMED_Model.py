import pandas as pd
import torch
import torch.nn as nn
import numpy as np
import matplotlib.pyplot as plt
import os

# OpenMP 런타임 중복 초기화 오류 해결
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

df = pd.read_csv('BMED_data_v5.csv')
df1 = df[df['exp'] == df['exp'].unique()[0]]

X = df1[['T','V','E','Ci','Ki']].iloc[0]
S = df1[['NF_LA','NA_LA','NF_K','NB_K','VF','VA','VB']].iloc[0]
Y = df1[['NF_LA','NA_LA','NF_K','NB_K','VF','VA','VB']]

# 데이터를 텐서로 변환
X = torch.FloatTensor(X.values)
S = torch.FloatTensor(S.values)
Y = torch.FloatTensor(Y.values)

# 시퀀스 길이를 Y 데이터의 길이로 설정
seq_length = Y.shape[0]

# 학습 데이터 준비
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
        
    def forward(self, x, s0, seq_length):
        batch_size = x.size(0)
        h0 = torch.zeros(self.num_layers, batch_size, self.hidden_size).to(x.device)
        c0 = torch.zeros(self.num_layers, batch_size, self.hidden_size).to(x.device)
        
        # 초기 상태 설정
        current_s = s0.unsqueeze(1)  # [batch_size, 1, s_size]
        outputs = []
        
        # 각 시퀀스 스텝마다
        for t in range(seq_length):
            # 현재 입력 x와 상태 s 결합
            combined_input = torch.cat([x[:, t:t+1], current_s], dim=2)
            
            # LSTM 순전파
            lstm_out, (h0, c0) = self.lstm(combined_input, (h0, c0))
            step_output = self.fc(lstm_out)
            
            # 현재 스텝의 출력을 상태에 더함
            current_s = current_s + step_output
            
            outputs.append(step_output)
        
        # 모든 출력을 시퀀스로 결합
        return torch.cat(outputs, dim=1)

# 모델 파라미터 설정
x_size = X.shape[0]
s_size = S.shape[0]
hidden_size = 64

# 모델 초기화
model = BMEDLSTM(x_size, s_size, hidden_size)

# 손실 함수와 옵티마이저 정의
criterion = nn.MSELoss()
optimizer = torch.optim.Adam(model.parameters(), lr=0.001)

# 학습
num_epochs = 100
for epoch in range(num_epochs):
    model.train()
    optimizer.zero_grad()
    
    # 순전파
    outputs = model(X_batch, S_batch, seq_length)
    loss = criterion(outputs, Y_batch)
    
    # 역전파
    loss.backward()
    optimizer.step()
    
    if (epoch + 1) % 10 == 0:
        print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.4f}')

# 예측
model.eval()
with torch.no_grad():
    predictions = model(X_batch, S_batch, seq_length)
    print("\n예측 결과:")
    print(predictions[0])

# 시간 데이터 준비
time = df1['t'].values

# Y 값들의 이름
y_names = ['NF_LA', 'NA_LA', 'NF_K', 'NB_K', 'VF', 'VA', 'VB']

# 예측값과 실제값 시각화
plt.figure(figsize=(15, 10))
for i in range(len(y_names)):
    plt.subplot(3, 3, i+1)
    plt.plot(time, Y_batch[0, :, i].numpy(), 'b-', label='실제값')
    plt.plot(time, predictions[0, :, i].numpy(), 'r--', label='예측값')
    plt.title(f'{y_names[i]}')
    plt.xlabel('시간 (t)')
    plt.ylabel('값')
    plt.legend()
    plt.grid(True)

plt.tight_layout()
plt.show()






