import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import matplotlib.pyplot as plt

# 데이터 생성 함수
def generate_data(num_samples=1000):
    x1 = np.random.uniform(0, 10, num_samples)
    x2 = np.random.uniform(0, 5, num_samples)
    x3 = np.random.uniform(0, 5, num_samples)
    x4 = np.random.uniform(0, 5, num_samples)
    x5 = np.random.uniform(0, 5, num_samples)
    
    # a1, a2, a3를 x2, x3, x4, x5의 함수로 정의 (예시)
    a1 = 2 + 0.5*x2 + 0.3*x3
    a2 = 0.5 + 0.1*x3 + 0.2*x4
    a3 = 1 + 0.2*x4 + 0.4*x5
    
    y = a1 * np.exp(-a2 * x1) + a3 + np.random.normal(0, 0.1, num_samples)
    
    return x1, x2, x3, x4, x5, y

# 신경망 모델 정의
class ComplexExponentialModel(nn.Module):
    def __init__(self):
        super().__init__()
        self.hidden = nn.Sequential(
            nn.Linear(4, 64),
            nn.ReLU(),
            nn.Linear(64, 32),
            nn.ReLU(),
            nn.Linear(32, 3)
        )

    def forward(self, x):
        x1 = x[:, 0].unsqueeze(1)
        x_rest = x[:, 1:]
        
        a1, a2, a3 = self.hidden(x_rest).chunk(3, dim=1)
        
        # a1과 a2는 항상 양수여야 하므로 exp 함수 사용
        a1 = torch.exp(a1)
        a2 = torch.exp(a2)
        
        return a1 * torch.exp(-a2 * x1) + a3

# 데이터 생성 및 전처리
x1, x2, x3, x4, x5, y = generate_data()
X = np.column_stack((x1, x2, x3, x4, x5))
X_tensor = torch.FloatTensor(X)
y_tensor = torch.FloatTensor(y).unsqueeze(1)

# 모델, 손실 함수, 옵티마이저 초기화
model = ComplexExponentialModel()
criterion = nn.MSELoss()
optimizer = optim.Adam(model.parameters(), lr=0.001)

# 학습 루프
num_epochs = 1000
losses = []

for epoch in range(num_epochs):
    # 순전파
    y_pred = model(X_tensor)
    loss = criterion(y_pred, y_tensor)
    
    # 역전파
    optimizer.zero_grad()
    loss.backward()
    optimizer.step()
    
    losses.append(loss.item())
    
    if (epoch + 1) % 100 == 0:
        print(f'Epoch [{epoch+1}/{num_epochs}], Loss: {loss.item():.4f}')

# 결과 시각화
plt.figure(figsize=(10, 6))
plt.scatter(x1, y, alpha=0.5, label='Data')
y_pred = model(X_tensor).detach().numpy()
plt.scatter(x1, y_pred, color='red', alpha=0.5, label='Predictions')
plt.legend()
plt.xlabel('x1')
plt.ylabel('y')
plt.title('Complex Exponential Function Fitting')
plt.show()

# 학습 곡선 시각화
plt.figure(figsize=(10, 6))
plt.plot(losses)
plt.title('Model Loss')
plt.ylabel('Loss')
plt.xlabel('Epoch')
plt.show()

# 몇 가지 샘플에 대한 a1, a2, a3 값 출력
sample_X = X_tensor[:5]
sample_outputs = model.hidden(sample_X[:, 1:])
a1, a2, a3 = [torch.exp(param) if i < 2 else param for i, param in enumerate(sample_outputs.chunk(3, dim=1))]

print("\n샘플 데이터에 대한 추정된 a1, a2, a3 값:")
for i in range(5):
    print(f"Sample {i+1}: a1 = {a1[i].item():.4f}, a2 = {a2[i].item():.4f}, a3 = {a3[i].item():.4f}")