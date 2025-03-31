import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import matplotlib.pyplot as plt


# 간단한 신경망 모델 정의
class SimpleNN(nn.Module):
    def __init__(self):
        super(SimpleNN, self).__init__()
        self.model = nn.Sequential(
            nn.Linear(1, 64),
            nn.ReLU(),
            nn.Linear(64, 64),
            nn.ReLU(),
            nn.Linear(64, 1)
        )
    
    def forward(self, x):
        return self.model(x)


# 예제 데이터 생성
def generate_data(samples=100):
    x = np.linspace(-5, 5, samples).reshape(-1, 1)
    y = 0.2 * x**3 + 0.5 * x**2 - x + 2 + np.random.randn(samples, 1) * 0.5
    return torch.FloatTensor(x), torch.FloatTensor(y)


# 모델 학습
def train_model(model, x, y, epochs=1000):
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=0.01)
    
    losses = []
    
    for epoch in range(epochs):
        optimizer.zero_grad()
        outputs = model(x)
        loss = criterion(outputs, y)
        loss.backward()
        optimizer.step()
        
        losses.append(loss.item())
        
        if epoch % 100 == 0:
            print(f"Epoch {epoch}/{epochs}, Loss: {loss.item():.4f}")
    
    return losses


# 실행 함수
def main():
    # 장치 설정
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")
    
    # 데이터 생성
    x, y = generate_data(200)
    x, y = x.to(device), y.to(device)
    
    # 모델 초기화
    model = SimpleNN().to(device)
    
    # 모델 학습
    losses = train_model(model, x, y, epochs=1000)
    
    # 결과 시각화
    plt.figure(figsize=(12, 5))
    
    # 손실 그래프
    plt.subplot(1, 2, 1)
    plt.plot(losses)
    plt.title('Training Loss')
    plt.xlabel('Epoch')
    plt.ylabel('MSE Loss')
    
    # 예측 결과
    plt.subplot(1, 2, 2)
    
    # CPU로 데이터 이동
    x_cpu = x.cpu().detach().numpy()
    y_cpu = y.cpu().detach().numpy()
    
    # 정렬 (시각화를 위해)
    idx = np.argsort(x_cpu.flatten())
    x_sorted = x_cpu[idx]
    y_sorted = y_cpu[idx]
    
    # 예측 수행
    with torch.no_grad():
        pred = model(x).cpu().numpy()
    pred_sorted = pred[idx]
    
    plt.scatter(x_sorted, y_sorted, label='Actual Data', alpha=0.5)
    plt.plot(x_sorted, pred_sorted, 'r-', linewidth=2, label='Model Prediction')
    plt.title('Regression Result')
    plt.xlabel('Input (x)')
    plt.ylabel('Output (y)')
    plt.legend()
    
    plt.tight_layout()
    plt.savefig('regression_result.png')
    
    # Docker 컨테이너에서는 show() 호출하지 않음
    # plt.show() 대신 파일로만 저장
    
    print("Training completed! Results saved to regression_result.png.")


if __name__ == "__main__":
    main() 
