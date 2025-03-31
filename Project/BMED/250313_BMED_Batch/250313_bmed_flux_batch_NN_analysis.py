import pandas as pd
import numpy as np
import torch
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from sklearn.preprocessing import StandardScaler
from bmed_flux_batch_NN_train import NNmodel, RawDataLoader

# StandardScaler를 안전한 전역 변수로 등록
torch.serialization.add_safe_globals([StandardScaler])

def analyze_model(model_path, data_path='250314_BMED_train data_v4.xlsx'):
    # 데이터 로드
    data_loader = RawDataLoader(path=data_path)
    X_data = data_loader.X_data
    Y_data = data_loader.Y_data
    
    # 모델 로드
    model = NNmodel()
    model.load_model(model_path, weights_only=False)
    
    # 데이터 정규화
    X_scaled = model.scaler_X.transform(X_data)
    X_tensor = torch.FloatTensor(X_scaled).to(model.device)
    
    # 예측 수행
    with torch.no_grad():
        predictions = model.model(X_tensor)
        predictions = predictions.cpu().numpy()
    
    # 예측값 역정규화
    predictions_original = model.scaler_Y.inverse_transform(predictions)
    
    # 성능 지표 계산
    output_names = ['dNLA', 'dNK', 'dVF', 'dVA', 'dVB']
    performance_metrics = {}
    
    for i, name in enumerate(output_names):
        r2 = r2_score(Y_data[:, i], predictions_original[:, i])
        mse = mean_squared_error(Y_data[:, i], predictions_original[:, i])
        rmse = np.sqrt(mse)
        mae = mean_absolute_error(Y_data[:, i], predictions_original[:, i])
        
        performance_metrics[name] = {
            'R2': r2,
            'RMSE': rmse,
            'MAE': mae
        }
    
    # 결과 시각화
    plot_results(Y_data, predictions_original, output_names, performance_metrics)
    
    # 성능 지표 출력
    print("\n=== 성능 지표 ===")
    for name, metrics in performance_metrics.items():
        print(f"\n{name}:")
        for metric_name, value in metrics.items():
            print(f"{metric_name}: {value:.4f}")
    
    return performance_metrics

def plot_results(Y_true, Y_pred, output_names, performance_metrics):
    for i, name in enumerate(output_names):
        plt.figure(figsize=(10, 6))
        
        # 산점도
        plt.scatter(Y_true[:, i], Y_pred[:, i], alpha=0.5, 
                   label=f'R² = {performance_metrics[name]["R2"]:.3f}')
        
        # 대각선
        min_val = min(min(Y_true[:, i]), min(Y_pred[:, i]))
        max_val = max(max(Y_true[:, i]), max(Y_pred[:, i]))
        plt.plot([min_val, max_val], [min_val, max_val], 'r--')
        
        plt.xlabel(f'실제 {name}')
        plt.ylabel(f'예측 {name}')
        plt.title(f'{name} 예측 결과')
        plt.legend()
        plt.grid(True)
        plt.savefig(f'analysis_plot_{name}.png')
        plt.close()

if __name__ == "__main__":
    model_path = 'bmed_flux_batch_NN_model_v0.pth'
    performance_metrics = analyze_model(model_path) 