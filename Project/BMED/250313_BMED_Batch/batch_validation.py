import pandas as pd
import numpy as np
import torch
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error
from bmed_flux_batch_NN_train import NNmodel, RawDataLoader

def validate_by_experiment(model_path, data_path='250314_BMED_train data_v4.xlsx'):
    # 데이터 로드
    data_loader = RawDataLoader(path=data_path)
    df = pd.read_excel(data_path)
    
    # 모델 로드
    model = NNmodel()
    model.load_model(model_path, weights_only=False)
    
    # 실험 번호 목록 추출
    exp_numbers = df['exp'].unique()
    output_names = ['dNLA', 'dNK', 'dVF', 'dVA', 'dVB']
    
    # 실험별 성능 저장을 위한 딕셔너리
    exp_metrics = {}
    
    for exp_num in exp_numbers:
        # 실험별 데이터 추출
        exp_mask = df['exp'] == exp_num
        X_exp = data_loader.X_data[exp_mask]
        Y_exp = data_loader.Y_data[exp_mask]
        
        # 데이터 정규화 및 예측
        X_scaled = model.scaler_X.transform(X_exp)
        X_tensor = torch.FloatTensor(X_scaled).to(model.device)
        
        with torch.no_grad():
            predictions = model.model(X_tensor)
            predictions = predictions.cpu().numpy()
        
        # 예측값 역정규화
        predictions_original = model.scaler_Y.inverse_transform(predictions)
        
        # 실험별 성능 지표 계산
        exp_metrics[exp_num] = {}
        for i, name in enumerate(output_names):
            r2 = r2_score(Y_exp[:, i], predictions_original[:, i])
            rmse = np.sqrt(mean_squared_error(Y_exp[:, i], predictions_original[:, i]))
            mae = mean_absolute_error(Y_exp[:, i], predictions_original[:, i])
            
            exp_metrics[exp_num][name] = {
                'R2': r2,
                'RMSE': rmse,
                'MAE': mae
            }
        
        # 실험별 결과 시각화
        plot_exp_results(Y_exp, predictions_original, output_names, exp_num, exp_metrics[exp_num])
    
    # 전체 실험 결과 요약
    print_summary(exp_metrics, exp_numbers, output_names)
    plot_summary(exp_metrics, exp_numbers, output_names)
    
    return exp_metrics

def plot_exp_results(Y_true, Y_pred, output_names, exp_num, metrics):
    for i, name in enumerate(output_names):
        plt.figure(figsize=(8, 6))
        
        plt.scatter(Y_true[:, i], Y_pred[:, i], alpha=0.6,
                   label=f'R² = {metrics[name]["R2"]:.3f}')
        
        min_val = min(min(Y_true[:, i]), min(Y_pred[:, i]))
        max_val = max(max(Y_true[:, i]), max(Y_pred[:, i]))
        plt.plot([min_val, max_val], [min_val, max_val], 'r--')
        
        plt.xlabel(f'실제 {name}')
        plt.ylabel(f'예측 {name}')
        plt.title(f'실험 {exp_num} - {name} 예측 결과')
        plt.legend()
        plt.grid(True)
        plt.savefig(f'exp_{exp_num}_{name}_validation.png')
        plt.close()

def print_summary(exp_metrics, exp_numbers, output_names):
    print("\n=== 실험별 성능 지표 요약 ===")
    
    # 각 출력 변수별로 실험 간 성능 비교
    for name in output_names:
        print(f"\n{name} 결과:")
        print("실험번호    R²      RMSE     MAE")
        print("-" * 35)
        
        for exp_num in exp_numbers:
            metrics = exp_metrics[exp_num][name]
            print(f"{exp_num:^8} {metrics['R2']:6.3f} {metrics['RMSE']:8.3f} {metrics['MAE']:8.3f}")

def plot_summary(exp_metrics, exp_numbers, output_names):
    # 출력 변수별 R² 값 비교 그래프
    plt.figure(figsize=(12, 6))
    x = np.arange(len(exp_numbers))
    width = 0.15
    
    for i, name in enumerate(output_names):
        r2_values = [exp_metrics[exp]['name']['R2'] for exp in exp_numbers]
        plt.bar(x + i*width, r2_values, width, label=name)
    
    plt.xlabel('실험 번호')
    plt.ylabel('R² Score')
    plt.title('실험별 R² 점수 비교')
    plt.xticks(x + width*2, exp_numbers)
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.savefig('exp_comparison_summary.png')
    plt.close()

if __name__ == "__main__":
    model_path = 'bmed_flux_batch_NN_model_v0.pth'
    exp_metrics = validate_by_experiment(model_path)
