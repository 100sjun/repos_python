# 전체 실험 데이터를 시각화하는 새로운 코드
# 8개의 개별 plot을 생성하여 각 plot은 모든 실험의 해당 feature를 subplot으로 표시

import matplotlib.pyplot as plt
import numpy as np
import torch
import pandas as pd
import math

def comprehensive_model_evaluation(evaluation_model, dataloader, ndf, range_mm, model_info, best_train_loss=None, best_epoch=None):
    """
    전체 실험 데이터에 대한 종합적인 모델 평가 시각화
    8개 feature별로 개별 figure를 생성하고, 각 figure는 모든 실험의 해당 feature를 subplot으로 표시
    
    Args:
        evaluation_model: 평가할 모델
        dataloader: 전체 데이터로더
        ndf: 정규화된 데이터프레임
        range_mm: 정규화 범위 정보
        model_info: 모델 정보 문자열
        best_train_loss: 최적 훈련 손실 (선택적)
        best_epoch: 최적 에포크 (선택적)
    """
    
    # Device 설정
    device = next(evaluation_model.parameters()).device
    
    def denormalize_data(normalized_data, feature_name):
        """정규화된 데이터를 원래 단위로 복원"""
        if feature_name in range_mm:
            min_val = range_mm[feature_name]['min']
            max_val = range_mm[feature_name]['max']
            return normalized_data * (max_val - min_val) + min_val
        return normalized_data

    # feature 정보 정의
    feature_names = ['VF', 'VA', 'VB', 'CFLA', 'CALA', 'CFK', 'CBK', 'I']
    feature_indices = [2, 3, 4, 5, 6, 7, 8, 9]
    feature_units = ['L', 'L', 'L', 'mol/L', 'mol/L', 'mol/L', 'mol/L', 'A']
    feature_titles = ['Feed Volume', 'Acid Volume', 'Base Volume', 'Feed LA Conc.', 'Acid LA Conc.', 'Feed K Conc.', 'Base K Conc.', 'Current']
    
    # 실험 번호와 데이터로더 인덱스 매핑 생성
    exp_num_to_dataloader_idx = {}
    dataloader_idx_to_exp_num = {}
    
    for i, exp_num in enumerate(sorted(ndf['exp'].unique())):
        exp_num_to_dataloader_idx[exp_num] = i
        dataloader_idx_to_exp_num[i] = exp_num
    
    available_exp_nums = sorted(ndf['exp'].unique())
    total_experiments = len(available_exp_nums)
    print(f"🔍 전체 {total_experiments}개 실험 데이터 분석 시작...")
    
    # 전체 실험 데이터 수집
    all_experiment_data = {}
    
    evaluation_model.eval()
    with torch.no_grad():
        for dataloader_idx, (input_seq, seq_lengths) in enumerate(dataloader):
            try:
                # 실제 실험 번호 확인
                actual_exp_num = dataloader_idx_to_exp_num[dataloader_idx]
                
                input_seq = input_seq.to(device)
                seq_lengths = seq_lengths.to(device)
                
                # Teacher forcing 데이터 준비 (기존 tf_data 함수 사용)
                inputs = input_seq[:, :-1, :-1]  # 전류 제외한 입력
                targets = input_seq[:, 1:, :]    # 전류 포함한 타겟
                target_seq_lengths = seq_lengths - 1
                
                # 모델 예측
                predictions = evaluation_model(inputs, target_seq_lengths)
                
                # CPU로 이동
                predictions_cpu = predictions.cpu().numpy()
                targets_cpu = targets.cpu().numpy()
                input_seq_cpu = input_seq.cpu().numpy()
                
                # 첫 번째 샘플 사용 (batch size = 4이지만 같은 실험 데이터)
                pred_sample = predictions_cpu[0]
                actual_sample = targets_cpu[0]
                initial_sample = input_seq_cpu[0]
                seq_len = target_seq_lengths[0].item()
                full_seq_len = seq_lengths[0].item()
                
                # 실제 시퀀스 길이만큼만 사용
                pred_sample = pred_sample[:seq_len]
                actual_sample = actual_sample[:seq_len]
                initial_sample = initial_sample[:full_seq_len]
                
                # 시간 축 생성
                time_points = np.arange(full_seq_len) * 0.25
                
                # 예측값을 초기값과 연결하여 전체 시퀀스 구성
                prediction_full = np.zeros((full_seq_len, 10))
                prediction_full[0] = initial_sample[0]  # 초기값
                prediction_full[1:seq_len+1] = pred_sample  # 예측값
                
                # 실험 조건 정보
                voltage = denormalize_data(initial_sample[0, 0], 'V')
                electrolyte = denormalize_data(initial_sample[0, 1], 'E')
                
                # 데이터 저장
                all_experiment_data[actual_exp_num] = {
                    'dataloader_idx': dataloader_idx,
                    'time_points': time_points,
                    'actual_data': initial_sample,
                    'predicted_data': prediction_full,
                    'seq_len': seq_len,
                    'full_seq_len': full_seq_len,
                    'voltage': voltage,
                    'electrolyte': electrolyte,
                    'actual_sample': actual_sample,
                    'pred_sample': pred_sample
                }
                
                print(f"✅ 실험 {actual_exp_num} 데이터 처리 완료 (시간: {full_seq_len*0.25:.2f}h, V: {voltage:.1f}V, E: {electrolyte:.3f}M)")
                
            except Exception as e:
                print(f"❌ 실험 {dataloader_idx} 처리 중 오류: {str(e)}")
                continue
    
    print(f"📊 총 {len(all_experiment_data)}개 실험 데이터 수집 완료!")
    
    # 서브플롯 그리드 크기 계산 (24개 실험을 위한 최적 배치)
    sqrt_exp = math.ceil(math.sqrt(total_experiments))
    if sqrt_exp * (sqrt_exp - 1) >= total_experiments:
        grid_rows = sqrt_exp - 1
        grid_cols = sqrt_exp
    else:
        grid_rows = sqrt_exp
        grid_cols = sqrt_exp
    
    print(f"📐 서브플롯 그리드: {grid_rows}x{grid_cols} ({grid_rows * grid_cols}개 위치)")
    
    # 8개 feature별로 개별 figure 생성
    for feature_idx, (feature_name, feature_index, unit, title) in enumerate(zip(feature_names, feature_indices, feature_units, feature_titles)):
        print(f"🎨 {feature_name} ({title}) 플롯 생성 중...")
        
        # 새 figure 생성
        fig, axes = plt.subplots(grid_rows, grid_cols, figsize=(4*grid_cols, 3*grid_rows))
        fig.suptitle(f'{title} ({feature_name}) - All Experiments\n{model_info}', fontsize=16, fontweight='bold')
        
        # axes를 1차원 배열로 변환 (단일 subplot인 경우 처리)
        if grid_rows == 1 and grid_cols == 1:
            axes = [axes]
        elif grid_rows == 1 or grid_cols == 1:
            axes = axes.flatten()
        else:
            axes = axes.flatten()
        
        # 전체 성능 지표 수집
        all_rmse = []
        all_r2 = []
        
        # 각 실험별 subplot 생성
        for plot_idx, exp_num in enumerate(sorted(all_experiment_data.keys())):
            if plot_idx >= len(axes):
                break
                
            ax = axes[plot_idx]
            exp_data = all_experiment_data[exp_num]
            
            # 데이터 denormalize
            actual_denorm = denormalize_data(exp_data['actual_data'][:, feature_index], feature_name)
            pred_denorm = denormalize_data(exp_data['predicted_data'][:exp_data['seq_len']+1, feature_index], feature_name)
            
            # 실제 vs 예측 플롯
            ax.plot(exp_data['time_points'], actual_denorm, 'b-', linewidth=1.5, label='Actual', marker='o', markersize=2)
            ax.plot(exp_data['time_points'][:exp_data['seq_len']+1], pred_denorm, 'r--', linewidth=1.5, label='Predicted', marker='s', markersize=2)
            
            # 초기값과 예측 시작점 표시
            ax.plot(exp_data['time_points'][0], actual_denorm[0], 'go', markersize=4)
            ax.axvline(x=exp_data['time_points'][1], color='gray', linestyle=':', alpha=0.5)
            
            # 성능 지표 계산 (예측 구간만)
            actual_pred_denorm = denormalize_data(exp_data['actual_sample'][:, feature_index], feature_name)
            pred_pred_denorm = denormalize_data(exp_data['pred_sample'][:, feature_index], feature_name)
            
            rmse = np.sqrt(np.mean((actual_pred_denorm - pred_pred_denorm)**2))
            r2 = 1 - (np.sum((actual_pred_denorm - pred_pred_denorm)**2) / np.sum((actual_pred_denorm - np.mean(actual_pred_denorm))**2))
            
            all_rmse.append(rmse)
            all_r2.append(r2)
            
            # 서브플롯 제목과 레이블 설정
            ax.set_title(f'Exp {exp_num}\nV:{exp_data["voltage"]:.1f}V, E:{exp_data["electrolyte"]:.2f}M\nRMSE:{rmse:.3f}', fontsize=8)
            ax.set_xlabel('Time (h)', fontsize=8)
            ax.set_ylabel(f'{feature_name} ({unit})', fontsize=8)
            ax.grid(True, alpha=0.3)
            ax.tick_params(labelsize=7)
            
            # 범례는 첫 번째 subplot에만 표시
            if plot_idx == 0:
                ax.legend(fontsize=7, loc='best')
        
        # 사용되지 않는 subplot 숨기기
        for plot_idx in range(len(all_experiment_data), len(axes)):
            axes[plot_idx].set_visible(False)
        
        # 전체 성능 요약 텍스트 추가
        avg_rmse = np.mean(all_rmse)
        avg_r2 = np.mean(all_r2)
        std_rmse = np.std(all_rmse)
        std_r2 = np.std(all_r2)
        
        summary_text = f"Overall Performance:\nRMSE: {avg_rmse:.4f}±{std_rmse:.4f} {unit}\nR²: {avg_r2:.3f}±{std_r2:.3f}"
        fig.text(0.02, 0.02, summary_text, fontsize=10, bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
        
        plt.tight_layout()
        plt.subplots_adjust(top=0.93, bottom=0.1)
        plt.show()
        
        print(f"✅ {feature_name} 플롯 완료 - 평균 RMSE: {avg_rmse:.4f} {unit}, 평균 R²: {avg_r2:.3f}")
    
    # 전체 요약 통계
    print("\n" + "="*80)
    print("COMPREHENSIVE MODEL EVALUATION SUMMARY")
    print(f"{model_info}")
    print("="*80)
    
    # Feature별 전체 성능 요약
    overall_performance = {}
    for feature_idx, (feature_name, feature_index, unit) in enumerate(zip(feature_names, feature_indices, feature_units)):
        all_rmse = []
        all_r2 = []
        all_mae = []
        
        for exp_num in sorted(all_experiment_data.keys()):
            exp_data = all_experiment_data[exp_num]
            
            # 성능 지표 계산
            actual_pred_denorm = denormalize_data(exp_data['actual_sample'][:, feature_index], feature_name)
            pred_pred_denorm = denormalize_data(exp_data['pred_sample'][:, feature_index], feature_name)
            
            rmse = np.sqrt(np.mean((actual_pred_denorm - pred_pred_denorm)**2))
            mae = np.mean(np.abs(actual_pred_denorm - pred_pred_denorm))
            r2 = 1 - (np.sum((actual_pred_denorm - pred_pred_denorm)**2) / np.sum((actual_pred_denorm - np.mean(actual_pred_denorm))**2))
            
            all_rmse.append(rmse)
            all_mae.append(mae)
            all_r2.append(r2)
        
        overall_performance[feature_name] = {
            'rmse_mean': np.mean(all_rmse),
            'rmse_std': np.std(all_rmse),
            'mae_mean': np.mean(all_mae),
            'mae_std': np.std(all_mae),
            'r2_mean': np.mean(all_r2),
            'r2_std': np.std(all_r2),
            'unit': unit
        }
        
        print(f"{feature_name:6s}: RMSE={np.mean(all_rmse):.4f}±{np.std(all_rmse):.4f} {unit}, "
              f"MAE={np.mean(all_mae):.4f}±{np.std(all_mae):.4f} {unit}, "
              f"R²={np.mean(all_r2):.3f}±{np.std(all_r2):.3f}")
    
    # 전체 평균 성능
    all_features_rmse = [overall_performance[fname]['rmse_mean'] for fname in feature_names]
    all_features_r2 = [overall_performance[fname]['r2_mean'] for fname in feature_names]
    
    print("-"*80)
    print(f"OVERALL AVERAGE:")
    print(f"  Average RMSE across features: {np.mean(all_features_rmse):.4f} (mixed units)")
    print(f"  Average R² across features: {np.mean(all_features_r2):.3f}")
    print(f"  Total experiments evaluated: {len(all_experiment_data)}")
    
    if best_train_loss is not None and best_epoch is not None:
        print(f"  Best training loss: {best_train_loss:.6f} (Epoch {best_epoch})")
    
    print("="*80)
    
    return overall_performance, all_experiment_data

# 사용 예시:
# comprehensive_model_evaluation(best_model, dataloader, ndf, range_mm, model_info, best_train_loss, best_epoch)