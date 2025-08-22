# 기존 notebook에 추가할 수 있는 comprehensive evaluation 코드
# 기존 변수들 (best_model, dataloader, ndf, range_mm 등)을 사용하여 전체 실험 데이터 시각화

# 모든 실험 데이터에 대한 comprehensive evaluation
print("🚀 전체 실험 데이터 comprehensive evaluation 시작...")

# Device 확인
device = next(best_model.parameters()).device if best_model is not None else next(model.parameters()).device
evaluation_model = best_model if best_model is not None else model

# feature 정보 재정의
feature_names = ['VF', 'VA', 'VB', 'CFLA', 'CALA', 'CFK', 'CBK', 'I']
feature_indices = [2, 3, 4, 5, 6, 7, 8, 9]
feature_units = ['L', 'L', 'L', 'mol/L', 'mol/L', 'mol/L', 'mol/L', 'A']
feature_titles = ['Feed Volume', 'Acid Volume', 'Base Volume', 'Feed LA Conc.', 'Acid LA Conc.', 'Feed K Conc.', 'Base K Conc.', 'Current']

# 전체 실험 데이터 수집
all_experiment_results = {}
available_exp_nums = sorted(ndf['exp'].unique())
total_experiments = len(available_exp_nums)

# 실험 번호 매핑 (기존 코드에서 가져옴)
exp_num_to_dataloader_idx = {}
for i, exp_num in enumerate(available_exp_nums):
    exp_num_to_dataloader_idx[exp_num] = i

print(f"📋 전체 {total_experiments}개 실험 데이터 수집 중...")

evaluation_model.eval()
with torch.no_grad():
    for dataloader_idx, (input_seq, seq_lengths) in enumerate(dataloader):
        try:
            # 해당 데이터로더 인덱스의 실험 번호 찾기
            actual_exp_num = available_exp_nums[dataloader_idx]
            
            input_seq = input_seq.to(device)
            seq_lengths = seq_lengths.to(device)
            
            # Teacher forcing 데이터 준비 (기존 함수 재사용)
            inputs, targets, target_seq_lengths = prepare_teacher_forcing_data(input_seq, seq_lengths)
            
            # 모델 예측
            predictions = evaluation_model(inputs, target_seq_lengths)
            
            # CPU로 이동
            predictions_cpu = predictions.cpu().numpy()
            targets_cpu = targets.cpu().numpy()
            input_seq_cpu = input_seq.cpu().numpy()
            
            # 첫 번째 샘플 사용
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
            
            # 예측값을 초기값과 연결
            prediction_full = np.zeros((full_seq_len, 10))
            prediction_full[0] = initial_sample[0]
            prediction_full[1:seq_len+1] = pred_sample
            
            # 실험 조건 정보
            voltage = denormalize_data(initial_sample[0, 0], 'V')
            electrolyte = denormalize_data(initial_sample[0, 1], 'E')
            
            # 결과 저장
            all_experiment_results[actual_exp_num] = {
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
            
            print(f"✅ 실험 {actual_exp_num}: {full_seq_len*0.25:.1f}h, {voltage:.1f}V, {electrolyte:.3f}M")
            
        except Exception as e:
            print(f"❌ 실험 {dataloader_idx} 오류: {str(e)}")
            continue

print(f"📊 총 {len(all_experiment_results)}개 실험 데이터 수집 완료!")

# 서브플롯 배치 계산
import math
sqrt_exp = math.ceil(math.sqrt(total_experiments))
if sqrt_exp * (sqrt_exp - 1) >= total_experiments:
    grid_rows, grid_cols = sqrt_exp - 1, sqrt_exp
else:
    grid_rows, grid_cols = sqrt_exp, sqrt_exp

print(f"📐 각 feature당 {grid_rows}x{grid_cols} 서브플롯 그리드 사용")

# 8개 feature별로 개별 figure 생성
for feature_idx, (feature_name, feature_index, unit, title) in enumerate(zip(feature_names, feature_indices, feature_units, feature_titles)):
    print(f"🎨 Feature {feature_idx+1}/8: {feature_name} ({title}) 생성 중...")
    
    # 새 figure 생성
    fig, axes = plt.subplots(grid_rows, grid_cols, figsize=(4*grid_cols, 3*grid_rows))
    if best_model is not None:
        suptitle = f'{title} ({feature_name}) - All Experiments\nBest Model (Epoch {best_epoch}, Loss: {best_train_loss:.6f})'
    else:
        suptitle = f'{title} ({feature_name}) - All Experiments\nCurrent Model'
    fig.suptitle(suptitle, fontsize=16, fontweight='bold')
    
    # axes 처리 (1차원으로 변환)
    if grid_rows * grid_cols == 1:
        axes = [axes]
    else:
        axes = axes.flatten() if hasattr(axes, 'flatten') else [axes]
    
    # 성능 지표 수집
    rmse_values = []
    r2_values = []
    
    # 각 실험별 subplot
    for plot_idx, exp_num in enumerate(sorted(all_experiment_results.keys())):
        if plot_idx >= len(axes):
            break
            
        ax = axes[plot_idx]
        exp_data = all_experiment_results[exp_num]
        
        # 데이터 denormalize
        actual_denorm = denormalize_data(exp_data['actual_data'][:, feature_index], feature_name)
        pred_denorm = denormalize_data(exp_data['predicted_data'][:exp_data['seq_len']+1, feature_index], feature_name)
        
        # 플롯 생성
        ax.plot(exp_data['time_points'], actual_denorm, 'b-', linewidth=1.5, label='Actual', alpha=0.8)
        ax.plot(exp_data['time_points'][:exp_data['seq_len']+1], pred_denorm, 'r--', linewidth=1.5, label='Predicted', alpha=0.8)
        
        # 시작점 표시
        ax.plot(exp_data['time_points'][0], actual_denorm[0], 'go', markersize=3)
        ax.axvline(x=exp_data['time_points'][1], color='gray', linestyle=':', alpha=0.5)
        
        # 성능 지표 계산
        actual_pred_denorm = denormalize_data(exp_data['actual_sample'][:, feature_index], feature_name)
        pred_pred_denorm = denormalize_data(exp_data['pred_sample'][:, feature_index], feature_name)
        
        rmse = np.sqrt(np.mean((actual_pred_denorm - pred_pred_denorm)**2))
        if np.sum((actual_pred_denorm - np.mean(actual_pred_denorm))**2) > 0:
            r2 = 1 - (np.sum((actual_pred_denorm - pred_pred_denorm)**2) / np.sum((actual_pred_denorm - np.mean(actual_pred_denorm))**2))
        else:
            r2 = 0.0
            
        rmse_values.append(rmse)
        r2_values.append(r2)
        
        # 서브플롯 설정
        condition_text = f'{exp_data["voltage"]:.1f}V, {exp_data["electrolyte"]:.2f}M'
        ax.set_title(f'Exp {exp_num}: {condition_text}\nRMSE: {rmse:.3f} {unit}', fontsize=9)
        ax.set_xlabel('Time (hours)', fontsize=8)
        ax.set_ylabel(f'{feature_name} ({unit})', fontsize=8)
        ax.grid(True, alpha=0.3)
        ax.tick_params(labelsize=7)
        
        # 첫 번째 subplot에만 범례 추가
        if plot_idx == 0:
            ax.legend(fontsize=8, loc='upper right')
    
    # 빈 subplot 숨기기
    for plot_idx in range(len(all_experiment_results), len(axes)):
        axes[plot_idx].set_visible(False)
    
    # 전체 성능 요약
    avg_rmse = np.mean(rmse_values)
    std_rmse = np.std(rmse_values)
    avg_r2 = np.mean(r2_values)
    std_r2 = np.std(r2_values)
    
    # 성능 요약 텍스트 박스
    summary_text = f'Overall Performance:\nRMSE: {avg_rmse:.4f} ± {std_rmse:.4f} {unit}\nR²: {avg_r2:.3f} ± {std_r2:.3f}\nExperiments: {len(rmse_values)}'
    fig.text(0.02, 0.02, summary_text, fontsize=10, 
             bbox=dict(boxstyle='round,pad=0.5', facecolor='lightblue', alpha=0.8))
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.92, bottom=0.12, left=0.05, right=0.98)
    plt.show()
    
    print(f"✅ {feature_name} 완료: 평균 RMSE = {avg_rmse:.4f} {unit}, 평균 R² = {avg_r2:.3f}")

# 전체 요약 출력
print("\n" + "="*80)
print("COMPREHENSIVE EVALUATION SUMMARY - ALL EXPERIMENTS")
print("="*80)

for feature_idx, (feature_name, feature_index, unit) in enumerate(zip(feature_names, feature_indices, feature_units)):
    # 각 feature별 전체 성능 통계
    all_rmse = []
    all_mae = []
    all_r2 = []
    
    for exp_num in sorted(all_experiment_results.keys()):
        exp_data = all_experiment_results[exp_num]
        
        actual_pred_denorm = denormalize_data(exp_data['actual_sample'][:, feature_index], feature_name)
        pred_pred_denorm = denormalize_data(exp_data['pred_sample'][:, feature_index], feature_name)
        
        rmse = np.sqrt(np.mean((actual_pred_denorm - pred_pred_denorm)**2))
        mae = np.mean(np.abs(actual_pred_denorm - pred_pred_denorm))
        
        if np.sum((actual_pred_denorm - np.mean(actual_pred_denorm))**2) > 0:
            r2 = 1 - (np.sum((actual_pred_denorm - pred_pred_denorm)**2) / np.sum((actual_pred_denorm - np.mean(actual_pred_denorm))**2))
        else:
            r2 = 0.0
        
        all_rmse.append(rmse)
        all_mae.append(mae)
        all_r2.append(r2)
    
    print(f"{feature_name:6s}: RMSE={np.mean(all_rmse):.4f}±{np.std(all_rmse):.4f} {unit}, "
          f"MAE={np.mean(all_mae):.4f}±{np.std(all_mae):.4f} {unit}, "
          f"R²={np.mean(all_r2):.3f}±{np.std(all_r2):.3f}")

print("-"*80)
print(f"📊 Total experiments analyzed: {len(all_experiment_results)}")
if best_model is not None:
    print(f"🏆 Best model used: Epoch {best_epoch}, Loss: {best_train_loss:.6f}")
else:
    print("📋 Current model used")
print("="*80)

print("🎉 Comprehensive evaluation 완료! 8개의 feature별 플롯이 모든 실험 데이터를 포함하여 생성되었습니다.")