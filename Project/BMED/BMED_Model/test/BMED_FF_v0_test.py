import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import torch
import pickle
from matplotlib.gridspec import GridSpec
import sys

# 상대 경로 임포트를 위한 설정
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

# BMED_FF_v0.py 모듈 직접 임포트
from BMED_FF_v0 import (
    MembraneSystemModel,
    PhysicsLayer,
    simulate_with_model,
    plot_simulation_results
)

# PyTorch 2.6 호환성을 위한 모델 로딩 함수 재정의
def load_model_with_scalers(save_dir='result', model_name='bmed_model', device='cuda' if torch.cuda.is_available() else 'cpu'):
    """
    저장된 모델과 스케일러를 불러오는 함수 (PyTorch 2.6 호환성 개선)
    
    Parameters:
    -----------
    save_dir : str
        저장된 디렉토리 경로
    model_name : str
        저장된 모델 이름
    device : str
        모델을 불러올 디바이스 ('cuda' 또는 'cpu')
        
    Returns:
    --------
    model : nn.Module
        불러온 모델
    scalers : dict
        불러온 스케일러 딕셔너리
    """
    # 경로가 상대 경로인 경우 절대 경로로 변환
    if not os.path.isabs(save_dir):
        save_dir = os.path.join(parent_dir, save_dir)
    
    # 파일 경로
    model_path = os.path.join(save_dir, f'{model_name}.pt')
    weights_path = os.path.join(save_dir, f'{model_name}_weights.pt')
    scalers_path = os.path.join(save_dir, f'{model_name}_scalers.pkl')
    
    print(f"모델 경로: {model_path}")
    print(f"가중치 경로: {weights_path}")
    print(f"스케일러 경로: {scalers_path}")
    
    # 파일 존재 확인
    if not os.path.exists(model_path) and not os.path.exists(weights_path):
        raise FileNotFoundError(f"모델 파일이 존재하지 않습니다: {model_path} 또는 {weights_path}")
    
    if not os.path.exists(scalers_path):
        raise FileNotFoundError(f"스케일러 파일이 존재하지 않습니다: {scalers_path}")
    
    # 모델 로딩 시도
    model = None
    print("PyTorch 2.6+ 호환 방식으로 모델 로딩 시도 중...")
    
    # 방법 1: 전체 모델 로딩 시도 (weights_only=False 사용)
    try:
        print("방법 1: 전체 모델 로딩 시도 (weights_only=False)")
        model = torch.load(model_path, map_location=device, weights_only=False)
        print("전체 모델 로딩 성공!")
    except Exception as e:
        print(f"전체 모델 로딩 실패: {str(e)}")
        
        # 방법 2: 모델 클래스 생성 후 가중치만 로딩
        try:
            print("방법 2: 모델 가중치만 로딩 시도")
            # 기본 모델 생성
            model = MembraneSystemModel(hidden_units=256, hidden_layers=20, time_step=0.1)
            model.to(device)
            
            # 가중치 로딩 시도
            if os.path.exists(weights_path):
                model.load_state_dict(torch.load(weights_path, map_location=device))
                print("모델 가중치 로딩 성공!")
            else:
                # 방법 3: 안전한 글로벌 목록에 MembraneSystemModel 추가 후 시도
                print("방법 3: safe_globals 사용 시도")
                from torch.serialization import safe_globals
                with safe_globals([MembraneSystemModel, PhysicsLayer]):
                    model = torch.load(model_path, map_location=device)
                print("safe_globals 방식으로 모델 로딩 성공!")
        except Exception as inner_e:
            print(f"모델 로딩 모든 방법 실패: {str(inner_e)}")
            raise RuntimeError("모델 로딩 실패")
    
    # 스케일러 불러오기
    with open(scalers_path, 'rb') as f:
        scalers = pickle.load(f)
    
    return model, scalers

def compare_experiment_and_simulation():
    """
    BMED 모델의 시뮬레이션 결과와 실험 데이터를 비교하는 함수
    """
    print("BMED 피드포워드 모델 테스트 - 실험 1번과 시뮬레이션 비교")
    
    # 경로 설정
    model_dir = 'result'
    data_path = 'data/BMED_data_v5.csv'
    output_dir = 'result_test'
    
    # 출력 디렉토리 생성
    os.makedirs(output_dir, exist_ok=True)
    
    # 모델 및 스케일러 불러오기
    print("모델 로딩 중...")
    try:
        model, scalers = load_model_with_scalers(save_dir=model_dir, model_name='bmed_model')
        print("모델 로딩 완료")
    except FileNotFoundError as e:
        print(f"파일을 찾을 수 없습니다: {str(e)}")
        return
    except Exception as e:
        print(f"모델 로딩 중 오류 발생: {str(e)}")
        return
    
    # 실험 데이터 로드
    print("실험 데이터 로딩 중...")
    try:
        df = pd.read_csv(data_path)
        exp1_data = df[df['exp'] == 1].sort_values('t')
        if len(exp1_data) == 0:
            raise ValueError("1번 실험 데이터가 존재하지 않습니다.")
        print(f"1번 실험 데이터 로딩 완료: {len(exp1_data)}개 데이터 포인트")
    except Exception as e:
        print(f"데이터 로드 중 오류 발생: {str(e)}")
        return
    
    # 실험 데이터의 초기 상태 가져오기
    print("실험 초기 상태 추출 중...")
    initial_state = exp1_data.iloc[0][['T', 'V', 'E', 'CF_LA', 'CF_K', 'CA_LA', 'CB_K', 'VF', 'VA', 'VB']].values
    print("초기 상태:", initial_state)
    
    # 시뮬레이션 파라미터 설정
    time_step = 0.1  # 시간 간격 (시간)
    max_time = exp1_data['t'].max()  # 실험의 최대 시간까지 시뮬레이션
    steps = int(max_time / time_step) + 1  # 시뮬레이션 단계 수 계산
    
    print(f"시뮬레이션 실행: time_step={time_step}, 최대 시간={max_time}, 시뮬레이션 단계={steps}")
    
    # 시뮬레이션 실행
    states_history, mol_changes_history, t_history = simulate_with_model(
        model=model,
        initial_state=initial_state,
        steps=steps,
        scalers=scalers,
        time_step=time_step
    )
    
    print(f"시뮬레이션 완료: {len(states_history)}개 시뮬레이션 데이터 포인트 생성")
    
    # 시뮬레이션 결과와 실험 데이터 비교 시각화
    print("결과 시각화 중...")
    
    # 1. 변수별 시간에 따른 변화 비교
    plot_simulation_results(
        states_history=states_history,
        t_history=t_history,
        experiment_data=exp1_data,
        save_path=os.path.join(output_dir, 'simulation_vs_experiment_exp1.png')
    )
    
    # 2. 몇 가지 중요 변수에 대한 좀 더 상세한 비교 그래프
    important_vars = ['CF_LA', 'CA_LA', 'CB_K', 'VF', 'VA', 'VB']
    var_indices = {
        'T': 0, 'V': 1, 'E': 2, 'CF_LA': 3, 'CF_K': 4,
        'CA_LA': 5, 'CB_K': 6, 'VF': 7, 'VA': 8, 'VB': 9
    }
    
    fig, axes = plt.subplots(3, 2, figsize=(16, 12))
    axes = axes.flatten()
    
    for i, var in enumerate(important_vars):
        ax = axes[i]
        idx = var_indices[var]
        
        # 실험 데이터
        ax.scatter(exp1_data['t'], exp1_data[var], color='red', marker='o', s=50, 
                   label='실험 데이터', zorder=3)
        
        # 시뮬레이션 결과
        ax.plot(t_history, states_history[:, idx], 'b-', linewidth=2, 
                label='시뮬레이션')
        
        ax.set_xlabel('시간 (시간)')
        ax.set_ylabel(var)
        ax.set_title(f'{var} - 실험 vs 시뮬레이션')
        ax.grid(True, alpha=0.3)
        ax.legend()
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'important_vars_comparison_exp1.png'), dpi=300)
    plt.close()

    # 3. 각 변수별 실험 데이터와 시뮬레이션 데이터의 R² 계산
    print("변수별 R² 점수 계산 중...")
    
    # 실험 데이터 시간과 가장 가까운 시뮬레이션 데이터 포인트 찾기
    var_names = ['T', 'V', 'E', 'CF_LA', 'CF_K', 'CA_LA', 'CB_K', 'VF', 'VA', 'VB']
    r2_scores = {}
    
    # 각 실험 데이터 시간에 가장 가까운 시뮬레이션 시간 인덱스 찾기
    exp_times = exp1_data['t'].values
    sim_indices = []
    
    for t in exp_times:
        idx = np.argmin(np.abs(t_history - t))
        sim_indices.append(idx)
    
    # 선택된 시뮬레이션 데이터 포인트 추출
    selected_sim_data = states_history[sim_indices]
    
    # 각 변수별 R² 계산
    from sklearn.metrics import r2_score
    
    plt.figure(figsize=(15, 12))
    for i, var in enumerate(var_names):
        # 실험 데이터
        exp_values = exp1_data[var].values
        # 대응하는 시뮬레이션 데이터
        sim_values = selected_sim_data[:, i]
        
        # R² 계산
        r2 = r2_score(exp_values, sim_values)
        r2_scores[var] = r2
        
        # 산점도 그리기 (Parity Plot)
        plt.subplot(3, 4, i+1)
        plt.scatter(exp_values, sim_values, alpha=0.7)
        
        # 이상적인 예측선 (y=x)
        min_val = min(np.min(exp_values), np.min(sim_values))
        max_val = max(np.max(exp_values), np.max(sim_values))
        plt.plot([min_val, max_val], [min_val, max_val], 'r--')
        
        plt.title(f'{var}: R² = {r2:.4f}')
        plt.xlabel('실험값')
        plt.ylabel('시뮬레이션값')
        plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'r2_parity_plots_exp1.png'), dpi=300)
    plt.close()
    
    # R² 점수 출력 및 저장
    print("\n각 변수의 R² 점수:")
    for var, r2 in r2_scores.items():
        print(f"{var}: {r2:.4f}")
    
    # R² 점수를 파일로 저장
    with open(os.path.join(output_dir, 'r2_scores_exp1.txt'), 'w') as f:
        f.write("변수,R²\n")
        for var, r2 in r2_scores.items():
            f.write(f"{var},{r2:.6f}\n")
    
    # 4. 몰 변화량 시각화
    mol_change_names = ['LA 몰 변화량', 'K+ 몰 변화량', 'Acid 방향 물 흐름', 'Base 방향 물 흐름']
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    axes = axes.flatten()
    
    for i, name in enumerate(mol_change_names):
        ax = axes[i]
        ax.plot(t_history[:-1], mol_changes_history[:, i], 'g-', linewidth=2)
        ax.set_xlabel('시간 (시간)')
        ax.set_ylabel(name)
        ax.set_title(f'{name} - 시뮬레이션')
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'mol_changes_exp1.png'), dpi=300)
    plt.close()
    
    print(f"\n모든 결과가 '{output_dir}' 디렉토리에 저장되었습니다.")

if __name__ == "__main__":
    compare_experiment_and_simulation() 