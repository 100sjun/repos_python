import optuna
from optuna.samplers import GPSampler
import numpy as np
import matplotlib.pyplot as plt

# 최적화하고자 하는 1차원 함수 (실제로는 모르는 함수라고 가정)
def objective_function(x):
    return (x - 2) ** 2 * np.sin(x) + 0.1 * x

# Optuna objective 함수
def objective(trial):
    x = trial.suggest_float('x', -5, 5)
    return objective_function(x)

# 시각화 함수
def visualize_optimization(study):
    plt.figure(figsize=(12, 6))
    
    # 실제 함수 그리기
    x = np.linspace(-5, 5, 1000)
    y = objective_function(x)
    plt.plot(x, y, 'b-', label='True Function')
    
    # 시도된 점들 그리기
    trials_x = [t.params['x'] for t in study.trials]
    trials_y = [t.value for t in study.trials]
    
    # 시도된 순서대로 컬러맵 설정
    colors = plt.cm.viridis(np.linspace(0, 1, len(trials_x)))
    
    for i, (x_, y_, c) in enumerate(zip(trials_x, trials_y, colors)):
        if i == 0:
            plt.scatter(x_, y_, color=c, s=100, label='Trials')
        else:
            plt.scatter(x_, y_, color=c, s=100)
    
    # 전역 최적값 표시
    best_trial = study.best_trial
    plt.scatter(best_trial.params['x'], best_trial.value, 
               color='red', s=200, marker='*', label='Best Trial')
    
    plt.xlabel('x')
    plt.ylabel('Objective Value')
    plt.title('Bayesian Optimization Process')
    plt.legend()
    plt.grid(True)
    plt.show()
    
    # 최적값 출력
    print(f"Best value: {best_trial.value:.4f} at x = {best_trial.params['x']:.4f}")

# Optuna study 생성 및 최적화 실행
sampler = GPSampler(n_startup_trials=10)  # 처음 5번은 랜덤 샘플링
study = optuna.create_study(sampler=sampler)
study.optimize(objective, n_trials=100)

# 최적화 과정 시각화
visualize_optimization(study)