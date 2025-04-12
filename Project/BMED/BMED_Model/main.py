import torch
import os
import time
import json

from src.util.hpOpt import hpOpt

def main(mode, datapath):
    # CUDA 가용성 확인
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    print(f'장치 사용: {device}')

    if device == 'cuda':
        # GPU 이름 출력
        gpu_name = torch.cuda.get_device_name(0)
        print(f'GPU 이름: {gpu_name}')

    # 결과 저장 디렉토리 생성
    timestamp = time.strftime('%Y%m%d%H%M%S')
    result_dir = f'./results/{timestamp}'
    os.makedirs(result_dir, exist_ok=True)

    # 모드 설정
    if mode == 1:
        print('Mode: Cross Validation')
    elif mode == 2:
        print('Mode: Train')
    elif mode == 3:
        print('Mode: Simulation')

    # Model Selection
    if mode == 1:
        hp_ranges = {
            'hidden_nodes': (16, 512),
            'hidden_layers': (3, 50),
            'lr': (1e-5, 1e-1),
            'rstop': (0.05, 0.3),
            'weight_decay': (1e-5, 1e-1),
            'epochs': (50, 3000)
        }
        hpOptimizer = hpOpt(n_trials=1000, mode=mode, datapath=datapath, hp_ranges=hp_ranges, result_dir=result_dir)
        results = hpOptimizer.optimize()

        save_data, json_path = hpOptimizer.save_json(results)
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(save_data, f, indent=4)

if __name__ == '__main__':
    '''mode
    1: CV
    2: train
    3: simulation
    '''
    mode = 1
    datapath = './data/BMED_data_v6+spline.xlsx'
    main(mode, datapath)