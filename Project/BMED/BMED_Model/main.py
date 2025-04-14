import torch
import os
import time
import torch
import json
import pickle

from src.util.hpOpt import hpOpt
from src.util.hpRes import hpRes
from src.model.model_for_train import model_for_train
from src.model.model_for_eval import model_for_eval

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

    # best params 저장 디렉토리 로드
    params_dir = './src/params'

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
        save_data, save_path = hpOptimizer.save_json(results)

    elif mode == 2:
        json_path = os.path.join(params_dir, 'best_trial.json')
        # best_trial.json에서 하이퍼파라미터 로드
        with open(json_path, 'r') as f:
            best_params = json.load(f)['best_params']
        
        # 하이퍼파라미터 추출
        hidden_nodes = best_params['hidden_nodes']
        hidden_layers = best_params['hidden_layers']
        lr = best_params['lr']
        rstop = best_params['rstop'] 
        weight_decay = best_params['weight_decay']
        epochs = best_params['epochs']

        r2_final,r2_min, worst_var, best_model_state, scalers = model_for_train(mode, datapath, fold_idx=0, hidden_nodes=hidden_nodes, hidden_layers=hidden_layers, lr=lr, rstop=rstop, weight_decay=weight_decay, epochs=epochs)
        
        # 모델 저장
        model_save_path = os.path.join(result_dir, 'best_model.pth')
        torch.save(best_model_state, model_save_path)

        # 스케일러 저장
        scaler_save_path = os.path.join(result_dir, 'scalers.pkl')
        with open(scaler_save_path, 'wb') as f:
            pickle.dump(scalers, f)
        
        # 학습 결과 출력
        print(f'R2: {r2_min:.6f}, Worst Var: {worst_var}')
        print("\n final R2 score")
        for var, score in r2_final.items():
            print(f'{var}: {score:.6f}')
    
    elif mode == 3:
        model_for_eval(mode, datapath, params_dir)     


    elif mode == 4:
        res = hpRes(res_dir=result_dir)
        res.execute()

if __name__ == '__main__':
    '''mode
    1: CV
    2: train
    3: simulation
    4: CV selection
    '''
    mode = 2
    datapath = './data/BMED_data_v6+spline.xlsx'
    main(mode, datapath)