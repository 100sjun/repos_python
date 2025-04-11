import torch
import os
import time

from src.model.model_with_kfold import model_with_kfold

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
    result_dir = f'./result/{timestamp}'
    # os.makedirs(result_dir, exist_ok=True)

    # 모드 설정
    if mode == 1:
        print('Mode: Cross Validation')
    elif mode == 2:
        print('Mode: Train')
    elif mode == 3:
        print('Mode: Simulation')

    # Model Selection
    if mode == 1:
        model = model_with_kfold(mode, datapath)
        print(model)

if __name__ == '__main__':
    '''mode
    1: CV
    2: train
    3: simulation
    '''
    mode = 1
    datapath = './data/BMED_data_v6+spline.xlsx'
    main(mode, datapath)