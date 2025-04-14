from src.data.BMEDDataset import BMEDDataset
from src.model.MembraneSystemModel import MembraneSystemModel
from torch.utils.data import DataLoader
from src.data.BatchSampler import BatchSampler
import os
import pickle
import json
import torch
import numpy as np
from sklearn.metrics import r2_score


def model_for_eval(mode, datapath, params_dir):

    # path 설정
    json_path = os.path.join(params_dir, 'best_trial.json')
    model_save_path = os.path.join(params_dir, 'best_model.pth')
    scaler_save_path = os.path.join(params_dir, 'scalers.pkl')

    # best_trial.json에서 하이퍼파라미터 로드
    with open(json_path, 'r') as f:
        best_params = json.load(f)['best_params']
    
    # 하이퍼파라미터 추출   
    hidden_nodes = best_params['hidden_nodes']
    hidden_layers = best_params['hidden_layers']

    # 스케일러 로드
    scaler_save_path = os.path.join(params_dir, 'scalers.pkl')
    with open(scaler_save_path, 'rb') as f:
        scalers = pickle.load(f) 

    # dataset 생성
    dataset = BMEDDataset(
        mode = mode,
        data_path=datapath,
        fold_idx=0,
        scalers=scalers,
        train=False
    )

    # 모델 생성
    model = MembraneSystemModel(
        hidden_nodes=hidden_nodes,
        hidden_layers=hidden_layers
    )        

    model_state = torch.load(model_save_path)
    model.load_state_dict(model_state)

    # 데이터 로더 생성
    eval_batch_sampler = BatchSampler(dataset)
    eval_loader = DataLoader(dataset, batch_sampler=eval_batch_sampler)

    # 모델 평가
    model.eval()
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model.to(device)

    # 결과 저장용 리스트
    all_preds = []
    all_actuals = []

    # 평가 실행
    with torch.no_grad():
        for features, migrations, states in eval_loader:
            # 데이터를 GPU로 이동
            features = features.to(device)
            states = states.to(device)
            
            # 예측 수행
            migration_pred, state_pred = model(features)
            
            # CPU로 이동하고 numpy로 변환
            state_pred = state_pred.cpu().numpy()
            states = states.cpu().numpy()
            
            # 결과 저장
            all_preds.append(state_pred)
            all_actuals.append(states)
            
        # 결과 합치기
    all_preds = np.concatenate(all_preds, axis=0)
    all_actuals = np.concatenate(all_actuals, axis=0)
    
    # 역정규화 수행
    predictions_denorm = scalers['state'].inverse_transform(all_preds)
    actuals_denorm = scalers['state'].inverse_transform(all_actuals)
    
    # R2 점수 계산
    var_names = ['T', 'V', 'E', 'CF_LA', 'CF_K', 'CA_LA', 'CB_K', 'VF', 'VA', 'VB']
    r2_scores = {}
    
    for i, var in enumerate(var_names):
        r2 = r2_score(actuals_denorm[:, i], predictions_denorm[:, i])
        r2_scores[var] = r2
    
    # 결과 출력
    print("\n평가 결과:")
    print("\n각 변수별 R2 점수:")
    for var, score in r2_scores.items():
        print(f"{var}: {score:.4f}")
    
    min_r2 = min(r2_scores.values())
    worst_var = min(r2_scores.keys(), key=lambda k: r2_scores[k])
    print(f"\n최소 R2 점수: {min_r2:.4f} ({worst_var})")
    
    return {
        'predictions': predictions_denorm,
        'actuals': actuals_denorm,
        'r2_scores': r2_scores,
        'min_r2': min_r2,
        'worst_var': worst_var
    }


