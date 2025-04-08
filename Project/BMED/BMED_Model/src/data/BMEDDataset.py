
from torch.utils.data import Dataset
from sklearn.preprocessing import StandardScaler
import numpy as np
import pandas as pd
import torch

class BMEDDataset(Dataset):
    def __init__(self, mode, data_path, fold_idx=0, train=True, scalers=None):
        self.mode = mode
        self.fold_idx = fold_idx
        self.train = train

        # 데이터 로드
        if mode == 1 or mode == 2:
            self.df = pd.read_excel(data_path, sheet_name='spline_data')
        elif mode == 3:
            self.df = pd.read_excel(data_path, sheet_name='raw_data')

        # 스케일러 초기화 혹은 전달받은 스케일러 사용
        if scalers is None:
            self.feature_scaler = StandardScaler()  
            self.migration_scaler = StandardScaler()
            self.state_scaler = StandardScaler()
        else:
            self.feature_scaler = scalers['feature']
            self.migration_scaler = scalers['migration']
            self.state_scaler = scalers['state']

        # 실험 번호별로 데이터 분리
        self.dict_spline = {}
        for exp in self.df['exp'].unique():
            self.dict_spline[exp] = self.df[self.df['exp'] == exp].sort_values('t')

        # k-fold를 위한 실험 분할
        all_exps = list(self.dict_spline.keys())
        np.random.seed(42)
        np.random.shuffle(all_exps)
        np.random.seed(None)

        if mode == 1:
            fold_size = len(all_exps) // 5 + 1
            test_start = fold_idx * fold_size
            test_end = (fold_idx + 1) * fold_size if fold_idx < 4 else len(all_exps)

            # 실험 훈련 데이터 분리
            if train:
                self.exps_to_use = all_exps[:test_start] + all_exps[test_end:]
            else:
                self.exps_to_use = all_exps[test_start:test_end]

        # 데이터 준비
        self.prepare_data()

    def get_scalers(self):
        """스케일러 반환"""
        return {
            'feature': self.feature_scaler,
            'migration': self.migration_scaler,
            'state': self.state_scaler
        }
    
    def prepare_data(self):
        features_list = []
        migrations_list = []
        states_list = []
        
        # 하나의 실험셋 로드
        for exp in self.exps_to_use:
            exp_data = self.dict_spline[exp]

            # 각 시간 단계에 대해 하나의 샘플 생성, 마지막 데이터 포인트는 다음 상태가 없으므로 패스
            for i in range(len(exp_data) - 1): 
                cur = exp_data.iloc[i]
                next = exp_data.iloc[i+1]

                # current feature vector
                feature = [
                    cur['T'], cur['V'], cur['E'],
                    cur['CF_LA'], cur['CF_K'], cur['CA_LA'], cur['CB_K'],
                    cur['VF'], cur['VA'], cur['VB']
                ]

                # time step
                dt = next['t'] - cur['t']

                # mole number (mol)
                cur_NA_LA = cur['CA_LA'] * cur['VA']
                cur_NB_K = cur['CB_K'] * cur['VB']
                next_NA_LA = next['CA_LA'] * next['VA']
                next_NB_K = next['CB_K'] * next['VB']

                # migration (mol/hr)
                dNLA = (next_NA_LA - cur_NA_LA) / dt
                dNK = (next_NB_K - cur_NB_K) / dt
                dVA = (next['VA'] - cur['VA']) / dt
                dVB = (next['VB'] - cur['VB']) / dt

                # next state vector
                state = [
                    next['T'], next['V'], next['E'],
                    next['CF_LA'], next['CF_K'], next['CA_LA'], next['CB_K'],
                    next['VF'], next['VA'], next['VB']
                ]

                features_list.append(feature)
                migrations_list.append([dNLA, dNK, dVA, dVB])
                states_list.append(state)
                
        # 데이터 배열로 변환
        self.features_raw = np.array(features_list, dtype=np.float32)
        self.migrations_raw = np.array(migrations_list, dtype=np.float32)
        self.states_raw = np.array(states_list, dtype=np.float32)

        # 데이터 정규화
        if self.train:
            self.feature_scaler.fit(self.features_raw)
            self.migration_scaler.fit(self.migrations_raw)
            self.state_scaler.fit(self.states_raw)

        # 정규화 적용
        self.features = self.feature_scaler.transform(self.features_raw)
        self.migrations = self.migration_scaler.transform(self.migrations_raw)
        self.states = self.state_scaler.transform(self.states_raw)

    def __len__(self):
        return len(self.features)
    
    def __getitem__(self, idx):
        return (
            torch.FloatTensor(self.features[idx]).unsqueeze(0),
            torch.FloatTensor(self.migrations[idx]).unsqueeze(0),
            torch.FloatTensor(self.states[idx]).unsqueeze(0)
        )
        

            




        
