
from torch.utils.data import Dataset
from sklearn.preprocessing import StandardScaler
import numpy as np
import pandas as pd
import torch

class BMEDDataset(Dataset):
    def __init__(self, mode, data_path, scalers=None):
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
        
