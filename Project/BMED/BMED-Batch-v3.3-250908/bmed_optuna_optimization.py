# BMED TF 모델의 Optuna 하이퍼파라미터 최적화
# K-fold Cross Validation과 함께 구현

import torch
import torch.nn as nn
import torch.nn.functional as F
import pandas as pd
from torch.nn.utils.rnn import pad_sequence
from torch.utils.data import TensorDataset, DataLoader, Subset
from sklearn.model_selection import KFold
import numpy as np
import optuna
from optuna.trial import TrialState
import copy
import pickle
import json
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# 기존 모델 클래스들을 임포트하거나 여기에 정의
# (기존 notebook의 모든 클래스들을 그대로 사용)

class LayerNormLSTM(nn.Module):
    """LSTM layer with layer normalization applied to gates"""
    def __init__(self, input_node, hidden_node):
        super().__init__()
        self.input_node = input_node
        self.hidden_node = hidden_node

        self.w_i = nn.Linear(input_node, 4 * hidden_node, bias=False)
        self.w_h = nn.Linear(hidden_node, 4 * hidden_node, bias=False)

        self.ln_i = nn.LayerNorm(hidden_node)
        self.ln_h = nn.LayerNorm(hidden_node)
        self.ln_g = nn.LayerNorm(hidden_node)
        self.ln_o = nn.LayerNorm(hidden_node)
        self.ln_c = nn.LayerNorm(hidden_node)

    def forward(self, input, hidden):
        h_prev, c_prev = hidden
        gi = self.w_i(input)
        gh = self.w_h(h_prev)
        i_i, i_f, i_g, i_o = gi.chunk(4, dim=-1)
        h_i, h_f, h_g, h_o = gh.chunk(4, dim=-1)

        i_g = torch.sigmoid(self.ln_i(i_i + h_i))
        f_g = torch.sigmoid(self.ln_h(i_f + h_f))
        g_g = torch.tanh(self.ln_g(i_g + h_g))
        o_g = torch.sigmoid(self.ln_o(i_o + h_o))

        c_new = f_g * c_prev + i_g * g_g
        c_new = self.ln_c(c_new)
        h_new = o_g * torch.tanh(c_new)

        return h_new, c_new

class StateExtr(nn.Module):
    def __init__(self, input_node, hidden_node, n_layer, dropout):
        super().__init__()
        self.hidden_node = hidden_node
        self.n_layer = n_layer
        self.input_node = input_node

        self.lstm_cells = nn.ModuleList()
        self.lstm_cells.append(LayerNormLSTM(input_node, hidden_node))
        for _ in range(n_layer - 1):
            self.lstm_cells.append(LayerNormLSTM(hidden_node, hidden_node))

        self.dropout = nn.Dropout(dropout)
        self.final_layer_norm = nn.LayerNorm(hidden_node)
        self.final_dropout = nn.Dropout(dropout)

    def forward(self, x, seq_len):
        batch_size, max_len, input_node = x.size()
        device = x.device

        h_states = []
        c_states = []
        for _ in range(self.n_layer):
            h_states.append(torch.zeros(batch_size, self.hidden_node, device=device))
            c_states.append(torch.zeros(batch_size, self.hidden_node, device=device))
        
        outputs = []
        for t in range(max_len):
            x_t = x[:, t, :]
            layer_input = x_t
            for layer_idx, lstm_cell in enumerate(self.lstm_cells):
                h_new, c_new = lstm_cell(layer_input, (h_states[layer_idx], c_states[layer_idx]))
                h_states[layer_idx] = h_new
                c_states[layer_idx] = c_new

                if layer_idx < len(self.lstm_cells) - 1:
                    layer_input = self.dropout(h_new)
                else:
                    layer_input = h_new
            outputs.append(layer_input)
        
        output_tensor = torch.stack(outputs, dim=1)
        seq_len_cpu = seq_len.detach().cpu().long()
        mask = torch.arange(max_len, device='cpu')[None, :] < seq_len_cpu[:, None]
        mask = mask.float().to(device).unsqueeze(-1)
        masked_output = output_tensor * mask
        normalized = self.final_layer_norm(masked_output)
        return self.final_dropout(normalized)

class PhysicalChangeDecoder(nn.Module):
    def __init__(self, input_node, output_node, n_layer, hidden_node, dropout):
        super().__init__()
        self.layers = nn.ModuleList()
        
        self.layers.append(nn.Linear(input_node, hidden_node))
        self.layers.append(nn.LayerNorm(hidden_node))
        self.layers.append(nn.ReLU())
        self.layers.append(nn.Dropout(dropout))

        for i in range(n_layer - 1):
            self.layers.append(nn.Linear(hidden_node, hidden_node))
            self.layers.append(nn.LayerNorm(hidden_node))
            self.layers.append(nn.ReLU())
            self.layers.append(nn.Dropout(dropout))

        self.layers.append(nn.Linear(hidden_node, output_node))
    
    def forward(self, hidden_states):
        x = hidden_states
        for layer in self.layers:
            x = layer(x)
        return x

class CurrentPredictor(nn.Module):
    def __init__(self, input_node, hidden_node, n_layer, dropout):
        super().__init__()
        self.layers = nn.ModuleList()
        
        self.layers.append(nn.Linear(input_node, hidden_node))
        self.layers.append(nn.LayerNorm(hidden_node))
        self.layers.append(nn.ReLU())
        self.layers.append(nn.Dropout(dropout))
        
        for i in range(n_layer - 1):
            self.layers.append(nn.Linear(hidden_node, hidden_node))
            self.layers.append(nn.LayerNorm(hidden_node))
            self.layers.append(nn.ReLU())
            self.layers.append(nn.Dropout(dropout))
        
        self.layers.append(nn.Linear(hidden_node, 1))
    
    def forward(self, new_state):
        x = new_state
        for layer in self.layers:
            x = layer(x)
        return x

class PhysicsConstraintLayer(nn.Module):
    def __init__(self, range_mm, current_predictor, eps=1e-2):
        super().__init__()
        self.sps = eps
        self.current_predictor = current_predictor
        self.register_buffer('range_mm_tensor', self._convert_range_to_tensor(range_mm))

    def _convert_range_to_tensor(self, range_mm):
        feature_names = ['V','E','VF','VA','VB','CFLA','CALA','CFK','CBK','I']
        ranges = torch.zeros(len(feature_names),2)
        for i, name in enumerate(feature_names):
            if name in range_mm:
                ranges[i, 0] = range_mm[name]['min']
                ranges[i, 1] = range_mm[name]['max']
        return ranges
    
    def normalize(self, data, feature_idx):
        min_val = self.range_mm_tensor[feature_idx, 0]
        max_val = self.range_mm_tensor[feature_idx, 1]
        return (data - min_val) / (max_val - min_val)

    def denormalize(self, data, feature_idx):
        min_val = self.range_mm_tensor[feature_idx, 0]
        max_val = self.range_mm_tensor[feature_idx, 1]
        return data * (max_val - min_val) + min_val

    def forward(self, physical_changes, current_state):
        V_idx, E_idx, VF_idx, VA_idx, VB_idx = 0, 1, 2, 3, 4
        CFLA_idx, CALA_idx, CFK_idx, CBK_idx, I_idx = 5, 6, 7, 8, 9

        VF = self.denormalize(current_state[..., 2:3], VF_idx)
        VA = self.denormalize(current_state[..., 3:4], VA_idx)
        VB = self.denormalize(current_state[..., 4:5], VB_idx)
        CFLA = self.denormalize(current_state[..., 5:6], CFLA_idx)
        CALA = self.denormalize(current_state[..., 6:7], CALA_idx)
        CFK = self.denormalize(current_state[..., 7:8], CFK_idx)
        CBK = self.denormalize(current_state[..., 8:9], CBK_idx)

        dVA = physical_changes[..., 0:1]
        dVB = physical_changes[..., 1:2]
        rratio = physical_changes[..., 2:3]
        dNBK = physical_changes[..., 3:4]

        ratio = torch.sigmoid(rratio)
        dNALA = ratio * dNBK

        NFLA = CFLA * VF
        NALA = CALA * VA
        NFK = CFK * VF
        NBK = CBK * VB

        if VF < dVA + dVB:
            dVA = 0
            dVB = 0
        if NFLA < dNALA:
            dNALA = 0
        if NFK < dNBK:
            dNBK = 0

        nVF = VF - dVA - dVB
        nVA = VA + dVA
        nVB = VB + dVB

        nVF = torch.clamp(nVF, min=self.sps)
        nVA = torch.clamp(nVA, min=self.sps)
        nVB = torch.clamp(nVB, min=self.sps)
        
        nNFLA = NFLA - dNALA
        nNALA = NALA + dNALA
        nNFK = NFK - dNBK
        nNBK = NBK + dNBK

        nCFLA = nNFLA / nVF
        nCALA = nNALA / nVA
        nCFK = nNFK / nVF
        nCBK = nNBK / nVB

        V = current_state[..., 0:1]
        E = current_state[..., 1:2]
        nVF_norm = self.normalize(nVF, VF_idx)
        nVA_norm = self.normalize(nVA, VA_idx)
        nVB_norm = self.normalize(nVB, VB_idx)
        nCFLA_norm = self.normalize(nCFLA, CFLA_idx)
        nCALA_norm = self.normalize(nCALA, CALA_idx)
        nCFK_norm = self.normalize(nCFK, CFK_idx)
        nCBK_norm = self.normalize(nCBK, CBK_idx)

        temp_state = torch.cat([
            V, E, nVF_norm, nVA_norm, nVB_norm, nCFLA_norm, nCALA_norm, nCFK_norm, nCBK_norm
        ], dim=-1)
        
        nI_pred_norm = self.current_predictor(temp_state)
        nI_real = self.denormalize(nI_pred_norm, I_idx)
        nI_real = torch.clamp(nI_real, min=0.0)
        nI_norm = self.normalize(nI_real, I_idx)

        next_state = torch.cat([
            V, E, nVF_norm, nVA_norm, nVB_norm, nCFLA_norm, nCALA_norm, nCFK_norm, nCBK_norm, nI_norm
        ], dim=-1)
        
        return next_state

class BMEDAutoregressiveModel(nn.Module):
    def __init__(self, state_extr_params, decoder_params, current_predictor_params, range_mm):
        super().__init__()
        self.state_extr = StateExtr(**state_extr_params)
        self.physical_decoder = PhysicalChangeDecoder(**decoder_params)
        self.current_predictor = CurrentPredictor(**current_predictor_params)
        self.physics_constraint = PhysicsConstraintLayer(range_mm, self.current_predictor)

    def forward(self, x, seq_len):
        hidden_states = self.state_extr(x, seq_len)
        physical_changes = self.physical_decoder(hidden_states)
        new_x = self.physics_constraint(physical_changes, x)
        return new_x

class NoamScheduler:
    def __init__(self, optimizer, model_size, warmup_epochs, factor=1.0):
        self.optimizer = optimizer
        self.model_size = model_size
        self.warmup_epochs = warmup_epochs
        self.factor = factor
        self.epoch_num = 0

    def step_epoch(self):
        self.epoch_num += 1
        lr = self.factor * (
            self.model_size ** (-0.5) *
            min(self.epoch_num ** (-0.5), self.epoch_num * self.warmup_epochs ** (-1.5))
        )
        for param_group in self.optimizer.param_groups:
            param_group['lr'] = lr
        return lr

# 유틸리티 함수들
def df_treat(name):
    df = pd.read_csv(name)
    ndf = pd.DataFrame()
    range_mm={
        'V': {'min':df['V'].min()*0.8, 'max': df['V'].max()*1.2},
        'E': {'min':df['E'].min()*0.8, 'max': df['E'].max()*1.2},
        'VF': {'min':df['VF'].min()*0.8, 'max': df['VF'].max()*1.2},
        'VA': {'min':df['VA'].min()*0.8, 'max': df['VA'].max()*1.2},
        'VB': {'min':df['VB'].min()*0.8, 'max': df['VB'].max()*1.2},
        'CFLA': {'min':0, 'max': df['CFLA'].max()*1.2},
        'CALA': {'min':0, 'max': df['CALA'].max()*1.2},
        'CFK': {'min':0, 'max': df['CFK'].max()*1.2},
        'CBK': {'min':0, 'max': df['CBK'].max()*1.2},
        'I': {'min':0, 'max': df['I'].max()*1.2},
    }
    ndf['exp'] = df['exp']; ndf['t'] = df['t']

    for col in ['V', 'E', 'VF', 'VA', 'VB', 'CFLA', 'CALA', 'CFK', 'CBK', 'I']:
        if col in range_mm:
            ndf[col] = (df[col] - range_mm[col]['min'])/(range_mm[col]['max'] - range_mm[col]['min'])
        else:
            ndf[col] = df[col]

    exp_num_list = sorted(ndf['exp'].unique())
    return df, ndf, range_mm, exp_num_list

def seq_data(ndf, exp_num_list):
    seq = []
    feature_cols = ['V', 'E', 'VF', 'VA', 'VB', 'CFLA', 'CALA', 'CFK', 'CBK', 'I']
    for exp in exp_num_list:
        exp_df = ndf[ndf['exp'] == exp]
        seq.append(exp_df[feature_cols].values)
    return seq

def pad_seq(seq):
    max_len = max([len(s) for s in seq])
    seq_len = [len(s) for s in seq]
    pad_seq = pad_sequence([torch.tensor(s) for s in seq], batch_first=True, padding_value=-1)
    return pad_seq, seq_len, max_len

def gen_dataset(pad_seq, seq_len):
    input_tensor = pad_seq.float()
    seq_len_tensor = torch.tensor(seq_len)
    dataset = TensorDataset(input_tensor, seq_len_tensor)
    return dataset

def masked_mse_loss(pred, target, seq_len):
    batch_size, max_len, features = pred.shape
    seq_len_cpu = seq_len.detach().cpu().long()
    mask = torch.arange(max_len, device='cpu')[None, :] < seq_len_cpu[:, None]
    mask = mask.float().to(pred.device)
    loss = F.mse_loss(pred, target, reduction='none')
    masked_loss = loss * mask.unsqueeze(-1)
    total_loss = masked_loss.sum()
    total_elements = mask.sum()
    masked_loss = total_loss / total_elements
    return masked_loss

def tf_data(input_seq, seq_len):
    inputs = input_seq[:, :-1, :-1]
    targets = input_seq[:, 1:, :]
    target_seq_len = seq_len - 1
    return inputs, targets, target_seq_len

# Optuna 목적 함수
def objective(trial):
    """
    Optuna trial을 위한 목적 함수
    K-fold cross validation을 사용하여 하이퍼파라미터 최적화
    """
    
    # 1. 하이퍼파라미터 제안
    # LSTM StateExtractor 파라미터
    lstm_hidden_size = trial.suggest_categorical('lstm_hidden_size', [16, 32, 48, 64,72,96])
    lstm_n_layers = trial.suggest_int('lstm_n_layers', 2, 4)
    lstm_dropout = trial.suggest_float('lstm_dropout', 0.1, 0.5, step=0.1)
    
    # PhysicalChangeDecoder 파라미터
    decoder_hidden_size = trial.suggest_categorical('decoder_hidden_size', [16, 32, 48, 64,72,96])
    decoder_n_layers = trial.suggest_int('decoder_n_layers', 2, 4)
    decoder_dropout = trial.suggest_float('decoder_dropout', 0.1, 0.6, step=0.1)
    
    # CurrentPredictor 파라미터
    current_hidden_size = trial.suggest_categorical('current_hidden_size', [16, 32, 48, 64,72,96])
    current_n_layers = trial.suggest_int('current_n_layers', 2, 4)
    current_dropout = trial.suggest_float('current_dropout', 0.1, 0.6, step=0.1)
    
    # NoamScheduler 파라미터
    noam_factor = trial.suggest_float('noam_factor', 0.5, 2.0, step=0.1)
    warmup_ratio = trial.suggest_float('warmup_ratio', 0.05, 0.3, step=0.05)
    
    # Batch size 파라미터
    batch_size = trial.suggest_categorical('batch_size', [4, 8, 16])
    
    # 2. K-fold Cross Validation
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    n_splits = 5
    kfold = KFold(n_splits=n_splits, shuffle=True, random_state=42)
    
    fold_losses = []
    
    # 데이터 로드 (global 변수 사용)
    indices = list(range(len(dataset)))
    
    for fold, (train_idx, val_idx) in enumerate(kfold.split(indices)):
        print(f"  🔄 Trial {trial.number}, Fold {fold+1}/{n_splits}")
        
        # 폴드별 데이터셋 준비
        train_subset = Subset(dataset, train_idx)
        val_subset = Subset(dataset, val_idx)
        
        train_loader = DataLoader(train_subset, batch_size=batch_size, shuffle=True)
        val_loader = DataLoader(val_subset, batch_size=batch_size, shuffle=False)
        
        # 3. 모델 파라미터 설정
        state_extr_params = {
            'input_node': 9,
            'hidden_node': lstm_hidden_size,
            'n_layer': lstm_n_layers,
            'dropout': lstm_dropout
        }
        
        decoder_params = {
            'input_node': lstm_hidden_size,
            'hidden_node': decoder_hidden_size,
            'n_layer': decoder_n_layers,
            'dropout': decoder_dropout,
            'output_node': 4
        }
        
        current_predictor_params = {
            'input_node': 9,
            'hidden_node': current_hidden_size,
            'n_layer': current_n_layers,
            'dropout': current_dropout
        }
        
        # 4. 모델 초기화
        model = BMEDAutoregressiveModel(state_extr_params, decoder_params, current_predictor_params, range_mm)
        model = model.to(device)
        
        # 5. 옵티마이저 및 스케줄러 설정
        optimizer = torch.optim.AdamW(model.parameters(), lr=1.0)
        
        # 총 에포크 수와 warmup 에포크 계산
        total_epochs = 100  # Optuna 최적화를 위해 에포크 수 감소
        warmup_epochs = int(total_epochs * warmup_ratio)
        
        scheduler = NoamScheduler(
            optimizer, 
            model_size=lstm_hidden_size,
            warmup_epochs=warmup_epochs,
            factor=noam_factor
        )
        
        # 6. 훈련
        best_total_loss = float('inf')
        patience = 20
        patience_counter = 0
        
        for epoch in range(total_epochs):
            # Learning rate 업데이트
            current_lr = scheduler.step_epoch()
            
            # 훈련
            model.train()
            train_loss = 0.0
            train_batches = 0
            
            for input_seq, seq_len in train_loader:
                try:
                    input_seq = input_seq.to(device)
                    seq_len = seq_len.to(device)
                    
                    inputs, targets, target_seq_len = tf_data(input_seq, seq_len)
                    
                    optimizer.zero_grad()
                    pred = model(inputs, target_seq_len)
                    loss = masked_mse_loss(pred, targets, target_seq_len)
                    
                    loss.backward()
                    torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                    optimizer.step()
                    
                    train_loss += loss.item()
                    train_batches += 1
                    
                except Exception as e:
                    continue
            
            if train_batches == 0:
                break
                
            train_loss = train_loss / train_batches
            
            # 검증
            model.eval()
            val_loss = 0.0
            val_batches = 0
            
            with torch.no_grad():
                for input_seq, seq_len in val_loader:
                    try:
                        input_seq = input_seq.to(device)
                        seq_len = seq_len.to(device)
                        
                        inputs, targets, target_seq_len = tf_data(input_seq, seq_len)
                        
                        pred = model(inputs, target_seq_len)
                        loss = masked_mse_loss(pred, targets, target_seq_len)
                        
                        val_loss += loss.item()
                        val_batches += 1
                        
                    except Exception as e:
                        continue
            
            if val_batches == 0:
                break
                
            val_loss = val_loss / val_batches
            
            # Calculate total loss
            total_loss = train_loss + val_loss
            
            # Early stopping
            if total_loss < best_total_loss:
                best_total_loss = total_loss
                patience_counter = 0
            else:
                patience_counter += 1
                
            if patience_counter >= patience:
                print(f"    Early stopping at epoch {epoch+1}")
                break
        
        fold_losses.append(best_total_loss)
        print(f"    Fold {fold+1} best total loss: {best_total_loss:.6f}")
        
        # 메모리 정리
        del model, optimizer, scheduler
        torch.cuda.empty_cache()
    
    # 7. K-fold 평균 손실 반환
    avg_loss = np.mean(fold_losses)
    std_loss = np.std(fold_losses)
    
    print(f"  📊 Trial {trial.number} - Average CV Loss: {avg_loss:.6f} (±{std_loss:.6f})")
    
    return avg_loss

# 메인 최적화 함수
def run_optuna_optimization():
    """Optuna를 사용한 하이퍼파라미터 최적화 실행"""
    
    print("🚀 BMED TF Model Hyperparameter Optimization with Optuna")
    print("="*80)
    
    # 전역 데이터 로드
    global dataset, range_mm
    
    print("📋 데이터 로드 중...")
    df, ndf, range_mm, exp_num_list = df_treat('BMED_DATA_AG.csv')
    seq = seq_data(ndf, exp_num_list)
    pad, seq_len, max_len = pad_seq(seq)
    dataset = gen_dataset(pad, seq_len)
    
    print(f"   - 총 실험 개수: {len(exp_num_list)}")
    print(f"   - 총 데이터 포인트: {len(dataset)}")
    print(f"   - 최대 시퀀스 길이: {max_len}")
    
    # SQLite 데이터베이스를 사용한 Optuna study 생성
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    # timestamp = "20250903_075457"
    db_url = f"sqlite:///bmed_optuna_study_{timestamp}.db"
    
    study = optuna.create_study(
        direction='minimize',
        study_name='bmed_tf_optimization',
        sampler=optuna.samplers.TPESampler(seed=42),
        storage=db_url,
        load_if_exists=True
    )
    
    # 최적화 실행
    n_trials = 100
    print(f"🔍 최적화 시작 (총 {n_trials} trials)")
    
    try:
        study.optimize(objective, n_trials=n_trials, timeout=None)
    except KeyboardInterrupt:
        print("\n⚠️ 최적화가 사용자에 의해 중단되었습니다.")
    
    # 결과 분석
    print("\n" + "="*80)
    print("📊 OPTIMIZATION RESULTS")
    print("="*80)
    
    print(f"✅ 완료된 trials: {len(study.trials)}")
    print(f"🏆 최고 성능 trial: {study.best_trial.number}")
    print(f"💯 최고 성능 값: {study.best_value:.6f}")
    
    print(f"\n🎯 최적 하이퍼파라미터:")
    for key, value in study.best_params.items():
        print(f"   {key}: {value}")
    
    # 상위 5개 trial 정보
    print(f"\n📈 상위 5개 Trials:")
    trials_df = study.trials_dataframe().sort_values('value').head(5)
    for idx, (_, trial) in enumerate(trials_df.iterrows()):
        print(f"   {idx+1}. Trial {int(trial['number'])}: {trial['value']:.6f}")
    
    # 결과 저장
    result_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    
    # Best parameters 저장
    best_params_file = f"bmed_optuna_best_params_{result_timestamp}.json"
    with open(best_params_file, 'w') as f:
        json.dump(study.best_params, f, indent=2)
    print(f"\n💾 최적 파라미터가 저장되었습니다: {best_params_file}")
    
    # Study 객체 저장
    study_file = f"bmed_optuna_study_{result_timestamp}.pkl"
    with open(study_file, 'wb') as f:
        pickle.dump(study, f)
    print(f"💾 Study 객체가 저장되었습니다: {study_file}")
    
    # Trials 결과 CSV로 저장
    trials_file = f"bmed_optuna_trials_{result_timestamp}.csv"
    trials_df = study.trials_dataframe()
    trials_df.to_csv(trials_file, index=False)
    print(f"💾 모든 trials 결과가 저장되었습니다: {trials_file}")
    
    # SQLite 데이터베이스 정보
    print(f"💾 SQLite 데이터베이스에 실시간 저장됨: {db_url}")
    print(f"   - 중단 후 재시작 시 자동으로 기존 결과를 불러옵니다")
    print(f"   - 다른 프로세스에서 진행상황 모니터링 가능합니다")
    
    print("="*80)
    print("🎉 하이퍼파라미터 최적화 완료!")
    
    return study

if __name__ == "__main__":
    study = run_optuna_optimization()