# BMED LSTM-MLP 하이퍼파라미터 최적화 워크플로우

## 📊 프로젝트 개요
**목표**: Bipolar Membrane Electrodialysis (BMED) 배치 실험 데이터를 이용한 LSTM-MLP 하이브리드 모델 개발  
**데이터**: 0.25시간 간격 측정, Feed/Acid/Base의 LA/K 농도, 부피, 전압, 전류 변화  
**정규화**: 20% 확장 MinMaxScaler (정보 보존 + 외삽 안정성)  
**최적화**: Optuna 기반 하이퍼파라미터 튜닝, 5-fold Cross Validation

## 데이터 구조 분석
**입력 변수**: exp, V, E, t, VF, VA, VB, CFLA, CALA, CBLA, CFK, CAK, CBK, I  
**시계열 특성**: 0.25시간 간격, 배치별 시퀀스 데이터  
**타겟 변수**: dCALA, dCBLA, dCAK, dCBK, dVA, dVB (변화량)

---

# Phase 1: 데이터 전처리 & 구조화 (2-3일)

## 1.1 데이터 로딩 & 탐색적 분석
```python
# 순차적 코드 작성 (함수화 나중에)
import pandas as pd
import numpy as np
import torch
import torch.nn as nn
from sklearn.model_selection import KFold
import optuna
import matplotlib.pyplot as plt

# 데이터 로딩
df = pd.read_csv('BMED_DATA_AG.csv')
experiments = df['exp'].unique()
time_steps = df['t'].unique()

print(f"실험 수: {len(experiments)}")
print(f"시간 단계: {len(time_steps)}")
print(f"데이터 형태: {df.shape}")
```

## 1.2 확장 MinMax Scaling 구현 (20% 외삽 허용)
```python
# 20% 확장 MinMax Scaling - 정보 보존 + 외삽 안정성
class ExtendedMinMaxScaler:
    def __init__(self, extend_ratio=0.2):
        self.extend_ratio = extend_ratio
        self.original_min_ = None
        self.original_max_ = None
        self.extended_min_ = None
        self.extended_max_ = None
    
    def fit(self, X):
        self.original_min_ = X.min(axis=0)
        self.original_max_ = X.max(axis=0)
        
        # 20% 확장 범위 계산 (외삽 허용)
        range_size = self.original_max_ - self.original_min_
        extension = range_size * self.extend_ratio
        
        self.extended_min_ = self.original_min_ - extension
        self.extended_max_ = self.original_max_ + extension
        
        # 농도 변수는 0 이하로 갈 수 없음 (물리적 제약)
        concentration_cols = [7, 8, 9, 10, 11]  # CFLA, CALA, CBLA, CFK, CAK, CBK 인덱스
        for i in concentration_cols:
            if i < len(self.extended_min_):
                self.extended_min_[i] = max(0, self.extended_min_[i])
        
        return self
    
    def transform(self, X):
        return (X - self.extended_min_) / (self.extended_max_ - self.extended_min_)
    
    def inverse_transform(self, X_scaled):
        return X_scaled * (self.extended_max_ - self.extended_min_) + self.extended_min_

# 변수별 맞춤형 확장 비율 설정
extension_ratios = {
    'V': 0.15,      # 전압: 15% (측정 안정적)
    'E': 0.20,      # 전기장: 20% (변동 있음)
    'VF': 0.10,     # Feed 부피: 10% (물리적 제약)
    'VA': 0.10,     # Acid 부피: 10% (물리적 제약)  
    'VB': 0.10,     # Base 부피: 10% (물리적 제약)
    'CFLA': 0.25,   # Feed LA농도: 25% (큰 변동)
    'CALA': 0.25,   # Acid LA농도: 25% (큰 변동)
    'CBLA': 0.25,   # Base LA농도: 25% (큰 변동)
    'CFK': 0.25,    # Feed K농도: 25% (큰 변동)
    'CAK': 0.25,    # Acid K농도: 25% (큰 변동)
    'CBK': 0.25,    # Base K농도: 25% (큰 변동)
    'I': 0.30       # 전류: 30% (가장 변동 심함)
}

feature_cols = ['V', 'E', 'VF', 'VA', 'VB', 'CFLA', 'CALA', 'CBLA', 'CFK', 'CAK', 'CBK', 'I']

# 변수별 개별 스케일링 적용
scalers = {}
df_scaled = df.copy()

for i, col in enumerate(feature_cols):
    ratio = extension_ratios.get(col, 0.2)  # 기본 20%
    scaler = ExtendedMinMaxScaler(extend_ratio=ratio)
    df_scaled[col] = scaler.fit_transform(df[[col]]).flatten()
    scalers[col] = scaler

print("변수별 스케일링 범위:")
for col in feature_cols:
    orig_min, orig_max = df[col].min(), df[col].max()
    ext_min, ext_max = scalers[col].extended_min_[0], scalers[col].extended_max_[0]
    print(f"{col}: [{orig_min:.3f}, {orig_max:.3f}] → [{ext_min:.3f}, {ext_max:.3f}]")
```

## 1.3 시퀀스 데이터 구조화
```python
# 실험별 시퀀스 분리 (teacher forcing용) - 스케일링된 데이터 사용
sequences = []
for exp in experiments:
    exp_data = df_scaled[df_scaled['exp'] == exp].sort_values('t')
    sequences.append(exp_data[feature_cols].values)

print(f"시퀀스 개수: {len(sequences)}")
print(f"평균 시퀀스 길이: {np.mean([len(seq) for seq in sequences])}")
print(f"스케일링된 데이터 범위 확인:")
for i, col in enumerate(feature_cols):
    seq_values = [seq[:, i] for seq in sequences]
    all_values = np.concatenate(seq_values)
    print(f"{col}: [{all_values.min():.3f}, {all_values.max():.3f}]")
```

## 1.4 변화량 계산 (Delta 타겟)
```python
# dCALA, dCBLA, dCAK, dCBK, dVA, dVB 계산 - 스케일링된 데이터로 계산
delta_targets = ['dCALA', 'dCBLA', 'dCAK', 'dCBK', 'dVA', 'dVB']
target_sequences = []
target_cols = ['CALA', 'CBLA', 'CAK', 'CBK', 'VA', 'VB']

for exp in experiments:
    exp_data = df_scaled[df_scaled['exp'] == exp].sort_values('t')
    
    # 스케일링된 데이터로 변화량 계산
    deltas = np.diff(exp_data[target_cols].values, axis=0)
    target_sequences.append(deltas)

print(f"타겟 시퀀스 형태: {[seq.shape for seq in target_sequences[:3]]}")

# 변화량 범위 확인
all_deltas = np.concatenate(target_sequences, axis=0)
print("변화량 통계:")
for i, col in enumerate(target_cols):
    delta_values = all_deltas[:, i]
    print(f"d{col}: mean={delta_values.mean():.4f}, std={delta_values.std():.4f}, "
          f"range=[{delta_values.min():.4f}, {delta_values.max():.4f}]")
```

---

# Phase 2: 모델 아키텍처 설계 (3-4일)

## 2.1 LSTM 인코더 설계
```python
# LSTM 하이퍼파라미터 범위 정의
lstm_hidden_sizes = [32, 64, 128, 256]
lstm_num_layers = [1, 2, 3]
lstm_dropout = [0.1, 0.2, 0.3]

class LSTMEncoder(nn.Module):
    def __init__(self, input_size, hidden_size, num_layers, dropout=0.1):
        super(LSTMEncoder, self).__init__()
        self.lstm = nn.LSTM(input_size, hidden_size, num_layers, 
                           batch_first=True, dropout=dropout)
        
    def forward(self, x):
        output, (hidden, cell) = self.lstm(x)
        return output, hidden, cell
```

## 2.2 MLP 디코더 설계
```python
# MLP 하이퍼파라미터 범위
mlp_hidden_sizes = [64, 128, 256, 512]
mlp_num_layers = [2, 3, 4]
mlp_dropout = [0.1, 0.2, 0.3]

class MLPDecoder(nn.Module):
    def __init__(self, input_size, hidden_size, output_size, num_layers, dropout=0.1):
        super(MLPDecoder, self).__init__()
        layers = []
        
        # 첫 번째 레이어
        layers.append(nn.Linear(input_size, hidden_size))
        layers.append(nn.ReLU())
        layers.append(nn.Dropout(dropout))
        
        # 중간 레이어들
        for _ in range(num_layers - 2):
            layers.append(nn.Linear(hidden_size, hidden_size))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(dropout))
        
        # 출력 레이어
        layers.append(nn.Linear(hidden_size, output_size))
        
        self.mlp = nn.Sequential(*layers)
    
    def forward(self, x):
        return self.mlp(x)
```

## 2.3 Mass Balance Layer 구현
```python
# 물리 법칙 기반 계산 레이어
class MassBalanceLayer(nn.Module):
    def __init__(self):
        super(MassBalanceLayer, self).__init__()
    
    def forward(self, deltas, current_state):
        """
        deltas: [dCALA, dCBLA, dCAK, dCBK, dVA, dVB]
        current_state: [CALA, CBLA, CAK, CBK, VA, VB]
        """
        # 농도 변화 적용
        next_cala = current_state[:, 0] + deltas[:, 0]
        next_cbla = current_state[:, 1] + deltas[:, 1]
        next_cak = current_state[:, 2] + deltas[:, 2]
        next_cbk = current_state[:, 3] + deltas[:, 3]
        
        # 부피 변화 적용
        next_va = current_state[:, 4] + deltas[:, 4]
        next_vb = current_state[:, 5] + deltas[:, 5]
        
        # 물리적 제약 조건 적용 (농도 >= 0, 부피 >= 0)
        next_cala = torch.clamp(next_cala, min=0)
        next_cbla = torch.clamp(next_cbla, min=0)
        next_cak = torch.clamp(next_cak, min=0)
        next_cbk = torch.clamp(next_cbk, min=0)
        next_va = torch.clamp(next_va, min=0)
        next_vb = torch.clamp(next_vb, min=0)
        
        return torch.stack([next_cala, next_cbla, next_cak, next_cbk, next_va, next_vb], dim=1)
```

## 2.4 통합 모델 설계
```python
class BMEDModel(nn.Module):
    def __init__(self, lstm_params, mlp_params):
        super(BMEDModel, self).__init__()
        self.lstm_encoder = LSTMEncoder(**lstm_params)
        self.mlp_decoder = MLPDecoder(**mlp_params)
        self.mass_balance = MassBalanceLayer()
        
    def forward(self, x, current_state, teacher_forcing_ratio=0.8):
        # LSTM으로 시퀀스 인코딩
        lstm_output, _, _ = self.lstm_encoder(x)
        
        # MLP로 변화량 예측
        deltas = self.mlp_decoder(lstm_output[:, -1, :])  # 마지막 타임스텝 사용
        
        # Mass Balance로 다음 상태 계산
        next_state = self.mass_balance(deltas, current_state)
        
        return deltas, next_state
```

## 2.5 Teacher Forcing 메커니즘
```python
# 훈련 시 실제 값 사용, 추론 시 예측 값 사용
def teacher_forcing_step(model, input_seq, target_seq, current_states, teacher_forcing_ratio=0.8):
    batch_size, seq_len, input_size = input_seq.shape
    outputs = []
    
    for t in range(seq_len - 1):
        # 현재까지의 시퀀스로 다음 상태 예측
        current_input = input_seq[:, :t+1, :]
        current_state = current_states[:, t, :]
        
        deltas, next_state = model(current_input, current_state)
        outputs.append(deltas)
        
        # Teacher forcing 적용
        if np.random.random() < teacher_forcing_ratio:
            # 실제 값 사용
            current_states[:, t+1, :] = target_seq[:, t, :]
        else:
            # 예측 값 사용
            current_states[:, t+1, :] = next_state.detach()
    
    return torch.stack(outputs, dim=1)
```

---

# Phase 3: Optuna 하이퍼파라미터 최적화 (4-5일)

## 3.1 Objective Function 설계
```python
def objective(trial):
    # LSTM 하이퍼파라미터
    lstm_hidden = trial.suggest_categorical('lstm_hidden', [32, 64, 128, 256])
    lstm_layers = trial.suggest_int('lstm_layers', 1, 3)
    lstm_dropout = trial.suggest_float('lstm_dropout', 0.1, 0.3)
    
    # MLP 하이퍼파라미터
    mlp_hidden = trial.suggest_categorical('mlp_hidden', [64, 128, 256, 512])
    mlp_layers = trial.suggest_int('mlp_layers', 2, 4)
    mlp_dropout = trial.suggest_float('mlp_dropout', 0.1, 0.3)
    
    # 학습률 & 배치 사이즈
    lr = trial.suggest_loguniform('lr', 1e-5, 1e-2)
    batch_size = trial.suggest_categorical('batch_size', [16, 32, 64])
    
    # 모델 파라미터 설정
    lstm_params = {
        'input_size': len(feature_cols),
        'hidden_size': lstm_hidden,
        'num_layers': lstm_layers,
        'dropout': lstm_dropout
    }
    
    mlp_params = {
        'input_size': lstm_hidden,
        'hidden_size': mlp_hidden,
        'output_size': 6,  # dCALA, dCBLA, dCAK, dCBK, dVA, dVB
        'num_layers': mlp_layers,
        'dropout': mlp_dropout
    }
    
    # 5-Fold Cross Validation
    kfold = KFold(n_splits=5, shuffle=True, random_state=42)
    cv_scores = []
    
    for fold, (train_idx, val_idx) in enumerate(kfold.split(sequences)):
        # 폴드별 데이터 분할
        train_sequences = [sequences[i] for i in train_idx]
        val_sequences = [sequences[i] for i in val_idx]
        train_targets = [target_sequences[i] for i in train_idx]
        val_targets = [target_sequences[i] for i in val_idx]
        
        # 모델 초기화
        model = BMEDModel(lstm_params, mlp_params)
        optimizer = torch.optim.Adam(model.parameters(), lr=lr)
        criterion = nn.MSELoss()
        
        # 훈련
        model.train()
        for epoch in range(50):  # 빠른 평가를 위해 50 에포크로 제한
            train_loss = 0
            for i, (seq, target) in enumerate(zip(train_sequences, train_targets)):
                if len(seq) < 2:  # 최소 2개 타임스텝 필요
                    continue
                    
                optimizer.zero_grad()
                
                # 배치 데이터 준비
                input_tensor = torch.FloatTensor(seq[:-1]).unsqueeze(0)
                target_tensor = torch.FloatTensor(target).unsqueeze(0)
                current_states = torch.FloatTensor(seq[1:, [7, 8, 10, 11, 5, 6]]).unsqueeze(0)  # CALA, CBLA, CAK, CBK, VA, VB
                
                # 예측
                predicted_deltas = teacher_forcing_step(model, input_tensor, target_tensor, current_states)
                
                # 손실 계산
                loss = criterion(predicted_deltas, target_tensor)
                loss.backward()
                optimizer.step()
                
                train_loss += loss.item()
        
        # 검증
        model.eval()
        val_loss = 0
        val_count = 0
        
        with torch.no_grad():
            for i, (seq, target) in enumerate(zip(val_sequences, val_targets)):
                if len(seq) < 2:
                    continue
                    
                input_tensor = torch.FloatTensor(seq[:-1]).unsqueeze(0)
                target_tensor = torch.FloatTensor(target).unsqueeze(0)
                current_states = torch.FloatTensor(seq[1:, [7, 8, 10, 11, 5, 6]]).unsqueeze(0)
                
                predicted_deltas = teacher_forcing_step(model, input_tensor, target_tensor, current_states, teacher_forcing_ratio=0.0)
                loss = criterion(predicted_deltas, target_tensor)
                
                val_loss += loss.item()
                val_count += 1
        
        if val_count > 0:
            cv_scores.append(val_loss / val_count)
    
    return np.mean(cv_scores) if cv_scores else float('inf')
```

## 3.2 최적화 실행
```python
# Optuna study 생성 및 실행
study = optuna.create_study(direction='minimize')
study.optimize(objective, n_trials=200)

# 최적 파라미터 출력
best_params = study.best_params
best_score = study.best_value

print("최적 하이퍼파라미터:")
for key, value in best_params.items():
    print(f"  {key}: {value}")
print(f"최적 CV 점수: {best_score}")
```

---

# Phase 4: 최종 모델 훈련 & 평가 (3-4일)

## 4.1 최적 모델 훈련
```python
# 최적 파라미터로 모델 구성
lstm_params = {
    'input_size': len(feature_cols),
    'hidden_size': best_params['lstm_hidden'],
    'num_layers': best_params['lstm_layers'],
    'dropout': best_params['lstm_dropout']
}

mlp_params = {
    'input_size': best_params['lstm_hidden'],
    'hidden_size': best_params['mlp_hidden'],
    'output_size': 6,
    'num_layers': best_params['mlp_layers'],
    'dropout': best_params['mlp_dropout']
}

# 최종 모델 훈련
final_model = BMEDModel(lstm_params, mlp_params)
optimizer = torch.optim.Adam(final_model.parameters(), lr=best_params['lr'])
criterion = nn.MSELoss()

# 전체 데이터로 훈련
epochs = 200
train_losses = []

for epoch in range(epochs):
    final_model.train()
    epoch_loss = 0
    batch_count = 0
    
    for i, (seq, target) in enumerate(zip(sequences, target_sequences)):
        if len(seq) < 2:
            continue
            
        optimizer.zero_grad()
        
        input_tensor = torch.FloatTensor(seq[:-1]).unsqueeze(0)
        target_tensor = torch.FloatTensor(target).unsqueeze(0)
        current_states = torch.FloatTensor(seq[1:, [7, 8, 10, 11, 5, 6]]).unsqueeze(0)
        
        predicted_deltas = teacher_forcing_step(final_model, input_tensor, target_tensor, current_states)
        loss = criterion(predicted_deltas, target_tensor)
        
        loss.backward()
        optimizer.step()
        
        epoch_loss += loss.item()
        batch_count += 1
    
    if batch_count > 0:
        avg_loss = epoch_loss / batch_count
        train_losses.append(avg_loss)
        
        if (epoch + 1) % 20 == 0:
            print(f"Epoch {epoch+1}/{epochs}, Loss: {avg_loss:.6f}")
```

## 4.2 결과 시각화
```python
import matplotlib.pyplot as plt

# 훈련 손실 그래프
plt.figure(figsize=(10, 6))
plt.plot(train_losses)
plt.title('Training Loss')
plt.xlabel('Epoch')
plt.ylabel('MSE Loss')
plt.grid(True)
plt.show()

# Optuna 결과 시각화
optuna.visualization.plot_param_importances(study).show()
optuna.visualization.plot_optimization_history(study).show()
```

## 4.3 물리 법칙 검증
```python
# Mass balance 검증
def validate_mass_balance(model, sequences, tolerance=1e-3):
    model.eval()
    violations = []
    
    with torch.no_grad():
        for seq in sequences:
            if len(seq) < 2:
                continue
                
            input_tensor = torch.FloatTensor(seq[:-1]).unsqueeze(0)
            current_states = torch.FloatTensor(seq[1:, [7, 8, 10, 11, 5, 6]]).unsqueeze(0)
            
            deltas, next_states = model(input_tensor, current_states[:, 0, :].unsqueeze(1))
            
            # 물리적 제약 위반 확인
            negative_concentrations = (next_states < 0).sum().item()
            if negative_concentrations > 0:
                violations.append(negative_concentrations)
    
    print(f"물리 법칙 위반 사례: {len(violations)}/{len(sequences)}")
    return violations

violations = validate_mass_balance(final_model, sequences)
```

---

# Phase 5: 모델 저장 & 배포 준비 (1-2일)

## 5.1 모델 저장
```python
# 모델과 관련 정보 저장
torch.save({
    'model_state_dict': final_model.state_dict(),
    'lstm_params': lstm_params,
    'mlp_params': mlp_params,
    'best_params': best_params,
    'scalers': scalers,  # 변수별 스케일러 딕셔너리
    'extension_ratios': extension_ratios,  # 확장 비율 정보
    'cv_score': best_score,
    'train_losses': train_losses,
    'feature_cols': feature_cols,
    'target_cols': target_cols
}, 'bmed_lstm_mlp_optimized.pth')

print("모델이 'bmed_lstm_mlp_optimized.pth'에 저장되었습니다.")
```

## 5.2 추론 파이프라인
```python
def predict_next_state(model, scalers, current_sequence, feature_cols, target_cols):
    """
    새로운 데이터에 대한 예측 함수 (20% 확장 MinMaxScaler 사용)
    """
    model.eval()
    
    # 변수별 개별 전처리
    normalized_data = {}
    for col in feature_cols:
        scaler = scalers[col]
        normalized_data[col] = scaler.transform(current_sequence[[col]]).flatten()
    
    # 정규화된 시퀀스 구성
    normalized_seq = np.column_stack([normalized_data[col] for col in feature_cols])
    input_tensor = torch.FloatTensor(normalized_seq).unsqueeze(0)
    
    # 현재 상태 추출 (CALA, CBLA, CAK, CBK, VA, VB)
    target_indices = [feature_cols.index(col) for col in target_cols]
    current_state = torch.FloatTensor(normalized_seq[-1, target_indices]).unsqueeze(0)
    
    with torch.no_grad():
        deltas, next_state = model(input_tensor, current_state)
    
    # 역정규화 (변수별 개별 처리)
    next_state_denorm = np.zeros_like(next_state.numpy())
    for i, col in enumerate(target_cols):
        scaler = scalers[col]
        next_state_denorm[:, i] = scaler.inverse_transform(
            next_state[:, i:i+1].numpy()
        ).flatten()
    
    return next_state_denorm

def load_model_and_predict(model_path, new_sequence_data):
    """
    저장된 모델 로딩 및 예측 함수
    """
    # 모델 로딩
    checkpoint = torch.load(model_path)
    
    # 모델 재구성
    lstm_params = checkpoint['lstm_params']
    mlp_params = checkpoint['mlp_params']
    model = BMEDModel(lstm_params, mlp_params)
    model.load_state_dict(checkpoint['model_state_dict'])
    
    # 스케일러 및 설정 로딩
    scalers = checkpoint['scalers']
    feature_cols = checkpoint['feature_cols']
    target_cols = checkpoint['target_cols']
    
    # 예측 수행
    prediction = predict_next_state(model, scalers, new_sequence_data, 
                                  feature_cols, target_cols)
    
    return prediction

# 사용 예시
# 새로운 실험 데이터로 예측
# new_prediction = load_model_and_predict('bmed_lstm_mlp_optimized.pth', new_experiment_data)
# print(f"예측된 다음 상태: {new_prediction}")

# 외삽 테스트 (20% 범위 내)
# 예를 들어, 훈련 데이터 전압 범위가 10-20V였다면
# 새로운 데이터에서 8-24V 범위까지 안정적으로 예측 가능
```

---

# 🔧 구현 전략 & 팁

## 개발 방식
- **순차적 개발**: 함수/클래스 없이 셀별로 코드 작성 → 나중에 리팩토링
- **점진적 테스트**: 각 Phase마다 중간 결과 확인
- **물리 법칙 검증**: Mass balance 제약 조건 지속적 모니터링

## 성능 최적화
- **GPU 활용**: CUDA 사용 가능 시 자동 활용
- **배치 처리**: 메모리 효율성을 위한 적절한 배치 크기
- **Early Stopping**: 과적합 방지를 위한 조기 종료

## 예상 결과
- **모델 정확도**: MSE < 0.01 (정규화된 데이터 기준)
- **물리 법칙 준수**: >95% 물리적 제약 조건 만족
- **최적화 효율**: 200회 시행으로 최적 하이퍼파라미터 발견

---

# 📚 필요 라이브러리
```python
pip install torch pandas scikit-learn optuna matplotlib numpy
```

**예상 소요 시간**: 2-3주  
**핵심 검증 지표**: MSE, 물리 법칙 준수율, Cross-validation 점수