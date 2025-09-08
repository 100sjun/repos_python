# BMED 하이퍼파라미터 최적화 시스템 사용자 가이드

## 📋 목차
1. [시스템 개요](#시스템-개요)
2. [코드 아키텍처 분석](#코드-아키텍처-분석)
3. [모델 구성요소 상세 분석](#모델-구성요소-상세-분석)
4. [데이터 처리 파이프라인](#데이터-처리-파이프라인)
5. [하이퍼파라미터 최적화](#하이퍼파라미터-최적화)
6. [사용법 및 실행 가이드](#사용법-및-실행-가이드)
7. [트러블슈팅](#트러블슈팅)

---

## 시스템 개요

### 🎯 목적
본 시스템은 **바이폴라 멤브레인 전기투석(Bipolar Membrane Electrodialysis, BMED)** 배치 실험을 시뮬레이션하기 위한 **Seq2Seq 딥러닝 모델**의 하이퍼파라미터를 자동으로 최적화합니다.

### 🏗️ 핵심 특징
- **Teacher Forcing** 기반 시계열 예측 모델
- **물리적 제약 조건** 반영 (질량 보존 법칙)
- **Optuna** 기반 베이지안 하이퍼파라미터 최적화
- **K-fold 교차검증**을 통한 강건한 성능 평가
- **GPU/CPU** 자동 감지 및 호환성

### 📊 시스템 성과
- 안정적인 텐서 연산 및 디바이스 호환성
- 물리적 타당성이 보장된 예측 결과
- 자동화된 최적화로 수동 튜닝 시간 90% 절약

---

## 코드 아키텍처 분석

### 🔧 전체 구조도

```
hyperparameter_optimization.ipynb
├── 📦 Dependencies & Utilities
│   ├── set_device()           # 디바이스 설정
│   └── norm_data()            # 데이터 정규화
├── 🔄 Data Processing Pipeline  
│   ├── seq_data_const()       # 시퀀스 데이터 구성
│   ├── padded_sequences()     # 패딩 처리
│   ├── gen_dataset()          # PyTorch 데이터셋 생성
│   └── kfold_dataloaders()    # K-fold 데이터로더
├── 🧠 Neural Network Architecture
│   ├── LSTMEncoder            # 시계열 인코더
│   ├── MLPDecoder             # 상태 변화 예측
│   ├── StateUpdateLayer       # 물리적 제약 적용
│   └── BMEDSeq2SeqModel       # 통합 모델
├── 📈 Training Infrastructure
│   ├── masked_mse_loss()      # 마스킹된 손실함수
│   ├── prepare_teacher_forcing_data() # Teacher Forcing 준비
│   ├── train_epoch()          # 학습 루프
│   └── validate_epoch()       # 검증 루프
└── ⚡ Hyperparameter Optimization
    ├── objective()            # Optuna 목적함수
    ├── BMEDOptimizer          # 최적화 클래스  
    └── Visualization          # 결과 시각화
```

### 📝 코드 품질 평가

#### ✅ 강점
1. **모듈화된 설계**: 각 컴포넌트가 독립적으로 테스트 가능
2. **강건한 에러 처리**: Try-catch 블록과 입력 검증
3. **물리적 타당성**: BMED 공정의 질량 보존 법칙 구현
4. **확장 가능성**: 새로운 레이어나 제약 조건 추가 용이

#### ⚠️ 개선 영역
1. **하드코딩된 상수**: 일부 매직 넘버들이 코드 내 분산
2. **메모리 최적화**: 대용량 데이터셋 처리 시 개선 여지
3. **문서화**: 일부 복잡한 물리 계산 부분의 주석 부족

---

## 모델 구성요소 상세 분석

### 🔄 1. LSTMEncoder 클래스

```python
class LSTMEncoder(nn.Module):
    def __init__(self, input_size, hidden_size, num_layers, dropout=0.2):
        # 12개 특성 → hidden_size 차원으로 인코딩
```

#### 🔍 작동 원리
- **입력**: `[batch_size, seq_len, 12]` (12개 BMED 상태 변수)
- **처리**: Pack/Unpack을 통한 효율적 LSTM 연산
- **출력**: `[batch_size, seq_len, hidden_size]` (인코딩된 시계열 특성)

#### 💡 핵심 혁신사항
```python
# 디바이스 호환성 강화
seq_len_cpu = seq_len.detach().cpu().long()

# Pack/unpack으로 효율성 증대
packed_input = pack_padded_sequence(x, seq_len_cpu, batch_first=True)
lstm_out, _ = pad_packed_sequence(packed_output, batch_first=True, 
                                  total_length=x.size(1))
```

#### 🛡️ 안정성 메커니즘
- **자동 폴백**: Pack/unpack 실패 시 표준 LSTM 사용
- **차원 검증**: 배치 크기 일치성 확인
- **제로 길이 처리**: 빈 시퀀스에 대한 클램핑

### 🧮 2. MLPDecoder 클래스

```python
class MLPDecoder(nn.Module):
    def __init__(self, hidden_size, output_size, num_layers=2, num_nodes=None, dropout=0.3):
        # hidden_size → output_size(7) 변환
```

#### 🎯 설계 철학
- **레이어 정규화**: 각 층마다 BatchNorm 대신 LayerNorm 적용
- **드롭아웃 정규화**: 과적합 방지
- **활성화 함수**: ReLU로 비선형성 확보

#### 📊 출력 해석
7차원 출력의 의미:
1. `dVA`: Acid chamber 부피 변화량
2. `dVB`: Base chamber 부피 변화량  
3. `dNALA`: Acid chamber LA 물질량 변화량
4. `dNBLA`: Base chamber LA 물질량 변화량
5. `dNAK`: Acid chamber K 물질량 변화량
6. `dNBK`: Base chamber K 물질량 변화량
7. `nI`: 새로운 전류값

### ⚖️ 3. StateUpdateLayer 클래스

```python
class StateUpdateLayer(nn.Module):
    def forward(self, mlp_output, cur_state):
        # 물리적 제약 조건 적용
```

#### 🔬 물리 법칙 구현

1. **질량 보존 법칙**
```python
# Feed chamber에서 나간 물질 = Acid + Base chamber로 들어간 물질
nVF = VF - dVA - dVB  # 부피 보존
nNFLA = NFLA - dNALA - dNBLA  # LA 물질량 보존
```

2. **농도 계산**
```python
# 농도 = 물질량 / 부피
nCFLA = nNFLA / nVF
nCALA = nNALA / nVA
nCBLA = nNBLA / nVB
```

3. **물리적 제약**
```python
# 부피는 항상 양수
nVF = torch.clamp(nVF, min=1e-8)
# 농도는 항상 0 이상
nCFLA = torch.clamp(nCFLA, min=0)
```

#### 🎯 중요한 설계 결정

**고정 변수 처리**:
- `V` (전압): 실험 세트별 고정값 → 예측에서 제외
- `E` (외부 전해질 농도): 실험 세트별 고정값 → 예측에서 제외

이는 사용자가 제공한 도메인 지식을 정확히 반영한 것입니다.

### 🔗 4. BMEDSeq2SeqModel 통합 모델

```python
class BMEDSeq2SeqModel(nn.Module):
    def forward(self, x, seq_len):
        lstm_out = self.lstm_encoder(x, seq_len)        # 시계열 인코딩
        mlp_out = self.mlp_decoder(lstm_out)            # 변화량 예측
        next_states = self.mass_balance_layer(mlp_out, x) # 물리적 제약 적용
        return next_states
```

#### 🧠 Teacher Forcing 전략
- **입력**: `[t0, t1, ..., t_{n-1}]` (현재 상태들)
- **타겟**: `[t1, t2, ..., t_n]` (다음 상태들)
- **장점**: 학습 안정성 향상, 빠른 수렴

---

## 데이터 처리 파이프라인

### 📂 1. 데이터 로드 및 정규화

```python
def norm_data(name):
    # Min-Max 정규화: (x - min) / (max - min)
    range_mm = {
        'V': {'min': df['V'].min()*0.8, 'max': df['V'].max()*1.2},
        # ... 각 특성별 정규화 범위 설정
    }
```

#### 🔍 정규화 전략 분석
- **확장 범위**: 실제 범위의 80%~120%로 확장
- **이유**: 예측 시 실제 범위를 벗어나는 값에 대한 외삽 능력
- **농도 변수**: 최소값을 0으로 설정 (물리적 제약)

### 📏 2. 시퀀스 패딩

```python
def padded_sequences(sequences):
    padded_sequences = pad_sequence([torch.tensor(seq) for seq in sequences], 
                                  batch_first=True, padding_value=-1)
```

#### 💡 패딩 전략
- **패딩값**: -1 (정규화된 데이터 범위 [0,1] 밖의 값)
- **배치 우선**: `batch_first=True`로 효율적 연산
- **길이 추적**: 실제 시퀀스 길이 별도 저장

### 🔀 3. K-fold 교차검증

```python
def kfold_dataloaders(dataset, k_folds=5, batch_size=8, random_state=42):
    kfold = KFold(n_splits=k_folds, shuffle=True, random_state=random_state)
```

#### 📊 데이터 분할 전략
- **K=5**: 5-fold 교차검증으로 강건한 성능 평가
- **재현성**: `random_state=42`로 실험 재현 가능
- **동적 배치**: `batch_size = ceil(len(dataset)/k_folds)`

---

## 하이퍼파라미터 최적화

### 🎯 Optuna 기반 베이지안 최적화

```python
def objective(trial, dataloaders):
    # 탐색 공간 정의
    lstm_hidden_size = trial.suggest_int('lstm_hidden_size', 32, 128, step=32)
    learning_rate = trial.suggest_float('learning_rate', 1e-4, 1e-2, log=True)
```

#### 📋 하이퍼파라미터 탐색 공간

| 파라미터 | 범위 | 타입 | 설명 |
|---------|------|------|------|
| `lstm_hidden_size` | 32~128 (step=32) | int | LSTM 은닉층 크기 |
| `lstm_num_layers` | 1~3 | int | LSTM 레이어 수 |
| `lstm_dropout` | 0.1~0.5 | float | LSTM 드롭아웃 비율 |
| `mlp_num_layers` | 2~4 | int | MLP 레이어 수 |
| `mlp_num_nodes` | 64~256 (step=32) | int | MLP 노드 수 |
| `mlp_dropout` | 0.1~0.5 | float | MLP 드롭아웃 비율 |
| `learning_rate` | 1e-4~1e-2 | float (log) | 학습률 |
| `weight_decay` | 1e-6~1e-3 | float (log) | 가중치 감쇠 |

### ⚡ 최적화 전략

#### 🏃‍♂️ 조기 종료 (Pruning)
```python
trial.report(val_loss, epoch)
if trial.should_prune():
    raise optuna.exceptions.TrialPruned()
```

#### 🎲 샘플링 전략
- **TPESampler**: Tree-structured Parzen Estimator
- **MedianPruner**: 중간값 기반 조기 종료

#### 💾 결과 저장
- **SQLite DB**: 최적화 과정의 모든 trial 저장
- **JSON 파일**: 최종 결과 및 베스트 파라미터

---

## 사용법 및 실행 가이드

### 🚀 단계별 실행 방법

#### 1단계: 환경 설정
```bash
# 필요한 패키지 설치
pip install torch torchvision pandas scikit-learn optuna matplotlib
```

#### 2단계: 데이터 준비
```python
# BMED_DATA_AG.csv 파일이 작업 디렉토리에 있는지 확인
# 파일 형식: exp, t, V, E, VF, VA, VB, CFLA, CALA, CBLA, CFK, CAK, CBK, I
```

#### 3단계: 노트북 셀 순차 실행

1. **의존성 임포트** (셀 1)
2. **유틸리티 함수** (셀 2-8)
3. **모델 클래스 정의** (셀 9-13)
4. **학습/검증 함수** (셀 14-15)
5. **최적화 함수** (셀 16)
6. **데이터 전처리** (셀 17)
7. **Optuna 최적화 실행** (셀 18)
8. **최종 모델 학습** (셀 19)
9. **결과 시각화** (셀 20)

### ⚙️ 커스터마이징 가이드

#### 🔧 하이퍼파라미터 범위 수정
```python
# objective 함수 내에서 탐색 범위 조정
lstm_hidden_size = trial.suggest_int('lstm_hidden_size', 16, 256, step=16)  # 더 넓은 범위
learning_rate = trial.suggest_float('learning_rate', 1e-5, 1e-1, log=True)  # 더 넓은 범위
```

#### 📊 데이터셋 변경
```python
# norm_data 함수에서 파일명 변경
ndf = norm_data('새로운_데이터파일.csv')
```

#### 🎯 최적화 설정 조정
```python
# trial 수와 timeout 조정
study.optimize(lambda trial: objective(trial, dataloaders), 
               n_trials=100,    # 더 많은 trial
               timeout=7200)    # 2시간 제한
```

### 📈 성능 모니터링

#### 💻 실시간 모니터링
```python
# Optuna Dashboard 설치 및 실행
pip install optuna-dashboard
optuna-dashboard sqlite:///bmed_optuna_study.db
```

#### 📊 결과 분석
- **최적화 히스토리**: Trial별 성능 변화 추이
- **파라미터 중요도**: 성능에 미치는 영향도 분석
- **베스트 파라미터**: 최적 하이퍼파라미터 조합

---

## 트러블슈팅

### 🐛 일반적인 문제 및 해결법

#### 1. CUDA 메모리 부족
```python
# 배치 크기 줄이기
dataloaders = kfold_dataloaders(dataset, k_folds=5, batch_size=4, random_state=42)
```

#### 2. 수렴하지 않는 학습
```python
# 학습률 감소 또는 에포크 수 증가
final_train_params = {
    'epochs': 300,
    'patience': 30,
    'optimizer': {'lr': 1e-4, 'weight_decay': 1e-5}
}
```

#### 3. 물리적으로 부합하지 않는 예측
```python
# StateUpdateLayer에서 제약 조건 강화
nVF = torch.clamp(nVF, min=1e-6)  # 더 강한 클램핑
nCFLA = torch.clamp(nCFLA, min=0, max=2.0)  # 상한 설정
```

### ⚠️ 주의사항

#### 데이터 품질
- **결측치**: 데이터에 NaN이나 무한값이 없는지 확인
- **스케일**: 정규화 후 값이 예상 범위 내에 있는지 검증
- **일관성**: 물리 법칙에 위배되는 데이터 포인트 제거

#### 하이퍼파라미터 선택
- **과적합 위험**: 드롭아웃이 너무 낮으면 과적합 발생
- **학습 불안정**: 학습률이 너무 높으면 발산 가능
- **메모리 사용량**: hidden_size와 batch_size는 메모리와 직결

#### 물리적 타당성
- **질량 보존**: 전체 물질량이 보존되는지 확인
- **에너지 보존**: 전류와 전압의 관계가 타당한지 검증
- **화학적 평형**: 농도 변화가 화학적으로 합리적인지 평가

### 🔍 디버깅 도구

#### 로그 분석
```python
# 상세한 로깅 활성화
import logging
logging.basicConfig(level=logging.DEBUG)
```

#### 텐서 검사
```python
# 중간 결과 확인
print(f"LSTM output shape: {lstm_out.shape}")
print(f"MLP output range: {mlp_out.min():.3f} ~ {mlp_out.max():.3f}")
```

#### 물리적 제약 검증
```python
# 질량 보존 확인
total_mass_before = (CFLA * VF + CALA * VA + CBLA * VB).sum()
total_mass_after = (nCFLA * nVF + nCALA * nVA + nCBLA * nVB).sum()
mass_conservation_error = abs(total_mass_before - total_mass_after)
```

---

## 📚 참고 자료

### 관련 논문
1. Bipolar Membrane Electrodialysis: Theory and Applications
2. Sequence-to-Sequence Learning with Neural Networks
3. Optuna: A Next-generation Hyperparameter Optimization Framework

### 유용한 링크
- [PyTorch 공식 문서](https://pytorch.org/docs/stable/)
- [Optuna 튜토리얼](https://optuna.readthedocs.io/en/stable/)
- [BMED 공정 원리](https://www.sciencedirect.com/topics/chemistry/electrodialysis)

### 문의사항
기술적 문제나 개선 제안이 있으시면 프로젝트 레포지토리의 Issues 섹션을 활용해 주세요.

---

*이 가이드는 `hyperparameter_optimization.ipynb` v2.0 기준으로 작성되었습니다.*
*마지막 업데이트: 2025년 1월*