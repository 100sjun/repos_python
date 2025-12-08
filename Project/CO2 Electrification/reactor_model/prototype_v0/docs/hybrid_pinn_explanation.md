# Hybrid Physics-Informed Neural Network with CasADi

## 개요

이 문서는 Neural Network(NN)와 CasADi FDM(Finite Difference Method) solver를 결합한 하이브리드 모델의 작동 원리를 설명합니다.

## 핵심 질문: "어떻게 동시에 작동하는가?"

### 답변: Forward Pass에서의 순차적 결합

```
실험 데이터 → NN (k 예측) → CasADi (물리 모델 풀이) → Loss 계산 → Backprop (NN 학습)
```

**중요**: CasADi는 미분 가능한 연산이므로 전체 파이프라인에서 gradient가 흐를 수 있습니다.

---

## 1. 전체 파이프라인

### 1.1 학습 과정 (Training Loop)

```python
# 한 번의 학습 iteration
for 실험_데이터 in 훈련_데이터셋:
    # Step 1: 특징 추출
    features = extract_features(실험_데이터.T_measured)
    # features = [T_avg, dT/dx, q, T_var, T_max]

    # Step 2: NN이 열전도도 예측
    k_pred = neural_network(features)  # NN forward pass

    # Step 3: CasADi로 물리 방정식 풀이
    T_pred = casadi_solver(k_pred, boundary_conditions, q_gen)
    # 모든 노드에서 열방정식을 만족하는 온도 분포 계산

    # Step 4: 예측과 측정값 비교
    loss = MSE(T_pred, 실험_데이터.T_measured)  # 모든 노드에서의 오차

    # Step 5: Backpropagation
    loss.backward()  # NN 파라미터에 대한 gradient 계산
    optimizer.step()  # NN 파라미터 업데이트
```

---

## 2. 각 구성요소의 역할

### 2.1 Neural Network (파라미터 예측기)

**입력**: 온도 프로파일의 특징 벡터
```python
features = [
    T_avg,      # 평균 온도
    dT_dx,      # 온도 구배
    q_gen,      # 발열량
    T_var,      # 온도 분산
    T_max       # 최대 온도
]
```

**출력**: 열전도도 `k` (scalar 값)

**네트워크 구조**:
```
Input(5) → Dense(32) → Tanh → Dense(16) → Tanh → Dense(8) → Tanh → Dense(1) → Softplus → k
```

**Softplus**: `k > 0` 보장 (물리적 제약조건)

---

### 2.2 CasADi FDM Solver (물리 방정식 풀이기)

#### 2.2.1 지배 방정식

1차원 정상상태 열전도 방정식:

```
d/dx(k * dT/dx) + q_gen = 0
```

k가 상수라고 가정하면:

```
k * d²T/dx² + q_gen = 0
```

#### 2.2.2 유한차분법 이산화

공간을 N개의 노드로 나눔 (예: N=50)

중앙차분법으로 2차 미분 근사:

```
d²T/dx² ≈ (T[i+1] - 2*T[i] + T[i-1]) / (Δx)²
```

각 내부 노드 i = 1, 2, ..., N-2에서:

```
k * (T[i+1] - 2*T[i] + T[i-1]) / (Δx)² + q_gen = 0
```

이를 residual 형태로:

```
residual[i] = k * (T[i+1] - 2*T[i] + T[i-1]) / (Δx)² + q_gen
```

#### 2.2.3 CasADi 최적화 문제

**변수**: `T = [T[0], T[1], ..., T[N-1]]` (N개 노드의 온도)

**목적함수**: 모든 내부 노드에서 residual 최소화
```python
objective = Σ (residual[i])²  for i = 1, 2, ..., N-2
```

**제약조건**: 경계조건
```python
constraints:
    T[0] = T_left      # 왼쪽 경계
    T[N-1] = T_right   # 오른쪽 경계
```

**CasADi가 푸는 문제**:
```
minimize  Σ (residual[i])²
subject to
    T[0] = T_left
    T[N-1] = T_right
```

#### 2.2.4 CasADi 코드 구조

```python
def solve_physics_model(k, T_left, T_right, q_gen):
    # 1. 변수 정의 (N개 노드의 온도)
    T = ca.MX.sym('T', n_nodes)

    # 2. 목적함수: residual 최소화
    obj = 0
    for i in range(1, n_nodes-1):
        d2T_dx2 = (T[i+1] - 2*T[i] + T[i-1]) / (dx**2)
        residual = k * d2T_dx2 + q_gen
        obj += residual**2

    # 3. 제약조건: 경계조건
    g = [T[0] - T_left, T[-1] - T_right]

    # 4. NLP 문제 정의
    nlp = {'x': T, 'f': obj, 'g': ca.vertcat(*g)}

    # 5. IPOPT solver로 풀기
    solver = ca.nlpsol('solver', 'ipopt', nlp)
    sol = solver(x0=initial_guess, lbg=0, ubg=0)

    return sol['x']  # 최적화된 온도 분포
```

---

## 3. 중요: IPOPT가 매번 실행됩니다!

### 3.1 IPOPT의 역할

**핵심**: 매 forward pass마다 **IPOPT가 물리 방정식을 풉니다**.

```python
# CasADi가 IPOPT를 호출
solver = ca.nlpsol('solver', 'ipopt', nlp)  # ← IPOPT 사용
sol = solver(x0=T0, lbg=0, ubg=0)
```

**IPOPT 호출 횟수**:
- 100 epochs × 70 training samples = **7,000번 IPOPT 실행**
- 각 IPOPT 실행: 50개 노드의 온도를 최적화
- 각 IPOPT는 내부적으로 수십~수백 번 iteration

**이것이 계산 비용이 높은 이유입니다!**

### 3.2 현재 코드의 Gradient 문제

현재 코드는 **gradient가 제대로 연결되지 않습니다**:

```python
# Step 1: NN이 k 예측 (PyTorch)
k_pred_tensor = nn_model(features)  # requires_grad=True
k_pred = k_pred_tensor.item()       # ← 여기서 gradient 연결 끊김!

# Step 2: IPOPT가 물리 방정식 풀이
T_pred = casadi_solver(k_pred, ...)  # ← numpy array 반환

# Step 3: Loss 계산
T_pred_tensor = torch.FloatTensor(T_pred)  # numpy → tensor
loss_data = MSE(T_pred_tensor, T_measured_tensor)

# Step 4: Backprop
loss_data.backward()  # ← gradient가 k_pred_tensor로 안 흐름!
```

**문제**:
- `.item()`을 호출하면 tensor가 scalar로 변환되면서 **computational graph가 끊김**
- CasADi는 PyTorch와 별개로 실행됨
- `loss_data.backward()`를 해도 `k_pred_tensor`로 gradient가 전달 안 됨

### 3.3 그럼 어떻게 학습이 되는가?

코드를 보면 **k를 직접 비교하는 loss**가 있습니다:

```python
# k 직접 비교 loss
loss_k = torch.mean((k_pred_tensor - k_true_tensor)**2)
total_loss = loss_data + 0.01 * loss_k  # ← 이 부분이 학습을 가능하게 함
```

**실제로 일어나는 일**:
1. `loss_data`: gradient 연결 안 됨 ❌ (온도 프로파일 fitting)
2. `loss_k`: gradient 연결 됨 ✅ (k 직접 비교)
3. **결과**: NN은 k_true를 **암기**하는 supervised learning

**즉, 현재 코드는 진짜 Physics-Informed가 아닙니다!**
- 온도 데이터로부터 k를 추론하는 것이 아님
- 단순히 k_true를 외우는 것

---

## 4. 완전 미분 가능 파이프라인이란?

### 4.1 목표

온도 데이터만으로 k를 학습하려면:

```python
# 이상적인 방식
loss = MSE(T_pred, T_measured)  # k_true 없이!
loss.backward()  # ← gradient가 NN까지 흘러야 함
```

이를 위해서는:
```
∂Loss/∂(NN_params) = ∂Loss/∂T × ∂T/∂k × ∂k/∂(NN_params)
                               ↑
                      이 부분이 핵심!
```

### 4.2 문제: ∂T/∂k 계산

IPOPT가 푼 결과 T는 k의 함수입니다:

```
T = IPOPT(k, boundary_conditions)
```

**질문**: k를 조금 변화시키면 T가 얼마나 변하는가? (= ∂T/∂k)

**방법 1: Numerical Differentiation (느림)**
```python
# k를 약간 변화시켜서 T 변화량 측정
T1 = IPOPT(k, ...)
T2 = IPOPT(k + ε, ...)  # ← IPOPT를 한 번 더 실행!
dT_dk ≈ (T2 - T1) / ε
```

**방법 2: Analytical Differentiation (빠름)**
- CasADi는 **implicit function theorem**을 사용해서 ∂T/∂k를 계산 가능
- 추가 IPOPT 실행 없이 gradient 계산 가능!

### 4.3 완전 미분 가능 파이프라인 구현

CasADi와 PyTorch를 연결하는 custom function:

```python
class CasADiPhysicsLayer(torch.autograd.Function):
    """
    CasADi IPOPT solver를 미분 가능한 PyTorch layer로 만들기
    """

    @staticmethod
    def forward(ctx, k_tensor, T_left, T_right, q_gen):
        """
        Forward: IPOPT로 물리 방정식 풀기

        입력: k (PyTorch tensor)
        출력: T (PyTorch tensor)
        """
        k = k_tensor.item()  # scalar로 변환

        # IPOPT 실행 (7,000번 중 1번)
        T_pred = solve_with_ipopt(k, T_left, T_right, q_gen)

        # Backward를 위해 정보 저장
        ctx.save_for_backward(k_tensor)
        ctx.T_pred = T_pred
        ctx.params = (T_left, T_right, q_gen)

        return torch.tensor(T_pred, dtype=torch.float32)

    @staticmethod
    def backward(ctx, grad_output):
        """
        Backward: ∂T/∂k 계산 (CasADi sensitivity)

        Chain rule:
        ∂Loss/∂k = ∂Loss/∂T × ∂T/∂k
                   ↑         ↑
              grad_output  여기서 계산
        """
        k_tensor = ctx.saved_tensors[0]
        T_pred = ctx.T_pred
        T_left, T_right, q_gen = ctx.params

        # CasADi로 ∂T/∂k 계산 (추가 IPOPT 없이!)
        dT_dk = compute_casadi_sensitivity(k_tensor.item(), T_left, T_right, q_gen)
        # dT_dk.shape = (50,)  # 각 노드에서 ∂T[i]/∂k

        # Chain rule 적용
        grad_k = grad_output @ dT_dk  # scalar

        return grad_k, None, None, None


# 사용 방법
physics_layer = CasADiPhysicsLayer.apply

for data in train_data:
    # NN이 k 예측
    k_pred_tensor = nn_model(features)  # requires_grad=True

    # CasADi physics layer (미분 가능!)
    T_pred_tensor = physics_layer(
        k_pred_tensor,
        data['T_left'],
        data['T_right'],
        data['q_gen']
    )
    # ← gradient가 k_pred_tensor로 흐름!

    # Loss: 온도 프로파일만으로!
    loss = MSE(T_pred_tensor, T_measured_tensor)

    # Backward: 제대로 된 gradient
    loss.backward()
    # ∂Loss/∂(NN params) 올바르게 계산됨!
```

### 4.4 CasADi Sensitivity 계산

```python
def compute_casadi_sensitivity(k, T_left, T_right, q_gen):
    """
    CasADi로 ∂T/∂k 계산

    IPOPT가 푼 최적화 문제:
        minimize  Σ (k * d²T/dx² + q_gen)²
        subject to T[0]=T_left, T[N-1]=T_right

    최적해 T*에서 k를 조금 변화시키면?
    → Implicit Function Theorem 사용
    """

    # 1. k를 파라미터로 하는 심볼릭 문제 정의
    k_sym = ca.MX.sym('k')
    T_sym = ca.MX.sym('T', 50)

    # 2. NLP 문제 (k는 파라미터)
    obj, g = build_heat_equation(T_sym, k_sym, T_left, T_right, q_gen)
    nlp = {'x': T_sym, 'p': k_sym, 'f': obj, 'g': g}

    solver = ca.nlpsol('solver', 'ipopt', nlp)

    # 3. 현재 k에서 최적해 구하기
    sol = solver(x0=T_init, p=k)
    T_optimal = sol['x']

    # 4. Sensitivity: ∂T*/∂k
    # CasADi가 자동으로 계산 (implicit differentiation)
    sensitivity = solver.factory('sens', ['x0', 'p'], ['x'])
    dT_dp = sensitivity(x0=T_init, p=k)

    return np.array(dT_dp['x']).flatten()
```

### 4.5 두 방식 비교

| 항목 | 현재 코드 (분리 방식) | 완전 미분 가능 방식 |
|------|---------------------|-------------------|
| **Gradient 연결** | ❌ 끊김 (`.item()`) | ✅ 연결됨 (Custom Function) |
| **학습 방식** | k_true를 암기 | 온도 데이터로부터 학습 |
| **Loss 항** | MSE(T) + MSE(k) | MSE(T)만으로 충분 |
| **Physics-Informed** | ❌ 부분적 | ✅ 완전함 |
| **일반화** | 낮음 (단순 암기) | 높음 (물리 이해) |
| **IPOPT 횟수** | 7,000번 | 7,000번 (동일) |
| **Sensitivity 계산** | ❌ 없음 | ✅ 있음 (추가 비용) |
| **구현 난이도** | 쉬움 | 어려움 |

---

## 5. 학습 과정 시각화

### 5.1 초기 상태 (Epoch 0)

```
실험 데이터: T_measured = [350, 355, 362, ..., 320] K
k_true = 50 W/m·K

NN 예측: k_pred = 25 W/m·K (무작위 초기값)
↓
CasADi: T_pred = [350, 358, 365, ..., 320] K
↓
Loss = MSE(T_pred, T_measured) = 큰 값
↓
Backprop: NN 파라미터 업데이트
```

### 5.2 중간 상태 (Epoch 50)

```
NN 예측: k_pred = 45 W/m·K (개선됨)
↓
CasADi: T_pred = [350, 356, 363, ..., 320] K
↓
Loss = MSE(T_pred, T_measured) = 중간 값
↓
Backprop: 계속 개선
```

### 5.3 수렴 상태 (Epoch 100)

```
NN 예측: k_pred = 49.5 W/m·K (거의 정확)
↓
CasADi: T_pred = [350, 355, 362, ..., 320] K
↓
Loss = MSE(T_pred, T_measured) = 작은 값
```

---

## 6. 코드 실행 흐름 상세 분석

### 6.1 Main Function

```python
def main():
    # 1. 데이터 생성 (실험 데이터 시뮬레이션)
    exp_data = generate_experimental_data(n_samples=100)
    # 각 샘플: {T_measured[50], T_left, T_right, q_gen, k_true}

    # 2. 모델 초기화
    physics_model = HybridPINNModel(n_nodes=50)
    nn_model = ThermalConductivityNN(input_dim=5)

    # 3. 학습
    trainer = HybridTrainer(physics_model, nn_model)
    history = trainer.train(train_data, val_data, n_epochs=100)

    # 4. 평가
    plot_results(history, test_data, physics_model, nn_model)
```

### 6.2 Training Loop 상세

```python
def train(self, train_data, val_data, n_epochs):
    for epoch in range(n_epochs):
        for data in train_data:
            # ========== Forward Pass ==========

            # Step 1: 특징 추출
            features = extract_features(
                data['T_measured'],  # [350, 355, 362, ..., 320]
                data['T_left'],      # 350
                data['T_right'],     # 320
                data['q_gen']        # 500
            )
            # features = [355.2, -30.0, 0.5, 12.4, 370.5]

            # Step 2: NN이 k 예측
            features_tensor = torch.FloatTensor(features)  # shape: (5,)
            k_pred_tensor = nn_model(features_tensor.unsqueeze(0))  # shape: (1, 1)
            # k_pred_tensor = [[45.3]] (with grad)

            k_pred = k_pred_tensor.item()  # 45.3 (scalar)

            # Step 3: CasADi로 물리 방정식 풀이
            T_pred = physics_model.solve_physics_model(
                k=45.3,
                T_left=350,
                T_right=320,
                q_gen=500
            )
            # T_pred = [350, 356, 363, ..., 320] (numpy array, shape: (50,))

            # Step 4: Loss 계산
            T_measured_tensor = torch.FloatTensor(data['T_measured'])  # (50,)
            T_pred_tensor = torch.FloatTensor(T_pred)  # (50,)

            loss_data = torch.mean((T_pred_tensor - T_measured_tensor)**2)
            # loss_data = 25.4 (tensor with grad)

            # Regularization: k도 맞아야 함
            k_true_tensor = torch.FloatTensor([data['k_true']])  # [50.0]
            loss_k = torch.mean((k_pred_tensor - k_true_tensor)**2)
            # loss_k = 22.09

            total_loss = loss_data + 0.01 * loss_k
            # total_loss = 25.4 + 0.2209 = 25.62

            # ========== Backward Pass ==========

            optimizer.zero_grad()
            total_loss.backward()  # NN 파라미터에 대한 gradient 계산
            optimizer.step()       # NN 파라미터 업데이트
```

### 6.3 CasADi Solver 내부

```python
def solve_physics_model(k=45.3, T_left=350, T_right=320, q_gen=500):
    # 공간 이산화: 50개 노드, dx = 1.0/49 = 0.0204 m

    # 변수: T[0], T[1], ..., T[49]
    T = ca.MX.sym('T', 50)

    # 목적함수: 내부 노드 i=1~48에서 열방정식 residual 최소화
    obj = 0
    for i in range(1, 49):
        # 2차 미분 근사
        d2T_dx2 = (T[i+1] - 2*T[i] + T[i-1]) / (0.0204**2)

        # 열방정식 residual
        # 45.3 * d2T_dx2 + 500 = 0
        residual = 45.3 * d2T_dx2 + 500

        obj += residual**2

    # 제약조건: 경계조건
    g = [T[0] - 350,    # T[0] = 350
         T[49] - 320]    # T[49] = 320

    # NLP 문제
    nlp = {
        'x': T,           # 50개 변수
        'f': obj,         # minimize Σ residual²
        'g': ca.vertcat(*g)  # subject to boundary conditions
    }

    # IPOPT로 풀기
    solver = ca.nlpsol('solver', 'ipopt', nlp)

    # 초기값: 선형 보간
    T0 = np.linspace(350, 320, 50)  # [350, 349.4, 348.8, ..., 320]

    # 최적화 실행
    sol = solver(x0=T0, lbg=0, ubg=0)

    # 결과: 열방정식을 만족하는 온도 분포
    return sol['x']  # [350, 356, 363, ..., 320]
```

---

## 7. 주요 특징 및 장단점

### 7.1 장점

1. **물리 법칙 준수**: CasADi가 물리 방정식을 정확히 풀므로 비물리적인 결과 방지
2. **데이터 효율성**: 순수 NN보다 적은 데이터로 학습 가능
3. **해석 가능성**: k 값이 명시적으로 예측되므로 물리적 의미 명확
4. **외삽 성능**: 물리 법칙 덕분에 훈련 범위 밖 데이터에도 robust

### 7.2 단점

1. **계산 비용**: 매 iteration마다 CasADi optimization 필요
2. **수렴 속도**: 순수 NN보다 학습이 느릴 수 있음
3. **Gradient 근사**: 현재 구현은 numerical gradient 사용
4. **수치 안정성**: CasADi solver가 실패하면 학습 중단 가능

---

## 8. 확장 가능성

### 8.1 더 복잡한 물리 모델

```python
# 온도 의존적 열전도도
k(T) = k0 + k1*T + k2*T²

# 비정상 상태
ρ*cp*∂T/∂t = ∂/∂x(k*∂T/∂x) + q_gen

# 2D/3D 확장
∇·(k∇T) + q_gen = 0
```

### 8.2 다중 파라미터 예측

```python
# NN이 여러 물성을 동시에 예측
output = [k, ρ, cp, q_gen]
```

### 8.3 실시간 역문제 해결

```python
# 실시간 센서 데이터 → NN → CasADi → 파라미터 추정
# 산업 공정 모니터링에 활용
```

---

## 9. 핵심 요약

| 항목 | 설명 |
|------|------|
| **NN 역할** | 온도 데이터 → 열전도도 k 예측 |
| **CasADi 역할** | k와 경계조건 → 모든 노드의 온도 계산 (FDM) |
| **학습 목표** | CasADi 예측 온도가 실험 데이터와 일치하도록 NN 학습 |
| **Loss** | MSE(모든 노드에서 예측 vs 측정) + k 정규화 |
| **Gradient** | PyTorch autograd (numerical) 또는 CasADi sensitivity |
| **이산화 횟수** | 매 forward pass마다 50개 노드에서 FDM 수행 |

---

## 10. 실행 및 결과

### 10.1 실행 방법

```bash
cd /home/sjbaek/repos_python/Project/CO2\ Electrification/reactor_model/prototype_v0/src
python hybrid_pinn_casadi.py
```

### 10.2 출력 파일

1. **hybrid_training_history.png**: 학습 곡선 (loss, k error)
2. **hybrid_test_predictions.png**: 테스트 샘플의 온도 프로파일 비교
3. **hybrid_k_prediction.png**: k 예측 정확도 (true vs predicted)

### 10.3 예상 성능

- **k 예측 MAE**: ~2-5 W/m·K
- **온도 예측 RMSE**: ~1-3 K
- **R² score**: >0.95

---

## 참고자료

- CasADi Documentation: https://web.casadi.org/
- Physics-Informed Neural Networks (PINNs): Raissi et al., 2019
- Hybrid Modeling: Combining data-driven and physics-based approaches

---

# 11. 두 코드 비교: 분리 방식 vs 완전 미분 가능 방식

## 11.1 코드 파일 비교

| 항목 | hybrid_pinn_casadi.py | fully_differentiable_pinn.py |
|------|----------------------|------------------------------|
| **파일명** | 분리 방식 (Detached Mode) | 완전 미분 가능 방식 (Fully Differentiable) |
| **Gradient 연결** | ❌ 끊김 | ✅ 완전 연결 |
| **학습 방식** | Supervised (k_true 필요) | Physics-Informed (k_true 불필요) |
| **구현 난이도** | 쉬움 | 중간 |

---

## 11.2 핵심 차이점 상세 비교

### 차이점 1: Gradient 연결 방식

#### 분리 방식 (hybrid_pinn_casadi.py)

```python
# NN 예측
k_pred_tensor = nn_model(features)  # requires_grad=True

# ❌ Gradient graph 단절!
k_pred = k_pred_tensor.item()  # tensor → scalar 변환

# CasADi 실행 (PyTorch와 분리)
T_pred = physics_model.solve_physics_model(k_pred, ...)  # numpy array

# Loss
T_pred_tensor = torch.FloatTensor(T_pred)  # numpy → tensor
loss_data = MSE(T_pred_tensor, T_measured_tensor)

# ❌ Backward 시 gradient가 NN에 전달 안 됨!
loss_data.backward()
```

**문제**:
- `.item()`으로 computational graph 끊김
- `loss_data`에서 backward 해도 `k_pred_tensor`로 gradient 안 흐름

**해결책**:
```python
# k를 직접 비교하는 loss 추가
loss_k = MSE(k_pred_tensor, k_true_tensor)  # ← gradient 연결됨!
total_loss = loss_data + 0.01 * loss_k
total_loss.backward()  # loss_k를 통해 NN 학습
```

**결과**: k_true를 암기하는 supervised learning

---

#### 완전 미분 가능 방식 (fully_differentiable_pinn.py)

```python
# NN 예측
k_pred_tensor = nn_model(features)  # requires_grad=True

# ✅ Custom PyTorch Function으로 gradient 연결 유지!
physics_layer = CasADiPhysicsLayer.apply
T_pred_tensor = physics_layer(
    k_pred_tensor,  # ← tensor 그대로 전달
    physics_solver,
    T_left, T_right, q_gen
)

# Loss: 온도 데이터만으로!
loss = MSE(T_pred_tensor, T_measured_tensor)

# ✅ Backward 시 gradient가 NN에 제대로 전달됨!
loss.backward()
```

**핵심**:
```python
class CasADiPhysicsLayer(torch.autograd.Function):
    @staticmethod
    def forward(ctx, k_tensor, ...):
        # IPOPT 실행
        k = k_tensor.item()
        T_optimal = solve_with_ipopt(k, ...)
        ctx.save_for_backward(k_tensor, ...)
        return torch.tensor(T_optimal)

    @staticmethod
    def backward(ctx, grad_output):
        # ✅ CasADi로 ∂T/∂k 계산!
        dT_dk = compute_sensitivity(...)
        grad_k = grad_output @ dT_dk
        return grad_k, None, ...
```

**결과**: 온도 프로파일로부터 k를 추론하는 진짜 Physics-Informed learning

---

### 차이점 2: Loss Function

#### 분리 방식

```python
# Loss 구성
loss_data = MSE(T_pred, T_measured)  # 온도 fitting (gradient 끊김)
loss_k = MSE(k_pred, k_true)         # k 직접 비교 (gradient 연결)

total_loss = loss_data + 0.01 * loss_k

# 문제점
# 1. k_true가 필요함 (실험에서 알 수 없는 경우 사용 불가)
# 2. NN이 k_true를 단순 암기
# 3. 물리적 관계를 학습하지 못함
```

**데이터 요구사항**:
```python
exp_data = {
    'T_measured': [...],  # 온도 측정값
    'k_true': 50.0        # ← 필수!
}
```

---

#### 완전 미분 가능 방식

```python
# Loss 구성
loss = MSE(T_pred, T_measured)  # 온도 fitting만! (gradient 연결됨)

# 장점
# 1. k_true 불필요! (실제 실험 데이터에 적용 가능)
# 2. NN이 온도-k 관계를 물리 법칙을 통해 학습
# 3. 일반화 성능 우수
```

**데이터 요구사항**:
```python
exp_data = {
    'T_measured': [...],  # 온도 측정값만!
    # k_true 불필요!
}
```

---

### 차이점 3: Sensitivity 계산

#### 분리 방식

**Sensitivity 없음**:
- ∂T/∂k를 계산하지 않음
- CasADi와 PyTorch가 완전히 분리되어 실행
- Gradient 근사 없음

---

#### 완전 미분 가능 방식

**CasADi Sensitivity 계산**:

```python
def compute_sensitivity(k, T_left, T_right, q_gen, T_optimal):
    """
    CasADi로 ∂T/∂k 계산

    Implicit Function Theorem:
    If F(T*, k) = 0 (optimality condition), then:
    ∂T*/∂k = -(∂F/∂T)^(-1) × (∂F/∂k)
    """

    # 1. 열방정식 residual 정의
    F = heat_equation_residuals(T_sym, k_sym, ...)

    # 2. Jacobian 계산
    dF_dT = ca.jacobian(F, T_sym)  # ∂F/∂T (50×50 행렬)
    dF_dk = ca.jacobian(F, k_sym)  # ∂F/∂k (50×1 벡터)

    # 3. 현재 최적해에서 평가
    dF_dT_val = evaluate(dF_dT, T_optimal, k, ...)
    dF_dk_val = evaluate(dF_dk, T_optimal, k, ...)

    # 4. 선형 시스템 풀이: dF_dT × dT_dk = -dF_dk
    dT_dk = -np.linalg.solve(dF_dT_val, dF_dk_val)

    return dT_dk  # shape: (50,)
```

**계산 비용**:
- Jacobian 계산: 자동 미분 (빠름)
- 선형 시스템 풀이: 50×50 행렬 역행렬 (~0.1ms)
- IPOPT 재실행 불필요!

---

### 차이점 4: 학습 과정 비교

#### 분리 방식의 학습

```
Epoch 1:
  Sample 1: k_true=50.0
    NN 예측: k_pred=25.0
    IPOPT: T_pred = solve(k=25.0) → [350, 358, ...]
    Loss: MSE(T_pred, T_measured) + MSE(25.0, 50.0)
           ↑ gradient 안 흐름      ↑ gradient 흐름

    Backward: ∂Loss/∂(NN params) ← MSE(k_pred, k_true)만 사용

  NN 업데이트: k_true=50.0을 향해 직접 이동

Epoch 50:
  k_pred → 48.5 (k_true에 근접, 단순 암기)
```

**학습 특성**:
- k_true를 목표로 직접 학습
- 온도 프로파일은 검증용으로만 사용
- Supervised learning과 동일

---

#### 완전 미분 가능 방식의 학습

```
Epoch 1:
  Sample 1: (k_true는 모름, 온도 데이터만 있음)
    NN 예측: k_pred=25.0
    IPOPT: T_pred = solve(k=25.0) → [350, 358, ...]
    Loss: MSE(T_pred, T_measured)
           ↑ gradient 흐름!

    Backward:
      ∂Loss/∂T = 2(T_pred - T_measured) / n  # PyTorch
      ∂T/∂k = compute_sensitivity(...)       # CasADi
      ∂k/∂(NN params) = NN.backward()        # PyTorch

      Chain rule: ∂Loss/∂(NN params) = ∂Loss/∂T × ∂T/∂k × ∂k/∂(NN params)

    NN 업데이트: 온도 예측을 개선하는 방향으로 k 조정

Epoch 50:
  k_pred → 49.2 (온도 프로파일을 통해 간접 학습)
```

**학습 특성**:
- k는 온도 프로파일을 맞추기 위한 수단
- 물리 법칙(열방정식)을 통해 T-k 관계 학습
- Physics-Informed learning

---

### 차이점 5: 일반화 성능

#### 분리 방식

**훈련 범위 내**:
- k_true를 암기했으므로 정확함
- 온도 예측도 정확

**훈련 범위 밖 (외삽)**:
```python
# 훈련: k ∈ [10, 80]
# 테스트: k = 100 (훈련 범위 밖)

# NN은 k=100을 본 적 없음
# → k를 잘못 예측 (예: k_pred=75)
# → 온도 예측 부정확
```

**문제**: 새로운 k 값에 대한 일반화 약함

---

#### 완전 미분 가능 방식

**훈련 범위 내**:
- 물리 법칙을 학습했으므로 정확함

**훈련 범위 밖 (외삽)**:
```python
# 훈련: k ∈ [10, 80]
# 테스트: k = 100 (훈련 범위 밖)

# NN은 열방정식을 통해 T-k 관계를 학습
# → 온도 패턴을 보고 k를 추론
# → k=100에 가까운 값 예측 가능 (예: k_pred=95)
# → 온도 예측 더 정확
```

**장점**: 물리 법칙 덕분에 외삽 성능 우수

---

## 11.3 계산 비용 비교

### IPOPT 호출 횟수 (동일)

두 방식 모두:
- Training: 100 epochs × 70 samples = **7,000번**
- Validation: 100 epochs × 15 samples = **1,500번**
- **총 8,500번 IPOPT 실행**

### 추가 계산 비용

#### 분리 방식
- Sensitivity 계산: **없음**
- 추가 비용: **0%**

#### 완전 미분 가능 방식
- Sensitivity 계산: **7,000번** (training만)
- 각 sensitivity: ~0.1-0.5ms (Jacobian + 선형 시스템)
- 추가 비용: **~10-20%**

### 총 계산 시간 비교

50 노드, 100 epochs 기준:

| 방식 | IPOPT | Sensitivity | 총 시간 | 상대 비용 |
|------|-------|-------------|---------|----------|
| 분리 방식 | ~5분 | 0초 | **~5분** | 100% |
| 완전 미분 가능 | ~5분 | ~30초 | **~5.5분** | 110% |

**결론**: 추가 비용은 적지만, 성능 향상이 큼!

---

## 11.4 사용 시나리오

### 분리 방식 (hybrid_pinn_casadi.py) 추천 상황

✅ **k_true를 알고 있는 경우**
- 실험에서 k를 직접 측정한 데이터
- Synthetic data로 proof-of-concept

✅ **빠른 프로토타이핑**
- 구현이 간단함
- 디버깅 쉬움

✅ **Supervised learning이 목적**
- k 예측이 주 목표
- 물리 법칙은 검증용

---

### 완전 미분 가능 방식 (fully_differentiable_pinn.py) 추천 상황

✅ **k_true를 모르는 실제 실험 데이터**
- 온도 센서 데이터만 있음
- k를 역산하고 싶은 경우 (inverse problem)

✅ **물리 법칙을 활용한 학습**
- 데이터가 적을 때
- 외삽 성능이 중요할 때

✅ **진짜 Physics-Informed NN**
- 물리 지식을 학습에 활용
- 학술 연구, 논문 작성

---

## 11.5 실행 및 비교

### 분리 방식 실행

```bash
cd src/
python hybrid_pinn_casadi.py
```

**출력**:
- `hybrid_training_history.png`
- `hybrid_test_predictions.png`
- `hybrid_k_prediction.png`

**예상 결과**:
- k 예측 MAE: ~2-3 W/m·K (k_true로 직접 학습)
- 온도 RMSE: ~1-2 K

---

### 완전 미분 가능 방식 실행

```bash
cd src/
python fully_differentiable_pinn.py
```

**출력**:
- `fully_diff_training_history.png`
- `fully_diff_test_predictions.png`
- `fully_diff_k_prediction.png`

**예상 결과**:
- k 예측 MAE: ~3-5 W/m·K (온도로부터 간접 학습)
- 온도 RMSE: ~1-2 K (더 정확한 물리 이해)

---

## 11.6 코드 구조 비교

### 분리 방식

```
hybrid_pinn_casadi.py
├── ThermalConductivityNN
├── HybridPINNModel
│   └── solve_physics_model()  ← IPOPT 실행, numpy 반환
├── HybridTrainer
│   └── train()
│       ├── k_pred = nn(features).item()  ← gradient 끊김
│       ├── T_pred = solve(k_pred)
│       └── loss = MSE(T) + MSE(k_pred, k_true)
└── generate_experimental_data()
```

### 완전 미분 가능 방식

```
fully_differentiable_pinn.py
├── ThermalConductivityNN
├── CasADiPhysicsSolver
│   ├── solve()                    ← IPOPT 실행
│   └── compute_sensitivity()      ← ∂T/∂k 계산
├── CasADiPhysicsLayer             ← Custom PyTorch Function
│   ├── forward()                  ← IPOPT 호출
│   └── backward()                 ← Sensitivity 사용
└── train_fully_differentiable()
    ├── k_pred_tensor = nn(features)  ← tensor 유지
    ├── T_pred_tensor = physics_layer(k_pred_tensor)
    └── loss = MSE(T_pred, T_measured)  ← k_true 불필요!
```

---

## 11.7 핵심 요약

| 측면 | 분리 방식 | 완전 미분 가능 방식 |
|------|----------|-------------------|
| **Gradient 연결** | ❌ 끊김 | ✅ 연결 |
| **k_true 필요** | ✅ 필수 | ❌ 불필요 |
| **학습 방식** | Supervised | Physics-Informed |
| **일반화 성능** | 낮음 (암기) | 높음 (물리 이해) |
| **계산 비용** | 100% | 110% |
| **구현 난이도** | 쉬움 | 중간 |
| **실제 활용** | Proof-of-concept | 실제 inverse problem |
| **외삽 성능** | 약함 | 강함 |

**권장사항**:
- 학습 목적, proof-of-concept → **분리 방식**
- 실제 데이터 분석, inverse problem → **완전 미분 가능 방식**
