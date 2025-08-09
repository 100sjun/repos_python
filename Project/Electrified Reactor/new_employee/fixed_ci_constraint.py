import pyomo.environ as pyo
import numpy as np

# 모델 생성
model = pyo.ConcreteModel()

# 전역 매개변수
R = 0.082057  # L·atm/(mol·K)
P = 1  # atm
T = 473.15  # K
rhoi = P/R/T  # mol/L

# 성분 정의
comp = ['CO2', 'H2', 'CH4', 'CO', 'H2O']

# 올바른 정의: Ci와 zi를 Var로, rhoi를 Param으로
model.zi = pyo.Var(comp, initialize=0)
model.Ci = pyo.Var(comp, initialize=0)  # Param이 아닌 Var로!
model.rhoi = pyo.Param(initialize=rhoi)

# 초기 zi 값 설정
model.zi['CO2'] = 0.2
model.zi['H2'] = 0.8
model.zi['CH4'] = 0.0
model.zi['CO'] = 0.0
model.zi['H2O'] = 0.0

# Ci 계산을 위한 Constraint 정의
def ci_rule(model, c):
    return model.Ci[c] == model.rhoi * model.zi[c]

model.ci_constraint = pyo.Constraint(comp, rule=ci_rule)

print("=== Constraint 적용 확인 ===")
print("1. 초기 상태 (솔버 실행 전):")
for c in comp:
    print(f"zi[{c}] = {model.zi[c].value}, Ci[{c}] = {model.Ci[c].value}")

print("\n2. 솔버 실행 후 (Constraint 적용됨):")
solver = pyo.SolverFactory('glpk')
result = solver.solve(model, tee=False)

for c in comp:
    print(f"zi[{c}] = {model.zi[c].value}, Ci[{c}] = {model.Ci[c].value}")

print("\n3. zi 값 변경 후:")
model.zi['CO2'] = 0.1
model.zi['H2'] = 0.7
model.zi['CH4'] = 0.1
model.zi['CO'] = 0.05
model.zi['H2O'] = 0.05

print("변경된 zi 값:")
for c in comp:
    print(f"zi[{c}] = {model.zi[c].value}")

print("\n4. 다시 솔버 실행 (Constraint 재적용):")
result = solver.solve(model, tee=False)

for c in comp:
    print(f"zi[{c}] = {model.zi[c].value}, Ci[{c}] = {model.Ci[c].value}")

print("\n5. 수동 계산과 비교:")
for c in comp:
    manual_calc = model.rhoi.value * model.zi[c].value
    constraint_value = model.Ci[c].value
    print(f"{c}: 수동={manual_calc:.6f}, Constraint={constraint_value:.6f}, 일치={abs(manual_calc-constraint_value)<1e-10}")



