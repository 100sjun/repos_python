import pandas as pd
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib as mpl

# 한글 폰트 설정
plt.rcParams['font.family'] = 'NanumGothic'
plt.rcParams['axes.unicode_minus'] = False  # 마이너스 기호 깨짐 방지


df = pd.read_csv('data.csv')
df1 = df[df['label']=='time']

# parameters initialization
Wi = 66.625
VL = 41.253968 # mL# h2 solubility parameter
p00 = 0.2638
p01 = 5.443e-4
p02 = -2.45e-7
p10 = -1.545e-3
p20 = 2.205e-6
p11 = 9.105e-7
mw_h2o = 18 # g/mol, water molweight
mw_h2 = 2 # g/mol, hydrogen molweight

# Concentration initialization
N0 = np.array([
    df1.iloc[0]['N_Gly [mmol]'],    # 글리세롤
    df1.iloc[0]['N_KOH [mmol]'],    # KOH
    df1.iloc[0]['N_H2 [mmol]'],     # 수소
    df1.iloc[0]['N_LA [mmol]'],     # 젖산
    df1.iloc[0]['N_PDO [mmol]'],    # 1,3-PDO
    df1.iloc[0]['N_FA [mmol]'],     # 포름산
    df1.iloc[0]['N_GA [mmol]'],     # 글리콜산
    df1.iloc[0]['N_ME [mmol]'],     # 메탄올
    df1.iloc[0]['N_EG [mmol]']      # 에틸렌글리콜
], dtype=np.float32)

N_exp = np.array([
    df1.iloc[0:]['N_Gly [mmol]'].values,    # 글리세롤 0
    df1.iloc[0:]['N_KOH [mmol]'].values,    # KOH 1
    df1.iloc[0:]['N_H2 [mmol]'].values,     # 수소 2
    df1.iloc[0:]['N_LA [mmol]'].values,     # 젖산 3
    df1.iloc[0:]['N_PDO [mmol]'].values,    # 1,3-PDO 4
    df1.iloc[0:]['N_FA [mmol]'].values,     # 포름산 5
    df1.iloc[0:]['N_GA [mmol]'].values,     # 글리콜산 6
    df1.iloc[0:]['N_ME [mmol]'].values,     # 메탄올 7
    df1.iloc[0:]['N_EG [mmol]'].values      # 에틸렌글리콜 8
], dtype=np.float32)

def ode_system(t, N, k_params):

    wcat = 0.3*0.02 # Pt 함량 0.3 g의 2%

    VL = 41.253968 # mL
    VG = 145 - VL
    C = N/VL # mmol/mL = M
    T = 200.0
    P_H2O = 10**(8.14019-1810.94/(244.485+T))/750.062 # Antoine equation
    vg = 8.314*(T+273.15)/(P_H2O*1e5)
    NG_H2O = 1/(vg*1e3)*VG

    P_H2 = (N[2]/1e3)/(VG/1e6)*8.314*(T+273.15)/1e5
    S_H2 = p00 + (T)*p10 + P_H2*p01 + (T)**2*p20 + P_H2**2*p02 + (T)*P_H2*p11
    WL = Wi - (N[2])*mw_h2/1000 - NG_H2O*mw_h2o/1000
    NL_H2 = S_H2*WL
    NLf_H2 = min(N[2], NL_H2)

    C[2] = NLf_H2/VL    
    # kg1 0 # KGly 1 # KH2 2 # KPDO 3 # kd1 4 # kd2 5 # kd3 6 # kp1 7 # kp2 8

    # 글리세롤 0 # KOH 1 # 수소 2 # 젖산 3 # 1,3-PDO 4 # 포름산 5 # 글리콜산 6 # 메탄올 7 # 에틸렌글리콜 8
    
    # parameter
    k = 10**k_params
    
    theta = 1 + k[1]*C[0] + (k[2]*C[2])**0.5 + k[3]*C[4]
    rg1= k[0]*k[1]*C[0]/(theta**3)*wcat # Glycerol -> DHA + H2
    rd1 = k[4]*wcat # DHA -> PRA + H2O
    rd2 = k[5]*(C[1]**2)*wcat # DHA + 2KOH -> FA + GA + 2H2
    rd3 = k[6]*(C[2]**2)*wcat # DHA + 2H2 -> ME + EG
    rp1 = k[7]*C[1]*wcat # PRA + KOH -> LA
    rp2 = k[8]*k[2]*C[2]/(theta**3)*wcat # PRA + H2 -> PDO

    fd1 = rd1/(rd1+rd2+rd3)
    fd2 = rd2/(rd1+rd2+rd3)
    fd3 = rd3/(rd1+rd2+rd3)

    fp1 = rp1/(rp1+rp2)
    fp2 = rp2/(rp1+rp2)

    dN0dt = -rg1 
    dN1dt = -2*rg1*fd2 - rg1*fd1*fp1
    dN2dt = rg1 + 2*rg1*fd2 - 2*rg1*fd3 - rg1*fd1*fp2
    dN3dt = rg1*fd1*fp1
    dN4dt = rg1*fd1*fp2
    dN5dt = rg1*fd2
    dN6dt = rg1*fd2
    dN7dt = rg1*fd3
    dN8dt = rg1*fd3

    return np.array([dN0dt, dN1dt, dN2dt, dN3dt, dN4dt, dN5dt, dN6dt, dN7dt, dN8dt])


t_eval = np.array(df1.iloc[0:]['Time [h]'].values, dtype=np.float32)

def solve_ode(k_params, t_eval, y0):
    sol = solve_ivp(
        ode_system,
        t_span=(t_eval[0], t_eval[-1]),
        y0=y0,
        t_eval=t_eval,
        args=(k_params,),
        method='RK45',
        rtol=1e-8,
        atol=1e-10
    )
    return sol.y

def objective_function(k_params):
    try:
        N_pred = solve_ode(k_params, t_eval, N0)

        error = np.mean((N_pred - N_exp)**2)
        return error
    except:
        return 1e10
    


best_loss = float('inf')
best_x = None
x0 = np.zeros(9)

print("최적화 시작...")
for i in range(100):  # 100번 반복
    result = minimize(
        objective_function,
        x0=x0,
        method='L-BFGS-B', 
    )
    
    current_loss = result.fun
    improvement = (best_loss - current_loss) / best_loss * 100 if best_loss != float('inf') else 0
    
    print(f"반복 {i+1}/100 - Loss: {current_loss:.6e}", end='')
    if current_loss < best_loss:
        print(f" (개선율: {improvement:.2f}%)")
        best_loss = current_loss
        best_x = result.x.copy()
        x0 = best_x + np.random.randn(9) * 0.1  # 작은 변동폭으로 탐색
    else:
        print(" (개선 없음)")
        x0 = best_x + np.random.randn(9) * 0.5

print("\n최종 최적화 수행 중...")
result = minimize(
    objective_function,
    x0=best_x,
    method='L-BFGS-B',
)
print(f"최종 Loss: {result.fun:.6e}")

# 예측값 계산을 위한 더 조밀한 시간 격자 생성
t_pred = np.linspace(t_eval[0], t_eval[-1], 100)
N_pred = solve_ode(result.x, t_pred, N0)

# 성분별 비교 그래프 그리기
fig, axes = plt.subplots(3, 3, figsize=(15, 15))
components = ['Gly', 'KOH', 'H2', 'LA', 'PDO', 'FA', 'GA', 'ME', 'EG']

for i, (ax, comp) in enumerate(zip(axes.flat, components)):
    ax.plot(t_eval, N_exp[i], 'o', label='실험값')
    ax.plot(t_pred, N_pred[i], '-', label='예측값')
    ax.set_xlabel('시간 [h]')
    ax.set_ylabel(f'N_{comp} [mmol]')
    ax.legend()
    ax.grid(True)

plt.tight_layout()
plt.show()