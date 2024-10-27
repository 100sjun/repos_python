import subprocess
import pandas as pd
import scipy.constants as const
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize
import re
import signal
import math
from scipy.stats import norm
from IPython.display import clear_output
import time

def to_fortran_double(number):
    exponent = math.floor(math.log10(number))
    mantissa = number / (10 ** exponent)

    return f"{mantissa:.3f}d{exponent}"

def update_kinet_inp(fortran_doubles, det):
    try:
        with open('kinet.inp', 'r') as file:
            lines = file.readlines()

        # 파라미터 라인을 찾아 수정
        pattern = re.compile(r'\$ double precision, parameter :: (f\d+)\s*=\s*([\d.d+-]+)')
        for i, line in enumerate(lines):
            match = pattern.match(line)
            if match:
                param_name = match.group(1)
                param_number = int(param_name[1:])
                if param_number < len(fortran_doubles):
                    new_value = fortran_doubles[param_number]
                    lines[i] = f'$ double precision, parameter :: {param_name}= {new_value}\n'

        # 수정된 내용을 파일에 다시 쓰기
        if det == 0:
            with open('kinet2.inp', 'w') as file:
                file.writelines(lines)
            print(f"kinet2.inp 파일의 파라미터가 성공적으로 업데이트되었습니다.")    
        if det == 1:
            with open('kinet.inp', 'w') as file:
                file.writelines(lines)
                print(f"kinet.inp 파일의 파라미터가 성공적으로 업데이트되었습니다.")   
        

    except FileNotFoundError:
        print(f"오류: kinet.inp 파일을 찾을 수 없습니다.")
    except IOError as e:
        print(f"파일 처리 중 오류가 발생했습니다: {e}")
    except Exception as e:
        print(f"예상치 못한 오류가 발생했습니다: {e}")

def run_processor(det):
    process = subprocess.Popen(
        'preprocessor.exe',
        stdout = subprocess.DEVNULL,
        stderr = subprocess.DEVNULL,
        stdin = subprocess.PIPE,
        universal_newlines= True
    )

    if det == 0:
        process.stdin.write('kinet2.inp\n')
        process.stdin.flush()
    else:
        process.stdin.write('.\n')
        process.stdin.flush()
        
    while process.poll() is None:
        try:
            process.stdin.write('.\n')
            process.stdin.flush()
        except OSError:
            break

def run_compile(version):
    command = f'gfortran -o run_plasRxn_r{version}.exe dvode_f90_m.F90 zdplaskin_m.F90 run_plasRxn_r{version}.F90 bolsig_x86_64_g.dll'
    result = subprocess.run(command)

def run_exe(exe_path):
    try:
        # 실행 파일 실행
        process = subprocess.Popen(
            exe_path,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
            bufsize=1
        )

        # 실시간으로 출력 읽기
        while True:
            output = process.stdout.readline()
            if "ZDPlasKin INFO: the density of species not configured for BOLSIG+ solver exceeds 1.00D+00" in output:
                print("지정된 문자열이 감지되었습니다. 프로그램을 종료합니다.")
                process.send_signal(signal.SIGTERM)
                break
            if "PRESS ENTER TO EXIT" in output:
                print("지정된 문자열이 감지되었습니다. 프로그램을 종료합니다.")
                process.send_signal(signal.SIGTERM)
                break
            if output:
                print(output.strip())
        clear_output()
        # 프로세스 종료 대기
        return_code = process.wait()

    except:
        pass

def cal_error(version):
    species = []
    with open('qt_species_list.txt','r') as file:
        i = 0
        for line in file:
            line = line.strip()
            if i < 9:
                line = line[2:]
            else:
                line = line[3:]
            species.append(line)
            i += 1
        file.close()

    df_sp = pd.read_csv('qt_densities.txt', sep=r'\s+', header=0, names=['Time [s]']+species)
    
    # mass concentration [g/mL]
    t = df_sp['Time [s]']
    CH4 = (df_sp['CH4'] + df_sp['CH4(V13)'] + df_sp['CH4(V24)'])/const.N_A
    C2H6 = (df_sp['C2H6'] + df_sp['C2H6(V13)'] + df_sp['C2H6(V24)'])/const.N_A
    C2H4 = (df_sp['C2H4'] + df_sp['C2H4(V1)'] + df_sp['C2H4(V2)'])/const.N_A
    C2H2 = (df_sp['C2H2'] + df_sp['C2H2(V13)'] + df_sp['C2H2(V2)'] + df_sp['C2H2(V5)'])/const.N_A
    C3H6 = (df_sp['C3H6'] + df_sp['C3H6(V)'])/const.N_A
    C3H8 = (df_sp['C3H8'] + df_sp['C3H8(V1)'] + df_sp['C3H8(V2)'])/const.N_A
    C4H10 = df_sp['C4H9H']/const.N_A
    H2 = (df_sp['H2'])/const.N_A

    conv = (CH4.iloc[0] - CH4.iloc[-1]) / CH4.iloc[0] * 100
    S_H2 = H2.iloc[-1] / (CH4.iloc[0] - CH4.iloc[-1]) / 2 * 100
    S_C2H6 = 2*C2H6.iloc[-1] / (CH4.iloc[0] - CH4.iloc[-1]) * 100
    S_C2H4 = 2*C2H4.iloc[-1] / (CH4.iloc[0] - CH4.iloc[-1]) * 100
    S_C2H2 = 2*C2H2.iloc[-1] / (CH4.iloc[0] - CH4.iloc[-1]) * 100
    S_C3H8 = 3*C3H8.iloc[-1] / (CH4.iloc[0] - CH4.iloc[-1]) * 100
    S_C3H6 = 3*C3H6.iloc[-1] / (CH4.iloc[0] - CH4.iloc[-1]) * 100 
    S_C4H10 = 4*C4H10.iloc[-1] / (CH4.iloc[0] - CH4.iloc[-1]) * 100     

    if version == 1:
        exp_conv = 10.59278
        exp_C2H6 = 33.69427
        exp_H2 = 19.07643
        exp_C2 = 19.07643
        exp_C4 = 12.7707
        exp_C3 = 0.955414
    else:
        exp_conv = 21.59794
        exp_C2H6 = 42.89809
        exp_H2 = 36.81529
        exp_C2 = 14.61783
        exp_C4 = 19.26752
        exp_C3 = 1.910828
    
    mse_conv = (conv - exp_conv)**2
    mse_H2 = (S_H2 - exp_H2)**2
    mse_C2H6 = (S_C2H6 - exp_C2H6)**2
    mse_C2 = (S_C2H4 + S_C2H2 - exp_C2)**2
    mse_C3 = (S_C3H8 + S_C3H6 - exp_C3)**2
    mse_C4 = (S_C4H10 - exp_C4)**2

    return 10*mse_conv + mse_H2 + mse_C2H6 + mse_C2 + mse_C3 + mse_C4

def objective_function(X,det,x0):
    err_arr = []
    for i in range(len(X)):
        fi = x0 * 10 ** X[i]
        fortran_doubles = list(map(to_fortran_double, fi))
        err = 0

        update_kinet_inp(fortran_doubles,det)
        run_processor(det)

        print(f'1 실행중!')
        run_compile(str(1))
        run_exe('run_plasRxn_r1.exe')
        err += cal_error(1)

        print(f'5 실행중!')
        run_compile(str(5))
        run_exe('run_plasRxn_r5.exe')
        err += cal_error(5)
        
        err_arr.append(err)

    return err_arr

class GaussianProcess:
    def __init__(self, kernel, noise=1e-8):
        self.kernel = kernel
        self.noise = noise

    def fit(self, X, y):
        self.X_train = X
        self.y_train = y
        self.K = self.kernel(self.X_train, self.X_train) + self.noise * np.eye(len(self.X_train))
        self.L = np.linalg.cholesky(self.K)
        self.alpha = np.linalg.solve(self.L.T, np.linalg.solve(self.L, self.y_train))

    def predict(self, X_test):
        K_s = self.kernel(self.X_train, X_test)
        K_ss = self.kernel(X_test, X_test)

        LK = np.linalg.solve(self.L, K_s)
        mu = K_s.T.dot(self.alpha)
        s2 = np.diag(K_ss) - np.sum(LK**2, axis=0)
        return mu, s2
    
def rbf_kernel(X1, X2, l=1.0, sigma_f=1.0):
    sqdist = np.sum(X1**2,1).reshape(-1,1) + np.sum(X2**2,1) - 2 * np.dot(X1, X2.T)
    return sigma_f**2 * np.exp(-0.5 / l**2 * sqdist)

def expected_improvement(X, X_sample, gp, y_sample, xi=0.01):
    mu, sigma = gp.predict(X)
    mu_sample = gp.predict(X_sample)[0]
    
    sigma = sigma.reshape(-1, 1)
    mu_sample_opt = np.max(mu_sample)
    
    with np.errstate(divide='warn'):
        imp = mu - mu_sample_opt - xi
        Z = imp / sigma
        ei = imp * norm.cdf(Z) + sigma * norm.pdf(Z)
        ei[sigma == 0.0] = 0.0
        
    return ei

def propose_location(acquisition, X_sample, Y_sample, gp, bounds, n_restarts=25):
    dim = X_sample.shape[1]
    min_val = 1
    min_x = None
    
    def min_obj(X):
        return -acquisition(X.reshape(1, -1), X_sample, gp, Y_sample).ravel()
    
    for x0 in np.random.uniform(bounds[:, 0], bounds[:, 1], size=(n_restarts, dim)):
        res = minimize(min_obj, x0=x0, bounds=bounds, method='L-BFGS-B')
        if res.fun < min_val:
            min_val = res.fun
            min_x = res.x
                  
    return min_x.reshape(1, -1)

def bayesian_optimization(n_iters, objective_function, bounds, x0):
    dim = bounds.shape[0]
    det = 0
    X_sample = np.random.uniform(bounds[:, 0], bounds[:, 1], size=(1, dim))
    Y_sample = objective_function(X_sample, det, x0)
    
    # Y_sample을 명확하게 numpy 배열로 변환
    if isinstance(Y_sample, list):
        Y_sample = np.array(Y_sample)
    
    gp = GaussianProcess(kernel=rbf_kernel)
    
    for i in range(n_iters):
        gp.fit(X_sample, Y_sample)
        
        # 다음 샘플링 위치 제안
        next_sample = propose_location(expected_improvement, X_sample, Y_sample, gp, bounds)
        
        # 새 샘플에 대한 목적 함수 계산
        
        next_sample_value = objective_function(next_sample, det, x0)[0]
        if next_sample_value < Y_sample.min():
            det = 1
            x0 *= 10**np.array(next_sample)[0]
            X_sample = np.array([[1,1,1,1,1,1,1,1,1,1,1]])
            Y_sample = np.array([objective_function(X_sample, det, x0)[0]])
        else:
            det = 0
            X_sample = np.vstack((X_sample, next_sample))
            Y_sample = np.concatenate([Y_sample, [next_sample_value]])
        print(Y_sample)
        print(f"Iteration {i+1}/{n_iters}, Best error: {Y_sample.min()}")
    with open('best_para.txt', 'w') as f:
        f.write(str(X_sample[Y_sample.argmin()]))
    return X_sample[Y_sample.argmin()]

best = []
x0 = [1e3, 1, 1e13, 1e7, 1e-3, 1, 1e7, 1e7, 1e7,1e-5,1e-5]

bounds = np.array([[-1, 1]] * 11)  # 각 파라미터의 범위를 -10에서 10으로 설정
best_params = bayesian_optimization(1000, objective_function, bounds, x0)
best = best_params
print("Best parameters found:", best_params)

