# import modules
import subprocess
import pandas as pd
import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel
import ast
from scipy.stats import norm
from IPython.display import clear_output

# MVs controller
step = 0.01
pre_size = 3
max_error = 10000

# compile zdplaskin
def run_prep(inp_path):
    try:
        process = subprocess.Popen(
            'preprocessor.exe',
            stdout = subprocess.DEVNULL,        # ignore outputs
            stderr = subprocess.DEVNULL,        # ignore errors
            stdin = subprocess.PIPE,            # recognize input
            universal_newlines=True
        )
        
        process.stdin.write(inp_path)
        process.stdin.flush()                   # send a data

        while process.poll() is None:           # check the program state, if None, program is still in the run
            process.stdin.write('.\n')
            process.stdin.flush()
    except:
        pass
    print('check the run of preprocessor')
    return process

# Compile exe
def compile_zdp(name):
    compile_command = [
        'gfortran', '-o', name, 'dvode_f90_m.F90', 'zdplaskin_m.F90',
        'run_plasRxn_v2.F90', 'bolsig_x86_64_g.dll'
    ]
    
    try:
        subprocess.run(compile_command)
    except:
        pass
    print('check the compiler')

# Run exe
def run_exe(exe_path):
    try:
        process = subprocess.Popen(
            exe_path,
            stdout = subprocess.PIPE, # Read standard outputs
            stderr = subprocess.PIPE, # Read standard errors
            universal_newlines=True,  # outputs to str variables
            bufsize = 1               # control the size of buffer
        )

        log_flag = False             # The flag for starting log after "Caculation Start!!"
        while True:
            output = process.stdout.readline()
            if not output:
                break
            if "Calculation Start" in output:
                log_flag = True

            if log_flag:
                print(f'\r{output.strip()}           ',end='',flush=True)

            if "PRESS ENTER TO EXIT" in output:
                process.kill()        # forced shutdown
                break
            if "WARNING: BOLSIG+ convergence failed" in output:
                process.kill()        # forced shutdown
                break
    except:
        pass
    return process

# Error calculation
def cal_error(exp_result):
    # Read a result
    species = []
    with open('qt_species_list.txt','r') as file:
        for line in file:
            line = line.rstrip()
            line = line[3:]
            species.append(line)
        file.close()

    df_sp = pd.read_csv('qt_densities.txt', sep=r'\s+', header=0, names=['Time [s]']+species)

    CH4 = (df_sp['CH4'] + df_sp['CH4(V13)'] + df_sp['CH4(V24)'])
    C2H2 = (df_sp['C2H2'] + df_sp['C2H2(V13)']+ df_sp['C2H2(V2)']+ df_sp['C2H2(V5)'])
    C2H4 = (df_sp['C2H4'] + df_sp['C2H4(V1)']+ df_sp['C2H4(V2)'])
    C2H6 = (df_sp['C2H6'] + df_sp['C2H6(V13)']+ df_sp['C2H6(V24)'])
    C3H6 = (df_sp['C3H6'] + df_sp['C3H6(V)'])
    C3H8 = (df_sp['C3H8'] + df_sp['C3H8(V1)'] + df_sp['C3H8(V2)'])
    C4H10 = (df_sp['C4H9H'])
    H2 = df_sp['H2']
    C = df_sp['C']
    H = df_sp['H']
    CH = df_sp['CH']
    CH2 = df_sp['CH2']
    CH3 = df_sp['CH3']
    C2H3 = df_sp['C2H3']
    C2H5 = df_sp['C2H5']
    C3H7 = df_sp['C3H7']
    CH3p = df_sp['CH3^+']
    C2H5p = df_sp['C2H5^+']

    exp = exp_result

    t1 = abs(df_sp['Time [s]']-5.65).argmin()
    t2 = abs(df_sp['Time [s]']-7.27).argmin()
    t3 = abs(df_sp['Time [s]']-10.18).argmin()
    t4 = abs(df_sp['Time [s]']-16.96).argmin()
    t = [t1, t2, t3, t4]

    for i in t:
        sim = []
        for j in range(18):
            sim_XCH4 = (CH4.iloc[0] - CH4.iloc[i])/CH4.iloc[0]*100
            sim_SH2 = 0.5*H2.iloc[i]/(CH4.iloc[0] - CH4.iloc[i])*100
            sim_SC2H6 = 2*C2H6.iloc[i]/(CH4.iloc[0] - CH4.iloc[i])*100
            sim_SC2H4 = 2*C2H4.iloc[i]/(CH4.iloc[0] - CH4.iloc[i])*100
            sim_SC2H2 = 2*C2H2.iloc[i]/(CH4.iloc[0] - CH4.iloc[i])*100
            sim_SC3H8 = 3*C3H8.iloc[i]/(CH4.iloc[0] - CH4.iloc[i])*100
            sim_SC3H6 = 3*C3H6.iloc[i]/(CH4.iloc[0] - CH4.iloc[i])*100
            sim_SC4H10 = 4*C4H10.iloc[i]/(CH4.iloc[0] - CH4.iloc[i])*100
            sim_SH = 0.25*H.iloc[i]/(CH4.iloc[0] - CH4.iloc[i])*100
            sim_SCH = CH.iloc[i]/(CH4.iloc[0] - CH4.iloc[i])*100
            sim_SCH2 = CH2.iloc[i]/(CH4.iloc[0] - CH4.iloc[i])*100
            sim_SCH3 = CH3.iloc[i]/(CH4.iloc[0] - CH4.iloc[i])*100
            sim_SC2H3 = 2*C2H3.iloc[i]/(CH4.iloc[0] - CH4.iloc[i])*100
            sim_SC2H5 = 2*C2H5.iloc[i]/(CH4.iloc[0] - CH4.iloc[i])*100
            sim_SC3H7 = 3*C3H7.iloc[i]/(CH4.iloc[0] - CH4.iloc[i])*100
            sim_SCH3p = CH3p.iloc[i]/(CH4.iloc[0] - CH4.iloc[i])*100
            sim_SC2H5p = 2*C2H5p.iloc[i]/(CH4.iloc[0] - CH4.iloc[i])*100
            newsim = [sim_XCH4,sim_SH2,sim_SC2H6,sim_SC2H4,sim_SC2H2,sim_SC3H8,sim_SC3H6,sim_SC4H10,sim_SH,sim_SCH,sim_SCH2,sim_SCH3,sim_SC2H3,sim_SC2H5,sim_SC3H7,sim_SCH3p,sim_SC2H5p]
            sim.append(newsim)

    err = 0 
    for i in range(len(exp)):
        for j in range(len(exp[i])):
            if j < 9:
                err += ((exp[i][j] - sim[i][j])/exp[i][j])**2
            else:
                err += (exp[i][j] - sim[i][j])**2
    
    return err, sim, df_sp['Time [s]'].iloc[-1] 

# initial variable sset sampling: latin_hypercube sampling
def latin_hypercube_sampling(bounds, n_samples):
    return np.random.uniform(
        low = np.array([b[0] for b in bounds]),
        high = np.array([b[1] for b in bounds]),
        size = (n_samples, len(bounds))
    )

# update kinet.inp
def update_kinet(samples,index):
    with open(f'kinet.inp','r') as file:
        lines = file.readlines()
        file.close()
    new_lines=[]
    for i in range(316):
        new_value = 10**samples[i]
        new_value2 = f'{new_value:.4e}'.replace('e-0','e-').replace('e+00','e0').replace('e+0','e+').replace('e','d')
        new_line = f'$ double precision, parameter :: f{i} = {new_value2}\n'
        new_lines.append(new_line)
    new_inp = lines[:18] + new_lines + lines[334:]
    with open(f'./kinet{index}.inp', 'w') as file:
        file.writelines(new_inp)

# Gaussian Surrogate Model
def GP_Model(X,y):
    # Gaussian Process Surrogate Model
    constant_kernel = ConstantKernel(constant_value=1.0, constant_value_bounds=(1e-5,1e5))
    RBF_kernel = RBF(1.0, length_scale_bounds = (1e-5,1e5))
    kernel = constant_kernel * RBF_kernel

    gp_model = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=10, alpha=1e-6)
    gp_model.fit(X,y)
    return gp_model

# Expected Improvement
def EI(X, gp, y_best, xi):
    mu, sigma = gp.predict(X, return_std=True)
    sigma = sigma
    #improvement = mu - y_best - xi
    improvement = y_best - mu - xi
    Z = improvement / sigma
    ei = improvement * norm.cdf(Z) + sigma * norm.pdf(Z)
    return ei

# Bayesian Optimization for minimization
def BO(gp, bounds, n_iter, X, y, exp_result):
    """Perform Bayesian Optimization using EI for minimization."""
    for i in range(n_iter):
        y_best = np.min(y)  # 현재까지의 최소값
        
        def random_search(bounds, gp, y_best, num_samples=1000, xi=0.01):
            candidate_points = np.random.uniform(
                low = np.array(bounds)[:,0],
                high = np.array(bounds)[:, 1],
                size = (num_samples, len(bounds))
            )

            ei_values = EI(candidate_points, gp, y_best, xi = xi)

            best_index = np.argmax(ei_values)
            return candidate_points[best_index]
        
        x_next = random_search(bounds, gp, y_best)

        # 새로운 관측값 추가 (사용자 정의 시뮬레이션 함수 호출 필요)
        para_index = [12,13,14,19,22,24,34,45,46,80,88,116,117,118,134,211,220,221,223,226,233,268,269,273,276,291,295,299,300,302,304,308,313,315]
        paraset = -21 * np.ones(316)
        paraset[list(range(12))] = 0
        paraset[list(range(90,102))] = 0
        for j in range(len(para_index)):
            paraset[para_index[j]] = x_next[j]

        update_kinet(paraset,i+pre_size)
        inp_path = f'kinet{i+pre_size}.inp\n'
        exe_path = 'run_plasRxn_v2.exe'
        run_prep(inp_path)
        compile_zdp(exe_path)
        run_exe(exe_path)
        y_next,_, t_time = cal_error(exp_result)  # 사용자의 실제 목적 함수

        # 데이터셋 업데이트
        if t_time > 16.96:
            X = np.vstack([X, x_next])
            y = np.append(y, y_next)
            
            if y_next < y_best:
                bounds = [[num * (1-step), num * (1+step)] for num in x_next]

            np.savetxt('X_final.csv', X, delimiter=',', fmt='%f')
            np.savetxt('y_final.csv', y, delimiter=',', fmt='%f')
        else:
            X = np.vstack([X, x_next])
            y = np.append(y, max_error)

        
        # Surrogate 모델 업데이트
        gp = GP_Model(X, y)
        clear_output()
        print(f'현재 iteration: {i+pre_size}, dataset 크기: {len(X)}, 현재 최소값: {y_best}, 이번 y: {y[-1]}\n',)
    return X, y

# initialize parameter samples
df_set = pd.read_csv('parameter_set.csv')
para_index = [12,13,14,19,22,24,34,45,46,80,88,116,117,118,134,211,220,221,223,226,233,268,269,273,276,291,295,299,300,302,304,308,313,315]
problem = {
    'num_bars': len(para_index),
    'bounds' : [[-step,step]] * len(para_index)
}

initial_samples = latin_hypercube_sampling(problem['bounds'],pre_size)

# set the experimental results
result_list = ['XCH4','SH2','SC2H6','SC2H4','SC2H2','SC3H8','SC3H6','SC4','SH','SCH','SCH2','SCH3','SC2H3','SC2H5','SC3H7','SCH3^+','SC2H5^+']
exp_result = [[7.37,60.84,21.1,3.07,2.16,7.83,0.67,1.23,2.09,0,0,0,0,0,0,0,0],
              [7.79,53.13,24.77,3.24,2.38,9.55,0.77,1.59,2.56,0,0,0,0,0,0,0,0],
              [8.73,62.08,19.51,1.59,1.38,8.64,0.56,1.84,2.40,0,0,0,0,0,0,0,0],
              [21.53,43.31,31.17,3.42,2.67,12.71,0.94,2.31,3.46,0,0,0,0,0,0,0,0]]

# pre-trained mode
db_error = []
db_paraset = []

for i in range(len(initial_samples)):
    inp_path = f'kinet{i}.inp\n'
    exe_path = 'run_plasRxn_v2.exe'

    print(f'run {i}')

    paraset = -21 * np.ones(316)
    paraset[list(range(12))] = 0
    paraset[list(range(90,102))] = 0
    for j in range(len(para_index)):
        paraset[para_index[j]] = initial_samples[i][j]

    update_kinet(paraset,i)
    prep_process = run_prep(inp_path)
    prep_process.wait()
    compile_zdp(exe_path)
    exe_process = run_exe(exe_path)
    exe_process.wait()
    error, sim, total_time = cal_error(exp_result)
    if float(total_time) > 16.96:
        db_error.append(error)
        db_paraset.append(initial_samples[i])
    clear_output()

# pre-trained data gathering
with open('db_error.txt', 'w') as file:
    for i in range(len(db_error)):
        file.writelines(str(db_error[i]))
        file.writelines('\n')
    file.close()
with open('db_paraset.txt', 'w') as file:
    for i in range(len(db_paraset)):
        file.writelines(str(db_paraset[i].tolist()))
        file.writelines('\n')
    file.close()

# Load Data
with open('db_error.txt', 'r') as file:
    lines_error = file.readlines()
    file.close()
with open('db_paraset.txt', 'r') as file:
    lines_paraset = file.readlines()
    file.close()
pre_X = np.array([ast.literal_eval(i) for i in lines_paraset])
pre_Y = np.array([float(i[:-1]) for i in lines_error])

# BO run
gp_model = GP_Model(pre_X,pre_Y)
gp_model.predict(pre_X)
bounds = [[-step,step]] * len(para_index)
X_final, y_final = BO(gp_model, bounds, 10,pre_X,pre_Y,exp_result)

np.savetxt('X_final.csv', X_final, delimiter=',', fmt='%f')
np.savetxt('y_final.csv', y_final, delimiter=',', fmt='%f')