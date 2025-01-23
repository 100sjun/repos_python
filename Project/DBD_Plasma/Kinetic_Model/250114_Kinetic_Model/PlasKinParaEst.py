# Modules
import subprocess
from IPython.display import clear_output
import pandas as pd
import re
import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel
import ast
from scipy.stats import norm

# Hyperparameters
step = 0.01
max_error = 1000
pre_size = 1

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
def compile_zdp1(name):
    compile_command = [
        'gfortran', '-o', name, 'dvode_f90_m.F90', 'zdplaskin_m.F90',
        'run_plasRxn1.F90', 'bolsig_x86_64_g.dll'
    ]
    
    try:
        subprocess.run(compile_command)
    except:
        pass
    print('check the compiler_run1')

def compile_zdp2(name):
    compile_command = [
        'gfortran', '-o', name, 'dvode_f90_m.F90', 'zdplaskin_m.F90',
        'run_plasRxn2.F90', 'bolsig_x86_64_g.dll'
    ]
    
    try:
        subprocess.run(compile_command)
    except:
        pass
    print('check the compiler_run2')

# Compile exe
def compile_zdp3(name):
    compile_command = [
        'gfortran', '-o', name, 'dvode_f90_m.F90', 'zdplaskin_m.F90',
        'run_plasRxn3.F90', 'bolsig_x86_64_g.dll'
    ]
    
    try:
        subprocess.run(compile_command)
    except:
        pass
    print('check the compiler_run3')

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
    conditions = []
    with open('qt_conditions_list.txt','r') as file:
        for line in file:
            line = line.strip()
            line = line[2:]
            conditions.append(line)
        file.close()

    species = []
    with open('qt_species_list.txt','r') as file:
        for line in file:
            line = line.rstrip()
            line = line[3:]
            species.append(line)
        file.close()

    reactions = []
    reaction_list = pd.read_csv('parameter_set.csv')
    reactions = reaction_list['Reaction'].to_list()
    df_cd = pd.read_csv('qt_conditions.txt', sep=r'\s+', header=0, names=['Time [s]']+conditions)
    df_sp = pd.read_csv('qt_densities.txt', sep=r'\s+', header=0, names=['Time [s]']+species)
    df_rx = pd.read_csv('qt_rates.txt', sep=r'\s+', header=0, names=['Time [s]']+reactions)
    top_rate = df_rx.iloc[:,1:].sum().sort_values(ascending=False)

    CH4 = (df_sp['CH4'] + df_sp['CH4(V13)'] + df_sp['CH4(V24)'])
    C2H2 = (df_sp['C2H2'] + df_sp['C2H2(V13)']+ df_sp['C2H2(V2)']+ df_sp['C2H2(V5)'])
    C2H4 = (df_sp['C2H4'] + df_sp['C2H4(V1)']+ df_sp['C2H4(V2)'])
    C2H6 = (df_sp['C2H6'] + df_sp['C2H6(V13)']+ df_sp['C2H6(V24)'])
    C3H6 = (df_sp['C3H6'] + df_sp['C3H6(V)'])
    C3H8 = (df_sp['C3H8'] + df_sp['C3H8(V1)'] + df_sp['C3H8(V2)'])
    C4H10 = (df_sp['C4H9H'])
    C5H12 = (df_sp['C5H12'])
    H2 = df_sp['H2']
    C = df_sp['C']
    H = df_sp['H']
    CH = df_sp['CH']
    CH2 = df_sp['CH2']
    CH3 = df_sp['CH3']
    C2H3 = df_sp['C2H3']
    C2H5 = df_sp['C2H5']
    C3H5 = df_sp['C3H5']
    C3H7 = df_sp['C3H7']
    C4H9 = df_sp['C4H9']

    exp = exp_result
    
    sim_XCH4 = (CH4.iloc[0] - CH4.iloc[-1])/CH4.iloc[0]*100
    sim_SH2 = 0.5*H2.iloc[-1]/(CH4.iloc[0] - CH4.iloc[-1])*100
    sim_SC2H6 = 2*C2H6.iloc[-1]/(CH4.iloc[0] - CH4.iloc[-1])*100
    sim_SC2H4 = 2*C2H4.iloc[-1]/(CH4.iloc[0] - CH4.iloc[-1])*100
    sim_SC2H2 = 2*C2H2.iloc[-1]/(CH4.iloc[0] - CH4.iloc[-1])*100
    sim_SC3H8 = 3*C3H8.iloc[-1]/(CH4.iloc[0] - CH4.iloc[-1])*100
    sim_SC3H6 = 3*C3H6.iloc[-1]/(CH4.iloc[0] - CH4.iloc[-1])*100
    sim_SC4H10 = 4*C4H10.iloc[-1]/(CH4.iloc[0] - CH4.iloc[-1])*100
    sim_SC5H12 = 5*C5H12.iloc[-1]/(CH4.iloc[0] - CH4.iloc[-1])*100
    sim_SH = 0.25*H.iloc[-1]/(CH4.iloc[0] - CH4.iloc[-1])*100
    sim_SCH = CH.iloc[-1]/(CH4.iloc[0] - CH4.iloc[-1])*100
    sim_SCH2 = CH2.iloc[-1]/(CH4.iloc[0] - CH4.iloc[-1])*100
    sim_SCH3 = CH3.iloc[-1]/(CH4.iloc[0] - CH4.iloc[-1])*100
    sim_SC2H3 = 2*C2H3.iloc[-1]/(CH4.iloc[0] - CH4.iloc[-1])*100
    sim_SC2H5 = 2*C2H5.iloc[-1]/(CH4.iloc[0] - CH4.iloc[-1])*100
    sim_SC3H5 = 3*C3H5.iloc[-1]/(CH4.iloc[0] - CH4.iloc[-1])*100
    sim_SC3H7 = 3*C3H7.iloc[-1]/(CH4.iloc[0] - CH4.iloc[-1])*100
    sim_SC4H9 = 4*C4H9.iloc[-1]/(CH4.iloc[0] - CH4.iloc[-1])*100

    sim = []
    sim.append(sim_XCH4)
    sim.append(sim_SH2)
    sim.append(sim_SC2H6)
    sim.append(sim_SC2H4)
    sim.append(sim_SC2H2)
    sim.append(sim_SC3H8)
    sim.append(sim_SC3H6)
    sim.append(sim_SC4H10)
    sim.append(sim_SH)
    sim.append(sim_SCH)
    sim.append(sim_SCH2)
    sim.append(sim_SCH3)
    sim.append(sim_SC2H3)
    sim.append(sim_SC2H5)
    sim.append(sim_SC3H5)
    sim.append(sim_SC3H7)
    sim.append(sim_SC4H9)
    
    err = 0
    
    for i in range(len(sim)):
        if i < 8:
            err += ((exp[i] - sim[i])/exp[i])**2
        else:
            err += (exp[i] - sim[i])**2
    
    return err, top_rate, sim, df_cd['Time [s]'].iloc[-1], sim_XCH4, sim_SC2H6, sim_SC2H4, sim_SC2H2

# initial variable sset sampling: latin_hypercube sampling
def latin_hypercube_sampling(bounds, n_samples):
    return np.random.uniform(
        low = np.array([b[0] for b in bounds]),
        high = np.array([b[1] for b in bounds]),
        size = (n_samples, len(bounds))
    )

# update kinet.inp
def update_kinet1(samples,index):
    with open(f'kinet_blank.inp','r') as file:
        lines = file.readlines()
        file.close()
    new_lines=[]
    for i in range(316):
        new_value = 10**samples[i]
        new_value2 = f'{new_value:.4e}'.replace('e-0','e-').replace('e+00','e0').replace('e+0','e+').replace('e','d')
        new_line = f'$ double precision, parameter :: f{i} = {new_value2}\n'
        new_lines.append(new_line)
    new_inp = lines[:18] + new_lines + lines[334:]
    with open(f'./kinet1_{index}.inp', 'w') as file:
        file.writelines(new_inp)

def update_kinet2(samples,index):
    with open(f'kinet_blank.inp','r') as file:
        lines = file.readlines()
        file.close()
    new_lines=[]
    for i in range(316):
        new_value = 10**samples[i]
        new_value2 = f'{new_value:.4e}'.replace('e-0','e-').replace('e+00','e0').replace('e+0','e+').replace('e','d')
        new_line = f'$ double precision, parameter :: f{i} = {new_value2}\n'
        new_lines.append(new_line)
    new_inp = lines[:18] + new_lines + lines[334:]
    with open(f'./kinet2_{index}.inp', 'w') as file:
        file.writelines(new_inp)

def update_kinet3(samples,index):
    with open(f'kinet_blank.inp','r') as file:
        lines = file.readlines()
        file.close()
    new_lines=[]
    for i in range(316):
        new_value = 10**samples[i]
        new_value2 = f'{new_value:.4e}'.replace('e-0','e-').replace('e+00','e0').replace('e+0','e+').replace('e','d')
        new_line = f'$ double precision, parameter :: f{i} = {new_value2}\n'
        new_lines.append(new_line)
    new_inp = lines[:18] + new_lines + lines[334:]
    with open(f'./kinet3_{index}.inp', 'w') as file:
        file.writelines(new_inp)

def GP_Model(X,y):
    avg_pre_Y = sum(y)/len(y)
    X_range = X.max(axis=0) - X.min(axis=0)
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

        x1 = []
        x2 = []
        x3 = []
        for j in range(int(len(x_next)/3)):
            para1 = x_next[3*j] + x_next[3*j+1] * 10.566 + x_next[3*j+2] * 2.8274
            para2 = x_next[3*j] + x_next[3*j+1] * 11.3 + x_next[3*i+2] * 4.7124
            para3 = x_next[3*j] + x_next[3*j+1] * 4.823 + x_next[3*i+2] * 3.7699
            x1.append(para1)
            x2.append(para2)
            x3.append(para3)
        para_index = [12,13,14,22,24,34,46,116,134,220,223,226,273,276,313,315]
        print(f'run {i+pre_size}')

        inp_path1 = f'kinet1_{i+pre_size}.inp\n'
        exe_path1 = f'run_plasRxn1.exe'

        paraset1 = -21 * np.ones(316)
        paraset1[list(range(12))] = 0
        paraset1[list(range(90,102))] = 0
        for j in range(len(para_index)):
            paraset1[para_index[j]] = x1[j]

        print(f'run {i+pre_size}')

        update_kinet1(paraset1,i+pre_size)
        prep_process = run_prep(inp_path1)
        prep_process.wait()
        compile_zdp1(exe_path1)
        exe_process = run_exe(exe_path1)
        exe_process.wait()
        y1,_,_,t_time1,_,_,_,_ = cal_error(exp_result[0])
        print(f'run {i+pre_size} set1 completion')

        inp_path2 = f'kinet2_{i+pre_size}.inp\n'
        exe_path2 = f'run_plasRxn2.exe'

        paraset2 = -21 * np.ones(316)
        paraset2[list(range(12))] = 0
        paraset2[list(range(90,102))] = 0
        for j in range(len(para_index)):
            paraset2[para_index[j]] = x2[j]

        update_kinet2(paraset2,i+pre_size)
        prep_process = run_prep(inp_path2)
        prep_process.wait()
        compile_zdp2(exe_path2)
        exe_process = run_exe(exe_path2)
        exe_process.wait()
        y2,_,_,t_time2,_,_,_,_ = cal_error(exp_result[1])
        print(f'run {i+pre_size} set2 completion')

        inp_path3 = f'kinet3_{i+pre_size}.inp\n'
        exe_path3 = f'run_plasRxn3.exe'

        paraset3 = -21 * np.ones(316)
        paraset3[list(range(12))] = 0
        paraset3[list(range(90,102))] = 0
        for j in range(len(para_index)):
            paraset3[para_index[j]] = x3[j]

        update_kinet3(paraset3,i+pre_size)
        prep_process = run_prep(inp_path3)
        prep_process.wait()
        compile_zdp3(exe_path3)
        exe_process = run_exe(exe_path3)
        exe_process.wait()
        y3,_,_,t_time3,_,_,_,_ = cal_error(exp_result[2])
        print(f'run {i+pre_size} set3 completion')
       
        y_next = y1 + y2 + y3
      
        # 데이터셋 업데이트
        if (t_time1 > 2.1) & (t_time2 > 3.5) & (t_time3 > 5.6):
            X = np.vstack([X, x_next])
            y = np.append(y, y_next)
            
            bnds_low = [i[0] for i in bounds]
            bnds_high = [i[1] for i in bounds]
            if y_next < y_best:
                for j in range(len(bounds)):
                    bnds_low[j] = x_next[j] * (1-step)
                    bnds_high[j] = x_next[j] + (1+step)
            bounds = [[bnds_low[j],bnds_high[j]] for j in range(len(bnds_low))]

            np.savetxt('X_final.csv', X, delimiter=',', fmt='%f')
            np.savetxt('y_final.csv', y, delimiter=',', fmt='%f')
        else:
            X = np.vstack([X, x_next])
            y = np.append(y, max_error)

        # Surrogate 모델 업데이트
        gp = GP_Model(X, y)
        clear_output()
        print(f'현재 iteration: {i+pre_size}, dataset 크기: {len(X)}, 현재 최소값: {y_best}, 이번 y: {y[-1]}\n')
    return X, y

# Optimization
df_set = pd.read_csv('parameter_set.csv')
para_index = [12,13,14,22,24,34,46,116,134,220,223,226,273,276,313,315]
para_init = [-0.3905,0.0223, 0.0066, -0.1598, 0.0242, -0.0527, 1.1040, 0.0960, -0.1287, 10.0715, -0.0049, -0.0076,
             0.2022, 0.0171, -0.1276, 1.5793, 0.0127, -0.0859, 1.4756, 0.0503, -0.0348, 12.1475, 0.0134, 0.0102,
             -0.0418, 0.0592, -0.0372, 11.3629, 0.0049, 0.0230, 10.6719, -0.0204, 0.0066, 14.9529, 0.0073, -0.0079,
             -0.9522, -0.0040, -0.0019, -2.6273, -0.0056, -0.0260, -0.6586, 0.0335, -0.0023, 18.8026, 0.0221, -0.0110]

problem = {
    'num_vars': len(para_index),
    'bounds': [[num * (1-step), num * (1+step)] for num in para_init]
}
initial_samples = latin_hypercube_sampling(problem['bounds'],pre_size)
initial_samples = np.vstack((np.array(para_init), initial_samples))
db_error = []
db_paraset = []
db_toprate = []
db_leng = []
result_list = ['XCH4','SH2','SC2H6','SC2H4','SC2H2','SC3H8','SC3H6','SC4+','SH','SCH','SCH2','SCH3','SC2H3','SC2H5','SC3H5','SC3H7','SC4H9']
exp_result = [[16.1,77.031,12.999,1.535,1.302,3.499,0.467,0.422,0,0,0,0,0,0,0,0,0],[18.873,74.645,14.318,1.739,1.556,3.968,0.555,0.508,0,0,0,0,0,0,0,0,0],[9.05848,58.23691,21.87466,2.52725,1.87428,8.83052,0.56969,1.45217,0,0,0,0,0,0,0,0,0]]

for i in range(len(initial_samples)):
    sample1 = []
    sample2 = []
    sample3 = []
    
    for j in range(int(len(initial_samples[i])/3)):
        para1 = initial_samples[i][3*j] + initial_samples[i][3*j+1] * 10.566 + initial_samples[i][3*j+2] * 2.8274
        para2 = initial_samples[i][3*j] + initial_samples[i][3*j+1] * 11.3 + initial_samples[i][3*j+2] * 4.7124
        para3 = initial_samples[i][3*j] + initial_samples[i][3*j+1] * 4.823 + initial_samples[i][3*j+2] * 3.7699
        sample1.append(para1)
        sample2.append(para2)
        sample3.append(para3)

    inp_path1 = f'kinet1_{i}.inp\n'
    exe_path1 = f'run_plasRxn1.exe'

    paraset1 = -21 * np.ones(316)
    paraset1[list(range(12))] = 0
    paraset1[list(range(90,102))] = 0
    for j in range(len(para_index)):
        paraset1[para_index[j]] = sample1[j]

    print(f'run {i}')

    update_kinet1(paraset1,i)
    prep_process = run_prep(inp_path1)
    prep_process.wait()
    compile_zdp1(exe_path1)
    exe_process = run_exe(exe_path1)
    exe_process.wait()
    error1,_,_,total_time1,_,_,_,_ = cal_error(exp_result[0])
    print(f'run {i} set1 completion')

    inp_path2 = f'kinet2_{i}.inp\n'
    exe_path2 = f'run_plasRxn2.exe'

    paraset2 = -21 * np.ones(316)
    paraset2[list(range(12))] = 0
    paraset2[list(range(90,102))] = 0
    for j in range(len(para_index)):
        paraset2[para_index[j]] = sample2[j]

    update_kinet2(paraset2,i)
    prep_process = run_prep(inp_path2)
    prep_process.wait()
    compile_zdp2(exe_path2)
    exe_process = run_exe(exe_path2)
    exe_process.wait()
    error2,_,_,total_time2,_,_,_,_ = cal_error(exp_result[1])
    print(f'run {i} set2 completion')

    inp_path3 = f'kinet3_{i}.inp\n'
    exe_path3 = f'run_plasRxn3.exe'

    paraset3 = -21 * np.ones(316)
    paraset3[list(range(12))] = 0
    paraset3[list(range(90,102))] = 0
    for j in range(len(para_index)):
        paraset3[para_index[j]] = sample3[j]

    update_kinet3(paraset3,i)
    prep_process = run_prep(inp_path3)
    prep_process.wait()
    compile_zdp3(exe_path3)
    exe_process = run_exe(exe_path3)
    exe_process.wait()
    error3,_,_,total_time3,_,_,_,_ = cal_error(exp_result[2])
    print(f'run {i} set3 completion')


    if (float(total_time1) > 2.1) & (float(total_time2) > 3.5) & (float(total_time3) > 5.6):
        db_error.append(error1+error2+error3)
        db_paraset.append(initial_samples[i])
    clear_output()

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

gp_model = GP_Model(pre_X,pre_Y)
gp_model.predict(pre_X)
bounds = [[num * (1-step), num * (1+step)] for num in para_init]
X_final, y_final = BO(gp_model, bounds, 1000,pre_X,pre_Y,exp_result)

np.savetxt('X_final.csv', X_final, delimiter=',', fmt='%f')
np.savetxt('y_final.csv', y_final, delimiter=',', fmt='%f')