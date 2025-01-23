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
from scipy.optimize import minimize

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
    sim = []
    for i in t:
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
            if j < 8:
                err += ((exp[i][j] - sim[i][j])/exp[i][j])**2
            else:
                err += (exp[i][j] - sim[i][j])**2
    
    return err, sim

# experimental results
result_list = ['XCH4','SH2','SC2H6','SC2H4','SC2H2','SC3H8','SC3H6','SC4','SH','SCH','SCH2','SCH3','SC2H3','SC2H5','SC3H7','SCH3^+','SC2H5^+']
exp_result = [[7.37,60.84,21.1,3.07,2.16,7.83,0.67,1.23,0,0,0,0,0,0,0,0,0],
              [7.79,53.13,24.77,3.24,2.38,9.55,0.77,1.59,0,0,0,0,0,0,0,0,0],
              [8.73,62.08,19.51,1.59,1.38,8.64,0.56,1.84,0,0,0,0,0,0,0,0,0],
              [21.53,43.31,31.17,3.42,2.67,12.71,0.94,2.31,0,0,0,0,0,0,0,0,0]]

# run a zdplaskin
inp_path = f'kinet.inp'
exe_path = 'run_plasRxn_v2.exe'
prep_process = run_prep(inp_path)
prep_process.wait()
compile_zdp(exe_path)
exe_process = run_exe(exe_path)
exe_process.wait()
error, sim = cal_error(exp_result)

# show the result dataframe
df_res = pd.DataFrame({
    'list': result_list,
    'exp_t1': exp_result[0],
    'sim_t1': [round(num,2) for num in sim[0]],
    'exp_t2': exp_result[1],
    'sim_t2': [round(num,2) for num in sim[1]],
    'exp_t3': exp_result[2],
    'sim_t3': [round(num,2) for num in sim[2]],
    'exp_t4': exp_result[3],
    'sim_t4': [round(num,2) for num in sim[3]]
})
print(df_res)
print(f'error: {error}')