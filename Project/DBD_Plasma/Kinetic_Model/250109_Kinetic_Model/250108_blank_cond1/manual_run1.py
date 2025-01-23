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
    # sim_SH2 = H2.iloc[-1] / (H2.iloc[-1] + C2H6.iloc[-1] + C2H4.iloc[-1] + C2H2.iloc[-1] + C3H8.iloc[-1] + C3H6.iloc[-1] + C4H10.iloc[-1] + C5H12.iloc[-1])*100
    # sim_SC2H6 = C2H6.iloc[-1] / (H2.iloc[-1] + C2H6.iloc[-1] + C2H4.iloc[-1] + C2H2.iloc[-1] + C3H8.iloc[-1] + C3H6.iloc[-1] + C4H10.iloc[-1] + C5H12.iloc[-1])*100
    # sim_SC2H4 = C2H4.iloc[-1] / (H2.iloc[-1] + C2H6.iloc[-1] + C2H4.iloc[-1] + C2H2.iloc[-1] + C3H8.iloc[-1] + C3H6.iloc[-1] + C4H10.iloc[-1] + C5H12.iloc[-1])*100
    # sim_SC2H2 = C2H2.iloc[-1] / (H2.iloc[-1] + C2H6.iloc[-1] + C2H4.iloc[-1] + C2H2.iloc[-1] + C3H8.iloc[-1] + C3H6.iloc[-1] + C4H10.iloc[-1] + C5H12.iloc[-1])*100
    # sim_SC3H8 = C3H8.iloc[-1] / (H2.iloc[-1] + C2H6.iloc[-1] + C2H4.iloc[-1] + C2H2.iloc[-1] + C3H8.iloc[-1] + C3H6.iloc[-1] + C4H10.iloc[-1] + C5H12.iloc[-1])*100
    # sim_SC3H6 = C3H6.iloc[-1] / (H2.iloc[-1] + C2H6.iloc[-1] + C2H4.iloc[-1] + C2H2.iloc[-1] + C3H8.iloc[-1] + C3H6.iloc[-1] + C4H10.iloc[-1] + C5H12.iloc[-1])*100
    # sim_SC4H10 = C4H10.iloc[-1] / (H2.iloc[-1] + C2H6.iloc[-1] + C2H4.iloc[-1] + C2H2.iloc[-1] + C3H8.iloc[-1] + C3H6.iloc[-1] + C4H10.iloc[-1] + C5H12.iloc[-1])*100
    # sim_SC5H12 = C5H12.iloc[-1] / (H2.iloc[-1] + C2H6.iloc[-1] + C2H4.iloc[-1] + C2H2.iloc[-1] + C3H8.iloc[-1] + C3H6.iloc[-1] + C4H10.iloc[-1] + C5H12.iloc[-1])*100
    # sim_SH = H.iloc[-1] / (H2.iloc[-1] + C2H6.iloc[-1] + C2H4.iloc[-1] + C2H2.iloc[-1] + C3H8.iloc[-1] + C3H6.iloc[-1] + C4H10.iloc[-1] + C5H12.iloc[-1])*100
    # sim_SCH = CH.iloc[-1] / (H2.iloc[-1] + C2H6.iloc[-1] + C2H4.iloc[-1] + C2H2.iloc[-1] + C3H8.iloc[-1] + C3H6.iloc[-1] + C4H10.iloc[-1] + C5H12.iloc[-1])*100
    # sim_SCH2 = CH2.iloc[-1] / (H2.iloc[-1] + C2H6.iloc[-1] + C2H4.iloc[-1] + C2H2.iloc[-1] + C3H8.iloc[-1] + C3H6.iloc[-1] + C4H10.iloc[-1] + C5H12.iloc[-1])*100
    # sim_SCH3 = CH3.iloc[-1] / (H2.iloc[-1] + C2H6.iloc[-1] + C2H4.iloc[-1] + C2H2.iloc[-1] + C3H8.iloc[-1] + C3H6.iloc[-1] + C4H10.iloc[-1] + C5H12.iloc[-1])*100
    # sim_SC2H3 = C2H3.iloc[-1] / (H2.iloc[-1] + C2H6.iloc[-1] + C2H4.iloc[-1] + C2H2.iloc[-1] + C3H8.iloc[-1] + C3H6.iloc[-1] + C4H10.iloc[-1] + C5H12.iloc[-1])*100
    # sim_SC2H5 = C2H5.iloc[-1] / (H2.iloc[-1] + C2H6.iloc[-1] + C2H4.iloc[-1] + C2H2.iloc[-1] + C3H8.iloc[-1] + C3H6.iloc[-1] + C4H10.iloc[-1] + C5H12.iloc[-1])*100
    # sim_SC3H5 = C3H5.iloc[-1] / (H2.iloc[-1] + C2H6.iloc[-1] + C2H4.iloc[-1] + C2H2.iloc[-1] + C3H8.iloc[-1] + C3H6.iloc[-1] + C4H10.iloc[-1] + C5H12.iloc[-1])*100
    # sim_SC3H7 = C3H7.iloc[-1] / (H2.iloc[-1] + C2H6.iloc[-1] + C2H4.iloc[-1] + C2H2.iloc[-1] + C3H8.iloc[-1] + C3H6.iloc[-1] + C4H10.iloc[-1] + C5H12.iloc[-1])*100
    # sim_SC4H9 = C4H9.iloc[-1] / (H2.iloc[-1] + C2H6.iloc[-1] + C2H4.iloc[-1] + C2H2.iloc[-1] + C3H8.iloc[-1] + C3H6.iloc[-1] + C4H10.iloc[-1] + C5H12.iloc[-1])*100


    sim = []
    sim.append(sim_XCH4)
    sim.append(sim_SH2)
    sim.append(sim_SC2H6)
    sim.append(sim_SC2H4)
    sim.append(sim_SC2H2)
    sim.append(sim_SC3H8)
    sim.append(sim_SC3H6)
    sim.append(sim_SC4H10)
    sim.append(sim_SC5H12)
    sim.append(sim_SH)
    sim.append(sim_SCH)
    sim.append(sim_SCH2)
    sim.append(sim_SCH3)
    sim.append(sim_SC2H3)
    sim.append(sim_SC2H5)
    sim.append(sim_SC3H5)
    sim.append(sim_SC3H7)
    sim.append(sim_SC4H9)

    sim2 = []
    sim2.append(sim_XCH4)
    sim2.append(sim_SH2)
    sim2.append(sim_SC2H6)
    sim2.append(sim_SC2H4)
    sim2.append(sim_SC2H2)
    sim2.append(sim_SC3H8)
    sim2.append(sim_SC3H6)
    sim2.append(sim_SC4H10)
    sim2.append(sim_SH)
    sim2.append(sim_SCH)
    sim2.append(sim_SCH2)
    sim2.append(sim_SCH3)
    sim2.append(sim_SC2H3)
    sim2.append(sim_SC2H5)
    sim2.append(sim_SC3H5)
    sim2.append(sim_SC3H7)
    sim2.append(sim_SC4H9)

    

    err =  float(np.sum((np.asarray(exp) - np.asarray(sim2))**2))
    
    return err, top_rate, sim, df_cd['Time [s]'].iloc[-1],df_sp.iloc[-1].sort_values(ascending=False)

result_list = ['XCH4','SH2','SC2H6','SC2H4','SC2H2','SC3H8','SC3H6','SC4','SC5+','SH','SCH','SCH2','SCH3','SC2H3','SC2H5','SC3H5','SC3H7','SC4H9']
exp_result = [16.1,77.031,12.999,1.535,1.302,3.499,0.467,0.422,0,0,0,0,0,0,0,0,0]
exp_result_expr = [16.1,77.031,12.999,1.535,1.302,3.499,0.467,0.422,2.741,0,0,0,0,0,0,0,0,0]

inp_path = f'kinet_blank_cond1_1.inp'
exe_path = 'run_plasRxn_v2.exe'
prep_process = run_prep(inp_path)
prep_process.wait()
compile_zdp(exe_path)
exe_process = run_exe(exe_path)
exe_process.wait()
error, top_rate, sim_result, total_time, total_species = cal_error(exp_result)

sim_result2 = [round(float(i),3) for i in sim_result]

df_res = pd.DataFrame({
    'res': result_list,
    'exp': exp_result_expr,
    'sim': sim_result2
    }
)
print(df_res)
print(total_time)
print(error)

