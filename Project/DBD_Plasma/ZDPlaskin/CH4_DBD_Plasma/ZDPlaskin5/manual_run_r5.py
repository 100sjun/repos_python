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

def run_processor():
    process = subprocess.Popen(
        'preprocessor.exe',
        stdout = subprocess.DEVNULL,
        stderr = subprocess.DEVNULL,
        stdin = subprocess.PIPE,
        universal_newlines= True
    )

    process.stdin.write('kinet_r5.inp\n')
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
        # 프로세스 종료 대기
        return_code = process.wait()

    except:
        pass

run_processor()
run_compile(str(5))
run_exe('run_plasRxn_r5.exe')

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

reactions = []
with open('qt_reactions_list.txt','r') as file:
    for line in file:        
        line = line.strip()
        line = line[2:]
        reactions.append(line)
    file.close()

df_sp_5 = pd.read_csv('qt_densities.txt', sep=r'\s+', header=0, names=['Time [s]']+species)
df_rx_5 = pd.read_csv('qt_rates.txt', sep=r'\s+', header=0, names=['Time [s]']+reactions)
# mass concentration [g/mL]
t5 = df_sp_5['Time [s]']
CH4_5 = (df_sp_5['CH4'] + df_sp_5['CH4(V13)'] + df_sp_5['CH4(V24)'])/const.N_A
C2H6_5 = (df_sp_5['C2H6'] + df_sp_5['C2H6(V13)'] + df_sp_5['C2H6(V24)'])/const.N_A
C2H4_5 = (df_sp_5['C2H4'] + df_sp_5['C2H4(V1)'] + df_sp_5['C2H4(V2)'])/const.N_A
C2H2_5 = (df_sp_5['C2H2'] + df_sp_5['C2H2(V13)'] + df_sp_5['C2H2(V2)'] + df_sp_5['C2H2(V5)'])/const.N_A
C3H6_5 = (df_sp_5['C3H6'] + df_sp_5['C3H6(V)'])/const.N_A
C3H8_5 = (df_sp_5['C3H8'] + df_sp_5['C3H8(V1)'] + df_sp_5['C3H8(V2)'])/const.N_A
C4H10_5 = df_sp_5['C4H9H']/const.N_A
H2_5 = (df_sp_5['H2'])/const.N_A

conv_5 = (CH4_5.iloc[0] - CH4_5.iloc[-1]) / CH4_5.iloc[0] * 100
S_H2_5 = H2_5.iloc[-1] / (CH4_5.iloc[0] - CH4_5.iloc[-1]) / 2 * 100
S_C2H6_5 = 2*C2H6_5.iloc[-1] / (CH4_5.iloc[0] - CH4_5.iloc[-1]) * 100
S_C2H4_5 = 2*C2H4_5.iloc[-1] / (CH4_5.iloc[0] - CH4_5.iloc[-1]) * 100
S_C2H2_5 = 2*C2H2_5.iloc[-1] / (CH4_5.iloc[0] - CH4_5.iloc[-1]) * 100
S_C3H8_5 = 3*C3H8_5.iloc[-1] / (CH4_5.iloc[0] - CH4_5.iloc[-1]) * 100
S_C3H6_5 = 3*C3H6_5.iloc[-1] / (CH4_5.iloc[0] - CH4_5.iloc[-1]) * 100 
S_C4H10_5 = 4*C4H10_5.iloc[-1] / (CH4_5.iloc[0] - CH4_5.iloc[-1]) * 100

pd.set_option('display.float_format', '{:.3f}'.format)
result = pd.DataFrame({
    'species' : ['Conv','H2','C2H6','C2','C4','C3'],
    'r5_exp':[21.59794,36.81529,42.89809,14.61783,40.573,1.910828],
    'r5_cal':[conv_5,S_H2_5,S_C2H6_5,S_C2H4_5+S_C2H2_5,S_C4H10_5,S_C3H8_5+S_C3H6_5]
})

print(result)

plt.plot(t5,C2H4_5*1000000,label='C2H4')
plt.plot(t5,C2H2_5*1000000,label='C2H2')
plt.show()