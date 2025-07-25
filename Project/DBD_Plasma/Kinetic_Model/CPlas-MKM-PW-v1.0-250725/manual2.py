import subprocess  # subprocess 모듈 임포트 - 외부 프로세스 실행을 위한 모듈
import pandas as pd
#import matplotlib.pyplot as plt


# load kinetic data
df_kinetic = pd.read_csv('kinetic_data.csv')
df_mol = pd.read_csv('kinetic_data.csv').filter(like='mol%')
df_time = pd.read_csv('kinetic_data.csv').filter(like='Residence')

def run_preprocessor(input_path):  # 입력 파일 경로를 매개변수로 받는 함수 정의
    process = subprocess.Popen(  # 외부 프로세스 실행을 위한 Popen 객체 생성
        './preprocessor.exe',    # 실행할 프로그램 경로 지정
        stdin=subprocess.PIPE,   # 표준 입력 파이프 설정 - 프로세스에 입력을 전달하기 위함
        stdout=subprocess.PIPE,  # 표준 출력 파이프 설정 - 프로세스의 출력을 받기 위함
        stderr=subprocess.PIPE,  # 표준 에러 파이프 설정 - 프로세스의 에러를 받기 위함
        universal_newlines=True) # 텍스트 모드로 파이프 처리 - 문자열을 자동으로 인코딩/디코딩

    process.stdin.write(f'{input_path}\n')  # 입력 파일 경로를 프로세스에 전달
    process.stdin.flush()  # 버퍼를 비워서 데이터가 즉시 전송되도록 함

    process.stdin.write('.\n')  # 현재 디렉토리 표시를 프로세스에 전달
    process.stdin.flush()  # 버퍼를 비워서 데이터가 즉시 전송되도록 함

    output, error = process.communicate()  # 프로세스 실행 완료를 기다리고 출력과 에러를 받음
    
    if "ERROR" in output or "ERROR" in error:  # 출력이나 에러에 "ERROR" 문자열이 있는지 확인
        raise Exception("zdplaskin_m.F90 컴파일 실패")  # 에러가 있으면 예외 발생
    
    print('check the run of preprocessor')  # 전처리기 실행 확인 메시지 출력   
    return output, error  # 프로세스의 출력과 에러를 반환

def compile_zdp(name):  # ZDPlaskin 컴파일 함수 정의
    compile_command = [  # 컴파일 명령어 리스트 생성
        'gfortran', '-o', name, 'dvode_f90_m.F90', 'zdplaskin_m.F90', 'run_15kV.F90', 'bolsig_x86_64_g.dll'
    ]
    result = subprocess.run(compile_command, capture_output=True, text=True)  # 컴파일 명령 실행
    
    if result.returncode != 0:  # 컴파일 결과 확인
        raise Exception(f"{name} 컴파일 실패")  # 컴파일 실패시 예외 발생
    print('check the compile_zdp')  # 컴파일 완료 메시지 출력

def run_exe(exe_path):  # ZDPlaskin 실행 함수 정의
    try:
        process = subprocess.Popen(  # 실행 파일 실행
           exe_path,    
            stdout=subprocess.PIPE,  
            stderr=subprocess.STDOUT,  
            stdin=subprocess.PIPE,   
            universal_newlines=True,
            bufsize=1,              
        )
        
        while True:
            output = process.stdout.readline()
            
            if not output:
                break
                
            print(f'\r{output.strip()}                                                                              ',end='',flush=True)
            
            if "PRESS ENTER TO EXIT" in output:
                print(end='\n')
                process.kill()
                break

            if "WARNING: BOLSIG+ convergence failed" in output:
                process.stdin.write('\n')
                process.stdin.flush()

    except:
        pass
    return process

kinet_path = 'kinet_ori_treat.inp'  # 입력 파일 경로 설정
run_preprocessor(kinet_path)  # 전처리기 실행
exe_name = 'run_15kV.exe'  # 실행 파일 이름 설정
compile_zdp(exe_name)  # ZDPlaskin 컴파일 실행
run_exe(exe_name)

exp_list = ['H2', 'CH4', 'C2H6', 'C2H4', 'C2H2', 'C3H8', 'C3H6', 'C4H10', 'C5H12', 'C']
exp = [df_mol.iloc[0].tolist(),
        df_mol.iloc[1].tolist(),
        df_mol.iloc[2].tolist(),
        df_mol.iloc[3].tolist()]


species = []
with open('qt_species_list.txt','r') as f:
    for line in f:
        comp = line[2:]
        species.append(comp.strip())
    f.close()
df_sp = pd.read_csv('qt_densities.txt', sep=r'\s+', header=0, names=['Time [s]']+species)
H2 = (df_sp['H2'])
CH4 = (df_sp['CH4'] + df_sp['CH4(V13)'] + df_sp['CH4(V24)'])
C2H2 = (df_sp['C2H2'] + df_sp['C2H2(V2)'] + df_sp['C2H2(V5)'] + df_sp['C2H2(V13)'])
C2H4 = (df_sp['C2H4'] + df_sp['C2H4(V1)'] + df_sp['C2H4(V2)'])
C2H6 = (df_sp['C2H6'] + df_sp['C2H6(V13)'] + df_sp['C2H6(V24)'])
C3H6 = (df_sp['C3H6'] + df_sp['C3H6(V)'])
C3H8 = (df_sp['C3H8'] + df_sp['C3H8(V1)'] + df_sp['C3H8(V2)'])
C4H10 = (df_sp['C4H9H'])
C5H12 = (df_sp['C5H12'])
C = (df_sp['C'])
CH = (df_sp['CH'])
CH2 = (df_sp['CH2'])
CH3 = (df_sp['CH3'])
C2H3 = (df_sp['C2H3'])
C2H5 = (df_sp['C2H5'])
C3H7 = (df_sp['C3H7'])
H = (df_sp['H'])
CH3_plus = (df_sp['CH3^+'])
CH4_plus = df_sp['CH4^+']
CH5_plus = df_sp['CH5^+']
C2H2_plus = df_sp['C2H2^+']
C2H4_plus = df_sp['C2H4^+']
C2H5_plus = df_sp['C2H5^+']
C2H6_plus = df_sp['C2H6^+']
C3H6_plus = df_sp['C3H6^+']
C3H8_plus = df_sp['C3H8^+']

t0 = abs(df_sp['Time [s]']-5.654867).argmin()
t1 = abs(df_sp['Time [s]']-7.270543).argmin()
t2 = abs(df_sp['Time [s]']-10.17876).argmin()

sim = []
for t in [t0, t1, t2]:
    sim_H2 = float(format(H2.iloc[t], '.6f'))
    sim_CH4 = float(format(CH4.iloc[t], '.6f'))
    sim_C2H2 = float(format(C2H2.iloc[t], '.6f'))
    sim_C2H4 = float(format(C2H4.iloc[t], '.6f'))
    sim_C2H6 = float(format(C2H6.iloc[t], '.6f'))
    sim_C3H6 = float(format(C3H6.iloc[t], '.6f'))
    sim_C3H8 = float(format(C3H8.iloc[t], '.6f'))
    sim_C4H10 = float(format(C4H10.iloc[t], '.6f'))
    sim_C5H12 = float(format(C5H12.iloc[t], '.6f'))
    sim_C = float(format(C.iloc[t], '.6f'))
    sim_CH = float(format(CH.iloc[t], '.6f'))
    sim_CH2 = float(format(CH2.iloc[t], '.6f'))
    sim_CH3 = float(format(CH3.iloc[t], '.6f'))
    sim_C2H3 = float(format(C2H3.iloc[t], '.6f'))
    sim_C2H5 = float(format(C2H5.iloc[t], '.6f'))
    sim_C3H7 = float(format(C3H7.iloc[t], '.6f'))
    sim_H = float(format(H.iloc[t], '.6f'))
    sim_CH3_plus = float(format(CH3_plus.iloc[t], '.6f'))
    sim_CH4_plus = float(format(CH4_plus.iloc[t], '.6f'))
    sim_CH5_plus = float(format(CH5_plus.iloc[t], '.6f'))
    sim_C2H2_plus = float(format(C2H2_plus.iloc[t], '.6f'))
    sim_C2H4_plus = float(format(C2H4_plus.iloc[t], '.6f'))
    sim_C2H5_plus = float(format(C2H5_plus.iloc[t], '.6f'))
    sim_C2H6_plus = float(format(C2H6_plus.iloc[t], '.6f'))
    sim_C3H6_plus = float(format(C3H6_plus.iloc[t], '.6f'))
    sim_C3H8_plus = float(format(C3H8_plus.iloc[t], '.6f'))

    out_H2 = sim_H2 + 0.5*sim_CH - sim_CH2 - 0.5*sim_CH3 - 0.5*sim_CH3_plus + 0.5*sim_CH5_plus - 0.5*sim_C2H3 + 0.5*sim_C2H5 + 0.5*sim_C2H5_plus + 0.5*sim_C3H7 + 0.5*sim_H
    out_CH4 = sim_CH4 + sim_CH2 + sim_CH3 + sim_CH3_plus + sim_CH4_plus + sim_CH5_plus
    out_C2H6 = sim_C2H6 + sim_C2H6_plus
    out_C2H4 = sim_C2H4 + sim_C2H3 + sim_C2H5 + sim_C2H4_plus + sim_C2H5_plus
    out_C2H2 = sim_C2H2 + sim_C2H2_plus
    out_C3H8 = sim_C3H8 + sim_C3H8_plus
    out_C3H6 = sim_C3H6 + sim_C3H7 + sim_C3H6_plus
    out_C4H10 = sim_C4H10
    out_C5H12 = sim_C5H12
    out_C = sim_C + sim_CH
    
    out_total = out_H2 + out_CH4 + out_C2H6 + out_C2H4 + out_C2H2 + out_C3H8 + out_C3H6 + out_C4H10 + out_C5H12 + out_C
    newsim = [out_H2/out_total*100, out_CH4/out_total*100, out_C2H6/out_total*100, out_C2H4/out_total*100, out_C2H2/out_total*100, out_C3H8/out_total*100, out_C3H6/out_total*100, out_C4H10/out_total*100, out_C5H12/out_total*100, out_C/out_total*100]
    sim.append(newsim)


# 결과 비교를 위한 함수 정의
def compare_results(exp_data, sim_data, residence_time):
    result = pd.DataFrame({
        'Species': exp_list,
        'Experimental (mol%)': exp_data,
        'Simulation (mol%)': sim_data,
        'Absolute Error (mol%)': [round(s - e, 4) for e, s in zip(exp_data, sim_data)],
        'Relative Error (%)': [round((s - e)/e * 100, 2) if e != 0 else float('inf') for e, s in zip(exp_data, sim_data)]
    })
    
    print(f'\nResults Comparison at Residence Time {residence_time:.2f}s')
    print('=' * 80)
    
    # Set display format for dataframe
    pd.set_option('display.float_format', '{:.4f}'.format)
    print(result.to_string(index=False))
    
    # Calculate MSE (excluding C5H12)
    mse = 0
    for i in range(len(exp_data)-2):  # Excluding C5H12 and C
        if exp_data[i] != 0:
            mse += ((exp_data[i] - sim_data[i])/exp_data[i])**2
    
    print(f'\nMSE (excluding C5H12, C): {mse:.6f}')
    return mse

t0 = abs(df_sp['Time [s]']-5.654867).argmin()
t1 = abs(df_sp['Time [s]']-7.270543).argmin()
t2 = abs(df_sp['Time [s]']-10.17876).argmin()


# 각 체류시간별 결과 비교
residence_times = [5.654867, 7.270543, 10.17876]
for i, (e, s, rt) in enumerate(zip(exp, sim, residence_times)):
    compare_results(e, s, rt)

# def plot_comparison(exp, sim, residence_times, exp_list):
#     # Set figure size and layout
#     fig = plt.figure(figsize=(15, 10))
#     fig.suptitle('Comparison of Experimental and Simulation Results vs. Residence Time', fontsize=16)
    
#     # Create 2x5 subplots
#     for i, species in enumerate(exp_list):
#         ax = plt.subplot(2, 5, i+1)
        
#         # Extract experimental and simulation values
#         exp_values = [exp_data[i] for exp_data in exp]
#         sim_values = [sim_data[i] for sim_data in sim]
        
#         # Reverse the order of residence times and corresponding values
#         reversed_times = residence_times[::-1]
#         reversed_exp = exp_values[::-1]
#         reversed_sim = sim_values[::-1]
        
#         # Plot lines
#         ax.plot(reversed_times, reversed_exp, 'o-', label='Experimental', color='blue')
#         ax.plot(reversed_times, reversed_sim, 's--', label='Simulation', color='red')
        
#         # Set graph properties
#         ax.set_title(f'{species}', fontsize=12)
#         ax.set_xlabel('Residence Time (s)')
#         ax.set_ylabel('Mole Fraction (%)')
#         ax.grid(True, linestyle='--', alpha=0.7)
#         ax.legend()
        
#         # Set x-axis limits (low to high)
#         ax.set_xlim(min(residence_times) - 1, max(residence_times) + 1)
    
#     # Adjust layout
#     plt.tight_layout()
#     plt.show()

# # 모든 데이터에 대한 그래프 생성
# plot_comparison(exp, sim, residence_times, exp_list)

# # 그래프 데이터를 CSV 파일로 저장
# def save_results_to_csv(exp, sim, residence_times, exp_list):
#     # 데이터 준비
#     data = {}
    
#     # 체류 시간 추가
#     data['Residence_Time'] = residence_times
    
#     # 실험 데이터 추가
#     for i, species in enumerate(exp_list):
#         exp_values = [exp_data[i] for exp_data in exp]
#         data[f'{species}_Experimental'] = exp_values
    
#     # 시뮬레이션 데이터 추가
#     for i, species in enumerate(exp_list):
#         sim_values = [sim_data[i] for sim_data in sim]
#         data[f'{species}_Simulation'] = sim_values
    
#     # 데이터프레임 생성 및 CSV 저장
#     df_results = pd.DataFrame(data)
#     csv_filename = 'comparison_results.csv'
#     df_results.to_csv(csv_filename, index=False)
    
#     print(f'그래프 데이터가 {csv_filename}에 저장되었습니다.')

# # 그래프 데이터 CSV 저장 함수 호출
# save_results_to_csv(exp, sim, residence_times, exp_list)


rx = []
with open('qt_reactions_list.txt','r') as f:
    for line in f:
        comp = line[3:]
        rx.append(comp.strip())
    f.close()

# Create column names using reaction numbers (001, 002, etc.)
rx_numbers = [f"{i:03d}" for i in range(1, len(rx)+1)]

df_rx = pd.read_csv('qt_rates.txt', sep=r'\s+', header=0, names=['Time [s]']+rx_numbers)

reaction_extents = {}
reaction_extents['Time [s]'] = df_rx['Time [s]']
for i in range(len(rx)):
    reaction = []
    for j in range(len(df_rx)):
        if j == 0:
            reaction.append(0)
        else:
            reaction.append(reaction[j-1] +df_rx[rx_numbers[i]].iloc[j-1]*(df_rx['Time [s]'].iloc[j]-df_rx['Time [s]'].iloc[j-1]))
    reaction_extents[rx[i]] = reaction

df_rex = pd.DataFrame(reaction_extents)
rx_end = df_rex.iloc[-1].sort_values(ascending=False)

print(rx_end.head(20))

print(df_sp.iloc[-1].sort_values(ascending=False).head(20))


