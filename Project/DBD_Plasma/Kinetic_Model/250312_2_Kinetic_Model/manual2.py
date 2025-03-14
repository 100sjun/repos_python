import subprocess  # subprocess 모듈 임포트 - 외부 프로세스 실행을 위한 모듈
import pandas as pd


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
        'gfortran', '-o', name, 'dvode_f90_m.F90', 'zdplaskin_m.F90', 'run.F90', 'bolsig_x86_64_g.dll'
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

kinet_path = 'kinet_opt.inp'  # 입력 파일 경로 설정
run_preprocessor(kinet_path)  # 전처리기 실행
exe_name = 'run.exe'  # 실행 파일 이름 설정
compile_zdp(exe_name)  # ZDPlaskin 컴파일 실행
run_exe(exe_name)

exp_list = ['H2', 'CH4', 'C2H6', 'C2H4', 'C2H2', 'C3H8', 'C3H6', 'C4H10', 'C5H12', 'C']
exp = [[3.208145, 95.21883, 0.996818, 0.130302, 0.0958, 0.256151, 0.020783, 0.031923, 0.041253, 0],
       [5.319082, 90.80403, 2.428802, 0.197735, 0.171795, 0.717088, 0.046734, 0.114829, 0.119573, 0]]


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

all_sp = df_sp.sum(axis=1) - df_sp['E']

t1 = abs(df_sp['Time [s]']-7.27).argmin()
t2 = abs(df_sp['Time [s]']-16.96).argmin()

sim = []  
for t in [t1, t2]:
    sim_H2 = float(format(H2.iloc[t]/all_sp.iloc[t]*100, '.6f'))
    sim_CH4 = float(format(CH4.iloc[t]/all_sp.iloc[t]*100, '.6f'))
    sim_C2H2 = float(format(C2H2.iloc[t]/all_sp.iloc[t]*100, '.6f'))
    sim_C2H4 = float(format(C2H4.iloc[t]/all_sp.iloc[t]*100, '.6f'))
    sim_C2H6 = float(format(C2H6.iloc[t]/all_sp.iloc[t]*100, '.6f'))
    sim_C3H6 = float(format(C3H6.iloc[t]/all_sp.iloc[t]*100, '.6f'))
    sim_C3H8 = float(format(C3H8.iloc[t]/all_sp.iloc[t]*100, '.6f'))
    sim_C4H10 = float(format(C4H10.iloc[t]/all_sp.iloc[t]*100, '.6f'))
    sim_C5H12 = float(format(C5H12.iloc[t]/all_sp.iloc[t]*100, '.6f'))
    sim_C = float(format(C.iloc[t]/all_sp.iloc[t]*100, '.6f'))
    sim_CH = float(format(CH.iloc[t]/all_sp.iloc[t]*100, '.6f'))
    sim_CH2 = float(format(CH2.iloc[t]/all_sp.iloc[t]*100, '.6f'))
    sim_CH3 = float(format(CH3.iloc[t]/all_sp.iloc[t]*100, '.6f'))
    sim_C2H3 = float(format(C2H3.iloc[t]/all_sp.iloc[t]*100, '.6f'))
    sim_C2H5 = float(format(C2H5.iloc[t]/all_sp.iloc[t]*100, '.6f'))
    sim_C3H7 = float(format(C3H7.iloc[t]/all_sp.iloc[t]*100, '.6f'))
    sim_H = float(format(H.iloc[t]/all_sp.iloc[t]*100, '.6f'))
    sim_CH3_plus = float(format(CH3_plus.iloc[t]/all_sp.iloc[t]*100, '.6f'))
    sim_CH4_plus = float(format(CH4_plus.iloc[t]/all_sp.iloc[t]*100, '.6f'))
    sim_CH5_plus = float(format(CH5_plus.iloc[t]/all_sp.iloc[t]*100, '.6f'))
    sim_C2H2_plus = float(format(C2H2_plus.iloc[t]/all_sp.iloc[t]*100, '.6f'))
    sim_C2H4_plus = float(format(C2H4_plus.iloc[t]/all_sp.iloc[t]*100, '.6f'))
    sim_C2H5_plus = float(format(C2H5_plus.iloc[t]/all_sp.iloc[t]*100, '.6f'))
    sim_C2H6_plus = float(format(C2H6_plus.iloc[t]/all_sp.iloc[t]*100, '.6f'))
    sim_C3H6_plus = float(format(C3H6_plus.iloc[t]/all_sp.iloc[t]*100, '.6f'))
    sim_C3H8_plus = float(format(C3H8_plus.iloc[t]/all_sp.iloc[t]*100, '.6f'))
    newsim = [sim_H2 - 0.5*sim_CH3 - sim_CH2 - 0.5*sim_C2H3 + 0.5*sim_C2H5 + 0.5*sim_H + 0.5*sim_C3H7, sim_CH4 + sim_CH3 + sim_CH2, sim_C2H6, sim_C2H4 + sim_C2H3 + sim_C2H5, sim_C2H2, sim_C3H8, sim_C3H6 + sim_C3H7, sim_C4H10, sim_C5H12, sim_C]
    sim.append(newsim)


result1 = pd.DataFrame({
    'species': exp_list,
    'exp': exp[0],
    'sim': sim[0]   
})

print(result1)

err1 = 0
for i in range(len(exp[0])):
    if i < 8:
        err1 += ((exp[0][i] - sim[0][i])/exp[1][i])**2
    else:
        err1 += 0
err1 -= (exp[0][8] - sim[0][8])**2
print(err1)

result2 = pd.DataFrame({
    'species': exp_list,
    'exp': exp[1],
    'sim': sim[1]   
})

print(result2)

err2 = 0
for i in range(len(exp[1])):
    if i < 8:
        err2 += ((exp[1][i] - sim[1][i])/exp[1][i])**2
    else:
        err2 += 0
err2 -= (exp[1][8] - sim[1][8])**2

print(err2)
