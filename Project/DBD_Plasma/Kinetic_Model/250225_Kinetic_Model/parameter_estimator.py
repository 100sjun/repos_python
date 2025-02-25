import re
import pandas as pd
import subprocess
import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C
from scipy.stats import norm

class ParameterOptmizer:
    def __init__(self):
        # define the experimental data
        self.exp_list = ['H2', 'CH4', 'C2H6', 'C2H4', 'C2H2', 'C3H8', 'C3H6', 'C4H10', 'C5H12', 
                        'C', 'CH', 'CH2', 'CH3', 'C2H3', 'C2H5', 'C3H7', 'H',
                        'CH3^+', 'CH4^+', 'CH5^+', 'C2H2^+', 'C2H4^+', 'C2H5^+', 'C2H6^+', 'C3H6^+', 'C3H8^+']
        self.exp_values = [[3.208145, 95.21883, 0.996818, 0.130302, 0.0958, 0.256151, 0.020783, 
                           0.031923, 0.041253, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                            [5.319082, 90.80403, 2.428802, 0.197735, 0.171795, 0.717088, 0.046734, 
                           0.114829, 0.119573, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
        
        # define the file path
        self.kinet_path = 'kinet.inp'
        self.exe_path = 'run2.exe'

        # dictionary to save the current parameters
        self.current_parameters = {}
        self.load_current_parameters(init=True)

        # define manipulated variables
        self.max_error = 100
        self.lr = 0.01 # bounds learning rate
        self.pre_iter = 20 # pre-train iteration
        self.iter = 100000 # BO iteration
        self.xi = 0.05 # EI parameters (0.05 is a balance between exploration and exploitation)
    
    def load_current_parameters(self, init=True):
        # 시작인 경우 kinet_ori.inp에서 parameter를 불러옴
        if init:
            with open('kinet_ori.inp', 'r') as f:
                content = f.readlines()
        else:
            # kinet.inp 파일에서 현재 파라미터 값들을 읽어옴
            with open(self.kinet_path, 'r') as f:
                content = f.readlines()

        for line in content:
            if "$ double precision, parameter :: f" in line:
                # f parameter와 값 추출
                pattern = r'f(\d+)\s*=\s*([\d\.\+\-d]+)'
                match = re.search(pattern, line)
                if match:
                    param_name = 'f' + match.group(1)
                    value = float(match.group(2).replace('d','e'))
                    self.current_parameters[param_name] = value

    def modify_parameter(self, param_name: str, new_value: float) -> None:
        # modify the parameter value in the kinet.inp file
        with open (self.kinet_path, 'r') as f:
            content = f.read()
        
        # modify the parameter
        pattern = rf'(parameter :: {param_name} = )([\d\.\+\-d]+)'
        new_value_str = f'{new_value:.4e}'.replace('e', 'd')
        content = re.sub(pattern, f'parameter :: {param_name} = {new_value_str}', content)

        with open(self.kinet_path, 'w') as f:
            f.write(content)
    
    def run_preprocessor(self):
        process = subprocess.Popen(  # 외부 프로세스 실행을 위한 Popen 객체 생성
        './preprocessor.exe',    # 실행할 프로그램 경로 지정
        stdin=subprocess.PIPE,   # 표준 입력 파이프 설정 - 프로세스에 입력을 전달하기 위함
        stdout=subprocess.PIPE,  # 표준 출력 파이프 설정 - 프로세스의 출력을 받기 위함
        stderr=subprocess.PIPE,  # 표준 에러 파이프 설정 - 프로세스의 에러를 받기 위함
        universal_newlines=False) # 텍스트 모드로 파이프 처리 - 문자열을 자동으로 인코딩/디코딩

        process.stdin.write(f'{self.kinet_path}\n'.encode())  # 입력 파일 경로를 프로세스에 전달
        process.stdin.flush()  # 버퍼를 비워서 데이터가 즉시 전송되도록 함

        process.stdin.write('.\n'.encode())  # 현재 디렉토리 표시를 프로세스에 전달
        process.stdin.flush()  # 버퍼를 비워서 데이터가 즉시 전송되도록 함

        output, error = process.communicate()  # 프로세스 실행 완료를 기다리고 출력과 에러를 받음
        output_str = output.decode('utf-8', errors='ignore') if output else ''
        error_str = error.decode('utf-8', errors='ignore') if error else ''
    
        print('check the run of preprocessor')  # 전처리기 실행 확인 메시지 출력 
        return output_str, error_str
    
    def compile_zdp(self):
        compile_command = [  # 컴파일 명령어 리스트 생성
        'gfortran', '-o', self.exe_path, 'dvode_f90_m.F90', 'zdplaskin_m.F90', 'run2.F90', 'bolsig_x86_64_g.dll'
        ]
        result = subprocess.run(compile_command, capture_output=True, text=True)  # 컴파일 명령 실행
    
        if result.returncode != 0:  # 컴파일 결과 확인
            raise Exception(f"{self.exe_path} 컴파일 실패")  # 컴파일 실패시 예외 발생
        print('check the compile_zdp')  # 컴파일 완료 메시지 출력

    def run_simulation(self):
        try:
            process = subprocess.Popen(  # 실행 파일 실행
                self.exe_path,    
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
                    print()
                    process.kill()
                    break

                if "WARNING: BOLSIG+ convergence failed" in output:
                    process.stdin.write('\n')
                    process.stdin.flush()

        except:
            pass
        return process
    
    def err_calculation(self):
        # read the result file
        species = []
        with open('qt_species_list.txt', 'r') as f:
            for line in f:
                comp = line[2:]
                species.append(comp.strip())
        
        df_sp = pd.read_csv('qt_densities.txt', sep=r'\s+', header=0, names=['Time [s]']+species)
        
        # calculate the concentration
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
            newsim = [sim_H2, sim_CH4, sim_C2H6, sim_C2H4, sim_C2H2, sim_C3H8, sim_C3H6, sim_C4H10, sim_C5H12,
                      sim_C, sim_CH, sim_CH2, sim_CH3, sim_C2H3, sim_C2H5, sim_C3H7, sim_H, sim_CH3_plus, sim_CH4_plus,
                      sim_CH5_plus, sim_C2H2_plus, sim_C2H4_plus, sim_C2H5_plus, sim_C2H6_plus, sim_C3H6_plus, sim_C3H8_plus]
            sim.append(newsim)
        
        err = 0
        for i in range(len(self.exp_values)):
            for j in range(len(self.exp_values[i])):
                if j < 8:
                    err += ((self.exp_values[i][j] - sim[i][j])/self.exp_values[i][j])**2
                else:
                    err += ((self.exp_values[i][j] - sim[i][j]))**2
            err -= (self.exp_values[i][8] - sim[i][8])**2
            err -= (self.exp_values[i][9] - sim[i][9])**2

        return err, df_sp['Time [s]'].iloc[-1]  

    def LH_sampling(self):
        # load a current parameter set
        paraset = self.current_parameters

        # define the bounds
        bounds = np.array([[-self.lr,self.lr]] * len(paraset))

        # sampling
        samples = np.random.uniform(bounds[:,0], bounds[:,1], (self.pre_iter, len(paraset)))

        return samples
    
    def GP_Model(self,X,y):
        # Gaussian Process 모델 생성
        kernel = C(6.3442, (1e-2, 1e3)) * RBF(length_scale=1.0, length_scale_bounds=(1e-2, 1e3))
        gp_model = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=10)
        gp_model.fit(X,y)
        return gp_model
    
    def Bayesian_Opt(self,pre_X,y,ori_X):
        bounds = np.array([[-self.lr,self.lr]] * len(ori_X))
        X = pre_X/ori_X - 1
        gp_model = self.GP_Model(X,y)

        for i in range(self.iter):
            y_best = np.min(y)
            candidate_points = np.random.uniform(
            low = np.array(bounds)[:,0],
            high = np.array(bounds)[:,1],
            size = (len(ori_X)*100,len(bounds))
            )
            mu, sigma = gp_model.predict(candidate_points,return_std=True)
            improvement = y_best - mu - 0.05 # 0.05는 탐색과 착쥐의 균형을 조절하는 값

            Z = improvement / (sigma + 1e-9) # 0으로 나누는 것을 방지하기 위해 1e-9를 더함
            ei = improvement * norm.cdf(Z) + sigma * norm.pdf(Z)

            best_index = np.argmax(ei)
            x_next = candidate_points[best_index]

            para_next = (10**x_next)*ori_X
            for j in range(len(ori_X)):
                self.modify_parameter(f'f{j}', para_next[j])
            self.load_current_parameters(init=False)
            self.run_preprocessor()
            self.compile_zdp()
            self.run_simulation()
            
            y_next, t_time = self.err_calculation()

            # update dataset
            if t_time > 16.96:
                X = np.vstack([X,x_next])
                y = np.append(y,y_next)

                if y_next < y_best:
                    bounds = [[num * (1-self.lr), num * (1+self.lr)] for num in x_next]
            else:
                X = np.vstack([X,x_next])
                y = np.append(y,self.max_error)
                
            gp_model = self.GP_Model(X,y)
            # db_set.csv 파일에서 f 파라미터에 해당하는 index를 불러오기
            df = pd.read_csv('db_set.csv')

            res_bo = {f'f{k}': para_next[k] for k in range(len(para_next))}
            res_bo['err'] = y[-1]
            res_bo['index'] = f'Bayesian_Opt {i}'
            print(f'State: Bayesian_Opt , iteration: {i}, error = {y_next}, t_time = {t_time}')
            df = pd.concat([df, pd.DataFrame([res_bo])], ignore_index=True)
            df.to_csv('db_set.csv', index=False)
        
if __name__ == "__main__":
    optimizer = ParameterOptmizer()
    # iniitialization
    optimizer.load_current_parameters(init=True)
    param = np.array(list(optimizer.current_parameters.values()))
    for i in range(len(param)):
        optimizer.modify_parameter(f'f{i}', param[i])
    optimizer.run_preprocessor()
    optimizer.compile_zdp()
    optimizer.run_simulation()

    res_init = optimizer.current_parameters.copy()
    res_init['err'] = optimizer.err_calculation()[0]
    res_init['index'] = 'initialization'

    df_set = pd.DataFrame([res_init])
    df_set.to_csv('db_set.csv', index=False)
    print(f'State: initialization Done')

    # LH sample pre-train data
    param_LH = 10**optimizer.LH_sampling()
    param_init = np.array(list(optimizer.current_parameters.values()))
    for j in range(len(optimizer.LH_sampling())):
        param = param_LH[j] * param_init
        for i in range(len(param)):
            optimizer.modify_parameter(f'f{i}', param[i])
        optimizer.load_current_parameters(init=False)
        optimizer.run_preprocessor()
        optimizer.compile_zdp()
        optimizer.run_simulation()

        res = optimizer.current_parameters.copy()
        total_time = optimizer.err_calculation()[1]
        
        if total_time > 16.96:
            res['err'] = optimizer.err_calculation()[0]  
        else:
            res['err'] = optimizer.max_error
        res['index'] = f'LH pre-train {j}'
    
        print(f'State: LH pre-train, iteration: {j}, error = {res['err']}, t_time = {total_time}')
        df_set = pd.concat([df_set, pd.DataFrame([res])], ignore_index=True)
        df_set.to_csv('db_set.csv', index=False)

    pre_tr_data = pd.read_csv('db_set.csv')
    # define the input and output data
    pre_X = pre_tr_data.iloc[:, :37].values
    pre_y = pre_tr_data['err'].values
    ori_X = pre_X[0]

    optimizer.Bayesian_Opt(pre_X,pre_y,ori_X)