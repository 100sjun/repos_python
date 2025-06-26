import re
import pandas as pd
import subprocess
import numpy as np
import optuna
from optuna.samplers import TPESampler
import sqlite3
import time

# load kinetic data
df_kinetic = pd.read_csv('kinetic_data.csv')
df_mol = pd.read_csv('kinetic_data.csv').filter(like='mol%')
df_time = pd.read_csv('kinetic_data.csv').filter(like='Residence')

class ParameterOptimizer:
    def __init__(self):
        # define the experimental data
        self.exp_list = ['H2', 'CH4', 'C2H6', 'C2H4', 'C2H2', 'C3H8', 'C3H6', 'C4H10', 'C5H12', 
                        'C', 'CH', 'CH2', 'CH3', 'C2H3', 'C2H5', 'C3H7', 'H',
                        'CH3^+', 'CH4^+', 'CH5^+', 'C2H2^+', 'C2H4^+', 'C2H5^+', 'C2H6^+', 'C3H6^+', 'C3H8^+']
        self.exp_values = [df_mol.iloc[0].tolist() + [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                           df_mol.iloc[1].tolist() + [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                           df_mol.iloc[2].tolist() + [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                           df_mol.iloc[3].tolist() + [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]
        
        # define the file path
        self.kinet_path = 'kinet.inp'
        self.exe_path = 'run_15kV.exe'

        # dictionary to save the current parameters
        self.current_parameters = {}
        self.original_parameters = {}
        self.load_current_parameters(init=True)
        self.original_parameters = self.current_parameters.copy()

        # define manipulated variables
        self.max_error = 100
        self.lr = 0.05  # bounds learning rate
        self.n_trials = 100000  # number of optimization trials
        
        # dynamic bounds for parameters
        self.current_bounds = {}
        self.initialize_bounds()
        
        # SQLite database setup (separate from Optuna's database)
        self.db_path = 'detailed_results3.db'  # 별도의 상세 결과 저장용 DB
        self.optuna_db_path = 'optuna_study3.db'  # Optuna 전용 DB
        self.setup_database()
        
        # 기존 최적화 결과가 있으면 bounds 업데이트
        self.load_best_bounds_from_db()
    
    def initialize_bounds(self):
        """Initialize parameter bounds based on original values and learning rate"""
        for param_name, value in self.original_parameters.items():
            self.current_bounds[param_name] = {
                'low': value * (1 - self.lr),
                'high': value * (1 + self.lr)
            }
    
    def update_bounds(self, best_params):
        """Update parameter bounds based on best parameters found so far"""
        for param_name in self.current_bounds.keys():
            if param_name in best_params:
                value = best_params[param_name]
                new_low = value * (1 - self.lr)
                new_high = value * (1 + self.lr)
                
                # 기존 범위와 새로운 범위를 비교하여 더 넓은 범위로 업데이트
                self.current_bounds[param_name] = {
                    'low': min(self.current_bounds[param_name]['low'], new_low),
                    'high': max(self.current_bounds[param_name]['high'], new_high)
                }
    
    def load_best_bounds_from_db(self):
        """기존 데이터베이스에서 최고 성능의 파라미터를 찾아서 bounds 업데이트"""
        try:
            conn = sqlite3.connect(self.db_path, timeout=30.0)
            cursor = conn.cursor()
            
            # 가장 낮은 에러값을 가진 trial 찾기
            cursor.execute("SELECT * FROM optimization_results WHERE error_value = (SELECT MIN(error_value) FROM optimization_results)")
            best_row = cursor.fetchone()
            
            if best_row:
                # 컬럼 이름 가져오기
                cursor.execute("PRAGMA table_info(optimization_results)")
                columns = [col[1] for col in cursor.fetchall()]
                
                # 파라미터 값들 추출
                best_params = {}
                for i, col_name in enumerate(columns):
                    if col_name.startswith('f') and col_name in self.current_parameters:
                        best_params[col_name] = best_row[i]
                
                if best_params:
                    print(f"기존 최고 성능 발견! 에러값: {best_row[columns.index('error_value')]}")
                    print("탐색 범위를 최고 성능 파라미터 기준으로 업데이트합니다.")
                    self.update_bounds(best_params)
            
            conn.close()
            
        except sqlite3.OperationalError:
            # 테이블이 없거나 첫 실행인 경우
            print("기존 최적화 결과를 찾을 수 없습니다. 초기 범위로 시작합니다.")
        except Exception as e:
            print(f"기존 결과 로드 중 오류: {e}")
        finally:
            try:
                conn.close()
            except:
                pass
    
    def setup_database(self):
        """Setup SQLite database for storing optimization results"""
        try:
            conn = sqlite3.connect(self.db_path, timeout=30.0)
            cursor = conn.cursor()
            
            # Create parameters columns dynamically
            param_columns = ', '.join([f'{param} REAL' for param in self.current_parameters.keys()])
            
            cursor.execute(f'''
            CREATE TABLE IF NOT EXISTS optimization_results (
                trial_number INTEGER PRIMARY KEY,
                {param_columns},
                error_value REAL,
                simulation_time REAL,
                status TEXT,
                timestamp DATETIME DEFAULT CURRENT_TIMESTAMP
            )
            ''')
            
            conn.commit()
            conn.close()
        except Exception as e:
            print(f"데이터베이스 설정 오류: {e}")
        finally:
            try:
                conn.close()
            except:
                pass
    
    def save_to_database(self, trial_number, params, error_value, sim_time, status):
        """Save optimization results to SQLite database"""
        max_retries = 3
        retry_delay = 0.1
        
        for attempt in range(max_retries):
            try:
                conn = sqlite3.connect(self.db_path, timeout=30.0)
                cursor = conn.cursor()
                
                columns = ['trial_number'] + list(params.keys()) + ['error_value', 'simulation_time', 'status']
                values = [trial_number] + list(params.values()) + [error_value, sim_time, status]
                
                placeholders = ', '.join(['?' for _ in values])
                columns_str = ', '.join(columns)
                
                # INSERT OR REPLACE를 사용하여 중복 trial_number 문제 해결
                cursor.execute(f'''
                INSERT OR REPLACE INTO optimization_results ({columns_str})
                VALUES ({placeholders})
                ''', values)
                
                conn.commit()
                conn.close()
                return  # 성공하면 함수 종료
                
            except sqlite3.OperationalError as e:
                if "database is locked" in str(e) and attempt < max_retries - 1:
                    print(f"데이터베이스 잠금, {retry_delay}초 후 재시도... (시도 {attempt + 1}/{max_retries})")
                    time.sleep(retry_delay)
                    retry_delay *= 2  # 지수적 백오프
                    continue
                else:
                    print(f"데이터베이스 저장 오류: {e}")
                    break
            except Exception as e:
                print(f"데이터베이스 저장 중 예외 발생: {e}")
                break
            finally:
                try:
                    conn.close()
                except:
                    pass
    
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
        process = subprocess.Popen(
        './preprocessor.exe',
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=False)

        process.stdin.write(f'{self.kinet_path}\n'.encode())
        process.stdin.flush()

        process.stdin.write('.\n'.encode())
        process.stdin.flush()

        output, error = process.communicate()
        output_str = output.decode('utf-8', errors='ignore') if output else ''
        error_str = error.decode('utf-8', errors='ignore') if error else ''
    
        return output_str, error_str
    
    def compile_zdp(self):
        compile_command = [
        'gfortran', '-o', self.exe_path, 'dvode_f90_m.F90', 'zdplaskin_m.F90', 'run_15kV.F90', 'bolsig_x86_64_g.dll'
        ]
        result = subprocess.run(compile_command, capture_output=True, text=True)
    
        if result.returncode != 0:
            raise Exception(f"{self.exe_path} 컴파일 실패")

    def run_simulation(self):
        try:
            process = subprocess.Popen(
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

        t0 = abs(df_sp['Time [s]']-5.654867).argmin()
        t1 = abs(df_sp['Time [s]']-7.270543).argmin()
        t2 = abs(df_sp['Time [s]']-10.17876).argmin()
        t3 = abs(df_sp['Time [s]']-16.9646).argmin()

        sim = []
        final_out_values = []  # Store output values for all time points
        
        time_points = [5.654867, 7.270543, 10.17876, 16.9646]
        for i, t in enumerate([t0, t1, t2, t3]):
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
            
            # Store output values for each time point
            final_out_values.append({
                'time': time_points[i],
                'out_H2': out_H2/out_total*100,
                'out_CH4': out_CH4/out_total*100, 
                'out_C2H6': out_C2H6/out_total*100,
                'out_C2H4': out_C2H4/out_total*100,
                'out_C2H2': out_C2H2/out_total*100,
                'out_C3H8': out_C3H8/out_total*100,
                'out_C3H6': out_C3H6/out_total*100,
                'out_C4H10': out_C4H10/out_total*100,
                'out_C5H12': out_C5H12/out_total*100,
                'out_C': out_C/out_total*100
            })
            
        w_factor = [0, 1, 10, 100, 10, 10, 10, 10, 1, 0]
        err = 0
        for i in range(len(self.exp_values)):
            for j in range(len(self.exp_values[i])):
                if j < 10:
                    #err += w_factor[j] * ((self.exp_values[i][j] - sim[i][j])/self.exp_values[i][j])**2
                    err += w_factor[j] * ((self.exp_values[i][j] - sim[i][j]))**2
                else:
                    err += 0

        # Print the output values and experimental values for all time points
        print("\n모든 시간대의 출력값 비교:")
        print(f"{'종':>10}", end='')
        for time_data in final_out_values:
            print(f" {time_data['time']:>10.6f}s", end='')
        print()  # New line
        print("-" * (10 + 11 * len(final_out_values)))

        species_names = ['H2', 'CH4', 'C2H6', 'C2H4', 'C2H2', 'C3H8', 'C3H6', 'C4H10', 'C5H12', 'C']
        
        for i, species in enumerate(species_names):
            print(f"{species:>10}", end='')
            for time_idx, time_data in enumerate(final_out_values):
                print(f" {time_data[f'out_{species}']:>10.6f}", end='')
            print()  # New line
            print(f"{'(실험값)':>10}", end='')
            for time_idx in range(len(final_out_values)):
                print(f" {self.exp_values[time_idx][i]:>10.6f}", end='')
            print()  # New line
            if i < len(species_names) - 1:
                print("-" * (10 + 11 * len(final_out_values)))

        return err, df_sp['Time [s]'].iloc[-1]

    def objective(self, trial):
        """Optuna objective function"""
        print(f"\n{'='*50}")
        print(f"Trial {trial.number} 시작")
        print(f"{'='*50}")
        
        # Suggest parameters within current bounds
        params = {}
        for param_name in self.current_parameters.keys():
            bounds = self.current_bounds[param_name]
            params[param_name] = trial.suggest_float(
                param_name, 
                bounds['low'], 
                bounds['high']
            )
        
        print("파라미터 설정 중...")
        # Modify parameters in the input file
        for param_name, value in params.items():
            self.modify_parameter(param_name, value)
        
        # Update current parameters
        self.load_current_parameters(init=False)
        
        try:
            print("전처리기 실행 중...")
            # Run simulation pipeline
            self.run_preprocessor()
            
            print("컴파일 중...")
            self.compile_zdp()
            
            # Run simulation with trial number for progress display
            self.run_simulation()
            
            print("결과 계산 중...")
            # Calculate error
            error, sim_time = self.err_calculation()
            
            # Check if simulation completed successfully
            if sim_time < 17.0:
                error = self.max_error
                status = "simulation_incomplete"
                print(f"시뮬레이션 미완료 (시간: {sim_time:.2f}s)")
            else:
                status = "success"
                print(f"시뮬레이션 완료 (시간: {sim_time:.2f}s)")
            
            print(f"Trial {trial.number} 에러값: {error:.6f}")
            
            # Save results to database
            self.save_to_database(trial.number, params, error, sim_time, status)
            
        except Exception as e:
            error = self.max_error
            sim_time = 0
            status = f"error: {str(e)}"
            print(f"Trial {trial.number} 실행 중 오류 발생: {e}")
            # 오류가 발생해도 데이터베이스 저장 시도
            try:
                self.save_to_database(trial.number, params, error, sim_time, status)
            except Exception as db_error:
                print(f"데이터베이스 저장 실패: {db_error}")
                # 데이터베이스 저장이 실패해도 optimization은 계속 진행
        
        print(f"Trial {trial.number} 완료\n")
        return error
    
    def update_bounds_callback(self, study, trial):
        """Callback to update bounds when a new best trial is found"""
        if study.best_trial.number == trial.number:
            # Update bounds based on the new best parameters
            self.update_bounds(study.best_params)
    
    def optimize(self):
        """Run optimization using Optuna"""
        # Create study with TPE sampler and SQLite storage
        # 새로운 study 이름으로 충돌 방지
        import time
        study_name = f'parameter_optimization_{int(time.time())}'
        
        study = optuna.create_study(
            direction='minimize',
            sampler=TPESampler(seed=42),
            storage=f"sqlite:///{self.optuna_db_path}",
            study_name=study_name,
            load_if_exists=True
        )
        
    def optimize(self):
        """Run optimization using Optuna"""
        # Create study with TPE sampler and SQLite storage
        study = optuna.create_study(
            direction='minimize',
            sampler=TPESampler(seed=42, n_startup_trials=10),
            storage=f"sqlite:///{self.optuna_db_path}",
            study_name='parameter_optimization',
            load_if_exists=True
        )
        
        # 기존 최적화가 있었는지 확인
        n_existing_trials = len(study.trials)
        print(f"기존 Optuna trials: {n_existing_trials}개")
        
        # 기존 study에서 best params가 있으면 bounds 업데이트
        if n_existing_trials > 0:
            try:
                best_params = study.best_params
                print(f"Optuna에서 최고 성능 발견! 에러값: {study.best_value}")
                print("탐색 범위를 Optuna 최고 성능 파라미터 기준으로 업데이트합니다.")
                self.update_bounds(best_params)
            except:
                print("기존 Optuna 결과에서 best params를 가져올 수 없습니다.")
        
        # # Run initial evaluation only if no trials exist
        # if n_existing_trials == 0:
        #     print("초기 평가를 실행합니다...")
        #     for param_name, value in self.original_parameters.items():
        #         self.modify_parameter(param_name, value)
            
        #     self.run_preprocessor()
        #     self.compile_zdp() 
        #     self.run_simulation()
            
        #     initial_error, sim_time = self.err_calculation()
        #     status = "success" if sim_time >= 17.0 else "simulation_incomplete"
            
        #     self.save_to_database(0, self.original_parameters, initial_error, sim_time, "initial_evaluation")
        #     print(f"초기 평가 완료: error = {initial_error}")
        # else:
        #     print("기존 결과가 있으므로 초기 평가를 건너뜁니다.")
        
        # Optimize with callback for dynamic bounds update
        study.optimize(
            self.objective, 
            n_trials=self.n_trials,
            callbacks=[self.update_bounds_callback]
        )
        
        return study

if __name__ == "__main__":
    optimizer = ParameterOptimizer()
    
    print("Starting parameter optimization with Optuna...")
    print(f"Detailed results database: {optimizer.db_path}")
    print(f"Optuna database: {optimizer.optuna_db_path}")
    print(f"Number of parameters: {len(optimizer.current_parameters)}")
    print(f"Maximum trials: {optimizer.n_trials}")
    print(f"Learning rate: {optimizer.lr}")
    
    study = optimizer.optimize()
    
    print("\nOptimization completed!")
    print(f"Best trial: {study.best_trial.number}")
    print(f"Best error: {study.best_value}")
    print("Best parameters saved in database.")