import pandas as pd
import re
import subprocess
import numpy as np
import optuna
from optuna.samplers import TPESampler
import sqlite3
import time

# load kinetic data
df_kinetic = pd.read_csv('kinetic_data.csv')
df_mol = pd.read_csv('kinetic_data.csv').filter(like='mol%')
df_power = pd.read_csv('kinetic_data.csv').filter(like='Power')

class ParameterOptimizer:
    def __init__(self):
        
        # define the experimental data
        self.exp_list = ['H2', 'CH4', 'C2H6', 'C2H4', 'C2H2', 'C3H8', 'C3H6', 'C4H10', 'C5H12', 'C']
        self.exp_values = [df_mol.iloc[0].tolist(),
                           df_mol.iloc[1].tolist(),
                           df_mol.iloc[2].tolist()]

        # dictionary to save the current parameters
        self.ori_para = {}
        self.init_para()  

        # define manipulated variables
        self.lr = 0.01 # bounds learning rate
        self.max_error = 1000 # maximum error when simulation is failed
        self.n_trials = 100 # number of optimization trials

        # dynamic bounds for parameters
        self.bnds = {}
        self.init_bnds()

        # SQLite database setup (separate from Optuna's database)
        self.db_path = 'detailed_results.db'  # detailed results
        self.optuna_db_path = 'optuna_study.db'  # Optuna study
        self.setup_database()

        # 기존 최적화 결과가 있으면 bounds 업데이트
        self.load_best_bnds_from_db()
        
    
    def init_bnds(self):
        '''initialize parameter bounds based on learning rate'''
        low = np.zeros(len(self.ori_para)) - self.lr
        high = np.zeros(len(self.ori_para)) + self.lr
        self.bnds = np.vstack((low, high)).T
    
    def update_bnds(self, best_params):
        '''update parameter bounds based on best parameters found so far'''
        for i in range(len(self.bnds)):
            value = best_params[i]
            new_low = value * (1 - self.lr)
            new_high = value * (1 + self.lr)

            # update the bounds with the new values compared to the current bounds
            self.bnds[i] = [min(self.bnds[i][0], new_low), max(self.bnds[i][1], new_high)]

    def load_best_bnds_from_db(self):
        '''update parameter bounds based on best parameters in the database'''
        try:
            conn = sqlite3.connect(self.db_path, timeout=30.0)
            cursor = conn.cursor()

            # looking for the trial with the lowest error
            cursor.execute("SELECT * FROM optimization_results WHERE error_value = (SELECT MIN(error_value) FROM optimization_results)")
            best_row = cursor.fetchone()

            if best_row:
                # get the column names
                cursor.execute("PRAGMA table_info(optimization_results)")
                columns = [col[1] for col in cursor.fetchall()]

                # extract parameter values
                best_params = {}
                for i, col_name in enumerate(columns):
                    if col_name.startswith('p') and col_name in self.ori_para:
                        best_params[col_name] = best_row[i]
                
                if best_params:
                    print(f"found the best performance! error value: {best_row[columns.index('error_value')]}")
                    print("update the bounds based on the best parameters")
                    self.update_bnds(best_params)

            conn.close()

        except sqlite3.OperationalError:
            print("no previous optimization results found. start with initial bounds")
        except Exception as e:
            print(f'error loading best bounds from datbase: {e}')
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
            param_columns = ', '.join([f'{param} REAL' for param in self.ori_para.keys()])
            
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
            
    def init_para(self):
        with open('kinet_ori.inp', 'r') as f:
            content = f.readlines()
       
        for line in content:
            if '$ double precision, parameter :: p' in line:
                # extract p parameters and values
                pattern = r'p(\d+)\s*=\s*([\d\.\+\-d]+)'
                match = re.search(pattern, line)
                if match:
                    param_name = 'p' + match.group(1)
                    value = float(match.group(2).replace('d','e'))
                    self.ori_para[param_name] = value
                    
    def modify_parameter(self, param_name: str, new_value: float, Voltage: float) -> None:
        # modify the parameter value in the kinet.inp file
        with open(f'kinet{Voltage}.inp', 'r') as f:
            content = f.read()

        # modify the parameter
        pattern = rf'(parameter :: {param_name} = )([\d\.\+\-d]+)'
        new_value_str = f'{new_value:.4e}'.replace('e', 'd')
        content = re.sub(pattern, f'parameter :: {param_name} = {new_value_str}', content)
        
        with open(f'kinet{Voltage}.inp', 'w') as f:
            f.write(content)

    def run_prep(self, Voltage):
        kinet_path = f'kinet{Voltage}.inp'
        process = subprocess.Popen(
            './preprocessor.exe',
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=False)
        
        process.stdin.write(f'{kinet_path}\n'.encode())
        process.stdin.flush()

        process.stdin.write('.\n'.encode())
        process.stdin.flush()

        output, error = process.communicate()
        output_str = output.decode('utf-8', errors='ignore') if output else ''
        error_str = error.decode('utf-8', errors='ignore') if error else ''

        return output_str, error_str
    
    def compile_zdp(self, Voltage):
        exe_path = f'run_{Voltage}kV.exe'
        compile_command= [
            'gfortran', '-o', exe_path, 'dvode_f90_m.F90', 'zdplaskin_m.F90', f'run_{Voltage}kV.F90', 'bolsig_x86_64_g.dll'
        ]
        result = subprocess.run(compile_command, capture_output=True, text=True)

        if result.returncode != 0:
            raise Exception(f"{exe_path} 컴파일 실패")

    def run_zdplaskin(self,Voltage):
        exe_path = f'run_{Voltage}kV.exe'
        try:
            process = subprocess.Popen(
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
                    process.kill()
                    break

                if "WARNING: BOLSIG+ convergence failed" in output:
                    process.stdin.write('\n')
                    process.stdin.flush()

        except:
            pass
        return process
    
    def err_cal(self, Voltage):
        #read the result file
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

        t = abs(df_sp['Time [s]']-16.9646).argmin()

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

        out_H2 = sim_H2 - sim_CH2 - 0.5*sim_CH3 - 0.5*sim_CH3_plus + 0.5*sim_CH5_plus - 0.5*sim_C2H3 + 0.5*sim_C2H5 + 0.5*sim_C2H5_plus + 0.5*sim_C3H7 + 0.5*sim_H
        out_CH4 = sim_CH4 + sim_CH2 + sim_CH3 + sim_CH3_plus + sim_CH4_plus + sim_CH5_plus
        out_C2H6 = sim_C2H6 + sim_C2H6_plus
        out_C2H4 = sim_C2H4 + sim_C2H3 + sim_C2H5 + sim_C2H4_plus + sim_C2H5_plus
        out_C2H2 = sim_C2H2 + sim_C2H2_plus
        out_C3H8 = sim_C3H8 + sim_C3H8_plus
        out_C3H6 = sim_C3H6 + sim_C3H7 + sim_C3H6_plus
        out_C4H10 = sim_C4H10
        out_C5H12 = sim_C5H12
        out_C = sim_C
        out_total = out_H2 + out_CH4 + out_C2H6 + out_C2H4 + out_C2H2 + out_C3H8 + out_C3H6 + out_C4H10 + out_C5H12 + out_C
        newsim = [out_H2/out_total*100, out_CH4/out_total*100, out_C2H6/out_total*100, out_C2H4/out_total*100, out_C2H2/out_total*100, out_C3H8/out_total*100, out_C3H6/out_total*100, out_C4H10/out_total*100, out_C5H12/out_total*100, out_C/out_total*100]
        
        w_factor = [0, 1, 10, 100, 10, 10, 10, 10, 1, 0]
        err = 0

        if Voltage == 10:
            id_exp = 2
        elif Voltage == 12.5:
            id_exp = 1
        elif Voltage == 15:
            id_exp = 0

        for i in range(len(self.exp_values[id_exp])):
            err += w_factor[i] * ((self.exp_values[id_exp][i] - newsim[i]))**2

        return err, df_sp['Time [s]'].iloc[-1]
    
    def objective(self, trial):
        '''optuna objective function'''
        print(f'\n{'='*50}')
        print(f'Trial {trial.number} started')
        print(f'{'='*50}')

        # suggest parameters within current bounds
        alpha = {}
        for i, param_name in enumerate(self.ori_para.keys()):
            alpha[param_name] = trial.suggest_float(param_name, self.bnds[i][0], self.bnds[i][1])

        PD0 = 4.317
        Volt = [10, 12.5]
        terr = 0
        runtime = []

        try:
            for v in Volt:
                if v == 10:
                    PD = 2.64
                elif v == 12.5:
                    PD = 3.273
                
                para = self.ori_para.copy()
                for i, para_name in enumerate(para.keys()):
                    para[para_name] = np.exp(alpha[param_name]*(PD0-PD))
                
                print(f'set the parameters: {v} kV case')
                # modify parameters in the input file
                for para_name, value in para.items():
                    self.modify_parameter(para_name, value, v)
                
                print(f'Run {v}kV')
                self.run_prep(v)
                self.compile_zdp(v)
                self.run_zdplaskin(v)
                err,time = self.err_cal(v)
                terr += err
                runtime.append(time)
                print('')

            # check if simulation completed successfully
            if min(runtime) < 17.0:
                terr = self.max_error
                status = 'simulation_incomplete'
                print(f'simulation is incomplete (min time: {min(runtime):.2f}s)')
            else:
                status = 'success'
                print(f'simulation is completed (min time: {min(runtime):.2f}s)')

            print(f'Trial {trial.number} error: {terr:.6f}')

            # save results to database
            self.save_to_database(trial.number, alpha, terr, min(runtime), status)

        except Exception as e:
            terr = self.max_error
            sim_time = 0
            status = f"error: {str(e)}"
            print(f"Error in running Trial {trial.number}: {e}")
            # try to save the results even if error occurs
            try:
                self.save_to_database(trial.number, alpha, terr, sim_time, status)
            except Exception as db_error:
                print(f"Failed to save results to database: {db_error}")

        print(f"Trial {trial.number} completed\n")
        return terr
    
    def update_bounds_callback(self, study, trial):
        """Callback to update bounds when a new best trial is found"""
        if study.best_trial.number == trial.number:
            # Update bounds based on the new best parameters
            self.update_bnds(study.best_params)

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

        # check if there is any existing optimization
        n_existing_trials = len(study.trials)
        print(f'existing Optuna trials: {n_existing_trials} trials')

        # if there is any existing optimization, update bounds based on the best parameters
        if n_existing_trials > 0:
            try:
                best_params = study.best_params
                print(f'found the best performance! error value: {study.best_value}')
                print('update the bounds based on the best parameters')
                self.update_bnds(best_params)
            except:
                print('failed to load best parameters from database')
        
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
    print(f"Number of parameters: {len(optimizer.ori_para)}")
    print(f"Maximum trials: {optimizer.n_trials}")
    print(f"Learning rate: {optimizer.lr}")
    
    study = optimizer.optimize()
    
    print("\nOptimization completed!")
    print(f"Best trial: {study.best_trial.number}")
    print(f"Best error: {study.best_value}")
    print("Best parameters saved in database.")        