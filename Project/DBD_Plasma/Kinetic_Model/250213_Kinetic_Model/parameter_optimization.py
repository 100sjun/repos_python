import subprocess
import pandas as pd
import numpy as np
import re
from typing import Dict, List, Tuple
import copy

class ParameterOptimizer:
    def __init__(self):
        # 실험 데이터 설정
        self.exp_list = ['H2', 'CH4', 'C2H6', 'C2H4', 'C2H2', 'C3H8', 'C3H6', 'C4H10', 'C5H12', 
                        'C', 'CH', 'CH2', 'CH3', 'C2H3', 'C2H5', 'C3H7', 'H',
                        'CH3^+', 'CH4^+', 'CH5^+', 'C2H2^+', 'C2H4^+', 'C2H5^+', 'C2H6^+', 'C3H6^+', 'C3H8^+']
        
        self.exp_values = [5.319082, 90.80403, 2.428802, 0.197735, 0.171795, 0.717088, 0.046734, 
                          0.114829, 0.119573, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        
        # 현재 파라미터 값들을 저장할 딕셔너리
        self.current_parameters = {}
        self.load_current_parameters()

    def load_current_parameters(self):
        """kinet.inp 파일에서 현재 파라미터 값들을 읽어옴"""
        with open('kinet.inp', 'r') as f:
            content = f.read()
            # 정규표현식으로 파라미터 값 추출
            pattern = r'parameter :: (f\d+) = ([\d.]+d[-+]?\d+)'
            matches = re.finditer(pattern, content)
            for match in matches:
                param_name, value = match.groups()
                self.current_parameters[param_name] = float(value.replace('d', 'e'))

    def modify_parameter(self, param_name: str, new_value: float) -> None:
        """kinet.inp 파일의 특정 파라미터 값을 수정"""
        with open('kinet.inp', 'r') as f:
            content = f.read()
        
        # 파라미터 값 수정
        pattern = f'(parameter :: {param_name} = )([\d.]+d[-+]?\d+)'
        new_value_str = f'{new_value:.4e}'.replace('e', 'd')
        content = re.sub(pattern, f'\\1{new_value_str}', content)
        
        with open('kinet.inp', 'w') as f:
            f.write(content)

    def run_simulation(self) -> pd.DataFrame:
        """시뮬레이션 실행 및 결과 반환"""
        # manual.py 실행
        subprocess.run(['python', 'manual.py'], capture_output=True)
        
        # 결과 파일 읽기
        species = []
        with open('qt_species_list.txt', 'r') as f:
            for line in f:
                comp = line[2:]
                species.append(comp.strip())
        
        df_sp = pd.read_csv('qt_densities.txt', sep=r'\s+', header=0, names=['Time [s]']+species)
        
        # 농도 계산
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

        t = abs(df_sp['Time [s]']-16.96).argmin()

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

        result_df = pd.DataFrame({
            'species': self.exp_list,
            'exp': self.exp_values,
            'sim': [sim_H2, sim_CH4, sim_C2H6, sim_C2H4, sim_C2H2, sim_C3H8, sim_C3H6, sim_C4H10, sim_C5H12,
                    sim_C, sim_CH, sim_CH2, sim_CH3, sim_C2H3, sim_C2H5, sim_C3H7, sim_H,
                    sim_CH3_plus, sim_CH4_plus, sim_CH5_plus, sim_C2H2_plus, sim_C2H4_plus,
                    sim_C2H5_plus, sim_C2H6_plus, sim_C3H6_plus, sim_C3H8_plus]
        })
        
        return result_df

    def calculate_error(self, sim_results: pd.DataFrame) -> float:
        """실험값과 시뮬레이션 결과의 차이 계산"""
        errors = []
        for species, exp_value in zip(self.exp_list, self.exp_values):
            if exp_value > 0:  # 실험값이 0인 경우 제외
                sim_value = sim_results.loc[sim_results['species'] == species, 'sim'].iloc[0]
                relative_error = abs(exp_value - sim_value) / exp_value
                errors.append(relative_error)
        return np.mean(errors)

    def optimize_parameters(self, iterations: int = 10):
        """파라미터 최적화 실행"""
        best_error = float('inf')
        best_parameters = copy.deepcopy(self.current_parameters)
        
        for i in range(iterations):
            print(f"Iteration {i+1}/{iterations}")
            
            # 각 파라미터에 대해 10배 증가/감소 시도
            for param_name in self.current_parameters:
                # 10배 증가 시도
                self.modify_parameter(param_name, self.current_parameters[param_name] * 10)
                results = self.run_simulation()
                error_increase = self.calculate_error(results)
                
                # 10배 감소 시도
                self.modify_parameter(param_name, self.current_parameters[param_name] / 10)
                results = self.run_simulation()
                error_decrease = self.calculate_error(results)
                
                # 가장 좋은 결과 선택
                min_error = min(best_error, error_increase, error_decrease)
                if min_error < best_error:
                    best_error = min_error
                    if min_error == error_increase:
                        best_parameters[param_name] *= 10
                    elif min_error == error_decrease:
                        best_parameters[param_name] /= 10
                
                # 원래 값으로 복구
                self.modify_parameter(param_name, self.current_parameters[param_name])
            
            # 최적의 파라미터 적용
            self.current_parameters = copy.deepcopy(best_parameters)
            for param_name, value in best_parameters.items():
                self.modify_parameter(param_name, value)
            
            print(f"Current best error: {best_error}")
            print("Current best parameters:")
            for param_name, value in best_parameters.items():
                print(f"{param_name}: {value:.2e}")

if __name__ == "__main__":
    optimizer = ParameterOptimizer()
    optimizer.optimize_parameters(iterations=5)