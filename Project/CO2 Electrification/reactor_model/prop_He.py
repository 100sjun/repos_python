"""헬륨 관련 정보를 저장하는 class"""

# libraries
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

class prop_He:
    def __init__(self):
        self.M = 4.0026 # g/mol, molar mass

    def linear(self, T, a, b):
        return a + b * T
    
    def inverse(self, T, a, b):
        return a / (T + b)

    def Tc(self, T):
        T_fit = np.concatenate(([20], np.arange(100, 701, 100)))
        factor_fit = np.array([151, 179, 212, 244, 274, 302, 329, 353])
        params, _ = curve_fit(self.linear, T_fit, factor_fit)
        a, b = params
        return self.linear(T, a, b)/1000 # W/m/K, thermal conductivity

    def rho(self, T):
        T_fit = np.array([26.9, 46.9, 66.9, 86.9, 127, 227, 327, 427, 527, 627, 727])
        rho_fit = np.array([0.1604, 0.1504, 0.1415, 0.1337, 0.1203, 0.096926, 0.08022, 0.06876, 0.06017, 0.05348, 0.04814])
        params, _ = curve_fit(self.inverse, T_fit, rho_fit)
        a, b = params
        return self.inverse(T, a, b)/1000 # g/mL, mass density

    def Cp(self, T):
        TK = T + 273.15
        A = 20.78603
        B = 4.850638e-10
        C = -1.582916e-10
        D = 1.525102e-11
        E = 3.196347e-11
        return A + B*TK + C*TK**2 + D*TK**3 + E/TK**2 # J/mol/K, specific heat capacity
