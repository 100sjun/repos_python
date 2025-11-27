"""반응기 관련 정보를 저장하는 class"""

# libraries
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

class reactor_config:
    def __init__(self):
        self.L = 50e-2 # m, length
        self.OD = 6.33e-3 # m, outer diameter
        self.ID = 5.03e-3 # m , inner diameter
        self.W = 57.449e-3 # kg, weight

    def linear(self, T, a, b):
        return a + b * T
    
    def log(self, T, a, b):
        return a + b * np.log(T)

    def rho(self):
        return self.W / (np.pi * (self.OD**2 - self.ID**2) / 4 * self.L) # kg/m3, density

    def Er(self, T):
        T_fit = np.concatenate(([20], np.arange(100, 1301, 100)))
        factor_fit = np.array([1., 1., 1.01, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.07, 1.07, 1.08, 1.08])
        params, _ = curve_fit(self.linear, T_fit, factor_fit)
        a, b = params
        return 1.35e-6 * self.linear(T, a, b) # ohm-m2/m, resistivity

    def Tc(self, T):
        T_fit = np.array([50, 600, 800, 1000, 1200])
        Tc_fit = np.array([11, 20, 22, 26, 27])
        params, _ = curve_fit(self.linear, T_fit, Tc_fit)
        a, b = params
        return self.linear(T, a, b) # W/m/K, thermal conductivity

    def Cp(self, T):
        T_fit = np.array([20, 200, 400, 800, 1000, 1200])
        Cp_fit = np.array([0.46, 0.56, 0.63, 0.71, 0.72, 0.74])
        params, _ = curve_fit(self.log, T_fit, Cp_fit)
        a, b = params
        return self.log(T, a, b) # kJ/kg/K, specific heat capacity