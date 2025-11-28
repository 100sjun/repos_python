'''
FDM solver for 1D Transient Coupled Gas-Wall Heat Transfer

Governing Equations:
- Gas (quasi-steady): mdot * Cp * dTg/dz = h * pi * Di * (Tw - Tg)
- Wall (transient): rho * A * Cp * dTw/dt = Qgen - h*pi*Di*(Tw - Tg) - qloss

Pytorch-based implementation for automatic differentiation
'''

# Libraries
import torch
import torch.nn as nn
import numpy as np
from Geometry import ReactorGeometry

# ====================================================================
# FDM Solver Class
# ====================================================================
class TheramlFDMSolver(nn.Module):
    '''
    1D Transient FDM solver for coupled gas-wall heat transfer

    Solves:
    - Gas: quasi-steady convection equation (marching from inlet)
    - Wall: trainsient energy balance with generation, convection, and loss
    '''

    def __init__(self, geometry: ReactorGeometry):
        