# 2D Forced Convection Model for He Gas Flow

## Overview

This module implements a 2D forced convection model for helium (He) gas flow in a cylindrical reactor using the Finite Difference Method (FDM) with CasADi for symbolic computation and optimization.

## Features

- **Cylindrical coordinates (r, z)**: Axisymmetric geometry
- **2nd order Forward Finite Difference (FFD)**: For z-axis derivatives (flow direction)
- **2nd order Central Finite Difference (CFD)**: For r-axis derivatives (radial direction)
- **CasADi integration**: Symbolic computation and efficient optimization
- **Reactor configuration**: Automatically loads geometry from `reactor_config.py`
- **He properties**: Temperature-dependent properties from `prop_He.py`

## Model Equations

### Governing Equation

Steady-state energy equation with forced convection:

```
ρ·Cp·uz·∂T/∂z = k·∇²T
```

where:
- ρ: He density (kg/m³)
- Cp: Specific heat capacity (J/kg/K)
- uz: Axial velocity (m/s)
- k: Thermal conductivity (W/m/K)
- T: Temperature (K)

### Laplacian in Cylindrical Coordinates

```
∇²T = ∂²T/∂r² + (1/r)·∂T/∂r + ∂²T/∂z²
```

### Velocity Profile

Laminar parabolic profile for pipe flow:

```
uz(r) = 2·u_avg·[1 - (r/R)²]
```

where:
- u_avg = Q/A: Average velocity
- Q: Volumetric flow rate (m³/s)
- A = πR²: Cross-sectional area
- R: Reactor inner radius

## Discretization Schemes

### Forward Finite Difference (z-axis)

**First derivative (2nd order FFD):**
```
∂T/∂z ≈ (-3T[i,j] + 4T[i+1,j] - T[i+2,j]) / (2Δz)
```

**Second derivative:**
```
∂²T/∂z² ≈ (T[i,j] - 2T[i+1,j] + T[i+2,j]) / Δz²
```

### Central Finite Difference (r-axis)

**First derivative (2nd order CFD):**
```
∂T/∂r ≈ (T[i,j+1] - T[i,j-1]) / (2Δr)
```

**Second derivative:**
```
∂²T/∂r² ≈ (T[i,j-1] - 2T[i,j] + T[i,j+1]) / Δr²
```

### Centerline Singularity Treatment

At r=0, the term (1/r)·∂T/∂r is singular. We apply L'Hôpital's rule:

```
lim(r→0) [(1/r)·∂T/∂r] = ∂²T/∂r²
```

## Boundary Conditions

1. **Inlet (z=0)**: Dirichlet BC
   - T(r, 0) = T_inlet

2. **Outlet (z=L)**: Neumann BC
   - ∂T/∂z = 0 (convective outflow)

3. **Centerline (r=0)**: Symmetry BC
   - ∂T/∂r = 0

4. **Wall (r=R)**: Dirichlet BC
   - T(R, z) = T_wall

## Installation

### Dependencies

```bash
pip install numpy scipy matplotlib casadi
```

### Required Files

- `reactor_config.py`: Reactor geometry and properties
- `prop_He.py`: Helium gas properties

## Usage

### Basic Example

```python
from convection_model_2d import ConvectionModel2D

# Create model
model = ConvectionModel2D(
    flow_rate_sccm=50,    # He flow rate in sccm
    T_inlet=298.15,        # Inlet temperature (K)
    T_wall=400.0,          # Wall temperature (K)
    Nr=20,                 # Number of radial grid points
    Nz=50                  # Number of axial grid points
)

# Solve steady-state problem
T_solution = model.solve_steady_state()

# Visualize results
model.visualize(T_solution, save_path='results.png')
```

### Test Script

Run the provided test script:

```bash
python test_convection.py
```

## Model Parameters

### ConvectionModel2D Class

**Parameters:**
- `flow_rate_sccm` (float): He flow rate in standard cubic centimeters per minute
- `T_inlet` (float): Inlet temperature in Kelvin
- `T_wall` (float): Wall temperature in Kelvin
- `Nr` (int): Number of grid points in radial direction
- `Nz` (int): Number of grid points in axial direction
- `P_atm` (float): Pressure in atmospheres

**Key Attributes:**
- `L`: Reactor length (m)
- `R`: Reactor inner radius (m)
- `Q_actual`: Actual volumetric flow rate (m³/s)
- `u_avg`: Average velocity (m/s)
- `uz`: Velocity profile uz(r)

## Methods

### create_grid()
Creates 2D spatial mesh in cylindrical coordinates (r, z).

### calculate_velocity_profile()
Calculates axial velocity profile from flow rate using laminar pipe flow assumption.

### get_properties(T)
Returns He properties (ρ, Cp, k) at specified temperature.

### discretize_laplacian(T_sym)
Discretizes Laplacian operator using CFD (r-axis) and FFD (z-axis).

### discretize_convection(T_sym)
Discretizes convection term using FFD for z-derivative.

### build_residual(T_avg)
Builds residual equations for steady-state heat equation.

### solve_steady_state(T_initial, T_avg_guess)
Solves steady-state problem using CasADi's IPOPT solver.

**Returns:** Temperature field T(r,z) as numpy array (Nr × Nz)

### visualize(T_solution, save_path)
Generates visualization:
- 2D contour plot of temperature field
- Radial temperature profiles at different axial positions

## Example Results

For He flow at 50 sccm with:
- Inlet: 298.15 K (25°C)
- Wall: 400.0 K (127°C)
- Grid: 15 × 30 points

**Results:**
- Solution converges in ~10 iterations (~0.04 seconds)
- Temperature range: 265 K - 444 K
- Centerline outlet: ~399 K (126°C)
- Temperature rise: ~101 K

## Physical Validation

### Reynolds Number

```python
Re = ρ·u_avg·D / μ
```

For typical conditions:
- ρ ≈ 0.16 kg/m³
- u_avg ≈ 0.046 m/s
- D = 2R ≈ 5 mm
- μ ≈ 2×10⁻⁵ Pa·s

**Re ≈ 18** → Laminar flow (Re < 2300) ✓

### Péclet Number

```python
Pe = Re·Pr = (ρ·Cp·u_avg·D) / k
```

Indicates relative importance of convection vs. diffusion.

## File Structure

```
reactor_model/
├── convection_model_2d.py      # Main model implementation
├── test_convection.py           # Test script
├── reactor_config.py            # Reactor geometry class
├── prop_He.py                   # Helium properties class
├── README_convection.md         # This file
└── test_convection_results.png  # Example visualization
```

## Numerical Considerations

### Grid Resolution

- **Radial (Nr)**: Typically 15-30 points
  - Captures velocity gradient near wall
  - Resolves thermal boundary layer

- **Axial (Nz)**: Typically 30-100 points
  - Captures axial temperature development
  - Depends on reactor length and flow velocity

### CFL Condition (for transient extension)

```
CFL = u·Δt/Δz < 1
```

### Solver Options

IPOPT settings:
- `tol`: 1e-4 (convergence tolerance)
- `max_iter`: 500 (maximum iterations)
- `mu_strategy`: 'adaptive' (barrier parameter update)

## Extensions

### Possible Enhancements

1. **Transient model**: Add time-dependent term ∂T/∂t
2. **Non-uniform wall temperature**: T_wall = f(z)
3. **Heat generation**: Add source term Q_gen(r,z)
4. **Coupled momentum**: Solve full Navier-Stokes equations
5. **Turbulent flow**: Implement turbulence models for higher Re
6. **Multi-species**: Add species transport equations
7. **Radiation**: Include radiative heat transfer

### Transient Formulation

```python
# ∂T/∂t = (k·∇²T - ρ·Cp·uz·∂T/∂z) / (ρ·Cp)
# Use ca.integrator() with 'idas' or 'cvodes'
```

## References

1. Incropera, F.P., & DeWitt, D.P. (2002). *Fundamentals of Heat and Mass Transfer*
2. Patankar, S.V. (1980). *Numerical Heat Transfer and Fluid Flow*
3. Andersson, J., et al. (2019). *CasADi: A software framework for nonlinear optimization and optimal control*

## Authors

- Implementation: Claude Code + sjbaek
- Date: 2025-11-27

## License

Part of CO2 Electrification reactor modeling project.
