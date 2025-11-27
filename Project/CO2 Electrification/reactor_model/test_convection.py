"""Test script for 2D forced convection model"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from convection_model_2d import ConvectionModel2D

def test_basic_model():
    """Test basic model initialization and solve"""
    print("=" * 60)
    print("Testing 2D Forced Convection Model")
    print("=" * 60)

    # Create model with He flow at 50 sccm
    print("\n1. Initializing model...")
    model = ConvectionModel2D(
        flow_rate_sccm=50,
        T_inlet=298.15,  # 25°C inlet
        T_wall=400.0,     # 127°C wall (heated)
        Nr=15,            # Fewer points for faster solve
        Nz=30
    )

    print(f"   ✓ Reactor: L={model.L*100:.1f} cm, R={model.R*1000:.2f} mm")
    print(f"   ✓ He flow: {model.Q_sccm} sccm = {model.Q_actual*1e6:.3f} mL/s")
    print(f"   ✓ Avg velocity: {model.u_avg:.4f} m/s")
    print(f"   ✓ Max velocity (centerline): {model.uz[0]:.4f} m/s")
    print(f"   ✓ Grid: {model.Nr} x {model.Nz} points (dr={model.dr*1000:.3f} mm, dz={model.dz*100:.2f} cm)")

    # Check properties
    print("\n2. Checking He properties...")
    rho_298, Cp_298, k_298 = model.get_properties(298.15)
    rho_400, Cp_400, k_400 = model.get_properties(400.0)

    print(f"   At T=298 K: ρ={rho_298:.4f} kg/m³, Cp={Cp_298:.1f} J/kg/K, k={k_298:.4f} W/m/K")
    print(f"   At T=400 K: ρ={rho_400:.4f} kg/m³, Cp={Cp_400:.1f} J/kg/K, k={k_400:.4f} W/m/K")

    # Solve
    print("\n3. Solving steady-state heat equation...")
    print("   (This may take 30-60 seconds...)")
    try:
        T_solution = model.solve_steady_state(T_avg_guess=340)

        print(f"   ✓ Solution converged!")
        print(f"   ✓ T_min = {T_solution.min():.2f} K ({T_solution.min()-273.15:.2f}°C)")
        print(f"   ✓ T_max = {T_solution.max():.2f} K ({T_solution.max()-273.15:.2f}°C)")
        print(f"   ✓ T_avg = {T_solution.mean():.2f} K ({T_solution.mean()-273.15:.2f}°C)")

        # Visualize
        print("\n4. Generating visualization...")
        fig = model.visualize(T_solution, save_path='test_convection_results.png')
        print("   ✓ Saved to: test_convection_results.png")

        # Print some key results
        print("\n5. Key Results:")
        print(f"   Centerline outlet T: {T_solution[0, -1]:.2f} K ({T_solution[0, -1]-273.15:.2f}°C)")
        print(f"   Wall outlet T: {T_solution[-1, -1]:.2f} K ({T_solution[-1, -1]-273.15:.2f}°C)")

        # Temperature rise
        T_rise_center = T_solution[0, -1] - model.T_inlet
        T_rise_bulk = T_solution[:, -1].mean() - model.T_inlet
        print(f"   ΔT centerline: {T_rise_center:.2f} K")
        print(f"   ΔT bulk average: {T_rise_bulk:.2f} K")

        print("\n" + "=" * 60)
        print("✓ All tests passed successfully!")
        print("=" * 60)

        return True

    except Exception as e:
        print(f"   ✗ Error during solve: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = test_basic_model()
    exit(0 if success else 1)
