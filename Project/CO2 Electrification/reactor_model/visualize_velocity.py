"""Visualize velocity profile for He gas flow in cylindrical reactor"""

import numpy as np
import matplotlib.pyplot as plt
from convection_model_2d import ConvectionModel2D

def visualize_velocity_profile():
    """Create comprehensive velocity profile visualizations"""

    # Create model
    print("Creating convection model...")
    model = ConvectionModel2D(
        flow_rate_sccm=50,
        T_inlet=298.15,
        T_wall=400.0,
        Nr=30,
        Nz=50
    )

    print(f"Reactor: L={model.L*100:.1f} cm, R={model.R*1000:.2f} mm")
    print(f"He flow: {model.Q_sccm} sccm = {model.Q_actual*1e6:.3f} mL/s")
    print(f"Average velocity: {model.u_avg*100:.2f} cm/s")
    print(f"Max velocity (centerline): {model.uz[0]*100:.2f} cm/s")
    print(f"Reynolds number: {calculate_reynolds(model):.1f}")

    # Create figure with multiple subplots
    fig = plt.figure(figsize=(16, 10))

    # 1. Radial velocity profile
    ax1 = plt.subplot(2, 3, 1)
    ax1.plot(model.r * 1000, model.uz * 100, 'b-', linewidth=2)
    ax1.axhline(y=model.u_avg * 100, color='r', linestyle='--',
                label=f'u_avg = {model.u_avg*100:.2f} cm/s')
    ax1.set_xlabel('Radial Position r (mm)', fontsize=11)
    ax1.set_ylabel('Axial Velocity u_z (cm/s)', fontsize=11)
    ax1.set_title('Radial Velocity Profile (Parabolic)', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    # 2. Normalized velocity profile
    ax2 = plt.subplot(2, 3, 2)
    r_normalized = model.r / model.R
    u_normalized = model.uz / model.uz[0]
    ax2.plot(r_normalized, u_normalized, 'g-', linewidth=2)
    ax2.set_xlabel('Normalized Radius r/R', fontsize=11)
    ax2.set_ylabel('Normalized Velocity u/u_max', fontsize=11)
    ax2.set_title('Normalized Velocity Profile', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim([0, 1])
    ax2.set_ylim([0, 1.05])

    # Add theoretical parabolic profile
    r_theory = np.linspace(0, 1, 100)
    u_theory = 1 - r_theory**2
    ax2.plot(r_theory, u_theory, 'r--', linewidth=1.5,
             label='Theory: 1 - (r/R)²', alpha=0.7)
    ax2.legend()

    # 3. 2D velocity field
    ax3 = plt.subplot(2, 3, 3)
    # uz_field has shape (Nz, Nr), need to transpose for contourf
    # Create meshgrid: Z_mesh and R_mesh have shape (Nz, Nr)
    contour = ax3.contourf(model.Z_mesh * 100, model.R_mesh * 1000,
                           model.uz_field * 100, levels=20, cmap='viridis')
    ax3.set_xlabel('Axial Position z (cm)', fontsize=11)
    ax3.set_ylabel('Radial Position r (mm)', fontsize=11)
    ax3.set_title('2D Velocity Field (cm/s)', fontsize=12, fontweight='bold')
    plt.colorbar(contour, ax=ax3, label='u_z (cm/s)')

    # 4. Volumetric flow rate distribution
    ax4 = plt.subplot(2, 3, 4)
    # Calculate flow rate through annular rings
    dr = model.dr
    dQ = np.zeros(model.Nr)
    for i in range(model.Nr):
        if i == 0:
            # Center point: treat as small circle
            r_avg = model.r[i] / 2
        else:
            r_avg = (model.r[i] + model.r[i-1]) / 2

        dA = 2 * np.pi * r_avg * dr
        dQ[i] = model.uz[i] * dA * 1e6  # Convert to mL/s

    ax4.bar(model.r * 1000, dQ, width=dr * 1000 * 0.8,
            color='steelblue', edgecolor='black', alpha=0.7)
    ax4.set_xlabel('Radial Position r (mm)', fontsize=11)
    ax4.set_ylabel('Flow Rate in Ring dQ (mL/s)', fontsize=11)
    ax4.set_title('Volumetric Flow Distribution', fontsize=12, fontweight='bold')
    ax4.grid(True, alpha=0.3, axis='y')

    total_Q = dQ.sum()
    ax4.text(0.98, 0.95, f'Total: {total_Q:.3f} mL/s\nExpected: {model.Q_actual*1e6:.3f} mL/s',
             transform=ax4.transAxes, fontsize=9,
             verticalalignment='top', horizontalalignment='right',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # 5. Velocity gradient (shear rate)
    ax5 = plt.subplot(2, 3, 5)
    # Calculate du/dr using central differences
    du_dr = np.zeros_like(model.uz)
    for i in range(1, model.Nr-1):
        du_dr[i] = (model.uz[i+1] - model.uz[i-1]) / (2 * model.dr)

    # Boundary conditions for gradient
    du_dr[0] = 0  # Symmetry at centerline
    du_dr[-1] = (model.uz[-1] - model.uz[-2]) / model.dr  # Forward diff at wall

    ax5.plot(model.r * 1000, -du_dr, 'r-', linewidth=2)  # Negative for positive shear
    ax5.set_xlabel('Radial Position r (mm)', fontsize=11)
    ax5.set_ylabel('Shear Rate |du/dr| (1/s)', fontsize=11)
    ax5.set_title('Velocity Gradient (Shear Rate)', fontsize=12, fontweight='bold')
    ax5.grid(True, alpha=0.3)

    # Highlight wall shear rate
    wall_shear = -du_dr[-1]
    ax5.axhline(y=wall_shear, color='orange', linestyle='--', alpha=0.7,
                label=f'Wall shear: {wall_shear:.2f} 1/s')
    ax5.legend()

    # 6. Flow statistics table
    ax6 = plt.subplot(2, 3, 6)
    ax6.axis('off')

    # Calculate statistics
    Re = calculate_reynolds(model)
    wall_shear_stress = calculate_wall_shear_stress(model)
    friction_factor = calculate_friction_factor(model)

    stats_text = f"""
    FLOW CHARACTERISTICS
    {'='*40}

    Geometry:
      Reactor length (L):        {model.L*100:.1f} cm
      Inner radius (R):          {model.R*1000:.2f} mm
      Cross-section area:        {np.pi*model.R**2*1e6:.2f} mm²

    Flow Rate:
      Standard (SCCM):           {model.Q_sccm:.1f} sccm
      Actual volumetric:         {model.Q_actual*1e6:.3f} mL/s
      Mass flow rate:            {model.Q_actual*0.1614:.6f} g/s

    Velocity:
      Average (u_avg):           {model.u_avg*100:.2f} cm/s
      Maximum (centerline):      {model.uz[0]*100:.2f} cm/s
      Ratio (u_max/u_avg):       {model.uz[0]/model.u_avg:.2f}

    Flow Regime:
      Reynolds number (Re):      {Re:.1f}
      Flow type:                 {'Laminar' if Re < 2300 else 'Turbulent'}

    Wall Effects:
      Wall shear rate:           {wall_shear:.2f} 1/s
      Wall shear stress:         {wall_shear_stress:.4f} Pa
      Friction factor (f):       {friction_factor:.4f}

    Residence Time:
      Average:                   {model.L/model.u_avg:.2f} s
    """

    ax6.text(0.1, 0.95, stats_text, transform=ax6.transAxes,
             fontsize=9, verticalalignment='top', family='monospace',
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))

    plt.suptitle(f'Velocity Profile Analysis - He Flow at {model.Q_sccm} sccm',
                 fontsize=14, fontweight='bold', y=0.995)

    plt.tight_layout()

    # Save figure
    save_path = 'velocity_profile_analysis.png'
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"\n✓ Saved to: {save_path}")

    plt.show()

    return model


def calculate_reynolds(model):
    """Calculate Reynolds number"""
    # Get He properties at inlet temperature
    rho, Cp, k = model.get_properties(model.T_inlet)

    # Dynamic viscosity for He at 298K (from literature)
    # μ ≈ 1.96e-5 Pa·s at 25°C
    mu = 1.96e-5  # Pa·s

    D = 2 * model.R  # Diameter
    Re = rho * model.u_avg * D / mu

    return Re


def calculate_wall_shear_stress(model):
    """Calculate wall shear stress"""
    # Get He properties
    rho, Cp, k = model.get_properties(model.T_inlet)
    mu = 1.96e-5  # Pa·s

    # For laminar pipe flow: τ_wall = 8μu_avg/R
    tau_wall = 8 * mu * model.u_avg / model.R

    return tau_wall


def calculate_friction_factor(model):
    """Calculate Darcy friction factor"""
    Re = calculate_reynolds(model)

    # For laminar flow: f = 64/Re
    if Re < 2300:
        f = 64 / Re
    else:
        # Blasius equation for turbulent flow
        f = 0.316 / Re**0.25

    return f


if __name__ == "__main__":
    print("="*60)
    print("Velocity Profile Visualization")
    print("="*60)

    model = visualize_velocity_profile()

    print("\n" + "="*60)
    print("✓ Velocity profile analysis complete!")
    print("="*60)
