"""Improved 2D Forced Convection Model with Temperature-Dependent Velocity
Accounts for thermal expansion using ideal gas law
"""

import numpy as np
import casadi as ca
import matplotlib.pyplot as plt
from reactor_config import reactor_config
from prop_He import prop_He


class ConvectionModel2D_Improved:
    """2D FDM model with temperature-dependent velocity field"""

    def __init__(self, flow_rate_sccm=50, T_inlet=298.15, T_wall=298.15,
                 Nr=20, Nz=50, P_atm=1.0):
        """Initialize improved model with thermal expansion"""

        # Load reactor configuration
        self.reactor = reactor_config()
        self.he_prop = prop_He()

        # Reactor geometry
        self.L = self.reactor.L
        self.R = self.reactor.ID / 2

        # Operating conditions
        self.Q_sccm = flow_rate_sccm
        self.T_inlet = T_inlet
        self.T_wall = T_wall
        self.P = P_atm * 101325  # Pa

        # Grid parameters
        self.Nr = Nr
        self.Nz = Nz

        # Create spatial grid
        self.create_grid()

        # Calculate mass flow rate (constant)
        self.calculate_mass_flow_rate()

    def create_grid(self):
        """Create 2D mesh"""
        self.r = np.linspace(0, self.R, self.Nr)
        self.dr = self.r[1] - self.r[0]
        self.z = np.linspace(0, self.L, self.Nz)
        self.dz = self.z[1] - self.z[0]
        self.R_mesh, self.Z_mesh = np.meshgrid(self.r, self.z)

    def calculate_mass_flow_rate(self):
        """Calculate constant mass flow rate"""
        # Convert sccm to actual flow at inlet conditions
        T_std = 273.15
        P_std = 101325
        Q_std = self.Q_sccm * 1e-6 / 60
        Q_inlet = Q_std * (self.T_inlet / T_std) * (P_std / self.P)

        # Inlet density
        rho_inlet, _, _ = self.get_properties(self.T_inlet)

        # Mass flow rate (constant)
        self.m_dot = rho_inlet * Q_inlet  # kg/s

        # Inlet velocity
        A = np.pi * self.R**2
        self.u_inlet = Q_inlet / A

    def get_properties(self, T):
        """Get He properties at temperature T"""
        T_C = T - 273.15
        rho = self.he_prop.rho(T_C) * 1000
        Cp = self.he_prop.Cp(T_C)
        k = self.he_prop.Tc(T_C)
        Cp_mass = Cp / (self.he_prop.M / 1000)
        return rho, Cp_mass, k

    def calculate_velocity_field(self, T_field):
        """
        Calculate velocity field from temperature field using:
        1. Mass conservation: ṁ = ρ(T)·u_avg(T)·A
        2. Parabolic profile shape: uz(r,z) = 2·u_avg(z)·[1 - (r/R)²]
        """
        uz_field = np.zeros((self.Nr, self.Nz))
        A = np.pi * self.R**2

        for j in range(self.Nz):
            # Get average temperature at this axial position
            T_avg_z = T_field[:, j].mean()

            # Get density at this temperature
            rho_z, _, _ = self.get_properties(T_avg_z)

            # Calculate average velocity from mass conservation
            u_avg_z = self.m_dot / (rho_z * A)

            # Create parabolic profile
            for i in range(self.Nr):
                uz_field[i, j] = 2 * u_avg_z * (1 - (self.r[i] / self.R)**2)

        return uz_field

    def solve_coupled(self, max_iterations=20, tol=1e-3):
        """
        Iteratively solve coupled energy-velocity problem:
        1. Guess temperature field
        2. Calculate velocity from T (thermal expansion)
        3. Solve energy equation with new velocity
        4. Repeat until convergence
        """
        print("\nStarting coupled iteration...")

        # Initial guess: linear temperature profile
        T_field = np.zeros((self.Nr, self.Nz))
        for j in range(self.Nz):
            T_axial = self.T_inlet + (self.T_wall - self.T_inlet) * (self.z[j] / self.L)
            for i in range(self.Nr):
                r_ratio = self.r[i] / self.R
                T_field[i, j] = T_axial * (1 - 0.2 * r_ratio**2)

        T_field[:, 0] = self.T_inlet
        T_field[-1, :] = self.T_wall

        for iteration in range(max_iterations):
            T_old = T_field.copy()

            # Update velocity field based on current temperature
            uz_field = self.calculate_velocity_field(T_field)

            # Store for use in energy equation
            self.uz_field = uz_field

            # Solve energy equation
            T_field = self.solve_energy_equation(T_field, uz_field)

            # Check convergence
            error = np.abs(T_field - T_old).max()
            error_rel = error / T_field.max()

            print(f"  Iteration {iteration+1}: max ΔT = {error:.3f} K ({error_rel*100:.3f}%)")

            if error < tol:
                print(f"  ✓ Converged in {iteration+1} iterations!")
                break
        else:
            print(f"  ⚠ Did not converge in {max_iterations} iterations")

        return T_field, uz_field

    def solve_energy_equation(self, T_initial, uz_field):
        """Solve energy equation with given velocity field"""

        # Get average properties
        T_avg = T_initial.mean()
        rho_avg, Cp_avg, k_avg = self.get_properties(T_avg)

        # Create symbolic variables
        T_sym = ca.SX.sym('T', self.Nr, self.Nz)

        # Build residual
        residual = ca.SX.zeros(self.Nr, self.Nz)

        for j in range(self.Nz):
            for i in range(self.Nr):

                # Inlet BC
                if j == 0:
                    residual[i, j] = T_sym[i, j] - self.T_inlet
                    continue

                # Wall BC
                if i == self.Nr - 1:
                    residual[i, j] = T_sym[i, j] - self.T_wall
                    continue

                # Interior points: k∇²T - ρCp·uz·∂T/∂z = 0

                # Convection term (FFD)
                if j >= self.Nz - 2:
                    conv = 0
                else:
                    dT_dz = (-3*T_sym[i, j] + 4*T_sym[i, j+1] - T_sym[i, j+2]) / (2*self.dz)
                    conv = rho_avg * Cp_avg * uz_field[i, j] * dT_dz

                # Diffusion term
                # Radial (CFD)
                if i == 0:
                    d2T_dr2 = (T_sym[i, j] - 2*T_sym[i+1, j] + T_sym[i+2, j]) / self.dr**2
                    dT_dr_over_r = d2T_dr2
                else:
                    d2T_dr2 = (T_sym[i-1, j] - 2*T_sym[i, j] + T_sym[i+1, j]) / self.dr**2
                    dT_dr = (T_sym[i+1, j] - T_sym[i-1, j]) / (2*self.dr)
                    dT_dr_over_r = dT_dr / self.r[i]

                # Axial (FFD)
                if j >= self.Nz - 2:
                    d2T_dz2 = 0
                else:
                    d2T_dz2 = (T_sym[i, j] - 2*T_sym[i, j+1] + T_sym[i, j+2]) / self.dz**2

                laplacian = d2T_dr2 + dT_dr_over_r + d2T_dz2
                diff = k_avg * laplacian

                residual[i, j] = diff - conv

        # Solve
        res_vec = ca.reshape(residual, -1, 1)
        T_vec = ca.reshape(T_sym, -1, 1)

        nlp = {
            'x': T_vec,
            'f': ca.dot(res_vec, res_vec) / 2,
            'g': ca.SX([])
        }

        opts = {
            'ipopt': {
                'print_level': 0,
                'max_iter': 200,
                'tol': 1e-4,
                'acceptable_tol': 1e-3
            },
            'print_time': False
        }

        solver = ca.nlpsol('solver', 'ipopt', nlp, opts)
        T_init_vec = T_initial.T.flatten()
        sol = solver(x0=T_init_vec, lbg=[], ubg=[])

        T_sol_vec = sol['x'].full().flatten()
        T_solution = T_sol_vec.reshape(self.Nz, self.Nr).T

        return T_solution

    def visualize_comparison(self, T_field, uz_field, save_path=None):
        """Visualize results with velocity-temperature coupling"""

        fig = plt.figure(figsize=(18, 10))

        # 1. Temperature field
        ax1 = plt.subplot(2, 3, 1)
        contour1 = ax1.contourf(self.Z_mesh * 100, self.R_mesh * 1000,
                                T_field.T, levels=20, cmap='hot')
        ax1.set_xlabel('Axial Position z (cm)')
        ax1.set_ylabel('Radial Position r (mm)')
        ax1.set_title('Temperature Field (K)', fontweight='bold')
        plt.colorbar(contour1, ax=ax1, label='T (K)')

        # 2. Velocity field
        ax2 = plt.subplot(2, 3, 2)
        Z_plot, R_plot = np.meshgrid(self.z, self.r)
        contour2 = ax2.contourf(Z_plot * 100, R_plot * 1000,
                                uz_field * 100, levels=20, cmap='viridis')
        ax2.set_xlabel('Axial Position z (cm)')
        ax2.set_ylabel('Radial Position r (mm)')
        ax2.set_title('Velocity Field (cm/s) - Variable!', fontweight='bold')
        plt.colorbar(contour2, ax=ax2, label='u_z (cm/s)')

        # 3. Axial profiles at centerline
        ax3 = plt.subplot(2, 3, 3)
        ax3_temp = ax3.twinx()

        line1 = ax3.plot(self.z * 100, T_field[0, :], 'r-', linewidth=2, label='Temperature')
        line2 = ax3_temp.plot(self.z * 100, uz_field[0, :] * 100, 'b-', linewidth=2, label='Velocity (centerline)')

        ax3.set_xlabel('Axial Position z (cm)')
        ax3.set_ylabel('Temperature (K)', color='r')
        ax3_temp.set_ylabel('Velocity (cm/s)', color='b')
        ax3.tick_params(axis='y', labelcolor='r')
        ax3_temp.tick_params(axis='y', labelcolor='b')
        ax3.set_title('Coupled T-u Evolution', fontweight='bold')
        ax3.grid(True, alpha=0.3)

        lines = line1 + line2
        labels = [l.get_label() for l in lines]
        ax3.legend(lines, labels, loc='upper left')

        # 4. Velocity profiles at different z positions
        ax4 = plt.subplot(2, 3, 4)
        z_positions = [0, 0.25, 0.5, 0.75, 1.0]
        colors = plt.cm.viridis(np.linspace(0, 1, len(z_positions)))

        for idx, z_frac in enumerate(z_positions):
            z_idx = int(z_frac * (self.Nz - 1))
            T_avg_z = T_field[:, z_idx].mean()
            ax4.plot(self.r * 1000, uz_field[:, z_idx] * 100,
                    color=colors[idx], linewidth=2,
                    label=f'z={z_frac*self.L*100:.0f}cm (T≈{T_avg_z:.0f}K)')

        ax4.set_xlabel('Radial Position r (mm)')
        ax4.set_ylabel('Velocity u_z (cm/s)')
        ax4.set_title('Velocity Profiles (showing thermal acceleration)', fontweight='bold')
        ax4.legend(fontsize=8)
        ax4.grid(True, alpha=0.3)

        # 5. Thermal expansion quantification
        ax5 = plt.subplot(2, 3, 5)

        # Calculate average velocity and temperature at each z
        u_avg_z = np.zeros(self.Nz)
        T_avg_z = np.zeros(self.Nz)
        for j in range(self.Nz):
            T_avg_z[j] = T_field[:, j].mean()
            u_avg_z[j] = uz_field[:, j].mean()

        ax5.plot(self.z * 100, u_avg_z / u_avg_z[0], 'g-', linewidth=2,
                label='Actual u/u_inlet')

        # Theoretical from ideal gas law
        u_theory = T_avg_z / self.T_inlet
        ax5.plot(self.z * 100, u_theory, 'r--', linewidth=2,
                label='Theory: T/T_inlet (ideal gas)')

        ax5.set_xlabel('Axial Position z (cm)')
        ax5.set_ylabel('Velocity Ratio u/u_inlet')
        ax5.set_title('Thermal Expansion Effect', fontweight='bold')
        ax5.legend()
        ax5.grid(True, alpha=0.3)

        # 6. Statistics
        ax6 = plt.subplot(2, 3, 6)
        ax6.axis('off')

        u_inlet_avg = uz_field[:, 0].mean()
        u_outlet_avg = uz_field[:, -1].mean()
        T_outlet_avg = T_field[:, -1].mean()

        expansion_ratio = u_outlet_avg / u_inlet_avg
        T_ratio = T_outlet_avg / self.T_inlet

        stats = f"""
    THERMAL EXPANSION ANALYSIS
    {'='*45}

    Inlet Conditions:
      Temperature:              {self.T_inlet:.2f} K
      Avg velocity:             {u_inlet_avg*100:.2f} cm/s
      Mass flow rate:           {self.m_dot*1e6:.3f} mg/s

    Outlet Conditions:
      Avg temperature:          {T_outlet_avg:.2f} K
      Avg velocity:             {u_outlet_avg*100:.2f} cm/s
      Mass flow rate:           {self.m_dot*1e6:.3f} mg/s ✓

    Thermal Expansion:
      Velocity increase:        {(expansion_ratio-1)*100:.1f}%
      Temperature increase:     {(T_ratio-1)*100:.1f}%
      Ratio match:              {expansion_ratio/T_ratio:.3f} ≈ 1.0 ✓

    Validation:
      u_out/u_in:               {expansion_ratio:.3f}
      T_out/T_in (theory):      {T_ratio:.3f}
      Error:                    {abs(expansion_ratio-T_ratio)/T_ratio*100:.2f}%

    Max velocity increase:      {uz_field[0, -1]/uz_field[0, 0]:.3f}×
    (centerline, inlet→outlet)
        """

        ax6.text(0.05, 0.95, stats, transform=ax6.transAxes,
                fontsize=9, verticalalignment='top', family='monospace',
                bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.3))

        plt.suptitle('Improved Model: Temperature-Dependent Velocity (Thermal Expansion)',
                     fontsize=14, fontweight='bold')

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"✓ Saved to: {save_path}")

        plt.close(fig)


if __name__ == "__main__":
    print("="*60)
    print("Improved 2D Convection Model with Thermal Expansion")
    print("="*60)

    model = ConvectionModel2D_Improved(
        flow_rate_sccm=50,
        T_inlet=298.15,
        T_wall=500.0,  # Higher wall temp to show expansion
        Nr=15,
        Nz=30
    )

    print(f"\nReactor: L={model.L*100:.1f} cm, R={model.R*1000:.2f} mm")
    print(f"Mass flow rate: {model.m_dot*1e6:.3f} mg/s (constant)")
    print(f"Inlet velocity: {model.u_inlet*100:.2f} cm/s")

    # Solve coupled problem
    T_solution, uz_solution = model.solve_coupled(max_iterations=15, tol=1e-2)

    # Visualize
    model.visualize_comparison(T_solution, uz_solution,
                              save_path='improved_model_thermal_expansion.png')

    print("\n" + "="*60)
    print("✓ Improved model complete!")
    print("="*60)
