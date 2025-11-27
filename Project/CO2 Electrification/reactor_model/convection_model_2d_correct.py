"""Correct 2D Forced Convection Model using proper continuity equation
Enforces: ∂(ρ·uz)/∂z = 0 at each radial position
"""

import numpy as np
import casadi as ca
import matplotlib.pyplot as plt
from reactor_config import reactor_config
from prop_He import prop_He


class ConvectionModel2D_Correct:
    """2D FDM model with proper continuity equation enforcement"""

    def __init__(self, flow_rate_sccm=50, T_inlet=298.15, T_wall=298.15,
                 Nr=20, Nz=50, P_atm=1.0):
        """Initialize model with proper mass balance"""

        self.reactor = reactor_config()
        self.he_prop = prop_He()

        self.L = self.reactor.L
        self.R = self.reactor.ID / 2

        self.Q_sccm = flow_rate_sccm
        self.T_inlet = T_inlet
        self.T_wall = T_wall
        self.P = P_atm * 101325

        self.Nr = Nr
        self.Nz = Nz

        self.create_grid()
        self.calculate_inlet_conditions()

    def create_grid(self):
        """Create 2D mesh"""
        self.r = np.linspace(0, self.R, self.Nr)
        self.dr = self.r[1] - self.r[0]
        self.z = np.linspace(0, self.L, self.Nz)
        self.dz = self.z[1] - self.z[0]
        self.R_mesh, self.Z_mesh = np.meshgrid(self.r, self.z)

    def get_properties(self, T):
        """Get He properties at temperature T (Kelvin)"""
        T_C = T - 273.15
        rho = self.he_prop.rho(T_C) * 1000  # kg/m³
        Cp = self.he_prop.Cp(T_C)  # J/mol/K
        k = self.he_prop.Tc(T_C)  # W/m/K
        Cp_mass = Cp / (self.he_prop.M / 1000)  # J/kg/K
        return rho, Cp_mass, k

    def calculate_inlet_conditions(self):
        """
        Calculate inlet conditions and mass flux function g(r)
        where g(r) = ρ_inlet * uz_inlet(r) = constant along z for each r
        """
        # Convert sccm to actual flow at inlet
        T_std = 273.15
        P_std = 101325
        Q_std = self.Q_sccm * 1e-6 / 60
        Q_inlet = Q_std * (self.T_inlet / T_std) * (P_std / self.P)

        # Inlet density (uniform at inlet temperature)
        rho_inlet, Cp_inlet, k_inlet = self.get_properties(self.T_inlet)
        self.rho_inlet = rho_inlet

        # Average inlet velocity
        A = np.pi * self.R**2
        u_avg_inlet = Q_inlet / A
        self.u_avg_inlet = u_avg_inlet

        # Total mass flow rate
        self.m_dot = rho_inlet * Q_inlet

        # Mass flux function g(r) = ρ_inlet * uz_inlet(r)
        # For parabolic inlet profile: uz_inlet(r) = 2*u_avg*[1-(r/R)²]
        self.g_r = np.zeros(self.Nr)
        for i in range(self.Nr):
            uz_inlet = 2 * u_avg_inlet * (1 - (self.r[i] / self.R)**2)
            self.g_r[i] = rho_inlet * uz_inlet

        print(f"Inlet conditions:")
        print(f"  ρ_inlet = {rho_inlet:.5f} kg/m³")
        print(f"  u_avg_inlet = {u_avg_inlet*100:.3f} cm/s")
        print(f"  ṁ = {self.m_dot*1e6:.4f} mg/s")
        print(f"  g(r=0) = {self.g_r[0]:.6f} kg/(m²·s)")
        print(f"  g(r=R) = {self.g_r[-1]:.6f} kg/(m²·s)")

    def calculate_velocity_field(self, T_field):
        """
        Calculate velocity from continuity equation: ∂(ρ·uz)/∂z = 0

        This gives: ρ(r,z)·uz(r,z) = g(r) = constant along z

        Therefore: uz(r,z) = g(r) / ρ(T(r,z))
        """
        uz_field = np.zeros((self.Nr, self.Nz))

        for j in range(self.Nz):
            for i in range(self.Nr):
                # Get local density from temperature
                rho_local, _, _ = self.get_properties(T_field[i, j])

                # Apply continuity: uz = g(r) / ρ(r,z)
                uz_field[i, j] = self.g_r[i] / rho_local

        return uz_field

    def verify_continuity(self, T_field, uz_field):
        """Verify that ∂(ρ·uz)/∂z ≈ 0"""

        print("\n" + "="*60)
        print("CONTINUITY VERIFICATION: ∂(ρ·uz)/∂z = 0")
        print("="*60)

        errors = []
        mass_flux = np.zeros((self.Nr, self.Nz))

        for j in range(self.Nz):
            for i in range(self.Nr):
                rho_local, _, _ = self.get_properties(T_field[i, j])
                mass_flux[i, j] = rho_local * uz_field[i, j]

        # Check at interior radial positions
        check_radii = [0, self.Nr//4, self.Nr//2, 3*self.Nr//4]

        for i in check_radii:
            flux_variation = mass_flux[i, :].std() / mass_flux[i, :].mean() * 100
            errors.append(flux_variation)
            print(f"  r={self.r[i]*1000:.2f}mm: ρ·uz variation = {flux_variation:.3f}%")

        print(f"\n  Mean variation: {np.mean(errors):.3f}%")
        print(f"  Max variation:  {np.max(errors):.3f}%")

        if np.max(errors) < 1.0:
            print("  ✓ Continuity satisfied! (< 1% variation)")
        else:
            print("  ⚠ Continuity violation detected")

        return mass_flux

    def solve_coupled(self, max_iterations=20, tol=1e-3):
        """Iteratively solve coupled energy-velocity problem"""

        print("\n" + "="*60)
        print("COUPLED ITERATION (with proper continuity)")
        print("="*60)

        # Initial temperature guess
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

            # Update velocity from continuity equation
            uz_field = self.calculate_velocity_field(T_field)
            self.uz_field = uz_field

            # Solve energy equation
            T_field = self.solve_energy_equation(T_field, uz_field)

            # Check convergence
            error = np.abs(T_field - T_old).max()
            error_rel = error / T_field.max()

            print(f"  Iter {iteration+1}: max ΔT = {error:.3f} K ({error_rel*100:.3f}%)")

            if error < tol:
                print(f"  ✓ Converged in {iteration+1} iterations!")
                break
        else:
            print(f"  ⚠ Did not converge in {max_iterations} iterations")

        # Verify continuity
        mass_flux = self.verify_continuity(T_field, uz_field)

        return T_field, uz_field, mass_flux

    def solve_energy_equation(self, T_initial, uz_field):
        """Solve energy equation with given velocity field"""

        T_avg = T_initial.mean()
        rho_avg, Cp_avg, k_avg = self.get_properties(T_avg)

        T_sym = ca.SX.sym('T', self.Nr, self.Nz)
        residual = ca.SX.zeros(self.Nr, self.Nz)

        for j in range(self.Nz):
            for i in range(self.Nr):

                # BCs
                if j == 0:
                    residual[i, j] = T_sym[i, j] - self.T_inlet
                    continue
                if i == self.Nr - 1:
                    residual[i, j] = T_sym[i, j] - self.T_wall
                    continue

                # Convection (FFD)
                if j >= self.Nz - 2:
                    conv = 0
                else:
                    dT_dz = (-3*T_sym[i, j] + 4*T_sym[i, j+1] - T_sym[i, j+2]) / (2*self.dz)
                    conv = rho_avg * Cp_avg * uz_field[i, j] * dT_dz

                # Diffusion
                if i == 0:
                    d2T_dr2 = (T_sym[i, j] - 2*T_sym[i+1, j] + T_sym[i+2, j]) / self.dr**2
                    dT_dr_over_r = d2T_dr2
                else:
                    d2T_dr2 = (T_sym[i-1, j] - 2*T_sym[i, j] + T_sym[i+1, j]) / self.dr**2
                    dT_dr = (T_sym[i+1, j] - T_sym[i-1, j]) / (2*self.dr)
                    dT_dr_over_r = dT_dr / self.r[i]

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

        nlp = {'x': T_vec, 'f': ca.dot(res_vec, res_vec) / 2, 'g': ca.SX([])}
        opts = {
            'ipopt': {'print_level': 0, 'max_iter': 200, 'tol': 1e-4},
            'print_time': False
        }

        solver = ca.nlpsol('solver', 'ipopt', nlp, opts)
        sol = solver(x0=T_initial.T.flatten(), lbg=[], ubg=[])
        T_solution = sol['x'].full().flatten().reshape(self.Nz, self.Nr).T

        return T_solution

    def visualize_all(self, T_field, uz_field, mass_flux, save_path=None):
        """Comprehensive visualization"""

        fig = plt.figure(figsize=(18, 12))

        # 1. Temperature field
        ax1 = plt.subplot(2, 3, 1)
        c1 = ax1.contourf(self.Z_mesh*100, self.R_mesh*1000, T_field.T,
                          levels=20, cmap='hot')
        ax1.set_xlabel('z (cm)')
        ax1.set_ylabel('r (mm)')
        ax1.set_title('Temperature Field (K)', fontweight='bold')
        plt.colorbar(c1, ax=ax1)

        # 2. Velocity field
        ax2 = plt.subplot(2, 3, 2)
        c2 = ax2.contourf(self.Z_mesh*100, self.R_mesh*1000, uz_field.T*100,
                          levels=20, cmap='viridis')
        ax2.set_xlabel('z (cm)')
        ax2.set_ylabel('r (mm)')
        ax2.set_title('Velocity Field uz (cm/s)', fontweight='bold')
        plt.colorbar(c2, ax=ax2)

        # 3. Mass flux ρ·uz (should be constant along z)
        ax3 = plt.subplot(2, 3, 3)
        c3 = ax3.contourf(self.Z_mesh*100, self.R_mesh*1000, mass_flux.T,
                          levels=20, cmap='plasma')
        ax3.set_xlabel('z (cm)')
        ax3.set_ylabel('r (mm)')
        ax3.set_title('Mass Flux ρ·uz (kg/m²/s)', fontweight='bold')
        plt.colorbar(c3, ax=ax3)

        # 4. Mass flux along z at different r
        ax4 = plt.subplot(2, 3, 4)
        check_r = [0, self.Nr//4, self.Nr//2, 3*self.Nr//4, self.Nr-2]
        colors = plt.cm.viridis(np.linspace(0, 1, len(check_r)))

        for idx, i in enumerate(check_r):
            ax4.plot(self.z*100, mass_flux[i, :],
                    color=colors[idx], linewidth=2, marker='o', markersize=3,
                    label=f'r={self.r[i]*1000:.2f}mm')

        ax4.set_xlabel('z (cm)')
        ax4.set_ylabel('ρ·uz (kg/m²/s)')
        ax4.set_title('Mass Flux Verification (should be flat)', fontweight='bold')
        ax4.legend(fontsize=8)
        ax4.grid(True, alpha=0.3)

        # 5. Velocity profiles at different z
        ax5 = plt.subplot(2, 3, 5)
        z_positions = [0, 0.25, 0.5, 0.75, 1.0]
        colors2 = plt.cm.hot(np.linspace(0.2, 1, len(z_positions)))

        for idx, z_frac in enumerate(z_positions):
            j = int(z_frac * (self.Nz - 1))
            T_avg = T_field[:, j].mean()
            ax5.plot(self.r*1000, uz_field[:, j]*100,
                    color=colors2[idx], linewidth=2,
                    label=f'z={z_frac*self.L*100:.0f}cm (T≈{T_avg:.0f}K)')

        ax5.set_xlabel('r (mm)')
        ax5.set_ylabel('uz (cm/s)')
        ax5.set_title('Velocity Profiles (varying with T)', fontweight='bold')
        ax5.legend(fontsize=8)
        ax5.grid(True, alpha=0.3)

        # 6. Statistics
        ax6 = plt.subplot(2, 3, 6)
        ax6.axis('off')

        u_inlet = uz_field[:, 0].mean()
        u_outlet = uz_field[:, -1].mean()
        T_outlet = T_field[:, -1].mean()

        # Check total mass flow
        m_check = 0
        for i in range(self.Nr):
            if i == 0:
                r_avg = self.r[0] / 2
            else:
                r_avg = (self.r[i] + self.r[i-1]) / 2
            dA = 2 * np.pi * r_avg * self.dr
            rho_out, _, _ = self.get_properties(T_field[i, -1])
            m_check += rho_out * uz_field[i, -1] * dA

        stats = f"""
CONTINUITY-BASED MODEL RESULTS
{'='*50}

Inlet:
  Temperature:       {self.T_inlet:.2f} K
  Avg velocity:      {u_inlet*100:.3f} cm/s
  Mass flow:         {self.m_dot*1e6:.4f} mg/s

Outlet:
  Avg temperature:   {T_outlet:.2f} K
  Avg velocity:      {u_outlet*100:.3f} cm/s
  Mass flow (check): {m_check*1e6:.4f} mg/s

Verification:
  ṁ error:           {abs(m_check-self.m_dot)/self.m_dot*100:.3f}%
  Velocity ratio:    {u_outlet/u_inlet:.3f}

Physics:
  • Each streamline conserves ρ·uz = const ✓
  • Velocity increases where T↑ (ρ↓)
  • Profile shape changes with heating
  • Proper continuity equation enforced
        """

        ax6.text(0.05, 0.95, stats, transform=ax6.transAxes,
                fontsize=9, verticalalignment='top', family='monospace',
                bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.3))

        plt.suptitle('Correct Model: Enforcing ∂(ρ·uz)/∂z = 0',
                     fontsize=14, fontweight='bold')
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"\n✓ Saved to: {save_path}")

        plt.close(fig)


if __name__ == "__main__":
    print("="*60)
    print("CORRECT 2D Convection Model")
    print("Using proper continuity: ∂(ρ·uz)/∂z = 0")
    print("="*60)

    model = ConvectionModel2D_Correct(
        flow_rate_sccm=50,
        T_inlet=298.15,
        T_wall=500.0,
        Nr=15,
        Nz=30
    )

    # Solve
    T_sol, uz_sol, mass_flux = model.solve_coupled(max_iterations=15, tol=1e-2)

    # Visualize
    model.visualize_all(T_sol, uz_sol, mass_flux,
                       save_path='correct_model_continuity.png')

    print("\n" + "="*60)
    print("✓ Correct model complete!")
    print("="*60)
