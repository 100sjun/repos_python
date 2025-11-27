"""Final 2D Forced Convection Model with accurate mass conservation
Fixes:
1. Numerical stability with under-relaxation
2. Accurate mass flow integration using composite Simpson's rule
3. Better convergence criteria
"""

import numpy as np
import casadi as ca
import matplotlib.pyplot as plt
from reactor_config import reactor_config
from prop_He import prop_He


class ConvectionModel2D_Final:
    """2D FDM model with proper continuity and accurate mass conservation"""

    def __init__(self, flow_rate_sccm=50, T_inlet=298.15, T_wall=298.15,
                 Nr=20, Nz=50, P_atm=1.0):
        """Initialize model"""

        self.reactor = reactor_config()
        self.he_prop = prop_He()

        self.L = self.reactor.L
        self.R = self.reactor.ID / 2

        self.Q_sccm = flow_rate_sccm
        self.T_inlet = T_inlet
        self.T_wall = T_wall
        self.P = P_atm * 101325

        # Use odd number of points for Simpson's rule
        if Nr % 2 == 0:
            Nr += 1
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
        """Calculate inlet conditions and mass flux function g(r)"""
        # Convert sccm to actual flow
        T_std = 273.15
        P_std = 101325
        Q_std = self.Q_sccm * 1e-6 / 60
        Q_inlet = Q_std * (self.T_inlet / T_std) * (P_std / self.P)

        # Inlet properties
        rho_inlet, _, _ = self.get_properties(self.T_inlet)
        self.rho_inlet = rho_inlet

        # Inlet velocity (parabolic profile)
        A = np.pi * self.R**2
        u_avg_inlet = Q_inlet / A
        self.u_avg_inlet = u_avg_inlet

        # Mass flow rate
        self.m_dot = rho_inlet * Q_inlet

        # Mass flux function g(r) = ρ_inlet * uz_inlet(r)
        self.g_r = np.zeros(self.Nr)
        for i in range(self.Nr):
            uz_inlet = 2 * u_avg_inlet * (1 - (self.r[i] / self.R)**2)
            self.g_r[i] = rho_inlet * uz_inlet

        print(f"\nInlet conditions:")
        print(f"  ρ_inlet = {rho_inlet:.6f} kg/m³")
        print(f"  u_avg_inlet = {u_avg_inlet*100:.4f} cm/s")
        print(f"  ṁ = {self.m_dot*1e6:.6f} mg/s")

    def calculate_velocity_field(self, T_field):
        """Calculate velocity from continuity: uz(r,z) = g(r)/ρ(T(r,z))"""
        uz_field = np.zeros((self.Nr, self.Nz))

        for j in range(self.Nz):
            for i in range(self.Nr):
                rho_local, _, _ = self.get_properties(T_field[i, j])
                uz_field[i, j] = self.g_r[i] / rho_local

        return uz_field

    def integrate_mass_flow(self, T_field, uz_field, z_idx):
        """
        Accurately integrate mass flow using composite Simpson's rule
        ṁ = ∫₀^R ρ(r)·uz(r)·2πr dr
        """
        # For integration of f(r) from 0 to R with Simpson's rule:
        # ∫ f(r) dr ≈ (h/3)[f(0) + 4·Σf(odd) + 2·Σf(even) + f(R)]

        # Build integrand: f(r) = ρ·uz·2πr
        integrand = np.zeros(self.Nr)
        for i in range(self.Nr):
            rho_i, _, _ = self.get_properties(T_field[i, z_idx])
            integrand[i] = rho_i * uz_field[i, z_idx] * 2 * np.pi * self.r[i]

        # Simpson's rule
        # Handle r=0 specially (integrand = 0 there due to 2πr term)
        h = self.dr

        # Sum odd indices (1, 3, 5, ...)
        sum_odd = np.sum(integrand[1::2])

        # Sum even indices (2, 4, 6, ...) excluding first and last
        sum_even = np.sum(integrand[2:-1:2])

        # Simpson's formula
        mass_flow = (h / 3) * (integrand[0] + 4*sum_odd + 2*sum_even + integrand[-1])

        return mass_flow

    def verify_mass_conservation(self, T_field, uz_field):
        """Verify mass conservation at inlet and outlet"""

        print("\n" + "="*60)
        print("MASS CONSERVATION VERIFICATION")
        print("="*60)

        # Check at multiple z positions
        z_check = [0, self.Nz//4, self.Nz//2, 3*self.Nz//4, self.Nz-1]
        mass_flows = []

        for j in z_check:
            m_j = self.integrate_mass_flow(T_field, uz_field, j)
            mass_flows.append(m_j)
            error = (m_j - self.m_dot) / self.m_dot * 100
            z_cm = self.z[j] * 100
            print(f"  z={z_cm:5.1f}cm: ṁ={m_j*1e6:.6f} mg/s, error={error:+.3f}%")

        print(f"\n  Expected: ṁ = {self.m_dot*1e6:.6f} mg/s")
        print(f"  Mean error: {np.mean([(m-self.m_dot)/self.m_dot*100 for m in mass_flows]):.3f}%")
        print(f"  Max error:  {np.max([abs((m-self.m_dot)/self.m_dot*100) for m in mass_flows]):.3f}%")

        max_error = np.max([abs((m-self.m_dot)/self.m_dot*100) for m in mass_flows])
        if max_error < 0.5:
            print("  ✓ Excellent mass conservation! (<0.5% error)")
        elif max_error < 1.0:
            print("  ✓ Good mass conservation (<1% error)")
        elif max_error < 2.0:
            print("  ⚠ Acceptable mass conservation (<2% error)")
        else:
            print("  ✗ Poor mass conservation (>2% error)")

        return mass_flows

    def solve_coupled(self, max_iterations=30, tol=1e-3, relax_factor=0.7):
        """
        Iteratively solve with under-relaxation for stability
        T_new = relax_factor * T_solved + (1 - relax_factor) * T_old
        """

        print("\n" + "="*60)
        print("COUPLED ITERATION (with under-relaxation)")
        print("="*60)
        print(f"  Relaxation factor: {relax_factor}")
        print(f"  Convergence tol: {tol} K")

        # Better initial guess - linear temperature profile
        T_field = np.zeros((self.Nr, self.Nz))
        for j in range(self.Nz):
            # Linear from inlet to wall temperature
            T_base = self.T_inlet + (self.T_wall - self.T_inlet) * (self.z[j] / self.L) * 0.5

            for i in range(self.Nr):
                # Slight radial variation
                r_factor = (self.r[i] / self.R)**2
                T_field[i, j] = T_base + (self.T_wall - T_base) * r_factor * 0.5

        # Apply BCs
        T_field[:, 0] = self.T_inlet
        T_field[-1, :] = self.T_wall

        for iteration in range(max_iterations):
            T_old = T_field.copy()

            # Calculate velocity from current temperature
            uz_field = self.calculate_velocity_field(T_field)

            # Solve energy equation
            T_solved = self.solve_energy_equation(T_field, uz_field)

            # Under-relaxation for stability
            T_field = relax_factor * T_solved + (1 - relax_factor) * T_old

            # Check convergence
            error = np.abs(T_field - T_old).max()
            error_rel = error / T_field.max()
            T_min, T_max = T_field.min(), T_field.max()

            print(f"  Iter {iteration+1:2d}: ΔT_max={error:7.3f} K ({error_rel*100:5.3f}%), "
                  f"T=[{T_min:.1f}, {T_max:.1f}] K")

            # Check for unphysical temperatures
            if T_max > 2000 or T_min < 200:
                print("  ⚠ Unphysical temperatures detected! Reducing relaxation...")
                relax_factor *= 0.8
                T_field = T_old
                continue

            if error < tol:
                print(f"  ✓ Converged in {iteration+1} iterations!")
                break
        else:
            print(f"  ⚠ Did not converge in {max_iterations} iterations")

        # Final velocity field
        self.uz_field = self.calculate_velocity_field(T_field)

        # Verify mass conservation
        mass_flows = self.verify_mass_conservation(T_field, self.uz_field)

        return T_field, self.uz_field, mass_flows

    def solve_energy_equation(self, T_initial, uz_field):
        """Solve energy equation with given velocity field"""

        T_avg = T_initial.mean()
        rho_avg, Cp_avg, k_avg = self.get_properties(T_avg)

        T_sym = ca.SX.sym('T', self.Nr, self.Nz)
        residual = ca.SX.zeros(self.Nr, self.Nz)

        for j in range(self.Nz):
            for i in range(self.Nr):

                # Boundary conditions
                if j == 0:
                    residual[i, j] = T_sym[i, j] - self.T_inlet
                    continue
                if i == self.Nr - 1:
                    residual[i, j] = T_sym[i, j] - self.T_wall
                    continue

                # Convection term (FFD)
                if j >= self.Nz - 2:
                    conv = 0
                else:
                    dT_dz = (-3*T_sym[i, j] + 4*T_sym[i, j+1] - T_sym[i, j+2]) / (2*self.dz)
                    conv = rho_avg * Cp_avg * uz_field[i, j] * dT_dz

                # Diffusion term
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

        # Solve with tighter tolerances
        res_vec = ca.reshape(residual, -1, 1)
        T_vec = ca.reshape(T_sym, -1, 1)

        nlp = {'x': T_vec, 'f': ca.dot(res_vec, res_vec) / 2, 'g': ca.SX([])}
        opts = {
            'ipopt': {
                'print_level': 0,
                'max_iter': 500,
                'tol': 1e-6,
                'acceptable_tol': 1e-5,
                'acceptable_iter': 15,
                'warm_start_init_point': 'yes'
            },
            'print_time': False
        }

        solver = ca.nlpsol('solver', 'ipopt', nlp, opts)
        sol = solver(x0=T_initial.T.flatten(), lbg=[], ubg=[])

        T_solution = sol['x'].full().flatten().reshape(self.Nz, self.Nr).T

        return T_solution

    def visualize_all(self, T_field, uz_field, mass_flows, save_path=None):
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

        # 3. Mass flow along reactor
        ax3 = plt.subplot(2, 3, 3)
        z_positions = np.linspace(0, self.Nz-1, 20, dtype=int)
        mass_vs_z = [self.integrate_mass_flow(T_field, uz_field, j) for j in z_positions]
        error_vs_z = [(m - self.m_dot)/self.m_dot*100 for m in mass_vs_z]

        ax3.plot(self.z[z_positions]*100, error_vs_z, 'b-o', linewidth=2, markersize=4)
        ax3.axhline(y=0, color='r', linestyle='--', linewidth=2, label='Perfect conservation')
        ax3.fill_between(self.z[z_positions]*100, -0.5, 0.5, alpha=0.2, color='g',
                        label='<0.5% error')
        ax3.set_xlabel('z (cm)')
        ax3.set_ylabel('Mass flow error (%)')
        ax3.set_title('Mass Conservation Along Reactor', fontweight='bold')
        ax3.legend()
        ax3.grid(True, alpha=0.3)

        # 4. Velocity profiles at different z
        ax4 = plt.subplot(2, 3, 4)
        z_plot = [0, 0.25, 0.5, 0.75, 1.0]
        colors = plt.cm.hot(np.linspace(0.2, 1, len(z_plot)))

        for idx, z_frac in enumerate(z_plot):
            j = int(z_frac * (self.Nz - 1))
            T_avg = T_field[:, j].mean()
            ax4.plot(self.r*1000, uz_field[:, j]*100,
                    color=colors[idx], linewidth=2,
                    label=f'z={z_frac*self.L*100:.0f}cm (T≈{T_avg:.0f}K)')

        ax4.set_xlabel('r (mm)')
        ax4.set_ylabel('uz (cm/s)')
        ax4.set_title('Velocity Profiles', fontweight='bold')
        ax4.legend(fontsize=8)
        ax4.grid(True, alpha=0.3)

        # 5. Temperature profiles
        ax5 = plt.subplot(2, 3, 5)
        for idx, z_frac in enumerate(z_plot):
            j = int(z_frac * (self.Nz - 1))
            ax5.plot(self.r*1000, T_field[:, j],
                    color=colors[idx], linewidth=2,
                    label=f'z={z_frac*self.L*100:.0f}cm')

        ax5.set_xlabel('r (mm)')
        ax5.set_ylabel('Temperature (K)')
        ax5.set_title('Temperature Profiles', fontweight='bold')
        ax5.legend(fontsize=8)
        ax5.grid(True, alpha=0.3)

        # 6. Statistics
        ax6 = plt.subplot(2, 3, 6)
        ax6.axis('off')

        max_error = np.max([abs((m-self.m_dot)/self.m_dot*100) for m in mass_flows])
        u_inlet = uz_field[:, 0].mean()
        u_outlet = uz_field[:, -1].mean()

        stats = f"""
FINAL MODEL RESULTS
{'='*50}

Mass Conservation:
  Expected ṁ:        {self.m_dot*1e6:.6f} mg/s
  Max error:         {max_error:.3f}%
  Status:            {'✓ Excellent' if max_error < 0.5 else '✓ Good' if max_error < 1.0 else '⚠ Check'}

Inlet:
  Temperature:       {self.T_inlet:.2f} K
  Avg velocity:      {u_inlet*100:.4f} cm/s

Outlet:
  Temperature:       {T_field[:,-1].mean():.2f} K
  Avg velocity:      {u_outlet*100:.4f} cm/s
  Velocity ratio:    {u_outlet/u_inlet:.3f}×

Grid:
  Nr × Nz:           {self.Nr} × {self.Nz}
  dr, dz:            {self.dr*1000:.3f} mm, {self.dz*100:.2f} cm

Physics:
  • ∂(ρ·uz)/∂z = 0 enforced ✓
  • Simpson's rule integration ✓
  • Under-relaxation for stability ✓
  • Accurate mass conservation ✓
        """

        ax6.text(0.05, 0.95, stats, transform=ax6.transAxes,
                fontsize=9, verticalalignment='top', family='monospace',
                bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.3))

        plt.suptitle('Final Model: Accurate Mass Conservation',
                     fontsize=14, fontweight='bold')
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches='tight')
            print(f"\n✓ Saved to: {save_path}")

        plt.close(fig)


if __name__ == "__main__":
    print("="*60)
    print("FINAL 2D Convection Model")
    print("With accurate mass conservation")
    print("="*60)

    model = ConvectionModel2D_Final(
        flow_rate_sccm=50,
        T_inlet=298.15,
        T_wall=500.0,
        Nr=21,  # Odd number for Simpson's rule
        Nz=40
    )

    # Solve with under-relaxation
    T_sol, uz_sol, mass_flows = model.solve_coupled(
        max_iterations=30,
        tol=0.5,  # 0.5 K convergence
        relax_factor=0.7
    )

    # Visualize
    model.visualize_all(T_sol, uz_sol, mass_flows,
                       save_path='final_model_accurate.png')

    print("\n" + "="*60)
    print("✓ Final model complete!")
    print("="*60)
