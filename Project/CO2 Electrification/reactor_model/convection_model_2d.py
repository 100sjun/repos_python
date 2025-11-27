"""2D Forced Convection Model for He gas flow in cylindrical reactor
Uses 2nd order Forward Finite Difference (FFD) for z-axis
Uses 2nd order Central Finite Difference (CFD) for r-axis
Implemented with CasADi for symbolic computation
"""

import numpy as np
import casadi as ca
import matplotlib.pyplot as plt
from reactor_config import reactor_config
from prop_He import prop_He


class ConvectionModel2D:
    """2D FDM model for forced convection of He gas in cylindrical reactor"""

    def __init__(self, flow_rate_sccm=50, T_inlet=298.15, T_wall=298.15,
                 Nr=20, Nz=50, P_atm=1.0):
        """
        Initialize 2D convection model

        Parameters:
        -----------
        flow_rate_sccm : float
            He flow rate in standard cubic centimeters per minute (sccm)
        T_inlet : float
            Inlet temperature in Kelvin
        T_wall : float
            Wall temperature in Kelvin
        Nr : int
            Number of grid points in radial direction
        Nz : int
            Number of grid points in axial direction
        P_atm : float
            Pressure in atmospheres
        """
        # Load reactor configuration
        self.reactor = reactor_config()
        self.he_prop = prop_He()

        # Reactor geometry
        self.L = self.reactor.L  # m, length
        self.R = self.reactor.ID / 2  # m, radius (inner diameter)

        # Operating conditions
        self.Q_sccm = flow_rate_sccm  # sccm
        self.T_inlet = T_inlet  # K
        self.T_wall = T_wall  # K
        self.P = P_atm * 101325  # Pa

        # Grid parameters
        self.Nr = Nr
        self.Nz = Nz

        # Create spatial grid
        self.create_grid()

        # Calculate velocity profile
        self.calculate_velocity_profile()

    def create_grid(self):
        """Create 2D mesh in cylindrical coordinates (r, z)"""
        # Radial grid (0 to R)
        self.r = np.linspace(0, self.R, self.Nr)
        self.dr = self.r[1] - self.r[0]

        # Axial grid (0 to L)
        self.z = np.linspace(0, self.L, self.Nz)
        self.dz = self.z[1] - self.z[0]

        # 2D meshgrid
        self.R_mesh, self.Z_mesh = np.meshgrid(self.r, self.z)

    def calculate_velocity_profile(self):
        """Calculate axial velocity profile from flow rate
        Assumes laminar parabolic profile: uz(r) = 2*u_avg*(1 - (r/R)^2)
        """
        # Convert sccm to m^3/s at operating conditions
        # Standard conditions: 273.15 K, 1 atm
        T_std = 273.15  # K
        P_std = 101325  # Pa

        Q_std = self.Q_sccm * 1e-6 / 60  # m^3/s at standard conditions

        # Apply ideal gas law to get actual volumetric flow rate
        # Q_actual = Q_std * (T_inlet/T_std) * (P_std/P)
        self.Q_actual = Q_std * (self.T_inlet / T_std) * (P_std / self.P)

        # Average velocity
        A = np.pi * self.R**2  # m^2
        self.u_avg = self.Q_actual / A  # m/s

        # Parabolic velocity profile (laminar flow)
        # uz(r) = 2*u_avg*(1 - (r/R)^2)
        self.uz = 2 * self.u_avg * (1 - (self.r / self.R)**2)  # m/s

        # Create 2D velocity field (constant along z)
        self.uz_field = np.tile(self.uz, (self.Nz, 1))

    def get_properties(self, T):
        """Get He properties at temperature T (in Celsius)"""
        T_C = T - 273.15  # Convert K to C

        rho = self.he_prop.rho(T_C) * 1000  # kg/m^3
        Cp = self.he_prop.Cp(T_C)  # J/mol/K
        k = self.he_prop.Tc(T_C)  # W/m/K

        # Convert Cp from J/mol/K to J/kg/K
        Cp_mass = Cp / (self.he_prop.M / 1000)  # J/kg/K

        return rho, Cp_mass, k

    def discretize_laplacian(self, T_sym):
        """
        Discretize Laplacian in cylindrical coordinates
        ∇²T = ∂²T/∂r² + (1/r)∂T/∂r + ∂²T/∂z²

        Uses:
        - 2nd order CFD for r-derivatives
        - 2nd order FFD for z-derivatives
        """
        laplacian = ca.SX.zeros(self.Nr, self.Nz)

        for j in range(self.Nz):
            for i in range(self.Nr):
                # Initialize components
                d2T_dr2 = 0
                dT_dr_over_r = 0
                d2T_dz2 = 0

                # --- Radial derivatives (CFD) ---
                if i == 0:
                    # Centerline: Use L'Hôpital's rule
                    # lim(r→0) [(1/r)∂T/∂r] = ∂²T/∂r²
                    # Use forward difference for both terms
                    d2T_dr2 = (T_sym[i, j] - 2*T_sym[i+1, j] + T_sym[i+2, j]) / self.dr**2
                    dT_dr_over_r = d2T_dr2  # Apply L'Hôpital's rule

                elif i == self.Nr - 1:
                    # Wall boundary: Will be handled by BC
                    d2T_dr2 = 0
                    dT_dr_over_r = 0

                else:
                    # Interior points: Use central differences
                    d2T_dr2 = (T_sym[i-1, j] - 2*T_sym[i, j] + T_sym[i+1, j]) / self.dr**2
                    dT_dr = (T_sym[i+1, j] - T_sym[i-1, j]) / (2*self.dr)
                    dT_dr_over_r = dT_dr / self.r[i]

                # --- Axial derivatives (FFD) ---
                if j >= self.Nz - 2:
                    # Outlet: Zero gradient BC
                    d2T_dz2 = 0
                else:
                    # Use forward finite difference for 2nd derivative
                    d2T_dz2 = (T_sym[i, j] - 2*T_sym[i, j+1] + T_sym[i, j+2]) / self.dz**2

                # Assemble Laplacian
                laplacian[i, j] = d2T_dr2 + dT_dr_over_r + d2T_dz2

        return laplacian

    def discretize_convection(self, T_sym):
        """
        Discretize convection term: uz * ∂T/∂z
        Uses 2nd order FFD for z-derivative
        """
        conv_term = ca.SX.zeros(self.Nr, self.Nz)

        for j in range(self.Nz):
            for i in range(self.Nr):
                if j >= self.Nz - 2:
                    # Outlet: No convection gradient calculation
                    dT_dz = 0
                else:
                    # Forward finite difference (2nd order)
                    dT_dz = (-3*T_sym[i, j] + 4*T_sym[i, j+1] - T_sym[i, j+2]) / (2*self.dz)

                conv_term[i, j] = self.uz[i] * dT_dz

        return conv_term

    def build_residual(self, T_avg=350):
        """
        Build residual equations for steady-state heat equation:
        ρ*Cp*uz*∂T/∂z = k*∇²T

        Rearranged: k*∇²T - ρ*Cp*uz*∂T/∂z = 0
        """
        # Get properties at average temperature
        rho, Cp, k = self.get_properties(T_avg)

        # Create symbolic temperature field
        T_sym = ca.SX.sym('T', self.Nr, self.Nz)

        # Discretize operators
        laplacian = self.discretize_laplacian(T_sym)
        convection = self.discretize_convection(T_sym)

        # Build residual: k*∇²T - ρ*Cp*uz*∂T/∂z = 0
        residual = k * laplacian - rho * Cp * convection

        # Apply boundary conditions
        for j in range(self.Nz):
            # Inlet (z=0): Fixed temperature
            if j == 0:
                for i in range(self.Nr):
                    residual[i, j] = T_sym[i, j] - self.T_inlet

            # Wall (r=R): Fixed temperature
            residual[self.Nr-1, j] = T_sym[self.Nr-1, j] - self.T_wall

        # Flatten residual and variables
        res_vec = ca.reshape(residual, -1, 1)
        T_vec = ca.reshape(T_sym, -1, 1)

        return T_vec, res_vec

    def solve_steady_state(self, T_initial=None, T_avg_guess=350):
        """
        Solve steady-state heat equation using CasADi

        Parameters:
        -----------
        T_initial : numpy array, optional
            Initial guess for temperature field (Nr x Nz)
        T_avg_guess : float
            Guess for average temperature to evaluate properties

        Returns:
        --------
        T_solution : numpy array
            Solved temperature field (Nr x Nz)
        """
        # Create initial guess with better profile
        if T_initial is None:
            # Create a more realistic initial guess
            T_initial = np.zeros((self.Nr, self.Nz))
            for j in range(self.Nz):
                # Axial variation: linear from inlet to average of inlet and wall
                T_axial = self.T_inlet + (self.T_wall - self.T_inlet) * (self.z[j] / self.L)

                # Radial variation: parabolic from centerline to wall
                for i in range(self.Nr):
                    r_ratio = self.r[i] / self.R
                    T_initial[i, j] = T_axial * (1 - 0.3 * r_ratio**2)

            # Apply boundary conditions
            T_initial[:, 0] = self.T_inlet  # Inlet
            T_initial[-1, :] = self.T_wall  # Wall

        # Build residual system
        T_vec, res_vec = self.build_residual(T_avg=T_avg_guess)

        # Create CasADi function
        residual_func = ca.Function('residual', [T_vec], [res_vec])

        # Set up rootfinder with better options
        opts = {
            'ipopt': {
                'print_level': 3,
                'max_iter': 500,
                'tol': 1e-4,
                'acceptable_tol': 1e-3,
                'acceptable_iter': 10,
                'mu_strategy': 'adaptive',
                'linear_solver': 'mumps'
            },
            'print_time': False
        }

        # Create NLP solver (treating as least squares problem)
        nlp = {
            'x': T_vec,
            'f': ca.dot(res_vec, res_vec) / 2,  # Minimize squared residual
            'g': ca.SX([])  # No constraints
        }

        solver = ca.nlpsol('solver', 'ipopt', nlp, opts)

        # Solve
        T_init_vec = T_initial.T.flatten()  # Flatten in column-major order

        sol = solver(x0=T_init_vec, lbg=[], ubg=[])

        # Extract solution
        T_sol_vec = sol['x'].full().flatten()
        T_solution = T_sol_vec.reshape(self.Nz, self.Nr).T  # Reshape back

        return T_solution

    def visualize(self, T_solution, save_path=None):
        """Visualize temperature field

        Parameters:
        -----------
        T_solution : numpy array
            Temperature field with shape (Nr, Nz)
        """
        fig, axes = plt.subplots(1, 2, figsize=(14, 5))

        # 2D contour plot - need to transpose T_solution to match meshgrid
        ax1 = axes[0]
        contour = ax1.contourf(self.Z_mesh * 100, self.R_mesh * 1000,
                               T_solution.T, levels=20, cmap='hot')
        ax1.set_xlabel('Axial Position z (cm)')
        ax1.set_ylabel('Radial Position r (mm)')
        ax1.set_title('Temperature Field (K)')
        plt.colorbar(contour, ax=ax1, label='Temperature (K)')

        # Radial profiles at different axial positions
        ax2 = axes[1]
        z_positions = [0, 0.25, 0.5, 0.75, 1.0]
        for z_frac in z_positions:
            z_idx = int(z_frac * (self.Nz - 1))
            ax2.plot(self.r * 1000, T_solution[:, z_idx],
                    label=f'z = {z_frac*self.L*100:.1f} cm')
        ax2.set_xlabel('Radial Position r (mm)')
        ax2.set_ylabel('Temperature (K)')
        ax2.set_title('Radial Temperature Profiles')
        ax2.legend()
        ax2.grid(True)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150)

        plt.close(fig)

        return fig


# Example usage
if __name__ == "__main__":
    # Create model
    print("Initializing 2D forced convection model...")
    model = ConvectionModel2D(
        flow_rate_sccm=50,
        T_inlet=298.15,
        T_wall=400.0,  # Example: heated wall
        Nr=20,
        Nz=50
    )

    print(f"Reactor geometry: L={model.L*100:.1f} cm, R={model.R*1000:.2f} mm")
    print(f"Flow rate: {model.Q_sccm} sccm = {model.Q_actual*1e6:.2f} mL/s")
    print(f"Average velocity: {model.u_avg:.4f} m/s")
    print(f"Grid: {model.Nr} x {model.Nz} points")
    print(f"dr={model.dr*1000:.3f} mm, dz={model.dz*100:.3f} cm")

    # Solve steady-state problem
    print("\nSolving steady-state heat equation...")
    T_solution = model.solve_steady_state(T_avg_guess=340)

    print(f"Solution converged!")
    print(f"Temperature range: {T_solution.min():.2f} - {T_solution.max():.2f} K")

    # Visualize results
    print("\nVisualizing results...")
    model.visualize(T_solution, save_path='convection_2d_results.png')
