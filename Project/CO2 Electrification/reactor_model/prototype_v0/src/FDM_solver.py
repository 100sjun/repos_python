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
from Geometry import ReactorGeometry, Mesh1D, OperatingConditions
from Property import Prop_He, Prop_kanthal

# ====================================================================
# FDM Solver Class
# ====================================================================
class ThermalFDMSolver(nn.Module):
    '''
    1D Transient FDM solver for coupled gas-wall heat transfer

    Solves:
    - Gas: quasi-steady convection equation (marching from inlet)
    - Wall: transient energy balance with generation, convection, and loss

    Args:
        geometry: ReactorGeometry instance
        mesh: Mesh1D instance
        operating_conditions: OperatingConditions instance
        device: torch device ('cpu' or 'cuda')
    '''

    def __init__(self,
                 geometry: ReactorGeometry,
                 mesh: Mesh1D,
                 operating_conditions: OperatingConditions,
                 device: str = 'cpu'):
        super(ThermalFDMSolver, self).__init__()

        # Store configuration
        self.geom = geometry
        self.mesh = mesh
        self.ops = operating_conditions
        self.device = torch.device(device)

        # Mesh parameters
        self.n_nodes = mesh.n_nodes
        self.dz = mesh.dz_actual

        # Initialize temperature fields as parameters for autograd
        self.register_buffer('Tg', torch.ones(self.n_nodes, device=self.device) * operating_conditions.Tin)
        self.register_buffer('Tw', torch.ones(self.n_nodes, device=self.device) * operating_conditions.Tin)

        # Time tracking
        self.t = 0.0
        self.dt = operating_conditions.dt

        # History storage
        self.history = {
            'time': [],
            'Tg': [],
            'Tw': []
        }

    def compute_heat_transfer_coefficient(self, Tg: torch.Tensor, Tw: torch.Tensor) -> torch.Tensor:
        """
        Calculate convective heat transfer coefficient [W/m2*K]

        Args:
            Tg: Gas temperature [C] (n_nodes,)
            Tw: Wall temperature [C] (n_nodes,)

        Returns:
            h: Heat transfer coefficient [W/m2*K] (n_nodes,)
        """
        # Film temperature
        T_film = 0.5 * (Tg + Tw)

        # Properties at film temperature
        rho = Prop_He.rho(T_film)
        cp_molar = Prop_He.cp(T_film)
        cp_mass = cp_molar / (Prop_He.Mw / 1000)  # J/kg/K
        mu = Prop_He.mu(T_film)
        k = Prop_He.k(T_film)

        # Prandtl number
        Pr = cp_mass * mu / k

        # Cross-sectional area and velocity
        A_cs = torch.pi * (self.geom.id ** 2) / 4
        velocity = self.ops.Fm / (rho * A_cs)

        # Reynolds number
        Re = rho * velocity * self.geom.id / mu

        # Nusselt number (laminar flow)
        L_D = self.geom.L / self.geom.id
        Nu_developing = 1.86 * Re * Pr / L_D
        Nu_developed = 4.36
        Nu = torch.maximum(Nu_developing, torch.ones_like(Nu_developing) * Nu_developed)

        # Heat transfer coefficient
        h = Nu * k / self.geom.id

        return h

    def solve_gas_temperature(self, Tw: torch.Tensor) -> torch.Tensor:
        """
        Solve quasi-steady gas energy equation by marching from inlet

        Equation: mdot * Cp * dTg/dz = h * pi * Di * (Tw - Tg)

        Args:
            Tw: Wall temperature [C] (n_nodes,)

        Returns:
            Tg: Gas temperature [C] (n_nodes,)
        """
        Tg = torch.zeros_like(Tw)
        Tg[0] = torch.clamp(torch.tensor(self.ops.Tin), min=1.0, max=1500.0)  # Inlet boundary condition

        # March downstream using upwind scheme
        for i in range(1, self.n_nodes):
            # Clamp Tg at previous node to ensure valid property evaluation
            Tg[i-1] = torch.clamp(Tg[i-1], min=1.0, max=1500.0)

            # Use properties at previous node for stability
            h = self.compute_heat_transfer_coefficient(Tg[i-1:i], Tw[i-1:i])

            # Gas properties at previous node
            cp_molar = Prop_He.cp(Tg[i-1])
            cp_mass = cp_molar / (Prop_He.Mw / 1000)  # J/kg/K

            # Source term: heat transfer from wall
            Q_conv = h[0] * torch.pi * self.geom.id * (Tw[i-1] - Tg[i-1])

            # Temperature gradient
            dTg_dz = Q_conv / (self.ops.Fm * cp_mass)

            # Check for numerical issues
            if torch.isnan(dTg_dz) or torch.isinf(dTg_dz):
                print(f"🔴 NaN/Inf in dTg_dz at node {i}!")
                print(f"   h={h[0]:.2e}, Q_conv={Q_conv:.2e}, cp_mass={cp_mass:.2e}")
                print(f"   Tg[{i-1}]={Tg[i-1]:.2f}, Tw[{i-1}]={Tw[i-1]:.2f}")

            # Update temperature
            Tg[i] = Tg[i-1] + dTg_dz * self.dz

            # Clamp immediately after update
            Tg[i] = torch.clamp(Tg[i], min=1.0, max=1500.0)

        return Tg

    def compute_wall_rhs(self, Tg: torch.Tensor, Tw: torch.Tensor) -> torch.Tensor:
        """
        Compute RHS of wall energy equation

        Equation: rho * A * Cp * dTw/dt = Qgen - h*pi*Di*(Tw - Tg) - Qloss

        Args:
            Tg: Gas temperature [C] (n_nodes,)
            Tw: Wall temperature [C] (n_nodes,)

        Returns:
            dTw_dt: Time derivative of wall temperature [C/s] (n_nodes,)
        """
        # DEBUG: Check for negative temperatures before property calculation
        if (Tw <= 0).any():
            print(f"⚠️  WARNING: Tw has non-positive values! Tw_min={Tw.min():.2f}°C at t={self.t:.2f}s")
            print(f"   This will cause log(Tw) = NaN in Prop_kanthal.cp()")

        # Wall properties
        rho_w = Prop_kanthal.rho(Tw)
        cp_w = Prop_kanthal.cp(Tw) * 1000  # Convert kJ/kg/K to J/kg/K

        # DEBUG: Check for NaN after property calculation
        if torch.isnan(cp_w).any():
            print(f"🔴 NaN detected in cp_w!")
            print(f"   Tw range: [{Tw.min():.2f}, {Tw.max():.2f}]")
            print(f"   rho_w range: [{rho_w.min():.2f}, {rho_w.max():.2f}]")

        # Heat transfer coefficient
        h = self.compute_heat_transfer_coefficient(Tg, Tw)

        # Heat generation (uniform along tube)
        q_gen = self.ops.Pw / self.geom.Vw  # W/m3
        Q_gen = q_gen * self.geom.Aw * self.dz  # W per cell

        # Convective heat loss to gas
        Q_conv = h * torch.pi * self.geom.id * self.dz * (Tw - Tg)

        # Radiation/natural convection loss (simplified model)
        # Assuming ambient cooling proportional to temperature difference
        h_amb = 10.0  # W/m2/K (natural convection + radiation)
        Q_loss = h_amb * torch.pi * self.geom.od * self.dz * (Tw - self.ops.Tamb)

        # Net energy balance
        Q_net = Q_gen - Q_conv - Q_loss

        # Mass per cell
        m_cell = rho_w * self.geom.Aw * self.dz

        # Temperature rate of change
        dTw_dt = Q_net / (m_cell * cp_w)

        # DEBUG: Check result
        if torch.isnan(dTw_dt).any() or torch.isinf(dTw_dt).any():
            print(f"🔴 NaN/Inf in dTw_dt!")
            print(f"   Q_gen={Q_gen:.2e}, Q_conv_mean={Q_conv.mean():.2e}, Q_loss_mean={Q_loss.mean():.2e}")
            print(f"   dTw_dt range: [{dTw_dt.min():.2e}, {dTw_dt.max():.2e}]")

        return dTw_dt

    def step(self, verbose: bool = False) -> tuple[torch.Tensor, torch.Tensor]:
        """
        Advance solution by one time step using Forward Euler

        Args:
            verbose: Print detailed debugging info

        Returns:
            Tg: Gas temperature [C] (n_nodes,)
            Tw: Wall temperature [C] (n_nodes,)
        """
        if verbose:
            print(f"\n🔍 Step start (t={self.t:.3f}s):")
            print(f"   Tw: [{self.Tw.min():.2f}, {self.Tw.max():.2f}]°C")
            print(f"   Tg: [{self.Tg.min():.2f}, {self.Tg.max():.2f}]°C")

        # Apply physical constraints BEFORE computations to ensure valid property evaluations
        # This prevents log(negative) in Prop_kanthal.cp()
        self.Tw = torch.clamp(self.Tw, min=1.0, max=1500.0)  # Min 1°C to prevent log(0)
        self.Tg = torch.clamp(self.Tg, min=1.0, max=1500.0)

        # Solve gas temperature (quasi-steady)
        self.Tg = self.solve_gas_temperature(self.Tw)

        if verbose:
            print(f"   After solve_gas: Tg=[{self.Tg.min():.2f}, {self.Tg.max():.2f}]°C")
            if torch.isnan(self.Tg).any():
                print(f"   🔴 NaN in Tg after solve_gas_temperature!")

        # Compute wall temperature rate of change
        dTw_dt = self.compute_wall_rhs(self.Tg, self.Tw)

        if verbose:
            print(f"   dTw_dt: [{dTw_dt.min():.2e}, {dTw_dt.max():.2e}] K/s")
            if torch.isnan(dTw_dt).any():
                print(f"   🔴 NaN in dTw_dt after compute_wall_rhs!")

        # Update wall temperature (Forward Euler)
        Tw_old = self.Tw.clone()
        self.Tw = self.Tw + dTw_dt * self.dt

        if verbose:
            delta_Tw = (self.Tw - Tw_old)
            print(f"   ΔTw: [{delta_Tw.min():.2f}, {delta_Tw.max():.2f}]°C")
            print(f"   Tw after update: [{self.Tw.min():.2f}, {self.Tw.max():.2f}]°C")
            if torch.isnan(self.Tw).any():
                print(f"   🔴 NaN in Tw after update!")

        # Apply physical constraints again after update
        self.Tw = torch.clamp(self.Tw, min=1.0, max=1500.0)
        self.Tg = torch.clamp(self.Tg, min=1.0, max=1500.0)

        # Check for NaN or Inf
        if torch.isnan(self.Tw).any() or torch.isnan(self.Tg).any():
            print(f"\n🔴 NaN DETECTED at t={self.t:.3f}s!")
            print(f"   Tw: {self.Tw}")
            print(f"   Tg: {self.Tg}")
            print(f"   dTw_dt: {dTw_dt}")
            raise ValueError(f"NaN detected at t={self.t:.3f}s. "
                           f"Tw range: [{self.Tw.min():.1f}, {self.Tw.max():.1f}], "
                           f"Tg range: [{self.Tg.min():.1f}, {self.Tg.max():.1f}]")

        if torch.isinf(self.Tw).any() or torch.isinf(self.Tg).any():
            raise ValueError(f"Inf detected at t={self.t:.3f}s. Reduce time step dt.")

        # Update time
        self.t += self.dt

        return self.Tg.clone(), self.Tw.clone()

    def solve(self, save_interval: int = 10) -> dict:
        """
        Run transient simulation until end time

        Args:
            save_interval: Save history every N steps

        Returns:
            history: Dictionary with time, Tg, Tw arrays
        """
        n_steps = int(self.ops.tend / self.dt)

        print(f"Starting simulation: {n_steps} steps, dt={self.dt}s")
        print(f"Power: {self.ops.Pw}W, Flow: {self.ops.Fv_sccm} sccm")

        for step in range(n_steps):
            # Verbose output for first 60 steps (0.6s) to catch the problem
            verbose = (step < 60) or (step % 10 == 0 and step < 100)
            Tg, Tw = self.step(verbose=verbose)

            # Save history
            if step % save_interval == 0:
                self.history['time'].append(self.t)
                self.history['Tg'].append(Tg.detach().cpu().numpy())
                self.history['Tw'].append(Tw.detach().cpu().numpy())

                if step % 100 == 0:
                    print(f"Step {step}/{n_steps}, t={self.t:.1f}s, "
                          f"Tg_out={Tg[-1]:.1f}C, Tw_max={Tw.max():.1f}C, Tw_min={Tw.min():.1f}C")

        # Convert lists to arrays
        self.history['time'] = torch.tensor(self.history['time'])
        self.history['Tg'] = torch.tensor(self.history['Tg'])
        self.history['Tw'] = torch.tensor(self.history['Tw'])

        print("Simulation completed!")
        return self.history

    def get_tc_temperatures(self) -> tuple[torch.Tensor, torch.Tensor]:
        """
        Extract temperatures at thermocouple locations

        Returns:
            tc_times: Time array [s]
            tc_temps: Temperature array [C] (n_times, n_tc)
        """
        tc_temps = self.history['Tw'][:, self.mesh.tc_indices]
        return self.history['time'], tc_temps

# ====================================================================
# Test / Example Usage
# ====================================================================
if __name__ == "__main__":
    from Geometry import ReactorGeometry, Mesh1D, OperatingConditions
    import matplotlib.pyplot as plt

    print("=" * 60)
    print("1D FDM Thermal Solver - Test Case")
    print("=" * 60)

    # Create geometry
    geom = ReactorGeometry()
    print(geom.summary())

    # Create mesh (1mm spacing)
    mesh = Mesh1D(geom, dz=1e-3)
    print(mesh.summary())

    # Operating conditions: 1.4W, 50 sccm He
    ops = OperatingConditions(
        Fv_sccm=50.0,
        Tin=25.0,
        P=1.01325,
        Pw=0.01,  # 1.4W as specified
        Tamb=25.0,
        tend=60.0,  # 1 minute for quick test
        dt=0.1  # Reduced from 0.1 to 0.01 for numerical stability
    )
    print(ops.summary())

    # Create solver
    solver = ThermalFDMSolver(geom, mesh, ops, device='cpu')

    # Run simulation
    history = solver.solve(save_interval=10)

    # Extract TC temperatures
    tc_times, tc_temps = solver.get_tc_temperatures()

    # Plot results
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 1. Temperature profiles at different times
    ax = axes[0, 0]
    time_indices = [0, len(history['time'])//4, len(history['time'])//2, -1]
    for idx in time_indices:
        t = history['time'][idx]
        ax.plot(mesh.z_mm, history['Tw'][idx], label=f't={t:.1f}s')
    ax.set_xlabel('Axial position [mm]')
    ax.set_ylabel('Wall temperature [°C]')
    ax.set_title('Wall Temperature Profiles')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 2. TC temperature evolution
    ax = axes[0, 1]
    for i in range(geom.n_tc):
        ax.plot(tc_times, tc_temps[:, i], label=f'TC{i+1}')
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Temperature [°C]')
    ax.set_title('Thermocouple Temperature Evolution')
    ax.legend(ncol=2, fontsize=8)
    ax.grid(True, alpha=0.3)

    # 3. Gas vs Wall temperature at final time
    ax = axes[1, 0]
    ax.plot(mesh.z_mm, history['Tg'][-1], 'b-', label='Gas')
    ax.plot(mesh.z_mm, history['Tw'][-1], 'r-', label='Wall')
    ax.set_xlabel('Axial position [mm]')
    ax.set_ylabel('Temperature [°C]')
    ax.set_title(f'Final Temperature Profile (t={history["time"][-1]:.1f}s)')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 4. Temperature difference
    ax = axes[1, 1]
    dT = history['Tw'][-1] - history['Tg'][-1]
    ax.plot(mesh.z_mm, dT, 'g-')
    ax.set_xlabel('Axial position [mm]')
    ax.set_ylabel('Temperature difference [°C]')
    ax.set_title('Wall-Gas Temperature Difference')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('fdm_results.png', dpi=150)
    print("\nResults saved to 'fdm_results.png'")
    plt.show()

    # Summary statistics
    print("\n" + "=" * 60)
    print("Simulation Summary")
    print("=" * 60)
    print(f"Final gas outlet temperature: {history['Tg'][-1, -1]:.1f} °C")
    print(f"Maximum wall temperature: {history['Tw'][-1].max():.1f} °C")
    print(f"Wall-gas temperature difference: {dT.mean():.2f} ± {dT.std():.2f} °C")
    print(f"Simulation time: {history['time'][-1]:.1f} s")
    print("=" * 60)