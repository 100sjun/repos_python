"""
Physical properties module for heat balance calculation.
- Helium gas: temperature-dependent properties
- KanthalD (FeCrAl alloy): heating element properties
"""

# Libraries
import torch
import json
import os

ArrayLike = float | torch.Tensor
BASE_DIR = os.path.dirname(__file__)
DATA_PATH = os.path.join(BASE_DIR, "..", "data", "property.json")

# ====================================================================
# Helium Gas Properties (Temperature-Dependent)
# ====================================================================

class Prop_He:
    """
    Helium gas properties as functions of temperature.
    Reference: Nist Chemistry WebBook
    """
    Mw = 4.0026

    @staticmethod
    def rho(T: torch.Tensor) -> torch.Tensor:
        """
        Mass density from property data [kg/m3]

        Args:
            T: Temperature [C]
        Returns:
            rho: Mass density [kg/m3]
        """
        a = json.load(open(DATA_PATH))['helium']['rho']['a']
        b = json.load(open(DATA_PATH))['helium']['rho']['b']
        function = json.load(open(DATA_PATH))['helium']['rho']['function']

        if function == 'inverse':
            return a / (T + b) * 1000 # g/ml to kg/m3
        elif function == 'linear':
            return a + b * T
        elif function == 'log':
            return a + b * torch.log(T)
        else:
            return a

    @staticmethod
    def cp(T: torch.Tensor) -> torch.Tensor:
        """
        Specific heat capacity at constant pressure [J/mol*K]
        
        Args:
            T: Temperature [C]
        Returns:
            cp: specific heat [J/mol*K]
        """
        a = json.load(open(DATA_PATH))['helium']['Cp']['a']
        b = json.load(open(DATA_PATH))['helium']['Cp']['b']
        function = json.load(open(DATA_PATH))['helium']['Cp']['function']

        if function == 'constant':
            return a
        elif function == 'linear':
            return a + b * T
        elif function == 'log':
            return a + b * torch.log(T)
        else:
            return a
    
    @staticmethod
    def mu(T: torch.Tensor) -> torch.Tensor:
        """
        Dynamic viscosity [Pa*s]
        
        Args:
            T: Temperature [C]
        Returns:
            mu: Dynamic viscosity [Pa*s]
        """
        a = json.load(open(DATA_PATH))['helium']['mu']['a']
        b = json.load(open(DATA_PATH))['helium']['mu']['b']
        function = json.load(open(DATA_PATH))['helium']['mu']['function']

        if function == 'constant':
            return a
        elif function == 'linear':
            return a + b * T
        elif function == 'log':
            return a + b * torch.log(T)
        else:
            return a
    
    @staticmethod
    def k(T: torch.Tensor) -> torch.Tensor:
        """
        Thermal conductivity [W/m/K]
        
        Args:
            T: Temperature [C]
        Returns:
            k: Thermal conductivity [W/m/K]
        """
        a = json.load(open(DATA_PATH))['helium']['k']['a']
        b = json.load(open(DATA_PATH))['helium']['k']['b']
        function = json.load(open(DATA_PATH))['helium']['k']['function']

        if function == 'constant':
            return a
        elif function == 'linear':
            return a + b * T
        elif function == 'log':
            return a + b * torch.log(T)
        else:
            return a
    
    @staticmethod
    def prandtl(T: torch.Tensor) -> torch.Tensor:
        """Prandtl number [-]"""
        cp = Prop_He.cp(T) / (Prop_He.Mw/1000)
        mu = Prop_He.mu(T)
        k = Prop_He.k(T)
        return cp * mu / k

    @classmethod
    def all_properties(cls, T: torch.Tensor) -> dict:
        """Return all properties as dictionary"""
        return {
            'rho': cls.rho(T),
            'cp': cls.cp(T),
            'mu': cls.mu(T),
            'k': cls.k(T),
            'Pr': cls.prandtl(T)

        }

# ====================================================================
# KanthalD Properties (Temperature-Dependent)
# ====================================================================

class Prop_kanthal:
    """
    Kanthal D heating element properties as functions of temperature

    Reference: Kanthal D datasheet
    """

    @staticmethod
    def rho(T: torch.Tensor) -> torch.Tensor:
        """
        Mass density from property data [kg/m3]

        Args:
            T: Temperature [C]
        Returns:
            rho: Mass density [kg/m3]
        """
        a = json.load(open(DATA_PATH))['kanthal']['rho']['a']
        b = json.load(open(DATA_PATH))['kanthal']['rho']['b']
        function = json.load(open(DATA_PATH))['kanthal']['rho']['function']

        if function == 'inverse':
            return a / (T + b)
        elif function == 'linear':
            return a + b * T
        elif function == 'log':
            return a + b * torch.log(T)
        else:
            return a
    
    @staticmethod
    def cp(T: torch.Tensor) -> torch.Tensor:
        """
        Specific heat capacity at constant pressure [kJ/kg/K]

        Args:
            T: Temperature [C]
        Returns:
            cp: Specific heat capacity at constant pressure [kJ/kg/K]
        """
        a = json.load(open(DATA_PATH))['kanthal']['Cp']['a']
        b = json.load(open(DATA_PATH))['kanthal']['Cp']['b']
        function = json.load(open(DATA_PATH))['kanthal']['Cp']['function']

        if function == 'inverse':
            return a / (T + b)
        elif function == 'linear':
            return a + b * T
        elif function == 'log':
            # Clamp temperature to avoid log(0) or log(negative)
            T_safe = torch.clamp(T, min=1.0)  # Minimum 1°C for log safety
            return a + b * torch.log(T_safe)
        else:
            return a

    @staticmethod
    def k(T: torch.Tensor) -> torch.Tensor:
        """
        Thermal conductivity [W/m/K]
        
        Args:
            T: Temperature [C]
        Returns:
            k: Thermal conductivity [W/m/K]
        """
        a = json.load(open(DATA_PATH))['kanthal']['k']['a']
        b = json.load(open(DATA_PATH))['kanthal']['k']['b']
        function = json.load(open(DATA_PATH))['kanthal']['k']['function']

        if function == 'inverse':
            return a / (T + b)
        elif function == 'linear':
            return a + b * T
        elif function == 'log':
            return a + b * torch.log(T)
        else:
            return a
    
    @staticmethod
    def er(T: torch.Tensor) -> torch.Tensor:
        """
        Resistivity [ohm*m2/m]
        
        Args:
            T: Temperature [C]
        Returns:
            er: Resistivity [ohm*m2/m]
        """
        a = json.load(open(DATA_PATH))['kanthal']['er']['a']
        b = json.load(open(DATA_PATH))['kanthal']['er']['b']
        function = json.load(open(DATA_PATH))['kanthal']['er']['function']

        if function == 'inverse':
            return a / (T + b)
        elif function == 'linear':
            return a + b * T
        elif function == 'log':
            return a + b * torch.log(T)
        else:
            return a
    
    @classmethod
    def all_properties(cls, T: torch.Tensor) -> dict:
        """Return all properties as dictionary"""
        return {
            'rho': cls.rho(T),
            'cp': cls.cp(T),
            'k': cls.k(T),
            'er': cls.er(T)
        }

# ====================================================================
# Heat Transfer Correlations
# ====================================================================

def Nu_laminar(Re: torch.Tensor, Pr: torch.Tensor, L_D: float) -> torch.Tensor:
    """
    Nusselt number for laminar flow in circular tube
    """
    # Sieder-Tate for developing region
    Nu_developing = 1.86 * Re * Pr / L_D

    # Fully developed (constant heat flux)
    Nu_developed = 4.36

    # Take maximum (developing region has higher Nu)
    return torch.maximum(torch.tensor(Nu_developing), torch.tensor(Nu_developed))

def h_coeff(T_gas: torch.Tensor, T_wall: torch.Tensor, D_inner: float, mass_flow: float, L: float = None) -> torch.Tensor:
    """
    Calculate convective heat transfer coefficient [W/m2*K]

    Args:
        T_gas: Gas temperature [C]
        T_wall: Wall temperature [C]
        D_inner: Inner diameter [m]
        mass_flow: Mass flow rate [kg/s]
        L: Tube length [m], optional for developing flow

    Returns:
        h: Heat transfer coefficient [W/m2*K]
    """
    # Film temperature for property evaluation
    T_film = 0.5 * (T_gas + T_wall)

    # Get properties at film temperature
    props_film = Prop_He.all_properties(T_film)

    # Cross-sectional area
    A_cs = torch.pi * (D_inner**2) / 4

    # Velocity and Reynolds number
    velocity = mass_flow / (props_film['rho'] * A_cs)
    Re = props_film['rho'] * velocity * D_inner / props_film['mu']
    Pr = props_film['Pr']

    # Nusselt number
    L_D = L / D_inner
    Nu = Nu_laminar(torch.tensor(Re), torch.tensor(Pr), L_D)

    # Heat transfer coefficient
    h = Nu * props_film['k'] / D_inner

    return h

# ====================================================================
# Utility Functions
# ====================================================================

def sccm_to_kg_per_s(sccm: float, M: float = 4.0026) -> float:
    """
    Convert sccm to kg/s

    Args:
        sccm: Flow rate in standard mL/min (at 0 C, 1atm)
        M: Molecular weight [g/mol]
    
    Returns:
        Mass flow rate [kg/s]
    """

    mL_per_s = sccm / 60.0 # mL/s
    m3_per_s = mL_per_s * 1e-6 # m3/s

    # Molar volume at STP
    V_m = 8.314 * 273.15 / 101325 # m3/mol = 0.02241 m3/mol
    mol_per_s = m3_per_s / V_m
    kg_per_s = mol_per_s * M / 1000 # kg/s

    return kg_per_s

# ====================================================================
# Test / Example Usage
# ====================================================================

if __name__ == "__main__":
    print("=" * 60)
    print("Physical Properties Test")
    print("=" * 60)
    
    # Test temperatures
    T_test = torch.tensor([300, 500, 700, 900, 1100])  # K
    
    print("\n--- Helium Properties ---")
    print(f"{'T [C]':>8} {'rho [kg/m³]':>12} {'cp [J/kgK]':>12} {'μ [Pa·s]':>12} {'k [W/mK]':>10} {'Pr':>8}")
    print("-" * 70)
    for T in T_test:
        props = Prop_He.all_properties(T)
        print(f"{T:>8.0f} {props['rho']:>12.4f} {props['cp']:>12.3f} {props['mu']:>12.2e} {props['k']:>10.4f} {props['Pr']:>8.4f}")
    
    print("\n--- KanthalD Properties ---")
    print(f"{'T [C]':>8} {'rho [kg/m³]':>12} {'cp [J/kgK]':>12} {'k [W/mK]':>10} {'er [Ω·m2/m]':>12}")
    print("-" * 60)
    for T in T_test:
        props = Prop_kanthal.all_properties(T)
        print(f"{T:>8.0f} {props['rho']:>12.1f} {props['cp']:>12.3f} {props['k']:>10.2f} {props['er']:>12.2e}")
    
    print("\n--- Flow Rate Conversion ---")
    sccm = 50.0
    mdot = sccm_to_kg_per_s(sccm)
    print(f"Flow rate: {sccm} sccm He = {mdot:.6e} kg/s = {mdot*1e6:.4f} mg/s")
    
    print("\n--- Heat Transfer Coefficient ---")
    D_i = 5.03e-3  # m
    L = 0.4312  # m
    T_gas = 500  # K
    T_wall = 800  # K
    
    h = h_coeff(T_gas, T_wall, D_i, mdot, L)
    print(f"At T_gas={T_gas}C, T_wall={T_wall}C:")
    print(f"  h = {h:.2f} W/(m²·K)")
    
    # Reynolds number check
    props = Prop_He.all_properties((T_gas + T_wall)/2)
    A_cs = torch.pi * D_i**2 / 4
    v = mdot / (props['rho'] * A_cs)
    Re = props['rho'] * v * D_i / props['mu']
    print(f"  Re = {Re:.1f} ({'laminar' if Re < 2300 else 'turbulent'})")
    print(f"  velocity = {v:.3f} m/s")

