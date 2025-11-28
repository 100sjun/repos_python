'''
Geometry module for tubular reactor heat balance calculation.

Reactor configurations:
- Kanthal D single tube reactor
- Inner diameter: 5.03 mm
- Outer diameter: 6.33 mm
- Total length: 430 mm
'''

# Libraries
import numpy as np
import matplotlib.pyplot as plt

# ====================================================================
# Reactor Geometry Class
# ====================================================================

class ReactorGeometry:
    def __init__(self, inner_diameter = 5.03e-3, outer_diameter = 6.33e-3, length = 430e-3):
        self.id = inner_diameter
        self.od = outer_diameter
        self.L = length
        self.tc_pos = np.array([35e-3, 48e-3, 64e-3, 83e-3, 105e-3, 130e-3, 158e-3, 189.6e-3, 241e-3, 272e-3, 300e-3,
        325e-3, 347e-3, 366e-3, 382e-3, 395e-3])

        self._calculate_derived()

    def _calculate_derived(self):
        '''Calculate derived geometric quuantities'''
        # Radii
        self.ri = self.id / 2
        self.ro = self.od / 2

        # Wall thickness
        self.tkn = self.ro - self.ri

        # Cross-sectional areas
        self.Ai = np.pi * self.ri**2
        self.Ao = np.pi * self.ro**2
        self.Aw = self.Ao - self.Ai

        # Perimeters
        self.pri = np.pi * self.id
        self.pro = np.pi * self.od

        # Surface ares (total tube)
        self.Si = self.pri * self.L
        self.So = self.pro * self.L

        # Volume
        self.Vi = self.Ai * self.L
        self.Vw = self.Aw * self.L

    @property
    def n_tc(self) -> int:
        '''Number of thermocouples'''
        return len(self.tc_pos)

    @property
    def tc_pos_mm(self) -> np.ndarray:
        '''TC positions in mm for convenience'''
        return self.tc_pos * 1000

    def summary(self) -> str:
        '''Return formatted summary of geometry'''
        lines =[
            '=' * 50,
            'Reactor Geometry Summary',
            '=' * 50,
            f'Inner diameter:   {self.id*1000:.2f} mm',
            f'Outer diameter:   {self.od*1000:.2f} mm',
            f'Wall thickness:   {self.tkn*1000:.2f} mm',
            f'Length:           {self.L*1000:.1f} mm',
            '',
            f'Inner area:       {self.Ai*1e6:.2f} mm²',
            f'Wall area:        {self.Aw*1e6:.2f} mm²',
            f'Inner surface:    {self.Si*1e6:.2f} mm²',
            f'Outer surface:    {self.So*1e6:.2f} mm²',
            f'Wall Volume:      {self.Vw*1e6:.2f} mm³',
            '',
            f'Thermocouples:    {self.n_tc}',
            f'Tc positions:     {self.tc_pos_mm[0]:.1f} - {self.tc_pos_mm[-1]:.1f} mm',
            '=' * 50
        ]
        return '\n'.join(lines)

# ====================================================================
# 1D Mesh Class
# ====================================================================

class Mesh1D:
    '''
    1D mesh for axial discretization

    Handles both uniform mesh and mesh aligned with TC positions
    '''

    def __init__(self, geometry: ReactorGeometry, dz: float = 1e-3):
        self.geometry = geometry
        self.dz = dz

        self._generate_uniform_mesh()
        self._find_tc_indices()

    def _generate_uniform_mesh(self):
        """Generate uniform mesh with specified spacing"""
        self.n_nodes = int(np.ceil(self.geometry.L / self.dz)) + 1
        self.z = np.linspace(0, self.geometry.L, self.n_nodes)
        self.dz_actual = self.z[1] - self.z[0]

    def _find_tc_indices(self):
        """Find mesh indices closest to TC positions"""
        self.tc_indices = np.zeros(self.geometry.n_tc, dtype = int)
        for i, tc_pos in enumerate(self.geometry.tc_pos):
            self.tc_indices[i] = np.argmin(np.abs(self.z - tc_pos))
        
        # Also store the actual z positiions at TC indices
        self.tc_z = self.z[self.tc_indices]

    @property
    def z_mm(self) -> np.ndarray:
        '''Mesh positions in mm'''
        return self.z * 1000
    
    def summary(self) -> str:
        lines = [
            '=' * 50,
            '1D Mesh Summary',
            '=' * 50,
            f'Number of nodes:  {self.n_nodes}',
            f'Spacing (dz):     {self.dz_actual*1000:.3f} mm',
            f'z range:          {self.z_mm[0]:.1f} - {self.z_mm[-1]:.1f} mm',
            '',
            f'TC positions and mesh mapping:',
        ]

        for i, (tc_pos, idx) in enumerate(zip(self.geometry.tc_pos,self.tc_indices)):
            z_mesh = self.z[idx]
            error = (z_mesh - tc_pos) * 1000
            lines.append(f'     TC{i+1:2d}: {tc_pos*1000:6.1f} mm -> '
                            f'node {idx:4d} (z={z_mesh*1000:6.2f} mm, '
                            f'err={error:+.2f} mm')
        lines.append('=' * 50)
        return '\n'.join(lines)

    def plot_tc_positions(self, ax=None, show: bool = True):
        '''Visualize TC positions on the reactor'''
        if ax is None:
            fig, ax = plt.subplots(figsize=(12, 3))

        # Draw tube outline
        L_mm = self.geometry.L * 1000
        ri_mm = self.geometry.id * 1000
        ro_mm = self.geometry.od * 1000

        # Outer wall
        ax.plot([0, L_mm], [ro_mm, ro_mm], 'k-', linewidth=2)
        ax.plot([0, L_mm], [-ro_mm, -ro_mm], 'k-', linewidth=2)

        # Inner Wall (dashed)
        ax.plot([0, L_mm], [ri_mm, ri_mm], 'k--', linewidth=1, alpha=0.5)
        ax.plot([0, L_mm], [-ri_mm, -ri_mm], 'k--', linewidth=1, alpha=0.5)

        # End caps
        ax.plot([0, 0], [-ro_mm, ro_mm], 'k-', linewidth=2)
        ax.plot([L_mm, L_mm], [-ro_mm, ro_mm], 'k-', linewidth=2)

        # TC positions
        tc_mm = self.geometry.tc_pos * 1000
        for i, tc in enumerate(tc_mm):
            ax.axvline(tc, color='r', linestyle='-', alpha=0.7, linewidth=1)
            ax.plot(tc, ro_mm + 0.5, 'rv', markersize=8)
            ax.text(tc, ro_mm + 1, f'{i+1}', ha='center', va='bottom', fontsize=8)
        
        ax.set_xlabel('Axial position [mm]')
        ax.set_ylabel('Radial [mm]')
        ax.set_title('Reactor Geometry with TC Positions')
        ax.set_xlim(-10, L_mm + 10)
        ax.set_ylim(-ro_mm - 2, ro_mm + 3)
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        
        if show:
            plt.tight_layout()
            plt.show()
        
        return ax

# ====================================================================
# Operating Conditions Class
# ====================================================================

class OperatingConditions:
    '''
    Operating conditions for the reactor
    '''

    def __init__(self, Fv_sccm: float = 50.0,
                Tin: float = 25.0,
                P: float = 1.01325,
                Pw: float = 140.0,
                Tamb: float = 25.0,
                tend: float = 215 * 60,
                dt: float = 1.0):
        self.Fv_sccm = Fv_sccm
        self.Tin = Tin
        self.P = P
        self.Pw = Pw
        self.Tamb = Tamb
        self.tend = tend
        self.dt = dt

        # Calculate mass flow rate
        from Property import sccm_to_kg_per_s
        self.Fm = sccm_to_kg_per_s(self.Fv_sccm)

    def heat_generation_per_length(self, geometry: ReactorGeometry) -> float:
        '''
        Volumetric heat generation rate [W/m3]
        '''
        return self.Pw / geometry.Vw

    def summary(self) -> str:
        '''Return formatted summary'''
        lines= [
            '=' * 50,
            'Operating Conditions',
            '=' * 50,
            f'He flow rate      {self.Fv_sccm:.1f} sccm',
            f'                  {self.Fm*1e6:.4f} mg/s',
            f'Inlet temp        {self.Tin:.1f} C',
            f'Pressure          {self.P:.1f} bar',
            f'Ambient temp      {self.Tamb:.1f} C',
            f'Total time        {self.tend/60:.1f} min',
            f'Time step         {self.dt:.1f} s',
            '=' * 50
        ]
        return '\n'.join(lines)

# =============================================================================
# Test / Example Usage
# =============================================================================

if __name__ == "__main__":

    # Create geometry
    geom = ReactorGeometry()
    print(geom.summary())

    # Create mesh
    mesh = Mesh1D(geom, dz=1e-3)
    print(mesh.summary())

    # Create operating conditions
    ops = OperatingConditions()
    print(ops.summary())

    # Plot TC positions
    mesh.plot_tc_positions()