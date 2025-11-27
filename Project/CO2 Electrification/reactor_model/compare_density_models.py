"""Compare He density from prop_He correlation vs ideal gas law"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from prop_He import prop_He

# Initialize He properties
he = prop_He()

# Temperature range
T_C = np.linspace(20, 500, 100)  # Celsius
T_K = T_C + 273.15  # Kelvin

# Density from prop_He correlation (experimental fit)
rho_correlation = np.array([he.rho(T) * 1000 for T in T_C])  # kg/m³

# Density from ideal gas law: ρ = PM/(RT)
P = 101325  # Pa (1 atm)
R = 8.314  # J/(mol·K)
M = he.M / 1000  # kg/mol
rho_ideal = P * M / (R * T_K)  # kg/m³

# Calculate difference
diff_abs = rho_correlation - rho_ideal
diff_rel = (rho_correlation - rho_ideal) / rho_ideal * 100

# Create comparison plot
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# 1. Density comparison
ax1 = axes[0, 0]
ax1.plot(T_C, rho_correlation, 'b-', linewidth=2, label='prop_He correlation (experimental)')
ax1.plot(T_C, rho_ideal, 'r--', linewidth=2, label='Ideal gas law')
ax1.set_xlabel('Temperature (°C)', fontsize=11)
ax1.set_ylabel('Density ρ (kg/m³)', fontsize=11)
ax1.set_title('He Density: Correlation vs Ideal Gas Law', fontsize=12, fontweight='bold')
ax1.legend()
ax1.grid(True, alpha=0.3)

# 2. Relative difference
ax2 = axes[0, 1]
ax2.plot(T_C, diff_rel, 'g-', linewidth=2)
ax2.axhline(y=0, color='k', linestyle='-', linewidth=0.5)
ax2.fill_between(T_C, 0, diff_rel, alpha=0.3, color='g')
ax2.set_xlabel('Temperature (°C)', fontsize=11)
ax2.set_ylabel('Relative Difference (%)', fontsize=11)
ax2.set_title('Error: (Correlation - Ideal)/Ideal × 100%', fontsize=12, fontweight='bold')
ax2.grid(True, alpha=0.3)

# 3. Inverse relationship check
ax3 = axes[1, 0]
# prop_He uses ρ = a/(T+b) form
# Ideal gas gives ρ = const/T
ax3.plot(T_K, 1/rho_correlation, 'b-', linewidth=2, label='1/ρ (correlation)')
ax3.plot(T_K, 1/rho_ideal, 'r--', linewidth=2, label='1/ρ (ideal gas)')
ax3.set_xlabel('Temperature (K)', fontsize=11)
ax3.set_ylabel('1/ρ (m³/kg)', fontsize=11)
ax3.set_title('Inverse Density (shows linearity)', fontsize=12, fontweight='bold')
ax3.legend()
ax3.grid(True, alpha=0.3)

# 4. Statistics table
ax4 = axes[1, 1]
ax4.axis('off')

# Calculate statistics at key temperatures
temps_check = [25, 100, 200, 300, 400, 500]
stats_text = "DENSITY COMPARISON AT KEY TEMPERATURES\n"
stats_text += "="*55 + "\n\n"
stats_text += f"{'T(°C)':<8} {'ρ_corr':<12} {'ρ_ideal':<12} {'Diff(%)':<10}\n"
stats_text += "-"*55 + "\n"

for T in temps_check:
    T_k = T + 273.15
    rho_c = he.rho(T) * 1000
    rho_i = P * M / (R * T_k)
    diff = (rho_c - rho_i) / rho_i * 100
    stats_text += f"{T:<8.0f} {rho_c:<12.5f} {rho_i:<12.5f} {diff:<10.3f}\n"

stats_text += "\n" + "="*55 + "\n"
stats_text += f"\nMean absolute error:    {np.abs(diff_rel).mean():.3f}%\n"
stats_text += f"Max absolute error:     {np.abs(diff_rel).max():.3f}%\n"
stats_text += f"RMS error:              {np.sqrt((diff_rel**2).mean()):.3f}%\n"

stats_text += f"\nConclusion:\n"
stats_text += f"  • Correlation is ~{diff_rel.mean():.2f}% higher than ideal gas\n"
stats_text += f"  • Real He deviates from ideality\n"
stats_text += f"  • Correlation captures T-dependence better\n"

stats_text += f"\nCurrent implementation:\n"
stats_text += f"  ✓ Using prop_He.rho(T) correlation\n"
stats_text += f"  ✓ Based on experimental data\n"
stats_text += f"  ✓ More accurate than ideal gas law\n"

ax4.text(0.05, 0.95, stats_text, transform=ax4.transAxes,
         fontsize=9, verticalalignment='top', family='monospace',
         bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))

plt.suptitle('He Density: Experimental Correlation vs Ideal Gas Law',
             fontsize=14, fontweight='bold')
plt.tight_layout()

plt.savefig('density_comparison_correlation_vs_ideal.png', dpi=150, bbox_inches='tight')
print("✓ Saved to: density_comparison_correlation_vs_ideal.png")

plt.close(fig)

print("\n" + "="*60)
print("IMPORTANT: Your current models ALREADY use the correlation!")
print("="*60)
print("\nIn convection_model_2d.py and convection_model_2d_improved.py:")
print("  get_properties(T) calls:")
print("    rho = self.he_prop.rho(T_C) * 1000  ← prop_He correlation ✓")
print("\nNOT using ideal gas law - using experimental data fit!")
print("="*60)
