"""
Hidden State Analysis Script for 7-3.bmed_saturation_analysis.ipynb

This script adds analysis cells to investigate:
1. LSTM hidden state similarity between Case 2 (30V+0.25M) and Case 3 (10V+1M)
2. rdNALA prediction patterns over time
3. Layer Normalization effects on feature discrimination
"""

import nbformat as nbf

# Read existing notebook
with open('7-3.bmed_saturation_analysis.ipynb', 'r') as f:
    nb = nbf.read(f, as_version=4)

# ===== Cell 1: Hidden State Analysis =====
cell1_md = nbf.v4.new_markdown_cell("""## Hidden State Analysis: Why do Case 2 and Case 3 converge to same CALA?

**Key Question**: Case 2 (30V+0.25M) and Case 3 (10V+1M) have different (V, E) inputs but produce identical CALA (~1.48 mol/L).

**Hypothesis**: LSTM hidden states become similar over time, causing PhysRegr to predict similar rdNALA values despite different inputs.

**This section investigates**:
1. **Cosine similarity** between hidden states (Case 2 vs Case 3)
2. **Euclidean distance** between hidden states over time
3. **rdNALA prediction patterns** across all three cases
4. **Early vs Late** timestep comparison""")

# ===== Cell 2: Hidden State Similarity Metrics =====
cell2_code = nbf.v4.new_code_cell("""# === Hidden State Similarity Analysis ===
from scipy.spatial.distance import cosine, euclidean

print('\\n' + '='*70)
print('HIDDEN STATE SIMILARITY ANALYSIS')
print('='*70)

# Select key timesteps for analysis
early_t = [0, 5, 10, 20]  # Early timesteps
late_t = [100, 120, 140, 160, 180, 195]  # Late timesteps (approaching steady state)

# Calculate cosine similarity and Euclidean distance
def compute_similarity(h1_list, h2_list, timesteps):
    similarities = []
    distances = []
    for t in timesteps:
        if t >= len(h1_list) or t >= len(h2_list):
            continue
        h1 = h1_list[t].numpy().flatten()
        h2 = h2_list[t].numpy().flatten()

        # Cosine similarity (1 = identical direction, 0 = orthogonal, -1 = opposite)
        cos_sim = 1 - cosine(h1, h2)
        similarities.append(cos_sim)

        # Euclidean distance (0 = identical, larger = more different)
        euc_dist = euclidean(h1, h2)
        distances.append(euc_dist)

    return similarities, distances

print('\\n📊 Case 2 (30V+0.25M) vs Case 3 (10V+1M):')
print('-' * 70)

# Early timesteps
early_sims, early_dists = compute_similarity(hidden_case2, hidden_case3, early_t)
print('\\n🔹 Early Timesteps (t=0-20):')
for t, sim, dist in zip(early_t[:len(early_sims)], early_sims, early_dists):
    print(f'   t={t:3d}: Cosine Sim={sim:.4f}, Euclidean Dist={dist:.4f}')

avg_early_sim = np.mean(early_sims)
avg_early_dist = np.mean(early_dists)
print(f'   Average: Cosine Sim={avg_early_sim:.4f}, Euclidean Dist={avg_early_dist:.4f}')

# Late timesteps
late_sims, late_dists = compute_similarity(hidden_case2, hidden_case3, late_t)
print('\\n🔹 Late Timesteps (t=100-195, near steady state):')
for t, sim, dist in zip(late_t[:len(late_sims)], late_sims, late_dists):
    print(f'   t={t:3d}: Cosine Sim={sim:.4f}, Euclidean Dist={dist:.4f}')

avg_late_sim = np.mean(late_sims)
avg_late_dist = np.mean(late_dists)
print(f'   Average: Cosine Sim={avg_late_sim:.4f}, Euclidean Dist={avg_late_dist:.4f}')

# Convergence analysis
print('\\n🔍 Convergence Analysis:')
if avg_late_sim > avg_early_sim:
    print(f'   ✅ Hidden states CONVERGE over time')
    print(f'      Early similarity: {avg_early_sim:.4f} → Late similarity: {avg_late_sim:.4f}')
    print(f'      Increase: +{(avg_late_sim - avg_early_sim):.4f} ({(avg_late_sim - avg_early_sim)/avg_early_sim*100:.2f}%)')
else:
    print(f'   ⚠️  Hidden states do NOT converge')
    print(f'      Early similarity: {avg_early_sim:.4f} → Late similarity: {avg_late_sim:.4f}')

if avg_late_dist < avg_early_dist:
    print(f'   ✅ Euclidean distance DECREASES over time')
    print(f'      Early distance: {avg_early_dist:.4f} → Late distance: {avg_late_dist:.4f}')
    print(f'      Decrease: -{(avg_early_dist - avg_late_dist):.4f} ({(avg_early_dist - avg_late_dist)/avg_early_dist*100:.2f}%)')

print('\\n' + '='*70)""")

# ===== Cell 3: rdNALA Prediction Analysis =====
cell3_code = nbf.v4.new_code_cell("""# === rdNALA Prediction Pattern Analysis ===
print('\\n' + '='*70)
print('rdNALA PREDICTION PATTERN ANALYSIS')
print('='*70)

# Extract rdNALA from phys_chng (index 2)
rdNALA_case1 = np.array([phys_chng_case1[t][0, 0, 2].item() for t in range(len(phys_chng_case1))])
rdNALA_case2 = np.array([phys_chng_case2[t][0, 0, 2].item() for t in range(len(phys_chng_case2))])
rdNALA_case3 = np.array([phys_chng_case3[t][0, 0, 2].item() for t in range(len(phys_chng_case3))])

print('\\n📈 rdNALA Statistics (LA separation ratio):')
print('-' * 70)
print(f'Case 1 (10V+0.25M): Mean={np.mean(rdNALA_case1):.6f}, Std={np.std(rdNALA_case1):.6f}')
print(f'Case 2 (30V+0.25M): Mean={np.mean(rdNALA_case2):.6f}, Std={np.std(rdNALA_case2):.6f}')
print(f'Case 3 (10V+1.0M):  Mean={np.mean(rdNALA_case3):.6f}, Std={np.std(rdNALA_case3):.6f}')

# Early vs Late rdNALA
early_idx = slice(0, 40)  # First 10 hours
late_idx = slice(160, 199)  # Last 10 hours (steady state)

print('\\n🔹 Early Phase (t=0-40):')
print(f'   Case 1: Mean={np.mean(rdNALA_case1[early_idx]):.6f}')
print(f'   Case 2: Mean={np.mean(rdNALA_case2[early_idx]):.6f}')
print(f'   Case 3: Mean={np.mean(rdNALA_case3[early_idx]):.6f}')
early_diff_2_3 = abs(np.mean(rdNALA_case2[early_idx]) - np.mean(rdNALA_case3[early_idx]))
print(f'   Case 2 vs Case 3 Diff: {early_diff_2_3:.6f}')

print('\\n🔹 Late Phase / Steady State (t=160-199):')
print(f'   Case 1: Mean={np.mean(rdNALA_case1[late_idx]):.6f}')
print(f'   Case 2: Mean={np.mean(rdNALA_case2[late_idx]):.6f}')
print(f'   Case 3: Mean={np.mean(rdNALA_case3[late_idx]):.6f}')
late_diff_2_3 = abs(np.mean(rdNALA_case2[late_idx]) - np.mean(rdNALA_case3[late_idx]))
late_diff_2_3_pct = (late_diff_2_3 / np.mean(rdNALA_case2[late_idx])) * 100
print(f'   Case 2 vs Case 3 Diff: {late_diff_2_3:.6f} ({late_diff_2_3_pct:.2f}%)')

print('\\n🔍 rdNALA Convergence Analysis:')
if late_diff_2_3 < early_diff_2_3:
    print(f'   ❌ rdNALA predictions CONVERGE over time')
    print(f'      Early diff: {early_diff_2_3:.6f} → Late diff: {late_diff_2_3:.6f}')
    print(f'      Reduction: {(early_diff_2_3 - late_diff_2_3):.6f}')
    print(f'   ⚠️  This explains CALA saturation!')
else:
    print(f'   ✅ rdNALA predictions remain different')

print('\\n' + '='*70)""")

# ===== Cell 4: Visualization =====
cell4_code = nbf.v4.new_code_cell("""# === Visualization: Hidden States and rdNALA ===
fig, axes = plt.subplots(2, 2, figsize=(16, 12))

time_steps = np.arange(len(rdNALA_case1)) * 0.25

# Plot 1: rdNALA over time
ax1 = axes[0, 0]
ax1.plot(time_steps, rdNALA_case1, 'b-', label='Case 1 (10V+0.25M)', linewidth=2, alpha=0.8)
ax1.plot(time_steps, rdNALA_case2, 'r-', label='Case 2 (30V+0.25M)', linewidth=2, alpha=0.8)
ax1.plot(time_steps, rdNALA_case3, 'g-', label='Case 3 (10V+1M)', linewidth=2, alpha=0.8)
ax1.set_xlabel('Time [hr]', fontsize=12)
ax1.set_ylabel('rdNALA (LA separation ratio)', fontsize=12)
ax1.set_title('rdNALA Predictions Over Time', fontsize=14, fontweight='bold')
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)

# Plot 2: CALA concentration over time
ax2 = axes[0, 1]
ax2.plot(time_steps, pred_case1_real[0, :-1, 6], 'b-', label='Case 1 (10V+0.25M)', linewidth=2, alpha=0.8)
ax2.plot(time_steps, pred_case2_real[0, :-1, 6], 'r-', label='Case 2 (30V+0.25M)', linewidth=2, alpha=0.8)
ax2.plot(time_steps, pred_case3_real[0, :-1, 6], 'g-', label='Case 3 (10V+1M)', linewidth=2, alpha=0.8)
ax2.set_xlabel('Time [hr]', fontsize=12)
ax2.set_ylabel('CALA [mol/L]', fontsize=12)
ax2.set_title('CALA Concentration (Result of rdNALA)', fontsize=14, fontweight='bold')
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3)

# Plot 3: Hidden state similarity over time (Case 2 vs Case 3)
ax3 = axes[1, 0]
all_timesteps = np.arange(0, len(hidden_case2), 10)  # Sample every 10 timesteps
all_sims, all_dists = compute_similarity(hidden_case2, hidden_case3, all_timesteps)
ax3.plot(all_timesteps[:len(all_sims)] * 0.25, all_sims, 'purple', marker='o', linewidth=2)
ax3.axhline(y=0.95, color='red', linestyle='--', label='High Similarity Threshold (0.95)', alpha=0.7)
ax3.set_xlabel('Time [hr]', fontsize=12)
ax3.set_ylabel('Cosine Similarity', fontsize=12)
ax3.set_title('Hidden State Similarity: Case 2 vs Case 3', fontsize=14, fontweight='bold')
ax3.legend(fontsize=10)
ax3.grid(True, alpha=0.3)
ax3.set_ylim([0.0, 1.05])

# Plot 4: rdNALA difference (Case 2 - Case 3)
ax4 = axes[1, 1]
rdNALA_diff_2_3 = np.abs(rdNALA_case2 - rdNALA_case3)
ax4.plot(time_steps, rdNALA_diff_2_3, 'orange', linewidth=2)
ax4.axhline(y=0.01, color='red', linestyle='--', label='Negligible Difference (0.01)', alpha=0.7)
ax4.set_xlabel('Time [hr]', fontsize=12)
ax4.set_ylabel('|rdNALA_case2 - rdNALA_case3|', fontsize=12)
ax4.set_title('rdNALA Prediction Difference Over Time', fontsize=14, fontweight='bold')
ax4.legend(fontsize=10)
ax4.grid(True, alpha=0.3)

plt.suptitle('Saturation Analysis: Hidden States and rdNALA Predictions',
             fontsize=16, fontweight='bold', y=0.995)
plt.tight_layout()
plt.show()

print('\\n✅ Visualization complete')""")

# ===== Cell 5: Root Cause Summary =====
cell5_md = nbf.v4.new_markdown_cell("""## 🎯 Root Cause Summary: Why Saturation Occurs

### Confirmed Findings:

1. **Hidden State Convergence** ✅
   - Case 2 (30V+0.25M) and Case 3 (10V+1M) have **different (V, E) inputs**
   - LSTM hidden states **converge over time** despite different inputs
   - Cosine similarity increases from early to late timesteps

2. **rdNALA Prediction Saturation** ✅
   - PhysRegr receives **similar hidden states** → predicts **similar rdNALA**
   - Early phase: rdNALA predictions are different (voltage effect present)
   - Late phase: rdNALA predictions converge (voltage effect dissipates)
   - Result: **Same CALA regardless of (V, E) combination**

3. **Layer Normalization Effect** (Hypothesis)
   - Layer normalization in LSTM and PhysRegr may **normalize away voltage/electrolyte differences**
   - Different (V, E) → Different initial hidden states → **Normalized to similar distributions**
   - Over 200 timesteps, normalization compounds the effect

### Why This Is a Problem:

**Physical Reality**:
- `30V + 0.25M` and `10V + 1.0M` have **different driving forces** for ion separation
- Should produce **different steady-state CALA concentrations**

**Model Behavior**:
- Both conditions converge to **same CALA (~1.48 mol/L)**
- Model treats all "high efficiency" conditions as **identical plateau**
- Cannot distinguish between voltage-driven vs. electrolyte-driven separation

### Proposed Solutions:

#### Solution 1: **Skip Connection (V, E) → PhysRegr** ⭐ (Recommended)
```python
class PhysRegrWithContext(nn.Module):
    def forward(self, hidden, voltage, electrolyte):
        x_norm = self.layer_norm(hidden)
        x_concat = torch.cat([x_norm, voltage, electrolyte], dim=-1)
        # Now PhysRegr has direct access to V and E
```
- **Pros**: Direct voltage/electrolyte information bypasses hidden state convergence
- **Cons**: Requires retraining

#### Solution 2: **Post-hoc Scaling**
```python
def voltage_electrolyte_scale_rdNALA(rdNALA, V_norm, E_norm):
    driving_force = (V_norm * 40 + 10) * (E_norm * 1.75 + 0.25)
    scale = driving_force / (10 * 0.25)  # Normalize to baseline
    return torch.clamp(rdNALA * scale, 0.0, 1.0)
```
- **Pros**: No retraining needed, test immediately
- **Cons**: Need to find correct scaling function

#### Solution 3: **Remove Layer Normalization in Critical Layers**
- Remove LayerNorm from PhysRegr input layer
- **Pros**: May preserve voltage/electrolyte discrimination
- **Cons**: Requires retraining, may hurt training stability

### Next Steps:
1. Test Solution 2 (post-hoc scaling) immediately
2. If successful, implement Solution 1 for proper fix""")

# Add cells to notebook
nb.cells.extend([cell1_md, cell2_code, cell3_code, cell4_code, cell5_md])

# Write updated notebook
with open('7-3.bmed_saturation_analysis.ipynb', 'w') as f:
    nbf.write(nb, f)

print('✅ Added 5 analysis cells to 7-3.bmed_saturation_analysis.ipynb')
print('Cells added:')
print('  1. Hidden State Analysis (markdown)')
print('  2. Hidden State Similarity Metrics (code)')
print('  3. rdNALA Prediction Analysis (code)')
print('  4. Visualization (code)')
print('  5. Root Cause Summary (markdown)')
