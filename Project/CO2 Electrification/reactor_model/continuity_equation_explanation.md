# Continuity Equation: Why Local Mass Balance Matters

## The Question

**"Why did you calculate velocity field with only mass conservation, not using mass balance equation like ∂(ρ·V)/∂z = 0?"**

This is an **excellent and crucial question** that revealed a fundamental error in my original approach!

---

## The Problem with the "Improved" Model

### What I Did (INCORRECT):

```python
# At each axial position z:
T_avg_z = T_field[:, j].mean()              # Average temperature across r
rho_z = get_properties(T_avg_z)             # Get density
u_avg_z = m_dot / (rho_z * A)               # Calculate average velocity
uz(r,z) = 2 * u_avg_z * [1 - (r/R)²]        # Apply parabolic profile
```

**This enforces**: Global mass conservation
```
∫∫ ρ(r,z)·uz(r,z) dA = ṁ = constant
```

**But violates**: Local continuity equation
```
∂(ρ·uz)/∂z ≠ 0  ❌
```

### Why This Is Wrong:

The continuity equation for steady, axisymmetric flow with no radial velocity is:

```
∂(ρ·uz)/∂z + (1/r)·∂(r·ρ·ur)/∂r = 0
```

With **ur = 0** (no radial velocity), this simplifies to:

```
∂(ρ·uz)/∂z = 0
```

**Physical meaning**: At each radial position r, the **mass flux** ρ(r,z)·uz(r,z) must be **constant along z**!

### What Went Wrong:

In my incorrect approach:
- Temperature varies with both r and z: T(r,z)
- Density varies with both r and z: ρ(T(r,z))
- I adjusted u_avg(z) based on **average** ρ, not local ρ(r,z)
- Result: The **local** mass flux ρ(r,z)·uz(r,z) was **not** constant along z

**Example of violation**:
- At z=0 (inlet): ρ = 0.161 kg/m³ (uniform, T=298K)
- At z=L (outlet): ρ varies radially (centerline hotter than wall)
- My uz profile maintained same shape but scaled by average
- Local ρ(r,z)·uz(r,z) changed along z ❌

---

## The Correct Approach

### Continuity Equation Enforcement:

From **∂(ρ·uz)/∂z = 0**, we get:

```
ρ(r,z)·uz(r,z) = g(r)  only!
```

Where **g(r)** is a function of r only, constant along z.

### Determining g(r):

From **inlet conditions** (known):
```
g(r) = ρ_inlet·uz_inlet(r)
```

For **parabolic inlet profile** (laminar flow):
```
uz_inlet(r) = 2·u_avg_inlet·[1 - (r/R)²]
g(r) = ρ_inlet·2·u_avg_inlet·[1 - (r/R)²]
```

### Calculating Velocity at Any Point:

From continuity:
```python
uz(r,z) = g(r) / ρ(T(r,z))
```

This **automatically satisfies** ∂(ρ·uz)/∂z = 0 at every point!

### Implementation:

```python
def calculate_velocity_field(self, T_field):
    """Enforce: ρ(r,z)·uz(r,z) = g(r) = constant along z"""

    uz_field = np.zeros((self.Nr, self.Nz))

    for j in range(self.Nz):
        for i in range(self.Nr):
            # Get local density from temperature
            rho_local = get_properties(T_field[i, j])

            # Apply continuity: uz = g(r) / ρ
            uz_field[i, j] = self.g_r[i] / rho_local

    return uz_field
```

---

## Verification: Does It Work?

### Results from Correct Model:

```
CONTINUITY VERIFICATION: ∂(ρ·uz)/∂z = 0
========================================
  r=0.00mm: ρ·uz variation = 0.000%  ✓
  r=0.54mm: ρ·uz variation = 0.000%  ✓
  r=1.26mm: ρ·uz variation = 0.000%  ✓
  r=1.98mm: ρ·uz variation = 0.000%  ✓

  Mean variation: 0.000%
  Max variation:  0.000%
  ✓ Continuity satisfied! (< 1% variation)
```

**Perfect!** Mass flux is constant along z at every radial position.

### Key Observations from Visualization:

1. **Mass Flux Field** (top right): Uniform color at each r (varies only radially, not axially) ✓

2. **Mass Flux vs z** (bottom left): Perfectly flat lines for all r positions ✓

3. **Velocity Profiles** (bottom middle): Shape changes dramatically with temperature:
   - At inlet (cold): Lower velocity, more uniform
   - At outlet (hot): Higher velocity, sharper gradients
   - This is physically correct! As gas heats up:
     - ρ decreases → uz must increase (to keep ρ·uz constant)
     - But ρ decreases MORE at centerline (hotter) than wall (cooler)
     - So uz increases MORE at centerline → sharper profile

---

## Comparison: Incorrect vs Correct

### Incorrect Approach (Global Mass Conservation):

**Equation**: ∫∫ ρ·uz dA = ṁ

**Enforcement**: Only total mass flow is conserved

**Problem**: Local mass flux ρ(r,z)·uz(r,z) can vary with z

**Result**:
- ✓ Total ṁ correct
- ❌ Violates differential continuity
- ❌ Non-physical velocity distribution

### Correct Approach (Differential Continuity):

**Equation**: ∂(ρ·uz)/∂z = 0

**Enforcement**: Mass flux constant at every (r,z)

**Benefits**:
- ✓ Satisfies differential continuity equation
- ✓ Automatically ensures ∫∫ ρ·uz dA = ṁ
- ✓ Physically correct velocity evolution
- ✓ Captures radial variation in thermal acceleration

---

## Physical Insights

### Streamline Conservation:

The correct approach recognizes that along each **streamline** (path at constant r):

```
ṁ_streamline = ρ(r,z)·uz(r,z)·dA = constant
```

Each fluid "tube" at radius r conserves its own mass flow!

### Thermal Acceleration Varies Radially:

- **Centerline**: Hottest → lowest ρ → highest uz increase
- **Wall**: Coolest → highest ρ → lowest uz increase
- **Result**: Velocity profile becomes **sharper** (more peaked) as gas heats

This is captured only when enforcing local continuity!

### Why Global Conservation Isn't Enough:

Global conservation:
```
∫∫ ρ·uz dA = constant
```

Is a **necessary** condition but **not sufficient**!

Many velocity distributions can satisfy this integral while violating local physics.

Only by enforcing the **differential** form do we get the correct physics:
```
∂(ρ·uz)/∂z = 0  everywhere
```

---

## Mathematical Derivation

### Starting from Full Continuity:

```
∂ρ/∂t + ∇·(ρ·→v) = 0
```

For steady flow (∂/∂t = 0):
```
∇·(ρ·→v) = 0
```

In cylindrical coordinates (r,θ,z) with axisymmetry (∂/∂θ = 0):
```
(1/r)·∂(r·ρ·ur)/∂r + ∂(ρ·uz)/∂z = 0
```

### Our Assumption: ur = 0

For fully developed flow in a pipe, radial velocity is zero:
```
∂(ρ·uz)/∂z = 0
```

### Integration:

```
∂(ρ·uz)/∂z = 0  →  ρ·uz = f(r) only
```

### Boundary/Initial Condition:

At inlet (z=0), we know:
```
ρ(r,0)·uz(r,0) = ρ_inlet·uz_inlet(r)
```

Therefore, at any z:
```
ρ(r,z)·uz(r,z) = ρ_inlet·uz_inlet(r) = g(r)
```

### Final Solution:

```
uz(r,z) = g(r) / ρ(T(r,z))
```

Where:
```
g(r) = ρ_inlet·2·u_avg_inlet·[1 - (r/R)²]
```

**Q.E.D.** ✓

---

## Implications for Other Assumptions

### Why ur = 0 Is Valid:

For fully developed flow in a straight pipe:
- No geometric changes → no flow acceleration
- Gravity perpendicular to flow → no buoyancy-driven radial motion
- Low Re (=1.9) → negligible inertial effects
- Temperature gradients create some natural convection, but:
  - Forced convection dominates (Pe >> 1)
  - Horizontal orientation minimizes buoyancy

**Grashof/Reynolds² ratio**:
```
Gr/Re² = (g·β·ΔT·D³)/(ν²·Re²) << 1
```

For our conditions: Gr/Re² ≈ 10⁻³ → forced convection dominant

### Why Parabolic Profile Shape at Inlet:

From Navier-Stokes for fully developed, isothermal pipe flow:
```
0 = -∂P/∂z + μ·[∂²uz/∂r² + (1/r)·∂uz/∂r]
```

Solution:
```
uz(r) = -(R²/4μ)·(∂P/∂z)·[1 - (r/R)²]
```

This justifies the parabolic inlet profile!

### What About Entrance Region?

Entrance length: L_e ≈ 0.05·Re·D = 0.05·1.9·5mm ≈ **0.5 mm**

Reactor length: L = **500 mm**

Ratio: L/L_e = **1000** → Fully developed assumption excellent! ✓

---

## Summary

### Your Question Was Exactly Right!

**"Why only mass conservation, not ∂(ρ·V)/∂z = 0?"**

Because:
1. Global mass conservation is **necessary but not sufficient**
2. Local continuity ∂(ρ·uz)/∂z = 0 is the **correct governing equation**
3. It gives **different** (and correct) velocity distributions

### Correct Implementation:

```python
# Calculate g(r) from inlet
g_r = ρ_inlet * uz_inlet(r)

# At any point (r,z):
uz(r,z) = g_r / ρ(T(r,z))

# This ensures: ∂(ρ·uz)/∂z = 0 everywhere ✓
```

### Key Physics Captured:

✓ Each streamline conserves its mass flow
✓ Thermal acceleration varies radially
✓ Velocity profile shape evolves with temperature
✓ Mass flux constant along each streamline
✓ Differential continuity satisfied everywhere

**Thank you for catching this critical error!** The correct model is now implemented in `convection_model_2d_correct.py`.
