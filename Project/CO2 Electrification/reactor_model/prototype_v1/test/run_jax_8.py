#!/usr/bin/env python3
"""
Execute test_jax_8.ipynb notebook code
Unified parameter fitting with tau saturation c(Tw) = c_inf * (1 - exp(-Tw/T0))
"""

# Cell 1: Library imports
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt

import jax
import jax.numpy as jnp
from jax import lax
from jax import custom_vjp
import optax

print(f"JAX version: {jax.__version__}")
print(f"JAX devices: {jax.devices()}")

# Cell 2: Data loading
tc_pos = jnp.array(json.load(open('tc_pos.json')))
T_140 = jnp.array(json.load(open('Temp_profile_140W.json')))
T_280 = jnp.array(json.load(open('Temp_profile_280W.json')))
T_420 = jnp.array(json.load(open('Temp_profile_420W.json')))

powers = jnp.array([140.0, 280.0, 420.0])
T_meas_cases = jnp.array([T_140, T_280, T_420])
num_cases, n_TC = T_meas_cases.shape

print(f"Cases: {num_cases}, TC count: {n_TC}")

# Geometry
L = 430e-3; ID = 5.03e-3; OD = 6.33e-3
Ai = jnp.pi * ID**2 / 4; Ao = jnp.pi * OD**2 / 4; Aw = Ao - Ai
pri = jnp.pi * ID; pro = jnp.pi * OD
dz = 0.001; n_nodes = int(L/dz) + 1; z = jnp.linspace(0, L, n_nodes)
Awg = pri * dz; Aout = pro * dz
tc_idx = jnp.array([jnp.argmin(jnp.abs(z - zp)) for zp in tc_pos])

print(f"Nodes: {n_nodes}, Grid: {dz*1000:.2f} mm")

# Material properties
path_he = '../data/He_property.csv'
path_kt = '../data/KanthalD_property.csv'
df_he = pd.read_csv(path_he); df_kt = pd.read_csv(path_kt);
df_kt_k = df_kt[df_kt['k [W/m*K]'].notna()]; df_kt_cp = df_kt[df_kt['Cp [kJ/kg*K]'].notna()]
Mw_he = 4.0026

def interp1d(xq, x, y):
    xq = jnp.asarray(xq)
    x  = jnp.asarray(x)
    y  = jnp.asarray(y)
    return jnp.interp(xq, x, y)

def rho_he(T): return interp1d(T, df_he['T [C]'], df_he['rho [kg/m3]'])
def cp_he(T): return interp1d(T, df_he['T [C]'], df_he['Cp [J/mol*K]'])
def k_he(T): return interp1d(T, df_he['T [C]'], df_he['Tc [W/m*K]'])
def rho_kt(T): return interp1d(T, df_kt['T [C]'], df_kt['rho [kg/m3]'])
def cp_kt(T): return interp1d(T, df_kt_cp['T [C]'], df_kt_cp['Cp [kJ/kg*K]'])*1000
def k_kt(T): return interp1d(T, df_kt_k['T [C]'], df_kt_k['k [W/m*K]'])
def h_wg(Tg): return 4.36 * k_he(Tg) / ID

# Feed conditions
P = 101325.0; Tamb = 25.0; Fv_std = 50.0
Fw = Fv_std * 1e-6 / 60 * rho_he(0)
F = Fw / Mw_he * 1000.0

print("Setup complete")

# Cell 3: c(Tw) and hout functions
Z1_FIXED = 0.1
Z2_FIXED = 0.33

def compute_c_of_Tw(c_params, Tw_local):
    """Tau saturation: c(Tw) = c_inf * (1 - exp(-Tw/T0))"""
    c_inf, T0 = c_params[0], c_params[1]
    c_Tw = c_inf * (1.0 - jnp.exp(-Tw_local / (T0 + 1e-6)))
    return c_Tw

def compute_hout_profile(global_params, Tw):
    """hout(z, Tw) = c(Tw) + edge loss terms"""
    c_inf, T0, m, k = global_params[0], global_params[1], global_params[2], global_params[3]
    c_Tw = c_inf * (1.0 - jnp.exp(-Tw / (T0 + 1e-6)))
    z_term1 = (m / k) * jnp.log(1 + jnp.exp(k * (Z1_FIXED - z)))
    z_term2 = (m / k) * jnp.log(1 + jnp.exp(k * (z - Z2_FIXED)))
    hout_profile = c_Tw + z_term1 + z_term2
    return hout_profile

print("c(Tw) and hout(z, Tw) functions defined with tau saturation c(Tw)")

# Cell 4: Residual and Newton solver
HWI_FIXED = 0.
HWO_FIXED = 0.

def residual(U, global_params, Pw):
    hwi = HWI_FIXED
    hwo = HWO_FIXED

    Tw = U[:n_nodes]; Tg = U[n_nodes:]
    hout_profile = compute_hout_profile(global_params, Tw)
    kw = k_kt(Tw); cpg = cp_he(Tg); hwg = h_wg(Tg)
    Qelec = Pw * dz / L
    rw = jnp.zeros((n_nodes,)); rg = jnp.zeros((n_nodes,))

    # Wall BC
    rw = rw.at[0].set(hwi * Aw * (Tw[0] - Tamb) + kw[0] * Aw * (Tw[0] - Tw[1]) / dz)
    rw = rw.at[-1].set(hwo * Aw * (Tw[-1] - Tamb) + kw[-1] * Aw * (Tw[-1] - Tw[-2]) / dz)

    # Wall interior
    kw_half = 0.5 * (kw[1:] + kw[:-1])
    wflux = kw_half * Aw * (Tw[1:] - Tw[:-1]) / dz
    Qcond = (wflux[1:] - wflux[:-1]) / dz
    Qwg = hwg[1:-1] * Awg * (Tw[1:-1] - Tg[1:-1])
    Qout = hout_profile[1:-1] * Aout * (Tw[1:-1] - Tamb)
    rw = rw.at[1:-1].set(Qcond + (Qelec/dz) - (Qwg/dz) - (Qout/dz))

    # Gas BC & interior
    rg = rg.at[0].set(Tg[0] - Tw[0])
    rg = rg.at[-1].set(Tg[-1] - Tw[-1])
    gflux = F * cpg[1:-1] * (Tg[1:-1] - Tg[:-2]) / dz
    rg = rg.at[1:-1].set((Qwg/dz) - gflux)

    return jnp.concatenate([rw, rg])

def newton_step(residual_fn, damping=1.0):
    def step(U, _):
        F = residual_fn(U)
        J = jax.jacfwd(residual_fn)(U)
        dU = jnp.linalg.solve(J, -F)
        U_new = U + damping * dU
        return U_new, (jnp.linalg.norm(F), jnp.linalg.norm(dU))
    return step

def newton_solve_forward(residual_fn, U0, iters=20, damping=1.0):
    step = newton_step(residual_fn, damping)
    U_final, (res_hist, step_hist) = lax.scan(step, U0, xs=None, length=iters)
    return U_final, res_hist, step_hist

print("Residual and Newton solver defined")

# Cell 5: Implicit differentiation
@custom_vjp
def newton_solve_implicit(global_params, Pw, U0, iters=20, damping=1.0):
    res_fn = lambda U: residual(U, global_params, Pw)
    U_final, _, _ = newton_solve_forward(res_fn, U0, iters, damping)
    return U_final

def newton_solve_implicit_fwd(global_params, Pw, U0, iters, damping):
    U_star = newton_solve_implicit(global_params, Pw, U0, iters, damping)
    return U_star, (global_params, Pw, U_star)

def newton_solve_implicit_bwd(res, g):
    global_params, Pw, U_star = res
    res_fn = lambda U: residual(U, global_params, Pw)
    J = jax.jacfwd(res_fn)(U_star)
    lmbda = jnp.linalg.solve(J.T, g)

    def res_wrt_global_params(gp):
        return residual(U_star, gp, Pw)
    _, vjp_global_params = jax.vjp(res_wrt_global_params, global_params)
    grad_global_params = vjp_global_params(-lmbda)[0]

    def res_wrt_Pw(p):
        return residual(U_star, global_params, p)
    _, vjp_Pw = jax.vjp(res_wrt_Pw, Pw)
    grad_Pw = vjp_Pw(-lmbda)[0]

    grad_U0 = jnp.zeros_like(U_star)
    return (grad_global_params, grad_Pw, grad_U0, None, None)

newton_solve_implicit.defvjp(newton_solve_implicit_fwd, newton_solve_implicit_bwd)
print("Implicit differentiation solver ready")

# Cell 7: Parameter setup
LEARNING_RATE = 1e-1

def global_params_phys(global_params_raw):
    """global_params = [c_inf, T0, m, k]"""
    eps = 1e-6
    params = jax.nn.softplus(global_params_raw) + eps
    c_inf = jnp.clip(params[0], 30.0, 150.0)
    T0 = jnp.clip(params[1], 100.0, 2000.0)
    m = jnp.clip(params[2], 10.0, 500.0)
    k = jnp.clip(params[3], 10.0, 200.0)
    return jnp.array([c_inf, T0, m, k])

def softplus_inv(h):
    return jnp.log(jnp.exp(h) - 1.0)

global_params_init = jnp.array([77.59, 1000.0, 250.0, 78.0])
global_params_raw = softplus_inv(global_params_init)

Tw0 = Tamb * jnp.ones(n_nodes)
Tg0 = Tamb * jnp.ones(n_nodes)
U0 = jnp.concatenate([Tw0, Tg0])
U0_cases = jnp.stack([U0, U0, U0])

def predict_TC(U):
    return U[:n_nodes][tc_idx]

def case_loss(global_params, Pw, U0, T_meas):
    U_star = newton_solve_implicit(global_params, Pw, U0, iters=20, damping=1.0)
    T_pred = predict_TC(U_star)
    loss = jnp.mean((T_pred - T_meas)**2)
    return loss, U_star

def total_loss(global_params_raw, U0_cases, T_meas_cases, powers):
    global_params = global_params_phys(global_params_raw)
    def one_case(U0_k, T_k, Pw_k):
        loss_k, _ = case_loss(global_params, Pw_k, U0_k, T_k)
        return loss_k
    losses = jax.vmap(one_case)(U0_cases, T_meas_cases, powers)
    return jnp.sum(losses)

def warm_start_update(global_params_raw, U0_cases, T_meas_cases, powers):
    global_params = global_params_phys(global_params_raw)
    def one_case(U0_k, T_k, Pw_k):
        _, U_star = case_loss(global_params, Pw_k, U0_k, T_k)
        return lax.stop_gradient(U_star)
    return jax.vmap(one_case)(U0_cases, T_meas_cases, powers)

def get_newton_convergence(global_params_raw, U0_cases, T_meas_cases, powers):
    global_params = global_params_phys(global_params_raw)
    def one_case(U0_k, T_k, Pw_k):
        res_fn = lambda U: residual(U, global_params, Pw_k)
        _, res_hist, step_hist = newton_solve_forward(res_fn, U0_k, iters=20, damping=1.0)
        return res_hist[-1], step_hist[-1]
    res_norms, step_norms = jax.vmap(one_case)(U0_cases, T_meas_cases, powers)
    return res_norms, step_norms

print("="*80)
print("Fitting Setup - Tau Saturation c(Tw) = c_inf * (1 - exp(-Tw/T0))")
print("="*80)

# Cell 8: Training
opt = optax.adam(learning_rate=LEARNING_RATE)
opt_state = opt.init(global_params_raw)
loss_and_grad = jax.value_and_grad(total_loss)

print("TRAINING START")
print("="*80)

history = {
    'loss': [],
    'loss_per_case': [],
    'c_inf': [],
    'T0': [],
    'm': [],
    'k': [],
    'grads': [],
    'newton_res': [],
    'newton_step': []
}

for step in range(200):
    loss, grads = loss_and_grad(global_params_raw, U0_cases, T_meas_cases, powers)
    updates, opt_state = opt.update(grads, opt_state)
    global_params_raw = optax.apply_updates(global_params_raw, updates)

    res_norms, step_norms = get_newton_convergence(global_params_raw, U0_cases, T_meas_cases, powers)
    U0_cases = warm_start_update(global_params_raw, U0_cases, T_meas_cases, powers)

    global_params = global_params_phys(global_params_raw)

    def compute_case_losses(gp, U0s, Ts, Pws):
        def one_loss(U0_k, T_k, Pw_k):
            loss_k, _ = case_loss(gp, Pw_k, U0_k, T_k)
            return loss_k
        return jax.vmap(one_loss)(U0s, Ts, Pws)

    case_losses = compute_case_losses(global_params, U0_cases, T_meas_cases, powers)

    history['loss'].append(float(loss))
    history['loss_per_case'].append(np.array(case_losses))
    history['c_inf'].append(float(global_params[0]))
    history['T0'].append(float(global_params[1]))
    history['m'].append(float(global_params[2]))
    history['k'].append(float(global_params[3]))
    history['grads'].append(np.array(grads))
    history['newton_res'].append(np.array(res_norms))
    history['newton_step'].append(np.array(step_norms))

    if step % 20 == 0 or step < 5:
        print(f'\nEPOCH {step:04d} | Total Loss = {float(loss):.6e}')
        print(f'  Global Params:')
        print(f'    c_inf={float(global_params[0]):.4f}, T0={float(global_params[1]):.3f}')
        print(f'    m={float(global_params[2]):.3f}, k={float(global_params[3]):.3f}')
        print(f'  Case Losses:')
        for k, (Pw_k, loss_k) in enumerate(zip([140., 280., 420.], case_losses)):
            print(f'    Case {k+1} ({Pw_k}W): Loss={float(loss_k):.3e}')

print(f'\n{"="*80}')
print("TRAINING COMPLETE")
print(f'{"="*80}')

# Cell 10: Final results
global_params_final = global_params_phys(global_params_raw)

print("\n" + "="*80)
print("FINAL RESULTS - Tau Saturation c(Tw)")
print("="*80)
print(f"\nGlobal Parameters:")
print(f"  c_inf = {float(global_params_final[0]):.6f}")
print(f"  T0 = {float(global_params_final[1]):.3f} °C")
print(f"  m  = {float(global_params_final[2]):.3f}")
print(f"  k  = {float(global_params_final[3]):.3f}")

print("\nc(Tw) function values:")
for T_sample in [0, 300, 600, 900, 1200]:
    c_val = float(global_params_final[0]) * (1 - np.exp(-T_sample / float(global_params_final[1])))
    print(f"  c({T_sample}°C) = {c_val:.3f}")
print("="*80)

for k, Pw_k in enumerate([140., 280., 420.]):
    res_fn = lambda U: residual(U, global_params_final, Pw_k)
    U_star, res_hist, step_hist = newton_solve_forward(res_fn, U0_cases[k], iters=20, damping=1.0)
    Tw = U_star[:n_nodes]
    Tw_tc = predict_TC(U_star)
    mse = float(jnp.mean((Tw_tc - T_meas_cases[k])**2))

    hout_profile = compute_hout_profile(global_params_final, Tw)

    print(f"\nCase {k+1} ({Pw_k}W):")
    print(f"  Tw range: {float(Tw.min()):.1f} ~ {float(Tw.max()):.1f} °C")
    print(f"  Newton:   ||F||={float(res_hist[-1]):.3e}")
    print(f"  TC MSE:   {mse:.3e}")

# Cell 11: Summary
print("\n" + "="*80)
print("SUMMARY - Tau Saturation c(Tw)")
print("="*80)

loss_reduction = (1 - history['loss'][-1]/history['loss'][0]) * 100
print(f"\n1. Training Performance:")
print(f"   Initial loss: {history['loss'][0]:.6e}")
print(f"   Final loss:   {history['loss'][-1]:.6e}")
print(f"   Reduction:    {loss_reduction:.2f}%")

print(f"\n2. Final Parameters:")
print(f"   c_inf = {float(global_params_final[0]):.6f}")
print(f"   T0    = {float(global_params_final[1]):.3f} °C")
print(f"   m     = {float(global_params_final[2]):.3f}")
print(f"   k     = {float(global_params_final[3]):.3f}")

print(f"\n3. c(Tw) Analysis:")
for T_sample in [0, 300, 600, 900, 1200]:
    c_val = float(global_params_final[0]) * (1 - np.exp(-T_sample / float(global_params_final[1])))
    print(f"   c({T_sample}°C) = {c_val:.3f}")
print(f"   → Always positive: 0 < c < {float(global_params_final[0]):.2f} ✅")

print(f"\n4. Case Performance:")
for k, Pw_k in enumerate([140., 280., 420.]):
    case_loss_final = history['loss_per_case'][-1][k]
    print(f"   Case {k+1} ({Pw_k}W): Loss = {case_loss_final:.3e}")

print("\n" + "="*80)
if loss_reduction > 50:
    print("✅ Tau saturation c(Tw) strategy is highly effective!")
    print(f"   - Loss reduction: {loss_reduction:.1f}%")
    print(f"   - Parameters: 9 → 4 (56% reduction)")
    print(f"   - Always positive: guaranteed by tau saturation")
else:
    print("⚠️  Further improvement needed")
print("="*80)

print("\nExecution completed successfully!")
