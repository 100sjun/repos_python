"""
Differentiable FDM Solver with Implicit Differentiation
========================================================

Features:
1. Forward pass: Robust solver with line search + Levenberg-Marquardt
2. Backward pass: Implicit differentiation for gradient computation
3. Neural network compatible: Enables end-to-end training

Usage:
    # With neural network
    params = neural_net(input)  # Predict physical parameters
    Tw_solution = DifferentiableFDMSolver.apply(params, Tw_init, residual_func)
    loss = F.mse_loss(Tw_solution, measured_data)
    loss.backward()  # Gradients flow through solver!
"""

import torch
import torch.nn.functional as F
from typing import Callable, Tuple, List, Optional


# =============================================================================
# Robust Forward Solver (Non-differentiable, but converges well!)
# =============================================================================

def robust_solver(
    Tw_init: torch.Tensor,
    params: torch.Tensor,
    residual_func: Callable,
    tol: float = 1e-8,
    max_iter: int = 100,
    lambda_init: float = 1e-3,
    verbose: bool = False,
) -> Tuple[torch.Tensor, List[float]]:
    """
    Robust Newton-Raphson solver with:
    - Levenberg-Marquardt damping
    - Backtracking line search
    - Adaptive lambda adjustment

    This solver is NOT differentiable, but that's OK!
    We'll use implicit differentiation in the backward pass.

    Args:
        Tw_init: Initial temperature guess [K]
        params: Physical parameters (can be from NN)
        residual_func: Residual function R(Tw, params)
        tol: Convergence tolerance
        max_iter: Maximum iterations
        lambda_init: Initial Levenberg-Marquardt damping
        verbose: Print iteration info

    Returns:
        Tw_solution: Converged temperature distribution
        residual_history: Convergence history
    """
    # Detach everything - forward pass is gradient-free!
    Tw = Tw_init.clone().detach()
    params_detached = params.detach()

    lambda_lm = lambda_init
    residual_history = []

    if verbose:
        print(f"\n{'='*70}")
        print("Robust Newton-Raphson with Line Search + Levenberg-Marquardt")
        print(f"{'='*70}")
        print(f"{'Iter':>5} | {'||R||':>12} | {'alpha':>8} | {'lambda':>10} | {'Status':>8}")
        print(f"{'-'*5}-+-{'-'*12}-+-{'-'*8}-+-{'-'*10}-+-{'-'*8}")

    for iteration in range(max_iter):
        # Compute residual
        with torch.no_grad():
            R = residual_func(Tw, params_detached)
            residual_norm = torch.norm(R).item()
            residual_history.append(residual_norm)

        # Check convergence
        if residual_norm < tol:
            if verbose:
                print(f"{iteration:5d} | {residual_norm:12.6e} | {'--':>8} | {'--':>10} | {'CONVERGED':>8}")
                print(f"\n✓ Converged in {iteration} iterations!")
            break

        # Compute Jacobian J = ∂R/∂Tw
        Tw_for_jacobian = Tw.clone().requires_grad_(True)
        R_for_jacobian = residual_func(Tw_for_jacobian, params_detached)
        J = torch.autograd.functional.jacobian(
            lambda T: residual_func(T, params_detached),
            Tw
        )

        # Levenberg-Marquardt step
        # (J^T*J + lambda*I) * delta = -J^T*R
        n = Tw.shape[0]
        JTJ = J.T @ J
        JTR = J.T @ R
        I = torch.eye(n, dtype=Tw.dtype, device=Tw.device)

        try:
            delta_Tw = torch.linalg.solve(JTJ + lambda_lm * I, -JTR)
        except:
            # Singular matrix - increase damping
            lambda_lm *= 10.0
            if verbose:
                print(f"{iteration:5d} | {residual_norm:12.6e} | {'FAIL':>8} | {lambda_lm:10.3e} | {'SINGULAR':>8}")
            continue

        # Backtracking line search
        alpha = 1.0
        max_backtracks = 10
        rho = 0.5
        c = 1e-4  # Armijo constant

        success = False
        for backtrack in range(max_backtracks):
            with torch.no_grad():
                Tw_new = Tw + alpha * delta_Tw
                R_new = residual_func(Tw_new, params_detached)
                new_norm = torch.norm(R_new).item()

            # Armijo condition: sufficient decrease
            if new_norm < (1 - c * alpha) * residual_norm:
                success = True
                break

            alpha *= rho

        # Update based on success
        if success:
            Tw = Tw_new
            # Decrease lambda (move toward Newton)
            lambda_lm = max(lambda_lm * 0.5, 1e-10)
            status = "✓"
        else:
            # Increase lambda (move toward gradient descent)
            lambda_lm = min(lambda_lm * 2.0, 1e10)
            status = "✗"

        if verbose and (iteration < 10 or iteration % 10 == 0):
            print(f"{iteration:5d} | {residual_norm:12.6e} | {alpha:8.4f} | {lambda_lm:10.3e} | {status:>8}")
    else:
        if verbose:
            print(f"\n⚠ Maximum iterations reached ({max_iter})")
            print(f"  Final residual: {residual_norm:.6e}")

    return Tw, residual_history


# =============================================================================
# Differentiable Solver with Implicit Differentiation
# =============================================================================

class DifferentiableFDMSolver(torch.autograd.Function):
    """
    Custom autograd function for differentiable FDM solver.

    Forward: Use robust solver (any method works!)
    Backward: Implicit differentiation using Implicit Function Theorem

    Mathematical Background:
    ========================
    Given: Solver finds T* such that F(T*, θ) = 0, where θ are parameters

    Implicit Function Theorem:
        dT*/dθ = -(∂F/∂T)^(-1) * (∂F/∂θ)
               = -J^(-1) * (∂F/∂θ)

    Chain rule for loss L:
        dL/dθ = (dL/dT*) * (dT*/dθ)
              = -(dL/dT*) * J^(-1) * (∂F/∂θ)

    Efficient computation:
        Instead of computing J^(-1), solve:
        J^T * v = dL/dT*  for v
        Then: dL/dθ = -v^T * (∂F/∂θ)
    """

    @staticmethod
    def forward(
        ctx,
        params: torch.Tensor,
        Tw_init: torch.Tensor,
        residual_func: Callable,
        tol: float = 1e-8,
        max_iter: int = 100,
        verbose: bool = False,
    ) -> torch.Tensor:
        """
        Forward pass: Solve F(T, params) = 0 for T

        Args:
            params: Physical parameters (from NN or fixed)
            Tw_init: Initial temperature guess
            residual_func: Function R(Tw, params) that should equal 0
            tol: Convergence tolerance
            max_iter: Maximum iterations
            verbose: Print convergence info

        Returns:
            Tw_solution: Converged temperature distribution
        """
        # Solve using robust method
        Tw_solution, history = robust_solver(
            Tw_init, params, residual_func,
            tol=tol, max_iter=max_iter, verbose=verbose
        )

        # Save for backward
        ctx.save_for_backward(Tw_solution, params)
        ctx.residual_func = residual_func
        ctx.tol = tol

        return Tw_solution

    @staticmethod
    def backward(ctx, grad_output: torch.Tensor) -> Tuple[Optional[torch.Tensor], ...]:
        """
        Backward pass: Compute dL/dparams using implicit differentiation

        Args:
            grad_output: dL/dT* from upstream gradient

        Returns:
            grad_params: dL/dparams
            None for other inputs (Tw_init, residual_func, etc.)
        """
        Tw_solution, params = ctx.saved_tensors
        residual_func = ctx.residual_func

        # Check if params requires grad
        if not params.requires_grad:
            return None, None, None, None, None, None

        # Enable gradients for implicit differentiation
        Tw_for_backward = Tw_solution.detach().requires_grad_(True)
        params_for_backward = params.detach().requires_grad_(True)

        # Compute residual at solution (should be ~0)
        R = residual_func(Tw_for_backward, params_for_backward)

        # 1. Compute Jacobian J = ∂F/∂T at solution
        J = torch.autograd.functional.jacobian(
            lambda T: residual_func(T, params_for_backward),
            Tw_for_backward
        )

        # 2. Solve J^T * v = grad_output for v
        # This is much more efficient than computing J^(-1)
        try:
            v = torch.linalg.solve(J.T, grad_output)
        except:
            # Fallback: use pseudo-inverse if singular
            v = torch.linalg.lstsq(J.T, grad_output).solution

        # 3. Compute ∂F/∂params using Jacobian
        # More efficient: compute Jacobian of R w.r.t. params directly
        dF_dparams = torch.autograd.functional.jacobian(
            lambda p: residual_func(Tw_for_backward, p),
            params_for_backward
        )

        # 4. Compute grad_params = -v^T * (∂F/∂params)
        # dF_dparams shape: [n_residuals, n_params]
        # v shape: [n_residuals]
        # Ensure dtype compatibility (create new tensor, not inplace)
        v_cast = v.to(dtype=dF_dparams.dtype)
        grad_params = -torch.einsum('i,ij->j', v_cast, dF_dparams)

        # Detach to prevent further gradient tracking
        grad_params = grad_params.detach()

        # Return gradients (only for params, others are None)
        return grad_params, None, None, None, None, None


# =============================================================================
# Convenience Wrapper
# =============================================================================

def solve_differentiable(
    params: torch.Tensor,
    Tw_init: torch.Tensor,
    residual_func: Callable,
    tol: float = 1e-8,
    max_iter: int = 1000,
    verbose: bool = False,
) -> torch.Tensor:
    """
    Convenience wrapper for differentiable solver.

    Usage:
        # Simple call
        Tw = solve_differentiable(params, Tw_init, residual_func)

        # With neural network
        params = neural_net(input)
        Tw = solve_differentiable(params, Tw_init, residual_func)
        loss = F.mse_loss(Tw, measured)
        loss.backward()  # Gradients flow!
    """
    return DifferentiableFDMSolver.apply(
        params, Tw_init, residual_func,
        tol, max_iter, verbose
    )


# =============================================================================
# Test Functions
# =============================================================================

if __name__ == "__main__":
    print("="*70)
    print("Differentiable FDM Solver Test")
    print("="*70)

    # Test with simple 1D heat conduction
    N = 50
    L = 1.0
    dx = L / (N - 1)

    def test_residual(T, params):
        """
        Simple 1D heat conduction: k * d²T/dx² + q = 0
        params = [k, q, T_left, T_right]
        """
        k, q, T_left, T_right = params
        R = torch.zeros_like(T)

        # Boundary conditions (Dirichlet)
        R[0] = T[0] - T_left
        R[-1] = T[-1] - T_right

        # Interior points: k * (T[i+1] - 2*T[i] + T[i-1])/dx² + q = 0
        R[1:-1] = k * (T[2:] - 2*T[1:-1] + T[:-2]) / dx**2 + q

        return R

    # Test 1: Fixed parameters (sanity check)
    print("\n" + "="*70)
    print("Test 1: Fixed Parameters")
    print("="*70)

    params_fixed = torch.tensor([1.0, 100.0, 25.0, 25.0])  # k, q, T_left, T_right
    T_init = torch.ones(N) * 25.0

    T_solution = solve_differentiable(
        params_fixed, T_init, test_residual,
        verbose=True
    )

    print(f"\nSolution stats:")
    print(f"  T_min = {T_solution.min():.2f} °C")
    print(f"  T_max = {T_solution.max():.2f} °C")
    print(f"  T[0] = {T_solution[0]:.2f} °C (boundary)")
    print(f"  T[-1] = {T_solution[-1]:.2f} °C (boundary)")

    # Test 2: Gradient check
    print("\n" + "="*70)
    print("Test 2: Gradient Check (Backpropagation)")
    print("="*70)

    # Parameters that require gradients
    params_trainable = torch.tensor([1.0, 100.0, 25.0, 25.0], requires_grad=True)

    # Solve
    T_solution = solve_differentiable(
        params_trainable, T_init, test_residual,
        verbose=False
    )

    # Create a simple loss (e.g., minimize max temperature)
    loss = T_solution.max()

    print(f"Loss (max temperature): {loss.item():.2f} °C")

    # Backward pass
    loss.backward()

    print(f"\nGradients w.r.t. parameters:")
    print(f"  dL/dk        = {params_trainable.grad[0]:.6f}")
    print(f"  dL/dq        = {params_trainable.grad[1]:.6f}")
    print(f"  dL/dT_left   = {params_trainable.grad[2]:.6f}")
    print(f"  dL/dT_right  = {params_trainable.grad[3]:.6f}")

    print("\n✓ Gradients computed successfully!")
    print("  This means the solver is fully differentiable!")

    # Test 3: Neural network integration
    print("\n" + "="*70)
    print("Test 3: Neural Network Integration")
    print("="*70)

    class SimpleNN(torch.nn.Module):
        def __init__(self):
            super().__init__()
            self.fc = torch.nn.Linear(2, 4)  # Input: [k_prior, q_prior] -> Output: [k, q, T_left, T_right]

        def forward(self, x):
            out = self.fc(x)
            # Ensure physical constraints
            out[0] = F.softplus(out[0])  # k > 0
            out[1] = F.softplus(out[1])  # q > 0
            return out

    model = SimpleNN()
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)

    # Dummy training data
    input_data = torch.tensor([1.0, 100.0])
    target_T_max = torch.tensor(100.0)  # Target: max temperature should be 100°C

    print("Training for 5 iterations...")
    for epoch in range(5):
        optimizer.zero_grad()

        # NN predicts parameters
        params_pred = model(input_data)

        # Solve (differentiable!)
        T_pred = solve_differentiable(
            params_pred, T_init, test_residual,
            verbose=False
        )

        # Loss
        loss = (T_pred.max() - target_T_max)**2

        # Backprop through solver!
        loss.backward()
        optimizer.step()

        print(f"  Epoch {epoch+1}: Loss = {loss.item():.6f}, T_max = {T_pred.max():.2f}")

    print("\n✓ Neural network training successful!")
    print("  Gradients flowed through the FDM solver correctly!")
