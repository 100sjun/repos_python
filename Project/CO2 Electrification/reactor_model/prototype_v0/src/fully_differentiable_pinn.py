"""
Fully Differentiable Physics-Informed Neural Network with CasADi
================================================================
Complete gradient flow: Temperature data → NN → CasADi (with sensitivity) → Loss
No need for k_true in training!
"""

import numpy as np
import casadi as ca
import matplotlib.pyplot as plt
from typing import Tuple, Dict, List
import torch
import torch.nn as nn
import torch.optim as optim


class ThermalConductivityNN(nn.Module):
    """Neural Network to predict thermal conductivity"""

    def __init__(self, input_dim: int = 5, hidden_sizes: list = [32, 16, 8]):
        super(ThermalConductivityNN, self).__init__()

        layers = []
        prev_size = input_dim

        for hidden_size in hidden_sizes:
            layers.append(nn.Linear(prev_size, hidden_size))
            layers.append(nn.Tanh())
            prev_size = hidden_size

        layers.append(nn.Linear(prev_size, 1))
        layers.append(nn.Softplus())  # k > 0

        self.network = nn.Sequential(*layers)

    def forward(self, x):
        return self.network(x)


class CasADiPhysicsSolver:
    """CasADi-based physics solver with sensitivity computation"""

    def __init__(self, length: float = 1.0, n_nodes: int = 50):
        self.L = length
        self.n = n_nodes
        self.dx = length / (n_nodes - 1)
        self.x = np.linspace(0, length, n_nodes)

        # Build symbolic problem once (for efficiency)
        self._build_symbolic_problem()

    def _build_symbolic_problem(self):
        """Build symbolic CasADi problem with k as parameter"""

        # Symbolic variables
        T_sym = ca.MX.sym('T', self.n)
        k_sym = ca.MX.sym('k')
        T_left_sym = ca.MX.sym('T_left')
        T_right_sym = ca.MX.sym('T_right')
        q_gen_sym = ca.MX.sym('q_gen')

        # Objective: minimize heat equation residual
        obj = 0
        for i in range(1, self.n - 1):
            d2T_dx2 = (T_sym[i+1] - 2*T_sym[i] + T_sym[i-1]) / (self.dx**2)
            residual = k_sym * d2T_dx2 + q_gen_sym
            obj += residual**2

        # Constraints: boundary conditions
        g = ca.vertcat(T_sym[0] - T_left_sym, T_sym[-1] - T_right_sym)

        # Parameters
        p = ca.vertcat(k_sym, T_left_sym, T_right_sym, q_gen_sym)

        # NLP problem
        nlp = {
            'x': T_sym,
            'p': p,
            'f': obj,
            'g': g
        }

        # Solver options
        opts = {
            'ipopt.print_level': 0,
            'print_time': 0,
            'ipopt.max_iter': 1000
        }

        # Create solver
        self.solver = ca.nlpsol('solver', 'ipopt', nlp, opts)

        print("✓ Symbolic CasADi problem built successfully")

    def solve(self, k: float, T_left: float, T_right: float,
             q_gen: float) -> np.ndarray:
        """
        Solve heat equation with given parameters

        Parameters:
        -----------
        k : float
            Thermal conductivity [W/m·K]
        T_left : float
            Left boundary temperature [K]
        T_right : float
            Right boundary temperature [K]
        q_gen : float
            Volumetric heat generation [W/m³]

        Returns:
        --------
        T : np.ndarray
            Temperature distribution [K]
        """
        # Initial guess
        T0 = np.linspace(T_left, T_right, self.n)

        # Parameters
        p = [k, T_left, T_right, q_gen]

        # Solve
        sol = self.solver(x0=T0, p=p, lbg=0, ubg=0)

        return np.array(sol['x']).flatten()

    def compute_sensitivity(self, k: float, T_left: float, T_right: float,
                           q_gen: float, T_optimal: np.ndarray) -> np.ndarray:
        """
        Compute sensitivity ∂T/∂k using CasADi

        Uses implicit function theorem:
        If F(T*, k) = 0 (optimality condition), then:
        ∂T*/∂k = -(∂F/∂T)^(-1) × (∂F/∂k)

        Parameters:
        -----------
        k : float
            Thermal conductivity [W/m·K]
        T_left, T_right : float
            Boundary temperatures [K]
        q_gen : float
            Heat generation [W/m³]
        T_optimal : np.ndarray
            Optimal temperature from forward solve

        Returns:
        --------
        dT_dk : np.ndarray
            Sensitivity ∂T[i]/∂k for each node i
        """
        # Build sensitivity problem
        T_sym = ca.MX.sym('T', self.n)
        k_sym = ca.MX.sym('k')
        T_left_sym = ca.MX.sym('T_left')
        T_right_sym = ca.MX.sym('T_right')
        q_gen_sym = ca.MX.sym('q_gen')

        # Heat equation residuals (KKT conditions)
        residuals = []
        for i in range(1, self.n - 1):
            d2T_dx2 = (T_sym[i+1] - 2*T_sym[i] + T_sym[i-1]) / (self.dx**2)
            residual = k_sym * d2T_dx2 + q_gen_sym
            residuals.append(residual)

        # Boundary conditions
        residuals.insert(0, T_sym[0] - T_left_sym)
        residuals.append(T_sym[-1] - T_right_sym)

        F = ca.vertcat(*residuals)

        # Compute Jacobians
        dF_dT = ca.jacobian(F, T_sym)  # ∂F/∂T
        dF_dk = ca.jacobian(F, k_sym)  # ∂F/∂k

        # Create functions
        jac_T_func = ca.Function('jac_T', [T_sym, k_sym, T_left_sym, T_right_sym, q_gen_sym], [dF_dT])
        jac_k_func = ca.Function('jac_k', [T_sym, k_sym, T_left_sym, T_right_sym, q_gen_sym], [dF_dk])

        # Evaluate at optimal point
        dF_dT_val = jac_T_func(T_optimal, k, T_left, T_right, q_gen)
        dF_dk_val = jac_k_func(T_optimal, k, T_left, T_right, q_gen)

        # Solve linear system: dF_dT × dT_dk = -dF_dk
        dF_dT_dense = np.array(dF_dT_val.full())
        dF_dk_dense = np.array(dF_dk_val.full()).flatten()

        # ∂T/∂k = -(∂F/∂T)^(-1) × (∂F/∂k)
        dT_dk = -np.linalg.solve(dF_dT_dense, dF_dk_dense)

        return dT_dk


class CasADiPhysicsLayer(torch.autograd.Function):
    """
    Custom PyTorch layer that wraps CasADi solver
    Enables gradient flow through physics solver
    """

    @staticmethod
    def forward(ctx, k_tensor, physics_solver, T_left, T_right, q_gen):
        """
        Forward pass: solve physics with IPOPT

        Parameters:
        -----------
        ctx : context
            PyTorch context for saving info
        k_tensor : torch.Tensor
            Thermal conductivity (requires_grad=True)
        physics_solver : CasADiPhysicsSolver
            Physics solver instance
        T_left, T_right, q_gen : float
            Boundary conditions and heat generation

        Returns:
        --------
        T_pred : torch.Tensor
            Predicted temperature distribution
        """
        k = k_tensor.item()

        # Solve with IPOPT
        T_optimal = physics_solver.solve(k, T_left, T_right, q_gen)

        # Save for backward
        ctx.save_for_backward(k_tensor)
        ctx.physics_solver = physics_solver
        ctx.T_optimal = T_optimal
        ctx.params = (T_left, T_right, q_gen)

        return torch.tensor(T_optimal, dtype=torch.float32)

    @staticmethod
    def backward(ctx, grad_output):
        """
        Backward pass: compute ∂Loss/∂k using chain rule

        Chain rule:
        ∂Loss/∂k = ∂Loss/∂T × ∂T/∂k
                   ↑         ↑
              grad_output  CasADi sensitivity

        Parameters:
        -----------
        ctx : context
            Saved context from forward
        grad_output : torch.Tensor
            Gradient from loss ∂Loss/∂T

        Returns:
        --------
        grad_k : torch.Tensor
            Gradient ∂Loss/∂k
        """
        k_tensor = ctx.saved_tensors[0]
        physics_solver = ctx.physics_solver
        T_optimal = ctx.T_optimal
        T_left, T_right, q_gen = ctx.params

        # Compute sensitivity ∂T/∂k using CasADi
        dT_dk = physics_solver.compute_sensitivity(
            k_tensor.item(), T_left, T_right, q_gen, T_optimal
        )

        # Chain rule: ∂Loss/∂k = (∂Loss/∂T) · (∂T/∂k)
        grad_output_np = grad_output.detach().numpy()
        grad_k = np.dot(grad_output_np, dT_dk)

        # Convert to tensor
        grad_k_tensor = torch.tensor(grad_k, dtype=torch.float32)

        return grad_k_tensor, None, None, None, None


def extract_features(T_data: np.ndarray, T_left: float,
                    T_right: float, q_gen: float, L: float = 1.0) -> np.ndarray:
    """Extract features from temperature data"""
    T_avg = np.mean(T_data)
    dT_dx_avg = (T_right - T_left) / L
    q_normalized = q_gen / 1000.0
    T_var = np.var(T_data)
    T_max = np.max(T_data)

    return np.array([T_avg, dT_dx_avg, q_normalized, T_var, T_max])


def generate_experimental_data(n_samples: int = 100, n_nodes: int = 50,
                               noise_level: float = 0.5) -> List[Dict]:
    """Generate synthetic experimental data"""
    physics_solver = CasADiPhysicsSolver(length=1.0, n_nodes=n_nodes)
    exp_data = []

    print(f"\nGenerating {n_samples} experimental datasets...")

    for i in range(n_samples):
        k_true = np.random.uniform(10, 80)
        T_left = np.random.uniform(350, 450)
        T_right = np.random.uniform(300, 400)
        q_gen = np.random.uniform(0, 800)

        T_true = physics_solver.solve(k_true, T_left, T_right, q_gen)
        T_measured = T_true + np.random.normal(0, noise_level, n_nodes)

        exp_data.append({
            'T_measured': T_measured,
            'T_left': T_left,
            'T_right': T_right,
            'q_gen': q_gen,
            'k_true': k_true
        })

        if (i + 1) % 20 == 0:
            print(f"  Generated {i + 1}/{n_samples} samples")

    return exp_data


def train_fully_differentiable(nn_model: ThermalConductivityNN,
                               physics_solver: CasADiPhysicsSolver,
                               train_data: List[Dict],
                               val_data: List[Dict],
                               n_epochs: int = 50,
                               lr: float = 0.01) -> Dict:
    """
    Train with fully differentiable pipeline

    Key difference: NO k_true needed in loss!
    Only temperature data is used.
    """
    optimizer = optim.Adam(nn_model.parameters(), lr=lr)
    physics_layer = CasADiPhysicsLayer.apply

    history = {
        'train_loss': [],
        'val_loss': [],
        'train_k_error': [],
        'val_k_error': []
    }

    print(f"\n{'='*70}")
    print("Training Fully Differentiable PINN (NO k_true in loss!)")
    print(f"{'='*70}\n")

    for epoch in range(n_epochs):
        # ========== Training ==========
        nn_model.train()
        train_loss = 0
        train_k_error = 0

        for data in train_data:
            # Extract features
            features = extract_features(
                data['T_measured'],
                data['T_left'],
                data['T_right'],
                data['q_gen']
            )
            features_tensor = torch.FloatTensor(features).unsqueeze(0)

            # NN predicts k
            k_pred_tensor = nn_model(features_tensor).squeeze()

            # Physics layer (fully differentiable!)
            T_pred_tensor = physics_layer(
                k_pred_tensor,
                physics_solver,
                data['T_left'],
                data['T_right'],
                data['q_gen']
            )

            # Loss: ONLY temperature fitting (no k_true!)
            T_measured_tensor = torch.FloatTensor(data['T_measured'])
            loss = torch.mean((T_pred_tensor - T_measured_tensor)**2)

            # Backward (gradient flows through CasADi!)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            train_loss += loss.item()
            train_k_error += abs(k_pred_tensor.item() - data['k_true'])

        train_loss /= len(train_data)
        train_k_error /= len(train_data)

        # ========== Validation ==========
        nn_model.eval()
        val_loss = 0
        val_k_error = 0

        with torch.no_grad():
            for data in val_data:
                features = extract_features(
                    data['T_measured'],
                    data['T_left'],
                    data['T_right'],
                    data['q_gen']
                )
                features_tensor = torch.FloatTensor(features).unsqueeze(0)
                k_pred_tensor = nn_model(features_tensor).squeeze()

                # Forward pass
                T_pred = physics_solver.solve(
                    k_pred_tensor.item(),
                    data['T_left'],
                    data['T_right'],
                    data['q_gen']
                )

                loss = np.mean((T_pred - data['T_measured'])**2)
                val_loss += loss
                val_k_error += abs(k_pred_tensor.item() - data['k_true'])

        val_loss /= len(val_data)
        val_k_error /= len(val_data)

        history['train_loss'].append(train_loss)
        history['val_loss'].append(val_loss)
        history['train_k_error'].append(train_k_error)
        history['val_k_error'].append(val_k_error)

        if (epoch + 1) % 10 == 0:
            print(f"Epoch [{epoch+1}/{n_epochs}]")
            print(f"  Train - Loss: {train_loss:.4f}, k Error: {train_k_error:.2f} W/m·K")
            print(f"  Val   - Loss: {val_loss:.4f}, k Error: {val_k_error:.2f} W/m·K")

    return history


def plot_results(history: Dict, test_data: List[Dict],
                nn_model: ThermalConductivityNN,
                physics_solver: CasADiPhysicsSolver):
    """Plot training results"""

    # Plot 1: Training history
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    axes[0].plot(history['train_loss'], label='Train Loss', linewidth=2)
    axes[0].plot(history['val_loss'], label='Val Loss', linewidth=2)
    axes[0].set_xlabel('Epoch', fontsize=12)
    axes[0].set_ylabel('MSE Loss (Temperature)', fontsize=12)
    axes[0].set_title('Training History: Temperature Fitting\n(NO k_true in loss!)', fontsize=13)
    axes[0].legend(fontsize=11)
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(history['train_k_error'], label='Train k Error', linewidth=2)
    axes[1].plot(history['val_k_error'], label='Val k Error', linewidth=2)
    axes[1].set_xlabel('Epoch', fontsize=12)
    axes[1].set_ylabel('k Error [W/m·K]', fontsize=12)
    axes[1].set_title('k Prediction Error\n(learned from temperature data only)', fontsize=13)
    axes[1].legend(fontsize=11)
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('fully_diff_training_history.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("\n✓ Saved: fully_diff_training_history.png")

    # Plot 2: Test predictions
    nn_model.eval()
    n_test_plots = min(4, len(test_data))
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    k_true_list = []
    k_pred_list = []

    for i in range(n_test_plots):
        data = test_data[i]

        features = extract_features(
            data['T_measured'],
            data['T_left'],
            data['T_right'],
            data['q_gen']
        )
        features_tensor = torch.FloatTensor(features).unsqueeze(0)

        with torch.no_grad():
            k_pred = nn_model(features_tensor).item()

        T_pred = physics_solver.solve(
            k_pred, data['T_left'], data['T_right'], data['q_gen']
        )

        axes[i].plot(physics_solver.x, data['T_measured'], 'o',
                    label='Measured', alpha=0.6, markersize=5)
        axes[i].plot(physics_solver.x, T_pred, '-',
                    label=f'Predicted', linewidth=2)
        axes[i].set_xlabel('Position [m]', fontsize=11)
        axes[i].set_ylabel('Temperature [K]', fontsize=11)
        axes[i].set_title(f'Test {i+1}: k_true={data["k_true"]:.1f}, k_pred={k_pred:.1f} W/m·K',
                         fontsize=12)
        axes[i].legend(fontsize=10)
        axes[i].grid(True, alpha=0.3)

        k_true_list.append(data['k_true'])
        k_pred_list.append(k_pred)

    plt.tight_layout()
    plt.savefig('fully_diff_test_predictions.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("✓ Saved: fully_diff_test_predictions.png")

    # Plot 3: k prediction accuracy (all test data)
    for data in test_data[n_test_plots:]:
        features = extract_features(
            data['T_measured'],
            data['T_left'],
            data['T_right'],
            data['q_gen']
        )
        features_tensor = torch.FloatTensor(features).unsqueeze(0)

        with torch.no_grad():
            k_pred = nn_model(features_tensor).item()

        k_true_list.append(data['k_true'])
        k_pred_list.append(k_pred)

    k_true_arr = np.array(k_true_list)
    k_pred_arr = np.array(k_pred_list)

    plt.figure(figsize=(8, 8))
    plt.scatter(k_true_arr, k_pred_arr, alpha=0.6, s=50)
    plt.plot([k_true_arr.min(), k_true_arr.max()],
            [k_true_arr.min(), k_true_arr.max()], 'r--', lw=2, label='Perfect prediction')
    plt.xlabel('True Thermal Conductivity [W/m·K]', fontsize=12)
    plt.ylabel('Predicted Thermal Conductivity [W/m·K]', fontsize=12)
    plt.title('k Prediction Accuracy\n(Fully Differentiable PINN)', fontsize=13)
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=11)

    r2 = 1 - np.sum((k_true_arr - k_pred_arr)**2) / np.sum((k_true_arr - np.mean(k_true_arr))**2)
    mae = np.mean(np.abs(k_true_arr - k_pred_arr))
    mape = np.mean(np.abs((k_true_arr - k_pred_arr) / k_true_arr)) * 100

    textstr = f'R² = {r2:.4f}\nMAE = {mae:.2f} W/m·K\nMAPE = {mape:.2f}%'
    plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes,
            verticalalignment='top', fontsize=11,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))

    plt.tight_layout()
    plt.savefig('fully_diff_k_prediction.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("✓ Saved: fully_diff_k_prediction.png")


def main():
    """Main demonstration"""

    print("\n" + "="*70)
    print("Fully Differentiable Physics-Informed Neural Network")
    print("="*70)
    print("\nKey Feature: Gradient flows through CasADi solver!")
    print("Training uses ONLY temperature data (no k_true in loss)")

    # Parameters
    n_nodes = 50
    n_samples = 100
    noise_level = 1.0

    # Step 1: Generate data
    print("\n[Step 1] Generating experimental data...")
    exp_data = generate_experimental_data(n_samples, n_nodes, noise_level)

    # Split
    n_train = int(0.7 * n_samples)
    n_val = int(0.15 * n_samples)

    train_data = exp_data[:n_train]
    val_data = exp_data[n_train:n_train+n_val]
    test_data = exp_data[n_train+n_val:]

    print(f"\nDataset split:")
    print(f"  Training:   {len(train_data)} samples")
    print(f"  Validation: {len(val_data)} samples")
    print(f"  Test:       {len(test_data)} samples")

    # Step 2: Initialize models
    print("\n[Step 2] Initializing models...")
    physics_solver = CasADiPhysicsSolver(length=1.0, n_nodes=n_nodes)
    nn_model = ThermalConductivityNN(input_dim=5, hidden_sizes=[32, 16, 8])

    print("✓ Neural network initialized")

    # Step 3: Train
    print("\n[Step 3] Training fully differentiable model...")
    history = train_fully_differentiable(
        nn_model,
        physics_solver,
        train_data,
        val_data,
        n_epochs=50,
        lr=0.01
    )

    # Step 4: Evaluate
    print("\n[Step 4] Evaluating and visualizing results...")
    plot_results(history, test_data, nn_model, physics_solver)

    # Final statistics
    print("\n" + "="*70)
    print("Final Results:")
    print("="*70)
    print(f"Final Train Loss: {history['train_loss'][-1]:.4f}")
    print(f"Final Val Loss:   {history['val_loss'][-1]:.4f}")
    print(f"Final Train k Error: {history['train_k_error'][-1]:.2f} W/m·K")
    print(f"Final Val k Error:   {history['val_k_error'][-1]:.2f} W/m·K")

    print("\n✓ All figures saved:")
    print("  - fully_diff_training_history.png")
    print("  - fully_diff_test_predictions.png")
    print("  - fully_diff_k_prediction.png")
    print("\n" + "="*70)


if __name__ == "__main__":
    try:
        import casadi
        import torch
        main()
    except ImportError as e:
        print(f"Error: Required package not found - {e}")
        print("\nPlease install: pip install casadi torch matplotlib")
