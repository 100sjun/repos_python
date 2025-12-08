"""
Fully Differentiable PINN with Full Temperature Profile Input
==============================================================
Uses entire temperature profile (50 nodes) as NN input instead of summary statistics
This should significantly improve k prediction accuracy!
"""

import numpy as np
import casadi as ca
import matplotlib.pyplot as plt
from typing import Tuple, Dict, List
import torch
import torch.nn as nn
import torch.optim as optim


class ThermalConductivityNN(nn.Module):
    """Neural Network that takes full temperature profile as input"""

    def __init__(self, n_nodes: int = 50, hidden_sizes: list = [64, 32, 16]):
        """
        Parameters:
        -----------
        n_nodes : int
            Number of spatial nodes (input dimension)
        hidden_sizes : list
            Hidden layer sizes
        """
        super(ThermalConductivityNN, self).__init__()

        layers = []
        prev_size = n_nodes

        for hidden_size in hidden_sizes:
            layers.append(nn.Linear(prev_size, hidden_size))
            layers.append(nn.Tanh())
            layers.append(nn.Dropout(0.1))
            prev_size = hidden_size

        layers.append(nn.Linear(prev_size, 1))
        layers.append(nn.Softplus())  # k > 0

        self.network = nn.Sequential(*layers)

    def forward(self, x):
        """
        x: temperature profile [batch_size, n_nodes]
        output: thermal conductivity [batch_size, 1]
        """
        return self.network(x)


class CasADiPhysicsSolver:
    """CasADi-based physics solver with sensitivity computation"""

    def __init__(self, length: float = 1.0, n_nodes: int = 50):
        self.L = length
        self.n = n_nodes
        self.dx = length / (n_nodes - 1)
        self.x = np.linspace(0, length, n_nodes)
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

        self.solver = ca.nlpsol('solver', 'ipopt', nlp, opts)
        print("✓ Symbolic CasADi problem built")

    def solve(self, k: float, T_left: float, T_right: float,
             q_gen: float) -> np.ndarray:
        """Solve heat equation with IPOPT"""
        T0 = np.linspace(T_left, T_right, self.n)
        p = [k, T_left, T_right, q_gen]
        sol = self.solver(x0=T0, p=p, lbg=0, ubg=0)
        return np.array(sol['x']).flatten()

    def compute_sensitivity(self, k: float, T_left: float, T_right: float,
                           q_gen: float, T_optimal: np.ndarray) -> np.ndarray:
        """
        Compute sensitivity ∂T/∂k using implicit function theorem

        Returns:
        --------
        dT_dk : np.ndarray, shape (n_nodes,)
            Sensitivity at each node
        """
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
        dF_dT = ca.jacobian(F, T_sym)
        dF_dk = ca.jacobian(F, k_sym)

        # Create functions
        jac_T_func = ca.Function('jac_T', [T_sym, k_sym, T_left_sym, T_right_sym, q_gen_sym], [dF_dT])
        jac_k_func = ca.Function('jac_k', [T_sym, k_sym, T_left_sym, T_right_sym, q_gen_sym], [dF_dk])

        # Evaluate at optimal point
        dF_dT_val = jac_T_func(T_optimal, k, T_left, T_right, q_gen)
        dF_dk_val = jac_k_func(T_optimal, k, T_left, T_right, q_gen)

        # Solve: dF_dT × dT_dk = -dF_dk
        dF_dT_dense = np.array(dF_dT_val.full())
        dF_dk_dense = np.array(dF_dk_val.full()).flatten()

        # ∂T/∂k
        dT_dk = -np.linalg.solve(dF_dT_dense, dF_dk_dense)

        return dT_dk


class CasADiPhysicsLayer(torch.autograd.Function):
    """Custom PyTorch layer with gradient through CasADi"""

    @staticmethod
    def forward(ctx, k_tensor, physics_solver, T_left, T_right, q_gen):
        """Forward: solve with IPOPT"""
        k = k_tensor.item()
        T_optimal = physics_solver.solve(k, T_left, T_right, q_gen)

        # Save for backward
        ctx.save_for_backward(k_tensor)
        ctx.physics_solver = physics_solver
        ctx.T_optimal = T_optimal
        ctx.params = (T_left, T_right, q_gen)

        return torch.tensor(T_optimal, dtype=torch.float32)

    @staticmethod
    def backward(ctx, grad_output):
        """Backward: compute ∂Loss/∂k using chain rule"""
        k_tensor = ctx.saved_tensors[0]
        physics_solver = ctx.physics_solver
        T_optimal = ctx.T_optimal
        T_left, T_right, q_gen = ctx.params

        # Compute sensitivity ∂T/∂k
        dT_dk = physics_solver.compute_sensitivity(
            k_tensor.item(), T_left, T_right, q_gen, T_optimal
        )

        # Chain rule: ∂Loss/∂k = (∂Loss/∂T) · (∂T/∂k)
        grad_output_np = grad_output.detach().numpy()
        grad_k = np.dot(grad_output_np, dT_dk)

        grad_k_tensor = torch.tensor(grad_k, dtype=torch.float32)

        return grad_k_tensor, None, None, None, None


def generate_experimental_data(n_samples: int = 200, n_nodes: int = 50,
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

        if (i + 1) % 50 == 0:
            print(f"  Generated {i + 1}/{n_samples} samples")

    return exp_data


def normalize_temperature(T: np.ndarray, T_left: float, T_right: float) -> np.ndarray:
    """
    Normalize temperature profile to [0, 1] range
    This helps NN training stability
    """
    T_min = min(T_left, T_right)
    T_max = max(T_left, T_right, T.max())
    T_range = T_max - T_min

    if T_range < 1e-6:
        return np.zeros_like(T)

    return (T - T_min) / T_range


def train_fully_differentiable(nn_model: ThermalConductivityNN,
                               physics_solver: CasADiPhysicsSolver,
                               train_data: List[Dict],
                               val_data: List[Dict],
                               n_epochs: int = 100,
                               lr: float = 0.001) -> Dict:
    """Train with full temperature profile as input"""

    optimizer = optim.Adam(nn_model.parameters(), lr=lr)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode='min',
                                                     factor=0.5, patience=10, verbose=True)
    physics_layer = CasADiPhysicsLayer.apply

    history = {
        'train_loss': [],
        'val_loss': [],
        'train_k_error': [],
        'val_k_error': []
    }

    print(f"\n{'='*70}")
    print("Training with FULL Temperature Profile (50 nodes)")
    print(f"{'='*70}\n")

    best_val_loss = float('inf')
    patience_counter = 0
    early_stop_patience = 20

    for epoch in range(n_epochs):
        # ========== Training ==========
        nn_model.train()
        train_loss = 0
        train_k_error = 0

        for data in train_data:
            # Normalize temperature profile
            T_norm = normalize_temperature(
                data['T_measured'],
                data['T_left'],
                data['T_right']
            )
            T_tensor = torch.FloatTensor(T_norm).unsqueeze(0)

            # NN predicts k from full temperature profile
            k_pred_tensor = nn_model(T_tensor).squeeze()

            # Physics layer (fully differentiable!)
            T_pred_tensor = physics_layer(
                k_pred_tensor,
                physics_solver,
                data['T_left'],
                data['T_right'],
                data['q_gen']
            )

            # Loss: ONLY temperature fitting
            T_measured_tensor = torch.FloatTensor(data['T_measured'])
            loss = torch.mean((T_pred_tensor - T_measured_tensor)**2)

            # Backward
            optimizer.zero_grad()
            loss.backward()

            # Gradient clipping for stability
            torch.nn.utils.clip_grad_norm_(nn_model.parameters(), max_norm=1.0)

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
                T_norm = normalize_temperature(
                    data['T_measured'],
                    data['T_left'],
                    data['T_right']
                )
                T_tensor = torch.FloatTensor(T_norm).unsqueeze(0)
                k_pred_tensor = nn_model(T_tensor).squeeze()

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

        # Learning rate scheduling
        scheduler.step(val_loss)

        history['train_loss'].append(train_loss)
        history['val_loss'].append(val_loss)
        history['train_k_error'].append(train_k_error)
        history['val_k_error'].append(val_k_error)

        # Early stopping
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            patience_counter = 0
        else:
            patience_counter += 1

        if (epoch + 1) % 10 == 0:
            print(f"Epoch [{epoch+1}/{n_epochs}]")
            print(f"  Train - Loss: {train_loss:.4f}, k Error: {train_k_error:.2f} W/m·K")
            print(f"  Val   - Loss: {val_loss:.4f}, k Error: {val_k_error:.2f} W/m·K")

        if patience_counter >= early_stop_patience:
            print(f"\nEarly stopping at epoch {epoch+1}")
            break

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
    axes[0].set_title('Training History\n(Full Temperature Profile Input)', fontsize=13)
    axes[0].legend(fontsize=11)
    axes[0].grid(True, alpha=0.3)
    axes[0].set_yscale('log')

    axes[1].plot(history['train_k_error'], label='Train k Error', linewidth=2)
    axes[1].plot(history['val_k_error'], label='Val k Error', linewidth=2)
    axes[1].set_xlabel('Epoch', fontsize=12)
    axes[1].set_ylabel('k Error [W/m·K]', fontsize=12)
    axes[1].set_title('k Prediction Error', fontsize=13)
    axes[1].legend(fontsize=11)
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('full_profile_training_history.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("\n✓ Saved: full_profile_training_history.png")

    # Plot 2: Test predictions
    nn_model.eval()
    n_test_plots = min(4, len(test_data))
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    k_true_list = []
    k_pred_list = []

    for i in range(n_test_plots):
        data = test_data[i]

        T_norm = normalize_temperature(
            data['T_measured'],
            data['T_left'],
            data['T_right']
        )
        T_tensor = torch.FloatTensor(T_norm).unsqueeze(0)

        with torch.no_grad():
            k_pred = nn_model(T_tensor).item()

        T_pred = physics_solver.solve(
            k_pred, data['T_left'], data['T_right'], data['q_gen']
        )

        axes[i].plot(physics_solver.x, data['T_measured'], 'o',
                    label='Measured', alpha=0.6, markersize=5)
        axes[i].plot(physics_solver.x, T_pred, '-',
                    label='Predicted', linewidth=2)
        axes[i].set_xlabel('Position [m]', fontsize=11)
        axes[i].set_ylabel('Temperature [K]', fontsize=11)

        error = abs(k_pred - data['k_true'])
        error_pct = error / data['k_true'] * 100
        axes[i].set_title(
            f'Test {i+1}: k_true={data["k_true"]:.1f}, k_pred={k_pred:.1f} W/m·K\n'
            f'Error: {error:.1f} W/m·K ({error_pct:.1f}%)',
            fontsize=11
        )
        axes[i].legend(fontsize=10)
        axes[i].grid(True, alpha=0.3)

        k_true_list.append(data['k_true'])
        k_pred_list.append(k_pred)

    plt.tight_layout()
    plt.savefig('full_profile_test_predictions.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("✓ Saved: full_profile_test_predictions.png")

    # Plot 3: k prediction accuracy (all test data)
    for data in test_data[n_test_plots:]:
        T_norm = normalize_temperature(
            data['T_measured'],
            data['T_left'],
            data['T_right']
        )
        T_tensor = torch.FloatTensor(T_norm).unsqueeze(0)

        with torch.no_grad():
            k_pred = nn_model(T_tensor).item()

        k_true_list.append(data['k_true'])
        k_pred_list.append(k_pred)

    k_true_arr = np.array(k_true_list)
    k_pred_arr = np.array(k_pred_list)

    plt.figure(figsize=(9, 9))
    plt.scatter(k_true_arr, k_pred_arr, alpha=0.6, s=60, edgecolors='black', linewidths=0.5)

    k_min = min(k_true_arr.min(), k_pred_arr.min())
    k_max = max(k_true_arr.max(), k_pred_arr.max())
    plt.plot([k_min, k_max], [k_min, k_max], 'r--', lw=2, label='Perfect prediction')

    plt.xlabel('True Thermal Conductivity [W/m·K]', fontsize=13)
    plt.ylabel('Predicted Thermal Conductivity [W/m·K]', fontsize=13)
    plt.title('k Prediction Accuracy\n(Full Temperature Profile Input)', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=12)

    r2 = 1 - np.sum((k_true_arr - k_pred_arr)**2) / np.sum((k_true_arr - np.mean(k_true_arr))**2)
    mae = np.mean(np.abs(k_true_arr - k_pred_arr))
    rmse = np.sqrt(np.mean((k_true_arr - k_pred_arr)**2))
    mape = np.mean(np.abs((k_true_arr - k_pred_arr) / k_true_arr)) * 100

    textstr = f'R² = {r2:.4f}\nMAE = {mae:.2f} W/m·K\nRMSE = {rmse:.2f} W/m·K\nMAPE = {mape:.2f}%'
    plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes,
            verticalalignment='top', fontsize=12, fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))

    plt.tight_layout()
    plt.savefig('full_profile_k_prediction.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("✓ Saved: full_profile_k_prediction.png")

    # Print statistics
    print(f"\n{'='*70}")
    print("Test Set Performance:")
    print(f"{'='*70}")
    print(f"R² Score:  {r2:.4f}")
    print(f"MAE:       {mae:.2f} W/m·K")
    print(f"RMSE:      {rmse:.2f} W/m·K")
    print(f"MAPE:      {mape:.2f}%")
    print(f"{'='*70}")


def main():
    """Main demonstration"""

    print("\n" + "="*70)
    print("Fully Differentiable PINN with Full Temperature Profile")
    print("="*70)
    print("\nKey Improvement: Uses all 50 temperature nodes as NN input!")
    print("Previous: 5 summary statistics → Poor k prediction")
    print("Now:      50 temperature values → Much better k prediction")

    # Parameters
    n_nodes = 50
    n_samples = 200  # More data for better learning
    noise_level = 0.5

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

    # NN takes full temperature profile (50 nodes) as input
    nn_model = ThermalConductivityNN(n_nodes=n_nodes, hidden_sizes=[64, 32, 16])

    print(f"✓ Neural network initialized (input: {n_nodes} nodes)")

    # Step 3: Train
    print("\n[Step 3] Training fully differentiable model...")
    history = train_fully_differentiable(
        nn_model,
        physics_solver,
        train_data,
        val_data,
        n_epochs=100,
        lr=0.001
    )

    # Step 4: Evaluate
    print("\n[Step 4] Evaluating and visualizing results...")
    plot_results(history, test_data, nn_model, physics_solver)

    print("\n✓ All figures saved:")
    print("  - full_profile_training_history.png")
    print("  - full_profile_test_predictions.png")
    print("  - full_profile_k_prediction.png")
    print("\n" + "="*70)


if __name__ == "__main__":
    try:
        import casadi
        import torch
        main()
    except ImportError as e:
        print(f"Error: Required package not found - {e}")
        print("\nPlease install: pip install casadi torch matplotlib")
