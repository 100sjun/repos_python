"""
Easier Problem for Testing Fully Differentiable PINN
====================================================
Makes the inverse problem easier by:
1. Larger temperature differences (wider T_left, T_right range)
2. Larger heat generation (more distinct patterns)
3. Wider k range (easier to distinguish)
4. Sensitivity validation
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

        T_sym = ca.MX.sym('T', self.n)
        k_sym = ca.MX.sym('k')
        T_left_sym = ca.MX.sym('T_left')
        T_right_sym = ca.MX.sym('T_right')
        q_gen_sym = ca.MX.sym('q_gen')

        obj = 0
        for i in range(1, self.n - 1):
            d2T_dx2 = (T_sym[i+1] - 2*T_sym[i] + T_sym[i-1]) / (self.dx**2)
            residual = k_sym * d2T_dx2 + q_gen_sym
            obj += residual**2

        g = ca.vertcat(T_sym[0] - T_left_sym, T_sym[-1] - T_right_sym)
        p = ca.vertcat(k_sym, T_left_sym, T_right_sym, q_gen_sym)

        nlp = {'x': T_sym, 'p': p, 'f': obj, 'g': g}
        opts = {'ipopt.print_level': 0, 'print_time': 0, 'ipopt.max_iter': 1000}

        self.solver = ca.nlpsol('solver', 'ipopt', nlp, opts)
        print("✓ Symbolic CasADi problem built")

    def solve(self, k: float, T_left: float, T_right: float, q_gen: float) -> np.ndarray:
        """Solve heat equation with IPOPT"""
        T0 = np.linspace(T_left, T_right, self.n)
        p = [k, T_left, T_right, q_gen]
        sol = self.solver(x0=T0, p=p, lbg=0, ubg=0)
        return np.array(sol['x']).flatten()

    def compute_sensitivity(self, k: float, T_left: float, T_right: float,
                           q_gen: float, T_optimal: np.ndarray) -> np.ndarray:
        """Compute sensitivity ∂T/∂k using implicit function theorem"""

        T_sym = ca.MX.sym('T', self.n)
        k_sym = ca.MX.sym('k')
        T_left_sym = ca.MX.sym('T_left')
        T_right_sym = ca.MX.sym('T_right')
        q_gen_sym = ca.MX.sym('q_gen')

        residuals = []
        for i in range(1, self.n - 1):
            d2T_dx2 = (T_sym[i+1] - 2*T_sym[i] + T_sym[i-1]) / (self.dx**2)
            residual = k_sym * d2T_dx2 + q_gen_sym
            residuals.append(residual)

        residuals.insert(0, T_sym[0] - T_left_sym)
        residuals.append(T_sym[-1] - T_right_sym)
        F = ca.vertcat(*residuals)

        dF_dT = ca.jacobian(F, T_sym)
        dF_dk = ca.jacobian(F, k_sym)

        jac_T_func = ca.Function('jac_T', [T_sym, k_sym, T_left_sym, T_right_sym, q_gen_sym], [dF_dT])
        jac_k_func = ca.Function('jac_k', [T_sym, k_sym, T_left_sym, T_right_sym, q_gen_sym], [dF_dk])

        dF_dT_val = jac_T_func(T_optimal, k, T_left, T_right, q_gen)
        dF_dk_val = jac_k_func(T_optimal, k, T_left, T_right, q_gen)

        dF_dT_dense = np.array(dF_dT_val.full())
        dF_dk_dense = np.array(dF_dk_val.full()).flatten()

        dT_dk = -np.linalg.solve(dF_dT_dense, dF_dk_dense)

        return dT_dk

    def compute_sensitivity_numerical(self, k: float, T_left: float, T_right: float,
                                     q_gen: float, epsilon: float = 0.1) -> np.ndarray:
        """Compute sensitivity numerically for validation"""
        T1 = self.solve(k, T_left, T_right, q_gen)
        T2 = self.solve(k + epsilon, T_left, T_right, q_gen)
        return (T2 - T1) / epsilon


def validate_sensitivity(physics_solver: CasADiPhysicsSolver):
    """Validate that analytical sensitivity matches numerical sensitivity"""

    print("\n" + "="*70)
    print("Sensitivity Validation Test")
    print("="*70)

    # Test cases with EASIER parameters
    test_cases = [
        {'k': 10.0, 'T_left': 600.0, 'T_right': 200.0, 'q_gen': 2000.0},
        {'k': 50.0, 'T_left': 700.0, 'T_right': 100.0, 'q_gen': 5000.0},
        {'k': 150.0, 'T_left': 800.0, 'T_right': 150.0, 'q_gen': 3000.0},
    ]

    for i, case in enumerate(test_cases):
        print(f"\nTest Case {i+1}: k={case['k']:.1f}, T_left={case['T_left']:.1f}, "
              f"T_right={case['T_right']:.1f}, q_gen={case['q_gen']:.1f}")

        # Solve for optimal temperature
        T_optimal = physics_solver.solve(case['k'], case['T_left'],
                                        case['T_right'], case['q_gen'])

        # Analytical sensitivity
        dT_dk_analytical = physics_solver.compute_sensitivity(
            case['k'], case['T_left'], case['T_right'], case['q_gen'], T_optimal
        )

        # Numerical sensitivity
        dT_dk_numerical = physics_solver.compute_sensitivity_numerical(
            case['k'], case['T_left'], case['T_right'], case['q_gen'], epsilon=0.1
        )

        # Compare
        max_error = np.max(np.abs(dT_dk_analytical - dT_dk_numerical))
        mean_error = np.mean(np.abs(dT_dk_analytical - dT_dk_numerical))
        relative_error = mean_error / (np.mean(np.abs(dT_dk_numerical)) + 1e-10) * 100

        print(f"  Max Error:      {max_error:.6f}")
        print(f"  Mean Error:     {mean_error:.6f}")
        print(f"  Relative Error: {relative_error:.2f}%")

        if max_error < 0.1:
            print("  ✓ PASS: Sensitivity calculation is accurate")
        else:
            print("  ✗ FAIL: Sensitivity calculation has significant error")

    print("\n" + "="*70)


class CasADiPhysicsLayer(torch.autograd.Function):
    """Custom PyTorch layer with gradient through CasADi"""

    @staticmethod
    def forward(ctx, k_tensor, physics_solver, T_left, T_right, q_gen):
        k = k_tensor.item()
        T_optimal = physics_solver.solve(k, T_left, T_right, q_gen)

        ctx.save_for_backward(k_tensor)
        ctx.physics_solver = physics_solver
        ctx.T_optimal = T_optimal
        ctx.params = (T_left, T_right, q_gen)

        return torch.tensor(T_optimal, dtype=torch.float32)

    @staticmethod
    def backward(ctx, grad_output):
        k_tensor = ctx.saved_tensors[0]
        physics_solver = ctx.physics_solver
        T_optimal = ctx.T_optimal
        T_left, T_right, q_gen = ctx.params

        dT_dk = physics_solver.compute_sensitivity(
            k_tensor.item(), T_left, T_right, q_gen, T_optimal
        )

        grad_output_np = grad_output.detach().numpy()
        grad_k = np.dot(grad_output_np, dT_dk)
        grad_k_tensor = torch.tensor(grad_k, dtype=torch.float32)

        return grad_k_tensor, None, None, None, None


def generate_easier_data(n_samples: int = 200, n_nodes: int = 50,
                        noise_level: float = 1.0) -> List[Dict]:
    """
    Generate EASIER experimental data with:
    - Larger temperature range (200-800K instead of 300-450K)
    - Larger heat generation (0-5000 instead of 0-800)
    - Wider k range (5-200 instead of 10-80)
    """
    physics_solver = CasADiPhysicsSolver(length=1.0, n_nodes=n_nodes)
    exp_data = []

    print(f"\nGenerating {n_samples} EASIER experimental datasets...")
    print("  Wider ranges for more distinct temperature patterns:")
    print("  - k:      5-200 W/m·K (vs 10-80)")
    print("  - T_left: 500-800 K (vs 350-450)")
    print("  - T_right: 100-300 K (vs 300-400)")
    print("  - q_gen:  0-5000 W/m³ (vs 0-800)")

    for i in range(n_samples):
        # MUCH WIDER RANGES!
        k_true = np.random.uniform(5, 200)       # 5-200 (was 10-80)
        T_left = np.random.uniform(500, 800)     # 500-800 (was 350-450)
        T_right = np.random.uniform(100, 300)    # 100-300 (was 300-400)
        q_gen = np.random.uniform(0, 5000)       # 0-5000 (was 0-800)

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
    """Normalize temperature profile to [0, 1] range"""
    T_min = min(T_left, T_right, T.min())
    T_max = max(T_left, T_right, T.max())
    T_range = T_max - T_min

    if T_range < 1e-6:
        return np.zeros_like(T)

    return (T - T_min) / T_range


def train_model(nn_model: ThermalConductivityNN,
               physics_solver: CasADiPhysicsSolver,
               train_data: List[Dict],
               val_data: List[Dict],
               n_epochs: int = 100,
               lr: float = 0.001) -> Dict:
    """Train with easier problem"""

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
    print("Training with EASIER Problem")
    print(f"{'='*70}\n")

    best_val_loss = float('inf')
    patience_counter = 0
    early_stop_patience = 20

    for epoch in range(n_epochs):
        # Training
        nn_model.train()
        train_loss = 0
        train_k_error = 0

        for data in train_data:
            T_norm = normalize_temperature(data['T_measured'], data['T_left'], data['T_right'])
            T_tensor = torch.FloatTensor(T_norm).unsqueeze(0)

            k_pred_tensor = nn_model(T_tensor).squeeze()

            T_pred_tensor = physics_layer(
                k_pred_tensor, physics_solver,
                data['T_left'], data['T_right'], data['q_gen']
            )

            T_measured_tensor = torch.FloatTensor(data['T_measured'])
            loss = torch.mean((T_pred_tensor - T_measured_tensor)**2)

            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(nn_model.parameters(), max_norm=1.0)
            optimizer.step()

            train_loss += loss.item()
            train_k_error += abs(k_pred_tensor.item() - data['k_true'])

        train_loss /= len(train_data)
        train_k_error /= len(train_data)

        # Validation
        nn_model.eval()
        val_loss = 0
        val_k_error = 0

        with torch.no_grad():
            for data in val_data:
                T_norm = normalize_temperature(data['T_measured'], data['T_left'], data['T_right'])
                T_tensor = torch.FloatTensor(T_norm).unsqueeze(0)
                k_pred_tensor = nn_model(T_tensor).squeeze()

                T_pred = physics_solver.solve(
                    k_pred_tensor.item(), data['T_left'], data['T_right'], data['q_gen']
                )

                loss = np.mean((T_pred - data['T_measured'])**2)
                val_loss += loss
                val_k_error += abs(k_pred_tensor.item() - data['k_true'])

        val_loss /= len(val_data)
        val_k_error /= len(val_data)

        scheduler.step(val_loss)

        history['train_loss'].append(train_loss)
        history['val_loss'].append(val_loss)
        history['train_k_error'].append(train_k_error)
        history['val_k_error'].append(val_k_error)

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
    """Plot results"""

    # Plot 1: Training history
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    axes[0].plot(history['train_loss'], label='Train Loss', linewidth=2)
    axes[0].plot(history['val_loss'], label='Val Loss', linewidth=2)
    axes[0].set_xlabel('Epoch', fontsize=12)
    axes[0].set_ylabel('MSE Loss', fontsize=12)
    axes[0].set_title('Training History\n(EASIER Problem)', fontsize=13, fontweight='bold')
    axes[0].legend(fontsize=11)
    axes[0].grid(True, alpha=0.3)
    axes[0].set_yscale('log')

    axes[1].plot(history['train_k_error'], label='Train k Error', linewidth=2)
    axes[1].plot(history['val_k_error'], label='Val k Error', linewidth=2)
    axes[1].set_xlabel('Epoch', fontsize=12)
    axes[1].set_ylabel('k Error [W/m·K]', fontsize=12)
    axes[1].set_title('k Prediction Error', fontsize=13, fontweight='bold')
    axes[1].legend(fontsize=11)
    axes[1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('easier_training_history.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("\n✓ Saved: easier_training_history.png")

    # Plot 2: k prediction accuracy
    nn_model.eval()
    k_true_list = []
    k_pred_list = []

    with torch.no_grad():
        for data in test_data:
            T_norm = normalize_temperature(data['T_measured'], data['T_left'], data['T_right'])
            T_tensor = torch.FloatTensor(T_norm).unsqueeze(0)
            k_pred = nn_model(T_tensor).item()

            k_true_list.append(data['k_true'])
            k_pred_list.append(k_pred)

    k_true_arr = np.array(k_true_list)
    k_pred_arr = np.array(k_pred_list)

    plt.figure(figsize=(10, 10))
    plt.scatter(k_true_arr, k_pred_arr, alpha=0.6, s=80, edgecolors='black', linewidths=0.5)

    k_min = min(k_true_arr.min(), k_pred_arr.min())
    k_max = max(k_true_arr.max(), k_pred_arr.max())
    plt.plot([k_min, k_max], [k_min, k_max], 'r--', lw=3, label='Perfect prediction')

    plt.xlabel('True Thermal Conductivity [W/m·K]', fontsize=14)
    plt.ylabel('Predicted Thermal Conductivity [W/m·K]', fontsize=14)
    plt.title('k Prediction Accuracy\n(EASIER Problem: Wider Ranges)', fontsize=15, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend(fontsize=12)

    r2 = 1 - np.sum((k_true_arr - k_pred_arr)**2) / np.sum((k_true_arr - np.mean(k_true_arr))**2)
    mae = np.mean(np.abs(k_true_arr - k_pred_arr))
    rmse = np.sqrt(np.mean((k_true_arr - k_pred_arr)**2))
    mape = np.mean(np.abs((k_true_arr - k_pred_arr) / k_true_arr)) * 100

    textstr = (f'R² = {r2:.4f}\n'
              f'MAE = {mae:.2f} W/m·K\n'
              f'RMSE = {rmse:.2f} W/m·K\n'
              f'MAPE = {mape:.2f}%')

    plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes,
            verticalalignment='top', fontsize=13, fontweight='bold',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))

    plt.tight_layout()
    plt.savefig('easier_k_prediction.png', dpi=150, bbox_inches='tight')
    plt.close()
    print("✓ Saved: easier_k_prediction.png")

    # Statistics
    print(f"\n{'='*70}")
    print("Test Set Performance (EASIER Problem):")
    print(f"{'='*70}")
    print(f"R² Score:  {r2:.4f}")
    print(f"MAE:       {mae:.2f} W/m·K")
    print(f"RMSE:      {rmse:.2f} W/m·K")
    print(f"MAPE:      {mape:.2f}%")
    print(f"{'='*70}")


def main():
    """Main test"""

    print("\n" + "="*70)
    print("EASIER Problem Test for Fully Differentiable PINN")
    print("="*70)

    # Step 1: Validate sensitivity
    print("\n[Step 1] Validating CasADi Sensitivity Calculation...")
    physics_solver = CasADiPhysicsSolver(length=1.0, n_nodes=50)
    validate_sensitivity(physics_solver)

    # Step 2: Generate easier data
    print("\n[Step 2] Generating EASIER experimental data...")
    n_samples = 200
    exp_data = generate_easier_data(n_samples, n_nodes=50, noise_level=1.0)

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

    # Step 3: Initialize and train
    print("\n[Step 3] Initializing and training model...")
    nn_model = ThermalConductivityNN(n_nodes=50, hidden_sizes=[64, 32, 16])

    history = train_model(
        nn_model, physics_solver,
        train_data, val_data,
        n_epochs=100, lr=0.001
    )

    # Step 4: Evaluate
    print("\n[Step 4] Evaluating results...")
    plot_results(history, test_data, nn_model, physics_solver)

    print("\n✓ All figures saved!")
    print("\n" + "="*70)


if __name__ == "__main__":
    try:
        import casadi
        import torch
        main()
    except ImportError as e:
        print(f"Error: {e}")
        print("\nInstall: pip install casadi torch matplotlib")
