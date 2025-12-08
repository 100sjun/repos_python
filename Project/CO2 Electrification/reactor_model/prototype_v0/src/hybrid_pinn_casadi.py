"""
Hybrid Physics-Informed Neural Network with CasADi
===================================================
Neural Network predicts thermal conductivity → CasADi solves physics model
→ Loss computed against experimental temperature data at all points
"""

import numpy as np
import casadi as ca
import matplotlib.pyplot as plt
from typing import Tuple, Dict, List
import torch
import torch.nn as nn
import torch.optim as optim


class ThermalConductivityNN(nn.Module):
    """Neural Network to predict thermal conductivity from input features"""

    def __init__(self, input_dim: int = 3, hidden_sizes: list = [32, 16, 8]):
        """
        Initialize neural network to predict thermal conductivity

        Parameters:
        -----------
        input_dim : int
            Number of input features (e.g., average temperature, gradient, etc.)
        hidden_sizes : list
            Hidden layer sizes
        """
        super(ThermalConductivityNN, self).__init__()

        layers = []
        prev_size = input_dim

        for hidden_size in hidden_sizes:
            layers.append(nn.Linear(prev_size, hidden_size))
            layers.append(nn.Tanh())  # Tanh for smoother gradients
            prev_size = hidden_size

        # Output layer: thermal conductivity (positive value)
        layers.append(nn.Linear(prev_size, 1))
        layers.append(nn.Softplus())  # Ensures k > 0

        self.network = nn.Sequential(*layers)

    def forward(self, x):
        """Forward pass"""
        return self.network(x)


class HybridPINNModel:
    """Hybrid model combining NN prediction with CasADi physics solver"""

    def __init__(self, length: float = 1.0, n_nodes: int = 50):
        """
        Initialize hybrid model

        Parameters:
        -----------
        length : float
            Domain length [m]
        n_nodes : int
            Number of spatial nodes
        """
        self.L = length
        self.n = n_nodes
        self.dx = length / (n_nodes - 1)
        self.x = np.linspace(0, length, n_nodes)

    def solve_physics_model(self, k: float, T_left: float, T_right: float,
                           q_gen: float = 0.0) -> np.ndarray:
        """
        Solve 1D heat conduction using CasADi with given thermal conductivity

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
        # Define temperature variables
        T = ca.MX.sym('T', self.n)

        # Objective: minimize heat equation residual
        obj = 0

        # Constraints: boundary conditions
        g = []
        g.append(T[0] - T_left)
        g.append(T[-1] - T_right)

        # Heat equation at interior nodes
        # k * d²T/dx² + q_gen = 0
        for i in range(1, self.n - 1):
            d2T_dx2 = (T[i+1] - 2*T[i] + T[i-1]) / (self.dx**2)
            residual = k * d2T_dx2 + q_gen
            obj += residual**2

        # Setup NLP
        nlp = {'x': T, 'f': obj, 'g': ca.vertcat(*g)}
        opts = {'ipopt.print_level': 0, 'print_time': 0}
        solver = ca.nlpsol('solver', 'ipopt', nlp, opts)

        # Solve
        T0 = np.linspace(T_left, T_right, self.n)
        sol = solver(x0=T0, lbg=0, ubg=0)

        return np.array(sol['x']).flatten()

    def extract_features(self, T_data: np.ndarray, T_left: float,
                        T_right: float, q_gen: float) -> np.ndarray:
        """
        Extract features from temperature data for NN input

        Parameters:
        -----------
        T_data : np.ndarray
            Measured temperature profile [n_nodes]
        T_left : float
            Left boundary temperature [K]
        T_right : float
            Right boundary temperature [K]
        q_gen : float
            Heat generation [W/m³]

        Returns:
        --------
        features : np.ndarray
            Feature vector for NN input
        """
        # Feature 1: Average temperature
        T_avg = np.mean(T_data)

        # Feature 2: Temperature gradient (approximate)
        dT_dx_avg = (T_right - T_left) / self.L

        # Feature 3: Normalized heat generation
        q_normalized = q_gen / 1000.0  # Normalize to [0, 1] range

        # Feature 4: Temperature variance
        T_var = np.var(T_data)

        # Feature 5: Max temperature
        T_max = np.max(T_data)

        features = np.array([T_avg, dT_dx_avg, q_normalized, T_var, T_max])

        return features


def generate_experimental_data(n_samples: int = 100, n_nodes: int = 50,
                               noise_level: float = 0.5) -> List[Dict]:
    """
    Generate synthetic experimental data (mimics real measurements)

    Parameters:
    -----------
    n_samples : int
        Number of experimental samples
    n_nodes : int
        Number of measurement points
    noise_level : float
        Measurement noise standard deviation [K]

    Returns:
    --------
    exp_data : list of dict
        Each dict contains: T_measured, T_left, T_right, q_gen, k_true
    """
    model = HybridPINNModel(length=1.0, n_nodes=n_nodes)
    exp_data = []

    print(f"Generating {n_samples} experimental datasets...")

    for i in range(n_samples):
        # Random true thermal conductivity
        k_true = np.random.uniform(10, 80)

        # Random boundary conditions
        T_left = np.random.uniform(350, 450)
        T_right = np.random.uniform(300, 400)

        # Random heat generation
        q_gen = np.random.uniform(0, 800)

        # Solve for "true" temperature profile
        T_true = model.solve_physics_model(k_true, T_left, T_right, q_gen)

        # Add measurement noise
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


class HybridTrainer:
    """Training loop for hybrid PINN model"""

    def __init__(self, model: HybridPINNModel, nn_model: ThermalConductivityNN):
        """
        Initialize trainer

        Parameters:
        -----------
        model : HybridPINNModel
            Physics model
        nn_model : ThermalConductivityNN
            Neural network model
        """
        self.physics_model = model
        self.nn_model = nn_model
        self.optimizer = optim.Adam(nn_model.parameters(), lr=0.01)

    def compute_loss(self, exp_data: Dict) -> Tuple[torch.Tensor, np.ndarray]:
        """
        Compute loss: MSE between physics model output and experimental data

        Parameters:
        -----------
        exp_data : dict
            Experimental data sample

        Returns:
        --------
        loss : torch.Tensor
            Total loss
        T_pred : np.ndarray
            Predicted temperature profile from physics model
        """
        # Extract features from experimental data
        features = self.physics_model.extract_features(
            exp_data['T_measured'],
            exp_data['T_left'],
            exp_data['T_right'],
            exp_data['q_gen']
        )

        # NN predicts thermal conductivity
        features_tensor = torch.FloatTensor(features).unsqueeze(0)
        k_pred = self.nn_model(features_tensor).item()

        # Solve physics model with predicted k
        T_pred = self.physics_model.solve_physics_model(
            k_pred,
            exp_data['T_left'],
            exp_data['T_right'],
            exp_data['q_gen']
        )

        # Compute loss: MSE between predicted and measured temperatures
        T_measured = exp_data['T_measured']
        loss_data = np.mean((T_pred - T_measured)**2)

        # Convert to torch tensor for backprop
        loss = torch.tensor(loss_data, requires_grad=False)

        # Add regularization: predicted k should be reasonable
        k_true = exp_data['k_true']
        loss_k = (k_pred - k_true)**2
        loss = loss + 0.01 * loss_k  # Small weight on k error

        return loss, T_pred, k_pred

    def train(self, train_data: List[Dict], val_data: List[Dict],
              n_epochs: int = 100) -> Dict:
        """
        Train the hybrid model

        Parameters:
        -----------
        train_data : list of dict
            Training data
        val_data : list of dict
            Validation data
        n_epochs : int
            Number of training epochs

        Returns:
        --------
        history : dict
            Training history
        """
        history = {
            'train_loss': [],
            'val_loss': [],
            'train_k_error': [],
            'val_k_error': []
        }

        print(f"\nTraining hybrid PINN model for {n_epochs} epochs...")

        for epoch in range(n_epochs):
            # Training
            self.nn_model.train()
            train_loss = 0
            train_k_error = 0

            for data in train_data:
                # Forward pass through NN → CasADi → loss
                features = self.physics_model.extract_features(
                    data['T_measured'], data['T_left'],
                    data['T_right'], data['q_gen']
                )
                features_tensor = torch.FloatTensor(features).unsqueeze(0)
                k_pred_tensor = self.nn_model(features_tensor)
                k_pred = k_pred_tensor.item()

                # Solve physics model
                T_pred = self.physics_model.solve_physics_model(
                    k_pred, data['T_left'], data['T_right'], data['q_gen']
                )

                # Data fitting loss
                T_measured_tensor = torch.FloatTensor(data['T_measured'])
                T_pred_tensor = torch.FloatTensor(T_pred)
                loss_data = torch.mean((T_pred_tensor - T_measured_tensor)**2)

                # Physics-informed loss (k should be reasonable)
                k_true_tensor = torch.FloatTensor([data['k_true']])
                loss_k = torch.mean((k_pred_tensor - k_true_tensor)**2)

                # Total loss
                loss = loss_data + 0.01 * loss_k

                # Backward pass
                self.optimizer.zero_grad()
                loss.backward()
                self.optimizer.step()

                train_loss += loss.item()
                train_k_error += abs(k_pred - data['k_true'])

            train_loss /= len(train_data)
            train_k_error /= len(train_data)

            # Validation
            self.nn_model.eval()
            val_loss = 0
            val_k_error = 0

            with torch.no_grad():
                for data in val_data:
                    features = self.physics_model.extract_features(
                        data['T_measured'], data['T_left'],
                        data['T_right'], data['q_gen']
                    )
                    features_tensor = torch.FloatTensor(features).unsqueeze(0)
                    k_pred = self.nn_model(features_tensor).item()

                    T_pred = self.physics_model.solve_physics_model(
                        k_pred, data['T_left'], data['T_right'], data['q_gen']
                    )

                    loss_data = np.mean((T_pred - data['T_measured'])**2)
                    val_loss += loss_data
                    val_k_error += abs(k_pred - data['k_true'])

            val_loss /= len(val_data)
            val_k_error /= len(val_data)

            history['train_loss'].append(train_loss)
            history['val_loss'].append(val_loss)
            history['train_k_error'].append(train_k_error)
            history['val_k_error'].append(val_k_error)

            if (epoch + 1) % 10 == 0:
                print(f"  Epoch [{epoch+1}/{n_epochs}]")
                print(f"    Train Loss: {train_loss:.4f}, k Error: {train_k_error:.2f}")
                print(f"    Val Loss:   {val_loss:.4f}, k Error: {val_k_error:.2f}")

        return history


def plot_results(history: Dict, test_data: List[Dict],
                physics_model: HybridPINNModel,
                nn_model: ThermalConductivityNN):
    """Plot training results and predictions"""

    # Plot 1: Training history
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    axes[0].plot(history['train_loss'], label='Train Loss')
    axes[0].plot(history['val_loss'], label='Val Loss')
    axes[0].set_xlabel('Epoch')
    axes[0].set_ylabel('MSE Loss')
    axes[0].set_title('Training History: Temperature Fitting Loss')
    axes[0].legend()
    axes[0].grid(True)

    axes[1].plot(history['train_k_error'], label='Train k Error')
    axes[1].plot(history['val_k_error'], label='Val k Error')
    axes[1].set_xlabel('Epoch')
    axes[1].set_ylabel('k Error [W/m·K]')
    axes[1].set_title('Training History: Thermal Conductivity Error')
    axes[1].legend()
    axes[1].grid(True)

    plt.tight_layout()
    plt.savefig('hybrid_training_history.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Plot 2: Test predictions
    nn_model.eval()
    n_test_plots = min(4, len(test_data))
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    k_true_list = []
    k_pred_list = []

    for i in range(n_test_plots):
        data = test_data[i]

        # Predict k
        features = physics_model.extract_features(
            data['T_measured'], data['T_left'],
            data['T_right'], data['q_gen']
        )
        features_tensor = torch.FloatTensor(features).unsqueeze(0)
        k_pred = nn_model(features_tensor).item()

        # Solve physics model
        T_pred = physics_model.solve_physics_model(
            k_pred, data['T_left'], data['T_right'], data['q_gen']
        )

        # Plot
        axes[i].plot(physics_model.x, data['T_measured'], 'o',
                    label='Measured (with noise)', alpha=0.6)
        axes[i].plot(physics_model.x, T_pred, '-',
                    label=f'Predicted (k={k_pred:.1f})', linewidth=2)
        axes[i].set_xlabel('Position [m]')
        axes[i].set_ylabel('Temperature [K]')
        axes[i].set_title(f'Test Case {i+1}: k_true={data["k_true"]:.1f}, k_pred={k_pred:.1f} W/m·K')
        axes[i].legend()
        axes[i].grid(True)

        k_true_list.append(data['k_true'])
        k_pred_list.append(k_pred)

    plt.tight_layout()
    plt.savefig('hybrid_test_predictions.png', dpi=150, bbox_inches='tight')
    plt.close()

    # Plot 3: k prediction accuracy
    for data in test_data:
        features = physics_model.extract_features(
            data['T_measured'], data['T_left'],
            data['T_right'], data['q_gen']
        )
        features_tensor = torch.FloatTensor(features).unsqueeze(0)
        k_pred = nn_model(features_tensor).item()
        k_true_list.append(data['k_true'])
        k_pred_list.append(k_pred)

    k_true_arr = np.array(k_true_list)
    k_pred_arr = np.array(k_pred_list)

    plt.figure(figsize=(8, 8))
    plt.scatter(k_true_arr, k_pred_arr, alpha=0.6)
    plt.plot([k_true_arr.min(), k_true_arr.max()],
            [k_true_arr.min(), k_true_arr.max()], 'r--', lw=2)
    plt.xlabel('True Thermal Conductivity [W/m·K]')
    plt.ylabel('Predicted Thermal Conductivity [W/m·K]')
    plt.title('Thermal Conductivity Prediction Accuracy')
    plt.grid(True)

    # Calculate R²
    r2 = 1 - np.sum((k_true_arr - k_pred_arr)**2) / np.sum((k_true_arr - np.mean(k_true_arr))**2)
    mae = np.mean(np.abs(k_true_arr - k_pred_arr))
    plt.text(0.05, 0.95, f'R² = {r2:.4f}\nMAE = {mae:.2f} W/m·K',
            transform=plt.gca().transAxes, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    plt.tight_layout()
    plt.savefig('hybrid_k_prediction.png', dpi=150, bbox_inches='tight')
    plt.close()


def main():
    """Main demonstration"""

    print("=" * 70)
    print("Hybrid Physics-Informed Neural Network with CasADi")
    print("=" * 70)

    # Parameters
    n_nodes = 50
    n_samples = 100
    noise_level = 1.0  # K

    # Step 1: Generate experimental data
    print("\n[Step 1] Generating experimental data...")
    exp_data = generate_experimental_data(n_samples, n_nodes, noise_level)

    # Split data
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
    physics_model = HybridPINNModel(length=1.0, n_nodes=n_nodes)

    # NN input: 5 features (T_avg, dT/dx, q, T_var, T_max)
    nn_model = ThermalConductivityNN(input_dim=5, hidden_sizes=[32, 16, 8])

    # Step 3: Train
    print("\n[Step 3] Training hybrid model...")
    trainer = HybridTrainer(physics_model, nn_model)
    history = trainer.train(train_data, val_data, n_epochs=100)

    # Step 4: Evaluate and visualize
    print("\n[Step 4] Evaluating and visualizing results...")
    plot_results(history, test_data, physics_model, nn_model)

    print("\n[Complete] All figures saved:")
    print("  - hybrid_training_history.png")
    print("  - hybrid_test_predictions.png")
    print("  - hybrid_k_prediction.png")


if __name__ == "__main__":
    try:
        import casadi
        import torch
        main()
    except ImportError as e:
        print(f"Error: Required package not found - {e}")
        print("\nPlease install: pip install casadi torch matplotlib")
