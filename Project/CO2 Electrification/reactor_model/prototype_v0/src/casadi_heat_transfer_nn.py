"""
1D Heat Transfer Simulation with CasADi and Neural Network
==========================================================
This module demonstrates:
1. 1D heat conduction using CasADi for optimization
2. Neural network to predict thermal conductivity from temperature data
3. Integration of physics-based model with data-driven approach
"""

import numpy as np
import casadi as ca
import matplotlib.pyplot as plt
from typing import Tuple, Dict
import torch
import torch.nn as nn
import torch.optim as optim


class HeatConduction1D:
    """1D Heat Conduction Solver using CasADi"""

    def __init__(self, length: float = 1.0, n_nodes: int = 50):
        """
        Initialize 1D heat conduction problem

        Parameters:
        -----------
        length : float
            Length of the domain [m]
        n_nodes : int
            Number of spatial nodes
        """
        self.L = length
        self.n = n_nodes
        self.dx = length / (n_nodes - 1)
        self.x = np.linspace(0, length, n_nodes)

    def solve_steady_state(self, k: float, T_left: float, T_right: float,
                          q_gen: float = 0.0) -> np.ndarray:
        """
        Solve steady-state 1D heat conduction equation using CasADi

        Equation: d/dx(k * dT/dx) + q_gen = 0

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
        # Define optimization variables (temperature at each node)
        T = ca.MX.sym('T', self.n)

        # Objective function (minimize residual of heat equation)
        obj = 0

        # Constraints list
        g = []

        # Boundary conditions
        g.append(T[0] - T_left)  # Left BC
        g.append(T[-1] - T_right)  # Right BC

        # Interior nodes - Finite difference discretization
        # d²T/dx² = -q_gen/k
        for i in range(1, self.n - 1):
            # Second derivative using central difference
            d2T_dx2 = (T[i+1] - 2*T[i] + T[i-1]) / (self.dx**2)
            # Heat equation residual
            residual = k * d2T_dx2 + q_gen
            obj += residual**2

        # Setup optimization problem
        nlp = {
            'x': T,
            'f': obj,
            'g': ca.vertcat(*g)
        }

        # Solver options
        opts = {
            'ipopt.print_level': 0,
            'print_time': 0,
            'ipopt.max_iter': 1000
        }

        solver = ca.nlpsol('solver', 'ipopt', nlp, opts)

        # Initial guess (linear interpolation)
        T0 = np.linspace(T_left, T_right, self.n)

        # Solve
        sol = solver(x0=T0, lbg=0, ubg=0)

        return np.array(sol['x']).flatten()

    def solve_transient(self, k: float, rho: float, cp: float,
                       T_init: np.ndarray, T_left: float, T_right: float,
                       t_final: float, dt: float, q_gen: float = 0.0) -> Tuple[np.ndarray, np.ndarray]:
        """
        Solve transient 1D heat conduction using explicit time integration

        Equation: rho * cp * dT/dt = d/dx(k * dT/dx) + q_gen

        Parameters:
        -----------
        k : float
            Thermal conductivity [W/m·K]
        rho : float
            Density [kg/m³]
        cp : float
            Specific heat capacity [J/kg·K]
        T_init : np.ndarray
            Initial temperature distribution [K]
        T_left : float
            Left boundary temperature [K]
        T_right : float
            Right boundary temperature [K]
        t_final : float
            Final time [s]
        dt : float
            Time step [s]
        q_gen : float
            Volumetric heat generation [W/m³]

        Returns:
        --------
        T_history : np.ndarray
            Temperature history [n_time, n_nodes]
        t_history : np.ndarray
            Time history [n_time]
        """
        # Thermal diffusivity
        alpha = k / (rho * cp)

        # Stability check (CFL condition)
        dt_crit = 0.5 * self.dx**2 / alpha
        if dt > dt_crit:
            print(f"Warning: dt={dt:.2e} > dt_crit={dt_crit:.2e}, reducing dt")
            dt = dt_crit * 0.9

        n_steps = int(t_final / dt)
        t_history = np.linspace(0, t_final, n_steps)
        T_history = np.zeros((n_steps, self.n))

        # Initial condition
        T = T_init.copy()
        T_history[0, :] = T

        # Time integration
        for step in range(1, n_steps):
            T_new = T.copy()

            # Interior nodes
            for i in range(1, self.n - 1):
                d2T_dx2 = (T[i+1] - 2*T[i] + T[i-1]) / (self.dx**2)
                dT_dt = alpha * d2T_dx2 + q_gen / (rho * cp)
                T_new[i] = T[i] + dt * dT_dt

            # Boundary conditions
            T_new[0] = T_left
            T_new[-1] = T_right

            T = T_new
            T_history[step, :] = T

        return T_history, t_history


class ThermalConductivityNN(nn.Module):
    """Neural Network to predict thermal conductivity from temperature profile"""

    def __init__(self, n_nodes: int, hidden_sizes: list = [64, 32, 16]):
        """
        Initialize neural network

        Parameters:
        -----------
        n_nodes : int
            Number of spatial nodes (input features)
        hidden_sizes : list
            Hidden layer sizes
        """
        super(ThermalConductivityNN, self).__init__()

        layers = []
        input_size = n_nodes

        for hidden_size in hidden_sizes:
            layers.append(nn.Linear(input_size, hidden_size))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(0.2))
            input_size = hidden_size

        # Output layer (single value: thermal conductivity)
        layers.append(nn.Linear(input_size, 1))

        self.network = nn.Sequential(*layers)

    def forward(self, x):
        """Forward pass"""
        return self.network(x)


def generate_training_data(n_samples: int = 1000, n_nodes: int = 50) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate training data by solving heat conduction for various thermal conductivities

    Parameters:
    -----------
    n_samples : int
        Number of samples to generate
    n_nodes : int
        Number of spatial nodes

    Returns:
    --------
    X : np.ndarray
        Temperature profiles [n_samples, n_nodes]
    y : np.ndarray
        Thermal conductivities [n_samples]
    """
    solver = HeatConduction1D(length=1.0, n_nodes=n_nodes)

    X = np.zeros((n_samples, n_nodes))
    y = np.zeros(n_samples)

    print(f"Generating {n_samples} training samples...")

    for i in range(n_samples):
        # Random thermal conductivity [1, 100] W/m·K
        k = np.random.uniform(1, 100)

        # Random boundary temperatures
        T_left = np.random.uniform(300, 500)
        T_right = np.random.uniform(300, 500)

        # Random heat generation [0, 1000] W/m³
        q_gen = np.random.uniform(0, 1000)

        # Solve heat equation
        T = solver.solve_steady_state(k, T_left, T_right, q_gen)

        X[i, :] = T
        y[i] = k

        if (i + 1) % 100 == 0:
            print(f"  Generated {i + 1}/{n_samples} samples")

    return X, y


def train_neural_network(X_train: np.ndarray, y_train: np.ndarray,
                        X_val: np.ndarray, y_val: np.ndarray,
                        n_epochs: int = 1000, batch_size: int = 32) -> ThermalConductivityNN:
    """
    Train neural network to predict thermal conductivity

    Parameters:
    -----------
    X_train : np.ndarray
        Training temperature profiles
    y_train : np.ndarray
        Training thermal conductivities
    X_val : np.ndarray
        Validation temperature profiles
    y_val : np.ndarray
        Validation thermal conductivities
    n_epochs : int
        Number of training epochs
    batch_size : int
        Batch size

    Returns:
    --------
    model : ThermalConductivityNN
        Trained model
    """
    # Convert to PyTorch tensors
    X_train_t = torch.FloatTensor(X_train)
    y_train_t = torch.FloatTensor(y_train).reshape(-1, 1)
    X_val_t = torch.FloatTensor(X_val)
    y_val_t = torch.FloatTensor(y_val).reshape(-1, 1)

    # Initialize model
    n_nodes = X_train.shape[1]
    model = ThermalConductivityNN(n_nodes)

    # Loss and optimizer
    criterion = nn.MSELoss()
    optimizer = optim.Adam(model.parameters(), lr=0.001)

    # Training loop
    train_losses = []
    val_losses = []

    print(f"\nTraining neural network for {n_epochs} epochs...")

    for epoch in range(n_epochs):
        model.train()

        # Shuffle training data
        perm = torch.randperm(X_train_t.size(0))

        epoch_loss = 0
        n_batches = 0

        for i in range(0, X_train_t.size(0), batch_size):
            indices = perm[i:i+batch_size]
            batch_X = X_train_t[indices]
            batch_y = y_train_t[indices]

            # Forward pass
            optimizer.zero_grad()
            outputs = model(batch_X)
            loss = criterion(outputs, batch_y)

            # Backward pass
            loss.backward()
            optimizer.step()

            epoch_loss += loss.item()
            n_batches += 1

        train_loss = epoch_loss / n_batches
        train_losses.append(train_loss)

        # Validation
        model.eval()
        with torch.no_grad():
            val_outputs = model(X_val_t)
            val_loss = criterion(val_outputs, y_val_t).item()
            val_losses.append(val_loss)

        if (epoch + 1) % 10 == 0:
            print(f"  Epoch [{epoch+1}/{n_epochs}], Train Loss: {train_loss:.4f}, Val Loss: {val_loss:.4f}")

    # Plot training history
    plt.figure(figsize=(10, 5))
    plt.plot(train_losses, label='Training Loss')
    plt.plot(val_losses, label='Validation Loss')
    plt.xlabel('Epoch')
    plt.ylabel('MSE Loss')
    plt.title('Training History')
    plt.legend()
    plt.grid(True)
    plt.savefig('training_history.png', dpi=150, bbox_inches='tight')
    plt.close()

    return model


def evaluate_model(model: ThermalConductivityNN, X_test: np.ndarray,
                  y_test: np.ndarray) -> Dict:
    """
    Evaluate trained model

    Parameters:
    -----------
    model : ThermalConductivityNN
        Trained model
    X_test : np.ndarray
        Test temperature profiles
    y_test : np.ndarray
        True thermal conductivities

    Returns:
    --------
    metrics : dict
        Evaluation metrics
    """
    model.eval()

    with torch.no_grad():
        X_test_t = torch.FloatTensor(X_test)
        y_pred = model(X_test_t).numpy().flatten()

    # Calculate metrics
    mse = np.mean((y_test - y_pred)**2)
    rmse = np.sqrt(mse)
    mae = np.mean(np.abs(y_test - y_pred))
    mape = np.mean(np.abs((y_test - y_pred) / y_test)) * 100
    r2 = 1 - np.sum((y_test - y_pred)**2) / np.sum((y_test - np.mean(y_test))**2)

    metrics = {
        'MSE': mse,
        'RMSE': rmse,
        'MAE': mae,
        'MAPE': mape,
        'R²': r2
    }

    print("\nModel Evaluation Metrics:")
    print(f"  MSE:  {mse:.4f}")
    print(f"  RMSE: {rmse:.4f}")
    print(f"  MAE:  {mae:.4f}")
    print(f"  MAPE: {mape:.2f}%")
    print(f"  R²:   {r2:.4f}")

    # Plot predictions vs true values
    plt.figure(figsize=(10, 5))

    plt.subplot(1, 2, 1)
    plt.scatter(y_test, y_pred, alpha=0.5)
    plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--', lw=2)
    plt.xlabel('True Thermal Conductivity [W/m·K]')
    plt.ylabel('Predicted Thermal Conductivity [W/m·K]')
    plt.title(f'Predictions vs True Values (R²={r2:.4f})')
    plt.grid(True)

    plt.subplot(1, 2, 2)
    errors = y_pred - y_test
    plt.hist(errors, bins=50, edgecolor='black')
    plt.xlabel('Prediction Error [W/m·K]')
    plt.ylabel('Frequency')
    plt.title(f'Error Distribution (MAE={mae:.4f})')
    plt.grid(True)

    plt.tight_layout()
    plt.savefig('model_evaluation.png', dpi=150, bbox_inches='tight')
    plt.close()

    return metrics


def demonstrate_example():
    """Demonstrate the complete workflow"""

    print("=" * 70)
    print("1D Heat Conduction with CasADi and Neural Network")
    print("=" * 70)

    # Parameters
    n_nodes = 50
    n_samples = 1000

    # Step 1: Generate training data
    print("\n[Step 1] Generating training data...")
    X, y = generate_training_data(n_samples=n_samples, n_nodes=n_nodes)

    # Split data
    n_train = int(0.7 * n_samples)
    n_val = int(0.15 * n_samples)

    X_train, y_train = X[:n_train], y[:n_train]
    X_val, y_val = X[n_train:n_train+n_val], y[n_train:n_train+n_val]
    X_test, y_test = X[n_train+n_val:], y[n_train+n_val:]

    print(f"\nDataset split:")
    print(f"  Training:   {len(X_train)} samples")
    print(f"  Validation: {len(X_val)} samples")
    print(f"  Test:       {len(X_test)} samples")

    # Step 2: Train neural network
    print("\n[Step 2] Training neural network...")
    model = train_neural_network(X_train, y_train, X_val, y_val, n_epochs=1000)

    # Step 3: Evaluate model
    print("\n[Step 3] Evaluating model...")
    metrics = evaluate_model(model, X_test, y_test)

    # Step 4: Test on a specific case
    print("\n[Step 4] Testing on a specific case...")
    solver = HeatConduction1D(length=1.0, n_nodes=n_nodes)

    # Known thermal conductivity
    k_true = 50.0
    T_left = 400.0
    T_right = 300.0
    q_gen = 500.0

    # Solve heat equation
    T = solver.solve_steady_state(k_true, T_left, T_right, q_gen)

    # Predict thermal conductivity using NN
    model.eval()
    with torch.no_grad():
        T_tensor = torch.FloatTensor(T).unsqueeze(0)
        k_pred = model(T_tensor).item()

    print(f"\nTest Case Results:")
    print(f"  True thermal conductivity:      {k_true:.2f} W/m·K")
    print(f"  Predicted thermal conductivity: {k_pred:.2f} W/m·K")
    print(f"  Error:                          {abs(k_true - k_pred):.2f} W/m·K ({abs(k_true - k_pred)/k_true*100:.2f}%)")

    # Plot temperature profile
    plt.figure(figsize=(10, 5))
    plt.plot(solver.x, T, 'b-', linewidth=2, label='Temperature Profile')
    plt.xlabel('Position [m]')
    plt.ylabel('Temperature [K]')
    plt.title(f'1D Heat Conduction (k={k_true:.1f} W/m·K, k_pred={k_pred:.1f} W/m·K)')
    plt.grid(True)
    plt.legend()
    plt.savefig('test_case_temperature.png', dpi=150, bbox_inches='tight')
    plt.close()

    print("\n[Complete] All figures saved to current directory")
    print("  - training_history.png")
    print("  - model_evaluation.png")
    print("  - test_case_temperature.png")


if __name__ == "__main__":
    # Check if required packages are available
    try:
        import casadi
        import torch
        demonstrate_example()
    except ImportError as e:
        print(f"Error: Required package not found - {e}")
        print("\nPlease install required packages:")
        print("  pip install casadi torch matplotlib")
