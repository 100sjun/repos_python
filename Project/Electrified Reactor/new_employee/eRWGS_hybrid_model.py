"""
First-Principle Hybrid Neural Network Model for eRWGS Reactor
- Input: Fv_st (1 variable)
- Output: Q[nx] + k0, k1, k2 (205 variables total)
- Physics: 1D RWGS reactor model integrated with neural network
"""

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from typing import Tuple, Dict, Any

class eRWGSPhysicsModel:
    """Physical model of eRWGS reactor"""
    
    def __init__(self):
        # Global Parameters
        self.R = 0.082057  # L·atm/(mol·K)
        self.cpCO2 = 37.1
        self.cpH2 = 28.8
        self.cpAr = 35.3
        self.cpCO = 29.1
        self.cpH2O = 33.6
        self.cpCH4 = 20.8
        
        # Reactor Configurations
        self.L = 40  # cm
        self.ID = 0.052  # cm
        self.A = np.pi * (self.ID/2)**2  # cm^2
        self.nx = 100
        self.z_span = (0, self.L)
        self.z_eval = np.linspace(0, self.L, self.nx)
        
        # Feed Conditions (fixed except Fv_st)
        self.P = 1  # atm
        self.Ti = 473.15  # K
        self.zi_CO2 = 0.24
        self.zi_H2 = 0.72
        self.zi_Ar = 0.04
        self.zi_CO = 0
        self.zi_H2O = 0
        self.zi_CH4 = 0
        
        # Reaction Parameters (fixed)
        self.Ea0 = 110000
        self.Ea1 = 97100
        self.Ea2 = 77820
        self.Hr0 = 41000
        self.Hr1 = -165000
        self.Hr2 = 247300
        
        # Base values for neural network correction factors
        self.k0_base = 100000000  # Base kinetic constant
        self.k1_base = 100000000
        self.k2_base = 100000000
        self.Q_base = 10000  # Base power density W/cm³
    
    def calculate_initial_conditions(self, Fv_st: float) -> np.ndarray:
        """Calculate initial conditions based on Fv_st"""
        Fi = self.P * Fv_st/1000/60 / self.R / 298.15  # mol/s
        rhoi = self.P / self.R / self.Ti  # mol/L
        
        # Initial concentrations
        Ci_CO2 = rhoi * self.zi_CO2
        Ci_H2 = rhoi * self.zi_H2
        Ci_Ar = rhoi * self.zi_Ar
        Ci_CO = rhoi * self.zi_CO
        Ci_H2O = rhoi * self.zi_H2O
        Ci_CH4 = rhoi * self.zi_CH4
        
        # Initial velocity
        Fvi = Fi / rhoi * 1000  # mL/s
        Vi = Fvi / self.A  # cm/s
        
        return np.array([Vi, Ci_CO2, Ci_H2, Ci_CO, Ci_CH4, Ci_H2O, Ci_Ar, self.Ti])
    
    def pfr_ode(self, z: float, y: np.ndarray, Q_factors: np.ndarray, k_factors: np.ndarray) -> np.ndarray:
        """PFR ODE system with neural network correction factors"""
        V, CO2, H2, CO, CH4, H2O, Ar, T = y
        
        # Total concentration
        total_conc = CO2 + H2 + CO + CH4 + H2O + Ar
        
        # Mole fractions
        zCO2 = CO2 / total_conc
        zH2 = H2 / total_conc
        zCO = CO / total_conc
        zCH4 = CH4 / total_conc
        zH2O = H2O / total_conc
        zAr = Ar / total_conc
        
        # Partial pressures
        pCO2 = zCO2 * self.P
        pH2 = zH2 * self.P
        pCO = zCO * self.P
        pCH4 = zCH4 * self.P
        pH2O = zH2O * self.P
        pAr = zAr * self.P
        
        # Heat capacity
        cp = (zCO2*self.cpCO2 + zH2*self.cpH2 + zCO*self.cpCO + 
              zCH4*self.cpCH4 + zH2O*self.cpH2O + zAr*self.cpAr)
        
        # Reaction rates with neural network correction factors
        k0_corrected = self.k0_base * k_factors[0]
        k1_corrected = self.k1_base * k_factors[1] 
        k2_corrected = self.k2_base * k_factors[2]
        
        r0 = k0_corrected * np.exp(-self.Ea0/8.314/T) * pCO2 * pH2        # CO2 + H2 → CO + H2O
        r1 = k1_corrected * np.exp(-self.Ea1/8.314/T) * pCO2 * pH2**4     # CO2 + 4H2 → CH4 + 2H2O
        r2 = k2_corrected * np.exp(-self.Ea2/8.314/T) * pCH4 * pCO2       # CH4 + CO2 → 2CO + 2H2
        
        # Ideal gas density
        rho = self.P / (self.R * T)
        
        # Total mole change rate
        total_mole_change_rate = (-2*r1 + 2*r2) / V
        
        # Heat generation with neural network correction factor
        Q_corrected = self.Q_base * Q_factors[int(z)]
        
        # Density change due to temperature
        drho_dT = -self.P / (self.R * T**2)
        dTdz = (Q_corrected - self.Hr0*r0 - self.Hr1*r1 - self.Hr2*r2) / (rho * cp * V)
        
        # Total density change
        drhodz = total_mole_change_rate + drho_dT * dTdz
        
        # Continuity equation: velocity change
        dVdz = -(V / rho) * drhodz
        
        # Mass balance equations with convection
        convection_term = (1/V) * dVdz
        
        dCO2dz = (-r0 - r1 - r2) / V - CO2 * convection_term
        dH2dz = (-r0 - 4*r1 + 2*r2) / V - H2 * convection_term
        dCOdz = (r0 + 2*r2) / V - CO * convection_term
        dCH4dz = (r1 - r2) / V - CH4 * convection_term
        dH2Odz = (r0 + 2*r1) / V - H2O * convection_term
        dArdz = 0 - Ar * convection_term
        
        return np.array([dVdz, dCO2dz, dH2dz, dCOdz, dCH4dz, dH2Odz, dArdz, dTdz])
    
    def solve_reactor(self, Fv_st: float, Q_factors: np.ndarray, k_factors: np.ndarray) -> Dict[str, Any]:
        """Solve reactor with given correction factors"""
        var0 = self.calculate_initial_conditions(Fv_st)
        
        def ode_func(z, y):
            return self.pfr_ode(z, y, Q_factors, k_factors)
        
        try:
            sol = solve_ivp(ode_func, self.z_span, var0, t_eval=self.z_eval, method='RK45')
            
            if sol.success:
                return {
                    'success': True,
                    'z': sol.t,
                    'V': sol.y[0],
                    'CO2': sol.y[1],
                    'H2': sol.y[2],
                    'CO': sol.y[3],
                    'CH4': sol.y[4],
                    'H2O': sol.y[5],
                    'Ar': sol.y[6],
                    'T': sol.y[7]
                }
            else:
                return {'success': False, 'message': 'ODE solver failed'}
        except Exception as e:
            return {'success': False, 'message': f'Error: {str(e)}'}


class HybridNeuralNetwork(nn.Module):
    """
    Hybrid Neural Network for eRWGS reactor
    Input: Fv_st (1D)
    Output: Q_factors[nx] + k_factors[3] (correction factors)
    """
    
    def __init__(self, nx: int = 100):
        super(HybridNeuralNetwork, self).__init__()
        self.nx = nx
        
        # Neural network architecture
        # Input layer
        self.input_layer = nn.Linear(1, 64)
        
        # Hidden layers
        self.hidden1 = nn.Linear(64, 128)
        self.hidden2 = nn.Linear(128, 256)
        self.hidden3 = nn.Linear(256, 512)
        self.hidden4 = nn.Linear(512, 256)
        self.hidden5 = nn.Linear(256, 128)
        
        # Output branches
        # Q correction factors branch (nx outputs for spatial distribution)
        self.q_factors_branch = nn.Sequential(
            nn.Linear(128, 64),
            nn.ReLU(),
            nn.Linear(64, nx),
            nn.Softplus()  # Ensure positive correction factors (typically 0.1 to 10)
        )
        
        # Kinetic correction factors branch (k0, k1, k2)
        self.k_factors_branch = nn.Sequential(
            nn.Linear(128, 32),
            nn.ReLU(),
            nn.Linear(32, 16),
            nn.ReLU(),
            nn.Linear(16, 3),
            nn.Softplus()  # Ensure positive correction factors (typically 0.1 to 10)
        )
        
        # Activation functions
        self.relu = nn.ReLU()
        self.dropout = nn.Dropout(0.1)
        
    def forward(self, x: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        Forward pass
        Args:
            x: Input tensor (batch_size, 1) - Fv_st values
        Returns:
            Q_factors: Heat correction factors tensor (batch_size, nx)
            k_factors: Kinetic correction factors tensor (batch_size, 3) - [k0_factor, k1_factor, k2_factor]
        """
        # Input processing
        x = self.relu(self.input_layer(x))
        
        # Hidden layers with residual connections
        h1 = self.relu(self.hidden1(x))
        h1 = self.dropout(h1)
        
        h2 = self.relu(self.hidden2(h1))
        h2 = self.dropout(h2)
        
        h3 = self.relu(self.hidden3(h2))
        h3 = self.dropout(h3)
        
        h4 = self.relu(self.hidden4(h3))
        h4 = self.dropout(h4)
        
        h5 = self.relu(self.hidden5(h4))
        
        # Output branches
        Q_factors = self.q_factors_branch(h5)  # (batch_size, nx)
        k_factors = self.k_factors_branch(h5)  # (batch_size, 3)
        
        return Q_factors, k_factors


class HybridTrainer:
    """Training class for the hybrid model"""
    
    def __init__(self, model: HybridNeuralNetwork, physics_model: eRWGSPhysicsModel):
        self.model = model
        self.physics_model = physics_model
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.model.to(self.device)
        
    def physics_loss(self, Fv_st: torch.Tensor, Q_factors: torch.Tensor, k_factors: torch.Tensor, 
                    targets: Dict[str, torch.Tensor]) -> torch.Tensor:
        """
        Calculate physics-based loss by running the physical model
        Args:
            Fv_st: Flow rate inputs (batch_size, 1)
            Q_factors: Predicted heat correction factors (batch_size, nx)
            k_factors: Predicted kinetic correction factors (batch_size, 3)
            targets: Target concentrations/temperatures
        """
        total_loss = 0.0
        batch_size = Fv_st.shape[0]
        
        for i in range(batch_size):
            fv_val = Fv_st[i, 0].item()
            Q_factors_val = Q_factors[i].detach().cpu().numpy()
            k_factors_val = k_factors[i].detach().cpu().numpy()
            
            # Run physics model with correction factors
            result = self.physics_model.solve_reactor(fv_val, Q_factors_val, k_factors_val)
            
            if result['success']:
                # Calculate loss based on final concentrations and temperature
                pred_CO2_final = result['CO2'][-1]
                pred_CO_final = result['CO'][-1]
                pred_T_final = result['T'][-1]
                
                # Convert to tensors
                pred_tensor = torch.tensor([pred_CO2_final, pred_CO_final, pred_T_final], 
                                         device=self.device, dtype=torch.float32)
                target_tensor = targets['final_state'][i]
                
                # MSE loss
                loss = nn.MSELoss()(pred_tensor, target_tensor)
                total_loss += loss
            else:
                # Penalty for failed solutions
                total_loss += 1000.0
        
        return total_loss / batch_size
    
    def train_step(self, Fv_st_batch: torch.Tensor, targets: Dict[str, torch.Tensor], 
                   optimizer: torch.optim.Optimizer) -> float:
        """Single training step"""
        self.model.train()
        optimizer.zero_grad()
        
        # Forward pass
        Q_factors_pred, k_factors_pred = self.model(Fv_st_batch)
        
        # Calculate physics-based loss
        loss = self.physics_loss(Fv_st_batch, Q_factors_pred, k_factors_pred, targets)
        
        # Add regularization
        reg_loss = 0.001 * sum(p.pow(2.0).sum() for p in self.model.parameters())
        total_loss = loss + reg_loss
        
        # Backward pass
        total_loss.backward()
        optimizer.step()
        
        return total_loss.item()


def create_synthetic_data(n_samples: int = 100) -> Tuple[np.ndarray, Dict[str, np.ndarray]]:
    """
    Create synthetic training data
    This would typically be replaced with experimental data
    """
    # Generate Fv_st values
    Fv_st_values = 89.4  # mL/min
    
    # Create synthetic targets (would come from experiments)
    # For now, use some reasonable values
    targets = {
        'final_state': np.random.rand(n_samples, 3)  # [CO2_final, CO_final, T_final]
    }
    
    # Normalize targets to reasonable ranges
    targets['final_state'][:, 0] *= 0.1  # CO2 concentration
    targets['final_state'][:, 1] *= 0.05  # CO concentration  
    targets['final_state'][:, 2] = targets['final_state'][:, 2] * 100 + 500  # Temperature
    
    return Fv_st_values, targets


def main():
    """Main training loop"""
    print("Initializing First-Principle Hybrid Neural Network for eRWGS...")
    
    # Initialize models
    physics_model = eRWGSPhysicsModel()
    neural_model = HybridNeuralNetwork(nx=100)
    trainer = HybridTrainer(neural_model, physics_model)
    
    # Create synthetic data (replace with real data)
    Fv_st_data, target_data = create_synthetic_data(200)
    
    # Convert to tensors and move to device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    Fv_st_tensor = torch.tensor(Fv_st_data.reshape(-1, 1), dtype=torch.float32).to(device)
    target_tensor = torch.tensor(target_data['final_state'], dtype=torch.float32).to(device)
    
    # Setup training
    optimizer = optim.Adam(neural_model.parameters(), lr=0.001)
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(optimizer, patience=10)
    
    # Training loop
    n_epochs = 100
    batch_size = 16
    
    print(f"Starting training for {n_epochs} epochs...")
    
    for epoch in range(n_epochs):
        epoch_loss = 0.0
        n_batches = 0
        
        # Batch training
        for i in range(0, len(Fv_st_tensor), batch_size):
            batch_Fv_st = Fv_st_tensor[i:i+batch_size].to(trainer.device)
            batch_targets = {'final_state': target_tensor[i:i+batch_size].to(trainer.device)}
            
            loss = trainer.train_step(batch_Fv_st, batch_targets, optimizer)
            epoch_loss += loss
            n_batches += 1
        
        avg_loss = epoch_loss / n_batches
        scheduler.step(avg_loss)
        
        if (epoch + 1) % 10 == 0:
            print(f"Epoch {epoch+1}/{n_epochs}, Loss: {avg_loss:.6f}")
    
    print("Training completed!")
    
    # Test the model
    test_Fv_st = 89.4  # Original value from eRWGS notebook
    test_input = torch.tensor([[test_Fv_st]], dtype=torch.float32).to(trainer.device)
    
    neural_model.eval()
    with torch.no_grad():
        Q_factors_pred, k_factors_pred = neural_model(test_input)
        
    print(f"\nTest Results for Fv_st = {test_Fv_st} mL/min:")
    print(f"Predicted kinetic correction factors:")
    print(f"k0_factor: {k_factors_pred[0, 0].item():.3f} (k0 = {physics_model.k0_base * k_factors_pred[0, 0].item():.2e})")
    print(f"k1_factor: {k_factors_pred[0, 1].item():.3f} (k1 = {physics_model.k1_base * k_factors_pred[0, 1].item():.2e})")
    print(f"k2_factor: {k_factors_pred[0, 2].item():.3f} (k2 = {physics_model.k2_base * k_factors_pred[0, 2].item():.2e})")
    print(f"Q correction factors shape: {Q_factors_pred.shape}")
    print(f"Q_factor mean: {Q_factors_pred[0].mean().item():.3f} (Q_base = {physics_model.Q_base:.2f} W/cm³)")
    print(f"Q_factor std: {Q_factors_pred[0].std().item():.3f}")
    
    # Save model
    torch.save(neural_model.state_dict(), 'eRWGS_hybrid_model.pth')
    print("\nModel saved as 'eRWGS_hybrid_model.pth'")


if __name__ == "__main__":
    main()