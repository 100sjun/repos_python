import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import torch
import pickle
from sklearn.metrics import r2_score

# Get parent directory path
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.append(parent_dir)

# Check for CUDA
device = 'cuda' if torch.cuda.is_available() else 'cpu'
print(f"Using device: {device}")

def load_model_with_scalers(save_dir='result', model_name='bmed_model', device='cuda' if torch.cuda.is_available() else 'cpu'):
    """
    Function to load saved model and scalers
    
    Parameters:
    -----------
    save_dir : str
        Directory path where model is saved
    model_name : str
        Name of the saved model
    device : str
        Device to load the model to ('cuda' or 'cpu')
        
    Returns:
    --------
    model : nn.Module
        Loaded model
    scalers : dict
        Loaded scalers dictionary
    """
    # Model path
    model_path = os.path.join(save_dir, f'{model_name}.pt')
    
    # Weights path (as fallback)
    weights_path = os.path.join(save_dir, f'{model_name}_weights.pt')
    
    # Scalers path
    scalers_path = os.path.join(save_dir, f'{model_name}_scalers.pkl')
    
    # Check if files exist
    if not os.path.exists(model_path):
        raise FileNotFoundError(f"Model file not found: {model_path}")
    
    if not os.path.exists(scalers_path):
        raise FileNotFoundError(f"Scalers file not found: {scalers_path}")
    
    # Try loading the model with weights_only=False for PyTorch 2.6+ compatibility
    try:
        print(f"Loading model from: {model_path}")
        model = torch.load(model_path, map_location=device, weights_only=False)
        print("Model loaded successfully!")
    except Exception as e:
        print(f"Error loading full model: {e}")
        
        # Try loading model weights as fallback
        try:
            print(f"Trying to load model weights from: {weights_path}")
            
            # Try to import the model class from various module paths
            try:
                from test.BMED_LSTM_v5 import MembraneSystemModel
                print("Imported model class from BMED_LSTM_v5")
            except ImportError:
                try:
                    from test.BMED_LSTM_v3 import MembraneSystemModel
                    print("Imported model class from BMED_LSTM_v3")
                except ImportError:
                    from test.BMED_FF_v0 import MembraneFFNetwork as MembraneSystemModel
                    print("Imported model class from BMED_FF_v0")
            
            # Create a new model instance
            # Note: This assumes model parameters are standard; adjust as needed
            model = MembraneSystemModel(
                input_size=10,
                hidden_size=64,
                lstm_layers=2,
                fc_layers=2,
                dropout=0.2,
                device=device
            ).to(device)
            
            # Load the state dictionary
            model.load_state_dict(torch.load(weights_path, map_location=device))
            print("Model weights loaded successfully!")
        except Exception as e2:
            raise RuntimeError(f"Failed to load model weights: {e2}")
    
    # Load scalers
    try:
        print(f"Loading scalers from: {scalers_path}")
        with open(scalers_path, 'rb') as f:
            scalers = pickle.load(f)
        print("Scalers loaded successfully!")
    except Exception as e:
        raise RuntimeError(f"Failed to load scalers: {e}")
    
    # Set model to evaluation mode
    model.eval()
    
    return model, scalers

def recursive_simulation_from_initial_state(model, initial_state, max_time=6.0, time_step=0.1, scalers=None):
    """
    Function to simulate the system recursively from initial state
    
    Parameters:
    -----------
    model : nn.Module
        Trained PyTorch model for simulation
    initial_state : numpy.ndarray
        Initial state vector [T, V, E, CF_LA, CF_K, CA_LA, CB_K, VF, VA, VB]
    max_time : float
        Maximum simulation time in hours
    time_step : float
        Time step for simulation in hours
    scalers : dict
        Dictionary of scalers for normalizing/denormalizing data
        
    Returns:
    --------
    dict
        Dictionary containing simulation results
    """
    # Set model to evaluation mode
    model.eval()
    
    # Initialize variables
    current_state = initial_state.copy()
    
    # Normalize initial state
    if scalers is not None:
        # Reshape for scaler
        state_reshaped = current_state.reshape(1, -1)
        current_state = scalers['state'].transform(state_reshaped)[0]
    
    # Convert to tensor
    device = next(model.parameters()).device
    
    # Lists to store history
    states_history = []
    mol_changes_history = []
    t_history = []
    
    # Store initial state
    states_history.append(current_state)
    t_history.append(0.0)
    
    # Number of steps
    n_steps = int(max_time / time_step)
    
    # Run simulation
    print(f"Starting recursive simulation for {n_steps} steps...")
    
    for step in range(n_steps):
        # Current time
        t_current = (step + 1) * time_step
        
        # Convert current state to tensor for model input
        current_state_tensor = torch.tensor(current_state, dtype=torch.float32).unsqueeze(0).to(device)
        
        # Forward pass - predict mole changes
        with torch.no_grad():
            mol_changes = model(current_state_tensor)
            
        # Convert to numpy array
        mol_changes = mol_changes.cpu().numpy()[0]
        
        # Store mole changes
        mol_changes_history.append(mol_changes)
        
        # Update state using mole changes and physical model
        if scalers is not None:
            # Denormalize mole changes
            mol_changes_denorm = scalers['mol_change'].inverse_transform(mol_changes.reshape(1, -1))[0]
            
            # Extract current state in original scale
            state_denorm = scalers['state'].inverse_transform(current_state.reshape(1, -1))[0]
            
            # Unpack state values
            T, V, E, CF_LA, CF_K, CA_LA, CB_K, VF, VA, VB = state_denorm
            
            # Unpack mole changes (rates per hour)
            d_NLA, d_NK, d_VA, d_VB = mol_changes_denorm
            
            # Calculate mole changes for this time step
            d_NLA_step = d_NLA * time_step
            d_NK_step = d_NK * time_step
            d_VA_step = d_VA * time_step
            d_VB_step = d_VB * time_step
            
            # Update volumes
            VF_new = VF
            VA_new = VA + d_VA_step
            VB_new = VB + d_VB_step
            
            # Update concentrations (C = N/V)
            # Feed compartment - constant
            # Acid compartment - LA concentration
            if VA_new > 0:
                CA_LA_new = (CA_LA * VA + d_NLA_step) / VA_new
            else:
                CA_LA_new = CA_LA
            
            # Base compartment - K concentration
            if VB_new > 0:
                CB_K_new = (CB_K * VB + d_NK_step) / VB_new
            else:
                CB_K_new = CB_K
                
            # Create updated state
            updated_state_denorm = np.array([
                T, V, E, CF_LA, CF_K, CA_LA_new, CB_K_new, VF_new, VA_new, VB_new
            ])
            
            # Normalize updated state for next iteration
            updated_state = scalers['state'].transform(updated_state_denorm.reshape(1, -1))[0]
        else:
            # Without scalers, use denormalized values directly
            updated_state = current_state
        
        # Store updated state and time
        states_history.append(updated_state)
        t_history.append(t_current)
        
        # Update current state for next iteration
        current_state = updated_state
    
    # Convert history lists to numpy arrays
    states_history = np.array(states_history)
    mol_changes_history = np.array(mol_changes_history)
    t_history = np.array(t_history)
    
    # Denormalize results if scalers are provided
    if scalers is not None:
        states_history_denorm = scalers['state'].inverse_transform(states_history)
        mol_changes_history_denorm = scalers['mol_change'].inverse_transform(mol_changes_history)
    else:
        states_history_denorm = states_history
        mol_changes_history_denorm = mol_changes_history
    
    # Return simulation results
    simulation_results = {
        'states': states_history_denorm,
        'mol_changes': mol_changes_history_denorm,
        'time': t_history
    }
    
    return simulation_results

def compare_with_experiment(exp_id=1, data_path=None, model_dir=None, save_dir=None):
    """
    Function to compare recursive simulation with experimental data
    
    Parameters:
    -----------
    exp_id : int
        Experiment ID to compare with
    data_path : str
        Path to the Excel file with experimental data
    model_dir : str
        Directory containing the saved model
    save_dir : str
        Directory to save the comparison results
        
    Returns:
    --------
    dict
        Dictionary containing simulation results
    """
    # Load data
    print(f"Loading data from: {data_path}")
    df = pd.read_excel(data_path, sheet_name='raw_data')
    
    # Filter data for specified experiment
    exp_data = df[df['exp'] == exp_id].copy()
    
    if len(exp_data) == 0:
        raise ValueError(f"No data found for experiment {exp_id}")
    
    print(f"Found {len(exp_data)} data points for experiment {exp_id}")
    
    # Load model and scalers
    print(f"Loading model from: {model_dir}")
    model, scalers = load_model_with_scalers(save_dir=model_dir, device=device)
    
    # Get variable names
    var_names = ['T', 'V', 'E', 'CF_LA', 'CF_K', 'CA_LA', 'CB_K', 'VF', 'VA', 'VB']
    mol_change_names = ['d_NLA', 'd_NK', 'd_VA', 'd_VB']
    
    # Extract initial state
    initial_state = exp_data.iloc[0][var_names].values
    
    # Get max time from experimental data
    max_time = exp_data['t'].max() + 0.1  # Add a small margin
    
    # Run recursive simulation
    print(f"Running recursive simulation for experiment {exp_id}...")
    simulation_results = recursive_simulation_from_initial_state(
        model=model,
        initial_state=initial_state,
        max_time=max_time,
        time_step=0.1,
        scalers=scalers
    )
    
    # Extract simulation results
    states_history = simulation_results['states']
    mol_changes_history = simulation_results['mol_changes']
    t_history = simulation_results['time']
    
    # Plot state variables
    print("Generating state variables comparison plot...")
    fig, axes = plt.subplots(3, 4, figsize=(20, 12))
    axes = axes.flatten()
    
    # Variables to plot (all 10 state variables)
    for i, var in enumerate(var_names):
        ax = axes[i]
        
        # Plot experimental data
        ax.plot(exp_data['t'], exp_data[var], 'ro', label='Experiment')
        
        # Plot simulation result
        ax.plot(t_history, states_history[:, i], 'b-', label='Simulation')
        
        ax.set_xlabel('Time [h]')
        ax.set_ylabel(var)
        ax.set_title(f'{var} - Exp {exp_id}')
        ax.grid(True, alpha=0.3)
        ax.legend()
    
    # Remove empty subplots
    for i in range(len(var_names), len(axes)):
        fig.delaxes(axes[i])
    
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, f'recursive_simulation_exp{exp_id}.png'), dpi=300)
    plt.close()
    
    # Plot mole change variables
    print("Generating mole changes plot...")
    
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    axes = axes.flatten()
    
    for i, name in enumerate(mol_change_names):
        ax = axes[i]
        ax.plot(t_history[:-1], mol_changes_history[:, i], 'g-', linewidth=2)
        ax.set_xlabel('Time [h]')
        ax.set_ylabel(name)
        ax.set_title(f'{name} - Simulation')
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, f'mol_changes_exp{exp_id}.png'), dpi=300)
    plt.close()
    
    # 3. Calculate R² scores
    # Find indices of simulation times closest to experimental times
    exp_times = exp_data['t'].values
    sim_indices = []
    
    for t in exp_times:
        idx = np.argmin(np.abs(t_history - t))
        sim_indices.append(idx)
    
    # Get selected simulation data points
    selected_sim_data = states_history[sim_indices]
    
    # Calculate R² scores
    r2_scores = {}
    for i, var in enumerate(var_names):
        r2 = r2_score(exp_data[var].values, selected_sim_data[:, i])
        r2_scores[var] = r2
    
    # Print R² scores
    print("\nRecursive simulation R² scores:")
    for var, score in r2_scores.items():
        print(f"  {var}: {score:.4f}")
    
    # 4. Print final values
    print("\nSimulation final values:")
    for i, var in enumerate(var_names):
        print(f"  {var}: {states_history[-1, i]:.6f}")
    
    return simulation_results

def main():
    """Main function"""
    # Set paths
    save_dir = os.path.join(parent_dir, 'result_after')
    data_path = os.path.join(parent_dir, 'data', 'BMED_data_v6+spline.xlsx')
    model_dir = os.path.join(parent_dir, 'result')
    
    # Experiment number (default: 1)
    exp_id = 1
    
    # Create results directory
    os.makedirs(save_dir, exist_ok=True)
    
    # Run recursive simulation and compare with experimental data
    simulation_results = compare_with_experiment(
        exp_id=exp_id,
        data_path=data_path,
        model_dir=model_dir,
        save_dir=save_dir
    )
    
    print(f"\nAll results saved to '{save_dir}' directory.")

if __name__ == "__main__":
    main() 