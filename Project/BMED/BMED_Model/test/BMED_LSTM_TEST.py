import torch
import os
import sys
import importlib.util
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score, mean_squared_error

# CUDA availability check
cuda_available = torch.cuda.is_available()
if cuda_available:
    print(f"CUDA available: {torch.cuda.get_device_name(0)}")
    device = torch.device("cuda")
else:
    print("CUDA not available, using CPU.")
    device = torch.device("cpu")

# Set model file paths
result_dir = "./result"
model_name = "bmed_model"
model_path = os.path.join(result_dir, f"{model_name}.pt")
weights_path = os.path.join(result_dir, f"{model_name}_weights.pt")
scalers_path = os.path.join(result_dir, f"{model_name}_scalers.pkl")

print("Loading model...")
print(f"Model file path: {model_path}")
print(f"Weights file path: {weights_path}")
print(f"Device: {device}")

# Import model classes from BMED_LSTM_v5.py
try:
    print("\nAttempting to import model classes from BMED_LSTM_v5.py...")
    
    # Regular import
    try:
        sys.path.append('./test')
        from BMED_LSTM_v5 import MembraneSystemModel, PhysicsLayer
        print("Import successful!")
    except ImportError:
        # Import using spec (direct file path)
        print("Regular import failed, trying spec import...")
        module_path = os.path.join('test', 'BMED_LSTM_v5.py')
        if os.path.exists(module_path):
            spec = importlib.util.spec_from_file_location("BMED_LSTM_v5", module_path)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            MembraneSystemModel = module.MembraneSystemModel
            PhysicsLayer = module.PhysicsLayer
            print(f"Spec import successful: {module_path}")
        else:
            raise ImportError(f"Module file not found: {module_path}")
    
    # Register classes in safe globals
    import torch.serialization
    torch.serialization.add_safe_globals([MembraneSystemModel, PhysicsLayer])
    print("Model classes registered in safe globals.")
    
except Exception as e:
    print(f"Model class import failed: {e}")
    print("Warning: Model loading may fail without class definitions.")

# Load model
try:
    print("\nAttempting to load model (weights_only=False)...")
    model = torch.load(model_path, map_location=device, weights_only=False)
    print("Model loaded successfully!")
    
    # Check if model is on GPU
    next_param = next(model.parameters(), None)
    if next_param is not None:
        print(f"Model parameter device: {next_param.device}")
        # Move model to device if necessary
        if next_param.device != device:
            model = model.to(device)
            print(f"Model moved to {device}.")
    
except Exception as e:
    print(f"Model loading failed: {e}")
    
    # Alternative: create model instance and load weights
    if os.path.exists(weights_path) and 'MembraneSystemModel' in globals():
        try:
            print("\nAlternative: Attempting to load from weights file...")
            # Create model instance (removed sequence_length parameter)
            model = MembraneSystemModel(lstm_units=64, lstm_layer=3, fc_units=64, fc_layer=5, time_step=0.1)
            # Load weights
            model.load_state_dict(torch.load(weights_path, map_location=device))
            model = model.to(device)
            print("Model successfully loaded from weights file!")
        except Exception as e2:
            print(f"Weight file loading failed: {e2}")
            print("All methods of model loading failed.")
    else:
        print("Weight file doesn't exist or model class not available.")

# Set model to evaluation mode if successfully loaded
if 'model' in locals() and isinstance(model, torch.nn.Module):
    model.eval()
    print("Model set to evaluation mode.")
    print("\nModel structure:")
    print(model)
else:
    print("Model loading failed.")
    exit(1)

# Load scalers
import pickle
try:
    with open(scalers_path, 'rb') as f:
        scalers = pickle.load(f)
    print("\nScalers loaded successfully!")
    print(f"Scaler keys: {scalers.keys()}")
except Exception as e:
    print(f"Scaler loading failed: {e}")
    print("Accurate predictions are difficult without scalers.")
    exit(1)

# Load data
print("\nLoading data...")
data_path = "./data/BMED_data_v6+spline.xlsx"
try:
    # Load data from raw_data sheet
    df_raw = pd.read_excel(data_path, sheet_name='raw_data')
    print(f"Data loaded: {len(df_raw)} rows")
    
    # Get unique experiment numbers
    exp_numbers = df_raw['exp'].unique()
    print(f"Found {len(exp_numbers)} experiments: {exp_numbers}")
    
    # Feature names for model input
    feature_names = ['T', 'V', 'E', 'CF_LA', 'CF_K', 'CA_LA', 'CB_K', 'VF', 'VA', 'VB']
    state_names = feature_names.copy()  # Same names for states
    mol_names = ['d_NLA', 'd_NK', 'd_VA', 'd_VB']
    
except Exception as e:
    print(f"Data loading failed: {e}")
    print("Check if the data path is correct.")
    exit(1)

# Recursive simulation function
def simulate_recursively(model, initial_state, max_time, time_step=0.1, scalers=None):
    """
    Recursively simulate the system from initial state using the model
    
    Parameters:
    -----------
    model : MembraneSystemModel
        Trained LSTM model
    initial_state : numpy.ndarray
        Initial state vector [T, V, E, CF_LA, CF_K, CA_LA, CB_K, VF, VA, VB]
    max_time : float
        Maximum simulation time (hours)
    time_step : float
        Simulation time interval (hours)
    scalers : dict
        Scalers for normalization/denormalization
    
    Returns:
    --------
    dict: Simulation results (states_history, mol_changes_history, t_history)
    """
    model.eval()
    
    # Get scaler keys
    feature_key = 'feature' if 'feature' in scalers else 'features'
    mol_key = 'mol_change' if 'mol_change' in scalers else 'mol_changes'
    state_key = 'state' if 'state' in scalers else 'states'
    
    # Normalize initial state
    initial_state_norm = scalers[feature_key].transform(initial_state.reshape(1, -1))
    
    # Convert current state to tensor (add batch and sequence dimensions)
    current_state = torch.FloatTensor(initial_state_norm).reshape(1, 1, -1).to(device)
    
    # Lists for storing results
    states_history = [initial_state.copy()]  # Store in original scale
    mol_changes_history = []
    t_history = [0.0]  # Start time recording
    current_time = 0.0
    
    # Initialize LSTM hidden state
    hidden = None
    
    # Calculate number of steps based on max_time
    steps = int(max_time / time_step) + 1
    
    # Run simulation
    with torch.no_grad():
        for step in range(steps):
            # Skip if we've exceeded max_time
            if current_time > max_time:
                break
                
            # Model prediction (maintain hidden state)
            mol_changes, new_state, hidden = model(current_state, hidden)
            
            # Denormalize results
            mol_changes_np = mol_changes.cpu().numpy().reshape(1, -1)
            new_state_np = new_state.cpu().numpy().reshape(1, -1)
            
            mol_changes_orig = scalers[mol_key].inverse_transform(mol_changes_np)[0]
            next_state_orig = scalers[state_key].inverse_transform(new_state_np)[0]
            
            # Store results
            mol_changes_history.append(mol_changes_orig)
            states_history.append(next_state_orig)
            
            # Use current prediction as input for next time step
            # Convert to normalized form
            next_state_norm = scalers[feature_key].transform(next_state_orig.reshape(1, -1))
            current_state = torch.FloatTensor(next_state_norm).reshape(1, 1, -1).to(device)
            
            # Update time
            current_time += time_step
            t_history.append(current_time)
    
    return {
        'states_history': np.array(states_history),
        'mol_changes_history': np.array(mol_changes_history),
        't_history': np.array(t_history)
    }

# Calculate metrics for evaluation
def calculate_metrics(predictions, actuals, variable_names):
    """Calculate R² and RMSE metrics for each variable"""
    metrics = {}
    
    for i, name in enumerate(variable_names):
        pred_values = predictions[:, i]
        actual_values = actuals[:, i]
        
        # Calculate metrics
        r2 = r2_score(actual_values, pred_values)
        rmse = np.sqrt(mean_squared_error(actual_values, pred_values))
        
        metrics[name] = {
            'R2': r2,
            'RMSE': rmse
        }
    
    # Calculate overall metrics
    overall_r2 = r2_score(actuals.flatten(), predictions.flatten())
    overall_rmse = np.sqrt(mean_squared_error(actuals.flatten(), predictions.flatten()))
    
    metrics['Overall'] = {
        'R2': overall_r2,
        'RMSE': overall_rmse
    }
    
    return metrics

# Run simulations for each experiment
for exp_num in exp_numbers:
    print(f"\n\n{'='*50}")
    print(f"Evaluating Experiment {exp_num}")
    print(f"{'='*50}")
    
    # Filter data for current experiment
    exp_data = df_raw[df_raw['exp'] == exp_num].sort_values('t')
    
    if len(exp_data) == 0:
        print(f"No data found for experiment {exp_num}. Skipping.")
        continue
    
    # Extract initial state and maximum time
    initial_row = exp_data.iloc[0]
    max_time = exp_data['t'].max()
    
    # Create initial feature vector
    initial_features = initial_row[feature_names].values
    
    print(f"Initial values for experiment {exp_num}:")
    for col in ['t', 'T', 'V', 'E', 'CF_LA', 'CF_K', 'CA_LA', 'CB_K', 'VF', 'VA', 'VB']:
        print(f"  {col}: {initial_row[col]}")
    
    print(f"Maximum time: {max_time} hours")
    
    # Run recursive simulation up to the last measurement time point
    print(f"\nStarting recursive simulation for experiment {exp_num}...")
    time_step = 0.1
    
    simulation_results = simulate_recursively(
        model=model,
        initial_state=initial_features,
        max_time=max_time,
        time_step=time_step,
        scalers=scalers
    )
    
    # Extract simulation results
    states_history = simulation_results['states_history']
    mol_changes_history = simulation_results['mol_changes_history']
    t_history = simulation_results['t_history']
    
    print(f"Simulation completed: {len(t_history)-1} steps, total {t_history[-1]:.1f} hours")
    print(f"States history shape: {states_history.shape}")
    
    # Create plots comparing simulation results with experimental data
    print("\nGenerating plots...")
    fig, axs = plt.subplots(5, 2, figsize=(15, 20))
    axs = axs.flatten()
    
    # Create concentration and volume plots
    for i, name in enumerate(state_names):
        ax = axs[i]
        
        # Plot simulation results
        ax.plot(t_history, states_history[:, i], 'b-', linewidth=2, label='Model Prediction')
        
        # Plot experimental data
        if name in feature_names:
            ax.scatter(exp_data['t'], exp_data[name], color='red', marker='o', label='Experimental Data')
        
        ax.set_xlabel('Time (hours)')
        ax.set_ylabel(name)
        ax.set_title(f'{name} vs Time')
        ax.grid(True, alpha=0.3)
        ax.legend()
    
    plt.tight_layout()
    plt.savefig(f'exp{exp_num}_simulation_result.png', dpi=300)
    
    # Calculate metrics by interpolating simulation results to experimental data times
    # For each time point in experimental data, find the closest simulated time point
    exp_times = exp_data['t'].values
    actual_states = exp_data[feature_names].values
    
    # Find nearest simulation time points to experimental times
    predicted_states = np.zeros_like(actual_states)
    
    for i, t in enumerate(exp_times):
        # Find closest time in simulation
        closest_idx = np.argmin(np.abs(np.array(t_history) - t))
        predicted_states[i] = states_history[closest_idx]
    
    # Calculate metrics
    metrics = calculate_metrics(predicted_states, actual_states, feature_names)
    
    # Print metrics
    print(f"\nPerformance Metrics for Experiment {exp_num}:")
    print(f"{'Variable':<10} {'R²':>10} {'RMSE':>10}")
    print("-" * 30)
    
    for var_name, var_metrics in metrics.items():
        print(f"{var_name:<10} {var_metrics['R2']:>10.4f} {var_metrics['RMSE']:>10.4f}")
    
    # Print final predicted values
    print(f"\nFinal predicted values for experiment {exp_num} at t={t_history[-1]:.1f} hours:")
    for i, name in enumerate(state_names):
        print(f"  {name}: {states_history[-1, i]:.6f}")

print("\nAll experiments evaluated successfully!")

