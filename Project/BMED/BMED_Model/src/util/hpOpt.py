import optuna
from optuna.samplers import GPSampler
from src.model.model_with_kfold import model_with_kfold
from datetime import datetime
import os
'''optuna를 활용하여 Feedforward Neural Network의 하이퍼파라미터를 최적화'''

class hpOpt:
    '''Hyperparameter optimization with gaussian process'''

    def __init__(self, n_trials, mode=1, datapath='./data/BMED_data_v6+spline.xlsx', hp_ranges=None, result_dir=None):
        self.n_trials = n_trials
        self.mode = mode
        self.datapath = datapath
        self.hp_ranges = hp_ranges
        self.storage_name = f'sqlite:///results/hpOpt.db'
        self.result_dir = os.path.join(result_dir, 'hpOpt')
        os.makedirs(self.result_dir, exist_ok=True)


        if self.hp_ranges is None:
            self.hp_ranges = {
                'hidden_nodes': (16, 512),
                'hidden_layers': (3, 20),
                'lr': (1e-5, 1e-1),
                'rstop': (0.05, 0.3),
                'weight_decay': (1e-5, 1e-1),
                'epochs': (50, 1000)
            }
        else:
            self.hp_ranges = hp_ranges


    def objective(self, trial):
        # set the hyperparameters
        hidden_nodes = trial.suggest_int('hidden_nodes', self.hp_ranges['hidden_nodes'][0], self.hp_ranges['hidden_nodes'][1])
        hidden_layers = trial.suggest_int('hidden_layers', self.hp_ranges['hidden_layers'][0], self.hp_ranges['hidden_layers'][1])
        lr = trial.suggest_float('lr', self.hp_ranges['lr'][0], self.hp_ranges['lr'][1], log=True)
        rstop = trial.suggest_float('rstop', self.hp_ranges['rstop'][0], self.hp_ranges['rstop'][1])
        weight_decay = trial.suggest_float('weight_decay', self.hp_ranges['weight_decay'][0], self.hp_ranges['weight_decay'][1])
        epochs = trial.suggest_int('epochs', self.hp_ranges['epochs'][0], self.hp_ranges['epochs'][1])

        # run 5-fold cross validation
        result = model_with_kfold(self.mode, self.datapath, hidden_nodes, hidden_layers, epochs, lr, rstop, weight_decay)

        avg_min_r2 = result['avg_min_r2']
        return avg_min_r2

    def optimize(self):
        # generate study with gaussian process sampler
        sampler = GPSampler(n_startup_trials=10, seed=42)
        study = optuna.create_study(
            study_name = 'hpOpt',
            storage = self.storage_name,
            direction='maximize',
            sampler=sampler,
            load_if_exists = True)

        # optimize the hyperparameters
        study.optimize(self.objective, n_trials=self.n_trials)

        # return the best hyperparameters
        best_params = study.best_trial.params
        best_value = study.best_value

        # optimize the model with the best hyperparameters
        results = {
            'best_params': best_params,
            'best_value': best_value,
            'study': study
        }

        return results
    
    def save_json(self, results):
        '''return json file from optimization results'''
        save_data = {
            'timestamp': datetime.now().strftime('%Y%m%d%H%M%S'),
            'best_params': results['best_params'],
            'best_r2_score': float(results['best_value']),
            'optimization_settings': {
                'n_trials': self.n_trials,
                'hp_ranges': self.hp_ranges
            }
        }

        filename = f'hpOpt_{save_data["timestamp"]}.json'

        save_path = os.path.join(self.result_dir, filename)

        return save_data, save_path
            

