import optuna
from optuna.samplers import GPSampler
from src.model.model_with_kfold import model_with_kfold
'''optuna를 활용하여 Feedforward Neural Network의 하이퍼파라미터를 최적화'''

class hpOpt:
    '''Hyperparameter optimization with gaussian process'''

    def __init__(self, n_trials, mode=1, datapath='./data/BMED_data_v6+spline.xlsx'):
        self.n_trials = n_trials
        self.mode = mode
        self.datapath = datapath

    def objective(self, trial):
        # set the hyperparameters
        hidden_nodes = trial.suggest_int('hidden_nodes', 16, 512)
        hidden_layers = trial.suggest_int('hidden_layers',3, 20)
        lr = trial.suggest_loguniform('lr', 1e-5, 1e-1, log=True)
        rstop = trial.suggest_float('rstop', 0.05, 0.3)
        weight_decay = trial.suggest_loguniform('weight_decay', 1e-5, 1e-1, log=True)
        epochs = trial.suggest_int('epochs', 50, 1000)

        # run 5-fold cross validation
        result = model_with_kfold(self.mode, self.datapath, hidden_nodes, hidden_layers, epochs, lr, rstop, weight_decay)

        avg_min_r2 = result['avg_min_r2']
        return -avg_min_r2  # Optuna는 기본적으로 최소화하므로 음수로 반환


