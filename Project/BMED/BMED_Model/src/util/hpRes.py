import optuna
import matplotlib.pyplot as plt
import os
import json
import datetime as dt

'''
학습된 hyperparameter 결과를 불러와서 결과를 확인하는 코드
'''

class hpRes:
    def __init__(self, res_dir):
        self.storage = 'sqlite:///results/hpOpt.db'
        self.study_name = 'hpOpt'
        self.res_dir = os.path.join(res_dir, 'hpOpt')
        os.makedirs(self.res_dir, exist_ok=True)

    def load_study(self):
        '''Load Optuna study'''
        self.study = optuna.load_study(study_name=self.study_name, storage=self.storage)
        return self.study
    
    def plot_history(self, image_path = 'study_history.png'):
        '''Plot history of Optuna study'''
        study = self.load_study()
        plt.figure(figsize=(10, 6))
        optuna.visualization.matplotlib.plot_optimization_history(study)
        plt.title('History of Hyperparameter Optimization')
        plt.xlabel('Number of Trials')
        plt.ylabel('Objective Value (r2 score)')
        plt.ylim(0, 1)
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(self.res_dir, image_path))
        plt.close()

    def plot_param_importance(self, image_path = 'param_importance.png'):
        '''Plot importance of hyperparameters'''
        study = self.load_study()
        plt.figure(figsize=(10, 6))
        optuna.visualization.matplotlib.plot_param_importances(study)
        plt.title('Importance of Hyperparameters')
        plt.xlabel('Importance')
        plt.ylabel('Hyperparameter')
        plt.tight_layout()
        plt.savefig(os.path.join(self.res_dir, image_path))
        plt.close()

    def save_best_trial_json(self, json_path = 'best_trial.json'):
        '''Save best trial to json file'''
        study = self.load_study()
        best_trial = study.best_trial

        # json config
        save_data = {
            'timestamp': dt.datetime.now().strftime('%Y%m%d%H%M%S'),
            'trial_number': best_trial.number,
            'best_value': best_trial.value,
            'best_params': best_trial.params,
            'study_settings': {
                'n_trials': len(study.trials),
                'datetime_start': study.trials[0].datetime_start.strftime('%Y-%m-%d %H:%M:%S'),
                'datetime_complete': study.trials[-1].datetime_complete.strftime('%Y-%m-%d %H:%M:%S')
            }
        }

        # save json
        save_path = os.path.join(self.res_dir, json_path)
        with open(save_path, 'w', encoding='utf-8') as f:
            json.dump(save_data, f, indent=4)
        print(f'Best trial saved to {save_path}')
        return save_data
    
    def execute(self):
        self.plot_history()
        self.plot_param_importance()
        self.save_best_trial_json()
