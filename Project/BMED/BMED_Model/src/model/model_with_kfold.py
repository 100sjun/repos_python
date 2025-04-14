from src.data.BMEDDataset import BMEDDataset
from src.data.BatchSampler import BatchSampler
from src.model.MembraneSystemModel import MembraneSystemModel
from src.model.MembraneSystemTrainer import MembraneSystemTrainer
from torch.utils.data import DataLoader
import numpy as np

def model_with_kfold(mode, datapath, hidden_nodes=64, hidden_layers=3, epochs=100, lr = 0.01, rstop = 0.1, weight_decay = 1e-2):
    fold_results = []
    min_r2_per_fold = []
    for fold in range(5):
        # 현재 fold의 train dataset 생성
        train_dataset = BMEDDataset(
            mode = mode, 
            data_path=datapath, 
            fold_idx=fold, 
            scalers=None, 
            train=True
        )

        # 스케일러 가져오기
        scalers = train_dataset.get_scalers()

        # 현재 fold의 test dataset 생성
        test_dataset = BMEDDataset(
            mode = mode,
            data_path=datapath,
            fold_idx=fold,
            scalers=scalers,
            train=False
        )

        # 데이터 로더 생성
        train_batch_sampler = BatchSampler(train_dataset)
        test_batch_sampler = BatchSampler(test_dataset)
        train_loader = DataLoader(train_dataset, batch_sampler=train_batch_sampler)
        test_loader = DataLoader(test_dataset, batch_sampler=test_batch_sampler)


        # 새로운 모델 초기화
        model = MembraneSystemModel(
            hidden_nodes=hidden_nodes,
            hidden_layers=hidden_layers
        )
        
        trainer = MembraneSystemTrainer(model, epochs=epochs, lr=lr, rstop=rstop, weight_decay=weight_decay)
        
        # 모델 훈련
        train_losses, val_losses, r2_final, r2_min, worst_var, best_model_state = trainer.train(
            train_loader = train_loader, 
            val_loader = test_loader
        )
        
        min_r2_per_fold.append(r2_min)

        # save results
        fold_results.append({
            'fold': fold,
            'r2_final': r2_final,
            'r2_min': r2_min,
            'worst_var': worst_var
        })

        print(f'Fold {fold+1} - min R2: {r2_min:.6f}, Worst Var: {worst_var}')

    # calcualte average min r2
    avg_min_r2 = np.mean(min_r2_per_fold)
    print(f'Average Min R2: {avg_min_r2:.6f}')

    return {
        'fold_results': fold_results,
        'avg_min_r2': avg_min_r2
    }