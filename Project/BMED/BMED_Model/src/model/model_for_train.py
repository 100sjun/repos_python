from src.data.BMEDDataset import BMEDDataset
from src.data.BatchSampler import BatchSampler
from torch.utils.data import DataLoader
from src.model.MembraneSystemModel import MembraneSystemModel
from src.model.MembraneSystemTrainer import MembraneSystemTrainer

def model_for_train(mode, datapath, fold_idx=0, hidden_nodes=64, hidden_layers=3, lr=0.01, rstop=0.1, weight_decay=1e-2, epochs=100):
    # train dataset 생성
    train_dataset = BMEDDataset(
        mode = mode,
        data_path=datapath,
        fold_idx=fold_idx,
        scalers=None,
        train=True
    )

    # 스케일러 가져오기
    scalers = train_dataset.get_scalers()

    # 현재 fold의 test dataset 생성
    test_dataset = BMEDDataset(
        mode = mode,
        data_path=datapath,
        fold_idx=fold_idx,
        scalers=scalers,
        train=False
    )

    # 데이터 로더 생성
    train_batch_sampler = BatchSampler(train_dataset)
    test_batch_sampler = BatchSampler(test_dataset)
    train_loader = DataLoader(train_dataset, batch_sampler=train_batch_sampler)
    test_loader = DataLoader(test_dataset, batch_sampler=test_batch_sampler)

    # 모델 생성
    model = MembraneSystemModel(
        hidden_nodes=hidden_nodes,
        hidden_layers=hidden_layers
    )

    # 트레이너 생성
    trainer = MembraneSystemTrainer(
        model=model,
        lr=lr,
        rstop=rstop,
        weight_decay=weight_decay,
        epochs=epochs
    )

    # 학습 시작
    train_losses, val_losses, r2_final, r2_min, worst_var, best_model_state = trainer.train(train_loader, test_loader,mode)

    return r2_final,r2_min, worst_var, best_model_state, scalers

