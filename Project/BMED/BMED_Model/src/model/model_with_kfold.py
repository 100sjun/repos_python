from src.data.BMEDDataset import BMEDDataset
from src.data.BatchSampler import BatchSampler
from torch.utils.data import DataLoader

def model_with_kfold(mode, datapath):
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

        print("train_loader:", train_loader)
        print("test_loader:", test_loader)

