from src.data.BMEDDataset import BMEDDataset

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
        print(scalers)
