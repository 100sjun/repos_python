from src.data.BMEDDataset import BMEDDataset

def model_with_kfold(mode, datapath):
    df = BMEDDataset(mode, datapath)
    print(df.df.head())
