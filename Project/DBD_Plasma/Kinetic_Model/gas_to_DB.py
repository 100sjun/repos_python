### Overview
# gas kinetic DB 생성 모듈

# Module
import pandas as pd

def run():
    # gas kinetic DB
    path_gas = './crossA_DB/gas.txt'

    with open(path_gas,'r') as f_gas:
        lines = f_gas.readlines()
        data = [line.replace('\n','').rstrip() for line in lines]
        f_gas.close()

    reaction = []
    rate = []
    for i, line in enumerate(data):
        reaction.append(line.split('\t')[0])
        rate.append(line.split('\t')[1])
    reaction = reaction[1:]
    rate = rate[1:]
    df = pd.DataFrame({
        'Num': list(range(len(reaction))),
        'Reaction': reaction,
        'Rate': rate,
    })

    df.to_csv('./kinetic_DB/reaction_gas.csv',index=False)

if __name__ == '__main__':
    run()