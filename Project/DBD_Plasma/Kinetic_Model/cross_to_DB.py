### Overview
# cross section are txt 파일을 kinetic_DB 폴더에 CSV 파일로 저장하는 모듈

# Modules
import pandas as pd
import re
import io

def cross_DB(num):
    # cross area data path
    index_ca = num    # cross area data index
    path = f'./crossA_DB/{index_ca}.txt'

    # Reading a text file
    with open(path, 'r') as file:
        lines = file.readlines()
        data = [line.replace('\n','').rstrip() for line in lines]
        file.close()

    # Meta Data
    type = data[0]
    reaction = data[1]
    eloss = float(data[2])
    species = re.search(r'SPECIES: (.+)', data[3]).group(1)
    formula = re.search(r'PROCESS: (.+?),',data[4]).group(1)

    reaction = re.sub(r'(\w+)\+', r'\1^+', reaction)
    species = re.sub(r'(\w+)\+', r'\1^+', species)
    formula = re.sub(r'(\w+)\+', r'\1^+', formula)

    # cross area data
    data_start = data.index('-----------------------------') + 1
    data_end = data[data_start:].index('-----------------------------') + data_start
    data_lines = data[data_start:data_end]

    # Dataframe
    data_io = io.StringIO('\n'.join(data_lines))
    df = pd.read_csv(data_io, sep=r'\s+', names=['Energy(eV)', 'Cross_section(m2)'])

    # Dataframe Metadata Insertion
    df['Num'] = index_ca
    df['Type'] = type
    df['Reaction'] = reaction
    df['ELoss'] = eloss
    df['Species'] = species
    df['Formula'] = formula

    df = df[['Num', 'Type', 'Reaction', 'ELoss', 'Species', 'Formula', 'Energy(eV)', 'Cross_section(m2)']]

    return df

def run():
    for num in range(183):
        df = cross_DB(num)
        df.to_csv(f'./kinetic_DB/reaction_{num}.csv', index=False)

if __name__ == '__main__':
    run()