### Overview
# bolsigdb 파일 생성 모듈

# Module
import pandas as pd

# Convert function
def to_pow(value, decimal_places):
    '''입력된 숫자를 소수점 자리수의 지수로 변환'''
    format_string = f'{{:.{decimal_places}e}}'
    return format_string.format(value)

def col_to_pow(df,column_name,decimal_places):
    '''Dataframe 특정 열의 data를 지수 표기법으로 변환'''
    df[column_name] = df[column_name].apply(lambda x: to_pow(x,decimal_places))
    return df

def run():
    # path
    file_path = './kinetic_DB/'
    i_list = list(range(183))

    # bolsigdb.dat
    with open('bolsigdb.dat','w') as file:
        for i in i_list:
            file_name = f'reaction_{i}.csv'
            f = pd.read_csv(file_path+file_name)
            f = col_to_pow(f,'Energy(eV)',4)
            f = col_to_pow(f,'Cross_section(m2)',4)

            type = f['Type'].iloc[0]
            reaction = f['Reaction'].iloc[0]
            if type == 'ELASTIC' or type == 'EFFECTIVE':
                eloss = str(f['ELoss'].iloc[0]) + ' / mass ratio'
            else:
                eloss = str(f['ELoss'].iloc[0]) + ' / threshold energy'
            line = '-----------------------------'

            file.write(type+'\n')
            file.write(reaction+'\n')
            file.write(eloss+'\n')
            file.write(line+'\n')
            for j in range(len(f)):
                file.write('  ' + str(f['Energy(eV)'].iloc[j]) + '  ' + str(f['Cross_section(m2)'].iloc[j]) + '\n')
            file.write(line + '\n')
            file.write('\n')
        file.close()

if __name__ == '__main__':
    run()