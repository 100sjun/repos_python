### Overview
# kinet.inp 파일 생성

# Module
import pandas as pd
import re

# Species 정렬함수
def custom_sort_key(s):
    # 기호가 가장 앞에 오도록 우선순위 부여
    if not re.match(r'[A-Za-z]', s):
        return(0, s)
    # C가 적은 순으로 우선순위 부여
    elif s in ['C','CH','CH^+']:
        return(1, s)
    elif re.match(r'^CH\d+',s):
        return(2, s)
    elif re.match(r'^C2H\d+',s):
        return(3, s)
    elif re.match(r'^C3H\d+',s):
        return(4, s)
    elif re.match(r'^C4H\d+',s):
        return(5, s)
    return (6, s)

def run():
    # Path
    path = './kinetic_DB/'
    f = pd.read_csv(path+'reaction_10.csv')

    # data treatment
    species = []
    bolsig = []
    formula_e = []
    bol_comp = []
    formula_g = []
    rate_g = []
    formula_i = []
    rate_i = []

    for i in range(183):
        file_name = f'reaction_{i}.csv'
        f1 = pd.read_csv(path+file_name)
        species += f1['Formula'].iloc[0].split(' ')
        bolsig += f1['Species'].iloc[0].split(' ')
        # formula 
        if f1['Type'].iloc[0] == 'EXCITATION':
            form = f1['Formula'].iloc[0].replace('E','e').replace('->','=>')
            b_comp = form.replace('e','').replace('+','').replace('=>','->')
            b_comp = b_comp.split(' ')
            b_comp = [item for item in b_comp if item != '']
            b_comp = ' '.join(b_comp)
            formula_e.append(form)
            bol_comp.append(b_comp)
            
        elif f1['Type'].iloc[0] == 'IONIZATION':
            form = f1['Formula'].iloc[0].replace('E','e').replace('->','=>')
            b_comp = form.replace('e','').replace('+','').replace('=>','->') + '+'
            b_comp = b_comp.split(' ')
            b_comp = [item for item in b_comp if item != '']
            b_comp = ' '.join(b_comp)
            formula_e.append(form)
            bol_comp.append(b_comp)
        else:
            pass

    f2 = pd.read_csv(path+'reaction_gas.csv')
    for i in range(len(f2)):
        species += f2['Reaction'].iloc[i].split(' ')
        formula_g.append(f2['Reaction'].iloc[i])
        rate_g.append(f2['Rate'].iloc[i])

    f3 = pd.read_csv(path+'reaction_ion.csv')
    for i in range(len(f3)):
        formula_i.append(f3['Reaction'].iloc[i])
        rate_i.append(f3['Rate'].iloc[i])

    species = sorted(list(set(species)),key=custom_sort_key)
    species.remove('+')
    species.remove('->')
    species.remove('=>')
    species.remove('E')
    species.remove('e')
    species = ['e'] + species

    if len(species) > 25:
        species1 = species[:25]
        species2 = species[25:]
    else:
        species1 = species
        species2 = []

    species1 = ' '.join(species1)
    species2 = ' '.join(species2)

    bolsig = list(set(bolsig))
    bolsig = sorted(list(set(bolsig)),key=custom_sort_key)
    bolsig.remove('/')
    bolsig.remove('e')
    bolsig = ' '.join(bolsig)

    # overlapping between ground and vibrational excitation
    t_rxn = formula_e + formula_i + formula_g
    df_rxn = pd.DataFrame({'Reaction': t_rxn})
    count_index = 0
    index = []
    for i in range(len(df_rxn)):
        if i > 0:
            curr_react, curr_prod = df_rxn['Reaction'].iloc[i].split('=>')
            prev_react, prev_prod = df_rxn['Reaction'].iloc[i-1].split('=>')
            
            curr_react_norm = re.sub(r'\(v\d*\)','',curr_react).strip()
            prev_react_norm = re.sub(r'\(v\d*\)','',prev_react).strip()

            if curr_react_norm == prev_react_norm and curr_prod == prev_prod:
                pass
            else:
                count_index += 1
        index.append(count_index)
    df_rxn['index'] = index

    df_rxn.to_csv('parameter_set.csv',index=False)

    # compensation parameter
    para = []
    for i in range(df_rxn['index'].iloc[-1]+1):
        para_sentence = f'$ double precision, parameter :: f{i} = 1.000d-21'
        para.append(para_sentence)

    kinet = ''

    kinet += 'ELEMENTS\n'
    kinet += 'e C H\n'
    kinet += 'END\n'
    kinet += '\n'
    kinet += 'SPECIES\n'
    kinet += species1+'\n'
    kinet += species2+'\n'
    kinet += 'END\n'
    kinet += '\n'
    kinet += 'BOLSIG\n'
    kinet += bolsig + '\n'
    kinet += 'END\n'
    kinet += '\n'
    kinet += 'REACTIONS\n'
    kinet += '# Gas Constant\n'
    kinet += '$ double precision, parameter :: R = 8.314d-3\n'
    kinet += '\n'
    kinet += '# Compensation Parameters\n'
    for i in range(len(para)):
        kinet += para[i] + '\n'
    kinet += '\n'
    kinet += '# Electron Collision Reaction\n'
    for i in range(len(formula_e)):
        kinet += formula_e[i] + '\t'*2 + f'! f{df_rxn['index'].iloc[i]} * Bolsig+ ' + bol_comp[i] + '\n'
    kinet += '\n'
    kinet += '# Ion Recombination Reaction\n'
    for i in range(len(formula_i)):
        kinet += formula_i[i] + '\t'*2 + f'! f{df_rxn['index'].iloc[i+len(formula_e)]} * ' + rate_i[i] + '\n'
    kinet += '\n'
    kinet += '# Gas Phase Reaction\n'

    for i in range(len(formula_g)):
        kinet += formula_g[i] + '\t'*2 + f'! f{df_rxn['index'].iloc[i+len(formula_e)+len(formula_i)]} * ' + rate_g[i] + '\n'
    kinet += 'END'


    with open('kinet.inp', 'w') as file:
        file.write(kinet)
        file.close()

if __name__ == '__main__':
    run()

