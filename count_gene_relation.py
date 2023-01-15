import pandas as pd
from pandas import ExcelWriter


#path = 'D:\\downloads\\Study\\laboratory\\EloE\data\downloaded\Ralstonia_complete_genome\\4quantiles\\quantile1_def.xlsx'
#path = 'D:\downloads\Study\laboratory\kor\\group_cogs_descr1.xlsx'
path = 'D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\10quantiles\group_cogs.xlsx'
df = pd.read_excel(path)

df = df.drop_duplicates()

nec = df['necator'].to_list()
soil = df['soil'].to_list()
plant = df['plant'].to_list()
'''

path = 'D:\\downloads\\Study\\laboratory\\EloE\data\downloaded\Ralstonia_complete_genome\\4quantiles\\quantile1_def.xlsx'
df = pd.read_excel(path, sheet_name = '1')
df = df.drop_duplicates()

nec = df['pseudosolanacearum'].to_list()
soil = df['solanacearum'].to_list()
plant = df['syzygii'].to_list()

types = []
for i in range(len(df)):
    if nec[i] < 0.2:
        if soil[i] > 0.8:
            if plant[i] > 0.8:
                #types.append('sp')
                types.append('so_sy')
            elif plant[i] < 0.2:
                #types.append('s')
                types.append('so')
            else:
                types.append(None)

        elif soil[i] < 0.2:
            if plant[i] > 0.8:
                #types.append('p')
                types.append('sy')
            elif plant[i] < 0.02:
                types.append('-')
            else:
                types.append(None)
        else:
            types.append(None)

    elif nec[i] > 0.8:
        if soil[i] > 0.8:
            if plant[i] < 0.2:
                #types.append('ns')
                types.append('ps_so')
            elif plant[i] > 0.8:
                #types.append('nps')
                types.append('ps_so_sy')
            else:
                types.append(None)

            
        elif soil[i] < 0.2:
            if plant[i] > 0.8:
                #types.append('np')
                types.append('ps_sy')
            elif plant[i] < 0.2:
                #types.append('n')
                types.append('ps')
            else:
                types.append(None)
        else:
            types.append(None)

    else:
        types.append(None)

'''
types = []
for i in range(len(df)):
    if nec[i] < 0.2:
        if soil[i] > 0.8:
            if plant[i] > 0.8:
                types.append('sp')

            elif plant[i] < 0.2:
                types.append('s')

            else:
                types.append(None)

        elif soil[i] < 0.2:
            if plant[i] > 0.8:
                types.append('p')

            elif plant[i] < 0.02:
                types.append('-')
            else:
                types.append(None)
        else:
            types.append(None)

    elif nec[i] > 0.8:
        if soil[i] > 0.8:
            if plant[i] < 0.2:
                types.append('ns')

            elif plant[i] > 0.8:
                types.append('nps')

            else:
                types.append(None)

            
        elif soil[i] < 0.2:
            if plant[i] > 0.8:
                types.append('np')

            elif plant[i] < 0.2:
                types.append('n')
            else:
                types.append(None)
        else:
            types.append(None)

    else:
        types.append(None)

df['types'] = types

#with ExcelWriter('D:\\downloads\\Study\\laboratory\\EloE\data\downloaded\Ralstonia_complete_genome\\4quantiles\\quantile1_pl_types.xlsx') as all_file:
#with ExcelWriter('D:\downloads\Study\laboratory\kor\\group_cogs_psn.xlsx') as all_file:
with ExcelWriter('D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\10quantiles\gene_name_psn.xlsx') as all_file:
    df.to_excel(all_file)