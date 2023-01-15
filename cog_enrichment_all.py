import pandas as pd
import os
import re
from collections import Counter

path = 'D:\my_downloads\Study\laboratory\kor'
descr = 'cog-20.def.tab.txt'
desc = pd.read_csv(path+'\\'+descr, sep='\t', header=None)
cog_dic = {}
#преобразуем в словарь для ускоренного доступа к данным
for i in range(len(desc)):
    key = desc.iloc[i][0]
    val = [char for char in desc.iloc[i][1]]
    cog_dic[key] = val

def to_group(o):
    org = o.split("_")[1]
    if org == 'insidiosa' or org == 'pickettii' or org == 'mannitolilytica':
        return 'S'

    elif org == 'necator':
        return 'N'

    elif org == 'syzygii' or org == 'solanacearum' or org == 'pseudosolanacearum':
        return 'P'

#посчитаем P S N для каждой из трех групп
plant = {}
soil = {}
necator = {}

for o in os.listdir(path+'\\eloe_res'):
    files = os.listdir(path+'\\eloe_res\\'+o)
    for f in files:
        if re.search('eei[1-5].txt', f):
            print(f)
            
            org = pd.read_csv(path+'\\eloe_res\\'+o+'\\'+f, sep='\t')
            group = to_group(o)
            cogs = []
            for cog in org['COG'].to_list():
                if cog != '-':
                    cogs = cogs +cog_dic[cog]
            org_stat = Counter(cogs)
            if group=="P":
                plant = dict(Counter(plant)+Counter(org_stat))
            elif group=="S":
                soil = dict(Counter(soil)+Counter(org_stat))
            elif group=="N":
                necator = dict(Counter(necator)+Counter(org_stat))
            

dicts = [plant, soil, necator]
res = pd.DataFrame(dicts)
res1 = res.T
res1.columns = ['plant', 'soil', 'necator']
with pd.ExcelWriter(path + '\\PSN_cogs_enrichment_all.xlsx') as myfile:
    res1.to_excel(myfile)



