import pandas as pd
import os
path = 'D:\downloads\Study\laboratory\MinPath\pathways'
necator = {}
plant = {}
soil = {}
n = 0
p = 0
s = 0

def cont_d(pw_file, org_dict):
    for pw in pw_file:
        try:
            org_dict[pw] += 1
        except:
            org_dict[pw] = 1
    return org_dict

for org in os.listdir(path):
    print(org) 
    pws = []
    with open(path+'\\'+org) as pw_file:
        for line in pw_file:
            pws.append(line[:-1])

    if 'necator' in org:
        necator = cont_d(pws, necator)
        n += 1
    elif 'insidiosa' in org or 'pickettii' in org or 'mannitolilytica' in org:
        soil = cont_d(pws, soil)
        s += 1 
    elif 'solanacearum' in org or 'pseudosolanacearum' in org or 'syzygii' in org:
        plant = cont_d(pws, plant)
        p += 1
    else:
        print(org)

res = 'D:\downloads\Study\laboratory\MinPath'

with open(res+'\\necator.txt', 'w') as nec:
    for key in necator:
        if (necator[key]/n) > 0.8:
            nec.write(key+'\n')

with open(res+'\\soil.txt', 'w') as f:
    for key in soil:
        if (soil[key]/s) > 0.8:
            f.write(key+'\n')

with open(res+'\\plant.txt', 'w') as f:
    for key in plant:
        if (plant[key]/p) > 0.8:
            f.write(key+'\n')
