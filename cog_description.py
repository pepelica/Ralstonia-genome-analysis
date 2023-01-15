import pandas as pd
import os

cog_desc = {}
section = None
letter_desc = {}
 
'''
расшифровка эукариотических KOG

with open("D:\downloads\Study\laboratory\kor\COG one letter code descriptions.txt", "r") as f:
    for line in f:
        if '<b>' in line:
            if '[' not in line :
                section = (line.split('</b>')[0]).split('>')[-1]
                print(section)
            else:
                letter = (line.split('</b>')[0]).split('>')[-1]
                desc = (line.split('</b>')[1]).split('<')[0]
                letter_desc[letter] = (section, desc)

            
with open('D:\downloads\Study\laboratory\kor\\COG function code descriptions.txt') as cog_file:
    for line in cog_file:
        cog = (line.split('<')[0]).split()[1]
        desc = " ".join((line.split('<')[0]).split()[2:])
        letter = (line.split('<')[0]).split()[0]
        cog_desc[cog] = (letter, desc)
'''
with open('D:\my_downloads\Study\laboratory\kor\cogs_letter.txt') as c:
    for line in c:
        letter = line.split('\t')[0]
        descr = line.split('\t')[1]
        letter_desc[letter] = descr

path = 'D:\my_downloads\Study\laboratory\kor\cogs'
files = os.listdir(path)
for myf in files:
    with open(path+'//'+myf) as f:
        for line in f:
            cog = line.split('\t')[2]
            letter = line.split('\t')[3]
            descr = line.split('\t')[4]
            cog_desc[cog] = (letter, descr)

df = pd.read_excel('D:\my_downloads\Study\laboratory\kor\\group_cogs.xlsx')

cog_todf = []
letter_todf = []
for i in range(len(df)):
    cogdf = df.iloc[i]['COG']
    try:
        cog_todf.append((cog_desc[cogdf])[1])
        letter_todf.append((cog_desc[cogdf])[0])
    except:
        letter_todf.append(None)
        cog_todf.append(None)

df['COG description'] = cog_todf
df["COG letter"] = letter_todf

section = []
l_desc = []
for i in range(len(df)):
    lt = df.iloc[i]['COG letter']
    l_d = []
    for l in lt:
        if l != ' ':
            l_d.append(letter_desc[l])
    l_desc.append(l_d)
    '''
    try:
        section.append((letter_desc[lt])[0])
        l_desc.append((letter_desc[lt])[1])
    except:
        section.append(None)
        l_desc.append(None)
    '''
df["COG letter descriprion"] = l_desc
#df["COG section"] = section

with pd.ExcelWriter('D:\my_downloads\Study\laboratory\kor\\group_cogs_descr2.xlsx') as w:
    df.to_excel(w)