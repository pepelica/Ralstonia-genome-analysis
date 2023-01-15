import pandas as pd
import requests
import json
import time
import os
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

#распарсить файл с путями
path = 'D:\downloads\Study\laboratory\MinPath'
data_path = 'D:\downloads\Study\laboratory\MinPath\mp_details'
files = os.listdir(data_path)

pw_file = pd.read_csv(path+'\\venn_result.txt', delimiter= '\t')
pathways = {}

report = open(path+'//enzyme_report.txt', 'w')

names = {}
name = None
pw_file = pw_file.fillna(0)
for i in range(len(pw_file)):
    n = pw_file['Names'][i]
    if n != 0:
        name = n
    pathways[pw_file['elements'][i]] = [name, [], []]
    try:
        names[name].append(pw_file['elements'][i])
    except:
        names[name] = [pw_file['elements'][i]]

enzymes = {}
det = 'D:\downloads\Study\laboratory\MinPath\minpath_reports'
#распарсить details, получить список enzymes
for det_file in os.listdir(det):
    if 'details' in det_file:
        with open(det+'//'+det_file) as details:
            pw = details.readline()
            for line in details:
                if 'path' in line:
                    pw = line.split(' ')[1]
                else:
                    enz = line.split(' ')[-1][:-1]
                    try:
                        pathways[pw][1].append(enz)
                        enzymes[enz]= []
                    except:pass
print(pathways)

#удалим дупликации
for key in pathways:
    pathways[key][1] = list( dict.fromkeys(pathways[key][1]) )

gbk_p = 'D:\downloads\Study\laboratory\kor\gbk_files'
#распарсить gbk, получить список генов
for gbk_file in os.listdir(gbk_p):
    genebank = os.listdir(gbk_p+'\\'+gbk_file)[0]  
    record = SeqIO.read(gbk_p+'\\'+gbk_file+'\\'+genebank, 'genbank')
    for feature in record.features:
        if feature.type == "CDS":
            try:
                enzyme = feature.qualifiers['EC_number'][0]
                gene = feature.qualifiers['gene'][0]
                try:
                    enzymes[enzyme].append(gene)
                except:
                    pass
            except:
                pass
                
#удалим дупликации
for key in pathways:
    pathways[key][2] = list( dict.fromkeys(pathways[key][2]) )

#сопоставим путям гены
for ptway in pathways:
    genes = []
    enzs = pathways[ptway][1]
    for enz in enzs:
        genes = genes + enzymes[enz]
    pathways[ptway][2] = genes


#сопоставить гены с таблицей, присвоить им групповую принадлежность
table_p = 'D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\10quantiles\gene_name_psn_descr.xlsx'
table = pd.read_excel(table_p, sheet_name='Sheet3')
pws_list= []
pws_type = []

for i in range(len(table)):
    g_pw = []
    t_pw = []
    gene = table['genes'][i]
    for pw in pathways:
        if gene in pathways[pw][2]:
            g_pw.append(pw)
            t_pw.append(pathways[pw][0]+'|')
    pws_list.append(g_pw)
    pws_type.append(t_pw)
table['pathways'] = pws_list
table['pathways type'] = pws_type
with pd.ExcelWriter(path+'\\genes_pathways_type_enz_all.xlsx') as res:
    table.to_excel(res)
