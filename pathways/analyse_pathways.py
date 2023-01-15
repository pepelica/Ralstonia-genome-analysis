import pandas as pd
import requests
import json
from xml.etree import ElementTree as ET
import time

#распарсить файл с путями
path = 'D:\my_downloads\Study\laboratory\MinPath'
pw_file = pd.read_csv(path+'\\venn_result.txt', delimiter= '\t')
pathways = {}

report = open(path+'//ecocyc_report1_2.txt', 'w')

name = None
pw_file = pw_file.fillna(0)
for i in range(len(pw_file)):
    n = pw_file['Names'][i]
    if n != 0:
        name = n
    pathways[pw_file['elements'][i]] = [name]

#распарсить файл rest api, получить список генов
def pw_genes(pw):
    #api_link = 'https://websvc.biocyc.org/apixml?fn=genes-of-pathway&id=PWY0-1479'
    #api_link = 'https://websvc.biocyc.org/apixml?fn=genes-of-pathway&id='+str(pw)
    api_link = 'https://websvc.biocyc.org/apixml?fn=genes-of-pathway&id=META:'+str(pw)
    gene_list = []
    headers = {'Accept': 'application/json;odata=verbose'}
    responce=requests.get(api_link,verify=True)
    root = ET.fromstring(responce.content) 
    print(root.tag)
    print(root.attrib)
    for child in root:

        for product in child.findall('common-name'):
            pass
            try:
                g_name = product.text
                print(g_name)
                gene_list.append(g_name)
                
                #try:
                #    for s in product.findall('synonym'):
                #        g_name.append(s.text)
                #except:
                #    pass


            except:
                pass
    return gene_list

for ptway in pathways:
    try:
        genes = pw_genes(ptway) 
        report.write(ptway+' получен\n')  
        time.sleep(15)
    except:
        genes = []
        print(ptway+' не содержится в biocyc\n')
        report.write(ptway+' не содержится в biocyc\n')
    pathways[ptway].append(genes)

report.close()

#сопоставить гены с таблицей, присвоить им групповую принадлежность
table_p = 'D:\my_downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\10quantiles\gene_name_psn_descr.xlsx'
table = pd.read_excel(table_p, sheet_name='Sheet1')
pws_list= []
pws_type = []

for i in range(len(table)):
    g_pw = []
    t_pw = []
    gene = table['genes'][i]
    for pw in pathways:
        if gene in pathways[pw][1]:
            g_pw.append(pw)
            t_pw.append(pathways[pw][0]+'|')
    pws_list.append(g_pw)
    pws_type.append(t_pw)
table['pathways'] = pws_list
table['pathways type'] = pws_type
with pd.ExcelWriter(path+'\\genes_pathways_type1_2.xlsx') as res:
    table.to_excel(res)
