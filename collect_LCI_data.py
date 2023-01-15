import os
import re
import pandas as pd

#файл со списком высокоэкспрессируемых генов

stats_df = pd.DataFrame(columns=["locus_tag", "lci1", "lci2", "tai"])
with open('D:\my_downloads\Study\laboratory\kor\\all_genes_efficiency1.csv', 'w') as res_file:
    res_file.write("locus_tag,lci1,lci2,tai\n")
    #папка с результатами елое
    #eloe_path = 'D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\10quantiles\eei_example'
    eloe_path = 'D:\my_downloads\Study\laboratory\kor\eloe_res'
    for folder in os.listdir(eloe_path):
        for myfile in os.listdir(os.path.join(eloe_path, folder)):
            if 'genes_and_flanks.txt' in myfile:
                with open(os.path.join(eloe_path, folder, myfile)) as lci_file:
                    for line in lci_file:
                        if line.startswith('>'):
                            locus = line.split('\t')[3]
                            lci1 =  (line.split('\t')[13]).split('=')[-1]
                            lci2 = (line.split('\t')[14]).split('=')[-1]
                            tai = (line.split('\t')[15]).split('=')[-1]
                            res_file.write(locus+','+lci1+','+lci2+','+tai+"\n")