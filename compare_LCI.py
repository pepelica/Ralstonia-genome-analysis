import os
import re
import pandas as pd

#файл со списком высокоэкспрессируемых генов
high_expressed_path = 'D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\10quantiles\high_expressed.txt'
high_ex = pd.read_csv(high_expressed_path, delimiter='\t')
print(len(high_ex))
eei_list = len(high_ex)*[None]
id_list = len(high_ex)*[None]
lci1 = len(high_ex)*[None]
lci2 = len(high_ex)*[None]

#папка с результатами елое
#eloe_path = 'D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\10quantiles\eei_example'
eloe_path = 'D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\EloE_res'
for folder in os.listdir(eloe_path):
    for myfile in os.listdir(os.path.join(eloe_path, folder)):
        #файл с ИЭЭ
        #получим идентификатор для каждого гена
        if re.search('_eei\d.txt', myfile):
            id_df = pd.read_csv(os.path.join(eloe_path, folder, myfile), delimiter='\t')
            for i in range(len(id_df)):
                for j in range(len(high_ex)):
                    if high_ex['genes'][j] == id_df['gene_name'][i]:
                        eei_list[j] = id_df['EEI'][i]
                        id_list[j] = id_df['locus_tag'][i]

            high_ex['eei'] = eei_list
            high_ex['locus_tag'] = id_list   
            print(id_df)
            pass

        #файл с LCI

        if 'genes_and_flanks.txt' in myfile:
            with open(os.path.join(eloe_path, folder, myfile)) as lci_file:
                for line in lci_file:
                    if line.startswith('>'):
                        locus = line.split('\t')[3]
                        for i in range(len(id_list)):
                            if id_list[i] == locus:
                                lci1[i] =  (line.split('\t')[13]).split('=')[-1]

                                lci2[i] = (line.split('\t')[14]).split('=')[-1]
            high_ex['LCI1'] = lci1
            high_ex['LCI2'] = lci2


    high_ex.to_csv(f'D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\10quantiles\\lci_all\\{folder}_lci1.txt', sep = '\t')

        
#каждому гену добавим соответствующие индексы lci1, lci2
#
#