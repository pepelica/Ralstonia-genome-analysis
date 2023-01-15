# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 10:09:12 2022

@author: Aleksandra
"""
from Bio import SeqIO
import os
import pandas as pd


org_def = []
orgs_new = []
ids = []
types = []
complete = []
#genomeFolder = genomeFolder.replace('_', ' ', genomeFolder.count('_')-1)  #использовать если названия геномов через пробел

#genomesFolder = 'D:\downloads\Study\laboratory\EloE\data\downloaded\genome_assemblies_genome_gb\merged'
#genomesFolder = 'D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\merged'
#genomesFolder = "D:\downloads\Study\laboratory\kor\gbk_files"
genomesFolder = "D:\my_downloads\Study\laboratory\EloE\data\downloaded\genome_assemblies_genome_gb\\ncbi-genomes-2020-08-06"
comp_genomes_path= "D:\my_downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\ncbi-genomes-2021-01-28"
comp_genomes=[genomeFolder.split('_')[0]+'_'+genomeFolder.split('_')[1][:-2] for genomeFolder in os.listdir(comp_genomes_path)]
phylotype = []
for genomeFolder in os.listdir(genomesFolder):
    try:
        contents = os.listdir(genomesFolder+'/'+genomeFolder)           
    except(FileNotFoundError):
        print('file not found')
        #читаем геномную карточку
    ids.append(genomeFolder.split('_')[0]+'_'+genomeFolder.split('_')[1])
    genome = contents[0]
    for record in SeqIO.parse(genomesFolder+'/'+genomeFolder+'/'+genome, 'genbank'):
        pass

    if ids[-1][:-2] in comp_genomes:
        complete.append('complete')
    else:
        complete.append('draft')


    org = record.annotations['organism'].split()[1]
    org_def.append("Ralstonia "+org)
    #использовать если названия геномов через пробел
    #org = record.annotations['organism'].split('_')[1] #для прокки не подходит
    print(org)
    if org == 'insidiosa':
        new_org = 'insidiosa'
        types.append('soil')
    elif org == 'pickettii':
        new_org = 'pickettii'
        types.append('soil')

    elif org == 'mannitolilytica':
        new_org = 'mannitolilytica'
        types.append('soil')

    elif org == 'necator':
        new_org = 'necator'
    elif org == 'syzygii':
        new_org = 'syzygii'
        types.append('phyto')

    elif org == 'solanacearum':
        types.append('phyto')
        print('determine species')
        new_org='sp'
        for record in SeqIO.parse(genomesFolder+'/'+genomeFolder+'/'+genome, 'genbank'):

            #ищем в ней сигнальные последовательности филотпов
            if 'CGTTGATGAGGCGCGCAATTT' in record:
                #phylotype.append('1')
                new_org = 'pseudosolanacearum'
    
            if 'AGTTATGGACGGTGGAAGTC' in record:
                #phylotype.append('2')
                new_org = 'solanacearum'
    
            if 'ATTACGAGAGCAATCGAAAGATT' in record:
                #phylotype.append('3')
                new_org = 'pseudosolanacearum'
    
            if 'ATTGCCAAGACGAGAGAAGTA' in record:
                #phylotype.append('4')
                new_org = 'syzygii'
    else:
        new_org = 'sp'
        types.append(None)
    print(new_org)
    orgs_new.append("Ralstonia "+new_org)
df = pd.DataFrame()
df['Assembly accession'] = ids
df['Assembly accuracy'] = complete
df['Species'] = org_def
df['New species'] = orgs_new
df['type'] = types
with pd.ExcelWriter('D:\my_downloads\Study\laboratory\EloE\data\downloaded\genome_assemblies_genome_gb\\describe2.xlsx') as all_file:
    df.to_excel(all_file)

   
