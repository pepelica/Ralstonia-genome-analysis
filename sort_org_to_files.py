import os
import shutil
from Bio import SeqIO

#path = 'D:\\downloads\\Study\\laboratory\\Results\\upper_quartile1'
path = 'D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\4quantiles'

def find_phylotype(genomeFolder):
    #genomesFolder = 'D:\downloads\Study\laboratory\EloE\data\downloaded\genome_assemblies_genome_gb\merged'
    genomesFolder = 'D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\merged'
    phylotype = []

    try:
        contents = os.listdir(genomesFolder+'/'+genomeFolder)           
    except(FileNotFoundError):
        return None 
        #читаем геномную карточку
    
    genome = contents[0]
    record = SeqIO.read(genomesFolder+'/'+genomeFolder+'/'+genome, 'genbank')

    org = record.annotations['organism'].split()[1]
    print(org)
    if org == 'insidiosa':
        return 'insidiosa'
    elif org == 'pickettii':
        return 'pickettii'

    elif org == 'mannitolilytica':
        return 'mannitolilytica'

    elif org == 'necator':
        return 'necator'

    else:

        #ищем в ней сигнальные последовательности филотпов
        if 'CGTTGATGAGGCGCGCAATTT' in record:
            #phylotype.append('1')
            return 1

        if 'AGTTATGGACGGTGGAAGTC' in record:
            #phylotype.append('2')
            return 2

        if 'ATTACGAGAGCAATCGAAAGATT' in record:
            #phylotype.append('3')
            return 3

        if 'ATTGCCAAGACGAGAGAAGTA' in record:
            #phylotype.append('4')
            return 4
        else:
            return None


for folder in os.listdir(path+'\\quantile'):
    
    if len(os.listdir(path + '\\quantile\\'+folder)) == 0:
        os.rmdir(path+'\\quantile\\'+folder)
    else:
        phyl = find_phylotype(folder)
        try:
            if phyl == 4:
                shutil.move(path+'\\quantile\\'+folder, path+'\\sorted\\syzygii')

            elif phyl == 3 or phyl == 1:
                shutil.move(path+'\\quantile\\'+folder, path+'\\sorted\\pseudosolanacearum')

            elif phyl == 2:
                shutil.move(path+'\\quantile\\'+folder, path+'\\sorted\\solanacearum')

            elif phyl == 'mannitolilytica':
                shutil.move(path+'\\quantile\\'+folder, path+'\\sorted\\mannitolilytica')

                
            elif phyl == 'pickettii':
                shutil.move(path+'\\quantile\\'+folder, path+'\\sorted\\pickettii')
            
            elif phyl == 'insidiosa':
                shutil.move(path+'\\quantile\\'+folder, path+'\\sorted\\insidiosa')
        except:
            continue



