from Bio import SeqIO
from Bio.Seq import Seq, UnknownSeq, Alphabet
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
import pandas as pd
import os
import re
from pandas import ExcelWriter
from Bio import pairwise2


def find_species(genomeFolder):
    genomeFolder = genomeFolder.replace('_', ' ', genomeFolder.count('_')-1)

    #genomesFolder = 'D:\downloads\Study\laboratory\EloE\data\downloaded\genome_assemblies_genome_gb\merged'
    genomesFolder = 'D:\my_downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\merged'
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
    elif org == 'syzygii':
        return 'syzygii'

    else:

        #ищем в ней сигнальные последовательности филотпов
        if 'CGTTGATGAGGCGCGCAATTT' in record:
            #phylotype.append('1')
            return 'pseudosolanacearum'

        if 'AGTTATGGACGGTGGAAGTC' in record:
            #phylotype.append('2')
            return 'solanacearum'

        if 'ATTACGAGAGCAATCGAAAGATT' in record:
            #phylotype.append('3')
            return 'pseudosolanacearum'

        if 'ATTGCCAAGACGAGAGAAGTA' in record:
            #phylotype.append('4')
            return 'syzygii'
        else:
            return None      

#находит файл eei и разбивает его на заданное количество частей
def find_quantile(EloE_output_path, output, sample, quantile):
    EloEresults = os.listdir(EloE_output_path+'\\'+sample)
    for EloE_file in EloEresults:
        if re.search('_eei\d.txt', EloE_file):
            indx_type = EloE_file[-5]
            df = pd.read_csv(os.path.join(EloE_output_path, sample, EloE_file), delimiter='\t', 
                                   usecols = ['locus_tag', 'protein_id', 'gene_id',
                                              'gene_name', 'COG', 'EEI'])
            df['EEI_type'] = [indx_type]*len(df)
            org = find_species(sample)
            if org == None:
                return None
            sample = sample.replace('_', ' ', sample.count('_')-1)
            strain = sample.split(' ')[0]+' '+org+' '+ ' '.join(sample.split(' ')[2:])
            start = 0
            end = 0
            for i in range(quantile):
                start = end
                end = len(df)-round(len(df)/quantile*(quantile - 1-i))
#                df1 = df.iloc[round(i*len(df)/4):len(df)-round(len(df)/4*(quantile - 1-i)), :]
                df1 = df.iloc[start:end, :]
                file_path = os.path.join(output, f'{strain}_quantile{i+1}.txt')
                df1.to_csv(file_path, sep = '\t', na_rep = '-')
                break
    return strain

class Gene():
    __slots__ = ('gene', 'locus_tag', 'species', 'samples', 'translation', 'definition', 'relative_gene', 'protein_id')



    def get_sequence(self, genomesFolder):        
        try:
            contents = os.listdir(genomesFolder+'/'+self.samples[0])

        except(FileNotFoundError):
            return None 
            #читаем геномную карточку
        genome = contents[0]
        record = SeqIO.read(genomesFolder+'/'+self.samples[0]+'/'+genome, 'genbank')
        for feature in record.features:
            if feature.type == "CDS":
                locus = feature.qualifiers['locus_tag']
                if locus[0] == self.locus_tag:
                    self.translation = feature.qualifiers['translation']
                    try:
                        self.definition = feature.qualifiers['product']
                        if 'Traceback' in self.definition:
                            raise Exception
                    except:
                        self.definition = None
                    try:
                        self.protein_id = feature.qualifiers['protein_id']
                        if 'Traceback' in self.protein_id:
                            raise Exception
                    except:
                        self.protein_id = None

    def __init__(self, gene, locus_tag, sample, translation, definition, protein_id):
        self.species = {'solanacearum':0, 'syzygii':0, 'pseudosolanacearum':0, 'necator':0, 'mannitolilytica':0, 'insidiosa':0, 'pickettii':0}
        self.samples = []
        self.relative_gene = []
        self.gene = gene
        self.locus_tag = locus_tag
        self.species[sample.split()[1]] += 1
        self.samples.append(sample)
        #self.get_sequence(genomesFolder)
        self.translation = translation
        self.definition = definition
        self.protein_id = protein_id

    def align_sequence(self, seq):
  
        alignments = pairwise2.align.globalxs(seq, self.translation, -9, -1)
        score = alignments[0][2]/(len(seq)+len(self.translation))*2
        return score


    def count_stats(self, args): #args = sol_amount, syz_amount, pseud_amount, nec_amount, mann_amount, insid_amount
        return [self.species['solanacearum']/args[0], self.species['syzygii']/args[1], self.species['pseudosolanacearum']/args[2], self.species['necator']/args[3], self.species['mannitolilytica']/args[4], self.species['insidiosa']/args[5]]
    
    def add_species(self, sample):
        self.samples.append(sample)
        self.species[sample.split()[1]] += 1

path_eei = 'D:\my_downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\EloE_res'
path_percentile = 'D:\my_downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\10quantiles\quantile'
genomesFolder = 'D:\my_downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\merged'
quantile_amount = 10
species_statistics = {'solanacearum':1, 'syzygii':1, 'pseudosolanacearum':1, 'necator':1, 'mannitolilytica':1, 'insidiosa':1, 'pickettii':1}




def define_tax_group(s):
    if 'solanacearum' in s or 'syzygii' in s:
        return 'P'
    elif  'mannitolilytica' in s or 'insidiosa' in s or 'pickettii' in s:
        return 'S'
    elif 'necator' in s:
        return 'N'
#path = 'D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\4quantiles\sorted'

#файлы елое гены верней квантили(PHEG) с COG
path_pheg='D:\my_downloads\Study\laboratory\kor\quantile_10'
#файлы елое с COG (все гены) - названия файлов отличаются
path_all='D:\my_downloads\Study\laboratory\kor\eloe_res'
path_all='D:\my_downloads\Study\laboratory\kor\eloe_res'

# читаем все образцы
df=pd.DataFrame({'tax':['P_PHEG', 'S_PHEG', 'N_PHEG', 'P_all', 'S_all', 'N_all']})
samples = os.listdir(path_pheg)
for s in samples:
    tax=define_tax_group(s)
    
    # подсчтываем коги верхней квантили
    file_df=pd.read_csv(path_pheg+'\\'+s, sep='\t')
    cogs=file_df['COG'].to_list()
    for cog in cogs:
        if cog in df.columns:
            #df.iloc['COG']+=1
            df.loc[df['tax'] == tax+'_PHEG', cog] +=1
        else:
            df[cog]=[0]*len(df)
            df.loc[df['tax'] == tax+'_PHEG', cog] +=1
print(df)


# добавим иноформацию по всем остальным когам
samples = os.listdir(path_all)
for s in samples:
    for EloE_file in os.listdir(path_all+'\\'+s):
        if re.search('_eei\d.txt', EloE_file):
            tax=define_tax_group(s)

            file_df=pd.read_csv(path_all+'\\'+s+'\\'+EloE_file, sep='\t')
            cogs=file_df['COG'].to_list()
            for cog in cogs:
                if cog in df.columns:
                    #df.iloc['COG']+=1
                    df.loc[df['tax'] == tax+'_all', cog] +=1
                else:
                    df[cog]=[0]*len(df)
                    df.loc[df['tax'] == tax+'_all', cog] +=1
print(df)


#список всех когов - столбцы таблицы
df_res = df.T

with ExcelWriter('D:\my_downloads\Study\laboratory\kor\cog_fisher\\all_cogs.xlsx') as all_file:
    df_res.to_excel(all_file)


