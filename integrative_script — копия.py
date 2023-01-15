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

path_eei = 'D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\EloE_res'
path_percentile = 'D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\10quantiles\quantile'
genomesFolder = 'D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\merged'
quantile_amount = 10
species_statistics = {'solanacearum':1, 'syzygii':1, 'pseudosolanacearum':1, 'necator':1, 'mannitolilytica':1, 'insidiosa':1, 'pickettii':1}

#path = 'D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\4quantiles\sorted'
samples = os.listdir(path_eei)
df = pd.DataFrame(columns = ['Org', 'Genes amount', 'Species'])
df_means = pd.DataFrame(columns=['genes'])
geneObjects = {}
geneName = 1
myfile = open('D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\10quantiles\\non-defying_genes.txt', 'w')
myfile.write('gene_name,strain,locus_tag,protein_id,translation\n')
myfile.close()
'''
for s in samples:
    df_genes = pd.DataFrame(columns=['Org'])
    org = os.listdir(path_eei + '\\'+s)
    species = find_quantile(path_eei, path_percentile, s, quantile_amount)
'''
genes_from_record = []
for qfile in os.listdir(path_percentile):        
    # + '\\' + f'{org}_quantile{i+1}.txt):
    species = qfile[:-14]
    print(species)
    species_name = species.split()[1]
    species_statistics[species_name] += 1
    if 'quantile1' in qfile:
        quantile = pd.read_csv(path_percentile + '\\'+qfile, delimiter = '\t')
        genes = quantile['gene_name'].tolist()

        locuses = quantile['locus_tag'].tolist()
        contents = os.listdir(genomesFolder+'/'+species)
            #читаем геномную карточку
        genome = contents[0]
        record = SeqIO.read(genomesFolder+'/'+species+'/'+genome, 'genbank')
        for feature in record.features:
            if feature.type == "CDS":
                locus = feature.qualifiers['locus_tag']
                
                for i in range(len(locuses)):
                    if locus[0] == locuses[i]:
                        
                        translation = feature.qualifiers['translation']
                        try:
                            definition = feature.qualifiers['product']
                            if 'Traceback' in definition:
                                raise Exception
                        except:
                            definition = None
                        try:
                            protein_id = feature.qualifiers['protein_id']
                            if 'Traceback' in protein_id:
                                raise Exception
                        except:
                            protein_id = None
                        new_gene = Gene(genes[i], locuses[i], species, translation[0], definition[0], protein_id[0])
                        genes_from_record.append(new_gene)
                        locuses.pop(i)
                        genes.pop(i)
                        break
            if len(genes) == 0:
                break

        del record

        
        #genes_dict = {'Org':species}
        

        #print(genes)
        
        #создадим объект для каждого гена организма
        for i in range(len(genes_from_record)):
            added = False
            #проверяем, записаны ли гены с таким названием ( если нет, создаем новый объект)
            if genes_from_record[i].gene != '-':
                for obj in geneObjects.keys():

                    #одинаковые названия гена
                    if genes_from_record[i].gene == geneObjects[obj].gene:
                        geneObjects[obj].add_species(species) #добавляем название штамма к гену
                        added = True
                    
                    #одинаковые описания гена
                    #плохой критерий - есть разные гены с одинаковым описанием (напр. 'hypothetical protein')
                    elif genes_from_record[i].definition == geneObjects[obj].definition:
                        if  genes_from_record[i].definition == "hypothetical protein":
                            geneObjects[obj].add_species(species) #добавляем название штамма к гену
                            added = True
                
                    
            
                    #копируем объект к списку генов
                if added == False:
                    geneObjects[genes_from_record[i].gene] = genes_from_record[i]
            #если у гена нет названия
            else:
                
                for obj in geneObjects.keys():
                    seq = genes_from_record[i].translation
                    score = geneObjects[obj].align_sequence(seq)
                    print(score)
                    if score > 0.7:
                        geneObjects[obj].add_species(species)
                        added = True
                
                                    
                if added == False:
                    geneObjects['gene'+str(geneName)] = genes_from_record[i]
                    with open('D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\10quantiles\\non-defying_genes.txt', 'a') as myfile:
                        myfile.write('gene'+str(geneName)+','+geneObjects['gene'+str(geneName)].samples[0]+','+geneObjects['gene'+str(geneName)].protein_id+','+geneObjects['gene'+str(geneName)].locus_tag+','+geneObjects['gene'+str(geneName)].translation+'\n')
                    geneName += 1
            
    break
df_res = pd.DataFrame()
for key in geneObjects.keys():
    #df_res[key] = geneObjects[key].count_stats(list(species_statistics.values()))+ [geneObjects[key].definition]
    df_res[geneObjects[key].gene] = geneObjects[key].count_stats(list(species_statistics.values()))+ [geneObjects[key].definition]

df_res = df_res.transpose()
df_res.columns = ['solanacearum', 'syzygii', 'pseudosolanacearum', 'necator', 'mannitolilytica', 'insidiosa', 'pickettii']
print(df_res)
with ExcelWriter(f'D:\\downloads\\Study\\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\10quantiles\\genes.xlsx') as sp_file:
    df_res.to_excel(sp_file)
'''
            for col in df_genes.columns[1:]:
            #   print(col)
                if col in genes:
                    genes_dict[col] = 1
                else:
                    genes_dict[col] = 0

            df_genes = df_genes.append(genes_dict, ignore_index=True)

            for gene in genes:                   
                if gene not in df_genes.columns:
                    df_genes[gene] = [0] *(len(df_genes)-1)+ [1]

        
            df = df.append({'Org': org, 'Genes amount':len(quantile), 'Species':sp}, ignore_index=True)
my_mean = df_genes.mean(axis = 0)
df_mean = my_mean.to_frame(name = sp)
df_mean['genes'] = df_genes.columns[1:]
print(df_mean, type(my_mean))
#df_mean = pd.DataFrame(data = [sp] + my_mean, columns=df_genes.columns)
df_means = df_means.merge(df_mean, how = 'outer', on='genes')

#    with ExcelWriter(f'D:\\downloads\\Study\\laboratory\\Results\\upper_quartile1\\{sp}.xlsx') as sp_file:
with ExcelWriter(f'D:\\downloads\\Study\\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\4quantiles\\{sp}.xlsx') as sp_file:
    df_genes.to_excel(sp_file)
                

df_means = df_means.fillna(0)
print(df)
with ExcelWriter('D:\\downloads\\Study\\laboratory\\EloE\data\downloaded\Ralstonia_complete_genome\\4quantiles\\quantile1.xlsx') as all_file:
    df.to_excel(all_file, sheet_name='All orgs')
    df_means.to_excel(all_file, sheet_name='Means')
    
'''

