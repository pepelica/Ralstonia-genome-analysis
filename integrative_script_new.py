from Bio import SeqIO
from Bio.Seq import Seq, UnknownSeq, Alphabet
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
import pandas as pd
import os
import re


def find_species(genomeFolder):
    #genomeFolder = genomeFolder.replace('_', ' ', genomeFolder.count('_')-1)  #использовать если названия геномов через пробел

    #genomesFolder = 'D:\downloads\Study\laboratory\EloE\data\downloaded\genome_assemblies_genome_gb\merged'
    #genomesFolder = 'D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\merged'
    genomesFolder = "D:\downloads\Study\laboratory\kor\gbk_files"
    phylotype = []

    try:
        contents = os.listdir(genomesFolder+'/'+genomeFolder)           
    except(FileNotFoundError):
        return None 
        #читаем геномную карточку
    
    genome = contents[0]
    record = SeqIO.read(genomesFolder+'/'+genomeFolder+'/'+genome, 'genbank')


    #org = record.annotations['organism'].split()[1]    #использовать если названия геномов через пробел
    #org = record.annotations['organism'].split('_')[1] #для прокки не подходит
    org = genomeFolder.split("_")[1]
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
            #sample = sample.replace('_', ' ', sample.count('_')-1)
            #strain = sample.split(' ')[0]+' '+org+' '+ ' '.join(sample.split(' ')[2:]) #использовать если названия геномов через пробел
            strain = sample.split('_')[0]+'_'+org+'_'+ '_'.join(sample.split('_')[2:])
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
    __slots__ = ('cog', 'gene', 'locus_tag', 'species', 'samples', 'translation', 'definition', 'relative_gene', 'protein_id')



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

    def __init__(self, cog, gene, locus_tag, sample, genomesFolder):
        self.species = {'solanacearum':0, 'syzygii':0, 'pseudosolanacearum':0, 'necator':0, 'mannitolilytica':0, 'insidiosa':0, 'pickettii':0}
        self.samples = []
        self.cog = cog
        self.relative_gene = []
        self.gene = gene
        self.locus_tag = locus_tag
        self.species[sample.split('_')[1]] += 1
        self.samples.append(sample)
        self.get_sequence(genomesFolder)

    def align_sequence(self, seq):
        print('Align')    
        score = 11 
        return score

    def count_stats(self, list_amount):
        sol_amount, syz_amount, pseud_amount, nec_amount, mann_amount, insid_amount, pick_amount = list_amount
        if sol_amount == 0:sol = 0
        else: sol = self.species['solanacearum']/sol_amount 

        if syz_amount == 0: syz = 0
        else: syz = self.species['syzygii']/syz_amount

        if pseud_amount == 0: pseud = 0
        else: pseud = self.species['pseudosolanacearum']/pseud_amount
        
        if nec_amount == 0: nec = 0
        else: nec = self.species['necator']/nec_amount

        if mann_amount == 0: mann = 0
        else: mann = self.species['mannitolilytica']/mann_amount

        if insid_amount == 0: insid = 0
        else: insid = self.species['insidiosa']/insid_amount

        if pick_amount == 0: pick = 0
        else: pick = self.species['pickettii']/insid_amount

        return [sol, syz, pseud, nec, mann, insid, pick]
    
    def count_group(self, amount_list):
        stats_list = self.count_stats(amount_list)
        plant = (stats_list[0]+stats_list[1]+stats_list[2])/3
        soil = (stats_list[4]+stats_list[5]+stats_list[6])/3
        return [plant, soil, stats_list[3]]

    def add_species(self, sample):
        self.samples.append(sample)
        self.species[sample.split('_')[1]] += 1

#path_eei = 'D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\EloE_res'
path_eei = "D:\downloads\Study\laboratory\kor\eloe_res"
#path_percentile = 'D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\10quantiles\quantile'
path_percentile = "D:\downloads\Study\laboratory\kor\quantile"
genomesFolder =  "D:\downloads\Study\laboratory\kor\gbk_files"
#genomesFolder = 'D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\merged'
quantile_amount = 10
species_statistics = {'solanacearum':0, 'syzygii':0, 'pseudosolanacearum':0, 'necator':0, 'mannitolilytica':0, 'insidiosa':0, 'pickettii':0}
groups_statistics = {'plant':0, 'necator':0, 'soil':0}

#path = 'D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\4quantiles\sorted'
samples = os.listdir(path_eei)
df = pd.DataFrame(columns = ['Org', 'Genes amount', 'Species'])
df_means = pd.DataFrame(columns=['genes'])
geneObjects = {}
geneName = 1
#myfile = open('D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\non-defying_genes.txt', 'w')
#myfile.close()
'''
for s in samples:
    df_genes = pd.DataFrame(columns=['Org'])
    org = os.listdir(path_eei + '\\'+s)
    species = find_quantile(path_eei, path_percentile, s, quantile_amount)
'''

for qfile in os.listdir(path_percentile):        
    # + '\\' + f'{org}_quantile{i+1}.txt):
    species = qfile[:-18]
    print(species)
    species_name = species.split('_')[1]
    species_statistics[species_name] += 1
    if 'quantile1' in qfile:
        quantile = pd.read_csv(path_percentile + '\\'+qfile, delimiter = '\t')
        genes = quantile['gene_name'].tolist()
        cogs = quantile["COG"].tolist()

        locuses = quantile['locus_tag'].tolist()


        #создадим объект для каждого гена организма
        for i in range(len(cogs)):
            added = False
            #проверяем, записаны ли гены с таким названием ( если нет, создаем новый объект)
            if cogs[i] != '-':
                for obj in geneObjects.keys():
                    if cogs[i] == obj:
                        geneObjects[obj].add_species(species)
                        added = True
                    #создаем новый объект
                if added == False:
                    geneObjects[cogs[i]] = Gene(cogs[i], genes[i], locuses[i], species, genomesFolder)
            '''
            #если у гена нет названия
            else:
                for obj in geneObjects.keys():
                    seq = geneObjects[obj].translation

                    score = geneObjects[obj].align_sequence(seq)
                    if score > 10:
                        geneObjects[obj].add_species(species)
                        added = True
                                    
                if added == False:
                    geneObjects['gene'+str(geneName)] = Gene('gene'+str(geneName), locuses[i], species, genomesFolder)
                    with open('D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\non-defying_genes.txt', 'a') as myfile:
                        myfile.write('gene'+str(geneName)+'\t'+geneObjects['gene'+str(geneName)].samples[0]+'\t'+geneObjects['gene'+str(geneName)].protein_id[0]+'\t'+locuses[i]+'\t'+geneObjects['gene'+str(geneName)].translation[0])
                    geneName += 1
            '''

df_res = pd.DataFrame()
df_int_res = pd.DataFrame()
for key in geneObjects.keys():
    df_res[key] = geneObjects[key].count_stats(list(species_statistics.values()))
    df_int_res[key] = geneObjects[key].count_group(list(species_statistics.values()))
df_res = df_res.transpose()
#df_res.columns = ["COG", "solanacearum", "syzygii", "pseudosolanacearum", "necator", "mannitolilytica", "insidiosa", "pickettii"]
print(df_res)

df_int_res = df_int_res.transpose()
#df_int_res.columns = ["COG", "plant", "soil", "necator"]
with pd.ExcelWriter(f'D:\downloads\Study\laboratory\kor\\species_cogs.xlsx') as sp_file:
    df_res.to_excel(sp_file)

with pd.ExcelWriter(f'D:\downloads\Study\laboratory\kor\\group_cogs.xlsx') as sp1_file:
    df_int_res.to_excel(sp1_file)



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

