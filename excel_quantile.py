import os
import pandas as pd
from pandas import ExcelWriter

#path = 'D:\\downloads\\Study\\laboratory\\Results\\upper_quartile1\\sorted'
path = 'D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\4quantiles\sorted'
species = os.listdir(path)
df = pd.DataFrame(columns = ['Org', 'Genes amount', 'Species'])
df_means = pd.DataFrame(columns=['genes'])

for sp in species:
    df_genes = pd.DataFrame(columns=['Org'])
    orgs = os.listdir(path + '\\'+sp)
    for org in orgs:
        for myfile in os.listdir(path + '\\' + sp + '\\' + org):
            if 'quantile1' in myfile:
                genes_dict = {'Org':org}
                quantile = pd.read_csv(path + '\\' + sp + '\\' + org+'\\'+myfile, delimiter = '\t')
                genes = quantile['gene_name'].tolist()
                #print(genes)
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
    
