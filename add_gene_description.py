import pandas as pd
from Bio import SeqIO
from pandas import ExcelWriter

#path = 'D:\\downloads\\Study\\laboratory\\EloE\data\downloaded\Ralstonia_complete_genome\\4quantiles\\quantile1.xlsx'
path = 'D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\10quantiles\gene_name_psn.xlsx'
#df = pd.read_excel(path, sheet_name ='Means')
df = pd.read_excel(path)
df_genes = pd.DataFrame({'genes':[None], 'function':[None], 'inference':[None], 'product':[None]})
genes = df['genes'].to_list()

genomeFolder = "D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\merged\\ralstonia pseudosolanacearum   ['NZ_CP021764']\\ralstonia pseudosolanacearum.gbk"
record = SeqIO.read(genomeFolder, 'genbank')
for feature in record.features:
    if feature.type == 'CDS':
        try:
            gene_is = feature.qualifiers.get("gene")[0]


            if gene_is in genes:
                function_is = feature.qualifiers.get("function")
                inference_is = feature.qualifiers.get("inference")
                product_is = feature.qualifiers.get("product")
                df_genes = pd.concat([df_genes, pd.DataFrame({'genes':[gene_is], 'function':[function_is], 'inference':[inference_is], 'product':[product_is]})])
        except:
            pass
#org = record.annotations['organism'].split()[1]
df = pd.merge(df, df_genes, how='left', on = 'genes')
#with ExcelWriter('D:\\downloads\\Study\\laboratory\\EloE\data\downloaded\Ralstonia_complete_genome\\4quantiles\\quantile1_def.xlsx') as all_file:
with ExcelWriter('D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\10quantiles\gene_name_psn_descr.xlsx') as all_file:
    df.to_excel(all_file)
