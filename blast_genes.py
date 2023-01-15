from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
import pandas as pd
import time
#help(NCBIWWW.qblast)

genes_path = 'D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\10quantiles\\non-defying_genes.txt'
seq_file = pd.read_csv(genes_path, sep = ',')
print(seq_file.head())
for i in range(len(seq_file)):
    seq = seq_file.iloc[i].at['translation']
    result_handle = NCBIWWW.qblast('blastp', 'nr', seq, hitlist_size = 10)
    blast_records = NCBIXML.parse(result_handle)
    
    for blast_record in blast_records:
        print(blast_record)
        desc = blast_record.descriptions
        for d in desc:
            print(d("title"))
            d.split('|')
    #res_h = result_handle.read()
    #print(res_h)
    

    print(blast_records)
    time.sleep(5)
    break
