import os
import pandas as pd
from scipy.stats import ttest_ind
import statistics
import seaborn as sns
import numpy as np
from scipy.stats import mannwhitneyu
#path = 'D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\10quantiles\\res\\'
path = 'D:\my_downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\10quantiles\lci_all'
#path = 'D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\10quantiles\lci_rals_test'
lci_list =[]
'''
#для 3 образцов

for org in os.listdir(path):
    df = pd.read_csv(path+'\\'+org, delimiter='\t')
    lci1 = df['LCI1'].to_list()
    lci_list.append(lci1)
    print(org)

    #lci2_list.append(lci2)
for org in os.listdir(path):
    df = pd.read_csv(path+'\\'+org, delimiter='\t')
    lci2 = df['LCI2'].to_list()
    lci_list.append(lci2)
    print(org)
    ax = sns.displot(lci2,
                  bins=100,
                  kde=True,
                  color='skyblue')
    ax.set(xlabel='Normal Distribution', ylabel='Frequency')
    import matplotlib.pyplot as plt
    plt.show()

print(statistics.mean(lci_list[4]), '\n', statistics.mean(lci_list[5]))
'''

'''
#для всех
p_lci1_df = pd.DataFrame(columns=["genes"])
p_lci2_df = pd.DataFrame(columns=["genes"])
n_lci1_df = pd.DataFrame(columns=["genes"])
n_lci2_df = pd.DataFrame(columns=["genes"])
s_lci1_df = pd.DataFrame(columns=["genes"])
s_lci2_df = pd.DataFrame(columns=["genes"])

for org in os.listdir(path):
    df1 = pd.read_csv(path+'\\'+org, delimiter='\t', usecols=['LCI1', 'genes'])
    df1 = df1.drop_duplicates()
    df2 = pd.read_csv(path+'\\'+org, delimiter='\t', usecols=['LCI2', 'genes'])
    df2 = df2.drop_duplicates()
    print(org)
    if 'necator' in org:
        n_lci1_df = pd.merge(df1, n_lci1_df, how="outer", on=["genes"])
        n_lci2_df = pd.merge(df2, n_lci2_df, how="outer", on=["genes"])

    elif 'insidiosa' in org or 'pickettii' in org or 'mannitolilytica' in org:
        s_lci1_df = pd.merge(df1, s_lci1_df, how="outer", on=["genes"])
        s_lci2_df = pd.merge(df2, s_lci2_df, how="outer", on=["genes"])  
    elif 'solanacearum' in org or 'pseudosolanacearum' in org or 'syzygii' in org:
        p_lci1_df = pd.merge(df1, p_lci1_df, how="outer", on=["genes"])
        p_lci2_df = pd.merge(df2, p_lci2_df, how="outer", on=["genes"])
    else:
        print(org)
df_lists = [n_lci1_df, s_lci1_df, p_lci1_df, n_lci2_df, s_lci2_df, p_lci2_df]
for df in df_lists:
    mymean = df.mean(axis = 1)
    lci_list.append(mymean.to_list())
    


for lci in lci_list:
    ax = sns.displot(lci,
                  bins=100,
                  kde=True,
                  color='skyblue')
    ax.set(xlabel='Normal Distribution', ylabel='Frequency')
    import matplotlib.pyplot as plt
    plt.show()
    try:
        print(f'mean = {statistics.mean(lci)}; len = {len(lci)}')
    except:
        print('Empty df')
        

#сохраним картинки  
     
import seaborn as sns
import matplotlib.pyplot as plt
box_res_path = 'D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\10quantiles\\'
lci1_df = pd.DataFrame(data = {'Group':(['Necator']* len(lci_list[0]))+(['Soil']*len(lci_list[1]))+(['Plant']*len(lci_list[2])), 'LCI1':lci_list[0]+lci_list[1]+lci_list[2]})
axn = sns.boxplot(x="Group", y="LCI1", data = lci1_df)
axn.figure.savefig(box_res_path+f'\\lci1.jpg', format='jpeg', dpi=500)
plt.show()

lci2_df = pd.DataFrame(data = {'Group':(['Necator']* len(lci_list[3]))+(['Soil']*len(lci_list[4]))+(['Plant']*len(lci_list[5])), 'LCI2':lci_list[3]+lci_list[4]+lci_list[5]})
axy = sns.boxplot(x="Group", y="LCI2", data = lci2_df)
axy.figure.savefig(box_res_path+f'\\lci2.jpg', format='jpeg', dpi=500)



res_path = 'D:\downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\10quantiles\\lci_stats1.txt'

#посчитаем и сохраним статистику
with open(res_path, 'w') as res:
    res.write('-\tNecator_lci1_t\tNecator_lci1_p\tsoil_lci1_t\tsoil_lci1_p\tplant_lci1_t\tplant_lci1_p\t\tNecator_lci2_t\tNecator_lci2_p\tsoil_lci2_t\tsoil_lci2_p\tplant_lci2_t\tplant_lci2_p\n')
for lci_i in lci_list:
    with open(res_path, 'a') as res:
        res.writelines('-\t')
    for lci_j in lci_list:
        t, p = ttest_ind(lci_i, lci_j, equal_var=False)
        
        with open(res_path, 'a') as res:
            res.writelines(str(t)+'\t'+str(p)+'\t')
    with open(res_path, 'a') as res:
        res.writelines('\n')
'''        
####################################################################
#то же самое, но LCI не усредняем, 1 точка - значение для 1 гена одного генома
#сохраним картинки  

#для всех
p_lci1_df = pd.DataFrame(columns=["genes"])
p_lci2_df = pd.DataFrame(columns=["genes"])
n_lci1_df = pd.DataFrame(columns=["genes"])
n_lci2_df = pd.DataFrame(columns=["genes"])
s_lci1_df = pd.DataFrame(columns=["genes"])
s_lci2_df = pd.DataFrame(columns=["genes"])

for org in os.listdir(path):
    df1 = pd.read_csv(path+'\\'+org, delimiter='\t', usecols=['LCI1', 'genes'])
    df1 = df1.drop_duplicates()
    df2 = pd.read_csv(path+'\\'+org, delimiter='\t', usecols=['LCI2', 'genes'])
    df2 = df2.drop_duplicates()
    print(org)
    if 'necator' in org:
        n_lci1_df = pd.concat([df1, n_lci1_df])
        n_lci2_df = pd.concat([df2, n_lci2_df])

    elif 'insidiosa' in org or 'pickettii' in org or 'mannitolilytica' in org:
        s_lci1_df = pd.concat([df1, s_lci1_df])
        s_lci2_df = pd.concat([df2, s_lci2_df])  
    elif 'solanacearum' in org or 'pseudosolanacearum' in org or 'syzygii' in org:
        p_lci1_df = pd.concat([df1, p_lci1_df])
        p_lci2_df = pd.concat([df2, p_lci2_df])
    else:
        print(org)
df_lists = [n_lci1_df, s_lci1_df, p_lci1_df, n_lci2_df, s_lci2_df, p_lci2_df]
for df in df_lists:
    try:
        lci_list.append(df['LCI1'].to_list())
    except:
        lci_list.append(df['LCI2'].to_list())
    


for lci in lci_list:
    ax = sns.displot(lci,
                  bins=100,
                  kde=True,
                  color='skyblue')
    ax.set(xlabel='Normal Distribution', ylabel='Frequency')
    import matplotlib.pyplot as plt
    plt.show()
    try:
        print(f'mean = {statistics.mean(lci)}; len = {len(lci)}')
    except:
        print('Empty df')
     
import seaborn as sns
import matplotlib.pyplot as plt
box_res_path = 'D:\my_downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\10quantiles\\'
lci1_df = pd.DataFrame(data = {'Group':(['Necator']* len(lci_list[0]))+(['Soil']*len(lci_list[1]))+(['Phyto']*len(lci_list[2])), 'LCI1':lci_list[0]+lci_list[1]+lci_list[2]})
#axn = sns.boxplot(x="Group", y="LCI1", data = lci1_df)
axn = sns.violinplot(x="Group", y="LCI1", data = lci1_df)
axn.figure.savefig(box_res_path+f'\\lci1_nomean_vio.jpg', format='jpeg', dpi=500)

plt.show()

lci2_df = pd.DataFrame(data = {'Group':(['Necator']* len(lci_list[3]))+(['Soil']*len(lci_list[4]))+(['Phyto']*len(lci_list[5])), 'LCI2':lci_list[3]+lci_list[4]+lci_list[5]})
#axy = sns.boxplot(x="Group", y="LCI2", data = lci2_df)
axy = sns.violinplot(x="Group", y="LCI2", data = lci2_df)
axy.figure.savefig(box_res_path+f'\\lci2_nomean_vio.jpg', format='jpeg', dpi=500)


res_path = 'D:\my_downloads\Study\laboratory\EloE\data\downloaded\Ralstonia_complete_genome\\10quantiles\\lci_stats1_no_mean.txt'
#calculate statistics
# Shapiro-Wilk Test
from numpy.random import seed
from numpy.random import randn
from scipy.stats import shapiro
# seed the random number generator
seed(1)
# generate univariate observations
data = 5 * randn(100) + 50

# normality test
def normal(data):
    stat, p = shapiro(data)
    print('Statistics=%.3f, p=%.3f' % (stat, p))
    # interpret
    alpha = 0.05
    if p > alpha:
    	print('Sample looks Gaussian (fail to reject H0)')
    else:
    	print('Sample does not look Gaussian (reject H0)')
normal(data)
with open(res_path, 'w') as res:
    res.write('-\tNecator_lci1_t\tNecator_lci1_p\tsoil_lci1_t\tsoil_lci1_p\tplant_lci1_t\tplant_lci1_p\tNecator_lci2_t\tNecator_lci2_p\tsoil_lci2_t\tsoil_lci2_p\tplant_lci2_t\tplant_lci2_p\n')
for lci_i in lci_list:
    with open(res_path, 'a') as res:
        res.writelines('-\t')
        normal(lci_i)
    for lci_j in lci_list:
        #непараметрический тест манна-уитни
        U1, p=mannwhitneyu(lci_i, lci_j, use_continuity=True, alternative='two-sided')
        #t, p = ttest_ind(lci_i, lci_j, equal_var=False)
        
        with open(res_path, 'a') as res:
            res.writelines(str(U1)+'\t'+str(p)+'\t')
    with open(res_path, 'a') as res:
        res.writelines('\n')


