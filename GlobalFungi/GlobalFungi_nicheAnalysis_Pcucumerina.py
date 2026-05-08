import matplotlib.pyplot as plt
from statsmodels.graphics.mosaicplot import mosaic
import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
from statsmodels.stats import multitest


## Parsing data
data=pd.read_csv('data_GlobalFungi_Pcucumerina.csv',sep=';').set_index('Niche').sort_values(by="All samples",ascending=False)[['All samples', 'Pc']]
data['without Pc']=data['All samples']-data['Pc']
data4plot=data.T
data4plot['sum']=data4plot.sum(axis=1)


## Compute test for overrepresentation of each niche in samples containing P. cucumerina
ps=[]
cat=[]
signif=[]
for c in data4plot.columns[:-1]:
    #prepare table for Fisher's exact test (one sided: tests enrichment but not depletion)
    table=[[data4plot.loc['Pc',c], data4plot.loc['without Pc',c]],[data4plot.loc['Pc','sum']-data4plot.loc['Pc',c], data4plot.loc['without Pc','sum']-data4plot.loc['without Pc',c]]]
    res=stats.fisher_exact(table, alternative='greater')
    ps.append(res.pvalue)
    cat.append(c)                 
                        
FDRs=multitest.multipletests(ps, alpha=0.05, method='fdr_bh') # Adjusted p-values
testResults=pd.concat([pd.DataFrame(cat).rename(columns={0:'Niche'}),pd.DataFrame(ps).rename(columns={0:'P-value'}),pd.DataFrame(FDRs[1]).rename(columns={0:'FDR'})],axis=1)
print(testResults.sort_values(by='FDR'))

## Calculate percentages and fold-change values for plotting
data4plot=data4plot.T
data4plot['Percentage of samples with P. cucumerina']=(data4plot['Pc']/data4plot['All samples'])*100
for niche in data4plot.index:
	data4plot.loc[niche,'FoldChange']=(data4plot.loc[niche,'Pc']/data4plot.loc[niche,'All samples'])/((data4plot.loc['sum','Pc']-data4plot.loc[niche,'Pc'])/(data4plot.loc['sum','All samples']-data4plot.loc[niche,'All samples']))

data4plot=data4plot.sort_values(by='FoldChange',ascending=False)
data4plot=data4plot.reset_index(drop=False)
data4plot=data4plot[:-1]


## Plotting:
fig,ax=plt.subplots(1,2,width_ratios=(1,0.7),sharey=True,figsize=(6,5))
sns.barplot(ax=ax[0],data=data4plot,x='Percentage of samples with P. cucumerina',y='Niche',color='black')
ax[0].set_xlim(0,30)
sns.barplot(ax=ax[1],data=data4plot,x='FoldChange',y='Niche',color='grey')
ax[1].axvline(x=1, ymin=0, ymax=12,color='grey')
ax[1].set_xlim(0,5)
N=0
for ind in data4plot.index:
    ax[0].text(data4plot.loc[ind,'Percentage of samples with P. cucumerina']+1,N,data4plot.loc[ind,'Pc'],va='center')
    if data4plot.loc[ind,'Niche'] in signif:
         ax[1].text(data4plot.loc[ind,'FoldChange']+0.25,N,'*',va='center')
    N+=1


plt.savefig('foldchange_niche.pdf')
plt.close()
