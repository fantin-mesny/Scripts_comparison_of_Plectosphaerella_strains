import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

data=pd.read_csv('data_GlobalFungi_4SHs.tsv',sep='\t').set_index('Niche').sort_values(by="All samples",ascending=False)
data4plot=data.T
data4plot['sum']=data4plot.sum(axis=1)


data4plot=data4plot.T
fcs=pd.DataFrame(columns=data4plot.columns[:-1],index=data4plot.index[:-2])
for niche in fcs.index:
    for sh in fcs.columns:
	    fcs.loc[niche,sh]=(data4plot.loc[niche,sh]/data4plot.loc[niche,'All samples'])/((data4plot.loc['sum',sh]-data4plot.loc[niche,sh])/(data4plot.loc['sum','All samples']-data4plot.loc[niche,'All samples']))

fig,ax=plt.subplots(1,1,figsize=(20,10))
sns.heatmap(fcs.sort_values(by='SH1002406.10FU').astype(float),center=1,cmap='coolwarm')
plt.savefig('foldchange_niche_4SHs.pdf')
plt.close()
