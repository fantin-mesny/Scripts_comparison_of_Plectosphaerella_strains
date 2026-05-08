import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import math
from numpy import std, mean, sqrt

def hedges_g(x,y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    SDpooled=sqrt(((nx-1)*std(x, ddof=1) ** 2 + (ny-1)*std(y, ddof=1) ** 2) / dof)
    return (mean(x) - mean(y)) / SDpooled


## Parse data
df=pd.read_csv('data_shootFreshWeight.csv')
df_metadata=df[['HostFamily','HostGroup','order','Treatment']].set_index('Treatment').drop_duplicates()
hostFam=df_metadata['HostFamily'].to_dict()
hostGroup=df_metadata['HostGroup'].to_dict()
order=df_metadata['order'].to_dict()

palette={'Brassicaceae':'#6699FF',
	'Not plant-isolated':'#BFBFBF',
	'Other dicots':'#79d279',
	'Solanaceae':'#FF6666',
	'Cucurbitaceae':'#FFBF80',
	'Monocots':'#FFFF66',
	'Rosaceae':'#D9B3FF',
	'Fabaceae':'#D9B38C',
	'Marchantiaceae':'#999966',
	'Unknown':'#FFFFFF'}
phylo_palette={'subsp1':'#0000CD',
               'subsp2':'#1E90FF',
            	'notInSubsp':'#8E00F6',
                'P. delsorboi':'#DF6767',
                'P. plurivora':'#38771D'}

hue_order=['Brassicaceae','Cucurbitaceae','Solanaceae','Fabaceae','Rosaceae','Other dicots','Monocots','Marchantiaceae','Not plant-isolated','Unknown']
phyloDic={'P64':'P. delsorboi','P30':'P. plurivora','P117':'P. plurivora','P47':'notInSubsp','P40':'notInSubsp','P49':'notInSubsp','P23':'notInSubsp','P6':'subsp1', 'P32':'subsp1', 'P33':'subsp1', 'P63':'subsp1', 'P61':'subsp1', 'P66':'subsp1', 'P42':'subsp1', 'P612':'subsp1', 'P421':'subsp1', 'P43':'subsp1', 'P19':'subsp1', 'P68':'subsp1', 'P18':'subsp1', 'P51':'subsp1', 'P35':'subsp1', 'P340':'subsp1', 'P25':'subsp1', 'P5':'subsp1', 'P34':'subsp1', 'P17':'subsp1', 'P4':'subsp1', 'P38':'subsp1', 'P36':'subsp1', 'P39':'subsp1', 'P455':'subsp1', 'P29':'subsp1', 'P24':'subsp1', 'P14':'subsp1', 'P41':'subsp1', 'P27':'subsp1', 'P26':'subsp1', 'P16':'subsp1', 'P45':'subsp1', 'P44':'subsp1', 'P28':'subsp1', 'P48':'subsp1', 'P143':'subsp1', 'P52':'subsp1', 'P15':'subsp1'}
for p in set(df.Treatment):
	if p not in phyloDic:
		phyloDic[p]='subsp2'

host_order=['Brassicaceae','Cucurbitaceae','Fabaceae','Rosaceae','Solanaceae','Other dicots','Monocots','Marchantiaceae', 'Not plant-isolated','Unknown']
phylo_order=['subsp1','subsp2','notInSubsp','P. delsorboi','P. plurivora']



## Calculate Hedges' g
df_h=[]
for f in set(df.Treatment):
    x=df[df['Treatment']==f]['SFW']
    y=df[df['Treatment']=='MOCK']['SFW']
    df_h.append({'Fungus':f,"Hedges' g":hedges_g(x,y),'mean SFW':x.mean()})
df_h=pd.DataFrame(df_h).sort_values(by="Hedges' g").set_index('Fungus')
df_h=df_h.dropna()
df_h['HostGroup']=df_h.index.map(hostGroup)
df_h['PhyloGroup']=df_h.index.map(phyloDic)
df_h[df_h.index!='MOCK'].to_csv('data_forStatisticalAnalysis.csv')


## Plot
fig,ax=plt.subplots(2,1,figsize=(5.2,5),sharex=True,height_ratios=(2,3.5))
sns.boxplot(ax=ax[1],y='HostGroup',x="Hedges' g",data=df_h,order=host_order,palette=palette,linewidth=0.7,fliersize=0)
sns.swarmplot(ax=ax[1],y='HostGroup',x="Hedges' g",order=host_order,data=df_h,size=3.25,legend=None,linewidth=0.1,color='black',edgecolor='white')
sns.boxplot(ax=ax[0],y='PhyloGroup',x="Hedges' g",data=df_h,order=phylo_order,palette=phylo_palette,linewidth=0.7,fliersize=0)
sns.swarmplot(ax=ax[0],y='PhyloGroup',x="Hedges' g",order=phylo_order,data=df_h,size=3.25,legend=None,linewidth=0.1,color='black',edgecolor='white')
ax[0].set_xlim(-2.25,0.6)
ax[1].set_xlim(-2.25,0.6)
plt.savefig('hedges.pdf')
plt.close()
