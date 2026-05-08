import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


## Parsing:
order='P64 P30 P117 P47 P40 P49 P23 P9 P65 P22 P20 P7 P46 P2 P8 P31 P1 P13 P50 P3 P12 P11 P10 P21 P56 P54 P67 P60 P59 P62 P58 P57 P55 P6 P32 P33 P63 P61 P66 P42 P612 P421 P43 P19 P68 P18 P51 P35 P340 P25 P5 P34 P17 P4 P38 P36 P39 P455 P29 P24 P14 P41 P27 P26 P16 P45 P44 P28 P48 P143 P52 P15'.split(' ')#[::-1]
order_mapper={p:order.index(p) for p in order}
df=pd.read_csv('data_temperature.csv').set_index('Unnamed: 0')
df['order']=df['strain'].map(order_mapper)
df=df.sort_values(by='order')
df=df[df['strain'].isin(order)]
metadata=pd.read_csv('metadata_phyloAndHost.csv').rename(columns={'Fungus':'strain'})
df=df.merge(metadata,on='strain')

## Plotting
fig,ax=plt.subplots(1,3,figsize=(10,10),sharey=True)
sns.boxplot(y='strain',x='25 degrees', data=df,ax=ax[0],color='lightgray')
sns.swarmplot(y='strain',x='25 degrees', data=df,ax=ax[0],size=2,color='black')
sns.boxplot(y='strain',x='30 degrees', data=df,ax=ax[1],color='lightgray')
sns.swarmplot(y='strain',x='30 degrees', data=df,ax=ax[1],size=2,color='black')
sns.barplot(y='strain',x='Growth difference',data=df,ax=ax[2],color='gray')
plt.savefig('experimentResults.pdf')
plt.close()

## Calculate mean values per strain
df_phylo=df[['strain','Growth difference']].groupby('strain').mean()
df_phylo=df_phylo.merge(metadata,on='strain')
df_host=df[['strain','Growth difference']].groupby('strain').mean()
df_host=df_host.merge(metadata,on='strain')

## plot difference per phylogenetic group, host phylogenetic group and continent of origin
fig,ax=plt.subplots(3,1,figsize=(5,12))
sns.boxplot(y='PhyloGroup',x='Growth difference', data=df_phylo,ax=ax[0],color='lightgray',fliersize=0,order=['subsp1','subsp2','notInSubsp','P. delsorboi','P. plurivora'])
sns.swarmplot(y='PhyloGroup',x='Growth difference', data=df_phylo,ax=ax[0],color='black',size=2.5,order=['subsp1','subsp2','notInSubsp','P. delsorboi','P. plurivora'])
sns.boxplot(y='HostGroup',x='Growth difference', data=df_host,ax=ax[1],color='lightgray',order=['Brassicaceae','Cucurbitaceae','Fabaceae','Rosaceae','Solanaceae','Other dicots','Monocots','Marchantiaceae','Not plant-isolated','Unknown'],fliersize=0)
sns.swarmplot(y='HostGroup',x='Growth difference', data=df_host,ax=ax[1],color='black',size=2.5,order=['Brassicaceae','Cucurbitaceae','Fabaceae','Rosaceae','Solanaceae','Other dicots','Monocots','Marchantiaceae','Not plant-isolated','Unknown'])
sns.boxplot(y='Continent of isolation',x='Growth difference', data=df_host,ax=ax[2],color='lightgray',order=['Europe','North America','Asia','Oceania','Africa','Unknown'],fliersize=0)
sns.swarmplot(y='Continent of isolation',x='Growth difference', data=df_host,ax=ax[2],color='black',size=2.5,order=['Europe','North America','Asia','Oceania','Africa','Unknown'])
plt.savefig('growthDifferencePerGroup.pdf')
plt.close()
