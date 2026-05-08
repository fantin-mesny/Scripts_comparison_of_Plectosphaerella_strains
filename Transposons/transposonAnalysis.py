import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import numpy as np
import matplotlib as mpl
from scipy import stats

order='P64 P30 P117 P47 P40 P49 P23 P9 P65 P22 P20 P7 P46 P2 P8 P31 P1 P13 P50 P3 P12 P11 P10 P21 P56 P54 P67 P60 P59 P62 P58 P57 P55 P6 P32 P33 P63 P61 P66 P42 P612 P421 P43 P19 P68 P18 P51 P35 P340 P25 P5 P34 P17 P4 P38 P36 P39 P455 P29 P24 P14 P41 P27 P26 P16 P45 P44 P28 P48 P143 P52 P15'.split(' ')
subsp2='P9 P65 P22 P20 P7 P46 P2 P8 P31 P1 P13 P50 P3 P12 P11 P10 P21 P56 P54 P67 P60 P59 P62 P58 P57 P55'.split(' ')
subsp1='P6 P32 P33 P63 P61 P66 P42 P612 P421 P43 P19 P68 P18 P51 P35 P340 P25 P5 P34 P17 P4 P38 P36 P39 P455 P29 P24 P14 P41 P27 P26 P16 P45 P44 P28 P48 P143 P52 P15'.split(' ')


#parsing genome sizes
data=pd.DataFrame(index=order)
for p in order:
    with open('assembly_stats/'+p+'.fasta_stats') as inp:
        inpLines=inp.readlines()
    data.loc[p,'Assembly length']=int(inpLines[1].split(' ')[2][:-1])

#parsing total bp covered by transposons
for p in order:
    with open('transposon_annotation/'+p+'.Statistics_FinalAnnotations.txt') as inp:
        inpLines=inp.readlines()
    for index_line in range(len(inpLines)):
        if '#BP_transposons' in inpLines[index_line] and 'classes' not in inpLines[index_line]:
            #print(p, inpLines[index_line+1].split('\t')[-1][3:-1])
            data.loc[p,'BP_transposons']=int(inpLines[index_line+1].split('\t')[-1][3:-1])
data['Coverage in transposons']=data['BP_transposons']/data['Assembly length']

#parsing transposon counts
allContents=[]
for p in order:
    tc=pd.read_csv('transposon_annotation/'+p+'.Statistics_FinalAnnotations.txt',skiprows=2,nrows=1,sep='\t',header=None).rename(index={0:p})
    tc=tc.drop(columns=[0,1,21])
    tc.columns='Number of transposons  1  1/1  1/1/1  1/1/2  1/1/3  1/2  1/2/1  1/2/2   2  2/1  2/1/1  2/1/2  2/1/3  2/1/4  2/1/5  2/1/6  2/2  2/3'.split('  ')
    allContents.append(tc)
allContents=pd.concat(allContents).astype(float)
allContents=allContents.drop(columns=['1','1/1','1/2',' 2','2/1'])
allContents=allContents.rename(columns={'1/1/1':'Copia','1/1/2':'Gypsy','1/1/3':'ERV','1/2/1':'LINE','1/2/2':'SINE','2/1/1':'Tc1-Mariner','2/1/2':'hAT','2/1/3':'CMC','2/1/4':'Sola','2/1/5':'Zator','2/1/6':'Novosib','2/2':'Helitron','2/3':'MITEs'})
data=data.merge(allContents,left_index=True,right_index=True).reset_index(drop=False)
data=data.rename(columns={'index':'Strain'})

#plot general transposon annotation
fig,ax=plt.subplots(1,3,figsize=(10,16),width_ratios=[1.2,1.2,5],sharey=True)
sns.barplot(y="Strain",x='Coverage in transposons',color='black',ax=ax[0],data=data)
sns.barplot(y="Strain",x='Number of transposons',color='grey',ax=ax[1],data=data)

classes=data.drop(columns=['Assembly length','BP_transposons','Coverage in transposons','Number of transposons'])
classes=pd.DataFrame(classes.set_index('Strain').stack()).reset_index(drop=False).rename(columns={'level_1':'x',0:'number'})
cmap = mpl.cm.Greys(np.linspace(0,1,20))
cmap = mpl.colors.ListedColormap(cmap[2:,:-1])
sns.scatterplot(y='Strain',x='x',hue='number',size='number',palette=cmap,ax=ax[2],data=classes,sizes=(20, 500))

for a in range(3):
    ax[a].xaxis.tick_top()
    ax[a].xaxis.set_label_position('top')
ax[2].tick_params(axis='x', rotation=90)

plt.savefig('transposons_figure.pdf')
plt.close()

#boxplots differences between subspecies
for ind in data.index:
    if data.loc[ind,'Strain'] in subsp1:
        data.loc[ind,'Subspecies']="1"
    elif data.loc[ind,'Strain'] in subsp2:
        data.loc[ind,'Subspecies']="2"
    else:
        data.loc[ind,'Subspecies']="Other"

data=data[data['Subspecies'].isin(["1","2"])]
data['Number of DNA transposons']=data[['Copia','Gypsy','ERV','LINE','SINE']].sum(axis=1)
data['Number of Retrotransposons']=data[['Tc1-Mariner','hAT','CMC','Sola','Zator','Novosib']].sum(axis=1)

S=3
fig,ax=plt.subplots(1,4,figsize=(9,3))
sns.boxplot(data=data,x='Subspecies',y='Coverage in transposons',ax=ax[0],color='lightgray',showfliers=False)
sns.swarmplot(size=S,data=data,x='Subspecies',y='Coverage in transposons',ax=ax[0],color='black')
sns.boxplot(data=data,x='Subspecies',y='Number of transposons',ax=ax[1],color='lightgray',showfliers=False)
sns.swarmplot(size=S,data=data,x='Subspecies',y='Number of transposons',ax=ax[1],color='black')
sns.boxplot(data=data,x='Subspecies',y='Number of DNA transposons',ax=ax[2],color='lightgray',showfliers=False)
sns.swarmplot(size=S,data=data,x='Subspecies',y='Number of DNA transposons',ax=ax[2],color='black')
sns.boxplot(data=data,x='Subspecies',y='Number of Retrotransposons',ax=ax[3],color='lightgray',showfliers=False)
sns.swarmplot(size=S, data=data,x='Subspecies',y='Number of Retrotransposons',ax=ax[3],color='black')
plt.tight_layout()
plt.savefig('diffSubsp1-2.pdf')
plt.close()

## Run statistical testing and print results:
print('Coverage in transposons: ',stats.mannwhitneyu(data[data['Strain'].isin(subsp1)]['Coverage in transposons'],data[data['Strain'].isin(subsp2)]['Coverage in transposons']))
print('Number of transposons:',stats.mannwhitneyu(data[data['Strain'].isin(subsp1)]['Number of transposons'],data[data['Strain'].isin(subsp2)]['Number of transposons']))
print('Number of DNA transposons:',stats.mannwhitneyu(data[data['Strain'].isin(subsp1)]['Number of DNA transposons'],data[data['Strain'].isin(subsp2)]['Number of DNA transposons']))
print('Number of Retrotransposons:',stats.mannwhitneyu(data[data['Strain'].isin(subsp1)]['Number of Retrotransposons'],data[data['Strain'].isin(subsp2)]['Number of Retrotransposons']))
