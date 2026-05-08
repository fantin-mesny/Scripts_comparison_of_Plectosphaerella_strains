import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from statsmodels.stats import multitest

## Parse orthology prediction results
og=pd.read_csv('data_Orthogroups.tsv',sep='\t').set_index('Orthogroup').fillna('')
og_unassigned=pd.read_csv('data_Orthogroups_UnassignedGenes.tsv',sep='\t').set_index('Orthogroup')
og=pd.concat([og,og_unassigned.fillna('')])
og_c=pd.read_csv('data_Orthogroups.GeneCount.tsv',sep='\t').set_index('Orthogroup')
og_unassigned_c=og_unassigned.where(pd.isna(og_unassigned), 1).fillna(0)
og_c=pd.concat([og_c,og_unassigned_c])
og_pa=og_c.drop(columns='Total')
og_pa[og_pa>1]=1
og_c=og_c.drop(columns='Total')

## Define phylogenetic groups
otherSp=['P64','P30','P117'] # P. plurivora and P. delsorboi
notInSubsp=['P47','P40','P49','P23'] # P. cucumerina outgroup
subsp1='P6 P32 P33 P63 P61 P66 P42 P612 P421 P43 P19 P68 P18 P51 P35 P340 P25 P5 P34 P17 P4 P38 P36 P39 P455 P29 P24 P14 P41 P27 P26 P16 P45 P44 P28 P48 P143 P52 P15'.split(' ')
subsp2=[p for p in og.columns if p not in otherSp+notInSubsp+subsp1]

og_pa=og_pa.drop(columns=otherSp) #Remove non-P. cucumerina strains from dataset


## Calculate percentages genomes with gene families
og_pa['pcOfGenomes']=100*(og_pa.sum(axis=1)/len(notInSubsp+subsp1+subsp2))
og_pa['pcOfGenomesInSubsp1']=100*(og_pa[subsp1].sum(axis=1)/len(subsp1))
og_pa['pcOfGenomesInSubsp2']=100*(og_pa[subsp2].sum(axis=1)/len(subsp2))
og_pa=og_pa[og_pa['pcOfGenomes']>0]

## Classify each gene family into categories based on their prevalence among phylogenetic groups:
coreOG=list(og_pa[og_pa['pcOfGenomes']==100].index) # present in all
broadlyConservedOG=list(og_pa[(og_pa['pcOfGenomes']>=75) & (og_pa['pcOfGenomes']<100)].index) #present in 75-99% genomes
subsp1OG=list(og_pa[(og_pa['pcOfGenomesInSubsp1']>=20) & (og_pa['pcOfGenomesInSubsp2']==0)].index) #Subsp1-specific
subsp2OG=list(og_pa[(og_pa['pcOfGenomesInSubsp2']>=20) & (og_pa['pcOfGenomesInSubsp1']==0)].index) #Subsp2-specific
lowConservationOG=list(og_pa[(og_pa['pcOfGenomes']<=20) & ~(og_pa.index.isin(subsp1+subsp2))].index) # present in less than 20% genomes
classifiedOG=coreOG+broadlyConservedOG+lowConservationOG+subsp1OG+subsp2OG
sparselyDistributedOG=[o for o in og_pa.index if o not in classifiedOG] #not classified in other groups
totalOG=sparselyDistributedOG+classifiedOG

def OGtoConservation(og):
	if og in coreOG:
		return 'core'
	elif og in broadlyConservedOG:
		return 'broadlyConserved'
	elif og in lowConservationOG:
		return 'lowConservation'
	elif og in sparselyDistributedOG:
		return 'sparseDistribution'
	elif og in subsp1OG:
		return 'subsp1'
	elif og in subsp2OG:
		return 'subsp2'
	else:
		return 'NOT FOUND'

def gen_colors(df):
    col_d = {'core':'#454444','broadlyConserved':'#807e7e','sparseDistribution':'#D3D3D3','subsp1':'#0000CD','subsp2':'#1E90FF','lowConservation':'#F5F5F5'}
    return [col_d[col] for col in df.columns]


## Number of genes at each conservation level for each Plectosphaerella cucumerina isolate
data=pd.DataFrame(index='P47 P40 P49 P23 P9 P65 P22 P20 P7 P46 P2 P8 P31 P1 P13 P50 P3 P12 P11 P10 P21 P56 P54 P67 P60 P59 P62 P58 P57 P55 P6 P32 P33 P63 P61 P66 P42 P612 P421 P43 P19 P68 P18 P51 P35 P340 P25 P5 P34 P17 P4 P38 P36 P39 P455 P29 P24 P14 P41 P27 P26 P16 P45 P44 P28 P48 P143 P52 P15'.split(' ')[::-1],columns=['core','broadlyConserved','sparseDistribution','subsp1','subsp2','lowConservation']).fillna(0)
for o in og_c.index:
	for p in data.index:
		cons=OGtoConservation(o)
		if cons!='NOT FOUND':
			data.loc[p,cons]+=og_c.loc[o,p]
			
## Plot content of each genome
data.plot(kind='barh', stacked=True,figsize=(4,15),color=gen_colors(data))
plt.savefig('barplot_geneConservation.pdf')
plt.close()


## Plot per-genome proportion of gene conservation levels
data['sum']=data.sum(axis=1)
for ind in data.index:
    for col in data.columns[:-1]:
        data.loc[ind,col]=data.loc[ind,col]/data.loc[ind,'sum']
data=data.drop(columns='sum').stack().reset_index(drop=False).rename(columns={'level_1':'Gene conservation',0:'Per-genome percentage of genes'})
plt.figure()
sns.boxplot(y='Gene conservation',x='Per-genome percentage of genes',data=data,fliersize=3)
plt.savefig('boxplots_geneConservation.pdf')
plt.close()



## Load functional annotations
annotCSEP=pd.read_csv('data_annotationEffectors.csv')
annotCAZ=pd.read_csv('data_annotationCazymes.csv')
annotProteases=pd.read_csv('data_annotationProteases.csv')
annotTotal=pd.DataFrame(index=range(0,len(totalOG)),columns=['Orthogroup'])
annotTotal['Orthogroup']=totalOG

## Link function to conservation level
def addConservation(annot):
    for ind in annot.index:
        annot.loc[ind,'Conservation']=OGtoConservation(annot.loc[ind,'Orthogroup'])
    annot=annot[['Orthogroup','Conservation']].drop_duplicates()
    return annot
annotCSEP=addConservation(annotCSEP)
annotCAZ=addConservation(annotCAZ)
annotProteases=addConservation(annotProteases)
annotTotal=addConservation(annotTotal)

## Get functional groups (CAzymes, Proteases, Effectors) of orthogroups per conservation levels
data=[annotTotal['Conservation'].value_counts().to_dict(),annotCAZ['Conservation'].value_counts().to_dict(),annotProteases['Conservation'].value_counts().to_dict(), annotCSEP['Conservation'].value_counts().to_dict()]
data[3]['genes']='CSEPs'
data[2]['genes']='Proteases'
data[1]['genes']='CAZymes'
data[0]['genes']='Total'

## Plot percentage of CAzymes, Proteases, Effectors in each conservation level
data=pd.DataFrame(data).set_index('genes')
data=data[['core','broadlyConserved','sparseDistribution','subsp1','subsp2','lowConservation']]
data['sum']=data.sum(axis=1)
for ind in data.index:
    for col in data.columns[:-1]:
        data.loc[ind,col]=data.loc[ind,col]/data.loc[ind,'sum']
data=data.drop(columns='sum')[['core','broadlyConserved','sparseDistribution','subsp1','subsp2','lowConservation']]
data.plot(kind='bar', stacked=True,figsize=(1.5,4),color=gen_colors(data))
plt.savefig('barplot.pdf')
plt.close()


## Compute Fisher's test for overrepresentation of CAzymes, Proteases, Effectors in each conservation level
def performTestForGroup(OGsOfInterest,OGsTotal,annotation):
	OGsNotOfInterest=[og for og in OGsTotal if og not in OGsOfInterest]

	annotatedOfInterest=[og for og in set(annotation['Orthogroup']) if og in OGsOfInterest]
	annotatedNotOfInterest=[og for og in set(annotation['Orthogroup']) if og in OGsNotOfInterest]
	notAnnotatedOfInterest=[og for og in OGsOfInterest if og not in set(annotation['Orthogroup'])]
	notAnnotatedNotOfInterest=[og for og in OGsNotOfInterest if og not in set(annotation['Orthogroup'])]

	statistic,pvalue=stats.fisher_exact([[len(annotatedOfInterest),len(annotatedNotOfInterest)],[len(notAnnotatedOfInterest),len(notAnnotatedNotOfInterest)]],alternative='two-sided')
	return statistic,pvalue,len(annotatedOfInterest)/(len(annotatedOfInterest)+len(notAnnotatedOfInterest)),len(annotatedNotOfInterest)/(len(annotatedNotOfInterest)+len(notAnnotatedNotOfInterest))

groups={'core':coreOG,'broadlyConserved':broadlyConservedOG,'sparselyDistributed':sparselyDistributedOG,'subsp1':subsp1OG,'subsp2':subsp2OG,'lowConservation':lowConservationOG}
output=[]
for group in groups:
	cazres=performTestForGroup(groups[group],totalOG,annotCAZ)
	output.append({'Function':'CAZymes','Conservation':group,'statistic':cazres[0],'pvalue':cazres[1],'percentage annotated in group':cazres[2],'percentage annotated in rest':cazres[3]})
	csepres=performTestForGroup(groups[group],totalOG,annotCSEP)
	output.append({'Function':'Effectors','Conservation':group,'statistic':csepres[0],'pvalue':csepres[1],'percentage annotated in group':csepres[2],'percentage annotated in rest':csepres[3]})
	protres=performTestForGroup(groups[group],totalOG,annotProteases)
	output.append({'Function':'Proteases','Conservation':group,'statistic':protres[0],'pvalue':protres[1],'percentage annotated in group':protres[2],'percentage annotated in rest':protres[3]})

output=pd.DataFrame(output)
output['FDR']=multitest.multipletests(output['pvalue'],method='fdr_bh')[1]
print(output[output['FDR']<0.05])
output.to_csv('FisherTests.output.tsv',sep='\t')
