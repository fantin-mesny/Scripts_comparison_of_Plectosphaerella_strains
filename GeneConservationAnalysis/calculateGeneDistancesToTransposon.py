import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns

def parseStatTable(fileName):
	mappingDic={}
	with open(fileName,'r') as inp:
		for line in inp:
			if line.startswith('seq'):
				mappingDic[line.split('\t')[0]]=line.split('\t')[1]
	return mappingDic

	
def getDistanceToClosestTransposon(transposonDf,geneStart,geneEnd):
	closestValueToStart=min(list(transposonDf[3])+list(transposonDf[4]), key=lambda x:abs(x-geneStart))
	distanceFromStart=abs(geneStart-closestValueToStart)
	closestValueToEnd=min(list(transposonDf[3])+list(transposonDf[4]), key=lambda x:abs(x-geneEnd))
	distanceFromEnd=abs(geneEnd-closestValueToEnd)
	return min([distanceFromStart,distanceFromEnd])


## For each genome, calculate the distance of each gene to the nearest transposon (in the same contig)
distanceToTransposons=[]	
for gff in os.listdir('data_gff'):
	p=gff.split('.')[0]

	genes=pd.read_csv('data_gff/'+gff,sep='\t',comment='#',header=None)
	transposons=pd.read_csv('../Transposons/'+p+'.FinalAnnotations_Transposons.gff3',sep='\t',comment='#',header=None)
	genes=genes[genes[2]=='gene']
	transposons[0]=transposons[0].map(parseStatTable('../Transposons/transposon_annotation/'+p+'.Statistics_FinalAnnotations.txt'))

	df=pd.concat([transposons,genes]).sort_values(by=[0,3]).reset_index(drop=True)
	
	for ind in df[df[2]=='gene'].index:
		transposonsInContig=df[(df[2]=='transposon') & (df[0]==df.loc[ind,0])]
		if len(transposonsInContig)>0:
			distanceToTransposons.append({'strain':p,'gene':df.loc[ind,8],'distance to transposon':getDistanceToClosestTransposon(transposonsInContig,df.loc[ind,3],df.loc[ind,4])})
		else:
			distanceToTransposons.append({'strain':p,'gene':df.loc[ind,8],'distance to transposon':np.nan})
distanceToTransposons=pd.DataFrame(distanceToTransposons)
distanceToTransposons.to_csv('data_distanceToTransposons.csv')

## Parse orthogroups to classify them by conservation levels
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

## distances to transposon within each orthogroup
meanGeneDistanceToTransposonPerOG={}
OGtoConservation={}
for o in og.index:
	if o in coreOG:
		OGtoConservation[o]='core'
	elif o in broadlyConservedOG:
		OGtoConservation[o]='broadlyConserved'
	elif o in lowConservationOG:
		OGtoConservation[o]='lowConservation'
	elif o in sparselyDistributedOG:
		OGtoConservation[o]='sparselyDistributedOG'
	elif o in subsp1OG:
		OGtoConservation[o]='subsp1'
	elif o in subsp2OG:
		OGtoConservation[o]='subsp2'


	distancesInOG=[]
	for p in og.columns: #iterate over each protein in each orthogroup to get individual distances and calculate a mean
		distanceToTransposons_p=distanceToTransposons[distanceToTransposons['strain']==p].set_index('gene')
		for prot in og.loc[o,p].split(', '):
			if prot!='':
				if p not in ['P16','P117']:
					protID='ID='+prot.split('|')[-1]
					distancesInOG.append(distanceToTransposons_p.loc[protID,'distance to transposon'])
				else:
					protID=';proteinId='+prot.split('|')[2]+';'
					if len(distanceToTransposons_p[distanceToTransposons_p.index.str.contains(protID)])!=1:
						print('ERROR',p,protID,distanceToTransposons_p[distanceToTransposons_p.index.str.contains(protID)])
					else:
						distancesInOG.append(list(distanceToTransposons_p[distanceToTransposons_p.index.str.contains(protID)]['distance to transposon'])[0])
	meanGeneDistanceToTransposonPerOG[o]=pd.Series(distancesInOG).dropna().mean()


meanGeneDistanceToTransposonPerOG_df=pd.DataFrame(index=og.index,columns=['mean distance to transposon','conservation'])
meanGeneDistanceToTransposonPerOG_df['mean distance to transposon']=meanGeneDistanceToTransposonPerOG_df.index.map(meanGeneDistanceToTransposonPerOG)
meanGeneDistanceToTransposonPerOG_df['conservation']=meanGeneDistanceToTransposonPerOG_df.index.map(OGtoConservation)
meanGeneDistanceToTransposonPerOG_df.to_csv('data_meanGeneDistanceToTransposonPerOG.csv')
				

eff=pd.read_csv('data_annotationEffectors.csv')
caz=pd.read_csv('data_annotationCazymes.csv')
prot=pd.read_csv('data_annotationProteases.csv')

meanGeneDistanceToTransposonPerOG_df['CAZymes']=meanGeneDistanceToTransposonPerOG_df['Orthogroup'].isin(caz['Orthogroup'])
meanGeneDistanceToTransposonPerOG_df['Effectors']=meanGeneDistanceToTransposonPerOG_df['Orthogroup'].isin(eff['Orthogroup'])
meanGeneDistanceToTransposonPerOG_df['Proteases']=meanGeneDistanceToTransposonPerOG_df['Orthogroup'].isin(eff['Orthogroup'])
meanGeneDistanceToTransposonPerOG_df['CAZymes']=meanGeneDistanceToTransposonPerOG_df['CAZymes'].map({True:'CAZymes',False:'Others'})
meanGeneDistanceToTransposonPerOG_df['Effectors']=meanGeneDistanceToTransposonPerOG_df['Effectors'].map({True:'Effectors',False:'Others'})
meanGeneDistanceToTransposonPerOG_df['Proteases']=meanGeneDistanceToTransposonPerOG_df['Proteases'].map({True:'Proteases',False:'Others'})
df2=pd.concat([meanGeneDistanceToTransposonPerOG_df[meanGeneDistanceToTransposonPerOG_df['CAZymes']=='CAZymes'].rename(columns={'CAZymes':'Group'})[['Orthogroup','mean distance to transposon','conservation','Group']],meanGeneDistanceToTransposonPerOG_df[meanGeneDistanceToTransposonPerOG_df['Proteases']=='Proteases'].rename(columns={'Proteases':'Group'})[['Orthogroup','mean distance to transposon','conservation','Group']],meanGeneDistanceToTransposonPerOG_df[meanGeneDistanceToTransposonPerOG_df['Effectors']=='Effectors'].rename(columns={'Effectors':'Group'})[['Orthogroup','mean distance to transposon','conservation','Group']],meanGeneDistanceToTransposonPerOG_df[(meanGeneDistanceToTransposonPerOG_df['CAZymes']=='Others') & (meanGeneDistanceToTransposonPerOG_df['Proteases']=='Others') & (meanGeneDistanceToTransposonPerOG_df['Effectors']=='Others')].rename(columns={'Effectors':'Group'})[['Orthogroup','mean distance to transposon','conservation','Group']]])

fig,ax=plt.subplots(2,1,figsize=(10,7),sharex=True)
sns.boxplot(y='conservation',x='mean distance to transposon',data=meanGeneDistanceToTransposonPerOG_df,ax=ax[0],order=['core','broadlyConserved','sparselyDistributedOG','subsp1','subsp2','lowConservation'],fliersize=2,flierprops={"marker": "o"})
sns.boxplot(y='Group',x='mean distance to transposon',data=df2,ax=ax[1],fliersize=2,flierprops={"marker": "o"},order=['CAZymes','Proteases','Effectors','Others'])
ax[0].set_xscale('log') 
ax[1].set_xscale('log') 

plt.savefig('meanGeneDistanceToTransposonPerOG_boxplot.pdf')

plt.close()