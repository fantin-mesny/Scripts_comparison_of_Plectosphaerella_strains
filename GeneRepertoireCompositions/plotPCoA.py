import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

metadata=pd.read_csv('metadata_phyloAndHost.csv').set_index('Fungus')
phylo_palette={'subsp1':'#0000CD','subsp2':'#1E90FF','notInSubsp':'#8E00F6','P. delsorboi':'#DF6767','P. plurivora':'#38771D'}
host_palette={'Brassicaceae':'#6699FF',	'Not plant isolated':'#BFBFBF',	'Other dicots':'#79d279','Solanaceae':'#FF6666','Cucurbitaceae':'#FFBF80','Monocots':'#FFFF66','Rosaceae':'#D9B3FF','Fabaceae':'#D9B38C','Marchantiaceae':'#999966','Unknown':'#FFFFFF'}


for dataset in ['only69Pcucumerina','all72Plectosphaerella']:
	N=0
	fig,ax=plt.subplots(2,3,figsize=(18,7))
	for cat in ['CAZymes','Proteases','Effectors']:
		pca=pd.read_csv(cat+'.'+dataset+'.pcoaCoordinates.csv').set_index('Unnamed: 0')
		pca=pca.merge(metadata,left_index=True,right_index=True,how='left')
		sns.scatterplot(x=pca.columns[0],y=pca.columns[1],data=pca,hue='PhyloGroup',palette=phylo_palette,ax=ax[0][N],legend=None,s=50,edgecolor='black')
		ax[0][N].set_title(cat)
		sns.scatterplot(x=pca.columns[0],y=pca.columns[1],data=pca,hue='HostGroup',palette=host_palette,ax=ax[1][N],legend=None,s=50,edgecolor='black')
		N+=1
	plt.savefig('CAZymesProteasesEffectors.'+dataset+'.pcoa.pdf')
	plt.close()

	fig,ax=plt.subplots(2,1,figsize=(6,7),sharex=True)
	pca=pd.read_csv('Orthogroups.'+dataset+'.pcoaCoordinates.csv').set_index('Unnamed: 0')
	pca=pca.merge(metadata,left_index=True,right_index=True,how='left')
	sns.scatterplot(x=pca.columns[0],y=pca.columns[1],data=pca,hue='PhyloGroup',palette=phylo_palette,ax=ax[0],legend=None,s=50,edgecolor='black')
	ax[0].set_title('Orthogroups')
	sns.scatterplot(x=pca.columns[0],y=pca.columns[1],data=pca,hue='HostGroup',palette=host_palette,ax=ax[1],legend=None,s=50,edgecolor='black')
	plt.savefig('Orthogroups.'+dataset+'.pdf')
	plt.close()
