import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import pandas as pd
import numpy as np


## Parsing data:
root=pd.read_csv('data_GlobalFungi_rootSamples.tsv',sep='\t')
rootAndRhizo=pd.read_csv('data_GlobalFungi_rootAndRhizosphereSamples.txt')
rhizo=pd.read_csv('data_GlobalFungi_rhizosphereSamples.txt')
data=pd.concat([root,rootAndRhizo,rhizo])
data=data[['permanent_id','host','sample_type','latitude','longitude']].sort_values(by='host')
data=data[data['host']!='?']

## Update some species name to the actual one:
species_rename={'Lycopodium volubile':'Pseudodiphasium volubile','Populous deltoids':'Populus deltoides','Lilaea scilloides':'Triglochin scilloides','Scirpus mucronatus':'Schoenoplectiella mucronata','Lophozonia obliqua':'Nothofagus obliqua'}
for host in list(data['host']):
	if host not in species_rename:
		species_rename[host]=host
data['host_renamed']=data['host'].map(species_rename)

## Match species to family
host_phylo=pd.read_csv('data_Taxallnomy_plantTaxonomy.tsv',sep='\t',skiprows=1)
host_to_family={}
for h in set(data['host_renamed']):
	genus=h.split(' ')[0]
	if genus in list(host_phylo['genus']):
		host_to_family[h]=list(set(host_phylo[host_phylo['genus']==genus]['family']))[0]
	else:
		print(h,' ERROR - no phylogeny')
data['host_family']=data['host_renamed'].map(host_to_family)

## Order of the phylogenetic tree obtained with Taxallnomy:
order="Lycopodiaceae Selaginellaceae Alismataceae Hydrocharitaceae Juncaginaceae Potamogetonaceae Zosteraceae Magnoliaceae Orchidaceae Arecaceae Cyperaceae Poaceae Musaceae Amaranthaceae Caryophyllaceae Montiaceae Crassulaceae Aquifoliaceae Asteraceae Boraginaceae Ericaceae Theaceae Rubiaceae Lamiaceae Oleaceae Plantaginaceae Solanaceae Brassicaceae Fabaceae Betulaceae Fagaceae Juglandaceae Nothofagaceae Geraniaceae Euphorbiaceae Salicaceae Cannabaceae Urticaceae Pinaceae Cupressaceae Ophioglossaceae Gleicheniaceae Polytrichaceae".split(' ')
order_index={o:order.index(o) for o in order}

## Prepare data for plotting:
Nsamples=pd.DataFrame(data.groupby('host_family')['sample_type'].value_counts()).rename(columns={'sample_type':'N_samples'})#.reset_index(drop=False)
Nsamples=Nsamples.unstack('sample_type').fillna(0).droplevel(0, axis=1)
Nsamples['order']=Nsamples.index.map(order_index)
Nsamples=Nsamples.sort_values(by='order',ascending=False).drop(columns='order')
Nsamples=Nsamples.rename(columns={'rhizosphere soil':'rhizosphere','root + rhizosphere soil':'rhizosphere + root'})[['rhizosphere','rhizosphere + root','root']]
Nspecies=data[['host_renamed','host_family']].drop_duplicates()
Nspecies=pd.DataFrame(Nspecies['host_family'].value_counts()).rename(columns={'count':'N_species'})
print(Nspecies)
Nspecies['order']=Nspecies.index.map(order_index)
Nspecies=Nspecies.sort_values(by='order').reset_index(drop=False)
print(Nspecies)



## Plotting:
min_val, max_val = 0.3,1.0
n = 10
orig_cmap = plt.cm.Greys
colors = orig_cmap(np.linspace(min_val, max_val, n))
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("mycmap", colors)
fig,ax=plt.subplots(1,2,figsize=(5,10))
sns.scatterplot(y='host_family',x=1,size='N_species',data=Nspecies,ax=ax[0],sizes=(20,200))
Nsamples.plot(kind='barh', stacked=True,ax=ax[1],colormap=cmap)
plt.savefig('detectionInDiversePlants.pdf')
plt.close()



