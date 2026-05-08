library('vegan')
library('RVAideMemoire')
library('multcompView')

## Get letters indicating significant differences between groups according to pairwise PERMANOVA
getLettersPPM<- function(ppm){
  p_mat=ppm$p.value
  lt_indices <-!is.na(p_mat)
  p_vector <- p_mat[lt_indices]
  names(p_vector) <- paste(rownames(p_mat)[row(p_mat)[lt_indices]], 
                          colnames(p_mat)[col(p_mat)[lt_indices]], sep = "-")
  print(multcompLetters(p_vector, compare = "<", threshold = 0.05))
}

## Parse metadata
meta0<-read.csv('metadata_phyloAndHost.csv',row.names='Fungus')
meta0 <- with(meta0, meta0[order(rownames(meta0)) , ])


catalogs<-c('CAZymes','Effectors','Proteases')
for (dataset in c('only69Pcucumerina','all72Plectosphaerella')){ # Each analysis will be run once on P.cucumerina-only, once on all Plectosphaerella genomes
    ## Run analysis for CAZymes, effectors and proteases
    for (cat in catalogs) {
        csv<-read.csv(paste('data_',cat,'.annotationTable.csv',sep=''),row.names='Fungus') #parse annotation table
        if (dataset=='only69Pcucumerina'){
            csv <- csv[!(row.names(csv) %in% c("P64","P117","P30")),]
            meta <- meta0[!(row.names(meta0) %in% c("P64","P117","P30")),]
        }
        else {
           meta <- meta0
        }
        csv <- with(csv,  csv[order(rownames(csv)) , ])
        dm<-vegdist(csv, method='jaccard') # Jaccard distance matrix
        pca <- prcomp(dm) # PCoA on Jaccard distances
        PCA <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2])
        pc1Name<-paste("PC1 (",as.character(round(summary(pca)$importance[2,'PC1']*100,digits=2)),"%)",sep='')
        pc2Name<-paste("PC2 (",as.character(round(summary(pca)$importance[2,'PC2']*100,digits=2)),"%)",sep='')
        names(PCA)[1] <- pc1Name
        names(PCA)[2] <- pc2Name
        write.csv(PCA, paste(cat,dataset,'pcoaCoordinates.csv',sep='.'))
        perm <- adonis2(dm~PhyloGroup+HostGroup, data=meta, permutations = 9999,by = "terms") # PERMANOVA test
        capture.output(perm, file=paste(cat,dataset,'permanova.txt',sep='.'))
        ppm=pairwise.perm.manova(dm, fact=meta$PhyloGroup) # Post-hoc pairwise PERMANOVA on phylogenetic group
        getLettersPPM(ppm)
        capture.output(ppm, file=paste(cat,dataset,'pairwisePermManova_phylo.txt',sep='.'))
        ppm=pairwise.perm.manova(dm, fact=meta$HostGroup) # Post-hoc pairwise PERMANOVA on host phylogenetic group
        getLettersPPM(ppm)
        capture.output(ppm, file=paste(cat,dataset,'pairwisePermManova_host.txt',sep=''))
    }
    ## Run analysis for orthogroups
    og<- read.csv('data_Orthogroups.GeneCount.tsv',row.names='Orthogroup',sep='\t')
    og <- og[, !names(og) %in% c("Total")]
    og<-t(og) #NOTE: rows are in the same order as meta already
    if (dataset=='only69Pcucumerina'){
        meta <- meta0[!(row.names(meta0) %in% c("P64","P117","P30")),]
        og <- og[!(row.names(og) %in% c("P64","P117","P30")),]
    } else {meta <- meta0}
    dm<-vegdist(og, method='jaccard') # Jaccard distance matrix
    pca <- prcomp(dm) # PCoA on Jaccard distances
    PCA <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2])
    pc1Name<-paste("PC1 (",as.character(round(summary(pca)$importance[2,'PC1']*100,digits=2)),"%)",sep='')
    pc2Name<-paste("PC2 (",as.character(round(summary(pca)$importance[2,'PC2']*100,digits=2)),"%)",sep='')
    names(PCA)[1] <- pc1Name
    names(PCA)[2] <- pc2Name
    write.csv(PCA, paste('Orthogroups',dataset,'pcoaCoordinates.csv',sep='.'))
    perm <- adonis2(dm~PhyloGroup+HostGroup, data=meta, permutations = 9999,by = "terms") # PERMANOVA test
    capture.output(perm, file=paste('Orthogroups',dataset,'permanova.txt',sep='.'))
    ppm=pairwise.perm.manova(dm, fact=meta$PhyloGroup)
    getLettersPPM(ppm)
    capture.output(ppm, file=paste('Orthogroups',dataset,'pairwisePermManova_phylo.txt',sep='.')) # Post-hoc pairwise PERMANOVA
    ppm=pairwise.perm.manova(dm, fact=meta$HostGroup)
    getLettersPPM(ppm)
    capture.output(ppm, file=paste('Orthogroups',dataset,'pairwisePermManova_host.txt',sep='')) # Post-hoc pairwise PERMANOVA

}