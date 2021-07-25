
######## sample similarity
library(PCAtools)

## VSD
sampleDists_vsd <- dist(t(assay(vsd_JIA)))
sampleDists_vsd

sampleDistMatrix_vsd_JIA <- as.matrix(sampleDists_vsd)
rownames(sampleDistMatrix_vsd_JIA) <- paste(vsd_JIA$ifActive)
colnames(sampleDistMatrix_vsd_JIA) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "BrBG")))(255)
sampleDists_vsd_JIA_plot <- pheatmap(sampleDistMatrix_vsd_JIA,
                                      clustering_distance_rows = sampleDists_vsd,
                                      clustering_distance_cols = sampleDists_vsd,
                                      col = colors)
sampleDists_vsd_JIA_plot 

## RLOG
sampleDists_rld <- dist(t(assay(rld_JIA)))
sampleDists_rld

sampleDistMatrix_rld_JIA <- as.matrix(sampleDists_rld)
rownames(sampleDistMatrix_rld_JIA) <- paste(rld_JIA$ifActive)
colnames(sampleDistMatrix_rld_JIA) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "BrBG")))(255)
sampleDists_rld_JIA_plot <- pheatmap(sampleDistMatrix_rld_JIA,
                                     clustering_distance_rows = sampleDists_rld,
                                     clustering_distance_cols = sampleDists_rld,
                                     col = colors)
sampleDists_rld_JIA_plot 

##### PCA #####
vsd_pca_JIA_plot <- plotPCA(vsd_JIA, intgroup = "ifActive")
vsd_pca_JIA_plot

rld_pca_JIA_plot <- plotPCA(rld_JIA, intgroup = "ifActive")
rld_pca_JIA_plot



##### PCA TOOLS #####
vsd_counts <- assay(vsd_JIA)
p <- PCAtools::pca(vsd_counts, metadata = colData(dds_JIA), removeVar = 0.1)

screeplot(p, axisLabSize = 18, titleLabSize = 22)

# biplot
biplot(p, showLoadings = TRUE, lab = NULL)
biplot(p, showLoadings = TRUE,
       labSize = 5, pointSize = 5, sizeLoadingsNames = 5)

# pairs plot
pairsplot(p)

# loadings plot
plotloadings(p, labSize = 3)

#eigen plot
eigencorplot(p,
             metavars = colDataVars)

horn <- parallelPCA(vsd_counts)
horn$n
elbow <- findElbowPoint(p$variance)
elbow
