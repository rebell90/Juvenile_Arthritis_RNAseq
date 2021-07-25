

###### PERFORM DESEQ #######

## Steps in DESEQ:
## estimation in size factors( controlling/adjusting for sequencing depth in samples)
## estimation of dispersion for each gene
## fiting a GLM to the data

dds_JIA <- DESeq(dds_JIA)

## Results (DESeqDataSet)
res_JIA_inactive <- DESeq2::results(dds_JIA) ## returns the last comparison unless otherwise specified (Inactive vs control)
res_JIA_inactive

# repeat for "ifActive_Active_vs_Control" via contrast argumnet
res_JIA_active <- DESeq2::results(dds_JIA, contrast  = c("ifActive", "Active", "Control")) 
res_JIA_active

 ## coefficient names that can be extracted
resultsNames(res_JIA_inactive)
# results metadata
mcols(res_JIA, use.names = TRUE) ## meaning of the cols
# baseMean:  mean of the normalized count values divided by size factors taken over all samples
# log2FoldChange: effect size estimate (how much a genes expression has changed compared to control samples )
# lfcSE:  standard error estimate for the log2 fold change parameter
# pval: probability that the results occured by chance
# padj: adjusted p-val

## summary
summary(res_JIA_inactive)
summary(res_JIA_active)

# the model is idcenfiying a large number of genes as significanlty differentially expressed, so,
# in order to be more discerning/strict about which genes are actuall significant we could
# lower the padj parameter/threshold, and/or increase the log2 fold change threshold (via the lcfThreshold argument)

###### Inactive vs control ######
## raise log2 threshold from 0 to 1 
res_JIA_inactive_LFC1 <- DESeq2::results(dds_JIA, lfcThreshold=1, contrast = c("ifActive", "Inactive", "Control"))
table(res_JIA_inactive_LFC1$padj < 0.1) 
res_JIA_inactive_LFC1.05 <- DESeq2::results(dds_JIA, lfcThreshold=2, 
                                            contrast = c("ifActive", "Inactive", "Control"),alpha = 0.5)
table(res_JIA_inactive_LFC1.05$padj < 0.05) ## must macth alpha (False Discovery Rate) to match p val threshold

summary(res_JIA_inactive_LFC1)
summary(res_JIA_inactive_LFC1.05)

###### Active vs control ######
## raise log2 threshold from 0 to 1 
res_JIA_active_LFC1 <- DESeq2::results(dds_JIA, lfcThreshold=1, contrast = c("ifActive", "Active", "Control"))
table(res_JIA_active_LFC1$padj < 0.1) 
res_JIA_active_LFC1.05 <- DESeq2::results(dds_JIA, lfcThreshold=2, 
                                            contrast = c("ifActive", "Active", "Control"),alpha = 0.5)
table(res_JIA_active_LFC1.05$padj < 0.05) ## must macth alpha (False Discovery Rate) to match p val threshold

summary(res_JIA_active_LFC1)
summary(res_JIA_active_LFC1.05)

coef(dds_JIA)
#### shrinkage estimators (for ranking and visualization)
# apeglm: adaptive t prior shrinkage estimator
# ashr:  adaptive shrinkage estimator to fit a mixture of normal distributions to form the prior

# normal
resNorm_inactive <- lfcShrink(dds_JIA, lfcThreshold=2, coef="ifActive_Inactive_vs_Control", type="normal")
resNorm_active <- lfcShrink(dds_JIA, lfcThreshold=2, coef="ifActive_Active_vs_Control", type="normal")

# ashr (note: lfcthreshold not used)
resAshr_inactive <- lfcShrink(dds_JIA,  coef="ifActive_Inactive_vs_Control", type="ashr")
resAshr_active <- lfcShrink(dds_JIA, coef="ifActive_Active_vs_Control", type="ashr")

# apeglm
resApe_inactive <- lfcShrink(dds_JIA, lfcThreshold=2, coef="ifActive_Inactive_vs_Control", type="apeglm")
resApe_active <- lfcShrink(dds_JIA, lfcThreshold=2, coef="ifActive_Active_vs_Control", type="apeglm")



par(mfrow=c(1,3), mar=c(4,4,2,1))
xlim <- c(1,1e5); ylim <- c(-3,3)
# inactive
plotMA(resApe_inactive, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm_inactive, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAshr_inactive, xlim=xlim, ylim=ylim, main="ashr")

#active
plotMA(resApe_active, xlim=xlim, ylim=ylim, main="apeglm")
plotMA(resNorm_active, xlim=xlim, ylim=ylim, main="normal")
plotMA(resAshr_active, xlim=xlim, ylim=ylim, main="ashr")


### use apeglm for remaining analysis
resSig_active <- subset(resApe_active, lfcThreshold=2)
head(resSig_active[ order(resSig_active$log2FoldChange), ])

resSig_inactive <- subset(resApe_inactive, lfcThreshold=2)
head(resSig_inactive[ order(resSig_inactive$log2FoldChange), ])

##### PLOTS ######

hist(resApe_active$log2FoldChange[resApe_active$baseMean > 1],  col = "grey50", border = "white")

par(mar=c(8,5,2,2))
boxplot(log10(assays(dds_JIA)[["cooks"]]), range=0, las=2)

library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd_JIA)), decreasing = TRUE), 20)
mat  <- assay(vsd_JIA)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(vsd_JIA$ifActive)
rownames(anno) <- rownames(colData(vsd_JIA))
pheatmap(mat, annotation_col = anno)



## top gene objects 
topGene_active <- rownames(resApe_active)[which.max(resApe_active$log2FoldChange)]
topGene_inactive <- rownames(resApe_inactive)[which.max(resApe_inactive$log2FoldChange)]

## plot results in genomic space
resGR_inactive <- lfcShrink(dds_JIA, coef="ifActive_Inactive_vs_Control", lfcThreshold = 2,
                   type="apeglm", format="GRanges")
resGR_inactive 

gene_ens_inactive <-names(resGR_inactive)
resGR_inactive$symbol <- mapIds(org.Hs.eg.db, gene_ens_inactive, "SYMBOL", "ENSEMBL")

window_inactive <- resGR_inactive[topGene_inactive] + 1e6
strand(window_inactive) <- "*"
resGRsub_inactive <- resGR_inactive[sp::over(resGR_inactive,window_inactive)]
naOrDup_inactive <- is.na(resGRsub_inactive$symbol) | duplicated(resGRsub_inactive$symbol)
resGRsub_inactive$group <- ifelse(naOrDup_inactive, names(resGRsub_inactive), resGRsub_inactive$symbol)

status_inactive <- factor(ifelse(resGRsub_inactive$padj < 0.05 & !is.na(resGRsub_inactive$padj),
                        "sig", "notsig"))
options(ucscChromosomeNames = FALSE)
g_inactive <- GenomeAxisTrack()
a_inactive <- AnnotationTrack(resGRsub_inactive, name = "gene ranges", feature = status)
d_inactive <- DataTrack(resGRsub_inactive, data = "log2FoldChange", baseline = 0,
               type = "h", name = "log2 fold change", strand = "+")
plotTracks(list(g_inactive, d_inactive, a_inactive), groupAnnotation = "group",
           notsig = "grey", sig = "hotpink")

?`trim,GenomicRanges-method`
