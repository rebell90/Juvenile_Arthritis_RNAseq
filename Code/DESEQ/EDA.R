######################### DESEQ OBJECT ###############################
## http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#pre-filtering-the-dataset

dds_JIA <- DESeqDataSetFromMatrix(countData =  as.matrix(JIA_counts) ,
                                    colData =  colData,
                                    rowRanges = TxByGns_JIA,
                                    design =  ~ 1 + ifActive)
dds_JIA$ifActive <-relevel(dds_JIA$ifActive, ref = "Control") ## revel var
table(dds_JIA$ifActive)

## Boxplot
boxplot(log10(counts(dds_JIA)+1), las=2) 

# eliminate zeros
keep_feature <- rowSums(counts(dds_JIA) > 20 ) > 9 ## keep genes with at least 20 counts in over 9 samples across all samples
dds_JIA<-  dds_JIA[keep_feature, ]


## look at counts
summary(assay(dds_JIA)) ## some genes seem very strongly expressed across ALL samples 

## boxplot
boxplot(assay(dds_JIA), las=2) 

### mean-var plots 
msd <- meanSdPlot(counts(dds_JIA), ranks = FALSE)
msd

log.cts.one <- log2(counts(dds_JIA) + 1)
log.cts.one_msd <- meanSdPlot(log.cts.one , ranks = FALSE)
log.cts.one_msd

## variance stabilizing with counts
vsd_JIA <- DESeq2::vst(dds_JIA, blind = TRUE)
head(vsd_JIA, 3)
## box plot
par(mar=c(12,5,4,2)+0.1)
boxplot(assay(vsd_JIA), xlab="", ylab="Log2 counts per million",las=2, main="Normalized Count Distributions")
abline(h=median(assay(vsd_JIA)), col="blue")


# rlog
rld_JIA <- rlog(dds_JIA, blind = TRUE)
head(assay(rld_JIA), 3)
## box plot
par(mar=c(12,5,4,2)+0.1)
boxplot(assay(rld_JIA), xlab="", ylab="Log2 counts per million",las=2, main="RLD Normalized Count Distributions")
abline(h=median(assay(rld_JIA)), col="blue")


##### plot comparison ####
dds_JIA <- estimateSizeFactors(dds_JIA)

df_JIA <- bind_rows(
  as_tibble(log2(counts(dds_JIA, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(vsd_JIA)[, 1:2]) %>% mutate(transformation = "vst"),
  as_tibble(assay(rld_JIA)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df_JIA)[1:2] <- c("x", "y")  

lvls <- c("log2(x + 1)", "vst", "rlog")
df_JIA$transformation <- factor(df_JIA$transformation, levels=lvls)

compare_transforms <- ggplot(df_JIA, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  
compare_transforms

### R LOG LOOKS MOST NORMALIZED 