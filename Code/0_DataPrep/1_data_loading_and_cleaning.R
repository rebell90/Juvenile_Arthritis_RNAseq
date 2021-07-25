
####### Load Packages -- includes automatic installation from CRAN/Bioconductor if package is not yet installed
#install.packages("pacman")
pacman::p_load(scRNAseq, Glimma, DropletUtils,iSEE,BiocParallel,mclust,dplyr,biomaRt, scater,robustbase, scran, BiocNeighbors,
               Seurat, patchwork,htmltools, raster, readr, data.table, stringr, ggplot2, SingleCellExperiment, sctransform, 
               AnnotationHub, limma, edgeR, tidyverse, cowplot, Matrix.utils, magrittr, Matrix, purrr, reshape2, S4Vectors, zoo,
               tibble, clustree, GenomicRanges, EnsDb.Hsapiens.v75, pheatmap, apeglm, png, DESeq2, ashr, sva, PCAtools, GEOquery, RColorBrewer, vsn, PoiClaClu,
               EnrichmentBrowser, regioneR, airway, ALL, hgu95av2.db, Hmisc, Gviz, BSgenome.Hsapiens.UCSC.hg19.masked, tximport, OUTRIDER, EDASeq) 

gc() #garbage cleanup
## Juvenile Idiopathic Arthritis(formerly known as Juvenile Rheumatoid Arthritis) is  a rare autoimmune 
## disorder in children that has similarities to adult RA.  
## a risk of JIA, especially systemic JIA, is macrophage activation syndrome , which is a 
## life-threatening episode of hyperinflammation driven by IFN-y.

## JIA is idiopathic, but there are genetic links, and it's thought to be a combination of genetics, 
## environmental factors, and triggering events. 
## Even within JIA, there are many sub-types , and unique patterns disease expression.

## This study looks at the molecular profile of Active sJIA (MAS, New Onset JIA) and Inactive JIA, 
## along with healthy controls 

## For more information on Juvenile Idiopathic Arthritis, see resources.

## GEO Link
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147608


## reference values for na handling in metadata
# https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781444345186.app2

######### systemic JIA #################
## download Raw files from geo via url and unzip
url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE147608&format=file"
utils::download.file(url, destfile = "GSE147608_RAW.tar", mode = "wb")
utils::untar("GSE147608_RAW.tar", exdir = ".")
files_gz <- Sys.glob("GSM*.tsv.gz")
files_txt <- for(f in files_gz)
  R.utils::gunzip(f, overwrite = TRUE)

## store list of file-names in an object
files_txt <- c( "GSM4435365_run0282_lane6_read1_index702_09FLW.tsv", 
                "GSM4435366_run0282_lane7_read1_index705_09F9D.tsv",
                "GSM4435367_run0294_lane1_read1_index704_09FEM.tsv", 
                "GSM4435368_run0297_lane1_read1_index709_09FK4.tsv",
                "GSM4435369_run0297_lane1_read1_index710_09FKF.tsv", 
                "GSM4435370_run0297_lane1_read1_index711_09FKG.tsv",
                "GSM4435371_run0297_lane2_read1_index701_09FLC.tsv", 
                "GSM4435372_run0818_lane1_read1_indexN706-S502_A_409_FLAR_09FSD_15073.tsv",
                "GSM4435373_run0282_lane6_read1_index703_09F9R.tsv", 
                "GSM4435374_run0294_lane2_read1_index708_09FJQ.tsv", 
                "GSM4435375_run0297_lane1_read1_index712_09FL1.tsv",
                "GSM4435376_run0818_lane1_read1_indexN702-S506_A_409_FLAR_09FMB_15069.tsv",
                "GSM4435377_run0828_lane2_read1_indexN703-S507_A_409_FLAR_09FNK_15070.tsv", 
                "GSM4435378_run0282_lane7_read1_index708_09FB5.tsv",
                "GSM4435379_run0282_lane8_read1_index710_09FBY.tsv", 
                "GSM4435380_run0294_lane1_read1_index701_09FCT.tsv",
                "GSM4435381_run1808_lane12_read1_indexN702-S502_C6_82_16232.tsv",
                "GSM4435382_run1810_lane12_read1_indexN711-S502_C14_90_16240.tsv", 
                "GSM4435383_run1810_lane12_read1_indexN712-S503_C15_91_2_16241.tsv",
                "GSM4435384_run1808_lane12_read1_indexN703-S503_C7_83_16233.tsv",
                "GSM4435385_run1808_lane12_read1_indexN704-S504_C8_84_16234.tsv", 
                "GSM4435386_run1810_lane12_read1_indexN701-S504_C16_US_16242.tsv", 
                "GSM4435387_run1810_lane12_read1_indexN705-S505_C9_85_16235.tsv",
                "GSM4435388_run1810_lane12_read1_indexN706-S506_C10_86_16236.tsv",
                "GSM4435389_run1810_lane12_read1_indexN707-S507_C11_87_16237.tsv", 
                "GSM4435390_run1810_lane12_read1_indexN708-S508_C12_88_16238.tsv",
                "GSM4435391_run1810_lane12_read1_indexN710-S517_C13_89_16239.tsv", 
                "GSM4435392_run0294_lane1_read1_index702_09FD6.tsv", 
                "GSM4435393_run0297_lane2_read1_index702_09FLV.tsv",
                "GSM4435394_run0818_lane1_read1_indexN705-S517_A_409_FLAR_09FRY_15072.tsv", 
                "GSM4435395_run0818_lane2_read1_indexN708-S504_A_409_FLAR_09FT6_15075.tsv", 
                "GSM4435396_run0818_lane2_read1_indexN709-S505_A_409_FLAR_09FTG_15076.tsv", 
                "GSM4435397_run0282_lane6_read1_index704_09FAK.tsv", 
                "GSM4435398_run0282_lane7_read1_index706_09FAL.tsv", 
                "GSM4435399_run0282_lane7_read1_index707_09FAM.tsv", 
                "GSM4435400_run0282_lane8_read1_index709_09FK9.tsv", 
                "GSM4435401_run0282_lane8_read1_index711_09FC6.tsv")

## rename files to only be sample name
names(files_txt) <- str_sub(files_txt, 1, 10)

# counts are on the transcript level, and must be converted to the gene-level 
#(easier to interpret and most algorithms are designed for gene-level data)
# cannot simply sum transcripts by gene or "rollup" per say.  
# Overlaps, lengths, locations, lib size, etc all factor in to determine the gene count 
# relative to the individual transcript counts within each gene


#### GENE DATA PREP/ANNOTATION #####
library(AnnotationHub)
ah <- AnnotationHub() ## load annotationHub data into object
query(ah, c("EnsDb", "Hsapiens", "v94"))  ## query to retrieve Ensembl data
ensembl94 <-  ah[["AH64923"]] ## record from  query-result 
columns(ensembl94) ## list of columns 

keytypes(ensembl94)

## create key for transcript id
k <- keys(ensembl94, keytype = "TXID")

## create object that links transcript data to gene data
tx2gene <- biomaRt::select(ensembl94, k, "GENEID", "TXID")
txi.kallisto.tsv <- tximport(files_txt, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion=TRUE)

## count object
JIA_counts <- txi.kallisto.tsv$counts
JIA_counts <- round(JIA_counts, 0)

## gene anno dataframe 
gene_id <- rownames(JIA_counts)
gene_anno <- biomaRt::select(ensembl94,
                             keys= gene_id,
                             columns=c("GENEID", "TXID", "SYMBOL", "DESCRIPTION", "GENEBIOTYPE"),
                             keytype="GENEID")

## Granges list object (rowRanges )
# TxByGns <- transcriptsBy(ensembl94, by = "gene", filter = TxIdFilter(tx2gene$TXID))

# remove pseudogenes, spikes, and mitochondiral genes: quality can be compromised in these types
# and/or expression is harder to accurately quantify. For simplicity, we will remove these genes

is.spike <- grepl("^ERCC", gene_anno$SYMBOL)
is.mito <- grepl("mitochon", gene_anno$DESCRIPTION)
is.pseudo <- grepl("pseudo", gene_anno$GENEBIOTYPE)

spike_genes <- gene_anno[is.spike,] # 102
mito_genes <- gene_anno[is.mito,] #2420
pseudo_genes <- gene_anno[is.pseudo,] #20582

## remove from counts
JIA_counts <- JIA_counts[which(rownames(JIA_counts) %nin% unique(spike_genes$GENEID)) , ]
JIA_counts <- JIA_counts[which(rownames(JIA_counts) %nin% unique(mito_genes$GENEID)) , ]
JIA_counts <- JIA_counts[which(rownames(JIA_counts) %nin% unique(pseudo_genes$GENEID)) , ]

## update gene_anno table 
gene_anno <- gene_anno[which(gene_anno$GENEID %in% rownames(JIA_counts)),]

# filter to only include transcripts in experiment
TxByGns_JIA <- transcriptsBy(ensembl94, by = "gene", filter = TxIdFilter(gene_anno$TXID))

## length (for normalization)
JIA_normMat <- txi.kallisto.tsv$length

## read in metadata
GSE147608_META <- read.table("jiameta.txt", sep = ",", header = TRUE)
row.names(GSE147608_META) <- GSE147608_META$Sample.Name
summary(GSE147608_META)
head(GSE147608_META)
view(GSE147608_META)

#Age -- change NAs to rounded mean (10)
mean(GSE147608_META$AGE, na.rm = TRUE)
GSE147608_META$AGE[is.na(GSE147608_META$AGE)] <- 10

#bundle New-onsets (NOS) with active and convert to Factor
GSE147608_META$SAMPLE_TYPE <- str_trim(GSE147608_META$SAMPLE_TYPE)
GSE147608_META$SAMPLE_TYPE <- as.factor(GSE147608_META$SAMPLE_TYPE)


## dplyr to look for NAs within group for measurement variables
library(dplyr)

### CRP = measures inflammation level
## change crp to numeric (remove < and convert)

GSE147608_META$crp <- gsub("<", "", GSE147608_META$crp, fixed = TRUE)
GSE147608_META$crp <- as.numeric(GSE147608_META$crp)

### NA fields in healthy controls to be populated with average/healthy ranges for children for each 
### respectiove test.
### ref data:  https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781444345186.app2
#              https://www.ucsfbenioffchildrens.org/medical-tests/esr

## create a sub-data frame of Healthy controls only, so the NA fields can be populated
control_GSE147608_META_na <- 
  GSE147608_META %>%
  mutate(
    ESR = ifelse(is.na(ESR) , 7.5, ESR),  ## average val in normal range for children
    hemoglobin = ifelse(is.na(hemoglobin) , 13.5, hemoglobin), ## average val in normal range for children
    crp = ifelse(is.na(crp)  , 0.29, crp),  ## < 1  in normal range for children, but this df uses <0.29 
    wbc = ifelse(is.na(wbc) , 7.5, wbc), ## average val in normal range for children
    Platelets = ifelse(is.na(Platelets) , 300, Platelets),  ## average val in normal range for children
    Fever = ifelse(str_trim(Fever) == '', 'N', str_trim(Fever)),
    arthritis = ifelse(str_trim(arthritis) == '', 'N', str_trim(arthritis)),
    rash = ifelse(str_trim(rash) == '', 'N', str_trim(rash)),
    hepatosplenomegaly = ifelse(str_trim(hepatosplenomegaly) == '', 'N', str_trim(hepatosplenomegaly)),
    lymphoadenopathy = ifelse(str_trim(lymphoadenopathy) == '', 'N', str_trim(lymphoadenopathy))) %>% 
  dplyr::filter(SAMPLE_TYPE == 'Control')

control_GSE147608_META_na <- distinct(control_GSE147608_META_na)

## convert NAs for non-healthy controls

na_num_cols <- c("hemoglobin", "Platelets", "crp", "ESR", "wbc")
GSE147608_META2 <- GSE147608_META %>% 
  group_by(SAMPLE_TYPE) %>% 
  mutate_at(na_num_cols, na.aggregate) %>%
  dplyr::filter(SAMPLE_TYPE != 'Control')


## categorical vars (populate empty-string vals with 'N')
GSE147608_META2 <- GSE147608_META2 %>%
  mutate(
    Fever = ifelse(str_trim(Fever) == '', 'N', str_trim(Fever)),
    arthritis = ifelse(str_trim(arthritis) == '', 'N', str_trim(arthritis)),
    rash = ifelse(str_trim(rash) == '', 'N', str_trim(rash)),
    hepatosplenomegaly = ifelse(str_trim(hepatosplenomegaly) == '', 'N', str_trim(hepatosplenomegaly)),
    lymphoadenopathy = ifelse(str_trim(lymphoadenopathy) == '', 'N', str_trim(lymphoadenopathy)))

## Combine dataframes together 
GSE147608_META <- rbind(GSE147608_META2, control_GSE147608_META_na)

GSE147608_META <- as.data.frame(GSE147608_META) ## convert to df
rownames(GSE147608_META) <- GSE147608_META$Sample.Name ## populate rownames

# check if any Nas left
sum(is.na(GSE147608_META)) ## 0 
## View metadata
head(GSE147608_META, 5)
summary(GSE147608_META)

## convert data-types
cat_vars <- c("Batch",  "sex", "Sample.Name",  "arthritis", "Fever", "hepatosplenomegaly", 
              "rash", "lymphoadenopathy")
GSE147608_META[,cat_vars] <- lapply(GSE147608_META[,cat_vars] , as.factor)
GSE147608_META<- GSE147608_META %<>% mutate_at(cat_vars, factor)

summary(GSE147608_META)
## order rownames
GSE147608_META <- GSE147608_META[ order(row.names(GSE147608_META)), ]
GSE147608_META$ifJIA <- ifelse(GSE147608_META$SAMPLE_TYPE == 'Control', 'Control', 'JIA')
GSE147608_META$ifJIA <- as.factor(GSE147608_META$ifJIA)

GSE147608_META$ifActive <- ifelse(GSE147608_META$SAMPLE_TYPE == 'Control', 'Control', 
                                  (ifelse(GSE147608_META$SAMPLE_TYPE == 'Inactive', 'Inactive', 'Active' )))
GSE147608_META$ifActive <- as.factor(GSE147608_META$ifActive)

### make clean colData dataframe 
colDataVars <- c("AGE", "sex", "SAMPLE_TYPE", "ifJIA", "ifActive", "Batch", "hemoglobin", "arthritis", "Fever", "Platelets", "lymphoadenopathy", "rash",
                 "hepatosplenomegaly", "ESR", "crp", "wbc")

colData <- GSE147608_META[,colDataVars]
colData
rownames(colData)

barplot(colSums(JIA_counts))
#GSM4435375 ## barely any counts for these two samples, so we will remove
#GSM4435393 

colData <- colData[which(rownames(colData) %nin% c("GSM4435375", "GSM4435393")), ]
## adjust counts
JIA_counts <- JIA_counts[, which(colnames(JIA_counts) %nin% c("GSM4435375", "GSM4435393"))]
