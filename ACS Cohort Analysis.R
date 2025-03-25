library(GEOquery)
library(tidyverse)
library(umap)
library(DESeq2)

#ACS cohorts
#Cleanup metadata
ACS.meta <- getGEO(GEO = "GSE79362")
ACS.meta <- pData(ACS.meta[[1]])

ACS.meta.v2 <- read.csv(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/DK Postdoc Data/Coding stuff/ACS Cohort/ACS_meta_v2.csv")
ACS.meta.v2.filt <- ACS.meta.v2 %>% filter(GSM != "")
ACS.meta.v2.filt <- ACS.meta.v2.filt[1:15]

ACS.counts.sheet1 <- read_csv(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/DK Postdoc Data/Coding stuff/ACS Cohort/ACS_counts_sheet1.csv")
ACS.counts.sheet2 <- read_csv(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/DK Postdoc Data/Coding stuff/ACS Cohort/ACS_counts_sheet2.csv")
ACS.combined.counts <- full_join(ACS.counts.sheet1, ACS.counts.sheet2, by = "entry")
ACS.combined.counts <- ACS.combined.counts[c(6:270,277:366)]

#A lot of genes had multiple rows per gene so I averaged reads per gene.
ACS.combined.counts <- ACS.combined.counts %>% group_by(gene.x) %>% summarise_all("mean") %>% na.omit()
ACS.combined.counts <- as.data.frame(ACS.combined.counts, check.names=F)
row.names(ACS.combined.counts) <- ACS.combined.counts$gene.x
ACS.combined.counts <- ACS.combined.counts[2:355]

#Use meta data to filter counts tibble to match prior to running DESeq2
ACS.meta <- ACS.meta %>% filter(geo_accession %in% ACS.meta.v2.filt$GSM)
#Split into TB progressor or non progessor then make the sample_name variable then bind the rows to merge back into one dataframe as how I am exctracting the sample name needs to be different due to different variable lengths
ACS.meta.prog <- ACS.meta %>% filter(source_name_ch1 == "Blood, case")
ACS.meta.non <- ACS.meta %>% filter(source_name_ch1 == "Blood, control")
ACS.meta.prog$sample_name <- str_sub(ACS.meta.prog$title,16, end = -2)
ACS.meta.non$sample_name <- str_sub(ACS.meta.non$title,17, end = -2)
ACS.meta <- rbind(ACS.meta.prog, ACS.meta.non)

ACS.meta.name <- ACS.meta %>% select(geo_accession, sample_name)
colnames(ACS.meta.name) <- c("GSM","sample_name")
ACS.meta.v2.filt <- left_join(ACS.meta.v2.filt, ACS.meta.name, by = "GSM")
ACS.meta.v2.filt <- ACS.meta.v2.filt %>% filter(QFT != "indeterminate")

#Filter counts to only include data for which there is metadata
ACS.combined.counts <- ACS.combined.counts %>% select(ACS.meta.v2.filt$sample_name)

#Finish formatting the metadata to prepare it for DESeq2
row.names(ACS.meta.v2.filt) <- ACS.meta.v2.filt$sample_name
ACS.meta.v2.filt <- ACS.meta.v2.filt[1:15]

#Generate dds object
dds <- DESeqDataSetFromMatrix(countData = round(ACS.combined.counts),
                              colData = ACS.meta.v2.filt,
                              design = ~ Group)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# run DESeq
dds <- DESeq(dds)

#Normalize and plot PCA colored by group
dds.pca <- vst(dds)
plotPCA(dds.pca,
        intgroup = c('Group'),
        returnData = FALSE)

resultsNames(dds)

#Get results from DESeq2
res <- results(dds)
res

#Normalized Counts
ACS.counts.norm <- counts(dds, normalized = T)
ACS.counts.norm <- as.data.frame(ACS.counts.norm)

saveRDS(dds, file = '/Users/dmitrikotov/Library/CloudStorage/Box-Box/Dmitri Personal/DK Postdoc Data and Analysis/Coding stuff/ACS Cohort/ACS_progressors')

IL1RN.progressors <- plotCounts(dds, gene="IL1RN", intgroup="TimeToTB", returnData=T)
IL1RN.progressors <- IL1RN.progressors %>% filter(TimeToTB != "-")
IL1RN.progressors$TimeToTB <- as.numeric(IL1RN.progressors$TimeToTB)
ggplot(IL1RN.progressors, aes(x = TimeToTB, y = count)) + geom_point() + geom_smooth() + geom_vline(xintercept = 200, linetype = "dotted", color = "red", size = 1.5)
