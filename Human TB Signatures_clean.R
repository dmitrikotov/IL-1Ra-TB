library(GEOquery)
library(tidyverse)
library(umap)
library(DESeq2)
library(readxl)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(sigQC)
library(sjmisc)
library(pROC)
library(hacksig)

#Human data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114192
#TB = active TB, IH = intermediate hyperglycemia, DM = diabetes

#Cleanup metadata
human.tb.db.meta <- getGEO(GEO = "GSE114192",GSEMatrix = TRUE)
human.tb.db.meta <- pData(human.tb.db.meta[[1]])
row.names(human.tb.db.meta) <- human.tb.db.meta$description
rownames(human.tb.db.meta)[rownames(human.tb.db.meta) == c("RSEQ341_p1","RSEQ460_p1")] <- c("RSEQ341","RSEQ460")
colnames(human.tb.db.meta)[colnames(human.tb.db.meta) == 'disease state (disease_category):ch1'] <- "disease.state"
human.tb.db.meta$TB.status <- ifelse(human.tb.db.meta$disease.state == "TB_DM" | human.tb.db.meta$disease.state == "TB_only" | human.tb.db.meta$disease.state == "TB_IH" , "TB","Uninfected")

#Clean up counts
folderfiles <- list.files("/Users/dmitrikotov/Downloads/GSE114192_RAW", pattern = "*.txt.gz", full.names = TRUE)
data_csv <- lapply(folderfiles, function(x){
  read.delim(file = x, header = FALSE, row.names = 1)
})
data <- as.data.frame(data_csv) %>% slice_head(n = 63677)
sample.names <- substr(folderfiles,55,61)
colnames(data) <- sample.names

#Generate dds object
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = human.tb.db.meta,
                              design = ~ disease.state)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$disease.state <- relevel(dds$disease.state, ref = "Healthy_Control")

# run DESeq
dds <- DESeq(dds)

saveRDS(dds, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Dmitri Data/Coding stuff/Human TB Gene Signatures/dds.TBDM")

#Normalize and plot PCA colored by group
dds.pca <- vst(dds)
plotPCA(dds.pca,
        intgroup = c('disease.state'),
        returnData = FALSE)
plotPCA(dds.pca,
        intgroup = c('TB.status'),
        returnData = FALSE)

#Normalized Counts
tb.counts.norm <- counts(dds, normalized = T)
tb.counts.norm <- as.data.frame(tb.counts.norm)

#Graph plots prettier just by entering ENSMUSG id
graph <- function(x){
  d <- plotCounts(dds, gene= x, intgroup="TB.status", returnData=TRUE)
  ggplot(d, aes(x = TB.status, y = count, color = TB.status)) + 
    geom_boxplot(coef=0, outlier.shape = NA) +
    geom_jitter(size = 3, alpha = 0.9) +
    theme_bw() +
    theme(
      legend.position = "none",
      plot.title = element_text(size=11)
    ) +
    ggtitle(x) +
    xlab("")
}

#Example for plotting graph
#Sp140
graph("ENSG00000079263")

#ISG20
graph("ENSG00000172183")

#IFITM2
graph("ENSG00000185201")

#IL1RN
graph("ENSG00000136689")


#DEG analysis

#Generate dds object
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = human.tb.db.meta,
                              design = ~ TB.status)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds$TB.status <- relevel(dds$TB.status, ref = "Uninfected")

# run DESeq
dds <- DESeq(dds)

res <- results(dds) #IL1RN TB vs uninfected padj is 1.17683e-21 using Wald test p-value
res.v2 <- results(dds, contrast = c("TB.status","TB","Uninfected"))

#TB only is 45 samples
#Healthy control is 35 samples

#convert ENSG gene names to gene symbol
tb.counts.norm$gene <- mapIds(org.Hs.eg.db, row.names(tb.counts.norm), keytype="ENSEMBL", column="SYMBOL", multiVals = "first")
tb.counts.norm <- tb.counts.norm[!is.na(tb.counts.norm$gene),]
tb.counts.norm$gene <- make.unique(tb.counts.norm$gene, sep="-")
row.names(tb.counts.norm) <- tb.counts.norm$gene

#Get IDs for the two classes to compare
TB.infected <- row.names(human.tb.db.meta[human.tb.db.meta$TB.status == "TB",])
Uninfected <- row.names(human.tb.db.meta[human.tb.db.meta$TB.status == "Uninfected",])

#Flip dataframe with normalized counts
counts <- rotate_df(tb.counts.norm)
counts <- counts[1:249,]
counts$sample <- row.names(counts) #23687
counts$Class <- ifelse(grepl(paste(TB.infected, collapse = "|"), counts$sample), "A","B")

#Import signatures - these are previously published by our group https://doi.org/10.1016/j.cell.2023.11.002
human.ifnb <- read.csv(file="/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/B6 vs Sp140 scRNA-seq 10-9-20/Plots 230320/human_ifnb_sig.csv")
human.ifng <- read.csv(file="/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/B6 vs Sp140 scRNA-seq 10-9-20/Plots 230320/human_ifng_sig.csv")
mouse.ifnb <- read.csv(file="/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/B6 vs Sp140 scRNA-seq 10-9-20/Plots 230320/mouse_ifnb_sig.csv")
mouse.ifng <- read.csv(file="/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/B6 vs Sp140 scRNA-seq 10-9-20/Plots 230320/mouse_ifng_sig.csv")

#Filter to genes in signature and have column 1 be sample name and column 2 be class
ifnb.score <- counts %>% dplyr::select("sample","Class",intersect(human.ifnb$x, colnames(counts)))
ifng.score <- counts %>% dplyr::select("sample","Class",intersect(human.ifng$x, colnames(counts)))
mouse.ifnb.score <- counts %>% dplyr::select("sample","Class",intersect(mouse.ifnb$x, colnames(counts)))
mouse.ifng.score <- counts %>% dplyr::select("sample","Class",intersect(mouse.ifng$x, colnames(counts)))

#Try roc for just IL1RN
IL1RN.score <- counts %>% dplyr::select(c("sample","Class","IL1RN"))
IL1RN.score$status <- ifelse(grepl(paste(TB.infected, collapse = "|"), IL1RN.score$sample), "TB","Uninfected")
IL1RN.score$IL1RN <- as.numeric(IL1RN.score$IL1RN)
IL1RN.roc <- roc(IL1RN.score$status, IL1RN.score$IL1RN)
plot(IL1RN.roc, print.auc = TRUE, col = "red", print.auc.y =0.47)


#gene symbol = row names and sample IDs = column names
tb.counts.norm <- tb.counts.norm[,1:249]
check_sig(tb.counts.norm, list(mouse.ifng$x))

#Methods to use for classifying based on score - "zscore","ssgsea","singscore"
scored.tb <- hack_sig(tb.counts.norm, signatures = c(list(human.ifng$x),list(human.ifnb$x),list(mouse.ifng$x),list(mouse.ifnb$x)), method = "zscore")
scored.tb$status <- ifelse(grepl(paste(TB.infected, collapse = "|"), scored.tb$sample_id), "TB","Uninfected")
ifng.roc <- roc(scored.tb$status, scored.tb$sig1, levels = c("TB","Uninfected")) #IFNg zscore AUC 83.8 ssgsea AUC 91.6 singscore AUC 88.7
ifnb.roc <- roc(scored.tb$status, scored.tb$sig2) #IFNb zscore AUC 87.4 ssgsea AUC 81.8 singscore AUC 77.5
plot(ifng.roc, print.auc = TRUE, col = "red", print.auc.y =0.47)
plot(ifnb.roc, add=TRUE, print.auc = TRUE, col = "blue",  print.auc.y =0.4)

mouse.ifng.roc <- roc(scored.tb$status, scored.tb$sig3, levels = c("TB","Uninfected")) #Mouse IFNg zscore AUC 76.9
mouse.ifnb.roc <- roc(scored.tb$status, scored.tb$sig4) #Mouse IFNb zscore AUC 79.7
plot(mouse.ifng.roc, print.auc = TRUE, col = "red", print.auc.y =0.47)
plot(mouse.ifnb.roc, add=TRUE, print.auc = TRUE, col = "blue",  print.auc.y =0.4)

#Plot human and mouse
plot(ifng.roc, print.auc = TRUE, col = "red", print.auc.y =0.47)
plot(ifnb.roc, add=TRUE, print.auc = TRUE, col = "blue",  print.auc.y =0.4)
plot(mouse.ifng.roc,add = TRUE, print.auc = TRUE, col = "chartreuse", print.auc.y =0.33)
plot(mouse.ifnb.roc, add=TRUE, print.auc = TRUE, col = "cadetblue2",  print.auc.y =0.27)

#Import human PBMC single cell data from Hao et al. (2021) Cell Satija lab
library(Seurat)
pbmc <- readRDS(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Human TB Gene Signatures/human PBMC.rds")
Idents(pbmc) <- "celltype.l2"
DimPlot(pbmc, reduction = "wnn.umap", raster = FALSE, label = TRUE)
subset.genes <- FindAllMarkers(pbmc, max.cells.per.ident = 1000)

#Find genes specific for each major cell type
cell.types <- unique(pbmc$celltype.l2)
pbmc.simplified <- subset(pbmc, idents = cell.types[!cell.types %in% c("ASDC","dnT", "Doublet","CD4 Proliferating","CD8 Proliferating","NK Proliferating")])
DimPlot(pbmc.simplified, reduction = "wnn.umap", raster = FALSE)
subset.genes.simplified <- FindAllMarkers(pbmc.simplified, max.cells.per.ident = 1000)

#change Ensembl ids to gene symbols
library(org.Hs.eg.db)
library(AnnotationDbi)

#convert ENSG gene names to gene symbol
subset.genes.simplified$symbol <- mapIds(org.Hs.eg.db, subset.genes.simplified$gene, keytype="ENSEMBL", column="SYMBOL", multiVals = "first")
subset.genes.simplified <- subset.genes.simplified[!is.na(subset.genes.simplified$symbol),]
pop.markers <- subset.genes.simplified[subset.genes.simplified$avg_log2FC > 1,]
subset.genes.simplified$symbol <- make.unique(subset.genes.simplified$symbol, sep="-")

#Score Immunosuppressive genes - complete signature AUC = 0.815 with zscore, 0.896 with ssgsea, 0.799 with singscore, 0.896 with default
#Genes impact (AUC without this gene) IL1RN = 0.789, IL18BP = 0.872, CD274 = 0.750, TGFB1 = 0.807, IDO1 = 0.818, IL10 = 0.802, ARG1 = 0.802
#"zscore","ssgsea","singscore"
immunosuppresion <- list(c("IL1RN","IL18BP","CD274","TGFB1","IDO1","IL10","ARG1"))
check_sig(tb.counts.norm, immunosuppresion)
immunosuppresion.scored <- hack_sig(tb.counts.norm, signatures = immunosuppresion, method = "ssgsea")
immunosuppresion.scored$status <- ifelse(grepl(paste(TB.infected, collapse = "|"), immunosuppresion.scored$sample_id), "TB","Uninfected")
immunosuppresion.roc <- roc(immunosuppresion.scored$status, immunosuppresion.scored$sig1, levels = c("TB","Uninfected")) #IFNg zscore AUC 83.8 ssgsea AUC 91.6 singscore AUC 88.7
plot(immunosuppresion.roc, print.auc = TRUE, col = "red", print.auc.y =0.47)

#Apply mouse signatures to human dataset
B6.strict <- read.csv(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/B6 Strict Genes.csv") %>% filter(change == "up")
Sp140.strict <- read.csv(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/Sp140 Strict Genes.csv") %>% filter(change == "up")
Irg1.strict <- read.csv(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/Irg1 Strict Genes.csv") %>% filter(change == "up")
iNOS.strict <- read.csv(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/iNOS Strict Genes.csv") %>% filter(change == "up")

#Getting rid of duplicate human genes
Sp140.strict <- filter(Sp140.strict, !(X %in% c("H2-Q6","Ifit3b")))
Irg1.strict <- filter(Irg1.strict, !(X %in% c("Gbp4","Gbp8")))
iNOS.strict <- filter(iNOS.strict, !(X %in% c("Gbp4","Gbp8")))

#HeatMap of Colors
B6.strict <- B6.strict %>% filter(avg_log2FC > 3)  %>%  arrange(avg_log2FC) %>% dplyr::select(avg_log2FC, gene)
Sp140.strict <- Sp140.strict %>% filter(avg_log2FC > 3)  %>%  arrange(avg_log2FC) %>% dplyr::select(avg_log2FC, gene)
Irg1.strict <- Irg1.strict %>% filter(avg_log2FC > 3)  %>%  arrange(avg_log2FC) %>% dplyr::select(avg_log2FC, gene)
iNOS.strict <- iNOS.strict %>% filter(avg_log2FC > 3)  %>%  arrange(avg_log2FC) %>% dplyr::select(avg_log2FC, gene)

length(unique(c(Sp140.strict$gene,Irg1.strict$gene,iNOS.strict$gene,B6.strict$gene))) #43
all.heatmap.genes <- data.frame(avg_log2FC = rep(-1, 43), gene = unique(c(Sp140.strict$gene,Irg1.strict$gene,iNOS.strict$gene,B6.strict$gene)))

B6.strict <- rbind(B6.strict, filter(all.heatmap.genes, gene %in% setdiff(all.heatmap.genes$gene, B6.strict$gene))) %>%  arrange(avg_log2FC)
B6.strict$genotype <- "B6"
Sp140.strict <- rbind(Sp140.strict, filter(all.heatmap.genes, gene %in% setdiff(all.heatmap.genes$gene, Sp140.strict$gene))) %>%  arrange(avg_log2FC)
Sp140.strict$genotype <- "Sp140"
Irg1.strict <- rbind(Irg1.strict, filter(all.heatmap.genes, gene %in% setdiff(all.heatmap.genes$gene, Irg1.strict$gene))) %>%  arrange(avg_log2FC)
Irg1.strict$genotype <- "Irg1"
iNOS.strict <- rbind(iNOS.strict, filter(all.heatmap.genes, gene %in% setdiff(all.heatmap.genes$gene, iNOS.strict$gene))) %>%  arrange(avg_log2FC)
iNOS.strict$genotype <- "iNOS"
gene.sig.combined <- rbind(B6.strict,Sp140.strict,Irg1.strict,iNOS.strict) %>% arrange(gene)

ggplot(gene.sig.combined, aes(genotype, gene, fill= avg_log2FC)) + 
  geom_tile() +
  scale_fill_gradient2(low = "black")


#ssgsea makes the point well but singscore is ever better!
mouse.tb <- hack_sig(tb.counts.norm, signatures = c(list(B6.strict$gene),list(Sp140.strict$gene),list(Irg1.strict$gene),list(iNOS.strict$gene)), method = "singscore")
mouse.tb$status <- ifelse(grepl(paste(TB.infected, collapse = "|"), mouse.tb$sample_id), "TB","Uninfected")
b6.roc <- roc(mouse.tb$status, mouse.tb$sig1, levels = c("TB","Uninfected"))
sp140.roc <- roc(mouse.tb$status, mouse.tb$sig2, levels = c("TB","Uninfected"))
irg1.roc <- roc(mouse.tb$status, mouse.tb$sig3, levels = c("TB","Uninfected"))
inos.roc <- roc(mouse.tb$status, mouse.tb$sig4, levels = c("TB","Uninfected"))
plot(b6.roc, print.auc = TRUE, col = "red", print.auc.y =0.47)
plot(sp140.roc, add=TRUE, print.auc = TRUE, col = "blue",  print.auc.y =0.4)
plot(irg1.roc,add = TRUE, print.auc = TRUE, col = "chartreuse", print.auc.y =0.33)
plot(inos.roc, add=TRUE, print.auc = TRUE, col = "cadetblue2",  print.auc.y =0.27)

#Use the Berry South Africa dataset for the comparison
berry.southafrica <- getGEO(GEO= "GSE107992")
berry.southafrica <- pData(berry.southafrica[[1]])
row.names(berry.southafrica) <- berry.southafrica$title
berry.southafrica <- berry.southafrica[2:43]
colnames(berry.southafrica)[colnames(berry.southafrica) == 'group:ch1'] <- "disease.state"
berry.southafrica.data <- read_excel("/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Stromal Cell IL-1 scRNA-seq 7-22-21/GSE107992_Raw_counts_Berry_SouthAfrica.xlsx")
berry.southafrica.data <- as.data.frame(berry.southafrica.data)
row.names(berry.southafrica.data) <- berry.southafrica.data$Genes
berry.southafrica.data <- berry.southafrica.data[4:50]

#Generate dds object
berry.southafrica.dds <- DESeqDataSetFromMatrix(countData = berry.southafrica.data,
                                           colData = berry.southafrica,
                                           design = ~ disease.state)
keep <- rowSums(counts(berry.southafrica.dds)) >= 10
berry.southafrica.dds <- berry.southafrica.dds[keep,]
berry.southafrica.dds$disease.state <- relevel(berry.southafrica.dds$disease.state, ref = "LTBI")

# run DESeq
berry.southafrica.dds <- DESeq(berry.southafrica.dds)

#Normalize and plot PCA colored by group
berry.southafrica.dds.vst <- vst(berry.southafrica.dds)
plotPCA(berry.southafrica.dds.vst,
        intgroup = c('disease.state'),
        returnData = FALSE)

#Normalized Counts
berry.southafrica.norm <- counts(berry.southafrica.dds, normalized = T)
berry.southafrica.norm <- as.data.frame(berry.southafrica.norm)

#convert ENSG gene names to gene symbol
berry.southafrica.norm$gene <- mapIds(org.Hs.eg.db, row.names(berry.southafrica.norm), keytype="ENSEMBL", column="SYMBOL", multiVals = "first")
berry.southafrica.norm <- berry.southafrica.norm[!is.na(berry.southafrica.norm$gene),]
berry.southafrica.norm$gene <- make.unique(berry.southafrica.norm$gene, sep="-")
row.names(berry.southafrica.norm) <- berry.southafrica.norm$gene

#Get IDs for the two classes to compare
berry.southafrica.tb <- row.names(berry.southafrica[berry.southafrica$disease.state == "Active_TB",])
berry.southafrica.ltbi <- row.names(berry.southafrica[berry.southafrica$disease.state == "LTBI",])

saveRDS(berry.southafrica.tb, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Dmitri Data/Coding stuff/Human TB Gene Signatures/berry.southafrica.tb")

#Perform ROC analysis on Berry South Africa
berry.southafrica.norm <- berry.southafrica.norm[,1:47]

saveRDS(berry.southafrica.norm, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Dmitri Data/Coding stuff/Human TB Gene Signatures/berry.southafrica.norm")

berry.southafrica.sig <- hack_sig(berry.southafrica.norm, signatures = c(list(B6.strict$gene),list(Sp140.strict$gene),list(Irg1.strict$gene),list(iNOS.strict$gene)), method = "singscore")
berry.southafrica.sig$status <- ifelse(grepl(paste(berry.southafrica.tb, collapse = "|"), berry.southafrica.sig$sample_id), "TB", "LTBI")
b6.berry.southafrica.roc <- roc(berry.southafrica.sig$status, berry.southafrica.sig$sig1, levels = c("TB","LTBI"))
sp140.berry.southafrica.tb.roc <- roc(berry.southafrica.sig$status, berry.southafrica.sig$sig2, levels = c("TB","LTBI"))
irg1.berry.southafrica.tb.roc <- roc(berry.southafrica.sig$status, berry.southafrica.sig$sig3, levels = c("TB","LTBI"))
inos.berry.southafrica.tb.roc <- roc(berry.southafrica.sig$status, berry.southafrica.sig$sig4, levels = c("TB","LTBI"))
plot(b6.berry.southafrica.roc, print.auc = TRUE, col = "red", print.auc.y =0.47)
plot(sp140.berry.southafrica.tb.roc, add=TRUE, print.auc = TRUE, col = "blue",  print.auc.y =0.4)
plot(irg1.berry.southafrica.tb.roc,add = TRUE, print.auc = TRUE, col = "chartreuse", print.auc.y =0.33)
plot(inos.berry.southafrica.tb.roc, add=TRUE, print.auc = TRUE, col = "cadetblue2",  print.auc.y =0.27)

library(sjmisc)
counts <- rotate_df(berry.southafrica.norm)
counts <- counts[1:47,]
counts$sample <- row.names(counts) #23687
counts$Class <- ifelse(grepl(paste(berry.southafrica.tb, collapse = "|"), berry.southafrica.sig$sample_id), "TB", "LTBI")

IL1RN.score.berrySA <- counts %>% dplyr::select(c("sample","Class","IL1RN"))
IL1RN.score.berrySA$status <- ifelse(grepl(paste(berry.southafrica.tb, collapse = "|"), berry.southafrica.sig$sample_id), "TB", "LTBI")
IL1RN.score.berrySA$IL1RN <- as.numeric(IL1RN.score.berrySA$IL1RN)
berrySA.IL1RN.roc <- roc(IL1RN.score.berrySA$status, IL1RN.score.berrySA$IL1RN)
plot(berrySA.IL1RN.roc, print.auc = TRUE, col = "red", print.auc.y =0.47)

#Look at data comparing Lung TB, adenocarcinoma, and sarcoidosis.
lung.inflammation<- getGEO(GEO= "GSE148036")
lung.inflammation <- pData(lung.inflammation[[1]])
lung.inflammation <- lung.inflammation %>% mutate(across(c('characteristics_ch1.3'), substr, 9, nchar(characteristics_ch1.3)))
row.names(lung.inflammation) <- lung.inflammation$title
colnames(lung.inflammation)[colnames(lung.inflammation) == 'characteristics_ch1.3'] <- "disease.state"
lung.inflammation.data <- read.delim("/Users/dmitrikotov/Downloads/GSE148036_PRJNA609278count_matrix.txt")
row.names(lung.inflammation.data) <- lung.inflammation.data$Geneid
lung.inflammation.data <- lung.inflammation.data[2:21]
lung.inflammation.samples <- row.names(lung.inflammation)
colnames(lung.inflammation.data) <- lung.inflammation.samples

#Generate dds object
lung.inflammation.dds <- DESeqDataSetFromMatrix(countData = lung.inflammation.data,
                                        colData = lung.inflammation,
                                        design = ~ disease.state)
keep <- rowSums(counts(lung.inflammation.dds)) >= 10
lung.inflammation.dds <- lung.inflammation.dds[keep,]
lung.inflammation.dds$disease.state <- relevel(lung.inflammation.dds$disease.state, ref = " Normal")

# run DESeq
lung.inflammation.dds <- DESeq(lung.inflammation.dds)

#Normalize and plot PCA colored by group
lung.inflammation.dds.vst <- vst(lung.inflammation.dds)
plotPCA(lung.inflammation.dds.vst,
        intgroup = c('disease.state'),
        returnData = FALSE)

#Normalized Counts
lung.inflammation.norm <- counts(lung.inflammation.dds, normalized = T)
lung.inflammation.norm <- as.data.frame(lung.inflammation.norm)

#convert ENSG gene names to gene symbol
lung.inflammation.norm$gene <- mapIds(org.Hs.eg.db, row.names(lung.inflammation.norm), keytype="ENSEMBL", column="SYMBOL", multiVals = "first")
lung.inflammation.norm <- lung.inflammation.norm[!is.na(lung.inflammation.norm$gene),]
lung.inflammation.norm$gene <- make.unique(lung.inflammation.norm$gene, sep="-")
row.names(lung.inflammation.norm) <- lung.inflammation.norm$gene
lung.inflammation.norm <- lung.inflammation.norm[,1:20]

#Get IDs for the three classes to compare
lung.inflammation.tb <- row.names(lung.inflammation[lung.inflammation$disease.state == " Tuberculosis",])
lung.inflammation.adeno <- row.names(lung.inflammation[lung.inflammation$disease.state == " Adenocarcinoma",])
lung.inflammation.sacro <- row.names(lung.inflammation[lung.inflammation$disease.state == " Sacrodosis",])

#Perform ROC analysis on Lung inflammation dataset 1
lung.inflammation.sig <- hack_sig(lung.inflammation.norm, signatures = c(list(B6.strict$gene),list(Sp140.strict$gene),list(Irg1.strict$gene),list(iNOS.strict$gene)), method = "singscore")
lung.inflammation.sig$status <- ifelse(grepl(paste(lung.inflammation.tb, collapse = "|"), lung.inflammation.sig$sample_id), "TB", ifelse(grepl(paste(lung.inflammation.adeno, collapse = "|"), lung.inflammation.sig$sample_id), "Adeno",ifelse(grepl(paste(lung.inflammation.sacro, collapse = "|"), lung.inflammation.sig$sample_id), "Sacro","Normal")))
lung.inflammation.sig.tb.adeno <- lung.inflammation.sig %>% filter(status == "TB" | status == "Adeno")
b6.lung.inflammation.sig.tb.adeno <- roc(lung.inflammation.sig.tb.adeno$status, lung.inflammation.sig.tb.adeno$sig1, levels = c("TB","Adeno"))
sp140.lung.inflammation.sig.tb.adeno.roc <- roc(lung.inflammation.sig.tb.adeno$status, lung.inflammation.sig.tb.adeno$sig2, levels = c("TB","Adeno"))
irg1.lung.inflammation.sig.tb.adeno.roc <- roc(lung.inflammation.sig.tb.adeno$status, lung.inflammation.sig.tb.adeno$sig3, levels = c("TB","Adeno"))
inos.lung.inflammation.sig.tb.adeno.roc <- roc(lung.inflammation.sig.tb.adeno$status, lung.inflammation.sig.tb.adeno$sig4, levels = c("TB","Adeno"))
plot(b6.lung.inflammation.sig.tb.adeno, print.auc = TRUE, col = "red", print.auc.y =0.47)
plot(sp140.lung.inflammation.sig.tb.adeno.roc, add=TRUE, print.auc = TRUE, col = "blue",  print.auc.y =0.4)
plot(irg1.lung.inflammation.sig.tb.adeno.roc,add = TRUE, print.auc = TRUE, col = "chartreuse", print.auc.y =0.33)
plot(inos.lung.inflammation.sig.tb.adeno.roc, add=TRUE, print.auc = TRUE, col = "cadetblue2",  print.auc.y =0.27)

lung.inflammation.sig.tb.sacro <- lung.inflammation.sig %>% filter(status == "TB" | status == "Sacro")
b6.lung.inflammation.sig.tb.sacro <- roc(lung.inflammation.sig.tb.sacro$status, lung.inflammation.sig.tb.sacro$sig1, levels = c("TB","Sacro"))
sp140.lung.inflammation.sig.tb.sacro.roc <- roc(lung.inflammation.sig.tb.sacro$status, lung.inflammation.sig.tb.sacro$sig2, levels = c("TB","Sacro"))
irg1.lung.inflammation.sig.tb.sacro.roc <- roc(lung.inflammation.sig.tb.sacro$status, lung.inflammation.sig.tb.sacro$sig3, levels = c("TB","Sacro"))
inos.lung.inflammation.sig.tb.sacro.roc <- roc(lung.inflammation.sig.tb.sacro$status, lung.inflammation.sig.tb.sacro$sig4, levels = c("TB","Sacro"))
plot(b6.lung.inflammation.sig.tb.sacro, print.auc = TRUE, col = "red", print.auc.y =0.47)
plot(sp140.lung.inflammation.sig.tb.sacro.roc, add=TRUE, print.auc = TRUE, col = "blue",  print.auc.y =0.4)
plot(irg1.lung.inflammation.sig.tb.sacro.roc,add = TRUE, print.auc = TRUE, col = "chartreuse", print.auc.y =0.33)
plot(inos.lung.inflammation.sig.tb.sacro.roc, add=TRUE, print.auc = TRUE, col = "cadetblue2",  print.auc.y =0.27)

#Look at data comparing lymph node TB and sarcoidosis #2.
ln.inflammation <- getGEO(GEO= "GSE157671")
ln.inflammation <- pData(ln.inflammation[[1]])
row.names(ln.inflammation) <- ln.inflammation$geo_accession
colnames(ln.inflammation)[colnames(ln.inflammation) == 'sample type:ch1'] <- "disease.state"
ln.inflammation$tissue.disease <- paste(ln.inflammation$source_name_ch1, ln.inflammation$disease.state, sep = "_")
ln.inflammation.data.2 <- read.table(file = "/Users/dmitrikotov/Downloads/GSE157671_raw_counts_GRCh38.p13_NCBI.tsv", header = T)
rownames(ln.inflammation.data.2) <- ln.inflammation.data.2$GeneID
ln.inflammation.data.2 <- ln.inflammation.data.2[2:32]

#Generate dds object
ln.inflammation.dds <- DESeqDataSetFromMatrix(countData = ln.inflammation.data.2,
                                                colData = ln.inflammation,
                                                design = ~ tissue.disease)
keep <- rowSums(counts(ln.inflammation.dds)) >= 10
ln.inflammation.dds <- ln.inflammation.dds[keep,]
ln.inflammation.dds$tissue.disease <- relevel(ln.inflammation.dds$tissue.disease, ref = "lymph nodes_healthy tissue")

# run DESeq
ln.inflammation.dds <- DESeq(ln.inflammation.dds)

#Normalize and plot PCA colored by group
ln.inflammation.dds.vst <- vst(ln.inflammation.dds)
plotPCA(ln.inflammation.dds.vst,
        intgroup = c('tissue.disease'),
        returnData = FALSE)

#Normalized Counts
ln.inflammation.norm <- counts(ln.inflammation.dds, normalized = T)
ln.inflammation.norm <- as.data.frame(ln.inflammation.norm)

#convert ENSG gene names to gene symbol
ln.inflammation.norm$gene <- mapIds(org.Hs.eg.db, row.names(ln.inflammation.norm), keytype="ENTREZID", column="SYMBOL", multiVals = "first")
ln.inflammation.norm <- ln.inflammation.norm[!is.na(ln.inflammation.norm$gene),]
ln.inflammation.norm$gene <- make.unique(ln.inflammation.norm$gene, sep="-")
row.names(ln.inflammation.norm) <- ln.inflammation.norm$gene
ln.inflammation.norm <- ln.inflammation.norm[,1:31]

#Get IDs for the two classes to compare
ln.inflammation.tb <- row.names(ln.inflammation[ln.inflammation$tissue.disease == "lymph nodes_tuberculosis granuloma",])
ln.inflammation.sacro <- row.names(ln.inflammation[ln.inflammation$tissue.disease == "lymph nodes_sarcoidosis granuloma",])

#Perform ROC analysis on Lung inflammation dataset 1
ln.inflammation.sig <- hack_sig(ln.inflammation.norm, signatures = c(list(B6.strict$gene),list(Sp140.strict$gene),list(Irg1.strict$gene),list(iNOS.strict$gene)), method = "singscore")
ln.inflammation.sig$status <- ifelse(grepl(paste(ln.inflammation.tb, collapse = "|"), ln.inflammation.sig$sample_id), "TB", ifelse(grepl(paste(ln.inflammation.sacro, collapse = "|"), ln.inflammation.sig$sample_id), "Sarco","healthy"))
ln.inflammation.sig.tb <- ln.inflammation.sig %>% filter(status == "TB" | status == "Sarco")
b6.ln.inflammation.sig.tb <- roc(ln.inflammation.sig.tb$status, ln.inflammation.sig.tb$sig1, levels = c("TB","Sarco"))
sp140.ln.inflammation.sig.tb <- roc(ln.inflammation.sig.tb$status, ln.inflammation.sig.tb$sig2, levels = c("TB","Sarco"))
irg1.ln.inflammation.sig.tb <- roc(ln.inflammation.sig.tb$status, ln.inflammation.sig.tb$sig3, levels = c("TB","Sarco"))
inos.ln.inflammation.sig.tb <- roc(ln.inflammation.sig.tb$status, ln.inflammation.sig.tb$sig4, levels = c("TB","Sarco"))
plot(b6.ln.inflammation.sig.tb, print.auc = TRUE, col = "red", print.auc.y =0.47)
plot(sp140.ln.inflammation.sig.tb, add=TRUE, print.auc = TRUE, col = "blue",  print.auc.y =0.4)
plot(irg1.ln.inflammation.sig.tb,add = TRUE, print.auc = TRUE, col = "chartreuse", print.auc.y =0.33)
plot(inos.ln.inflammation.sig.tb, add=TRUE, print.auc = TRUE, col = "cadetblue2",  print.auc.y =0.27)

#Do QC on mouse signatures using the TANDEM cohort.
mRNA_expr_matrix = list()
mRNA_expr_matrix[["TANDEM"]] = tb.counts.norm
mRNA_expr_matrix[["Berry SA"]] = berry.southafrica.norm

gene_sigs_list = list()
gene_sigs_list[['B6']] = as.matrix(B6.strict$gene)
gene_sigs_list[['Nos2']] = as.matrix(iNOS.strict$gene)
gene_sigs_list[['Acod1']] = as.matrix(Irg1.strict$gene)
gene_sigs_list[['Sp140']] = as.matrix(Sp140.strict$gene)

names_sigs = c("B6","Nos2","Acod1","Sp140")
names_datasets = c("TANDEM","Berry SA")

showResults <- FALSE # we do not want to show the reuslts in R graphics windows
doNegativeControl <- FALSE # we do not want to compute the negative or permutation controls for time purposes

make_all_plots(gene_sigs_list, mRNA_expr_matrix, showResults = showResults, names_sigs = names_sigs, names_datasets = names_datasets , doNegativeControl = doNegativeControl, out_dir = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Human TB Gene Signatures/Gene signature QC")

