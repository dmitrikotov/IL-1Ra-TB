library(GEOquery)
library(tidyverse)
library(umap)
library(DESeq2)

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

#change Ensembl ids to gene symbols
library(org.Hs.eg.db)
library(AnnotationDbi)

#convert ENSG gene names to gene symbol
tb.counts.norm$gene <- mapIds(org.Hs.eg.db, row.names(tb.counts.norm), keytype="ENSEMBL", column="SYMBOL", multiVals = "first")
tb.counts.norm <- tb.counts.norm[!is.na(tb.counts.norm$gene),]
tb.counts.norm$gene <- make.unique(tb.counts.norm$gene, sep="-")
row.names(tb.counts.norm) <- tb.counts.norm$gene

#Get IDs for the two classes to compare
TB.infected <- row.names(human.tb.db.meta[human.tb.db.meta$TB.status == "TB",])
Uninfected <- row.names(human.tb.db.meta[human.tb.db.meta$TB.status == "Uninfected",])

#Flip dataframe with normalized counts
library(sjmisc)
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

library(pROC)
library(hacksig)
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
B6.mtb.pos.strict <- read.csv(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/B6 Mtb Pos Strict.csv") %>% filter(change == "up")
B6.bystander.strict <- read.csv(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/B6 Bystander Strict.csv") %>% filter(change == "up")
Sp140.mtb.pos.strict <- read.csv(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/Sp140 Mtb Pos Strict.csv") %>% filter(change == "up")
Irg1.mtb.pos.strict <- read.csv(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/Irg1 Mtb Pos Strict.csv") %>% filter(change == "up")
iNOS.mtb.pos.strict <- read.csv(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/iNOS Mtb Pos Strict.csv") %>% filter(change == "up")

#Getting rid of duplicate human genes
Sp140.strict <- filter(Sp140.strict, !(X %in% c("H2-Q6","Ifit3b")))
Irg1.strict <- filter(Irg1.strict, !(X %in% c("Gbp4","Gbp8")))
iNOS.strict <- filter(iNOS.strict, !(X %in% c("Gbp4","Gbp8")))

#HeatMap of Colors
B6.strict <- B6.strict %>% filter(avg_log2FC > 3)  %>%  arrange(avg_log2FC) %>% select(avg_log2FC, gene)
Sp140.strict <- Sp140.strict %>% filter(avg_log2FC > 3)  %>%  arrange(avg_log2FC) %>% select(avg_log2FC, gene)
Irg1.strict <- Irg1.strict %>% filter(avg_log2FC > 3)  %>%  arrange(avg_log2FC) %>% select(avg_log2FC, gene)
iNOS.strict <- iNOS.strict %>% filter(avg_log2FC > 3)  %>%  arrange(avg_log2FC) %>% select(avg_log2FC, gene)

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

#Only looking at Mtb positive cells
mouse.mtb.pos <- hack_sig(tb.counts.norm, signatures = c(list(B6.mtb.pos.strict$gene),list(Sp140.mtb.pos.strict$gene),list(Irg1.mtb.pos.strict$gene),list(iNOS.mtb.pos.strict$gene)), method = "singscore")
mouse.mtb.pos$status <- ifelse(grepl(paste(TB.infected, collapse = "|"), mouse.mtb.pos$sample_id), "TB","Uninfected")
b6.roc <- roc(mouse.tb$status, mouse.tb$sig1, levels = c("TB","Uninfected"))
sp140.roc <- roc(mouse.tb$status, mouse.tb$sig2, levels = c("TB","Uninfected"))
irg1.roc <- roc(mouse.tb$status, mouse.tb$sig3, levels = c("TB","Uninfected"))
inos.roc <- roc(mouse.tb$status, mouse.tb$sig4, levels = c("TB","Uninfected"))
plot(b6.roc, print.auc = TRUE, col = "red", print.auc.y =0.47)
plot(sp140.roc, add=TRUE, print.auc = TRUE, col = "blue",  print.auc.y =0.4)
plot(irg1.roc,add = TRUE, print.auc = TRUE, col = "chartreuse", print.auc.y =0.33)
plot(inos.roc, add=TRUE, print.auc = TRUE, col = "cadetblue2",  print.auc.y =0.27)


#TB infected vs bystander cells for B6
B6.Mtb.pos <- read.csv(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/Mtb Pos Genes.csv")
B6.Mtb.neg <- read.csv(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/Bystander Genes.csv")
B6.mtb.pos.strict
B6.bystander.strict
B6.Mtb.pos <- B6.Mtb.pos %>% filter(change == "up") %>% filter(avg_log2FC >= 2)
B6.Mtb.neg <-  B6.Mtb.neg %>% filter(change == "up")  %>% filter(avg_log2FC >= 2)
TB.infection.status <- hack_sig(tb.counts.norm, signatures = c(list(B6.Mtb.pos$gene),list(B6.Mtb.neg$gene)), method = "ssgsea")
TB.infection.status$status <- ifelse(grepl(paste(TB.infected, collapse = "|"), TB.infection.status$sample_id), "TB","Uninfected")
B6.Mtb.pos.roc <- roc(TB.infection.status$status, TB.infection.status$sig1, levels = c("TB","Uninfected"))
B6.Mtb.neg.roc <- roc(TB.infection.status$status, TB.infection.status$sig2, levels = c("TB","Uninfected"))
plot(B6.Mtb.pos.roc, print.auc = TRUE, col = "red", print.auc.y =0.47)
plot(B6.Mtb.neg.roc, add=TRUE, print.auc = TRUE, col = "blue",  print.auc.y =0.4)