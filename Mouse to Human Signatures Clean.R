library(Seurat)
library(tidyverse)
library(patchwork)
library(cowplot)
library(nichenetr)
library(GEOquery)
library(umap)
library(DESeq2)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(hacksig)
library(pROC)

#Define signatures for B6 infected versus naive
super.integrated <- readRDS(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/b6 sp140 irg1 inos v2")
super.integrated$infected <- ifelse(super.integrated$state == "Mtb+" | super.integrated$state == "Mtb-", "infected","naive")
super.integrated$cell_type.infected <- paste(super.integrated$cell_type, super.integrated$infected, sep = "_")
Idents(super.integrated) <- "cell_type.infected"
B6 <- FindMarkers(super.integrated, ident.1 = "B6_infected", ident.2 = "B6_naive")
Sp140 <- FindMarkers(super.integrated, ident.1 = "Sp140_infected", ident.2 = "Sp140_naive")
Irg1 <- FindMarkers(super.integrated, ident.1 = "Irg1_infected", ident.2 = "Irg1_naive")
iNOS <- FindMarkers(super.integrated, ident.1 = "iNOS_infected", ident.2 = "iNOS_naive")
B6$gene <- rownames(B6)
Sp140$gene <- rownames(Sp140)
Irg1$gene <- rownames(Irg1)
iNOS$gene <- rownames(iNOS)

#Filtering Steps
B6.strict <- B6 %>% filter(avg_log2FC > 3) %>% mutate(gene = convert_mouse_to_human_symbols(gene)) %>% drop_na()
Sp140.strict <- Sp140 %>% filter(avg_log2FC > 3) %>% mutate(gene = convert_mouse_to_human_symbols(gene)) %>% drop_na()
Irg1.strict <- Irg1 %>% filter(avg_log2FC > 3) %>% mutate(gene = convert_mouse_to_human_symbols(gene)) %>% drop_na()
iNOS.strict <- iNOS %>% filter(avg_log2FC > 3) %>% mutate(gene = convert_mouse_to_human_symbols(gene)) %>% drop_na()

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
        intgroup = c('TB.status'),
        returnData = FALSE)

#Normalized Counts
tb.counts.norm <- counts(dds, normalized = T)
tb.counts.norm <- as.data.frame(tb.counts.norm)

#convert ENSG gene names to gene symbol
tb.counts.norm$gene <- mapIds(org.Hs.eg.db, row.names(tb.counts.norm), keytype="ENSEMBL", column="SYMBOL", multiVals = "first")
tb.counts.norm <- tb.counts.norm[!is.na(tb.counts.norm$gene),]
tb.counts.norm$gene <- make.unique(tb.counts.norm$gene, sep="-")
row.names(tb.counts.norm) <- tb.counts.norm$gene

#Get IDs for the two classes to compare
TB.infected <- row.names(human.tb.db.meta[human.tb.db.meta$TB.status == "TB",])
Uninfected <- row.names(human.tb.db.meta[human.tb.db.meta$TB.status == "Uninfected",])

saveRDS(TB.infected, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Dmitri Data/Coding stuff/Human TB Gene Signatures/TB")

#gene symbol = row names and sample IDs = column names
tb.counts.norm <- tb.counts.norm[,1:249]
saveRDS(tb.counts.norm, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Dmitri Data/Coding stuff/Human TB Gene Signatures/tb.counts.norm")

#Score Immunosuppressive genes
immunosuppresion <- list(c("IL1RN","IL18BP","CD274","TGFB1","IDO1","IL10","ARG1"))
check_sig(tb.counts.norm, immunosuppresion)
immunosuppresion.scored <- hack_sig(tb.counts.norm, signatures = immunosuppresion, method = "ssgsea")
immunosuppresion.scored$status <- ifelse(grepl(paste(TB.infected, collapse = "|"), immunosuppresion.scored$sample_id), "TB","Uninfected")
immunosuppresion.roc <- roc(immunosuppresion.scored$status, immunosuppresion.scored$sig1, levels = c("TB","Uninfected")) #IFNg zscore AUC 83.8 ssgsea AUC 91.6 singscore AUC 88.7
plot(immunosuppresion.roc, print.auc = TRUE, col = "red", print.auc.y =0.47)

#ssgsea makes the point well but singscore is ever better!
mouse.tb <- hack_sig(tb.counts.norm, signatures = c(list(unique(B6.strict$gene)),list(unique(Sp140.strict$gene)),list(unique(Irg1.strict$gene)),list(unique(iNOS.strict$gene))) , method = "singscore")
mouse.tb$status <- ifelse(grepl(paste(TB.infected, collapse = "|"), mouse.tb$sample_id), "TB","Uninfected")
b6.roc <- roc(mouse.tb$status, mouse.tb$sig1, levels = c("TB","Uninfected"))
sp140.roc <- roc(mouse.tb$status, mouse.tb$sig2, levels = c("TB","Uninfected"))
irg1.roc <- roc(mouse.tb$status, mouse.tb$sig3, levels = c("TB","Uninfected"))
inos.roc <- roc(mouse.tb$status, mouse.tb$sig4, levels = c("TB","Uninfected"))
plot(b6.roc, print.auc = TRUE, col = "red", print.auc.y =0.47)
plot(sp140.roc, add=TRUE, print.auc = TRUE, col = "blue",  print.auc.y =0.4)
plot(irg1.roc,add = TRUE, print.auc = TRUE, col = "chartreuse", print.auc.y =0.33)
plot(inos.roc, add=TRUE, print.auc = TRUE, col = "cadetblue2",  print.auc.y =0.27)

write.csv(B6.strict, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Human TB Gene Signatures/Mouse Gene Signatures/B6.csv")
write.csv(Sp140.strict, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Human TB Gene Signatures/Mouse Gene Signatures/Sp140.csv")
write.csv(Irg1.strict, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Human TB Gene Signatures/Mouse Gene Signatures/Irg1.csv")
write.csv(iNOS.strict, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Human TB Gene Signatures/Mouse Gene Signatures/iNOS.csv")