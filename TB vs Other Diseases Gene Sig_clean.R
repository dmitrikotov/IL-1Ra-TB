library(tidyverse)
library(hacksig)
library(pROC)
library(readxl)

#Read in mouse signatures
B6.strict <- read.csv(file = '/Users/dmitrikotov/Library/CloudStorage/Box-Box/DK Postdoc Data/Coding stuff/Human TB Gene Signatures/Mouse Gene Signatures/B6.csv')
Sp140.strict <- read.csv(file = '/Users/dmitrikotov/Library/CloudStorage/Box-Box/DK Postdoc Data/Coding stuff/Human TB Gene Signatures/Mouse Gene Signatures/Sp140.csv')
iNOS.strict <- read.csv(file = '/Users/dmitrikotov/Library/CloudStorage/Box-Box/DK Postdoc Data/Coding stuff/Human TB Gene Signatures/Mouse Gene Signatures/iNOS.csv')
Irg1.strict <- read.csv(file = '/Users/dmitrikotov/Library/CloudStorage/Box-Box/DK Postdoc Data/Coding stuff/Human TB Gene Signatures/Mouse Gene Signatures/Irg1.csv')

#Analysis of GSE42834 starting with the processed data available from GEO with per gene signal averaged when there are multiple probes for a given gene - 
GSE42834 <- read_xlsx("/Users/dmitrikotov/Library/CloudStorage/Box-Box/Dmitri_Kotov/GSE42834_cleanedmatrix.xlsx",col_names = F)
#Filter down to the whole sample dataset - this strips out the cell separation dataset
counts.GSE42834 <- GSE42834[c(3,4,111:391)]
counts.GSE42834 <- counts.GSE42834[c(10,12:nrow(counts.GSE42834)),]
counts.GSE42834[1,1] <- "Symbol"
counts.GSE42834[1,2] <- "Entrez_Gene_ID"
counts.GSE42834 <- counts.GSE42834 %>% na.omit %>% as.data.frame()
colnames(counts.GSE42834) <- counts.GSE42834[1,]
counts.GSE42834 <- counts.GSE42834[2:nrow(counts.GSE42834),]

meta.42834 <- GSE42834[c(1,4,6,7,10),]
meta.42834 <- meta.42834[111:391]
meta.42834 <- t(meta.42834)
meta.42834 <- data.frame(meta.42834[,1:4], row.names = meta.42834[,5])
colnames(meta.42834) <- c("disease_state","Tissue","diagnosis","Tissue_2")


#gene symbol = row names and sample IDs = column names
counts.GSE42834$Symbol <- make.unique(counts.GSE42834$Symbol, sep="-")
row.names(counts.GSE42834) <- counts.GSE42834$Symbol
counts.GSE42834 <- counts.GSE42834[3:283]

human.disease.42834 <- hack_sig(counts.GSE42834, signatures = c(list(unique(B6.strict$gene)),list(unique(Sp140.strict$gene)),list(unique(Irg1.strict$gene)),list(unique(iNOS.strict$gene))) , method = "singscore")
human.disease.42834$status <- ifelse(grepl(paste(row.names(filter(meta.42834, disease_state == "TB")), collapse = "|"), human.disease.42834$sample_id), "TB",
                                      ifelse(grepl(paste(row.names(filter(meta.42834, disease_state == "Cancer")), collapse = "|"), human.disease.42834$sample_id), "Cancer",
                                             ifelse(grepl(paste(row.names(filter(meta.42834, disease_state == "Sarcoid")), collapse = "|"), human.disease.42834$sample_id), "Sarcoid",
                                                    ifelse(grepl(paste(row.names(filter(meta.42834, disease_state == "Pneumonia")), collapse = "|"), human.disease.42834$sample_id), "Pneumonia",
                                                           ifelse(grepl(paste(row.names(filter(meta.42834, disease_state == "Pneumonia_Pre-treatment")), collapse = "|"), human.disease.42834$sample_id), "Pneumonia",
                                                                  ifelse(grepl(paste(row.names(filter(meta.42834, disease_state == "Control")), collapse = "|"), human.disease.42834$sample_id), "Control","Other"))))))

b6.roc <- roc(human.disease.42834$status, human.disease.42834$sig1, levels = c("TB","Pneumonia"))
sp140.roc <- roc(human.disease.42834$status, human.disease.42834$sig2, levels = c("TB","Pneumonia"))
irg1.roc <- roc(human.disease.42834$status, human.disease.42834$sig3, levels = c("TB","Pneumonia"))
inos.roc <- roc(human.disease.42834$status, human.disease.42834$sig4, levels = c("TB","Pneumonia"))
plot(b6.roc, print.auc = TRUE, col = "red", print.auc.y =0.47)
plot(sp140.roc, add=TRUE, print.auc = TRUE, col = "blue",  print.auc.y =0.4)
plot(irg1.roc,add = TRUE, print.auc = TRUE, col = "chartreuse", print.auc.y =0.33)
plot(inos.roc, add=TRUE, print.auc = TRUE, col = "cadetblue2",  print.auc.y =0.27)

b6.roc <- roc(human.disease.42834$status, human.disease.42834$sig1, levels = c("TB","Cancer"))
sp140.roc <- roc(human.disease.42834$status, human.disease.42834$sig2, levels = c("TB","Cancer"))
irg1.roc <- roc(human.disease.42834$status, human.disease.42834$sig3, levels = c("TB","Cancer"))
inos.roc <- roc(human.disease.42834$status, human.disease.42834$sig4, levels = c("TB","Cancer"))
plot(b6.roc, print.auc = TRUE, col = "red", print.auc.y =0.47)
plot(sp140.roc, add=TRUE, print.auc = TRUE, col = "blue",  print.auc.y =0.4)
plot(irg1.roc,add = TRUE, print.auc = TRUE, col = "chartreuse", print.auc.y =0.33)
plot(inos.roc, add=TRUE, print.auc = TRUE, col = "cadetblue2",  print.auc.y =0.27)

#Save the files
saveRDS(counts.GSE42834, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Dmitri Personal/DK Postdoc Data and Analysis/Coding stuff/TB vs Other Diseases Gene Sig/counts_GSE42834")
saveRDS(meta.42834, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Dmitri Personal/DK Postdoc Data and Analysis/Coding stuff/TB vs Other Diseases Gene Sig/meta_GSE42834")
