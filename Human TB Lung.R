library(Seurat)
library(tidyverse)
library(patchwork)
library(nichenetr)
library(EnhancedVolcano)
library(UCell)

# Load the human cell datasets from PMID 37690670
SP019H <- Read10X_h5("/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Human TB resected lung/GSE192483_RAW/GSM5747736_SP019H.filtered_feature_bc_matrix.h5")
SP019L <- Read10X_h5("/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Human TB resected lung/GSE192483_RAW/GSM5747737_SP019L.filtered_feature_bc_matrix.h5")
SP020L <- Read10X_h5("/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Human TB resected lung/GSE192483_RAW/GSM5747738_SP020L.filtered_feature_bc_matrix.h5")
SP021H <- Read10X_h5("/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Human TB resected lung/GSE192483_RAW/GSM5747739_SP021H.filtered_feature_bc_matrix.h5")
SP021L <- Read10X_h5("/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Human TB resected lung/GSE192483_RAW/GSM5747740_SP021L.filtered_feature_bc_matrix.h5")
SP023H <- Read10X_h5("/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Human TB resected lung/GSE192483_RAW/GSM5747741_SP023H.filtered_feature_bc_matrix.h5")
SP023L <- Read10X_h5("/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Human TB resected lung/GSE192483_RAW/GSM5747742_SP023L.filtered_feature_bc_matrix.h5")
SP024H <- Read10X_h5("/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Human TB resected lung/GSE192483_RAW/GSM5747743_SP024H.filtered_feature_bc_matrix.h5")
SP024L <- Read10X_h5("/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Human TB resected lung/GSE192483_RAW/GSM5747744_SP024L.filtered_feature_bc_matrix.h5")
SP025H <- Read10X_h5("/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Human TB resected lung/GSE192483_RAW/GSM5747745_SP025H.filtered_feature_bc_matrix.h5")
SP025L <- Read10X_h5("/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Human TB resected lung/GSE192483_RAW/GSM5747746_SP025L.filtered_feature_bc_matrix.h5")

#Create Seurat object
SP019H <- CreateSeuratObject(counts = SP019H, project = "SP019")
SP019L <- CreateSeuratObject(counts = SP019L, project = "SP019")
SP020L <- CreateSeuratObject(counts = SP020L, project = "SP020")
SP021H <- CreateSeuratObject(counts = SP021H, project = "SP021")
SP021L <- CreateSeuratObject(counts = SP021L, project = "SP021")
SP023H <- CreateSeuratObject(counts = SP023H, project = "SP023")
SP023L <- CreateSeuratObject(counts = SP023L, project = "SP023")
SP024H <- CreateSeuratObject(counts = SP024H, project = "SP024")
SP024L <- CreateSeuratObject(counts = SP024L, project = "SP024")
SP025H <- CreateSeuratObject(counts = SP025H, project = "SP025")
SP025L <- CreateSeuratObject(counts = SP025L, project = "SP025")

# store mitochondrial percentage in object meta data
SP019H <- PercentageFeatureSet(SP019H, pattern = "^MT-", col.name = "percent.mt")
SP019L <- PercentageFeatureSet(SP019L, pattern = "^MT-", col.name = "percent.mt")
SP020L <- PercentageFeatureSet(SP020L, pattern = "^MT-", col.name = "percent.mt")
SP021H <- PercentageFeatureSet(SP021H, pattern = "^MT-", col.name = "percent.mt")
SP021L <- PercentageFeatureSet(SP021L, pattern = "^MT-", col.name = "percent.mt")
SP023H <- PercentageFeatureSet(SP023H, pattern = "^MT-", col.name = "percent.mt")
SP023L <- PercentageFeatureSet(SP023L, pattern = "^MT-", col.name = "percent.mt")
SP024H <- PercentageFeatureSet(SP024H, pattern = "^MT-", col.name = "percent.mt")
SP024L <- PercentageFeatureSet(SP024L, pattern = "^MT-", col.name = "percent.mt")
SP025H <- PercentageFeatureSet(SP025H, pattern = "^MT-", col.name = "percent.mt")
SP025L <- PercentageFeatureSet(SP025L, pattern = "^MT-", col.name = "percent.mt")

#Add Metadata for patient information
SP019H[["Patient"]] <- "SP019"
SP019H[["FDG.signal"]] <- "High"
SP019L[["Patient"]] <- "SP019"
SP019L[["FDG.signal"]] <- "Low"
SP020L[["Patient"]] <- "SP020"
SP020L[["FDG.signal"]] <- "Low"
SP021H[["Patient"]] <- "SP021"
SP021H[["FDG.signal"]] <- "High"
SP021L[["Patient"]] <- "SP021"
SP021L[["FDG.signal"]] <- "Low"
SP023H[["Patient"]] <- "SP023"
SP023H[["FDG.signal"]] <- "High"
SP023L[["Patient"]] <- "SP023"
SP023L[["FDG.signal"]] <- "Low"
SP024H[["Patient"]] <- "SP024"
SP024H[["FDG.signal"]] <- "High"
SP024L[["Patient"]] <- "SP024"
SP024L[["FDG.signal"]] <- "Low"
SP025H[["Patient"]] <- "SP025"
SP025H[["FDG.signal"]] <- "High"
SP025L[["Patient"]] <- "SP025"
SP025L[["FDG.signal"]] <- "Low"

#standard filters - change nFeature_RNA to have over 750
VlnPlot(SP019H, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
SP019H <- subset(SP019H, subset = nFeature_RNA > 750 & nFeature_RNA < 3000 & percent.mt < 10)
VlnPlot(SP019L, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
SP019L <- subset(SP019L, subset = nFeature_RNA > 750 & nFeature_RNA < 3000 & percent.mt < 10)
VlnPlot(SP020L, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
SP020L <- subset(SP020L, subset = nFeature_RNA > 750 & nFeature_RNA < 3000 & percent.mt < 10)
VlnPlot(SP021H, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
SP021H <- subset(SP021H, subset = nFeature_RNA > 750 & nFeature_RNA < 3000 & percent.mt < 10)
VlnPlot(SP021L, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
SP021L <- subset(SP021L, subset = nFeature_RNA > 750 & nFeature_RNA < 3000 & percent.mt < 10)
VlnPlot(SP023H, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
SP023H <- subset(SP023H, subset = nFeature_RNA > 750 & nFeature_RNA < 3000 & percent.mt < 10)
VlnPlot(SP023L, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
SP023L <- subset(SP023L, subset = nFeature_RNA > 750 & nFeature_RNA < 3000 & percent.mt < 10)
VlnPlot(SP024H, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
SP024H <- subset(SP024H, subset = nFeature_RNA > 750 & nFeature_RNA < 3000 & percent.mt < 10)
VlnPlot(SP024L, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
SP024L <- subset(SP024L, subset = nFeature_RNA > 750 & nFeature_RNA < 3000 & percent.mt < 10)
VlnPlot(SP025H, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
SP025H <- subset(SP025H, subset = nFeature_RNA > 750 & nFeature_RNA < 3000 & percent.mt < 10)
VlnPlot(SP025L, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
SP025L <- subset(SP025L, subset = nFeature_RNA > 750 & nFeature_RNA < 3000 & percent.mt < 10)

# run sctransform
SP019H <- SCTransform(SP019H, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)
SP019L <- SCTransform(SP019L, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)
SP020L <- SCTransform(SP020L, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)
SP021H <- SCTransform(SP021H, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)
SP021L <- SCTransform(SP021L, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)
SP023H <- SCTransform(SP023H, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)
SP023L <- SCTransform(SP023L, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)
SP024H <- SCTransform(SP024H, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)
SP024L <- SCTransform(SP024L, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)
SP025H <- SCTransform(SP025H, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)
SP025L <- SCTransform(SP025L, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)

#Integration using pearson residuals
sample.list <- list(SP019H = SP019H, SP019L = SP019L, SP020L = SP020L, SP021H = SP021H, SP021L = SP021L, SP023H = SP023H, SP023L = SP023L, SP024H = SP024H, SP024L = SP024L, SP025H = SP025H, SP025L = SP025L)
features <- SelectIntegrationFeatures(object.list = sample.list, nfeatures = 3000)
sample.list <- PrepSCTIntegration(object.list = sample.list, anchor.features = features)
immune.anchors <- FindIntegrationAnchors(object.list = sample.list, normalization.method = "SCT",
                                         anchor.features = features)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")

#Analysis of integrated data
immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30, verbose = FALSE)
immune.combined.sct <- FindNeighbors(immune.combined.sct, reduction = "pca", dims = 1:30)

DefaultAssay(immune.combined.sct) <- "integrated"
immune.combined.sct <- FindClusters(immune.combined.sct, resolution = 0.3)
immune.combined.sct <- FindClusters(immune.combined.sct, resolution = 0.8)
immune.combined.sct <- FindClusters(immune.combined.sct, resolution = 1.5)
immune.combined.sct <- FindClusters(immune.combined.sct, resolution = 2)
immune.combined.sct <- FindClusters(immune.combined.sct, resolution = 3)

#Work with the SCT corrected values rather than the integrated values
DefaultAssay(immune.combined.sct) <- "SCT"
immune.combined.sct <- PrepSCTFindMarkers(immune.combined.sct)

saveRDS(immune.combined.sct, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Human TB resected lung/Human_lungv2.RDS")

human.tb.lung.markers <- FindAllMarkers(immune.combined.sct)

#V2 Resolution of 3
Idents(immune.combined.sct) <- "integrated_snn_res.3"

#Following clusters have low gene numbers: 2, 11, 16, 17, 23, 26, 35, 38 and 44 is extra low.

#Apply AM and TREM2 IM signatures from NHP to human data
mac.sig <- list()
mac.sig$nhp_AM <- c("MARCO", "MCEMP1", "CD109", "RNASE2", "SCD","PDLIM1", "SLC1A5","NGFRAP1","FHL1","PPARG")
mac.sig$nhp_TREM2 <- c("FAM26F", "RGS1", "TMEM176A", "TSPAN33", "SLC9A9","PADI2", "GUCY1A3","GRASP","CLDN1","FOLR2","MEF2C","SECTM1","TLR7","STAB1","AXL","MS4A4A","F13A1","EVL","CEBPD","ITM2C","ITGAD")
mac.sig$nhp_TREM2v2 <- c("TREM2", "C1QA", "C1QB", "C1QC", "AXL","APOE", "APOC2","APOC4","CCL18","CD101","MMP2")

immune.combined.sct <- AddModuleScore_UCell(immune.combined.sct, features = mac.sig)
#Populations high for AM are 14, 47, 56
#Populations high but mixed for AM are 16, ~25, 38, 45, 55, ~62

#Populations very high for TREM2 IM v2 are 25, 47, 56; high are 14, 16, 38, 45, 55; intermediate/mixed are 7, 44, 62
#AM - 14, 16, 47
#TREM2 IM - 55, 56, 25, 29, 62
#FOLR2 IM - 7
#27 and 43 are IM
#57 is proliferating macs
#44 is cDC2
#53 is cDC1
#70 is mregDC or DC3?
#35 is monocytes
#68 is pDCs
#69 is plasma cells
#B cells are 60 and 0
#67 is mast
#

#52 is proliferating cells - T cells?

#Tregs - 71, 46, 8

immune.combined.sct <- readRDS(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Human TB resected lung/Human_lungv2.RDS")
human.tb <- RenameIdents(immune.combined.sct, '60' = "Mast", '37' = "ILC", '47' = "Cycling T", '4' = "Treg",'45' = "CD8 T",'39' = "CD8 T",'30' = "CD4 T", '57' = "CD4 T", '15' = "CD4 T",'53' = "CD4 T",
                         '31' = "CD8 T",'20' = "CD8 T",'22' = "CD8 T",'36' = "CD8 T",'18' = "CD8 T",'9' = "CD8 T",'33'="CD4 T",'55'="CD4 T",'6'="CD4 T",'2'="CD4 T",'1'="CD4 T",'10'="CD4 T",'29'="CD4 T",
                         '42' = "CD4 T",'3' = "CD4 T",'8'="CD4 T",'11'="CD4 T",'38'="CD4 T",'28'="CD4 T",'7'="CD4 T",'5'="CD4 T",'43'="CD4 T",'62' = "plasma cell",'61' = "pDCs",'59'="stromal",'26'="stromal",
                         '34'="stromal",'56'="stromal",'54'="B cell",'0'="B cell", '35' = "monocytes",'63'="mregDC",'48'="cDC1",'41'="cDC2",'52'="Cycling M",'19'="IM",'27'="IM",'46'="SPP1 IM",'13'="FOLR2 IM",
                         '58'="IM",'49'="TREM2 IM",'40'="TREM2 IM",'14'="TREM2 IM",'16'="TREM2 IM",'24'="TREM2 IM",'50'="TREM2 IM",'44'="stromal",'12'="NK",'32'="NK",'17'= "gd T",'51' = "gd T",'25'="gd T",'21'="gd T",'23'="gd T")


human.tb <- RenameIdents(immune.combined.sct, '52' = "Proliferating T", '14' = "AM",'16' = "AM", '47' = "AM",'55' = "TREM2 IM",'56' = "TREM2 IM",'25' = "TREM2 IM", '29' = "TREM2 IM",'62' = "TREM2 IM",'7' = "FOLR2 IM",
                         '27' = "IM",'43' = "IM",'57' = "Proliferating Mac", '44' = "cDC2", '53' = "cDC1", '70' = "mregDC",'35' = "monocytes", '68' = "pDCs", '69' = "plasma cell",'60' = "B cell",'0' = "B cell", '67' = "mast cell",
                         '71' = "Treg", '46' = "Treg",'8' = "Treg", '54' = "SPP1 IM")

human.tb <- RenameIdents(immune.combined.sct, '52' = "Proliferating T", '14' = "TREM2 IM",'16' = "TREM2 IM", '47' = "TREM2 IM",'55' = "TREM2 IM",'56' = "TREM2 IM",'25' = "TREM2 IM", '29' = "TREM2 IM",'62' = "TREM2 IM",'7' = "FOLR2 IM",
                         '27' = "IM",'43' = "IM",'57' = "Proliferating Mac", '44' = "cDC2", '53' = "cDC1", '70' = "mregDC",'35' = "monocytes", '68' = "pDCs", '69' = "plasma cell",'60' = "B cell",'0' = "B cell", '67' = "mast cell",
                         '71' = "Treg", '46' = "Treg",'8' = "Treg", '54' = "SPP1 IM")


#Super simplified
simple.human.tb <- RenameIdents(immune.combined.sct, '52' = "Proliferating T", '14' = "TREM2 IM",'16' = "TREM2 IM", '47' = "TREM2 IM",'55' = "TREM2 IM",'56' = "TREM2 IM",'25' = "TREM2 IM", '29' = "TREM2 IM",'62' = "TREM2 IM",'7' = "FOLR2 IM",
                         '27' = "IM",'43' = "IM",'57' = "Proliferating Mac", '44' = "cDC2", '53' = "cDC1", '70' = "mregDC",'35' = "monocytes", '68' = "pDCs", '69' = "plasma cell",'60' = "B cell",'0' = "B cell", '67' = "mast cell",
                         '71' = "Treg", '46' = "Treg",'8' = "Treg", '54' = "SPP1 IM", '45' = "stromal", '66' = "stromal", '48' = "stromal", '49' = "stromal", '33' = "stromal", '39' = "ILC", '9' = "NK",'30' = "NK",'36'="NK",'13' = "gd T", '28' = "gd T",
                         '37' = "gd T",'20' = "gd T",'51' = "gd T",'31' = "CD4 T",'65' = "CD4 T",'11' = "CD4 T",'58' = "CD4 T",'32' = "CD8 T",'42' = "CD8 T",'19' = "CD8 T",'18' = "CD8 T",'15' = "CD8 T",'12' = "CD8 T",'50' = "CD8 T",'41' = "CD8 T",
                         '2' = "CD4 T",'23' = "CD4 T",'59' = "CD4 T",'24' = "CD4 T",'22' = "CD4 T",'17' = "CD4 T",'4' = "CD4 T",'72' = "CD4 T",'26' = "CD4 T",'21' = "CD4 T",'61' = "CD4 T",'40' = "CD4 T",'3' = "CD4 T",'5' = "CD4 T",'1' = "CD4 T",
                         '34' = "CD4 T",'63' = "CD4 T",'6' = "CD4 T",'10' = "CD4 T",'38' = "CD4 T",'64' = "stromal")

levels(human.tb) <- c("Cycling T","TREM2 IM","ILC","B cell","NK","FOLR2 IM","Mast","monocytes","Treg","pDCs","cDC2","CD8 T","IM","plasma cell","Cycling M","CD4 T","SPP1 IM","mregDC","gd T","stromal","cDC1")
DimPlot(simple.human.tb, label = T)


#CD8 T cell - 12, 15, 19, 24, 32, 42, 59
#CD4 T cell

#31, 65, 11, 58 are RORC+ CCR6+ SELL- CD4+ T cells

#previous analysis:

#48 is pDCs
#10 is Tregs
#6 is B cells
#45 is plasma cells
#46 is mast cells
#24 is stromal cells
#31 is cDC2
#39 is cDC1
#36 is ILC - Maybe ILC3
#15 has some ILC markers like KIT, TTLL10, TRDC, KLRB1 but not IL1R1 - ILC1
#22 NK cell
#25 is gd T cells


#34 is monocytes
#12 is SPP1+ Macs
#16 is IM
#11, 9 TREM2+ IM
#13 and 30AM


#Probably some kind of NK cells - 15, 22, 8, 25, 27, 42 and maybe 36

#Resolution of 3
#High for C1QA = 18, 19, 34, 41 and intermediate at 32, 45, 47
#TREM2 consistently high at 41

#32 and 19 are AM
#41 is FOLR2 IM
#34 18 and 5? are TREM2+ IM
#48 is SPP1 IM
#40 is Mono
 
human.tb <- RenameIdents(immune.combined.sct, '48' = "pDC", '10' = "Treg", '6' = "CXCR4 B", '47' = "Follicular B", '45' = "Plasma", '46' = "Mast", '24' = "stromal", '31' = "cDC2", '39' = "cDC1",'36' = "ILC3", '15' = "ILC1",
                         '22' = "NK", '25' = "gd T", '34' = "mono", '12' = "SPP1 IM", '16' = "IM", '11' = "TREM2 IM", '9' = "TREM2 IM", '13' = "AM", '30' = "AM")

human.tb$celltype <- Idents(human.tb)
#selectLab = c("Cxcl9","Spp1","Il1rn","C1qa","C1qc","C1qb","Trem2","Ccl5","Mmp2","Tnfsf13b","Lyz2","Apoe","Ubd","Dnase1l3","Saa3","C3","Cd9","Inhba","Lgals3","Cd53","Upp1","Il7r","Cxcl3","Ccl4","Ccl2","Ccr1","Phlda1","Slfn4","Clec4a1","C1ra","Clec4n","Edn1","Cd24a","Basp1")
#selectLab = c("SPP1","TREM2","C1QC","C1QB","AXL","IDO1","IL18BP","CD274","IL6","IL1RN","EGLN3","NDRG1","CLEC6A","ADAM8","CHST11", "MMP2","APOC2","CD101","LY86","GPR155","ITGB5","ITGA9","GLO1","CCL18")

saveRDS(human.tb, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Human TB resected lung/human.tb")

#Compare SPP1 IM and TREM2 IM
Idents(human.tb) <- "FDG.Signal"
human.low <- subset(human.tb, idents = "Low")
human.high <- subset(human.tb, idents = "High")

DefaultAssay(human.low) <- "SCT"
human.low <- PrepSCTFindMarkers(human.low)
DefaultAssay(human.high) <- "SCT"
human.high <- PrepSCTFindMarkers(human.high)

Idents(human.low) <- "celltype"
Idents(human.high) <- "celltype"
human.high.SPP1.TREM2.genes <- FindMarkers(human.high, ident.1 = "SPP1 IM", ident.2 = "TREM2 IM")
human.high.SPP1.TREM2.volcano <- human.high.SPP1.TREM2.genes[c(2,5)]
EnhancedVolcano(human.high.SPP1.TREM2.volcano,
                lab = rownames(human.high.SPP1.TREM2.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "SPP1 IM vs TREM2 IM",
                subtitle = "",
                drawConnectors = T,
                arrowheads = F,
                selectLab = c("SPP1", "TREM2","IL1RN","IL18BP","IDO1","C1QB","C1QC","C1QA","APOC1","APOE","CCL18","ALCAM","IL6","CD274","AXL","IL10","AREG","OSM","IL1R2"),
                raster = T)

human.low.SPP1.TREM2.genes <- FindMarkers(human.low, ident.1 = "SPP1 IM", ident.2 = "TREM2 IM")
human.low.SPP1.TREM2.volcano <- human.low.SPP1.TREM2.genes[c(2,5)]
EnhancedVolcano(human.low.SPP1.TREM2.volcano,
                lab = rownames(human.low.SPP1.TREM2.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "SPP1 IM vs TREM2 IM",
                subtitle = "",
                drawConnectors = T,
                arrowheads = F,
                selectLab = c("SPP1", "TREM2","IL1RN","IL18BP","IDO1","C1QB","C1QC","C1QA","APOC1","APOE","CCL18","ALCAM","IL6","CD274","AXL","IL10","AREG","OSM","IL1R2"),
                raster = T)


#IL1R2 OSM AREG IL10 IL1RN SPP1 IL6

EnhancedVolcano(human.SPP1.TREM2.volcano,
                lab = rownames(human.SPP1.TREM2.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "SPP1 IM vs TREM2 IM",
                subtitle = "",
                drawConnectors = T,
                arrowheads = F,
                selectLab = c("SPP1", "TREM2","IL1RN","IDO1","C1QB","C1QC","C1QA","APOC1","APOE","CCL18","CD274","AXL","IL10","OSM","IL1R2","CD9","FCN1","APOBEC3A","S100A9","CD93","MPEG1","SELL","FABP4","ALDH2","MRC1"),
                raster = T)

#Add FCN1 APOBEC3A S100A9 CD93 MPEG1 SELL
#Add FABP4 CD9 ALDH2 MRC1 

human.SPP1.TREM2.genes <- FindMarkers(human.tb, ident.1 = "SPP1 IM", ident.2 = "TREM2 IM")
human.SPP1.TREM2.volcano <- human.SPP1.TREM2.genes[c(2,5)]
EnhancedVolcano(human.SPP1.TREM2.volcano,
                lab = rownames(human.SPP1.TREM2.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "SPP1 IM vs TREM2 IM",
                subtitle = "",
                drawConnectors = T,
                arrowheads = F,
                selectLab = c("SPP1", "TREM2","IL1RN","IDO1","C1QB","C1QC","C1QA","APOC1","APOE","CCL18","CD274","AXL","IL10","OSM","IL1R2","FCN1","APOBEC3A","MPEG1","SRGN","METRNL","SCIMP","LILRA5","CD93",
                              "CCR2","IGF1","MARCO","CXCL16","LGALS3","IFITM3","CD59","ALDH2","BASP1"),
                raster = T)

EnhancedVolcano(human.SPP1.TREM2.volcano,
                lab = rownames(human.SPP1.TREM2.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "SPP1 IM vs TREM2 IM",
                subtitle = "",
                drawConnectors = T,
                arrowheads = F,
                raster = T)

human.SPP1.FOLR2.genes <- FindMarkers(human.tb, ident.1 = "SPP1 IM", ident.2 = "FOLR2 IM")
human.SPP1.FOLR2.volcano <- human.SPP1.FOLR2.genes[c(2,5)]
EnhancedVolcano(human.SPP1.FOLR2.volcano,
                lab = rownames(human.SPP1.FOLR2.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "SPP1 IM vs FOLR2 IM",
                subtitle = "",
                drawConnectors = T,
                arrowheads = F,
                raster = T)

human.SPP1.IM.genes <- FindMarkers(human.tb, ident.1 = "SPP1 IM", ident.2 = "IM")
human.SPP1.IM.volcano <- human.SPP1.IM.genes[c(2,5)]
EnhancedVolcano(human.SPP1.IM.volcano,
                lab = rownames(human.SPP1.IM.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "SPP1 IM vs IM",
                subtitle = "",
                drawConnectors = T,
                arrowheads = F,
                selectLab = c("SPP1", "TREM2","IL1RN","IDO1","C1QB","C1QC","C1QA","APOC1","APOE","CCL18","ALCAM","CD274","AXL","IL10","OSM","IL1R2","CD9"),
                raster = T)


selectLab = c("SPP1", "TREM2","IL1RN","IDO1","C1QB","C1QC","C1QA","APOC1","APOE","CCL18","ALCAM","CD274","AXL","IL10","OSM","IL1R2","CD9")

#subsetting
myeloid <- subset(human.tb, idents = c("IM", "SPP1 IM", "TREM2 IM", "FOLR2 IM","monocytes"))
levels(myeloid) <- c("SPP1 IM","FOLR2 IM","TREM2 IM","IM","monocytes")
VlnPlot(myeloid, features = "IL1RN", idents = c("monocytes","IM","TREM2 IM", "FOLR2 IM","SPP1 IM"), split.by = "FDG.signal")
DefaultAssay(myeloid) <- "SCT"
myeloid <- PrepSCTFindMarkers(myeloid)
myeloid$celltype_FDG <- paste(myeloid$celltype, myeloid$FDG.signal, sep = "_")
Idents(myeloid) <- "celltype_FDG"
DotPlot(myeloid, features = c("SPP1","TREM2","C1QA","FOLR2","SELENOP","IL1RN","IL18BP","IDO1","CD274","IL1R2","IL10","CD63","CD9"))


#3 cells clustered weirdly and look like background noise so they were removed from the analysis.
plot <- DimPlot(myeloid)
select.cells <- CellSelector(plot = plot) #"Array6_29618_CGGTGTGGCAGC" "Array2_29818_ACGATTGGAACT" "Array2_29818_TGGCTTGGAGAT"
select.cells <- c("Array6_29618_CGGTGTGGCAGC","Array2_29818_ACGATTGGAACT","Array2_29818_TGGCTTGGAGAT")
myeloid_filtered <- myeloid[,!colnames(myeloid) %in% select.cells]
DimPlot(myeloid_filtered)

human.markers <- FindAllMarkers(human.tb)
levels(human.tb) <- c("stromal","Mast","plasma cell","B cell","pDCs","mregDC","cDC2","cDC1","gd T","ILC","NK","CD8 T","Treg","CD4 T","Cycling T","SPP1 IM","FOLR2 IM","TREM2 IM","IM","monocytes","Cycling M")
DotPlot(human.tb, features = c("MKI67","MAFB","FCGR3A","CD14","ITGAX","TREM2","FOLR2","SPP1","CD3G","CD4","CD8A","FOXP3","KLRD1","IL7R","TRGV9","FLT3","CLEC9A","CLEC10A","CD200","CLEC4C","CD79B","MZB1","FCER1A","EPCAM"))
levels(human.tb) <- c("stromal","Mast","plasma cell","B cell","pDCs","mregDC","cDC1","cDC2","monocytes","IM","FOLR2 IM","TREM2 IM","SPP1 IM","gd T","ILC","NK","CD8 T","Treg","CD4 T","Cycling T","Cycling M")
