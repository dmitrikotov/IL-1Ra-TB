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
immune.combined.sct <- FindClusters(immune.combined.sct, resolution = 3)

#Work with the SCT corrected values rather than the integrated values
DefaultAssay(immune.combined.sct) <- "SCT"
immune.combined.sct <- PrepSCTFindMarkers(immune.combined.sct)

saveRDS(immune.combined.sct, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Human TB resected lung/Human_lungv2.RDS")

human.tb.lung.markers <- FindAllMarkers(immune.combined.sct)

#V2 Resolution of 3
Idents(immune.combined.sct) <- "integrated_snn_res.3"

immune.combined.sct <- readRDS(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Human TB resected lung/Human_lungv2.RDS")

human.tb <- RenameIdents(immune.combined.sct, '60' = "Mast", '37' = "ILC", '47' = "Cycling T", '4' = "Treg",'45' = "CD8 T",'39' = "CD8 T",'30' = "CD4 T", '57' = "CD4 T", '15' = "CD4 T",'53' = "CD4 T",
                         '31' = "CD8 T",'20' = "CD8 T",'22' = "CD8 T",'36' = "CD8 T",'18' = "CD8 T",'9' = "CD8 T",'33'="CD4 T",'55'="CD4 T",'6'="CD4 T",'2'="CD4 T",'1'="CD4 T",'10'="CD4 T",'29'="CD4 T",
                         '42' = "CD4 T",'3' = "CD4 T",'8'="CD4 T",'11'="CD4 T",'38'="CD4 T",'28'="CD4 T",'7'="CD4 T",'5'="CD4 T",'43'="CD4 T",'62' = "plasma cell",'61' = "pDCs",'59'="stromal",'26'="stromal",
                         '34'="stromal",'56'="stromal",'54'="B cell",'0'="B cell", '35' = "Non-classical",'63'="mregDC",'48'="cDC1",'41'="cDC2",'52'="Cycling M",'19'="monocyte",'27'="IM",'46'="SPP1 IM",'13'="FOLR2 IM",
                         '58'="IM",'49'="AM",'40'="AM",'14'="AM",'16'="AM",'24'="TREM2 IM",'50'="AM",'44'="stromal",'12'="NK",'32'="NK",'17'= "gd T",'51' = "gd T",'25'="gd T",'21'="gd T",'23'="gd T")

levels(human.tb) <- c("AM","Cycling T","TREM2 IM","ILC","B cell","NK","FOLR2 IM","Mast","Non-classical","Treg","pDCs","cDC2","CD8 T","monocyte","plasma cell","Cycling M","CD4 T","SPP1 IM","mregDC","gd T","stromal","cDC1")


human.tb$celltype <- Idents(human.tb)

saveRDS(human.tb, file = "/Users/dmitrikotov/Downloads/Box Data - Postdoc/Old Stuff/Postdoc Files/Coding stuff/Human TB resected lung/human.tb")
human.tb <- readRDS(file = '/Users/dmitrikotov/Library/CloudStorage/Box-Box/Dmitri Data/Coding stuff/Human TB resected lung/Human_lungv2.RDS')

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
                selectLab = c("SPP1", "TREM2","IL1RN","IDO1","C1QB","C1QC","C1QA","APOC1","APOE","CCL18","AXL","IL27","TGFB1","OSM","IL1R2","FCN1","APOBEC3A","MPEG1","SRGN","METRNL","SCIMP","LILRA5","CD93",
                              "CCR2","IGF1","MARCO","CXCL16","LGALS3","IFITM3","CD59","ALDH2","IL10","CD274"),
                raster = T,
                pCutoff = 0.05)


#subsetting
myeloid <- subset(human.tb, idents = c("monocyte", "SPP1 IM", "TREM2 IM", "FOLR2 IM","Non-classical","AM"))
levels(myeloid) <- c("AM","Non-classical","monocyte","TREM2 IM","FOLR2 IM","SPP1 IM")
VlnPlot(myeloid, features = "IL1RN", idents = c("AM","Non-classical","monocyte","TREM2 IM", "FOLR2 IM","SPP1 IM"), split.by = "FDG.signal")
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
levels(human.tb) <- c("stromal","Mast","plasma cell","B cell","pDCs","mregDC","cDC2","cDC1","gd T","ILC","NK","CD8 T","Treg","CD4 T","Cycling T","AM","SPP1 IM","FOLR2 IM","TREM2 IM","monocyte","Non-classical","Cycling M")
DotPlot(human.tb, features = c("MKI67","MAFB","FCGR3A","CSF1R","CD14","ITGAX","TREM2","FOLR2","SPP1","FABP4","CD3G","CD4","CD8A","FOXP3","KLRD1","IL7R","TRGV9","FLT3","CLEC9A","CLEC10A","CD200","CLEC4C","CD79B","MZB1","FCER1A","EPCAM"))
levels(human.tb) <- c("stromal","Mast","plasma cell","B cell","pDCs","mregDC","cDC1","cDC2","monocytes","IM","FOLR2 IM","TREM2 IM","SPP1 IM","AM","gd T","ILC","NK","CD8 T","Treg","CD4 T","Cycling T","Cycling M")