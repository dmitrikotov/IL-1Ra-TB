library(Seurat)
library(tidyverse)
library(patchwork)
library(cowplot)
library(EnhancedVolcano)

# Load the infected cell dataset
il1.mtb.rna <- Read10X(data.dir = "/Users/DmitriKotov/Documents/Sequencing Data/scRNA-seq 7-22-21/RVDK004A/outs/filtered_feature_bc_matrix", strip.suffix = TRUE)
wt.mtb.rna <- Read10X(data.dir = "/Users/DmitriKotov/Documents/Sequencing Data/scRNA-seq 7-22-21/RVDK004C/outs/filtered_feature_bc_matrix", strip.suffix = TRUE)
il1.naive.rna  <- Read10X(data.dir = "/Users/DmitriKotov/Documents/Sequencing Data/scRNA-seq 7-22-21/RVDK004B/outs/filtered_feature_bc_matrix", strip.suffix = TRUE)
wt.naive.rna <- Read10X(data.dir = "/Users/DmitriKotov/Documents/Sequencing Data/scRNA-seq 7-22-21/RVDK004D/outs/filtered_feature_bc_matrix", strip.suffix = TRUE)

# Initialize the Seurat object with the raw (non-normalized data).
wt.mtb <- CreateSeuratObject(counts = wt.mtb.rna, project = "stromal")
wt.naive <- CreateSeuratObject(counts = wt.naive.rna, project = "stromal")
il1.mtb <- CreateSeuratObject(counts = il1.mtb.rna, project = "stromal")
il1.naive <- CreateSeuratObject(counts = il1.naive.rna, project = "stromal")
wt.mtb #check seurat object

#check for dying cells based on mitochondrial percentage
wt.mtb[["percent.mt"]] <- PercentageFeatureSet(wt.mtb, pattern = "^mt-")
wt.naive[["percent.mt"]] <- PercentageFeatureSet(wt.naive, pattern = "^mt-")
il1.mtb[["percent.mt"]] <- PercentageFeatureSet(il1.mtb, pattern = "^mt-")
il1.naive[["percent.mt"]] <- PercentageFeatureSet(il1.naive, pattern = "^mt-")

#standard filters are 
VlnPlot(wt.mtb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# standard log-normalization
wt.mtb <- NormalizeData(wt.mtb)
wt.naive <- NormalizeData(wt.naive)
il1.mtb <- NormalizeData(il1.mtb)
il1.naive <- NormalizeData(il1.naive)

# choose ~1k variable features
wt.mtb <- FindVariableFeatures(wt.mtb)
wt.naive <- FindVariableFeatures(wt.naive)
il1.mtb <- FindVariableFeatures(il1.mtb)
il1.naive <- FindVariableFeatures(il1.naive)

# standard scaling (no regression)
wt.mtb <- ScaleData(wt.mtb)
wt.naive <- ScaleData(wt.naive)
il1.mtb <- ScaleData(il1.mtb)
il1.naive <- ScaleData(il1.naive)

#Filter to exclude cells that have unique feature counts more than 2500 (multiplets) and less than 200 (empty drops) and removing cells that have >5% mitochondrial counts (dead/dying cells)
wt.mtb <- subset(wt.mtb, subset = nFeature_RNA < 2500 & nFeature_RNA > 200 & percent.mt < 5)
wt.naive <- subset(wt.naive, subset = nFeature_RNA < 2500 & nFeature_RNA > 200 & percent.mt < 5)
il1.mtb <- subset(il1.mtb, subset = nFeature_RNA < 2500 & nFeature_RNA > 200 & percent.mt < 5)
il1.naive <- subset(il1.naive, subset = nFeature_RNA < 2500 & nFeature_RNA > 200 & percent.mt < 5)

#Add metadata to each group
wt.mtb[["genotype"]] <- "WT"
wt.naive[["genotype"]] <- "WT"
il1.mtb[["genotype"]] <- "IL1KO"
il1.naive[["genotype"]] <- "IL1KO"
wt.mtb[["state"]] <- "Mtb"
wt.naive[["state"]] <- "naive"
il1.mtb[["state"]] <- "Mtb"
il1.naive[["state"]] <- "naive"

#continue with integrated analysis here: #what is appropriate dims for integration?
object.list <-c(wt.mtb, wt.naive, il1.mtb, il1.naive)
immune.anchors <- FindIntegrationAnchors(object.list = object.list, dims = 1:30)
stromal.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:30)
DefaultAssay(stromal.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
stromal.combined <- ScaleData(stromal.combined, verbose = FALSE)
stromal.combined <- RunPCA(stromal.combined, npcs = 30, verbose = FALSE)

# t-SNE and Clustering
stromal.combined <- RunUMAP(stromal.combined, reduction = "pca", dims = 1:30)
stromal.combined <- FindNeighbors(stromal.combined, reduction = "pca", dims = 1:30)
stromal.combined <- FindClusters(stromal.combined, resolution = 1.5)

#Visualize integrated dataset
p1 <- DimPlot(stromal.combined, reduction = "umap", group.by = "state")
p2 <- DimPlot(stromal.combined, reduction = "umap", label = TRUE)
p3 <- DimPlot(stromal.combined, reduction = "umap", group.by = "genotype")
plot_grid(p1, p2, p3)

DimPlot(stromal.combined, reduction = "umap", split.by = "genotype")
DimPlot(stromal.combined, reduction = "umap", split.by = "state")

#Find markers for each cluster
DefaultAssay(stromal.combined) <- "RNA"
all.markers <- FindAllMarkers(stromal.combined)


#Identifying markers for each cluster:
DotPlot(stromal.combined, features = c("Lyve1","Pecam1","Epcam","Col1a1","Ptprc","Msln","Sftpc","Akap5","Scgb1a1","Col14a1","Col13a1","Tgfbi","Notch3","Car4","Gpihbp1","Bmx","Slc6a2","Esm1","Ly6g","Ccr2","Fcgr1","Itgax","Ear1","Cd3e","Sell","Cd4","Cd8a","Cd79a","Gzma"))

#re-label clusters
#Reference for cluster labels is Hurskainen M Nature Communications 2021
stromal.combined <- RenameIdents(stromal.combined, `0` = "gCap 2", `1` = "gCap 1", `2` = "Col13a1+ fibroblast 1", 
                                 `3` = "Myofibroblast 2", `4` = "Myofibroblast 3", `5` = "Club 2", `6` = "Club 1", `7` = "Col13a1+ fibroblast 2", `8` = "gCap 3", `9` = "B", 
                                 `10` = "Art", `11` = "Vein", `12` = "Myofibroblast 1", `13` = "Naive T",`14` = "Effector CD4 T",`15` = "aCap 1",
                                 `16` = "Col13a1+ fibroblast 3",`17` = "Effector CD8 T",`18` = "Col13a1+ fibroblast 4",`19` = "Col14a1 fibroblast 1",`20` = "Monocyte",`21` = "Pericyte",
                                 `22` = "Fibromyo/SMC",`23` = "Lymph",`24` = "Esm1 Endothelial",`25` = "Neutrophil",`26` = "Col14a1 fibroblast 2",`27` = "AM",`28` = "aCap 2",`29` = "AT2",
                                 `30` = "NK",`31` = "IM",`32` = "AT1",`33` = "Scgb3a1+ Club",`34` = "Stickler Syndrome",`35` = "Myofibroblast 4")

#Ciliated cells might be unclustered offshoot of cluster 6.
#Col13a1+ fibroblast 2 looks different from rest of Col13a1+ fibroblasts

stromal.combined$old.ident <- Idents(stromal.combined)
stromal.combined$state_genotype <- paste(stromal.combined$state, stromal.combined$genotype, sep = "_")
stromal.combined$celltype_state_genotype <- paste(Idents(stromal.combined),stromal.combined$state, stromal.combined$genotype, sep = "_")

Idents(stromal.combined) <- "old.ident"
levels(stromal.combined) <- c("gCap 1","gCap 2", "gCap 3","aCap 1","aCap 2","Art","Vein","Esm1 Endothelial",
                              "AT1","Club 1","Club 2","AT2","Scgb3a1+ Club","Stickler Syndrome","B","NK","Neutrophil","Monocyte","IM","AM","Naive T","Effector CD4 T","Effector CD8 T",
                              "Col13a1+ fibroblast 1","Col13a1+ fibroblast 2","Col13a1+ fibroblast 3","Col13a1+ fibroblast 4","Col14a1 fibroblast 1","Col14a1 fibroblast 2",
                              "Pericyte","Fibromyo/SMC","Myofibroblast 1","Myofibroblast 2","Myofibroblast 3","Myofibroblast 4","Lymph")
DimPlot(stromal.combined, label = TRUE, repel = TRUE)

simplified.stromal <- stromal.combined
simplified.stromal <- RenameIdents(simplified.stromal, 'Club 1' = "Club",'Club 2' = "Club",'Myofibroblast 1' = "Myofibroblast",'Myofibroblast 2' = "Myofibroblast",
                                   'Myofibroblast 3' = "Myofibroblast",'Myofibroblast 4' = "Myofibroblast",'Col14a1 fibroblast 1' = "Col14a1 fibroblast", 'Col14a1 fibroblast 2' = "Col14a1 fibroblast",
                                   'Col13a1+ fibroblast 1' = "Col13a1+ fibroblast",'Col13a1+ fibroblast 2' = "Col13a1+ fibroblast",'Col13a1+ fibroblast 3' = "Col13a1+ fibroblast",
                                   'Col13a1+ fibroblast 4' = "Col13a1+ fibroblast",'gCap 1' = "gCap",'gCap 2' = "gCap",'gCap 3' = "gCap",'aCap 1' = "aCap",'aCap 2' = "aCap")
DimPlot(simplified.stromal, label = TRUE, repel = TRUE)
simplified.stromal$old.ident <- Idents(simplified.stromal)
simplified.stromal$celltype_state_genotype <- paste(Idents(simplified.stromal),simplified.stromal$state, simplified.stromal$genotype, sep = "_")

#Do DEG of various populations
Idents(simplified.stromal) <- "celltype_state_genotype"
AT2.il1.markers <- FindMarkers(simplified.stromal, assay = "RNA", ident.1 = "AT2_Mtb_WT", ident.2 = "AT2_Mtb_IL1KO")
AT2.il1.volcano <- AT2.il1.markers[c(2,5)]
EnhancedVolcano(AT2.il1.volcano,
                lab = rownames(AT2.il1.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Mtb infected WT vs IL-1R KO AT2 cells",
                subtitle = "",
                pCutoff = 0.05)

AT2.mtb.markers <- FindMarkers(simplified.stromal, assay = "RNA", ident.1 = "AT2_Mtb_WT", ident.2 = "AT2_naive_WT")
AT2.mtb.volcano <- AT2.mtb.markers[c(2,5)]
EnhancedVolcano(AT2.mtb.volcano,
                lab = rownames(AT2.mtb.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Mtb infected vs uninfected WT AT2 cells",
                subtitle = "",
                pCutoff = 0.05)

Col14a1.il1.markers <- FindMarkers(simplified.stromal, assay = "RNA", ident.1 = "Col14a1 fibroblast_Mtb_WT", ident.2 = "Col14a1 fibroblast_Mtb_IL1KO")
Col14a1.il1.volcano <- Col14a1.il1.markers[c(2,5)]
EnhancedVolcano(Col14a1.il1.volcano,
                lab = rownames(Col14a1.il1.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Mtb infected WT vs IL-1R KO Col14a1 fibroblast cells",
                subtitle = "",
                pCutoff = 0.05)
Col14a1.il1.volcano <- Col14a1.il1.markers[c(1,2)]
EnhancedVolcano(Col14a1.il1.volcano,
                lab = rownames(Col14a1.il1.volcano),
                x = 'avg_log2FC',
                y = 'p_val',
                title = "Mtb infected WT vs IL-1R KO Col14a1 fibroblast cells",
                subtitle = "",
                pCutoff = 0.05,
                selectLab = c("Ccl7","Ccl2","Cxcl1","Scgb1a1","Rpl29","Icam1","Nars","Cryab","Cxcl10","Csf1","Sod2"))


Col14a1.mtb.markers <- FindMarkers(simplified.stromal, assay = "RNA", ident.1 = "Col14a1 fibroblast_Mtb_WT", ident.2 = "Col14a1 fibroblast_naive_WT")
Col14a1.mtb.volcano <- Col14a1.mtb.markers[c(2,5)]
EnhancedVolcano(Col14a1.mtb.volcano,
                lab = rownames(Col14a1.mtb.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Mtb infected vs naive WT Col14a1 fibroblast cells",
                subtitle = "",
                pCutoff = 0.05)

Col13a1.il1.markers <- FindMarkers(simplified.stromal, assay = "RNA", ident.1 = "Col13a1+ fibroblast_Mtb_WT", ident.2 = "Col13a1+ fibroblast_Mtb_IL1KO")
Col13a1.il1.volcano <- Col13a1.il1.markers[c(2,5)]
EnhancedVolcano(Col13a1.il1.volcano,
                lab = rownames(Col13a1.il1.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Mtb infected WT vs IL-1R KO Col13a1 fibroblast cells",
                subtitle = "",
                pCutoff = 0.05)
Col13a1.il1.volcano <- Col13a1.il1.markers[c(1,2)]
EnhancedVolcano(Col13a1.il1.volcano,
                lab = rownames(Col13a1.il1.volcano),
                x = 'avg_log2FC',
                y = 'p_val',
                title = "Mtb infected WT vs IL-1R KO Col13a1 fibroblast cells",
                subtitle = "",
                pCutoff = 0.05)


Col13a1.mtb.markers <- FindMarkers(simplified.stromal, assay = "RNA", ident.1 = "Col13a1+ fibroblast_Mtb_WT", ident.2 = "Col13a1+ fibroblast_naive_WT")
Col13a1.mtb.volcano <- Col13a1.mtb.markers[c(2,5)]
EnhancedVolcano(Col13a1.mtb.volcano,
                lab = rownames(Col13a1.mtb.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Mtb infected vs naive WT Col13a1 fibroblast cells",
                subtitle = "",
                pCutoff = 0.05)

Col13a1.mtb.il1.markers <- FindMarkers(simplified.stromal, assay = "RNA", ident.1 = "Col13a1+ fibroblast_Mtb_IL1KO", ident.2 = "Col13a1+ fibroblast_naive_IL1KO")
Col13a1.mtb.il1.volcano <- Col13a1.mtb.il1.markers[c(2,5)]
EnhancedVolcano(Col13a1.mtb.il1.volcano,
                lab = rownames(Col13a1.mtb.il1.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Mtb infected vs naive IL-1R KO Col13a1 fibroblast cells",
                subtitle = "",
                pCutoff = 0.05)

Myofibroblast.il1.markers <- FindMarkers(simplified.stromal, assay = "RNA", ident.1 = "Myofibroblast_Mtb_WT", ident.2 = "Myofibroblast_Mtb_IL1KO")
Myofibroblast.il1.volcano <- Myofibroblast.il1.markers[c(2,5)]
EnhancedVolcano(Myofibroblast.il1.volcano,
                lab = rownames(Myofibroblast.il1.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Mtb infected WT vs IL-1R KO Myofibroblast cells",
                subtitle = "",
                pCutoff = 0.05)
Myofibroblast.il1.volcano <- Myofibroblast.il1.markers[c(1,2)]
EnhancedVolcano(Myofibroblast.il1.volcano,
                lab = rownames(Myofibroblast.il1.volcano),
                x = 'avg_log2FC',
                y = 'p_val',
                title = "Mtb infected WT vs IL-1R KO Myofibroblast cells",
                subtitle = "",
                pCutoff = 0.05)

Myofibroblast.mtb.markers <- FindMarkers(simplified.stromal, assay = "RNA", ident.1 = "Myofibroblast_Mtb_WT", ident.2 = "Myofibroblast_naive_WT")
Myofibroblast.mtb.volcano <- Myofibroblast.mtb.markers[c(2,5)]
EnhancedVolcano(Myofibroblast.mtb.volcano,
                lab = rownames(Myofibroblast.mtb.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Mtb infected vs naive WT Myofibroblast cells",
                subtitle = "",
                pCutoff = 0.05)

gCap.il1.markers <- FindMarkers(simplified.stromal, assay = "RNA", ident.1 = "gCap_Mtb_WT", ident.2 = "gCap_Mtb_IL1KO")
gCap.il1.volcano <- gCap.il1.markers[c(2,5)]
EnhancedVolcano(gCap.il1.volcano,
                lab = rownames(gCap.il1.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Mtb infected WT vs IL-1R KO gCap cells",
                subtitle = "",
                pCutoff = 0.05)
gCap.il1.volcano <- gCap.il1.markers[c(1,2)]
EnhancedVolcano(gCap.il1.volcano,
                lab = rownames(gCap.il1.volcano),
                x = 'avg_log2FC',
                y = 'p_val',
                title = "Mtb infected WT vs IL-1R KO gCap cells",
                subtitle = "",
                pCutoff = 0.05)

gCap.mtb.markers <- FindMarkers(simplified.stromal, assay = "RNA", ident.1 = "gCap_Mtb_WT", ident.2 = "gCap_naive_WT")
gCap.mtb.volcano <- gCap.mtb.markers[c(2,5)]
EnhancedVolcano(gCap.mtb.volcano,
                lab = rownames(gCap.mtb.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Mtb infected vs naive WT gCap cells",
                subtitle = "",
                pCutoff = 0.05)

aCap.il1.markers <- FindMarkers(simplified.stromal, assay = "RNA", ident.1 = "aCap_Mtb_WT", ident.2 = "aCap_Mtb_IL1KO")
aCap.il1.volcano <- aCap.il1.markers[c(2,5)]
EnhancedVolcano(aCap.il1.volcano,
                lab = rownames(aCap.il1.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Mtb infected WT vs IL-1R KO aCap cells",
                subtitle = "",
                pCutoff = 0.05)
aCap.il1.volcano <- aCap.il1.markers[c(1,2)]
EnhancedVolcano(aCap.il1.volcano,
                lab = rownames(aCap.il1.volcano),
                x = 'avg_log2FC',
                y = 'p_val',
                title = "Mtb infected WT vs IL-1R KO aCap cells",
                subtitle = "",
                pCutoff = 0.05)

aCap.mtb.markers <- FindMarkers(simplified.stromal, assay = "RNA", ident.1 = "aCap_Mtb_WT", ident.2 = "aCap_naive_WT")
aCap.mtb.volcano <- aCap.mtb.markers[c(2,5)]
EnhancedVolcano(aCap.mtb.volcano,
                lab = rownames(aCap.mtb.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Mtb infected vs naive WT aCap cells",
                subtitle = "",
                pCutoff = 0.05)

club.il1.markers <- FindMarkers(simplified.stromal, assay = "RNA", ident.1 = "Club_Mtb_WT", ident.2 = "Club_Mtb_IL1KO")
club.il1.volcano <- club.il1.markers[c(2,5)]
EnhancedVolcano(club.il1.volcano,
                lab = rownames(club.il1.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Mtb infected WT vs IL-1R KO Club cells",
                subtitle = "",
                pCutoff = 0.05)

#Combine IL-1 receptor expressing cells genes up > 2 fold.
super.simp.stromal <- simplified.stromal
super.simp.stromal <- RenameIdents(super.simp.stromal, 'Pericyte' = "IL1_responsive", 'Myofibroblast' ="IL1_responsive", 'Col13a1+ fibroblast' = "IL1_responsive", 'Col14a1 fibroblast' = "IL1_responsive",
                                   'AT2' = "IL1_responsive", 'aCap' = "IL1_responsive", 'Art' = "IL1_responsive", 'gCap' = "IL1_responsive", 'Vein' = "IL1_responsive", 'Naive T' = "nonresponsive",
                                   'B' = "nonresponsive", 'NK' = "nonresponsive", 'Monocyte' = "nonresponsive", 'Neutrophil' = "nonresponsive", 'Club' = "nonresponsive", 'Scgb3a1+ Club' = "nonresponsive")
DimPlot(super.simp.stromal, label = TRUE, repel = TRUE)
super.simp.stromal$celltype_state_genotype <- paste(Idents(super.simp.stromal),super.simp.stromal$state, super.simp.stromal$genotype, sep = "_")

Idents(super.simp.stromal) <- "celltype_state_genotype"

il1.Mtb.markers <- FindMarkers(super.simp.stromal, assay = "RNA", ident.1 = "IL1_responsive_Mtb_WT", ident.2 = "IL1_responsive_Mtb_IL1KO")
il1.Mtb.volcano <- il1.Mtb.markers[c(2,5)]
EnhancedVolcano(il1.Mtb.volcano,
                lab = rownames(il1.Mtb.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Mtb infected WT vs IL-1R KO cells",
                subtitle = "",
                pCutoff = 0.05)

il1.responders.markers <- FindMarkers(super.simp.stromal, assay = "RNA", ident.1 = "IL1_responsive_Mtb_WT", ident.2 = "nonresponsive_Mtb_WT")
il1.responders.markers.2 <- FindMarkers(super.simp.stromal, assay = "RNA", ident.1 = "IL1_responsive_Mtb_IL1KO", ident.2 = "nonresponsive_Mtb_IL1KO")
il1.responders.markers$gene <- rownames(il1.responders.markers)
il1.responders.markers.2$gene <- rownames(il1.responders.markers.2)
il1.responders.markers <- anti_join(il1.responders.markers,il1.responders.markers.2, by = "gene")

il1.responders.volcano <- il1.responders.markers[c(2,5)]
EnhancedVolcano(il1.responders.volcano,
                lab = rownames(il1.responders.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Mtb infected IL-1R1 expressing vs non-expressing cells",
                subtitle = "",
                pCutoff = 0.05)

saveRDS(all.markers, file = "/Users/DmitriKotov/Box Sync/Coding stuff/Stromal Cell IL-1 scRNA-seq 7-22-21/allmarkers")
saveRDS(stromal.combined, file = "/Users/DmitriKotov/Box Sync/Coding stuff/Stromal Cell IL-1 scRNA-seq 7-22-21/stromalcombined")
stromal.combined <- readRDS(file = "/Users/DmitriKotov/library/CloudStorage/Box-Box/Coding stuff/Stromal Cell IL-1 scRNA-seq 7-22-21/stromalcombined")
all.markers <- readRDS(file = "/Users/DmitriKotov/Box Sync/Coding stuff/Stromal Cell IL-1 scRNA-seq 7-22-21/allmarkers")

#Add scores for TNF, IFNg, IFN sig - Looks good for IFNb and IFNg but too strict for tnf and tgfb.
library(UCell)
cytokine <- read.csv(file = "/Users/DmitriKotov/Box Sync/Coding stuff/B6 vs Sp140 scRNA-seq 10-9-20/exvivo_markers.csv")
cytokine.ifn <- na.omit(cytokine) %>% filter(avgExpr > 2) %>% filter(gsi > 1.5) #gsi > 1.3 for TGFb and TNF?
cytokine.ifnb <- filter(cytokine.ifn, group == "IFNB")  
cytokine.ifng <- filter(cytokine.ifn, group == "IFNG")  
cytokine.tnf <- filter(cytokine.ifn, group == "TNFA")  
cytokine.tgfb <- filter(cytokine.ifn, group == "TGFB") 
cytokine.null <- filter(cytokine.ifn, group == "null")

cyto.sig <- list()
cyto.sig$ifnb <- setdiff(cytokine.ifnb$mouse_symbol, cytokine.ifng$mouse_symbol) %>% setdiff(cytokine.tnf$mouse_symbol) %>% setdiff(cytokine.tgfb$mouse_symbol) %>% setdiff(cytokine.null$mouse_symbol)
cyto.sig$ifng <- setdiff(cytokine.ifng$mouse_symbol, cytokine.ifnb$mouse_symbol) %>% setdiff(cytokine.tnf$mouse_symbol) %>% setdiff(cytokine.tgfb$mouse_symbol) %>% setdiff(cytokine.null$mouse_symbol)

stromal.combined <- AddModuleScore_UCell(stromal.combined, features = cyto.sig)

Idents(stromal.combined) <- "old.ident"
VlnPlot(stromal.combined, features = c("ifnb_UCell","ifng_UCell"), split.by = "state", ncol = 2)
FeaturePlot(stromal.combined, features = c("ifnb_UCell","ifng_UCell"), split.by = "state")

#gene signature for TNF and TGFb
cytokine.2 <- na.omit(cytokine) %>% filter(avgExpr > 2) %>% filter(gsi > 1.3)
cytokine.ifnb.2 <- filter(cytokine.2, group == "IFNB")  
cytokine.ifng.2 <- filter(cytokine.2, group == "IFNG")  
cytokine.tnf.2 <- filter(cytokine.2, group == "TNFA")  
cytokine.tgfb.2 <- filter(cytokine.2, group == "TGFB") 
cytokine.null.2 <- filter(cytokine.2, group == "null")

cyto.sig.2 <- list()
cyto.sig.2$tnf <- setdiff(cytokine.tnf.2$mouse_symbol, cytokine.ifnb.2$mouse_symbol) %>% setdiff(cytokine.ifng.2$mouse_symbol) %>% setdiff(cytokine.tgfb.2$mouse_symbol) %>% setdiff(cytokine.null.2$mouse_symbol)
cyto.sig.2$tgfb <- setdiff(cytokine.tgfb.2$mouse_symbol, cytokine.ifnb.2$mouse_symbol) %>% setdiff(cytokine.ifng.2$mouse_symbol) %>% setdiff(cytokine.tnf.2$mouse_symbol) %>% setdiff(cytokine.null.2$mouse_symbol)

stromal.combined <- AddModuleScore_UCell(stromal.combined, features = cyto.sig.2)

Idents(stromal.combined) <- "old.ident"
VlnPlot(stromal.combined, features = c("tgfb_UCell"), split.by = "state", ncol = 2)
FeaturePlot(stromal.combined, features = c("tgfb_UCell"), split.by = "state")

#Save file after adding cytokine signaling analysis
saveRDS(stromal.combined, file = "/Users/DmitriKotov/Box Sync/Coding stuff/Stromal Cell IL-1 scRNA-seq 7-22-21/stromalcombined2")
stromal.combined <- readRDS(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Stromal Cell IL-1 scRNA-seq 7-22-21/stromalcombined2")

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
library(nichenetr)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

m.s.genes <-convert_human_to_mouse_symbols(s.genes)
m.g2m.genes <- convert_human_to_mouse_symbols(g2m.genes)

DefaultAssay(stromal.combined) <-"RNA"
stromal.combined <- CellCycleScoring(stromal.combined, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
Idents(stromal.combined) <- "Phase"
DimPlot(stromal.combined, label = TRUE)
Idents(stromal.combined) <- "old.ident"


#Interesting result - It looks like there is a type I interferon response in Club cells and Col14a1+ Fibroblasts, but a IFN-g response in capillary cells.
#Cutoffs for positive signal might be too stringent for this dataset. Play with it a bit.
stromal.combined$ifn_sig <- ifelse(stromal.combined$ifnb_UCell >= 0.35 & stromal.combined$ifng_UCell >= 0.4, "both", 
                                   ifelse(stromal.combined$ifnb_UCell >= 0.35 & stromal.combined$ifng_UCell < 0.4, "ifnb", 
                                          ifelse(stromal.combined$ifnb_UCell < 0.35 & stromal.combined$ifng_UCell >= 0.4, "ifng","neither")))

#Simplified Stromal Cell
simplified.stromal <- RenameIdents(stromal.combined,'Myofibroblast 1' = "Myofibroblast",'Myofibroblast 2' = "Myofibroblast", 'Myofibroblast 3' = "Myofibroblast", 'Myofibroblast 4' = "Myofibroblast",
                                   'Fibromyo/SMC' = "Myofibroblast", 'Pericyte' = "Myofibroblast", 'Col13a1+ fibroblast 1' = "Fibroblast", 'Col13a1+ fibroblast 2' = "Fibroblast",
                                   'Col13a1+ fibroblast 3' = "Fibroblast",'Col13a1+ fibroblast 4' = "Fibroblast", 'Col14a1 fibroblast 1' = "Fibroblast", 'Col14a1 fibroblast 2' = "Fibroblast",
                                   'aCap 1' = "Endothelial", 'aCap 2' = "Endothelial", 'gCap 3' = "Endothelial", 'Esm1 Endothelial' = "Endothelial", 'Art' = "Endothelial", 'gCap 2' = "Endothelial",
                                   'gCap 1' = "Endothelial", 'Vein' = "Endothelial", 'Club 1' = "Club", 'Club 2' = "Club")
simplified.stromal$celltype_state_genotype <- paste(Idents(simplified.stromal),simplified.stromal$state, simplified.stromal$genotype, sep = "_")

#Look at differentially expressed genes
potential.ligands <- unique(lr_network$from) %>% convert_human_to_mouse_symbols() %>% na.omit

Idents(simplified.stromal) <- "celltype_state_genotype"
Myo.il1.markers <- FindMarkers(simplified.stromal, assay = "RNA", ident.1 = "Myofibroblast_Mtb_WT", ident.2 = "Myofibroblast_Mtb_IL1KO")
Myo.il1.volcano <- Myo.il1.markers[c(2,5)]
regulated.Myo <- Myo.il1.volcano %>% filter(avg_log2FC > 1 | avg_log2FC < -1)  %>% row.names()
Myo.ligands <- intersect(potential.ligands, regulated.Myo)

p1 <- EnhancedVolcano(Myo.il1.volcano,
                lab = rownames(Myo.il1.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Mtb infected WT vs IL-1R KO Myofibroblast cells",
                subtitle = "",
                pCutoff = 0.05,
                xlim = c(-5,5),
                ylim = c(0,60),
                selectLab = Myo.ligands,
                drawConnectors = TRUE,
                arrowheads = FALSE,
                labSize = 10)

Fibro.il1.markers <- FindMarkers(simplified.stromal, assay = "RNA", ident.1 = "Fibroblast_Mtb_WT", ident.2 = "Fibroblast_Mtb_IL1KO")
Fibro.il1.volcano <- Fibro.il1.markers[c(2,5)]
regulated.Fibro <- Fibro.il1.volcano %>% filter(avg_log2FC > 1 | avg_log2FC < -1)  %>% row.names()
Fibro.ligands <- intersect(potential.ligands, regulated.Fibro)

p2 <- EnhancedVolcano(Fibro.il1.volcano,
                      lab = rownames(Fibro.il1.volcano),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      title = "Mtb infected WT vs IL-1R KO Fibroblast cells",
                      subtitle = "",
                      pCutoff = 0.05,
                      xlim = c(-5,5),
                      ylim = c(0,60),
                      selectLab = Fibro.ligands,
                      drawConnectors = TRUE,
                      arrowheads = FALSE,
                      labSize = 10)

Endo.il1.markers <- FindMarkers(simplified.stromal, assay = "RNA", ident.1 = "Endothelial_Mtb_WT", ident.2 = "Endothelial_Mtb_IL1KO")
Endo.il1.volcano <- Endo.il1.markers[c(2,5)]
regulated.Endo <- Endo.il1.volcano %>% filter(avg_log2FC > 1 | avg_log2FC < -1)  %>% row.names()
Endo.ligands <- intersect(potential.ligands, regulated.Endo)

p3 <- EnhancedVolcano(Endo.il1.volcano,
                      lab = rownames(Endo.il1.volcano),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      title = "Mtb infected WT vs IL-1R KO Endothelial cells",
                      subtitle = "",
                      pCutoff = 0.05,
                      xlim = c(-5,5),
                      ylim = c(0,60),
                      selectLab = Endo.ligands,
                      drawConnectors = TRUE,
                      arrowheads = FALSE,
                      labSize = 10)

p1 + p2 + p3
