library(Seurat)
library(tidyverse)
library(patchwork)
library(cowplot)
library(ggplot2)
library(nichenetr)
library(UCell)
library(EnhancedVolcano)
library(VennDiagram)
library(RColorBrewer)
library(ggpubr)

# Load the infected cell dataset
naive.rna <- Read10X(data.dir = "/Users/DmitriKotov/Documents/Sequencing Data/scRNA-seq 2-24-22/cDNA/RVDK005A/filtered_feature_bc_matrix", strip.suffix = TRUE)
naive.hto <- Read10X(data.dir = "/Users/DmitriKotov/Documents/Sequencing Data/scRNA-seq 2-24-22/ADT and HTO/RVDK005G HTO/umi_count", gene.column = 1)
naive.adt <- Read10X(data.dir = "/Users/DmitriKotov/Documents/Sequencing Data/scRNA-seq 2-24-22/ADT and HTO/RVDK005D ADT/umi_count", gene.column = 1)
infected.rna <- Read10X(data.dir = "/Users/DmitriKotov/Documents/Sequencing Data/scRNA-seq 2-24-22/cDNA/RVDK005B/filtered_feature_bc_matrix", strip.suffix = TRUE)
infected.hto <- Read10X(data.dir = "/Users/DmitriKotov/Documents/Sequencing Data/scRNA-seq 2-24-22/ADT and HTO/RVDK005H HTO/umi_count", gene.column = 1)
infected.adt <- Read10X(data.dir = "/Users/DmitriKotov/Documents/Sequencing Data/scRNA-seq 2-24-22/ADT and HTO/RVDK005E ADT/umi_count", gene.column = 1)
uninfected.rna <- Read10X(data.dir = "/Users/DmitriKotov/Documents/Sequencing Data/scRNA-seq 2-24-22/cDNA/RVDK005C/filtered_feature_bc_matrix", strip.suffix = TRUE)
uninfected.hto <- Read10X(data.dir = "/Users/DmitriKotov/Documents/Sequencing Data/scRNA-seq 2-24-22/ADT and HTO/RVDK005I HTO/umi_count", gene.column = 1)
uninfected.adt <- Read10X(data.dir = "/Users/DmitriKotov/Documents/Sequencing Data/scRNA-seq 2-24-22/ADT and HTO/RVDK005F ADT/umi_count", gene.column = 1)

#Keep cells with data for RNA, ADT, and HTO - infected
joint.infected = Reduce("intersect", list(colnames(infected.rna), colnames(infected.hto), colnames(infected.adt)))
infected.rna = infected.rna[, joint.infected]
infected.hto = as.matrix(infected.hto[-nrow(infected.hto), joint.infected])
infected.adt = as.matrix(infected.adt[-nrow(infected.adt), joint.infected])

#Keep cells with data for RNA, ADT, and HTO - naive
joint.naive = Reduce("intersect", list(colnames(naive.rna), colnames(naive.hto), colnames(naive.adt)))
naive.rna = naive.rna[, joint.naive]
naive.hto = as.matrix(naive.hto[-nrow(naive.hto), joint.naive])
naive.adt = as.matrix(naive.adt[-nrow(naive.adt), joint.naive])

#Keep cells with data for RNA, ADT, and HTO - uninfected
joint.uninfected = Reduce("intersect", list(colnames(uninfected.rna), colnames(uninfected.hto), colnames(uninfected.adt)))
uninfected.rna = uninfected.rna[, joint.uninfected]
uninfected.hto = as.matrix(uninfected.hto[-nrow(uninfected.hto), joint.uninfected])
uninfected.adt = as.matrix(uninfected.adt[-nrow(uninfected.adt), joint.uninfected])

# Confirm that the HTO have the correct names
rownames(infected.hto)
rownames(infected.hto) <- c("HTO_A", "HTO_B", "HTO_C", "HTO_D","HTO_E", "HTO_F")
rownames(naive.hto) <- c("HTO_A", "HTO_B", "HTO_C", "HTO_D","HTO_E", "HTO_F")
rownames(uninfected.hto) <- c("HTO_A", "HTO_B", "HTO_C", "HTO_D","HTO_E", "HTO_F")

# Confirm that the ADT have the correct names
rownames(infected.adt)
rownames(infected.adt) <- c("Ly6C", "CD44", "PD-L1", "SiglecF", "Csf1r", "Ly6G", "CD11b", "CD11c", "CD86", "MHCII", "CX3CR1", "CCR2", "CD62L", "CD45-2")
rownames(naive.adt) <- c("Ly6C", "CD44", "PD-L1", "SiglecF", "Csf1r", "Ly6G", "CD11b", "CD11c", "CD86", "MHCII", "CX3CR1", "CCR2", "CD62L", "CD45-2")
rownames(uninfected.adt) <- c("Ly6C", "CD44", "PD-L1", "SiglecF", "Csf1r", "Ly6G", "CD11b", "CD11c", "CD86", "MHCII", "CX3CR1", "CCR2", "CD62L", "CD45-2")

# Initialize the Seurat object with the raw (non-normalized data).
infected <- CreateSeuratObject(counts = infected.rna, project = "infected")
infected #check seurat object

# Add HTO data as a new assay independent from RNA
infected[["HTO"]] <- CreateAssayObject(counts = infected.hto)

# Add ADT data as a new assay independent from RNA
infected[["ADT"]] <- CreateAssayObject(counts = infected.adt)

#Create Seurat object for naive and uninfected
uninfected <- CreateSeuratObject(counts = uninfected.rna, project = "uninfected")
uninfected[["HTO"]] <- CreateAssayObject(counts = uninfected.hto)
uninfected[["ADT"]] <- CreateAssayObject(counts = uninfected.adt)
naive <- CreateSeuratObject(counts = naive.rna, project = "naive")
naive[["HTO"]] <- CreateAssayObject(counts = naive.hto)
naive[["ADT"]] <- CreateAssayObject(counts = naive.adt)

#check for dying cells based on mitochondrial percentage
infected[["percent.mt"]] <- PercentageFeatureSet(infected, pattern = "^mt-")
uninfected[["percent.mt"]] <- PercentageFeatureSet(uninfected, pattern = "^mt-")
naive[["percent.mt"]] <- PercentageFeatureSet(naive, pattern = "^mt-")

#standard filters are 
VlnPlot(infected, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# standard log-normalization
infected <- NormalizeData(infected)
uninfected <- NormalizeData(uninfected)
naive <- NormalizeData(naive)

# choose ~1k variable features
infected <- FindVariableFeatures(infected)
uninfected <- FindVariableFeatures(uninfected)
naive <- FindVariableFeatures(naive)

# standard scaling (no regression)
infected <- ScaleData(infected)
uninfected <- ScaleData(uninfected)
naive <- ScaleData(naive)

#HTO demultiplexing

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
infected <- NormalizeData(infected, assay = "HTO", normalization.method = "CLR")
uninfected <- NormalizeData(uninfected, assay = "HTO", normalization.method = "CLR")
naive <- NormalizeData(naive, assay = "HTO", normalization.method = "CLR")

#Run HTO demux
infected <- HTODemux(infected, assay = "HTO", positive.quantile = 0.99)
uninfected <- HTODemux(uninfected, assay = "HTO", positive.quantile = 0.99)
naive <- HTODemux(naive, assay = "HTO", positive.quantile = 0.99)

# Global classification results
table(infected$HTO_classification.global)

# Group cells based on the max HTO signal
Idents(infected) <- "HTO_maxID"

#Visualize multiplets
RidgePlot(infected, assay = "HTO", features = rownames(infected[["HTO"]])[1:2], ncol = 2)
FeatureScatter(infected, feature1 = "hto_HTO-A", feature2 = "hto_HTO-B")
Idents(infected) <- "HTO_classification.global"
VlnPlot(infected, features = "nCount_RNA", pt.size = 0.1, log = TRUE)

# First, we will remove negative cells from the object
infected.subset <- subset(infected, idents = "Negative", invert = TRUE)
uninfected.subset <- subset(uninfected, idents = "Negative", invert = TRUE)
naive.subset <- subset(naive, idents = "Negative", invert = TRUE)

# Calculate a tSNE embedding of the HTO data
DefaultAssay(infected.subset) <- "HTO"
infected.subset <- ScaleData(infected.subset, features = rownames(infected.subset), 
                             verbose = FALSE)
infected.subset <- RunPCA(infected.subset, features = rownames(infected.subset), approx = FALSE)
infected.subset <- RunTSNE(infected.subset, dims = 1:6, perplexity = 100, check_duplicates = FALSE) #check why there are duplicates
DimPlot(infected.subset)

# Extract the singlets
Idents(infected) <- "HTO_classification.global"
infected.singlet <- subset(infected, idents = "Singlet")
Idents(uninfected) <- "HTO_classification.global"
uninfected.singlet <- subset(uninfected, idents = "Singlet")
Idents(naive) <- "HTO_classification.global"
naive.singlet <- subset(naive, idents = "Singlet")

#Filter to exclude cells that have unique feature counts more than 4500 (multiplets) and less than 200 (empty drops) and removing cells that have >5% mitochondrial counts (dead/dying cells)
#Macrophages were heavily removed when removing cells with feature counts over 2500, hence why I raised it to 4500.
infected.singlet <- subset(infected.singlet, subset = nFeature_RNA < 4500 & nFeature_RNA > 200 & percent.mt < 5)
naive.singlet <- subset(naive.singlet, subset = nFeature_RNA < 4500 & nFeature_RNA > 200 & percent.mt < 5)
uninfected.singlet <- subset(uninfected.singlet, subset = nFeature_RNA < 4500 & nFeature_RNA > 200 & percent.mt < 5)

#for dividing the datasets into Irg1 or iNOS
infected.Irg1 <- subset(infected.singlet, subset = HTO_maxID == "HTO-A" | HTO_maxID == "HTO-B" | HTO_maxID == "HTO-C")
infected.iNOS <- subset(infected.singlet, subset = HTO_maxID == "HTO-D" | HTO_maxID == "HTO-E" | HTO_maxID == "HTO-F")
uninfected.Irg1 <- subset(uninfected.singlet, subset = HTO_maxID == "HTO-A" | HTO_maxID == "HTO-B" | HTO_maxID == "HTO-C")
uninfected.iNOS <- subset(uninfected.singlet, subset = HTO_maxID == "HTO-D" | HTO_maxID == "HTO-E" | HTO_maxID == "HTO-F")
naive.Irg1 <- subset(naive.singlet, subset = HTO_maxID == "HTO-A" | HTO_maxID == "HTO-B" | HTO_maxID == "HTO-C")
naive.iNOS <- subset(naive.singlet, subset = HTO_maxID == "HTO-D" | HTO_maxID == "HTO-E" | HTO_maxID == "HTO-F")

#Add metadata about genotype and infection status (Mtb+ = cell was infected; Mtb- = cell was not infected but came from an infected mouse)
infected.Irg1[["state"]] <- "Mtb+"
infected.iNOS[["state"]] <- "Mtb+"
infected.Irg1[["cell_type"]] <- "Irg1"
infected.iNOS[["cell_type"]] <- "iNOS"
uninfected.Irg1[["state"]] <- "Mtb-"
uninfected.iNOS[["state"]] <- "Mtb-"
uninfected.Irg1[["cell_type"]] <- "Irg1"
uninfected.iNOS[["cell_type"]] <- "iNOS"
naive.Irg1[["state"]] <- "naive"
naive.iNOS[["state"]] <- "naive"
naive.Irg1[["cell_type"]] <- "Irg1"
naive.iNOS[["cell_type"]] <- "iNOS"

#continue with integrated analysis here: #what is appropriate dims for integration?
object.list <-c(infected.Irg1, infected.iNOS, uninfected.Irg1, uninfected.iNOS, naive.Irg1, naive.iNOS)
immune.anchors <- FindIntegrationAnchors(object.list = object.list, dims = 1:30) #resulted in cells with duplicat naming. Look into subsetting by hash.ID instead of HTO_maxID
irg1.inos <- IntegrateData(anchorset = immune.anchors, dims = 1:30)
DefaultAssay(irg1.inos) <- "integrated"

# Run the standard workflow for visualization and clustering
irg1.inos <- ScaleData(irg1.inos, verbose = FALSE)
irg1.inos <- RunPCA(irg1.inos, npcs = 30, verbose = FALSE)

# t-SNE and Clustering
irg1.inos <- RunUMAP(irg1.inos, reduction = "pca", dims = 1:30)
irg1.inos <- FindNeighbors(irg1.inos, reduction = "pca", dims = 1:30)
irg1.inos <- FindClusters(irg1.inos, resolution = 0.8)

#Visualize integrated dataset
p1 <- DimPlot(irg1.inos, reduction = "umap", group.by = "state")
p2 <- DimPlot(irg1.inos, reduction = "umap", label = TRUE)
p3 <- DimPlot(irg1.inos, reduction = "umap", group.by = "cell_type")
plot_grid(p1, p2, p3)

DimPlot(irg1.inos, reduction = "umap", split.by = "cell_type")
DimPlot(irg1.inos, reduction = "umap", split.by = "state")

#Improve clustering through Weighted nearest neighbor analysis

# Process ADT data and  set a dimensional reduction name to avoid overwriting the RNA PCA
DefaultAssay(irg1.inos) <- "ADT"
VariableFeatures(irg1.inos) <- rownames(irg1.inos[["ADT"]])
irg1.inos <- NormalizeData(irg1.inos, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

irg1.inos <- FindMultiModalNeighbors(
  irg1.inos, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:10), modality.weight.name = "RNA.weight"
)

irg1.inos <- RunUMAP(irg1.inos, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
irg1.inos <- FindClusters(irg1.inos, graph.name = "wsnn", algorithm = 3, resolution = 1.5, verbose = FALSE)

plot <- DimPlot(irg1.inos, reduction = "wnn.umap", label = TRUE)

saveRDS(irg1.inos, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/irg1 inos")
irg1.inos <- readRDS(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/irg1 inos")

DefaultAssay(irg1.inos) <- "integrated"
FeaturePlot(irg1.inos, reduction = "wnn.umap", features = c("adt_Ly6G", "adt_Ly6C", "adt_CX3CR1", "adt_CCR2","adt_Csf1r", "adt_MHCII", "adt_SiglecF", "adt_CD11b", "adt_CD45-2",
                                          "Ly6g", "Ly6c1", "Cx3cr1","Csf1r","H2-Ab1","Siglecf"), min.cutoff = "q05", max.cutoff = "q95", ncol = 5)

#Try integrating iNOS and Irg1 datasets with the B6 and Sp140 datasets
Idents(irg1.inos) <- "cell_type"
Irg1 <-subset(irg1.inos, idents = "Irg1")
iNOS <-subset(irg1.inos, idents = "iNOS")
Irg1.list <- SplitObject(Irg1, split.by = "state")
names(Irg1.list) <- c("Irg1.Mtb+","Irg1.Mtb-","Irg1.naive")
iNOS.list <- SplitObject(iNOS, split.by = "state")
names(iNOS.list) <- c("iNOS.Mtb+","iNOS.Mtb-","iNOS.naive")

immune.combined <- readRDS(file = "/Users/DmitriKotov/Box Sync/Coding stuff/B6 vs Sp140 scRNA-seq 10-9-20/immunecombined3")
Idents(immune.combined) <- "cell_type"
b6 <- subset(immune.combined, idents = "B6")
sp140 <- subset(immune.combined, idents = "Sp140")
b6.list <- SplitObject(b6, split.by = "state")
names(b6.list) <- c("b6.Mtb+","b6.Mtb-","b6.naive")
sp140.list <- SplitObject(sp140, split.by = "state")
names(sp140.list) <- c("sp140.Mtb+","sp140.Mtb-","sp140.naive")
super.combined.list <- c(b6.list, sp140.list, Irg1.list, iNOS.list)

super.combined.list <- lapply(X = super.combined.list, FUN = function(x) {
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = super.combined.list)

super.combined.list <- lapply(X = super.combined.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
anchors <- FindIntegrationAnchors(object.list = super.combined.list, reference = c(1:6), reduction = "rpca",
                                  dims = 1:50)
super.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
super.integrated <- ScaleData(super.integrated, verbose = FALSE)
super.integrated <- RunPCA(super.integrated, verbose = FALSE)
super.integrated <- RunUMAP(super.integrated, dims = 1:50)
super.integrated[["cell_type.state"]] <- paste(super.integrated$cell_type, super.integrated$state, sep = "_")

super.integrated <- FindNeighbors(super.integrated, reduction = "pca", dims = 1:30)
super.integrated <- FindClusters(super.integrated, graph.name= "integrated_snn", resolution = 1.5)

DimPlot(super.integrated, group.by = "cell_type.state")
p1 <- DimPlot(super.integrated, reduction = "umap", label = TRUE)
Idents(super.integrated) <- "celltype"
old.labels <- subset(super.integrated, idents = c("Basophil", "B cell", "DC", "AM", "Cycling AM", "Mono", "Int Mono", "Inflam Mono", "CD16-2 Mono", "Lyve1 IM", "IM", "IFN IM", "Neut", "SigF Neut", "Mmp8 Neut", "Stfa Neut", "mNeut", "Act Neut", "Aged Neut", "IFN Aged Neut", "Nos2 Neut", "Itpr2 Nreg", "Mif Nreg"))
p2 <- DimPlot(old.labels, label = TRUE)
p1 + p2

#Figuring out neutrophil subsets
FeaturePlot(super.integrated, features = c("Mmp8","Cd274", "Cxcr2","Cxcr4","Isg15","Cd63","Camp","Sell")) + p1 + p2
#Based on the Neutrotime paper from Immgen (https://www.nature.com/articles/s41467-021-22973-9#Sec46), Camp and Ltf are early genes followed by Mmp8. Mature Neut's lacked Mmp8 expression but still had CXCR2.
# Aged Neuts had lower Sell and Cxcr2 while increasing in Cxcr4 expression. ISG Young and ISG Old had the strongest ISG signature, but it was still present in the ISG Activated and ISG Mature cells.
# CD63 were a unique population, with high CD63 expression, no CXCR2, no CD62L, and low CXCR4. Nos2 Neut had various activation signatures but lacked a strong ISG signature.

#Find markers for each cluster
DefaultAssay(super.integrated) <- "RNA"
Idents(super.integrated) <- "seurat_clusters"
markers <- FindAllMarkers(super.integrated, only.pos = TRUE, min.pct = 0.25, verbose = TRUE)

#Label clusters and save labeling in metadata
super.integrated <- RenameIdents(super.integrated, `0` = "ISG Mature Neut", `1` = "ISG Activated Neut", `2` = "Aged Neut", 
                                `3` = "Mature Neut", `4` = "ISG Old Neut", `5` = "ISG Young Neut", `6` = "CD63 Neut", `7` = "Mono", `8` = "Aged Neut", `9` = "AM", 
                                `10` = "Mmp8 Neut", `11` = "IM", `12` = "ISG IM", `13` = "Nos2 Neut",`14` = "Activated IM",`15` = "CD16-2 Mono",
                                `16` = "Activated AM",`17` = "Nos2 Neut",`18` = "Camp Mmp8 Neut",`19` = "Stromal?",`20` = "gCap",`21` = "Cycling AM",`22` = "DC",
                                `23` = "Lipofibro?",`24` = "B cell",`25` = "fibroblast",`26` = "Basophil",`27`="aCap")
super.integrated$celltype <- Idents(super.integrated)

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

m.s.genes <-convert_human_to_mouse_symbols(s.genes)
m.g2m.genes <- convert_human_to_mouse_symbols(g2m.genes)

super.integrated <- CellCycleScoring(super.integrated, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
# view cell cycle scores and phase assignments
head(super.integrated[[]])
#visualize different cell cycle phases
DimPlot(super.integrated)
#switch ident back to clusters
Idents(super.integrated) <- "celltype"

#Remove Lipofibro because they don't belong to any single cluster and don't appear to be proliferating cells or doublets.
Lipofibro <- WhichCells(super.integrated, idents = "Lipofibro?")
DimPlot(super.integrated, label = T, cells.highlight = Lipofibro)
super.integrated <- subset(super.integrated, subset = celltype != "Lipofibro?")

#Save the dataset
saveRDS(super.integrated, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/b6 sp140 irg1 inos")
super.integrated <- readRDS(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/b6 sp140 irg1 inos")

#Try to perform weighted nearest neighbor with combined dataset - doesn't work because of differences in ADT data
DefaultAssay(super.integrated) <- "ADT"
VariableFeatures(super.integrated) <- rownames(super.integrated[["ADT"]])
super.integrated <- NormalizeData(super.integrated, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

super.integrated <- FindMultiModalNeighbors(
  super.integrated, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:10), modality.weight.name = "RNA.weight"
)

super.integrated <- RunUMAP(super.integrated, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
super.integrated <- FindClusters(super.integrated, graph.name = "wsnn", algorithm = 3, resolution = 1.5, verbose = FALSE)

plot <- DimPlot(super.integrated, reduction = "wnn.umap", label = TRUE)

#Plot naive IL1Rn expression in neutrophils for all four genotypes
neutrophil <- super.integrated
neutrophil <- RenameIdents(neutrophil, 'ISG Mature Neut' = "Neut",'ISG Activated Neut'="Neut",'Aged Neut'="Neut",'Mature Neut'="Neut",'ISG Old Neut'="Neut",
                           'ISG Young Neut'="Neut",'CD63 Neut'="Neut",'Mmp8 Neut'="Neut","Nos2 Neut"="Neut",'Camp Mmp8 Neut'="Neut")
neutrophil$celltype <- Idents(neutrophil)
Idents(neutrophil) <- "state"
neutrophil.naive <- subset(neutrophil, idents = "naive")
Idents(neutrophil.naive) <- "celltype"
VlnPlot(neutrophil.naive, features = "Il1rn",idents = c("Neut","AM"), split.by = "cell_type")

#Plot IFNb expressing cells.
infected.mice <- subset(super.integrated, idents = c("Mtb+","Mtb-"))
FeaturePlot(infected.mice, features = "Ifnb1", split.by = "cell_type", order = T)
FeaturePlot(super.integrated, features = "Ifnb1", split.by = "state", order = T)
ifnb.pos <- subset(super.integrated, Ifnb1 > 0)

#Add statistics to the Ifnb1+ expression level plot
vp_case1 <- function(gene_signature, test_sign, y_max){
  plot_case1 <- function(signature){
    VlnPlot(ifnb.pos, features = signature,
            pt.size = 0.1, 
            group.by = "cell_type", 
            y.max = y_max, # add the y-axis maximum value - otherwise p-value hidden
    ) + stat_compare_means(comparisons = test_sign, label = "p.signif")
  }
  purrr::map(gene_signature, plot_case1) %>% cowplot::plot_grid(plotlist = .)
}
gene_sig <- c("Ifnb1")
comparisons <- list(c("B6","Sp140"))
vp_case1(gene_signature = gene_sig, test_sign = comparisons, y_max = 3)

#Reorder celltype identities.
levels(super.integrated) <- c("gCap","aCap","fibroblast","Stromal?","DC","B cell","Basophil","Camp Mmp8 Neut","Mmp8 Neut","Mature Neut",
                              "Aged Neut","ISG Mature Neut","ISG Activated Neut","ISG Young Neut","ISG Old Neut","Nos2 Neut","CD63 Neut",
                              "AM","Activated AM","Cycling AM","CD16-2 Mono","Mono","IM","Activated IM","ISG IM")

super.integrated$celltype <- factor(super.integrated$celltype, levels = c("gCap","aCap","fibroblast","Stromal?","DC","B cell","Basophil","Camp Mmp8 Neut","Mmp8 Neut","Mature Neut",
                                                                          "Aged Neut","ISG Mature Neut","ISG Activated Neut","ISG Young Neut","ISG Old Neut","Nos2 Neut","CD63 Neut",
                                                                          "AM","Activated AM","Cycling AM","CD16-2 Mono","Mono","IM","Activated IM","ISG IM"))
DimPlot(super.integrated, split.by = "state")

#Full example genes for all subsets
DotPlot(super.integrated, features = c("Ptprc","adt_Ly6G","Cxcr2","Cxcr4","Camp","Mmp8","Sell","adt_SiglecF","Mertk","Fcgr1","Ccr2","adt_CX3CR1","Fcgr4","Cd19","Itgax","adt_MHCII","adt_CD11b","Nos2","Cd63","Isg15","Ch25h","Fcer1a"))
#Subset of genes for general classification
DotPlot(super.integrated, features = c("Ptprc","adt_Ly6G","adt_SiglecF","Mertk","Fcgr1","Ccr2","Fcgr4","Cd19","Fcer1a"))

#Demonstrating markers for subsetting Neutorphils
Neut.plot <- subset(super.integrated, idents = c("Camp Mmp8 Neut","Mmp8 Neut","Mature Neut",
                                                 "Aged Neut","ISG Mature Neut","ISG Activated Neut","ISG Young Neut","ISG Old Neut","Nos2 Neut","CD63 Neut"))
DotPlot(Neut.plot, features = c("adt_CD11b","adt_Ly6G","Ly6g","Cxcr2","Cxcr4","Camp","Mmp8","Sell","Cd63","Isg15","Nos2","Hypoxia"))

#Markers for differentiating between the various monocyte/macrophage populations
Mac.plot <- subset(super.integrated, idents = c("AM","Activated AM","Cycling AM","CD16-2 Mono","Mono","IM","Activated IM","ISG IM"))
DotPlot(Mac.plot, features = c("Csf1r","adt_CD11b","adt_SiglecF","Fcgr1","Mertk","adt_CD11c","Ccr2","adt_CX3CR1",
                               "Fcgr4","Ly6c2","Sell","adt_MHCII","Nos2","Cd63","Isg15","adt_PD-L1","adt_CD86","Mki67"))

#Add Mouse Identity Data for infected mice
super.integrated$cell_type.state.hash <- paste(super.integrated$cell_type.state, super.integrated$hash.ID, sep= "_")
combos <- unique(super.integrated$cell_type.state.hash)
super.cells <- lapply(X = combos, FUN = function(x){
  Idents(super.integrated) <- "cell_type.state.hash"
  WhichCells(super.integrated, idents = x)
})
names(super.cells) <- combos

super.integrated$mouse_id <- ifelse(colnames(super.integrated) %in% c(super.cells$`B6_Mtb+_HTO-E`,super.cells$`B6_Mtb-_HTO-E`,
                                                                      super.cells$`Sp140_Mtb+_HTO-A`,super.cells$`Sp140_Mtb-_HTO-A`,
                                                                      super.cells$`Irg1_Mtb+_HTO-A`,super.cells$`Irg1_Mtb-_HTO-A`,
                                                                      super.cells$`iNOS_Mtb+_HTO-E`,super.cells$`iNOS_Mtb-_HTO-E`), 1,
                                    ifelse(colnames(super.integrated) %in% c(super.cells$`B6_Mtb+_HTO-F`,super.cells$`B6_Mtb-_HTO-F`,
                                                                             super.cells$`Sp140_Mtb+_HTO-B`,super.cells$`Sp140_Mtb-_HTO-B`,
                                                                             super.cells$`Irg1_Mtb+_HTO-B`,super.cells$`Irg1_Mtb-_HTO-B`,
                                                                             super.cells$`iNOS_Mtb+_HTO-F`,super.cells$`iNOS_Mtb-_HTO-F`), 2,
                                           ifelse(colnames(super.integrated) %in% c(super.cells$`B6_Mtb+_HTO-D`,super.cells$`B6_Mtb-_HTO-D`,
                                                                                    super.cells$`Sp140_Mtb+_HTO-C`,super.cells$`Sp140_Mtb-_HTO-C`,
                                                                                    super.cells$`Irg1_Mtb+_HTO-C`,super.cells$`Irg1_Mtb-_HTO-C`,
                                                                                    super.cells$`iNOS_Mtb+_HTO-D`,super.cells$`iNOS_Mtb-_HTO-D`), 3, 
                                                  ifelse(colnames(super.integrated) %in% c(super.cells$`B6_naive_HTO-A`, super.cells$`Sp140_naive_HTO-D`,
                                                                                          super.cells$`iNOS_naive_HTO-D`, super.cells$`Irg1_naive_HTO-A`), 4,
                                                         ifelse(colnames(super.integrated) %in% c(super.cells$`B6_naive_HTO-B`, super.cells$`Sp140_naive_HTO-E`,
                                                                                                 super.cells$`iNOS_naive_HTO-E`, super.cells$`Irg1_naive_HTO-B`), 5,
                                                                ifelse(colnames(super.integrated) %in% c(super.cells$`B6_naive_HTO-B`, super.cells$`Sp140_naive_HTO-E`,
                                                                                                        super.cells$`iNOS_naive_HTO-F`, super.cells$`Irg1_naive_HTO-C`), 6,NA))))))

super.integrated$cell_type.mouse_id <- paste(super.integrated$cell_type,super.integrated$mouse_id, sep = "_")

#Neutrophil B6 vs Sp140 volcano plot, IM B6 vs Sp140 volcano plot - both from infected mice
super.simple.neut <- subset(super.simple, idents = "Neutrophil")
Idents(super.simple.neut) <- "cell_type"
Neut.B6.Sp140.markers <- FindMarkers(super.simple.neut, ident.1 = "B6", ident.2 = "Sp140")
Neut.B6.Sp140.volcano <- Neut.B6.Sp140.markers[c(2,5)]
EnhancedVolcano(Neut.B6.Sp140.volcano,
                lab = rownames(Neut.B6.Sp140.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "B6 vs Sp140 Neutrophils",
                subtitle = "",
                raster = T)

super.simple.IM <- subset(super.simple, idents = "IM")
Idents(super.simple.IM) <- "cell_type"
IM.B6.Sp140.markers <- FindMarkers(super.simple.IM, ident.1 = "B6", ident.2 = "Sp140")
IM.B6.Sp140.volcano <- IM.B6.Sp140.markers[c(2,5)]
EnhancedVolcano(IM.B6.Sp140.volcano,
                lab = rownames(IM.B6.Sp140.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "B6 vs Sp140 IM",
                subtitle = "",
                raster = T)

#Check whether any genes are differentially expressed in naive mice
Idents(super.integrated) <- "state"
naive.integrated <- subset(super.integrated,idents = "naive")
Idents(naive.integrated) <- "celltype"
naive.integrated <- RenameIdents(naive.integrated, 'ISG Mature Neut' = "Neutrophil", 'ISG Activated Neut' = "Neutrophil", 'Aged Neut' = "Neutrophil", 'Mature Neut' = "Neutrophil",
                                  'ISG Old Neut' = "Neutrophil", 'ISG Young Neut' = "Neutrophil", 'CD63 Neut' = "Neutrophil", 'Mmp8 Neut' = "Neutrophil",'Nos2 Neut' = "Neutrophil",
                                  'Camp Mmp8 Neut' = "Neutrophil", 'CD16-2 Mono' = "Mono", 'ISG IM' = "IM", 'Activated IM' = "IM", 'Activated AM' = "AM")
DimPlot(naive.integrated, split.by = "cell_type")

#Make Venn Diagrams for B6 vs Sp140; B6 vs Nos2; B6 vs Acod1 for IM, AM, Mono, and Neutrophils for naive mice
myCol <- brewer.pal(3, "Pastel2")

#Neutrophils
naive.neutrophil <- subset(naive.integrated, idents = "Neutrophil")
Idents(naive.neutrophil) <- "cell_type"
naive.neutrophils.sp140 <- FindMarkers(naive.neutrophil, ident.1 = "B6",ident.2 = "Sp140")
naive.neutrophils.nos2 <- FindMarkers(naive.neutrophil, ident.1 = "B6",ident.2 = "iNOS")
naive.neutrophils.acod1 <- FindMarkers(naive.neutrophil, ident.1 = "B6",ident.2 = "Irg1")
naive.neutrophils.sp140 <- filter(naive.neutrophils.sp140, avg_log2FC > 1 | avg_log2FC < -1 & p_val_adj < 0.05)
naive.neutrophils.nos2 <- filter(naive.neutrophils.nos2, avg_log2FC > 1 | avg_log2FC < -1 & p_val_adj < 0.05)
naive.neutrophils.acod1 <- filter(naive.neutrophils.acod1, avg_log2FC > 1 | avg_log2FC < -1 & p_val_adj < 0.05)
venn.diagram(x = list(rownames(naive.neutrophils.sp140),rownames(naive.neutrophils.nos2),rownames(naive.neutrophils.acod1)),
            category.names = c("Sp140", "Nos2","Acod1"), filename = "naive neutrophils.tiff", resolution = 600, 
            col=myCol,
            fill = c(alpha(myCol[1],0.3), alpha(myCol[2],0.3), alpha(myCol[3],0.3)),
            euler.d = FALSE,
            scaled = FALSE,
            cat.cex = 0)

#AM
naive.AM <- subset(naive.integrated, idents = "AM")
Idents(naive.AM) <- "cell_type"
naive.AM.sp140 <- FindMarkers(naive.AM, ident.1 = "B6",ident.2 = "Sp140")
naive.AM.nos2 <- FindMarkers(naive.AM, ident.1 = "B6",ident.2 = "iNOS")
naive.AM.acod1 <- FindMarkers(naive.AM, ident.1 = "B6",ident.2 = "Irg1")
naive.AM.sp140 <- filter(naive.AM.sp140, avg_log2FC > 1 | avg_log2FC < -1 & p_val_adj < 0.05)
naive.AM.nos2 <- filter(naive.AM.nos2, avg_log2FC > 1 | avg_log2FC < -1 & p_val_adj < 0.05)
naive.AM.acod1 <- filter(naive.AM.acod1, avg_log2FC > 1 | avg_log2FC < -1 & p_val_adj < 0.05)
venn.diagram(x = list(rownames(naive.AM.sp140),rownames(naive.AM.nos2),rownames(naive.AM.acod1)),
             category.names = c("Sp140", "Nos2","Acod1"), filename = "naive AM.tiff", resolution = 600, 
             col=myCol,
             fill = c(alpha(myCol[1],0.3), alpha(myCol[2],0.3), alpha(myCol[3],0.3)),
             euler.d = FALSE,
             scaled = FALSE,
             cat.cex = 0)

#IM
naive.IM <- subset(naive.integrated, idents = "IM")
Idents(naive.IM) <- "cell_type"
naive.IM.sp140 <- FindMarkers(naive.IM, ident.1 = "B6",ident.2 = "Sp140")
naive.IM.nos2 <- FindMarkers(naive.IM, ident.1 = "B6",ident.2 = "iNOS")
naive.IM.acod1 <- FindMarkers(naive.IM, ident.1 = "B6",ident.2 = "Irg1")
naive.IM.sp140 <- filter(naive.IM.sp140, avg_log2FC > 1 | avg_log2FC < -1 & p_val_adj < 0.05)
naive.IM.nos2 <- filter(naive.IM.nos2, avg_log2FC > 1 | avg_log2FC < -1 & p_val_adj < 0.05)
naive.IM.acod1 <- filter(naive.IM.acod1, avg_log2FC > 1 | avg_log2FC < -1 & p_val_adj < 0.05)
venn.diagram(x = list(rownames(naive.IM.sp140),rownames(naive.IM.nos2),rownames(naive.IM.acod1)),
             category.names = c("Sp140", "Nos2","Acod1"), filename = "naive IM.tiff", resolution = 600, 
             col=myCol,
             fill = c(alpha(myCol[1],0.3), alpha(myCol[2],0.3), alpha(myCol[3],0.3)),
             euler.d = FALSE,
             scaled = FALSE,
             cat.cex = 0)

#Monocytes
naive.mono <- subset(naive.integrated, idents = "Mono")
Idents(naive.mono) <- "cell_type"
naive.mono.sp140 <- FindMarkers(naive.mono, ident.1 = "B6",ident.2 = "Sp140")
naive.mono.nos2 <- FindMarkers(naive.mono, ident.1 = "B6",ident.2 = "iNOS")
naive.mono.acod1 <- FindMarkers(naive.mono, ident.1 = "B6",ident.2 = "Irg1")
naive.mono.sp140 <- filter(naive.mono.sp140, avg_log2FC > 1 | avg_log2FC < -1 & p_val_adj < 0.05)
naive.mono.nos2 <- filter(naive.mono.nos2, avg_log2FC > 1 | avg_log2FC < -1 & p_val_adj < 0.05)
naive.mono.acod1 <- filter(naive.mono.acod1, avg_log2FC > 1 | avg_log2FC < -1 & p_val_adj < 0.05)
venn.diagram(x = list(rownames(naive.mono.sp140),rownames(naive.mono.nos2),rownames(naive.mono.acod1)),
             category.names = c("Sp140", "Nos2","Acod1"), filename = "naive monocytes.tiff", resolution = 600, 
             col=myCol,
             fill = c(alpha(myCol[1],0.3), alpha(myCol[2],0.3), alpha(myCol[3],0.3)),
             euler.d = FALSE,
             scaled = FALSE,
             cat.cex = 0)

#Define signatures for B6 infected versus naive
super.integrated <- readRDS(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/b6 sp140 irg1 inos v2")
super.integrated$infected <- ifelse(super.integrated$state == "Mtb+" | super.integrated$state == "Mtb-", "infected","naive")
super.integrated$cell_type.infected <- paste(super.integrated$cell_type, super.integrated$infected, sep = "_")
Idents(super.integrated) <- "cell_type.infected"
B6 <- FindMarkers(super.integrated, ident.1 = "B6_infected", ident.2 = "B6_naive")
Sp140 <- FindMarkers(super.integrated, ident.1 = "Sp140_infected", ident.2 = "Sp140_naive")
Irg1 <- FindMarkers(super.integrated, ident.1 = "Irg1_infected", ident.2 = "Irg1_naive")
iNOS <- FindMarkers(super.integrated, ident.1 = "iNOS_infected", ident.2 = "iNOS_naive")
Sp140.b6 <- FindMarkers(super.integrated, ident.1 = "Sp140_infected", ident.2 = "B6_infected")
Irg1.b6 <- FindMarkers(super.integrated, ident.1 = "Irg1_infected", ident.2 = "B6_infected")
iNOS.b6 <- FindMarkers(super.integrated, ident.1 = "iNOS_infected", ident.2 = "B6_infected")
Idents(super.integrated) <- "cell_type.state"
b6.mtb.pos <- FindMarkers(super.integrated, ident.1 = "B6_Mtb+", ident.2 = "B6_Mtb-")
B6$gene <- rownames(B6)
Sp140$gene <- rownames(Sp140)
Irg1$gene <- rownames(Irg1)
iNOS$gene <- rownames(iNOS)
Sp140.b6$gene <- rownames(Sp140.b6)
Irg1.b6$gene <- rownames(Irg1.b6)
iNOS.b6$gene <- rownames(iNOS.b6)
b6.mtb.pos$gene <- rownames(b6.mtb.pos)

Idents(super.integrated) <- "cell_type.state"
b6.mtb.pos.neg <- FindMarkers(super.integrated, ident.1 = "B6_Mtb+", ident.2 = "B6_naive")
b6.bystander.neg <- FindMarkers(super.integrated, ident.1 = "B6_Mtb-", ident.2 = "B6_naive")
b6.mtb.pos.neg$gene <- rownames(b6.mtb.pos.neg)
b6.bystander.neg$gene <- rownames(b6.bystander.neg)
Sp140.mtb.pos <- FindMarkers(super.integrated, ident.1 = "Sp140_Mtb+", ident.2 = "Sp140_naive")
Irg1.mtb.pos <- FindMarkers(super.integrated, ident.1 = "Irg1_Mtb+", ident.2 = "Irg1_naive")
iNOS.mtb.pos <- FindMarkers(super.integrated, ident.1 = "iNOS_Mtb+", ident.2 = "iNOS_naive")
Sp140.mtb.pos$gene <- rownames(Sp140.mtb.pos)
Irg1.mtb.pos$gene <- rownames(Irg1.mtb.pos)
iNOS.mtb.pos$gene <- rownames(iNOS.mtb.pos)

#Filtering steps
Sp140.infected <- Sp140.b6 %>% filter(avg_log2FC > 1 | avg_log2FC < -1) %>% inner_join(filter(Sp140, avg_log2FC > 1), by = "gene")
Irg1.infected <- Irg1.b6 %>% filter(avg_log2FC > 1 | avg_log2FC < -1) %>% inner_join(filter(Irg1, avg_log2FC > 1), by = "gene")
iNOS.infected <- iNOS.b6 %>% filter(avg_log2FC > 1 | avg_log2FC < -1) %>% inner_join(filter(iNOS, avg_log2FC > 1), by = "gene")
conserved <- Sp140.infected %>% inner_join(Irg1.infected, by = "gene") %>% inner_join(iNOS.infected, by = "gene") %>% mutate(change = ifelse(avg_log2FC.x > 1, "up", "down")) %>% select(gene, change) %>% mutate(gene = convert_mouse_to_human_symbols(gene)) %>% drop_na()
Sp140.unique <- Sp140.infected %>% anti_join(Irg1.infected, by = "gene") %>% anti_join(iNOS.infected, by = "gene") %>% mutate(change = ifelse(avg_log2FC.x > 1, "up", "down")) %>% select(gene, change) %>% mutate(gene = convert_mouse_to_human_symbols(gene)) %>% drop_na()
B6.general <- B6 %>% filter(avg_log2FC > 1 | avg_log2FC < -1) %>% mutate(change = ifelse(avg_log2FC > 1, "up", "down")) %>% mutate(gene = convert_mouse_to_human_symbols(gene)) %>% drop_na()
B6.Mtb.Pos <- b6.mtb.pos %>% filter(avg_log2FC > 1 | avg_log2FC < -1) %>% inner_join(filter(B6, avg_log2FC > 1 | avg_log2FC < -1), by = "gene") %>% mutate(change = ifelse(avg_log2FC.x > 1, "up", "down")) %>% select(gene, change) %>% mutate(gene = convert_mouse_to_human_symbols(gene)) %>% drop_na()

#Filtering Steps Round 2 - incorporate feedback from Sanjana
Sp140.general <- Sp140 %>% filter(avg_log2FC > 1 | avg_log2FC < -1) %>% mutate(change = ifelse(avg_log2FC > 1, "up", "down")) %>% mutate(gene = convert_mouse_to_human_symbols(gene)) %>% drop_na()
Irg1.general <- Irg1 %>% filter(avg_log2FC > 1 | avg_log2FC < -1) %>% mutate(change = ifelse(avg_log2FC > 1, "up", "down")) %>% mutate(gene = convert_mouse_to_human_symbols(gene)) %>% drop_na()
iNOS.general <- iNOS %>% filter(avg_log2FC > 1 | avg_log2FC < -1) %>% mutate(change = ifelse(avg_log2FC > 1, "up", "down")) %>% mutate(gene = convert_mouse_to_human_symbols(gene)) %>% drop_na()
b6.mtb.pos.neg.1 <- b6.mtb.pos.neg %>% filter(avg_log2FC > 1 | avg_log2FC < -1) %>% mutate(change = ifelse(avg_log2FC > 1, "up", "down")) %>% mutate(gene = convert_mouse_to_human_symbols(gene)) %>% drop_na()
b6.bystander.neg <- b6.bystander.neg %>% filter(avg_log2FC > 1 | avg_log2FC < -1) %>% mutate(change = ifelse(avg_log2FC > 1, "up", "down")) %>% mutate(gene = convert_mouse_to_human_symbols(gene)) %>% drop_na()

#Filtering Steps Round 3 - more stringent to reduce the type I IFN signature from iNOS and Irg1
B6.strict <- B6 %>% filter(avg_log2FC > 2 | avg_log2FC < -2) %>% mutate(change = ifelse(avg_log2FC > 1, "up", "down")) %>% mutate(gene = convert_mouse_to_human_symbols(gene)) %>% drop_na()
Sp140.strict <- Sp140 %>% filter(avg_log2FC > 2 | avg_log2FC < -2) %>% mutate(change = ifelse(avg_log2FC > 1, "up", "down")) %>% mutate(gene = convert_mouse_to_human_symbols(gene)) %>% drop_na()
Irg1.strict <- Irg1 %>% filter(avg_log2FC > 2 | avg_log2FC < -2) %>% mutate(change = ifelse(avg_log2FC > 1, "up", "down")) %>% mutate(gene = convert_mouse_to_human_symbols(gene)) %>% drop_na()
iNOS.strict <- iNOS %>% filter(avg_log2FC > 2 | avg_log2FC < -2) %>% mutate(change = ifelse(avg_log2FC > 1, "up", "down")) %>% mutate(gene = convert_mouse_to_human_symbols(gene)) %>% drop_na()
b6.mtb.pos.neg.strict <- b6.mtb.pos.neg %>% filter(avg_log2FC > 2 | avg_log2FC < -2) %>% mutate(change = ifelse(avg_log2FC > 1, "up", "down")) %>% mutate(gene = convert_mouse_to_human_symbols(gene)) %>% drop_na()
b6.bystander.neg.strict <- b6.bystander.neg %>% filter(avg_log2FC > 2 | avg_log2FC < -2) %>% mutate(change = ifelse(avg_log2FC > 1, "up", "down")) %>% mutate(gene = convert_mouse_to_human_symbols(gene)) %>% drop_na()
Sp140.mtb.pos.strict <- Sp140.mtb.pos %>% filter(avg_log2FC > 2 | avg_log2FC < -2) %>% mutate(change = ifelse(avg_log2FC > 1, "up", "down")) %>% mutate(gene = convert_mouse_to_human_symbols(gene)) %>% drop_na()
Irg1.mtb.pos.strict <- Irg1.mtb.pos %>% filter(avg_log2FC > 2 | avg_log2FC < -2) %>% mutate(change = ifelse(avg_log2FC > 1, "up", "down")) %>% mutate(gene = convert_mouse_to_human_symbols(gene)) %>% drop_na()
iNOS.mtb.pos.strict <- iNOS.mtb.pos %>% filter(avg_log2FC > 2 | avg_log2FC < -2) %>% mutate(change = ifelse(avg_log2FC > 1, "up", "down")) %>% mutate(gene = convert_mouse_to_human_symbols(gene)) %>% drop_na()

b6.mtb.pos.neg <- b6.mtb.pos.neg %>% filter(avg_log2FC > 1 | avg_log2FC < -1) %>% mutate(change = ifelse(avg_log2FC > 1, "up", "down")) %>% mutate(gene = convert_mouse_to_human_symbols(gene)) %>% drop_na()
b6.mtb.pos.neg.strict <- b6.mtb.pos.neg %>% filter(avg_log2FC > 2 | avg_log2FC < -2) %>% mutate(change = ifelse(avg_log2FC > 1, "up", "down")) %>% mutate(gene = convert_mouse_to_human_symbols(gene)) %>% drop_na()

#Export try 2 for gene signatures
write.csv(B6.general, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/B6 Genes.csv")
write.csv(Sp140.general, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/Sp140 Genes.csv")
write.csv(Irg1.general, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/Irg1 Genes.csv")
write.csv(iNOS.general, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/iNOS Genes.csv")
write.csv(b6.mtb.pos.neg, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/Mtb Pos Genes.csv")
write.csv(b6.bystander.neg, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/Bystander Genes.csv")

#Export try 2 for gene signatures
write.csv(B6.strict, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/B6 Strict Genes.csv")
write.csv(Sp140.strict, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/Sp140 Strict Genes.csv")
write.csv(Irg1.strict, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/Irg1 Strict Genes.csv")
write.csv(iNOS.strict, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/iNOS Strict Genes.csv")
write.csv(b6.mtb.pos.neg.strict, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/B6 Mtb Pos Strict.csv")
write.csv(b6.bystander.neg.strict, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/B6 Bystander Strict.csv")
write.csv(Sp140.mtb.pos.strict, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/Sp140 Mtb Pos Strict.csv")
write.csv(Irg1.mtb.pos.strict, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/Irg1 Mtb Pos Strict.csv")
write.csv(iNOS.mtb.pos.strict, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/iNOS Mtb Pos Strict.csv")

#Overview of total cells by infection status and the infected mice by genotype
Idents(super.integrated) <- "infected"
super.infected <- subset(super.integrated, idents = "infected")
Idents(super.infected) <- "celltype"
Idents(super.integrated) <- "celltype"
DimPlot(super.infected, split.by = "cell_type")
DimPlot(super.integrated, split.by = "state")

#Blow up of monocytes and IMs for classification. Example CD63, CD9, and CD16.2 expression
infected.mono <- subset(super.infected, idents = c("CD16-2 Mono","Mono","IM","ISG IM", "Activated IM"))

#5 cells clustered weirdly and look like background noise so they were removed from the analysis.
CellSelector(DimPlot(infected.mono))
cells.remove <- c("GATGATCTCCTAGCCT_4_3","GTCGTAACACTTGGCG_4_3","CACAGGCTCTGAGCAT_6_6","GTAGGTTGTGGATCAG_6_6","ATAGAGACAAGGCGTA_1_7")
infected.mono_filtered <- infected.mono[,!colnames(infected.mono) %in% cells.remove]
DimPlot(infected.mono_filtered)

#Plot marker expression corresopnding to new flow cytometry gating strategy
FeaturePlot(infected.mono_filtered, features = c("Cd63","Cd9","Fcgr4","Ly6c2"))

#General markers for subsets
Idents(super.integrated) <- "celltype"
levels(super.integrated) <- c("gCap","aCap","fibroblast","Stromal?","DC","B cell","Basophil","Camp Mmp8 Neut","Mmp8 Neut","Mature Neut",
                              "Aged Neut","ISG Mature Neut","ISG Activated Neut","ISG Young Neut","ISG Old Neut","Nos2 Neut","CD63 Neut",
                              "AM","Activated AM","Cycling AM","CD16-2 Mono","Mono","IM","Activated IM","ISG IM")
DotPlot(super.integrated, features = c("Ptprc","Cd19","Fcer1a","adt_Ly6G","adt_MHCII","Itgax","adt_SiglecF","Fcgr1","Mertk","Ccr2","Mki67"))

#Markers for subsetting neutrophils
super.neut <- subset(super.integrated, idents = c("Camp Mmp8 Neut","Mmp8 Neut","Mature Neut","Aged Neut","ISG Mature Neut",
                                                  "ISG Activated Neut","ISG Young Neut","ISG Old Neut","Nos2 Neut","CD63 Neut"))
DotPlot(super.neut, features = c("Camp","Mmp8","Cxcr2","Cxcr4","Isg15","Nos2","Cd63"))

#Demonstrate that the scRNA-seq consists of different mice.
Idents(super.integrated) <- "infected"
super.naive <- subset(super.integrated, idents = "naive")
Idents(super.infected) <- "cell_type"
Idents(super.naive) <- "cell_type"
super.infected.inos.irg1 <- subset(super.infected, idents = c("iNOS","Irg1"))
super.naive.inos.irg1 <- subset(super.naive, idents = c("iNOS","Irg1"))
Idents(super.infected.inos.irg1) <- "mouse_id"
Idents(super.naive.inos.irg1) <- "mouse_id"
DimPlot(super.infected.inos.irg1, split.by = "cell_type")
DimPlot(super.naive.inos.irg1, split.by = "cell_type")

#Save current version.
saveRDS(super.integrated, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/b6 sp140 irg1 inos v2")
super.integrated <- readRDS(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/b6 sp140 irg1 inos v2")

simple.integrated <- super.integrated
Idents(simple.integrated) <- "celltype"
simple.integrated <- RenameIdents(simple.integrated, 'ISG Mature Neut' = "Neutrophil", 'ISG Activated Neut' = "Neutrophil", 'Aged Neut' = "Neutrophil", 'Mature Neut' = "Neutrophil",
                                  'ISG Old Neut' = "Neutrophil", 'ISG Young Neut' = "Neutrophil", 'CD63 Neut' = "Neutrophil", 'Mmp8 Neut' = "Neutrophil",'Nos2 Neut' = "Neutrophil",
                                  'Camp Mmp8 Neut' = "Neutrophil", 'CD16-2 Mono' = "Mono", 'ISG IM' = "IM", 'Activated IM' = "IM", 'Activated AM' = "AM")
super.simple <- subset(simple.integrated, idents = c("Neutrophil","Mono","IM","AM","DC"))
super.simple$old.celltype <- super.simple$celltype
super.simple$celltype <- Idents(super.simple)

Idents(super.simple) <- "cell_type"
simp.inos <- subset(super.simple, idents = "iNOS")
simp.irg1 <- subset(super.simple, idents = "Irg1")
simp.sp140 <- subset(super.simple, idents = "Sp140")
Idents(simp.inos) <- "celltype"
Idents(simp.irg1) <- "celltype"
Idents(simp.sp140) <- "celltype"
simp.inos.neut <- subset(simp.inos, idents = "Neutrophil")
simp.irg1.neut <- subset(simp.irg1, idents = "Neutrophil")
simp.sp140.neut <- subset(simp.sp140, idents = "Neutrophil")
Idents(simp.inos.neut) <- "state"
inos.neut.genes <- FindMarkers(simp.inos.neut, ident.1 = "Mtb+", ident.2 = "naive")
Idents(simp.irg1.neut) <- "state"
irg1.neut.genes <- FindMarkers(simp.irg1.neut, ident.1 = "Mtb+", ident.2 = "naive")
Idents(simp.sp140.neut) <- "state"
sp140.neut.genes <- FindMarkers(simp.sp140.neut, ident.1 = "Mtb+", ident.2 = "Mtb-")

Idents(super.simple) <- "celltype"
mtb.infected.neutrophil <- subset(super.simple, idents = "Neutrophil")
Idents(mtb.infected.neutrophil) <- "state"
mtb.infected.neutrophil <- subset(mtb.infected.neutrophil, idents = "Mtb+")
Idents(mtb.infected.neutrophil) <- "cell_type"
sp140.b6.neut.infected.genes <- FindMarkers(mtb.infected.neutrophil, ident.1 = "Sp140", ident.2 = "B6")
irg1.b6.neut.infected.genes <- FindMarkers(mtb.infected.neutrophil, ident.1 = "Irg1", ident.2 = "B6")
inos.b6.neut.infected.genes <- FindMarkers(mtb.infected.neutrophil, ident.1 = "iNOS", ident.2 = "B6")

sp140.b6.neut.infected.volcano <- sp140.b6.neut.infected.genes[c(2,5)]
EnhancedVolcano(sp140.b6.neut.infected.volcano,
                lab = rownames(sp140.b6.neut.infected.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Mtb-infected B6 vs Sp140 Neutrophils",
                subtitle = "",
                raster = T,
                xlim = c(-5,5))

irg1.b6.neut.infected.volcano <- irg1.b6.neut.infected.genes[c(2,5)]
EnhancedVolcano(irg1.b6.neut.infected.volcano,
                lab = rownames(irg1.b6.neut.infected.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Mtb-infected B6 vs Irg1 Neutrophils",
                subtitle = "",
                raster = T,
                xlim = c(-5,5))

inos.b6.neut.infected.volcano <- inos.b6.neut.infected.genes[c(2,5)]
EnhancedVolcano(inos.b6.neut.infected.volcano,
                lab = rownames(inos.b6.neut.infected.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Mtb-infected B6 vs iNOS Neutrophils",
                subtitle = "",
                raster = T,
                xlim = c(-5,5))

Idents(super.integrated) <- "infected"
super.infected <- subset(super.integrated, idents = "infected")
Idents(super.infected) <- "celltype"
super.infected.IM <- subset(super.infected, idents = c("Activated IM","ISG IM","IM","Mono","CD16-2 Mono"))
DimPlot(super.infected.IM)

ISG.IM.infected.genes <- FindMarkers(super.infected.IM, ident.1 = "ISG IM", ident.2 = "Activated IM")
ISG.IM.infected.volcano <- ISG.IM.infected.genes[c(2,5)]
EnhancedVolcano(ISG.IM.infected.volcano,
                lab = rownames(ISG.IM.infected.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Mtb-infected ISG IM vs Activated IM",
                subtitle = "",
                raster = T)

#Get signatures of upregulated genes.
B6.Mtb.Pos.less.stringent

#remove weirdo outliers
cells.remove <- c("GATGATCTCCTAGCCT_4_3","GTCGTAACACTTGGCG_4_3","CACAGGCTCTGAGCAT_6_6","GTAGGTTGTGGATCAG_6_6","ATAGAGACAAGGCGTA_1_7")
super.infected.IM_filtered <- super.infected.IM[,!colnames(super.infected.IM) %in% cells.remove]
Idents(super.infected.IM_filtered) <- "celltype"
DimPlot(super.infected.IM_filtered, split.by = "cell_type")
DimPlot(super.infected.IM_filtered, split.by = "state")
Idents(super.infected.IM_filtered) <- "state"
super.infected.IM_filtered.mtbpos <- subset(super.infected.IM_filtered, idents = "Mtb+")
Idents(super.infected.IM_filtered.mtbpos) <- "celltype"
DimPlot(super.infected.IM_filtered.mtbpos, split.by = "cell_type")

#Volcano plot comparing Spp1 and Cxcl9 IMs
Idents(super.infected.IM_filtered.mtbpos) <- "celltype"
Spp1.IM.genes <- FindMarkers(super.infected.IM_filtered.mtbpos, ident.1 = "ISG IM", ident.2 = "Activated IM")
Spp1.IM.volcano <- Spp1.IM.genes[c(2,5)]
EnhancedVolcano(Spp1.IM.volcano,
                lab = rownames(Spp1.IM.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Mtb-infected Spp1 IM vs CXCL9 IM",
                subtitle = "",
                drawConnectors = T,
                arrowheads = F,
                raster = T)

EnhancedVolcano(Spp1.IM.volcano,
                lab = rownames(Spp1.IM.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "Mtb-infected Spp1 IM vs CXCL9 IM",
                subtitle = "",
                drawConnectors = T,
                arrowheads = F,
                selectLab = c("Cxcl9","Spp1","Il1rn","C1qa","C1qc","C1qb","Trem2","Ccl5","Mmp2","Tnfsf13b","Lyz2","Apoe","Ubd","Dnase1l3","Saa3","C3","Cd9","Inhba","Lgals3","Cd53","Upp1","Il7r","Cxcl3","Ccl4","Ccl2","Ccr1","Phlda1","Slfn4","Clec4a1","C1ra","Clec4n","Edn1","Cd24a","Basp1"),
                raster = T)

#Plot IL-1Ra expression
FeaturePlot(super.infected.IM_filtered.mtbpos, features = "Il1rn")
levels(super.infected.IM_filtered.mtbpos) <- c("ISG IM","Activated IM","IM","CD16-2 Mono","Mono")

#check ifnb expressors
ifnb.pos <- subset(super.integrated, Ifnb1 > 0)

#Defining IFNg and IFNb responders using the mouse signatures
cyto.sig <- list()
cyto.sig$mouse_ifnb <- specific.ifnb.mouse$gene
cyto.sig$mouse_ifng <- specific.ifng.mouse$gene

super.integrated <- AddModuleScore_UCell(super.integrated, features = cyto.sig)

super.integrated$mouse_ifn_sig <- ifelse(super.integrated$mouse_ifnb_UCell >= 0.24 & super.integrated$mouse_ifng_UCell >= 0.32, "both", 
                                        ifelse(super.integrated$mouse_ifnb_UCell >= 0.24 & super.integrated$mouse_ifng_UCell < 0.32, "ifnb", 
                                               ifelse(super.integrated$mouse_ifnb_UCell < 0.24 & super.integrated$mouse_ifng_UCell >= 0.32, "ifng","neither")))

Idents(super.integrated) <- "mouse_ifn_sig"
DimPlot(super.integrated, split.by = "state",cols = c("grey","purple","red","blue"), order = c("ifnb","ifng","both","neither"))
DimPlot(super.integrated, split.by = "cell_type",cols = c("grey","purple","red","blue"), order = c("ifnb","ifng","both","neither"))
Idents(super.integrated) <- "state"
super.mtb <- subset(super.integrated, idents = "Mtb+")
Idents(super.mtb) <- "mouse_ifn_sig"
DimPlot(super.mtb, split.by = "cell_type",cols = c("grey","purple","red","blue"), order = c("ifnb","ifng","both","neither"))

Idents(super.integrated) <- "state"
super.naive <- subset(super.integrated, idents = "naive")
Idents(super.naive) <- "mouse_ifn_sig"
DimPlot(super.naive, cols = c("grey","purple","red","blue"), order = c("ifnb","ifng","both","neither"))

#The problem is that not every group has every cell type. Mono, IM, Neutrophil are in every sample. DC and AM are in 5 out of 6 samples (not sample 2).
Idents(super.mtb) <- "celltype"
simple.mtb <- RenameIdents(super.mtb, 'ISG Mature Neut' = "Neutrophil", 'ISG Activated Neut' = "Neutrophil", 'Aged Neut' = "Neutrophil", 'Mature Neut' = "Neutrophil",
                                  'ISG Old Neut' = "Neutrophil", 'ISG Young Neut' = "Neutrophil", 'CD63 Neut' = "Neutrophil", 'Mmp8 Neut' = "Neutrophil",'Nos2 Neut' = "Neutrophil",
                                  'Camp Mmp8 Neut' = "Neutrophil", 'CD16-2 Mono' = "Mono", 'ISG IM' = "IM", 'Activated IM' = "IM", 'Activated AM' = "AM")
simple.mtb <- subset(simple.mtb, idents = c("Neutrophil","Mono","IM","AM","DC"))
simple.mtb$old.celltype <- simple.mtb$celltype
simple.mtb$celltype <- Idents(simple.mtb)

simp.genotype <- list()
Idents(simple.mtb) <- "mouse_id"

populations <- as.character(unique(simple.mtb$celltype))
cell_type <- rep(unique(simple.mtb$cell_type),3)

for(z in 1:3){
  a <- subset(simple.mtb, idents = z)
  Idents(a) <- "cell_type"
  for(i in 1:4){
    y <- assign(paste("simp.", cell_type[i], sep = ""), subset(a, idents = cell_type[i]))
    simp.genotype[[length(simp.genotype)+1]] <- y
    Idents(simp.genotype[[length(simp.genotype)]]) <- "celltype"
  }
}


#Missing populations
#1 - AM, DC
#3 - DC
#7 - DC
#11 - DC
cell_type.pop <- list()
for(i in 1:3){
  for(j in 1:12){
    y <- assign(paste("simp", cell_type[j], populations[i], sep = "."), subset(simp.genotype[[j]], idents = populations[i]))
    cell_type.pop[[length(cell_type.pop)+1]] <- y
  }
}

cyt.sig.total <- lapply(cell_type.pop, function(x) {(prop.table(table(Idents(x), x$mouse_ifn_sig)))})
cyt.sig.total <- lapply(cyt.sig.total, function(x) {as.data.frame.matrix(x)})
cyt.sig.total <- data.table:::rbindlist(cyt.sig.total, fill=TRUE)
cyt.sig.total[is.na(cyt.sig.total)] = 0
cyt.sig.total <- cyt.sig.total[]*100
cyt.sig.total$mouse.id <- rep(c(rep(1,4),rep(2,4),rep(3,4)),3)
cyt.sig.total$genotype <- c(rep(cell_type,3))
cyt.sig.total$celltype <- c(rep("IM",12),rep("Neutrophil",12),rep("Monocyte",12))
write.csv(cyt.sig.total, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/Plots 230806/IFN Response IM Neut Mono.csv")

#For AM exclude 1
cell_type.pop <- list()
for(i in 4){
  for(j in c(2:12)){
    y <- assign(paste("simp", cell_type[j], populations[i], sep = "."), subset(simp.genotype[[j]], idents = populations[i]))
    cell_type.pop[[length(cell_type.pop)+1]] <- y
  }
}

cyt.sig.total <- lapply(cell_type.pop, function(x) {(prop.table(table(Idents(x), x$mouse_ifn_sig)))})
cyt.sig.total <- lapply(cyt.sig.total, function(x) {as.data.frame.matrix(x)})
cyt.sig.total <- data.table:::rbindlist(cyt.sig.total, fill=TRUE)
cyt.sig.total[is.na(cyt.sig.total)] = 0
cyt.sig.total <- cyt.sig.total[]*100
cyt.sig.total$mouse.id <- c(rep(1,3),rep(2,4),rep(3,4))
cyt.sig.total$genotype <- c("Sp140","Irg1","iNOS",rep(c("B6","Sp140","Irg1","iNOS"),2))
cyt.sig.total$celltype <- c(rep("AM",11))
write.csv(cyt.sig.total, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/Plots 230806/IFN Response AM.csv")

#Frequency of the IM populations by genotype for Mtb-infected cells:
Idents(super.infected.IM_filtered.mtbpos) <- "celltype"
simple.IM.mtb <- subset(super.infected.IM_filtered.mtbpos, idents = c("Activated IM","ISG IM","IM"))
simple.IM.mtb$old.celltype <- simple.IM.mtb$celltype
simple.IM.mtb$celltype <- Idents(simple.IM.mtb)

simp.genotype <- list()
Idents(simple.IM.mtb) <- "mouse_id"

populations <- as.character(unique(simple.IM.mtb$celltype))
cell_type <- rep(unique(simple.IM.mtb$cell_type),3)

for(z in 1:3){
  a <- subset(simple.IM.mtb, idents = z)
  Idents(a) <- "cell_type"
  for(i in 1:4){
    y <- assign(paste("simp.", cell_type[i], sep = ""), subset(a, idents = cell_type[i]))
    simp.genotype[[length(simp.genotype)+1]] <- y
    Idents(simp.genotype[[length(simp.genotype)]]) <- "celltype"
  }
}

im.total <- lapply(simp.genotype, function(x) {(prop.table(table(Idents(x), x$celltype)))})
im.total <- lapply(im.total, function(x) {as.data.frame.matrix(x)})
im.total <- data.table:::rbindlist(im.total, fill=TRUE)
im.total[is.na(im.total)] = 0
im.total <- im.total[]*100
im.total$mouse.id <- rep(c(rep(1,12),rep(2,12),rep(3,12)))
im.total$genotype <- rep(c(rep("B6",3),rep("Sp140",3),rep("Irg1",3),rep("iNOS",3)),3)
write.csv(im.total, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/Plots 240116/Mtb-infected IM Subsets by Genotype.csv")

#Plot Il1rn and IL1r2 levels by state and genotype
levels(super.integrated) <- c("B6_naive","B6_Mtb-","B6_Mtb+","Sp140_naive","Sp140_Mtb-","Sp140_Mtb+","Irg1_naive","Irg1_Mtb-","Irg1_Mtb+","iNOS_naive","iNOS_Mtb-","iNOS_Mtb+")
VlnPlot(super.integrated, features = c("Il1r2","Il1rn"))