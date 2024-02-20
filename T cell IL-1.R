library(Seurat)
library(dplyr)
library(patchwork)
library(cowplot)
library(multtest)
library(metap)
library(ggplot2)
library(EnhancedVolcano)

# Load the infected cell dataset - use the 30000 expected cell HTO and ADT mapping -remap HTO for B6 and Sp140 using 5 tags.
T.rna <- Read10X(data.dir = "/Users/DmitriKotov/Documents/Sequencing Data/scRNA-seq 12-30-20/RVDK003A/outs/filtered_feature_bc_matrix", strip.suffix = TRUE)
T.hto <- Read10X(data.dir = "/Users/DmitriKotov/Documents/Sequencing Data/scRNA-seq 12-30-20/HTO RVDK003E/umi_count", gene.column = 1)
T.adt <- Read10X(data.dir = "/Users/DmitriKotov/Documents/Sequencing Data/scRNA-seq 12-30-20/ADT RVDK003C/umi_count", gene.column = 1)

#Keep cells with data for RNA, ADT, and HTO - T cell
joint.T = Reduce("intersect", list(colnames(T.rna), colnames(T.hto), colnames(T.adt)))
T.rna = T.rna[, joint.T]
T.hto = as.matrix(T.hto[-nrow(T.hto), joint.T])
T.adt = as.matrix(T.adt[-nrow(T.adt), joint.T])

# Confirm that the HTO have the correct names
rownames(T.hto)
rownames(T.hto) <- c("HTO_D", "HTO_E", "HTO_F")

# Confirm that the ADT have the correct names
rownames(T.adt)
rownames(T.adt) <- c("Ly6C", "CD44", "TCRb","PD-L1","CX3CR1","CD49a","NKp46","CD4","ST2","CD25","CD45-1","CD45-2","CD8a","CXCR5","CD62L","CCR6")

# Initialize the Seurat object with the raw (non-normalized data).
Tcell <- CreateSeuratObject(counts = T.rna, project = "IL-1")
Tcell #check seurat object

# Add HTO data as a new assay independent from RNA
Tcell[["HTO"]] <- CreateAssayObject(counts = T.hto)

# Add ADT data as a new assay independent from RNA
Tcell[["ADT"]] <- CreateAssayObject(counts = T.adt)

#check for dying cells based on mitochondrial percentage
Tcell[["percent.mt"]] <- PercentageFeatureSet(Tcell, pattern = "^mt-")

#standard filters are 
VlnPlot(Tcell, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# standard log-normalization
Tcell <- NormalizeData(Tcell)

# choose ~1k variable features
Tcell <- FindVariableFeatures(Tcell)

# standard scaling (no regression)
Tcell <- ScaleData(Tcell)

Tcell <- RunPCA(Tcell, features = VariableFeatures(object = Tcell))
ElbowPlot(Tcell) #dims = 17
#HTO demultiplexing

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
Tcell <- NormalizeData(Tcell, assay = "HTO", normalization.method = "CLR")

#Run HTO demux
Tcell <- HTODemux(Tcell, assay = "HTO", positive.quantile = 0.99)

# Global classification results
table(Tcell$HTO_classification.global)

# Group cells based on the max HTO signal
Idents(Tcell) <- "HTO_maxID"

#Visualize multiplets
RidgePlot(Tcell, assay = "HTO", features = rownames(Tcell[["HTO"]])[1:2], ncol = 2)
FeatureScatter(Tcell, feature1 = "hto_HTO-D", feature2 = "hto_HTO-E")
Idents(Tcell) <- "HTO_classification.global"
VlnPlot(Tcell, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
VlnPlot(Tcell, features = "nFeature_RNA", pt.size = 0.1, log = TRUE)

# First, we will remove negative cells from the object
Tcell.subset <- subset(Tcell, idents = "Negative", invert = TRUE)

# Calculate a tSNE embedding of the HTO data
DefaultAssay(Tcell.subset) <- "HTO"
Tcell.subset <- ScaleData(Tcell.subset, features = rownames(Tcell.subset), 
                          verbose = FALSE)
Tcell.subset <- RunPCA(Tcell.subset, features = rownames(Tcell.subset), approx = FALSE)
Tcell.subset <- RunTSNE(Tcell.subset, dims = 1:3, perplexity = 100, check_duplicates = FALSE) #check why there are duplicates
DimPlot(Tcell.subset)

# Extract the singlets for T cells
Tcell.singlet <- subset(Tcell, idents = "Singlet")

#Filter to exclude cells that have unique feature counts more than 2500 (multiplets) and less than 200 (empty drops) and removing cells that have >5% mitochondrial counts (dead/dying cells)
VlnPlot(Tcell.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Tcell.singlet <- subset(Tcell.singlet, subset = nFeature_RNA < 2500 & nFeature_RNA > 200 & percent.mt < 5)

# Process ADT data and  set a dimensional reduction name to avoid overwriting the RNA PCA
DefaultAssay(Tcell.singlet) <- "ADT"
VariableFeatures(Tcell.singlet) <- rownames(Tcell.singlet[["ADT"]])
Tcell.singlet <- NormalizeData(Tcell.singlet, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

Tcell.singlet <- FindMultiModalNeighbors(
  Tcell.singlet, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:15), modality.weight.name = "RNA.weight"
)

Tcell.singlet <- RunUMAP(Tcell.singlet, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
Tcell.singlet <- FindClusters(Tcell.singlet, graph.name = "wsnn", algorithm = 3, resolution = 1, verbose = FALSE)

T.plot <- DimPlot(Tcell.singlet, reduction = "wnn.umap", label = TRUE, repel = TRUE)

#RNA contribution to each cluster
VlnPlot(Tcell.singlet, features = "RNA.weight", sort = TRUE, pt.size = 0.1) +
  NoLegend()

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
library(nichenetr)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

m.s.genes <-convert_human_to_mouse_symbols(s.genes)
m.g2m.genes <- convert_human_to_mouse_symbols(g2m.genes)

Tcell.singlet <- CellCycleScoring(Tcell.singlet, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
# view cell cycle scores and phase assignments
head(Tcell.singlet[[]])
#visualize different cell cycle phases
DimPlot(Tcell.singlet, reduction = "wnn.umap")
#switch ident back to clusters
Idents(Tcell.singlet) <- "old.ident"

#Markers for T cell populations
DotPlot(Tcell.singlet, features = c("Cd3e","Mki67","Cd4","Cd8a","Cx3cr1","Cxcr3","Sell","adt_CD44","adt_CXCR5","Il17a","Foxp3","Isg15"))

#Name clusters
Tcell.singlet <- RenameIdents(Tcell.singlet, `0` = "Th1", `1` = "naive CD4", `2` = "CX3CR1 Th1", 
                              `3` = "CX3CR1 CD8", `4` = "naive CD8", `5` = "CXCR3 CD8", `6` = "Tcm CD8", `7` = "CX3CR1 Th1", `8` = "Int Th1", `9` = "Int Th1", 
                              `10` = "Treg", `11` = "Th1", `12` = "CX3CR1 CD8", `13` = "Th17",`14` = "Cycling T",`15` = "Cycling T",
                              `16` = "naive CD4",`17` = "Tfh",`18` = "naive CD8",`19` = "dying cells",`20` = "ISG+ CD4")

#Split Th17 cluster into Th17 and CD4 effector memory?

#Re-order clusters to group related cell types
levels(x = Tcell.singlet) <- c("Cycling T","dying cells","ISG+ CD4","naive CD4","Th1","Int Th1","CX3CR1 Th1","Treg","Th17","Tfh","naive CD8","Tcm CD8","CXCR3 CD8","CX3CR1 CD8")

Tcell.singlet$celltype <- Idents(Tcell.singlet)

#Figure out how to distinguish CD45.1 and CD45.2 cells - CE45.1+ seems to be at least > 0.5
CD45.plot <- FeatureScatter(Tcell.singlet, feature1 = "adt_CD45-1", feature2 = "adt_CD45-2")

#select the KO cells (CD45.1- CD45.2+)
select.cells <- CellSelector(plot = CD45.plot)
IL1RKO.cells <- select.cells

#select the KO cells (CD45.1+ CD45.2+)
select.cells <- CellSelector(plot = CD45.plot)
WT.cells <- select.cells
CD45.meta <- data.frame("BM" = 1:c(length(IL1RKO.cells) + length(WT.cells)))
CD45.meta[1:length(IL1RKO.cells),] <- "IL-1RKO" 
CD45.meta[c(length(IL1RKO.cells) + 1):c(length(IL1RKO.cells) + length(WT.cells)),] <- "WT"
rownames(CD45.meta) <- c(IL1RKO.cells,WT.cells)
Tcell.singlet <- AddMetaData(Tcell.singlet, CD45.meta, col.name = "BM")
Idents(Tcell.singlet) <- "BM"
DimPlot(Tcell.singlet, reduction = "wnn.umap")
Idents(Tcell.singlet) <-"celltype"

#add metadata to immune.combined to specify infection status and cluster in one value.
Tcell.singlet$celltype.BM <- paste(Idents(Tcell.singlet), Tcell.singlet$BM, sep = "_")

#Save Tcell.singlet file for later recall
saveRDS(Tcell.singlet, file = "/Users/DmitriKotov/Box Sync/Coding stuff/T cell ILC IL-1 scRNA-seq 12-30-20/Tcellsinglet")

#load saved ILC.singlet file
Tcell.singlet <- readRDS(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/T cell ILC IL-1 scRNA-seq 12-30-20/Tcellsinglet")

#Visualize differences between WT and IL-1R KO cells for various clusters.
potential.ligands <- unique(lr_network$from) %>% convert_human_to_mouse_symbols() %>% na.omit
unregulated.Th17 <- Th17.volcano %>% filter(avg_log2FC < 1) %>% filter(avg_log2FC > -1) %>% row.names()
potential.ligands.Th17 <- setdiff(potential.ligands, unregulated.Th17)

regulated.Th17 <- Th17.volcano %>% filter(avg_log2FC > 1 | avg_log2FC < -1)  %>% row.names()
Th17.ligands <- intersect(potential.ligands, regulated.Th17)


Idents(Tcell.singlet) <- "celltype.BM"
Th17.volcano <- FindMarkers(Tcell.singlet, ident.1 = "Th17_WT", ident.2 = "Th17_IL-1RKO", logfc.threshold = 0)
Th17.volcano <- Th17.volcano[c(1,2)]
EnhancedVolcano(Th17.volcano,
                lab = rownames(Th17.volcano),
                x = 'avg_log2FC',
                y = 'p_val',
                title = "Th17 WT vs IL-1R KO",
                subtitle = "",
                pCutoff = 0.05,
                FCcutoff = 1,
                selectLab = Th17.ligands,
                drawConnectors = TRUE,
                arrowheads = FALSE,
                labSize = 10)

Th1.volcano <- FindMarkers(Tcell.singlet, ident.1 = "Th1_WT", ident.2 = "Th1_IL-1RKO")
Th1.volcano <- Th1.volcano[c(1,2)]
EnhancedVolcano(Th1.volcano,
                lab = rownames(Th1.volcano),
                x = 'avg_log2FC',
                y = 'p_val',
                title = "Th1 WT vs IL-1R KO",
                subtitle = "",
                pCutoff = 0.05,
                FCcutoff = 0.75)

CX3CR1.Th1.volcano <- FindMarkers(Tcell.singlet, ident.1 = "CX3CR1 Th1_WT", ident.2 = "CX3CR1 Th1_IL-1RKO")
CX3CR1.Th1.volcano <- CX3CR1.Th1.volcano[c(1,2)]
EnhancedVolcano(CX3CR1.Th1.volcano,
                lab = rownames(CX3CR1.Th1.volcano),
                x = 'avg_log2FC',
                y = 'p_val',
                title = "CX3CR1 Th1 WT vs IL-1R KO",
                subtitle = "",
                pCutoff = 0.05,
                FCcutoff = 0.75)

CX3CR1.Th1.volcano <- FindMarkers(Tcell.singlet, ident.1 = "CX3CR1 Th1_WT", ident.2 = "CX3CR1 Th1_IL-1RKO",assay = "ADT")
CX3CR1.Th1.volcano <- CX3CR1.Th1.volcano[c(1,2)]
EnhancedVolcano(CX3CR1.Th1.volcano,
                lab = rownames(CX3CR1.Th1.volcano),
                x = 'avg_log2FC',
                y = 'p_val',
                title = "CX3CR1 Th1 WT vs IL-1R KO",
                subtitle = "",
                pCutoff = 0.05,
                FCcutoff = 0.75)

CX3CR1.CD8.volcano <- FindMarkers(Tcell.singlet, ident.1 = "CX3CR1 CD8_WT", ident.2 = "CX3CR1 CD8_IL-1RKO")
CX3CR1.CD8.volcano <- CX3CR1.CD8.volcano[c(1,2)]
EnhancedVolcano(CX3CR1.CD8.volcano,
                lab = rownames(CX3CR1.CD8.volcano),
                x = 'avg_log2FC',
                y = 'p_val',
                title = "CX3CR1 CD8 WT vs IL-1R KO",
                subtitle = "",
                pCutoff = 0.05,
                FCcutoff = 0.75)

CX3CR1.CD8.volcano <- FindMarkers(Tcell.singlet, ident.1 = "CX3CR1 CD8_WT", ident.2 = "CX3CR1 CD8_IL-1RKO", assay = "ADT")
CX3CR1.CD8.volcano <- CX3CR1.CD8.volcano[c(1,2)]
EnhancedVolcano(CX3CR1.CD8.volcano,
                lab = rownames(CX3CR1.CD8.volcano),
                x = 'avg_log2FC',
                y = 'p_val',
                title = "CX3CR1 CD8 WT vs IL-1R KO",
                subtitle = "",
                pCutoff = 0.05,
                FCcutoff = 0.75)

CXCR3.CD8.volcano <- FindMarkers(Tcell.singlet, ident.1 = "CXCR3 CD8_WT", ident.2 = "CXCR3 CD8_IL-1RKO")
CXCR3.CD8.volcano <- CXCR3.CD8.volcano[c(1,2)]
EnhancedVolcano(CXCR3.CD8.volcano,
                lab = rownames(CXCR3.CD8.volcano),
                x = 'avg_log2FC',
                y = 'p_val',
                title = "CXCR3 CD8 WT vs IL-1R KO",
                subtitle = "",
                pCutoff = 0.05,
                FCcutoff = 0.75)

Tfh.volcano <- FindMarkers(Tcell.singlet, ident.1 = "Tfh_WT", ident.2 = "Tfh_IL-1RKO")
Tfh.volcano <- Tfh.volcano[c(1,2)]
EnhancedVolcano(Tfh.volcano,
                lab = rownames(Tfh.volcano),
                x = 'avg_log2FC',
                y = 'p_val',
                title = "Tfh WT vs IL-1R KO",
                subtitle = "",
                pCutoff = 0.05,
                FCcutoff = 0.75)

Treg.volcano <- FindMarkers(Tcell.singlet, ident.1 = "Treg_WT", ident.2 = "Treg_IL-1RKO")
Treg.volcano <- Treg.volcano[c(1,2)]
EnhancedVolcano(Treg.volcano,
                lab = rownames(Treg.volcano),
                x = 'avg_log2FC',
                y = 'p_val',
                title = "Treg WT vs IL-1R KO",
                subtitle = "",
                pCutoff = 0.05,
                FCcutoff = 0.75)

Idents(Tcell.singlet) <- "BM"
wt.volcano <- FindMarkers(Tcell.singlet, ident.1 = "WT", ident.2 = "IL-1RKO")
wt.volcano <- wt.volcano[c(1,2)]
EnhancedVolcano(wt.volcano,
                lab = rownames(wt.volcano),
                x = 'avg_log2FC',
                y = 'p_val',
                title = "WT vs IL-1R KO",
                subtitle = "",
                pCutoff = 0.05,
                FCcutoff = 0.75)

#Do nichenet analysis and compare with the LogFC on senders graph using genes up in WT over IL1R KO with Log2FC > 0.3 (1.25X in WT over IL-1R KO)
library(nichenetr)
library(tidyverse)

#load in myeloid / granulocyte seurat object
mac.neut <- readRDS(file = "/Users/DmitriKotov/Box Sync/Coding stuff/B6 vs Sp140 scRNA-seq 10-9-20/immunecombined2")

#Nichenet setup
ligand_target_matrix = readRDS(url("https://zenodo.org/record/3260758/files/ligand_target_matrix.rds"))
lr_network = readRDS(url("https://zenodo.org/record/3260758/files/lr_network.rds"))
weighted_networks = readRDS(url("https://zenodo.org/record/3260758/files/weighted_networks.rds"))
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))

#eliminate predicted protein interaction datasets from the NicheNet analysis
lr_network = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")

#convert the loaded in networks to mouse symbols
lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()

ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

weighted_networks_lr = weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()

#Nichenet analysis on IM interactions with various CD4 and CD8 T cell senders with IM vs AM DEG.

## receiver
receiver = "IM"
expressed_genes_receiver = get_expressed_genes(receiver, assay_oi = "RNA", mac.neut, pct = 0.10)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## sender
sender_celltypes = c("Th17")
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, Tcell.singlet, 0.10, assay_oi = "RNA") # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

#Define a genset of interest - 
DE_table_receiver = FindMarkers(object = mac.neut, assay = "RNA", ident.1 = "IM", ident.2 = "AM", min.pct = 0.10) %>% rownames_to_column("gene")
geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#Define a set of potential ligands
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

#Preform NicheNet analysis to define ligand activity
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))

#use pearons = 0.1 as cutoff for downstream analysis
best_upstream_ligands = ligand_activities %>% filter(pearson > 0.1) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()

cluster.ligand <- DotPlot(Tcell.singlet, assay = "RNA", idents = sender_celltypes, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis() + coord_flip()

#active target inference
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
p_ligand_target_network

#receptors of top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()

p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network

#expression of receptors for prioritized ligands
p_receptor_expression = DotPlot(mac.neut, idents = c("IFN IM","IM","AM","Mono") , assay = "RNA", features = rownames(vis_ligand_receptor_network), cols = "RdYlBu") + RotatedAxis()
p_receptor_expression

#Had to change the get_lfc_celltype function because new version of Seurat labels column as avg_log2FC rather than avg_logFC
get_lfc_celltype = function(celltype_oi, seurat_obj, condition_colname, condition_oi, condition_reference, expression_pct = 0.10){
  requireNamespace("Seurat")
  requireNamespace("dplyr")
  seurat_obj_celltype = SetIdent(seurat_obj, value = seurat_obj[["celltype"]])
  seuratObj_sender = subset(seurat_obj_celltype, idents = celltype_oi)
  seuratObj_sender = SetIdent(seuratObj_sender, value = seuratObj_sender[[condition_colname]])
  DE_table_sender = FindMarkers(object = seuratObj_sender, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = expression_pct, logfc.threshold = 0.05) %>% rownames_to_column("gene")
  DE_table_sender = DE_table_sender %>% as_tibble() %>% dplyr::select(-p_val) %>% dplyr::select(gene, avg_log2FC)
  colnames(DE_table_sender) = c("gene",celltype_oi)
  return(DE_table_sender)
}

DE_table_all = Idents(Tcell.singlet) %>% levels() %>% intersect(sender_celltypes) %>% lapply(get_lfc_celltype, seurat_obj = Tcell.singlet, condition_colname = "BM", condition_oi = "WT", condition_reference = "IL-1RKO", expression_pct = 0.10) %>% reduce(full_join)

DE_table_all[is.na(DE_table_all)] = 0

# Combine ligand activities with DE information
ligand_activities_de = ligand_activities %>% dplyr::select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0

# make LFC heatmap
lfc_matrix = ligand_activities_de  %>% dplyr::select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = as.matrix(lfc_matrix[order_ligands,])

colnames(vis_ligand_lfc) = "Th17"

p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
p_ligand_lfc

#IL-1R1 expression by the sender populations
p_Il1r1_sender <- DotPlot(Tcell.singlet, idents = sender_celltypes, features = c("Il1r1")) + coord_flip()
p_Il1r1_sender

#Alternative visualization that I like less
#p_ligand_lfc = p_ligand_lfc + scale_fill_gradientn(colors = c("midnightblue","blue", "grey95", "grey99","firebrick1","red"),values = c(0,0.1,0.2,0.25, 0.40, 0.7,1), limits = c(vis_ligand_lfc %>% min() - 0.1, vis_ligand_lfc %>% max() + 0.1))
#p_ligand_lfc

#display ligand pearson alognside p_ligand_lfc
# ligand activity heatmap
ligand_pearson_matrix = ligand_activities %>% dplyr::select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))

p_ligand_pearson + cluster.ligand + p_ligand_lfc + p_ligand_receptor_network + p_receptor_expression

#Based on this analysis replot ligand expression and receptor expresssion for Th17 and IM populations for ligands with positive LFC for WT vs IL-1R KO.
#Tnf, Ccl5, Icam1, Ifng, Cd274, Itgb2, Tgfb1, Il17a
lfc_ligands = c("Tnf","Ccl5","Icam1","Cd274","Itgb2","Tgfb1","Il17a")
lfc.cluster.ligand <- DotPlot(Tcell.singlet, assay = "RNA", idents = c("Th17","Th1","CX3CR1 Th1","naive CD4"), features = lfc_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis() + coord_flip()

lfc_receptors = c("Tnfrsf1a","Tnfrsf1b","Ccr1","Ccr5","Itgax","Cd80","Icam1","Itgal","Tgfbr1","Tgfbr2","Il17ra")
p_lfc_receptor_expression = DotPlot(mac.neut, idents = c("IFN IM","IM","AM","Mono") , assay = "RNA", features = lfc_receptors, cols = "RdYlBu") + RotatedAxis()
p_lfc_receptor_expression

lfc.cluster.ligand + p_lfc_receptor_expression

#Further filtered ligands based on receptor expression pattern
lfc_receptor_ligands = c("Tnf","Ccl5","Cd274","Itgb2","Il17a")
lfc.receptor.cluster.ligand <- DotPlot(Tcell.singlet, assay = "RNA", idents = c("Th17","Th1","CX3CR1 Th1","naive CD4"), features = lfc_receptor_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis() + coord_flip()

#Rerun nichenet with IM vs Mono DEG, IM vs Lyve1 IM DEG

#Run nichenet on IFN IM recievers and IFN IM vs Lyve1 IM DEG

## receiver
receiver = "IFN IM"
expressed_genes_receiver = get_expressed_genes(receiver, assay_oi = "RNA", mac.neut, pct = 0.10)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## sender
sender_celltypes = c("Th17")
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, Tcell.singlet, 0.10, assay_oi = "RNA") # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

#Define a genset of interest - 
DE_table_receiver = FindMarkers(object = mac.neut, assay = "RNA", ident.1 = "IFN IM", ident.2 = "Lyve1 IM", min.pct = 0.10) %>% rownames_to_column("gene")
geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.25) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#Define a set of potential ligands
ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

#Preform NicheNet analysis to define ligand activity
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))

#use pearons = 0.1 as cutoff for downstream analysis
best_upstream_ligands = ligand_activities %>% filter(pearson > 0.1) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()

cluster.ligand <- DotPlot(Tcell.singlet, assay = "RNA", idents = sender_celltypes, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis() + coord_flip()

#active target inference
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
p_ligand_target_network

#receptors of top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()

p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network

#expression of receptors for prioritized ligands
p_receptor_expression = DotPlot(mac.neut, idents = c("IFN IM","IM","AM","Mono") , assay = "RNA", features = rownames(vis_ligand_receptor_network), cols = "RdYlBu") + RotatedAxis()
p_receptor_expression

#Had to change the get_lfc_celltype function because new version of Seurat labels column as avg_log2FC rather than avg_logFC
get_lfc_celltype = function(celltype_oi, seurat_obj, condition_colname, condition_oi, condition_reference, expression_pct = 0.10){
  requireNamespace("Seurat")
  requireNamespace("dplyr")
  seurat_obj_celltype = SetIdent(seurat_obj, value = seurat_obj[["celltype"]])
  seuratObj_sender = subset(seurat_obj_celltype, idents = celltype_oi)
  seuratObj_sender = SetIdent(seuratObj_sender, value = seuratObj_sender[[condition_colname]])
  DE_table_sender = FindMarkers(object = seuratObj_sender, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = expression_pct, logfc.threshold = 0.05) %>% rownames_to_column("gene")
  DE_table_sender = DE_table_sender %>% as_tibble() %>% dplyr::select(-p_val) %>% dplyr::select(gene, avg_log2FC)
  colnames(DE_table_sender) = c("gene",celltype_oi)
  return(DE_table_sender)
}

DE_table_all = Idents(Tcell.singlet) %>% levels() %>% intersect(sender_celltypes) %>% lapply(get_lfc_celltype, seurat_obj = Tcell.singlet, condition_colname = "BM", condition_oi = "WT", condition_reference = "IL-1RKO", expression_pct = 0.10) %>% reduce(full_join)

DE_table_all[is.na(DE_table_all)] = 0

# Combine ligand activities with DE information
ligand_activities_de = ligand_activities %>% dplyr::select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0

# make LFC heatmap
lfc_matrix = ligand_activities_de  %>% dplyr::select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = as.matrix(lfc_matrix[order_ligands,])

colnames(vis_ligand_lfc) = "Th17"

p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
p_ligand_lfc

#IL-1R1 expression by the sender populations
p_Il1r1_sender <- DotPlot(Tcell.singlet, idents = sender_celltypes, features = c("Il1r1")) + coord_flip()
p_Il1r1_sender

#Alternative visualization that I like less
#p_ligand_lfc = p_ligand_lfc + scale_fill_gradientn(colors = c("midnightblue","blue", "grey95", "grey99","firebrick1","red"),values = c(0,0.1,0.2,0.25, 0.40, 0.7,1), limits = c(vis_ligand_lfc %>% min() - 0.1, vis_ligand_lfc %>% max() + 0.1))
#p_ligand_lfc

#display ligand pearson alognside p_ligand_lfc
# ligand activity heatmap
ligand_pearson_matrix = ligand_activities %>% dplyr::select(pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

rownames(ligand_pearson_matrix) = rownames(ligand_pearson_matrix) %>% make.names()
colnames(ligand_pearson_matrix) = colnames(ligand_pearson_matrix) %>% make.names()

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)") + theme(legend.text = element_text(size = 9))

p_ligand_pearson + cluster.ligand + p_ligand_lfc + p_ligand_receptor_network + p_receptor_expression

#Plot IL-1 regulated genes expressed by Th17 cells.
library(tidyverse)
library(nichenetr)

Idents(Tcell.singlet) <- "celltype.BM"
Th17.IL1 <- FindMarkers(Tcell.singlet, ident.1 = "Th17_WT", ident.2 = "Th17_IL-1RKO")

inter_Th17_IL1 = Th17.IL1 %>% filter(avg_log2FC > 0.5 & pct.1 > 0.1)
