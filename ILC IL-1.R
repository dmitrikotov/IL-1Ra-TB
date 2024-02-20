library(Seurat)
library(dplyr)
library(patchwork)
library(cowplot)
library(multtest)
library(metap)
library(ggplot2)
library(EnhancedVolcano)

# Load the infected cell dataset - use the 30000 expected cell HTO and ADT mapping -remap HTO for B6 and Sp140 using 5 tags.
ILC.rna <- Read10X(data.dir = "/Users/DmitriKotov/Documents/Sequencing Data/scRNA-seq 12-30-20/RVDK003B/outs/filtered_feature_bc_matrix", strip.suffix = TRUE)
ILC.hto <- Read10X(data.dir = "/Users/DmitriKotov/Documents/Sequencing Data/scRNA-seq 12-30-20/HTO RVDK003F/umi_count", gene.column = 1)
ILC.adt <- Read10X(data.dir = "/Users/DmitriKotov/Documents/Sequencing Data/scRNA-seq 12-30-20/ADT RVDK003D/umi_count", gene.column = 1)

#Keep cells with data for RNA, ADT, and HTO - ILC
joint.ILC = Reduce("intersect", list(colnames(ILC.rna), colnames(ILC.hto), colnames(ILC.adt)))
ILC.rna = ILC.rna[, joint.ILC]
ILC.hto = as.matrix(ILC.hto[-nrow(ILC.hto), joint.ILC])
ILC.adt = as.matrix(ILC.adt[-nrow(ILC.adt), joint.ILC])

# Confirm that the HTO have the correct names
rownames(ILC.hto)
rownames(ILC.hto) <- c("HTO_A", "HTO_B", "HTO_C")

# Confirm that the ADT have the correct names
rownames(ILC.adt)
rownames(ILC.adt) <- c("Ly6C", "CD44", "TCRb","PD-L1","CX3CR1","CD49a","NKp46","CD4","ST2","CD25","CD45-1","CD45-2","CD8a","CXCR5","CD62L","CCR6")

# Initialize the Seurat object with the raw (non-normalized data).
ILC <- CreateSeuratObject(counts = ILC.rna, project = "IL-1")
ILC #check seurat object

# Add HTO data as a new assay independent from RNA
ILC[["HTO"]] <- CreateAssayObject(counts = ILC.hto)

# Add ADT data as a new assay independent from RNA
ILC[["ADT"]] <- CreateAssayObject(counts = ILC.adt)

#check for dying cells based on mitochondrial percentage
ILC[["percent.mt"]] <- PercentageFeatureSet(ILC, pattern = "^mt-")

#standard filters are 
VlnPlot(ILC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# standard log-normalization
ILC <- NormalizeData(ILC)

# choose ~1k variable features
ILC <- FindVariableFeatures(ILC)

# standard scaling (no regression)
ILC <- ScaleData(ILC)

ILC <- RunPCA(ILC, features = VariableFeatures(object = ILC))
ElbowPlot(ILC) #dims = 18

#HTO demultiplexing

# Normalize HTO data, here we use centered log-ratio (CLR) transformation
ILC <- NormalizeData(ILC, assay = "HTO", normalization.method = "CLR")

#Run HTO demux
ILC <- HTODemux(ILC, assay = "HTO", positive.quantile = 0.99)

# Global classification results
table(ILC$HTO_classification.global)

# Group cells based on the max HTO signal
Idents(ILC) <- "HTO_maxID"

#Visualize multiplets
RidgePlot(ILC, assay = "HTO", features = rownames(ILC[["HTO"]])[1:2], ncol = 2)
FeatureScatter(ILC, feature1 = "hto_HTO-A", feature2 = "hto_HTO-B")
Idents(ILC) <- "HTO_classification.global"
VlnPlot(ILC, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
VlnPlot(ILC, features = "nFeature_RNA", pt.size = 0.1, log = TRUE)

# First, we will remove negative cells from the object
ILC.subset <- subset(ILC, idents = "Negative", invert = TRUE)

# Calculate a tSNE embedding of the HTO data
DefaultAssay(ILC.subset) <- "HTO"
ILC.subset <- ScaleData(ILC.subset, features = rownames(ILC.subset), 
                        verbose = FALSE)
ILC.subset <- RunPCA(ILC.subset, features = rownames(ILC.subset), approx = FALSE)
ILC.subset <- RunTSNE(ILC.subset, dims = 1:3, perplexity = 100, check_duplicates = FALSE) #check why there are duplicates
DimPlot(ILC.subset)

#We will remove doublets for ILC
ILC.singlet <- subset(ILC, idents = "Doublet", invert = TRUE)

#Filter to exclude cells that have unique feature counts more than 2500 (multiplets) and less than 200 (empty drops) and removing cells that have >5% mitochondrial counts (dead/dying cells)
VlnPlot(ILC.singlet, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ILC.singlet <- subset(ILC.singlet, subset = nFeature_RNA < 2500 & nFeature_RNA > 200 & percent.mt < 5)

#Improve clustering through Weighted nearest neighbor analysis

# Process ADT data and  set a dimensional reduction name to avoid overwriting the RNA PCA
DefaultAssay(ILC.singlet) <- "ADT"
VariableFeatures(ILC.singlet) <- rownames(ILC.singlet[["ADT"]])
ILC.singlet <- NormalizeData(ILC.singlet, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')

ILC.singlet <- FindMultiModalNeighbors(
  ILC.singlet, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:15), modality.weight.name = "RNA.weight"
)

ILC.singlet <- RunUMAP(ILC.singlet, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
ILC.singlet <- FindClusters(ILC.singlet, graph.name = "wsnn", algorithm = 3, resolution = 1, verbose = FALSE)

ILC.plot <- DimPlot(ILC.singlet, reduction = "wnn.umap", label = TRUE, repel = TRUE)

#RNA contribution to each cluster
VlnPlot(ILC.singlet, features = "RNA.weight", sort = TRUE, pt.size = 0.1) +
  NoLegend()

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
library(nichenetr)
library(tidyverse)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

m.s.genes <-convert_human_to_mouse_symbols(s.genes)
m.g2m.genes <- convert_human_to_mouse_symbols(g2m.genes)

ILC.singlet <- CellCycleScoring(ILC.singlet, s.features = m.s.genes, g2m.features = m.g2m.genes, set.ident = TRUE)
# view cell cycle scores and phase assignments
head(ILC.singlet[[]])
#visualize different cell cycle phases
DimPlot(ILC.singlet, reduction = "wnn.umap")
#switch ident back to clusters
Idents(ILC.singlet) <- "old.ident"

#Name clusters
ILC.singlet <- RenameIdents(ILC.singlet, `0` = "ST2 ILC2", `1` = "NK", `2` = "IL18ra ILC2", 
                      `3` = "Act ILC2", `4` = "ILC3", `5` = "Effector CD8", `6` = "NKT17", `7` = "Th1", `8` = "naive CD4", `9` = "naive CD8", 
                      `10` = "Cycling T", `11` = "NKT1", `12` = "DC", `13` = "Mono/Mac",`14` = "Lyve1 IM",`15` = "B cell",`16` = "Neutrophil")

#Re-order clusters to group related cell types
levels(x = ILC.singlet) <- c("B cell","DC","Mono/Mac","Lyve1 IM","Neutrophil","Cycling T","naive CD4","Th1","naive CD8","Effector CD8",
                             "NKT17","NKT1","NK","ST2 ILC2","IL18ra ILC2","Act ILC2","ILC3")

ILC.singlet$celltype <- Idents(ILC.singlet)

#Figure out how to distinguish CD45.1 and CD45.2 cells - CE45.1+ seems to be at least > 0.5
CD45.plot <- FeatureScatter(ILC.singlet, feature1 = "adt_CD45-1", feature2 = "adt_CD45-2")

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
ILC.singlet <- AddMetaData(ILC.singlet, CD45.meta, col.name = "BM")
Idents(ILC.singlet) <- "BM"
DimPlot(ILC.singlet, reduction = "wnn.umap")
Idents(ILC.singlet) <-"celltype"

#add metadata to immune.combined to specify infection status and cluster in one value.
ILC.singlet$celltype.BM <- paste(Idents(ILC.singlet), ILC.singlet$BM, sep = "_")

#Save ILC.singlet file for later recall
saveRDS(ILC.singlet, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/T cell ILC IL-1 scRNA-seq 12-30-20/ILCsinglet")

#load saved ILC.singlet file
ILC.singlet <- readRDS(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/T cell ILC IL-1 scRNA-seq 12-30-20/ILCsinglet")

#Visualize WT vs IL-1R KO gene expression for various clusters
Idents(ILC.singlet) <- "celltype.BM"
ILC3.volcano <- FindMarkers(ILC.singlet, ident.1 = "ILC3_WT", ident.2 = "ILC3_IL-1RKO", logfc.threshold = 0)
ILC3.volcano <- ILC3.volcano[c(1,2)]
p1 <- EnhancedVolcano(ILC3.volcano,
                lab = rownames(ILC3.volcano),
                x = 'avg_log2FC',
                y = 'p_val',
                title = "ILC3 WT vs IL-1R KO",
                subtitle = "",
                pCutoff = 0.05,
                FCcutoff = 0.75,
                xlim = c(-5,5),
                ylim = c(0,10))

NKT17.volcano <- FindMarkers(ILC.singlet, ident.1 = "NKT17_WT", ident.2 = "NKT17_IL-1RKO", logfc.threshold = 0)
NKT17.volcano <- NKT17.volcano[c(1,2)]
p2 <- EnhancedVolcano(NKT17.volcano,
                lab = rownames(NKT17.volcano),
                x = 'avg_log2FC',
                y = 'p_val',
                title = "NKT17 WT vs IL-1R KO",
                subtitle = "",
                pCutoff = 0.05,
                FCcutoff = 0.75,
                xlim = c(-5,5),
                ylim = c(0,10))

ST2.ILC2.volcano <- FindMarkers(ILC.singlet, ident.1 = "ST2 ILC2_WT", ident.2 = "ST2 ILC2_IL-1RKO", logfc.threshold = 0)
ST2.ILC2.volcano <- ST2.ILC2.volcano[c(1,2)]
p3 <- EnhancedVolcano(ST2.ILC2.volcano,
                lab = rownames(ST2.ILC2.volcano),
                x = 'avg_log2FC',
                y = 'p_val',
                title = "ST2 ILC2 WT vs IL-1R KO",
                subtitle = "",
                pCutoff = 0.05,
                FCcutoff = 1)

p3 + p1 + p2

#Label all the ligands by getting potential ligands from NicheNet
potential.ligands <- unique(lr_network$from) %>% convert_human_to_mouse_symbols() %>% na.omit
regulated.ILC3 <- ILC3.volcano %>% filter(avg_log2FC > 1 | avg_log2FC < -1)  %>% row.names()
ILC3.ligands <- intersect(potential.ligands, regulated.ILC3)

regulated.NKT17 <- NKT17.volcano %>% filter(avg_log2FC > 1 | avg_log2FC < -1)  %>% row.names()
NKT17.ligands <- intersect(potential.ligands, regulated.NKT17)

p1 <- EnhancedVolcano(ILC3.volcano,
                      lab = rownames(ILC3.volcano),
                      x = 'avg_log2FC',
                      y = 'p_val',
                      title = "ILC3 WT vs IL-1R KO",
                      subtitle = "",
                      pCutoff = 0.05,
                      FCcutoff = 1,
                      xlim = c(-5,5),
                      ylim = c(0,10),
                      selectLab = ILC3.ligands,
                      drawConnectors = TRUE,
                      arrowheads = FALSE,
                      labSize = 10)

p2 <- EnhancedVolcano(NKT17.volcano,
                      lab = rownames(NKT17.volcano),
                      x = 'avg_log2FC',
                      y = 'p_val',
                      title = "NKT17 WT vs IL-1R KO",
                      subtitle = "",
                      pCutoff = 0.05,
                      FCcutoff =1,
                      xlim = c(-5,5),
                      ylim = c(0,10),
                      selectLab = NKT17.ligands,
                      drawConnectors = TRUE,
                      arrowheads = FALSE,
                      labSize = 10)

#Perform Nichenet on myleoid / granulocyte interactions with ILCs.

#load in myeloid / granulocyte seurat object
mac.neut <- readRDS(file = "/Users/DmitriKotov/Box Sync/Coding stuff/B6 vs Sp140 scRNA-seq 10-9-20/immunecombined2")
mac.neut <- immune.combined

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

#Nichenet analysis on IM interactions with ILC and NKT senders with IM vs AM DEG.

## receiver
receiver = "IM"
expressed_genes_receiver = get_expressed_genes(receiver, assay_oi = "RNA", mac.neut, pct = 0.10)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## sender
sender_celltypes = c("ILC3","NKT17")
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, ILC.singlet, 0.10, assay_oi = "RNA") # lapply to get the expressed genes of every sender cell type separately here
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

cluster.ligand <- DotPlot(ILC.singlet, assay = "RNA", idents = sender_celltypes, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis() + coord_flip()
mocha <- ILC.singlet
mocha$celltype <- factor(x = mocha$celltype, levels = c("NKT17","ILC3","NK","ST2 ILC2","Effector CD8", "Th1","NKT1","Cycling T","naive CD4","IL18ra ILC2","DC","Act ILC2","naive CD8","B cell","Lyve1 IM","Mono/Mac","Neutrophil"))

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
p_receptor_expression = DotPlot(mac.neut, idents = c("IFN IM","IM","Mono","AM") , assay = "RNA", features = rownames(vis_ligand_receptor_network), cols = "RdYlBu") + RotatedAxis()
p_receptor_expression

DE_table_all = Idents(ILC.singlet) %>% levels() %>% intersect(sender_celltypes) %>% lapply(get_lfc_celltype, seurat_obj = ILC.singlet, condition_colname = "BM", condition_oi = "WT", condition_reference = "IL-1RKO", expression_pct = 0.10) %>% reduce(full_join)

DE_table_all[is.na(DE_table_all)] = 0

# Combine ligand activities with DE information
ligand_activities_de = ligand_activities %>% dplyr::select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0

# make LFC heatmap
lfc_matrix = ligand_activities_de  %>% dplyr::select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = lfc_matrix[order_ligands,]

colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()

p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
p_ligand_lfc

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

#Focus on IL-1 regulated genes for NKT17 and ILC3
lfc_ligands = c("Il22","Tnf","Sema4d","Cd28","Csf2","Cd274","Chad","Itgb1","Tgfb1","Icam1","Il17a","Hsp90b1","Arf1","Vegfa","Gpi1")
lfc.cluster.ligand <- DotPlot(ILC.singlet, assay = "RNA", idents = c("ST2 ILC2","NK","NKT17","ILC3"), features = lfc_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis() + coord_flip()

lfc_receptors = c("Il22ra1","Il10rb","Tnfrsf1a","Tnfrsf1b","Met","Cd80","Csf2ra","Csf2rb","Itgb1","Itgax","Tgfbr1","Tgfbr2","Il17ra","Pld2","Nrp1","Nrp2","Amfr")
p_lfc_receptor_expression = DotPlot(mac.neut, idents = c("IFN IM","IM","AM","Mono") , assay = "RNA", features = lfc_receptors, cols = "RdYlBu") + RotatedAxis()
p_lfc_receptor_expression

lfc.cluster.ligand + p_lfc_receptor_expression

#Focus on IL-1 regulated genes for NKT17 and ILC3
lfc_receptor_ligands = c("Il22","Tnf","Sema4d","Csf2","Tgfb1","Icam1","Il17a","Arf1","Vegfa","Gpi1")
lfc.receptor.cluster.ligand <- DotPlot(ILC.singlet, assay = "RNA", idents = c("ST2 ILC2","NK","NKT17","ILC3"), features = lfc_receptor_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis() + coord_flip()

#provide unique names for various NicheNet analysis output for IM recievers and ILC/NK senders with IM vs AM DEG to enable cross reference.
IM.ligand_activities <- ligand_activities
IM.lfc_matrix <- lfc_matrix
IM.lr_network_top_matrix <- lr_network_top_matrix #For finding ligand - receptor interactions

#Repeat IM interaction with ILC and NKT Nichenet analysis, but with IM vs Mono DEG
## receiver
receiver = "IM"
expressed_genes_receiver = get_expressed_genes(receiver, assay_oi = "RNA", mac.neut, pct = 0.10)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

## sender
sender_celltypes = c("ST2 ILC2","IL18ra ILC2","Act ILC2","ILC3","NK","NKT17","NKT1")
list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, ILC.singlet, 0.10, assay_oi = "RNA") # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

#Define a genset of interest - 
DE_table_receiver = FindMarkers(object = mac.neut, assay = "RNA", ident.1 = "IM", ident.2 = "Mono", min.pct = 0.10) %>% rownames_to_column("gene")
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

cluster.ligand <- DotPlot(ILC.singlet, assay = "RNA", idents = sender_celltypes, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis() + coord_flip()

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
DE_table_all = Idents(ILC.singlet) %>% levels() %>% intersect(sender_celltypes) %>% lapply(get_lfc_celltype, seurat_obj = ILC.singlet, condition_colname = "BM", condition_oi = "WT", condition_reference = "IL-1RKO", expression_pct = 0.10) %>% reduce(full_join)

DE_table_all[is.na(DE_table_all)] = 0

# Combine ligand activities with DE information
ligand_activities_de = ligand_activities %>% dplyr::select(test_ligand, pearson) %>% rename(ligand = test_ligand) %>% left_join(DE_table_all %>% rename(ligand = gene))
ligand_activities_de[is.na(ligand_activities_de)] = 0

# make LFC heatmap
lfc_matrix = ligand_activities_de  %>% dplyr::select(-ligand, -pearson) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities_de$ligand)
rownames(lfc_matrix) = rownames(lfc_matrix) %>% make.names()

order_ligands = order_ligands[order_ligands %in% rownames(lfc_matrix)]
vis_ligand_lfc = lfc_matrix[order_ligands,]

colnames(vis_ligand_lfc) = vis_ligand_lfc %>% colnames() %>% make.names()

p_ligand_lfc = vis_ligand_lfc %>% make_threecolor_heatmap_ggplot("Prioritized ligands","LFC in Sender", low_color = "midnightblue",mid_color = "white", mid = median(vis_ligand_lfc), high_color = "red",legend_position = "top", x_axis_position = "top", legend_title = "LFC") + theme(axis.text.y = element_text(face = "italic"))
p_ligand_lfc

#IL-1R1 expression by the sender populations
p_Il1r1_sender <- DotPlot(ILC.singlet, idents = sender_celltypes, features = c("Il1r1")) + coord_flip()
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

p_ligand_pearson + cluster.ligand + p_ligand_lfc + p_ligand_receptor_network + p_receptor_expression + p_Il1r1_sender

#provide unique names for various NicheNet analysis output for IM recievers and ILC/NK senders with IM vs Mono DEG to enable cross reference.
IM.mono.ligand_activities <- ligand_activities
IM.mono.lfc_matrix <- lfc_matrix
IM.mono.lr_network_top_matrix <- lr_network_top_matrix #For finding ligand - receptor interactions

#Plot IL-1 regulated genes expressed by NKT17 and ILC3 cells.
library(tidyverse)
library(nichenetr)

Idents(ILC.singlet) <- "celltype.BM"
NKT17.IL1 <- FindMarkers(ILC.singlet, ident.1 = "NKT17_WT", ident.2 = "NKT17_IL-1RKO")
ILC3.IL1 <- FindMarkers(ILC.singlet, ident.1 = "ILC3_WT", ident.2 = "ILC3_IL-1RKO")

inter_NKT17_IL1 = NKT17.IL1 %>% filter(avg_log2FC > 0.5 & pct.1 > 0.1)
inter_ILC3_IL1 = ILC3.IL1 %>% filter(avg_log2FC > 0.5 & pct.1 > 0.1)

NKT17.ILC3 <- intersect(rownames(inter_NKT17_IL1),rownames(inter_ILC3_IL1))
NKT17.Th17 <- intersect(rownames(inter_NKT17_IL1),rownames(inter_Th17_IL1))
ILC3.Th17 <- intersect(rownames(inter_ILC3_IL1),rownames(inter_Th17_IL1))

NKT17.ILC3.Th17 <- intersect(ILC3.Th17, NKT17.ILC3)

vis_NKT17_IL1 = NKT17.IL1 %>% filter(avg_log2FC > 0.5 & pct.1 > 0.1) %>% select(avg_log2FC) %>% arrange(avg_log2FC) %>% as.matrix()

p_NKT17_IL1 = vis_NKT17_IL1 %>% make_heatmap_ggplot("IL-1 regulated genes","Avg log2 Fold Change", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Avg log2 Fold Change") + theme(legend.text = element_text(size = 9))

#Combine NKT17, ILC3, and Th17 datasets into one big dataset
combined.IL1 <- rbind(inter_NKT17_IL1, inter_ILC3_IL1, inter_Th17_IL1)
length(unique(rownames(combined.IL1))) #returns 389 more conservative, 1111 less conservative
combined.IL1$genes <- rownames(combined.IL1)
with.receptor <- intersect(unique(rownames(combined.IL1)), lr_network$from) #28 with ppi and 23 without

with.receptor.dataset <- inner_join(combined.IL1, lr_network, by = c("genes" = "from")) #56 unique receptors.
length(unique(with.receptor.dataset$to)) #43
