library(Seurat)
library(tidyverse)
library(EnhancedVolcano)

#Download data from the Broad single cell portal
nhp.data <- Read10X(data.dir = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Flynn Shalek CD8 Depletion NHP data/Unprocessed data/", gene.column = 1)

# Initialize the Seurat object with the raw (non-normalized data).
cd8.nhp <- CreateSeuratObject(counts = nhp.data, project = "CD8 NHP")
cd8.nhp #check seurat object

#Metadata
nhp.metadata <- read.csv(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Flynn Shalek CD8 Depletion NHP data/formatted_metadata.csv", header = T)
nhp.metadata <- nhp.metadata[2:41977,]

cd8.nhp <- AddMetaData(object = cd8.nhp, metadata = nhp.metadata, col.name = colnames(nhp.metadata))

# store mitochondrial percentage in object meta data
cd8.nhp <- PercentageFeatureSet(cd8.nhp, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
cd8.nhp <- SCTransform(cd8.nhp, vars.to.regress = "percent.mt", verbose = FALSE)
cd8.nhp <- RunPCA(cd8.nhp, verbose = FALSE)
cd8.nhp <- RunUMAP(cd8.nhp, dims = 1:30, verbose = FALSE)

cd8.nhp <- FindNeighbors(cd8.nhp, dims = 1:30, verbose = FALSE)
cd8.nhp <- FindClusters(cd8.nhp, verbose = FALSE)
DimPlot(cd8.nhp, label = TRUE)

saveRDS(cd8.nhp, file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Flynn Shalek CD8 Depletion NHP data/nhp.rds")

FeaturePlot(cd8.nhp, features = c("CD163","MRC1","ITGAM","ITGAX","CD14","CD68"), order = T)

#Annotating the dataset
Idents(cd8.nhp) <- "treatment"
cd8.nhp.igg <- subset(cd8.nhp, idents = "IgG")
FeaturePlot(cd8.nhp.igg, features = c("TREM2","IL1RN","SPP1"))
Idents(cd8.nhp.igg) <- "seurat_clusters"
DimPlot(cd8.nhp.igg, label = T)
cd8.nhp.igg <- RenameIdents(cd8.nhp.igg, `0` = "TREM2 IM",`1` = "T cell", `2` = "AM",`3` = "Proliferating T cells", `4` = "SPP1 IM", `5` = "Monocytes",`6` = "Mast",`7` = "T cell",`8` = "NK",`9` = "T cell", `10` = "IM", `11` = "TREM2 IM", `12` = "TREM2 IM", `13` = "B Cells", `14` = "cDC2",`15` = "NKT", 
                            `16` = "T cell",`17` = "Plasma Cells", `18` = "T2P",`19` = "IM", `20` = "cDC1", `21` = "Endothelial", `22` = "Proliferating Macs", `23` = "Fibroblasts",`24` = "T cell", `25` = "pDC", `26` = "Proliferating Mast",`27` = "Neutrophils", `28` = "T1P")

#Compare SPP1 IM and TREM2 IM
nhp.SPP1.TREM2.genes <- FindMarkers(cd8.nhp.igg, ident.1 = "SPP1 IM", ident.2 = "TREM2 IM")
nhp.SPP1.TREM2.volcano <- nhp.SPP1.TREM2.genes[c(2,5)]
EnhancedVolcano(nhp.SPP1.TREM2.volcano,
                lab = rownames(nhp.SPP1.TREM2.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "SPP1 IM vs TREM2 IM",
                subtitle = "",
                drawConnectors = T,
                arrowheads = F,
                raster = T)

nhp.SPP1.TREM2.genes <- FindMarkers(cd8.nhp.igg, ident.1 = "SPP1 IM", ident.2 = "TREM2 IM")
nhp.SPP1.TREM2.volcano <- nhp.SPP1.TREM2.genes[c(2,5)]
EnhancedVolcano(nhp.SPP1.TREM2.volcano,
                lab = rownames(nhp.SPP1.TREM2.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "SPP1 IM vs TREM2 IM",
                subtitle = "",
                drawConnectors = T,
                arrowheads = F,
                selectLab = c("SPP1","TREM2","C1QC","C1QB","AXL","IDO1","IL18BP","CD274","IL6","IL1RN","EGLN3","NDRG1","CLEC6A","ADAM8","CHST11", "MMP2","APOC2","CD101","LY86","GPR155","ITGB5","ITGA9","GLO1","CCL18"),
                raster = T)

nhp.SPP1.TREM2.genes$genes <- row.names(nhp.SPP1.TREM2.genes)

EnhancedVolcano(nhp.SPP1.TREM2.volcano,
                lab = rownames(nhp.SPP1.TREM2.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "SPP1 IM vs TREM2 IM",
                subtitle = "",
                drawConnectors = T,
                arrowheads = F,
                selectLab = c("SPP1","TREM2","IL1RN","C1QC","INHBA","AXL"),
                raster = T)

nhp.IL1RN.TREM2.genes <- FindMarkers(cd8.nhp, ident.1 = "IL1RN IM", ident.2 = "TREM2 IM")
nhp.IL1RN.TREM2.volcano <- nhp.IL1RN.TREM2.genes[c(2,5)]
EnhancedVolcano(nhp.IL1RN.TREM2.volcano,
                lab = rownames(nhp.IL1RN.TREM2.volcano),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = "IL1RN IM vs TREM2 IM",
                subtitle = "",
                drawConnectors = T,
                arrowheads = F,
                raster = T)

#Scatter plot for myeloid cells
cd8.nhp.igg$celltypes <- Idents(cd8.nhp.igg)
Idents(cd8.nhp.igg) <- "celltypes"
nhp.mac <- subset(cd8.nhp.igg, idents = c("TREM2 IM","AM","SPP1 IM","Monocytes","IM"))
FeatureScatter(nhp.mac, feature1 = "SPP1", feature2 = "IL1RN")
DimPlot(nhp.mac)
FeaturePlot(nhp.mac, features = c("SPP1", "IL1RN"), blend = TRUE)


#3 cells clustered weirdly and look like background noise so they were removed from the analysis.
plot <- DimPlot(nhp.mac)
select.cells <- CellSelector(plot = plot) #"Array6_29618_CGGTGTGGCAGC" "Array2_29818_ACGATTGGAACT" "Array2_29818_TGGCTTGGAGAT"
select.cells <- c("Array6_29618_CGGTGTGGCAGC","Array2_29818_ACGATTGGAACT","Array2_29818_TGGCTTGGAGAT")
nhp.mac_filtered <- nhp.mac[,!colnames(nhp.mac) %in% select.cells]
DimPlot(nhp.mac_filtered)

