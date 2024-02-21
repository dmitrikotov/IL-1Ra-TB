library(Seurat)
library(tidyverse)
library(patchwork)
library(nichenetr)

super.integrated <- readRDS(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Coding stuff/Irg1 and iNOS scRNAseq 220315/b6 sp140 irg1 inos v2")

#Nichenet setup
options(timeout=600)
organism = "mouse"

if(organism == "human"){
  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
  weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
} else if(organism == "mouse"){
  lr_network = readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
  ligand_target_matrix = readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
  weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
  
}

lr_network = lr_network %>% distinct(from, to)
head(lr_network)
ligand_target_matrix[1:5,1:5]

weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network, by = c("from","to"))
head(weighted_networks$lr_sig)


#Define ligands regulated by IL-1
il1.ligands.myofibroblast <- c("Icam1","Il6","Ccl2","Wnt11")
il1.ligands.fibroblast <- c("Cxcl1","Il6","Ccl2","Cxcl9")
il1.ligands.endothelial <- c("Icam1")
il1.ligands.ilc3 <- c("Aimp1","Cxcl2")
il1.ligands.nkt17 <- c("Tnf","Cxcl2","Il17a","Itgb2","Rtn4")
il1.ligands.th17 <- c("Tnf","Gzmb","Ifng")
il1.ligands <- c(il1.ligands.myofibroblast, il1.ligands.fibroblast, il1.ligands.endothelial, il1.ligands.ilc3, il1.ligands.nkt17, il1.ligands.th17)
il1.ligands <- unique(il1.ligands)

#Simplify the integrated data to find the genes expressed by the receiver cells
Idents(super.integrated) <- "state"
simple.integrated <- subset(super.integrated, idents = c("Mtb+","Mtb-"))
Idents(simple.integrated) <- "celltype"
simple.integrated <- RenameIdents(simple.integrated, 'ISG Mature Neut' = "Neutrophil", 'ISG Activated Neut' = "Neutrophil", 'Aged Neut' = "Neutrophil", 'Mature Neut' = "Neutrophil",
                                  'ISG Old Neut' = "Neutrophil", 'ISG Young Neut' = "Neutrophil", 'CD63 Neut' = "Neutrophil", 'Mmp8 Neut' = "Neutrophil",'Nos2 Neut' = "Neutrophil",
                                  'Camp Mmp8 Neut' = "Neutrophil", 'CD16-2 Mono' = "Mono", 'ISG IM' = "IM", 'Activated IM' = "IM", 'Activated AM' = "AM")

#Define reciever expressed genes
receiver = "IM"
expressed_genes_receiver = get_expressed_genes(receiver, assay_oi = "RNA", simple.integrated, pct = 0.10)
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

#get the expressed ligand-receptor pairs
ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands,il1.ligands)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)

lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)

potential_ligands = lr_network_expressed %>% pull(from) %>% unique()
head(potential_ligands)

#Define a geneset of interest - comparing Activated IM (Trem2+ IM) vs IM
Idents(super.integrated) <- "celltype"
DE_table_receiver = FindMarkers(object = super.integrated, assay = "RNA", ident.1 = "Activated IM", ident.2 = "IM", min.pct = 0.10) %>% rownames_to_column("gene")
geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 1) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#Preform NicheNet analysis to define ligand activity - TREM2+ IM
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities %>% arrange(-aupr_corrected)
ligand_activities_activated <- ligand_activities

# show histogram of ligand activity scores - TREM2+ IM
p_hist_lig_activity = ggplot(ligand_activities, aes(x=aupr_corrected)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(30, aupr_corrected) %>% pull(aupr_corrected))), color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity
p_hist_lig_activity_Activated <- p_hist_lig_activity
p_hist_lig_activity_Activated

#active target inference - TREM2 IM
best_upstream_ligands = ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)
head(best_upstream_ligands)

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)
nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network_Activated = vis_ligand_target %>% make_heatmap_ggplot("IL-1 Regulated Ligands","Genes upregulated in TREM2+ IM vs IM", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))

p_ligand_target_network_Activated

# get the ligand-receptor network of the top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

# get the weights of the ligand-receptor interactions as used in the NicheNet model
weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

# convert to a matrix
lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

# perform hierarchical clustering to order the ligands and receptors
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_ligands_receptor <- rev(best_upstream_ligands)
vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("IL-1 Regulated Ligands","Receptors expressed by IMs", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network

#Show Ligand Activity - TREM2 IM
ligand_aupr_matrix = ligand_activities %>% select(aupr_corrected) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
ligand_aupr_matrix = ligand_activities_activated %>% select(aupr_corrected) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected)
order_ligands = intersect(ligand_activities$test_ligand, colnames(active_ligand_target_links)) %>% rev()

vis_ligand_aupr = ligand_aupr_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("AUPR")
p_ligand_aupr = vis_ligand_aupr %>% make_heatmap_ggplot("IL-1 Regulated Ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "AUPR\n(target gene prediction ability)")
p_ligand_aupr_Activated <- p_ligand_aupr


p_ligand_aupr_Activated


#Redo the analysis as above but comparing Spp1+ IM to IM


#Define a geneset of interest - comparing ISG IM (Spp1+ IM) vs IM
DE_table_receiver = FindMarkers(object = super.integrated, assay = "RNA", ident.1 = "ISG IM", ident.2 = "IM", min.pct = 0.10) %>% rownames_to_column("gene")
geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 1) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

#Preform NicheNet analysis to define ligand activity - SPP1+ IM
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities %>% arrange(-aupr_corrected)
ligand_activities_ISG <- ligand_activities

# show histogram of ligand activity scores - TREM2+ IM
p_hist_lig_activity = ggplot(ligand_activities, aes(x=aupr_corrected)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(30, aupr_corrected) %>% pull(aupr_corrected))), color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity
p_hist_lig_activity_ISG <- p_hist_lig_activity

p_hist_lig_activity_ISG

#active target inference - SPP1 IM
best_upstream_ligands = ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)
head(best_upstream_ligands)

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)
nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network_ISG = vis_ligand_target %>% make_heatmap_ggplot("IL-1 Regulated Ligands","Genes upregulated in SPP1+ IM vs IM", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential") + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.005,0.01)) + theme(axis.text.x = element_text(face = "italic"))

p_ligand_target_network_ISG

# get the ligand-receptor network of the top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

# get the weights of the ligand-receptor interactions as used in the NicheNet model
weighted_networks = readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

# convert to a matrix
lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

# perform hierarchical clustering to order the ligands and receptors
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_ligands_receptor <- rev(best_upstream_ligands)
vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("IL-1 Regulated Ligands","Receptors expressed by IMs", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network

#Show Ligand Activity - SPP1 IM
ligand_aupr_matrix = ligand_activities %>% select(aupr_corrected) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)
ligand_aupr_matrix = ligand_activities_activated %>% select(aupr_corrected) %>% as.matrix() %>% magrittr::set_rownames(ligand_activities$test_ligand)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected)
order_ligands = intersect(ligand_activities$test_ligand, colnames(active_ligand_target_links)) %>% rev()

vis_ligand_aupr = ligand_aupr_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("AUPR")
p_ligand_aupr = vis_ligand_aupr %>% make_heatmap_ggplot("IL-1 Regulated Ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "AUPR\n(target gene prediction ability)")
p_ligand_aupr_ISG <- p_ligand_aupr

p_ligand_aupr_ISG