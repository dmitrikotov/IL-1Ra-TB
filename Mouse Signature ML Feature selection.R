library(tidyverse)
library(Seurat)
library(caret)
library(hacksig)
library(pROC)
library(nichenetr)

# ==============================================================================
# 1. Filter scRNA-seq Differentially Expressed Genes (DEGs)
# ==============================================================================
# Input: DEG results table from scRNA-seq (e.g., Seurat FindMarkers output)
super.integrated <- readRDS(file = '/Users/dmitrikotov/Library/CloudStorage/Box-Box/Dmitri Personal/DK Postdoc Data and Analysis/Coding stuff/Irg1 and iNOS scRNAseq 220315/b6 sp140 irg1 inos v2')
Idents(super.integrated) <- "cell_type.infected"
b6 <- FindMarkers(super.integrated, ident.1 = "B6_infected", ident.2 = "B6_naive")
sp140 <- FindMarkers(super.integrated, ident.1 = "Sp140_infected", ident.2 = "Sp140_naive")
irg1 <- FindMarkers(super.integrated, ident.1 = "Irg1_infected", ident.2 = "Irg1_naive")
nos2 <- FindMarkers(super.integrated, ident.1 = "iNOS_infected", ident.2 = "iNOS_naive")

b6$gene <- row.names(b6)
sp140$gene <- row.names(sp140)
irg1$gene <- row.names(irg1)
nos2$gene <- row.names(nos2)

# Filter by significance (adj p < 0.05) and effect size (|log2FC| > 1.0)
#tried changing abs(avg_log2FC) to just avg_log2FC
filtered_b6 <- b6 %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 1.0) %>%
  mutate(human_gene = convert_mouse_to_human_symbols(gene)) %>%
  filter(!is.na(human_gene)) %>% pull(human_gene)

cat("Candidate scRNA-seq DEGs passed filter:", length(filtered_b6), "\n")

filtered_sp140 <- sp140 %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 1.0) %>%
  mutate(human_gene = convert_mouse_to_human_symbols(gene)) %>%
  filter(!is.na(human_gene)) %>% pull(human_gene)

cat("Candidate scRNA-seq DEGs passed filter:", length(filtered_sp140), "\n")

filtered_irg1 <- irg1 %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 1.0) %>%
  mutate(human_gene = convert_mouse_to_human_symbols(gene)) %>%
  filter(!is.na(human_gene)) %>% pull(human_gene)

cat("Candidate scRNA-seq DEGs passed filter:", length(filtered_irg1), "\n")

filtered_nos2 <- nos2 %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 1.0) %>%
  mutate(human_gene = convert_mouse_to_human_symbols(gene)) %>%
  filter(!is.na(human_gene)) %>% pull(human_gene)

cat("Candidate scRNA-seq DEGs passed filter:", length(filtered_nos2), "\n")

# ==============================================================================
# 2. Feature Selection & Model Training with caret
# ==============================================================================
# Load primary dataset for model training (Samples in rows, Genes + 'class' in columns)
tb.counts.norm <- readRDS(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Dmitri Personal/DK Postdoc Data and Analysis/Coding stuff/Human TB Gene Signatures/tb_counts_norm")

TB.infected <- readRDS(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Dmitri Personal/DK Postdoc Data and Analysis/Coding stuff/TB_infected")

tb.counts.norm_transposed <- tb.counts.norm %>%
  # 1. Convert existing rownames into an explicit column
  tibble::rownames_to_column(var = "row_names") %>%

  # 2. Collapse all data columns into long format
  pivot_longer(cols = -row_names, names_to = "variable", values_to = "value") %>%

  # 3. Re-expand the data spreading the original row names across columns
  pivot_wider(names_from = row_names, values_from = value)

tb.counts.norm_transposed <- column_to_rownames(tb.counts.norm_transposed, var = "variable")

# Independent test expression dataset (Genes in rows, Samples in columns)
berry_southafrica_tb <- readRDS(file = "/Users/dmitrikotov/Library/CloudStorage/Box-Box/Dmitri Personal/DK Postdoc Data and Analysis/Coding stuff/Human TB Gene Signatures/berry_southafrica_tb")
berry_southafrica_norm <-readRDS(file = '/Users/dmitrikotov/Library/CloudStorage/Box-Box/Dmitri Personal/DK Postdoc Data and Analysis/Coding stuff/Human TB Gene Signatures/berry_southafrica_norm')

berry_meta <- data.frame("sample_id" = colnames(berry_southafrica_norm))
berry_meta <- berry_meta %>% mutate(label = if_else(sample_id %in% berry_southafrica_tb, "TB", "Control"))

#Identify genes shared between training and test dataset
common_genes <- intersect(row.names(berry_southafrica_norm), colnames(tb.counts.norm_transposed))
tb.counts.norm_transposed <- tb.counts.norm_transposed[, common_genes]
berry_southafrica_norm <- berry_southafrica_norm[common_genes,]

# Add the new column based on vector matching
tb.counts.norm_transposed <- tb.counts.norm_transposed %>%
  mutate(condition = if_else(row.names(tb.counts.norm_transposed) %in% TB.infected, "TB", "Control"))

# Intersect filtered DEGs with available features in training data
candidate_features_b6 <- intersect(filtered_b6, colnames(tb.counts.norm_transposed))
X_train_b6 <- tb.counts.norm_transposed[, candidate_features_b6]

candidate_features_sp140 <- intersect(filtered_sp140, colnames(tb.counts.norm_transposed))
X_train_sp140 <- tb.counts.norm_transposed[, candidate_features_sp140]

candidate_features_irg1 <- intersect(filtered_irg1, colnames(tb.counts.norm_transposed))
X_train_irg1 <- tb.counts.norm_transposed[, candidate_features_irg1]

candidate_features_nos2 <- intersect(filtered_nos2, colnames(tb.counts.norm_transposed))
X_train_nos2 <- tb.counts.norm_transposed[, candidate_features_nos2]

y_train <- factor(tb.counts.norm_transposed$condition, levels = c("Control", "TB"))

# Configure Recursive Feature Elimination (RFE) for feature selection
rfe_ctrl <- rfeControl(
  functions = rfFuncs,
  method = "repeatedcv",
  number = 5,
  repeats = 3,
  verbose = FALSE
)

set.seed(42)
rfe_fit_b6 <- rfe(
  x = X_train_b6,
  y = y_train,
  sizes = c(5, 10, 15, 20, 25),
  rfeControl = rfe_ctrl
)

rfe_fit_sp140 <- rfe(
  x = X_train_sp140,
  y = y_train,
  sizes = c(5, 10, 15, 20, 25),
  rfeControl = rfe_ctrl
)

rfe_fit_nos2 <- rfe(
  x = X_train_nos2,
  y = y_train,
  sizes = c(5, 10, 15, 20, 25),
  rfeControl = rfe_ctrl
)

rfe_fit_irg1 <- rfe(
  x = X_train_irg1,
  y = y_train,
  sizes = c(5, 10, 15, 20, 25),
  rfeControl = rfe_ctrl
)

# Extract optimal selected gene signature
selected_genes_b6 <- predictors(rfe_fit_b6)
cat("Selected signature genes:", paste(selected_genes_b6, collapse = ", "), "\n")

selected_genes_sp140 <- predictors(rfe_fit_sp140)
cat("Selected signature genes:", paste(selected_genes_sp140, collapse = ", "), "\n")

selected_genes_irg1 <- predictors(rfe_fit_irg1)
cat("Selected signature genes:", paste(selected_genes_irg1, collapse = ", "), "\n")

selected_genes_nos2 <- predictors(rfe_fit_nos2)
cat("Selected signature genes:", paste(selected_genes_nos2, collapse = ", "), "\n")

# Fit final Random Forest model on selected features
train_ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary
)

final_model_b6 <- caret::train(
  x = X_train_b6[, selected_genes_b6],
  y = y_train,
  method = "rf",
  metric = "ROC",
  trControl = train_ctrl
)

final_model_sp140 <- caret::train(
  x = X_train_sp140[, selected_genes_sp140],
  y = y_train,
  method = "rf",
  metric = "ROC",
  trControl = train_ctrl
)


final_model_irg1 <- caret::train(
  x = X_train_irg1[, selected_genes_irg1],
  y = y_train,
  method = "rf",
  metric = "ROC",
  trControl = train_ctrl
)

final_model_nos2 <- caret::train(
  x = X_train_nos2[, selected_genes_nos2],
  y = y_train,
  method = "rf",
  metric = "ROC",
  trControl = train_ctrl
)

#Features used in final model
b6_final_features <- varImp(final_model_b6, scale = FALSE)
print(b6_final_features)

sp140_final_features <- varImp(final_model_sp140, scale = FALSE)
print(sp140_final_features)

irg1_final_features <- varImp(final_model_irg1, scale = FALSE)
print(irg1_final_features)

nos2_final_features <- varImp(final_model_nos2, scale = FALSE)
print(nos2_final_features)


# ==============================================================================
# 3. AUC & Performance Evaluation using pROC
# ==============================================================================

#Evaluate caret Machine Learning Classifier Predictions
X_test_b6 <- t(berry_southafrica_norm[selected_genes_b6, ])
test_preds_prob_b6 <- predict(final_model_b6, newdata = X_test_b6, type = "prob")

roc_ml_b6 <- roc(
  response = factor(berry_meta$label, levels = c("Control", "TB")),
  predictor = test_preds_prob_b6$TB,
  levels = c("Control", "TB"),
  direction = "<"
)

auc_ml_b6 <- auc(roc_ml_b6)
ci_ml_b6  <- ci.auc(roc_ml_b6)

cat(sprintf("B6 caret ML Classifier AUC: %.3f (95%% CI: %.3f - %.3f)\n",
            auc_ml_b6, ci_ml_b6[1], ci_ml_b6[3]))

X_test_sp140 <- t(berry_southafrica_norm[selected_genes_sp140, ])
test_preds_prob_sp140 <- predict(final_model_sp140, newdata = X_test_sp140, type = "prob")

roc_ml_sp140 <- roc(
  response = factor(berry_meta$label, levels = c("Control", "TB")),
  predictor = test_preds_prob_sp140$TB,
  levels = c("Control", "TB"),
  direction = "<"
)

auc_ml_sp140 <- auc(roc_ml_sp140)
ci_ml_sp140  <- ci.auc(roc_ml_sp140)

cat(sprintf("Sp140 caret ML Classifier AUC: %.3f (95%% CI: %.3f - %.3f)\n",
            auc_ml_sp140, ci_ml_sp140[1], ci_ml_sp140[3]))

X_test_irg1 <- t(berry_southafrica_norm[selected_genes_irg1, ])
test_preds_prob_irg1 <- predict(final_model_irg1, newdata = X_test_irg1, type = "prob")

roc_ml_irg1 <- roc(
  response = factor(berry_meta$label, levels = c("Control", "TB")),
  predictor = test_preds_prob_irg1$TB,
  levels = c("Control", "TB"),
  direction = "<"
)

auc_ml_irg1 <- auc(roc_ml_irg1)
ci_ml_irg1  <- ci.auc(roc_ml_irg1)

cat(sprintf("Irg1 caret ML Classifier AUC: %.3f (95%% CI: %.3f - %.3f)\n",
            auc_ml_irg1, ci_ml_irg1[1], ci_ml_irg1[3]))

X_test_nos2 <- t(berry_southafrica_norm[selected_genes_nos2, ])
test_preds_prob_nos2 <- predict(final_model_nos2, newdata = X_test_nos2, type = "prob")

roc_ml_nos2 <- roc(
  response = factor(berry_meta$label, levels = c("Control", "TB")),
  predictor = test_preds_prob_nos2$TB,
  levels = c("Control", "TB"),
  direction = "<"
)

auc_ml_nos2 <- auc(roc_ml_nos2)
ci_ml_nos2  <- ci.auc(roc_ml_nos2)

cat(sprintf("Nos2 caret ML Classifier AUC: %.3f (95%% CI: %.3f - %.3f)\n",
            auc_ml_nos2, ci_ml_nos2[1], ci_ml_nos2[3]))


#overview
cat("B6 Candidate scRNA-seq DEGs passed filter:", length(filtered_b6), "\n")
cat("B6 Selected signature genes:", paste(selected_genes_b6, collapse = ", "), "\n")
cat(sprintf("B6 caret ML Classifier AUC: %.3f (95%% CI: %.3f - %.3f)\n",
            auc_ml_b6, ci_ml_b6[1], ci_ml_b6[3]))

cat("Sp140 Candidate scRNA-seq DEGs passed filter:", length(filtered_sp140), "\n")
cat("Sp140 Selected signature genes:", paste(selected_genes_sp140, collapse = ", "), "\n")
cat(sprintf("Sp140 caret ML Classifier AUC: %.3f (95%% CI: %.3f - %.3f)\n",
            auc_ml_sp140, ci_ml_sp140[1], ci_ml_sp140[3]))

cat("Irg1 Candidate scRNA-seq DEGs passed filter:", length(filtered_irg1), "\n")
cat("Irg1 Selected signature genes:", paste(selected_genes_irg1, collapse = ", "), "\n")
cat(sprintf("Irg1 caret ML Classifier AUC: %.3f (95%% CI: %.3f - %.3f)\n",
            auc_ml_irg1, ci_ml_irg1[1], ci_ml_irg1[3]))

cat("Nos2 Candidate scRNA-seq DEGs passed filter:", length(filtered_nos2), "\n")
cat("Nos2 Selected signature genes:", paste(selected_genes_nos2, collapse = ", "), "\n")
cat(sprintf("Nos2 caret ML Classifier AUC: %.3f (95%% CI: %.3f - %.3f)\n",
            auc_ml_nos2, ci_ml_nos2[1], ci_ml_nos2[3]))

#Genes used in final model:
model_genes_b6 <- predictors(final_model_b6)
cat("B6 Final model gene count:", length(model_genes_b6), "\n")
print(model_genes_b6)

model_genes_sp140 <- predictors(final_model_sp140)
cat("Sp140 Final model gene count:", length(model_genes_sp140), "\n")
print(model_genes_sp140)

model_genes_irg1 <- predictors(final_model_irg1)
cat("Irg1 Final model gene count:", length(model_genes_irg1), "\n")
print(model_genes_irg1)

model_genes_nos2 <- predictors(final_model_nos2)
cat("Nos2 Final model gene count:", length(model_genes_nos2), "\n")
print(model_genes_nos2)

# Plot comparative ROC curves
pROC::plot.roc(roc_ml_b6, print.auc = TRUE, col = "red", print.auc.y =0.47)
pROC::plot.roc(roc_ml_sp140, add=TRUE, print.auc = TRUE, col = "blue",  print.auc.y =0.4)
pROC::plot.roc(roc_ml_irg1,add = TRUE, print.auc = TRUE, col = "chartreuse", print.auc.y =0.33)
pROC::plot.roc(roc_ml_nos2, add=TRUE, print.auc = TRUE, col = "cadetblue2",  print.auc.y =0.27)

#roc.test using Boostrap test
roc.test(roc_ml_sp140, roc_ml_b6) #p-value = 3.061e-06
roc.test(roc_ml_irg1, roc_ml_b6) #p-value = 4.049e-06
roc.test(roc_ml_nos2, roc_ml_b6) #p-value = 2.439e-05
