library(glmnet); library(randomForest); library(xgboost); library(e1071)
source("msvmRFE.R")

# LASSO
set.seed(3)
x <- as.matrix(data_roc[, -1]); y <- as.matrix(data_roc[, 1])
cvfit <- cv.glmnet(x, y, family = "binomial")
coef_min <- coef(cvfit, s = cvfit$lambda.min)
lassoGene <- setdiff(rownames(coef_min)[which(coef_min != 0)], "(Intercept)")
lassoGene

# RF
set.seed(123)
rf  <- randomForest(as.factor(group) ~ ., data = data_roc, ntree = 500, importance = TRUE)
rf2 <- randomForest(as.factor(group) ~ ., data = data_roc,
                    ntree = which.min(rf$err.rate[, 1]), importance = TRUE)
imp_rf <- importance(rf2)[, "MeanDecreaseGini"]
rfGenes <- names(sort(imp_rf[imp_rf > 3], decreasing = TRUE))
rfGenes

# SVM-RFE
set.seed(12345)
input <- data_roc
nfold <- 10
folds <- rep(1:nfold, len = nrow(input))[sample(nrow(input))]
folds <- lapply(1:nfold, function(x) which(folds == x))
results <- lapply(folds, svmRFE.wrap, input, k = 10, halve.above = 100)
top.features <- WriteFeatures(results, input, save = FALSE) %>% arrange(AvgRank)
featsweep <- lapply(1:(ncol(input) - 1), FeatSweep.wrap, results, input)
errors  <- sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))
top <- top.features[1:which.min(errors), "FeatureName"]
top

# XGBoost
set.seed(123)
labels <- as.numeric(factor(data_roc$group)) - 1
dtrain <- xgb.DMatrix(data = as.matrix(data_roc[, -1]), label = labels)
model  <- xgb.train(list(objective = "binary:logistic", 
                         eval_metric = "logloss",
                         max_depth = 6, 
                         eta = 0.3, 
                         nthread = 4), 
                    data = dtrain, nrounds = 50)
imp_xgb <- xgb.importance(feature_names = colnames(data_roc)[-1], model = model)
imp_xgb


# Single cell
library(harmony)
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(ggsci)
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(pheatmap)
library(scplotter)
library(RColorBrewer)
clusterCols <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488",
                 "#F39B7F", "#8491B4")

sample_PATH <- list.files("./GSE196638/data/")
sample <- c("NL1","NL2","NL3","TXP24", "TXP25", "TXP26")

dataList <- list()
for (i in 1:length(sample_PATH)) {
  print(sample_PATH[[i]])
  x = Read10X(data.dir = paste0("./GSE196638/data/", sample_PATH[[i]]))
  dataList[[i]] <- CreateSeuratObject(counts = x, 
                                      project = paste0(sample[i]),
                                      min.cells = 3, 
                                      min.features = 200)
  dataList[[i]][["percent.mt"]] <- PercentageFeatureSet(dataList[[i]], pattern = "^MT-")
}

data_merge <- merge(dataList[[1]], y = dataList[-1], add.cell.ids = sample)

p1 <- VlnPlot(object = data_merge,
              features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
              cols = clusterCols,
              ncol = 3, group.by = "orig.ident", pt.size = 0) +
  theme_bw() + theme(panel.grid = element_blank());p1

dat_sub <- subset(data_merge, subset = nCount_RNA > 1000 & 
                    nCount_RNA < 25000 &
                    nFeature_RNA > 300 & 
                    nFeature_RNA < 5000 & 
                    percent.mt < 10)

p2 <- VlnPlot(object = dat_sub,
              features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
              cols = clusterCols,
              ncol = 3, group.by = "orig.ident", pt.size = 0) +
  theme_bw() + theme(panel.grid = element_blank());p2

data_sub <- NormalizeData(dat_sub, 
                          normalization.method = "LogNormalize",
                          scale.factor = 10000)
data_sub <- FindVariableFeatures(data_sub, 
                                 selection.method = "vst", 
                                 nfeatures = 2000)
data_sub <- ScaleData(data_sub, features = VariableFeatures(data_sub))
data_sub <- Seurat::RunPCA(data_sub, features = VariableFeatures(object = data_sub))

ElbowPlot(data_sub, ndims = 50)
dat_harmony <- RunHarmony(data_sub, 
                          group.by.vars = "orig.ident",
                          max_iter = 20, 
                          reduction.save = "harmony")
ElbowPlot(dat_harmony, ndims = 50, reduction = "harmony")
dat_harmony <- FindNeighbors(dat_harmony, reduction = "harmony", dims = 1:20)

dat_harmony <- FindClusters(dat_harmony, resolution = 0.2, verbose = FALSE)
dat_harmony <- RunUMAP(dat_harmony, reduction = "harmony", dims = 1:20)
dat_harmony <- RunTSNE(dat_harmony, reduction = "harmony", dims = 1:20)
Idents(dat_harmony) <- "RNA_snn_res.0.2"

marker_list <- list(
  "CD4+ T cell"    = c("CD3D", "CD3E", "CD4", "IL7R", "CD40LG", "TRAC"),
  "NK/CD8+ T cell" = c("CD8A", "NKG7", "GNLY", "PRF1", "KLRD1"),
  "B cell"        = c("MS4A1", "CD79A", "BANK1"),
  "Plasma cell"   = c("MZB1", "JCHAIN", "DERL3", "XBP1"),
  "Monocyte"     = c("FCN1", "S100A8", "S100A9", "VCAN", "S100A12", "LST1", "FCGR3A"),
  "Macrophage"    = c("C1QA", "CD68", "CD163", "MSR1", "MARCO", "FOLR2", "MRC1"),
  "pDC"      = c("CLEC4C", "LILRA4", "IRF7", "PTCRA"),
  "Mast cell"     = c("TPSAB1", "CPA3", "MS4A2", "HDC"),
  "Endothelial cell"    = c("PECAM1", "CDH5", "CLDN5", "VWF", "KDR"),
  "Fibroblast"    = c("DCN", "LUM", "COL1A1", "FBLN1", "PDGFRA"),
  "Smooth muscle cell"      = c("ACTA2", "MYH11", "NOTCH3", "TAGLN", "DES"),
  "AT1"      = c("AGER", "PDPN", "CAV1", "RTKN2"),
  "AT2"      = c("ABCA3", "LAMP3", "MUC1", "PIGR"),
  "Basal"    = c("KRT15", "KRT17", "TRIM29", "DSC3", "TP63"))

marker_list_flt <- lapply(marker_list, function(g) g[g %in% rownames(dat_harmony)])
marker_list_flt <- marker_list_flt[lengths(marker_list_flt) > 0]
marker_list_flt <- marker_list_flt[order(names(marker_list_flt))]

pal <- rev(RColorBrewer::brewer.pal(n = 8, name = "RdBu"))
Idents(dat_harmony) <- "RNA_snn_res.0.2"
p_marker <- DotPlot(dat_harmony, 
                    features = marker_list_flt,
                    assay = "RNA") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
        axis.text.y = element_text(colour = "black"),
        text = element_text(colour = "black", size = 11),
        strip.text.x = element_text(size = 9, margin = margin(2, 1, 2, 1)),
        panel.grid = element_blank()) +
  scale_colour_gradientn(colors = pal) +
  labs(x = "", y = "");p_marker

new.cluster.ids <- c(
  "0"="Monocyte", 
  "1"="CD4+ T cell",
  "2"="NK/CD8+ T cell",
  "3"="Fibroblast",
  "4"="Macrophage",
  "5"="Endothelial cell",
  "6"="AT2",
  "7"="Endothelial cell",
  "8"="B cell",
  "9"="Mast cell",
  "10"="Smooth muscle cell",
  "11"="Plasma cell",
  "12"="Endothelial cell",
  "13"="AT1",
  "14"="Basal",
  "15"="Unkown",
  "16"="Unkown",
  "17"="pDC",
  "18"="Unkown")              

Idents(dat_harmony) <- "RNA_snn_res.0.2"
dat_harmony$celltype <- factor(unname(new.cluster.ids[as.character(dat_harmony$seurat_clusters)]))
stopifnot(!any(is.na(dat_harmony$celltype)))


keep_cells <- colnames(dat_harmony)[dat_harmony$celltype != "Unkown"]
dat_harmony <- subset(dat_harmony, cells = keep_cells)

dat_harmony$celltype <- droplevels(dat_harmony$celltype)
table(dat_harmony$celltype)
stopifnot(!"Unkown" %in% levels(dat_harmony$celltype))

Idents(dat_harmony) <- "celltype"

dat_harmony$group <- "Normal"
dat_harmony$group[dat_harmony$orig.ident %in% c("TXP24", "TXP25", "TXP26")] <- "COPD"

celltype_order <- c("CD4+ T cell", 
                    "NK/CD8+ T cell", 
                    "B cell", 
                    "Plasma cell",
                    "Monocyte", "Macrophage", 
                    "pDC", "Mast cell",
                    "Endothelial cell",
                    "Fibroblast", "Smooth muscle cell", "AT1", "AT2",
                    "Basal") %>% sort()
celltype_order <- intersect(celltype_order, levels(dat_harmony$celltype))
dat_harmony$celltype <- factor(dat_harmony$celltype, levels = celltype_order)
Idents(dat_harmony) <- "celltype"
set.seed(1)
pal16 <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(celltype_order))
celltypeCols <- setNames(pal16, celltype_order);celltypeCols

p <- DimPlot(dat_harmony, reduction = "tsne", label = F,
             cols = celltypeCols, raster = FALSE, repel = TRUE) +
  theme_bw() +
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        legend.position = "right",
        legend.title = element_text(colour = "black", size = 10),
        legend.text = element_text(colour = "black", size = 9.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        panel.grid = element_blank());p
