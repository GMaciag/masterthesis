# IMPORTS AND SETTINGS -----------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(reshape2)
library(limma)
library(edgeR)
library(ggrepel)
library(AnnotationDbi)
library(ggpubr)
library(GSVA)
library(pheatmap)
library(caret)
library(pROC)

set.seed(2019)

theme_set(theme_pubr())

options(stringsAsFactors = FALSE)

annotation <- read_csv("~/Google Drive/KU/Thesis/R/counts/counting_annotation_entrez.csv")
annotation$GeneID <- mapIds(org.Hs.eg.db::org.Hs.eg.db,
                            keys=as.character(annotation$GeneID), # The IDs
                            column="SYMBOL", 
                            keytype="ENTREZID")
annotation$GeneID[is.na(annotation$GeneID)] <- "NO" 
annotation$GeneID <- make.unique(annotation$GeneID)

# CB friendly palette for plotting.
cbPalette <- c(tumor="#999999", normal="#E69F00", caco="#0072B2", org="#009E73", Ctrl="#F0E442", IBD="#56B4E9", "#D55E00", "#CC79A7")

EM <- read_csv("~/Google Drive/KU/Thesis/R/counts/counts_matrix_entrez.csv")
label <- read_csv("~/Google Drive/KU/Thesis/R/label/goodlabelfinal.csv")

EM_IBD <- read_tsv("/Users/xvk205/Google Drive/KU/Thesis/IBDdata/counts/genesCounts.tsv")
label_IBD <- read_tsv("/Users/xvk205/Google Drive/KU/Thesis/IBDdata/label/design.tsv")

# Remove genderless sample.
noG <- label_IBD$Sample[label_IBD$Gender=="NAN"]
EM_IBD <- EM_IBD[(colnames(EM_IBD)!=noG)]
label_IBD <- label_IBD[label_IBD$Sample!=noG,]

EM_normal <- read_tsv("~/Google Drive/KU/Thesis/R/PCA/data/GSM1697009_06_01_15_TCGA_24_normal_Rsubread_FeatureCounts.txt")
label_normal <- read_tsv("~/Google Drive/KU/Thesis/R/PCA/data/GSE62944_06_01_15_TCGA_24_Normal_CancerType_Samples.txt", col_names = FALSE)
EM_tumor <- read_tsv("~/Google Drive/KU/Thesis/R/PCA/data/GSM1536837_06_01_15_TCGA_24_tumor_Rsubread_FeatureCounts.txt")
label_tumor <- read_tsv("~/Google Drive/KU/Thesis/R/PCA/data/GSE62944_01_27_15_TCGA_20_CancerType_Samples.txt", col_names = FALSE)

# PREPARE LABEL AND EM -----------------------------------------------------------

label <- as.data.frame(label)
label <- label[order(match(label$libid,colnames(EM))),]
rownames(label) <- label$SampleID
label$confluence <- factor(paste0(label$confluence, "%"), levels = c("0%", "50%", "100%"))
label$time[which(label$time==6)] <- 0
label$time <- factor(paste0(label$time, "h"))
label$time <- factor(label$time, levels = c("0h", "4h", "24h"))
label$batch <- factor(paste0("Batch", label$batch))
label$celltype <- factor(label$celltype)

label_caco <- label[label$celltype == "caco",]
label_org <- label[label$celltype == "org",]
label_caco$group <- factor(label_caco$group)
label_org$group <- factor(label_org$group)
label_org$patient <- factor(paste0("Patient", label_org$patient))
caco_status <- ifelse(label$time[label$celltype == "caco"] == "0h", "control", "inflammation/cancer")
org_status <- ifelse(label$time[label$celltype == "org"] == "0h", "control", "inflammation/cancer")

EM_df <- EM[-1] %>% 
  as.data.frame() %>% 
  set_rownames(EM$GeneID) %>% 
  set_colnames(label$SampleID)

genes_EM_df <- mapIds(org.Hs.eg.db::org.Hs.eg.db,
                      keys=as.character(rownames(EM_df)), # The IDs
                      column="SYMBOL", 
                      keytype="ENTREZID") %>% as.vector()
genes_EM_df[is.na(genes_EM_df)] <- "NO"
genes_EM_df <- make.unique(genes_EM_df)
rownames(EM_df) <- genes_EM_df

EM_caco <- EM_df[, colnames(EM_df) %in% label_caco$SampleID]
EM_org <- EM_df[, colnames(EM_df) %in% label_org$SampleID]

# IBD - PREPARE LABEL AND EM -----------------------------------------------------------

label_IBD <- as.data.frame(label_IBD[1:4])
label_IBD <- label_IBD[order(match(label_IBD$Sample,colnames(EM_IBD))),]
rownames(label_IBD) <- label_IBD$Sample
label_IBD$Class <- factor(label_IBD$Class)
label_IBD$Batch <- factor(label_IBD$Batch)
label_IBD$Gender <- factor(label_IBD$Gender)
label_IBD$group <- ifelse(label_IBD$Class == "con", "Ctrl", "IBD")
label_IBD$group <- factor(label_IBD$group)
IBD_status <- ifelse(label_IBD$Class == "con", "control", "inflammation/cancer")

EM_IBD <- EM_IBD[-1] %>% 
  as.data.frame() %>% 
  set_rownames(EM_IBD$GeneID) %>% 
  set_colnames(label_IBD$Sample)

genes_EM_IBD <- mapIds(org.Hs.eg.db::org.Hs.eg.db,
                      keys=as.character(rownames(EM_IBD)), # The IDs
                      column="SYMBOL", 
                      keytype="ENTREZID") %>% as.vector()
genes_EM_IBD[is.na(genes_EM_IBD)] <- "NO"
genes_EM_IBD <- make.unique(genes_EM_IBD)
rownames(EM_IBD) <- genes_EM_IBD

# TCGA - PREPARE LABEL AND EM -----------------------------------------------------------

EM_normal <- EM_normal[-1] %>% 
  as.data.frame() %>% 
  set_rownames(EM_normal$X1)
label_normal <- label_normal[order(match(label_normal$X1,colnames(EM_normal))),]
EM_normal <- EM_normal[, label_normal$X2=="COAD"]

EM_tumor <- EM_tumor[-1] %>% 
  as.data.frame() %>% 
  set_rownames(EM_tumor$X1)
more_names <- colnames(EM_tumor)[!(colnames(EM_tumor)%in%label_tumor$X1)]
label_tumor <- rbind(label_tumor, cbind(X1=more_names, X2="NO"))
label_tumor <- label_tumor[order(match(label_tumor$X1, colnames(EM_tumor))),]
EM_tumor <- EM_tumor[, label_tumor$X2=="COAD"]

EM_TCGA <- cbind(EM_normal, EM_tumor)

label_TCGA <- cbind(sample=colnames(EM_TCGA), group=c(rep("normal", length(EM_normal[1,])), rep("tumor", length(EM_tumor[1,]))))
label_TCGA <- as.data.frame(label_TCGA)
label_TCGA$group <- factor(label_TCGA$group)
TCGA_status <- ifelse(label_TCGA$group == "normal", "control", "inflammation/cancer")

# PREPARE THE DESIGN ------------------------------------------------------

# Add batch to the design.
design_caco <- model.matrix(~ 0 + group, data = label_caco)
batch1_caco <- label_caco$batch
design_org <- model.matrix(~ 0 + group, data = label_org)
batch1_org <- label_org$batch
batch2_org <- label_org$patient
colnames(design_caco) <- unique(levels(label_caco$group))
colnames(design_org) <- unique(levels(label_org$group))

# FILTERING AND NORMALISATION ---------------------------------------------

dge_caco <- DGEList(counts=EM_caco, group = as.numeric(factor(label_caco$group)))
# Remove genes that consistently have zero or very low counts.
keep_caco <- filterByExpr(dge_caco, design_caco)
dge_caco <- dge_caco[keep_caco,,keep.lib.sizes=FALSE]
# Apply scale normalization.
dge_caco <- calcNormFactors(dge_caco)
# Add gene length information. 
dge_caco$genes <- annotation$Length[match(rownames(dge_caco), annotation$GeneID)]
# logRPKM normalisation.
logRPKM_caco <- rpkm(dge_caco, gene.length = dge_caco$genes, log=TRUE, prior.count=1)
# logCPM normalization.
logCPM_caco <- cpm(dge_caco, log=TRUE, prior.count=1)
# Remove variation coming from bacth.  
mat_caco <- removeBatchEffect(logCPM_caco, batch = batch1_caco, design = design_caco)

dge_org <- DGEList(counts=EM_org, group = as.numeric(factor(label_org$group)))
# Remove genes that consistently have zero or very low counts.
keep_org <- filterByExpr(dge_org, design_org)
dge_org <- dge_org[keep_org,,keep.lib.sizes=FALSE]
# Apply scale normalization.
dge_org <- calcNormFactors(dge_org)
# Add gene length information. 
dge_org$genes <- annotation$Length[match(rownames(dge_org), annotation$GeneID)]
# logRPKM normalisation.
logRPKM_org <- rpkm(dge_org, gene.length = dge_org$genes, log=TRUE, prior.count=1)
# logCPM normalization.
logCPM_org <- cpm(dge_org, log=TRUE, prior.count=1)
# Remove variation coming from batch and patient. 
mat_org <- removeBatchEffect(logCPM_org, batch = batch1_org, batch2 = batch2_org, design = design_org)

# IBD - PREPARE THE DESIGN ------------------------------------------------------

# Add batch to the design.
design_ibd <- model.matrix(~ 0 + group, data = label_IBD)
batch1_ibd <- label_IBD$Batch
colnames(design_ibd) <- unique(levels(label_IBD$group))

# IBD - FILTERING AND NORMALISATION ---------------------------------------------

dge_IBD <- DGEList(counts=EM_IBD, group = as.numeric(factor(label_IBD$group)))
# Remove genes that consistently have zero or very low counts.
keep_IBD <- filterByExpr(dge_IBD, design_ibd)
dge_IBD <- dge_IBD[keep_IBD,,keep.lib.sizes=FALSE]
# Apply scale normalization.
dge_IBD <- calcNormFactors(dge_IBD)
# Add gene length information. 
dge_IBD$genes <- annotation$Length[match(rownames(dge_IBD), annotation$GeneID)]
# logRPKM normalisation.
logRPKM_IBD <- rpkm(dge_IBD, gene.length = dge_IBD$genes, log=TRUE, prior.count=1)
# logCPM normalization.
logCPM_IBD <- cpm(dge_IBD, log=TRUE, prior.count=1)
# Remove variation coming from bacth.  
mat_IBD <- removeBatchEffect(logCPM_IBD, batch = batch1_ibd, design = design_ibd)

# TCGA - PREPARE THE DESIGN ------------------------------------------------------

# Add batch to the design.
design_TCGA <- model.matrix(~ 0 + group, data = label_TCGA)
colnames(design_TCGA) <- unique(levels(label_TCGA$group))

# TCGA - FILTERING AND NORMALISATION ---------------------------------------------

dge_TCGA <- DGEList(counts=EM_TCGA, group = as.numeric(factor(label_TCGA$group)))
# Remove genes that consistently have zero or very low counts.
keep_TCGA <- filterByExpr(dge_TCGA, design_TCGA)
dge_TCGA <- dge_TCGA[keep_TCGA,,keep.lib.sizes=FALSE]
# Apply scale normalization.
dge_TCGA <- calcNormFactors(dge_TCGA)
# Add gene length information. 
dge_TCGA$genes <- annotation$Length[match(rownames(dge_TCGA), annotation$GeneID)]
# logRPKM normalisation.
logRPKM_TCGA <- rpkm(dge_TCGA, gene.length = dge_TCGA$genes, log=TRUE, prior.count=1)
# logCPM normalization.
logCPM_TCGA <- cpm(dge_TCGA, log=TRUE, prior.count=1)
# Remove variation coming from bacth.  
mat_TCGA <- logCPM_TCGA

# COMMON GENES ------------------------------------------------------------

keep.genes <- intersect(rownames(mat_caco), rownames(mat_org)) %>% intersect(rownames(mat_IBD)) %>% intersect(rownames(mat_TCGA))

mat_caco <- mat_caco[rownames(mat_caco) %in% keep.genes,]
mat_org <- mat_org[rownames(mat_org) %in% keep.genes,]
mat_IBD <- mat_IBD[rownames(mat_IBD) %in% keep.genes,]
mat_TCGA <- mat_TCGA[rownames(mat_TCGA) %in% keep.genes,]

mat_org <- mat_org[order(match(rownames(mat_org), rownames(mat_caco))),]
mat_IBD <- mat_IBD[order(match(rownames(mat_IBD), rownames(mat_caco))),]
mat_TCGA <- mat_TCGA[order(match(rownames(mat_TCGA), rownames(mat_caco))),]

# PCA ---------------------------------------------------------------------

pca_TCGA <- prcomp(t(mat_TCGA), scale. = TRUE)

TCGA_in_TCGA <- pca_TCGA$x %>% as.tibble()
caco_in_TCGA <- predict(pca_TCGA, t(mat_caco)) %>% as.tibble()
org_in_TCGA <- predict(pca_TCGA, t(mat_org)) %>% as.tibble()
IBD_in_TCGA <- predict(pca_TCGA, t(mat_IBD)) %>% as.tibble()

plotPCA <- function(i, j) {
  comp1 <- paste0("PC", i)
  comp2 <- paste0("PC", j)
  pc1var <- round(summary(pca_TCGA)$importance[2,i]*100, digits=1)
  pc2var <- round(summary(pca_TCGA)$importance[2,j]*100, digits=1)
  ggplot() +
    geom_point(data = TCGA_in_TCGA, aes_string(x=comp1, y=comp2, color=label_TCGA$group), alpha=0.6, size = 4) +
    geom_point(data = IBD_in_TCGA, aes_string(x=comp1, y=comp2, color=label_IBD$group), alpha=0.6, size = 4) +
    geom_point(data = caco_in_TCGA, aes_string(x=comp1, y=comp2, color=label_caco$celltype), alpha=0.6, size = 4) +
    geom_point(data = org_in_TCGA, aes_string(x=comp1, y=comp2, color=label_org$celltype), alpha=0.6, size = 4) +
    labs(x = sprintf('PC%s (%s%%)', i, pc1var), y = sprintf('PC%s (%s%%)', j, pc2var)) +
    guides(col=guide_legend(title="dataset")) +
    scale_color_manual(values = cbPalette, breaks = c("normal", "tumor", "caco", "org", "Ctrl", "IBD")) +  
    theme(text = element_text(size=20))
}
p <- vector("list", 4)
n <- 1
for (i in c(1,3,5,7)) {
    p[[n]] <- plotPCA(i = i, j = i+1)
    n <- n+1
}

figure <- ggarrange(plotlist = p,
                    ncol = 2, nrow = 2,
                    common.legend = TRUE, legend = "none")
figure

# KNN --------------------------------------------------------

cumsum(summary(pca_TCGA)$importance[2,1:60]*100)

# 51.0% - PC 1:6
# 61.8% - PC 1:14
# 70.3% - PC 1:26
# 80.1% - PC 1:56

## Prepare sets.

ALL_set <- TCGA_in_TCGA

control_validation <- sample(41, 10)
tumor_validation <- sample(seq(42, 509), 117)

VALIDATE_set <- ALL_set[c(control_validation, tumor_validation),]
VALIDATE_label <- label_TCGA$group[c(control_validation, tumor_validation)]
VALIDATE_label_predict <- VALIDATE_label
levels(VALIDATE_label_predict) <- c(0, 1)
TRAIN_set <- ALL_set[-c(control_validation, tumor_validation),]
TRAIN_label <- label_TCGA$group[-c(control_validation, tumor_validation)]

TEST_set <- rbind(TCGA_in_TCGA, caco_in_TCGA, org_in_TCGA, IBD_in_TCGA)


## Train.

trControl <- trainControl(method  = "cv",
                          number  = 5,
                          summaryFunction = twoClassSummary,
                          classProbs = TRUE,
                          savePredictions='all')
knn_plots <- list()
for (i in c(1, 6, 14, 26, 56)) {
  
  fit <- train(TRAIN_set[,1:i], TRAIN_label,
               method     = "knn",
               tuneGrid   = expand.grid(k = 1:20),
               trControl  = trControl)
  print(sprintf("####################    No. of PCs == %s    ####################", i))
  print(fit)
  knn_plots[[i]] <- plot(fit, cex=1.7)
  knn_plots[[i+1]] <- plot.roc(fit$pred$obs,
           fit$pred$tumor)
  
}

  # 60% of variance (with 14 PCs) gives Sensitivity of 1. Specificity is all the same for all k. 
# I will choose moderate k, k = 10.

## Validate.

knn_fit <- train(x=TRAIN_set[,1:14], TRAIN_label,
                 method     = "knn",
                 tuneGrid   = expand.grid(k = 17),
                 trControl  = trControl,
                 metric     = "ROC")

knn_predict <- predict(knn_fit, VALIDATE_set[,1:14])

# summarize results
confusionMatrix(knn_predict, VALIDATE_label)
plot.roc(knn_predict,
         as.numeric(VALIDATE_label_predict), cex.lab=1.7, cex.axis=1.7, cex.main=1.7, cex.sub=1.7)


## Test.

knn_predict_test <- predict(knn_fit, TEST_set[,1:14])

results_caco_knn <- cbind(prediction=knn_predict_test[510:526], label_caco)
results_org_knn <- cbind(prediction=knn_predict_test[527:566], label_org)
results_IBD_knn <- cbind(prediction=knn_predict_test[567:639], label_IBD)

results_caco_knn <- table(results_caco_knn$prediction, status=caco_status) %>% melt() %>% cbind(dataset="caco")
results_org_knn <- table(results_org_knn$prediction, status=org_status) %>% melt() %>% cbind(dataset="org")
results_IBD_knn <- table(results_IBD_knn$prediction, status=IBD_status) %>% melt() %>% cbind(dataset="IBD")

results_knn <- rbind(results_caco_knn, results_org_knn) %>% rbind(results_IBD_knn)
results_knn <-results_knn %>% cbind(fraction=c(results_knn$value[1:2]/6, results_knn$value[3:4]/11, results_knn$value[5:6]/16, results_knn$value[7:8]/24, results_knn$value[9:10]/29, results_knn$value[11:12]/44))

ggplot(data=results_knn, aes(x = status, y = fraction, fill = Var1)) + facet_wrap(~factor(dataset, levels = c("caco", "org", "IBD"))) +
  geom_bar(stat = "identity", position = position_stack()) +
  xlab("dataset") +
  ylab("fraction") +
  labs(fill = "prediction") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=15),
        plot.title = element_text(hjust = 0.5, size = 10),
        text = element_text(size=20)) +
  scale_fill_manual(labels = c("normal", "tumor"), values = c("#D3D3D3", "#696969")) +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  scale_x_discrete(labels = c("control", "inflammation"))

# LDA --------------------------------------------------------

## Train.

trControl <- trainControl(method  = "cv",
                          number  = 5,
                          summaryFunction = twoClassSummary,
                          classProbs = TRUE,
                          savePredictions='all')
lda_plots <- list()
for (i in c(6, 14, 26, 56)) {
  
  fit <- train(TRAIN_set[,1:i], TRAIN_label,
               method     = "lda",
               trControl  = trControl)
  print(sprintf("####################    No. of PCs == %s    ####################", i))
  print(fit)
  lda_plots[[i+1]] <- plot.roc(fit$pred$obs,
                               fit$pred$tumor)
  
}
# 60% of variance (with 14 PCs) gives Sensitivity of 1. Specificity is 1 for 56 PCs, but very high (0.997) for lower #PCs.
# I'll choose to use 14 PCs.

## Validate.

lda_fit <- train(x=TRAIN_set[,1:14], TRAIN_label,
                 method     = "lda",
                 trControl  = trControl,
                 metric     = "ROC")

lda_predict <- predict(lda_fit, VALIDATE_set[,1:14])

# summarize results
confusionMatrix(lda_predict, VALIDATE_label)
plot.roc(lda_predict,
         as.numeric(VALIDATE_label_predict), cex.lab=1.7, cex.axis=1.7, cex.main=1.7, cex.sub=1.7)


## Test.

lda_predict_test <- predict(lda_fit, TEST_set[,1:14])

results_caco_lda <- cbind(prediction=lda_predict_test[510:526], label_caco)
results_org_lda <- cbind(prediction=lda_predict_test[527:566], label_org)
results_IBD_lda <- cbind(prediction=lda_predict_test[567:639], label_IBD)

results_caco_lda <- table(results_caco_lda$prediction, status=caco_status) %>% melt() %>% cbind(dataset="caco")
results_org_lda <- table(results_org_lda$prediction, status=org_status) %>% melt() %>% cbind(dataset="org")
results_IBD_lda <- table(results_IBD_lda$prediction, status=IBD_status) %>% melt() %>% cbind(dataset="IBD")

results_lda <- rbind(results_caco_lda, results_org_lda) %>% rbind(results_IBD_lda)
results_lda <- results_lda %>% cbind(fraction=c(results_lda$value[1:2]/6, results_lda$value[3:4]/11, results_lda$value[5:6]/16, results_lda$value[7:8]/24, results_lda$value[9:10]/29, results_lda$value[11:12]/44))

ggplot(data=results_lda, aes(x = status, y = fraction, fill = Var1)) + facet_wrap(~factor(dataset, levels = c("caco", "org", "IBD"))) +
  geom_bar(stat = "identity", position = position_stack()) +
  xlab("dataset") +
  ylab("fraction") +
  labs(fill = "prediction") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size=15),
        plot.title = element_text(hjust = 0.5, size = 10),
        text = element_text(size=20)) +
  scale_fill_manual(labels = c("normal", "tumor"), values = c("#D3D3D3", "#696969")) +
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  scale_x_discrete(labels = c("control", "inflammation"))
# GSVA --------------------------------------------------------------------

keep.genes.gsva <- intersect(rownames(logRPKM_caco), rownames(logRPKM_org)) %>% intersect(rownames(logRPKM_IBD)) %>% intersect(rownames(logRPKM_TCGA))

logRPKM_caco <- logRPKM_caco[rownames(logRPKM_caco) %in% keep.genes.gsva,]
logRPKM_org <- logRPKM_org[rownames(logRPKM_org) %in% keep.genes.gsva,]
logRPKM_IBD <- logRPKM_IBD[rownames(logRPKM_IBD) %in% keep.genes.gsva,]
logRPKM_TCGA <- logRPKM_TCGA[rownames(logRPKM_TCGA) %in% keep.genes.gsva,]

logRPKM_org <- logRPKM_org[order(match(rownames(logRPKM_org), rownames(logRPKM_caco))),]
logRPKM_IBD <- logRPKM_IBD[order(match(rownames(logRPKM_IBD), rownames(logRPKM_caco))),]
logRPKM_TCGA <- logRPKM_TCGA[order(match(rownames(logRPKM_TCGA), rownames(logRPKM_caco))),]

logRPKM_data <- cbind(logRPKM_caco, logRPKM_org, logRPKM_IBD, logRPKM_TCGA)

HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION <- read_tsv(
  "genesets/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.txt", col_names = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION") %>% as.list()
HALLMARK_MYC_TARGETS_V2 <- read_tsv(
  "genesets/HALLMARK_MYC_TARGETS_V2.txt", col_names = "HALLMARK_MYC_TARGETS_V2") %>% as.list()
GRADE_COLON_CANCER_UP <- read_tsv(
  "genesets/GRADE_COLON_CANCER_UP.txt", col_names = "GRADE_COLON_CANCER_UP") %>% as.list()
CAMPS_COLON_CANCER_COPY_NUMBER_UP <- read_tsv(
  "genesets/CAMPS_COLON_CANCER_COPY_NUMBER_UP.txt", col_names = "CAMPS_COLON_CANCER_COPY_NUMBER_UP") %>% as.list()
SABATES_COLORECTAL_ADENOMA_DN <- read_tsv(
  "genesets/SABATES_COLORECTAL_ADENOMA_DN.txt", col_names = "SABATES_COLORECTAL_ADENOMA_DN") %>% as.list()
IBD_GENES <- read_tsv("genesets/IBD_UP.txt", 
                      col_names = c("IBD_UP", "desc"), 
                      col_types = cols(col_character(), col_character()))[,-2] %>% as.list()
IBD_GENES[[1]] <- mapIds(org.Hs.eg.db::org.Hs.eg.db,
                         keys=IBD_GENES[[1]], # The IDs
                         column="SYMBOL", 
                         keytype="ENTREZID")
TNF_GENES <- read_tsv("genesets/TNF_UP.txt", 
                      col_names = c("TNF_UP", "desc"), 
                      col_types = cols(col_character(), col_character()))[,-2] %>% as.list()
TNF_GENES[[1]] <- mapIds(org.Hs.eg.db::org.Hs.eg.db,
                         keys=TNF_GENES[[1]], # The IDs
                         column="SYMBOL", 
                         keytype="ENTREZID")
INFLAMMATORY_RESPONSE_UP_GENES <- read_tsv("genesets/INFLAMMATORY_RESPONSE_UP.txt", 
                                           col_names = c("INFLAMMATORY_RESPONSE_UP"), 
                                           col_types = cols(col_character())) %>% as.list()
IMMUNE_RESPONSE_UP_GENES <- read_tsv("genesets/IMMUNE_RESPONSE_UP.txt", 
                                     col_names = c("IMMUNE_RESPONSE_UP"), 
                                     col_types = cols(col_character())) %>% as.list()

plotGSVA <- function(gset) {
  
  all_gsva <- gsva(logRPKM_data, gset.idx.list = gset, method = "ssgsea") %>% 
    t() %>% 
    as.numeric() %>% 
    cbind(X1=., group=c(rep("caco", length(label_caco[,1])), rep("org", length(label_org[,1])), rep("IBD", length(label_IBD[,1])), rep("TCGA", length(label_TCGA[,1]))),
          status=c(caco_status, org_status, IBD_status, TCGA_status)) %>% 
    as.data.frame()
  
  ggplot(all_gsva, aes(
    x = factor(group, levels = c("caco", "org", "IBD", "TCGA")),
    y = as.numeric(X1),
    min = 0,
    fill=factor(status, levels = c("control", "inflammation/cancer")))) +
    geom_boxplot(alpha = .75) +
    ylab("Enrichment score") +
    xlab("dataset") +
    ggtitle(names(gset)) +
    theme(plot.title = element_text(hjust = 0.5, size = 10),
          text = element_text(size=20),
          legend.title = element_blank()) +
    scale_fill_manual(labels = c("control", "inflammation/cancer"), values = c("#D3D3D3", "#696969")) %>% return()
}

p1 <- plotGSVA(SABATES_COLORECTAL_ADENOMA_DN)
p2 <- plotGSVA(HALLMARK_MYC_TARGETS_V2)
p3 <- plotGSVA(GRADE_COLON_CANCER_UP)
p4 <- plotGSVA(CAMPS_COLON_CANCER_COPY_NUMBER_UP)

ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, common.legend = TRUE, legend = "top") %>% ggsave(filename = "plots/gsva/cancer_withTGCA_ssgseaMethod.pdf")

p5 <- plotGSVA(IBD_GENES)
p6 <- plotGSVA(TNF_GENES)
p7 <- plotGSVA(INFLAMMATORY_RESPONSE_UP_GENES)
p8 <- plotGSVA(IMMUNE_RESPONSE_UP_GENES)

ggarrange(p5, p6, p7, p8, ncol = 2, nrow = 2, common.legend = TRUE, legend = "top") %>% ggsave(filename = "plots/gsva/inflammation_withTGCA_ssgseaMethod.pdf")

# MULTI GSVA --------------------------------------------------------------

load("~/Google Drive/KU/Thesis/R/genesets/whole_sets/human_c2_v5p2.rdata")
load("~/Google Drive/KU/Thesis/R/genesets/whole_sets/human_c6_v5p2.rdata")
load("~/Google Drive/KU/Thesis/R/genesets/whole_sets/human_c7_v5p2.rdata")

IDmapper <- function(x) {
  y <- mapIds(org.Hs.eg.db::org.Hs.eg.db,
              keys=x, # The IDs
              column="SYMBOL", 
              keytype="ENTREZID")
  return(as.character(y))
}

c2_CANCER <- Hs.c2[grep("CANCER", names(Hs.c2))]
c2_CANCER <- lapply(c2_CANCER, IDmapper)
c2_CANCER_UP <- c2_CANCER[grep("UP", names(c2_CANCER))]
c2_CANCER_DOWN <- c2_CANCER[grep("DN", names(c2_CANCER))]

c2_INFLA <- Hs.c2[grep("INFLA", names(Hs.c2))]
c2_INFLA <- lapply(c2_INFLA, IDmapper)
c2_INFLA_UP <- c2_INFLA[grep("UP", names(c2_INFLA))]
c2_INFLA_DOWN <- c2_INFLA[grep("DN", names(c2_INFLA))]

c2_INFECT <- Hs.c2[grep("INFECT", names(Hs.c2))]
c2_INFECT <- lapply(c2_INFECT, IDmapper)
c2_INFECT_UP <- c2_INFECT[grep("UP", names(c2_INFECT))]
c2_INFECT_DOWN <- c2_INFECT[grep("DN", names(c2_INFECT))]

c7_TNF <- Hs.c7[grep("TNF", names(Hs.c7))]
c7_TNF <- lapply(c7_TNF, IDmapper)
c7_TNF_UP <- c7_TNF[grep("UP", names(c7_TNF))]
c7_TNF_DOWN <- c7_TNF[grep("DN", names(c7_TNF))]

plotMultiGSVA <- function(gset, name, method) {
  
  all_gsva <- gsva(logRPKM_data, gset.idx.list = gset, method = method)
  all_gsva <- all_gsva %>% t() %>% as.data.frame()
  all_gsva <- data.frame(apply(all_gsva, 2, function(x) as.numeric(as.character(x)))) %>% 
    cbind(group=c(rep("caco", length(label_caco[,1])), rep("org", length(label_org[,1])), rep("IBD", length(label_IBD[,1])), rep("TCGA", length(label_TCGA[,1]))),
          status=c(caco_status, org_status, IBD_status, TCGA_status), .)
  all_gsva$status[which(all_gsva$status=="control")] <- "con"
  all_gsva$status[which(all_gsva$status=="inflammation/cancer")] <- "infl"
  
  all_gsva <- all_gsva %>% group_by(group, status) %>% summarise_all(median)
  
  all_gsva <- data.frame(
    rbind(
      cbind(score=all_gsva[1,-c(1,2)], group=all_gsva[1,1], status=all_gsva[1,2], geneset=colnames(all_gsva[-c(1,2)])),
      cbind(score=all_gsva[2,-c(1,2)], group=all_gsva[2,1], status=all_gsva[2,2], geneset=colnames(all_gsva[-c(1,2)])),
      cbind(score=all_gsva[3,-c(1,2)], group=all_gsva[3,1], status=all_gsva[3,2], geneset=colnames(all_gsva[-c(1,2)])),
      cbind(score=all_gsva[4,-c(1,2)], group=all_gsva[4,1], status=all_gsva[4,2], geneset=colnames(all_gsva[-c(1,2)])),
      cbind(score=all_gsva[5,-c(1,2)], group=all_gsva[5,1], status=all_gsva[5,2], geneset=colnames(all_gsva[-c(1,2)])),
      cbind(score=all_gsva[6,-c(1,2)], group=all_gsva[6,1], status=all_gsva[6,2], geneset=colnames(all_gsva[-c(1,2)])),
      cbind(score=all_gsva[7,-c(1,2)], group=all_gsva[7,1], status=all_gsva[7,2], geneset=colnames(all_gsva[-c(1,2)])),
      cbind(score=all_gsva[8,-c(1,2)], group=all_gsva[8,1], status=all_gsva[8,2], geneset=colnames(all_gsva[-c(1,2)]))))
  
  all_gsva$score <- as.numeric(all_gsva$score)
  all_gsva$group <- as.character(all_gsva$group)
  all_gsva$status <- as.character(all_gsva$status)
  all_gsva$geneset <- as.character(all_gsva$geneset)
  
  all_gsva <- all_gsva %>% 
    tidyr::spread(status, score) %>%
    dplyr::mutate(is_increasing = infl > con) %>%
    tidyr::gather("status", "score", 3:4)
  
  ggplot(all_gsva, aes(x = status, y = score)) + facet_wrap(~group) +
    geom_boxplot(aes(fill = status), alpha = 0.75, col = "black") +
    geom_point(col = "black") +
    geom_line(aes(group = geneset, col = is_increasing), alpha = 0.5, size = 2) +
    scale_colour_manual(values = c("blue", "orange"), labels = c("Decrease", "Increase")) +
    ylab("Enrichment score") +
    xlab("dataset") +
    ggtitle(sprintf("%s genesets (n = %s)", name, length(gset))) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    theme(plot.title = element_text(hjust = 0.5, size = 25),
          text = element_text(size=20),
          legend.title = element_blank(),
          legend.position = "top",
          legend.justification = "center",
          strip.background = element_blank(),
          axis.text.x = element_blank()) +
    scale_fill_manual(labels = c("control", "inflammation/cancer"), values = c("#D3D3D3", "#696969")) %>% return()
}

plotMultiGSVA(c2_CANCER_UP, "CANCER UP", "ssgsea") %>% ggsave(filename = "plots/gsva/c2_CANCER_UP_ssgsea.pdf")

plotMultiGSVA(c2_CANCER_DOWN, "CANCER DOWN", "ssgsea") %>% ggsave(filename = "plots/gsva/c2_CANCER_DOWN_ssgsea.pdf")

plotMultiGSVA(c2_INFLA_UP, "INFLAMMATION UP", "ssgsea") %>% ggsave(filename = "plots/gsva/c2_INFLAMMATION_UP_ssgsea.pdf")

plotMultiGSVA(c2_INFLA_DOWN, "INFLAMMATION DOWN", "ssgsea") %>% ggsave(filename = "plots/gsva/c2_INFLAMMATION_DOWN_ssgsea.pdf")

plotMultiGSVA(c2_INFECT_UP, "INFECTION UP", "ssgsea") %>% ggsave(filename = "plots/gsva/c2_INFECTION_UP_ssgsea.pdf")

plotMultiGSVA(c2_INFECT_DOWN, "INFECTION DOWN", "ssgsea") %>% ggsave(filename = "plots/gsva/c2_INFECTION_DOWN_ssgsea.pdf")

plotMultiGSVA(c7_TNF_UP, "TNF UP", "ssgsea") %>% ggsave(filename = "plots/gsva/c7_TNF_UP_ssgsea.pdf")

plotMultiGSVA(c7_TNF_DOWN, "TNF DOWN", "ssgsea") %>% ggsave(filename = "plots/gsva/c7_TNF_DOWN_ssgsea.pdf")
