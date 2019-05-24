# IMPORTS AND SETTINGS -----------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(reshape2)
library(limma)
library(edgeR)
library(pathview)
library(ggrepel)
library(pheatmap)
library(rowr)
library(Vennerable)
library(ggpubr)
library(GSVA)
library(parallel)
library(viridis)

options(stringsAsFactors = FALSE)

# CB friendly palette for plotting.
cbPalette <- c("#56B4E9", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette_resp <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#000000", "#0072B2", "#D55E00", "#CC79A7")

EM <- read_csv("~/Google Drive/KU/Thesis/R/counts/counts_matrix_entrez.csv")
label <- read_csv("~/Google Drive/KU/Thesis/R/label/goodlabelfinal.csv")
annotation <- read_csv("~/Google Drive/KU/Thesis/R/counts/counting_annotation_entrez.csv")

# PREPARE DESIGN -----------------------------------------------------------

label <- as.data.frame(label)
label <- label[order(match(label$libid,colnames(EM))),]
rownames(label) <- label$SampleID
label$confluence <- factor(paste0(label$confluence, "%"), levels = c("0%", "50%", "100%"))
label$time[which(label$time==6)] <- 0
label$time <- factor(paste0(label$time, "h"))
label$time <- factor(label$time, levels = c("0h", "4h", "24h"))
label$batch <- factor(paste0("Batch", label$batch))
label$patient <- factor(paste0("Patient", label$patient))
label$celltype <- factor(label$celltype)
label$group <- factor(label$group)

# Add batch to the design.
design <- model.matrix(~ 0 + group + batch, data = label)
colnames(design) <- c(unique(levels(label$group)), unique(levels(label$batch)[-1]))
# Make comparisons between the experimental conditions.
contrasts <- makeContrasts(
  "caco50_at_0to4"=caco.50.4-caco.50.0,
  "caco50_at_4to24"=caco.50.24-caco.50.4,
  "caco50_at_0to24"=caco.50.24-caco.50.0,
  "caco100_at_0to4"=caco.100.4-caco.100.0,
  "caco100_at_4to24"=caco.100.24-caco.100.4,
  "caco100_at_0to24"=caco.100.24-caco.100.0,
  "org_at_0to4"=org.4-org.0,
  "org_at_4to24"=org.24-org.4,
  "org_at_0to24"=org.24-org.0,
  "caco50_vs_caco100_at_0"=caco.100.0-caco.50.0,
  "caco50_vs_caco100_at_4"=caco.100.4-caco.50.4,
  "caco50_vs_caco100_at_24"=caco.100.24-caco.50.24,
  "caco50_vs_org_at_0"=org.0-caco.50.0,
  "caco50_vs_org_at_4"=org.4-caco.50.4,
  "caco50_vs_org_at_24"=org.24-caco.50.24,
  "caco100_vs_org_at_0"=org.0-caco.100.0,
  "caco100_vs_org_at_4"=org.4-caco.100.4,
  "caco100_vs_org_at_24"=org.24-caco.100.24,
  
  "caco50_vs_caco100_at_0to4"=(caco.50.4-caco.50.0)-(caco.100.4-caco.100.0),
  "caco50_vs_caco100_at_4to24"=(caco.50.24-caco.50.4)-(caco.100.24-caco.100.4),
  "caco50_vs_caco100_at_0to24"=(caco.50.24-caco.50.0)-(caco.100.24-caco.100.0),
  
  "caco50_vs_org_at_0to4"=(caco.50.4-caco.50.0)-(org.4-org.0),
  "caco50_vs_org_at_4to24"=(caco.50.24-caco.50.4)-(org.24-org.4),
  "caco50_vs_org_at_0to24"=(caco.50.24-caco.50.0)-(org.24-org.0),
  
  "caco100_vs_org_at_0to4"=(caco.100.4-caco.100.0)-(org.4-org.0),
  "caco100_vs_org_at_4to24"=(caco.100.24-caco.100.4)-(org.24-org.4),
  "caco100_vs_org_at_0to24"=(caco.100.24-caco.100.0)-(org.24-org.0),
  levels=design
)

# FILTERING AND NORMALIZATION ---------------------------------------------

EM_df <- EM[-1] %>% 
  as.data.frame() %>% 
  set_rownames(EM$GeneID) %>% 
  set_colnames(label$SampleID)
dge <- DGEList(counts=EM_df, group = as.numeric(factor(label$group)))
# Remove genes that consistently have zero or very low counts.
keep <- filterByExpr(dge, design)
dge <- dge[keep,,keep.lib.sizes=FALSE]
# Apply scale normalization.
dge <- calcNormFactors(dge)
# logCPM normalization.
logCPM <- cpm(dge, log=TRUE, prior.count=1)
# Estimate the correlation between measurements made on the same subject.
corfit <- duplicateCorrelation(logCPM, design, block=label$patient)
# Fit a linear model with inter-subject correlation.
fit <- lmFit(logCPM, design, block=label$patient, correlation=corfit$consensus)
# Compute the contrasts.
fit2 <- contrasts.fit(fit, contrasts)
# Apply empirical Bayes smoothing to the standard errors.
eb <- eBayes(fit2, trend=TRUE, robust = TRUE)
# Add gene information. 
eb$genes <- as_data_frame(annotation)[match(rownames(eb), annotation$GeneID),]
# Top genes differential expressed.
results <- decideTests(eb)
# DE GENES ----------------------------------------------------------------

topgenes <- list()
for (i in seq(length(contrasts[1,]))) {
  name <- names(contrasts[1,])[i]
  t <- topTable(eb, coef = i, number = 10)
  topgenes[[name]] <- t
  t$GeneID <- mapIds(org.Hs.eg.db::org.Hs.eg.db,
                     keys=as.character(t$GeneID), # The IDs
                     column="SYMBOL", 
                     keytype="ENTREZID")
  colnames(t)[1] <- "GeneSymbol"
  write.csv(t, sprintf('DEgenes/%s.csv', name), row.names = TRUE)
}

# For STRING. 
topgenes <- list()
for (i in seq(length(contrasts[1,]))) {
  name <- names(contrasts[1,])[i]
  t <- topTable(eb, coef = i, number = 2000, p.value=0.05)
  topgenes[[name]] <- t
  t$GeneID <- mapIds(org.Hs.eg.db::org.Hs.eg.db,
                     keys=as.character(t$GeneID), # The IDs
                     column="SYMBOL", 
                     keytype="ENTREZID")
  colnames(t)[1] <- ""
  t_string <- t[1]
  write.csv(t_string, sprintf('DEgenes/STRING/%s.csv', name), row.names = FALSE, col.names = FALSE)
}


# DE COUNTS -----------------------------------------------------------------

# Plots of number of differentially expressed genes.
res.plot <- results %>% 
  summary() %>% 
  as_tibble(.name_repair = c("minimal")) %>% 
  set_colnames(c("expression", "group", "count")) %>% 
  cbind(comparison=c(rep("intra", 27), rep("inter", 27), rep("complx", 27)))
res.plot$group <- factor(res.plot$group, levels = unique(res.plot$group))
res.plot$comparison <- factor(res.plot$comparison, levels = c("intra", "inter", "complx"))
res.plot.intra <- filter(res.plot, comparison == "intra")
res.plot.inter <- filter(res.plot, comparison == "inter")
res.plot.complx <- filter(res.plot, comparison == "complx")

plotDEGnumbers <- function(res.plot) {
  ggplot(res.plot, aes(group)) +
    geom_bar(data =  subset(res.plot, expression == "Up"),
             aes(y = count, fill = expression), stat = "identity", position = "dodge") +
    geom_bar(data = subset(res.plot, expression == "Down"), 
             aes(y = -count, fill = expression), stat = "identity", position = "dodge") + 
    geom_hline(yintercept = 0,colour = "grey90") +
    theme_pubr() + theme(axis.text.x = element_text(angle = 60, hjust = 1, size=15),
                       strip.background = element_blank(),
                       strip.text.x = element_blank(), 
                       text = element_text(size=20),
                       legend.title = element_blank()) +
    ylab("No. of genes") + xlab("Interaction") + 
    guides(fill=guide_legend(title="DE")) + 
    scale_fill_manual(values=cbPalette[1:2], breaks = c("Up", "Down"), labels = c("DE Up", "DE Down")) + 
    scale_y_continuous(labels=abs) +
    geom_text(data = subset(res.plot, expression == "Up"), 
              aes(x = group, y = count, group = expression, label = count),
              position = position_dodge(width=0.9), vjust = -.2, size=6) +
    geom_text(data = subset(res.plot, expression == "Down"), 
              aes(x = group, y = -count, group = expression, label = count),
              position = position_dodge(width=0.9), vjust = 1.2, size=6)
}

plotDEGnumbers(res.plot.intra)
plotDEGnumbers(res.plot.inter)
plotDEGnumbers(res.plot.complx)

# DE COMPARISON -----------------------------------------------------------

plotInteraction <- function(group1, group2, interaction) {
  r <- data.frame(
    eb$coefficients[,group1, drop = FALSE], 
    eb$coefficients[,group2, drop = FALSE], 
    Interaction=eb$coefficients[,interaction])
  r$gene.names <- mapIds(org.Hs.eg.db::org.Hs.eg.db,
                        keys=rownames(r), # The IDs
                        column="SYMBOL", 
                        keytype="ENTREZID")
  p <- ggplot(r, aes(x = r[,1], y = r[,2])) + 
    geom_point(aes(colour = Interaction), size = 3) + 
    geom_abline() +
    coord_fixed(xlim = c(-8, 8), ylim = c(-8, 8)) +
    theme_pubr() +
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size=20)) +
    scale_color_gradient2() +
    xlab(group1) +
    ylab(group2)
  ggsave(filename = sprintf("plots/interactions/%s.pdf", interaction), plot = p, device = "pdf")
}

plotInteraction(group2="caco50_at_0to4", group1="caco100_at_0to4", interaction="caco50_vs_caco100_at_0to4")
plotInteraction(group2="caco50_at_4to24", group1="caco100_at_4to24", interaction="caco50_vs_caco100_at_4to24")
plotInteraction(group2="caco50_at_0to24", group1="caco100_at_0to24", interaction="caco50_vs_caco100_at_0to24")

plotInteraction(group2="caco50_at_0to4", group1="org_at_0to4", interaction="caco50_vs_org_at_0to4")
plotInteraction(group2="caco50_at_4to24", group1="org_at_4to24", interaction="caco50_vs_org_at_4to24")
plotInteraction(group2="caco50_at_0to24", group1="org_at_0to24", interaction="caco50_vs_org_at_0to24")

plotInteraction(group2="caco100_at_0to4", group1="org_at_0to4", interaction="caco100_vs_org_at_0to4")
plotInteraction(group2="caco100_at_4to24", group1="org_at_4to24", interaction="caco100_vs_org_at_4to24")
plotInteraction(group2="caco100_at_0to24", group1="org_at_0to24", interaction="caco100_vs_org_at_0to24")

# GO & KEGG ---------------------------------------------------------------

# Gene ontology analysis.
for (i in seq(1:length(eb$coefficients[1,]))) {
  name <- colnames(eb$coefficients)[i]
  g <- goana(eb, coef=i, trend = TRUE)
  g <- g %>% filter(N<2000) %>% topGO(n=20,truncate="50") 
  write.csv(g, sprintf('GOterms/%s.csv', name))
}
# KEGG pathways analysis.
for (i in seq(1:length(eb$coefficients[1,]))) {
  name <- colnames(eb$coefficients)[i]
  k <- kegga(eb, coef=i, trend = TRUE)
  k <- k %>% filter(N<2000) %>% topKEGG(n=20,truncate="50")
  write.csv(k, sprintf('KEGGterms/%s.csv', name))
}
# KEGG pathways with IBD universe.
gene.pathway <- rbind(cbind(IBD_UP, "1"), cbind(IBD_DOWN, "2"))
pathway.names <- cbind(c("1", "2"), c("IBD UP", "IBD DOWN"))

for (i in seq(1:length(eb$coefficients[i,]))) {
  name <- colnames(eb$coefficients)[i]
  k <- kegga(eb, coef=i, trend = TRUE, 
             gene.pathway = gene.pathway, 
             pathway.names = pathway.names)
  #k <- k %>% topKEGG()
  write.csv(k, sprintf('KEGGterms/IBD/%s.csv', name))
}

# Do the pathview analysis on top120 genes

top120_caco50_0to4_UP <- topTable(eb, coef = "caco50_at_0to4", number = Inf)
top120_caco50_0to4_UP <- top120_caco50_0to4_UP[which(top120_caco50_0to4_UP$logFC > 0), ] [1:120,]
results_caco50_0to4_UP <- results@.Data[,"caco50_at_0to4"]
results_caco50_0to4_UP <- results_caco50_0to4_UP[names(results_caco50_0to4_UP) %in% top120_caco50_0to4_UP$GeneID]

top120_caco100_0to4_UP <- topTable(eb, coef = "caco100_at_0to4", number = Inf)
top120_caco100_0to4_UP <- top120_caco100_0to4_UP[which(top120_caco100_0to4_UP$logFC > 0), ] [1:120,]
results_caco100_0to4_UP <- results@.Data[,"caco100_at_0to4"]
results_caco100_0to4_UP <- results_caco100_0to4_UP[names(results_caco100_0to4_UP) %in% top120_caco100_0to4_UP$GeneID]

top120_org_0to4_UP <- topTable(eb, coef = "org_at_0to4", number = Inf)
top120_org_0to4_UP <- top120_org_0to4_UP[which(top120_org_0to4_UP$logFC > 0), ] [1:120,]
results_org_0to4_UP <- results@.Data[,"org_at_0to4"]
results_org_0to4_UP <- results_org_0to4_UP[names(results_org_0to4_UP) %in% top120_org_0to4_UP$GeneID]

top120_caco50_0to24_UP <- topTable(eb, coef = "caco50_at_0to24", number = Inf)
top120_caco50_0to24_UP <- top120_caco50_0to24_UP[which(top120_caco50_0to24_UP$logFC > 0), ] [1:120,]
results_caco50_0to24_UP <- results@.Data[,"caco50_at_0to24"]
results_caco50_0to24_UP <- results_caco50_0to24_UP[names(results_caco50_0to24_UP) %in% top120_caco50_0to24_UP$GeneID]

top120_caco100_0to24_UP <- topTable(eb, coef = "caco100_at_0to24", number = Inf)
top120_caco100_0to24_UP <- top120_caco100_0to24_UP[which(top120_caco100_0to24_UP$logFC > 0), ] [1:120,]
results_caco100_0to24_UP <- results@.Data[,"caco100_at_0to24"]
results_caco100_0to24_UP <- results_caco100_0to24_UP[names(results_caco100_0to24_UP) %in% top120_caco100_0to24_UP$GeneID]

top120_org_0to24_UP <- topTable(eb, coef = "org_at_0to24", number = Inf)
top120_org_0to24_UP <- top120_org_0to24_UP[which(top120_org_0to24_UP$logFC > 0), ] [1:120,]
results_org_0to24_UP <- results@.Data[,"org_at_0to24"]
results_org_0to24_UP <- results_org_0to24_UP[names(results_org_0to24_UP) %in% top120_org_0to24_UP$GeneID]



# Visualise the IBD pathway.
pathway <- "05321"
pathview.data <- cbind(caco50_at_0to24 = results_caco50_0to24_UP, 
                       caco100_at_0to24 = results_caco100_0to24_UP, 
                       org_at_0to24 = results_org_0to24_UP)
pv.out <- pathview(gene.data = pathview.data, pathway.id = pathway,
                   species = "hsa", out.suffix = "IBD_0to24_120UP", kegg.native = T, 
                   dsicrete=list(gene=TRUE, cpd=FALSE), 
                   bins = list(gene = 4, cpd= 10), 
                   low = list(gene = "#1e90ff", cpd = "blue"), 
                   mid = list(gene = "gray", cpd = "gray"), 
                   high = list(gene = "orange", cpd = "yellow"))

# HEATMAPS -----------------------------------------------------------

#, Patient = factor(label$patient)

makeHeat <- function(group1, group2, interaction, patient=F) {
  samples <- label$group == group1 | label$group == group2
  DEgenes <- topTable(eb, coef = interaction, number = 5000, p.value = 0.05)$GeneID
  SelectGenes <- rownames(logCPM) %in% DEgenes
  plotData <- logCPM[SelectGenes,samples]
  #plotData <- plotData[,order(colnames(plotData))]
  if (patient) {
    group_annotation <- data.frame(Group = factor(label$group, levels = c(group1, group2)),
                                   Patient = factor(label$patient)) %>%  
      set_rownames(label$SampleID)
  } else {
    group_annotation <- data.frame(Group = factor(label$group, levels = c(group1, group2))) %>%  
      set_rownames(label$SampleID)
  }
  
  pheatmap(plotData, show_rownames = FALSE, show_colnames = FALSE,
           width = 8, height = 8, scale = "row", fontsize = 14,
           annotation_col = group_annotation,
           cluster_cols = FALSE,
           filename=sprintf("plots/heatmaps/%s.pdf", interaction))
}

makeHeatcomplx <- function(group1, group2, group3, group4, interaction, patient=F) {
  samples <- label$group == group1 | label$group == group2 | label$group == group3 | label$group == group4
  DEgenes <- topTable(eb, coef = interaction, number = 5000, p.value = 0.05)$GeneID
  SelectGenes <- rownames(logCPM) %in% DEgenes
  plotData <- logCPM[SelectGenes,samples]
  #plotData <- plotData[,order(colnames(plotData))]
  if (patient) {
    group_annotation <- data.frame(Group = factor(label$group, levels = c(group1, group3, group2, group4)),
                                   Patient = factor(label$patient)) %>%  
      set_rownames(label$SampleID)
  } else {
    group_annotation <- data.frame(Group = factor(label$group, levels = c(group1, group3, group2, group4))) %>%  
      set_rownames(label$SampleID)
  }
  
  pheatmap(plotData, show_rownames = FALSE, show_colnames = FALSE,
           width = 8, height = 8, scale = "row", fontsize = 14,
           annotation_col = group_annotation,
           cluster_cols = FALSE,
           filename=sprintf("plots/heatmaps/%s.pdf", interaction))
}

makeHeat(group1="caco.50.0", group2="caco.50.4", interaction="caco50_at_0to4")
makeHeat(group1="caco.50.4", group2="caco.50.24", interaction="caco50_at_4to24")
makeHeat(group1="caco.50.0", group2="caco.50.24", interaction="caco50_at_0to24")
makeHeat(group1="caco.100.0", group2="caco.100.4", interaction="caco100_at_0to4")
makeHeat(group1="caco.100.4", group2="caco.100.24", interaction="caco100_at_4to24")
makeHeat(group1="caco.100.0", group2="caco.100.24", interaction="caco100_at_0to24")
makeHeat(group1="org.0", group2="org.4", interaction="org_at_0to4", patient = T)
makeHeat(group1="org.4", group2="org.24", interaction="org_at_4to24", patient = T)
makeHeat(group1="org.0", group2="org.24", interaction="org_at_0to24", patient = T)

makeHeat(group1="caco.50.0", group2="caco.100.0", interaction="caco50_vs_caco100_at_0")
makeHeat(group1="caco.50.4", group2="caco.100.4", interaction="caco50_vs_caco100_at_4")
makeHeat(group1="caco.50.24", group2="caco.100.24", interaction="caco50_vs_caco100_at_24")
makeHeat(group1="caco.50.0", group2="org.0", interaction="caco50_vs_org_at_0", patient = T)
makeHeat(group1="caco.50.4", group2="org.4", interaction="caco50_vs_org_at_4", patient = T)
makeHeat(group1="caco.50.24", group2="org.24", interaction="caco50_vs_org_at_24", patient = T)
makeHeat(group1="caco.100.0", group2="org.0", interaction="caco100_vs_org_at_0", patient = T)
makeHeat(group1="caco.100.4", group2="org.4", interaction="caco100_vs_org_at_4", patient = T)
makeHeat(group1="caco.100.24", group2="org.24", interaction="caco100_vs_org_at_24", patient = T)

makeHeatcomplx(group1="caco.50.0", group2="caco.100.0", group3="caco.50.4", group4="caco.100.4", interaction="caco50_vs_caco100_at_0to4")
makeHeatcomplx(group1="caco.50.4", group2="caco.100.4", group3="caco.50.24", group4="caco.100.24", interaction="caco50_vs_caco100_at_4to24")
makeHeatcomplx(group1="caco.50.0", group2="caco.100.0", group3="caco.50.24", group4="caco.100.24", interaction="caco50_vs_caco100_at_0to24")
makeHeatcomplx(group1="caco.50.0", group2="org.0", group3="caco.50.4", group4="org.4", interaction="caco50_vs_org_at_0to4", patient = T)
makeHeatcomplx(group1="caco.50.4", group2="org.4", group3="caco.50.24", group4="org.24", interaction="caco50_vs_org_at_4to24", patient = T)
makeHeatcomplx(group1="caco.50.0", group2="org.0", group3="caco.50.24", group4="org.24", interaction="caco50_vs_org_at_0to24", patient = T)
makeHeatcomplx(group1="caco.100.0", group2="org.0", group3="caco.100.4", group4="org.4", interaction="caco100_vs_org_at_0to4", patient = T)
makeHeatcomplx(group1="caco.100.4", group2="org.4", group3="caco.100.24", group4="org.24", interaction="caco100_vs_org_at_4to24", patient = T)
makeHeatcomplx(group1="caco.100.0", group2="org.0", group3="caco.100.24", group4="org.24", interaction="caco100_vs_org_at_0to24", patient = T)

# GENAS -------------------------------------------------------------------

pdf(file = sprintf("plots/genas/%s.pdf", colnames(contrasts)[19]), width = 5, height = 5)
genas(eb, coef=c("caco100_at_0to4","caco50_at_0to4"), plot=TRUE, alpha=0.4, subset = "Fpval")
dev.off()
pdf(file = sprintf("plots/genas/%s.pdf", colnames(contrasts)[20]), width = 5, height = 5)
genas(eb, coef=c("caco100_at_4to24","caco50_at_4to24"), plot=TRUE, alpha=0.4, subset = "Fpval")
dev.off()
pdf(file = sprintf("plots/genas/%s.pdf", colnames(contrasts)[21]), width = 5, height = 5)
genas(eb, coef=c("caco100_at_0to24","caco50_at_0to24"), plot=TRUE, alpha=0.4, subset = "Fpval")
dev.off()

pdf(file = sprintf("plots/genas/%s.pdf", colnames(contrasts)[22]), width = 5, height = 5)
genas(eb, coef=c("org_at_0to4","caco50_at_0to4"), plot=TRUE, alpha=0.4, subset = "Fpval")
dev.off()
pdf(file = sprintf("plots/genas/%s.pdf", colnames(contrasts)[23]), width = 5, height = 5)
genas(eb, coef=c("org_at_4to24","caco50_at_4to24"), plot=TRUE, alpha=0.4, subset = "Fpval")
dev.off()
pdf(file = sprintf("plots/genas/%s.pdf", colnames(contrasts)[24]), width = 5, height = 5)
genas(eb, coef=c("org_at_0to24","caco50_at_0to24"), plot=TRUE, alpha=0.4, subset = "Fpval")
dev.off()

pdf(file = sprintf("plots/genas/%s.pdf", colnames(contrasts)[25]), width = 5, height = 5)
genas(eb, coef=c("org_at_0to4","caco100_at_0to4"), plot=TRUE, alpha=0.4, subset = "Fpval")
dev.off()
pdf(file = sprintf("plots/genas/%s.pdf", colnames(contrasts)[26]), width = 5, height = 5)
genas(eb, coef=c("org_at_4to24","caco100_at_4to24"), plot=TRUE, alpha=0.4, subset = "Fpval")
dev.off()
pdf(file = sprintf("plots/genas/%s.pdf", colnames(contrasts)[27]), width = 5, height = 5)
genas(eb, coef=c("org_at_0to24","caco100_at_0to24"), plot=TRUE, alpha=0.4, subset = "Fpval")
dev.off()
# INDIVIDUAL RESPONSE -----------------------------------------------------

plotResponse <- function(data, geneSymbol) {
  p <- ggplot(data, aes(x = group, y = logCPM, min = 0)) +
    geom_jitter(aes(color = patient), shape=16, position=position_jitter(0.2), size = 4) +
    geom_boxplot(outlier.shape = NA, alpha = 0) +
    
    ggtitle(sprintf("%s", geneSymbol)) +
    theme_pubr() +
    theme(plot.title = element_text(hjust = 0.5, size = 30),
          axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size=20),
          legend.title = element_blank()) +
    scale_fill_grey() + 
    scale_color_manual(values=cbPalette_resp, breaks=c("Patient1", "Patient2", "Patient3", "Patient4"))
  return(p)
}
indvResponse <- function(gene) {
  geneSymbol <- mapIds(org.Hs.eg.db::org.Hs.eg.db,
                       keys=gene, # The IDs
                       column="SYMBOL", 
                       keytype="ENTREZID")
  respall <- tibble()
  for (group in c("caco.50", "caco.100", "org")) {
    group1 <- paste0(group, ".0")
    group2 <- paste0(group, ".4")
    group3 <- paste0(group, ".24")
    sample1 <- label$group == group1
    sample2 <- label$group == group2
    sample3 <- label$group == group3
    resp1 <- cbind(tibble(logCPM=logCPM[gene, sample1]), 
                          group = group1, 
                          biggroup = group, 
                          patient = label$patient[match(names(logCPM[gene, sample1]), label$SampleID)])
    resp2 <- cbind(tibble(logCPM=logCPM[gene, sample2]), 
                   group = group2, 
                   biggroup = group, 
                   patient = label$patient[match(names(logCPM[gene, sample2]), label$SampleID)])
    resp3 <- cbind(tibble(logCPM=logCPM[gene, sample3]), 
                   group = group3, 
                   biggroup = group, 
                   patient = label$patient[match(names(logCPM[gene, sample3]), label$SampleID)])
    respall <- rbind(respall, resp1, resp2, resp3)
  }
  respall$group <- factor(respall$group, 
                          levels = stringr::str_sort(unique(respall$group), numeric = TRUE))
  ggsave(filename = sprintf("plots/Responses/%s_all.pdf", geneSymbol), 
         plot = plotResponse(respall, geneSymbol), device = "pdf", scale = 1)
  ggsave(filename = sprintf("plots/Responses/%s_caco50.pdf", geneSymbol), 
         plot = plotResponse(filter(respall, biggroup == "caco.50"), geneSymbol), device = "pdf", scale = .65)
  ggsave(filename = sprintf("plots/Responses/%s_caco100.pdf", geneSymbol), 
         plot = plotResponse(filter(respall, biggroup == "caco.100"), geneSymbol), device = "pdf", scale = .65)
  ggsave(filename = sprintf("plots/Responses/%s_org.pdf", geneSymbol), 
         plot = plotResponse(filter(respall, biggroup == "org"), geneSymbol), device = "pdf", scale = .65)
  
}

indvResponse("3569")
indvResponse("7124")

# BARCODE PLOTS -----------------------------------------------------------

IBD_stats <- data_frame(genes=names(eb_IBD$coefficients[,1]), weights=eb_IBD$coefficients[,1])

for (i in seq(1:length(eb$t[1,]))) {
  name <- colnames(eb)[i]
  t_stats <- eb$t[,i]
  genes <- names(eb$t[,i])
  plot_stats <- data_frame(genes, t_stats) %>% 
    left_join(., IBD_stats, by = "genes") %>% 
    replace_na(list(weights = 0))
  pdf(file = sprintf(sprintf("plots/barcode/%s.pdf", name)))
  barcodeplot(statistics = plot_stats$t_stats, 
            weights.label = "logFC", 
            gene.weights = plot_stats$weights,
            main = sprintf("%s", name))
  dev.off()

}

# VENN DIAGRAMES ----------------------------------------------------------

for (i in seq(1:3)) {
  name <- colnames(eb$coefficients)[i] %>% substr(11, nchar(.))
  caco_name <- colnames(eb$coefficients)[i]
  org_name <- colnames(eb$coefficients)[i+6]
  caco_genes_UP <- names(which(results[,i] == 1))
  caco_genes_DOWN <- names(which(results[,i] == -1))
  org_genes_UP <- names(which(results[,i+6] == 1))
  org_genes_DOWN <- names(which(results[,i+6] == -1))
  
  data_list_UP <- vector(mode="list", length=3)
  names(data_list_UP) <- paste0(c(caco_name, org_name, "IBD"), "_UP")
  data_list_DOWN <- vector(mode="list", length=3)
  names(data_list_DOWN) <- paste0(c(caco_name, org_name, "IBD"), "_DOWN")
  data_list_UP[[1]] <- caco_genes_UP
  data_list_UP[[2]] <- org_genes_UP
  data_list_UP[[3]] <- IBD_UP
  data_list_DOWN[[1]] <- caco_genes_DOWN
  data_list_DOWN[[2]] <- org_genes_DOWN
  data_list_DOWN[[3]] <- IBD_DOWN
  venn_UP <- Venn(data_list_UP)
  venn_DOWN <- Venn(data_list_DOWN)
  
  pdf(file = sprintf(sprintf("plots/venn/genes_%s_UP.pdf", name)))
  plot(venn_UP, doWeights = TRUE, doEuler = TRUE)
  dev.off()
  pdf(file = sprintf(sprintf("plots/venn/genes_%s_DOWN.pdf", name)))
  plot(venn_DOWN, doWeights = TRUE, doEuler = TRUE)
  dev.off()
}

for (i in seq(1:3)) {
  name <- colnames(eb$coefficients)[i] %>% substr(11, nchar(.))
  caco_name <- colnames(eb$coefficients)[i]
  org_name <- colnames(eb$coefficients)[i+6]
  
  g_caco <- goana(eb, coef=i, trend = TRUE)
  g_caco_up_adj <- p.adjust(g_caco$P.Up)
  g_caco_UP <- rownames(g_caco[g_caco_up_adj<0.05,])
  g_caco_down_adj <- p.adjust(g_caco$P.Down)
  g_caco_DOWN <- rownames(g_caco[g_caco_down_adj<0.05,])
  
  g_org <- goana(eb, coef=i+6, trend = TRUE)
  g_org_up_adj <- p.adjust(g_org$P.Up)
  g_org_UP <- rownames(g_org[g_org_up_adj<0.05,])
  g_org_down_adj <- p.adjust(g_org$P.Down)
  g_org_DOWN <- rownames(g_org[g_org_down_adj<0.05,])
  
  data_list_UP <- vector(mode="list", length=3)
  names(data_list_UP) <- paste0(c(caco_name, org_name, "IBD"), "_UP")
  data_list_DOWN <- vector(mode="list", length=3)
  names(data_list_DOWN) <- paste0(c(caco_name, org_name, "IBD"), "_DOWN")
  data_list_UP[[1]] <- g_caco_UP
  data_list_UP[[2]] <- g_org_UP
  data_list_UP[[3]] <- IBD_UP_GO
  data_list_DOWN[[1]] <- g_caco_DOWN
  data_list_DOWN[[2]] <- g_org_DOWN
  data_list_DOWN[[3]] <- IBD_DOWN_GO
  
  if(!isEmpty(g_caco_UP) & !isEmpty(g_org_UP)) {
    venn_UP <- Venn(data_list_UP)
    pdf(file = sprintf(sprintf("plots/venn/GO1_%s_UP.pdf", name)))
    plot(venn_UP, doWeights = TRUE, doEuler = TRUE)
    dev.off()
  }
  if(!isEmpty(g_caco_DOWN) & !isEmpty(g_org_DOWN)) {
    venn_DOWN <- Venn(data_list_DOWN)
    pdf(file = sprintf(sprintf("plots/venn/GO1_%s_DOWN.pdf", name)))
    plot(venn_DOWN, doWeights = TRUE, doEuler = TRUE)
    dev.off()
  }
}

for (i in seq(1:3)) {
  name <- colnames(eb$coefficients)[i] %>% substr(11, nchar(.))
  caco_name <- colnames(eb$coefficients)[i]
  org_name <- colnames(eb$coefficients)[i+6]
  
  k_caco <- kegga(eb, coef=i, trend = TRUE)
  k_caco_up_adj <- p.adjust(k_caco$P.Up)
  k_caco_UP <- rownames(k_caco[k_caco_up_adj<0.05,])
  k_caco_down_adj <- p.adjust(k_caco$P.Down)
  k_caco_DOWN <- rownames(k_caco[k_caco_down_adj<0.05,])
  
  k_org <- kegga(eb, coef=i+6, trend = TRUE)
  k_org_up_adj <- p.adjust(k_org$P.Up)
  k_org_UP <- rownames(k_org[k_org_up_adj<0.05,])
  k_org_down_adj <- p.adjust(k_org$P.Down)
  k_org_DOWN <- rownames(k_org[k_org_down_adj<0.05,])
  
  data_list_UP <- vector(mode="list", length=3)
  names(data_list_UP) <- paste0(c(caco_name, org_name, "IBD"), "_UP")
  data_list_DOWN <- vector(mode="list", length=3)
  names(data_list_DOWN) <- paste0(c(caco_name, org_name, "IBD"), "_DOWN")
  data_list_UP[[1]] <- k_caco_UP
  data_list_UP[[2]] <- k_org_UP
  data_list_UP[[3]] <- IBD_UP_KEGG
  data_list_DOWN[[1]] <- k_caco_DOWN
  data_list_DOWN[[2]] <- k_org_DOWN
  data_list_DOWN[[3]] <- IBD_DOWN_KEGG
  
  if(!isEmpty(k_caco_UP) & !isEmpty(k_org_UP)) {
    venn_UP <- Venn(data_list_UP)
    pdf(file = sprintf(sprintf("plots/venn/KEGG1_%s_UP.pdf", name)))
    plot(venn_UP, doWeights = TRUE, doEuler = TRUE)
    dev.off()
  }
  if(!isEmpty(k_caco_DOWN) & !isEmpty(k_org_DOWN)) {
    venn_DOWN <- Venn(data_list_DOWN)
    pdf(file = sprintf(sprintf("plots/venn/KEGG1_%s_DOWN.pdf", name)))
    plot(venn_DOWN, doWeights = TRUE, doEuler = TRUE)
    dev.off()
  }
}

# FISHER TESTS -------------------------------------------------------------

all_genes <- rownames(eb)
all_go <- rownames(goana(eb, coef=1))
all_kegg <- rownames(kegga(eb, coef=1))
IBD_UP_in_all_genes <- all_genes %in% IBD_UP
IBD_DOWN_in_all_genes <- all_genes %in% IBD_DOWN
IBD_UP_in_all_go <- all_go %in% IBD_UP_GO
IBD_DOWN_in_all_go <- all_go %in% IBD_DOWN_GO
IBD_UP_in_all_kegg <- all_kegg %in% IBD_UP_KEGG
IBD_DOWN_in_all_kegg <- all_kegg %in% IBD_DOWN_KEGG

# Fisher tests for gene set enrichments.
f_tests <- vector("list", 18)
for (i in seq(1:9)) {
  name <- colnames(eb$coefficients)[i]
  DE_UP <- names(which(results[,i] == 1))
  DE_UP <- all_genes %in% DE_UP
  DE_DOWN <- names(which(results[,i] == -1))
  DE_DOWN <- all_genes %in% DE_DOWN
  t_UP <- table(DE_UP, IBD_UP_in_all_genes)
  t_DOWN <- table(DE_DOWN, IBD_DOWN_in_all_genes)
  f_UP <- fisher.test(t_UP)
  f_UP$data.name <- name
  capture.output(f_UP, file = sprintf('fishers/%s_UP_gene.txt', name)) 
  f_DOWN <- fisher.test(t_DOWN)
  f_DOWN$data.name <- name
  capture.output(f_DOWN, file = sprintf('fishers/%s_DOWN_gene.txt', name))
  f_tests[[i]] <- c(f_UP$data.name, 
                    "UP", 
                    f_UP$estimate[[1]], 
                    f_UP$conf.int[[1]],
                    f_UP$conf.int[[2]],
                    f_UP$p.value)
  f_tests[[19-i]] <- c(f_DOWN$data.name, 
                    "DOWN", 
                    f_DOWN$estimate[[1]], 
                    f_DOWN$conf.int[[1]],
                    f_DOWN$conf.int[[2]],
                    f_DOWN$p.value)
  
}

f_tests_go <- vector("list", 18)
for (i in seq(1:9)) {
  name <- colnames(eb$coefficients)[i]
  g <- goana(eb, coef=i, trend = TRUE)
  g_up_adj <- p.adjust(g$P.Up)
  g_UP <- rownames(g[g_up_adj<0.05,])
  g_UP <- all_go %in% g_UP
  g_down_adj <- p.adjust(g$P.Down)
  g_DOWN <- rownames(g[g_down_adj<0.05,])
  g_DOWN <- all_go %in% g_DOWN
  t_UP <- table(g_UP, IBD_UP_in_all_go)
  t_DOWN <- table(g_DOWN, IBD_DOWN_in_all_go)
  if (len(t_DOWN)!=2) {
    t_DOWN <- rbind(t_DOWN, c(0,0))
  }
  if (len(t_UP)!=2) {
    t_UP <- rbind(t_UP, c(0,0))
  }
  f_UP <- fisher.test(t_UP)
  f_UP$data.name <- name
  capture.output(f_UP, file = sprintf('fishers/GO/%s_UP.txt', name)) 
  f_DOWN <- fisher.test(t_DOWN)
  f_DOWN$data.name <- name
  capture.output(f_DOWN, file = sprintf('fishers/GO/%s_DOWN.txt', name))
  f_tests_go[[i]] <- c(f_UP$data.name, 
                    "UP", 
                    f_UP$estimate[[1]], 
                    f_UP$conf.int[[1]],
                    f_UP$conf.int[[2]],
                    f_UP$p.value)
  f_tests_go[[19-i]] <- c(f_DOWN$data.name, 
                       "DOWN", 
                       f_DOWN$estimate[[1]], 
                       f_DOWN$conf.int[[1]],
                       f_DOWN$conf.int[[2]],
                       f_DOWN$p.value)
  
}

f_tests_kegg <- vector("list", 18)
for (i in seq(1:9)) {
  name <- colnames(eb$coefficients)[i]
  k <- kegga(eb, coef=i, trend = TRUE)
  k_up_adj <- p.adjust(k$P.Up)
  k_UP <- rownames(k[k_up_adj<0.05,])
  k_UP <- all_kegg %in% k_UP
  k_down_adj <- p.adjust(k$P.Down)
  k_DOWN <- rownames(k[k_down_adj<0.05,])
  k_DOWN <- all_kegg %in% k_DOWN
  t_UP <- table(k_UP, IBD_UP_in_all_kegg)
  t_DOWN <- table(k_DOWN, IBD_DOWN_in_all_kegg)
  if (len(t_DOWN)!=2) {
    t_DOWN <- rbind(t_DOWN, c(0,0))
  }
  if (len(t_UP)!=2) {
    t_UP <- rbind(t_UP, c(0,0))
  }
  f_UP <- fisher.test(t_UP)
  f_UP$data.name <- name
  capture.output(f_UP, file = sprintf('fishers/KEGG/%s_UP.txt', name)) 
  f_DOWN <- fisher.test(t_DOWN)
  f_DOWN$data.name <- name
  capture.output(f_DOWN, file = sprintf('fishers/KEGG/%s_DOWN.txt', name))
  f_tests_kegg[[i]] <- c(f_UP$data.name, 
                    "UP", 
                    f_UP$estimate[[1]], 
                    f_UP$conf.int[[1]],
                    f_UP$conf.int[[2]],
                    f_UP$p.value)
  f_tests_kegg[[19-i]] <- c(f_DOWN$data.name, 
                       "DOWN", 
                       f_DOWN$estimate[[1]], 
                       f_DOWN$conf.int[[1]],
                       f_DOWN$conf.int[[2]],
                       f_DOWN$p.value)
  
}

plot_fisher <- function(list) {
  f_tests <- do.call(rbind.data.frame, list)
  names(f_tests) <- c("group", "change", "OR", "CIup", "CIdown", "p")
  f_tests <- f_tests %>% mutate(group = factor(group, levels = c("caco50_at_0to4", "caco50_at_4to24", "caco50_at_0to24", "caco100_at_0to4", "caco100_at_4to24", "caco100_at_0to24", "org_at_0to4", "org_at_4to24", "org_at_0to24")),
                                change = factor(change, levels = c("UP", "DOWN")) ,
                                OR = as.numeric(OR),
                                CIup = as.numeric(CIup),
                                CIdown = as.numeric(CIdown),
                                p = cut(as.numeric(p), 
                                        breaks = c(-Inf, 0.05, Inf), 
                                        labels = c("Significant", "Not-Significant")))
  pd <- position_dodge(0.1)
  breaks <- c("caco50_at_0to4", "caco50_at_4to24", "caco50_at_0to24", "org_at_0to4", "org_at_4to24", "org_at_0to24")
  ggplot(f_tests[f_tests$group %in% breaks,], aes(x=group, y=OR, fill=p)) + facet_wrap(~change, ncol = 1) +
    geom_errorbar(aes(ymin=CIdown, ymax=CIup), colour="black", width=.1, position=pd) + 
    geom_point(position=pd, size=3, shape=21) + # 21 is filled circle
    geom_abline(slope = 0, intercept = 1, color = "red", linetype = 4) +
    xlab("Contrast") +
    ylab("Odds Ratio (OR)") +
    scale_fill_manual(values = cbPalette_resp) +  
    guides(fill=guide_legend(title="p-value")) +
    #scale_y_continuous(breaks = seq(0, 5, by = 1),
    #                  limits = c(0, 5)) +
    scale_y_continuous(limits = c(0, 150)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          text = element_text(size=20),
          legend.position = "top")
}

plot_fisher(f_tests)
plot_fisher(f_tests_go)
plot_fisher(f_tests_kegg)

# VENN DIAGRAMES 150 ----------------------------------------------------------

for (i in seq(1:3)) {
  name <- colnames(eb$coefficients)[i] %>% substr(11, nchar(.))
  caco_name <- colnames(eb$coefficients)[i]
  org_name <- colnames(eb$coefficients)[i+6]
  
  top150_caco <- topTable(eb, coef = i, number = Inf)
  top150_caco <- top150_caco[which(top150_caco$logFC > 0), ] [1:150,]
  top150_org <- topTable(eb, coef = i+6, number = Inf)
  top150_org <- top150_org[which(top150_org$logFC > 0), ] [1:150,]
  
  g_caco <- goana(top150_caco$GeneID, trend = TRUE, universe = rownames(eb))
  g_caco_adj <- p.adjust(g_caco$P.DE)
  g_caco <- rownames(g_caco[g_caco_adj<0.05,])
  
  g_org <- goana(top150_org$GeneID, trend = TRUE, universe = rownames(eb))
  g_org_adj <- p.adjust(g_org$P.DE)
  g_org <- rownames(g_org[g_org_adj<0.05,])
  
  data_list <- vector(mode="list", length=3)
  names(data_list) <- paste0(c(caco_name, org_name, "IBD"), "_UP")
  data_list[[1]] <- g_caco
  data_list[[2]] <- g_org
  data_list[[3]] <- IBD_UP_GO
  
  if(!isEmpty(g_caco) & !isEmpty(g_org)) {
    venn <- Venn(data_list)
    pdf(file = sprintf(sprintf("plots/venn/150/GO_%s_UP.pdf", name)))
    plot(venn, doWeights = TRUE, doEuler = TRUE)
    dev.off()
  }
}

for (i in seq(1:3)) {
  name <- colnames(eb$coefficients)[i] %>% substr(11, nchar(.))
  caco_name <- colnames(eb$coefficients)[i]
  org_name <- colnames(eb$coefficients)[i+6]
  
  top150_caco <- topTable(eb, coef = i, number = Inf)
  top150_caco <- top150_caco[which(top150_caco$logFC > 0), ] [1:150,]
  top150_org <- topTable(eb, coef = i+6, number = Inf)
  top150_org <- top150_org[which(top150_org$logFC > 0), ] [1:150,]
  
  k_caco <- kegga(top150_caco$GeneID, trend = TRUE, universe = rownames(eb))
  k_caco_adj <- p.adjust(k_caco$P.DE)
  k_caco <- rownames(k_caco[k_caco_adj<0.05,])
  
  k_org <- kegga(top150_org$GeneID, trend = TRUE, universe = rownames(eb))
  k_org_adj <- p.adjust(k_org$P.DE)
  k_org <- rownames(k_org[k_org_adj<0.05,])
  
  data_list <- vector(mode="list", length=3)
  names(data_list) <- paste0(c(caco_name, org_name, "IBD"), "_UP")
  data_list[[1]] <- k_caco
  data_list[[2]] <- k_org
  data_list[[3]] <- IBD_UP_KEGG
  
  if(!isEmpty(k_caco) & !isEmpty(k_org)) {
    venn <- Venn(data_list)
    pdf(file = sprintf(sprintf("plots/venn/150/KEGG1_%s_UP.pdf", name)))
    plot(venn, doWeights = FALSE, doEuler = FALSE)
    dev.off()
  }
}

# FISHER TESTS 150 -------------------------------------------------------------

all_go <- rownames(goana(eb, coef=1))
all_kegg <- rownames(kegga(eb, coef=1))
IBD_UP_in_all_go <- all_go %in% IBD_UP_GO
IBD_UP_in_all_kegg <- all_kegg %in% IBD_UP_KEGG

# Fisher tests for gene set enrichments.

f_tests_go_150 <- vector("list", 18)
for (i in c(1,3,4,6,7,9)) {
  name <- colnames(eb$coefficients)[i]
  top150 <- topTable(eb, coef = i, number = Inf)
  top150 <- top150[which(top150$logFC > 0), ] [1:150,]
  g <- goana(top150$GeneID, trend = TRUE, universe = rownames(eb))
  g_adj <- p.adjust(g$P.DE)
  g <- rownames(g[g_adj<0.05,])
  g <- all_go %in% g
  t <- table(g, IBD_UP_in_all_go)
  if (len(t)!=2) {
    t <- rbind(t, c(0,0))
  }
  f <- fisher.test(t)
  f$data.name <- name
  capture.output(f, file = sprintf('fishers/150/GO/%s_UP.txt', name)) 
  f_tests_go_150[[i]] <- c(f$data.name, 
                       "UP", 
                       f$estimate[[1]], 
                       f$conf.int[[1]],
                       f$conf.int[[2]],
                       f$p.value)
  
}

f_tests_kegg_150 <- vector("list", 18)
for (i in c(1,3,4,6,7,9)) {
  name <- colnames(eb$coefficients)[i]
  top150 <- topTable(eb, coef = i, number = Inf)
  top150 <- top150[which(top150$logFC > 0), ] [1:150,]
  k <- kegga(top150$GeneID, trend = TRUE, universe = rownames(eb))
  k_adj <- p.adjust(k$P.DE)
  k <- rownames(k[k_adj<0.05,])
  k <- all_kegg %in% k
  t <- table(k, IBD_UP_in_all_kegg)
  if (len(t)!=2) {
    t <- rbind(t, c(0,0))
  }
  f <- fisher.test(t)
  f$data.name <- name
  capture.output(f, file = sprintf('fishers/150/KEGG/%s_UP.txt', name)) 
  f_tests_kegg_150[[i]] <- c(f$data.name, 
                       "UP", 
                       f$estimate[[1]], 
                       f$conf.int[[1]],
                       f$conf.int[[2]],
                       f$p.value)
  
}

plot_fisher(f_tests_go_150)
plot_fisher(f_tests_kegg_150)

# HISTOGRAM OF log2FC -----------------------------------------------------

toptab_caco50_at_0to4 <- topTable(eb, coef = "caco50_at_0to4", number = Inf)
p1 <- ggplot(toptab_caco50_at_0to4, aes(x=logFC)) + 
  geom_histogram(binwidth=0.1) +
  theme_bw() + theme(text = element_text(size=15)) +
  ggtitle("caco50_at_0to4")

toptab_org_at_0to4 <- topTable(eb, coef = "org_at_0to4", number = Inf)
p2 <- ggplot(toptab_org_at_0to4, aes(x=logFC)) + 
  geom_histogram(binwidth=0.1) +
  theme_bw() + theme(text = element_text(size=15)) +
  ggtitle("org_at_0to4")

toptab_caco50_at_0to4_top150 <- toptab_caco50_at_0to4[which(toptab_caco50_at_0to4$logFC > 0), ] [1:150,]
p3 <- ggplot(toptab_caco50_at_0to4_top150, aes(x=logFC)) + 
  geom_histogram(binwidth=0.1) +
  theme_bw() + theme(text = element_text(size=15)) +
  ggtitle("caco50_at_0to4_top150")

toptab_org_at_0to4_top150 <- toptab_org_at_0to4[which(toptab_org_at_0to4$logFC > 0), ] [1:150,]
p4 <- ggplot(toptab_org_at_0to4_top150, aes(x=logFC)) + 
  geom_histogram(binwidth=0.1) +
  theme_bw() + theme(text = element_text(size=15)) +
  ggtitle("org_at_0to4_top150")

ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D")) %>% ggsave(filename = "plots/FCcomparisons/logFC_histograms.pdf")



# LOG2 FC comparison ------------------------------------------------------

plotFCcomparison <- function(group1, group2) {
  logFC_caco <- eb$coefficients[, group1, drop = FALSE]
  logFC_org <- eb$coefficients[, group2, drop = FALSE]
  logFC_IBD <- eb_IBD$coefficients[, 1, drop = FALSE]
  keep.genes <- intersect(rownames(logFC_org), rownames(logFC_IBD))
  logFC_caco <- as.data.frame(logFC_caco[rownames(logFC_caco) %in% keep.genes,])
  logFC_org <- as.data.frame(logFC_org[rownames(logFC_org) %in% keep.genes,])
  logFC_IBD <- as.data.frame(logFC_IBD[rownames(logFC_IBD) %in% keep.genes,])
  logFC_IBD <- as.data.frame(logFC_IBD[order(match(rownames(logFC_IBD), rownames(logFC_org))),, drop = FALSE])
  logFC_caco <<- logFC_caco
  logFC_org <<- logFC_org
  logFC_IBD <<- logFC_IBD
  data <- data.frame(caco=logFC_caco[,1], org=logFC_org[,1], IBD=logFC_IBD[,1])
  
  p1 <- ggplot(data, aes(x = caco, y = org)) + 
    scale_fill_viridis(begin = 0, end = 1, option = "C", limit = range(c(10,3500)), breaks=c(10,1000,2000,3500)) +
    geom_hex(aes(fill = ..count..), bins = 20) +
    geom_abline() +
    stat_cor(method = "spearman", label.x = -3, label.y = 8, size = 4.5) +
    coord_fixed(xlim = c(-3, 8), ylim = c(-5, 8)) +
    theme_pubr() +
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size=20),
          legend.text = element_text(size=10)) +
    xlab(group1) +
    ylab(group2)
  p2 <- ggplot(data, aes(x = caco, y = IBD)) + 
    scale_fill_viridis(begin = 0, end = 1, option = "C", limit = range(c(10,3500)), breaks=c(10,1000,2000,3500)) +
    geom_hex(aes(fill = ..count..), bins = 20) +
    geom_abline() +
    stat_cor(method = "spearman", label.x = -3, label.y = 8, size = 4.5) +
    coord_fixed(xlim = c(-3, 8), ylim = c(-5, 8)) +
    theme_pubr() +
    theme(plot.title = element_text(hjust = 0.75),
          text = element_text(size=20),
          legend.text = element_text(size=10)) +
    xlab(group1) +
    ylab("IBD")
  p3 <- ggplot(data, aes(x = org, y = IBD)) + 
    scale_fill_viridis(begin = 0, end = 1, option = "C", limit = range(c(10,3500)), breaks=c(10,1000,2000,3500)) +
    geom_hex(aes(fill = ..count..), bins = 20) +
    geom_abline() +
    stat_cor(method = "spearman", label.x = -3, label.y = 8, size = 4.5) +
    coord_fixed(xlim = c(-3, 8), ylim = c(-5, 8)) +
    theme_pubr() +
    theme(plot.title = element_text(hjust = 0.75),
          text = element_text(size=20),
          legend.text = element_text(size=10)) +
    xlab(group2) +
    ylab("IBD")
  
  ggarrange(p1,                                                 
            ggarrange(p2, p3, ncol = 2, labels = c("B", "C"), legend = "none"), 
            nrow = 2,
            labels = "A", common.legend = TRUE, legend = "top") %>% return()

}

plotFCcomparison("caco50_at_0to4", "org_at_0to4") %>% ggsave(filename = "plots/FCcomparisons/logFC_comparisons_0to4.pdf")

p1 <- ggqqplot(logFC_caco[,1]) + 
  ggtitle("caco50_at_0to4") + 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30)) 
p2 <- ggqqplot(logFC_org[,1]) + 
  ggtitle("org_at_0to4") + 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30))
p3 <- ggqqplot(logFC_IBD[,1]) + 
  ggtitle("IBD") + 
  theme(plot.title = element_text(hjust = 0.5),text = element_text(size=30))
ggarrange(p1, p2, p3, ncol = 2, nrow = 2) %>% ggsave(filename = "plots/FCcomparisons/GaussianCheck.pdf", width = 8, height = 8)
