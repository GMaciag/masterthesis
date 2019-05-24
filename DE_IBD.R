# IMPORTS AND SETTINGS -----------------------------------------------------------------

library(tidyverse)
library(magrittr)
library(reshape2)
library(limma)
library(edgeR)
library(pheatmap)

options(stringsAsFactors = FALSE)

# CB friendly palette for plotting.
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

setwd("../IBDdata")

EM <- read_tsv("/Users/xvk205/Google Drive/KU/Thesis/IBDdata/counts/genesCounts.tsv")
label <- read_tsv("/Users/xvk205/Google Drive/KU/Thesis/IBDdata/label/design.tsv")
#annotation <- read_csv("~/Google Drive/KU/Thesis/R/counts/counting_annotation_entrez.csv")

# Remove genderless sample.
noG <- label$Sample[label$Gender=="NAN"]
EM <- EM[(colnames(EM)!=noG)]
label <- label[label$Sample!=noG,]

# PREPARE DESIGN -----------------------------------------------------------

label <- as.data.frame(label[1:4])
label <- label[order(match(label$Sample,colnames(EM))),]
rownames(label) <- label$Sample
label$Class <- factor(label$Class)
label$Batch <- factor(label$Batch)
label$Gender <- factor(label$Gender)
label$group <- ifelse(label$Class == "con", "Ctrl", "IBD")
label$group <- factor(label$group)

# Add batch to the design.
design <- model.matrix(~ 0 + group + Batch + Gender, data = label)
colnames(design) <- c(unique(levels(label$group)), unique(levels(label$Batch)[-1]), unique(levels(label$Gender)[-1]))
# Make comparisons between the experimental conditions.
contrasts <- makeContrasts(
  "Ctrl_vs_IBD"=IBD-Ctrl,
  levels=design
)

# FILTERING AND NORMALIZATION ---------------------------------------------

EM_df <- EM[-1] %>% 
  as.data.frame() %>% 
  set_rownames(EM$GeneID) %>% 
  set_colnames(label$Sample)
dge <- DGEList(counts=EM_df, group = as.numeric(factor(label$group)))
# Remove genes that consistently have zero or very low counts.
keep <- filterByExpr(dge, design)
dge <- dge[keep,,keep.lib.sizes=FALSE]
# Apply scale normalization.
dge <- calcNormFactors(dge)
# logCPM normalization.
logCPM <- cpm(dge, log=TRUE, prior.count=1)
# Fit a linear model.
fit <- lmFit(logCPM, design)
# Compute the contrasts.
fit2 <- contrasts.fit(fit, contrasts)
# Apply empirical Bayes smoothing to the standard errors.
eb <- eBayes(fit2, trend=TRUE, robust = TRUE)

# Top genes differential expressed.
results <- decideTests(eb)
results_IBD <- results
eb_IBD <- eb
logCPM_IBD <- logCPM
label_IBD <- label
save(results_IBD, file = "results_IBD.RData")
save(eb_IBD, file = "eb_IBD.RData")
save(logCPM_IBD, label_IBD, file = "logCPM_IBD.RData")

# DE GENES ----------------------------------------------------------------

topgenes <- list()

t <- topTable(eb, coef = 1, number = 10)
t <- add_column(t, GeneSymbol = mapIds(org.Hs.eg.db::org.Hs.eg.db,
                     keys=as.character(rownames(t)), # The IDs
                     column="SYMBOL", 
                     keytype="ENTREZID"), .before = 1)
write.csv(t, 'DEgenes/Ctrl_vs_IBD.csv', row.names = TRUE)

# PREPARE UNIVERSE --------------------------------------------------------

g <- goana(eb, coef=1, trend = TRUE)
k <- kegga(eb, coef=1, trend = TRUE)

g_up_adj <- p.adjust(g$P.Up)
g_down_adj <- p.adjust(g$P.Down)
k_up_adj <- p.adjust(k$P.Up)
k_down_adj <- p.adjust(k$P.Down)

IBD_UP_GO <- rownames(g[g_up_adj<0.05,])
IBD_DOWN_GO <- rownames(g[g_down_adj<0.05,])

IBD_UP_KEGG <- rownames(k[k_up_adj<0.05,])
IBD_DOWN_KEGG <- rownames(k[k_down_adj<0.05,])

IBD_UP <- rownames(results)[results@.Data==1]
IBD_DOWN <- rownames(results)[results@.Data==-1]

save(IBD_UP, IBD_DOWN, IBD_UP_GO, IBD_DOWN_GO, IBD_UP_KEGG, IBD_DOWN_KEGG, file = "IBD_UNIVERSE.RData")
