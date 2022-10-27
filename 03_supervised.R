# Differential Analysis by Status
# Laste update: 10/26/2022

rm(list=ls())
library(limma)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)

DIR <- "~/Documents/datasets/Bariatric/"
P_CUTOFF <- 0.05 / 5000 #Bonferroni

load(paste0(DIR,"calculated_profiles/post_combat.RData"))
dim(expr)

meta$subject_unique <- factor( paste0("s", gsub(" ", "", meta$subject_unique)) )
stopifnot(identical(meta$Accession, colnames(expr))) #important checkpoint

BLOCK <- meta$subject_unique
BLOCK

design <- model.matrix( ~ Status, data=meta)
design

dupcor <- duplicateCorrelation(expr, design, block=BLOCK)
dupcor$consensus.correlation #0.1068558

fit <- lmFit(expr, design=design, correlation=dupcor$consensus.correlation, block=BLOCK, method="robust")
fit <- eBayes(fit, robust=TRUE)

res_dge <- topTable(fit, number=Inf, coef="Status1", adjust.method="bonferroni", sort.by="p")
res_dge$Direction <- NA
res_dge$Direction[res_dge$adj.P.Val < 0.05 & res_dge$logFC > 0] <- "Up"
res_dge$Direction[res_dge$adj.P.Val < 0.05 & res_dge$logFC < 0] <- "Down"
table(res_dge$Direction)
res_dge$color <- factor(res_dge$Direction)

keyvals <- ifelse(res_dge$Direction=="Up", "red", "royalblue")
keyvals[is.na(keyvals)] <- "dimgray"
names(keyvals) <- res_dge$Direction

labSele <- rownames(res_dge)[(res_dge$logFC>0.8 & res_dge$adj.P.Val<0.05) | 
                             (res_dge$logFC < -0.9 & res_dge$adj.P.Val<0.01)]
options(ggrepel.max.overlaps=Inf) #default: 10
EnhancedVolcano(
  res_dge,
  rownames(res_dge),
  "logFC", 
  "P.Value", 
  pCutoff = P_CUTOFF,
  FCcutoff = 0,
  selectLab = labSele,
  colCustom = keyvals,
  colAlpha = 0.4,
  pointSize = 3,
  legendLabSize = 20,
  legendIconSize = 6,
  xlim = c(-1.05, 1.05), 
  ylim = c(0, 16.5),
  title = "",
  subtitle = ""
)

# write.csv(res_dge, file=paste0(DIR,"calculated_profiles/results/dgeBinary.csv"), row.names=TRUE, quote=FALSE)
