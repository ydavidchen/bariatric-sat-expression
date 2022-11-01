# Processing the Joined Data Matrix
# Laste update: 10/26/2022

rm(list=ls())
library(sva)
library(ggplot2)
library(ggfortify)
library(pheatmap)

DIR <- "~/Documents/datasets/Bariatric/"

THEME_PC <- theme_bw() +
  theme(axis.text=element_text(size=15,color="black"), 
        axis.title=element_text(size=20,color="black"),
        legend.text=element_text(size=15,color="black"),
        legend.title=element_text(size=15,color="black"),
        legend.position="top")

winsorize <- function(mat, upper, lower) {
  mat[mat > upper] <- upper
  mat[mat < lower] <- lower
  return(mat)
}

## Sample metadata:
meta <- read.csv(paste0(DIR,"calculated_profiles/samplesheet.csv"))
meta <- subset(meta, SAT)
meta$Status <- factor(1 * (meta$Timepoint > 0))
meta$Dataset[meta$Dataset %in% c("GSE29409","GSE29410")] <- "GSE29411" #superseries

table(meta$Dataset, meta$Status)
table(meta$Dataset, meta$Timepoint)

## Expression:
mGse <- read.csv(paste0(DIR,"calculated_profiles/merged_expression.csv"), row.names=1)
mGse <- data.matrix(mGse)
mGse <- mGse[ , colnames(mGse) != "GSM1599916"]
qplot(as.numeric(mGse), xlab="Value", ylab="Freq") + THEME_PC

meta <- subset(meta, Accession %in% colnames(mGse))
mGse <- mGse[ , colnames(mGse) %in% meta$Accession]
mGse <- mGse[ , match(meta$Accession, colnames(mGse))]
stopifnot(identical(meta$Accession, colnames(mGse)))
stopifnot(! anyNA(mGse) )

## PCA by Dataset & Status before ComBat:
pca <- prcomp(t(mGse), center=TRUE, scale.=TRUE)
autoplot(pca, data=meta, colour="Dataset", size=4, alpha=0.5) + THEME_PC #requires ggfortify
autoplot(pca, data=meta, colour="Status", size=4, alpha=0.5) + THEME_PC

## ComBat correction:
stopifnot(identical(meta$Accession, colnames(mGse))) #checkpoint
batch <- meta$Dataset
modcombat <- model.matrix(~ 1, data=meta)

expr <- ComBat(dat=mGse, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
qplot(as.numeric(expr), xlab="Value", ylab="Freq") + THEME_PC

## Scaling across all samples & winsorize:
expr <- t( scale(t(expr)) )
qplot(as.numeric(expr), xlab="Value", ylab="Freq") + THEME_PC
expr <- winsorize(expr, 5, -5)
qplot(as.numeric(expr), xlab="Normalized Expression", ylab="Freq") + THEME_PC

## Remove small number of least variant genes across all datasets, post-ComBat
expr <- expr[order(rowSds(expr), decreasing=TRUE)[1:5000], ]
dim(expr)

## Export:
# save(list=c("expr","meta"), file=paste0(DIR,"calculated_profiles/post_combat.RData"), compress=TRUE)
# write.csv(expr, paste0(DIR,"calculated_profiles/merged_expression_final.csv"), row.names=TRUE, quote=FALSE)
# write.csv(meta, paste0(DIR,"calculated_profiles/samplesheet_final.csv"), row.names=FALSE, quote=FALSE)
