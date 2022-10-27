# Unsupervised Analyses & Visualizations
# Laste update: 10/26/2022

rm(list=ls())
library(matrixStats)
library(ggplot2)
library(ggfortify)
library(gridExtra)
library(reshape2)
library(pheatmap)

DIR <- "~/Documents/datasets/Bariatric/"
COL_ATTR <- c("Accession", "Dataset", "Status")
CL_PARAMS <- c("manhattan", "ward.D")

HEAT_COLS <- colorRampPalette(c("blue","lightgray","red"))(1024)

ANN_COLOR <- list(
  Status = c(`1`="black", `0`="lightgray"),
  Dataset = c(GSE29411="lightgray",GSE42715="gray90", GSE65540="gray60",GSE66921="dimgray",GSE199063="gray30")
)

THEME_PCA <- theme(axis.text=element_text(size=15,color="black"), 
                   axis.title=element_text(size=20,color="black"),
                   legend.text=element_text(size=20,color="black"),
                   legend.title=element_text(size=20,color="black"),
                   legend.position="top")

THEME_BOX <- theme(axis.text=element_text(size=20,color="black"), 
                   axis.title.y=element_text(size=20,color="black"),
                   axis.title.x=element_blank(),
                   legend.text=element_text(size=20,color="black"),
                   legend.title=element_text(size=20,color="black"),
                   legend.position="top")

run2x2 <- function(var1, var2, data, flip1=FALSE, flip2=FALSE, test=FALSE) {
  cTab <- table(data[ , var1], data[ , var2])
  if(flip1) cTab <- cTab[c(2,1), ]
  if(flip2) cTab <- cTab[ , c(2,1)]
  print(cTab)
  if(test) print(fisher.test(cTab))
}

load(paste0(DIR,"calculated_profiles/post_combat.RData"))

# ----------------- Principal Component Analysis -----------------
wrapper_pca <- function(dat, pltSlope=NULL) {
  ## Execute PCA & format results:
  pca_combat <- prcomp(dat, center=TRUE, scale.=TRUE)
  resPca <- pca_combat$x[ , c(1,2)]
  resPca <- scale(resPca) #columnwise
  resPca <- as.data.frame(resPca)
  resPca$Accession <- rownames(resPca)
  resPca <- merge(meta[,COL_ATTR], resPca, by="Accession")
  
  ## 2D visualization:
  print(
    ggplot(resPca, aes(PC1, PC2, color=Dataset)) + 
    geom_point(size=3.5, alpha=0.75) + 
    theme_bw() + 
    THEME_PCA
  )
  
  pltS <- ggplot(resPca, aes(PC1, PC2, color=Status)) + 
    geom_point(size=3.5, alpha=0.75) + 
    scale_color_brewer(palette="Dark2") +
    theme_bw() + 
    THEME_PCA
  
  if(!is.null(pltSlope)) pltS <- pltS+geom_abline(slope=pltSlope, intercept=0, linetype="dashed")
  
  ## Boxplot:
  mResPca <- melt(resPca)
  pltB <- ggplot(mResPca, aes(variable, value, fill=Status)) +
    geom_boxplot(outlier.colour=NA, alpha=0.5) +
    geom_jitter(aes(color=Status), alpha=0.5, position=position_jitterdodge()) +
    scale_color_brewer(palette="Dark2") +
    scale_fill_brewer(palette="Dark2") +
    theme_bw() +
    THEME_BOX

  print( gridExtra::grid.arrange(grobs=list(pltS,pltB), widths=c(0.62,0.38)) )
  
  print( t.test(resPca$PC1 ~ resPca$Status) )
  print( t.test(resPca$PC2 ~ resPca$Status) )
  return(resPca)
}

pcaRes <- wrapper_pca(t(expr), -1)

## Decision function: PC2 = -PC1
pcaRes$Group <- as.integer(pcaRes$PC2 > -pcaRes$PC1)
table(pcaRes$Group)
run2x2("Group","Status", pcaRes, FALSE, TRUE, TRUE)

# ------------------------- Heatmap & Association Tests on Most Variable -------------------------
## Variance distribution:
sampSds <- rowSds(expr)
summary(sampSds)
plot(sort(sampSds, decreasing=TRUE), ylab="Inter-sample SD")
thresh <- sort(sampSds, decreasing=TRUE)[500]
thresh

eSub <- expr[sampSds >= thresh, ]

hm_annot <- data.frame(
  row.names = meta$Accession,
  Status = meta$Status,
  Dataset = meta$Dataset
)

wrapper_hclust <- function(mat, num_cl=2, rowname="Accession") {
  ph <- pheatmap(
    mat,
    clustering_distance_rows = CL_PARAMS[1],
    clustering_distance_cols = CL_PARAMS[1],
    clustering_method = CL_PARAMS[2],
    cutree_cols = num_cl,
    show_rownames = FALSE,
    show_colnames = FALSE,
    fontsize = 15,
    annotation_col = hm_annot, 
    annotation_colors = ANN_COLOR,
    color = HEAT_COLS,
    border_color = NA
  )
  
  ## Tree only:
  myHcl <- ph$tree_col #hclust object
  plot(myHcl, labels=FALSE)
  rect.hclust(myHcl, k=num_cl)
  
  membership <- data.frame(cutree(myHcl, k=num_cl))
  colnames(membership) <- "Cluster"
  if(! is.null(rowname)) {
    membership[rowname] <- rownames(membership)
    rownames(membership) <- NULL
  }
  
  ## Heatmap:
  plot.new()
  print(ph)
  return(membership[ , c(2,1)])
}

## Test cluster membership:
res_cl <- wrapper_hclust(eSub, 2)
res_cl <- merge(res_cl, meta, by="Accession")

table(res_cl$Cluster) #rename cluster if necessary
res_cl$ClusterRenamed <-  ifelse(res_cl$Cluster==2, "High", "Low")
run2x2("ClusterRenamed","Status", res_cl, TRUE, TRUE, TRUE)
