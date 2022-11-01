# Cleaning Individual Datasets & First Joining
# Execution: Rscript <this script>
# Laste update: 10/26/2022

rm(list=ls())
library(data.table)
library(matrixStats)
library(biomaRt)

## Constants:
DIR <- "~/Documents/datasets/Bariatric/"
CHROMS <- c(paste0("chr", 1:22), "chrX")
STRNA <- c("", "---")
MART <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
BM_ATTR <- c("ensembl_gene_id", "hgnc_symbol", "refseq_mrna", "chromosome_name")

## Helpers & callbacks:
intersect_all <- function(a,b,...) Reduce(intersect, list(a,b,...))
cbMean <- function(x) mean(as.numeric(x), na.rm=TRUE)

## Wrappers:
custom_prep <- function(mat, rowname="ID_REF", pLeastVar=0.001) {
  #'@description Wrapper to format matrix & filter least variable per dataset
  #'@param pLeastVar Proportion of least variable genes to filter
  
  ## Format-Standardize:
  rownames(mat) <- as.character(unlist(mat[ , rowname]))
  mat <- mat[ , colnames(mat) != rowname]
  mat <- data.matrix(mat)
  hist(mat, main="", xlab="",  ylab="")
  
  ## Filter flat-zero genes:
  print(paste("Number of annotated autosomal genes:", nrow(mat)))
  
  mat <- mat[rowSums(mat) > 0, ]
  print(paste("Number of genes after filtering all 0:", nrow(mat)))
  
  ## (Optional) Filter least variable genes:
  if(pLeastVar > 0) {
    sampSds <- rowSds(mat)
    sdThreh <- quantile(sampSds, probs=pLeastVar)
    mat <- mat[sampSds >= sdThreh, ]
    print(paste("Number of genes after filtering additional", 100*pLeastVar, "% invariant:", nrow(mat)))
  }
  return(mat)
}

custom_proc <- function(mat, features=NULL, normBySamp=TRUE, bound=0) {
  if(! is.null(features)) mat <- mat[features, ]
  if(normBySamp) mat <- t( scale(t(mat)) )
  if(!is.null(bound) & bound > 0) {
    mat[mat < -bound] <- -bound
    mat[mat > bound] <- bound
  }
  hist(mat, main="", xlab="", ylab="")
  return(mat)
}

custom_cbind <- function(x, xnew) {
  #'@description Wrapper to iteratively bind datasets into one
  if(!is.null(x)) stopifnot(identical(rownames(x), rownames(xnew)))
  if(any(colnames(xnew) %in% colnames(x))) stop("Some columns already in X!")
  return(cbind(x, xnew))
}

# ---------------------------- Process Microarray Platform Annotation Files ---------------------------- 
## GSE29411 superseries:
gp7020 <- fread(paste0(DIR,"platforms/GPL7020.annot.gz"), data.table=FALSE, na.strings=STRNA)
gp7020 <- subset(gp7020, !is.na(`Gene symbol`) & !grepl("///", `Gene symbol`))

gp7020$seqname <- gsub("Chromosome ", "chr", gp7020$`Chromosome annotation`)
gp7020$seqname <- lapply(strsplit(gp7020$seqname, ","), function(x) x[1])
gp7020$seqname <- gsub("///.*", "", gp7020$seqname)
gp7020 <- subset(gp7020, seqname %in% CHROMS)

## GSE42715
gp5715 <- fread(paste0(DIR,"platforms/GPL5175-3188.txt"), data.table=FALSE, na.strings=STRNA)
gp5715$ID <- as.character(gp5715$ID)
gp5715 <- subset(gp5715, ! is.na(gene_assignment))
gp5715 <- subset(gp5715, seqname %in% CHROMS)

gp5715$gene_assign_id <- gp5715$gene_assign_symbol <- NA
for(k in 1:nrow(gp5715)) {
  rk <- strsplit(gp5715$gene_assignment[k], " // ", fixed=TRUE)[[1]]
  gp5715$gene_assign_id[k] <- rk[1]
  gp5715$gene_assign_symbol[k] <- rk[2]
}
rm(k, rk)

## GSE65540 RNAseq: use ENSEMBL
ensembl <- getBM(attributes=BM_ATTR, mart=MART)
ensembl[ensembl==""] <- NA
ensembl$chromosome_name <- paste0("chr", ensembl$chromosome_name)

ensembl <- subset(ensembl, ! is.na(hgnc_symbol))
ensembl <- subset(ensembl, chromosome_name %in% CHROMS)

## GSE66921
gp13607 <- fread(paste0(DIR,"platforms/GPL13607-20416.txt"), data.table=FALSE, na.strings=STRNA)
gp13607$ID <- as.character(gp13607$ID)
gp13607$spotid_symbol <- gsub(".*:", "", gp13607$SPOT_ID)
gp13607 <- subset(gp13607, ControlType==0 & chr_coord!="unmapped" & !is.na(accessions))
gp13607$seqname <- gsub(":.*", "", gsub("hs|","",gp13607$chr_coord,fixed=TRUE))
gp13607 <- subset(gp13607, seqname %in% CHROMS)

## GSE199063
gp23126 <- fread(paste0(DIR,"platforms/GPL23126-131.txt"), data.table=FALSE, na.strings=STRNA)
gp23126$spotid_symbol <- gsub(" // .*", "", gp23126$SPOT_ID)
gp23126$spotid_category <- gsub(".* // ", "", gp23126$SPOT_ID)
gp23126 <- subset(gp23126, spotid_category=="RefSeq")
gp23126 <- subset(gp23126, seqname %in% CHROMS)

gp23126$gene_assign_id <- gp23126$gene_assign_symbol <- NA
for(k in 1:nrow(gp23126)) {
  rk <- strsplit(gp23126$gene_assignment[k], " // ", fixed=TRUE)[[1]]
  gp23126$gene_assign_id[k] <- rk[1]
  gp23126$gene_assign_symbol[k] <- rk[2]
}
rm(k, rk)

## Investigate common features in annotation files:
SYMSHARD <- intersect_all(ensembl$hgnc_symbol, gp13607$GeneName, gp23126$gene_assign_symbol, gp5715$gene_assign_symbol, gp7020$`Gene symbol`)
length(SYMSHARD) #12558

# ---------------------------- Expression Dataset Cleaning ----------------------------
## GSE29411:
gse29411 <- merge(
  fread(paste0(DIR,"raw_mat/GSE29409_series_matrix.txt"), data.table=FALSE),
  fread(paste0(DIR,"raw_mat/GSE29410_series_matrix.txt"), data.table=FALSE),
  by = "ID_REF"
)
gse29411 <- subset(gse29411, ID_REF %in% gp7020$ID)
gp7020 <- subset(gp7020, ID %in% gse29411$ID_REF)
gse29411 <- gse29411[match(gp7020$ID, gse29411$ID_REF), ]
stopifnot(identical(gse29411$ID_REF, gp7020$ID))
gse29411$ID_REF <- gp7020$`Gene symbol`

mean(duplicated(gse29411$ID_REF)) #0.1001395
gse29411 <- aggregate(. ~ ID_REF, data=gse29411, FUN=cbMean)

## GSE42715:
gse42715 <- fread(paste0(DIR,"raw_mat/GSE42715_series_matrix.txt"), data.table=FALSE)
gse42715$ID_REF <- as.character(gse42715$ID_REF)
gse42715 <- subset(gse42715, ID_REF %in% gp5715$ID)
gp5715 <- subset(gp5715, ID %in% gse42715$ID_REF)
gse42715 <- gse42715[match(gp5715$ID, gse42715$ID_REF), ]
stopifnot(identical(gse42715$ID_REF, gp5715$ID))
gse42715$ID_REF <- gp5715$gene_assign_symbol

mean(duplicated(gse42715$ID_REF)) #0.01750796
gse42715 <- aggregate(. ~ ID_REF, data=gse42715, FUN=cbMean)

## GSE65540 RNA-seq:
gse65540 <- fread(paste0(DIR,"raw_mat/GSE65540_SAT_processed_counts.txt.gz"), data.table=FALSE)
colnames(gse65540)[1] <- "ID_REF"
gse65540 <- subset(gse65540, `Ensembl Gene ID` %in% ensembl$ensembl_gene_id)
stopifnot(! anyDuplicated(gse65540$ID_REF)) #confirm no need to aggregate
gse65540$`Ensembl Gene ID` <- NULL

## GSE66921:
gse66921 <- fread(paste0(DIR,"raw_mat/GSE66921_series_matrix.txt"), data.table=FALSE)
gse66921$ID_REF <- as.character(gse66921$ID_REF)
gse66921 <- subset(gse66921, ID_REF %in% gp13607$ID)
gp13607 <- subset(gp13607, ID %in% gse66921$ID_REF)
stopifnot(identical(gp13607$ID, gse66921$ID_REF))
gse66921$ID_REF <- gp13607$GeneName

mean(duplicated(gse66921$ID_REF)) #0.4650961
gse66921 <- aggregate(. ~ ID_REF, data=gse66921, FUN=cbMean)

## GSE199063:
gse199063 <- fread(paste0(DIR,"raw_mat/GSE199063_series_matrix.txt"), header=TRUE, data.table=FALSE)
gse199063 <- subset(gse199063, ID_REF %in% gp23126$ID)
gp23126 <- subset(gp23126, ID %in% gse199063$ID_REF)

stopifnot(identical(gp23126$ID, gse199063$ID_REF)) #if fail, match
gse199063$ID_REF <- gp23126$gene_assign_symbol

mean(duplicated(gse199063$ID_REF)) #0.005534541
gse199063 <- aggregate(. ~ ID_REF, data=gse199063, FUN=cbMean)

# ---------------------------- Commonize & Process Expression Datasets ----------------------------
par(mfrow=c(2,5))

## Phase 1: Format matrices & impose variability filter:
gse29411 <- custom_prep(gse29411)
gse42715 <- custom_prep(gse42715) 
gse65540 <- custom_prep(gse65540)
gse66921 <- custom_prep(gse66921)
gse199063 <- custom_prep(gse199063) #avg count

## Ensure all sample names are GEO accessions:
meta <- read.csv(paste0(DIR,"calculated_profiles/samplesheet.csv"), stringsAsFactors=FALSE)
meta <- subset(meta, Dataset=="GSE65540")
stopifnot(all(colnames(gse65540) %in% meta$SampleID))
gse65540 <- gse65540[ , match(meta$SampleID, colnames(gse65540))]
stopifnot(identical(meta$SampleID, colnames(gse65540)))
colnames(gse65540) <- meta$Accession

## Phase 2: Commonize & Normalize on common set
commonSet <- intersect_all(rownames(gse199063), rownames(gse66921), rownames(gse42715), rownames(gse29411), rownames(gse65540))
length(commonSet)

gse29411 <- custom_proc(gse29411, commonSet)
gse42715 <- custom_proc(gse42715, commonSet)
gse65540 <- custom_proc(gse65540, commonSet)
gse66921 <- custom_proc(gse66921, commonSet, FALSE) #already normalized
gse199063 <- custom_proc(gse199063, commonSet)

## Phase 3: Join
expr <- NULL
expr <- custom_cbind(expr, gse29411)
expr <- custom_cbind(expr, gse42715)
expr <- custom_cbind(expr, gse65540)
expr <- custom_cbind(expr, gse66921)
expr <- custom_cbind(expr, gse199063)
dim(expr)

# write.csv(expr, paste0(DIR,"calculated_profiles/merged_expression.csv"), row.names=TRUE, quote=FALSE)
