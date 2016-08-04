setwd("~/Desktop/presentation")
source("sup.moa2.R")
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(EnrichmentBrowser))
suppressPackageStartupMessages(library(hgu95av2))
suppressPackageStartupMessages(library(ComplexHeatmap))



### Read in Normalized RNA Seq
se <- readRDS("se_breastMAE.rds")
se <- se$RNAseqNorm

###  clinical tables
subtypes <- readRDS("subtypes.rds")
pdata <- readRDS("pheno.rds")


### Diagnosis
# Generate diagnosis table
specimen <- na.exclude(readr::read_tsv("brca_specimen.txt"))
data_idc_ilc <- readr::read_tsv("BRCA_IDC_ILC.txt")
specimen <- as.data.frame(sapply(specimen, make.true.NA))
specimen2 <- specimen
for(i in 1:ncol(specimen)){
  if(i %in% 7:14){specimen[,i] <- as.numeric(as.character(specimen[,i]))}
  else{specimen[, i] <- as.character(specimen[,i])}
}

stype <- subtypes$breast
stype <- stype %>% dplyr::select(1:2) %>%
  dplyr::mutate(bcr_sample_barcode = gsub(".", "-", substr(sample, 1, 16), fixed = TRUE),
                colname = gsub(".", "-", substr(sample, 1, 15), fixed = TRUE),
                sample = gsub(".", "-", substr(sample, 1, 12), fixed = TRUE)) %>%
  dplyr::left_join(data_idc_ilc, by="sample")
stype$Diagnosis[is.na(stype$Diagnosis)] <- "not_available"


# Filter tables
#stype <- as.data.frame(stype[stype$subtype == "LumA" & stype$Diagnosis %in% c("ILC", "IDC") & stype$colname %in% intersect(colnames(se), stype$colname), ])
#rownames(stype) <- stype$colname
#diagnosis <- stype[, c("colname", "Diagnosis")]
#diagnosis$Diagnosis <- as.factor(diagnosis$Diagnosis)
#diagnosis$GROUP <- ifelse(diagnosis$Diagnosis == "ILC", 1,0)
#RNAseqMat <- se[, stype$colname]

# Fix Rownames with Gene Symbol
#RNAseqMat <- RNAseqMat[16:nrow(RNAseqMat),] # Get rid of rows without gene symbols
#x <- rownames(RNAseqMat)
#id_x <- cbind(orgid=x,reshape2:::colsplit(x ,"\\|", c("symbol", "entrezid")))
#RNAseqMat <- RNAseqMat[!duplicated(id_x$symbol), ]
#rownames(RNAseqMat) <- id_x[!duplicated(id_x$symbol), "symbol"]


### Generate expressionSet
#pdata <- AnnotatedDataFrame(data = stype)
#ExprSet <- new("ExpressionSet", exprs = RNAseqMat, phenoData = pdata, annotation = "hgu95av2")
ExprSet <- readRDS("ExprSet.rds")

### Normalization
#before.norm <- exprs(ExprSet)
#boxplot(before.norm[1:4000,])
#ExprSet <- EnrichmentBrowser::normalize(ExprSet, norm.method="quantile")
#after.norm <- exprs(ExprSet)
#boxplot(after.norm[1:4000,])
ExprSet <- readRDS("after.norm.ExprSet.rds")

### Differential expression
# Make Groups
#pData(ExprSet)$GROUP <- ifelse(ExprSet$Diagnosis == "ILC", 1,0) # ILC is coded as 1 and IDC is coded as 0
#table(pData(ExprSet)$GROUP)
# Differential expresssion using edgeR
#diffexp <- de.ana(expr = ExprSet, de.method = "edgeR")

## Starting point
#ExprSet <- readRDS("ExprSet.rds")
#diffexp <- readRDS("BRCA_diffexp.rds")

### Get Genesets
# File path
msigdb_immune_filepath <- file.path("GeneSets/MsigDB/c7.all.v5.1.symbols.gmt")
genesigdb_filepath <- file.path("GeneSets/GeneSigDB/ALL_SIGSv4.gmt")
# Parse gene sets
msigdb_immune_gs <- parse.genesets.from.GMT(msigdb_immune_filepath)
genesigdb_gs <- parse.genesets.from.GMT(genesigdb_filepath)
bindea_gs <- readRDS("GeneSets/BindeaDB/BindeaSigGS_GeneSymbol.rds")

## GSEA
#GSEA_immune <- sbea(method = "gsea", eset = diffexp, gs = msigdb_immune_gs, alpha = 0.05)
#GSEA_genesig <- sbea(method = "gsea", eset = diffexp, gs = genesigdb_gs, alpha = 0.05)
#GSEA_bindea <- sbea(method = "gsea", eset = diffexp, gs = bindea_gs, alpha = 0.05)
GSEA_immune <- readRDS("GSEA_immune.rds")
GSEA_genesig <- readRDS("GSEA_genesig.rds")
GSEA_bindea <- readRDS("GSEA_bindea.rds")


#### gs ranking
#gs.ranks.bindea <- gs.ranking(GSEA_bindea)
gs.ranks.bindea <- readRDS("bindea_sig_df.rds")
gs.ranks.genesig <- readRDS("genesig_sig_df.rds")
gs.ranks.immune <- readRDS("immune_sig_df.rds")


# Make Heatmaps
# 