###########################################################################

# egsea.R
# Author: Makena, Aedin
# Date: 05/08/16
# Description: Based off Makena_Enrichment.R >> Run GSA using MsigDB/ImmunoDB, 
############## GeneSigDB and BindeaDB on LumA IDC v ILC data. 
############## Includes quantiile normalization step

###########################################################################

#### Library ####

## dplyr, readr, Biobase, EnrichmentBrowser, hgu95av2, ComplexHeatmap, limma, sup.moa2.R (additional functions)

source(file.path(srcdir, "sup.moa2.R"))
suppressPackageStartupMessages((library(dplyr)))             #  flexible grammar of data manipulation
suppressPackageStartupMessages((library(readr)))             #  Read flat/tabular text files from disk
suppressPackageStartupMessages((library(Biobase)))           #  Functions that are needed by many other packages or which
                                                             #  replace R functions
suppressPackageStartupMessages((library(EnrichmentBrowser))) #  essential functionality for the enrichment analysis of 
                                                             #  gene expression data. The analysis combines the advantages of
                                                             #  set-based and network-based enrichment analysis in order to derive 
                                                             #  high-confidence gene sets and biological pathways that are differentially 
                                                             #  regulated in the expression data under investigation. 
                                                             #  Besides, the package facilitates the visualization and exploration of such sets and pathways.
suppressPackageStartupMessages((library(hgu95av2)))          #  Affymetrix Human Genome U95 Set annotation data (hgu95av2) assembled using data from public data repositories
suppressPackageStartupMessages((library(ComplexHeatmap)))    #  Visualize associations between different sources of data sets and reveal potential structures. 
suppressPackageStartupMessages((library(limma)))             #  Data analysis, linear models and differential expression for microarray data.

#### Files ####

## Directories
setwd("/Users/makybincos/Desktop/Git/DFCI/")
wkdir ="~/Desktop/Git/DFCI/"
srcdir= file.path(wkdir, "src")
datadir = file.path(wkdir, "data")

## Normalized RNA Seq
nRNAseq <- readRDS(file.path(datadir, "original_data", "se_breastMAE.rds")) 
#file.path joins text as ".."/".."/"..."
nRNAseq <- nRNAseq$RNAseqNorm #filter for normalized data

## Clinical tables
subtypes <- readRDS(file.path(datadir, "original_data", "subtypes.rds")) # data on subtypes
pdata <- readRDS(file.path(datadir, "original_data", "pheno.rds"))       # data on phenotypes

## Gene Sets
#File location
msigdb_immune_filepath <- file.path("~/Desktop/presentation/GeneSets/MsigDB/c7.all.v5.1.symbols.gmt")
genesigdb_filepath <- file.path("~/Desktop/presentation/GeneSets/GeneSigDB/ALL_SIGSv4.gmt")
bindea_gs <- readRDS("~/Desktop/presentation/GeneSets/BindeaDB/BindeaSigGS_GeneSymbol.rds")
#Parse Gene Sets
msigdb_immune_gs <- parse.genesets.from.GMT(msigdb_immune_filepath)
genesigdb_gs <- parse.genesets.from.GMT(genesigdb_filepath)

###########################################################################

#####Diagnosis##########

## Generate Expression Set

ExprSet <- readRDS(file.path(datadir,"ExprSet.rds")) # Loads Expression Set

normSet <- as.data.frame(voom(ExprSet)$E) #estimates the mean-variance relationship of the log-counts, generates a precision weight for each observation and enters these into the limma empirical Bayes analysis pipeline

annot <- readRDS(file.path(datadir, "Makena_clinical_table.rds"))
annot <- merge(pData(ExprSet), annot) 
rownames(annot) <- annot$colname
annot <- annot[colnames(normSet),]
identical(rownames(annot), colnames(normSet)) #Checks if rownames of annot and colnames of normSet are the same
normSet<-makeEset(normSet, annot) #makeEset in sup.moa2

normSet$highTumorFac<-factor(normSet$percent_tumor_cells>70)

lmRes<-eBayes(lmFit((normSet), model.matrix(~Diagnosis, data=pData(ExprSet))))

normSet<-normSet[,!is.na(normSet$highTumorFac)] #we don't want high tumor fac
lmRes <- eBayes(lmFit(normSet, model.matrix(~ Diagnosis + cut(percent_tumor_cells,3) , data = pData(normSet))))

summary(decideTests(lmRes, p.value=0.05,lfc=1))
lmResSig<-topTable(lmRes, coef="DiagnosisILC", p.value=0.05, lfc=1, n=20000)

paste("Number genes ") #what's this for?
nrow(lmResSig[lmResSig$t < 0,])
nrow(lmResSig[lmResSig$t > 0,])

require(made4)
made4::heatplot(normSet[rownames(lmResSig),], classvec=normSet$Diagnosis)

-require(ComplexHeatmap)
ha = HeatmapAnnotation(df = pData(normSet)[,c ("Diagnosis", "percent_tumor_cells")], points = anno_points(normSet$percent_tumor_cells))
Heatmap(exprs(normSet[rownames(lmResSig),]), top_annotation = ha, top_annotation_height = unit(2, "cm"),clustering_distance_columns = "pearson",clustering_distance_rows = "pearson" , show_column_names = FALSE,column_names_gp = gpar(fontsize = 5))


# Get EntrezIds of genes
require(org.Hs.eg.db)
entID<-mapIds(org.Hs.eg.db, keys=rownames(normSet), column="ENTREZID", keytype="SYMBOL")
entID<-entID[!(is.na(entID))] #deletes NAs
print(paste(length(entID), "mapped to EntrezID"))
eset<-ExprSet[names(entID),]
rownames(eset)<-as.character(entID)
eset$percent_tumor_cells =annot$percent_tumor_cells
eset<-eset[,complete.cases(eset$Diagnosis, eset$percent_tumor_cells)]

eset$TumFac<-factor(cut(eset$percent_tumor_cells,3), labels=c("low", "med", "high"))

m=model.matrix(~ Diagnosis + TumFac, pData(eset) )
v= voom(eset,design=m)

#################
# Map Bindea GS
##################
#1. Map Symbols to EntrezGene IDs
bindea_gsENTREZ<-lapply(bindea_gs, function(x) as.character(mapIds(org.Hs.eg.db, keys= x,  keytype="SYMBOL",column="ENTREZID")))

#2. Build Index.
gs.bindea<-buildCustomIdx(rownames(v$E), bindea_gsENTREZ, anno = NULL, label = "bindea_gs",name = "Bindea Gene Sets", species = "Human", min.size = 3)


#################
# GeneSIgDB
################
genesigdb_Entrez<-GSEABase::getGmt(file.path(datadir, "GENESIGDBv4Entrez.gmt"), geneIdType=EntrezIdentifier())



######
## Build Index
####
gs.annots = buildIdx(entrezIDs=rownames(v$E), species="human",
                     msigdb.gsets="none",
                     kegg.updated=FALSE, kegg.exclude = "all")
gs.annots= buildMSigDBIdx(rownames(v$E), geneSets = "all", species = "Homo sapiens", min.size = 5, rdata.dir = NULL)
gs.annots<-gs.annots[c("c6", "c7")]

v$genes= data.frame(cbind(FeatureID = as.character(entID), Symbols=names(entID)))


res2<-egsea(v, makeContrasts(ILC = DiagnosisILC-TumFacmed -TumFachigh, levels=v$design), symbolsMap=v$genes,gs.annots=gs.annots[[1]], baseGSEAs=c("globaltest"),display.top = 100, sort.by="avg.rank",egsea.dir="./EGSEA_MSigDB-report")


# Get per tumor gene set scores
# For this we will use gsva.
# I had hope these scores would be provided by egsea but I can't see any argument to permit this
gene_set_surivival_analysis<-function(){}

# # Requie Globaltest
# require(globaltest)
# p.value(gt(normSet, Y= PT, genesets=geneset, sampling=TRUE))
#
# ## GSEA
# #GSEA_immune <- sbea(method = "gsea", eset = diffexp, gs = msigdb_immune_gs, alpha = 0.05)
# #GSEA_genesig <- sbea(method = "gsea", eset = diffexp, gs = genesigdb_gs, alpha = 0.05)
# GSEA_bindea <- sbea(method = "gsea", eset = diffexp, gs = bindea_gs, alpha = 0.05)
# GSEA_immune <- readRDS("GSEA_immune.rds")
# GSEA_genesig <- readRDS("GSEA_genesig.rds")
# GSEA_bindea <- readRDS(file.path(datadir,"GSEA_bindea.rds"))
#
#
# #### gs ranking
# #gs.ranks.bindea <- gs.ranking(GSEA_bindea)
# gs.ranks.bindea <- readRDS("bindea_sig_df.rds")
# gs.ranks.genesig <- readRDS("genesig_sig_df.rds")
# gs.ranks.immune <- readRDS("immune_sig_df.rds")
#
#
# # Make Heatmaps
# #
