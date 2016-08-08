

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(Biobase))
suppressPackageStartupMessages(library(EnrichmentBrowser))
#suppressPackageStartupMessages(library(hgu95av2))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(EGSEA))
suppressPackageStartupMessages(library(survHD))
suppressPackageStartupMessages(library(survival))

#setwd("~/Desktop/presentation")
wkdir ="~/git_src/Makena/"
srcdir= file.path(wkdir, "src")
datadir = file.path(wkdir, "data")

source(file.path(srcdir,"sup.moa2.R"))

## Files To be Read

### Read in Normalized RNA Seq
se <- readRDS(file.path(datadir, "original_data", "se_breastMAE.rds"))
se <- se$RNAseqNorm  # Class matrix

###  clinical tables
subtypes <- readRDS(file.path(datadir, "original_data", "subtypes.rds"))
pdata <- readRDS(file.path(datadir, "original_data", "pheno.rds"))


### Genesets

msigdb_immune_filepath <- file.path("~/Dropbox/DFCI (1)/GeneSets/MsigDB/c7.all.v5.1.symbols.gmt")
genesigdb_filepath <- file.path("~/Dropbox/DFCI (1)/GeneSets/GeneSigDB/ALL_SIGSv4.gmt")
bindea_gs <- readRDS("~/Dropbox/DFCI (1)/GeneSets/BindeaDB/BindeaSigGS_GeneSymbol.rds")

# Parse gene sets
msigdb_immune_gs <- parse.genesets.from.GMT(msigdb_immune_filepath)
genesigdb_gs <- parse.genesets.from.GMT(genesigdb_filepath)


## CLinical annotation (Ignore)
makeAnnot<-function(ExprSet=ExprSet) {
  annot <- readRDS(file.path(datadir,"Makena_clinical_table.rds"))
  annot<-merge(pData(ExprSet), annot)


  clin<-read.table(file.path(datadir, "BRCA.clin.merged.picked.txt"), sep="\t", header = TRUE,as.is=TRUE)
  clin<-t(clin)
  colnames(clin) = clin[1,]
  clin<-clin[-1,]
  clin<-clin[,!colnames(clin)%in%"Composite Element REF"]
  rownames(clin)<-toupper(gsub("\\.", "-", rownames(clin)))
  annot<-cbind(annot,clin[annot$sample,])
  #annot[,grep("days", colnames(annot))][1:5,]
  annot$days_to_last_followup<-as.numeric(as.character(annot$days_to_last_followup))
  annot$days_to_death<-as.numeric(as.character(annot$days_to_death))
  annot$TIME<-annot$days_to_last_followup
  annot$TIME[which(is.na(annot$TIME))]<-annot$days_to_death[which(is.na(annot$TIME))]
  #cbind(annot$TIME,annot[,grep("days", colnames(annot))])[1:5,]
}


### Diagnosis
###############################
## Depreciated Code
#-----------------------#
# Generate diagnosis table
#specimen <- na.exclude(readr::read_tsv(file.path(datadir,"original_data", "brca_specimen.txt")))

#specimen <- as.data.frame(sapply(specimen, make.true.NA))
#specimen2 <- specimen

#for(i in 1:ncol(specimen)){
#  if(i %in% 7:14){specimen[,i] <- as.numeric(as.character(specimen[,i]))}
#  else{specimen[, i] <- as.character(specimen[,i])}
#}
#--------------------------#

# specimen <- read.table(file.path(datadir,"original_data", "brca_specimen.txt"), sep="\t", header=TRUE, as.is=TRUE, na.strings ="[Not Available]")
#
# data_idc_ilc <- readr::read_tsv(file.path(datadir,"original_data", "BRCA_IDC_ILC.txt"))
#
# stype <- subtypes$breast
# stype <- stype %>% dplyr::select(1:2) %>%
#   dplyr::mutate(bcr_sample_barcode = gsub(".", "-", substr(sample, 1, 16), fixed = TRUE),
#                 colname = gsub(".", "-", substr(sample, 1, 15), fixed = TRUE),
#                 sample = gsub(".", "-", substr(sample, 1, 12), fixed = TRUE)) %>%
#   dplyr::left_join(data_idc_ilc, by="sample")
# stype$Diagnosis[is.na(stype$Diagnosis)] <- "not_available"
#

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
ExprSet <- readRDS(file.path(datadir, "ExprSet.rds"))  # Expression Set

normSet<-  as.data.frame(voom(ExprSet)$E)

annot<-annot[colnames(normSet),]
identical(rownames(annot), colnames(normSet))

normSet<-makeEset(normSet, annot)

normSet$highTumorFac<-factor(normSet$percent_tumor_cells>70)

lmRes<-eBayes(lmFit((normSet), model.matrix(~Diagnosis, data=pData(ExprSet) ) ))


### Normalization
#before.norm <- exprs(ExprSet)
#boxplot(before.norm[1:4000,])
#ExprSet <- EnrichmentBrowser::normalize(ExprSet, norm.method="quantile")
#after.norm <- exprs(ExprSet)
#boxplot(after.norm[1:4000,])
#ExprSet <- readRDS(file.path(datadir,"after.norm.ExprSet.rds"))

### Differential expression
# Make Groups
#pData(ExprSet)$GROUP <- ifelse(ExprSet$Diagnosis == "ILC", 1,0) # ILC is coded as 1 and IDC is coded as 0
#table(pData(ExprSet)$GROUP)
# Differential expresssion using edgeR
#diffexp <- de.ana(expr = as.matrix(normSet+6), grp=pData(ExprSet)$GROUP, de.method = "edgeR")

## Starting point
#ExprSet <- readRDS("ExprSet.rds")
#diffexp <- readRDS(file.path(datadir,"BRCA_diffexp.rds"))

normSet <- normSet[,!is.na(normSet$highTumorFac)]
lmRes <- eBayes(lmFit(normSet, model.matrix(~ Diagnosis + cut(percent_tumor_cells,3) , data = pData(normSet) ) ))

summary(decideTests(lmRes, p.value=0.05,lfc=1))
lmResSig<-topTable(lmRes, coef="DiagnosisILC", p.value=0.05, lfc=1, n=20000)

paste("Number genes ")
nrow(lmResSig[lmResSig$t < 0,])
nrow(lmResSig[lmResSig$t > 0,])

require(made4)
heatplot(normSet[rownames(lmResSig),], classvec=normSet$Diagnosis)

require(ComplexHeatmap)
ha = HeatmapAnnotation(df = pData(normSet)[,c ("Diagnosis", "percent_tumor_cells")], points = anno_points(normSet$percent_tumor_cells))
Heatmap(exprs(normSet[rownames(lmResSig),]), top_annotation = ha, top_annotation_height = unit(2, "cm"),clustering_distance_columns = "pearson",clustering_distance_rows = "pearson" , show_column_names = FALSE,column_names_gp = gpar(fontsize = 8))


# Get EntrezIds of genes
require(org.Hs.eg.db)
entID<-mapIds(org.Hs.eg.db, keys=rownames(normSet), column="ENTREZID", keytype="SYMBOL")
entID<-entID[!(is.na(entID))]
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


## limma

#lmFit(v) -> lmRes
#contrasts.fit(lmRes, makeContrasts(ILC = DiagnosisILC-TumFacmed -TumFachigh, levels=v$design))
#topTable(eBayes(lmRes))

# Use default baseGSEA
bindea_egsea<-egsea(v, makeContrasts(ILC = DiagnosisILC-TumFacmed -TumFachigh, levels=v$design), symbolsMap=v$genes,gs.annots=gs.bindea,display.top = 100, sort.by="avg.rank",egsea.dir="./EGSEA_bindea-report")

topSets(bindea_egsea,names.only=FALSE)
require(GSVA)
bindea_gsva<-gsva(v$E, bindea_gsENTREZ, mx.diff=1)$es.obs

# Clinical
# downloaded from http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/BRCA/20160128/gdac.broadinstitute.org_BRCA.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz

annot- read.csv(file.path(datadir, file='updated_annot.csv'))


print("Number of events")
table(annot$vital_status, annot$Diagnosis)


# Survival in TCGA data, probably not enough

plotKM(y=Surv(annot$TIME/365*12,annot$vital_status==1), strata=factor(annot$Diagnosis))

survcomp::km.coxph.plot(Surv(TIME/365*12,vital_status==1)~factor(Diagnosis), data.s=annot,x.label="Time (months)", y.label="Probability of survival", main.title="",show.n.risk=TRUE, n.risk.step=10, n.risk.cex=.3)

# Code to test genesets
a<-annot[colnames(bindea_gsva),]
bindea_surv<-list()

bindea_surv$CoxPH_Models<-apply(bindea_gsva,1, function(x) coxph(Surv(annot[colnames(bindea_gsva),"TIME"]/365*12,as.character(annot[colnames(bindea_gsva),"vital_status"])==1)~x))

bindea_surv$CoxPH_Summary<-sapply(bindea_surv$CoxPH_Models, broom::tidy)
nsig= which(bindea_surv$CoxPH_Summary["p.value",]<0.05)
bindea_surv$CoxPH_Summary[,nsig]

sapply(nsig,function(x) plotKMStratifyBy("median",y =Surv(annot[colnames(bindea_gsva),"TIME"]/365*12,as.character(annot[colnames(bindea_gsva),"vital_status"])==1),linearriskscore=bindea_gsva[x,], main=rownames(bindea_gsva)[x]))

