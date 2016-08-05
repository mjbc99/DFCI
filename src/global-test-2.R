################################
# global-test-2.R
# Author: Makena, Aedin
# Date: Aug 4th 2016
# Run global test/ancova GSA on genesets that discriminate between BRCA TCGA cases of IDC and ILC but correct for tumor purity
###################################

#setwd("~/Desktop/presentation")
wkdir ="~/git_src/Makena/"
srcdir= file.path(wkdir, "src")
datadir = file.path(wkdir, "data")

#Do our genesets predict histology (idc vs ilc) and not tumor purity?

#Load libraries
library(globaltest)
library(GlobalAncova)

#Tumor Purity Data
annot <- readRDS(file.path(datadir,"Makena_clinical_table.rds"))


runGT<-function(eset, Histology="Diagnosis", PT= "percent_tumor_cells", lumA=
                  eset$subtype=="LumA",geneset=probeID) {

  # eset is a ExpressionSet,
  # Histology, PT are varLabels of eset
  # lumA is a logical vector of TRUE/FALSE <- Is this necessary?
  # geneset is a list of probeIDs which are contained in featureNames(eset)

  res= vector(length=5)
  res[1]= p.value(globaltest(eset, Y= Histology, genesets=geneset, sampling=TRUE))
  res[2]= p.value(globaltest(eset, Y= PT, genesets=geneset, sampling=TRUE))
  res[3]= p.value(globaltest(eset, Y=Histology, adjust=PT, genesets=geneset, sampling=TRUE))
  res[4]= p.value(globaltest(eset, Y= PT, adjust=Histology, genesets=geneset, sampling=TRUE))
  #lumA= pData
  res[5]= p.value(globaltest(eset[,lumA], Y= Histology, genesets=geneset, sampling=TRUE))
  res = signif(res, 5)
  res = cbind(test=c("Geneset_Diagnosis", "Percent_tumor", "Geneset_diagnosis_wo_Percent_tumor",
                     "MolecularSubtype_wo_Site", "LumA_Geneset"), res)
  return(res)
}
