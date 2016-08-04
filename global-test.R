#Do our genesets predict histology (idc vs ilc) and not tumor purity?

#Load libraries
library(globaltest)
library(GlobalAncova)

#Tumor Purity Data
eset <- readRDS("Makena_clinical_table.rds")


runGT<-function(eset, Histology="Diagnosis", PT= "percent_tumor_cells", basal=
                  eSetClin82$MolecularSubtype=="Basal", geneset=probeID) {
 
   # eset is a ExpressionSet,
  # Histology, PT are varLabels of eset
  # basal is a logical vector of TRUE/FALSE
  # geneset is a list of probeIDs which are contained in featureNames(eset)

  res= vector(length=5)
  res[1]= p.value(globaltest(eset, Y= Histology, genesets=geneset, sampling=TRUE))
  res[2]= p.value(globaltest(eset, Y= PT, genesets=geneset, sampling=TRUE))
  res[3]= p.value(globaltest(eset, Y=Histology, adjust=PT, genesets=geneset, sampling=TRUE))
  res[4]= p.value(globaltest(eset, Y= PT, adjust=Histology, genesets=geneset, sampling=TRUE))
  #basal= pData
  res[5]= p.value(globaltest(eset[,basal], Y= Site, genesets=geneset, sampling=TRUE))
  res = signif(res, 5)
  res = cbind(test=c("Met_Site", "MolecularSubtype", "MetSite_wo_MolSubtype",
                     "MolecularSubtype_wo_Site", "BasalSet"), res)
  return(res)
}
