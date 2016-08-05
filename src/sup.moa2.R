sup.moa2 <- function (X, sup, nf = 2, ks.stat = FALSE, ks.B = 1000, ks.cores = NULL)
{
  N <- length(sup)
  fs <- X@fac.scr
  repr <- sapply(X@data, dim)[1, ]
  nn <- names(X@data)
  load <- split(X@loading, f = rep(nn, repr))
  load <- load[nn]
  w <- rowSums(sapply(sup, colSums))
  if (any(w == 0))
    stop("unrelated pathways involved")
  normsup <- sup
  GSCoordinate_sep <- mapply(SIMPLIFY = FALSE, function(load,
                                                        sup, A) {
    a <- t(sup * A) %*% as.matrix(load[, 1:nf, drop = FALSE])
    colnames(a) <- paste("PC", 1:nf, sep = "")
    return(a)
  }, load = load, sup = normsup, A = split(X@w.data, names(X@w.data))[nn])
  GSCoordinate_comb <- Reduce("+", GSCoordinate_sep)
  contribution <- lapply(GSCoordinate_sep, function(supcor,
                                                    score) {
    a <- lapply(1:nf, function(i) {
      r <- outer(supcor[, i], score[, i])
      colnames(r) <- rownames(score)
      return(r)
    })
    a[is.na(a)] <- 0
    names(a) <- paste("PC", 1:nf, sep = "")
    return(a)
  }, score = fs)
  contribution_dataset <- lapply(contribution, function(x) {
    Reduce("+", x)
  })
  contribution_pc <- lapply(1:nf, function(i, cont) {
    a <- lapply(cont, function(x) x[[i]])
    Reduce("+", a)
  }, cont = contribution)
  names(contribution_pc) <- paste("PC", 1:nf, sep = "")
  contribution_total <- Reduce("+", contribution_dataset)
  csup <- do.call("rbind", sup)
  if (!ks.stat) {
    pmat <- .signifGS(X = X, sup = csup, A = X@w.data, score = contribution_total,
                      nf = nf)
    attr(pmat, "method") <- "zscore"
  }
  else {
    if (is.null(ks.cores))
      ks.cores <- getOption("mc.cores", 2L)
    cat("running bootstrapping for p values of KS.stat ...\n")
    pmat <- .ks.pval(X, sup, ks.B = ks.B, A = X@w.data, nf = nf,
                     mc.cores = ks.cores)
    attr(pmat, "method") <- "KS.stat"
  }
  res <- new("moa.sup", sup = sup, coord.comb = GSCoordinate_comb,
             coord.sep = GSCoordinate_sep, score = contribution_total,
             score.data = contribution_dataset, score.pc = contribution_pc,
             score.sep = contribution, p.val = pmat)
  return(res)
}

.signifGS <- function(X, sup, A, score, nf) {

  # define function
  ff <- function(x, n, score, infinite=FALSE) {
    lx <- length(x)
    if (infinite)
      sf <- 1 else
        sf <- sqrt((lx-n)/(lx-1))
      sum_sd <- sf * sd(x)/sqrt(n) * n
      sum_mean <- mean(x) * n
      pp <- abs(pnorm(score, mean = sum_mean, sd = sum_sd))
      2 * Biobase::rowMin(cbind(pp, 1-pp))
  }
  # reconstuct matrix using nf PCs
  U <- as.matrix(X@loading[, 1:nf, drop=FALSE])
  D <- diag(sqrt(X@eig[1:nf]), nrow = nf)
  V <- as.matrix(X@eig.vec[, 1:nf, drop=FALSE])
  rec <- (U %*% D %*% t(V)) * A
  # the number of feature in each GS
  supn <- colSums(sup != 0)
  # calculate the P value
  pmat <- sapply(1:ncol(rec), function(i) ff(rec[, i], supn, score = score[, i]))
  colnames(pmat) <- colnames(score)
  rownames(pmat) <- rownames(score)
  return(pmat)
}

matchgs <- function(genes, geneSets) {
  mat <- sapply(geneSets, function(x) match(genes, x, nomatch = 0))
  mat[mat > 0] <- 1
  return(mat)
}

make.true.NA <- function(x) if(is.character(x)||is.factor(x)){
  is.na(x) <- x=="[Not Available]"; x} else {
    x}

coxReport<-function(coxph.model) {
  fit = summary(coxph.model)
  n=nrow(fit$conf.int)
  res=sapply(seq_along(1:n), function(x) c(hazard.ratio=fit$conf.int[x] ,coef=fit$coefficients[x,"coef"],se_coef=fit$coefficients[1,3] ,lower95=fit$conf.int[x,"lower .95"] ,upper95=fit$conf.int[x,"upper .95"], coef.pvalue= fit$coefficients[x,"Pr(>|z|)"] ,Score_logrank.p.value =fit$sctest[3],   Wald.p.value = fit$waldtest[3], LikelihoodRatio.p.value = fit$logtest[3]))
  colnames(res) = rownames(fit$coefficients)
  #res$Wald.p.value = fit$waldtest[3]
  #res$LikelihoodRatio.p.value = fit$logtest[3]
  return(res)
}

require(mogsa)
makemgsa2 <- function(mfa, sup)
{
  return(c(moa=mfa, sup=sup))
}

getmgsa2 <- function (mgsa, value)
{
  if (value %in% c("call", "moa", "sup"))
    r <- slot(mgsa, value)
  else if (value %in% c("eig", "tau", "partial.eig", "eig.vec",
                        "loading", "fac.scr", "partial.fs", "ctr.obs", "ctr.var",
                        "ctr.tab", "RV"))
    r <- slot(mgsa@moa, value)
  else if (value %in% c("data", "coord.sep", "coord.comb",
                        "score", "score.data", "score.pc", "score.sep", "p.val"))
    r <- slot(mgsa@sup, value)
  else stop("unknown value selected.")
  return(r)
}


makeEset<-function(eSet, annt){
  metadata <- data.frame(labelDescription = colnames(annt), row.names=colnames(annt))
  phenoData<-new("AnnotatedDataFrame", data=annt, varMetadata=metadata)
  # pData(eSet) = pData(phenoData)
  if (inherits(eSet, "data.frame")) eSet= as.matrix(eSet)
  if (inherits(eSet, "ExpressionSet")) eSet=exprs(eSet)
  data.eSet<-new("ExpressionSet", exprs=eSet, phenoData=phenoData)
  print(varLabels(data.eSet))
  return(data.eSet)
}

