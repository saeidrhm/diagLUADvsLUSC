invisible(lapply(c("scales","intervals","modeest", "bioDist", "Hmisc", "party", "nnet", "caret", "klaR",  "e1071", "rpart", "data.table", "abind", "plyr", "raster", "gplots", "ggplot2", "pheatmap", "reshape", "parallel", "zoo", "directlabels",  "Biobase", "GEOquery", "limma"), require, character.only = T))
options(mc.cores=4)
ai.algs <- c(nnet = function(...) nnet(...,size=20,maxit=300,trace=F), rpart=rpart, svm = svm, naiveBayes=naiveBayes)
ai.algs.pred <- c(rep("class",3), rep("none", length(ai.algs)-4))
factorize.data.frame <- function(x) do.call(cbind.data.frame, lapply(x, factor))

echo <- function(quiet, ...) if (! quiet) cat(paste(...))

Classify.init <- function(bm, target.col, cols=colnames(bm), class.num=2, cut.quantiles=F, alg=ai.algs, name=names(alg), pred.type=ai.algs.pred, print.table = F, return.model=F, parallel=T, noise.sd = 0, quiet=F) {
  bm <- data.frame(bm,check.names = F)
  echo(quiet, "Starting classification\n")
  cols <- unique(cols)
  if (any(! cols %in% colnames(bm))) {
    cat("Following columns do not exist: ", setdiff(cols, colnames(bm)))
    return(NULL)
  }
  cols <- colnames(bm[,unique(cols), drop=F])
  if (! target.col %in% cols) cols <- c(cols, target.col)
  bm <- na.omit(bm[,cols])
  cols.test <- cols[-which(cols==target.col)]
  if (class(bm[,target.col]) != "factor")
  {
    if (class(bm[,target.col])=="character") bm[,target.col]=factor(bm[,target.col]) else {
      print("Target column is not of type factor. Shaping it as factor\n")
      if (cut.quantiles) bm[,target.col] <- cut2(bm[,target.col], g=class.num) 
      else bm[,target.col] <- as.numeric(cut2(bm[,target.col], seq(min(bm[,target.col],na.rm=T)-1e-6, max(bm[,target.col],na.rm=T)+1e-6, length.out=class.num+1)))
      cat("Number of samples per classes:\n", table(bm[,target.col]), "\n") 
    }
  }
  if (class(alg) != "list") {
    alg <- list(alg)
    pred.type <- list(pred.type)
  } 
  if (! parallel) mclapply <- lapply
  list(bm=bm, target.col=target.col, cols=cols, cols.test = cols.test, alg=alg, pred.type=pred.type, name=name, mclapply=mclapply, print.table=print.table, noise.sd=noise.sd, return.model=return.model, quiet=quiet)
}

Classify.sample <- function(bm, target.col, trainSize=nrow(bm)/2, testSize=nrow(bm)/2, balance.factors=F) {
  if (balance.factors) prob <- 1/ddply(bm,target.col,nrow)[as.numeric(bm[,target.col]),2] 
  else prob <- NULL 
  all<- sample(nrow(bm), trainSize+testSize, prob=prob, replace=balance.factors || (trainSize + testSize) > nrow(bm))
  x <- sample(seq_along(all), trainSize)
  tr <- all[x]
  te <- all[-x]
  list(tr=tr, te=te)
}

Test.sample <- function(i, model, classify.init, classify.sample, test.data) {
  #     x <- switch(pred.type[[i]], none=predict(model[[i]], bm[te, cols.test,drop=F]), predict(model[[i]], bm[te, cols.test,drop=F], type=pred.type[[i]]))
  with(c(classify.init, classify.sample), {
    y <- test.data[,cols.test, drop=F]
    if (noise.sd > 0) y <- y + rnorm(nrow(y) * ncol(y), sd=noise.sd)
    x <- predict(model[[i]], y)
    if(class(x) == "list") x <- x[[1]]
    if (! is.null(dim(x))) { 
      x <- factor(apply(x,1,which.max))
      levels(x)=levels(bm[,target.col])
    }
    return(x)
  })
  
}

Classify.Accuracy <- function(classify.init, classify.sample, ...){
  with(c(classify.init, classify.sample), {
    model <- setNames(mclapply(seq(alg), function(i) alg[[i]](as.formula(paste0(target.col,"~.")),data=bm[tr,cols,drop=F], ...)), name)
    echo(quiet, "Model made\n")
    test.data <- na.omit(bm[te, c(cols.test, target.col),drop=F])
    test <- setNames(mclapply(seq(alg), Test.sample, model=model, classify.init=classify.init, classify.sample=classify.sample, test.data=test.data), name)
    if (print.table)
      lapply(seq(test), function(i) print(c(names(test)[i], list(table(Pred=test[[i]], True=test.data[,target.col])))))
    echo(quiet, "Finished classification\n")
    if (class(bm[,target.col]) == "factor") res <- sapply(test, function(x) mean(x==test.data[,target.col])) 
    else res <- sort(sapply(test, function(x) mean(abs(x - test.data[, target.col]))))
    #   hyb <- do.call(cbind, test)[,order(-res)]
    #   hyb <- factor(apply(hyb[,seq(hyb.algs)], 1, Mode))
    #   levels(hyb) <- levels(bm[,target.col])
    ##res <- list(accuracy=sort(res, decreasing=T))
    # if (return.test) res <- c(res, test = data.frame(pred=hyb, true=bm[te,target.col]))
    if (return.model) res <- c(res, model = model)
    return(res)
  })
}



ClassifyKFold <- function(bm, target.col, cols=colnames(bm), trainSize=nrow(bm)/2, testSize=nrow(bm)/2, alg=ai.algs, name=names(alg), pred.type=ai.algs.pred, print.table = F, balance.factors=F,quiet=F,hyb.algs=3,parallel=T,return.test=F,return.model=F, class.num=2, cut.quantiles=F, noise.sd = 0, Fold=4, ...) {
  ClassFoldRes <- lapply(seq(1,Fold),function(x)Classify(bm, target.col, cols, trainSize, testSize, alg, name, pred.type, print.table, balance.factors,quiet,hyb.algs,parallel,return.test,return.model, class.num, cut.quantiles, noise.sd))
  #browser()
  print(cols)
  print(ClassFoldRes)
  MeanRes <- mean(unlist(ClassFoldRes))
  print(target.col)
  print(MeanRes)
  return(MeanRes)
}

# bm: big matrix of data
# cols: all variables to use in building model
# target.col: the response
# trainSize: number of rows (samples) of bm to use in training model
# testSize: number of rows (samples) of bm to use while testing models
# algs: list of algorithms to use for classification, can be c(nnet = function(...) nnet(...,size=30,maxit=300,trace=F), knn3 = knn3, NaiveBayes = NaiveBayes, etc.
# balance.factors: balance the number of samples used for training and testing over different types of response - after factorization of response if needed
# quiet: don't print
# parallel: run in parallel
# return.model: return the model as a part of result
# class.num: number of classes used for factorization of the response
# cut.quantiles: if set, the response variable is factorized to class.num levels with almost the same number of samples in each level
#                otherwise the whole range of response is broken to class.num levels with equal sizes
Classify <- function(bm, target.col, cols=colnames(bm), trainSize=nrow(bm)/2, testSize=nrow(bm)/2, alg=ai.algs, name=names(alg), pred.type=ai.algs.pred, print.table = F, balance.factors=F,quiet=F,hyb.algs=3,parallel=T,return.test=F,return.model=F, class.num=2, cut.quantiles=F, noise.sd = 0, ...) {
  init.res <- Classify.init(bm, target.col, cols, class.num, cut.quantiles, alg, name, pred.type, print.table=print.table, return.model=return.model, parallel=parallel, noise.sd=noise.sd, quiet)
  bm <- init.res$bm
  samp.res <- Classify.sample(bm, target.col, trainSize, testSize, balance.factors)
  temp <- Classify.Accuracy(init.res, samp.res,...)
  #browser()
  return(temp)
}
ClassifyCombinations <- function(bm, target, cols = colnames(bm), comb.num=1, iterations=100, balance.factors=T, parallel=T, alg=c(svm=svm),...) {
  cols <- cols[cols != target]
  if (! parallel) mclapply <- lapply
  if (comb.num == 1) z <- lapply(seq(iterations),function(i) unlist(mclapply(cols, function(x) Classify(bm,target,as.character(x),alg=alg,quiet=T, parallel = F, balance.factors=balance.factors,...))))
  else z <- mclapply(seq(iterations),function(i) as.numeric(combn(cols,comb.num,function(x) Classify(bm,target,as.character(x),alg=alg,quiet=T, parallel = F, balance.factors=balance.factors,...))))
  z <- do.call(cbind, z)
  if (class(z) == "numeric") z <- matrix(z,ncol=iterations)
  z <- setNames(as.numeric(rowMeans(z)), apply(combn(cols,comb.num),2,paste,collapse=" "))
  sort(z, decreasing=T)
}

ClassifyCrossCombinations <- function(bm, target, cols1 = colnames(bm), cols2 = colnames(bm), comb.num1=1, comb.num2=1, iterations=100, balance.factors=T, parallel=T, alg=c(svm=svm), ...) {
  cols1 <- setdiff(cols1, target)
  cols2 <- setdiff(cols2, target)
  if (! parallel) mclapply <- lapply
  cols1.comb <- combn(cols1, comb.num1,simplify=F)
  cols2.comb <- combn(cols2, comb.num2,simplify=F)
  comb <- expand.grid(cols1.comb, cols2.comb)
  comb <- t(apply(comb, 1, function(x) sort(unlist(strsplit(unlist(x), " ")))))
  comb <- comb[!apply(comb, 1, function(x) any(duplicated(x))),,drop=F]
  comb <- data.frame(t(comb[! duplicated(comb),]))
  
  z <- lapply(seq(iterations),function(i) unlist(mclapply(comb, function(x) Classify(bm,target,as.character(x),alg=alg,quiet=T, parallel = F, balance.factors=balance.factors,...))))
  z <- do.call(cbind, z)
  if (class(z) == "numeric") z <- matrix(z,ncol=iterations)
  z <- setNames(as.numeric(rowMeans(z)), apply(comb,2,paste,collapse=" "))
  sort(z, decreasing=T)
}

SelectFeatures <- function(bm, target, cols = colnames(bm), max.comb.num=1, candidates.num=10, iterations=5, balance.factors=T, parallel=T, alg=c(svm=svm),...) {
  print(alg)
  print(max.comb.num)
  candidates <- ClassifyCombinations(bm, target, cols, parallel = parallel, balance.factors=balance.factors, iterations=iterations, ...)
  results <- list(candidates)
  if (max.comb.num > 1) {
    for (i in 2:max.comb.num) {
      cat("Round ", i-1, "Done.\n")
      print(head(candidates, candidates.num))
      candidates <- ClassifyCrossCombinations(bm, target, head(names(candidates), candidates.num), cols, 1, 1, iterations=iterations, balance.factors = balance.factors, parallel=parallel,alg=alg, ... )
      results <- c(results, list(candidates))
    }
  }
  setNames(results, paste("Set",seq(max.comb.num)))
}

library("miRBaseConverter")
library(plyr)
library(e1071)

exp_AD = readRDS("D:\\Drmowla\\LungCancer\\LUAD\\TCGA_LUAD_miRNA2.rds")
exp_SC = readRDS("D:\\Drmowla\\LungCancer\\LUSC\\TCGA_LUSC_miRNA.rds")

rownames(exp_AD) <- miRNA_AccessionToName(rownames(exp_AD),targetVersion = "v21")$TargetName
rownames(exp_SC) <- miRNA_AccessionToName(rownames(exp_SC),targetVersion = "v21")$TargetName

cntSCC <- exp_SC[,grepl("01A",colnames(exp_SC))|grepl("01B",colnames(exp_SC))|grepl("02A",colnames(exp_SC))|grepl("01C",colnames(exp_SC))]
cntADC <- exp_AD[,grepl("01A",colnames(exp_AD))|grepl("01B",colnames(exp_AD))]

intersectmiRs <- intersect(rownames(cntADC),rownames(cntSCC))
mergedcnt <- cbind(as.matrix(cntADC[intersectmiRs,]), as.matrix(cntSCC[intersectmiRs,]))
gr <- factor(c(rep("ADC", ncol(cntADC)),rep("SCC", ncol(cntSCC))))
label <- as.data.frame(gr)
colnames(label)= 'label'
mergedcnt <- cbind(as.data.frame(t(mergedcnt)), label)

bm <- mergedcnt
target <- 'label'



SelectFeaturesRes <- (SelectFeatures(bm, target, candidates.num = 100, max.comb.num = 2, alg =c(nnet = function(...) nnet(...,size=20,maxit=300,trace=F)), parallel = T))
saveRDS(SelectFeaturesRes, file = "SelectFeaturesRes_nnet_3_v1.Rds")
data.copy <- readRDS(file = "SelectFeaturesRes_nnet_3_v1.Rds")
