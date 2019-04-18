# Test2 droplasso on single-cell RNA-seq from conquer
# Same as testjp.R, but perform correct parameter selection for elasticnet and dropout/droplasso
# JP Vert
# March 23, 2019

#### Initialization ###

# Load packages
library(MultiAssayExperiment)
library(parallel)
library(glmnet)
library(droplasso)
library(ROCR)
library(xtable)

# Random number initialization
set.seed(1234)

# Number of cores to use
numCores <- detectCores()-1

# Data to extract from conquer
conquerdata = list(
  # list(
  #   name = "EMTAB2805",
  #   classid = "cell_cycle_stage",
  #   class = c("G1",
  #             "G2M"),
  #   displayname = "EMTAB2805, G1 vs G2M"
  # ),
  # list(
  #   name = "EMTAB2805",
  #   classid = "cell_cycle_stage",
  #   class = c("G1",
  #             "S"),
  #   displayname = "EMTAB2805, G1 vs S"
  # ),
  # list(
  #   name = "EMTAB2805",
  #   classid = "cell_cycle_stage",
  #   class = c("S",
  #             "G2M"),
  #   displayname = "EMTAB2805, S vs G2M"
  # ),
  # list(
  #   name = "GSE45719",
  #   classid = "source_name_ch1",
  #   class = c("16-cell stage blastomere",
  #             "Mid blastocyst cell (92-94h post-fertilization)"),
  #   displayname = "GSE45719, 16-cell vs Mid blastocyst"
  # ),
  # list(
  #   name = "GSE45719",
  #   classid = "source_name_ch1",
  #   class = c("16-cell stage blastomere",
  #             "8-cell stage blastomere"),
  #   displayname = "GSE45719, 16-cell vs 8-cell"
  # ),
  list(
    name = "GSE48968-GPL13112",
    classid = "source_name_ch1",
    class = c("BMDC (1h LPS Stimulation)",
              "BMDC (4h LPS Stimulation)"),
    displayname = "GSE48968, 1h vs 4h"
  ),
  list(
    name = "GSE48968-GPL13112",
    classid = "source_name_ch1",
    class = c("BMDC (4h LPS Stimulation)",
              "BMDC (6h LPS Stimulation)"),
    displayname = "GSE48968, 4h vs 6h"
  ),
  list(
    name = "GSE74596",
    classid = "source_name_ch1",
    class = c("Single_cell_RNA-seq_NKT0",
              "Single_cell_RNA-seq_NKT17"),
    displayname = "GSE74596, NKT0 vs NKT17"
  ),
  list(
    name = "GSE74596",
    classid = "source_name_ch1",
    class = c("Single_cell_RNA-seq_NKT0",
              "Single_cell_RNA-seq_NKT1"),
    displayname = "GSE74596, NKT0 vs NKT1"
  ),
  list(
    name = "GSE74596",
    classid = "source_name_ch1",
    class = c("Single_cell_RNA-seq_NKT1",
              "Single_cell_RNA-seq_NKT2"),
    displayname = "GSE74596, NKT1 vs NKT2"
  )
)

# Fraction of cells that need to have nonzero counts to keep a gene
keepgenethreshold = 0.1

# Threshold to count non-zero weights in final models
epsilon <- 1e-8

# Number of repeats and folds for cross-validation
nrepeats <- 2
nfolds <- 5

# Number of lambda values to test (regularization parameter for glmnet and droplasso)
nlambda <- 10

# Grid of dropout probabilities (10 values)
dl_proba <- 1/(1+exp(seq(-5,4)))

# Grid of alpha parameters for elasticnet
alphalist <- seq(0,1,0.1)

# Number of passes during droplasso optimization
n_passes <- 5000

test_droplasso <- function(xtrain, ytrain, xvalid, yvalid, xtest, ytest, dl_proba, nlambda) {
  resauc <- list(valid = matrix(data=0, nrow=nlambda, ncol=length(dl_proba)) , test=matrix(data=0, nrow=nlambda, ncol=length(dl_proba)), nfeat=matrix(data=0, nrow=nlambda, ncol=length(dl_proba)))
  for (idlp in seq(length(dl_proba))) {
    dlp <- dl_proba[idlp]
    # train
    m <- droplasso(x=xtrain, y=ytrain, family="binomial", keep_prob=dlp, nlambda=nlambda, n_passes = n_passes)
    # Number of features
    resauc$nfeat[,idlp] <- apply(predict(m,s=m$lambda,type="coefficients"), 2, function(u) {sum(abs(u)>epsilon)})
    # prediction on valid
    yval <- predict(m, xvalid)
    for (ipar in seq(ncol(yval))) {
      pred <- prediction(yval[,ipar], yvalid)
      resauc$valid[ipar,idlp] <- performance(pred, "auc")@y.values[[1]]
    }
    #prediction on test
    yval <- predict(m, xtest)
    for (ipar in seq(ncol(yval))) {
      pred <- prediction(yval[,ipar], ytest)
      resauc$test[ipar,idlp] <- performance(pred, "auc")@y.values[[1]]
    }
  }
  return(resauc)
}



test_dropout <- function(xtrain, ytrain, xvalid, yvalid, xtest, ytest, dl_proba) {
  resauc <- list(valid=matrix(data=0, nrow=1, ncol=length(dl_proba)), test=matrix(data=0, nrow=1, ncol=length(dl_proba)), nfeat=matrix(data=0, nrow=1, ncol=length(dl_proba)))
  for (idlp in seq(length(dl_proba))) {
    dlp <- dl_proba[idlp]
    # train
    m <- droplasso(x=xtrain, y=ytrain, family="binomial", keep_prob=dlp, lambda=0, n_passes = n_passes)
    # Number of features
    resauc$nfeat[,idlp] <- apply(predict(m,s=m$lambda,type="coefficients"), 2, function(u) {sum(abs(u)>epsilon)})
    
    # prediction on valid and test sets
    yval <- predict(m, xvalid)
    for (ipar in seq(ncol(yval))) {
      pred <- prediction(yval[,ipar], yvalid)
      resauc$valid[ipar,idlp] <- performance(pred, "auc")@y.values[[1]]
    }
    
    yval <- predict(m, xtest)
    for (ipar in seq(ncol(yval))) {
      pred <- prediction(yval[,ipar], ytest)
      resauc$test[ipar,idlp] <- performance(pred, "auc")@y.values[[1]]
    }
  }
  return(resauc)
}


test_elasticnet <- function(xtrain, ytrain, xvalid, yvalid, xtest, ytest, alphalist, nlambda) {
  resauc <- list(valid=matrix(data=0, nrow=nlambda, ncol=length(alphalist)), test=matrix(data=0, nrow=nlambda, ncol=length(alphalist)), nfeat=matrix(data=0, nrow=nlambda, ncol=length(alphalist)))
  for (ialpha in seq(length(alphalist))) {
    myalpha <- alphalist[ialpha]
    # train
    m <- glmnet(x=xtrain, y=ytrain, family="binomial", alpha=myalpha, nlambda=nlambda,  intercept = F,standardize = F)
    # Number of features
    resauc$nfeat[,ialpha] <- apply(predict(m,s=m$lambda,type="coefficients"), 2, function(u) {sum(abs(u)>epsilon)})
    # prediction on valid and test sets
    yval <- predict(m, xvalid)
    for (ipar in seq(ncol(yval))) {
      pred <- prediction(yval[,ipar], yvalid)
      resauc$valid[ipar,ialpha] <- performance(pred, "auc")@y.values[[1]]
    }
    
    yval <- predict(m, xtest)
    for (ipar in seq(ncol(yval))) {
      pred <- prediction(yval[,ipar], ytest)
      resauc$test[ipar,ialpha] <- performance(pred, "auc")@y.values[[1]]
    }
  }
  return(resauc)
}



### Main loop to run experiments ###
result = list()
idata = 1

for (dataset in conquerdata) {
  print(
    paste(
      "Now working on dataset ",
      dataset$name,
      " with classes ",
      dataset$class[1],
      " and ",
      dataset$class[2],
      sep = ""
    )
  )
  
  # Read data from conquer
  print("Loading data")
  d = readRDS(gzcon(url(
    paste(
      "http://imlspenticton.uzh.ch/robinson_lab/conquer/data-mae/",
      dataset$name,
      ".rds",
      sep = ""
    )
  )))
  
  # Extract count matrix
  cts = assays(experiments(d)[["gene"]])[["count_lstpm"]]
  print(paste(ncol(cts), " cells, ", nrow(cts), " genes.", sep = ""))
  
  # Extract positive and negative cells
  pdata <- colData(d)
  pos = which(pdata[[dataset$classid]] == dataset$class[1])
  neg = which(pdata[[dataset$classid]] == dataset$class[2])
  npos = length(pos)
  nneg = length(neg)
  n = npos + nneg
  x <- t(cts)[c(pos, neg),]
  # Labels
  y = c(rep(0, npos), rep(1, nneg))
  # Remove genes with positive counts in enough cells
  nposcell <- apply(x, 2, function(i) {
    sum(i > 0)
  })
  x = x[, nposcell >= n * keepgenethreshold]
  # Take log(count+1)
  #x = log(x + 1)
  print(paste(
    npos,
    " positives, ",
    nneg,
    " negatives, ",
    ncol(x),
    " genes kept.",
    sep = ""
  ))
  
  # Make folds
  folds <- list()
  for (i in seq(nrepeats)) {
    folds <- c(folds, split(sample(seq(n)), rep(1:nfolds, length = n)))
  }
  
  # Main loop over methods
  resultloc <-
    parallel::mclapply(seq(length(folds)), function(ifold) {
      cat('.')
      itrain <- folds[[ifold]]
      ivalid <- folds[[ifold + 1 - nfolds*(ifold%%nfolds==0)]]
      itest <- seq(n)[-c(itrain,ivalid)]
      
      auc <- list()
      nfeat <- list()
      bestp <- list()
      # Dropout
      r <- test_dropout(xtrain=x[itrain,], ytrain=y[itrain], xvalid=x[ivalid,], yvalid=y[ivalid], xtest=x[itest,], ytest=y[itest], dl_proba = dl_proba)
      bestp[["dropout"]] <- which(r$valid == max(r$valid), arr.ind = T)[1,]
      auc[["dropout"]] <- r$test[bestp[["dropout"]][1],bestp[["dropout"]][2]]
      nfeat[["dropout"]] <- r$nfeat[bestp[["dropout"]][1],bestp[["dropout"]][2]]
      
      # Droplasso
      r <- test_droplasso(xtrain=x[itrain,], ytrain=y[itrain], xvalid=x[ivalid,], yvalid=y[ivalid], xtest=x[itest,], ytest=y[itest], dl_proba = dl_proba, nlambda=nlambda)
      bestp[["droplasso"]] <- which(r$valid == max(r$valid), arr.ind = T)[1,]
      auc[["droplasso"]] <- r$test[bestp[["droplasso"]][1],bestp[["droplasso"]][2]]
      nfeat[["droplasso"]] <- r$nfeat[bestp[["droplasso"]][1],bestp[["droplasso"]][2]]
      
      # Elasticnet
      r <- test_elasticnet(xtrain=x[itrain,], ytrain=y[itrain], xvalid=x[ivalid,], yvalid=y[ivalid], xtest=x[itest,], ytest=y[itest], alphalist = alphalist, nlambda=nlambda)
      bestp[["elasticnet"]] <- which(r$valid == max(r$valid), arr.ind = T)[1,]
      auc[["elasticnet"]] <- r$test[bestp[["elasticnet"]][1],bestp[["elasticnet"]][2]]
      nfeat[["elasticnet"]] <- r$nfeat[bestp[["elasticnet"]][1],bestp[["elasticnet"]][2]]
      
      # Ridge
      r <- test_elasticnet(xtrain=x[itrain,], ytrain=y[itrain], xvalid=x[ivalid,], yvalid=y[ivalid], xtest=x[itest,], ytest=y[itest], alphalist = c(0), nlambda=nlambda)
      bestp[["ridge"]] <- which(r$valid == max(r$valid), arr.ind = T)[1,]
      auc[["ridge"]] <- r$test[bestp[["ridge"]][1],bestp[["ridge"]][2]]
      nfeat[["ridge"]] <- r$nfeat[bestp[["ridge"]][1],bestp[["ridge"]][2]]
      
      # Lasso
      r <- test_elasticnet(xtrain=x[itrain,], ytrain=y[itrain], xvalid=x[ivalid,], yvalid=y[ivalid], xtest=x[itest,], ytest=y[itest], alphalist = c(1), nlambda=nlambda)
      bestp[["lasso"]] <- which(r$valid == max(r$valid), arr.ind = T)[1,]
      auc[["lasso"]] <- r$test[bestp[["lasso"]][1],bestp[["lasso"]][2]]
      nfeat[["lasso"]] <- r$nfeat[bestp[["lasso"]][1],bestp[["lasso"]][2]]
      
      return(list(auc=auc, nfeat=nfeat, bestp=bestp))
    }, mc.cores = numCores)
  
  # Compute mean AUC for each method
  auc <- sapply(resultloc, function(v) unlist(v[["auc"]]))
  print("Mean AUC")
  print(apply(auc, 1, mean))
  print("Standard deviation")
  print(apply(auc, 1, sd))
  # Compute P-values
  nmethods <- nrow(auc)
  pval <- matrix(NA,nrow=nmethods,ncol=nmethods,dimnames = list(rownames(auc),rownames(auc)))
  for (i in seq(nmethods)) {
    for (j in seq(nmethods)) {
      if (sd(auc[i,])+sd(auc[j,])>0) {
        pval[i,j] = t.test(auc[i,],auc[j,],alternative="greater")[["p.value"]]
      }
    }
  }
  # Number of selected features
  nfeat <- sapply(resultloc, function(v) unlist(v[["nfeat"]]))
  print("Mean Nfeat")
  print(apply(nfeat, 1, mean))
  print("Standard deviation")
  print(apply(nfeat, 1, sd))
  bestp <- sapply(resultloc, function(v) unlist(v[["bestp"]]))
  print("Mean bestparam")
  print(apply(bestp, 1, mean))
  print("bestparam Standard deviation")
  print(apply(bestp, 1, sd))
  result[[idata]] = list(dataset=dataset, npos=npos, nneg=nneg, ngenes=ncol(x), auc=auc, auc.mean=apply(auc, 1, mean), auc.sd=apply(auc, 1, sd), pval=pval, nfeat=nfeat, nfeat.mean=apply(nfeat, 1, mean), nfeat.sd=apply(nfeat, 1, sd),bestp= bestp, bestp.mean=apply(bestp, 1, mean))
  idata = idata+1
  save(result, file="result_new2.Rdata")
}

### Print and plot results ###

datasetnames <- sapply(conquerdata, function(v) v$displayname)
# Print mean AUC table
meanAUC <- t(sapply(result, function(v) v[["auc.mean"]]))
rownames(meanAUC) <- datasetnames
print(xtable(meanAUC), file="meanAUC2.tex")
# Print mean nfeat table
meanfeat <- t(sapply(result, function(v) v[["nfeat.mean"]]))
rownames(meanfeat) <- datasetnames
print(xtable(meanfeat), file="meanNfeat2.tex")
# Print pvalue table
smallpval <- Reduce('+', lapply(result, function(v) {pp <- v[["pval"]]; pp[is.na(pp)]=0.5; return(pp<0.05)}))
print(xtable(smallpval, caption=paste("Number of experiments, out of a total of ",length(result),", where each method (in row) is significantly better than each other method (in column).",sep="")),file="pval2.tex")
# Plot how many times each method is significantly better than another
pdf("pvalwins2.pdf",width=6,height=4)
par(mar = c(4, 7, 0, 0) + 0.2)
barplot(apply(smallpval, 1, sum),horiz=T,las=2,xlab="Number of wins")
dev.off()