
#### Initialization ###
load("~/result_new.Rdata")
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
  list(
    name = "EMTAB2805",
    classid = "cell_cycle_stage",
    class = c("G1",
              "G2M"),
    displayname = "EMTAB2805, G1 vs G2M"
  ),
  list(
    name = "EMTAB2805",
    classid = "cell_cycle_stage",
    class = c("G1",
              "S"),
    displayname = "EMTAB2805, G1 vs S"
  ),
  list(
    name = "EMTAB2805",
    classid = "cell_cycle_stage",
    class = c("S",
              "G2M"),
    displayname = "EMTAB2805, S vs G2M"
  ),
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

models = list()

idata=1

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
  
  m_dl <- droplasso(x=x, y=y, family="binomial", keep_prob=dl_proba[round(result[[idata]]$bestp.mean[4])], nlambda=nlambda, n_passes = n_passes)
  m_el <- glmnet(x=x, y=y, family="binomial", alpha=alphalist[round(result[[idata]]$bestp.mean[6])], nlambda=nlambda,  intercept = F,standardize = F)
  m_lasso <- glmnet(x=x, y=y, family="binomial", alpha=c(1), nlambda=nlambda,  intercept = F,standardize = F)
  
  list_dl <- colnames(x)[which(abs(m_dl$beta[,round(result[[idata]]$bestp.mean[3])])>epsilon)]
  list_el <- colnames(x)[which(abs(m_el$beta[,round(result[[idata]]$bestp.mean[5])])>epsilon)]
  list_lasso <- colnames(x)[which(abs(m_lasso$beta[,round(result[[idata]]$bestp.mean[9])])>epsilon)]
  
  models[[idata]] = list(dataset=dataset,npos=npos, nneg=nneg, ngenes=ncol(x), dl=m_dl,el=m_el,lasso=m_lasso, list_all=colnames(x), list_dl=list_dl,list_el=list_el,list_lasso=list_lasso)
  idata = idata+1 }

save(models, file="models_conquer.Rdata")


