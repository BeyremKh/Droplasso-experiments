#dependecies 
library(parallel)
ncores <- detectCores() -1 
library(mvtnorm)
library(glmnet)
library(ROCR)
library(droplasso)


#generating function
generate_data <- function(n=100, d=20, d1=2, pi=1, w=10 , q=1) {
  # The samples 
  mu <- c(rep(0,d))
  Sigma <-  diag(d)
  rawvars <- rmvnorm(n, mean=mu, sigma=Sigma)
  pvars <- pnorm(rawvars)  
  z  <- qpois(pvars, pi)
  
  w <-c(rep(w,d1),rep(-w,d1),rep(0,d-2*d1))
  # The labels y
  y <- rbinom(n, 1, 1/(1+exp(-z %*%  w)) )
  # The corrupted samples x
  x <-sapply(1:d, function(i) z[,i] * rbinom(n,1,q))
  
  return(list(z=z, x=x, y=y, w=w))
}


#Comparaison function
compare <- function(nruns=10 , nlambda=10,   alphalist= seq(0,1,0.1), n_passes=1000, q=1) {  
  
  #preparing accuracies tables
  auc_tot_el=matrix(0,nrow=11,ncol=10)
  auc_tot_dl=matrix(0,nrow=11,ncol=10)
  auc_tot_drop=matrix(0,nrow=11,ncol=1)
  
  
  tot_auc <- mclapply(1:nruns,function(iexp){
    data_train <- generate_data(q=q)
    data_valid <- generate_data(n=10000,q=q)
    data_test <- generate_data(n=10000,q=q)
    d=nrow(data_train)
    d1=length(which(abs(data_train$w)>0))
    # The values of alpha to interpolate between lasso (alpha=1) and ridge (alpha=0)
    firstalpha <- TRUE
    for (ind in seq(length(alphalist))) {
      # Train glmnet model
      m_glm <- glmnet(data_train$x, data_train$y, family="binomial", intercept=F, alpha=alphalist[ind], standardize = F, nlambda=nlambda)
      # Predict on the validation set
      ypred <- data_valid$x %*% m_glm$beta
      pred <- sapply(seq(ncol(ypred)), function(j) {prediction(ypred[,j], data_valid$y)})
      auc <- unlist(lapply(pred, function(p) {performance(p, "auc")@y.values[[1]]}))
      auc_tot_el[ind,] <-  c(auc,rep(0.5, nlambda-length(auc))) 
    }
    # Assess AUC of the best lambda on the test set for ElasticNet
    lambdas=m_glm$lambda
    bestindex <- which(auc_tot_el == max(auc_tot_el), arr.ind = TRUE)[1,]
    m_glm <- glmnet(data_train$x, data_train$y, family="binomial", intercept=F, alpha=alphalist[bestindex[1]], standardize = F, nlambda=10)
    cat(m_glm$beta[,max(bestindex[2],length(m_glm$lambda))])
    ypred <- data_test$x %*% m_glm$beta[,max(bestindex[2],length(m_glm$lambda))]
    auc <- performance(prediction(ypred, data_test$y), "auc")@y.values[[1]]
    auclist_el <- auc
    # Assess AUC of the best lambda on the test set for LASSO 
    bestindex <- which(auc_tot_el[length(alphalist),] == max(auc_tot_el[length(alphalist),]))
    
    #train droplasso
    # The values of p we want to test
    probalist <- 0.6^seq(0,10)
    # The number of lambda values we want to test to form the regularization path
    # The number of epoch of the optimization procedure
    
    for (ind in seq(length(probalist))) {
      # Train on the train set
      m <- droplasso(data_train$x, data_train$y, family="binomial", keep_prob=probalist[ind], nlambda=nlambda, n_passes = n_passes,decay=0.1)
      # Pick lambda with the best AUC on the validation set
      ypred <- predict(m, data_valid$x)
      pred <- sapply(seq(ncol(ypred)), function(j) {prediction(ypred[,j], data_valid$y)})
      auc <- unlist(lapply(pred, function(p) {performance(p, "auc")@y.values[[1]]}))
      auc_tot_dl[ind,] <-  auc 
      m_drop <- droplasso(data_train$x, data_train$y, family="binomial", keep_prob=probalist[ind], lambda=0 , n_passes = n_passes,decay=0.1)
      ypred <- predict(m_drop, data_valid$x)
      pred <- sapply(seq(ncol(ypred)), function(j) {prediction(ypred[,j], data_valid$y)})
      auc_tot_drop[ind] <- unlist(lapply(pred, function(p) {performance(p, "auc")@y.values[[1]]}))
    }
    # Assess AUC of the best lambda on the test set
    bestindex <- which(auc_tot_dl == max(auc_tot_dl), arr.ind = TRUE)[1,]
    m <- droplasso(data_train$x, data_train$y, family="binomial", keep_prob=probalist[bestindex[1]], nlambda=nlambda, n_passes = n_passes,decay=0.1)
    cat(m$beta[,bestindex[2]])
    ypred <- predict(m, data_test$x,s=m$lambda[max(bestindex[2])])
    auc <- performance(prediction(ypred, data_test$y), "auc")@y.values[[1]]
    auclist_dl <- auc
    bestdrop <- which.max(auc_tot_drop)
    m<- droplasso(data_train$x, data_train$y, family="binomial", keep_prob=probalist[bestdrop], lambda=0, n_passes = n_passes,decay=0.1)
    ypred <- predict(m, data_test$x)
    pred <- sapply(seq(ncol(ypred)), function(j) {prediction(ypred[,j], data_valid$y)})
    auc <- performance(prediction(ypred, data_test$y), "auc")@y.values[[1]]
    auclist_drop <- auc
    return(c(auclist_el,auclist_dl,auclist_drop))
  },mc.cores=ncores)
  return(sapply(tot_auc, unlist))
}



#processing results:  
proc_compare <- function(nruns=10,q=1 nlambda=10,   alphalist= seq(0,1,0.1), n_passes=1000) {
  tot_auc = compare(nruns=nruns,q=q,  nlambda=nlambda,   alphalist= alphalist, n_passes=n_passes)  
  auclist_dl=tot_auc[1,]
  auclist_el=tot_auc[2,]
  auclist_drop=tot_auc[3,]
  
  
  p_auc= t.test(auclist_el, auclist_dl,var.equal=(var.test(auclist_el, auclist_dl)$p.value>0.05),paired=T)$p.value
  
  mean_acc_el=mean(auclist_el)
  mean_acc_dl=mean(auclist_dl)
  mean_acc_drop=mean(auclist_drop)
  
  
  
  sd_acc_el=sd(auclist_el)
  sd_acc_dl=sd(auclist_dl)
  sd_acc_drop=sd(auclist_drop)
  
  print(paste("the best average AUC for Elasticnet is ",mean_acc_el," with a standard deviation of ",sd_acc_el))
  print(paste("the best average AUC for DropLasso is ",mean_acc_dl," with a standard deviation of ",sd_acc_dl))
  print(paste("the best average AUC for dropout is ", mean_acc_drop," with a standard deviation of ",sd_acc_drop))
  print(paste("the p-value between DropLasso and ElasticNet is ", p_auc))
  
  
  boxplot(cbind(auclist_dl,auclist_drop,auclist_el),ylab="AUC",xlab="lambda index", col=rainbow(3),main=paste("Test AUC, noise rate:",q))
}
