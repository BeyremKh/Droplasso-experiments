#load models 
load("~/models_conquer.Rdata")

for (k in (1:length(models))){
  a=list(sapply(1:length(models[[k]]$list_lasso) , function(i) strsplit(models[[k]]$list_lasso,"[.]")[[i]][1]))
  file.create(paste(paste(models[[k]]$dataset$displayname,"lasso",sep="_"),"txt",sep="."))
  fileConn<-file(paste(paste(models[[k]]$dataset$displayname,"lasso",sep="_"),"txt",sep="."))
  writeLines(unlist(lapply(a[[1]], paste, collapse=" ")),fileConn)
  close(fileConn)
  a=list(sapply(1:length(models[[k]]$list_el) , function(i) strsplit(models[[k]]$list_el,"[.]")[[i]][1]))
  file.create(paste(paste(models[[k]]$dataset$displayname,"el",sep="_"),"txt",sep="."))
  fileConn<-file(paste(paste(models[[k]]$dataset$displayname,"el",sep="_"),"txt",sep="."))
  writeLines(unlist(lapply(a[[1]], paste, collapse=" ")),fileConn)
  close(fileConn)
  a=list(sapply(1:length(models[[k]]$list_dl) , function(i) strsplit(models[[k]]$list_dl,"[.]")[[i]][1]))
  file.create(paste(paste(models[[k]]$dataset$displayname,"dl",sep="_"),"txt",sep="."))
  fileConn<-file(paste(paste(models[[k]]$dataset$displayname,"dl",sep="_"),"txt",sep="."))
  writeLines(unlist(lapply(a[[1]], paste, collapse=" ")),fileConn)
  close(fileConn)
  a=list(sapply(1:length(models[[k]]$list_all) , function(i) strsplit(models[[k]]$list_all,"[.]")[[i]][1]))
  file.create(paste(paste(models[[k]]$dataset$displayname,"all",sep="_"),"txt",sep="."))
  fileConn<-file(paste(paste(models[[k]]$dataset$displayname,"all",sep="_"),"txt",sep="."))
  writeLines(unlist(lapply(a[[1]], paste, collapse=" ")),fileConn)
  close(fileConn)    
}

legend("bottomleft", legend=leg, col=1:5, pch=c(4,15:18), bty="n")
meanfeat_norm=sapply(1:nrow(meanfeat),function(i) meanfeat[i,]/ result[[i]]$ngenes)
#sparsity vs accuracy plot
matplot(meanAUC, 1- t(meanfeat_norm), ylab="Model sparsity",pch = c(4,15:18),bty="l")
