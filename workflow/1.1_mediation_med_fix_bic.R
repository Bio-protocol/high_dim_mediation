setwd('/common/jyanglab/zhikaiyang/projects/high_dim_mediation')   #set working environment to the git repo



####################################### load required packages and functions #######################################

library(data.table)
library(glmnet)
library(MASS)
library(rrBLUP)
library(parallel)
library(doParallel)
library(CompQuadForm)
source('lib/highmed2019.r')
source('lib/fromSKAT.R')



###' New wrapper functions for each method, they source the function script highmed2019.r

adlassoMixWrapper <- function(Z,X,X0,y,kernel='linear',ncores=1,pmax=length(y)-2){  
  #This function is a wrapper function that runs the adaptive lasso for linear mixed model for a sequence of penalty parameter lambda, and output the result that minimizes BIC
  if(kernel=='shrink_EJ'){
    K = A.mat(Z,shrink=list(method="EJ"))
  }
  if(kernel=='linear'){
    K = A.mat(Z,shrink=FALSE)
  }
  n = length(y)
  eigK = eigen(K+sqrt(n)*diag(n))
  Qeig = eigK$vectors
  thetaeig = eigK$values-sqrt(n)
  Xt = t(Qeig)%*%cbind(cbind(1,X0),X)
  yt = t(Qeig)%*%y
  p = ncol(X)
  lambda.seq = getLambda(Xt,yt,nlambda=100,intercept=F,penalty.factor = c(rep(0,1+ncol(X0)),rep(1,p)))
  #    getDoParWorkers()
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  results.all = foreach(lambda=lambda.seq) %dopar% {
    library(glmnet)
    library(MASS)
    source('lib/highmed2019.r')
    try(adlassoMix(X,y,K,X0=X0,lambda=lambda,init=list(tau=var(y)/2,ve=var(y)/2),pmax = pmax,
                   err.max=1e-5,iter.max=200,iter.lasso.max=1e4,method.var='REML')) ## 'MLE') ## method.var=MLE or REML
  }
  names(results.all) = lambda.seq
  negll.all = sapply(results.all, getElement,name='negll')
  s0.all = sapply(results.all, getElement,name='s0')
  lambda.all = sapply(results.all, getElement,name='lambda')
  bic.all = sapply(results.all, getElement,name='bic')
  mod.final = results.all[[which.min(bic.all)]]
  return(list(model=mod.final,bic=bic.all,negll=negll.all,s0=s0.all,lambda=lambda.all,results.all=results.all)) 
}

###' 
###' 
adlasso2typeWrapper <- function(y,X0,X,Z,pz=seq(0.01,0.99,length.out=20),pmax=length(y)-2,ncores=16){ ## wrapper of medfix with tuning parameter selection based on BIC
  ## the fixed model should NOT include the intercept in X0
  n = length(y)
  p0 = ncol(X0)
  p = ncol(X)
  q = ncol(Z)
  getDoParWorkers()
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  out.all = foreach(ppzz=pz) %dopar% {
    library(glmnet)
    library(MASS)
    source('lib/highmed2019.r')
    adlasso.2type(y,X,Z,ppzz,X0,pmax=pmax)
  }
  stopCluster(cl)
  bic = sapply(out.all, getElement,name='bic.min')
  id.opt = which.min(bic)[1]
  out.final = out.all[[id.opt]]
  lambda.opt = out.final[['lambda.seq']][which.min(out.final[['bic']])[1]]
  mod.final = out.final[['model']]
  mod.final[['X0']] = X0
  mod.final[['X']] = X
  mod.final[['Z']] = Z
  mod.final[['y']] = y
  mod.eq= adlasso.2type(y,X,Z,pz=0.5,X0,pmax=pmax)[['model']]
  mod.eq[['X0']] = X0
  mod.eq[['X']] = X
  mod.eq[['Z']] = Z
  mod.eq[['y']] = y    
  return(list(model=mod.final,mod.eq=mod.eq,pz.opt=pz[id.opt],lambda.opt=lambda.opt,bic.prop=bic,bic.lambda=out.final[['bic']],pz.seq=pz,lambda.seq=out.final[['lambda.seq']]))
}

e2mFixed <- function(mod,ncores=8){ ## calculate the exposure to mediator effect using the fixed model
  ### mod = mod.fixed
  X = mod[['X']]
  n = nrow(X)
  X0 = mod[['X0']]#[,-1]
  if(is.null(X0)){
    p0=1
  }else{
    p0=ncol(X0)+1
  }
  Z = mod[['Z']]
  b = as.matrix(mod[['b']])
  ixnz = (b!=0)
  ns = sum(ixnz)
  #    print(ns)
  if(ns>0){
    bnz = b[ixnz]
    Xnz = matrix(X[,ixnz],nrow=n)
    #    print('good so far')
    getDoParWorkers()
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    Ball = foreach(ii=1:ns) %dopar% {
      library(glmnet)
      library(MASS)
      source('lib/highmed2019.r')
      ## ii = 1
      vs.lasso = adLasso(X=Z,y=Xnz[,ii],X0=X0)
      all.lasso = bic.glmnet(vs.lasso,cbind(1,cbind(X0,Z)),Xnz[,ii],p0=p0)
      mod.lasso = all.lasso$model
      list(coef=mod.lasso$b,cov=mod.lasso$cov.coef,dv=mod.lasso$dev,s0=mod.lasso$s0)#,bic=all.lasso$bic,df=all.lasso$df,lambda=all.lasso$lambda.seq)
    }
    stopCluster(cl)
  }else{
    Ball = NULL
  }
  return(Ball)
}






###' a function extracting the mediators making the cut
reportDirectSNPs <- function(mod=mod.fixed, tst=tests.fixed){
  coef.outcome = mod[['cov.coef']][['u.nz']]
  p.direct = tst$direct
  output <- data.frame()
  if(!is.null(coef.outcome)){
    output = data.frame(snp=names(coef.outcome), pval=p.direct, coef=coef.outcome)
  }
  return(output)
}

reportMediator <- function(mod=mod.fixed, tst=tests.fixed){
  coef.outcome = mod[['cov.coef']][['b.nz']]
  med <- data.frame()
  
  if(!is.null(coef.outcome)){
    med = tst[['med.individual']]
    # ix.output = (med$padj <= pval.cut)
    # output = list(pval=med[ix.output,],coef.outcome=coef.outcome[ix.output])
    med$coef <- coef.outcome
    med$id <- row.names(med)
  }
  return(med)
}

reportINDirectSNPs <- function(mod=mod.fixed, e2m=e2m.fixed, tst=tests.fixed, 
                               pval.cut=pval.cut, Z=Z, colnames_X=colnames(X)){
  med = tst[['med.individual']]
  ix.output = med$padj <= pval.cut
  idx <- c(1:length(ix.output))[ix.output]
  
  output <- data.frame()
  if(length(idx)>=1){
    e2m.sig <- e2m[idx]
    medid <- med$id[ix.output]
    colnames_X <- colnames_X
    mednames <- colnames_X[medid]
    Znames <- colnames(Z)
    for(i in 1:length(med$id[ix.output])){
      idx <- which(e2m.sig[[i]]$coef != 0)
      snpsnames <- Znames[idx]
      med_snps <- data.frame(medi=mednames[i], snps_for_medi=snpsnames, coef=e2m.sig[[i]]$coef[idx])
      output <- rbind(output, med_snps)
    }
  }
  return(output)
}


###' a function extracting the mediators making the cut
reportDirectSNPs <- function(mod=mod.fixed, tst=tests.fixed){
  coef.outcome = mod[['cov.coef']][['u.nz']]
  p.direct = tst$direct
  output <- data.frame()
  if(!is.null(coef.outcome)){
    output = data.frame(snp=names(coef.outcome), pval=p.direct, coef=coef.outcome)
  }
  return(output)
}

reportMediator <- function(mod=mod.fixed, tst=tests.fixed){
  coef.outcome = mod[['cov.coef']][['b.nz']]
  med <- data.frame()
  
  if(!is.null(coef.outcome)){
    med = tst[['med.individual']]
    # ix.output = (med$padj <= pval.cut)
    # output = list(pval=med[ix.output,],coef.outcome=coef.outcome[ix.output])
    med$coef <- coef.outcome
    med$id <- row.names(med)
  }
  return(med)
}

reportINDirectSNPs <- function(mod=mod.fixed, e2m=e2m.fixed, tst=tests.fixed, pval.cut=pval.cut, Z=Z, colnames_X=colnames(X)){
  med = tst[['med.individual']]
  ix.output = med$padj <= pval.cut
  idx <- c(1:length(ix.output))[ix.output]
  
  output <- data.frame()
  if(length(idx)>=1){
    e2m.sig <- e2m[idx]
    medid <- med$id[ix.output]
    colnames_X <- colnames_X
    mednames <- colnames_X[medid]
    Znames <- colnames(Z)
    for(i in 1:length(med$id[ix.output])){
      idx <- which(e2m.sig[[i]]$coef != 0)
      snpsnames <- Znames[idx]
      med_snps <- data.frame(medi=mednames[i], snps_for_medi=snpsnames, coef=e2m.sig[[i]]$coef[idx])
      output <- rbind(output, med_snps)
    }
  }
  return(output)
}










####################################### Load example  data and formating #######################################

y <- fread("input/y_matrix.txt", header=T,data.table=FALSE)
y = as.matrix(y)
Z <- fread("input/Z_matrix.txt", header=T,data.table=FALSE)
Z = as.matrix(Z)
#X0 = prcomp(Z)$x[,1:3] # use this line of code to calculate principal components if no X0_matrix.txt file provided
X0 <- fread("input/X0_matrix.txt", header=T,data.table=FALSE)
X0 = as.matrix(X0)
X <- fread("input/X_matrix.txt", header=T,data.table=FALSE)
X = as.matrix(X)

ncores = 4

trait = colnames(y)[1]


####################################### run mediation analysis #######################################

####################################### run the estimation step for Med_Fix_BIC 

vs.fixed = adlasso2typeWrapper(y,X0,X,Z,ncores=ncores) 

mod.fixed = vs.fixed[['model']]  # extract the model that minimizes BIC
e2m.fixed= e2mFixed(mod.fixed,ncores=ncores) # for mod.fixed, run the mediator models for the mediators with non-zero estimated coefficient in the outcome model

###' perform the hypotheses testing step that control fasle discovery proportion (FDP) below gamma
###' P(FDP>gamma)<alpha where gamma = 0.1, and alpha could be 0.05
p.adj.method = 'holmr' # I wrote an implementation of Holm test that controls FDP
pval.cut=0.05
tests.fixed=testMedFix(mod.fixed,e2m.fixed,p.adj.method=p.adj.method)
(res.fixed.pcut = medH.L2fixed(mod.fixed,e2m.fixed,tests.fixed,pval.cut=pval.cut))# report the results including the PVM calculation for MedFix_BIC
pathname = paste("output/res_fixed_bic_trait_", ".csv", sep = trait )
fwrite(as.list(res.fixed.pcut), pathname, sep=",", row.names=FALSE, quote=FALSE)



directsnps.fixed = reportDirectSNPs(mod=mod.fixed, tst=tests.fixed)
pathname = paste("output/dsnps_fixed_bic_trait_", ".csv", sep = trait)
if(nrow(directsnps.fixed)>=1){fwrite(directsnps.fixed, pathname , sep=",", row.names=FALSE, quote=FALSE)}


mediators.fixed = reportMediator(mod=mod.fixed, tst=tests.fixed)
pathname = paste("output/mediators_fixed_bic_trait_", ".csv", sep = trait)
if(nrow(mediators.fixed) >= 1){fwrite(mediators.fixed, pathname, sep=",", row.names=FALSE, quote=FALSE)}

indirect_snps.fixed <- reportINDirectSNPs(mod=mod.fixed, e2m=e2m.fixed, tst=tests.fixed, pval.cut=pval.cut, Z=Z, colnames_X=colnames(X))
pathname = paste("output/isnps_fixed_bic_trait_", ".csv", sep = trait)
if(nrow(indirect_snps.fixed)>=1){fwrite(indirect_snps.fixed, pathname, sep=",", row.names=FALSE, quote=FALSE)}



