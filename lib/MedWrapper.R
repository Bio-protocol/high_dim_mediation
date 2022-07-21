#' 
run_GMA <- function(y,X0,X,Z, ncores, model="MedFix_eq", output_folder="output/"){
  trait = colnames(y)[1]
  out_dsnp <- paste0(output_folder, "/dSNP_", model, "_trait_", trait, ".csv")
  out_isnp <- paste0(output_folder, "/iSNP_", model, "_trait_", trait, ".csv")
  out_med <- paste0(output_folder, "/mediator_", model, "_trait_", trait, ".csv")
  
  ###' perform the hypotheses testing step that control fasle discovery proportion (FDP) below gamma
  ###' P(FDP>gamma)<alpha where gamma = 0.1, and alpha could be 0.05
  p.adj.method = 'holmr'
  pval.cut=0.05
  message(sprintf("###>>> running [ %s ] model with [ %s ] cores ...", model, ncores))
  
  if(model == "MedFix_eq"){
    vs.fixed = adlasso2typeWrapper(y,X0,X,Z,pz=seq(0.01,0.99,length.out=20),pmax=length(y)-2,ncores=ncores) 
    
    mod.eq = vs.fixed[['mod.eq']] # extract the model that assign equal penalty on the two data types.
    # for mod.eq, run the mediator models for the mediators with non-zero estimated coefficient in the outcome model
    e2m.eq= e2mFixed(mod.eq,ncores=ncores) 
    tests.eq=testMedFix(mod.eq,e2m.eq,p.adj.method=p.adj.method)
    
    directsnps.eq = reportDirectSNPs(mod=mod.eq, tst=tests.eq)
    if(nrow(directsnps.eq) >=1){
      fwrite(directsnps.eq, out_dsnp, sep=",", row.names=FALSE, quote=FALSE)}
    
    mediators.eq = reportMediator(mod=mod.eq, tst=tests.eq)
    if(nrow(mediators.eq) >= 1){
      fwrite(mediators.eq, out_med, sep=",", row.names=FALSE, quote=FALSE)}
    
    indirect_snps.eq <- reportINDirectSNPs(mod=mod.eq, e2m=e2m.eq, tst=tests.eq, pval.cut=pval.cut, Z=Z, colnames_X=colnames(X))
    if(nrow(indirect_snps.eq)>=1){
      fwrite(indirect_snps.eq, out_isnp, sep=",", row.names=FALSE, quote=FALSE)}
  }
  if(model == "MedFix_fixed"){
    vs.fixed = adlasso2typeWrapper(y,X0,X,Z,ncores=ncores) 
    mod.fixed = vs.fixed[['model']]  # extract the model that minimizes BIC
    # for mod.fixed, run the mediator models for the mediators with non-zero estimated coefficient in the outcome model
    e2m.fixed= e2mFixed(mod.fixed,ncores=ncores) 
    
    tests.fixed=testMedFix(mod.fixed,e2m.fixed,p.adj.method=p.adj.method)
    directsnps.fixed = reportDirectSNPs(mod=mod.fixed, tst=tests.fixed)
    if(nrow(directsnps.fixed)>=1){fwrite(directsnps.fixed, out_dsnp, sep=",", row.names=FALSE, quote=FALSE)}
    
    mediators.fixed = reportMediator(mod=mod.fixed, tst=tests.fixed)
    if(nrow(mediators.fixed) >= 1){fwrite(mediators.fixed, out_med, sep=",", row.names=FALSE, quote=FALSE)}
    
    indirect_snps.fixed <- reportINDirectSNPs(mod=mod.fixed, e2m=e2m.fixed, tst=tests.fixed, pval.cut=pval.cut, Z=Z, colnames_X=colnames(X))
    if(nrow(indirect_snps.fixed)>=1){fwrite(indirect_snps.fixed, out_isnp, sep=",", row.names=FALSE, quote=FALSE)}
  }
  if(model == "MedMix_linear"){
    ###' run the estimation step for MedMix using linear kernel, and extract the model that minimizes BIC
    vs.mix.linear = adlassoMixWrapper(y,X0,X,Z,kernel='linear',ncores=ncores)
    mod.mix.linear = vs.mix.linear[['model']]
    tests.mix.linear = testMedMix(mod.mix.linear,p.adj.method=p.adj.method)
    ###' Report the mediators selected by each method
    mediators.mix.linear = reportMediator(mod=mod.mix.linear, tst=tests.mix.linear)
    if(nrow(mediators.mix.linear) >= 1){
      fwrite(mediators.mix.linear, out_med, sep=",", row.names=FALSE, quote=FALSE)}
  }
  if(model == "MedMix_shrink"){
    ###' run the estimation step for MedMix using shrink_EJ kernel, and extract the model that minimizes BIC
    vs.mix.shrink = adlassoMixWrapper(y,X0,X,Z,kernel='shrink_EJ',ncores=ncores)
    mod.mix.shrink = vs.mix.shrink[['model']]
    tests.mix.shrink = testMedMix(mod.mix.shrink,p.adj.method=p.adj.method)
    ###' Report the mediators selected by each method
    mediators.mix.shrink = reportMediator(mod=mod.mix.shrink, tst=tests.mix.shrink)
    if(nrow(mediators.mix.shrink) >= 1){
      fwrite(mediators.mix.shrink, out_med, sep=",", row.names=FALSE, quote=FALSE)}
  }
  message(sprintf("###>>> Analysis finished! Find results in folder [ %s ].", output_folder))
}



###' New wrapper functions for each method, they source the function script highmed2019.r
adlassoMixWrapper <- function(y,X0,X,Z,kernel='linear',ncores=1,pmax=length(y)-2){  
  # This function is a wrapper function that runs the adaptive lasso for linear mixed model 
  # for a sequence of penalty parameter lambda, and output the result that minimizes BIC
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
    try(adlassoMix(X,y,K,X0=X0,lambda=lambda,init=list(tau=var(y)/2,ve=var(y)/2),pmax = pmax,err.max=1e-5,iter.max=200,iter.lasso.max=1e4,method.var='REML')) ## 'MLE') ## method.var=MLE or REML
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
adlasso2typeWrapper <- function(y,X0,X,Z,pz=seq(0.01,0.99,length.out=20),pmax=length(y)-2,ncores=16){ 
  ## wrapper of medfix with tuning parameter selection based on BIC
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

e2mFixed <- function(mod,ncores=8){ 
  ## calculate the exposure to mediator effect using the fixed model
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
      list(coef=mod.lasso$b,cov=mod.lasso$cov.coef,dv=mod.lasso$dev,s0=mod.lasso$s0)
      #,bic=all.lasso$bic,df=all.lasso$df,lambda=all.lasso$lambda.seq)
    }
    stopCluster(cl)
  }else{
    Ball = NULL
  }
  return(Ball)
}

