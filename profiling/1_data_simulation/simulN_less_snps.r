nQTL=20000
# we know that we have 20 chromosomes with 1000 QTLs each + 1 id
#nQTLperM=20
#nQTLar=1000
#sigma2m=2
#h2m=0.75
#nM=500
#h2=0.20
#h2ar=0.04
#varalpha=0.1
#c2m=(h2-h2ar)/h2m

library(data.table)
library(dplyr)


data_sim <- function(nQTLperM=100, nQTLar=100, h2_as=0.5, seed=1234 ){
  
  print(paste(nQTLperM, nQTLar, h2_as, seed))
  # Values from Weishaar et al paper, RFI
  h2=0.47*(h2_as/0.47)
  h2ar=0.20*(h2_as/0.47)
  h2m=0.61*(h2_as/0.47)
  # more realistic values, also use 5000 putative QTLs (1 every 4 markers) and 500 QTL per omic  and for the "ar" trait
  nloci=20000
  inter=4
  nQTL=nloci/inter
  sigma2m=2
  nM=1200
  sigma2alpha=0.01
  c2m=(h2-h2ar)/h2m
  
  cat("c2m = ",c2m)
  if(c2m<0 | c2m>1){
    cat("impossible values\n")
  }
  
  # from here we obtain
  sigma2am=h2m*nM*sigma2m*sigma2alpha
  sigma2em=(1-h2m)*nM*sigma2m*sigma2alpha
  cat("sigma2am= ",sigma2am,"\n")
  cat("sigma2em= ",sigma2em,"\n")
  sigma2p=sigma2am+sigma2em
  cat("sigma2p=sigma2em+sigma2am= ",sigma2p,"\n")
  sigma2y=h2m*sigma2p/(h2-h2ar) 
  cat("sigma2y= ",sigma2y,"\n")
  sigma2ar=sigma2y*h2ar
  cat("sigma2ar= ",sigma2ar,"\n")
  sigma2epsilon=sigma2y*(1-c2m-h2ar)
  cat("sigma2epsilon= ",sigma2epsilon,"\n")
  
  
  
  
  
  #a=scan("kk.txt") # file with genotypes
  #geno=matrix(a,ncol=nloci+1,byrow=TRUE)
  a = fread("/common/jyanglab/zhikaiyang/projects/mediation/largedata/geno/gddtosilk_22ksnps_matrix_new.txt",header = T, data.table = F)
  geno = cbind(1:nrow(a), a[,3:ncol(a)])
  set.seed(seed)
  id_sample = sort(sample(2:ncol(geno),20000))
  pathname = paste("id_sample_seed_",seed,".txt",sep = "")
  fwrite(as.list(id_sample), pathname, sep = "\n", row.names = F, col.names = F)
  geno = geno[,c(1,id_sample)]
  colnames(geno)[1] = "id"
  geno = as.matrix(geno)
  geno[,2:ncol(geno)]=geno[,2:ncol(geno)]+1
  # reduce to matrix of QTL keeping every 4 locus
  pos=which(1:nloci %% inter ==0)
  # add 1st col to keep id
  pos=c(1,pos+1)
  geno=geno[,pos]
  
  # we compute scaling factors using all markers for simplicity
  # remmber there is the id in geno
  freq=apply(as.matrix(geno[,-1]),2,mean)/2
  sc = sqrt(nQTLperM*2*sum(freq*(1-freq))/nQTL)
  
  
  id=geno[,1]
  geno=as.matrix(geno[,-1]) # remove id
  #rm(a)
  qtl=matrix(NA,ncol=nM,nrow=nQTLperM)
  u  =matrix(NA,ncol=nM,nrow=nQTLperM)
  
  ##########################
  set.seed(seed) # so that same residuals are not sampled over and over
  
  alpha=rnorm(nM,0,sqrt(sigma2alpha))
  alpha_mat = matrix(0,ncol = 4, nrow = length(alpha))
  alpha_mat[,1]=alpha
  alpha_mat[,2]=abs(alpha)
  alpha_mat[,3]=1:length(alpha)
  alpha_mat = as.data.frame(alpha_mat)
  colnames(alpha_mat) = c("alpha", "abs", "id", "order")
  
  
  alpha_sorted = alpha_mat %>% arrange(desc(abs))
  alpha_sorted$order = 1:length(alpha)
  alpha_sorted = alpha_sorted %>% arrange(id)
  
  
  #aid = c(which(alpha < quantile(alpha, c(0.005,0.995))[1]),which(alpha > quantile(alpha, c(0.005,0.995))[2]))
  #anid = c(1:1200)[-aid]
  #alpha[anid] = 0
  
  
  for (t in 12:1) {
    topi = t
    anid = which(alpha_sorted$order > topi)
    alpha[anid] = 0
    
    for (i in 1:nM){
      qtl[,i]=sample(1:nQTL,nQTLperM) # perhaps with overlap across Ms but we assume overlap across Ms to be small
      u[,i]=rnorm(nQTLperM)
    }
    
    qtlar=sample(1:nQTL,nQTLar)
    uar=rnorm(nQTLar)
    
    
    pathname = paste("mediationh2/alpha_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".txt",sep = "")
    
    write(alpha,file=pathname,ncolumns=1)
    
    pathname = paste("mediationh2/qtl_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".txt",sep = "")
    
    write(qtl,file=pathname,ncolumns=nM)
    
    pathname = paste("mediationh2/u_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".txt",sep = "")
    write(u,file=pathname,ncolumns=nM)
    
    pathname = paste("mediationh2/qtlar_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".txt",sep = "")
    write(qtlar,file = pathname,ncolumns = 1)
    
    pathname = paste("mediationh2/uar_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".txt",sep = "")
    write(uar,file=pathname,ncolumns=1)
    
    
    ############################
    lastgen = 1
    
    
    # 1st layer
    G=matrix(NA,ncol=nM,nrow=nrow(geno)) # here I keep only the genetic part, what Ole calls G
    E=matrix(NA,ncol=nM,nrow=nrow(geno)) 
    for (i in 1:nM){
      W=geno[,qtl[,i]]
      G[,i]=W%*%as.matrix(u[,i])
      # rescale to h2m so that on output var(M)=sigma2m
      # scale constant is based on the 1st 1100 animals (generation 1) so scaling is constant
      # across all calls
      G[,i]=sqrt(sigma2m)*sqrt(h2m)*G[,i]/sc
      varG=var(G[,i])
      varE=sigma2m*(1-h2m)
      E[,i]=rnorm(nrow(geno),0,sd=sqrt(varE))
      varE=var(E[,i])
      cat("nM=",i,"varG= ",varG,"varE= ",varE,"varM= ",varG+varE,"h2m ",varG/(varG+varE),"\n")
    }
    
    
    
    pathname = paste("mediationh2/M_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".txt",sep = "")
    fwrite(as.data.frame(G+E),file=pathname, quote = F, row.names = F, col.names = F)
    
    tbv_m = as.vector(G%*%alpha)
    e_m   = as.vector(E%*%alpha)
    p = tbv_m + e_m
    cat("var(tbv_m) sampled ",var(tbv_m),"\n")
    cat("var(e_m) sampled ",var(e_m),"\n")
    cat("var(p) sampled ",var(p),"\n")
    
    
    
    cat("computing required numbers\n")
    s2p=var(p)
    cat("s2p= ",s2p,"\n")
    s2y=h2m*s2p/(h2-h2ar)
    cat("s2y= ",s2y,"\n")
    var_ar=h2ar*s2y
    cat("var_ar= ",var_ar,"\n")
    c2m=(h2-h2ar)/h2m
    cat("c2m= ",c2m,"\n")
    var_epsilon=(1-c2m-h2ar)*s2y
    cat("var_epsilon= ",var_epsilon,"\n")
    W=geno[,qtlar]
    ar=as.vector(W%*%as.matrix(uar))
    ar=sqrt(var_ar)*ar/sd(ar)
    
    epsilon=rnorm(nrow(geno),0,sqrt(var_epsilon))
    
    y=tbv_m+e_m+ar+epsilon
    tbv=tbv_m+ar
    res=e_m+epsilon
    
    cat("after simulation: \n")
    cat("vary= ",var(y),"\n")
    cat("c2m= ",var(tbv_m+e_m)/var(y),"\n")
    cat("h2= ",var(tbv_m+ar)/var(y),"\n")
    cat("h2ar= ",var(ar)/var(y),"\n")
    
    # these are total genetic values , part of them are TBV because we don't multiply here
    pathname = paste("mediationh2/tbv_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".txt",sep = "")
    write(t(cbind(tbv_m,ar,tbv)),file=pathname,ncolumns=3,append=(lastgen!=1))
    
    pathname = paste("mediationh2/components_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".txt",sep = "")
    write(t(cbind(tbv_m,e_m,p,ar,epsilon)),file=pathname,ncolumns=5,append=(lastgen!=1))
    
    pathname = paste("mediationh2/y_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".txt",sep = "")
    write(y,file=pathname,ncolumns=1,append=(lastgen!=1))
    
    
  }
  
  
}







data_sim <- function(nQTLperM=100, nQTLar=100, h2_as=0.75, seed=1234 ){
  
  print(paste(nQTLperM, nQTLar, h2_as, seed))
  # Values from Weishaar et al paper, RFI
  h2=h2_as
  h2ar=0.2
  h2m=0.75
  # more realistic values, also use 5000 putative QTLs (1 every 4 markers) and 500 QTL per omic  and for the "ar" trait
  nloci=20000
  inter=4
  nQTL=nloci/inter
  sigma2m=2
  nM=1200
  sigma2alpha=0.01
  c2m=(h2-h2ar)/h2m
  
  cat("c2m = ",c2m)
  if(c2m<0 | c2m>1){
    cat("impossible values\n")
  }
  
  # from here we obtain
  sigma2am=h2m*nM*sigma2m*sigma2alpha
  sigma2em=(1-h2m)*nM*sigma2m*sigma2alpha
  cat("sigma2am= ",sigma2am,"\n")
  cat("sigma2em= ",sigma2em,"\n")
  sigma2p=sigma2am+sigma2em
  cat("sigma2p=sigma2em+sigma2am= ",sigma2p,"\n")
  sigma2y=h2m*sigma2p/(h2-h2ar) 
  cat("sigma2y= ",sigma2y,"\n")
  sigma2ar=sigma2y*h2ar
  cat("sigma2ar= ",sigma2ar,"\n")
  sigma2epsilon=sigma2y*(1-c2m-h2ar)
  cat("sigma2epsilon= ",sigma2epsilon,"\n")
  
  
  
  
  
  #a=scan("kk.txt") # file with genotypes
  #geno=matrix(a,ncol=nloci+1,byrow=TRUE)
  a = fread("/common/jyanglab/zhikaiyang/projects/mediation/largedata/geno/gddtosilk_22ksnps_matrix_new.txt",header = T, data.table = F)
  geno = cbind(1:nrow(a), a[,3:ncol(a)])
  set.seed(seed)
  id_sample = sort(sample(2:ncol(geno),20000))
  pathname = paste("id_sample_seed_",seed,".txt",sep = "")
  fwrite(as.list(id_sample), pathname, sep = "\n", row.names = F, col.names = F)
  geno = geno[,c(1,id_sample)]
  colnames(geno)[1] = "id"
  geno = as.matrix(geno)
  geno[,2:ncol(geno)]=geno[,2:ncol(geno)]+1
  # reduce to matrix of QTL keeping every 4 locus
  pos=which(1:nloci %% inter ==0)
  # add 1st col to keep id
  pos=c(1,pos+1)
  geno=geno[,pos]
  
  # we compute scaling factors using all markers for simplicity
  # remmber there is the id in geno
  freq=apply(as.matrix(geno[,-1]),2,mean)/2
  sc = sqrt(nQTLperM*2*sum(freq*(1-freq))/nQTL)
  
  
  id=geno[,1]
  geno=as.matrix(geno[,-1]) # remove id
  #rm(a)
  qtl=matrix(NA,ncol=nM,nrow=nQTLperM)
  u  =matrix(NA,ncol=nM,nrow=nQTLperM)
  
  ##########################
  set.seed(seed) # so that same residuals are not sampled over and over
  
  alpha=rnorm(nM,0,sqrt(sigma2alpha))
  alpha_mat = matrix(0,ncol = 4, nrow = length(alpha))
  alpha_mat[,1]=alpha
  alpha_mat[,2]=abs(alpha)
  alpha_mat[,3]=1:length(alpha)
  alpha_mat = as.data.frame(alpha_mat)
  colnames(alpha_mat) = c("alpha", "abs", "id", "order")
  
  
  alpha_sorted = alpha_mat %>% arrange(desc(abs))
  alpha_sorted$order = 1:length(alpha)
  alpha_sorted = alpha_sorted %>% arrange(id)
  
  
  #aid = c(which(alpha < quantile(alpha, c(0.005,0.995))[1]),which(alpha > quantile(alpha, c(0.005,0.995))[2]))
  #anid = c(1:1200)[-aid]
  #alpha[anid] = 0
  
  
  for (t in 12:1) {
    topi = t
    anid = which(alpha_sorted$order > topi)
    alpha[anid] = 0
    
    for (i in 1:nM){
      qtl[,i]=sample(1:nQTL,nQTLperM) # perhaps with overlap across Ms but we assume overlap across Ms to be small
      u[,i]=rnorm(nQTLperM)
    }
    
    qtlar=sample(1:nQTL,nQTLar)
    uar=rnorm(nQTLar)
    
    
    pathname = paste("mediationh2/alpha_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".txt",sep = "")
    
    write(alpha,file=pathname,ncolumns=1)
    
    pathname = paste("mediationh2/qtl_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".txt",sep = "")
    
    write(qtl,file=pathname,ncolumns=nM)
    
    pathname = paste("mediationh2/u_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".txt",sep = "")
    write(u,file=pathname,ncolumns=nM)
    
    pathname = paste("mediationh2/qtlar_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".txt",sep = "")
    write(qtlar,file = pathname,ncolumns = 1)
    
    pathname = paste("mediationh2/uar_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".txt",sep = "")
    write(uar,file=pathname,ncolumns=1)
    
    
    ############################
    lastgen = 1
    
    
    # 1st layer
    G=matrix(NA,ncol=nM,nrow=nrow(geno)) # here I keep only the genetic part, what Ole calls G
    E=matrix(NA,ncol=nM,nrow=nrow(geno)) 
    for (i in 1:nM){
      W=geno[,qtl[,i]]
      G[,i]=W%*%as.matrix(u[,i])
      # rescale to h2m so that on output var(M)=sigma2m
      # scale constant is based on the 1st 1100 animals (generation 1) so scaling is constant
      # across all calls
      G[,i]=sqrt(sigma2m)*sqrt(h2m)*G[,i]/sc
      varG=var(G[,i])
      varE=sigma2m*(1-h2m)
      E[,i]=rnorm(nrow(geno),0,sd=sqrt(varE))
      varE=var(E[,i])
      cat("nM=",i,"varG= ",varG,"varE= ",varE,"varM= ",varG+varE,"h2m ",varG/(varG+varE),"\n")
    }
    
    
    
    pathname = paste("mediationh2/M_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".txt",sep = "")
    fwrite(as.data.frame(G+E),file=pathname, quote = F, row.names = F, col.names = F)
    
    tbv_m = as.vector(G%*%alpha)
    e_m   = as.vector(E%*%alpha)
    p = tbv_m + e_m
    cat("var(tbv_m) sampled ",var(tbv_m),"\n")
    cat("var(e_m) sampled ",var(e_m),"\n")
    cat("var(p) sampled ",var(p),"\n")
    
    
    
    cat("computing required numbers\n")
    s2p=var(p)
    cat("s2p= ",s2p,"\n")
    s2y=h2m*s2p/(h2-h2ar)
    cat("s2y= ",s2y,"\n")
    var_ar=h2ar*s2y
    cat("var_ar= ",var_ar,"\n")
    c2m=(h2-h2ar)/h2m
    cat("c2m= ",c2m,"\n")
    var_epsilon=(1-c2m-h2ar)*s2y
    cat("var_epsilon= ",var_epsilon,"\n")
    W=geno[,qtlar]
    ar=as.vector(W%*%as.matrix(uar))
    ar=sqrt(var_ar)*ar/sd(ar)
    
    epsilon=rnorm(nrow(geno),0,sqrt(var_epsilon))
    
    y=tbv_m+e_m+ar+epsilon
    tbv=tbv_m+ar
    res=e_m+epsilon
    
    cat("after simulation: \n")
    cat("vary= ",var(y),"\n")
    cat("c2m= ",var(tbv_m+e_m)/var(y),"\n")
    cat("h2= ",var(tbv_m+ar)/var(y),"\n")
    cat("h2ar= ",var(ar)/var(y),"\n")
    
    # these are total genetic values , part of them are TBV because we don't multiply here
    pathname = paste("mediationh2/tbv_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".txt",sep = "")
    write(t(cbind(tbv_m,ar,tbv)),file=pathname,ncolumns=3,append=(lastgen!=1))
    
    pathname = paste("mediationh2/components_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".txt",sep = "")
    write(t(cbind(tbv_m,e_m,p,ar,epsilon)),file=pathname,ncolumns=5,append=(lastgen!=1))
    
    pathname = paste("mediationh2/y_corn_s0_top_", topi, "_qtl_",nQTLperM,"_qtlar_",nQTLar,"_h2_",h2_as,"_seed_",seed,".txt",sep = "")
    write(y,file=pathname,ncolumns=1,append=(lastgen!=1))
    
    
  }
  
  
}









for (nQTLperM in c(10,50,100)) {
  for (nQTLar in c(10,50,100)) {
    for (h2_as in c(0.25,0.5,0.75)) {
      for (seed in 1231:1240) {
        data_sim(nQTLperM, nQTLar, h2_as, seed)
      }
    }
  }
}











for (nQTLperM in c(2,4,8,16,32,64,128)) {
  for (nQTLar in c(2,4,8,16,32,64,128)) {
    for (h2_as in c(0.25,0.5,0.75)) {
      for (seed in 1231:1240) {
        data_sim(nQTLperM, nQTLar, h2_as, seed)
      }
    }
  }
}


for (seed in 1231:1240) {
  a = fread("/common/jyanglab/zhikaiyang/projects/mediation/largedata/geno/gddtosilk_22ksnps_matrix_new.txt",header = T, data.table = F)
  geno = cbind(1:nrow(a), a[,3:ncol(a)])
  set.seed(seed)
  id_sample = sort(sample(2:ncol(geno),20000))
  pathname = paste("id_sample_seed_",seed,".txt",sep = "")
  fwrite(as.list(id_sample), pathname, sep = "\n", row.names = F, col.names = F)
  geno = geno[,c(1,id_sample)]
  colnames(geno)[1] = "id"
  geno = as.matrix(geno)
  geno[,2:ncol(geno)]=geno[,2:ncol(geno)]+1
  # reduce to matrix of QTL keeping every 4 locus
  pos=which(1:nloci %% inter ==0)
  # add 1st col to keep id
  pos=c(1,pos+1)
  geno=geno[,pos]
}

a = fread("/common/jyanglab/zhikaiyang/projects/mediation/largedata/geno/gddtosilk_22ksnps_matrix_new.txt",header = T, data.table = F)
geno = cbind(1:nrow(a), a[,3:ncol(a)])
set.seed(seed)
id_sample = sort(sample(2:ncol(geno),20000))
pathname = paste("id_sample_seed_",seed,".txt",sep = "")
fwrite(as.list(id_sample), pathname, sep = "\n", row.names = F, col.names = F)
geno = geno[,c(1,id_sample)]
colnames(geno)[1] = "id"
geno = as.matrix(geno)
geno[,2:ncol(geno)]=geno[,2:ncol(geno)]+1
# reduce to matrix of QTL keeping every 4 locus
pos=which(1:nloci %% inter ==0)
# add 1st col to keep id
pos=c(1,pos+1)
geno=geno[,pos]


