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

#' \code{Obtain the first several PCs}
#'
#'
#' @param Z The `n x m` genotype matrix of data. [matrix, NULL]
#' @param p The first n number of PCs. [num, p=3]
#'
#' @return return a `n x p` matrix.
#'
#' @examples
#' Z <- fread("input/Z_matrix.txt", header=T,data.table=FALSE)
#' Z = as.matrix(Z)
#' X0 <- getpca(Z, p=3)
#'
#' @export
getpca <- function(Z, p=3) {
  X0 = prcomp(Z)$x[,1:p]
  return(X0)
}

