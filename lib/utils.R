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

#GWAS_TRAIT_FILE_PROCESSING

y <- fread("./input/y_matrix.txt", header=T,data.table=FALSE)
y = as.matrix(y)
Z <- fread("./input/Z_matrix.txt", header=T,data.table=FALSE)
Z = as.matrix(Z)


#' \code{Conduct quick GWAS using rrBLUP}
#' 
#' Performs genome-wide association analysis based on the mixed model (Yu et al. 2006).
#'
#' @param y 
#' @param Z 
#'
#' @return return GWAS results.
#'
#' @examples
#'
#' @export
qGWAS <- function(y, Z, ...) {
  Zt = t(Z)
  
  Ztd = data.frame(marker=rownames(Zt),chrom=as.integer(as.character(gsub("-.*", "",rownames(Zt)))),
                   pos=as.integer(as.character(gsub(".*-", "",rownames(Zt)))),Zt,check.names=FALSE)
  rownames(Ztd) = 1:nrow(Ztd)
  yd = data.frame(line = 1:nrow(y), y=y)
  scores <- GWAS(pheno = yd, geno = Ztd, ...)
  
  gwas <- data.frame(Chr = paste0("Chr", scores$chrom), Start = scores$pos, End = scores$pos, p = scores$V1)  ## qq from Ropt
  return(gwas)
}


