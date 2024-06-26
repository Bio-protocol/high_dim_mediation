

#' \code{Obtain the first several PCs}
#'
#'
#' @param Z The `n x m` genotype matrix of data. [matrix, NULL]
#' @param p The first n number of PCs. [num, p=3]
#' @param chrlen The input file name of chromosome length. [chr, ="input/Chromosome_v4.txt"] 
#' @param med_res The input file name of mediation results. [chr, ="output/res_fixed_bic_trait_V1.csv"] 
#' @param tiff_out The plot output file name. [chr, ="graphs/circos.tiff"] 
#'
#' @return return a `n x p` matrix.
#'
#' @examples
#' Z <- fread("input/Z_matrix.txt", header=T,data.table=FALSE)
#' Z = as.matrix(Z)
#' X0 <- getpca(Z, p=3)
#'
#' @export
circos_med <- function(gwas_res="output/gwas_results.csv",
                       med_res="output/mediators_fixed_bic_trait_V1.csv", 
                       dsnp_res="output/dsnps_fixed_bic_trait_V1.csv", 
                       isnp_res="output/isnps_fixed_bic_trait_V1.csv",
                       chrlen="input/Chromosome_v4.txt", 
                       gene_position= "input/gene_pos.csv",
                       out_tiff = "graphs/circos.tiff"
                       ){
  

  tiff(out_tiff, res=600, units = "mm", height = 120, width = 120)
  par(mar=c(0,0,1,0))
  ### Set up the chrs
  data <- read.table(chrlen,head=T,stringsAsFactors=FALSE,sep='\t') 
  p_chr <- vector()
  p_start <- vector()
  p_end <- vector()
  for (i in 1:nrow(data)) {
    p_chr <- c(p_chr, rep(data[i,1],50))
    t <- sample(c(data[i,2]:(data[i,3] - 4000000)), 50)
    p_start <- c(p_start, t ) 
    p_end <- c(p_end, (t + 4000000))
  }
  p_lines <- data.frame(chr= p_chr, start = p_start, end = p_end, p_5 = rep(5, 50*nrow(data)), p_10 = rep(10, 50*nrow(data)) )
  
  ####################################################################################################################################
  # visual GWAS results
  circos.genomicInitialize(data,plotType="NULL")
  #chr_col=colours()[c(12,41,46,52,60,79,125,414,429,190)]
  circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
    chr = CELL_META$sector.index
    xlim = CELL_META$xlim
    ylim = CELL_META$ylim
    circos.rect(xlim[1], 0, xlim[2], 1, border = NA,col = colours()[407])
    circos.text(mean(xlim), mean(ylim), chr, cex = 1, col = "black",font = 2,
                facing = "inside", niceFacing = TRUE)
  }, track.height = 0.1, bg.border = NA)
  
  
  bg.col <- rep(colours()[c(407,140)], 5)
  
  gwas <- read.csv(gwas_res)
  ch=gwas[,1]
  ch=gsub("Chr","",ch)
  col=ifelse(as.numeric(ch)%%2==1,"darkblue", "darkred")
  gwas$Chr = as.factor(gwas$Chr)
  circos.genomicTrackPlotRegion(gwas,panel.fun = function(region, value, ...){
    #circos.genomicLines(region, value, type = "h")
    circos.genomicPoints(region, value, pch = 16, cex = 0.25, ...)}
    ,bg.col =bg.col, bg.border = "white",track.height = 0.40
  ) 
  

    #dSNPs
    dsnps_fixed_bic <- fread(dsnp_res, header = T, data.table=FALSE)
    dsnps_fixed_bic_d = data.frame(Chr = paste0("Chr",gsub("-.*","", dsnps_fixed_bic$snp)), 
                                   Start = as.integer(gsub(".*-","", dsnps_fixed_bic$snp)), 
                                   End = as.integer(gsub(".*-","", dsnps_fixed_bic$snp)), 
                                   value = 10.5)
    message(sprintf("###>>> ploting [ %s ] dSNPs ...", nrow(dsnps_fixed_bic_d)))
    dsnps_fixed_bic_d$Chr = as.factor(dsnps_fixed_bic_d$Chr)
    
    dsnps_fixed_bic_d <- rbind(dsnps_fixed_bic_d[1,], dsnps_fixed_bic_d)
    dsnps_fixed_bic_d$value[1] = 1
    circos.genomicTrack(dsnps_fixed_bic_d, 
                        panel.fun = function(region, value, ...) {
                          circos.genomicLines(region, value, type = "h", lty=2, col="red")
                        }, track.index =2)
    
    
    circos.genomicTrack(p_lines, stack = TRUE, 
                        panel.fun = function(region, value, ...) {
                          i = getI(...)
                          circos.genomicLines(region, value, col = "#0000FF", ...)
                        }, track.index =2)
    
  
  ####################################################################################################################################

    mediators_fixed_bic <- fread(med_res, header = T , data.table=FALSE)
    # mediators_fixed_bic = subset(mediators_fixed_bic, padj <0.05)
    message(sprintf("###>>> ploting [ %s ] mediators ...", nrow(mediators_fixed_bic)))
    
    color = brewer.pal(n = nrow(mediators_fixed_bic), name = "Set1")
    
    ### gene position
    genepos <- read.csv(gene_position)
    mediators_fixed_bic_pos <- merge(genepos, mediators_fixed_bic, by="id")
    mediators_fixed_bic_pos$col <- color
    
    circos.track(ylim = c(0, 0.05), track.height = 0.05, bg.border = "black")
    
    for(i in 1:nrow(mediators_fixed_bic_pos))
    {
      circos.rect((mediators_fixed_bic_pos$start[i]-1000000),0,(mediators_fixed_bic_pos$end[i] +1000000),0.05,
                  sector.index=mediators_fixed_bic_pos$chr[i], col= mediators_fixed_bic_pos$col[i], 
                  border= mediators_fixed_bic_pos$col[i],track.index =3)
    }
    
    
    isnps_fixed_bic <- fread(isnp_res, data.table=FALSE)
    isnp <- isnps_fixed_bic[order(isnps_fixed_bic$medi, abs(isnps_fixed_bic$coef), decreasing=TRUE),]
    
    # plot top 3 iSNPs
    isnp3 <- data.frame()
    tb <- data.frame(table(isnp$medi))
    for(i in 1:nrow(tb)){
      sub <- subset(isnp, medi %in% tb$Var1[i])
      if(nrow(sub) > 3){
        sub3 <- sub[1:3, ]
      }else{
        sub3 <- sub
      }
      isnp3 <- rbind(isnp3, sub3)
    }
    
    mediators_isnps_fixed_bic = data.frame(medi = isnp3$medi, isnps_for_medi = isnp3$snps_for_medi, 
                                           snp_chr = paste0("Chr", as.integer(gsub("-.*","",isnp3$snps_for_medi))), 
                                           snp_pos = as.integer(gsub(".*-","",isnp3$snps_for_medi)))
    #mediators_fixed_bic_pos_n = cbind(medi= mediators_fixed_bic$id, mediators_fixed_bic_pos[,c(1,2,3,5)])
    mediators_isnps_fixed_bic = merge(mediators_isnps_fixed_bic, mediators_fixed_bic_pos, by.x = "medi", by.y="id")
    
    for (i in 1 : nrow(mediators_isnps_fixed_bic)) {
      circos.link(mediators_isnps_fixed_bic$chr[i], c(mediators_isnps_fixed_bic$start[i], mediators_isnps_fixed_bic$end[i]), 
                  mediators_isnps_fixed_bic$snp_chr[i], mediators_isnps_fixed_bic$snp_pos[i], 
                  col = mediators_isnps_fixed_bic$col[i],  border = mediators_isnps_fixed_bic$col[i])
    }
  dev.off()
}

