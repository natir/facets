# heterozygous and keep flags of the SNPs
procSnps <- function(rcmat, ndepth=35, het.thresh=0.25, snp.nbhd=250, nX=23, unmatched=FALSE, ndepthmax=1000) {
    gc(verbose=TRUE, reset=TRUE, full=TRUE)
    cat("procSNPS: BEGIN\n")

    # keep only chromsomes 1-22 & X for humans and 1-19, X for mice
    # for other genomes (gbuild = udef) nX is number of autosomes plus 1
    chromlevels <- c(1:(nX-1),"X")
    # keep only snps with normal read depth between ndepth and 1000
    # reduce the data frame to these snps
    rcmat <- rcmat[(rcmat$Chromosome %in% chromlevels) & ((rcmat$NOR.DP >= ndepth) & (rcmat$NOR.DP < ndepthmax)),]
    
    gc(verbose=TRUE, reset=TRUE, full=TRUE)
    print("In procSNPS 1")
    
    # output data frame
    out <- list()
    out$chrom <- rcmat$Chromosome
    out$maploc <- rcmat$Position
    out$rCountT <- rcmat$TUM.DP
    out$rCountN <- rcmat$NOR.DP
    out$vafT <- 1 - rcmat$TUM.RD/rcmat$TUM.DP
    out$vafN <- 1 - rcmat$NOR.RD/rcmat$NOR.DP

    
    rm(rcmat)
    print("In procSNPS 2")
    gc(verbose=TRUE, reset=TRUE, full=TRUE)
    
    # make chromosome ordered and numeric
    out$chrom <- as.numeric(ordered(out$chrom, levels=chromlevels))
    print("OUT CHROM")
    gc(verbose=TRUE, reset=TRUE, full=TRUE)
    
    # call a snp heterozygous if min(vafN, 1-mafN) > het.thresh
    if (unmatched) {
        if (het.thresh == 0.25) het.thresh <- 0.1
        out$het <- 1*(pmin(out$vafT, 1-out$vafT) > het.thresh & out$rCountT >= 50)
    } else {
        out$het <- 1*(pmin(out$vafN, 1-out$vafN) > het.thresh)
    }
    print("OUT HET")
    gc(verbose=TRUE, reset=TRUE, full=TRUE)
    
    # scan maploc for snps that are close to one another (within snp.nbhd bases)
    # heep all the hets (should change if too close) and only one from a nbhd
    out <- as.data.frame(out)
    out <- out[scanSnp(out$maploc, out$het, snp.nbhd)==1 & out$rCountT>0,]
    
    print("OUT KEEP")
    gc(verbose=TRUE, reset=TRUE, full=TRUE)
    
    out
}

scanSnp <- function(maploc, het, nbhd) {
    n <- length(maploc)
    zzz <- .Fortran("scansnp",
                    as.integer(n),
                    as.double(maploc),
                    as.double(het),
                    keep=double(n),
                    as.double(nbhd))
    zzz$keep
}

# obtain logR and logOR from read counts and GC-correct logR
counts2logROR <- function(mat, gbuild, unmatched=FALSE, ugcpct=NULL, f=0.2) {
    gc(verbose=TRUE, reset=TRUE, full=TRUE)
    cat("counts2logROR: BEGIN\n")

    nchr <- max(mat$chrom) # IMPACT doesn't have X so only 22
    
    # gc percentage
    mat$gcpct <- rep(NA_real_, nrow(mat))

    gc(verbose=TRUE, reset=TRUE, full=TRUE)
    cat("counts2logROR: END GC allocation\n")

    # get GC percentages from pctGCdata package
    # loop thru chromosomes
    
    for (i in 1:nchr) {
        ii <- which(mat$chrom==i)
        # allow for chromosomes with no SNPs i.e. not targeted
        if (length(ii) > 0) {
            if (gbuild == "udef") {
                mat$gcpct[ii] <- getGCpct(i, mat$maploc[ii], gbuild, ugcpct)
            } else {
                mat$gcpct[ii] <- getGCpct(i, mat$maploc[ii], gbuild)
            }
        }
        rm(ii)
        gc(verbose=TRUE, reset=TRUE, full=TRUE)
        cat("\tcounts2logROR: in chr loop\n")
    }

    ##### log-ratio with gc correction and maf log-odds ratio steps
    chrom <- mat$chrom
    maploc_len <- length(mat$maploc)
    rCountN <- mat$rCountN
    rCountT <- mat$rCountT
    vafT <- mat$vafT
    vafN <- mat$vafN
    het <- mat$het
    gcpct <- mat$gcpct

    # compute gc bias
    ncount <- tapply(rCountN, gcpct, sum)
    tcount <- tapply(rCountT, gcpct, sum)
    pctgc <- as.numeric(names(ncount))
    tscl <- sum(ncount)/sum(tcount)
    gcb <- lowess(pctgc, log2(tcount*tscl)-log2(ncount), f=f)
    jj <- match(gcpct, gcb$x)
    gcbias <- gcb$y[jj]

    rm(ncount)
    gc(verbose=TRUE, reset=TRUE, full=TRUE)
    cat("counts2logROR: clean ncount\n")

    rm(tcount)
    gc(verbose=TRUE, reset=TRUE, full=TRUE)
    cat("counts2logROR: clean tcount\n")

    rm(pctgc)
    gc(verbose=TRUE, reset=TRUE, full=TRUE)
    cat("counts2logROR: clean pctgc\n")

    rm(gcb)
    gc(verbose=TRUE, reset=TRUE, full=TRUE)
    cat("counts2logROR: clean gcb\n")

    rm(jj)
    gc(verbose=TRUE, reset=TRUE, full=TRUE)
    cat("counts2logROR: clean jj\n")
    
    # compute cn log-ratio (gc corrected) and baf log odds-ratio
    cnlr <- log2(1+rCountT*tscl) - log2(1+rCountN) - gcbias

    rm(tscl)
    gc(verbose=TRUE, reset=TRUE, full=TRUE)
    cat("counts2logROR: clean tscl\n")
    
    # minor allele log-odds ratio and weights
    lorvar <- valor <- rep(NA_real_, maploc_len)
    
    if (unmatched) {
        # read count matrix for odds ratio etc
        rcmat <- round(cbind(vafT[het==1]*rCountT[het==1], (1-vafT[het==1])*rCountT[het==1]))
        # folded log of Tukey (with 1/6 correction)
        valor[het==1] <- log(rcmat[,1]+1/6) - log(rcmat[,2]+1/6)
        # variance - approximation using delta method
        lorvar[het==1] <- 1/(rcmat[,1]+1/6) + 1/(rcmat[,2]+1/6)
    } else {
        # read count matrix for odds ratio etc
        rcmat <- round(cbind(vafT[het==1]*rCountT[het==1], (1-vafT[het==1])*rCountT[het==1], vafN[het==1]*rCountN[het==1], (1-vafN[het==1])*rCountN[het==1]))
        # log-odds-ratio (Haldane correction)
        valor[het==1] <- log(rcmat[,1]+0.5) - log(rcmat[,2]+0.5) - log(rcmat[,3]+0.5) + log(rcmat[,4]+0.5)
        # variance of log-odds-ratio (Haldane; Gart & Zweifel Biometrika 1967)
        lorvar[het==1] <- (1/(rcmat[,1]+0.5) + 1/(rcmat[,2]+0.5) + 1/(rcmat[,3]+0.5) + 1/(rcmat[,4]+0.5))
    }
    # put them together
    mat$lorvar <- mat$valor <- mat$cnlr <- mat$gcbias <- rep(NA_real_, nrow(mat))
    mat$gcbias <- gcbias
    mat$cnlr <- cnlr
    mat$valor <- valor
    mat$lorvar <- lorvar

    gc(verbose=TRUE, reset=TRUE, full=TRUE)
    cat("END COUNTS2LOGROR\n")
    
    mat
}
