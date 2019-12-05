

# This file reads in gvcf data,
# Among the outputs are two files written to
# ./peaks/ that include the estimate of the peaks
# and the counts of heterozygous positions in each window.


vcf_files <- list.files("../gvcf/", pattern = "2n.g.vcf.gz$")

samples <- unlist(lapply(strsplit(vcf_files, ".", fixed = TRUE), function(x){x[1]}))


#clean <- TRUE
#
clean <- FALSE

if( clean == TRUE ){
  unlink( paste("hists/", samples[i], "*", sep="") )
  unlink( paste("log/", samples[i], "*", sep="") )
  unlink( paste("peaks/", samples[i], "*", sep="") )
  unlink( paste("quantiles/", samples[i], "*", sep="") )
}


library(vcfR)


#
for(i in 1:length(samples)){
  if( !file.exists( paste("./peaks/", samples[i], "_peaks.csv", sep="") ) ){
    
    # cat( paste("\nSample:", samples[i], "\n\n") )
    
    time1 <- Sys.time()
    
    logFile <- paste("./log/", samples[i], "_log.txt", sep = "")
    cat(c(), file = logFile)
    
    vcf <- read.vcfR(paste("../gvcf/", vcf_files[i], sep = ""))
    
    # The gVCF band records include a homozygous for the
    # reference allele genotype and zeros for subsequent data.
    # These records are therefore uninformative and should be omitted.
    vcf <- vcf[grep("END=", getINFO(vcf), invert = TRUE),]
    
    ad <- extract.gt(vcf, element = 'AD')
    allele1 <- masplit(ad, record = 1)
    allele2 <- masplit(ad, record = 2)
    allele3 <- masplit(ad, record = 3)

    quants <- quantile(allele1[,1], probs=seq(0,1, by=0.05), na.rm = TRUE)
    write.table(quants, paste("./quantiles/", samples[i], "_quantiles.csv", sep=""), sep=",", col.names = FALSE)

    # Filter on depth

    # Allele 1
    #
    vcf@gt[allele1[,1] < quants['10%'] & !is.na(allele1[,1]),2] <- NA
    #vcf@gt[allele1[,1] < quants['15%'] & !is.na(allele1[,1]),2] <- NA
    #
    vcf@gt[allele1[,1] > quants['90%'] & !is.na(allele1[,1]),2] <- NA
    #vcf@gt[allele1[,1] > quants['95%'] & !is.na(allele1[,1]),2] <- NA

    # Allele 2
    #
    vcf@gt[allele2[,1] < quants['10%'] & !is.na(allele2[,1]),2] <- NA
    #vcf@gt[allele2[,1] < quants['15%'] & !is.na(allele2[,1]),2] <- NA
    #
    vcf@gt[allele2[,1] > quants['90%'] & !is.na(allele2[,1]),2] <- NA
    #vcf@gt[allele2[,1] > quants['95%'] & !is.na(allele2[,1]),2] <- NA
    
    # Create histograms for alllele depth
    pdf(paste("./hists/", samples[i], "_allele_hist.pdf", sep=""))
    par(mfrow=c(3,1))
    hist(allele1, breaks = seq(0,max(allele1, na.rm = TRUE), by=1), col="#beaed4", xlim = c(0,100))
    abline(v=quants[c('10%', '95%')], col=2)

    hist(allele2[ allele2 > 0 ], breaks = seq(0,max(allele1, na.rm = TRUE), by=1), col="#beaed4", xlim = c(0,100))
    abline(v=quants[c('10%', '95%')], col=2)

    hist(allele3[ allele3 > 0 ], breaks = seq(0,max(allele1, na.rm = TRUE), by=1), col="#beaed4", xlim = c(0,100))
    abline(v=quants[c('10%', '95%')], col=2)

    title(main=samples[i], outer=TRUE, line=-1)

    dev.off()
    par(mfrow=c(1,1))

    ## Heterozygous positions
    gt <- extract.gt(vcf)
    hets <- is_het(gt)
    vcf <- vcf[hets,]
    
#
    for(j in 1:4921){
      
        scont <- vcf[getCHROM(vcf) == paste('Supercontig_1.', j, sep=""),]
        
        cat(paste('Supercontig_1.', j, sep=""), file = logFile, append = TRUE)
        cat("\n", file = logFile, append = TRUE)
        
        # Extract allele depths
        ad <- extract.gt(scont, element = 'AD')
        allele1 <- masplit(ad, record = 1)
        allele2 <- masplit(ad, record = 2)
        ad1 <- allele1 / (allele1 + allele2)
        ad2 <- allele2 / (allele1 + allele2)
        
        # Parameters
        #winsize <- 1e5
        #
        winsize <- 2e5
        #
        bin_width <- 0.02
        
        freq1 <- ad1/(ad1+ad2)
        freq2 <- ad2/(ad1+ad2)
        myPeaks1 <- freq_peak(freq1, getPOS(scont), winsize = winsize, bin_width = bin_width)
        myPeaks2 <- freq_peak(freq2, getPOS(scont), winsize = winsize, bin_width = bin_width, lhs = FALSE)
        
        out <- cbind(myPeaks1$wins[,1:2, drop = FALSE], myPeaks1$peaks, myPeaks2$peaks)
        cnts <- cbind(myPeaks1$wins[,1:2, drop = FALSE], myPeaks1$counts, myPeaks2$counts)

        colnames(out)[3] <- paste(colnames(out)[3], "_allele1", sep="")
        colnames(out)[4] <- paste(colnames(out)[4], "_allele2", sep="")
        colnames(cnts)[3] <- paste(colnames(cnts)[3], "_allele1", sep="")
        colnames(cnts)[4] <- paste(colnames(cnts)[4], "_allele2", sep="")
        
        if( nrow(out) > 0 ){
          rownames(out) <- paste("Supercontig_1.", j, "_", rownames(out), sep="")
          rownames(cnts) <- paste("Supercontig_1.", j, "_", rownames(cnts), sep="")
        }
        
        if(j == 1){
          write.table(out, file = paste("./peaks/", samples[i], "_peaks.csv", sep=""), sep=",", append = FALSE, col.names=TRUE)
          write.table(cnts, file = paste("./peaks/", samples[i], "_counts.csv", sep=""), sep=",", append = FALSE, col.names=TRUE)
        } else if( nrow(out) > 0 ){
          write.table(out, file = paste("./peaks/", samples[i], "_peaks.csv", sep=""), sep=",", append = TRUE, col.names=FALSE)
          write.table(cnts, file = paste("./peaks/", samples[i], "_counts.csv", sep=""), sep=",", append = TRUE, col.names=FALSE)
        }
    }
    
    time2 <- Sys.time()

    cat("Start time: ", file = logFile, append = TRUE)
    cat(time1, file = logFile, append = TRUE)
    cat('\n', file = logFile, append = TRUE)

    cat("End time: ", file = logFile, append = TRUE)
    cat(time2, file = logFile, append = TRUE)
    cat('\n', file = logFile, append = TRUE)

    cat("Elapsed time (s): ", file = logFile, append = TRUE)
    cat(time2 - time1, file = logFile, append = TRUE)
    cat('\n', file = logFile, append = TRUE)

    cat("Elapsed time: ", file = logFile, append = TRUE)
    cat(format(time2 - time1, digits= 3), file = logFile, append = TRUE)
    cat('\n', file = logFile, append = TRUE)
  }
#
}


