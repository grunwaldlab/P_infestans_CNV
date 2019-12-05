
# This file processes the peaks from 
# ../proc_gvcf.
#
# The resulting output includes:
# pinf_freq_wins.csv
# pinf_cnts_wins.csv


wins <- read.table("../pinf_ref/pinf_ref_wins.csv", header=TRUE, sep=",")
rownames(wins) <- paste(wins$CHROM, wins$start, wins$end, sep = "_")
wins2 <- wins


myPeaks <- list.files("../proc_gvcf/peaks/", pattern = "peaks")
myPeaks[1:4]

myCnts <- list.files("../proc_gvcf/peaks/", pattern = "counts")

#i <- 1

for(i in 1:length(myPeaks)){
  peak <- read.csv(paste("../proc_gvcf/peaks/", myPeaks[i], sep = ""))
  cnt  <- read.csv(paste("../proc_gvcf/peaks/", myCnts[i], sep = ""))
  
  rownames(peak) <- paste(unlist(lapply(strsplit(rownames(peak), "_"), function(x){paste(x[1], x[2], sep="_")})), peak$START, peak$END, sep="_")
  rownames(cnt) <- paste(unlist(lapply(strsplit(rownames(cnt), "_"), function(x){paste(x[1], x[2], sep="_")})), cnt$START, cnt$END, sep="_")
  
  myThreshold <- 20
  
  is.na(peak[,3:4][cnt[,3:4] < myThreshold]) <- TRUE
  
  wins <- cbind(wins, NA)
  colnames(wins)[ncol(wins)] <- colnames(peak)[3]
  wins[ rownames(peak), ncol(wins) ] <- peak[,3]
  
  wins2 <- cbind(wins2, NA)
  colnames(wins2)[ncol(wins2)] <- colnames(cnt)[3]
  wins2[ rownames(cnt), ncol(wins2) ] <- cnt[,3]
  
  wins <- cbind(wins, NA)
  colnames(wins)[ncol(wins)] <- colnames(peak)[4]
  wins[ rownames(peak), ncol(wins) ] <- peak[,4]
  
  wins2 <- cbind(wins2, NA)
  colnames(wins2)[ncol(wins2)] <- colnames(cnt)[4]
  wins2[ rownames(cnt), ncol(wins2) ] <- cnt[,4]
}


write.csv(wins, file = "pinf_freq_wins.csv")
write.csv(wins2, file = "pinf_cnts_wins.csv")

