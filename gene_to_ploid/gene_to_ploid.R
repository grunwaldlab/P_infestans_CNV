
# This file reads in several data files including:
# ../peak_proc/pinf_freq_wins.csv
#
# Output includes:
# gene_ploidy.csv
#


# Annotations
#genes <- read.table("../../broad/phytophthora_infestans_t30-4_1_genome_summary_bjk1.csv", header=TRUE, sep=",")
#genes <- read.table("../../pinf/broad/phytophthora_infestans_t30-4_1_genome_summary_bjk.txt", header=TRUE, sep="\t")
genes <- read.table("../pinf_ref/phytophthora_infestans_t30-4_1_genome_summary_bjk.txt", header=TRUE, sep="\t")

genes$CHROMOSOME <- paste("Supercontig_1.", genes$CHROMOSOME, sep="")
genes <- genes[!duplicated(genes$LOCUS), ]

genes[1:4, 1:6]

# Peaks
peaks <- read.csv("../peak_proc/pinf_freq_wins.csv", row.names = 1)

peaks$CHROM <- as.character(peaks$CHROM)

peaks[1:4, 14:18]


# Initialize a result data structure.
gmat <- matrix(nrow=nrow(genes), ncol=length(grep("allele1", colnames(peaks))))
colnames(gmat) <- grep("allele1", colnames(peaks), value = TRUE)
rownames(gmat) <- genes$LOCUS

gmat[1:4, 1:6]



#
# Set to TRUE to actually process the data.
#
#proc <- TRUE
#
proc <- FALSE


##### ##### ##### ##### #####
# Don't run!!!

if( proc == TRUE){

  for(i in 1:nrow(genes)){
    myPeak <- peaks[genes$CHROMOSOME[i] == peaks$CHROM, , drop = FALSE]
    myPeak <- myPeak[genes$START[i] >= myPeak$start, , drop = FALSE]
    myPeak <- myPeak[genes$STOP[i] <= myPeak$end, , drop = FALSE]
  
    myPeak <- myPeak[,grep("allele1", colnames(myPeak))]
  
    gmat[i,] <- colMeans(myPeak, na.rm = FALSE)
  }

#
  write.table(gmat, file = "gene_peak.csv", sep = ",")
  
  # Bin to ploidy
#  critical <- 1/4 - (1/3-1/4)/2
#  critical <- c(critical, 1/4 + (1/3-1/4)/2)
#  critical <- c(critical, 1/2 - (2/3 - 1/2)/2)
#  critical <- c(critical, 1/2 + (2/3 - 1/2)/2)
#  critical <- c(critical, 3/4 - (1/3-1/4)/2)
#  critical <- c(critical, 3/4 + (1/3-1/4)/2)

  critical <- c( 9/40, 7/24, 5/12, 7/12, 17/24, 31/40)
  
  
  gmat[ gmat <= 1 & gmat > critical[6] ] <- 5
  gmat[ gmat <= critical[6] & gmat > critical[5] ] <- 4
  gmat[ gmat <= critical[5] & gmat > critical[4] ] <- 3
  gmat[ gmat <= critical[4] & gmat >= critical[3] ] <- 2
  gmat[ gmat < critical[3] & gmat >= critical[2] ] <- 3
  gmat[ gmat < critical[2] & gmat >= critical[1] ] <- 4
  gmat[ gmat < critical[1] & gmat >= 0 ] <- 5

  ##### ##### ##### ##### #####
  # Write to file.
#  
  write.table(gmat, file = "gene_ploidy.csv", sep = ",")
}

gmat <- read.table("gene_ploidy.csv", header = TRUE, sep = ",", row.names = 1)
gmat <- as.matrix(gmat)
nrow(gmat)

##### ##### ##### ##### #####
# Visualize

library(vcfR)
library(RColorBrewer)

heatmap.bp(gmat, col.ramp = brewer.pal(n=4, name = "Dark2"), rlabels = FALSE)


##### ##### ##### ##### #####


phas  <- gmat[,grep("F18", colnames(gmat)), drop = FALSE]
gmat1 <- gmat[,grep("F18", colnames(gmat), invert = TRUE), drop = FALSE]
ipo   <- gmat1[,grep("00Ip5|Ipom_1.2|Ipom_2.4|Pipo5|P3001|PIC99167", colnames(gmat1))]
gmat1 <- gmat1[,grep("00Ip5|Ipom_1.2|Ipom_2.4|Pipo5|P3001", colnames(gmat1), invert = TRUE)]
mir   <- gmat1[,grep("00M_410|CBS.678.85|PIC99167|Pmir5|PIC99114|P7722", colnames(gmat1))]
gmat1 <- gmat1[,grep("00M_410|CBS.678.85|PIC99167|Pmir5|PIC99114|P7722", colnames(gmat1), invert = TRUE)]
para  <- gmat1[,grep("CJ01A1|INRA.310|P10297|P1569|P1976", colnames(gmat1))]
gmat1 <- gmat1[,grep("CJ01A1|INRA.310|P10297|P1569|P1976", colnames(gmat1), invert = TRUE)]
adna  <-  gmat1[,grep("Kew|KM|M.", colnames(gmat1))]
gmat1 <-  gmat1[,grep("Kew|KM|M.", colnames(gmat1), invert = TRUE)]
and   <-  gmat1[,grep("PaX|EC3425|EC3394|P13803", colnames(gmat1))]
gmat1 <-  gmat1[,grep("PaX|EC3425|EC3394|P13803", colnames(gmat1), invert = TRUE)]


gmat1 <- gmat1[, sort.int(colSums(gmat1, na.rm = TRUE)/apply(gmat1, MARGIN=2, function(x){sum(!is.na(x))}), index.return = TRUE)$ix]


# Mexico
regex <- "PIC97136|PIC97146|PIC97335|PIC97442|PIC97750|PIC97785" # Grunwald
#regex <- paste(regex, "P1362|P3681|P3683|P3685|P6629|P6634|PIC97207|PIC97605|PIC97630|P99189|PIC98372", spe="|") # Martin 2016 Mexico
#regex <- paste(regex, "P6635", sep = "|") # Martin 2016 Central Mexico
regex <- paste(regex, "P7036|P8140|P8141|P8143|P8144|P10650|P6633", sep = "|") # Martin 2016 Toluca, Mexico
mex   <-  gmat1[,grep(regex, colnames(gmat1))]


# Raffaele samples
gmat1 <- gmat1[,grep("90128|PIC99189|t30.4", colnames(gmat1), invert = TRUE)]


#heatmap.bp(phas, col.ramp = brewer.pal(n=4, name = "Dark2"), rlabels = FALSE)
#heatmap.bp(ipo, col.ramp = brewer.pal(n=4, name = "Dark2"), rlabels = FALSE)
#
heatmap.bp(mir, col.ramp = brewer.pal(n=4, name = "Dark2"), rlabels = FALSE)

heatmap.bp(para, col.ramp = brewer.pal(n=4, name = "Dark2"), rlabels = FALSE)
#heatmap.bp(adna, col.ramp = brewer.pal(n=4, name = "Dark2"), rlabels = FALSE)
heatmap.bp(and, col.ramp = brewer.pal(n=4, name = "Dark2"), rlabels = FALSE)


heatmap.bp(gmat1, col.ramp = brewer.pal(n=4, name = "Dark2"), rlabels = FALSE)

heatmap.bp(mex, col.ramp = brewer.pal(n=4, name = "Dark2"), rlabels = FALSE)


##### ##### ##### ##### #####
# Constant ploid?


mex[apply(mex,MARGIN = 1, function(x){ sum(x == 2, na.rm = TRUE) == length(x) }),]


regex <- paste(rownames(mex[apply(mex,MARGIN = 1, function(x){ sum(x == 2, na.rm = TRUE) == length(x) }),]), collapse = "$|^")

regex <- paste("^", regex, "$", sep="")


genes$NAME[grep(regex, genes$LOCUS)]


#
mex[apply(mex,MARGIN = 1, function(x){ sum(x == 3, na.rm = TRUE) / length(x) }) > 0.8,]

regex <- paste(rownames(mex[apply(mex,MARGIN = 1, function(x){ sum(x == 3, na.rm = TRUE) / length(x) }) > 0.8,]), collapse = "$|^")
regex <- paste("^", regex, "$", sep="")

genes$NAME[grep(regex, genes$LOCUS)]


#
mex[apply(mex,MARGIN = 1, function(x){ sum(x == 4, na.rm = TRUE) / length(x) }) > 0.45 ,]


#gmat1[apply(gmat1,MARGIN = 1, function(x){ sum(x == 2, na.rm = TRUE) == length(x) }),]

##### ##### ##### ##### #####
#


nrow(gmat)

sum( apply(gmat, MARGIN = 1, function(x){ length(grep("2|3", x)) / length(na.omit(x)) }) == 1, na.rm = TRUE )

sum( apply(gmat1, MARGIN = 1, function(x){ length(grep("2|3", x)) / length(x) }) == 1 )

##### ##### ##### ##### #####

# Slow!
# heatmap(gmat1)

##### ##### ##### ##### #####


pmat <- matrix(ncol=ncol(gmat), nrow = 5)
colnames(pmat) <- colnames(gmat)

#rownames(pmat) <- c("5","4","3","2","NA")
rownames(pmat) <- c("NA","2","3","4","5")

pmat["5",] <- colSums(gmat == 5, na.rm = TRUE)
pmat["4",] <- colSums(gmat == 4, na.rm = TRUE)
pmat["3",] <- colSums(gmat == 3, na.rm = TRUE)
pmat["2",] <- colSums(gmat == 2, na.rm = TRUE)
pmat["NA",] <- apply(gmat, MARGIN = 2, function(x){ sum(is.na(x)) })



phas  <- pmat[,grep("F18", colnames(pmat)), drop = FALSE]
pmat <- pmat[,grep("F18", colnames(pmat), invert = TRUE), drop = FALSE]
ipo   <- pmat[,grep("00Ip5|Ipom_1.2|Ipom_2.4|Pipo5|P3001|PIC99167", colnames(pmat))]
pmat <- pmat[,grep("00Ip5|Ipom_1.2|Ipom_2.4|Pipo5|P3001", colnames(pmat), invert = TRUE)]
mir   <- pmat[,grep("00M_410|CBS.678.85|PIC99167|Pmir5|PIC99114|P7722", colnames(pmat))]
pmat <- pmat[,grep("00M_410|CBS.678.85|PIC99167|Pmir5|PIC99114|P7722", colnames(pmat), invert = TRUE)]
para  <- pmat[,grep("CJ01A1|INRA.310|P10297|P1569|P1976", colnames(pmat))]
pmat <- pmat[,grep("CJ01A1|INRA.310|P10297|P1569|P1976", colnames(pmat), invert = TRUE)]
adna  <-  pmat[,grep("Kew|KM|M.", colnames(pmat))]
pmat <-  pmat[,grep("Kew|KM|M.", colnames(pmat), invert = TRUE)]
and   <-  pmat[,grep("PaX|EC3425|EC3394|P13803", colnames(pmat))]
pmat <-  pmat[,grep("PaX|EC3425|EC3394|P13803", colnames(pmat), invert = TRUE)]


pmat <- pmat[,sort.int(pmat['2',], index.return = TRUE)$ix]

# Mexico
regex <- "PIC97136|PIC97146|PIC97335|PIC97442|PIC97750|PIC97785" # Grunwald
#regex <- paste(regex, "P1362|P3681|P3683|P3685|P6629|P6634|PIC97207|PIC97605|PIC97630|P99189|PIC98372", spe="|") # Martin 2016 Mexico
#regex <- paste(regex, "P6635", sep = "|") # Martin 2016 Central Mexico
regex <- paste(regex, "P7036|P8140|P8141|P8143|P8144|P10650|P6633", sep = "|") # Martin 2016 Toluca, Mexico
mex   <-  pmat[,grep(regex, colnames(pmat))]


sort(colnames(pmat))


# Raffaele samples
pmat <- pmat[,grep("90128|PIC99189|t30.4", colnames(pmat), invert = TRUE)]


pdf('ploidy_bp.pdf')

par(mar=c(9,4,4,2))

#barplot(pmat, las = 3, col=brewer.pal(n=5, name = "Dark2"))
barplot(pmat[-1,], las = 3, col=brewer.pal(n=5, name = "Dark2"))
title(main = "P. infestans")
#abline(h=500)

barplot(ipo[-1,], las = 3, col=brewer.pal(n=5, name = "Dark2"))
title(main = "P. ipomoeae")
#abline(h=500)

barplot(mir[-1,], las = 3, col=brewer.pal(n=5, name = "Dark2"))
title(main = "P. mirabilis")
#abline(h=500)

barplot(para[-1,], las = 3, col=brewer.pal(n=5, name = "Dark2"))
title(main = "P. parasitica")
#abline(h=500)

#barplot(adna[-1,], las = 3, col=brewer.pal(n=5, name = "Dark2"))
#abline(h=500)

barplot(and[-1,], las = 3, col=brewer.pal(n=5, name = "Dark2"))
title(main = "P. andina")
#abline(h=500)

barplot(mex[-1,], las = 3, col=brewer.pal(n=5, name = "Dark2"))
title(main = "P. infestans, MX")
abline(h=500)

par(mar=c(5,4,4,2))

dev.off()


##### ##### ##### ##### #####


png("ploidy_bp.png", width = 800, height = 800)

par(mar=c(9,4,4,2))

barplot(pmat[-1,], las = 3, col=brewer.pal(n=5, name = "Dark2"), space = 0)
title(main = "P. infestans")

par(mar=c(5,4,4,2))

dev.off()



##### ##### ##### ##### #####
