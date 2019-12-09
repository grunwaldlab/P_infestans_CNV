
FIR_ggplot <- function(x, bins=30, legend=TRUE, xlab = TRUE, ylab = TRUE, ...){
  require(ggplot2)
  require(hexbin)
  x <- x[!is.na(x$fiveprime),]
  x <- x[!is.na(x$threeprime),]
  x <- x[!x$fiveprime == 0, ]
  x <- x[!x$threeprime == 0, ]
  #x <- ggplot(x,aes(log(fiveprime),log(threeprime)))
  x <- ggplot(x,aes(fiveprime, threeprime))
  x <- x + scale_y_continuous(trans=scales::log_trans(), breaks=10^(1:5))
  x <- x + scale_x_continuous(trans=scales::log_trans(), breaks=10^(1:5))
  x <- x + stat_binhex(bins = bins)
  #x <- x + scale_fill_gradientn(colours = c("white", "darkblue", "forestgreen", "goldenrod1", "orangered", "red3", "darkred"))
  x <- x + scale_fill_gradientn(colours = c("white","slategray1","darkseagreen2","lightgoldenrod1","indianred1","indianred3"))
#  x <- x + xlim(c(0,12.5)) + ylim(c(0,12.5))

  x <- x + labs(x = "5' distance (log bp)")
  x <- x + labs(y = "3' distance (log bp)")
  x <- x + theme_bw()

  if( legend == FALSE ){
    x <- x + theme(legend.position="none")
  }
  if( xlab == FALSE ){
    x <- x + theme(axis.title.x=element_blank(), 
                   axis.text.x=element_blank())
  }
  if( ylab == FALSE ){
    x <- x + theme(axis.title.y=element_blank(), 
                   axis.text.y=element_blank())
  }

  return(x)
}


enrich_FIR <- function(boc=boc, fir=fir, label=NULL, label_size = 10, threshold=0, legend = TRUE){

  # Isolate deleted genes.
  deleted_boc <- apply(boc, MARGIN=1, function(x){ sum(x <= threshold) == length(x) })
  deleted_boc <- deleted_boc[deleted_boc == TRUE]
  deleted_FIR <- fir[names(deleted_boc),]
  deleted_FIR <- deleted_FIR[!is.na(deleted_FIR$fiveprime) & !is.na(deleted_FIR$threeprime),]
  
  # Fisher's exact test.
  subFir <- fir[!is.na(fir$fiveprime) & !is.na(fir$threeprime), ]

  genes <- matrix(
    c(nrow(subFir[subFir$fiveprime <= 1500 & subFir$threeprime <= 1500,]),
      nrow(subFir[subFir$fiveprime > 1500 & subFir$threeprime > 1500,]),
      nrow(deleted_FIR[deleted_FIR$fiveprime <= 1500 & deleted_FIR$threeprime <= 1500,]),
      nrow(deleted_FIR[deleted_FIR$fiveprime > 1500 & deleted_FIR$threeprime > 1500,])
    ), nrow = 2, dimnames = list(c('Present','Absent'), c('Dense','Sparse'))
  )
  feP <- round(fisher.test(genes, alternative = "greater")$p.value, digits = 4)
  
  # Plot
  if( legend == TRUE ){
    p <- FIR_ggplot(fir)
  } else {
    p <- FIR_ggplot(fir, legend = FALSE, xlab = FALSE, ylab = FALSE)
  }
  p <- p + geom_point( data = deleted_FIR, aes(fiveprime,threeprime), colour="#00000088",size=1)
  
  if( !is.null(label) ){
    p <- p + geom_rect(aes(xmin = 1e0, xmax=1e6, ymin=1, ymax=1e1), colour="#808080", fill="#808080", alpha = 0.003)
    p <- p + geom_text(aes(x=1e3,y=3),label=label,color="white", size = label_size)
  }
#p <- p + geom_text(aes(x=3e0,y=17),label="n=14",color="black", size = 6)
  p <- p + geom_text(aes(x=4e0,y=30),label=paste("n=", nrow(deleted_FIR), sep = ""),color="black", size = 3)
  p <- p + geom_text(aes(x=1e5,y=30),label=paste("p-value=", feP, sep = ""),color="black", size = 3)
  return(p)
}




