
#' A function to calculate binding strengh by summing weighted chip-seq peaks.
#'
#' @description This function is intended to calculate chip-seq binding score and
#'	associated probability for each gene. The input chip-seq list
#'	contains all peaks found in up/down stream within certain range
#'	whose weights contribute to final score is assigned based on the
#'	its distance to transcription start site.
#' @param ChipList A matrix contains all peaks for each genes with each row
#'		represent one peak and the first column contains ref-seq ID,
#'		second column contain entrez ID, third column contains
#'		chromosome ID, fourth column contains distance from the
#'		center of the peak to transcription start site and the fifth
#'		column contains the peak's strength.
#' @param verbose If TRUE, prints out results of every iteration.
#' @param genome Transcript definition table created by
#'		makeTranscriptDbFromUCSC function and used to create ChipList
#'
#' @return This function returns a matrix with first column contains the
#'	entrez ID, the second column contains the score and the third
#'	column contains the associated probability.
#'
#' @details This function is intended to calculate chip-seq binding score and
#'	associated probability for each gene.
#'
#' @author Mehdi Fazel-Najafabadi, Mario Medvedovic
#'
#' @references ...
#' @seealso ...
#' @keywords annotation
<<<<<<< HEAD
#' @importFrom GenomicFeatures transcripts
#' @importFrom IRanges IRanges findOverlaps
=======
# importFrom GenomicFeatures transcripts
# importFrom IRanges IRanges findOverlaps
>>>>>>> 78f437d6e9df8be758b38ccc911f460bf8ad6de2
#' @export
#' @examples
#' ## not run
#' data(erAlpha)
#' refTable <- GenomicFeatures::makeTxDbFromUCSC(genome="hg19",tablename="refGene")
#' ChipSeq <- annotateChipSeqPeaks(chip.seq=erAlpha[[3]], transcriptDB=refTable, distanceRange=c(-1e+06,1e+06))
#' tregBindingERalpha <- chipSeqWeightedSum(ChipSeq, verbose=TRUE, genome=refTable)
#' ## not run

chipSeqWeightedSum <- function(ChipList,verbose=FALSE, genome=NULL) {
	emExpUniform <- function(data,p=0.1, lamda=0.1, steps=1000,stopCond=0.01) {
		maxValue <- max(data)
		minValue <- min(data)
		likelihood1 <- sum(log(p*lamda*exp(-lamda*data)+(1-p)/(maxValue-minValue)))
		likelihood2 <- 9999
		if(verbose)cat("Log likelihood ",likelihood1, "\nIteration\tLikelihood\tWeights\tLamda\n")
		iter <- 1
		while(abs(likelihood2-likelihood1) > stopCond & iter < steps)
		{
			if(iter > 1) likelihood1 <- likelihood2
			tao <- p*lamda*exp(-lamda*data)/(p*lamda*exp(-lamda*data)+(1-p)/(maxValue-minValue))
			p <- sum(tao)/length(data)
			lamda <- sum(tao)/sum(tao*data)
			likelihood2 <- sum(log(p*lamda*exp(-lamda*data)+(1-p)/(maxValue-minValue)))
			iter <- iter+1
			if(verbose)cat(iter,likelihood2,p,lamda,"\n")
		}
		res <- list(p=p,lamda=lamda)
		res
<<<<<<< HEAD
	}
=======
	} 
>>>>>>> 78f437d6e9df8be758b38ccc911f460bf8ad6de2

    emExpNormal <- function(data,p=0.1,lamda=0.1,mean=2.0,var=1.0,steps=1000,stopCond=0.01) {
	      likelihood1 <- sum(log(p*lamda*exp(-lamda*data)+(1-p)*dnorm(data,mean=mean,sd=sqrt(var))))
	      likelihood2 <- 9999
	     if(verbose) cat("Log likelihood ",likelihood1, "\nIteration\tLikelihood\tWeights\tLamda\tMean\tVariance\n")
	      iter <- 1
	      while(abs(likelihood2-likelihood1) > stopCond & iter < steps)
	      {
		      if(iter > 1) likelihood1 <- likelihood2
		      tao1 <- p*lamda*exp(-lamda*data)/(p*lamda*exp(-lamda*data)+(1-p)*dnorm(data,mean=mean,sd=sqrt(var)))
		      tao2 <- 1-tao1
		      p <- sum(tao1)/length(data)
		      lamda <- sum(tao1)/sum(tao1*data)
		      mean <- sum(tao2*data)/sum(tao2)
		      var <- sum(tao2*(data-mean)^2)/sum(tao2)
		      likelihood2 <- sum(log(p*lamda*exp(-lamda*data)+(1-p)*dnorm(data,mean=mean,sd=sqrt(var))))
		      iter <- iter+1
		      if(verbose)cat(iter,likelihood2,p,lamda,mean,var,"\n")
	      }
	      res <- list(p=p,lamda=lamda,mean=mean,var=var)
	      res
<<<<<<< HEAD
      }
=======
      } 
>>>>>>> 78f437d6e9df8be758b38ccc911f460bf8ad6de2

    emExpNormal.prior <- function(data,p=0.1,lamda=0.1,mean=2.0,var=1.0,steps=20,stopCond=0.01,background=3,probs=0.5) {
	    if(background <= 0) background <- 0.1
	    lamda <- -log(probs)/background
	    mean <- background+1
	    likelihood1 <- sum(log(p*lamda*exp(-lamda*data)+(1-p)*dnorm(data,mean=mean,sd=sqrt(var))))
	    likelihood2 <- 9999
	    if(verbose) cat("Log likelihood ",likelihood1, "\nIteration\tLikelihood\tWeights\tLamda\tMean\tVariance\n")
	    iter <- 1
	    res <- list()
	    while(abs(likelihood2 - likelihood1) > stopCond & iter < steps) {
		    if(iter > 1) likelihood1 <- likelihood2
		    tao1 <- p*lamda*exp(-lamda*data)/(p*lamda*exp(-lamda*data)+(1-p)*dnorm(data,mean=mean,sd=sqrt(var)))
		    tao2 <- 1-tao1
		    p <- sum(tao1)/length(data)
		    mean <- sum(tao2*data)/sum(tao2)
		    var <- sum(tao2*(data-mean)^2)/sum(tao2)
		    likelihood2 <- sum(log(p*lamda*exp(-lamda*data)+(1-p)*dnorm(data,mean=mean,sd=sqrt(var))))
		    res[[iter]] <- list(p=p, lamda=lamda, mean=mean, var=var)
		    if(verbose) cat(iter,likelihood2,p,lamda,mean,var,"\n")
		    if (likelihood2 == Inf) return(res[[iter-1]])
		    iter <- iter+1
	    }
	    res <- list(p=p, lamda=lamda, mean=mean, var=var)
	    res
<<<<<<< HEAD
    }
=======
    } 
>>>>>>> 78f437d6e9df8be758b38ccc911f460bf8ad6de2

    estimate.background <- function(data,para) {
	    minDist <- min(as.numeric(data[,"TSSDist"]))
	    maxDist <- max(as.numeric(data[,"TSSDist"]))
<<<<<<< HEAD

=======
	    
>>>>>>> 78f437d6e9df8be758b38ccc911f460bf8ad6de2
 	    windowsize <- table(data[,"RefID"])
	    windowsize <- max(windowsize)
	    windowsize <- as.integer((maxDist-minDist)/windowsize)
	    res <- 0
	    iter <- minDist
	    bin.start <- seq(minDist,maxDist-windowsize,by=windowsize)
	    bin.end <- seq(windowsize,maxDist,by=windowsize)
	    bins <- IRanges::IRanges(start=bin.start,end=bin.end)
	    query <- IRanges::IRanges(as.numeric(data[,"TSSDist"]),as.numeric(data[,"TSSDist"])+0.5)
	    tables <- as.matrix(IRanges::findOverlaps(query,bins))
	    bin <- unique(tables[,2])
	    peaks <- unlist(sapply(bin,function(x) mean(as.numeric(data[tables[tables[,2]==x,1],"Score"]))))
	    Weights <- para$p*para$lamda*exp(-para$lamda*(bin.start+windowsize/2))
	    Weights <- Weights/(Weights+(1-para$p)/(maxDist-minDist))
	    res <- Weights*peaks
	    res <- sum(res)
<<<<<<< HEAD

=======
	    
>>>>>>> 78f437d6e9df8be758b38ccc911f460bf8ad6de2
	    res
    }

    WeightedSum <- function(data,TSSDist.Bound=10000) {
	    data <- data[as.numeric(data[,"TSSDist"])<=TSSDist.Bound,]
	    para <- emExpUniform(as.numeric(data[,"TSSDist"]),p=0.1,lamda=0.1)
	    if(verbose) cat("Estimate background...")
	    uniform <- estimate.background(data,para)
 	    if(verbose) cat(log(uniform+1),"\n")
	    Weights <- para$p*para$lamda*exp(-para$lamda*as.numeric(data[,"TSSDist"]))
	    Weights <- Weights/(Weights+(1-para$p)/(max(as.numeric(data[,"TSSDist"]))-min(as.numeric(data[,"TSSDist"]))))
	    data <- cbind(data,WeightScore=Weights*as.numeric(data[,"Score"]))
  	    Weighted.Sum <- matrix(0,nrow=length(unique(data[,"RefID"])),ncol=3)
 	    colnames(Weighted.Sum) <- c("RefID","EntrezID","Score")
	    if(verbose) cat("Summarize Chip.Seq peaks...\n")
	    temp <- split(data,data[,"RefID"])
	    Weighted.Sum[,1] <- names(temp)
 	    Weighted.Sum[,2] <- sapply(temp,function(x)x[1,"GeneID"])
 	    Weighted.Sum[,3] <- sapply(temp,function(x) sum(x[,"WeightScore"]))
	    rownames(Weighted.Sum) <- names(temp)
	    Weighted.Sum <- data.frame(Weighted.Sum,stringsAsFactors=F)

	    res <- list(Score=Weighted.Sum,Uniform=uniform)
	    res
    }
<<<<<<< HEAD

=======
  
>>>>>>> 78f437d6e9df8be758b38ccc911f460bf8ad6de2
    WeightedSum.Prob <- function(data,prob=0.005,prior=TRUE,log=TRUE) {
	data.2 <- data[[1]]
	if(log) {
		data.2[,"Score"] <- log(as.numeric(data.2[,"Score"])+1)
		backgrounds <- log(data[[2]]+1)
	} else backgrounds <- data[[2]]
<<<<<<< HEAD

=======
	
>>>>>>> 78f437d6e9df8be758b38ccc911f460bf8ad6de2
	if(prior) {para <- emExpNormal.prior(as.numeric(data.2[,"Score"]), background=backgrounds,probs=prob)
	} else para <- emExpNormal(as.numeric(data.2[,"Score"]))
	Prob <- (1-para$p)*dnorm(as.numeric(data.2[,"Score"]), mean=para$mean,sd=sqrt(para$var))
	Prob <- Prob/(Prob+para$p*para$lamda*exp(-para$lamda*as.numeric(data.2[,"Score"])))
 	data.2 <- cbind(data.2, Prob=Prob)
	rownames(data.2) <- data.2[,1]
	data.2 <- data.2[,-1]
	data.2
      }
	getRefSeqs <- function(reftable, keys="TXNAME", columns=c("GENEID","TXNAME"), keytype="TXNAME") {
	refSeqs <- select(reftable, keys(reftable, keys), columns=columns ,keytype=keytype)
	refSeqs
	}

    if(verbose) cat("Sum Chip-seq peaks...\n")
    Chip.Weighted <- WeightedSum(ChipList,TSSDist.Bound=1000000)
    if(verbose) cat("Get binding probs...\n")
    Chip.Weighted.Prob <- WeightedSum.Prob(Chip.Weighted, prob=0.001, prior=TRUE)
    if(!is.null(genome)){
      print("Assigning zeros to non-bound genes")
 	refseqs <- getRefSeqs(genome)
      matchingAllReseqs <- match(refseqs[,1],rownames(Chip.Weighted.Prob))
      Chip.Weighted.Prob <- Chip.Weighted.Prob[matchingAllReseqs,]
      Chip.Weighted.Prob$Score[is.na(matchingAllReseqs)] <- 0
      Chip.Weighted.Prob$Prob[is.na(matchingAllReseqs)] <- 0
      Chip.Weighted.Prob$EntrezID <- refseqs[, 2]
     rownames(Chip.Weighted.Prob) <- refseqs[, 1]
    }
    Chip.Weighted.Prob
}
