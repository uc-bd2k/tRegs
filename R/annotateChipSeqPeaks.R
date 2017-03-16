
#' A function to match each peak to associated genes in terms of RefID.
#'
#' @description This function allows you to annotate chip-sequencing reads to a reference.
#'	It is intended to search all peaks within up and down
#'	distances for each gene. Users could specify their own platforms
#'	or the function will download platform from UCSC as default.
#' @param chip.seq A matrix listing all CHIP-Seq peaks in which the chromosome
#'		ID for all peaks are listed in the first column, peaks' start
#'		positions are listed in the second column, peaks' end
#'		position are listed in the third column, peaks' ID are listed
#'		in the forth column and the fifth column lists peaks' score.
#' @param transcriptDB This is Transcript definition table created by
#'		makeTranscriptDbFromUCSC or makeTxDbFromUCSC functions.
#' @param distanceRange The up- and downstream window aroung gene TSS to be used
#'		in associating peaks to a gene. The default is c(-1e+06,1e+06).
#'
#' @return This function returns a matrix in which the first column contains
#'	the RefSeqID, the second column contains Entrez ID , the third
#'	column contains chromosome ID, the fourth column contains the
#'	absolute distance from each peak to transcription start site and
#'	the fifth column contains chip-seq score.
#'
#' @details This function is intended to search all peaks within up and down distances for each gene.
#'
#' @author Mehdi Fazel-Najafabadi, Mario Medvedovic
#'
#' @export
#' @examples
#' ## not run
#' data(erAlpha)
#' refTable <- GenomicFeatures::makeTxDbFromUCSC(genome="hg19",tablename="refGene")
#' ChipSeq <- annotateChipSeqPeaks(chip.seq=erAlpha[[3]], transcriptDB=refTable, distanceRange=c(-1e+06,1e+06))
#' ## not run

annotateChipSeqPeaks <- function(chip.seq, transcriptDB=NULL, distanceRange=c(-1e+06,1e+06)) {
	colnames(chip.seq) <- c("Chromosome","Start","End","Des","Score")
	if(is.null(transcriptDB) | class(transcriptDB) != "TxDb") stop("Please provide refGenome platform\n") else {
	    genomeDB <- transcriptDB
	}
	getRefSeqs <- function(reftable, keys="TXNAME", columns=c("GENEID","TXNAME"), keytype="TXNAME") {
	refSeqs <- select(reftable, keys(reftable, keys), columns=columns ,keytype=keytype)
	refSeqs
	}
	refSeqs <- getRefSeqs(genomeDB)
# importFrom GenomicFeatures transcripts 
# importFrom IRanges IRanges findOverlaps NCList values
	allchr <- unique(chip.seq[,1])
	exps = matrix(NA,ncol=5)
 	colnames(exps) <- c("RefID","GeneID","Chromosome","TSSDist","Score")
	for(chr in allchr)
	{
	  cat("Annotating ",chr,"\n")
	  genome <- GenomicFeatures::transcripts(genomeDB, filter = list(tx_name=refSeqs[,"TXNAME"],tx_chrom=chr))
	  genome.starts <- start(genome)
	  genome.ends <- end(genome)
	  genome.strands <- as.vector(strand(genome))
	  genome.0.starts <- genome.starts
	  genome.0.ends <- genome.ends
	  genome.0.starts[genome.strands=="+"] <- genome.starts[genome.strands=="+"]+distanceRange[1]
	  genome.0.starts[genome.0.starts<0] <- 0
	  genome.0.ends[genome.strands=="+"] <- genome.starts[genome.strands=="+"]+distanceRange[2]
	  genome.0.starts[genome.strands=="-"] <- genome.ends[genome.strands=="-"]+distanceRange[1]
	  genome.0.starts[genome.0.starts<0] <- 0
	  genome.0.ends[genome.strands=="-"] <- genome.ends[genome.strands=="-"]+distanceRange[2]
	  gene.id <- refSeqs$GENEID[match(as.character(IRanges::values(genome)[,"tx_name"]), refSeqs$TXNAME)]

	  tree <- IRanges::NCList(IRanges::IRanges(genome.0.starts,genome.0.ends))
	  temp <- chip.seq[chip.seq$Chromosome==chr,]
	  tables <- as.matrix(IRanges::findOverlaps(query=IRanges::IRanges(temp$Start,temp$End),tree,type="within",select="all"))
	  distances <- genome.starts[tables[,2]]
	  distances[is.element(tables[,2],which(genome.strands=="-"))] <- genome.ends[tables[is.element(tables[,2],which(genome.strands=="-")),2]]
	  starts1 <- temp[tables[,1],"Start"]
	  starts2 <- temp[tables[,1],"End"]
	  distances <- as.integer(abs((starts1+starts2)/2-distances))
	  tables<-cbind(tables,Dist=distances,Score=temp[tables[,1],"Score"])
	  rowID <- paste(IRanges::values(genome)[tables[,2],"tx_name"], IRanges::values(genome)[tables[,2],"tx_id"],tables[,3],tables[,1],sep="-")
	  if(length(rowID)>length(unique(rowID))) cat(chr,"\n")
 	  tables <- data.frame(RefID=IRanges::values(genome)[tables[,2],"tx_name"],GeneID=gene.id[tables[,2]],Chromosome=as.character(rep(chr,dim(tables)[1])),
	    TSSDist=tables[,3],Score=tables[,4],row.names=rowID,stringsAsFactors=FALSE)
	  exps <- rbind(exps,tables)
	  rm(tables)
	}
	exps <- exps[!is.na(exps[,1]),]
	exps
}
