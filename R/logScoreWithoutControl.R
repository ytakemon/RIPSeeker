# Function Name: 	logScoreWithoutControl
# Description: 		x
# Input: 			y
# Output:			z
#
# Author: Yue Li
###############################################################################

logScoreWithoutControl <- function(nbhGRRIP, padjMethod="BH")
{	
	
	stopifnot(!missing(nbhGRRIP))
			
	logOddScore <- computeLogOdd(nbhGRRIP)
		
	values(nbhGRRIP) <- cbind(as.data.frame(values(nbhGRRIP)), logOddScore)
	
#	if(existsFunction("mcols")) values(nbhGRRIP) <- cbind(mcols(nbhGRRIP), logOddScore)
	
	enrichIdx <- which(values(nbhGRRIP)$viterbi_state == 2)
	
	# assign 1 as dummy index if no enriched state present
	if(length(enrichIdx) == 0) enrichIdx <- 1
		
	unmergedRIP <- nbhGRRIP[enrichIdx]	
	
	mergedRIP <- reduce(unmergedRIP, 
			min.gapwidth = median(width(unmergedRIP) ))
	
	
	overlapIdx <- findOverlaps(mergedRIP, unmergedRIP)
	
	
	mergedRIPList <- lapply(split(overlapIdx, queryHits(overlapIdx)),
			scoreMergedBins, unmergedRIP, mergedRIP)
	
	
	names(mergedRIPList) <- NULL
	
	
	mergedRIP <- do.call(c, mergedRIPList)
	
		
	# compute p-value of the averaged logOddScore
	pval <- pnorm(values(mergedRIP)$logOddScore, mean(logOddScore), 
			sd(logOddScore), lower.tail=FALSE)
			
	pvalAdj <- p.adjust(pval, padjMethod)
	
	values(mergedRIP) <-  cbind(as.data.frame(values(mergedRIP)), pval, pvalAdj)
	
#	if(existsFunction("mcols")) values(mergedRIP) <-  cbind(mcols(mergedRIP), pval, pvalAdj)
	
	return(mergedRIP)	
}
