# Function Name: 	galp2gal
# Description: 		Convert GappedAlignmentPairs to GappedAlignments using CIGAR to mark flanked portion of the pairs as 'N'					
# Input: 			GappedAlignmentPairs
# Output:			GappedAlignments
#
# Author: Yue Li
###############################################################################

galp2gal <- function(galp)
{	
	betweenPairCigar <- paste(start(right(galp)) - end(left(galp)) + 1, "N", sep="")
	
	galcigar <- paste(cigar(left(galp)), betweenPairCigar, cigar(right(galp)), sep="")
	
	gal <- GappedAlignments(
				seqnames = seqnames(galp),
			
				pos = start(left(galp)),
			
				cigar = galcigar,

				strand = strand(left(galp)),
			
				names = names(left(galp)),
				
				seqlengths = seqlengths(galp))
	return(gal)
}
