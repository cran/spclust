#' Generate NAM founder genotypes
#'
#' Generates founder genotypes for a Nested Association Mapping population
#' @export
#' @param nparents Total number of parents
#' @param nmarkers Number of markers on all chromosomes
#' @param map Marker map
#' @param option Process for generating genotypes. See details below.
#' @note If option==1 (default), the common parent (1) is assumed to have a different value from all the other parents. This would be achieved if markers were chosen for which that parent had a rare allele. If option==2, then the minor allele frequency for each markers is uniformly distributed on [0, .5] and parents are randomly assigned one of two alleles (0 or 1) according to the allele frequency. Note that in all cases "1" is the minor allele.  The marker map is currently only used for the marker names, but alternative methods of generating founder genotypes could depend on it. If it is not input a default marker naming will be used. 
#' @return A matrix with nparents rows and nmarkers columns, with values in {0, 1}. The attribute "MAF" will give a vector containing minor allele frequencies for each marker (set at 1/nparents if option 1 is selected). 
generateNAMfounders <- function(nparents, nmarkers, map, option=1) 
{
  if (missing(map)) {
	cat("Warning: Marker map not input; will use default marker names\n")
	mrknam <- paste("D1M", 1:nmarkers, sep="")
  } else mrknam <- unlist(lapply(map, names))


  if (option==1) {
	foumat <- rbind(rep(1, nmarkers), matrix(0, nrow=nparents-1, ncol=nmarkers))
	maf <- 1/nparents
  } else if (option==2) {
	maf <- runif(nmarkers, 0, .5)
	foumat <- matrix(runif(nmarkers*nparents, 0, 1), nrow=nparents, ncol=nmarkers)
	foumat <- t(apply(foumat, 1, function(x) as.numeric(x<maf)))
  }
  colnames(foumat) <- mrknam
  attr(foumat, "MAF") <- maf
  return(foumat)
}

