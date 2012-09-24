#' Generates NAM pedigree
#'
#' Generated a pedigree for a Nested Association Mapping population 
#' @export
#' @param nparents Total number of parents
#' @param nperfam Number of lines per family
#' @param nssdgen Number of generations of selfing
#' @return A four column matrix with IDs for each individual, maternal parent ID, paternal parent ID, and a flag for whether the individual was observed.  
#' @note The first founder will be the common parent, assumed for all F1s to be the maternal parent. Default parameters are those for the original maize NAM population (ref).  

generateNAMpedigree <- 
function(nparents=26, nperfam=200, nssdgen=6)
{
  obs <- vector()
  # start with founders
  ped <- cbind(1:nparents, matrix(0, nrow=nparents, ncol=3))
  
  ## create F1 generation
  ped <- rbind(ped, cbind(nparents+1:(nparents-1), rep(1, nparents-1), 2:nparents, rep(0, nparents-1)))


  ## create families - self each F1 for nssdgen generations
  if (nssdgen>0) {
  for (i in 1:(nparents-1)) {
     ## first selfing event
     ped <- rbind(ped, cbind(nrow(ped)+1:nperfam, rep(i+nparents, nperfam), rep(i+nparents, nperfam), rep(nssdgen==1, nperfam)))

     ## remaining selfing events
     if (nssdgen>1)
     for (j in 2:nssdgen) 
	ped <- rbind(ped, cbind(nrow(ped)+1:nperfam, nrow(ped)-nperfam+1:nperfam, nrow(ped)-nperfam+1:nperfam, rep(nssdgen==j, nperfam)))

    } 
  }  else { stop("Must have at least one selfing event\n")}  

  ped <- as.data.frame(ped)
  names(ped) <- c("id", "m", "f", "obs")
  return(ped)
}

