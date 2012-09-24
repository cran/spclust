#' Convert a NAM mpcross object to a set of cross objects
#' 
#' Takes a Nested Association Mapping population formatted as an mpcross object and converts it to a list of RIL crosses with (potentially) different sets of segregating markers
#' @export
#' @param mpcross NAM mpcross object 
#' @return A list of cross objects - one for each family in the NAM. Potentially each cross has a different number of lines and/or markers. 

NAM2cross <- function(mpcross) 
{
  output <- list()

  ## Record which family each line belongs to
  fam <- vector(length=nrow(mpcross$finals))
  for (i in 1:nrow(mpcross$finals)) {
     p <- mpcross$pedigree[mpcross$id[i],3]
     while(mpcross$pedigree[p,3]!=0)
       p <- mpcross$pedigree[p,3]
     fam[i] <- p-1
  } 
  
  nfamilies <- sum(mpcross$pedigree[,2]==0 & mpcross$pedigree[,3]==0)-1

  ## For each family need to do separate processing to remove non-polymorphic
  ## markers, the appropriate lines, etc. 
  for (i in 1:nfamilies) {
    ### export to R/qtl format
    index <- which(fam==i)
    if (!is.null(mpcross$pheno)) ph <- mpcross$pheno[index,1] else ph <- rnorm(length(index))
    mat <- data.frame(name=c("ph", unlist(lapply(mpcross$map, names))), chr=c("", rep(names(mpcross$map), unlist(lapply(mpcross$map, length)))), pos=c("", unlist(mpcross$map)))
    dat <- mpcross$finals[index,]
    dat <- apply(dat, 1, function(x) 2*(x==mpcross$founders[1+i,])+(x==mpcross$founders[1,]))
    dat[dat==1] <- "AA"
    dat[dat==2] <- "AB"
    dat <- rbind(ph, dat)

    ### should also recode values so that 1 = common parent; 2 = other
    keep <- which(mpcross$founders[1,] != mpcross$founders[1+i,])
    out <- cbind(mat, dat)[c(1, 1+keep), ]

    write.csv(out, file="tmp.csv", row.names=FALSE, quote=FALSE)  

    cr <- read.cross(format="csvr", file="tmp.csv", genotypes=c("AA", "AB"), skip=1)
    attr(cr, "linenames") <- rownames(mpcross$finals)[index]
    output[[i]] <- cr
  }
  names(output) <- paste("F", 1:nfamilies, sep="")

  attr(output, "map") <- mpcross$map
  attr(output, "type") <- "nam"
  return(output)
}

