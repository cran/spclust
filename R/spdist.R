#' Compute genetic distances for experimental crosses
#'
#' This function computes the expected percentage of markers shared IBD across
#' the genome for experimental crosses of type backcross, doubled haploid, 
#' recombinant inbred line ("bc"), F2 ("f2"), Multiparent Advanced 
#' Generation Inter-Cross with 4 or 8 parents ("magic"), or Nested Association
#' Mapping ("nam"). 
#'
#' @useDynLib spclust
#' @param object cross or mpcross object containing genetic data
#' @param type type of experimental cross
#' @export spdist
#' @seealso \code{\link[spclust]{spclust}}
#' @return distance object for use in spclust
#' @examples
#' # Simulate a map and data using qtl package
#' map <- sim.map(len=100, n.mar=101, eq.spacing=TRUE, include.x=FALSE)
#' dat1 <- sim.cross(map, n.ind=200, type="bc")
#' dat2 <- sim.cross(map, n.ind=500, type="f2")
#' # Calculate the distances between lines in datasets
#' dist1 <- spdist(dat1, "bc")
#' dist2 <- spdist(dat2, "f2")
#' # Heatmap of distances
#' heatmap(as.matrix(dist1), Rowv=NA, Colv=NA, scale="none", col=topo.colors(10))
#' par(mfrow=c(2, 1))
#' # Histogram of minimum distances between lines for full datasets
#' d2 <- as.matrix(dist2) 
#' diag(d2) <- NA
#' mind <- apply(d2, 1, function(x) min(x, na.rm=TRUE))
#' hist(mind, col="tomato", main="Minimum distances between all F2 lines", xlab="Distance")
#' # Histogram of minimum distances between selected lines
#' sp <- spclust(dat2, 100, method="ward")
#' hist(sp$mind[, 2], col="tomato", main="Minimum distances between 100 selected F2 lines", xlab="Distance")

spdist <- function(object, type=c("bc", "f2", "magic", "nam")) 
{
	if (type=="bc")
	{
		if (!inherits(object, "cross")) stop("Must input object of class cross to compute BC distance.\n")
		crp <- calc.genoprob(object)
		prmat <- lapply(crp$geno, function(x) x$prob[,,1])
		prmat <- do.call("cbind", prmat)
		rownames(prmat) <- paste("L", 1:nrow(prmat), sep="")
		distmat <- .C("distbc", nrow(prmat), ncol(prmat), as.numeric(t(prmat)), as.numeric(1-t(prmat)), distmat=numeric(nrow(prmat)*nrow(prmat)), PACKAGE="spclust")$distmat
	}

	if (type=="f2")
	{
		if (!inherits(object, "cross")) stop("Must input object of class cross to compute F2 distance.\n")
		crp <- calc.genoprob(object)
		prmat <- lapply(crp$geno, function(x) x$prob[,,1])
		prmat <- do.call("cbind", prmat)
		prmat2 <- lapply(crp$geno, function(x) x$prob[,,2])
		prmat2 <- do.call("cbind", prmat2)
		rownames(prmat) <- rownames(prmat2) <- paste("L", 1:nrow(prmat), sep="")
		distmat <- .C("distf2", nrow(prmat), ncol(prmat), as.numeric(t(prmat)), as.numeric(t(prmat2)), as.numeric(1-t(prmat)-t(prmat2)), distmat=numeric(nrow(prmat)*nrow(prmat)), PACKAGE="spclust")$distmat
	}

	if (type=="magic")
	{
		if (!inherits(object, "mpcross")) stop("Must input object of class mpcross to compute MAGIC distance.\n")
		mpp <- mpprob(object, program="qtl", est=FALSE)
		prmat <- lapply(mpp$prob, function(x) x)
		prmat <- do.call("cbind", prmat)
		rownames(prmat) <- rownames(mpp$finals)
		distmat <- .C("distmagic", nrow(prmat), ncol(prmat), nrow(mpp$founders), as.numeric(t(prmat)), distmat=numeric(nrow(prmat)*nrow(prmat)), PACKAGE="spclust")$distmat
	}

	if (type=="nam")
	{
		if (!inherits(object, "nam")) stop("Must input object of class nam to compute NAM distance.\n")
		## think about how to do this - do we want to treat each 
		## founder allele separately or consider just as 1s/0s??
		## need to make sure same full map is input for all 
		nfamilies <- length(object)
		prmat <- vector()
		for (i in 1:nfamilies) {
		  crp <- calc.genoprob2(object[[i]], object$map)
		  pr <- lapply(crp$geno, function(x) x$prob[,,1])
		  pr <- do.call("cbind", pr)
		  rownames(pr) <- attr(object[[i]], "linenames")
		  prmat <- rbind(prmat, pr)
		}
		### this gives the distance matrix for lines within a family
		### what about between families?
		  distmat <- .C("distbc", nrow(prmat), ncol(prmat), as.numeric(t(prmat)), as.numeric(1-t(prmat)), distmat=numeric(nrow(prmat)*nrow(prmat)), PACKAGE="spclust")$distmat
	}

	output <- matrix(distmat, nrow=nrow(prmat), ncol=nrow(prmat), byrow=TRUE)
	output <- as.dist(output)
	attr(output, "linenames") <- rownames(prmat)
	output
}

