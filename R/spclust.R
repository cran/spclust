#' Perform single-stage selective phenotyping clustering
#'
#' This function implements the SPCLUST algorithm to perform selective
#' phenotyping in experimental crosses by maximizing the genetic diversity
#' in the selected subsample. The plot function plots the clusters from spclust with some summary information.
#' For all clustering methods the silhouette plot is displayed for clusters 
#' containing more than one individual. For hierarchical clustering methods
#' the dendrogram is displayed with clusters and selected individuals marked.
#'
#' @useDynLib spclust
#' @rdname spclust1-all
#' @aliases spclust plot.spclust summary.spclust print.spclust
#' @param object Cross or mpcross object containing genetic data; for summary, spclust object
#' @param x spclust object to plot
#' @param nlines Size of selected subsample
#' @param file Optional argument, filename for outputting clusters to a file
#' @param method Clustering method
#' @param ... Additional argument passed through to plclust 
#' @export spclust 
#' @usage spclust (object, nlines, file, method=c("average", "ward" , "pam"))
#' \method{plot}{spclust} (x, ...)
#' \method{print}{spclust} (x, ...)
#' \method{summary}{spclust} (object, ...)
#' @return list with components:
#' \item{lines}{indices of selected lines}
#' \item{mind}{For each selected line, minimum distance to other lines in sample}
#' \item{tree}{Hierarchical clustering tree}
#' \item{clusters}{Assignment of all lines to clusters}
#' @seealso \code{\link[spclust]{spdist}}, \code{\link[cluster]{plot.silhouette}}, \code{\link[spclust]{plclust_in_colour}}
#' @note The SPCLUST algorithm performs the following steps in order to 
#' select a subsample with high genetic diversity. First, genetic distances
#' are estimated between all lines in the sample, based on the expected
#' proportion of alleles not shared IBD across the genome. Second, lines are
#' clustered based on the genetic distance, with the number of clusters 
#' matching the number of lines desired for selection. Third, a representative
#' line is selected from each cluster as the one most similar to other lines 
#' in the cluster. 
#' @examples
#' # Simulate a map and data using the qtl package
#' map <- sim.map(len=rep(100, 5), n.mar=21, eq.spacing=TRUE, include.x=FALSE)
#' dat <- sim.cross(map, n.ind=500, type="f2")
#' # Selection of 100 lines
#' sp <- spclust(dat, 100, method="ward")
#' summary(sp)
#' plot(sp)


spclust <- function(object, nlines, file, method=c("average", "ward", "pam"))
{
	output <- list()
	if (!(inherits(object, "cross") | inherits(object, "mpcross"))) stop("Need to have object of type cross or mpcross for selection. \n")
	if (missing(nlines)) stop("Need to know how many lines to select.\n")
	if (inherits(object, "mpcross")) type <- "magic"
	else if (inherits(object, "cross")) type <- class(object)[1]
	
	if (!(type %in% c("bc", "f2", "magic"))) stop("Invalid type of experimental cross has been input.\n")

	if (inherits(object, "cross"))  nmrk <- sum(nmar(object)) else nmrk <- ncol(object$finals)

	distmat <- spdist(object, type)
	linenames <- attr(distmat, "linenames")
	d <- as.matrix(distmat)

	if (method == "pam") {
		hc <- pam(d, k=nlines, diss=TRUE, keep.diss=TRUE)
		keep <- as.numeric(hc$medoids)
		cl <- hc$clustering
	}
	else if (method == "ward" | method == "average") {
	 	hc <- hclust(distmat, method=method)
		cl <- cutree(hc, k=nlines)
  		keep <- unlist(by(1:nrow(d), cl, function(x) if (length(x)==1) return(x) else x[which.min(apply(d[x,x], 2, sum))]))
		hc$clustering <- cl
	}
	else stop("Clustering method is not available. \n")

	dkeep <- d[keep, keep]
  	diag(dkeep) <- NA
  	mind <- apply(dkeep, 2, function(x) min(x, na.rm=T))
  	output$lines <- linenames[keep]
  	output$mind <- data.frame(lines=linenames[keep], min.dist=mind)
  	output$result <- hc
	cl2 <- cbind(linenames, cl)
	colnames(cl2) <- c("LineName", "Cluster")
  	cl2 <- as.data.frame(cl2)
	cl2 <- cl2[order(cl2[,2]),]
  	output$clusters <- cl2
	
  	if (!missing(file)) write.csv(cbind(1:nlines, output$lines), file=file, row.names=FALSE)

	output$dmatrix <- d
	attr(output, "type") <- type
	attr(output, "method") <- method
	attr(output, "nmrk") <- nmrk
	class(output) <- "spclust"
	output
}

#' @S3method plot spclust
#' @method plot spclust
plot.spclust <- function(x, ...)
{
	pp <- par(ask=TRUE)
	on.exit(par(pp))

	if (attr(x, "method") %in% c("ward", "average")) {
		lab <- rep("", nrow(x$clusters))
		lab[x$lines] <- x$mind[,1]
		lab.col <- rep("black", nrow(x$clusters))
		lab.col[x$lines] <- "red"
		plclust_in_colour(x$result, lab=lab, lab.col=lab.col, ...)
		rect.hclust(x$result, k=length(x$lines), border="blue")
	}
	## silhouette plot
	sil <- silhouette(x$result$clustering, dmatrix=x$dmatrix)
	## modify sil to remove clusters with only one individual
	a <- which(sil[,1] %in% names(table(sil[,1]))[table(sil[,1])>1])
	sil1 <- silhouette(sil[a,1], x$dmatrix[a,a])
	clnames <- names(table(sil1[,1]))
	clcol <- rep("firebrick", nrow(sil1))
	clcol[which(sil1[,1] %in% clnames[seq(2,length(clnames), 2)])] <- "seagreen4"
	plot(sil1, col=clcol, main="Silhouette plot of clusters with size > 1")
}

#' @S3method summary spclust
#' @method summary spclust
summary.spclust <-
function(object, ...)
{
  cat("-------------------------------------------------------\n")
  cat("Summary of spclust object\n")
  cat("-------------------------------------------------------\n")
  cat("Type of experimental cross: ", attr(object, "type"), "\n")
  cat("Number of genetic markers used: ", attr(object, "nmrk"), "\n")
  cat("Original population size: ", nrow(object$clusters), "\n")
  cat("Selected sample size: ", length(object$lines), "\n")
  cat("Clustering method: ", attr(object, "method"), "\n")
  cat("Range of minimum distances between selected individuals: ", round(range(object$mind[,2]),3), "\n")
}

#' @method print spclust
#' @S3method print spclust

print.spclust <-
function(x, ...)
{
  cat(" This is an object of class \"spclust\".\n")
  cat(" We provide the following summary. \n")
  summary(x)
}


