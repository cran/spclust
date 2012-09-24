#' Perform selective phenotyping clustering
#'
#' This function implements the SPCLUST algorithm to perform selective
#' phenotyping in experimental crosses by maximizing the genetic diversity in 
#' the selected subsample. Selection can be done in one or multiple stages. 
#' The plot function plots the clusters with some summary information. 
#' Graphical genotypes are displayed for individuals selected with maximal
#' recombinations. For hierarchical clustering methods the dendrogram is 
#' displayed with clusters and selected individuals marked.
#'
#' @useDynLib spclust
#' @rdname spclust-all
#' @aliases spclust plot.spclust summary.spclust print.spclust
#' @param object Cross or mpcross object containing genetic data; for summary function, spclust object
#' @param x spclust object input to plot function
#' @param nlines Number of lines to be selected (note: does not include the number of input lines)
#' @param method Selection method - options include hierarchical clustering (average, ward), partitioning around medoids (pam), or based on the maximal number of recombinations (maxrec) 
#' @param inputlines Names of lines which must be included in the selected sample. See details below. 
#' @param file Optional argument, filename for outputting clusters to a file
#' @param step Step size used in estimating recombinations (default=5 cM)
#' @param threshold Threshold used in estimating recombinations (default=0.7)
#' @param type Style of plot to draw; 1=Silhouette; 2=Dendrogram; 3=Recombinations; 4=All (that are appropriate)
#' @param ... Additional arguments to be passed on to plot functions
#' @usage spclust (object, nlines, method=c("average", "ward", "pam", "maxrec"), inputlines=NULL, file, step=5, threshold=.7) 
#' \method{plot}{spclust} (x, type=4, ...) 
#' \method{print}{spclust} (x, ...)
#' \method{summary}{spclust} (object, ...)
#' @export spclust
#' @return list with components:
#' \item{numlines}{indices of selected lines from both stages}
#' \item{lines}{names of selected lines from both stages}
#' \item{mind}{For each selected line, minimum distance to other lines in sample}
#' \item{tree}{Hierarchical clustering tree}
#' \item{clusters}{Assignment of all lines to clusters}
#' \item{recmat}{If method="maxrec", returns matrix of recombinations for genomic region}
#' @seealso \code{\link[spclust]{spdist}}, \code{\link[spclust]{plot.spclust}}, \code{\link[spclust]{spclust}}, \code{\link[spclust]{plclust_in_colour}}
#' @note This function can perform both single-stage or multi-stage 
#' selective phenotyping clustering. In a single stage, the SPCLUST 
#' algorithm performs the following steps in order to 
#' select a subsample with high genetic diversity. First, genetic distances
#' are estimated between all lines in the sample, based on the expected
#' proportion of alleles not shared IBD across the genome. Second, lines are
#' clustered based on the genetic distance, with the number of clusters 
#' matching the number of lines desired for selection. Third, a representative
#' line is selected from each cluster as the one most similar to other lines 
#' in the cluster. 
#' 
#' If the inputlines argument is used, SPCLUST performs the following steps 
#' in order to select a sample with high genetic diversity while accounting 
#' for the input lines. First, genetic distances are estimated
#' between all lines in the sample, based on the expected proportion of 
#' alleles not shared IBD across the genome. Second, if the "maxrec" method
#' is selected, nlines lines are selected with the highest number of estimated
#' recombinations, excluding those which have already been selected in the 
#' first  stage. Otherwise, all lines are clustered based on the
#' genetic distance, with the number of clusters equal to
#' the sum of the number of input lines and the nlines argument.
#' All input lines selected in stage 1 are included in
#' the final sample, and clusters containing these lines are
#' excluded from further selection. From nlines of the
#' remaining clusters, a representative is selected as the
#' one most similar to other lines in the cluster.
#' @examples
#' # Simulate a map and data using qtl package
#' map <- sim.map(len=rep(100, 5), n.mar=21, eq.spacing=TRUE, include.x=FALSE)
#' dat <- sim.cross(map, n.ind=500, type="bc")
#' # Select two samples of size 100 in two stages
#' sp <- spclust(dat, 100, method="ward")
#' sp2 <- spclust(dat, 100, method="maxrec", inputlines=sp$lines)
#' summary(sp2)
#' plot(sp2)

spclust <- function(object, nlines, method=c("average", "ward", "pam", "maxrec"), inputlines=NULL, file, step=5, threshold=.7)
{
  if (missing(object)) stop("Cannot select lines without input genetic data\n")
      
  output <- list()
#  output$stage1 <- spobj
#  output$stage2 <- list()

#  if (method=="stage1") attr(output$stage2, "method") <- attr(output$stage1, "method") else attr(output$stage2, "method") <- method

  if (!is.null(attr(object, "type")) && attr(object, "type")=="nam") type <- "nam" 
  else if (inherits(object, "mpcross")) type <- "magic"
  else if (inherits(object, "cross")) type <- class(object)[1]

  if (!(type %in% c("bc", "f2", "magic", "nam"))) stop("Invalid type of experimental cross has been input.\n")

  if (method=="maxrec" & type=="f2") stop("Cannot currently use this method with this type of cross\n")

  if (inherits(object, "cross"))  nmrk <- sum(nmar(object)) 
  else if (attr(object, "type")=="nam") nmrk <- length(unlist(attr(object, "map"))) 
  else nmrk <- ncol(object$finals)

  attr(output, "nmrk") <- nmrk

  distmat <- spdist(object, type)
  linenames <- attr(distmat, "linenames")
  d <- as.matrix(distmat)

  if (is.null(inputlines)) {
        cat("No required lines input; will only select a single-stage sample\n")
        ss1 <- 0
  } else {
	ss1 <- length(inputlines)
  	if (is.character(inputlines)) inputlines <- match(inputlines, linenames)
  	if (any(is.na(inputlines))) stop("Some lines to be included are not in dataset\n")
  }

  if (missing(nlines)) {
        cat("Will select equal number of lines to those input\n")
        nlines <- ss1
	if (nlines==0) stop("No sample size input for selection\n")
  }

  if (method=="maxrec") {
	## need to separate out magic from other types
	if (type=="magic") {
  	  if (!inherits(object, "mpprob")) {
	  mpp <- mpprob(object, step=step, program="qtl", threshold=threshold)
	  cat("mpprob object not input; computing multi-point probabilities with input step and threshold values\n")
	  } else mpp <- object
	gen <- do.call("cbind", lapply(mpp$estfnd, function(x) cbind(x, rep(0, nrow(x)))))
	nrec <- do.call("cbind", lapply(mpp$estfnd, function(x) apply(x, 1, function(y) return(sum(diff(y[!is.na(y)])!=0)))))
	}
	
	if (type=="bc") {
	   gen <- lapply(object$geno, function(x) x$data)
	   nrec <- do.call("cbind", lapply(gen, function(x) apply(x, 1, function(y) return(sum(diff(y[!is.na(y)])!=0)))))
	   gen <- do.call("cbind", lapply(gen, function(x) cbind(x, rep(0, nrow(x)))))
	}

	if (type=="nam") {
	   nfamilies <- length(object)
	   gen <- list()
	   nrec <- list()
	   for (i in 1:nfamilies) {
		gen[[i]] <- lapply(object$geno, function(x) x$data)
	  	nrec[[i]] <- lapply(gen[[i]], function(x) apply(x, 1, function(y) return(sum(diff(y[!is.na(y)])!=0))))
	   	gen[[i]] <- do.call("cbind", lapply(gen[[i]], function(x) cbind(x, rep(0, nrow(x)))))
		nrec[[i]] <- do.call("cbind", nrec[[i]])
	   }
	   gen <- do.call("rbind", gen)
	   nrec <- do.call("rbind", nrec)
 	}

  	totnrec <- apply(nrec, 1, sum)
	ord <- data.frame(LineName=linenames, No.Rec=totnrec) 
  	ord <- ord[order(ord[,2], decreasing=TRUE),]

	keep <- inputlines
  	start <- 1

   	### need to select most recombinant lines
  	while (length(keep) < ss1+nlines) {
    	  rem <- ss1+nlines-length(keep)
     	  keep <- unique(c(keep, ord[start:(start+rem-1),1]))
     	  start <- start+rem
  	}
	allkeep <- keep
	keep <- setdiff(keep, inputlines)
	sel <- rep(0, nrow(ord))
	sel[match(linenames[inputlines], ord[,1])] <- 1
	sel[setdiff(keep, inputlines)] <- 2
	ord$StageSelected <- sel
	output$result <- ord
	# last column is all 0s, omit
	output$recmat <- gen[, 1:(ncol(gen)-1)]
	attr(output, "step") <- step
	attr(output, "threshold") <- threshold
	if (type=="bc") n.founders <- 2 
	else if (type=="nam") n.founders <- length(object)+1
	else n.founders <- nrow(mpp$founders)
	attr(output, "n.founders") <- n.founders 
  }

  ## need to modify this so that it works properly when stage1 method = pam
  if (method %in% c("ward", "average", "pam")) {
    if (method == "pam") {
	hc <- pam(d, k=nlines, diss=TRUE, keep.diss=TRUE)

	reps <- hc$medoids
	clust <- vector(length=ss1+nlines)
	if (ss1>0)
	for (i in 1:(ss1+nlines)) clust[i] <- length(intersect(inputlines, which(hc$clustering==i)))

	ncl <- sum(clust==0)
     	### sample from clusters which do not contain input lines
     	clust2 <- sample(which(clust==0), nlines)
     	keep <- as.numeric(reps[clust2])

     	allkeep <- c(keep, inputlines)

     	sel <- rep(0, nrow(d))
     	sel[keep] <- 2
     	sel[inputlines] <- 1
     	clall <- data.frame(LineName=linenames, Cluster=hc$clustering, StageSelected=sel)
     	clall <- clall[order(clall[,2]),]

       	cl2 <- clall[clall$Cluster%in%clust2,]
   
     	output$result <- hc
     	output$clusters <- clall
    }
    if (method != "pam") {
      hc <- hclust(distmat, method=method)
      cl <- cutree(hc, k=ss1+nlines)
      hc$clustering <- as.numeric(factor(cl, levels=unique(cl[hc$order])))

      reps <- unlist(by(1:nrow(d), cl, function(x) if (length(x)==1) return(x) else x[which.min(apply(d[x,x], 2, sum))]))
     ### count how many clusters do not contain lines from stage 1
     clust <- vector(length=ss1+nlines)
     if (ss1>0)
     for (i in 1:(ss1+nlines)) clust[i] <- length(intersect(inputlines, which(cl == i)))
     ncl <- sum(clust==0)
     ### sample from clusters which do not contain lines from stage 1
     clust2 <- sample(which(clust==0), nlines)
     keep <- reps[clust2]

     allkeep <- c(keep, inputlines)

     sel <- rep(0, nrow(d))
     sel[keep] <- 2
     sel[inputlines] <- 1
     clall <- data.frame(LineName=linenames, Cluster=cl, StageSelected=sel)
     clall <- clall[order(clall[,2]),]

     cl2 <- clall[clall$Cluster%in%clust2,]

     ## can also produce a clustering matrix for stage 2 lines
     
#     output$stage2$clusters <- cl2
     output$result <- hc
     output$clusters <- clall
     } 
  }

  dkeep <- d[keep, keep]
  diag(dkeep) <- NA
  mind <- apply(dkeep, 2, function(x) min(x, na.rm=TRUE))

  dallkeep <- d[allkeep, allkeep]
  diag(dallkeep) <- NA
  mindall <- apply(dallkeep, 2, function(x) min(x, na.rm=TRUE))

#  output$stage2$mind <- data.frame(lines=linenames[keep], min.dist=mind)
#  output$stage2$lines <- linenames[keep]
#  output$stage2$numlines <- keep
  output$mind <- data.frame(lines=linenames[allkeep], min.dist=mindall) 
  output$lines <- linenames[allkeep]
  output$numlines <- allkeep
  output$inputlines <- inputlines
  output$dmatrix <- d

  if (!missing(file)) write.csv(cbind(1:(ss1+nlines), output$lines), file=file, row.names=FALSE)

  attr(output, "method") <- method
  class(output) <- "spclust"
  output
}

#' @method plot spclust
#' @S3method plot spclust

plot.spclust <- function(x, type=4, ...)
{
	if (!type %in% c(1:4)) stop("Type argument must be 1, 2, 3 or 4\n")

	pp <- par(ask=TRUE)
	on.exit(par(pp))

	## Alter this so that stage1 and stage2 colours are different
	if (attr(x, "method") %in% c("ward", "average") & type %in% c(2,4)) {
		lab <- rep("", nrow(x$clusters))
		lab.col <- rep("blue", nrow(x$clusters))
#		lab[match(x$lines, x$clusters[,1])] <- x$lines
		lab[match(x$lines, x$clusters[,1])] <- "*"
		lab.col[match(x$lines, x$clusters[,1])] <- "red"

		if (!is.null(x$inputlines)) 
		  lab.col[match(x$inputlines, x$clusters[,1])] <- "black"
		
		plclust_in_colour(x$result, lab=lab, lab.col=lab.col, ...)
		rect.hclust(x$result, k=length(x$lines), border="blue")
	}

	## silhouette plot
	if (attr(x, "method")!="maxrec" & type %in% c(1, 4)) {
	 	sil1 <- silhouette(x$result$clustering[x$result$order], dmatrix=x$dmatrix[x$result$order, x$result$order])
	 	## modify sil to remove clusters with only one individual
#	 	a <- which(sil[,1] %in% names(table(sil[,1]))[table(sil[,1])>1])
#	 	sil1 <- silhouette(sil[a,1], x$dmatrix[a,a])
	 	clnames <- unique(sil1[,1])
	 	clcol <- rep("firebrick", nrow(sil1))
	 	clcol[which(sil1[,1] %in% clnames[seq(2,length(clnames), 2)])] <- "seagreen4"
	 	plot(sil1, col=clcol, main="Silhouette plot of clusters", ...)
	}
	
	## Recombination diagram
	if (attr(x, "method")=="maxrec" & type %in% c(3, 4)) {
	  	n.founders <- attr(x, "n.founders")
		eightcol <- c("lightblue2", "royalblue1", "darkseagreen1", "seagreen3", "lemonchiffon", "darkgoldenrod1", "indianred1", "firebrick")
		if (n.founders<=8) col <- eightcol else col <- rainbow(n.founders)
		mat <- x$recmat[match(x$lines, x$result[,1]),]
	  	nrec <- apply(mat, 1, function(x) return(sum(diff(x[!is.na(x)])!=0)))
	  	rownames(mat) <- x$lines
	  	heatmap(mat[order(nrec, apply(mat, 1, function(x) max(x, na.rm=TRUE)), decreasing=FALSE),], col=col[1:n.founders], main="Recombinations for selected lines", Rowv=NA, Colv=NA, scale="none", ylab="Lines", xlab="Positions", ...)
	}

}

#' @method summary spclust
#' @S3method summary spclust

### change this if necessary. 
summary.spclust <-
function(object, ...)
{
  cat("-------------------------------------------------------\n")
  cat("Summary of spclust object\n")
  cat("-------------------------------------------------------\n")
  cat("Type of experimental cross: ", attr(object, "type"), "\n")
  cat("Original population size: ", nrow(object$clusters), "\n")
  cat("Total sample size: ", length(object$lines), "\n")
  cat("Number of input lines: ", length(object$inputlines), "\n")
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

