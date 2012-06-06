#' Perform multi-stage selective phenotyping clustering
#'
#' This function performs a second stage of selective phenotyping clustering
#' where the second stage may aim to maximize diversity in regions of the 
#' genome highlighted in results from the first stage of phenotyping. The plot
#' function plots the clusters from spclust2 with some summary information.
#' For hierarchical clustering methods
#' the dendrogram is displayed with clusters and individuals selected in 
#' each stage are marked. Graphical genotypes are displayed for individuals
#' selected with maximal recombinations. 
#'
#' @useDynLib spclust
#' @rdname spclust2-all
#' @aliases spclust2 plot.spclust2 summary.spclust2 print.spclust2
#' @param object Cross or mpcross object containing genetic data; for summary function, spclust object
#' @param x spclust object input to plot function
#' @param spobj spclust object with output from the first stage of SPCLUST
#' @param ss2 Sample size to be selected in second stage
#' @param file Optional argument, filename for outputting clusters to a file
#' @param method Selection method in second stage - may either match that in first stage or use maximum recombinations 
#' @param step Step size used in estimating recombinations (default=5 cM)
#' @param threshold Threshold used in estimating recombinations (default=0.7)
#' @param ... Additional arguments to be passed on to plot functions
#' @usage spclust2 (object, spobj, ss2, file, method=c("stage1", "maxrec"), step=5, threshold=.7) 
#' \method{plot}{spclust2} (x, ...) 
#' \method{print}{spclust2} (x, ...)
#' \method{summary}{spclust2} (object, ...)
#' @export spclust2
#' @return list with components:
#' \item{stage1}{Input spclust object containing stage 1 selections}
#' \item{stage2}{spclust object containing stage 2 selections}
#' \item{lines}{indices of selected lines from both stages}
#' \item{mind}{For each selected line, minimum distance to other lines in sample}
#' \item{tree}{Hierarchical clustering tree}
#' \item{clusters}{Assignment of all lines to clusters}
#' @seealso \code{\link[spclust]{spdist}}, \code{\link[spclust]{plot.spclust}}, \code{\link[spclust]{spclust}}, \code{\link[spclust]{plclust_in_colour}}
#' @note Stage 2 SPCLUST performs the following steps in order to select a
#' second stage of selected lines with high genetic diversity after accounting
#' for the originally selected lines. First, genetic distances are estimated
#' between all lines in the sample, based on the expected proportion of 
#' alleles not shared IBD across the genome. Second, if the "maxrec" method
#' is selected, ss2 lines are selected with the highest number of estimated
#' recombinations, excluding those which have already been selected in the 
#' first  stage. Otherwise, all lines are clustered based on the
#'  genetic distance, with the number of clusters equal to
#'   the sum of sample sizes indicated for selection in stages
#'  1 and 2. All lines selected in stage 1 are included in
#'  the final sample, and clusters containing these lines are
#'   excluded from further selection. From ss2 of the
#'  remaining clusters, a representative is selected as the
#'  one most similar to other lines in the cluster.
#' @examples
#' # Simulate a map and data using qtl package
#' map <- sim.map(len=rep(100, 5), n.mar=21, eq.spacing=TRUE, include.x=FALSE)
#' dat <- sim.cross(map, n.ind=500, type="f2")
#' # Select two samples of size 100 in two stages
#' sp <- spclust(dat, 100, method="ward")
#' sp2 <- spclust2(dat, sp, 100, method="stage1")
#' summary(sp2)
#' plot(sp2)

spclust2 <- function(object, spobj, ss2, file, method=c("stage1", "maxrec"), step=5, threshold=.7)
{
  if (missing(object)) stop("Cannot select lines without input genetic data\n")
      
  output <- list()
  output$stage1 <- spobj
  output$stage2 <- list()

  if (method=="stage1") attr(output$stage2, "method") <- attr(output$stage1, "method") else attr(output$stage2, "method") <- method

  if (inherits(object, "mpcross")) type <- "magic"
  else if (inherits(object, "cross")) type <- class(object)[1]

  if (!(type %in% c("bc", "f2", "magic"))) stop("Invalid type of experimental cross has been input.\n")

  if (method=="maxrec" & type=="f2") stop("Cannot currently use this method with this type of cross\n")

  if (inherits(object, "cross"))  nmrk <- sum(nmar(object)) else nmrk <- ncol(object$finals)

  attr(output$stage2, "nmrk") <- nmrk

  if (missing(spobj)) {
        cat("No lines from stage 1 input; will only select stage 2 sample\n")
	cat("Method used will be Ward clustering\n")
	attr(output$stage1, "method") <- "ward"
        ss1 <- 0
  } else ss1 <- length(spobj$lines)

  if (missing(ss2)) {
        cat("Will select equal number of lines to stage 1\n")
        ss2 <- ss1
  }

  distmat <- spdist(object, type)
  linenames <- attr(distmat, "linenames")
  d <- as.matrix(distmat)

  if (method=="maxrec") {
	## need to separate out magic from other types
	if (type=="magic") {
  	  if (!inherits(object, "mpprob")) {
	  mpp <- mpprob(object, step=step, program="qtl", threshold=threshold)
	  cat("mpprob object not input; computing multi-point probabilities with input step and threshold values\n")
	  } else mpp <- object
	gen <- do.call("cbind", lapply(mpp$estfnd, function(x) cbind(x, rep(0, nrow(x)))))
	nrec <- lapply(mpp$estfnd, function(x) apply(x, 1, function(y) return(sum(diff(y[!is.na(y)])!=0))))
	}
	
	if (type=="bc") {
	   gen <- lapply(object$geno, function(x) x$data)
	   nrec <- lapply(gen, function(x) apply(x, 1, function(y) return(sum(diff(y[!is.na(y)])!=0))))
	   gen <- do.call("cbind", lapply(gen, function(x) cbind(x, rep(0, nrow(x)))))
	}

  	totnrec <- apply(do.call("cbind", nrec), 1, sum)
	ord <- data.frame(LineName=linenames, No.Rec=totnrec) 
  	ord <- ord[order(ord[,2], decreasing=TRUE),]

  	keep <- match(spobj$lines, linenames)
  	start <- 1

   	### need to select most recombinant lines
  	while (length(keep) < ss1+ss2) {
    	  rem <- ss1+ss2-length(keep)
     	  keep <- unique(c(keep, ord[start:(start+rem-1),1]))
     	  start <- start+rem
  	}
	allkeep <- keep
	keep <- setdiff(keep, match(spobj$lines, linenames))
	sel <- rep(0, nrow(ord))
	sel[match(spobj$lines, ord[,1])] <- 1
	sel[setdiff(keep, match(spobj$lines, linenames))] <- 2
	ord$StageSelected <- sel
	output$result <- ord
	output$stage2$recmat <- gen
	attr(output$stage2, "step") <- step
	attr(output$stage2, "threshold") <- threshold
	if (type=="bc") n.founders <- 2 else n.founders <- nrow(mpp$founders)
	attr(output$stage2, "n.founders") <- n.founders 
  }

  ## need to modify this so that it works properly when stage1 method = pam
  if (method=="stage1") {

    if (attr(output$stage1, "method") != "pam") {
      hc <- hclust(distmat, method=attr(output$stage1, "method"))
      cl <- cutree(hc, k=ss1+ss2)
      hc$clustering <- cl
      reps <- unlist(by(1:nrow(d), cl, function(x) if (length(x)==1) return(x) else x[which.min(apply(d[x,x], 2, sum))]))
     ### count how many clusters do not contain lines from stage 1
     clust <- vector(length=ss1+ss2)
     if (ss1>0)
     for (i in 1:(ss1+ss2)) clust[i] <- length(intersect(match(spobj$lines, linenames), which(cl == i)))
     ncl <- sum(clust==0)
     ### sample from clusters which do not contain lines from stage 1
     clust2 <- sample(which(clust==0), ss2)
     keep <- reps[clust2]

     allkeep <- c(keep, match(output$stage1$lines, linenames))

     sel <- rep(0, nrow(d))
     sel[keep] <- 2
     sel[match(output$stage1$lines, linenames)] <- 1
     clall <- data.frame(LineName=linenames, Cluster=cl, StageSelected=sel)
     clall <- clall[order(clall[,2]),]

     cl2 <- clall[clall$Cluster%in%clust2,]

     ## can also produce a clustering matrix for stage 2 lines
     
     output$stage2$clusters <- cl2
     output$result <- hc
     output$clusters <- clall
     } else cat("Cannot currently accommodate pam clustering for stage 2\n\n")
  }

  dkeep <- d[keep, keep]
  diag(dkeep) <- NA
  mind <- apply(dkeep, 2, function(x) min(x, na.rm=TRUE))

  dallkeep <- d[allkeep, allkeep]
  diag(dallkeep) <- NA
  mindall <- apply(dallkeep, 2, function(x) min(x, na.rm=TRUE))

  output$stage2$mind <- data.frame(lines=linenames[keep], min.dist=mind)
  output$stage2$lines <- linenames[keep]
  output$mind <- data.frame(lines=linenames[allkeep], min.dist=mindall) 
  output$lines <- linenames[allkeep]

  if (!missing(file)) write.csv(cbind(1:(ss1+ss2), output$lines), file=file, row.names=FALSE)

  class(output) <- "spclust2"
  output
}

#' @method plot spclust2
#' @S3method plot spclust2

plot.spclust2 <- function(x, ...)
{
	pp <- par(ask=TRUE)
	on.exit(par(pp))

	## Alter this so that stage1 and stage2 colours are different
	if (attr(x$stage2, "method") %in% c("ward", "average")) {
		lab <- rep("", nrow(x$clusters))
		lab[match(x$stage1$lines, x$clusters[,1])] <- x$stage1$lines
		lab[match(x$stage2$lines, x$clusters[,1])] <- x$stage2$lines
		lab.col <- rep("blue", nrow(x$clusters))
		lab.col[match(x$stage1$lines, x$clusters[,1])] <- "black"
		lab.col[match(x$stage2$lines, x$clusters[,1])] <- "red"
		plclust_in_colour(x$result, lab=lab, lab.col=lab.col, ...)
		rect.hclust(x$result, k=length(x$lines), border="blue")
	}
	## silhouette plot
#	sil <- silhouette(x$result$clustering, dmatrix=x$dmatrix)
	## modify sil to remove clusters with only one individual
#	a <- which(sil[,1] %in% names(table(sil[,1]))[table(sil[,1])>1])
#	sil1 <- silhouette(sil[a,1], x$dmatrix[a,a])
#	clnames <- names(table(sil1[,1]))
#	clcol <- rep("firebrick", nrow(sil1))
#	clcol[which(sil1[,1] %in% clnames[seq(2,length(clnames), 2)])] <- "seagreen4"
#	plot(sil1, col=clcol, main="Silhouette plot of clusters with size > 1")

	## Recombination diagram
	if (attr(x$stage2, "method")=="maxrec") {
	  n.founders <- attr(x$stage2, "n.founders")
	  mat <- x$stage2$recmat[match(x$stage2$lines, x$result[,1]),]
	  rownames(mat) <- x$stage2$lines
	  heatmap(mat, col=c("lightblue2", "royalblue1", "darkseagreen1", "seagreen3", "lemonchiffon", "darkgoldenrod1", "indianred1", "firebrick")[1:n.founders], main="Recombinations for lines selected in stage 2", Rowv=NA, Colv=NA, scale="none", ylab="Lines", xlab="Positions")
	}

}

#' @method summary spclust2
#' @S3method summary spclust2

summary.spclust2 <-
function(object, ...)
{
  cat("-------------------------------------------------------\n")
  cat("Summary of spclust2 object\n")
  cat("-------------------------------------------------------\n")
  cat("Type of experimental cross: ", attr(object, "type"), "\n")
  cat("Original population size: ", nrow(object$clusters), "\n")
  cat("Total selected sample size: ", length(object$lines), "\n")
  cat("-------------------------------------------------------\n")
  cat("First Stage\n\n")
  cat("Selected sample size: ", length(object$stage1$lines), "\n")
  cat("Number of genetic markers used: ", attr(object$stage1, "nmrk"), "\n")
  cat("Clustering method: ", attr(object$stage1, "method"), "\n")
  cat("Range of minimum distances between selected individuals: ", round(range(object$stage1$mind[,2]),3), "\n")
  cat("-------------------------------------------------------\n")
  cat("Second Stage\n\n")
  cat("Selected sample size: ", length(object$stage2$lines), "\n")
  cat("Number of genetic markers: ", attr(object$stage2, "nmrk"), "\n")
  cat("Clustering method: ", attr(object$stage2, "method"), "\n")
  cat("Range of minimum distances between selected individuals: ", round(range(object$stage2$mind[,2]),3), "\n")
}

#' @method print spclust2
#' @S3method print spclust2

print.spclust2 <-
function(x, ...)
{
  cat(" This is an object of class \"spclust2\".\n")
  cat(" We provide the following summary. \n")
  summary(x)
}

