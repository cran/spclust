#' Hierarchical clustering *in colour*
#' 
#' Modification of plclust for plotting hclust objects in colour
#' @export
#' @param hclust hclust object
#' @param lab Character vector of labels of leaves of the tree
#' @param lab.col Colour for the labels; NA gives the defaul device foreground colour
#' @param hang Fraction of the plot height by which labels hang belore the rest of the plot, as in plclust
#' @param ... Additional arguments passed through to plclust
#' @return Display of dendrogram with coloured leaf labels
#' @seealso \code{\link[spclust]{plot.spclust}}
#' @author originally written by Eva KF Chan, modification by Emma Huang
#' @examples
#' map <- sim.map(len=rep(100, 5), n.mar=21, eq.spacing=TRUE, include.x=FALSE)
#' dat <- sim.cross(map, n.ind=500, type="f2")
#' distmat <- spdist(dat, type="f2")
#' hc <- hclust(distmat, method="average")
#' lab <- rep("", 500)
#' lab[seq(1, 500, 50)] <- "RED"
#' lab[seq(10, 500, 50)] <- "BLACK"
#' lab.col <- rep("red", 500)
#' lab.col[seq(1,500, 50)] <- "blue"
#' plclust_in_colour(hc, lab, lab.col)

plclust_in_colour <- function( hclust, lab=hclust$labels, lab.col=rep(1,length(hclust$labels)), hang=0.1, ... )
{
## modification of plclust for plotting hclust objects *in colour*!
## Copyright Eva KF Chan 2009
## Arguments:
##hclust:hclust object
##lab:a character vector of labels of the leaves of the tree
##lab.col:colour for the labels; NA=default device foreground colour
##hang: as in hclust & plclust
## Side effect:
##A display of hierarchical cluster with coloured leaf labels. 
     y <- rep(hclust$height,2)
     x <- as.numeric(hclust$merge)
     y <- y[which(x<0)]
     x <- x[which(x<0)]
     x <- abs(x)
     y <- y[order(x)]
     x <- x[order(x)]
     plot( hclust, labels=F, hang=hang, ... )
text( x=x+.5, y=y[hclust$order]-diff(range(hclust$height))*hang, labels=lab[hclust$order], col=lab.col[hclust$order], srt=90, adj=c(1,0.5), xpd=NA, ... )
}


