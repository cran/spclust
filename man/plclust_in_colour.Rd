\name{plclust_in_colour}
\alias{plclust_in_colour}
\title{Hierarchical clustering *in colour*}
\usage{
  plclust_in_colour(hclust, lab = hclust$labels,
    lab.col = rep(1, length(hclust$labels)), hang = 0.1,
    ...)
}
\arguments{
  \item{hclust}{hclust object}

  \item{lab}{Character vector of labels of leaves of the
  tree}

  \item{lab.col}{Colour for the labels; NA gives the defaul
  device foreground colour}

  \item{hang}{Fraction of the plot height by which labels
  hang belore the rest of the plot, as in plclust}

  \item{...}{Additional arguments passed through to
  plclust}
}
\value{
  Display of dendrogram with coloured leaf labels
}
\description{
  Modification of plclust for plotting hclust objects in
  colour
}
\examples{
map <- sim.map(len=rep(100, 5), n.mar=21, eq.spacing=TRUE, include.x=FALSE)
dat <- sim.cross(map, n.ind=500, type="f2")
distmat <- spdist(dat, type="f2")
hc <- hclust(distmat, method="average")
lab <- rep("", 500)
lab[seq(1, 500, 50)] <- "RED"
lab[seq(10, 500, 50)] <- "BLACK"
lab.col <- rep("red", 500)
lab.col[seq(1,500, 50)] <- "blue"
plclust_in_colour(hc, lab, lab.col)
}
\author{
  originally written by Eva KF Chan, modification by Emma
  Huang
}
\seealso{
  \code{\link[spclust]{plot.spclust}}
}

