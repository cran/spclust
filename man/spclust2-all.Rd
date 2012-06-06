\name{spclust2}
\alias{plot.spclust2}
\alias{print.spclust2}
\alias{spclust2}
\alias{summary.spclust2}
\title{Perform multi-stage selective phenotyping clustering}
\usage{
spclust2 (object, spobj, ss2, file, method=c("stage1", "maxrec"), step=5, threshold=.7)
\method{plot}{spclust2} (x, ...)
\method{print}{spclust2} (x, ...)
\method{summary}{spclust2} (object, ...)
}
\arguments{
  \item{object}{Cross or mpcross object containing genetic
  data; for summary function, spclust object}

  \item{x}{spclust object input to plot function}

  \item{spobj}{spclust object with output from the first
  stage of SPCLUST}

  \item{ss2}{Sample size to be selected in second stage}

  \item{file}{Optional argument, filename for outputting
  clusters to a file}

  \item{method}{Selection method in second stage - may
  either match that in first stage or use maximum
  recombinations}

  \item{step}{Step size used in estimating recombinations
  (default=5 cM)}

  \item{threshold}{Threshold used in estimating
  recombinations (default=0.7)}

  \item{...}{Additional arguments to be passed on to plot
  functions}
}
\value{
  list with components: \item{stage1}{Input spclust object
  containing stage 1 selections} \item{stage2}{spclust
  object containing stage 2 selections}
  \item{lines}{indices of selected lines from both stages}
  \item{mind}{For each selected line, minimum distance to
  other lines in sample} \item{tree}{Hierarchical
  clustering tree} \item{clusters}{Assignment of all lines
  to clusters}
}
\description{
  This function performs a second stage of selective
  phenotyping clustering where the second stage may aim to
  maximize diversity in regions of the genome highlighted
  in results from the first stage of phenotyping. The plot
  function plots the clusters from spclust2 with some
  summary information. For hierarchical clustering methods
  the dendrogram is displayed with clusters and individuals
  selected in each stage are marked. Graphical genotypes
  are displayed for individuals selected with maximal
  recombinations.
}
\note{
  Stage 2 SPCLUST performs the following steps in order to
  select a second stage of selected lines with high genetic
  diversity after accounting for the originally selected
  lines. First, genetic distances are estimated between all
  lines in the sample, based on the expected proportion of
  alleles not shared IBD across the genome. Second, if the
  "maxrec" method is selected, ss2 lines are selected with
  the highest number of estimated recombinations, excluding
  those which have already been selected in the first
  stage. Otherwise, all lines are clustered based on the
  genetic distance, with the number of clusters equal to
  the sum of sample sizes indicated for selection in stages
  1 and 2. All lines selected in stage 1 are included in
  the final sample, and clusters containing these lines are
  excluded from further selection. From ss2 of the
  remaining clusters, a representative is selected as the
  one most similar to other lines in the cluster.
}
\examples{
# Simulate a map and data using qtl package
map <- sim.map(len=rep(100, 5), n.mar=21, eq.spacing=TRUE, include.x=FALSE)
dat <- sim.cross(map, n.ind=500, type="f2")
# Select two samples of size 100 in two stages
sp <- spclust(dat, 100, method="ward")
sp2 <- spclust2(dat, sp, 100, method="stage1")
summary(sp2)
plot(sp2)
}
\seealso{
  \code{\link[spclust]{spdist}},
  \code{\link[spclust]{plot.spclust}},
  \code{\link[spclust]{spclust}},
  \code{\link[spclust]{plclust_in_colour}}
}

