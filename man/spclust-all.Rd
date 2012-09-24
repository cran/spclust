\name{spclust}
\alias{plot.spclust}
\alias{print.spclust}
\alias{spclust}
\alias{summary.spclust}
\title{Perform selective phenotyping clustering}
\usage{
  spclust (object, nlines, method=c("average", "ward",
    "pam", "maxrec"), inputlines=NULL, file, step=5,
    threshold=.7) 
  \method{plot}{spclust} (x, type=4, ...)
  \method{print}{spclust} (x, ...)
  \method{summary}{spclust} (object, ...)
}
\arguments{
  \item{object}{Cross or mpcross object containing genetic
  data; for summary function, spclust object}

  \item{x}{spclust object input to plot function}

  \item{nlines}{Number of lines to be selected (note: does
  not include the number of input lines)}

  \item{method}{Selection method - options include
  hierarchical clustering (average, ward), partitioning
  around medoids (pam), or based on the maximal number of
  recombinations (maxrec)}

  \item{inputlines}{Names of lines which must be included
  in the selected sample. See details below.}

  \item{file}{Optional argument, filename for outputting
  clusters to a file}

  \item{step}{Step size used in estimating recombinations
  (default=5 cM)}

  \item{threshold}{Threshold used in estimating
  recombinations (default=0.7)}

  \item{type}{Style of plot to draw; 1=Silhouette;
  2=Dendrogram; 3=Recombinations; 4=All (that are
  appropriate)}

  \item{...}{Additional arguments to be passed on to plot
  functions}
}
\value{
  list with components: \item{numlines}{indices of selected
  lines from both stages} \item{lines}{names of selected
  lines from both stages} \item{mind}{For each selected
  line, minimum distance to other lines in sample}
  \item{tree}{Hierarchical clustering tree}
  \item{clusters}{Assignment of all lines to clusters}
  \item{recmat}{If method="maxrec", returns matrix of
  recombinations for genomic region}
}
\description{
  This function implements the SPCLUST algorithm to perform
  selective phenotyping in experimental crosses by
  maximizing the genetic diversity in the selected
  subsample. Selection can be done in one or multiple
  stages. The plot function plots the clusters with some
  summary information. Graphical genotypes are displayed
  for individuals selected with maximal recombinations. For
  hierarchical clustering methods the dendrogram is
  displayed with clusters and selected individuals marked.
}
\note{
  This function can perform both single-stage or
  multi-stage selective phenotyping clustering. In a single
  stage, the SPCLUST algorithm performs the following steps
  in order to select a subsample with high genetic
  diversity. First, genetic distances are estimated between
  all lines in the sample, based on the expected proportion
  of alleles not shared IBD across the genome. Second,
  lines are clustered based on the genetic distance, with
  the number of clusters matching the number of lines
  desired for selection. Third, a representative line is
  selected from each cluster as the one most similar to
  other lines in the cluster.

  If the inputlines argument is used, SPCLUST performs the
  following steps in order to select a sample with high
  genetic diversity while accounting for the input lines.
  First, genetic distances are estimated between all lines
  in the sample, based on the expected proportion of
  alleles not shared IBD across the genome. Second, if the
  "maxrec" method is selected, nlines lines are selected
  with the highest number of estimated recombinations,
  excluding those which have already been selected in the
  first stage. Otherwise, all lines are clustered based on
  the genetic distance, with the number of clusters equal
  to the sum of the number of input lines and the nlines
  argument. All input lines selected in stage 1 are
  included in the final sample, and clusters containing
  these lines are excluded from further selection. From
  nlines of the remaining clusters, a representative is
  selected as the one most similar to other lines in the
  cluster.
}
\examples{
# Simulate a map and data using qtl package
map <- sim.map(len=rep(100, 5), n.mar=21, eq.spacing=TRUE, include.x=FALSE)
dat <- sim.cross(map, n.ind=500, type="bc")
# Select two samples of size 100 in two stages
sp <- spclust(dat, 100, method="ward")
sp2 <- spclust(dat, 100, method="maxrec", inputlines=sp$lines)
summary(sp2)
plot(sp2)
}
\seealso{
  \code{\link[spclust]{spdist}},
  \code{\link[spclust]{plot.spclust}},
  \code{\link[spclust]{spclust}},
  \code{\link[spclust]{plclust_in_colour}}
}

