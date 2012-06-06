\name{spclust}
\alias{plot.spclust}
\alias{print.spclust}
\alias{spclust}
\alias{summary.spclust}
\title{Perform single-stage selective phenotyping clustering}
\usage{
spclust (object, nlines, file, method=c("average", "ward", "pam")) 
\method{plot}{spclust} (x, ...)
\method{print}{spclust} (x, ...)
\method{summary}{spclust} (object, ...)
}
\arguments{
  \item{object}{Cross or mpcross object containing genetic
  data; for summary, spclust object}

  \item{x}{spclust object to plot}

  \item{nlines}{Size of selected subsample}

  \item{file}{Optional argument, filename for outputting
  clusters to a file}

  \item{method}{Clustering method}

  \item{...}{Additional argument passed through to plclust}
}
\value{
  list with components: \item{lines}{indices of selected
  lines} \item{mind}{For each selected line, minimum
  distance to other lines in sample}
  \item{tree}{Hierarchical clustering tree}
  \item{clusters}{Assignment of all lines to clusters}
}
\description{
  This function implements the SPCLUST algorithm to perform
  selective phenotyping in experimental crosses by
  maximizing the genetic diversity in the selected
  subsample. The plot function plots the clusters from
  spclust with some summary information. For all clustering
  methods the silhouette plot is displayed for clusters
  containing more than one individual. For hierarchical
  clustering methods the dendrogram is displayed with
  clusters and selected individuals marked.
}
\note{
  The SPCLUST algorithm performs the following steps in
  order to select a subsample with high genetic diversity.
  First, genetic distances are estimated between all lines
  in the sample, based on the expected proportion of
  alleles not shared IBD across the genome. Second, lines
  are clustered based on the genetic distance, with the
  number of clusters matching the number of lines desired
  for selection. Third, a representative line is selected
  from each cluster as the one most similar to other lines
  in the cluster.
}
\examples{
# Simulate a map and data using the qtl package
map <- sim.map(len=rep(100, 5), n.mar=21, eq.spacing=TRUE, include.x=FALSE)
dat <- sim.cross(map, n.ind=500, type="f2")
# Selection of 100 lines
sp <- spclust(dat, 100, method="ward")
summary(sp)
plot(sp)
}
\seealso{
  \code{\link[spclust]{spdist}},
  \code{\link[cluster]{plot.silhouette}},
  \code{\link[spclust]{plclust_in_colour}}
}

