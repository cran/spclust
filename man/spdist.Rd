\name{spdist}
\alias{spdist}
\title{Compute genetic distances for experimental crosses}
\usage{
  spdist(object, type = c("bc", "f2", "magic", "nam"))
}
\arguments{
  \item{object}{cross or mpcross object containing genetic
  data}

  \item{type}{type of experimental cross}
}
\value{
  distance object for use in spclust
}
\description{
  This function computes the expected percentage of markers
  shared IBD across the genome for experimental crosses of
  type backcross, doubled haploid, recombinant inbred line
  ("bc"), F2 ("f2"), Multiparent Advanced Generation
  Inter-Cross with 4 or 8 parents ("magic"), or Nested
  Association Mapping ("nam").
}
\examples{
# Simulate a map and data using qtl package
map <- sim.map(len=100, n.mar=101, eq.spacing=TRUE, include.x=FALSE)
dat1 <- sim.cross(map, n.ind=200, type="bc")
dat2 <- sim.cross(map, n.ind=500, type="f2")
# Calculate the distances between lines in datasets
dist1 <- spdist(dat1, "bc")
dist2 <- spdist(dat2, "f2")
# Heatmap of distances
heatmap(as.matrix(dist1), Rowv=NA, Colv=NA, scale="none", col=topo.colors(10))
par(mfrow=c(2, 1))
# Histogram of minimum distances between lines for full datasets
d2 <- as.matrix(dist2)
diag(d2) <- NA
mind <- apply(d2, 1, function(x) min(x, na.rm=TRUE))
hist(mind, col="tomato", main="Minimum distances between all F2 lines", xlab="Distance")
# Histogram of minimum distances between selected lines
sp <- spclust(dat2, 100, method="ward")
hist(sp$mind[, 2], col="tomato", main="Minimum distances between 100 selected F2 lines", xlab="Distance")
}
\seealso{
  \code{\link[spclust]{spclust}}
}

