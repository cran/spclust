\name{generateNAMpedigree}
\alias{generateNAMpedigree}
\title{Generates NAM pedigree}
\usage{
  generateNAMpedigree(nparents = 26, nperfam = 200,
    nssdgen = 6)
}
\arguments{
  \item{nparents}{Total number of parents}

  \item{nperfam}{Number of lines per family}

  \item{nssdgen}{Number of generations of selfing}
}
\value{
  A four column matrix with IDs for each individual,
  maternal parent ID, paternal parent ID, and a flag for
  whether the individual was observed.
}
\description{
  Generated a pedigree for a Nested Association Mapping
  population
}
\note{
  The first founder will be the common parent, assumed for
  all F1s to be the maternal parent. Default parameters are
  those for the original maize NAM population (ref).
}

