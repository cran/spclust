\name{generateNAMfounders}
\alias{generateNAMfounders}
\title{Generate NAM founder genotypes}
\usage{
  generateNAMfounders(nparents, nmarkers, map, option = 1)
}
\arguments{
  \item{nparents}{Total number of parents}

  \item{nmarkers}{Number of markers on all chromosomes}

  \item{map}{Marker map}

  \item{option}{Process for generating genotypes. See
  details below.}
}
\value{
  A matrix with nparents rows and nmarkers columns, with
  values in {0, 1}. The attribute "MAF" will give a vector
  containing minor allele frequencies for each marker (set
  at 1/nparents if option 1 is selected).
}
\description{
  Generates founder genotypes for a Nested Association
  Mapping population
}
\note{
  If option==1 (default), the common parent (1) is assumed
  to have a different value from all the other parents.
  This would be achieved if markers were chosen for which
  that parent had a rare allele. If option==2, then the
  minor allele frequency for each markers is uniformly
  distributed on [0, .5] and parents are randomly assigned
  one of two alleles (0 or 1) according to the allele
  frequency. Note that in all cases "1" is the minor
  allele.  The marker map is currently only used for the
  marker names, but alternative methods of generating
  founder genotypes could depend on it. If it is not input
  a default marker naming will be used.
}

