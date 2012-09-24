\name{NAM2cross}
\alias{NAM2cross}
\title{Convert a NAM mpcross object to a set of cross objects}
\usage{
  NAM2cross(mpcross)
}
\arguments{
  \item{mpcross}{NAM mpcross object}
}
\value{
  A list of cross objects - one for each family in the NAM.
  Potentially each cross has a different number of lines
  and/or markers.
}
\description{
  Takes a Nested Association Mapping population formatted
  as an mpcross object and converts it to a list of RIL
  crosses with (potentially) different sets of segregating
  markers
}

