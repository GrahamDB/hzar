\name{hzar.multiFitRequest}
\alias{hzar.multiFitRequest}
\title{
  Collect hzar.fitRequest objects to simplify automation 
}
\description{
  This method manipulates a list of \code{hzar.fitRequest} objects to
  ensure each object has an independent seed.  If requested, it can
  replicate each \code{hzar.fitRequest} to create independent chains. 
}
\usage{
hzar.multiFitRequest(fitL, each = 1, baseSeed = c(1234, 2345, 3456, 4567, 5678, 6789, 7890, 8901, 9012), rotateSeed = TRUE, baseChannel = 50, adjChannel = 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitL}{
    A single \code{\link{hzar.fitRequest}} objects or a list of
    \code{\link{hzar.fitRequest}} objects 
  }
  \item{each}{
    %%     ~~Describe \code{\link{each}} here~~
  }
  \item{baseSeed}{
    %%     ~~Describe \code{\link{baseSeed}} here~~
  }
  \item{rotateSeed}{
    %%     ~~Describe \code{\link{rotateSeed}} here~~
  }
  \item{baseChannel}{
    %%     ~~Describe \code{\link{baseChannel}} here~~
  }
  \item{adjChannel}{
    Amount to increment the seed channel
  }
}
\details{
  The behavior of this method is complex.
}
\value{
   Returns a list of \code{\link{hzar.fitRequest}} object, suitable for
   \code{\link{hzar.doFit.multi}}. 

}

\author{
  Graham Derryberry \email{asterion@alum.mit.edu}
}
\seealso{
  \code{\link{hzar.fitRequest}}
  \code{\link{hzar.doFit.multi}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.


}