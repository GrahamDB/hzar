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
hzar.multiFitRequest(fitL, each = 1, baseSeed = c(1234, 2345, 3456,
4567, 5678, 6789, 7890, 8901, 9012, 123), rotateSeed = TRUE, baseChannel
= 50, adjChannel = 50) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitL}{
    A single \code{\link{hzar.fitRequest}} objects or a list of
    \code{\link{hzar.fitRequest}} objects 
  }
  \item{each}{
    How many times to replicate each \code{\link{hzar.fitRequest}}
    object. 
  }
  \item{baseSeed}{
    The pool of values from which to draw seeds.

    If \code{NULL}, do not change the seed.
    
    If \code{rotateSeed} is \code{TRUE} this pool is automatically
    reduced to unique values.  
  }
  \item{rotateSeed}{
    If \code{TRUE}, a unique set of six values is drawn from
    \code{baseSeed}.
    If \code{FALSE}, the first 6 values of \code{baseSeed} are used. The
    method \code{\link{rep}} is used to expand \code{baseSeed} to 6
    values if needed. 
  }
  \item{baseChannel}{
    The initial stream channel to set for each element of \code{fitL}.
    The method \code{\link{rep}} is used to expand \code{baseChannel} to
    the length of \code{fitL} if needed. 
    If \code{NULL}, use the original stream channel from each element of
    \code{fitL}. 
  }
  \item{adjChannel}{
    Amount to increment the stream channel
  }
}
\details{
  This method assumes that the user wishes to generate independent
  chains unless instructed other wise.

  By default, this method will use a unique seed for each element of
  \code{fitL} and increment the stream channel for each replication of
  each element of \code{fitL}.

  If \code{rotateSeed} is \code{TRUE}, \code{baseSeed} is numeric, and
  \code{adjChannel} is not numeric, then this method will use a unique
  seed for every element of the result.

  If \code{rotateSeed} is \code{FALSE}, \code{baseChannel} is a numeric
  of length 1, \code{adjChannel} is numeric, and \code{each} is greater
  than one, then this method will increment the stream channel for every
  element of the result. 
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