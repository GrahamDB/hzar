\name{hzar.doFit.multi}
\alias{hzar.doFit.multi}
\alias{hzar.doChain.multi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Run hzar fit commands on a list of hzar.fitRequest objects
}
\description{
  These methods simplify repeated calling of \code{\link{hzar.doFit}} or
  \code{\link{hzar.chain.doSeq}} while taking advantage of
  \code{\link{\%dopar\%}} if requested. 
}
\usage{
hzar.doFit.multi(mFitR, doPar = FALSE, inOrder = TRUE)
hzar.doChain.multi(mFitR, doPar = FALSE, inOrder = TRUE, ...)
}
\arguments{
  \item{mFitR}{
    Provide a list of \code{\link{hzar.fitRequest}} objects. Use
    \code{\link{hzar.multiFitRequest}} to ensure independent seeds and
    to request independent chains.
  }
  \item{doPar}{
    Use \code{\link{\%dopar\%}}?
  }
  \item{inOrder}{
    Should the results be returned in order? If \code{FALSE}, the
    results are returned in the order of completion. See
    \code{\link{foreach}} for more information.
  }
  \item{\dots}{
    Additional arguments to pass to \code{\link{hzar.chain.doSeq}}
  }
}
\value{
  A list of the fitted \code{\link{hzar.fitRequest}} objects.
}
\author{
  Graham Derryberry \email{asterion@alum.mit.edu}
}
\seealso{
  \code{\link{hzar.chain.doSeq}}
  \code{\link{hzar.doFit}}
  \code{\link{hzar.fitRequest}}
  \code{\link{hzar.multiFitRequest}} 
  \code{\link{\%dopar\%}}
  \code{\link{foreach}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
