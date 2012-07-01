\name{hzar.next.fitRequest}
\alias{hzar.next.fitRequest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate a new fitRequest using data from another fitRequest.
%%  ~~function to do ... ~~
}
\description{
  The method is the glue for parallel runs of sequential chains.
  It returns a ready to run \code{\link{hzar.fitRequest}} object based
  on the \code{\link{hzar.fitRequest}} supplied.  If
  \code{oldFitRequest} had already been successfully run, this method's
  result will be dependent on the prior run.  If not, this method's
  result will be a request for an independent run.
}
\usage{
hzar.next.fitRequest(oldFitRequest)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{oldFitRequest}{
    A \code{\link{hzar.fitRequest}} object.
  }
}
\details{
  If \code{oldFitRequest} describes a successful run, a new covariance
  matrix is generated, modelParam$init is updated to the covariance
  matrix center, and the lecuyer seed channel incremented by 1. If the
  mersenne twister was used previously, the lecuyer random number
  generator is requested on channel 2.

  If \code{oldFitRequest} does not describe a successful run, everything
  is copied, except the lecuyer seed channel incremented by 10 to
  prevent overlapping. If the mersenne twister was used previously, the
  lecuyer random number generator is requested on channel 11. 
  
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A \code{\link{hzar.fitRequest}} object.
}
% \references{
% %% ~put references to the literature/web site here ~
% }
\author{
  Graham Derryberry \email{asterion@alum.mit.edu}
}
% \note{
%   The comment about overlapping seed channels is a result of the
%   following scenario:

%   \code{ 
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{hzar.fitRequest}}
  \code{\link{hzar.chain.doSeq}}
  \code{\link{hzar.cov.mcmc}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(manakinMolecular);
mknAdaA <-
  hzar.doMolecularData1DPops(manakinMolecular$distance,
                             manakinMolecular$ada.A,
                             manakinMolecular$ada.nSamples);
mknAdaAmodel <-
  hzar.makeCline1DFreq(mknAdaA, scaling="fixed",tails="none");

mknAdaAmodelFit <-
   hzar.first.fitRequest.old.ML(model=mknAdaAmodel ,
                                mknAdaA,
                                verbose=FALSE);

mknAdaAmodelFit$mcmcParam$chainLength <- 1e4;
str(hzar.next.fitRequest(mknAdaAmodelFit))
\dontrun{
mknAdaAinitialFit <- hzar.doFit(mknAdaAmodelFit);
str(hzar.next.fitRequest(mknAdaAinitialFit))
}
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
