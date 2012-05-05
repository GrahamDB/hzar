\name{hzar.doFit}
\alias{hzar.doFit}
\title{
  Run the optimizer.
}
\description{
 Run the optimizer using the parameters listed in the
 \code{\link{hzar.fitRequest}} given.
}
\usage{
hzar.doFit(fitRequest)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitRequest}{
    The \code{\link{hzar.fitRequest}} object to be processed
  }
}
\value{
  An updated \code{\link{hzar.fitRequest}} object.
}

\author{
  Graham Derryberry \email{asterion@alum.mit.edu}
}

\seealso{
  \code{\link{hzar.fitRequest}}
  \code{\link{MCMCmetrop1R}} 
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
mknAdaAmodelFit$mcmcParam$chainLength <- 1e5;
str(hzar.doFit(mknAdaAmodelFit));
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
