\name{hzar.mcmc.bindLL}
\alias{hzar.mcmc.bindLL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate a mcmc object with sampled parameters and log likelihoods.
}
\description{
  This function returns the mcmc data with an added a log likelihood
  column. 
}
\usage{
hzar.mcmc.bindLL(fitRequest,
                 dataGroup = hzar.fit2DataGroup(fitRequest),
                 mcmcData =
                   ifelse(inherits(fitRequest,"hzar.fitRequest"),
                          fitRequest$mcmcRaw,
                          as.mcmc(dataGroup$data.mcmc)),
                 llData = dataGroup$data.LL,
                 t0 = start(mcmcData),
                 tF = thin(mcmcData))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitRequest}{
    The \code{\link{hzar.fitRequest}} or \code{\link{hzar.dataGroup}} object to use. 
  }
  \item{dataGroup}{
    The \code{\link{hzar.dataGroup}} object to use.
  }
  \item{mcmcData}{
    The \code{\link{mcmc}} object with the series parameter values. 
  }
  \item{llData}{
    The series of log likelihoods of the parameter values.
  }
  \item{t0}{
    The \code{\link{start.mcmc}} attribute of the result.
  }
  \item{tF}{
    The \code{\link{thin.mcmc}} attribute of the result.
  }
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  A \code{\link{mcmc}} object, with columns for each free parameter and
  the log likelihood of each row. 
}

\author{
Graham Derryberry \email{asterion@alum.mit.edu}

}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{hzar.fitRequest}}
  \code{\link{hzar.dataGroup}}
  \code{\link{mcmc}}
  
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
plot(hzar.mcmc.bindLL(hzar.doFit(mknAdaAmodelFit)));
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
