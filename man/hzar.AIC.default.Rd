\name{hzar.AIC.default}
\alias{hzar.AIC.default}
\alias{hzar.AICc.default}
\alias{hzar.AIC.hzar.cline}
\alias{hzar.AICc.hzar.cline}
\alias{hzar.AIC.hzar.dataGroup}
\alias{hzar.AICc.hzar.dataGroup}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the AIC score.
}
\description{
  Calculate the AIC or the corrected AIC (AICc) for the given
  likelihood, number of parameters and number of observations.

  Extracts the parameters as needed when passed the correct hzar
  object. 
}
\usage{
hzar.AIC.default(maxLL, param.count)
hzar.AICc.default(maxLL, param.count, nObs)
hzar.AIC.hzar.cline(cline)
hzar.AICc.hzar.cline(cline,nObs)
hzar.AIC.hzar.dataGroup(dataGroup)
hzar.AICc.hzar.dataGroup(dataGroup)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{maxLL}{
    The maximum log likelihood value.
  }
  \item{param.count}{
    The number of free parameters, also known as the number of degrees of freedom.
  }
  \item{nObs}{
    The number of samples observed.
  }
  \item{cline}{
    A \code{\link{hzar.cline}} object.
  }
  \item{dataGroup}{
    A \code{\link{hzar.dataGroup}} object.
  }
}

\details{
  The formula for AIC used is 2 * (param.count - maxLL).

  The formula for AICc used is:
  AIC
  + 2 * param.count * (param.count + 1)
    / (nObs - param.count - 1)
}
\value{
  The AIC or AICc score calculated.
}

% \references{
% %% ~put references to the literature/web site here ~
% }
\author{
  Graham Derryberry \email{asterion@alum.mit.edu}
}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{AIC}}
  \code{\link{hzar.AIC.hzar.obsDataGroup}}
}
\examples{

print(hzar.AIC.default(-8,3))
print(hzar.AICc.default(-8,3,30))

data(manakinMolecular);
mknAdaA <-
  hzar.doMolecularData1DPops(manakinMolecular$distance,
                             manakinMolecular$ada.A,
                             manakinMolecular$ada.nSamples);
hzar.plot.obsData(mknAdaA);
mknAdaAmodel <-
  hzar.makeCline1DFreq(mknAdaA, scaling="fixed",tails="none");
mknAdaAmodel <-
  hzar.model.addBoxReq(mknAdaAmodel,-30,600);
mknAdaAmodelFitR <-
   hzar.first.fitRequest.old.ML(model=mknAdaAmodel ,
                                mknAdaA,
                                verbose=FALSE);
print(hzar.AIC.hzar.dataGroup(hzar.fit2DataGroup(mknAdaAmodelFitR)))

mknAdaAcline <- hzar.gen.cline(list(center=300,width=10),
                               mknAdaAmodelFitR);

print(hzar.AIC.hzar.cline(mknAdaAcline));
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
