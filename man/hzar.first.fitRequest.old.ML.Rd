\name{hzar.first.fitRequest.old.ML}
\alias{hzar.first.fitRequest.old.ML}
\title{
  Generate a ML based hzar.fitRequest using a meta model structure.
}
\description{
  This method generates a \code{\link{hzar.fitRequest}} object suitable
  for \code{\link{hzar.doFit}}.
}
\usage{
hzar.first.fitRequest.old.ML(model, obsData, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
    A \code{\link{clineMetaModel}} object. 
  }
  \item{obsData}{
    The \code{\link{hzar.obsData}} object to which the meta model is to be fit.
  }
  \item{verbose}{
    Should \code{\link{MCMCmetrop1R}} be verbose?
  }
}

% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
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
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{hzar.fitRequest}}
  \code{\link{hzar.obsData}}
  \code{\link{clineMetaModel}}
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
mknAdaAmodelFit$mcmcParam$chainLength <- 1e5;
str(mknAdaAmodelFit);
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
