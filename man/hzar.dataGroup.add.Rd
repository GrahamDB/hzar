\name{hzar.dataGroup.add}
\alias{hzar.dataGroup.add}
\alias{hzar.fit2DataGroup}
\title{
  Prepare optimizer output for analysis.
}
\description{
  Intended to group multiple fits of the same model and the same
  observation data into a single object.  As it works with
  \code{\link{hzar.fitRequest}} objects, it is simpler to use than
  \code{\link{hzar.make.dataGroup}}. 
}
\usage{
hzar.dataGroup.add(dataGroup, fitRequestL = list(), doPar = FALSE)
hzar.fit2DataGroup(fitRequest, doPar = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataGroup}{
    A single \code{\link{hzar.dataGroup}} object to update.  If
    \code{fitRequestL} is a list of length 0, this argument may also be
    a list of \code{\link{hzar.fitRequest}} or
    \code{\link{hzar.dataGroup}} objects.
  }
  \item{fitRequestL}{
    A \code{\link{hzar.fitRequest}} or \code{\link{hzar.dataGroup}}
    object, a list of such objects, or a deep list of such objects. 
  }
  \item{fitRequest}{
    A single \code{\link{hzar.fitRequest}} object. A
    \code{\link{hzar.dataGroup}} object may also be used.
  }
  \item{doPar}{
    This is argument is passed to \code{\link{hzar.eval.clineLL}}.
  }
}
% \details{
  
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  A \code{\link{hzar.dataGroup}} object.
}
% \references{
% %% ~put references to the literature/web site here ~
% }
\author{
Graham Derryberry \email{asterion@alum.mit.edu}
}
\note{
   A deep list of \code{T} is a list of length greater > 1 that contains
   only deep lists of \code{T} or objects of class \code{T}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{hzar.fitRequest}}
  \code{\link{hzar.dataGroup}}
  \code{\link{hzar.obsDataGroup}}
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
mknAdaAFit <- hzar.doFit(mknAdaAmodelFit);
str(hzar.fit2DataGroup(mknAdaAFit));
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
