\name{hzar.dataGroup.null}
\alias{hzar.dataGroup.null}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Datagroup placeholder for the null model (frequency independent of location)
%%  ~~function to do ... ~~
}
\description{
  Generates a \code{\link{hzar.dataGroup}} object representing a fit of
  the null model to a \code{\link{hzar.obsData}} object.
}
\usage{
hzar.dataGroup.null(obsData)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obsData}{
    The \code{\link{hzar.obsData}} object for which to generate a
    \code{\link{hzar.dataGroup}} object. 
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
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{hzar.dataGroup}}
  \code{\link{hzar.obsData}}
  \code{\link{hzar.make.LLfunc.null}}
}
\examples{
data(manakinMolecular);
mknAdaA <-
  hzar.doMolecularData1DPops(manakinMolecular$distance,
                             manakinMolecular$ada.A,
                             manakinMolecular$ada.nSamples);
mkn.AdaA.null <- hzar.dataGroup.null(mknAdaA)
hzar.plot.cline(mkn.AdaA.null);
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
