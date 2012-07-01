\name{hzar.sameModel}
\alias{hzar.sameModel}
\alias{hzar.sameObsData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Do the hzar objects share the same model or data?
}
\description{
  Test hzar objects for identical associated cline models or identical
  hzar.obsData objects.
}
\usage{
hzar.sameModel(fitA, fitB)
hzar.sameObsData(fitA, fitB)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitA}{
    An object to compare.
  }
  \item{fitB}{
    An object to compare.
  }
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  FALSE if the objects can not be compared (either does not have an
  associated model or observed data).

  For hzar.sameModel:
  
  TRUE if the model associated with fitA is identical (same equation,
  parameters, free parameters, and fixed parameter values).

  FALSE otherwise.

  For hzar.sameObsData:

  TRUE if the [["frame"]] of the hzar.obsData object associated with
  fitA is identical to the [["frame"]] of the hzar.obsData object
  associated with fitB.
  
  FALSE otherwise.

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

% \seealso{
% %% ~~objects to See Also as \code{\link{help}}, ~~~
% }
% \examples{
% ##TODO
% }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
