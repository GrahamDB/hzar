\name{hzar.qScores.getScores}
\alias{hzar.qScores.getScores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate credibility scores.
}
\description{
  Determine the series of weights on a series of intervals from a series
  of weighted samples. 
}
\usage{
hzar.qScores.getScores(xSeries, dSeries, raw.x, raw.wt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xSeries}{
    Locations to divide integration intervals (must be in ascending order).
  }
  \item{dSeries}{
    Widths of locations.
  }
  \item{raw.x}{
    Data samples of the distribution.
  }
  \item{raw.wt}{
    Log of the weight to apply to each sample.
  }
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  The weights on the intervals given.
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
% ##
% }
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
\keyword{ internal }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
