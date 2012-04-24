\name{hzar.qScores.getScores}
\alias{hzar.qScores.getScores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate credibility scores.
}
\description{
  To be written.

  Currently for internal use.
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
    Data samples for numerical integration.
  }
  \item{raw.wt}{
    Log likelihood of weight to apply.
  }
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  The values of the weighted integration performed.
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
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
