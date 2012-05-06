\name{hzar.getCredParamRed}
\alias{hzar.getCredParamRed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Create a 95\% credibility hzar.fzCline object
}
\description{
  Generate a \code{\link{hzar.fzCline}} of the parameter subset selected
  from a \code{\link{hzar.dataGroup}} by \code{\link{hzar.getCredParam}}
  with a \code{rejectionPercent} of \code{0.05}.
}
\usage{
hzar.getCredParamRed(dataGroup)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataGroup}{
    The \code{\link{hzar.dataGroup}} to analyze.
  }
}
% \details{
%   %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  A \code{\link{hzar.fzCline}} object for the samples selected by
  \code{\link{hzar.getCredParam}}. 
}
% \references{
% %% ~put references to the literature/web site here ~
% }
\author{
  Graham Derryberry \email{asterion@alum.mit.edu}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{hzar.make.fzCline}}
  \code{\link{hzar.getCredParam}}
  \code{\link{hzar.dataGroup}}
  \code{\link{hzar.fzCline}}
}
% \examples{
% ##TODO
% }
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
