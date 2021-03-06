\name{hzar.extract.obsData}
\alias{hzar.extract.obsData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Extract the observation data used by the optimizer.
%%  ~~function to do ... ~~
}
\description{
  Most hzar objects have at least an indirectly if not directly
  associated hzar.obsData object. This function returns that
  hzar.obsData object.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.extract.obsData(fitRequest)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitRequest}{
    A \code{\link{hzar.dataGroup}}, \code{\link{hzar.obsDataGroup}},
    \code{\link{hzar.fitRequest}} or \code{\link{hzar.obsData}} object,
    or a likelihood function generated by hzar.
  }
}

% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }

\value{
  The associated hzar.obsData object.
}
% \references{
% %% ~put references to the literature/web site here ~
% }
\author{
  Graham Derryberry \email{asterion@alum.mit.edu}

}
\note{
  This function if passed a hzar.obsData objects returns the same
  object.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

% \seealso{
%   % \code{\link{hzar.}}
% }
\examples{
##TODO
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
