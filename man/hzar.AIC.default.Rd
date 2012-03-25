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
    A hzar.cline object.
  }
  \item{dataGroup}{
    A hzar.dataGroup object.
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
  Graham Derryberry <email: asterion@alum.mit.edu>

}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{AIC}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

% ## The function is currently defined as
% function (maxLL, param.count) 
% {
%     return(2 * (param.count - maxLL))
%   }

% ## The function is currently defined as
% function (maxLL, param.count, nObs) 
% {
%     return(2 * (param.count - maxLL + (param.count * (param.count + 
%         1)/(nObs - param.count - 1))))
%   }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
