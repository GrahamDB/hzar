\name{hzar.gen.cline}
\alias{hzar.gen.cline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Make a \code{hzar.cline} object using the given parameters and model.
%%  ~~function to do ... ~~
}
\description{
  Make a \code{\link{hzar.cline}} object using the given parameters and model.

}
\usage{
hzar.gen.cline(free.parameters, fitRequest)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{free.parameters}{
    A named list of free parameter values.
  }
  \item{fitRequest}{
    An object refering to the model that should be used, such as a
    \code{\link{hzar.fitRequest}} or \code{\link{hzar.dataGroup}}
    object. 
  }
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  The \code{\link{hzar.cline}} object requested.
}

\author{
Graham Derryberry \email{asterion@alum.mit.edu}
}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{hzar.cline}}
}
\examples{
##TODO
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
