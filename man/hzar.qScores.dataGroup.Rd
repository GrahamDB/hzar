\name{hzar.qScores.dataGroup}
\alias{hzar.qScores.dataGroup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate quantiles for the cline model parameters.
%%  ~~function to do ... ~~
}
\description{
  To be written.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.qScores.dataGroup(dataGroup, probs = c(0.025, 0.5, 0.975))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataGroup}{
%%     ~~Describe \code{dataGroup} here~~
}
  \item{probs}{
%%     ~~Describe \code{probs} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Graham Derryberry \email{asterion@alum.mit.edu}

%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

% ## The function is currently defined as
% function (dataGroup, probs = c(0.025, 0.5, 0.975)) 
% {
%     return(hzar.qScores(as.data.frame(dataGroup$data.mcmc), dataGroup$data.LL$model.LL, 
%         probs = probs))
%   }
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
