\name{hzar.dataGroup.null}
\alias{hzar.dataGroup.null}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Datagroup placeholder for the null model (frequency independent of location)
%%  ~~function to do ... ~~
}
\description{
  To be written.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.dataGroup.null(obsData)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obsData}{
%%     ~~Describe \code{obsData} here~~
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
% function (obsData) 
% {
%     clineLLfunc <- hzar.make.LLfunc.null(obsData)
%     nPS <- 100
%     pSeries <- (1:(nPS - 1))/(nPS)
%     pExpected <- sum(obsData$frame$n * obsData$frame$obsFreq)/sum(obsData$frame$n)
%     pSeries <- c(pSeries, pExpected)
%     data.mcmc <- mcmc(data = data.frame(pVal = pSeries))
%     return(hzar.make.dataGroup(data.mcmc = data.mcmc, llFunc = clineLLfunc, 
%         ML.cline = hzar.gen.cline(list(pVal = pExpected), clineLLfunc)))
%   }
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
