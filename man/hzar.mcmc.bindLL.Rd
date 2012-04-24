\name{hzar.mcmc.bindLL}
\alias{hzar.mcmc.bindLL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate a data frame with sampled parameters and log likelihoods.
}
\description{
  To be written. Primarily for internal use.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.mcmc.bindLL(fitRequest,
                 dataGroup = hzar.fit2DataGroup(fitRequest),
                 mcmcData = fitRequest$mcmcRaw,
                 llFunc = dataGroup$llFunc,
                 llData = dataGroup$data.LL,
                 t0 = start(mcmcData),
                 tF = thin(mcmcData))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitRequest}{
    
}
  \item{dataGroup}{
%%     ~~Describe \code{dataGroup} here~~
}
  \item{mcmcData}{
%%     ~~Describe \code{mcmcData} here~~
}
  \item{llFunc}{
%%     ~~Describe \code{llFunc} here~~
}
  \item{llData}{
%%     ~~Describe \code{llData} here~~
}
  \item{t0}{
%%     ~~Describe \code{t0} here~~
}
  \item{tF}{
%%     ~~Describe \code{tF} here~~
}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  A data.frame, with columns for each free parameter and the log
  likelihood of each row.
}

\author{
Graham Derryberry \email{asterion@alum.mit.edu}

}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{hzar.dataGroup}}
}
\examples{
##TODO
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
