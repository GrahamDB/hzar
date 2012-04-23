\name{hzar.make.dataGroup}
\alias{hzar.make.dataGroup}
\alias{hzar.dataGroup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Make a hzar.dataGroup object with given parameters.
}
\description{
  To be written.
}
\usage{
hzar.make.dataGroup(data.mcmc,
                    llFunc,
                    ML.cline = NULL,
                    doPar = FALSE,
                    data.LL = hzar.eval.clineLL(llFunc = llFunc,
                                                data = data.mcmc,
                                                doPar = doPar),
                    data.param = as.data.frame(data.mcmc),
                    obsData = hzar.extract.obsData(llFunc))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.mcmc}{
    %%     ~~Describe \code{data.mcmc} here~~
  }
  \item{llFunc}{
    %%     ~~Describe \code{llFunc} here~~
  }
  \item{ML.cline}{
    %%     ~~Describe \code{ML.cline} here~~
  }
  \item{doPar}{
    %%     ~~Describe \code{doPar} here~~
  }
  \item{data.LL}{
    %%     ~~Describe \code{data.LL} here~~
  }
  \item{data.param}{
    %%     ~~Describe \code{data.param} here~~
  }
  \item{obsData}{
    %%     ~~Describe \code{obsData} here~~
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  Object of class \code{hzar.dataGroup}
\item{llFunc}{Description of 'comp2'}
\item{data.mcmc}{Description of 'comp1'}
\item{data.param}{Description of 'comp2'}
\item{data.LL}{Description of 'comp2'}
\item{ML.cline}{Description of 'comp2'}
\item{obsData}{Description of 'comp2'}
%% ...
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

\seealso{
  \code{\link{hzar.eval.clineLL}}
  \code{\link{hzar.extract.obsData}}
  \code{\link{hzar.dataGroup.add}}
}
\examples{
##TODO
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
