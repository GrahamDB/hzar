\name{hzar.make.LLfunc.null}
\alias{hzar.make.LLfunc.null}
\title{
  Likelihood function for the null model
}
\description{
  To be written.
  Likelihood function for the null model (frequency independent of
  location). 
}
\usage{
hzar.make.LLfunc.null(obsData,
                      model.LL = obsData$model.LL,
                      LLrejectedModel = -1e+08) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obsData}{
    A \code{\link{hzar.obsData}} object.
  }
  \item{model.LL}{
    The likelihood function for the observed data.
  }
  \item{LLrejectedModel}{
    A finite likelihood to return when the likelihood calculated is
    \code{-INF}. 
  }
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  A function of \code{pEst}, that returns the likelihood that
  \code{pEst} is the mean frequency at all observed populations.  
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
  \code{\link{hzar.obsData}}
}
\examples{
##TODO
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
