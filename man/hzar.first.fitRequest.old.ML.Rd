\name{hzar.first.fitRequest.old.ML}
\alias{hzar.first.fitRequest.old.ML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate a ML based hzar.fitRequest using a meta model structure.
}
\description{
  This method generates a \code{\link{hzar.fitRequest}} object suitable
  for \code{\link{hzar.doFit}}.
}
\usage{
hzar.first.fitRequest.old.ML(model, obsData, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
    A \code{\link{clineMetaModel}} object. 
    %%     ~~Describe \code{model} here~~
  }
  \item{obsData}{
    The \code{\link{hzar.obsData}} object to which the meta model is to be fit.
    %%     ~~Describe \code{obsData} here~~
  }
  \item{verbose}{
    Should \code{\link{MCMCmetrop1R}} be verbose?
    %%     ~~Describe \code{verbose} here~~
  }
}

% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  A \code{\link{hzar.fitRequest}} object.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
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
  \code{\link{hzar.fitRequest}}
  \code{\link{hzar.obsData}}
  \code{\link{clineMetaModel}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

% ## The function is currently defined as
% function (model, obsData, verbose = TRUE) 
% {
%     if (verbose) {
%         mcmcParam <- hzar.default.mcmc
%     }
%     else {
%         mcmcParam <- hzar.quiet.mcmc
%     }
%     modelParam <- splitParameters(model$parameterTypes)
%     clineLLfunc <- hzar.make.clineLLfunc.old.ML(names(modelParam$init), 
%         modelParam$fixed, model$req, model$func, obsData$model.LL)
%     covMatrix <- NULL
%     try(covMatrix <- hzar.cov.rect(clineLLfunc, modelParam$lower, 
%         modelParam$upper, random = 10000))
%     return(hzar.make.fitRequest(modelParam, covMatrix, clineLLfunc, 
%         mcmcParam))
%   }
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
