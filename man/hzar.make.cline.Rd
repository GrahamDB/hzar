\name{hzar.make.cline}
\alias{hzar.make.cline}
\alias{hzar.cline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Make a \code{hzar.cline} object.
}
\description{
  This method creates a \code{hzar.cline} object, which describes a
  cline model with specific parameter values.  A log likelihood of the
  parameters can be assigned to this object. See
  \code{\link{hzar.gen.cline}} for a simpler method which just requires
  the free parameters and a cline model reference.
}
\usage{
hzar.make.cline(free.parameters, parameters, func, LL, isValid = is.function(func))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{free.parameters}{
    The optimized parameter values for this cline.
  }
  \item{parameters}{
    All of the parameter values for this cline.
  }
  \item{func}{
    The estimator function for this cline.
  }
  \item{LL}{
    The log likelihood of this cline.
  }
  \item{isValid}{
    Is this cline valid?
  }
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  A \code{hzar.cline} object.
  A list with the values:
  \item{param.free}{See \code{free.parameters} above}
  \item{param.all}{See \code{parameters} above}
  \item{clineFunc}{See \code{func} above}
  \item{logLike}{See \code{LL} above}
  \item{isValid}{See \code{isValid} above}
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

% %% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{hzar.gen.cline}}
}
% \examples{
% ##TODO
% }
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
