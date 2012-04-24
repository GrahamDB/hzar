\name{hzar.make.cline}
\alias{hzar.make.cline}
\alias{hzar.cline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Make a \code{hzar.cline} object.
%%  ~~function to do ... ~~
}
\description{
  To be written.
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
 \item{param.free}{\code{free.parameters}}
 \item{param.all}{\code{parameters}}
 \item{clineFunc}{\code{func}}
 \item{logLike}{\code{LL}}
 \item{isValid}{\code{isValid}}
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

% \seealso{
% %% ~~objects to See Also as \code{\link{help}}, ~~~
% }
\examples{
##TODO
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
