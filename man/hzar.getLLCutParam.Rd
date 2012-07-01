\name{hzar.getLLCutParam}
\alias{hzar.getLLCutParam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Get the region of parameter space close to the maximum likelihood
}
\description{
  This function returns the range of parameter values that are within
  two log likelihood units of the maximum likelihood for a provided
  character vector of parameters. 
}
\usage{
hzar.getLLCutParam(dataGroups, params, cutValue = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataGroups}{
    Either a \code{\link{hzar.dataGroup}} object, or a list of
    \code{\link{hzar.dataGroup}} objects.
}
\item{params}{
  The parameters to report.
}
\item{cutValue}{
  The number of log likelihood units to retain.
}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  A \code{\link{data.frame}}, with 2 columns for each parameter
  requested, containing the maximum and minimum parameter value within
  \code{cutValue} log likelihood units of the maximum likelihood
  observed. 
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
\code{\link{hzar.dataGroup}}
}
\examples{
##TODO
}
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
