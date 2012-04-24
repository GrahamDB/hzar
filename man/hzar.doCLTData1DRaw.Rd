\name{hzar.doCLTData1DRaw}
\alias{hzar.doCLTData1DRaw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Create an hzar.obsData object using a table of individual traits.
%%  ~~function to do ... ~~
}
\description{
  To be written.

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hzar.doCLTData1DRaw(distance, traitValue)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{distance}{
%%     ~~Describe \code{distance} here~~
}
  \item{traitValue}{
%%     ~~Describe \code{traitValue} here~~
}
}
\details{
  If for any population site, there is only a small number of samples
  taken, warnings will be issued.

  If at any population site, the sample variance is 0, a warning is
  issued, and additional variance is included by estimating the amount
  of variance ignored due to measurement error.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 A \code{\link{hzar.obsData}} object.
}
% \references{
% %% ~put references to the literature/web site here ~
% }
\author{
Graham Derryberry \email{asterion@alum.mit.edu}
}

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
