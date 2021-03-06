\name{hzar.eval.clineLL}
\alias{hzar.eval.clineLL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculate the Log Likehoods of the table of parameters provided.
}
\description{
  Using the likelihood function and the table of parameter values
  provided, calculate the likelihood of each row of parameter values.  
}
\usage{
hzar.eval.clineLL(data, llFunc, doPar = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
    A \code{\link{data.frame}} of the free parameter values.  Each
    column name should match the corresponding parameter name. 
  }
  \item{llFunc}{
    The log likelihood function to use.
  }
  \item{doPar}{
    If TRUE, use \code{\link{\%dopar\%}} to iterate over the rows of
    \code{data}.  
  }
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  A data.frame with a single column (\code{model.LL}) containing the log
  likelihoods for each row of \code{data}.
}
% \references{
% %% ~put references to the literature/web site here ~
% }
\author{
Graham Derryberry \email{asterion@alum.mit.edu}
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{\%dopar\%}}
}
% \examples{
% ##TODO
% }
% % Add one or more standard keywords, see file 'KEYWORDS' in the
% % R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
