\name{hzar.meta.init}
\alias{hzar.meta.init}
\alias{hzar.meta.tune}
\alias{hzar.meta.fix}
\alias{hzar.meta.lower}
\alias{hzar.meta.upper}
\alias{hzar.meta.init<-}
\alias{hzar.meta.tune<-}
\alias{hzar.meta.fix<-}
\alias{hzar.meta.lower<-}
\alias{hzar.meta.upper<-}

\title{
  Observe and Alter the model parameters in the clineMetaModel
}
\description{
  This is a collection of methods to get or set attributes of the
  various model parameters.
}
\usage{
hzar.meta.init(x)
hzar.meta.init(x) <- value
hzar.meta.tune(x)
hzar.meta.tune(x) <- value
hzar.meta.fix(x)
hzar.meta.fix(x) <- value
hzar.meta.lower(x)
hzar.meta.lower(x) <- value
hzar.meta.upper(x)
hzar.meta.upper(x) <- value
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    The \code{\link{clineMetaModel}} to use.
  }
  \item{value}{
    The new value or values to set.
  }
}
%% \details{
%% %%  ~~ If necessary, more details than the description above ~~
%% }
\value{
  Returns a list, with one numeric or boolean value per cline parameter.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%% \references{
%% %% ~put references to the literature/web site here ~
%% }
\author{
Graham Derryberry \email{asterion@alum.mit.edu}
%%  ~~who you are~~
}
%% \note{
%% %%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{clineMetaModel}}
}
%% \examples{

%% }
