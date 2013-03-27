\name{hzar.map.dms2deg}
\alias{hzar.map.dms2deg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Convert degrees, minutes, seconds, direction to degrees.
}
\description{
  This method converts angular measurements from degree, minute, second,
  direction notation to decimal notation.
}
\usage{
hzar.map.dms2deg(deg, min, sec, dir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{deg}{
    The degrees of the angular measurement.
  }
  \item{min}{
    The minutes of the angular measurement.
  }
  \item{sec}{
    The seconds of the angular measurement.
  }
  \item{dir}{
    A character vector indicating the direction of the angular
    measurement. 
    The case insensitive values \code{"N"}, \code{"S"}, \code{"E"},
    \code{"W"}, \code{"West"}, \code{"North"}, \code{"South"} and
    \code{"East"} are all valid, any others will result in an error. 
  }
}
% \details{
%   %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  A numeric vector of angular measurements.
}
\author{
  Graham Derryberry \email{asterion@alum.mit.edu}
}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{hzar.map.latLong.dms}}
  \code{\link{hzar.map.latLongSites.dms}}
}
\examples{

hzar.map.dms2deg(9,24,0,"N")

hzar.map.dms2deg(82,33,50,"W")

}