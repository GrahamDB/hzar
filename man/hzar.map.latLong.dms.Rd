\name{hzar.map.latLong.dms}
\alias{hzar.map.latLong.dms}
\title{
  Convert D-M-S C strings to degrees
}
\description{
  Translate a character vector of typical map coordinate(s) to a matrix
  of numeric values.
}
\usage{
hzar.map.latLong.dms(coordinates)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coordinates}{
    A character vector, with each string containing one or more
    latitude / longitude measurements, seperated by the cardinal
    direction of the measurement (N/S for latitude, E/W for longitude).  
  }
}

\value{
  A NxM matrix of numeric values, where N is the number of strings in
  \code{coordinates} and M is the greatest number of measurements in a
  single string.
  Each row contains the measurements from each string in
  \code{coordinates}, filling one column from left to right for each
  measurement.  The remaining columns (if any) in a row have the value
  \code{NA}. 
  %% ...
}
\author{
  Graham Derryberry \email{asterion@alum.mit.edu}
}
% \note{
% %%  ~~further notes~~
% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{hzar.map.dms2deg}}
}
\examples{

# A single value:
hzar.map.latLong.dms(c("9-52 N"))
# A 1:2 matrix:
hzar.map.latLong.dms(c("9-52 N 82 W"))
# A 2:2 matrix with NA for the upper right corner:
hzar.map.latLong.dms(c("9-52-34 E" ,"9-52 N 82 W"))
}