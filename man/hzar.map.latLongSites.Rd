\name{hzar.map.latLongSites}
\alias{hzar.map.latLongSites}
\alias{hzar.map.latLongSites.dms}
\alias{hzar.map.distanceFromSite}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate a table of site locations.
}
\description{
  Methods for processing tables of site location:
  Generate a table of latitude and longitude of site locations from
  either precalculated values or GPS coordinates.
  Generate a table of distances from an included site, using a table of
  latitude and longitude of site locations. 
  
}
\usage{
hzar.map.latLongSites(siteIDs, site.lat, site.long, degrees = TRUE)
hzar.map.latLongSites.dms(siteIDs, coordinates)
hzar.map.distanceFromSite(latLongSites, site0, units = "Km")

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{siteIDs}{
    A character vector used to identify each site.
  }
  \item{site.lat}{
    A numeric vector of site longitudes.
    
  }
  \item{site.long}{
    A numeric vector of site longitudes.
  }
  \item{degrees}{
    Are \code{site.lat} and \code{site.long} in degrees? If
    \code{FALSE}, \code{site.lat} and \code{site.long} should be in radians.
  }
  
  \item{coordinates}{
    A character vector the same length as \code{siteIDs}. Each value
    must contain the latitude and longitude in DMS format, such as:
    \code{9-22 N 82-34-50 W}
  }
  \item{latLongSites}{
    The result of either the \code{hzar.map.latLongSites} or
    \code{hzar.map.latLongSites.dms} method. 
  }
  \item{site0}{
    The ID string for the site place at 0.
  }
  \item{units}{
    The units of distances returned.  Only the case sensitive values
    \code{"Km"} for kilometers, \code{"miles"} for miles, and
    \code{"nautical"} for nautical miles are valid.
  }
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  A \code{\link{data.frame}}:
  For all methods:
  \item{site }{A character vector used to identify each site.}
  For \code{hzar.map.latLongSites} and \code{hzar.map.latLongSites.dms}:
  \item{lat.rad}{ The site latitude in radians. }
  \item{long.rad}{The site longitude in radians. }
  \item{lat.deg}{ The site latitude in degrees. }
  \item{long.deg}{The site longitude in degrees. }
  For \code{hzar.map.distanceFromSite}:
  \item{distance }{ The distance to each site from a common origin. }
  
}

\author{
  Graham Derryberry \email{asterion@alum.mit.edu}
}
\note{
  Distances in \code{hzar.map.distanceFromSite} are calcutated using the
  method \code{\link{hzar.map.greatCircleDistance}}.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{hzar.map.greatCircleDistance}}
}
\examples{

a=hzar.map.latLongSites(c("Norway"),60.4,11)
b=hzar.map.latLongSites.dms(c("Norway","Sweden"),c("60-24 N 11 E","58 N 15 E"))
hzar.map.distanceFromSite(b,"Norway")
}
