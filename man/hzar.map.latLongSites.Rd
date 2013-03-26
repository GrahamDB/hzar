\name{hzar.map.latLongSites}
\alias{hzar.map.latLongSites}
\alias{hzar.map.latLongSites.dms}
\alias{hzar.map.distanceFromSite}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Generate a table of site locations.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
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
    A character vector the same length as \code(siteIDs}. Each value
  must contain the latitude and longitude in DMS format, such as:
  \code{9-22 N 82-34-50 W}
  }
  \item{latLongSites}{
    
  }
  \item{site0}{
    %%     ~~Describe \code{site0} here~~
  }
  \item{units}{
    %%     ~~Describe \code{units} here~~
  }
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
"hzar.map.latLongSites"
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
