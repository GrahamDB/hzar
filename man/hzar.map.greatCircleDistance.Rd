\name{hzar.map.greatCircleDistance}
\alias{hzar.map.greatCircleDistance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  The distance between two points on the Earth's surface.
}
\description{
  The distance along a great circle between two points on a spheroid
  approximation of the Earth's surface.
}
\usage{
hzar.map.greatCircleDistance(lat1, long1, lat2, long2, units = "Km", degrees = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lat1}{
    The latitude of point 1.
  }
  \item{long1}{
    The latitude of point 1.
  }
  \item{lat2}{
    The latitude of point 2.
  }
  \item{long2}{
    The latitude of point 2.
  }
  \item{units}{
    The units of distances returned.  Only the case sensitive values
    \code{"Km"} for kilometers, \code{"miles"} for miles, and
    \code{"nautical"} for nautical miles are valid.
  }
  \item{degrees}{
    Are the latitude and longitude in degrees?
  }
}
\details{
  The Lambert formula is the approximation used to calculate the
  distance, due to its high accuracy and robustness.

  Geometry of Earth
  \code{R = 6371.009 #Earth radius in kilometers}
  \code{earthSphd.r = 298.257223563 #WGS84}
  \code{earthSphd.ep= (2*earthSphd.r -1)/(earthSphd.r-1)^2}
  
  \code{dLat=p2$lat.rad-p1$lat.rad;}
  \code{dLong=p2$long.rad-p1$long.rad;}
  \code{dLong=ifelse( dLong>pi, dLong-2*pi, ifelse(-dLong>pi, dLong+2*pi, dLong));}
  
  
  \code{mLat=(p2$lat.rad+p1$lat.rad)/2;}
  
  \code{reLat1=atan((earthSphd.r -1)*tan(p1$lat.rad) /earthSphd.r )}
  \code{reLat2=atan((earthSphd.r -1)*tan(p2$lat.rad) /earthSphd.r )}
  \code{cenNum=sqrt((cos(reLat2)*sin(dLong))^2+(cos(reLat1)*sin(reLat2)-cos(reLat2)*sin(reLat1)*cos(dLong))^2);}
  \code{cenDen=sin(reLat1)*sin(reLat2)+cos(reLat2)*cos(reLat1)*cos(dLong);}
  \code{central <- atan2(cenNum,cenDen);}
  \code{lFP <- (reLat1+ reLat2)/2 ;}
  \code{lFQ <- (-reLat1+ reLat2)/2 ;}
  \code{lFX <- (central-sin(central))*sin(lFP)^2*cos(lFQ)^2/cos(central/2)^2;}
  \code{lFY <- (central+sin(central))*cos(lFP)^2*sin(lFQ)^2/sin(central/2)^2;}
  
  \code{lambertFormulaeD <- ifelse(central==0,0,R*(central-(lFX+lFY)/(2*earthSphd.r)));}
  
}
\value{
  The great circle distance between points 1 and 2.
}
\references{
  Lambert, W. D. 1942. The distance between two widely separated points on the surface of the earth. J. Washington Academy of Sciences 32 (5): 125-130.
}
\author{
  Graham Derryberry \email{asterion@alum.mit.edu}
}

\examples{
hzar.map.greatCircleDistance(89.5,60,89.5,390)
# 28.87587
hzar.map.greatCircleDistance(-89.5,-90,89.5,90)
# 19981.56
hzar.map.greatCircleDistance(-89.5,90,89.5,90)
# 19869.99
hzar.map.greatCircleDistance(0,90,89.5,90)
# 9934.996
hzar.map.greatCircleDistance(0,90,0,-90)
# 20015.12
hzar.map.greatCircleDistance(0,90,0,180)
# 10007.56
}