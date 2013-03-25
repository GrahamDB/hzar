## Calculate the distance between two points on the earth.

hzar.map.greatCircleDistance <-  function( lat1, long1, lat2, long2,units="Km",degrees=TRUE){
  latLong.d(latLong.p(lat1,long1,degrees),
            latLong.p(lat2,long2,degrees),
            units);
}

latLongRE.NE <- "n(orth)?|e(ast)?"
latLongRE.SW <- "w(est)?|s(outh)?"
latLongRE.NESW <- "n(orth)?|e(ast)?|w(est)?|s(outh)?"

hzar.map.dms2deg <- function(deg,min,sec,dir){
  
  if(!all(grepl(latLongRE.NESW,dir,ignore.case=TRUE)))
    stop(paste("Direction",dir[which(!grepl(latLongRE.NESW,dir,ignore.case=TRUE))],"not recognized."))
  
  res <- deg+min/60+sec/3600
  if(length(dir)==1){
    if(grepl(latLongRE.NE,dir,ignore.case=TRUE))
     return(res)
    return(-res)
  }
  ifelse(grepl(latLongRE.NE,dir,ignore.case=TRUE),res,-res)
}

pdms2deg <- function(x){
  if(length(x)==2)
    return(hzar.dms2deg(as.numeric(x[1]),0,0,x[2]))
  if(length(x)==3)
    return(hzar.dms2deg(as.numeric(x[1]),as.numeric(x[2]),0,x[3]))
  if(length(x)==4)
    return(hzar.dms2deg(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[3]),x[4]))
  stop(paste(paste(x,collapse=" "),"not understood."))
}
  

hzar.map.latLong.dms <- function(coordinates){
  res <- lapply(strsplit(test.coor,"(?<=[NSEW]) ",perl=TRUE), strsplit,"-| ")
  result <- matrix(NA,nrow=length(coordinates),max(sapply(res,length)))
  res <- lapply(res,lapply,pdms2deg)
  for(iter in 1:length(coordinates)){
    result[iter,1:length(res[[iter]])] <- as.numeric(res[[iter]])
  }
  return(result)
  
}
hzar.map.latLongSites <- function(siteIDs,site.lat,site.long,degrees=TRUE){
  cbind(site=siteIDs,latLong.p(site.lat,site.long,degrees))
}
hzar.map.latLongSites.dms <- function(siteIDs,coordinates){
  res <- hzar.latLong.dms(coordinates)
##   if(length(siteIDs)!=nrow(res))
##     stop(
  cbind(site=siteIDs,latLong.p(site.lat,site.long,TRUE))
}
hzar.map.distanceFromSite <-  function(latLongSites,site0,units="Km"){
  if(is.character(site0)&&site0[[1]] %in% latLongSites$site){
    i=which(latLongSites$site %in% site0[1]);
    return(latLong.d(latLongSites[i[1],],latLongSites,units))
  }else{
    stop(paste(site0,"not listed in {",paste(latLongSites$site,collapse=", "),"}"))
  }
}

## Useful values
earthD.Km   = 6371.009 #kilometers
earthD.mile = 3958.761 #statute miles
earthD.ntMl = 3440.069 #nautical miles
earthSphd.r = 298.257223563 #WGS84
earthSphd.ep= (2*earthSphd.r -1)/(earthSphd.r-1)^2
## latLongCoors <- function(lat,long){
##   return(cbind(lat=lat,long=long));
## }

latLong.p <- function(lat,long,degrees=TRUE){
  if(degrees){
    latD=lat;longD=long;
    lat=pi*lat/180
    long=pi*long/180
    return(data.frame(lat.rad=lat,long.rad=long,lat.deg=latD,long.deg=longD));
  }
  
  latD=180*lat/pi
  longD=180*long/pi
  return(data.frame(lat.rad=lat,long.rad=long,lat.deg=latD,long.deg=longD));
}


latLong.d <- function(p1,p2,units="Km"){
  if(identical(units,"Km")){ 
    R=earthD.Km;
  }else{
    if(identical(units,"miles")){
       R=earthD.mile;
     } else {
       if(identical(units,"nautical")){
         R=earthD.ntMl;
       }else{
         stop(units +"not identifiable");
       }
     }
  }
  dLat=p2$lat.rad-p1$lat.rad;
  dLong=p2$long.rad-p1$long.rad;
  dLong=ifelse(
    dLong>pi,
    dLong-2*pi,
    ifelse(-dLong>pi,
           dLong+2*pi,
           dLong));
  
  reLat1=atan((earthSphd.r -1)*tan(p1$lat.rad) /earthSphd.r )
  reLat2=atan((earthSphd.r -1)*tan(p2$lat.rad) /earthSphd.r )
  cenNum=sqrt((cos(reLat2)*sin(dLong))^2+(cos(reLat1)*sin(reLat2)-cos(reLat2)*sin(reLat1)*cos(dLong))^2);
  cenDen=sin(reLat1)*sin(reLat2)+cos(reLat2)*cos(reLat1)*cos(dLong);
  central <- atan2(cenNum,cenDen);
  lFP <- (reLat1+ reLat2)/2 ;
  lFQ <- (-reLat1+ reLat2)/2 ;
  lFX <- (central-sin(central))*sin(lFP)^2*cos(lFQ)^2/cos(central/2)^2;
  lFY <- (central+sin(central))*cos(lFP)^2*sin(lFQ)^2/sin(central/2)^2;
  
  lambertFormulaeD <- ifelse(central==0,0,R*(central-(lFX+lFY)/(2*earthSphd.r)));
  return( lambertFormulaeD);
  
}
