## Intended for plotting functions using the new hzar class
## structures.

hzar.yLabel<-function(x) {
  oD<-hzar.extract.obsData(x);
  if(inherits(oD, "clineSampleData1D")){
    ## Add name extractor
    return("Frequency");
  }
  if(inherits(oD, "clineSampleData1DCLT")){
    ## Add name extractor
    return("Trait Value");
  }
  if( is.list(oD)){
    if(prod(as.logical(lapply(oD,inherits,what="clineSampleData1D")))==1)
      return("Frequency");
    if(prod(as.logical(lapply(oD,inherits,what="clineSampleData1DCLT")))==1)
      return("Trait Value");
    
  }
  return("");
}

hzar.plot.obsData<-function(x,type="p",pch="+",xlab="Distance",ylab=hzar.yLabel(x),add=FALSE,...){
  oD<-hzar.extract.obsData(x);
  if(inherits(oD,c("hzar.obsData",
                   "clineSampleData1D",
                   "clineSampleData1DCLT"))){
    oDF<-oD$frame;
    yData<-numeric(length(oDF$dist));
    if(inherits(oD,"clineSampleData1D"))
      yData<-oDF$obsFreq
    if(inherits(oD,"clineSampleData1DCLT"))
      yData<-oDF$obsMean
    if(add){
      points(y=yData,
           x=oDF$dist,
           type=type,pch=pch,...);
    }else{
      plot(y=yData,
           x=oDF$dist,
           type=type,pch=pch,xlab=xlab,ylab=ylab,...);
    }
  }
}


hzar.plot.cline<-function(cline,add=FALSE,...){
  if(inherits(cline,"hzar.cline"))
    curve(cline$clineFunc(x),add=add,...);
  if(inherits(cline,c("hzar.dataGroup","hzar.fitRequest"))){
    dataGroup<-hzar.fit2DataGroup(cline);
    hzar.plot.obsData(dataGroup,add=add,...);
    hzar.plot.cline(hzar.get.ML.cline(dataGroup),add=TRUE,...);
  }
  if(inherits(cline,c("hzar.obsDataGroup"))){
    hzar.plot.obsData(cline,add=add,...);
    lapply(cline$data.groups,function(dataGroup) hzar.plot.cline(hzar.get.ML.cline(dataGroup),add=TRUE,... ));
  }
}


hzar.mcmc.bindLL <-
  function(fitRequest,
           dataGroup=hzar.fit2DataGroup(fitRequest),
           mcmcData=ifelse(
             inherits(fitRequest,"hzar.fitRequest"),
             fitRequest$mcmcRaw,
             as.mcmc(dataGroup$data.mcmc)),
           llData=dataGroup$data.LL,
           t0=start(mcmcData),
           tF=thin(mcmcData)){
  data<-cbind(mcmcData,llData);
## print( data[1:50,]);
  result<-mcmc(data=as.matrix(data),start=t0,thin=tF);
  return(result);
}

              

hzar.plot.fzCline<-function(dataGroup,fzCline=hzar.getCredParamRed(dataGroup) ,type="p",pch="+",col="black",fzCol="gray",...){
  hzar.plot.obsData(dataGroup,col="transparent",...);
  xSeries<-seq(from=par("usr")[1],to=par("usr")[2],length.out=109)
  if(par("xaxs")=="r")
    xSeries<-xSeries[2:108];
  
  fzCoor<-fzCline$fzCline(xSeries);
  polygon(x=c(fzCoor$x,rev(fzCoor$x)),y=c(fzCoor$yMin,rev(fzCoor$yMax)),border=fzCol,col=fzCol);
  lines(x=xSeries,y=dataGroup$ML.cline$clineFunc(xSeries),col=col);
  hzar.plot.obsData(dataGroup,col=col,type=type,pch=pch,add=TRUE);
}

hzar.overPlot.fzCline<-function(dataGroupSet,
                                fzClineSet=sapply(dataGroupSet,
                                  hzar.getCredParamRed, simplify=FALSE),
                                type="p",
                              ##  pch="+", col="black",
                                fzDens=8,
                                fzShadeAngle=((1:length(dataGroupSet))*180
                                              ) %/% (1+length(dataGroupSet)),
                                ...){
  
}
