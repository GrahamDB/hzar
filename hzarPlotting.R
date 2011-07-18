
plotClineFunc <- function(clineFrame){

  attach(clineFrame$data$frame);
  x<-seq(min(dist),max(dist),length.out=101);
  detach();
  plot(obsFreq ~ dist, data=clineFrame$data$frame);
  lines(clineFrame$maxLL.cline(x)~x);

}




plotLLspace<-function(clineFrame){
  attach(clineFrame$allClines);
credibleLLspace<-data.frame(LL=sort(model.LL),
                            percentile=cumsum(exp(sort(model.LL)))/
                              sum(exp(sort(model.LL))));
detach(clineFrame$allClines)
plot(percentile~LL,data=credibleLLspace,log="y",xlim=c(max(LL)-6,max(LL)),ylim=c(1e-4,1))


lines(c(625e-8,1)~c(max(LL)-8,max(LL)),data=credibleLLspace,col="blue")
lines(c(125e-6,1)~c(max(LL)-8,max(LL)-2),data=credibleLLspace,col="blue")
lines(c(25e-4,1)~c(max(LL)-8,max(LL)-4),data=credibleLLspace,col="blue")
lines(c(5e-2,1)~c(max(LL)-7,max(LL)-5),data=credibleLLspace,col="blue")
lines(c(25e-4,1)~c(max(LL)-7,max(LL)-3),data=credibleLLspace,col="blue")
lines(c(125e-6,1)~c(max(LL)-7,max(LL)-1),data=credibleLLspace,col="blue")
}

plotLLspaceDelta<-function(clineFrame){
  attach(clineFrame$allClines);
  credibleLLspace<-data.frame(LL=sort(model.LL)-max(model.LL),
                              percentile=cumsum(exp(sort(model.LL)))/
                              sum(exp(sort(model.LL))));
  detach(clineFrame$allClines)
  attach(clineFrame$allClines.old);
  credibleLLspace.old<-data.frame(LL=sort(model.LL)-max(clineFrame$allClines$model.LL),
                              percentile=cumsum(exp(sort(model.LL)))/
                              sum(exp(sort(model.LL))));
  detach(clineFrame$allClines.old)
  attach(clineFrame$allClines.new);
  credibleLLspace.new<-data.frame(LL=sort(model.LL)-max(clineFrame$allClines$model.LL),
                              percentile=cumsum(exp(sort(model.LL)))/
                              sum(exp(sort(model.LL))));
  detach(clineFrame$allClines.new)
  plot(-log(percentile)+3*LL/2~LL,data=credibleLLspace,xlim=c(-6,0),ylim=c(-3,3))
  
  points(-log(percentile)+3*LL/2~LL,data=credibleLLspace.old,pch="X",col="orange")
  points(-log(percentile)+3*LL/2~LL,data=credibleLLspace.new,pch="+",col="green")
  llMarkers=data.frame(cSpace=log(c(625e-8,125e-6,25e-4,5e-2,1)) ,
    lSpace=seq(-8,0,by=2) );
  lines(-cSpace[4]+3*lSpace/2~lSpace,data=llMarkers,col="red");
  lines(c(0,0)~c(-7,1),col="blue")
  lines(-cSpace[3]+3*lSpace/2~lSpace,data=llMarkers,col="red");
  
  lines(-log(0.01)+3*lSpace/2~lSpace,data=llMarkers,col="red");
## lines(log(c(125e-6,1))~c(max(LL)-8,max(LL)-2),data=credibleLLspace,col="blue")
## lines(log(c(25e-4,1))~c(max(LL)-8,max(LL)-4),data=credibleLLspace,col="blue")
## lines(log(c(5e-2,1))~c(max(LL)-7,max(LL)-5),data=credibleLLspace,col="blue")
## lines(log(c(25e-4,1))~c(max(LL)-7,max(LL)-3),data=credibleLLspace,col="blue")
## lines(log(c(125e-6,1))~c(max(LL)-7,max(LL)-1),data=credibleLLspace,col="blue")
}
