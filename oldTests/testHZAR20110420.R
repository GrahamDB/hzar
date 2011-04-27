source("hzarMultiModel.R")

 mannequin.data<-read.csv("molecular-data-robb2.csv",header=TRUE)

test.data<-doMolecularData1DPops(distance=mannequin.data$distance,pObs=mannequin.data$L5.b,nSamples=mannequin.data$L5.n)

library(doMC)
 registerDoMC()


test.mModel<-makeMultiCline1D(test.data)

test.multiCline<-fitMultipleModels(test.mModel,test.data);


plotClineFunc(test.multiCline$clineFrames[[1]])


attach(test.lastCline$allClines);
credibleLLspace<-data.frame(LL=sort(model.LL),
                            percentile=cumsum(exp(sort(model.LL)))/
                              sum(exp(sort(model.LL))));
detach(test.lastCline$allClines)
plot(percentile~LL,data=credibleLLspace,log="y",xlim=c(max(LL)-6,max(LL)),ylim=c(1e-4,1))


lines(c(125e-6,1)~c(max(LL)-6,max(LL)),data=credibleLLspace,col="blue")
lines(c(25e-4,1)~c(max(LL)-6,max(LL)-2),data=credibleLLspace,col="blue")
lines(c(5e-2,1)~c(max(LL)-6,max(LL)-4),data=credibleLLspace,col="blue")
lines(c(5e-2,1)~c(max(LL)-7,max(LL)-5),data=credibleLLspace,col="blue")
lines(c(25e-4,1)~c(max(LL)-7,max(LL)-3),data=credibleLLspace,col="blue")
lines(c(125e-6,1)~c(max(LL)-7,max(LL)-1),data=credibleLLspace,col="blue")
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
credibleLLspace<-data.frame(LL=sort(model.LL),
                            percentile=cumsum(exp(sort(model.LL)))/
                              sum(exp(sort(model.LL))));
detach(clineFrame$allClines)
plot(log(percentile)~LL,data=credibleLLspace,xlim=c(max(LL)-6,max(LL)),ylim=c(-9,0))


lines(log(c(625e-8,1))~c(max(LL)-8,max(LL)),data=credibleLLspace,col="blue")
lines(log(c(125e-6,1))~c(max(LL)-8,max(LL)-2),data=credibleLLspace,col="blue")
lines(log(c(25e-4,1))~c(max(LL)-8,max(LL)-4),data=credibleLLspace,col="blue")
lines(log(c(5e-2,1))~c(max(LL)-7,max(LL)-5),data=credibleLLspace,col="blue")
lines(log(c(25e-4,1))~c(max(LL)-7,max(LL)-3),data=credibleLLspace,col="blue")
lines(log(c(125e-6,1))~c(max(LL)-7,max(LL)-1),data=credibleLLspace,col="blue")
}
