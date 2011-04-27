#Grab data to test against.

mannequin.data<- read.csv("molecular-data-robb2.csv");
source("hzarClasses.R");
attach(mannequin.data );
obs.data<-list(L5.b=doMolecularData1DPops(
                 distance=distance,
                 pObs=L5.b,
                 nSamples=L5.n));


cline.bothTails.L5.b<-makeTailedCline1D(obs.data$L5.b,scaling="free")

fit.cline.bothTails.L5.b<-fitClineModel(cline.bothTails.L5.b,obs.data$L5.b,mcmc=1e4);
refit.cline.bothTails.L5.b<-reFitClineFunc(fit.cline.bothTails.L5.b);
plotClineFunc(refit.cline.bothTails.L5.b)

attach(refit.cline.bothTails.L5.b$allClines);
credibleLLspace<-data.frame(LL=sort(model.LL),
                            percentile=cumsum(exp(sort(model.LL)))/
                              sum(exp(sort(model.LL))));
detach(refit.cline.bothTails.L5.b$allClines)
plot(percentile~LL,data=credibleLLspace,log="y",xlim=c(max(LL)-6,max(LL)))

refit.cline.bothTails.L5.b<-reFitClineFunc(refit.cline.bothTails.L5.b);
plotClineFunc(refit.cline.bothTails.L5.b)

attach(refit.cline.bothTails.L5.b$allClines);
credibleLLspace<-data.frame(LL=sort(model.LL),
                            percentile=cumsum(exp(sort(model.LL)))/
                              sum(exp(sort(model.LL))));
detach(refit.cline.bothTails.L5.b$allClines)
plot(percentile~LL,data=credibleLLspace,log="y",xlim=c(max(LL)-6,max(LL)),ylim=c(1e-4,1))


lines(c(125e-6,1)~c(max(LL)-6,max(LL)),data=credibleLLspace,col="blue")
lines(c(25e-4,1)~c(max(LL)-6,max(LL)-2),data=credibleLLspace,col="blue")
lines(c(5e-2,1)~c(max(LL)-6,max(LL)-4),data=credibleLLspace,col="blue")
lines(c(5e-2,1)~c(max(LL)-7,max(LL)-5),data=credibleLLspace,col="blue")
lines(c(25e-4,1)~c(max(LL)-7,max(LL)-3),data=credibleLLspace,col="blue")
lines(c(125e-6,1)~c(max(LL)-7,max(LL)-1),data=credibleLLspace,col="blue")

refit.cline.bothTails.L5.b<-reFitClineFunc(refit.cline.bothTails.L5.b,mcmc=1e7,verbose=1e6);
attach(refit.cline.bothTails.L5.b$allClines);
credibleLLspace<-data.frame(LL=sort(model.LL),
                            percentile=cumsum(exp(sort(model.LL)))/
                              sum(exp(sort(model.LL))));
detach(refit.cline.bothTails.L5.b$allClines)
plot(percentile~LL,data=credibleLLspace,log="y",xlim=c(max(LL)-6,max(LL)),ylim=c(1e-4,1))


lines(c(125e-6,1)~c(max(LL)-6,max(LL)),data=credibleLLspace,col="blue")
lines(c(25e-4,1)~c(max(LL)-6,max(LL)-2),data=credibleLLspace,col="blue")
lines(c(5e-2,1)~c(max(LL)-6,max(LL)-4),data=credibleLLspace,col="blue")
lines(c(5e-2,1)~c(max(LL)-7,max(LL)-5),data=credibleLLspace,col="blue")
lines(c(25e-4,1)~c(max(LL)-7,max(LL)-3),data=credibleLLspace,col="blue")
lines(c(125e-6,1)~c(max(LL)-7,max(LL)-1),data=credibleLLspace,col="blue")

refit.cline.bothTails.L5.b<-reFitClineFunc(refit.cline.bothTails.L5.b,mcmc=1e6,verbose=1e4);

plot(refit.cline.bothTails.L5.b$mcmc)



LLregions=2*(-4:0)
LLcols=c("black","red","blue","green");
clineRegions<-list();
for(iter in 1:4){
clineRegions[[iter]]<-subset(refit.cline.bothTails.L5.b$allClines,((refit.cline.bothTails.L5.b$allClines$model.LL-refit.cline.bothTails.L5.b$maxLL)>LLregions[[iter]])&((refit.cline.bothTails.L5.b$allClines$model.LL-refit.cline.bothTails.L5.b$maxLL)<=LLregions[[iter+1]]));
}
plot(center~width,data=clineRegions[[1]],col=LLcols[[1]],pch="+")
for(iter in 2:4){
points(center~width,data=clineRegions[[iter]],col=LLcols[[iter]],pch="+")
}



attach(lastClineFrame$allClines);
credibleLLspace<-data.frame(LL=sort(model.LL),
                            percentile=cumsum(exp(sort(model.LL)))/
                              sum(exp(sort(model.LL))));
detach(lastClineFrame$allClines)
plot(percentile~LL,data=credibleLLspace,log="y",xlim=c(max(LL)-6,max(LL)),ylim=c(1e-4,1))


lines(c(125e-6,1)~c(max(LL)-6,max(LL)),data=credibleLLspace,col="blue")
lines(c(25e-4,1)~c(max(LL)-6,max(LL)-2),data=credibleLLspace,col="blue")
lines(c(5e-2,1)~c(max(LL)-6,max(LL)-4),data=credibleLLspace,col="blue")
lines(c(5e-2,1)~c(max(LL)-7,max(LL)-5),data=credibleLLspace,col="blue")
lines(c(25e-4,1)~c(max(LL)-7,max(LL)-3),data=credibleLLspace,col="blue")
lines(c(125e-6,1)~c(max(LL)-7,max(LL)-1),data=credibleLLspace,col="blue")
