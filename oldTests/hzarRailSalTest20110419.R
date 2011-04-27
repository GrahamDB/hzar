                                            
railData<-read.csv("RailGenData.csv",header=TRUE)
                                        #prepare the data (pull out observations and associate with appropriate likelihood function)

obs.p.F <- doMolecularData1DPops(railData$dist, railData$p.F, railData$p.sal.N)

obs.p.F

#choose Model 0

#starting with 2 degrees of freedom - estimate cline center and width, scaling none

model0.cline.p.F <- makeSimpleCline1D(obs.p.F, scaling="none")

#fitting model to cline
model0.fit.cline.p.F <- fitClineModel(model0.cline.p.F, obs.p.F, mcmc=1e6)

plotClineFunc(model0.fit.cline.p.F);

plot(model0.fit.cline.p.F$mcmc)

#restart fitting model using the previous mcmc chain using the maximum likelihood for all the fitted parameters and estimates the covariance matrix using the 95% credibility of the chain

refit.model0.cline.p.F <- reFitClineFunc(model0.fit.cline.p.F)
plotClineFunc(refit.model0.cline.p.F);
#plot MCMC traces

plot(refit.model0.cline.p.F$mcmc)



attach(refit.model0.cline.p.F$allClines);
credibleLLspace<-data.frame(LL=sort(model.LL),
                            percentile=cumsum(exp(sort(model.LL)))/
                              sum(exp(sort(model.LL))));
detach(refit.model0.cline.p.F$allClines)
plot(percentile~LL,data=credibleLLspace,log="y",xlim=c(max(LL)-6,max(LL)),ylim=c(1e-4,1))


lines(c(125e-6,1)~c(max(LL)-6,max(LL)),data=credibleLLspace,col="blue")
lines(c(25e-4,1)~c(max(LL)-6,max(LL)-2),data=credibleLLspace,col="blue")
lines(c(5e-2,1)~c(max(LL)-6,max(LL)-4),data=credibleLLspace,col="blue")
lines(c(5e-2,1)~c(max(LL)-7,max(LL)-5),data=credibleLLspace,col="blue")
lines(c(25e-4,1)~c(max(LL)-7,max(LL)-3),data=credibleLLspace,col="blue")
lines(c(125e-6,1)~c(max(LL)-7,max(LL)-1),data=credibleLLspace,col="blue")
