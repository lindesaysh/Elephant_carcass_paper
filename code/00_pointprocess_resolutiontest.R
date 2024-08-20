require(spatstat)
require(ggplot2)
devtools::load_all('C:/Users/lass/Documents/GitHub/MRSea')

# load in pan data (2.5km internal buffer)
largepan<-read.csv('data/largepan.csv')
largepan<-largepan/1000
largepan.list<-list(x=largepan$x, y=largepan$y)

# ENP boundary plus 20km
fence<-read.csv('data/databoundary.csv')
fence<-fence[nrow(fence):1,]/1000
fence.list<-list(x=fence$x.pos, y=fence$y.pos)

# read in points data
datcomb<- read.csv("data/datcomb.csv", header=TRUE)

#subset just to get the presences
dat.pres<-datcomb[datcomb$response==1,]
dat<-datcomb[datcomb$response==1,c('x.pos', 'y.pos')]/1000

ggplot(dat.pres) + geom_point(aes(x.pos, y.pos), alpha=1/2) + facet_wrap(~year)
ggplot(dat.pres) + geom_point(aes(x.pos, y.pos), alpha=1/2) + facet_wrap(~month)

# create boundary listing so that quad points created in correct area
Z<-owin(poly=list(fence.list, largepan.list))

x<-dat$x.pos
y<-dat$y.pos

P<-ppp(x, y, poly=list(fence.list, largepan.list))
plot(P)



spacings<-c(5, 4, 3, 2, 1.5, 1.25, 1)#, 0.5, 0.25)
startknots<-c(10, 20, 30, 40)
basis<-c("exponential", "gaussian")
outputs<-NULL
betaall<-NULL
betaseall<-NULL
for(b in basis){
  for(sk in startknots){
    beta<-matrix(NA, nrow=(sk+1), ncol=length(spacings))
    betase<-matrix(NA, nrow=(sk+1), ncol=length(spacings))
    for(i in 1:length(spacings)){
      #i=1
      dummyQpts<-quadscheme(P, eps=spacings[i])$dummy
      dpts<-data.frame(x.pos=dummyQpts$x, y.pos=dummyQpts$y)
      plot(dummyQpts)
      plot(P, add = TRUE, cex = 0.5)
      
      analysisdat<-rbind(data.frame(dat, response=1), data.frame(dpts, response=0))
      kg<- cbind(dat$x.pos, dat$y.pos)
      duplicatedid<-which(duplicated(kg)==TRUE)
      if(length(duplicatedid)>0){ kg<-kg[-duplicatedid,]}
      distMats<-makeDists(cbind(analysisdat$x.pos, analysisdat$y.pos), kg)
      
      #
      p.wt<-rep(1e-6, nrow(analysisdat))
      p.wt[analysisdat$response==0]<-summary(dummyQpts)$window$area/sum(analysisdat$response==0)
      
      initialModel<- glm(response/p.wt ~ 1 , family='poisson', data=analysisdat, weights=p.wt)
      
      salsa2dlist<-list(fitnessMeasure = 'BIC',
                        knotgrid = kg,
                        startKnots=sk,
                        minKnots=sk,
                        maxKnots=sk,
                        cv.opts=list(K=5, cv.gamMRSea.seed=1))
      
      r_seq<-getRadiiChoices(10, distMats$dataDist, basis=b)
      
      salsa2dOutput<-runSALSA2D(initialModel, salsa2dlist, d2k=distMats$dataDist, k2k=distMats$knotDist, splineParams=NULL,suppress.printout = TRUE, basis=b, chooserad = FALSE)
      
      beta[,i]<-as.numeric(summary(salsa2dOutput$bestModel)$coeff[,1])
      betase[,i] <- as.numeric(summary(salsa2dOutput$bestModel)$coeff[,2])
      
      quilt.plot(analysisdat$x.pos, analysisdat$y.pos, fitted(salsa2dOutput$bestModel), asp=1, nrow=64, ncol=64)
      points(dat$x.pos, dat$y.pos, pch=20)
      
      loglik<-sum(p.wt*(((analysisdat$response/p.wt)*log(fitted(salsa2dOutput$bestModel))) - fitted(salsa2dOutput$bestModel)))
      
      outputs<-rbind(outputs, 
                     data.frame(Basis = b, 
                                startknots = sk,
                                gridspacing = spacings[i],
                                numabs = nrow(dpts),
                                AIC = salsa2dOutput$bestModel$aic,
                                Loglik = logLik(salsa2dOutput$bestModel), 
                                loglik = loglik, 
                                numparams = length(coefficients(salsa2dOutput$bestModel))))
      print(outputs)
    }
    betaall<-rbind(betaall, cbind(sk, spacings[i], beta))
    betaseall<-rbind(betaall, cbind(sk, spacings[i], betase))
  }
}


outputs$startknots<-as.factor(outputs$startknots)

saveRDS(outputs, file='results/pp_resolutions_160621')

output<-readRDS((file='results/pp_resolutions_160621'))

#require(ggplot2)
require(dplyr)

p1 <- mutate(output, 
       group=paste0(Basis, startknots)) %>%
ggplot() + 
  geom_line(aes(numabs, loglik, group = group, colour=startknots, linetype=Basis)) + 
  geom_point(aes(numabs, loglik, group = group, colour=startknots)) + 
  theme_bw() + 
  geom_vline(aes(xintercept=9690), linetype=2) + 
  xlab('Number Qudrature Points') + 
  ylab('Log-Likelihood') +
  scale_color_discrete("Start Knots") + 
  scale_linetype_discrete(labels=c("Exponential", "Gaussian"))

png('results/ppconverge.png', res = 300,height=900, width=1500)
p1
dev.off()

p2 <- mutate(output, 
       group=paste0(Basis, startknots)) %>%
  group_by(group) %>% 
  mutate(scaledLL = scale(loglik)) %>%
  ggplot() + 
  geom_line(aes(numabs, scaledLL, group = group, colour=startknots, linetype=Basis)) + 
  geom_point(aes(numabs, scaledLL, group = group, colour=startknots)) + 
  theme_bw() + 
  geom_vline(aes(xintercept=9690), linetype=2) + 
  xlab('Number Qudrature Points') + 
  ylab('Scaled and Centered Log-Likelihood') +
  scale_color_discrete("Start Knots") + 
  scale_linetype_discrete(labels=c("Exponential", "Gaussian"))

png('results/ppconverge_scaled.png', res = 300,height=900, width=1500)
p2
dev.off()

