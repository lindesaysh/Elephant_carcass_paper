setwd("C:/Users/lass/OneDrive - University of St Andrews/papers/papers/CarcassSALSApaper/August2020")
devtools::load_all('C:/Users/lass/Documents/GitHub/MRSea')

require(dplyr)
#require(MRSea)
require(MuMIn)

set.seed(1)
res<- c()
rad<- FALSE # choose radii

# # load in pan data (2.5km internal buffer)
# largepan<-read.csv('data/largepan.csv')
# largepan<-largepan/1000
# largepan.list<-list(x=largepan$x, y=largepan$y)
# # ENP boundary plus 20km
# fence<-read.csv('data/databoundary.csv')
# fence<-fence[nrow(fence):1,]/1000
# fence.list<-list(x=fence$x.pos, y=fence$y.pos)

# # read in points data
# datcomb<- read.csv("data/datcomb.csv", header=TRUE)
# #subset just to get the presences
# dat.pres<-datcomb[datcomb$response==1,]
# dat<-datcomb[datcomb$response==1,c('x.pos', 'y.pos')]/1000
# 
# # create boundary listing so that quad points created in correct area
# Z<-spatstat::owin(poly=list(fence.list, largepan.list))
# P<-spatstat::ppp(dat$x.pos, dat$y.pos, poly=list(fence.list, largepan.list))
# dummyQpts<-spatstat::quadscheme(P, eps=2)$dummy
# dpts<-data.frame(x.pos=dummyQpts$x, y.pos=dummyQpts$y)
# 
# analysisdat<-rbind(data.frame(dat, response=1), data.frame(dpts, response=0))
# p.wt<-rep(1e-6, nrow(analysisdat))
# p.wt[analysisdat$response==0]<-summary(dummyQpts)$window$area/sum(analysisdat$response==0)
# # must have column in data called pp.wts
# analysisdat$pp.wts<-p.wt
# write.csv(analysisdat, file='data/analysisdat.csv', row.names=FALSE)

analysisdat<-read.csv(file='data/analysisdat.csv')
require(rgdal)
proj4string = CRS("+proj=utm +zone=33 +units=km")


kg<-read.csv(file='data/kg.csv')

distMats<-makeDists(cbind(analysisdat$x.pos, analysisdat$y.pos), kg)

# load geodesic matrix and reduce to remove duplicated knots
load('data/knot2knotDist_geo.Robj')
load('data/data2knotDist_geo.Robj')


distMats_Geo<-list(dataDist= d2k, knotDist=k2k)
rm(d2k, k2k)

require(fields)
par(mfrow=c(2,1))
quilt.plot( analysisdat$x.pos, analysisdat$y.pos, distMats$dataDist[,1], asp=1, nrow=170, ncol=77)
quilt.plot( analysisdat$x.pos, analysisdat$y.pos, distMats_Geo$dataDist[,1], asp=1, nrow=170, ncol=77)

analysisdat$y <- analysisdat$response/analysisdat$pp.wts

initialModel<- glm(response/pp.wts ~ 1 , family='poisson', data=analysisdat, weights=pp.wts)

type=c('EucE', 'EucG', 'GeoE', 'GeoG')

listoftypeandk<- data.frame(type=rep(type,4), sk=rep(seq(65,70, by=5), each=2))

# load folddata
allfolds<-readRDS('data/folds_10_100.rds')


require(parallel)
cat('Code in parallel')
nCores<-5
computerCores <- getOption("cl.cores", detectCores())
if(nCores>(computerCores-2)){nCores<-(computerCores-2)}
myCluster <- makeCluster(nCores) ; myCluster # initialise the cluster
clusterExport(myCluster, ls(), envir=environment()) # export all objects to each cluster
# export directory and functions to each cluster
clusterEvalQ(myCluster, {
  require(splines)
  require(dplyr)
  devtools::load_all('C:/Users/lass/Documents/GitHub/MRSea')
  require(MuMIn)
})

# only do parametric bootsrap if no data re-sampling and no nhats provided
Routputs<-pbapply::pblapply(cl=myCluster, X=1:nrow(listoftypeandk), FUN=function(i){
  
  sk<-listoftypeandk[i,2]
  type<-listoftypeandk[i,1]
  
  salsa2dlist<-list(fitnessMeasure = 'BIC', 
                    knotgrid = kg, 
                    startKnots=sk,
                    minKnots=2, 
                    maxKnots=100, 
                    #cv.opts=list(K=5, cv.gamMRSea.seed=1),
                    modelType='pointProcess')
  
  starttime<-proc.time()
  if(type=='EucG'){
    salsa2dOutput<-runSALSA2D(initialModel, salsa2dlist, d2k=distMats$dataDist, k2k=distMats$knotDist, splineParams=NULL,suppress.printout = FALSE, basis='gaussian', chooserad = rad)
  }
  
  if(type=='GeoG'){
    salsa2dOutput<-runSALSA2D(initialModel, salsa2dlist, d2k=distMats_Geo$dataDist, k2k=distMats_Geo$knotDist, splineParams=NULL,suppress.printout = FALSE,basis='gaussian',  chooserad=rad)
  }
  
  if(type=='EucE'){
    salsa2dOutput<-runSALSA2D(initialModel, salsa2dlist, d2k=distMats$dataDist, k2k=distMats$knotDist, splineParams=NULL,suppress.printout = FALSE, basis='exponential', chooserad = rad)
  }
  
  if(type=='GeoE'){
    salsa2dOutput<-runSALSA2D(initialModel, salsa2dlist, d2k=distMats_Geo$dataDist, k2k=distMats_Geo$knotDist, splineParams=NULL,suppress.printout = FALSE,basis='exponential',  chooserad=rad)
  }
  endtime<-proc.time() - starttime
  
  bestModel<- salsa2dOutput$bestModel
  splist<-bestModel$splineParams
  

  ll<-sum(analysisdat$pp.wts * (((analysisdat$response/analysisdat$pp.wts) * log(fitted(bestModel))) - fitted(bestModel)))

  res<- data.frame(StartKnots=sk, 
                   type=type, 
                   noKnots=(length(coef(bestModel))-1),
                   ll=ll,
                   AICc=AICc(bestModel), 
                   BIC=BIC(bestModel),
                   timing_min =  as.numeric(endtime[3]/60))
  
  save(res, file=paste('results/res_bic_pp', sk, type, '_140621.RData' , sep=''), compress='bzip2')
  
  return(list(res=res, splist=splist))
})

save(Routputs, file='results/cressmrsea11062121_pp3.RData', compress='bzip2')

stopCluster(myCluster)
# 
# cat("Creating Outputs...")
# 

res1 = data.frame(t(sapply(Routputs, '[[','res')))
typedf<-data.frame(type.num=c(1:4), type=c('EucE', 'EucG', 'GeoE', 'GeoG'))

fulldat<- res1
fulldat$index<-1:nrow(fulldat)
fulldat$StartKnots<-as.numeric(fulldat$StartKnots)
fulldat$type<-as.character(fulldat$type)
fulldat$noKnots<-as.numeric(fulldat$noKnots)
fulldat$ll<-as.numeric(fulldat$ll)
fulldat$LCL<-as.numeric(fulldat$LCL)
fulldat$UCL<-as.numeric(fulldat$UCL)
fulldat$AICc<-as.numeric(fulldat$AICc)
fulldat$BIC<-as.numeric(fulldat$BIC)
fulldat$timing_min<-as.numeric(fulldat$timing_min)
fulldat <- left_join(typedf, fulldat)
# 
# 
# 
 splistall = sapply(Routputs, '[[','splist')
# # 

# require(dplyr)
# 

fulldat <- fulldat %>% mutate(mybic = (-2*ll) + (log(nrow(analysisdat)) * noKnots)) 
 
arrange(fulldat, desc(ll))
str(fulldat)
# 
require(ggplot2)
ggplot(fulldat) + 
   geom_line(aes(x=StartKnots, y=-ll, group=type, colour=type)) + theme_bw() + 
  ylab('Negative Log Likelihood')

# 
 ggplot(fulldat) + 
   geom_line(aes(x=noKnots, y=-ll, colour=type, group=type)) + theme_bw() +
   ylab('Negative Log Likelihood')
# 

modid<-arrange(fulldat, desc(ll))$index[1]
modelspec<-filter(fulldat, index==modid)
modelspec
basis_ll<- LRF.e(radiusIndices = splistall[[modid]]$radiusIndices, dists = splistall[[modid]]$dist, radii = splistall[[modid]]$radii, aR = splistall[[modid]]$knotPos)
  
mymodll<- gamMRSea(response/pp.wts ~ basis_ll , family='poisson', data=analysisdat, weights=pp.wts)

require(ggplot2)
analysisdat$fitsll<-fitted(mymodll)
p1<-ggplot() + 
  geom_tile(data=filter(analysisdat, response==0), aes(x=x.pos, y=y.pos, fill=fitsll), height=2, width=2) +
scale_fill_distiller(palette = "Spectral",name="Intensity", limits=c(0, 0.81)) + 
  coord_equal() + 
  theme_bw() + ggtitle(paste0("Best LL: Type = ", modelspec$type, ", No. Knots = ", modelspec$noKnots)) +
  geom_point(data=filter(analysisdat, response==1), aes(x=x.pos, y=y.pos), shape=1)

# ~~~~~~~~~

modid2<-arrange(fulldat, AICc)$index[1]
modelspec<-filter(fulldat, index==modid2)
basis_aic<- LRF.e(radiusIndices = splistall[[modid2]]$radiusIndices, 
                  dists = splistall[[modid2]]$dist, 
                  radii = splistall[[modid2]]$radii, 
                  aR = splistall[[modid2]]$knotPos)

mymodaic<- gamMRSea(response/pp.wts ~ basis_aic , family='poisson', data=analysisdat, weights=pp.wts)

analysisdat$fitsaic<-fitted(mymodaic)
p2<-ggplot() + 
  geom_tile(data=filter(analysisdat, response==0), aes(x=x.pos, y=y.pos, fill=fitsaic), height=2, width=2) +
  scale_fill_distiller(palette = "Spectral",name="Intensity", limits=c(0, 0.81)) + 
  coord_equal() + 
  theme_bw() + ggtitle(paste0("Best AIC: Type = ", modelspec$type, ", No. Knots = ", modelspec$noKnots))  +
  geom_point(data=filter(analysisdat, response==1), aes(x=x.pos, y=y.pos), shape=1)

# ~~~~~~~~~~~~~

modid3<-arrange(fulldat, BIC)$index[1]
modelspec<-filter(fulldat, index==modid3)
modelspec
basis_bic<- LRF.e(radiusIndices = splistall[[modid3]]$radiusIndices, 
                  dists = splistall[[modid3]]$dist, 
                  radii = splistall[[modid3]]$radii, 
                  aR = splistall[[modid3]]$knotPos)

mymodbic<- gamMRSea(response/pp.wts ~ basis_bic , family='poisson', data=analysisdat, weights=pp.wts)

analysisdat$fitsbic<-fitted(mymodbic)
p3<-ggplot() + 
  geom_tile(data=filter(analysisdat, response==0), aes(x=x.pos, y=y.pos, fill=fitsbic), height=2, width=2) +
  scale_fill_distiller(palette = "Spectral",name="Intensity", limits=c(0, 0.81)) + 
  coord_equal() + 
  theme_bw() + ggtitle(paste0("Best BIC: Type = ", modelspec$type, ", No. Knots = ", modelspec$noKnots))  +
  geom_point(data=filter(analysisdat, response==1), aes(x=x.pos, y=y.pos), shape=1)

# ~~~~~~~~~~~~



gridExtra::grid.arrange(p1, p2, p3, ncol=1)
