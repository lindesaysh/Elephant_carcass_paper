setwd("C:/Users/lass/OneDrive - University of St Andrews/papers/papers/CarcassSALSApaper/August2020")
devtools::load_all('C:/Users/lass/Documents/GitHub/MRSea')

require(dplyr)


set.seed(1)
res<- c()
rad<- FALSE # choose radii


analysisdat<-read.csv(file='data/analysisdat.csv')
require(rgdal)
proj4string = CRS("+proj=utm +zone=33 +units=km")

kg<-read.csv(file='data/kg.csv')

kseq<-seq(5, 60, by = 5)
sktime <- proc.time()
for(k in kseq){
  knotgrid<-getKnotgrid(kg, numKnots = k, plot=FALSE)
  save(knotgrid, file=paste('kg', k, '.RData', sep=''))
}
endktime <- proc.time() - sktime

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

require(MuMIn)
require(boot)
res<- c()
preds<- c()

starttime<-proc.time()

#for(k in c(50, 55, 60)){
for(k in kseq){
  print(k)
  
  load(paste('data/kg', k, '.RData', sep=''))
  
  for(r in 1:10){
    print(r)
    
    stime<-proc.time()
    
    knotid<-attributes(knotgrid)$points.selected
    
    r_seqG<-getRadiiChoices(10, distMats$dataDist,
                            basis="gaussian")
    r_seqE<-getRadiiChoices(10, distMats$dataDist,
                            basis="exponential")
    
    r_seqG_geo<-getRadiiChoices(10, distMats_Geo$dataDist,
                                basis="gaussian")
    r_seqE_geo<-getRadiiChoices(10, distMats_Geo$dataDist,
                                basis="exponential")
    
    
    
    b_G<-LRF.g(rep(r,times=k), distMats$dataDist, radii=r_seqG,
               aR=attr(knotgrid, "points.selected"))
    b_E<-LRF.e(rep(r,times=k), distMats$dataDist, radii=r_seqE,
               aR=attr(knotgrid, "points.selected"))
    
     b_G_geo<-LRF.g(rep(r,times=k), distMats_Geo$dataDist, 
                   radii=r_seqG_geo, aR=attr(knotgrid, "points.selected"))
    b_E_geo<-LRF.e(rep(r,times=k), distMats_Geo$dataDist,
                   radii=r_seqE_geo, aR=attr(knotgrid, "points.selected"))
    
   # 
    dattest_G<- data.frame(response=analysisdat$response,
                           pp.wts = analysisdat$pp.wts, b_G)
    dattest_E<- data.frame(response=analysisdat$response,
                           pp.wts = analysisdat$pp.wts, b_E)
    
    dattest_G_geo<- data.frame(response=analysisdat$response,
                               pp.wts = analysisdat$pp.wts, b_G_geo)
    dattest_E_geo<- data.frame(response=analysisdat$response,
                               pp.wts = analysisdat$pp.wts, b_E_geo)
    
    modFit_G<- try(glm(response/pp.wts ~ ., family='poisson', data=dattest_G, weights=pp.wts), silent = TRUE)
    modFit_E<- try(glm(response/pp.wts ~ ., family='poisson', data=dattest_E, weights=pp.wts), silent = TRUE)
   # 
    
    modFit_G_geo<- try(glm(response/pp.wts ~ ., family='poisson', data=dattest_G_geo, weights=pp.wts), silent = TRUE)
    modFit_E_geo<- try(glm(response/pp.wts ~ ., family='poisson', data=dattest_E_geo, weights=pp.wts), silent = TRUE)
    
    if(class(modFit_E)=="try-error"){
      AICc_E=NA
      BIC_E=NA
      predE<-NA
    }else{
      AICc_E=AICc(modFit_E)
      BIC_E=BIC(modFit_E)
      predE<-fitted(modFit_E)[analysisdat$response==0]
    }
   
    if(class(modFit_G)=="try-error"){
      AICc_G=NA
      BIC_G=NA
      predG<-NA
    }else{
      AICc_G=AICc(modFit_G)
      BIC_G=BIC(modFit_G)
      predG<-fitted(modFit_G)[analysisdat$response==0]
    }
    
    if(class(modFit_E_geo)=="try-error"){
      AICc_E_geo=NA
      BIC_E_geo=NA
      predE_geo<-NA
    }else{
      AICc_E_geo=AICc(modFit_E_geo)
      BIC_E_geo=BIC(modFit_E_geo)
      predE_geo<-fitted(modFit_E_geo)[analysisdat$response==0]
    }
    
    if(class(modFit_G_geo)=="try-error"){
      AICc_G_geo=NA
      BIC_G_geo=NA
      predG_geo<-NA
    }else{
      AICc_G_geo=AICc(modFit_G_geo)
      BIC_G_geo=BIC(modFit_G_geo)
      predG_geo<-fitted(modFit_G_geo)[analysisdat$response==0]
    }
    
    etime = (proc.time() - stime)[3]
    
    res<- rbind(res, data.frame(k, r, 
                                AICc_G,
                                AICc_E,
                                AICc_G_geo, 
                                AICc_E_geo, 
                                BIC_G,
                                BIC_E,
                                BIC_G_geo,
                                BIC_E_geo,
                                time=etime))
    write.csv(res, "results/resultsmodavg.csv", row.names=F)
    
    preds<- rbind(preds, data.frame(k, r, 
                                    x.pos=analysisdat$x.pos[analysisdat$response==0],
                                    y.pos=analysisdat$y.pos[analysisdat$response==0],
                                    predE, predE_geo, predG, predG_geo))
    write.csv(preds, "results/predsmodavg.csv", row.names=F)
    
  }
}

endtime = proc.time() - starttime

save.image("results/modelavgwkspace_180621.RData")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


res<-read.csv("results/resultsmodavg.csv")
rownames(res)<-NULL
preds<-read.csv("results/predsmodavg.csv")

res$diffAICc_G<- res$AICc_G-min(res$AICc_G, na.rm=TRUE)
res$diffAICc_E<- res$AICc_E-min(res$AICc_E, na.rm=TRUE)
res$diffAICc_G_geo<- res$AICc_G_geo-min(res$AICc_G_geo, na.rm=TRUE)
res$diffAICc_E_geo<- res$AICc_E_geo-min(res$AICc_E_geo, na.rm=TRUE)


res$w_E<- exp(-1/2*res$diffAICc_E)/sum(exp(-1/2*res$diffAICc_E), na.rm=TRUE)
res$w_G<- exp(-1/2*res$diffAICc_G)/sum(exp(-1/2*res$diffAICc_G), na.rm=TRUE)
res$w_E_geo<- exp(-1/2*res$diffAICc_E_geo)/sum(exp(-1/2*res$diffAICc_E_geo), na.rm=TRUE)
res$w_G_geo<- exp(-1/2*res$diffAICc_G_geo)/sum(exp(-1/2*res$diffAICc_G_geo), na.rm=TRUE)


res_E<-res[order(res$diffAICc_E),]
res_G<-res[order(res$diffAICc_G),]
res_E_geo<-res[order(res$diffAICc_E_geo),]
res_G_geo<-res[order(res$diffAICc_G_geo),]


head(res_E[,c('k', 'r', 'w_E', 'diffAICc_E')])
head(res_G[,c('k', 'r', 'w_G', 'diffAICc_G')])
head(res_E_geo[,c('k', 'r', 'w_E_geo', 'diffAICc_E_geo')])
head(res_G_geo[,c('k', 'r', 'w_G_geo', 'diffAICc_G_geo')])



makeAVGpreds<-function(result, scorename, preds, predname){
  result$delta<-as.vector(res[,scorename]-min(res[,scorename],na.rm=TRUE))
  avg.set<-which(result$delta<=10)
  res.avg<-result[avg.set,]
  res.avg$w<-AICc.weights(res.avg$delta)
  res.avg$keepid<-1
  preds$index<-1:nrow(preds)
  avg.preds<-merge(preds, res.avg)
  avg.preds<-arrange(avg.preds, k, r, index)
  avg.preds$wpreds<-avg.preds[,predname]*avg.preds$w
  wpreds.mat<-matrix(NA, 9690, nrow(res.avg))
  for(i in 1:nrow(res.avg)){
    predid<-which(avg.preds$k==res.avg$k[i] & avg.preds$r==res.avg$r[i])
    wpreds.mat[,i]<-avg.preds$wpreds[predid]
  }
  finalpred<-apply(wpreds.mat, 1, sum)
  return(finalpred)
}

AICc.weights<-function(delta){
  exp((-1/2)*delta)/sum(exp(-(1/2)*delta))
}


require(dplyr)
preds.exp.euc<-makeAVGpreds(res, "AICc_E", preds, "predE")
preds.gau.euc<-makeAVGpreds(res, "AICc_G", preds, "predG")
preds.exp.geo<-makeAVGpreds(res, "AICc_E_geo", preds, "predE_geo")
preds.gau.geo<-makeAVGpreds(res, "AICc_G_geo", preds, "predG_geo")

predgrid<-analysisdat[analysisdat$response==0,]
require(fields)
par(mfrow=c(2,2))
quilt.plot(predgrid$x.pos, predgrid$y.pos, preds.exp.euc, nrow=170, ncol=77, main="Exp, Euc", asp=1)
quilt.plot(predgrid$x.pos, predgrid$y.pos, preds.gau.euc, nrow=170, ncol=77, main="Gaussian, Euc", asp=1)
quilt.plot(predgrid$x.pos, predgrid$y.pos, preds.exp.geo, nrow=170, ncol=77, main="Exp, Geo", asp=1)
quilt.plot(predgrid$x.pos, predgrid$y.pos, preds.gau.geo, nrow=170, ncol=77, main="Gaussian, Geo", asp=1)


score<-"diffAICc_G"
avg.set<-which(res[,score]<=10)
res.avg<-res[avg.set,]
w<-AICc.weights(res.avg[,score])
modelout_Gaus_Euc<-data.frame(k=res.avg$k, r=res.avg$r, w=w)

score<-"diffAICc_E"
avg.set<-which(res[,score]<=10)
res.avg<-res[avg.set,]
w<-AICc.weights(res.avg[,score])
modelout_Exp_Euc<-data.frame(k=res.avg$k, r=res.avg$r, w=w)

score<-"diffAICc_G_geo"
avg.set<-which(res[,score]<=10)
res.avg<-res[avg.set,]
w<-AICc.weights(res.avg[,score])
modelout_Gaus_Geo<-data.frame(k=res.avg$k, r=res.avg$r, w=w)

score<-"diffAICc_E_geo"
avg.set<-which(res[,score]<=10)
res.avg<-res[avg.set,]
w<-AICc.weights(res.avg[,score])
modelout_Exp_Geo<-data.frame(k=res.avg$k, r=res.avg$r, w=w)

outmodels<-list(expeuc=modelout_Exp_Euc, expgeo=modelout_Exp_Geo, gauseuc=modelout_Gaus_Euc, gausgeo=modelout_Gaus_Geo)
save(outmodels, file='results/resultsmodelavg_models.RData')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
require(MuMIn)
analysisdat$y <- analysisdat$response/analysisdat$pp.wts
scores<-NULL
for(m in 1:4){
 #m=1
  #for(cvs in 1:100){
  #  myfold <- allfolds[[cvs]]
  #  cv=0

    #for(f in 1:10){
      res<-NULL
      preds<-NULL

      res.in<-outmodels[[m]]

      if(m==1 | m==2){
        basis='exponential'
      }else{
        basis='gaussian'
      }

      if(m==1 | m==3){
        disttype='euc'
        dMat <- distMats
      }else{
        disttype='geo'
        dMat <- distMats_Geo
      }

      for(i in 1:nrow(res.in)){

        k<-res.in[i,'k']
        r<-res.in[i,'r']

        load(paste('data/kg', k, '.RData', sep=''))

        knotid<-attributes(knotgrid)$points.selected

        r_seq<-getRadiiChoices(10, dMat$dataDist, basis=basis)

        if(basis=='gaussian'){
          b<-LRF.g(rep(r,times=k), dMat$dataDist, radii=r_seq, aR=attr(knotgrid, "points.selected"))
        }else{
          b<-LRF.e(rep(r,times=k), dMat$dataDist, radii=r_seq, aR=attr(knotgrid, "points.selected"))

        }

        dattest<- data.frame(response=analysisdat$response, pp.wts = analysisdat$pp.wts, b)
        modFit<- glm(response/pp.wts ~ ., family='poisson', data=dattest, weights=pp.wts)

        # get predicted intensity for all points
        pred =  fitted(modFit)

        res<- rbind(res, data.frame(k, r, AICc=AICc(modFit),                             BIC=BIC(modFit), w=res.in$w[i]))

        preds<- cbind(preds,  pred)

      } # end i loop (single model run)

      if(res$w==1){
        preds.avg <- preds
      }else{
        preds.avg<-preds%*%res$w  
      }
      
      # calculate loglik
      LL<- sum(analysisdat$pp.wts * ((analysisdat$y * log(preds.avg)) - preds.avg))
      
      mscore<-data.frame(Model=m,
                         Basis=basis,
                         DistanceType=disttype,
                         No.Models=nrow(res.in),
                         LogLik = LL)
      scores<-rbind(scores, mscore)

} # end m loop

scores<-arrange(scores, desc(LogLik))
write.csv(scores, file=paste('results/ModAvgScores.csv', sep=''), row.names = FALSE)

# ~~~~~~~~~~
plotdat <- filter(analysisdat, response==0)
plotdat$preds.exp.euc<-preds.exp.euc
plotdat$preds.gau.euc<-preds.gau.euc
plotdat$preds.exp.geo<-preds.exp.geo
plotdat$preds.gau.geo<-preds.gau.geo

require(ggplot2)
p1<-ggplot() + 
  geom_tile(data=plotdat, aes(x=x.pos, y=y.pos, fill=preds.exp.euc), height=2, width=2) +
  scale_fill_distiller(palette = "Spectral",name="Intensity", limits=c(0, 0.3)) + 
  coord_equal() + 
  theme_bw() + ggtitle("Exponential, Euclidean") +
  geom_point(data=filter(analysisdat, response==1), aes(x=x.pos, y=y.pos), shape=1)

p2<-ggplot() + 
  geom_tile(data=plotdat, aes(x=x.pos, y=y.pos, fill=preds.gau.euc), height=2, width=2) +
  scale_fill_distiller(palette = "Spectral",name="Intensity", limits=c(0, 0.3)) + 
  coord_equal() + 
  theme_bw() + ggtitle("Gaussian, Euclidean") +
  geom_point(data=filter(analysisdat, response==1), aes(x=x.pos, y=y.pos), shape=1)

p3<-ggplot() + 
  geom_tile(data=plotdat, aes(x=x.pos, y=y.pos, fill=preds.exp.geo), height=2, width=2) +
  scale_fill_distiller(palette = "Spectral",name="Intensity", limits=c(0, 0.3)) + 
  coord_equal() + 
  theme_bw() + ggtitle("Exponential, Geodesic") +
  geom_point(data=filter(analysisdat, response==1), aes(x=x.pos, y=y.pos), shape=1)

p4<-ggplot() + 
  geom_tile(data=plotdat, aes(x=x.pos, y=y.pos, fill=preds.gau.geo), height=2, width=2) +
  scale_fill_distiller(palette = "Spectral",name="Intensity", limits=c(0, 0.56)) + 
  coord_equal() + 
  theme_bw() + ggtitle("Gaussian, Geodesic") +
  geom_point(data=filter(analysisdat, response==1), aes(x=x.pos, y=y.pos), shape=1)

require(gridExtra)
grid.arrange(p3, p4, p1, p2, nrow=2)

png('results/modelavgpredplots.png', height=600, width=1000)
grid.arrange(p3, p4, p1, p2, nrow=2)
dev.off()

write.csv(plotdat, 'results/modelavgplotdata_140621.csv', row.names=FALSE)



# getting knot locations for final model with radii.

res.in<-outmodels[[2]]
basis='exponential'
disttype='geo'
dMat <- distMats_Geo
r_seq<-getRadiiChoices(10, dMat$dataDist, basis=basis)
kgnew<-NULL

for(i in 1:nrow(res.in)){
  
  k<-res.in[i,'k']
  r<-res.in[i,'r']
  
  load(paste('data/kg', k, '.RData', sep=''))
  
  
  knotid<-attributes(knotgrid)$points.selected
  b<-LRF.e(rep(r,times=k), dMat$dataDist, radii=r_seq, aR=attr(knotgrid, "points.selected"))
  dattest<- data.frame(response=analysisdat$response, pp.wts = analysisdat$pp.wts, b)
  modFit<- glm(response/pp.wts ~ ., family='poisson', data=dattest, weights=pp.wts)
  
  kgnew<-rbind(kgnew, data.frame(knotgrid, 
                                 radii = r_seq[r],
                                 coefsign=as.factor(ifelse(coef(modFit)[-1]<0, 0, 1))))
  
} # end i loop (single model run)


library(viridis)
largepan<-read.csv('data/largepan.csv')[,2:3]/1000
fence<-read.csv('data/fencepoly.csv')/1000

modavgknots<-ggplot() +
  geom_point(data=kgnew,
             aes(x=x.pos, y=y.pos, colour=coefsign, size=radii), alpha=1/5) +
  theme_bw() + coord_equal() +
  labs(x='Easting (Km)', y='Northing (Km)') +
  geom_point(data=filter(analysisdat, response==1), aes(x=x.pos, y=y.pos), alpha=1/5, size=0.5) +
  geom_polygon(data=data.frame(largepan), aes(x=x, y=y), fill='lightblue')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  geom_polygon(data=fence, aes(x=x.pos, y=y.pos), fill=NA, colour='black') +
  scale_size(guide=FALSE) +
  scale_color_viridis(discrete=TRUE) +
  guides(colour="none") 

png('results/modavgknots.png', height=900, width=1500, res=300)
modavgknots
dev.off()




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

