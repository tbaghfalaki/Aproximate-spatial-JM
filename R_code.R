rm(list=ls())
library(mvtnorm)
library(INLA);library(extraDistr);library(splines);
library(SpatialEpi)
library(rgdal)
library(cplm)
library(survival)
library(ggplot2)
#library(dplyr)
library(ggfortify)
setwd("/Users/taban/Downloads/SJM")

###########
joint=list()
longdat <-list()# joint$longitudinal
survdat <-list()# joint$survival
set.seed(1)
Data=read.table("Data_new.txt",header=TRUE)
names(Data)
attach(Data)
head(Data)

library(CorrMixed)
Spaghetti.Plot(Dataset = Data, Outcome = cd4, Id=patient, Time = timecd4,
               Add.Mean = FALSE, Add.Median = FALSE,Col="gray4")

n2=length(unique(patient))
Nobs=length(patient)
########################
Fox_Lattice <- readOGR("Map/Brasil.shp")
plot(Fox_Lattice)

ZONE_CODE=c()
Surv0=Surv1=Cen=Gender=Race=Age=POI=c()
m=c()
TIME=CD4=matrix(NA,n2,17)
r=0
for(k in unique(patient)){
  r=r+1
  m[r]=length(states[patient==k])
  ZONE_CODE[r]=as.numeric(states[patient==k][1])
  Surv0[r]=as.numeric(dtini[patient==k][1])
  Surv1[r]=as.numeric(dtend[patient==k][1])
  Cen[r]=as.numeric(max(Censure[patient==k]))
  Gender[r]=as.numeric(max(sex[patient==k]))
  Race[r]=as.numeric(max(race[patient==k]))
  Age[r]=as.numeric(max(age50[patient==k]))
  POI[r]=as.numeric(max(prevoi[patient==k]))
  TIME[r,1:length(timecd4[patient==k])]=timecd4[patient==k]
  CD4[r,1:length(timecd4[patient==k])]=cd4[patient==k]
}
Surv=Surv1-Surv0
ID=rep(c(1:n2),times=m)
table(ZONE_CODE)
require(spdep)  # a package that can tabulate contiguity in spatial objects, i.e., the state of bordering or being in contact with something
Lattice_Temp <- poly2nb(Fox_Lattice)  # construct the neighbour list
nb2INLA("Lattice.graph", Lattice_Temp) # create the adjacency matrix in INLA format
Lattice.adj <- paste(getwd(),"/Lattice.graph",sep="") # name the object

inla.setOption(scale.model.default = F)
H <- inla.read.graph(filename = "Lattice.graph")  # and save it as a graph

# Plot adjacency matrix
image(inla.graph2matrix(H), xlab = "", ylab = "")

knots=quantile(timecd4/365,prob=c(.25,.50,.75)) # B-spline
SS=bsp(timecd4/365, k=15)$Z

##########
longdat1=survdat1=list()
longdat1$y=sqrt(cd4)
longdat1$TIME =timecd4/365#/max(timecd4)
longdat1$Gender =sex
longdat1$Race =race
longdat1$Age =age50
longdat1$POI =prevoi
longdat1$ID=ID

#####################################
survdat1$CENSOR=Cen
survdat1$SURVTIME=Surv/365#/max(timecd4)
survdat1$Gender=Gender
survdat1$Race =Race
survdat1$Age =Age
survdat1$POI =POI

longdat=longdat1
survdat=survdat1
n1 <- length(longdat1$y)
n2 <- length(survdat1$CENSOR)


################## ################## ################## ###############
y.long <- c(longdat$y)
y.surv <- inla.surv(time = c(rep(NA, n1), survdat$SURVTIME), event = c(rep(NA, n1), survdat$CENSOR))
Yjoint <- list(y.long, y.surv)

km_trt_fit <- survfit(Surv(survdat1$SURVTIME,survdat1$CENSOR) ~ 1)
ggplot2::autoplot(km_trt_fit)


linear.covariate <- data.frame(mu = as.factor(c(rep(1, n1), rep(2, n2))), 
                               l.TIME = c(longdat$TIME, rep(0, n2)), 
                               l.Age = c(longdat$Age,  rep(0, n2)), 
                               l.Gender = c(longdat$Gender,  rep(0, n2)),
                               l.POI = c(longdat$POI,  rep(0, n2)),
                               s.Age = c(rep(0, n1),  survdat$Age),
                               s.Gender  = c(rep(0, n1),  survdat$Gender ),
                               s.POI = c(rep(0, n1),  survdat$POI)
)

SS2=matrix(0,n2,length(SS[1,]))

linear.covariate$l.SS=rbind(SS,SS2)


random.covariate <- list(U11 = c(longdat1$ID, rep(NA, n2)),
                         U21 = c(longdat1$ID+n2, rep(NA, n2)), 
                         U12 = c(rep(NA, n1), 1:n2),
                         U22 = c(rep(NA, n1),  n2+(1:n2)),
                         ZONE_CODE=c(rep(NA, n1), ZONE_CODE))



joint.data <- c(linear.covariate,random.covariate)
joint.data$Y <- Yjoint


formula = Y~ mu-1 +  l.Age +l.Gender+l.POI+ l.SS+s.Age +s.Gender+s.POI+ f(U11, model="iid2d",n=2*n2) +
  f(U21, l.TIME,  copy="U11") +
  f(U12, copy = "U11", fixed = FALSE) +
  f(U22, copy = "U11", fixed = FALSE)+ f(ZONE_CODE, model = "besag",       # spatial effect: ZONE_CODE is a numeric identifier for each area in the lattice  (does not work with factors)
                                         graph = Lattice.adj)  

#joint.inla <- inla(formula, family = c("gaussian","exponential.surv"),
#              data = joint.data, control.compute=list(dic=TRUE,waic = TRUE,cpo = TRUE),control.family = list(list(),list()))
#summary(joint.inla)
#LPML1 <- sum(log(joint.inla$cpo$cpo))
##########################################
joint.inla<- inla(formula, family = c("gaussian","weibullsurv"),
                  data = joint.data, control.compute=list(dic=TRUE,waic = TRUE,cpo = TRUE),control.family = list(list(),list()))
#summary(joint.inla)


LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1

##################    #####
predl=joint.inla$summary.fitted.values$mean[1:n1]
preds=joint.inla$summary.fitted.values$mean[(1+n1):(n2+n1)]



predL=matrix(NA,n2,17)
r=0
for(k in unique(patient)){
  r=r+1
  predL[r,1:length(timecd4[patient==k])]=predl[patient==k]
  
}

patients <- c(12, 230,489)
par(mfrow = c(3, 2))
par(mar = c(4, 4, 4, 4))
alpha <- joint.inla$summary.hyperpar$mean[1]
j_est <- joint.inla$summary.hyperpar$mean[2:6]
f_est <- joint.inla$summary.fixed$mean

for (patientnr in patients) {
  dataHi <- Data[Data$patient == patientnr, ]
  plot(TIME[patientnr,], sqrt(CD4[patientnr,]),
       ylab = "CD4 count", xlab = "Time (months)", type = "l",
       xlim=c(0,1820),ylim=c(0,45),
       main = paste("CD4 trajectory - patient", patientnr))
  
  lines(TIME[patientnr,],predL[patientnr,], col = "blue", lty = 2)
  legend("topright",c("Observed","Prediction"),lty=1:2,col=c(1,"blue"))
  
  
  r= joint.inla$summary.hyperpar$mean[2]
  
  X=c(1,linear.covariate$s.Age[n1+patientnr],linear.covariate$s.Gender[n1+patientnr],linear.covariate$s.POI[n1+patientnr])
  llll=length(joint.inla$summary.fixed$mean)
  lambda=c(joint.inla$summary.fixed$mean[2],joint.inla$summary.fixed$mean[(llll-2):llll])%*%X
  
  
  pred_time <- seq(0, 5, by = 0.01)
  bbbb=(as.numeric(pred_time)/as.numeric(-lambda))^r
  plot(pred_time, exp(-bbbb),
       type = "l", ylab = "Survival probability", xlab = "Time (months)",
       main = paste("Survival curve - patient", patientnr))
  abline(h = 0.5, col = "red")
}

######################
par(mfrow=c(1,2))
plot(joint.inla$summary.fitted.values$mean[1:n1],y.long,
     ylab="Observed",xlab="Fitted",pch=19,main="(a)")
lines(1:42,1:42,col=2,lwd=2)

#plot(joint.inla$summary.fitted.values$mean[1:n1]-y.long,pch=19,
#   ylab="Residual",xlab="Index")
sigma2= joint.inla$summary.hyperpar$mean[1]^-1
d11= joint.inla$summary.hyperpar$mean[2]^-1
d22= joint.inla$summary.hyperpar$mean[3]^-1
sigmau= joint.inla$summary.hyperpar$mean[5]^-1

V=sigma2+d11+d22*longdat$TIME^2
Z=(joint.inla$summary.fitted.values$mean[1:n1]-y.long)/sqrt(V)
plot(Z,pch=19,ylab="Standardized Residual",xlab="Index",main="(b)")

#qqnorm(Z, pch=19)
#qqline(Z, col =2, lwd = 2
############################# 
####### cox -snell res
#############################
r1=joint.inla$summary.hyperpar$mean[7]
r2=joint.inla$summary.hyperpar$mean[8]
D=matrix(c(joint.inla$summary.hyperpar$mean[2]^-1,
           joint.inla$summary.hyperpar$mean[2]^-1*joint.inla$summary.hyperpar$mean[3]^-1*joint.inla$summary.hyperpar$mean[4],
           joint.inla$summary.hyperpar$mean[2]^-1*joint.inla$summary.hyperpar$mean[3]^-1*joint.inla$summary.hyperpar$mean[4],
           joint.inla$summary.hyperpar$mean[3]^-1),2,2)

a=1
D1=a*diag(1,2)+(1-a)*D
b=rmvnorm(n2,rep(0,2),D1)

tau=joint.inla$summary.hyperpar$mean[5]
Q=matrix(tau,H$n,H$n)
for(kk in 1:H$n){
  Q[kk,kk]=tau*(H$nnbs[kk])
}
a=0.01
Q1=a*diag(1,H$n)+(1-a)*Q
u=c()
u=rmvnorm(1,rep(0,H$n),solve(Q1))

lambda=c()
for(patientnr in 1:n2){ 
  X=c(1,linear.covariate$s.Age[n1+patientnr],linear.covariate$s.Gender[n1+patientnr],linear.covariate$s.POI[n1+patientnr])
  llll=length(joint.inla$summary.fixed$mean)
  lambda[patientnr]=c(joint.inla$summary.fixed$mean[2],joint.inla$summary.fixed$mean[(llll-2):llll])%*%X+
    r1*b[patientnr,1]+r2*b[patientnr,2]+u[ZONE_CODE[patientnr]]
}


#lambda=joint.inla$summary.linear.predictor$mean[(n1+1):(n1+n2)]

par(mfrow=c(1,1))

#r=1
surv2=1-pweibull(survdat1$SURVTIME,r,exp(-lambda/r))
#surv2=1-pexp(survdat1$SURVTIME,exp(lambda))

res.coxsnell <- -log(surv2)#apply(mui.t, 2, mean)*(data$tee) # tee is the time of observation
# either by death or by censure;
# in a Weibull(alpha,mut[i]) the Cox snell Residuals are r_ci(t)= mu_i(t)*t^alpha
km <- survfit(Surv(res.coxsnell,survdat1$CENSOR,type='right')~1)
plot(km, xlab="Cox-Snell Residuals", xlim=c(0,1.2),main="Survival Function of Cox-Snell Residuals",
     mark.time=F, conf.int=F,lwd=2)
#Plot the survival function for an exponential distribution $Exp(1)$
x <- seq(0,10,0.01)
lines(x,1-pexp(x),"l",lwd=1,lty=2)



joint.inla$summary.fixed[,c(1:3,5)]
joint.inla$summary.hyperpar[,c(1:3,5)]


#################################################
formula = Y~ mu-1 +  l.Age +l.Gender+l.POI+ l.SS+s.Age +s.Gender+s.POI+ f(U11, model="iid2d",n=2*n2) +
  f(U21, l.TIME,  copy="U11") +
  f(U12, copy = "U11", fixed = FALSE) +
  f(U22, copy = "U11", fixed = FALSE)
joint.inla<- inla(formula, family = c("gaussian","weibullsurv"),
                  data = joint.data, control.compute=list(dic=TRUE,waic = TRUE,cpo = TRUE),control.family = list(list(),list()))
summary(joint.inla)
LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1
#################################################

formula = Y~ mu-1 +  l.Age +l.Gender+l.POI+ l.SS+s.Age +s.Gender+s.POI+ f(U11, model="iid2d",n=2*n2) +
  f(U21, l.TIME,  copy="U11") + f(ZONE_CODE, model = "besag",       # spatial effect: ZONE_CODE is a numeric identifier for each area in the lattice  (does not work with factors)
                                  graph = Lattice.adj)  

joint.inla<- inla(formula, family = c("gaussian","weibullsurv"),
                  data = joint.data, control.compute=list(dic=TRUE,waic = TRUE,cpo = TRUE),control.family = list(list(),list()))
summary(joint.inla)
LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1

#################################################

formula = Y~ mu-1 +  l.Age +l.Gender+l.POI+ l.SS+s.Age +s.Gender+s.POI+ f(U11, model="iid2d",n=2*n2) +
  f(U21, l.TIME,  copy="U11") 

joint.inla<- inla(formula, family = c("gaussian","weibullsurv"),
                  data = joint.data, control.compute=list(dic=TRUE,waic = TRUE,cpo = TRUE),control.family = list(list(),list()))
summary(joint.inla)
LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1

#################################################
formula = Y~ mu-1 +  l.Age +l.Gender+l.POI+ l.SS+s.Age +s.Gender+s.POI

joint.inla<- inla(formula, family = c("gaussian","weibullsurv"),
                  data = joint.data, control.compute=list(dic=TRUE,waic = TRUE,cpo = TRUE),control.family = list(list(),list()))
summary(joint.inla)
LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1
################################################# b1 0 0 
formula = Y~ mu-1 +  l.Age +l.Gender+l.POI+ l.SS+s.Age +s.Gender+s.POI+ f(U11) 

joint.inla<- inla(formula, family = c("gaussian","weibullsurv"),
                  data = joint.data, control.compute=list(dic=TRUE,waic = TRUE,cpo = TRUE),control.family = list(list(),list()))
summary(joint.inla)
LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1
#################################################
formula = Y~ mu-1 +  l.Age +l.Gender+l.POI+ l.SS+s.Age +s.Gender+s.POI+ f(U11) +
  f(U12, copy = "U11", fixed = FALSE) 

joint.inla<- inla(formula, family = c("gaussian","weibullsurv"),
                  data = joint.data, control.compute=list(dic=TRUE,waic = TRUE,cpo = TRUE),control.family = list(list(),list()))
summary(joint.inla)
LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1
#################################################
formula = Y~ mu-1 +  l.Age +l.Gender+l.POI+ l.SS+s.Age +s.Gender+s.POI+ f(U11) +
  f(U12, copy = "U11", fixed = FALSE)  + f(ZONE_CODE, model = "besag",      
                                           graph = Lattice.adj) 

joint.inla<- inla(formula, family = c("gaussian","weibullsurv"),
                  data = joint.data, control.compute=list(dic=TRUE,waic = TRUE,cpo = TRUE),control.family = list(list(),list()))
summary(joint.inla)
LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1
#################################################
formula = Y~ mu-1 +  l.Age +l.Gender+l.POI+ l.SS+s.Age +s.Gender+s.POI+ f(ZONE_CODE, model = "besag",      
                                                                          graph = Lattice.adj)  

joint.inla<- inla(formula, family = c("gaussian","weibullsurv"),
                  data = joint.data, control.compute=list(dic=TRUE,waic = TRUE,cpo = TRUE),control.family = list(list(),list()))
summary(joint.inla)
LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1

################################################# b1+b2t, gamma b1
formula = Y~ mu-1 +  l.Age +l.Gender+l.POI+ l.SS+s.Age +s.Gender+s.POI+ f(U11, model="iid2d",n=2*n2) +
  f(U21, l.TIME,  copy="U11") +
  f(U12, copy = "U11", fixed = FALSE) + f(ZONE_CODE, model = "besag",      
                                          graph = Lattice.adj)  

joint.inla<- inla(formula, family = c("gaussian","weibullsurv"),
                  data = joint.data, control.compute=list(dic=TRUE,waic = TRUE,cpo = TRUE),control.family = list(list(),list()))
summary(joint.inla)
LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1
################################################# b1+b2t, gamma b2
formula = Y~ mu-1 +  l.Age +l.Gender+l.POI+ l.SS+s.Age +s.Gender+s.POI+ f(U11, model="iid2d",n=2*n2) +
  f(U21, l.TIME,  copy="U11") +
  f(U22, copy = "U11", fixed = FALSE) + f(ZONE_CODE, model = "besag",      
                                          graph = Lattice.adj)  

joint.inla<- inla(formula, family = c("gaussian","weibullsurv"),
                  data = joint.data, control.compute=list(dic=TRUE,waic = TRUE,cpo = TRUE),control.family = list(list(),list()))
summary(joint.inla)


LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1

################################################# b1+b2t, gamma b1
formula = Y~ mu-1 +  l.Age +l.Gender+l.POI+ l.SS+s.Age +s.Gender+s.POI+ f(U11, model="iid2d",n=2*n2) +
  f(U21, l.TIME,  copy="U11") +
  f(U12, copy = "U11", fixed = FALSE) 

joint.inla<- inla(formula, family = c("gaussian","weibullsurv"),
                  data = joint.data, control.compute=list(dic=TRUE,waic = TRUE,cpo = TRUE),control.family = list(list(),list()))
summary(joint.inla)
LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1
################################################# b1+b2t, gamma b2
formula = Y~ mu-1 +  l.Age +l.Gender+l.POI+ l.SS+s.Age +s.Gender+s.POI+ f(U11, model="iid2d",n=2*n2) +
  f(U21, l.TIME,  copy="U11") +
  f(U22, copy = "U11", fixed = FALSE) 

joint.inla<- inla(formula, family = c("gaussian","weibullsurv"),
                  data = joint.data, control.compute=list(dic=TRUE,waic = TRUE,cpo = TRUE),control.family = list(list(),list()))
summary(joint.inla)


LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1

#################################################
formula = Y~ mu-1 +  l.Age +l.Gender+l.POI+ l.SS+s.Age +s.Gender+s.POI+ f(U11, model="iid2d",n=2*n2) +
  f(U21, l.TIME,  copy="U11") +
  f(U12, copy = "U11", fixed = FALSE) +
  f(U22, copy = "U11", fixed = FALSE) + f(ZONE_CODE, model = "besag",      
                                          graph = Lattice.adj)  

joint.inla<- inla(formula, family = c("gaussian","weibullsurv"),
                  data = joint.data, control.compute=list(dic=TRUE,waic = TRUE,cpo = TRUE),control.family = list(list(),list()))
summary(joint.inla)


LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1

#################################################
formula = Y~ mu-1 +  l.Age +l.Gender+l.POI+ l.SS+s.Age +s.Gender+s.POI+ f(U11, model="iid2d",n=2*n2) +
  f(U21, l.TIME,  copy="U11") +
  f(U12, copy = "U11", fixed = FALSE) +
  f(U22, copy = "U11", fixed = FALSE) + f(ZONE_CODE, model = "besag",      
                                          graph = Lattice.adj)  

joint.inla<- inla(formula, family = c("gaussian","weibullsurv"),
                  data = joint.data, control.compute=list(dic=TRUE,waic = TRUE,cpo = TRUE),control.family = list(list(),list()))
summary(joint.inla)


LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1
#################################################
formula = Y~ mu-1 +  l.Age +l.Gender+l.POI+ l.SS+s.Age +s.Gender+s.POI+ f(U11, model="iid2d",n=2*n2) +
  f(U21, l.TIME,  copy="U11") +
  f(U12, copy = "U11", fixed = FALSE) +
  f(U22, copy = "U11", fixed = FALSE) + f(ZONE_CODE, model = "besag",      
                                          graph = Lattice.adj)  

joint.inla<- inla(formula, family = c("gaussian","weibullsurv"),
                  data = joint.data, control.compute=list(dic=TRUE,waic = TRUE,cpo = TRUE),control.family = list(list(),list()))
summary(joint.inla)


LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1

#################################################
formula = Y~ mu-1 +  l.Age +l.Gender+l.POI+ l.SS+s.Age +s.Gender+s.POI+ f(U11, model="iid2d",n=2*n2) +
  f(U21, l.TIME,  copy="U11") +
  f(U12, copy = "U11", fixed = FALSE) +
  f(U22, copy = "U11", fixed = FALSE) + f(ZONE_CODE, model = "besag",      
                                          graph = Lattice.adj)  

joint.inla<- inla(formula, family = c("gaussian","weibullsurv"),
                  data = joint.data, control.compute=list(dic=TRUE,waic = TRUE,cpo = TRUE),control.family = list(list(),list()))
summary(joint.inla)


LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1




###################### Map 
formula = Y~ mu-1 +  l.Age +l.Gender+l.POI+ l.SS+s.Age +s.Gender+s.POI+ f(U11, model="iid2d",n=2*n2) +
  f(U21, l.TIME,  copy="U11") +
  f(U12, copy = "U11", fixed = FALSE) +
  f(U22, copy = "U11", fixed = FALSE)+ f(ZONE_CODE, model = "besag",       # spatial effect: ZONE_CODE is a numeric identifier for each area in the lattice  (does not work with factors)
                                         graph = Lattice.adj)  

#joint.inla <- inla(formula, family = c("gaussian","exponential.surv"),
#              data = joint.data, control.compute=list(dic=TRUE,waic = TRUE,cpo = TRUE),control.family = list(list(),list()))
#summary(joint.inla)
#LPML1 <- sum(log(joint.inla$cpo$cpo))
##########################################
joint.inla<- inla(formula, family = c("gaussian","weibullsurv"),
                  data = joint.data, control.compute=list(dic=TRUE,waic = TRUE,cpo = TRUE),control.family = list(list(),list()))
#summary(joint.inla)
LPML1 <- sum(log(joint.inla$cpo$cpo))
LPML1


####################### Map #########################
library(sf)
library(ggplot2)
library(viridis)

nc <- st_read("Map/Brasil.shp",
              quiet = TRUE
)

plot(nc )
SID74=joint.inla$summary.random$ZONE_CODE[,2]
ggplot(data = nc, aes(fill = SID74)) + geom_sf() +
  scale_fill_viridis() + theme_bw()

