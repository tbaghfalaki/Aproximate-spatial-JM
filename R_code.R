rm(list=ls()) 
library(mvtnorm)
library(INLA);library(extraDistr);library(splines);
library(SpatialEpi)
library(rgdal)
library(cplm)
library(survival)
library(ggplot2)
library(ggfortify)
###########
joint=list()
longdat <-list()# joint$longitudinal
survdat <-list()# joint$survival
set.seed(1)
Data=read.table("/Users/taban/Desktop/Github/Data_new.txt",header=TRUE)
#Data=read.table("C:/Users/taban/Desktop/Github/Data_new.txt",header=TRUE)
names(Data)
attach(Data)
head(Data)

library(CorrMixed)
#Spaghetti.Plot(Dataset = Data, Outcome = cd4, Id=patient, Time = timecd4,
     #        Add.Mean = FALSE, Add.Median = FALSE,Col="gray4")

n2=length(unique(patient))
Nobs=length(patient)
########################
Fox_Lattice <- readOGR("/Users/taban/Desktop/Github/Map/Brasil.shp")
#Fox_Lattice <- readOGR("C:/Users/taban/Desktop/Github/Map/Brasil.shp")
#plot(Fox_Lattice)

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

set.seed(1)
joint.inla<- inla(formula, family = c("gaussian","weibullsurv"),
              data = joint.data, control.compute=list(dic=TRUE,waic = TRUE),control.family = list(list(),list()))


summary(joint.inla)



##########3