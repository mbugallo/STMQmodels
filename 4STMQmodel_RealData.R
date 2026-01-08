
######################################################
###   Fitting MQ, TMQ and STMQ models. Real data   ###
######################################################

# Remove all objects from the current workspace to avoid conflicts
rm(list=ls())

# Load required functions for M-quantile regression
source("QRLM.R")
source("QRLMweights.R")

# Load required libraries
suppressWarnings(suppressMessages({
  library(MASS)
  library(dplyr)
  library(forecast)
  }))

# Data loading and preprocessing: harmonisation of identifiers and 
# ordering of sample and auxiliary datasets

weights.df <- read.csv('county_weights_matrix_haversine.csv', 
                       sep=',', row.names=1)
colnames(weights.df)  <- gsub("^X", "", colnames(weights.df))

data <- read.csv("data.csv")
data$GEOID <- sprintf("%05d", data$GeoFIPS) 
data$State.Code <- sprintf("%02d", data$State.Code) 
data$id <- 1:dim(data)[1]

data.aux <- read.csv("data_aux.csv")
data.aux$GEOID <- sprintf("%05d", data.aux$GeoFIPS) 
data.aux$State.Code <- sprintf("%02d", data.aux$State.Code) 
data.aux$id <- 1:dim(data.aux)[1]

# Read mapping for counties missing from the analysis sample
missing.counties.dist <- read.csv("missing_counties_dist.csv", 
                                  row.names=NULL)[, 2:5]
missing.counties.dist[, 1:2] <- apply(missing.counties.dist[, 1:2], 2, 
                                      function(x) sprintf("%05d", x))
missing.counties.dist[, 3:4] <- apply(missing.counties.dist[, 3:4], 2, 
                                      function(x) sprintf("%02d", x))

data.aux <- merge(data.aux, missing.counties.dist, by.x='GEOID', 
                  by.y='missing_positions.name', all.x=T, sort = FALSE)
data.aux$positions.name.min[is.na(data.aux$positions.name.min)] =
  data.aux$GEOID[is.na(data.aux$positions.name.min)]

data.aux$positions.name.min.ST[is.na(data.aux$positions.name.min.ST)] = 
  data.aux$missing_positions.name.ST[is.na(data.aux$missing_positions.name.ST)] = 
  data.aux$State.Code[is.na(data.aux$positions.name.min.ST)]

data.aux <- data.aux %>% arrange(id)
data <- data %>% arrange(id)
attach(data)

# Compute sample sizes and sampling fractions
kd <- table(State.Code, Year)
Kd <- table(data.aux$State.Code, data.aux$Year)
fd.100 <- ifelse(is.na(100 * kd / Kd), 0, 100 * kd / Kd)

y.s <- Median.AQI

# Standardized auxiliary data to fit the models
# Build design matrices for sample and prediction domains
x.s <- data.frame('Intercept' =  1, Poverty.Proportion.All.Ages, 
                  Poverty.Proportion.Age.0.17,  Median.Household.Income, 
                  Real.GDP, Current.Dollar.GDP, Total.Pop.County, 
                  Dens.Pop.County, Total.Pop.State, Dens.Pop.State, 
                  ID.Coast, ID.Metrop)

x.r <- data.frame('Intercept' =  1, data.aux$Poverty.Proportion.All.Ages, 
                  data.aux$Poverty.Proportion.Age.0.17, 
                  data.aux$Median.Household.Income, 
                  data.aux$Real.GDP, data.aux$Current.Dollar.GDP, 
                  data.aux$Total.Pop.County, data.aux$Dens.Pop.County, 
                  data.aux$Total.Pop.State, data.aux$Dens.Pop.State, 
                  data.aux$ID.Coast, data.aux$ID.Metrop)

media.r <- colMeans(x.r[, 4:10], na.rm = TRUE)
sd.r <- apply(x.r[, 4:10], 2, sd, na.rm = TRUE)
for (i in 1:7){
  x.s[, (3+i)] <- (x.s[, (3+i)] - media.r[i])/sd.r[i]
  x.r[, (3+i)] <- (x.r[, (3+i)] - media.r[i])/sd.r[i] }

correlations <- lapply(apply(x.s[, 2:10], 2, function(x) cor.test(x, y.s)), 
                       function(x) c(estimate = round(x$estimate,3), 
                                     conf.int = round(x$conf.int, 3)))

names(x.r) <- names(x.s)

# Delete variables: Median.Household.Income, Dens.Pop.County, Dens.Pop.State
x.s <- as.matrix(x.s[, !(names(x.s) %in% c('Median.Household.Income', 
                                           'Dens.Pop.County', 'Dens.Pop.State', 
                                           'Total.Pop.State'))])
x.r <- as.matrix(x.r[, !(names(x.r) %in% c('Median.Household.Income', 
                                           'Dens.Pop.County', 'Dens.Pop.State', 
                                           'Total.Pop.State'))])

regioncode.s <- Postal.Code
timecode.s <- Year
GEOID.s <- GEOID

regioncode.r <- data.aux$Postal.Code
timecode.r <- data.aux$Year
GEOID.r <- data.aux$GEOID
GEOID.r.in.s <- data.aux$positions.name.min

states.sr <- sort(unique(regioncode.s))
GEOID.sr <- sort(unique(GEOID.s))
years <- sort(unique(timecode.s))

n <- dim(data)[1]
m <- length(unique(regioncode.s))
T <- length(unique(timecode.s))
nit <- table(regioncode.s, timecode.s)
nu.t <- table(timecode.s) + c(0,table(timecode.s)[-T])
n.t <- colSums(nit)


#############################
###     Area-level MQ     ###
#############################  

# Area-level MQ estimation using a two-step quantile regression approach

mod50 <- QRLM(x=x.s, y=y.s, q=0.50, maxit=35, k = 1.345)
coef50 <- mod50$coef
res50 <- mod50$residuals
s50 <- fun.MAD(res50)
p <- length(coef50)

var.beta50 <- n^2*s50^2/(n-p)*sum(hub.psi(res50/s50, k=1.345)^2)/
  (sum(der.hub.psi(res50/s50, k=1.345))^2)*solve(t(x.s)%*%x.s)
sd.beta.diag50 <- sqrt(diag(var.beta50))

(IC50 <-round(cbind(coef50 + qnorm(0.025)*sd.beta.diag50/sqrt(n), 
                   coef50 - qnorm(0.025)*sd.beta.diag50/sqrt(n)), 3))


# 1º Quantile regression
tau <- sort(unique(c(seq(0.006, 0.99, 0.045), 1- seq(0.006,0.99,0.045), 0.50)))
mod <- QRLM(x=x.s, y=y.s, q=tau, maxit=35, k = 1.345)

## Linear Interpolation
qo <- matrix(c(gridfitinter(y.s,mod$fitted.values,mod$q.values)),nrow=n,ncol=1)
qmat <- matrix(c(qo,regioncode.s),nrow=n,ncol=2)
mqo <- aggregate(as.numeric(qmat[,1]),by=list(qmat[,2]),mean)[, 2]

# 2º Quantile regression
mod.SAE <- QRLM(x=x.s, y=y.s, q=mqo, maxit=35, k = 1.345)
beta.SAE <- mod.SAE$coef


#############################
###    Predicting MQ      ###
#############################

# Prediction for sampled and non-sampled units under the MQ-SAE framework
pred.st<-pred.rt<-id.list.s<-id.list.r<-list()
for (j in 1:m){  pred.st[[j]]<-id.list.s[[j]]<-
  pred.rt[[j]]<-id.list.r[[j]]<-list()
	for (t in 1:T){ 
		pred.st[[j]][[t]] <- x.s[regioncode.s==states.sr[j] & 
		                     timecode.s==years[t],]%*%beta.SAE[,j]  
		pred.rt[[j]][[t]] <- x.r[regioncode.r==states.sr[j] & 
		                     timecode.r==years[t],]%*%beta.SAE[,j]
		id.list.s[[j]][[t]]<- data$id[regioncode.s==states.sr[j] &  
		                                timecode.s==years[t]]
		id.list.r[[j]][[t]]<- data.aux$id[regioncode.r==states.sr[j] & 
		                                    timecode.r==years[t]]
	}	
} 

pred.st.id <- data.frame(unlist(id.list.s), unlist(pred.st))
pred.st.id <- pred.st.id %>% arrange(unlist.id.list.s.)
names(pred.st.id) <- c('id', 'pred.st')
pred.st.id <- merge(data, pred.st.id, by='id')

pred.rt.id <- data.frame(unlist(id.list.r), unlist(pred.rt))
pred.rt.id <- pred.rt.id %>% arrange(unlist.id.list.r.)
names(pred.rt.id) <- c('id', 'pred.rt')
pred.rt.id <- merge(data.aux, pred.rt.id, by='id')

pred.st.id <- merge(pred.st.id, pred.rt.id[, c('GEOID', 'Year', 'pred.rt')], 
                    by=c('GEOID', 'Year'), sort=F)

res.st <- y.s- pred.st.id $pred.st
res.st <- res.st/sd(res.st)


##################################
###    Temporal MQ models      ###
##################################

# Temporal MQ (TMQ): borrowing strength across neighbouring time periods
wt.results <- compute.wt(x.s, y.s, regioncode.s, timecode.s, mod.SAE, x.r)

ww <- wt.results$ww
P <- wt.results$P

years.incl.t <- list()
for(t in 1:T){
  if ((t == 1)|(P==0)){ incl.t = t } else { incl.t = t-P+1; 
  if (incl.t < 1){ incl.t=1 }}
  years.incl.t[[t]] <- years[incl.t:t] 
}

# 1º Quantile regression

mod.TMQ <- lapply(1:T, function(t) {
    cond.t <- which(timecode.s %in% years.incl.t[[t]])
    return(QRLMweights(x = x.s[cond.t , ], y = y.s[cond.t], q = tau, 
                       maxit = 35, k = 1.345,
                       w1 = diag(ww[[t]])[cond.t, cond.t] ))
})

fitted.values <-list()
for (t in 1:T){ 
  cond.tt <- which(timecode.s == years[t])-which(timecode.s %in% 
                                                   years.incl.t[[t]])[1]+1 
  fitted.values[[t]] <- mod.TMQ[[t]]$fitted.values[cond.tt, ] 
}

fitted.values <- Reduce('rbind', fitted.values)

## Linear Interpolation

qo.TMQ <- matrix(c(gridfitinter(y.s, fitted.values, tau)), nrow=n, ncol=1)
qmat.TMQ<-matrix(c(qo.TMQ, regioncode.s), nrow=n, ncol=2)
mqo.TMQ.L <- aggregate(as.numeric(qmat.TMQ[,1]), by=list(qmat.TMQ[,2]),mean)
mqo.TMQ <- mqo.TMQ.L[, 2]

# 2º Quantile regression

modSAE.TMQ <- lapply(1:T, function(t) {
  cond.t <- which(timecode.s %in% years.incl.t[[t]])
  return(QRLMweights(x = x.s[cond.t , ], y = y.s[cond.t], q = mqo.TMQ, 
                     maxit = 35, k = 1.345,
                     w1 = diag(ww[[t]])[cond.t, cond.t] ))
})


######################################
###    Predicting Temporal MQ      ###
######################################

# Predicted values
pred.TMQ.st<-pred.TMQ.rt<-id.list.TMQ.s<-id.list.TMQ.r<-list()
for (j in 1:m){  pred.TMQ.st[[j]]<-id.list.TMQ.s[[j]]<-
  pred.TMQ.rt[[j]]<-id.list.TMQ.r[[j]]<-list()
for (t in 1:T){ 
  pred.TMQ.st[[j]][[t]] <- x.s[regioncode.s==states.sr[j] & 
                             timecode.s==years[t],]%*%modSAE.TMQ[[t]]$coefficients[,j]  
  pred.TMQ.rt[[j]][[t]] <- x.r[regioncode.r==states.sr[j] & 
                             timecode.r==years[t],]%*%modSAE.TMQ[[t]]$coefficients[,j]
  id.list.TMQ.s[[j]][[t]]<- data$id[regioncode.s==states.sr[j] &  
                                  timecode.s==years[t]]
  id.list.TMQ.r[[j]][[t]]<- data.aux$id[regioncode.r==states.sr[j] & 
                                      timecode.r==years[t]]
}	
} 

pred.TMQ.st.id <- data.frame(unlist(id.list.TMQ.s), unlist(pred.TMQ.st))
pred.TMQ.st.id <- pred.TMQ.st.id %>% arrange(unlist.id.list.TMQ.s.)
names(pred.TMQ.st.id) <- c('id', 'pred.st')
pred.TMQ.st.id <- merge(data, pred.TMQ.st.id, by='id')

pred.TMQ.rt.id <- data.frame(unlist(id.list.TMQ.r), unlist(pred.TMQ.rt))
pred.TMQ.rt.id <- pred.TMQ.rt.id %>% arrange(unlist.id.list.TMQ.r.)
names(pred.TMQ.rt.id) <- c('id', 'pred.rt')
pred.TMQ.rt.id <- merge(data.aux, pred.TMQ.rt.id, by='id')

pred.TMQ.st.id <- merge(pred.TMQ.st.id, pred.TMQ.rt.id[, c('GEOID', 'Year', 'pred.rt')], 
                    by=c('GEOID', 'Year'), sort=F)

res.TMQ.st <- y.s- pred.TMQ.st.id $pred.st
res.TMQ.st <- res.TMQ.st/sd(res.TMQ.st)


######################################
###    Temporal MQ-GWR models      ###
######################################

# Spatio-temporal MQ-GWR (STMQ): incorporation of spatial weights 
# based on geographic distance
weights.df.counties <- lapply(lapply(1:nrow(weights.df), function(j) {
  county_data <- unlist(lapply(years, function(i) {
    weights.df[j, colnames(weights.df) %in% data$GEOID[data$Year == i]]
  }))
  county_data[!is.na(county_data)]
}), as.numeric)


# 1º Quantile regression

mod.STMQ <- lapply(1:length(weights.df.counties), function(j) {
  lapply(1:T, function(t) {
    cond.t <- which(timecode.s %in% years.incl.t[[t]])
    return(QRLMweights(x = x.s[cond.t , ], y = y.s[cond.t], q = tau, 
              maxit = 35, k = 1.345,
              w1 = diag(weights.df.counties[[j]] * ww[[t]])[cond.t, cond.t] ))
  })
})

fitted.values <- matrix(nrow=n, ncol=length(tau))
for (j in 1:n){ 
  t <-which(timecode.s[j] == years)
  cond.t <- which(timecode.s %in% years.incl.t[[t]])
  id.tj <- which(id[cond.t]==j)
  fitted.values[j, ] <- mod.STMQ[[which(GEOID.s[j] == 
                          GEOID.sr)]][[t]]$fitted.values[id.tj , ] 
  }

## Linear Interpolation

qo.STMQ<-matrix(c(gridfitinter(y.s, fitted.values, tau)),nrow=n,ncol=1)
qmat.STMQ<-matrix(c(qo.STMQ, regioncode.s), nrow=n, ncol=2)
mqo.STMQ <- aggregate(as.numeric(qmat.STMQ[,1]),by=list(qmat.STMQ[,2]),mean)
names(mqo.STMQ) <- c('Postal.Code', 'mqo.STMQ')
data <- merge(data, mqo.STMQ, by='Postal.Code') %>% arrange(id)
data.aux <- merge(data.aux, mqo.STMQ, by='Postal.Code') %>% arrange(id)

mqo.GEOID.STMQ <- unique(data[, c('GEOID', 'mqo.STMQ')]) %>% arrange(GEOID)

rm(mod.STMQ)

# 2º Quantile regression

modSAE.STMQ <- lapply(1:length(weights.df.counties), function(j) {
  lapply(1:T, function(t) {
    cond.t <- which(timecode.s %in% years.incl.t[[t]])
    return(QRLMweights(x = x.s[cond.t , ], y = y.s[cond.t], q = mqo.GEOID.STMQ[j,2], 
                   maxit = 35, k = 1.345, 
                   w1 = diag(weights.df.counties[[j]] * ww[[t]])[cond.t, cond.t] ))
  })
})

beta.SAE.STMQ <- data.frame(rep(GEOID.sr, each=T), rep(years, times=length(GEOID.sr)), 
                             Reduce(rbind, lapply(modSAE.STMQ, function(x){ t(sapply(x, 
                             function(y){ return(y$'coefficients') })) })) )
colnames(beta.SAE.STMQ) <- c('GEOID', 'Year', paste('beta', 1:p, sep=''))


mod50.STMQ <- lapply(1:length(weights.df.counties), function(j) {
  lapply(1:T, function(t) {
    cond.t <- which(timecode.s %in% years.incl.t[[t]])
    return(QRLMweights(x = x.s[cond.t , ], y = y.s[cond.t], q = 0.50, 
                       maxit = 35, k = 1.345, 
                       w1 = diag(weights.df.counties[[j]] * ww[[t]])[cond.t, cond.t] ))
  })
})


##########################################
###     Predicting Temporal MQ-GWR     ###
##########################################

STMQ.pred.results <- STMQ.pred(y.s, x.s, x.r, timecode.r, GEOID.sr, 
                      beta.SAE.STMQ, years, years.incl.t)

pred.st.N.GWR <- STMQ.pred.results[[1]]
# pred.st.GWR <- STMQ.pred.results[[2]]
pred.rt.N.GWR <- STMQ.pred.results[[3]]
# pred.rt.GWR <- STMQ.pred.results[[4]]
id.list.s.GWR <- STMQ.pred.results[[5]] 
id.list.r.GWR <- STMQ.pred.results[[6]] 

pred.GWR <- data.frame(unlist(id.list.s.GWR), unlist(pred.st.N.GWR))
names(pred.GWR) <- c('id', 'pred.st')
pred.GWR <- merge(data, pred.GWR, by='id') %>% arrange(id)

predr.GWR <- data.frame(unlist(id.list.r.GWR), unlist(pred.rt.N.GWR))
names(predr.GWR) <- c('id', 'pred.rt')
predr.GWR <- merge(data.aux, predr.GWR, by='id') %>% arrange(id)

pred.GWR <- merge(pred.GWR, predr.GWR[, c('GEOID', 'Year', 'pred.rt')], 
                        by=c('GEOID', 'Year'), sort=F)

res.st.GWR <- y.s - pred.GWR$pred.st
res.st.GWR <- res.st.GWR/sd(res.st.GWR)


##########################################
###   MSE estimation Temporal MQ-GWR  ###
##########################################

MSE <- STMQ.mse(x.s, y.s, x.r, modSAE.STMQ, mod50.STMQ, 
                   regioncode.s, timecode.s, regioncode.r, timecode.r)

predN.mse <- MSE$pred.rt.N.GWR
predBC.mse <- MSE$pred.rt.BC.GWR
MSE.id <- MSE$id.list.r.GWR
c.phi <- MSE$c.phi

bias  <- MSE$bias
rmse1 <- MSE$rmse1
rmse2 <- MSE$rmse2
rmse3 <- MSE$rmse3

biasBC  <- MSE$biasBC
rmse1BC <- MSE$rmse1BC
rmse2BC <- MSE$rmse2BC
rmse3BC <- MSE$rmse3BC

MSE <- merge(data.aux, data.frame(MSE.id, predN.mse, predBC.mse, 
                       c.phi, bias, rmse1, rmse2, rmse3,  rmse1BC, 
                       rmse2BC, rmse3BC, biasBC), by.x='id', 
                       by.y='MSE.id') %>% arrange(id)

data$dom <- paste(data[['GEOID']], data[['Year']], sep = ":")
MSE$dom <- paste(MSE[['GEOID']], MSE[['Year']], sep = ":")

MSE$id.in.s <- ifelse(MSE$dom %in% data$dom, 1, 0)
MSE$rrmse1 <- MSE$rmse1/MSE$predN.mse
MSE$rrmse2 <- MSE$rmse2/MSE$predN.mse
MSE$rrmse3 <- MSE$rmse3/MSE$predN.mse

MSE$rrmse1BC <- MSE$rmse1BC/MSE$predBC.mse
MSE$rrmse2BC <- MSE$rmse2BC/MSE$predBC.mse
MSE$rrmse3BC <- MSE$rmse3BC/MSE$predBC.mse

MSEm <- MSE %>% group_by(Postal.Code, Year) %>%
  summarise(c.phi = mean(c.phi, na.rm = TRUE), .groups = "drop")

c.phi2 <- reshape2::acast(MSE, GEOID ~ Year, value.var = "c.phi")
c.phi3 <- reshape2::acast(MSEm, Postal.Code ~ Year, value.var = "c.phi")

#############################
###   Outlier detection   ###
#############################
# Detection of spatial and temporal outliers based on MSE diagnostics

friedman <- data.frame('respuesta'=MSE$c.phi, 'area'=MSE$Postal.Code, 
                       'time'=MSE$Year)

# Area-level outliers         
Tukey <- agricolae::HSD.test(aov(friedman$respuesta~friedman$time+friedman$area),
                             "friedman$area", group=TRUE,alpha=0.05); Tukey

# Time-level outliers 
Tukey <- agricolae::HSD.test(aov(friedman$respuesta ~ friedman$time+friedman$area),
                             "friedman$time",  group=TRUE,alpha=0.05); Tukey

mqo.TMQ.L[mqo.TMQ.L$Group.1 %in% names(which(apply(t(c.phi3), 
                                             2, max, na.rm=T)>3)), ]

# save.image('STMQpred.Rdata')
