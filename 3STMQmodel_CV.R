
#######################################################
###   Cross validation for STMQ models. Real data   ###
#######################################################

# INPUT DATASETS:
#  - data.csv : analysis dataset with target variable (AQI) and covariates
#  - data_aux.csv : auxiliary-only dataset for out-of-sample prediction 
#  - missing_counties_dist.csv : mapping of missing counties to nearest available counties
#  - county_distance_matrix_haversine.csv : county-to-county distance matrix

# EXTERNAL CODE:
#  - QRLM.R          : M-quantile regression fitting functions
#  - QRLMweights.R   : weighted M-quantile regression functions

# OUTPUT OBJECTS:
#  - Cross-validation results for the spatial decay parameter b

# NOTES:
#  - Auxiliary covariates are standardized using statistics from the auxiliary sample
#  - Spatial and temporal borrowing of strength are combined in STMQ models
#  - Leave-one-county cross validation is performed for each candidate value of b

# Clear the workspace
rm(list=ls())

# Load required functions for M-quantile regression
source("QRLM.R")
source("QRLMweights.R")

# Load required libraries
suppressWarnings(
  suppressMessages({
    library(MASS)
    library(dplyr)
    library(forecast)
    library(ggplot2)
}))

# Read target and auxiliary datasets and format identifiers
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

# Match missing counties to nearest available counties
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

# Compute sample sizes by state and year
kd <- table(State.Code, Year)
Kd <- table(data.aux$State.Code, data.aux$Year)
fd.100 <- ifelse(is.na(100 * kd / Kd), 0, 100 * kd / Kd)

# Define the response variable (median AQI)
y.s <- Median.AQI

# Standardized auxiliary data to fit the models
# Build design matrices for sample and prediction domains
x.s <- data.frame('Intercept' =  1, Poverty.Proportion.All.Ages, 
                  Poverty.Proportion.Age.0.17,  Median.Household.Income, 
                  Real.GDP, Current.Dollar.GDP,Total.Pop.County, 
                  Dens.Pop.County, Total.Pop.State, Dens.Pop.State, 
                  ID.Coast, ID.Metrop)

x.r <- data.frame('Intercept' =  1, data.aux$Poverty.Proportion.All.Ages, 
                  data.aux$Poverty.Proportion.Age.0.17, 
                  data.aux$Median.Household.Income, 
                  data.aux$Real.GDP, data.aux$Current.Dollar.GDP, 
                  data.aux$Total.Pop.County, data.aux$Dens.Pop.County, 
                  data.aux$Total.Pop.State, data.aux$Dens.Pop.State, 
                  data.aux$ID.Coast, data.aux$ID.Metrop)

# Standardize continuous covariates using auxiliary-domain statistics
media.r <- colMeans(x.r[, 4:10], na.rm = TRUE)
sd.r <- apply(x.r[, 4:10], 2, sd, na.rm = TRUE)
for (i in 1:7){
  x.s[, (3+i)] <- (x.s[, (3+i)] - media.r[i])/sd.r[i]
  x.r[, (3+i)] <- (x.r[, (3+i)] - media.r[i])/sd.r[i] }

# Compute correlations between covariates and the target variables
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

# Define region, time and spatial identifiers
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
n.t <- colSums(table(regioncode.s, timecode.s))
ni. <- rowSums(table(regioncode.s, timecode.s))
n.ts <- c(0,cumsum(n.t))


#############################
###     Area-level MQ     ###
#############################  

# 1º Quantile regression
tau <- sort(unique(c(seq(0.006, 0.99, 0.045), 1- seq(0.006,0.99,0.045), 0.50)))
mod <- QRLM(x=x.s, y=y.s, q=tau, maxit=35, k = 1.345)

## Linear Interpolation
qo<-matrix(c(gridfitinter(y.s,mod$fitted.values,mod$q.values)),nrow=n,ncol=1)
qmat<-matrix(c(qo,regioncode.s),nrow=n,ncol=2)
mqo<-aggregate(as.numeric(qmat[,1]),by=list(qmat[,2]),mean)
names(mqo) <- c('Postal.Code', 'mqo')
data <- merge(data, mqo, by='Postal.Code') %>% arrange(id)
mqo <- mqo$mqo

# 2º Quantile regression
mod.SAE<-QRLM(x=x.s, y=y.s, q=mqo, maxit=35, k = 1.345)
beta.SAE <- mod.SAE$coef


#############################
###    Predicting MQ      ###
#############################

# Generate spatio-temporal predictions for sample and auxiliary units
pred.st<-pred.rt<-id.list.s<-id.list.r<-list()
for (j in 1:m){  pred.st[[j]] <- id.list.s[[j]]<-
  pred.rt[[j]] <- id.list.r[[j]] <- list()
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

# Compute standardized residuals
res.st <- y.s- pred.st.id $pred.rt
res.st <- res.st/sd(res.st)


###########################################
###    Temporal MQ-GWR models - CV      ###
###########################################

# Compute temporal weights to borrow strength across years
wt.results <- compute.wt(x.s, y.s, regioncode.s, timecode.s, mod.SAE, x.r)

ww <- wt.results$ww
P <- wt.results$P

years.incl.t <- list()
for(t in 1:T){
  if ((t == 1)|(P==0)){ incl.t = t } else { incl.t = t-P+1; if (incl.t < 1){ incl.t=1 }}
  years.incl.t[[t]] <- years[incl.t:t] 
 }

# Borrowing strength from spatial locations
# Read spatial distance matrix between counties
distance.df <- read.csv("county_distance_matrix_haversine.csv", sep=',', row.names=1)
county_names <- colnames(distance.df)  <- gsub("^X", "", colnames(distance.df))

# Cross validation
mod.STMQ50 <- CV.dif2 <- list()
iter <- 0
# Define grid of spatial decay parameters for cross validation
grid.b <- seq(0, 50, by=.5)

# Perform leave-one-county-out cross validation
for(b in grid.b){
  weights <- function(x, b){ exp( -b*x) }
  weights.df <- apply(distance.df, 1, weights,  b=b)
  rownames(weights.df) <- colnames(weights.df) <- county_names
  
  weights.df.counties <- lapply(lapply(1:nrow(weights.df), function(j) {
    county_data <- unlist(lapply(years, function(i) {
      weights.df[j, colnames(weights.df) %in% data$GEOID[data$Year == i]]
    }))
    county_data[!is.na(county_data)]
  }), as.numeric)
  
  iter <- iter +1
  mod.STMQ50[[iter]] <- lapply(1:length(weights.df.counties), function(j) {
    lapply(1:T, function(t) {
      cond.t <- which(timecode.s %in% years.incl.t[[t]])
      delete.CV <- which(data[cond.t, ]$GEOID%in% colnames(weights.df)[j])
      
      if(length(delete.CV) == 0){
        # Fit weighted M-quantile models excluding the target county
        model <- QRLMweights(x = (x.s[cond.t , ]), y = (y.s[cond.t]), 
                             q = data$mqo[GEOID==GEOID.sr[j]][1], 
                             maxit = 35, k = 1.345,
                             w1 = (diag(weights.df.counties[[j]]*ww[[t]])[cond.t,
                                                                          cond.t]))
      } else {
        model <- QRLMweights(x = (x.s[cond.t , ])[-delete.CV, ], 
                             y = (y.s[cond.t])[-delete.CV], 
                             q = data$mqo[GEOID==GEOID.sr[j]][1], 
                             maxit = 35, k = 1.345,
                             w1 = (diag(weights.df.counties[[j]]*ww[[t]])[cond.t, 
                                              cond.t])[-delete.CV, -delete.CV] )
      }
      return(model=model)
    })
  })
  
  # Compute cross-validation loss function
  fitted.values.CV <- sapply(1:n, function(j) {
    t <- which(timecode.s[j] == years)
    cond.t <- which(timecode.s %in% years.incl.t[[t]])
    fitted.values.aux <- (x.s[cond.t, ])[data[cond.t, ]$GEOID %in% GEOID[j], ] %*% 
      mod.STMQ50[[iter]][[which(GEOID.s[j] == GEOID.sr)]][[t]]$coefficients
    tail(fitted.values.aux, 1)  
  })
  
  CV.dif2[[iter]] <- sum((y.s-fitted.values.CV)^2)
  # Identify the optimal spatial decay parameter
  print(sprintf("For b=%d, the value of CV.dif is %f", b, CV.dif2[[iter]]))
}


## The results are: 

values <- c(538406.9, 538405.9, 538404.8, 538404.2, 538403.6, 538402.8, 538402.0, 
            538401.1, 538400.1, 538399.2, 538398.3, 538397.4, 538396.4, 538395.7, 
            538395.0, 538394.2, 538393.5, 538392.7, 538391.8, 538391.2, 538390.5, 
            538389.8, 538389.0, 538388.3, 538387.6, 538387.1, 538386.5, 538386.1, 
            538385.7, 538385.6, 538385.5, 538385.2, 538384.9, 538384.4, 538383.9, 
            538383.5, 538383.1, 538382.6, 538382.2, 538381.7, 538381.2, 538380.8, 
            538380.3, 538379.9, 538379.6, 538379.3, 538379.0, 538378.7, 538378.4, 
            538378.2, 538378.0, 538377.8, 538377.5, 538377.2, 538377.0, 538376.7, 
            538376.4, 538376.2, 538376.0, 538375.8, 538375.7, 538375.6, 538375.4, 
            538375.1, 538374.7, 538374.2, 538373.8, 538373.6, 538373.4, 538373.3, 
            538373.2, 538373.6, 538374.0, 538373.9, 538373.9, 538374.1, 538374.2, 
            538374.6, 538374.9, 538374.9, 538375.0, 538375.2, 538375.4, 538375.9, 
            538376.5, 538376.9, 538377.4, 538377.8, 538378.2, 538378.6, 538378.9, 
            538379.2, 538379.6, 538379.9, 538380.2, 538380.6, 538381.0, 538381.4, 
            538381.9)

# Plot cross-validation loss as a function of b
quartz(width = 6, height = 6) 
df <- data.frame(x = seq(1, 50, by=0.5), y = sqrt(values))
ggplot(df, aes(x, y)) +
  geom_line(color='blue', size=1) +
  geom_point(aes(x = x[which.min(y)], y = min(y)), 
             color='red', size=3) +
  geom_hline(yintercept=min(df$y), linetype="dashed", color="red") +
  labs(title=" ", x="b", y='Square root of the loss function') +
  theme_minimal() +  
  theme(legend.position = c(0.9, 0.25), 
        axis.title.x = element_text(size = 14),    
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14))

# Best result: b=1/36; 1/b=36



