
### QUANTILE REGRESSION FITTING ALGORITHM WITH WEIGHTS - IRLS

QRLMweights <- function (x, y, case.weights = rep(1, nrow(x)), var.weights = rep(1, nrow(x)), ...,
                         w = rep(1, nrow(x)), init = "ls", psi = psi.huber, scale.est = c("MAD", 
                         "Huber", "proposal 2"), k2 = 1.345, method = c("M", "MM"), maxit = 20, 
                         acc = 1e-04, test.vec = "resid", q = 0.5, w1)
{
  irls.delta <- function(old, new) sqrt(sum((old - new)^2)/max(1e-20, sum(old^2)))
  irls.rrxwr <- function(x, w, r) {
    w <- sqrt(w)
    max(abs((matrix(r*w,1,length(r)) %*% x)/sqrt(matrix(w,1,length(r)) %*% (x^2))))/sqrt(sum(w*r^2))
  }
  method <- match.arg(method)
  nmx <- deparse(substitute(x))
  if (is.null(dim(x))) {
    x <- as.matrix(x)
    colnames(x) <- nmx
  }
  else x <- as.matrix(x)
  if (is.null(colnames(x))) 
    colnames(x) <- paste("X", seq(ncol(x)), sep = "")
  if (qr(x)$rank < ncol(x)) 
    stop("x is singular: singular fits are not implemented in rlm")
  if (!(any(test.vec == c("resid", "coef", "w", "NULL")) || is.null(test.vec))) 
    stop("invalid testvec")
  if (length(var.weights) != nrow(x)) 
    stop("Length of var.weights must equal number of observations")
  if (any(var.weights < 0)) 
    stop("Negative var.weights value")
  if (length(case.weights) != nrow(x)) 
    stop("Length of case.weights must equal number of observations")
  w <- (w * case.weights)/var.weights
  if (method == "M") {
    scale.est <- match.arg(scale.est)
    if (!is.function(psi)) 
      psi <- get(psi, mode = "function")
    arguments <- list(...)
    if (length(arguments)) {
      pm <- pmatch(names(arguments), names(formals(psi)), nomatch = 0)
      if (any(pm == 0)) 
        warning(paste("some of ... do not match"))
      pm <- names(arguments)[pm > 0]
      formals(psi)[pm] <- unlist(arguments[pm])
    }
    if (is.character(init)) {
      if (init == "ls") 
        temp <- lm.wfit(x, y, w, method = "qr")
      else if (init == "lts") 
        temp <- lqs.default(x, y, intercept = FALSE, nsamp = 200)
      else stop("init method is unknown")
      coef <- temp$coef
      resid <- temp$resid
    }
    else {
      if (is.list(init)) 
        coef <- init$coef
      else coef <- init
      resid <- y - x %*% coef
    }
  }
  else if (method == "MM") {
    scale.est <- "MM"
    temp <- lqs.default(x, y, intercept = FALSE, method = "S", k0 = 1.548)
    coef <- temp$coef
    resid <- temp$resid
    psi <- psi.bisquare
    if (length(arguments <- list(...))) 
      if (match("c", names(arguments), nomatch = FALSE)) {
        c0 <- arguments$c
        if (c0 > 1.548) {
          psi$c <- c0
        }
        else warning("c must be at least 1.548 and has been ignored")
      }
    scale <- temp$scale
  }
  else stop("method is unknown")
  done <- FALSE
  conv <- NULL
  n1 <- nrow(x) - ncol(x)
  if (scale.est != "MM") 
    scale <- mad(resid/sqrt(var.weights), 0)
  theta <- 2 * pnorm(k2) - 1
  gamma <- theta + k2^2 * (1 - theta) - 2 * k2 * dnorm(k2)
  qest <- matrix(0, nrow = ncol(x), ncol = length(q))
  qwt <- matrix(0, nrow = nrow(x), ncol = length(q))
  qfit <- matrix(0, nrow = nrow(x), ncol = length(q))
  qres <- matrix(0, nrow = nrow(x), ncol = length(q))
  for(i in 1:length(q)) {
    for (iiter in 1:maxit) {
      if (!is.null(test.vec)) 
        testpv <- get(test.vec)
      if (scale.est != "MM") {
        if (scale.est == "MAD") 
          scale <- median(abs(resid/sqrt(var.weights)))/0.6745
        else scale <- sqrt(sum(pmin(resid^2/var.weights,(k2*scale)^2))/(n1*gamma))
        if (scale == 0) {
          done <- TRUE
          break
        }
      }
      w <- psi(resid/(scale * sqrt(var.weights))) * case.weights
      ww <- 2 * (1 - q[i]) * w
      ww[resid > 0] <- 2 * q[i] * w[resid > 0]
      w <- ww*diag(w1)
      temp <- lm.wfit(x, y, w, method = "qr")
      coef <- temp$coef
      resid <- temp$residuals
      if (!is.null(test.vec)) 
        convi <- irls.delta(testpv, get(test.vec))
      else convi <- irls.rrxwr(x, wmod, resid)
      conv <- c(conv, convi)
      done <- (convi <= acc)
      if (done) 
        break
    }
    if (!done) 
      warning(paste("rlm failed to converge in", maxit, "steps at q = ", q[i]))
    qest[, i] <- coef
    qwt[, i] <- w
    qfit[, i] <- temp$fitted.values
    qres[,i] <- resid
  }
  list(fitted.values = qfit, residuals = qres, q.values = q, q.weights = qwt, coefficients = qest)
}

# COMPUTING OF THE QUANTILE-ORDERS
"zerovalinter"<-function(y, x)
{
  if(min(y) > 0) {
    xmin <- x[y == min(y)]
    if(length(xmin) > 0)
      xmin <- xmin[length(xmin)]
    xzero <- xmin
  }
  else {
    if(max(y) < 0) {
      xmin <- x[y == max(y)]
      if(length(xmin) > 0)
        xmin <- xmin[1]
      xzero <- xmin
    }
    else {
      y1 <- min(y[y > 0])
      if(length(y1) > 0)
        y1 <- y1[length(y1)]
      y2 <- max(y[y < 0])
      if(length(y2) > 0)
        y2 <- y2[1]
      x1 <- x[y == y1]
      if(length(x1) > 0)
        x1 <- x1[length(x1)]
      x2 <- x[y == y2]
      if(length(x2) > 0)
        x2 <- x2[1]
      xzero <- (x2 * y1 - x1 * y2)/(y1 - y2)
      xmin <- x1
      if(abs(y2) < y1)
        xmin <- x2
    }
  }
  resu <-  xzero
  resu
}

# LINEAR INTERPOLATION FUNCTION
# It assumes that the "zerovalinter" function has been already loaded

"gridfitinter"<-function(y,expectile,Q)
  # Computing of the expectile-order of each observation of y by linear interpolation
{
  nq<-length(Q)
  diff <- y %*% t(as.matrix(rep(1, nq))) - expectile        
  vectordest <- apply(diff, 1, zerovalinter,Q)    
}

# BORROWING STRENGTH FROM TIME: TWMQ WEIGHTS

compute.wt <- function(x.s, y.s, regioncode.s, timecode.s, mod.SAE, x.r){
  m <- length(unique(regioncode.s))
  T <- length(unique(timecode.s))
  p <- dim(x.s)[2]
  n <- dim(x.s)[1]
  
  # Predicted values
  pred.st<-list()
  for (j in 1:m){ pred.st[[j]]<-list()
  for (t in 1:T){ pred.st[[j]][[t]] <- as.matrix(x.s[regioncode.s==states.sr[j] & 
                                                       timecode.s==years[t],])%*%beta.SAE[,j]}
  }
  
  # Time series analysis
  serie<- ts(aggregate(y.s-unlist(pred.st), by=list(regioncode.s,timecode.s), mean)[,3],frequency=m)
  serie<-serie+abs(min(serie))+1
  
  fit.Arima <- auto.arima(serie, max.q=0, max.Q=0, max.P=0, 
                          seasonal=TRUE, approximation=FALSE, stepwise=TRUE)
  P<-fit.Arima$arma[1]; P
  
  # Definition of time weights
  if (P>=1){
    coefARabs<-as.numeric(abs(fit.Arima$coef[P:1])/sum(abs(fit.Arima$coef[P:1])))
    w<-ww<-list()
    for (t in 1:T){
      if(t <= P){ w[[t]]<-c(coefARabs[1:t],rep(0, T-t))/sum(coefARabs[1:t]) } 
      else {	w[[t]]<-c(rep(0,t-P),coefARabs,rep(0, T-t)) }
      ww[[t]]<-list()
      for (tt in 1:T){ ww[[t]][[tt]]<-rep(w[[t]][tt], each=n.t[tt]) }
      ww[[t]]<-unlist(ww[[t]]) }
    
    order.s<-with(data, order(timecode.s))
    order.r<-with(data.aux, order(timecode.r))
    
    x.s<-x.s[order.s, ]
    y.s<-y.s[order.s]
    timecode.s<-timecode.s[order.s]
    regioncode.s<-regioncode.s[order.s]
    
    x.r<-x.r[order.r, ]
    timecode.r<-timecode.r[order.r]
    regioncode.r<-regioncode.r[order.r]
  }
  
  return(list(ww = ww, P = P))
}

# Predicting Temporal MQ-GWR

STMQ.pred <- function(y.s, x.s, x.r, timecode.r, GEOID.sr, betaSAE.STMQ, 
                      years, years.incl.t){
  
  res.modSAE <- res.ST.modSAE <- s.MAD <- pred.st.GWR <- pred.rt.GWR <-
    id.list.s.GWR <- id.list.r.GWR <- pred.rt.N.GWR <- pred.st.N.GWR <- list() 
  
  for(j in 1:dim(betaSAE.STMQ)[1]){
    
    res.modSAE[[j]] <- as.numeric(y.s - x.s %*% t(as.matrix(betaSAE.STMQ[j, -c(1,2)])))
    s.MAD[[j]] <- fun.MAD(res.modSAE[[j]])
    res.ST.modSAE[[j]] <- res.modSAE[[j]]/s.MAD[[j]]
  }
  
  for (j in 1:length(GEOID.sr)){  
    
    pred.st.N.GWR[[j]] <- pred.st.GWR[[j]] <- id.list.s.GWR[[j]] <- 
      pred.rt.N.GWR[[j]]  <- pred.rt.GWR[[j]] <- id.list.r.GWR[[j]] <- list()
    for (t in 1:T){ 
      
      select.jt <- which((betaSAE.STMQ[, 'GEOID'] == GEOID.sr[j]))[t]
      
      beta.jt <- (betaSAE.STMQ)[select.jt, -c(1,2)]
      cond.jt <- (GEOID.s==GEOID.sr[j] &  timecode.s==years[t])
      
      BC.term.jt <- s.MAD[[select.jt]]*mean(hub.psi(
        res.ST.modSAE[[select.jt]][timecode.s %in% years.incl.t[t][[1]]], k=3))
      
      pred.st.N.GWR[[j]][[t]] <- x.s[cond.jt,]%*%t(beta.jt)
      pred.st.GWR[[j]][[t]] <-  pred.st.N.GWR[[j]][[t]] + BC.term.jt
      
      pred.rt.N.GWR[[j]][[t]] <- x.r[data.aux$positions.name.min==GEOID.sr[j] & 
                                       timecode.r==years[t],]%*% t(beta.jt)
      pred.rt.GWR[[j]][[t]] <-  pred.rt.N.GWR[[j]][[t]]  + BC.term.jt 
      
      id.list.s.GWR[[j]][[t]] <- data$id[(cond.jt)]
      id.list.r.GWR[[j]][[t]] <- which(data.aux$positions.name.min==GEOID.sr[j] &  
                                         timecode.r==years[t])
    }	
  } 
  
  return(list(pred.st.N.GWR, pred.st.GWR, pred.rt.N.GWR, pred.rt.GWR,
              id.list.s.GWR, id.list.r.GWR))
}


STMQ.mse<-function(x.s, y.s, x.r, modSAE.STMQ, mod50.STMQ, 
                   regioncode.s, timecode.s, regioncode.r, timecode.r){
  
  grid <- seq(0, 10, by=.01)
  
  m <- length(unique(regioncode.s))
  T <- length(unique(timecode.s))
  p <- dim(x.s)[2]
  n <- dim(x.s)[1]
  GEOID.sr <- sort(unique(GEOID.s))
  years <- sort(unique(timecode.s))
  
  betaSAE.STMQ <- data.frame(rep(GEOID.sr, each=T), rep(years, times=length(GEOID.sr)), 
                              Reduce(rbind, lapply(modSAE.STMQ, function(x){ t(sapply(x, 
                              function(y){ return(y$'coefficients') })) })) )
  colnames(betaSAE.STMQ) <- c('GEOID', 'Year', paste('beta', 1:p, sep=''))
  
  
  beta50.STMQ <- data.frame(rep(GEOID.sr, each=T), rep(years, times=length(GEOID.sr)), 
                            Reduce(rbind, lapply(mod50.STMQ, function(x){ t(sapply(x, 
                            function(y){ return(y$'coefficients') })) })) )
  colnames(beta50.STMQ) <- c('GEOID', 'Year', paste('beta', 1:p, sep=''))
  
  
  preds <- id.order <- list()
  for (j in 1:length(GEOID.sr)){  
    preds[[j]] <- id.order[[j]] <- list()
    for (t in 1:T){ 
      cond.jt <- (GEOID.s==GEOID.sr[j] &  timecode.s==years[t])
      beta.jt <- (betaSAE.STMQ)[which((betaSAE.STMQ[, 'GEOID'] == 
                                          GEOID.sr[j]))[t], -c(1,2)]
      preds[[j]][[t]] <- x.s[cond.jt,]%*%t(beta.jt)
      id.order[[j]][[t]] <- data$id[(cond.jt)]
    }	
  } 
  
  preds <- data.frame('preds'=unlist(preds), 'id'=unlist(id.order))
  preds <- preds[order(preds$id), ]$preds
  
  
  pred.rt.N.GWR <- pred.rt.BC.GWR <- var1 <- var2 <- var3 <- 
    bias1 <- bias1BC <-  var1BC <- var2BC <- rmse1 <- rmse2 <- 
    rmse3 <- rmse1BC <- rmse2BC <- rmse3BC <- 
    id.list.r.GWR <- c.phi <- list() 
  
  for (j in 1:length(GEOID.sr)){
    
    print(j)
    pred.rt.N.GWR[[j]] <- pred.rt.BC.GWR[[j]] <- var1[[j]] <- var2[[j]] <- 
      var3[[j]] <- bias1[[j]] <- rmse1[[j]] <- rmse2[[j]] <- rmse3[[j]] <- 
      bias1BC[[j]] <- rmse1BC[[j]] <- rmse2BC[[j]] <- rmse3BC[[j]] <-
      id.list.r.GWR[[j]] <- c.phi[[j]] <-  list() 
    
    for (t in 1:T){
      if ((t == 1)|(P==0)){ incl.t = t } else { 
        incl.t = t-P+1; if (incl.t < 1){ incl.t=1 } }
      
      regioncode.s.t<-regioncode.s[timecode.s %in% years[incl.t:t]]
      timecode.s.t<-timecode.s[timecode.s %in% years[incl.t:t]]
      x.s.t<-x.s[timecode.s %in% years[incl.t:t], ]
      y.s.t<-y.s[timecode.s %in% years[incl.t:t]]
      GEOID.s.t <- GEOID.s[timecode.s %in% years[incl.t:t]]
      preds.t <- preds[timecode.s %in% years[incl.t:t]]
      
      regioncode.r.t<-regioncode.r[timecode.r %in% years[incl.t:t]]
      timecode.r.t<-timecode.r[timecode.r %in% years[incl.t:t]]
      x.r.t<-x.r[timecode.r %in% years[incl.t:t], ]
      GEOID.r.t <- GEOID.r[timecode.r %in% years[incl.t:t]]
      data.aux.r <- data.aux[timecode.r %in% years[incl.t:t],]
      
      ndt.matrix <- table(regioncode.s, timecode.s)
      nu.t<-length(timecode.s.t)
      weig.dt<-diag(as.numeric(modSAE.STMQ[[j]][[t]]$q.weights))
      
      x.r.t.j <- as.matrix(x.r.t[data.aux.r$positions.name.min[timecode.r.t %in% 
                     years[incl.t:t]]==GEOID.sr[j] & timecode.r.t==years[t], ])
      
      if(dim(x.r.t.j)[2] ==1 ){
        a.dtk <- weig.dt%*%x.s.t%*%solve(t(x.s.t)%*%weig.dt%*%x.s.t)%*% x.r.t.j
      } else if(dim(x.r.t.j)[1] ==0){ 
        next  } else {
        a.dtk <- weig.dt%*%x.s.t%*%solve(t(x.s.t)%*%weig.dt%*%x.s.t)%*% t(x.r.t.j)
      }
      
      index.r.t.k <- data.aux.r$positions.name.min==GEOID.sr[j] & timecode.r.t==years[t]
      
      var.ydtk.out <- rep(1/nu.t, nu.t) 
      
      lambda.dt <- a.dtk^2
      L <- ncol(as.matrix(a.dtk))
      
      pred.rt.N.GWR[[j]][[t]] <- t(a.dtk)%*%y.s.t
      id.list.r.GWR[[j]][[t]] <- which(data.aux$positions.name.min==GEOID.sr[j] &  
                                         timecode.r==years[t])
      
      res.d.t <- y.s.t - x.s.t%*%modSAE.STMQ[[j]][[t]]$coefficients
      s.d.t <- fun.MAD(res.d.t)
      
      var.beta <- nu.t^2*s.d.t^2/(nu.t-p) * sum(hub.psi(res.d.t/s.d.t, k=1.345)^2)/
        (sum(der.hub.psi(res.d.t/s.d.t, k=1.345))^2) * solve(t(x.s.t)%*%x.s.t)
      
      var.aux <- rep(NA, L)
      x.r.aux <- x.r[data.aux$positions.name.min==GEOID.sr[j] & 
                       timecode.r==years[t],]
      
      for(m in 1:L){
        if(ncol(lambda.dt) ==1 ){ 
          var.aux[m] <- x.r.aux %*% var.beta %*% x.r.aux
        } else 
          var.aux[m] <- x.r.aux[m, ] %*% var.beta %*% x.r.aux[m, ]
      }
      
      pos <- which(colnames(weights.df)==GEOID.sr[j])[1]
      suma.var1k = suma.var3k = rep(NA, dim(a.dtk)[2])
      
      for(k in 1:L){
        suma.var1 = suma.var3 = 0
        for(d.k in sort(unique(regioncode.s.t))){
          index <- (regioncode.s.t==d.k)
          suma.var1=suma.var1+sum(lambda.dt[index,k]*(y.s.t[index]-
                    x.s.t[index,]%*%mod50.STMQ[[pos]][[t]]$coefficients)^2)
          
          suma.var3=suma.var3+sum(var.ydtk.out[index]*(y.s.t[index]-
                    x.s.t[index,]%*%modSAE.STMQ[[pos]][[t]]$coefficients)^2)
          
        }
        suma.var1k[k] <- suma.var1
        suma.var3k[k] <- suma.var3
      }
      
      var1[[j]][[t]]  <- suma.var1k
      var2[[j]][[t]]  <- t(lambda.dt)%*%(preds.t-y.s.t)^2 
      
      bias1[[j]][[t]] <- t(a.dtk)%*%preds.t - pred.rt.N.GWR[[j]][[t]]  
      rmse1[[j]][[t]] <- sqrt( bias1[[j]][[t]]^2 + var1[[j]][[t]] )
      rmse2[[j]][[t]] <- sqrt( bias1[[j]][[t]]^2 + var2[[j]][[t]] )
      rmse3[[j]][[t]] <- sqrt( bias1[[j]][[t]]^2 + var.aux + suma.var3)
      
      pos.i <- apply(matrix(sapply(grid, function(k) {
        bias_BC <- s.d.t / nu.t * sum(hub.psi(res.d.t / s.d.t, k = k))
        var_BC <- (s.d.t / nu.t)^2 * sum(hub.psi(res.d.t / s.d.t, k = k)^2)
        
        return( (bias1[[j]][[t]] + bias_BC)^2 + var_BC )
      }),nrow=L), 1, which.min)
      
      c.phi[[j]][[t]] <- grid[pos.i]
      
      bias.BC <- varBC <- rep(NA, L)
      for(k in 1:L){
        bias.BC[k] <- s.d.t / nu.t * sum(hub.psi(res.d.t / s.d.t, k = c.phi[[j]][[t]][k]))
        varBC[k] <- (s.d.t / nu.t)^2 * sum(hub.psi(res.d.t / s.d.t, k = c.phi[[j]][[t]][k])^2)
      }
      
      pred.rt.BC.GWR[[j]][[t]] <- pred.rt.N.GWR[[j]][[t]] + bias.BC
      
      bias1BC[[j]][[t]] <- bias1[[j]][[t]] + bias.BC
      rmse1BC[[j]][[t]] <- sqrt( bias1BC[[j]][[t]]^2 + var1[[j]][[t]] + varBC)
      rmse2BC[[j]][[t]] <- sqrt( bias1BC[[j]][[t]]^2 + var2[[j]][[t]] + varBC)
      rmse3BC[[j]][[t]] <- sqrt( bias1BC[[j]][[t]]^2 + var.aux + suma.var3 + varBC)
      
    }
  }
  
  return(list(pred.rt.N.GWR = unlist(pred.rt.N.GWR),
              pred.rt.BC.GWR = unlist(pred.rt.BC.GWR),
              c.phi = unlist(c.phi),
              id.list.r.GWR = unlist(id.list.r.GWR), 
              bias=unlist(bias1), biasBC=unlist(bias1BC), 
              rmse1 = unlist(rmse1), rmse2 = unlist(rmse2), 
              rmse3 = unlist(rmse3), rmse1BC = unlist(rmse1BC),  
              rmse2BC = unlist(rmse2BC), rmse3BC = unlist(rmse3BC)
  ))
}  

