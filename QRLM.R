
###############################################
###    External functions from MQ models    ###
###############################################


### HUBER FUNCTION

hub.psi <- function(x, k){ ifelse(abs(x) <= k, x, sign(x) * k) }
der.hub.psi <- function(x, k){ ifelse(abs(x) <= k, 1, 0) }


### MEDIAN ABSOLUTE DEVIATION

fun.MAD<-function(x){ median(abs( x-median(x)) )/0.6745 }


### QUANTILE REGRESSION FITTING ALGORITHM - IRLS

QRLM <-function (x, y, case.weights = rep(1, nrow(x)), k=1.345, var.weights = rep(1, nrow(x)), ...,
                 w = rep(1, nrow(x)), init = "ls", psi = psi.huber, scale.est = c("MAD", "Huber", "proposal 2"), 
                 k2 = 1.345, method = c("M", "MM"), maxit = 20, acc = 1e-04, test.vec = "resid", q = 0.5)
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
  
  qest <- matrix(0, nrow = ncol(x), ncol = length(q))
  qwt <- matrix(0, nrow = nrow(x), ncol = length(q))
  qfit <- matrix(0, nrow = nrow(x), ncol = length(q))
  qres <- matrix(0, nrow = nrow(x), ncol = length(q))
  qvar.matrix <- array(rep(0,ncol(x)*ncol(x)),dim=c(ncol(x),ncol(x),length(q)))
  qscale <- NULL
  for(i in 1:length(q)) {
    for (iiter in 1:maxit) {
      if (!is.null(test.vec))
        testpv <- get(test.vec)
      if (scale.est != "MM") {
        if (scale.est == "MAD")
          scale <- median(abs(resid/sqrt(var.weights)))/0.6745
        else {gamma<- 4*k2^2*(1-pnorm(k2))*((1-q[i])^2+q[i]^2) - 4*k2*dnorm(k2)*((1-q[i])^2+q[i]^2) + 
          4*(1-q[i])^2*(pnorm(0)-(1-pnorm(k2))) + 4*q[i]^2*(pnorm(k2)-pnorm(0))
        scale <- sqrt(sum(pmin(resid^2/var.weights,(k2*scale)^2))/(n1*gamma))
        }
        if (scale == 0) {
          done <- TRUE
          break
        }
      }
      w <- psi(resid/(scale * sqrt(var.weights)),k=k) * case.weights
      ww <- 2 * (1 - q[i]) * w
      ww[resid > 0] <- 2 * q[i] * w[resid > 0]
      w <- ww
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
    qscale[i]<-scale
    qwt[, i] <- w
    qfit[, i] <- temp$fitted.values
    qres[,i] <- resid
    
    tmp.res.mq<-qres[,i]/qscale[i]
    Epsi2<-(sum((qwt[,i]*tmp.res.mq)^2)/(nrow(x)-ncol(x)))
    Epsi<-(1/qscale[i])*(sum(2*(q[i]*(0<=tmp.res.mq & tmp.res.mq<= k)+
                                  (1-q[i])*(-k <=tmp.res.mq & tmp.res.mq<0)))/nrow(x))
    qvar.matrix[,,i]<- (((Epsi2)/Epsi^2)*solve(t(x)%*%x))
    
  }
  list(fitted.values = qfit, residuals = qres, q.values = q, q.weights = qwt, coef= qest,
       qscale=qscale,var.beta=qvar.matrix)
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

