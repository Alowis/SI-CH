mexDependence <- function (x, which, dqu, margins = "laplace", constrain = TRUE, 
                           v = 10, maxit = 1e+06, start = c(0.01, 0.01), marTransform = "mixture", 
                           nOptim = 1, PlotLikDo = FALSE, PlotLikRange = list(a = c(-1, 
                                                                                    1), b = c(-3, 1)), PlotLikTitle = NULL) 
{
  theCall <- match.call()
  if (class(x) != "migpd") 
    stop("you need to use an object created by migpd")
  margins <- list(casefold(margins), p2q = switch(casefold(margins), 
                                                  gumbel = function(p) -log(-log(p)),
                                                  frechet = function(p) -1/log(p),
                                                  laplace = function(p) ifelse(p < 
                                                                                 0.5, log(2 * p), -log(2 * (1 - p)))), q2p = switch(casefold(margins), 
                                                                                                                                    gumbel = function(q) exp(-exp(-q)),  frechet = function(p) -1/log(p), laplace = function(q) ifelse(q < 0, exp(q)/2, 1 - 0.5 * exp(-q))))
  x <- mexTransform(x, margins = margins, method = marTransform)
  if (margins[[1]] == "gumbel" & constrain) {
    warning("With Gumbel margins, you can't constrain, setting constrain=FALSE")
    constrain <- FALSE
  }
  if (missing(which)) {
    message("Missing 'which'. Conditioning on", dimnames(x$transformed)[[2]][1], 
            "\n")
    which <- 1
  }
  else if (length(which) > 1) 
    stop("which must be of length 1")
  else if (is.character(which)) 
    which <- match(which, dimnames(x$transformed)[[2]])
  if (missing(dqu)) {
    warning("Assuming same quantile for dependence thesholding as was used\n     to fit corresponding marginal model...\n")
    dqu <- x$mqu[which]
  }
  dth <- quantile(x$transformed[, which], dqu)
  dependent <- (1:(dim(x$data)[[2]]))[-which]
  if (length(dqu) < length(dependent)) 
    dqu <- rep(dqu, length = length(dependent))
  aLow <- ifelse(margins[[1]] == "gumbel", 10^(-10), -1 + 10^(-10))
  if (missing(start)) {
    start <- c(0.01, 0.01)
  }
  else if (class(start) == "mex") {
    start <- start$dependence$coefficients[1:2, ]
  }
  if (length(start) == 2) {
    start <- matrix(rep(start, length(dependent)), nrow = 2)
  }
  if (length(start) != 2 * length(dependent)) {
    stop("start should be of type 'mex' or be a vector of length 2, or be a matrix with 2 rows and ncol equal to the number of dependence models to be estimated")
  }
  if (!missing(PlotLikRange)) {
    PlotLikDo <- TRUE
  }
  qfun <- function(X, yex, wh, aLow, margins, constrain, v, 
                   maxit, start) {
    Qpos <- function(param, yex, ydep, constrain, v, aLow) {
      a <- param[1]
      b <- param[2]
      res <- PosGumb.Laplace.negProfileLogLik(yex, ydep, 
                                              a, b, constrain, v, aLow)
      res$profLik
    }
    o <- try(optim(par = start, fn = Qpos, control = list(maxit = maxit), 
                   yex = yex[wh], ydep = X[wh], constrain = constrain, 
                   v = v, aLow = aLow), silent = TRUE)
    if (class(o) == "try-error") {
      warning("Error in optim call from mexDependence")
      o <- as.list(o)
      o$par <- rep(NA, 6)
      o$value <- NA
    }
    else if (o$convergence != 0) {
      warning("Non-convergence in mexDependence")
      o <- as.list(o)
      o$par <- rep(NA, 6)
    }
    else if (nOptim > 1) {
      for (i in 2:nOptim) {
        o <- try(optim(par = o$par, fn = Qpos, control = list(maxit = maxit), 
                       yex = yex[wh], ydep = X[wh], constrain = constrain, 
                       v = v, aLow = aLow), silent = TRUE)
        if (class(o) == "try-error") {
          warning("Error in optim call from mexDependence")
          o <- as.list(o)
          o$par <- rep(NA, 6)
          o$value <- NA
          (break)()
        }
        else if (o$convergence != 0) {
          warning("Non-convergence in mexDependence")
          o <- as.list(o)
          o$par <- rep(NA, 6)
          (break)()
        }
      }
    }
    if (PlotLikDo) {
      nGridPlotLik <- 50
      a.grid <- seq(PlotLikRange$a[1], PlotLikRange$a[2], 
                    length = nGridPlotLik)
      b.grid <- seq(PlotLikRange$b[1], PlotLikRange$b[2], 
                    length = nGridPlotLik)
      NegProfLik <- matrix(0, nrow = nGridPlotLik, ncol = nGridPlotLik)
      for (i in 1:nGridPlotLik) {
        for (j in 1:nGridPlotLik) {
          NegProfLik[i, j] <- PosGumb.Laplace.negProfileLogLik(yex = yex[wh], 
                                                               ydep = X[wh], a = a.grid[i], b = b.grid[j], 
                                                               constrain = constrain, v = v, aLow = aLow)$profLik
        }
      }
      NegProfLik[NegProfLik > 10^10] <- NA
      if (sum(!is.na(NegProfLik))) {
        filled.contour(a.grid, b.grid, -NegProfLik, main = paste("Profile likelihood", 
                                                                 PlotLikTitle), color.palette = terrain.colors, 
                       xlab = "a", ylab = "b", plot.axes = {
                         axis(1)
                         axis(2)
                         points(o$par[1], o$par[2])
                       })
      }
    }
    if (!is.na(o$par[1])) {
      if (margins == "gumbel" & o$par[1] <= 10^(-5) & o$par[2] < 
          0) {
        lo <- c(10^(-10), -Inf, -Inf, 10^(-10), -Inf, 
                10^(-10))
        Qneg <- function(yex, ydep, param) {
          param <- param[-1]
          b <- param[1]
          cee <- param[2]
          d <- param[3]
          m <- param[4]
          s <- param[5]
          obj <- function(yex, ydep, b, cee, d, m, s) {
            mu <- cee - d * log(yex) + m * yex^b
            sig <- s * yex^b
            log(sig) + 0.5 * ((ydep - mu)/sig)^2
          }
          res <- sum(obj(yex, ydep, b, cee, d, m, s))
          res
        }
        o <- try(optim(c(0, 0, 0, 0, 0, 1), Qneg, method = "L-BFGS-B", 
                       lower = lo, upper = c(1, 1 - 10^(-10), Inf, 
                                             1 - 10^(-10), Inf, Inf), yex = yex[wh], ydep = X[wh]), 
                 silent = TRUE)
        if (class(o) == "try-error" || o$convergence != 
            0) {
          warning("Non-convergence in mexDependence")
          o <- as.list(o)
          o$par <- rep(NA, 6)
        }
      }
      else {
        Z <- (X[wh] - yex[wh] * o$par[1])/(yex[wh]^o$par[2])
        o$par <- c(o$par[1:2], 0, 0, mean(Z), sd(Z))
      }
    }
    c(o$par[1:6], o$value)
  }
  yex <- c(x$transformed[, which])
  wh <- yex > unique(dth)
  res <- sapply(1:length(dependent), function(X, dat, yex, 
                                              wh, aLow, margins, constrain, v, maxit, start) qfun(dat[, 
                                                                                                      X], yex, wh, aLow, margins, constrain, v, maxit, start[, 
                                                                                                                                                             X]), dat = as.matrix(x$transformed[, dependent]), yex = yex, 
                wh = wh, aLow = aLow, margins = margins[[1]], constrain = constrain, 
                v = v, maxit = maxit, start = start)
  loglik <- -res[7, ]
  res <- matrix(res[1:6, ], nrow = 6)
  dimnames(res)[[1]] <- c(letters[1:4], "m", "s")
  dimnames(res)[[2]] <- dimnames(x$transformed)[[2]][dependent]
  gdata <- as.matrix(x$transformed[wh, -which])
  tfun <- function(i, data, yex, a, b, cee, d) {
    data <- data[, i]
    a <- a[i]
    b <- b[i]
    cee <- cee[i]
    d <- d[i]
    if (is.na(a)) 
      rep(NA, length(data))
    else {
      if (a < 10^(-5) & b < 0) 
        a <- cee - d * log(yex)
      else a <- a * yex
      (data - a)/(yex^b)
    }
  }
  z <- try(sapply(1:(dim(gdata)[[2]]), tfun, data = gdata, 
                  yex = yex[wh], a = res[1, ], b = res[2, ], cee = res[3, 
                                                                       ], d = res[4, ]))
  if (class(z) %in% c("Error", "try-error")) {
    z <- matrix(nrow = 0, ncol = dim(x$data)[[2]] - 1)
  }
  else if (is.R()) {
    if (!is.array(z)) {
      z <- matrix(nrow = 0, ncol = dim(x$data)[[2]] - 1)
    }
  }
  dimnames(z) <- list(NULL, dimnames(x$transformed)[[2]][dependent])
  res2 <- list(coefficients = res, Z = z, dth = unique(dth), 
               dqu = unique(dqu), which = which, conditioningVariable = colnames(x$data)[which], 
               loglik = loglik, margins = margins, constrain = constrain, 
               v = v)
  oldClass(res2) <- "mexDependence"
  output <- list(margins = x, dependence = res2, call = theCall)
  oldClass(output) <- "mex"
  output
}




mexTransform<-function( x ,
                        method = "mixture",
                        divisor = "n+1",
                        na.rm = TRUE,
                        margins = "laplace"){
  
  if ( !is.element( method, c( "mixture", "empirical" ) ) )
    stop( "method should be either 'mixture' or 'empirical'" )
  
  if ( divisor == "n" )
    divisor <- dim( x$data )[[ 1 ]]
  else if ( divisor == "n+1" )
    divisor <- dim( x$data )[[ 1 ]] + 1
  else stop( "divisor can be 'n' or 'n+1'" )
  
  transFun <- function( i , x , mod , th, divisor, method ){
    x <- x[ , i ]
    mod <- mod[[ i ]]
    th <- th[ i ]
    
    ox <- order( x )
    names( x ) <- 1:length( x )
    
    x <- sort( x )
    run <- rle( x )
    p <- cumsum( run$lengths )/ divisor
    p <- rep( p, run$lengths )
    names(p)
    p <- p[ order(as.character( names( x )))]
    x <- x[ order( as.character( names( x ) ) ) ]
    sort(p)
    
    Femp <- p
    if ( method == "mixture" ){
      sigma <- exp( mod$coefficients[ 1 ] )
      xi <- mod$coefficients[ 2 ]
      
      Para <- ( 1 + xi * ( x - th) / sigma ) ^ ( -1 / xi )
      a<-x>th
      sum(a)
      Para <- 1 - mean( x > th ) * Para
      res <- ifelse( x <= th , Femp , Para )
    }
    else res <- Femp
    
    res[ ox ] <- sort( res )
    res
  } # Close transfun
  
  res <- sapply( 1:( dim( x1$data )[[ 2 ]] ), transFun,
                 x = x1$data, mod = x1$models, th = x1$mth,
                 divisor = divisor, method=method
  )
  dimnames( res ) <- list( NULL, names( x$models ) )
  
  if (margins == "gumbel"){
    x$transformed <- -log( -log( res ) )
  }
  else if (margins == "laplace"){
    #x$transformed <- sign(res - .5) * log(1 - 2*abs(res - .5))
    
    x$transformed <- ifelse(res < .5,  log(2 * res), -log(2 * (1 - res) ))
  }
  else if (margins == "frechet"){ 
    
    x$transformed<- (-log(res))^(-1)}
  else { stop("margins need to be gumbel, frechet or laplace") }
  invisible(x)
}


