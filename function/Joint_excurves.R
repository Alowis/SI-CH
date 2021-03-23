
JointExceedanceCurve <- function(Sample, ExceedanceProb,...) {
  theCall <- match.call()
  UseMethod("JointExceedanceCurve", Sample)
}

#' @rdname JointExceedanceCurve
#' @export
JointExceedanceCurve.default <- function(Sample, ExceedanceProb,n=50,x=NULL,...) {
  s <- calcJointExceedanceCurve(Sample,ExceedanceProb,n,x)
  names(s) <- names(Sample)
  attr(s,"ExceedanceProb") <- ExceedanceProb
  s
}

#' @rdname JointExceedanceCurve
#' @export
#' @method JointExceedanceCurve mexMC
#' @param which Vector length two identifying which margins to use for joint exceedance curve estimation. Can be integer vector, giving column numbers of original data matrix, or character vector identifying variables by name (these must match column names in original data).
#' @param ... Further aguments to be passed to methods
JointExceedanceCurve.mexMC <- function(Sample, ExceedanceProb,n=50,x=NULL,which=1:2,...) {
  S <- Sample$MCsample[,which]
  Sample <- as.matrix(S)
  Sample <- Sample[!is.na(Sample[,1]),]
  Sample <- Sample[!is.na(Sample[,2]),]
  
  s <- calcJointExceedanceCurve(Sample,ExceedanceProb,n,x)
  names(s) <- names(S)
  attr(s,"ExceedanceProb") <- ExceedanceProb
  s
}

#' @rdname JointExceedanceCurve
#' @export
#' @method JointExceedanceCurve predict.mex
JointExceedanceCurve.predict.mex <- function(Sample, ExceedanceProb,n=200,x=NULL,boots=F,which=1:2,meth="J",...) {
  CondExceedanceProb <- 1-Sample$pqu
  if(is.numeric(which)){ # since column order of predict order sample is not original data column order
    d <- dim(Sample$data$simulated)[2]
    w <- Sample$which

    which <- order(c(w,c(1:d)[-w]))[which]
  }
  if (boots==T){
    st<-list()
    for (ids in (1:length(Sample$replicates))){
    S <- Sample$replicate[[ids]][,which]
    Sample1 <- as.matrix(S)
    if(ExceedanceProb > CondExceedanceProb) stop("ExceedanceProb must be less than the probability of exceeding the threshold used for importance sampling in the call to predict")
    
    s <- calcJointExceedanceCurve(Sample1,ExceedanceProb/CondExceedanceProb,n,x)
    names(s) <- names(S)
    attr(s,"ExceedanceProb") <- ExceedanceProb
    st<-c(st,list(s))
    }
    return(st)
  }

  if(boots==F){
  S <- Sample$data$simulated[,which]
  Sample <- as.matrix(S)
  condex<-CondExceedanceProb
  if(ExceedanceProb > CondExceedanceProb) stop("ExceedanceProb must be less than the probability of exceeding the threshold used for importance sampling in the call to predict")
  if (meth=="J"){
  s <- calcJointExceedanceCurve(Sample,ExceedanceProb/CondExceedanceProb,n,x)}
  if (meth=="C"){
  s <- calcCondExceedanceCurve(Sample,CondExceedanceProb,ExceedanceProb/CondExceedanceProb,n,x,w)}
  names(s) <- names(S)
  attr(s,"ExceedanceProb") <- ExceedanceProb
  s
  }

}
  # Sample <- as.matrix(S)
  # if(ExceedanceProb > CondExceedanceProb) stop("ExceedanceProb must be less than the probability of exceeding the threshold used for importance sampling in the call to predict")
  # 
  # s <- calcJointExceedanceCurve(Sample,ExceedanceProb/CondExceedanceProb,n,x)
  # names(s) <- names(S)
  # attr(s,"ExceedanceProb") <- ExceedanceProb
  # s


calcCondExceedanceCurve  <- function(Sample,condex,ExceedanceProb,n=50,x=NULL,w) {
  # mx, my are marginal upper limits
  # px, py are plotting points
  
  mx <- quantile(Sample[,1],1-ExceedanceProb)
  my <- quantile(Sample[,2],1-ExceedanceProb)
  if(is.null(x)){
    px <- seq(min(Sample[,1],na.rm=T),mx,length=n)
    py <- seq(min(Sample[,2],na.rm=T),my,length=n)
  } else {
    px <- x
  }
  mxx <- quantile(Sample[,w],1-ExceedanceProb)
  pxx<-seq(min(Sample[,w]),mxx,length=n)
  as<-Sample[,w]
  if(w==1)condexprob<-1-condex+(sapply(pxx,function(pxx){(mean(as < pxx))})*condex)
  if(w==2)condexprob<-sapply(pxx,function(pxx){(mean(as > pxx))})
  ixd<-seq(1,length(pxx))
  
  f <- function(z,x,y) {
    sapply(z,function(z){
      g <- function(w) mean(x > z & y > w)/(1-condex+mean(x<z)*condex) - ExceedanceProb
      if(g(range(y)[1]) * g(range(y)[2]) < 0) {
        out <- uniroot(g,lower=range(y)[1],upper=range(y)[2])$root
      } else {
        out <- min(y[x>z])
      }
      out
    }
    )
  }
  
  #calculate curve values at plotting points
  if(w==1)cx <- f(px,Sample[,1],Sample[,2])
  if(w==2)cx <- f(px,Sample[,2],Sample[,1])
  # print(w)
  # if(is.null(x)){
  #   cy <- f(py,Sample[,2],Sample[,1])
  #   res <- data.frame(x=c(px,cy),y=c(cx,py))
  #   # sort results
  #   res <- res[order(res[,1]),]
  #   row.names(res) <- NULL
  # } else {
  if(w==1)res <- data.frame(x=px,y=cx)
  if(w==2)res <- data.frame(x=px,y=cx)
  # }
  oldClass(res) <- "jointExcCurve"
  res
}

calcJointExceedanceCurve  <- function(Sample,ExceedanceProb,n=50,x=NULL) {
  # mx, my are marginal upper limits
  # px, py are plotting points
  
  mx <- quantile(Sample[,1],1-ExceedanceProb)
  my <- quantile(Sample[,2],1-ExceedanceProb)
  if(is.null(x)){
    px <- seq(min(Sample[,1],na.rm=T),mx,length=n)
    py <- seq(min(Sample[,2],na.rm=T),my,length=n)
  } else {
    px <- x
  }
  
  f <- function(z,x,y) {
    sapply(z,function(z){
      g <- function(w) mean(x > z & y > w) - ExceedanceProb
      if(g(range(y)[1]) * g(range(y)[2]) < 0) {
        out <- uniroot(g,lower=range(y)[1],upper=range(y)[2])$root
      } else {
        out <- min(y[x>z])
      }
      out
    }
    )
  }
  
  #calculate curve values at plotting points
  cx <- f(px,Sample[,1],Sample[,2])
  if(is.null(x)){
    cy <- f(py,Sample[,2],Sample[,1])
    res <- data.frame(x=c(px,cy),y=c(cx,py))
    # sort results
    res <- res[order(res[,1]),]
    row.names(res) <- NULL
  } else {
    res <- data.frame(x=px,y=cx)
   }
  oldClass(res) <- "jointExcCurve"
  res
}

print.jointExcCurve <- function(x, ...){
  P <- attributes(x)$ExceedanceProb
  cat("\n Estimated curve with constant joint exceedance probability equal to",P,"\n")
  res <- as.data.frame(cbind(x[[1]],x[[2]]))
  colnames(res) <- attributes(x)$names
  print(res)
  invisible(x)
}


geom_jointExcCurve <- function(x,...){
  dat <- as.data.frame(cbind(x[[1]],x[[2]]))
  colnames(dat) <- attributes(x)$names
  geom_line(data=dat,...)
}



