#' runMinMax
#'
#' \code{runMinMax} Compute the running Min Max or Both on a matrix-like object. 
#'
#' An extremely fast method of calculating running min and max windows in C. 
#'
#'
#' @param x a vector or matrix-like object.   
#' @param win A scaler integer for determining the window size.
#' @param dropAttr A bool which if set to false attempts to copy appropriate attributes.  
#'
#' @examples 
#'  x <- sapply(1:2,function(i) arima.sim(list(ar=c(.95)),n=500))
#'  colnames(x) <- paste0("k",1:2)
#'
#'  k1 <- runmin(x[,-1,drop=FALSE],win=10,dropAttr=FALSE)
#'  str(k1)
#'  k2 <- runmax(x[,-1,drop=FALSE],win=10,dropAttr=FALSE)
#'  str(k2)
#'  matplot(cbind(x[,1],k1,k2),,'l')
#'
#'  # multivariate
#'  k3 <- runmin(x,win=10,dropAttr=FALSE)
#'  k4 <- runmax(x,win=10,dropAttr=FALSE)
#'
#'  matplot(cbind(x,k3,k4),,'l')
#'
#'
#' @references
#' \url{http://richardhartersworld.com/cri/2001/slidingmin.html}.
#'


runMinMax <- function(x,win,dropAttr=FALSE){
  if(missing(win)) stop("missing win\n")
  .Call(run_min_max_C, x, win,3L,dropAttr)
}

#'
#' @rdname runMinMax
#'

runmin <- function(x,win,dropAttr=FALSE){
  if(missing(win)) stop("missing win\n")
  .Call(run_min_max_C, x, win,1L,dropAttr)
}

#'
#' @rdname runMinMax
#'

runmax <- function(x,win,dropAttr=FALSE){
  if(missing(win)) stop("missing win\n")
  .Call(run_min_max_C, x, win,2L,dropAttr)
}



