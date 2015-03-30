#' runZnorm runMeanSD
#'
#' \code{runZnorm} Compute a running Z-normalization on a numeric vector or on
#' each column of a matrix-like object. 
#'
#' 
#' runMeanSD ruturns a list of the running mean and standard deviation 
#' components in which runZnorm(x) should be equal to (x - runMeanSD(x)[[1]]) / runMeanSD(x)[[2]]
#'
#' For detailed information see
#' \url{http://www.cs.ucr.edu/~eamonn/UCRsuite.html}.
#'
#' @param x The orginal data to be transformed 1 or 2 dim non list or data
#' frame.  @param win A scaler integer for determining the sliding window size
#' @param epoch A scaler integer stating when to restart the normalizaion
#' process so we don't have numerical instability.  
#' @param backfill Setting backfill=TRUE simply takes first normalized 
#' value and back-propogates it from
#' win:1 otherwise 1:win will contain NAs.
#'
#' @return 
#' A Z-normalization vector or matrix-like object. Attributes are duplicated so it 
#' should be xts-friendly.
#' 
#'  @examples 
#'  set.seed(1001)
#'  x <- sapply(1:6,function(i) arima.sim(list(ar=c(.95)),n=500))
#'  colnames(x) <- paste0("k",1:6)
#'  
#'  mfrow <- par("mfrow")
#'  par("mfrow"=c(2,1))
#'
#'  k1 <- runZnorm(x[,1],win=10)
#'  matplot(k1,,'l')
#'
#'  k2 <- runZnorm(x,win=10)
#'  matplot(k2,,'l')
#'
#'  ms <- runMeanSD(x,win=10)
#'  all.equal(k2,((x-ms$mu)/ms$sd))
#'
#'  par("mfrow"=mfrow)
#'
#'
#' @author William Pleasant adapted from UCR Suite
#' \url{http://www.cs.ucr.edu/~eamonn/UCRsuite.html}
#'
#' @references Rakthanmanon, T., Campana, B., Mueen, A., Batista, G., Westover,
#' B., Zhu, Q., ... & Keogh, E. (2012, August). Searching and mining trillions
#' of time series subsequences under dynamic time warping. In Proceedings of the
#' 18th ACM SIGKDD international conference on Knowledge discovery and data
#' mining (pp. 262-270). ACM.
#'



runZnorm <- function(x,win,epoch=100000L, backfill=FALSE){
  if(is.list(x)) stop("x can not be a list or data.frame\n")
  .Call("run_znorm_C",x,win,epoch,backfill)
}

#' @rdname runZnorm
#'


runMeanSD <- function(x,win,epoch=100000L,backfill=FALSE){
  if(is.list(x)) stop("x can not be a list or data.frame\n")
  .Call("run_MeanSD_C",x,win,epoch,backfill)
}
