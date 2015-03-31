#' LB_Keogh
#'
#' \code{lbKeogh} Compute the univariate or multivariate LB_Keogh on a
#' matrix-like object. 
#'
#' A lower bound for Dynamic Time Warping introduced by Keogh and
#' Ratanamahatana.  lbKeogh is a c-level utility that returns object of
#' different dimensions depending to the number of columns in x and the value of
#' returnSum.  If y has one dimension and x is a matrix-like
#' object then we compute the lbKeogh on each column in x using a window (win)
#' created on y. If y more than one column then the the ncols(y)==ncols(x) and a
#' lbKeogh will computed on each column of x with an envelope created on the matching
#' column of y.  With returnSum=FALSE we return an object with the same
#' dimensions as x otherwise we return a numeric vector length equal the number of 
#' columns in x.  The default setting for dropAttr is FALSE.  Setting drapAttr=TRUE
#' attempts to duplicate all of the objects, if returnSum=FALSE or just the colnames 
#' with returnSum=TRUE.
#'
#' The distance used is Euclidean without the sqrt() calculation ie (x-y)^2. To
#' get the true distance simply take sqrt(lbKeogh(x,y,returnSum=TRUE)) or for the
#' multivariate version sqrt(sum(lbKeogh(x,y,returnSum=TRUE))).
#'
#' envelope calculates the upper and lower bounds used in the lbKeogh. The
#' implementation idea was introduced by Danial Lemire in 2009. It returns a
#' list containing the upper and lower envelopes for each column in x.
#'
#' For detailed information see:
#'
#' \url{http://www.cs.ucr.edu/~eamonn/LB_Keogh.htm}.
#'
#' \url{http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.14.6347&rep=rep1&type=pdf}.
#' 
#'  @examples 
#'  set.seed(1001)
#'  x <- sapply(1:6,function(i) arima.sim(list(ar=c(.95)),n=500))
#'  colnames(x) <- paste0("k",1:6)
#'
#'  k1 <- lbKeogh(x[,-1,drop=FALSE],x[,1],win=10,dropAttr=FALSE)
#'  matplot(k1,,'l')
#'  str(k1)
#'  k2 <- lbKeogh(x[,-1,drop=FALSE],x[,1],win=10,returnSum=TRUE)
#'
#'  all.equal(colSums(k1),k2,check.attributes=FALSE)
#'
#'  # multivariate
#'  k3 <- lbKeogh(x[,1:3,drop=FALSE],x[,4:6],win=10,dropAttr=FALSE)
#'  matplot(k3,,'l')
#'  str(k3)
#'
#'  k4 <- lbKeogh(x[,1:3,drop=FALSE],x[,4:6],win=10,returnSum=TRUE)
#'  matplot(k4,,'l')
#'  str(k4)
#'
#'  all.equal(colSums(k3),k4,check.attributes=FALSE)

#'
#'
#' @author William Pleasant adapted from UCR Suite
#' \url{http://www.cs.ucr.edu/~eamonn/UCRsuite.html}
#'
#' @param x The reference data which can be a 1 or 2 dimensional. x must have the 
#'  same number of rows as y.
#'   vector or matrix-like object.   
#' @param y The query data which can be a 1 or 2 dimensional 
#'   vector or matrix-like object, see details.
#' @param win A scaler integer for determining the window size for the upper 
#'   and lower bounds.
#' @param returnSum  A bool which if TRUE sums the columns of the reference set x.
#' @param cSum A bool to determine if the column of the return value should be
#' cumulatively summed. cSum is over ridden if retunSum is TRUE. The only
#' expected use for cSum=TRUE is for the early abandoning of later index
#' searches.
#' @param dropAttr A bool which if set to false attempts to copy appropriate attributes.  
#'
#' @references
#' Keogh, E. (2002). Exact indexing of dynamic time warping. In 28th
#' International Conference on Very Large Data Bases. Hong Kong. pp 406-417.
#'
#' Rath, Toni M., and R. Manmatha. "Lower-Bounding of Dynamic Time Warping
#' Distances for Multivariate Time Series." 
#'
#' "Faster Retrieval with a Two-Pass Dynamic-Time-Warping Lower Bound", Pattern
#' Recognition 42(9), 2009.


lbKeogh <- function(x,y,win=20,returnSum=FALSE,cSum=FALSE,dropAttr=TRUE){
  .Call("lbKeogh_C",x,y,as.integer(win),returnSum,cSum,dropAttr)
}

#'
#' @rdname lbKeogh
#'

envelope <- function(x,win){
  if(missing(win)) stop("missing win\n")
  .Call("lower_upper_env_C",x,win)
}



