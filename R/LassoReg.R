#' @title LASSO Regression
#' @description a function that computes the LASSO regression coefficient estimates given X, y and lambda.
#'
#' @param X an n*p design matrix
#' @param y an n*1 response vector
#' @param lambda a pre-specified tuning parameter
#' @param soft logical. If TRUE (the default) use soft-thresholding operator
#' @param thresh convergence threshold for coordinate descent
#' @param max_iter max number of iterations for coordinate descent
#'
#' @author Rongqian Zhang
#'
#' @return an p*1 LASSO regression coefficients vector
#'
#' @examples
#' #example 1 n>p
#' set.seed(123)
#' X=matrix(rnorm(100*12),nrow=100,ncol=12)
#' beta=c(rnorm(2,2,1),rnorm(3,-1,1),rep(0,7))
#' y=X%*%beta+rnorm(100,sd=1)
#' beta_lasso<-LassoReg(X,y,lambda=0.01,thresh = 1e-7,max_iter = 5000,soft = FALSE)
#' plot(beta)
#' points(beta_lasso,col='red')
#'
#'
#'
#' #example 2 n<p
#' set.seed(123)
#' X=matrix(rnorm(10*20),nrow=10,ncol=20)
#' beta=c(rnorm(2,2,1),rnorm(3,-1,1),rep(0,7),rnorm(4,3,2),rnorm(1,-4,1),rep(0,3))
#' y=X%*%beta+rnorm(10,sd=0.5)
#' beta_lasso<-LassoReg(X,y,lambda=.1,thresh = 1e-7,max_iter = 5000,soft = FALSE)
#' plot(beta)
#' points(beta_lasso,col='red')


#'
#' @export
#'

LassoReg <- function(
  X,                   # model matrix
  y,                   # target
  lambda  = .1,        # penalty parameter
  soft    = FALSE,      # soft vs. hard thresholding
  thresh     = 1e-6,      # tolerance
  max_iter    = 10000     # number of max iterations
  ) {

  # soft thresholding function
  soft_thresh <- function(a, b) {
    out = rep(0, length(a))
    out[a >  b] = a[a > b] - b
    out[a < -b] = a[a < -b] + b
    out
  }

  w = solve(crossprod(X) + diag(lambda, ncol(X))) %*% crossprod(X,y)
  thresh_curr = 1
  J = ncol(X)
  a = rep(0, J)
  c_ = rep(0, J)
  i = 1

  while (thresh < thresh_curr && i < max_iter) {
    w_old = w
    a = colSums(X^2)
    l = length(y)*lambda  # for consistency with glmnet approach
    c_ = sapply(1:J, function(j)  sum( X[,j] * (y - X[,-j] %*% w_old[-j]) ))
    if (soft) {
      for (j in 1:J) {
        w[j] = soft_thresh(c_[j]/a[j], l/a[j])
      }
    }else {
      w = w_old
      w[c_< l & c_ > -l] = 0
    }

   thresh_curr = crossprod(w - w_old)
    i = i + 1
    #if (verbose && i%%10 == 0) message(i)
  }
  if (i == max_iter) message('the algorithm cannot reach convergence within max_iter, please increase max_iter')
  return(w)
}

