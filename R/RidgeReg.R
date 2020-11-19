#' @title Ridge Regression
#' @description a function that computes the ridge regression coefficient estimates given X, y and lambda.
#'
#' @param X an n*p design matrix
#' @param y an n*1 response vector
#' @param lambda a pre-specified tuning parameter
#'
#' @author Rongqian Zhang
#'
#' @return an p*1 ridge regression coefficients vector
#'
#' @examples
#' #example 1 n>p
#' set.seed(123)
#' X=matrix(rnorm(100*80),nrow=100,ncol=80)
#' beta=rt(80,df=2)
#' y=X%*%beta+rnorm(100,sd=2)
#' RidgeReg(X,y,1)
#' plot(beta)
#' points(RidgeReg(X,y,1),col='red')
#'
#' #example 2 n<p
#' set.seed(123)
#' X=matrix(rnorm(100*80),nrow=80,ncol=100)
#' beta=rt(100,df=2)
#' y=X%*%beta+rnorm(80,sd=2)
#' RidgeReg(X,y,1)
#' plot(beta)
#' points(RidgeReg(X,y,1),col='red')
#'
#' @export
#'

RidgeReg<-function(X,y,lambda)
{
  num_x<-dim(X)
  if (num_x[1]>=num_x[2])
  {
    res<-svd.lm.fit1(X,y,lambda)
  }else{
    res<-svd.lm.fit2(X,y,lambda)
  }
  return(res)
}

svd.lm.fit1 <- function(x,y,z){
  res_svd <- svd(x,LINPACK=TRUE)
  tuy <- crossprod(res_svd$u,y)
  de <- res_svd$d/(res_svd$d*res_svd$d + rep(z,length(res_svd$d)))

  return(res_svd$v%*%(tuy*de))
}

svd.lm.fit2 <- function(x,y,z){
  res_svd <- svd(t(x),LINPACK=TRUE)
  tuy <- crossprod(res_svd$v,y)
  de <- res_svd$d/(res_svd$d*res_svd$d + rep(z,length(res_svd$d)))

  return(res_svd$u%*%(tuy*de))
}

