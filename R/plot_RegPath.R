#' @title plot_RegPath
#' @description a function that produces a coefficient profile plot of the coefficient paths given lambda
#'
#' @param beta_hat an p*m estimated coefficient matrix
#' @param lambda an m*1 regularization parameter vector
#'
#'
#' @author Rongqian Zhang
#'
#' @return a image showing the coefficient paths given lambda
#'
#' @examples
#' data(longley)
#' lambda=seq(from=0,to=0.1,by=0.01)
#' beta_hat=matrix(nrow=6,ncol = length(lambda))
#' for (i in 1:length(lambda))
#' {
#' beta_hat[,i]=RidgeReg(scale(longley[,2:7]),longley$y,lambda = lambda[i])
#' }
#' plot_RegPath(lambda,beta_hat)
#'
#' @importFrom("graphics","plot","points")
#' @export


plot_RegPath<-function(lambda,beta_hat)
{
  plot(lambda,beta_hat[1,],type = 'l',ylim = c(min(beta_hat)-0.1,max(beta_hat)+0.1),ylab = 'Coefficients',xlab = 'Lambda')
for(i in 2:nrow(beta_hat))
{
  points(lambda,beta_hat[i,],type = 'l',col=i)
}
}


