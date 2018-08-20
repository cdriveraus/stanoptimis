#' Diagnostics for importance sampling from stanoptim
#'
#' @param fit Output from stanoptimis when isloops > 0
#'
#' @return Nothing. Plots convergence of parameter mean estimates from initial Hessian based distribution to final sampling distribution.
#' @export
#'
#' @examples
#' \dontrun{
#' library(rstan)
#' scode <- "
#' parameters {
#'   real y[2];
#' }
#' model {
#'   y[1] ~ normal(0, 1);
#'   y[2] ~ double_exponential(0, 2);
#' }
#' "
#'
#' sm <- stan_model(model_code=scode)
#'
#' fit1 <- sampling(object = sm, iter = 10, verbose = FALSE)
#' print(fit1)
#' fit2 <- stan(fit = fit1, iter = 10000, verbose = FALSE)
#'
#' ## extract samples as a list of arrays
#' e2 <- extract(fit2, permuted = TRUE)
#'
#' ## using as.array on the stanfit object to get samples
#' a2 <- as.array(fit2)
#'
#'
#' optimfit <- optimstan(standata = list(),sm = sm,isloops=10,issamples = 1000,cores=3)
#'
#' apply(optimfit$posterior,2,mean)
#' apply(optimfit$posterior,2,sd)
#' isdiag(optimfit)
#'
#' plot(density(optimfit$posterior))
#' points(density(e2$y))
#' }

isdiag <- function(fit,wait=TRUE){
  iter=length(fit$isdiags$cov)
  mcov <- fit$isdiags$cov
  samplecov <- cov(fit$posterior)
  means <- simplify2array(fit$isdiags$means)
  means <- (means - means[,ncol(means)])
  means <- t(means / sqrt(diag(samplecov)))

  # smeans <- matrix(apply(means,2,function(x) x))),byrow=TRUE,ncol=iter)
matplot(means,type='l',main='Mean convergence',xlab='Sampling loop',ylab=' Z divergence relative to finish',xlim=c(0,iter*1.2))

if(wait) readline(prompt = 'Press a key to go to next plot')

legend('topright',bty='n',legend = paste0('par',1:ncol(means)),lty = 1:5,col=1:6,text.col=1:6,cex = .7)

sds <- simplify2array(mcov)
sds <- apply(sds,3,function(x) sqrt(diag(x)))
sds <- t((sds - sds[,iter]) / sqrt(diag(samplecov)))
matplot(sds,type='l',main='SD convergence',xlab='Sampling loop',ylab='Z divergence relative to finish',xlim=c(0,iter*1.2))

legend('topright',bty='n',legend = paste0('par',1:ncol(means)),lty = 1:5,col=1:6,text.col=1:6,cex = .7)


}
