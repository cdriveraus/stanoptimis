#' Optimize / importance sample a stan model.
#'
#' @param standata list object conforming to rstan data standards.
#' @param sm compiled stan model object.
#' @param init init argument conforming to rstan init standards.
#' @param deoptim Do first pass optimization using differential evolution? Slower, but better for cases with multiple
#' minima / difficult optimization.
#' @param decontrol List of control parameters for differential evolution step, to pass to \code{\link[DEoptim]{DEoptim.control}}.
#' @param cores Number of cpu cores to use.
#' @param isloops Number of iterations of adaptive importance sampling to perform after optimization.
#' @param isloopsize Number of samples per iteration of importance sampling.
#' @param issamples Number of samples to use for final results of importance sampling.
#' @param adjtransform Logical. If TRUE, the MAP estimate used as the basis for sampling will be based on the *unconstrained*
#' parameters. If FALSE, it will use the estimate for *constrained* parameters.
#' @param verbose Logical. Extra output while fitting?
#'
#' @return list containing: optimfit, the result from optimizing with stans optimizer;
#' stanfit, the stanfit object used;
#' posterior, the sampled posterior;
#' transformedpars, the full model output based on the posterior;
#' transformedpars_old, a simple summary approach based on the minimum and Hessian;
#' isdiags, a list containing elements used by \code{\link{isdiag}} for importance sampling diagnostics.
#' @export
#' @import DEoptim
#'
#' @examples
#' #' \dontrun{
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
#' fit <- sampling(sm, iter = 10000)
#' summary(fit)$summary
#'
#' ## extract samples as a list of arrays
#' e <- extract(fit, permuted = TRUE)
#'
#' optimis <- stanoptimis(standata = list(),sm = sm,isloops=30,issamples = 3000,cores=1)
#'
#'
#' apply(optimis$posterior,2,mean)
#' apply(optimis$posterior,2,sd)
#' isdiag(optimis)
#'
#' plot(density(optimis$posterior))
#' points(density(e$y),type='l',col=2)
#' }
stanoptimis <- function(standata, sm, init=0, verbose=FALSE,adjtransform=TRUE,
  deoptim=FALSE,
  decontrol=list(),
  isloops=5, isloopsize=500, issamples=500, cores=1){

  if(is.null(decontrol$steptol)) decontrol$steptol=5
  if(is.null(decontrol$reltol)) decontrol$reltol=1e-4
  if(is.null(decontrol$NP)) decontrol$NP='auto'
  if(is.null(decontrol$CR)) decontrol$CR=.9
  if(is.null(decontrol$trace)) decontrol$trace =as.integer(verbose)


  message('Optimizing...')
  betterfit<-TRUE
  try2 <- FALSE
  while(betterfit){ #repeat loop if importance sampling improves on optimized max
    betterfit <- FALSE

    #create fit object for use with rstan functions
    suppressWarnings(suppressOutput(smf<-sampling(sm,iter=1,chains=1,data=standata,check_data=FALSE, control=list(max_treedepth=0))))

    npars=get_num_upars(smf)

    if(deoptim){ #init with DE
      if(decontrol$NP=='auto') NP=min(c(40,10*npars)) else NP = decontrol$NP

      decontrollist <- c(decontrol,DEoptim.control())
      decontrollist <- decontrollist[unique(names(decontrollist))]

      lp2 = function(parm) {
        out<-try(log_prob(smf, upars=parm,adjust_transform=adjtransform,gradient=FALSE),silent = TRUE)
        if(class(out)=='try-error') {
          out=-1e20
        }
        return(-out)
      }
      deinit <- matrix(rnorm(npars*NP,0,2),nrow = NP)
      deinit[2,] <- 0
      if(length(init)>1) {
        deinit[1,] <- unconstrain_pars(smf,init)
        if(NP > 10) deinit[3:9,] =  matrix( rnorm(npars*(7),rep(deinit[1,],each=7),.1), nrow = 7)
      }
      decontrollist$initialpop=deinit
      decontrollist$NP = NP
      optimfitde <- suppressWarnings(DEoptim(fn = lp2,lower = rep(-1e10, npars), upper=rep(1e10, npars),
        control = decontrollist))
      init=constrain_pars(object = smf,optimfitde$optim$bestmem)
    }

    suppressWarnings(suppressOutput(optimfit <- optimizing(sm,standata, hessian=FALSE, iter=40000, init=init,as_vector=FALSE,
      tol_obj=1e-12, tol_rel_obj=0,init_alpha=.001, tol_grad=0,tol_rel_grad=1e1,tol_param=1e-12,history_size=100),verbose=verbose))


    estc=optimfit$par
    bestfit <-optimfit$value
    estu=unconstrain_pars(smf, estc)



    lp<-function(parm) {
      out<-try(log_prob(smf, upars=parm,adjust_transform=adjtransform,gradient=FALSE),silent = TRUE)
      if(class(out)=='try-error') {
        out=-Inf
      }
      return(out)
    }

    grf<-function(parm,...) {
      out=try(grad_log_prob(smf, upars=parm, adjust_transform = adjtransform))
      if(class(out)=='try-error') {
        out=rep(NA,length(parm))
      }
      return(out)
    }

    grmat<-function(func,pars,step=1e-8){
      gradout<-matrix(NA,nrow=length(pars),ncol=length(pars))
      for(i in 1:length(pars)){
        stepsize <- step * 10
        while(any(is.na(gradout[i,])) && stepsize > 1e-14){
          stepsize <- stepsize * .1
          uppars<-pars
          downpars<-pars
          uppars[i]<-pars[i]+stepsize
          downpars[i]<-pars[i]-stepsize
          gradout[i,]<-((func(uppars)) - (func(downpars)))/stepsize/2
        }
      }
      return(gradout)
    }

    # A more numerically stable way of calculating log( sum( exp( x ))) Source:
    # http://r.789695.n4.nabble.com/logsumexp-function-in-R-td3310119.html
    log_sum_exp <- function(x) {
      xmax <- which.max(x)
      log1p(sum(exp(x[-xmax] - x[xmax]))) + x[xmax]
    }


    hess=grmat(func=grf,pars=estu)
    if(any(is.na(hess))) stop(paste0('Hessian could not be computed for pars ', which(apply(hess,1,function(x) any(is.na(x)))), ' -- consider reparameterising.'))
    hess = (hess/2) + t(hess/2)
    mchol=try(t(chol(solve(-hess))),silent=TRUE)
    if(class(mchol)=='try-error') message('Hessian not positive-definite -- check importance sampling convergence with isdiag')
    mcov=MASS::ginv(-hess)
    mcov=as.matrix(Matrix::nearPD(mcov)$mat)

    mcovl <- list()
    mcovl[[1]]=mcov
    delta=list()
    delta[[1]]=estu
    samples <-c()
    resamples <- c()
    prop_dens <-c()
    target_dens<-c()
    sample_prob<-c()
    counter <- 0
    ess <- 0
    qdiag<-0

    cl <- parallel::makeCluster(cores, type = "PSOCK")
    parallel::clusterExport(cl, c('sm','standata'),environment())

    if(isloops == 0) {
      nresamples = issamples
      resamples <- matrix(unlist(lapply(1:5000,function(x){
        delta[[1]] + t(chol(mcovl[[1]])) %*% t(matrix(rnorm(length(delta[[1]])),nrow=1))
      } )),byrow=TRUE,ncol=length(delta[[1]]))
      message('Importance sampling not done -- interval estimates via Hessian based sampling only')
    }

    if(isloops > 0){
      message('Adaptive importance sampling, loop:')
      for(j in 1:isloops){
        message(paste0('  ', j, ' / ', isloops, '...'))
        if(j==1){
          df=2
          samples <- mvtnorm::rmvt(isloopsize, delta = delta[[j]], sigma = mcovl[[j]],   df = df)
        } else {
          # if(j>5) df <- 3
          delta[[j]]=colMeans(resamples)
          mcovl[[j]] = cov(resamples)+diag(1e-6,ncol(samples))
          samples <- rbind(samples,mvtnorm::rmvt(isloopsize, delta = delta[[j]], sigma = mcovl[[j]],   df = df))
        }
        prop_dens <- mvtnorm::dmvt(tail(samples,isloopsize), delta[[j]], mcovl[[j]], df = df)

        parallel::clusterExport(cl, c('samples'),environment())

        target_dens[[j]] <- unlist(parallel::parLapply(cl, parallel::clusterSplit(cl,1:isloopsize), function(x){
          eval(parse(text=paste0('library(rstan)')))

          #have to regenerate stan fit object... inefficient
          smf<-sampling(sm,iter=1,chains=1,data=standata,check_data=FALSE,control=list(max_treedepth=0,adapt_engaged=FALSE))

          lp<-function(parm) { #recreate function in here because scoping
            out<-try(log_prob(smf, upars=parm, adjust_transform = adjtransform, gradient=FALSE),silent = TRUE)
            if(class(out)=='try-error') {
              out=-Inf
            }
            return(out)
          }
          out <- apply(tail(samples,isloopsize)[x,],1,lp)

          #unload old rstan dlls (attempt to reduce crashes)
          try(dyn.unload(file.path(tempdir(), paste0(smf@stanmodel@dso@dso_filename, .Platform$dynlib.ext))),silent = TRUE)

          return(out)

        }))

        if(all(target_dens[[j]] < -1e29)) stop('Could not sample from optimum! Try reparamaterizing?')
        if(any(target_dens[[j]] > bestfit && (j < isloops && !try2))){
          oldfit <- bestfit
          try2 <- TRUE
          bestfit<-max(target_dens[[j]],na.rm=TRUE)
          betterfit<-TRUE
          init = rstan::constrain_pars(object = smf, samples[which(unlist(target_dens) == bestfit),])
          message('Improved fit found - ', bestfit,' vs ', oldfit,' - restarting optimization')
          break
        }
        nresamples = ifelse(j==isloops,issamples,5000)

        target_dens2 <- target_dens[[j]] + (0-max(target_dens[[j]])) #adjustment to get in decent range
        target_dens2[!is.finite(target_dens[[j]])] <- -1e30
        weighted_dens <- target_dens2 - prop_dens

        sample_prob <- c(sample_prob,exp((weighted_dens - log_sum_exp(weighted_dens)))) #sum to 1 for each iteration, normalise later
        sample_prob[!is.finite(sample_prob)] <- 0
        resample_i <- sample(1:nrow(samples), size = nresamples, replace = ifelse(j == isloops+1,FALSE,TRUE),
          prob = sample_prob / sum(sample_prob))
        resamples <- samples[resample_i, , drop = FALSE]

        ess[j] <- (sum(sample_prob[resample_i]))^2 / sum(sample_prob[resample_i]^2)
        qdiag[j]<-mean(unlist(lapply(sample(x = 1:length(sample_prob),size = 500,replace = TRUE),function(i){
          (max(sample_prob[resample_i][1:i])) / (sum(sample_prob[resample_i][1:i]) )
        })))

      }
    }
  }#end after isloops when no better fit found by IS

  if(isloops==0) lpsamples <- NA else lpsamples <- unlist(target_dens)[resample_i]


  message('Computing quantities...') #at least, I think this is what is taking so long...

  relistarrays <- function(flesh, skeleton){
    skelnames <- names(skeleton)
    skelstruc <- lapply(skeleton,dim)
    count=1
    npars <- length(flesh)
    out <- list()
    for(ni in skelnames){
      if(!is.null(skelstruc[[ni]])){
        out[[ni]] <- array(flesh[count:(count+prod(skelstruc[[ni]]))],dim = skelstruc[[ni]])
        count <- count + prod(skelstruc[[ni]])
      } else {
        out[[ni]] <- flesh[count]
        count <- count + 1
      }
    }
    return(out)
  }


  parallel::clusterExport(cl, c('relistarrays','resamples','sm','standata','optimfit'),environment())

  transformedpars <- parallel::parLapply(cl, parallel::clusterSplit(cl,1:nresamples), function(x){
    Sys.sleep(.1)
    smf<-sampling(sm,iter=1,chains=1,data=standata,check_data=FALSE,control=list(max_treedepth=0,adapt_engaged=FALSE))
    Sys.sleep(.1)

    out <- list()
    for(li in 1:length(x)){
      Sys.sleep(.01)
      flesh = unlist(rstan::constrain_pars(smf, resamples[x[li],]))
      names(flesh) <- c()
      skeleton=optimfit$par
      out[[li]] <-relistarrays(flesh, skeleton)
    }
    return(out)
  })
  parallel::stopCluster(cl)
  transformedpars<-unlist(transformedpars,recursive = FALSE)



  tostanarray <- function(flesh, skeleton){
    skelnames <- names(skeleton)
    skelstruc <- lapply(skeleton,dim)
    count=1
    npars <- ncol(flesh)
    niter=nrow(flesh)
    out <- list()
    for(ni in skelnames){
      if(prod(skelstruc[[ni]])>0){
        if(!is.null(skelstruc[[ni]])){
          out[[ni]] <- array(flesh[,count:(count+prod(skelstruc[[ni]])-1)],dim = c(niter,skelstruc[[ni]]))
          count <- count + prod(skelstruc[[ni]])
        } else {
          out[[ni]] <- array(flesh[,count],dim = c(niter))
          count <- count + 1
        }
      }
    }
    return(out)
  }

  transformedpars=tostanarray(flesh = matrix(unlist(transformedpars),byrow=TRUE, nrow=nresamples), skeleton = optimfit$par)

  #compute simplistic sd's
  sds=try(suppressWarnings(sqrt(diag(mcov))))
  if(class(sds)=='try-error') sds <- rep(NA,length(estu))
  lest= estu - 1.96 * sds
  uest= estu + 1.96 * sds

  transformedpars_old=NA
  try(transformedpars_old<-cbind(unlist(constrain_pars(smf, lest)),
    unlist(constrain_pars(smf, estu)),
    unlist(constrain_pars(smf, uest))),silent=TRUE)
  try(colnames(transformedpars_old)<-c('2.5%','mean','97.5%'),silent=TRUE)

  stanoptimfit=list(optimfit=optimfit,stanfit=smf, posterior = resamples, transformedpars=transformedpars,transformedpars_old=transformedpars_old,
    isdiags=list(cov=mcovl,means=delta,ess=ess,qdiag=qdiag,lpsamples=lpsamples ))
  return(stanoptimfit)
}
