

#based on rstan function, very cut down, may fail in some cases...
#' @importFrom Rcpp cpp_object_initializer
getcxxfun <- function(object) {
  if (length(object@dso_saved) == 0){
    return(function(...) stop("this function should not be called"))
  }  else  return(object@.CXXDSOMISC$cxxfun)
}


stan_reinitsf <- function(model, data,fast=FALSE){
  if(fast) sf <- new(model@mk_cppmodule(model),data,0L,getcxxfun(model@dso))

  if(!fast) suppressMessages(suppressWarnings(suppressOutput(sf<-sampling(model,iter=0,chains=0,init=0,data=data,check_data=FALSE,
    control=list(max_treedepth=0),save_warmup=FALSE,test_grad=FALSE))))

  return(sf)
}



flexsapply <- function(cl, X, fn,cores=1){
  if(cores > 1) parallel::parSapply(cl,X,fn) else sapply(X, fn)
}


stan_constrainsamples<-function(sm,standata, samples,cores=2){
  smfull <- stan_reinitsf(model = sm,data = standata)
  message('Computing quantities for ', nrow(samples),' samples...')
  est1=NA
  class(est1)<-'try-error'
  i=0
  while(class(est1)=='try-error'){
    i=i+1
    est1=try(constrain_pars(smfull, upars=samples[i,]),silent=TRUE)
  }
  if(class(est1)=='try-error') stop('All samples generated errors! Respecify, try stochastic optimizer, try again?')

  cl2 <- parallel::makeCluster(cores, type = "PSOCK")
  on.exit(parallel::stopCluster(cl2))
  parallel::clusterExport(cl2, c('sm','standata','samples','est1'),environment())
  parallel::clusterApply(cl2,1:cores,function(x) require(ctsem))

  transformedpars <- try(flexsapply(cl2, parallel::clusterSplit(cl2,1:nrow(samples)), function(x){ #could pass smaller samples
    if(!is.null(standata$savescores) && !standata$savescores) standata$dokalmanrows <- as.integer(c(1,standata$subject[-1] - standata$subject[-standata$ndatapoints]))
    smfull <- stan_reinitsf(sm,standata)
    out <- list()
    skeleton=est1
    for(li in 1:length(x)){
      out[[li]] <- try(constrain_pars(smfull, upars=samples[x[li],]))
      if(any(sapply(out[[li]], function(x) any(c(is.nan(x),is.infinite(x),is.na(x)))))) class(out[[li]]) <- c(class(out[[li]]),'try-error')
    }
    return(out)
  },cores=cores))
  transformedpars=unlist(transformedpars,recursive = FALSE)
  missingsamps <-sapply(transformedpars, function(x) 'try-error' %in% class(x))
  nasampscount <- sum(missingsamps)

  transformedpars <- transformedpars[!missingsamps]
  nresamples <- nrow(samples) - nasampscount
  if(nasampscount > 0) {
    message(paste0(nasampscount,' NAs generated during final sampling of ', nrow(samples), '. Biased estimates may result -- consider importance sampling, respecification, or full HMC sampling'))
  }
  #this seems inefficient and messy, should be a better way...
  transformedpars=try(tostanarray(flesh=matrix(unlist(transformedpars),byrow=TRUE, nrow=nresamples), skeleton = est1))

  return(transformedpars)
}



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





#' Optimize / importance sample a stan or ctStan model.
#'
#' @param standata list object conforming to rstan data standards.
#' @param sm compiled stan model object.
#' @param init vector of unconstrained parameter values, or character string 'random' to initialise with
#' random values very close to zero.
#' @param initsd positive numeric specifying sd of normal distribution governing random sample of init parameters,
#' if init='random' .
#' @param sampleinit either NA, or an niterations * nparams matrix of samples to initialise importance sampling.
#' @param deoptim Do first pass optimization using differential evolution? Slower, but better for cases with multiple
#' minima / difficult optimization.
#' @param stochastic Logical. Use stochastic gradient descent instead of mize (bfgs) optimizer.
#' Still experimental, worth trying for either robustness checks or problematic, high dimensional, nonlinear, problems.
#' @param plotsgd Logical. If TRUE, plot iteration details when using stochastic optimizer.
#' @param estonly if TRUE,just return point estimates under $rawest subobject.
#' @param verbose Integer from 0 to 2. Higher values print more information during model fit -- for debugging.
#' @param decontrol List of control parameters for differential evolution step, to pass to \code{DEoptim.control}.
#' @param tol objective tolerance.
#' @param cores Number of cpu cores to use.
#' @param is Logical. Use importance sampling, or just return map estimates?
#' @param isloopsize Number of samples of approximating distribution per iteration of importance sampling.
#' @param finishsamples Number of samples to draw for final results of importance sampling.
#' @param finishmultiply Importance sampling stops once available samples reach \code{finishsamples * finishmultiply , then the final samples are drawn
#' without replacement from this set.
#' @param tdf degrees of freedom of multivariate t distribution. Higher (more normal) generally gives more efficent
#' importance sampling, at risk of truncating tails.
#' @param chancethreshold drop iterations of importance sampling where any samples are chancethreshold times more likely to be drawn than expected.
#'
#' @return list containing fit elements
#' @importFrom mize mize
#' @importFrom utils head tail
#' @importFrom Rcpp evalCpp
#' @import rstan
#' @export
#' @examples
#' #' \donttest{
#'
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
#' #for ml or map estimates
#' optimis <- stanoptimis(standata = list(),sm = sm,finishsamples = 3000,cores=2)
#' optimis$optimfit
#'
#' #for posterior distributions
#' optimis <- stanoptimis(standata = list(),sm = sm,finishsamples = 3000,cores=6,tdf=2)
#'
#' apply(optimis$rawposterior,2,mean)
#' apply(optimis$rawposterior,2,sd)
#' isdiag(optimis)
#'
#' plot(density(optimis$rawposterior[,2],bw=.05))
#' points(density(e$y[,2],bw=.05),type='l',col=2)
#' }
stanoptimis <- function(standata, sm, init='random',initsd=.01,sampleinit=NA,
  deoptim=FALSE, estonly=FALSE,tol=1e-14,
  decontrol=list(),
  stochastic = FALSE, #'auto',
  plotsgd=FALSE,
  is=TRUE, isloopsize=1000, finishsamples=500, tdf=2,chancethreshold=100,finishmultiply=5,
  verbose=0,cores=2){

  standata$verbose=as.integer(verbose)

  if(is.null(decontrol$steptol)) decontrol$steptol=5
  if(is.null(decontrol$reltol)) decontrol$reltol=1e-4
  if(is.null(decontrol$NP)) decontrol$NP='auto'
  if(is.null(decontrol$CR)) decontrol$CR=.9
  if(is.null(decontrol$trace)) decontrol$trace =ifelse(verbose>0,1,0)




  sgd <- function(init,stepbase=1e-4,gmeminit=ifelse(is.na(startnrows),.9,.9),gmemmax=.95,maxparchange = .5,
    startnrows=NA,gsmoothness = 1,roughnessmemory=.95,groughnesstarget=.5,lproughnesstarget=.1,
    gsmoothroughnesstarget=.2,
    warmuplength=50,
    minparchange=1e-20,maxiter=50000,nconvergeiter=20, itertol=1e-3, deltatol=1e-5){
    pars=init
    bestpars = pars
    maxpars=pars
    minpars=pars
    changepars=pars
    step=rep(stepbase,length(init))
    g=try(smf$grad_log_prob(upars=init,adjust_transform=TRUE),silent = TRUE) #rnorm(length(init),0,.001)
    if(class(g)=='try-error') {
      i = 0
      message('Problems initialising, trying random values...')
      while(i < 50 && class(g)=='try-error'){
        if(i %%5 == 0) init = rep(0,length(init))
        init=init+rnorm(length(init),0,abs(init)+ .1)
        g=try(smf$grad_log_prob(upars=init,adjust_transform=TRUE),silent = TRUE)
        i = i + 1
      }
    }
    g=sign(g)#*sqrt(abs(g))
    gsmooth=g
    oldg=g
    groughness = rep(groughnesstarget,length(g))
    gsmoothroughness = rep(gsmoothroughnesstarget,length(g))
    deltasmoothsq=.01
    lproughness=lproughnesstarget
    gmemory <- gmeminit
    oldgmemory <- gmemory
    oldlpdif <- 0
    lpdif <- 0
    maxlp <- -Inf
    i=0
    lp<-c()
    oldlp <- -Inf
    converged <- FALSE
    nrows <- ifelse(is.na(startnrows),standata$ndatapoints, min(standata$ndatapoints, startnrows))

    while(!converged && i < maxiter){
      print
      i = i + 1
      accepted <- FALSE
      lproughnesstarget2 = ifelse(nrows==standata$ndatapoints,lproughnesstarget,.49)
      while(!accepted){
        newpars = bestpars
        delta =   step  * (gsmooth*(gsmoothness) + g*(1-(gsmoothness))) * exp((rnorm(length(g),0,.01)))
        delta[abs(delta) > maxparchange] <- maxparchange*sign(delta[abs(delta) > maxparchange])
        newpars = newpars + delta

        #sub sampling
        if(!is.na(startnrows) || (nrows!=standata$ndatapoints)){
          subjects <- sample(1:standata$ndatapoints,nrows,replace = FALSE)
          standata$dokalmanrows <- as.integer(standata$subject %in% subjects) #rep(1L,standata$ndatapoints) #
          # smf<-stan_reinitsf(sm,standata,fast=TRUE)
        }

        # lpg = try(smf$log_prob(upars=newpars,adjust_transform=TRUE,gradient=TRUE),silent = FALSE)
        lpg= -neglpgf(newpars)
        if(class(lpg) !='try-error' && !is.nan(lpg[1]) && all(!is.nan(attributes(lpg)$gradient))) accepted <- TRUE else step <- step * .1
        if(is.na(startnrows) && i < warmuplength && i > 1 && lpg[1] < lp[1]) {
          accepted <- FALSE
          if(plotsgd) print('not accepted!')
          # step = step * .5
          gsmooth=gsmooth*.5
        }
      }
      lp[i]=lpg[1]
      pars <- newpars

      oldg=g
      g=attributes(lpg)$gradient
      g=sign(g)*sqrt(abs(g))
      oldgsmooth = gsmooth
      gmemory2 = gmemory * min(i/warmuplength,1)^(1/8)
      gsmooth= gsmooth*gmemory2 + (1-gmemory2)*g#^2 #should it really be squared? sgd algorithms do so
      roughnessmemory2 = roughnessmemory * min(i/warmuplength,1)^(1/8)

      stdgdifold = (g-oldg) * step
      stdgdifsmooth = (g-gsmooth) * step
      groughness = groughness * (roughnessmemory2) + (1-(roughnessmemory2)) * as.numeric(sign(g)!=sign(oldg))
      gsmoothroughness = gsmoothroughness * (roughnessmemory2) + (1-(roughnessmemory2)) * as.numeric(sign(gsmooth)!=sign(oldgsmooth))
      if(i > 1) lproughness = lproughness * (roughnessmemory2) + (1-(roughnessmemory2)) * as.numeric(lp[i-1] >= (lp[i] + sd(tail(lp,min(i,3)))))

      # print(stdgdif)
      # step=exp(mean(log(step))+(.99*(log(step)-mean(log(step)))))
      # step[oldsigng == signg] = step[which(oldsigng == signg)] * sqrt(2-gmemory) #exp((1-gmemory)/2)
      # step[oldsigng != signg] = step[which(oldsigng != signg)] / sqrt(2-gmemory) #ifelse(nrows == standata$ndatapoints, (2-gmemory),1.1) #1.2 #exp((1-gmemory)/2)

      signdifmod = step
      signdifmod[sign(oldg) == sign(g)] =  .1 #/ (1.5-inv_logit(abs(stdgdif[oldsigng == signg])))^4 #* (1/ ( ( (roughness*.05+.95)^2) ))
      signdifmod[sign(oldg) != sign(g)]  = -.1 #10* ((1.5-inv_logit(abs(stdgdifold[sign(oldg) != sign(g)])))-1) #* ( ( (roughness*.05+.95)^2) )
      signdifmod[is.nan(signdifmod)] <- .5 #oldstep[is.nan(step)] #because of overflow in some cases

      deltasmoothsq = deltasmoothsq * gmemory + (1-gmemory)*delta^2
      lproughnessmod= 2 * ( ( (1/(-lproughness-lproughnesstarget2)) / (1/-lproughnesstarget2) + .5) -1) #balanced eq for any centre / target
      gmemoryupd = min(gmemmax,max(.1,gmemory /  ( (1/(-mean(groughness)-groughnesstarget)) / (1/-groughnesstarget) + .5) ))
      gsmoothroughnessmod =  .5 *(( ( (1/(-(gsmoothroughness)-gsmoothroughnesstarget)) / (1/-gsmoothroughnesstarget) + .5) ) -1)
      groughnessmod = .5 *( ( ( (1/(-(groughness)-groughnesstarget)) / (1/-groughnesstarget) + .5) ) -1)
      rmsstepmod = sqrt(abs(gsmooth+1e-7))/step -1 #like adagrad but with decaying gradient

      step = (step
        # + sqrt(deltasmoothsq)/abs(gsmooth) /2
        # + step
        + step*signdifmod #* min(sqrt(deltasmoothsq),1)
        + step*lproughnessmod
        + step* gsmoothroughnessmod #* min(sqrt(deltasmoothsq),1)
        + step* groughnessmod# * min(sqrt(deltasmoothsq),1)
        # + step * rmsstepmod
      )

      # gmemory = gmemory * roughnessmemory2 + (1-roughnessmemory2)*gmemoryupd

      if(lp[i] >= max(lp)) {
        # step = step * sqrt(2-gmemory) #exp((1-gmemory)/8)
        if(i > warmuplength/2) {
          gsmooth[pars>maxpars | pars < minpars] <- gsmooth[pars>maxpars | pars < minpars]  + .2*delta[pars>maxpars | pars < minpars] /step[pars>maxpars | pars < minpars]
          step[pars>maxpars | pars < minpars] <- step[pars>maxpars | pars < minpars] * 1.5  #+ pars[pars>maxpars | pars < minpars]
          changepars=pars
          changepars[!(pars>maxpars | pars < minpars)] <- NA
          lproughness = lproughness * .9
        }
        # pars <- newpars
        bestpars <- pars

        maxpars[pars>maxpars] <-pars[pars>maxpars]
        minpars[pars<minpars] <-pars[pars<minpars]
      }
      if(i > 1 && lp[i] < lp[i-1]) {
        # step[sign(oldgsmooth) != signg] = .5 * step[sign(oldgsmooth) != signg]
        # gsmooth[sign(oldgsmooth) != signg] = .05 * g[sign(oldgsmooth) != signg]
        # gmemory=min(.8,gmemory)
        # signg <- oldsigng
        # gsmooth = oldgsmooth
        # pars <- bestpars
        # if(nrows == standata$ndatapoints) {
        # gsmooth[sign(gsmooth) != signg] = gsmooth[sign(gsmooth) != signg] #* .5 + g[sign(gsmooth) != signg] * .5
        # step = step/ (2-gmemory)#step  / max( 1.5, (-10*(lp[i] - lp[i-1]) / sd(head(tail(lp,20),10)))) #exp((1-gmemory)/4)
        # } else {
        #   step = step / 1.06
        # }
      }
      # if(i %%10 ==0) gmemory = min(gmemory+.1,gmemmax)# * exp(mean(sign(diff(tail(lp,20)))))
      # if(i %%20 ==0) gmemory =  max(gmeminit, min(gmemmax, 1.6*(1-(log(sd(tail(lp,20)) ) -log(itertol)) / (log(sd(head(lp,20)))-log(itertol)))* (1-gmeminit) + gmeminit))

      if(i > 30 && i %% 20 == 0) {
        lpdif <- sum(diff(tail(lp,10)))
        oldlpdif <- sum(diff(head(tail(lp,10),20)))
        if(oldlpdif >= lpdif) gmemory <- oldgmemory
        proposal = gmemory*2-oldgmemory
        gmemory <- min(gmemmax, max(0, proposal + runif(1,-.05,.1)))
        oldgmemory <- gmemory
      }

      step[step > maxparchange] <- maxparchange
      step[step < minparchange] <- minparchange

      if(plotsgd){
        par(mfrow=c(2,3),mgp=c(2,.8,0),mar=c(2,3,1,0)+.2)
        plot(pars)
        points(changepars,pch=17,col='red')

        plot(log(step))

        plot(groughness,col='red',ylim=c(0,1))
        abline(h=mean(gsmoothroughness),col='blue',lty=2)
        abline(h=(gsmoothroughnesstarget),col='blue',lty=1,lwd=2)
        points(gsmoothroughness,ylim=c(0,1),col='blue')
        abline(h=mean(groughness),col='red',lty=2)
        # abline(h=(groughnesstarget),col='red',lty=1)

        abline(h=lproughnesstarget,lty=1,col='green')
        abline(h=lproughness, col='green',lty=2)

        plot(tail(log(-(lp-max(lp)-1)),500),type='l')
        plot(gsmooth,ylim= c(-max(abs(gsmooth)),max(abs(gsmooth))))

        matplot(cbind(signdifmod,gsmoothroughnessmod),col=c('black','blue'),pch=1,ylim=c(-1,1))
        points(groughnessmod,col='red')
        abline(h=lproughnessmod,col='green')

        # message(paste0('Iter = ',i, '   Best LP = ', max(lp),'   grad = ', sqrt(sum(g^2)), '   gmem = ', gmemory))
      }

      #check convergence
      if(i > 30){
        if(max(tail(lp,nconvergeiter)) - min(tail(lp,nconvergeiter)) < itertol) converged <- TRUE
        # print(max(tail(lp,nconvergeiter)) - min(tail(lp,nconvergeiter)))
        if(max(diff(tail(lp,nconvergeiter))) < deltatol) converged <- TRUE
        if(nrows < standata$ndatapoints && (length(lp) - match(min(lp),lp)) > nconvergeiter) converged <- TRUE
      }
      if(converged & nrows != standata$ndatapoints){
        converged <- FALSE
        nrows <- min(standata$ndatapoints, nrows * 4)
        if(nrows > standata$ndatapoints/2){
          nrows <- standata$ndatapoints
          i=0
          lp=c()
        }
        message(paste0('nrows now ',nrows, ' out of ',standata$ndatapoints),' total')

      }
    }
    return(list(itervalues = lp, value = max(lp),par=bestpars) )
  }








  message('Optimizing...')

  betterfit<-TRUE
  bestfit <- -9999999999
  try2 <- FALSE
  while(betterfit){ #repeat loop if importance sampling improves on optimized max
    betterfit <- FALSE

    smfull <- stan_reinitsf(sm,standata)
    smf <- stan_reinitsf(sm,standata,fast=TRUE)
    npars=smf$num_pars_unconstrained()

    if(all(init %in% 'random')) init <- rnorm(npars, 0, initsd)
    if(all(init == 0)) init <- rep(0,npars)

    if(is.na(sampleinit[1])){

      if(deoptim){ #init with DE
        if(requireNamespace('DEoptim',quietly = TRUE)) {
          if(decontrol$NP=='auto') NP=min(c(40,10*npars)) else NP = decontrol$NP

          decontrollist <- c(decontrol,DEoptim::DEoptim.control())
          decontrollist <- decontrollist[unique(names(decontrollist))]

          lp2 = function(parm) {
            out<-try(smf$log_prob(upars=parm,adjust_transform=TRUE,gradient=FALSE),silent = TRUE)
            if(class(out)=='try-error') {
              out=-1e200
            }
            return(-out)
          }

          deinit <- matrix(rnorm(npars*NP,0,2),nrow = NP)
          deinit[2,] <- rnorm(npars,0,.0002)
          if(length(init)>1 & try2) {
            deinit[1,] <- unconstrain_pars(smfull,init)
            if(NP > 10) deinit[3:9,] =  matrix( rnorm(npars*(7),rep(deinit[1,],each=7),.1), nrow = 7)
          }
          decontrollist$initialpop=deinit
          decontrollist$NP = NP
          optimfitde <- suppressWarnings(DEoptim::DEoptim(fn = lp2,lower = rep(-1e10, npars), upper=rep(1e10, npars),
            control = decontrollist))
          # init=constrain_pars(object = smf,optimfitde$optim$bestmem)
          init=optimfitde$optim$bestmem
          bestfit = -optimfitde$optim$bestval
        } else stop(paste0('use install.packages(\"DEoptim\") to use deoptim')) #end require deoptim
      }

      gradout <- c()
      bestlp <- -Inf

      neglpgf<-function(parm) {
          out<-try(smf$log_prob(upars=parm,adjust_transform=TRUE,gradient=TRUE),silent = FALSE)
        if(class(out)=='try-error' || is.nan(out)) {
          out=-99999999
          gradout <<- rep(NaN,length(parm))
          attributes(out) <- list(gradient=rep(0,length(parm)))
        } else {
          if(out[1] > bestlp) {
            bestlp <<- out[1]
            gradout <<- attributes(out)$gradient
          }
        }
          if(verbose > 0) message('target = ', out)
        return(-out)
      }


      parbase=par()


      # require(mize)
      mizelpg=list(
        fg=function(pars){
          r=neglpgf(pars)
          r=list(fn=r[1],gr= -attributes(r)$gradient)
          return(r)
        },
        fn=neglpgf,
        gr=function(pars) -attributes(neglpgf(pars))$gradient
      )


      if(stochastic=='auto' && npars > 100){
        message('> 100 parameters and stochastic="auto" so stochastic gradient descent used -- try disabling if slow!')
        stochastic <- TRUE
      } else if(stochastic=='auto') stochastic <- FALSE

      if(!stochastic) {
        # optimfit <- ucminf(init,fn = neglpgf,gr = grffromlp,control=list(grtol=1e-99,xtol=tol,maxeval=10000),hessian=2)

        optimfit <- mize(init, fg=mizelpg, max_iter=99999,
          # method = 'NAG', nest_q = .01, #nest_convex_approx=TRUE,
          method="L-BFGS",memory=100,
          # method='SR1',
          line_search='Schmidt',c1=1e-10,c2=.9,step0='schmidt',ls_max_fn=999,
          abs_tol=tol,grad_tol=0,rel_tol=0,step_tol=0,ginf_tol=0)
        optimfit$value = optimfit$f
        init = optimfit$par
        bestfit <- -optimfit$value
        optimfit$value <- -optimfit$value

      }

      if(stochastic || is.infinite(bestfit)){
        if(is.infinite(bestfit)) {
          message('Switching to stochastic optimizer -- failed initialisation with bfgs')
        }
        optimfit <- sgd(init)
        bestfit <-optimfit$value
      }


      est2=optimfit$par #unconstrain_pars(smf, est1)
    }

    if(!estonly){

      cl <- parallel::makeCluster(cores, type = "PSOCK") #should vary this depending on os for memory management
      parallel::clusterExport(cl = cl, varlist = c('sm','smf','standata'),envir = environment())
      parallel::clusterApply(cl,1:cores,function(x) {
        library(Rcpp)
        library(stanoptimis)
        })
      on.exit(parallel::stopCluster(cl))

      lpg<-function(parm) {
        out<-try(smf$log_prob(upars=parm,adjust_transform=TRUE,gradient=TRUE),silent = FALSE)
        if(class(out)=='try-error' || is.nan(out)) {
          out=Inf
          gradout <<- rep(NaN,length(parm))
        } else {
          gradout <<- attributes(out)$gradient
        }
        return(out)
      }

      grmat<-function(pars,step=1e-5,lpdifmin=1e-8, lpdifmax=1e-3){

        hessout <- flexsapply(cl=cl, X = 1:length(pars), function(i) {
        # hessout <- sapply(X = 1:length(pars), function(i) {
          smf <- stan_reinitsf(sm,standata,fast=TRUE)
          # for(i in 1:length(pars)){
          stepsize <- step *10
          colout <- NA
          dolpchecks <- FALSE #set to true to try the log prob checks again...
          while(any(is.na(colout)) && stepsize > 1e-16){
            stepsize <- stepsize * .1
            lpdifok<-FALSE
            lpdifcount <- 0
            lpdifdirection <- 0
            lpdifmultiplier <- 1
            # message('par',i)
            while(!lpdifok & lpdifcount < 15){
              # message(stepsize)
              # message(paste(i,'  col=',colout,'  lpdifmultiplier=',lpdifmultiplier, '  stepsize=',stepsize))
              lpdifok <- TRUE
              lpdifcount <- lpdifcount + 1
              uppars<-pars
              downpars<-pars
              uppars[i]<-pars[i]+stepsize
              downpars[i]<-pars[i]-stepsize
              uplp= try(smf$log_prob(upars=uppars,adjust_transform=TRUE,gradient=TRUE)) #lpg(uppars)
              downlp = try(smf$log_prob(upars=downpars,adjust_transform=TRUE,gradient=TRUE)) #lpg(downpars)
              if(class(uplp)=='try-error' || class(downlp)=='try-error'){
                lpdifok <- TRUE
                upgrad <- rep(NA,length(pars))
                downgrad <- rep(NA,length(pars))
                dolpchecks <- FALSE
                # if(stepsize < 1e-12) browser()
              } else{
                upgrad= attributes(uplp)$gradient
                downgrad = attributes(downlp)$gradient

                if(dolpchecks){
                  if(abs(uplp-downlp) > lpdifmax) {
                    # message(paste0('decreasing step for ', i))
                    lpdifok <- FALSE
                    if(lpdifdirection== 1) {
                      lpdifmultiplier = lpdifmultiplier * .5
                    }
                    stepsize = stepsize * 1e-2 * lpdifmultiplier
                    lpdifdirection <- -1
                  }
                  if(abs(uplp-downlp) < lpdifmin) {
                    # browser()
                    # message(paste0('increasing step for ', i))
                    lpdifok <- FALSE
                    if(lpdifdirection== -1) {
                      lpdifmultiplier = lpdifmultiplier * .5
                    }
                    stepsize = stepsize * 100 * lpdifmultiplier
                    lpdifdirection <- 1
                  }
                  if(any(is.na(c(uplp,downlp)))) stepsize = stepsize * .1
                }
              }
            }
            # hessout[i,]<- (upgrad-downgrad) /stepsize/2
            colout<- (upgrad-downgrad) /stepsize/2
            # print(colout)
          }
          return(rbind(colout))
        },cores=cores)
        return(t(hessout))
      }


      hess1s<-function(pars,direction=1,step=1e-5,lpdifmin=1e-6, lpdifmax=1e-1){
        hessout<-matrix(NA,nrow=length(pars),ncol=length(pars))
        bestlp=lpg(pars)
        basegrad=gradout
        for(i in 1:length(pars)){
          stepsize <- step
          lpdifok<-FALSE
          lpdifcount <- 0
          lpdifdirection <- 0
          lpdifmultiplier <- 1
          while(!lpdifok & lpdifcount < 15){
            lpdifok <- TRUE
            lpdifcount <- lpdifcount + 1
            uppars<-pars
            uppars[i]<-pars[i]+stepsize*direction
            uplp=lpg(uppars)
            upgrad=gradout
            if(abs(uplp-bestlp) > lpdifmax) {
              # message(paste0('decreasing step for ', i))
              lpdifok <- FALSE
              if(lpdifdirection== 1) {
                lpdifmultiplier = lpdifmultiplier * .5
              }
              stepsize = stepsize * 1e-2 * lpdifmultiplier
              lpdifdirection <- -1
            }
            if(abs(uplp-bestlp) < lpdifmin) {
              # message(paste0('increasing step for ', i))
              lpdifok <- FALSE
              if(lpdifdirection== -1) {
                lpdifmultiplier = lpdifmultiplier * .5
              }
              stepsize = stepsize * 100 * lpdifmultiplier
              lpdifdirection <- 1
            }
            hessout[i,]<- (upgrad-basegrad) /stepsize*direction
          }

          # }
        }
        return(t(hessout))
      }


      # A more numerically stable way of calculating log( sum( exp( x ))) Source:
      # http://r.789695.n4.nabble.com/logsumexp-function-in-R-td3310119.html
      log_sum_exp <- function(x) {
        xmax <- which.max(x)
        log1p(sum(exp(x[-xmax] - x[xmax]))) + x[xmax]
      }


      if(is.na(sampleinit[1])){
        # hessup=hess1s(pars = est2,direction = 1,step = 1e-4,lpdifmin = 1e-4,lpdifmax = 1e-3)
        # hessdown=hess1s(pars = est2,direction = -1,step = 1e-4,lpdifmin = 1e-4,lpdifmax = 1e-3)
        # hess=(hessup+hessdown)/2
        hess=grmat(pars=est2,step=1e-6)
        if(any(is.na(hess))) warning(paste0('Hessian could not be computed for pars ', paste0(which(apply(hess,1,function(x) any(is.na(x)))),collapse=', '), ' -- standard errors will be nonsense, model adjustment may be needed.',collapse=''))
        diag(hess)[is.na(diag(hess))]<- -1
        hess[is.na(hess)] <- 0
        hess = (hess/2) + t(hess/2)
        # neghesschol = try(chol(-hess),silent=TRUE)

        mchol=try(t(chol(solve(-hess))),silent=TRUE)
        if(class(mchol)=='try-error') {
          message('Hessian not positive-definite so approximating, treat SE\'s with caution, consider respecification / priors.')
          npd <- TRUE
        } else npd <- FALSE
        # if(class(mchol)=='try-error') {
        mcov=MASS::ginv(-hess) #-optimfit$hessian)
        mcov=as.matrix(Matrix::nearPD(mcov)$mat)
      }

      if(!is.na(sampleinit[1])){
        mcov = cov(sampleinit)*1.5+diag(1e-6,ncol(sampleinit))
        est2 = apply(sampleinit,2,mean)
        bestfit = 9e100
        optimfit <- suppressWarnings(list(par=sampling(sm,standata,iter=2,control=list(max_treedepth=1),chains=1,show_messages = FALSE,refresh=0)@inits[[1]]))
      }

      mcovl <- list()
      mcovl[[1]]=mcov
      delta=list()
      delta[[1]]=est2
      samples <-matrix(NA)
      resamples <- c()
      prop_dens <-c()
      target_dens<-c()
      sample_prob<-c()
      ess <- 0
      qdiag<-0

      if(!is) {
        nresamples = finishsamples
        resamples <- matrix(unlist(lapply(1:nresamples,function(x){
          delta[[1]] + t(chol(mcovl[[1]])) %*% t(matrix(rnorm(length(delta[[1]])),nrow=1))
        } )),byrow=TRUE,ncol=length(delta[[1]]))
        message('Importance sampling not done -- interval estimates via hessian based sampling only')
      }

      if(is){
        targetsamples <- finishsamples * finishmultiply
        # message('Adaptive importance sampling, loop:')
        j <- 0
        while(nrow(samples) < targetsamples){
          j<- j+1
          if(j==1){
            # if(!npd)
            samples <- mvtnorm::rmvt(isloopsize, delta = delta[[j]], sigma = mcovl[[j]],   df = tdf)

          } else {
            delta[[j]]=colMeans(resamples)
            mcovl[[j]] = as.matrix(Matrix::nearPD(cov(resamples))$mat) #+diag(1e-12,ncol(samples))
            samples <- rbind(samples,mvtnorm::rmvt(isloopsize, delta = delta[[j]], sigma = mcovl[[j]],   df = tdf))
            # samples <- rbind(samples, MASS::mvrnorm(isloopsize, delta[[j]],mcovl[[j]]))
          }
          # if(j > 1 || !npd)
          prop_dens <- mvtnorm::dmvt(tail(samples,isloopsize), delta[[j]], mcovl[[j]], df = tdf,log = TRUE)
          # prop_dens <- mvtnorm::dmvnorm(tail(samples,isloopsize), delta[[j]], mcovl[[j]],log = TRUE)
          # prop_dens <- ctdmvnorm(tail(samples,isloopsize), delta[[j]], mcovl[[j]])

          parallel::clusterExport(cl, c('samples'),environment())

          target_dens[[j]] <- c(flexsapply(cl, parallel::clusterSplit(cl,1:isloopsize), function(x){
            # eval(parse(text=paste0('library(rstan)')))

            smf <- stan_reinitsf(sm,standata, fast=TRUE)

            lp<-function(parm) {
              out<-try(smf$log_prob(upars=parm, adjust_transform = TRUE, gradient=FALSE),silent = TRUE)
              if(class(out)=='try-error') {
                out=-1e60
              }
              return(out)
            }
            out <- apply(tail(samples,isloopsize)[x,],1,lp)

            try(dyn.unload(file.path(tempdir(), paste0(smf@stanmodel@dso@dso_filename, .Platform$dynlib.ext))),silent = TRUE)
            return(out)

          },cores=cores))
          target_dens[[j]][is.na(target_dens[[j]])] <- -1e60
          if(all(target_dens[[j]] < -1e100)) stop('Could not sample from optimum! Try reparamaterizing?')
          if(any(target_dens[[j]] > bestfit)){
            oldfit <- bestfit
            try2 <- TRUE
            bestfit<-max(target_dens[[j]],na.rm=TRUE)
            betterfit<-TRUE
            init = samples[which(unlist(target_dens) == bestfit),]
            message('Improved fit found - ', bestfit,' vs ', oldfit,' - restarting optimization')
            break
          }

          if(nrow(samples) >= targetsamples) nresamples <- finishsamples else nresamples = min(5000,length(samples)/5)


          target_dens2 <- target_dens[[j]] -max(target_dens[[j]],na.rm=TRUE) + max(prop_dens) #adjustment to get in decent range, doesnt change to prob
          target_dens2[!is.finite(target_dens[[j]])] <- -1e30
          weighted_dens <- target_dens2 - prop_dens
          # psis_dens <- psis(matrix(target_dens2,ncol=length(target_dens2)),r_eff=NA)
          # sample_prob <- weights(psis_dens,normalize = TRUE,log=FALSE)
          # plot(target_dens2,prop_dens)

          newsampleprob <- exp((weighted_dens - log_sum_exp(weighted_dens)))
          counter <- 1
          while(counter < 30 &&
              length(unique(sample(x=1:length(newsampleprob),size=length(newsampleprob),prob=newsampleprob,replace=TRUE))
                ) < 50) {
           if(counter==1) message ('Sampling problematic -- trying to recover... ')
            counter = counter + 1
            if(counter == 30) stop('Importance sampling failed -- either posterior mode not found, or mode is inadequate starting point for sampling')
            weighted_dens <- weighted_dens /2
            newsampleprob <- exp((weighted_dens - log_sum_exp(weighted_dens)))
            # plot(newsampleprob)
          }
          sample_prob <- c(sample_prob,newsampleprob) #sum to 1 for each iteration, normalise later

          sample_prob[!is.finite(sample_prob)] <- 0
          sample_prob[is.na(sample_prob)] <- 0
          # points(target_dens2[sample_prob> (1/isloopsize * 10)], prop_dens[sample_prob> (1/isloopsize * 10)],col='red')
          resample_i <- sample(1:nrow(samples), size = nresamples, replace = ifelse(nrow(samples) > targetsamples,FALSE,TRUE),
            prob = sample_prob / sum(sample_prob))
          # resample_i <- sample(tail(1:nrow(samples),isloopsize), size = nresamples, replace = ifelse(j == isloops+1,FALSE,TRUE),
          #   prob = tail(sample_prob,isloopsize) / sum(tail(sample_prob,isloopsize) ))

            message(paste0('Importance sample loop ',j,', ',length(unique(resample_i)), ' unique samples, from ', nresamples,' resamples of ', nrow(samples),' actual, prob sd = ', round(sd(sample_prob),4),
              ', max chance = ',max(sample_prob) * isloopsize))
            if(length(unique(resample_i)) < 100) {
              message('Sampling ineffective, unique samples < 100 -- try increasing samples per step (isloopsize), or use HMC (non optimizing) approach.')
              # return(est)
          }

          resamples <- samples[resample_i, , drop = FALSE]

          #check max resample probability and drop earlier samples if too high
          dropuntil <- ceiling(max(c(0,which(sample_prob > (chancethreshold / isloopsize)) / isloopsize),na.rm=TRUE))*isloopsize
          if((isloopsize - dropuntil) > isloopsize) dropuntil <- dropuntil -isloopsize

          if(dropuntil > 0){
            resamples <- resamples[-(0:dropuntil),,drop=FALSE]
            sample_prob <- sample_prob[-(0:dropuntil)]
            samples <- samples[-(0:dropuntil),,drop=FALSE]
          }

          ess[j] <- (sum(sample_prob[resample_i]))^2 / sum(sample_prob[resample_i]^2)
          qdiag[j]<-mean(unlist(lapply(sample(x = 1:length(sample_prob),size = 500,replace = TRUE),function(i){
            (max(sample_prob[resample_i][1:i])) / (sum(sample_prob[resample_i][1:i]) )
          })))

        }
      }
    }
  }#end while no better fit



  if(!estonly){
    if(!is) lpsamples <- NA else lpsamples <- unlist(target_dens)[resample_i]

    transformedpars=stan_constrainsamples(sm = sm,standata = standata,samples=resamples,cores=cores)

    # quantile(sapply(transformedpars, function(x) x$rawpopcorr[3,2]),probs=c(.025,.5,.975))
    # quantile(sapply(transformedpars, function(x) x$DRIFT[1,2,2]),probs=c(.025,.5,.975))

    sds=try(suppressWarnings(sqrt(diag(mcov))))  #try(sqrt(diag(solve(optimfit$hessian))))
    if(class(sds)=='try-error') sds <- rep(NA,length(est2))
    lest= est2 - 1.96 * sds
    uest= est2 + 1.96 * sds

    transformedpars_old=NA
    try(transformedpars_old<-cbind(unlist(constrain_pars(smfull, upars=lest)),
      unlist(constrain_pars(smfull, upars= est2)),
      unlist(constrain_pars(smfull, upars= uest))),silent=TRUE)
    try(colnames(transformedpars_old)<-c('2.5%','mean','97.5%'),silent=TRUE)

    stanfit=list(optimfit=optimfit,stanfit=smfull, rawest=est2, rawposterior = resamples, transformedpars=transformedpars,transformedpars_old=transformedpars_old,
      isdiags=list(cov=mcovl,means=delta,ess=ess,qdiag=qdiag,lpsamples=lpsamples ))
  }
  if(estonly) stanfit=list(optimfit=optimfit,stanfit=smf, rawest=est2)
  suppressWarnings(do.call(par,parbase)) #reset par in case plots done
  return(stanfit)
}

