#' @export
#' @title forecasting using an ltvar model
#' @param obj a fitted ltvar model
#' @param n.ahead forecast horizon
#' @param ... currently not used

forecast <- function(obj,n.ahead,...) UseMethod("forecast")

#' @export
#' @title forecasting using an ltvar model
#' @param obj a fitted ltvar model
#' @param n.ahead forecast horizon
#' @param ... currently not used

forecast.ltvar <- function(obj,n.ahead,...){

  # Preliminaries
  y = obj$data_info$data
  p = obj$model_info$p
  intercept = obj$model_info$intercept

  nreps <- dim(obj$mcmc_draws$mb)[1]
  nTotal <- dim(obj$mcmc_draws$mb)[2]

  nb <- dim(obj$mcmc_draws$mb)[3]
  na <- dim(obj$mcmc_draws$ma)[3]
  nk <- dim(obj$mcmc_draws$mh)[3]

  forecast_storage <- array(0,dim=c(n.ahead,nk,nreps))


  for(ii in 1:nreps){

    mb    <- obj$mcmc_draws$mb[ii,nTotal,]
    vmub  <- obj$mcmc_draws$vmub[ii,,]
    mSigb <- obj$mcmc_draws$mSigb[ii,,]
    mPhib <- obj$mcmc_draws$mPhib[ii,,]
    vdb   <- obj$mcmc_draws$vdb[ii,,]

    ma    <- obj$mcmc_draws$ma[ii,nTotal,]
    vmua  <- obj$mcmc_draws$vmua[ii,,]
    mSiga <- obj$mcmc_draws$mSiga[ii,,]
    mPhia <- obj$mcmc_draws$mPhia[ii,,]
    vda   <- obj$mcmc_draws$vda[ii,,]

    mPhih <- obj$mcmc_draws$mPhih[ii,,]
    mSigh <- obj$mcmc_draws$mSigh[ii,,]

    vbf  <- obj$mcmc_draws$mb[ii, nTotal - 1,]
    vaf  <- obj$mcmc_draws$ma[ii, nTotal - 1,]
    vhf  <- obj$mcmc_draws$mh[ii, nTotal - 1,]

    vgam <- obj$mcmc_draws$vgam[ii,,]

    myf <- y[nTotal:(nTotal - p + 1),]

    for(jj in 1:(n.ahead)){

      # Calculate time-varying parameters

      vbf  <- vmub + mPhib %*% (vbf - vmub) + t(chol(mSigb)) %*% rnorm(nb)
      vbfd <- vbf
      vbfd[abs(vbf) < vdb] <- 0

      vaf  <- vmua + mPhia %*% (vaf - vmua) + t(chol(mSiga)) %*% rnorm(na)
      vafd <- vaf
      vafd[abs(vaf) < vda] <- 0
      vhf  <- mPhih %*% vhf + t(chol(mSigh)) %*% rnorm(nk)

      mOmsq <- solve(fAt(vafd,nk)) %*% diag(exp(as.vector(vhf)/2)) * sqrt(vgam)

      vb <- matrix(vbfd,ncol = nk)
      if(intercept){

        tmp <-c(1,myf)
      }
      else{

        tmp <- myf

      }
      vyf <- t(tmp %*% vb) + t(chol(mOmsq)) %*% rnorm(nk)

      # Geet new variables
      if(p > 1){

        myf[1:(p-1),] <- myf[2:p,]
        myf[p,] <- vyf

      }
      else{

        myf <- vyf

      }
      # Store forecasts
      forecast_storage[jj,,ii] <- vyf

    }
  }
  return(forecast_storage)
}
