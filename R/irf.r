#' @export
#' @title impulse-responses for latent threshold VARs
#' @param obj A fitted ltvar model
#' @param n.ahead Integer specifying the steps
impulse_response <- function(obj,n.ahead) UseMethod("impulse_response")

#' Impulse-responses for latent threshold VARs.
#'
#' This function computes impulse-response functions for latent-threshold VARS. The user supplies an estimated latent-threshold model and the number of steps for the IRFs. The function returns an object of the class ltirf with the computed mean and the 5% and 95% quantil for the IRFs.
#'
#' @param obj A fitted ltvar model
#' @param n.ahead Integer specifying the steps
#' @importFrom stats quantile
#' @export
impulse_response.ltvar <- function(obj,n.ahead){

  dims  <- dim(obj$mcmc_draws$mb)
  nreps <- dims[1]
  nT    <- dims[2]
  dims  <- dim(obj$mcmc_draws$vda)
  nvar  <- dims[2]

  nolags <- obj$model_info$p
  intercept <- obj$model_info$intercept

  irf_storage <- array(0,dim=c(nvar,nvar,n.ahead,nT,3))

  # Outer loop is over all periods, while inner loop is over all draws.
  # This way we can also get quantiles of the impulse-response function

  for(ii in 1:nT){ # Loop over all periods

    # temporary storage of impulse-response functions
    tmp_irf <- array(0,dim=c(nvar,nvar,n.ahead,nreps))

    for(jj in 1:nreps){ # loop over all draws

      # Prepare draw for beta

      Alpha <- obj$mcmc_draws$mb[jj,ii,]
      threshold_beta <- obj$mcmc_draws$vdb[jj,,1]
      Alpha[abs(Alpha) < threshold_beta] <- 0
      Alpha <- matrix(Alpha,ncol=nvar)

      if(intercept == TRUE){

        Alpha <- Alpha[-c(1),]

      }

      # Prepare draw for variance-covariance matrix
      vda <- obj$mcmc_draw$ma[jj,ii,]
      ma  <- fAt(vda,nvar)
      mh  <- diag(obj$mcmc_draw$mh[jj,ii,])
      Sigma <- solve(ma) %*% (mh %*% mh) %*% solve(t(ma))

      # Get Impulse-Response
      tmp_irf[,,,jj] <- compute_impulse_responses(Alpha = Alpha, Sigma = Sigma, nolags = nolags, nhor = n.ahead)

    }
    for(kk in 1:nvar){
      for(ll in 1:nvar){
        for(mm in 1:n.ahead){

          irf_storage[kk,ll,mm,ii,1] <- mean(tmp_irf[kk,ll,mm,])
          irf_storage[kk,ll,mm,ii,2] <- stats::quantile(tmp_irf[kk,ll,mm,],probs=c(0.05))
          irf_storage[kk,ll,mm,ii,3] <- stats::quantile(tmp_irf[kk,ll,mm,],probs=c(0.95))

        }
      }
    }

  }
  retlist = structure(list(irf = irf_storage),class="ltirf")
  return(retlist)

}

#' Function to draw one single path for the impulse-response functions.
#' @param Alpha the estimate of the VAR-coefficients
#' @param Sigma the estimate of the Variance-covariance matrix
#' @param nolags Number of lags in the model
#' @param nhor the horizon of the impulse-response functions
compute_impulse_responses <- function(Alpha, Sigma, nolags,nhor){

  K <- nrow(Sigma)
  bigj <- matrix(0,K,K*nolags)
  bigj[1:K,1:K] <- diag(K)

  PHI_mat <- companion(Alpha,nolags)
  biga <- PHI_mat
  bigai <- biga

  #Prepare storage variales for Impulse-Responses
  shock <- t(chol(Sigma))
  impresp <- matrix(0,K,K*nhor)
  impresp[1:K,1:K] <- shock

  for(ii in 1:(nhor - 1)){

    impresp[,(ii*K+1):((ii + 1) * K)] <- bigj %*% bigai %*% t(bigj) %*% shock
    bigai <- bigai %*% biga

  }

  # Reshapre impulse-response matrix
  imp <- array(0,dim=c(K,K,nhor))
  for(ii in 1:K){
    for(ij in 1:nhor){

      jj <- ii + K  * (ij - 1 )
      imp[ii,,ij] <- impresp[,jj]
    }

  }

  return(imp)
}


companion <- function(betadraw,nolags){

  K <- ncol(betadraw)

  comp <- array(0, dim=c( K * nolags, K * nolags ))

  Ai <- array(0,dim=c( K, K, nolags ))

  if(nolags > 1){

    for(jj in 1:nolags){

      indxmin  <- ( jj - 1 ) * K + 1
      indxmax  <- ( jj - 1 ) * K + K

      Ai[,,jj] <- betadraw[indxmin:indxmax,]

      comp[1:K,indxmin:indxmax] <- betadraw[indxmin:indxmax,]
    }

    indxmin <- K + 1
    indxmax <- K * nolags
    indxrep <- indxmax - indxmin
    for(jj in 0:indxrep){

      comp[indxmin + jj, jj + 1] <- 1

    }
  }
  else{

    comp <- betadraw

  }

  return(comp)

}
