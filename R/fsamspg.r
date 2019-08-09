#' @export
#' @title samples sigma, phi and gamma for stochastic volatility sampler
#' @param my current draw of stochastic volatility
#' @param mh previous draw of stochastic volatility
#' @param mPhi Phi-matrix from previous draw
#' @param vgam gamma draw from previous round
#' @param dnu
#' @param dg0,dG0 prior parameters on gamma, inverse gamma distribution
#' @param da0,db0 prior parameters on phi, shape parameters for beta distribution
#' @param dnu,dV0 prior parameters on sigma, inverse gamma distribution

fSampSPG <- function(my,mh,mPhi,vgam,dnu,dV0,da0,db0,dg0,dG0){

  # Preliminaries

  ns <- nrow(my)
  nk <- ncol(my)


  mSign <- diag(1,nk)
  mPhin <- mPhi
  vgamn <- vgam

  # Sampling Sigma
  sse <- 0
  for(jj in 2:ns){

    sse <- sse + ( mh[jj,] - mh[jj-1,] %*% mPhi )^2

  }

  vV <- dV0 + mh[1,]^2 %*% ( diag(1,nk) - mPhi^2 ) + sse

  for(jj in 1:nk){

    mSign[jj,jj] <- 1/rgamma(1,dnu/2,vV[jj]/2)

  }

  # sampling gamma
  sums <- array(0,dim=c(nk))

  for(jj in 1:nk){

    sums[jj] <- sum(my[,jj,jj]^2/exp(mh[,jj]))

  }

  vG <- dG0 + sums

  for(jj in 1:nk){

    vgamn[jj] <- 1/rgamma(1,(dg0+ns),vG[jj]/2)

  }

  # Sampling Phi
  sumv <- 0
  sumv2 <- 0
  for(jj in 2:ns){
    sumv <- sumv + ( mh[jj,] )^2
    sumv2 <- sumv2 + ( mh[jj,] ) * ( mh[jj-1,] )
  }

  vphii <- sumv2 / sumv
  vsigi <- diag(mSign) / sumv

  for(jj in 1:nk){
    docounter <- 0
    repeat{

      dphin <- vphii[jj] + rnorm(1) * sqrt(vsigi[jj])

      docounter <- docounter + 1
      if( docounter > 100 ){ break }
      if( dphin < 1){ break }

    }

    if( (docounter < 100) & (dphin < 1) ){

      dphio <- mPhi[jj,jj]
      dfrac <- ( dbeta((dphin+1)/2,da0,db0) / dbeta((dphio+1)/2,da0,db0) ) * ( sqrt( 1 - dphin^2 ) / sqrt( 1 - dphio^2 ) )
      if(is.na(dfrac)){ dfrac <- 0 }

      if( runif(1) < dfrac ){
        mPhin[jj,jj] <- dphin
      }

    }
  }



  # return results

  return( list( mPhi = mPhin,vgam = vgamn, mSig = mSign ) )

}
