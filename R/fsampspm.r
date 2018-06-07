#fSampSPM <- function(mp,mSig,mPhi,vmu,dnu,dV0,da0,db0,dm0,ds0,vd,dk0){
#
#  ns <- nrow(mp)
#  np <- ncol(mp)
#
#  mSign <- diag(1,np)
#  mPhin <- mPhi
#  vmun  <- vmu
#
#  # Sampling Sigma
#
#  sse <- 0
#
#  for(jj in 1:(ns-1)){
#
#    sse <- sse + (( mp[jj+1,] - mp[jj,] %*% mPhi - t(vmu) %*% (diag(1,np)-mPhi))^2 )^2
#
#  }
#
#  vV <- dV0 + ( mp[1,] - t(vmu))^2 %*% ( diag(1,np) - mPhi^2 ) + sse
#
#  for(ii in 1:np){
#
#    docounter <- 0
#
#    repeat{
#
#      dsign <- 1 / rgamma(1,dnu/2,vV[ii]/2)
#      dup <- fk0( vmu[ii], dsign, mPhi[ii,ii], dk0)
#
#      if( docounter > 100 ){ break }
#      if( vd[ii] < dup ){ break }
#
#      docounter <- docounter + 1
#
#    }
#
#    if( docounter < 100 ){
#
#      dup2 <- fk0(vmu[ii], mSig[ii,ii], mPhi[ii,ii], dk0)
#      dfrac <- dup2 / dup
#
#      if(runif(1)<dfrac){
#
#        mSign[ii,ii] <- dsign
#
#      }
#
#    }
#
#  }
#
#  # Sammpling Phi
#
#  vsum <- array(0,dim=c(np))
#  vphii <- array(0,dim=c(np))
#
#  for(jj in 1:np){
#
#    vsum[jj] <- sum( ( mp[1:(ns-1),jj] - vmu[jj] )^2 )
#    vphii[jj] <- sum( ( mp[2:ns,jj] - vmu[jj] ) * (mp[1:(ns-1),jj] - vmu[jj]) ) / vsum[jj]
#
#  }
#
#  vsigi <- diag(mSign)/vsum
#
#  for(ii in 1:np){
#
#    doCounter <- 0
#
#    repeat{
#
#      dphin <- vphii[ii] + rnorm(1) * sqrt(vsigi[ii])
#
#      if( abs(dphin) < 1 ){ break }
#      if( doCounter > 100 ){ break }
#
#      doCounter <- doCounter + 1
#
#    }
#    if(doCounter < 100){
#
#      dup <- fk0(vmu[ii],mSign[ii,ii],dphin,dk0)
#      dphio <- mPhi[ii,ii]
#      x1 <- dbeta((dphin+1)/2,da0,db0)
#      x2 <- dbeta((dphio+1)/2,da0,db0)
#
#      dfrac <- x1 / x2 * (( sqrt(1 - dphin^2 )/sqrt( 1 - dphio^2 )) * fk0(vmu[ii],mSign[ii,ii],dphio,dk0) )
#
#      if(runif(1) < dfrac){
#        mPhin[ii,ii] <- dphin
#      }
#
#    }
#  }
#
#  mSig0 <- sqrt(mSign%*%(diag(1,np)-mPhin^2))
#
#  # sampling mu
#
#  x1 <- 1 / ds0^2 + ( 1 - diag( mPhin^2 ) + ( ns - 1 ) * (1-diag(mPhin))^2)/diag(mSign)
#  vsigi <- 1/x1
#
#  sumc <- 0
#
#  for(jj in 1:(ns-1)){
#
#    sumc <- sumc + mp[jj+1,] - t(mPhin %*% mp[jj,])
#
#  }
#  vmui <- vsigi * ( dm0 / ds0^2 + ( mp[1,] * ( 1 - diag(mPhin)^2 ) +( 1 - diag(mPhin)) * sumc ) / diag(mSign) )
#
#  for(jj in 1:np){
#
#    docounter <- 0
#
#    repeat{
#
#      dmun <- abs(vmu[ii]) + rnorm(1) * sqrt(vsigi[ii])
#      dup  <- dmun+dk0 * mSig0[ii,ii]
#
#      if(docounter > 100){ break }
#      if(vd[ii] < dup){ break }
#      docounter <- docounter + 1
#
#    }
#
#    if(docounter < 100){
#
#      dup2 <- (abs(vmu[ii])+dk0*mSig0[ii,ii])
#      dfrac <- dup2/dup
#
#      if(runif(1)<dfrac){
#        vmun[ii] <- dmun
#      }
#    }
#
#  }
#  return(list(mSig=mSign,mPhi=mPhin,vmu=vmun))
#}

fk0 <- function(dmu,dsig2,dphi,dk0){

  return( abs(dmu) + dk0 * sqrt( dsig2 / (1-dphi^2) ))

}

fxh <- function(vyh,nk,na){
  mxh <- array(0,dim=c(nk,na))
  for(ii in 1:(nk-1)){
    ind1 <- ii*(ii-1)/2+1
    ind2 <- (ii+1)*ii/2
    mxh[ii+1,ind1:ind2] <- -vyh[1:ii]
  }
  return(mxh)
}
