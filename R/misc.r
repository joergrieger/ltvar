
#
#
# Function to create lagged data
#
#
#' @importFrom stats embed
lagdata <-function(y,lags,intercept=FALSE){

  T<-nrow(y)
  K<-ncol(y)

  obs <- T-lags

  x  <- stats::embed(y,dimension=lags+1)[,-(1:K)]

  if(intercept==TRUE){

    x<-cbind(1,x)

  }
  yi <- y[(lags+1):T,]

  return(list(y=yi,x=x,obs=obs,T=T,K=K));

}

fAt <- function(va,nk){

  mAt <- diag(1,nk)

  for(ii in 2:nk){

    ind1 <- (ii - 1) * (ii - 2)/2 + 1
    ind2 <- ii * (ii-1) / 2

    mAt[ii,1:(ii-1)] <- va[ind1:ind2]

  }
  return(mAt)
}

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
