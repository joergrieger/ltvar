fSampLT <- function(my,amX,mp,amG2inv,vd,vdc){
  ns <- nrow(mp)
  nk <- ncol(mp)
  vdo <- vd

  for(ii in 1:nk){
    vdn <- vdo
    vdn[ii] <- runif(1) * vdc[ii]
    mdn <- mp
    mdo <- mp
    
    for(jj in 1:ns){
      
      if(abs(mp[jj,ii]) < vdn[ii] ){ mdn[jj,ii] <- 0 }
      if(abs(mp[jj,ii]) < vdo[ii] ){ mdo[jj,ii] <- 0 }
      
    }
    
    dlo <- 0
    dln <- 0
    
    for(jj in 1:ns){
      
      ve  <- my[jj,] - mdn[jj,] %*% amX[jj,,]
      dln <- dln + ve %*% amG2inv[jj,,] %*% t(ve)
      
      ve  <- my[jj,] - mdo[jj,] %*% amX[jj,,]
      dlo <- dlo + ve %*% amG2inv[jj,,] %*% t(ve)
      
    }

    if(runif(1)<exp(-0.5*(dln-dlo))){
      
      vdo <- vdn
      
    }
    
  }
  return(vdo)
}