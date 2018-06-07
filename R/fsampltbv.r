#fSampLTBv <- function(my,amX,mPhi,vmu,mSig,amG2inv,vd,mbo){
#  .Call('fSampLTBv',package='ltvar',my,amX,mPhi,vmu,mSig,amG2inv,vd,mbo)
#
#  ns <- nrow(my)
#  nk <- ncol(my)
#  np <- ncol(mPhi)
#  mba <- mbo
#
#  mSigi <- solve(mSig)
#  mSigp <- mSigi %*% ( diag(1,np) + mPhi^2 )
#  vpm1  <- ( diag(1,np) - mPhi ) %*% vmu
#  vpm2  <- ( diag(1,np) -2 * mPhi + mPhi^2) %*% vmu
#
#  mSig0i <- (diag(1,np)-mPhi^2) %*% mSigi
#
#  for(ii in 1:ns){
#
#    if(ii==1){
#
#      mSighi <- amX[ii,,] %*% amG2inv[ii,,] %*% t(amX[ii,,]) + mSig0i + mSigi %*% mPhi^2
#      mSigh  <- ginv(mSighi)
#      vbh    <- mSigh %*% ( amX[ii,,] %*% amG2inv[ii,,] %*% my[ii,] + mSig0i %*% vmu + mSigi %*% mPhi %*% ( mba[ii+1,] - vpm1 ))
#
#    }
#    else if(ii < ns){
#
#      mSighi <- amX[ii,,] %*% amG2inv[ii,,] %*% t(amX[ii,,]) + mSigp
#      mSigh  <- ginv(mSighi)
#      vbh    <- mSigh %*% ( amX[ii,,] %*% amG2inv[ii,,] %*% my[ii,] + mSigi %*% (mPhi %*% (vba1 + mba[ii+1,]) + vpm2))
#
#    }
#    else if(ii == ns){
#
#      mSighi <- amX[ii,,] %*% amG2inv[ii,,] %*% t(amX[ii,,]) + mSigi
#      mSigh  <- ginv(mSighi)
#      vbh    <- amX[ii,,] %*% amG2inv[ii,,] %*% my[ii,] + mSigi %*% ( mPhi %*% vba1 + vpm1)
#
#    }
#
#    vba1 <- vbh
#    vbo  <- mba[ii,]
#    vbn  <- mvrnorm(1,mu=vbh,Sigma=mSigh)
#
#    dhn  <- -0.5 * log(det(mSigh) + 1e-6)
#    dho  <-  dhn
#
#    dhn  <- dhn - 0.5 * t(vbn - vbh) %*% mSighi %*% (vbn - vbh)
#    dho  <- dho - 0.5 * t(vbo - vbh) %*% mSighi %*% (vbo - vbh)
#
#    if(sum(abs(vbn)<vd)>0){
#
#      mxh <- amX[ii,,]
#
#      for(jj in 1:nk){
#
#        mxh[abs(vbn)<vd,jj] <- 0
#
#      }
#
#      if(ii == 1){
#
#        mSighi <- mxh %*% amG2inv[ii,,] %*% t(mxh) + mSig0i + mSigi %*% mPhi^2
#        mSigh  <- ginv(mSighi)
#        vbh    <- mSigh %*% ( mxh %*% amG2inv[ii,,] %*% my[ii,] + mSig0i %*% vmu + mSigi %*% mPhi %*% ( mba[ii+1,] - vpm1 ))
#
#      }
#      else if(ii < ns){
#
#        mSighi <- mxh %*% amG2inv[ii,,] %*% t(mxh) + mSigp
#        mSigh  <- ginv(mSighi)
#        vbh    <- mSigh %*% ( mxh %*% amG2inv[ii,,] %*% my[ii,] + mSigi %*% ( mPhi %*% ( vba1 + mba[ii+1,] ) + vpm2 ))
#
#      }
#      else if(ii == ns){
#
#        mSighi <- mxh %*% amG2inv[ii,,] %*% t(mxh) + mSigi
#        mSigh  <- ginv(mSighi)
#        vbh    <- mSigh %*% ( mxh %*% amG2inv[ii,,] %*% my[ii,] + mSigi %*% ( mPhi %*% vba1 + vpm1 ))
#
#      }
#
#      dln <- -0.5 * log( det(mSigh) + 1e-6) - 0.5 * t(vbn - vbh) %*% ginv(mSigh) %*% ( vbn - vbh )
#
#    }
#    else{
#      dln <- dhn
#    }
#
#    if(sum(abs(vbo)<vd)>0){
#
#      s1 <- as.vector(abs(vbo)>vd)
#      mxh <- amX[ii,,]*s1
#
#      if(ii == 1){
#
#        mSighi <- mxh %*% amG2inv[ii,,] %*% t(mxh) + mSig0i + mSigi %*% mPhi^2
#        mSigh  <- ginv(mSighi)
#        vbh    <- mSigh %*% (mxh %*% amG2inv[ii,,] %*% my[ii,] + mSig0i %*% vmu + mSigi %*% mPhi %*% ( mba[ii+1,] - vpm1 ))
#
#      }
#      else if(ii < ns){
#
#        mSighi <- mxh %*% amG2inv[ii,,] %*% t(mxh) + mSigp
#        mSigh  <- ginv(mSighi)
#        vbh    <- mSigh %*% ( mxh %*% amG2inv[ii,,] %*% my[ii,] + mSigi %*% ( mPhi %*% ( vba1 + mba[ii+1,]) + vpm2 ))
#
#      }
#      else if(ii == ns){
#
#        mSighi <- mxh %*% amG2inv[ii,,] %*% t(mxh) + mSigi
#        mSigh  <- ginv(mSighi)
#        vbh    <- mSigh %*% ( mxh %*% amG2inv[ii,,] %*% my[ii,] + mSigi %*% ( mPhi %*% vba1 + vpm1 ))
#
#      }
#
#      dlo <- -0.5 * log( det(mSigh) + 1e-6) - 0.5 * t(vbo-vbh) %*% ginv(mSigh) %*% (vbo-vbh)
#
#    }
#    else{
#
#      dlo <- dho
#
#    }
#
#    dfrac <- min(1,exp(dln-dhn-dlo+dho))
#    if(is.na(dfrac)){ dfrac = 0}
#
#    if(runif(1) < dfrac){
#
#      mba[ii,] <- vbn
#      vba1     <- vbn
#
#    }
#    else{
#
#      vba1 <- mba[ii,]
#
#    }
#
#  }
#
#  return(list(mba=mba))
#
#}
