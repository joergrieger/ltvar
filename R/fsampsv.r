fSVSamp <- function(vy,vh,dphi,dsig2,dh00,dsig02,nk){

  ns   <- length(vy)
  nite <- 10

  vhs <- vh

  repeat{

    v1 <- t(floor(ns * ( seq(1:nk) + runif(nk) ) / (nk+2)))
    vk <- as.vector(cbind(1, v1, ns + 1))
    dvk <- diff(vk)

    if(!sum(dvk<2)){ break }

  }

  for(jj in 1:(nk+1)){

    ir  <- vk[jj]
    id  <- vk[jj+1] - vk[jj]
    ird <- vk[jj+1] - 1

    vyi <- vy[ir:ird]
    vho <- vh[ir:ird]
    vhn <- array(0,dim=c(id,1))

    if(jj < nk ){

      dhrd1 <- vh[ird+1]

    }

    for(ii in 0:(nite+1)){

      if(ii==1){

        vhh <-  vhn

      }
      else{

        vhh <- vho

      }

      vgder2 <- -0.5 * vyi^2 / exp(vhh)
      vgder1 <- -0.5 - vgder2
      vsiga2 <- -1 / vgder2
      vha <- vhh + vsiga2 * vgder1

      if( jj <= nk ){

        vsiga2[id] <- 1 / ( -vgder2[id-1] + dphi^2 / dsig2 )
        vha[id]    <- vsiga2[id] * ( vgder1[id] - vgder2[id] * vhh[id] + dphi*dhrd1/dsig2)

      }



      if( jj == 0){

        da  <- dphi * vhs[ir-1]
        dho <- dphi * vhs[ir-1]

      }
      else{

        da  <- dh00
        dh0 <- dh00

      }

      dH2 <- dsig2

      if( jj == 0){

        dH20 <- dsig2
        dP   <- dH20

      }
      else{

        dH20 <- dsig02
        dP   <- dH20

      }



      ve    <- array(0,dim=c(id,1))
      vDinv <- array(0,dim=c(id,1))
      vK    <- array(0,dim=c(id,1))
      vL    <- array(0,dim=c(id,1))
      vu    <- array(0,dim=c(id,1))

      for(tt in 1:id){

        ve[tt]    <- vha[tt] - da
        vDinv[tt] <- 1 / (dP + vsiga2[tt] )
        vK[tt]    <- dphi * dP * vDinv[tt]
        vL[tt]    <- dphi - vK[tt]

        da <- dphi * da + vK[tt] * ve[tt]
        dP <- dphi * dP + vL[tt] + dH2

      }

      if( ii < nite ){

        dr <- 0
        dU <- 0

        for(tt in id:1){

          dC <- dH2 * ( 1 - dU * dH2 )
          deps <- 0
          vu[tt] <- dH2 * dr + deps
          dV <- dH2 * dU * vL[tt]

          dCinv <- 1 / dC
          dr <- vDinv[tt] * ve[tt] + vL[tt] * dr - dV * dCinv * deps
          dU <- vDinv[tt] + vL[tt]^2 * dU + dV^2 * dCinv

        }

        dC  <- max(0, dH20 * ( 1 - dU * dH20 ))
        du0 <- dH20 * dr + sqrt(dC) * rnorm(1,1)
        vhn[1] <- dh0 + du0

        for(tt in 2:id){

          vhn[tt] <- dphi * vhn[tt-1] + vu[tt-1]

        }

      }
      else{

        docounter <- 1
        repeat{

          dr <- 0
          dU <- 0
          for(tt in id:1){

            dC   <- max(0, dH2 * ( 1 - dU * dH2 ))
            deps <- sqrt(dC) * rnorm(1)
            vu[tt] <- dH2 * dr + deps
            dV <- dH2 * dU * vL[tt]

            dCinv <- 1/dC
            dr <- vDinv[tt] * ve[tt] + vL[tt] * dr - dV * dCinv * deps
            dU <- vDinv[tt] + vL[tt]^2 * dU + dV^2 * dCinv

          }

          dC  <- max(0,dH20 * ( 1 - dU * dH20 ))
          du0 <- dH20 * dr + sqrt(dC) * rnorm(1)

          vhn[1] <- dh0 + du0
          for(tt in 1:(id-1)){
            vhn[tt+1] <- dphi * vhn[tt] + vu[tt]
          }

          dpron <- sum( -0.5 * ( vhh + vyi^2 / exp(vhh) ) + vgder1 * (vhn - vhh) + 0.5 * vgder2 * ( vhn - vhh )^2 )
          dposn <- sum( - 0.5 * ( vhn + vyi^2 / exp(vhn) ))
          docounter <- docounter + 1

          dfrac <- exp(dposn - dpron)
          if(is.na(dfrac)){ dfrac <- 0}

          if( docounter > 100 ){ break }
          if( runif(1) < dfrac ){ break }
        }

        if( docounter < 100 ){

          dproo <- sum(-0.5*( vhh + vyi^2 / exp(vhh) ) + vgder1 * ( vho - vhh ) + 0.5 * vgder2 * ( vho - vhh )^2)
          dposo <- sum(-0.5 * ( vho + vyi^2 / exp(vho) ))

          dfrac <- exp( dposn + min(dposo,dproo) - dposo - min(dposn,dpron) )
          if(is.na(dfrac)){ dfrac <- 0 }

          if( runif(1) < dfrac ){

            vhs[ir:ird] <- vhn

          }

        }

      }

    }

  }

  return(vhs)

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
