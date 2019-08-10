##
##  vhs = fsvsamp(vy, vh,dphi, dsig2, dh00, dsig02, nK)
##
##  "svsamp" implements multi-move sampler for SV model
##  by Shephard & Pitt (1997) and Omori & Watanabe (2004)
##
##  [model]
##    y_t = exp(h_t/2)*eps_t,  eps_t ~ N(0, 1)
##    h_{t+1} = \phi*h_t + ets_t,  eta_t ~ N(0, sig^2)
##    h_0 = h00,  eta_0 ~ N(0, sig0^2)
##
##  [input]
##      vy:     response (ns*1 vector)
##      vh:     current point of h (ns*1 vector)
##      dsig2, dsig02:
##              parameter (scalar)
##      nK:     # of blocks for multi-move sampler
##
##  [output]
##      vhs:  sampled stochastic volatility (ns*1 vector)
##

#' @export
#' @title implements multi-move sampler for SV model
#' @param vy response (ns*1 vector)
#' @param vh current point of h (ns*1 vector)
#' @param nK number of blocks for multi-move sampler
#' @param dphi AR(1) coefficient for transition of volatility
#' @param dsig2,dsig02 parameters (scalar)
#' @param dh00 mean of cofficients

fSVSamp <- function(vy, vh, dphi, dsig2, dh00, dsig02, nK=NULL){



  ##--- set variables ---##
  #dphi = 1

  ns <- length(vy)
  nite <- 5# # of iteration for h_hat
  if(is.null(nK)){nK = 1/3*ns}

  vhs <- vh# current point

  vk <- c(1,2)

  while (sum(diff(vk)<2) > 0){

    vk <- c(1,floor(ns * (c(1:nK) + runif(nK)) / (nK+2)),ns+1)# stochastic knots

  }

  ##--- sampling start ---##

  for (i in 1 : (nK+1) ){
    ir  <- vk[i]
    id  <- vk[i+1] - vk[i]
    ird <- vk[i+1] - 1

    vyi <- vy[ir:ird]
    vho <- vh[ir:ird]# current (old) point
    vhn <- array(0,dim=c(id,1))# new point

    if (i <= (nK)){

        dhrd1 <- vh[ird+1]# h(r+d+1)

    }

    ##--- finding model & draw candidate ---#

    for (j in 1 : (nite + 1)){

        if (j == 1){

          vhh <- vho# h_hat

        }
        else {
          vhh <- vhn
        }

        vgder2 <- -0.5 * vyi^2 / exp(vhh)# g''(h)
        vgder1 <- -0.5 - vgder2# g'(h)
        vsiga2 <- -1 / vgder2# sig2_ast
        vha <- vhh + vsiga2 * vgder1# h_ast

        if (i <= nK){

          vsiga2[id] <- 1 / (-vgder2[id] + 1/dsig2)
          vha[id] <- vsiga2[id] * (vgder1[id] - vgder2[id]*vhh[id] + dphi * dhrd1/dsig2)

        }

        ##--- simulation smoother ---##

        if (i == 1){

          dh0 <- dh00
          dH20 <- dsig02

        } else {

          dh0 <- dphi * vhs[ir-1]
          dH20 <- dsig2

        }

        da <- dh0
        dP <- dH20
        dH2 <- dsig2
        ve <- array(0,dim=c(id, 1))
        vDinv <- array(0,dim=c(id, 1))
        vK <- array(0,dim=c(id, 1))
        vL <- array(0,dim=c(id, 1))
        vu <- array(0,dim=c(id, 1))

        for (t in 1 : id){
          ve[t] <- vha[t] - da# Kalman filter
          vDinv[t] <- 1 / (dP + vsiga2[t])
          vK[t] <- dphi * dP * vDinv[t]
          vL[t] <- dphi - vK[t]

          da <- dphi * da + vK[t] * ve[t]
          dP <- dphi * dP * vL[t] + dH2
        }

        if (j <= nite){
                    # finding mode
          dr <- 0
          dU <- 0

          t <- id# simulation smoother
          while (t >= 1){

            dC <- dH2 * (1 - dU * dH2)
            deps <- 0
            vu[t] <- dH2 * dr + deps
            dV <- dH2 * dU * vL[t]

            dCinv <- 1 / dC
            dr <- vDinv[t] * ve[t] + vL[t] * dr - dV * dCinv * deps
            dU <- vDinv[t] + vL[t]^2 * dU + dV^2 * dCinv
            t <- t - 1

          }

          du0 <- dH20 * dr
          vhn[1] <- dh0 + du0

          for (t in 1 : (id-1)){

            vhn[t+1] <- dphi * vhn[t] + vu[t]

          }

        } else {
                    # draw candidate
          fl <- 0
          icyc <- 0
          while ((fl  ==  0) && (icyc < 100)){

            dr <- 0
            dU <- 0

            t <- id# simulation smoother
            while (t >=  1){
              dC <- dH2 * (1 - dU * dH2)
              deps <- sqrt(dC) * rnorm(1)
              vu[t] <- dH2 * dr + deps
              dV <- dH2 * dU * vL[t]

              dCinv <- 1 / dC
              dr <- vDinv[t] * ve[t] + vL[t] * dr - dV * dCinv * deps
              dU <- vDinv[t] + vL[t]^2 * dU + dV^2 * dCinv
              t <- t - 1

            }

            dC <- dH20 * (1 - dU * dH20)
            du0 <- dH20 * dr + sqrt(dC) * rnorm(1,0,1)
            vhn[1] <- dh0 + du0


            for (t in 1 : (id-1)){

              vhn[t+1] <- dphi * vhn[t] + vu[t]

            }

          ##--- AR step ---##

            dpron <- sum( -0.5 * (vhh + vyi^2 / exp(vhh)) # g(h_hat)
                          + vgder1 * (vhn - vhh) # g'(h_hat)
                          + 0.5 * vgder2 * (vhn - vhh)^2)# g''(h_hat)

            dposn <- sum(-0.5 * (vhn + vyi^2 / exp(vhn)))

            dfrac <- exp(dposn-dpron)
            if(is.na(dfrac)){dfrac <- 0}


            if (runif(1) < dfrac){
                fl <- 1
            }

            icyc <- icyc + 1

          }

        }

    }

    if (icyc < 100){

      ##--- MH step ---##

        dproo <- sum( -0.5 * (vhh + vyi^2 / exp(vhh)) # g(h_hat)
            + vgder1 * (vho - vhh) # g'(h_hat)
            + 0.5 * vgder2 * (vho - vhh)^2)# g''(h_hat)

        dposo <- sum(-0.5 * (vho + vyi^2 / exp(vho)))

        dfrac <- exp(  dposn + min(c(dposo, dproo))
                    - dposo - min(c(dposn, dpron) ))

        if (runif(1) < dfrac){

            vhs[ir:ird] <- vhn

        }
    }


  }
  return(vhs)
}
