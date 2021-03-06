
#' MCMC computation of latent-threshold models.
#'
#' Estimate a latent-threshold VAR model using a single-move Gibbs sampler as in Nakajima and West (2013)
#'
#' @param y A TxK matrix with the data
#' @param p number of lags
#' @param Intercept Logical flag whether the model contains an intercept
#' @param nreps number of total mcmc draws
#' @param burnin number of burn-in draws
#' @param nKnots number of blocks in the stochastic volatility sampler
#' @param dvb0,dVb0 prior on volatility of betas
#' @param dva0,dVa0 prior on volatiltiy of covariances
#' @param dvh0,dVh0 prior on volatility of variances
#' @param dm0,ds0 prior on intercept
#' @param da0 prior on phi for covariances
#' @param db0 prior on phi for betas
#' @param dg0,dG0 prior on parameters for stochastic volatility
#' @param dk0 latent threshold prior
#' @references Nakajima, J. and M. West (2013) Bayesian Analysis of Latent Threshold Dynamic Models; Journal of Business & Economic Statistics 31 (2), 151-164
#' @export
#' @importFrom MASS ginv
ltvar <- function(y,p=2, Intercept=TRUE,nreps=100,burnin=10,
                  # Priors
                  dvb0 = 20, dVb0 = 0.002, dva0 = 2, dVa0 = 0.002, dvh0 = 2, dVh0 = 0.002,
                  dm0 = 0, ds0 = 1, da0 = 1, db0 = 1, dg0 = 4, dG0 = 0.1, dk0 = 3,
                  nKnots = NULL){
  # Lag data
  xy       <- lagdata(y,p,intercept=Intercept)
  ylagged  <- xy$y
  xlagged  <- xy$x

  # preliminaries
  ns <- nrow(ylagged)
  nk <- ncol(ylagged)
  nl <- p

  constant <- 0
  if(Intercept == TRUE){constant = 1}

  nb <- nk*(nk*nl+constant) # number of beta-coefficients
  na <- nk*(nk-1)/2         # number of parameters in a


  mxLagged <- array(0,dim=c(nrow(xlagged),nb,nk))
  am0minv  <- array(0,dim=c(nrow(xlagged),nk,nk))
  amG2inv  <- array(0,dim=c(nrow(xlagged),nk,nk))

  for(ii in 1:nrow(xlagged)){

    xi <- diag(1,nk) %x% xlagged[ii,]
    mxLagged[ii,,] <- xi
    am0minv[ii,,]  <- diag(1,nk)
  }

  mb <- array(0,dim=c(ns,nb))
  mbd <- mb

  ma <- array(0,dim=c(ns,na))
  mad <- ma

  mh <- array(0,dim=c(ns,nk))

  mSigb <- diag(1,nb) * 0.01
  mSiga <- diag(1,na) * 0.01
  mSigh <- diag(1,nk) * 0.01

  vmub <- array(0,dim=c(nb,1))
  vmua <- array(0,dim=c(na,1))
  vgam <- array(1,dim=c(1,nk))

  mPhib <- diag(1,nb) * 0.95
  mPhia <- diag(1,na) * 0.95
  mPhih <- diag(1,nk) * 0.95

  vdb <- array(0.0,dim=c(nb,1))
  vda <- array(0.0,dim=c(na,1))

  vidB <- seq(1:nb)-1

  if( Intercept == TRUE ){

    x1 <- seq(1:nk)-1
    vidi <- x1 * ( nk * nl + 1 )

  }

  dnub <- dvb0 + ns - nl
  dnua <- dva0 + ns - nl
  dnuh <- dvh0 + ns - nl


  # Declare Variables for Storage

  mbSave  <- array(0, dim=c( nreps - burnin, ns, nb ))
  maSave  <- array(0, dim=c( nreps - burnin, ns, na ))
  mhSave  <- array(0, dim=c( nreps - burnin, ns, nk ))

  mPhiasave <- array(0, dim=c( nreps - burnin, na, na ))
  mPhibsave <- array(0, dim=c( nreps - burnin, nb, nb ))
  mPhihsave <- array(0, dim=c( nreps - burnin, nk, nk ))

  vdbSave <- array(0, dim=c( nreps - burnin, nb, 1 ))
  vdaSAve <- array(0, dim=c( nreps - burnin, na, 1 ))

  vmuasave <- array(0, dim=c( nreps - burnin, na, 1 ))
  vmubsave <- array(0, dim=c( nreps - burnin, nb, 1 ))

  vgamsave <- array(0, dim=c( nreps - burnin, 1, nk ))

  mSigasave <- array(0, dim = c( nreps - burnin, na, na ))
  mSigbsave <- array(0, dim = c( nreps - burnin, nb, nb ))
  mSighsave <- array(0, dim = c( nreps - burnin, nk, nk))


  # Start sampling
  for(ii in 1:nreps){
    print(ii)

    # Sample beta
    xtmp <- fSampLTBv(my = ylagged, amX = mxLagged, mPhi = mPhib, vmu = vmub,
                      mSig = mSigb, amG2inv = am0minv, vd = vdb, mbo = mb,
                      stabletest = FALSE, Intercept = FALSE, nl = nl)

    mb <- xtmp
    mbd  <- mb

    # Sample Phib, Sigmab, mub

    xtmp <- fSampSPM(mp = mb, mSig = mSigb, mPhi = mPhib, vmu = vmub,
                     dnu = dnub, dV0 = dVb0, da0 = da0, db0 = db0,
                     dm0 = dm0, ds0 =ds0, vd = vdb, dk0 = dk0)

    vmub <- xtmp$vmu
    mPhib <- xtmp$mPhi
    mSigb <- xtmp$mSig

    # Sample threshold for b

    vdbc <- abs(vmub) + ( dk0 * sqrt( diag( mSigb %*% MASS::ginv( diag(1,nb) - mPhib^2 ))))

    xtmp <- fSampLT( my = ylagged, amX = mxLagged, mp = mb,am0minv, vd = vdb, vdc = vdbc)
    vdb  <- xtmp


    for(jj in 1:ns){

      mbd[jj,abs(mbd[jj,])<vdb] <- 0

    }


    # sampling a

    myh   <- array(0,dim=c(ns,nk))
    amXhi <- array(0,dim=c(ns,na,nk))

    for(jj in 1:ns){

      myh[jj,]      <- ylagged[jj,]-mbd[jj,]%*%mxLagged[jj,,]
      amXhi[jj,,]   <- t(fxh(myh[jj,],nk,na))

      xx <- as.vector(exp(-mh[jj,])/vgam)
      amG2inv[jj,,] <- diag(xx)

    }

    xtmp <- fSampLTBv(myh,amXhi,mPhia,vmua,mSiga,amG2inv,vda,ma,
                      stabletest = FALSE, Intercept = FALSE, nl = nl)
    ma   <- xtmp
    mad  <- ma

    for(jj in 1:ns){

      mad[jj,abs(mad[jj,])<vda] <- 0

    }

    # Sampling sigmaa, phia, mua

    xtmp  <- fSampSPM(ma,mSiga,mPhia,vmua,dnua,dVa0,da0,db0,dm0,ds0,vda,dk0)
    mSiga <- xtmp$mSig
    mPhia <- xtmp$mPhi
    vmua  <- xtmp$vmu


    # Sampling threshold of a

    vdac <- abs(vmua) + dk0 * sqrt( diag( mSiga %*% MASS::ginv( diag(1,na) - mPhia )))

    xtmp <- fSampLT(myh,amXhi,ma,amG2inv,vda,vdac)
    vda  <- xtmp

    # Sampling h
    mya <- array(0,dim=c(ns,nk,nk))

    for(jj in 1:ns){

      mya[jj,,] <- myh[jj,] %*% t(fAt(mad[jj,],nk))

    }

    mSh0 <- mSigh %*% MASS::ginv(diag(1,nk)-mPhih^2)

    for( jj in 1:nk){

      xtmp <- fSVSamp(mya[,jj,jj]/sqrt(vgam[jj]),mh[,jj],mPhih[jj,jj],mSigh[jj,jj],0,mSh0[jj,jj],nKnots)
      mh[,jj] <- xtmp

    }

    # Sampling Sigmah, phih, gammah

    xtmp <- fSampSPG(mya, mh, mPhih, vgam, dnuh, dVh0, da0, db0, dg0, dG0)
    mSigh <- xtmp$mSig
    mPhih <- xtmp$mPhi
    vgam  <- xtmp$vgam


    for( jj in 1:ns){

      mA <- fAt(mad[jj,],nk)
      mAinv <- solve(mA)
      xtmp <- mA %*% diag( as.vector( exp( - mh[jj,] ) / vgam ) ) %*% t(mA)
      am0minv[jj,,] <- xtmp

    }

    # Store Variables

    if( ii > burnin){
      vmubsave[ii - burnin,,] <- vmub
      vmuasave[ii - burnin,,] <- vmua

      mPhiasave[ii - burnin,,] <- mPhia
      mPhibsave[ii - burnin,,] <- mPhib
      mPhihsave[ii - burnin,,] <- mPhih

      mSigasave[ii - burnin,,] <- mSiga
      mSigbsave[ii - burnin,,] <- mSigb
      mSighsave[ii - burnin,,] <- mSigh

      mbSave[ii - burnin,,] <- mb
      maSave[ii - burnin,,] <- ma
      mhSave[ii - burnin,,] <- mh

      vgamsave[ii - burnin,,] <- vgam

      vdbSave[ii - burnin,,] <- vdb
      vdaSAve[ii - burnin,,] <- vda

    }
  }
  draws = list(mb = mbSave, ma = maSave, mh = mhSave, vdb = vdbSave,vda = vdaSAve,
               mPhia = mPhiasave, mPhib = mPhibsave, mPhih = mPhihsave,
               vmua = vmuasave, vmub = vmubsave,
               mSiga = mSigasave, mSigb = mSigbsave, mSigh = mSighsave,
               vgam = vgamsave)

  model_info = list(p=p,intercept=Intercept)
  data_info  = list(data = y)
  retlist <- structure(list(mcmc_draws=draws,
                            data_info = data_info,
                            model_info=model_info),class="ltvar")
  return(retlist)

}
