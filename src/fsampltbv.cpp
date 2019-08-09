// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <cmath>
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <Rcpp.h>


using namespace std;
using namespace Rcpp;

//' @export
//' @title samples the parameter
//' This function samples the VARS-parameter of a latent-threshold VAR using a single move algorithm.
//' @param my LHS of data
//' @param amX RHS of data
//' @param mPhi transition matrix for VAR-Parameters
//' @param vmu mean of transition matrix
//' @param mSig variance-covariance matrix for standard deviation of parameters
//' @param amG2inv a KxKxT array with the stochastic volatility
//' @param vd latent threshold
//' @param mbo old betas
//' @param stabletest Test if candidate draw is stable (currently not implemented)
//' @param Intercept whether the model contains an intercept or not.
//' @param nl number of lags

// [[Rcpp::depends("RcppArmadillo")]]
// [[Rcpp::export]]

arma::mat fSampLTBv(arma::mat my,arma::cube amX,arma::mat mPhi,arma::colvec vmu,arma::mat mSig,arma::cube amG2inv,arma::vec vd, arma::mat mbo,int stabletest,int Intercept,int nl){

  // Declaring variables
  int np = mPhi.n_cols;
  int ns = my.n_rows;
  int nk = my.n_cols;
  int nb = mbo.n_cols;
  arma::mat mba=mbo;
  arma::mat mSigi = pinv(mSig);
  arma::mat xx = arma::eye<arma::mat>(np,np);
  arma::mat mSigp = mSigi*(xx+mPhi*mPhi);
  arma::colvec vpm1 = (xx-mPhi)*vmu;
  arma::colvec vpm2 = (xx-2*mPhi+mPhi*mPhi)*vmu;
  arma::mat mSig0i = (xx-mPhi*mPhi)*mSigi;
  arma::colvec vbh=arma::zeros(nb);
  arma::colvec vba1=arma::zeros(nb);
  arma::mat mSighi,mSigh;


  for(int i=0;i<(ns-1);i++){
    arma::mat amxtmp   = amX(arma::span(i),arma::span::all,arma::span::all);
    arma::mat amgtmp   = amG2inv(arma::span(i),arma::span::all,arma::span::all);
    arma::rowvec ytmp  = my.row(i);
    arma::rowvec mbtmp = mba.row(i+1);
    arma::rowvec vbo   = mba.row(i);
    // Calculate expected mean of candidate
    if(i==0){

      mSighi = amxtmp*amgtmp*amxtmp.t()+mSig0i+mSigi*mPhi*mPhi;
      mSigh  = pinv(mSighi);
      vbh = mSigh*(amxtmp*amgtmp*ytmp.t()+mSig0i*vmu+mSigi*mPhi*(arma::trans(mbtmp)-vpm1));

    }
    else if(i < (ns-1)){

      mSighi = amxtmp*amgtmp*amxtmp.t()+mSigp;
      mSigh  = pinv(mSighi);
      vbh    = mSigh*(amxtmp*amgtmp*ytmp.t()+mSigi*(mPhi*(vba1+arma::trans(mbtmp))+vpm2));

    }
    else{
      mSighi = amxtmp*amgtmp*amxtmp.t()+mSigi;
      mSigh  = pinv(mSighi);
      vbh = mSigh*(amxtmp*amgtmp*ytmp.t()+mSigi*(mPhi*vba1+vpm1));
    }
    vba1=vbh;
    // draw candidate
    arma::vec vbn = vbh + arma::chol(mSigh).t() * as<arma::colvec>(rnorm(nb));

    arma::colvec dho,dhn;
    dho = -0.5*log(arma::det(mSigh));
    dhn = dho;
    dhn = dhn - 0.5 * arma::trans(vbn-vbh) * mSighi * (vbn-vbh);
    dho = dho - 0.5 * arma::trans(arma::trans(vbo)-vbh) * mSighi * (arma::trans(vbo)-vbh);

    arma::mat mxh = amX(arma::span(i),arma::span::all,arma::span::all);
    // test if candidate draws are smaller than threshold
    for(int j=0;j<nk;j++){
      for(int k=0;k<nb;k++){
        double vbntest = vbn(k);
        double vdtest = vd(k);
        if(abs(vbntest) < vdtest){
          mxh(k,j)=0;
        }
      }
    }
    if(i==0){
      mSighi = mxh*amgtmp*mxh.t()+mSig0i+mSigi*mPhi*mPhi;
      mSigh  = pinv(mSighi);
      vbh = mSigh*(mxh*amgtmp*ytmp.t()+mSig0i*vmu+mSigi*mPhi*(arma::trans(mbtmp)-vpm1));
    }
    else if(i<(ns-1)){
      mSighi = mxh*amgtmp*mxh.t()+mSigp;
      mSigh  = pinv(mSighi);
      vbh    = mSigh*(mxh*amgtmp*ytmp.t()+mSigi*(mPhi*(vba1+arma::trans(mbtmp))+vpm2));
    }
    else{
      mSighi = mxh*amgtmp*mxh.t()+mSigi;
      mSigh  = pinv(mSighi);
      vbh    = mSigh*(mxh*amgtmp*ytmp.t()+mSigi*(mPhi*vba1+vpm1));
    }
    arma::colvec dln = -0.5*log(arma::det(mSigh))-0.5*arma::trans(vbn-vbh)*mSighi*(vbn-vbh);

    mxh = amX(arma::span(i),arma::span::all,arma::span::all);
    // test if candidate draws are smaller than threshold
    for(int j=0;j<nk;j++){
      for(int k=0;k<nb;k++){
        double vbntest = vbo(k);
        double vdtest = vd(k);
        if(abs(vbntest) < vdtest){
          mxh(k,j)=0;
        }
      }
    }
    if(i==0){
      mSighi = mxh*amgtmp*mxh.t()+mSig0i+mSigi*mPhi*mPhi;
      mSigh  = pinv(mSighi);
      vbh = mSigh*(mxh*amgtmp*ytmp.t()+mSig0i*vmu+mSigi*mPhi*(arma::trans(mbtmp)-vpm1));
    }
    else if(i<(ns-1)){
      mSighi = mxh*amgtmp*mxh.t()+mSigp;
      mSigh  = pinv(mSighi);
      vbh    = mSigh*(mxh*amgtmp*ytmp.t()+mSigi*(mPhi*(vba1+arma::trans(mbtmp))+vpm2));
    }
    else{
      mSighi = mxh*amgtmp*mxh.t()+mSigi;
      mSigh  = pinv(mSighi);
      vbh    = mSigh*(mxh*amgtmp*ytmp.t()+mSigi*(mPhi*vba1+vpm1));
    }
    arma::colvec dlo = -0.5*log(arma::det(mSigh))-0.5*arma::trans(arma::trans(vbo)-vbh)*mSighi*(arma::trans(vbo)-vbh);
    double dfrac = (exp(dln-dhn-dlo+dho)).eval()(0,0);
    double dftest = (arma::randu(1)).eval()(0,0);

    if(dftest < dfrac){
      mba(i,arma::span::all) = arma::trans(vbn);
      vba1 = vbn;
    }
    else{
      vba1 = mba(i,arma::span::all).t();
    }

  }
  return(mba);
}


