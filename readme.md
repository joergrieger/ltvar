# Overview

The LTVAR package estimates the Latent Threshold Vector Autoregressive Model of Nakajima and West (2013)

# Model

The LTVAR package estimates the following model:

\[
y_t = c_t +B_{1t}y_{t-1}+B_{2t}y_{t-2}+B_{pt}Y_{t-p}+u_t\quad u_t\sim N(0,\Sigma_t)
\]
Construct the vector $b_t$ by stacking the set of $c_t$ and $B_j$ by rows and by order $j=1,\ldots,p$. Define the matrix $X_t=I\otimes(1,y'_{t-1},\ldots,y'_{t-p})$ we have
\[
y_t=X_tb_t+u_t.
\]
With the elements of $b_t$ as $b_{it}=\beta_{it}I(|\beta_{it}\geq d_i)$. The evolution of $\beta_{t}$ is given below.


The time-varying covariane matrices are modelled as follows:
\[
\Sigma_t = A_t^{-1}\Lambda^2_t(A'_t)^{-1}, A_t=\left(\begin{array}{cccc}
1&0&\cdots&0\\
a_{21,t}&\ddots&\ddots&\vdots\\
\vdots&\ddots&\ddots&0\\
a_{m1,t}&\cdots&a_{m,m-1,t}&1
\end{array}\right),
\Lambda_t=\left(\begin{array}{cccc}
\sigma_{1t}&0&\cdots&0\\
0&\ddots&\ddots&\vdots\\
\vdots&\ddots&\ddots&0\\
0&\cdots&0&\sigma_{mt}
\end{array}\right)
\]

Let $a_t$ be the vector of the strictly lower triangular elements of $A_t$ and define $h_t=(h_{1t},\ldots,h_{mt})$ where $h_{jt}=\log(\sigma^2_{jt})$. Then dynamics of covariances and variances are specified jointly with the Time-Varying coefficient matrix $\beta_t$ as
$$\begin{aligned}
\beta_t=\mu_{\beta}+\Phi_\beta(\beta_{t-1}-\mu_\beta)+\eta_{\beta,t}\\
a_t=\mu_a+\Phi_a(a_{t-1}-\mu_a)+\eta_{a,t}\\
h_t=\mu_h+\Phi_h(h_{t-1}-\mu_h)+\eta_{h,t}
\end{aligned}
$$

# Installation

To install the package you need the devtools package. If you don't have the devtools package, you can install it with

    install.packages("devtools")
    
Once you have installed the devtools package you can install the ltvar package with

    library(devtools)
    devtools::install_github('joergrieger/ltvar')
