# load library and data
library(ltvar)
data(jpdata)

# Test fitting a latent threshold model

ltfit <- ltvar(y=jpdata,p=3)

# Test computing impulse-responses

ltirf <- impulse_response(ltfit,n.ahead = 24)
