# load library and data
library(ltvar)
data(jpdata)

# Test fitting a latent threshold model

ltfit1 <- ltvar(y=jpdata,p=1)
ltfit2 <- ltvar(y=jpdata,p=2)

# Test computing impulse-responses

ltirf1 <- impulse_response(ltfit1,n.ahead = 12)
ltirf2 <- impulse_response(ltfit2,n.ahead = 12)

# Test forecasting

ltfc1 <- forecast(ltfit1,n.ahead = 6)
ltfc2 <- forecast(ltfit2,n.ahead = 6)
