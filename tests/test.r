# load library and data
library(ltvar)
data(jpdata)

# Test fitting a latent threshold model
# Models 1 and 2 with intercept
# Models 3 and 4 without intercept


ltfit1 <- ltvar(y=jpdata,p=1)
ltfit2 <- ltvar(y=jpdata,p=2)
ltfit3 <- ltvar(y=jpdata,p=1,Intercept = FALSE)
ltfit4 <- ltvar(y=jpdata,p=2,Intercept = FALSE)

# Test computing impulse-responses

ltirf1 <- impulse_response(ltfit1,n.ahead = 12)
ltirf2 <- impulse_response(ltfit2,n.ahead = 12)
ltirf3 <- impulse_response(ltfit3,n.ahead = 12)
ltirf4 <- impulse_response(ltfit4,n.ahead = 12)

# Test forecasting

ltfc1 <- forecast(ltfit1,n.ahead = 6)
ltfc2 <- forecast(ltfit2,n.ahead = 6)
ltfc3 <- forecast(ltfit3,n.ahead = 6)
ltfc4 <- forecast(ltfit4,n.ahead = 6)
