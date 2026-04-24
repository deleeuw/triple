yold <- matrix(rnorm(12), 4, 3)
library(RSpectra)

set.seed(12345)
ymat <- matrix(rnorm(12), 4, 3)
rmat <- crossprod(matrix(rnorm(40), 10, 4)) / 10
cmat <- crossprod(matrix(rnorm(30), 10, 3)) / 10
wmat <- matrix(rnorm(12)^2, 4, 3)

rckr <- kronecker(cmat, rmat)
yvec <- as.vector(ymat)
wvec <- as.vector(wmat)
wdig <- diag(wvec)

hmat <- rckr * outer(wvec, wvec)
# print(matrix(hmat %*% yvec, 4, 3))
# print((wmat^2) * rmat %*% ymat %*% cmat)
# print(wmat * (rmat %*% (wmat * ymat) %*% cmat))