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

checker <- function(x, y, rmat, cmat, wmat) {
  sx <- loss(x, rmat, cmat)
  so <- loss(y, rmat, cmat)
  lb <- eigs_sym(rmat, 1)$values * eigs_sym(cmat, 1)$values
  go <- rmat %*% y %*% cmat
  vo <- y - go/lb
  s1 <- so + 2 * sum(go * (x - y)) + lb * sum((x - y)^2)
  s2 <- so + lb * sum((x - vo)^2) - sum(go^2)/lb
  print(c(sx, s1, s2))
  print(c(loss(x - y, rmat, cmat), lb * sum((x - y)^2)))
}