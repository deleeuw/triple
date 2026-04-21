
library(RSpectra)

set.seed(12345)
ymat <- matrix(rnorm(12), 4, 3)
yold <- matrix(rnorm(12), 4, 3)
rmat <- crossprod(matrix(rnorm(40), 10, 4)) / 10
cmat <- crossprod(matrix(rnorm(30), 10, 3)) / 10
wmat <- matrix(rnorm(12)^2, 4, 3)

rckr <- kronecker(cmat, rmat)
yvec <- as.vector(ymat)
wvec <- as.vector(wmat)
wdig <- diag(wvec)

print(matrix(rckr %*% yvec, 4, 3))
print(rmat %*% ymat %*% cmat)

print(matrix((rckr %*% (wvec * yvec)) / wvec, 4, 3))
print((rmat %*% (wmat * ymat) %*% cmat) / wmat)

checker <- function(x, y, rmat, cmat, wmat) {
  sx <- loss(x, rmat, cmat, wmat)
  so <- loss(y, rmat, cmat, wmat)
  lb <- eigs_sym(rmat, 1)$values * eigs_sym(cmat, 1)$values
  go <- rmat %*% y %*% cmat
  vo <- y - go/lb
  s1 <- so + 2 * sum(go * (x - y)) + lb * sum((x - y)^2)
  s2 <- so + lb * sum((x - vo)^2) - sum(go^2)/lb
  print(c(sx, s1, s2))
  print(c(loss(x - y, rmat, cmat), lb * sum((x - y)^2)))
}

loss <- function(x, rmat, cmat) {
  return(sum((rmat %*% x %*% cmat) * x))
}

triple <- function(ymat,
                   wmat = array(1, dim(ymat)),
                   rmat = diag(nrow(ymat)),
                   cmat = diag(ncol(ymat)),
                   p = 2,
                   zini = NULL,
                   itmax = 100,
                   inmax = 1,
                   eps = 1e-6,
                   verbose = TRUE) {
  nr <- nrow(ymat)
  nc <- ncol(ymat)
  labd <- eigs_sym(rmat, 1)$values * eigs_sym(cmat, 1)$values
  if (is.null(zini)) {
    zold <- eckart_young(ymat, p)
  } else {
    zold <- zinit
  }
  rold <- wmat * (ymat - zold)
  sold <- loss(rold, rmat, cmat)
  itel <- 1
  repeat {
    gold <- rmat %*% zold %*% cmat
    ynew <- zold - (gold / wmat) / labd
    znew <- eckart_young(ynew, p)
    rnew <- wmat * (ymat - znew)
    snew <- loss(rnew, rmat, cmat)
    epsi <- max(abs(zold - znew))
    if (verbose) {
      cat("itel", formatC(itel, digits = 3, format = "d"),
          "sold", formatC(sold, digits = 10, width = 15, format = "f"),
          "snew", formatC(snew, digits = 10, width = 15, format = "f"),
          "epsi", formatC(epsi, digits = 10, width = 15, format = "f"),
          "\n")
    }
    if ((itel == itmax) || (epsi < eps)) {
      break
    }
    itel <- itel + 1
    zold <- znew
    sold <- snew
  }
}

eckart_young <- function(x, p) {
  h <- svd(x, nu = p, nv = p)
  return(tcrossprod(h$u, h$v %*% diag(h$d[1:p])))
}



