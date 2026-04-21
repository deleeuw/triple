
library(RSpectra)

set.seed(12345)
ymat <- matrix(rnorm(12), 4, 3)




triple <- function(ymat,
                   wmat = array(1, dim(ymat)),
                   rmat = diag(nrow(ymat)),
                   cmat = diag(ncol(ymat)),
                   p = 2,
                   xini = NULL,
                   itmax = 100,
                   inmax = 1,
                   eps = 1e-6,
                   verbose = TRUE) {
  nr <- nrow(ymat)
  nc <- ncol(ymat)
  wvec <- as.vector(wmat) 
  hmat <- kronecker(cmat, rmat) * outer(wvec, wvec)
  labd <- eigs_sym(hmat, 1)$values
  if (is.null(xini)) {
    xold <- eckart_young(ymat, p)
  } else {
    xold <- xini
  }
  rold <- wmat * (ymat - xold)
  sold <- loss(rold, rmat, cmat)
  itel <- 1
  repeat {
    gold <- (wmat ^ 2) * (rmat %*% (ymat - xold) %*% cmat) 
    ynew <- xold + 2 * gold / labd
    xnew <- eckart_young(ynew, p)
    rnew <- wmat * (ymat - xnew)
    snew <- loss(rnew, rmat, cmat)
    epsi <- max(abs(xold - xnew))
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
    xold <- xnew
    sold <- snew
  }
}

loss <- function(x, rmat, cmat) {
  return(sum((rmat %*% x %*% cmat) * x))
}

eckart_young <- function(x, p) {
  h <- svd(x, nu = p, nv = p)
  return(tcrossprod(h$u, h$v %*% diag(h$d[1:p])))
}



