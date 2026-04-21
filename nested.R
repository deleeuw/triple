
library(RSpectra)




nested <- function(ymat,
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
  labd <- eigs_sym(kronecker(cmat, rmat), 1)$values
  if (is.null(xini)) {
    xold <- eckart_young(ymat, p)
  } else {
    xold <- xini
  }
  rold <- wmat * (ymat - xold)
  sold <- loss(rold, rmat, cmat)
  itel <- 1
  repeat {
    gold <- (rmat %*% (wmat * (ymat - xold)) %*% cmat) / wmat
    ynew <- xold + gold / labd
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



