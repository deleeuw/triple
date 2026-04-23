
library(RSpectra)

source("auxiliary.R")

triple <- function(ymat,
                   wmat = array(1, dim(ymat)),
                   rmat = diag(nrow(ymat)),
                   cmat = diag(ncol(ymat)),
                   func = eckart_young,
                   p = 2,
                   xini = NULL,
                   itmax = 100000,
                   eps = 1e-6,
                   verbose = TRUE) {
  nr <- nrow(ymat)
  nc <- ncol(ymat)
  wvec <- as.vector(wmat) 
  hmat <- kronecker(cmat, rmat) * outer(wvec, wvec)
  labd <- eigs_sym(hmat, 1)$values
  if (is.null(xini)) {
    if (length(formals(func)) == 1) {
      xold <- func(ymat)
    } else {
      xold <- func(ymat, p)
    }
  } else {
    xold <- xini
  }
  rold <- wmat * (ymat - xold)
  sold <- loss(rold, rmat, cmat)
  itel <- 1
  repeat {
    gold <- wmat * (rmat %*% (wmat * (ymat - xold)) %*% cmat) 
    ynew <- xold + gold / labd
    if (length(formals(func)) == 1) {
      xnew <- func(ynew)
    } else {
      xnew <- func(ynew, p)
    }
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
  return(list(x = xnew, loss = snew, itel = itel))
}

