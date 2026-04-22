
library(RSpectra)

source("auxiliary.R")

nested <- function(ymat,
                   wmat = array(1, dim(ymat)),
                   rmat = diag(nrow(ymat)),
                   cmat = diag(ncol(ymat)),
                   func = eckart_young,
                   p = 2,
                   xini = NULL,
                   itmax = 100000,
                   inmax = 10,
                   eps = 1e-6,
                   ips = 1e-6,
                   verbose = TRUE) {
  nr <- nrow(ymat)
  nc <- ncol(ymat)
  labd <- eigs_sym(kronecker(cmat, rmat), 1)$values
  wmax <- max(wmat^2)
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
    gold <- (rmat %*% (wmat * (ymat - xold)) %*% cmat) / wmat
    ymid <- xold + gold / labd
    xmid <- xold
    ntel <- 1
    repeat {
      ydob <- xmid + ((wmat^2) / wmax) * (ymid - xmid)
      if (length(formals(func)) == 1) {
        xnew <- func(ydob)
      } else {
        xnew <- func(ydob, p)
      }
      ipsi <- max(abs(xmid - xnew))
      if ((ntel == inmax) || (ipsi < ips)) {
        break
      }
      ntel <- ntel + 1
      xmid <- xnew
    }
    rnew <- wmat * (ymat - xnew)
    snew <- loss(rnew, rmat, cmat)
    epsi <- max(abs(xold - xnew))
    if (verbose) {
      cat("itel", formatC(itel, digits = 3, format = "d"),
          "ntel", formatC(ntel, digits = 3, format = "d"),
          "sold", formatC(sold, digits = 10, width = 12, format = "f"),
          "snew", formatC(snew, digits = 10, width = 12, format = "f"),
          "epsi", formatC(epsi, digits = 10, width = 12, format = "f"),
          "ipsi", formatC(epsi, digits = 10, width = 12, format = "f"),
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


