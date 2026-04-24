loss <- function(x, rmat, cmat) {
  return(sum((rmat %*% x %*% cmat) * x))
}

decompose <- function(xmat, rmat, cmat) {
  n <- nrow(xmat)
  m <- ncol(xmat)
  pmat <- matrix(colSums(rmat), n, n, byrow = TRUE) / sum(rmat)
  qmat <- matrix(rowSums(cmat), m, m) / sum(cmat)
  px <- pmat %*% xmat
  xq <- xmat %*% qmat
  mu <- px %*% qmat
  ae <- px - mu
  eb <- xq - mu
  ab <- xmat - ae - eb - mu
  print()
}

eckart_young <- function(x, p) {
  h <- svd(x, nu = p, nv = p)
  if (p == 1) {
    return(h$d[1] * outer(drop(h$u), drop(h$v)))
  } else {
    return(tcrossprod(h$u, h$v %*% diag(h$d[1:p])))
  }
}

column_adjust <- function(x, p) {
  x <- t(apply(x, 1, function(z)
    z - mean(z)))
  return(eckart_young(x, p))
}

row_adjust <- function(x, p) {
  x <- apply(x, 2, function(z)
    z - mean(z))
  return(eckart_young(x, p))
}

double_center <- function(x, p) {
  rm <- apply(x, 1, mean)
  cm <- apply(x, 2, mean)
  mm <- mean(x)
  x <- x - outer(rm, cm, "+") + mm
  return(eckart_young(x, p))
}

hankel <- function(x, p) {
  n <- nrow(x)
  m <- 1:n
  z <- array(0, dim(x))
  for (i in (1 - n):(n - 1)) {
    y <- ifelse(outer(m, rev(m), "-") == i, 1, 0)
    z <- z + (sum(y * x) / sum(y)) * y
  }
  return(z)
}

toeplitz <- function(x) {
  n <- nrow(x)
  m <- 1:n
  z <- array(0, dim(x))
  for (i in (1 - n):(n - 1)) {
    y <- ifelse(outer(m, m, "-") == i, 1, 0)
    z <- z + (sum(y * x) / sum(y)) * y
  }
  return(z)
}

symmetric <- function(x) {
  return((x + t(x)) / 2)
}

anti_symmetric <- function(x) {
  return((x - t(x)) / 2)
}

additive <- function(x) {
  rm <- apply(x, 1, mean)
  cm <- apply(x, 2, mean)
  mm <- mean(x)
  return(outer(rm, cm, "+") - mm)
}