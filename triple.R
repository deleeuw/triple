# example 1

set.seed(12345)
rmat <- crossprod(matrix(rnorm(400), 100, 4))/100
cmat <- crossprod(matrix(rnorm(300), 100, 3))/100
a <- rnorm(4)^2
b <- rnorm(3)^2
c <- rnorm(4)^2
d <- rnorm(3)^2
m <- 3
z <- as.vector(m + outer(a, b, "+") + outer(c, d))
ymat <- matrix(rpois(12, z), 4, 3)
wmat <- 1 / sqrt(xmat)
aini <- 1:4
bini <- 1:3
cini <- 1:4
dini <- 1:3
mini <- 1
ztilde <- mini + outer(aini, bini, "+") + outer(cini, dini)
xtilde <- wmat * (ymat - ztilde)
gtilde <- rmat %*% xtilde %*% cmat
vtilde <- xtilde - vtilde
krx <- matrix(kr %*% x, 4, 3)
l <- sum(x * (kr %*% x))
