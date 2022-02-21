speedTest <- function(n = 1000, nsim = 1000, seed = 1234, Mfac = 2, f.idx){
  set.seed(seed)


  x <- rnorm(n)
  nTimes2 <- 2^(floor(log2(n)) + 2)
  M <- nTimes2 * Mfac

  if (is.null(k)){
    k <- 0:(M-1)
  }

  start <- Sys.time()
  for (i in 1:nsim){
    a <- goertzel.R(x, k)
  }
  times <- list(Goertzel.R = Sys.time() - start)

  start <- Sys.time()
  for (i in 1:nsim){
    a <- goertzel(x, k)
  }
  times[["Goertzel.F"]] <- Sys.time() - start

  times
}

test.zp <- function(){
  x <- rnorm(100)
  N <- length(x)
  M.fft <- 2^(floor(log2(N))+2)
  f.idx <- 57 # original

  x.zp <- c(x, rep(0, M.fft-N))
  x.fft <- fft(x.zp)[f.idx]
  fac <- 10
  f.idx.r <- (f.idx - 1) * fac + 1

  x.zp.r <- c(x, rep(0, fac*M.fft - N))

  f.grd <- (f.idx.r - (fac - 1)):(f.idx.r + (fac-1))
  x.fft.r <- fft(x.zp.r)[f.idx.r]
  x.g <- goertzel(x.zp, k = f.idx-1)

  start <- Sys.time()
  for (i in 1:1000){
    x.fft.grd <- fft(x.zp.r)[f.grd]
  }
  Sys.time() - start

  start2 <- Sys.time()
  for (i in 1:1000){
    x.g.r <- goertzel(x.zp, k = (f.grd - 1) / fac)
  }
  Sys.time() - start2


  aa <- cbind(x.fft.grd, x.g.r)
  aa[, 2] - aa[, 1]
}

test <- function(){
  set.seed(1234)
  N <- 4
  xx <- rnorm(N)
  xx <- 1:N
  # M <- 2^(floor(log2(N)) + 2)
  M <- 8
  x.zp <- c(xx, rep(0, M-N))

  x.fft <- fft(xx)
  k.idx <- 0:3
  x.gft <- goertzel(xx, k.idx)
  x.gftr <- goertzel.R(xx, k.idx)
  round(matrix(c(x.fft[k.idx+1], x.gft, x.gftr), ncol = 3), 3)
  # round(matrix(c(x.fft[k.idx+1], x.gft, abs(x.fft[k.idx+1] - x.gft)), ncol = 3), 3)

  k.idx <- 0:(M-1)
  x.fft.zp <- fft(x.zp)
  x.gft.zp <- goertzel(x.zp, k.idx)
  x.gftr.zp <- goertzel.R(x.zp, k.idx)
  m.zp <- matrix(c(x.fft.zp[k.idx+1], x.gft.zp, x.gftr.zp), ncol = 3)
  round(cbind(m.zp, m.zp[, 1] - m.zp[, 2], m.zp[, 1] - m.zp[, 3]), 3)

  # round(matrix(c(x.fft.zp[k.idx+1], x.gft.zp, abs(x.fft.zp[k.idx+1] - Conj(x.gft.zp))), ncol = 3), 2)


  fft(1:8)
  x.fft.zp[1:10]
  goertzel(1:8, 1.5)

  x.zp2 <- c(xx, rep(0, 4*M - N))
  x.fft.zp2 <- fft(x.zp2)
  x.gft.zp2 <- goertzel(x.zp, seq(0, M-0.25, by = 0.25))

  round(x.fft.zp2 - x.gft.zp2  , 3)

}


test2N <- function(){
  set.seed(1234)
  N <- floor(runif(1, 100, 200))
  x <- rnorm(N)
  M <- 2^(floor(log2(N))+2)
  x.zp <- c(x, rep(0, M - N))
  x.gzp <- c(x, rep(0, N))

  x.fft <- abs(fft(x.zp))
  x.gft <- abs(goertzel(x.gzp, k = 0:(2*N-1)))

  f.f <- seq(0, 1, by = 1/M)
  f.g <- seq(0, 1, by = 1/(2*N))

  plot(f.f[1:M], x.fft, type='l')
  lines(f.g[1:length(x.gft)], x.gft, col = 'blue')
}
