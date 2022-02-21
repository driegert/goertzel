
#' Generatlized Goertzel Algorithm
#'
#' Goertzel algorithm generalized to non-integer multiples of fundamental frequency
#'
#' @param x Vector containing real-valued data.
#' @param k Vector of real-valued elements containing the frequency components desired,
#'   where k = i/length(x)  (k start at 0 and should be less than length(x)).
#'
#' @details fft(x)[1] = goertzel(x, 0)
#'
#' @references Sysel, Petr, and Pavel Rajmic. "Goertzel algorithm generalized to non-integer multiples of
#'   fundamental frequency." EURASIP Journal on Advances in Signal Processing 2012.1 (2012): 1-8.
#' @export
#'
#' @useDynLib goertzel
goertzel <- function(x, k){
  N <- length(x)
  nk <- length(k)
  out <- .Fortran("goertzel"
                  , x = as.double(x)
                  , N = as.integer(N)
                  , k = as.double(k)
                  , nk = as.integer(nk)
                  , y = complex(nk))
  out$y
}

#' @export
goertzel.R <- function(x, k){
  N <- length(x)
  A = 2*pi*k/N
  B <- 2*cos(A)
  C = exp(-(0+1i) * A)
  D <- exp(-(0+1i)*A * (N-1))

  s0 <- 0; s1 <- 0; s2 <- 0

  for (i in 0:(N-2)){
    s0 = x[i+1] + B * s1 - s2
    s2 <- s1
    s1 <- s0
  }

  s0 <- x[N] + B*s1 - s2
  y <- s0 - s1*C
  y = y*D

  y
}
