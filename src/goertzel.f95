! Figure 4 of:
! Sysel, Petr, and Pavel Rajmic. "Goertzel algorithm generalized to 
! non-integer multiples of fundamental frequency." EURASIP Journal 
! on Advances in Signal Processing 2012.1 (2012): 1-8.
subroutine goertzel(x, N, k, nk, y)
  implicit none
  integer :: N, nk, i
  real(kind=8) :: x(0:(N-1)), k(nk), A(nk), B(nk), s0(nk), s1(nk), s2(nk)
  complex(kind=8) :: y(nk), C(nk), D(nk)
  complex(kind=8), parameter :: eye = cmplx(0.d0, 1.d0)
  real(kind=8), parameter :: pi = 4.0d0*atan(1.d0)
  
  A = 2*pi*k/N
  B = 2*cos(A)
  C = exp(-eye*A)
  D = exp(-eye * A * (N-1))
  
  s0 = 0.d0
  s1 = 0.d0
  s2 = 0.d0
  
  do i = 0, N-2
    s0 = x(i) + B*s1 - s2
    s2 = s1
    s1 = s0
  end do
  
  s0 = x(N-1) + B*s1 - s2
  y = s0 - s1*C
  y = y*D
end subroutine goertzel
