      subroutine fft(n,f,m,g)

      implicit none

      integer,intent(in) :: m, n

      real,intent(in) :: f(n)

      complex, intent(out) :: g(m)

      complex, parameter :: i = (0,1)

      real, parameter :: pi = 3.1415927

      real :: ac(m,m),as(m,m), xx(m)

      integer :: j, k

      !$OMP PARALLEL DO
      do k = 1,m
        do j = 1,m
          ac(j,k) = cos(-2*pi*(k-1)*(j-1)/m)
          as(j,k) = sin(-2*pi*(k-1)*(j-1)/m)
        end do
      end do
      !$OMP PARALLEL END DO      
  
      if (n.ge.m) then
        xx = f(1:m)
        g = matmul(ac,xx)+i*matmul(as,xx)
      else 
        xx(1:n) = f
        g = matmul(ac,xx)+i*matmul(as,xx)
      end if        

      end subroutine
