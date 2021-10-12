      subroutine hanning(n,g)

      implicit none

      integer, intent(in) :: n

      real, intent(out) :: g(n)

      real, parameter :: pi = 3.1415927

      !$OMP PARALLEL DO
      do i = 1,n
        g(i)  = 0.5*(1+cos(pi*i/n))
      end do
      !$OMP PARALLEL END DO

      end subroutine
