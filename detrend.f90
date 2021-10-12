      subroutine detrend(n,x,f,g)

      implicit none

      integer, intent(in) :: n

      real, intent(in) :: x(n), f(n)

      real, intent(out) :: f

      real :: m, b, sumxy, sumx, sumy, sumxx, y(n)

      sumx = sum(x)
 
      sumy = sum(f)
 
      sumxx = sum(x*x)

      sumxy = sum(x*f)

      m = (n*sumxy-sumx*sumy)/(n*sumxx-(sumx)**2)

      b = (sumy*sumxx-sumx*sumxy)/(n*sumxx-sumx**2)

      y = m*x+b

      g = f-y

      end subroutine
