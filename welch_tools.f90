      module welch_tools

      implicit none

      real, parameter :: pi = 4.0 *ATAN(1.0)
      
      contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       1) hanning        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function hanning(n)

        integer, intent(in) :: n

        real, dimension(n) :: hanning

        integer :: i

        !$OMP PARALLEL DO
        do i = 1,n
          hanning(i)  = 0.5*(1-cos(2*pi*i/(n+1)))
        end do
        !$OMP PARALLEL END DO
                
        end function hanning

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       2) detrend        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function detrend(n,x,f)

        integer, intent(in) :: n

        real, dimension(n), intent(in) :: x, f

        real, dimension(n) :: detrend,y

        real :: m, b, sumxy, sumx, sumy, sumxx

        sumx = sum(x)
        sumy = sum(f)
        sumxx = sum(x*x)
        sumxy = sum(x*f)

        m = (n*sumxy-sumx*sumy)/(n*sumxx-(sumx)**2)
        b = (sumy*sumxx-sumx*sumxy)/(n*sumxx-sumx**2)

        y = m*x+b
        detrend = f-y

        end function detrend

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       3) prealocation   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function pre_allocation(n,window,numzeros)

        integer,intent(in) :: n, window

        real, intent(in) :: numzeros

        integer :: pre_allocation ! -> M

        integer :: L, D, LL

        L = window
        LL = numzeros*L
        pre_allocation = (LL+L) ! -> M

        end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       4) FFT            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function fft(n,f,m)

        integer,intent(in) :: m, n

        real, dimension(n), intent(in) :: f

        complex, dimension(m) :: fft

        complex, parameter :: i = (0,1)

        real, dimension(m,m) :: ac, as

        real, dimension(m) :: xx

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
          fft = matmul(ac,xx)+i*matmul(as,xx)
        else 
          xx(1:n) = f
          fft = matmul(ac,xx)+i*matmul(as,xx)
        end if        

        end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       5) spatial welch  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!
        function spatial_welch(n,x,y,window,overlap,delta,numzeros,M)

        integer, intent(in) :: n, window, M
  
        real, dimension(n), intent(in) :: x, y

        real, intent(in) :: overlap, delta, numzeros

        integer :: L, D, LL, k

        real :: spatial_welch(int(M/2),2), dy(n), U               
  
        complex :: ak(M)

        dy = detrend(n,x,y) 

        L = window

        D = overlap*L

        LL = numzeros*L
       
        U = sum(hanning(L)**2)/L

        ak = fft(n,detrend(L,x(1:L),y(1:L))*hanning(L),M)

        !$OMP PARALLEL DO
        do k = 2,int(n/window)          
          ak = ak + fft(L+D,detrend(L+D,x(1+L*(k-1)-D:k*L),y(1+L*(k-1)-D:k*L)) &
              *hanning(L+D),M)
        end do
        !$OMP PARALLEL END DO

        ak = (abs(ak/M)*M**2)/(L*U)
   
        !$OMP PARALLEL DO
        do k = 1,int(M/2)
          spatial_welch(k,1) = k/(M*delta)
          spatial_welch(k,2) = ak(k)
        end do
        !$OMP PARALLEL END DO
                       
        end function

      end module welch_tools
