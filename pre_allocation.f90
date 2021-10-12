      subroutine pre_allocation(n,window,numzeros,M)

      implicit none

      integer,intent(in) :: n, window
  
      real, intent(in) :: overlap, numzeros

      integer, intent(out) :: M

      integer :: L, D, LL
    
      L = window

      LL = numzeros*L

      M = (LL+L)/2

      end subroutine
