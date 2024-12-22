! Rochester Potential Functions (rochester.f08)
!
! Author:  D. Younis
!          University of Rochester
!          Department of Physics
!
! Written: 6/19/2020
! Revised: 7/3/2023

module rochester

use prec
use math

implicit none

real(num), parameter, private :: half3 = 1.50_num

contains

pure real(num) function Vsc(Z,x,s)
    real(num), intent(in) :: Z(2), x(:), s
    Vsc = Z(1)*Z(2)/sqrt(sum(x**2) + s**2)
end function Vsc

pure real(num) function DVsc(Z,x,s,j)
    real(num), intent(in) :: Z(2), x(:), s
    integer, intent(in) :: j
    DVsc = -x(j)*Z(1)*Z(2)/(sum(x**2) + s**2)**half3
end function DVsc

end module rochester
