! Rochester Potential Functions (rochester.f08)
!
! Author:  D. Younis
!          University of Rochester
!          Department of Physics
!
! Written: 6/19/2020
! Revised: 8/23/2022

module rochester

use prec
use math

implicit none

contains

pure real(num) function Vroc(Z1, Z2, x, y, ac)
    implicit none
    real(num), intent(in) :: Z1, Z2, x, y, ac

    Vroc = (Z1*Z2)/sqrt(x**2 + y**2 + ac**2)
    return
end function Vroc

real(num) function DVroc(Z1, Z2, x, y, ac, cmp)
    implicit none
    real(num), intent(in) :: Z1, Z2, x, y, ac
    integer, intent(in) :: cmp

    select case (cmp)
    case (1)
        DVroc = -x*(Z1*Z2)/((x**2 + y**2 + ac**2)**(3.0/2.0_num))
    case (2)
        DVroc = -y*(Z1*Z2)/((x**2 + y**2 + ac**2)**(3.0/2.0_num))
    end select

    return
end function DVroc

end module rochester
