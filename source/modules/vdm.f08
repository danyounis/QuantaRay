! Virtual Detector Module (vdm.f08)
!
! Author:  D. Younis
!          University of Rochester
!          Department of Physics
!
! Written: 6/2/2020
! Revised: 8/23/2022
!
! References:
! [1] B. Feuerstein and U. Thumm, J. Phys. B 36, 707-716 (2003).
! [2] X. Wang, J. Tian, and J. H. Eberly, J. Phys. B 51, 084002 (2018).

module vdm

use prec
use math
use quantum

implicit none

type, public :: vdet
    real(num) :: xl, yl
    integer :: ixl, iyl
    real(num), dimension(:), allocatable :: rho, phase
    real(num), dimension(:,:), allocatable :: Krt, Jrt
contains
    procedure :: init => vdet_initialize
    generic, public :: trigger => vdet_calculate_current_1d, vdet_calculate_current_2d
    procedure, private :: vdet_calculate_current_1d, vdet_calculate_current_2d
end type vdet

type, public :: edet
    integer :: nde
    real(num), dimension(2) :: bfwt
    real(num), dimension(:,:), allocatable :: data
contains
    procedure :: escan => edet_detect
end type edet

type, public :: particle_electron
    real(num) :: x, y
    real(num) :: px, py
    real(num) :: weight, phase
    integer :: ix, iy
    logical :: propagate
contains
    generic, public :: apush => electron_dt_propagate_analytic_1d, electron_dt_propagate_analytic_2d
    procedure, private :: electron_dt_propagate_analytic_1d, electron_dt_propagate_analytic_2d
    procedure :: npush => electron_dt_propagate_numeric
end type particle_electron

type, public :: grid_electron
    real(num), dimension(:), allocatable :: x, y
    real, dimension(2) :: x_lim, y_lim
    real(num) :: dr(2)
    integer :: nr(2)
end type grid_electron

contains

subroutine vdet_initialize(this, R0, Nv, x, y, n, nt, sdim)
    implicit none
    class(vdet), intent(inout) :: this

    real(num), intent(in) :: x(:), y(:)
    real, intent(in) :: R0
    integer, intent(in) :: Nv, n, nt, sdim

    select case (sdim)
    case (2)
        ! set virtual detector (VD) position
        this%xl = R0*cos(2.*pi*(n-1)/real(Nv,num))
        this%yl = R0*sin(2.*pi*(n-1)/real(Nv,num))

        ! get index of nearest upper-left grid point
        this%ixl = minloc(abs(x-this%xl),1)
        this%iyl = minloc(abs(y-this%yl),1)
        if (x(this%ixl) >= this%xl) this%ixl = this%ixl - 1
        if (y(this%iyl) >= this%yl) this%iyl = this%iyl - 1
    case (1)
        this%xl = R0*cos(2.*pi*(n-1)/real(Nv,num))
        this%yl = 0.0_num

        this%ixl = minloc(abs(x-this%xl),1)
        this%iyl = 1
        if (x(this%ixl) >= this%xl) this%ixl = this%ixl - 1
    end select

    ! initialize blank data arrays
    allocate(this%Krt(sdim,nt), this%Jrt(sdim,nt))
    allocate(this%rho(nt), this%phase(nt))
    this%Krt = 0.0_num; this%Jrt = 0.0_num;
    this%rho = 0.0_num; this%phase = 0.0_num;

    return
end subroutine vdet_initialize

subroutine vdet_calculate_current_1d(this, wavefn, x, k)
    implicit none
    class(vdet), intent(inout) :: this
    class(SchrodingerWavefunction1D), intent(in) :: wavefn
    real(num), intent(in) :: x(:)
    integer, intent(in) :: k

    ! temporary interpolation variables
    real(num) :: xa, xb
    complex(num) :: psi, grad_psi

    ! linearly-interpolate wavefunction to VD position
    xa = x(this%ixl); xb = x(this%ixl+1);

    psi = ((this%xl-xa)*wavefn%psi(this%ixl+1) &
        + (xb-this%xl)*wavefn%psi(this%ixl))/(xb-xa)

    grad_psi = ((this%xl-xa)*wavefn%grad_psi(this%ixl+1) &
        + (xb-this%xl)*wavefn%grad_psi(this%ixl))/(xb-xa)

    ! PROBABILITY CURRENT !
    this%Jrt(1,k) = real((i/2.0)*(psi*conjg(grad_psi) - conjg(psi)*grad_psi), num)

    ! MOMENTUM !
    this%rho(k) = abs(psi)**2
    this%Krt(1,k) = this%Jrt(1,k)/this%rho(k)

    ! PHASE !
    this%phase(k) = atan2(aimag(psi),real(psi))

    return
end subroutine vdet_calculate_current_1d

subroutine vdet_calculate_current_2d(this, wavefn, IDM, x, y, k)
    implicit none
    class(vdet), intent(inout) :: this
    class(SchrodingerWavefunction2D), intent(in) :: wavefn

    integer, intent(in) :: k
    real(num), intent(in) :: x(:), y(:), IDM(:,:,:)

    ! temporary interpolation variables
    complex(num) :: psi, psix1, psix2
    real(num), dimension(4) :: by, by1, by2, by12
    real(num) :: ansy, ansy1, ansy2

    ! indexing variables
    integer :: j, l, lp, c(4,2)

    ! grid corners surrounding VD position
    c(1,:) = (/this%ixl+1, this%iyl/)
    c(2,:) = (/this%ixl+1, this%iyl+1/)
    c(3,:) = (/this%ixl, this%iyl+1/)
    c(4,:) = (/this%ixl, this%iyl/)

    ! PROBABILITY CURRENT !
    do j=1,2
        ! set IDM lookup index
        if (j == 1) lp = 5
        if (j == 2) lp = 9

        ! get IDM data for J1/J2
        do l=1,4
            by(l)   = IDM(c(l,1), c(l,2), lp)
            by1(l)  = IDM(c(l,1), c(l,2), lp+1)
            by2(l)  = IDM(c(l,1), c(l,2), lp+2)
            by12(l) = IDM(c(l,1), c(l,2), lp+3)
        end do

        ! interpolate J1/J2 to VD position
        call bcuint(by, by1, by2, by12, x(this%ixl), y(this%iyl), &
            this%xl, this%yl, (/x(2)-x(1),y(2)-y(1)/), ansy, ansy1, ansy2)

        ! record probability current
        this%Jrt(j,k) = ansy
    end do

    ! MOMENTUM !
    ! get IDM data for rho
    do l=1,4
        by(l)   = IDM(c(l,1), c(l,2), 1)
        by1(l)  = IDM(c(l,1), c(l,2), 2)
        by2(l)  = IDM(c(l,1), c(l,2), 3)
        by12(l) = IDM(c(l,1), c(l,2), 4)
    end do

    ! interpolate rho to VD position
    call bcuint(by, by1, by2, by12, x(this%ixl), y(this%iyl), &
        this%xl, this%yl, (/x(2)-x(1),y(2)-y(1)/), ansy, ansy1, ansy2)
    this%rho(k) = ansy

    ! record momentum
    do j=1,2
        this%Krt(j,k) = this%Jrt(j,k)/this%rho(k)
    end do

    ! PHASE !
    ! bilinear interpolation of wavefunction to VD position

    psix1 = ((x(this%ixl+1)-this%xl)*wavefn%psi(c(4,1),c(4,2)) &
        + (this%xl-x(this%ixl))*wavefn%psi(c(1,1),c(1,2)))/(x(this%ixl+1)-x(this%ixl))

    psix2 = ((x(this%ixl+1)-this%xl)*wavefn%psi(c(3,1),c(3,2)) &
        + (this%xl-x(this%ixl))*wavefn%psi(c(2,1),c(2,2)))/(x(this%ixl+1)-x(this%ixl))

    psi = ((y(this%iyl+1)-this%yl)*psix1 + (this%yl-y(this%iyl))*psix2)/(y(this%iyl+1)-y(this%iyl))

    this%phase(k) = atan2(aimag(psi),real(psi))

    return
end subroutine vdet_calculate_current_2d

subroutine vdet_calculate_current_from_phase(this, psi, dr, nr, k)
    implicit none
    class(vdet),  intent(inout) :: this
    integer,      intent(in)    :: nr(2), k
    real(num),    intent(in)    :: dr(2)
    complex(num), intent(in)    :: psi(nr(1),nr(2))

    real(num) :: phase_psi(nr(1),nr(2)), grad_phase_psi(nr(1),nr(2))

    integer :: j

    ! compute wavefunction phase
    phase_psi = atan2(aimag(psi),real(psi))
    call phase_unwrap(phase_psi)

    do j=1,2
        ! wavefunction gradient along axis j
        call deriv(phase_psi,dr(j),j,grad_phase_psi)
        ! momentum
        this % Krt(j,k) = grad_phase_psi(this%ixl,this%iyl)
        ! probability current
        this % Jrt(j,k) = abs(psi(this%ixl,this%iyl))**2 * this%Krt(j,k)
    end do

    return
end subroutine vdet_calculate_current_from_phase

subroutine edet_detect(this, electron, sdim)
    implicit none
    class(edet),              intent(inout) :: this
    class(particle_electron), intent(in)    :: electron(:)
    integer,                  intent(in)    :: sdim

    integer :: ine

    ! total number of electrons
    this % nde = size(electron,1)

    select case (sdim)
    case (1) ! 1D trajectories
        allocate(this % data(this%nde,4))
        do ine=1,this%nde
            this % data(ine,:) = (/ &
                electron(ine) % x,  &
                electron(ine) % px, &
                electron(ine) % phase, &
                electron(ine) % weight &
                /)
        end do

    case (2) ! 2D trajectories
        allocate(this % data(this%nde,6))
        do ine=1,this%nde
            this % data(ine,:) = (/ &
                electron(ine) % x,  &
                electron(ine) % y,  &
                electron(ine) % px, &
                electron(ine) % py, &
                electron(ine) % phase, &
                electron(ine) % weight &
                /)
        end do
    end select

    return
end subroutine edet_detect

subroutine electron_dt_propagate_analytic_1d(this, pdot, dt)
    implicit none
    class(particle_electron), intent(inout) :: this
    real(num), intent(in) :: dt, pdot

    ! check propagation flag
    ! if (.not.this%propagate) return

    ! advance momentum
    this%px = this%px + pdot*dt

    ! advance position
    this%x = this%x + (this%px)*dt

    return
end subroutine electron_dt_propagate_analytic_1d

subroutine electron_dt_propagate_analytic_2d(this, pdot, dt)
    implicit none
    class(particle_electron), intent(inout) :: this
    real(num), intent(in) :: dt, pdot(2)

    ! check propagation flag
    ! if (.not.this%propagate) return

    ! advance momenta
    this%px = this%px + pdot(1)*dt
    this%py = this%py + pdot(2)*dt

    ! advance positions
    this%x = this%x + (this%px)*dt
    this%y = this%y + (this%py)*dt

    return
end subroutine electron_dt_propagate_analytic_2d

subroutine electron_dt_propagate_numeric(this, V, x, y, dr, dt, nr)
    implicit none
    class(particle_electron), intent(inout) :: this
    integer,   intent(in) :: nr(2)
    real(num), intent(in) :: x(nr(1)), y(nr(2)), dr(2), dt

    real(num), intent(in) :: V(nr(1),nr(2))
    real(num) :: grad_V(nr(1),nr(2))

    ! check propagation flag
    if (.not.this%propagate) return

    ! advance momenta
    call deriv(V,dr(1),1,grad_V)
    this%px = this%px - grad_V(this%ix,this%iy)*dt
    call deriv(V,dr(2),2,grad_V)
    this%py = this%py - grad_V(this%ix,this%iy)*dt

    ! advance positions
    this%x = this%x + (this%px)*dt
    this%y = this%y + (this%py)*dt

    ! update position indices
    this%ix = nint((this%x-x(1))/dr(1)) + 1
    this%iy = nint((this%y-y(1))/dr(2)) + 1

    ! disable pushing if electron leaves simulation grid
    if ((this%ix <= 2).or.(this%ix >= nr(1)-1)) this%propagate = .false.
    if ((this%iy <= 2).or.(this%iy >= nr(2)-1)) this%propagate = .false.

    return
end subroutine electron_dt_propagate_numeric

end module vdm
