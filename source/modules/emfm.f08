! Electromagnetic Field Module (emfm.f08)
!
! Author:  D. Younis
!          University of Rochester
!          Department of Physics
!
! Written: 6/2/2020
! Revised: 8/23/2022

module emfm

use prec
use math

implicit none

type, public :: emf
    real(num) :: E0, omg0, eps, ch1, Ncyc_rf, Ncyc_pl, CEP, t_on, t_off, T0, Tp, Tpk, Tfwhm
    integer :: it_on, it_off
    character(len=:), allocatable :: profile
    real(num), dimension(:), allocatable :: Ex, Ey
    real(num), dimension(:), allocatable :: Ax, Ay
    real(num), dimension(:), allocatable :: Cx, Cy
contains
    procedure :: init_trapz => emf_trapezoidal_pulse
    procedure :: init_sine2 => emf_sine_squared_pulse
    procedure :: init_gauss => emf_gaussian_pulse
    procedure :: init_gauss_l => emf_gaussian_l_pulse
end type emf

contains

subroutine emf_trapezoidal_pulse(this, t)
    implicit none
    class(emf), intent(inout) :: this
    real(num), intent(in) :: t(:)

    real(num), allocatable :: ft(:)
    real(num) :: dt
    integer :: nt, it, idx(4)

    dt = t(2) - t(1); nt = size(t);
    allocate(this%Ex(nt), this%Ey(nt))
    allocate(this%Ax(nt), this%Ay(nt))
    allocate(this%Cx(nt), this%Cy(nt))

    read(this%profile(7:),*) this%Ncyc_rf ! no. cycles rise/fall
    this%T0 = (2.*pi)/this%omg0 ! period
    this%Tp = (2*this%Ncyc_rf + this%Ncyc_pl)*this%T0 ! total duration

    ! set time indices
    idx(1) = nint(this%t_on/dt) + 1 ! pulse start
    idx(2) = nint((this%t_on + this%Ncyc_rf*this%T0)/dt) + 1 ! plateau start
    idx(3) = nint((this%t_on + this%Ncyc_rf*this%T0 + this%Ncyc_pl*this%T0)/dt) + 1 ! plateau end
    idx(4) = nint((this%t_on + this%Tp)/dt) + 1 ! pulse end

    ! construct trapezoidal function with Ncyc_rf-cycle turn-on/turn-off
    allocate(ft(nt)); ft = 0.0_num;
    do it=1,nt
        if ((it >= idx(1)).and.(it < idx(2))) then
            ft(it) = (t(it) - this%t_on)/(this%Ncyc_rf*this%T0)
        else if ((it >= idx(2)).and.(it <= idx(3))) then
            ft(it) = 1.0_num
        else if ((it > idx(3)).and.(it <= idx(4))) then
            ft(it) = (this%Tp + this%t_on - t(it))/(this%Ncyc_rf*this%T0)
        end if
    end do

    ! construct laser electric field
    this%Ex = this%E0/sqrt(1.0+this%eps**2) * ft &
        * sin(this%omg0*(t-this%t_on) + this%ch1*(t-this%t_on)**2 + this%CEP)
    this%Ey = this%E0/sqrt(1.0+this%eps**2) * ft &
        * this%eps*cos(this%omg0*(t-this%t_on) + this%ch1*(t-this%t_on)**2 + this%CEP)

    ! compute vector potential
    this%Ax = -1.0 * simint(this%Ex, 0.0_num, dt)
    this%Ay = -1.0 * simint(this%Ey, 0.0_num, dt)

    ! compute excursion
    this%Cx = simint(this%Ax, 0.0_num, dt)
    this%Cy = simint(this%Ay, 0.0_num, dt)

    ! set pulse start/end times
    this%t_off = this%t_on + this%Tp
    this%it_on = nint(this%t_on/dt) + 1
    this%it_off = nint(this%t_off/dt) + 1

    ! ! check that the pulse ends before the simulation (annoying)
    ! if (this%t_off > t(nt)) then
    !     print '(A)', 'Error: Pulse duration exceeds simulation time.'
    !     print '(A)', 'Terminating execution.'
    !     call exit(1)
    ! end if

    return
end subroutine emf_trapezoidal_pulse

subroutine emf_sine_squared_pulse(this, t)
    implicit none
    class(emf), intent(inout) :: this
    real(num), intent(in) :: t(:)

    real(num), allocatable :: ft(:)
    real(num) :: dt
    integer :: nt, it

    dt = t(2) - t(1); nt = size(t);
    allocate(this%Ex(nt), this%Ey(nt))
    allocate(this%Ax(nt), this%Ay(nt))
    allocate(this%Cx(nt), this%Cy(nt))

    this%T0 = (2.*pi)/this%omg0 ! period

    ! set pulse start/end times
    this%t_off = this%t_on + this%Tp
    this%it_on = nint(this%t_on/dt) + 1
    this%it_off = nint(this%t_off/dt) + 1

    ! construct sine-squared function
    allocate(ft(nt)); ft = 0.0_num;
    do it=1,nt
        if ((t(it) >= this%t_on).and.(t(it) <= this%t_off)) &
        ft(it) = sin(pi*(t(it)-this%t_on)/this%Tp)**2
    end do

    ! construct laser vector potential
    this%Ax = this%E0/this%omg0/sqrt(1.0+this%eps**2) * ft &
        * sin(this%omg0*(t-this%t_on) + this%ch1*(t-this%t_on)**2 + this%CEP)
    this%Ay = this%E0/this%omg0/sqrt(1.0+this%eps**2) * ft &
        * this%eps*cos(this%omg0*(t-this%t_on) + this%ch1*(t-this%t_on)**2 + this%CEP)

    ! compute electric field
    call deriv(-this%Ax, dt, this%Ex)
    call deriv(-this%Ay, dt, this%Ey)

    ! compute excursion
    this%Cx = simint(this%Ax, 0.0_num, dt)
    this%Cy = simint(this%Ay, 0.0_num, dt)

    ! ! check that the pulse ends before the simulation (annoying)
    ! if (this%t_off > t(nt)) then
    !     print '(A)', 'Error: Pulse duration exceeds simulation time.'
    !     print '(A)', 'Terminating execution.'
    !     call exit(1)
    ! end if

    return
end subroutine emf_sine_squared_pulse

subroutine emf_gaussian_pulse(this, t)
    implicit none
    class(emf), intent(inout) :: this
    real(num), intent(in) :: t(:)

    real(num), allocatable :: ft(:)
    real(num) :: dt, w0t
    integer :: nt, idx(2)

    dt = t(2) - t(1); nt = size(t);
    allocate(this%Ex(nt), this%Ey(nt))
    allocate(this%Ax(nt), this%Ay(nt))
    allocate(this%Cx(nt), this%Cy(nt))

    w0t = this%Tfwhm/sqrt(2.0*log(2.0)) ! 1/e^2 duration from FWHM
    this%T0 = (2.*pi)/this%omg0 ! period

    ! construct gaussian function
    allocate(ft(nt)); ft = exp(-(t-this%Tpk)**2/w0t**2);

    ! set time indices, start/end ~ time_center +/- 3*pulse_radius (3*w0t)
    ! approximate since gaussian is of infinite extent
    idx(1) = nint((this%Tpk - 3.0*w0t)/dt) + 1
    idx(2) = nint((this%Tpk + 3.0*w0t)/dt) + 1

    ! set pulse start/end times
    this%it_on = 1
    this%it_off = nt
    if (idx(1) > this%it_on) this%it_on = idx(1)
    if (idx(2) < this%it_off) this%it_off = idx(2)
    this%t_on = t(this%it_on)
    this%t_off = t(this%it_off)
    this%Tp = this%t_off - this%t_on

    ! construct laser electric field
    this%Ex = this%E0/sqrt(1.0+this%eps**2) * ft &
        * sin(this%omg0*(t-this%t_on) + this%ch1*(t-this%t_on)**2 + this%CEP)
    this%Ey = this%E0/sqrt(1.0+this%eps**2) * ft &
        * this%eps*cos(this%omg0*(t-this%t_on) + this%ch1*(t-this%t_on)**2 + this%CEP)

    ! compute vector potential
    this%Ax = -1.0 * simint(this%Ex, 0.0_num, dt)
    this%Ay = -1.0 * simint(this%Ey, 0.0_num, dt)

    ! compute excursion
    this%Cx = simint(this%Ax, 0.0_num, dt)
    this%Cy = simint(this%Ay, 0.0_num, dt)

    return
end subroutine emf_gaussian_pulse

subroutine emf_gaussian_l_pulse(this, t)
    implicit none
    class(emf), intent(inout) :: this
    real(num), intent(in) :: t(:)

    real(num), allocatable :: ft(:)
    real(num) :: dt, w0t
    integer :: nt, it, idx(4)

    dt = t(2) - t(1); nt = size(t);
    allocate(this%Ex(nt), this%Ey(nt))
    allocate(this%Ax(nt), this%Ay(nt))
    allocate(this%Cx(nt), this%Cy(nt))

    w0t = this%Tfwhm/sqrt(2.0*log(2.0)) ! 1/e^2 duration from FWHM
    this%T0 = (2.*pi)/this%omg0 ! period

    ! construct gaussian function
    allocate(ft(nt)); ft = exp(-(t-this%Tpk)**2/w0t**2);

    ! set time indices
    idx(1) = nint((this%Tpk - 3.0*w0t)/dt) + 1 ! rise start
    idx(2) = nint((this%Tpk - 2.0*w0t)/dt) + 1 ! rise end
    idx(3) = nint((this%Tpk + 2.0*w0t)/dt) + 1 ! fall start
    idx(4) = nint((this%Tpk + 3.0*w0t)/dt) + 1 ! fall end

    ! check that the pulse fits in the simulation time domain
    if ((idx(1) < 1).or.(idx(4) > nt)) then
        print '(A)', 'Error: Gaussian pulse does not fit in simulation time domain.'
        print '(A,F10.2,A,F10.2,A)', 'Required window: [TIME_CENTER +/- 3*PULSE_RADIUS] = [', &
            this%Tpk-3.0*w0t, ', ', this%Tpk+3.0*w0t, '] au'
        print '(A)', 'Terminating execution.'
        call exit(1)
    end if

    ft(1:idx(1)) = 0.0_num
    ft(idx(4):nt) = 0.0_num

    ! linear cutoffs between [2,3] pulse radius (w0t)
    do it=idx(1),idx(2)
        ft(it) = (t(it) - t(idx(1)))/(t(idx(2)) - t(idx(1))) * ft(idx(2))
    end do
    do it=idx(3),idx(4)
        ft(it) = (t(it) - t(idx(4)))/(t(idx(3)) - t(idx(4))) * ft(idx(3))
    end do

    ! set pulse start/end times
    this%it_on = idx(1)
    this%it_off = idx(4)
    this%t_on = t(this%it_on)
    this%t_off = t(this%it_off)
    this%Tp = this%t_off - this%t_on

    ! construct laser electric field
    this%Ex = this%E0/sqrt(1.0+this%eps**2) * ft &
        * sin(this%omg0*(t-this%t_on) + this%ch1*(t-this%t_on)**2 + this%CEP)
    this%Ey = this%E0/sqrt(1.0+this%eps**2) * ft &
        * this%eps*cos(this%omg0*(t-this%t_on) + this%ch1*(t-this%t_on)**2 + this%CEP)

    ! compute vector potential
    this%Ax = -1.0 * simint(this%Ex, 0.0_num, dt)
    this%Ay = -1.0 * simint(this%Ey, 0.0_num, dt)

    ! compute excursion
    this%Cx = simint(this%Ax, 0.0_num, dt)
    this%Cy = simint(this%Ay, 0.0_num, dt)

    return
end subroutine emf_gaussian_l_pulse

end module emfm
