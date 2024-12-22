! Quantum Mechanics Module (quantum.f08)
!
! Author:  D. Younis
!          University of Rochester
!          Department of Physics
!
! Written: 4/30/2020
! Revised: 12/22/2024

module quantum

use prec
use math
use omp_lib

implicit none

! __________________________________________________________________________________________________
!
! SCHRODINGER WAVEFUNCTION TYPE-STRUCTURES
! __________________________________________________________________________________________________
!

type, public :: SchrodingerWavefunction2D
    complex(num), dimension(:,:),   allocatable :: psi
    complex(num), dimension(:,:,:), allocatable :: grad_psi
    real(num),    dimension(:,:),   allocatable :: phase
    real(num),    dimension(:),     allocatable :: norm, energy
    real(num),    dimension(:),     allocatable :: D2x, D2y, M2x, M2y
    integer,      dimension(:,:),   allocatable :: stabx, staby
contains
    procedure :: init_vars => psi2d_initialize_vars
    procedure :: init_form => psi2d_initialize_form
    procedure :: make_cnn_mats => psi2d_make_cnn_matrices
    procedure :: propagate_fft => psi2d_dt_propagate_fft
    procedure :: propagate_cnn => psi2d_dt_propagate_cnn
    procedure :: destroy => psi2d_destructor
end type SchrodingerWavefunction2D

type, public :: SchrodingerWavefunction1D
    complex(num), dimension(:), allocatable :: psi, grad_psi
    real(num),    dimension(:), allocatable :: phase, norm, energy
    real(num),    dimension(:), allocatable :: D2x, M2x
    integer,    dimension(:,:), allocatable :: stab
contains
    procedure :: init_vars => psi1d_initialize_vars
    procedure :: init_form => psi1d_initialize_form
    procedure :: make_cnn_mats => psi1d_make_cnn_matrices
    procedure :: propagate_fft => psi1d_dt_propagate_fft
    procedure :: propagate_cnn => psi1d_dt_propagate_cnn
    procedure :: destroy => psi1d_destructor
end type SchrodingerWavefunction1D

type, public :: SchrodingerWavefunction1DR
    complex(num), dimension(:,:), allocatable :: phi
    real(num),    dimension(:,:), allocatable :: V0, dV0_dr
    complex(num), dimension(:),   allocatable :: Va
    real(num),    dimension(:),   allocatable :: norm, energy, clm
    real(num),    dimension(:),   allocatable :: D1r, D2r, M1r, M2r
    integer,      dimension(:,:), allocatable :: stab
    real(num) :: Z, D211, M211
    integer   :: m
contains
    procedure :: init_vars => psi1dr_initialize_vars
    procedure :: init_prop => psi1dr_initialize_propagators
    procedure :: init_form => psi1dr_initialize_form
    procedure :: prep_atom => psi1dr_prepare_atomic_state
    procedure :: destroy => psi1dr_destructor
    generic, public :: propagate => psi1dr_dt_propagate_free, psi1dr_dt_propagate_full
    procedure, private :: psi1dr_dt_propagate_free, psi1dr_dt_propagate_full
end type SchrodingerWavefunction1DR

! __________________________________________________________________________________________________
!
! MISCELLANEOUS TYPE-STRUCTURES
! __________________________________________________________________________________________________
!

type, public :: tSURFF2D
    complex(num), dimension(:,:), allocatable :: p_dist
    real(num),    dimension(:),   allocatable :: kx, ky
    real(num),    dimension(:),   allocatable :: xl, yl
    integer,      dimension(:),   allocatable :: ixl, iyl
    character(len=2) :: interp
    real(num) :: R0, dti, dphi, kx_lim(2), ky_lim(2)
    integer   :: iti, Ns, nk(2)
    logical   :: enable
contains
    procedure :: init => tSURFF2D_initialize
    procedure :: dt_step => tSURFF2D_dt_step
    procedure :: destroy => tSURFF2D_destructor
end type tSURFF2D

type, public :: pconst_mks
    real(dble_t) :: c = 2.99792458d+8 ! speed of light (m/s, exact)
    real(dble_t) :: e = 1.602176634d-19 ! elementary charge (C, exact)
    real(dble_t) :: m_e = 9.1093837015d-31 ! electron mass (kg, +/- 2.8d-40 kg)
    real(dble_t) :: hbar = 1.0545718176461565d-34 ! reduced Planck constant (J.s, exact)
    real(dble_t) :: alpha = 7.2973525693d-3 ! fine-structure constant (+/- 1.1d-12)
    real(dble_t) :: lambda_Compton = 2.42631023867d-12 ! Compton wavelength (m, +/- 7.3d-22 m)
end type pconst_mks

type, public :: pconst_cgs
    real(dble_t) :: c = 2.99792458d+10 ! speed of light (cm/s, exact)
    real(dble_t) :: e = 4.803204712570263e-10 ! elementary charge (statC, exact)
    real(dble_t) :: m_e = 9.1093837015d-28 ! electron mass (g, +/- 2.8d-37 g)
    real(dble_t) :: hbar = 1.0545718176461565d-27 ! reduced Planck constant (erg.s, exact)
    real(dble_t) :: alpha = 7.2973525693d-3 ! fine-structure constant (+/- 1.1d-12)
    real(dble_t) :: lambda_Compton = 2.42631023867d-10 ! Compton wavelength (cm, +/- 7.3d-20 cm)
end type pconst_cgs

! __________________________________________________________________________________________________
!
! ENERGY/RADIATION CALCULATION PROCEDURES
! __________________________________________________________________________________________________
!

interface calc_energy
    module procedure expectE_1D_fft
    module procedure expectE_1D_cnn
    module procedure expectE_1DR
    module procedure expectE_2D_fft
    module procedure expectE_2D_cnn
end interface calc_energy

interface calc_spectrum
    module procedure photoe_spectrum_winop_1D
    module procedure photoe_spectrum_winop_1DR
end interface calc_spectrum

interface calc_radintens
    module procedure radiative_intensity_1DR
    module procedure radiative_intensity_2D1e
    module procedure radiative_intensity_1D2e
end interface calc_radintens

contains

! __________________________________________________________________________________________________
!
! WAVEFUNCTION INITIALIZATION ROUTINES
! __________________________________________________________________________________________________
!

subroutine psi1d_initialize_vars(this, nx, nt)
    implicit none
    class(SchrodingerWavefunction1D), intent(inout) :: this
    integer, intent(in) :: nx, nt

    allocate(this%psi(nx))
    allocate(this%grad_psi(nx), this%phase(nx))
    allocate(this%norm(nt), this%energy(nt))

    return
end subroutine psi1d_initialize_vars

subroutine psi1dr_initialize_vars(this, nr, nt, l_max)
    implicit none
    class(SchrodingerWavefunction1DR), intent(inout) :: this
    integer, intent(in) :: nr, nt, l_max

    allocate(this%phi(l_max+1,nr), this%Va(nr))
    allocate(this%norm(nt), this%energy(nt))

    this%Va = (0.0_num,0.0_num)

    return
end subroutine psi1dr_initialize_vars

subroutine psi2d_initialize_vars(this, nr, nt)
    implicit none
    class(SchrodingerWavefunction2D), intent(inout) :: this
    integer, intent(in) :: nr(2), nt

    allocate(this%psi(nr(1),nr(2)))
    allocate(this%grad_psi(nr(1),nr(2),2), this%phase(nr(1),nr(2)))
    allocate(this%norm(nt), this%energy(nt))

    return
end subroutine psi2d_initialize_vars

subroutine psi1d_initialize_form(this, parity, dx)
    implicit none
    class (SchrodingerWavefunction1D), intent(inout) :: this
    integer, intent(in) :: parity
    real(num), intent(in) :: dx

    real(num), dimension(:), allocatable :: Re_psi, Im_psi
    integer :: nx

    nx = size(this%psi)
    allocate(Re_psi(nx), Im_psi(nx))

    ! randomize wavefunction components
    call init_RNG()
    call random_number(Re_psi)
    call random_number(Im_psi)
    this%psi = Re_psi + i*Im_psi

    ! impose parity condition
    if (parity /= 0) then
        this%psi = this%psi + sign(1,parity)*fliplr(this%psi)
    end if

    ! normalize
    this%psi = this%psi/sqrt(dx*trapz(abs(this%psi)**2))

    return
end subroutine psi1d_initialize_form

subroutine psi1dr_initialize_form(this, dr)
    implicit none
    class (SchrodingerWavefunction1DR), intent(inout) :: this
    real(num), intent(in) :: dr

    real(num), dimension(:), allocatable :: Re_phi, Im_phi
    integer :: l, l_max, nr

    l_max = size(this%phi,1) - 1
    nr = size(this%phi,2)

    allocate(Re_phi(nr), Im_phi(nr))

    ! randomize wavefunction components
    call init_RNG()
    !$OMP PARALLEL DO SHARED(this) PRIVATE(Re_phi, Im_phi)
    do l=0,l_max
        call random_number(Re_phi)
        call random_number(Im_phi)
        this%phi(l+1,:) = Re_phi + i*Im_phi

        ! normalize
        this%phi(l+1,:) = this%phi(l+1,:)/sqrt(dr*trapz(abs(this%phi(l+1,:))**2))
    end do
    !$OMP END PARALLEL DO

    return
end subroutine psi1dr_initialize_form

subroutine psi2d_initialize_form(this, parity, dr)
    implicit none
    class (SchrodingerWavefunction2D), intent(inout) :: this
    integer, intent(in) :: parity
    real(num), intent(in) :: dr(2)

    real(num), dimension(:,:), allocatable :: Re_psi, Im_psi
    integer :: nr(2)

    nr(1) = size(this%psi,1); nr(2) = size(this%psi,2);
    allocate(Re_psi(nr(1),nr(2)), Im_psi(nr(1),nr(2)))

    ! randomize wavefunction components
    call init_RNG()
    call random_number(Re_psi)
    call random_number(Im_psi)
    this%psi = Re_psi + i*Im_psi

    ! impose parity condition
    select case (parity)
    case default ! indefinite parity
        continue
    case (-2) ! anti-symmetric under pi rotation
        this%psi = this%psi + sign(1,parity)*fliplr(this%psi,ax=0)
    case (-1,1) ! anti-symmetric under pi/2 rotation or centro-symmetric
        this%psi = this%psi + sign(1,parity)*fliplr(this%psi,ax=1)
        this%psi = this%psi + sign(1,parity)*fliplr(this%psi,ax=2)
    end select

    ! normalize
    this%psi = this%psi/sqrt(dr(1)*dr(2)*trapz(abs(this%psi)**2))

    return
end subroutine psi2d_initialize_form

subroutine psi1d_make_cnn_matrices(this, nx, dx)
    implicit none
    class(SchrodingerWavefunction1D), intent(inout) :: this
    real(num), intent(in) :: dx
    integer, intent(in) :: nx
    integer :: j, k

    ! D2x/M2x are tri-diagonal
    ! nx diag + 2*(nx-1) off-diag elts = 3*nx-2 total elts
    ! row-major order dense-to-sparse convention
    allocate(this%D2x(3*nx-2), this%M2x(3*nx-2))

    ! sparse tri-diagonal element sequences:
    !   0 diag: 1, 3*nx-2, 3
    !  +1 diag: 2, 3*nx-4, 3
    !  -1 diag: 3, 3*nx-3, 3

    ! 2nd derivative centered-difference matrix
    this%D2x = 1.0
    do j=1, 3*nx-2, 3
        this%D2x(j) = -2.0
    end do
    this%D2x = this%D2x/(dx**2)

    ! Numerov matrix
    this%M2x = (dx**2/12.0)*this%D2x
    do j=1, 3*nx-2, 3
        this%M2x(j) = 1.0 + this%M2x(j)
    end do

    ! set the populated indices of sparse tri-diagonal matrices
    allocate(this%stab(3*nx-2,2))
    ! key: stab(i,:) = (j,k) is the i-th tri-diagonal element at position (j,k)
    this%stab(1,:) = (/1,1/)
    this%stab(2,:) = (/1,2/)
    k = 2
    do j=3, 3*nx-6, 3
        this%stab(j,:)   = (/k, k-1/)
        this%stab(j+1,:) = (/k, k  /)
        this%stab(j+2,:) = (/k, k+1/)
        k = k + 1
    end do
    this%stab(3*nx-3,:) = (/nx, nx-1/)
    this%stab(3*nx-2,:) = (/nx, nx  /)

    return
end subroutine psi1d_make_cnn_matrices

subroutine psi2d_make_cnn_matrices(this, nr, dr)
    implicit none
    class(SchrodingerWavefunction2D), intent(inout) :: this
    real(num), intent(in) :: dr(2)
    integer, intent(in) :: nr(2)
    integer :: j, k

    ! D2(x,y) and M2(x,y) are tri-diagonal
    ! n diag + 2*(n-1) off-diag elts = 3*n-2 total elts
    ! row-major order dense-to-sparse convention
    allocate(this%D2x(3*nr(1)-2), this%M2x(3*nr(1)-2))
    allocate(this%D2y(3*nr(2)-2), this%M2y(3*nr(2)-2))

    ! sparse tri-diagonal element sequences:
    !   0 diag: 1, 3*n-2, 3
    !  +1 diag: 2, 3*n-4, 3
    !  -1 diag: 3, 3*n-3, 3

    ! 2nd derivative centered-difference matrix
    this%D2x = 1.0
    do j=1, 3*nr(1)-2, 3
        this%D2x(j) = -2.0
    end do
    this%D2x = this%D2x/(dr(1)**2)

    this%D2y = 1.0
    do j=1, 3*nr(2)-2, 3
        this%D2y(j) = -2.0
    end do
    this%D2y = this%D2y/(dr(2)**2)

    ! Numerov matrix
    this%M2x = (dr(1)**2/12.0)*this%D2x
    do j=1, 3*nr(1)-2, 3
        this%M2x(j) = 1.0 + this%M2x(j)
    end do

    this%M2y = (dr(2)**2/12.0)*this%D2y
    do j=1, 3*nr(2)-2, 3
        this%M2y(j) = 1.0 + this%M2y(j)
    end do

    ! set the populated indices of sparse tri-diagonal matrices
    allocate(this%stabx(3*nr(1)-2,2), this%staby(3*nr(2)-2,2))
    ! key: stab(i,:) = (j,k) is the i-th tri-diagonal element at position (j,k)
    this%stabx(1,:) = (/1,1/); this%staby(1,:) = (/1,1/);
    this%stabx(2,:) = (/1,2/); this%staby(2,:) = (/1,2/);
    k = 2
    do j=3, 3*nr(1)-6, 3
        this%stabx(j,:)   = (/k, k-1/)
        this%stabx(j+1,:) = (/k, k  /)
        this%stabx(j+2,:) = (/k, k+1/)
        k = k + 1
    end do
    k = 2
    do j=3, 3*nr(2)-6, 3
        this%staby(j,:)   = (/k, k-1/)
        this%staby(j+1,:) = (/k, k  /)
        this%staby(j+2,:) = (/k, k+1/)
        k = k + 1
    end do
    this%stabx(3*nr(1)-3,:) = (/nr(1), nr(1)-1/)
    this%stabx(3*nr(1)-2,:) = (/nr(1), nr(1)  /)
    this%staby(3*nr(2)-3,:) = (/nr(2), nr(2)-1/)
    this%staby(3*nr(2)-2,:) = (/nr(2), nr(2)  /)

    return
end subroutine psi2d_make_cnn_matrices

subroutine psi1dr_initialize_propagators(this, r, dr)
    implicit none
    class(SchrodingerWavefunction1DR), intent(inout) :: this
    real(num), intent(in) :: r(:), dr
    integer :: l_max, nr
    integer :: j, k, l

    l_max = size(this%phi,1) - 1
    nr = size(this%phi,2)

    ! RADIAL COORDINATE SPACE !

    ! atomic potential with centrifugal term
    allocate(this%V0(l_max+1,nr), this%dV0_dr(l_max+1,nr))
    do l=0,l_max
        this % V0(l+1,:) = -this%Z/r + real(l*(l+1),num)/2.0_num/(r**2)
        this % dV0_dr(l+1,:) = this%Z/(r**2) - real(l*(l+1),num)/(r**3)
    end do

    ! r-matrices are tri-diagonal
    ! nr diag + 2*(nr-1) off-diag elts = 3*nr-2 total elts
    ! row-major order dense-to-sparse convention
    allocate(this%D1r(3*nr-2), this%D2r(3*nr-2))
    allocate(this%M1r(3*nr-2), this%M2r(3*nr-2))

    ! sparse tri-diagonal element sequences:
    !   0 diag: 1, 3*nr-2, 3
    !  +1 diag: 2, 3*nr-4, 3
    !  -1 diag: 3, 3*nr-3, 3

    ! 1st derivative Hermitian-corrected centered-difference matrix
    this%D1r = 0.0
    do j=2, 3*nr-4, 3
        this%D1r(j) = 1.0
    end do
    do j=3, 3*nr-3, 3
        this%D1r(j) = -1.0
    end do
    this%D1r(1) = sqrt(3.0_num) - 2.0
    this%D1r(3*nr-2) = -sqrt(3.0_num) + 2.0
    this%D1r = this%D1r/(2.0*dr)

    ! 2nd derivative centered-difference matrix
    this%D2r = 1.0
    do j=1, 3*nr-2, 3
        this%D2r(j) = -2.0
    end do
    this%D2r = this%D2r/(dr**2)

    ! 1st- and 2nd-derivative Muller matrices
    this%M1r = 1.0/6.0_num
    this%M2r = -1.0/6.0_num
    do j=1, 3*nr-2, 3
        this%M1r(j) = 4.0/6.0_num
        this%M2r(j) = -10.0/6.0_num
    end do
    this%M1r(1) = this%M1r(1) + (sqrt(3.0_num)-2.0)/6.0
    this%M1r(3*nr-2) = this%M1r(3*nr-2) + (sqrt(3.0_num)-2.0)/6.0

    ! upper-element corrections for l=m=0
    this%D211 = -(2.0/dr**2)*(1.0 - this%Z*dr/(12.0 - 10.0*this%Z*dr))
    this%M211 = -2.0*(1.0 + (dr**2)*this%D211/12.0)

    ! set the populated indices of sparse tri-diagonal matrices
    allocate(this%stab(3*nr-2,2))
    ! key: stab(i,:) = (j,k) is the i-th tri-diagonal element at position (j,k)
    this%stab(1,:) = (/1,1/)
    this%stab(2,:) = (/1,2/)
    k = 2
    do j=3, 3*nr-6, 3
        this%stab(j,:)   = (/k, k-1/)
        this%stab(j+1,:) = (/k, k  /)
        this%stab(j+2,:) = (/k, k+1/)
        k = k + 1
    end do
    this%stab(3*nr-3,:) = (/nr, nr-1/)
    this%stab(3*nr-2,:) = (/nr, nr  /)

    ! ORBITAL ANGULAR MOMENTUM SPACE !

    allocate(this%clm(l_max))

    ! elements of Clm matrix
    do l=0,l_max-1
        this%clm(l+1) = sqrt(real((l+1)**2-this%m**2,num)/real(2*l+1,num)/real(2*l+3,num))
    end do

    return
end subroutine psi1dr_initialize_propagators

! __________________________________________________________________________________________________
!
! WAVEFUNCTION PROPAGATION ROUTINES
! __________________________________________________________________________________________________
!

subroutine psi1d_dt_propagate_fft(this, T, V, dx, dt, j)
    implicit none
    class(SchrodingerWavefunction1D), intent(inout) :: this
    real(num), intent(in) :: dx, dt
    complex, intent(in) :: j
    real(num), intent(in) :: T(:)
    complex(num), intent(in) :: V(:)

    complex(num), allocatable, save :: phi(:)

    if (.not.allocated(phi)) &
        allocate(phi(size(V)))

    ! SPLIT-OPERATOR METHOD !
    ! Apply 1st half potential energy part of Trotter expansion.
    this%psi = exp(-(i*j)*V*dt/2.0) * this%psi

    ! Fourier transform wavefunction; momentum-space representation.
    phi = fft(this%psi)

    ! Apply kinetic energy part of Trotter expansion.
    phi = exp(-(i*j)*T*dt) * phi

    ! Inverse Fourier transform wavefunction.
    this%psi = ifft(phi)

    ! Apply 2nd half potential energy part of Trotter expansion.
    this%psi = exp(-(i*j)*V*dt/2.0) * this%psi

    return
end subroutine psi1d_dt_propagate_fft

subroutine psi1d_dt_propagate_cnn(this, V, dt, j)
    implicit none
    class(SchrodingerWavefunction1D), intent(inout) :: this
    complex(num), intent(in) :: V(:)
    real(num), intent(in) :: dt
    complex, intent(in) :: j

    ! propagation variables
    complex(num), dimension(:), allocatable, save :: Ap, Am, yt
    integer :: nx
    nx = size(V)

    ! propagator Cayley form terms (tri-diagonal)
    if (.not.allocated(Ap)) then
        allocate(Ap(3*nx-2), Am(3*nx-2), yt(nx))
    end if

    ! Ap = M2x - (i*j)*(dt/2.0)*(-0.5*D2x + matmul(M2x,diag(V)))
    ! Am = M2x + (i*j)*(dt/2.0)*(-0.5*D2x + matmul(M2x,diag(V)))

    call tridiag_matmul_dmat(Ap, this%M2x, V)

    Ap = -0.5*this%D2x + Ap; Am = Ap;
    Ap = this%M2x - (i*j)*(dt/2.0)*Ap
    Am = this%M2x + (i*j)*(dt/2.0)*Am

    ! advance wavefunction: solve (A-).psi(t+dt) = (A+).psi(t)

    ! calculate yt = (A+).psi(t)
    call tridiag_matmul_cvec(yt, Ap, this%psi)

    ! solve (A-).psi(t+dt) = yt
    call tridiag_fbwd_subs(Am, this%psi, yt)

    return
end subroutine psi1d_dt_propagate_cnn

subroutine psi1dr_prepare_atomic_state(this, lp, ntau, dr, dt, pure)
    implicit none
    class(SchrodingerWavefunction1DR), intent(inout) :: this
    integer, intent(in) :: lp, ntau
    real(num), intent(in) :: dr, dt
    logical, intent(in) :: pure

    ! propagation variables
    real(num), dimension(:), allocatable :: Wp, Wm
    complex(num), dimension(:), allocatable :: rt

    integer :: l_max, nr
    integer :: l, it

    l_max = size(this%phi,1) - 1
    nr = size(this%phi,2)

    ! propagator Cayley form terms (tri-diagonal)
    allocate(Wp(3*nr-2), Wm(3*nr-2), rt(nr))

    if (pure) then

        !$OMP PARALLEL DO DEFAULT(SHARED)
        do l=0,l_max
            if (l /= lp) this%phi(l+1,:) = 0.0
        end do
        !$OMP END PARALLEL DO

        ! sparse calculation of Crank-Nicolson propagator components
        ! W(p/m) = M2r +/- (dt/2)*(D2r + matmul(M2r,diag(V0)))

        call tridiag_matmul_dmat(Wp, this%M2r, this%V0(lp+1,:))

        Wp = this%D2r + Wp; Wm = Wp;
        Wp = this%M2r + (dt/2.0)*Wp
        Wm = this%M2r - (dt/2.0)*Wm

        ! accuracy boost for Coulomb-like potentials
        if ((lp == 0).and.(this%m == 0)) then
            Wp(1) = this%M211 + (dt/2.0)*(this%D211 + this%M211*this%V0(1,1))
            Wm(1) = this%M211 - (dt/2.0)*(this%D211 + this%M211*this%V0(1,1))
        end if

        ! propagate the lp-state through imaginary time
        ! repeatedly solve (W+).phi(t+dt) = (W-).phi(t)
        do it=1,ntau
            ! calculate rt = (W-).phi(t)
            call tridiag_matmul_cvec(rt, Wm, this%phi(lp+1,:))

            ! solve (W+).phi(t+dt) = rt
            call tridiag_fbwd_subs(Wp, this%phi(lp+1,:), rt)

            ! normalize
            this%phi(lp+1,:) = this%phi(lp+1,:)/sqrt(dr*trapz(abs(this%phi(lp+1,:))**2))
        end do

    else

        do l=0,l_max
            ! sparse calculation of Crank-Nicolson propagator components
            ! W(p/m) = M2r +/- (dt/2)*(D2r + matmul(M2r,diag(V0)))

            call tridiag_matmul_dmat(Wp, this%M2r, this%V0(l+1,:))

            Wp = this%D2r + Wp; Wm = Wp;
            Wp = this%M2r + (dt/2.0)*Wp
            Wm = this%M2r - (dt/2.0)*Wm

            ! accuracy boost for Coulomb-like potentials
            if ((l == 0).and.(this%m == 0)) then
                Wp(1) = this%M211 + (dt/2.0)*(this%D211 + this%M211*this%V0(1,1))
                Wm(1) = this%M211 - (dt/2.0)*(this%D211 + this%M211*this%V0(1,1))
            end if

            ! propagate current l-state through imaginary time
            ! repeatedly solve (W+).phi(t+dt) = (W-).phi(t)
            do it=1,ntau
                ! calculate rt = (W-).phi(t)
                call tridiag_matmul_cvec(rt, Wm, this%phi(l+1,:))

                ! solve (W+).phi(t+dt) = rt
                call tridiag_fbwd_subs(Wp, this%phi(l+1,:), rt)

                ! normalize
                this%phi(l+1,:) = this%phi(l+1,:)/sqrt(dr*trapz(abs(this%phi(l+1,:))**2))
            end do
        end do

    end if

    return
end subroutine psi1dr_prepare_atomic_state

subroutine psi1dr_dt_propagate_free(this, dt)
    implicit none
    class(SchrodingerWavefunction1DR), intent(inout) :: this
    real(num), intent(in) :: dt

    ! propagation variables
    complex(num), dimension(:), allocatable, save :: Wp, Wm, rt

    integer :: l_max, nr
    integer :: l

    l_max = size(this%phi,1) - 1
    nr = size(this%phi,2)

    ! propagator Cayley form terms (tri-diagonal)
    if (.not.allocated(Wp)) then
        allocate(Wp(3*nr-2), Wm(3*nr-2), rt(nr))
    end if

    ! advance angular momentum components in-parallel
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(Wp, Wm, rt)
    do l=0,l_max
        ! Crank-Nicolson propagator components
        ! W(p/m) = M2r +/- (i*dt/2)*(D2r + matmul(M2r,diag(V0+Va)))

        call tridiag_matmul_dmat(Wp, this%M2r, this%V0(l+1,:)+this%Va)

        Wp = this%D2r + Wp; Wm = Wp;
        Wp = this%M2r + (i*dt/2.0)*Wp
        Wm = this%M2r - (i*dt/2.0)*Wm

        ! accuracy boost for Coulomb-like potentials
        if ((l == 0).and.(this%m == 0)) then
            Wp(1) = this%M211 + (i*dt/2.0)*(this%D211 + this%M211*this%V0(1,1))
            Wm(1) = this%M211 - (i*dt/2.0)*(this%D211 + this%M211*this%V0(1,1))
        end if

        ! advance wavefunction: solve (W+).phi(t+dt) = (W-).phi(t)

        ! calculate rt = (W-).phi(t)
        call tridiag_matmul_cvec(rt, Wm, this%phi(l+1,:))

        ! solve (W+).phi(t+dt) = rt
        call tridiag_fbwd_subs(Wp, this%phi(l+1,:), rt)
    end do
    !$OMP END PARALLEL DO

    return
end subroutine psi1dr_dt_propagate_free

subroutine psi1dr_dt_propagate_full(this, A, r, dt)
    implicit none
    class(SchrodingerWavefunction1DR), intent(inout) :: this
    real(num), intent(in) :: A, r(:), dt

    ! propagation variables
    real(num), allocatable, save :: Ylm(:,:), B(:,:), B_Rlm(:,:), Rlm_Bt(:,:)
    complex(num), allocatable, save :: rt(:), phi_l_vec(:), phi_lr_block(:,:)
    real(num) :: beta

    integer :: l_max, nr
    integer :: l, j, k

    l_max = size(this%phi,1) - 1
    nr = size(this%phi,2)

    if (.not.allocated(Ylm)) then
        allocate(Ylm(2,3*nr-2), B(2,2), B_Rlm(2,2), Rlm_Bt(2,2))
        allocate(rt(nr), phi_l_vec(2), phi_lr_block(2,nr))
        ! l-space B matrix
        B(1,:) = (/ 1.0, 1.0/)
        B(2,:) = (/-1.0, 1.0/)
        B = B/sqrt(2.0_num)
    end if

    ! STEP 1: apply 1st series of pure- and mixed-angular matrices !
    do l=l_max-1, 0, -1
        ! apply purely-angular 2x2 B_Rlm matrix to neighboring l-vectors, across r-space
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(beta, B_Rlm, phi_l_vec, j)
        do k=1,nr
            beta = -A*dt/(4.0*r(k))*real(l+1,num)*this%clm(l+1)
            B_Rlm(1,1) = (-beta**2 - 2*beta + 1)/sqrt(2.0_num)/(1 + beta**2)
            B_Rlm(1,2) = (-beta**2 + 2*beta + 1)/sqrt(2.0_num)/(1 + beta**2)
            B_Rlm(2,1) = (+beta**2 - 2*beta - 1)/sqrt(2.0_num)/(1 + beta**2)
            B_Rlm(2,2) = (-beta**2 - 2*beta + 1)/sqrt(2.0_num)/(1 + beta**2)

            phi_l_vec = matmul(B_Rlm, (/this%phi(l+1,k),this%phi(l+2,k)/))
            ForAll(j=1:2) this%phi(l+j,k) = phi_l_vec(j)
        end do
        !$OMP END PARALLEL DO

        ! create mixed-angular (Y-)lm matrix
        Ylm(1,:) = this%M1r - (A*dt*this%clm(l+1)/4.0)*this%D1r ! operates on l-th phi(r)
        Ylm(2,:) = this%M1r + (A*dt*this%clm(l+1)/4.0)*this%D1r ! operates on (l+1)th phi(r)

        do j=1,2
            ! apply (Y-)lm to (l+j)th phi(r)
            call tridiag_matmul_cvec(phi_lr_block(j,:), Ylm(j,:), this%phi(l+j,:))
            ! copy back to eigenmode matrix
            this%phi(l+j,:) = phi_lr_block(j,:)
        end do

        ! apply (Y+)lm inverse matrix
        ! now Ylm(2,:) operates on l-th phi'(r)
        ! and Ylm(1,:) operates on (l+1)th phi'(r)

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(rt)
        do j=1,2
            rt = this%phi(l+j,:)
            ! solve (Y+).phi'(r) = phi(r)
            call tridiag_fbwd_subs(Ylm(mod(j,2)+1,:), this%phi(l+j,:), rt)
        end do
        !$OMP END PARALLEL DO

        ! apply B-transpose to neighboring l-vectors
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(phi_l_vec,j)
        do k=1,nr
            phi_l_vec = matmul(transpose(B), (/this%phi(l+1,k),this%phi(l+2,k)/))
            ForAll(j=1:2) this%phi(l+j,k) = phi_l_vec(j)
        end do
        !$OMP END PARALLEL DO
    end do

    ! STEP 2: apply intermediate r-space Crank-Nicolson matrices !
    ! this operation is equivalent to field-free dt-propagation
    call psi1dr_dt_propagate_free(this,dt)

    ! STEP 3: apply last series of pure- and mixed-angular matrices !
    do l=0,l_max-1
        ! apply B to neighboring l-vectors
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(phi_l_vec,j)
        do k=1,nr
            phi_l_vec = matmul(B, (/this%phi(l+1,k),this%phi(l+2,k)/))
            ForAll(j=1:2) this%phi(l+j,k) = phi_l_vec(j)
        end do
        !$OMP END PARALLEL DO

        ! create mixed-angular (Y-)lm matrix
        Ylm(1,:) = this%M1r - (A*dt*this%clm(l+1)/4.0)*this%D1r ! operates on l-th phi(r)
        Ylm(2,:) = this%M1r + (A*dt*this%clm(l+1)/4.0)*this%D1r ! operates on (l+1)th phi(r)

        do j=1,2
            ! apply (Y-)lm to (l+j)th phi(r)
            call tridiag_matmul_cvec(phi_lr_block(j,:), Ylm(j,:), this%phi(l+j,:))
            ! copy back to eigenmode matrix
            this%phi(l+j,:) = phi_lr_block(j,:)
        end do

        ! apply (Y+)lm inverse matrix
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(rt)
        do j=1,2
            rt = this%phi(l+j,:)
            ! solve (Y+).phi'(r) = phi(r)
            call tridiag_fbwd_subs(Ylm(mod(j,2)+1,:), this%phi(l+j,:), rt)
        end do
        !$OMP END PARALLEL DO

        ! apply purely-angular 2x2 Rlm_Bt matrix to neighboring l-vectors, across r-space
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(beta, Rlm_Bt, phi_l_vec, j)
        do k=1,nr
            beta = -A*dt/(4.0*r(k))*real(l+1,num)*this%clm(l+1)
            Rlm_Bt(1,1) = (-beta**2 + 2*beta + 1)/sqrt(2.0_num)/(1 + beta**2)
            Rlm_Bt(1,2) = (+beta**2 + 2*beta - 1)/sqrt(2.0_num)/(1 + beta**2)
            Rlm_Bt(2,1) = (-beta**2 - 2*beta + 1)/sqrt(2.0_num)/(1 + beta**2)
            Rlm_Bt(2,2) = (-beta**2 + 2*beta + 1)/sqrt(2.0_num)/(1 + beta**2)

            phi_l_vec = matmul(Rlm_Bt, (/this%phi(l+1,k),this%phi(l+2,k)/))
            ForAll(j=1:2) this%phi(l+j,k) = phi_l_vec(j)
        end do
        !$OMP END PARALLEL DO
    end do

    return
end subroutine psi1dr_dt_propagate_full

subroutine psi2d_dt_propagate_fft(this, T, V, dr, dt, j)
    implicit none
    class(SchrodingerWavefunction2D), intent(inout) :: this
    real(num), intent(in) :: dr(2), dt
    complex, intent(in) :: j
    real(num), intent(in) :: T(:,:)
    complex(num), intent(in) :: V(:,:)

    complex(num), allocatable, save :: phi(:,:)

    if (.not.allocated(phi)) &
        allocate(phi(size(V,1),size(V,2)))

    ! SPLIT-OPERATOR METHOD !
    ! Apply 1st half potential energy part of Trotter expansion.
    this%psi = exp(-(i*j)*V*dt/2.0) * this%psi

    ! Fourier transform wavefunction; momentum-space representation.
    phi = fft(this%psi)

    ! Apply kinetic energy part of Trotter expansion.
    phi = exp(-(i*j)*T*dt) * phi

    ! Inverse Fourier transform wavefunction.
    this%psi = ifft(phi)

    ! Apply 2nd half potential energy part of Trotter expansion.
    this%psi = exp(-(i*j)*V*dt/2.0) * this%psi

    return
end subroutine psi2d_dt_propagate_fft

subroutine psi2d_dt_propagate_cnn(this, V, dt, j)
    implicit none
    class(SchrodingerWavefunction2D), intent(inout) :: this
    complex(num), intent(in) :: V(:,:)
    real(num), intent(in) :: dt
    complex, intent(in) :: j

    ! propagation variables
    complex(num), allocatable, save :: Apx(:), Amx(:), Apy(:), Amy(:), yt(:,:)
    integer :: nr(2), ix, iy

    nr(1) = size(V,1); nr(2) = size(V,2);

    ! propagator Cayley form terms (tri-diagonal)
    if (.not.allocated(Apx)) then
        allocate(Apx(3*nr(1)-2), Amx(3*nr(1)-2), Apy(3*nr(2)-2), Amy(3*nr(2)-2), yt(nr(1),nr(2)))
    end if

    ! Ap(x,y) = M2(x,y) - (i*j)*(dt/4.0)*(-D2(x,y) + matmul(M2(x,y),diag(V)))
    ! Am(x,y) = M2(x,y) + (i*j)*(dt/4.0)*(-D2(x,y) + matmul(M2(x,y),diag(V)))

    ! PEACEMAN-RACHFORD ALTERNATING DIRECTION IMPLICIT METHOD !
    ! unitary propagation for commutative x- and y-Hamiltonians

    ! apply (1 - i*dt*Hx/2)
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(Apx)
    do iy=1,nr(2)
        call tridiag_matmul_dmat(Apx, this%M2x, V(:,iy))
        Apx = -this%D2x + Apx; Apx = this%M2x - (i*j)*(dt/4.0)*Apx;
        call tridiag_matmul_cvec(yt(:,iy), Apx, this%psi(:,iy))
    end do
    !$OMP END PARALLEL DO

    ! invert (1 + i*dt*Hy/2)
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(Amy)
    do ix=1,nr(1)
        call tridiag_matmul_dmat(Amy, this%M2y, V(ix,:))
        Amy = -this%D2y + Amy; Amy = this%M2y + (i*j)*(dt/4.0)*Amy;
        call tridiag_fbwd_subs(Amy, this%psi(ix,:), yt(ix,:))
    end do
    !$OMP END PARALLEL DO

    ! apply (1 - i*dt*Hy/2)
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(Apy)
    do ix=1,nr(1)
        call tridiag_matmul_dmat(Apy, this%M2y, V(ix,:))
        Apy = -this%D2y + Apy; Apy = this%M2y - (i*j)*(dt/4.0)*Apy;
        call tridiag_matmul_cvec(yt(ix,:), Apy, this%psi(ix,:))
    end do
    !$OMP END PARALLEL DO

    ! invert (1 + i*dt*Hx/2)
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(Amx)
    do iy=1,nr(2)
        call tridiag_matmul_dmat(Amx, this%M2x, V(:,iy))
        Amx = -this%D2x + Amx; Amx = this%M2x + (i*j)*(dt/4.0)*Amx;
        call tridiag_fbwd_subs(Amx, this%psi(:,iy), yt(:,iy))
    end do
    !$OMP END PARALLEL DO

    return
end subroutine psi2d_dt_propagate_cnn

! __________________________________________________________________________________________________
!
! TYPE DESTRUCTORS
! __________________________________________________________________________________________________
!

subroutine psi1d_destructor(this)
    implicit none
    class(SchrodingerWavefunction1D), intent(inout) :: this

    if (allocated(this%psi)) deallocate(this%psi)
    if (allocated(this%grad_psi)) deallocate(this%grad_psi)
    if (allocated(this%phase)) deallocate(this%phase)
    if (allocated(this%norm)) deallocate(this%norm)
    if (allocated(this%energy)) deallocate(this%energy)
    if (allocated(this%D2x)) deallocate(this%D2x)
    if (allocated(this%M2x)) deallocate(this%M2x)
    if (allocated(this%stab)) deallocate(this%stab)

    return
end subroutine psi1d_destructor

subroutine psi1dr_destructor(this)
    implicit none
    class(SchrodingerWavefunction1DR), intent(inout) :: this

    if (allocated(this%V0)) deallocate(this%V0)
    if (allocated(this%Va)) deallocate(this%Va)
    if (allocated(this%phi)) deallocate(this%phi)
    if (allocated(this%clm)) deallocate(this%clm)
    if (allocated(this%D1r)) deallocate(this%D1r)
    if (allocated(this%D2r)) deallocate(this%D2r)
    if (allocated(this%M1r)) deallocate(this%M1r)
    if (allocated(this%M2r)) deallocate(this%M2r)
    if (allocated(this%norm)) deallocate(this%norm)
    if (allocated(this%stab)) deallocate(this%stab)
    if (allocated(this%energy)) deallocate(this%energy)
    if (allocated(this%dV0_dr)) deallocate(this%dV0_dr)

    return
end subroutine psi1dr_destructor

subroutine psi2d_destructor(this)
    implicit none
    class(SchrodingerWavefunction2D), intent(inout) :: this

    if (allocated(this%psi)) deallocate(this%psi)
    if (allocated(this%grad_psi)) deallocate(this%grad_psi)
    if (allocated(this%phase)) deallocate(this%phase)
    if (allocated(this%norm)) deallocate(this%norm)
    if (allocated(this%energy)) deallocate(this%energy)
    if (allocated(this%D2x)) deallocate(this%D2x)
    if (allocated(this%D2y)) deallocate(this%D2y)
    if (allocated(this%M2x)) deallocate(this%M2x)
    if (allocated(this%M2y)) deallocate(this%M2y)
    if (allocated(this%stabx)) deallocate(this%stabx)
    if (allocated(this%staby)) deallocate(this%staby)

    return
end subroutine psi2d_destructor

! __________________________________________________________________________________________________
!
! CALCULATE HAMILTONIAN EXPECTATION VALUE
! __________________________________________________________________________________________________
!

real(num) function expectE_1D_fft(this, T, V, dx, dp)
    implicit none
    class(SchrodingerWavefunction1D), intent(in) :: this
    real(num), intent(in), dimension(:) :: T, V
    real(num), intent(in) :: dx, dp

    complex(num), dimension(:), allocatable, save :: phi
    integer :: nx
    nx = size(V)

    if (.not.allocated(phi)) then
        allocate(phi(nx))
    end if

    phi = fft(this%psi)
    phi = phi/sqrt(dp*trapz(abs(phi)**2))

    expectE_1D_fft = dx*trapz(V*abs(this%psi)**2) + dp*trapz(T*abs(phi)**2)

    return
end function expectE_1D_fft

real(num) function expectE_1D_cnn(this, V, dx)
    implicit none
    class(SchrodingerWavefunction1D), intent(in) :: this
    real(num), intent(in) :: V(:)
    real(num), intent(in) :: dx

    complex(num), dimension(:), allocatable, save :: lt, rt
    integer :: nx
    nx = size(V)

    if (.not.allocated(lt)) then
        allocate(lt(nx), rt(nx))
    end if

    ! solve M2x.(psi')* = (psi)*
    rt = conjg(this%psi)
    call tridiag_fbwd_subs(this%M2x, lt, rt)

    ! calculate D2x.(psi)
    call tridiag_matmul_cvec(rt, this%D2x, this%psi)

    ! energy = -1/2 (psi,inv(M2).D2.psi) + (psi,V(x).psi)
    expectE_1D_cnn = -0.5*dx*trapz(real(lt*rt,num)) + dx*trapz(real(conjg(this%psi)*V*this%psi,num))

    return
end function expectE_1D_cnn

real(num) function expectE_1DR(this, dr)
    implicit none
    class(SchrodingerWavefunction1DR), intent(in) :: this
    real(num), intent(in) :: dr

    real(num), dimension(:), allocatable, save :: els, Ut
    complex(num), dimension(:), allocatable, save :: lt, rt

    integer :: l_max, nr
    integer :: l

    l_max = size(this%phi,1) - 1
    nr = size(this%phi,2)

    if (.not.allocated(els)) then
        allocate(els(l_max+1), Ut(3*nr-2), lt(nr), rt(nr))
    end if

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(Ut,lt,rt)
    do l=0,l_max
        ! solve M2r.(phi')* = (phi)*
        rt = conjg(this%phi(l+1,:))
        call tridiag_fbwd_subs(this%M2r, lt, rt)

        ! construct D2r + matmul(M2r,diag(V0))
        call tridiag_matmul_dmat(Ut, this%M2r, this%V0(l+1,:))
        Ut = this%D2r + Ut

        ! calculate Ut.phi
        call tridiag_matmul_cvec(rt, Ut, this%phi(l+1,:))

        ! calculate energy of current l-component
        els(l+1) = dr*trapz(real(lt*rt,num))
    end do
    !$OMP END PARALLEL DO

    ! calculate final energy
    expectE_1DR = sum(els)

    return
end function expectE_1DR

real(num) function expectE_2D_fft(this, T, V, dr, dp)
    implicit none
    class(SchrodingerWavefunction2D), intent(in) :: this
    real(num), intent(in), dimension(:,:) :: T, V
    real(num), intent(in) :: dr(2), dp(2)

    complex(num), dimension(:,:), allocatable, save :: phi
    integer :: nr(2)

    nr(1) = size(V,1); nr(2) = size(V,2);

    if (.not.allocated(phi)) then
        allocate(phi(nr(1),nr(2)))
    end if

    phi = fft(this%psi)
    phi = phi/sqrt(dp(1)*dp(2)*trapz(abs(phi)**2))

    expectE_2D_fft = dr(1)*dr(2)*trapz(V*abs(this%psi)**2) + dp(1)*dp(2)*trapz(T*abs(phi)**2)

    return
end function expectE_2D_fft

real(num) function expectE_2D_cnn(this, V, dr)
    implicit none
    class(SchrodingerWavefunction2D), intent(in) :: this
    real(num), intent(in) :: V(:,:), dr(2)

    complex(num), dimension(:,:), allocatable, save :: ltx, lty, rtx, rty
    integer :: nr(2), ix, iy

    nr(1) = size(V,1); nr(2) = size(V,2);

    if (.not.allocated(ltx)) then
        allocate(ltx(nr(1),nr(2)), lty(nr(1),nr(2)))
        allocate(rtx(nr(2),nr(1)), rty(nr(2),nr(1)))
    end if

    ! solve <psi'|.M2(x,y) = <psi|
    rtx = transpose(conjg(this%psi))
    !$OMP PARALLEL DO DEFAULT(SHARED)
    do iy=1,nr(2)
        call tridiag_fbwd_subs(this%M2x, ltx(:,iy), rtx(iy,:))
    end do
    !$OMP END PARALLEL DO

    rty = transpose(conjg(this%psi))
    !$OMP PARALLEL DO DEFAULT(SHARED)
    do ix=1,nr(1)
        call tridiag_fbwd_subs(this%M2y, lty(ix,:), rty(:,ix))
    end do
    !$OMP END PARALLEL DO

    ! calculate D2(x,y).(psi)
    !$OMP PARALLEL DO DEFAULT(SHARED)
    do iy=1,nr(2)
        call tridiag_matmul_cvec(rtx(iy,:), this%D2x, this%psi(:,iy))
    end do
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO DEFAULT(SHARED)
    do ix=1,nr(1)
        call tridiag_matmul_cvec(rty(:,ix), this%D2y, this%psi(ix,:))
    end do
    !$OMP END PARALLEL DO

    ! energy = -1/2 <psi|inv(M2x).D2x|psi> - 1/2 <psi|inv(M2y).D2y|psi> + <V(x,y)>
    expectE_2D_cnn = -0.5*dr(1)*dr(2)*trapz(real(ltx*transpose(rtx),num)) &
        - 0.5*dr(1)*dr(2)*trapz(real(lty*transpose(rty),num)) &
        + dr(1)*dr(2)*trapz(real(conjg(this%psi)*V*this%psi,num))

    return
end function expectE_2D_cnn

! __________________________________________________________________________________________________
!
! CALCULATE ELECTRON ENERGY SPECTRUM
! __________________________________________________________________________________________________
!

subroutine photoe_spectrum_winop_1D(psi, V0, E, dx, n, W)
    implicit none
    complex(num), intent(in) :: psi(:)
    real(num),    intent(in) :: V0(:), E(:), dx
    integer,      intent(in) :: n
    real(num),   intent(out) :: W(:)

    complex(num), dimension(:), allocatable, save :: OP, ft, rt

    real(num) :: gamma, q_nk
    integer   :: nx, nbins, v, k, j

    nx = size(V0); nbins = size(E);

    if (.not.allocated(OP)) &
        allocate(OP(3*nx-2), ft(nx), rt(nx))

    ! energy resolution
    gamma = E(2) - E(1)

    select case (n)
    ! bypass condition
    case (0)
    W = 0.0_num
    return

    ! first order
    case (1)
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(OP,ft,rt,j)
    do v=1,nbins
        rt = gamma*psi

        ! make window operator
        OP = -0.5/(dx**2)
        ForAll(j=1:nx) &
            OP(1 + 3*(j-1)) = 1.0/(dx**2) + V0(j) - (E(v) + i*gamma)

        ! solve for amplitude
        ! (H0 - Ev - i*gamma).Amp = gamma*psi
        call tridiag_fbwd_subs(OP, ft, rt)

        ! set spectrum value
        W(v) = dx*trapz(real(conjg(ft)*ft,num))
    end do
    !$OMP END PARALLEL DO

    ! general order
    case default
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(OP,ft,rt,j,k,q_nk)
    do v=1,nbins
        rt = (gamma**(2**(n-1)))*psi

        do k=1,2**(n-2)
            q_nk = real((2*k-1)*pi,num)/real(2**n,num)

            ! make window operator
            OP = -0.5/(dx**2)
            ForAll(j=1:nx) &
                OP(1 + 3*(j-1)) = 1.0/(dx**2) + V0(j) - E(v) - exp(i*q_nk)*gamma

            ! solve for amplitude
            ! (H0 - Ev - exp(i*q_nk)*gamma).(Amp') = Amp
            call tridiag_fbwd_subs(OP, ft, rt)
            rt = ft

            ! repeat for partner factorization
            OP = -0.5/(dx**2)
            ForAll(j=1:nx) &
                OP(1 + 3*(j-1)) = 1.0/(dx**2) + V0(j) - E(v) + exp(i*q_nk)*gamma

            ! iterate partial amplitude
            ! (H0 - Ev + exp(i*q_nk)*gamma).(Amp'') = Amp'
            call tridiag_fbwd_subs(OP, ft, rt)
            rt = ft
        end do

        ! set spectrum value
        W(v) = dx*trapz(real(conjg(ft)*ft,num))
    end do
    !$OMP END PARALLEL DO
    end select

    return
end subroutine photoe_spectrum_winop_1D

subroutine photoe_spectrum_winop_1DR(this, E, dr, n, W)
    implicit none
    class(SchrodingerWavefunction1DR), intent(in) :: this
    real(num), intent(in)  :: E(:), dr
    integer,   intent(in)  :: n
    real(num), intent(out) :: W(:)

    complex(num), dimension(:), allocatable, save :: Wls, OP, ft, rt
    real(num) :: gamma, q_nk

    integer :: l_max, nr, nbins
    integer :: l, j, v, k

    l_max = size(this%phi,1) - 1
    nr = size(this%phi,2)
    nbins = size(E)

    if (.not.allocated(Wls)) &
        allocate(Wls(l_max+1), OP(3*nr-2), ft(nr), rt(nr))

    ! energy resolution
    gamma = E(2) - E(1)

    select case (n)
    ! bypass condition
    case (0)
    W = 0.0_num
    return

    ! first order
    case (1)
    do v=1,nbins
        ! calculate partial spectra at current energy
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(OP,ft,rt,j)
        do l=0,l_max
            rt = gamma*this%phi(l+1,:)

            ! make window operator
            OP = -0.5/(dr**2)
            ForAll(j=1:nr) &
                OP(1 + 3*(j-1)) = 1.0/(dr**2) + this%V0(l+1,j) - (E(v) + i*gamma)

            ! solve for amplitude
            ! (H0 - Ev - i*gamma).Amp = gamma*psi
            call tridiag_fbwd_subs(OP, ft, rt)

            ! set partial spectrum value
            Wls(l+1) = dr*trapz(real(conjg(ft)*ft,num))
        end do
        !$OMP END PARALLEL DO

        ! set total spectrum of current energy
        W(v) = sum(Wls)
    end do

    ! general order
    case default
    do v=1,nbins
        ! calculate partial spectra at current energy
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(OP,ft,rt,j,k,q_nk)
        do l=0,l_max
            rt = (gamma**(2**(n-1)))*this%phi(l+1,:)

            do k=1,2**(n-2)
                q_nk = real((2*k-1)*pi,num)/real(2**n,num)

                ! make window operator
                OP = -0.5/(dr**2)
                ForAll(j=1:nr) &
                    OP(1 + 3*(j-1)) = 1.0/(dr**2) + this%V0(l+1,j) - E(v) - exp(i*q_nk)*gamma

                ! solve for amplitude
                ! (H0 - Ev - exp(i*q_nk)*gamma).(Amp') = Amp
                call tridiag_fbwd_subs(OP, ft, rt)
                rt = ft

                ! repeat for partner factorization
                OP = -0.5/(dr**2)
                ForAll(j=1:nr) &
                    OP(1 + 3*(j-1)) = 1.0/(dr**2) + this%V0(l+1,j) - E(v) + exp(i*q_nk)*gamma

                ! iterate partial amplitude
                ! (H0 - Ev + exp(i*q_nk)*gamma).(Amp'') = Amp'
                call tridiag_fbwd_subs(OP, ft, rt)
                rt = ft
            end do

            ! set partial spectrum value
            Wls(l+1) = dr*trapz(real(conjg(ft)*ft,num))
        end do
        !$OMP END PARALLEL DO

        ! set total spectrum of current energy
        W(v) = sum(Wls)
    end do
    end select

    return
end subroutine photoe_spectrum_winop_1DR

! __________________________________________________________________________________________________
!
! CALCULATE ELECTRON RADIATION EMISSION INTENSITY
! __________________________________________________________________________________________________
!

subroutine radiative_intensity_1DR(this, Et, dr, S)
    implicit none
    class(SchrodingerWavefunction1DR), intent(in) :: this
    real(num), intent(in) :: Et, dr
    complex(num), intent(out) :: S

    integer :: l, l_max
    l_max = size(this%phi,1) - 1

    S = 0.0_num

    do l=0,l_max
    if (l == 0) then
        S = S + dr*trapz(this%dV0_dr(l+1,:) * conjg(this%phi(l+1,:)) &
            * this%clm(l+2)*this%phi(l+2,:))
    else if (l == l_max) then
        S = S + dr*trapz(this%dV0_dr(l+1,:) * conjg(this%phi(l+1,:)) &
            * this%clm(l)*this%phi(l,:))
    else
        S = S + dr*trapz(this%dV0_dr(l+1,:) * conjg(this%phi(l+1,:)) &
            * (this%clm(l)*this%phi(l,:) + this%clm(l+2)*this%phi(l+2,:)))
    end if
    end do

    ! expectation value of acceleration in direction of polarization
    S = -(S + Et)

    return
end subroutine radiative_intensity_1DR

subroutine radiative_intensity_2D1e(this, aV, Et, dr, S)
    implicit none
    class(SchrodingerWavefunction2D), intent(in) :: this
    real(num), intent(in) :: aV(:,:,:), Et(2), dr(2)
    complex(num), intent(out) :: S
    S = -dr(1)*dr(2)*trapz(conjg(this%psi)*aV(:,:,1)*this%psi) &
        -dr(1)*dr(2)*trapz(conjg(this%psi)*aV(:,:,2)*this%psi) - sum(Et)
    return
end subroutine radiative_intensity_2D1e

subroutine radiative_intensity_1D2e(this, aV, Et, dr, S)
    implicit none
    class(SchrodingerWavefunction2D), intent(in) :: this
    real(num), intent(in) :: aV(:,:,:), Et, dr(2)
    complex(num), intent(out) :: S(2)
    S(1) = -dr(1)*dr(2)*trapz(conjg(this%psi)*aV(:,:,1)*this%psi) - Et
    S(2) = -dr(1)*dr(2)*trapz(conjg(this%psi)*aV(:,:,2)*this%psi) - Et
    return
end subroutine radiative_intensity_1D2e

! __________________________________________________________________________________________________
!
! TIME-DEPENDENT SURFACE-FLUX METHOD PROCEDURES
! __________________________________________________________________________________________________
!

subroutine tSURFF2D_initialize(S, x, y, dt)
    implicit none
    class(tSURFF2D), intent(inout) :: S
    real(num), intent(in) :: x(:), y(:), dt
    integer :: n

    if (S%enable) then

        ! set momentum distribution bins
        allocate(S%kx(S%nk(1)), S%ky(S%nk(2)))
        S%kx = linspace(S%kx_lim(1), S%kx_lim(2), S%nk(1))
        S%ky = linspace(S%ky_lim(1), S%ky_lim(2), S%nk(2))

        ! set surface integration coordinates
        allocate(S%xl(S%Ns), S%yl(S%Ns))
        allocate(S%ixl(S%Ns), S%iyl(S%Ns))

        do n=1,S%Ns
            S%xl(n) = S%R0 * cos(2.*pi*(n-1)/real(S%Ns,num))
            S%yl(n) = S%R0 * sin(2.*pi*(n-1)/real(S%Ns,num))

            ! index in spatial grid is nearest upper-left point
            S%ixl(n) = minloc(abs(x-S%xl(n)),1)
            S%iyl(n) = minloc(abs(y-S%yl(n)),1)
            if (x(S%ixl(n)) >= S%xl(n)) S%ixl(n) = S%ixl(n) - 1
            if (y(S%iyl(n)) >= S%yl(n)) S%iyl(n) = S%iyl(n) - 1
        end do

        ! angular spacing between surface points
        S%dphi = 2.*pi/real(S%Ns,num)

        ! blank amplitude
        allocate(S%p_dist(S%nk(1),S%nk(2)))
        S%p_dist = 0.0_num

        ! set integration step frequency
        S%iti = nint(S%dti/dt)

    else

        ! nullify unused variables
        allocate(S%kx(1), S%ky(1))
        allocate(S%xl(1), S%yl(1))
        allocate(S%ixl(1), S%iyl(1))
        allocate(S%p_dist(1,1))
        S%kx = 0.0_num; S%ky = 0.0_num;
        S%xl = 0.0_num; S%yl = 0.0_num;
        S%ixl = 0; S%iyl = 0; S%Ns = 1;
        S%dphi = 0.0_num; S%p_dist = 0.0_num;
        S%iti = huge(1_4); S%nk = (/1,1/);

    end if

    return
end subroutine tSURFF2D_initialize

subroutine tSURFF2D_dt_step(S, wavefn, A, C, x, y, t)
    implicit none
    class(tSURFF2D), intent(inout) :: S
    class(SchrodingerWavefunction2D), intent(in) :: wavefn
    real(num), intent(in) :: A(2), C(2), x(:), y(:), t

    ! integration variables
    complex(num) :: term(S%Ns), volkovc
    real(num) :: r(2), k(2), r_hat(2), dx(2), dy(2), r_mag, k_para, k_sq
    integer :: j, l, ikx, iky, ir(2)

    ! temporary interpolation variables
    complex(num), dimension(:,:), allocatable, save :: grad12_psi
    complex(num), dimension(4) :: g, g1, g2, g12
    complex(num) :: FR, F1R, F2R
    integer :: m(4,2)

    if (.not.S%enable) return

    if (S%interp == 'b3') then
        if (.not.allocated(grad12_psi)) &
        allocate(grad12_psi(size(x),size(y)))
        call deriv(wavefn%grad_psi(:,:,1), y(2)-y(1), 2, grad12_psi)
    end if

    ! dt-surface-integrate each (kx,ky) cell
    !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(S,wavefn,A,C,x,y,t,grad12_psi)
    do ikx=1,S%nk(1)
    do iky=1,S%nk(2)
        ! current wave vector
        k = (/S%kx(ikx),S%ky(iky)/); k_sq = sum(k*k);

        ! evaluate surface integration terms
        do j=1,S%Ns

            ! current position vector
            r = (/S%xl(j),S%yl(j)/); ir = (/S%ixl(j),S%iyl(j)/);
            r_mag = sqrt(sum(r*r)); r_hat = r/r_mag; k_para = sum(k*r_hat);

            ! (complex-conjugate of) associated Volkov wave
            volkovc = exp(i*k_sq*t/2.)*exp(-i*sum(k*(r-C)))/(2.*pi);

            ! interpolate wavefunction to this surface point
            select case (S%interp)
            ! nearest-neighbor
            case ('nn')
                FR = wavefn % psi(ir(1),ir(2))
                F1R = wavefn % grad_psi(ir(1),ir(2),1)
                F2R = wavefn % grad_psi(ir(1),ir(2),2)
            ! linear
            case ('li')
                dx = (/ x(ir(1)+1)-r(1) , r(1)-x(ir(1)) /)
                dy = (/ y(ir(2)+1)-r(2) , r(2)-y(ir(2)) /)

                FR = dx(1)*dy(1)*wavefn%psi(ir(1),ir(2)) &
                   + dx(1)*dy(2)*wavefn%psi(ir(1),ir(2)+1) &
                   + dx(2)*dy(1)*wavefn%psi(ir(1)+1,ir(2)) &
                   + dx(2)*dy(2)*wavefn%psi(ir(1)+1,ir(2)+1);

                F1R = dx(1)*dy(1)*wavefn%grad_psi(ir(1),ir(2),1) &
                    + dx(1)*dy(2)*wavefn%grad_psi(ir(1),ir(2)+1,1) &
                    + dx(2)*dy(1)*wavefn%grad_psi(ir(1)+1,ir(2),1) &
                    + dx(2)*dy(2)*wavefn%grad_psi(ir(1)+1,ir(2)+1,1);

                F2R = dx(1)*dy(1)*wavefn%grad_psi(ir(1),ir(2),2) &
                    + dx(1)*dy(2)*wavefn%grad_psi(ir(1),ir(2)+1,2) &
                    + dx(2)*dy(1)*wavefn%grad_psi(ir(1)+1,ir(2),2) &
                    + dx(2)*dy(2)*wavefn%grad_psi(ir(1)+1,ir(2)+1,2);

                FR = FR/(x(2)-x(1))/(y(2)-y(1))
                F1R = F1R/(x(2)-x(1))/(y(2)-y(1))
                F2R = F2R/(x(2)-x(1))/(y(2)-y(1))
            ! bicubic
            case ('b3')
                m(1,:) = (/ir(1)+1, ir(2)/); m(2,:) = (/ir(1)+1, ir(2)+1/);
                m(3,:) = (/ir(1), ir(2)+1/); m(4,:) = (/ir(1), ir(2)/);

                do l=1,4
                    g(l) = wavefn % psi(m(l,1),m(l,2))
                    g1(l) = wavefn % grad_psi(m(l,1),m(l,2),1)
                    g2(l) = wavefn % grad_psi(m(l,1),m(l,2),2)
                    g12(l) = grad12_psi(m(l,1),m(l,2))
                end do

                call bcuint(g, g1, g2, g12, x(ir(1)), y(ir(2)), &
                    r(1), r(2), (/x(2)-x(1),y(2)-y(1)/), FR, F1R, F2R)
            ! default: nearest-neighbor
            case default
                FR = wavefn % psi(ir(1),ir(2))
                F1R = wavefn % grad_psi(ir(1),ir(2),1)
                F2R = wavefn % grad_psi(ir(1),ir(2),2)
            end select

            ! evaluate differential-angle t-SURFF probability amplitude
            term(j) = r_mag * volkovc &
                * (sum(r_hat*A)*FR - (i/2.0)*(sum(r_hat*(/F1R,F2R/)) + i*sum(k*r_hat)*FR))

        end do

        ! evaluate t-SURFF probability amplitude
        S%p_dist(ikx,iky) = S%p_dist(ikx,iky) + (S%dti/2.0)*S%dphi*trapz(term)
    end do
    end do
    !$OMP END PARALLEL DO

    return
end subroutine tSURFF2D_dt_step

subroutine tSURFF2D_destructor(this)
    implicit none
    class(tSURFF2D), intent(inout) :: this
    if (allocated(this%p_dist)) deallocate(this%p_dist)
    if (allocated(this%kx)) deallocate(this%kx)
    if (allocated(this%ky)) deallocate(this%ky)
    if (allocated(this%xl)) deallocate(this%xl)
    if (allocated(this%yl)) deallocate(this%yl)
    if (allocated(this%ixl)) deallocate(this%ixl)
    if (allocated(this%iyl)) deallocate(this%iyl)
    return
end subroutine tSURFF2D_destructor

! __________________________________________________________________________________________________
!
! MISCELLANEOUS PROCEDURES
! __________________________________________________________________________________________________
!

subroutine calc_bohm_velocity(bv, Jxt, psi, dx, nx, nt, k)
    implicit none
    integer,      intent(in)                    :: nx, nt, k
    real(num),    intent(in)                    :: dx
    complex(num), intent(in),  dimension(nt,nx) :: psi

    real(num),    intent(out), dimension(nt,nx) :: bv, Jxt
    real(num),                 dimension(nx)    :: rho
    complex(num),              dimension(nx)    :: dpsi_dx, dpsiC_dx

    integer :: j

    rho = abs(psi(k,:))**2

    ! compute derivative of wavefunction
    dpsi_dx = grad(psi(k,:),dx)
    dpsiC_dx = conjg(dpsi_dx)

    ! compute probability current
    ! note: Im(J) is zero within numerical error.
    Jxt(k,:) = real((i/2.0)*(psi(k,:)*dpsiC_dx - conjg(psi(k,:))*dpsi_dx), num)

    ! compute Bohmian velocity
    bv(k,:) = Jxt(k,:)/rho
    return
end subroutine calc_bohm_velocity

subroutine calc_bohm_velocity_from_phase(bv, psi, dx, nx, nt, k)
    implicit none
    integer,      intent(in)                    :: nx, nt, k
    real(num),    intent(in)                    :: dx
    complex(num), intent(in),  dimension(nt,nx) :: psi
    real(num),    intent(out), dimension(nt,nx) :: bv
    real(num),                 dimension(nx)    :: S

    real(num) :: delta
    integer   :: j

    ! compute phase of wavefunction
    S = atan2(aimag(psi(k,:)),real(psi(k,:)))

    ! phase-unwrap S(x,t)
    call phase_unwrap(S)

    ! compute Bohmian velocity
    bv(k,:) = grad(S,dx)
    return
end subroutine calc_bohm_velocity_from_phase

subroutine calc_bohm_trajectories(bx, bv, x, dt, nx, nt)
    implicit none
    integer,   intent(in) :: nx, nt
    real(num), intent(in) :: dt
    real(num), intent(in),  dimension(nx)    :: x
    real(num), intent(in),  dimension(nt,nx) :: bv
    real(num), intent(out), dimension(nt,nx) :: bx

    integer :: j

    ! numerically integrate the velocity field at each position
    do j=1,nx
        bx(:,j) = x(j) + simint(bv(:,j), bv(1,j), dt)
    end do
    return
end subroutine calc_bohm_trajectories

subroutine chk_continuity_eqn(cty, psi, Jxt, dx, dt, nx, nt)
    implicit none
    integer,      intent(in)                    :: nx, nt
    real(num),    intent(in)                    :: dx, dt
    real(num),    intent(in),  dimension(nt,nx) :: Jxt
    complex(num), intent(in),  dimension(nt,nx) :: psi
    real(num),    intent(out), dimension(nt,nx) :: cty

    real(num), dimension(nt,nx) :: rho
    real(num), dimension(nt,nx) :: drho_dt
    real(num), dimension(nt)    :: divJ

    integer :: j,k

    rho = abs(psi)**2

    ! calculate d(rho)/dt
    call deriv(rho,dt,1,drho_dt)

    ! calculate div(J)
    do k=1,nt
        divJ(k) = sum(grad(Jxt(k,:),dx))
    end do

    ! evaluate d(rho)/dt + div(J)
    do k=1,nt
        cty(k,:) = drho_dt(k,:) + divJ(k)*ones(nx)
    end do
    return
end subroutine chk_continuity_eqn

function E_hydrogen(n,l) result(En)
    implicit none
    integer, intent(in) :: n, l
    real(dble_t) :: En, c, j
    type(pconst_mks) :: mks

    c = 1.0/mks%alpha
    j = l + 0.5

    En = (c**2)/sqrt(1.0 + (mks%alpha**2)/(real(n,num)-j-0.5+sqrt((j+0.5)**2-(mks%alpha**2)))**2)
    En = En - (c**2) ! subtract rest-mass energy

    return
end function E_hydrogen

end module quantum
