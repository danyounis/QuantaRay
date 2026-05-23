! QED Monte-Carlo Particle Module (qed_mcpm.F08)
!
! Author:  D. Younis
!          University of Rochester
!          Department of Physics
!
! Written: 1/24/2025
! Revised: 3/20/2026

module qed_mcpm

use prec
use math, only: pi, operator(.cross.), operator(.dot.)
use quantum, only: pconst_mks, pconst_cgs
use optimize, only: BinarySearch

implicit none

! Units
real(num) :: lam0_cm ! defined
real(num) :: k0_icm, w0_Hz, hbar_, q_, F_ ! derived

! Parameters
real(num) :: Ep_min = 1.0_num
real(num) :: dt_max = 1.0_num
real(num) :: dt_scl = real(1.d-2, num)
logical :: photon_dynamics = .true.

type(pconst_mks) :: mks
type(pconst_cgs) :: cgs

! __________________________________________________________________________________________________
!
! PARTICLE & LINKED-LIST TYPES
! __________________________________________________________________________________________________
!

type, public :: particle
    real :: m, q
    real(num) :: x(0:3), p(0:3), dt
    real(num) :: qp = 0.0_num, x0_, p0_
    integer :: id, stat
contains
    generic, public :: push => push_particle_rk4, push_photon
    procedure, private :: push_particle_rk4, push_photon
end type particle

type, public :: pnode_t
    type(particle) :: par
    type(pnode_t), pointer :: prev=>null()
    type(pnode_t), pointer :: next=>null()
end type pnode_t

type, public :: plist_t
    type(pnode_t), pointer :: head=>null()
    type(pnode_t), pointer :: tail=>null()
    integer :: N=0
contains
! core
    procedure :: create => create_particle
    procedure :: delete => delete_particle
    procedure :: reset => reset_list
    procedure :: push => list_push
    procedure :: emit => list_emit
! support
    procedure :: count => count_particles
    procedure :: has_stat => exists_particle_wstat
    procedure :: get_id => get_particle_ids
    procedure :: get_charge => get_particle_charges
    procedure :: get_energy => get_particle_energies
    procedure :: get_qparam => get_particle_quantum_parameters
    procedure :: get_ivalue => get_particle_initial_quantities
    procedure :: get_fvalue => get_particle_final_quantities
end type plist_t

type(plist_t) :: annihilated_photons

type, public :: diag_t
    real(num) :: min_qp(2) = huge(1.d0)
    real(num) :: max_qp(2) = 0.0_num
    real(num) :: min_prob(2) = huge(1.d0)
    real(num) :: max_prob(2) = 0.0_num
end type diag_t

type(diag_t) :: diag

! __________________________________________________________________________________________________
!
! FIELD, FORCE, AND RATE-TABLE TYPES
! __________________________________________________________________________________________________
!

type, public :: field
    procedure(emfunc), pointer, nopass :: Ex=>null(), Ey=>null(), Ez=>null()
    procedure(emfunc), pointer, nopass :: Bx=>null(), By=>null(), Bz=>null()
contains
! core
    procedure :: set => set_field_components
    procedure :: unset => unset_field_components
! support
    procedure :: E => get_field_electric
    procedure :: B => get_field_magnetic
end type field

type, public :: force
    procedure(frfunc), pointer, nopass :: eval=>null()
contains
    procedure :: set => set_force_function
    procedure :: unset => unset_force_function
end type force

type, public :: qmc_table_t
    real(kind=dble_t), dimension(:), allocatable :: eta0, chi1, frac1, eta2, h2, chi3, g3
    real(kind=dble_t), dimension(:,:), allocatable :: chi0, Py0, Pf1
    character(len=:), allocatable :: path
contains
    procedure :: load => load_qmc_tables
end type qmc_table_t

abstract interface
    ! field function template
    function emfunc(r,t)
        use prec, only: num
        implicit none
        real(num), intent(in) :: r(3), t
        real(num) :: emfunc
    end function emfunc

    ! force function template
    function frfunc(m,q,r,p,t,F)
        use prec, only: num
        import field
        implicit none
        real, intent(in) :: m, q
        real(num), intent(in) :: r(3), p(3), t
        class(field), intent(in) :: F
        real(num) :: frfunc(3)
    end function frfunc
end interface

interface Norm
    module procedure Norm3
    module procedure Norm4
end interface Norm

contains

! __________________________________________________________________________________________________
!
! MODULE PARAMETER ROUTINES
! __________________________________________________________________________________________________
!

subroutine set_qedmcpm_units(lam0_cm_)
    implicit none
    real(num), intent(in) :: lam0_cm_
    lam0_cm = lam0_cm_ ! normalization length [cm]
    k0_icm = (2*pi)/lam0_cm ! wavenumber [1/cm]
    w0_Hz = cgs%c*k0_icm ! frequency [Hz]
    hbar_ = cgs%hbar*k0_icm/cgs%m_e/cgs%c ! reduced Planck constant [sim.u.]
    q_ = sqrt(cgs%alpha*hbar_) ! electric charge [sim.u.]
    F_ = (w0_Hz*cgs%m_e*cgs%c)/cgs%e ! E/B-field conversion factor [CGS]->[sim.u.]
    return
end subroutine set_qedmcpm_units

subroutine set_qedmcpm_params(Ep_min_,dt_max_,dt_scl_,photon_dynamics_)
    implicit none
    real(num), intent(in), optional :: Ep_min_, dt_max_, dt_scl_
    logical, intent(in), optional :: photon_dynamics_
    if (present(Ep_min_)) Ep_min = Ep_min_
    if (present(dt_max_)) dt_max = dt_max_
    if (present(dt_scl_)) dt_scl = dt_scl_
    if (present(photon_dynamics_)) photon_dynamics = photon_dynamics_
    return
end subroutine set_qedmcpm_params

! __________________________________________________________________________________________________
!
! PARTICLE PROPAGATION ROUTINES
! __________________________________________________________________________________________________
!

subroutine push_particle_rk4(par,G,F)
    implicit none
    class(particle), intent(inout) :: par
    class(force), intent(in) :: G
    class(field), intent(in) :: F

    real(num) :: p_dot(3), rel
    integer :: j

    ! 4th-order Runge-Kutta (RK4) coefficients
    real(num) :: r_(3), p_(3), t_
    real(num) :: x1_(4), x2_(4), x3_(4)
    real(num) :: p1_(4), p2_(4), p3_(4)

    ! compute RK4 coefficients
    do j=1,4
        select case (j)
        case (1)
            r_ = par%x(1:3)
            p_ = par%p(1:3)
            t_ = par%x(0)
        case (2:3)
            r_ = [par%x(1) + 0.5*x1_(j-1), par%x(2) + 0.5*x2_(j-1), par%x(3) + 0.5*x3_(j-1)]
            p_ = [par%p(1) + 0.5*p1_(j-1), par%p(2) + 0.5*p2_(j-1), par%p(3) + 0.5*p3_(j-1)]
            t_ = par%x(0) + 0.5*par%dt
        case (4)
            r_ = [par%x(1)+x1_(j-1), par%x(2)+x2_(j-1), par%x(3)+x3_(j-1)]
            p_ = [par%p(1)+p1_(j-1), par%p(2)+p2_(j-1), par%p(3)+p3_(j-1)]
            t_ = par%x(0) + par%dt
        end select

        p_dot = G%eval(par%m,par%q,r_,p_,t_,F)
        rel = Lorentz(p_,par%m)

        x1_(j) = par%dt*p_(1)/(rel*par%m)
        x2_(j) = par%dt*p_(2)/(rel*par%m)
        x3_(j) = par%dt*p_(3)/(rel*par%m)
        p1_(j) = par%dt*p_dot(1)
        p2_(j) = par%dt*p_dot(2)
        p3_(j) = par%dt*p_dot(3)
    end do

    ! dt-advance 4-position and 4-momentum
    par%x(1) = par%x(1) + (x1_(1) + 2*x1_(2) + 2*x1_(3) + x1_(4))/6.0
    par%x(2) = par%x(2) + (x2_(1) + 2*x2_(2) + 2*x2_(3) + x2_(4))/6.0
    par%x(3) = par%x(3) + (x3_(1) + 2*x3_(2) + 2*x3_(3) + x3_(4))/6.0
    par%p(1) = par%p(1) + (p1_(1) + 2*p1_(2) + 2*p1_(3) + p1_(4))/6.0
    par%p(2) = par%p(2) + (p2_(1) + 2*p2_(2) + 2*p2_(3) + p2_(4))/6.0
    par%p(3) = par%p(3) + (p3_(1) + 2*p3_(2) + 2*p3_(3) + p3_(4))/6.0
    par%x(0) = par%x(0) + par%dt
    par%p(0) = par%m*Lorentz(par%p(1:3),par%m)
    return
end subroutine push_particle_rk4

subroutine push_particle_vay(par,G,F)
    implicit none
    class(particle), intent(inout) :: par
    class(force), intent(in) :: G
    class(field), intent(in) :: F

    real(num) :: dt_h, qdt_hm, E(3), B(3)
    real(num) :: us, sig, s
    real(num) :: u_(2,3), up(3), tau(3), gam(2), tt(3)

    ! Note: This is a leapfrog scheme; momentum (p) half-advanced by dt/2 over position (x). !

    dt_h = 0.5*par%dt
    qdt_hm = dt_h*par%q/par%m
    E = F%E(par%x(1:3), par%x(0)+dt_h) ! E(r^{i+1/2}, t+dt/2)
    B = F%B(par%x(1:3), par%x(0)+dt_h) ! B(r^{i+1/2}, t+dt/2)

    ! u^{i+1/2}
    u_(1,:) = par%p(1:3)/par%m &
        + qdt_hm * (E + (par%p(1:3).cross.B)/Lorentz(par%p(1:3),par%m)/par%m)

    up = u_(1,:) + qdt_hm*E ! u-prime
    tau = qdt_hm*B
    us = up.dot.tau ! u*

    gam(1) = sqrt(1.0 + sum(up**2)) ! gamma-prime
    sig = gam(1)**2 - sum(tau**2) ! sigma
    gam(2) = sqrt(sig + sqrt(sig**2 + 4*(sum(tau**2) + us**2)))/sqrt(2.) ! gamma^{i+1/2}
    tt = tau/gam(2) ! t-vector
    s = (1.0 + sum(tt**2))**(-1)

    u_(2,:) = s*(up + (up.dot.tt)*tt + (up.cross.tt)) ! u^{i+1}

    ! dt-advance 4-position and 4-momentum
    par%x(1:3) = par%x(1:3) + par%dt*par%p(1:3)/Lorentz(par%p(1:3),par%m)/par%m
    par%x(0) = par%x(0) + par%dt
    par%p(1:3) = par%m*u_(2,:)
    par%p(0) = par%m*Lorentz(par%p(1:3),par%m)
    return
end subroutine push_particle_vay

subroutine push_photon(par)
    implicit none
    class(particle), intent(inout) :: par
    integer :: j
    par%x(0) = par%x(0) + par%dt
    ForAll(j=1:3) &
        par%x(j) = par%x(j) + par%dt*par%p(j)/par%p(0)
    return
end subroutine push_photon

! __________________________________________________________________________________________________
!
! PARTICLE CREATION & ANNIHILATION PROCEDURES
! __________________________________________________________________________________________________
!

subroutine create_particle(lst,par)
    implicit none
    class(plist_t), intent(inout) :: lst
    class(particle), intent(in) :: par
    class(pnode_t), pointer :: new_node
    ! instantiate node
    allocate(new_node)
    new_node%par = par
    new_node%prev => lst%tail
    new_node%next => null()
    if (.not.associated(lst%head)) then
        ! the list is empty
        lst%head => new_node
        lst%tail => new_node
        lst%N = 1 ! initialize particle count
    else
        ! append to the list
        lst%tail%next => new_node
        lst%tail => new_node
        lst%N = lst%N+1 ! increment particle count
    end if
    return
end subroutine create_particle

subroutine delete_particle(lst,node)
    implicit none
    class(plist_t), intent(inout) :: lst
    class(pnode_t), intent(inout), pointer :: node
    class(pnode_t), pointer :: curr
    curr => lst%head

    outer: do while (associated(curr))
        if (associated(curr,node)) then
            ! singular node
            if (.not.associated(curr%prev).and.&
                .not.associated(curr%next)) then
                lst%head => null()
                lst%tail => null()
                exit outer

            ! head node
            else if (.not.associated(curr%prev)) then
                lst%head => curr%next
                curr%next%prev => null()
                exit outer

            ! tail node
            else if (.not.associated(curr%next)) then
                lst%tail => curr%prev
                curr%prev%next => null()
                exit outer
            end if

            ! middle node
            curr%prev%next => curr%next
            curr%next%prev => curr%prev
            exit outer
        end if
        curr => curr%next
    end do outer

    if (associated(curr)) then
        deallocate(curr) ! release the memory
        nullify(node) ! no dangling pointer
        lst%N = lst%N-1 ! decrement particle count
    end if
    return
end subroutine delete_particle

subroutine reset_list(lst)
    implicit none
    class(plist_t), intent(inout) :: lst
    class(pnode_t), pointer :: curr, next

    curr => lst%head
    do while (associated(curr))
        next => curr%next
        deallocate(curr)
        curr => next
    end do

    nullify(lst%head)
    nullify(lst%tail)
    lst%N = 0 ! reset particle count

    return
end subroutine reset_list

function new_particle(m,q,id,x,p,E,N,dt,stat) result(par)
    type(particle) :: par
    real, intent(in) :: m, q
    integer, intent(in) :: id
    real(num), intent(in) :: x(0:3)
    real(num), intent(in), optional :: p(3), E, N(3), dt
    integer, intent(in), optional :: stat
    real(num) :: rel, N_(3)
    logical :: C(3)

#ifndef UNSAFE
    C(1) = .not.(present(p).or.present(E).or.present(N)) ! ~p ^ ~E ^ ~N
    C(2) = (.not.present(p)).and.((.not.present(E)).or.(.not.present(N))) ! ~p ^ (~E v ~N)
    C(3) = present(p).and.(present(E).or.present(N)) ! p ^ (E v N)

    if (C(1).or.C(2).or.C(3)) &
        call ERROR(message='Particle energy/momentum under or over-specified.', &
          fix='Set the 3-momentum `p` /or/ relativistic energy `E` and direction-vector `N`.', &
          in='function new_particle')
#endif

    par%m = m
    par%q = q
    par%id = id
    par%x = x

    if (present(dt)) then; par%dt = dt; else; par%dt = 1.d-3; end if; ! particle-push time-step
    if (present(stat)) then; par%stat = stat; else; par%stat = 1; end if; ! default to normal initial status
    if ((m==0.).and.(q==0.).and.(.not.photon_dynamics)) par%stat = 0; ! freeze newly-created photon

    if (present(p)) then
        par%p(1:3) = p
        if (m /= 0.0) then; par%p(0) = m*Lorentz(p,m); ! massive particle
        else; par%p(0) = Norm(p); end if; ! photon
        return
    end if

    if (present(E).and.present(N)) then
#ifndef UNSAFE
        if (Norm(N) == 0.0_num) &
            call ERROR('Particle direction vector `N` is ill-defined (zero norm).', in='function new_particle')
#endif
        rel = E/m ! relativistic factor
        N_ = N/Norm(N) ! normalized direction
        par%p(0) = E
        if (m /= 0.0) then
            par%p(1:3) = E*N_*sqrt(1.0 - rel**(-2))
        else
            par%p(1:3) = E*N_
        end if
    end if

    ! massive particle initialized at rest
    if (par%p(0) == real(m,num)) par%p(1:3) = 0.0_num

#ifndef UNSAFE
    if (par%p(0) < real(m,num)) &
        call ERROR('Particle energy cannot be less than rest-mass (E !< m).', in='function new_particle')
#endif

    ! store time/energy at birth
    par%x0_ = par%x(0)
    par%p0_ = par%p(0)

end function new_particle

! __________________________________________________________________________________________________
!
! PARTICLE-LIST PHYSICS ROUTINES
! __________________________________________________________________________________________________
!

subroutine list_push(lst,t_max,enrg_min,dt,G,F)
    implicit none
    class(plist_t), intent(inout) :: lst
    real(num), intent(in), optional :: t_max, enrg_min, dt
    class(force), intent(in), optional :: G
    class(field), intent(in), optional :: F

    class(pnode_t), pointer :: curr
    class(particle), pointer :: par
    real(num) :: maxrate

    logical, save :: suppress_warning(2) = .false.

    curr => lst%head
    if (.not.associated(curr)) return ! the list is empty

#ifndef UNSAFE
    if ((.not.present(dt)).and.(.not.present(F))) &
        call ERROR(message='The adaptive time-step calculator(s) `dt_QED_` require a field object `F`.', &
          fix='If `F` is omitted, then a fixed time-step value `dt` must be provided.', &
          in='subroutine list_push')
#endif

    !-----/ MASSIVE-PARTICLE PUSHER /-----!
    if (curr%par%m /= 0.0) then

#ifndef UNSAFE
        if (.not.(present(G).and.present(F))) &
            call ERROR(message='Missing variables required to push particles of nonzero mass.', &
              fix='Provide a force `G` and a field `F` object.', &
              in='subroutine list_push')
#endif

        do while (associated(curr))
            par => curr%par

            if (par%stat == 1) then
              ! check particle freeze condition(s) !
              if (present(t_max)) then; if (par%x(0)>=t_max) then; par%stat=0; goto 11; end if; end if;
              if (present(enrg_min)) then; if (par%p(0)<=enrg_min) then; par%stat=0; goto 11; end if; end if;
              ! set particle-push time-step !
              if (present(dt)) then; par%dt = dt; ! fixed
              else ! adaptive
                maxrate = maxrate_QED_SYNCH(F%E(par%x(1:3),par%x(0)), F%B(par%x(1:3),par%x(0)))
                par%dt = dt_scl/maxrate
                if (par%dt > dt_max) then
                  par%dt = dt_max
#ifndef NOWARN
                  if (.not.suppress_warning(1)) then
                  call WARNING(message="A massive particle's local time-step has exceeded `dt_max`.", &
                    effect='The particle will propagate at `dt_max` for this cycle.', &
                    in='subroutine list_push')
                  print '(A)', '(Subsequent warnings of this kind will be suppressed.)'
                  suppress_warning(1) = .true.
                  end if
#endif
                end if
              end if

              ! push particle !
              call par%push(G,F)
            end if

11          curr => curr%next
        end do

    !-----/ PHOTON PUSHER /-----!
    else

        if (photon_dynamics) then
        do while (associated(curr))
            par => curr%par

            if (par%stat == 1) then
              ! check particle freeze condition(s) !
              if (present(t_max)) then; if (par%x(0)>=t_max) then; par%stat=0; goto 22; end if; end if;
              if (present(enrg_min)) then; if (par%p(0)<=enrg_min) then; par%stat=0; goto 22; end if; end if;
              ! set particle-push time-step !
              if (present(dt)) then; par%dt = dt; ! fixed
              else ! adaptive
                maxrate = maxrate_QED_SYNCH(F%E(par%x(1:3),par%x(0)), F%B(par%x(1:3),par%x(0)))
                par%dt = dt_scl/maxrate
                if (par%dt > dt_max) then
                  par%dt = dt_max
#ifndef NOWARN
                  if (.not.suppress_warning(2)) then
                  call WARNING(message="A photon's local time-step has exceeded `dt_max`.", &
                    effect='The photon will propagate at `dt_max` for this cycle.', &
                    in='subroutine list_push')
                  print '(A)', '(Subsequent warnings of this kind will be suppressed.)'
                  suppress_warning(2) = .true.
                  end if
#endif
                end if
              end if

              ! push particle !
              call par%push()
            end if

22          curr => curr%next
        end do
        end if

    end if

    return
end subroutine list_push

subroutine list_emit(emit_lst,recv_lst,F,L)
    implicit none
    class(plist_t), intent(inout) :: emit_lst ! emitter
    class(plist_t), intent(inout) :: recv_lst ! receiver
    class(field), intent(in) :: F
    class(qmc_table_t), intent(in) :: L

    class(pnode_t), pointer :: curr, next
    class(particle), pointer :: par

    real(num) :: rate, frac, prob
    real(num) :: rand(0:1), enrg(0:1)

    curr => emit_lst%head
    if (.not.associated(curr)) return ! the list is empty

    if (curr%par%m /= 0.0) then

      !-----/ PHOTON EMISSION /-----!
      do while (associated(curr))
        par => curr%par
        if (par%stat == 1) then
          par%qp = Calc_QP_eta(par%x(1:3), par%p(1:3), par%x(0), F)
          rate = Calc_Rate_Synchrotron(par%qp, par%p(0), L)
          prob = rate*par%dt

          call random_number(rand(0))
          call update_diag(par%qp,prob,1)

          if (prob > rand(0)) then
            check1: do
              enrg(0) = Calc_newPhoton_Energy(par%qp,par,F,L)
              if (enrg(0) <= par%p(0)-par%m) exit check1
            end do check1

            if (enrg(0) >= Ep_min) then
              ! create photon
              call recv_lst%create(new_particle(m=0., q=0., id=par%id+1, x=par%x, E=enrg(0), N=par%p))
              ! update lepton energy/momentum (radiation reaction)
              par%p(1:3) = par%p(1:3) - enrg(0)*par%p(1:3)/Norm(par%p(1:3))
              par%p(0) = par%m*Lorentz(par%p(1:3),par%m)
            end if
          end if
        end if
        curr => curr%next
      end do

    else

      !-----/ e-/+ PAIR-PRODUCTION /-----!
      if (photon_dynamics) then
      do while (associated(curr))
        par => curr%par
        next => curr%next

        if (par%p(0) < 2.0_num) &
          par%stat = 0 ! photon energy below pair-creation threshold

        if (par%stat == 1) then
          par%qp = Calc_QP_chi(par%x(1:3), par%p(1:3), par%x(0), F)
          rate = Calc_Rate_PairProduction(par%qp, curr%par%p(0), L)
          prob = rate*par%dt

          call random_number(rand)
          call update_diag(par%qp,prob,2)

          if (prob > rand(0)) then
            ! assign energy fractions to e-/+ pair
            ! energy-conserving scheme, /not/ momentum-conserving
            check2: do
              frac = Calc_newPair_Fraction(par%qp,L)
              enrg = par%p(0)*[frac, (1.-frac)]
              if (all(enrg >= 1.0_num)) exit check2
            end do check2

            ! create e-/+ pair
            call recv_lst%create(new_particle(m=1., q=-1., id=par%id+1, x=par%x, E=enrg(int(2*rand(1))), N=par%p))
            call recv_lst%create(new_particle(m=1., q=+1., id=par%id+1, x=par%x, E=enrg(1-int(2*rand(1))), N=par%p))
            ! add photon to list of annihilated
            call annihilated_photons%create(par)
            ! annihilate photon (in-situ deletion is safe)
            call emit_lst%delete(curr)
            ! mark parent photon for annihilation (deprecated method)
            ! par%stat = -1
          end if
        end if

        curr => next
      end do
      end if

      ! Note: The `next` pointer is used in the e-/+ pair-production section,
      ! but not in the photon emission section, because of the possibility
      ! of photon annihilation (the `delete` subroutine nullifies `curr`).

      ! !-----/ PHOTON ANNIHILATION /-----! (deprecated method)
      ! curr => emit_lst%head
      ! do while (associated(curr))
      !   next => curr%next
      !   if (curr%par%stat == -1) &
      !     call emit_lst%delete(curr)
      !   curr => next
      ! end do

    end if

    return
contains

! update diagnostic variables
subroutine update_diag(qp,prob,i)
    real(num), intent(in) :: qp ! quantum parameter (eta or chi)
    real(num), intent(in) :: prob ! photon or e-/+ emission probability
    integer, intent(in) :: i ! storage index (1/2)
    if (qp < diag%min_qp(i)) diag%min_qp(i) = qp
    if (qp > diag%max_qp(i)) diag%max_qp(i) = qp
    if (prob < diag%min_prob(i)) diag%min_prob(i) = prob
    if (prob > diag%max_prob(i)) diag%max_prob(i) = prob
end subroutine update_diag

end subroutine list_emit

! __________________________________________________________________________________________________
!
! PARTICLE-LIST DIAGNOSTIC ROUTINES
! __________________________________________________________________________________________________
!

integer function count_particles(lst)
    implicit none
    class(plist_t), intent(in) :: lst
    class(pnode_t), pointer :: curr
    count_particles = 0
    curr => lst%head
    do while (associated(curr))
        count_particles = count_particles + 1
        curr => curr%next
    end do
    return
end function count_particles

function exists_particle_wstat(lst,stat) result(ans)
    implicit none
    class(plist_t), intent(in) :: lst
    integer, intent(in) :: stat
    class(pnode_t), pointer :: curr

    logical :: ans
    ans=.false.

    curr => lst%head
    do while (associated(curr))
        if (curr%par%stat == stat) then; ans=.true.; return; end if;
        curr => curr%next
    end do
    return
end function exists_particle_wstat

function get_particle_ids(lst) result(ans)
    implicit none
    class(plist_t), intent(in) :: lst
    class(pnode_t), pointer :: curr
    integer, allocatable :: ans(:)
    integer :: j

    allocate(ans(lst%N))
    if (lst%N == 0) return

    curr => lst%head; j = 0;
    do while (associated(curr))
        j = j+1
        ans(j) = curr%par%id
        curr => curr%next
    end do
    return
end function get_particle_ids

function get_particle_charges(lst) result(ans)
    implicit none
    class(plist_t), intent(in) :: lst
    class(pnode_t), pointer :: curr
    real, allocatable :: ans(:)
    integer :: j

    allocate(ans(lst%N))
    if (lst%N == 0) return

    curr => lst%head; j = 0;
    do while (associated(curr))
        j = j+1
        ans(j) = curr%par%q
        curr => curr%next
    end do
    return
end function get_particle_charges

function get_particle_energies(lst) result(ans)
    implicit none
    class(plist_t), intent(in) :: lst
    class(pnode_t), pointer :: curr
    real(num), allocatable :: ans(:)
    integer :: j

    allocate(ans(lst%N))
    if (lst%N == 0) return

    curr => lst%head; j = 0;
    do while (associated(curr))
        j = j+1
        ans(j) = curr%par%p(0)
        curr => curr%next
    end do
    return
end function get_particle_energies

function get_particle_quantum_parameters(lst) result(ans)
    implicit none
    class(plist_t), intent(in) :: lst
    class(pnode_t), pointer :: curr
    real(num), allocatable :: ans(:)
    integer :: j

    allocate(ans(lst%N))
    if (lst%N == 0) return

    curr => lst%head; j = 0;
    do while (associated(curr))
        j = j+1
        ans(j) = curr%par%qp
        curr => curr%next
    end do
    return
end function get_particle_quantum_parameters

function get_particle_initial_quantities(lst,which) result(ans)
    implicit none
    class(plist_t), intent(in) :: lst
    character(len=*), intent(in) :: which
    class(pnode_t), pointer :: curr
    real(num), allocatable :: ans(:)
    integer :: j

    allocate(ans(lst%N))
    if (lst%N == 0) return

    curr => lst%head; j = 0;

    select case (which)
    case ('time','T','t')
    do while (associated(curr))
        j = j+1
        ans(j) = curr%par%x0_
        curr => curr%next
    end do
    case ('energy','E','e')
    do while (associated(curr))
        j = j+1
        ans(j) = curr%par%p0_
        curr => curr%next
    end do
    case default
        ans = 0.0_num
    end select
    return
end function get_particle_initial_quantities

function get_particle_final_quantities(lst,which,ic) result(ans)
    implicit none
    class(plist_t), intent(in) :: lst
    character(len=*), intent(in) :: which
    integer, intent(in) :: ic
    class(pnode_t), pointer :: curr
    real(num), allocatable :: ans(:)
    integer :: j

    allocate(ans(lst%N))
    if (lst%N == 0) return
    if (ic<0 .or. ic>3) then
        call ERROR("Component index `ic` must be in {0,1,2,3}.", &
            in='function get_particle_final_quantities')
    end if

    curr => lst%head; j = 0;

    select case (which)
    case ('position','X','x')
    do while (associated(curr))
        j = j+1
        ans(j) = curr%par%x(ic)
        curr => curr%next
    end do
    case ('momentum','P','p')
    do while (associated(curr))
        j = j+1
        ans(j) = curr%par%p(ic)
        curr => curr%next
    end do
    case default
        ans = 0.0_num
    end select
    return
end function get_particle_final_quantities

! __________________________________________________________________________________________________
!
! FIELD, FORCE, AND RATE-TABLE MANAGEMENT PROCEDURES
! __________________________________________________________________________________________________
!

subroutine set_field_components(F,Ex,Ey,Ez,Bx,By,Bz)
    implicit none
    class(field), intent(inout) :: F
    procedure(emfunc) :: Ex,Ey,Ez,Bx,By,Bz
    F%Ex=>Ex; F%Ey=>Ey; F%Ez=>Ez;
    F%Bx=>Bx; F%By=>By; F%Bz=>Bz;
    return
end subroutine set_field_components

subroutine unset_field_components(F)
    implicit none
    class(field), intent(inout) :: F
    if (associated(F%Ex)) nullify(F%Ex)
    if (associated(F%Ey)) nullify(F%Ey)
    if (associated(F%Ez)) nullify(F%Ez)
    if (associated(F%Bx)) nullify(F%Bx)
    if (associated(F%By)) nullify(F%By)
    if (associated(F%Bz)) nullify(F%Bz)
    return
end subroutine unset_field_components

function get_field_electric(F,r,t) result(E)
    class(field), intent(in) :: F
    real(num), intent(in) :: r(3), t
    real(num) :: E(3)
    E = [F%Ex(r,t), F%Ey(r,t), F%Ez(r,t)]
end function get_field_electric

function get_field_magnetic(F,r,t) result(B)
    class(field), intent(in) :: F
    real(num), intent(in) :: r(3), t
    real(num) :: B(3)
    B = [F%Bx(r,t), F%By(r,t), F%Bz(r,t)]
end function get_field_magnetic

subroutine set_force_function(G,eval)
    implicit none
    class(force), intent(inout) :: G
    procedure(frfunc), optional :: eval
    G%eval => Force_Lorentz
    if (present(eval)) G%eval => eval
    return
end subroutine set_force_function

subroutine unset_force_function(G)
    implicit none
    class(force), intent(inout) :: G
    if (associated(G%eval)) nullify(G%eval)
    return
end subroutine unset_force_function

subroutine load_qmc_tables(L,path)
    implicit none
    class(qmc_table_t), intent(inout) :: L
    character(len=*), intent(in) :: path
    integer :: j, iu, n(2)
    logical :: exist

    character(len=7), dimension(10) :: filenames
    filenames = [character(len=7) :: 'eta_0i','chi_0ji','Py_0ij','chi_1i',&
        'frac_1j','Pf_1ij','eta_2i','h_2i','chi_3i','g_3i']

    allocate(character(len=len(trim(adjustl(path)))) :: L%path)
    L%path = trim(adjustl(path))

    do j=1,size(filenames)
        ! ensure data file exists
        inquire(file=L%path//'/'//trim(adjustl(filenames(j)))//'.bin', exist=exist)
        if (.not.exist) then
            call ERROR("Could not find '"//trim(adjustl(filenames(j)))//".bin' in path"&
            //NEW_LINE('A')//'  <'//L%path//'>', in='subroutine load_qmc_tables')
        end if

        ! load data file
        open(newunit=iu, file=L%path//'/'//trim(adjustl(filenames(j)))//'.bin', form='unformatted', access='stream')
        select case (j)
        case (1)
            read(iu) n(1); allocate(L%eta0(n(1))); read(iu) L%eta0;
        case (2)
            read(iu) n(1), n(2); allocate(L%chi0(n(1),n(2))); read(iu) L%chi0;
        case (3)
            read(iu) n(1), n(2); allocate(L%Py0(n(1),n(2))); read(iu) L%Py0;
        case (4)
            read(iu) n(1); allocate(L%chi1(n(1))); read(iu) L%chi1;
        case (5)
            read(iu) n(1); allocate(L%frac1(n(1))); read(iu) L%frac1;
        case (6)
            read(iu) n(1), n(2); allocate(L%Pf1(n(1),n(2))); read(iu) L%Pf1;
        case (7)
            read(iu) n(1); allocate(L%eta2(n(1))); read(iu) L%eta2;
        case (8)
            read(iu) n(1); allocate(L%h2(n(1))); read(iu) L%h2;
        case (9)
            read(iu) n(1); allocate(L%chi3(n(1))); read(iu) L%chi3;
        case (10)
            read(iu) n(1); allocate(L%g3(n(1))); read(iu) L%g3;
        end select
        close(iu)
    end do

    ! ensure lookup data are correctly shaped
    if (size(L%eta0)/=size(L%chi0,1) .or. size(L%eta0)/=size(L%Py0,1) .or. size(L%chi0,2)/=size(L%Py0,2)) &
        call ERROR('data shape mismatch (eta_0i.bin, chi_0ji.bin, and Py_0ij.bin).', in='subroutine load_qmc_tables')

    if (size(L%chi1)/=size(L%Pf1,1) .or. size(L%frac1)/=size(L%Pf1,2)) &
        call ERROR('data shape mismatch (chi_1i.bin, frac_1j.bin, and Pf_1ij.bin).', in='subroutine load_qmc_tables')

    if (size(L%eta2)/=size(L%h2)) &
        call ERROR('data shape mismatch (eta_2i.bin and h_2i.bin).', in='subroutine load_qmc_tables')

    if (size(L%chi3)/=size(L%g3)) &
        call ERROR('data shape mismatch (chi_3i.bin and g_3i.bin).', in='subroutine load_qmc_tables')

    return
end subroutine load_qmc_tables

! __________________________________________________________________________________________________
!
! QED RATE & PARAMETER CALCULATION AND INTERPOLATION FUNCTIONS
! __________________________________________________________________________________________________
!

real(num) function Calc_QP_eta(r,p,t,F)
    real(num), intent(in) :: r(3), p(3), t
    class(field), intent(in) :: F
    real(num) :: v(3), x1(3), x2
    ! mass := 1.0 implied throughout
    v = p/Lorentz(p,1.)
    x1 = F%E(r,t) + (v.cross.F%B(r,t))
    x2 = (v.dot.F%E(r,t))**2
    Calc_QP_eta = hbar_*Lorentz(p,1.)*sqrt((x1.dot.x1)-x2)
end function Calc_QP_eta

real(num) function Calc_QP_chi(r,p,t,F)
    real(num), intent(in) :: r(3), p(3), t
    class(field), intent(in) :: F
    real(num) :: n(3), x1(3), x2
    ! Norm(p) = photon energy / m_e / (c^2)
    n = p/Norm(p) ! unit k-vector
    x1 = F%E(r,t) + (n.cross.F%B(r,t))
    x2 = (n.dot.F%E(r,t))**2
    Calc_QP_chi = 0.5*hbar_*Norm(p)*sqrt((x1.dot.x1)-x2)
end function Calc_QP_chi

function Calc_Rate_Synchrotron(eta,enrg,L) result(rate)
    real(num), intent(in) :: eta, enrg
    class(qmc_table_t), intent(in) :: L
    real(num) :: rate, synch, x_(2), y_(2)
    integer :: k

    ! Search quantum-parameter table & interpolate synchrotron function
    if (eta <= L%eta2(1)) then; synch = L%h2(1);
    else if (eta >= L%eta2(size(L%eta2))) then; synch = L%h2(size(L%h2));
    else
        k = BinarySearch(L%eta2, 1, size(L%eta2), eta)
        x_ = L%eta2(k:k+1)
        y_ = L%h2(k:k+1)
        synch = y_(1) + (eta-x_(1))*(y_(2)-y_(1))/(x_(2)-x_(1))
    end if

    ! Calculate synchrotron rate (1/time [sim.u.])
    ! Note: In the rate equation, `enrg` is the Lorentz factor (mass := 1.0 implied).
    rate = sqrt(3.)*cgs%alpha/(2*pi*hbar_)*(eta/enrg)*synch
end function Calc_Rate_Synchrotron

function Calc_Rate_PairProduction(chi,enrg,L) result(rate)
    real(num), intent(in) :: chi, enrg
    class(qmc_table_t), intent(in) :: L
    real(num) :: rate, emiss, x_(2), y_(2)
    integer :: k

    ! Search quantum-parameter table & interpolate emissivity function
    if (chi <= L%chi3(1)) then; emiss = L%g3(1);
    else if (chi >= L%chi3(size(L%chi3))) then; emiss = L%g3(size(L%g3));
    else
        k = BinarySearch(L%chi3, 1, size(L%chi3), chi)
        x_ = L%chi3(k:k+1)
        y_ = L%g3(k:k+1)
        emiss = y_(1) + (chi-x_(1))*(y_(2)-y_(1))/(x_(2)-x_(1))
    end if

    ! Calculate pair-production rate (1/time [sim.u.])
    ! Note: In the rate equation, `enrg` is the photon energy.
    rate = (cgs%alpha/hbar_)*(chi/enrg)*emiss
end function Calc_Rate_PairProduction

function Calc_newPhoton_Energy(eta,par,F,L) result(enrg)
    real(num), intent(in) :: eta ! of generating lepton `par`
    class(particle), intent(in) :: par
    class(field), intent(in) :: F
    class(qmc_table_t), intent(in) :: L
    real(num) :: enrg ! of new photon

    ! interpolation variables
    real(num) :: Py, chi, denom, Py_(2,2), chi_(2,2), eta_(2)
    real(num) :: ku(3), x1(3), x2
    integer :: i, j(2)

    ! Find lepton quantum parameter in table
    if (eta < L%eta0(1)) then; enrg = 0.0_num; return; ! below range
    else if (eta > L%eta0(size(L%eta0))) then; i = size(L%eta0)-1; ! above range
    else
        i = BinarySearch(L%eta0, 1, size(L%eta0), eta)
        if (i == size(L%eta0)) i = i-1
    end if

    ! Generate and search table for Py, the probability
    ! of emitting a chi(i,j) photon by an eta(i) lepton.
    call random_number(Py)
    j(1) = BinarySearch(L%Py0(i,:), 1, size(L%Py0,2), Py)
    j(2) = BinarySearch(L%Py0(i+1,:), 1, size(L%Py0,2), Py)

    ! Interpolate new-photon quantum parameter
    eta_ = L%eta0(i:i+1)

    chi_(1,:) = L%chi0(i, j(1):j(1)+1)
    chi_(2,:) = L%chi0(i+1, j(2):j(2)+1)

    Py_(1,:) = L%Py0(i, j(1):j(1)+1)
    Py_(2,:) = L%Py0(i+1, j(2):j(2)+1)

    chi = (eta_(2)-eta_(1))*Py &
        + (eta-eta_(1))*(chi_(2,1)*Py_(2,2)-chi_(2,2)*Py_(2,1))/(chi_(2,2)-chi_(2,1)) &
        + (eta_(2)-eta)*(chi_(1,1)*Py_(1,2)-chi_(1,2)*Py_(1,1))/(chi_(1,2)-chi_(1,1))

    denom = (eta-eta_(1))*(Py_(2,2)-Py_(2,1))/(chi_(2,2)-chi_(2,1)) &
        + (eta_(2)-eta)*(Py_(1,2)-Py_(1,1))/(chi_(1,2)-chi_(1,1))

    chi = chi/denom

    ! Calculate new-photon energy from `chi`
    ku = par%p(1:3)/Norm(par%p(1:3)) ! photon parallel to lepton

    x1 = F%E(par%x(1:3),par%x(0)) + (ku.cross.F%B(par%x(1:3),par%x(0)))
    x2 = (ku.dot.F%E(par%x(1:3),par%x(0)))**2

    enrg = (2*chi/hbar_)/sqrt((x1.dot.x1)-x2)
end function Calc_newPhoton_Energy

function Calc_newPair_Fraction(chi,L) result(frac)
    real(num), intent(in) :: chi ! of annihilating photon
    class(qmc_table_t), intent(in) :: L
    real(num) :: frac ! of energy split

    ! interpolation variables
    real(num) :: Pf, denom, chi_(2), frac_(2,2), Pf_(2,2)
    integer :: i, j(2)

    ! Find photon quantum parameter in table
    if (chi <= L%chi1(1)) then; i = 1; ! below range; `frac` will be ~1/2.
    else if (chi > L%chi1(size(L%chi1))) then; i = size(L%chi1)-1; ! above range
    else
        i = BinarySearch(L%chi1, 1, size(L%chi1), chi)
        if (i == size(L%chi1)) i = i-1
    end if

    ! Generate and search table for Pf, the probability of a
    ! chi-photon providing pair-energy fractions f and (1-f).
    call random_number(Pf)
    j(1) = BinarySearch(L%Pf1(i,:), 1, size(L%Pf1,2), Pf)
    j(2) = BinarySearch(L%Pf1(i+1,:), 1, size(L%Pf1,2), Pf)

    ! Interpolate pair energy share-fraction
    chi_ = L%chi1(i:i+1)

    frac_(1,:) = L%frac1(j(1):j(1)+1)
    frac_(2,:) = L%frac1(j(2):j(2)+1)

    Pf_(1,:) = L%Pf1(i, j(1):j(1)+1)
    Pf_(2,:) = L%Pf1(i+1, j(2):j(2)+1)

    frac = (chi_(2)-chi_(1))*Pf &
        + (chi-chi_(1))*(frac_(2,1)*Pf_(2,2)-frac_(2,2)*Pf_(2,1))/(frac_(2,2)-frac_(2,1)) &
        + (chi_(2)-chi)*(frac_(1,1)*Pf_(1,2)-frac_(1,2)*Pf_(1,1))/(frac_(1,2)-frac_(1,1))

    denom = (chi-chi_(1))*(Pf_(2,2)-Pf_(2,1))/(frac_(2,2)-frac_(2,1)) &
        + (chi_(2)-chi)*(Pf_(1,2)-Pf_(1,1))/(frac_(1,2)-frac_(1,1))

    frac = frac/denom
end function Calc_newPair_Fraction

! __________________________________________________________________________________________________
!
! UTILITY FUNCTIONS
! __________________________________________________________________________________________________
!

pure real(num) function Norm3(a)
    real(num), intent(in) :: a(3)
    Norm3 = sqrt(sum(a*a))
end function Norm3

pure real(num) function Norm4(a,sign)
    real(num), intent(in) :: a(0:3)
    integer, intent(in) :: sign
    Norm4 = sign*(a(0)**2 - sum(a(1:3)**2))
end function Norm4

pure real(num) function Lorentz(p,m)
    real(num), intent(in) :: p(3)
    real, intent(in) :: m
    Lorentz = sqrt(1.0 + sum(p*p)/(m**2))
end function Lorentz

function Force_Lorentz(m,q,r,p,t,F) result(force)
    real, intent(in) :: m, q
    real(num), intent(in) :: r(3), p(3), t
    class(field), intent(in) :: F
    real(num) :: force(3)
    force = q*F%E(r,t) + (q/m)*(p.cross.F%B(r,t))/Lorentz(p,m)
end function Force_Lorentz

real(num) function maxrate_QED_SYNCH(E,B)
    real(num), intent(in) :: E(3), B(3)
    ! Maximum photon emission rate (frequency units, [sim.u.])
    maxrate_QED_SYNCH = 5.24*sqrt(3.)*cgs%alpha/pi * max(maxval(abs(E)), maxval(abs(B)))
end function maxrate_QED_SYNCH

real(num) function maxrate_QED_PPROD(E,B)
    real(num), intent(in) :: E(3), B(3)
    ! Maximum pair-production rate (frequency units, [sim.u.])
    maxrate_QED_PPROD = 0.2*cgs%alpha * max(maxval(abs(E)), maxval(abs(B)))
end function maxrate_QED_PPROD

real(num) function dt_QED_SYNCH(E,B)
    real(num), intent(in) :: E(3), B(3)
    ! Inverse of maximum photon emission rate (time units, [sim.u.])
    dt_QED_SYNCH = maxrate_QED_SYNCH(E,B)**(-1)
end function dt_QED_SYNCH

real(num) function dt_QED_PPROD(E,B)
    real(num), intent(in) :: E(3), B(3)
    ! Inverse of maximum pair-production rate (time units, [sim.u.])
    dt_QED_PPROD = maxrate_QED_PPROD(E,B)**(-1)
end function dt_QED_PPROD

! __________________________________________________________________________________________________
!
! ERROR / WARNING UTILITIES
! __________________________________________________________________________________________________
!

subroutine ERROR(message, in, fix)
    character(len=*), intent(in) :: message
    character(len=*), intent(in), optional :: in, fix
    if (present(in)) print '(A)', '[In `'//in//'`]'
    print '(A)', 'Error: '//message
    if (present(fix)) print '(A)', 'Fix: '//fix
    ! print '(A)', 'Terminating execution.'
    call exit(1)
end subroutine ERROR

subroutine WARNING(message, in, effect)
    character(len=*), intent(in) :: message
    character(len=*), intent(in), optional :: in, effect
    if (present(in)) print '(A)', '[In `'//in//'`]'
    print '(A)', 'Warning: '//message
    if (present(effect)) print '(A)', 'Effect: '//effect
    return
end subroutine WARNING

end module qed_mcpm
