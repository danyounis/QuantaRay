! QuantaRay (qray)
!
! Main Classical Calculation Routine (CCR)
!
! Author:  D. Younis
!          University of Rochester
!          Department of Physics
!
! Written: 6/21/2020
! Revised: 8/23/2022
!
! DESCRIPTION
!   Propagates classical electron trajectories whose initial conditions are
!   determined by the wavefunction probability current and phase gradient.
!
! INPUT
!   Input deck defining the simulation parameters (refer to qr-stdin.deck)
!   and HDF5 dataset containing virtual detector data (vd_data.h5).
!
! OUTPUT
!   HDF5 dataset containing final electron trajectory data.

program main
    use prec
    use math
    use quantum
    use emfm
    use vdm
    use rochester
    use dataio
    use omp_lib

    implicit none

    ! code revision date
    character(len=8), parameter :: crv = '20220823'
    ! execution date / time
    character(len=8) :: exec_date
    character(len=6) :: exec_time
    ! total computation time
    real, dimension(2) :: cpu_t, wall_t

    ! input deck variables
    character(len=24) :: deck_name
    character(len=9) :: mode
    character(len=3) :: int_method, dtp_unit
    integer :: nt_update, iu

    ! global simulation time parameters
    real(num) :: t_max, dt, dtp

    ! time variables
    real(num), dimension(:), allocatable :: t
    integer :: nt, it, itp, its

    ! em-field variables
    type(emf) :: field
    real(num) :: fldt(2,3)

    ! virtual detector container variables
    real(num), dimension(:,:),   allocatable :: VD_loc, VD_rho, VD_phase
    real(num), dimension(:,:,:), allocatable :: VD_Krt, VD_Jrt
    integer :: N_VD, n, k, j

    ! end detector variables
    type(edet) :: ED

    ! electron variables
    type(particle_electron), allocatable :: electron(:)
    real(num) :: erkc(4,4), epdot(2), elagr(2)
    integer   :: ne = 0, tne, ine

    ! kinetic/potential variables
    real(num) :: ac(3), Z1, Z2

    ! particle tracking variables
    real(num), dimension(:,:), allocatable :: tt_data
    integer, dimension(:), allocatable :: idl, itl, ipt
    integer :: max_tt_per_vd
    real(num) :: wgt_perc

    ! start computation clock
    call cpu_time(cpu_t(1))
    wall_t(1) = OMP_get_wtime()
    ! get date/time of execution
    call date_and_time(date=exec_date, time=exec_time)
    ! utilize all available threads
    call OMP_set_num_threads(OMP_get_max_threads())
    call OMP_set_nested(.false.)

    ! parse input deck
    call get_command_argument(1,deck_name)
    call read_deck

    ! make time domain
    nt = nint(t_max/dt) + 1
    t = linspace(0.0_num, t_max, dt)

    ! initialize em-field
    field%ch1 = 0.0; field%CEP = field%CEP * pi;
    if (mode(1:4) /= '2d1e') field%eps = 0.0
    if (field%profile(1:5) == 'trapz') call field % init_trapz(t)
    if (field%profile == 'sine2') call field % init_sine2(t)
    if (field%profile == 'gauss') call field % init_gauss(t)
    if (field%profile == 'gauss-l') call field % init_gauss_l(t)

    ! aligned 2-electron model
    if (mode(1:4) == '1d2e') then
        field%Ey = field%Ex
        field%Ay = field%Ax
        field%Cy = field%Cx
    end if

    ! set integration time-step
    if (dtp_unit == 'cyc') dtp = dtp*field%T0
    its = nint(dtp/dt)

    ! LOAD VIRTUAL DETECTOR DATA !
    select case (mode(1:4))
    case ('2d1e','1d2e')
        call READ_DATA_CCR(VD_loc, VD_rho, VD_phase, VD_Krt, VD_Jrt, sdim=2)
        N_VD = size(VD_loc,1) ! dim(VD_loc) = (N_VD,2)
    case ('1d1e')
        call READ_DATA_CCR(VD_loc, VD_rho, VD_phase, VD_Krt, VD_Jrt, sdim=1)
        N_VD = size(VD_loc,2) ! dim(VD_loc) = (1,N_VD)
        int_method = 'fwe'
    end select

    ! PARTICLE TRACKING !
    if (mode(6:9) == 'trac') call particle_tracking

    ! static memory allocation of electron trajectories
    tne = N_VD * (field%it_off - field%it_on + 1) ! total no. electrons
    allocate(electron(tne))

    ! dynamic memory allocation of electron trajectories
    ! allocate(temp_electron(ne))
    ! temp_electron(1:ne-1) = electron
    ! call move_alloc(temp_electron, electron)

    ! initialize runtime log
    call init_runtime_log

    ! TIME PROPAGATION BRANCH !
    ! 1D
    if (mode(1:4) == '1d1e') call qrayC1d1e
    ! 2D
    select case (int_method)
    case ('rk4')
        call tprop_rk4
    case ('fwe')
        call tprop_fwe
    end select

    ! record all electron trajectory information
    call remove_bound_trajectories
    call ED % escan(electron, sdim=2)

    ! output end detector data
    open(newunit=iu, file='runtime_ccr.log', status='old', access='append')
    write(iu,'(A)') NEW_LINE('A')//'Writing data...'
    close(iu)

    call WRITE_DATA_CCR(ED)

    ! stop computation clock
    call cpu_time(cpu_t(2))
    wall_t(2) = OMP_get_wtime()

    ! finalize runtime log
    call end_runtime_log
    call exit()

contains

! RUNGE-KUTTA PARTICLE PROPAGATION ROUTINE !
subroutine tprop_rk4
    do it = field%it_on, field%it_off
        if (mod(it,nt_update) == 0) call update_runtime_log

        ! emit batch of electrons
        do n=1,N_VD
            ne = ne + 1
            electron(ne) % x = VD_loc(n,1)
            electron(ne) % y = VD_loc(n,2)
            electron(ne) % px = VD_Krt(1,it,n)
            electron(ne) % py = VD_Krt(2,it,n)
            electron(ne) % phase = VD_phase(n,it)
            electron(ne) % weight = VD_rho(n,it)
            electron(ne) % propagate = .true.
        end do

        ! push electrons in batch
        do itp = it, field%it_off, its
            call cache_field(fldt,itp,its)
            !$OMP PARALLEL DO SHARED(electron) PRIVATE(erkc,elagr)
            do ine = ne-N_VD+1, ne
                ! calculate RK4 coefficients
                call calc_coeff_rk4(electron(ine),fldt,erkc)

                ! calculate lagrangian
                call calc_lagr(electron(ine),fldt(:,1),elagr(1))

                ! dt-propagate electron
                call dt_prop(electron(ine),erkc)

                ! recalculate lagrangian
                call calc_lagr(electron(ine),fldt(:,3),elagr(2))

                ! calculate dt-accumulated phase
                call dt_phase(electron(ine),elagr)
            end do
            !$OMP END PARALLEL DO
        end do

        ! calculate final phase
        !$OMP PARALLEL DO SHARED(electron)
        do ine = ne-N_VD+1, ne
            electron(ine)%phase = electron(ine)%phase &
                - electron(ine)%px*electron(ine)%x - electron(ine)%py*electron(ine)%y
            electron(ine)%propagate = .false.
        end do
        !$OMP END PARALLEL DO
    end do
    return
end subroutine tprop_rk4

! FORWARD-EULER PARTICLE PROPAGATION ROUTINE !
subroutine tprop_fwe
    do it = field%it_on, field%it_off
        if (mod(it,nt_update) == 0) call update_runtime_log

        ! emit batch of electrons
        do n=1,N_VD
            ne = ne + 1
            electron(ne) % x = VD_loc(n,1)
            electron(ne) % y = VD_loc(n,2)
            electron(ne) % px = VD_Krt(1,it,n)
            electron(ne) % py = VD_Krt(2,it,n)
            electron(ne) % phase = VD_phase(n,it)
            electron(ne) % weight = VD_rho(n,it)
            electron(ne) % propagate = .true.
        end do

        ! push electrons in batch
        do itp = it, field%it_off, its
            call cache_field(fldt,itp,its)
            !$OMP PARALLEL DO SHARED(electron) PRIVATE(epdot,elagr)
            do ine = ne-N_VD+1, ne
                ! calculate force components
                select case (mode(1:4))
                case ('2d1e')
                    epdot(1) = -DVroc(Z1,Z2,electron(ine)%x,electron(ine)%y,ac(1),1) + Z2*field%Ex(itp)
                    epdot(2) = -DVroc(Z1,Z2,electron(ine)%x,electron(ine)%y,ac(1),2) + Z2*field%Ey(itp)
                case ('1d2e')
                    epdot(1) = -DVroc(Z1, Z2, electron(ine)%x, 0.0_num, ac(1), 1) + Z2*field%Ex(itp) &
                        - DVroc(Z2, Z2, electron(ine)%x-electron(ine)%y, 0.0_num, ac(3), 1)
                    epdot(2) = -DVroc(Z1, Z2, 0.0_num, electron(ine)%y, ac(2), 2) + Z2*field%Ex(itp) &
                        - DVroc(Z2, Z2, electron(ine)%y-electron(ine)%x, 0.0_num, ac(3), 1)
                end select

                ! calculate lagrangian
                call calc_lagr(electron(ine),fldt(:,1),elagr(1))

                ! dt-propagate electron
                call electron(ine) % apush(epdot,dtp)

                ! recalculate lagrangian
                call calc_lagr(electron(ine),fldt(:,3),elagr(2))

                ! calculate dt-accumulated phase
                call dt_phase(electron(ine),elagr)
            end do
            !$OMP END PARALLEL DO
        end do

        ! calculate final phase
        !$OMP PARALLEL DO SHARED(electron)
        do ine = ne-N_VD+1, ne
            electron(ine)%phase = electron(ine)%phase &
                - electron(ine)%px*electron(ine)%x - electron(ine)%py*electron(ine)%y
            electron(ine)%propagate = .false.
        end do
        !$OMP END PARALLEL DO
    end do
    return
end subroutine tprop_fwe

! MAIN COMPUTATION ROUTINE FOR A ONE-DIMENSIONAL WAVEFUNCTION !
subroutine qrayC1d1e
    do it = field%it_on, field%it_off
        if (mod(it,nt_update) == 0) call update_runtime_log

        ! emit batch of electron trajectories
        do n=1,N_VD
            ne = ne + 1
            electron(ne) % x = VD_loc(1,n)
            electron(ne) % px = VD_Krt(n,it,1)
            electron(ne) % phase = VD_phase(n,it)
            electron(ne) % weight = VD_rho(n,it)
            electron(ne) % propagate = .true.
        end do

        ! push electron trajectories in batch
        do itp = it, field%it_off, its
            !$OMP PARALLEL DO SHARED(electron) PRIVATE(epdot,elagr)
            do ine = ne-N_VD+1, ne
                ! calculate force
                epdot(1) = -DVroc(Z1,Z2,electron(ine)%x,0.0_num,ac(1),1) + Z2*field%Ex(itp)

                ! calculate lagrangian
                elagr(1) = (electron(ine)%px**2)/2.0 - Vroc(Z1,Z2,electron(ine)%x,0.0_num,ac(1)) &
                    + Z2*electron(ine)%x*field%Ex(itp)

                ! dt-propagate electron
                call electron(ine) % apush(epdot(1),dtp)

                ! recalculate lagrangian
                elagr(2) = (electron(ine)%px**2)/2.0 - Vroc(Z1,Z2,electron(ine)%x,0.0_num,ac(1))
                if ((itp+its) < field%it_off) elagr(2)=elagr(2) + Z2*electron(ine)%x*field%Ex(itp+its)

                ! calculate dt-accumulated phase
                electron(ine)%phase = electron(ine)%phase + 0.5*(elagr(1)+elagr(2))*dtp
            end do
            !$OMP END PARALLEL DO
        end do

        ! calculate final phase
        !$OMP PARALLEL DO SHARED(electron)
        do ine = ne-N_VD+1, ne
            electron(ine)%phase = electron(ine)%phase - electron(ine)%px*electron(ine)%x
            electron(ine)%propagate = .false.
        end do
        !$OMP END PARALLEL DO
    end do

    ! record all electron trajectory information
    call ED % escan(electron, sdim=1)

    ! output end detector data
    open(newunit=iu, file='runtime_ccr.log', status='old', access='append')
    write(iu,'(A)') NEW_LINE('A')//'writing data...'
    close(iu)

    call WRITE_DATA_CCR(ED)

    ! stop computation clock
    call cpu_time(cpu_t(2))
    wall_t(2) = OMP_get_wtime()

    ! finalize runtime log
    call end_runtime_log
    call exit()
end subroutine qrayC1d1e

! PARTICLE TRACKING ROUTINE !
subroutine particle_tracking

    ! cache tracking info
    ! tne: total no. electrons
    ! (n,i) = (idl(k),itl(k)): electron k launches from detector n at time-step i
    call ttr_initialize

    ! initialize propagation variables
    allocate(electron(tne))

    !$OMP PARALLEL DO SHARED(electron)
    do k=1,tne
        electron(k) % x = VD_loc(idl(k),1)
        electron(k) % y = VD_loc(idl(k),2)
        electron(k) % px = VD_Krt(1,itl(k),idl(k))
        electron(k) % py = VD_Krt(2,itl(k),idl(k))
        electron(k) % phase = VD_phase(idl(k),itl(k))
        electron(k) % weight = VD_rho(idl(k),itl(k))
        electron(k) % propagate = .true.

        ! initialize data output
        call ttr_store_data(electron(k),k,itl(k),'new')
    end do
    !$OMP END PARALLEL DO

    ! RK4 PROPAGATION with TRACKING !
    !$OMP PARALLEL DO SHARED(electron) PRIVATE(itp,erkc,elagr,fldt)
    do k=1,tne
    do itp = itl(k), field%it_off, its
        ! calculate RK4 coefficients
        call cache_field(fldt,itp,its)
        call calc_coeff_rk4(electron(k),fldt,erkc)

        ! calculate lagrangian
        call calc_lagr(electron(k),fldt(:,1),elagr(1))

        ! dt-propagate electron
        call dt_prop(electron(k),erkc)

        ! recalculate lagrangian
        call calc_lagr(electron(k),fldt(:,3),elagr(2))

        ! calculate dt-accumulated phase
        call dt_phase(electron(k),elagr)

        ! record dynamical variables
        call ttr_store_data(electron(k),k,itp+its,'old')
    end do
    end do
    !$OMP END PARALLEL DO

    ! get partitioning indices
    call ttr_list_partition

    ! fetch data from written files
    call ttr_compile_data('t')
    call ttr_compile_data('x')
    call ttr_compile_data('y')
    call ttr_compile_data('px')
    call ttr_compile_data('py')
    call ttr_compile_data('fz')

    select case (mode(1:4))
    case ('2d1e')
        call ttr_compile_data('en')
    case ('1d2e')
        call ttr_compile_data('en1')
        call ttr_compile_data('en2')
    end select

    ! output trajectory data
    call WRITE_DATA_CCR_TT(tt_data, ipt)
    call exit()
end subroutine particle_tracking

! calculate RK4 propagation coefficients
subroutine calc_coeff_rk4(e, fldt, erkc)
    type(particle_electron), intent(in) :: e
    real(num), intent(in) :: fldt(2,3)
    real(num), intent(out) :: erkc(4,4)

    ! RK4 PROPAGATION !
    !   Let x-dot = f(p) and p-dot = g(t,x).
    !
    ! 4th order Runge-Kutta coefficients
    !   K1x = dt * f(p)
    !   K1p = dt * g(t,x)
    !   K2x = dt * f(p + K1p/2)
    !   K2p = dt * g(t + dt/2, x + K1x/2)
    !   K3x = dt * f(p + K2p/2)
    !   K3p = dt * g(t + dt/2, x + K2x/2)
    !   K4x = dt * f(p + K3p)
    !   K4p = dt * g(t + dt, x + K3x)
    !
    ! dt-advanced position and momentum
    !   x(n+1) = x(n) + (K1x + 2*K2x + 2*K3x + K4x)/6
    !   p(n+1) = p(n) + (K1p + 2*K2p + 2*K3p + K4p)/6

    select case (mode(1:4))
    case ('2d1e')
    erkc(1,1) = dtp*e%px ! K1x
    erkc(1,2) = dtp*e%py ! K1y
    erkc(1,3) = dtp*(Z2*fldt(1,1) - DVroc(Z1,Z2,e%x,e%y,ac(1),1)) ! K1px
    erkc(1,4) = dtp*(Z2*fldt(2,1) - DVroc(Z1,Z2,e%x,e%y,ac(1),2)) ! K1py

    erkc(2,1) = dtp*(e%px + 0.5*erkc(1,3)) ! K2x
    erkc(2,2) = dtp*(e%py + 0.5*erkc(1,4)) ! K2y
    erkc(2,3) = dtp*(Z2*fldt(1,2) - DVroc(Z1,Z2,e%x+0.5*erkc(1,1),e%y+0.5*erkc(1,2),ac(1),1)) ! K2px
    erkc(2,4) = dtp*(Z2*fldt(2,2) - DVroc(Z1,Z2,e%x+0.5*erkc(1,1),e%y+0.5*erkc(1,2),ac(1),2)) ! K2py

    erkc(3,1) = dtp*(e%px + 0.5*erkc(2,3)) ! K3x
    erkc(3,2) = dtp*(e%py + 0.5*erkc(2,4)) ! K3y
    erkc(3,3) = dtp*(Z2*fldt(1,2) - DVroc(Z1,Z2,e%x+0.5*erkc(2,1),e%y+0.5*erkc(2,2),ac(1),1)) ! K3px
    erkc(3,4) = dtp*(Z2*fldt(2,2) - DVroc(Z1,Z2,e%x+0.5*erkc(2,1),e%y+0.5*erkc(2,2),ac(1),2)) ! K3py

    erkc(4,1) = dtp*(e%px + erkc(3,3)) ! K4x
    erkc(4,2) = dtp*(e%py + erkc(3,4)) ! K4y
    erkc(4,3) = dtp*(Z2*fldt(1,3) - DVroc(Z1,Z2,e%x+erkc(3,1),e%y+erkc(3,2),ac(1),1)) ! K4px
    erkc(4,4) = dtp*(Z2*fldt(2,3) - DVroc(Z1,Z2,e%x+erkc(3,1),e%y+erkc(3,2),ac(1),2)) ! K4py

    case ('1d2e')
    erkc(1,1) = dtp*e%px ! K1x1
    erkc(1,2) = dtp*e%py ! K1x2
    erkc(1,3) = dtp*(-DVroc(Z1,Z2,e%x,0.0_num,ac(1),1) + Z2*fldt(1,1) & ! K1p1
        - DVroc(Z2,Z2,e%x-e%y,0.0_num,ac(3),1))
    erkc(1,4) = dtp*(-DVroc(Z1,Z2,0.0_num,e%y,ac(2),2) + Z2*fldt(1,1) & ! K1p2
        - DVroc(Z2,Z2,e%y-e%x,0.0_num,ac(3),1))

    erkc(2,1) = dtp*(e%px + 0.5*erkc(1,3)) ! K2x1
    erkc(2,2) = dtp*(e%py + 0.5*erkc(1,4)) ! K2x2
    erkc(2,3) = dtp*(-DVroc(Z1,Z2,e%x+0.5*erkc(1,1),0.0_num,ac(1),1) + Z2*fldt(1,2) &
        - DVroc(Z2,Z2,e%x+0.5*erkc(1,1)-e%y-0.5*erkc(1,2),0.0_num,ac(3),1)) ! K2p1
    erkc(2,4) = dtp*(-DVroc(Z1,Z2,0.0_num,e%y+0.5*erkc(1,2),ac(2),2) + Z2*fldt(1,2) &
        - DVroc(Z2,Z2,e%y+0.5*erkc(1,2)-e%x-0.5*erkc(1,1),0.0_num,ac(3),1)) ! K2p2

    erkc(3,1) = dtp*(e%px + 0.5*erkc(2,3)) ! K3x1
    erkc(3,2) = dtp*(e%py + 0.5*erkc(2,4)) ! K3x2
    erkc(3,3) = dtp*(-DVroc(Z1,Z2,e%x+0.5*erkc(2,1),0.0_num,ac(1),1) + Z2*fldt(1,2) &
        - DVroc(Z2,Z2,e%x+0.5*erkc(2,1)-e%y-0.5*erkc(2,2),0.0_num,ac(3),1)) ! K3p1
    erkc(3,4) = dtp*(-DVroc(Z1,Z2,0.0_num,e%y+0.5*erkc(2,2),ac(2),2) + Z2*fldt(1,2) &
        - DVroc(Z2,Z2,e%y+0.5*erkc(2,2)-e%x-0.5*erkc(2,1),0.0_num,ac(3),1)) ! K3p2

    erkc(4,1) = dtp*(e%px + erkc(3,3)) ! K4x1
    erkc(4,2) = dtp*(e%py + erkc(3,4)) ! K4x2
    erkc(4,3) = dtp*(-DVroc(Z1,Z2,e%x+erkc(3,1),0.0_num,ac(1),1) + Z2*fldt(1,3) &
        - DVroc(Z2,Z2,e%x+erkc(3,1)-e%y-erkc(3,2),0.0_num,ac(3),1)) ! K4p1
    erkc(4,4) = dtp*(-DVroc(Z1,Z2,0.0_num,e%y+erkc(3,2),ac(2),2) + Z2*fldt(1,3) &
        - DVroc(Z2,Z2,e%y+erkc(3,2)-e%x-erkc(3,1),0.0_num,ac(3),1)) ! K4p2
    end select

    return
end subroutine calc_coeff_rk4

! evaluate instantaneous lagrangian
subroutine calc_lagr(e, fldt, elagr)
    type(particle_electron), intent(in) :: e
    real(num), intent(in) :: fldt(2)
    real(num), intent(out) :: elagr

    select case (mode(1:4))
    case ('2d1e')
        elagr = 0.5*(e%px**2 + e%py**2) &
            - Vroc(Z1,Z2,e%x,e%y,ac(1)) &
            + Z2*e%x*fldt(1) + Z2*e%y*fldt(2)
    case ('1d2e')
        elagr = 0.5*(e%px**2 + e%py**2) &
            - Vroc(Z1, Z2, e%x, 0.0_num, ac(1)) &
            - Vroc(Z1, Z2, 0.0_num, e%y, ac(2)) &
            - Vroc(Z2, Z2, e%x-e%y, 0.0_num, ac(3)) &
            + Z2*(e%x + e%y)*fldt(1)
    end select

    return
end subroutine calc_lagr

! cache field components at t, t+dt/2, and t+dt.
subroutine cache_field(fldt, it0, it1)
    real(num), intent(out) :: fldt(2,3)
    integer, intent(in) :: it0, it1
    fldt = 0.0_num

    select case (mode(1:4))
    case ('2d1e')
        fldt(1,1) = field%Ex(it0)
        fldt(2,1) = field%Ey(it0)

        if ((it0+it1/2) < field%it_off) then
            fldt(1,2) = field%Ex(it0+it1/2)
            fldt(2,2) = field%Ey(it0+it1/2)
        end if

        if ((it0+it1) < field%it_off) then
            fldt(1,3) = field%Ex(it0+it1)
            fldt(2,3) = field%Ey(it0+it1)
        end if

    case ('1d2e')
        fldt(1,1) = field%Ex(it0)
        if ((it0+it1/2) < field%it_off) fldt(1,2) = field%Ex(it0+it1/2)
        if ((it0+it1) < field%it_off) fldt(1,3) = field%Ex(it0+it1)
    end select

    return
end subroutine cache_field

! dt-propagate an electron (RK4 scheme)
subroutine dt_prop(e, erkc)
    type(particle_electron), intent(inout) :: e
    real(num), intent(in) :: erkc(4,4)

    e%x  = e%x  + (erkc(1,1) + 2*erkc(2,1) + 2*erkc(3,1) + erkc(4,1))/6.0
    e%y  = e%y  + (erkc(1,2) + 2*erkc(2,2) + 2*erkc(3,2) + erkc(4,2))/6.0
    e%px = e%px + (erkc(1,3) + 2*erkc(2,3) + 2*erkc(3,3) + erkc(4,3))/6.0
    e%py = e%py + (erkc(1,4) + 2*erkc(2,4) + 2*erkc(3,4) + erkc(4,4))/6.0

    return
end subroutine dt_prop

! dt-accumulate phase from action integral
subroutine dt_phase(e, elagr)
    type(particle_electron), intent(inout) :: e
    real(num), intent(in) :: elagr(2)

    e%phase = e%phase + 0.5*(elagr(1)+elagr(2))*dtp

    return
end subroutine dt_phase

! remove bound trajectories by making them weightless
subroutine remove_bound_trajectories
    logical :: cnd(2)

    select case (mode(1:4))

    case ('2d1e')
    !$OMP PARALLEL DO SHARED(electron) PRIVATE(cnd)
    do ine=1,tne
        cnd(1) = 0.5*(electron(ine)%px**2 + electron(ine)%py**2) &
            + Vroc(Z1,Z2,electron(ine)%x,electron(ine)%y,ac(1)) < 0.0

        if (cnd(1)) electron(ine)%weight = 0.0_num
    end do
    !$OMP END PARALLEL DO

    case ('1d2e')
    !$OMP PARALLEL DO SHARED(electron) PRIVATE(cnd)
    do ine=1,tne
        cnd(1) = (electron(ine)%px**2)/2.0 + Vroc(Z1,Z2,electron(ine)%x,0.0_num,ac(1)) &
            + Vroc(Z2,Z2,electron(ine)%x-electron(ine)%y,0.0_num,ac(3)) < 0.0
        cnd(2) = (electron(ine)%py**2)/2.0 + Vroc(Z1,Z2,0.0_num,electron(ine)%y,ac(2)) &
            + Vroc(Z2,Z2,electron(ine)%y-electron(ine)%x,0.0_num,ac(3)) < 0.0

        if (cnd(1).or.cnd(2)) electron(ine)%weight = 0.0_num
    end do
    !$OMP END PARALLEL DO

    end select

    ! calculate bound/free weight totals
    ED % bfwt = 0.0_num
    do ine=1,tne
        ED%bfwt(2) = ED%bfwt(2) + electron(ine)%weight ! free
    end do
    ED%bfwt(1) = sum(VD_rho) - ED%bfwt(2) ! bound = total - free

    return
end subroutine remove_bound_trajectories

! initialize particle tracking
subroutine ttr_initialize
    real(num) :: wgt_trac, wgt_lims(2)
    integer :: iostat, cnt

    ! initialize tab files
    open(newunit=iu, file='idl.dat', status='new', action='write'); close(iu); ! launching detector indices
    open(newunit=iu, file='itl.dat', status='new', action='write'); close(iu); ! launching time indices

    ! count and record desired particle indices
    ! i.e. those whose weight is +/- 0.5% within wgt_perc of each detector's max weight
    ! max_tt_per_vd: no. tracked particles per detector

    tne = 0 ! initialize total no. electrons

    do n=1,N_VD
        wgt_trac = (wgt_perc*1.d-2)*maxval(VD_rho(n,:)); cnt = 0;
    do it=field%it_on,field%it_off
        wgt_lims(1) = VD_rho(n,it) - 5.d-3*VD_rho(n,it)
        wgt_lims(2) = VD_rho(n,it) + 5.d-3*VD_rho(n,it)

        ! increment total no. trajectories and record launch indices
        if ((wgt_trac >= wgt_lims(1)).and.(wgt_trac <= wgt_lims(2))) then
            tne = tne + 1; cnt = cnt + 1;

            open(newunit=iu, file='idl.dat', status='old', action='write', position='append')
            write(iu,*) n
            close(iu)

            open(newunit=iu, file='itl.dat', status='old', action='write', position='append')
            write(iu,*) it
            close(iu)
        end if

        if (cnt == max_tt_per_vd) exit ! break
    end do
    end do

    ! store desired particle launch indices
    allocate(idl(tne))
    open(newunit=iu, file='idl.dat', status='old', action='read'); k = 1;
    do
        read(iu,*,iostat=iostat) n
        if (iostat > 0) then
            ! unknown error
            print '(A)', 'Error: File I/O illegal data'
            print '(A)', 'Terminating execution.'
            call exit(1)
        else if (iostat < 0) then
            ! done reading, delete file
            close(iu, status='delete')
            exit ! break
        else
            ! store index
            idl(k) = n
            k = k + 1
        end if
    end do

    allocate(itl(tne))
    open(newunit=iu, file='itl.dat', status='old', action='read'); k = 1;
    do
        read(iu,*,iostat=iostat) it
        if (iostat > 0) then
            ! unknown error
            print '(A)', 'Error: File I/O illegal data'
            print '(A)', 'Terminating execution.'
            call exit(1)
        else if (iostat < 0) then
            ! done reading, delete file
            close(iu, status='delete')
            exit ! break
        else
            ! store index
            itl(k) = it
            k = k + 1
        end if
    end do

    return
end subroutine ttr_initialize

! write single-particle dynamical variables to temporary file
subroutine ttr_store_data(e, k, it, status)
    type(particle_electron), intent(in) :: e
    integer, intent(in) :: k, it ! trajectory no., time index
    character(len=3), intent(in) :: status ! new/old dat file
    character(len=9) :: k_str
    real(num) :: temp

    write(k_str,'(I9)') k

    ! append data to temporary k-th trajectory file

    ! time-step
    open(k, file='t-'//trim(adjustl(k_str))//'.dat', status=status, action='write', position='append')
    write(k,*) t(it)
    close(k)

    ! position
    open(k, file='x-'//trim(adjustl(k_str))//'.dat', status=status, action='write', position='append')
    write(k,*) e%x
    close(k)

    open(k, file='y-'//trim(adjustl(k_str))//'.dat', status=status, action='write', position='append')
    write(k,*) e%y
    close(k)

    ! momenta
    open(k, file='px-'//trim(adjustl(k_str))//'.dat', status=status, action='write', position='append')
    write(k,*) e%px
    close(k)

    open(k, file='py-'//trim(adjustl(k_str))//'.dat', status=status, action='write', position='append')
    write(k,*) e%py
    close(k)

    ! phase
    open(k, file='fz-'//trim(adjustl(k_str))//'.dat', status=status, action='write', position='append')
    write(k,*) e%phase
    close(k)

    ! energy
    select case (mode(1:4))
    case ('2d1e')
        temp = 0.5*(e%px**2 + e%py**2) &
            + Vroc(Z1,Z2,e%x,e%y,ac(1))! - Z2*(e%x*field%Ex(it) + e%y*field%Ey(it))

        open(k, file='en-'//trim(adjustl(k_str))//'.dat', status=status, action='write', position='append')
        write(k,*) temp
        close(k)

    case ('1d2e')
        temp = 0.5*e%px**2 + Vroc(Z1,Z2,e%x,0.0_num,ac(1)) &
            + Vroc(Z2,Z2,e%x-e%y,0.0_num,ac(3))! - Z2*e%x*field%Ex(it)

        open(k, file='en1-'//trim(adjustl(k_str))//'.dat', status=status, action='write', position='append')
        write(k,*) temp
        close(k)

        temp = 0.5*e%py**2 + Vroc(Z1,Z2,e%y,0.0_num,ac(2)) &
            + Vroc(Z2,Z2,e%y-e%x,0.0_num,ac(3))! - Z2*e%y*field%Ex(it)

        open(k, file='en2-'//trim(adjustl(k_str))//'.dat', status=status, action='write', position='append')
        write(k,*) temp
        close(k)
    end select

    return
end subroutine ttr_store_data

! record indices that partition trajectories in master list
subroutine ttr_list_partition
    character(len=9) :: k_str

    ! data i/o variables
    real(num) :: datln
    integer :: iostat

    allocate(ipt(tne)); ipt(1) = 0;

    ! iterate through time files
    do k=1,tne
        write(k_str,'(I9)') k
        open(newunit=iu, file='t-'//trim(adjustl(k_str))//'.dat', status='old', action='read')
        if (k > 1) ipt(k) = ipt(k-1)

        do
            read(iu,*,iostat=iostat) datln
            if (iostat > 0) then
                ! unknown error
                print '(A)', 'Error: File I/O illegal data'
                print '(A)', 'Terminating execution.'
                call exit(1)
            else if (iostat < 0) then
                ! done reading
                close(iu, status='keep')
                exit ! break
            else
                ! increment partition counter
                ipt(k) = ipt(k) + 1
            end if
        end do
    end do

    return
end subroutine ttr_list_partition

! compile trajectory data files into one and then transfer them to master variable (tt_data)
subroutine ttr_compile_data(var)
    character(len=*), intent(in) :: var
    character(len=9) :: k_str

    ! data i/o variables
    real(num) :: datln
    integer :: iostat, iu(2)

    ! compile single-trajectory files into one
    open(newunit=iu(1), file=var//'.dat', status='new', action='write')

    do k=1,tne
        write(k_str,'(I9)') k
        open(newunit=iu(2), file=var//'-'//trim(adjustl(k_str))//'.dat', status='old', action='read')

        ! iterate through current trajectory file, copying data
        do
            read(iu(2),*,iostat=iostat) datln
            if (iostat > 0) then
                ! unknown error
                print '(A)', 'Error: File I/O illegal data'
                print '(A)', 'Terminating execution.'
                call exit(1)
            else if (iostat < 0) then
                ! done reading, delete single-trajectory file
                close(iu(2), status='delete')
                exit ! break
            else
                ! transfer data line
                write(iu(1),*) datln
            end if
        end do
    end do
    close(iu(1))

    ! ipt(tne) = total no. of trajectory recordings
    if (.not.allocated(tt_data)) then
        allocate( tt_data(8,ipt(tne)) ); tt_data = 0.0_num;
    end if

    ! store and delete compiled dat files
    open(newunit=iu(1), file=var//'.dat', status='old', action='read')
    do j=1,ipt(tne)
        select case (var)
        case ('t')
            read(iu(1),*) tt_data(1,j)
        case ('x')
            read(iu(1),*) tt_data(2,j)
        case ('y')
            read(iu(1),*) tt_data(3,j)
        case ('px')
            read(iu(1),*) tt_data(4,j)
        case ('py')
            read(iu(1),*) tt_data(5,j)
        case ('fz')
            read(iu(1),*) tt_data(6,j)
        case ('en','en1')
            read(iu(1),*) tt_data(7,j)
        case ('en2')
            read(iu(1),*) tt_data(8,j)
        end select
    end do
    close(iu(1), status='delete')

    return
end subroutine ttr_compile_data

! read input deck
subroutine read_deck
    character(len=128) :: label, temp_profile
    open(newunit=iu, file=deck_name, status='old', action='read')
    read(iu,*) label, mode
    read(iu,*) ! skip ccr
    read(iu,*)
    read(iu,*) ! global simulation time
    read(iu,*) label, t_max
    read(iu,*) label, dt
    read(iu,*)
    read(iu,*) ! particle propagation
    read(iu,*) label, int_method
    read(iu,*) label, dtp, dtp_unit
    read(iu,*)
    read(iu,*) ! particle tracking
    read(iu,*) label, max_tt_per_vd
    read(iu,*) label, wgt_perc
    read(iu,*)
    read(iu,*) ! wavefunction domain
    read(iu,*) ! skip x range
    read(iu,*) ! skip y range
    read(iu,*) ! skip nx
    read(iu,*) ! skip ny
    read(iu,*) ! skip dm
    read(iu,*)
    read(iu,*) ! atomic properties
    read(iu,*) label, Z1
    read(iu,*) label, Z2
    read(iu,*) ! skip Pa
    read(iu,*) ! skip En_str
    read(iu,*) ! skip tau
    read(iu,*) label, ac(1), ac(2), ac(3)
    read(iu,*)
    read(iu,*) ! field parameters
    read(iu,*) label, temp_profile
    allocate( character(len=len(trim(adjustl(temp_profile)))) :: field%profile )
    field%profile = trim(adjustl(temp_profile))
    read(iu,*) label, field % E0
    read(iu,*) label, field % omg0
    read(iu,*) label, field % eps
    if (field%profile(1:5) == 'trapz') read(iu,*) label, field % Ncyc_pl
    if (field%profile(1:5) == 'sine2') read(iu,*) label, field % Tp
    if (field%profile(1:5) == 'gauss') read(iu,*) label, field % Tfwhm
    read(iu,*) label, field % CEP
    if (field%profile(1:5) == 'trapz') read(iu,*) label, field % t_on
    if (field%profile(1:5) == 'sine2') read(iu,*) label, field % t_on
    if (field%profile(1:5) == 'gauss') read(iu,*) label, field % Tpk
    read(iu,*)
    read(iu,*) ! virtual detectors
    read(iu,*) ! skip Nv
    read(iu,*) ! skip R0
    read(iu,*) ! skip unwrap_phase
    read(iu,*)
    read(iu,*) ! t-SURFF distribution
    read(iu,*) ! skip enable
    read(iu,*) ! skip interp
    read(iu,*) ! skip dti
    read(iu,*) ! skip Ns
    read(iu,*) ! skip R0
    read(iu,*) ! skip kx
    read(iu,*) ! skip ky
    read(iu,*) ! skip nk
    read(iu,*)
    read(iu,*) ! qm grid output
    read(iu,*) ! skip dt_save
    read(iu,*)
    read(iu,*) ! qm proj output
    read(iu,*) ! skip dt_save
    read(iu,*)
    read(iu,*) ! update runtime log every N time steps
    read(iu,*) label, nt_update
    read(iu,*)
    read(iu,*) ! echo this input deck?
    read(iu,*) ! skip choice
    close(iu)
    return
end subroutine read_deck

! initialize runtime log
subroutine init_runtime_log
    open(newunit=iu, file='runtime_ccr.log', status='replace')
    write(iu,'(A,I2)') 'QuantaRay version: ccr'//mode(1:4)//'-'//crv//'_x',8*num
    write(iu,'(A)') 'executed: '//exec_date//' '//exec_time
    write(iu,'(A,I2)') 'OMP threads: ', OMP_get_max_threads()
    write(iu,'(A)') NEW_LINE('A')//'integration method:'
    select case (int_method)
    case ('rk4')
        write(iu,'(A,F4.2,A)') 'Runge-Kutta, dtp = ',dtp,' au'
    case ('fwe')
        write(iu,'(A,F4.2,A)') 'Forward-Euler, dtp = ',dtp,' au'
    end select
    write(iu,'(A,I9)') NEW_LINE('A')//'No. electron trajectories: ', tne
    write(iu,'(A)') NEW_LINE('A')//'Starting propagation...'
    write(iu,'(A,F9.2,A,I6,A)') 't = ', t(field%it_on), ' au    nt: ', field%it_on-1, '    0.0%'
    close(iu)
    return
end subroutine init_runtime_log

! update runtime log
subroutine update_runtime_log
    open(newunit=iu, file='runtime_ccr.log', status='old', access='append')
    write(iu,'(A,F9.2,A,I6,A,F5.1,A)') &
    't = ', t(it)+dt, ' au    nt: ', it,'  ',1e2*(it-field%it_on)/(field%it_off-field%it_on),'%'
    close(iu)
    return
end subroutine update_runtime_log

! finalize runtime log
subroutine end_runtime_log
    open(newunit=iu, file='runtime_ccr.log', status='replace')
    write(iu,'(A,I2)') 'QuantaRay version: ccr'//mode(1:4)//'-'//crv//'_x',8*num
    write(iu,'(A)') 'executed: '//exec_date//' '//exec_time
    write(iu,'(A,I2)') 'OMP threads: ', OMP_get_max_threads()
    write(iu,'(A)') NEW_LINE('A')//'Propagation complete.'
    write(iu,'(A,F9.2,A)') 'total execution time:', (cpu_t(2)-cpu_t(1))/3.60d+3, ' core-hrs'
    write(iu,'(F9.2,A)') (wall_t(2)-wall_t(1))/3.60d+3, ' wall clock time (hrs)'
    close(iu)
    return
end subroutine end_runtime_log

end program main
