! QuantaRay (qray)
!
! Main Quantum Calculation Routine (QCR)
!
! Author:  D. Younis
!          University of Rochester
!          Department of Physics
!
! Written: 6/2/2020
! Revised: 8/23/2022
!
! DESCRIPTION
!   Numerically solves the time-dependent Schrodinger equation for an atomic
!   electron in an electromagnetic field and records information using numerical
!   detectors situated around the nucleus.
!
! INPUT
!   Refer to the input deck template file: stdin.deck.
!
! OUTPUT
!   HDF5 datasets containing wavefunction and virtual detector data.
!
! REFERENCES
!   [1] D. Bauer, Computational Strong-Field Quantum Dynamics.
!   [2] B. Feuerstein and U. Thumm, J. Phys. B 36, 707-716 (2003).
!   [3] X. Wang, J. Tian, and J. H. Eberly, J. Phys. B 51, 084002 (2018).
!   [4] S. L. Haan, P. S. Wheeler, R. Panfili, and J. H. Eberly, Phys. Rev. A 66, 061402(R) (2002).
!   [5] C. A. Moyer, Am. J. Phys. 72, 3 (2004).

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
    ! current working directory
    character(len=255) :: CWD
    ! execution date / time
    character(len=8) :: exec_date
    character(len=6) :: exec_time
    ! total computation time
    real, dimension(2) :: cpu_t, wall_t

    ! input deck variables
    character(len=24) :: deck_name
    character(len=6) :: En_str
    character(len=9) :: mode
    character(len=4) :: pp_method
    logical :: run_ccr, echo_input
    integer :: nt_update

    ! global simulation time parameters
    real(num) :: t_max, tau, dt

    ! wavefunction domain parameters
    real(num) :: x_lim(2), y_lim(2), r_abs, d_msk

    ! grid variables
    real(num), dimension(:), allocatable :: t, x, y, px, py
    real(num) :: dr(2), dp(2)
    integer :: nr(2), nt, ntau, it, ix, iy, iu, j

    ! wavefunction variables
    type(SchrodingerWavefunction1D) :: wavefn1
    type(SchrodingerWavefunction2D) :: wavefn2

    ! em-field variables
    type(emf) :: field

    ! t-SURFF variables
    type(tSURFF2D) :: S

    ! virtual detector variables
    type(vdet), dimension(:),     allocatable :: VD
    real(num),  dimension(:,:,:), allocatable :: IDM
    logical :: unwrap_VD_phase
    integer :: N_VD, n
    real :: R0

    ! kinetic/potential variables
    integer :: Pa
    real(num) :: ac(3), Z1, Z2, En
    real(num),    dimension(:,:), allocatable :: Tk, Vp0, fmsk
    complex(num), dimension(:,:), allocatable :: Vp, Va, psi_e

    ! output info for qm data
    real(num) :: dtg(2)
    integer :: ig(2,3)
    character(len=3) :: go_unit(2)
    ! output containers - full grid / projection data
    complex(num), dimension(:,:,:), allocatable :: psif_out, psip_out
    real(num),    dimension(:,:,:), allocatable :: Vp_out
    real(num),    dimension(:),     allocatable :: tf_out, tp_out

    ! start computation clock
    call cpu_time(cpu_t(1))
    wall_t(1) = OMP_get_wtime()
    ! get current working directory
    call get_environment_variable('PWD',CWD)
    ! get date/time of execution
    call date_and_time(date=exec_date, time=exec_time)
    ! utilize all available threads
    call OMP_set_num_threads(OMP_get_max_threads())
    call OMP_set_nested(.false.)

    ! parse input deck
    call get_command_argument(1,deck_name)
    call read_deck

    nt = nint(t_max/dt) + 1
    ntau = nint(tau/dt) + 1

    ! make time domain
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

    ! set pre-propagation (field-free) method
    ! for finding the desired electronic state
    if (trim(adjustl(En_str)) == 'ground') then
        pp_method = 'imag'
        Pa = 0
    else
        pp_method = 'eigs'
        read(En_str,*) En
    end if

    if ((mode == '2d1e').or.(mode == '2d1e-trac')) mode = '2d1e-cnn'
    if ((mode == '1d2e').or.(mode == '1d2e-trac')) mode = '1d2e-cnn'
    if ((mode == '1d1e').or.(mode == '1d1e-trac')) mode = '1d1e-cnn'
    if (mode(1:4) == '1d1e') call qrayQ1d1e

    ! construct wavefunction domain
    allocate(x(nr(1)), y(nr(2)), px(nr(1)), py(nr(2)))
    x = linspace(x_lim(1), x_lim(2), nr(1)); dr(1) = x(2) - x(1);
    y = linspace(y_lim(1), y_lim(2), nr(2)); dr(2) = y(2) - y(1);
    px = 2.*pi*fft_freq(nr(1),dr(1),.false.); dp(1) = px(2) - px(1);
    py = 2.*pi*fft_freq(nr(2),dr(2),.false.); dp(2) = py(2) - py(1);

    ! allocate variables
    allocate(Tk(nr(1),nr(2)), Vp(nr(1),nr(2)), Vp0(nr(1),nr(2)))
    if (run_ccr) then
        allocate(VD(N_VD), IDM(nr(1),nr(2),12))
    else
        allocate(VD(1), IDM(1,1,1))
    end if

    ! initialize wavefunction
    call wavefn2 % init_vars(nr,nt)
    call wavefn2 % init_form(Pa,dr)
    if (mode(6:8) == 'cnn') &
    call wavefn2 % make_cnn_mats(nr,dr)

    ! construct kinetic energy operator
    ForAll(ix=1:nr(1),iy=1:nr(2)) Tk(ix,iy) = (px(ix)**2 + py(iy)**2)/2.0

    ! construct atomic + masking function
    select case (mode(1:4))
    case ('2d1e')
        ForAll(ix=1:nr(1),iy=1:nr(2)) &
            Vp0(ix,iy) = Vroc(Z1, Z2, x(ix), y(iy), ac(1))
    case ('1d2e')
        ForAll(ix=1:nr(1),iy=1:nr(2)) &
            Vp0(ix,iy) = Vroc(Z1, Z2, x(ix), 0.0_num, ac(1)) & ! electron-nucleus
                       + Vroc(Z1, Z2, y(iy), 0.0_num, ac(2)) &
                       + Vroc(Z2, Z2, x(ix)-y(iy), 0.0_num, ac(3)) ! electron-electron
    end select
    Vp = Vp0
    call mk_mask_function

    ! initialize virtual detectors
    if (run_ccr) then
    !$OMP PARALLEL DO SHARED(VD)
    do n=1,N_VD
        call VD(n) % init(R0, N_VD, x, y, n, nt, sdim=2)
    end do
    !$OMP END PARALLEL DO
    end if

    ! initialize runtime log
    call init_runtime_log

    ! ATOMIC STATE PREPARATION !
    select case (pp_method)
    case ('imag') ! ground state
        do it=1,ntau
            if (mod(it,nt_update) == 0) call update_runtime_log('pre')
            ! advance wavefunction in imaginary time
            select case (mode(6:8))
            case ('cnn')
                call wavefn2 % propagate_cnn(cmplx(Vp0,kind=num), dt, -i)
            case ('fft')
                call wavefn2 % propagate_fft(Tk, cmplx(Vp0,kind=num), dr, dt, -i)
            end select
            ! normalize
            wavefn2%psi = wavefn2%psi/sqrt(dr(1)*dr(2)*trapz(abs(wavefn2%psi)**2))
        end do

    case ('eigs') ! excited state
        allocate(psi_e(nr(1),nr(2))); psi_e = 0.0_num;

        do it=1,ntau
            if (mod(it,nt_update) == 0) call update_runtime_log('pre')
            ! advance wavefunction
            select case (mode(6:8))
            case ('cnn')
                call wavefn2 % propagate_cnn(cmplx(Vp0,kind=num), dt, cmplx(1.0))
            case ('fft')
                call wavefn2 % propagate_fft(Tk, cmplx(Vp0,kind=num), dr, dt, cmplx(1.0))
            end select
            ! integrate with Hanning window
            psi_e = psi_e + wavefn2%psi*winHann((it-1)*dt,tau)*exp(i*En*(it-1)*dt)*dt
        end do

        ! set normalized eigenstate
        wavefn2%psi = psi_e/sqrt(dr(1)*dr(2)*trapz(abs(psi_e)**2))
        deallocate(psi_e)
    end select

    ! compute initial wavefunction gradient
    if (run_ccr) then
    do j=1,2
        call deriv(wavefn2%psi,dr(j),j,wavefn2%grad_psi(:,:,j))
    end do
    end if

    ! calculate time-dependent wavefunction quantities
    wavefn2 % norm(1) = dr(1)*dr(2)*trapz(abs(wavefn2%psi)**2)

    select case (mode(6:8))
    case ('cnn')
        wavefn2 % energy(1) = calc_energy(wavefn2, real(Vp,num), dr)
    case ('fft')
        wavefn2 % energy(1) = calc_energy(wavefn2, Tk, real(Vp,num), dr, dp)
    end select

    ! set output info for qm data
    call set_qm_data_output_info

    allocate(psif_out(nr(1),nr(2),ig(1,1)), Vp_out(nr(1),nr(2),ig(1,1)), tf_out(ig(1,1)))
    allocate(psip_out(ig(2,1),nr(1),2), tp_out(ig(2,1)))

    ! save initial qm data
    call store_qm_grid_data(1)
    call store_qm_proj_data(1)
    tf_out(1) = t(1)
    tp_out(1) = t(1)

    ! output relaxation energy
    open(newunit=iu, file='runtime_qcr.log', status='old', access='append')
    write(iu,'(A)') NEW_LINE('A')//'Atomic state preparation complete.'
    write(iu,'(A,F6.3,A)') 'En = ', wavefn2%energy(1), ' au'//NEW_LINE('A')
    write(iu,'(A)') 'Starting propagation...'//NEW_LINE('A')
    write(iu,'(A,F9.2,A,I6,A)') 't = ', t(1), ' au    nt: ', 1, '    0.0%'
    close(iu)

    ! echo input deck
    if (echo_input) call echo_deck

    ! initialize t-SURFF calculation
    call S%init(x,y,dt)
    call S%dt_step(wavefn2,(/field%Ax(1),field%Ay(1)/),(/field%Cx(1),field%Cy(1)/),x,y,t(1))

    ! FULL ATOM + FIELD PROPAGATION !
    do it=1,nt
        if (mod(it,nt_update) == 0) call update_runtime_log('ord')

        ! update total potential
        ForAll (ix=1:nr(1),iy=1:nr(2)) &
            Vp(ix,iy) = Vp0(ix,iy) - Z2*x(ix)*field%Ex(it) - Z2*y(iy)*field%Ey(it)

        ! advance wavefunction
        select case (mode(6:8))
        case ('cnn')
            call wavefn2 % propagate_cnn(Vp, dt, cmplx(1.0))
        case ('fft')
            call wavefn2 % propagate_fft(Tk, Vp, dr, dt, cmplx(1.0))
        end select

        wavefn2%psi = fmsk*wavefn2%psi
        wavefn2%norm(it) = dr(1)*dr(2)*trapz(abs(wavefn2%psi)**2)

        select case (mode(6:8))
        case ('cnn')
            wavefn2 % energy(it) = calc_energy(wavefn2, real(Vp,num), dr)
        case ('fft')
            wavefn2 % energy(it) = calc_energy(wavefn2, Tk, real(Vp,num), dr, dp)
        end select

        if (run_ccr) then
            ! compute wavefunction phase
            wavefn2 % phase = atan2(aimag(wavefn2%psi),real(wavefn2%psi))
            ! compute wavefunction gradient
            do j=1,2
                call deriv(wavefn2%psi,dr(j),j,wavefn2%grad_psi(:,:,j))
            end do
            ! compute interpolation data matrix
            call POP_IDM
            ! trigger each virtual detector
            !$OMP PARALLEL DO SHARED(VD)
            do n=1,N_VD
                call VD(n) % trigger(wavefn2, IDM, x, y, it)
            end do
            !$OMP END PARALLEL DO
        end if

        ! t-SURFF calculation
        if (mod(it,S%iti) == 0) &
        call S%dt_step(wavefn2,(/field%Ax(it),field%Ay(it)/),(/field%Cx(it),field%Cy(it)/),x,y,t(it))

        if (mod(it,ig(1,2)) == 0) call store_qm_grid_data(it)
        if (mod(it,ig(2,2)) == 0) call store_qm_proj_data(it)
    end do

    ! unwrap virtual detector phase information
    if (run_ccr.and.unwrap_VD_phase) then
        !$OMP PARALLEL DO SHARED(VD)
        do n=1,N_VD
            call phase_unwrap(VD(n)%phase)
        end do
        !$OMP END PARALLEL DO
    end if

    ! output data
    open(newunit=iu, file='runtime_qcr.log', status='old', access='append')
    write(iu,'(A)') NEW_LINE('A')//'Writing data...'
    close(iu)

    call WRITE_DATA_QCR_2D(t, x, y, VD, field, wavefn2, &
        tf_out, psif_out, tp_out, psip_out, Vp_out, S, fmsk, run_ccr)

    ! stop computation clock
    call cpu_time(cpu_t(2))
    wall_t(2) = OMP_get_wtime()

    ! finalize runtime log
    call end_runtime_log
    ! free memory
    call cleanup
    ! execute the CCR
    if (run_ccr) &
    call execute_command_line(trim(adjustl(CWD))//'/qrayC '//deck_name, wait=.true.)
    ! terminate qray
    call exit()

contains

! MAIN COMPUTATION ROUTINE FOR A ONE-DIMENSIONAL WAVEFUNCTION !
subroutine qrayQ1d1e
    ! construct wavefunction domain
    allocate(x(nr(1)), px(nr(1)))
    x = linspace(x_lim(1), x_lim(2), nr(1)); dr(1) = x(2) - x(1);
    px = 2.*pi*fft_freq(nr(1),dr(1),.false.); dp(1) = px(2) - px(1);
    N_VD = 2

    ! allocate variables
    allocate(Tk(1,nr(1)), Vp(1,nr(1)), Vp0(1,nr(1)), VD(N_VD))

    ! initialize wavefunction
    call wavefn1 % init_vars(nr(1),nt)
    call wavefn1 % init_form(Pa,dr(1))
    if (mode(6:8) == 'cnn') &
    call wavefn1 % make_cnn_mats(nr(1),dr(1))

    ! construct kinetic energy operator
    Tk(1,:) = (px**2)/2.0

    ! construct atomic + masking function
    ForAll(ix=1:nr(1)) Vp0(1,ix) = Vroc(Z1, Z2, x(ix), 0.0_num, ac(1))
    Vp(1,:) = Vp0(1,:)
    call mk_mask_function

    ! initialize virtual detectors
    if (run_ccr) then
    !$OMP PARALLEL DO SHARED(VD)
    do n=1,N_VD
        call VD(n) % init(R0, N_VD, x, zeros(2), n, nt, sdim=1)
    end do
    !$OMP END PARALLEL DO
    end if

    ! initialize runtime log
    call init_runtime_log

    ! ATOMIC STATE PREPARATION !
    select case (pp_method)
    case ('imag') ! ground state
        do it=1,ntau
            if (mod(it,nt_update) == 0) call update_runtime_log('pre')
            ! advance wavefunction in imaginary time
            select case (mode(6:8))
            case ('cnn')
                call wavefn1 % propagate_cnn(cmplx(Vp0(1,:),kind=num), dt, -i)
            case ('fft')
                call wavefn1 % propagate_fft(Tk(1,:), cmplx(Vp0(1,:),kind=num), dr(1), dt, -i)
            end select
            ! normalize
            wavefn1%psi = wavefn1%psi/sqrt(dr(1)*trapz(abs(wavefn1%psi)**2))
        end do

    case ('eigs') ! excited state
        allocate(psi_e(1,nr(1))); psi_e = 0.0_num;

        do it=1,ntau
            if (mod(it,nt_update) == 0) call update_runtime_log('pre')
            ! advance wavefunction
            select case (mode(6:8))
            case ('cnn')
                call wavefn1 % propagate_cnn(cmplx(Vp0(1,:),kind=num), dt, cmplx(1.0))
            case ('fft')
                call wavefn1 % propagate_fft(Tk(1,:), cmplx(Vp0(1,:),kind=num), dr(1), dt, cmplx(1.0))
            end select
            ! integrate with Hanning window
            psi_e(1,:) = psi_e(1,:) + wavefn1%psi*winHann((it-1)*dt,tau)*exp(i*En*(it-1)*dt)*dt
        end do

        ! set normalized eigenstate
        wavefn1%psi = psi_e(1,:)/sqrt(dr(1)*trapz(abs(psi_e(1,:))**2))
        deallocate(psi_e)
    end select

    ! compute initial wavefunction gradient
    if (run_ccr) call deriv(wavefn1%psi,dr(1),wavefn1%grad_psi)

    ! calculate time-dependent wavefunction quantities
    wavefn1 % norm(1) = dr(1)*trapz(abs(wavefn1%psi)**2)

    select case (mode(6:8))
    case ('cnn')
        wavefn1 % energy(1) = calc_energy(wavefn1, real(Vp(1,:),num), dr(1))
    case ('fft')
        wavefn1 % energy(1) = calc_energy(wavefn1, Tk(1,:), real(Vp(1,:),num), dr(1), dp(1))
    end select

    ! set output info for qm data
    call set_qm_data_output_info

    allocate(psif_out(ig(1,1),nr(1),1), Vp_out(ig(1,1),nr(1),1), tf_out(ig(1,1)))

    ! save initial qm data
    call store_qm_grid_data(1)
    tf_out(1) = t(1)

    ! output relaxation energy
    open(newunit=iu, file='runtime_qcr.log', status='old', access='append')
    write(iu,'(A)') NEW_LINE('A')//'Atomic state preparation complete.'
    write(iu,'(A,F6.3,A)') 'En = ', wavefn1%energy(1), ' au'//NEW_LINE('A')
    write(iu,'(A)') 'Starting propagation...'//NEW_LINE('A')
    write(iu,'(A,F9.2,A,I6,A)') 't = ', t(1), ' au    nt: ', 1, '    0.0%'
    close(iu)

    ! echo input deck
    if (echo_input) call echo_deck

    ! FULL ATOM + FIELD PROPAGATION !
    do it=1,nt
        if (mod(it,nt_update) == 0) call update_runtime_log('ord')

        ! update total potential
        Vp(1,:) = Vp0(1,:) - Z2*x*field%Ex(it)

        ! advance wavefunction
        select case (mode(6:8))
        case ('cnn')
            call wavefn1 % propagate_cnn(Vp(1,:), dt, cmplx(1.0))
        case ('fft')
            call wavefn1 % propagate_fft(Tk(1,:), Vp(1,:), dr(1), dt, cmplx(1.0))
        end select

        wavefn1%psi = fmsk(1,:)*wavefn1%psi
        wavefn1%norm(it) = dr(1)*trapz(abs(wavefn1%psi)**2)

        select case (mode(6:8))
        case ('cnn')
            wavefn1 % energy(it) = calc_energy(wavefn1, real(Vp(1,:),num), dr(1))
        case ('fft')
            wavefn1 % energy(it) = calc_energy(wavefn1, Tk(1,:), real(Vp(1,:),num), dr(1), dp(1))
        end select

        if (run_ccr) then
            ! compute wavefunction gradient
            call deriv(wavefn1%psi,dr(1),wavefn1%grad_psi)
            ! trigger each virtual detector
            !$OMP PARALLEL DO SHARED(VD)
            do n=1,N_VD
                call VD(n) % trigger(wavefn1, x, it)
            end do
            !$OMP END PARALLEL DO
        end if

        if (mod(it,ig(1,2)) == 0) call store_qm_grid_data(it)
    end do

    ! unwrap virtual detector phase information
    if (run_ccr.and.unwrap_VD_phase) then
        !$OMP PARALLEL DO SHARED(VD)
        do n=1,N_VD
            call phase_unwrap(VD(n)%phase)
        end do
        !$OMP END PARALLEL DO
    end if

    ! output data
    open(newunit=iu, file='runtime_qcr.log', status='old', access='append')
    write(iu,'(A)') NEW_LINE('A')//'Writing data...'
    close(iu)

    call WRITE_DATA_QCR_1D(t, x, VD, field, &
        wavefn1, tf_out, psif_out(:,:,1), Vp_out(:,:,1), fmsk(1,:), run_ccr)

    ! stop computation clock
    call cpu_time(cpu_t(2))
    wall_t(2) = OMP_get_wtime()

    ! finalize runtime log
    call end_runtime_log
    ! free memory
    call cleanup
    ! execute the CCR
    if (run_ccr) &
    call execute_command_line(trim(adjustl(CWD))//'/qrayC '//deck_name, wait=.true.)
    ! terminate qray
    call exit()
end subroutine qrayQ1d1e

! construct a quartic cylindrical absorbing boundary of a specified radius (r_abs)
subroutine mk_abs_potential

    ! select case (mode(1:4))
    ! case ('2d1e','1d2e')
    !     allocate(Va(nr(1),nr(2))); Va = 0.0_num;
    !     do ix=1,nr(1)
    !     do iy=1,nr(2)
    !         if (sqrt(x(ix)**2 + y(iy)**2) > r_abs) &
    !             Va(ix,iy) = -i*abs(sqrt(x(ix)**2 + y(iy)**2) - r_abs)**4
    !     end do
    !     end do
    ! case ('1d1e')
    !     allocate(Va(1,nr(1))); Va = 0.0_num;
    !     do ix=1,nr(1)
    !         if (abs(x(ix)) > r_abs) &
    !             Va(1,ix) = -i*abs(abs(x(ix)) - r_abs)**4
    !     end do
    ! end select

    return
end subroutine mk_abs_potential

! construct a cos-1/8 masking function of a specified distance (d_msk)
subroutine mk_mask_function

    select case (mode(1:4))
    case ('2d1e','1d2e')
        allocate(fmsk(nr(1),nr(2))); fmsk = 1.0_num;
        do ix=1,nr(1)
        do iy=1,nr(2)
            if (abs(y(iy)) > y(nr(2))-d_msk) &
                fmsk(ix,iy) = cos(pi*(abs(y(iy))-y(nr(2))+d_msk)/2.0/d_msk)**(1.0/8.0_num)
        end do
        end do
        fmsk = fmsk * rotM90(fmsk)
    case ('1d1e')
        allocate(fmsk(1,nr(1))); fmsk = 1.0_num;
        do ix=1,nr(1)
            if (abs(x(ix)) > x(nr(1))-d_msk) &
                fmsk(1,ix) = cos(pi*(abs(x(ix))-x(nr(1))+d_msk)/2.0/d_msk)**(1.0/8.0_num)
        end do
    end select

    return
end subroutine mk_mask_function

! populate the data matrix used for VD interpolation
subroutine POP_IDM
    select case (mode(1:4))
    case ('2d1e','1d2e')
        ! IDM(..,1): rho - wavefunction square modulus
        ! IDM(..,2): rho,1 = d(rho)/dr1
        ! IDM(..,3): rho,2 = d(rho)/dr2
        ! IDM(..,4): rho,12 = d2(rho)/dr1,dr2
        IDM(:,:,1) = abs(wavefn2%psi)**2
        call deriv(IDM(:,:,1),dr(1),1,IDM(:,:,2))
        call deriv(IDM(:,:,1),dr(2),2,IDM(:,:,3))
        call deriv(IDM(:,:,2),dr(2),2,IDM(:,:,4))
        ! IDM(..,5): J1 - probability current density along x1
        ! IDM(..,6): J1,1 = d(J1)/dr1
        ! IDM(..,7): J1,2 = d(J1)/dr2
        ! IDM(..,8): J1,12 = d2(J1)/dr1,dr2
        IDM(:,:,5) = real((i/2.0) * &
            (wavefn2%psi*conjg(wavefn2%grad_psi(:,:,1)) - conjg(wavefn2%psi)*wavefn2%grad_psi(:,:,1)), num)
        call deriv(IDM(:,:,5),dr(1),1,IDM(:,:,6))
        call deriv(IDM(:,:,5),dr(2),2,IDM(:,:,7))
        call deriv(IDM(:,:,6),dr(2),2,IDM(:,:,8))
        ! IDM(..,9):  J2 - probability current density along x2
        ! IDM(..,10): J2,1 = d(J2)/dr1
        ! IDM(..,11): J2,2 = d(J2)/dr2
        ! IDM(..,12): J2,12 = d2(J2)/dr1,dr2
        IDM(:,:,9) = real((i/2.0) * &
            (wavefn2%psi*conjg(wavefn2%grad_psi(:,:,2)) - conjg(wavefn2%psi)*wavefn2%grad_psi(:,:,2)), num)
        call deriv(IDM(:,:,9),dr(1),1,IDM(:,:,10))
        call deriv(IDM(:,:,9),dr(2),2,IDM(:,:,11))
        call deriv(IDM(:,:,10),dr(2),2,IDM(:,:,12))
    end select
    return
end subroutine POP_IDM

! set save parameters for qm data
subroutine set_qm_data_output_info
    integer :: k

    ! ig(n,:) = (total, frequency, current index) of saves

    ! ig(1,:) - full grid outputs
    select case (go_unit(1))
    case default
        ig(1,2) = nint(dtg(1)/dt)
    case ('au')
        ig(1,2) = nint(dtg(1)/dt)
    case ('cyc')
        ig(1,2) = nint(dtg(1)*field%T0/dt)
    end select

    ! ig(2,:) - wavefunction projection outputs
    if ((dtg(2)==0.0_num).or.(mode(1:4)=='1d1e')) then
        ig(2,:) = (/1, huge(1_4), 1/) ! nullify
    else
        select case (go_unit(2))
        case default
            ig(2,2) = nint(dtg(2)/dt)
        case ('au')
            ig(2,2) = nint(dtg(2)/dt)
        case ('cyc')
            ig(2,2) = nint(dtg(2)*field%T0/dt)
        end select
    end if

    ig(1,1) = 1; ig(2,1) = 1;
    ig(1,3) = 1; ig(2,3) = 1;

    do k=1,nt
        if (mod(k,ig(1,2)) == 0) ig(1,1) = ig(1,1) + 1
        if (mod(k,ig(2,2)) == 0) ig(2,1) = ig(2,1) + 1
    end do

    return
end subroutine set_qm_data_output_info

! store qm grid data for time step k
subroutine store_qm_grid_data(k)
    integer, intent(in) :: k

    complex(num), allocatable :: psif_dyn(:,:,:)
    real(num), allocatable :: Vp_dyn(:,:,:)
    real(num), allocatable :: tf_dyn(:)

    ! catch overflow from save indexing mismatch
    ! re-allocate _out arrays to avoid segfault
    if (ig(1,3) > ig(1,1)) then
    select case (mode(1:4))
    case ('2d1e','1d2e')
        allocate(psif_dyn(nr(1),nr(2),ig(1,3)), &
            Vp_dyn(nr(1),nr(2),ig(1,3)), &
            tf_dyn(ig(1,3)))

        psif_dyn(:,:,1:ig(1,1)) = psif_out
        Vp_dyn(:,:,1:ig(1,1)) = Vp_out
        tf_dyn(1:ig(1,1)) = tf_out

        call move_alloc(psif_dyn, psif_out)
        call move_alloc(Vp_dyn, Vp_out)
        call move_alloc(tf_dyn, tf_out)
        ig(1,1) = ig(1,1) + 1

    case ('1d1e')
        allocate(psif_dyn(ig(1,3),nr(1),1), &
            Vp_dyn(ig(1,3),nr(1),1), &
            tf_dyn(ig(1,3)))

        psif_dyn(1:ig(1,1),:,:) = psif_out
        Vp_dyn(1:ig(1,1),:,:) = Vp_out
        tf_dyn(1:ig(1,1)) = tf_out

        call move_alloc(psif_dyn, psif_out)
        call move_alloc(Vp_dyn, Vp_out)
        call move_alloc(tf_dyn, tf_out)
        ig(1,1) = ig(1,1) + 1
    end select
    end if

    select case (mode(1:4))
    case ('2d1e','1d2e')
        psif_out(:,:,ig(1,3)) = wavefn2 % psi
        Vp_out(:,:,ig(1,3)) = real(Vp,num)
    case ('1d1e')
        psif_out(ig(1,3),:,1) = wavefn1 % psi
        Vp_out(ig(1,3),:,1) = real(Vp(1,:),num)
    end select
    tf_out(ig(1,3)) = t(k) + dt
    ig(1,3) = ig(1,3) + 1

    return
end subroutine store_qm_grid_data

! store qm projection data for time step k
subroutine store_qm_proj_data(k)
    integer, intent(in) :: k

    complex(num), allocatable :: psip_dyn(:,:,:)
    real(num), allocatable :: tp_dyn(:)

    if (ig(2,3) > ig(2,1)) then
        allocate(psip_dyn(ig(2,3),nr(1),2), tp_dyn(ig(2,3)))

        psip_dyn(1:ig(2,1),:,:) = psip_out
        tp_dyn(1:ig(2,1)) = tp_out

        call move_alloc(psip_dyn, psip_out)
        call move_alloc(tp_dyn, tp_out)
        ig(2,1) = ig(2,1) + 1
    end if

    psip_out(ig(2,3),:,1) = dr(1)*trapz(wavefn2%psi,1)
    psip_out(ig(2,3),:,2) = dr(2)*trapz(wavefn2%psi,2)
    tp_out(ig(2,3)) = t(k) + dt
    ig(2,3) = ig(2,3) + 1

    return
end subroutine store_qm_proj_data

! read input deck
subroutine read_deck
    character(len=128) :: label, temp_profile
    open(newunit=iu, file=deck_name, status='old', action='read')
    read(iu,*) label, mode
    read(iu,*) label, run_ccr
    read(iu,*)
    read(iu,*) ! global simulation time
    read(iu,*) label, t_max
    read(iu,*) label, dt
    read(iu,*)
    read(iu,*) ! particle propagation
    read(iu,*) ! skip int
    read(iu,*) ! skip dtp
    read(iu,*)
    read(iu,*) ! trajectory tracking
    read(iu,*) ! skip no
    read(iu,*) ! skip wgt
    read(iu,*)
    read(iu,*) ! wavefunction domain
    read(iu,*) label, x_lim(1), x_lim(2)
    read(iu,*) label, y_lim(1), y_lim(2)
    read(iu,*) label, nr(1)
    read(iu,*) label, nr(2)
    read(iu,*) label, d_msk
    read(iu,*)
    read(iu,*) ! atomic properties
    read(iu,*) label, Z1
    read(iu,*) label, Z2
    read(iu,*) label, Pa
    read(iu,*) label, En_str
    read(iu,*) label, tau
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
    read(iu,*) label, N_VD
    read(iu,*) label, R0
    read(iu,*) label, unwrap_VD_phase
    read(iu,*)
    read(iu,*) ! t-SURFF distribution
    read(iu,*) label, S%enable
    read(iu,*) label, S%interp
    read(iu,*) label, S%dti
    read(iu,*) label, S%Ns
    read(iu,*) label, S%R0
    read(iu,*) label, S%kx_lim(1), S%kx_lim(2)
    read(iu,*) label, S%ky_lim(1), S%ky_lim(2)
    read(iu,*) label, S%nk(1), S%nk(2)
    read(iu,*)
    read(iu,*) ! qm grid output
    read(iu,*) label, dtg(1), go_unit(1)
    read(iu,*)
    read(iu,*) ! qm proj output
    read(iu,*) label, dtg(2), go_unit(2)
    read(iu,*)
    read(iu,*) ! update runtime log every N time steps
    read(iu,*) label, nt_update
    read(iu,*)
    read(iu,*) ! echo this input deck?
    read(iu,*) label, echo_input
    close(iu)
    return
end subroutine read_deck

! echo simulation parameters to text file
subroutine echo_deck
    open(newunit=iu, file='input-echo.txt', status='replace')
    write(iu,'(A,I2)') 'QuantaRay version: qcr'//trim(adjustl(mode))//'-'//crv//'_x',8*num
    write(iu,'(A)') 'executed: '//exec_date//' '//exec_time
    write(iu,'(A,I2)') 'OMP threads: ', OMP_get_max_threads()
    write(iu,'(A)')
    write(iu,'(A)') 'global simulation time'
    write(iu,'(A,F7.2,A)') 't  = ', t(nt), ' au'
    write(iu,'(A,F5.3,A)') 'dt = ', dt, ' au'
    write(iu,'(A,I6,A)') 'nt = ', nt, ' pts'//NEW_LINE('A')
    write(iu,'(A)') 'wavefunction domain'
    write(iu,'(A,I4,A,I3,A)') 'x: [',nint(x_lim(1)),', ',nint(x_lim(2)),'] au'
    if (mode(1:4) /= '1d1e') &
    write(iu,'(A,I4,A,I3,A)') 'y: [',nint(y_lim(1)),', ',nint(y_lim(2)),'] au'
    write(iu,'(A,E9.3,A)') 'dx = ', dr(1), ' au'
    if (mode(1:4) /= '1d1e') &
    write(iu,'(A,E9.3,A)') 'dy = ', dr(2), ' au'
    write(iu,'(A,I4,A)') 'nx = ', nr(1), ' pts'
    if (mode(1:4) /= '1d1e') &
    write(iu,'(A,I4,A)') 'ny = ', nr(2), ' pts'
    write(iu,'(A,F4.1,A)') 'dm = ', d_msk, ' au'//NEW_LINE('A')
    write(iu,'(A)') 'atomic properties'
    write(iu,'(A,F5.2)') 'Z1 = ', Z1
    write(iu,'(A,F5.2)') 'Z2 = ', Z2
    write(iu,'(A,I2)') 'Pa = ', Pa
    write(iu,'(A)') 'En = '//trim(adjustl(En_str))//' au'
    write(iu,'(A,F6.2,A)') 'tau = ', tau, ' au'
    write(iu,'(A,F4.2,A,F4.2,A,F4.2,A)') 'aps = ',ac(1),' ',ac(2),' ',ac(3),' au'//NEW_LINE('A')
    write(iu,'(A)') 'field parameters'
    write(iu,'(A)') 'type: '//field%profile
    write(iu,'(A,F6.4,A)') 'E0   = ', field%E0, ' au'
    write(iu,'(A,F6.4,A)') 'omg0 = ', field%omg0, ' au'
    if (mode(1:4) == '2d1e') &
    write(iu,'(A,F4.2)') 'eps  = ', field%eps
    select case (field%profile(1:5))
    case ('trapz')
        write(iu,'(A,F4.1,A)') 'time = ', field%Ncyc_pl, ' cycles (plateau)'
    case ('sine2')
        write(iu,'(A,F7.2,A)') 'time = ', field%Tp, ' au (total)'
    case ('gauss')
        write(iu,'(A,F7.2,A)') 'time = ', field%Tfwhm, ' au (intensity FWHM)'
    end select
    write(iu,'(A,F5.3,A)') 'CEP  = ', field%CEP/pi, ' pi'
    select case (field%profile(1:5))
    case default
        write(iu,'(A,F4.1,A)') 't_on = ', field%t_on, ' au'//NEW_LINE('A')
    case ('gauss')
        write(iu,'(A,F7.2,A)') 'Tpk = ', field%Tpk, ' au'//NEW_LINE('A')
    end select
    if (run_ccr) then
        write(iu,'(A)') 'virtual detectors'
        write(iu,'(A,I4)') 'Nv: ', N_VD
        write(iu,'(A,F5.1,A)') 'R0: ', R0, ' au'//NEW_LINE('A')
    end if
    if (S%enable) then
        write(iu,'(A)') 't-SURFF distribution'
    select case (S%interp)
    case ('nn')
        write(iu,'(A)') 'int: nearest-neighbor'
    case ('b3')
        write(iu,'(A)') 'int: bicubic'
    case default
        write(iu,'(A)') 'int: nearest-neighbor'
    end select
        write(iu,'(A,F4.2)') 'dti: ', S%dti
        write(iu,'(A,I4)') 'Ns: ', S%Ns
        write(iu,'(A,F5.1)') 'R0: ', S%R0
        write(iu,'(A,I3,A,I2,A)') 'kx: [',nint(S%kx_lim(1)),', ',nint(S%kx_lim(2)),'] au'
        write(iu,'(A,I3,A,I2,A)') 'ky: [',nint(S%ky_lim(1)),', ',nint(S%ky_lim(2)),'] au'
        write(iu,'(A,I3,A,I3,A)') 'nk: [',S%nk(1),', ',S%nk(2),'] pts'//NEW_LINE('A')
    end if
    write(iu,'(A)') 'qm grid output'
    write(iu,'(A,E9.3,A)') 'dt save: ', dtg(1), ' '//go_unit(1)
    write(iu,'(A,I6)') 'nt saves: ', ig(1,1)
    if (mode(1:4) /= '1d1e') then
        write(iu,'(A)') NEW_LINE('A')//'qm proj output'
        write(iu,'(A,E9.3,A)') 'dt save: ', dtg(2), ' '//go_unit(2)
        write(iu,'(A,I6)') 'nt saves: ', ig(2,1)
    end if
    close(iu)
    return
end subroutine echo_deck

! deallocate variables
subroutine cleanup
    if (allocated(t)) deallocate(t)
    if (allocated(x)) deallocate(x)
    if (allocated(y)) deallocate(y)
    if (allocated(px)) deallocate(px)
    if (allocated(py)) deallocate(py)

    if (allocated(VD)) deallocate(VD)
    if (allocated(IDM)) deallocate(IDM)

    if (allocated(Tk)) deallocate(Tk)
    if (allocated(Vp0)) deallocate(Vp0)
    if (allocated(fmsk)) deallocate(fmsk)

    if (allocated(Vp)) deallocate(Vp)
    if (allocated(Va)) deallocate(Va)
    if (allocated(psi_e)) deallocate(psi_e)

    if (allocated(psif_out)) deallocate(psif_out)
    if (allocated(psip_out)) deallocate(psip_out)
    if (allocated(Vp_out)) deallocate(Vp_out)
    if (allocated(tf_out)) deallocate(tf_out)
    if (allocated(tp_out)) deallocate(tp_out)

    call S % destroy()
    call wavefn1 % destroy()
    call wavefn2 % destroy()

    return
end subroutine cleanup

! initialize runtime log
subroutine init_runtime_log
    open(newunit=iu, file='runtime_qcr.log', status='replace')
    write(iu,'(A,I2)') 'QuantaRay version: qcr'//trim(adjustl(mode))//'-'//crv//'_x',8*num
    write(iu,'(A)') 'executed: '//exec_date//' '//exec_time
    write(iu,'(A,I2)') 'OMP threads: ', OMP_get_max_threads()
    write(iu,'(A)')
    write(iu,'(A)') 'Preparing atomic state...'//NEW_LINE('A')
    select case (mode(1:4))
    case ('2d1e','1d1e')
        write(iu,'(A,F4.2,A)') 'aps = ', ac(1), ' au'
    case ('1d2e')
        write(iu,'(A,F4.2,A,F4.2,A,F4.2,A)') 'aps = (',ac(1),', ',ac(2),', ',ac(3),') au'
    end select
    write(iu,'(A)') 'En = '//trim(adjustl(En_str))//' au'
    write(iu,'(A,I2)') 'Pa = ', Pa
    write(iu,'(A,F9.2,A,I6,A)') NEW_LINE('A')//'t = ', t(1), ' au    nt: ', 1, '    0.0%'
    close(iu)
    return
end subroutine init_runtime_log

! update runtime log
subroutine update_runtime_log(stage)
    character(len=3), intent(in) :: stage
    open(newunit=iu, file='runtime_qcr.log', status='old', access='append')
    select case (stage)
    case ('pre') ! pre-propagation
        write(iu,'(A,F9.2,A,I6,A,F5.1,A)') 't = ',it*dt,' au    nt: ',it,'  ',1e2*it/ntau,'%'
    case ('ord') ! full propagation
        write(iu,'(A,F9.2,A,I6,A,F5.1,A)') 't = ',t(it)+dt,' au    nt: ',it,'  ',1e2*it/nt,'%'
    end select
    close(iu)
    return
end subroutine update_runtime_log

! finalize runtime log
subroutine end_runtime_log
    open(newunit=iu, file='runtime_qcr.log', status='replace')
    write(iu,'(A,I2)') 'QuantaRay version: qcr'//trim(adjustl(mode))//'-'//crv//'_x',8*num
    write(iu,'(A)') 'executed: '//exec_date//' '//exec_time
    write(iu,'(A,I2)') 'OMP threads: ', OMP_get_max_threads()
    write(iu,'(A)') NEW_LINE('A')//'Simulation complete.'
    write(iu,'(A,F9.2,A)') 'total execution time:', (cpu_t(2)-cpu_t(1))/3.60d+3, ' core-hrs'
    write(iu,'(F9.2,A)') (wall_t(2)-wall_t(1))/3.60d+3, ' wall clock time (hrs)'
    close(iu)
    return
end subroutine end_runtime_log

end program main
