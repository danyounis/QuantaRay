! Rabi-split (main)
!
! Author:  D. Younis
!          University of Rochester
!          Department of Physics
!
! Written: 1/26/2021
! Revised: 8/23/2022

program main
    use prec
    use math
    use quantum
    use emfm
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
    integer :: nt_update

    ! global simulation time parameters
    real(num) :: t_max, tau, dt, dtau

    ! wavefunction domain parameters
    real(num) :: r_max, r_abs, r_msk, dr
    integer :: n_init, l_init, l_max

    ! grid variables
    real(num), dimension(:), allocatable :: t, r
    integer :: nt, nr, ntau, it, ir, iu, l, j

    ! wavefunction variables
    type(SchrodingerWavefunction1DR) :: wavefn, wavefn_temp
    complex(num), allocatable :: phi_e(:), eigs_b(:,:), eigs_f(:,:)
    real(num), allocatable :: fmsk(:)
    real(num) :: En

    ! em-field variables
    type(emf) :: pump, probe

    ! photo-electron spectrum variables
    real(num) :: e_lim(2), de
    real(num), allocatable :: E(:), W(:,:)
    integer :: winop_ord, npe_bins

    ! eigenstate projection variables
    real(num), allocatable :: Pbt(:,:), Pft(:,:), Ec(:)
    integer, allocatable :: lc(:)
    integer :: nc

    ! data output info
    real(num) :: dtg
    integer :: ig(3)
    character(len=3) :: go_unit
    ! output containers
    complex(num), dimension(:,:,:), allocatable :: phif_out
    real(num),    dimension(:),     allocatable :: tf_out

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

    nt = nint(t_max/dt) + 1
    ntau = nint(tau/dtau) + 1

    ! construct wavefunction domain
    t = linspace(0.0_num, t_max, dt)
    r = linspace(dr, r_max, dr); nr = size(r);

    ! make photo-electron energy bins
    E = linspace(e_lim(1), e_lim(2), de)
    npe_bins = size(E)

    ! initialize wavefunction
    call wavefn % init_vars(nr,nt,l_max)
    call wavefn % init_prop(r,dr)
    call wavefn % init_form(dr)
    allocate(eigs_b(l_max+1,nr))

    ! initialize pump-probe fields
    pump%eps = 0.0; pump%ch1 = 0.0;
    pump%CEP = pump%CEP * pi
    if (pump%profile(1:5) == 'trapz') call pump % init_trapz(t)
    if (pump%profile == 'sine2') call pump % init_sine2(t)
    if (pump%profile == 'gauss') call pump % init_gauss(t)
    if (pump%profile == 'gauss-l') call pump % init_gauss_l(t)

    probe%eps = 0.0; probe%ch1 = 0.0;
    probe%CEP = probe%CEP * pi
    if (probe%profile(1:5) == 'trapz') call probe % init_trapz(t)
    if (probe%profile == 'sine2') call probe % init_sine2(t)
    if (probe%profile == 'gauss') call probe % init_gauss(t)
    if (probe%profile == 'gauss-l') call probe % init_gauss_l(t)

    ! construct absorbing potential
    call mk_abs_potential
    ! construct masking function
    call mk_mask_function

    ! initialize runtime log
    call init_runtime_log

    ! PREPARE INITIAL ATOMIC STATE !
    call wavefn % prep_atom(l_init, ntau, dr, dtau, pure=.false.)
    do l=0,l_max
        eigs_b(l+1,:) = wavefn%phi(l+1,:)
        if (l /= l_init) wavefn%phi(l+1,:) = 0.0_num
    end do

    if (n_init /= (l_init+1)) then
        allocate(phi_e(nr)); phi_e = 0.0_num;
        En = E_hydrogen(n_init,l_init)

        call wavefn % init_form(dr)

        do l=0,l_max
            if (l /= l_init) wavefn%phi(l+1,:) = 0.0_num
        end do

        do it=1,ntau
            call wavefn % propagate(dtau)
            phi_e = phi_e &
                + wavefn%phi(l_init+1,:)*winHann((it-1)*dtau,tau)*exp(i*En*(it-1)*dtau)*dtau
        end do

        wavefn%phi(l_init+1,:) = phi_e/sqrt(dr*trapz(abs(phi_e)**2))
        deallocate(phi_e)
    end if

    ! cache continuum states of interest
    if (nc > 0) then
        allocate(eigs_f(nc,nr)); eigs_f = 0.0_num;
        allocate(Pft(nc,nt)); Pft = 0.0_num;

        wavefn_temp = wavefn
        call wavefn_temp % init_form(dr)

        do it=1,ntau
            call wavefn_temp % propagate(dtau)
            do j=1,nc
                eigs_f(j,:) = eigs_f(j,:) &
                    + wavefn_temp%phi(lc(j)+1,:)*winHann((it-1)*dtau,tau)*exp(i*Ec(j)*(it-1)*dtau)*dtau
            end do
        end do

        do j=1,nc
            ! normalize
            eigs_f(j,:) = eigs_f(j,:)/sqrt(dr*trapz(abs(eigs_f(j,:))**2))
            ! calculate state energies, over-writing input targets
            wavefn_temp % phi = 0.0_num
            wavefn_temp % phi(lc(j)+1,:) = eigs_f(j,:)
            Ec(j) = calc_energy(wavefn_temp,dr)
        end do

        call wavefn_temp % destroy()
    end if

    wavefn % norm(1) = sum(dr*trapz(abs(wavefn%phi)**2, ax=1))
    wavefn % energy(1) = calc_energy(wavefn,dr)

    ! eigenstate projection
    allocate(Pbt(l_max+1,nt))
    !$OMP PARALLEL DO DEFAULT(SHARED)
    do l=0,l_max
        Pbt(l+1,1) = abs(dr*trapz(conjg(eigs_b(l+1,:))*wavefn%phi(l+1,:)))**2
    end do
    !$OMP END PARALLEL DO

    if (nc > 0) then
        !$OMP PARALLEL DO DEFAULT(SHARED)
        do j=1,nc
            Pft(j,1) = abs(dr*trapz(conjg(eigs_f(j,:))*wavefn%phi(lc(j)+1,:)))**2
        end do
        !$OMP END PARALLEL DO
    end if

    ! set data output info
    call set_data_output_info

    allocate(phif_out(l_max+1,nr,ig(1)), tf_out(ig(1)), W(ig(1),npe_bins))

    ! store initial data
    call store_data(1)
    tf_out(1) = t(1)

    ! output relaxation energy
    open(newunit=iu, file='runtime.log', status='old', access='append')
    write(iu,'(A)') 'Atomic state preparation complete.'
    write(iu,'(A,F6.3,A)') 'En = ', wavefn%energy(1), ' au'//NEW_LINE('A')
    if (nc > 0) then
        write(iu,'(A)') 'Continuum state energies.'
        do j=1,nc
            write(iu,'(A,I1,A,E11.5,A)') 'Ec(',j,') = ',Ec(j),' au'
        end do
        write(iu,'(A)')
    end if
    write(iu,'(A)') 'Starting propagation...'//NEW_LINE('A')
    write(iu,'(A,F9.2,A,I7,A)') 't = ', t(1), ' au    nt: ', 1, '    0.0%'
    close(iu)

    ! FULL ATOM + PUMP-PROBE FIELD PROPAGATION !
    do it=2,nt
        if (mod(it,nt_update) == 0) call update_runtime_log

        ! advance wavefunction
        call wavefn % propagate(pump%Ax(it) + probe%Ax(it), r, dt)

        ! mask wavefunction in l-r space
        wavefn%phi(l_max+1,:) = 0.0_num
        do l=0,l_max
            wavefn%phi(l+1,:) = fmsk*wavefn%phi(l+1,:)
        end do

        wavefn % norm(it) = sum(dr*trapz(abs(wavefn%phi)**2, ax=1))
        wavefn % energy(it) = calc_energy(wavefn,dr)

        ! eigenstate projection
        !$OMP PARALLEL DO DEFAULT(SHARED)
        do l=0,l_max
            Pbt(l+1,it) = abs(dr*trapz(conjg(eigs_b(l+1,:))*wavefn%phi(l+1,:)))**2
        end do
        !$OMP END PARALLEL DO

        if (nc > 0) then
            !$OMP PARALLEL DO DEFAULT(SHARED)
            do j=1,nc
                Pft(j,it) = abs(dr*trapz(conjg(eigs_f(j,:))*wavefn%phi(lc(j)+1,:)))**2
            end do
            !$OMP END PARALLEL DO
        end if

        ! store grid data
        if (mod(it,ig(2)) == 0) call store_data(it)
    end do

    ! normalize differential-energy probability density
    W = W*real(winop_ord,num)*dsin(pi/2.0/real(winop_ord,num))/(pi*de)

    ! nullify unused variables
    if (nc == 0) then
        allocate(Pft(1,1),Ec(1),lc(1))
        Pft = 0.0_num; Ec = 0.0_num; lc = 0;
    end if

    ! output data
    open(newunit=iu, file='runtime.log', status='old', access='append')
    write(iu,'(A)') NEW_LINE('A')//'Writing data...'
    close(iu)

    call WRITE_DATA(t, r, pump, probe, wavefn, tf_out, &
        phif_out, E, W, Pbt, Pft, Ec, lc)

    ! stop computation clock
    call cpu_time(cpu_t(2))
    wall_t(2) = OMP_get_wtime()

    ! finalize runtime log
    call end_runtime_log
    call exit()

contains

! construct a quartic absorbing boundary of a specified radius (r_abs)
subroutine mk_abs_potential

    ! do ir=1,nr
    !     if (r(ir) > r(nr)-r_abs) &
    !         wavefn%Va(ir) = -i*abs(r(ir)-r(nr)+r_abs)**4
    ! end do

    return
end subroutine mk_abs_potential

! construct a cos-1/8 masking function of a specified radius (r_msk)
subroutine mk_mask_function
    allocate(fmsk(nr)); fmsk = 1.0_num;

    do ir=1,nr
        if (r(ir) > r(nr)-r_msk) &
            fmsk(ir) = cos(pi*(r(ir)-r(nr)+r_msk)/2.0/r_msk)**(1.0/8.0_num)
    end do

    return
end subroutine mk_mask_function

! set data save parameters
subroutine set_data_output_info
    ! ig = (total, frequency, index) of saves
    select case (go_unit)
    case default
        ig(1) = nint(t_max/dtg) + 1
        ig(2) = nint(dtg/dt)
    case ('au')
        ig(1) = nint(t_max/dtg) + 1
        ig(2) = nint(dtg/dt)
    case ('cyc')
        ig(1) = floor(t_max/pump%T0/dtg) + 1
        ig(2) = floor(dtg*pump%T0/dt)
    end select
    ig(3) = 1

    return
end subroutine set_data_output_info

! store data for time step k
subroutine store_data(k)
    integer, intent(in) :: k

    complex(num), allocatable :: phif_dyn(:,:,:)
    real(num), allocatable :: tf_dyn(:), W_dyn(:,:)

    ! catch overflow from save indexing mismatch
    ! re-allocate _out arrays to avoid segfault
    if (ig(3) > ig(1)) then
        allocate(phif_dyn(l_max+1,nr,ig(3)), tf_dyn(ig(3)), W_dyn(ig(3),npe_bins))

        phif_dyn(:,:,1:ig(1)) = phif_out
        tf_dyn(1:ig(1)) = tf_out
        W_dyn(1:ig(1),:) = W

        call move_alloc(phif_dyn, phif_out)
        call move_alloc(tf_dyn, tf_out)
        call move_alloc(W_dyn, W)
        ig(1) = ig(1) + 1
    end if

    phif_out(:,:,ig(3)) = wavefn % phi
    tf_out(ig(3)) = t(k) + dt

    ! calculate photo-electron spectrum
    call calc_spectrum(wavefn, E, dr, winop_ord, W(ig(3),:))
    ig(3) = ig(3) + 1

    return
end subroutine store_data

! read input deck
subroutine read_deck
    character(len=128) :: label, temp_profile
    open(newunit=iu, file=deck_name, status='old', action='read')
    read(iu,*) ! global simulation time
    read(iu,*) label, t_max
    read(iu,*) label, dt
    read(iu,*)
    read(iu,*) ! wavefunction domain
    read(iu,*) label, r_max
    read(iu,*) label, l_max
    read(iu,*) label, dr
    read(iu,*) label, r_msk
    read(iu,*)
    read(iu,*) ! atomic properties
    read(iu,*) label, wavefn%Z
    read(iu,*) label, n_init, l_init, wavefn%m
    read(iu,*) label, tau
    read(iu,*) label, dtau
    read(iu,*)
    read(iu,*) ! p-electron spectrum
    read(iu,*) label, e_lim(1), e_lim(2)
    read(iu,*) label, de
    read(iu,*) label, winop_ord
    read(iu,*)
    read(iu,*) ! c-level occupations
    read(iu,'(A)') label
    nc = count([(label(j:j),j=1,len(label))] == '.')
    if (nc > 0) then
        allocate(Ec(nc),lc(nc))
        read(label(4:),*) Ec
    end if
    read(iu,'(A)') label
    if (nc > 0) read(label(4:),*) lc
    read(iu,*)
    read(iu,*) ! pump field parameters
    read(iu,*) label, temp_profile
    allocate( character(len=len(trim(adjustl(temp_profile)))) :: pump%profile )
    pump%profile = trim(adjustl(temp_profile))
    read(iu,*) label, pump % E0
    read(iu,*) label, pump % omg0
    if (pump%profile(1:5) == 'trapz') read(iu,*) label, pump % Ncyc_pl
    if (pump%profile(1:5) == 'sine2') read(iu,*) label, pump % Tp
    if (pump%profile(1:5) == 'gauss') read(iu,*) label, pump % Tfwhm
    read(iu,*) label, pump % CEP
    if (pump%profile(1:5) == 'trapz') read(iu,*) label, pump % t_on
    if (pump%profile(1:5) == 'sine2') read(iu,*) label, pump % t_on
    if (pump%profile(1:5) == 'gauss') read(iu,*) label, pump % Tpk
    read(iu,*)
    read(iu,*) ! probe field parameters
    read(iu,*) label, temp_profile
    allocate( character(len=len(trim(adjustl(temp_profile)))) :: probe%profile )
    probe%profile = trim(adjustl(temp_profile))
    read(iu,*) label, probe % E0
    read(iu,*) label, probe % omg0
    if (probe%profile(1:5) == 'trapz') read(iu,*) label, probe % Ncyc_pl
    if (probe%profile(1:5) == 'sine2') read(iu,*) label, probe % Tp
    if (probe%profile(1:5) == 'gauss') read(iu,*) label, probe % Tfwhm
    read(iu,*) label, probe % CEP
    if (probe%profile(1:5) == 'trapz') read(iu,*) label, probe % t_on
    if (probe%profile(1:5) == 'sine2') read(iu,*) label, probe % t_on
    if (probe%profile(1:5) == 'gauss') read(iu,*) label, probe % Tpk
    read(iu,*)
    read(iu,*) ! data output
    read(iu,*) label, dtg, go_unit
    read(iu,*)
    read(iu,*) ! update runtime log every N time steps
    read(iu,*) label, nt_update
    close(iu)
    return
end subroutine read_deck

! initialize runtime log
subroutine init_runtime_log
    open(newunit=iu, file='runtime.log', status='replace')
    write(iu,'(A,I2)') 'Rabi-split version: '//crv//'_x',8*num
    write(iu,'(A)') 'executed: '//exec_date//' '//exec_time
    write(iu,'(A,I2)') 'OMP threads: ', OMP_get_max_threads()
    write(iu,'(A)')
    write(iu,'(A)') 'Preparing atomic state...'
    write(iu,'(A,I1,A,I1,A,I1,A)') '(n,l,m) = (',n_init,',',l_init,',',wavefn%m,')'//NEW_LINE('A')
    close(iu)
    return
end subroutine init_runtime_log

! update runtime log
subroutine update_runtime_log
    open(newunit=iu, file='runtime.log', status='old', access='append')
    write(iu,'(A,F9.2,A,I7,A,F5.1,A)') 't = ',t(it)+dt,' au    nt: ',it,'  ',1e2*it/nt,'%'
    close(iu)
    return
end subroutine update_runtime_log

! finalize runtime log
subroutine end_runtime_log
    open(newunit=iu, file='runtime.log', status='replace')
    write(iu,'(A,I2)') 'Rabi-split version: '//crv//'_x',8*num
    write(iu,'(A)') 'executed: '//exec_date//' '//exec_time
    write(iu,'(A,I2)') 'OMP threads: ', OMP_get_max_threads()
    write(iu,'(A)') NEW_LINE('A')//'Simulation complete.'
    write(iu,'(A,F9.2,A)') 'total execution time:', (cpu_t(2)-cpu_t(1))/3.60d+3, ' core-hrs'
    write(iu,'(F9.2,A)') (wall_t(2)-wall_t(1))/3.60d+3, ' wall clock time (hrs)'
    close(iu)
    return
end subroutine end_runtime_log

end program main
