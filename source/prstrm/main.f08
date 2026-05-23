! PairStorm (main)
!
! Author:  D. Younis
!          University of Rochester
!          Department of Physics
!
! Written: 1/24/2025
! Revised: 3/20/2026
!
! USAGE
!   mpirun -np <number of processes> prstrm <input.deck>
!
! DEPENDS
!   mpifort/13.3.0, mpirun/4.1.6

program main
    use prec
    use qed_mcpm
    use math, only: pi, linspace, init_RNG
    use mpi_f08
    ! use omp_lib

    implicit none

    ! code revision date
    character(len=8), parameter :: crv = '20260320'
    ! execution date / time
    character(len=8) :: exec_date
    character(len=6) :: exec_time
    ! total computation time
    real :: wall_t(2)

    ! input deck variables
    integer :: iu, N_seed
    real(num) :: lam0_cm_, Ep_min_, dt_max_, dt_scl_
    real(num) :: t_max, E_, r_(3), k_(3), Ro, I_A
    logical :: photon_dynamics_
    character(len=:), allocatable :: dt_alg, tables_path
    character(len=24) :: deck_name

    ! MPI variables
    type(MPI_Status) :: MPI_Status_
    type(MPI_Datatype) :: MPI_R_NUM
    type(MPI_Datatype) :: MPI_C_NUM
    character(len=3) :: process_rank_str
    integer :: size_cluster, process_rank, ierr
    integer, parameter :: root = 0

    integer :: nn, N_ens ! no. ensemble members assigned to this MPI process

    type(plist_t) :: leptons
    type(plist_t) :: photons

    type(force) :: force_
    type(field) :: field_
    type(qmc_table_t) :: table_

    ! get date/time of execution
    call date_and_time(date=exec_date, time=exec_time)

    ! parse input deck
    call get_command_argument(1,deck_name)
    if (deck_name=='-v' .or. deck_name=='--version') then
        print '(A,I2)', 'PairStorm v.'//crv//'_x',8*num
        call exit()
    end if
    call init_RNG()
    call read_deck()

    ! initialize qed_mcpm module
    call set_qedmcpm_units(lam0_cm_)
    call set_qedmcpm_params(Ep_min_, dt_max_, dt_scl_, photon_dynamics_)

    call force_%set()
    call field_%set(Ex,Ey,Ez,Bx,By,Bz)
    call table_%load(tables_path)

    ! MPI initialization
    call MPI_Start_Interface()
    if (process_rank == root) wall_t(1) = MPI_Wtime()
    if (process_rank == root) call inquire_diag()

    ! distribute ensemble members to MPI processes
    N_ens = N_seed*(process_rank+1)/size_cluster - N_seed*process_rank/size_cluster ! implied integer division (floor)
    write(process_rank_str,'(I3)') process_rank

do nn=1,N_ens
    ! initialize seed lepton (id:=0)
    call leptons%create(new_particle(m=1., q=-1., id=0, &
        x=[0.0_num, 2*pi*r_(1), 2*pi*r_(2), 2*pi*r_(3)], &
        E=E_*mks%e/mks%m_e/mks%c**2, N=k_))

    !-----/ PAIRSTORM MAIN CYCLE /-----!
    do while (leptons%has_stat(1).or.photons%has_stat(1))
        !~ LEPTON PROPAGATION ~!
        ! if ((process_rank == root).and.(nn == 1)) print *, leptons%head%par%x
        call leptons%push(t_max=t_max, F=field_, G=force_)

        !~ PHOTON EMISSION ~!
        call leptons%emit(recv_lst=photons, F=field_, L=table_)

        !~ PHOTON PROPAGATION ~!
        call photons%push(t_max=t_max, F=field_, G=force_)

        !~ PAIR PRODUCTION ~!
        call photons%emit(recv_lst=leptons, F=field_, L=table_)
    end do

    call output_diag()
    call output_data()

    call leptons%reset()
    call photons%reset()
    call annihilated_photons%reset()
end do

    if (process_rank == root) then
        wall_t(2) = MPI_Wtime()
        print '(A,F7.3,A)', 'Execution time: ', (wall_t(2)-wall_t(1))/3.6e+3, ' hours'
    end if

    call MPI_Stop_Interface()
    call exit()

contains

real(num) function Ex(r,t)
    real(num), intent(in) :: r(3), t
    Ex = 0.0_num
end function Ex

real(num) function Ey(r,t)
    real(num), intent(in) :: r(3), t
    Ey = 0.0_num
end function Ey

real(num) function Ez(r,t)
    real(num), intent(in) :: r(3), t
    Ez = 0.0_num
end function Ez

real(num) function Bx(r,t)
    real(num), intent(in) :: r(3), t
    real(num) :: r_xy
    r_xy = 1.d-2*sqrt(r(1)**2 + r(2)**2)/k0_icm ! [m]
    ! r->[MKS], then B [MKS]->[CGS]->[sim.u.]

    ! infinitesimally-thin wire !
    ! Bx = -1.d+4*(mks%mu_0*I_A)/(2*pi) * sin(atan2(r(2),r(1))) / r_xy / F_

    ! finite-thickness wire (linear current density) !
    if (r_xy < Ro) then
        Bx = -1.d+4*(mks%mu_0*I_A)/(2*pi) * sin(atan2(r(2),r(1))) * (r_xy**2)/(Ro**3) / F_
    else
        Bx = -1.d+4*(mks%mu_0*I_A)/(2*pi) * sin(atan2(r(2),r(1))) / r_xy / F_
    end if
end function Bx

real(num) function By(r,t)
    real(num), intent(in) :: r(3), t
    real(num) :: r_xy
    r_xy = 1.d-2*sqrt(r(1)**2 + r(2)**2)/k0_icm ! [m]
    ! r->[MKS], then B [MKS]->[CGS]->[sim.u.]

    ! infinitesimally-thin wire !
    ! By = +1.d+4*(mks%mu_0*I_A)/(2*pi) * cos(atan2(r(2),r(1))) / r_xy / F_

    ! finite-thickness wire (linear current density) !
    if (r_xy < Ro) then
        By = +1.d+4*(mks%mu_0*I_A)/(2*pi) * cos(atan2(r(2),r(1))) * (r_xy**2)/(Ro**3) / F_
    else
        By = +1.d+4*(mks%mu_0*I_A)/(2*pi) * cos(atan2(r(2),r(1))) / r_xy / F_
    end if
end function By

real(num) function Bz(r,t)
    real(num), intent(in) :: r(3), t
    Bz = 0.0_num
end function Bz

! read input deck
subroutine read_deck()
    integer :: iu
    character(len=64) :: label
    character(len=128) :: label_long
    open(newunit=iu, file=deck_name, status='old', action='read')
    read(iu,*) label, lam0_cm_ ! normalization length [cm]
    read(iu,*)
    read(iu,*) label, N_seed ! number of seed electrons
    read(iu,*)
    read(iu,*) label, label ! time-stepping algorithm
    allocate(character(len=len(trim(adjustl(label)))) :: dt_alg)
    dt_alg = trim(adjustl(label))
    read(iu,*) label, dt_max_ ! maximum time-step [sim.u.]
    read(iu,*) label, dt_scl_ ! adaptive time-step scale factor, in (0,1]
    read(iu,*)
    read(iu,*) label, Ep_min_ ! cutoff energy for photon emission [sim.u.]
    read(iu,*) label, photon_dynamics_ ! toggle photon propagation and e-/e+ pair production
    read(iu,*)
    read(iu,*) label ! path to QED lookup tables
    read(iu,*) label_long
    allocate(character(len=len(trim(adjustl(label_long)))) :: tables_path)
    tables_path = trim(adjustl(label_long))
    read(iu,*)
    read(iu,*) ! Electron parameters
    read(iu,*) label, E_ ! energy [eV]
    read(iu,*) label, r_(1), r_(2), r_(3) ! injection point [sim.u.]
    read(iu,*) label, k_(1), k_(2), k_(3) ! propagation direction (un-normalized)
    read(iu,*)
    read(iu,*) ! Field parameters
    read(iu,*) label, Ro ! wire thickness [m]
    read(iu,*) label, I_A ! current [A] (+Z oriented)
    read(iu,*)
    read(iu,*) ! End criterion
    read(iu,*) label, t_max
    close(iu)
    ! Error-checking !
    if (dt_alg == 'fixed') &
        call ERROR('Fixed time-stepping has not been implemented yet.', &
            fix='Change `dt_alg` from `fixed` to `adaptive`.', in=trim(adjustl(deck_name)))
    if (dt_alg/='fixed' .and. dt_alg/='adaptive') &
        call ERROR('Unknown time-stepping method `'//dt_alg//'`.', &
            fix='Choose between `fixed` and `adaptive`.', in=trim(adjustl(deck_name)))
    if (Ep_min_ <= 0.0_num) &
        call ERROR('Cutoff energy for photon emission (Ep_min) must be positive.', in=trim(adjustl(deck_name)))
    if (dt_scl_<=0.0_num .or. dt_scl_>1.0_num) &
        call ERROR('Adaptive time-step scale factor must be 0 < dt_scl <= 1.', in=trim(adjustl(deck_name)))
    return
end subroutine read_deck

! check if `diagnostics.txt` file already exists
! parse contents to preserve old min/max values
subroutine inquire_diag()
    integer :: iu, j
    logical :: b_reload
    character(len=32) :: label
    inquire(file='diagnostics.txt', exist=b_reload)
    if (b_reload) then
        open(newunit=iu, file='diagnostics.txt', status='old', action='read')
        do j=1,6; read(iu,*); end do; ! skip 6 lines
        read(iu,*) label, diag%min_qp(1), diag%min_qp(2)
        read(iu,*) label, diag%max_qp(1), diag%max_qp(2)
        do j=1,4; read(iu,*); end do; ! skip 4 lines
        read(iu,*) label, diag%min_prob(1), diag%min_prob(2)
        read(iu,*) label, diag%max_prob(1), diag%max_prob(2)
        close(iu)
    end if
    return
end subroutine inquire_diag

! output diagnostic data
subroutine output_diag()
    real(num), allocatable, dimension(:,:), save :: min_qp_, max_qp_
    real(num), allocatable, dimension(:,:), save :: min_prob_, max_prob_
    integer :: irank, l

    ! Gather min/max quantum parameter and
    ! emission probability values to Rank 0
    if ((process_rank == root).and.(size_cluster > 1)) then
        if (.not.allocated(min_qp_)) then
            allocate(min_qp_(2,0:size_cluster-1))
            allocate(max_qp_(2,0:size_cluster-1))
            allocate(min_prob_(2,0:size_cluster-1))
            allocate(max_prob_(2,0:size_cluster-1))
        end if

        ! insert Rank 0 data
        min_qp_(:,0) = diag%min_qp
        max_qp_(:,0) = diag%max_qp
        min_prob_(:,0) = diag%min_prob
        max_prob_(:,0) = diag%max_prob

        ! loop over other ranks
        do irank=1,size_cluster-1
            call MPI_Recv(min_qp_(:,irank),   2, MPI_R_NUM, irank, 1, MPI_COMM_WORLD, MPI_Status_, ierr)
            call MPI_Recv(max_qp_(:,irank),   2, MPI_R_NUM, irank, 2, MPI_COMM_WORLD, MPI_Status_, ierr)
            call MPI_Recv(min_prob_(:,irank), 2, MPI_R_NUM, irank, 3, MPI_COMM_WORLD, MPI_Status_, ierr)
            call MPI_Recv(max_prob_(:,irank), 2, MPI_R_NUM, irank, 4, MPI_COMM_WORLD, MPI_Status_, ierr)
        end do

        ! Overwrite Rank-0 data with global min/max values
        ForAll(l=1:2) diag%min_qp(l) = minval(min_qp_(l,:))
        ForAll(l=1:2) diag%max_qp(l) = maxval(max_qp_(l,:))
        ForAll(l=1:2) diag%min_prob(l) = minval(min_prob_(l,:))
        ForAll(l=1:2) diag%max_prob(l) = maxval(max_prob_(l,:))

    else if ((process_rank /= root).and.(size_cluster > 1)) then
        call MPI_Send(diag%min_qp,   2, MPI_R_NUM, root, 1, MPI_COMM_WORLD, ierr)
        call MPI_Send(diag%max_qp,   2, MPI_R_NUM, root, 2, MPI_COMM_WORLD, ierr)
        call MPI_Send(diag%min_prob, 2, MPI_R_NUM, root, 3, MPI_COMM_WORLD, ierr)
        call MPI_Send(diag%max_prob, 2, MPI_R_NUM, root, 4, MPI_COMM_WORLD, ierr)
    end if

    if (.not.photon_dynamics) then
        diag%min_qp(2) = 0.0_num
        diag%max_qp(2) = 0.0_num
        diag%min_prob(2) = 0.0_num
        diag%max_prob(2) = 0.0_num
    end if

    if (process_rank == root) then
        open(newunit=iu, file='diagnostics.txt', status='replace')
        write(iu,'(A,I2)') 'PairStorm version '//crv//'_x',8*num
        write(iu,'(A)') 'executed: '//exec_date//' '//exec_time//NEW_LINE('A')
        write(iu,'(A)') 'Quantum parameters (*1)'
        write(iu,'(A)') '------------------------------'
        write(iu,'(A)') '    Leptons     Photons'
        write(iu,'(A,ES9.3,A,ES9.3)') 'Min ', diag%min_qp(1), '   ', diag%min_qp(2)
        write(iu,'(A,ES9.3,A,ES9.3,A)') 'Max ', diag%max_qp(1), '   ', diag%max_qp(2), NEW_LINE('A')
        write(iu,'(A)') 'Emission probabilities (*1)'
        write(iu,'(A)') '------------------------------'
        write(iu,'(A)') '    e + nγ →γ   γ + nγ→e±'
        write(iu,'(A,ES9.3,A,ES9.3)') 'Min ', diag%min_prob(1), '   ', diag%min_prob(2)
        write(iu,'(A,ES9.3,A,ES9.3,A)') 'Max ', diag%max_prob(1), '   ', diag%max_prob(2), NEW_LINE('A')
        write(iu,'(A)') '(*1) Values encountered at some time-step.'
        close(iu)
    end if

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    return
end subroutine output_diag

! output linked-list data
subroutine output_data()
    ! Rank-0 variables
    real, allocatable :: charge_lepton(:)
    integer, allocatable :: id_lepton(:), id_photon(:), id_annihilated_photon(:)
    real(num), allocatable :: energy_lepton(:), energy_photon(:), energy_annihilated_photon(:)
    real(num), allocatable :: qparam_lepton(:), qparam_photon(:), qparam_annihilated_photon(:)
    real(num), allocatable :: x0_lepton(:), x0_photon(:), x0_annihilated_photon(:)
    real(num), allocatable :: xf_lepton(:,:), pf_lepton(:,:)
    real(num), allocatable :: xf_annihilated_photon(:,:), pf_annihilated_photon(:,:)
    real(num), allocatable :: p0_lepton(:)
    ! temporary variables
    integer, allocatable :: id(:)
    real, allocatable :: charge(:)
    real(num), allocatable :: energy(:), qparam(:), x0(:), p0(:), xf(:,:), pf(:,:)
    integer :: j, irank, isize, old_size, new_size

    ! Initialize data on Rank 0
    if (process_rank == root) then
        allocate(id_lepton(leptons%N), source=leptons%get_id())
        allocate(energy_lepton(leptons%N), source=leptons%get_energy())
        allocate(qparam_lepton(leptons%N), source=leptons%get_qparam())
        allocate(x0_lepton(leptons%N), source=leptons%get_ivalue('t'))
        allocate(p0_lepton(leptons%N), source=leptons%get_ivalue('E'))
        allocate(charge_lepton(leptons%N), source=leptons%get_charge())

        allocate(id_photon(photons%N), source=photons%get_id())
        allocate(energy_photon(photons%N), source=photons%get_energy())
        allocate(qparam_photon(photons%N), source=photons%get_qparam())
        allocate(x0_photon(photons%N), source=photons%get_ivalue('t'))

        allocate(id_annihilated_photon(annihilated_photons%N), source=annihilated_photons%get_id())
        allocate(energy_annihilated_photon(annihilated_photons%N), source=annihilated_photons%get_energy())
        allocate(qparam_annihilated_photon(annihilated_photons%N), source=annihilated_photons%get_qparam())
        allocate(x0_annihilated_photon(annihilated_photons%N), source=annihilated_photons%get_ivalue('t'))

        allocate(xf_lepton(leptons%N,3))
        allocate(pf_lepton(leptons%N,3))
        do j=1,3
            xf_lepton(:,j) = leptons%get_fvalue('x',j)
            pf_lepton(:,j) = leptons%get_fvalue('p',j)
        end do

        allocate(xf_annihilated_photon(annihilated_photons%N,3))
        allocate(pf_annihilated_photon(annihilated_photons%N,3))
        do j=1,3
            xf_annihilated_photon(:,j) = annihilated_photons%get_fvalue('x',j)
            pf_annihilated_photon(:,j) = annihilated_photons%get_fvalue('p',j)
        end do
    end if

    ! Gather all data to Rank 0
    if ((process_rank == root).and.(size_cluster > 1)) then

        do irank=1,size_cluster-1
            !--- lepton data ---!
            call MPI_Recv(isize, 1, MPI_INTEGER, irank, 0, MPI_COMM_WORLD, MPI_Status_, ierr)

            old_size = size(id_lepton)
            new_size = size(id_lepton) + isize

            allocate(id(new_size)); id(1:old_size) = id_lepton;
            allocate(energy(new_size)); energy(1:old_size) = energy_lepton;
            allocate(qparam(new_size)); qparam(1:old_size) = qparam_lepton;
            allocate(x0(new_size)); x0(1:old_size) = x0_lepton;
            allocate(p0(new_size)); p0(1:old_size) = p0_lepton;
            allocate(charge(new_size)); charge(1:old_size) = charge_lepton;
            allocate(xf(new_size,3)); xf(1:old_size,:) = xf_lepton;
            allocate(pf(new_size,3)); pf(1:old_size,:) = pf_lepton;

            ! id
            call MPI_Recv(id(old_size+1:new_size), isize, MPI_INTEGER, irank, 1, MPI_COMM_WORLD, MPI_Status_, ierr)
            call move_alloc(id, id_lepton)
            ! energy
            call MPI_Recv(energy(old_size+1:new_size), isize, MPI_R_NUM, irank, 2, MPI_COMM_WORLD, MPI_Status_, ierr)
            call move_alloc(energy, energy_lepton)
            ! quantum parameter
            call MPI_Recv(qparam(old_size+1:new_size), isize, MPI_R_NUM, irank, 3, MPI_COMM_WORLD, MPI_Status_, ierr)
            call move_alloc(qparam, qparam_lepton)
            ! creation time
            call MPI_Recv(x0(old_size+1:new_size), isize, MPI_R_NUM, irank, 4, MPI_COMM_WORLD, MPI_Status_, ierr)
            call move_alloc(x0, x0_lepton)
            ! creation energy
            call MPI_Recv(p0(old_size+1:new_size), isize, MPI_R_NUM, irank, 5, MPI_COMM_WORLD, MPI_Status_, ierr)
            call move_alloc(p0, p0_lepton)
            ! charge
            call MPI_Recv(charge(old_size+1:new_size), isize, MPI_REAL, irank, 6, MPI_COMM_WORLD, MPI_Status_, ierr)
            call move_alloc(charge, charge_lepton)
            ! final positions
            call MPI_Recv(xf(old_size+1:new_size,1), isize, MPI_R_NUM, irank, 7, MPI_COMM_WORLD, MPI_Status_, ierr)
            call MPI_Recv(xf(old_size+1:new_size,2), isize, MPI_R_NUM, irank, 8, MPI_COMM_WORLD, MPI_Status_, ierr)
            call MPI_Recv(xf(old_size+1:new_size,3), isize, MPI_R_NUM, irank, 9, MPI_COMM_WORLD, MPI_Status_, ierr)
            call move_alloc(xf, xf_lepton)
            ! final momenta
            call MPI_Recv(pf(old_size+1:new_size,1), isize, MPI_R_NUM, irank, 10, MPI_COMM_WORLD, MPI_Status_, ierr)
            call MPI_Recv(pf(old_size+1:new_size,2), isize, MPI_R_NUM, irank, 11, MPI_COMM_WORLD, MPI_Status_, ierr)
            call MPI_Recv(pf(old_size+1:new_size,3), isize, MPI_R_NUM, irank, 12, MPI_COMM_WORLD, MPI_Status_, ierr)
            call move_alloc(pf, pf_lepton)

            !--- photon data ---!
            call MPI_Recv(isize, 1, MPI_INTEGER, irank, 13, MPI_COMM_WORLD, MPI_Status_, ierr)

            old_size = size(id_photon)
            new_size = size(id_photon) + isize

            allocate(id(new_size)); id(1:old_size) = id_photon;
            allocate(energy(new_size)); energy(1:old_size) = energy_photon;
            allocate(qparam(new_size)); qparam(1:old_size) = qparam_photon;
            allocate(x0(new_size)); x0(1:old_size) = x0_photon;

            ! id
            call MPI_Recv(id(old_size+1:new_size), isize, MPI_INTEGER, irank, 14, MPI_COMM_WORLD, MPI_Status_, ierr)
            call move_alloc(id, id_photon)
            ! energy
            call MPI_Recv(energy(old_size+1:new_size), isize, MPI_R_NUM, irank, 15, MPI_COMM_WORLD, MPI_Status_, ierr)
            call move_alloc(energy, energy_photon)
            ! quantum parameter
            call MPI_Recv(qparam(old_size+1:new_size), isize, MPI_R_NUM, irank, 16, MPI_COMM_WORLD, MPI_Status_, ierr)
            call move_alloc(qparam, qparam_photon)
            ! creation time
            call MPI_Recv(x0(old_size+1:new_size), isize, MPI_R_NUM, irank, 17, MPI_COMM_WORLD, MPI_Status_, ierr)
            call move_alloc(x0, x0_photon)

            !--- annihilated-photon data ---!
            call MPI_Recv(isize, 1, MPI_INTEGER, irank, 18, MPI_COMM_WORLD, MPI_Status_, ierr)

            old_size = size(id_annihilated_photon)
            new_size = size(id_annihilated_photon) + isize

            allocate(id(new_size)); id(1:old_size) = id_annihilated_photon;
            allocate(energy(new_size)); energy(1:old_size) = energy_annihilated_photon;
            allocate(qparam(new_size)); qparam(1:old_size) = qparam_annihilated_photon;
            allocate(x0(new_size)); x0(1:old_size) = x0_annihilated_photon;
            allocate(xf(new_size,3)); xf(1:old_size,:) = xf_annihilated_photon;
            allocate(pf(new_size,3)); pf(1:old_size,:) = pf_annihilated_photon;

            ! id
            call MPI_Recv(id(old_size+1:new_size), isize, MPI_INTEGER, irank, 19, MPI_COMM_WORLD, MPI_Status_, ierr)
            call move_alloc(id, id_annihilated_photon)
            ! energy
            call MPI_Recv(energy(old_size+1:new_size), isize, MPI_R_NUM, irank, 20, MPI_COMM_WORLD, MPI_Status_, ierr)
            call move_alloc(energy, energy_annihilated_photon)
            ! quantum parameter
            call MPI_Recv(qparam(old_size+1:new_size), isize, MPI_R_NUM, irank, 21, MPI_COMM_WORLD, MPI_Status_, ierr)
            call move_alloc(qparam, qparam_annihilated_photon)
            ! creation time
            call MPI_Recv(x0(old_size+1:new_size), isize, MPI_R_NUM, irank, 22, MPI_COMM_WORLD, MPI_Status_, ierr)
            call move_alloc(x0, x0_annihilated_photon)
            ! final positions
            call MPI_Recv(xf(old_size+1:new_size,1), isize, MPI_R_NUM, irank, 23, MPI_COMM_WORLD, MPI_Status_, ierr)
            call MPI_Recv(xf(old_size+1:new_size,2), isize, MPI_R_NUM, irank, 24, MPI_COMM_WORLD, MPI_Status_, ierr)
            call MPI_Recv(xf(old_size+1:new_size,3), isize, MPI_R_NUM, irank, 25, MPI_COMM_WORLD, MPI_Status_, ierr)
            call move_alloc(xf, xf_annihilated_photon)
            ! final momenta
            call MPI_Recv(pf(old_size+1:new_size,1), isize, MPI_R_NUM, irank, 26, MPI_COMM_WORLD, MPI_Status_, ierr)
            call MPI_Recv(pf(old_size+1:new_size,2), isize, MPI_R_NUM, irank, 27, MPI_COMM_WORLD, MPI_Status_, ierr)
            call MPI_Recv(pf(old_size+1:new_size,3), isize, MPI_R_NUM, irank, 28, MPI_COMM_WORLD, MPI_Status_, ierr)
            call move_alloc(pf, pf_annihilated_photon)
        end do

    else if ((process_rank /= root).and.(size_cluster > 1)) then

        !--- lepton data ---!
        call MPI_Send(leptons%N, 1, MPI_INTEGER, root, 0, MPI_COMM_WORLD, ierr)
        call MPI_Send(leptons%get_id(), leptons%N, MPI_INTEGER, root, 1, MPI_COMM_WORLD, ierr)
        call MPI_Send(leptons%get_energy(), leptons%N, MPI_R_NUM, root, 2, MPI_COMM_WORLD, ierr)
        call MPI_Send(leptons%get_qparam(), leptons%N, MPI_R_NUM, root, 3, MPI_COMM_WORLD, ierr)
        call MPI_Send(leptons%get_ivalue('t'), leptons%N, MPI_R_NUM, root, 4, MPI_COMM_WORLD, ierr)
        call MPI_Send(leptons%get_ivalue('E'), leptons%N, MPI_R_NUM, root, 5, MPI_COMM_WORLD, ierr)
        call MPI_Send(leptons%get_charge(), leptons%N, MPI_REAL, root, 6, MPI_COMM_WORLD, ierr)
        call MPI_Send(leptons%get_fvalue('x',1), leptons%N, MPI_R_NUM, root, 7, MPI_COMM_WORLD, ierr)
        call MPI_Send(leptons%get_fvalue('x',2), leptons%N, MPI_R_NUM, root, 8, MPI_COMM_WORLD, ierr)
        call MPI_Send(leptons%get_fvalue('x',3), leptons%N, MPI_R_NUM, root, 9, MPI_COMM_WORLD, ierr)
        call MPI_Send(leptons%get_fvalue('p',1), leptons%N, MPI_R_NUM, root, 10, MPI_COMM_WORLD, ierr)
        call MPI_Send(leptons%get_fvalue('p',2), leptons%N, MPI_R_NUM, root, 11, MPI_COMM_WORLD, ierr)
        call MPI_Send(leptons%get_fvalue('p',3), leptons%N, MPI_R_NUM, root, 12, MPI_COMM_WORLD, ierr)

        !--- photon data ---!
        call MPI_Send(photons%N, 1, MPI_INTEGER, root, 13, MPI_COMM_WORLD, ierr)
        call MPI_Send(photons%get_id(), photons%N, MPI_INTEGER, root, 14, MPI_COMM_WORLD, ierr)
        call MPI_Send(photons%get_energy(), photons%N, MPI_R_NUM, root, 15, MPI_COMM_WORLD, ierr)
        call MPI_Send(photons%get_qparam(), photons%N, MPI_R_NUM, root, 16, MPI_COMM_WORLD, ierr)
        call MPI_Send(photons%get_ivalue('t'), photons%N, MPI_R_NUM, root, 17, MPI_COMM_WORLD, ierr)

        !--- annihilated-photon data ---!
        call MPI_Send(annihilated_photons%N, 1, MPI_INTEGER, root, 18, MPI_COMM_WORLD, ierr)
        call MPI_Send(annihilated_photons%get_id(), annihilated_photons%N, MPI_INTEGER, root, 19, MPI_COMM_WORLD, ierr)
        call MPI_Send(annihilated_photons%get_energy(), annihilated_photons%N, MPI_R_NUM, root, 20, MPI_COMM_WORLD, ierr)
        call MPI_Send(annihilated_photons%get_qparam(), annihilated_photons%N, MPI_R_NUM, root, 21, MPI_COMM_WORLD, ierr)
        call MPI_Send(annihilated_photons%get_ivalue('t'), annihilated_photons%N, MPI_R_NUM, root, 22, MPI_COMM_WORLD, ierr)
        call MPI_Send(annihilated_photons%get_fvalue('x',1), annihilated_photons%N, MPI_R_NUM, root, 23, MPI_COMM_WORLD, ierr)
        call MPI_Send(annihilated_photons%get_fvalue('x',2), annihilated_photons%N, MPI_R_NUM, root, 24, MPI_COMM_WORLD, ierr)
        call MPI_Send(annihilated_photons%get_fvalue('x',3), annihilated_photons%N, MPI_R_NUM, root, 25, MPI_COMM_WORLD, ierr)
        call MPI_Send(annihilated_photons%get_fvalue('p',1), annihilated_photons%N, MPI_R_NUM, root, 26, MPI_COMM_WORLD, ierr)
        call MPI_Send(annihilated_photons%get_fvalue('p',2), annihilated_photons%N, MPI_R_NUM, root, 27, MPI_COMM_WORLD, ierr)
        call MPI_Send(annihilated_photons%get_fvalue('p',3), annihilated_photons%N, MPI_R_NUM, root, 28, MPI_COMM_WORLD, ierr)

    end if

    ! Write / overwrite data to file
    if (process_rank == root) then
        open(newunit=iu, file='id_lepton.dat', action='write', form='unformatted', &
            access='stream', position='append', status='unknown')
        write(iu) id_lepton
        close(iu)
        open(newunit=iu, file='id_photon.dat', action='write', form='unformatted', &
            access='stream', position='append', status='unknown')
        write(iu) id_photon
        close(iu)
        open(newunit=iu, file='id_annihilated-photon.dat', action='write', form='unformatted', &
            access='stream', position='append', status='unknown')
        write(iu) id_annihilated_photon
        close(iu)
        open(newunit=iu, file='energy_lepton.dat', action='write', form='unformatted', &
            access='stream', position='append', status='unknown')
        write(iu) energy_lepton
        close(iu)
        open(newunit=iu, file='energy_photon.dat', action='write', form='unformatted', &
            access='stream', position='append', status='unknown')
        write(iu) energy_photon
        close(iu)
        open(newunit=iu, file='energy_annihilated-photon.dat', action='write', form='unformatted', &
            access='stream', position='append', status='unknown')
        write(iu) energy_annihilated_photon
        close(iu)
        open(newunit=iu, file='qparam_lepton.dat', action='write', form='unformatted', &
            access='stream', position='append', status='unknown')
        write(iu) qparam_lepton
        close(iu)
        open(newunit=iu, file='qparam_photon.dat', action='write', form='unformatted', &
            access='stream', position='append', status='unknown')
        write(iu) qparam_photon
        close(iu)
        open(newunit=iu, file='qparam_annihilated-photon.dat', action='write', form='unformatted', &
            access='stream', position='append', status='unknown')
        write(iu) qparam_annihilated_photon
        close(iu)
        open(newunit=iu, file='x0_lepton.dat', action='write', form='unformatted', &
            access='stream', position='append', status='unknown')
        write(iu) x0_lepton
        close(iu)
        open(newunit=iu, file='x0_photon.dat', action='write', form='unformatted', &
            access='stream', position='append', status='unknown')
        write(iu) x0_photon
        close(iu)
        open(newunit=iu, file='x0_annihilated-photon.dat', action='write', form='unformatted', &
            access='stream', position='append', status='unknown')
        write(iu) x0_annihilated_photon
        close(iu)
        open(newunit=iu, file='p0_lepton.dat', action='write', form='unformatted', &
            access='stream', position='append', status='unknown')
        write(iu) p0_lepton
        close(iu)
        open(newunit=iu, file='charge_lepton.dat', action='write', form='unformatted', &
            access='stream', position='append', status='unknown')
        write(iu) charge_lepton
        close(iu)
        open(newunit=iu, file='x1f_lepton.dat', action='write', form='unformatted', &
            access='stream', position='append', status='unknown')
        write(iu) xf_lepton(:,1)
        close(iu)
        open(newunit=iu, file='x2f_lepton.dat', action='write', form='unformatted', &
            access='stream', position='append', status='unknown')
        write(iu) xf_lepton(:,2)
        close(iu)
        open(newunit=iu, file='x3f_lepton.dat', action='write', form='unformatted', &
            access='stream', position='append', status='unknown')
        write(iu) xf_lepton(:,3)
        close(iu)
        open(newunit=iu, file='p1f_lepton.dat', action='write', form='unformatted', &
            access='stream', position='append', status='unknown')
        write(iu) pf_lepton(:,1)
        close(iu)
        open(newunit=iu, file='p2f_lepton.dat', action='write', form='unformatted', &
            access='stream', position='append', status='unknown')
        write(iu) pf_lepton(:,2)
        close(iu)
        open(newunit=iu, file='p3f_lepton.dat', action='write', form='unformatted', &
            access='stream', position='append', status='unknown')
        write(iu) pf_lepton(:,3)
        close(iu)
        open(newunit=iu, file='x1f_annihilated-photon.dat', action='write', form='unformatted', &
            access='stream', position='append', status='unknown')
        write(iu) xf_annihilated_photon(:,1)
        close(iu)
        open(newunit=iu, file='x2f_annihilated-photon.dat', action='write', form='unformatted', &
            access='stream', position='append', status='unknown')
        write(iu) xf_annihilated_photon(:,2)
        close(iu)
        open(newunit=iu, file='x3f_annihilated-photon.dat', action='write', form='unformatted', &
            access='stream', position='append', status='unknown')
        write(iu) xf_annihilated_photon(:,3)
        close(iu)
        open(newunit=iu, file='p1f_annihilated-photon.dat', action='write', form='unformatted', &
            access='stream', position='append', status='unknown')
        write(iu) pf_annihilated_photon(:,1)
        close(iu)
        open(newunit=iu, file='p2f_annihilated-photon.dat', action='write', form='unformatted', &
            access='stream', position='append', status='unknown')
        write(iu) pf_annihilated_photon(:,2)
        close(iu)
        open(newunit=iu, file='p3f_annihilated-photon.dat', action='write', form='unformatted', &
            access='stream', position='append', status='unknown')
        write(iu) pf_annihilated_photon(:,3)
        close(iu)
    end if

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    return
end subroutine output_data

subroutine MPI_Start_Interface()
    implicit none
    call MPI_Init(ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, size_cluster, ierr)
    call MPI_Comm_rank(MPI_COMM_WORLD, process_rank, ierr)
    ! set MPI floating-point precision
    select case (num)
    case (4)
        MPI_R_NUM = MPI_REAL
        MPI_C_NUM = MPI_COMPLEX
    case (8)
        MPI_R_NUM = MPI_DOUBLE_PRECISION
        MPI_C_NUM = MPI_DOUBLE_COMPLEX
    end select
    return
end subroutine MPI_Start_Interface

subroutine MPI_Stop_Interface()
    implicit none
    call MPI_Finalize(ierr)
    return
end subroutine MPI_Stop_Interface

end program main
