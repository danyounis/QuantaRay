! 2-electron Correlation Degree Calculator (correl2e.f08)
!
! Author: D. Younis
!         University of Rochester
!         Department of Physics
!
! Written: 2/22/2023
! Revised: 8/31/2023
!
! DESCRIPTION
!   Given a time-set of 2 one-dimensional electron wavefunctions,
!   calculates the single-particle density operator
!   $\rho(x,x')=\int\!dx_2\,\Psi(x,x_2)\Psi^*(x',x_2)$,
!   its eigenvalues, and the 2e degree of correlation as defined in Ref. [1].
!
! REFERENCES
!   [1] R. Grobe, K. Rzazewski, and J. H. Eberly, J. Phys. B 27, L503 (1994).
!       doi.org/10.1088/0953-4075/27/16/001
!
! INPUT
!   QuantaRay (qray) dataset containing full volumetric wavefunction data.
!
! OUTPUT
!   Raw binary files containing (1) the single-particle density operator (rho),
!   (2) its eigenvalues (eigs), and (3) the 2e correlation degree vs. time (Kt).
!
! USAGE
!   ./correl2e full/part qm_data.h5
!
! DEPENDS
!   LAPACK, BLAS, HDF5
!
! NOTES
!   - full: Compute everything.
!     part: Compute only the single-particle density operator.
!
!   - To read the output datasets in python, use, for example:
!       from scipy.io import FortranFile
!       with FortranFile('Re_rho.dat','r') as f:
!           Re_rho = f.read_reals('float64').reshape((nx,nx,nt), order='F')
!
!   - Depending on the wavefunction resolution (nx,ny), diagonalizing the
!     density operator can be an expensive task. One can instead call correl2e
!     with the "part" argument to bypass those subroutine calls and compute only
!     a few of the largest eigenvalues in python, using, for example:
!       from scipy.sparse.linalg import eigs
!       (w,v) = eigs(rho, k=?, which='LR')
!     where k is the no. of (w,v) = (eigenvalue, eigenvector) pairs to compute,
!     see the documentation for scipy.sparse.linalg.eigs.

program correl2e
    use prec
    use math
    use omp_lib

    implicit none

    character(len=4) :: mode
    character(len=24) :: dataset_name
    real(num) :: dx, dy
    integer :: nx, ny, nt

    real(num), allocatable :: x(:), y(:), Kt(:)
    complex(num), allocatable :: psi(:,:,:), rho(:,:,:), eigs(:,:)

    ! utilize all available threads
    call OMP_set_num_threads(OMP_get_max_threads())
    call OMP_set_nested(.false.)

    ! load qray output data
    call get_command_argument(1,mode)
    call get_command_argument(2,dataset_name)
    call load_data()
    ! calculate single-particle density operator
    call calc_sp_dens_op()
if (mode == 'full') then
    ! its eigenvalues
    call calc_sp_dens_op_eigenvals()
    ! and the 2e degree of correlation
    call calc_deg_correl()
end if
    ! output the calculated quantities
    call write_data()
    ! terminate correl2e
    call exit()

contains

! Load the qray dataset specified by 'dataset_name'
SUBROUTINE load_data()
    USE HDF5
    IMPLICIT NONE
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: Re_psi, Im_psi
    INTEGER(HID_T) :: H5T_REAL_PRECX, file_id, dspace_id, dset_id
    INTEGER(HSIZE_T) :: dims2(2), dims3(3)
    INTEGER :: error

    ! initialize HDF5 interface and open file
    CALL H5OPEN_F(error)
    CALL H5FOPEN_F(dataset_name, H5F_ACC_RDONLY_F, file_id, error)

    ! set H5 precision
    IF (num.EQ.4) H5T_REAL_PRECX = H5T_IEEE_F32LE
    IF (num.EQ.8) H5T_REAL_PRECX = H5T_IEEE_F64LE

    ! load grid data
    CALL H5DOPEN_F(file_id, 'x', dset_id, error)
    CALL H5DGET_SPACE_F(dset_id, dspace_id, error)
    CALL H5SGET_SIMPLE_EXTENT_DIMS_F(dspace_id, dims2, dims2, error)
    ALLOCATE(x(dims2(2))) ! dims2 = (/1,nx/)
    CALL H5DREAD_F(dset_id, H5T_REAL_PRECX, x, dims2, error)

    CALL H5DOPEN_F(file_id, 'y', dset_id, error)
    CALL H5DGET_SPACE_F(dset_id, dspace_id, error)
    CALL H5SGET_SIMPLE_EXTENT_DIMS_F(dspace_id, dims2, dims2, error)
    ALLOCATE(y(dims2(2))) ! dims2 = (/1,ny/)
    CALL H5DREAD_F(dset_id, H5T_REAL_PRECX, y, dims2, error)

    ! load wavefunction data
    CALL H5DOPEN_F(file_id, 'Re(psi)', dset_id, error)
    CALL H5DGET_SPACE_F(dset_id, dspace_id, error)
    CALL H5SGET_SIMPLE_EXTENT_DIMS_F(dspace_id, dims3, dims3, error)
    ALLOCATE(Re_psi(dims3(1),dims3(2),dims3(3))) ! dims3 = (/nx,ny,nt/)
    CALL H5DREAD_F(dset_id, H5T_REAL_PRECX, Re_psi, dims3, error)

    CALL H5DOPEN_F(file_id, 'Im(psi)', dset_id, error)
    CALL H5DGET_SPACE_F(dset_id, dspace_id, error)
    CALL H5SGET_SIMPLE_EXTENT_DIMS_F(dspace_id, dims3, dims3, error)
    ALLOCATE(Im_psi(dims3(1),dims3(2),dims3(3))) ! dims3 = (/nx,ny,nt/)
    CALL H5DREAD_F(dset_id, H5T_REAL_PRECX, Im_psi, dims3, error)

    ! close file and end HDF5 interface
    CALL H5FCLOSE_F(file_id, error)
    CALL H5CLOSE_F(error)

    ! store wavefunction and grid parameters
    ALLOCATE(psi(dims3(1),dims3(2),dims3(3)))
    psi = Re_psi + i*Im_psi
    nx = SIZE(psi,1); ny = SIZE(psi,2); nt = SIZE(psi,3);
    dx = x(2)-x(1); dy = y(2)-y(1);
    DEALLOCATE(Re_psi,Im_psi)

    RETURN
END SUBROUTINE load_data

! Calculate the single-particle density operator
! $\rho(x,x')=\int\!dx_2\,\Psi(x,x_2)\Psi^*(x',x_2)$
subroutine calc_sp_dens_op()
    implicit none
    integer :: it, ix1, ix2
    allocate(rho(nx,nx,nt), source=(0.0_num,0.0_num))
    !$OMP PARALLEL DO SHARED(rho) PRIVATE(ix1,ix2)
    do it=1,nt
    do ix2=1,nx
    do ix1=1,nx
        rho(ix1,ix2,it) = dx*dy*trapz(psi(ix1,:,it)*conjg(psi(ix2,:,it)))
    end do
    end do
    end do
    !$OMP END PARALLEL DO
    return
end subroutine calc_sp_dens_op

! Calculate the eigenvalues of the single-particle density operator
subroutine calc_sp_dens_op_eigenvals()
    implicit none
    complex(num), allocatable :: rho_temp(:,:,:)
    integer :: it
    ! temporary variables, LAPACK eigenvalue subroutine
    EXTERNAL :: CGEES, ZGEES, SELECT
    INTEGER :: INFO, SDIM, LWORK
    LOGICAL, ALLOCATABLE :: BWORK(:)
    ! for CGEES
    REAL, ALLOCATABLE :: RWORKc(:)
    COMPLEX, ALLOCATABLE :: Wc(:), VSc(:,:), WORKc(:)
    ! for ZGEES
    DOUBLE PRECISION, ALLOCATABLE :: RWORKz(:)
    COMPLEX*16, ALLOCATABLE :: Wz(:), VSz(:,:), WORKz(:)

    LWORK = 2*nx
    allocate(eigs(nx,nt), source=(0.0_num,0.0_num))
    allocate(rho_temp(nx,nx,nt), source=rho)

    select case (num)
    case (4) ! single precision

        allocate(BWORK(nx), RWORKc(nx), Wc(nx), VSc(1,1), WORKc(LWORK))
        !$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(eigs,rho_temp)
        do it=1,nt
            CALL CGEES('N','N',UNUSED,nx,rho_temp(:,:,it),nx,SDIM,Wc,&
                VSc,1,WORKc,LWORK,RWORKc,BWORK,INFO)
            eigs(:,it) = Wc
        end do
        !$OMP END PARALLEL DO

    case (8) ! double precision

        allocate(BWORK(nx), RWORKz(nx), Wz(nx), VSz(1,1), WORKz(LWORK))
        !$OMP PARALLEL DO DEFAULT(FIRSTPRIVATE) SHARED(eigs,rho_temp)
        do it=1,nt
            CALL ZGEES('N','N',UNUSED,nx,rho_temp(:,:,it),nx,SDIM,Wz,&
                VSz,1,WORKz,LWORK,RWORKz,BWORK,INFO)
            eigs(:,it) = Wz
        end do
        !$OMP END PARALLEL DO

    end select

    deallocate(rho_temp)

    return
end subroutine calc_sp_dens_op_eigenvals

! Calculate the 2-electron degree of correlation
subroutine calc_deg_correl()
    implicit none
    allocate(Kt(nt), source=0.0_num)
    Kt = 1.0_num/sum(abs(2.*real(eigs,num))**2, dim=1)
    return
end subroutine calc_deg_correl

! Output calculated quantities to unformatted binary files
subroutine write_data()
    implicit none
    integer :: iu
    open(newunit=iu, file='Re_rho.dat', form='unformatted')
    write(iu) real(rho,num)
    close(iu)
    open(newunit=iu, file='Im_rho.dat', form='unformatted')
    write(iu) aimag(rho)
    close(iu)
if (mode == 'full') then
    open(newunit=iu, file='Re_eigs.dat', form='unformatted')
    write(iu) real(eigs,num)
    close(iu)
    open(newunit=iu, file='Im_eigs.dat', form='unformatted')
    write(iu) aimag(eigs) ! Im(eigs) ~ 0.0, rho is Hermitian
    close(iu)
    open(newunit=iu, file='Kt.dat', form='unformatted')
    write(iu) Kt
    close(iu)
end if
    return
end subroutine write_data

! Dummy function for LAPACK's eigenvalue subroutine
LOGICAL FUNCTION UNUSED(X,Y)
    IMPLICIT NONE
    REAL(num), INTENT(IN) :: X,Y
    UNUSED = .TRUE.
    RETURN
END FUNCTION UNUSED

end program correl2e
