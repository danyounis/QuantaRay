MODULE DATAIO

USE prec
USE hdf5
USE quantum
USE emfm
USE vdm

IMPLICIT NONE

CONTAINS

SUBROUTINE WRITE_DATA_QCR_2D(t, x, y, VD, field, wavefn, tf, psif, tp, psip, Vp, S, fmsk, run_ccr)
    IMPLICIT NONE
    REAL(num),    INTENT(IN) :: t(:), x(:), y(:), tf(:), tp(:), Vp(:,:,:), fmsk(:,:)
    COMPLEX(num), INTENT(IN) :: psif(:,:,:), psip(:,:,:)
    LOGICAL,      INTENT(IN) :: run_ccr
    TYPE(vdet),   INTENT(IN) :: VD(:)
    TYPE(emf),    INTENT(IN) :: field
    TYPE(SchrodingerWavefunction2D), INTENT(IN) :: wavefn
    TYPE(tSURFF2D), INTENT(IN) :: S

    INTEGER :: nt, nx, ny, ntf, ntp, N_VD

    REAL(num), DIMENSION(:,:),   ALLOCATABLE :: VD_loc, VD_rho, VD_phase
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: VD_Krt, VD_Jrt
    REAL(num), DIMENSION(:,:),   ALLOCATABLE :: Ert, Art, Crt

    INTEGER, PARAMETER :: nvars1 = 6, nvars2 = 20
    CHARACTER(LEN=10)  :: varnames1(nvars1), varnames2(nvars2)

    INTEGER :: n, i, j, error
    INTEGER, PARAMETER :: rank2 = 2, rank3 = 3
    INTEGER(HSIZE_T)   :: dims1(4,2), dims2(13,2), dims_vd(3), dims_psif(3), dims_psip(3)

    CHARACTER(LEN=10), PARAMETER :: filename1 = 'vd_data.h5'
    CHARACTER(LEN=10), PARAMETER :: filename2 = 'qm_data.h5'

    INTEGER(HID_T) :: file_id, dspace_id, dset_id
    INTEGER(HID_T) :: H5T_REAL_PRECX

    nt = SIZE(t)
    nx = SIZE(x)
    ny = SIZE(y)
    ntf = SIZE(tf)
    ntp = SIZE(tp)
    N_VD = SIZE(VD)

    varnames1(1) = 't';        dims1(1,:) = (/1,nt/)
    varnames1(2) = 'VD_loc';   dims1(2,:) = (/N_VD,2/)
    varnames1(3) = 'VD_rho';   dims1(3,:) = (/N_VD,nt/)
    varnames1(4) = 'VD_phase'; dims1(4,:) = (/N_VD,nt/)
    varnames1(5) = 'VD_Krt';   dims_vd    = (/2,nt,N_VD/)
    varnames1(6) = 'VD_Jrt'; ! dims_vd

    varnames2(1)  = 't';            dims2(1,:)  = (/1,nt/)
    varnames2(2)  = 'tf';           dims2(2,:)  = (/1,ntf/)
    varnames2(3)  = 'tp';           dims2(3,:)  = (/1,ntp/)
    varnames2(4)  = 'x';            dims2(4,:)  = (/1,nx/)
    varnames2(5)  = 'y';            dims2(5,:)  = (/1,ny/)
    varnames2(6)  = 'mask';         dims2(6,:)  = (/nx,ny/)
    varnames2(7)  = 'Ert';          dims2(7,:)  = (/2,nt/)
    varnames2(8)  = 'Art';          dims2(8,:)  = (/2,nt/)
    varnames2(9)  = 'Crt';          dims2(9,:)  = (/2,nt/)
    varnames2(10) = 'kx';           dims2(10,:) = (/1,S%nk(1)/)
    varnames2(11) = 'ky';           dims2(11,:) = (/1,S%nk(2)/)
    varnames2(12) = 'Re(tSURFF)';   dims2(12,:) = S%nk
    varnames2(13) = 'Im(tSURFF)';   dims2(13,:) = S%nk
    varnames2(14) = 'Re(psi)';      dims_psif   = (/nx,ny,ntf/)
    varnames2(15) = 'Im(psi)';    ! dims_psif
    varnames2(16) = 'Vp';         ! dims_psif
    varnames2(17) = 'Re(p-psi)';    dims_psip  = (/ntp,nx,2/)
    varnames2(18) = 'Im(p-psi)';  ! dims_psip
    varnames2(19) = 'norm(psi)';  ! dims2(1,:)
    varnames2(20) = 'enrg(psi)';  ! dims2(1,:)

    IF (run_ccr) THEN
    ALLOCATE(VD_loc(N_VD,2), VD_rho(N_VD,nt), VD_phase(N_VD,nt), VD_Krt(2,nt,N_VD), VD_Jrt(2,nt,N_VD))

    ! get virtual detector data
    DO n=1,N_VD
        ! locations
        VD_loc(n,:) = (/VD(n)%xl, VD(n)%yl/)
        ! probability density
        VD_rho(n,:) = VD(n)%rho
        ! phase
        VD_phase(n,:) = VD(n)%phase
        ! momentum and current density
        VD_Krt(:,:,n) = VD(n)%Krt
        VD_Jrt(:,:,n) = VD(n)%Jrt
    END DO
    END IF

    ALLOCATE(Ert(2,nt), Art(2,nt), Crt(2,nt))

    ! get electric field data
    Ert(1,:) = field%Ex
    Ert(2,:) = field%Ey
    ! get vector potential data
    Art(1,:) = field%Ax
    Art(2,:) = field%Ay
    ! get excursion data
    Crt(1,:) = field%Cx
    Crt(2,:) = field%Cy

    ! initialize HDF5 interface
    CALL H5OPEN_F(error)

    ! set H5 precision
    IF (num.EQ.4) H5T_REAL_PRECX = H5T_IEEE_F32LE
    IF (num.EQ.8) H5T_REAL_PRECX = H5T_IEEE_F64LE

    ! OUTPUT VIRTUAL DETECTOR DATA !
    IF (run_ccr) THEN
    CALL H5FCREATE_F(filename1, H5F_ACC_EXCL_F, file_id, error)
    DO i=1,nvars1
        SELECT CASE (i)
        CASE (1:4)
            CALL H5SCREATE_SIMPLE_F(rank2, dims1(i,:), dspace_id, error)
        CASE (5:6)
            CALL H5SCREATE_SIMPLE_F(rank3, dims_vd, dspace_id, error)
        END SELECT
        CALL H5DCREATE_F(file_id, varnames1(i), H5T_REAL_PRECX, dspace_id, dset_id, error)
        SELECT CASE (i)
        CASE (1)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, t, dims1(i,:), error)
        CASE (2)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, VD_loc, dims1(i,:), error)
        CASE (3)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, VD_rho, dims1(i,:), error)
        CASE (4)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, VD_phase, dims1(i,:), error)
        CASE (5)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, VD_Krt, dims_vd, error)
        CASE (6)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, VD_Jrt, dims_vd, error)
        END SELECT
        CALL H5DCLOSE_F(dset_id, error)
        CALL H5SCLOSE_F(dspace_id, error)
    END DO
    CALL H5FCLOSE_F(file_id, error)
    END IF

    ! OUTPUT QUANTUM MECHANICS DATA !
    CALL H5FCREATE_F(filename2, H5F_ACC_EXCL_F, file_id, error)
    DO i=1,nvars2
        SELECT CASE (i)
        CASE (1:13)
            CALL H5SCREATE_SIMPLE_F(rank2, dims2(i,:), dspace_id, error)
        CASE (14:16)
            CALL H5SCREATE_SIMPLE_F(rank3, dims_psif, dspace_id, error)
        CASE (17:18)
            CALL H5SCREATE_SIMPLE_F(rank3, dims_psip, dspace_id, error)
        CASE (19:20)
            CALL H5SCREATE_SIMPLE_F(rank2, dims2(1,:), dspace_id, error)
        END SELECT
        CALL H5DCREATE_F(file_id, varnames2(i), H5T_REAL_PRECX, dspace_id, dset_id, error)
        SELECT CASE (i)
        CASE (1)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, t, dims2(i,:), error)
        CASE (2)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, tf, dims2(i,:), error)
        CASE (3)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, tp, dims2(i,:), error)
        CASE (4)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, x, dims2(i,:), error)
        CASE (5)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, y, dims2(i,:), error)
        CASE (6)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, fmsk, dims2(i,:), error)
        CASE (7)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, Ert, dims2(i,:), error)
        CASE (8)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, Art, dims2(i,:), error)
        CASE (9)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, Crt, dims2(i,:), error)
        CASE (10)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, S%kx, dims2(i,:), error)
        CASE (11)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, S%ky, dims2(i,:), error)
        CASE (12)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, REAL(S%p_dist,num), dims2(i,:), error)
        CASE (13)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, AIMAG(S%p_dist), dims2(i,:), error)
        CASE (14)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, REAL(psif,num), dims_psif, error)
        CASE (15)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, AIMAG(psif), dims_psif, error)
        CASE (16)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, Vp, dims_psif, error)
        CASE (17)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, REAL(psip,num), dims_psip, error)
        CASE (18)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, AIMAG(psip), dims_psip, error)
        CASE (19)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, wavefn%norm, dims2(1,:), error)
        CASE (20)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, wavefn%energy, dims2(1,:), error)
        END SELECT
        CALL H5DCLOSE_F(dset_id, error)
        CALL H5SCLOSE_F(dspace_id, error)
    END DO
    CALL H5FCLOSE_F(file_id, error)

    ! end HDF5 interface
    CALL H5CLOSE_F(error)

    RETURN
END SUBROUTINE WRITE_DATA_QCR_2D

SUBROUTINE WRITE_DATA_QCR_1D(t, x, VD, field, wavefn, tf, psi, Vp, fmsk, run_ccr)
    IMPLICIT NONE
    REAL(num),    INTENT(IN) :: t(:), x(:), tf(:), Vp(:,:), fmsk(:)
    COMPLEX(num), INTENT(IN) :: psi(:,:)
    LOGICAL,      INTENT(IN) :: run_ccr
    TYPE(vdet),   INTENT(IN) :: VD(:)
    TYPE(emf),    INTENT(IN) :: field
    TYPE(SchrodingerWavefunction1D), INTENT(IN) :: wavefn

    INTEGER :: nt, nx, ntf, N_VD

    REAL(num), DIMENSION(:),   ALLOCATABLE :: VD_loc
    REAL(num), DIMENSION(:,:), ALLOCATABLE :: VD_rho, VD_phase, VD_Krt, VD_Jrt

    INTEGER, PARAMETER :: nvars1 = 6, nvars2 = 11
    CHARACTER(LEN=10)  :: varnames1(nvars1), varnames2(nvars2)

    INTEGER :: n, i, j, error
    INTEGER, PARAMETER :: rank2 = 2, rank3 = 3
    INTEGER(HSIZE_T)   :: dims1(2,2), dims2(3,2), dims_vd(2), dims_psi(2)

    CHARACTER(LEN=10), PARAMETER :: filename1 = 'vd_data.h5'
    CHARACTER(LEN=10), PARAMETER :: filename2 = 'qm_data.h5'

    INTEGER(HID_T) :: file_id, dspace_id, dset_id
    INTEGER(HID_T) :: H5T_REAL_PRECX

    nt = SIZE(t)
    nx = SIZE(x)
    ntf = SIZE(tf)
    N_VD = SIZE(VD)

    varnames1(1) = 't';        dims1(1,:) = (/1,nt/)
    varnames1(2) = 'VD_loc';   dims1(2,:) = (/1,N_VD/)
    varnames1(3) = 'VD_rho'; ! dims_vd
    varnames1(4) = 'VD_phase'; dims_vd = (/N_VD,nt/)
    varnames1(5) = 'VD_Krt'; ! dims_vd
    varnames1(6) = 'VD_Jrt'; ! dims_vd

    varnames2(1)  = 't';            dims2(1,:) = (/1,nt/)
    varnames2(2)  = 'tf';           dims2(2,:) = (/1,ntf/)
    varnames2(3)  = 'x';            dims2(3,:) = (/1,nx/)
    varnames2(4)  = 'Et';         ! dims2(1,:)
    varnames2(5)  = 'At';         ! dims2(1,:)
    varnames2(6)  = 'mask';       ! dims2(3,:)
    varnames2(7)  = 'Re(psi)';      dims_psi = (/ntf,nx/)
    varnames2(8)  = 'Im(psi)';    ! dims_psi
    varnames2(9)  = 'Vp';         ! dims_psi
    varnames2(10) = 'norm(psi)';  ! dims2(1,:)
    varnames2(11) = 'enrg(psi)';  ! dims2(1,:)

    IF (run_ccr) THEN
    ALLOCATE(VD_loc(N_VD), VD_rho(N_VD,nt), VD_phase(N_VD,nt), VD_Krt(N_VD,nt), VD_Jrt(N_VD,nt))

    ! get virtual detector data
    DO n=1,N_VD
        ! locations
        VD_loc(n) = VD(n)%xl
        ! probability density
        VD_rho(n,:) = VD(n)%rho
        ! phase
        VD_phase(n,:) = VD(n)%phase
        ! momentum and current density
        VD_Krt(n,:) = VD(n)%Krt(1,:)
        VD_Jrt(n,:) = VD(n)%Jrt(1,:)
    END DO
    END IF

    ! initialize HDF5 interface
    CALL H5OPEN_F(error)

    ! set H5 precision
    IF (num.EQ.4) H5T_REAL_PRECX = H5T_IEEE_F32LE
    IF (num.EQ.8) H5T_REAL_PRECX = H5T_IEEE_F64LE

    ! OUTPUT VIRTUAL DETECTOR DATA !
    IF (run_ccr) THEN
    CALL H5FCREATE_F(filename1, H5F_ACC_EXCL_F, file_id, error)
    DO i=1,nvars1
        SELECT CASE (i)
        CASE (1:2)
            CALL H5SCREATE_SIMPLE_F(rank2, dims1(i,:), dspace_id, error)
        CASE (3:6)
            CALL H5SCREATE_SIMPLE_F(rank2, dims_vd, dspace_id, error)
        END SELECT
        CALL H5DCREATE_F(file_id, varnames1(i), H5T_REAL_PRECX, dspace_id, dset_id, error)
        SELECT CASE (i)
        CASE (1)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, t, dims1(i,:), error)
        CASE (2)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, VD_loc, dims1(i,:), error)
        CASE (3)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, VD_rho, dims_vd, error)
        CASE (4)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, VD_phase, dims_vd, error)
        CASE (5)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, VD_Krt, dims_vd, error)
        CASE (6)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, VD_Jrt, dims_vd, error)
        END SELECT
        CALL H5DCLOSE_F(dset_id, error)
        CALL H5SCLOSE_F(dspace_id, error)
    END DO
    CALL H5FCLOSE_F(file_id, error)
    END IF

    ! OUTPUT QUANTUM MECHANICS DATA !
    CALL H5FCREATE_F(filename2, H5F_ACC_EXCL_F, file_id, error)
    DO i=1,nvars2
        SELECT CASE (i)
        CASE (1:3)
            CALL H5SCREATE_SIMPLE_F(rank2, dims2(i,:), dspace_id, error)
        CASE (4:5)
            CALL H5SCREATE_SIMPLE_F(rank2, dims2(1,:), dspace_id, error)
        CASE (6)
            CALL H5SCREATE_SIMPLE_F(rank2, dims2(3,:), dspace_id, error)
        CASE (7:9)
            CALL H5SCREATE_SIMPLE_F(rank2, dims_psi, dspace_id, error)
        CASE (10:11)
            CALL H5SCREATE_SIMPLE_F(rank2, dims2(1,:), dspace_id, error)
        END SELECT
        CALL H5DCREATE_F(file_id, varnames2(i), H5T_REAL_PRECX, dspace_id, dset_id, error)
        SELECT CASE (i)
        CASE (1)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, t, dims2(i,:), error)
        CASE (2)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, tf, dims2(i,:), error)
        CASE (3)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, x, dims2(i,:), error)
        CASE (4)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, field%Ex, dims2(1,:), error)
        CASE (5)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, field%Ax, dims2(1,:), error)
        CASE (6)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, fmsk, dims2(3,:), error)
        CASE (7)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, REAL(psi,num), dims_psi, error)
        CASE (8)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, AIMAG(psi), dims_psi, error)
        CASE (9)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, Vp, dims_psi, error)
        CASE (10)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, wavefn%norm, dims2(1,:), error)
        CASE (11)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, wavefn%energy, dims2(1,:), error)
        END SELECT
        CALL H5DCLOSE_F(dset_id, error)
        CALL H5SCLOSE_F(dspace_id, error)
    END DO
    CALL H5FCLOSE_F(file_id, error)

    ! end HDF5 interface
    CALL H5CLOSE_F(error)

    RETURN
END SUBROUTINE WRITE_DATA_QCR_1D

SUBROUTINE WRITE_DATA_CCR(ED)
    IMPLICIT NONE
    TYPE(edet), INTENT(IN) :: ED

    INTEGER, PARAMETER :: nvars = 2
    CHARACTER(LEN=4) :: varnames(nvars)

    INTEGER :: i, error
    INTEGER, PARAMETER :: rank2 = 2
    INTEGER(HSIZE_T)   :: dims(nvars,2)

    CHARACTER(LEN=10) :: filename = 'ed_data.h5'

    INTEGER(HID_T) :: file_id, dspace_id, dset_id
    INTEGER(HID_T) :: H5T_REAL_PRECX

    varnames(1) = 'data'
    varnames(2) = 'bfwt'

    dims(1,:) = (/SIZE(ED%data,1),SIZE(ED%data,2)/)
    dims(2,:) = (/1,2/)

    ! initialize HDF5 interface and create file
    CALL H5OPEN_F(error)
    CALL H5FCREATE_F(filename, H5F_ACC_EXCL_F, file_id, error)

    ! set H5 precision
    IF (num.EQ.4) H5T_REAL_PRECX = H5T_IEEE_F32LE
    IF (num.EQ.8) H5T_REAL_PRECX = H5T_IEEE_F64LE

    ! OUTPUT END DETECTOR DATA !
    DO i=1,nvars
        CALL H5SCREATE_SIMPLE_F(rank2, dims(i,:), dspace_id, error)
        CALL H5DCREATE_F(file_id, varnames(i), H5T_REAL_PRECX, dspace_id, dset_id, error)
        SELECT CASE (i)
        CASE (1)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, ED%data, dims(i,:), error)
        CASE (2)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, ED%bfwt, dims(i,:), error)
        END SELECT
        CALL H5DCLOSE_F(dset_id, error)
        CALL H5SCLOSE_F(dspace_id, error)
    END DO

    ! close file and end HDF5 interface
    CALL H5FCLOSE_F(file_id, error)
    CALL H5CLOSE_F(error)

    RETURN
END SUBROUTINE WRITE_DATA_CCR

SUBROUTINE WRITE_DATA_CCR_TT(tt_data, ipt)
    IMPLICIT NONE
    REAL(num), INTENT(IN) :: tt_data(:,:)
    INTEGER, INTENT(IN) :: ipt(:)

    INTEGER, PARAMETER :: nvars = 8
    CHARACTER(LEN=2) :: varnames(nvars)

    INTEGER :: i, error
    INTEGER, PARAMETER :: rank2 = 2
    INTEGER(HSIZE_T)   :: dims(3,2)

    CHARACTER(LEN=10), PARAMETER :: filename = 'tt_data.h5'

    INTEGER(HID_T) :: file_id, dspace_id, dset_id
    INTEGER(HID_T) :: H5T_REAL_PRECX

    varnames(1) = 't';    dims(1,:) = (/1,SIZE(tt_data,2)/)
    varnames(2) = 'x';  ! dims(1,:)
    varnames(3) = 'y';  ! dims(1,:)
    varnames(4) = 'px'; ! dims(1,:)
    varnames(5) = 'py'; ! dims(1,:)
    varnames(6) = 'fz'; ! dims(1,:)
    varnames(7) = 'en';   dims(2,:) = (/2,SIZE(tt_data,2)/)
    varnames(8) = 'id';   dims(3,:) = (/1,SIZE(ipt)/)

    ! initialize HDF5 interface and create file
    CALL H5OPEN_F(error)
    CALL H5FCREATE_F(filename, H5F_ACC_EXCL_F, file_id, error)

    ! set H5 precision
    IF (num.EQ.4) H5T_REAL_PRECX = H5T_IEEE_F32LE
    IF (num.EQ.8) H5T_REAL_PRECX = H5T_IEEE_F64LE

    ! OUTPUT TRAJECTORY TRACKING DATA !
    DO i=1,nvars
        SELECT CASE (i)
        CASE DEFAULT
            CALL H5SCREATE_SIMPLE_F(rank2, dims(1,:), dspace_id, error)
            CALL H5DCREATE_F(file_id, varnames(i), H5T_REAL_PRECX, dspace_id, dset_id, error)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, tt_data(i,:), dims(1,:), error)
        CASE (7)
            CALL H5SCREATE_SIMPLE_F(rank2, dims(2,:), dspace_id, error)
            CALL H5DCREATE_F(file_id, varnames(i), H5T_REAL_PRECX, dspace_id, dset_id, error)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, tt_data(i:i+1,:), dims(2,:), error)
        CASE (8)
            CALL H5SCREATE_SIMPLE_F(rank2, dims(3,:), dspace_id, error)
            CALL H5DCREATE_F(file_id, varnames(i), H5T_NATIVE_INTEGER, dspace_id, dset_id, error)
            CALL H5DWRITE_F(dset_id, H5T_NATIVE_INTEGER, ipt, dims(3,:), error)
        END SELECT
        CALL H5DCLOSE_F(dset_id, error)
        CALL H5SCLOSE_F(dspace_id, error)
    END DO

    ! close file and end HDF5 interface
    CALL H5FCLOSE_F(file_id, error)
    CALL H5CLOSE_F(error)

    RETURN
END SUBROUTINE WRITE_DATA_CCR_TT

SUBROUTINE READ_DATA_CCR(VD_loc, VD_rho, VD_phase, VD_Krt, VD_Jrt, sdim)
    IMPLICIT NONE
    CHARACTER(LEN=10), PARAMETER :: filename = 'vd_data.h5'
    REAL(num), INTENT(OUT), DIMENSION(:,:),   ALLOCATABLE :: VD_loc, VD_rho, VD_phase
    REAL(num), INTENT(OUT), DIMENSION(:,:,:), ALLOCATABLE :: VD_Krt, VD_Jrt
    INTEGER, INTENT(IN) :: sdim

    ! HDF5 file variables
    INTEGER(HID_T) :: file_id, dspace_id, dset_id
    INTEGER(HID_T) :: H5T_REAL_PRECX
    INTEGER(HSIZE_T) :: dims2(2), dims3(3)
    INTEGER :: error

    ! initialize HDF5 interface and open file
    CALL H5OPEN_F(error)
    CALL H5FOPEN_F(filename, H5F_ACC_RDONLY_F, file_id, error)

    ! set H5 precision
    IF (num.EQ.4) H5T_REAL_PRECX = H5T_IEEE_F32LE
    IF (num.EQ.8) H5T_REAL_PRECX = H5T_IEEE_F64LE

    ! load VD location data
    CALL H5DOPEN_F(file_id, 'VD_loc', dset_id, error)
    CALL H5DGET_SPACE_F(dset_id, dspace_id, error)
    CALL H5SGET_SIMPLE_EXTENT_DIMS_F(dspace_id, dims2, dims2, error)
    ALLOCATE(VD_loc(dims2(1), dims2(2)))
    CALL H5DREAD_F(dset_id, H5T_REAL_PRECX, VD_loc, dims2, error)

    ! load VD probability density data
    CALL H5DOPEN_F(file_id, 'VD_rho', dset_id, error)
    CALL H5DGET_SPACE_F(dset_id, dspace_id, error)
    CALL H5SGET_SIMPLE_EXTENT_DIMS_F(dspace_id, dims2, dims2, error)
    ALLOCATE(VD_rho(dims2(1), dims2(2)))
    CALL H5DREAD_F(dset_id, H5T_REAL_PRECX, VD_rho, dims2, error)

    ! load VD phase data
    CALL H5DOPEN_F(file_id, 'VD_phase', dset_id, error)
    CALL H5DGET_SPACE_F(dset_id, dspace_id, error)
    CALL H5SGET_SIMPLE_EXTENT_DIMS_F(dspace_id, dims2, dims2, error)
    ALLOCATE(VD_phase(dims2(1), dims2(2)))
    CALL H5DREAD_F(dset_id, H5T_REAL_PRECX, VD_phase, dims2, error)

    ! load VD momentum data
    CALL H5DOPEN_F(file_id, 'VD_Krt', dset_id, error)
    CALL H5DGET_SPACE_F(dset_id, dspace_id, error)
    CALL H5SGET_SIMPLE_EXTENT_DIMS_F(dspace_id, dims3, dims3, error)
    IF (sdim.EQ.1) dims3(3) = 1
    ALLOCATE(VD_Krt(dims3(1), dims3(2), dims3(3)))
    CALL H5DREAD_F(dset_id, H5T_REAL_PRECX, VD_Krt, dims3, error)

    ! load VD current density data
    CALL H5DOPEN_F(file_id, 'VD_Jrt', dset_id, error)
    CALL H5DGET_SPACE_F(dset_id, dspace_id, error)
    CALL H5SGET_SIMPLE_EXTENT_DIMS_F(dspace_id, dims3, dims3, error)
    IF (sdim.EQ.1) dims3(3) = 1
    ALLOCATE(VD_Jrt(dims3(1), dims3(2), dims3(3)))
    CALL H5DREAD_F(dset_id, H5T_REAL_PRECX, VD_Jrt, dims3, error)

    ! close file and end HDF5 interface
    CALL H5FCLOSE_F(file_id, error)
    CALL H5CLOSE_F(error)

    RETURN
END SUBROUTINE READ_DATA_CCR

END MODULE DATAIO
