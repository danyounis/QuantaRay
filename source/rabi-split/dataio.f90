MODULE DATAIO

USE prec
USE hdf5
USE quantum
USE emfm

IMPLICIT NONE

CONTAINS

SUBROUTINE WRITE_DATA(t, r, pump, probe, wavefn, tf, phif, E, W, Pbt, Pft, Ec, lc)
    IMPLICIT NONE
    REAL(num),    INTENT(IN) :: t(:), r(:), tf(:), Ec(:), E(:), W(:,:), Pbt(:,:), Pft(:,:)
    INTEGER,      INTENT(IN) :: lc(:)
    COMPLEX(num), INTENT(IN) :: phif(:,:,:)
    TYPE(emf),    INTENT(IN) :: pump, probe
    TYPE(SchrodingerWavefunction1DR), INTENT(IN) :: wavefn

    INTEGER :: nt, nr, na(2), ntf, l_max, npe_bins

    INTEGER, PARAMETER :: nvars = 17
    CHARACTER(LEN=10)  :: varnames(nvars)

    INTEGER :: i, error
    INTEGER, PARAMETER :: rank2 = 2, rank3 = 3
    INTEGER(HSIZE_T)   :: dims(9,2), dims_phif(3)

    CHARACTER(LEN=9), PARAMETER :: filename = 'output.h5'

    INTEGER(HID_T) :: file_id, dspace_id, dset_id
    INTEGER(HID_T) :: H5T_REAL_PRECX

    nt = SIZE(t)
    nr = SIZE(r)
    ntf = SIZE(tf)
    npe_bins = SIZE(E)
    na(1) = size(Pft,1)
    na(2) = size(Pft,2)
    l_max = SIZE(wavefn%phi,1) - 1

    varnames(1)  = 't';            dims(1,:) = (/1,nt/)
    varnames(2)  = 'tf';           dims(2,:) = (/1,ntf/)
    varnames(3)  = 'r';            dims(3,:) = (/1,nr/)
    varnames(4)  = 'E';            dims(4,:) = (/1,npe_bins/)
    varnames(5)  = 'dW_dE';        dims(5,:) = (/ntf,npe_bins/)
    varnames(6)  = 'Pb(t)';        dims(6,:) = (/l_max+1,nt/)
    varnames(7)  = 'Pf(t)';        dims(7,:) = (/na(1),na(2)/)
    varnames(8)  = 'E(c)';         dims(8,:) = (/1,na(1)/)
    varnames(9)  = 'l(c)';         dims(9,:) = (/1,na(1)/)
    varnames(10) = 'E_u(t)';     ! dims(1,:)
    varnames(11) = 'A_u(t)';     ! dims(1,:)
    varnames(12) = 'E_r(t)';     ! dims(1,:)
    varnames(13) = 'A_r(t)';     ! dims(1,:)
    varnames(14) = 'norm(phi)';  ! dims(1,:)
    varnames(15) = 'enrg(phi)';  ! dims(1,:)
    varnames(16) = 'Re(phi)';      dims_phif  = (/l_max+1,nr,ntf/)
    varnames(17) = 'Im(phi)';    ! dims_phif

    ! initialize HDF5 interface
    CALL H5OPEN_F(error)

    ! set H5 precision
    IF (num.EQ.4) H5T_REAL_PRECX = H5T_IEEE_F32LE
    IF (num.EQ.8) H5T_REAL_PRECX = H5T_IEEE_F64LE

    ! OUTPUT QUANTUM MECHANICS DATA !
    CALL H5FCREATE_F(filename, H5F_ACC_EXCL_F, file_id, error)
    DO i=1,nvars
        SELECT CASE (i)
        CASE (1:9)
            CALL H5SCREATE_SIMPLE_F(rank2, dims(i,:), dspace_id, error)
        CASE (10:15)
            CALL H5SCREATE_SIMPLE_F(rank2, dims(1,:), dspace_id, error)
        CASE (16:17)
            CALL H5SCREATE_SIMPLE_F(rank3, dims_phif, dspace_id, error)
        END SELECT
        CALL H5DCREATE_F(file_id, varnames(i), H5T_REAL_PRECX, dspace_id, dset_id, error)
        SELECT CASE (i)
        CASE (1)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, t, dims(i,:), error)
        CASE (2)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, tf, dims(i,:), error)
        CASE (3)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, r, dims(i,:), error)
        CASE (4)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, E, dims(i,:), error)
        CASE (5)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, W, dims(i,:), error)
        CASE (6)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, Pbt, dims(i,:), error)
        CASE (7)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, Pft, dims(i,:), error)
        CASE (8)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, Ec, dims(i,:), error)
        CASE (9)
            CALL H5DWRITE_F(dset_id, H5T_STD_I32LE, lc, dims(i,:), error)
        CASE (10)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, pump%Ex, dims(1,:), error)
        CASE (11)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, pump%Ax, dims(1,:), error)
        CASE (12)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, probe%Ex, dims(1,:), error)
        CASE (13)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, probe%Ax, dims(1,:), error)
        CASE (14)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, wavefn%norm, dims(1,:), error)
        CASE (15)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, wavefn%energy, dims(1,:), error)
        CASE (16)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, REAL(phif,num), dims_phif, error)
        CASE (17)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, AIMAG(phif), dims_phif, error)
        END SELECT
        CALL H5DCLOSE_F(dset_id, error)
        CALL H5SCLOSE_F(dspace_id, error)
    END DO
    CALL H5FCLOSE_F(file_id, error)

    ! end HDF5 interface
    CALL H5CLOSE_F(error)

    RETURN
END SUBROUTINE WRITE_DATA

END MODULE DATAIO
