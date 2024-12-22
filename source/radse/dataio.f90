MODULE DATAIO

USE prec
USE hdf5
USE quantum
USE emfm

IMPLICIT NONE

CONTAINS

SUBROUTINE WRITE_DATA(t, r, field, wavefn, fmsk, tf, phif, E, W, S)
    IMPLICIT NONE
    REAL(num),    INTENT(IN) :: t(:), r(:), tf(:), fmsk(:), E(:), W(:,:)
    COMPLEX(num), INTENT(IN) :: S(:), phif(:,:,:)
    TYPE(emf),    INTENT(IN) :: field
    TYPE(SchrodingerWavefunction1DR), INTENT(IN) :: wavefn

    INTEGER :: nt, nr, ntf, l_max, npe_bins

    INTEGER, PARAMETER :: nvars = 14
    CHARACTER(LEN=10)  :: varnames(nvars)

    INTEGER :: n, i, j, error
    INTEGER, PARAMETER :: rank2 = 2, rank3 = 3
    INTEGER(HSIZE_T)   :: dims(6,2), dims_phif(3)

    CHARACTER(LEN=9), PARAMETER :: filename = 'output.h5'

    INTEGER(HID_T) :: file_id, dspace_id, dset_id
    INTEGER(HID_T) :: H5T_REAL_PRECX

    nt = SIZE(t)
    nr = SIZE(r)
    ntf = SIZE(tf)
    l_max = SIZE(wavefn%phi,1) - 1
    npe_bins = SIZE(E)

    varnames(1)  = 't';            dims(1,:) = (/1,nt/)
    varnames(2)  = 'tf';           dims(2,:) = (/1,ntf/)
    varnames(3)  = 'r';            dims(3,:) = (/1,nr/)
    varnames(4)  = 'mask';         dims(4,:) = (/1,nr/)
    varnames(5)  = 'E';            dims(5,:) = (/1,npe_bins/)
    varnames(6)  = 'dW_dE';        dims(6,:) = (/ntf,npe_bins/)
    varnames(7)  = 'Re(S)';      ! dims(1,:)
    varnames(8)  = 'Im(S)';      ! dims(1,:)
    varnames(9)  = 'E(t)';       ! dims(1,:)
    varnames(10) = 'A(t)';       ! dims(1,:)
    varnames(11) = 'norm(phi)';  ! dims(1,:)
    varnames(12) = 'enrg(phi)';  ! dims(1,:)
    varnames(13) = 'Re(phi)';      dims_phif  = (/l_max+1,nr,ntf/)
    varnames(14) = 'Im(phi)';    ! dims_phif

    ! initialize HDF5 interface
    CALL H5OPEN_F(error)

    ! set H5 precision
    IF (num.EQ.4) H5T_REAL_PRECX = H5T_IEEE_F32LE
    IF (num.EQ.8) H5T_REAL_PRECX = H5T_IEEE_F64LE

    ! OUTPUT QUANTUM MECHANICS DATA !
    CALL H5FCREATE_F(filename, H5F_ACC_EXCL_F, file_id, error)
    DO i=1,nvars
        SELECT CASE (i)
        CASE (1:6)
            CALL H5SCREATE_SIMPLE_F(rank2, dims(i,:), dspace_id, error)
        CASE (7:12)
            CALL H5SCREATE_SIMPLE_F(rank2, dims(1,:), dspace_id, error)
        CASE (13:14)
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
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, fmsk, dims(i,:), error)
        CASE (5)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, E, dims(i,:), error)
        CASE (6)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, W, dims(i,:), error)
        CASE (7)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, REAL(S,num), dims(1,:), error)
        CASE (8)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, AIMAG(S), dims(1,:), error)
        CASE (9)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, field%Ex, dims(1,:), error)
        CASE (10)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, field%Ax, dims(1,:), error)
        CASE (11)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, wavefn%norm, dims(1,:), error)
        CASE (12)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, wavefn%energy, dims(1,:), error)
        CASE (13)
            CALL H5DWRITE_F(dset_id, H5T_REAL_PRECX, REAL(phif,num), dims_phif, error)
        CASE (14)
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
