! w2_output.hdf5
! Branches 
! Segments
! Time Series
! Flow and Mass Balance
! WQ Pathway Fluxes
! Withdrawal Outlets  
!
! Zhonglong Zhang (PSU)
! Updated 1/2024
!===========================================================================================================================   
MODULE OUTPUTHDF5
  USE HDF5; USE HDF5MOD
  USE ISO_C_BINDING  
  USE MAIN; USE GLOBAL 
  USE GEOMC,     ONLY: DLX, ELWS, H, BI, B
  USE TVDC,      ONLY: NAC, NACD, NACTR, CONSTITUENTS
  USE SELWC,     ONLY: NSTR
  USE SURFHE,    ONLY: PHI0
  USE NAMESC,    ONLY: TITLE
  USE LOGICC,    ONLY: PRINT_SEDIMENT
  USE EDDY,      ONLY: FRIC
  USE SCREENC,   ONLY: MINDLT 
  USE KINETIC,   ONLY: NAF
  USE GDAYC
  USE CEMAVars,  ONLY: SEDIMENT_DIAGENESIS
  IMPLICIT NONE
  
  CHARACTER(LEN=15), PARAMETER  :: filename = "W2_output.hdf5" 
  INTEGER, PARAMETER            :: real_kind_15 = SELECTED_REAL_KIND(15,307)
  INTEGER, PARAMETER            :: rank1 = 1
  INTEGER, PARAMETER            :: rank2 = 2
  INTEGER, PARAMETER            :: rank3 = 3
  INTEGER, PARAMETER            :: rank4 = 4
  INTEGER                       :: error
  INTEGER                       :: FLX_I, PRF_I  
  INTEGER, ALLOCATABLE, DIMENSION(:) :: INI_KTWB
  INTEGER(HID_T) :: file_id     
  INTEGER(HID_T) :: group_id4, group_id4_1_2, group_id4_1_3, group_id4_1_3_1, group_id4_1_3_2, group_id4_1_3_2_1, group_id4_1_3_2_2, group_id4_1_3_2_3,   &
                    group_id4_1_4, group_id4_1_5, group_id4_1_7, group_id4_1_8, group_id4_1_10    
  INTEGER(HID_T) :: groupname4_1_3_1_1, groupname4_1_3_1_2, groupname4_1_3_1_3, groupname4_1_3_1_4, groupname4_1_3_1_5, groupname4_1_3_1_6
  INTEGER(HID_T), DIMENSION(99) :: group_id1_2_1, group_id4_1_3_1_1, group_id4_1_5_2
  !
  CHARACTER(50), TARGET :: dataset_name4_1_8a, dataset_name4_1_10a, dataset_name4_1_10b, dataset_name4_1_10c, dataset_name4_1_10d, dataset_name4_1_4a, dataset_name4_1_4b                                         
  CHARACTER(50), DIMENSION(99), TARGET  :: dataset_name4_1_2, dataset_name4_1_3_1ts, dataset_name4_1_3_1cs, dataset_name4_1_3_1ds, dataset_name4_1_3_1s, &
                                           dataset_name4_1_3_2_2t, dataset_name4_1_3_2_3c, dataset_name4_1_3_2_3d, dataset_name4_1_5a, &                                          
                                           dataset_name4_1_7a, dataset_name4_1_7b, dataset_name4_1_8, dataset_name4_1_3_2_1q
  CHARACTER(50), DIMENSION(999), TARGET :: dataset_name4_1_5_2  
  
  ! Dataset dimensions    
  INTEGER(HSIZE_T), DIMENSION(1)               :: maxdims1, adims1, adims2, adims3, adims4 
  INTEGER(HSIZE_T), DIMENSION(1:2)             :: maxdims2 
  INTEGER(HSIZE_T), DIMENSION(1:3)             :: maxdims3
  INTEGER(HSIZE_T), DIMENSION(1:4)             :: maxdims4
  INTEGER(HSIZE_T), DIMENSION(999,9999,1:2)    :: dimN
  INTEGER(HSIZE_T), DIMENSION(999,1:2)         :: dims
  INTEGER(HSIZE_T), DIMENSION(999,9999,1:3)    :: dimssX, dimssY, dimssZ
      
  ! Dataset parameters
  TYPE DSETPARS
    SEQUENCE
    INTEGER                          :: rank
    INTEGER(HSIZE_T), DIMENSION(1:1) :: inidim, dim, maxdim
    INTEGER(HSIZE_T), DIMENSION(1:2) :: inidims, dims, maxdims
    INTEGER(HSIZE_T), DIMENSION(1:3) :: inidimss, dimss, maxdimss
    INTEGER(HSIZE_T), DIMENSION(1:4) :: inidimsss, dimsss, maxdimsss
    INTEGER(HID_T)                   :: group
    INTEGER(HID_T)                   :: parent_group
    INTEGER(HID_T)                   :: dataset
    INTEGER(HID_T), DIMENSION(1:3)   :: timeset
    INTEGER(HID_T)                   :: dataspace
    INTEGER(HID_T), DIMENSION(1:3)   :: timespace
    INTEGER(HID_T)                   :: crplist
    INTEGER(HID_T)                   :: datatype
    CHARACTER(50), POINTER           :: name
    LOGICAL, POINTER                 :: status
    INTEGER(HID_T)                   :: aspace
    INTEGER(HID_T)                   :: aID
    INTEGER(HID_T)                   :: time_group
  END TYPE DSETPARS
      
  INTERFACE HDF5_OPEN_DATASET
    MODULE PROCEDURE HDF5_OPEN_DATASET1C
    MODULE PROCEDURE HDF5_OPEN_DATASET2
    MODULE PROCEDURE HDF5_OPEN_DATASET2C
    MODULE PROCEDURE HDF5_OPEN_DATASETS
  END INTERFACE
     
  INTERFACE HDF5_ADD_ATTRIBUTE
    MODULE PROCEDURE HDF5_ADD_ATTRIBUTEI
    MODULE PROCEDURE HDF5_ADD_ATTRIBUTEI1
    MODULE PROCEDURE HDF5_ADD_ATTRIBUTEI_JW
    MODULE PROCEDURE HDF5_ADD_ATTRIBUTEIS
    MODULE PROCEDURE HDF5_ADD_ATTRIBUTEC
    MODULE PROCEDURE HDF5_ADD_ATTRIBUTEC1
    MODULE PROCEDURE HDF5_ADD_ATTRIBUTEC_JW
    MODULE PROCEDURE HDF5_ADD_ATTRIBUTECS
    MODULE PROCEDURE HDF5_ADD_ATTRIBUTECS_LENGTH
    MODULE PROCEDURE HDF5_ADD_ATTRIBUTER
    MODULE PROCEDURE HDF5_ADD_ATTRIBUTER1
    MODULE PROCEDURE HDF5_ADD_ATTRIBUTERS
    MODULE PROCEDURE HDF5_ADD_ATTRIBUTET1
    MODULE PROCEDURE HDF5_ADD_ATTRIBUTET_JW
  END INTERFACE
     
  INTERFACE HDF5_ADD_DATA
    MODULE PROCEDURE HDF5_ADD_DATA1C
    MODULE PROCEDURE HDF5_ADD_DATA2
    MODULE PROCEDURE HDF5_ADD_DATA2I
    MODULE PROCEDURE HDF5_ADD_DATA2C
    MODULE PROCEDURE HDF5_ADD_DATA2N
    MODULE PROCEDURE HDF5_ADD_DATA3
  END INTERFACE
     
  INTERFACE HDF5_ADD_DATE
    MODULE PROCEDURE HDF5_ADD_DATE1
    MODULE PROCEDURE HDF5_ADD_DATE_JW
  END INTERFACE  
     
  INTERFACE XINT
    MODULE PROCEDURE XINT1
    MODULE PROCEDURE XINT2
  END INTERFACE
     
  CONTAINS
  
  !=========================================================================================================================== 
  ! SET REAL PRECISIONS
  SUBROUTINE XINT1(X, N, Y)
    IMPLICIT NONE  
    INTEGER     :: N, M
    REAL        :: X, Y
    
    IF (X<=-1.0/(10**N) .OR. X>=1.0/(10**N)) THEN
      IF (N/=0) THEN
        Y = (ANINT(X*10**N))/10**N
      ELSE
        Y = IFIX(X)
      END IF
    ELSE IF (X/= 0.0) THEN
      DO M = 1, 100
        IF(X*(10**M)>=1.0 .OR. X*(10**M)<=-1.0) EXIT
      END DO
      IF (M<=6) THEN 
        Y = ANINT(X*(10**(M+N)))/10**(M+N)
      ELSE IF (M>100) THEN
        Y = 0.0
      ELSE
        Y = X
      END IF
    ELSE
      Y = X
    END IF
    return
  END SUBROUTINE XINT1
     
  SUBROUTINE XINT2(X, N, Y)
    IMPLICIT NONE  
    INTEGER :: X, N
    REAL    :: Y
    Y = X
    return
  END SUBROUTINE XINT2
  
  !=========================================================================================================================== 
  ! GROUPS AND PARAMETERS
  SUBROUTINE OUTPUTHDF5INI
    IMPLICIT NONE
      
    ! GROUPS
    CHARACTER(LEN=11), PARAMETER :: groupname4 = "W2_H5output"   
    CHARACTER(LEN=8),  PARAMETER :: groupname4_1_2 = "Segments"
    !
    CHARACTER(LEN=18), PARAMETER :: groupname4_1_3 = "Withdrawal Outlets"
    CHARACTER(LEN=10), PARAMETER :: groupname4_1_3_1 = "Structures"
    CHARACTER(LEN=10), PARAMETER :: groupname4_1_3_1_1 = "Segment"
    CHARACTER(LEN=4),  PARAMETER :: groupname4_1_3_2_1 = "Flow"
    CHARACTER(LEN=11), PARAMETER :: groupname4_1_3_2_2 = "Temperature"
    CHARACTER(LEN=12), PARAMETER :: groupname4_1_3_2_3 = "Constituents" 
    !
    CHARACTER(LEN=22), PARAMETER :: groupname4_1_4 = "Flow and Mass Balances"
    !
    CHARACTER(LEN=20), PARAMETER :: groupname4_1_7 = "Water Quality Fluxes"

    CHARACTER(LEN=11), PARAMETER :: groupname4_1_8 = "Time Series"
    !
    CHARACTER(LEN=4),  PARAMETER :: groupname4_1_9_2_1 = "Flow"
    CHARACTER(LEN=11), PARAMETER :: groupname4_1_9_2_2 = "Temperature"
    CHARACTER(LEN=12), PARAMETER :: groupname4_1_9_2_3 = "Constituents" 
    !
    CHARACTER(LEN=8),  PARAMETER :: groupname4_1_10 = "Branches"
             
    ! PARAMETERS
    FLX_i = 0
    PRF_i = 0
    DO I = 1, 99
      WRITE(HDF5_JWBS(I), '(I2)') I
    END DO
    ALLOCATE(DATAS2(IMX, 9999))
    adims1 = (/1/)
    adims2 = (/2/)
    adims3 = (/3/)
    
    DO I = 1, 999
      dims(I, 1:2) = (/I, 1/)
      DO J = 1, 9999
        dimN(I, J, 1:2)   = (/I, J/)
        dimssX(I, J, 1:3) = (/1, I, J/)
        dimssY(I, J, 1:3) = (/I, 1, J/)
        dimssZ(I, J, 1:3) = (/I, J, 1/)
      END DO
    END DO
     
    ALLOCATE(HDF5_JTS(NWB))
    ALLOCATE(HDF5_NTR(NTR))
    ALLOCATE(HDF5_A(NWB))
    ALLOCATE(HDF5_MSEG(NWB))
    ALLOCATE(HDF5_WBNSEG(NWB))
    ALLOCATE(SEGMTFN(NWB), LAYRFN(NWB))
    ALLOCATE(HDF5_NBR(NBR))
    ALLOCATE(HDF5_KKSEG(KMX, NWB))
    ALLOCATE(INI_KTWB(NWB))
    
    HDF5_NN = 0
    DO JW = 1, NWB
      HDF5_JTS(JW) = 0
      DO JB = BS(JW), BE(JW)
        DO I = US(JB), DS(JB)
          HDF5_NN = HDF5_NN + 1
        END DO
      END DO
      DO K = 1, KMX
        HDF5_KKSEG(K, JW) = 0
        DO I = US(BS(JW))-1, DS(BE(JW))+1
          IF(B(K,I)>0.0) HDF5_KKSEG(K, JW) = HDF5_KKSEG(K, JW) + 1
        END DO
      END DO
      HDF5_WBNSEG(JW) = DS(BE(JW))-US(BS(JW)) + 3
      SEGMTFN(JW) = 3600 + JW-1
      LAYRFN(JW) = 37000 + (JW-1)*500
    END DO
      
    HDF5_KSEG = 0
    DO JT = 1, NTR
      HDF5_NTR(JT) = 2 + NACTR(JT)
      DO J = 1, NACTR(JT)
        HDF5_KSEG = HDF5_KSEG + 1
      END DO
    END DO
    DO JB = 1, NBR
      HDF5_NBR(JB) = 0
    END DO
    DO JW = 1, NWB
      IF (CONSTITUENTS) THEN
        HDF5_A(JW) = 18 + NAC+NEP+NMC+NACD(JW)+NAF(JW) + 3*NAL
        IF(SEDIMENT_CALC(JW)) HDF5_A(JW) = HDF5_A(JW) + 4
      ELSE
        HDF5_A(JW) = 18
      END IF
      IF(ICE_COMPUTATION)   HDF5_A(JW) = HDF5_A(JW) + 1
    END DO
    DO JW = 1, NWB
      IF (CONSTITUENTS) THEN
        HDF5_MSEG(JW) = 11 + NAC+NEP+NMC+NACD(JW)+NAF(JW) + 4*NAL + 4*NEP + 3*NMC
        IF(SEDIMENT_CALC(JW)) HDF5_MSEG(JW) = HDF5_MSEG(JW) + 4
      ELSE
        HDF5_MSEG(JW) = 11
      END IF
    END DO
    DO JW =1 , NWB
      DO J = 1, NIKTSR
        IF(ITSR(J)>= US(BS(JW)) .AND. ITSR(J)<= DS(BE(JW))) HDF5_JTS(JW) = HDF5_JTS(JW) + 1 
      END DO
    END DO
      
    HDF5_JWB = 1
    HDF5_NST = 0
    HDF5_N = 0
    !HDF5_M = 14 + NAC
    HDF5_M = 11
    DO J = 1, NBR
      IF(NSTR(J)>0) HDF5_NST = HDF5_NST + NSTR(J)
      HDF5_N = HDF5_N + DS(J)-US(J) + 1
    END DO
    IF(NWD>0) HDF5_NST = HDF5_NST + NWD
    IF(NSP>0) HDF5_NST = HDF5_NST + NSP
    IF(NPU>0) HDF5_NST = HDF5_NST + NPU
    IF(NPI>0) HDF5_NST = HDF5_NST + NPI
    IF(NGT>0) HDF5_NST = HDF5_NST + NGT
    
    IF(CONSTITUENTS) ALLOCATE(TEMPSETC2(HDF5_NST*NAC,1), TEMPSETD2(HDF5_NST*NACD(1),1))
          
    ! HDF5
    CALL h5open_f(error)
    CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)            
    CALL h5gcreate_f(file_id,       groupname4,       group_id4, error)
    CALL h5gcreate_f(group_id4,     groupname4_1_2,   group_id4_1_2, error)
    CALL h5gcreate_f(group_id4,     groupname4_1_3,   group_id4_1_3, error)
    DO I = 1, NIWDO
      WRITE(CHARFROMWD, '(I4)') IWDO(I)
      CALL h5gcreate_f(group_id4_1_3, groupname4_1_3_1_1//CHARFROMWD, group_id4_1_3_1_1(I), error)
    END DO   
    CALL h5gcreate_f(group_id4_1_3, groupname4_1_3_2_1, group_id4_1_3_2_1, error)
    CALL h5gcreate_f(group_id4_1_3, groupname4_1_3_2_2, group_id4_1_3_2_2, error)
    CALL h5gcreate_f(group_id4_1_3, groupname4_1_3_2_3, group_id4_1_3_2_3, error)
    CALL h5gcreate_f(group_id4,     groupname4_1_4,     group_id4_1_4, error)  
    !
    CALL h5gcreate_f(group_id4,     groupname4_1_7,     group_id4_1_7, error)
    CALL h5gcreate_f(group_id4,     groupname4_1_8,     group_id4_1_8, error)
    CALL h5gcreate_f(group_id4,     groupname4_1_10,    group_id4_1_10, error)   
  END SUBROUTINE OUTPUTHDF5INI
       
  SUBROUTINE  ENDOUTPUTHDF5
    IF(CONSTITUENTS) DEALLOCATE(TEMPSETC2, TEMPSETD2) 
    DEALLOCATE(DATAS2) 
    DEALLOCATE(HDF5_NTR, HDF5_A, HDF5_WBNSEG, HDF5_MSEG, HDF5_NBR, SEGMTFN, LAYRFN, HDF5_JTS, HDF5_KKSEG)
    
    CALL h5gclose_f(group_id4, error)
    CALL h5gclose_f(group_id4_1_2, error)
    CALL h5gclose_f(group_id4_1_3, error)
    DO I = 1, NIWDO
      CALL h5gclose_f(group_id4_1_3_1_1(I), error)
    END DO   
    CALL h5gclose_f(group_id4_1_3_2_1, error)
    CALL h5gclose_f(group_id4_1_3_2_2, error)
    CALL h5gclose_f(group_id4_1_3_2_3, error)
    CALL h5gclose_f(group_id4_1_4, error)
    !
    CALL h5gclose_f(group_id4_1_7, error)
    CALL h5gclose_f(group_id4_1_8, error)
    CALL h5gclose_f(group_id4_1_10, error)
    !
    CALL h5fclose_f(file_id, error)
    CALL h5close_f(error)
  END SUBROUTINE ENDOUTPUTHDF5

  !===========================================================================================================================   
  ! DATASET PARAMETERS
  SUBROUTINE FIND_ID(FILEID, SETNAME, PARX)
    USE CEMAVars, ONLY: SEDIMENT_DIAGENESIS  
    IMPLICIT NONE
    
    INTEGER(HID_T)    :: FILEID
    CHARACTER(*)      :: SETNAME
    TYPE(DSETPARS)    :: PARX
                    
    IF (HDF5_J>0 .AND. ALLOCATED(WDO)) THEN  
      IF (FILEID == WDO(HDF5_J,1)) THEN
        PARX%group      = group_id4_1_3_2_1
        PARX%time_group = group_id4_1_3
        PARX%rank = rank2
        PARX%inidims = dims(HDF5_NST+1,1:2)
        CALL H5Tcopy_f(H5T_IEEE_F64LE, PARX%datatype, error)
        IF(SETNAME /= 'N') dataset_name4_1_3_2_1q(HDF5_J) = ADJUSTL(SETNAME)
        PARX%name => dataset_name4_1_3_2_1q(HDF5_J)
        PARX%maxdims = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
        PARX%status => dataset_status4_1_3_2_1q(HDF5_J)
      END IF
            
      IF (FILEID ==  WDO(HDF5_J,2)) THEN
        PARX%group = group_id4_1_3_2_2
        PARX%rank = rank2
        PARX%inidims = dims(HDF5_NST+1,1:2)
        CALL H5Tcopy_f(H5T_IEEE_F64LE, PARX%datatype, error)
        IF(SETNAME /= 'N') dataset_name4_1_3_2_2t(HDF5_J) = ADJUSTL(SETNAME)
        PARX%name => dataset_name4_1_3_2_2t(HDF5_J)
        PARX%maxdims = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
        PARX%status => dataset_status4_1_3_2_2t(HDF5_J)
      END IF
      
IF (CONSTITUENTS) THEN     
      IF (FILEID ==  WDO(HDF5_J,3)) THEN
        PARX%group = group_id4_1_3_2_3
        PARX%rank = rank2
        PARX%inidims = dims(NAC*(1+HDF5_NST),1:2)
        PARX%inidims = dims(NAC,1:2)
        CALL H5Tcopy_f(H5T_IEEE_F64LE, PARX%datatype, error)
        IF(SETNAME /= 'N') dataset_name4_1_3_2_3c(HDF5_J) = ADJUSTL(SETNAME)
        PARX%name => dataset_name4_1_3_2_3c(HDF5_J)
        PARX%maxdims = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
        PARX%status => dataset_status4_1_3_2_3c(HDF5_J)
      END IF
    
      IF (FILEID ==  WDO(HDF5_J,4)) THEN
        PARX%group = group_id4_1_3_2_3
        PARX%rank = rank2
        PARX%inidims = dims(HDF5_NACD*(1+HDF5_NST),1:2)
        PARX%inidims = dims(HDF5_NACD,1:2)
        CALL H5Tcopy_f(H5T_IEEE_F64LE, PARX%datatype, error)
        IF(SETNAME /= 'N') dataset_name4_1_3_2_3d(HDF5_J) = ADJUSTL(SETNAME)
        PARX%name => dataset_name4_1_3_2_3d(HDF5_J)
        PARX%maxdims = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
        PARX%status => dataset_status4_1_3_2_3d(HDF5_J)
      END IF
END IF      
    
    END IF
        
    IF (JFILE > 0 .AND.  ALLOCATED(WDO2).AND.HDF5_JWD>0) THEN
      IF (FILEID == WDO2(JFILE,1)) THEN
        PARX%group = group_id4_1_3_1_1(HDF5_JWD)
        PARX%rank = rank2
        PARX%inidims = dims(2+NAC+HDF5_NACD,1:2)
        CALL H5Tcopy_f(H5T_IEEE_F64LE, PARX%datatype, error)
        IF(SETNAME /= 'N') dataset_name4_1_3_1s(JFILE) = ADJUSTL(SETNAME)
        PARX%name => dataset_name4_1_3_1s(JFILE)
        PARX%maxdims = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
        PARX%status => dataset_status4_1_3_1s(JFILE)
      END IF
    END IF
               
    ! BRANCH  
    IF (NBR>0) THEN
      ! FLOW
      IF (FILEID == BRANCHFN0 + 1) THEN
        PARX%group      = group_id4_1_10
        PARX%time_group = group_id4_1_10
        PARX%rank = rank3
        PARX%inidimss = dimssZ(NBR, KMX+NST+6, 1:3)
        CALL H5Tcopy_f(H5T_IEEE_F64LE, PARX%datatype, error)
        IF(SETNAME /= 'N') dataset_name4_1_10a = ADJUSTL(SETNAME)
        PARX%name => dataset_name4_1_10a
        PARX%maxdimss = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
        PARX%status => dataset_status4_1_10a
      END IF      
      ! VOLUME BALANCE
      IF (FILEID == BRANCHFN0 + 2) THEN
        PARX%group = group_id4_1_10
        PARX%rank = rank3
        PARX%inidimss = dimssZ(NBR, 4, 1:3)
        CALL H5Tcopy_f(H5T_IEEE_F64LE, PARX%datatype, error)
        IF(SETNAME /= 'N') dataset_name4_1_10b = ADJUSTL(SETNAME)
        PARX%name => dataset_name4_1_10b
        PARX%maxdimss = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
        PARX%status => dataset_status4_1_10b
      END IF   
      ! ENERGY BALNACE
      IF (FILEID == BRANCHFN0 + 3) THEN
        PARX%group = group_id4_1_10
        PARX%rank = rank3
        PARX%inidimss = dimssZ(NBR, 4, 1:3)
        CALL H5Tcopy_f(H5T_IEEE_F64LE, PARX%datatype, error)
        IF(SETNAME /= 'N') dataset_name4_1_10c = ADJUSTL(SETNAME)
        PARX%name => dataset_name4_1_10c
        PARX%maxdimss = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
        PARX%status => dataset_status4_1_10c
      END IF
      !
      IF (NMC+NAC>0) THEN
        IF (FILEID == BRANCHFN0 + 4) THEN
          PARX%group = group_id4_1_10
          PARX%rank = rank3 
          PARX%inidimss = dimssZ(NBR, 4*NAC+4*NMC, 1:3)
          CALL H5Tcopy_f(H5T_IEEE_F64LE, PARX%datatype, error)
          IF(SETNAME /= 'N') dataset_name4_1_10d = ADJUSTL(SETNAME)
          PARX%name => dataset_name4_1_10d
          PARX%maxdimss = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
          PARX%status => dataset_status4_1_10d
        END IF
      END IF
    END IF        
    
    ! SEGMENTS
    IF (HDF5_M>0 .AND. HDF5_WBNSEG(HDF5_JWB)>0) THEN
      IF (FILEID ==  SEGMTFN(HDF5_JWB)) THEN
        PARX%group      = group_id4_1_2
        PARX%time_group = group_id4_1_2
        PARX%rank = rank3
        PARX%inidimss = dimssZ(HDF5_M, HDF5_WBNSEG(HDF5_JWB), 1:3)
        CALL H5Tcopy_f(H5T_IEEE_F64LE, PARX%datatype, error)
        IF (SETNAME /= 'N') dataset_name4_1_2(HDF5_JWB) = ADJUSTL(SETNAME)
        PARX%name => dataset_name4_1_2(HDF5_JWB)
        PARX%maxdimss = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
        PARX%status => dataset_status4_1_2(HDF5_JWB)
      END IF
    END IF
    
    ! TIME SERIES TSR(J)
    IF (HDF5_A(HDF5_JWB)>0 .AND. HDF5_JTS(HDF5_JWB)>0) THEN
      IF (FILEID == TSRFN0+HDF5_JWB) THEN
        PARX%group = group_id4_1_8
        PARX%rank = rank3
        PARX%inidimss = dimssZ(HDF5_JTS(HDF5_JWB), HDF5_A(HDF5_JWB), 1:3)
        CALL H5Tcopy_f(H5T_IEEE_F64LE, PARX%datatype, error)
        IF(SETNAME /= 'N') dataset_name4_1_8(HDF5_JWB) = ADJUSTL(SETNAME)
        PARX%name => dataset_name4_1_8(HDF5_JWB) 
        PARX%maxdimss = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
        PARX%status => dataset_status4_1_8(HDF5_JWB) 
      END IF
    END IF
        
    ! TIME_SERIES WLFN
    IF (HDF5_N>0) THEN
      IF (FILEID == WLFN) THEN
        PARX%group      = group_id4_1_8
        PARX%time_group = PARX%group
        PARX%rank = rank2
        PARX%inidims = dims(HDF5_N, 1:2)
        CALL H5Tcopy_f(H5T_IEEE_F64LE, PARX%datatype, error)
        IF(SETNAME /= 'N') dataset_name4_1_8a = ADJUSTL(SETNAME)
        PARX%name => dataset_name4_1_8a
        PARX%maxdims = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
        PARX%status => dataset_status4_1_8a
      END IF
    END IF

    ! FLUX FLX(JW)
    IF (HDF5_JWB>0) THEN
      IF (FILEID == FLX(HDF5_JWB)) THEN
        PARX%group      = group_id4_1_7
        PARX%time_group = PARX%group 
        PARX%rank = rank3        
        PARX%inidimss = dimssX(NISNP(HDF5_JWB), KBR(HDF5_JWB)-KTWB(HDF5_JWB)+1, 1:3) 
        IF(NAF(HDF5_JWB)>0) PARX%inidimss = dimssX(NISNP(HDF5_JWB), NAF(HDF5_JWB)*(KBR(HDF5_JWB)-KTWB(HDF5_JWB)+1), 1:3)
        CALL H5Tcopy_f(H5T_IEEE_F64LE, PARX%datatype, error)
        IF(SETNAME /= 'N') dataset_name4_1_7a(HDF5_JWB) = ADJUSTL(SETNAME)
        PARX%name => dataset_name4_1_7a(HDF5_JWB)
        PARX%maxdimss = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
        PARX%status => dataset_status4_1_7a(HDF5_JWB)
      END IF
      
      ! FLUX OUTPUT FLX2(JW)
      IF (FILEID == FLX2(HDF5_JWB)) THEN
        PARX%group = group_id4_1_7
        PARX%rank = rank2
        PARX%inidims = dims(1+NAF(HDF5_JWB), 1:2)
        CALL H5Tcopy_f(H5T_IEEE_F64LE, PARX%datatype, error)
        IF(SETNAME /= 'N') dataset_name4_1_7b(HDF5_JWB) = ADJUSTL(SETNAME)
        PARX%name => dataset_name4_1_7b(HDF5_JWB)
        PARX%maxdims = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
        PARX%status => dataset_status4_1_7b(HDF5_JWB)
      END IF
    END IF 
            
    ! FLOW BALANCE
    IF (FILEID == FLOWBFN) THEN
      PARX%group      = group_id4_1_4
      PARX%time_group = PARX%group
      PARX%rank = rank2
      ALLOCATENO = 9
      IF(VOLUME_BALANCE(HDF5_JWB)) ALLOCATENO = 10
      PARX%inidims = dims(ALLOCATENO, 1:2)
      CALL H5Tcopy_f(H5T_IEEE_F64LE, PARX%datatype, error)
      IF(SETNAME /= 'N') dataset_name4_1_4a = ADJUSTL(SETNAME)
      PARX%name => dataset_name4_1_4a
      PARX%maxdims = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
      PARX%status => dataset_status4_1_4a
    END IF
         
    ! MASS BALANCE
    IF (FILEID == MASSBFN) THEN
      PARX%group = group_id4_1_4
      PARX%rank = rank2
      ALLOCATENO = 26
      IF(SEDIMENT_DIAGENESIS) ALLOCATENO = 29
      PARX%inidims = dims(ALLOCATENO, 1:2)
      CALL H5Tcopy_f(H5T_IEEE_F64LE, PARX%datatype, error)
      IF(SETNAME /= 'N') dataset_name4_1_4b = ADJUSTL(SETNAME)
      PARX%name => dataset_name4_1_4b
      PARX%maxdims = (/H5S_UNLIMITED_F, H5S_UNLIMITED_F/)
      PARX%status => dataset_status4_1_4b
    END IF 
  END SUBROUTINE FIND_ID
    
  !===========================================================================================================================   
  ! 1D STRING DATASET  
  SUBROUTINE HDF5_OPEN_DATASET1C(SETNO, SETNAME, LENS)
    IMPLICIT NONE  
    
    INTEGER            :: SETNO
    CHARACTER(*)       :: SETNAME
    CHARACTER          :: STATUS
    INTEGER            :: LENS
    TYPE(DSETPARS)     :: pard    
     
    CALL FIND_ID(SETNO, SETNAME, pard)
    IF (STATUS=='U') THEN
      CALL h5dopen_f(pard%group, pard%name, pard%dataset, error)
      IF (error==-1) THEN
        pard%status = .TRUE.
        CALL H5Tset_size_f(pard%datatype, LENS, error)
        CALL h5screate_simple_f(pard%rank, pard%inidim, pard%dataspace, error, pard%maxdim)
        CALL h5pcreate_f(H5P_DATASET_CREATE_F, pard%crplist, error)
        CALL h5pset_chunk_f(pard%crplist, pard%rank, pard%inidim, error)
        CALL h5pset_deflate_f(pard%crplist, 6, error) 
        CALL h5dcreate_f(pard%group, pard%name, pard%datatype, pard%dataspace, pard%dataset, error, pard%crplist)
        CALL h5pclose_f(pard%crplist, error)
        CALL h5dclose_f(pard%dataset, error)
        CALL h5sclose_f(pard%dataspace, error)
      END IF
    ELSE IF (STATUS=='A') THEN
      CALL h5dopen_f(pard%group, pard%name, pard%dataset, error)
    ELSE
      WRITE(w2err, '(A,I0)') 'ERROR-- WRONG HDF5_OPEN_TABLE INPUT STATUS PARAMETER'
      ERROR_OPEN = .TRUE.
    END IF
  END SUBROUTINE HDF5_OPEN_DATASET1C
  
  !=========================================================================================================================== 
  ! 2D STRING DATASET  
  SUBROUTINE HDF5_OPEN_DATASET2C(SETNO, SETNAME, LENS, STATUS, NEWSET)
    IMPLICIT NONE      
    
    INTEGER            :: SETNO
    CHARACTER(*)       :: SETNAME
    CHARACTER          :: STATUS
    INTEGER            :: LENS
    LOGICAL            :: NEWSET
    TYPE(DSETPARS)     :: pard    
     
    CALL FIND_ID(SETNO, SETNAME, pard)
    IF (STATUS=='U') THEN
      CALL h5dopen_f(pard%group, pard%name, pard%dataset, error)
      IF (error==-1) THEN
        pard%status = .TRUE.
        CALL H5Tset_size_f(pard%datatype, LENS, error)
        CALL h5screate_simple_f(pard%rank, pard%inidims, pard%dataspace, error, pard%maxdims)
        CALL h5pcreate_f(H5P_DATASET_CREATE_F, pard%crplist, error)
        CALL h5pset_chunk_f(pard%crplist, pard%rank, pard%inidims, error)
        CALL h5pset_deflate_f(pard%crplist, 6, error) 
        CALL h5dcreate_f(pard%group, pard%name, pard%datatype, pard%dataspace, pard%dataset, error, pard%crplist)
        CALL h5pclose_f(pard%crplist, error)
        CALL h5dclose_f(pard%dataset, error)
        CALL h5sclose_f(pard%dataspace, error)
      END IF
    ELSE IF (STATUS=='A') THEN
      CALL h5dopen_f(pard%group, pard%name, pard%dataset, error)
    ELSE
      WRITE(w2err, '(A,I0)') 'ERROR-- WRONG HDF5_OPEN_TABLE INPUT STATUS PARAMETER'
      ERROR_OPEN = .TRUE.
    END IF
  END SUBROUTINE HDF5_OPEN_DATASET2C
       
  SUBROUTINE HDF5_OPEN_DATASET2(SETNO, SETNAME, STATUS)
    IMPLICIT NONE    
    
    INTEGER            :: SETNO
    CHARACTER(*)       :: SETNAME
    CHARACTER          :: STATUS
    TYPE(DSETPARS)     :: pard    
     
    CALL FIND_ID(SETNO, SETNAME, pard)
    IF (STATUS=='U') THEN
      CALL h5dopen_f(pard%group, pard%name, pard%dataset, error)
      IF (error==-1) THEN
        pard%status = .TRUE.
        CALL h5screate_simple_f(pard%rank, pard%inidims, pard%dataspace, error, pard%maxdims)
        CALL h5pcreate_f(H5P_DATASET_CREATE_F, pard%crplist, error)
        CALL h5pset_chunk_f(pard%crplist, pard%rank, pard%inidims, error)
        CALL h5pset_deflate_f(pard%crplist, 6, error) 
        CALL h5dcreate_f(pard%group, pard%name, pard%datatype, pard%dataspace, pard%dataset, error, pard%crplist)
        CALL h5pclose_f(pard%crplist, error)
        CALL h5dclose_f(pard%dataset, error)
        CALL h5sclose_f(pard%dataspace, error)
      END IF
    ELSE IF (STATUS == 'A') THEN
      CALL h5dopen_f(pard%group, pard%name, pard%dataset, error)
    ELSE
      WRITE(w2err, '(A,I0)') 'ERROR-- WRONG HDF5_OPEN_TABLE INPUT STATUS PARAMETER'
      ERROR_OPEN = .TRUE.
    END IF
  END SUBROUTINE HDF5_OPEN_DATASET2
  
  !=========================================================================================================================== 
  ! 3D DATASET
  SUBROUTINE HDF5_OPEN_DATASETS(SETNO, SETNAME, STATUS, DATADIMS)
    IMPLICIT NONE       
    
    INTEGER            :: SETNO, DATADIMS
    CHARACTER(*)       :: SETNAME
    CHARACTER          :: STATUS
    TYPE(DSETPARS)     :: pard    
     
    CALL FIND_ID(SETNO, SETNAME, pard)
    IF (STATUS=='U') THEN
      CALL h5dopen_f(pard%group, pard%name, pard%dataset, error)
      IF (error==-1) THEN
        pard%status = .TRUE.
        IF(DATADIMS==3) CALL h5screate_simple_f(pard%rank, pard%inidimss, pard%dataspace, error, pard%maxdimss)
        IF(DATADIMS==4) CALL h5screate_simple_f(pard%rank, pard%inidimsss, pard%dataspace, error, pard%maxdimsss)
        CALL h5pcreate_f(H5P_DATASET_CREATE_F, pard%crplist, error)
        IF(DATADIMS==3) CALL h5pset_chunk_f(pard%crplist, pard%rank, pard%inidimss, error)
        IF(DATADIMS==4) CALL h5pset_chunk_f(pard%crplist, pard%rank, pard%inidimsss, error)
        CALL h5pset_deflate_f(pard%crplist, 6, error) 
        CALL h5dcreate_f(pard%group, pard%name, pard%datatype, pard%dataspace, pard%dataset, error, pard%crplist)
        CALL h5pclose_f(pard%crplist, error)
        CALL h5dclose_f(pard%dataset, error)
        CALL h5sclose_f(pard%dataspace, error)
      END IF
    ELSE IF (STATUS=='A') THEN
      CALL h5dopen_f(pard%group, pard%name, pard%dataset, error)
    ELSE
      WRITE(w2err, '(A,I0)') 'ERROR-- WRONG HDF5_OPEN_TABLE INPUT STATUS PARAMETER'
      ERROR_OPEN = .TRUE.
    END IF
  END SUBROUTINE HDF5_OPEN_DATASETS
   
  !=========================================================================================================================== 
  ! DATE/TIME
  SUBROUTINE HDF5_ADD_DATE1(TABNO, DAY) 
    IMPLICIT NONE 
    
    INTEGER                               :: TABNO
    REAL                                  :: JDAY, DAY
    INTEGER(HSIZE_T), DIMENSION(1:1)      :: offset, data_size
    INTEGER(HID_T)                        :: stamptype, memspace, dim0
    CHARACTER(19), POINTER                :: STAMP1 
    CHARACTER(19), TARGET                 :: STAMP11(1)    
    TYPE(DSETPARS)                        :: pard
    TYPE(C_PTR)                           :: f_ptr
    INTEGER(HSIZE_T)                      :: DIM1
    
    STAMP1 => STAMP11(1)     
    CALL FIND_ID(TABNO, 'N', pard)
    data_size(1:1) = (/1/)
    pard%inidim(1:1) = (/1/)
    pard%maxdim = (/H5S_UNLIMITED_F/)
    CALL XINT(DAY, 3, JDAY)
    CALL JDAYTOTS(DAY, STAMP1)
    
    CALL h5dopen_f(pard%time_group, 'Time', pard%timeset(1), error)  
    IF (error==-1) THEN
      ! Time
      CALL h5screate_simple_f(1, pard%inidim, pard%timespace(1), error, pard%maxdim)
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, pard%crplist, error)
      CALL h5pset_chunk_f(pard%crplist, 1,  pard%inidim, error)
      CALL h5pset_deflate_f(pard%crplist, 6, error) 
      CALL h5dcreate_f(pard%time_group, 'Time', H5T_IEEE_F64LE, pard%timespace(1), pard%timeset(1), error, pard%crplist)
      CALL h5dwrite_f(pard%timeset(1), H5T_IEEE_F64LE, JDAY, data_size, error)
      CALL h5pclose_f(pard%crplist, error)
      CALL h5sclose_f(pard%timespace(1), error)
      ! 
      CALL h5screate_simple_f(1, adims1, pard%aspace, error)
      CALL h5screate_f(H5S_SCALAR_F, pard%aspace, error)
      CALL h5acreate_f(pard%timeset(1), 'Number of actual Time Steps', H5T_NATIVE_INTEGER, pard%aspace, pard%aID, error)
      CALL h5awrite_f(pard%aID, H5T_NATIVE_INTEGER, 1, adims1, error)
      CALL h5aclose_f(pard%aID, error)
      CALL h5sclose_f(pard%aspace, error)
      CALL h5dclose_f(pard%timeset(1), error)  
    ELSE
      CALL h5dget_space_f(pard%timeset(1), pard%timespace(1), error)
      CALL h5sget_simple_extent_dims_f(pard%timespace(1), pard%dim, pard%maxdim, error) 
      dim0 = pard%dim(1)
      pard%dim(1:1) = pard%dim(1:1) + 1
      offset(1:1) = (/dim0/)
      CALL h5dset_extent_f(pard%timeset(1), pard%dim, error)
      CALL h5screate_simple_f (1, data_size, memspace, error) 
      CALL h5dget_space_f(pard%timeset(1), pard%timespace(1), error)
      CALL h5sselect_hyperslab_f(pard%timespace(1), H5S_SELECT_SET_F, offset, data_size, error)
      CALL h5dwrite_f(pard%timeset(1), H5T_IEEE_F64LE, JDAY, data_size, error, memspace, pard%timespace(1))
      CALL h5pclose_f(pard%crplist, error)
      CALL h5sclose_f(memspace, error)
      ! 
      CALL h5aopen_f(pard%timeset(1), 'Number of actual Time Steps', pard%aID, error)
      CALL h5awrite_f(pard%aID, H5T_NATIVE_INTEGER, DIM0+1, adims1, error)
      CALL h5aclose_f(pard%aID, error)
      CALL h5sclose_f(pard%timespace(1), error)
      CALL h5dclose_f(pard%timeset(1), error) 
    END IF

    ! Time Date Stamp
    CALL H5Tcopy_f(H5T_FORTRAN_S1, stamptype, error)
    CALL H5Tset_size_f(stamptype, 19, error)
    CALL h5dopen_f(pard%time_group, 'Time Data Stamp', pard%timeset(2), error)
    IF (error == -1) THEN
      CALL h5screate_simple_f(1, pard%inidim, pard%timespace(2), error, pard%maxdim)
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, pard%crplist, error)
      CALL h5pset_chunk_f(pard%crplist, 1, pard%inidim, error)
      CALL h5pset_deflate_f(pard%crplist, 6, error) 
      CALL h5dcreate_f(pard%time_group, 'Time Data Stamp', stamptype, pard%timespace(2), pard%timeset(2), error, pard%crplist)            
      f_ptr = C_LOC(STAMP1(1:1))
      CALL h5dwrite_f(pard%timeset(2), stamptype, f_ptr, error)
      CALL h5pclose_f(pard%crplist, error)
      CALL h5dclose_f(pard%timeset(2), error)
      CALL h5sclose_f(pard%timespace(2), error)      
    ELSE 
      CALL h5dget_space_f(pard%timeset(2), pard%timespace(2), error)
      CALL h5sget_simple_extent_dims_f(pard%timespace(2), pard%dim, pard%maxdim, error)
      dim0 = pard%dim(1)
      pard%dim(1:1) = pard%dim(1:1) + 1
      offset(1:1) = (/dim0/)
      CALL h5dset_extent_f(pard%timeset(2), pard%dim, error)
      CALL h5screate_simple_f (1, data_size, memspace, error) 
      CALL h5dget_space_f(pard%timeset(2), pard%timespace(2), error)
      CALL h5sselect_hyperslab_f(pard%timespace(2), H5S_SELECT_SET_F, offset, data_size, error)
      f_ptr = C_LOC(STAMP1(1:1))
      CALL h5dwrite_f(pard%timeset(2), stamptype, f_ptr, error, memspace, pard%timespace(2))
      CALL h5sclose_f(pard%timespace(2), error)
      CALL h5sclose_f(memspace, error)
      CALL h5dclose_f(pard%timeset(2), error)
    END IF
    
    CALL H5Tclose_f(stamptype, error)
  END SUBROUTINE HDF5_ADD_DATE1   
     
  SUBROUTINE HDF5_ADD_DATE_JW(JW, TABNO, DAY) 
    IMPLICIT NONE    
    
    INTEGER                               :: JW, TABNO
    REAL                                  :: JDAY, DAY
    INTEGER(HSIZE_T), DIMENSION(1:1)      :: offset, data_size
    INTEGER(HID_T)                        :: stamptype, memspace, dim0
    CHARACTER(19), POINTER                :: STAMP1 
    CHARACTER(19), TARGET                 :: STAMP11(1)   
    TYPE(DSETPARS)                        :: pard
    TYPE(C_PTR)                           :: f_ptr  
    
    STAMP1 => STAMP11(1)
    CALL FIND_ID(TABNO, 'N', pard)
    data_size(1:1) = (/1/)
    pard%inidim(1:1) = (/1/)
    pard%maxdim = (/H5S_UNLIMITED_F/)
    CALL XINT(DAY, 3, JDAY)
    CALL JDAYTOTS(DAY, STAMP1)
    
    CALL h5dopen_f(pard%time_group, 'Time'//HDF5_JWBS(JW), pard%timeset(1), error)   
    IF (error==-1) THEN
      ! Time
      CALL h5screate_simple_f(1, pard%inidim, pard%timespace(1), error, pard%maxdim)
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, pard%crplist, error)
      CALL h5pset_chunk_f(pard%crplist, 1,  pard%inidim, error)
      CALL h5pset_deflate_f(pard%crplist, 6, error) 
      CALL h5dcreate_f(pard%time_group, 'Time'//HDF5_JWBS(JW), H5T_IEEE_F64LE, pard%timespace(1), pard%timeset(1), error, pard%crplist)
      CALL h5dwrite_f(pard%timeset(1), H5T_IEEE_F64LE, JDAY, data_size, error)
      CALL h5pclose_f(pard%crplist, error)
      CALL h5sclose_f(pard%timespace(1), error)   
      ! 
      CALL h5screate_simple_f(1, adims1, pard%aspace, error)
      CALL h5screate_f(H5S_SCALAR_F, pard%aspace, error)
      CALL h5acreate_f(pard%timeset(1), 'Number of actual Time Steps', H5T_NATIVE_INTEGER, pard%aspace, pard%aID, error)
      CALL h5awrite_f(pard%aID, H5T_NATIVE_INTEGER, 1, adims1, error)
      CALL h5aclose_f(pard%aID, error)
      CALL h5sclose_f(pard%aspace, error)
      CALL h5dclose_f(pard%timeset(1), error)  
    ELSE
      CALL h5dget_space_f(pard%timeset(1), pard%timespace(1), error)
      CALL h5sget_simple_extent_dims_f(pard%timespace(1), pard%dim, pard%maxdim, error) 
      dim0 = pard%dim(1)
      pard%dim(1:1) = pard%dim(1:1) + 1
      offset(1:1) = (/dim0/)
      CALL h5dset_extent_f(pard%timeset(1), pard%dim, error)
      CALL h5screate_simple_f(1, data_size, memspace, error) 
      CALL h5dget_space_f(pard%timeset(1), pard%timespace(1), error)
      CALL h5sselect_hyperslab_f(pard%timespace(1), H5S_SELECT_SET_F, offset, data_size, error)
      CALL h5dwrite_f(pard%timeset(1), H5T_IEEE_F64LE, JDAY, data_size, error, memspace, pard%timespace(1))
      CALL h5pclose_f(pard%crplist, error)
      CALL h5sclose_f(memspace, error)
      ! 
      CALL h5aopen_f(pard%timeset(1), 'Number of actual Time Steps', pard%aID, error)
      CALL h5awrite_f(pard%aID, H5T_NATIVE_INTEGER, DIM0+1, adims1, error)
      CALL h5aclose_f(pard%aID, error)
      CALL h5sclose_f(pard%timespace(1), error)
      CALL h5dclose_f(pard%timeset(1), error)
    END IF

    ! Time Date Stamp
    CALL H5Tcopy_f(H5T_FORTRAN_S1, stamptype, error)
    CALL H5Tset_size_f(stamptype, 19, error)
    CALL h5dopen_f(pard%time_group, 'Time Data Stamp'//HDF5_JWBS(JW), pard%timeset(2), error)
    IF (error == -1) THEN
      CALL h5screate_simple_f(1, pard%inidim, pard%timespace(2), error, pard%maxdim)
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, pard%crplist, error)
      CALL h5pset_chunk_f(pard%crplist, 1, pard%inidim, error)
      CALL h5pset_deflate_f(pard%crplist, 6, error) 
      CALL h5dcreate_f(pard%time_group, 'Time Data Stamp'//HDF5_JWBS(JW), stamptype, pard%timespace(2), pard%timeset(2), error, pard%crplist)            
      f_ptr = C_LOC(STAMP1(1:1))
      CALL h5dwrite_f(pard%timeset(2), stamptype, f_ptr, error)
      CALL h5pclose_f(pard%crplist, error)
      CALL h5dclose_f(pard%timeset(2), error)
      CALL h5sclose_f(pard%timespace(2), error)    
    ELSE 
      CALL h5dget_space_f(pard%timeset(2), pard%timespace(2), error)
      CALL h5sget_simple_extent_dims_f(pard%timespace(2), pard%dim, pard%maxdim, error)
      dim0 = pard%dim(1)
      pard%dim(1:1) = pard%dim(1:1) + 1
      offset(1:1) = (/dim0/)
      CALL h5dset_extent_f(pard%timeset(2), pard%dim, error)
      CALL h5screate_simple_f (1, data_size, memspace, error) 
      CALL h5dget_space_f(pard%timeset(2), pard%timespace(2), error)
      CALL h5sselect_hyperslab_f(pard%timespace(2), H5S_SELECT_SET_F, offset, data_size, error)
      f_ptr = C_LOC(STAMP1(1:1))
      CALL h5dwrite_f(pard%timeset(2), stamptype, f_ptr, error, memspace, pard%timespace(2))
      CALL h5sclose_f(pard%timespace(2), error)
      CALL h5sclose_f(memspace, error)
      CALL h5dclose_f(pard%timeset(2), error)
    END IF
    
    CALL H5Tclose_f(stamptype, error)
  END SUBROUTINE HDF5_ADD_DATE_JW  

  !===========================================================================================================================   
  ! 1D CHARACTER DATA
  SUBROUTINE HDF5_ADD_DATA1C(TABNO, DATAX, DATAS, LENS)  
    IMPLICIT NONE      
    
    INTEGER                               :: TABNO, LENS
    INTEGER                               :: DATAX
    CHARACTER(LENS), POINTER              :: DATAS
    CHARACTER(LENS), TARGET               :: DATAS0(DATAX)
    INTEGER(HSIZE_T), DIMENSION(1:1)      :: offset, data_size
    INTEGER(HID_T)                        :: memspace, dim0
    TYPE(DSETPARS)                        :: pard
    TYPE(C_PTR)                           :: f_ptr
     
    DATAS => DATAS0(1)
    data_size(1:1) = (/1/)
    CALL FIND_ID(TABNO, 'N', pard)
    CALL H5Tset_size_f(pard%datatype, LENS, error)
    
    CALL h5dopen_f(pard%group, pard%name, pard%dataset, error)    
    IF (error==-1) THEN
      CALL h5screate_simple_f(pard%rank, pard%inidim, pard%dataspace, error, pard%maxdim)
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, pard%crplist, error)
      CALL h5pset_chunk_f(pard%crplist, pard%rank, pard%inidim, error)
      CALL h5pset_deflate_f(pard%crplist, 6, error) 
      CALL h5dcreate_f(pard%group, pard%name, pard%datatype, pard%dataspace, pard%dataset, error, pard%crplist) 
      f_ptr = C_LOC(DATAS(1:1))
      CALL h5dwrite_f(pard%dataset, pard%datatype, f_ptr, error)
      CALL h5pclose_f(pard%crplist, error)
      CALL h5dclose_f(pard%dataset, error)
      CALL h5sclose_f(pard%dataspace, error)    
    ELSE    
      CALL h5dget_space_f(pard%dataset, pard%dataspace, error)
      CALL h5sget_simple_extent_dims_f(pard%dataspace, pard%dim, pard%maxdim, error)
      dim0 = pard%dim(1)
      pard%dim(1:1) = pard%dim(1:1) + 1
      offset(1:1) = (/dim0/)
      CALL h5dset_extent_f(pard%dataset, pard%dim, error)
      CALL h5screate_simple_f (1, data_size, memspace, error) 
      CALL h5dget_space_f(pard%dataset, pard%dataspace, error)
      CALL h5sselect_hyperslab_f(pard%dataspace, H5S_SELECT_SET_F, offset, data_size, error)
      f_ptr = C_LOC(DATAS(1:1))
      CALL h5dwrite_f(pard%dataset, pard%datatype, f_ptr, error, memspace, pard%dataspace)
      CALL h5sclose_f(pard%dataspace, error)
      CALL h5sclose_f(memspace, error)
      CALL h5dclose_f(pard%dataset, error)
    END IF
  END SUBROUTINE HDF5_ADD_DATA1C
  
  !=========================================================================================================================== 
  ! 2D DATA
  SUBROUTINE HDF5_ADD_DATA2(TABNO,  DATAX, DATAY, DATAS)
    IMPLICIT NONE  
    
    INTEGER                               :: TABNO
    INTEGER                               :: DATAX, DATAY
    REAL, DIMENSION(DATAX,DATAY)          :: DATAS
    TYPE(DSETPARS)                        :: pard
    INTEGER(HSIZE_T), DIMENSION(1:2)      :: offset, data_size
    INTEGER(HID_T)                        :: memspace
    INTEGER(HSIZE_T)                      :: dim1, dim2
    
    CALL FIND_ID(TABNO, 'N', pard)
    data_size(1:2)  = (/DATAX, DATAY/)
    dim1 = 0
    
    CALL h5dopen_f(pard%group, pard%name, pard%dataset, error)  
    IF (error==-1) THEN
      pard%status = .TRUE.
      CALL h5screate_simple_f(pard%rank, pard%inidims, pard%dataspace, error, pard%maxdims)
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, pard%crplist, error)
      CALL h5pset_chunk_f(pard%crplist, pard%rank, pard%inidims, error)
      CALL h5pset_deflate_f(pard%crplist, 6, error) 
      CALL h5dcreate_f(pard%group, pard%name, pard%datatype, pard%dataspace, pard%dataset, error, pard%crplist)
      CALL h5dwrite_f(pard%dataset, pard%datatype, DATAS, data_size, error)
      CALL h5pclose_f(pard%crplist, error)
      CALL h5dclose_f(pard%dataset, error)
      CALL h5sclose_f(pard%dataspace, error)  
    ELSE
      IF (pard%status) THEN
        CALL h5dwrite_f(pard%dataset, pard%datatype, DATAS, data_size, error)
        CALL h5dclose_f(pard%dataset, error)
        pard%status = .FALSE.
      ELSE
        CALL h5dget_space_f(pard%dataset, pard%dataspace, error)
        CALL h5sget_simple_extent_dims_f(pard%dataspace, pard%dims, pard%maxdims, error)
        dim2 = pard%dims(2)
        offset(1:2) = (/dim1, dim2/)
        pard%dims(2) = pard%dims(2) + DATAY 
        CALL h5dset_extent_f(pard%dataset, pard%dims, error)
        CALL h5screate_simple_f (pard%rank, data_size, memspace, error) 
        CALL h5dget_space_f(pard%dataset, pard%dataspace, error)
        CALL h5sselect_hyperslab_f(pard%dataspace, H5S_SELECT_SET_F, offset, data_size, error)
        CALL h5dwrite_f(pard%dataset, pard%datatype, DATAS, data_size, error, memspace, pard%dataspace)
        CALL h5sclose_f(pard%dataspace, error)
        CALL h5sclose_f(memspace, error)
        CALL h5dclose_f(pard%dataset, error)
      END IF
    END IF
  END SUBROUTINE HDF5_ADD_DATA2
     
  ! 2D INTEGER DATA
  SUBROUTINE HDF5_ADD_DATA2I(TABNO, DATAX, DATAY, DATAS)
    IMPLICIT NONE    
    
    INTEGER                               :: TABNO
    INTEGER                               :: DATAX, DATAY
    INTEGER, DIMENSION(DATAX,DATAY)       :: DATAS
    TYPE(DSETPARS)                        :: pard
    INTEGER(HSIZE_T), DIMENSION(1:2)      :: offset, data_size
    INTEGER(HID_T)                        :: memspace
    INTEGER(HSIZE_T)                      :: dim1, dim2
    
    CALL FIND_ID(TABNO, 'N', pard)
    data_size(1:2)  = (/DATAX, DATAY/)
    dim1 = 0
    
    CALL h5dopen_f(pard%group, pard%name, pard%dataset, error)  
    IF (error==-1) THEN
      pard%status = .TRUE.
      CALL h5screate_simple_f(pard%rank, pard%inidims, pard%dataspace, error, pard%maxdims)
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, pard%crplist, error)
      CALL h5pset_chunk_f(pard%crplist, pard%rank, pard%inidims, error)
      CALL h5pset_deflate_f(pard%crplist, 6, error) 
      CALL h5dcreate_f(pard%group, pard%name, pard%datatype, pard%dataspace, pard%dataset, error, pard%crplist)
      CALL h5dwrite_f(pard%dataset, pard%datatype, DATAS, data_size, error)
      CALL h5pclose_f(pard%crplist, error)
      CALL h5dclose_f(pard%dataset, error)
      CALL h5sclose_f(pard%dataspace, error)     
    ELSE
      IF (pard%status) THEN
        CALL h5dwrite_f(pard%dataset, pard%datatype, DATAS, data_size, error)
        CALL h5dclose_f(pard%dataset, error)
        pard%status = .FALSE.
      ELSE
        CALL h5dget_space_f(pard%dataset, pard%dataspace, error)
        CALL h5sget_simple_extent_dims_f(pard%dataspace, pard%dims, pard%maxdims, error)
        dim2 = pard%dims(2)
        offset(1:2) = (/dim1, dim2/)
        pard%dims(2) = pard%dims(2) + DATAY    
        CALL h5dset_extent_f(pard%dataset, pard%dims, error)
        CALL h5screate_simple_f(pard%rank, data_size, memspace, error) 
        CALL h5dget_space_f(pard%dataset, pard%dataspace, error)
        CALL h5sselect_hyperslab_f(pard%dataspace, H5S_SELECT_SET_F, offset, data_size, error)
        CALL h5dwrite_f(pard%dataset, pard%datatype, DATAS, data_size, error, memspace, pard%dataspace)
        CALL h5sclose_f(pard%dataspace, error)
        CALL h5sclose_f(memspace, error)
        CALL h5dclose_f(pard%dataset, error)
      END IF
    END IF
  END SUBROUTINE HDF5_ADD_DATA2I
     
  ! 2D CHARACTER DATA
  SUBROUTINE HDF5_ADD_DATA2C(TABNO, DATAX, DATAY, DATAS, LENS)
    USE ISO_C_BINDING  
    IMPLICIT NONE    
    
    INTEGER                               :: TABNO, LENS
    INTEGER                               :: DATAX, DATAY
    CHARACTER(3), POINTER                 :: DATAS(:,:)
    CHARACTER(3), POINTER                 :: DATAS2(:,:)
    CHARACTER(3), TARGET                  :: DATAS00(DATAX,DATAY)
    TYPE(DSETPARS)                        :: pard
    TYPE(C_PTR)                           :: f_ptr
    
    DATAS => DATAS00
    DATAS2 => DATAS00

    CALL FIND_ID(TABNO, 'N', pard)
    CALL H5Tset_size_f(pard%datatype, LENS, error)
    
    CALL h5dopen_f(pard%group, pard%name, pard%dataset, error)   
    IF (error==-1) THEN
      pard%status = .FALSE.
      CALL h5screate_simple_f(pard%rank, pard%inidims, pard%dataspace, error, pard%maxdims)
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, pard%crplist, error)
      CALL h5pset_chunk_f(pard%crplist, pard%rank, pard%inidims, error)
      CALL h5pset_deflate_f(pard%crplist, 6, error) 
      CALL h5dcreate_f(pard%group, pard%name, pard%datatype, pard%dataspace, pard%dataset, error, pard%crplist) 
      f_ptr = C_LOC(DATAS(1,1)(1:1))
      CALL h5dwrite_f(pard%dataset, pard%datatype, f_ptr, error)
      CALL h5pclose_f(pard%crplist, error)
      CALL h5dclose_f(pard%dataset, error)
      CALL h5sclose_f(pard%dataspace, error)
      DATAS2(1:DATAX,1:DATAY) = DATAS(1:DATAX,1:DATAY)     
    ELSE
      CALL h5dget_space_f(pard%dataset, pard%dataspace, error)
      CALL h5sget_simple_extent_dims_f(pard%dataspace, pard%dims, pard%maxdims, error)
      IF (pard%status) THEN
        DATAS2(1:DATAX,1:DATAY) = DATAS
      ELSE
        DATAS2(1:DATAX,pard%dims(2)+1:pard%dims(2)+DATAY) = DATAS
        pard%dims(2) =  pard%dims(2) + DATAY 
        CALL h5dset_extent_f(pard%dataset, pard%dims, error)
      END IF
      f_ptr = C_LOC(DATAS2(1,1)(1:1))
      CALL h5dwrite_f(pard%dataset, pard%datatype, f_ptr, error)
      CALL h5sclose_f(pard%dataspace, error)
      CALL h5dclose_f(pard%dataset, error)
      pard%status = .FALSE.
    END IF
  END SUBROUTINE HDF5_ADD_DATA2C
     
  SUBROUTINE HDF5_ADD_DATA2N(TABNO, SETNAME, DATAX, DATAY, DATAS)
    IMPLICIT NONE  
    
    INTEGER                               :: TABNO
    CHARACTER(*)                          :: SETNAME
    INTEGER                               :: DATAX, DATAY
    REAL, DIMENSION(DATAX,DATAY)          :: DATAS
    TYPE(DSETPARS)                        :: pard
    INTEGER(HSIZE_T), DIMENSION(1:2)      :: offset, data_size   
    INTEGER(HSIZE_T)                      :: dim1, dim2
    INTEGER(HID_T)                        :: memspace
    
    CALL FIND_ID(TABNO, SETNAME, pard)
    data_size(1:2)  = (/DATAX, DATAY/)
    dim1 = 0
    
    CALL h5dopen_f(pard%group, pard%name, pard%dataset, error) 
    IF (error==-1) THEN
      pard%status = .TRUE.
      CALL h5screate_simple_f(pard%rank, pard%inidims, pard%dataspace, error, pard%maxdims)
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, pard%crplist, error)
      CALL h5pset_chunk_f(pard%crplist, pard%rank, pard%inidims, error)
      CALL h5pset_deflate_f(pard%crplist, 6, error) 
      CALL h5dcreate_f(pard%group, pard%name, pard%datatype, pard%dataspace, pard%dataset, error, pard%crplist)
      CALL h5dwrite_f(pard%dataset, pard%datatype, DATAS, data_size, error)
      CALL h5pclose_f(pard%crplist, error)
      CALL h5dclose_f(pard%dataset, error)
      CALL h5sclose_f(pard%dataspace, error) 
    ELSE
      IF (pard%status) THEN
        CALL h5dwrite_f(pard%dataset, pard%datatype, DATAS, data_size, error)
        CALL h5dclose_f(pard%dataset, error)
        pard%status = .FALSE.
      ELSE
        CALL h5dget_space_f(pard%dataset, pard%dataspace, error)
        CALL h5sget_simple_extent_dims_f(pard%dataspace, pard%dims, pard%maxdims, error)
        dim2 = pard%dims(2)
        offset(1:2) = (/dim1, dim2/)
        pard%dims(2) = pard%dims(2) + DATAY 
        CALL h5dset_extent_f(pard%dataset, pard%dims, error)
        CALL h5screate_simple_f(pard%rank, data_size, memspace, error) 
        CALL h5dget_space_f(pard%dataset, pard%dataspace, error)
        CALL h5sselect_hyperslab_f(pard%dataspace, H5S_SELECT_SET_F, offset, data_size, error)
        CALL h5dwrite_f(pard%dataset, pard%datatype, DATAS, data_size, error, memspace, pard%dataspace)
        CALL h5sclose_f(pard%dataspace, error)
        CALL h5sclose_f(memspace, error)
        CALL h5dclose_f(pard%dataset, error)
      END IF
    END IF
  END SUBROUTINE HDF5_ADD_DATA2N
  
  !=========================================================================================================================== 
  ! 3D DATA
  SUBROUTINE HDF5_ADD_DATA3(TABNO, DATAX, DATAY, DATAZ, DATAS, DIMX)
    IMPLICIT NONE       
    
    INTEGER                               :: TABNO
    INTEGER                               :: DATAX, DATAY, DATAZ
    REAL, DIMENSION(DATAX,DATAY, DATAZ)   :: DATAS
    CHARACTER                             :: DIMX
    TYPE(DSETPARS)                        :: pard
    INTEGER(HSIZE_T), DIMENSION(1:3)      :: offset, data_size
    INTEGER(HID_T)                        :: memspace
    INTEGER(HSIZE_T)                      :: dim1, dim2, dim3
    
    CALL FIND_ID(TABNO, 'N', pard)
    data_size(1:3)  = (/DATAX, DATAY, DATAZ/)
    
    CALL h5dopen_f(pard%group, pard%name, pard%dataset, error) 
    IF (error==-1) THEN 
      pard%status = .TRUE.
      CALL h5screate_simple_f(pard%rank, pard%inidimss, pard%dataspace, error, pard%maxdimss)
      CALL h5pcreate_f(H5P_DATASET_CREATE_F, pard%crplist, error)
      CALL h5pset_chunk_f(pard%crplist, pard%rank, pard%inidimss, error)
      CALL h5pset_deflate_f(pard%crplist, 6, error) 
      CALL h5dcreate_f(pard%group, pard%name, pard%datatype, pard%dataspace, pard%dataset, error, pard%crplist)
      CALL h5dwrite_f(pard%dataset, pard%datatype, DATAS, data_size, error)
      CALL h5pclose_f(pard%crplist, error)
      CALL h5dclose_f(pard%dataset, error)
      CALL h5sclose_f(pard%dataspace, error)    
    ELSE
      IF (pard%status) THEN
        CALL h5dwrite_f(pard%dataset, pard%datatype, DATAS, data_size, error)
        CALL h5dclose_f(pard%dataset, error)
        pard%status = .FALSE.
      ELSE
        CALL h5dget_space_f(pard%dataset, pard%dataspace, error)
        CALL h5sget_simple_extent_dims_f(pard%dataspace, pard%dimss, pard%maxdimss, error)
        IF (DIMX== 'X') THEN
        dim1 = pard%dimss(1)
        dim2 = 0
        dim3 = 0
        pard%dimss(1) =  pard%dimss(1) + DATAX
        pard%dimss(2) = pard%inidimss(2)
        pard%dimss(3) = pard%inidimss(3)
        END IF
        IF (DIMX== 'Y') THEN
          dim1 = 0
          dim2 = pard%dimss(2)
          dim3 = 0
          pard%dimss(2) =  pard%dimss(2) + DATAY
          pard%dimss(1) = pard%inidimss(1)
          pard%dimss(3) = pard%inidimss(3)
        END IF
        IF (DIMX== 'Z') THEN
          dim1 = 0
          dim2 = 0
          dim3 = pard%dimss(3)
          pard%dimss(3) =  pard%dimss(3) + DATAZ
          pard%dimss(1) = pard%inidimss(1)
          pard%dimss(2) = pard%inidimss(2)
        END IF
        offset(1:3) = (/dim1, dim2, dim3/)    
        CALL h5dset_extent_f(pard%dataset, pard%dimss, error)
        CALL h5screate_simple_f(pard%rank, data_size, memspace, error) 
        CALL h5dget_space_f(pard%dataset, pard%dataspace, error)
        CALL h5sselect_hyperslab_f(pard%dataspace, H5S_SELECT_SET_F, offset, data_size, error)
        CALL h5dwrite_f(pard%dataset, pard%datatype, DATAS, data_size, error, memspace, pard%dataspace)
        CALL h5sclose_f(pard%dataspace, error)
        CALL h5sclose_f(memspace, error)
        CALL h5dclose_f(pard%dataset, error)
      END IF
    END IF
  END SUBROUTINE HDF5_ADD_DATA3
   
  !=========================================================================================================================== 
  ! INTEGER ATTRIBUTE
  SUBROUTINE HDF5_ADD_ATTRIBUTEI(ITEMNO, ATTRINAME, ATTRIDATA)
    IMPLICIT NONE
    
    INTEGER        :: ITEMNO
    CHARACTER(*)   :: ATTRINAME
    INTEGER        :: ATTRIDATA
    INTEGER(HID_T) :: atype
    TYPE(DSETPARS) :: para
    
    CALL FIND_ID(ITEMNO, 'N', para)
    CALL h5dopen_f(para%group, para%name, para%dataset, error)
    CALL h5screate_simple_f(1, adims1, para%aspace, error)
    CALL h5screate_f(H5S_SCALAR_F, para%aspace, error)
    CALL h5acreate_f(para%dataset, ATTRINAME, H5T_NATIVE_INTEGER, para%aspace, para%aID, error)
    CALL h5awrite_f(para%aID, H5T_NATIVE_INTEGER, ATTRIDATA, adims1, error)
    CALL h5aclose_f(para%aID, error)
    CALL h5sclose_f(para%aspace, error)
    CALL h5dclose_f(para%dataset, error)      
  END SUBROUTINE HDF5_ADD_ATTRIBUTEI
     
  ! CHARACTER ATTRIBUTE
  SUBROUTINE HDF5_ADD_ATTRIBUTEC(ITEMNO, ATTRINAME, ATTRIDATA)
    IMPLICIT NONE
    
    INTEGER         :: ITEMNO
    CHARACTER(*)    :: ATTRINAME
    CHARACTER(*)    :: ATTRIDATA
    INTEGER(HID_T)  :: atype
    TYPE(DSETPARS)  :: para
    
    CALL FIND_ID(ITEMNO, 'N', para)
    CALL h5dopen_f(para%group, para%name, para%dataset, error)
    CALL h5screate_simple_f(1, adims1, para%aspace, error)
    CALL h5screate_f(H5S_SCALAR_F, para%aspace, error)
    CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype, error)
    CALL h5tset_size_f(atype, LEN_TRIM(ATTRIDATA), error)
    CALL h5acreate_f(para%dataset, ATTRINAME, atype, para%aspace, para%aID, error)
    CALL h5awrite_f(para%aID, atype, ATTRIDATA, adims1, error)
    CALL h5aclose_f(para%aID, error)
    CALL h5sclose_f(para%aspace, error)
    CALL h5tclose_f(atype, error)
    CALL h5dclose_f(para%dataset, error)
  END SUBROUTINE HDF5_ADD_ATTRIBUTEC
     
  ! REAL ATTRIBUTE
  SUBROUTINE HDF5_ADD_ATTRIBUTER(ITEMNO, ATTRINAME, PRECISIONS, ATTRIDATA)
    IMPLICIT NONE  
    
    INTEGER        :: ITEMNO
    CHARACTER(*)   :: ATTRINAME
    INTEGER        :: PRECISIONS
    REAL           :: ATTRIDATA
    REAL           :: RDATA
    INTEGER(HID_T) :: atype
    TYPE(DSETPARS) :: para
    
    CALL FIND_ID(ITEMNO, 'N', para)
    CALL XINT(ATTRIDATA, PRECISIONS, RDATA)
    CALL h5dopen_f(para%group, para%name, para%dataset, error)
    CALL h5screate_simple_f(1, adims1, para%aspace, error)
    CALL h5screate_f(H5S_SCALAR_F, para%aspace, error)
    CALL h5tcopy_f(H5T_NATIVE_DOUBLE, atype, error)
    CALL h5acreate_f(para%dataset, ATTRINAME, atype, para%aspace, para%aID, error)
    CALL h5awrite_f(para%aID, atype, RDATA, adims1, error)
    CALL h5aclose_f(para%aID, error)
    CALL h5sclose_f(para%aspace, error)
    CALL h5tclose_f(atype, error)
    CALL h5dclose_f(para%dataset, error)
  END SUBROUTINE HDF5_ADD_ATTRIBUTER
     
  ! 1D INTEGER ATTRIBUTE 
  SUBROUTINE HDF5_ADD_ATTRIBUTEIS(ITEMNO, ATTRINAME, ATTRISIZE, ATTRIDATA)
    IMPLICIT NONE     
    
    INTEGER                           :: ITEMNO
    CHARACTER(*)                      :: ATTRINAME
    INTEGER                           :: ATTRISIZE
    INTEGER, DIMENSION(ATTRISIZE)     :: ATTRIDATA
    INTEGER(HID_T)                    :: atype
    INTEGER(HSIZE_T), DIMENSION(1)    :: adims
    TYPE(DSETPARS)                    :: para
    
    CALL FIND_ID(ITEMNO, 'N', para) 
    adims = ATTRISIZE
    CALL h5dopen_f(para%group, para%name, para%dataset, error)
    CALL h5screate_simple_f(1, adims, para%aspace, error)
    CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype, error)
    CALL h5acreate_f(para%dataset, ATTRINAME, atype, para%aspace, para%aID, error)
    CALL h5awrite_f(para%aID, atype, ATTRIDATA, adims, error)
    CALL h5aclose_f(para%aID, error)
    CALL h5sclose_f(para%aspace, error)
    CALL h5tclose_f(atype, error)
    CALL h5dclose_f(para%dataset, error)
  END SUBROUTINE HDF5_ADD_ATTRIBUTEIS
     
  ! 1D CHARACTER ATTRIBUTE 
  SUBROUTINE HDF5_ADD_ATTRIBUTECS(ITEMNO, ATTRINAME, ATTRISIZE, ATTRIDATA)
    IMPLICIT NONE      
    
    INTEGER                               :: ITEMNO
    CHARACTER(*)                          :: ATTRINAME
    INTEGER                               :: ATTRISIZE
    CHARACTER(*),DIMENSION(ATTRISIZE)     :: ATTRIDATA
    INTEGER(HID_T)                        :: atype
    INTEGER(HSIZE_T), DIMENSION(1)        :: adims
    TYPE(DSETPARS)                        :: para
    
    CALL FIND_ID(ITEMNO, 'N', para)
    adims = ATTRISIZE
    CALL h5dopen_f(para%group, para%name, para%dataset, error)
    CALL h5screate_simple_f(1, adims, para%aspace, error)
    CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype, error)
    CALL h5tset_size_f(atype, 99, error)
    CALL h5acreate_f(para%dataset, ATTRINAME, atype, para%aspace, para%aID, error)
    CALL h5awrite_f(para%aID, atype, ATTRIDATA, adims, error)
    CALL h5aclose_f(para%aID, error)
    CALL h5sclose_f(para%aspace, error)
    CALL h5tclose_f(atype, error)
    CALL h5dclose_f(para%dataset, error)
  END SUBROUTINE HDF5_ADD_ATTRIBUTECS
     
  ! 1D CHARACTER ATTRIBUTE 
  SUBROUTINE HDF5_ADD_ATTRIBUTECS_LENGTH(ITEMNO, ATTRINAME, ATTRISIZE, ATTRIDATA, LENGTH)
    IMPLICIT NONE  
    
    INTEGER                               :: ITEMNO, LENGTH
    CHARACTER(*)                          :: ATTRINAME
    INTEGER                               :: ATTRISIZE
    CHARACTER(*),DIMENSION(ATTRISIZE)     :: ATTRIDATA
    INTEGER(HID_T)                        :: atype
    INTEGER(HSIZE_T), DIMENSION(1)        :: adims
    TYPE(DSETPARS)                        :: para
    
    CALL FIND_ID(ITEMNO, 'N', para)
    adims = ATTRISIZE
    CALL h5dopen_f(para%group, para%name, para%dataset, error)
    CALL h5screate_simple_f(1, adims, para%aspace, error)
    CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype, error)
    CALL h5tset_size_f(atype, LENGTH, error)
    CALL h5acreate_f(para%dataset, ATTRINAME, atype, para%aspace, para%aID, error)
    CALL h5awrite_f(para%aID, atype, ATTRIDATA, adims, error)
    CALL h5aclose_f(para%aID, error)
    CALL h5sclose_f(para%aspace, error)
    CALL h5tclose_f(atype, error)
    CALL h5dclose_f(para%dataset, error)
  END SUBROUTINE HDF5_ADD_ATTRIBUTECS_LENGTH
     
  ! 1D REAL ATTRIBUTE 
  SUBROUTINE HDF5_ADD_ATTRIBUTERS(ITEMNO, ATTRINAME, ATTRISIZE, PRECISIONS, ATTRIDATA)
    IMPLICIT NONE  
    
    INTEGER                           :: ITEMNO
    CHARACTER(*)                      :: ATTRINAME
    INTEGER                           :: PRECISIONS
    INTEGER                           :: ATTRISIZE
    REAL, DIMENSION(ATTRISIZE)        :: ATTRIDATA
    REAL, DIMENSION(ATTRISIZE)        :: RDATA
    INTEGER(HID_T)                    :: atype
    INTEGER(HSIZE_T), DIMENSION(1)    :: adims
    TYPE(DSETPARS)                    :: para
    
    CALL FIND_ID(ITEMNO, 'N', para)
    ATTRISIZE = SIZE(ATTRIDATA)
    DO I = 1, ATTRISIZE
      CALL XINT(ATTRIDATA(I), PRECISIONS, RDATA(I))
    END DO
    adims = ATTRISIZE
    CALL h5dopen_f(para%group, para%name, para%dataset, error)
    CALL h5screate_simple_f(1, adims, para%aspace, error)
    CALL h5tcopy_f(H5T_NATIVE_DOUBLE, atype, error)
    CALL h5acreate_f(para%dataset, ATTRINAME, atype, para%aspace, para%aID, error)
    CALL h5awrite_f(para%aID, atype, RDATA, adims, error)
    CALL h5aclose_f(para%aID, error)
    CALL h5sclose_f(para%aspace, error)
    CALL h5tclose_f(atype, error)
    CALL h5dclose_f(para%dataset, error)
  END SUBROUTINE HDF5_ADD_ATTRIBUTERS
     
  ! INTEGER ATTRIBUTE TO TIME (JDAY)
  SUBROUTINE HDF5_ADD_ATTRIBUTET1(NO, ITEMNO, ATTRINAME)
    IMPLICIT NONE  
    
    INTEGER           :: NO, ITEMNO
    CHARACTER(*)      :: ATTRINAME
    INTEGER(HSIZE_T)  :: DIM1
    INTEGER           :: DIM0
    TYPE(DSETPARS)    :: para
    
    IF (NO==1) THEN
      CALL FIND_ID(ITEMNO, 'N', para)
      CALL h5dopen_f(para%time_group, 'Time', para%timeset(1), error)
      CALL h5aopen_f(para%timeset(1), ATTRINAME, para%aID, error)
      IF (error==-1) THEN
        CALL h5screate_simple_f(1, adims1, para%aspace, error)
        CALL h5screate_f(H5S_SCALAR_F, para%aspace, error)
        CALL h5acreate_f(para%timeset(1), ATTRINAME, H5T_NATIVE_INTEGER, para%aspace, para%aID, error)
        CALL h5awrite_f(para%aID, H5T_NATIVE_INTEGER, 1, adims1, error)
        CALL h5aclose_f(para%aID, error)
        CALL h5sclose_f(para%aspace, error)
        CALL h5dclose_f(para%timeset(1), error)
      ELSE
        CALL h5dget_space_f(para%timeset(1), para%timespace(1), error)
        CALL h5sget_simple_extent_npoints_f(para%timespace(1), DIM1, error)
        DIM0 = DIM1
        CALL h5awrite_f(para%aID, H5T_NATIVE_INTEGER, DIM0, adims1, error)
        CALL h5aclose_f(para%aID, error)
        CALL h5sclose_f(para%aspace, error)
        CALL h5sclose_f(para%timespace(1), error)
        CALL h5dclose_f(para%timeset(1), error)
      END IF
    END IF
  END SUBROUTINE HDF5_ADD_ATTRIBUTET1
     
  ! INTEGER ATTRIBUTE TO TIME (JDAY)
  SUBROUTINE HDF5_ADD_ATTRIBUTET_JW(NO, JW, ITEMNO, ATTRINAME)
    IMPLICIT NONE  
    
    INTEGER           :: NO, JW, ITEMNO
    CHARACTER(*)      :: ATTRINAME
    INTEGER(HSIZE_T)  :: DIM1
    INTEGER           :: DIM0
    TYPE(DSETPARS)    :: para
    
    IF (NO==1) THEN
      CALL FIND_ID(ITEMNO, 'N', para)
      CALL h5dopen_f(para%time_group, 'Time'//HDF5_JWBS(JW), para%timeset(1), error)
      CALL h5aopen_f(para%timeset(1), ATTRINAME, para%aID, error)
      IF (error==-1) THEN
        CALL h5screate_simple_f(1, adims1, para%aspace, error)
        CALL h5screate_f(H5S_SCALAR_F, para%aspace, error)
        CALL h5acreate_f(para%timeset(1), ATTRINAME, H5T_NATIVE_INTEGER, para%aspace, para%aID, error)
        CALL h5awrite_f(para%aID, H5T_NATIVE_INTEGER, 1, adims1, error)
        CALL h5aclose_f(para%aID, error)
        CALL h5sclose_f(para%aspace, error)
        CALL h5dclose_f(para%timeset(1), error)
      ELSE
        CALL h5dget_space_f(para%timeset(1), para%timespace(1), error)
        CALL h5sget_simple_extent_npoints_f(para%timespace(1), DIM1, error)
        DIM0 = DIM1
        CALL h5awrite_f(para%aID, H5T_NATIVE_INTEGER, DIM0, adims1, error)
        CALL h5aclose_f(para%aID, error)
        CALL h5sclose_f(para%aspace, error)
        CALL h5sclose_f(para%timespace(1), error)
        CALL h5dclose_f(para%timeset(1), error)
      END IF
    END IF
  END SUBROUTINE HDF5_ADD_ATTRIBUTET_JW
     
  ! STRING ATTRIBUTE 
  SUBROUTINE HDF5_ADD_ATTRIBUTEC1(NO, ITEMNO, ATTRINAME, ATTRIDATA)
    IMPLICIT NONE 
    
    INTEGER           :: NO, ITEMNO
    CHARACTER(*)      :: ATTRINAME
    CHARACTER(*)      :: ATTRIDATA
    INTEGER(HID_T)    :: atype
    TYPE(DSETPARS)    :: para
       
    IF (NO==1) THEN
      CALL FIND_ID(ITEMNO, 'N', para)
      CALL h5dopen_f(para%time_group, 'Time', para%timeset(1), error)
      CALL h5screate_simple_f(1, adims1, para%aspace, error)
      CALL h5screate_f(H5S_SCALAR_F, para%aspace, error)
      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype, error)
      CALL h5tset_size_f(atype, LEN_TRIM(ATTRIDATA), error)
      CALL h5acreate_f(para%timeset(1), ATTRINAME, atype, para%aspace, para%aID, error)
      CALL h5awrite_f(para%aID, atype, ATTRIDATA, adims1, error)
      CALL h5aclose_f(para%aID, error)
      CALL h5sclose_f(para%aspace, error)
      CALL h5tclose_f(atype, error)
      CALL h5dclose_f(para%timeset(1), error)
    END IF
    
    IF (NO==0) THEN
      CALL FIND_ID(ITEMNO, 'Name', para)
      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype, error)
      CALL h5tset_size_f(atype, LEN_TRIM(ATTRIDATA), error)
      CALL h5screate_simple_f(1, adims1, para%aspace, error)
      CALL h5acreate_f(para%group, ATTRINAME, atype, para%aspace, para%aID, error)
      CALL h5awrite_f(para%aID, atype, ATTRIDATA, adims1, error)
      CALL h5aclose_f(para%aID, error)
      CALL h5sclose_f(para%aspace, error)
      CALL h5tclose_f(atype, error)
    END IF
  END SUBROUTINE HDF5_ADD_ATTRIBUTEC1
     
  ! STRING ATTRIBUTE 
  SUBROUTINE HDF5_ADD_ATTRIBUTEC_JW(NO, JW, ITEMNO, ATTRINAME, ATTRIDATA)
    IMPLICIT NONE     
    
    INTEGER           :: NO, JW, ITEMNO
    CHARACTER(*)      :: ATTRINAME
    CHARACTER(*)      :: ATTRIDATA
    INTEGER(HID_T)    :: atype
    TYPE(DSETPARS)    :: para
       
    IF (NO==1) THEN
      CALL FIND_ID(ITEMNO, 'N', para)
      CALL h5dopen_f(para%group, 'Time'//HDF5_JWBS(JW), para%timeset(1), error)
      CALL h5screate_simple_f(1, adims1, para%aspace, error)
      CALL h5screate_f(H5S_SCALAR_F, para%aspace, error)
      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype, error)
      CALL h5tset_size_f(atype, LEN_TRIM(ATTRIDATA), error)
      CALL h5acreate_f(para%timeset(1), ATTRINAME, atype, para%aspace, para%aID, error)
      CALL h5awrite_f(para%aID, atype, ATTRIDATA, adims1, error)
      CALL h5aclose_f(para%aID, error)
      CALL h5sclose_f(para%aspace, error)
      CALL h5tclose_f(atype, error)
      CALL h5dclose_f(para%timeset(1), error)
    END IF
    
    IF (NO==0) THEN
      CALL FIND_ID(ITEMNO, 'Name', para)
      CALL h5tcopy_f(H5T_NATIVE_CHARACTER, atype, error)
      CALL h5tset_size_f(atype, LEN_TRIM(ATTRIDATA), error)
      CALL h5screate_simple_f(1, adims1, para%aspace, error)
      CALL h5acreate_f(para%group, ATTRINAME, atype, para%aspace, para%aID, error)
      CALL h5awrite_f(para%aID, atype, ATTRIDATA, adims1, error)
      CALL h5aclose_f(para%aID, error)
      CALL h5sclose_f(para%aspace, error)
      CALL h5tclose_f(atype, error)
    END IF
  END SUBROUTINE HDF5_ADD_ATTRIBUTEC_JW
     
  ! INTEGER ATTRIBUTE 
  SUBROUTINE HDF5_ADD_ATTRIBUTEI1(NO, ITEMNO, ATTRINAME, ATTRIDATA)
    IMPLICIT NONE   
    
    INTEGER           :: NO, ITEMNO
    CHARACTER(*)      :: ATTRINAME
    INTEGER           :: ATTRIDATA
    INTEGER(HID_T)    :: atype
    TYPE(DSETPARS)    :: para
       
    IF (NO==0) THEN
      CALL FIND_ID(ITEMNO, 'Name', para)
      CALL h5screate_simple_f(1, adims1, para%aspace, error)
      CALL h5screate_f(H5S_SCALAR_F, para%aspace, error)
      CALL h5acreate_f(para%group, ATTRINAME, H5T_NATIVE_INTEGER, para%aspace, para%aID, error)
      CALL h5awrite_f(para%aID, H5T_NATIVE_INTEGER, ATTRIDATA, adims1, error)
      CALL h5aclose_f(para%aID, error)
      CALL h5sclose_f(para%aspace, error)
    END IF
  END SUBROUTINE HDF5_ADD_ATTRIBUTEI1
     
  ! INTEGER ATTRIBUTE 
  SUBROUTINE HDF5_ADD_ATTRIBUTEI_JW(NO, JW, ITEMNO, ATTRINAME, ATTRIDATA)
    IMPLICIT NONE
       
    INTEGER           :: NO, JW, ITEMNO
    CHARACTER(*)      :: ATTRINAME
    INTEGER           :: ATTRIDATA
    INTEGER(HID_T)    :: atype
    TYPE(DSETPARS)    :: para
       
    IF (NO==1) THEN
      CALL FIND_ID(ITEMNO, 'N', para)
      CALL h5dopen_f(para%group, 'Time'//HDF5_JWBS(JW), para%timeset(1), error)
      CALL h5screate_simple_f(1, adims1, para%aspace, error)
      CALL h5screate_f(H5S_SCALAR_F, para%aspace, error)
      CALL h5tcopy_f(H5T_NATIVE_INTEGER, atype, error)
      CALL h5acreate_f(para%timeset(1), ATTRINAME, atype, para%aspace, para%aID, error)
      CALL h5awrite_f(para%aID, atype, ATTRIDATA, adims1, error)
      CALL h5aclose_f(para%aID, error)
      CALL h5sclose_f(para%aspace, error)
      CALL h5tclose_f(atype, error)
      CALL h5dclose_f(para%timeset(1), error)
    END IF
    
    IF (NO==0) THEN
      CALL FIND_ID(ITEMNO, 'Name', para)
      CALL h5screate_simple_f(1, adims1, para%aspace, error)
      CALL h5screate_f(H5S_SCALAR_F, para%aspace, error)
      CALL h5acreate_f(para%group, ATTRINAME, H5T_NATIVE_INTEGER, para%aspace, para%aID, error)
      CALL h5awrite_f(para%aID, H5T_NATIVE_INTEGER, ATTRIDATA, adims1, error)
      CALL h5aclose_f(para%aID, error)
      CALL h5sclose_f(para%aspace, error)
    END IF
  END SUBROUTINE HDF5_ADD_ATTRIBUTEI_JW
       
  SUBROUTINE HDF5_ADD_ATTRIBUTER1(NO, ITEMNO, ATTRINAME, PRECISIONS, ATTRIDATA)
    IMPLICIT NONE  
    
    INTEGER           :: NO, ITEMNO
    CHARACTER(*)      :: ATTRINAME
    INTEGER           :: PRECISIONS
    REAL              :: ATTRIDATA
    REAL              :: RDATA
    INTEGER(HID_T)    :: atype
    TYPE(DSETPARS)    :: para
       
    IF (NO==0) THEN
      CALL FIND_ID(ITEMNO, 'N', para)
      CALL XINT(ATTRIDATA, PRECISIONS, RDATA)
      CALL h5screate_simple_f(1, adims1, para%aspace, error)
      CALL h5screate_f(H5S_SCALAR_F, para%aspace, error)
      CALL h5tcopy_f(H5T_NATIVE_DOUBLE, atype, error)
      CALL h5acreate_f(para%group, ATTRINAME, atype, para%aspace, para%aID, error)
      CALL h5awrite_f(para%aID, atype, RDATA, adims1, error)
      CALL h5aclose_f(para%aID, error)
      CALL h5sclose_f(para%aspace, error)
      CALL h5tclose_f(atype, error)
    END IF
  END SUBROUTINE HDF5_ADD_ATTRIBUTER1
          
  !****************************************************************************************************************************************!    
  ! TIME STAMP
  SUBROUTINE JDAYTOTS(JDAY, STAMP)
    IMPLICIT NONE  
    
    REAL         :: JDAY
    CHARACTER(19):: STAMP
    CHARACTER(2) :: GDAYC, HOURC, MINC, SECC
    CHARACTER(3) :: MSECC, MONTHC
    CHARACTER(4) :: YEARC
    
    CALL JDAYTODATE(JDAY, GDAYC, MONTHC, YEARC, HOURC, MINC, SECC)
    STAMP = GDAYC//MONTHC//' '//YEARC//' '//HOURC//':'//MINC//':'//SECC
  END SUBROUTINE JDAYTOTS 
          
  SUBROUTINE JDAYTODATE(JDAY, GDAYCC, MONTHC, YEARC, HOURC, MINC, SECC)
    USE GDAYC
    IMPLICIT NONE 
    
    REAL         :: JDAY
    CHARACTER(2) :: GDAYCC, HOURC, MINC, SECC
    CHARACTER(3) :: MONTHC
    CHARACTER(4) :: YEARC
    INTEGER      :: GDAY1, YEAR1, JDAYG1, INCR1
    CHARACTER(9) :: MONTH1
    LOGICAL      :: LEAP_YEAR1
             
    WRITE(YEARC, '(I4)') YEAR
    MONTHC = ADJUSTL(MONTH)
    WRITE(GDAYCC, '(I2)') GDAY
    WRITE(HOURC, '(I2)') INT((JDAY-INT(JDAY))*24)
    IF(HOURC(1:1)==' ') HOURC(1:1) = '0'
    WRITE(MINC,  '(I2)') INT(((JDAY-INT(JDAY))*24-INT((JDAY-INT(JDAY))*24))*60)
    IF(MINC(1:1)==' ') MINC(1:1) = '0'
    WRITE(SECC,  '(I2)') INT((((JDAY-INT(JDAY))*24-INT((JDAY-INT(JDAY))*24))*60-INT(((JDAY-INT(JDAY))*24-INT((JDAY-INT(JDAY))*24))*60))*60)
    IF(SECC(1:1)==' ') SECC(1:1) = '0'
  END SUBROUTINE JDAYTODATE
     
  !****************************************************************************************************************************************!     
  ! BALANCES
  SUBROUTINE HDF5_BALANCES1(JW, VOLINJW, VOLPRJW, VOLOUTJW, VOLWDJW, VOLEVJW, VOLDTJW, VOLTRBJW, VOLICEJW, FLOWBFN, JDAY)
    USE HDF5MOD; USE MAIN, ONLY: VOLUME_BALANCE, DLVR    
    IMPLICIT NONE   
    
    INTEGER       :: JW, FLOWBFN
    REAL          :: VOLINJW, VOLPRJW, VOLOUTJW, VOLWDJW, VOLEVJW, VOLDTJW, VOLTRBJW, VOLICEJW, JDAY
        
    HDF5_JWB = JW
    ALLOCATENO = 9
    IF(VOLUME_BALANCE(JW)) ALLOCATENO = 10
    ALLOCATE(DATASET2(ALLOCATENO,1))
    CALL XINT(JW, 0, DATASET2(1,1))
    CALL XINT(VOLINJW, 3, DATASET2(2,1))
    CALL XINT(VOLPRJW, 3, DATASET2(3,1))
    CALL XINT(VOLOUTJW, 3, DATASET2(4,1))
    CALL XINT(VOLWDJW, 3, DATASET2(5,1))
    CALL XINT(VOLEVJW, 3, DATASET2(6,1))
    CALL XINT(VOLDTJW, 3, DATASET2(7,1))
    CALL XINT(VOLTRBJW, 3, DATASET2(8,1))
    CALL XINT(VOLICEJW, 3, DATASET2(9,1))
    IF(VOLUME_BALANCE(JW)) DATASET2(10,1) = DLVR(JW)
    CALL HDF5_ADD_DATA(FLOWBFN, ALLOCATENO, 1, DATASET2)
    CALL HDF5_ADD_DATE(FLOWBFN, JDAY)
    DEALLOCATE(DATASET2)               
  END SUBROUTINE HDF5_BALANCES1
     
  SUBROUTINE HDF5_BALANCES2(JW, TPWB, TPSED, TPPLANT, TNWB, TNSED, TNPLANT, MASSBFN)
    USE HDF5MOD; USE MAIN, ONLY: VOLUME_BALANCE 
    USE GLOBAL,   ONLY: TPOUT, TPTRIB, TPDTRIB, TPWD, TPPR, TPIN, TP_SEDSOD_PO4, TNOUT,  &
                        TNTRIB, TNDTRIB, TNWD, TNPR, TNIN, TN_SEDSOD_NH4, ATMDEP_P, ATMDEP_N, NH3GASLOSS
    USE KINETIC,  ONLY: PFLUXIN, NFLUXIN
    USE CEMAVars, ONLY: SEDIMENT_DIAGENESIS, SDPFLUX, SDNH4FLUX, SDNO3FLUX
    IMPLICIT NONE     
    
    INTEGER   :: JW, MASSBFN
    REAL      :: TPWB, TPSED, TPPLANT, TNWB, TNSED, TNPLANT 
     
    IF (SEDIMENT_DIAGENESIS) THEN
      ALLOCATENO = 29 
      ALLOCATE(DATASET2(ALLOCATENO,1))
    ELSE
      ALLOCATENO = 26
      ALLOCATE(DATASET2(ALLOCATENO,1))
    END IF
    CALL XINT(JW, 0, DATASET2(1,1)) 
    CALL XINT(TPWB, 3, DATASET2(2,1)) 
    CALL XINT(TPSED, 3, DATASET2(3,1)) 
    CALL XINT(TPPLANT, 3, DATASET2(4,1)) 
    CALL XINT(TPOUT(JW), 3, DATASET2(5,1)) 
    CALL XINT(TPTRIB(JW), 3, DATASET2(6,1)) 
    CALL XINT(TPDTRIB(JW), 3, DATASET2(7,1)) 
    CALL XINT(TPWD(JW), 3, DATASET2(8,1)) 
    CALL XINT(TPPR(JW), 3, DATASET2(9,1)) 
    CALL XINT(TPIN(JW), 3, DATASET2(10,1)) 
    CALL XINT(ATMDEP_P(JW), 3, DATASET2(11,1))  
    CALL XINT(TP_SEDSOD_PO4(JW), 3, DATASET2(12,1)) 
    CALL XINT(PFLUXIN(JW), 3, DATASET2(13,1)) 
    
    IF (SEDIMENT_DIAGENESIS) THEN
      CALL XINT(SDPFLUX(JW), 3, DATASET2(14,1)) 
      CALL XINT(TNWB, 3, DATASET2(15,1)) 
      CALL XINT(TNSED, 3, DATASET2(16,1)) 
      CALL XINT(TNPLANT, 3, DATASET2(17,1))
      CALL XINT(TNOUT(JW), 3, DATASET2(18,1)) 
      CALL XINT(TNTRIB(JW), 3, DATASET2(19,1)) 
      CALL XINT(TNDTRIB(JW), 3, DATASET2(20,1))
      CALL XINT(TNWD(JW), 3, DATASET2(21,1)) 
      CALL XINT(TNPR(JW), 3, DATASET2(22,1)) 
      CALL XINT(TNIN(JW), 3, DATASET2(23,1)) 
      CALL XINT(ATMDEP_N(JW), 3, DATASET2(24,1))   
      CALL XINT(NH3GASLOSS(JW), 3, DATASET2(25,1))   
      CALL XINT(TN_SEDSOD_NH4(JW), 3, DATASET2(26,1))
      CALL XINT(NFLUXIN(JW), 3, DATASET2(27,1)) 
      CALL XINT(SDNH4FLUX(JW), 3, DATASET2(28,1)) 
      CALL XINT(SDNO3FLUX(JW), 3, DATASET2(29,1))
    ELSE
      CALL XINT(TNWB, 3, DATASET2(14,1)) 
      CALL XINT(TNSED, 3, DATASET2(15,1)) 
      CALL XINT(TNPLANT, 3, DATASET2(16,1)) 
      CALL XINT(TNOUT(JW), 3, DATASET2(17,1)) 
      CALL XINT(TNTRIB(JW), 3, DATASET2(18,1))
      CALL XINT(TNDTRIB(JW), 3, DATASET2(19,1)) 
      CALL XINT(TNWD(JW), 3, DATASET2(20,1))
      CALL XINT(TNPR(JW), 3, DATASET2(21,1)) 
      CALL XINT(TNIN(JW), 3, DATASET2(22,1)) 
      CALL XINT(ATMDEP_N(JW), 3, DATASET2(23,1))   
      CALL XINT(NH3GASLOSS(JW), 3, DATASET2(24,1))   
      CALL XINT(TN_SEDSOD_NH4(JW), 3, DATASET2(25,1)) 
      CALL XINT(NFLUXIN(JW), 3, DATASET2(26,1)) 
    END IF
    CALL HDF5_ADD_DATA(MASSBFN, ALLOCATENO, 1, DATASET2)
    DEALLOCATE(DATASET2)      
  END SUBROUTINE HDF5_BALANCES2
            
  !=========================================================================================================================== 
  ! OUTPUTINITW2TOOLS
  SUBROUTINE HDF5_OUTPUTINIT1
    USE HDF5MOD
    USE MAIN, ONLY: HDF, TSRFN1, ITSR, ICE_COMPUTATION, SEDIMENT_CALC, KFNAME2, CDN, NISNP, ISNP
    USE TVDC, ONLY: NAC, NACD, CN, CONSTITUENTS; USE KINETIC, ONLY: NAF, KFCN
    USE NAMESC, ONLY: CNAME2, CDNAME2, CUNIT, CUNIT3
    IMPLICIT NONE
           
    DO JW = 1, NWB
      HDF5_JWB = JW
      IF (HDF5_JTS(JW)>0) THEN
        HDF5_ITS = 0
        WRITE(CHARFROMI, '(I2)') JW
        CALL HDF5_OPEN_DATASET(TSRFN0+JW, TSRFN1(1:L1-1)//'_wb'//CHARFROMI, 'U',3)  
        ALLOCATE(STRING1(HDF5_A(JW)),INTS(HDF5_JTS(JW)))
        
        DO J = 1,NIKTSR
          IF (ITSR(J)>=US(BS(JW)) .AND. ITSR(J)<=DS(BE(JW))) THEN
            HDF5_ITS = HDF5_ITS +1
            INTS(HDF5_ITS) = ITSR(J)
          END IF
        END DO
        STRING1(1) = 'DLT(S)'
        STRING1(2) = 'ELWS(m)'
        STRING1(3) = 'T2(C)'
        STRING1(4) = 'U(m/s)'
        STRING1(5) = 'Q(m3/s)'
        STRING1(6) = 'SRON(W/m2)'
        STRING1(7) = 'EXT(m-1)'
        STRING1(8) = 'DEPTH(m)'
        STRING1(9) = 'WIDTH(m)'
        STRING1(10) = 'SHADE'
        
        IF (ICE_COMPUTATION) THEN
          STRING1(11) = 'ICETH(m)'
          STRING1(12) = 'Tvolag(C)'
          STRING1(13) = 'NetRad(W/m2)'
          STRING1(14) = 'SWSolar(W/m2)'
          STRING1(15) = 'LWRad(W/m2)'
          STRING1(16) = 'BackRad(W/m2)'
          STRING1(17) = 'EvapF(W/m2)'
          STRING1(18) = 'ConducF(W/m2)'
          STRING1(19) = 'ReaerationCoeff(day-1)'  
          
        IF (CONSTITUENTS) THEN
          DO JC = 1, NAC
            STRING1(19+JC) = CNAME2(CN(JC))//' ('//CUNIT(CN(JC))//")"
          END DO
          DO JE = 1, NEP
            WRITE(CHARFROMI, '(I4)')JE
            STRING1(19+NAC+JE) = 'EPI '//CHARFROMI
          END DO
          DO JM = 1, NMC
            WRITE(CHARFROMI, '(I4)')JM
            STRING1(19+NAC+NEP+JM) = 'MAC '//CHARFROMI
          END DO
          !
          IF (SEDIMENT_CALC(JW)) THEN 
            STRING1(19+NAC+NEP+NMC+1) = 'SED'//' (g/m3)'
            STRING1(19+NAC+NEP+NMC+2) = 'SEDP'//' (g/m3)'
            STRING1(19+NAC+NEP+NMC+3) = 'SEDN'//' (g/m3)'
            STRING1(19+NAC+NEP+NMC+4) = 'SEDC'//' (g/m3)'
            DO JD = 1, NACD(JW)
              STRING1(19+NAC+NEP+NMC+4+JD) = CDNAME2(CDN(JD,JW))//' ('//CUNIT3(CDN(JD,JW))//")"
            END DO
            DO JF = 1, NAF(JW)
              STRING1(19+NAC+NEP+NMC+4+NACD(JW)+JF) = KFNAME2(KFCN(JF,JW))
            END DO
            !
            DO JA = 1, NAL
              WRITE(CHARFROMI, '(I4)')JA
              STRING1(19+NAC+NEP+NMC+4+NACD(JW)+NAF(JW)+JA) = 'APLIM'//CHARFROMI
            END DO
            DO JA = 1, NAL
              WRITE(CHARFROMI, '(I4)')JA
              STRING1(19+NAC+NEP+NMC+4+NACD(JW)+NAF(JW)+NAL+JA) = 'ANLIM'//CHARFROMI
            END DO
            DO JA = 1, NAL
              WRITE(CHARFROMI, '(I4)')JA
              STRING1(19+NAC+NEP+NMC+4+NACD(JW)+NAF(JW)+2*NAL+JA) = 'ALLIM'//CHARFROMI
            END DO
            
          ELSE
            DO JD = 1, NACD(JW)
              STRING1(19+NAC+NEP+NMC+JD) = CDNAME2(CDN(JD,JW))//' ('//CUNIT3(CDN(JD,JW))//")"
            END DO
            DO JF = 1, NAF(JW)
              STRING1(19+NAC+NEP+NMC+NACD(JW)+JF) = KFNAME2(KFCN(JF,JW))
            END DO
            !
            DO JA = 1, NAL
              WRITE(CHARFROMI, '(I4)')JA
              STRING1(19+NAC+NEP+NMC+NACD(JW)+NAF(JW)+JA) = 'APLIM'//CHARFROMI
            END DO
            DO JA = 1, NAL
              WRITE(CHARFROMI, '(I4)')JA
              STRING1(19+NAC+NEP+NMC+NACD(JW)+NAF(JW)+NAL+JA) = 'ANLIM'//CHARFROMI
            END DO
            DO JA = 1, NAL
              WRITE(CHARFROMI, '(I4)')JA
              STRING1(19+NAC+NEP+NMC+NACD(JW)+NAF(JW)+2*NAL+JA) = 'ALLIM'//CHARFROMI
            END DO
          END IF
        END IF
          
        ELSE
          STRING1(11) = 'Tvolag(C)'
          STRING1(12) = 'NetRad(W/m2)'
          STRING1(13) = 'SWSolar(W/m2)'
          STRING1(14) = 'LWRad(W/m2)'
          STRING1(15) = 'BackRad(W/m2)'
          STRING1(16) = 'EvapF(W/m2)'
          STRING1(17) = 'ConducF(W/m2)'
          STRING1(18) = 'ReaerationCoeff(day-1)' 
          
        IF (CONSTITUENTS) THEN
          DO JC = 1, NAC
            STRING1(18+JC) = CNAME2(CN(JC))//' ('//CUNIT(CN(JC))//")"
          END DO
          DO JE = 1, NEP
            WRITE(CHARFROMI, '(I4)')JE
            STRING1(18+NAC+JE) = 'EPI '//CHARFROMI
          END DO
          DO JM = 1, NMC
            WRITE(CHARFROMI, '(I4)')JM
            STRING1(18+NAC+NEP+JM) = 'MAC '//CHARFROMI
          END DO
          !
          IF (SEDIMENT_CALC(JW)) THEN
            STRING1(18+NAC+NEP+NMC+1) = 'SED'//' (g/m3)'
            STRING1(18+NAC+NEP+NMC+2) = 'SEDP'//' (g/m3)'
            STRING1(18+NAC+NEP+NMC+3) = 'SEDN'//' (g/m3)'
            STRING1(18+NAC+NEP+NMC+4) = 'SEDC'//' (g/m3)'
            DO JD = 1, NACD(JW)
              STRING1(18+NAC+NEP+NMC+4+JD) = CDNAME2(CDN(JD,JW))//' ('//CUNIT3(CDN(JD,JW))//")"
            END DO
            DO JF = 1, NAF(JW)
              STRING1(18+NAC+NEP+NMC+4+NACD(JW)+JF) = KFNAME2(KFCN(JF,JW))
            END DO
            !
            DO JA = 1, NAL
              WRITE(CHARFROMI, '(I4)')JA
              STRING1(18+NAC+NEP+NMC+4+NACD(JW)+NAF(JW)+JA) = 'APLIM'//CHARFROMI
            END DO
            DO JA = 1, NAL
              WRITE(CHARFROMI, '(I4)')JA
              STRING1(18+NAC+NEP+NMC+4+NACD(JW)+NAF(JW)+NAL+JA) = 'ANLIM'//CHARFROMI
            END DO
            DO JA = 1, NAL
              WRITE(CHARFROMI, '(I4)')JA
              STRING1(18+NAC+NEP+NMC+4+NACD(JW)+NAF(JW)+2*NAL+JA) = 'ALLIM'//CHARFROMI
            END DO
            
          ELSE
            DO JD = 1, NACD(JW)
              STRING1(18+NAC+NEP+NMC+JD) = CDNAME2(CDN(JD,JW))//' ('//CUNIT3(CDN(JD,JW))//")"
            END DO
            DO JF = 1, NAF(JW)
              STRING1(18+NAC+NEP+NMC+NACD(JW)+JF) = KFNAME2(KFCN(JF,JW))
            END DO
            !
            DO JA = 1, NAL
              WRITE(CHARFROMI, '(I4)')JA
              STRING1(18+NAC+NEP+NMC+NACD(JW)+NAF(JW)+JA) = 'APLIM'//CHARFROMI
            END DO
            DO JA = 1, NAL
              WRITE(CHARFROMI, '(I4)')JA
              STRING1(18+NAC+NEP+NMC+NACD(JW)+NAF(JW)+NAL+JA) = 'ANLIM'//CHARFROMI
            END DO
            DO JA = 1, NAL
              WRITE(CHARFROMI, '(I4)')JA
              STRING1(18+NAC+NEP+NMC+NACD(JW)+NAF(JW)+2*NAL+JA) = 'ALLIM'//CHARFROMI
            END DO
          END IF
        END IF
        END IF
        
        CALL HDF5_ADD_ATTRIBUTE(TSRFN0+JW,'Time Series Dimension 2: Variable Names', HDF5_A(JW), STRING1)
        CALL HDF5_ADD_ATTRIBUTE(TSRFN0+JW, 'Time Series Dimension 3: Segment No', HDF5_JTS(JW), INTS)      
        CALL HDF5_ADD_ATTRIBUTE(TSRFN0+JW, 'Time Series Dimension 1:JDAY', 'As Time Stamp')
        DEALLOCATE(STRING1,INTS)
      END IF
    END DO    
  END SUBROUTINE HDF5_OUTPUTINIT1
  
  !=========================================================================================================================== 
  ! WDO
  SUBROUTINE HDF5_OUTPUTINIT4(JWD)
    USE HDF5MOD; USE MAIN, ONLY: WDO, IWDO
    IMPLICIT NONE
        
    INTEGER  :: JWD
        
    HDF5_J = JWD    
    CALL HDF5_OPEN_DATASET(WDO(JWD,1), 'qwo_'//SEGNUM(1:L), 'U')
    CALL HDF5_OPEN_DATASET(WDO(JWD,2), 'two_'//SEGNUM(1:L), 'U')
  
    ALLOCATE(STRING1(1+HDF5_NST))
    
    STRING1(1) = 'The sum of flows of Q (m3/s)'
    DO I = 1, HDF5_NST  
      WRITE(CHARFROMI, '(I4)') I
      STRING1(I+1) = 'Q '// CHARFROMI // ' (m3/s)'
    END DO
    CALL HDF5_ADD_ATTRIBUTE(WDO(JWD,1), 'Flow file for segment ', IWDO(JWD))
    CALL HDF5_ADD_ATTRIBUTE(WDO(JWD,1), 'Output Variable names', 1+HDF5_NST,STRING1)
    
    STRING1(1) = 'The sum of temperatures T (C)'
    DO I = 1, HDF5_NST 
      WRITE(CHARFROMI, '(I4)') I
      STRING1(I+1) = 'Temperature '// CHARFROMI //' T (C)'
    END DO
      
    CALL HDF5_ADD_ATTRIBUTE(WDO(JWD,2), '$Temperature file for segment ', IWDO(JWD))
    CALL HDF5_ADD_ATTRIBUTE(WDO(JWD,2), 'Output Variable names', 1+HDF5_NST,STRING1)
    DEALLOCATE(STRING1)
  END SUBROUTINE HDF5_OUTPUTINIT4
     
  SUBROUTINE HDF5_OUTPUTINIT5(JWD)
    USE HDF5MOD
    USE MAIN, ONLY: WDO, IWDO; USE TVDC, ONLY: NAC, CN; USE NAMESC, ONLY: CNAME2 
    IMPLICIT NONE
    
    INTEGER  :: JWD
        
    CALL HDF5_OPEN_DATASET(WDO(JWD,3), 'cwo_'//SEGNUM(1:L), 'U')
    ALLOCATE(STRING1(NAC)) 
    DO I = 1, NAC
      STRING1(I) = CNAME2(CN(I))//' (g/m3)'
    END DO
    CALL HDF5_ADD_ATTRIBUTE(WDO(JWD,3), '$Concentration file for segment ', IWDO(JWD))
    CALL HDF5_ADD_ATTRIBUTE(WDO(JWD,3), 'Output Variable names', NAC,STRING1)
    DEALLOCATE(STRING1)  
  END SUBROUTINE HDF5_OUTPUTINIT5
     
  SUBROUTINE HDF5_OUTPUTINIT6(JWD, JW)
    USE HDF5MOD
    USE TVDC, ONLY: NACD
    USE MAIN, ONLY: CDN, WDO, IWDO
    USE NAMESC, ONLY: CDNAME2
    IMPLICIT NONE
    
    INTEGER  :: JWD, JW
        
    HDF5_NACD = NACD(JW)
    CALL HDF5_OPEN_DATASET(WDO(JWD,4), 'dwo_'//SEGNUM(1:L), 'U')
    ALLOCATE(STRING1(NACD(JW)))
    DO I = 1, NACD(JW)
      STRING1(I) = CDNAME2(CDN(I,JW))
    END DO
    CALL HDF5_ADD_ATTRIBUTE(WDO(JWD,4), 'Derived Constituent file for segment ', IWDO(JWD))
    CALL HDF5_ADD_ATTRIBUTE(WDO(JWD,4), 'Output Variable names', NACD(JW),STRING1)
    DEALLOCATE(STRING1)
  END SUBROUTINE HDF5_OUTPUTINIT6
    
  SUBROUTINE HDF5_OUTPUTINIT7(JWD, JW, JS, STRING_TEMP1, STRING_TEMP2)
    USE HDF5MOD
    USE TVDC, ONLY: NAC, NACD, CN, CONSTITUENTS
    USE MAIN, ONLY: CDN, WDO2, IWDO, JFILE, DERIVED_CALC
    USE NAMESC, ONLY: CNAME2, CDNAME2
    IMPLICIT NONE 
    
    INTEGER      :: JWD, JW, JS
    CHARACTER(*) :: STRING_TEMP1
    CHARACTER(*) :: STRING_TEMP2
        
    HDF5_JWD = JWD
    ALLOCATE(STRING1(2 + NAC + NACD(JW)))
    HDF5_NACD = NACD(JW)
    CALL HDF5_OPEN_DATASET(WDO2(JFILE,1), STRING_TEMP1//SEGNUM(1:L), 'U')
    STRING1(1) = 'QWD (m3/s)'
    STRING1(2) = 'T (C)'
    
    IF (CONSTITUENTS) THEN
      DO I = 1, NAC
        STRING1(I+2) = CNAME2(CN(I))//' (g/m3)'
      END DO
    END IF
    IF (DERIVED_CALC) THEN
      DO I = 1, NACD(JW)
        STRING1(I+NAC+2) = CDNAME2(CDN(I,JW))
      END DO
    END IF
    
    CALL HDF5_ADD_ATTRIBUTE(WDO2(JFILE,1),'Output Variable names', 2+NAC+NACD(JW), STRING1)
    CALL HDF5_ADD_ATTRIBUTE(WDO2(JFILE,1), STRING_TEMP2//' No ', JS)
    CALL HDF5_ADD_ATTRIBUTE(WDO2(JFILE,1), 'Segment No ', IWDO(JWD))
    DEALLOCATE(STRING1)  
  END SUBROUTINE HDF5_OUTPUTINIT7
    
  !=========================================================================================================================== 
  ! WSEL
  SUBROUTINE HDF5_OUTPUTINIT8
    USE HDF5MOD
    USE MAIN, ONLY: WLFN
    IMPLICIT NONE
        
    CALL HDF5_OPEN_DATASET(WLFN, 'Water Level', 'U')
    ALLOCATE(STRING1(HDF5_N))
    HDF5_NSEG = 0
    DO JB = 1, NBR
      DO I = CUS(JB), DS(JB)
        WRITE(CHARFROMI,"(I4)") I
        HDF5_NSEG = HDF5_NSEG + 1
        STRING1(HDF5_NSEG) = 'Segment '// CHARFROMI// ' '
      END DO
    END DO
    CALL HDF5_ADD_ATTRIBUTE(WLFN, 'Variable Name', HDF5_N, STRING1)
    CALL HDF5_ADD_ATTRIBUTE(WLFN, 'Unit', '(m)')
    DEALLOCATE(STRING1)   
  END SUBROUTINE HDF5_OUTPUTINIT8
    
  !===========================================================================================================================  
  ! BRANCHES/SEGMENTS/FLUXES
  SUBROUTINE HDF5_OUTPUTINIT9
    USE HDF5MOD
    USE MAIN, ONLY: CDN, WDO2, IWDO, SNAPSHOT, PROFILE, IPRF, FLUX, FLXFN, NISNP, ISNP, KFNAME2, KBR, PRFFN, NIPRF, NDSP, &
                    CPRWBC, CDWBC, VOLUME_BALANCE, MASS_BALANCE, CONTOUR, FLOWBFN, DERIVED_CALC, MASSBFN
    USE TVDC, ONLY: NACTR, TRCN, NAC, NACD, CN, CONSTITUENTS; 
    USE NAMESC,  ONLY: CNAME, CDNAME, CNAME1, CNAME2, CDNAME2, CUNIT, CUNIT2
    USE KINETIC, ONLY: NAF, KFCN
    USE GEOMC,   ONLY: B, H; USE SCREENC, ONLY: JDAY;
    USE GDAYC,   ONLY: GDAY, YEAR, MONTH
    USE RSTART,  ONLY: PRFDP
    USE CEMAVars, ONLY: SEDIMENT_DIAGENESIS 
    IMPLICIT NONE
         
    ! BRANCHES                               
    CALL HDF5_OPEN_DATASET(BRANCHFN0+ 1, 'Branch Flow', 'U',3)
    CALL HDF5_OPEN_DATASET(BRANCHFN0+ 2, 'Branch Volume Balance', 'U',3)
    CALL HDF5_OPEN_DATASET(BRANCHFN0+ 3, 'Branch Energy Balance', 'U',3) 
    CALL HDF5_OPEN_DATASET(BRANCHFN0+ 4, 'Branch Mass Balance', 'U',3) 
    CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 1, 'Branch Dimension 1: JDAY',  'JDAY as Time Stamp (s)')
    CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 1, 'Branch Dimension 2: Variable 1 Column 0',  'Branch Flow [QSUM] (m3/s)')
    CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 1, 'Branch Dimension 2: Variable 2 Column 1',  'Evaporation Rate [EVBR] (m3/s)')
    CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 1, 'Branch Dimension 2: Variable 3 Column 2',  'Cumulative Evaporation [-VOLEV] (m3)')
    CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 1, 'Branch Dimension 2: Variable 4 Column 3',  'Precipitation [PR] (m/s)')
    CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 1, 'Branch Dimension 2: Variable 5 Column 4',  'Upstream Elevation [ELUH] (m)')
    CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 1, 'Branch Dimension 2: Variable 6 Column 5',  'Downstream Elevation [ELDH] (m)')
    WRITE(CHARFROMI, '(I4)') 5+KMX
    CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 1, 'Branch Dimension 2: Variable 7 Column 6 to'//CHARFROMI,' [QOUT] (m3/s)')
    WRITE(CHARFROMI, '(I4)') 6+KMX
    WRITE(CHARFROMJ, '(I4)') 5+KMX+NST
    CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 1, 'Branch Dimension 2: Variable 8 Column '//TRIM(CHARFROMI)//'to'//CHARFROMJ,' Structure Flow Discharge [QSTR] (m3/s)')
    CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 1, 'Branch Dimension 3: Number of Branch ', NBR)
    CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 2, 'Branch Dimension 1: JDAY',  'JDAY as Time Stamp (s)')
    CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 2, 'Branch Dimension 2: Variable 1 Column 0', 'Spatial Change [VOLSBR] (m3)')
    CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 2, 'Branch Dimension 2: Variable 2 Column 1', 'Temporal Change [VOLTBR] (m3)')
    CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 2, 'Branch Dimension 2: Variable 3 Column 2', 'Volume Error [VOLTBR-VOLSBR] (m3)')
    CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 2, 'Branch Dimension 2: Variable 4 Column 3', 'Percent Error [DLVBR*100.0] (%)')
    CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 2, 'Branch Dimension 3: Number of Branch ', NBR)
    CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 3, 'Branch Dimension 1: JDAY',  'JDAY as Time Stamp (s)')
    CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 3, 'Branch Dimension 2: Variable 1 Column 0', 'Spatially Integrated Energy [ESBR] (kJ)')
    CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 3, 'Branch Dimension 2: Variable 2 Column 1', 'Temporally Integrated Energy [ETBR] (kJ)')
    CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 3, 'Branch Dimension 2: Variable 3 Column 2', 'Energy Error [ESBR-ETBR] (kJ)')
    CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 3, 'Branch Dimension 2: Variable 4 Column 3', 'Percent Error [DLE*100.0] (%)')
    CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 3, 'Branch Dimension 3: Number of Branch ', NBR)
    CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 4, 'Branch Dimension 1: JDAY',  'JDAY as Time Stamp (s)')
    
    IF (CONSTITUENTS) THEN
    DO JC = 1, NAC
      WRITE(CHARFROMII, '(I4)') JC
      WRITE(CHARFROMH, '(I4)') (JC-1)*4
      WRITE(CHARFROMI, '(I4)') (JC-1)*4+1
      WRITE(CHARFROMJ, '(I4)') (JC-1)*4+2
      WRITE(CHARFROMK, '(I4)') (JC-1)*4+3
      CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 4, 'Branch Dimension 2: Variable'//CHARFROMII//TRIM(CNAME2(CN(JC)))//' Column '//CHARFROMH, 'Spatially Integrated Mass [CMBRS] (g/m3)')
      CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 4, 'Branch Dimension 2: Variable'//CHARFROMII//TRIM(CNAME2(CN(JC)))//' Column '//CHARFROMI, 'Temporally Integrated Mass [CMBRT] (g/m3)')
      CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 4, 'Branch Dimension 2: Variable'//CHARFROMII//TRIM(CNAME2(CN(JC)))//' Column '//CHARFROMJ, 'Mass Error [CMBRT-CMBRS] (g/m3)')
      CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 4, 'Branch Dimension 2: Variable'//CHARFROMII//TRIM(CNAME2(CN(JC)))//' Column '//CHARFROMK, 'Percent Error [DLMR*100.0] (%)')
    END DO
      
    DO JM = 1, NMC
      WRITE(CHARFROMII, '(I4)') NAC+JM
      WRITE(CHARFROMH, '(I4)') NAC*4+(JM-1)*4
      WRITE(CHARFROMI, '(I4)') NAC*4+(JM-1)*4+1
      WRITE(CHARFROMJ, '(I4)') NAC*4+(JM-1)*4+2
      WRITE(CHARFROMK, '(I4)') NAC*4+(JM-1)*4+3
      CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 4, 'Branch Dimension 2: Variable'//CHARFROMII//' Column '//CHARFROMH, 'Spatially Integrated Mass [MACMBRS] (g)')
      CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 4, 'Branch Dimension 2: Variable'//CHARFROMII//' Column '//CHARFROMI, 'Temporally Integrated Mass [MACMBRT] (g)')
      CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 4, 'Branch Dimension 2: Variable'//CHARFROMII//' Column '//CHARFROMJ, 'Mass Error [MACMBRT-MACMBRS] (g)')
      CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 4, 'Branch Dimension 2: Variable'//CHARFROMII//' Column '//CHARFROMK, 'Percent Error [DLMR*100.0] (%)')
    END DO 
    END IF
    CALL HDF5_ADD_ATTRIBUTE(BRANCHFN0+ 4, 'Branch Dimension 3: Number of Branch ', NBR)
      
    ! SEGMENTS
    DO JW = 1, NWB
      HDF5_JWB = JW
      IF (SNAPSHOT(JW)) THEN
        CALL HDF5_OPEN_DATASET(SEGMTFN(JW), 'Segments Output'//' wb'//HDF5_JWBS(JW), 'U', 3)
        ALLOCATE(STRING1(HDF5_WBNSEG(JW)))
        HDF5_NSEG = 0
        DO JB = BS(JW), BE(JW)
          DO I = CUS(JB), DS(JB)
            WRITE(CHARFROMJ,"(I4)") I
            HDF5_NSEG = HDF5_NSEG + 1
            STRING1(HDF5_NSEG) = 'Segment '// CHARFROMJ// ' '
          END DO
        END DO
        CALL HDF5_ADD_ATTRIBUTE(SEGMTFN(JW), 'Segments Dimension 1: JDAY',  'JDAY as Time Stamp (s)')
        CALL HDF5_ADD_ATTRIBUTE(SEGMTFN(JW), 'Segments Dimension 2: Segments', HDF5_WBNSEG(JW), STRING1)
        DEALLOCATE(STRING1)
        
        ALLOCATE(STRING1(HDF5_M))
        STRING1(1) = 'Water Level (ELWS) (m)'
        STRING1(2) = 'Water Surface Deviations (Z) (m)'
        STRING1(3) = 'Flow (Q) (m3/s)'
        STRING1(4) = 'Flow (QC) (m3/s)'
        STRING1(5) = 'X1'
        STRING1(6) = 'RN'
        STRING1(7) = 'RS'
        STRING1(8) = 'RB'
        STRING1(9) = 'RE'
        STRING1(10) = 'RC'
        STRING1(11) = 'KTI'
        CALL HDF5_ADD_ATTRIBUTE(SEGMTFN(JW), 'Segments Dimension 3: Variable Names', HDF5_M, STRING1)
        DEALLOCATE(STRING1)
      END IF 
    
      ! HDF5 FLUX OUTPUT                                                       
      IF (FLUX(JW)) THEN
        WRITE (SEGNUM,'(I0)') JW  
        SEGNUM = ADJUSTL(SEGNUM)  
        L      = LEN_TRIM(SEGNUM)
        IF (NAF(JW)>0) THEN 
          CALL HDF5_OPEN_DATASET(FLX(JW), FLXFN(JW)(1:3)//'_wb'//SEGNUM(1:L), 'U', 3)
          ALLOCATE (STRING12(11), STRING11(NISNP(JW)))
          DO I = 1, 11
            STRING12(I) = TITLE(I)
          END DO
          DO I = 1, NISNP(JW)
            WRITE(CHARFROMI,"(I4)") ISNP(I,JW)
            STRING11(I) = CHARFROMI
          END DO
          
          CALL HDF5_ADD_ATTRIBUTE(FLX(JW), 'JDAY', 3, JDAY)
          WRITE(CHARFROMI, "(I4)") GDAY
          WRITE(CHARFROMJ, "(I4)") YEAR
          CALL HDF5_ADD_ATTRIBUTE(FLX(JW), 'New date', MONTH//' '//CHARFROMI//', '//CHARFROMJ)
          WRITE(CHARFROMI, "(I4)") INT(JDAY)
          WRITE(CHARFROMJ, "(F0.2)") (JDAY-INT(JDAY))*24.0
          CALL HDF5_ADD_ATTRIBUTE(FLX(JW), 'FLUX Title', 11, STRING12)
          CALL HDF5_ADD_ATTRIBUTE(FLX(JW), 'FLUX Julian Date', CHARFROMI//'days '//CHARFROMJ//' hours')
          CALL HDF5_ADD_ATTRIBUTE(FLX(JW), 'FLUX Dimension 2: Segments', NISNP(JW), STRING11)
          !
          DO I = 1, NAF(JW)
            WRITE(CHARFROMKK(I),"(I3)") I
            INTS0(1) = (I-1)*(KBR(JW)-KTWB(JW)+1)
            INTS0(2) = I*(KBR(JW)-KTWB(JW)+1)-1
            CALL HDF5_ADD_ATTRIBUTE(FLX(JW), 'FLUX Dimension 1: KFN'//CHARFROMKK(I)//' '//TRIM(KFNAME(KFCN(I,JW)))//' No. Begin & End', 2, INTS0)
          END DO
          
          CALL HDF5_ADD_ATTRIBUTE(FLX(JW), 'FLUX Dimension 3: JDAY', 'As Time Stamp')
          CALL HDF5_ADD_ATTRIBUTE(FLX(JW), 'No of KFN', NAF(JW))
          CALL HDF5_OPEN_DATASET(FLX2(JW),'kflux_wb'//SEGNUM(1:L), 'U')
          DEALLOCATE(STRING12,STRING11)
          !  
          ALLOCATE(STRING1(1+NAF(JW)))
          STRING1(1) = 'ELTM (days)'
          DO JF = 1, NAF(JW)
            STRING1(1+JF) = KFNAME2(KFCN(JF,JW))
          END DO
          CALL HDF5_ADD_ATTRIBUTE(FLX2(JW), 'FLUX2_Variable Names', 1+NAF(JW), STRING1)              
          DEALLOCATE(STRING1)
        END IF
      END IF
      ! END FLUX OUTPUT
    END DO
    
    ! FLOW BALANCE                                        
    DO JW = 1, NWB  
      IF(VOLUME_BALANCE(JW) .AND. CONTOUR(JW)) THEN   
        CALL HDF5_OPEN_DATASET(FLOWBFN,'flowbal','U') 
        ALLOCATENO = 9
        IF(VOLUME_BALANCE(JW)) ALLOCATENO = 10
        ALLOCATE(STRING1(ALLOCATENO))
        STRING1(1) = 'WB'
        STRING1(2) = 'VOLIN(m3)'
        STRING1(3) = 'VOLPR(m3)'
        STRING1(4) = 'VOLOUT(m3)'
        STRING1(5) = 'VOLWD(m3)'
        STRING1(6) = 'VOLEV(m3)'
        STRING1(7) = 'VOLDT(m3)'
        STRING1(8) = 'VOLTRB(m3)'
        STRING1(9) = 'VOLICE(m3)'
        IF(VOLUME_BALANCE(JW)) STRING1(10) = '%VOLerror'
        CALL HDF5_ADD_ATTRIBUTE(FLOWBFN, 'Variable Names', ALLOCATENO, STRING1)
        DEALLOCATE(STRING1)
        EXIT  
      END IF  
    END DO 
    !        
    DO JW = 1, NWB  
      IF (MASS_BALANCE(JW) .AND. CONTOUR(JW) .AND. DERIVED_CALC) THEN  
        CALL HDF5_OPEN_DATASET(MASSBFN,'massbal','U')
        ALLOCATENO = 26
        IF(SEDIMENT_DIAGENESIS) ALLOCATENO = 29
        ALLOCATE(STRING1(ALLOCATENO))
        STRING1(1) = 'WB'
        STRING1(2) = 'TP-Waterbody(kg)'
        STRING1(3) = 'TP-Sediment(kg)'
        STRING1(4) = 'TP-Plants(kg)'
        STRING1(5) = 'OutflowTP(kg)'
        STRING1(6) = 'TributaryTP(kg)'
        STRING1(7) = 'DistributedTributaryTP(kg)'
        STRING1(8) = 'WithdrawalTP(kg)'
        STRING1(9) = 'PrecipitationTP(kg)'
        STRING1(10) = 'InflowTP(kg)'
        STRING1(11) = 'AtmosphericDepositionP(kg)'   
        STRING1(12) = 'SED+SOD_PRelease(kg)'
        STRING1(13) = 'PFluxtoSediments(kg)'
        IF (SEDIMENT_DIAGENESIS) THEN
          STRING1(14) = 'SedimentDiagenesisPFlux(kg)'
          STRING1(15) = 'TN-Waterbody(kg)'
          STRING1(16) = 'TN-Sediment(kg)'
          STRING1(17) = 'TN-Plants(kg)'
          STRING1(18) = 'OutflowTN(kg)'
          STRING1(19) = 'TributaryTN(kg)'
          STRING1(20) = 'DistributedTributaryTN(kg)'
          STRING1(21) = 'WithdrawalTN(kg)'
          STRING1(22) = 'PrecipitationTN(kg)'
          STRING1(23) = 'InflowTN(kg)'
          STRING1(24) = 'AtmosphericDepositionN(kg)'
          STRING1(25) = 'NH3GasLoss(kg)'
          STRING1(26) = 'SED+SOD_NRelease(kg)'
          STRING1(27) = 'NFluxtoSediments(kg)'
          STRING1(28) = 'SedimentDiagenesisNH4Flux(kg)'
          STRING1(29) = 'SedimentDiagenesisNO3Flux(kg)'
        ELSE
          STRING1(14) = 'TN-Waterbody(kg)'
          STRING1(15) = 'TN-Sediment(kg)'
          STRING1(16) = 'TN-Plants(kg)'
          STRING1(17) = 'OutflowTN(kg)'
          STRING1(18) = 'TributaryTN(kg)'
          STRING1(19) = 'DistributedTributaryTN(kg)'
          STRING1(20) = 'WithdrawalTN(kg)'
          STRING1(21) = 'PrecipitationTN(kg)'
          STRING1(22) = 'InflowTN(kg)'
          STRING1(23) = 'AtmosphericDepositionN(kg)'
          STRING1(24) = 'NH3GasLoss(kg)'
          STRING1(25) = 'SED+SOD_NRelease(kg)'
          STRING1(26) = 'NFluxtoSediments(kg)'
        END IF
        CALL HDF5_ADD_ATTRIBUTE(MASSBFN, 'Variable Names', ALLOCATENO, STRING1)
        DEALLOCATE(STRING1)
        EXIT  
      END IF  
    END DO     
  END SUBROUTINE HDF5_OUTPUTINIT9    
    
  !=========================================================================================================================== 
  ! WSL
  SUBROUTINE HDF5_OUTPUTA3(JDAY)
    USE HDF5MOD
    USE GEOMC, ONLY: ELWS; USE MAIN, ONLY: WLFN
    IMPLICIT NONE
    
    REAL :: JDAY
    
    ALLOCATE(DATASET2(HDF5_N,1))
    HDF5_NSEG = 0
    DO JB = 1, NBR
      DO I = US(JB), DS(JB)
        HDF5_NSEG = HDF5_NSEG + 1
        CALL XINT(ELWS(I),3, DATASET2(HDF5_NSEG,1))
      END DO
    END DO
    CALL HDF5_ADD_DATA(WLFN, HDF5_N, 1, DATASET2)
    CALL HDF5_ADD_DATE(WLFN, JDAY)   
    DEALLOCATE(DATASET2)
  END SUBROUTINE HDF5_OUTPUTA3
     
  !=========================================================================================================================== 
  ! TIME SERIES
  SUBROUTINE HDF5_OUTPUTA4(TVOLAVG)
    USE HDF5MOD
    USE MAIN, ONLY: ITSR, ETSR, ICE_COMPUTATION, RN, RS, RANLW, SEDIMENT_CALC, CDN, NISNP, ISNP, DERIVED_CALC
    USE GEOMC, ONLY: DEPTHB, ELWS, BI; 
    USE KINETIC, ONLY: REAER, QC, GAMMA, EPD, SED, SEDC, SEDN, SEDP, NAF, KFCN
    USE TVDC, ONLY: NAC, NACD, SRON, CN, CONSTITUENTS
    USE SURFHE, ONLY: SHADE, RB, RE, RC
    USE MACROPHYTEC, ONLY: MAC       
    IMPLICIT NONE
    
    REAL :: TVOLAVG
    REAL :: JDAY
      
    DO JW = 1, NWB
      HDF5_JWB = JW
      IF (HDF5_JTS(JW)>0) THEN
        HDF5_ITS = 0
        ALLOCATE(DATASET3(HDF5_JTS(JW), HDF5_A(JW),1))
        DATASET3(:,:,:) = -99.0
        
        DO J = 1,NIKTSR
          I = ITSR(J)
          ETSR(J) = 0  
          IF (ETSR(J) < 0) THEN  
            K = INT(ABS(ETSR(J)))  
          ELSE  
            DO K=KTWB(JW),KB(I)  
              IF (DEPTHB(K,I) > ETSR(J)) EXIT  
            END DO  
            IF(K > KB(I)) CYCLE  
          END IF
            
          IF (I >= US(BS(JW)) .AND. I <= DS(BE(JW))) THEN
            HDF5_ITS = HDF5_ITS+1 
            CALL XINT(DLT, 3,     DATASET3(HDF5_ITS,1,1))
            CALL XINT(ELWS(I), 3, DATASET3(HDF5_ITS,2,1))
            CALL XINT(T1(K,I), 3, DATASET3(HDF5_ITS,3,1))
            CALL XINT(U(K,I), 3,  DATASET3(HDF5_ITS,4,1))
            CALL XINT(QC(I), 3,   DATASET3(HDF5_ITS,5,1))
            CALL XINT(SRON(JW)*1.06, 3,   DATASET3(HDF5_ITS,6,1))
            CALL XINT(GAMMA(K,I), 3,      DATASET3(HDF5_ITS,7,1))
            CALL XINT(DEPTHB(KB(I),I), 3, DATASET3(HDF5_ITS,8,1))
            CALL XINT(BI(KTWB(JW),I), 3,  DATASET3(HDF5_ITS,9,1))
            CALL XINT(SHADE(I), 3,        DATASET3(HDF5_ITS,10,1))
                   
            IF (ICE_COMPUTATION) THEN
              CALL XINT(ICETH(I), 3, DATASET3(HDF5_ITS,11,1))
              CALL XINT(TVOLAVG, 3,  DATASET3(HDF5_ITS,12,1))
              CALL XINT(RN(I), 3,    DATASET3(HDF5_ITS,13,1))
              CALL XINT(RS(I), 3,    DATASET3(HDF5_ITS,14,1))
              CALL XINT(RANLW(JW), 3,DATASET3(HDF5_ITS,15,1))
              CALL XINT(RB(I), 3,    DATASET3(HDF5_ITS,16,1))
              CALL XINT(RE(I), 3,    DATASET3(HDF5_ITS,17,1))
              CALL XINT(RC(I), 3,    DATASET3(HDF5_ITS,18,1))
              CALL XINT(REAER(I)*86400., 3, DATASET3(HDF5_ITS,19,1))   
              
            IF (CONSTITUENTS) THEN
              DO JC = 1, NAC
                CALL XINT(C2(K,I,CN(JC))*CMULT(CN(JC)), 3, DATASET3(HDF5_ITS,19+JC,1))
              END DO
              DO JE = 1, NEP
                CALL XINT(EPD(K,I,JE), 3, DATASET3(HDF5_ITS,19+NAC+JE,1))
              END DO
              DO JM = 1, NMC
                CALL XINT(MAC(K,I,JM), 3, DATASET3(HDF5_ITS,19+NAC+NEP+JM,1))
              END DO
                
              IF (SEDIMENT_CALC(JW)) THEN 
                CALL XINT(SED(K,I), 3,  DATASET3(HDF5_ITS,19+NAC+NEP+NMC+1,1))
                CALL XINT(SEDP(K,I), 3, DATASET3(HDF5_ITS,19+NAC+NEP+NMC+2,1))
                CALL XINT(SEDN(K,I), 3, DATASET3(HDF5_ITS,19+NAC+NEP+NMC+3,1))
                CALL XINT(SEDC(K,I), 3, DATASET3(HDF5_ITS,19+NAC+NEP+NMC+4,1))
                DO JD = 1, NACD(JW)
                  CALL XINT(CD(K,I,CDN(JD,JW))*CDMULT(CDN(JD,JW)), 3, DATASET3(HDF5_ITS,19+NAC+NEP+NMC+4+JD,1))
                END DO
                DO JF = 1, NAF(JW)
                  CALL XINT(KF(K,I,(KFCN(JF,JW)))*VOL(K,I)/1000./DAY, 3, DATASET3(HDF5_ITS,19+NAC+NEP+NMC+4+NACD(JW)+JF,1))
                END DO
                !
                DO JA = 1, NAL
                  CALL XINT(APLIM(K,I,JA), 3, DATASET3(HDF5_ITS,19+NAC+NEP+NMC+4+NACD(JW)+NAF(JW)+JA,1))
                END DO
                DO JA = 1, NAL
                  CALL XINT(ANLIM(K,I,JA), 3, DATASET3(HDF5_ITS,19+NAC+NEP+NMC+4+NACD(JW)+NAF(JW)+NAL+JA,1))
                END DO
                DO JA = 1, NAL
                  CALL XINT(ALLIM(K,I,JA), 3, DATASET3(HDF5_ITS,19+NAC+NEP+NMC+4+NACD(JW)+NAF(JW)+2*NAL+JA,1))
                END DO
                
              ELSE
                DO JD = 1, NACD(JW)
                  CALL XINT(CD(K,I,CDN(JD,JW))*CDMULT(CDN(JD,JW)), 3, DATASET3(HDF5_ITS,19+NAC+NEP+NMC+JD,1))
                END DO
                DO JF = 1, NAF(JW)
                  CALL XINT(KF(K,I,(KFCN(JF,JW)))*VOL(K,I)/1000./DAY, 3, DATASET3(HDF5_ITS,19+NAC+NEP+NMC+NACD(JW)+JF,1))
                END DO
                !
                DO JA = 1, NAL
                  CALL XINT(APLIM(K,I,JA), 3, DATASET3(HDF5_ITS,19+NAC+NEP+NMC+NACD(JW)+NAF(JW)+JA,1))
                END DO
                DO JA = 1, NAL
                  CALL XINT(ANLIM(K,I,JA), 3, DATASET3(HDF5_ITS,19+NAC+NEP+NMC+NACD(JW)+NAF(JW)+NAL+JA,1))
                END DO
                DO JA = 1, NAL
                  CALL XINT(ALLIM(K,I,JA), 3, DATASET3(HDF5_ITS,19+NAC+NEP+NMC+NACD(JW)+NAF(JW)+2*NAL+JA,1))
                END DO
              END IF
            END IF
            !
            ELSE
              CALL XINT(TVOLAVG, 3,  DATASET3(HDF5_ITS,11,1))
              CALL XINT(RN(I), 3,    DATASET3(HDF5_ITS,12,1))
              CALL XINT(RS(I), 3,    DATASET3(HDF5_ITS,13,1))
              CALL XINT(RANLW(JW), 3,DATASET3(HDF5_ITS,14,1))
              CALL XINT(RB(I), 3,    DATASET3(HDF5_ITS,15,1))
              CALL XINT(RE(I), 3,    DATASET3(HDF5_ITS,16,1))
              CALL XINT(RC(I), 3,    DATASET3(HDF5_ITS,17,1))
              CALL XINT(REAER(I)*86400., 3, DATASET3(HDF5_ITS,18,1))
              
            IF (CONSTITUENTS) THEN
              DO JC = 1, NAC
                CALL XINT(C2(K,I,CN(JC))*CMULT(CN(JC)), 3, DATASET3(HDF5_ITS,18+JC,1))
              END DO
              DO JE = 1, NEP
                CALL XINT(EPD(K,I,JE), 3, DATASET3(HDF5_ITS,18+NAC+JE,1))
              END DO
              DO JM = 1, NMC
                CALL XINT(MAC(K,I,JM), 3, DATASET3(HDF5_ITS,18+NAC+NEP+JM,1))
              END DO
                  
              IF (SEDIMENT_CALC(JW)) THEN 
                CALL XINT(SED(K,I), 3,  DATASET3(HDF5_ITS,18+NAC+NEP+NMC+1,1))
                CALL XINT(SEDP(K,I), 3, DATASET3(HDF5_ITS,18+NAC+NEP+NMC+2,1))
                CALL XINT(SEDN(K,I), 3, DATASET3(HDF5_ITS,18+NAC+NEP+NMC+3,1))
                CALL XINT(SEDC(K,I), 3, DATASET3(HDF5_ITS,18+NAC+NEP+NMC+4,1))
                DO JD = 1, NACD(JW)
                  CALL XINT(CD(K,I,CDN(JD,JW))*CDMULT(CDN(JD,JW)), 3, DATASET3(HDF5_ITS,18+NAC+NEP+NMC+4+JD,1))
                END DO
                DO JF = 1, NAF(JW)
                  CALL XINT(KF(K,I,(KFCN(JF,JW)))*VOL(K,I)/1000./DAY, 3, DATASET3(HDF5_ITS,18+NAC+NEP+NMC+4+NACD(JW)+JF,1))
                END DO
                !
                DO JA = 1, NAL
                  CALL XINT(APLIM(K,I,JA), 3, DATASET3(HDF5_ITS,18+NAC+NEP+NMC+4+NACD(JW)+NAF(JW)+JA,1))
                END DO
                DO JA = 1, NAL
                  CALL XINT(ANLIM(K,I,JA), 3, DATASET3(HDF5_ITS,18+NAC+NEP+NMC+4+NACD(JW)+NAF(JW)+NAL+JA,1))
                END DO
                DO JA = 1, NAL
                  CALL XINT(ALLIM(K,I,JA), 3, DATASET3(HDF5_ITS,18+NAC+NEP+NMC+4+NACD(JW)+NAF(JW)+2*NAL+JA,1))
                END DO
                
              ELSE
                DO JD = 1, NACD(JW)
                  CALL XINT(CD(K,I,CDN(JD,JW))*CDMULT(CDN(JD,JW)), 3, DATASET3(HDF5_ITS,18+NAC+NEP+NMC+JD,1))
                END DO
                DO JF = 1, NAF(JW)
                  CALL XINT(KF(K,I,(KFCN(JF,JW)))*VOL(K,I)/1000./DAY, 3, DATASET3(HDF5_ITS,18+NAC+NEP+NMC+NACD(JW)+JF,1))
                END DO
                !
                DO JA = 1, NAL
                  CALL XINT(APLIM(K,I,JA), 3, DATASET3(HDF5_ITS,18+NAC+NEP+NMC+NACD(JW)+NAF(JW)+JA,1))
                END DO
                DO JA = 1, NAL
                  CALL XINT(ANLIM(K,I,JA), 3, DATASET3(HDF5_ITS,18+NAC+NEP+NMC+NACD(JW)+NAF(JW)+NAL+JA,1))
                END DO
                DO JA = 1, NAL
                  CALL XINT(ALLIM(K,I,JA), 3, DATASET3(HDF5_ITS,18+NAC+NEP+NMC+NACD(JW)+NAF(JW)+2*NAL+JA,1))
                END DO
              END IF
            END IF
            !
            END IF
          END IF
        END DO 
        CALL HDF5_ADD_DATA(TSRFN0+JW, HDF5_JTS(JW), HDF5_A(JW), 1, DATASET3, 'Z') 
        DEALLOCATE(DATASET3)
      END IF
    END DO
  END SUBROUTINE HDF5_OUTPUTA4
     
  !=========================================================================================================================== 
  ! SEGMENTS
  SUBROUTINE HDF5_OUTPUTA5(JW, JDAY)
    USE HDF5MOD
    USE MAIN, ONLY: Q, X1, RN, RS, ICE_COMPUTATION
    USE GEOMC, ONLY: ELWS, Z
    USE KINETIC, ONLY: QC; USE SURFHE, ONLY: SHADE, RB, RE, RC; 
    USE TVDC,    ONLY: NAC, NACD, SRON, CN     
    IMPLICIT NONE
    
    INTEGER    :: JW
    REAL       :: JDAY
      
    HDF5_JWB = JW
    ALLOCATE(DATASET3(HDF5_M, HDF5_WBNSEG(JW), 1)) 
    DATASET3(:,:,:) = -99.0
    HDF5_NSEG = 0
    DO I = US(BS(JW))-1, DS(BE(JW))+1
      HDF5_NSEG = HDF5_NSEG + 1
      CALL XINT(ELWS(I),3,DATASET3(1, HDF5_NSEG, 1))
      CALL XINT(Z(I), 3,  DATASET3(2, HDF5_NSEG, 1))
      CALL XINT(Q(I), 3,  DATASET3(3, HDF5_NSEG, 1))
      CALL XINT(QC(I), 3, DATASET3(4, HDF5_NSEG, 1))
      CALL XINT(X1(I), 3, DATASET3(5, HDF5_NSEG, 1))
      CALL XINT(RN(I), 3, DATASET3(6, HDF5_NSEG, 1))
      CALL XINT(RS(I), 3, DATASET3(7, HDF5_NSEG, 1))
      CALL XINT(RB(I), 3, DATASET3(8, HDF5_NSEG, 1))
      CALL XINT(RE(I), 3, DATASET3(9, HDF5_NSEG, 1))
      CALL XINT(RC(I), 3, DATASET3(10, HDF5_NSEG, 1)) 
      CALL XINT(KTI(I),3, DATASET3(11, HDF5_NSEG, 1))
    END DO     
    CALL HDF5_ADD_DATA(SEGMTFN(JW), HDF5_M, HDF5_WBNSEG(JW), 1, DATASET3, 'Z')
    CALL HDF5_ADD_DATE(JW, SEGMTFN(JW), JDAY)
    DEALLOCATE(DATASET3)
  END SUBROUTINE HDF5_OUTPUTA5
     
  !===========================================================================================================================
  ! FLUXES
  SUBROUTINE HDF5_OUTPUTA7(JW, JDAY)  
    USE HDF5MOD
    USE MAIN, ONLY: NISNP, KFJW
    USE KINETIC, ONLY: NAF, KFCN
    USE RSTART,  ONLY: ELTMF
    IMPLICIT NONE
    
    INTEGER   :: JW
    REAL      :: JDAY
    INTEGER   :: JAF
     
    HDF5_JWB = JW
    ALLOCATE(DATASET3(1, NISNP(JW), NAF(JW)*(KBR(JW)-KTWB(JW)+1)))
    DATASET3(:,:,:) = -99.0
    DO JAF = 1, NAF(JW)
      DO I=1,NISNP(JW)  
        DO K=KTWB(JW),KB(ISNP(I,JW))  
          CALL XINT(KFS(K,ISNP(I,JW),KFCN(JAF,JW)), 3, DATASET3(1,I, (JAF-1)*(KBR(JW)-KTWB(JW)+1)+K-KTWB(JW)+1))   
        END DO  
      END DO 
    END DO
    
    ! HDF5 FLUX DATA 
    IF (NAF(JW)>0 .AND. NISNP(JW)>0 .AND. KBR(JW)-KTWB(JW)+1>0) THEN 
      CALL HDF5_ADD_DATA(FLX(JW), 1, NISNP(JW), NAF(JW)*(KBR(JW)-KTWB(JW)+1), DATASET3, 'X')
      CALL HDF5_ADD_DATE(JW, FLX(JW), JDAY)
    END IF
    DEALLOCATE(DATASET3)
    !
    ALLOCATE(DATASET2(1+NAF(JW),1))
    CALL XINT(ELTMF(JW)/DAY, 3, DATASET2(1,1))
    DO K = 1, NAF(JW)
      CALL XINT(KFJW(JW,KFCN(K,JW)), 3, DATASET2(1+K,1))
    END DO
    CALL HDF5_ADD_DATA(FLX2(JW), 1+NAF(JW), 1, DATASET2) 
    DEALLOCATE(DATASET2)
  END SUBROUTINE HDF5_OUTPUTA7
             
  !===========================================================================================================================
  ! WITHDRAWAL
  SUBROUTINE HDF5_OUTPUTA8(J, JDAY, QOUTLET, TOUTLET)
    USE HDF5MOD
    USE MAIN, ONLY: IWD, IWDO, CDN, WDO, WDO2, DERIVED_CALC, PUMPON, QWDO, CWDO, CDWDO, jfile
    USE TVDC, ONLY: NAC, NACD, CONSTITUENTS, CN, QWD, TWDO 
    USE SELWC, ONLY: QSTR, NSTR, TAVG, TAVGW, CAVG, CDAVG, CAVGW, CDAVGW
    USE STRUCTURES, ONLY: JBUSP, JBUPU, JBUPI, JBUGT, IUSP, IUPU, IUPI, IUGT, &
                          JWUSP, JWUPU, JWUPI, JWUGT, QSP, QPU, QPI, QGT,     &
                          LATERAL_SPILLWAY, LATERAL_PIPE, LATERAL_GATE, LATERAL_PUMP
    IMPLICIT NONE
    
    INTEGER  :: J, JJ
    REAL     :: JDAY
    INTEGER  :: JSSS(100)
    REAL     :: QOUTLET(100), TOUTLET(100)
        
    HDF5_JFILE = 0
    HDF5_JWD = J 
    HDF5_E = 0                     
    HDF5_T = 0
    JWD = 0
    DO JW = 1,NWB  
      DO JB = BS(JW),BE(JW)
        IF (DS(JB) == IWDO(J)) THEN  
          DO JS = 1,NSTR(JB)
            HDF5_JFILE = HDF5_JFILE+1 
            jfile = hdf5_jfile
            ALLOCATE(TEMPSET2(2+NAC+NACD(JW),1))                           
            TEMPSET2(1,1) = QSTR(JS,JB)                              
            TEMPSET2(2,1) = TAVG(JS,JB)                         
            IF (JB==1) THEN
              HDF5_T = JS
            ELSE IF (JB >= 2) THEN
              HDF5_T = 0
              DO I = 1, JB-1
                HDF5_T = HDF5_T + NSTR(I)
              END DO
              HDF5_T = HDF5_T + JS
            END IF
            !
            IF (CONSTITUENTS) THEN
              DO JC = 1, NAC
                TEMPSETC2((HDF5_T-1)*NAC+JC,1) = CAVG(JS,JB,CN(JC))
                TEMPSET2(2+JC,1) = CAVG(JS,JB,CN(JC))                            
              END DO
            END IF
            IF (DERIVED_CALC) THEN
              DO JD = 1, NACD(JW)
                TEMPSETD2((HDF5_T-1)*NACD(JW)+JD,1) = CDAVG(JS,JB,CDN(JD,JW))
                TEMPSET2(2+NAC+JD,1) = CDAVG(JS,JB,CDN(JD,JW))                    
              END DO
            END IF
            
            IF(ALLOCATED(DATASET2)) DEALLOCATE(DATASET2)  
            ALLOCATE(DATASET2(2+NAC+NACD(JW),1))
            DO I = 1, 2+NAC+NACD(JW)
              CALL XINT(TEMPSET2(I,1), 3, DATASET2(I,1))
            END DO
            CALL HDF5_ADD_DATA(WDO2(HDF5_JFILE,1), 2+NAC+NACD(JW), 1, DATASET2)
            IF(ALLOCATED(TEMPSET2)) DEALLOCATE(TEMPSET2)
            IF(ALLOCATED(DATASET2)) DEALLOCATE(DATASET2)
          END DO
          JSSS(JB) = NSTR(JB)
        END IF  
      END DO  
    END DO 
    
    !-----------------------------------------------------------------------------  
    HDF5_E = HDF5_E + HDF5_T 
    DO JW = 1,NWB                   
      if(iwdo(j)>=us(bs(jw)) .and. iwdo(j)<=ds(be(JW))) exit  
    END DO
    DO JJ = 1,NWD  
      JWD = JWD+1
      IF (IWD(JWD) == IWDO(J)) THEN  
        HDF5_JFILE = HDF5_JFILE+1
        jfile = hdf5_jfile
        IF(ALLOCATED(TEMPSET2)) DEALLOCATE(TEMPSET2)
        ALLOCATE(TEMPSET2(2+NAC+NACD(JW),1))                       
        HDF5_T = HDF5_E + JJ 
        TEMPSET2(1,1) = QWD(JWD)                                     
        TEMPSET2(2,1) = TAVGW(JWD)                        
        IF (CONSTITUENTS) THEN
          DO JC = 1, NAC
            TEMPSETC2((HDF5_T-1)*NAC+JC,1) = CAVGW(JWD,CN(JC))
            TEMPSET2(2+JC,1) = CAVGW(JWD,CN(JC))                         
          END DO
        END IF
        IF (DERIVED_CALC) THEN
          DO JD = 1, NACD(JW)
            TEMPSETD2((HDF5_T-1)*NACD(JW)+JD,1) = CDAVGW(JWD,CDN(JD,JW))
            TEMPSET2(2+NAC+JD,1) = CDAVGW(JWD,CDN(JD,JW))                  
          END DO
        END IF
        
        IF(ALLOCATED(DATASET2)) DEALLOCATE(DATASET2)  
        ALLOCATE (DATASET2(2+NAC+NACD(JW),1))
        DO I = 1, 2+NAC+NACD(JW)
          CALL XINT(TEMPSET2(I,1), 3, DATASET2(I,1))
        END DO
        CALL HDF5_ADD_DATA(WDO2(HDF5_JFILE,1), 2+NAC+NACD(JW), 1, DATASET2)
        IF(ALLOCATED(TEMPSET2)) DEALLOCATE(TEMPSET2) 
        IF(ALLOCATED(DATASET2)) DEALLOCATE(DATASET2)
      END IF
    END DO
      
    !-----------------------------------------------------------------------------    
    HDF5_E = HDF5_E + NWD
    DO JS = 1,NSP  
      IF (LATERAL_SPILLWAY(JS)) THEN  
        JWD = JWD+1
      ELSE
        JSSS(JBUSP(JS)) = JSSS(JBUSP(JS))+1
      END IF
        
      IF (IWDO(J) == IUSP(JS)) THEN 
        HDF5_JFILE = HDF5_JFILE+1 
        jfile = hdf5_jfile
        IF(ALLOCATED(TEMPSET2)) DEALLOCATE(TEMPSET2)
        ALLOCATE(TEMPSET2(2+NAC+NACD(JWUSP(JS)),1))                     
        HDF5_T = HDF5_E + JS 
        TEMPSET2(1,1) = QSP(JS)
        IF (LATERAL_SPILLWAY(JS)) THEN     
          TEMPSET2(2,1) = TAVGW(JWD)                                 
          IF (CONSTITUENTS) THEN
            DO JC = 1, NAC
              TEMPSETC2((HDF5_T-1)*NAC+JC,1) = CAVGW(JWD,CN(JC))
              TEMPSET2(2+JC,1) = CAVGW(JWD,CN(JC))                            
            END DO
          END IF           
          IF (DERIVED_CALC) THEN
            DO JD = 1, NACD(JWUSP(JS))
              TEMPSETD2((HDF5_T-1)*NACD(JWUSP(JS))+JD,1) = CDAVGW(JWD,CDN(JD,JWUSP(JS)))
              TEMPSET2(2+NAC+JD,1) = CDAVGW(JWD,CDN(JD,JWUSP(JS)))                             
            END DO
          END IF
          !
        ELSE
          TEMPSET2(2,1) = TAVG(JSSS(JBUSP(JS)),JBUSP(JS))
          IF (CONSTITUENTS) THEN
            DO JC = 1, NAC
              TEMPSETC2((HDF5_T-1)*NAC+JC,1) = CAVG(JSSS(JBUSP(JS)),JBUSP(JS),CN(JC))
              TEMPSET2(2+JC,1) = CAVG(JSSS(JBUSP(JS)),JBUSP(JS),CN(JC))                          
            END DO
          END IF
          IF (DERIVED_CALC) THEN
            DO JD = 1, NACD(JWUSP(JS))
              TEMPSETD2((HDF5_T-1)*NACD(JWUSP(JS))+JD,1) = CDAVG(JSSS(JBUSP(JS)),JBUSP(JS),CDN(JD,JWUSP(JS)))
              TEMPSET2(2+NAC+JD,1) = CDAVG(JSSS(JBUSP(JS)),JBUSP(JS),CDN(JD,JWUSP(JS)))                  
            END DO
          END IF
        END IF
        
        IF(ALLOCATED(DATASET2)) DEALLOCATE(DATASET2)  
        ALLOCATE(DATASET2(2+NAC+NACD(JWUSP(JS)),1))
        DO I = 1, 2+NAC+NACD(JWUSP(JS))
          CALL XINT(TEMPSET2(I,1), 3, DATASET2(I,1))
        END DO
        CALL HDF5_ADD_DATA(WDO2(HDF5_JFILE,1), 2+NAC+NACD(JWUSP(JS)), 1, DATASET2)
        IF(ALLOCATED(TEMPSET2)) DEALLOCATE(TEMPSET2)
        IF(ALLOCATED(DATASET2)) DEALLOCATE(DATASET2)
      END IF
    END DO
    
    !-----------------------------------------------------------------------------    
    HDF5_E = HDF5_E + NSP 
    DO JS = 1,NPU  
      IF (LATERAL_PUMP(JS)) THEN  
        JWD = JWD+1 
      ELSE
        JSSS(JBUPU(JS)) = JSSS(JBUPU(JS))+1       
      END IF
        
      IF (IWDO(J) == IUPU(JS)) THEN  
        HDF5_JFILE = HDF5_JFILE+1  
        jfile = hdf5_jfile
        IF(ALLOCATED(TEMPSET2)) DEALLOCATE(TEMPSET2)
        ALLOCATE(TEMPSET2(2+NAC+NACD(JWUPU(JS)),1))
        HDF5_T = HDF5_E + JS 
        IF (PUMPON(JS)) THEN  
          TEMPSET2(1,1) = QPU(JS)                      
        ELSE  
          TEMPSET2(1,1) = 0.0                          
        END IF 
          
        IF (LATERAL_PUMP(JS)) THEN
          TEMPSET2(2,1) = TAVGW(JWD)          
          IF (CONSTITUENTS) THEN
            DO JC = 1, NAC
              TEMPSETC2((HDF5_T-1)*NAC+JC,1) = CAVGW(JWD,CN(JC))
              TEMPSET2(2+JC,1) = CAVGW(JWD,CN(JC))                      
            END DO
          END IF
          IF (DERIVED_CALC) THEN
            DO JD = 1, NACD(JWUPU(JS))
              TEMPSETD2((HDF5_T-1)*NACD(JWUPU(JS))+JD,1) = CDAVGW(JWD,CDN(JD,JWUPU(JS)))
              TEMPSET2(2+NAC+JD,1) = CDAVGW(JWD,CDN(JD,JWUPU(JS)))                    
            END DO
          END IF
          !
        ELSE
          TEMPSET2(2,1) = TAVG(JSSS(JBUPU(JS)),JBUPU(JS))      
          IF (CONSTITUENTS) THEN
            DO JC = 1, NAC
              TEMPSETC2((HDF5_T-1)*NAC+JC,1) = CAVG(JSSS(JBUPU(JS)),JBUPU(JS),CN(JC))
              TEMPSET2(2+JC,1) = CAVG(JSSS(JBUPU(JS)),JBUPU(JS),CN(JC))                       
            END DO
          END IF
          IF (DERIVED_CALC) THEN
            TEMPSETD2((HDF5_T-1)*NACD(JWUPU(JS))+JD,1) = CDAVG(JSSS(JBUPU(JS)),JBUPU(JS),CDN(JD,JWUPU(JS)))
            TEMPSET2(2+NAC+JD,1) = CDAVG(JSSS(JBUPU(JS)),JBUPU(JS),CDN(JD,JWUPU(JS)))   
          END IF  
        END IF
        
        IF(ALLOCATED(DATASET2)) DEALLOCATE(DATASET2)  
        ALLOCATE(DATASET2(2+NAC+NACD(JWUPU(JS)),1))
        DO I = 1, 2+NAC+NACD(JWUPU(JS))
          CALL XINT(TEMPSET2(I,1), 3, DATASET2(I,1))
        END DO
        CALL HDF5_ADD_DATA(WDO2(HDF5_JFILE,1), 2+NAC+NACD(JWUPU(JS)), 1, DATASET2)
        IF(ALLOCATED(TEMPSET2)) DEALLOCATE(TEMPSET2)
        IF(ALLOCATED(DATASET2)) DEALLOCATE(DATASET2)
      END IF
    END DO

    !-----------------------------------------------------------------------------    
    HDF5_E = HDF5_E + NPU 
    DO JS = 1,NPI 
      IF (LATERAL_PIPE(JS)) THEN  
        JWD = JWD+1  
      ELSE
        JSSS(JBUPI(JS)) = JSSS(JBUPI(JS))+1
      END IF
        
      IF (IWDO(J) == IUPI(JS)) THEN  
        HDF5_JFILE = HDF5_JFILE+1 
        jfile = hdf5_jfile
        IF(ALLOCATED(TEMPSET2)) DEALLOCATE(TEMPSET2)
        ALLOCATE(TEMPSET2(2+NAC+NACD(JWUPI(JS)),1))              
        HDF5_T = HDF5_E + JS 
        TEMPSET2(1,1) = QPI(JS)
        IF (LATERAL_PIPE(JS)) THEN        
          TEMPSET2(2,1) = TAVGW(JWD)                  
          IF (CONSTITUENTS) THEN
            DO JC = 1, NAC
              TEMPSETC2((HDF5_T-1)*NAC+JC,1) = CAVGW(JWD,CN(JC))
              TEMPSET2(2+JC,1) = CAVGW(JWD,CN(JC))                       
            END DO
          END IF
          IF (DERIVED_CALC) THEN
            DO JD = 1, NACD(JWUPI(JS))
              TEMPSETD2((HDF5_T-1)*NACD(JWUPI(JS))+JD,1) = CDAVGW(JWD,CDN(JD,JWUPI(JS)))
              TEMPSET2(2+NAC+JD,1) = CDAVGW(JWD,CDN(JD,JWUPI(JS)))                    
            END DO
          END IF
          !
        ELSE
          TEMPSET2(2,1) = TAVG(JSSS(JBUPI(JS)),JBUPI(JS))           
          IF (CONSTITUENTS) THEN
            DO JC = 1, NAC
              TEMPSETC2((HDF5_T-1)*NAC+JC,1) = CAVG(JSSS(JBUPI(JS)),JBUPI(JS),CN(JC))
              TEMPSET2(2+JC,1) = CAVG(JSSS(JBUPI(JS)),JBUPI(JS),CN(JC))                  
            END DO
          END IF
          IF (DERIVED_CALC) THEN
            DO JD = 1, NACD(JWUPI(JS))
              TEMPSETD2((HDF5_T-1)*NACD(JWUPI(JS))+JD,1) = CDAVG(JSSS(JBUPI(JS)),JBUPI(JS),CDN(JD,JWUPI(JS)))
              TEMPSET2(2+NAC+JD,1) = CDAVG(JSSS(JBUPI(JS)),JBUPI(JS),CDN(JD,JWUPI(JS)))             
            END DO
          END IF
        END IF
        
        IF(ALLOCATED(DATASET2)) DEALLOCATE(DATASET2)  
        ALLOCATE(DATASET2(2+NAC+NACD(JWUPI(JS)),1))
        DO I = 1, 2+NAC+NACD(JWUPI(JS))
          CALL XINT(TEMPSET2(I,1), 3, DATASET2(I,1))
        END DO
        CALL HDF5_ADD_DATA(WDO2(HDF5_JFILE,1), 2+NAC+NACD(JWUPI(JS)), 1, DATASET2)
        IF(ALLOCATED(TEMPSET2)) DEALLOCATE(TEMPSET2)
        IF(ALLOCATED(DATASET2)) DEALLOCATE(DATASET2)
      END IF
    END DO

    !-----------------------------------------------------------------------------    
    HDF5_E = HDF5_E + NPI
    DO JS = 1,NGT  
      IF (LATERAL_GATE(JS)) THEN  
        JWD = JWD+1 
      ELSE
        JSSS(JBUGT(JS)) = JSSS(JBUGT(JS))+1 
      END IF
        
      IF (IWDO(J) == IUGT(JS)) THEN 
        HDF5_JFILE = HDF5_JFILE+1  
        jfile = hdf5_jfile
        IF(ALLOCATED(TEMPSET2)) DEALLOCATE(TEMPSET2)
        ALLOCATE(TEMPSET2(2+NAC+NACD(JWUGT(JS)),1))               
        HDF5_T = HDF5_E + JS
        !HDF5_T = HDF5_E + 1
        TEMPSET2(1,1) = QGT(JS)
        IF (LATERAL_GATE(JS)) THEN             
          TEMPSET2(2,1) = TAVGW(JWD)                   
          IF (CONSTITUENTS) THEN
            DO JC = 1, NAC
              TEMPSETC2((HDF5_T-1)*NAC+JC,1) = CAVGW(JWD,CN(JC))
              TEMPSET2(2+JC,1) = CAVGW(JWD,CN(JC))                      
            END DO
          END IF
          IF (DERIVED_CALC) THEN
            DO JD = 1, NACD(JWUGT(JS))
              TEMPSETD2((HDF5_T-1)*NACD(JWUGT(JS))+JD,1) = CDAVGW(JWD,CDN(JD,JWUGT(JS)))
              TEMPSET2(2+NAC+JD,1) = CDAVGW(JWD,CDN(JD,JWUGT(JS)))                  
            END DO
          END IF
          !
        ELSE
          TEMPSET2(2,1) = TAVG(JSSS(JBUGT(JS)),JBUGT(JS))            
          IF (CONSTITUENTS) THEN
            DO JC = 1, NAC
              TEMPSETC2((HDF5_T-1)*NAC+JC,1) = CAVG(JSSS(JBUGT(JS)),JBUGT(JS),CN(JC))
              TEMPSET2(2+JC,1) = CAVG(JSSS(JBUGT(JS)),JBUGT(JS),CN(JC))
            END DO
          END IF
          IF (DERIVED_CALC) THEN
            DO JD = 1, NACD(JWUGT(JS))
              TEMPSETD2((HDF5_T-1)*NACD(JWUGT(JS))+JD,1) = CDAVG(JSSS(JBUGT(JS)),JBUGT(JS),CDN(JD,JWUGT(JS)))
              TEMPSET2(2+NAC+JD,1) = CDAVG(JSSS(JBUGT(JS)),JBUGT(JS),CDN(JD,JWUGT(JS)))
            END DO
          END IF
        END IF
        
        IF(ALLOCATED(DATASET2)) DEALLOCATE(DATASET2)  
        ALLOCATE(DATASET2(2+NAC+NACD(JWUGT(JS)),1))
        DO I = 1, 2+NAC+NACD(JWUGT(JS))
          CALL XINT(TEMPSET2(I,1), 3, DATASET2(I,1))
        END DO
        CALL HDF5_ADD_DATA(WDO2(HDF5_JFILE,1), 2+NAC+NACD(JWUGT(JS)), 1, DATASET2)
        IF(ALLOCATED(TEMPSET2)) DEALLOCATE(TEMPSET2)
        IF(ALLOCATED(DATASET2)) DEALLOCATE(DATASET2)
      END IF
    END DO
    
    !-----------------------------------------------------------------------------    
    HDF5_E = HDF5_E + NGT 
    HDF5_J = J
    IF(ALLOCATED(DATASET2)) DEALLOCATE(DATASET2)
    ALLOCATE(DATASET2(HDF5_NST+1,1))
    CALL XINT(QWDO(J), 3, DATASET2(1,1))
    
    DO I = 1, HDF5_NST
      CALL XINT(QOUTLET(I), 3, DATASET2(1+I,1))
    END DO
    CALL HDF5_ADD_DATA(WDO(J,1), 1+HDF5_NST, 1, DATASET2)
    CALL HDF5_ADD_DATE(WDO(J,1), JDAY)
    
    IF(ALLOCATED(DATASET2)) DEALLOCATE(DATASET2)
    ALLOCATE(DATASET2(HDF5_NST+1, 1))
    CALL XINT(TWDO(J), 3, DATASET2(1,1)) 
    DO I = 1, HDF5_NST
      CALL XINT(TOUTLET(I), 3, DATASET2(1+I,1))
    END DO
    CALL HDF5_ADD_DATA(WDO(J,2), HDF5_NST+1, 1, DATASET2)
      
    IF (CONSTITUENTS) THEN
      IF(ALLOCATED(DATASET2)) DEALLOCATE(DATASET2)
      ALLOCATE(DATASET2(NAC*(HDF5_NST+1), 1))
      DO I = 1, NAC
        CALL XINT(CWDO(CN(I),J), 3, DATASET2(I,1))
        DO JJ = 1, HDF5_E
          CALL XINT(TEMPSETC2(NAC*(JJ-1)+I,1), 3, DATASET2(NAC+NAC*(JJ-1)+I,1))
        END DO
      END DO
      CALL HDF5_ADD_DATA(WDO(J,3), NAC*(HDF5_NST+1), 1, DATASET2)
      IF(ALLOCATED(DATASET2)) DEALLOCATE(DATASET2)
    END IF
      
    IF (DERIVED_CALC) THEN
      IF(ALLOCATED(DATASET2)) DEALLOCATE(DATASET2)
      ALLOCATE(DATASET2(NACD(JW)*(HDF5_NST+1),1))
      DO I = 1, NACD(JW)
        CALL XINT(CDWDO(CDN(I,JW),J), 3, DATASET2(I,1))
        DO JJ = 1, HDF5_E
          CALL XINT(TEMPSETD2(NACD(JW)*(JJ-1)+I,1), 3, DATASET2(NACD(JW)+NACD(JW)*(JJ-1)+I,1))
        END DO
      END DO             
      CALL HDF5_ADD_DATA(WDO(J,4), NACD(JW)*(HDF5_NST+1), 1, DATASET2)
      IF(ALLOCATED(DATASET2)) DEALLOCATE(DATASET2)
    END IF   
  END SUBROUTINE HDF5_OUTPUTA8
    
  !===========================================================================================================================
  ! BRANCHEAS
  SUBROUTINE HDF5_OUTPUTA9(JDAY)
    USE HDF5MOD
    USE TVDC, ONLY: NAC, QTR, TTR, NACTR, CTR, TRCN, QSUM, PR, ELUH, ELDH, QOUT, CN
    USE MAIN, ONLY: EVBR, CMBRS; 
    USE RSTART, ONLY: VOLEV, VOLSBR, VOLTBR, ESBR, ETBR, CMBRT
    USE SELWC,  ONLY: QSTR
    USE MACROPHYTEC, ONLY: MACMBRS, MACMBRT
    IMPLICIT NONE
    
    REAL       :: JDAY
             
    ! FLOW 
    ALLOCATE(DATASET3(NBR, 6+KMX+NST, 1))
    DATASET3(:,:,:) = -99.0
    DO JB = 1, NBR
      CALL XINT(QSUM(JB), 2,  DATASET3(JB,1,1))
      CALL XINT(EVBR(JB), 3,  DATASET3(JB,2,1)) 
      CALL XINT(-VOLEV(JB), 1,DATASET3(JB,3,1))
      CALL XINT(PR(JB), 6,    DATASET3(JB,4,1))
      CALL XINT(ELUH(JB), 3,  DATASET3(JB,5,1))
      CALL XINT(ELDH(JB), 3,  DATASET3(JB,6,1))
      DO K = 1, KMX
        CALL XINT(QOUT(K,JB), 3, DATASET3(JB,6+K,1))
      END DO
      DO JS = 1, JSS(JB)
        CALL XINT(QSTR(JS,JB), 3, DATASET3(JB,6+KMX+JS,1))
      END DO
    END DO
    CALL HDF5_ADD_DATA(BRANCHFN0+ 1, NBR, 6+KMX+NST, 1, DATASET3, 'Z')
    CALL HDF5_ADD_DATE(BRANCHFN0+ 1, JDAY)
      
    ! VOLUME BALANCE
    IF(ALLOCATED(DATASET3)) DEALLOCATE(DATASET3)
    ALLOCATE(DATASET3(NBR,4,1))
    DO JB = 1, NBR
      CALL XINT(VOLSBR(JB), 3, DATASET3(JB,1,1))
      CALL XINT(VOLTBR(JB), 3, DATASET3(JB,2,1))
      CALL XINT(VOLTBR(JB)-VOLSBR(JB), 3, DATASET3(JB,3,1))
      IF (VOLSBR(JB)/=0.0) THEN
        CALL XINT((VOLTBR(JB)-VOLSBR(JB))/VOLSBR(JB)*100.0, 3, DATASET3(JB,4,1))
      ELSE
        DATASET3(JB,4,1) = 0.0
      END IF
    END DO
    CALL HDF5_ADD_DATA(BRANCHFN0+ 2, NBR, 4, 1, DATASET3, 'Z')
      
    ! ENERGY BALANCE 
    IF(ALLOCATED(DATASET3)) DEALLOCATE(DATASET3)
    ALLOCATE(DATASET3(NBR,4,1))
    DO JB = 1, NBR
      CALL XINT(ESBR(JB)*4.184E3, 3, DATASET3(JB,1,1))
      CALL XINT(ETBR(JB)*4.184E3, 3, DATASET3(JB,2,1))
      CALL XINT((ESBR(JB)-ETBR(JB))*4.184E3, 3, DATASET3(JB,3,1))
      IF (ESBR(JB)/=0.0) THEN
        CALL XINT((ESBR(JB)-ETBR(JB))/ESBR(JB)*100.0, 3, DATASET3(JB,4,1))
      ELSE
        DATASET3(JB,4,1) = 0.0
      END IF
    END DO
    CALL HDF5_ADD_DATA(BRANCHFN0+ 3, NBR, 4, 1, DATASET3, 'Z')
      
    ! MASS BALANCE 
    IF ((NMC+NAC)>0) THEN
      IF(ALLOCATED(DATASET3)) DEALLOCATE(DATASET3)
      ALLOCATE(DATASET3(NBR,4*NAC+4*NMC,1))
      DO JB = 1, NBR
        DO JC = 1, NAC
          CALL XINT(CMBRS(CN(JC),JB), 3, DATASET3(JB,(JC-1)*4+1,1))
          CALL XINT(CMBRT(CN(JC),JB), 3, DATASET3(JB,(JC-1)*4+2,1))
          CALL XINT(CMBRT(CN(JC),JB)-CMBRS(CN(JC),JB), 3, DATASET3(JB,(JC-1)*4+3,1))
          IF (CMBRS(CN(JC),JB)/=0.0) THEN
            CALL XINT((CMBRT(CN(JC),JB)-CMBRS(CN(JC),JB))/(CMBRS(CN(JC),JB)+NONZERO)*100.0, 3, DATASET3(JB,(JC-1)*4+4,1))
          ELSE
            DATASET3(JB,(JC-1)*4+4,1) = 0.0
          END IF
        END DO
          
        DO JM = 1, NMC
          CALL XINT(MACMBRS(JB,JM), 3, DATASET3(JB,4*NAC+(JM-1)*4+1,1))
          CALL XINT(MACMBRT(JB,JM), 3, DATASET3(JB,4*NAC+(JM-1)*4+2,1))
          CALL XINT(MACMBRT(JB,JM)-MACMBRS(JB,JM), 3, DATASET3(JB,4*NAC+(JM-1)*4+3,1))
          IF (MACMBRS(JB,JM)/=0.0) THEN
            CALL XINT((MACMBRT(JB,JM)-MACMBRS(JB,JM))/(MACMBRS(JB,JM)+NONZERO)*100.0, 3, DATASET3(JB,4*NAC+(JM-1)*4+4,1))
          ELSE
            DATASET3(JB, 4*NAC+(JM-1)*4+4, 1) = 0.0
          END IF
        END DO
      END DO
      CALL HDF5_ADD_DATA(BRANCHFN0+4, NBR, 4*NMC+4*NAC, 1, DATASET3, 'Z')
      IF(ALLOCATED(DATASET3)) DEALLOCATE(DATASET3)
    END IF
  END SUBROUTINE HDF5_OUTPUTA9
    
END MODULE OUTPUTHDF5