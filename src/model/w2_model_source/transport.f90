
!***********************************************************************************************************************************
!**                                           S U B R O U T I N E   T R A N S P O R T                                             **
!***********************************************************************************************************************************

SUBROUTINE TRANSPORT
  USE GLOBAL; USE GEOMC; USE TVDC; USE TRANS; USE LOGICC; USE STRUCTURES; USE PREC
  IMPLICIT NONE

! Type declarations

  REAL(R8),     SAVE, ALLOCATABLE, DIMENSION(:)     :: RATD,   CURX1,  CURX2,  CURX3
!  REAL,     SAVE, ALLOCATABLE, DIMENSION(:,:)   :: RATZ,   CURZ1,  CURZ2,  CURZ3
  REAL(R8),     SAVE, ALLOCATABLE, DIMENSION(:,:)   :: SF1X,   SF1Z,ALFA,RATS,CURS1,CURS2,CURS3,ALFAZ,RATSZ,CURS1Z,CURS2Z,CURS3Z
  REAL(R8),     SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: SF12X,  SF13X
  REAL(R8),     SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: SF2X,   SF3X,   SF4X,   SF5X,   SF6X,   SF7X,   SF8X,   SF9X,   SF10X,  SF11X
  REAL(R8),     SAVE, ALLOCATABLE, DIMENSION(:,:,:) :: SF2Z,   SF3Z,   SF4Z,   SF5Z,   SF6Z,   SF7Z,   SF8Z,   SF9Z,   SF10Z
  REAL(R8), SAVE, ALLOCATABLE, DIMENSION(:,:)   :: DX1,    DX2,    DX3
  REAL(R8), SAVE, ALLOCATABLE, DIMENSION(:,:)   :: AD1X,   AD2X,   AD3X,   AD1Z,   AD2Z,   AD3Z

  REAL(R8)                                      :: HT,HM,HB,HMIN,C2X,C3X,C1X,COUR,CART,CALF
  REAL(R8)                                      :: FLUX,FTEMP,CREF,CMAX1,CMIN1,C2Z,C1Z,C3Z,RATZI
  REAL(R8)                                      :: DELC,ADELC,ACURZ,RATDI,DLTDLXR
  REAL(R8)                                      :: DLXT,DLXMIN
  INTEGER                                       :: K

! Allocation declarations

  ALLOCATE (RATD(IMX),       CURX1(IMX),      CURX2(IMX),       CURX3(IMX))
  ALLOCATE (SF1X(KMX,IMX),   SF1Z(KMX,NWB))
!  ALLOCATE (RATZ(KMX,NWB),   CURZ1(KMX,NWB),  CURZ2(KMX,NWB),   CURZ3(KMX,NWB))
  ALLOCATE (DX1(KMX,IMX),    DX2(KMX,IMX),    DX3(KMX,IMX))
  ALLOCATE (AD1X(KMX,IMX),   AD2X(KMX,IMX),   AD3X(KMX,IMX))
  ALLOCATE (AD1Z(KMX,IMX),   AD2Z(KMX,IMX),   AD3Z(KMX,IMX))
  ALLOCATE (RATS(KMX,IMX),   CURS1(KMX,IMX), CURS2(KMX,IMX), CURS3(KMX,IMX), ALFA(KMX,IMX))
  ALLOCATE (RATSZ(KMX,IMX),   CURS1Z(KMX,IMX), CURS2Z(KMX,IMX), CURS3Z(KMX,IMX), ALFAZ(KMX,IMX))
  ALLOCATE (SF2X(KMX,IMX,2), SF3X(KMX,IMX,2), SF4X(KMX,IMX,2),  SF5X(KMX,IMX,2),  SF6X(KMX,IMX,2),  SF7X(KMX,IMX,2))
  ALLOCATE (SF8X(KMX,IMX,2), SF9X(KMX,IMX,2), SF10X(KMX,IMX,2), SF11X(KMX,IMX,2), SF12X(KMX,IMX,2), SF13X(KMX,IMX,2))
  ALLOCATE (SF2Z(KMX,2,NWB), SF3Z(KMX,2,NWB), SF4Z(KMX,2,NWB),  SF5Z(KMX,2,NWB),  SF6Z(KMX,2,NWB),  SF7Z(KMX,2,NWB))
  ALLOCATE (SF8Z(KMX,2,NWB), SF9Z(KMX,2,NWB), SF10Z(KMX,2,NWB))

! Variable initialization

  CT   = 0.0; AT   = 0.0; VT   = 0.0; DT   = 0.0; DX1  = 0.0; DX2  = 0.0; DX3  = 0.0; ADZ  = 0.0; ADX  = 0.0;  AD1X = 0.0
  AD2X = 0.0; AD3X = 0.0; AD1Z = 0.0; AD2Z = 0.0; AD3Z = 0.0
RETURN

!***********************************************************************************************************************************
!**                                        I N T E R P O L A T I O N  M U L T I P L I E R S                                       **
!***********************************************************************************************************************************

ENTRY INTERPOLATION_MULTIPLIERS

! Positive horizontal flows

  DO I=2,IMX-1
    DO K=2,KMX-1
      DLXT = DLX(I-1)
      IF (K > KB(I-1) .OR. INTERNAL_WEIR(K,I)) DLXT = DLX(I)
      DLXMIN       =  DMIN1(DLX(I+1),DLX(I))
      SF1X(K,I)    = (DLX(I+1)+DLX(I))*0.5D0
      SF2X(K,I,1)  =  DLX(I)/(DLX(I)+DLX(I+1))
      !SF3X(K,I,1)  =  DLX(I)**2
      SF3X(K,I,1)  =  DLX(I)*DLX(I)          ! SW 4/20/16 SPEED
      SF4X(K,I,1)  =  DLX(I+1)/(DLX(I)+DLX(I+1))
      SF5X(K,I,1)  =  0.25D0*(DLXT+2.0D0*DLX(I)+DLX(I+1))*(DLXT+DLX(I))
      SF6X(K,I,1)  = -0.25D0*(         DLX(I)+DLX(I+1))*(DLXT+DLX(I))
      SF7X(K,I,1)  =  0.25D0*(DLX(I)+DLX(I+1))*(DLXT+2.0D0*DLX(I)+DLX(I+1))
      SF8X(K,I,1)  =  0.50D0*(         DLX(I)-DLX(I+1))*DLXMIN
      SF9X(K,I,1)  =  0.50D0*(DLXT+2.0D0*DLX(I)-DLX(I+1))*DLXMIN
      SF10X(K,I,1) =  0.50D0*(DLXT+3.0D0*DLX(I))         *DLXMIN
      !SF11X(K,I,1) =  SF8X(K,I,1) /SF5X(K,I,1)/SF1X(K,I)           
      !SF12X(K,I,1) =  SF9X(K,I,1) /SF6X(K,I,1)/SF1X(K,I)
      !SF13X(K,I,1) =  SF10X(K,I,1)/SF7X(K,I,1)/SF1X(K,I)
      SF11X(K,I,1) =  SF8X(K,I,1) /(SF5X(K,I,1)*SF1X(K,I))             ! SW 4/20/16 SPEED
      SF12X(K,I,1) =  SF9X(K,I,1) /(SF6X(K,I,1)*SF1X(K,I))
      SF13X(K,I,1) =  SF10X(K,I,1)/(SF7X(K,I,1)*SF1X(K,I))
    END DO
  END DO

! Negative horizontal flows

  DO I=2,IMX-2
    DO K=2,KMX-1
      DLXT = DLX(I+2)
      IF (K > KB(I+2)) DLXT = DLX(I+1)
      DLXMIN       =  DMIN1(DLX(I),DLX(I+1))
      SF1X(K,I)    = (DLX(I+1)+DLX(I))*0.5D0
      SF2X(K,I,2)  =  DLX(I+1)/(DLX(I)+DLX(I+1))
      !SF3X(K,I,2)  =  DLX(I+1)**2
      SF3X(K,I,2)  =  DLX(I+1)*DLX(I+1)            ! SW 4/20/16 SPEED
      SF4X(K,I,2)  =  DLX(I)/(DLX(I)+DLX(I+1))
      SF5X(K,I,2)  =  0.25D0*(DLX(I)+2.0D0*DLX(I+1)+DLXT)*(DLX(I)+DLX(I+1))
      SF6X(K,I,2)  = -0.25D0*(           DLX(I+1)+DLXT)*(DLX(I)+DLX(I+1))
      SF7X(K,I,2)  =  0.25D0*(DLX(I)+2.0D0*DLX(I+1)+DLXT)*(DLX(I+1)+DLXT)
      SF8X(K,I,2)  = -0.50D0*(       3.0D0*DLX(I+1)+DLXT)*DLXMIN
      SF9X(K,I,2)  =  0.50D0*(DLX(I)-2.0D0*DLX(I+1)-DLXT)*DLXMIN
      SF10X(K,I,2) =  0.50D0*(DLX(I)-DLX(I+1))*DLXMIN
      !SF11X(K,I,2) =  SF8X(K,I,2) /SF5X(K,I,2)/SF1X(K,I)
      !SF12X(K,I,2) =  SF9X(K,I,2) /SF6X(K,I,2)/SF1X(K,I)
      !SF13X(K,I,2) =  SF10X(K,I,2)/SF7X(K,I,2)/SF1X(K,I)
      SF11X(K,I,2) =  SF8X(K,I,2) /(SF5X(K,I,2)*SF1X(K,I))         ! SW 4/20/16 SPEED
      SF12X(K,I,2) =  SF9X(K,I,2) /(SF6X(K,I,2)*SF1X(K,I))
      SF13X(K,I,2) =  SF10X(K,I,2)/(SF7X(K,I,2)*SF1X(K,I))
    END DO
  END DO

! Ultimate multipliers

  IF (ULTIMATE(JW)) THEN
    DO JB=BS(JW),BE(JW)
      DO I=US(JB),DS(JB)    !CONCURRENT(I=US(JB):DS(JB))      !FORALL                                                        !DO I=US(JB),DS(JB)
        RATD(I)  =  DLXR(I-1)/DLXR(I)
        !CURX1(I) =  2.0D0*DLX(I)**2/(DLXR(I)+DLXR(I-1))/DLXR(I-1)   ! code speed improvement SW 4/20/16
        !CURX2(I) = -2.0D0*DLX(I)**2/(DLXR(I)*DLXR(I-1))
        !CURX3(I) =  2.0D0*DLX(I)**2/(DLXR(I)+DLXR(I-1))/DLXR(I)
        CURX1(I) =  2.0D0*DLX(I)*DLX(I)/((DLXR(I)+DLXR(I-1))*DLXR(I-1))
        CURX2(I) = -2.0D0*DLX(I)*DLX(I)/(DLXR(I)*DLXR(I-1))
        CURX3(I) =  2.0D0*DLX(I)*DLX(I)/((DLXR(I)+DLXR(I-1))*DLXR(I))
      END DO
    END DO
  END IF

! Vertical positive flows

  DO K=2,KMX-1
    HT            =  H(K-1,JW)
    HM            =  H(K,JW)
    HB            =  H(K+1,JW)
    HMIN          =  DMIN1(HB,HM)
    SF1Z(K,JW)    = (HB+HM)*0.5D0
    !SF2Z(K,1,JW)  =  HM**2
    SF2Z(K,1,JW)  =  HM*HM       ! SW 4/20/16
    SF3Z(K,1,JW)  =  HM/(HM+HB)
    SF4Z(K,1,JW)  =  HB/(HM+HB)
    SF5Z(K,1,JW)  =  0.25D0*(HT+2.0D0*HM+HB)*(HT+HM)
    SF6Z(K,1,JW)  = -0.25D0*(HM+HB)*(HT+HM)
    SF7Z(K,1,JW)  =  0.25D0*(HM+HB)*(HT+2.0D0*HM+HB)
    SF8Z(K,1,JW)  =  0.50D0*(HM-HB)*HMIN
    SF9Z(K,1,JW)  =  0.50D0*(HT+2.0D0*HM-HB)*HMIN
    SF10Z(K,1,JW) =  0.50D0*(HT+3.0D0*HM)*HMIN
  END DO

! Vertical negative flows

  DO K=2,KMX-2
    HT            =  H(K,JW)
    HM            =  H(K+1,JW)
    HB            =  H(K+2,JW)
    HMIN          =  DMIN1(HT,HM)
    SF1Z(K,JW)    = (HM+HT)*0.5D0
    !SF2Z(K,2,JW)  =  HM**2
    SF2Z(K,2,JW)  =  HM*HM      ! SW 4/20/16
    SF3Z(K,2,JW)  =  HM/(HT+HM)
    SF4Z(K,2,JW)  =  HT/(HT+HM)
    SF5Z(K,2,JW)  =  0.25D0*(HT+2.0D0*HM+HB)*(HT+HM)
    SF6Z(K,2,JW)  = -0.25D0*(HM+HB)*(HT+HM)
    SF7Z(K,2,JW)  =  0.25D0*(HT+2.0D0*HM+HB)*(HM+HB)
    SF8Z(K,2,JW)  = -0.50D0*(3.0D0*HM+HB)*HMIN
    SF9Z(K,2,JW)  =  0.50D0*(HT-2.0D0*HM-HB)*HMIN
    SF10Z(K,2,JW) =  0.50D0*(HT-HM)*HMIN
  END DO

! Ultimate multipliers

  IF (ULTIMATE(JW)) THEN    !also called during UPDATE since surface layer properties change 
    DO K=2,KMX
      RATZ(K,JW)  =  AVH2(K-1,DS(BE(JW)))/AVH2(K,DS(BE(JW)))                                         ! SW 5/20/05
      !CURZ1(K,JW) =  2.0*H(K,JW)**2/(AVH2(K-1,DS(BE(JW)))+AVH2(K,DS(BE(JW))))/AVH2(K-1,DS(BE(JW)))   ! SW 5/20/05
      !CURZ2(K,JW) = -2.0*H(K,JW)**2/(AVH2(K-1,DS(BE(JW)))*AVH2(K,DS(BE(JW))))                        ! SW 5/20/05
      !CURZ3(K,JW) =  2.0*H(K,JW)**2/(AVH2(K-1,DS(BE(JW)))+AVH2(K,DS(BE(JW))))/AVH2(K,DS(BE(JW)))     ! SW 5/20/05
      CURZ1(K,JW) =  2.0*H(K,JW)*H(K,JW)/((AVH2(K-1,DS(BE(JW)))+AVH2(K,DS(BE(JW))))*AVH2(K-1,DS(BE(JW))))  ! SW 5/20/05   SPEED 4/20/16
      CURZ2(K,JW) = -2.0*H(K,JW)*H(K,JW)/(AVH2(K-1,DS(BE(JW)))*AVH2(K,DS(BE(JW))))                         ! SW 5/20/05
      CURZ3(K,JW) =  2.0*H(K,JW)*H(K,JW)/((AVH2(K-1,DS(BE(JW)))+AVH2(K,DS(BE(JW))))*AVH2(K,DS(BE(JW))))    ! SW 5/20/05
    END DO
  END IF
RETURN

!***********************************************************************************************************************************
!**                                          H O R I Z O N T A L  M U L T I P L I E R S                                           **
!***********************************************************************************************************************************

ENTRY HORIZONTAL_MULTIPLIERS1

! Horizontal advection and diffusion multipliers FIRST PASS

  IF (UPWIND(JW)) THEN
    DO I=IU,ID-1
      DO K=KT,KB(I)
        IF (U(K,I) >= 0.0) THEN
          DX2(K,I) = -DX(K,I)/SF1X(K,I)
          DX3(K,I) =  DX(K,I)/SF1X(K,I)
        ELSE
          DX1(K,I) = -DX(K,I)/SF1X(K,I)
          DX2(K,I) =  DX(K,I)/SF1X(K,I)
        END IF
      END DO
    END DO
  ELSE
!!$OMP PARALLEL DO
    DO I=IU,ID-1
      DLTDLXR=DLT/DLXR(I)
      DO K=KT,KB(I)
        COUR = U(K,I)*DLTDLXR
        IF (U(K,I) >= 0.0) THEN
          RATS(K,I)      =  RATD(I)
          CURS1(K,I)     =  CURX1(I)
          CURS2(K,I)     =  CURX2(I)
          CURS3(K,I)     =  CURX3(I)
          DX1(K,I)  =  DX(K,I)*SF11X(K,I,1)
          DX2(K,I)  =  DX(K,I)*SF12X(K,I,1)
          DX3(K,I)  =  DX(K,I)*SF13X(K,I,1)
          ALFA(K,I)      =  2.0D0*(DX(K,I)*DLT/(SF1X(K,I)*SF1X(K,I))-(1.0D0-COUR*COUR)*0.1666667D0)*SF3X(K,I,1)     !/6.0
          AD1X(K,I) = (ALFA(K,I)-COUR*SF8X(K,I,1)*0.5D0)/SF5X(K,I,1)
          AD2X(K,I) =  SF4X(K,I,1)+(ALFA(K,I)-COUR*SF9X(K,I,1) *0.5D0)/SF6X(K,I,1)
          AD3X(K,I) =  SF2X(K,I,1)+(ALFA(K,I)-COUR*SF10X(K,I,1)*0.5D0)/SF7X(K,I,1)
        ELSE
          RATS(K,I)      =  RATD(I+1)
          CURS1(K,I)     =  CURX1(I+1)
          CURS2(K,I)     =  CURX2(I+1)
          CURS3(K,I)     =  CURX3(I+1)
          DX1(K,I)  =  DX(K,I)*SF11X(K,I,2)
          DX2(K,I)  =  DX(K,I)*SF12X(K,I,2)
          DX3(K,I)  =  DX(K,I)*SF13X(K,I,2)
          ALFA(K,I)      =  2.0D0*(DX(K,I)*DLT/(SF1X(K,I)*SF1X(K,I))-(1.0D0-COUR*COUR)*0.1666667D0)*SF3X(K,I,2)    !/6.0
          AD1X(K,I) =  SF2X(K,I,2)+(ALFA(K,I)-COUR*SF8X(K,I,2)*0.5D0)/SF5X(K,I,2)
          AD2X(K,I) =  SF4X(K,I,2)+(ALFA(K,I)-COUR*SF9X(K,I,2)*0.5D0)/SF6X(K,I,2)
          AD3X(K,I) = (ALFA(K,I)-COUR*SF10X(K,I,2)*0.5D0)/SF7X(K,I,2)
        END IF
      END DO
    END DO
!!$OMP END PARALLEL DO
  END IF
RETURN

ENTRY HORIZONTAL_MULTIPLIERS

! Horizontal advection and diffusion multipliers

  IF (UPWIND(JW)) THEN
    DO I=IU,ID-1
      DO K=KT,KB(I)
        IF (U(K,I) >= 0.0) THEN
          C2X      =  COLD(K,I)
          C3X      =  COLD(K,I+1)
          ADX(K,I) = (DX2(K,I)-U(K,I))*C2X+DX3(K,I)*C3X
        ELSE
          C1X      =  COLD(K,I)
          C2X      =  COLD(K,I+1)
          ADX(K,I) =  DX1(K,I)*C1X+(DX2(K,I)-U(K,I))*C2X
        END IF
      END DO
    END DO
  ELSE
!!$OMP PARALLEL DO Private(K,COUR,C1X,C2X,C3X,CART,CALF,RATDI,DELC,ADELC,ACURZ,FLUX,FTEMP,CREF,CMAX1,CMIN1,DLTDLXR,IU,ID)
    DO I=IU,ID-1
    DLTDLXR=DLT/DLXR(I)
      DO K=KT,KB(I)
      COUR = U(K,I)*DLTDLXR
        IF (U(K,I) >= 0.0) THEN
          C1X = COLD(K,I-1)
          C2X = COLD(K,I)
          C3X = COLD(K,I+1)
          IF (U(K,I-1) <= 0.0 .OR. K > KB(I-1) .OR. INTERNAL_WEIR(K,I-1)) C1X = COLD(K,I)
          IF (INTERNAL_WEIR(K,I)) C3X = COLD(K,I)
          CART      =  C3X
          CALF      =  C1X
        ELSE
          C1X = COLD(K,I)
          C2X = COLD(K,I+1)
          C3X = COLD(K,I+2)
          IF (U(K,I+2) >= 0.0 .OR. K > KB(I+2) .OR. I == ID-1 .OR. INTERNAL_WEIR(K,I+1)) C3X = COLD(K,I+1)
          IF (INTERNAL_WEIR(K,I)) THEN
            C2X = COLD(K,I)
            C3X = COLD(K,I)
          END IF
          CART      =  C1X
          CALF      =  C3X
        END IF
        if( .not. ultimate(jw))then
        ADX(K,I) = (DX1(K,I)-U(K,I)*AD1X(K,I))*C1X+(DX2(K,I)-U(K,I)*AD2X(K,I))*C2X+(DX3(K,I)-U(K,I)*AD3X(K,I))*C3X
      !  IF (ULTIMATE(JW)) THEN    ! SW Code speedup 6/16/13
        else
          RATDI = 1.0/RATS(K,I)
          DELC  = RATS(K,I)*C3X+(RATDI-RATS(K,I))*C2X-RATDI*C1X
          DELC  = DSIGN(1.0,U(K,I))*DELC
          ADELC = DABS(DELC)
          ACURZ = DABS(CURS3(K,I)*C3X+CURS2(K,I)*C2X+CURS1(K,I)*C1X)
          IF (ACURZ <= 0.6*ADELC) THEN
            FLUX = AD1X(K,I)*C1X+AD2X(K,I)*C2X+AD3X(K,I)*C3X
          ELSE IF (ACURZ >= ADELC) THEN
            FLUX = C2X
          ELSE IF (ABS(COUR) > 0.0) THEN
            FTEMP = AD1X(K,I)*C1X+AD2X(K,I)*C2X+AD3X(K,I)*C3X
            CREF  = CALF+(C2X-CALF)/ABS(COUR)
            IF (DELC > 0.0) THEN
              CMAX1 = DMIN1(CREF,CART)
              IF (CREF < C2X) CMAX1 = CART
              FLUX = 0.5D0*(C2X+CMAX1)
              IF (FTEMP <= CMAX1 .AND. FTEMP >= C2X) FLUX = FTEMP
            ELSE
              CMIN1 = DMAX1(CREF,CART)
              IF (CREF > C2X) CMIN1 = CART
              IF (FTEMP >= CMIN1 .AND. FTEMP <= C2X) THEN
                FLUX = FTEMP
              ELSE IF (FTEMP > 0.0) THEN
                FLUX = 0.5D0*(C2X+CMIN1)
              ELSE
                FLUX = 0.0D0
              END IF
            END IF
          ELSE
            FLUX = 0.0D0
          END IF
          ADX(K,I) = (DX1(K,I)*C1X+DX2(K,I)*C2X+DX3(K,I)*C3X)-U(K,I)*FLUX
        END IF
      END DO
    END DO
!!$OMP END PARALLEL DO
  END IF
RETURN

!***********************************************************************************************************************************
!**                                            V E R T I C A L  M U L T I P L I E R S                                             **
!***********************************************************************************************************************************

ENTRY VERTICAL_MULTIPLIERS1    ! FIRST PASS

! Vertical advection multipliers

  IF (.NOT.UPWIND(JW)) THEN
!!$OMP PARALLEL DO
    DO I=IU,ID
      DO K=KT,KB(I)-1
        IF (W(K,I) >= 0.0) THEN
          RATSZ(K,I)  = RATZ(K,JW)
          CURS1Z(K,I) = CURZ1(K,JW)
          CURS2Z(K,I) = CURZ2(K,JW)
          CURS3Z(K,I) = CURZ3(K,JW)
          IF (K <= KT+1) THEN
            HT    =  H1(KT,I)
            HM    =  H1(K,I)
            HB    =  H1(K+1,I)
            RATSZ(K,I)  =  AVH1(KT,I)/AVH1(K,I)
            !CURS1Z(K,I) =  2.0D0*HM*HM/(AVH1(KT,I)+AVH1(K,I))/AVH1(KT,I)
            CURS1Z(K,I) =  2.0D0*HM*HM/((AVH1(KT,I)+AVH1(K,I))*AVH1(KT,I))      ! SW 4/20/16 SPEED
            CURS2Z(K,I) = -2.0D0*HM*HM/(AVH1(KT,I)*AVH1(K,I))
            !CURS3Z(K,I) =  2.0D0*HM*HM/(AVH1(KT,I)+AVH1(K,I))/AVH1(K,I)
            CURS3Z(K,I) =  2.0D0*HM*HM/((AVH1(KT,I)+AVH1(K,I))*AVH1(K,I))       ! SW 4/20/16 SPEED
            IF (K == KT) THEN
              HM    =  H1(KT,I)
              RATSZ(K,I)  =  1.0D0
              CURS3Z(K,I) =  1.0D0
              CURS2Z(K,I) = -2.0D0
              CURS1Z(K,I) =  1.0D0
            END IF
            HMIN          =  DMIN1(HB,HM)
            SF1Z(K,JW)    = (HB+HM)*0.5D0
            !SF2Z(K,1,JW)  =  HM**2
            SF2Z(K,1,JW)  =  HM*HM        ! SW 4/20/16
            SF3Z(K,1,JW)  =  HM/(HM+HB)
            SF4Z(K,1,JW)  =  HB/(HM+HB)
            SF5Z(K,1,JW)  =  0.25D0*(HT+2.0D0*HM+HB)*(HT+HM)
            SF6Z(K,1,JW)  = -0.25D0*(HM+HB)*(HT+HM)
            SF7Z(K,1,JW)  =  0.25D0*(HM+HB)*(HT+2.0D0*HM+HB)
            SF8Z(K,1,JW)  =  0.5D0*(HM-HB)*HMIN
            SF9Z(K,1,JW)  =  0.5D0*(HT+2.0D0*HM-HB)*HMIN
            SF10Z(K,1,JW) =  0.5D0*(HT+3.0D0*HM)*HMIN
          END IF
          COUR      =  W(K,I)*DLT/SF1Z(K,JW)
          ALFAZ(K,I) =  2.0D0*(DZQ(K,I)*DLT/(SF1Z(K,JW)*SF1Z(K,JW))-(1.0D0-COUR*COUR)*0.16666667D0)*SF2Z(K,1,JW)     !/6.0
          AD1Z(K,I) = (ALFAZ(K,I)-COUR*SF8Z(K,1,JW)*0.5D0)/SF5Z(K,1,JW)
          AD2Z(K,I) =  SF4Z(K,1,JW)+(ALFAZ(K,I)-COUR*SF9Z(K,1,JW) *0.5D0)/SF6Z(K,1,JW)
          AD3Z(K,I) =  SF3Z(K,1,JW)+(ALFAZ(K,I)-COUR*SF10Z(K,1,JW)*0.5D0)/SF7Z(K,1,JW)
        ELSE
          CURS3Z(K,I) = CURZ3(K+1,JW)
          CURS2Z(K,I) = CURZ2(K+1,JW)
          CURS1Z(K,I) = CURZ1(K+1,JW)
          RATSZ(K,I)  = AVH1(K,I)/AVH1(K+1,I)
          IF (K == KT) THEN
            HT            =  H1(KT,I)
            HM            =  H1(KT+1,I)
            HB            =  H1(KT+2,I)
            HMIN          =  DMIN1(HT,HM)
            RATSZ(K,I)          =  AVH1(KT,I)/AVH1(K,I)
            !CURS1Z(K,I)         =  2.0D0*HM*HM/(AVH1(KT,I)+AVH1(K,I))/AVH1(KT,I)        ! SW 4/20/16 SPEED
            CURS1Z(K,I)         =  2.0D0*HM*HM/((AVH1(KT,I)+AVH1(K,I))*AVH1(KT,I))       
            CURS2Z(K,I)         = -2.0D0*HM*HM/(AVH1(KT,I)*AVH1(K,I))
            !CURS3Z(K,I)         =  2.0D0*HM*HM/(AVH1(KT,I)+AVH1(K,I))/AVH1(K,I)          ! SW 4/20/16 SPEED
            CURS3Z(K,I)         =  2.0D0*HM*HM/((AVH1(KT,I)+AVH1(K,I))*AVH1(K,I))   
            SF1Z(K,JW)    = (HM+HT)*0.5D0
            !SF2Z(K,2,JW)  =  HM**2
            SF2Z(K,2,JW)  =  HM*HM      ! SW 4/20/16
            SF3Z(K,2,JW)  =  HM/(HT+HM)
            SF4Z(K,2,JW)  =  HT/(HT+HM)
            SF5Z(K,2,JW)  =  0.25D0*(HT+2.0D0*HM+HB)*(HT+HM)
            SF6Z(K,2,JW)  = -0.25D0*(HM+HB)*(HT+HM)
            SF7Z(K,2,JW)  =  0.25D0*(HT+2.0D0*HM+HB)*(HM+HB)
            SF8Z(K,2,JW)  = -0.5D0*(3.0D0*HM+HB)*HMIN
            SF9Z(K,2,JW)  =  0.5D0*(HT-2.0D0*HM-HB)*HMIN
            SF10Z(K,2,JW) =  0.5D0*(HT-HM)*HMIN
          END IF
          COUR      =  W(K,I)*DLT/SF1Z(K,JW)
          ALFAZ(K,I) =  2.0D0*(DZQ(K,I)*DLT/(SF1Z(K,JW)*SF1Z(K,JW))-(1.0D0-COUR*COUR)*0.16666667D0)*SF2Z(K,2,JW)    !/6.0
          AD1Z(K,I) =  SF3Z(K,2,JW)+(ALFAZ(K,I)-COUR*SF8Z(K,2,JW)*0.5D0)/SF5Z(K,2,JW)
          AD2Z(K,I) =  SF4Z(K,2,JW)+(ALFAZ(K,I)-COUR*SF9Z(K,2,JW)*0.5D0)/SF6Z(K,2,JW)
          AD3Z(K,I) = (ALFAZ(K,I)-COUR*SF10Z(K,2,JW)*0.5D0)/SF7Z(K,2,JW)
        END IF
      END DO
    END DO
!$END PARALLEL DO
  END IF
RETURN

ENTRY VERTICAL_MULTIPLIERS

! Vertical advection multipliers

  IF (UPWIND(JW)) THEN
    DO I=IU,ID
      DO K=KT,KB(I)-1
        C2Z = COLD(K+1,I)
        IF (W(K,I) >= 0.0) C2Z = COLD(K,I)
        ADZ(K,I) = -W(K,I)*C2Z
      END DO
    END DO
  ELSE
!!$OMP PARALLEL DO
    DO I=IU,ID
      DO K=KT,KB(I)-1
        IF (W(K,I) >= 0.0) THEN
          C1Z   = COLD(K-1,I)
          C2Z   = COLD(K,I)
          C3Z   = COLD(K+1,I)
          CART  = C3Z
          CALF  = C1Z
          IF (K <= KT+1) THEN
            C1Z   =  COLD(KT,I)
            CALF  =  C1Z
          END IF
        ELSE
          C1Z = COLD(K,I)
          C2Z = COLD(K+1,I)
          C3Z = COLD(K+2,I)
          IF (K == KB(I)-1) C3Z = COLD(K+1,I)
          CART  = C1Z
          CALF  = C3Z
        END IF
        if(.not.ultimate(jw))then
        ADZ(K,I) = -W(K,I)*(AD1Z(K,I)*C1Z+AD2Z(K,I)*C2Z+AD3Z(K,I)*C3Z)
        else
    !    IF (ULTIMATE(JW)) THEN    ! SW code speedup 6/16/13
          COUR  =  W(K,I)*DLT/SF1Z(K,JW)
          RATZI = 1.0D0/RATSZ(K,I)
          DELC  = RATSZ(K,I)*C3Z+(RATZI-RATSZ(K,I))*C2Z-RATZI*C1Z
          DELC  = DSIGN(1.0,W(K,I))*DELC
          ADELC = DABS(DELC)
          ACURZ = DABS(CURS3Z(K,I)*C3Z+CURS2Z(K,I)*C2Z+CURS1Z(K,I)*C1Z)
          IF (ACURZ <= 0.6*ADELC) THEN
            FLUX = AD1Z(K,I)*C1Z+AD2Z(K,I)*C2Z+AD3Z(K,I)*C3Z
          ELSE IF (ACURZ >= ADELC) THEN
            FLUX = C2Z
          ELSE IF (DABS(COUR) > 0.0) THEN
            FTEMP = AD1Z(K,I)*C1Z+AD2Z(K,I)*C2Z+AD3Z(K,I)*C3Z
            CREF  = CALF+(C2Z-CALF)/DABS(COUR)
            IF (DELC > 0.0) THEN
              CMAX1 = CART
              IF (CREF >= C2Z) CMAX1 = DMIN1(CREF,CART)
              FLUX = 0.5*(C2Z+CMAX1)
              IF (FTEMP <= CMAX1 .AND. FTEMP >= C2Z) FLUX = FTEMP
            ELSE
              CMIN1 = DMAX1(CREF,CART)
              IF (CREF > C2Z) CMIN1 = CART
              IF (FTEMP >= CMIN1 .AND. FTEMP <= C2Z) THEN
                FLUX = FTEMP
              ELSE IF (FTEMP > 0.0) THEN
                FLUX = 0.5D0*(C2Z+CMIN1)
              ELSE
                FLUX = 0.0D0
              END IF
            END IF
          ELSE
            FLUX = 0.0D0
          END IF
          ADZ(K,I) = -W(K,I)*FLUX
        END IF
      END DO
    END DO
!$END PARALLEL DO
  END IF
RETURN

!***********************************************************************************************************************************
!**                                           H O R I Z O N T A L  T R A N S P O R T                                              **not used since incorporated in Temp and WQ routines directly
!***********************************************************************************************************************************
!
!ENTRY HORIZONTAL_TRANSPORT
!  IF (CONSTITUENTS) THEN
!    DO I=IU,ID
!      DO K=KT,KB(I)     !CONCURRENT(K=KT:KB(I))    !FORALL                                                        !DO K=KT,KB(I)
!        CNEW(K,I) = (COLD(K,I)*BH2(K,I)/DLT+(ADX(K,I)*BHR1(K,I)-ADX(K,I-1)*BHR1(K,I-1))/DLX(I)+(1.0D0-THETA(JW))                     &
!                    *(ADZ(K,I)*BB(K,I)-ADZ(K-1,I)*BB(K-1,I))+SSB(K,I)/DLX(I))*DLT/BH1(K,I)+SSK(K,I)*DLT
!      END DO                                                    
!    END DO                                            
!  ELSE
!    DO I=IU,ID
!      DO K=KT,KB(I)     !CONCURRENT(K=KT:KB(I))      !FORALL                                         !DO K=KT,KB(I)
!        CNEW(K,I) = (COLD(K,I)*BH2(K,I)/DLT+(ADX(K,I)*BHR1(K,I)-ADX(K,I-1)*BHR1(K,I-1))/DLX(I)+(1.0D0-THETA(JW))                     &
!                    *(ADZ(K,I)*BB(K,I)-ADZ(K-1,I)*BB(K-1,I))+SSB(K,I)/DLX(I))*DLT/BH1(K,I)
!      END DO
!    END DO                                           
!  END IF
!RETURN
ENTRY DEALLOCATE_TRANSPORT
  DEALLOCATE (RATD, CURX1, CURX2, CURX3, SF1X,  SF1Z, RATZ, CURZ1, CURZ2, CURZ3, CT,   AT,   VT, DT,RATS,CURS1,CURS2,CURS3,ALFA)
  DEALLOCATE (DX1,  DX2,   DX3,   AD1X,  AD2X,  AD3X,  AD1Z, AD2Z, AD3Z,  SF2X,  SF3X, SF4X, SF5X, SF6X, SF7X)
  DEALLOCATE (RATSZ,CURS1Z,CURS2Z,CURS3Z,ALFAZ)
  DEALLOCATE (SF8X, SF9X,  SF10X, SF11X, SF12X, SF13X, SF2Z, SF3Z, SF4Z,  SF5Z,  SF6Z, SF7Z, SF8Z, SF9Z, SF10Z)
RETURN
END SUBROUTINE TRANSPORT

!***********************************************************************************************************************************
!*                                              S U B R O U T I N E    T R I D I A G                                              **
!***********************************************************************************************************************************

SUBROUTINE TRIDIAG(A,V,C,D,S,E,N,U)
  USE TRIDIAG_V
  INTEGER,                             INTENT(IN)  :: S, E, N
  REAL(R8),              DIMENSION(:), INTENT(IN)  :: A(E),V(E),C(E),D(E)
  REAL(R8),              DIMENSION(:), INTENT(OUT) :: U(N)
 ! REAL(R8), ALLOCATABLE, DIMENSION(:)              :: BTA, GMA
 ! REAL(R8), DIMENSION(1000)              :: BTA, GMA
  INTEGER                                          :: I
!  ALLOCATE (BTA(N),GMA(N))

  BTA1(S) = V(S)
  GMA1(S) = D(S)
  DO I=S+1,E
    BTA1(I) = V(I)-A(I)/BTA1(I-1)*C(I-1)
    GMA1(I) = D(I)-A(I)/BTA1(I-1)*GMA1(I-1)
  END DO
  U(E) = GMA1(E)/BTA1(E)
  DO I=E-1,S,-1
    U(I) = (GMA1(I)-C(I)*U(I+1))/BTA1(I)
  END DO
!  DEALLOCATE (BTA, GMA)                                                                                             ! SW 10/17/05
END SUBROUTINE TRIDIAG
