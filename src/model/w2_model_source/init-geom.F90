!***********************************************************************************************************************************
!**                                                    Task 1.1.3: Geometry                                                       **
!***********************************************************************************************************************************

! Layer elevations
SUBROUTINE INITGEOM
USE MAIN
USE GLOBAL;USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE POROSITYC; USE MACROPHYTEC
  USE RSTART
  IMPLICIT NONE

  INTEGER :: NNBP,NCBP,NINTERNAL,NUP,KTMAX,JJW,IEXIT,KUP,KDN,K1,ICON,ISEG
  REAL    :: ELL1,ELL2,EL1,EL2,B11,ERR1,ERR2,ELR,ELL,ELR2
  CHARACTER(2) :: ICHAR2

  NPOINT = 0
  DO JW=1,NWB
    IF (ZERO_SLOPE(JW)) THEN
      DO I=US(BS(JW))-1,DS(BE(JW))+1
        EL(KMX,I) = ELBOT(JW)
        DO K=KMX-1,1,-1
          EL(K,I) = EL(K+1,I)+H(K,JW)
        END DO
      END DO
    ELSE
      EL(KMX,DS(JBDN(JW))+1) = ELBOT(JW)
      JB                     = JBDN(JW)
      NPOINT(JB)             = 1
      NNBP                   = 1
      NCBP                   = 0
      NINTERNAL              = 0
      NUP                    = 0
      DO WHILE (NNBP <= (BE(JW)-BS(JW)+1))
        NCBP = NCBP+1
        IF (NINTERNAL == 0) THEN
          IF (NUP == 0) THEN
            DO I=DS(JB),US(JB),-1
              IF (I /= DS(JB)) THEN
                EL(KMX,I) = EL(KMX,I+1)+SINA(JB)*(DLX(I)+DLX(I+1))*0.5
              ELSE
                EL(KMX,I) = EL(KMX,I+1)
              END IF
              DO K=KMX-1,1,-1
                EL(K,I) = EL(K+1,I)+H(K,JW)*COSA(JB)
              END DO
            END DO
          ELSE
            DO I=US(JB),DS(JB)
              IF (I /= US(JB)) THEN
                EL(KMX,I) = EL(KMX,I-1)-SINA(JB)*(DLX(I)+DLX(I-1))*0.5
              ELSE
                EL(KMX,I) = EL(KMX,I-1)
              END IF
              DO K=KMX-1,1,-1
                EL(K,I) = EL(K+1,I)+H(K,JW)*COSA(JB)
              END DO
            END DO
            NUP = 0
          END IF
          DO K=KMX,1,-1
            IF (UP_HEAD(JB)) THEN
              EL(K,US(JB)-1) = EL(K,US(JB))+SINA(JB)*DLX(US(JB))
            ELSE
              EL(K,US(JB)-1) = EL(K,US(JB))
            END IF
            IF (DN_HEAD(JB)) THEN
              EL(K,DS(JB)+1) = EL(K,DS(JB))-SINA(JB)*DLX(DS(JB))
            ELSE
              EL(K,DS(JB)+1) = EL(K,DS(JB))
            END IF
          END DO
        ELSE
          DO K=KMX-1,1,-1
            EL(K,UHS(JJB)) = EL(K+1,UHS(JJB))+H(K,JW)*COSA(JB)
          END DO
          DO I=UHS(JJB)+1,DS(JB)
            EL(KMX,I) = EL(KMX,I-1)-SINA(JB)*(DLX(I)+DLX(I-1))*0.5
            DO K=KMX-1,1,-1
              EL(K,I) = EL(K+1,I)+H(K,JW)*COSA(JB)
            END DO
          END DO
          DO I=UHS(JJB)-1,US(JB),-1
            EL(KMX,I) = EL(KMX,I+1)+SINA(JB)*(DLX(I)+DLX(I+1))*0.5
            DO K=KMX-1,1,-1
              EL(K,I) = EL(K+1,I)+H(K,JW)*COSA(JB)
            END DO
          END DO
          NINTERNAL = 0
        END IF
        IF (NNBP == (BE(JW)-BS(JW)+1)) EXIT

!****** Find next branch connected to furthest downstream branch

        DO JB=BS(JW),BE(JW)
          IF (NPOINT(JB) /= 1) THEN
            DO JJB=BS(JW),BE(JW)
              IF (DHS(JB) >= US(JJB) .AND. DHS(JB) <=DS (JJB) .AND. NPOINT(JJB) == 1) THEN
                NPOINT(JB)       = 1
                EL(KMX,DS(JB)+1) = EL(KMX,DHS(JB))+SINA(JB)*(DLX(DS(JB))+DLX(DHS(JB)))*0.5
                NNBP             = NNBP+1; EXIT
              END IF
              IF (UHS(JJB) == DS(JB) .AND. NPOINT(JJB) == 1) THEN
                NPOINT(JB)       = 1
                EL(KMX,DS(JB)+1) = EL(KMX,US(JJB))+(SINA(JJB)*DLX(US(JJB))+SINA(JB)*DLX(DS(JB)))*0.5
                NNBP             = NNBP+1; EXIT
              END IF
              IF (UHS(JJB) >= US(JB) .AND. UHS(JJB) <= DS(JB) .AND. NPOINT(JJB)==1) THEN
                NPOINT(JB)       = 1
                EL(KMX,UHS(JJB)) = EL(KMX,US(JJB))+SINA(JJB)*DLX(US(JJB))*0.5
                NNBP             = NNBP+1
                NINTERNAL        = 1; EXIT
              END IF
              IF (UHS(JB) >= US(JJB) .AND. UHS(JB) <= DS(JJB) .AND. NPOINT(JJB) == 1) THEN
                NPOINT(JB)       = 1
                EL(KMX,US(JB)-1) = EL(KMX,UHS(JB))-SINA(JB)*DLX(US(JB))*0.5
                NNBP             = NNBP+1
                NUP              = 1; EXIT
              END IF
            END DO
            IF (NPOINT(JB)==1) EXIT
          END IF
        END DO
      END DO
    END IF
  END DO

! Minimum/maximum layer heights

  DO JW=1,NWB
    DO K=KMX-1,1,-1
      HMIN = DMIN1(H(K,JW),HMIN)
      HMAX = DMAX1(H(K,JW),HMAX)
    END DO
  END DO
  HMAX2 = HMAX**2

! Water surface and bottom layers

  DO JW=1,NWB
    DO JB=BS(JW),BE(JW)
      DO I=US(JB)-1,DS(JB)+1
        IF (.NOT. RESTART_IN) THEN
          KTI(I) = 2
          DO WHILE (EL(KTI(I),I) > ELWS(I))
            KTI(I) = KTI(I)+1
            if(kti(i) == kmx)then  ! cb 7/7/2010 if elws below grid, setting to elws to just within grid ! SW 5/27/17 DELETED
              kti(i)=kmx-1  ! 2
              elws(i)=el(kti(i),i)
              exit
            end if
          END DO
      
           Z(I)     = (EL(KTI(I),I)-ELWS(I))/COSA(JB)

          
          ZMIN(JW) = DMAX1(ZMIN(JW),Z(I))
          !KTI(I)   =  MAX(KTI(I)-1,2)   ! MOVED SW 5/27/17 
          KTMAX    =  MAX(2,KTI(I))
          KTWB(JW) =  MAX(KTMAX,KTWB(JW))
          KTI(I)   =  MAX(KTI(I)-1,2)    ! original    IF(KTI(I) /= KMX)
          IF (Z(I) > ZMIN(JW)) IZMIN(JW) = I
        END IF
        K = 2
        DO WHILE (B(K,I) > 0.0)
          KB(I) = K
          K     = K+1
        END DO
        KBMAX(JW) = MAX(KBMAX(JW),KB(I))
      END DO
      KB(US(JB)-1) = KB(US(JB))
      KB(DS(JB)+1) = KB(DS(JB))
    END DO


!** Correct for water surface going over several layers

    IF (.NOT. RESTART_IN) THEN
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        DO I=US(JB)-1,DS(JB)+1
          H2(KT,I) = H(KT,JW)-Z(I)
          K        = KTI(I)+1
          DO WHILE (KT > K)
            Z(I)     = Z(I)-H(K,JW)
            H2(KT,I) = H(KT,JW)-Z(I)
            K        = K+1
          END DO
        END DO
      END DO
    END IF
    ELKT(JW) = EL(KTWB(JW),DS(BE(JW)))-Z(DS(BE(JW)))*COSA(BE(JW))
  END DO
  DO JW=1,NWB
    KT = KTWB(JW)
    DO JB=BS(JW),BE(JW)
      IU = US(JB)
      ID = DS(JB)

!**** Boundary bottom layers

      IF (UH_EXTERNAL(JB)) KB(IU-1) = KB(IU)                                        ! CB 6/12/07
      IF (DH_EXTERNAL(JB)) KB(ID+1)  = KB(ID)

!**** Branch numbers corresponding to tributaries, withdrawals, and head

      IF (TRIBUTARIES) THEN
        DO JT=1,NTR
          IF (ITR(JT) >= US(JB) .AND. ITR(JT) <= DS(JB)) JBTR(JT) = JB
        END DO
      END IF
      IF (WITHDRAWALS) THEN
        DO JWD=1,NWD
          IF (IWD(JWD) >= US(JB) .AND. IWD(JWD) <= DS(JB)) JBWD(JWD) = JB
        END DO
      END IF
      IF (UH_INTERNAL(JB)) THEN
        DO JJB=1,NBR
          JBUH(JB) = JJB
          IF (UHS(JB) >= US(JJB) .AND. UHS(JB) <= DS(JJB)) EXIT
        END DO
        DO JJW=1,NWB
          JWUH(JB) = JJW
          IF (JBUH(JB) >= BS(JJW) .AND. JBUH(JB) <= BE(JJW)) EXIT
        END DO
      END IF
      IF (INTERNAL_FLOW(JB)) THEN
        DO JJB=1,NBR
          JBUH(JB) = JJB
          IF (UHS(JB) >= US(JJB) .AND. UHS(JB) <= DS(JJB)) EXIT
        END DO
        DO JJW=1,NWB
          JWUH(JB) = JJW
          IF (JBUH(JB) >= BS(JJW) .AND. JBUH(JB) <= BE(JJW)) EXIT
        END DO
      END IF
      IF (DH_INTERNAL(JB)) THEN
        DO JJB=1,NBR
          JBDH(JB) = JJB
          IF (DHS(JB) >= US(JJB) .AND. DHS(JB) <= DS(JJB)) EXIT
        END DO
        DO JJW=1,NWB
          JWDH(JB) = JJW
          IF (JBDH(JB) >= BS(JJW) .AND. JBDH(JB) <= BE(JJW)) EXIT
        END DO
      END IF

!**** Bottom boundary cells

      IF (UH_INTERNAL(JB)) THEN
        IF (JBUH(JB) >= BS(JW) .AND. JBUH(JB) <= BE(JW)) THEN
          KB(IU-1) = MIN(KB(UHS(JB)),KB(IU))                                 ! CB 6/12/07
        ELSE
          IF (EL(KB(IU),IU) >= EL(KB(UHS(JB)),UHS(JB))) THEN                 ! CB 6/12/07
            KB(IU-1) = KB(IU)                                                ! CB 6/12/07
          ELSE
            DO K=KT,KB(IU)                                                   ! CB 6/12/07
              IF (EL(KB(UHS(JB)),UHS(JB)) >= EL(K,IU)) THEN                  ! CB 6/12/07
                KB(IU-1) = K; EXIT                                           ! CB 6/12/07
              END IF
            END DO
          END IF
        END IF
      END IF
      IF (DH_INTERNAL(JB)) THEN
        IF (JBDH(JB) >= BS(JW) .AND. JBDH(JB) <= BE(JW)) THEN
          KB(ID+1) = MIN(KB(DHS(JB)),KB(ID))
        ELSE
          IF (EL(KB(ID),ID) >= EL(KB(DHS(JB)),DHS(JB))) THEN
            KB(ID+1) = KB(ID)
          ELSE
            DO K=KT,KB(ID)
              IF (EL(KB(DHS(JB)),DHS(JB)) >= EL(K,ID)) THEN
                KB(ID+1) = K; EXIT
              END IF
            END DO
          END IF
        END IF
      END IF

!**** Boundary segment lengths

      DLX(IU-1) = DLX(IU)
      DLX(ID+1) = DLX(ID)

!**** Minimum bottom layers and average segment lengths

      DO I=IU-1,ID
        KBMIN(I) =  MIN(KB(I),KB(I+1))
        DLXR(I)  = (DLX(I)+DLX(I+1))*0.5
      END DO
      KBMIN(ID+1) = KBMIN(ID)
      DLXR(ID+1)  = DLX(ID)

!**** Minimum/maximum segment lengths

      DO I=IU,ID
        DLXMIN = DMIN1(DLXMIN,DLX(I))
        DLXMAX = DMAX1(DLXMAX,DLX(I))
      END DO
    END DO
  END DO

  
! Constrictions                       ! SW 6/26/2018
    CONSTRICTION=.FALSE.   ! INIITALIZE VARIABLES
    BCONSTRICTION=0.0
    INQUIRE(FILE='w2_constriction.csv',EXIST=Constriction(1,1))     !  Just use the first element for file existance
    IF(CONSTRICTION(1,1))THEN
        OPEN(CON,FILE='w2_constriction.csv',STATUS='OLD')
        READ(CON,*)
        READ(CON,*)ICHAR2,ICON    ! NUMBER OF CONSTRICTIONS
        IF(ICHAR2/='ON')THEN
            CONSTRICTION=.FALSE.
        ELSE
        READ(CON,*)
        DO J=1,ICON
            READ(CON,*)ISEG,BCONSTRICTION(ISEG)
            CONSTRICTION(:,ISEG)=.TRUE.
        ENDDO
        ENDIF
    ENDIF
    
  
! Boundary widths

  DO JW=1,NWB
    KT = KTWB(JW)
    DO JB=BS(JW),BE(JW)
      IU = US(JB)
      ID = DS(JB)
      DO I=IU-1,ID+1
        B(1,I) = B(2,I)
        DO K=KB(I)+1,KMX
          B(K,I) = B(KB(I),I)
        END DO
      END DO
    END DO
  END DO

  DO JW=1,NWB
    KT = KTWB(JW)
    DO JB=BS(JW),BE(JW)
      IU    = US(JB)
      ID    = DS(JB)
      IEXIT = 0
      DO K=1,KMX-1
        B(K,IU-1) = B(K,IU)
        IF (UH_INTERNAL(JB) .OR. HEAD_FLOW(JB)) THEN
          IF (JBUH(JB) >= BS(JW) .AND. JBUH(JB) <= BE(JW)) THEN
            B(K,IU-1) = B(K,UHS(JB))
          ELSE
            ELR = EL(K,IU)+SINA(JB)*DLX(IU)*0.5
            ELL = EL(2,UHS(JB))-SINA(JBUH(JB))*DLX(UHS(JB))*0.5
            IF (ELR >= ELL) THEN
              B(K,IU-1) = B(2,UHS(JB))
            ELSE
              DO KUP=2,KMX-1
                ELL1 = EL(KUP,UHS(JB))-SINA(JBUH(JB))*DLX(UHS(JB))*0.5
                ELL2 = EL(KUP+1,UHS(JB))-SINA(JBUH(JB))*DLX(UHS(JB))*0.5
                IF (ELL1 > ELR .AND. ELL2 <= ELR) THEN
                  IF (KUP > KB(UHS(JB))) THEN
                    KB(IU-1)    = K-1
                    KBMIN(IU-1) = MIN(KB(IU),KB(IU-1))
                    IEXIT       = 1
                    EXIT
                  END IF
                  ELR2 = EL(K+1,IU)+SINA(JB)*DLX(IU)*0.5
                  IF (ELR2 >= ELL2) THEN
                    B(K,IU-1) = B(KUP,UHS(JB)); EXIT
                  ELSE
                    K1 = KUP+1
                    IF (K1 > KMX) EXIT
                    B11 = 0.0
                    EL1 = ELR
                    EL2 = EL(K1,UHS(JB))-SINA(JBUH(JB))*DLX(IU)*0.5
                    DO WHILE (ELR2 <= EL2)
                      B11 = B11+(EL1-EL2)*B(K1-1,UHS(JB))
                      EL1 = EL2
                      K1  = K1+1
                      IF (K1 >= KMX+1 .OR. EL2 == ELR2) EXIT
                      EL2 = EL(K1,UHS(JB))-SINA(JBUH(JB))*DLX(UHS(JB))*0.5
                      IF (EL2 <= ELR2) EL2 = ELR2
                    END DO
                    B(K,IU-1) = B11/H(K,JW); EXIT
                  END IF
                END IF
              END DO
              IF (EL(KMX,UHS(JB)) > EL(K,IU)) B(K,IU-1) = B(K-1,IU-1)
              IF (B(K,IU-1) == 0.0) B(K,IU-1) = B(K-1,IU-1)
              IF (IEXIT == 1) EXIT
            END IF
          END IF
        END IF
      END DO
      IEXIT = 0
      DO K=1,KMX-1
        B(K,ID+1) = B(K,ID)
        IF (DH_INTERNAL(JB)) THEN
          IF (JBDH(JB) >= BS(JW) .AND. JBDH(JB) <= BE(JW)) THEN
            B(K,ID+1) = B(K,DHS(JB))
          ELSE
            ELL = EL(K,ID)-SINA(JB)*DLX(ID)*0.5
            ELR = EL(2,DHS(JB))+SINA(JBDH(JB))*DLX(DHS(JB))*0.5
            IF (ELL >= ELR) THEN
              B(K,ID+1) = B(2,DHS(JB))
            ELSE
              DO KDN=2,KMX-1
                ERR1 = EL(KDN,DHS(JB))  +SINA(JBDH(JB))*DLX(DHS(JB))*0.5
                ERR2 = EL(KDN+1,DHS(JB))+SINA(JBDH(JB))*DLX(DHS(JB))*0.5
                IF (ERR1 >= ELL .AND. ERR2 < ELL) THEN
                  IF (KDN > KB(DHS(JB))) THEN
                    KB(ID+1)  = K-1
                    KBMIN(ID) = MIN(KB(ID),KB(ID+1))
                    IEXIT     = 1
                    EXIT
                  END IF
                  ELL2 = EL(K+1,ID)-SINA(JB)*DLX(ID)*0.5
                  IF (ELL2 >= ERR2) THEN
                    B(K,ID+1) = B(KDN,DHS(JB)); EXIT
                  ELSE
                    K1  = KDN+1
                    IF (K1 > KMX) EXIT
                    B11 = 0.0
                    EL2 = ELL
                    EL1 = EL(K1,DHS(JB))+SINA(JBDH(JB))*DLX(DHS(JB))*0.5
                    DO WHILE (ELL2 <= EL1)
                      B11 = B11+(EL2-EL1)*B(K1-1,DHS(JB))
                      EL2 = EL1
                      K1  = K1+1
                      IF (K1 >= KMX+1 .OR. EL1 == ELL2) EXIT
                      EL1 = EL(K1,DHS(JB))+SINA(JBDH(JB))*DLX(DHS(JB))*0.5
                      IF (EL1 <= ELL2) EL1 = ELL2
                    END DO
                    B(K,ID+1) = B11/H(K,JW); EXIT
                  END IF
                END IF
              END DO
              IF (EL(KMX,DHS(JB)) > EL(K,ID)) B(K,ID+1) = B(K-1,ID+1)
              IF (B(K,ID+1) == 0.0) B(K,ID+1) = B(K-1,ID+1)
              IF (IEXIT == 1) EXIT
            END IF
          END IF
        END IF
      END DO
    END DO   ! SW 1/23/06
   END DO    ! SW 1/23/06
   BNEW=B    ! SW 1/23/06
   KBI=KB    ! SW 10/29/2010

!**** Upstream active segment and single layer  ! 1/23/06 entire section moved SW
  DO JW=1,NWB
    KT = KTWB(JW)
    DO JB=BS(JW),BE(JW)
      IU    = US(JB)
      ID    = DS(JB)
      IEXIT = 0
      IF (SLOPE(JB) /= 0.0) THEN
        DO I=US(JB)-1,DS(JB)+1
          IF (KBi(I) < KT ) THEN     ! SW 1/23/06
            do k=kbi(i)+1,kt
            Bnew(K,I) = 0.000001      ! SW 1/23/06
            end do
            KB(I)   = KT
          END IF
        END DO
      END IF
      IUT = IU
      DO I=IU,ID
        IF (KB(I)-KT < NL(JB)-1) IUT = I+1
        ONE_LAYER(I) = KT == KB(I)
      END DO
      DO I=IU-1,ID   ! recompute kbmin after one_layer computation   11/12/07
      KBMIN(I) =  MIN(KB(I),KB(I+1))
      END DO
      KBMIN(ID+1) = KBMIN(ID)

      CUS(JB) = IUT
      IF(IUT>=DS(JB))BR_INACTIVE(JB)= .TRUE.    ! SW 6/12/2017
          

!**** Areas and bottom widths

      IF (.NOT. TRAPEZOIDAL(JW)) THEN                                                                                  !SW 07/16/04
        DO I=IU-1,ID+1
          DO K=1,KMX-1
            BH2(K,I) = B(K,I)*H(K,JW)
            BH(K,I)  = B(K,I)*H(K,JW)
            BB(K,I)  = B(K,I)-(B(K,I)-B(K+1,I))/(0.5*(H(K,JW)+H(K+1,JW)))*H(K,JW)*0.5                                  !SW 08/02/04
          END DO
          BH(KMX,I) = BH(KMX-1,I)
        END DO

! column widths
        DO I=IU-1,ID+1
          CW(KB(I),I)=B(KB(I),I)
          DO K=1,KB(I)-1
              CW(K,I)=B(K,I)-B(K+1,I)
          END DO
        END DO

!****** Derived geometry

        DO I=IU-1,ID+1
          BH2(KT,I) = B(KTI(I),I)*(EL(KT,I)-EL(KTI(I)+1,I)-Z(I)*COSA(JB))/COSA(JB)
          IF (KT == KTI(I)) BH2(KT,I) = H2(KT,I)*B(KT,I)
          DO K=KTI(I)+1,KT
            BH2(KT,I) = BH2(KT,I)+Bnew(K,I)*H(K,JW)                                      ! sw 1/23/06    BH(K,I)
          END DO
          BKT(I)   = BH2(KT,I)/H2(KT,I)
          BI(KT,I) = B(KTI(I),I)
        END DO
      ELSE                                                                                                             !SW 07/16/04
        DO I=IU-1,ID+1
          DO K=1,KMX-1
            BB(K,I)  = B(K,I)-(B(K,I)-B(K+1,I))/(0.5*(H(K,JW)+H(K+1,JW)))*H(K,JW)*0.5
          END DO
          BB(KB(I),I) = B(KB(I),I)*0.5
          BH2(1,I)    = B(1,I)*H(1,JW)
          BH(1,I)     = BH2(1,I)
          DO K=2,KMX-1
            BH2(K,I) = 0.25*H(K,JW)*(BB(K-1,I)+2.*B(K,I)+BB(K,I))
            BH(K,I)  = BH2(K,I)
          END DO
          BH(KMX,I) = BH(KMX-1,I)
        END DO
        DO I=IU-1,ID+1
          CALL GRID_AREA1(EL(KT,I)-Z(I),EL(KT+1,I),BH2(KT,I),BI(KT,I))
          BKT(I) = BH2(KT,I)/H2(KT,I)
        END DO
      END IF
      DO I=IU-1,ID
        DO K=1,KMX-1
          AVH2(K,I) = (H2(K,I) +H2(K+1,I)) *0.5
          AVHR(K,I) =  H2(K,I)+(H2(K,I+1)-H2(K,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                                  !SW 07/29/04
        END DO
        AVH2(KMX,I) = H2(KMX,I)
        DO K=1,KMX
          BR(K,I)   = B(K,I)  +(B(K,I+1)  -B(K,I))  /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                                !SW 07/29/04
          BHR(K,I)  = BH(K,I) +(BH(K,I+1) -BH(K,I)) /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                                !SW 07/29/04
          BHR2(K,I) = BH2(K,I)+(BH2(K,I+1)-BH2(K,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                                !SW 07/29/04
          IF(CONSTRICTION(K,I))THEN    ! SW 6/26/2018
              IF(BR(K,I) > BCONSTRICTION(I))          BR(K,I)  = BCONSTRICTION(I)
              IF(BHR(K,I) > BCONSTRICTION(I)*H(K,JW)) BHR(K,I) = BCONSTRICTION(I)*H(K,JW)
              IF(BHR2(K,I) > BCONSTRICTION(I)*H(K,JW))BHR2(K,I)= BCONSTRICTION(I)*H(K,JW)
          ENDIF
          
        END DO
      END DO
      DO K=1,KMX-1
        AVH2(K,ID+1) = (H2(K,ID+1)+H2(K+1,ID+1))*0.5
        BR(K,ID+1)   =   B(K,ID+1)
        BHR(K,ID+1)  =   BH(K,ID+1)
      END DO
      AVH2(KMX,ID+1) = H2(KMX,ID+1)
      AVHR(KT,ID+1)  = H2(KT,ID+1)
      BHR2(KT,ID+1)  = BH2(KT,ID+1)
      IUT = IU
      IF (UP_HEAD(JB)) IUT = IU-1
      DO I=IUT,ID
        DO K=1,KMX-1
          VOL(K,I) = B(K,I)*H2(K,I)*DLX(I)
        END DO
        VOL(KT,I)    = BH2(KT,I)*DLX(I)
        DEPTHB(KT,I) = H2(KT,I)
        DEPTHM(KT,I) = H2(KT,I)*0.5
        DO K=KT+1,KMX
          DEPTHB(K,I) = DEPTHB(K-1,I)+ H2(K,I)
          DEPTHM(K,I) = DEPTHM(K-1,I)+(H2(K-1,I)+H2(K,I))*0.5
        END DO
      END DO
    END DO
  END DO
  H1   = H2
  BH1  = BH2
  BHR1 = BHR2
  AVH1 = AVH2

! Temporary downstream head segment

  DO JB=1,NBR
   IF (DHS(JB).GT.0) THEN
     DO JJB=1,NBR
       IF (DHS(JB) >= US(JJB) .AND. DHS(JB) <= DS(JJB)) EXIT
     END DO
     IF (CUS(JJB) > DHS(JB)) CDHS(JB) = CUS(JJB)
   END IF
  END DO

! Total active cells

  DO JW=1,NWB
    DO JB=BS(JW),BE(JW)
      IF(BR_INACTIVE(JB))CYCLE   ! SW 6/12/2017
      DO I=CUS(JB),DS(JB)
        DO K=KTWB(JW),KB(I)
          NTAC = NTAC+1
        END DO
      END DO
      NTACMX = NTAC
      NTACMN = NTAC

!**** Wind fetch lengths

      DO I=US(JB),DS(JB)
        FETCHD(I,JB) = FETCHD(I-1,JB)+DLX(I)
      END DO
      DO I=DS(JB),US(JB),-1
        FETCHU(I,JB) = FETCHU(I+1,JB)+DLX(I)
      END DO
    END DO
  END DO

! Segment heights

  DO JW=1,NWB
    DO JB=BS(JW),BE(JW)
      DO I=US(JB)-1,DS(JB)+1
        DO K=MIN(KMX-1,KB(I)),2,-1
          HSEG(K,I) = HSEG(K+1,I)+H2(K,I)
        END DO
      END DO
    END DO
  END DO

! Beginning and ending segments/layers for snapshots

  DO JW=1,NWB
    DO I=1,NISNP(JW)
      KBR(JW) = MAX(KB(ISNP(I,JW)),KBR(JW))
    END DO
  END DO

  RETURN
  END SUBROUTINE INITGEOM
