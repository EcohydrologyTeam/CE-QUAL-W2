! System Total Dissolved Gas (SYSTDG)
! Spillbay TDG production is calculated based on the regression equations included in SYSTDG
! SYSTDG equations were provided by Tammy Threadgill and Ron Thompson.
! Sponsored by USACE Portland District through the CRSO EIS   
!  
! Three new ON/OFF controls are included in this version's control file
! 1) SYSTDG      - calculate Spillway TDG production (1 - 5 equations) and mixing downstream TDG
! 2) N2BDN/DOBND - specify saturation (%) for all inflow N2 and DO boundary conditions
! 3) TDGTA       - specify downstream or spillflow TDG target to reallocate spillflow into powerhouse to meet TDG targets
!
! If SYSTDG = ON, a SYSTDG control input file (w2_systdg.npt) is required.
! If TDGTA  = ON, a TDG target control input file (w2_tdgtarget.npt) is required.
!
! Developed from PSU W2V4.1 dated Sept 2, 2018 
! Drs Zhonglong Zhang (PSU) and Hongda Wang (UCD); Revised by SWells Sept 2019 for inclusion in updated 4.2 code; August 2020 csv file format added SW
!  
! 
MODULE modSYSTDG
 USE PREC,  ONLY: R8
 USE TDGAS; USE STRUCTURES; USE GLOBAL; USE MAIN, ONLY: EA, GTTYP, GTPC, Q, WBSEG, TEXT, ERROR_OPEN, SYSTDGC, N2BNDC, DOBNDC, TDGTAC,TMSTRT, CONTDG, TDG2BNDC
 USE TVDC, ONLY: TAIR, TDEW; USE KINETIC, ONLY: TDG
  IMPLICIT NONE
  REAL(R8), ALLOCATABLE, DIMENSION(:) :: TDG_PHS, TDG_FLS, TDG_TDP 
  REAL(R8), ALLOCATABLE, DIMENSION(:) :: BAYC, QBAY
  REAL(R8)                            :: qs, TDG_ROSP, TDG_TDG, TDG_REL
  REAL(R8)                            :: TDGP1, TDGP2, TDGP3, TDGP4, TDGP12, TDGP22, TDGP32, TDGP42, TDGE1, TDGE2, TDGE12, TDGE22, ROP1, ROP2, ROP3, ROP4
  REAL(R8)                            :: TWE, TWCE, TWE_TS, FBE, QSPILL, TDGSPMN
  INTEGER, ALLOCATABLE, DIMENSION(:)  :: POWGTNO, FLGTNO
  INTEGER                             :: POWNO, FLNO, BEGNO, ENDNO, NRO, NBAY
  INTEGER                             :: TWEMOD, TDGEQ, TDGROEQ, TDGENTEQ
  INTEGER                             :: ig, ip, ifl, ib
  LOGICAL,  ALLOCATABLE, DIMENSION(:) :: GTNAME
  CHARACTER(72)                       :: TITLESYSTDG(10)
  CHARACTER(8)                        :: TWETSC, TDGLOC
  CHARACTER(72)                       :: TWEFN 
  REAL                                :: NXTSPLIT3
  !
  CONTAINS
  !===========================================================================================================================
  SUBROUTINE INPUT_SYSTDG
  !USE MAIN, ONLY: TMSTRT, CONTDG       !, GTTYP, GTPC, SYSTDGC, N2BNDC, DOBNDC, TDGTAC
   IMPLICIT NONE
   CHARACTER(8) :: AID1
   LOGICAL      :: CSVFORMAT
   INTEGER      :: I, IG, N_POW, N_FLD, N_SPB
      NXTSPLIT3=TMSTRT
      open  (88888, FILE='TDG_output.csv', status='unknown')
      WRITE (88888, '(A, <NGT>("QGT-",I2,","))')'JDAY,TDG_TDG,SUM_QGT2,',(IG, IG = 1, NGT)
      NRO=0                                                 
      POWNO = 0
      FLNO  = 0
      NBAY  = 0
      
     CSVFORMAT=.FALSE.   
     READ (CONTDG,'(A)')TITLESYSTDG(1)
     IF(TITLESYSTDG(1)(1:1)=='$')CSVFORMAT=.TRUE.
     
     IF(.NOT.CSVFORMAT)THEN
     READ (CONTDG,'(A)')TITLESYSTDG(1)    ! READ NEXT LINE - IF COMMAS IN FIRST FEW FIELDS IT IS IN CSV FORMAT
     DO I=1,7
         IF(TITLESYSTDG(1)(I:I)==',')THEN
             CSVFORMAT=.TRUE.
             EXIT
         ENDIF
     ENDDO
     REWIND(CONTDG)
     ENDIF
      
      IF(CSVFORMAT)THEN
      READ (CONTDG,*)
      READ (CONTDG,*)
      READ (CONTDG,*)
      DO I=1,10
      READ (CONTDG,*) AID1,TITLESYSTDG(I)
      ENDDO
      READ (CONTDG,*)
      READ (CONTDG,*)
      READ (CONTDG,*)AID1,SYSTDGC, N2BNDC, DOBNDC,  TDG2BNDC,TDGTAC;SYSTDGC=ADJUSTR(SYSTDGC); N2BNDC=ADJUSTR(N2BNDC); DOBNDC=ADJUSTR(DOBNDC); TDGTAC=ADJUSTR(TDGTAC)
      TDG2BNDC=ADJUSTR(TDG2BNDC)
      READ (CONTDG,*)
      READ (CONTDG,*)
      DO IG=1,NGT
      READ (CONTDG,*) AID1,GTTYP(IG), GTPC(IG)
      ENDDO
      GTTYP=ADJUSTR(GTTYP)

      DO IG = 1, NGT
         IF (GTTYP(IG)=='     POW') POWNO = POWNO + 1
         IF (GTTYP(IG)=='     FLD') FLNO  = FLNO  + 1
         IF (GTTYP(IG)=='      RO')  NRO   = NRO  + 1
         IF (GTTYP(IG)=='     SPB') NBAY  = NBAY  + 1
      END DO       
      READ (CONTDG,*)
      READ (CONTDG,*)

      READ (CONTDG,*) AID1,FBE, TWCE, TWEMOD, TWE, TWETSC, TDGLOC, QSPILL, TDGSPMN;  TWETSC=ADJUSTR(TWETSC); TDGLOC=ADJUSTR(TDGLOC)
      READ (CONTDG,*)
      READ (CONTDG,*)
      READ (CONTDG,*) AID1,TDGEQ, TDGP1, TDGP2, TDGP3, TDGP4, TDGP12, TDGP22, TDGP32, TDGP42   
      IF (NRO>0) READ (CONTDG,*) AID1,TDGROEQ, ROP1, ROP2, ROP3, ROP4     
      READ (CONTDG,*)
      READ (CONTDG,*)

      READ (CONTDG,*) AID1,TDGENTEQ, TDGE1, TDGE2, TDGE12, TDGE22   
      READ (CONTDG,*)
      READ (CONTDG,*)

      READ (CONTDG,*) AID1, TWEFN                                                             
      CLOSE(CONTDG)       
ELSE
      READ (CONTDG, '(///(8X,A72))') (TITLESYSTDG(i), i=1,10)
      READ (CONTDG, '(//8x,5A8)')SYSTDGC, N2BNDC, DOBNDC, TDG2BNDC,TDGTAC
      READ (CONTDG,'(//(:8X,A8,F8.2))') (GTTYP(IG), GTPC(IG), IG=1,NGT) 

      DO IG = 1, NGT
         IF (GTTYP(IG)=='     POW') POWNO = POWNO + 1
         IF (GTTYP(IG)=='     FLD') FLNO  = FLNO  + 1
         IF (GTTYP(IG)=='      RO')  NRO   = NRO  + 1
         IF (GTTYP(IG)=='     SPB') NBAY  = NBAY  + 1
      END DO                                                                                             

      READ (CONTDG,'(//8X,2F8.3,I8,F8.3,2A8,2F8.3)') FBE, TWCE, TWEMOD, TWE, TWETSC, TDGLOC, QSPILL, TDGSPMN      
      READ (CONTDG,'(//8X,I8,8F8.3)') TDGEQ, TDGP1, TDGP2, TDGP3, TDGP4, TDGP12, TDGP22, TDGP32, TDGP42    
      IF (NRO>0) READ (CONTDG,'(8X,I8,4F8.5)') TDGROEQ, ROP1, ROP2, ROP3, ROP4                             
      READ (CONTDG,'(//8X,I8,4F8.3)') TDGENTEQ, TDGE1, TDGE2, TDGE12, TDGE22                               
      READ (CONTDG,'(//(8X,A72))')  TWEFN                                                             
      CLOSE(CONTDG)
ENDIF
IF(SYSTDGC == '     OFF')GO TO 100          ! DO NOT ALLOCATE ARRAYS IF WE ARE NOT USING SYSTDG
      IF (POWNO>0) THEN 
          ALLOCATE (POWGTNO(POWNO), TDG_PHS(POWNO))
      ELSE
          ALLOCATE (POWGTNO(1), TDG_PHS(1))
          POWGTNO(1)=0
      END IF
      IF (FLNO>0) THEN 
          ALLOCATE (FLGTNO(FLNO), TDG_FLS(FLNO))
      ELSE
         ALLOCATE (FLGTNO(1), TDG_FLS(1))
         FLGTNO(1)=0
      END IF
      N_POW = 0
      N_FLD = 0
      DO ig = 1, NGT
         IF (GTTYP(ig)=='     POW') THEN
            N_POW = N_POW + 1
            POWGTNO(N_POW)=ig
         END IF
         IF (GTTYP(ig)=='     FLD') THEN
            N_FLD = N_FLD + 1
            FLGTNO(N_FLD)=ig
         END IF
      END DO
      DO ig = 1, NGT
         IF (GTTYP(ig)=='     SPB') THEN
            BEGNO = ig
            EXIT
         END IF
      END DO
      ENDNO = BEGNO + NBAY -1
      IF (NBAY/=0) THEN
         ALLOCATE (BAYC(NBAY), QBAY(NBAY))
      ELSE
        ALLOCATE (BAYC(1), QBAY(1))
      END IF
      N_SPB = 0
      DO ig = 1, NGT
        IF (GTTYP(ig)=='     SPB') THEN
           N_SPB = N_SPB + 1
           BAYC(N_SPB)=GTPC(ig)
        END IF
      END DO
      !
      ALLOCATE (GTNAME(NGT), TDG_TDP(NGT))
      GTNAME=.FALSE.
      DO ig=1, NGT
          IF (GTTYP(ig)=='      RO') GTNAME(ig)=.TRUE.          ! C IN RO IS UPDATE FROM TDG
          IF (ig>=BEGNO .AND. ig<=ENDNO) GTNAME(ig)=.TRUE.      ! GTNAME = GATE IS A SPILLBAY== .TRUE.
      END DO
100      RETURN
  END SUBROUTINE
  
  !===========================================================================================================================
  ! allocate and initialize all input parameter
  SUBROUTINE SYSTDG_qs
    IMPLICIT NONE
    REAL(R8)          :: SUMQ, SUMQ1
    SUMQ = 0.0
    SUMQ1 = 0.0
    DO ib = 1, NBAY
      SUMQ = SUMQ + QBAY(ib)**BAYC(ib)
      SUMQ1 = SUMQ1 + QBAY(ib)**(BAYC(ib) - 1.0)
    END DO
    IF (SUMQ1 /= 0.0) qs = SUMQ / SUMQ1 *35.3147/1000.0    ! CMS TO KCFS 
  END SUBROUTINE SYSTDG_qs
  
  !===========================================================================================================================
  ! TDG production calculation in SYSTDG
  SUBROUTINE UPDATE_TDGC (NSAT, P, N, T, TDGC)
      IMPLICIT NONE
      INTEGER          :: NSAT, N
      REAL(R8)         :: P, T, TDGC
      REAL(R8)         :: SAT
      ! CALCULATE SAT
      IF (NSAT == 0) THEN   ! O2 saturation
        SAT = EXP(7.7117 - 1.31403 * (LOG(T + 45.93))) * P
      ELSE IF (NSAT == 1) THEN            ! N2 saturation
        EA = DEXP(2.3026D0 * (7.5D0 * TDEW(WBSEG(IUGT(N))) / (TDEW(WBSEG(IUGT(N))) + 237.3D0) + 0.6609D0)) * 0.001316   ! mmHg     
        SAT = (1.5568D06 * 0.79 * (P - EA) * (1.8816D-5 - 4.116D-7 * T + 4.6D-9 * T*T)) 
      ELSE
        SAT = P  
      END IF
      TDGC = TDG_TDP(N) * SAT / 100.0
      IF (POWNO>0 .AND. TDGLOC=='     REL') THEN
         TDGC = TDG_TDG * SAT / 100.0
      END IF
    END SUBROUTINE UPDATE_TDGC
    
    SUBROUTINE SYSTDG_TDG
    USE SCREENC, ONLY:JDAY; 
      IMPLICIT NONE
      REAL(R8)         :: P1, P2, P3, P4, E1, E2
      REAL(R8)         :: qs_RO, W2FBE, Q_SUM, TEMP_TW, SUM_TDGPHK
      REAL(R8)         :: Q_ROSP, QRO, QSP, QPH, QTOT, TDG_QROSP, TDG_QPH, TDG_QTOT, TDG_QENT
      REAL(R8)         :: SUM_TDG_ROS, SUM_TDG_SPS, SUM_TDG_PHS
      REAL(R8)         :: TDG_RO, TDG_SP, TDG_PH
      INTEGER*8        :: SUM_K, IK
      REAL             :: SUM_QGT2
      Q_SUM = 0.0
      ! ADD QRO AND QSP
      DO ig = 1, NGT
         IF (GTTYP(ig)=='      RO') THEN
            Q_SUM = Q_SUM + QGT(ig)
         ELSE
            IF (ig>=BEGNO.AND. ig<=ENDNO)  Q_SUM = Q_SUM + QGT(ig)
         END IF
      END DO
      TDG_ROSP = 0.0
      IF (Q_SUM/=0.0) THEN
          TDG_TDP(:)  = 0.0                                                ! INITIAL TDG FOR GATES
          !
          SUM_TDG_ROS = 0.0                                                ! INITIAL SUM TDG FOR RO
          QRO         = 0.0                                                ! INITIAL SUM Q FOR RO
          DO ig = 1, NGT
             IF (GTTYP(ig)=='      RO') THEN                               
                 QRO = QRO + QGT(ig)
             END IF 
          END DO
          qs_RO = QRO * 35.3147/1000.0                                     ! CMS TO KCFS 
          DO ig = 1, NGT
             IF (GTTYP(ig)=='      RO' .AND. QGT(ig)>0.0) THEN
                TDG_TDP(ig)= (ROP1 * (1 - EXP(ROP3 * qs_RO)) + PALT(IUGT(ig))*760.0)/(PALT(IUGT(ig))*760.0) *100.0     ! RO USE EQ 1     ! (mmgh) TO TDG (%)
                IF (TDG_TDP(ig) > 145.0) TDG_TDP(ig) = 145.0               ! TDG <= 145.0
                SUM_TDG_ROS = SUM_TDG_ROS + TDG_TDP(ig)*QGT(ig)            ! SUM TDG FOR RO
             END IF 
          END DO
          IF (QRO /= 0.0) THEN
             TDG_RO = SUM_TDG_ROS/QRO                                      ! FLOW AVERAGED TDG_RO (%)
          ELSE
             TDG_RO = 0.0                                                  ! ZERO FLOW --> RO TDG IS 0.0
          END IF
          ! SET E1 E2 INCASE ONLY RO FLOW TO CALCULATE TDG_REL
          IF (POWNO>0 .AND. TDGLOC=='     REL') THEN
             DO ip = 1, POWNO
                IF (QGT(POWGTNO(ip))/=0.0) THEN
                   W2FBE=Q(IUGT(ip)-US(JBUGT(ip))+1) ! FLOW DISCHARGE BEFORE DAM
                   EXIT
                END IF
             END DO
             IF (FBE < 0.0) THEN
                E1 = TDGE1
                E2 = TDGE2
             ELSE IF (W2FBE < FBE) THEN
                E1 = TDGE1
                E2 = TDGE2
             ELSE IF (W2FBE>= FBE) THEN
                E1 = TDGE12
                E2 = TDGE22
             END IF
          END IF
          !
          IF (NBAY/=0) QBAY(1:NBAY)=QGT(BEGNO:ENDNO)                       ! BAY FLOW QBAY FROM GATE FLOW QGT
          CALL SYSTDG_qs                                                   ! CALCULATE qs                              
          ! TWE Recalculation
          IF (TWETSC == '      ON') TWE=TWE_TS                             ! UPDATE TWE TO TWE_TS ACCORDING TO CONTROL VARIABLE TWETSC 
          IF (TWEMOD == 1) TWE=TWE * 0.934 + 4.94                          ! UPDATE TWE ACCORDING TO TWEMOD
          SUM_TDG_SPS = 0.0                                                ! INITIAL TDG*QSP FOR SPILL
          QSP         = 0.0                                                ! INITIAL SUM Q FOR SPILL
          TDG_SP      = 0.0                                                ! INITIAL TDG FOR SPILL
          ! CALCULATE TDG_TDP FOR EACH BAY FROM BEGNO TO ENDNO
          DO ig = BEGNO, ENDNO
             IF (QGT(ig)/=0.0) THEN
                ! READ W2FBE
                IF (IUGT(ig)>=US(JBUGT(ig))+1) THEN                        ! GATE IS NOT LOCATED IN THE 1ST SEGMENT                 
                    W2FBE=Q(IUGT(ig)-US(JBUGT(ig))+1)                      ! FLOW DISCHARGE BEFORE DAM
                ELSE 
                    W2FBE=0.0
                END IF
                IF (FBE < 0.0) THEN
                   P1 = TDGP1
                   P2 = TDGP2
                   P3 = TDGP3
                   P4 = TDGP4
                   E1 = TDGE1
                   E2 = TDGE2
                ELSE IF (W2FBE < FBE) THEN
                   P1 = TDGP1
                   P2 = TDGP2
                   P3 = TDGP3
                   P4 = TDGP4
                   E1 = TDGE1
                   E2 = TDGE2
                ELSE IF (W2FBE>= FBE) THEN
                   P1 = TDGP12
                   P2 = TDGP22
                   P3 = TDGP32
                   P4 = TDGP42
                   E1 = TDGE12
                   E2 = TDGE22
                END IF
                ! TDG EQUATIONS 1/2/3/4/5
                IF (TDGEQ == 2) THEN
                   TDG_TDP(ig) = P1 *((TWE - TWCE)**P2 ) * (1.0 - EXP(P3 * qs)) + P4 + PALT(IUGT(ig))*760.0  ! mmHg
                !
                ELSE IF (TDGEQ == 3) THEN
                   TDG_TDP(ig) = P1* ((TWE - TWCE)**P2) * (qs**P3) + P4 + PALT(IUGT(ig))*760.0               ! mmHg
                !
                ELSE IF (TDGEQ == 4) THEN
                   TDG_TDP(ig) = P1 * (TWE - TWCE) + P2 *(qs**P3) + P4 + PALT(IUGT(ig))*760.0                ! mmHg
                !
                ELSE IF (TDGEQ == 5) THEN                   
                   SUM_K = 0
                   TEMP_TW = 0.0
                   DO IK = KT, KB(IUGT(ig)+1)
                      IF (T2(IK,IUGT(ig)+1)/=-99.0 .AND. T2(IK,IUGT(ig)+1)>0.0) THEN                         ! WATER TEMP IS >0.0
                         SUM_K=SUM_K+1
                         TEMP_TW = TEMP_TW + T2(IK,IUGT(ig)+1)                                               ! SUM EVERY LAYER OF TAIL WATER
                      END IF
                   END DO
                IF (SUM_K/=0) TEMP_TW = TEMP_TW/SUM_K                                                        ! AVERAGED TEMP OF TAIL WATER
                   TDG_TDP(ig) = P1 * (1.0 - EXP(P2 * qs)) + P3 * (TEMP_TW - P4) + PALT(IUGT(ig))*760.0      ! mmHg
                !
                ELSE IF (TDGEQ == 1) THEN
                   TDG_TDP(ig) = P1 * (1.0 - EXP(P3 * qs)) + PALT(IUGT(ig))*760.0 ! mmHg
                !
                ELSE 
                TEXT='TDGEQ INPUT ERROR'
                ERROR_OPEN=.TRUE. 
                !
                END IF
                TDG_TDP(ig)= TDG_TDP(ig)/(PALT(IUGT(ig))*760.0) *100.0                                       ! (mmgh) TO TDG (%)
                IF (TDG_TDP(ig) > 145.0) TDG_TDP(ig) = 145.0                                                 ! TDG <= 145.0
                SUM_TDG_SPS = SUM_TDG_SPS + TDG_TDP(ig) * QGT(ig)                                            ! SUM TDG_SP*QGT 
                QSP = QSP + QGT(ig)                                                                          ! SUM Q OF SPILL
                ! If total spill <= 50 kcfs, TDG % saturation = 110 % in the spillway outlets.  
                IF ((QSP*35.3147/1000.0)>0.0 .and. (QSP*35.3147/1000.0)<=QSPILL) TDG_TDP(ig) =  TDGSPMN      
             END IF
          END DO
          IF (QSP/=0.0) TDG_SP= SUM_TDG_SPS/QSP                                                              ! FLOW AVERAGED TDG FOR SPILL
          Q_ROSP = QRO + QSP                                                                                 ! ADD QSP INTO Q_ROSP
          TDG_ROSP = (SUM_TDG_ROS+SUM_TDG_SPS)/Q_ROSP                                                        ! FLOW AVERAGED TDG FOR SPILL WITH RO
          TDG_TDG = TDG_ROSP                                                                                 ! OUTPUT TDG = TDG_SP
          ! If total spill <= 50 kcfs, TDG % saturation = 110 % in the spillway outlets.  
          IF ((QSP*35.3147/1000.0)>0.0 .and. (QSP*35.3147/1000.0)<=QSPILL) TDG_TDG =  TDGSPMN
          !
          ! QENT CALCULATIONS UPDATE TDG_TDG TO TDG_REL
          IF (POWNO>0 .AND. TDGLOC=='     REL') THEN
             ! TDG POWER HOUSE(S) AND Q POWER HOUSE(S)
             QPH  = 0.0                                                                                       ! INITIAL Q OF POWER HOUSE(S)
             TDG_PHS(:)   = 0.0                                                                               ! INITIAL TDG OF POWER HOUSE(S)                               
             SUM_TDG_PHS  = 0.0                                                                               ! INITIAL SUM TDG*QPH
             DO ip = 1, POWNO
                IF (QGT(POWGTNO(ip))/=0.0) THEN
                    QPH = QPH + QGT(POWGTNO(ip))                                                                 ! Q POWER HOUSE(S)
                    SUM_K=0                                                                                      ! INITIAL SUM K
                    SUM_TDGPHK=0.0                                                                               ! INITIAL SUM TDG OF POWER HOUSE(S) SEGMENT
                    DO IK=1, KMX
                        IF (TDG(IK, IUGT(POWGTNO(ip)))>0.0) THEN                                                 ! TDG > 0.0 LAYER 
                            SUM_K=SUM_K+1                                                                        ! SUM LAYER COUNT
                            SUM_TDGPHK=SUM_TDGPHK+ TDG(IK, IUGT(POWGTNO(ip)))                                    ! SUM TDG OF POWER HOUSE(S)
                        END IF
                    END DO
                    IF (SUM_K/=0) THEN
                        TDG_PHS(ip)=SUM_TDGPHK/SUM_K                                                             ! LAYER AVERAGED TDG FOR POWER HOUSE(S)
                    ELSE
                        TDG_PHS(ip)=0.0
                    END IF
                    SUM_TDG_PHS = SUM_TDG_PHS + TDG_PHS(ip)*QGT(POWGTNO(ip))                                     ! SUM TDG*QPH
                END IF
             END DO
             IF(QPH/=0.0) TDG_PH = SUM_TDG_PHS/QPH                                                                ! FLOW AVERAGED TDG FOR POWER HOUSE(S)
             QTOT = 0.0                                                                                           ! INITIAL QTOT
             DO ig = 1, NGT
                QTOT = QTOT + QGT(ig)                                                                             ! SUM QTOT
             END DO
             IF (FLNO>0) THEN
                DO ifl = 1, FLNO
                   QTOT = QTOT - QGT(FLGTNO(ifl))                                                                 ! QTOT = QTOT - Q FISH LADDER
                END DO 
             END IF
             TDG_QROSP = Q_ROSP *35.3147/1000.0                                                                   ! CMS TO KCFS
             TDG_QPH   = QPH    *35.3147/1000.0                                                                   ! CMS TO KCFS
             TDG_QTOT  = QTOT   *35.3147/1000.0                                                                   ! CMS TO KCFS
             ! CALCULATE QENT
             IF (TDGENTEQ == 1) THEN
                TDG_QENT=E1*TDG_QROSP+E2
                !
             ELSE IF (TDGENTEQ == 2) THEN
                TDG_QENT=MIN((TDG_QTOT/60),1.0)*E1*TDG_QROSP+E2
                !
             ELSE IF (TDGENTEQ== 3) THEN
                TDG_QENT=MIN((TDG_QROSP/20),1.0)*E1*TDG_QROSP+E2
                !
             ELSE
                TEXT='TDGENTEQ INPUT ERROR'
                ERROR_OPEN=.TRUE.
             END IF
             TDG_QENT=MIN(TDG_QENT,(TDG_QTOT-TDG_QROSP))
             TDG_REL= (TDG_ROSP*(TDG_QROSP+TDG_QENT) + TDG_PH*(TDG_QPH-TDG_QENT))/(TDG_QPH+TDG_QROSP)           ! TDG RELEASE CALCULATION
             IF (TDG_REL>145.0) TDG_REL = 145.0                                                                 ! TDG RELEASE <=145.0
             TDG_TDG=TDG_REL                                                                                    ! OUTPUT TDG = TDG_REL
          END IF                                                                                                ! END IF POWNO >0 
          IF (TDG_TDG > 145.0) TDG_TDG = 145.0                                                                  ! TDG <=145.0
          IF (JDAY>=NXTSPLIT3) THEN
             DO ig = 1, NGT
                SUM_QGT2=SUM_QGT2+QGT(ig)
             END DO
             WRITE (88888, '(A, F10.3, 2A, F10.3, A, F9.3, A, <NGT>(F9.3,","))')' ',JDAY,',  ', ', ',TDG_TDG,',  ',SUM_QGT2,',  ',(QGT(ig), ig = 1, NGT)
             NXTSPLIT3 = NXTSPLIT3 + 1.0
          END IF
       END IF                                                                                                   ! END IF Q_SUM/=0.0
    END SUBROUTINE SYSTDG_TDG
    !===========================================================================================================================
    SUBROUTINE DEALLOCATE_SYSTDG
     IMPLICIT NONE 
      DEALLOCATE (TDG_PHS, TDG_FLS, TDG_TDP)
      DEALLOCATE (BAYC, QBAY)
      DEALLOCATE (POWGTNO, FLGTNO)
      DEALLOCATE (GTNAME)
    END SUBROUTINE
    !
END MODULE modSYSTDG