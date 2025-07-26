!***********************************************************************************************************************************
!**    S U B R O U T I N E    I N I T I A L    W A T E R    L E V E L                                                            **
!**********************************************************************************************************************************

subroutine initial_water_level

USE MAIN
USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  USE MACROPHYTEC; USE POROSITYC; USE ZOOPLANKTONC 
  USE INITIALVELOCITY
  IMPLICIT NONE
  EXTERNAL RESTART_OUTPUT
  INTEGER :: JBU,JBD, JJW
  REAL :: QGATE,WSUP,BRLEN,WLSLOPE,DIST, WLDIFF
  
  LOOP_BRANCH=.FALSE.
    
! estimating initial flows in each segment
     QSSI=0.0
            
! first considering specified flows: upstream inflows, tributaries, distributed tribs, structural withdrawals and withdrawals    
    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        IU = CUS(JB)
        ID = DS(JB)
        DO I=IU,ID
          IF(I == IU .AND. UP_FLOW(JB))THEN
             QSSI(I)=QIN(JB)+QSSI(I)
          END IF
! if downstream of structure          
          !IF (i == iu .and. UHS(JB) < 0) then
          IF (I == IU .AND. DAM_INFLOW(JB)) THEN    ! CB 4/27/2011
            DO JJW=1,NWB
              DO JJB=BS(JJW),BE(JJW)
                IF(DS(JJB) == ABS(UHS(JB)))THEN
                  DO JS=1,NSTR(JJB)              
                    QSSI(I)=QSSI(I)+QSTR(JS,JJB)
                  END DO
                END IF
              END DO
            END DO
          END IF

          IF (TRIBUTARIES) THEN
            DO JT=1,NTR
              IF(ITR(JT) == I)QSSI(I)=QSSI(I)+QTR(JT)
            END DO
          END IF
          IF (DIST_TRIBS(JB)) THEN
            QSSI(I)=QSSI(I)+QDTR(JB)/REAL(ID-IU+1)    ! SINCE INITIAL WL UNKNOWN, DISTRIBUTING FLOW EVENLY BTW. SEGS.              
          END IF
          IF (WITHDRAWALS) THEN
            DO JWD=1,NWD
              IF(IWD(JWD) == I)QSSI(I)=QSSI(I)-QWD(JWD)            
            END DO
          END IF
          IF(I == ID)THEN
            DO JS=1,NSTR(JB)              
              QSSI(I)=QSSI(I)-QSTR(JS,JB)
            END DO
          END IF
        END DO
      END DO
    END DO
    
! INCLUDING FLOWS WITHIN BRANCH UPSTEAM OF SEGMENT
    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        IU = CUS(JB)
        ID = DS(JB)
        DO I=IU+1,ID
          QSSI(I)=QSSI(I)+QSSI(I-1)
        END DO
      END DO
    END DO

    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        IU = CUS(JB)
        ID = DS(JB)
        DO I=IU,ID          
! DETERMINING IF SEGMENT IS DOWNSTREAM INTERNAL HEAD BOUNDARY OF ANOTHER *UPSTREAM* BRANCH, AND ADDING FLOW TO SEGMENT AND SEGMENTS DOWNSTREAM
          DO JJW=1,NWB
            DO JJB=BS(JJW),BE(JJW)
              IF(DHS(JJB) == I)THEN
                DO II=I,ID
                  QSSI(II)=QSSI(II)+QSSI(DS(JJB))
                END DO
              END IF
            END DO
          END DO                    
! DETERMINING IF SEGMENT IS DOWNSTREAM OF SPILLWAY BELOW ANOTHER BRANCH
          DO JS=1,NSP
            IF(ESP(JS) < EL(2,I))THEN  ! DISREGARDING IF CREST ABOVE GRID
              IF(I == IDSP(JS))THEN
                DO II=I,ID
                  QSSI(II)=QSSI(II)+QSSI(IUSP(JS))
                END DO
              END IF
            END IF
          END DO
! DETERMINING IF SEGMENT IS DOWNSTREAM OF GATE BELOW ANOTHER BRANCH
          DO JG=1,NGT
            IF(EGT(JG) < EL(2,I))THEN  ! DISREGARDING IF CREST ABOVE GRID
              IF(I == IDGT(JG))THEN
                IF(DYNGTC(JG) == '    FLOW')THEN
                  QGATE = BGT(JG)
                  DO II=I,ID
                    QSSI(II)=QSSI(II)+QGATE
                  END DO
                ELSE
                  DO II=I,ID
                    QSSI(II)=QSSI(II)+QSSI(IUGT(JG))
                  END DO
                END IF
              END IF
            END IF
          END DO                    
        END DO
      END DO
    END DO

! DETERMINING IF BRANCH IS A SECONDARY BRANCH THAT LOOPS AROUND AN ISLAND WITH INTERNAL HEAD BOUNDARIES AT UPSTREAM AND DOWNSTREAM BC
! ATTACHED TO A SINGLE BRANCH
    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        IF(UH_INTERNAL(JB) .AND. DH_INTERNAL(JB))THEN
          DO JJW=1,NWB
              DO JJB=BS(JJW),BE(JJW)
                IF(UHS(JB) > US(JJB) .AND. UHS(JB) < DS(JJB))THEN
                  JBU=JJB
                END IF
              END DO
          END DO
          DO JJW=1,NWB
              DO JJB=BS(JJW),BE(JJW)
                IF(DHS(JB) > US(JJB) .AND. DHS(JB) < DS(JJB))THEN       ! WW 8/19/2013
                  JBD=JJB
                END IF
              END DO
          END DO
          LOOP_BRANCH(JB) = JBU == JBD
        END IF
      END DO
    END DO
              
              
! GIVEN ESTIMATED FLOWS FOR EACH SEGMENT, ESTIMATING WL WITH NORMAL DEPTH EQUATION
    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
! ONLY CONSIDERING BRANCHES WITH SLOPES > 0      
        IF (SLOPE(JB) > 0.0 .AND. .NOT. LOOP_BRANCH(JB))THEN
          IU = CUS(JB)
          ID = DS(JB)    
          DO I=IU,ID
            CALL NORMAL_DEPTH(QSSI(I))
          END DO
        END IF
      END DO
    END DO
    
! SMOOTHING WATER SURFACE ELEVATIONS, FIRST WITHIN BRANCHES
    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        IF (SLOPE(JB) > 0.0 .AND. .NOT. LOOP_BRANCH(JB))THEN
          IU = CUS(JB)
          ID = DS(JB)
          DO I=ID-1,IU,-1
            IF(ELWS(I+1) > ELWS(I))THEN
              ELWS(I) = ELWS(I+1)
            END IF          
          END DO
        END IF
      END DO
    END DO

! SMOOTHING WATER SURFACE ELEVATIONS AT INTERNAL HEAD BOUNDARIES
    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        IF (SLOPE(JB) > 0.0 .AND. .NOT. LOOP_BRANCH(JB) .AND. DH_INTERNAL(JB))THEN
          IU = CUS(JB)
          ID = DS(JB)
          IF(ELWS(DHS(JB)) > ELWS(ID))THEN
            ELWS(ID)=ELWS(DHS(JB))
            DO I=ID-1,IU,-1
              IF(ELWS(I+1) > ELWS(I))THEN
                ELWS(I) = ELWS(I+1)
              END IF   
            END DO
          END IF
        END IF
      END DO
    END DO

! IF SPILLWAY OR GATE AT DOWNSTREAM END OF BRANCH, MAKING SURE WATER LEVEL IS ABOVE CREST ELEVATION
    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        IF (SLOPE(JB) > 0.0 .AND. .NOT. LOOP_BRANCH(JB) .AND. NSTR(JB) == 0)THEN
          IU = CUS(JB)
          ID = DS(JB)
! SPILLWAYS          
          DO JS=1,NSP
            IF(ESP(JS) < EL(2,I))THEN  ! DISREGARDING IF CREST ABOVE GRID
              IF(ID == IUSP(JS))THEN
                 WSUP=ESP(JS)+(QSSI(ID)/A1SP(JS))**(1.0/B1SP(JS))  ! ESTIMATING UPSTREAM WS ELEV
                 IF(ELWS(ID) < WSUP)THEN
                   IF(IDSP(JS) .NE. 0)THEN    ! CB 8/10/10
                   IF(ELWS(IDSP(JS)) > WSUP) WSUP = ELWS(IDSP(JS))   ! CHECKING TO SEE IF DOWNSTREAM WS ELEVATION ISN'T ALREADY 'HIGH'
                   ELWS(ID)=WSUP
                   END IF                     ! CB 8/10/10
                   DO I=ID-1,IU,-1
                     IF(ELWS(I+1) > ELWS(I))THEN
                       ELWS(I) = ELWS(I+1)
                     END IF   
                   END DO
                 END IF
               END IF
            END IF
          END DO
! GATES
          DO JG=1,NGT
            IF(EGT(JG) < EL(2,I) .AND. DYNGTC(JG) .NE. '    FLOW')THEN  ! DISREGARDING IF CREST ABOVE GRID
              IF(ID == IUGT(JG))THEN
                IF(DYNGTC(JG) == '       B')THEN
                  WSUP=EGT(JG)+(QSSI(ID)/(A1GT(JG)*BGT(JG)**G1GT(JG)))**(1.0/B1GT(JG))
                END IF
                IF(DYNGTC(JG) == '     ZGT')THEN
                  WSUP=EGT(JG)+(QSSI(ID)/A1GT(JG))**(1.0/B1GT(JG))
                END IF
                IF(ELWS(ID) < WSUP)THEN
                   IF(ELWS(IDGT(JG)) > WSUP) WSUP = ELWS(IDGT(JG))   ! CHECKING TO SEE IF DOWNSTREAM WS ELEVATION ISN'T ALREADY 'HIGH'  WX 8/21/13
                   ELWS(ID)=WSUP
                   DO I=ID-1,IU,-1
                     IF(ELWS(I+1) > ELWS(I))THEN
                       ELWS(I) = ELWS(I+1)
                     END IF   
                   END DO
                END IF
              END IF
            END IF
          END DO
        END IF
      END DO
    END DO

! SMOOTHING WATER LEVEL AROUND LOOP BRANCHES
    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        IU = CUS(JB)
        ID = DS(JB)
        IF (LOOP_BRANCH(JB))THEN
          WLDIFF=ELWS(UHS(JB))-ELWS(UHS(JB))          
          BRLEN=0.0
          DO I=IU,ID
            BRLEN=BRLEN+DLX(I)
          END DO
          WLSLOPE=WLDIFF/BRLEN
          DIST=DLX(IU)/2.0
          ELWS(IU)=ELWS(UHS(JB))-WLSLOPE*DIST
          DO I=IU+1,ID
            DIST=DIST+(DLX(I-1)+DLX(I))/2.0
            ELWS(IU)=ELWS(UHS(JB))-WLSLOPE*DIST
          END DO
        END IF
      END DO
    END DO

    RETURN
    END SUBROUTINE INITIAL_WATER_LEVEL

!***********************************************************************************************************************************
!**        S U B R O U T I N E    N O R M A L    D E P T H                                                                        **
!***********************************************************************************************************************************

      SUBROUTINE NORMAL_DEPTH(FLOW)
      
      USE GLOBAL; USE GEOMC
                  
      INTEGER, PARAMETER ::JMAX=40    
      REAL(R8)   :: FLOW
      REAL :: X1, X2, FUNCVAL1, FUNCVAL2, XACC, FMID, FUNC1, RTBIS, DX, XMID     
      INTEGER :: JJ,J  
                  
! FIRST, BRACKETING ROOT
      X1=0.001
      X2=1.0
      CALL MANNINGS_EQN(FLOW,X1,FUNCVAL1)
      CALL MANNINGS_EQN(FLOW,X2,FUNCVAL2)
                  
      DO JJ=1,JMAX
        IF(FUNCVAL1*FUNCVAL2 > 0.0)THEN
          IF(ABS(FUNCVAL1)  < ABS(FUNCVAL2))THEN
            X1=X1/2.0
            CALL MANNINGS_EQN(FLOW,X1,FUNCVAL1)            
          ELSE
            X2=X2+1.5*(X2-X1)
            CALL MANNINGS_EQN(FLOW,X2,FUNCVAL2)
          END IF                                          
        ELSE
          EXIT
        END IF      
      END DO                 
      
! FINDING ROOT BY BISECTION      
      XACC=0.01
      CALL MANNINGS_EQN(FLOW,X2,FMID)
      CALL MANNINGS_EQN(FLOW,X1,FUNC1)      
  !    IF(FUNC1*FMID.GE.0.) PAUSE 'ROOT MUST BE BRACKETED IN RTBIS'
      IF(FUNC1.LT.0.)THEN
        RTBIS=X1
        DX=X2-X1
      ELSE
        RTBIS=X2
        DX=X1-X2
      ENDIF
      DO J=1,JMAX
        DX=DX*.5
        XMID=RTBIS+DX
        CALL MANNINGS_EQN(FLOW,XMID,FMID)        
        IF(FMID.LE.0.)RTBIS=XMID
        IF(ABS(DX).LT.XACC .OR. FMID.EQ.0.)THEN          
          ELWS(I)=RTBIS+EL(KB(I)+1,I)                                      ! SW 4/5/13
          RETURN
        END IF
      END DO
  !    PAUSE 'TOO MANY BISECTIONS IN RTBIS'
      END SUBROUTINE NORMAL_DEPTH

      
!***********************************************************************************************************************************
!**        S U B R O U T I N E    M A N N I N G S    E Q U A T I O N                                                              **
!***********************************************************************************************************************************

      SUBROUTINE MANNINGS_EQN(FLOW,DEPTH,FUNCVALUE)
      
      USE GLOBAL; USE GEOMC; USE EDDY; USE LOGICC
      
      REAL(R8) :: FLOW
      REAL     :: WSURF, DEPTH, XAREA, WPER, HRAD,FMANN, FUNCVALUE

 !     WSURF=EL(KB(I)-1,I)+DEPTH
      WSURF=EL(KB(I)+1,I)+DEPTH   ! CB 7/7/10
      CALL XSECTIONAL_AREA(WSURF,XAREA)
      WPER=B(KTI(I),I)+2.0*DEPTH
      HRAD=XAREA/WPER
      IF(MANNINGS_N(JW))THEN
        FMANN=FRIC(I)
      ELSE
        FMANN=HRAD**0.166666667/FRIC(I)
      END IF
      FUNCVALUE=FLOW-XAREA*HRAD**0.6667*SLOPEC(JB)**0.5/FMANN            ! SW 4/5/2013
      
      RETURN
      END SUBROUTINE MANNINGS_EQN
      
!***********************************************************************************************************************************
!**        S U B R O U T I N E    C R O S S    S E C T I O N A L    A R E A                                                       **
!***********************************************************************************************************************************

      SUBROUTINE XSECTIONAL_AREA(WSURF,XAREA)
      
      USE GLOBAL; USE GEOMC; USE MAIN; USE INITIALVELOCITY
      
      REAL :: XAREA, WSURF     ! 4/5/13 SW
      INTEGER :: KTTOP         ! 4/5/13 SW

 !     KTTOP = 2
      DO K=2,KMX-1
        IF(EL(K,I) < WSURF)THEN
          KTTOP=K-1
          EXIT
        END IF
        KTTOP=K     ! CB 8/10/10
      END DO
 !     DO WHILE (EL(KTTOP,I) > WSURF)
 !        KTTOP = KTTOP+1
 !     END DO            
      XAREA=(WSURF-EL(KTTOP+1,I)) * BSAVE(KTTOP,I)
      DO K=KTTOP+1,KBI(I)
         XAREA = XAREA+BSAVE(K,I)*H(K,JW)
      END DO      
      
      RETURN
      END  SUBROUTINE XSECTIONAL_AREA
      
      
!***********************************************************************************************************************************
!**        S U B R O U T I N E    I N I T I A L    H O R I Z O N T A L    V E L O C I T Y                                         **
!***********************************************************************************************************************************

      SUBROUTINE INITIAL_U_VELOCITY
      
      USE GLOBAL; USE GEOMC
      USE INITIALVELOCITY
      
      REAL :: XAREA, WSURF
      INTEGER :: K

      DO JW=1,NWB
        KT = KTWB(JW)
        
        DO JB=BS(JW),BE(JW)
          IF (SLOPE(JB) > 0.0 .AND. .NOT. LOOP_BRANCH(JB)) THEN
            IU = CUS(JB)
            ID = DS(JB)
            DO I=IU,ID
              WSURF=ELWS(I)              
              CALL XSECTIONAL_AREA(WSURF,XAREA)
              UAVG(I)=QSSI(I)/XAREA
              DO K=KT,KB(I)
                U(K,I)=UAVG(I)
              END DO
            END DO
          END IF
        END DO
      END DO
                
      RETURN
      END SUBROUTINE INITIAL_U_VELOCITY