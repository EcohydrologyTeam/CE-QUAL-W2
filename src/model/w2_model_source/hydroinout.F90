SUBROUTINE HYDROINOUT

USE MAIN
USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  USE MACROPHYTEC; USE POROSITYC; USE ZOOPLANKTONC  
  USE modSYSTDG, ONLY: GTNAME, SYSTDG_TDG, UPDATE_TDGC, TDG_TDG, TDG_ROSP, POWNO, POWGTNO, FLNO, FLGTNO, TDGLOC, ip, ifl          ! systdg
  IMPLICIT NONE
  EXTERNAL RESTART_OUTPUT
  INTEGER   :: JBU, JBD, JLAT, JWU
  REAL(R8)  :: ELW, CGAS, TM, VPTG, DTVL, RHOTR, VQTR, VQTRI, QTRFR, AKBR, FW
  REAL(R8)  :: TSUM,QSUMM
  REAL(R8)  :: dosat, n2sat   ! systdg - update powerhouse release tdg
!***********************************************************************************************************************************
!**                                            Task 2.1: Hydrodynamic sources/sinks                                               **
!***********************************************************************************************************************************

    QINSUM = 0.0; TINSUM = 0.0; CINSUM = 0.0; UXBR = 0.0; UYBR = 0.0;  tdgon=.false.        ! cb 1/16/13
    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        IF(BR_INACTIVE(JB))CYCLE
        IU   = CUS(JB)
        ID   = DS(JB)
        TSUM = 0.0; CSUM = 0.0; QSUM(JB) = 0.0; QOUT(:,JB) = 0.0; TOUT(JB)=0.0; COUT(:,JB)=0.0    

!****** Densities

        DO I=IU-1,ID+1
          DO K=KT,KB(I)
            TISS(K,I) = 0.0
            DO JS=1,NSS
              TISS(K,I) = TISS(K,I)+SS(K,I,JS)
            END DO
            RHO(K,I) = DENSITY(T2(K,I),DMAX1(TDS(K,I),0.0D0),DMAX1(TISS(K,I),0.0D0))
          END DO
        END DO
! v3.5 deleted pumpback code from v3.2
        DO JS=1,NSTR(JB)
          IF (QSTR(JS,JB) /= 0.0) THEN
            CALL DOWNSTREAM_WITHDRAWAL (JS)
          END IF
        END DO
        DO K=KT,KB(ID)
          QSUM(JB)        = QSUM(JB)       +QOUT(K,JB)
          TSUM            = TSUM           +QOUT(K,JB)*T2(K,ID)
          CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QOUT(K,JB)*C2(K,ID,CN(1:NAC))
        END DO
        IF (QSUM(JB) /= 0.0) THEN
          TOUT(JB)           = TSUM           /QSUM(JB)
          COUT(CN(1:NAC),JB) = CSUM(CN(1:NAC))/QSUM(JB)
        END IF
        IF (QSUM(JB) /= 0.0 .AND. DAM_OUTFLOW(JB)) THEN                                                                !TC 08/03/04
          TINSUM(JBDAM(JB))           = (TSUM           +QINSUM(JBDAM(JB))*TINSUM(JBDAM(JB)))          /(QSUM(JB)                  &
                                        +QINSUM(JBDAM(JB)))
          CINSUM(CN(1:NAC),JBDAM(JB)) = (CSUM(CN(1:NAC))+QINSUM(JBDAM(JB))*CINSUM(CN(1:NAC),JBDAM(JB)))/(QSUM(JB)                  &
                                        +QINSUM(JBDAM(JB)))
          QINSUM(JBDAM(JB))           =  QINSUM(JBDAM(JB))+QSUM(JB)
        END IF
      END DO
    END DO
    ILAT = 0
    JWW  = NWD
    withdrawals = jww > 0
    if(nwdt>nwd)qwd(nwd+1:nwdt)=0.0    ! SW 10/30/2017
    JTT  = NTR
    tributaries = jtt > 0
    if(ntrt>ntr)qtr(ntr+1:ntrt)=0.0    ! SW 10/30/2017
    JSS  = NSTR
    IF (SPILLWAY) THEN
      CALL SPILLWAY_FLOW
      DO JS=1,NSP

!****** Positive flows

        JLAT = 0
        JBU  = JBUSP(JS)
        JBD  = JBDSP(JS)
        tdgon=.false.        ! cb 1/16/13
        jsg=js
        nnsg=0
        if (cac(ndo) == '      ON' .and. gasspc(js) == '      ON')tdgon=.true.
        IF (QSP(JS) >= 0.0) THEN
          IF (LATERAL_SPILLWAY(JS)) THEN
            JWW       = JWW+1
            IWD(JWW)  = IUSP(JS)
            QWD(JWW)  = QSP(JS)
            KTWD(JWW) = KTUSP(JS)
            KBWD(JWW) = KBUSP(JS)
            EWD(JWW)  = ESP(JS)
            JBWD(JWW) = JBU
            I         = MAX(CUS(JBWD(JWW)),IWD(JWW))
            JB        = JBWD(JWW)
            JW        = JWUSP(JS)
            KT        = KTWB(JW)
            jwd=jww
            CALL LATERAL_WITHDRAWAL !(JWW)
            DO K=KTW(JWW),KBW(JWW)
              QSS(K,I) = QSS(K,I)-QSW(K,JWW)
            END DO
                IF (IDSP(JS) /= 0)then   ! cb 9/11/13
            JTT  = JTT+1
              QTR(JTT)         = QSP(JS)
              ITR(JTT)         = IDSP(JS)
              PLACE_QTR(JTT)   = PDSPC(JS) == ' DENSITY'
              SPECIFY_QTR(JTT) = PDSPC(JS) == ' SPECIFY'
              IF (SPECIFY_QTR(JTT)) THEN
                ELTRT(JTT) = ETDSP(JS)
                ELTRB(JTT) = EBDSP(JS)
              END IF
              JBTR(JTT) = JBD
                end if  ! cb 9/11/13
            IF (IDSP(JS) /= 0 .AND. QSP(JS)>0.0) THEN
            TSUM  =  0.0; QSUMM = 0.0; CSUM = 0.0
            DO K=KTW(JWW),KBW(JWW)
              QSUMM           = QSUMM          +QSW(K,JWW)
              TSUM            = TSUM           +QSW(K,JWW)*T2(K,IWD(JWW))
              CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QSW(K,JWW)*C2(K,IWD(JWW),CN(1:NAC))
            END DO
            TTR(JTT) = TSUM/QSUMM
            DO JC=1,NAC
              CTR(CN(JC),JTT) = CSUM(CN(JC))/QSUMM
              IF (CN(JC) == NDO .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN
                TDG_SPILLWAY(JWW,JS) = .TRUE.
                CALL TOTAL_DISSOLVED_GAS (0,PALT(I),0,JS,TTR(JTT),CTR(CN(JC),JTT))   ! DO
              END IF
              IF (CN(JC) == NN2 .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN
                TDG_SPILLWAY(JWW,JS) = .TRUE.
                CALL TOTAL_DISSOLVED_GAS (1,PALT(I),0,JS,TTR(JTT),CTR(CN(JC),JTT))   ! N2
              END IF
              IF (CN(JC) == NDGP .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN      !8/2020 TDGP
                TDG_SPILLWAY(JWW,JS) = .TRUE.
                CALL TOTAL_DISSOLVED_GAS (2,PALT(I),0,JS,TTR(JTT),CTR(CN(JC),JTT))   
              END IF
            END DO
          ELSE IF (CAC(NDO) == '      ON' .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN
            TDG_SPILLWAY(JWW,JS) = .TRUE.
          END IF
            
          ELSE
            JSS(JBU)                 =  JSS(JBU)+1
            KTSW(JSS(JBU),JBU)       =  KTUSP(JS)
            KBSW(JSS(JBU),JBU)       =  KBUSP(JS)
            JB                       =  JBU
            POINT_SINK(JSS(JBU),JBU) = .TRUE.
            ID                       =  IUSP(JS)
            QSTR(JSS(JBU),JBU)       =  QSP(JS)
            ESTR(JSS(JBU),JBU)       =  ESP(JS)
            KT                       =  KTWB(JWUSP(JS))
            JW                       =  JWUSP(JS)
            CALL DOWNSTREAM_WITHDRAWAL(JSS(JBU))           
              QSUM(JB) = 0.0; TSUM = 0.0; CSUM = 0.0
            DO K=KT,KB(ID)
              QSUM(JB) = QSUM(JB)+QOUT(K,JB)
              TSUM     = TSUM+QOUT(K,JB)*T2(K,ID)
              DO JC=1,NAC
                IF (CN(JC)==NDO .AND. CAC(NDO) == '      ON' .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN     ! MM 5/21/2009
                  T2R4=T2(K,ID)
                  CGAS=C2(K,ID,CN(JC))                                                                                      ! MM 5/21/2009
                  CALL TOTAL_DISSOLVED_GAS (0,PALT(ID),0,JS,T2R4,CGAS)    ! O2
                  CSUM(CN(JC)) = CSUM(CN(JC))+QOUT(K,JB)*CGAS
                !ELSEIF (CN(JC)==NGN2 .AND. CAC(NGN2) == '      ON' .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN     ! SW 10/27/15
                ELSEIF (CN(JC)==NN2 .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN     ! SW 10/27/15
                  if(CAC(NN2) == '      ON')then                                                                                ! cb 1/13/16
                    T2R4=T2(K,ID)
                    CGAS=C2(K,ID,CN(JC))                                                                                      ! 
                    CALL TOTAL_DISSOLVED_GAS (1,PALT(ID),0,JS,T2R4,CGAS)   ! N2
                    CSUM(CN(JC)) = CSUM(CN(JC))+QOUT(K,JB)*CGAS
                  end if
                ELSEIF (CN(JC)==NDGP .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN     !8/2020 TDGP
                  if(CAC(NDGP) == '      ON')then                                                                                
                    T2R4=T2(K,ID)
                    CGAS=C2(K,ID,CN(JC))                                                                                      
                    CALL TOTAL_DISSOLVED_GAS (2,PALT(ID),0,JS,T2R4,CGAS)   
                    CSUM(CN(JC)) = CSUM(CN(JC))+QOUT(K,JB)*CGAS
                  end if  
                ELSE
                  CSUM(CN(JC)) = CSUM(CN(JC))+QOUT(K,JB)*C2(K,ID,CN(JC))
                END IF
              END DO
            END DO
            IF (QSUM(JB) /= 0.0) THEN
              TOUT(JB)           = TSUM           /QSUM(JB)
              COUT(CN(1:NAC),JB) = CSUM(CN(1:NAC))/QSUM(JB)
            END IF
            IF (IDSP(JS) /= 0 .AND. US(JBD) == IDSP(JS)) THEN
              QSUMM = 0.0; TSUM  = 0.0;  CSUM  = 0.0
              DO K=KT,KB(ID)
                QSUMM = QSUMM+QNEW(K)
                TSUM  = TSUM+QNEW(K)*T2(K,ID)
                DO JC=1,NAC
                  IF (CN(JC)==NDO .AND. CAC(NDO) == '      ON' .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN             ! MM 5/21/2009
                    T2R4=T2(K,ID)
                    CGAS=C2(K,ID,CN(JC))                                                                                            ! MM 5/21/2009
                    CALL TOTAL_DISSOLVED_GAS (0,PALT(ID),0,JS,T2R4,CGAS)
                    CSUM(CN(JC)) = CSUM(CN(JC))+QNEW(K)*CGAS
                  !ELSEIF (CN(JC)==NGN2 .AND. CAC(NGN2) == '      ON' .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN             ! SW 10/27/15
                  ELSEIF (CN(JC)==NN2 .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN             ! SW 10/27/15
                    if(CAC(NN2) == '      ON')then                                                                                ! cb 1/13/16
                      T2R4=T2(K,ID)
                      CGAS=C2(K,ID,CN(JC))                                                                                            ! 
                      CALL TOTAL_DISSOLVED_GAS (1,PALT(ID),0,JS,T2R4,CGAS)
                      CSUM(CN(JC)) = CSUM(CN(JC))+QNEW(K)*CGAS
                    end if
                  ELSEIF (CN(JC)==NDGP .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN    !8/2020 TDGP        
                    if(CAC(NDGP) == '      ON')then                                                                       
                      T2R4=T2(K,ID)
                      CGAS=C2(K,ID,CN(JC))                                                                                            
                      CALL TOTAL_DISSOLVED_GAS (2,PALT(ID),0,JS,T2R4,CGAS)
                      CSUM(CN(JC)) = CSUM(CN(JC))+QNEW(K)*CGAS
                    end if  
                  ELSE
                    CSUM(CN(JC)) = CSUM(CN(JC))+QNEW(K)*C2(K,ID,CN(JC))
                  END IF
                END DO
              END DO
              IF (QSUMM /= 0.0) THEN
                TINSUM(JBD)           = (TSUM           +QINSUM(JBD)*TINSUM(JBD))          /(QSUMM+QINSUM(JBD))
                CINSUM(CN(1:NAC),JBD) = (CSUM(CN(1:NAC))+QINSUM(JBD)*CINSUM(CN(1:NAC),JBD))/(QSUMM+QINSUM(JBD))
                QINSUM(JBD)           =  QINSUM(JBD)+QSUMM
              END IF
        ELSEIF (IDSP(JS) /= 0) THEN
              JTT              = JTT+1
              QTR(JTT)         = QSP(JS)
              ITR(JTT)         = IDSP(JS)
              PLACE_QTR(JTT)   = PDSPC(JS) == ' DENSITY'
              SPECIFY_QTR(JTT) = PDSPC(JS) == ' SPECIFY'
              IF (SPECIFY_QTR(JTT)) THEN
                ELTRT(JTT) = ETDSP(JS)
                ELTRB(JTT) = EBDSP(JS)
              END IF
              JBTR(JTT) = JBD
                  TTR(JTT) =  TOUT(JB)
                DO JC=1,NAC
                  CTR(CN(JC),JTT) = COUT(CN(JC),JB)
                  IF (CN(JC) == NDO .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN
                    CALL TOTAL_DISSOLVED_GAS (0,PALT(ITR(JTT)),0,JS,TTR(JTT),CTR(CN(JC),JTT))
                  END IF
                  IF (CN(JC) == NN2 .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN     ! SW 10/27/15
                    CALL TOTAL_DISSOLVED_GAS (1,PALT(ITR(JTT)),0,JS,TTR(JTT),CTR(CN(JC),JTT))
                  END IF
                  IF (CN(JC) == NDGP .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN     
                    CALL TOTAL_DISSOLVED_GAS (2,PALT(ITR(JTT)),0,JS,TTR(JTT),CTR(CN(JC),JTT))
                  END IF
                END DO
               END IF
          END IF
        ELSE IF (QSP(JS) < 0.0) THEN
          JTT              =  JTT+1
          JWW              =  JWW+1
          IWD(JWW)         =  IDSP(JS)
          ITR(JTT)         =  IUSP(JS)
          QTR(JTT)         = -QSP(JS)
          QWD(JWW)         = -QSP(JS)
          KTWD(JWW)        =  KTDSP(JS)
          KBWD(JWW)        =  KBDSP(JS)
          EWD(JWW)         =  ESP(JS)
          PLACE_QTR(JTT)   =  PUSPC(JS) == ' DENSITY'
          SPECIFY_QTR(JTT) =  PUSPC(JS) == ' SPECIFY'
          IF (SPECIFY_QTR(JTT)) THEN
            ELTRT(JTT) = ETUSP(JS)
            ELTRB(JTT) = EBUSP(JS)
          END IF
          JBTR(JTT) = JBU
          JBWD(JWW) = JBD
          I         = MAX(CUS(JBWD(JWW)),IWD(JWW))
          JB        = JBWD(JWW)
          JW        = JWDSP(JS)
          KT        = KTWB(JW)
          jwd=jww
          CALL LATERAL_WITHDRAWAL !(JWW)
          DO K=KTW(JWW),KBW(JWW)
            QSS(K,I) = QSS(K,I)-QSW(K,JWW)
          END DO
          IF (IDSP(JS) /= 0) THEN
            TSUM  =  0.0; QSUMM = 0.0; CSUM = 0.0
            DO K=KTW(JWW),KBW(JWW)
              QSUMM           = QSUMM          +QSW(K,JWW)
              TSUM            = TSUM           +QSW(K,JWW)*T2(K,IWD(JWW))
              CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QSW(K,JWW)*C2(K,IWD(JWW),CN(1:NAC))
            END DO
            TTR(JTT) = TSUM/QSUMM
            DO JC=1,NAC
              CTR(CN(JC),JTT) = CSUM(CN(JC))/QSUMM
              IF (CN(JC) == NDO .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN
                TDG_SPILLWAY(JWW,JS) = .TRUE.
                CALL TOTAL_DISSOLVED_GAS (0,PALT(I),0,JS,TTR(JTT),CTR(CN(JC),JTT))    ! O2
              END IF
              IF (CN(JC) == NN2 .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN
                TDG_SPILLWAY(JWW,JS) = .TRUE.
                CALL TOTAL_DISSOLVED_GAS (1,PALT(I),0,JS,TTR(JTT),CTR(CN(JC),JTT))    ! N2
              END IF
              IF (CN(JC) == NDGP .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN     !8/2020 TDGP
                TDG_SPILLWAY(JWW,JS) = .TRUE.
                CALL TOTAL_DISSOLVED_GAS (2,PALT(I),0,JS,TTR(JTT),CTR(CN(JC),JTT))    
              END IF
            END DO
          ELSE IF (CAC(NDO) == '      ON' .AND. GASSPC(JS) == '      ON' .AND. QSP(JS) > 0.0) THEN
            TDG_SPILLWAY(JWW,JS) = .TRUE.
          END IF
        END IF
      END DO
    END IF
    IF (PUMPS) THEN
      DO JP=1,NPU
        JLAT = 0
        JWU  = JWUPU(JP)
        JBU  = JBUPU(JP)
        JBD  = JBDPU(JP)
        tdgon=.false.        ! cb 1/16/13
          IF (LATERAL_PUMP(JP)) THEN
            IF(PUMP_DOWNSTREAM(JP))THEN
                ELW = EL(KTWB(JWDPU(JP)),IDPU(JP))-Z(IDPU(JP))*COSA(JBD)
                JWW       = JWW+1      ! SW 10/30/2017
                JBWD(JWW) = JBU  
                IWD(JWW)  = IUPU(JP)
            ELSE
            ELW = EL(KTWB(JWU),IUPU(JP))-Z(IUPU(JP))*COSA(JBU)
            JWW       = JWW+1      ! SW 10/30/2017
            JBWD(JWW) = JBU  
            IWD(JWW)  = IUPU(JP)
            ENDIF
                ELSE
                    IF(PUMP_DOWNSTREAM(JP))THEN
                        ELW = EL(KTWB(JWDPU(JP)),IDPU(JP))-Z(IDPU(JP))*COSA(JBD)-SINA(JBD)*DLX(IDPU(JP))*0.5
                        JSS(JBU)                 =  JSS(JBU)+1     ! SW 10/30/2017
                    ELSE
            ELW = EL(KTWB(JWU),IUPU(JP))-Z(IUPU(JP))*COSA(JBU)-SINA(JBU)*DLX(IUPU(JP))*0.5
            JSS(JBU)                 =  JSS(JBU)+1     ! SW 10/30/2017
                   ENDIF
          END IF
          
        IF (JDAY >= ENDPU(JP)) PUMPON(JP) = .FALSE.                                                        !  CB 1/13/06
        IF (JDAY >= STRTPU(JP) .AND. JDAY < ENDPU(JP)) THEN
            IF(PUMP_DOWNSTREAM(JP))THEN    ! IF BASED ON DOWNSTREAM WATER LEVEL AND NOT UPSTREAM
            IF (ELW >= EOFFPU(JP)) PUMPON(JP) = .FALSE.                                                       ! CB 1/13/06
            IF (ELW < EOFFPU(JP) .AND. QPU(JP) > 0.0) THEN
            IF (ELW <= EONPU(JP)) PUMPON(JP) = .TRUE.
            IF (PUMPON(JP)) THEN
              IF (LATERAL_PUMP(JP)) THEN
                JLAT      = 1
                QWD(JWW)  = QPU(JP)
                KTWD(JWW) = KTPU(JP)
                KBWD(JWW) = KBPU(JP)
                EWD(JWW)  = EPU(JP)
                I         = MAX(CUS(JBWD(JWW)),IWD(JWW))
                JB        = JBWD(JWW)
                JW        = JWU
                KT        = KTWB(JW)
                jwd=jww
                CALL LATERAL_WITHDRAWAL         
                DO K=KTW(JWW),KBW(JWW)
                  QSS(K,I) = QSS(K,I)-QSW(K,JWW)
                END DO
                IF (IDPU(JP) /= 0) THEN           
                  JTT              = JTT+1
                  QTR(JTT)         = QPU(JP)
                  ITR(JTT)         = IDPU(JP)
                  PLACE_QTR(JTT)   = PPUC(JP) == ' DENSITY'
                  SPECIFY_QTR(JTT) = PPUC(JP) == ' SPECIFY'
                  IF (SPECIFY_QTR(JTT)) THEN
                    ELTRT(JTT) = ETPU(JP)
                    ELTRB(JTT) = EBPU(JP)
                  END IF
                  JBTR(JTT) = JBD
                    TSUM = 0.0; QSUMM = 0.0; CSUM(CN(1:NAC)) = 0.0
                    DO K=KTW(JWW),KBW(JWW)
                      QSUMM           = QSUMM          +QSW(K,JWW)
                      TSUM            = TSUM           +QSW(K,JWW)*T2(K,IWD(JWW))
                      CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QSW(K,JWW)*C2(K,IWD(JWW),CN(1:NAC))
                    END DO
                    IF(QSUMM > 0.0)THEN
                    TTR(JTT)           = TSUM           /QSUMM
                    CTR(CN(1:NAC),JTT) = CSUM(CN(1:NAC))/QSUMM
                    ENDIF
                 ENDIF        
              ELSE
                !JSS(JBU)                 =  JSS(JBU)+1     ! SW 9/25/13
                KTSW(JSS(JBU),JBU)       =  KTPU(JP)
                KBSW(JSS(JBU),JBU)       =  KBPU(JP)
                JB                       =  JBU
                POINT_SINK(JSS(JBU),JBU) = .TRUE.
                ID                       =  IUPU(JP)
                QSTR(JSS(JBU),JBU)       =  QPU(JP)
                ESTR(JSS(JBU),JBU)       =  EPU(JP)
                KT                       =  KTWB(JWU)
                JW                       =  JWU
                CALL DOWNSTREAM_WITHDRAWAL (JSS(JBU))
                IF (IDPU(JP) /= 0 .AND. US(JBD) == IDPU(JP)) THEN
                  QSUMM = 0.0; TSUM  = 0.0; CSUM  = 0.0
                  DO K=KT,KB(ID)
                    QSUMM           = QSUMM          +QNEW(K)
                    TSUM            = TSUM           +QNEW(K)*T2(K,ID)
                    CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QNEW(K)*C2(K,ID,CN(1:NAC))
                  END DO
                  IF (QSUMM /= 0.0) THEN
                    TINSUM(JBD)           = (TSUM           +TINSUM(JBD)          *QINSUM(JBD))/(QSUMM+QINSUM(JBD))
                    CINSUM(CN(1:NAC),JBD) = (CSUM(CN(1:NAC))+CINSUM(CN(1:NAC),JBD)*QINSUM(JBD))/(QSUMM+QINSUM(JBD))
                    QINSUM(JBD)           =  QINSUM(JBD)    +QSUMM
                  END IF
                END IF
                QSUM(JB) = 0.0; TSUM = 0.0; CSUM = 0.0
                DO K=KT,KB(ID)
                  QSUM(JB)        = QSUM(JB)       +QOUT(K,JB)
                  TSUM            = TSUM           +QOUT(K,JB)*T2(K,ID)
                  CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QOUT(K,JB)*C2(K,ID,CN(1:NAC))
                END DO
                IF (QSUM(JB) /= 0.0) THEN
                  TOUT(JB)           = TSUM           /QSUM(JB)
                  COUT(CN(1:NAC),JB) = CSUM(CN(1:NAC))/QSUM(JB)
                END IF
                 IF (IDPU(JP) /= 0) THEN     ! SW 9/25/13 Moved code start
                IF (US(JBD) /= IDPU(JP) .OR. HEAD_FLOW(JBD) .OR. UP_HEAD(JBD)) THEN
                  JTT              = JTT+1
                  QTR(JTT)         = QPU(JP)
                  ITR(JTT)         = IDPU(JP)
                  PLACE_QTR(JTT)   = PPUC(JP) == ' DENSITY'
                  SPECIFY_QTR(JTT) = PPUC(JP) == ' SPECIFY'
                  IF (SPECIFY_QTR(JTT)) THEN
                    ELTRT(JTT) = ETPU(JP)
                    ELTRB(JTT) = EBPU(JP)
                  END IF
                  JBTR(JTT) = JBD
                    TTR(JTT)          = TOUT(JB)
                    CTR(CN(1:NAC),JTT)= COUT(CN(1:NAC),JB)
                  END IF
                ENDIF                     
              END IF
            END IF
            ENDIF
            
                
                
                ELSE
          IF (ELW <= EOFFPU(JP)) PUMPON(JP) = .FALSE.                                                       ! CB 1/13/06
          IF (ELW > EOFFPU(JP) .AND. QPU(JP) > 0.0) THEN
            IF (ELW >= EONPU(JP)) PUMPON(JP) = .TRUE.
            IF (PUMPON(JP)) THEN
              IF (LATERAL_PUMP(JP)) THEN
                JLAT      = 1
                !JWW       = JWW+1               ! SW 9/25/13
                !JBWD(JWW) = JBU  
                !IWD(JWW)  = IUPU(JP)
                QWD(JWW)  = QPU(JP)
                KTWD(JWW) = KTPU(JP)
                KBWD(JWW) = KBPU(JP)
                EWD(JWW)  = EPU(JP)
                I         = MAX(CUS(JBWD(JWW)),IWD(JWW))
                JB        = JBWD(JWW)
                JW        = JWU
                KT        = KTWB(JW)
                jwd=jww
                CALL LATERAL_WITHDRAWAL         ! (JWW)
                DO K=KTW(JWW),KBW(JWW)
                  QSS(K,I) = QSS(K,I)-QSW(K,JWW)
                END DO
                IF (IDPU(JP) /= 0) THEN           ! MOVED CODE SW 9/25/13
                  JTT              = JTT+1
                  QTR(JTT)         = QPU(JP)
                  ITR(JTT)         = IDPU(JP)
                  PLACE_QTR(JTT)   = PPUC(JP) == ' DENSITY'
                  SPECIFY_QTR(JTT) = PPUC(JP) == ' SPECIFY'
                  IF (SPECIFY_QTR(JTT)) THEN
                    ELTRT(JTT) = ETPU(JP)
                    ELTRB(JTT) = EBPU(JP)
                  END IF
                  JBTR(JTT) = JBD
                    TSUM = 0.0; QSUMM = 0.0; CSUM(CN(1:NAC)) = 0.0
                    DO K=KTW(JWW),KBW(JWW)
                      QSUMM           = QSUMM          +QSW(K,JWW)
                      TSUM            = TSUM           +QSW(K,JWW)*T2(K,IWD(JWW))
                      CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QSW(K,JWW)*C2(K,IWD(JWW),CN(1:NAC))
                    END DO
                    IF(QSUMM > 0.0)THEN
                    TTR(JTT)           = TSUM           /QSUMM
                    CTR(CN(1:NAC),JTT) = CSUM(CN(1:NAC))/QSUMM
                    ENDIF
                 ENDIF        ! SW 9/25/13 END MOVED CODE
              ELSE
                !JSS(JBU)                 =  JSS(JBU)+1     ! SW 9/25/13
                KTSW(JSS(JBU),JBU)       =  KTPU(JP)
                KBSW(JSS(JBU),JBU)       =  KBPU(JP)
                JB                       =  JBU
                POINT_SINK(JSS(JBU),JBU) = .TRUE.
                ID                       =  IUPU(JP)
                QSTR(JSS(JBU),JBU)       =  QPU(JP)
                ESTR(JSS(JBU),JBU)       =  EPU(JP)
                KT                       =  KTWB(JWU)
                JW                       =  JWU
                CALL DOWNSTREAM_WITHDRAWAL (JSS(JBU))
                IF (IDPU(JP) /= 0 .AND. US(JBD) == IDPU(JP)) THEN
                  QSUMM = 0.0; TSUM  = 0.0; CSUM  = 0.0
                  DO K=KT,KB(ID)
                    QSUMM           = QSUMM          +QNEW(K)
                    TSUM            = TSUM           +QNEW(K)*T2(K,ID)
                    CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QNEW(K)*C2(K,ID,CN(1:NAC))
                  END DO
                  IF (QSUMM /= 0.0) THEN
                    TINSUM(JBD)           = (TSUM           +TINSUM(JBD)          *QINSUM(JBD))/(QSUMM+QINSUM(JBD))
                    CINSUM(CN(1:NAC),JBD) = (CSUM(CN(1:NAC))+CINSUM(CN(1:NAC),JBD)*QINSUM(JBD))/(QSUMM+QINSUM(JBD))
                    QINSUM(JBD)           =  QINSUM(JBD)    +QSUMM
                  END IF
                END IF
                QSUM(JB) = 0.0; TSUM = 0.0; CSUM = 0.0
                DO K=KT,KB(ID)
                  QSUM(JB)        = QSUM(JB)       +QOUT(K,JB)
                  TSUM            = TSUM           +QOUT(K,JB)*T2(K,ID)
                  CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QOUT(K,JB)*C2(K,ID,CN(1:NAC))
                END DO
                IF (QSUM(JB) /= 0.0) THEN
                  TOUT(JB)           = TSUM           /QSUM(JB)
                  COUT(CN(1:NAC),JB) = CSUM(CN(1:NAC))/QSUM(JB)
                END IF
                 IF (IDPU(JP) /= 0) THEN     ! SW 9/25/13 Moved code start
                IF (US(JBD) /= IDPU(JP) .OR. HEAD_FLOW(JBD) .OR. UP_HEAD(JBD)) THEN
                  JTT              = JTT+1
                  QTR(JTT)         = QPU(JP)
                  ITR(JTT)         = IDPU(JP)
                  PLACE_QTR(JTT)   = PPUC(JP) == ' DENSITY'
                  SPECIFY_QTR(JTT) = PPUC(JP) == ' SPECIFY'
                  IF (SPECIFY_QTR(JTT)) THEN
                    ELTRT(JTT) = ETPU(JP)
                    ELTRB(JTT) = EBPU(JP)
                  END IF
                  JBTR(JTT) = JBD
                    TTR(JTT)          = TOUT(JB)
                    CTR(CN(1:NAC),JTT)= COUT(CN(1:NAC),JB)
                  END IF
                ENDIF                     ! Moved code end SW 9/25/13
              END IF
            END IF
          END IF
          ENDIF
        END IF
      END DO
    END IF
    IF (PIPES) THEN
      YSS   = YS
      VSS   = VS
      VSTS  = VST
      YSTS  = YST
      DTPS  = DTP
      QOLDS = QOLD
      CALL PIPE_FLOW        ! (NIT)
      DO JP=1,NPI
       
        if(dynpipe(jp) == '      ON')then                     ! SW 5/10/10
        qpi(jp)=qpi(jp)*bp(jp)
        endif


!****** Positive flows

        JLAT = 0
        JBU  = JBUPI(JP)
        JBD  = JBDPI(JP)
        tdgon=.false.        ! cb 1/16/13
        IF (QPI(JP) >= 0.0) THEN
          IF (LATERAL_PIPE(JP)) THEN
            JLAT      = 1
            JWW       = JWW+1
            IWD(JWW)  = IUPI(JP)
            QWD(JWW)  = QPI(JP)
            KTWD(JWW) = KTUPI(JP)
            KBWD(JWW) = KBUPI(JP)
            EWD(JWW)  = EUPI(JP)
            JBWD(JWW) = JBU
            I         = MAX(CUS(JBWD(JWW)),IWD(JWW))
            JB        = JBWD(JWW)
            JW        = JWUPI(JP)
            KT        = KTWB(JW)
            jwd=jww
            CALL LATERAL_WITHDRAWAL    !(JWW)
            DO K=KTW(JWW),KBW(JWW)
              QSS(K,I) = QSS(K,I)-QSW(K,JWW)
            END DO
          ELSE
            JSS(JBU)                 =  JSS(JBU)+1
            KTSW(JSS(JBU),JBU)       =  KTDPI(JP)
            KBSW(JSS(JBU),JBU)       =  KBDPI(JP)
            JB                       =  JBU
            POINT_SINK(JSS(JBU),JBU) = .TRUE.
            ID                       =  IUPI(JP)
            QSTR(JSS(JBU),JBU)       =  QPI(JP)
            ESTR(JSS(JBU),JBU)       =  EUPI(JP)
            KT                       =  KTWB(JWUPI(JP))
            JW                       =  JWUPI(JP)
            CALL DOWNSTREAM_WITHDRAWAL(JSS(JBU))
            IF (IDPI(JP) /= 0 .AND. US(JBD) == IDPI(JP)) THEN
              QSUMM = 0.0; TSUM  = 0.0; CSUM  = 0.0
              DO K=KT,KB(ID)
                QSUMM           = QSUMM          +QNEW(K)
                TSUM            = TSUM           +QNEW(K)*T2(K,ID)
                CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QNEW(K)*C2(K,ID,CN(1:NAC))
              END DO
              IF (QSUMM /= 0.0) THEN
                TINSUM(JBD)           = (TSUM           +QINSUM(JBD)*TINSUM(JBD))          /(QSUMM+QINSUM(JBD))
                CINSUM(CN(1:NAC),JBD) = (CSUM(CN(1:NAC))+QINSUM(JBD)*CINSUM(CN(1:NAC),JBD))/(QSUMM+QINSUM(JBD))
                QINSUM(JBD)           =  QINSUM(JBD)    +QSUMM
              END IF
            END IF
            QSUM(JB) = 0.0; TSUM = 0.0; CSUM = 0.0
            DO K=KT,KB(ID)
              QSUM(JB)        = QSUM(JB)       +QOUT(K,JB)
              TSUM            = TSUM           +QOUT(K,JB)*T2(K,ID)
              CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QOUT(K,JB)*C2(K,ID,CN(1:NAC))
            END DO
            IF (QSUM(JB) /= 0.0) THEN
              TOUT(JB)           = TSUM           /QSUM(JB)
              COUT(CN(1:NAC),JB) = CSUM(CN(1:NAC))/QSUM(JB)
            END IF
          END IF
          IF (IDPI(JP) /= 0) THEN
            IF (US(JBD) /= IDPI(JP) .OR. HEAD_FLOW(JBD) .OR. UP_HEAD(JBD)) THEN
              JTT              = JTT+1
              QTR(JTT)         = QPI(JP)
              ITR(JTT)         = IDPI(JP)
              PLACE_QTR(JTT)   = PDPIC(JP) == ' DENSITY'
              SPECIFY_QTR(JTT) = PDPIC(JP) == ' SPECIFY'
              IF (SPECIFY_QTR(JTT)) THEN
                ELTRT(JTT) = ETDPI(JP)
                ELTRB(JTT) = EBDPI(JP)
              END IF
              JBTR(JTT) = JBD
              IF (JLAT == 1) THEN
                TSUM = 0.0; QSUMM = 0.0; CSUM(CN(1:NAC)) = 0.0
                DO K=KTW(JWW),KBW(JWW)
                  QSUMM           = QSUMM          +QSW(K,JWW)
                  TSUM            = TSUM           +QSW(K,JWW)*T2(K,IWD(JWW))
                  CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QSW(K,JWW)*C2(K,IWD(JWW),CN(1:NAC))
                END DO
                TTR(JTT)           = TSUM           /QSUMM
                CTR(CN(1:NAC),JTT) = CSUM(CN(1:NAC))/QSUMM
              ELSE
                TTR(JTT)           = TOUT(JB)
                CTR(CN(1:NAC),JTT) = COUT(CN(1:NAC),JB)
              END IF
            ELSE
              IF (LATERAL_PIPE(JP)) THEN
                TSUM      = 0.0; QSUMM = 0.0; CSUM = 0.0
                ILAT(JWW) = 1
                JB        = JBD
                DO K=KTW(JWW),KBW(JWW)
                  QSUMM           = QSUMM          +QSW(K,JWW)
                  TSUM            = TSUM           +QSW(K,JWW)*T2(K,IWD(JWW))
                  CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QSW(K,JWW)*C2(K,IWD(JWW),CN(1:NAC))
                END DO
                TINSUM(JB)           = (TINSUM(JB)          *QINSUM(JB)+TSUM)           /(QSUMM+QINSUM(JB))
                CINSUM(CN(1:NAC),JB) = (CINSUM(CN(1:NAC),JB)*QINSUM(JB)+CSUM(CN(1:NAC)))/(QSUMM+QINSUM(JB))
                QINSUM(JB)           =  QSUMM               +QINSUM(JB)
              END IF
            END IF
          END IF
        ELSE
          JTT              =  JTT+1
          JWW              =  JWW+1
          IWD(JWW)         =  IDPI(JP)
          ITR(JTT)         =  IUPI(JP)
          QTR(JTT)         = -QPI(JP)
          QWD(JWW)         = -QPI(JP)
          KTWD(JWW)        =  KTDPI(JP)
          KBWD(JWW)        =  KBDPI(JP)
          EWD(JWW)         =  EDPI(JP)
          PLACE_QTR(JTT)   =  PUPIC(JP) == ' DENSITY'
          SPECIFY_QTR(JTT) =  PUPIC(JP) == ' SPECIFY'
          IF (SPECIFY_QTR(JTT)) THEN
            ELTRT(JTT) = ETUPI(JP)
            ELTRB(JTT) = EBUPI(JP)
          END IF
          JBTR(JTT) = JBU
          JBWD(JWW) = JBD
          I         = MAX(CUS(JBWD(JWW)),IWD(JWW))
          JB        = JBWD(JWW)
          JW        = JWDPI(JP)
          KT        = KTWB(JW)
          jwd=jww
          CALL LATERAL_WITHDRAWAL   !(JWW)
          DO K=KTW(JWW),KBW(JWW)
            QSS(K,I) = QSS(K,I)-QSW(K,JWW)
          END DO
          IF (IDPI(JP) /= 0) THEN
            TSUM  = 0.0; QSUMM = 0.0; CSUM  = 0.0
            DO K=KTW(JWW),KBW(JWW)
              QSUMM           = QSUMM          +QSW(K,JWW)
              TSUM            = TSUM           +QSW(K,JWW)*T2(K,IWD(JWW))
              CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QSW(K,JWW)*C2(K,IWD(JWW),CN(1:NAC))
            END DO
            TTR(JTT)           = TSUM           /QSUMM
            CTR(CN(1:NAC),JTT) = CSUM(CN(1:NAC))/QSUMM
          END IF
        END IF
      END DO
    END IF
    IF (GATES) THEN
      CALL GATE_FLOW
      IF(TDGTA) CALL TDGtarget                   ! tdgtarget 
      IF (SYSTDG) CALL SYSTDG_TDG  ! SYSTDG - CALCULATE TDG_TDG (%)
      DO JG=1,NGT

!****** Positive flows

        JLAT = 0
        JBU  = JBUGT(JG)
        JBD  = JBDGT(JG)
        tdgon=.false.        ! cb 1/16/13
        jsg=jg
        nnsg=1
        if (cac(ndo) == '      ON' .and. gasgtc(jg) == '      ON')tdgon=.true.
        IF (QGT(JG) >= 0.0) THEN
          IF (LATERAL_GATE(JG)) THEN
 !           JLAT      = 1
            JWW       = JWW+1
            IWD(JWW)  = IUGT(JG)
            QWD(JWW)  = QGT(JG)
            KTWD(JWW) = KTUGT(JG)
            KBWD(JWW) = KBUGT(JG)
            EWD(JWW)  = EGT(JG)
            IF(DYNGTC(JG) == '     ZGT' .AND. GT2CHAR == 'EGT2ELEV')THEN                    ! SW 2/25/11
              IF(EGT2(JG) /=  0.0)THEN
              EWD(JWW) = EGT2(JG)
              ENDIF
            ENDIF
            JBWD(JWW) = JBU
            I         = MAX(CUS(JBWD(JWW)),IWD(JWW))
            JW        = JWUGT(JG)
            JB        = JBWD(JWW)
            KT        = KTWB(JW)
            jwd=jww
            CALL LATERAL_WITHDRAWAL    !(JWW)
            DO K=KTW(JWW),KBW(JWW)
              QSS(K,I) = QSS(K,I)-QSW(K,JWW)
            END DO
            IF (IDGT(JG) /= 0) THEN
              CSUM(CN(1:NAC)) = 0.0; TSUM = 0.0; QSUMM = 0.0
              JTT              = JTT+1                     !  SW 4/1/09
              QTR(JTT)         = QGT(JG)                   !  SW 4/1/09
              ITR(JTT)         = IDGT(JG)                  !  SW 4/1/09
              PLACE_QTR(JTT)   = PDGTC(JG) == ' DENSITY'   !  SW 4/1/09
              SPECIFY_QTR(JTT) = PDGTC(JG) == ' SPECIFY'   !  SW 4/1/09
              IF (SPECIFY_QTR(JTT)) THEN                   !  SW 4/1/09
                ELTRT(JTT) = ETDGT(JG)                     !  SW 4/1/09
                ELTRB(JTT) = EBDGT(JG)                     !  SW 4/1/09
              END IF
              JBTR(JTT) = JBD                              !  SW 4/1/09
              DO K=KTW(JWW),KBW(JWW)
                QSUMM           = QSUMM          +QSW(K,JWW)
                TSUM            = TSUM           +QSW(K,JWW)*T2(K,IWD(JWW))
                CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QSW(K,JWW)*C2(K,IWD(JWW),CN(1:NAC))
              END DO
              IF(QSUMM==0.0)THEN
              TTR(JTT)=0.0
              CTR(:,JTT)=0.0
              ELSE
              TTR(JTT) = TSUM/QSUMM
              DO JC=1,NAC
                CTR(CN(JC),JTT) = CSUM(CN(JC))/QSUMM
                IF (CN(JC) == NDO .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN
                  TDG_GATE(JWW,JG) = .TRUE.
                  !
                  ! systdg 
                  IF (SYSTDG) THEN 
                    IF(GTNAME(JG)) THEN                       
                      CALL  UPDATE_TDGC(0,PALT(ID),JG,TTR(JTT),CTR(CN(JC),JTT))   ! O2        
                    ELSE
                  CALL TOTAL_DISSOLVED_GAS(0,PALT(ID),1,JG,TTR(JTT),CTR(CN(JC),JTT))   ! O2
                END IF
                  ELSE
                      CALL TOTAL_DISSOLVED_GAS(0,PALT(ID),1,JG,TTR(JTT),CTR(CN(JC),JTT))   ! O2
                  END IF
                  !
                END IF
                IF (CN(JC) == NN2 .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN
                  TDG_GATE(JWW,JG) = .TRUE.
                  !
                  ! systdg 
                  IF (SYSTDG) THEN                                 
                    IF (GTNAME(JG)) THEN                       
                      CALL  UPDATE_TDGC(1,PALT(ID),JG,TTR(JTT),CTR(CN(JC),JTT))   ! N2         
                    ELSE
                      CALL TOTAL_DISSOLVED_GAS(1,PALT(ID),1,JG,TTR(JTT),CTR(CN(JC),JTT))   ! N2
                    END IF
                  ELSE 
                  CALL TOTAL_DISSOLVED_GAS(1,PALT(ID),1,JG,TTR(JTT),CTR(CN(JC),JTT))   ! N2
                END IF
                  !
                END IF
                IF (CN(JC) == NDGP .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN   !8/2020 TDGP
                  TDG_GATE(JWW,JG) = .TRUE.
                  !
                  IF (SYSTDG) THEN                                 
                    IF (GTNAME(JG)) THEN                       
                      CALL  UPDATE_TDGC(2,PALT(ID),JG,TTR(JTT),CTR(CN(JC),JTT))          
                    ELSE
                      CALL TOTAL_DISSOLVED_GAS(2,PALT(ID),1,JG,TTR(JTT),CTR(CN(JC),JTT))   
                    END IF
                  ELSE 
                  CALL TOTAL_DISSOLVED_GAS(2,PALT(ID),1,JG,TTR(JTT),CTR(CN(JC),JTT))   
                END IF
                  !
                END IF
              END DO
              ENDIF
            ELSE IF (CAC(NDO) == '      ON' .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN    
              TDG_GATE(JWW,JG) = .TRUE.
            END IF
          ELSE
            JSS(JBU)                 =  JSS(JBU)+1
            KTSW(JSS(JBU),JBU)       =  KTUGT(JG)
            KBSW(JSS(JBU),JBU)       =  KBUGT(JG)
            JB                       =  JBU
            POINT_SINK(JSS(JBU),JBU) = .TRUE.
            ID                       =  IUGT(JG)
            ESTR(JSS(JBU),JBU)       =  EGT(JG)
            IF(DYNGTC(JG) == '     ZGT' .AND. GT2CHAR == 'EGT2ELEV')THEN                    ! SW 2/25/11
              IF(EGT2(JG) /=  0.0)THEN
              ESTR(JSS(JBU),JBU)       =  EGT2(JG)
              ENDIF
            ENDIF
            QSTR(JSS(JBU),JBU)       =  QGT(JG)
            KT                       =  KTWB(JWUGT(JG))
            JW                       =  JWUGT(JG)
            CALL DOWNSTREAM_WITHDRAWAL (JSS(JBU))
            QSUM(JB) = 0.0; TSUM = 0.0; CSUM = 0.0
            DO K=KT,KB(ID)
              QSUM(JB) = QSUM(JB)+QOUT(K,JB)
              TSUM     = TSUM+QOUT(K,JB)*T2(K,ID)
              DO JC=1,NAC
                IF (CN(JC) == NDO .AND. CAC(NDO) == '      ON' .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN   ! MM 5/21/2009
                  T2R4=T2(K,ID)
                  CGAS=C2(K,ID,CN(JC))                                                                                    ! MM 5/21/2009                  
                  !
                  ! systdg 
                  IF (SYSTDG) THEN                      
                    IF(GTNAME(JG)) THEN                
                      CALL  UPDATE_TDGC(0,PALT(ID),JG,T2R4,CGAS)                       
                    ELSE
                      CALL TOTAL_DISSOLVED_GAS(0,PALT(ID),1,JG,T2R4,CGAS)
                    END IF
                  ELSE                                                                                    ! MM 5/21/2009                  
                  CALL TOTAL_DISSOLVED_GAS(0,PALT(ID),1,JG,T2R4,CGAS)
                  END IF
                  !
                  CSUM(CN(JC)) = CSUM(CN(JC))+QOUT(K,JB)*CGAS
                 !ELSEIF (CN(JC) == NGN2 .AND. CAC(NGN2) == '      ON' .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN   ! SW 10/27/15
                ELSEIF (CN(JC) == NN2 .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN   ! SW 10/27/15
                  if(CAC(NN2) == '      ON')then                                                                                ! cb 1/13/16
                      T2R4=T2(K,ID)
                      CGAS=C2(K,ID,CN(JC))   
                      !
                      ! systdg 
                      IF (SYSTDG) THEN                    
                        IF(GTNAME(JG)) THEN              
                          CALL  UPDATE_TDGC(1,PALT(ID),JG,T2R4,CGAS)               
                        ELSE
                          CALL TOTAL_DISSOLVED_GAS(1,PALT(ID),1,JG,T2R4,CGAS)
                        END IF
                      ELSE                                                                                                   
                      CALL TOTAL_DISSOLVED_GAS(1,PALT(ID),1,JG,T2R4,CGAS)
                      END IF
                      !
                      CSUM(CN(JC)) = CSUM(CN(JC))+QOUT(K,JB)*CGAS
                    end if
                ELSEIF (CN(JC)==NDGP .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN   !8/2020 TDGP
                  if(CAC(NDGP) == '      ON')then                                                                              
                      T2R4=T2(K,ID)
                      CGAS=C2(K,ID,CN(JC))   
                      !               
                      IF (SYSTDG) THEN                    
                        IF(GTNAME(JG)) THEN              
                          CALL  UPDATE_TDGC(2,PALT(ID),JG,T2R4,CGAS)               
                        ELSE
                          CALL TOTAL_DISSOLVED_GAS(2,PALT(ID),1,JG,T2R4,CGAS)
                        END IF
                      ELSE                                                                                                   
                      CALL TOTAL_DISSOLVED_GAS(2,PALT(ID),1,JG,T2R4,CGAS)
                      END IF
                      !
                      CSUM(CN(JC)) = CSUM(CN(JC))+QOUT(K,JB)*CGAS
                    end if  
                ELSE
                  CSUM(CN(JC)) = CSUM(CN(JC))+QOUT(K,JB)*C2(K,ID,CN(JC))
                END IF
              END DO
            END DO
            IF (QSUM(JB) /= 0.0) THEN
              TOUT(JB)           = TSUM           /QSUM(JB)
              COUT(CN(1:NAC),JB) = CSUM(CN(1:NAC))/QSUM(JB)
            END IF
            IF (IDGT(JG) /= 0 .AND. US(JBD) == IDGT(JG)) THEN
              QSUMM = 0.0
              TSUM  = 0.0
              CSUM  = 0.0
              DO K=KT,KB(ID)
                QSUMM = QSUMM+QNEW(K)
                TSUM  = TSUM+QNEW(K)*T2(K,ID)
                DO JC=1,NAC
                  IF (CN(JC) == NDO .AND. CAC(NDO) == '      ON' .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN   ! MM 5/21/2009
                    T2R4=T2(K,ID)
                    CGAS=C2(K,ID,CN(JC))                                                                                    ! MM 5/21/2009
                    !
                    ! systdg 
                    IF (SYSTDG) THEN                    
                        IF (GTNAME(JG)) THEN          
                          CALL  UPDATE_TDGC(0,PALT(ID),JG,T2R4,CGAS)   ! O2       
                        ELSE
                    CALL TOTAL_DISSOLVED_GAS(0,PALT(ID),1,JG,T2R4,CGAS)   ! O2
                        END IF
                    ELSE                                                                                   ! MM 5/21/2009
                        CALL TOTAL_DISSOLVED_GAS(0,PALT(ID),1,JG,T2R4,CGAS)   ! O2
                    END IF
                    !
                    CSUM(CN(JC)) = CSUM(CN(JC))+QNEW(K)*CGAS
                  !ELSEIF (CN(JC) == NGN2 .AND. CAC(NGN2) == '      ON' .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN   ! SW 10/27/15
                  ELSEIF (CN(JC) == NN2 .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN   ! SW 10/27/15
                    if(CAC(NN2) == '      ON')then
                      T2R4=T2(K,ID)
                      CGAS=C2(K,ID,CN(JC)) 
                      !
                      ! systdg 
                       IF (SYSTDG) THEN                 
                        IF(GTNAME(JG)) THEN           
                          CALL  UPDATE_TDGC(1,PALT(ID),JG,T2R4,CGAS)   ! N2           
                        ELSE
                          CALL TOTAL_DISSOLVED_GAS(1,PALT(ID),1,JG,T2R4,CGAS)   ! N2
                        END IF
                      ELSE                                                                               
                      CALL TOTAL_DISSOLVED_GAS(1,PALT(ID),1,JG,T2R4,CGAS)   ! N2
                      END IF
                      !
                      CSUM(CN(JC)) = CSUM(CN(JC))+QNEW(K)*CGAS
                    END IF
                  ELSEIF (CN(JC)==NDGP .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN   !8/2020 TDGP
                    if(CAC(NDGP) == '      ON')then
                      T2R4=T2(K,ID)
                      CGAS=C2(K,ID,CN(JC)) 
                      !
                       IF (SYSTDG) THEN                 
                        IF(GTNAME(JG)) THEN           
                          CALL  UPDATE_TDGC(2,PALT(ID),JG,T2R4,CGAS)            
                        ELSE
                          CALL TOTAL_DISSOLVED_GAS(2,PALT(ID),1,JG,T2R4,CGAS)   
                        END IF
                      ELSE                                                                               
                      CALL TOTAL_DISSOLVED_GAS(2,PALT(ID),1,JG,T2R4,CGAS) 
                      END IF
                      !
                      CSUM(CN(JC)) = CSUM(CN(JC))+QNEW(K)*CGAS
                    END IF
                  ELSE
                    CSUM(CN(JC)) = CSUM(CN(JC))+QNEW(K)*C2(K,ID,CN(JC))
                  END IF
                END DO
              END DO
              IF (QSUMM /= 0.0) THEN
                TINSUM(JBD)           = (TSUM           +QINSUM(JBD)*TINSUM(JBD))          /(QSUMM+QINSUM(JBD))
                CINSUM(CN(1:NAC),JBD) = (CSUM(CN(1:NAC))+QINSUM(JBD)*CINSUM(CN(1:NAC),JBD))/(QSUMM+QINSUM(JBD))
                QINSUM(JBD)           =  QINSUM(JBD)    +QSUMM
              END IF
           ELSEIF (IDGT(JG) /= 0) THEN
              JTT              = JTT+1
              QTR(JTT)         = QGT(JG)
              ITR(JTT)         = IDGT(JG)
              PLACE_QTR(JTT)   = PDGTC(JG) == ' DENSITY'
              SPECIFY_QTR(JTT) = PDGTC(JG) == ' SPECIFY'
              IF (SPECIFY_QTR(JTT)) THEN
                ELTRT(JTT) = ETDGT(JG)
                ELTRB(JTT) = EBDGT(JG)
              END IF
              JBTR(JTT) = JBD
                TTR(JTT) =  TOUT(JB)
                DO JC=1,NAC
                  CTR(CN(JC),JTT) = COUT(CN(JC),JB)
                  IF (CN(JC) == NDO .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN
                    !
                    ! systdg 
                    IF (SYSTDG) THEN                    
                        IF (GTNAME(JG)) THEN          
                          CALL  UPDATE_TDGC (0,PALT(ID),JG,TTR(JTT),CTR(CN(JC),JTT))    ! O2          
                        ELSE
                          CALL TOTAL_DISSOLVED_GAS (0,PALT(ID),0,JS,TTR(JTT),CTR(CN(JC),JTT))    ! O2
                        END IF
                      ELSE  
                    CALL TOTAL_DISSOLVED_GAS (0,PALT(ID),0,JS,TTR(JTT),CTR(CN(JC),JTT))    ! O2
                  END IF
                      !
                  END IF
                  IF (CN(JC) == NN2 .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN
                    ! 
                    ! systdg 
                    IF (SYSTDG) THEN                     
                        IF (GTNAME(JG)) THEN          
                          CALL  UPDATE_TDGC(1,PALT(ID),JG,TTR(JTT),CTR(CN(JC),JTT))    ! N2       
                        ELSE
                    CALL TOTAL_DISSOLVED_GAS(1,PALT(ID),0,JS,TTR(JTT),CTR(CN(JC),JTT))    ! N2
                  END IF
                      ELSE
                          CALL TOTAL_DISSOLVED_GAS(1,PALT(ID),0,JS,TTR(JTT),CTR(CN(JC),JTT))    ! N2
                      END IF
                      !
                  END IF
                  IF (CN(JC)==NDGP .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN    !8/2020 TDGP
                    ! 
                    IF (SYSTDG) THEN                     
                      IF (GTNAME(JG)) THEN          
                        CALL  UPDATE_TDGC(2,PALT(ID),JG,TTR(JTT),CTR(CN(JC),JTT))        
                        ELSE
                        CALL TOTAL_DISSOLVED_GAS(2,PALT(ID),0,JS,TTR(JTT),CTR(CN(JC),JTT))    
                  END IF
                      ELSE
                      CALL TOTAL_DISSOLVED_GAS(2,PALT(ID),0,JS,TTR(JTT),CTR(CN(JC),JTT))    
                      END IF
                      !
                  END IF
                END DO
            END IF
          END IF
        ELSE IF (QGT(JG) < 0.0) THEN
          JTT              =  JTT+1
          JWW              =  JWW+1
          IWD(JWW)         =  IDGT(JG)
          ITR(JTT)         =  IUGT(JG)
          QTR(JTT)         = -QGT(JG)
          QWD(JWW)         = -QGT(JG)
          KTWD(JWW)        =  KTDGT(JG)
          KBWD(JWW)        =  KBDGT(JG)
          EWD(JWW)         =  EGT(JG)
          PLACE_QTR(JTT)   =  PUGTC(JG) == ' DENSITY'
          SPECIFY_QTR(JTT) =  PUGTC(JG) == ' SPECIFY'
          IF (SPECIFY_QTR(JTT)) THEN
            ELTRT(JTT) = ETUGT(JG)
            ELTRB(JTT) = EBUGT(JG)
          END IF
          JBTR(JTT) = JBU
          JBWD(JWW) = JBD
          I         = MAX(CUS(JBWD(JWW)),IWD(JWW))
          JW        = JWDGT(JG)
          JB        = JBWD(JWW)
          KT        = KTWB(JW)
          jwd=jww
          CALL LATERAL_WITHDRAWAL !(JWW)
          DO K=KTW(JWW),KBW(JWW)
            QSS(K,I) = QSS(K,I)-QSW(K,JWW)
          END DO
          IF (IDGT(JG) /= 0) THEN
            CSUM(CN(1:NAC)) = 0.0; TSUM = 0.0; QSUMM = 0.0
            DO K=KTW(JWW),KBW(JWW)
              QSUMM           = QSUMM          +QSW(K,JWW)
              TSUM            = TSUM           +QSW(K,JWW)*T2(K,IWD(JWW))
              CSUM(CN(1:NAC)) = CSUM(CN(1:NAC))+QSW(K,JWW)*C2(K,IWD(JWW),CN(1:NAC))
            END DO
            TTR(JTT) = TSUM/QSUMM
            DO JC=1,NAC
              CTR(CN(JC),JTT) = CSUM(CN(JC))/QSUMM
              IF (CN(JC) == NDO .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN
                TDG_GATE(JWW,JG) = .TRUE.
                !
                ! systdg 
                 IF (SYSTDG) THEN                  
                     IF (GTNAME(JG)) THEN         
                     CALL  UPDATE_TDGC(0,PALT(ID),JG,TTR(JTT),CTR(CN(JC),JTT))            
                     ELSE
                CALL TOTAL_DISSOLVED_GAS(0,PALT(ID),1,JG,TTR(JTT),CTR(CN(JC),JTT))
              END IF
                 ELSE
                     CALL TOTAL_DISSOLVED_GAS(0,PALT(ID),1,JG,TTR(JTT),CTR(CN(JC),JTT))
                 END IF
                 !
              END IF
              IF (CN(JC) == NN2 .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN
                TDG_GATE(JWW,JG) = .TRUE.
                !
                ! systdg 
                 IF (SYSTDG) THEN                
                    IF (GTNAME(JG)) THEN         
                     CALL  UPDATE_TDGC(1,PALT(ID),JG,TTR(JTT),CTR(CN(JC),JTT))           
                    ELSE
                     CALL TOTAL_DISSOLVED_GAS(1,PALT(ID),1,JG,TTR(JTT),CTR(CN(JC),JTT))
                    END IF
                 ELSE
                CALL TOTAL_DISSOLVED_GAS(1,PALT(ID),1,JG,TTR(JTT),CTR(CN(JC),JTT))
              END IF
                 !
              END IF
              IF (CN(JC)==NDGP .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN      !8/2020 TDGP
                TDG_GATE(JWW,JG) = .TRUE.
                !
                IF (SYSTDG) THEN                
                  IF (GTNAME(JG)) THEN         
                    CALL UPDATE_TDGC(2,PALT(ID),JG,TTR(JTT),CTR(CN(JC),JTT))           
                  ELSE
                    CALL TOTAL_DISSOLVED_GAS(2,PALT(ID),1,JG,TTR(JTT),CTR(CN(JC),JTT))
                  END IF
                ELSE
                  CALL TOTAL_DISSOLVED_GAS(2,PALT(ID),1,JG,TTR(JTT),CTR(CN(JC),JTT))
              END IF
                 !
              END IF
            END DO
          ELSE IF (CAC(NDO) == '      ON' .AND. GASGTC(JG) == '      ON' .AND. QGT(JG) > 0.0) THEN
            TDG_GATE(JWW,JG) = .TRUE.
          END IF
        END IF
      END DO
      ! systdg - Add power house release tdg update
      IF (SYSTDG) THEN
         IF (POWNO>0 .AND. TDGLOC=='     REL') THEN
             DO ip = 1, POWNO
                 IF (QGT(POWGTNO(ip))>0.0 .AND. TDG_ROSP>0.0) THEN
                   IF (LATERAL_GATE(POWGTNO(ip))) THEN
                      IF(CONSTITUENTS) THEN
                         CALL UPDATE_TDGC(0, palt(IUGT(POWGTNO(ip))), POWGTNO(ip), tavgw(POWGTNO(ip)), CAVGW(POWGTNO(ip),NDO))
                         CALL UPDATE_TDGC(1, palt(IUGT(POWGTNO(ip))), POWGTNO(ip), tavgw(POWGTNO(ip)), CAVGW(POWGTNO(ip),NN2))
                         CALL UPDATE_TDGC(2, palt(IUGT(POWGTNO(ip))), POWGTNO(ip), tavgw(POWGTNO(ip)), CAVGW(POWGTNO(ip),NDGP))
                      END IF
                      IF(DERIVED_CALC) THEN
                         CDAVGW(POWGTNO(ip),O2DG_DER)=tdg_tdg
                         CDAVGW(POWGTNO(ip),TDG_DER)=tdg_tdg
                      END IF
                   ELSE
                      IF(CONSTITUENTS) THEN
                         CALL UPDATE_TDGC(0,PALT(IUGT(POWGTNO(ip))), POWGTNO(ip), Tavg(POWGTNO(ip),JBUGT(POWGTNO(ip))), CAVG(POWGTNO(ip),JBUGT(POWGTNO(ip)),NDO))
                         CALL UPDATE_TDGC(1,PALT(IUGT(POWGTNO(ip))), POWGTNO(ip), Tavg(POWGTNO(ip),JBUGT(POWGTNO(ip))), CAVG(POWGTNO(ip),JBUGT(POWGTNO(ip)),NN2))
                         CALL UPDATE_TDGC(2,PALT(IUGT(POWGTNO(ip))), POWGTNO(ip), Tavg(POWGTNO(ip),JBUGT(POWGTNO(ip))), CAVG(POWGTNO(ip),JBUGT(POWGTNO(ip)),NDGP))
                      END IF
                      IF(DERIVED_CALC) THEN
                         CDAVG(POWGTNO(ip),JBUGT(POWGTNO(ip)),O2DG_DER)=tdg_tdg
                         CDAVG(POWGTNO(ip),JBUGT(POWGTNO(ip)),TDG_DER)=tdg_tdg
                      END IF
                   END IF
                 ELSE IF(QGT(POWGTNO(ip))<0.0 .AND. TDG_ROSP>0.0) THEN! QGT<0.0
                      IF(CONSTITUENTS) THEN
                         CALL UPDATE_TDGC(0, palt(IUGT(POWGTNO(ip))), POWGTNO(ip), tavgw(POWGTNO(ip)), CAVGW(POWGTNO(ip),NDO))
                         CALL UPDATE_TDGC(1, palt(IUGT(POWGTNO(ip))), POWGTNO(ip), tavgw(POWGTNO(ip)), CAVGW(POWGTNO(ip),NN2))
                         CALL UPDATE_TDGC(2, palt(IUGT(POWGTNO(ip))), POWGTNO(ip), tavgw(POWGTNO(ip)), CAVGW(POWGTNO(ip),NDGP))
                      END IF
                      IF(DERIVED_CALC) THEN
                         CDAVGW(POWGTNO(ip),O2DG_DER)=tdg_tdg
                         CDAVGW(POWGTNO(ip),TDG_DER)=tdg_tdg
                      END IF
                 END IF
              END DO
          END IF
          IF(POWNO>0 .AND. TDGLOC=='     REL') THEN
          DO ifl = 1, FLNO
            IF (QGT(FLGTNO(ifl))>0.0 .AND. TDG_ROSP>0.0) THEN
              IF(CONSTITUENTS) THEN
                     CALL UPDATE_TDGC(0,PALT(IUGT(FLGTNO(ifl))), FLGTNO(ifl), Tavg(FLGTNO(ifl),JBUGT(FLGTNO(ifl))), CAVG(FLGTNO(ifl),JBUGT(FLGTNO(ifl)),NDO))
                     CALL UPDATE_TDGC(1,PALT(IUGT(FLGTNO(ifl))), FLGTNO(ifl), Tavg(FLGTNO(ifl),JBUGT(FLGTNO(ifl))), CAVG(FLGTNO(ifl),JBUGT(FLGTNO(ifl)),NN2))
                     CALL UPDATE_TDGC(2,PALT(IUGT(FLGTNO(ifl))), FLGTNO(ifl), Tavg(FLGTNO(ifl),JBUGT(FLGTNO(ifl))), CAVG(FLGTNO(ifl),JBUGT(FLGTNO(ifl)),NDGP))
                  END IF
                  IF(DERIVED_CALC) THEN
                     CDAVG(FLGTNO(ifl),JBUGT(FLGTNO(ifl)),O2DG_DER)=tdg_tdg
                     CDAVG(FLGTNO(ifl),JBUGT(FLGTNO(ifl)),TDG_DER)=tdg_tdg
                  END IF 
            ELSE IF (QGT(FLGTNO(ifl))<0.0 .AND. TDG_ROSP>0.0) THEN
                 IF(CONSTITUENTS) THEN
                     CALL UPDATE_TDGC(0, palt(IUGT(FLGTNO(ifl))), FLGTNO(ifl), tavgw(FLGTNO(ifl)), CAVGW(FLGTNO(ifl),NDO))
                     CALL UPDATE_TDGC(1, palt(IUGT(FLGTNO(ifl))), FLGTNO(ifl), tavgw(FLGTNO(ifl)), CAVGW(FLGTNO(ifl),NN2))
                     CALL UPDATE_TDGC(2, palt(IUGT(FLGTNO(ifl))), FLGTNO(ifl), tavgw(FLGTNO(ifl)), CAVGW(FLGTNO(ifl),NDGP))
                  END IF
                  IF(DERIVED_CALC) THEN
                     CDAVGW(FLGTNO(ifl),O2DG_DER)=tdg_tdg
                     CDAVGW(FLGTNO(ifl),TDG_DER)=tdg_tdg
                  END IF
            END IF
          END DO
          END IF
      END IF
      ! systdg - Add power house release tdg update
    END IF
    
    tdgon=.false.                         ! cb 1/17/13
    tributaries = jtt > 0
    withdrawals = jww > 0
    
    
    DO JW=1,NWB
        DO JB=BS(JW),BE(JW)
            IF(BR_INACTIVE(JB))THEN
                ! CONVERT INFLOWS TO TRIBS SET TO THE CUS(1) LOCATION
                JTT=JTT+1
                !ITR(JTT)         =  CUS(1)   ! HARDWIRED TO FIRST BRANCH
                ITR(JTT)         =  CUS(jbdn(jw))   ! changed HARDWIRE TO JBDN BRANCH  ! cb 11/20/19
                  QTR(JTT)         = QIN(JB)
                  TTR(JTT)         = TIN(JB)
                    DO JC=1,NAC
                    CTR(CN(JC),JTT) = CIN(CN(JC),JB)
                    END DO
               !  PLACE_QTR(JTT)   =  '   DISTR'
                  PLACE_QTR(JTT)   =  .FALSE.         !SR 01/22/2018
                  JBTR(JTT) = 1
            ENDIF        
        ENDDO
    ENDDO
    
    
    
    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
      IF(BR_INACTIVE(JB))CYCLE  ! SW 6/12/2017
        IU = CUS(JB)
        ID = DS(JB)
        IF (EVAPORATION(JW)) THEN
          EVBR(JB) = 0.0
          DO I=IU,ID
            FW = AFW(JW)+BFW(JW)*WIND2(I)**CFW(JW)
            IF (RH_EVAP(JW)) THEN
              EA = EXP(2.3026*(7.5*TDEW(JW)/(TDEW(JW)+237.3)+0.6609))
              ES = EXP(2.3026*(7.5*T2(KT,I)/(T2(KT,I)+237.3)+0.6609))
              IF (TDEW(JW) < 0.0) EA = EXP(2.3026*(9.5*TDEW(JW)/(TDEW(JW)+265.5)+0.6609))
              IF (T2(KT,I) < 0.0) ES = EXP(2.3026*(9.5*T2(KT,I)/(T2(KT,I)+265.5)+0.6609))
              TAIRV = (TAIR(JW)+273.0)/(1.0-0.378*EA/760.0)
              DTV   = (T2(KT,I)+273.0)/(1.0-0.378*ES/760.0)-TAIRV
              DTVL  =  0.0084*WIND2(I)**3
              IF (DTV < DTVL) DTV = DTVL
              FW = (3.59*DTV**0.3333333+4.26*WIND2(I))
            END IF
            TM    = (T2(KT,I)+TDEW(JW))*0.5
            VPTG  =  0.35+0.015*TM+0.0012*TM*TM
            EV(I) =  VPTG*(T2(KT,I)-TDEW(JW))*FW*BI(KT,I)*DLX(I)/2.45E9
            IF (EV(I) < 0.0 .OR. ICE(I)) EV(I) = 0.0
            QSS(KT,I) = QSS(KT,I)-EV(I)
            EVBR(JB)  = EVBR(JB)+EV(I)
          END DO
        END IF
        IF (PRECIPITATION(JW)) THEN
          QPRBR(JB) = 0.0
          DO I=IU,ID
            QPR(I)    = PR(JB)*BI(KT,I)*DLX(I)
            QPRBR(JB) = QPRBR(JB)+QPR(I)
            QSS(KT,I) = QSS(KT,I)+QPR(I)
          END DO
        END IF
        IF (TRIBUTARIES) THEN
          DO JT=1,JTT

!********** Inflow fractions

            IF (JB == JBTR(JT)) THEN
              I = MAX(ITR(JT),IU)
              QTRF(KT:KB(I),JT) = 0.0
              IF (PLACE_QTR(JT)) THEN

!************** Inflow layer

                SSTOT = 0.0
                DO J=NSSS,NSSE
                  SSTOT = SSTOT+CTR(J,JT)
                END DO
                RHOTR = DENSITY(TTR(JT),CTR(NTDS,JT),SSTOT)
                K     = KT
                DO WHILE (RHOTR > RHO(K,I) .AND. K < KB(I))
                  K = K+1
                END DO
                KTTR(JT) = K
                KBTR(JT) = K

!************** Layer inflows

                VQTR  =  QTR(JT)*DLT
                VQTRI =  VQTR
                QTRFR =  1.0
                INCR  = -1
                DO WHILE (QTRFR > 0.0)
                  IF (K <= KB(I)) THEN
                    V1 = VOL(K,I)
                    IF (VQTR > 0.5*V1) THEN
                      QTRF(K,JT) = 0.5*V1/VQTRI
                      QTRFR      = QTRFR-QTRF(K,JT)
                      VQTR       = VQTR-QTRF(K,JT)*VQTRI
                      IF (K == KT) THEN
                        K    = KBTR(JT)
                        INCR = 1
                      END IF
                    ELSE
                      QTRF(K,JT) = QTRFR
                      QTRFR      = 0.0
                    END IF
                    IF (INCR < 0) KTTR(JT) = K
                    IF (INCR > 0) KBTR(JT) = MIN(KB(I),K)
                    K = K+INCR
                  ELSE
                    QTRF(KT,JT) = QTRF(KT,JT)+QTRFR
                    QTRFR       = 0.0
                  END IF
                END DO
              ELSE
                IF (SPECIFY_QTR(JT)) THEN
                  KTTR(JT) = 2
     !             DO WHILE (EL(KTTR(JT),I) > ELTRT(JT))
                  DO WHILE (EL(KTTR(JT),I) > ELTRT(JT)  .and. EL(KTTR(JT)+1,I) > ELTRT(JT))    ! SW 10/3/13
                    KTTR(JT) = KTTR(JT)+1
                  END DO
                  KBTR(JT) = KMX-1
                  DO WHILE (EL(KBTR(JT),I) < ELTRB(JT))
                    KBTR(JT) = KBTR(JT)-1
                  END DO
                ELSE
                  KTTR(JT) = KT
                  KBTR(JT) = KB(I)
                END IF
                KTTR(JT) = MAX(KT,KTTR(JT))
                KBTR(JT) = MIN(KB(I),KBTR(JT))
                IF (KBTR(JT) < KTTR(JT)) KBTR(JT) = KTTR(JT)
                BHSUM = 0.0
                DO K=KTTR(JT),KBTR(JT)
                  BHSUM = BHSUM+BH2(K,I)
                END DO
                DO K=KTTR(JT),KBTR(JT)
                  QTRF(K,JT) = BH2(K,I)/BHSUM
                END DO
              END IF
              DO K=KTTR(JT),KBTR(JT)
                QSS(K,I) = QSS(K,I)+QTR(JT)*QTRF(K,JT)
              END DO
            END IF
          END DO
        END IF
        IF (DIST_TRIBS(JB)) THEN
          AKBR = 0.0
          DO I=IU,ID
            AKBR = AKBR+BI(KT,I)*DLX(I)
          END DO
          DO I=IU,ID
            QDT(I)    = QDTR(JB)*BI(KT,I)*DLX(I)/AKBR
            QSS(KT,I) = QSS(KT,I)+QDT(I)
          END DO
        END IF
        IF (WITHDRAWALS) THEN
          DO JWD=1,NWD
            IF (JB == JBWD(JWD)) THEN
              I = MAX(CUS(JBWD(JWD)),IWD(JWD))
              CALL LATERAL_WITHDRAWAL    !(JWD)
              DO K=KTW(JWD),KBW(JWD)
                QSS(K,I) = QSS(K,I)-QSW(K,JWD)
              END DO
            END IF
          END DO
        END IF
        IF (UH_INTERNAL(JB)) THEN
          IF (UHS(JB) /= DS(JBUH(JB)) .OR. DHS(JBUH(JB)) /= US(JB)) THEN
            IF (JBUH(JB) >= BS(JW) .AND. JBUH(JB) <= BE(JW)) THEN
              DO K=KT,KB(IU-1)
                QSS(K,UHS(JB)) = QSS(K,UHS(JB))-VOLUH2(K,JB)/DLT
              END DO
            ELSE
              CALL UPSTREAM_FLOW
            END IF
          END IF
        END IF
        IF (DH_INTERNAL(JB)) THEN
          IF (DHS(JB) /= US(JBDH(JB)) .OR. UHS(JBDH(JB)) /= DS(JB)) THEN
            IF (JBDH(JB) >= BS(JW) .AND. JBDH(JB) <= BE(JW)) THEN
              DO K=KT,KB(ID+1)
                QSS(K,CDHS(JB)) = QSS(K,CDHS(JB))+VOLDH2(K,JB)/DLT
              END DO
            ELSE
              CALL DOWNSTREAM_FLOW
            END IF
          END IF
        END IF
      END DO
    END DO

! including tributary flows for inactive branches
    DO JW=1,NWB           ! cb 11/20/19
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        IF(BR_INACTIVE(JB))then
          IU = CUS(JB)
          ID = DS(JB)
          IF (TRIBUTARIES) THEN
           DO JT=1,JTT

!********** Inflow fractions

            IF (JB == JBTR(JT)) THEN
              I = cus(jbdn(jw))        ! placing tributary flows in upstream end of main branch
              QTRF(KT:KB(I),JT) = 0.0                
              KTTR(JT) = KT
              KBTR(JT) = KB(I)                
              KTTR(JT) = MAX(KT,KTTR(JT))
              KBTR(JT) = MIN(KB(I),KBTR(JT))
              IF (KBTR(JT) < KTTR(JT)) KBTR(JT) = KTTR(JT)
              BHSUM = 0.0
              DO K=KTTR(JT),KBTR(JT)
                BHSUM = BHSUM+BH2(K,I)
              END DO
              DO K=KTTR(JT),KBTR(JT)
                 QTRF(K,JT) = BH2(K,I)/BHSUM
              END DO
              DO K=KTTR(JT),KBTR(JT)
                QSS(K,I) = QSS(K,I)+QTR(JT)*QTRF(K,JT)
              END DO
            END IF
          END DO
         END IF
        end if
      end do
    end do

!** Compute tributary contribution to cross-shear

    IF (TRIBUTARIES) THEN
      DO JW=1,NWB
        DO JB=BS(JW),BE(JW)
          DO JT=1,JTT
            IF (JB == JBTR(JT)) THEN
              I = MAX(CUS(JB),ITR(JT))
              DO K=KTWB(JW),KBMIN(I)
                UYBR(K,I) = UYBR(K,I)+ABS(QTR(JT))*QTRF(K,JT)
              END DO
            END IF
          END DO
        END DO
      END DO
    END IF

    return
    end subroutine hydroinout
