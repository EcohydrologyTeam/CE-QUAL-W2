SUBROUTINE BALANCES

USE MAIN
USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  USE MACROPHYTEC; USE POROSITYC; USE ZOOPLANKTONC
  Use CEMAVars
  IMPLICIT NONE
  EXTERNAL RESTART_OUTPUT
  REAL VOLINJW,VOLPRJW,VOLOUTJW,VOLWDJW,VOLEVJW,VOLDTJW,VOLTRBJW,VOLICEJW,TPWB,TPSED,TNWB,TNSED,TPPLANT,TNPLANT

!***********************************************************************************************************************************
!*                                                    TASK 2.6: BALANCES                                                          **
!***********************************************************************************************************************************

    !QINT  = 0.0
    !QOUTT = 0.0
    VOLSR = 0.0
    VOLTR = 0.0
    
    DO JW=1,NWB
    VOLINJW=0.0
    VOLPRJW=0.0
    VOLOUTJW=0.0
    VOLWDJW=0.0
    VOLEVJW=0.0
    VOLDTJW=0.0
    VOLTRBJW=0.0
    VOLICEJW=0.0
    
    TPWB=0.0
    TPSED=0.0
    TNWB=0.0
    TNSED=0.0
    TNPLANT=0.0
    TPPLANT=0.0
    
      KT = KTWB(JW)
        IF (VOLUME_BALANCE(JW)) THEN
         DO JB=BS(JW),BE(JW)
          IF(.NOT.BR_INACTIVE(JB))THEN    ! SW 8/8/2018 
          VOLSBR(JB) = VOLSBR(JB)+DLVOL(JB)
          VOLTBR(JB) = VOLEV(JB)+VOLPR(JB)+VOLTRB(JB)+VOLDT(JB)+VOLWD(JB)+VOLUH(JB)+VOLDH(JB)+VOLIN(JB)+VOLOUT(JB)+VOLICE(JB)
          if(sediment_diagenesis)then
            If(CEMARelatedCode .and. IncludeBedConsolidation)Then
              VOLTBR(JB) = VOLTBR(JB)+ VOLCEMA(JB)
            End If
          ENDIF
          VOLSR(JW)  = VOLSR(JW)+VOLSBR(JB)
          VOLTR(JW)  = VOLTR(JW)+VOLTBR(JB)
          ENDIF
          !QINT(JW)   = QINT(JW) +VOLIN(JB)+VOLTRB(JB)+VOLDT(JB)+VOLPR(JB)
          !QOUTT(JW)  = QOUTT(JW)-VOLEV(JB)-VOLWD(JB) -VOLOUT(JB)
          IF (ABS(VOLSBR(JB)-VOLTBR(JB)) > VTOL .AND. VOLTBR(JB) > 100.0*VTOL) THEN
            IF (VOLUME_WARNING) THEN
              WRITE (WRN,'(A,F0.4,/A,I0,A,I0,A,I0,3(:/A,E15.8,A))')                                                                &
                            'COMPUTATIONAL WARNING AT JULIAN DAY = ',JDAY,                                                         &
                            'WATERBODY=', JW, ', BRANCH=', JB, ', KT=', KT,                                                        &
                            'SPATIAL CHANGE  =', VOLSBR(JB), ' M^3',                                                               &
                            'TEMPORAL CHANGE =', VOLTBR(JB), ' M^3',                                                               &
                            'VOLUME ERROR    =', VOLSBR(JB)-VOLTBR(JB), ' M^3'                                          !SR 11/16/19
              WRITE(WRN,*)'LAYER CHANGE:',LAYERCHANGE(JW)
              WRITE (WRN,*) 'SZ', SZ(CUS(JB):DS(JB)), 'Z', Z(CUS(JB):DS(JB)), 'H2KT', H2(KT,CUS(JB):DS(JB)),                       &
                            'H1KT', H1(KT,CUS(JB):DS(JB)), 'WSE', ELWS(CUS(JB):DS(JB)), 'Q', Q(CUS(JB):DS(JB)),                    &
                            'QC', QC(CUS(JB):DS(JB)), 'T1', T1(KT,CUS(JB):DS(JB)), 'T2', T2(KT,CUS(JB):DS(JB)),                    &
                            'SUKT', SU(KT,CUS(JB):DS(JB)), 'UKT', U(KT,CUS(JB):DS(JB)), 'QIN', QINSUM(JB),                         &
                            'QTR', QTR, 'QWD', QWD                                                                      !SR 11/16/19
              WARNING_OPEN   = .TRUE.
              VOLUME_WARNING = .FALSE.
            END IF
          END IF
          IF (VOLSR(JW) /= 0.0) DLVR(JW) = (VOLTR(JW)-VOLSR(JW))/VOLSR(JW)*100.0
            VOLINJW=VOLINJW+VOLIN(JB)
            VOLPRJW=VOLPRJW+VOLPR(JB)
            VOLOUTJW=VOLOUTJW+VOLOUT(JB)
            VOLWDJW=VOLWDJW+VOLWD(JB)
            VOLEVJW=VOLEVJW+VOLEV(JB)
            VOLDTJW=VOLDTJW+VOLDT(JB)
            VOLTRBJW=VOLTRBJW+VOLTRB(JB)
            VOLICEJW=VOLICEJW+VOLICE(JB)
         END DO

        IF (FLOWBALC=='      ON') THEN
        IF(JDAY.GE.NXFLOWBAL)THEN
        !NXFLOWBAL = NXFLOWBAL+FLOWBALF  

            IF(VOLUME_BALANCE(JW))THEN
            WRITE(FLOWBFN,'(F10.3,",",1X,I3,",",11(E16.8,",",1X))')JDAY,JW,VOLINJW,VOLPRJW,VOLOUTJW,VOLWDJW,VOLEVJW,VOLDTJW,VOLTRBJW,VOLICEJW,DLVR(JW)
            ELSE
            WRITE(FLOWBFN,'(F10.3,",",1X,I3,",",10(E16.8,",",1X))')JDAY,JW,VOLINJW,VOLPRJW,VOLOUTJW,VOLWDJW,VOLEVJW,VOLDTJW,VOLTRBJW,VOLICEJW
            ENDIF
        ENDIF
        
        END IF  ! CONTOUR INTERVAL FOR WRITING OUT FLOW BALANCE 
        ENDIF   ! VOLUME BALANCE
      IF (ENERGY_BALANCE(JW)) THEN
        ESR(JW) = 0.0
        ETR(JW) = 0.0
        DO JB=BS(JW),BE(JW)
        IF(BR_INACTIVE(JB))CYCLE    ! SW 8/8/2018 
          ETBR(JB) = EBRI(JB)+TSSEV(JB)+TSSPR(JB)+TSSTR(JB)+TSSDT(JB)+TSSWD(JB)+TSSUH(JB)+TSSDH(JB)+TSSIN(JB)+TSSOUT(JB)+TSSS(JB)  &
                     +TSSB(JB)+TSSICE(JB)
          ESBR(JB) = 0.0
          DO I=CUS(JB),DS(JB)
            DO K=KT,KB(I)
              ESBR(JB) = ESBR(JB)+T1(K,I)*DLX(I)*BH1(K,I)
            END DO
          END DO
          ETR(JW) = ETR(JW)+ETBR(JB)
          ESR(JW) = ESR(JW)+ESBR(JB)
        END DO
      END IF
      IF (MASS_BALANCE(JW)) THEN
        DO JB=BS(JW),BE(JW)
        IF(BR_INACTIVE(JB))CYCLE    ! SW 8/8/2018 
          DO JC=1,NAC
            CMBRS(CN(JC),JB) = 0.0
            DO I=CUS(JB),DS(JB)
              DO K=KT,KB(I)
                CMBRS(CN(JC),JB) = CMBRS(CN(JC),JB)+C1(K,I,CN(JC))*DLX(I)*BH1(K,I)
                CMBRT(CN(JC),JB) = CMBRT(CN(JC),JB)+(CSSB(K,I,CN(JC))+CSSK(K,I,CN(JC))*BH1(K,I)*DLX(I))*DLT
              END DO
            END DO
          END DO
          IF(DERIVED_CALC)THEN
              DO I=CUS(JB),DS(JB)
              DO K=KT,KB(I)
                TPWB  =  TPWB  + TP(K,I)    *VOL(K,I)*0.001   !/1000.   ! kg
                TPSED =  TPSED + SEDP(K,I)  *VOL(K,I)*0.001   !/1000.   ! kg        
                TNWB  =  TNWB  + TN(K,I)    *VOL(K,I)*0.001   !/1000.  ! kg
                TNSED =  TNSED + SEDN(K,I)  *VOL(K,I)*0.001   !/1000.  ! kg     
                PFLUXIN(JW) = PFLUXIN(JW) + SEDPINFLUX(K,I)*VOL(K,I)*0.001   !/1000.   ! kg
                NFLUXIN(JW) = NFLUXIN(JW) + SEDNINFLUX(K,I)*VOL(K,I)*0.001   !/1000.   ! kg
                DO M=1,NMC
                TNPLANT=TNPLANT+ MAC(K,I,M)*VOL(K,I)*MN(M)*0.001   !/1000.
                TPPLANT=TPPLANT+ MAC(K,I,M)*VOL(K,I)*MP(M)*0.001   !/1000.
                ENDDO
                DO M=1,NEP
                TNPLANT=TNPLANT+ EPM(K,I,M)*EN(M)*0.001   !/1000.
                TPPLANT=TPPLANT+ EPM(K,I,M)*EP(M)*0.001   !/1000.
                ENDDO
              END DO
            END DO
          ENDIF
          
! MACROPHYTES
          DO M=1,NMC
            IF(MACROPHYTE_CALC(JW,M))THEN
              MACMBRS(JB,M) = 0.0
              DO I=CUS(JB),DS(JB)
                IF(KTICOL(I))THEN
                  JT=KTI(I)
                ELSE
                  JT=KTI(I)+1
                END IF
                JE=KB(I)
                DO J=JT,JE
                  IF(J.LT.KT)THEN
                    COLB=EL(J+1,I)
                  ELSE
                     COLB=EL(KT+1,I)
                  END IF
                  !COLDEP=ELWS(I)-COLB
                  coldep=EL(KT,i)-Z(i)*COSA(JB)-colb  ! cb 3/7/16
                  MACMBRS(JB,M) = MACMBRS(JB,M)+MACRM(J,KT,I,M)
                  MACMBRT(JB,M) = MACMBRT(JB,M)+(MACSS(J,KT,I,M)*COLDEP*CW(J,I)*DLX(I))*DLT
                END DO
                DO K=KT+1,KB(I)
                  JT=K
                  JE=KB(I)
                  DO J=JT,JE
                    MACMBRS(JB,M) =MACMBRS(JB,M)+MACRM(J,K,I,M)
!                    MACMBRT(JB,M) = MACMBRT(JB,M)+(MACSS(J,K,I,M)*H2(K,I)*CW(J,I)*DLX(I))*DLT
                    MACMBRT(JB,M) = MACMBRT(JB,M)+(MACSS(J,K,I,M)*(CW(J,I)/B(K,I))*BH1(K,I)*DLX(I))*DLT
                  END DO
                END DO
              END DO
            END IF
          END DO
! END MACROPHYTES
        END DO
        
        IF (NPBALC=='      ON') THEN
        IF(JDAY.GE.NXNPBAL)THEN
            IF(SEDIMENT_DIAGENESIS)THEN
            WRITE(MASSBFN,'(F10.3,",",1X,I3,",",41(E16.8,",",1X))')JDAY,JW,TPWB,TPSED,TPPLANT,TPOUT(JW),TPTRIB(JW),TPDTRIB(JW),TPWD(JW),TPPR(JW),TPIN(JW),ATMDEP_P(JW),TP_SEDSOD_PO4(JW),PFLUXIN(JW),SDPFLUX(JW),TNWB,TNSED,TNPLANT,TNOUT(JW),TNTRIB(JW),TNDTRIB(JW),TNWD(JW),TNPR(JW),TNIN(JW),ATMDEP_N(JW),NH3GASLOSS(JW),TN_SEDSOD_NH4(JW),NFLUXIN(JW),SDNH4FLUX(JW),SDNO3FLUX(JW)    ! TP_SEDBURIAL(JW),TN_SEDBURIAL(JW),
                ELSE
            WRITE(MASSBFN,'(F10.3,",",1X,I3,",",31(E16.8,",",1X))')JDAY,JW,TPWB,TPSED,TPPLANT,TPOUT(JW),TPTRIB(JW),TPDTRIB(JW),TPWD(JW),TPPR(JW),TPIN(JW),ATMDEP_P(JW),TP_SEDSOD_PO4(JW),PFLUXIN(JW),TNWB,TNSED,TNPLANT,TNOUT(JW),TNTRIB(JW),TNDTRIB(JW),TNWD(JW),TNPR(JW),TNIN(JW),ATMDEP_N(JW),NH3GASLOSS(JW),TN_SEDSOD_NH4(JW),NFLUXIN(JW)    ! TP_SEDBURIAL(JW),TN_SEDBURIAL(JW),
                ENDIF
        ENDIF
        
        END IF  ! CONTOUR INTERVAL FOR WRITING OUT FLOW BALANCE 
        
        
      END IF  ! MASS BALANCE
    END DO
    
    IF(JDAY.GE.NXNPBAL)NXNPBAL = NXNPBAL+NPBALF  
    IF (FLOWBALC=='      ON') THEN  ! cb 8/22/21
        IF(JDAY.GE.NXFLOWBAL)THEN
          NXFLOWBAL = NXFLOWBAL+FLOWBALF 
        end if
    end if
    
    RETURN
    END SUBROUTINE BALANCES
