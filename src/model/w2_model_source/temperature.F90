subroutine temperature

USE MAIN
USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  USE MACROPHYTEC; USE POROSITYC; USE ZOOPLANKTONC
  Use CEMAVars
  IMPLICIT NONE
  EXTERNAL RESTART_OUTPUT
  
  REAL(R8) :: BTA1(1000),GMA1(1000)   ! places a limit of 1000 vertical layers
  REAL     :: RN1    

DO JW=1,NWB
      IF (READ_EXTINCTION(JW))GAMMA(:,US(BS(JW)):DS(BE(JW))) = EXH2O(JW)      ! SW 1/28/13
      KT = KTWB(JW)
      IF (.NOT. NO_HEAT(JW)) THEN
        IF (.NOT. READ_RADIATION(JW)) CALL SHORT_WAVE_RADIATION (JDAY)
        IF (TERM_BY_TERM(JW))THEN                                      ! SW 1/25/05
           IF(TAIR(JW).GE.5.0)THEN
           !RANLW(JW) = 5.31D-13*(273.15D0+TAIR(JW))**6*(1.0D0+0.0017D0*CLOUD(JW)**2)*0.97D0
           RANLW(JW) = 5.31D-13*(273.15D0+TAIR(JW))**6*(1.0D0+0.0017D0*CLOUD(JW)*CLOUD(JW))*0.97D0    ! SW 4/20/16 SPEED
           ELSE
           !RANLW(JW) = 5.62D-8*(273.15D0+TAIR(JW))**4*(1.D0-0.261D0*DEXP(-7.77D-4*TAIR(JW)**2))*(1.0D0+0.0017D0*CLOUD(JW)**2)*0.97D0
           RANLW(JW) = 5.62D-8*(273.15D0+TAIR(JW))**4*(1.D0-0.261D0*DEXP(-7.77D-4*TAIR(JW)*TAIR(JW)))*(1.0D0+0.0017D0*CLOUD(JW)*CLOUD(JW))*0.97D0     ! SW 4/20/16 SPEED
           ENDIF
        ENDIF
      END IF
      DO JB=BS(JW),BE(JW)
      IF(BR_INACTIVE(JB))CYCLE
        IU = CUS(JB)
        ID = DS(JB)

!****** Heat exchange

        IF (.NOT. NO_HEAT(JW)) THEN
          DO I=IU,ID
            IF (DYNAMIC_SHADE(I)) CALL SHADING

!********** Surface

            IF (.NOT. ICE(I)) THEN
              IF (TERM_BY_TERM(JW)) THEN
                CALL SURFACE_TERMS (T2(KT,I))
                RS(I)     = SRON(JW)*SHADE(I)
                RN(I)     = RS(I)+RANLW(JW)-RB(I)-RE(I)-RC(I)
                HEATEX    = RN(I)/RHOWCP*BI(KT,I)*DLX(I)
              ELSE
                CALL EQUILIBRIUM_TEMPERATURE
                HEATEX = (ET(I)-T2(KT,I))*CSHE(I)*BI(KT,I)*DLX(I)
              END IF
              TSS(KT,I) =  TSS(KT,I)+HEATEX
              TSSS(JB)  =  TSSS(JB) +HEATEX*DLT
              SROOUT    = (1.0D0-BETA(JW))*(SRON(JW)*SHADE(I)/RHOWCP)*BI(KT,I)*DLX(I)*DEXP(-GAMMA(KT,I)*DEPTHB(KT,I))
              TSS(KT,I) =  TSS(KT,I)-SROOUT
              TSSS(JB)  =  TSSS(JB) -SROOUT*DLT
              IF(KT == KB(I))THEN    ! SW 4/18/07
              SROSED    =  SROOUT*TSEDF(JW)
              ELSE
              SROSED    =  SROOUT*(1.0D0-BI(KT+1,I)/BI(KT,I))*TSEDF(JW)
              ENDIF
              TSS(KT,I) =  TSS(KT,I)+SROSED
              TSSS(JB)  =  TSSS(JB) +SROSED*DLT
              SROIN     =  SROOUT*B(KT+1,I)/BI(KT,I)
              DO K=KT+1,KB(I)
                SROOUT   = SROIN*DEXP(-GAMMA(K,I)*(H1(K,I)))
                SRONET   = SROIN-SROOUT
                IF(K /= KB(I))THEN                                         ! SW 1/18/08
                SROSED   = SROOUT*(1.0D0-BI(K+1,I)/BI(K,I))*TSEDF(JW)
                ELSE
                SROSED   = SROOUT*TSEDF(JW)
                ENDIF
                TSS(K,I) = TSS(K,I)+ SRONET+SROSED
                TSSS(JB) = TSSS(JB)+(SRONET+SROSED)*DLT
                SROIN    = SROOUT*B(K+1,I)/B(K,I)
              END DO
            END IF

!********** Sediment/water

            DO K=KT,KB(I)
              IF(K==KB(I))THEN                ! SW 4/18/07
              TFLUX    = CBHE(JW)/RHOWCP*(TSED(JW)-T2(K,I))*BI(K,I)*DLX(I)
              ELSE
              TFLUX    = CBHE(JW)/RHOWCP*(TSED(JW)-T2(K,I))*(BI(K,I)-BI(K+1,I))*DLX(I)
              ENDIF
              TSS(K,I) = TSS(K,I)+TFLUX
              TSSB(JB) = TSSB(JB)+TFLUX*DLT
            END DO
          END DO

!******** Ice cover

          IF (ICE_CALC(JW)) THEN
 !           HIA = 0.2367*CSHE(I)/5.65E-8    ! SW 10/20/09 Duplicate line of code
            DO I=IU,ID
              ALLOW_ICE(I) = .TRUE.
              IF (T2(KT,I) > ICET2(JW)) ALLOW_ICE(I) = .FALSE.                        ! RC/SW 4/28/11
 !             DO K=KT,KB(I)                                                          ! RC 4/28/11: no reason to loop over all layers
 !               IF (T2(K,I) > ICET2(JW)) ALLOW_ICE(I) = .FALSE.
 !             END DO
            END DO
  !          ICE_IN(JB) = .TRUE.                                                      ! RC/SW 4/28/11 eliminate ICE_IN
  !          DO I=IU,ID
  !            IF (ICETH(I) < ICEMIN(JW)) ICE_IN(JB) = .FALSE.
  !          END DO
            DO I=IU,ID
              IF(SALT_WATER(JW))THEN        ! SW/RC 4/28/11
                IF(TDS(KT,I) < 35.)THEN
                RIMT = -0.0545*TDS(KT,I)                                         ! REGRESSION FOR TDS BETWEEN 0 AND 35 PPT
                ELSE
                RIMT=-0.31462-0.04177*TDS(KT,I)-0.000166*TDS(KT,I)*TDS(KT,I)     ! REGRESSION EQN FOR TDS>35 PPT
                ENDIF                                                                                
              ELSE
              RIMT=0.0
              ENDIF
              IF (DETAILED_ICE(JW)) THEN
                IF (T2(KT,I) < 0.0) THEN
                  IF (.NOT. ICE(I)) THEN
                    ICETH2 = -T2(KT,I)*RHO(KT,I)*CP*H2(KT,I)/RHOIRL1
                    IF (ICETH2 < ICE_TOL) THEN
                      ICETH2 = 0.0D0
                    ELSE
                      TFLUX      = T2(KT,I)*RHO(KT,I)*CP*H2(KT,I)*BI(KT,I)/(RHOWCP*DLT)*DLX(I)
                      TSS(KT,I)  = TSS(KT,I) -TFLUX
                      TSSICE(JB) = TSSICE(JB)-TFLUX*DLT
                    END IF
                  END IF
                END IF

!************** Ice balance

                IF (ICE(I)) THEN
                  TICE = TAIR(JW)
                  DEL  = 2.0D0
                  J    = 1
                    IF(TAIR(JW).GE.5.0)THEN
                    RANLW(JW) = 5.31D-13*(273.15D0+TAIR(JW))**6*(1.0D0+0.0017D0*CLOUD(JW)**2)*0.97D0
                    ELSE
                    RANLW(JW) = 5.62D-8*(273.15D0+TAIR(JW))**4*(1.D0-0.261D0*DEXP(-7.77D-4*TAIR(JW)**2))*(1.0D0+0.0017D0*CLOUD(JW)**2)*0.97D0
                    ENDIF
                    RN1=SRON(JW)/REFL*SHADE(I)*(1.0D0-ALBEDO(JW))*BETAI(JW)+RANLW(JW)               ! SW 4/19/10 eliminate spurious divsion of SRO by RHOCP
                  DO WHILE (ABS(DEL) > 1.0 .AND. J < 500)                                         ! SW 4/21/10 Should have been ABS of DEL
                    CALL SURFACE_TERMS (TICE)      
                    RN(I) = RN1-RB(I)-RE(I)-RC(I)    ! 4/19/10 
!                    RN(I) = SRON(JW)/(REFL*RHOWCP)*SHADE(I)*(1.0-ALBEDO(JW))*BETAI(JW)+RANLW(JW)-RB(I)-RE(JW)-RC(I)
                    DEL   = RN(I)+RK1*(RIMT-TICE)/ICETH(I)     ! RK1 is ice conductivity 2.12 W/m/oC
                    IF (ABS(DEL) > 1.0) TICE = TICE+DEL/500.0D0
                    J = J+1
                  END DO

!**************** Solar radiation attenuation

                  TFLUX      = DLX(I)*SRON(JW)/(RHOWCP*REFL)*SHADE(I)*(1.0D0-ALBEDO(JW))*(1.0D0-BETAI(JW))                             &   ! SW 4/21/10 Eliminate spurious divide by RHOCP
                               *DEXP(-GAMMAI(JW)*ICETH(I))*BI(KT,I)
                  TSS(KT,I)  = TSS(KT,I) +TFLUX
                  TSSICE(JB) = TSSICE(JB)+TFLUX*DLT
                  IF (TICE > 0.0) THEN
                    HICE   =  RHOICP*0.5D0*TICE*0.5D0*ICETH(I)*BI(KT,I)/(RHOWCP*DLT)
                    ICETHU = -DLT*HICE/B(KTI(I),I)*RHOWCP/RHOIRL1
                    TICE   =  0.0D0
                  END IF

!**************** Ice growth

                  IF (TICE < 0.0) ICETH1 = DLT*(RK1*(RIMT-TICE)/ICETH(I))/RHOIRL1

!**************** Ice melt from water-ice interface

                  IF (T2(KT,I) > 0.0) THEN
                    ICETH2     = -DLT*HWI(JW)*(T2(KT,I)-RIMT)/RHOIRL1
                    TFLUX      =  2.392D-7*HWI(JW)*(RIMT-T2(KT,I))*BI(KT,I)*DLX(I)
                    TSS(KT,I)  =  TSS(KT,I) +TFLUX
                    TSSICE(JB) =  TSSICE(JB)+TFLUX*DLT
                  END IF
                END IF

!************** Ice thickness

                ICETH(I) = ICETH(I)+ICETHU+ICETH1+ICETH2
                ! SW 9/4/15
                IF(ICEC(JW) == '    ONWB')THEN
                  IceThicknessChange = ICETHU+ICETH1+ICETH2
                  IF(IceThicknessChange >= 0.0  .and. ICETH(I) > ICE_TOL)THEN
                  IceQSS(I) =-IceThicknessChange*BI(KT,I)*DLX(I)*0.917/DLT    ! CEMA: removal of 0.917 m3 of water for every 1 m3 of ice formed
                  IceBank(I)=IceBank(I)+IceThicknessChange*BI(KT,I)*DLX(I)    ! Icebank - always + - volume of ice in m3
                  ELSE
                    IF (ICETH(I) < ICE_TOL)THEN
                        ICETH(I) = 0.0D0
                        ICEQSS(I)=ICEBANK(I)*0.917/DLT
                        ICEBANK(I)=0.0
                        IF(I==IU .AND. US(JB) < IU)THEN    ! CHECK SUBTRACTED SEGMENTS THAT MAY HAVE ICE, MELT THEM AND PUT WATER IN SEGMENT IU
                            DO II=US(JB),IU-1
                            ICETH(II)=0.0D0
                            ICEQSS(IU)=ICEQSS(IU)+ICEBANK(II)*0.917/DLT
                            ICEBANK(II)=0.0
                            ICE(II) = .FALSE.
                            ENDDO
                        ENDIF
                        
                    ELSE
                        ICEQSS(I)=-(ICETHICKNESSCHANGE/ICETH(I))*ICEBANK(I)*0.917/DLT     ! SW 9/29/15
                        ICEBANK(I)=(1.+ICETHICKNESSCHANGE/ICETH(I))*ICEBANK(I)            ! Note: ICETHICKNESSCHANGE is negative
                        IF(I==IU .AND. US(JB) < IU)THEN    ! CHECK SUBTRACTED SEGMENTS THAT MAY HAVE ICE, MELT THEM AND PUT WATER IN SEGMENT IU
                            DO II=US(JB),IU-1
                                ICEQSS(IU)=ICEQSS(IU)-(ICETHICKNESSCHANGE/ICETH(IU))*ICEBANK(II)*0.917/DLT 
                                ICEBANK(II)=(1.+ICETHICKNESSCHANGE/ICETH(IU))*ICEBANK(II)
                            ENDDO
                        ENDIF
                    ENDIF
                  ENDIF
                  !VolIce(jb)=VolIce(jb)+iceqss(i)*dlt   ! since this flow is not exercised until the next time step - moved code to main program calculation of qss
                 ENDIF 
                  
                IF (ICETH(I) < ICE_TOL) ICETH(I) = 0.0D0
     !           IF (WINTER .AND. (.NOT. ICE_IN(JB))) THEN            ! RC 4/28/11 No reason for this
     !             IF (.NOT. ALLOW_ICE(I)) ICETH(I) = 0.0
     !           END IF
                ICE(I)   = ICETH(I) > 0.0
                IF (ICE(I))THEN    ! 3/27/08 SW
                 ICESW(I) = 0.0
                 ELSE
                 ICESW(I) = 1.0
                ENDIF
                ICETHU   = 0.0
                ICETH1   = 0.0
                ICETH2   = 0.0
                IF (ICETH(I) < ICE_TOL .AND. ICETH(I) > 0.0) ICETH(I) = ICE_TOL
              ELSE                                                         ! IF no ice the preceding time step
                IF(TERM_BY_TERM(JW))CALL EQUILIBRIUM_TEMPERATURE           ! SW 10/20/09 Must call this first otherwise ET and CSHE are 0
                HIA      = 0.2367D0*CSHE(I)/5.65D-8                          ! JM 11/08 convert SI units of m/s to English (btu/ft2/d/F) and then back to SI W/m2/C
!                ICETH(I) = MAX(0.0,ICETH(I)+DLT*((RIMT-ET(I))/(ICETH(I)/RK1+1.0/HIA)-(T2(KT,I)-RIMT))/RHOIRL1)
                ICETH(I) = MAX(0.0,ICETH(I)+DLT*((RIMT-ET(I))/(ICETH(I)/RK1+1.0D0/HIA)-HWI(JW)*(T2(KT,I)-RIMT))/RHOIRL1)   ! SW 10/20/09 Revised missing HWI(JW)
                ICE(I)   = ICETH(I) > 0.0
                ICESW(I) = 1.0
                IF (ICE(I)) THEN
!                  TFLUX      = 2.392E-7*(RIMT-T2(KT,I))*BI(KT,I)*DLX(I)
                  TFLUX      = 2.392D-7*HWI(JW)*(RIMT-T2(KT,I))*BI(KT,I)*DLX(I)      ! SW 10/20/09 Revised missing HWI(JW)
                  TSS(KT,I)  = TSS(KT,I) +TFLUX
                  TSSICE(JB) = TSSICE(JB)+TFLUX*DLT
                  ICESW(I) = 0.0
                END IF
              END IF
            END DO
          END IF
        END IF
          !DO I=IU,ID
          !  icebank_all=icebank(i)+icebank_all
          !END DO
          !write(25600,*)jday,icebank_all
!****** Heat sources/sinks and total inflow/outflow

        IF (EVAPORATION(JW)) THEN
          DO I=IU,ID
            TSS(KT,I) = TSS(KT,I)-EV(I)*T2(KT,I)
            TSSEV(JB) = TSSEV(JB)-EV(I)*T2(KT,I)*DLT
            VOLEV(JB) = VOLEV(JB)-EV(I)         *DLT
          END DO
        END IF
        IF (PRECIPITATION(JW)) THEN
          DO I=IU,ID
            TSS(KT,I) = TSS(KT,I)+QPR(I)*TPR(JB)
            TSSPR(JB) = TSSPR(JB)+QPR(I)*TPR(JB)*DLT
            VOLPR(JB) = VOLPR(JB)+QPR(I)        *DLT
          END DO
        END IF
        IF (TRIBUTARIES) THEN
          DO JT=1,JTT
            IF (JB == JBTR(JT)) THEN
              I = ITR(JT)
              IF (I < CUS(JB)) I = CUS(JB)
              DO K=KTTR(JT),KBTR(JT)
                IF (QTR(JT) < 0) THEN
                  TSS(K,I)  = TSS(K,I) +T2(K,I)*QTR(JT)*QTRF(K,JT)
                  TSSTR(JB) = TSSTR(JB)+T2(K,I)*QTR(JT)*QTRF(K,JT)*DLT
                ELSE
                  TSS(K,I)  = TSS(K,I) +TTR(JT)*QTR(JT)*QTRF(K,JT)
                  TSSTR(JB) = TSSTR(JB)+TTR(JT)*QTR(JT)*QTRF(K,JT)*DLT
                END IF
              END DO
              VOLTRB(JB) = VOLTRB(JB)+QTR(JT)*DLT
            END IF
          END DO
        END IF
        IF (DIST_TRIBS(JB)) THEN
          DO I=IU,ID
            IF (QDT(I) < 0) THEN
              TSS(KT,I) = TSS(KT,I)+T2(KT,I)*QDT(I)
              TSSDT(JB) = TSSDT(JB)+T2(KT,I)*QDT(I)*DLT
            ELSE
              TSS(KT,I) = TSS(KT,I)+TDTR(JB)*QDT(I)
              TSSDT(JB) = TSSDT(JB)+TDTR(JB)*QDT(I)*DLT
            END IF
            VOLDT(JB) = VOLDT(JB)+QDT(I)*DLT
          END DO
        END IF
        IF (WITHDRAWALS) THEN
          DO JWD=1,JWW
            IF (QWD(JWD) /= 0.0) THEN
              IF (JB == JBWD(JWD)) THEN
                I = MAX(CUS(JBWD(JWD)),IWD(JWD))
                DO K=KTW(JWD),KBW(JWD)
                  TSS(K,I)  = TSS(K,I) -T2(K,I)*QSW(K,JWD)
                  TSSWD(JB) = TSSWD(JB)-T2(K,I)*QSW(K,JWD)*DLT
                END DO
                VOLWD(JB) = VOLWD(JB)-QWD(JWD)*DLT
              END IF
            END IF
          END DO
        END IF
        IF (UP_FLOW(JB)) THEN
            DO K=KT,KB(IU)
              IF (.NOT. HEAD_FLOW(JB)) THEN
                TSS(K,IU) = TSS(K,IU)+QINF(K,JB)*QIN(JB)*TIN(JB)
                TSSIN(JB) = TSSIN(JB)+QINF(K,JB)*QIN(JB)*TIN(JB)*DLT
              ELSE
                IF (U(K,IU-1) >= 0.0) THEN
                  TSS(K,IU) = TSS(K,IU)+U(K,IU-1)*BHR1(K,IU-1)*T1(K,IU-1)
                  TSSIN(JB) = TSSIN(JB)+U(K,IU-1)*BHR1(K,IU-1)*T1(K,IU-1)*DLT
                ELSE
                  TSS(K,IU) = TSS(K,IU)+U(K,IU-1)*BHR1(K,IU-1)*T1(K,IU)
                  TSSIN(JB) = TSSIN(JB)+U(K,IU-1)*BHR1(K,IU-1)*T1(K,IU)*DLT
                END IF
              END IF
            END DO
          VOLIN(JB) = VOLIN(JB)+QIN(JB)*DLT
        END IF
        IF (DN_FLOW(JB)) THEN
          DO K=KT,KB(ID)
            TSS(K,ID)  = TSS(K,ID) -QOUT(K,JB)*T2(K,ID+1)
            TSSOUT(JB) = TSSOUT(JB)-QOUT(K,JB)*T2(K,ID+1)*DLT
            VOLOUT(JB) = VOLOUT(JB)-QOUT(K,JB)           *DLT
          END DO
        END IF
        IF (UP_HEAD(JB)) THEN
          DO K=KT,KB(IU)
            IUT = IU
            IF (QUH1(K,JB) >= 0.0) IUT = IU-1
            TSSUH1(K,JB) = T2(K,IUT)*QUH1(K,JB)
            TSS(K,IU)    = TSS(K,IU)+TSSUH1(K,JB)
            TSSUH(JB)    = TSSUH(JB)+TSSUH1(K,JB)*DLT
            VOLUH(JB)    = VOLUH(JB)+QUH1(K,JB)  *DLT
          END DO
        END IF
        IF (UH_INTERNAL(JB)) THEN
          IF (UHS(JB) /= DS(JBUH(JB)) .OR. DHS(JBUH(JB)) /= US(JB)) THEN
            IF (JBUH(JB) >= BS(JW) .AND. JBUH(JB) <= BE(JW)) THEN
              DO K=KT,KB(IU-1)
                TSS(K,UHS(JB))  = TSS(K,UHS(JB)) -TSSUH2(K,JB)/DLT
                TSSUH(JBUH(JB)) = TSSUH(JBUH(JB))-TSSUH2(K,JB)
                VOLUH(JBUH(JB)) = VOLUH(JBUH(JB))-VOLUH2(K,JB)
              END DO
            ELSE
              CALL UPSTREAM_CONSTITUENT(T2,TSS)
              DO K=KT,KB(IU-1)
                TSSUH(JBUH(JB)) = TSSUH(JBUH(JB))-TSSUH2(K,JB)
                VOLUH(JBUH(JB)) = VOLUH(JBUH(JB))-VOLUH2(K,JB)
              END DO
            END IF
          END IF
        END IF
        IF (DN_HEAD(JB)) THEN
          DO K=KT,KB(ID+1)
            IDT = ID+1
            IF (QDH1(K,JB) >= 0.0) IDT = ID
            TSSDH1(K,JB) = T2(K,IDT)*QDH1(K,JB)
            TSS(K,ID)    = TSS(K,ID)-TSSDH1(K,JB)
            TSSDH(JB)    = TSSDH(JB)-TSSDH1(K,JB)*DLT
            VOLDH(JB)    = VOLDH(JB)-QDH1(K,JB)  *DLT
          END DO
        END IF
        IF (DH_INTERNAL(JB)) THEN
          IF (DHS(JB) /= US(JBDH(JB)) .OR. UHS(JBDH(JB)) /= DS(JB)) THEN
            IF (JBDH(JB) >= BS(JW) .AND. JBDH(JB) <= BE(JW)) THEN
              DO K=KT,KB(ID+1)
                TSS(K,CDHS(JB)) = TSS(K,CDHS(JB))+TSSDH2(K,JB)/DLT
                TSSDH(JBDH(JB)) = TSSDH(JBDH(JB))+TSSDH2(K,JB)
                VOLDH(JBDH(JB)) = VOLDH(JBDH(JB))+VOLDH2(K,JB)
              END DO
            ELSE
              CALL DOWNSTREAM_CONSTITUENT(T2,TSS)
              DO K=KT,KB(ID+1)
                TSSDH(JBDH(JB)) = TSSDH(JBDH(JB))+TSSDH2(K,JB)
                VOLDH(JBDH(JB)) = VOLDH(JBDH(JB))+VOLDH2(K,JB)
              END DO
            END IF
          END IF
        END IF
      END DO
    END DO

IF(NIT==0)THEN
    IF (DERIVED_CALC)THEN
        DO JW=1,NWB
        KT=KTWB(JW)
        DO JB=BS(JW),BE(JW)
            IU=CUS(JB)
            ID=DS(JB)
        CALL TEMPERATURE_RATES
        CALL KINETIC_RATES
        IF(CDWBC(PH_DER,JW)=='      ON')CALL PH_CO2
        ENDDO        
        CALL DERIVED_CONSTITUENTS
        ENDDO
    ENDIF
    CALL OUTPUTA    ! INITIALIZE OUTPUT FOR FIRST TIME STEP IF OUTPUT IS SET FOR FIRST TIME STEP ONLY ISSUE IS VELOCITY AND FLOW IS AT T+DT RATHER THAN AT T=0
ENDIF

!** Temperature transport
    COLD => HYD(:,:,4)
    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
      IF(BR_INACTIVE(JB))CYCLE
        IU =  CUS(JB)
        ID =  DS(JB)
        !COLD => HYD(:,:,4)
        CALL HORIZONTAL_MULTIPLIERS1
        CALL VERTICAL_MULTIPLIERS1
        CALL HORIZONTAL_MULTIPLIERS
        CALL VERTICAL_MULTIPLIERS
     !   CNEW => T1(:,:)
     !   SSB  => TSS(:,:)
     !   SSK  => CSSB(:,:,1)
     !   CALL HORIZONTAL_TRANSPORT
        
        DO I=IU,ID    !CONCURRENT(I=IU:ID)   !
        DO K=KT,KB(I)     !CONCURRENT(K=KT:KB(I))      !FORALL                                         !DO K=KT,KB(I)
        DT(K,I) = (COLD(K,I)*BH2(K,I)/DLT+(ADX(K,I)*BHR1(K,I)-ADX(K,I-1)*BHR1(K,I-1))/DLX(I)+(1.0D0-THETA(JW))                     &
                    *(ADZ(K,I)*BB(K,I)-ADZ(K-1,I)*BB(K-1,I))+TSS(K,I)/DLX(I))*DLT/BH1(K,I)
       END DO
       END DO
        
        
        DO I=IU,ID   !CONCURRENT(I=IU:ID)   !
          DO K=KT,KB(I)
            AT(K,I) = 0.0D0; CT(K,I) = 0.0D0; VT(K,I) = 0.0D0    !; DT(:,I) = 0.0D0    SW CODE SPEEDUP 6/15/13
          ENDDO
          DO K=KT,KB(I)   !CONCURRENT(K=KT:KB(I))          !FORALL(K=KT:KB(I))                                                 !DO K=KT,KB(I)
            AT(K,I) = -DLT/BH1(K,I)*(BB(K-1,I)*(DZ(K-1,I)/AVH1(K-1,I)+THETA(JW)*0.5D0*W(K-1,I)))
            CT(K,I) =  DLT/BH1(K,I)*(BB(K,I)*(THETA(JW)*0.5D0*W(K,I)-DZ(K,I)/AVH1(K,I)))
            VT(K,I) =  1.0D0+DLT/BH1(K,I)*(BB(K,I)*(DZ(K,I)/AVH1(K,I)+THETA(JW)*0.5D0*W(K,I))+BB(K-1,I)*(DZ(K-1,I)/AVH1(K-1,I)         &
                      -THETA(JW)*0.5D0*W(K-1,I)))
       !     DT(K,I) =  CNEW(K,I)
          END DO
         ! CALL TRIDIAG(AT(:,I),VT(:,I),CT(:,I),DT(:,I),KT,KB(I),KMX,CNEW(:,I))
           BTA1(KT)=VT(KT,I)
                    GMA1(KT)=DT(KT,I)
                  DO K=KT+1,KB(I)
                    BTA1(K) = VT(K,I)-AT(K,I)/BTA1(K-1)*CT(K-1,I)
                    GMA1(K) = DT(K,I)-AT(K,I)/BTA1(K-1)*GMA1(K-1)
                  END DO
                  T1(KB(I),I) = GMA1(KB(I))/BTA1(KB(I))
                  DO K=KB(I)-1,KT,-1
                    T1(K,I) = (GMA1(K)-CT(K,I)*T1(K+1,I))/BTA1(K)
                  END DO
        END DO
      END DO
    END DO


    END SUBROUTINE TEMPERATURE
