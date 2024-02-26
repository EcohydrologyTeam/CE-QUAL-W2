SUBROUTINE WQCONSTITUENTS

USE MAIN
USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  USE MACROPHYTEC; USE POROSITYC; USE ZOOPLANKTONC;USE TRIDIAG_V
  Use CEMAVars
  Use CEMASedimentDiagenesis, only: SedimentFlux; USE ALGAE_TOXINS
  
  IMPLICIT NONE
  EXTERNAL RESTART_OUTPUT
  REAL :: TPALG,TNALG,TPZ,TNZ,TPBOD,TNBOD

      IF(MACROPHYTE_ON.AND.UPDATE_KINETICS)CALL POROSITY 
      DO JW=1,NWB
        KT = KTWB(JW)
        DO JB=BS(JW),BE(JW)
        IF(BR_INACTIVE(JB))CYCLE    ! SW 6/12/2017
          IU = CUS(JB)
          ID = DS(JB)

!******** Kinetic sources/sinks

          IF (SEDIMENT_CALC(JW))then
            CALL SEDIMENT
            CALL SEDIMENTP
            CALL SEDIMENTN
            CALL SEDIMENTC
            IF(DYNSEDK(JW) == '      ON')CALL SEDIMENT_DECAY_RATE
          END IF
! Standing biomass decay
        IF(STANDING_BIOMASS_DECAY)THEN  ! SW 5/26/15
          IF (SEDIMENT_CALC1(JW))then
            CALL SEDIMENT1            
          end if
          IF (SEDIMENT_CALC2(JW))then
            CALL SEDIMENT2            
          end if
        ENDIF

          DO M=1,NMC
            IF (MACROPHYTE_CALC(JW,M))THEN
              CALL MACROPHYTE(M)
            END IF
          END DO

          IF (UPDATE_KINETICS) THEN
            IF (UPDATE_RATES) THEN
              CALL TEMPERATURE_RATES
              CALL KINETIC_RATES
            END IF
            DO JAC=1,NAC
              JC = CN(JAC)
              IF (JC == NPO4)                    CALL PHOSPHORUS
              IF (JC == NWAGE)                   CALL WATER_AGE
              IF (JC == NBACT)                   CALL BACTERIA
              IF (JC == NDGP)                    CALL DISSOLVED_GAS
              IF (JC == NN2)                     CALL DISSOLVED_N2
              IF (JC == NH2S)                    CALL SULFIDE
              IF (JC == NCH4)                    CALL METHANE
              IF (JC == NSO4)                    CALL SULFATE
              IF (JC == NFEII)                   CALL FERROUS
              IF (JC == NFEOOH)                  CALL OXIDIZEDFE
              IF (JC == NMNII)                   CALL BIVALENTMN
              IF (JC == NMNO2)                   CALL OXIDIZEDMN
              
              IF (JC == NNH4)                    CALL AMMONIUM
              IF (JC == NNO3)                    CALL NITRATE
              IF (JC == NDSI)                    CALL DISSOLVED_SILICA
              IF (JC == NPSI)                    CALL PARTICULATE_SILICA
              !
              IF(ORGC_CALC)THEN
                IF (JC == NLDOMC)                CALL LABILE_DOM_C
                IF (JC == NRDOMC)                CALL REFRACTORY_DOM_C
                IF (JC == NLPOMC)                CALL LABILE_POM_C
                IF (JC == NRPOMC)                CALL REFRACTORY_POM_C
              ELSE
              IF (JC == NLDOM)                   CALL LABILE_DOM
              IF (JC == NRDOM)                   CALL REFRACTORY_DOM
              IF (JC == NLPOM)                   CALL LABILE_POM
              IF (JC == NRPOM)                   CALL REFRACTORY_POM
              END IF
              IF (JC == NDO)                     CALL DISSOLVED_OXYGEN              
              IF (JC >= NGCS  .AND. JC <= NGCE)  CALL GENERIC_CONST(JC-NGCS+1)              
              IF (JC >= NSSS  .AND. JC <= NSSE)  CALL SUSPENDED_SOLIDS(JC-NSSS+1)
              IF (JC >= NAS   .AND. JC <= NAE)THEN
                IF(ALG_CALC(JC-NAS+1))CALL ALGAE(JC-NAS+1)
              ENDIF
              IF (JC >= NBODS .AND. JC <= NBODE)THEN
                DO JCB=1,NBOD       ! VARIABLE STOICHIOMETRY FOR CBOD, CB 6/6/10
                  IF(BOD_CALC(JCB))THEN
                    IF(JC == NBODC(JCB))CALL BIOCHEMICAL_O2_DEMAND(JCB)
                    IF(JC == NBODP(JCB) .AND. BOD_CALCP(JCB))CALL BIOCHEMICAL_O2_DEMAND_P(JCB)         ! CB 5/19/2011
                    IF(JC == NBODN(JCB) .AND. BOD_CALCN(JCB))CALL BIOCHEMICAL_O2_DEMAND_N(JCB)         ! CB 5/19/2011
                  END IF
                END DO       
              ENDIF
              IF (JC >= NZOOS  .AND. JC <= NZOOE .AND.ZOOPLANKTON_CALC)CALL ZOOPLANKTON  		
              IF (JC == NLDOMP)                CALL LABILE_DOM_P
              IF (JC == NRDOMP)                CALL REFRACTORY_DOM_P
              IF (JC == NLPOMP)                CALL LABILE_POM_P
              IF (JC == NRPOMP)                CALL REFRACTORY_POM_P
              IF (JC == NLDOMN)                CALL LABILE_DOM_N
              IF (JC == NRDOMN)                CALL REFRACTORY_DOM_N
              IF (JC == NLPOMN)                CALL LABILE_POM_N
              IF (JC == NRPOMN)                CALL REFRACTORY_POM_N
              !IF (JC == NALK .and. NONCON_ALKALINITY)                  CALL alkalinity
              IF (JC == NALK)                  CALL alkalinity    ! NW 2/11/16
              IF (JC >= NATS    .AND. JC <= NATE .AND. ALGAE_TOXIN)THEN
                  CALL INTRACELLULAR_TOXIN(JC-NATS+1)
                  CALL EXTRACELLULAR_TOXIN(JC-NATS+1)
              ENDIF
              
            END DO
            IF (PH_CALC(JW)) CALL INORGANIC_CARBON
            IF (PH_CALC(JW)) THEN
              if(ph_buffering)then  ! enhanced pH buffering                
                call pH_CO2_new
              else
                CALL PH_CO2
              end if
            END IF
            If(CEMARelatedCode .and. IncludeCEMASedDiagenesis) Call SedimentFlux
            
          END IF          
          DO JE=1,NEP   ! sw 5/16/06
            IF (EPIPHYTON_CALC(JW,JE)) CALL EPIPHYTON(JE)
          END DO

!******** External sources/sinks

            IF(AERATEC == "      ON")CALL AERATEMASS
            IF(EVAPORATION(JW) .AND. WATER_AGE_ACTIVE)THEN    ! CORRECT WATER AGE FOR EVAPORATION SR 7/27/2017
                DO I=IU,ID
                    !JC=NGCS+JG_AGE-1
                    CSSB(KT,I,NWAGE)=CSSB(KT,I,NWAGE)-EV(I)*WAGE(KT,I)   ! SW 10/17/2019 CG(KT,I,JC)           
                ENDDO
            ENDIF

          DO JAC=1,NAC
            JC = CN(JAC)
            IF (TRIBUTARIES) THEN
              DO JT=1,JTT
                IF (JB == JBTR(JT)) THEN
                  I = ITR(JT)
                  IF (I < CUS(JB)) I = CUS(JB)
                  DO K=KTTR(JT),KBTR(JT)
                    IF (QTR(JT) < 0.0) THEN
                      CSSB(K,I,JC) = CSSB(K,I,JC)+C1(K,I,JC)*QTR(JT)*QTRF(K,JT)
                    ELSE
                      CSSB(K,I,JC) = CSSB(K,I,JC)+CTR(JC,JT)*QTR(JT)*QTRF(K,JT)
                    END IF
                  END DO
                END IF
              END DO
            END IF
            IF (DIST_TRIBS(JB)) THEN
              DO I=IU,ID
                IF (QDT(I) < 0.0) THEN
                  CSSB(KT,I,JC) = CSSB(KT,I,JC)+C1(KT,I,JC)*QDT(I)
                ELSE
                  CSSB(KT,I,JC) = CSSB(KT,I,JC)+CDTR(JC,JB)*QDT(I)
                END IF
              END DO
            END IF
            IF (WITHDRAWALS) THEN
              DO JWD=1,JWW
                IF (QWD(JWD) /= 0.0) THEN
                  IF (JB == JBWD(JWD)) THEN
                    I = MAX(CUS(JBWD(JWD)),IWD(JWD))
                    DO K=KTW(JWD),KBW(JWD)      !CONCURRENT(K=KTW(JWD):KBW(JWD))                 ! FORALL
                      CSSB(K,I,JC) = CSSB(K,I,JC)-C1S(K,I,JC)*QSW(K,JWD)
                    END DO
                  END IF
                END IF
              END DO
            END IF
            IF (PRECIPITATION(JW)) THEN
              DO I=IU,ID    !CONCURRENT (I=IU:ID)                                  !FORALL
                CSSB(KT,I,JC) = CSSB(KT,I,JC)+CPR(JC,JB)*QPR(I)
              END DO
            END IF
            IF (UP_FLOW(JB)) THEN
              DO K=KT,KB(IU)
                IF (.NOT. HEAD_FLOW(JB)) THEN
                  CSSB(K,IU,JC) = CSSB(K,IU,JC)+QINF(K,JB)*QIN(JB)*CIN(JC,JB)
                ELSE
                  IF (U(K,IU-1) >= 0.0) THEN
                    CSSB(K,IU,JC) = CSSB(K,IU,JC)+U(K,IU-1)*BHR1(K,IU-1)*C1S(K,IU-1,JC)
                  ELSE
                    CSSB(K,IU,JC) = CSSB(K,IU,JC)+U(K,IU-1)*BHR1(K,IU-1)*C1S(K,IU,JC)
                  END IF
                END IF
              END DO
            END IF
            IF (DN_FLOW(JB)) CSSB(KT:KB(ID),ID,JC) = CSSB(KT:KB(ID),ID,JC)-QOUT(KT:KB(ID),JB)*C1S(KT:KB(ID),ID,JC)
            IF (UP_HEAD(JB)) THEN
                DO K=KT,KB(IU)
                IUT = IU
                IF (QUH1(K,JB) >= 0.0) IUT = IU-1
                CSSUH1(K,JC,JB) = C1S(K,IUT,JC)*QUH1(K,JB)
                CSSB(K,IU,JC)   = CSSB(K,IU,JC)+CSSUH1(K,JC,JB)
              END DO
              IF (UH_INTERNAL(JB)) THEN
                IF (UHS(JB) /= DS(JBUH(JB)) .OR. DHS(JBUH(JB)) /= US(JB)) THEN
                  IF (JBUH(JB) >= BS(JW) .AND. JBUH(JB) <= BE(JW)) THEN
                    I = UHS(JB)
                    !dir$ ivdep
                    DO K=KT,KB(IU)
                      CSSB(K,I,JC) = CSSB(K,I,JC)-CSSUH2(K,JC,JB)/DLT
                    END DO
                  ELSE
                    CALL UPSTREAM_CONSTITUENT(C2(:,:,JC),CSSB(:,:,JC))
                  END IF
                END IF
              END IF
            END IF
            IF (DN_HEAD(JB)) THEN
              DO K=KT,KB(ID+1)
                IDT = ID+1
                IF (QDH1(K,JB) >= 0.0) IDT = ID
                CSSDH1(K,JC,JB) = C1S(K,IDT,JC)*QDH1(K,JB)
                CSSB(K,ID,JC)   = CSSB(K,ID,JC)-CSSDH1(K,JC,JB)
              END DO
              IF (DH_INTERNAL(JB)) THEN
                IF (DHS(JB) /= US(JBDH(JB)) .OR. UHS(JBDH(JB)) /= DS(JB)) THEN
                  IF (JBDH(JB) >= BS(JW) .AND. JBDH(JB) <= BE(JW)) THEN
                    I = DHS(JB)
                    DO K=KT,KB(ID+1)
                      CSSB(K,I,JC) = CSSB(K,I,JC)+CSSDH2(K,JC,JB)/DLT
                    END DO
                  ELSE
                    CALL DOWNSTREAM_CONSTITUENT(C2(:,:,JC),CSSB(:,:,JC))
                  END IF
                END IF
              END IF
            END IF
          END DO
          
        IF (NPBALC=='      ON') THEN           !IF(MASS_BALANCE(JW).AND.CONTOUR(JW).AND.DERIVED_CALC)THEN     ! TO COMPUTE TP AND TN INFLOWS AND OUTFLOWS FOR MASSBAL.OPT FILE
            IF (TRIBUTARIES) THEN
              DO JT=1,JTT
                IF (JB == JBTR(JT)) THEN
                  I = ITR(JT)
                  IF (I < CUS(JB)) I = CUS(JB)
                  DO K=KTTR(JT),KBTR(JT)
                      TPALG=0.0;TNALG=0.0;TPZ=0.0;TNZ=0.0;TPBOD=0.0;TNBOD=0.0
                    IF (QTR(JT) < 0.0) THEN
                        DO JA=1,NAL
                            TPALG=TPALG+ALG(K,I,JA)*AP(JA)
                            TNALG=TNALG+ALG(K,I,JA)*AN(JA)
                        ENDDO
                        DO JCB=1,NBOD
                            TPBOD=TPBOD+CBODP(K,I,JCB)
                            TNBOD=TNBOD+CBODN(K,I,JCB)
                        ENDDO
                        DO JZ = 1,NZP
                            TPZ=TPZ+ZOO(K,I,JZ)*ZP(JZ)
                            TNZ=TNZ+ZOO(K,I,JZ)*ZN(JZ)
                        ENDDO
                      TPOUT = TPOUT-(TPALG+TPBOD+TPZ+PO4(K,I)+LDOMP(K,I)+LPOMP(K,I)+RDOMP(K,I)+RPOMP(K,I))*QTR(JT)*QTRF(K,JT)*DLT/1000.
                      TNOUT = TNOUT-(TNALG+TNBOD+TNZ+NO3(K,I)+NH4(K,I)+LDOMN(K,I)+LPOMN(K,I)+RDOMN(K,I)+RPOMN(K,I))*QTR(JT)*QTRF(K,JT)*DLT/1000.
                    ELSE
                        DO JC=NAS,NAE
                            TPALG=TPALG+CTR(JC,JT)*AP(JC-NAS+1)
                            TNALG=TNALG+CTR(JC,JT)*AN(JC-NAS+1)
                        ENDDO
                        DO JC=NBODS,NBODE,3
                            TPBOD=TPBOD+CTR(JC+1,JT)
                            TNBOD=TNBOD+CTR(JC+2,JT)
                        ENDDO
                        DO JC=NZOOS,NZOOE
                            TPZ=TPZ+CTR(JC,JT)*ZP(JC-NZOOS+1)
                            TNZ=TNZ+CTR(JC,JT)*ZN(JC-NZOOS+1)
                        ENDDO
                      TPTRIB(JW) = TPTRIB(JW)+ (TPALG+TPBOD+TPZ+CTR(NPO4,JT)+CTR(NLDOMP,JT)+CTR(NRDOMP,JT)+CTR(NLPOMP,JT)+CTR(NRPOMP,JT))*QTR(JT)*QTRF(K,JT)*DLT/1000.
                      TNTRIB(JW) = TNTRIB(JW)+ (TNALG+TNBOD+TNZ+CTR(NNH4,JT)+CTR(NNO3,JT)+CTR(NLDOMN,JT)+CTR(NRDOMN,JT)+CTR(NLPOMN,JT)+CTR(NRPOMN,JT))*QTR(JT)*QTRF(K,JT)*DLT/1000.
                    END IF
                  END DO
                END IF
              END DO
            END IF
            IF (DIST_TRIBS(JB)) THEN
              DO I=IU,ID
                TPALG=0.0;TNALG=0.0;TPZ=0.0;TNZ=0.0;TPBOD=0.0;TNBOD=0.0
                IF (QDT(I) < 0.0) THEN
                        DO JA=1,NAL
                            TPALG=TPALG+ALG(KT,I,JA)*AP(JA)
                            TNALG=TNALG+ALG(KT,I,JA)*AN(JA)
                        ENDDO
                        DO JCB=1,NBOD
                            TPBOD=TPBOD+CBODP(KT,I,JCB)
                            TNBOD=TNBOD+CBODN(KT,I,JCB)
                        ENDDO
                        DO JZ = 1,NZP
                            TPZ=TPZ+ZOO(KT,I,JZ)*ZP(JZ)
                            TNZ=TNZ+ZOO(KT,I,JZ)*ZN(JZ)
                        ENDDO
                      TPOUT = TPOUT-(TPALG+TPBOD+TPZ+PO4(KT,I)+LDOMP(KT,I)+LPOMP(KT,I)+RDOMP(KT,I)+RPOMP(KT,I))*QDT(I)*DLT/1000.
                      TNOUT = TNOUT-(TNALG+TNBOD+TNZ+NO3(KT,I)+NH4(KT,I)+LDOMN(KT,I)+LPOMN(KT,I)+RDOMN(KT,I)+RPOMN(KT,I))*QDT(I)*DLT/1000.
                ELSE
                        DO JC=NAS,NAE
                            TPALG=TPALG+CDTR(JC,JB)*AP(JC-NAS+1)
                            TNALG=TNALG+CDTR(JC,JB)*AN(JC-NAS+1)
                        ENDDO
                        DO JC=NBODS,NBODE,3
                            TPBOD=TPBOD+CDTR(JC+1,JB)
                            TNBOD=TNBOD+CDTR(JC+2,JB)
                        ENDDO
                        DO JC=NZOOS,NZOOE
                            TPZ=TPZ+CDTR(JC,JB)*ZP(JC-NZOOS+1)
                            TNZ=TNZ+CDTR(JC,JB)*ZN(JC-NZOOS+1)
                        ENDDO
                      TPDTRIB(JW) = TPDTRIB(JW)+ (TPALG+TPBOD+TPZ+CDTR(NPO4,JB)+CDTR(NLDOMP,JB)+CDTR(NRDOMP,JB)+CDTR(NLPOMP,JB)+CDTR(NRPOMP,JB))*QDT(I)*DLT/1000.
                      TNDTRIB(JW) = TNDTRIB(JW)+ (TNALG+TNBOD+TNZ+CDTR(NNH4,JB)+CDTR(NNO3,JB)+CDTR(NLDOMN,JB)+CDTR(NRDOMN,JB)+CDTR(NLPOMN,JB)+CDTR(NRPOMN,JB))*QDT(I)*DLT/1000.
                END IF
              END DO
            END IF
            IF (WITHDRAWALS) THEN
              DO JWD=1,JWW
                IF (QWD(JWD) /= 0.0) THEN
                  IF (JB == JBWD(JWD)) THEN
                    I = MAX(CUS(JBWD(JWD)),IWD(JWD))
                    DO K=KTW(JWD),KBW(JWD)      !CONCURRENT(K=KTW(JWD):KBW(JWD))                 ! FORALL
                         TPALG=0.0;TNALG=0.0;TPZ=0.0;TNZ=0.0;TPBOD=0.0;TNBOD=0.0
                        DO JC=NAS,NAE
                            TPALG=TPALG+C1S(K,I,JC)*AP(JC-NAS+1)
                            TNALG=TNALG+C1S(K,I,JC)*AN(JC-NAS+1)
                        ENDDO
                        DO JC=NBODS,NBODE,3
                            TPBOD=TPBOD+C1S(K,I,JC)
                            TNBOD=TNBOD+C1S(K,I,JC)
                        ENDDO
                        DO JC=NZOOS,NZOOE
                            TPZ=TPZ+C1S(K,I,JC)*ZP(JC-NZOOS+1)
                            TNZ=TNZ+C1S(K,I,JC)*ZN(JC-NZOOS+1)
                        ENDDO
                      TPWD(JW) = TPWD(JW) + (TPALG+TPBOD+TPZ+C1S(K,I,NPO4)+C1S(K,I,NLDOMP)+C1S(K,I,NRDOMP)+C1S(K,I,NLPOMP)+C1S(K,I,NRPOMP))*QSW(K,JWD)*DLT/1000.
                      TNWD(JW) = TNWD(JW) + (TNALG+TNBOD+TNZ+C1S(K,I,NNH4)+C1S(K,I,NNO3)+C1S(K,I,NLDOMN)+C1S(K,I,NRDOMN)+C1S(K,I,NLPOMN)+C1S(K,I,NRPOMN))*QSW(K,JWD)*DLT/1000.
                    END DO
                  END IF
                END IF
              END DO
            END IF
            IF (PRECIPITATION(JW)) THEN
              DO I=IU,ID    !CONCURRENT (I=IU:ID)                                  !FORALL
                TPALG=0.0;TNALG=0.0;TPZ=0.0;TNZ=0.0;TPBOD=0.0;TNBOD=0.0
                        DO JC=NAS,NAE
                            TPALG=TPALG+CPR(JC,JB)*AP(JC-NAS+1)
                            TNALG=TNALG+CPR(JC,JB)*AN(JC-NAS+1)
                        ENDDO
                        DO JC=NBODS,NBODE,3
                            TPBOD=TPBOD+CPR(JC,JB)
                            TNBOD=TNBOD+CPR(JC,JB)
                        ENDDO
                        DO JC=NZOOS,NZOOE
                            TPZ=TPZ+CPR(JC,JB)*ZP(JC-NZOOS+1)
                            TNZ=TNZ+CPR(JC,JB)*ZN(JC-NZOOS+1)
                        ENDDO
                      TPPR(JW) = TPPR(JW) + (TPALG+TPBOD+TPZ+CPR(NPO4,JB)+CPR(NLDOMP,JB)+CPR(NRDOMP,JB)+CPR(NLPOMP,JB)+CPR(NRPOMP,JB))*QPR(I)*DLT/1000.
                      TNPR(JW) = TNPR(JW) + (TNALG+TNBOD+TNZ+CPR(NNH4,JB)+CPR(NNO3,JB)+CPR(NLDOMN,JB)+CPR(NRDOMN,JB)+CPR(NLPOMN,JB)+CPR(NRPOMN,JB))*QPR(I)*DLT/1000.
              END DO
            END IF
            IF (UP_FLOW(JB)) THEN
              DO K=KT,KB(IU)
              TPALG=0.0;TNALG=0.0;TPZ=0.0;TNZ=0.0;TPBOD=0.0;TNBOD=0.0
                IF (.NOT. HEAD_FLOW(JB)) THEN
                        DO JC=NAS,NAE
                            TPALG=TPALG+CIN(JC,JB)*AP(JC-NAS+1)
                            TNALG=TNALG+CIN(JC,JB)*AN(JC-NAS+1)
                        ENDDO
                        DO JC=NBODS,NBODE,3
                            TPBOD=TPBOD+CIN(JC,JB)
                            TNBOD=TNBOD+CIN(JC,JB)
                        ENDDO
                        DO JC=NZOOS,NZOOE
                            TPZ=TPZ+CIN(JC,JB)*ZP(JC-NZOOS+1)
                            TNZ=TNZ+CIN(JC,JB)*ZN(JC-NZOOS+1)
                        ENDDO
                      TPIN(JW) = TPIN(JW) + (TPALG+TPBOD+TPZ+CIN(NPO4,JB)+CIN(NLDOMP,JB)+CIN(NRDOMP,JB)+CIN(NLPOMP,JB)+CIN(NRPOMP,JB))*QINF(K,JB)*QIN(JB)*DLT/1000.
                      TNIN(JW) = TNIN(JW) + (TNALG+TNBOD+TNZ+CIN(NNH4,JB)+CIN(NNO3,JB)+CIN(NLDOMN,JB)+CIN(NRDOMN,JB)+CIN(NLPOMN,JB)+CIN(NRPOMN,JB))*QINF(K,JB)*QIN(JB)*DLT/1000.
                END IF
              END DO
            END IF
            IF (DN_FLOW(JB)) THEN
                DO K=KT,KB(ID)
                        TPALG=0.0;TNALG=0.0;TPZ=0.0;TNZ=0.0;TPBOD=0.0;TNBOD=0.0
                        DO JC=NAS,NAE
                            TPALG=TPALG+C1S(K,ID,JC)*AP(JC-NAS+1)
                            TNALG=TNALG+C1S(K,ID,JC)*AN(JC-NAS+1)
                        ENDDO
                        DO JC=NBODS,NBODE,3
                            TPBOD=TPBOD+C1S(K,ID,JC)
                            TNBOD=TNBOD+C1S(K,ID,JC)
                        ENDDO
                        DO JC=NZOOS,NZOOE
                            TPZ=TPZ+C1S(K,ID,JC)*ZP(JC-NZOOS+1)
                            TNZ=TNZ+C1S(K,ID,JC)*ZN(JC-NZOOS+1)
                        ENDDO
                      TPOUT(JW) = TPOUT(JW) + (TPALG+TPBOD+TPZ+C1S(K,ID,NPO4)+C1S(K,ID,NLDOMP)+C1S(K,ID,NRDOMP)+C1S(K,ID,NLPOMP)+C1S(K,ID,NRPOMP))*QOUT(K,JB)*DLT/1000.   ! C1S(KT:KB(ID),ID,JC)
                      TNOUT(JW) = TNOUT(JW) + (TNALG+TNBOD+TNZ+C1S(K,ID,NNH4)+C1S(K,ID,NNO3)+C1S(K,ID,NLDOMN)+C1S(K,ID,NRDOMN)+C1S(K,ID,NLPOMN)+C1S(K,ID,NRPOMN))*QOUT(K,JB)*DLT/1000.
                ENDDO
            ENDIF
            
            !IF (UP_HEAD(JB)) THEN
            !    DO K=KT,KB(IU)
            !    IUT = IU
            !    IF (QUH1(K,JB) >= 0.0) IUT = IU-1
            !    !CSSUH1(K,JC,JB) = C1S(K,IUT,JC)*QUH1(K,JB)
            !
            !  END DO
            !  IF (UH_INTERNAL(JB)) THEN
            !    IF (UHS(JB) /= DS(JBUH(JB)) .OR. DHS(JBUH(JB)) /= US(JB)) THEN
            !      IF (JBUH(JB) >= BS(JW) .AND. JBUH(JB) <= BE(JW)) THEN
            !        I = UHS(JB)
            !        DO K=KT,KB(IU)
            !          !CSSB(K,I,JC) = CSSB(K,I,JC)-CSSUH2(K,JC,JB)/DLT
            !        END DO
            !      ELSE
            !        !CALL UPSTREAM_CONSTITUENT(C2(:,:,JC),CSSB(:,:,JC))
            !      END IF
            !    END IF
            !  END IF
            !END IF
            !IF (DN_HEAD(JB)) THEN
            !  DO K=KT,KB(ID+1)
            !    IDT = ID+1
            !    IF (QDH1(K,JB) >= 0.0) IDT = ID
            !    !CSSDH1(K,JC,JB) = C1S(K,IDT,JC)*QDH1(K,JB)
            !
            !  END DO
            !  IF (DH_INTERNAL(JB)) THEN
            !    IF (DHS(JB) /= US(JBDH(JB)) .OR. UHS(JBDH(JB)) /= DS(JB)) THEN
            !      IF (JBDH(JB) >= BS(JW) .AND. JBDH(JB) <= BE(JW)) THEN
            !        I = DHS(JB)
            !        DO K=KT,KB(ID+1)
            !          !CSSB(K,I,JC) = CSSB(K,I,JC)+CSSDH2(K,JC,JB)/DLT
            !        END DO
            !      ELSE
            !        !CALL DOWNSTREAM_CONSTITUENT(C2(:,:,JC),CSSB(:,:,JC))
            !      END IF
            !    END IF
            !  END IF
            !END IF
         ENDIF  ! END OF TP AND TN MASS BALANCES   
         
          
        END DO    ! JB loop
        
        ! Atmospheric Depsition original unit kg/km2/year
        IF(ATM_DEPOSITION(JW))THEN
            DO JB=BS(JW),BE(JW)
                DO I=CUS(JB),DS(JB)
                    DO JAC=1,NACATD(JW)
                            CSSB(KT,I,ATMDCN(JAC,JW))=CSSB(KT,I,ATMDCN(JAC,JW))+ATM_DEP_LOADING(ATMDCN(JAC,JW),JW)*BI(KT,I)*DLX(I)*3.17098E-11   ! Conversion: kg/km2/year to g/m2/s 1000/(365*86400*1000*1000)=3.17098E-11
                            IF(ATMDCN(JAC,JW)==NPO4.OR.ATMDCN(JAC,JW)==NLPOMP.OR.ATMDCN(JAC,JW)==NRPOMP)THEN 
                            ATMDEP_P(JW)=ATMDEP_P(JW)+ATM_DEP_LOADING(ATMDCN(JAC,JW),JW)*BI(KT,I)*DLX(I)*3.17098E-11*DLT/1000.    ! P MASS BALANCE IN KG - CUMULATIVE
                            ELSEIF(ATMDCN(JAC,JW)==NNO3.OR.ATMDCN(JAC,JW)==NLPOMN.OR.ATMDCN(JAC,JW)==NRPOMN.OR.ATMDCN(JAC,JW)==NNH4)THEN 
                            ATMDEP_N(JW)=ATMDEP_N(JW)+ATM_DEP_LOADING(ATMDCN(JAC,JW),JW)*BI(KT,I)*DLX(I)*3.17098E-11*DLT/1000.    ! N MASS BALANCE IN KG - CUMULATIVE
                            ENDIF
                    ENDDO
                ENDDO
            ENDDO
        ENDIF  

        
      END DO      ! JW Loop

!**** Kinetic fluxes

      DO JW=1,NWB
        KT = KTWB(JW)    ! SW 10/25/2017
        IF (FLUX(JW)) CALL KINETIC_FLUXES
      END DO

    !SP CEMA
    IF (UPDATE_KINETICS) THEN
      If(CEMARelatedCode .and. IncludeBedConsolidation) Call CEMASedimentModel
      !If(CEMARelatedCode .and. IncludeCEMASedDiagenesis)Call SedimentFlux
      If(CEMARelatedCode .and. IncludeCEMASedDiagenesis .AND. BUBBLES_CALCULATION)Call CEMACalculateRiseVelocity
      If(CEMARelatedCode .and. IncludeCEMASedDiagenesis .and. ApplyBubbTurb .AND. BUBBLES_CALCULATION)Call CEMABubblesReleaseTurbulence
      If(CEMARelatedCode .and. IncludeCEMASedDiagenesis .AND. BUBBLES_CALCULATION)Call CEMABubblesTransport
      If(CEMARelatedCode .and. IncludeCEMASedDiagenesis .AND. BUBBLES_CALCULATION)Call CEMABubbWatTransfer
      If(CEMARelatedCode .and. IncludeCEMASedDiagenesis .AND. BUBBLES_CALCULATION)Call CEMABubblesRelease
      If(IncludeFFTLayer)Call CEMAFFTLayerCode
    END IF
    !End SP CEMA

!**** Constituent transport

!!$OMP PARALLEL DO PRIVATE(I,JC,KT,JB,JW,DT,K,IU,ID,BTA1,GMA1)    !I,JC,KT,JW,JB,CNEW,SSB,SSK,COLD,AT,VT,CT,DT) 
    
    DO JAC=1,NAC    !CONCURRENT(JAC=1:NAC)              !JAC=1,NAC
    JC   =  CN(JAC)
    COLD => C1S(:,:,JC)
      DO JW=1,NWB
        KT = KTWB(JW)
        DO JB=BS(JW),BE(JW)
        IF(BR_INACTIVE(JB))CYCLE   ! SW 6/12/2017
          IU = CUS(JB)
          ID = DS(JB)
      !    DO JAC=1,NAC
      !      JC   =  CN(JAC)
            !COLD => C1S(:,:,JC)
            CALL HORIZONTAL_MULTIPLIERS       
            CALL VERTICAL_MULTIPLIERS           
            
     !       CNEW => C1(:,:,JC)
     !       SSB  => CSSB(:,:,JC)
     !       SSK  => CSSK(:,:,JC)
     !       CALL HORIZONTAL_TRANSPORT
            DO I=IU,ID
              DO K=KT,KB(I)      
              DT(K,I) = (C1S(K,I,JC)*BH2(K,I)/DLT+(ADX(K,I)*BHR1(K,I)-ADX(K,I-1)*BHR1(K,I-1))/DLX(I)+(1.0D0-THETA(JW))                     &
                    *(ADZ(K,I)*BB(K,I)-ADZ(K-1,I)*BB(K-1,I))+CSSB(K,I,JC)/DLX(I))*DLT/BH1(K,I)+CSSK(K,I,JC)*DLT
              END DO                                                    
            END DO
            DO I=IU,ID
        !      CALL TRIDIAG(AT(:,I),VT(:,I),CT(:,I),DT(:,I),KT,KB(I),KMX,CNEW(:,I))
                    BTA1(KT)=VT(KT,I)
                    GMA1(KT)=DT(KT,I)
                  DO K=KT+1,KB(I)
                           BTA1(K) = VT(K,I)-AT(K,I)/BTA1(K-1)*CT(K-1,I)
                    GMA1(K) = DT(K,I)-AT(K,I)/BTA1(K-1)*GMA1(K-1)
                  END DO
                  C1(KB(I),I,JC) = GMA1(KB(I))/BTA1(KB(I))
                  DO K=KB(I)-1,KT,-1
                    C1(K,I,JC) = (GMA1(K)-CT(K,I)*C1(K+1,I,JC))/BTA1(K)
                  END DO
            END DO
          END DO
        END DO
    END DO
!!$OMP END PARALLEL DO
      IF (DERIVED_CALC) CALL DERIVED_CONSTITUENTS

    END SUBROUTINE WQCONSTITUENTS
