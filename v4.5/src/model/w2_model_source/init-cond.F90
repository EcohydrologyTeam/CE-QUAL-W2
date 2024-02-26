!***********************************************************************************************************************************
!**                                              Task 1.4.4: Initial conditions                                                   **
!***********************************************************************************************************************************
SUBROUTINE INITCOND
USE MAIN
USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  USE MACROPHYTEC; USE POROSITYC; USE ZOOPLANKTONC
  Use CEMAVars, ONLY:CEMARelatedCode,IncludeCEMASedDiagenesis                         ! cb 07/23/18
  USE CEMASedimentDiagenesis
  IMPLICIT NONE
  EXTERNAL RESTART_OUTPUT

  REAL          :: TMAC,XSAR
  CHARACTER(1)  :: ICHAR
  CHARACTER(8)  :: IBLANK

   BIC=B

  DO JW=1,NWB
    KT = KTWB(JW)
    IF (VERT_PROFILE(JW)) THEN

!**** Temperature and water quality

      OPEN (VPR(JW),FILE=VPRFN(JW),STATUS='OLD')
      READ (VPR(JW),'(A1)')ICHAR
 
      IF(ICHAR=='$')THEN
         READ( VPR(JW),'(/)')
         IF (VERT_TEMP(JW)) READ (VPR(JW),*) IBLANK, (TVP(K,JW),K=KT,KBMAX(JW))
         IF (CONSTITUENTS) THEN
          DO JC=1,NCT
           IF (VERT_CONC(JC,JW))      READ (VPR(JW),*) IBLANK, (CVP(K,JC,JW),  K=KT,KBMAX(JW))
          END DO
          DO JE=1,NEP
            IF (VERT_EPIPHYTON(JW,JE)) READ (VPR(JW),*) IBLANK, (EPIVP(K,JW,JE),K=KT,KBMAX(JW))
          END DO
          DO m=1,nmc                                             ! cb 8/21/15
            IF (VERT_macrophyte(JW,m)) READ (VPR(JW),*) IBLANK, (macrcvp(K,JW,m),K=KT,KBMAX(JW))
          END DO
          IF (VERT_SEDIMENT(JW))       READ (VPR(JW),*) IBLANK,(SEDVP(K,JW),   K=KT,KBMAX(JW))
          END IF
      ELSE
         IF (VERT_TEMP(JW)) READ (VPR(JW),'(//(8X,9F8.0))') (TVP(K,JW),K=KT,KBMAX(JW))
         IF (CONSTITUENTS) THEN
           DO JC=1,NCT
            IF (VERT_CONC(JC,JW))      READ (VPR(JW),'(//(8X,9F8.0))') (CVP(K,JC,JW),  K=KT,KBMAX(JW))
           END DO
           DO JE=1,NEP
            IF (VERT_EPIPHYTON(JW,JE)) READ (VPR(JW),'(//(8X,9F8.0))') (EPIVP(K,JW,JE),K=KT,KBMAX(JW))
           END DO
           DO m=1,nmc                                              ! cb 8/21/15
            IF (VERT_macrophyte(JW,m)) READ (VPR(JW),'(//(8X,9F8.0))') (macrcvp(K,JW,m),K=KT,KBMAX(JW))
           END DO
           IF (VERT_SEDIMENT(JW))       READ (VPR(JW),'(//(8X,9F8.0))') (SEDVP(K,JW),   K=KT,KBMAX(JW))
          END IF
      ENDIF
    END IF

!** Longitudinal/vertical initial profiles

    IF (LONG_PROFILE(JW)) THEN
       OPEN (LPR(JW),FILE=LPRFN(JW),STATUS='OLD')
       READ (LPR(JW),'(A1)')ICHAR
       IF(ICHAR=='$')READ(LPR(JW),*)
    END IF

!** Branch related variables

    IF (.NOT. RESTART_IN) THEN
        IF(LONG_TEMP(JW).AND.ICHAR=='$')READ (LPR(JW),*)
      DO JB=BS(JW),BE(JW)

!****** Temperature

        DO I=CUS(JB),DS(JB)
          IF (LONG_TEMP(JW)) THEN
              IF(ICHAR=='$')THEN
              READ (LPR(JW),*) IBLANK,(T1(K,I),K=KT,KB(I))
              ELSE
              READ (LPR(JW),'(//(8X,9F8.0))') (T1(K,I),K=KT,KB(I))
              ENDIF
          ENDIF
          DO K=KT,KB(I)
            IF (ISO_TEMP(JW))  T1(K,I) = T2I(JW)
            IF (VERT_TEMP(JW)) T1(K,I) = TVP(K,JW)
            T2(K,I) = T1(K,I)
          END DO
        END DO
      END DO

!**** Constituents

      DO JC=1,NAC
        IF (LONG_CONC(CN(JC),JW).AND.ICHAR=='$')READ (LPR(JW),*)
        DO JB=BS(JW),BE(JW)
          DO I=CUS(JB),DS(JB)
            JAC = CN(JC)
            IF (LONG_CONC(JAC,JW))THEN
                IF(ICHAR=='$')THEN
                READ (LPR(JW),*) IBLANK,(C2(K,I,JAC),K=KT,KB(I))
                ELSE
                READ (LPR(JW),'(//(8X,9F8.0))') (C2(K,I,JAC),K=KT,KB(I))
                ENDIF
            ENDIF
            DO K=KT,KB(I)
              IF (ISO_CONC(JAC,JW))  C2(K,I,JAC) = C2I(JAC,JW)
              IF (VERT_CONC(JAC,JW)) C2(K,I,JAC) = CVP(K,JAC,JW)
              C1(K,I,JAC)  = C2(K,I,JAC)
              C1S(K,I,JAC) = C1(K,I,JAC)
            END DO
          END DO
        END DO
      END DO

!**** Epiphyton


    DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE)) THEN
        IF (LONG_EPIPHYTON(JW,JE).AND.ICHAR=='$')READ (LPR(JW),*)
           DO JB=BS(JW),BE(JW) 
           DO I=CUS(JB),DS(JB)
              IF (LONG_EPIPHYTON(JW,JE))THEN
                  IF(ICHAR=='$')THEN
                    READ (LPR(JW),*) IBLANK,(EPD(K,I,JE),K=KT,KB(I))
                      ELSE
                    READ (LPR(JW),'(//(8X,9F8.0))') (EPD(K,I,JE),K=KT,KB(I))
                ENDIF
              ENDIF
              IF (ISO_EPIPHYTON(JW,JE))  EPD(:,I,JE) = EPICI(JW,JE)
              IF (VERT_EPIPHYTON(JW,JE)) EPD(:,I,JE) = EPIVP(:,JW,JE)                                    ! CB 5/16/2009
            END DO
           END DO
        END IF
    END DO
    
!**** macrophytes - added 8/21/15


    DO m=1,nmc
        IF (macrophyte_CALC(JW,m)) THEN
        IF (LONG_macrophyte(JW,m).AND.ICHAR=='$')READ (LPR(JW),*)
           DO JB=BS(JW),BE(JW) 
           DO I=CUS(JB),DS(JB)
              IF (LONG_macrophyte(JW,m))THEN
                  IF(ICHAR=='$')THEN
                    READ (LPR(JW),*) IBLANK,(macrclp(K,I,m),K=KT,KB(I))
                      ELSE
                    READ (LPR(JW),'(//(8X,9F8.0))') (macrclp(K,I,m),K=KT,KB(I))
                ENDIF
              ENDIF
              !IF (ISO_macrophyte(JW,m))  macrc(:,I,m) = macwbci(JW,m)
              !IF (VERT_macrophyte(JW,m)) macrc(:,I,m) = macrcvp(:,JW,m)    
            END DO
           END DO
        END IF
      END DO

!**** Sediments

      DO JB=BS(JW),BE(JW)
 !     SDKV(:,US(JB):DS(JB))=SDK(JW)
        SDKV(:,US(JB)-1:DS(JB)+1)=SDK(JW)    ! SW 9/28/13
        IF (SEDIMENT_CALC(JW)) THEN
            IF(LONG_SEDIMENT(JW).AND.JB==BS(JW))READ (LPR(JW),*)
          DO I=CUS(JB),DS(JB)
            IF (LONG_SEDIMENT(JW))THEN
                IF(ICHAR=='$')THEN
                    READ (LPR(JW),*)IBLANK, (SED(K,I),K=KT,KB(I)) 
                    ELSE
                    READ (LPR(JW),'(//(8X,9F8.0))') (SED(K,I),K=KT,KB(I))
                ENDIF
            ENDIF
            DO K=KT,KB(I)
              IF (ISO_SEDIMENT(JW))  SED(K,I) = SEDCI(JW)
              IF (VERT_SEDIMENT(JW)) SED(K,I) = SEDVP(K,JW)
            END DO
            SED(KT,I)         = SED(KT,I)/H2(KT,I)
            If(CEMARelatedCode .and. IncludeCEMASedDiagenesis)then   ! cb 07/23/18
              SED(KT+1:KB(I)-1,I) = SED(KT+1:KB(I)-1,I)/H2(KT+1:KB(I)-1,I)             
              sed(kb(i),i)=0.0
            else              
            SED(KT+1:KB(I),I) = SED(KT+1:KB(I),I)/H2(KT+1:KB(I),I)             
            end if
          END DO
        END IF
      END DO
        DO JB=BS(JW),BE(JW)
          IF (SEDIMENT_CALC(JW)) THEN
            DO I=CUS(JB),DS(JB)
              DO K=KT,KB(I)
                IF (ISO_SEDIMENT(JW))SEDP(K,I) = ORGP(JW)*SEDCI(JW)
                IF (VERT_SEDIMENT(JW))SEDP(K,I) = SEDVP(K,JW)*ORGP(JW)
                IF (LONG_SEDIMENT(JW)) SEDP(K,I)=ORGP(JW)*SED(K,I)
              END DO
              SEDP(KT,I)         = SEDP(KT,I)/H2(KT,I)
              If(CEMARelatedCode .and. IncludeCEMASedDiagenesis)then   ! cb 07/23/18
                SEDp(KT+1:KB(I)-1,I) = SEDp(KT+1:KB(I)-1,I)/H2(KT+1:KB(I)-1,I)             
                sedp(kb(i),i)=0.0
              else                
              SEDP(KT+1:KB(I),I) = SEDP(KT+1:KB(I),I)/H2(KT+1:KB(I),I)
              end if
            END DO
          END IF
        END DO
        DO JB=BS(JW),BE(JW)
          IF (SEDIMENT_CALC(JW)) THEN
            DO I=CUS(JB),DS(JB)
              DO K=KT,KB(I)
                IF (ISO_SEDIMENT(JW))SEDn(K,I) = orgn(JW)*sedci(jw)
                IF (VERT_SEDIMENT(JW))SEDn(K,I) = SEDVP(K,JW)*orgn(jw)
                IF (LONG_SEDIMENT(JW)) sedn(k,i)=orgn(jw)*sed(k,i)
              END DO
              SEDn(KT,I)         = SEDn(KT,I)/H2(KT,I)
              If(CEMARelatedCode .and. IncludeCEMASedDiagenesis)then   ! cb 07/23/18                
                SEDn(KT+1:KB(I)-1,I) = SEDn(KT+1:KB(I)-1,I)/H2(KT+1:KB(I)-1,I)
                sedn(kb(i),i)=0.0
              else
              SEDn(KT+1:KB(I),I) = SEDn(KT+1:KB(I),I)/H2(KT+1:KB(I),I)
              end if
          END DO
        END IF
      END DO
        DO JB=BS(JW),BE(JW)
          IF (SEDIMENT_CALC(JW)) THEN
            DO I=CUS(JB),DS(JB)
              DO K=KT,KB(I)
                IF (ISO_SEDIMENT(JW))SEDc(K,I) = SEDCI(JW)*orgc(jw)
                IF (VERT_SEDIMENT(JW))SEDc(K,I) = SEDVP(K,JW)*orgc(jw)
                IF (LONG_SEDIMENT(JW)) sedc(k,i)=orgc(jw)*sed(k,i)
              END DO
              SEDc(KT,I)         = SEDc(KT,I)/H2(KT,I)
              If(CEMARelatedCode .and. IncludeCEMASedDiagenesis)then   ! cb 07/23/18                
                SEDc(KT+1:KB(I)-1,I) = SEDc(KT+1:KB(I)-1,I)/H2(KT+1:KB(I)-1,I)
                sedc(kb(i),i)=0.0
              else
              SEDc(KT+1:KB(I),I) = SEDc(KT+1:KB(I),I)/H2(KT+1:KB(I),I)
              end if
            END DO
          END IF
        END DO

      SED(:,US(BS(JW)):DS(BE(JW))) = SED(:,US(BS(JW)):DS(BE(JW)))*FSED(JW)
      SEDp(:,US(BS(JW)):DS(BE(JW))) = SEDp(:,US(BS(JW)):DS(BE(JW)))*FSED(JW)
      SEDn(:,US(BS(JW)):DS(BE(JW))) = SEDn(:,US(BS(JW)):DS(BE(JW)))*FSED(JW)
      SEDc(:,US(BS(JW)):DS(BE(JW))) = SEDc(:,US(BS(JW)):DS(BE(JW)))*FSED(JW)

!  Amaila start Additional sediment compartments
      DO JB=BS(JW),BE(JW)      
        IF (SEDIMENT_CALC1(JW)) THEN
            !IF(LONG_SEDIMENT(JW).AND.JB==BS(JW))READ (LPR(JW),*)
            IF(LONG_SEDIMENT(JW).AND.ICHAR=='$')READ (LPR(JW),*)      ! cb 6/10/13
          !DO I=CUS(JB),DS(JB)
           DO I=US(JB),DS(JB)              ! cb 6/17/17
            IF (LONG_SEDIMENT1(JW))THEN
                IF(ICHAR=='$')THEN
                    !READ (LPR(JW),*)IBLANK, (SED1(K,I),K=KT,KB(I)) 
                    READ (LPR(JW),*)IBLANK, (SED1(K,I),K=2,KB(I))   ! cb 6/17/17
                    ELSE
                    !READ (LPR(JW),'(//(8X,9F8.0))') (SED1(K,I),K=KT,KB(I))
                    READ (LPR(JW),'(//(8X,9F8.0))') (SED1(K,I),K=2,KB(I))  ! cb 6/17/17
                ENDIF
            ENDIF
            DO K=KT,KB(I)
              IF (ISO_SEDIMENT1(JW))  SED1(K,I) = SEDCI1(JW)
              IF (VERT_SEDIMENT1(JW)) SED1(K,I) = SEDVP1(K,JW)
            END DO
            !SED1(KT,I)         = SED1(KT,I)/H2(KT,I)   ! intial conditions for "tree" sediment compartments are given in g/m^3
            !SED1(KT+1:KB(I),I) = SED1(KT+1:KB(I),I)/H2(KT+1:KB(I),I)             
          END DO
        END IF
      END DO
      
       DO JB=BS(JW),BE(JW)      
        IF (SEDIMENT_CALC2(JW)) THEN
            !IF(LONG_SEDIMENT(JW).AND.JB==BS(JW))READ (LPR(JW),*)
            IF(LONG_SEDIMENT(JW).AND.ICHAR=='$')READ (LPR(JW),*)      ! cb 6/10/13
          !DO I=CUS(JB),DS(JB)
          DO I=US(JB),DS(JB)              ! cb 6/17/17
            IF (LONG_SEDIMENT2(JW))THEN
                IF(ICHAR=='$')THEN
                    !READ (LPR(JW),*)IBLANK, (SED2(K,I),K=KT,KB(I)) 
                    READ (LPR(JW),*)IBLANK, (SED2(K,I),K=2,KB(I))   ! cb 6/17/17
                    ELSE
                    !READ (LPR(JW),'(//(8X,9F8.0))') (SED2(K,I),K=KT,KB(I))
                    READ (LPR(JW),'(//(8X,9F8.0))') (SED2(K,I),K=2,KB(I))  ! cb 6/17/17
                ENDIF
            ENDIF
            DO K=KT,KB(I)
              IF (ISO_SEDIMENT2(JW))  SED2(K,I) = SEDCI2(JW)
              IF (VERT_SEDIMENT2(JW)) SED2(K,I) = SEDVP2(K,JW)
            END DO
            !SED2(KT,I)         = SED2(KT,I)/H2(KT,I)      ! intial conditions for "tree" sediment compartments are given in g/m^3
            !SED2(KT+1:KB(I),I) = SED2(KT+1:KB(I),I)/H2(KT+1:KB(I),I)             
          END DO
        END IF
      END DO

       DO JB=BS(JW),BE(JW)    ! 9/3/17      
         DO I=US(JB),DS(JB)
           do k=kt,kb(i)
             sdfirstadd(k,i)=.false.   
           end do
         end do
       end do

       SED1(:,US(BS(JW)):DS(BE(JW))) = SED1(:,US(BS(JW)):DS(BE(JW)))*FSEDc1(JW)  ! cb 6/7/17
       SED2(:,US(BS(JW)):DS(BE(JW))) = SED2(:,US(BS(JW)):DS(BE(JW)))*FSEDc2(JW)       
       sed1ic=sed1     ! cb 6/17/17
       sed2ic=sed2     ! cb 6/17/17


! Amaila end

      DO JB=BS(JW),BE(JW)
        DO M=1,NMC
          IF (MACROPHYTE_CALC(JW,M)) THEN

!C DISTRIBUTING INITIAL MACROPHYTE CONC TO BOTTOM COLUMN CELLS; MACWBCI = G/M^3
            DO I=CUS(JB),DS(JB)

              DEPKTI=ELWS(I)-EL(KTI(I)+1,I)

              IF(DEPKTI.GE.THRKTI)THEN
                KTICOL(I)=.TRUE.
                JT=KTI(I)
              ELSE
                KTICOL(I)=.FALSE.
                JT=KTI(I)+1
              END IF

              JE=KB(I)
              DO J=JT,JE
                IF(J.LE.KT)THEN
                  K=KT
                ELSE
                  K=J
                END IF
                !MACRC(J,K,I,M) = MACWBCI(JW,M)
                !SMACRC(J,K,I,M) = MACWBCI(JW,M)                
                IF (ISO_macrophyte(JW,m))  macrc(j,k,I,m) = macwbci(JW,m)     ! cb 8/24/15
                IF (VERT_macrophyte(JW,m)) macrc(j,k,I,m) = macrcvp(k,JW,m)    
                IF (long_macrophyte(JW,m)) macrc(j,k,I,m) = macrclp(K,I,m)
                SMACRC(J,K,I,M) = macrc(j,k,I,m)
              END DO
            END DO

            DO I=CUS(JB),DS(JB)
              TMAC=0.0
              XSAR=0.0
              DO K=KTI(I),KT
                JT=K
                JE=KB(I)
                COLB=EL(K+1,I)
                COLDEP=ELWS(I)-COLB
                DO J=JT,JE
                  TMAC=TMAC+MACRC(J,KT,I,M)*CW(J,I)*COLDEP
                  XSAR=XSAR+CW(J,I)*COLDEP
                END DO
              END DO
              MAC(KT,I,M)=TMAC/XSAR
              SMAC(KT,I,M)=MAC(KT,I,M)

              DO K=KT+1,KB(I)
                JT=K
                JE=KB(I)
                TMAC=0.0
                DO J=JT,JE
                  TMAC=TMAC+MACRC(J,K,I,M)*CW(J,I)
                END DO
                MAC(K,I,M)=TMAC/B(K,I)
                SMAC(K,I,M)=MAC(K,I,M)
              END DO
            END DO

            DO I=CUS(JB),DS(JB)
              JT=KTI(I)
              JE=KB(I)
              DO J=JT,JE
                IF(J.LT.KT)THEN
                  COLB=EL(J+1,I)
                ELSE
                  COLB=EL(KT+1,I)
                END IF
                COLDEP=ELWS(I)-COLB
                MACRM(J,KT,I,M)=MACRC(J,KT,I,M)*COLDEP*CW(J,I)*DLX(I)
                SMACRM(J,KT,I,M)=MACRM(J,KT,I,M)
              END DO

              DO K=KT+1,KB(I)

                JT=K
                JE=KB(I)

                DO J=JT,JE

                  MACRM(J,K,I,M)=MACRC(J,K,I,M)*H2(K,I)*CW(J,I)*DLX(I)
                  SMACRM(J,K,I,M)=MACRM(J,K,I,M)
                END DO

              END DO
            END DO

          END IF
        END DO
      END DO
! V3.5 END

!**** ENERGY

      DO JB=BS(JW),BE(JW)
        IF(BR_INACTIVE(JB))CYCLE    ! SW 12/18/2018
        DO I=CUS(JB),DS(JB)
          IF (ENERGY_BALANCE(JW)) THEN
            DO K=KT,KB(I)
              EBRI(JB) = EBRI(JB)+T2(K,I)*DLX(I)*BH2(K,I)
            END DO
          END IF
          DO K=KT,KB(I)
            CMBRT(CN(1:NAC),JB) = CMBRT(CN(1:NAC),JB)+C2(K,I,CN(1:NAC))*DLX(I)*BH2(K,I)
          END DO
        END DO

! V3.5 START
!C   INITIALIZING MACROPHYTE TEMPORAL MASS BALANCE TERM....
          DO M=1,NMC
            IF(MACROPHYTE_CALC(JW,M))THEN
              DO I=CUS(JB),DS(JB)
                IF(KTICOL(I))THEN
                  JT=KTI(I)
                ELSE
                  JT=KTI(I)+1
                END IF
                JE=KB(I)
                DO J=JT,JE
                  MACMBRT(JB,M) = MACMBRT(JB,M)+MACRM(J,KT,I,M)
                END DO
                DO K=KT+1,KB(I)
                  JT=K
                  JE=KB(I)
                  DO J=JT,JE
                    MACMBRT(JB,M) = MACMBRT(JB,M)+MACRM(J,K,I,M)
                  END DO
                END DO
              END DO
            END IF
          END DO

!****** Ice cover

        IF (ICE_CALC(JW)) THEN
          ICETH(CUS(JB):DS(JB)) = ICETHI(JW)               ! SW 9/29/15 only initialize CUS to DS
          ICE(US(JB):DS(JB))   = ICETH(US(JB):DS(JB)) > 0.0
          DO I=CUS(JB),DS(JB)                              ! SW 9/29/15
          ICEBANK(I)= ICETH(I)*BI(KT,I)*DLX(I)             ! Initial volume of ice in m3
          ENDDO
        END IF

!****** Vertical eddy viscosity

        IUT = CUS(JB)
        IDT = DS(JB)-1
        IF (UP_HEAD(JB)) IUT = IU-1
        IF (DN_HEAD(JB)) IDT = ID
        DO I=IUT,IDT
          DO K=KT,KB(I)-1
            AZ(K,I)    = AZMIN
            TKE(K,I,1) = 1.25E-7
            TKE(K,I,2) = 1.0E-9
          END DO
        END DO
        DO JWR=1,NIW
          IF (WEIR_CALC) AZ(MAX(KT,KTWR(JWR)-1):KBWR(JWR),IWR(JWR)) = 0.0
        END DO
      END DO
    END IF

!** Horizontal diffusivities

    DO JB=BS(JW),BE(JW)
      DO I=CUS(JB),DS(JB)-1
        DO K=KT,KBMIN(I)
          DX(K,I) = ABS(DXI(JW))    ! SW 8/2/2017 FIRST TIME STEP EVEN IF NEGATIVE USE AS ABS OF DX SINCE IT WILL ALWAYS BE LESS THAN 1     
          IF (INTERNAL_WEIR(K,I)) DX(K,I) = 0.0
        END DO
      END DO
    END DO
    IF (VERT_PROFILE(JW)) CLOSE (VPR(JW))
    IF (LONG_PROFILE(JW)) CLOSE (LPR(JW))   
  END DO

! Atmospheric pressure
    IF (CONSTITUENTS)THEN
        PALT(:) = (1.0-ELWS(:)/1000.0/44.3)**5.25       ! SW 2/3/08
        YEAROLD=YEAR
        IF(CO2YEARLYPPM=='      ON')THEN
            IF(YEAR<1980)THEN
             PCO2 = (0.000041392*REAL(YEAR*YEAR*YEAR) - 0.231409975*REAL(YEAR*YEAR) + 430.804190829*REAL(YEAR) - 266735.857433224)*PALT(DS(BE(1)))*1.0E-6      ! PPM CO2 AND ALTITUDE CORRECTION 
            ELSE
             PCO2  = (0.015903*YEAR**2 - 61.799598*YEAR + 60357.055057)*PALT(DS(BE(1)))*1.0E-6
            ENDIF
        ELSE
        PCO2=PCO2ATMPPM*PALT(DS(BE(1)))*1.0E-6    ! IN ATM
        ENDIF
    ENDIF
    
    ELWS_INI(:) = ELWS(:)         ! systdg - Add ELWS_INI
  RETURN

  END SUBROUTINE INITCOND
