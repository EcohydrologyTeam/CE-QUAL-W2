
!***********************************************************************************************************************************
!**                                            S U B R O U T I N E   K I N E T I C S                                              **
!***********************************************************************************************************************************

SUBROUTINE KINETICS
  USE SCREENC; USE GLOBAL; USE KINETIC; USE GEOMC; USE TVDC; USE LOGICC; USE SURFHE
  USE MACROPHYTEC; USE ZOOPLANKTONC; USE MAIN, ONLY:NPBALC, EPIPHYTON_CALC, BOD_CALC, &
      ALG_CALC, BOD_CALCN, BOD_CALCP, PO4_CALC, N_CALC, DSI_CALC, STANDING_BIOMASS_DECAY, NH3_DER, & 
      CDWBC,KF_NH4_SR,KF_NH4_SD,KF_PO4_SR,KF_PO4_SD,NLDOM, NRDOM, NLPOM, NRPOM, NDGP, ORGC_CALC, CO2_DER, HCO3_DER, CO3_DER,  &
      CBODU_DER,TOTSS_DER,O2DG_DER,TURB_DER,SECCHI_DER, CHLA_DER, GAS_TRANSFER_UPDATE
  USE ALGAE_TOXINS
  Use CEMAVars

! Type declarations
  IMPLICIT NONE
  
  REAL                                :: LAM1,   LAM2,   NH4PR,  NO3PR,  LIMIT,  LIGHT,  L, L0, L1, EA, N2SAT  ! SW 10/17/15
  REAL                                :: KW,     INCR,   OH,     K1,     K2, bicart, DLT13
  REAL                                :: CART,ALKT,T1K,S2,SQRS2,DH1,DH2,H2CO3T,CO3T,PHT,F,HION,HCO3T
  REAL                                :: LTCOEFM, LAVG,  MACEXT, TMAC,MACEXT1         ! CB 4/20/11
  REAL                                :: FETCH, U2, COEF1,COEF2,COEF3,COEF4,HS,TS,COEF,UORB,TAU
  REAL                                :: EPSILON, CBODSET, DOSAT,O2EX,CO2EX,SEDSI,SEDEM, SEDSO,SEDSIP
  REAL                                :: SEDSOP,SEDSON,SEDSOC,SEDSIC,SEDSIDK,SEDSUM,SEDSUMK,XDUM
  REAL                                :: BLIM, SEDSIN, COLB,COLDEP,BMASS,BMASSTEST,CVOL
  REAL                                :: ALGEX, SSEXT, TOTMAC, ZOOEXT, TOTSS0, FDPO4, ZMINFAC, SSR
  REAL                                :: ZGZTOT,CBODCT,CBODNT,CBODPT,BODTOT  ! CB 6/6/10
  REAL                                :: ALGP,ALGN,ZOOP,ZOON,TPSS,XX,KHCO2
  REAL                                :: RGAS=0.00008206, KH_NH3, K_NH3_NH4, K_NH3     ! atm m3/oK/mole for NH3 gas calculation
  ! enhanced pH buffering start
  real                                :: ammt,phost,omct,dh3,dhh,po4t,ht,hpo4t,h2po4t,oht
  real                                :: nh4t,nh3t,h3po4t,kamm,kp1,kp2,kp3
  ! enhanced pH buffering end
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: OMTRM,  SODTRM, NH4TRM, NO3TRM, BIBH2
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: DOM,    POM,    PO4BOD, NH4BOD, TICBOD
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: LAM2M  
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ATRM,   ATRMR,  ATRMF
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ETRM,   ETRMR,  ETRMF
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ASETTLE, DEN_AVG, DENP, DEN1, DEN2, ALLIM_OLD  ! CO 6/9/2019
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: DEN  ! CO 6/9/2019
  REAL, ALLOCATABLE, DIMENSION(:)     :: AMP, PHASE, C_COEFF_EXT, RAD, MIND, MAXD, DENSI, DENBI, T_DEC, C_DENINC, C_DENDEC, DEPTH_LIM, LOSS_FRAC  ! CO 6/3/2019
  REAL, ALLOCATABLE, DIMENSION(:)     :: I_C, C_DENINC_1, C_DENINC_2, C_DENDEC_1, C_DENDEC_2, DENP_MINS, DENP_MINB, DENP_MIN, DEN_COR, EXP_DEPTH    ! CO 6/10/2019
  REAL, ALLOCATABLE, DIMENSION(:,:)   :: MIGON, MIGOFF, LOLD, TWQ    ! CO 6/12/2019
  INTEGER, ALLOCATABLE, DIMENSION(:)  :: MIGRATE_GROUP, MIGRATE_MODEL, TS_DEC, DEPTH_LIM_ONOFF, NMINT, DEN_USE, DEPTH_CALC_ONOFF    ! CO 6/3/2019
  LOGICAL, ALLOCATABLE, DIMENSION(:)  :: ALGAE_SETTLING    ! CO 6/5/2019  
  INTEGER                             :: K, JA, JE, M, JS, JT, JJ, JJZ, JG, JCB, JBOD, LLM,J,JD,LL
  INTEGER                             :: MI,JAF,N,ITER,IBOD,ISETTLE

  integer                             :: jcg
  INTEGER                             :: ALGMIGRATION_DEBUG, NMIG, NITWQ, MIGI, KK ! CO 6/9/2019
  real                                :: ticch4, DK1_BACT, SET_BACT, PHOTO_BACT, KADG, SODTRMDDO3, OMTRMDO3                                              !h2sex,ch4ex,ticch4
  real                                :: sdalgc,sdepc,sdbodc,sdalgn,sdepn,sdbodn,sdalgp,sdepp,sdbodp
  REAL                                :: AIN, AINA, AINB, AOUT, AOUTA, AOUTB, AVG_LIGHT, LAM3, LAM4, VISCK, TOLD, iday    ! co 6/10/2019
  logical                             :: FeMn, ZOOP_SETTLING_EXIST 
  CHARACTER(2)                        :: MIGRATION  ! CO 6/3/2019
  SAVE

! Allocation declarations

  ALLOCATE (OMTRM(KMX,IMX),    SODTRM(KMX,IMX),    NH4TRM(KMX,IMX),    NO3TRM(KMX,IMX), DOM(KMX,IMX), POM(KMX,IMX))
  ALLOCATE (PO4BOD(KMX,IMX),   NH4BOD(KMX,IMX),    TICBOD(KMX,IMX))
  ALLOCATE (ATRM(KMX,IMX,NAL), ATRMR(KMX,IMX,NAL), ATRMF(KMX,IMX,NAL))
  ALLOCATE (ETRM(KMX,IMX,NEP), ETRMR(KMX,IMX,NEP), ETRMF(KMX,IMX,NEP))
  ALLOCATE (lam2m(KMX,kmx),    BIBH2(KMX,IMX), ALGAE_SETTLING(NAL))
  ALLOCATE (FE(KMX,IMX))       
  TICBOD=0.0; FE=0.0
  
  !ZS=0.0    ! SW 1/29/2019
  !ZSR=0.0
  !ZOOP_SETTLING_EXIST=.FALSE.
  !INQUIRE(FILE='zoop_settling.csv',EXIST=ZOOP_SETTLING_EXIST)    ! SW 1/28/2019
  !IF(ZOOP_SETTLING_EXIST)THEN
  !    OPEN(2450,FILE='zoop_settling.csv',STATUS='OLD')
  !    READ(2450,*)  ! SKIP HEADER
  !    READ(2450,*)(ZS(JZ),JZ=1,NZP)   ! zooplankton settling rate in m/day
  !    DO JZ=1,NZP
  !        ZS(JZ)=ZS(JZ)/86400.   ! CONVERT TO M/S
  !    ENDDO
  !    CLOSE(2450)
  ! ENDIF
  !ALGAE_SETTLING_EXIST=.FALSE.  ! CO 6/4/2019
  ALGAE_SETTLING(:)=.FALSE.  ! CO 6/5/2019
  NITWQ=0
  !INQUIRE(FILE='w2_AlgaeMigration.csv',EXIST=ALGAE_SETTLING_EXIST)    
  IF(ALGAE_SETTLING_EXIST)THEN
      OPEN(2450,FILE='w2_AlgaeMigration.csv',STATUS='OLD')
      READ(2450,*)  ! SKIP HEADER
      READ(2450,*)MIGRATION     ! '(A2)'
      IF(MIGRATION /= 'ON')GO TO 100
      READ(2450,*)
      READ(2450,*)NMIG,ALGMIGRATION_DEBUG
      IF(ALGMIGRATION_DEBUG==1)THEN
          OPEN(2451,file='algae_migration_debug.csv',status='unknown')
          WRITE(2451,*)'Method,K,I,JA,NITWQ,JDAY,Asettle,Dens2'
      ENDIF
      
      READ(2450,*)
      ALLOCATE(MIGRATE_GROUP(NMIG),MIGRATE_MODEL(NMIG),AMP(NMIG),PHASE(NMIG),C_COEFF_EXT(NMIG),RAD(NMIG),MIND(NMIG),MAXD(NMIG),DENSI(NMIG),DENBI(NMIG),T_DEC(NMIG),TS_DEC(NMIG),C_DENINC(NMIG),C_DENDEC(NMIG),&
          DEPTH_LIM_ONOFF(NMIG),DEPTH_LIM(NMIG),LOSS_FRAC(NMIG),I_C(NMIG),C_DENINC_1(NMIG),C_DENINC_2(NMIG),C_DENDEC_1(NMIG),C_DENDEC_2(NMIG),DENP_MINS(NMIG),DENP_MINB(NMIG),DEN_COR(MIGI),&
          NMINT(NMIG),MIGON(NMIG,12),MIGOFF(NMIG,12),DEPTH_CALC_ONOFF(NMIG),EXP_DEPTH(NMIG))
      DO I=1,NMIG
        READ(2450,*)MIGRATE_GROUP(I)
        READ(2450,*)
        READ(2450,*)NMINT(I)
        READ(2450,*)
        READ(2450,*)MIGON(I,1:NMINT(I))
        READ(2450,*)
        READ(2450,*)MIGOFF(I,1:NMINT(I))
        READ(2450,*)
        READ(2450,*)MIGRATE_MODEL(I)
        READ(2450,*)
        IF(MIGRATE_MODEL(I)==1 .OR. MIGRATE_MODEL(I)==2)THEN
            READ(2450,*)AMP(I),PHASE(I),C_COEFF_EXT(I),DEPTH_CALC_ONOFF(I),EXP_DEPTH(I),DEPTH_LIM_ONOFF(I),DEPTH_LIM(I),LOSS_FRAC(I)
            READ(2450,*)
            READ(2450,*)
            READ(2450,*)
            READ(2450,*)
            READ(2450,*)
        ELSEIF(MIGRATE_MODEL(I)==3)THEN
            READ(2450,*)
            READ(2450,*)
            READ(2450,*)RAD(I),MIND(I),MAXD(I),DENSI(I),DENBI(I),T_DEC(I),TS_DEC(I),C_DENINC(I),C_DENDEC(I),DEPTH_LIM_ONOFF(I),DEPTH_LIM(I),LOSS_FRAC(I)
            READ(2450,*)
            READ(2450,*)
            READ(2450,*)
        ELSE
            READ(2450,*)
            READ(2450,*)            
            READ(2450,*)
            READ(2450,*)            
            READ(2450,*)RAD(I),MIND(I),MAXD(I),DENSI(I),DENBI(I),I_C(I),C_DENINC_1(I),C_DENINC_2(I),C_DENDEC_1(I),C_DENDEC_2(I),DENP_MINS(I),DENP_MINB(I),DEN_COR(I),DEPTH_LIM_ONOFF(I),&
                DEPTH_LIM(I),LOSS_FRAC(I)
            READ(2450,*)            
        ENDIF
        ALGAE_SETTLING(MIGRATE_GROUP(I)) = .TRUE.  ! CO 6/5/2019
      ENDDO
      ALLOCATE(DEN_AVG(KMX,IMX,NMIG),DEN(KMX,IMX,MAXVAL(TS_DEC),NMIG),DEN1(KMX,IMX,NMIG),DEN2(KMX,IMX,NMIG),DENP(KMX,IMX,NMIG),DENP_MIN(KMX),ASETTLE(KMX,IMX,NAL),ALLIM_OLD(KMX,IMX,NMIG),&
          TWQ(MAXVAL(TS_DEC),NMIG),LOLD(IMX,NMIG))
100   CLOSE(2450)
  ENDIF

!ALGAE_TOXIN=.FALSE.
!INQUIRE(FILE='w2_Algae_Toxin.csv',EXIST=ALGAE_TOXIN_FILE)    
  IF(ALGAE_TOXIN)THEN
      OPEN(2450,FILE='w2_Algae_Toxin.csv',STATUS='OLD')
      READ(2450,*)  ! SKIP HEADER
      READ(2450,*)ATOX,ATOX_DEBUG     ! '(A2)'
      IF(ATOX == 'ON')THEN
      ALLOCATE(CTP(NUMATOXINS,NAL),CTB(NUMATOXINS,NAL),IN_TOXIN(KMX,IMX,NUMATOXINS))
      
      READ(2450,*)
      READ(2450,*)(CTP(1,JA),JA=1,NAL)
      READ(2450,*)(CTB(1,JA),JA=1,NAL)
      READ(2450,*) CTREL(1)      !(CTL(1,JA),JA=1,NAL)
      !READ(2450,*) CTDI(1)     !(CTDI(1,JA),JA=1,NAL)
      !READ(2450,*) CTA(1)      !(CTA(1,JA),JA=1,NAL)
      READ(2450,*) CTD(1)     !(CTDE(1,JA),JA=1,NAL)
      READ(2450,*)
      READ(2450,*)(CTP(2,JA),JA=1,NAL)
      READ(2450,*)(CTB(2,JA),JA=1,NAL)
      READ(2450,*) CTREL(2)     !(CTL(2,JA),JA=1,NAL)
      !READ(2450,*) CTDI(2)    !(CTDI(2,JA),JA=1,NAL)
      !READ(2450,*) CTA(2)     !(CTA(2,JA),JA=1,NAL)
      READ(2450,*) CTD(2)    !(CTDE(2,JA),JA=1,NAL)
      READ(2450,*)
      READ(2450,*)(CTP(3,JA),JA=1,NAL)
      READ(2450,*)(CTB(3,JA),JA=1,NAL)
      READ(2450,*) CTREL(3)     !(CTL(3,JA),JA=1,NAL)
      !READ(2450,*) CTDI(3)    !(CTDI(3,JA),JA=1,NAL)
      !READ(2450,*) CTA(3)     !(CTA(3,JA),JA=1,NAL)
      READ(2450,*) CTD(3)    !(CTDE(3,JA),JA=1,NAL)
      READ(2450,*)
      READ(2450,*)(CTP(4,JA),JA=1,NAL)
      READ(2450,*)(CTB(4,JA),JA=1,NAL)
      READ(2450,*) CTREL(4)     !(CTL(4,JA),JA=1,NAL)
      !READ(2450,*) CTDI(4)    !(CTDI(4,JA),JA=1,NAL)
      !READ(2450,*) CTA(4)     !(CTA(4,JA),JA=1,NAL)
      READ(2450,*) CTD(4)    !(CTDE(4,JA),JA=1,NAL)
      
      ! 
      DO JJ=1,NUMATOXINS
          CTREL(JJ)=CTREL(JJ)/86400.
          CTD(JJ)=CTD(JJ)/86400.
      ENDDO
      
      ELSE
      ALGAE_TOXIN=.FALSE.   
      ENDIF
  ENDIF
      
!!**** Cyanotoxin Constants
!CTP(J) =fraction of algae concentration producing toxin
!CTB(J) =ratio of intracellular toxin to dry weight OM
!
!!**** Cyanotoxin Rates
!NOT USED CTLR(J) = CTL(J) !leakage from the cell
!CTR release of toxin from inside the cell
!NOT USED CTDIR(J) = CTDI(J) !intracellular decay
!CTD !extracellular decay

      
RETURN

!***********************************************************************************************************************************
!**                                      T E M P E R A T U R E  R A T E  M U L T I P L I E R S                                    **
!***********************************************************************************************************************************

ENTRY TEMPERATURE_RATES
  DO I=IU,ID
    DO K=KT,KB(I)
      LAM1        = FR(T1(K,I),NH4T1(JW),NH4T2(JW),NH4K1(JW),NH4K2(JW))
      NH4TRM(K,I) = LAM1/(1.0+LAM1-NH4K1(JW))
      LAM1        = FR(T1(K,I),NO3T1(JW),NO3T2(JW),NO3K1(JW),NO3K2(JW))
      NO3TRM(K,I) = LAM1/(1.0+LAM1-NO3K1(JW))
      LAM1        = FR(T1(K,I),OMT1(JW),OMT2(JW),OMK1(JW),OMK2(JW))
      OMTRM(K,I)  = LAM1/(1.0+LAM1-OMK1(JW))
      LAM1        = FR(T1(K,I),SODT1(JW),SODT2(JW),SODK1(JW),SODK2(JW))
      SODTRM(K,I) = LAM1/(1.0+LAM1-SODK1(JW))
      DO JA=1,NAL
        IF(ALG_CALC(JA))THEN
        LAM1          = FR(T1(K,I),AT1(JA),AT2(JA),AK1(JA),AK2(JA))
        LAM2          = FF(T1(K,I),AT3(JA),AT4(JA),AK3(JA),AK4(JA))
        ATRMR(K,I,JA) = LAM1/(1.0+LAM1-AK1(JA))
        ATRMF(K,I,JA) = LAM2/(1.0+LAM2-AK4(JA))
        ATRM(K,I,JA)  = ATRMR(K,I,JA)*ATRMF(K,I,JA)
        ENDIF
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))then
        LAM1          = FR(T1(K,I),ET1(JE),ET2(JE),EK1(JE),EK2(JE))
        LAM2          = FF(T1(K,I),ET3(JE),ET4(JE),EK3(JE),EK4(JE))
        ETRMR(K,I,JE) = LAM1/(1.0+LAM1-EK1(JE))
        ETRMF(K,I,JE) = LAM2/(1.0+LAM2-EK4(JE))
        ETRM(K,I,JE)  = ETRMR(K,I,JE)*ETRMF(K,I,JE)
        endif
      END DO
      DO M=1,NMC
      IF(MACROPHYTE_CALC(JW,M))THEN
        LAM1    = FR(T1(K,I),MT1(M),MT2(M),MK1(M),MK2(M))
        LAM2    = FF(T1(K,I),MT3(M),MT4(M),MK3(M),MK4(M))
        MACTRMR(K,I,M) = LAM1/(1.0+LAM1-MK1(M))
        MACTRMF(K,I,M) = LAM2/(1.0+LAM2-MK4(M))
        MACTRM(K,I,M)  = MACTRMR(K,I,M)*MACTRMF(K,I,M)
      endif
      end do
      IF(ZOOPLANKTON_CALC)THEN
	    DO JZ = 1, NZP
          LAM1       = FR(T1(K,I),ZT1(JZ),ZT2(JZ),ZK1(JZ),ZK2(JZ))
          LAM2       = FF(T1(K,I),ZT3(JZ),ZT4(JZ),ZK3(JZ),ZK4(JZ))
          ZOORMR(K,I,JZ)= LAM1/(1.+LAM1-ZK1(JZ))
          ZOORMF(K,I,JZ)= LAM2/(1.+LAM2-ZK4(JZ))
          ZOORM(K,I,JZ) = ZOORMR(K,I,JZ)*ZOORMF(K,I,JZ)
        END DO
	  end if
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                 K I N E T I C   R A T E S                                                     **
!***********************************************************************************************************************************

ENTRY KINETIC_RATES

! Decay rates
!!$OMP PARALLEL DO
  DO I=IU,ID
      
      ! Gas Transfer
      IF(GAS_TRANSFER_UPDATE)THEN
          CALL GAS_TRANSFER
      ENDIF

    DO K=KT,KB(I)
      DO1(K,I)          = O2(K,I)/(O2(K,I)+KDO)                  
      DO2(K,I)          = 1.0 - DO1(K,I)                         !O2(K,I)/(O2(K,I)+KDO)
      DO3(K,I)          = (1.0+SIGN(1.0,O2(K,I)-1.E-10)) *0.5
      SODTRMDDO3        =   SODTRM(K,I)*SDKV(K,I)*DO3(K,I)
      SEDD(K,I)         =   SODTRMDDO3*SED(K,I)    !SODTRM(K,I) *SDKV(K,I)   *SED(K,I) *DO3(K,I)   !CB 10/22/06
      SEDDP(K,I)        =   SODTRMDDO3*SEDP(K,I)   !SODTRM(K,I) *SDKV(K,I)   *SEDP(K,I) *DO3(K,I)
      SEDDN(K,I)        =   SODTRMDDO3*SEDN(K,I)   !SODTRM(K,I) *SDKV(K,I)   *SEDN(K,I) *DO3(K,I)
      SEDDC(K,I)        =   SODTRMDDO3*SEDC(K,I)   !SODTRM(K,I) *SDKV(K,I)   *SEDC(K,I) *DO3(K,I)

      IF(STANDING_BIOMASS_DECAY)THEN
      SEDD1(K,I)         =   SODTRM(K,I) *SDK1(jw)   *SED1(K,I) *DO3(K,I)
      SEDD2(K,I)         =   SODTRM(K,I) *SDK2(jw)   *SED2(K,I) *DO3(K,I)
      ENDIF

      SEDBR(K,I)         =  SEDB(JW)    *SED(K,I)                           !CB 11/30/06
      SEDBRP(K,I)        =  SEDB(JW)    *SEDP(K,I)                          !CB 11/30/06
      SEDBRN(K,I)        =  SEDB(JW)    *SEDN(K,I)                          !CB 11/30/06
      SEDBRC(K,I)        =  SEDB(JW)    *SEDC(K,I)                          !CB 11/30/06
      NH4D(K,I)         =  NH4TRM(K,I) *NH4DK(JW) *NH4(K,I) *DO1(K,I)
      
      IF(CDWBC(NH3_DER,JW)=='      ON')THEN              ! if nh3 is ON as a derived variable
          KH_NH3= 0.00000001270002*T2(K,I)*T2(K,I) + 0.00000016769649*T2(K,I) + 0.000004794059   ! Henry's Law constant for ammonia, atm/m3/mole
          K_NH3_NH4=10**(-0.09018-2729.92/(T2(K,I)+273.15))                                       ! Equilibrium constant [-]
          F_NH3(K,I)=1./(1.+10**(-PH(K,I))/K_NH3_NH4)     ! FRACTION NH3
          NH4D(K,I)    =   NH4D(K,I)*(1.-F_NH3(K,I))     ! NH4 FRACTION FOR NITRIFICATION
          IF(K==KT)THEN
              K_NH3=KG_H2O_CONSTANT(JW)*WIND2(I)*KH_NH3/(RGAS*(T2(K,I)+273.15))*5.86806E-6       !0.5*1.014/86400.=5.86806E-6         ! Gas transfer coefficient for NH3: based on water vapor, gas film controls, 1.014 is ratio of (MW H2o/MW NH3) ^0.25, 0.5 is wind reduction based on wind measuring height, units: m/s
              NH3GAS(K,I)=K_NH3*NH4(K,I)*F_NH3(K,I)*BI(KT,I)/BH2(KT,I)
         ENDIF
      ENDIF
      
      NO3D(K,I)         =  NO3TRM(K,I) *NO3DK(JW) *NO3(K,I) *DO2(K,I)
      IF(ORGC_CALC)THEN
        OMTRMDO3        =  OMTRM(K,I)*DO3(K,I)/ORGC(JW)
        LDOMD(K,I)      =  OMTRMDO3*LDOMCDK(JW)*LDOMC(K,I) !OMTRM(K,I)  *LDOMCDK(JW)*LDOMC(K,I)/ORGC(JW)*DO3(K,I)
        RDOMD(K,I)      =  OMTRMDO3*RDOMCDK(JW)*RDOMC(K,I) !OMTRM(K,I)  *RDOMCDK(JW)*RDOMC(K,I)/ORGC(JW)*DO3(K,I)
        LPOMD(K,I)      =  OMTRMDO3*LPOMCDK(JW)*LPOMC(K,I) !OMTRM(K,I)  *LPOMCDK(JW)*LPOMC(K,I)/ORGC(JW)*DO3(K,I)
        RPOMD(K,I)      =  OMTRMDO3*RPOMCDK(JW)*RPOMC(K,I) !OMTRM(K,I)  *RPOMCDK(JW)*RPOMC(K,I)/ORGC(JW)*DO3(K,I)
        LRDOMD(K,I)     =  OMTRMDO3*LRDDK(JW) *LDOMC(K,I)  !OMTRM(K,I)  *LRDDK(JW) *LDOMC(K,I)/ORGC(JW)*DO3(K,I)
        LRPOMD(K,I)     =  OMTRMDO3*LRPDK(JW) *LPOMC(K,I)  !OMTRM(K,I)  *LRPDK(JW) *LPOMC(K,I)/ORGC(JW)*DO3(K,I)
        LPOMHD(K,I)     =  OMTRMDO3*LPOMHK(JW)*LPOMC(K,I)  !OMTRM(K,I)  *LPOMHK(JW)*LPOMC(K,I)/ORGC(JW)*DO3(K,I)     
        RPOMHD(K,I)     =  OMTRMDO3*RPOMHK(JW)*RPOMC(K,I)  !OMTRM(K,I)  *RPOMHK(JW)*RPOMC(K,I)/ORGC(JW)*DO3(K,I)
      ELSE
        OMTRMDO3        =  OMTRM(K,I)*DO3(K,I)
        LDOMD(K,I)      =  OMTRMDO3*LDOMDK(JW)*LDOM(K,I)   !OMTRM(K,I)  *LDOMDK(JW)*LDOM(K,I)*DO3(K,I)
        RDOMD(K,I)      =  OMTRMDO3*RDOMDK(JW)*RDOM(K,I)   !OMTRM(K,I)  *RDOMDK(JW)*RDOM(K,I)*DO3(K,I)
        LPOMD(K,I)      =  OMTRMDO3*LPOMDK(JW)*LPOM(K,I)   !OMTRM(K,I)  *LPOMDK(JW)*LPOM(K,I)*DO3(K,I)
        RPOMD(K,I)      =  OMTRMDO3*RPOMDK(JW)*RPOM(K,I)   !OMTRM(K,I)  *RPOMDK(JW)*RPOM(K,I)*DO3(K,I)
        LRDOMD(K,I)     =  OMTRMDO3*LRDDK(JW) *LDOM(K,I)   !OMTRM(K,I)  *LRDDK(JW) *LDOM(K,I)*DO3(K,I)
        LRPOMD(K,I)     =  OMTRMDO3*LRPDK(JW) *LPOM(K,I)   !OMTRM(K,I)  *LRPDK(JW) *LPOM(K,I)*DO3(K,I)
        LPOMHD(K,I)     =  OMTRMDO3*LPOMHK(JW)*LPOM(K,I)   !OMTRM(K,I)  *LPOMHK(JW)*LPOM(K,I)*DO3(K,I)           
        RPOMHD(K,I)     =  OMTRMDO3*RPOMHK(JW)*RPOM(K,I)   !OMTRM(K,I)  *RPOMHK(JW)*RPOM(K,I)*DO3(K,I)
      END IF
      !
      IF(CAC(NLDOMP) == '      ON')THEN
        LDOMPD(K,I)     =  OMTRM(K,I)  *LDOMPDK(JW) *LDOMP(K,I)*DO3(K,I)
        LRDOMPD(K,I)    =  OMTRM(K,I)  *LRDOMPDK(JW)*LDOMP(K,I)*DO3(K,I) 
      ELSE
        LDOMPD(K,I)     =  ORGP(JW)    *LDOMD(K,I)
        LRDOMPD(K,I)    =  ORGP(JW)    *LRDOMD(K,I)
      END IF
      IF(CAC(NRDOMP) == '      ON')THEN
        RDOMPD(K,I)     =  OMTRM(K,I)  *RDOMPDK(JW) *RDOMP(K,I)*DO3(K,I)
      ELSE
        RDOMPD(K,I)     =  ORGP(JW)    *RDOMD(K,I)
      END IF
      IF(CAC(NLPOMP) == '      ON')THEN
        LPOMPD(K,I)     =  OMTRM(K,I)  *LPOMPDK(JW) *LPOMP(K,I)*DO3(K,I)
        LRPOMPD(K,I)    =  OMTRM(K,I)  *LRPOMPDK(JW)*LPOMP(K,I)*DO3(K,I)
        LPOMPHD(K,I)    =  OMTRM(K,I)  *LPOMHK(JW)  *LPOMP(K,I)
      ELSE
        LPOMPD(K,I)     =  ORGP(JW)    *LPOMD(K,I)
        LRPOMPD(K,I)    =  ORGP(JW)    *LRPOMD(K,I)
        LPOMPHD(K,I)    =  ORGP(JW)    *LPOMHD(K,I)
      END IF
      IF(CAC(NRPOMP) == '      ON')THEN
        RPOMPD(K,I)     =  OMTRM(K,I)  *RPOMPDK(JW) *RPOMP(K,I)*DO3(K,I)
        RPOMPHD(K,I)    =  OMTRM(K,I)  *RPOMHK(JW)  *RPOMP(K,I)
      ELSE
        RPOMPD(K,I)     =  ORGP(JW)    *RPOMD(K,I)
        RPOMPHD(K,I)    =  ORGP(JW)    *RPOMHD(K,I)
      END IF
      IF(CAC(NLDOMN) == '      ON')THEN
        LDOMND(K,I)     =  OMTRM(K,I)  *LDOMNDK(JW) *LDOMN(K,I)*DO3(K,I)
        LRDOMND(K,I)    =  OMTRM(K,I)  *LRDOMNDK(JW)*LDOMN(K,I)*DO3(K,I)
      ELSE
        LDOMND(K,I)     =  ORGN(JW)    *LDOMD(K,I)
        LRDOMND(K,I)    =  ORGN(JW)    *LRDOMD(K,I)
      END IF
      IF(CAC(NRDOMN) == '      ON')THEN
        RDOMND(K,I)     =  OMTRM(K,I)  *RDOMNDK(JW) *RDOMN(K,I)*DO3(K,I)
      ELSE
        RDOMND(K,I)     =  ORGN(JW)    *RDOMD(K,I)
      END IF
      IF(CAC(NLPOMN) == '      ON')THEN
        LPOMND(K,I)     =  OMTRM(K,I)  *LPOMNDK(JW) *LPOMN(K,I)*DO3(K,I)
        LRPOMND(K,I)    =  OMTRM(K,I)  *LRPOMNDK(JW)*LPOMN(K,I)*DO3(K,I)
        LPOMNHD(K,I)    =  OMTRM(K,I)  *LPOMHK(JW)  *LPOMN(K,I)
      ELSE
        LPOMND(K,I)     =  ORGN(JW)    *LPOMD(K,I)
        LRPOMND(K,I)    =  ORGN(JW)    *LRPOMD(K,I)
        LPOMNHD(K,I)    =  ORGN(JW)    *LPOMHD(K,I)
      END IF
      IF(CAC(NRPOMN) == '      ON')THEN
        RPOMND(K,I)     =  OMTRM(K,I)  *RPOMNDK(JW) *RPOMN(K,I)*DO3(K,I)
        RPOMNHD(K,I)    =  OMTRM(K,I)  *RPOMHK(JW)  *RPOMN(K,I)
      ELSE
        RPOMND(K,I)     =  ORGN(JW)    *RPOMD(K,I)
        RPOMNHD(K,I)    =  ORGN(JW)    *RPOMHD(K,I)
      END IF
      !
      IF(ORGC_CALC)THEN
        OMTRMDO3      =  OMTRM(K,I)*DO3(K,I)
        LDOMCD(K,I)   =  OMTRMDO3*LDOMCDK(JW) *LDOMC(K,I)  !OMTRM(K,I)  *LDOMCDK(JW) *LDOMC(K,I)*DO3(K,I)
        LRDOMCD(K,I)  =  OMTRMDO3*LRDOMCDK(JW)*LDOMC(K,I)  !OMTRM(K,I)  *LRDOMCDK(JW)*LDOMC(K,I)*DO3(K,I)
        RDOMCD(K,I)   =  OMTRMDO3*RDOMCDK(JW) *RDOMC(K,I)  !OMTRM(K,I)  *RDOMCDK(JW) *RDOMC(K,I)*DO3(K,I)
        LPOMCD(K,I)   =  OMTRMDO3*LPOMCDK(JW) *LPOMC(K,I)  !OMTRM(K,I)  *LPOMCDK(JW) *LPOMC(K,I)*DO3(K,I)
        LRPOMCD(K,I)  =  OMTRMDO3*LRPOMCDK(JW)*LPOMC(K,I)  !OMTRM(K,I)  *LRPOMCDK(JW)*LPOMC(K,I)*DO3(K,I)
        LPOMCHD(K,I)  =  OMTRMDO3*LPOMHK(JW)  *LPOMC(K,I)  !OMTRM(K,I)  *LPOMHK(JW)  *LPOMC(K,I)*DO3(K,I)
        RPOMCD(K,I)   =  OMTRMDO3*RPOMCDK(JW) *RPOMC(K,I)  !OMTRM(K,I)  *RPOMCDK(JW) *RPOMC(K,I)*DO3(K,I)
        RPOMCHD(K,I)  =  OMTRMDO3*RPOMHK(JW)  *RPOMC(K,I)  !OMTRM(K,I)  *RPOMHK(JW)  *RPOMC(K,I)*DO3(K,I)
      ELSE
        LDOMCD(K,I)   =  ORGC(JW)    *LDOMD(K,I)
        LRDOMCD(K,I)  =  ORGC(JW)    *LRDOMD(K,I)
        RDOMCD(K,I)   =  ORGC(JW)    *RDOMD(K,I)
        LPOMCD(K,I)   =  ORGC(JW)    *LPOMD(K,I)
        LRPOMCD(K,I)  =  ORGC(JW)    *LRPOMD(K,I)
        LPOMCHD(K,I)  =  ORGC(JW)    *LPOMHD(K,I)
        RPOMCD(K,I)   =  ORGC(JW)    *RPOMD(K,I)
        RPOMCHD(K,I)  =  ORGC(JW)    *RPOMHD(K,I)
      END IF
      !
      CBODD(K,I,1:NBOD) =  KBOD(1:NBOD)*TBOD(1:NBOD)**(T1(K,I)-20.0)*DO3(K,I)
        IF(K == KB(I))THEN     ! SW 4/18/07
	  SODD(K,I)         =  SOD(I)/BH2(K,I)*SODTRM(K,I)*BI(K,I)
	    ELSE
      SODD(K,I)         =  SOD(I)/BH2(K,I)*SODTRM(K,I)*(BI(K,I)-BI(K+1,I))
	    ENDIF

! Inorganic suspended solids settling rates - P adsorption onto SS and Fe
    IF(PARTP(JW) > 0.0)THEN    ! SW 3/2019
      FPSS(K,I) = PARTP(JW)         /(PARTP(JW)*TISS(K,I)+PARTP(JW)*FE(K,I)*DO1(K,I)+1.0)
        FPFE(K,I) = PARTP(JW)*FE(K,I)*DO1(K,I) /(PARTP(JW)*TISS(K,I)+PARTP(JW)*FE(K,I)*DO1(K,I)+1.0)   !8/2020 corrected
    ENDIF
    
      IF(K.NE.KT)THEN
          SSSI(K,I) = SSSO(K-1,I)*BI(K,I)/BI(K-1,I)   ! SR 3/2019
      ELSE
          SSSI(K,I)=0.0                               ! SW 3/2019
      ENDIF
      
      TOTSS0    = 0.0
      DO JS=1,NSS
        TOTSS0 = TOTSS0+SSS(JS)*FPSS(K,I)*SS(K,I,JS)
      END DO
      SSSO(K,I) =  TOTSS0 + FeSetVel(JW)*FPFE(K,I)   !8/2020 corrected                                                 !(TOTSS0+FES(JW)*FPFE(K,I))*BI(K,I)/BH2(K,I)*DO1(K,I)                ! SW 11/7/07
      FPSS(K,I) =  FPSS(K,I)*TISS(K,I)

    ! OM stoichiometry
    !8/2020 add 12 temp variables
        IF(CAC(NLDOMP) == '      ON')THEN
      LDOP(K,I) = LDOMP(K,I)
    ELSE IF(ORGC_CALC)THEN
      LDOP(K,I) = LDOMC(K,I)/ORGC(JW) * ORGP(JW)
    ELSE 
      LDOP(K,I) = LDOM(K,I) * ORGP(JW)  
        END IF
        IF(CAC(NRDOMP) == '      ON')THEN
      RDOP(K,I) = RDOMP(K,I)
    ELSE IF(ORGC_CALC)THEN
      RDOP(K,I) = RDOMC(K,I)/ORGC(JW) * ORGP(JW)  
          ELSE
      RDOP(K,I) = RDOM(K,I) * ORGP(JW)
        END IF
        IF(CAC(NLPOMP) == '      ON')THEN
      LPOP(K,I) = LPOMP(K,I)
    ELSE IF(ORGC_CALC)THEN
      LPOP(K,I) = LPOMC(K,I)/ORGC(JW) * ORGP(JW)  
        ELSE
      LPOP(K,I) = LPOM(K,I) * ORGP(JW)
        END IF
        IF(CAC(NRPOMP) == '      ON')THEN
      RPOP(K,I) = RPOMP(K,I)
    ELSE IF(ORGC_CALC)THEN
      RPOP(K,I) = RPOMC(K,I)/ORGC(JW) * ORGP(JW)  
          ELSE
      RPOP(K,I) = RPOM(K,I) * ORGP(JW)
        END IF
        IF(CAC(NLDOMN) == '      ON')THEN
      LDON(K,I) = LDOMN(K,I)
    ELSE IF(ORGC_CALC)THEN
      LDON(K,I) = LDOMC(K,I)/ORGC(JW) * ORGN(JW)  
          ELSE
      LDON(K,I) = LDOM(K,I) * ORGN(JW)
        END IF
        IF(CAC(NRDOMN) == '      ON')THEN
      RDON(K,I) = RDOMN(K,I)
    ELSE IF(ORGC_CALC)THEN
      RDON(K,I) = RDOMC(K,I)/ORGC(JW) * ORGN(JW)  
        ELSE
      RDON(K,I) = RDOM(K,I) * ORGN(JW)
        END IF
        IF(CAC(NLPOMN) == '      ON')THEN
      LPON(K,I) = LPOMN(K,I)
    ELSE IF(ORGC_CALC)THEN
      LPON(K,I) = LPOMC(K,I)/ORGC(JW) * ORGN(JW)  
          ELSE
      LPON(K,I) = LPOM(K,I) * ORGN(JW)
    END IF

    IF(CAC(NRPOMN) == '      ON')THEN
      RPON(K,I) = RPOMN(K,I)
    ELSE IF(ORGC_CALC)THEN
      RPON(K,I) = RPOMC(K,I)/ORGC(JW) * ORGN(JW)  
        ELSE
      RPON(K,I) = RPOM(K,I) * ORGN(JW)
        END IF
    !
    IF(ORGC_CALC)THEN  
      LDOC(K,I) = LDOMC(K,I)
      RDOC(K,I) = RDOMC(K,I)
      LPOC(K,I) = LPOMC(K,I)
      RPOC(K,I) = RPOMC(K,I)
          ELSE
      LDOC(K,I) = LDOM(K,I) * ORGC(JW)
      RDOC(K,I) = RDOM(K,I) * ORGC(JW)
      LPOC(K,I) = LPOM(K,I) * ORGC(JW)
      RPOC(K,I) = RPOM(K,I) * ORGC(JW)  
          END IF        

! Light Extinction Coefficient
      IF (.NOT. READ_EXTINCTION(JW)) THEN
      ALGEX = 0.0; SSEXT = 0.0; ZOOEXT = 0.0                                                     ! SW 11/8/07
        DO JA=1,NAL
          IF(ALG_CALC(JA))ALGEX = ALGEX+EXA(JA)*ALG(K,I,JA)
        END DO
        DO JS=1,NSS
          SSEXT = SSEXT+EXSS(JW)*SS(K,I,JS)
        END DO
 !       TOTMAC=0.0                                                                ! SW 4/20/11 Delete this section?
 !       DO M=1,NMC
 !         IF(MACROPHYTE_CALC(JW,M))THEN
 !           JT=KTI(I)
 !           JE=KB(I)
 !           DO JJ=JT,JE
 !             TOTMAC = EXM(M)*MACRM(JJ,K,I,M)+TOTMAC
 !           END DO
 !         END IF
 !       END DO
 !       MACEXT=TOTMAC/(BH2(K,I)*DLX(I))

	    IF(ZOOPLANKTON_CALC)THEN
	        DO JZ = 1,NZP
	        ZOOEXT = ZOOEXT + ZOO(K,I,JZ)*EXZ(JZ)
	        END DO
	    ENDIF
      !
      GAMMA(K,I) = EXH2O(JW)+SSEXT+ALGEX+ZOOEXT
      IF(ORGC_CALC)THEN
        GAMMA(K,I) = GAMMA(K,I)+EXOM(JW)*(LPOMC(K,I)+RPOMC(K,I))/ORGC(JW)
      ELSE
        GAMMA(K,I) = GAMMA(K,I)+EXOM(JW)*(LPOM(K,I)+RPOM(K,I))
      END IF
		
	    IF(NMC>0)THEN    ! cb 4/20/11
	      MACEXT1=0.0    ! cb 4/20/11
          IF(KTICOL(I))THEN
            JT=KTI(I)
          ELSE
            JT=KTI(I)+1
          END IF
          JE=KB(I)
          DO JJ=JT,JE
            TOTMAC=0.0
            DO M=1,NMC
              IF(MACROPHYTE_CALC(JW,M))THEN
                TOTMAC = EXM(M)*MACRM(JJ,K,I,M)+TOTMAC
              END IF
            END DO
            IF(CW(JJ,I).GT.0.0)THEN
              MACEXT=TOTMAC/(CW(JJ,I)*DLX(I)*H2(K,I))
            ELSE
              MACEXT=0.0
            END IF
			GAMMAJ(JJ,K,I) = GAMMA(K,I)+MACEXT       ! SW 4/20/11
            MACEXT1 = MACEXT*CW(JJ,I)+MACEXT1    ! cb 4/20/11
          END DO
          GAMMA(K,I) = GAMMA(K,I) + MACEXT1/B(JT,I)                                      ! SW 4/21/11
        end if
      ELSE
        GAMMA(K,I) = EXH2O(JW)
      END IF

! Zooplankton Rates
   IF(ZOOPLANKTON_CALC)THEN
      DO JZ=1,NZP
        IF(ORGC_CALC)THEN 
          TGRAZE(K,I,JZ) = PREFP(JZ)*LPOMC(K,I)/ORGC(JW)
        ELSE
          TGRAZE(K,I,JZ) = PREFP(JZ)*LPOM(K,I)
        END IF
        DO JJZ = 1, NZP
          TGRAZE(K,I,JZ) = TGRAZE(K,I,JZ) + PREFZ(JJZ,JZ)*ZOO(K,I,JJZ)          !CB 5/17/2007
      END DO
        DO JA=1,NAL
          IF(ALG_CALC(JA))TGRAZE(K,I,JZ)=PREFA(JA,JZ)*ALG(K,I,JA)+TGRAZE(K,I,JZ)
        END DO
        ZMINFAC  = (1.0+SIGN(1.0,ZOO(K,I,JZ)-ZOOMIN(JZ)))*0.5
        ZRT(K,I,JZ) =  ZOORMR(K,I,JZ)*ZR(JZ)*ZMINFAC*DO3(K,I)
        IF (TGRAZE(K,I,JZ) <= 0.0 .OR. O2(K,I) < 2.0) THEN
          ZMU(K,I,JZ)       = 0.0
          AGZ(K,I,1:NAL,JZ) = 0.0
		  ZGZ(K,I,JZ,:) = 0.0
          IF (O2(K,I) < 2.0) ZMINFAC = 2*ZMINFAC
        ELSE
          ZMU(K,I,JZ) = MAX(ZOORM(K,I,JZ)*ZG(JZ)*(TGRAZE(K,I,JZ)-ZOOMIN(JZ))/(TGRAZE(K,I,JZ)+ZS2P(JZ)), 0.0)
          DO JA=1,NAL
          IF(ALG_CALC(JA))AGZ(K,I,JA,JZ) = ZMU(K,I,JZ)*ZOO(K,I,JZ)*(ALG(K,I,JA)*PREFA(JA,JZ)/TGRAZE(K,I,JZ))                      !  KV 5/26/2007
          END DO
          DO JJZ = 1,NZP ! OMNIVOROUS ZOOPLANKTON
          ZGZ(K,I,JJZ,JZ)  = ZMU(K,I,JZ)*ZOO(K,I,JZ)*(ZOO(K,I,JJZ)*PREFZ(JJZ,JZ)/TGRAZE(K,I,JZ))         !KV 5/26/2007
          END DO
        END IF
        ZMT(K,I,JZ) = MAX(1.0-ZOORMF(K,I,JZ),0.02)*ZM(JZ)*ZMINFAC
        ! zooplankton settling - adapted from SR 01/12/2004 Hagg Lake Model - in prep for dynamic vertical motion calculation SW 1/28/2019
        IF (ZS(JZ) >= 0.0) THEN                                                                                            
              IF (K == KT) THEN
              ZSR(K,I,JZ) = -ZS(JZ)*ZOO(K,I,JZ)*BI(K,I)/BH2(K,I)
              ELSE
              ZSR(K,I,JZ) =  ZS(JZ)*(ZOO(K-1,I,JZ)-ZOO(K,I,JZ))*BI(K,I)/BH2(K,I)  
              ENDIF
        ELSE                                                                                                           
          IF (K == KT) THEN
              ZSR(K,I,JZ) = -ZS(JZ)*ZOO(K+1,I,JZ)*BI(K+1,I)*DLX(I)/VOL(K,I)   
          ELSEIF(K == KB(I))THEN
              ZSR(K,I,JZ) =  ZS(JZ)*ZOO(K,I,JZ)*BI(K,I)/BH2(K,I) 
          ELSE
              ZSR(K,I,JZ) = -ZS(JZ)*(ZOO(K+1,I,JZ)*BI(K+1,I)/BH2(K,I)-ZOO(K,I,JZ)*BI(K,I)/BH2(K,I))    
          ENDIF                     
        END IF                                                                                                                   
      END DO   ! ZOOP LOOP 
   ENDIF

    END DO ! K LOOP
  END DO   ! I LOOP
!!$OMP END PARALLEL DO

! Algal rates
  if(ALGAE_SETTLING_EXIST .AND. jday>iday) NITWQ=NITWQ+1
     DO JA=1,NAL
      IF(ALG_CALC(JA))THEN
          IF(ALGAE_SETTLING(JA))THEN           ! If there is variable velocity - call algae migration subroutine   ! COMPUTE ALGAE MIGRATION RATE
            DO MIGI=1,NMIG                             ! CO 6/5/2019
                IF(MIGRATE_GROUP(MIGI) == JA)THEN
                    EXIT
                ENDIF
            ENDDO
            ISETTLE=0
            DO KK=1,NMINT(MIGI)
                IF(JDAY >= MIGON(MIGI,KK) .AND. JDAY <= MIGOFF(MIGI,KK))THEN
                    ISETTLE=1
                    EXIT
                ENDIF
            ENDDO
        ENDIF
      do i=iu,id
!**** Limiting factor
      LIGHT = (1.0-BETA(JW))*SRON(JW)*SHADE(I)/ASAT(JA)
      LAM1  =  LIGHT
      LAM2  =  LIGHT
      DO K=KT,KB(I)

!****** Limiting factor
        LAM1           = LAM2
        LAM2           = LAM1*EXP(-GAMMA(K,I)*H2(K,I))
        FDPO4          = 1.0-FPSS(K,I)-FPFE(K,I)
        ALLIM(K,I,JA)  = 2.718282*(EXP(-LAM2)-EXP(-LAM1))/(GAMMA(K,I)*H2(K,I))
        IF (AHSP(JA)  /= 0.0 .and. po4_calc) APLIM(K,I,JA) =  FDPO4*PO4(K,I)/(FDPO4*PO4(K,I)+AHSP(JA)+NONZERO)       ! cb 10/12/11
        IF (AHSN(JA)  /= 0.0 .and. n_calc) ANLIM(K,I,JA) = (NH4(K,I)+NO3(K,I))/(NH4(K,I)+NO3(K,I)+AHSN(JA)+NONZERO)  ! cb 10/12/11
        IF (AHSSI(JA) /= 0.0 .and. DSI_CALC) ASLIM(K,I,JA) =  DSI(K,I)/(DSI(K,I)+AHSSI(JA)+NONZERO)                  ! cb 10/12/11
        LIMIT          = MIN(APLIM(K,I,JA),ANLIM(K,I,JA),ASLIM(K,I,JA),ALLIM(K,I,JA))

!****** Algal rates
        AGR(K,I,JA) =  ATRM(K,I,JA)*AG(JA)*LIMIT
        ARR(K,I,JA) =  ATRM(K,I,JA)*AR(JA)*DO3(K,I)
        AMR(K,I,JA) = (ATRMR(K,I,JA)+1.0-ATRMF(K,I,JA))*AM(JA)
        AER(K,I,JA) =  MIN((1.0-ALLIM(K,I,JA))*AE(JA)*ATRM(K,I,JA),AGR(K,I,JA))
                    
        IF(ALGAE_SETTLING(JA) .AND. ISETTLE==1)THEN
            IF(MIGRATE_MODEL(MIGI) == 1)THEN     !   TIME VARYING VELOCTY
                ASETTLE(K,I,JA) = AMP(MIGI)*(2.*PI/86400.)*COS(2.*PI*JDAY + PHASE(MIGI))
                IF(ALGMIGRATION_DEBUG==1)WRITE(2451,'("Model1:,",3(I3,","),3(F15.5,","))')K,I,JA,JDAY,ASETTLE(K,I,JA)
            ELSEIF(MIGRATE_MODEL(MIGI) == 2)THEN   !  TIME AND SPACE VARYING VELOCITY
                IF(DEPTH_CALC_ONOFF(MIGI)==1)THEN  ! CALCULATE DEPTH LIMIT FOR INCREASING MIGRATION BASED ON LIGHT
                    LIGHT=(1.0-BETA(JW))*SRON(JW)*SHADE(I)
                    LAM3=LIGHT
                    KK=KT
                    DO WHILE(LAM3 > 0.01*LIGHT) 
                        KK=KK+1
                        IF(KK==KB(I))EXIT  
                        LAM3=LAM3*EXP(-GAMMA(KK,I)*H2(KK,I))
                    ENDDO
                    IF(DEPTHM(K,I)<=DEPTHM(KK,I))THEN
                        IF(LIGHT>0.0)THEN
                            ASETTLE(K,I,JA) = AMP(MIGI)*(2.*PI/86400.)*EXP(-C_COEFF_EXT(MIGI)*GAMMA(K,I)*(DEPTHM(KK,I)-DEPTHM(K,I)))*COS(2.*PI*JDAY+PHASE(MIGI)) 
                        ELSE 
                            ASETTLE(K,I,JA) = AMP(MIGI)*(2.*PI/86400.)*COS(2.*PI*JDAY+PHASE(MIGI))
                        ENDIF
                    ELSE
                        ASETTLE(K,I,JA) = AMP(MIGI)*(2.*PI/86400.)*COS(2.*PI*JDAY+PHASE(MIGI))
                    ENDIF
                ELSE    ! SET DEPTH LIMIT FOR INCREASING MIGRATION
                    IF(DEPTHM(K,I)<=EXP_DEPTH(MIGI))THEN
                        IF(LIGHT>0.0)THEN
                            ASETTLE(K,I,JA) = AMP(MIGI)*(2.*PI/86400.)*EXP(-C_COEFF_EXT(MIGI)*GAMMA(K,I)*(EXP_DEPTH(MIGI)-DEPTHM(K,I)))*COS(2.*PI*JDAY+PHASE(MIGI)) 
                        ELSE 
                            ASETTLE(K,I,JA) = AMP(MIGI)*(2.*PI/86400.)*COS(2.*PI*JDAY+PHASE(MIGI))
                        ENDIF
                    ELSE
                        ASETTLE(K,I,JA) = AMP(MIGI)*(2.*PI/86400.)*COS(2.*PI*JDAY+PHASE(MIGI))
                    ENDIF  
                ENDIF
                IF(ALGMIGRATION_DEBUG==1)WRITE(2451,'("Model2:,",3(I3,","),3(F15.5,","))')K,I,JA,JDAY,ASETTLE(K,I,JA)
            ELSEIF(MIGRATE_MODEL(MIGI) == 3)THEN    ! DENSITY CHANGE VELOCITY
                VISCK = DEXP((T2(K,I)+495.691)/(-37.3877)) ! dynamic viscosity of water
                IF (T2(K,I) > 30.0)  VISCK = DEXP((T2(K,I)+782.190)/(-57.7600))
                IF(NITWQ == 1)THEN ! SET INITIAL DENSITY
                    IF(K==KT)THEN
                        DEN(K,I,NITWQ,MIGI) = DENSI(MIGI)
                        IF(I==IU) TWQ(NITWQ,MIGI)=JDAY
                    ELSEIF(K==KB(I))THEN
                        DEN(K,I,NITWQ,MIGI) = DENBI(MIGI)
                    ELSE
                        DEN(K,I,NITWQ,MIGI) = DENSI(MIGI)+(DENBI(MIGI)-DENSI(MIGI))*(1.-EXP(-DEPTHM(K,I)));
                    ENDIF
                    RHO(K,I) = DENSITY(T2(K,I),DMAX1(TDS(K,I),0.0D0),DMAX1(TISS(K,I),0.0D0))
                    ASETTLE(K,I,JA) = 2.*G*(RAD(MIGI)**2)*(DEN(K,I,NITWQ,MIGI)/RHO(K,I)-1.)/(9.*VISCK) ! stoke's settling velocity  
                    ALLIM_OLD(K,I,MIGI) = ALLIM(K,I,JA)      
                    IF(ALGMIGRATION_DEBUG==1)WRITE(2451,'("Model3:,",4(I6,","),4(F15.5,","))')K,I,JA,nitwq,JDAY,ASETTLE(K,I,JA),DEN_avg(K,I,MIGI),den(k,i,min(nitwq,ts_dec(migi)),migi)
                ELSEIF(NITWQ <= TS_DEC(MIGI))THEN  ! STORE ALL DENSITY VALUES UNTIL MAXIMUM NUMBER IS REACHED
                    IF(K==KT .AND. I==IU) TWQ(NITWQ,MIGI)=JDAY
                    if(den(k,i,nitwq-1,migi)<=0)then
                        den(k,i,nitwq-1,migi) = den(k+1,i,nitwq-1,migi)
                        RHO(K,I) = DENSITY(T2(K,I),DMAX1(TDS(K,I),0.0D0),DMAX1(TISS(K,I),0.0D0))
                        allim_old(k,i,migi)=allim(k,i,migi)
                    endif
                    DEN(K,I,NITWQ,MIGI) = (C_DENINC(MIGI)*ALLIM_OLD(K,I,MIGI)-C_DENDEC(MIGI))*(JDAY-TWQ(NITWQ-1,MIGI))*86400. + DEN(K,I,NITWQ-1,MIGI) ! new colony density
                    DEN(K,I,NITWQ,MIGI) = MIN(DEN(K,I,NITWQ,MIGI),MAXD(MIGI)) ! maximum allowable colony density
                    DEN(K,I,NITWQ,MIGI) = MAX(DEN(K,I,NITWQ,MIGI),MIND(MIGI)) ! minimum allowable colony density
                        ALLOCATE(DEN_USE(COUNT(DEN(K,I,1:NITWQ,MIGI)>0)))
                        JJ=1
                        DO LL=1,NITWQ
                            IF (DEN(K,I,LL,MIGI).GT.0) THEN
                                DEN_USE(JJ) = LL
                                JJ = JJ+1
                            END IF
                        ENDDO
                    DEN_AVG(K,I,MIGI) = SUM(DEN(K,I,den_use,MIGI)*EXP(-T_DEC(MIGI)*(JDAY-TWQ(den_use,MIGI))))/SUM(EXP(-T_DEC(MIGI)*(JDAY-TWQ(den_use,MIGI))))
                    ASETTLE(K,I,JA) = 2.*G*(RAD(MIGI)**2)*(DEN_AVG(K,I,MIGI)/RHO(K,I)-1.)/(9.*VISCK) ! stoke's settling velocity
                        DEALLOCATE(DEN_USE)
                    ALLIM_OLD(K,I,MIGI) = ALLIM(K,I,JA)      
                    IF(ALGMIGRATION_DEBUG==1)WRITE(2451,'("Model3:,",4(I6,","),4(F15.5,","))')K,I,JA,nitwq,JDAY,ASETTLE(K,I,JA),DEN_avg(K,I,MIGI),den(k,i,min(nitwq,ts_dec(migi)),migi)
                ELSE  ! STORE SPECIFIED MAXIMUM NUMBER OF PAST DENSITY VALUES 
                    IF(K==KT .AND. I==IU)THEN
                        TWQ(1:TS_DEC(MIGI)-1,MIGI) = TWQ(2:TS_DEC(MIGI),MIGI)
                        TWQ(TS_DEC(MIGI),MIGI) = JDAY
                    ENDIF
                    DEN(K,I,1:TS_DEC(MIGI)-1,MIGI) = DEN(K,I,2:TS_DEC(MIGI),MIGI)                 
                    if(den(k,i,ts_dec(migi)-1,migi)<=0)then
                        den(k,i,ts_dec(migi)-1,migi) = den(k+1,i,ts_dec(migi)-1,migi)
                        RHO(K,I) = DENSITY(T2(K,I),DMAX1(TDS(K,I),0.0D0),DMAX1(TISS(K,I),0.0D0))
                        allim_old(k,i,migi)=allim(k,i,migi)
                    endif
                    DEN(K,I,TS_DEC(MIGI),MIGI) = (C_DENINC(MIGI)*ALLIM_OLD(K,I,MIGI)-C_DENDEC(MIGI))*(JDAY-TWQ(ts_dec(migi)-1,MIGI))*86400. + DEN(K,I,TS_DEC(MIGI)-1,MIGI) ! new colony density
                    DEN(K,I,TS_DEC(MIGI),MIGI) = MIN(DEN(K,I,TS_DEC(MIGI),MIGI),MAXD(MIGI)) ! maximum allowable colony density
                    DEN(K,I,TS_DEC(MIGI),MIGI) = MAX(DEN(K,I,TS_DEC(MIGI),MIGI),MIND(MIGI)) ! minimum allowable colony density
                        ALLOCATE(DEN_USE(COUNT(DEN(K,I,:,MIGI)>0)))
                        JJ=1
                        DO LL=1,TS_DEC(MIGI)
                            IF (DEN(K,I,LL,MIGI).GT.0) THEN
                                DEN_USE(JJ) = LL
                                JJ = JJ+1
                            END IF
                        ENDDO
                    DEN_AVG(K,I,MIGI) = SUM(DEN(K,I,DEN_USE,MIGI)*EXP(-T_DEC(MIGI)*(JDAY-TWQ(DEN_USE,MIGI))))/SUM(EXP(-T_DEC(MIGI)*(JDAY-TWQ(DEN_USE,MIGI)))) ! weighted density with time decay
                    ASETTLE(K,I,JA) = 2.*G*(RAD(MIGI)**2)*(DEN_AVG(K,I,MIGI)/RHO(K,I)-1.)/(9.*VISCK) ! stoke's settling velocity
                        DEALLOCATE(DEN_USE)
                    ALLIM_OLD(K,I,MIGI) = ALLIM(K,I,JA)      
                    IF(ALGMIGRATION_DEBUG==1)WRITE(2451,'("Model3:,",4(I6,","),4(F15.5,","))')K,I,JA,nitwq,JDAY,ASETTLE(K,I,JA),DEN_avg(K,I,MIGI),den(k,i,min(nitwq,ts_dec(migi)),migi)
                ENDIF
            ELSE     ! DENSITY CHANGE VELOCITY (VISSER)
                VISCK = DEXP((T2(K,I)+495.691)/(-37.3877)) ! dynamic viscosity of water
                IF (T2(K,I) > 30.0)  VISCK = DEXP((T2(K,I)+782.190)/(-57.7600))
                DENP_MIN(K) = DENP_MINS(MIGI) + (DEPTHM(K,I)/DEPTHB(KB(I),I))*(DENP_MINB(MIGI)-DENP_MINS(MIGI))
                IF(NITWQ == 1)THEN ! SET INITIAL DENSITY
                    IF(K==KT)THEN
                        DEN1(K,I,MIGI) = DENSI(MIGI)
                    ELSEIF(K==KB(I))THEN
                        DEN1(K,I,MIGI) = DENBI(MIGI)
                    ELSE
                        DEN1(K,I,MIGI) = DENSI(MIGI)+(DENBI(MIGI)-DENSI(MIGI))*(1-EXP(-DEPTHM(K,I)));
                    ENDIF
                    RHO(K,I) = DENSITY(T2(K,I),DMAX1(TDS(K,I),0.0D0),DMAX1(TISS(K,I),0.0D0))
                    ASETTLE(K,I,JA) = 2.*G*(RAD(MIGI)**2)*(DEN1(K,I,MIGI)/RHO(K,I)-1.)/(9.*VISCK) ! stoke's settling velocity
                    if(abs(asettle(k,i,ja))>1000.) asettle(k,1,ja)=0.0
                    DENP(K,I,MIGI) = DENP_MIN(K)
                    IF(ALGMIGRATION_DEBUG==1)WRITE(2451,'("Model4:,",4(I6,","),3(F15.5,","))')K,I,JA,NITWQ,JDAY,ASETTLE(K,I,JA),DEN1(K,I,MIGI)
                    IF(k==kb(i) .and. i==id .AND. jb==nbr .and. jw==nwb) TOLD = JDAY
                    IF(K==KB(I)) LOLD(I,MIGI) = LIGHT*ASAT(JA)
                ELSE
                    if(den1(k,i,migi)<=0)then
                        den1(k,i,migi) = den1(k+1,i,migi)
                        denp(k,i,migi) = denp(k+1,i,migi)
                        RHO(K,I) = DENSITY(T2(K,I),DMAX1(TDS(K,I),0.0D0),DMAX1(TISS(K,I),0.0D0))
                    endif
                    LAM3=LOLD(I,MIGI)
                    LAM4 = LAM3
                    DO KK=KT,K
                        LAM3 = LAM4
                        LAM4 = LAM3*EXP(-GAMMA(KK,I)*H2(KK,I))
                    ENDDO
                    AVG_LIGHT = LAM3*(EXP(-GAMMA(K,I)*H2(K,I))-1.)/(-GAMMA(K,I)*H2(K,I))
                    IF(AVG_LIGHT >= I_C(MIGI))THEN              
                        DEN2(K,I,MIGI) = (C_DENINC_1(MIGI)*AVG_LIGHT*EXP(-AVG_LIGHT/ASAT(JA))+C_DENINC_2(MIGI))*(JDAY-TOLD)*86400. + DEN1(K,I,MIGI) ! new colony density
                        DEN2(K,I,MIGI) = MIN(DEN2(K,I,MIGI),MAXD(MIGI)) ! maximum allowable colony density
                        DEN2(K,I,MIGI) = MAX(DEN2(K,I,MIGI),MIND(MIGI)) ! minimum allowable colony density
                        DENP(K,I,MIGI) = MAX(DEN2(K,I,MIGI),DENP_MIN(K))
                    ELSE
                        DEN2(K,I,MIGI) = (-C_DENDEC_1(MIGI)*(DENP(K,I,MIGI) + DEN_COR(MIGI)) + C_DENDEC_2(MIGI))*(JDAY-TOLD)*86400. + DEN1(K,I,MIGI) ! new colony density
                        DEN2(K,I,MIGI) = MIN(DEN2(K,I,MIGI),MAXD(MIGI)) ! maximum allowable colony density
                        DEN2(K,I,MIGI) = MAX(DEN2(K,I,MIGI),MIND(MIGI)) ! minimum allowable colony density
                        DENP(K,I,MIGI) = MAX(DENP(K,I,MIGI),DENP_MIN(K))
                    ENDIF 
                    ASETTLE(K,I,JA) = 2.*G*(RAD(MIGI)**2)*(DEN2(K,I,MIGI)/RHO(K,I)-1.)/(9.*VISCK) ! stoke's settling velocity
                    if(abs(asettle(k,i,ja))>1000.) asettle(k,1,ja)=0.
                    DEN1(K,I,MIGI) = DEN2(K,I,MIGI)
                    IF(ALGMIGRATION_DEBUG==1)WRITE(2451,'("Model4:,",4(I6,","),3(F15.5,","))')K,I,JA,NITWQ,JDAY,ASETTLE(K,I,JA),DEN2(K,I,MIGI)
                    IF(k==kb(i) .and. i==id .AND. jb==nbr .and. jw==nwb) TOLD = JDAY
                    IF(K==KB(I)) LOLD(I,MIGI) = LIGHT*ASAT(JA)
                ENDIF
            ENDIF
            if(depth_lim_onoff(migi) == 1 .and. depthb(k,i)<depth_lim(migi)) asettle(k,i,ja)=0
        ELSE ! ORIGINAL SETTLING EQUATIONS
        IF (AS(JA) >= 0.0) THEN
          IF(K == KT)THEN
          ASR(K,I,JA) =  AS(JA)*(-ALG(K,I,JA))*BI(K,I)/BH2(K,I)
          ELSE
          ASR(K,I,JA) =  AS(JA)*(ALG(K-1,I,JA)-ALG(K,I,JA))*BI(K,I)/BH2(K,I)
          ENDIF
        ELSE
          IF(K == KB(I))THEN
            ASR(K,I,JA) = -AS(JA)*(-ALG(K,I,JA)  *BI(K,I)/BH2(K,I))                                           !SW 11/8/07
          ELSEIF(K == KT)THEN
            ASR(K,I,JA) = -AS(JA)* ALG(K+1,I,JA)*BI(K+1,I)*DLX(I)/VOL(K,I)                                   !SW 11/8/07
          ELSE
            ASR(K,I,JA) = -AS(JA)*(ALG(K+1,I,JA)*BI(K+1,I)/BH2(K,I)-ALG(K,I,JA)*BI(K,I)/BH2(K,I))             !SP 8/27/07
          END IF
            END IF 
        ENDIF
      ENDDO   ! K LOOP

      IF(ALGAE_SETTLING(JA) .AND. ISETTLE==1)THEN
         do k=kt,kb(i)
            IF(K==KT)THEN
                IF(ASETTLE(K+1,I,JA) >= 0)THEN    ! incoming velocity from cell below
                    AIN = 0                
                ELSE 
                    AIN = -ASETTLE(K+1,I,JA)
                ENDIF             
                IF(ASETTLE(K,I,JA) >= 0)THEN    ! outgoing velocity for surface layer cell
                    AOUT = ASETTLE(K,I,JA)
                ELSE 
                    AOUT = 0
                ENDIF             
                ASR(K,I,JA) =  -AOUT*(ALG(K,I,JA))*BI(K,I)/BH2(K,I) + AIN*ALG(K+1,I,JA)*BI(K+1,I)*DLX(I)/VOL(K,I)                  
            ELSEIF(K==KB(I))THEN
                IF(ASETTLE(K-1,I,JA) >= 0)THEN    ! incoming velocity from cell ABOVE
                    AIN = ASETTLE(K-1,I,JA)              
                ELSE 
                    AIN = 0
                ENDIF    
                IF(ASETTLE(K,I,JA) >= 0)THEN    ! outgoing velocity for BOTTOM layer cell
                    AOUTB = ASETTLE(K,I,JA)
                ELSE 
                    AOUTA = -ASETTLE(K,I,JA)
                ENDIF      
                ASR(K,I,JA) = (AIN*ALG(K-1,I,JA) - (LOSS_FRAC(MIGI)*AOUTB + AOUTA)*ALG(K,I,JA))*BI(K,I)/BH2(K,I)      
            ELSE
                IF(ASETTLE(K-1,I,JA) >= 0)THEN    ! incoming velocity from cell ABOVE
                    AINA = ASETTLE(K-1,I,JA)              
                ELSE 
                    AINA = 0
                ENDIF    
                IF(ASETTLE(K+1,I,JA) >= 0)THEN    ! incoming velocity from cell below
                    AINB = 0                
                ELSE 
                    AINB = -ASETTLE(K+1,I,JA)
                ENDIF  
                AOUT = ABS(ASETTLE(K,I,JA))
                ASR(K,I,JA) =  (AINA*ALG(K-1,I,JA)-AOUT*ALG(K,I,JA))*BI(K,I)/BH2(K,I) + AINB*ALG(K+1,I,JA)*BI(K+1,I)/BH2(K,I)    ! NOT NECESSARY TO DO THE DIVISION - JUST DIVIDE BY H !SP 8/27/07
            ENDIF
      end do
      ENDIF
    end do
    ENDIF
  END DO    ! ALGAE LOOP
   if(algae_settling_exist) iday=jday
      
! Macrophyte Light/Nutrient Limitation and kinetic rates
  do m=1,nmc
  mGR(:,:,iu:id,m)=0.0; mRR(:,iu:id,m)=0.0; mmR(:,iu:id,m)=0.0  ! cb 3/8/16
  if(macrophyte_calc(jw,m))then
    DO I=IU,ID
      LTCOEFm = (1.0-BETA(jw))*SRON(jw)*SHADE(I)
      if(kticol(i))then
        jt=kti(i)
      else
        jt=kti(i)+1
      end if
      je=kb(i)
      do jj=jt,je
        lam1=ltcoefm
        lam2m(jj,kt)=lam1*exp(-gammaj(jj,kt,i)*h2(kt,i))
        lavg=(lam1-lam2m(jj,kt))/(GAMMAj(jj,kt,i)*H2(kt,i))
        mLLIM(jj,kt,I,m) = lavg/(lavg+msat(m))
        IF (mHSP(m)  /= 0.0.and.psed(m) < 1.0)then
          mPLIM(kt,I,m) =  FDPO4*PO4(kt,I)/(FDPO4*PO4(kt,I)+mHSP(m)+nonzero)
        else
          mPLIM(kt,I,m)=1.0
        end if
        IF (mHSN(m)  /= 0.0.and.nsed(m) < 1.0)then
          mNLIM(kt,I,m) = NH4(kt,I)/(NH4(kt,I)+mHSN(m)+nonzero)
        else
          mNLIM(kt,I,m)=1.0
        end if
        IF (mHSc(m) /= 0.0)then
          mcLIM(kt,i,m) = co2(kt,I)/(co2(kt,I)+mHSc(m)+NONZERO)
        end if
        LIMIT          = MIN(mPLIM(kt,I,m),mNLIM(kt,I,m),mcLIM(kt,I,m),mLLIM(jj,kt,I,m))

!************* sources/sinks

        mGR(jj,Kt,I,m) = macTRM(Kt,I,m)*mG(m)*LIMIT

      end do

      mRR(Kt,I,m) = macTRM(Kt,I,m)*mR(m)*DO3(Kt,I)
      mMR(Kt,I,m) = (macTRMR(Kt,I,m)+1.0-mAcTRMF(Kt,I,m))*mM(m)

      DO K=KT+1,KB(I)
        jt=k
        je=kb(i)
        do jj=jt,je
          lam1=lam2m(jj,k-1)
          lam2m(jj,k)=lam1*exp(-gammaj(jj,k,i)*h2(k,i))
          lavg=(lam1-lam2m(jj,k))/(GAMMAj(jj,k,i)*H2(k,i))
          mLLIM(jj,K,I,m) = lavg/(lavg+msat(m))
          IF (mHSP(m)  /= 0.0.and.psed(m) < 1.0)then
            mPLIM(K,I,m) =  FDPO4*PO4(K,I)/(FDPO4*PO4(K,I)+mHSP(m)+nonzero)
          else
            mPLIM(K,I,m)=1.0
          end if
          IF (mHSN(m)  /= 0.0.and.nsed(m) < 1.0)then
            mNLIM(K,I,m) = NH4(K,I)/(NH4(K,I)+mHSN(m)+nonzero)
          else
             mNLIM(K,I,m)=1.0
          end if
          IF (mHSc(m) /= 0.0)then
            mcLIM(k,i,m) = co2(K,I)/(co2(K,I)+mHSc(m)+NONZERO)
          end if
          LIMIT          = MIN(mPLIM(K,I,m),mNLIM(K,I,m),mcLIM(K,I,m),mLLIM(jj,K,I,m))

!************* sources/sinks

          mGR(jj,K,I,m) = macTRM(K,I,m)*mG(m)*LIMIT

        end do

        mRR(K,I,m) = macTRM(K,I,m)*mR(m)*DO3(K,I)
        mMR(K,I,m) = (macTRMR(K,I,m)+1.0-mAcTRMF(K,I,m))*mM(m)
      end do
    END DO
    ENDIF
  END DO

RETURN

!***********************************************************************************************************************************
!**                                             G E N E R I C   C O N S T I T U E N T                                             **
!***********************************************************************************************************************************

ENTRY GENERIC_CONST (JG)
    
  XX=0.0
    DO I=IU,ID
    LIGHT =  (1.0-BETA(JW))*SRON(JW)*SHADE(I)                  !LCJ 2/26/15
    LAM1  =  LIGHT
    LAM2  =  LIGHT

      DO K=KT,KB(I)
        LAM1           = LAM2
        LAM2           = LAM1*EXP(-GAMMA(K,I)*H2(K,I))
        LIGHT          = LAM1*(1.-EXP(-GAMMA(K,I)*H2(K,I)))/(GAMMA(K,I)*H2(K,I))     ! SW 10/17/15

        IF (CGS(JG) >= 0.0) THEN
            IF(K == KT)THEN
            XX =  CGS(JG)*(-CG(K,I,JG))*BI(K,I)/BH2(K,I)    ! AS(JA)*(-ALG(K,I,JA))*BI(K,I)/BH2(K,I)
            ELSE
            XX =  CGS(JG)*(CG(K-1,I,JG)-CG(K,I,JG))*BI(K,I)/BH2(K,I)     !AS(JA)*(ALG(K-1,I,JA)-ALG(K,I,JA))*BI(K,I)/BH2(K,I)
          END IF
        ELSE IF(CGS(JG)<0.0)then
          IF(K == KB(I))THEN
            XX = -CGS(JG)*(-CG(K,I,JG))*BI(K,I)/BH2(K,I)    !-AS(JA)*(-ALG(K,I,JA)  *BI(K,I)/BH2(K,I))                                           !SW 11/8/07
          ELSE IF(K == KT)THEN
            XX = -CGS(JG)*CG(K+1,I,JG)*BI(K+1,I)*DLX(I)/VOL(K,I)    !-AS(JA)* ALG(K+1,I,JA)*BI(K+1,I)*DLX(I)/VOL(K,I)                                   !SW 11/8/07
          ELSE
            XX = -CGS(JG)*(CG(K+1,I,JG)*BI(K+1,I)/BH2(K,I)-CG(K,I,JG)*BI(K,I)/BH2(K,I))    !-AS(JA)*(ALG(K+1,I,JA)*BI(K+1,I)/BH2(K,I)-ALG(K,I,JA)*BI(K,I)/BH2(K,I))             !SP 8/27/07
          END IF
        END IF
        !
         IF (CGQ10(JG) /= 0.0) THEN
          CGSS(K,I,JG) = -CG0DK(JG)*CGQ10(JG)**(T1(K,I)-20.0)-CG1DK(JG)*CGQ10(JG)**(T1(K,I)-20.0)*CG(K,I,JG)-CGLDK(JG)*LIGHT*CG(K,I,JG)+XX
         ELSE
          CGSS(K,I,JG) = -CG0DK(JG)-CG1DK(JG)*CG(K,I,JG)-CGLDK(JG)*LIGHT*CG(K,I,JG)+XX
         ENDIF
         IF(CGR(JG) /= 0.0)CGSS(K,I,JG) = CGSS(K,I,JG) + CGR(JG)*SODD(K,I)*DO2(K,I)
     END DO
    !
     IF(CGKLF(JG) /= 0.0)THEN
         IF (.NOT. ICE(I)) THEN                                                                             ! REAER in units of s-1
         CGSS(KT,I,JG) = CGSS(KT,I,JG)+REAER(I)*CGKLF(JG)*(CGCS(JG)-CG(KT,I,JG))                            ! this calculation performed in Gas-transfer.f90: *BI(KT,I)/BH2(KT,I)
        END IF    
     ENDIF
    !CGSS(K,I,JG) = CGSS(K,I,JG) + CGR(JG)*SODD(K,I)*DO2(K,I)
  END DO

RETURN

!***********************************************************************************************************************************
!**                                               S U S P E N D E D   S O L I D S                                                 **
!***********************************************************************************************************************************

ENTRY SUSPENDED_SOLIDS (J)
    
    If(IncludeBedConsolidation)Then
        !All resuspension done in CEMA code
        SEDIMENT_RESUSPENSION(J) = .FALSE.
    End If
    
  DO I=IU,ID
    SSR = 0.0
    IF (SEDIMENT_RESUSPENSION(J)) THEN
      FETCH = FETCHD(I,JB)
      IF (COS(PHI(JW)-PHI0(I)) < 0.0) FETCH = FETCHU(I,JB)
      FETCH = MAX(FETCH,BI(KT,I),DLX(I))
      U2    = WIND(JW)*WSC(I)*WIND(JW)*WSC(I)+NONZERO
      COEF1 = 0.53  *(G*DEPTHB(KT,I)/U2)**0.75
      COEF2 = 0.0125*(G*FETCH/U2)**0.42
      COEF3 = 0.833* (G*DEPTHB(KT,I)/U2)**0.375
      COEF4 = 0.077* (G*FETCH/U2)**0.25
      HS    = U2/G*0.283*TANH(COEF1)*TANH(COEF2/TANH(COEF1))    !8/2020 corrected
      !TS    = 2.0*PI*U2/G*1.2*  TANH(COEF3)*TANH(COEF4/TANH(COEF3))
      TS    = 2.0*PI*sqrt(U2)/G*1.2*  TANH(COEF3)*TANH(COEF4/TANH(COEF3))   ! cb 7/15/14
      L0    = G*TS*TS/(2.0*PI)
        L1 = L0             ! SW 6/28/2018 Allow for resuspension of surface layer
        L  = L0*TANH(2.0*PI*DEPTHB(KT,I)/L1)
        DO WHILE (ABS(L-L1) > 0.001)
          L1 = L
          L  = L0*TANH(2.0*PI*DEPTHB(KT,I)/L1)
        END DO
        COEF = MIN(710.0,2.0*PI*DEPTHB(KT,I)/L)
        UORB = PI*HS/TS*100.0/SINH(COEF)
        TAU  = 0.003*UORB*UORB
        IF (TAU-TAUCR(J) > 0.0) EPSILON = MAX(0.0,0.008/49.0*(TAU-TAUCR(J))**3*10000.0/DLT)
        SSR = EPSILON*DLX(I)*(BI(KT,I)-BI(KT+1,I))/VOL(KT,I)
      END IF
    SSSS(KT,I,J) = -SSS(J)*SS(KT,I,J)*BI(KT,I)/BH2(KT,I)+SSR
 !   DO K=KT-1,KB(I)-1                                             ! SW 4/3/09   KT,KB
     DO K=KT+1,KB(I)  !-1                 ! sw 12/2020 cb 9/29/14
      IF (SEDIMENT_RESUSPENSION(J)) THEN
        L1 = L0
        L  = L0*TANH(2.0*PI*DEPTHB(K,I)/L1)
        DO WHILE (ABS(L-L1) > 0.001)
          L1 = L
          L  = L0*TANH(2.0*PI*DEPTHB(K,I)/L1)
        END DO
        COEF = MIN(710.0,2.0*PI*DEPTHB(K,I)/L)
        UORB = PI*HS/TS*100.0/SINH(COEF)
        TAU  = 0.003*UORB*UORB
        IF (TAU-TAUCR(J) > 0.0) EPSILON = MAX(0.0,0.008/49.0*(TAU-TAUCR(J))**3*10000.0/DLT)
		if(k == kb(i))then   ! SW 4/18/07
		SSR = EPSILON*DLX(I)*BI(K,I)/VOL(K,I)
		else
        SSR = EPSILON*DLX(I)*(BI(K,I)-BI(K+1,I))/VOL(K,I)
		endif
      END IF
      SSSS(K,I,J) = SSS(J)*(SS(K-1,I,J)-SS(K,I,J))*BI(K,I)/BH2(K,I)+SSR
    END DO
    !IF (SEDIMENT_RESUSPENSION(J)) SSR = EPSILON*DLX(I)*BI(KB(I),I)/VOL(KB(I),I)

    If(IncludeFFTLayer .and. FFTActive)Then
        SSSS(KB(I),I,J) = (SSS(J)*SS(KB(I)-1,I,J)-FFTLayerSettVel*SS(KB(I),I,J))/H(KB(I),JW)+SSR
    End If
    If(IncludeFFTLayer .and. .NOT. FFTActive)Then
        SSSS(KB(I),I,J) = 0.d0
    End If
    If(.NOT. IncludeFFTLayer)Then
       SSSS(KB(I),I,J) = SSS(J)*(SS(KB(I)-1,I,J)-SS(KB(I),I,J))/H(KB(I),JW)+SSR
  ! Flocculation              !SR                                                      !New section on flocculation          !SR 04/21/13
    !DO K=KT,KB(I)
    !  SSF = 0.0
    !  IF (J > 1 .AND. SSFLOC(J-1) > 0.0) THEN
    !    IF (FLOCEQN(J-1) == 0) THEN
    !      SSF = MIN(SSFLOC(J-1), SS(K,I,J-1)/DLT)
    !    ELSE IF (FLOCEQN(J-1) == 1) THEN
    !      SSF = SSFLOC(J-1)*SS(K,I,J-1)
    !    ELSE IF (FLOCEQN(J-1) == 2) THEN
    !      SSF = SSFLOC(J-1)*SS(K,I,J-1)*SS(K,I,J-1)
    !    END IF
    !  END IF
    !  IF (J < NSS .AND. SSFLOC(J) > 0.0) THEN
    !    IF (FLOCEQN(J) == 0) THEN
    !      SSF = SSF - MIN(SSFLOC(J), SS(K,I,J)/DLT)
    !    ELSE IF (FLOCEQN(J) == 1) THEN
    !      SSF = SSF - SSFLOC(J)*SS(K,I,J)
    !    ELSE IF (FLOCEQN(J) == 2) THEN
    !      SSF = SSF - SSFLOC(J)*SS(K,I,J)*SS(K,I,J)
    !    END IF
    !  END IF
    !  SSSS(K,I,J) = SSSS(K,I,J) + SSF
    !END DO                                                                        !End new section on flocculation      !SR 04/21/13
    End If
    
  END DO
RETURN

!***********************************************************************************************************************************
!**                                               WATER AGE                                                                       **
!***********************************************************************************************************************************
ENTRY WATER_AGE
  DO I=IU,ID
    DO K=KT,KB(I)
      AGESS(K,I) = 1.0/DAY
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                               BACTERIA                                                                        **
!***********************************************************************************************************************************
ENTRY BACTERIA
  DO I=IU,ID
    LIGHT = (1.0-BETA(JW))*SRON(JW)*SHADE(I)                
    LAM1  = LIGHT
    LAM2  = LIGHT
    DO K=KT,KB(I)
      LAM1  = LAM2
      LAM2  = LAM1*EXP(-GAMMA(K,I)*H2(K,I))
      LIGHT = LAM1*(1.-EXP(-GAMMA(K,I)*H2(K,I)))/(GAMMA(K,I)*H2(K,I))  
      ! SETTELING
      IF(BACTS(JW) >= 0.0)THEN   
        IF(K == KT)THEN
          SET_BACT =  BACTS(JW)*(-BACT(K,I))*BI(K,I)/BH2(K,I) 
        ELSE
          SET_BACT =  BACTS(JW)*(BACT(K-1,I)-BACT(K,I))*BI(K,I)/BH2(K,I)
        END IF
      ELSE IF(BACTS(JW)<0.0)THEN
        IF(K == KB(I))THEN
          SET_BACT = -BACTS(JW)*(-BACT(K,I))*BI(K,I)/BH2(K,I)                                         
        ELSE IF(K == KT)THEN
          SET_BACT = -BACTS(JW)*BACT(K+1,I)*BI(K+1,I)*DLX(I)/VOL(K,I)                                    
        ELSE
          SET_BACT = -BACTS(JW)*(BACT(K+1,I)*BI(K+1,I)/BH2(K,I)-BACT(K,I)*BI(K,I)/BH2(K,I))   
        END IF
      END IF
      !FIRST-ORDER DECAY
      IF(BACTQ10(JW) /= 0.0)THEN
        DK1_BACT = -BACT1DK(JW)*BACTQ10(JW)**(T1(K,I)-20.0)*BACT(K,I)*DO3(K,I) 
      ELSE
        DK1_BACT = -BACT1DK(JW)*BACT(K,I)*DO3(K,I) 
      END IF
      !PHOTODEGRADATION
      PHOTO_BACT = -BACTLDK(JW) * LIGHT * BACT(K,I) 
      !
      BACTSS(K,I) = SET_BACT + DK1_BACT + PHOTO_BACT
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  DISSOLVED GAS PRESSURE (atm)                                                 **
!***********************************************************************************************************************************
ENTRY DISSOLVED_GAS
  DISGSS(:,IU:ID)=0.0
  DO I=IU,ID
    IF (.NOT. ICE(I)) THEN
        ! REAER in units of s-1
      DISGSS(KT,I) = -REAER(I)*DGPO2(JW)*(DGP(KT,I)-PALT(I))      ! this computation done in gas-transfer.f90 :*BI(KT,I)/BH2(KT,I)   
    END IF
    !
    DO K=KT,KB(I)
      TDG(K,I) = 100.*DGP(K,I)/PALT(I)
    END DO 
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      DISSOLVED NITROGEN GAS (N2)                                              **
!***********************************************************************************************************************************
ENTRY DISSOLVED_N2
  N2SS(:,IU:ID)=0.0
  DO I=IU,ID
     IF(.NOT. ICE(I))THEN
      EA = DEXP(2.3026D0*(7.5D0*TDEW(JW)/(TDEW(JW)+237.3D0)+0.6609D0))*0.001316   ! EA (atm)   0.0098692atm=7.5006151mmHg   
      N2SAT = 1.5568D06*0.79*(PALT(I)-EA)*(1.8816D-5 - 4.116D-7 * T1(KT,I) + 4.6D-9 * T1(KT,I)**2)
      N2SS(KT,I) = 1.034*REAER(I)*(N2SAT-N2(KT,I))      !This computation is done in gas-transfer.f90 *BI(KT,I)/BH2(KT,I)             ! KLN2=1.034*KLO2    ! ! REAER in units of s-1
      !
      DO K=KT,KB(I)
        DOSAT = SATO(T1(K,I),TDS(K,I),PALT(I),SALT_WATER(JW))
        TDG(K,I) = 100.*((0.79*N2(K,I)/N2SAT) + O2(K,I)/DOSAT*0.21)
      END DO 
    END IF
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      SULFIDE (H2S)                                                            **
!***********************************************************************************************************************************
ENTRY SULFIDE
  H2SD(:,IU:ID) = 0.0; H2SREAER(:,IU:ID) = 0.0; H2SSR(:,IU:ID) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I) 
      IF(H2SQ10(JW) /= 0.0)THEN
        H2SD(K,I) = -H2S1DK(JW)*H2SQ10(JW)**(T1(K,I)-20.0)*H2S(K,I)*DO3(K,I) 
      ELSE
        H2SD(K,I) = -H2S1DK(JW)*H2S(K,I)*DO3(K,I) 
      END IF
      H2SSR(K,I)  = H2SR(JW)*SODD(K,I)*DO2(K,I)
      H2SSS(K,I)  = H2SD(K,I) + H2SSR(K,I)
    END DO
    !
    IF(.NOT. ICE(I))THEN
      H2SREAER(KT,I) = 0.984*REAER(I)*(-H2S(KT,I))                  !This computation is done in gas-transfer.f90: *BI(KT,I)/BH2(KT,I)      REAER in units of s-1
      H2SSS(KT,I) = H2SSS(KT,I)+H2SREAER(KT,I)
    END IF
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      METHANE (CH4)                                                            **
!***********************************************************************************************************************************
ENTRY METHANE
  CH4D(:,IU:ID) = 0.0; CH4REAER(:,IU:ID) = 0.0; CH4SR(:,IU:ID) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      IF(CH4Q10(JW) /= 0.0)THEN
        CH4D(K,I) = -CH41DK(JW)*CH4Q10(JW)**(T1(K,I)-20.0)*CH4(K,I)  *DO3(K,I) 
      ELSE
        CH4D(K,I) = -CH41DK(JW)*CH4(K,I) *DO3(K,I) 
      END IF
      CH4SR(K,I)  = CH4R(JW)*SODD(K,I)*DO2(K,I)
      CH4SS(K,I)  = CH4D(K,I) + CH4SR(K,I)
    END DO
    !
    IF(.NOT. ICE(I))THEN
      CH4REAER(KT,I) = 1.188*REAER(I)*(-CH4(KT,I))              !This computation is done in gas-transfer.f90:*BI(KT,I)/BH2(KT,I)      REAER in units of s-1
      CH4SS(KT,I) = CH4SS(KT,I)+CH4REAER(KT,I)
    END IF
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      SULFATE (SO4)                                                            **
!***********************************************************************************************************************************
ENTRY SULFATE
  DO I=IU,ID
    DO K=KT,KB(I)
      SO4SS(K,I) = -H2SD(K,I) + SO4R(JW)*SODD(K,I)*DO2(K,I)   ! sulfate production from sulfide decay
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      FERROUS (FEII)                                                           **
!***********************************************************************************************************************************
ENTRY FERROUS
  FE2D(:,IU:ID) = 0.0; FEIISR(:,IU:ID) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      
      IF(FEII(K,I) > 1.0E-7)THEN   ! skip calculations if FE2 is negligible anyway
      IF(.NOT. PH_CALC(JW))THEN
          FE2D(K,I)    = -KFE_OXID(JW)*O2(K,I)*FEII(K,I)    ! assuming pH(K,I) = 7.0
      ELSEIF(PH(K,I) <= 4.0)THEN
          FE2D(K,I)    = -KFE_OXID(JW)*O2(K,I)*(1.0E-6)*FEII(K,I) 
      ELSEIF(PH(K,I) >= 8.0)THEN
          FE2D(K,I)    = -KFE_OXID(JW)*O2(K,I)*(1.0E+2)*FEII(K,I) 
      ELSE   
          FE2D(K,I)    = -KFE_OXID(JW)*O2(K,I)*10**(2.0*(pH(K,I)-7.0))*FEII(K,I) 
      ENDIF
      ENDIF
      
      FEIISR(K,I)  = FEIIR(JW)*SODD(K,I)*DO2(K,I)

      FEIISS(K,I)  = FEIISR(K,I) + FE2D(K,I) + KFE_RED(JW)*(KFEOOH_HalfSat(JW)/(O2(K,I)+KFEOOH_HalfSat(JW)))*FEOOH(K,I)
    END DO 
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      OXIDIZEDFE (FEOOH)                                                       **
!***********************************************************************************************************************************
ENTRY OXIDIZEDFE
  SDINFEOOH(:,IU:ID)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)  
      IF(K == KT)THEN
        SDINFEOOH(K,I) = FeSetVel(JW)*(-FEOOH(K,I))*BI(K,I)/BH2(K,I)      !(FEOOH(K-1,I)-FEOOH(K,I))*BI(K,I)/BH2(K,I) 
      ELSE
        SDINFEOOH(K,I) = FeSetVel(JW)*(FEOOH(K-1,I)-FEOOH(K,I))*BI(K,I)/BH2(K,I) 
      END IF
      FEOOHSS(K,I) = -FE2D(K,I) - KFE_RED(JW)*(KFEOOH_HalfSat(JW)/(O2(K,I)+KFEOOH_HalfSat(JW)))*FEOOH(K,I) + SDINFEOOH(K,I)
    END DO 
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                     BIVALENTMN (MNII)                                                         **
!***********************************************************************************************************************************
ENTRY BIVALENTMN
  MN2D(:,IU:ID) = 0.0; MNIISR(:,IU:ID) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      IF(MNII(K,I) > 1.0E-7) THEN   ! SKIP CALCS IF MNII IS NEGLIGIBLE ANYWAY
      IF(.NOT. PH_CALC(JW))THEN
          MN2D(K,I)   = -KMN_OXID(JW)*O2(K,I)*MNII(K,I)    ! assuming pH(K,I) = 7.0
      ELSEIF(PH(K,I) <= 4.0)THEN
          MN2D(K,I)   =  -KMN_OXID(JW)*O2(K,I)*(1.0E-6)*MNII(K,I)
      ELSEIF(PH(K,I) >= 8.0)THEN
          MN2D(K,I)   =   -KMN_OXID(JW)*O2(K,I)*(1.0E+2)*MNII(K,I)
      ELSE   
          MN2D(K,I)   = -KMN_OXID(JW)*O2(K,I)*10**(2.0*(pH(K,I)-7.0))*MNII(K,I)
      ENDIF
      ENDIF
      MNIISR(K,I) = MNIIR(JW)*SODD(K,I)*DO2(K,I)
      MNIISS(K,I) = MNIISR(K,I) + MN2D(K,I) + KMN_RED(JW)*(KMNO2_HalfSat(JW)/(O2(K,I)+KMNO2_HalfSat(JW)))*MNO2(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      OXIDIZEDMN (MNO2)                                                        **
!***********************************************************************************************************************************
ENTRY OXIDIZEDMN
  SDINMNO2(:,IU:ID)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      IF(K == KT)THEN
        SDINMNO2(K,I) = MnSetVel(JW)*(-MNO2(K,I))*BI(K,I)/BH2(K,I) 
      ELSE
        SDINMNO2(K,I) = MnSetVel(JW)*(MNO2(K-1,I)-MNO2(K,I))*BI(K,I)/BH2(K,I) 
      END IF
      MNO2SS(K,I) = -MN2D(K,I) - KMN_RED(JW)*(KMNO2_HalfSat(JW)/(O2(K,I)+KMNO2_HalfSat(JW)))*MNO2(K,I) + SDINMNO2(K,I)     
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      P H O S P H O R U S                                                      **
!***********************************************************************************************************************************

ENTRY PHOSPHORUS
  PO4AR(:,IU:ID) = 0.0; PO4AG(:,IU:ID) = 0.0; PO4ER(:,IU:ID) = 0.0; PO4EG(:,IU:ID) = 0.0; PO4BOD(:,IU:ID) = 0.0
  PO4MR(:,IU:ID) = 0.0; PO4MG(:,IU:ID) = 0.0; PO4ZR(:,IU:ID)=0.0  

  DO I=IU,ID
    DO K=KT,KB(I)
      DO JCB=1,NBOD
!        IF(BOD_CALC(JCB))PO4BOD(K,I) = PO4BOD(K,I)+CBODD(K,I,JCB)*CBOD(K,I,JCB)*BODP(JCB)
         IF(BOD_CALCp(JCB))then                                                ! cb 5/19/11
           PO4BOD(K,I) = PO4BOD(K,I)+CBODD(K,I,JCB)*CBODp(K,I,JCB)    
         else
           PO4BOD(K,I) = PO4BOD(K,I)+CBODD(K,I,JCB)*CBOD(K,I,JCB)*BODP(JCB)
         end if
      END DO
      DO JA=1,NAL
        IF(ALG_CALC(JA))THEN
        PO4AG(K,I) = PO4AG(K,I)+AGR(K,I,JA)*ALG(K,I,JA)*AP(JA)
        PO4AR(K,I) = PO4AR(K,I)+ARR(K,I,JA)*ALG(K,I,JA)*AP(JA)
        ENDIF
      END DO
      DO JE=1,NEP
      IF (EPIPHYTON_CALC(JW,JE))then
        PO4EG(K,I) = PO4EG(K,I)+EGR(K,I,JE)*EPC(K,I,JE)*EP(JE)
        PO4ER(K,I) = PO4ER(K,I)+ERR(K,I,JE)*EPC(K,I,JE)*EP(JE)
      endif
      END DO
      PO4EP(K,I)  = PO4ER(K,I)-PO4EG(K,I)
      PO4AP(K,I)  = PO4AR(K,I)-PO4AG(K,I)
      !
      PO4POM(K,I) = LPOMPD(K,I)+RPOMPD(K,I)     
      PO4DOM(K,I) = LDOMPD(K,I)+RDOMPD(K,I)     
      PO4OM(K,I)  = PO4POM(K,I)+PO4DOM(K,I)
            IF(STANDING_BIOMASS_DECAY)THEN  ! SW 5/26/15
           ! PO4SD(K,I)  = SEDDp(K,I)+sedd1(k,i)*orgp(jw) + sedd2(k,i)*orgp(jw)   ! Amaila
            PO4SD(K,I)  = SEDDp(K,I)+sedd1(k,i)*pbiom(jw) + sedd2(k,i)*pbiom(jw)   ! Amaila, cb 6/7/17
            ELSE
            PO4SD(K,I)  = SEDDp(K,I)
            ENDIF
      PO4SR(K,I)  = PO4R(JW)*SODD(K,I)*DO2(K,I)
      PO4NS(K,I)  = (SSSI(K,I)*PO4(K-1,I)-SSSO(K,I)*PO4(K,I))*BI(K,I)/BH2(K,I)   !8/2020 corrected
      
      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          IF(K.EQ.KT)THEN
            JT=KTI(I)
          ELSE
            JT=K
          END IF
          JE=KB(I)
          DO JJ=JT,JE
            PO4MG(K,I)= PO4MG(K,I)+MGR(JJ,K,I,M)*MACRM(JJ,K,I,M)*MP(M)*(1.0-PSED(M))
            PO4MR(K,I)= PO4MR(K,I)+MRR(K,I,M)*MACRM(JJ,K,I,M)*MP(M)
          END DO
        END IF
      END DO
      PO4MR(K,I) = PO4MR(K,I)/(DLX(I)*BH2(K,I)) !8/2020 corrected
      PO4MG(K,I) = PO4MG(K,I)/(DLX(I)*BH2(K,I)) !
      IF(ZOOPLANKTON_CALC)THEN
      DO JZ = 1,NZP
        PO4ZR(K,I) = PO4ZR(K,I) + ZRT(K,I,JZ)*ZOO(K,I,JZ)*ZP(JZ)
	  END DO
	  ENDIF

      PO4SS(K,I)  = PO4AP(K,I)+PO4EP(K,I)+PO4OM(K,I)+PO4SD(K,I)+PO4SR(K,I)+PO4NS(K,I)+PO4BOD(K,I)  &
                    +PO4MR(K,I)-PO4MG(K,I) +PO4ZR(K,I)    

    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                        A M M O N I U M                                                        **
!***********************************************************************************************************************************

ENTRY AMMONIUM
  NH4AG(:,IU:ID) = 0.0; NH4AR(:,IU:ID) = 0.0; NH4ER(:,IU:ID) = 0.0; NH4EG(:,IU:ID) = 0.0; NH4BOD(:,IU:ID) = 0.0
  NH4MG(:,IU:ID) = 0.0; NH4MR(:,IU:ID) = 0.0; NH4ZR(:,IU:ID)=0.0   
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JCB=1,NBOD
!        IF(BOD_CALC(JCB))NH4BOD(K,I) =  NH4BOD(K,I)+CBODD(K,I,JCB)*CBOD(K,I,JCB)*BODN(JCB)
         IF(BOD_CALCn(JCB))then                                                ! cb 5/19/11
           NH4BOD(K,I) =  NH4BOD(K,I)+CBODD(K,I,JCB)*CBODn(K,I,JCB)
         else
           NH4BOD(K,I) =  NH4BOD(K,I)+CBODD(K,I,JCB)*CBOD(K,I,JCB)*BODN(JCB)
         end if
      END DO
      DO JA=1,NAL
      IF(ALG_CALC(JA))THEN
        IF (ANEQN(JA).EQ.2) THEN
        NH4PR      = NH4(K,I)*NO3(K,I)/((ANPR(JA)+NH4(K,I))*(ANPR(JA)+NO3(K,I)))+NH4(K,I)*ANPR(JA)/((NO3(K,I)  &
                                        +NH4(K,I)+NONZERO)*(ANPR(JA)+NO3(K,I)))
        ELSE
        NH4PR = NH4(K,I)/(NH4(K,I)+NO3(K,I)+NONZERO)
        ENDIF
        IF (AHSN(JA) > 0.0) NH4AG(K,I) = NH4AG(K,I)+AGR(K,I,JA)*ALG(K,I,JA)*AN(JA)*NH4PR
        NH4AR(K,I) = NH4AR(K,I)+ARR(K,I,JA)*ALG(K,I,JA)*AN(JA)
      ENDIF
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))then
        IF (ENEQN(JE) == 2)THEN
        NH4PR = NH4(K,I)*NO3(K,I)/((ENPR(JE)+NH4(K,I))*(ENPR(JE)+NO3(K,I)))+NH4(K,I)*ENPR(JE)/((NO3(K,I)  &
                                        +NH4(K,I)+NONZERO)*(ENPR(JE)+NO3(K,I)))
        ELSE
        NH4PR = NH4(K,I)/(NH4(K,I)+NO3(K,I)+NONZERO)
        ENDIF
        IF (EHSN(JE) > 0.0) NH4EG(K,I) = NH4EG(K,I)+EGR(K,I,JE)*EPC(K,I,JE)*EN(JE)*NH4PR
        NH4ER(K,I) = NH4ER(K,I)+ERR(K,I,JE)*EPC(K,I,JE)*EN(JE)
        endif
      END DO
      NH4EP(K,I)  =  NH4ER(K,I) -NH4EG(K,I)
      NH4AP(K,I)  =  NH4AR(K,I) -NH4AG(K,I)
      !
      NH4DOM(K,I) = LDOMND(K,I)+RDOMND(K,I)  
      NH4POM(K,I) = LPOMND(K,I)+RPOMND(K,I)  
      NH4OM(K,I)  =  NH4DOM(K,I)+NH4POM(K,I)

            IF(STANDING_BIOMASS_DECAY)THEN  ! SW 5/26/15
            !NH4SD(K,I)  =  SEDDn(K,I) +sedd1(k,i)*orgn(jw) + sedd2(k,i)*orgn(jw)   ! Amaila
            NH4SD(K,I)  =  SEDDn(K,I) +sedd1(k,i)*nbiom(jw) + sedd2(k,i)*nbiom(jw)   ! Amaila, cb 6/7/17
            ELSE
            NH4SD(K,I)  =  SEDDn(K,I)
            ENDIF

      NH4SR(K,I)  =  NH4R(JW) *SODD(K,I)*DO2(K,I)

      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          IF(K.EQ.KT)THEN
            JT=KTI(I)
          ELSE
            JT=K
          END IF
          JE=KB(I)
          DO JJ=JT,JE
            NH4MR(K,I)= NH4MR(K,I)+MRR(K,I,M)*MACRM(JJ,K,I,M)*MN(M)
            NH4MG(K,I)= NH4MG(K,I)+MGR(JJ,K,I,M)*MACRM(JJ,K,I,M)*MN(M)*(1.0-NSED(M))
          END DO
        END IF
      END DO
      NH4MR(K,I) = NH4MR(K,I)/(DLX(I)*BH2(K,I))  !8/2020 corrected
      NH4MG(K,I) = NH4MG(K,I)/(DLX(I)*BH2(K,I))  ! 
	  IF(ZOOPLANKTON_CALC)THEN
	  DO JZ = 1,NZP
	    NH4ZR(K,I) = NH4ZR(K,I) + ZRT(K,I,JZ)*ZOO(K,I,JZ)*ZN(JZ) 
	  END DO
	  ENDIF
      NH4SS(K,I)  =  NH4AP(K,I)+NH4EP(K,I)+NH4OM(K,I)+NH4SD(K,I)+NH4SR(K,I)+NH4BOD(K,I)-NH4D(K,I)  &
         +NH4MR(K,I)-NH4MG(K,I) +NH4ZR(K,I)     
      IF(K==KT.AND.CDWBC(NH3_DER,JW)=='      ON')NH4SS(K,I)  = NH4SS(K,I)  - NH3GAS(K,I)     ! NH3 VOLATILIZATION IF PH IS ON AND NH3 DER ON
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                          N I T R A T E                                                        **
!***********************************************************************************************************************************

ENTRY NITRATE
  NO3AG(:,IU:ID) = 0.0; NO3EG(:,IU:ID) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
      IF(ALG_CALC(JA))THEN
        NO3PR = 1.0-NH4(K,I)/(NH4(K,I)+NO3(K,I)+NONZERO)
        IF (ANEQN(JA).EQ.2)  NO3PR      = 1.0-(NH4(K,I)*NO3(K,I)/((ANPR(JA)+NH4(K,I))*(ANPR(JA)+NO3(K,I)))+NH4(K,I)*ANPR(JA)       &
                                          /((NO3(K,I)+NH4(K,I)+NONZERO)*(ANPR(JA)+NO3(K,I))))
        IF (AHSN(JA).GT.0.0) NO3AG(K,I) = NO3AG(K,I)+AGR(K,I,JA)*ALG(K,I,JA)*NO3PR*AN(JA)
      ENDIF
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))then
        NO3PR = 1.0-NH4(K,I)/(NH4(K,I)+NO3(K,I)+NONZERO)
        IF (ENEQN(JE).EQ.2)  NO3PR      = 1.0-(NH4(K,I)*NO3(K,I)/((ENPR(JE)+NH4(K,I))*(ENPR(JE)+NO3(K,I)))+NH4(K,I)*ENPR(JE)       &
                                          /((NO3(K,I)+NH4(K,I)+NONZERO)*(ENPR(JE)+NO3(K,I))))
        IF (EHSN(JE).GT.0.0) NO3EG(K,I) = NO3EG(K,I)+EGR(K,I,JE)*EPC(K,I,JE)*NO3PR*EN(JE)
        ENDIF
      END DO
      IF(K == KB(I)) THEN      ! SW 4/18/07
      NO3SED(K,I) = NO3(K,I)*NO3S(JW)*NO3TRM(K,I)*(BI(K,I))/BH2(K,I)
	  ELSE
      NO3SED(K,I) = NO3(K,I)*NO3S(JW)*NO3TRM(K,I)*(BI(K,I)-BI(K+1,I))/BH2(K,I)
	  ENDIF
      NO3SS(K,I)  = NH4D(K,I)-NO3D(K,I)-NO3AG(K,I)-NO3EG(K,I)-NO3SED(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  D I S S O L V E D   S I L I C A                                              **
!***********************************************************************************************************************************

ENTRY DISSOLVED_SILICA
  DSIAG(:,IU:ID) = 0.0; DSIEG(:,IU:ID) = 0.0                          !; DSIBOD = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
      IF(ALG_CALC(JA))THEN
        DSIAG(K,I) = DSIAG(K,I)+AGR(K,I,JA)*ALG(K,I,JA)*ASI(JA)
      ENDIF
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))DSIEG(K,I) = DSIEG(K,I)+EGR(K,I,JE)*EPC(K,I,JE)*ESI(JE)
      END DO
      DSID(K,I)  =  PSIDK(JW)*PSI(K,I)
      DSISD(K,I) =  SEDD(K,I)*ORGSI(JW)
      DSISR(K,I) =  DSIR(JW)*SODD(K,I)*DO2(K,I)
      DSIS(K,I)  = (SSSI(K,I)*DSI(K-1,I)-SSSO(K,I)*DSI(K,I))*BI(K,I)/BH2(K,I)*PARTSI(JW)    !8/2020 corrected
      DSISS(K,I) =  DSID(K,I)+DSISD(K,I)+DSISR(K,I)+DSIS(K,I)-DSIAG(K,I)-DSIEG(K,I)    !+DSIBOD
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                P A R T I C U L A T E   S I L I C A                                            **
!***********************************************************************************************************************************

ENTRY PARTICULATE_SILICA
  PSIAM(:,IU:ID) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
      IF(ALG_CALC(JA))THEN
        PSIAM(K,I) = PSIAM(K,I)+AMR(K,I,JA)*ALG(K,I,JA)*ASI(JA)     !   PSI(K,I)   HA-Z  12/2016
      ENDIF
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE)) PSIEM(K,I) = PSIEM(K,I)+EMR(K,I,JE)*EPC(K,I,JE)*ESI(JE)                  
      END DO
      PSID(K,I)  = PSIDK(JW)*PSI(K,I)
      PSINS(K,I) = PSIS(JW)*(PSI(K-1,I)*DO1(K-1,I)-PSI(K,I)*DO1(K,I))*BI(K,I)/BH2(K,I)
      PSISS(K,I) = PSIAM(K,I)-PSID(K,I)+PSINS(K,I) + PSIEM(K,I)   !8/2020 corrected
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                       L A B I L E   D O M                                                     **
!***********************************************************************************************************************************

ENTRY LABILE_DOM
  LDOMAP(:,IU:ID) = 0.0; LDOMEP(:,IU:ID) = 0.0; LDOMMAC(:,IU:ID)= 0.0  
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        IF(ALG_CALC(JA))LDOMAP(K,I) = LDOMAP(K,I)+(AER(K,I,JA)+(1.0-APOM(JA))*AMR(K,I,JA))*ALG(K,I,JA)
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))LDOMEP(K,I) = LDOMEP(K,I)+(EER(K,I,JE)+(1.0-EPOM(JE))*EMR(K,I,JE))*EPC(K,I,JE)
      END DO

      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          IF(K.EQ.KT)THEN
            JT=KTI(I)
          ELSE
            JT=K
          END IF
          JE=KB(I)
          DO JJ=JT,JE
            LDOMMAC(K,I)=LDOMMAC(K,I)+(1.0-MPOM(M))*MMR(K,I,M)*MACRM(JJ,K,I,M)
          END DO
        END IF
      END DO
      LDOMMAC(K,I) = LDOMMAC(K,I)/(DLX(I)*BH2(K,I))               !8/2020 
      LDOMSS(K,I) = LDOMAP(K,I)+LDOMEP(K,I)-LDOMD(K,I)-LRDOMD(K,I)+LDOMMAC(K,I)+LPOMHD(K,I)

    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  R E F R A C T O R Y   D O M                                                  **
!***********************************************************************************************************************************

ENTRY REFRACTORY_DOM
  DO I=IU,ID
    DO K=KT,KB(I)
      RDOMSS(K,I) = LRDOMD(K,I)-RDOMD(K,I)+RPOMHD(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      L A B I L E   P O M                                                      **
!***********************************************************************************************************************************

ENTRY LABILE_POM
  LPOMAP(:,IU:ID) = 0.0; LPOMEP(:,IU:ID) = 0.0;   LPOMMAC(:,IU:ID) = 0.0; LPZOOIN(:,IU:ID)=0.0;LPZOOOUT(:,IU:ID)=0.0   ! cb 5/19/06
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        IF(ALG_CALC(JA))LPOMAP(K,I) = LPOMAP(K,I)+APOM(JA)*(AMR(K,I,JA)*ALG(K,I,JA))
      END DO
      DO JE=1,NEP                                                          ! cb 5/19/06
        IF (EPIPHYTON_CALC(JW,JE))LPOMEP(K,I) = LPOMEP(K,I)+EPOM(JE)*(EMR(K,I,JE)*EPC(K,I,JE))       ! cb 5/19/06
      END DO                                                               ! cb 5/19/06
      LPOMNS(K,I) = POMS(JW)*(LPOM(K-1,I)-LPOM(K,I))*BI(K,I)/BH2(K,I)
      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          JT=K
          JE=KB(I)
          DO JJ=JT,JE
            LPOMMAC(K,I)=LPOMMAC(K,I)+MPOM(M)*LRPMAC(M)*MMR(K,I,M)*MACRM(JJ,K,I,M)
          END DO
        END IF
      END DO
      LPOMMAC(K,I) = LPOMMAC(K,I)/(DLX(I)*BH2(K,I))          !8/2020
      IF(ZOOPLANKTON_CALC)THEN
      DO JZ = 1,NZP
        IF(TGRAZE(K,I,JZ) > 0.0)THEN
          LPZOOOUT(K,I)=LPZOOOUT(K,I)+ZOO(K,I,JZ)*(ZMT(K,I,JZ)+(ZMU(K,I,JZ)-(ZMU(K,I,JZ)*ZEFF(JZ))))
          LPZOOIN(K,I)=LPZOOIN(K,I)+ZOO(K,I,JZ)*PREFP(JZ)*ZMU(K,I,JZ)*LPOM(K,I)/TGRAZE(K,I,JZ)
        ELSE
          LPZOOOUT(K,I)=LPZOOOUT(K,I)+ZOO(K,I,JZ)*(ZMT(K,I,JZ)+(ZMU(K,I,JZ)-(ZMU(K,I,JZ)*ZEFF(JZ))))
          LPZOOIN(K,I)=0.0
        END IF
      END DO
      ENDIF
      LPOMSS(K,I) = LPOMAP(K,I)+LPOMEP(K,I)-LPOMD(K,I)+LPOMNS(K,I)-LRPOMD(K,I)+LPOMMAC(K,I)+LPZOOOUT(K,I)-LPZOOIN(K,I)-LPOMHD(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  R E F R A C T O R Y   P O M                                                  **
!***********************************************************************************************************************************

ENTRY REFRACTORY_POM
  RPOMMAC(:,IU:ID) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      RPOMNS(K,I) = POMS(JW)*(RPOM(K-1,I)-RPOM(K,I))*BI(K,I)/BH2(K,I)
      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          JT=K
          JE=KB(I)
          DO JJ=JT,JE
            RPOMMAC(K,I)=RPOMMAC(K,I)+MPOM(M)*(1.0-LRPMAC(M))*MMR(K,I,M)*MACRM(JJ,K,I,M)
          END DO
        END IF
      END DO
      RPOMMAC(K,I) = RPOMMAC(K,I)/(DLX(I)*BH2(K,I))              !8/2020   
      RPOMSS(K,I) = LRPOMD(K,I)+RPOMNS(K,I)-RPOMD(K,I)+RPOMMAC(K,I)-RPOMHD(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                         A L G A E                                                             **
!***********************************************************************************************************************************

ENTRY ALGAE (J)
  AGZT(:,IU:ID,J) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      IF(ZOOPLANKTON_CALC)THEN
      DO JZ = 1,NZP
	  AGZT(K,I,J) = AGZT(K,I,J) + AGZ(K,I,J,JZ)                       ! CB 5/26/07
	  END DO
	  ENDIF
      ASS(K,I,J) = ASR(K,I,J)+(AGR(K,I,J)-AER(K,I,J)-AMR(K,I,J)-ARR(K,I,J))*ALG(K,I,J)-AGZT(K,I,J)	
    END DO
  END DO
RETURN

RETURN



ENTRY INTRACELLULAR_TOXIN (J)
IN_TOXIN(KT:KMX-1,IU:ID,J)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL 
       IF(ALG_CALC(JA))THEN
       IN_TOXIN(K,I,J)=IN_TOXIN(K,I,J) + CTP(J,JA)*CTB(J,JA)*ALG(K,I,JA)                               !CTISS(K,I,J) = CTISS(K,I,J) + CTP(J,JA)*CTB(J,JA)*(AGR(K,I,JA)-AER(K,I,JA)-AMR(K,I,JA)-ARR(K,I,JA))*ALG(K,I,JA) 
       END IF
      END DO
      !CTISS(K,I,J) = CTISS(K,I,J) + (-CTL(J)-CTA(J)-CTDI(J))*INTOXIN(K,I,J) 
    END DO
  END DO
RETURN

        
ENTRY EXTRACELLULAR_TOXIN (J)
CTESS(KT:KMX-1,IU:ID,J)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        IF(ALG_CALC(JA))THEN
        CTESS(K,I,J) = CTESS(K,I,J) + CTP(J,JA)*CTB(J,JA)*AMR(K,I,JA)*ALG(K,I,JA) 
        END IF
      END DO
      CTESS(K,I,J) = CTESS(K,I,J) + CTREL(J)*IN_TOXIN(K,I,J)-CTD(J)*EX_TOXIN(K,I,J) 
    END DO
  END DO
RETURN


!***********************************************************************************************************************************
!**                                            B I O C H E M I C A L   O 2   D E M A N D                                          **
!***********************************************************************************************************************************

ENTRY BIOCHEMICAL_O2_DEMAND(JBOD)
  IF(JBOD == 1)CBODNS(:,IU:ID) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      CBODSET = CBODS(JBOD)*(CBOD(K-1,I,JBOD)-CBOD(K,I,JBOD))*BI(K,I)/BH2(K,I)
      CBODNS(K,I)=CBODNS(K,I)+CBODSET
      CBODSS(K,I,JBOD) = -CBODD(K,I,JBOD)*CBOD(K,I,JBOD)+CBODSET
    END DO
  END DO
RETURN

! VARIABLE STOCHIOMETRY FOR CBOD SECTION ! CB 6/6/10
!***********************************************************************************************************************************
!**                                            B I O C H E M I C A L   O 2   D E M A N D   P H O S P H O R U S                    **
!***********************************************************************************************************************************

ENTRY BIOCHEMICAL_O2_DEMAND_P(JBOD)
  IF(JBOD == 1)CBODNSP(:,IU:ID) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      CBODSET = CBODS(JBOD)*(CBODP(K-1,I,JBOD)-CBODP(K,I,JBOD))*BI(K,I)/BH2(K,I)
      CBODNSP(K,I)=CBODNSP(K,I)+CBODSET
      CBODPSS(K,I,JBOD) = -CBODD(K,I,JBOD)*CBODP(K,I,JBOD)+CBODSET
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                            B I O C H E M I C A L   O 2   D E M A N D   N I T R O G E N                        **
!***********************************************************************************************************************************

ENTRY BIOCHEMICAL_O2_DEMAND_N(JBOD)
  IF(JBOD == 1)CBODNSN(:,IU:ID) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      CBODSET = CBODS(JBOD)*(CBODN(K-1,I,JBOD)-CBODN(K,I,JBOD))*BI(K,I)/BH2(K,I)
      CBODNSN(K,I)=CBODNSN(K,I)+CBODSET
      CBODNSS(K,I,JBOD) = -CBODD(K,I,JBOD)*CBODN(K,I,JBOD)+CBODSET
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                D I S S O L V E D   O X Y G E N                                                **
!***********************************************************************************************************************************

ENTRY DISSOLVED_OXYGEN
  DOAP(:,IU:ID) = 0.0; DOAR(:,IU:ID) = 0.0; DOEP(:,IU:ID) = 0.0; DOER(:,IU:ID) = 0.0; DOBOD(:,IU:ID) = 0.0
  DOMP(:,IU:ID) = 0.0; DOMR(:,IU:ID) = 0.0; DOZR(:,IU:ID)=0.0    
  DOH2S(:,IU:ID)= 0.0; DOCH4(:,IU:ID)= 0.0; DOFE2(:,IU:ID) = 0.0; DOMN2(:,IU:ID) = 0.0 
  
  DO I=IU,ID
    DOSS(KT,I) = 0.0
    DO K=KT,KB(I)
      DO JCB=1,NBOD
        IF(BOD_CALC(JCB))DOBOD(K,I) = DOBOD(K,I)+RBOD(JCB)*CBODD(K,I,JCB)*CBOD(K,I,JCB)
      END DO
      DO JA=1,NAL
      IF(ALG_CALC(JA))THEN
        DOAP(K,I) = DOAP(K,I)+AGR(K,I,JA)*ALG(K,I,JA)*O2AG(JA)
        DOAR(K,I) = DOAR(K,I)+ARR(K,I,JA)*ALG(K,I,JA)*O2AR(JA)
      ENDIF
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))THEN
        DOEP(K,I) = DOEP(K,I)+EGR(K,I,JE)*EPC(K,I,JE)*O2EG(JE)
        DOER(K,I) = DOER(K,I)+ERR(K,I,JE)*EPC(K,I,JE)*O2ER(JE)
        ENDIF
      END DO

      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          IF(K.EQ.KT)THEN
            JT=KTI(I)
          ELSE
            JT=K
          END IF
          JE=KB(I)
          DO JJ=JT,JE
            DOMP(K,I)=DOMP(K,I)+MGR(JJ,K,I,M)*MACRM(JJ,K,I,M)*O2MG(M)
            DOMR(K,I)=DOMR(K,I)+MRR(K,I,M)*MACRM(JJ,K,I,M)*O2MR(M)
          END DO
        END IF
      END DO
      DOMP(K,I)=DOMP(K,I)/(DLX(I)*BH2(K,I))
      DOMR(K,I)=DOMR(K,I)/(DLX(I)*BH2(K,I))
      DOPOM(K,I) = (LPOMD(K,I)+RPOMD(K,I))*O2OM(JW)
      DODOM(K,I) = (LDOMD(K,I)+RDOMD(K,I))*O2OM(JW)
      DOOM(K,I)  =  DOPOM(K,I)+DODOM(K,I)+DOBOD(K,I)      
      DONIT(K,I) =  NH4D(K,I)*O2NH4(JW)
            IF(STANDING_BIOMASS_DECAY)THEN  ! SW 5/26/15
            DOSED(K,I) =  SEDD(K,I)*O2OM(JW) +SEDD1(K,I)*O2OM(JW)+SEDD2(K,I)*O2OM(JW)   !Amaila
            ELSE
            DOSED(K,I) =  SEDD(K,I)*O2OM(JW)
            ENDIF
      DOSOD(K,I) =  SODD(K,I)*DO3(K,I)
      
      DOCH4(K,I) = CH4D(K,I)*O2CH4      
      DOH2S(K,I) = H2SD(K,I)*O2H2S
      DOFE2(K,I) = FE2D(K,I)*O2FE2
      DOMN2(K,I) = MN2D(K,I)*O2MN2

     IF(ZOOPLANKTON_CALC)THEN
     DO JZ = 1, NZP
      DOZR(K,I)  = DOZR(K,I)+ZRT(K,I,JZ)*ZOO(K,I,JZ)*O2ZR(JZ)
	 END DO
	 ENDIF
    DOSS(K,I)  =  DOAP(K,I)+DOEP(K,I)-DOAR(K,I)-DOER(K,I)-DOOM(K,I)-DONIT(K,I)-DOSOD(K,I)-DOSED(K,I)+DOMP(K,I)-DOMR(K,I)-DOZR(K,I)+doch4(k,i)+doh2s(k,i)+dofe2(k,i)+doMn2(k,i)     ! doch4, doh2s,dofe2 already negative...
    END DO
    IF (.NOT. ICE(I)) THEN
      DOSAT = SATO(T1(KT,I),TDS(KT,I),PALT(I),SALT_WATER(JW))
      O2EX       =  REAER(I)                     ! in s-1
      DOAE(KT,I) = (DOSAT-O2(KT,I))*O2EX         !*BI(KT,I)/BH2(KT,I)  THIS IS NOW DONE IN gas-transfer.f90 - used to convert from m/d to 1/d then to 1/s
      DOSS(KT,I) =  DOSS(KT,I)+DOAE(KT,I)
    END IF
  END DO
RETURN

!***********************************************************************************************************************************
!**                                              I N O R G A N I C   C A R B O N                                                  **
!***********************************************************************************************************************************

ENTRY INORGANIC_CARBON
  TICAP(:,IU:ID) = 0.0; TICEP(:,IU:ID) = 0.0; TICBOD(:,IU:ID) = 0.0
  ticmc(:,iu:id) = 0.0; ticzr(:,iu:id)=0.0  !v3.5
  ticch4=0.0  ! CEMA
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JCB=1,NBOD
        IF(BOD_CALC(JCB))TICBOD(K,I) = TICBOD(K,I)+CBODD(K,I,JCB)*CBOD(K,I,JCB)*BODC(JCB)
      END DO
      DO JA=1,NAL
        IF(ALG_CALC(JA))TICAP(K,I) = TICAP(K,I)+AC(JA)*(ARR(K,I,JA)-AGR(K,I,JA))*ALG(K,I,JA)
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))TICEP(K,I) = TICEP(K,I)+EC(JE)*(ERR(K,I,JE)-EGR(K,I,JE))*EPC(K,I,JE)
      END DO
      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          IF(K.EQ.KT)THEN
            JT=KTI(I)
          ELSE
            JT=K
          END IF
          JE=KB(I)
          DO JJ=JT,JE
            TICMC(K,I)=TICMC(K,I)+(MRR(K,I,M)-MGR(JJ,K,I,M))*MACRM(JJ,K,I,M)*MC(M)
          END DO
        END IF
      END DO
      TICMC(K,I) = TICMC(K,I)/(DLX(I)*BH2(K,I))              !8/2020
      IF(ZOOPLANKTON_CALC)THEN
      DO JZ = 1,NZP
        TICZR(K,I)=TICZR(K,I)+ZRT(K,I,JZ)*ZOO(K,I,JZ)*ZC(JZ) !MLM
      END DO
      
      ticCH4=ch4d(k,i)
      
      ENDIF

            IF(STANDING_BIOMASS_DECAY)THEN  ! SW 5/26/15
               TICSS(K,I) = TICAP(K,I)+TICEP(K,I)+SEDDC(K,I)+ORGC(JW)*(LPOMD(K,I)+RPOMD(K,I)+LDOMD(K,I)+RDOMD(K,I))                          &
                   +CO2R(JW)*SODD(K,I)*DO3(K,I)+TICBOD(K,I)+TICMC(K,I)+TICZR(K,I)  + ticch4               &
                   +sedd1(k,i)*cbiom(jw) + sedd2(k,i)*cbiom(jw)             ! Amaila, cb 6/7/17

            ELSE
                TICSS(K,I) = TICAP(K,I)+TICEP(K,I)+SEDDC(K,I)+ORGC(JW)*(LPOMD(K,I)+RPOMD(K,I)+LDOMD(K,I)+RDOMD(K,I))                          &
              +CO2R(JW)*SODD(K,I)*DO3(K,I)+TICBOD(K,I)+TICMC(K,I)+TICZR(K,I) + ticch4      !!8/2020   
            ENDIF
      
    END DO
    IF (.NOT. ICE(I)) THEN
     ! IF (REAER(I) == 0.0) CALL GAS_TRANSFER
      CO2EX       = REAER(I)*0.923
      KHCO2=12000.*10**(2385.73/(T2(KT,I)+273.15)-14.0184+0.0152642*(T2(KT,I)+273.15))   ! KH CO2 in mg/l/atm as C    SW 8/16/2020

     ! CO2REAER(KT,I)=CO2EX*(0.286*EXP(-0.0314*(T2(KT,I))*PALT(I))-CO2(KT,I))*BI(KT,I)/BH2(KT,I)     ! SW 8/16/2020 
      
      CO2REAER(KT,I)=CO2EX*(PCO2*KHCO2-CO2(KT,I))         !This computation is done in gas-transfer.f90: *BI(KT,I)/BH2(KT,I)   REAER in units of s-1
      
      TICSS(KT,I) = TICSS(KT,I)+CO2REAER(KT,I)
    END IF
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      S E D I M E N T                                                          **
!***********************************************************************************************************************************

ENTRY SEDIMENT
  SEDAS(:,IU:ID) = 0.0; LPOMEP(:,IU:ID) = 0.0; SEDCB(:,IU:ID) = 0.0
  DO I=IU,ID
    SEDSI=0.0
    DO K=KT,KB(I)
    IF(K == KB(I))THEN
    BIBH2(K,I)=BI(K,I)/BH2(K,I)
    ELSE
    BIBH2(K,I)=BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
    ENDIF
      DO JA=1,NAL
        IF(ALG_CALC(JA))SEDAS(K,I) = SEDAS(K,I)+MAX(AS(JA),0.0)*ALG(K,I,JA)*BIBH2(K,I)                !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      END DO
      SEDEM = 0.0   ! CB 5/19/06
      DO JE=1,NEP
!        LPOMEP(K,I) = LPOMEP(K,I)+EPOM(JE)*(EMR(K,I,JE)*EPC(K,I,JE))
        IF (EPIPHYTON_CALC(JW,JE))SEDEM = SEDEM+EBR(K,I,JE)*EPC(K,I,JE)    ! SW 3/2019 SEDEM = SEDEM+EBR(K,I,JE)/H1(K,I)*EPC(K,I,JE)    ! cb 5/19/06     
      END DO
      DO JD=1,NBOD
        IF(BOD_CALC(JD))SEDCB(K,I) = SEDCB(K,I)+MAX(CBODS(JD),0.0)*CBOD(K,I,JD)*BIBH2(K,I)/O2OM(JW)           !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      END DO
      !
      IF(ORGC_CALC)THEN
        SEDOMS(K,I) = POMS(JW)*LPOMC(K,I)/ORGC(JW)*BIBH2(K,I)+POMS(JW)*RPOMC(K,I)/ORGC(JW)*BIBH2(K,I) 
      ELSE
        SEDOMS(K,I) = POMS(JW)*LPOM(K,I)*BIBH2(K,I)+POMS(JW)*RPOM(K,I)*BIBH2(K,I)
      END IF  
      IF(K==KB(I))THEN
      SEDSO       = 0.0
      ELSE
      SEDSO       = SEDS(JW)*SED(K,I)*BI(K+1,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      ENDIF
      SEDNS(K,I)  = SEDSI-SEDSO
      SEDSI       = SEDSO
      if(k < kb(i))then   ! CEMA sediment in kb layer goes to sediment diagenesis model
        SED(K,I)    = MAX(SED(K,I)+(SEDEM+SEDAS(K,I)+SEDCB(K,I)+SEDOMS(K,I)+SEDNS(K,I)-SEDD(K,I)-SEDBR(K,I))*DLT,0.0)   ! cb 11/30/06
      else if(k == kb(i))then
        SED(K,I)    = MAX(SED(K,I)+(SEDEM+SEDAS(K,I)+SEDCB(K,I)+SEDOMS(K,I)+SEDNS(K,I)-SEDD(K,I)-SEDBR(K,I))*DLT,0.0)
      end if
    END DO
  END DO
RETURN


!***********************************************************************************************************************************
!**                                                      S E D I M E N T   P H O S P H O R U S                                    **
!***********************************************************************************************************************************

ENTRY SEDIMENTP
  SEDASP(:,IU:ID) = 0.0; LPOMEPP(:,IU:ID) = 0.0; SEDCBP(:,IU:ID) = 0.0
  DO I=IU,ID
    SEDSIP=0.0
    DO K=KT,KB(I)
      DO JA=1,NAL
        IF(ALG_CALC(JA))then 
          SEDASP(K,I) = SEDASP(K,I)+MAX(AS(JA),0.0)*AP(JA)*ALG(K,I,JA)*BIBH2(K,I)          !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
        end if
      END DO
      SEDEM = 0.0   
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))then
          !LPOMEPP(K,I) = LPOMEPP(K,I)+EPOM(JE)*EP(JE)*(EMR(K,I,JE)*EPC(K,I,JE))
          SEDEM = SEDEM+EBR(K,I,JE)*EPC(K,I,JE)*EP(JE)   ! SW 3/2019  
        end if
      END DO
      DO JD=1,NBOD
!        IF(BOD_CALC(JD))SEDCBP(K,I)=SEDCBP(K,I)+MAX(CBODS(JD),0.0)*BODP(JD)*CBOD(K,I,JD)*BIBH2(K,I)      !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
        !IF(BOD_CALC(JD))then 
        !    SEDCBP(K,I)=SEDCBP(K,I)+MAX(CBODS(JD),0.0)*CBODP(K,I,JD)*BIBH2(K,I)    ! CB 6/6/10
        !end if
        IF(BOD_CALCP(JD))THEN
          SEDCBP(K,I) = SEDCBP(K,I) + MAX(CBODS(JD), 0.0)*CBODP(K,I,JD)*RBOD(JD)*BIBH2(K,I)
        ELSE
          SEDCBP(K,I) = SEDCBP(K,I) + MAX(CBODS(JD), 0.0)*CBOD(K,I,JD)*RBOD(JD)*BODP(JD)*BIBH2(K,I)
        END IF  
      END DO
      !SEDOMSP(K,I) = POMS(JW)*(LPOMP(K,I)+RPOMP(K,I))*BIBH2(K,I)
      SEDOMSP(K,I) = POMS(JW)*(LPOP(K,I)+RPOP(K,I))*BIBH2(K,I)
      !
      IF(K == KB(I))THEN
      SEDSOP       = 0.0
      ELSE
      SEDSOP       = SEDS(JW)*SEDP(K,I)*BI(K+1,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      ENDIF
      SEDNSP(K,I)  = SEDSIP-SEDSOP
      SEDSIP       = SEDSOP
! CEMA start
      SEDPINFLUX(K,I)=(SEDEM+SEDASP(K,I)+SEDOMSP(K,I)+SEDCBP(K,I))*DLT             !(LPOMEPP(K,I)+SEDASP(K,I)+SEDOMSP(K,I)+SEDCBP(K,I))*DLT
      if(k < kb(i))then
        !SEDP(K,I)    = MAX(SEDP(K,I)+(LPOMEPP(K,I)+SEDASP(K,I)+SEDOMSP(K,I)+SEDCBP(K,I)+SEDNSP(K,I)-SEDDP(K,I)   &   ! SW 4/8/16
         SEDP(K,I)    = MAX(SEDP(K,I)+SEDPINFLUX(K,I)+(SEDNSP(K,I)-SEDDP(K,I)   &
                     -SEDBRP(K,I))*DLT,0.0)                                                                 !cb 11/30/06
      else if(k == kb(i))then
      !SEDP(K,I)    = MAX(SEDP(K,I)+(LPOMEPP(K,I)+SEDASP(K,I)+SEDOMSP(K,I)+SEDCBP(K,I)+SEDNSP(K,I)-SEDDP(K,I)   &
        SEDP(K,I)    = MAX(SEDP(K,I)+SEDPINFLUX(K,I)+(SEDNSP(K,I)-SEDDP(K,I)   &
                     -SEDBRP(K,I))*DLT,0.0)                                                                 !cb 11/30/06
      end if
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      S E D I M E N T   N I T R O G E N                                        **
!***********************************************************************************************************************************

ENTRY SEDIMENTN
  SEDASN(:,IU:ID) = 0.0; LPOMEPN(:,IU:ID) = 0.0; SEDCBN(:,IU:ID) = 0.0
  DO I=IU,ID
    SEDSIN=0.0
    DO K=KT,KB(I)
      DO JA=1,NAL
        IF(ALG_CALC(JA))then 
            SEDASN(K,I) = SEDASN(K,I)+MAX(AS(JA),0.0)*AN(JA)*ALG(K,I,JA)*BIBH2(K,I)            !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))   
        end if
      END DO
      SEDEM=0.0
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))then
          !LPOMEPN(K,I) = LPOMEPN(K,I)+EPOM(JE)*EN(JE)*(EMR(K,I,JE)*EPC(K,I,JE))
            SEDEM = SEDEM+EBR(K,I,JE)*EPC(K,I,JE)*EN(JE)   ! SW 3/2019
        end if
      END DO
      DO JD=1,NBOD
!        IF(BOD_CALC(JD))SEDCBN(K,I)=SEDCBN(K,I)+MAX(CBODS(JD),0.0)*BODN(JD)*CBOD(K,I,JD)*BIBH2(K,I)        !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
        !IF(BOD_CALC(JD))then  
        !    SEDCBN(K,I)=SEDCBN(K,I)+MAX(CBODS(JD),0.0)*CBODN(K,I,JD)*BIBH2(K,I)    ! CB 6/6/10
        !end if
        IF(BOD_CALCN(JD))then
          SEDCBN(K,I) = SEDCBN(K,I) + MAX(CBODS(JD), 0.0)*CBODN(K,I,JD)*RBOD(JD)*BIBH2(K,I)    
        ELSE
          SEDCBN(K,I) = SEDCBN(K,I) + MAX(CBODS(JD), 0.0)*CBOD(K,I,JD)*BODN(JD)*RBOD(JD)*BIBH2(K,I)     
        END IF  
      END DO 
      !SEDOMSN(K,I) = POMS(JW)*(LPOMN(K,I)+RPOMN(K,I))*BIBH2(K,I)                           !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))  !CB 10/22/06 
      SEDOMSN(K,I) = POMS(JW)*(LPON(K,I) + RPON(K,I))*BIBH2(K,I)
      IF(K == KB(I)) THEN      ! SW 12/16/07
      SEDSON       = 0.0
	  ELSE
      SEDNO3(K,I)  = FNO3SED(JW)*NO3(K,I)*NO3S(JW)*NO3TRM(K,I)*(BI(K,I)-BI(K+1,I))/BH2(K,I)
      SEDSON       = SEDS(JW)*SEDN(K,I)*BI(K+1,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
	  ENDIF
      SEDNSN(K,I)  = SEDSIN-SEDSON
      SEDSIN       = SEDSON
! CEMA start
      SEDNINFLUX(K,I)=(SEDEM+SEDASN(K,I)+SEDOMSN(K,I)+SEDCBN(K,I))*DLT                 !SW 3/2019 (LPOMEPN(K,I)+SEDASN(K,I)+SEDOMSN(K,I)+SEDCBN(K,I))*DLT
      if(k < kb(i))then      
      SEDN(K,I)    = MAX(SEDN(K,I)+SEDNINFLUX(K,I)+(SEDNSN(K,I)+SEDNO3(K,I)   &          
                     -SEDDN(K,I)-SEDBRN(K,I))*DLT,0.0)  !CB 11/30/06                    
      else if(k == kb(i))then
        SEDN(K,I)    = MAX(SEDN(K,I)+SEDNINFLUX(K,I)+(SEDNSN(K,I)+SEDNO3(K,I)   &
                     -SEDDN(K,I)-SEDBRN(K,I))*DLT,0.0)  !CB 11/30/06
      end if
! CEMA end
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      S E D I M E N T   C A R B O N                                            **
!***********************************************************************************************************************************

ENTRY SEDIMENTC
  SEDASC(:,IU:ID) = 0.0; LPOMEPC(:,IU:ID) = 0.0; SEDCBC(:,IU:ID) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        IF(ALG_CALC(JA))then
          SEDASC(K,I) = SEDASC(K,I)+MAX(AS(JA),0.0)*AC(JA)*ALG(K,I,JA)*BIBH2(K,I)             !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
        end if
      END DO
      SEDEM=0.0
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))then
         ! LPOMEPC(K,I) = LPOMEPC(K,I)+EPOM(JE)*EC(JE)*(EMR(K,I,JE)*EPC(K,I,JE))
            SEDEM=SEDEM+EBR(K,I,JE)*EPC(K,I,JE)*EC(JE)   ! SW 3/2019  
        end if
      END DO
      DO JD=1,NBOD
        IF(BOD_CALC(JD))then
          SEDCBC(K,I)=SEDCBC(K,I)+MAX(CBODS(JD),0.0)*BODC(JD)*CBOD(K,I,JD)*BIBH2(K,I)         !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
        end if
      END DO
      !SEDOMSC(K,I) = POMS(JW)*ORGC(JW)*(LPOM(K,I)+RPOM(K,I))*BIBH2(K,I)                     !BI(K,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))   !CB 10/22/06
      SEDOMSC(K,I) = POMS(JW)*(LPOC(K,I) + RPOC(K,I))*BIBH2(K,I)
      IF(K == KB(I))THEN
      SEDSOC       = 0.0
      ELSE
      SEDSOC       = SEDS(JW)*SEDC(K,I)*BI(K+1,I)/BH2(K,I)*(1.0-BI(K+1,I)/BI(K,I))
      ENDIF
      SEDNSC(K,I)  = SEDSIC-SEDSOC
      SEDSIC       = SEDSOC
! CEMA start
      if(k < kb(i))then
      SEDC(K,I)    = MAX(SEDC(K,I)+(SEDEM+SEDASC(K,I)+SEDOMSC(K,I)+SEDCBC(K,I)+SEDNSC(K,I)-SEDDC(K,I)    &            !(LPOMEPC(K,I)+SEDASC(K,I)+SEDOMSC(K,I)+SEDCBC(K,I)+SEDNSC(K,I)-SEDDC(K,I)    &
                     -SEDBRC(K,I))*DLT,0.0)           
      else if(k == kb(i) .and. .not. sediment_diagenesis)then
        SEDC(K,I)    = MAX(SEDC(K,I)+(SEDEM+SEDASC(K,I)+SEDOMSC(K,I)+SEDCBC(K,I)+SEDNSC(K,I)-SEDDC(K,I)    &   ! (LPOMEPC(K,I)+SEDASC(K,I)+SEDOMSC(K,I)+SEDCBC(K,I)+SEDNSC(K,I)-SEDDC(K,I)    &
                     -SEDBRC(K,I))*DLT,0.0)
      end if
! CEMA end
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      S E D I M E N T   D E C A Y    R A T E                                   **
!***********************************************************************************************************************************

ENTRY SEDIMENT_DECAY_RATE
  DO I=IU,ID
    SEDSIDK=0.0
    DO K=KT,KB(I)
      SEDSUM=0.0
      SEDSUMK=0.0
      
      DO JA=1,NAL
        IF(ALG_CALC(JA))THEN
        XDUM=MAX(AS(JA),0.0)*ALG(K,I,JA)*BIBH2(K,I)
        IF(ORGC_CALC)THEN
            SEDSUMK = SEDSUMK + XDUM * LPOMCDK(JW)
        ELSE
            SEDSUMK = SEDSUMK + XDUM * LPOMDK(JW)
        END IF
        SEDSUM  = SEDSUM  + XDUM
        ENDIF
      END DO
      
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))THEN
        XDUM=EPOM(JE)*(EMR(K,I,JE)*EPC(K,I,JE))
            IF(ORGC_CALC)THEN
            SEDSUMK = SEDSUMK + XDUM * LPOMCDK(JW)
            ELSE
            SEDSUMK = SEDSUMK + XDUM * LPOMDK(JW)
            ENDIF
        SEDSUM  = SEDSUM  + XDUM
        ENDIF
      END DO
      
      DO JD=1,NBOD
        IF(BOD_CALC(JD))THEN
        XDUM=MAX(CBODS(JD),0.0)*CBOD(K,I,JD)*BIBH2(K,I)*RBOD(JD)/O2OM(JW)
        SEDSUMK = SEDSUMK+XDUM*CBODD(K,I,JD)               
        SEDSUM  = SEDSUM + XDUM
        ENDIF
      END DO
      !
      IF(ORGC_CALC)THEN
          SEDSUMK = SEDSUMK + POMS(JW)*LPOMC(K,I)/ORGC(JW)*LPOMCDK(JW)*BIBH2(K,I)+POMS(JW)*RPOMC(K,I)/ORGC(JW)*RPOMCDK(JW)*BIBH2(K,I)
          SEDSUM  = SEDSUM  + POMS(JW)*LPOMC(K,I)/ORGC(JW)*BIBH2(K,I)+POMS(JW)*RPOMC(K,I)/ORGC(JW)*BIBH2(K,I)      
      ELSE
          SEDSUMK = SEDSUMK + POMS(JW)*LPOM(K,I)*LPOMDK(JW)*BIBH2(K,I)+POMS(JW)*RPOM(K,I)*RPOMDK(JW)*BIBH2(K,I)
          SEDSUM  = SEDSUM  + POMS(JW)*LPOM(K,I)*BIBH2(K,I)+POMS(JW)*RPOM(K,I)*BIBH2(K,I)
      END IF
      
      SEDSUMK = SEDSUMK*DLT
      SEDSUM  = SEDSUM*DLT  
    
      IF((SEDSUM+SED(K,I)) > 0.0)THEN
      SDKV(K,I)    = (SEDSUMK+SED(K,I) * SDKV(K,I))/(SEDSUM+ SED(K,I))
      ELSE
      SDKV(K,I)=0.0
      ENDIF
            
    END DO
  END DO
RETURN

! 
! additional sediment compartments simulate slow and fast decaying OM left in standing trees
!***********************************************************************************************************************************
!**                                                      S E D I M E N T  1                                                       **
!***********************************************************************************************************************************

ENTRY SEDIMENT1  
  DO I=IU,ID    
    DO K=KT,KB(I)    
      SED1(K,I)    = MAX(SED1(K,I)+(-SEDD1(K,I))*DLT,0.0)      
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      S E D I M E N T  2                                                       **
!***********************************************************************************************************************************

ENTRY SEDIMENT2  
  DO I=IU,ID    
    DO K=KT,KB(I)    
      SED2(K,I)    = MAX(SED2(K,I)+(-SEDD2(K,I))*DLT,0.0)      
    END DO
  END DO
RETURN


!***********************************************************************************************************************************
!*                                                         E P I P H Y T O N                                                      **
!***********************************************************************************************************************************

ENTRY EPIPHYTON (J)
  DO I=IU,ID

!** Limiting factor

    LIGHT = (1.0-BETA(JW))*SRON(JW)*SHADE(I)/ESAT(J)
    LAM2  =  LIGHT
    LAM1  =  LIGHT
    DO K=KT,KB(I)

!**** Limiting factor

      LAM1          = LAM2
      LAM2          = LAM1*EXP(-GAMMA(K,I)*H1(K,I))
      FDPO4         = 1.0-FPSS(K,I)-FPFE(K,I)
      ELLIM(K,I,J)  = 2.718282*(EXP(-LAM2)-EXP(-LAM1))/(GAMMA(K,I)*H1(K,I))
      IF (EHSP(J)  /= 0.0) EPLIM(K,I,J) =  FDPO4*PO4(K,I)/(FDPO4*PO4(K,I)+EHSP(J)+NONZERO)
      IF (EHSN(J)  /= 0.0) ENLIM(K,I,J) = (NH4(K,I)+NO3(K,I))/(NH4(K,I)+NO3(K,I)+EHSN(J)+NONZERO)
      IF (EHSSI(J) /= 0.0) ESLIM(K,I,J) =  DSI(K,I)/(DSI(K,I)+EHSSI(J)+NONZERO)
      LIMIT         =  MIN(EPLIM(K,I,J),ENLIM(K,I,J),ESLIM(K,I,J),ELLIM(K,I,J))
      BLIM          =  1.0-(EPD(K,I,J)/(EPD(K,I,J)+EHS(J)))

!**** Sources/sinks

      EGR(K,I,J) =  MIN(ETRM(K,I,J)*EG(J)*LIMIT*BLIM,PO4(K,I)/(EP(J)*DLT*EPD(K,I,J)/H1(KT,I)+NONZERO),(NH4(K,I)+NO3(K,I))/(EN(J)   &
                    *DLT*EPD(K,I,J)/H1(K,I)+NONZERO))
      ERR(K,I,J) =  ETRM(K,I,J)*ER(J)*DO3(K,I)
      EMR(K,I,J) = (ETRMR(K,I,J)+1.0-ETRMF(K,I,J))*EM(J)
      EER(K,I,J) =  MIN((1.0-ELLIM(K,I,J))*EE(J)*ETRM(K,I,J),EGR(K,I,J))
!      EPD(K,I,J) =  MAX(EPD(K,I,J)+EPD(K,I,J)*(EGR(K,I,J)-ERR(K,I,J)-EMR(K,I,J)-EER(K,I,J)-EBR(K,I,J)/(H1(K,I)*0.0025))*DLT,0.0)
!      EPD(K,I,J) =  MAX(EPD(K,I,J)+EPD(K,I,J)*(EGR(K,I,J)-ERR(K,I,J)-EMR(K,I,J)-EER(K,I,J)-EBR(K,I,J)/H1(K,I))*DLT,0.00)   ! cb 5/18/06
      EPD(K,I,J) =  MAX(EPD(K,I,J)+EPD(K,I,J)*(EGR(K,I,J)-ERR(K,I,J)-EMR(K,I,J)-EER(K,I,J)-EBR(K,I,J))*DLT,0.00)   ! SW 3/2019
      if(k == kb(i)) then      ! SW 12/16/07
      EPM(K,I,J) =  EPD(K,I,J)*(BI(K,I)+2.0*H1(K,I))*DLX(I)
	  else
      EPM(K,I,J) =  EPD(K,I,J)*(BI(K,I)-BI(K+1,I)+2.0*H1(K,I))*DLX(I)
	  endif
      EPC(K,I,J) =  EPM(K,I,J)/VOL(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                       L A B I L E   D O M   P H O S P H O R U S                               **
!***********************************************************************************************************************************

ENTRY LABILE_DOM_P
  LDOMPAP(:,IU:ID) = 0.0; LDOMPEP(:,IU:ID) = 0.0; LDOMPMP(:,IU:ID)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        IF(ALG_CALC(JA))LDOMPAP(K,I) = LDOMPAP(K,I)+(AER(K,I,JA)+(1.0-APOM(JA))*AMR(K,I,JA))*ALG(K,I,JA)*AP(JA)
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))LDOMPEP(K,I) = LDOMPEP(K,I)+(EER(K,I,JE)+(1.0-EPOM(JE))*EMR(K,I,JE))*EPC(K,I,JE)*EP(JE)
      END DO
      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          IF(K.EQ.KT)THEN
            JT=KTI(I)
          ELSE
            JT=K
          END IF
          JE=KB(I)
          DO JJ=JT,JE
            LDOMPMP(K,I)=LDOMPMP(K,I)+(1.0-MPOM(M))*MMR(K,I,M)*MACRM(JJ,K,I,M)*MP(M)
          END DO
        END IF
      END DO
      LDOMPMP(K,I) = LDOMPMP(K,I)/(DLX(I)*BH2(K,I))       !8/2020
      LDOMPSS(K,I) = LDOMPAP(K,I)+LDOMPEP(K,I)+LDOMPMP(K,I)- LDOMPD(K,I) - LRDOMPD(K,I) + LPOMPHD(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  R E F R A C T O R Y   D O M   P H O S P H O R U S                            **
!***********************************************************************************************************************************

ENTRY REFRACTORY_DOM_P
  DO I=IU,ID
    DO K=KT,KB(I)
      RDOMPSS(K,I) = LRDOMPD(K,I) - RDOMPD(K,I) + RPOMPHD(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      L A B I L E   P O M   P H O S P H O R U S                                **
!***********************************************************************************************************************************

ENTRY LABILE_POM_P
  LPOMPAP(:,IU:ID) = 0.0;LPOMPMP(:,IU:ID)=0.0;LPZOOINP(:,IU:ID)=0.0; LPZOOOUTP(:,IU:ID)=0.0; LPOMEPP(:,IU:ID) = 0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        IF(ALG_CALC(JA))LPOMPAP(K,I) = LPOMPAP(K,I)+APOM(JA)*(AMR(K,I,JA)*ALG(K,I,JA))*AP(JA)
      END DO
      DO JE=1,NEP                                                          
        IF (EPIPHYTON_CALC(JW,JE))LPOMEPP(K,I) = LPOMEPP(K,I)+EPOM(JE)*EP(JE)*(EMR(K,I,JE)*EPC(K,I,JE))    
      END DO                                                              

      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          JT=K
          JE=KB(I)
          DO JJ=JT,JE
            LPOMPMP(K,I)=LPOMPMP(K,I)+MPOM(M)*LRPMAC(M)*MMR(K,I,M)*MACRM(JJ,K,I,M)*MP(M)
          END DO
        END IF
      END DO
      LPOMPMP(K,I) = LPOMPMP(K,I)/(DLX(I)*BH2(K,I))          !8/2020
	IF(ZOOPLANKTON_CALC)THEN
	DO JZ = 1,NZP
      IF(TGRAZE(K,I,JZ) > 0.0)THEN
        LPZOOOUTP(K,I)=LPZOOOUTP(K,I) + ZOO(K,I,JZ)*(ZMT(K,I,JZ)+(ZMU(K,I,JZ)-(ZMU(K,I,JZ)*ZEFF(JZ))))*ZP(JZ)
        IF(ORGC_CALC)THEN
          LPZOOINP(K,I) = LPZOOINP(K,I) + ZOO(K,I,JZ)*ZMU(K,I,JZ)*PREFP(JZ)*LPOMC(K,I)/ORGC(JW)/TGRAZE(K,I,JZ)*ZP(JZ)
        ELSE    
          LPZOOINP(K,I) = LPZOOINP(K,I) + ZOO(K,I,JZ)*ZMU(K,I,JZ)*PREFP(JZ)*LPOM(K,I)/TGRAZE(K,I,JZ)*ZP(JZ)
        END IF 
      ELSE
        LPZOOOUTP(K,I)=LPZOOOUTP(K,I)+ZOO(K,I,JZ)*(ZMT(K,I,JZ)+(ZMU(K,I,JZ)-(ZMU(K,I,JZ)*ZEFF(JZ))))*ZP(JZ)
        LPZOOINP(K,I)=0.0
      END IF
    END DO
    ENDIF
	  LPOMPNS(K,I) = POMS(JW) * (LPOMP(K-1,I)-LPOMP(K,I)) *BI(K,I)/BH2(K,I)                                                
    LPOMPSS(K,I) = LPOMPAP(K,I) + LPOMPEP(K,I) + LPOMPMP(K,I) - LPOMPD(K,I) + LPOMPNS(K,I) - LRPOMPD(K,I) - LPOMPHD(K,I)    
	  IF(ZOOPLANKTON_CALC)THEN
	!  DO JZ = 1,NZP                                           ! KV 4/24/12
	   LPOMPSS(K,I) =LPOMPSS(K,I) + LPZOOOUTP(K,I)-LPZOOINP(K,I)
	!  END DO                                                  ! KV 4/24/12
	  ENDIF

	END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  R E F R A C T O R Y   P O M   P H O S P H O R U S                            **
!***********************************************************************************************************************************

ENTRY REFRACTORY_POM_P
  RPOMPMP(:,IU:ID)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          JT=K
          JE=KB(I)
          DO JJ=JT,JE
            RPOMPMP(K,I)=RPOMPMP(K,I)+MPOM(M)*(1.0-LRPMAC(M))*MMR(K,I,M)*MACRM(JJ,K,I,M)*MP(M)
          END DO
        END IF
      END DO
      RPOMPMP(K,I) = RPOMPMP(K,I)/(DLX(I)*BH2(K,I))       !8/2020
      RPOMPNS(K,I) = POMS(JW) * (RPOMP(K-1,I) - RPOMP(K,I)) * BI(K,I)/BH2(K,I)     
      RPOMPSS(K,I) = LRPOMPD(K,I) + RPOMPNS(K,I) - RPOMPD(K,I) + RPOMPMP(K,I) - RPOMPHD(K,I)  
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                       L A B I L E   D O M   N I T R O G E N                                   **
!***********************************************************************************************************************************

ENTRY LABILE_DOM_N
  LDOMNAP(:,IU:ID) = 0.0; LDOMNEP(:,IU:ID) = 0.0; LDOMNMP(:,IU:ID)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        IF(ALG_CALC(JA))LDOMNAP(K,I) = LDOMNAP(K,I)+(AER(K,I,JA)+(1.0-APOM(JA))*AMR(K,I,JA))*ALG(K,I,JA)*AN(JA)
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE))LDOMNEP(K,I) = LDOMNEP(K,I)+(EER(K,I,JE)+(1.0-EPOM(JE))*EMR(K,I,JE))*EPC(K,I,JE)*EN(JE)
      END DO
      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          IF(K.EQ.KT)THEN
            JT=KTI(I)
          ELSE
            JT=K
          END IF
          JE=KB(I)
          DO JJ=JT,JE
            LDOMNMP(K,I)=LDOMNMP(K,I)+(1.0-MPOM(M))*MMR(K,I,M)*MACRM(JJ,K,I,M)*MN(M)
          END DO
        END IF
      END DO
      LDOMNMP(K,I) = LDOMNMP(K,I)/(DLX(I)*BH2(K,I))        !8/2020
      LDOMNSS(K,I) = LDOMNAP(K,I) + LDOMNEP(K,I) + LDOMNMP(K,I) - LDOMND(K,I) - LRDOMND(K,I) + LPOMNHD(K,I) 
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  R E F R A C T O R Y   D O M   N I T R O G E N                                **
!***********************************************************************************************************************************

ENTRY REFRACTORY_DOM_N
  DO I=IU,ID
    DO K=KT,KB(I)
      !RDOMNSS(K,I) = LRDOMD(K,I)*ORGNLD(K,I)-RDOMD(K,I)*ORGNRD(K,I)
       RDOMNSS(K,I) = LRDOMND(K,I) - RDOMND(K,I) + RPOMNHD(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      L A B I L E   P O M   N I T R O G E N                                    **
!***********************************************************************************************************************************

ENTRY LABILE_POM_N
  LPOMNAP(:,IU:ID) = 0.0;LPOMNMP(:,IU:ID)=0.0;LPZOOINN(:,IU:ID)=0.0; LPZOOOUTN(:,IU:ID)=0.0; LPOMEPN(:,IU:ID)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        IF(ALG_CALC(JA))LPOMNAP(K,I) = LPOMNAP(K,I)+APOM(JA)*(AMR(K,I,JA)*ALG(K,I,JA))*AN(JA)
      END DO
       DO JE=1,NEP                                                          
        IF (EPIPHYTON_CALC(JW,JE))LPOMEPN(K,I) = LPOMEPN(K,I)+EPOM(JE)*EN(JE)*(EMR(K,I,JE)*EPC(K,I,JE))    
      END DO                                                              

      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          JT=K
          JE=KB(I)
          DO JJ=JT,JE
            LPOMNMP(K,I)=LPOMNMP(K,I)+MPOM(M)*LRPMAC(M)*MMR(K,I,M)*MACRM(JJ,K,I,M)*MN(M)
          END DO
        END IF
      END DO
      LPOMNMP(K,I)=LPOMNMP(K,I)/(DLX(I)*BH2(K,I))
	IF(ZOOPLANKTON_CALC)THEN
	DO JZ = 1,NZP
      IF(TGRAZE(K,I,JZ) > 0.0)THEN
        LPZOOOUTN(K,I)=LPZOOOUTN(K,I)+ZOO(K,I,JZ)*(ZMT(K,I,JZ)+(ZMU(K,I,JZ)-(ZMU(K,I,JZ)*ZEFF(JZ))))*ZN(JZ)
        IF(ORGC_CALC)THEN
            LPZOOINN(K,I)=LPZOOINN(K,I)+ZOO(K,I,JZ)*PREFP(JZ)*ZMU(K,I,JZ)*LPOMC(K,I)/ORGC(JW)/TGRAZE(K,I,JZ)*ZN(JZ)
        ELSE
        LPZOOINN(K,I)=LPZOOINN(K,I)+ZOO(K,I,JZ)*PREFP(JZ)*ZMU(K,I,JZ)*LPOM(K,I)/TGRAZE(K,I,JZ)*ZN(JZ)
        END IF
      ELSE
        LPZOOOUTN(K,I)=LPZOOOUTN(K,I)+ZOO(K,I,JZ)*(ZMT(K,I,JZ)+(ZMU(K,I,JZ)-(ZMU(K,I,JZ)*ZEFF(JZ))))*ZN(JZ)
        LPZOOINN(K,I)=0.0
      END IF
	END DO
	ENDIF
    LPOMNNS(K,I) = POMS(JW) * (LPOMN(K-1,I) - LPOMN(K,I)) * BI(K,I)/BH2(K,I)      
    LPOMNSS(K,I) = LPOMNAP(K,I) + LPOMNEP(K,I) + LPOMNMP(K,I) - LPOMND(K,I) + LPOMNNS(K,I) - LRPOMND(K,I) - LPOMNHD(K,I)   &
                   + LPZOOOUTN(K,I) - LPZOOINN(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  R E F R A C T O R Y   P O M   N I T R O G E N                                **
!***********************************************************************************************************************************

ENTRY REFRACTORY_POM_N
  RPOMNMP(:,IU:ID)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          JT=K
          JE=KB(I)
          DO JJ=JT,JE
            RPOMNMP(K,I)=RPOMNMP(K,I)+MPOM(M)*(1.0-LRPMAC(M))*MMR(K,I,M)*MACRM(JJ,K,I,M)*MN(M)
          END DO
        END IF
      END DO
      RPOMNMP(K,I) = RPOMNMP(K,I)/(DLX(I)*BH2(K,I))
      RPOMNNS(K,I) = POMS(JW) * (RPOMN(K-1,I) - RPOMN(K,I)) * BI(K,I)/BH2(K,I)             
      RPOMNSS(K,I) = LRPOMND(K,I) + RPOMNNS(K,I) - RPOMND(K,I) + RPOMNMP(K,I) - RPOMNHD(K,I)   
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                       L A B I L E   D O C
!***********************************************************************************************************************************
ENTRY LABILE_DOM_C
  LDOMCAP(:,IU:ID) = 0.0; LDOMCEP(:,IU:ID) = 0.0; LDOMCMP(:,IU:ID)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        IF(ALG_CALC(JA)) LDOMCAP(K,I) = LDOMCAP(K,I) + (AER(K,I,JA)+(1.0-APOM(JA))*AMR(K,I,JA))*ALG(K,I,JA)*AC(JA)
      END DO
      DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE)) LDOMCEP(K,I) = LDOMCEP(K,I) + (EER(K,I,JE)+(1.0-EPOM(JE))*EMR(K,I,JE))*EPC(K,I,JE)*EC(JE)
      END DO
      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          IF(K.EQ.KT)THEN
            JT=KTI(I)
          ELSE
            JT=K
          END IF
          JE=KB(I)
          DO JJ=JT,JE
            LDOMCMP(K,I) = LDOMCMP(K,I) + (1.0-MPOM(M))*MMR(K,I,M)*MACRM(JJ,K,I,M)*MC(M)
          END DO
        END IF
      END DO
      LDOMCMP(K,I) = LDOMCMP(K,I) / (DLX(I)*BH2(K,I)) 
      LDOMCSS(K,I) = LDOMCAP(K,I) + LDOMCEP(K,I) + LDOMCMP(K,I) - LDOMCD(K,I) - LRDOMCD(K,I) + LPOMCHD(K,I) 
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  R E F R A C T O R Y   D O C
!***********************************************************************************************************************************
ENTRY REFRACTORY_DOM_C
  DO I=IU,ID
    DO K=KT,KB(I)
      RDOMCSS(K,I) = LRDOMCD(K,I) - RDOMCD(K,I) + RPOMCHD(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                      L A B I L E   P O C
!***********************************************************************************************************************************
ENTRY LABILE_POM_C
  LPOMCAP(:,IU:ID) = 0.0; LPOMCEP(:,IU:ID)=0.0; LPOMCMP(:,IU:ID)=0.0; LPZOOINC(:,IU:ID)=0.0; LPZOOOUTC(:,IU:ID)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO JA=1,NAL
        IF(ALG_CALC(JA)) LPOMCAP(K,I) = LPOMCAP(K,I) + APOM(JA)*(AMR(K,I,JA)*ALG(K,I,JA))*AC(JA)
      END DO
      DO JE=1,NEP
        IF(EPIPHYTON_CALC(JW,JE)) LPOMCEP(K,I) = LPOMCEP(K,I) + EPOM(JE)*EMR(K,I,JE)*EPC(K,I,JE)*EC(JE)
      END DO
      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          IF(K.EQ.KT)THEN
            JT=KTI(I)
          ELSE
            JT=K
          END IF
          JE=KB(I)
          DO JJ=JT,JE
            LPOMCMP(K,I) = LPOMCMP(K,I) + MPOM(M)*LRPMAC(M)*MMR(K,I,M)*MACRM(JJ,K,I,M)*MC(M)
          END DO
        END IF
      END DO
      LPOMCMP(K,I) = LPOMCMP(K,I)/(DLX(I)*BH2(K,I))
	    IF(ZOOPLANKTON_CALC)THEN
	      DO JZ = 1,NZP
          IF(TGRAZE(K,I,JZ) > 0.0)THEN
            LPZOOOUTC(K,I) = LPZOOOUTC(K,I) + ZOO(K,I,JZ)*(ZMT(K,I,JZ)+(ZMU(K,I,JZ)-(ZMU(K,I,JZ)*ZEFF(JZ))))*ZC(JZ)
            IF(CAC(NLPOM) == '      ON')THEN
              LPZOOINC(K,I) = LPZOOINC(K,I) + ZOO(K,I,JZ)*PREFP(JZ)*ZMU(K,I,JZ)*LPOM(K,I)/TGRAZE(K,I,JZ)*ZC(JZ)
            ELSE
              LPZOOINC(K,I) = LPZOOINC(K,I) + ZOO(K,I,JZ)*PREFP(JZ)*ZMU(K,I,JZ)*LPOMC(K,I)/ORGC(JW)/TGRAZE(K,I,JZ)*ZC(JZ)
            END IF
          ELSE
            LPZOOOUTC(K,I) = LPZOOOUTC(K,I) + ZOO(K,I,JZ)*(ZMT(K,I,JZ)+(ZMU(K,I,JZ)-(ZMU(K,I,JZ)*ZEFF(JZ))))*ZC(JZ)
            LPZOOINC(K,I) = 0.0
          END IF
	      END DO
	    ENDIF
      LPOMCNS(K,I) = POMS(JW) * (LPOMC(K-1,I)-LPOMC(K,I)) * BI(K,I)/BH2(K,I)
      LPOMCSS(K,I) = LPOMCAP(K,I) + LPOMCEP(K,I) + LPOMCMP(K,I) - LPOMCD(K,I) + LPOMCNS(K,I) - LRPOMCD(K,I) - LPOMCHD(K,I)   &
                   + LPZOOOUTC(K,I) - LPZOOINC(K,I)
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                  R E F R A C T O R Y   P O C
!***********************************************************************************************************************************
ENTRY REFRACTORY_POM_C
  RPOMCMP(:,IU:ID)=0.0
  DO I=IU,ID
    DO K=KT,KB(I)
      DO M=1,NMC
        IF(MACROPHYTE_CALC(JW,M))THEN
          IF(K.EQ.KT)THEN
            JT=KTI(I)
          ELSE
            JT=K
          END IF
          JE=KB(I)
          DO JJ=JT,JE
            RPOMCMP(K,I) = RPOMCMP(K,I) + MPOM(M)*(1.0-LRPMAC(M))*MMR(K,I,M)*MACRM(JJ,K,I,M)*MC(M)
          END DO
        END IF
      END DO
      RPOMCMP(K,I) = RPOMCMP(K,I) / (DLX(I)*BH2(K,I))
      RPOMCNS(K,I) = POMS(JW) * (RPOMC(K-1,I)-RPOMC(K,I)) * BI(K,I)/BH2(K,I)
      RPOMCSS(K,I) = LRPOMCD(K,I) + RPOMCNS(K,I) - RPOMCD(K,I) + RPOMCMP(K,I) - RPOMCHD(K,I) 
    END DO
  END DO
RETURN

!************************************************************************
!**                          M A C R O P H Y T E                       **
!************************************************************************

ENTRY MACROPHYTE(LLM)
  M=LLM
  DO I=IU,ID
    IF(KTICOL(I))THEN
      JT=KTI(I)
    ELSE
      JT=KTI(I)+1
    END IF
    JE=KB(I)
    DO JJ=JT,JE
      IF(JJ.LT.KT)THEN
        COLB=EL(JJ+1,I)
      ELSE
        COLB=EL(KT+1,I)
      END IF
      !COLDEP=ELWS(I)-COLB
       coldep=EL(KT,i)-Z(i)*COSA(JB)-colb  ! cb 3/7/16
      IF(MACRC(JJ,KT,I,M).GT.MMAX(M))THEN
        MGR(JJ,KT,I,M)=0.0
      END IF
      MACSS(JJ,KT,I,M) = (MGR(JJ,KT,I,M)-MMR(KT,I,M)-MRR(KT,I,M))*MACRC(JJ,KT,I,M)
      MACRM(JJ,KT,I,M)   = MACRM(JJ,KT,I,M)+MACSS(JJ,KT,I,M)*DLT*COLDEP*CW(JJ,I)*DLX(I)
    END DO

    DO K=KT+1,KB(I)
      JT=K
      JE=KB(I)
      DO JJ=JT,JE
        IF(MACRC(JJ,K,I,M).GT.MMAX(M))THEN
          MGR(JJ,K,I,M)=0.0
        END IF
        MACSS(JJ,K,I,M) = (MGR(JJ,K,I,M)-MMR(K,I,M)-MRR(K,I,M))*MACRC(JJ,K,I,M)
        IF(MACT(JJ,K,I).GT.MBMP(M).AND.MACT(JJ,K-1,I).LT.MBMP(M).AND.MACSS(JJ,K,I,M).GT.0.0)THEN
          IF(K-1.EQ.KT)THEN
            BMASS=MACSS(JJ,K,I,M)*DLT*H2(K,I)*CW(JJ,I)*DLX(I)
            MACRM(JJ,K-1,I,M)=MACRM(JJ,K-1,I,M)+BMASS
            COLB=EL(KT+1,I)
            !COLDEP=ELWS(I)-COLB
             coldep=EL(KT,i)-Z(i)*COSA(JB)-colb  ! cb 3/7/16
            MACSS(JJ,K-1,I,M)=BMASS/DLT/(COLDEP*CW(JJ,I)*DLX(I)) + MACSS(JJ,K-1,I,M)
          ELSE
            BMASS=MACSS(JJ,K,I,M)*DLT*H2(K,I)*CW(JJ,I)*DLX(I)
            MACRM(JJ,K-1,I,M)=MACRM(JJ,K-1,I,M)+BMASS
            MACSS(JJ,K-1,I,M)=BMASS/DLT/(H2(K-1,I)*CW(JJ,I)*DLX(I))+ MACSS(JJ,K-1,I,M)
          END IF
          MACSS(JJ,K,I,M)=0.0
        ELSE
          BMASSTEST=MACRM(JJ,K,I,M)+MACSS(JJ,K,I,M)*DLT*H2(K,I)*CW(JJ,I)*DLX(I)
          IF(BMASSTEST.GE.0.0)THEN
            MACRM(JJ,K,I,M)   = BMASSTEST
          ELSE
            MACSS(JJ,K,I,M)=-MACRM(JJ,K,I,M)/DLT/(H2(K,I)*CW(JJ,I)*DLX(I))
            MACRM(JJ,K,I,M)=0.0
          END IF
        END IF
      END DO
    END DO
  END DO
  DO I=IU,ID
    TMAC=0.0
    CVOL=0.0
    IF(KTICOL(I))THEN
      JT=KTI(I)
    ELSE
      JT=KTI(I)+1
    END IF
    JE=KB(I)

    DO JJ=JT,JE
      IF(JJ.LT.KT)THEN
        COLB=EL(JJ+1,I)
      ELSE
        COLB=EL(KT+1,I)
      END IF
      !COLDEP=ELWS(I)-COLB
       coldep=EL(KT,i)-Z(i)*COSA(JB)-colb  ! cb 3/7/16
      IF(CW(JJ,I).GT.0.0)THEN
        MACRC(JJ,KT,I,M)=MACRM(JJ,KT,I,M)/(CW(JJ,I)*COLDEP*DLX(I))
      ELSE
        MACRC(JJ,KT,I,M)=0.0
      END IF
      TMAC=TMAC+MACRM(JJ,KT,I,M)
      CVOL=CVOL+CW(JJ,I)*COLDEP*DLX(I)
    END DO

    MAC(KT,I,M)=TMAC/CVOL

    DO K=KT+1,KB(I)
      JT=K
      JE=KB(I)
      TMAC=0.0
      CVOL=0.0
      DO JJ=JT,JE
        IF(CW(JJ,I).GT.0.0)THEN
          MACRC(JJ,K,I,M)=MACRM(JJ,K,I,M)/(CW(JJ,I)*H2(K,I)*DLX(I))
        ELSE
          MACRC(JJ,K,I,M)=0.0
        END IF
        TMAC=TMAC+MACRM(JJ,K,I,M)
        CVOL=CVOL+CW(JJ,I)*H2(K,I)*DLX(I)
      END DO
      MAC(K,I,M)=TMAC/CVOL
    END DO
  END DO

  DO I=IU,ID
    TMAC=0.0
    CVOL=0.0
    DO K=KT,KB(I)
      IF(K.EQ.KT)THEN
        JT=KTI(I)
      ELSE
        JT=K
      END IF
      JE=KB(I)
      DO JJ=JT,JE
        MACT(JJ,K,I)=0.0
        DO MI=1,NMC
          IF(MACROPHYTE_CALC(JW,MI))THEN
            MACT(JJ,K,I)=MACRC(JJ,K,I,MI)+MACT(JJ,K,I)
          END IF
        END DO
      END DO
    END DO
  END DO
  RETURN

!***********************************************************************************************************************************
!*                                                  K I N E T I C   F L U X E S                                                   **
!***********************************************************************************************************************************

ENTRY KINETIC_FLUXES
  DO JAF=1,NAF(JW)
    DO JB=BS(JW),BE(JW)                ! SW 3/9/16
    DO I=CUS(JB),DS(JB)
      DO K=KT,KB(I)
        KFS(K,I,KFCN(JAF,JW)) = KFS(K,I,KFCN(JAF,JW))+KF(K,I,KFCN(JAF,JW))*VOL(K,I)*DLT      ! KF IN G/M3/S x VOL M3 x DT S == G
      END DO
    END DO
    END DO
  ENDDO

  IF(NPBALC=='      ON') THEN      
    DLT13=DLT*0.001     !/1000.
    DO JB=BS(JW),BE(JW)                ! SW 3/9/16
    DO I=CUS(JB),DS(JB)
      DO K=KT,KB(I)
            TN_SEDSOD_NH4(JW)= TN_SEDSOD_NH4(JW)+KF(K,I,KF_NH4_SR)*VOL(K,I)*DLT13+KF(K,I,KF_NH4_SD)*VOL(K,I)*DLT13
            TP_SEDSOD_PO4(JW)= TP_SEDSOD_PO4(JW)+KF(K,I,KF_PO4_SR)*VOL(K,I)*DLT13+KF(K,I,KF_PO4_SD)*VOL(K,I)*DLT13 
            IF(K==KT.AND.CDWBC(NH3_DER,JW)=='      ON')NH3GASLOSS(JW)=NH3GASLOSS(JW)+NH3GAS(K,I)*VOL(K,I)*DLT13    ! NH3GAS in g/m3/s, convert to kg, cumulative      
      END DO
    END DO
    END DO
 ENDIF

RETURN

!***********************************************************************************************************************************
!**                                                       p H   C O 2                                                             **
!***********************************************************************************************************************************

ENTRY PH_CO2

! pH and carbonate species

  DO I=IU,ID
    DO K=KT,KB(I)
      CART = TIC(K,I)/12000.0                ! CART=equivalents/liter of C    TIC=mg/l C (MW=12g/mole)
      ALKT = ALK(K,I)/5.0E+04                ! ALK=mg/l as CaCO3 (MW=50 g/mole; EQ=50g/eq))      ALKT=equivalents/l
      T1K  = T1(K,I)+273.15

!**** Ionic strength

      IF (FRESH_WATER(JW)) S2 = 2.5E-05*TDS(K,I)
      IF (SALT_WATER(JW))  S2 = 1.47E-3+1.9885E-2*TDS(K,I)+3.8E-5*TDS(K,I)*TDS(K,I)

!**** Debye-Huckel terms and activity coefficients

      SQRS2  =  SQRT(S2)
      DH1    = -0.5085*SQRS2/(1.0+1.3124*SQRS2)+4.745694E-03+4.160762E-02*S2-9.284843E-03*S2*S2
      DH2    = -2.0340*SQRS2/(1.0+1.4765*SQRS2)+1.205665E-02+9.715745E-02*S2-2.067746E-02*S2*S2
      H2CO3T =  10.0**(0.0755*S2)
      HCO3T  =  10.0**DH1
      CO3T   =  10.0**DH2
      OH     =  HCO3T

!**** Temperature adjustment

      KW = 10.0**(-283.971-0.05069842*T1K+13323.0/T1K+102.24447*LOG10(T1K)-1119669.0/(T1K*T1K))/OH
      K1 = 10.0**(-3404.71/T1K+14.8435-0.032786*T1K)*H2CO3T/HCO3T
      K2 = 10.0**(-2902.39/T1K+ 6.4980-0.023790*T1K)*HCO3T/CO3T

!**** pH evaluation

      PHT = -PH(K,I)-2.1
      IF (PH(K,I) <= 0.0) PHT = -14.0
      INCR = 10.0
      DO N=1,3
        F    = 1.0
        INCR = INCR/10.0
        ITER = 0
        DO WHILE (F > 0.0 .AND. ITER < 12)
          PHT    = PHT+INCR
          HION   = 10.0**PHT
          BICART = CART*K1*HION/(K1*HION+K1*K2+HION*HION)
          F      = BICART*(HION+2.0*K2)/HION+KW/HION-ALKT-HION/OH
          ITER   = ITER+1
        END DO
        PHT = PHT-INCR
      END DO

!**** pH, carbon dioxide, bicarbonate, and carbonate concentrations

      HION      =  10.0**PHT
      PH(K,I)   = -PHT
      CO2(K,I)  =  TIC(K,I)/(1.0+K1/HION+K1*K2/(HION*HION))          ! mg/l as C
      IF(CDWBC(HCO3_DER,JW)=='      ON')HCO3(K,I) =  TIC(K,I)/(1.0+HION/K1+K2/HION)                    ! mg/l as C
      IF(CDWBC(CO3_DER,JW)=='      ON') CO3(K,I)  =  TIC(K,I)/((HION*HION)/(K1*K2)+HION/K2+1.0)        ! mg/l as C
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                       p H   C O 2   N E W                                                     **
!***********************************************************************************************************************************

ENTRY PH_CO2_NEW ! Enhancements added for buffering by ammonia, phosphate, and OM ! SR 01/01/12
! pH and carbonate species
 DO I=IU,ID
   DO K=KT,KB(I)
     T1K = T1(K,I)+273.15
     CART = TIC(K,I)/12011. ! SR 01/01/12
     ALKT = ALK(K,I)/50044. ! SR 01/01/12
     AMMT = NH4(K,I)/14006.74 ! SR 01/01/12
     PHOST = PO4(K,I)/30973.762 ! SR 01/01/12
     OMCT = (LDOM(K,I)+RDOM(K,I))*ORGC(JW)/12011. ! moles carbon per liter from DOM ! SR 01/01/12
     IF (POM_BUFFERING) OMCT = OMCT + (LPOM(K,I)+RPOM(K,I))*ORGC(JW)/12011. ! SR 01/01/12
 !**** Ionic strength
     IF (FRESH_WATER(JW)) S2 = 2.5E-05*TDS(K,I)
     IF (SALT_WATER(JW)) S2 = 1.47E-3+1.9885E-2*TDS(K,I)+3.8E-5*TDS(K,I)*TDS(K,I)
!**** Debye-Huckel terms and activity coefficients
     SQRS2 = SQRT(S2)
     DH1 = -0.5085*SQRS2/(1.0+1.3124*SQRS2)+4.745694E-03+4.160762E-02*S2-9.284843E-03*S2*S2
     DH2 = -2.0340*SQRS2/(1.0+1.4765*SQRS2)+1.205665E-02+9.715745E-02*S2-2.067746E-02*S2*S2
     DH3 = -4.5765*SQRS2/(1.0+1.3124*SQRS2) ! extended Debye-Huckel for PO4 ! SR 01/01/12
     DHH = -0.5085*SQRS2/(1.0+2.9529*SQRS2) ! extended Debye-Huckel for H+ ion ! SR 01/01/12
     H2CO3T = 10.0**(0.0755*S2)
     HCO3T = 10.0**DH1
     CO3T = 10.0**DH2
     PO4T = 10.0**DH3 ! SR 01/01/12
     HT = 10.0**DHH ! activity coefficient for H+ ! SR 01/01/12
     HPO4T = CO3T ! tabled values similar to those for carbonate ! SR 01/01/12
     OHT = HCO3T ! tabled values similar to those for bicarbonate ! SR 01/01/12
     H2PO4T = HCO3T ! tabled values similar to those for bicarbonate ! SR 01/01/12
     NH4T = HCO3T ! tabled values similar to those for bicarbonate ! SR 01/01/12
     NH3T = H2CO3T ! neutral species, set coefficient to same as that for carbonic acid ! SR 01/01/12
     H3PO4T = H2CO3T ! neutral species, set coefficient to same as that for carbonic acid ! SR 01/01/12
!**** Temperature adjustment
     KW = 10.0**(-283.971 -0.05069842*T1K +13323.0/T1K +102.24447*LOG10(T1K) -1119669.0/(T1K*T1K))/OHT
     K1 = 10.0**(-356.3094 -0.06091964*T1K +21834.37/T1K +126.8339 *LOG10(T1K) -1684915 /(T1K*T1K))*H2CO3T/HCO3T
     K2 = 10.0**(-107.8871 -0.03252849*T1K + 5151.79/T1K + 38.92561*LOG10(T1K) - 563713.9/(T1K*T1K))*HCO3T/CO3T
     KAMM = 10.0**(-0.09018 -2729.92/T1K)*NH4T/NH3T ! SR 01/01/12
     KP1 = 10.0**(4.5535 -0.013486*T1K -799.31/T1K)*H3PO4T/H2PO4T ! Bates (1951) ! SR 01/21/12
     KP2 = 10.0**(5.3541 -0.019840*T1K -1979.5/T1K)*H2PO4T/HPO4T ! Bates and Acree (1943) ! SR 01/21/12
     KP3 = 10.0**(-12.38) *HPO4T/PO4T ! Dean (1985) ! SR 01/01/12
!**** pH evaluation
     PHT = -PH(K,I)-2.1
     IF (PH(K,I) <= 0.0) PHT = -14.0
     INCR = 10.0
     DO N=1,3
        F = 1.0
        INCR = INCR/10.0
        ITER = 0
        DO WHILE (F > 0.0 .AND. ITER < 12)
          PHT = PHT+INCR
          HION = 10.0**PHT
          F = CART*K1*(HION+2.0*K2)/(HION*HION+K1*HION+K1*K2)+KW/HION-ALKT-HION/HT ! SR 01/01/12
          IF (AMMONIA_BUFFERING) THEN ! SR 01/01/12
            F = F + AMMT*KAMM/(HION+KAMM) ! SR 01/01/12
          END IF ! SR 01/01/12
          IF (PHOSPHATE_BUFFERING) THEN ! SR 01/01/12
            F = F + PHOST*( KP1*KP2*HION + 2*KP1*KP2*KP3 - HION*HION*HION ) &
                /( HION*HION*HION + KP1*HION*HION + KP1*KP2*HION + KP1*KP2*KP3) ! SR 01/01/12
          END IF ! SR 01/01/12
          IF (OM_BUFFERING) THEN ! SR 01/01/12
            DO JA=1,NAG ! SR 01/01/12
              F = F + OMCT*SDEN(JA)*( 1.0/(1.0+HION*(10.0**PK(JA))) - 1.0/(1.0+(10.0**(PK(JA)-4.5))) ) ! SR 01/01/12
            END DO ! SR 01/01/12
          END IF ! SR 01/01/12
          ITER = ITER+1
        END DO
        PHT = PHT-INCR
     END DO
!**** pH, carbon dioxide, bicarbonate, and carbonate concentrations
     HION = 10.0**PHT
     PH(K,I) = -PHT
     CO2(K,I) = TIC(K,I)/(1.0+K1/HION+K1*K2/(HION*HION))
     IF(CDWBC(HCO3_DER,JW)=='      ON')HCO3(K,I) = TIC(K,I)/(1.0+HION/K1+K2/HION)
     IF(CDWBC(CO3_DER,JW)=='      ON') CO3(K,I) = TIC(K,I)/((HION*HION)/(K1*K2)+HION/K2+1.0)
   END DO
 END DO
RETURN


!**********************************************************
!**           SUBROUTINE ZOOPLANKTON                     **
!**********************************************************

ENTRY ZOOPLANKTON
  DO I=IU,ID
    DO K=KT,KB(I)
	  DO JZ = 1, NZP
            ZGZTOT=0.0                                                                                                   ! KV 5/9/2007
	        DO JJZ = 1,NZP
!             ZGZTOT=ZGZTOT+ZGZ(K,I,JZ,JJZ)*ZOO(K,I,JZ)                                                                   ! KV 5/9/2007
            ZGZTOT=ZGZTOT+ZGZ(K,I,JZ,JJZ)                                                                             ! CB 5/26/07
            END DO
        ZOOSS(K,I,JZ)= ZSR(K,I,JZ)+ (ZMU(K,I,JZ)*ZEFF(JZ)-ZRT(K,I,JZ)-ZMT(K,I,JZ))*ZOO(K,I,JZ) - ZGZTOT   ! OMNIVOROUS ZOOPLANKTON    ! KV 5/9/2007  ! SW 1/28/2019 SETTLING/RISING
	  END DO
    END DO
  END DO
RETURN

!***********************************************************************************************************************************
!**                                              D E R I V E D   C O N S T I T U E N T S                                          **
!***********************************************************************************************************************************

ENTRY DERIVED_CONSTITUENTS
  APR = 0.0; ATOT = 0.0; TOTSS = 0.0; CHLA = 0.0; CBODU=0.0
  DO JW=1,NWB
    DO JB=BS(JW),BE(JW)
      DO I=CUS(JB),DS(JB)
        DO K=KTWB(JW),KB(I)
          DO JA=1,NAL
            IF(ALG_CALC(JA))APR(K,I) = APR(K,I)+(AGR(K,I,JA)-ARR(K,I,JA))*ALG(K,I,JA)*H2(K,I)*DAY
          END DO
        END DO
        DO K=KTWB(JW),KB(I)
          CBODCT = 0.0; CBODNT = 0.0; CBODPT = 0.0; BODTOT = 0.0; ALGP = 0.0; ALGN = 0.0; ZOOP=0.0; ZOON=0.0; POC(K,I)=0.0
          DO JA=1,NAL
            IF(ALG_CALC(JA))ATOT(K,I) = ATOT(K,I)+ALG(K,I,JA)
          END DO
          DO IBOD=1,NBOD
          IF(BOD_CALC(IBOD))THEN      
            !CBODCt  = CBODCt+CBOD(K,I,IBOD)*BODC(IBOD)    ! cb 6/6/10
            !CBODNt  = CBODNt+CBODn(K,I,IBOD)              ! cb 6/6/10
            !CBODPt  = CBODPt+CBODp(K,I,IBOD)              ! cb 6/6/10
            !BODTOT = BODTOT+CBOD(K,I,IBOD)
            IF(BOD_CALC(IBOD))THEN      
              CBODCT  = CBODCT+CBOD(K,I,IBOD)*RBOD(IBOD)*BODC(IBOD)                          
              BODTOT  = BODTOT+CBOD(K,I,IBOD)*RBOD(IBOD)                                             
            ENDIF
            IF(BOD_CALCP(IBOD)) THEN                                                
              CBODPT  = CBODPT+CBODP(K,I,IBOD)*RBOD(IBOD)                                    
            ELSE 
              CBODPT  = CBODPT+CBOD(K,I,IBOD)*RBOD(IBOD)*BODP(IBOD)                        
            END IF
            IF(BOD_CALCN(IBOD)) THEN                                                
              CBODNT  = CBODNT+CBODN(K,I,IBOD)*RBOD(IBOD)                                   
            ELSE 
              CBODNT  = CBODNT+CBOD(K,I,IBOD)*RBOD(IBOD)*BODN(IBOD)                          
            END IF
            IF(CBODS(IBOD)>0.0)TOTSS(K,I) = TOTSS(K,I)+CBOD(K,I,IBOD)/O2OM(JW)               ! SW 9/5/13  Added particulate CBOD to TSS computation
          ENDIF
          END DO
          IF(ORGC_CALC)THEN
            DOM(K,I) = (LDOMC(K,I)+RDOMC(K,I))/ORGC(JW)
            POM(K,I) = (LPOMC(K,I)+RPOMC(K,I))/ORGC(JW)
          ELSE
          DOM(K,I) = LDOM(K,I)+RDOM(K,I)
          POM(K,I) = LPOM(K,I)+RPOM(K,I)
          END IF
          !
          !DOC(K,I) = DOM(K,I)*ORGC(JW)+CBODCt             ! cb 6/6/10
          !POC(K,I) = POM(K,I)*ORGC(JW)
          DOC(K,I) = LDOC(K,I) + RDOC(K,I) + CBODCT          
          POC(K,I) = LPOC(K,I) + RPOC(K,I)                 
          DO JA=1,NAL
          IF(ALG_CALC(JA))THEN
            POC(K,I) = POC(K,I)+ALG(K,I,JA)*AC(JA)
            ALGP     = ALGP+ALG(K,I,JA)*AP(JA)
            ALGN     = ALGN+ALG(K,I,JA)*AN(JA)
          ENDIF
          END DO
          IF(ZOOPLANKTON_CALC)THEN
            DO JZ=1,NZP
                POC(K,I)=POC(K,I)+ZC(JZ)*ZOO(K,I,JZ) !MLM BAULK
                ZOOP=ZOOP+ZOO(K,I,JZ)*ZP(JZ) !SW 3/25/2021
                ZOON=ZOON+ZOO(K,I,JZ)*ZN(JZ) !SW 3/25/2021
                IF(CDWBC(CBODU_DER,JW)=='      ON')CBODU(K,I) = CBODU(K,I) + O2OM(JW)*ZOO(K,I,JZ)
                IF(CDWBC(TOTSS_DER,JW)=='      ON')TOTSS(K,I) = TOTSS(K,I)+ZOO(K,I,JZ)               ! SW 9/5/13  Added zooplankton to TSS computation
	        END DO
	      ENDIF
          TOC(K,I)   = DOC(K,I)+POC(K,I)
          DOP(K,I)   = LDOP(K,I)+RDOP(K,I)+CBODPT       
          DON(K,I)   = LDON(K,I)+RDON(K,I)+CBODNT
          POP(K,I)   = LPOP(K,I)+RPOP(K,I)+ALGP+ZOOP
          PON(K,I)   = LPON(K,I)+RPON(K,I)+ALGN+ZOON
          !DOP(K,I)   = LDOM(K,I)*ORGPLD(K,I)+RDOM(K,I)*ORGPRD(K,I)+CBODPT    ! CB 6/6/10
          !DON(K,I)   = LDOM(K,I)*ORGNLD(K,I)+RDOM(K,I)*ORGNRD(K,I)+CBODNT    ! CB 6/6/10
          !POP(K,I)   = LPOM(K,I)*ORGPLP(K,I)+RPOM(K,I)*ORGPRP(K,I)+ALGP+ZOOP
          !PON(K,I)   = LPOM(K,I)*ORGNLP(K,I)+RPOM(K,I)*ORGNRP(K,I)+ALGN+ZOON   !SW 1/29/2019 ZOOP
          TOP(K,I)   = DOP(K,I)+POP(K,I)
          TON(K,I)   = DON(K,I)+PON(K,I)
          TKN(K,I)   = TON(K,I)+NH4(K,I)
          IF(CDWBC(NH3_DER,JW)=='      ON')NH3(K,I)   = F_NH3(K,I)*NH4(K,I)
          IF(CDWBC(CBODU_DER,JW)=='      ON')CBODU(K,I) = CBODU(K,I)+O2OM(JW)*(DOM(K,I)+POM(K,I)+ATOT(K,I))+BODTOT
          !TPSS       = 0.0    PO4 ALREADY INCLUDES PARTP  SR 3/17/2019
          !DO JS=1,NSS
          !  TPSS = TPSS+SS(K,I,JS)*PARTP(JW)
          !END DO
          TP(K,I)   =  TOP(K,I)+PO4(K,I)     !+TPSS   SR 3/17/2019
          TN(K,I)   =  TON(K,I)+NH4(K,I)+NO3(K,I)   ! note nh4 is total ammonia including nh3 if ON
          IF(CDWBC(O2DG_DER,JW)=='      ON')O2DG(K,I) = (O2(K,I)/SATO(T1(K,I),TDS(K,I),PALT(I),SALT_WATER(JW)))*100.0          
          IF(CDWBC(CHLA_DER,JW)=='      ON')THEN
            DO JA=1,NAL
               IF(ALG_CALC(JA))THEN
                 CHLA(K,I)  = CHLA(K,I) +ALG(K,I,JA)/ACHLA(JA)
               ENDIF
            END DO
          ENDIF
          IF(CDWBC(TOTSS_DER,JW)=='      ON')THEN
            DO JA=1,NAL
              IF(ALG_CALC(JA))THEN
              TOTSS(K,I) = TOTSS(K,I)+ALG(K,I,JA)
              ENDIF
            END DO
            TOTSS(K,I) = TOTSS(K,I)+TISS(K,I)+POM(K,I)
          ENDIF
          IF(CDWBC(TURB_DER,JW)=='      ON')TURB(K,I)    = EXP(CoeffA_Turb(JW)*LOG(TOTSS(K,I)) + CoeffB_Turb(JW))
          IF(CDWBC(SECCHI_DER,JW)=='      ON')SECCHID(K,I) = SECC_PAR(JW)/GAMMA(K,I)        ! Secchi Disk
          FE(K,I)      = FEII(K,I) + FEOOH(K,I)         ! Total Fe    
          !IF(CAC(NDGP)== '      ON') TDG(K,I) = 100.*DGP(K,I)/PALT(I)       
        END DO
      END DO
    END DO
  END DO
RETURN

!********************************************************************************************************************
!**                                             A L K A L I N I T Y                                                **
!********************************************************************************************************************
ENTRY ALKALINITY ! entire subroutine added ! SR 01/01/12
! According to Stumm and Morgan (1996), table 4.5 on page 173:
! Utilization of ammonium during photosynthesis results in an alkalinity decrease: 14 eq. alk per 16 moles ammonium
! Utilization of nitrate during photosynthesis results in an alkalinity increase: 18 eq. alk per 16 moles nitrate
! Production of ammonium during respiration results in an alkalinity increase: 14 eq. alk per 16 moles ammonium
! Nitrification of ammonium results in an alkalinity decrease: 2 eq. alk per 1 mole ammonium
! Denitrification of nitrate (to nitrogen gas) results in an alkalinity increase: 1 eq. alk per 1 mole nitrate
! Alkalinity is represented as mg/L CaCO3 (MW=100.088). CaCO3 has 2 equivalents of alk per mole.
! Nitrogen has an atomic mass of 14.00674. These numbers account for the factor of 50.044/14.00674 used below.

 
 DO I=IU,ID
   DO K=KT,KB(I)
       if(noncon_alkalinity)then
     ALKSS(K,I) = (50.044/14.00674) * ( 14./16.*(NH4AP(K,I)+NH4EP(K,I)+NH4ZR(K,I)+NH4MR(K,I)-NH4MG(K,I)) &
                  + 18./16.*(NO3AG(K,I)+NO3EG(K,I)) &
                  - 2.*NH4D(K,I) + NO3D(K,I) + NO3SED(K,I)*(1-FNO3SED(JW)) )
       else
         alkss(k,i)=0.0   ! NW 2/11/16
       end if
   END DO
 END DO
 
RETURN

ENTRY DEALLOCATE_KINETICS
  DEALLOCATE (OMTRM,  SODTRM, NH4TRM, NO3TRM, DOM, POM, PO4BOD, NH4BOD, TICBOD, ATRM,   ATRMR,  ATRMF, ETRM,   ETRMR,  ETRMF, BIBH2)
  DEALLOCATE (LAM2M, ALGAE_SETTLING,FE)       
  IF(MIGRATION == 'ON')DEALLOCATE (ASETTLE, DEN_AVG, DENP, DEN1, DEN2, ALLIM_OLD, DEN, AMP, PHASE, C_COEFF_EXT, RAD, MIND, MAXD, DENSI, DENBI, T_DEC, C_DENINC, C_DENDEC, DEPTH_LIM, LOSS_FRAC, TWQ, &
      I_C, C_DENINC_1, C_DENINC_2, C_DENDEC_1, C_DENDEC_2, DENP_MINS, DENP_MINB, DENP_MIN, DEN_COR, MIGRATE_GROUP, MIGRATE_MODEL, TS_DEC, DEPTH_LIM_ONOFF,&
      NMINT, MIGON, MIGOFF, LOLD, DEPTH_CALC_ONOFF, EXP_DEPTH)    ! CO 6/12/2019
  RETURN
END SUBROUTINE KINETICS

