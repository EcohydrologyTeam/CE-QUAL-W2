  !===========================================================================================================================
  ! Read csv input control file and additional parameter inputs are required only if the control switch is turned on.
  ! Updated 9/2020
  !===========================================================================================================================
  Subroutine CEMA_W2_Input
    Use MAIN 
    Use GLOBAL
    Use KINETIC
    Use GEOMC
    Use CEMAVars
    Use SCREENC, ONLY: JDAY
    
    ! Type declarations
    IMPLICIT NONE
    
    Logical SkipLoop         !, file_exists
    Character(256) MessageTemp
    CHARACTER(20) :: ADUMMY   ! SW 2/2019
    
    integer monzz,dayzz,yearzz,ninp,JSKIP
            
    SD_global           = .FALSE.
    IncludeIron         = .FALSE.
    IncludeManganese    = .FALSE.
    IncludeDynamicpH    = .FALSE.
    IncludeAlkalinity   = .FALSE.
    Bubbles_Calculation = .FALSE.    
    
    !INQUIRE(FILE="W2_diagenesis.npt", EXIST=file_exists)   ! file_exists will be TRUE if the file
    IF(SED_DIAG /='      ON')THEN
	      CEMARelatedCode = .FALSE.
	      IncludeBedConsolidation = .FALSE.
	      Return
    ENDIF
    CEMARelatedCode = .TRUE. 
	
    CEMAFilN=NUNIT; NUNIT=NUNIT+1   ! SW 7/8/2019
	  Open(CEMAFilN, File = "W2_diagenesis.npt", STATUS='OLD')
	  Open(CEMALogFilN, File = "DiagenesisLogFile.opt",STATUS='UNKNOWN')
	
	  !Read Header
	  SkipLoop = .FALSE.
	  Do While(.NOT. SkipLoop)
		  Read(CEMAFilN,'(a)')MessageTemp
		  If(index(MessageTemp, "$") == 0)SkipLoop = .TRUE.
	  End Do
	  BackSpace(CEMAFilN)
    !
    ! GROUP 1: Global Control
    Read(CEMAFilN,*)MessageTemp, SD_global
    If(.NOT. SD_global) Then         
	    CEMARelatedCode = .FALSE.
      IncludeFFTLayer = .FALSE.
	  IncludeBedConsolidation  = .FALSE.
      IncludeCEMASedDiagenesis = .FALSE.
	    Return
    End If
    
    ! GROUP 2: FFT Layer
    Read(CEMAFilN,*)MessageTemp, IncludeFFTLayer
    IF(IncludeFFTLayer) THEN
        FirstTimeInFFTCode = .TRUE.
        Read(CEMAFilN,*)MessageTemp, NumFFTActivePrds
        Allocate(FFTActPrdSt(NumFFTActivePrds), FFTActPrdEn(NumFFTActivePrds))
        Allocate(FFTLayConc(IMX))
        Read(CEMAFilN,*)MessageTemp, (FFTActPrdSt(i), i = 1, NumFFTActivePrds) 
        Read(CEMAFilN,*)MessageTemp, (FFTActPrdEn(i), i = 1, NumFFTActivePrds)
        Read(CEMAFilN,*)MessageTemp, InitFFTLayerConc
        Read(CEMAFilN,*)MessageTemp, FFTLayerSettVel
        FFTLayerSettVel = FFTLayerSettVel/DAY !m/d --> m/s
        FFTLayConc = 0.d00
        FFTActPrd = 1
        MoveFFTLayerDown = .FALSE.
        Read(CEMAFilN,*)MessageTemp, MoveFFTLayerDown
    ELSE
        FirstTimeInFFTCode = .FALSE.
        MoveFFTLayerDown   = .FALSE.
        DO JSKIP=1,6
            READ(CEMAFilN,*)
        ENDDO
    END IF
    !
    ! GROUP 3: Bed Consolidation
    Read(CEMAFilN,*)MessageTemp, IncludeBedConsolidation
    IF(IncludeBedConsolidation) THEN 
	    Read(CEMAFilN,*)MessageTemp, LayerAddThkFrac
	    Read(CEMAFilN,*)MessageTemp, NumConsolidRegns
	    Allocate(ConsolidationType(NumConsolidRegns),ConstConsolidRate(NumConsolidRegns))
	    Allocate(ConstPoreWtrRate(NumConsolidRegns),ConsolidRateTemp(NumConsolidRegns))
	    Allocate(ConsRegSegSt(NumConsolidRegns), ConsRegSegEn(NumConsolidRegns))
	    Read(CEMAFilN,*)MessageTemp, (ConsRegSegSt(i), i = 1, NumConsolidRegns)
	    Read(CEMAFilN,*)MessageTemp, (ConsRegSegEn(i), i = 1, NumConsolidRegns)
	    Read(CEMAFilN,*)MessageTemp, (ConsolidationType(i), i = 1, NumConsolidRegns)
	    Read(CEMAFilN,*)MessageTemp, (ConstConsolidRate(i), i = 1, NumConsolidRegns)
	    Read(CEMAFilN,*)MessageTemp, ConsolidRateRegnFil
        Read(CEMAFilN,*)MessageTemp, WriteBESnp
        Read(CEMAFilN,*)MessageTemp, WritePWSnp
    ELSE
        CEMASedimentProcessesInc = .FALSE.
        DO JSKIP=1,9
            READ(CEMAFilN,*)
        ENDDO
    ENDIF
    ! GROUP 4: Sediment Diagenesis
    Read(CEMAFilN,*)MessageTemp, IncludeCEMASedDiagenesis
    IF(IncludeCEMASedDiagenesis) THEN
        SEDIMENT_CALC = .FALSE.
        SOD = 0.0
    END IF
    Read(CEMAFilN,*)MessageTemp, BedElevationInit
	Read(CEMAFilN,*)MessageTemp, BedPorosityInit
    Read(CEMAFilN,*)MessageTemp, CEMAParticleSize
	    CEMAParticleSize = 1.d-6*CEMAParticleSize   !Microns to m
    Read(CEMAFilN,*)MessageTemp, CEMASedimentType
	Read(CEMAFilN,*)MessageTemp, CEMASedimentDensity
	Read(CEMAFilN,*)MessageTemp, CEMASedimentSVelocity
	    CEMASedimentSVelocity = CEMASedimentSVelocity/DAY  !m/d to m/s
    Read(CEMAFilN,*)MessageTemp, CEMASedimentProcessesInc
      !
      Allocate(BedElevation(IMX), BedElevationLayer(IMX), BedPorosity(IMX))
      Allocate(ConsolidRegnNum(IMX), BedConsolidRate(IMX), PorewaterRelRate(IMX))
      Allocate(CEMASedConc(IMX,KMX))
      Allocate(CEMACumPWRelease(IMX), CEMALayerAdded(IMX), CEMASSApplied(IMX))
      Allocate(CEMACumPWToRelease(IMX),CEMACumPWReleased(IMX))
      Allocate(NumCEMAPWInst(IMX))
      Allocate(ApplyCEMAPWRelease(IMX))
      Allocate(CEMACumPWReleaseRate(IMX))
      Allocate(EndBedConsolidation(IMX),BedConsolidationSeg(IMX))  
      Allocate(CEMATSSCopy(KMX,IMX))
      Allocate(VOLCEMA(NBR))
    !
    !
    IF(IncludeCEMASedDiagenesis) THEN
        sediment_diagenesis=.true.
        FirstTimeinCEMAMFTSedDiag = .TRUE.
        Read(CEMAFilN,*)MessageTemp, Bubbles_Calculation
        
            ! GROUP 5: Bubbles
    !IF (.NOT. IncludeCEMASedDiagenesis) Bubbles_Calculation = .FALSE.
    IF(Bubbles_Calculation) THEN    
        Read(CEMAFilN,*)MessageTemp, GasDiff_Sed    ! in m^2/s
        Read(CEMAFilN,*)MessageTemp, CalibParam_R1
        Read(CEMAFilN,*)MessageTemp, YoungModulus
        Read(CEMAFilN,*)MessageTemp, CritStressIF
        Read(CEMAFilN,*)MessageTemp, BubbRelScale
        Read(CEMAFilN,*)MessageTemp, CrackCloseFraction
        Read(CEMAFilN,*)MessageTemp, LimBubbSize
        Read(CEMAFilN,*)MessageTemp, MaxBubbRad
        Read(CEMAFilN,*)MessageTemp, UseReleaseFraction
        Read(CEMAFilN,*)MessageTemp, BubbRelFraction
        Read(CEMAFilN,*)MessageTemp, BubbAccFraction
        Read(CEMAFilN,*)MessageTemp, NumBubRelArr
        Read(CEMAFilN,*)MessageTemp, BubbRelFractionAtm
        Read(CEMAFilN,*)MessageTemp, BubbWatGasExchRate
        Read(CEMAFilN,*)MessageTemp, ApplyBubbTurb
        Read(CEMAFilN,*)MessageTemp, CEMATurbulenceScaling
        Allocate(BubblesCarried(IMX,NumBubRelArr), BubblesRadius(IMX,NumBubRelArr))
	      Allocate(BubblesLNumber(IMX,NumBubRelArr), BubblesStatus(IMX,NumBubRelArr))
	      Allocate(BubblesRiseV(IMX,NumBubRelArr))
	      Allocate(BubblesGasConc(IMX,NumBubRelArr,NumGas))
	      Allocate(BRVoluAGas(IMX,NumBubRelArr,NumGas), BRRateAGas(IMX,NumBubRelArr,NumGas))
	      !Allocate(FirstBubblesRelease(IMX,NumBubRelArr), BubblesReleaseAllValue(IMX,NumBubRelArr))    ! SW 7/1/2017
	      Allocate(BubblesAtSurface(IMX,NumBubRelArr))
    ELSE
        DO JSKIP=1,16
            READ(CEMAFilN,*)
        ENDDO
        LimBubbSize = .FALSE.
        UseReleaseFraction = .FALSE.
        ApplyBubbTurb = .FALSE.
    END IF   
       
        Read(CEMAFilN,*)MessageTemp, CEMA_POM_Resuspension
        
        IF(CEMA_POM_Resuspension) THEN
          Read(CEMAFilN,*)MessageTemp, TAUCRPOM
          Read(CEMAFilN,*)MessageTemp, crshields
          Read(CEMAFilN,*)MessageTemp, cao_method
          Read(CEMAFilN,*)MessageTemp, spgrav_POM
          Read(CEMAFilN,*)MessageTemp, dia_POM
        ELSE
            DO JSKIP=1,5
            READ(CEMAFilN,*)
            ENDDO
        END IF        
        
        Read(CEMAFilN,*)MessageTemp, IncludeAlkalinity
        Read(CEMAFilN,*)MessageTemp, IncludeIron        
        Read(CEMAFilN,*)MessageTemp, IncludeManganese
        !
        IF(IncludeAlkalinity) IncludeDynamicpH = .TRUE.
        !
        !IF(.NOT. IncludeBedConsolidation) THEN
        !    Read(CEMAFilN,*)MessageTemp, BedElevationInit
        !    Read(CEMAFilN,*)MessageTemp, BedPorosityInit
        !    Read(CEMAFilN,*)MessageTemp, CEMASedimentDensity
        !    Allocate(BedElevation(IMX), BedPorosity(IMX),EndBedConsolidation(IMX), PorewaterRelRate(IMX))
        !END IF
        !
        Read(CEMAFilN,*)MessageTemp, NumRegnsSedimentBedComposition
        Allocate(SDRegnPOC_T(NumRegnsSedimentBedComposition), SDRegnPON_T(NumRegnsSedimentBedComposition), SDRegnSul_T(NumRegnsSedimentBedComposition))
        Allocate(SDRegnPOP_T(NumRegnsSedimentBedComposition))
        Allocate(SDRegnH2S_T(NumRegnsSedimentBedComposition), SDRegnNH3_T(NumRegnsSedimentBedComposition), SDRegnCH4_T(NumRegnsSedimentBedComposition))
        Allocate(SDRegnTIC_T(NumRegnsSedimentBedComposition), SDRegnPO4_T(NumRegnsSedimentBedComposition), SDRegnNO3_T(NumRegnsSedimentBedComposition))
        IF(IncludeAlkalinity) Allocate(SDRegnALK_T(NumRegnsSedimentBedComposition))
        IF(IncludeIron)       Allocate(SDRegnFe2_T(NumRegnsSedimentBedComposition),SDRegnFeOOH_T(NumRegnsSedimentBedComposition))
        IF(IncludeManganese)  Allocate(SDRegnMn2_T(NumRegnsSedimentBedComposition),SDRegnMnO2_T(NumRegnsSedimentBedComposition))
        Allocate(SDRegnT_T(NumRegnsSedimentBedComposition))
        IF(.NOT. IncludeDynamicpH) Allocate(SDRegnpH(NumRegnsSedimentBedComposition))
        Allocate(SedBedInitRegSegSt(NumRegnsSedimentBedComposition), SedBedInitRegSegEn(NumRegnsSedimentBedComposition))
        !
        Read(CEMAFilN,*)   ! skip line for header
        Read(CEMAFilN,*)MessageTemp, (SedBedInitRegSegSt(i), i = 1, NumRegnsSedimentBedComposition)
        Read(CEMAFilN,*)MessageTemp, (SedBedInitRegSegEn(i), i = 1, NumRegnsSedimentBedComposition)
        Read(CEMAFilN,*)MessageTemp, (SDRegnT_T(i),   i = 1, NumRegnsSedimentBedComposition)
        IF(.NOT. IncludeDynamicpH) THEN
            Read(CEMAFilN,*)MessageTemp, (SDRegnpH(i), i = 1, NumRegnsSedimentBedComposition) 
        ELSE
            Read(CEMAFilN,*)
        ENDIF
        Read(CEMAFilN,*)MessageTemp, (SDRegnPOC_T(i), i = 1, NumRegnsSedimentBedComposition)
        Read(CEMAFilN,*)MessageTemp, (SDRegnPON_T(i), i = 1, NumRegnsSedimentBedComposition)
        Read(CEMAFilN,*)MessageTemp, (SDRegnPOP_T(i), i = 1, NumRegnsSedimentBedComposition)
        Read(CEMAFilN,*)MessageTemp, (SDRegnSul_T(i), i = 1, NumRegnsSedimentBedComposition)
        Read(CEMAFilN,*)MessageTemp, (SDRegnNH3_T(i), i = 1, NumRegnsSedimentBedComposition)
        Read(CEMAFilN,*)MessageTemp, (SDRegnNO3_T(i), i = 1, NumRegnsSedimentBedComposition)
        Read(CEMAFilN,*)MessageTemp, (SDRegnPO4_T(i), i = 1, NumRegnsSedimentBedComposition)
        Read(CEMAFilN,*)MessageTemp, (SDRegnH2S_T(i), i = 1, NumRegnsSedimentBedComposition)
        Read(CEMAFilN,*)MessageTemp, (SDRegnCH4_T(i), i = 1, NumRegnsSedimentBedComposition)
        Read(CEMAFilN,*)MessageTemp, (SDRegnTIC_T(i), i = 1, NumRegnsSedimentBedComposition)
        IF(IncludeAlkalinity)THEN
            Read(CEMAFilN,*)MessageTemp, (SDRegnALK_T(i), i = 1, NumRegnsSedimentBedComposition)
        ELSE
            READ(CEMAFilN,*)
        ENDIF
        !
        IF(IncludeIron) THEN
          Read(CEMAFilN,*)MessageTemp, (SDRegnFe2_T(i),   i = 1, NumRegnsSedimentBedComposition)
          Read(CEMAFilN,*)MessageTemp, (SDRegnFeOOH_T(i), i = 1, NumRegnsSedimentBedComposition)
        ELSE
            DO JSKIP=1,2
            READ(CEMAFilN,*)
            ENDDO            
        END IF
        IF(IncludeManganese) THEN
          Read(CEMAFilN,*)MessageTemp, (SDRegnMn2_T(i),  i = 1, NumRegnsSedimentBedComposition)
          Read(CEMAFilN,*)MessageTemp, (SDRegnMnO2_T(i), i = 1, NumRegnsSedimentBedComposition)
        ELSE
            DO JSKIP=1,2
            READ(CEMAFilN,*)
            ENDDO            
        END IF
        !
        Read(CEMAFilN,*)MessageTemp, NumRegnsSedimentDiagenesis
        Read(CEMAFilN,*)    ! SKIP LINE FOR HEADER
        Allocate(SDRegnPOC_L_Fr(NumRegnsSedimentDiagenesis),           SDRegnPOC_R_Fr(NumRegnsSedimentDiagenesis),         SDRegnPON_L_Fr(NumRegnsSedimentDiagenesis))
        Allocate(SDRegnPON_R_Fr(NumRegnsSedimentDiagenesis),           SDRegnPW_DiffCoeff(NumRegnsSedimentDiagenesis),     SDRegnOx_Threshold(NumRegnsSedimentDiagenesis))
        Allocate(SDRegnPOP_L_Fr(NumRegnsSedimentDiagenesis),           SDRegnPOP_R_Fr(NumRegnsSedimentDiagenesis))
        Allocate(SDRegnAe_NH3_NO3_L(NumRegnsSedimentDiagenesis),       SDRegnAe_NH3_NO3_H(NumRegnsSedimentDiagenesis),     SDRegnAe_NO3_N2_L(NumRegnsSedimentDiagenesis))
        Allocate(SDRegnAe_NO3_N2_H(NumRegnsSedimentDiagenesis),        SDRegnAn_NO3_N2(NumRegnsSedimentDiagenesis),        SDRegnAe_CH4_CO2(NumRegnsSedimentDiagenesis))
        Allocate(SDRegnAe_HS_NH4_Nit(NumRegnsSedimentDiagenesis),      SDRegnAe_HS_O2_Nit(NumRegnsSedimentDiagenesis),     SDRegn_Theta_PW(NumRegnsSedimentDiagenesis),SDRegn_Theta_PM(NumRegnsSedimentDiagenesis))
        Allocate(SDRegn_Theta_NH3_NO3(NumRegnsSedimentDiagenesis),     SDRegn_Theta_NO3_N2(NumRegnsSedimentDiagenesis),    SDRegn_Theta_CH4_CO2(NumRegnsSedimentDiagenesis))
        Allocate(SDRegn_Sulfate_CH4_H2S(NumRegnsSedimentDiagenesis),   SDRegnAe_H2S_SO4(NumRegnsSedimentDiagenesis),       SDRegn_Theta_H2S_SO4(NumRegnsSedimentDiagenesis))
        Allocate(SDRegn_NormConst_H2S_SO4(NumRegnsSedimentDiagenesis), SDRegn_MinRate_PON_Lab(NumRegnsSedimentDiagenesis), SDRegn_MinRate_PON_Ref(NumRegnsSedimentDiagenesis))
        Allocate(SDRegn_MinRate_PON_Ine(NumRegnsSedimentDiagenesis),   SDRegn_MinRate_POC_Lab(NumRegnsSedimentDiagenesis), SDRegn_MinRate_POC_Ref(NumRegnsSedimentDiagenesis))
        Allocate(SDRegn_MinRate_POC_Ine(NumRegnsSedimentDiagenesis),   SDRegn_Theta_PON_Lab(NumRegnsSedimentDiagenesis),   SDRegn_Theta_PON_Ref(NumRegnsSedimentDiagenesis))
        Allocate(SDRegn_Theta_PON_Ine(NumRegnsSedimentDiagenesis),     SDRegn_Theta_POC_Lab(NumRegnsSedimentDiagenesis),   SDRegn_Theta_POC_Ref(NumRegnsSedimentDiagenesis))
        Allocate(SDRegn_Theta_POC_Ine(NumRegnsSedimentDiagenesis),     SDRegn_CH4CompMethod(NumRegnsSedimentDiagenesis),   SDRegn_POMResuspMethod(NumRegnsSedimentDiagenesis))
        Allocate(SDRegn_Theta_POP_Lab(NumRegnsSedimentDiagenesis),     SDRegn_Theta_POP_Ref(NumRegnsSedimentDiagenesis),   SDRegn_Theta_POP_Ine(NumRegnsSedimentDiagenesis))
        Allocate(SDRegn_MinRate_POP_Lab(NumRegnsSedimentDiagenesis),   SDRegn_MinRate_POP_Ref(NumRegnsSedimentDiagenesis), SDRegn_MinRate_POP_Ine(NumRegnsSedimentDiagenesis))
        Allocate(SedBedDiaRCRegSegSt(NumRegnsSedimentDiagenesis),      SedBedDiaRCRegSegEn(NumRegnsSedimentDiagenesis))
        Allocate(Kdp2(NumRegnsSedimentDiagenesis),KdNH31(NumRegnsSedimentDiagenesis), KdNH32(NumRegnsSedimentDiagenesis)) 
        Allocate(delta_kpo41(NumRegnsSedimentDiagenesis),DOcr(NumRegnsSedimentDiagenesis))
        Allocate(KsOxch(NumRegnsSedimentDiagenesis))
        Allocate(KdH2S1(NumRegnsSedimentDiagenesis),KdH2S2(NumRegnsSedimentDiagenesis))
        Allocate(KdFe1(NumRegnsSedimentDiagenesis), KdFe2(NumRegnsSedimentDiagenesis), KdMn1(NumRegnsSedimentDiagenesis), KdMn2(NumRegnsSedimentDiagenesis))
	      Allocate(PartMixVel(NumRegnsSedimentDiagenesis),BurialVel(NumRegnsSedimentDiagenesis),POCr(NumRegnsSedimentDiagenesis))
        !
        DYNAMIC_SD=.FALSE.
        Read(CEMAFilN,*)MessageTemp, (SedBedDiaRCRegSegSt(i),     i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SedBedDiaRCRegSegEn(i),     i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegnPOC_L_Fr(i),          i = 1, NumRegnsSedimentDiagenesis)
        IF(SDRegnPOC_L_Fr(1) < 0.0)THEN
            DYNAMIC_SD=.TRUE.
            SDRegnPOC_L_Fr(1)=ABS(SDRegnPOC_L_Fr(1))
        ENDIF
        Read(CEMAFilN,*)MessageTemp, (SDRegnPOC_R_Fr(i),          i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegnPON_L_Fr(i),          i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegnPON_R_Fr(i),          i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegnPOP_L_Fr(i),          i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegnPOP_R_Fr(i),          i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegnPW_DiffCoeff(i),      i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (PartMixVel(i),              i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (BurialVel(i),               i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (POCr(i),                    i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegn_CH4CompMethod(i),    i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegnOx_Threshold(i),      i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegnAe_NH3_NO3_L(i),      i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegnAe_NH3_NO3_H(i),      i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegnAe_NO3_N2_L(i),       i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegnAe_NO3_N2_H(i),       i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegnAn_NO3_N2(i),         i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegnAe_CH4_CO2(i),        i = 1, NumRegnsSedimentDiagenesis)   !Eq. 10.35
        Read(CEMAFilN,*)MessageTemp, (KsOxch(i),                  i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegnAe_HS_NH4_Nit(i),     i = 1, NumRegnsSedimentDiagenesis)   !Eq. 3.3
        Read(CEMAFilN,*)MessageTemp, (SDRegnAe_HS_O2_Nit(i),      i = 1, NumRegnsSedimentDiagenesis)   !Eq. 3.3
        Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_PW(i),         i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_PM(i),         i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_NH3_NO3(i),    i = 1, NumRegnsSedimentDiagenesis)   
        Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_NO3_N2(i),     i = 1, NumRegnsSedimentDiagenesis)  
        Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_CH4_CO2(i),    i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegn_Sulfate_CH4_H2S(i),  i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegnAe_H2S_SO4(i),        i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_H2S_SO4(i),    i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegn_NormConst_H2S_SO4(i),i = 1, NumRegnsSedimentDiagenesis)   !Eq. 9.6
        Read(CEMAFilN,*)MessageTemp, (SDRegn_MinRate_POC_Lab(i),  i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegn_MinRate_POC_Ref(i),  i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegn_MinRate_POC_Ine(i),  i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegn_MinRate_PON_Lab(i),  i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegn_MinRate_PON_Ref(i),  i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegn_MinRate_PON_Ine(i),  i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegn_MinRate_POP_Lab(i),  i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegn_MinRate_POP_Ref(i),  i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegn_MinRate_POP_Ine(i),  i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_POC_Lab(i),    i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_POC_Ref(i),    i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_POC_Ine(i),    i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_PON_Lab(i),    i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_PON_Ref(i),    i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_PON_Ine(i),    i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_POP_Lab(i),    i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_POP_Ref(i),    i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (SDRegn_Theta_POP_Ine(i),    i = 1, NumRegnsSedimentDiagenesis)    
        Read(CEMAFilN,*)MessageTemp, (Kdp2(i),                    i = 1, NumRegnsSedimentDiagenesis) 
        Read(CEMAFilN,*)MessageTemp, (delta_kpo41(i),             i = 1, NumRegnsSedimentDiagenesis) 
        Read(CEMAFilN,*)MessageTemp, (DOcr(i),                    i = 1, NumRegnsSedimentDiagenesis) 
        Read(CEMAFilN,*)MessageTemp, (KdNH31(i),                  i = 1, NumRegnsSedimentDiagenesis)    
        Read(CEMAFilN,*)MessageTemp, (KdNH32(i),                  i = 1, NumRegnsSedimentDiagenesis) 
        Read(CEMAFilN,*)MessageTemp, (KdH2S1(i),                  i = 1, NumRegnsSedimentDiagenesis) 
        Read(CEMAFilN,*)MessageTemp, (KdH2S2(i),                  i = 1, NumRegnsSedimentDiagenesis) 
        Read(CEMAFilN,*)MessageTemp, (SDRegn_POMResuspMethod(i),  i = 1, NumRegnsSedimentDiagenesis)

        Read(CEMAFilN,*)MessageTemp, (KdFe1(i),               i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (KdFe2(i),               i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (KdMn1(i),               i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*)MessageTemp, (KdMn2(i),               i = 1, NumRegnsSedimentDiagenesis)
    ! GROUP 6: Output
        Read(CEMAFilN,*)MessageTemp, WriteCEMAMFTSedFlx
        READ(CEMAFilN,*)MessageTemp, SEDIAGFREQ   ! FREQUENCY OF OUTPUT SW 5/25/2017
    ELSE
        IncludeDynamicpH      = .FALSE.
        IncludeAlkalinity     = .FALSE.
        CEMA_POM_Resuspension = .FALSE.
        IncludeIron           = .FALSE.       
        IncludeManganese      = .FALSE.
        cao_method            = .FALSE. 
    END IF

    close(CEMAFilN)
    !
    IF(IncludeCEMASedDiagenesis .and. WriteCEMAMFTSedFlx) THEN
                    IF(RESTART_IN)THEN
                    Open(CEMASedFlxFilN4, File = 'DiagenesisSOD.csv', POSITION='APPEND')
                    JDAY1=0.0
                    REWIND (CEMASedFlxFilN4)  
                    READ   (CEMASedFlxFilN4,'(/)',END=101)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN4,'(A,F12.0)',END=101)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN4)  
101                JDAY1 = 0.0  
                    Open(CEMASedFlxFilN5, File = 'Diagenesis_POCG1.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN5)  
                    READ   (CEMASedFlxFilN5,'(/)',END=102)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN5,'(A,F12.0)',END=102)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN5)  
102                JDAY1 = 0.0  
                    Open(CEMASedFlxFilN6, File = 'Diagenesis_POCG2.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN6)  
                    READ   (CEMASedFlxFilN6,'(/)',END=103)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN6,'(A,F12.0)',END=103)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN6)  
103                JDAY1 = 0.0              
                    Open(CEMASedFlxFilN7, File = 'Diagenesis_JC.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN7)  
                    READ   (CEMASedFlxFilN7,'(/)',END=104)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN7,'(A,F12.0)',END=104)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN7)  
104                JDAY1 = 0.0             
                   Open(CEMASedFlxFilN8, File = 'Diagenesis_JN.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN8)  
                    READ   (CEMASedFlxFilN8,'(/)',END=105)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN8,'(A,F12.0)',END=105)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN8)  
105                JDAY1 = 0.0             
                    Open(CEMASedFlxFilN9, File = 'Diagenesis_PONG1.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN9)  
                    READ   (CEMASedFlxFilN9,'(/)',END=106)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN9,'(A,F12.0)',END=106)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN9)  
106                JDAY1 = 0.0             
                    Open(CEMASedFlxFilN10, File = 'Diagenesis_PONG2.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN10)  
                    READ   (CEMASedFlxFilN10,'(/)',END=107)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN10,'(A,F12.0)',END=107)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN10)  
107                JDAY1 = 0.0             
                    Open(CEMASedFlxFilN11, File = 'Diagenesis_SD_JCH4.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN11)  
                    READ   (CEMASedFlxFilN11,'(/)',END=108)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN11,'(A,F12.0)',END=108)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN11)  
108                JDAY1 = 0.0             
                    Open(CEMASedFlxFilN12, File = 'Diagenesis_SD_JNH4.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN12)  
                    READ   (CEMASedFlxFilN12,'(/)',END=109)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN12,'(A,F12.0)',END=109)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN12)  
109                 JDAY1 = 0.0                 
                    Open(CEMASedFlxFilN13, File = 'Diagenesis_SD_JNO3.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN13)  
                    READ   (CEMASedFlxFilN13,'(/)',END=110)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN13,'(A,F12.0)',END=110)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN13)  
110                JDAY1 = 0.0                 
                    Open(CEMASedFlxFilN14, File = 'Diagenesis_SD_JPO4.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN14)  
                    READ   (CEMASedFlxFilN14,'(/)',END=111)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN14,'(A,F12.0)',END=111)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN14)  
111                JDAY1 = 0.0                 
                    Open(CEMASedFlxFilN15, File = 'Diagenesis_POPG1.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN15)  
                    READ   (CEMASedFlxFilN15,'(/)',END=112)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN15,'(A,F12.0)',END=112)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN15)  
112                JDAY1 = 0.0                 
                    Open(CEMASedFlxFilN16, File = 'Diagenesis_POPG2.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN16)  
                    READ   (CEMASedFlxFilN16,'(/)',END=113)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN16,'(A,F12.0)',END=113)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN16)  
113                JDAY1 = 0.0                 
                    Open(CEMASedFlxFilN17, File = 'DiagenesisCSOD.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN17)  
                    READ   (CEMASedFlxFilN17,'(/)',END=114)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN17,'(A,F12.0)',END=114)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN17)  
114                JDAY1 = 0.0                 
                    Open(CEMASedFlxFilN18, File = 'DiagenesisNSOD.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN18)  
                    READ   (CEMASedFlxFilN18,'(/)',END=115)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN18,'(A,F12.0)',END=115)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN18)  
115                JDAY1 = 0.0                 
                    Open(CEMASedFlxFilN19, File = 'Diagenesis_JP.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN19)  
                    READ   (CEMASedFlxFilN19,'(/)',END=116)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN19,'(A,F12.0)',END=116)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN19)  
116                JDAY1 = 0.0                 
                    Open(CEMASedFlxFilN20, File = 'DiagenesisAerobicLayer.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN20)  
                    READ   (CEMASedFlxFilN20,'(/)',END=117)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN20,'(A,F12.0)',END=117)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN20)  
117                 JDAY1 = 0.0    
                    Open(CEMASedFlxFilN21, File = 'Diagenesis_TemperatureAerobicLayer.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN21)  
                    READ   (CEMASedFlxFilN21,'(/)',END=118)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN21,'(A,F12.0)',END=118)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN21)  
118                 JDAY1 = 0.0   
                    Open(CEMASedFlxFilN22, File = 'Diagenesis_TemperatureAnaerobicLayer.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN22)  
                    READ   (CEMASedFlxFilN22,'(/)',END=119)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN22,'(A,F12.0)',END=119)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN22)  
119                 JDAY1 = 0.0 
                    
                    Open(CEMASedFlxFilN23, File = 'Diagenesis_NO3AerobicLayer.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN23)  
                    READ   (CEMASedFlxFilN23,'(/)',END=120)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN23,'(A,F12.0)',END=120)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN23)  
120                 JDAY1 = 0.0   
                    Open(CEMASedFlxFilN24, File = 'Diagenesis_NO3AnaerobicLayer.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN24)  
                    READ   (CEMASedFlxFilN24,'(/)',END=121)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN24,'(A,F12.0)',END=121)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN24)  
121                 JDAY1 = 0.0 
                    Open(CEMASedFlxFilN25, File = 'Diagenesis_NH3AerobicLayer.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN25)  
                    READ   (CEMASedFlxFilN25,'(/)',END=122)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN25,'(A,F12.0)',END=122)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN25)  
122                 JDAY1 = 0.0   
                    Open(CEMASedFlxFilN26, File = 'Diagenesis_NH3AnaerobicLayer.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN26)  
                    READ   (CEMASedFlxFilN26,'(/)',END=123)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN26,'(A,F12.0)',END=123)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN26)  
123                 JDAY1 = 0.0 
                    Open(CEMASedFlxFilN27, File = 'Diagenesis_PO4AerobicLayer.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN27)  
                    READ   (CEMASedFlxFilN27,'(/)',END=124)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN27,'(A,F12.0)',END=124)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN27)  
124                 JDAY1 = 0.0   
                    Open(CEMASedFlxFilN28, File = 'Diagenesis_PO4AnaerobicLayer.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN28)  
                    READ   (CEMASedFlxFilN28,'(/)',END=125)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN28,'(A,F12.0)',END=125)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN28)  
125                 JDAY1 = 0.0 
                    Open(CEMASedFlxFilN29, File = 'Diagenesis_SO4AerobicLayer.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN29)  
                    READ   (CEMASedFlxFilN29,'(/)',END=126)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN29,'(A,F12.0)',END=126)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN29)  
126                 JDAY1 = 0.0   
                    Open(CEMASedFlxFilN30, File = 'Diagenesis_SO4AnaerobicLayer.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN30)  
                    READ   (CEMASedFlxFilN30,'(/)',END=127)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN30,'(A,F12.0)',END=127)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN30)  
127                 JDAY1 = 0.0           
                    
                    Open(CEMASedFlxFilN31, File = 'Diagenesis_FeIIAerobicLayer.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN31)  
                    READ   (CEMASedFlxFilN31,'(/)',END=128)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN31,'(A,F12.0)',END=128)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN31)  
128                 JDAY1 = 0.0           
                    Open(CEMASedFlxFilN32, File = 'Diagenesis_FeIIAnaerobicLayer.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN32)  
                    READ   (CEMASedFlxFilN32,'(/)',END=129)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN32,'(A,F12.0)',END=129)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN32)  
129                 JDAY1 = 0.0           
                    Open(CEMASedFlxFilN33, File = 'Diagenesis_MnIIAerobicLayer.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN33)  
                    READ   (CEMASedFlxFilN33,'(/)',END=130)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN33,'(A,F12.0)',END=130)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN33)  
130                 JDAY1 = 0.0           
                    Open(CEMASedFlxFilN34, File = 'Diagenesis_MnIIAnaerobicLayer.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN34)  
                    READ   (CEMASedFlxFilN34,'(/)',END=131)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN34,'(A,F12.0)',END=131)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN34)  
131                 JDAY1 = 0.0           
                    Open(CEMASedFlxFilN35, File = 'Diagenesis_CH4AerobicLayer.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN35)  
                    READ   (CEMASedFlxFilN35,'(/)',END=132)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN35,'(A,F12.0)',END=132)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN35)  
132                 JDAY1 = 0.0           
                    Open(CEMASedFlxFilN36, File = 'Diagenesis_CH4AnaerobicLayer.csv', POSITION='APPEND')
                    REWIND (CEMASedFlxFilN36)  
                    READ   (CEMASedFlxFilN36,'(/)',END=133)  
                    DO WHILE (JDAY1 < JDAY)  
                     READ (CEMASedFlxFilN36,'(A,F12.0)',END=133)ADUMMY, JDAY1  
                    END DO  
                    BACKSPACE (CEMASedFlxFilN36)  
133                 JDAY1 = 0.0           
        
            ELSE
        
        Open(CEMASedFlxFilN4, File = "DiagenesisSOD.csv", STATUS='unknown')
        Write(CEMASedFlxFilN4,'("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))')(SegNumI, SegNumI = 1, IMX),(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN5, File = "Diagenesis_POCG1.csv", STATUS='unknown')
        Write(CEMASedFlxFilN5,'("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))')(SegNumI, SegNumI = 1, IMX),(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN6, File = "Diagenesis_POCG2.csv", STATUS='unknown')
        Write(CEMASedFlxFilN6,'("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))')(SegNumI, SegNumI = 1, IMX),(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN7, File = "Diagenesis_JC.csv", STATUS='unknown')
        Write(CEMASedFlxFilN7,'("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))')(SegNumI, SegNumI = 1, IMX),(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN8, File = "Diagenesis_JN.csv", STATUS='unknown')
        Write(CEMASedFlxFilN8,'("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))')(SegNumI, SegNumI = 1, IMX),(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN9, File = "Diagenesis_PONG1.csv", STATUS='unknown')
        Write(CEMASedFlxFilN9,'("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))')(SegNumI, SegNumI = 1, IMX),(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN10, File = "Diagenesis_PONG2.csv", STATUS='unknown')
        Write(CEMASedFlxFilN10,'("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))')(SegNumI, SegNumI = 1, IMX),(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN11, File = "Diagenesis_SD_JCH4.csv", STATUS='unknown')
        Write(CEMASedFlxFilN11,'("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))')(SegNumI, SegNumI = 1, IMX),(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN12, File = "Diagenesis_SD_JNH4.csv", STATUS='unknown')
        Write(CEMASedFlxFilN12,'("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))')(SegNumI, SegNumI = 1, IMX),(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN13, File = "Diagenesis_SD_JNO3.csv", STATUS='unknown')
        Write(CEMASedFlxFilN13,'("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))')(SegNumI, SegNumI = 1, IMX),(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN14, File = "Diagenesis_SD_JPO4.csv", STATUS='unknown')
        Write(CEMASedFlxFilN14,'("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))')(SegNumI, SegNumI = 1, IMX),(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN15, File = "Diagenesis_POPG1.csv", STATUS='unknown')
        Write(CEMASedFlxFilN15,'("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))')(SegNumI, SegNumI = 1, IMX),(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN16, File = "Diagenesis_POPG2.csv", STATUS='unknown')
        Write(CEMASedFlxFilN16,'("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))')(SegNumI, SegNumI = 1, IMX),(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN17, File = "DiagenesisCSOD.csv", STATUS='unknown')
        Write(CEMASedFlxFilN17,'("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))')(SegNumI, SegNumI = 1, IMX),(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN18, File = "DiagenesisNSOD.csv", STATUS='unknown')
        Write(CEMASedFlxFilN18,'("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))')(SegNumI, SegNumI = 1, IMX),(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN19, File = "Diagenesis_JP.csv", STATUS='unknown')
        Write(CEMASedFlxFilN19,'("Variable,JDAY,",<IMX>(i5,","),<IMX>(i6,","))')(SegNumI, SegNumI = 1, IMX),(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN20, File = 'DiagenesisAerobicLayer.csv', STATUS='unknown')
        WRITE(CEMASedFlxFilN20,'("Variable,JDAY,",<IMX>(i5,","))')(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN21, File = 'Diagenesis_TemperatureAerobicLayer.csv', STATUS='unknown')
        Write(CEMASedFlxFilN21,'("Variable,JDAY,",<IMX>(i5,","))')(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN22, File = 'Diagenesis_TemperatureAnaerobicLayer.csv', STATUS='unknown')  
        Write(CEMASedFlxFilN22,'("Variable,JDAY,",<IMX>(i5,","))')(SegNumI, SegNumI = 1, IMX)
        
        Open(CEMASedFlxFilN23, File = 'Diagenesis_NO3AerobicLayer.csv', STATUS='unknown')
        Write(CEMASedFlxFilN23,'("Variable,JDAY,",<IMX>(i5,","))')(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN24, File = 'Diagenesis_NO3AnaerobicLayer.csv', STATUS='unknown')  
        Write(CEMASedFlxFilN24,'("Variable,JDAY,",<IMX>(i5,","))')(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN25, File = 'Diagenesis_NH3AerobicLayer.csv', STATUS='unknown')
        Write(CEMASedFlxFilN25,'("Variable,JDAY,",<IMX>(i5,","))')(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN26, File = 'Diagenesis_NH3AnaerobicLayer.csv', STATUS='unknown')  
        Write(CEMASedFlxFilN26,'("Variable,JDAY,",<IMX>(i5,","))')(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN27, File = 'Diagenesis_PO4AerobicLayer.csv', STATUS='unknown')
        Write(CEMASedFlxFilN27,'("Variable,JDAY,",<IMX>(i5,","))')(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN28, File = 'Diagenesis_PO4AnaerobicLayer.csv', STATUS='unknown')  
        Write(CEMASedFlxFilN28,'("Variable,JDAY,",<IMX>(i5,","))')(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN29, File = 'Diagenesis_SO4AerobicLayer.csv', STATUS='unknown')
        Write(CEMASedFlxFilN29,'("Variable,JDAY,",<IMX>(i5,","))')(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN30, File = 'Diagenesis_SO4AnaerobicLayer.csv', STATUS='unknown')  
        Write(CEMASedFlxFilN30,'("Variable,JDAY,",<IMX>(i5,","))')(SegNumI, SegNumI = 1, IMX)
        
        Open(CEMASedFlxFilN31, File = 'Diagenesis_FeIIAerobicLayer.csv', STATUS='unknown')
        Write(CEMASedFlxFilN31,'("Variable,JDAY,",<IMX>(i5,","))')(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN32, File = 'Diagenesis_FeIIAnaerobicLayer.csv', STATUS='unknown')  
        Write(CEMASedFlxFilN32,'("Variable,JDAY,",<IMX>(i5,","))')(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN33, File = 'Diagenesis_MnIIAerobicLayer.csv', STATUS='unknown')
        Write(CEMASedFlxFilN33,'("Variable,JDAY,",<IMX>(i5,","))')(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN34, File = 'Diagenesis_MnIIAnaerobicLayer.csv', STATUS='unknown')  
        Write(CEMASedFlxFilN34,'("Variable,JDAY,",<IMX>(i5,","))')(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN35, File = 'Diagenesis_CH4AerobicLayer.csv', STATUS='unknown')
        Write(CEMASedFlxFilN35,'("Variable,JDAY,",<IMX>(i5,","))')(SegNumI, SegNumI = 1, IMX)
        Open(CEMASedFlxFilN36, File = 'Diagenesis_CH4AnaerobicLayer.csv', STATUS='unknown')  
        Write(CEMASedFlxFilN36,'("Variable,JDAY,",<IMX>(i5,","))')(SegNumI, SegNumI = 1, IMX)

            ENDIF
            
    END IF
    !
    ! Fix the values here
    NH4_NH3_Eqb_Const = 9.1
    HS_H2S_Eqb_Const  = 9.0
    HenryConst_NH3 = 0.0179
    HenryConst_CH4 = 469.0 
    HenryConst_H2S = 10.0
    HenryConst_CO2 = 29.0
    !
    !Allocate other variablesnd 
    ALLOCATE(CellArea(KMX,IMX))
    IF(IncludeCEMASedDiagenesis) THEN
      Allocate(CEMAMFT_RandC_RegN(IMX), CEMAMFT_InCond_RegN(IMX), MFTSedFlxVars(KMX,IMX,59), CEMA_SD_Vars(KMX,IMX,22))
      Allocate(SD_NO3p2(2), SD_NH3p2(2), SD_NH3Tp2(2), SD_CH4p2(2), SD_PO4p2(2), SD_PO4Tp2(2),SD_PO4(2))
	    Allocate(SD_HSp2(2),  SD_HSTp2(2))
	    Allocate(SD_poc2(3),  SD_pon2(3), SD_pop2(3), SD_NH3Tp(2), SD_NO3p(2), SD_PO4Tp(2), SD_HSTp(2))
	    Allocate(SD_fpon(3),  SD_fpoc(3), SD_kdiaPON(3), SD_ThtaPON(3), SD_kdiaPOC(3), SD_ThtaPOC(3))
	    Allocate(SD_JPOC(3),  SD_JPON(3), SD_JPOP(3), SD_TDS(2))
      Allocate(SD_EPOC(3),  SD_EPON(3), SD_EPOP(3))
      Allocate(SD_Denit(2), SD_JDenit(2), SD_JO2NO3(2),SD_HS(2))   
      IF(IncludeIron)         Allocate(SD_Fe2(2))
      IF(IncludeManganese)    Allocate(SD_Mn2(2))
	    Allocate(SD_kdiaPOP(3), SD_ThtaPOP(3), SD_NH3T(2), SD_FPOP(3))
      Allocate(SD_pHValue(IMX))   
	    Allocate(SD_AerLayerThick(IMX))
	    IF(Bubbles_Calculation) THEN
        Allocate(H2SDis(IMX), H2SGas(IMX), CH4Dis(IMX), CH4Gas(IMX))
	      Allocate(NH4Dis(IMX), NH4Gas(IMX), CO2Dis(IMX), CO2Gas(IMX))
        Allocate(BubbleRadiusSed(IMX),PresBubbSed(IMX), PresCritSed(IMX))
	      Allocate(CgSed(IMX), C0Sed(IMX), CtSed(IMX))
	      Allocate(TConc(NumGas,KMX,IMX), TConcP(NumGas,KMX,IMX), SConc(NumGas,KMX,IMX))
	      Allocate(DissolvedGasSediments(NumGas,KMX,IMX))
	      Allocate(CrackOpen(IMX), MFTBubbReleased(IMX), LastDiffVolume(IMX))
	      Allocate(BottomTurbulence(IMX))
        Allocate(FirstBubblesRelease(IMX,NumBubRelArr), BubblesReleaseAllValue(IMX,NumBubRelArr),BubbleRelWB(NWB,NumGas))    ! SW 7/1/2017
        Allocate(BRRateAGasNet(IMX, NumGas))
      END IF
	    Allocate(SDPFLUX(NWB),SDNH4FLUX(NWB),SDNO3FLUX(NWB))
    END IF    
    
    Return
  End Subroutine
    
    
  SUBROUTINE INIT_CEMA
    USE CEMAVars; USE MAIN
    IMPLICIT NONE
    
    IF(IncludeCEMASedDiagenesis) THEN
        SD_NO3p2   = 0.d00; SD_NH3p2   = 0.d00; SD_NH3Tp2  = 0.d00;  SD_CH4p2 = 0.d00 
        SD_PO4p2   = 0.d00; SD_PO4Tp2  = 0.d00; SD_HSp2    = 0.d00;  SD_HSTp2 = 0.d00 
        SD_POC2    = 0.d00; SD_PON2    = 0.d00; SD_POP2    = 0.d00;  SD_NH3Tp = 0.d00 
        SD_NO3p    = 0.d00; SD_PO4Tp   = 0.d00; SD_HSTp    = 0.d00 
        SD_FPON    = 0.d00; SD_FPOC    = 0.d00; SD_kdiaPON = 0.d00;  SD_ThtaPON = 0.d00 
        SD_kdiaPOC = 0.d00; SD_ThtaPOC = 0.d00 
        SD_JPOC    = 0.d00; SD_JPON    = 0.d00; SD_JPOP    = 0.d00 
        SD_EPOC    = 0.d00; SD_EPON    = 0.d00; SD_EPOP    = 0.d00
        SD_Denit   = 0.d00; SD_JDenit  = 0.d00; SD_JO2NO3  = 0.d00
        SD_PO4     = 0.d00; SD_FPOP    = 0.d00; SD_HS      = 0.d00 
        !
        IF(IncludeIron) THEN
            SD_Fe2 = 0.d00
        END IF
        IF(IncludeManganese) THEN
            SD_Mn2 = 0.d00
        END IF
        SD_kdiaPOP = 0.d00; SD_ThtaPOP = 0.d00; SD_NH3T = 0.d00
        SD_AerLayerThick = 0.d00
        !
        IF(Bubbles_Calculation) THEN
            H2SDis = 0.d00; H2SGas = 0.d00; CH4Dis = 0.d00; CH4Gas = 0.d00 
            NH4Dis = 0.d00; NH4Gas = 0.d00; CO2Dis = 0.d00; CO2Gas = 0.d00 
            BubbleRadiusSed = 0.d00; PresBubbSed = 0.d00; PresCritSed = 0.d00 
            CgSed = 0.d00; C0Sed = 0.d00; CtSed = 0.d00; TConcP = 0.d00
            LastDiffVolume = 0.d00
            BubblesCarried = 0; BubblesLNumber = 0; BubblesStatus = 0
            BubblesRadius = 0.d00; BubblesRiseV = 0.d00; BubblesGasConc = 0.d00
            BubblesReleaseAllValue = 0.d00
            BRVoluAGas = 0.d00; BRRateAGas = 0.d00; BRRateAGasNet = 0.d00
            BottomTurbulence = 0.d00
            DissolvedGasSediments = 0.d00
            FirstTimeInBubbles  = .TRUE. 
            FirstBubblesRelease = .TRUE. 
            BubblesAtSurface    = .FALSE.
        END IF
        CEMAMFT_RandC_RegN = 0
        CEMA_SD_Vars = 0.d00
        SDPFLUX=0.0
        SDNH4FLUX=0.0
        SDNO3FLUX=0.0
    END IF
    !
    IF(IncludeBedConsolidation) THEN
        BedElevationLayer = 0.d00
        BedConsolidRate = 0.d00
        PorewaterRelRate = 0.d00
        CEMASedConc = 0.d00
        CEMACumPWRelease = 0.d00
        CEMACumPWReleaseRate = 0.d00
        CEMACumPWToRelease = 0.d00
        CEMACumPWReleased = 0.d00
        EndBedConsolidation = .FALSE.
        BedConsolidationSeg = .FALSE.  ! cb 6/28/18
        VOLCEMA = 0.d00
        NumCEMAPWInst = 0
        ApplyCEMAPWRelease = .FALSE.
    END IF
    !
    IF(IncludeCEMASedDiagenesis .and. (.not. IncludeBedConsolidation)) THEN
        EndBedConsolidation = .FALSE.
        PorewaterRelRate = 0.d00
    END IF
    IF(IncludeCEMASedDiagenesis .OR. IncludeBedConsolidation) BedElevation = BedElevationInit
    !
    IF(.NOT.RESTART_IN)THEN
        CellArea=0.0
        IF(IncludeCEMASedDiagenesis) THEN
            MFTSedFlxVars = 0.d00
            BedPorosity = BedPorosityInit
            IF(Bubbles_Calculation) THEN 
                MFTBubbReleased = 0 
                TConc = 0.d00; SConc = 0.d00
                CrackOpen = .FALSE.; BubbleRelWB=0.0
                GasReleaseCH4=0.0
            END IF
        END IF
    ENDIF
  END SUBROUTINE INIT_CEMA
    
  Subroutine Deallocate_CEMA
    Use CEMAVars
    IMPLICIT NONE
    !
    DEALLOCATE (CellArea)
    IF(IncludeBedConsolidation) THEN
        DEALLOCATE(ConsolidationType,ConstConsolidRate)
        DEALLOCATE(ConstPoreWtrRate, ConsolidRateTemp)
	    DEALLOCATE(ConsRegSegSt, ConsRegSegEn)   
    ENDIF        
        DEALLOCATE(ConsolidRegnNum, BedConsolidRate, PorewaterRelRate)   
        DEALLOCATE(CEMASedConc)
        DEALLOCATE(CEMACumPWRelease, CEMALayerAdded, CEMASSApplied)
        DEALLOCATE(CEMACumPWToRelease,CEMACumPWReleased)
        DEALLOCATE(NumCEMAPWInst)
        DEALLOCATE(ApplyCEMAPWRelease)
        DEALLOCATE(CEMACumPWReleaseRate)
        DEALLOCATE(EndBedConsolidation)
        DEALLOCATE(BedConsolidationSeg)  ! cb 6/28/18
        DEALLOCATE(CEMATSSCopy)
        DEALLOCATE(VOLCEMA)
        DEALLOCATE(BedElevationLayer)
    !END IF
    !
    IF(IncludeCEMASedDiagenesis .OR. IncludeBedConsolidation) DEALLOCATE(BedElevation, BedPorosity)
    IF(IncludeFFTLayer) THEN
      DEALLOCATE(FFTActPrdSt, FFTActPrdEn)
      DEALLOCATE(FFTLayConc)
    END IF
    IF(IncludeCEMASedDiagenesis) THEN
        DEALLOCATE(SDRegnPOC_T, SDRegnPON_T, SDRegnSul_T)
        DEALLOCATE(SDRegnPOP_T)
        DEALLOCATE(SDRegnH2S_T, SDRegnNH3_T, SDRegnCH4_T, SDRegnNO3_T)
        IF(IncludeAlkalinity) DEALLOCATE(SDRegnALK_T)
        IF(.NOT. IncludeDynamicpH)  DEALLOCATE(SDRegnpH)
        DEALLOCATE(SDRegnTIC_T, SDRegnPO4_T)
        IF(IncludeIron) DEALLOCATE(SDRegnFe2_T,SDRegnFeOOH_T, SD_Fe2, KdFe1, KdFe2)
        IF(IncludeManganese) DEALLOCATE(SDRegnMn2_T,SDRegnMnO2_T, SD_Mn2, KdMn1, KdMn2)
        DEALLOCATE(SDRegnT_T)
        DEALLOCATE(SedBedInitRegSegSt, SedBedInitRegSegEn)
        DEALLOCATE(SDRegnPOC_L_Fr,          SDRegnPOC_R_Fr,         SDRegnPON_L_Fr)
        DEALLOCATE(SDRegnPON_R_Fr,          SDRegnPW_DiffCoeff,     SDRegnOx_Threshold)
        DEALLOCATE(SDRegnPOP_L_Fr,          SDRegnPOP_R_Fr)
        DEALLOCATE(SDRegnAe_NH3_NO3_L,      SDRegnAe_NH3_NO3_H,     SDRegnAe_NO3_N2_L)
        DEALLOCATE(SDRegnAe_NO3_N2_H,       SDRegnAn_NO3_N2,        SDRegnAe_CH4_CO2)
        DEALLOCATE(SDRegnAe_HS_NH4_Nit,     SDRegnAe_HS_O2_Nit,     SDRegn_Theta_PW,SDRegn_Theta_PM)
        DEALLOCATE(SDRegn_Theta_NH3_NO3,    SDRegn_Theta_NO3_N2,    SDRegn_Theta_CH4_CO2)
        DEALLOCATE(SDRegn_Sulfate_CH4_H2S,  SDRegnAe_H2S_SO4,       SDRegn_Theta_H2S_SO4)
        DEALLOCATE(SDRegn_NormConst_H2S_SO4,SDRegn_MinRate_PON_Lab, SDRegn_MinRate_PON_Ref)
        DEALLOCATE(SDRegn_MinRate_PON_Ine,  SDRegn_MinRate_POC_Lab, SDRegn_MinRate_POC_Ref)
        DEALLOCATE(SDRegn_MinRate_POC_Ine,  SDRegn_Theta_PON_Lab,   SDRegn_Theta_PON_Ref)
        DEALLOCATE(SDRegn_Theta_PON_Ine,    SDRegn_Theta_POC_Lab,   SDRegn_Theta_POC_Ref)
        DEALLOCATE(SDRegn_Theta_POC_Ine,    SDRegn_CH4CompMethod,   SDRegn_POMResuspMethod)
        DEALLOCATE(SDRegn_Theta_POP_Lab,    SDRegn_Theta_POP_Ref,   SDRegn_Theta_POP_Ine)
        DEALLOCATE(Kdp2, KdNH31, KdNH32, KdH2S1, KdH2S2, delta_kpo41, DOcr)
        DEALLOCATE(PartMixVel,BurialVel,POCr,KsOxch)
        DEALLOCATE(SDRegn_MinRate_POP_Lab, SDRegn_MinRate_POP_Ref, SDRegn_MinRate_POP_Ine)
        DEALLOCATE(SedBedDiaRCRegSegSt, SedBedDiaRCRegSegEn)
        DEALLOCATE(CEMAMFT_RandC_RegN, CEMAMFT_InCond_RegN, MFTSedFlxVars, CEMA_SD_Vars)
        DEALLOCATE(SD_NO3p2, SD_NH3p2, SD_NH3Tp2, SD_CH4p2, SD_PO4p2, SD_PO4Tp2,SD_PO4)
	      DEALLOCATE(SD_HSp2, SD_HSTp2)
	      DEALLOCATE(SD_poc2, SD_pon2, SD_pop2, SD_NH3Tp, SD_NO3p, SD_PO4Tp, SD_HSTp)
	      DEALLOCATE(SD_fpon, SD_fpoc, SD_kdiaPON, SD_ThtaPON, SD_kdiaPOC, SD_ThtaPOC)
	      DEALLOCATE(SD_JPOC, SD_JPON, SD_JPOP,SD_TDS)
        DEALLOCATE(SD_EPOC, SD_EPON, SD_EPOP)
        DEALLOCATE(SD_Denit, SD_JDenit, SD_JO2NO3,  SD_HS)   ! cb 7/26/18
	      DEALLOCATE(SD_kdiaPOP, SD_ThtaPOP, SD_NH3T, SD_FPOP)
        DEALLOCATE(SD_pHValue)  
	      DEALLOCATE(SD_AerLayerThick)
	      IF(Bubbles_Calculation) THEN
          DEALLOCATE(H2SDis, H2SGas, CH4Dis, CH4Gas)
	        DEALLOCATE(NH4Dis, NH4Gas, CO2Dis, CO2Gas)
	        DEALLOCATE(BubbleRadiusSed,PresBubbSed, PresCritSed)
	        DEALLOCATE(CgSed, C0Sed, CtSed)
	        DEALLOCATE(TConc, TConcP, SConc)
	        DEALLOCATE(DissolvedGasSediments)
	        DEALLOCATE(CrackOpen, MFTBubbReleased, LastDiffVolume)
	        DEALLOCATE(BubblesCarried, BubblesRadius)
	        DEALLOCATE(BubblesLNumber,BubblesStatus)
	        DEALLOCATE(BubblesRiseV)
          DEALLOCATE(BubbleRelWB)
	        DEALLOCATE(BubblesGasConc)
	        DEALLOCATE(BRVoluAGas, BRRateAGas)
	        DEALLOCATE(FirstBubblesRelease, BubblesReleaseAllValue)
	        DEALLOCATE(BRRateAGasNet)
	        DEALLOCATE(BubblesAtSurface)
	        DEALLOCATE(BottomTurbulence)
        END IF
    END IF

    RETURN
  End Subroutine Deallocate_CEMA