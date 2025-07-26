!===========================================================================================================================
! Removed    1) water column processes/inputs and general BOD
! Add/Modify 2) sediment-water fluxes for all water column layers 
!            3) matrix solutions of CH4
!            4) initial conditions of NO3
!            5) partitioning of PO4 for sediment layer 1
!            6) particle mixing between two layers
!            7) segment averaged outputs
! Updated 9/2020
!===========================================================================================================================  
Module CEMASedimentDiagenesis   
	Use MAIN 
  Use GLOBAL
  Use GEOMC
  Use SCREENC
  Use RSTART
  Use PREC
  Use EDDY
  Use LOGICC
  Use TVDC
  Use KINETIC
  Use SURFHE
  Use CEMAVars
  Use TRANS
  Use RSTART
    
  ! Type declarations
  IMPLICIT NONE
	!Local Variables
	Integer(2)	iTemp, RegnNum, iter, InitConRegn
	Real(R8) SD_POC_L_Fr, SD_POC_R_Fr, SD_POC_I_Fr, SD_PON_L_Fr         
  Real(R8) SD_PON_R_Fr, SD_PON_I_Fr, SD_PW_DiffCoeff, SD_Ox_Threshol      
  Real(R8) SD_POP_L_Fr, SD_POP_R_Fr, SD_POP_I_Fr
  Real(R8) SD_Ae_NH3_NO3_L, SD_Ae_NH3_NO3_H, SD_Ae_NO3_N2_L, SD_Ae_NO3_N2_H      
  Real(R8) SD_An_NO3_N2, SD_Ae_CH4_CO2, SD_Ae_HS_NH4_Nit, SD_Ae_HS_O2_Nit     
  Real(R8) SD_Theta_PW, SD_Theta_NH3_NO3, SD_Theta_NO3_N2, SD_Theta_CH4_CO2, SD_Theta_PM, SD_PartMix    
  Real(R8) SD_Sulfate_CH4_H2S, SD_Ae_H2S_SO4, SD_Theta_H2S_SO4, SD_NormConst_H2S_SO4
  Real(R8) SD_MinRate_PON_Lab, SD_MinRate_PON_Ref, SD_MinRate_PON_Ine, SD_MinRate_POC_Lab   
  Real(R8) SD_MinRate_POC_Ref, SD_MinRate_POC_Ine, SD_Theta_PON_Lab, SD_Theta_PON_Ref    
  Real(R8) SD_Theta_PON_Ine, SD_Theta_POC_Lab, SD_Theta_POC_Ref, SD_Theta_POC_Ine
  Real(R8) SD_MinRate_POP_Lab, SD_MinRate_POP_Ref, SD_MinRate_POP_Ine
  Real(R8) SD_Theta_POP_Lab, SD_Theta_POP_Ref, SD_Theta_POP_Ine
  Real(R8) SD_tc, CellThickness
  Real(R8) SD_Jcin, SD_Jpin, SD_Jnin, SD_O20, SD_JFeOOHin, SD_JMnO2in
  Real(R8) SD_Depth, SD_Tw, SD_NH30, SD_NO30, SD_TIC0, SD_ALK0    
  Real(R8) SD_PO40, SD_CH40, SD_SOD, SD_SO40, SD_H2S0
  Real(R8) SD_Fe20, SD_FeOOH0, SD_Mn20, SD_MnO20 
  Real(R8) SD_JNH4, SD_JNO3, SD_JPO4, SD_JTIC, SD_JALK, SD_JFe2, SD_JFeOOH
  Real(R8) :: SD_JMn2, SD_JT, SD_Ksw, SD_rhowcp=4.186D6, SD_tsed          ! rhowcp:4.186*1.0e6   ! units J g-1 C-1 * g m-3= J C-1 m-3
  Real(R8) :: SD_Ae_NH3_NO3, SD_Ae_NO3_N2
  Real(R8) :: SD_BEN_STRp, SD_BEN_STRp2
  Real(R8) :: SD_POCT2, SD_PONT2, SD_POPT2, SDPOCT1
  Real(R8) :: SD_H1, SD_H2, SD_JC, SD_JN, SD_JP, SD_JC1, SD_JN1, SD_JP1
  Real(R8) :: maxit, SD_es, SD_KL12, SD_SODold, SD_JSOD
  Real(R8) SD_CH4SAT, SD_s, SD_NH3toNO3, SD_fd1, SD_fp1, SD_fd2, SD_fp2
  Real(R8) SD_fdn1, SD_fpn1, SD_fdn2, SD_fpn2
  Real(R8) SD_Mn2toMnO2, SD_MnO2toMn2
  Real(R8) SD_m1, SD_m2, SD_KdNH3, SD_a11, SD_a12, SD_b1, SD_a21, SD_a22, SD_b2
  Real(R8) SD_w12, SD_w2, SD_nh3conv, SD_NSOD, SD_JDenitT, SD_JO2NO3T, SD_JC_O2equiv, SD_POCr
  Real(R8) SD_CSODmax, SD_SECH_ARG, SD_CSOD, SD_CH4toCO2
  Real(R8) SD_xappd1, SD_xappp1, SD_k1h1d, SD_k1h1p, SD_k2h2d, SD_k2h2p, SD_F12
  Real(R8) SD_F21, SD_xk1, SD_xk2, SD_KappaH2Sp1, SD_ea
  Real(R8) SD_Rho, SD_Porosity
  Real(R8) SD_SJCH4, SD_JCH4g, SD_JCH4, SD_JHS, SD_BEN_STR
  Real(R8) SD_Ammonia, SD_Ammonium, SD_SulfiMinus, SD_Sulfide
  Real(R8) MW_Constituent
  Real(R8) Dissolved_CH4_Src, Dissolved_NH3_Src, Dissolved_H2S_Src
  Real(R8) Dissolved_NO3_Src, Dissolved_SO4_Src, Dissolved_CO2_Src
  Real(R8) Dissolved_ALK_Src, Dissolved_PO4_Src, Dissolved_Fe2_Src
  Real(R8) Dissolved_Mn2_Src, Sediment_Heat_Src
  Real(R8) Dissolved_O2_Snk,  CO2ProducedSrc1L1, SO4ProducedSrc1L1
  Real(R8) LPOM_Resuspension, RPOM_Resuspension, LPOMP_Resuspension, RPOMP_Resuspension
  Real(R8) LPOMN_Resuspension,RPOMN_Resuspension
  Real(R8) SO4ConsumedSnk1L2, CO2ProducedSrc2L2, SD1_Ammonia, SD1_Ammonium
	Real(R8) SD2_Ammonia, SD2_Ammonium, SedTemp, VolWater, SedTemp1, SedTemp2
	Real(R8) AmmoniaG_SD1, AmmoniaD_SD1, AmmoniaG_SD2, AmmoniaD_SD2
	Real(R8) SD1_SulfiMinus, SD1_Sulfide, SD2_SulfiMinus, SD2_Sulfide
	Real(R8) SulfideG_SD1, SulfideD_SD1, SulfideG_SD2, SulfideD_SD2
	Real(R8) CO2ProducedCon1L1, CO2ProducedCon2L2, SD_JSO4
  Integer  SD_CH4CompMethod, SD_POMResuspMethod
	Real(R8) PW_RelRate, PW_RelRate1, PW_RelRate2
  Real(R8) SD_POM, SD_E
  REAL(R8) FETCHW, U2, COEF1, COEF2, COEF3, COEF4, HS, TS, COEF, UORB, TAU, epsilon
  REAL(R8) LW, LW0, LW1, Vscour, shields, molvisc_h2o
  REAL(R8) SD_taubot, reyn_resusp, c_bottom, c_bottom2, CO2ProducedSrc1L2, SD_Fe2toFeOOH
  REAL(R8) SD_FeOOHtoFe2, SD_Jctest, SD_JinFeOOH, SD_POCT1   !, SedimentHeat_Src
  REAL(R8) SD_Kdp1, SD_Kdp2, SD_KdFe1, SD_KdFe2, SD_KdMn1, SD_KdMn2, SD_KdNH31, SD_KdNH32, SD_KdH2S1, SD_KdH2S2
  REAL(R8) SD_delta_kpo41, SD_DOcr, SD_KsOxch, FOxch, con_cox, FNH4, vsss, fdp
  !
  ! CEMA testing variables start
  Real(R8) aO2n,o2ss,xnh4,z1ss,z2ss,z3ss,z4ss,z5ss,z5ass,z5bss,z6ss,xjn,acn,jnin,jpin
  ! CEMA testing end
  !
  ! sediment pH start
  REAL(R8) CART,ALKT,T1K,S2,SQRS2,DH1,DH2,H2CO3T,CO3T,PHT,F1,HION,HCO3T,AMMT,DH3,DHH
  REAL(R8) KW, INCR1, OH, K1, K2, bicart, H2PO4T, H3PO4T, HPO4T, HT
  REAL(R8) KAMM, KP1, KP2, KP3, NH3T, NH4T, OHT, OMCT, PHOST, PO4T
  Real(R8) t1sed,ticsed,alksed,nh4sed,po4sed,pocsed,tdssed, phsed, CO2sed, HCO3sed, CO3sed
  Real(R8) FOxna, con_nit, a12_TNH4, a21_TNH4, a22_TNH4, b2_TNH4
  Real(R8) a12_NO3, a21_NO3, a22_NO3, b2_NO3, CH42_prev
  Real(R8) a12_CH4, a21_CH4, a22_CH4, b2_CH4
    
  Real(R8) :: BIBH22, DELTABI
  Real(R8), pointer :: SD_NH31, SD_NO31, SD_PO4T1, SD_SO41, SD_TIC1, SD_pH1, SD_HST1, SD_CH41, SD_T1
  Real(R8), pointer :: SD_NH32, SD_NO32, SD_PO4T2, SD_SO42, SD_TIC2, SD_pH2, SD_HST2, SD_CH42, SD_T2
  Real(R8), pointer, Dimension(:) :: SD_POC22, SD_PON22, SD_POP22
  Real(R8), pointer :: SD_ALK1, SD_Fe2T1, SD_FeOOH1, SD_Mn2T1, SD_MnO21
  Real(R8), pointer :: SD_ALK2, SD_Fe2T2, SD_FeOOH2, SD_Mn2T2, SD_MnO22
        
  INTEGER N,ITER1
  CHARACTER(20) :: ADUMMY   ! SW 2/2019
  ! sediment pH end
  !
  real(R8), target, allocatable, dimension(:,:,:) :: C2SF     ! sediment diagenesis state variables
  real(R8),         allocatable, dimension(:,:,:) :: KFSF     ! sediment diagenesis kinetic flux
  real(R8),         allocatable, dimension(:,:)   :: KFSFAV  
  real(R8),         allocatable, dimension(:,:)   :: sdinc1,sdinn1,sdinp1
  real(R8)                                        :: sum_ave(16) 
  !
    contains
    !  INPUTS 
    !  SD_Jcin = flux to sediments from settling organic carbon 
    !         from phytoplankton and detritus in oxygen equivalent units (gO2/m2/d) 
    !         (NOTE: gO2/m2/d = gC/m2/d * 2.67 gO2/gC)                            [SDINC] in gC/m2/d from WATER QUALITY SEDIMENTC Subroutine
    !  SD_Jnin = nitrogen flux in settling phytoplankton and detritus (gN/m2/d)   [SDINN]
    !  SD_Jpin = phosphorus flux in settling phytoplankton and detritus (gP/m2/d) [SDINP]
    !  SD_O20 = dissolved oxygen in water overlying the sediment (mgO2/L) 
    !  SD_depth = total water depth overlying the sediment (m) (used to calculate methane saturation concentration at in situ pressure)

    !  SD_Tw = temperature in water overlying the sediment (deg C) 
    !  SD_NH30 = ammonia N in water overlying the sediment (mgN/L) 
    !  SD_NO30 = nitrate N in water overlying the sediment (mgN/L) 
    !  SD_PO40 = soluble reactive P in water overlying the sediment (mgP/L) 
    !  SD_CH40 = fast reacting dissolved organic carbon and CBODu in the water overlying the sediment 
    !         in oxygen equivalent units (mgO2/L) 
    !         (NOTE: mgO2/L = mC/L * 2.67 mgO2/mgC) ...error 5.33 g O2/g C (see DiToro p. 197)
    !  SD_SALw = salinity in the water overlying the sediment (ppt) 
    ! 
    !  OUTPUTS 
    !  SOD = sediment oxygen demand flux of dissolved oxygen between the water and sediment (gO2/m2/d) 
    !        (positive is loss of O2 from water column) 
    !  Jnh4 = flux of ammonia N between the water and sediment (gN/m2/d) 
    !        (positive is source of NH4-N to water column) 
    !  Jno3 = flux of nitrate N between the water and sediment (gN/m2/d) 
    !        (positive is source of NO3-N to water column) 
    !  Jch4 = flux of dissolved methane, fast reacting C, and CBODu between water and sediment in O2 equivalent units (gO2/m2/d)
    !
    !        (positive is source of CBOD to water column) 
    !        (NOTE: gO2/m2/d = gC/m2/d * 2.67 gO2/gC) 
    !        (methane is not produced in salt water) 
    !  Jch4g = flux of methane gas bubbles between the water and sediment in O2 equivalent units (gO2/m2/d) 
    !        (positive is source of CH4 bubbles to water column) 
    !        (NOTE: gO2/m2/d = gC/m2/d * 2.67 gO2/gC) 
    !        (methane is not produced in salt water) 
    !  Jhs = flux of dissolved hydrogen sulfide (COD) between water and sediment in O2 equivalent units (gO2/m2/d) 
    !        (positive is source of COD to water column) 
    !        (hydrogen sulfide is not produced in freshwater) 
    !  Jpo4 = flux of soluble reactive P between the water and sedmiment (gP/m2/d) 
    !        (positive is source of PO4-P to water column) 
    !  NH3(1) and NH3(2) = ammonia N in the sediment layers 1 and 2 (mgN/L) 
    !  NO3(1) and NO3(2) = nitrate N in the sediment layers 1 and 2 (mgN/L) 
    !  CH4(1) = dissolved methane in the aerobic sediment layer 1 (O2 equivalent units mgO2/L) 
    !  HS(1) and HS(2) = dissolved sulfide in the sediment layers 1 and 2 (O2 equivalent units mgO2/L) 
    !  PO4(1) and PO4(2) = soluble reactive P in the sediment layers 1 and 2 (mgP/L) 
  !===========================================================================================================================
  !===========================================================================================================================
  Subroutine InitCond_SedFlux
    implicit none
    character(256) :: ParameterDesc
    allocate(C2SF(KMX,IMX,37),KFSF(KMX,IMX,17),KFSFAV(IMX,16),sdinc1(KMX,IMX),sdinn1(KMX,IMX),sdinp1(KMX,IMX))
    
    KFSF=0.0;KFSFAV=0.0;C2SF=0.0;SDINC1=0.0;SDINN1=0.0;SDINP1=0.0
        
    CEMAMFT_RandC_RegN = 1
	  Do RegnNum = 1, NumRegnsSedimentDiagenesis
	    Do SegNumI = SedBedDiaRCRegSegSt(RegnNum), SedBedDiaRCRegSegEn(RegnNum)
	        CEMAMFT_RandC_RegN(SegNumI) = RegnNum
	    End Do !SegNumI
    End Do !RegnNum
        
    !Get initial condition Region Numbers for each cell
	  CEMAMFT_InCond_RegN = 1
	  Do RegnNum = 1, NumRegnsSedimentBedComposition
	    Do SegNumI = SedBedInitRegSegSt(RegnNum), SedBedInitRegSegEn(RegnNum)
	        CEMAMFT_InCond_RegN(SegNumI) = RegnNum
	    End Do !SegNumI
    End Do !RegnNum
        
    !Initialize variables
    Do JW=1, NWB
      KT = KTWB(JW)
      Do JB=BS(JW),BE(JW)
        IU = US(JB)
        ID = DS(JB)
        Do SegNumI = IU, ID
                    
		      RegnNum = CEMAMFT_RandC_RegN(SegNumI)
		      InitConRegn = CEMAMFT_InCond_RegN(SegNumI)
    		        
		      Call CEMAMFTRatesandConstants
                    
          DO LayerNum = KT, KB(SegNumI)
    		                
            IF(.NOT.RESTART_IN)THEN                           
                C2SF(LayerNum,SegNumI,1)    =  SDRegnNH3_T(InitConRegn)
                C2SF(LayerNum,SegNumI,2)    =  SDRegnNH3_T(InitConRegn)
                C2SF(LayerNum,SegNumI,3)    =  SDRegnNO3_T(InitConRegn)
                C2SF(LayerNum,SegNumI,4)    =  SDRegnNO3_T(InitConRegn)
                C2SF(LayerNum,SegNumI,5)    =  SDRegnPO4_T(InitConRegn)
                C2SF(LayerNum,SegNumI,6)    =  SDRegnPO4_T(InitConRegn)
                C2SF(LayerNum,SegNumI,7)    =  SDRegnCH4_T(InitConRegn)
                C2SF(LayerNum,SegNumI,8)    =  SDRegnCH4_T(InitConRegn)
                C2SF(LayerNum,SegNumI,9)    =  SDRegnSul_T(InitConRegn)
                C2SF(LayerNum,SegNumI,10)   =  SDRegnSul_T(InitConRegn)
                C2SF(LayerNum,SegNumI,11)   =  SDRegnH2S_T(InitConRegn)
                C2SF(LayerNum,SegNumI,12)   =  SDRegnH2S_T(InitConRegn)
                C2SF(LayerNum,SegNumI,13)   =  SDRegnTIC_T(InitConRegn)
                C2SF(LayerNum,SegNumI,14)   =  SDRegnTIC_T(InitConRegn)
                C2SF(LayerNum,SegNumI,15)   =  SDRegnT_T(InitConRegn)
                C2SF(LayerNum,SegNumI,16)   =  SDRegnT_T(InitConRegn)
                !
                IF(.NOT. IncludeDynamicpH) THEN
                  C2SF(LayerNum,SegNumI,17)  =   SDRegnpH(InitConRegn)
                  C2SF(LayerNum,SegNumI,18)  =   SDRegnpH(InitConRegn)
                ELSE
                  C2SF(LayerNum,SegNumI,17)  =   7.0
                  C2SF(LayerNum,SegNumI,18)  =   7.0
                END IF
                !
                C2SF(LayerNum,SegNumI,19)   =  SD_POC_L_Fr*SDRegnPOC_T(InitConRegn)
                C2SF(LayerNum,SegNumI,20)   =  SD_POC_R_Fr*SDRegnPOC_T(InitConRegn) 
                C2SF(LayerNum,SegNumI,21)   =  SD_POC_I_Fr*SDRegnPOC_T(InitConRegn)
                C2SF(LayerNum,SegNumI,22)   =  SD_PON_L_Fr*SDRegnPON_T(InitConRegn)
                C2SF(LayerNum,SegNumI,23)   =  SD_PON_R_Fr*SDRegnPON_T(InitConRegn) 
                C2SF(LayerNum,SegNumI,24)   =  SD_PON_I_Fr*SDRegnPON_T(InitConRegn)
                C2SF(LayerNum,SegNumI,25)   =  SD_POP_L_Fr*SDRegnPOP_T(InitConRegn)
                C2SF(LayerNum,SegNumI,26)   =  SD_POP_R_Fr*SDRegnPOP_T(InitConRegn) 
                C2SF(LayerNum,SegNumI,27)   =  SD_POP_I_Fr*SDRegnPOP_T(InitConRegn)
                !
                IF(IncludeAlkalinity) THEN
                  C2SF(LayerNum,SegNumI,28)   =  SDRegnALK_T(InitConRegn)
                  C2SF(LayerNum,SegNumI,29)   =  SDRegnALK_T(InitConRegn)
                ELSE
                  C2SF(LayerNum,SegNumI,28)   =  0.0
                  C2SF(LayerNum,SegNumI,29)   =  0.0
                END IF
                !
                IF(IncludeIron) THEN
                  C2SF(LayerNum,SegNumI,30)   =  SDRegnFe2_T(InitConRegn)
                  C2SF(LayerNum,SegNumI,31)   =  SDRegnFe2_T(InitConRegn)
                  C2SF(LayerNum,SegNumI,32)   =  SDRegnFeOOH_T(InitConRegn)
                  C2SF(LayerNum,SegNumI,33)   =  SDRegnFeOOH_T(InitConRegn)
                ELSE
                  C2SF(LayerNum,SegNumI,30)   =  0.0
                  C2SF(LayerNum,SegNumI,31)   =  0.0
                  C2SF(LayerNum,SegNumI,32)   =  0.0
                  C2SF(LayerNum,SegNumI,33)   =  0.0
                END IF
                !
                IF(IncludeManganese) THEN
                  C2SF(LayerNum,SegNumI,34)   =  SDRegnMn2_T(InitConRegn)
                  C2SF(LayerNum,SegNumI,35)   =  SDRegnMn2_T(InitConRegn)
                  C2SF(LayerNum,SegNumI,36)   =  SDRegnMnO2_T(InitConRegn)
                  C2SF(LayerNum,SegNumI,37)   =  SDRegnMnO2_T(InitConRegn)
                ELSE
                  C2SF(LayerNum,SegNumI,34)   =  0.0
                  C2SF(LayerNum,SegNumI,35)   =  0.0
                  C2SF(LayerNum,SegNumI,36)   =  0.0
                  C2SF(LayerNum,SegNumI,37)   =  0.0
                END IF
            ENDIF
                        
          END DO !Layernum
        End Do !SegNumI
      End Do !JB
    End Do !JW 
    
  End Subroutine
    
  !===========================================================================================================================
  Subroutine SedimentFlux
	  implicit none	
	
	  SD_tc = dlt*CUF/DAY     !Time step in days
	
    DO SegNumI = IU, ID
        RegnNum = CEMAMFT_RandC_RegN(SegNumI)                        
        SD_H2 = BedElevation(SegNumI)
                        
        Call CEMAMFTRatesandConstants
                        
        sum_ave = 0.0
        DO LayerNum = KT, KB(SegNumI)
          IF(LayerNum==KB(SegNumI)) THEN
            BIBH22  =  BI(LayerNum,SegNumI) / BH2(LayerNum,SegNumI)
            DELTABI =  BI(LayerNum,SegNumI) / BI(KT,SegNumI) 
          ELSE
            BIBH22  = (BI(LayerNum,SegNumI) - BI(LayerNum+1,SegNumI)) / BH2(LayerNum,SegNumI)
            DELTABI = (BI(LayerNum,SegNumI) - BI(LayerNum+1,SegNumI)) / BI(KT,SegNumI) 
          END IF
          !
          SD_Depth      = DEPTHB(LayerNum,SegNumI)
          CellThickness = H1(LayerNum,SegNumI)
          CellArea(LayerNum,SegNumI) = BI(LayerNum,SegNumI)*DLX(SegNumI)
          SD_Tsed     =   TSED(JW)
          SD_Ksw      =   CBHE(JW)
          SD_Taubot   =   SB(LayerNum, SegNumI)/B(LayerNum,SegNumI)
          SD_Rho      =   CEMASedimentDensity
          SD_Porosity =   BedPorosity(SegnumI)
          !
          vsss  = SSSO(LayerNum,SegNumI)
          fdp   = 1.0 - FPSS(LayerNum,SegNumI) - FPFE(LayerNum,SegNumI)
          !
          !Obtain properties from the water column
          SD_Tw       =   MAX(T1(LayerNum,SegNumI),     0.0)
		      SD_O20      =   MAX(C2(LayerNum,SegNumI,NDO),0.01)     
		      SD_NH30     =   MAX(C2(LayerNum,SegNumI,NNH4),0.0)       
		      SD_NO30     =   MAX(C2(LayerNum,SegNumI,NNO3),0.0)        
		      SD_PO40     =   MAX(C2(LayerNum,SegNumI,NPO4),0.0)
          SD_CH40     =   MAX(C2(LayerNum,SegNumI,NCH4),0.0)*5.33   ! SW 5/26/2022  Convert from C to O2 units
          SD_H2S0     =   MAX(C2(LayerNum,SegNumI,NH2S),0.0)*1.88   ! SW 5/26/2022 Convert from S to O2 units
          SD_SO40     =   MAX(C2(LayerNum,SegNumI,NSO4),0.0)
          SD_TIC0     =   MAX(C2(LayerNum,SegNumI,NTIC),0.0)
          !
          IF(IncludeAlkalinity) SD_ALK0   =   MAX(C2(LayerNum,SegNumI,NALK),0.0)
          !
          IF(IncludeIron) THEN
            SD_Fe20   =   MAX(C2(LayerNum,SegNumI,NFEII),0.0)
            SD_FeOOH0 =   MAX(C2(LayerNum,SegNumI,NFEOOH),0.0)
          END IF
          !
          IF(IncludeManganese) THEN
            SD_Mn20   =   MAX(C2(LayerNum,SegNumI,NMNII),0.0)    !8/2020 Corrected
            SD_MnO20  =   MAX(C2(LayerNum,SegNumI,NMNO2),0.0)
          END IF
          !
          SD_NH31  => C2SF(LayerNum,SegNumI,1);  SD_NH32  => C2SF(LayerNum,SegNumI,2) 
          SD_NO31  => C2SF(LayerNum,SegNumI,3);  SD_NO32  => C2SF(LayerNum,SegNumI,4)
          SD_PO4T1 => C2SF(LayerNum,SegNumI,5);  SD_PO4T2 => C2SF(LayerNum,SegNumI,6)
          SD_CH41  => C2SF(LayerNum,SegNumI,7);  SD_CH42  => C2SF(LayerNum,SegNumI,8)
          SD_SO41  => C2SF(LayerNum,SegNumI,9);  SD_SO42  => C2SF(LayerNum,SegNumI,10)
          SD_HST1  => C2SF(LayerNum,SegNumI,11); SD_HST2  => C2SF(LayerNum,SegNumI,12)
          SD_TIC1  => C2SF(LayerNum,SegNumI,13); SD_TIC2  => C2SF(LayerNum,SegNumI,14)
          SD_T1    => C2SF(LayerNum,SegNumI,15); SD_T2    => C2SF(LayerNum,SegNumI,16)
          SD_pH1   => C2SF(LayerNum,SegNumI,17); SD_pH2   => C2SF(LayerNum,SegNumI,18)
          SD_POC22 => C2SF(LayerNum,SegNumI,19:21)
          SD_PON22 => C2SF(LayerNum,SegNumI,22:24)
          SD_POP22 => C2SF(LayerNum,SegNumI,25:27)
          !
          IF(IncludeAlkalinity) THEN
            SD_ALK1  => C2SF(LayerNum,SegNumI,28)
            SD_ALK2  => C2SF(LayerNum,SegNumI,29)
          END IF
          !
          IF(IncludeIron) THEN
            SD_Fe2T1  => C2SF(LayerNum,SegNumI,30)
            SD_Fe2T2  => C2SF(LayerNum,SegNumI,31)
            SD_FeOOH1 => C2SF(LayerNum,SegNumI,32)
            SD_FeOOH2 => C2SF(LayerNum,SegNumI,33)
          END IF
          !
          IF(IncludeManganese) THEN
            SD_Mn2T1 => C2SF(LayerNum,SegNumI,34)
            SD_Mn2T2 => C2SF(LayerNum,SegNumI,35)
            SD_MnO21 => C2SF(LayerNum,SegNumI,36)
            SD_MnO22 => C2SF(LayerNum,SegNumI,37)
          END IF
          !
          ! Sediment POM  
          call SedimentPOM
          !
          ! Compute SOD
          call ComputeSOD
          !
          ! Sediment flux pathways
          call SedimentReaction
          !
		      !Update source/sink terms
		      DOSS(LayerNum,SegNumI)       =     DOSS(LayerNum,SegNumI)     - Dissolved_O2_Snk      !DO
          DOSEDIA(LayerNum,SegNumI)    =     -Dissolved_O2_Snk                                 ! cb 3/19/13 storing for flux output
		      NH4SS(LayerNum,SegNumI)      =     NH4SS(LayerNum,SegNumI)    + Dissolved_NH3_Src     !NH4
          NO3SS(LayerNum,SegNumI)      =     NO3SS(LayerNum,SegNumI)    + Dissolved_NO3_Src     !NO3
		      H2SSS(LayerNum,SegNumI)      =     H2SSS(LayerNum,SegNumI)    + Dissolved_H2S_Src     !H2S
		      CH4SS(LayerNum,SegNumI)      =     CH4SS(LayerNum,SegNumI)    + Dissolved_CH4_Src     !CH4
          SO4SS(LayerNum,SegNumI)      =     SO4SS(LayerNum,SegNumI)    + Dissolved_SO4_Src     !SO4
          IF(IncludeIron)      FEIISS(LayerNum,SegNumI)     =     FEIISS(LayerNum,SegNumI)   + Dissolved_Fe2_Src     !Fe(II)
          IF(IncludeManganese) MNIISS(LayerNum,SegNumI)     =     MNIISS(LayerNum,SegNumI)   + Dissolved_Mn2_Src     !Mn(II)
          TICSS(LayerNum,SegNumI)      =     TICSS(LayerNum,SegNumI)    + Dissolved_CO2_Src     !CO2
          IF(IncludeAlkalinity) ALKSS(LayerNum,SegNumI)     =     ALKSS(LayerNum,SegNumI)   + Dissolved_ALK_Src      !ALkalinity
          PO4SS(LayerNum,SegNumI)      =     PO4SS(LayerNum,SegNumI)    + Dissolved_PO4_Src     !PO4     
          SDPFLUX(JW)                  =     Dissolved_PO4_Src*sd_tc*VOL(LayerNum,SEGNUMI)/1000.  + SDPFLUX(JW)      ! SW 8/31/2017 kg
          SDNH4FLUX(JW)                =     Dissolved_NH3_Src*sd_tc*VOL(LayerNum,SEGNUMI)/1000.  + SDNH4FLUX(JW)
          SDNO3FLUX(JW)                =     Dissolved_NO3_Src*sd_tc*VOL(LayerNum,SEGNUMI)/1000.  + SDNO3FLUX(JW)
          !
          IF(CEMA_POM_Resuspension) THEN
            IF(ORGC_CALC) THEN
              LPOMCSS(LayerNum,SegNumI) =  LPOMCSS(LayerNum,SegNumI)  + LPOM_Resuspension*ORGC(JW)
              RPOMCSS(LayerNum,SegNumI) =  RPOMCSS(LayerNum,SegNumI)  + RPOM_Resuspension*ORGC(JW)
            ELSE    
              LPOMSS(LayerNum,SegNumI)  =  LPOMSS(LayerNum,SegNumI)  + LPOM_Resuspension
              RPOMSS(LayerNum,SegNumI)  =  RPOMSS(LayerNum,SegNumI)  + RPOM_Resuspension
            END IF
            LPOMPSS(LayerNum,SegNumI)   =  LPOMPSS(LayerNum,SegNumI) + LPOMP_Resuspension
            RPOMPSS(LayerNum,SegNumI)   =  RPOMPSS(LayerNum,SegNumI) + RPOMP_Resuspension
            LPOMNSS(LayerNum,SegNumI)   =  LPOMNSS(LayerNum,SegNumI) + LPOMN_Resuspension   !8/2020 Corrected
            RPOMNSS(LayerNum,SegNumI)   =  RPOMNSS(LayerNum,SegNumI) + RPOMN_Resuspension 
          END IF
         ! TSS(LayerNum,SegNumI) = TSS(LayerNum,SegNumI)   + SedimentHeat_Src     !Heat   ! Old code error
          !
          ! Output kinetic flux
          KFSF(LayerNum,SegNumI,1)  = SD_JPOC(1) + SD_JPOC(2) + SD_JPOC(3)
          KFSF(LayerNum,SegNumI,2)  = SD_JPON(1) + SD_JPON(2) + SD_JPON(3)
          KFSF(LayerNum,SegNumI,3)  = SD_JPOP(1) + SD_JPOP(2) + SD_JPOP(3)
          KFSF(LayerNum,SegNumI,4)  = SD_SOD
          KFSF(LayerNum,SegNumI,5)  = SD_CSOD
          KFSF(LayerNum,SegNumI,6)  = SD_NSOD
          KFSF(LayerNum,SegNumI,7)  = SD_JC
          KFSF(LayerNum,SegNumI,8)  = SD_JN
          KFSF(LayerNum,SegNumI,9)  = SD_JP
          KFSF(LayerNum,SegNumI,10) = SD_JCH4
          KFSF(LayerNum,SegNumI,11) = SD_JHS
          KFSF(LayerNum,SegNumI,12) = SD_JSO4
          KFSF(LayerNum,SegNumI,13) = SD_JTIC
          KFSF(LayerNum,SegNumI,14) = SD_JNH4
          KFSF(LayerNum,SegNumI,15) = SD_JNO3
          KFSF(LayerNum,SegNumI,16) = SD_JPO4
          KFSF(LayerNum,SegNumI,17) = SD_JT
          !
          sum_ave(1) = sum_ave(1)   + SD_SOD      * DELTABI
          sum_ave(2) = sum_ave(2)   + SD_POC22(1) * DELTABI
          sum_ave(3) = sum_ave(3)   + SD_POC22(2) * DELTABI
          sum_ave(4) = sum_ave(4)   + SD_JC       * DELTABI
          sum_ave(5) = sum_ave(5)   + SD_JN       * DELTABI
          sum_ave(6) = sum_ave(6)   + SD_PON22(1) * DELTABI
          sum_ave(7) = sum_ave(7)   + SD_PON22(2) * DELTABI
          sum_ave(8) = sum_ave(8)   + SD_JCH4     * DELTABI
          sum_ave(9) = sum_ave(9)   + SD_JNH4     * DELTABI
          sum_ave(10) = sum_ave(10) + SD_JNO3     * DELTABI
          sum_ave(11) = sum_ave(11) + SD_JPO4     * DELTABI
          sum_ave(12) = sum_ave(12) + SD_POP22(1) * DELTABI
          sum_ave(13) = sum_ave(13) + SD_POP22(2) * DELTABI
          sum_ave(14) = sum_ave(14) + SD_CSOD     * DELTABI
          sum_ave(15) = sum_ave(15) + SD_NSOD     * DELTABI
          sum_ave(16) = sum_ave(16) + SD_JP       * DELTABI
        END DO !LayerNum
        
        KFSFAV(SegNumI,:) = sum_ave(:)
    End Do !SegNumI
    
    If(FirstTimeInBubbles) FirstTimeInBubbles = .FALSE. 
    
  End Subroutine
    
  Subroutine SedimentPOM
    implicit none
    !
    real(R8) :: FAP3, FEP3, FCBOD3, AlgaeSed, EPBurial, CLABILE,NLABILE,PLABILE,CCC,NNN,PPP
    integer  :: JJJ, ig
    !
    SDINC1(LayerNum,SegNumi) = 0.0
    SDINN1(LayerNum,SegNumi) = 0.0
    SDINP1(LayerNum,SegNumi) = 0.0
    
    IF(DYNAMIC_SD)THEN    ! SW 1/3/2022
        CLABILE=0.0
        NLABILE=0.0
        PLABILE=0.0
        
         DO JJJ=1,NAL
      IF(ALG_CALC(JJJ))THEN
        CCC=MAX(AS(JJJ)*DAY,0.0)*AC(JJJ)*ALG(LayerNum,SegNumi,JJJ)
        NNN=MAX(AS(JJJ)*DAY,0.0)*AN(JJJ)*ALG(LayerNum,SegNumi,JJJ)
        PPP=MAX(AS(JJJ)*DAY,0.0)*AP(JJJ)*ALG(LayerNum,SegNumi,JJJ)
        CLABILE=CLABILE+CCC      ! ASSUMING 100% LABILE
        NLABILE=NLABILE+NNN
        PLABILE=PLABILE+PPP
        
        SDINC1(LayerNum,SegNumi) = SDINC1(LayerNum,SegNumi) + CCC
        SDINN1(LayerNum,SegNumi) = SDINN1(LayerNum,SegNumi) + NNN
        SDINP1(LayerNum,SegNumi) = SDINP1(LayerNum,SegNumi) + PPP      
        
      END IF
    END DO
    !
    DO JJJ=1,NEP
      IF (EPIPHYTON_CALC(JW,JJJ))THEN
        IF (LayerNum == KB(SegNumi)) THEN
          EPBurial = max(EB(JJJ)*DAY, 0.0) / H1(LayerNum,SegNumi) * EPD(LayerNum,SegNumi,JJJ) * (BI(LayerNum,SegNumi)+2.0*H1(LayerNum,SegNumi)) / BI(LayerNum,SegNumi)
        ELSE
          EPBurial = max(EB(JJJ)*DAY, 0.0) / H1(LayerNum,SegNumi) * EPD(LayerNum,SegNumi,JJJ) * (BI(LayerNum,SegNumi)-BI(LayerNum+1,SegNumi)+2.0*H1(LayerNum,SegNumi)) / (BI(LayerNum,SegNumi)-BI(LayerNum+1,SegNumi))
        END IF
        CCC=EBR(LayerNum,SegNumi,JJJ)*EPC(LayerNum,SegNumi,JJJ)*EC(JJJ)
        NNN=EBR(LayerNum,SegNumi,JJJ)*EPC(LayerNum,SegNumi,JJJ)*EN(JJJ)
        PPP=EBR(LayerNum,SegNumi,JJJ)*EPC(LayerNum,SegNumi,JJJ)*EP(JJJ)
        CLABILE=CLABILE+CCC      ! ASSUMING 100% LABILE
        NLABILE=NLABILE+NNN
        PLABILE=PLABILE+PPP

        SDINC1(LayerNum,SegNumi) = SDINC1(LayerNum,SegNumi) + CCC
        SDINN1(LayerNum,SegNumi) = SDINN1(LayerNum,SegNumi) + NNN 
        SDINP1(LayerNum,SegNumi) = SDINP1(LayerNum,SegNumi) + PPP 
      END IF
    END DO
    !
    DO JJJ=1,NBOD
     IF(CBODS(JJJ) > 0.0)THEN
      IF(BOD_CALCP(JJJ)) THEN
          PPP=MAX(CBODS(JJJ)*DAY,0.0)*RBOD(JJJ)*CBODP(LayerNum,SegNumi,JJJ)
          IF(KBOD(JJJ) > 5.8E-7)THEN                 ! IF LESS THAN 0.05 DAY-1, THEN REFRACTORY
              PLABILE=PLABILE+PPP
          ENDIF
          SDINP1(LayerNum,SegNumi) = SDINP1(LayerNum,SegNumi) + PPP
      ELSE 
          PPP=MAX(CBODS(JJJ)*DAY,0.0)*RBOD(JJJ)*BODP(JJJ)*CBOD(LayerNum,SegNumi,JJJ)
          IF(KBOD(JJJ) > 5.8E-7)THEN                 ! IF LESS THAN 0.05 DAY-1, THEN REFRACTORY
              PLABILE=PLABILE+PPP
          ENDIF
          SDINP1(LayerNum,SegNumi) = SDINP1(LayerNum,SegNumi) + PPP
      END IF    
      IF(BOD_CALCN(JJJ)) THEN
          NNN=MAX(CBODS(JJJ)*DAY,0.0)*RBOD(JJJ)*CBODN(LayerNum,SegNumi,JJJ)
          IF(KBOD(JJJ) > 5.8E-7)THEN                 ! IF LESS THAN 0.05 DAY-1, THEN REFRACTORY
              NLABILE=NLABILE+NNN
          ENDIF
          SDINN1(LayerNum,SegNumi) = SDINN1(LayerNum,SegNumi) + NNN
      ELSE    
          NNN=MAX(CBODS(JJJ)*DAY,0.0)*RBOD(JJJ)*BODN(JJJ)*CBOD(LayerNum,SegNumi,JJJ)
          IF(KBOD(JJJ) > 5.8E-7)THEN                 ! IF LESS THAN 0.05 DAY-1, THEN REFRACTORY
              NLABILE=NLABILE+NNN
          ENDIF
          SDINN1(LayerNum,SegNumi) = SDINN1(LayerNum,SegNumi) + NNN
      END IF    
      IF(BOD_CALC(JJJ)) THEN
          CCC=MAX(CBODS(JJJ)*DAY,0.0)*RBOD(JJJ)*BODC(JJJ)*CBOD(LayerNum,SegNumi,JJJ)
          IF(KBOD(JJJ) > 5.8E-7)THEN                 ! IF LESS THAN 0.05 DAY-1, THEN REFRACTORY
              CLABILE=CLABILE+CCC
          ENDIF
          SDINC1(LayerNum,SegNumi) = SDINC1(LayerNum,SegNumi) + CCC
      END IF
     ENDIF
    END DO
    !
    CLABILE=CLABILE+POMS(JW)*DAY*LPOC(LayerNum,SegNumi)
    NLABILE=NLABILE+POMS(JW)*DAY*LPON(LayerNum,SegNumi)
    PLABILE=PLABILE+POMS(JW)*DAY*LPOP(LayerNum,SegNumi)
    SDINC1(LayerNum,SegNumi) = SDINC1(LayerNum,SegNumi) + POMS(JW)*DAY*(LPOC(LayerNum,SegNumi)+RPOC(LayerNum,SegNumi))
    SDINN1(LayerNum,SegNumi) = SDINN1(LayerNum,SegNumi) + POMS(JW)*DAY*(LPON(LayerNum,SegNumi)+RPON(LayerNum,SegNumi))
    SDINP1(LayerNum,SegNumi) = SDINP1(LayerNum,SegNumi) + POMS(JW)*DAY*(LPOP(LayerNum,SegNumi)+RPOP(LayerNum,SegNumi))
    !
    IF(IncludeIron)      SD_JFeOOHin = SdinFeOOH(LayerNum,Segnumi)
    IF(IncludeManganese) SD_JMnO2in  = SdinMnO2(LayerNum,Segnumi)
    !
    SD_Jcin = SDINC1(LayerNum,SegNumi)    
    SD_Jnin = SDINN1(LayerNum,SegNumi) 
    SD_Jpin = SDINP1(LayerNum,SegNumi)
		!
    IF(IncludeIron)       SD_JinFeOOH	=  SD_JinFeOOH*DAY
    IF(IncludeManganese)  SD_JMnO2in  =  SD_JMnO2in*DAY
    !
    !gp assign constants for G class 1 and 2 PON, POC, and POP and calculate G class 3 as 1-fpox1-fpox2 
		SD_FPON(1) = NLABILE/SDINN1(LAYERNUM,SEGNUMI)         !SD_PON_L_Fr
		SD_FPON(2) = 1.0-SD_FPON(1)                           !SD_PON_R_Fr
		SD_FPON(3) = 0.0                                      !SD_PON_I_Fr 
		SD_FPOC(1) = CLABILE/SDINC1(LAYERNUM,SEGNUMI)         !SD_POC_L_Fr
		SD_FPOC(2) = 1.0-SD_FPOC(1)                           !SD_POC_R_Fr
		SD_FPOC(3) = 0.0                                      !SD_POC_I_Fr
        SD_FPOP(1) = PLABILE/SDINP1(LAYERNUM,SEGNUMI)         !SD_POP_L_Fr
		SD_FPOP(2) = 1.0-SD_FPOP(1)                           !SD_POP_R_Fr
		SD_FPOP(3) = 0.0                                      !SD_POP_I_Fr
		!
		!gp assign constants for G class 1, 2, and 3 mineralization of PON, POC, POP 
		SD_kdiaPON(1)   = SD_MinRate_PON_Lab
		SD_ThtaPON(1)   = SD_Theta_PON_Lab
		SD_kdiaPON(2)   = SD_MinRate_PON_Ref
		SD_ThtaPON(2)   = SD_Theta_PON_Ref
		SD_kdiaPON(3)   = SD_MinRate_PON_Ine
		SD_ThtaPON(3)   = SD_Theta_PON_Ine
		SD_kdiaPOC(1)   = SD_MinRate_POC_Lab
		SD_ThtaPOC(1)   = SD_Theta_POC_Lab
		SD_kdiaPOC(2)   = SD_MinRate_POC_Ref
		SD_ThtaPOC(2)   = SD_Theta_POC_Ref
		SD_kdiaPOC(3)   = SD_MinRate_POC_Ine
		SD_ThtaPOC(3)   = SD_Theta_POC_Ine
        SD_kdiaPOP(1)   = SD_MinRate_POP_Lab
		SD_ThtaPOP(1)   = SD_Theta_POP_Lab
		SD_kdiaPOP(2)   = SD_MinRate_POP_Ref
		SD_ThtaPOP(2)   = SD_Theta_POP_Ref
		SD_kdiaPOP(3)   = SD_MinRate_POP_Ine
		SD_ThtaPOP(3)   = SD_Theta_POP_Ine
		!
		!Compute input fluxes 
		Do iTemp = 1, 3 
			SD_JPOC(iTemp) = SD_Jcin * SD_FPOC(iTemp) 
			SD_JPON(iTemp) = SD_Jnin * SD_FPON(iTemp) 
			SD_JPOP(iTemp) = SD_Jpin * SD_FPOP(iTemp) 
        End Do 
        
    ELSE
        
    !
    DO JJJ=1,NAL
      IF(ALG_CALC(JJJ))THEN
        SDINC1(LayerNum,SegNumi) = SDINC1(LayerNum,SegNumi) + MAX(AS(JJJ)*DAY,0.0)*AC(JJJ)*ALG(LayerNum,SegNumi,JJJ)
        SDINN1(LayerNum,SegNumi) = SDINN1(LayerNum,SegNumi) + MAX(AS(JJJ)*DAY,0.0)*AN(JJJ)*ALG(LayerNum,SegNumi,JJJ)
        SDINP1(LayerNum,SegNumi) = SDINP1(LayerNum,SegNumi) + MAX(AS(JJJ)*DAY,0.0)*AP(JJJ)*ALG(LayerNum,SegNumi,JJJ)
      END IF
    END DO
    !
    DO JJJ=1,NEP
      IF (EPIPHYTON_CALC(JW,JJJ))THEN
        IF (LayerNum == KB(SegNumi)) THEN
          EPBurial = max(EB(JJJ)*DAY, 0.0) / H1(LayerNum,SegNumi) * EPD(LayerNum,SegNumi,JJJ) * (BI(LayerNum,SegNumi)+2.0*H1(LayerNum,SegNumi)) / BI(LayerNum,SegNumi)
        ELSE
          EPBurial = max(EB(JJJ)*DAY, 0.0) / H1(LayerNum,SegNumi) * EPD(LayerNum,SegNumi,JJJ) * (BI(LayerNum,SegNumi)-BI(LayerNum+1,SegNumi)+2.0*H1(LayerNum,SegNumi)) / (BI(LayerNum,SegNumi)-BI(LayerNum+1,SegNumi))
        END IF
        SDINC1(LayerNum,SegNumi) = SDINC1(LayerNum,SegNumi) + EBR(LayerNum,SegNumi,JJJ)*EPC(LayerNum,SegNumi,JJJ)*EC(JJJ)
        SDINN1(LayerNum,SegNumi) = SDINN1(LayerNum,SegNumi) + EBR(LayerNum,SegNumi,JJJ)*EPC(LayerNum,SegNumi,JJJ)*EN(JJJ) 
        SDINP1(LayerNum,SegNumi) = SDINP1(LayerNum,SegNumi) + EBR(LayerNum,SegNumi,JJJ)*EPC(LayerNum,SegNumi,JJJ)*EP(JJJ) 
      END IF
    END DO
    !
    DO JJJ=1,NBOD
     IF(CBODS(JJJ) > 0.0)THEN
      IF(BOD_CALCP(JJJ)) THEN
          SDINP1(LayerNum,SegNumi) = SDINP1(LayerNum,SegNumi) + MAX(CBODS(JJJ)*DAY,0.0)*RBOD(JJJ)*CBODP(LayerNum,SegNumi,JJJ)
      ELSE 
          SDINP1(LayerNum,SegNumi) = SDINP1(LayerNum,SegNumi) + MAX(CBODS(JJJ)*DAY,0.0)*RBOD(JJJ)*BODP(JJJ)*CBOD(LayerNum,SegNumi,JJJ)
      END IF    
      IF(BOD_CALCN(JJJ)) THEN
          SDINN1(LayerNum,SegNumi) = SDINN1(LayerNum,SegNumi) + MAX(CBODS(JJJ)*DAY,0.0)*RBOD(JJJ)*CBODN(LayerNum,SegNumi,JJJ)
      ELSE    
          SDINN1(LayerNum,SegNumi) = SDINN1(LayerNum,SegNumi) + MAX(CBODS(JJJ)*DAY,0.0)*RBOD(JJJ)*BODN(JJJ)*CBOD(LayerNum,SegNumi,JJJ)
      END IF    
      IF(BOD_CALC(JJJ)) THEN
          SDINC1(LayerNum,SegNumi) = SDINC1(LayerNum,SegNumi) + MAX(CBODS(JJJ)*DAY,0.0)*RBOD(JJJ)*BODC(JJJ)*CBOD(LayerNum,SegNumi,JJJ)
      END IF
     ENDIF
    END DO
    !
    SDINC1(LayerNum,SegNumi) = SDINC1(LayerNum,SegNumi) + POMS(JW)*DAY*(LPOC(LayerNum,SegNumi)+RPOC(LayerNum,SegNumi))
    SDINN1(LayerNum,SegNumi) = SDINN1(LayerNum,SegNumi) + POMS(JW)*DAY*(LPON(LayerNum,SegNumi)+RPON(LayerNum,SegNumi))
    SDINP1(LayerNum,SegNumi) = SDINP1(LayerNum,SegNumi) + POMS(JW)*DAY*(LPOP(LayerNum,SegNumi)+RPOP(LayerNum,SegNumi))
    !
    IF(IncludeIron)      SD_JFeOOHin = SdinFeOOH(LayerNum,Segnumi)
    IF(IncludeManganese) SD_JMnO2in  = SdinMnO2(LayerNum,Segnumi)
    !
    SD_Jcin = SDINC1(LayerNum,SegNumi)    
    SD_Jnin = SDINN1(LayerNum,SegNumi) 
    SD_Jpin = SDINP1(LayerNum,SegNumi)
		!
    IF(IncludeIron)       SD_JinFeOOH	=  SD_JinFeOOH*DAY
    IF(IncludeManganese)  SD_JMnO2in  =  SD_JMnO2in*DAY
    !
    !gp assign constants for G class 1 and 2 PON, POC, and POP and calculate G class 3 as 1-fpox1-fpox2 
		SD_FPON(1) = SD_PON_L_Fr
		SD_FPON(2) = SD_PON_R_Fr
		SD_FPON(3) = SD_PON_I_Fr 
		SD_FPOC(1) = SD_POC_L_Fr
		SD_FPOC(2) = SD_POC_R_Fr
		SD_FPOC(3) = SD_POC_I_Fr
    SD_FPOP(1) = SD_POP_L_Fr
		SD_FPOP(2) = SD_POP_R_Fr
		SD_FPOP(3) = SD_POP_I_Fr
		!
		!gp assign constants for G class 1, 2, and 3 mineralization of PON, POC, POP 
		SD_kdiaPON(1)   = SD_MinRate_PON_Lab
		SD_ThtaPON(1)   = SD_Theta_PON_Lab
		SD_kdiaPON(2)   = SD_MinRate_PON_Ref
		SD_ThtaPON(2)   = SD_Theta_PON_Ref
		SD_kdiaPON(3)   = SD_MinRate_PON_Ine
		SD_ThtaPON(3)   = SD_Theta_PON_Ine
		SD_kdiaPOC(1)   = SD_MinRate_POC_Lab
		SD_ThtaPOC(1)   = SD_Theta_POC_Lab
		SD_kdiaPOC(2)   = SD_MinRate_POC_Ref
		SD_ThtaPOC(2)   = SD_Theta_POC_Ref
		SD_kdiaPOC(3)   = SD_MinRate_POC_Ine
		SD_ThtaPOC(3)   = SD_Theta_POC_Ine
        SD_kdiaPOP(1)   = SD_MinRate_POP_Lab
		SD_ThtaPOP(1)   = SD_Theta_POP_Lab
		SD_kdiaPOP(2)   = SD_MinRate_POP_Ref
		SD_ThtaPOP(2)   = SD_Theta_POP_Ref
		SD_kdiaPOP(3)   = SD_MinRate_POP_Ine
		SD_ThtaPOP(3)   = SD_Theta_POP_Ine
		!
		!Compute input fluxes 
		Do iTemp = 1, 3 
			SD_JPOC(iTemp) = SD_Jcin * SD_FPOC(iTemp) 
			SD_JPON(iTemp) = SD_Jnin * SD_FPON(iTemp) 
			SD_JPOP(iTemp) = SD_Jpin * SD_FPOP(iTemp) 
        End Do 
        ENDIF
		!
    ! computing resuspended POM
    SD_POCT2 = SD_POC22(1)+SD_POC22(2)+SD_POC22(3)
    SD_POM   = SD_POCT2/ORGC(JW)
    if(CEMA_POM_Resuspension)then
      if(SD_POMResuspMethod == 0)then
        Call CEMAWindInducedSedimentResuspension
      else
        Call CEMABottomScourResuspension
      end if
      !
      SD_EPOC = 0.0
      SD_EPON = 0.0
      SD_EPOP = 0.0
      Do iTemp = 1, 3
        if( SD_POM > NONZERO .and. SD_POCT2 > Nonzero) SD_EPOC(itemp)= SD_POCT2/SD_POM * SD_POC22(itemp)/SD_POCT2 * SD_E            
        if( SD_POM > NONZERO .and. SD_PONT2 > Nonzero) SD_EPON(itemp)= SD_PONT2/SD_POM * SD_PON22(itemp)/SD_PONT2 * SD_E
        if( SD_POM > NONZERO .and. SD_POPT2 > Nonzero) SD_EPOP(itemp)= SD_POPT2/SD_POM * SD_POP22(itemp)/SD_POPT2 * SD_E            
      end do                        
    end if
    !
		!Compute particulate organic forms 
		SD_POCT2 = 0.0
		SD_PONT2 = 0.0
		SD_POPT2 = 0.0
		!Equation 13.32 DiToro. See also Equation 12.2 - also includes resuspension
		Do iTemp = 1, 3
			SD_POC22(iTemp) = (SD_POC22(iTemp) + SD_JPOC(iTemp) * SD_tc / SD_H2 - SD_EPOC(iTemp) * SD_tc / SD_H2 ) /    &
                    (1. + SD_kdiaPOC(iTemp) * SD_ThtaPOC(iTemp) ** (SD_T2 - 20.) * SD_tc + SD_W2 * SD_tc / SD_H2)
			SD_PON22(iTemp) = (SD_PON22(iTemp) + SD_JPON(iTemp) * SD_tc / SD_H2 - SD_EPON(iTemp) * SD_tc / SD_H2 ) /   &
                    (1. + SD_kdiaPON(iTemp) * SD_ThtaPON(iTemp) ** (SD_T2 - 20.) * SD_tc + SD_W2 * SD_tc / SD_H2)
			SD_POP22(iTemp) = (SD_POP22(iTemp) + SD_JPOP(iTemp) * SD_tc / SD_H2 - SD_EPOP(iTemp) * SD_tc / SD_H2 ) /    &
                    (1. + SD_kdiaPOP(iTemp) * SD_ThtaPOP(iTemp) ** (SD_T2 - 20.) * SD_tc + SD_W2 * SD_tc / SD_H2)
      if (isnan(SD_POC22(iTemp))) SD_POC22(iTemp) = 0.0
      if (isnan(SD_PON22(iTemp))) SD_PON22(iTemp) = 0.0
      if (isnan(SD_POP22(iTemp))) SD_POP22(iTemp) = 0.0
      SD_POC22(iTemp) = max(SD_POC22(iTemp), 0.0) 
      SD_PON22(iTemp) = max(SD_PON22(iTemp), 0.0) 
      SD_POP22(iTemp) = max(SD_POP22(iTemp), 0.0) 
			SD_POCT2 = SD_POCT2 + SD_POC22(iTemp)
			SD_PONT2 = SD_PONT2 + SD_PON22(iTemp)
			SD_POPT2 = SD_POPT2 + SD_POP22(iTemp)
    End Do
    !
		!Compute diagenesis fluxes 
		!Equation 13.31 diagenesis term only. See also Equation 12.2, 12.5 and 12.6
		SD_Jc = 0
		SD_Jn = 0
		SD_Jp = 0 
		Do iTemp = 1,3 
			SD_Jc = SD_Jc + SD_kdiaPOC(iTemp) * SD_ThtaPOC(iTemp)**(SD_T2 - 20.)*SD_POC22(iTemp)*SD_H2 
			SD_Jn = SD_Jn + SD_kdiaPON(iTemp) * SD_ThtaPON(iTemp)**(SD_T2 - 20.)*SD_PON22(iTemp)*SD_H2 
			SD_Jp = SD_Jp + SD_kdiaPOP(iTemp) * SD_ThtaPOP(iTemp)**(SD_T2 - 20.)*SD_POP22(iTemp)*SD_H2 
    End Do 
    
  End Subroutine
    
  Subroutine CEMAMFTRatesandConstants    
    SD_POC_L_Fr             =     SDRegnPOC_L_Fr(RegnNum)
    SD_POC_R_Fr             =     SDRegnPOC_R_Fr(RegnNum)
    SD_POC_I_Fr             =     1 - SD_POC_L_Fr - SD_POC_R_Fr
    SD_PON_L_Fr             =     SDRegnPON_L_Fr(RegnNum)
    SD_PON_R_Fr             =     SDRegnPON_R_Fr(RegnNum)
    SD_PON_I_Fr             =     1 - SD_PON_L_Fr - SD_PON_R_Fr
    SD_POP_L_Fr             =     SDRegnPOP_L_Fr(RegnNum)
    SD_POP_R_Fr             =     SDRegnPOP_R_Fr(RegnNum)
    SD_POP_I_Fr             =     1 - SD_POP_L_Fr - SD_POP_R_Fr
    SD_PW_DiffCoeff         =     SDRegnPW_DiffCoeff(RegnNum)    ! m^2/d
    SD_Ox_Threshol          =     SDRegnOx_Threshold(RegnNum)
    SD_Ae_NH3_NO3_L         =     SDRegnAe_NH3_NO3_L(RegnNum)
    SD_Ae_NH3_NO3_H         =     SDRegnAe_NH3_NO3_H(RegnNum)
    SD_Ae_NO3_N2_L          =     SDRegnAe_NO3_N2_L(RegnNum)
    SD_Ae_NO3_N2_H          =     SDRegnAe_NO3_N2_H(RegnNum)
    SD_An_NO3_N2            =     SDRegnAn_NO3_N2(RegnNum)
    SD_Ae_CH4_CO2           =     SDRegnAe_CH4_CO2(RegnNum)
    SD_KsOxch               =     KsOxch(RegnNum)
    SD_Ae_HS_NH4_Nit        =     SDRegnAe_HS_NH4_Nit(RegnNum)
    SD_Ae_HS_O2_Nit         =     SDRegnAe_HS_O2_Nit(RegnNum)
    SD_Theta_PW             =     SDRegn_Theta_PW(RegnNum)
    SD_Theta_PM             =     SDRegn_Theta_PM(RegnNum)
    SD_Theta_NH3_NO3        =     SDRegn_Theta_NH3_NO3(RegnNum)   
    SD_Theta_NO3_N2         =     SDRegn_Theta_NO3_N2(RegnNum)  
    SD_Theta_CH4_CO2        =     SDRegn_Theta_CH4_CO2(RegnNum)
    SD_Sulfate_CH4_H2S      =     SDRegn_Sulfate_CH4_H2S(RegnNum)
    SD_Ae_H2S_SO4           =     SDRegnAe_H2S_SO4(RegnNum)
    SD_Theta_H2S_SO4        =     SDRegn_Theta_H2S_SO4(RegnNum)
    SD_NormConst_H2S_SO4    =     SDRegn_NormConst_H2S_SO4(RegnNum)
    SD_MinRate_PON_Lab      =     SDRegn_MinRate_PON_Lab(RegnNum)
    SD_MinRate_PON_Ref      =     SDRegn_MinRate_PON_Ref(RegnNum)
    SD_MinRate_PON_Ine      =     SDRegn_MinRate_PON_Ine(RegnNum)
    SD_MinRate_POC_Lab      =     SDRegn_MinRate_POC_Lab(RegnNum)
    SD_MinRate_POC_Ref      =     SDRegn_MinRate_POC_Ref(RegnNum)
    SD_MinRate_POC_Ine      =     SDRegn_MinRate_POC_Ine(RegnNum)
    SD_MinRate_POP_Lab      =     SDRegn_MinRate_POP_Lab(RegnNum)
    SD_MinRate_POP_Ref      =     SDRegn_MinRate_POP_Ref(RegnNum)
    SD_MinRate_POP_Ine      =     SDRegn_MinRate_POP_Ine(RegnNum)
    SD_Theta_PON_Lab        =     SDRegn_Theta_PON_Lab(RegnNum)
    SD_Theta_PON_Ref        =     SDRegn_Theta_PON_Ref(RegnNum)
    SD_Theta_PON_Ine        =     SDRegn_Theta_PON_Ine(RegnNum)
    SD_Theta_POC_Lab        =     SDRegn_Theta_POC_Lab(RegnNum)
    SD_Theta_POC_Ref        =     SDRegn_Theta_POC_Ref(RegnNum)
    SD_Theta_POC_Ine        =     SDRegn_Theta_POC_Ine(RegnNum)
    SD_Theta_POP_Lab        =     SDRegn_Theta_POP_Lab(RegnNum)
    SD_Theta_POP_Ref        =     SDRegn_Theta_POP_Ref(RegnNum)
    SD_Theta_POP_Ine        =     SDRegn_Theta_POP_Ine(RegnNum)
    SD_CH4CompMethod        =     SDRegn_CH4CompMethod(RegnNum)
    SD_POMResuspMethod      =     SDRegn_POMResuspMethod(RegnNum)
    SD_Kdp2                 =     Kdp2(RegnNum)
    SD_DOcr                 =     DOcr(RegnNum)
    SD_delta_kpo41          =     delta_kpo41(RegnNum)
    SD_KdNH31               =     KdNH31(RegnNum)
    SD_KdNH32               =     KdNH32(RegnNum)
    SD_KdH2S1               =     KdH2S1(RegnNum)
    SD_KdH2S2               =     KdH2S2(RegnNum)
    SD_PartMix              =     PartMixVel(RegnNum)
    SD_w2                   =     BurialVel(RegnNum)
    SD_POCr                 =     POCr(RegnNum)
    !
    IF(IncludeIron) THEN
        SD_KdFe1            =     KdFe1(RegnNum)   
        SD_KdFe2            =     KdFe2(RegnNum)
    END IF
    IF(IncludeManganese) THEN
        SD_KdMn1            =     KdMn1(RegnNum)   
        SD_KdMn2            =     KdMn2(RegnNum)
    END IF 
    
  End Subroutine
     
  !===========================================================================================================================
  ! Compute benthic sediment oxygen demand
  Subroutine ComputeSOD
    USE SCREENC, ONLY:JDAY
    integer  :: it
    !
    If(SD_O20 > SD_Ox_Threshol)Then
        SD_Ae_NH3_NO3   =   SD_Ae_NH3_NO3_H
        SD_Ae_NO3_N2    =   SD_Ae_NO3_N2_H
    Else
        SD_Ae_NH3_NO3   =   SD_Ae_NH3_NO3_L
        SD_Ae_NO3_N2    =   SD_Ae_NO3_N2_L
    End If
		!
		!
		SD_SOD = SD_Jc * 32.0 / 12.0 + 1.714 * SD_Jn
        !
		!Calculation of final benthic particle mixing velocity Equation 13.6 in DoToro (2001)
		!See also Equation 2.56 and 4.48
		SD_KL12 = SD_PW_DiffCoeff * (SD_Theta_PW ** (SD_T2-20.)) / (SD_H2/2.)
    if (isnan(SD_KL12)) SD_KL12 = 0.0
    !
    SD_W12  = SD_PartMix * (SD_Theta_PM ** (SD_T2-20.)) / (SD_H2/2.) * SD_POC22(1) / (SD_POCr*SD_Rho*(1.0-SD_Porosity))  
		if (isnan(SD_W12)) SD_W12 = 0.0
    !
    SD_S   = SD_SOD / SD_O20
    if (isnan(SD_S) .or. SD_S == 0.0) SD_S = 1.0E-8
    !
    SD_H1 = SD_KL12 * SD_H2 / SD_s
		If(SD_H1 > SD_H2)Then   
        SD_H1 = SD_H2
		End If
		SD_AerLayerThick(SegNumI) = SD_H1
    !
    !
    !CH4
    !SD_CH4SAT = 100.0D+00*(1.0D+00 + SD_depth/10.0D+00)*(1.024**(20.0D+00 - SD_T2))        ![gmO*/m3]
    SD_CH4SAT = 100.0D+00*(1.0D+00 + SD_depth/10.0D+00)*(1.024**(20.0D+00 - SD_Tw))     !Saturation conc. of methane in oxygen equivalent units (Equation 10.51) [gmO*/m3]  
    If(SD_CH4CompMethod == 1) CH42_prev = SD_CH42
    !------------------------------------------------------------------------------------------------------------------
    ! compute SOD
		maxit = 500   !1000 
		SD_es = 0.01  !0.001
    !
    DO it = 1, maxit
      !
      CO2ProducedSrc1L1 = 0.d0
			SO4ProducedSrc1L1 = 0.d0
			SO4ConsumedSnk1L2 = 0.d0
			CO2ProducedSrc2L2 = 0.d0
      !
      !-----------------------------------------------------------------------------------------------------------------------
      ! TNH41 and TNH42   !Calculate dissolved and particulate (sorbed) fractions 
      SD_fdn1 = 1.0/(1.0 + SD_KdNH31*SD_Rho*(1.0-SD_Porosity))
      SD_fpn1 = 1.0 - SD_fdn1             != ((m1*KdNH3)/(1 + m1*KdNH3)) 
      SD_fdn2 = 1.0/(1.0 + SD_KdNH32*SD_Rho*(1.0-SD_Porosity))
      SD_fpn2 = 1.0 - SD_fdn2             != ((m2*KdNH3)/(1 + m2*KdNH3)) 
      SD_NH3T(1) = SD_NH31/sd_fdn1
      SD_NH3T(2) = SD_NH32/sd_fdn2
      FOxna    = SD_O20 / (SD_Ae_HS_O2_Nit * 2.0 + SD_O20) 
      if (isnan(FOxna)) FOxna = 0.0
      con_nit  = ((SD_Ae_NH3_NO3*(SD_Theta_NH3_NO3**(SD_T1-20.)))**2.0)* FOxna * sd_fdn1
      !
      a12_TNH4 = SD_fdn2*SD_KL12 + SD_fpn2*SD_w12
      a21_TNH4 = SD_fdn1*SD_KL12 + SD_fpn1*SD_w12 + SD_w2
      a22_TNH4 = -SD_fdn2 * SD_KL12 - SD_fpn2 * SD_w12 - SD_w2 - SD_H2 / SD_tc
      b2_TNH4  = -SD_Jn - SD_H2 / SD_tc * SD_NH3T(2)
      !Equation 4.51
      ! 
      IF(SD_Ae_HS_NH4_Nit > 0.0) THEN
          FNH4 = SD_Ae_HS_NH4_Nit/(SD_Ae_HS_NH4_Nit + SD_NH3T(1)*SD_fdn1)
      ELSE 
          FNH4 = 1.0
      END IF
      !
      SD_a11 = -SD_fdn1*SD_KL12 - SD_fpn1*SD_w12 - SD_w2- con_nit * FNH4 / SD_s - SD_fdn1*SD_s  
			SD_b1  = -SD_s*SD_NH30             ![m/d]*[mg/m3]
                   
      Call Lin_Sys(SD_a11, a12_TNH4,a21_TNH4, a22_TNH4, SD_b1,b2_TNH4, SD_NH3T(1), SD_NH3T(2)) 
			SD_NH3T(1) = max(SD_NH3T(1),0.0)
      SD_NH3T(2) = max(SD_NH3T(2),0.0)
      !
      ! NO31 and NO32
      a12_NO3 = SD_KL12 
      a21_NO3 = SD_KL12    
      a22_NO3 = -SD_KL12 - SD_An_NO3_N2*(SD_Theta_NO3_N2**(SD_T2-20.)) - SD_H2 / SD_tc
      b2_NO3  = -SD_H2 / SD_tc * SD_NO32
      SD_a11  = -SD_KL12 - ((SD_Ae_NO3_N2*(SD_Theta_NO3_N2**(SD_T1-20.)))**2.0)/SD_s - SD_s 
			SD_b1   = -SD_s*SD_NO30 - con_nit * FNH4/SD_s*SD_NH3T(1)
      Call Lin_Sys(SD_a11, a12_NO3, a21_NO3, a22_NO3, SD_b1, b2_NO3, SD_NO31, SD_NO32)
      SD_NO31 = max(SD_NO31,0.0)
      SD_NO32 = max(SD_NO32,0.0)
      !
      !Denitrification in layers 1 and 2 (Equation 4.55) 
			SD_Denit(1) = ((SD_Ae_NO3_N2*SD_Theta_NO3_N2**(SD_T1-20.))**2.0)/SD_s 
			SD_Denit(2) = SD_An_NO3_N2*SD_Theta_NO3_N2**(SD_T2-20.) 
      !Denitrification Flux [mgN/m2d] 
			SD_JDenit(1) = SD_Denit(1) * SD_NO31 
			SD_JDenit(2) = SD_Denit(2) * SD_NO32 
			SD_JDenitT   = SD_JDenit(1) + SD_JDenit(2) 
			!    
			!Methane consumption due to denitrification (Equation 9.16) 
			SD_JO2NO3(1) = (32.0D+00 / 12.0D+00) * (10.0D+00 / 8.0D+00) * (12.0D+00 / 14.0D+00) * SD_JDenit(1) 
			SD_JO2NO3(2) = (32.0D+00 / 12.0D+00) * (10.0D+00 / 8.0D+00) * (12.0D+00 / 14.0D+00) * SD_JDenit(2) 
			!
			!Sum 
			SD_JO2NO3T = SD_JO2NO3(1) + SD_JO2NO3(2) 
			!
			!Calculate methane flux in oxygen equivalent units, adjusted for 
			!the methane consumed in denitrification                            
			!gp also used if sulfide is produced 
			SD_JC_O2equiv = SD_Jc * 32.0 / 12.0 - SD_JO2NO3T 
			SD_JC_O2equiv = max(SD_JC_O2equiv,1.0E-10)
      !
      ! CH41 and CH42
      If (SD_SO42 <= SD_Sulfate_CH4_H2S) Then   
				!gp freshwater methane production, no changes to original code 
				!CSODMAX Equations 10.28 and 10.30 
				!SD_CSODmax = DMin1((2.0D+00 * SD_KL12 * SD_CH4SAT * SD_JC_O2equiv)**2.0D+00, SD_JC_O2equiv)    ![gmO*/m2-d] = sqr([m/d] * [gmO*/m3] * [gmO*/m2-d])   ! SW 10/10/2017 MAJOR ERROR
        SD_CSODmax = DMin1((2.0D+00 * SD_KL12 * SD_CH4SAT * SD_JC_O2equiv)**0.5D+00, SD_JC_O2equiv)    ![gmO*/m2-d] = sqr([m/d] * [gmO*/m3] * [gmO*/m2-d])   ! SW 10/10/2017
				If(SD_CH4CompMethod == 0) Then
					!***********************************************************************
					!Analytical solution for methane
					SD_SECH_ARG = (SD_Ae_CH4_CO2 * SD_Theta_CH4_CO2 ** ((SD_T1-20.) / 2.0)) / SD_s
					!CSOD Equation 10.35
					!The hyperbolic secant is defined as HSec(X) = 2 / (Exp(X) + Exp(-X))
					If (SD_SECH_ARG < 400.0) Then !This is the usual case
						SD_CSOD = SD_CSODmax * (1.0 - (2.0 / (Exp(SD_SECH_ARG) + Exp(-SD_SECH_ARG))))
					Else !HSec(SECH_ARG) < 3.8E-174 ~ 0
						SD_CSOD = SD_CSODmax
					End If
					!***********************************************************************
        Else If(SD_CH4CompMethod == 1) Then	
					!
					!NumericalSolution for CH4
					!SD_CH4toCO2 = (SD_Ae_CH4_CO2 ** 2.0D+00 * SD_Theta_CH4_CO2**((SD_T(1) - 20.0D+00) / 2.0D+00)) / SD_s 
					!SD_CH4(1) = (SD_CSODmax + SD_s * SD_CH40) / (SD_CH4toCO2 + SD_s) 
					!SD_CSOD = SD_CH4toCO2 * SD_CH4(1) 
          !  
          ! CH41 and CH42 
          !CH42_prev = SD_CH42
          FOxch = SD_O20 / (SD_O20 + SD_KsOxch * 2.0)
          if (isnan(Foxch)) FOxch = 0.0
          con_cox = ((SD_Ae_CH4_CO2 * (SD_Theta_CH4_CO2**(SD_T1-20.)))**2.0) * FOxch
          a12_CH4 = SD_KL12
          a21_CH4 = SD_KL12
          a22_CH4 = -SD_KL12 - SD_H2 / SD_tc
          SD_a11 = -SD_KL12 - con_cox / SD_S - SD_S
          SD_b1  = -SD_S * SD_CH40
          SD_b2  = -SD_JC_O2equiv - CH42_prev * SD_H2 / SD_tc
                    
          !write(199,*) jday,SD_Ae_CH4_CO2,SD_Theta_CH4_CO2
          !if(segnumi==23 .and. layernum==kb(segnumi)) then
          !    write(199,'(11F10.4)') jday, sd_a11, a12_CH4, a21_CH4, a22_CH4,sd_b1, SD_b2,SD_JC_O2equiv,CH42_prev,SD_H2,SD_tc
          !end if
                    
          Call Lin_Sys(SD_a11, a12_CH4, a21_CH4, a22_CH4, SD_b1, SD_b2, SD_CH41, SD_CH42)
          SD_CH41 = max(SD_CH41,0.0)
          SD_CH42 = max(SD_CH42,0.0)
          IF (SD_CH42 > SD_CH4SAT) THEN
              SD_CH42 = SD_CH4SAT
              SD_CH41 = (SD_b1 - a12_CH4 * SD_CH42) / SD_a11
          END IF
          !if(segnumi==23 .and. layernum==kb(segnumi)) then
          !    write(199,'(2F10.4)') jday, SD_CH4SAT
          !end if
          SD_CH41  = max(SD_CH41, 0.0)
          SD_CH42  = max(SD_CH42, 0.0)
          SD_CSOD    = con_cox / SD_S * SD_CH41
				End IF
				
				!0.5CH4 + O2 --> 0.5 CO2 + H2O
				!CO2 Produced = 0.5*(12+32)/32 = 0.6875
				CO2ProducedSrc1L1 = SD_CSOD*0.6875    !g CO2/m�/d
				CO2ProducedCon1L1 = CO2ProducedSrc1L1*SD_tc/SD_H1   !g CO2/m�/d*d/m = g CO2/m�
				
			Else 
				!
				!gp saltwater sulfide production by C diagenesis based on DiToro (2001) Appendix B 
				!***** Calculate dissolved and particulate (sorbed) fractions for sulfide 
        SD_fd1=1.0/(1.0+SD_KdH2S1*SD_Rho*(1.0-SD_Porosity))
        !SD_fp1=(SD_KdH2S1*SD_Rho*(1.0-SD_Porosity))/(SD_Porosity+SD_KdH2S1*SD_Rho*(1.0-SD_Porosity))
        SD_fp1=1.0-SD_fd1                    != ((m1*KdH2S1)/(1 + m1*KdH2S1))
        SD_fd2=1.0/(1.0+SD_KdH2S2*SD_Rho*(1.0-SD_Porosity))
        !SD_fp2=(SD_KdH2S2*SD_Rho*(1.0-SD_Porosity))/(SD_Porosity+SD_KdH2S2*SD_Rho*(1.0-SD_Porosity)) 
        SD_fp2=1.0-SD_fd2                    != ((m2*KdH2S2)/(1 + m2*KdH2S2)) 
				!
				!***** Temperature adjusted reaction velocities 
				SD_xappd1 = SD_Ae_H2S_SO4 * SD_Theta_H2S_SO4 **((SD_T1-20.) / 2.0D+00) 
				SD_xappp1 = SD_KappaH2Sp1 * SD_Theta_H2S_SO4 **((SD_T1-20.) / 2.0D+00) 
				!
				!***** Transport and Decay terms 
				!Equation B.19
				SD_k1h1d = SD_xappd1**2.0D+00/SD_s*(SD_O20/SD_NormConst_H2S_SO4) + SD_s 
				SD_k1h1p = SD_xappp1**2.0D+00/SD_s*(SD_O20/SD_NormConst_H2S_SO4) 
				SD_k2h2d = 0.0d+00 
				SD_k2h2p = 0.0d+00
				SD_F12 = SD_w12 * SD_fp1 + SD_KL12 * SD_fd1 
				SD_F21 = SD_w12 * SD_fp2 + SD_KL12 * SD_fd2 
				SD_xk1 = SD_k1h1d * SD_fd1 + SD_k1h1p * SD_fp1 
				SD_xk2 = SD_k2h2d * SD_fd2 + SD_k2h2p * SD_fp2 
				!
				!***** Matrix and forcing function 
				SD_a11 = -SD_F12 - SD_xk1 - SD_w2                             !note: -fd1 * s is included in DiToro's -xk1 term 
        !SD_a11 = -SD_H1 / SD_tc - SD_F12 - SD_xk1 - SD_w2            !note: -fd1 * s is included in DiToro's -xk1 term 
				SD_a21 = SD_F12 + SD_w2 
				SD_a12 = SD_F21 
				SD_b1 = 0.0D+00 
				SD_a22 = -SD_F21 - SD_xk2 - SD_w2 - SD_H2 / SD_tc
				SD_b2  = -SD_JC_O2equiv - SD_H2 / SD_tc * SD_HST2
				
				Call Lin_Sys(SD_a11, SD_a12, SD_a21, SD_a22, SD_b1, SD_b2, SD_HST1, SD_HST2)   !, NFLog, NFCle) 
				SD_HST1 = max(SD_HST1,0.0)
        SD_HST2 = max(SD_HST2,0.0)
        !
				!***** dissolved concentrations 
				SD_HS(1) = SD_fd1*SD_HST1 
				SD_HS(2) = SD_fd2*SD_HST2 
				SD_CSOD = (SD_xappd1**2.0/SD_s*SD_fd1 + SD_xappp1**2.0/SD_s*SD_fp1)*(SD_O20/SD_NormConst_H2S_SO4)*SD_HST1 
				
				!H2S + 2 O2 --> 2 H+ + SO42-
				!SO42- Produced = (32+16*4)/(2*32) = 1.5
				SO4ProducedSrc1L1 = SD_CSOD*1.5    !g SO42-/m�/d
				
				!CH2O + 2 H+ + SO42- --> 2 CO2 + H2S + 2 H2O    !CH2O is represented in O2 equivalent in SD_JC_O2equiv
				!SO42- Consumed = (32+16*4)/(2*16) = 3.0
				!SO4ConsumedSnk1L2 = SD_JC_O2equiv*3.0    !g SO42-/m�/d
        SO4ConsumedSnk1L2 = -SD_JC_O2equiv*3.0    !g SO42-/m�/d   ! cb 7/26/18
				
				!CH2O + 2 H+ + SO42- --> 2 CO2 + H2S + 2 H2O    !CH2O is represented in O2 equivalent in SD_JC_O2equiv
				!CO2 Produced = 2*(12+16*2)/(2*16) = 2.75
				CO2ProducedSrc2L2 = SD_JC_O2equiv*2.75    !g CO2/m�/d
				CO2ProducedCon2L2 = CO2ProducedSrc2L2*SD_tc/SD_H2   !g CO2/m�/d*d/m = g CO2/m�
      End If
      !
      ! Metals Start
      if(includeIron)then
        SD_fd1 = 1.0/(1.0 + SD_KdFe1*SD_Rho*(1.0-SD_Porosity))
        SD_fp1 = 1.0 - SD_fd1
        SD_fd2 = 1.0/(1.0+SD_KdFe2*SD_Rho*(1.0-SD_Porosity))          
        SD_fp2 = 1.0 - SD_fd2
        SD_Fe2toFeOOH = kfe_oxid(JW)*SD_O20*10**(2.0*(sd_ph1-7.0))*SD_fd1*SD_Fe2T1
        SD_CSOD = SD_CSOD + (0.25*2.0*16.0/55.845)*SD_Fe2toFeOOH  ! lumping in DO consumed by Fe(II)>FeOOH into CSOD
        SD_FeOOHtoFe2 = kfe_red(JW)*SD_FeOOH1
        !FeOOH + 0.25CH2O + 2H+ > Fe(II) + 0.25CO2 + 1.75H20
        ! CO2 Produced = 0.25 * (12 + 2*16) / 55.845 = 0.197         !g CO2/m?d
        CO2ProducedSrc1L2 = CO2ProducedSrc1L2 + SD_FeOOHtoFe2 * 0.197
        !Write linear system of equations around total ferrous iron SD_Fe2T
			  !Equation 5.1
			  !
			  !Layer 1
			  SD_a11 = -SD_fd1*SD_KL12 - SD_fp1*SD_w12 - SD_fd1*SD_s - SD_w2 - SD_Fe2toFeOOH*SD_H1
			  SD_a12 = SD_fd2*SD_KL12 + SD_fp2*SD_w12 
			  SD_b1 = -SD_s*SD_Fe20
			  !
			  !Layer 2
			  SD_a21 = SD_fd1*SD_KL12 + SD_fp1*SD_w12 + SD_w2   
			  SD_a22 = -SD_fd2 * SD_KL12 - SD_fp2 * SD_w12 - SD_w2 - SD_H2 / SD_tc
			  SD_b2  = - SD_H2 / SD_tc * SD_Fe2T2 - SD_H2 * SD_FeOOHtoFe2
			!			
			  Call Lin_Sys(SD_a11, SD_a12, SD_a21, SD_a22, SD_b1, SD_b2, SD_Fe2T1, SD_Fe2T2) !, NFLog, NFCle) 
        SD_Fe2T1 = max(SD_Fe2T1,0.0)
        SD_Fe2T2 = max(SD_Fe2T2,0.0)
        !Write linear system of equations around total ferrous iron SD_FeOOH
			  !Equation 5.1
			  !
			  !Layer 1
			  SD_a11 = -SD_w12 -  SD_w2 
			  SD_a12 = SD_w12 
			  SD_b1 = -SD_JFeOOHin- SD_Fe2toFeOOH*SD_H1
			  !
			  !Layer 2
			  SD_a21 = SD_w12 + SD_w2    
			  SD_a22 = - SD_w12 - SD_w2 - SD_H2*SD_FeOOHtoFe2  - SD_H2 / SD_tc
			  SD_b2  = - SD_H2 / SD_tc * SD_FeOOH2
			!			
			  Call Lin_Sys(SD_a11, SD_a12, SD_a21, SD_a22, SD_b1, SD_b2, SD_FeOOH2, SD_FeOOH2)  !, NFLog, NFCle)                         		
        SD_FeOOH1 = max(SD_FeOOH1,0.0)
        SD_FeOOH2 = max(SD_FeOOH2,0.0)
      end if
      !
      if(includeManganese)then
        !calculating dissolved and particulate forms of Mn(II) (Chapra, eqn. 25.89)
        SD_fd1 = 1.0/(1.0 + SD_KdMn1*SD_Rho*(1.0-SD_Porosity))
        SD_fp1 = 1.0 - SD_fd1
        SD_fd2 = 1.0/(1.0 + SD_KdMn2*SD_Rho*(1.0-SD_Porosity))
        SD_fp2 = 1.0 - SD_fd2
        SD_Mn2toMnO2 = kMn_oxid(JW)*SD_O20*10**(2.0*(sd_ph1-7.0))*SD_fd1*SD_Mn2T1
        SD_CSOD = SD_CSOD + (16.0/54.94)*SD_Mn2toMnO2  ! lumping in DO consumed by Mn(II)>MnO2 into CSOD
        SD_MnO2toMn2 = kMn_red(JW)*SD_MnO21
        !MnO2 + 0.5CH2O + 2H+ > Mn(II) + 0.5CO2 + 1.5H20
        ! CO2 Produced = 0.5 * (12 + 2*16) / 54.94 = 0.400         !g CO2/m?d
        CO2ProducedSrc1L2 = CO2ProducedSrc1L2 + SD_MnO2toMn2 * 0.400
        !Write linear system of equations around total Mn(II) SD_Mn2T
			  !Equation 5.1
			  !
			  !Layer 1
			  SD_a11 = -SD_fd1*SD_KL12 - SD_fp1*SD_w12 - SD_fd1*SD_s - SD_w2 - SD_Mn2toMnO2*SD_H1
			  SD_a12 = SD_fd2*SD_KL12 + SD_fp2*SD_w12 
			  SD_b1 = -SD_s*SD_Mn20
			  !
			  !Layer 2
			  SD_a21 = SD_fd1*SD_KL12 + SD_fp1*SD_w12 + SD_w2 
			  SD_a22 = -SD_fd2 * SD_KL12 - SD_fp2 * SD_w12 - SD_w2 - SD_H2 / SD_tc
			  SD_b2  = - SD_H2 / SD_tc * SD_Mn2T2 - SD_H2 * SD_MnO2toMn2
			!			
			  Call Lin_Sys(SD_a11, SD_a12, SD_a21, SD_a22, SD_b1, SD_b2, SD_Mn2T1, SD_Mn2T2)   !, NFLog, NFCle) 
        SD_Mn2T1 = max(SD_Mn2T1,0.0)
        SD_Mn2T2 = max(SD_Mn2T2,0.0)        
        !
        !Write linear system of equations around manganese dioxide SD_MnO2
			  !Equation 5.1
			  !
			  !Layer 1
			  SD_a11 = -SD_w12 -  SD_w2 
			  SD_a12 = SD_w12 
			  SD_b1  = -SD_JMnO2in- SD_Mn2toMnO2*SD_H1
			  !
			  !Layer 2
			  SD_a21 = SD_w12 + SD_w2
			  SD_a22 = - SD_w12 - SD_w2 - SD_H2*SD_MnO2toMn2  - SD_H2 / SD_tc
			  SD_b2  = - SD_H2 / SD_tc * SD_MnO22
			!			
			  Call Lin_Sys(SD_a11, SD_a12, SD_a21, SD_a22, SD_b1, SD_b2, SD_MnO21, SD_MnO22)   !, NFLog, NFCle)                         		
        SD_MnO21 = max(SD_MnO21,0.0)
        SD_MnO22 = max(SD_MnO22,0.0)
      end if
      !
      ! new SOD
      SD_NSOD = 2.0*32.0/14.0 * con_nit * FNH4 / SD_S * SD_NH3T(1)
      SD_SODold = SD_SOD
			SD_SOD = (SD_SOD + SD_CSOD + SD_NSOD)/2.0  
			SD_ea = Abs((SD_SOD - SD_SODold)/SD_SOD)*100.0D+00 
			If (SD_ea <= SD_es) Exit 
      SD_s = SD_SOD/SD_O20
      if (isnan(SD_S) .or. SD_S == 0.0) SD_S = 1.0E-8

      If (it >= maxit-1) THEN
          IF(IT==MAXIT-1)THEN
              Write(WRN,'(A,f12.3,A,f10.6,A,f10.6,A,I3,A,I3)') 'Sediment Diagenesis: SOD iterations almost exceeded on JDAY:',JDAY,' SD_SOD:',SD_SOD,' SD_SODold:',SD_SODold,' LayerNum:',LayerNum,' SegNumI:',SegNumI
          ELSE
              Write(WRN,'(A,f12.3,A,f10.6,A,f10.6,A,I3,A,I3)') 'Sediment Diagenesis: SOD iterations exceeded on JDAY:',JDAY,' SD_SOD:',SD_SOD,' SD_SODold:',SD_SODold,' LayerNum:',LayerNum,' SegNumI:',SegNumI
          ENDIF
      ENDIF

    End Do
    !
    ! check if SOD solution is converged
	!	If (it >= maxit) Write(WRN,*) 'Sediment Diangenesis: SOD iterations exceeded on JDAY:',JDAY
  End Subroutine
    
  Subroutine SedimentReaction
    !
    SD_JSOD = SD_SOD
    SD_s = SD_SOD/SD_O20
    !
    SD_H1 = SD_KL12 * SD_H2 / SD_s
    If(SD_H1 > SD_H2)Then   
        SD_H1 = SD_H2
    End If
    SD_AerLayerThick(SegNumI) = SD_H1
    ! 
    ! pathways of TNH41/2, NO31/2, CH41/2, SO41/2, TH2S1/2, DIC1/2, TIP1/2, Si1/2
    ! TNH41 and TNH42
    !Dissolved Concentrations 
		SD_NH31 = SD_fdn1*SD_NH3T(1) 
		SD_NH32 = SD_fdn2*SD_NH3T(2)
    SD_JNH4 = SD_s * (SD_NH31 - SD_NH30)
    !
    !Calculate NH4+ and NH3 concentrations
!
    SD1_Ammonia = SD_NH31/(1. + 10.**(-SD_PH1)/10.**(-NH4_NH3_Eqb_Const))
    SD1_Ammonium =  SD_NH31 - SD1_Ammonia
    !SD2_Ammonia = SD_NH3(2)/(1. + 10.**(-SD_pHValue(SegNumI))/10.**(-NH4_NH3_Eqb_Const))
    SD2_Ammonia = SD_NH32/(1. + 10.**(-SD_PH2)/10.**(-NH4_NH3_Eqb_Const))
    SD2_Ammonium =  SD_NH32 - SD2_Ammonia
    !Calculate dissolved <--> gaseous phase distribution
    SedTemp1 = SD_T1 + 273.15
    SedTemp2 = SD_T2 + 273.15
    MW_Constituent = 17 !N = 14, H = 1
    !gp estimated thickness of the aerobic sediment layer 1 (DiToro Appendix B Page 576)  MOVED from below, cb 9/6/13
    !Layer 1
    VolWater = CellArea(LayerNum,SegNumI)*SD_H1
    Call CEMADisGasPhaseDistribution(SD1_Ammonia, MW_Constituent, SedTemp1, HenryConst_NH3, VolWater, AmmoniaG_SD1, AmmoniaD_SD1)
    !Layer 2
    VolWater = CellArea(LayerNum,SegNumI)*SD_H2
    Call CEMADisGasPhaseDistribution(SD2_Ammonia, MW_Constituent, SedTemp2, HenryConst_NH3, VolWater, AmmoniaG_SD2, AmmoniaD_SD2)
    !
    ! NO31 and NO32
    SD_JNO3 = SD_s * (SD_NO31 - SD_NO30)
    !
    !If (SD_SO4 <= SD_Sulfate_CH4_H2S) Then
    If (SD_SO42 <= SD_Sulfate_CH4_H2S) Then         ! 7/26/18
            
    Else
      !Calculate H2S and HS- concentrations
			!
			!SD1_SulfiMinus = SD_HS(1)/(1 + 10**(-SD_pHValue(SegNumI))/10**(-HS_H2S_Eqb_Const))
      SD1_SulfiMinus = SD_HS(1)/(1 + 10**(-SD_PH1)/10**(-HS_H2S_Eqb_Const))
			SD1_Sulfide =  SD_HS(1) - SD1_SulfiMinus
			!SD2_SulfiMinus = SD_HS(2)/(1 + 10**(-SD_pHValue(SegNumI))/10**(-HS_H2S_Eqb_Const))
      SD2_SulfiMinus = SD_HS(2)/(1 + 10**(-SD_PH2)/10**(-HS_H2S_Eqb_Const))
			SD2_Sulfide =  SD_HS(2) - SD2_SulfiMinus
    End If
		!
		!gp   methane or sulfide fluxes produced from C diagenesis 
		!If (SD_SO4 < SD_Sulfate_CH4_H2S) Then
    If (SD_SO42 <= SD_Sulfate_CH4_H2S) Then                 ! cb 7/26/18
			!gp freshwater sediment fluxes - methane
			If (SD_CH4CompMethod == 0) Then
				!***********************************************************************
				!Aqueous methane flux to water column
				SD_SJCH4 = SD_CSODmax - SD_CSOD
				SD_JCH4  = SD_SJCH4
				!Gaseous methane flux to water column
				SD_JCH4g = SD_JC_O2equiv - SD_JCH4 - SD_CSOD           ! flux in gO2/m2/day
				!***********************************************************************
				SD_CH41 = SD_JCH4/SD_s + SD_CH40
			Else
				!***********************************************************************
				!numerical solution for methane
				SD_JCH4 = SD_s * (SD_CH41 - SD_CH40)
				SD_JCH4g = SD_JC_O2equiv - SD_JCH4 - SD_CSOD           ! flux in gO2/m2/day
				!***********************************************************************
			End If
			SD_JHS = 0     !gp
		Else !gp
			!gp marine sediment fluxes - sulfide
			SD_JCH4  = 0
			SD_JCH4g = 0
			SD_JHS   = SD_s * (SD_HS(1) - SD_H2S0)
    End If
    !
    ! SO4
    !Write linear system of equations around sulfate   cb 7/26/18
    !Layer 1
		SD_a11 = -SD_KL12 -SD_s            
		SD_a12 = SD_KL12
    ! SO4ProducedSrc1L1 converted from  g SO4/m?d to g S/m?d ; 32/(4*16+32)=0.33333
    !H2S + 2 O2 --> 2 H+ + SO42-            
		SD_b1 = -SD_s*SD_SO40 - SO4ProducedSrc1L1 * 0.3333333
		!
		!Layer 2
		SD_a21 = SD_KL12
		SD_a22 = -SD_KL12 -  SD_H2 / SD_tc
		SD_b2  = - SD_H2 / SD_tc * SD_so42 - SO4ConsumedSnk1L2 * 0.33333333
		!					
		Call Lin_Sys(SD_a11, SD_a12, SD_a21, SD_a22, SD_b1, SD_b2, SD_so41, SD_so42)  !, NFLog, NFCle)
    SD_SO41 = max(SD_so41,0.0)
    SD_SO42 = max(SD_so42,0.0)           
    ! calculating diffusive flux between layer 1 and water column
    SD_JSO4= SD_s * (SD_so41-SD_so40)
    !
    !TIC
    !Write linear system of equations around Total inorganic carbon
    !Layer 1
		SD_a11 = -SD_KL12 -SD_s            
		SD_a12 = SD_KL12
    ! COCO2ProducedSrc1L1 converted from  g CO2/m?d to g C/m?d ; 12/(2*16+12)=0.272
		SD_b1 = -SD_s*SD_TIC0 - CO2ProducedSrc1L1 * 0.272
		!
		!Layer 2
		SD_a21 = SD_KL12
		SD_a22 = -SD_KL12 -  SD_H2 / SD_tc
		SD_b2  = - SD_H2 / SD_tc * SD_tic2 - CO2ProducedSrc2L2 * 0.272
		!					
		Call Lin_Sys(SD_a11, SD_a12, SD_a21, SD_a22, SD_b1, SD_b2, SD_tic1, SD_tic2)  !, NFLog, NFCle)
    SD_tic1 = max(SD_tic1,0.0)
    SD_tic2 = max(SD_tic2,0.0)
    ! calculating diffusive flux between layer 1 and water column
    SD_JTIC= SD_s * (SD_tic1-SD_TIC0)
    !
    !ALK
    if(IncludeAlkalinity)then
    !Write linear system of equations around Total Alkalinity
    ! Nitrification of ammonium results in an alkalinity decrease: 2 eq. alk per 1 mole ammonium
    ! Denitrification of nitrate (to nitrogen gas) results in an alkalinity increase: 1 eq. alk per 1 mole nitrate
    !Layer 1
		SD_a11 = -SD_KL12 -SD_s
    !SD_a11 = -SD_H1/SD_tc -SD_KL12 -SD_s
		SD_a12 = SD_KL12
		SD_b1 = -SD_s*SD_ALK0 +2.0*SD_NH3toNO3*SD_NH31 - SD_NO31*SD_Denit(1)                   ![m/d]*[mg/m3] 
		!
		!Layer 2
		SD_a21 = SD_KL12   
		SD_a22 = -SD_KL12 -  SD_H2 / SD_tc
		SD_b2  = - SD_H2 / SD_tc * SD_alk2 - SD_NO32*SD_Denit(2)
		!					
		Call Lin_Sys(SD_a11, SD_a12, SD_a21, SD_a22, SD_b1, SD_b2, SD_alk1, SD_alk2)  !, NFLog, NFCle)
    SD_alk1 = max(SD_alk1,0.0)
    SD_alk2 = max(SD_alk2,0.0)
    ! calculating diffusive flux between layer 1 and water column
    SD_JALK= SD_s * (SD_alk1-SD_alk0)
    end if
    !
    !PO4
    if (SD_O20 >= SD_DOcr) then
        SD_Kdp1 = SD_Kdp2 * SD_delta_kpo41
    else
        SD_Kdp1 = SD_Kdp2 * SD_delta_kpo41 ** (SD_O20 / SD_DOcr)
    end if  
    !calculating dissolved and particulate forms of phosphorus (Chapra, eqn. 25.89)
    SD_fd1 = 1.0/(1.0 + SD_Kdp1*SD_Rho*(1.0-SD_Porosity))
    SD_fp1 = 1.0 - SD_fd1
    SD_fd2 = 1.0/(1.0 + SD_Kdp2*SD_Rho*(1.0-SD_Porosity))
    SD_fp2 = 1.0 - SD_fd2
    !    
    !Write linear system of equations around total phosphate SD_PO4T
    !Layer 1
		SD_a11 = -SD_fd1*SD_KL12 - SD_fp1*SD_w12 - SD_fd1*SD_s - SD_w2  
    SD_a12 = SD_fd2*SD_KL12 + SD_fp2*SD_w12 
    SD_b1 = -SD_s*fdp*SD_PO40                   ![m/d]*[mg/m3] 
    !
    !Layer 2
    SD_a21 = SD_fd1*SD_KL12 + SD_fp1*SD_w12 + SD_w2
    SD_a22 = -SD_fd2 * SD_KL12 - SD_fp2 * SD_w12 - SD_w2 - SD_H2 / SD_tc
    SD_b2  = -SD_Jp - SD_H2 / SD_tc * SD_PO4T2 - vsss*SD_PO40
!					
    Call Lin_Sys(SD_a11, SD_a12, SD_a21, SD_a22, SD_b1, SD_b2, SD_PO4T1, SD_PO4T2)   !, NFLog, NFCle) 
    SD_PO4T1 = max(SD_PO4T1,0.0)
    SD_PO4T2 = max(SD_PO4T2,0.0)
    ! dissolved PO4
    SD_PO4(1) = SD_fd1*SD_PO4T1
    SD_PO4(2) = SD_fd2*SD_PO4T2		
    SD_JPO4   = SD_s * (SD_PO4(1)-fdp*SD_PO40)
    !
    ! Temperature
    !Write linear system of equations for temperature
    !SD_rhowcp = 4.186e6   !4.186*1.0e6   ! units J g-1 C-1 * g m-3= J C-1 m-3
		!Equation 5.1
		!
    !Layer 1
    !SD_a11 = -SD_KL12 -SD_s
    SD_a11 = -SD_KL12 -SD_s  -SD_H1/SD_tc   ! SW 8/29/2021
    SD_a12 = SD_KL12          
    !SD_b1 = -SD_s*SD_Tw 
    SD_b1 = -SD_s*SD_Tw -SD_H1/SD_tc * SD_T1  ! SW 8/29/2021

    !
    !Layer 2
    SD_a21 = SD_KL12 * SD_rhowcp
    !SD_a22 = -SD_KL12 * SD_rhowcp -  SD_H2 / SD_tc * SD_rhowcp   + SD_Ksw        ! SD_tc is time step in days=CUF*DLT/DAY
    SD_a22 = -SD_KL12 * SD_rhowcp -  SD_H2 / SD_tc * SD_rhowcp   - SD_Ksw*DAY   ! SW 8/28/2021  ! SD_tc is time step in days=CUF*DLT/DAY

    !SD_b2  = - SD_H2 * SD_rhowcp / SD_tc * SD_T2 - SD_Ksw * SD_Tsed
    SD_b2  = - SD_H2 * SD_rhowcp / SD_tc * SD_T2 - SD_Ksw * SD_Tsed*DAY
		!					
		Call Lin_Sys(SD_a11, SD_a12, SD_a21, SD_a22, SD_b1, SD_b2, SD_T1, SD_T2)   !, NFLog, NFCle)
    SD_T1 = max(SD_T1,0.0)
    SD_T2 = max(SD_T2,0.0)
    ! calculating diffusive heat flux between layer 1 and water column
    SD_JT= SD_s * SD_rhowcp* (SD_T1-SD_Tw)
    !
    if(IncludeDynamicpH)then
    ! Sediment pH for layers 1 and 2
      SD_TDS=0.0  !sediments not simulating tds at the moment        
      SD_POCT1=0.0 ! POC not predicted for aerobic layer.
      call PH_SEDIMENTS(SD_T1,sd_tic1,sd_alk1,sd_nh31,sd_po4(1),SD_POCT1,SD_TDS(1),SD_PH1)  !layer 1
      call PH_SEDIMENTS(SD_T2,sd_tic2,sd_alk2,sd_nh32,sd_po4(2),SD_POCT2,SD_TDS(2),SD_PH2)  !layer 2
    end if
    !
    !Bubbles formation
		!Aerobic layer is thin and so ignore gas formation from aerobic layer
    !H2S
		!Calculate dissolved <--> gaseous phase distribution
    IF(Bubbles_Calculation) THEN
		  SedTemp1 = SD_T1 + 273.15
		  MW_Constituent = 36. !S = 34, H = 1
		  !Calculate dissolved <--> gaseous phase distribution
		  !VolWater = CellArea(SegNumI)*SD_H2
		  SedTemp2 = SD_T2 + 273.15
		  !Layer 1
      VolWater = CellArea(LayerNum,SegNumI)*SD_H1
		  Call CEMADisGasPhaseDistribution(SD1_Sulfide, MW_Constituent, SedTemp1, HenryConst_H2S, VolWater, SulfideG_SD1, SulfideD_SD1)
		  !Layer 2
		  VolWater = CellArea(LayerNum,SegNumI)*SD_H2
		  Call CEMADisGasPhaseDistribution(SD2_Sulfide, MW_Constituent, SedTemp2, HenryConst_H2S, VolWater, SulfideG_SD2, SulfideD_SD2)
		  !Old TConc(1,SegNumI) = SD2_Sulfide
		  !Old SConc(1,SegNumI) = (TConc(1,SegNumI) - TConcP(1,SegNumI))/dlt    !gm/m�/s
		  TConc(1,LayerNum,SegNumI) = SulfideG_SD2*BubbAccFraction + TConcP(1,LayerNum,SegNumI)
		  SConc(1,LayerNum,SegNumI) = (SulfideG_SD2*BubbAccFraction)/sd_tc    !gm/m�/s
		  !Old TConcP(1,SegNumI) = SD2_Sulfide
		  TConcP(1,LayerNum,SegNumI) = TConc(1,LayerNum,SegNumI)
		
		  !CH4
		  TConc(2,LayerNum,SegNumI) = SD_CH41*BubbAccFraction + TConcP(2,LayerNum,SegNumI)
		  SConc(2,LayerNum,SegNumI) = (SD_CH41*BubbAccFraction)/sd_tc    !gm/m�/s
		  !Old TConcP(1,SegNumI) = SD2_Sulfide
		  TConcP(2,LayerNum,SegNumI) = TConc(2,LayerNum,SegNumI)
		
		  !NH3
		  !Old TConc(3,SegNumI) = SD2_Ammonia
		  !Old SConc(3,SegNumI) = (TConc(3,SegNumI) - TConcP(3,SegNumI))/dlt    !gm/m�/s
		  TConc(3,LayerNum,SegNumI) = SD_NH32 !SD2_Ammonia*BubbAccFraction + TConcP(3,SegNumI)
		  SConc(3,LayerNum,SegNumI) = (SD_NH32-TConcP(3,LayerNum,SegNumI))/sd_tc    !gm/m�/s
		  !Old TConcP(3,SegNumI) = SD2_Ammonia
		  TConcP(3,LayerNum,SegNumI) = TConc(3,LayerNum,SegNumI)
		
		  !CO2
		  !Old TConc(4,SegNumI) = CO2ProducedCon2L2
		  !Old SConc(4,SegNumI) = (TConc(4,SegNumI) - TConcP(4,SegNumI))/dlt    !gm/m�/s
		  TConc(4,LayerNum,SegNumI) = CO2ProducedCon2L2*BubbAccFraction + TConcP(4,LayerNum,SegNumI)
		  SConc(4,LayerNum,SegNumI) = (CO2ProducedCon2L2*BubbAccFraction)/sd_tc    !gm/m�/s
		  !Old TConcP(4,SegNumI) = CO2ProducedCon2L2
		  TConcP(4,LayerNum,SegNumI) = TConc(4,LayerNum,SegNumI)
		
		  !H2S
		  DissolvedGasSediments(1, LayerNum, SegNumI) = SulfideD_SD2
      !CH4
		  DissolvedGasSediments(2, LayerNum, SegNumI) = SD_CH41
		  !NH4
		  DissolvedGasSediments(3, LayerNum, SegNumI) = AmmoniaD_SD2
		  !CO2
		  DissolvedGasSediments(4, LayerNum, SegNumI) = 0.d0
    !
    END IF
		!
		IF(BUBBLES_CALCULATION) Call GasBubblesFormation(BubbleRadiusSed(SegNumI), sd_tc, VolWater)
    !
    !2
		!Porewater release
		!Flux of CH4d, NH3d + NH4d, H2Sd + HSd, SO42-d, NO3d, CO2d 
		!Volume of porewater
		VolWater = CellArea(LayerNum,SegNumI)*SD_H2
    !Aerobic Layer
    PW_RelRate1 = PorewaterRelRate(SegNumI)*SD_H1/(SD_H1 + SD_H2)*DAY  !m�/d
    !Anaerobic Layer
    PW_RelRate2 = PorewaterRelRate(SegNumI)*SD_H2/(SD_H1 + SD_H2)*DAY  !m�/d
		!Dissolved_CH4_Src = SD_JCH4*SD_tc/(SD_H1 + SD_H2)   !CHECK CODE This is incorrect original code since SD_tc is in units of s*d/s or d ==> g/m3 Why didvide by SD_H1 and SD_H2 inncorrect
    !Dissolved_CH4_Src = SD_JCH4*SD_tc/(SD_H1 + SD_H2)   !SW 10/10/2017    gO2/m2/d ==> gO2/m3/d
		!Dissolved_NH3_Src = (AmmoniaD_SD1*SD_H1 + AmmoniaD_SD2*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !gm/m�/s
		!Dissolved_NH3_Src = Dissolved_NH3_Src + (SD1_Ammonium*SD_H1 + SD2_Ammonium*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))  !gm/m�/s
    Dissolved_NH3_Src = (AmmoniaD_SD1*PW_RelRate1 + AmmoniaD_SD2*PW_RelRate2)/(CellThickness*CellArea(LayerNum,SegNumI))                      !gm/m�/d
		Dissolved_NH3_Src = Dissolved_NH3_Src + (SD1_Ammonium*PW_RelRate1 + SD2_Ammonium*PW_RelRate2)/(CellThickness*CellArea(LayerNum,SegNumI))  !gm/m�/d
    !Dissolved_H2S_Src = (SulfideD_SD1*SD_H1 + SulfideD_SD2*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !gm/m�/s
		!Dissolved_H2S_Src = Dissolved_H2S_Src + (SD1_SulfiMinus*SD_H1 + SD2_SulfiMinus*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))  !gm/m�/s
    Dissolved_H2S_Src = (SulfideD_SD1*PW_RelRate1 + SulfideD_SD2*PW_RelRate2)/(CellThickness*CellArea(LayerNum,SegNumI))                      !gm/m�/d
		Dissolved_H2S_Src = Dissolved_H2S_Src + (SD1_SulfiMinus*PW_RelRate1 + SD2_SulfiMinus*PW_RelRate2)/(CellThickness*CellArea(LayerNum,SegNumI))  !gm/m�/d
		!Dissolved_NO3_Src = (SD_NO3(1)*SD_H1 + SD_NO3(2)*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !gm/m�/s
    Dissolved_NO3_Src = (SD_NO31*PW_RelRate1 + SD_NO32*PW_RelRate2)/(CellThickness*CellArea(LayerNum,SegNumI))                      !gm/m�/d
		!Dissolved_SO4_Src = SD_SO4 * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !gm/m�/s
    !Dissolved_SO4_Src = SD_SO4 * (PW_RelRate1+PW_RelRate2)/(CellThickness*CellArea(SegNumI))                      !gm/m�/s
    Dissolved_SO4_Src = (SD_SO41*PW_RelRate1 + SD_SO42*PW_RelRate2)/(CellThickness*CellArea(LayerNum,SegNumI))     ! cb 7/26/18		
		!Dissolved_CO2_Src = (CO2ProducedCon1L1*SD_H1 + CO2ProducedCon2L2*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !gm/m�/s
    Dissolved_CO2_Src = (CO2ProducedCon1L1*PW_RelRate1 + CO2ProducedCon2L2*PW_RelRate2)/(CellThickness*CellArea(LayerNum,SegNumI))                      !gm/m�/d
    !Dissolved_ALK_Src = (SD_ALK(1)*SD_H1 + SD_ALK(2)*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !gm/m�/s
    if(IncludeAlkalinity)then
        Dissolved_ALK_Src = (SD_ALK1*PW_RelRate1 + SD_ALK2*PW_RelRate2)/(CellThickness*CellArea(LayerNum,SegNumI))                      !gm/m�/d
    else
        Dissolved_ALK_Src = 0.0
    end if
    !Dissolved_PO4_Src = (SD_PO4(1)*SD_H1 + SD_PO4(2)*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !gm/m�/s
    Dissolved_PO4_Src = (SD_PO4(1)*PW_RelRate1 + SD_PO4(2)*PW_RelRate2)/(CellThickness*CellArea(LayerNum,SegNumI))                      !gm/m�/d
    !Dissolved_Fe2_Src = (SD_Fe2(1)*SD_H1 + SD_Fe2(2)*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !gm/m�/s
    if(IncludeIron)then
        Dissolved_Fe2_Src = (SD_Fe2(1)*PW_RelRate1 + SD_Fe2(2)*PW_RelRate2)/(CellThickness*CellArea(LayerNum,SegNumI))                      !gm/m�/d
    else
        Dissolved_Fe2_Src = 0.0
    end if
    !Dissolved_Mn2_Src = (SD_Mn2(1)*SD_H1 + SD_Mn2(2)*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !gm/m�/s
    if(IncludeManganese)then
        Dissolved_Mn2_Src = (SD_Mn2(1)*PW_RelRate1 + SD_Mn2(2)*PW_RelRate2)/(CellThickness*CellArea(LayerNum,SegNumI))                      !gm/m�/d
    else
        Dissolved_Mn2_Src = 0.0
    end if
    !Sediment_Heat_Src = SD_rhowcp*(SD_T(1)*SD_H1 + SD_T(2)*SD_H2)/(SD_H1 + SD_H2) * PorewaterRelRate(SegNumI)/(CellThickness*CellArea(SegNumI))                      !J/m�/s
    Sediment_Heat_Src = SD_rhowcp*(SD_T1*PW_RelRate1 + SD_T2*PW_RelRate2)/(CellThickness*CellArea(LayerNum,SegNumI))                      !J/m�/d
    !3
    !Diffusive flux of NH3, NO3, CH4, SO4, H2S, CO2
    Dissolved_NH3_Src = Dissolved_NH3_Src + SD_JNH4*BIBH22   !g/m�/d SD_JNH4/CellThickness = g/m�/d/m = g/m�/d
    Dissolved_NO3_Src = Dissolved_NO3_Src + SD_JNO3*BIBH22   !g/m�/d SD_JNO3/CellThickness = g/m�/d/m = g/m�/d 
    Dissolved_CH4_Src = Dissolved_CH4_Src + SD_JCH4*BIBH22   !g/m�/d SD_JCH4/CellThickness = g/m�/d/m = g/m�/d commented out because double counting... ! SW 10/10/2017 added back because eliminated the problem above
    Dissolved_SO4_Src = Dissolved_SO4_Src + SD_JSO4*BIBH22   !g/m�/d SD_JSO4/CellThickness = g/m�/d/m = g/m�/d	
    Dissolved_H2S_Src = Dissolved_H2S_Src + SD_JHS*BIBH22
    Dissolved_CO2_Src = Dissolved_CO2_Src + SD_JTIC*BIBH22
    IF(IncludeAlkalinity) Dissolved_Alk_Src = Dissolved_ALK_Src + SD_JALK*BIBH22
    Dissolved_PO4_Src   = Dissolved_PO4_Src + SD_JPO4*BIBH22
    IF(IncludeIron)       Dissolved_Fe2_Src = Dissolved_Fe2_Src + SD_JFe2*BIBH22
    IF(IncludeManganese)  Dissolved_Mn2_Src = Dissolved_Mn2_Src + SD_JMn2*BIBH22
    !
    !GasReleaseCH4 = GasReleaseCH4+SD_JCH4g*B(LayerNum,SegNumI)*DLX(SegNumI)*SD_tc/2.67  ! Flux of gas CH4 in gC   SW 10/10/2017
    GasReleaseCH4 = GasReleaseCH4+SD_JCH4g*B(LayerNum,SegNumI)*DLX(SegNumI)*SD_tc/5.33  ! Flux of gas CH4 in gC   SW 5/26/2022  5.33 g O2/g CH4-C

    !
    Sediment_Heat_Src = Sediment_Heat_Src + SD_JT*BIBH22  ! Heat J/m3/d  SD_JT/CellThickness = J/m3/d/m= J/m3/d
		!Flux to CSOD and NSOD
		Dissolved_O2_Snk = SD_JSOD*BIBH22   !g/m�/d SD_JSOD/CellThickness = g/m�/d/m = g/m�/d
    !	
		!Convert all source/sink to g/m�/s from g/m�/d
		!NH3, NO3, CH4, SO4, DO, CO2, H2S
		Dissolved_NH3_Src   = Dissolved_NH3_Src/DAY      !g/m�/d --> g/m�/s
    Dissolved_NO3_Src   = Dissolved_NO3_Src/DAY      !g/m�/d --> g/m�/s
    !Dissolved_CH4_Src   = Dissolved_CH4_Src/DAY      !g/m�/d --> g/m�/s
    Dissolved_CH4_Src   = Dissolved_CH4_Src/DAY/5.33      !gO2/m�/d --> gC/m�/s   SW 5/26/2022

    Dissolved_SO4_Src   = Dissolved_SO4_Src/DAY      !g/m�/d --> g/m�/s
		Dissolved_CO2_Src   = Dissolved_CO2_Src/DAY      !g/m�/d --> g/m�/s
		Dissolved_H2S_Src   = Dissolved_H2S_Src/DAY/1.88      !gO2/m�/d --> gS/m�/s   SW 5/26/2022
		Dissolved_O2_Snk    = Dissolved_O2_Snk/DAY       !g/m�/d --> g/m�/s
    Dissolved_ALK_Src   = Dissolved_ALK_Src/DAY      !g/m�/d --> g/m�/s
    Dissolved_PO4_Src   = Dissolved_PO4_Src/DAY      !g/m�/d --> g/m�/s
    Dissolved_Fe2_Src   = Dissolved_Fe2_Src/DAY      !g/m�/d --> g/m�/s
    Dissolved_Mn2_Src   = Dissolved_Mn2_Src/DAY      !g/m�/d --> g/m�/s
    Sediment_Heat_Src   = Sediment_Heat_Src/DAY      !J/m�/d --> J/m�/s
    ! resuspension of POM, POP, and PON
    IF(CEMA_POM_Resuspension) THEN
      DO iTemp = 1, 2  ! only labile and refractory, not including inert for now 
        LPOM_Resuspension  = LPOM_Resuspension  + SD_EPOC(itemp)*BIBH22/ORGC(JW)
        RPOM_Resuspension  = RPOM_Resuspension  + SD_EPOC(itemp)*BIBH22/ORGC(JW)
        LPOMN_Resuspension = LPOMN_Resuspension + SD_EPON(itemp)*BIBH22
        RPOMN_Resuspension = RPOMN_Resuspension + SD_EPON(itemp)*BIBH22
        LPOMP_Resuspension = LPOMP_Resuspension + SD_EPOP(itemp)*BIBH22
        RPOMP_Resuspension = RPOMP_Resuspension + SD_EPOP(itemp)*BIBH22
      END DO
    END IF        
  End Subroutine
    
  Subroutine PH_SEDIMENTS(t1sed,ticsed,alksed,nh4sed,po4sed,pocsed,tdssed,phsed) ! Enhancements added for buffering by ammonia, phosphate, and OM ! SR 01/01/12
    ! pH and carbonate species
    !
    implicit none
    real(R8):: t1sed,ticsed,alksed,nh4sed,po4sed,pocsed,tdssed,phsed
    
     T1K = t1sed + 273.15
     CART = ticsed/12011. ! SR 01/01/12
     ALKT = alksed/50044. ! SR 01/01/12
     AMMT = NH4sed/14006.74 ! SR 01/01/12
     PHOST = PO4sed/30973.762 ! SR 01/01/12
     !OMCT = (LDOM(K,I)+RDOM(K,I))*ORGC(JW)/12011. ! moles carbon per liter from DOM ! SR 01/01/12
     omct=0.0  ! DOM is not simulated in the sediments yet...
     !IF (POM_BUFFERING) OMCT = OMCT + (LPOM(K,I)+RPOM(K,I))*ORGC(JW)/12011. ! SR 01/01/12
     IF (POM_BUFFERING) OMCT = OMCT + pocsed/12011. ! SR 01/01/12
     omct=0.0    ! 
 !**** Ionic strength
     IF (FRESH_WATER(JW)) S2 = 2.5E-05*TDSsed
     IF (SALT_WATER(JW))  S2 = 1.47E-3+1.9885E-2*TDSsed+3.8E-5*TDSsed*TDSsed
!**** Debye-Huckel terms and activity coefficients
     SQRS2 = SQRT(S2)
     DH1 = -0.5085*SQRS2/(1.0+1.3124*SQRS2)+4.745694E-03+4.160762E-02*S2-9.284843E-03*S2*S2
     DH2 = -2.0340*SQRS2/(1.0+1.4765*SQRS2)+1.205665E-02+9.715745E-02*S2-2.067746E-02*S2*S2
     DH3 = -4.5765*SQRS2/(1.0+1.3124*SQRS2) ! extended Debye-Huckel for PO4 ! SR 01/01/12
     DHH = -0.5085*SQRS2/(1.0+2.9529*SQRS2) ! extended Debye-Huckel for H+ ion ! SR 01/01/12
     H2CO3T= 10.0**(0.0755*S2)
     HCO3T = 10.0**DH1
     CO3T  = 10.0**DH2
     PO4T  = 10.0**DH3 ! SR 01/01/12
     HT    = 10.0**DHH ! activity coefficient for H+ ! SR 01/01/12
     HPO4T = CO3T   ! tabled values similar to those for carbonate ! SR 01/01/12
     OHT   = HCO3T  ! tabled values similar to those for bicarbonate ! SR 01/01/12
     H2PO4T= HCO3T  ! tabled values similar to those for bicarbonate ! SR 01/01/12
     NH4T  = HCO3T  ! tabled values similar to those for bicarbonate ! SR 01/01/12
     NH3T  = H2CO3T ! neutral species, set coefficient to same as that for carbonic acid ! SR 01/01/12
     H3PO4T= H2CO3T ! neutral species, set coefficient to same as that for carbonic acid ! SR 01/01/12
!**** Temperature adjustment
     KW   = 10.0**(-283.971 -0.05069842*T1K +13323.0/T1K +102.24447*LOG10(T1K) -1119669.0/(T1K*T1K))/OHT
     K1   = 10.0**(-356.3094 -0.06091964*T1K +21834.37/T1K +126.8339 *LOG10(T1K) -1684915 /(T1K*T1K))*H2CO3T/HCO3T
     K2   = 10.0**(-107.8871 -0.03252849*T1K + 5151.79/T1K + 38.92561*LOG10(T1K) - 563713.9/(T1K*T1K))*HCO3T/CO3T
     KAMM = 10.0**(-0.09018 -2729.92/T1K)*NH4T/NH3T ! SR 01/01/12
     KP1  = 10.0**(4.5535 -0.013486*T1K -799.31/T1K)*H3PO4T/H2PO4T ! Bates (1951) ! SR 01/21/12
     KP2  = 10.0**(5.3541 -0.019840*T1K -1979.5/T1K)*H2PO4T/HPO4T ! Bates and Acree (1943) ! SR 01/21/12
     KP3  = 10.0**(-12.38) *HPO4T/PO4T ! Dean (1985 )! SR 01/01/12
!**** pH evaluation
     !PHT = -PH(K,I)-2.1
     !IF (PH(K,I) <= 0.0) PHT = -14.0
     PHT = -phsed-2.1
     IF (phsed <= 0.0) PHT = -14.0
     INCR1 = 10.0
     DO N=1,3
        F1 = 1.0
        INCR1 = INCR1/10.0
        ITER1 = 0
        DO WHILE (F1 > 0.0 .AND. ITER1 < 12)
          PHT = PHT+INCR1
          HION = 10.0**PHT
          F1 = CART*K1*(HION+2.0*K2)/(HION*HION+K1*HION+K1*K2)+KW/HION-ALKT-HION/HT ! SR 01/01/12
          IF (AMMONIA_BUFFERING) THEN ! SR 01/01/12
            F1 = F1 + AMMT*KAMM/(HION+KAMM) ! SR 01/01/12
          END IF ! SR 01/01/12
          IF (PHOSPHATE_BUFFERING) THEN ! SR 01/01/12
            F1 = F1+ PHOST*( KP1*KP2*HION + 2*KP1*KP2*KP3 - HION*HION*HION ) &
                /( HION*HION*HION + KP1*HION*HION + KP1*KP2*HION + KP1*KP2*KP3) ! SR 01/01/12
          END IF ! SR 01/01/12
          IF (OM_BUFFERING) THEN ! SR 01/01/12
            DO JA=1,NAG ! SR 01/01/12
              F1 = F1 + OMCT*SDEN(JA)*( 1.0/(1.0+HION*(10.0**PK(JA))) - 1.0/(1.0+(10.0**(PK(JA)-4.5))) ) ! SR 01/01/12
            END DO ! SR 01/01/12
          END IF ! SR 01/01/12
          ITER1 = ITER1+1
        END DO
        PHT = PHT-INCR1
     END DO
!**** pH, carbon dioxide, bicarbonate, and carbonate concentrations
     HION   = 10.0**PHT
     PHsed  = -PHT
     CO2sed = TICsed/(1.0 + K1/HION+K1*K2/(HION*HION))
     HCO3sed= TICsed/(1.0 + HION/K1+K2/HION)
     CO3sed = TICsed/((HION*HION)/(K1*K2)+HION/K2 + 1.0)
  End Subroutine

  Subroutine CEMAWindInducedSedimentResuspension
    epsilon=0.0
    FETCHW = FETCHD(SegNumI,JB)
    IF (COS(PHI(JW)-PHI0(SegNumI)) < 0.0) FETCHW = FETCHU(SegNumI,JB)
    FETCHW = MAX(FETCHW,BI(KT,SegNumI),DLX(SegNumI))
    U2    = WIND(JW)*WSC(SegNumI)*WIND(JW)*WSC(SegNumI)+NONZERO
    COEF1 = 0.53  *(G*DEPTHB(LayerNum,SegNumI)/U2)**0.75
    COEF2 = 0.0125*(G*FETCHW/U2)**0.42
    COEF3 = 0.833* (G*DEPTHB(LayerNum,SegNumI)/U2)**0.375
    COEF4 = 0.077* (G*FETCHW/U2)**0.25
    HS    = 0.283 *U2/G*0.283*TANH(COEF1)*TANH(COEF2/TANH(COEF1))
    !TS    = 2.0*PI*U2/G*1.2*  TANH(COEF3)*TANH(COEF4/TANH(COEF3))
    TS    = 2.0*PI*sqrt(U2)/G*1.2*  TANH(COEF3)*TANH(COEF4/TANH(COEF3))   ! cb 5/9/14
    LW0    = G*TS*TS/(2.0*PI)  
 
    LW1 = LW0
    LW  = LW0*TANH(2.0*PI*DEPTHB(LayerNum,SegNumI)/LW1)
    DO WHILE (ABS(LW-LW1) > 0.001)
      LW1 = LW
      LW  = LW0*TANH(2.0*PI*DEPTHB(LayerNum,SegNumI)/LW1)
    END DO
    COEF = MIN(710.0,2.0*PI*DEPTHB(LayerNum,SegNumI)/LW)
    UORB = PI*HS/TS*100.0/SINH(COEF)
    TAU  = 0.003*UORB*UORB
    IF (TAU-TAUCRPOM > 0.0) EPSILON = MAX(0.0,0.008/49.0*(TAU-TAUCRPOM)**3*10000.0/sd_tc)						        
    SD_E = EPSILON*DLX(SegNumI)*BI(LayerNum,SegNumI)/VOL(LayerNum,SegNumI)  ! SD_E: g/m^2/s
  End Subroutine
        
  Subroutine CEMABottomScourResuspension
    if(cao_method)then
      reyn_resusp = dia_POM * sqrt(spgrav_POM * g * dia_POM)
      if(reyn_resusp < 6.0)then
        crshields = 0.1414 * reyn_resusp ** (-0.2306)
      else if(reyn_resusp >= 6.0 .and. reyn_resusp <= 282.8)then
        crshields = ( 1.0 + (0.0223 * reyn_resusp)**2.8358)**0.3542 / (3.0946 * reyn_resusp**0.6769)
      else if(reyn_resusp > 282.8)then
        crshields = 0.045
      end if
    end if  
    molvisc_h2o = 1.79e-6 * exp(-0.0266 * SD_T1)     ! LB 3/2019 SW 3/2019
    shields = SD_taubot / (g*(spgrav_POM-1.0) * dia_POM)
  
    Vscour = 0.00033*(shields/crshields - 1.0)*(spgrav_POM-1.0)**0.6 * g**0.6 * dia_POM**0.8/molvisc_h2o
    c_bottom = (c2(LayerNum,SegNumI,NLPOM) + c2(LayerNum,SegNumI,NRPOM)) * dexp( poms(jw) * h(LayerNum,jw) / DZ(LayerNum -1, SegNumI))
    
    if(spgrav_POM < 1.2)then
      c_bottom2=1.0
    else if(spgrav_POM >= 1.2 .and. spgrav_POM < 1.8)then
       c_bottom2= 1.0 * (1.8 - spgrav_POM)/0.6 + 3.0 * (spgrav_POM - 1.2)/0.6
    else if(spgrav_POM >= 1.8 .and. spgrav_POM <= 2.2)then
        c_bottom2= 3.0 * (2.2 - spgrav_POM)/0.4 + 5.0 * (spgrav_POM - 1.8)/0.4
    else if(spgrav_POM > 2.2)then
      c_bottom2=5.0
    end if 
  
    c_bottom=dmin1(c_bottom,c_bottom2)

    if(Vscour > 0.0)then
    SD_E = c_bottom*Vscour
    else
      SD_E = 0.0
    end if
  End Subroutine   

  Subroutine Lin_Sys(a11, a12, a21, a22, b1, b2, x1, x2)   !,NFLog,NFCle) 
	USE SCREENC, ONLY:JDAY  
  Real(8) a11, a12, a21, a22, b1, b2, x1, x2
	  !Byte NFLog, NFCle
	  !from 03-Nov-2003 version of Q2KMaster 
	  !This subroutine solves a linear system of 2 equations and 2 unknowns 
	  If (a11 * a22 - a12 * a21 == 0) Then 
		  !MsgBox "The sediment flux solution matrix is singular: " & a11 & ", " & a12 & ", " & a21 & ", " & a22 
		  !Write(NFLog,'(a)') 'The sediment flux solution matrix is singular: '
		  !Write(NFLog,*) 'a11  == ', a11, 'a12 = ', a12, 'a21 = ', a21, 'a22 = ',a22 
		  Write(w2err,'(a,F10.3)') 'The sediment flux solution matrix is singular on JDAY:', JDAY
		  Write(w2err,*) 'a11  = ', a11, 'a12 = ', a12, 'a21 = ', a21, 'a22 = ',a22 
          WRITE(W2ERR,'(A,E15.8,A,E15.8)') 'b1 = ', b1, 'b2  = ', b2
		  Write(w2err,*) 'Error in the solution of linear system of equations used in sediment diagenesis model'
		  !Write(NFLog,*) 'Error in the solution of linear system of equations used in sediment diagenesis model'
      write(w2err,'(A)')'Error in Sediment Diagenesis, review output files and review the sediment diagensis parameters'
		  Stop 'Please review the sediment diagensis parameters' 
	  End If 
	  x1 = (a22 * b1 - a12 * b2) / (a11 * a22 - a12 * a21) 
	  x2 = (a11 * b2 - a21 * b1) / (a11 * a22 - a12 * a21) 
  End Subroutine 


  Subroutine CEMADisGasPhaseDistribution(CTotal, MolWt, GasTemp, HenryConst, Vwtr, CGasPh, CLiqPh)
    Use CEMAVars
    
    Real(R8) CTotal, MolWt, GasTemp, HenryConst
    Real(R8) Vwtr, CGasPh, CLiqPh
    
    CGasPh  =   CTotal/(1 + GasConst_R*GasTemp/HenryConst)
    CLiqPh  =   CTotal*(GasConst_R*GasTemp/HenryConst)/(1 + GasConst_R*GasTemp/HenryConst)
  End Subroutine

End Module CEMASedimentDiagenesis

