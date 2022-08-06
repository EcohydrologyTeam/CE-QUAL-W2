!***********************************************************************************************************************************
!**                                                                                                                               **
!**                                                         CE-QUAL-W2-PRE                                                        **
!**                                                          Version 4.5                                                          **
!**                                                                                                                               **
!**                                                     A preprocessor code for                                                   **
!**                                                        CE-QUAL-W2 v4.5                                                        **
!**                                                                                                                               **
!**                                                 Currently maintained by:                                                      **
!**                                                        Scott A. Wells                                                         **
!**                                       Department of Civil and Environmental Engineering                                       **
!**                                                  Portland State University                                                    **
!**                                                         PO Box 751                                                            **
!**                                                 Portland, Oregon USA 97207-0751                                                  **
!**                                                 phone number: (503) 725-4276                                                  **
!**                                                 e-mail: wellss@pdx.edu                                                        **                 **
!**                                                                                                                               **
!**  Main contributors to CE-QUAL-W2 are shown in User Manual. Primary developer at the Corps has been                            **
!**                                                  Thomas M. Cole, Retired                                                      **
!**                                                Water Quality Modeling Group                                                   **
!**                                                U.S. Army Corps of Engineers                                                   **
!**                                                Waterways Experiment Station                                                   **
!**                                                Vicksburg, Mississippi 39180                                                   **
!**                                                                                                                               **
!**                                                                                                                               **
!**                                                                                                                               **
!***********************************************************************************************************************************

! Module definitions
MODULE PREC
  INTEGER, PARAMETER                               :: I2=SELECTED_INT_KIND(3)
  INTEGER, PARAMETER                               :: R8=SELECTED_REAL_KIND(15)
END MODULE PREC
MODULE FILE_ERROR
 CONTAINS
   SUBROUTINE DUPLICATE_FILENAME(ERR,FILE)
     INTEGER                                       :: ERR
     CHARACTER(2)                                  :: L
     CHARACTER(5)                                  :: FMTA
     CHARACTER(72)                                 :: FILE
     WRITE (L,'(I0)') LEN_TRIM(FILE)
     FMTA = '(A)'
     CALL ERRORS
     WRITE (ERR,FMTA) 'Duplicate file name in control file for file '//TRIM(FILE)
   END SUBROUTINE DUPLICATE_FILENAME
END MODULE FILE_ERROR
MODULE CONTROL
  REAL,    TARGET,   ALLOCATABLE, DIMENSION(:,:,:) :: C2,     CD
  INTEGER                                          :: IMX,    KMX,    UNIT
  INTEGER                                          :: NERR=0, NWRN=0, ERR=370
  INTEGER                                          :: WIN=6
  INTEGER                                          :: NASS,   NAAL, NBOD
  INTEGER,           ALLOCATABLE, DIMENSION(:)     :: SSAC,   ALAC
  LOGICAL                                          :: DELETE_ERR, DELETE_WRN,DOWNSTREAMDAM=.TRUE.  ! POLYGONS
END MODULE CONTROL
MODULE GRIDC
  REAL,              ALLOCATABLE, DIMENSION(:)     :: Z, XF
  REAL,              ALLOCATABLE, DIMENSION(:)     :: CMULT,  DMULT, HMULT, CDMULT
  REAL,              ALLOCATABLE, DIMENSION(:)     :: ICETH,  ALPHA
  REAL,              ALLOCATABLE, DIMENSION(:,:)   :: H,      EL
  REAL,              POINTER,     DIMENSION(:,:)   :: T2
  INTEGER                                          :: JW,     NAC,    INI=10
  INTEGER,           ALLOCATABLE, DIMENSION(:)     :: IEPR,   KEPR,   NADC,   NACTR,  NACIN,  NACPR,  NACDT
  INTEGER,           ALLOCATABLE, DIMENSION(:)     :: KTWB,   KB,     CUS,    US,     DS,     BS,     BE
  INTEGER,           ALLOCATABLE, DIMENSION(:)     :: CN
  INTEGER,           ALLOCATABLE, DIMENSION(:,:)   :: IPR,    CDN
  LOGICAL,           ALLOCATABLE, DIMENSION(:)     :: ICE_CALC
  CHARACTER(9),      ALLOCATABLE, DIMENSION(:)     :: HFMT
  CHARACTER(72),     ALLOCATABLE, DIMENSION(:)     :: TITLE
  CHARACTER(8),      ALLOCATABLE, DIMENSION(:,:)   :: HPRWBC, CDWBC,  CPRWBC
  CHARACTER(10),     ALLOCATABLE, DIMENSION(:,:)   :: CONV
  CHARACTER(43)                                    :: NAME
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: LNAME                      ! SW 1/16/04 added section
  !CHARACTER(8),      ALLOCATABLE, DIMENSION(:)       :: CDNAME2
  CHARACTER(6),      ALLOCATABLE, DIMENSION(:)       :: CUNIT,  CUNIT2
  CHARACTER(9),      ALLOCATABLE, DIMENSION(:)       :: FMTH,   FMTC,   FMTCD
  CHARACTER(19),     ALLOCATABLE, DIMENSION(:)       :: CNAME1
  CHARACTER(43),     ALLOCATABLE, DIMENSION(:)       :: CNAME,  CDNAME, HNAME
  REAL,          ALLOCATABLE, DIMENSION(:)   :: CMIN,   CMAX,   HYMIN,  HYMAX,  CDMIN,  CDMAX
END MODULE GRIDC
MODULE KINETICS
  REAL,    ALLOCATABLE, DIMENSION(:)             :: CoeffA_Turb, CoeffB_Turb,SECC_PAR,CGR,SSCS 
  REAL,    ALLOCATABLE, DIMENSION(:)             :: H2SQ10, H2S1DK, CH4Q10, CH41DK
  REAL,    ALLOCATABLE, DIMENSION(:)             :: CH4R, H2SR, FEIIR, MNIIR, SO4R
  REAL,    ALLOCATABLE, DIMENSION(:)             :: KFE_OXID, KFE_RED, KFEOOH_HalfSat, FeSetVel, KMN_OXID, KMN_RED, KMNO2_HalfSat, MnSetVel
  REAL,    ALLOCATABLE, DIMENSION(:)             :: BACTQ10, BACT1DK, BACTLDK, BACTS  !A_DISG, B_DISG, C_DISG              
  REAL,    POINTER,           DIMENSION(:,:)     :: N2, H2S, CH4, SO4, FEII, FEOOH, MNII, MNO2
  REAL,    POINTER,           DIMENSION(:,:)     :: WAGE, BACT, DGP
  REAL,    ALLOCATABLE,           DIMENSION(:)     :: ACHLA,  ORGP,   ORGN,   ORGC,   ORGSI,  O2OM,   ALGC,   PARTP
  REAL,    POINTER,               DIMENSION(:,:)   :: TN,     TP,     TKN
  REAL,    POINTER,               DIMENSION(:,:)   :: DON,    DOP,    DOC
  REAL,    POINTER,               DIMENSION(:,:)   :: PON,    POP,    POC
  REAL,    POINTER,               DIMENSION(:,:)   :: TON,    TOP,    TOC, APR
  REAL,    POINTER,               DIMENSION(:,:)   :: CHLA,   ATOT, NH3, TURB, TDG, pH, CO2, SECCHID, HCO3
  REAL,    POINTER,               DIMENSION(:,:)   :: O2DG,   CBODU,  TISS,   TOTSS, CO3
  REAL,    POINTER,               DIMENSION(:,:)   :: TDS,    COL,    NH4,    NO3,    PO4,    FE,     DSI,    PSI
  REAL,    POINTER,               DIMENSION(:,:)   :: LDOM,   RDOM,   LPOM,   RPOM,   O2,     TIC,    ALK
  REAL,    POINTER,               DIMENSION(:,:,:) :: SS,     ALG,    CG,     CBOD
  CHARACTER(8),      ALLOCATABLE, DIMENSION(:)     :: SROC, AVERTM
  REAL,              ALLOCATABLE, DIMENSION(:)       :: CBODS,sedcip,sedcin, sedcic,sedcis,volkti
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: ORGPLD,  ORGPRD,   ORGPLP,    ORGPRP,  ORGNLD,  ORGNRD, ORGNLP, ORGNRP,por,voli,nbodtot,pbodtot
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: LDOMPMP, LDOMNMP,  LPOMPMP,   LPOMNMP, RPOMPMP, RPOMNMP
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: LPZOOINP,LPZOOINN, LPZOOOUTP, LPZOOOUTN
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: SEDCAR,    SEDN, SEDP
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: SEDVPC,  SEDVPP, SEDVPN
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: SDKV,    SEDDKTOT
  REAL,    POINTER,               DIMENSION(:,:)     :: LDOMP,  RDOMP,  LPOMP,  RPOMP,  LDOMN,  RDOMN,  LPOMN,  RPOMN
  REAL,    POINTER,               DIMENSION(:,:)     :: LDOMPSS,  RDOMPSS, LPOMPSS, RPOMPSS, LDOMNSS, RDOMNSS
  REAL,    POINTER,               DIMENSION(:,:)     :: LPOMNSS,  RPOMNSS
  REAL,    POINTER,               DIMENSION(:,:)     :: LDOMPAP,  LDOMPEP, LPOMPAP, LPOMPNS, RPOMPNS
  REAL,    POINTER,               DIMENSION(:,:)     :: LDOMNAP,  LDOMNEP, LPOMNAP, LPOMNNS, RPOMNNS
  REAL,    POINTER,               DIMENSION(:,:)     :: SEDDP,    SEDASP,  SEDOMSP, SEDNSP,  LPOMEPP, SED
  REAL,    POINTER,               DIMENSION(:,:)     :: SEDDN,    SEDASN,  SEDOMSN, SEDNSN,  LPOMEPN, SEDNO3
  REAL,    POINTER,               DIMENSION(:,:)     :: SEDDC,    SEDASC,  SEDOMSC, SEDNSC,  LPOMEPC
  REAL,    POINTER,               DIMENSION(:,:)     :: CBODNS,   SEDCB,   SEDCBP,  SEDCBN,  SEDCBC
  REAL,    POINTER,               DIMENSION(:,:,:)   :: cbodp,    cbodn
  REAL,    POINTER,               DIMENSION(:,:,:)   :: cbodpss,    cbodnss
  REAL,    POINTER,               DIMENSION(:,:)     :: CBODNSp,  CBODNSn
  REAL,          ALLOCATABLE, DIMENSION(:)     :: CGQ10,  CG0DK,  CG1DK,  CGS,CGLDK,CGKLF,CGCS      ! SW 10/26/15
  REAL,          ALLOCATABLE, DIMENSION(:)     :: AG,     AR,     AE,     AM,     AS,     EXA,    ASAT
  REAL,          ALLOCATABLE, DIMENSION(:)     :: AHSP,   AHSN,   AHSSI,  AT1,    AT2,    AT3,    AT4
  REAL,          ALLOCATABLE, DIMENSION(:)     :: AK1,    AK2,    AK3,    AK4
  REAL,          ALLOCATABLE, DIMENSION(:)     :: ALGP,   ALGN,   ALGSI
  REAL,          ALLOCATABLE, DIMENSION(:)     :: EG,     ER,     EE,     EM,     EB
  REAL,          ALLOCATABLE, DIMENSION(:)     :: EHSN,   EHSP,   EHSSI,  ESAT,   EHS,    ENPR
  REAL,          ALLOCATABLE, DIMENSION(:)     :: EP,     EN,     EC,     ESI,    ECHLA,  EPOM
  REAL,          ALLOCATABLE, DIMENSION(:)     :: ET1,    ET2,    ET3,    ET4
  REAL,          ALLOCATABLE, DIMENSION(:)     :: EK1,    EK2,    EK3,    EK4
  REAL,          ALLOCATABLE, DIMENSION(:)     :: BODP,   BODN,   BODC                                                 !TC 01/15/02
  REAL,          ALLOCATABLE, DIMENSION(:)     :: LDOMDK, RDOMDK, LRDDK
  REAL,          ALLOCATABLE, DIMENSION(:)     :: LPOMDK, RPOMDK, LRPDK,  POMS,   APOM,   ANPR
  REAL,          ALLOCATABLE, DIMENSION(:)     :: OMT1,   OMT2,   OMK1,   OMK2
  REAL,          ALLOCATABLE, DIMENSION(:)     :: KBOD,   TBOD,   RBOD
  REAL,          ALLOCATABLE, DIMENSION(:)     :: PO4R
  REAL,          ALLOCATABLE, DIMENSION(:)     :: NH4R,   NH4DK
  REAL,          ALLOCATABLE, DIMENSION(:)     :: NH4T1,  NH4T2,  NH4K1,  NH4K2
  REAL,          ALLOCATABLE, DIMENSION(:)     :: NO3DK,  NO3S,FNO3SED
  REAL,          ALLOCATABLE, DIMENSION(:)     :: NO3T1,  NO3T2,  NO3K1,  NO3K2
  REAL,          ALLOCATABLE, DIMENSION(:)     :: DSIR,   PSIS,   PSIDK,  PARTSI
  REAL,          ALLOCATABLE, DIMENSION(:)     :: FER,    FES
  REAL,          ALLOCATABLE, DIMENSION(:)     :: CO2R
  REAL,          ALLOCATABLE, DIMENSION(:)     :: O2NH4,  O2AR,   O2AG
  REAL,          ALLOCATABLE, DIMENSION(:)     :: O2EG,   O2ER, KG_H2O_CONSTANT
  REAL, ALLOCATABLE, DIMENSION(:,:)            :: MIGON, MIGOFF, LOLD, TWQ   
  INTEGER, ALLOCATABLE, DIMENSION(:)           :: MIGRATE_GROUP, MIGRATE_MODEL, TS_DEC, DEPTH_LIM_ONOFF, NMINT, DEN_USE, DEPTH_CALC_ONOFF   
  INTEGER                                      :: ALGMIGRATION_DEBUG, NMIG, NITWQ, MIGI, KK
  CHARACTER(2)                                 :: MIGRATION  
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ASETTLE, DEN_AVG, DENP, DEN1, DEN2, ALLIM_OLD  ! CO 6/9/2019
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: DEN  ! CO 6/9/2019
  REAL, ALLOCATABLE, DIMENSION(:)     :: AMP, PHASE, C_COEFF_EXT, RAD, MIND, MAXD, DENSI, DENBI, T_DEC, C_DENINC, C_DENDEC, DEPTH_LIM, LOSS_FRAC  ! CO 6/3/2019
  REAL, ALLOCATABLE, DIMENSION(:)     :: I_C, C_DENINC_1, C_DENINC_2, C_DENDEC_1, C_DENDEC_2, DENP_MINS, DENP_MINB, DENP_MIN, DEN_COR, EXP_DEPTH    ! CO 6/10/2019
 
    END MODULE KINETICS
    MODULE EPIPHYTON
      REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: EPD
      CHARACTER(8),  ALLOCATABLE, DIMENSION(:,:)   :: EPIC
      INTEGER       :: NEP,    NEPT
    END MODULE EPIPHYTON

MODULE MACROPHYTEC
  REAL,    POINTER,               DIMENSION(:,:)     :: NH4MR,  NH4MG,  LDOMMAC, RPOMMAC, LPOMMAC, DOMP, DOMR, TICMC
  REAL,    POINTER,               DIMENSION(:,:)     :: PO4MR,  PO4MG
  REAL,              ALLOCATABLE, DIMENSION(:)       :: MG,     MR,     MM, MMAX,   MBMP
  REAL,              ALLOCATABLE, DIMENSION(:)       :: MT1,    MT2,    MT3,    MT4,    MK1,    MK2,    MK3,    MK4
  REAL,              ALLOCATABLE, DIMENSION(:)       :: MP,     MN,     MC
  REAL,              ALLOCATABLE, DIMENSION(:)       :: PSED,   NSED,   MHSP,   MHSN,   MHSC,   msat,   exm
  REAL,              ALLOCATABLE, DIMENSION(:)       :: CDSTEM, dmv, dwsa,anorm
  REAL,              ALLOCATABLE, DIMENSION(:)       :: ARMAC , sarea, iwind
  REAL,              ALLOCATABLE, DIMENSION(:)       :: O2MG,   O2MR,   LRPMAC,  MPOM
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: MACMBRS,MACMBRT,SSMACMB
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: CW,     BIC, macwbci,vstemkt
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: MACTRMR,MACTRMF,MACTRM,macrcvp,macrclp
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: MMR,    MRR
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: MAC,    MACT
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: MPLIM,  MNLIM, MCLIM
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: SMAC,   SMACT
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: GAMMAJ,vstem
  REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: MGR
  REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: MACRC,  MACRM
  REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: MLLIM
  REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: MACSS
  REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: SMACRC, SMACRM
  INTEGER                                            :: NMC, nmct
  LOGICAL,           ALLOCATABLE, DIMENSION(:)       :: KTICOL
  LOGICAL,           ALLOCATABLE, DIMENSION(:,:)     :: PRINT_MACROPHYTE, MACROPHYTE_CALC
  LOGICAL                                            :: MACROPHYTE_ON
  CHARACTER(8),      ALLOCATABLE, DIMENSION(:,:)     :: mprwbc, macwbc
  CHARACTER(10),      ALLOCATABLE, DIMENSION(:,:)    :: CONV2
  CHARACTER(10),     ALLOCATABLE, DIMENSION(:,:,:,:) :: MLFPR
  DATA                                                  SAVOLRAT /9000.0/, DEN /6.0E4/
END MODULE MACROPHYTEC
MODULE ZOOPLANKTONC
  LOGICAL                                            :: ZOOPLANKTON_CALC
  REAL,              ALLOCATABLE, DIMENSION(:)       :: zg,zm,zeff,PREFP,zr,ZOOMIN,ZS2P,EXZ,ZS
  REAL,              ALLOCATABLE, DIMENSION(:)       :: Zt1,Zt2,Zt3,Zt4,Zk1,Zk2,Zk3,Zk4
  REAL,              ALLOCATABLE, DIMENSION(:)       :: ZP,ZN,ZC,o2zr
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: PREFA, PREFZ ! OMNIVOROUS ZOOPLANKTON
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: po4zr,NH4ZR,DOZR,TICZR,LPZOOOUT,LPZOOIN
  REAL,    POINTER,               DIMENSION(:,:,:)   :: ZOO, ZOOSS
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: ZMU,TGRAZE,ZRT,ZMT
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: ZOORM,ZOORMR,ZOORMF
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: agzt
  REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: AGZ, ZGZ ! OMNIVOROUS ZOOPLANKTON
    END MODULE ZOOPLANKTONC
MODULE TDG_SYSTDG
  REAL                            :: TDGP1, TDGP2, TDGP3, TDGP4, TDGP12, TDGP22, TDGP32, TDGP42, TDGE1, TDGE2, TDGE12, TDGE22, ROP1, ROP2, ROP3, ROP4
  REAL                            :: TWE, TWCE, TWE_TS, FBE, QSPILL, TDGSPMN, GTPC(120)
  INTEGER                           :: POWNO, FLNO, NRO, NBAY
  INTEGER                           :: TWEMOD, TDGEQ, TDGROEQ, TDGENTEQ
  INTEGER                           :: ig
  CHARACTER(72)                       :: TITLESYSTDG(10)
  CHARACTER(8)                        :: TWETSC, TDGLOC, SYSTDGC, N2BNDC, DOBNDC, TDGTAC, GTTYP(120)
  CHARACTER(72)                       :: TWEFN 
  LOGICAL                             :: SYSTDG
ENDMODULE TDG_SYSTDG
MODULE HABITAT
integer :: ifish, n, nseg,iopenfish,kkmax,kseg,jjw
character*80, allocatable, dimension(:) :: fishname
character*80 :: conhab,conavg,consurf,consod
real, allocatable, dimension(:) :: fishtemph,fishtempl,fishdo,habvol,phabvol,cdo,cpo4,cno3,cnh4,cchla,ctotp,cdos,cpo4s,cno3s,cnh4s,cchlas,ctotps,cgamma,ssedd
integer, allocatable, dimension (:) :: isegvol
END MODULE HABITAT  
MODULE GRID
  REAL,          ALLOCATABLE, DIMENSION(:)     :: LAT, LONG, YNORTHING, XEASTING
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: UTMZONE, NPOINT, UHS, DHS, JBDN
  REAL,          ALLOCATABLE, DIMENSION(:)     :: DLX, PHI0
  INTEGER                                      :: NWB, NBR, NUP,NNBR,NCBP,NINTERNAL
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: B
  REAL*8, ALLOCATABLE, DIMENSION(:)            :: XC,YC,XNW,YNW,XW,YW,XSW,YSW,XSE,YSE,XNE,YNE,XE,YE,XCN,YCN,XCS,YCS
    END MODULE GRID

    Module CEMAVars
    USE PREC
    Integer(4), Allocatable, Dimension(:) :: ConsolidationType, ConstPoreWtrRate, NumCEMAPWInst
    Integer(4), Allocatable, Dimension(:) :: ConsRegSegSt, ConsRegSegEn, ConsolidRegnNum
    Character(256)                        :: ConsolidRateRegnFil
    Real(R8), allocatable, Dimension(:)   :: ConsolidRateTemp													  
    Real(R8), Allocatable, Dimension(:)   :: BedElevation, BedElevationLayer, BedPorosity
    Real(R8), Allocatable, Dimension(:,:) :: CellArea
    Real(R8), Allocatable, Dimension(:)   :: BedConsolidRate, PorewaterRelRate, ConstConsolidRate
    Real(R8), Allocatable, Dimension(:)   :: CEMACumPWRelease, CEMACumPWReleaseRate, CEMACumPWToRelease, CEMACumPWReleased
    Real(R8), Allocatable, Dimension(:,:) :: CEMASedConc,BubbleRelWB    ! SW 7/1/2017
    Real(R8), Allocatable, Dimension(:)   :: VOLCEMA
    Logical, Allocatable, Dimension(:)    :: CEMALayerAdded, CEMASSApplied
    Logical, Allocatable, Dimension(:) :: EndBedConsolidation, BedConsolidationSeg   ! cb 6/28/17
    Logical, Allocatable, Dimension(:) :: ApplyCEMAPWRelease
    
    Real(R8), Allocatable, Dimension(:) :: SDRegnPOC_T, SDRegnPON_T, SDRegnPOP_T, SDRegnSul_T
    Real(R8), Allocatable, Dimension(:) :: SDRegnPOC_L_Fr, SDRegnPOC_R_Fr, SDRegnPON_L_Fr
    Real(R8), Allocatable, Dimension(:) :: SDRegnPON_R_Fr, SDRegnPW_DiffCoeff, SDRegnOx_Threshold
    Real(R8), Allocatable, Dimension(:) :: SDRegnPOP_L_Fr, SDRegnPOP_R_Fr
    Real(R8), Allocatable, Dimension(:) :: SDRegnAe_NH3_NO3_L, SDRegnAe_NH3_NO3_H, SDRegnAe_NO3_N2_L
    Real(R8), Allocatable, Dimension(:) :: SDRegnAe_NO3_N2_H, SDRegnAn_NO3_N2, SDRegnAe_CH4_CO2
    Real(R8), Allocatable, Dimension(:) :: SDRegnAe_HS_NH4_Nit, SDRegnAe_HS_O2_Nit, SDRegn_Theta_PW,SDRegn_Theta_PM
    Real(R8), Allocatable, Dimension(:) :: SDRegn_Theta_NH3_NO3, SDRegn_Theta_NO3_N2, SDRegn_Theta_CH4_CO2
    Real(R8), Allocatable, Dimension(:) :: SDRegn_Sulfate_CH4_H2S, SDRegnAe_H2S_SO4, SDRegn_Theta_H2S_SO4
    Real(R8), Allocatable, Dimension(:) :: SDRegn_NormConst_H2S_SO4, SDRegn_MinRate_PON_Lab, SDRegn_MinRate_PON_Ref
    Real(R8), Allocatable, Dimension(:) :: SDRegn_MinRate_PON_Ine, SDRegn_MinRate_POC_Lab, SDRegn_MinRate_POC_Ref
    Real(R8), Allocatable, Dimension(:) :: SDRegn_MinRate_POC_Ine, SDRegn_Theta_PON_Lab, SDRegn_Theta_PON_Ref
    Real(R8), Allocatable, Dimension(:) :: SDRegn_Theta_PON_Ine, SDRegn_Theta_POC_Lab, SDRegn_Theta_POC_Ref
    Real(R8), Allocatable, Dimension(:) :: SDRegn_Theta_POC_Ine        ! cb 10/8/13    !Real(8), Allocatable, Dimension(:) :: SDRegn_Theta_POC_Ine, SDRegn_CH4CompMethod
    Real(R8), Allocatable, Dimension(:) :: SDRegn_Theta_POP_Lab, SDRegn_Theta_POP_Ref, SDRegn_Theta_POP_Ine
    Real(R8), Allocatable, Dimension(:) :: Kdp1, Kdp2, KdFe1, KdFe2, KdMn1, KdMn2, KdNH31, KdNH32, KdH2S1, KdH2S2
    Real(R8), Allocatable, Dimension(:) :: PartMixVel,BurialVel,POCr,delta_kpo41,DOcr,KsOxch
    Real(R8), Allocatable, Dimension(:) :: SDRegn_MinRate_POP_Lab, SDRegn_MinRate_POP_Ref, SDRegn_MinRate_POP_Ine
    Real(R8), Allocatable, Dimension(:) :: SD_NO3p2, SD_NH3p2, SD_NH3Tp2, SD_CH4p2, SD_PO4p2, SD_PO4Tp2
	  Real(R8), Allocatable, Dimension(:) :: SD_HSp2, SD_HSTp2
    Real(R8), Allocatable, Dimension(:) :: SD_poc2, SD_pon2, SD_pop2, SD_NH3Tp, SD_NO3p, SD_PO4Tp, SD_HSTp
    Real(R8), Allocatable, Dimension(:) :: SD_fpon, SD_fpoc, SD_kdiaPON, SD_ThtaPON, SD_kdiaPOC, SD_ThtaPOC
    Real(R8), Allocatable, Dimension(:) :: SD_kdiaPOP, SD_ThtaPOP, SD_NH3T, SD_PO4, SD_FPOP
    Real(R8), Allocatable, Dimension(:) :: SD_JPOC, SD_JPON, SD_JPOP,  SD_TDS
    
    Real(R8), Allocatable, Dimension(:) :: SD_Denit, SD_JDenit, SD_JO2NO3,  SD_HS   ! cb 7/26/18
    Real(R8), Allocatable, Dimension(:) :: SD_Fe2 
    Real(R8), Allocatable, Dimension(:) :: SD_Mn2 
    Real(R8), Allocatable, Dimension(:) :: SD_pHValue   ! cb 7/26/18    !Real(8), Allocatable, Dimension(:) :: SD_SO4Conc, SD_pHValue
    Real(R8), Allocatable, Dimension(:) :: SD_EPOC, SD_EPON, SD_EPOP
    Real(R8), Allocatable, Dimension(:) :: SD_AerLayerThick
        
    Real(R8), Allocatable, Dimension(:,:,:) :: MFTSedFlxVars, CEMA_SD_Vars
    Real(R8), Allocatable, Dimension(:,:) :: CEMATSSCopy
    
    Integer(4), Allocatable, Dimension(:) :: CEMAMFT_RandC_RegN, CEMAMFT_InCond_RegN
    Integer(4), Allocatable, Dimension(:) :: SedBedInitRegSegSt, SedBedInitRegSegEn
    Integer(4), Allocatable, Dimension(:) :: SedBedDiaRCRegSegSt, SedBedDiaRCRegSegEn
    
    Integer(4), Allocatable, Dimension(:) :: FFTActPrdSt, FFTActPrdEn
    Integer, Allocatable, Dimension(:) :: SDRegn_CH4CompMethod, SDRegn_POMResuspMethod
    Real(R8), Allocatable, Dimension(:) :: FFTLayConc
    
    Real(R8), Allocatable, Dimension(:) :: H2SDis, H2SGas, CH4Dis, CH4Gas, NH4Dis, NH4Gas, CO2Dis, CO2Gas
    Real(R8), Allocatable, Dimension(:) :: BubbleRadiusSed, PresBubbSed, PresCritSed
    Real(R8), Allocatable, Dimension(:) :: CgSed, C0Sed, CtSed
    Real(R8), Allocatable, Dimension(:,:,:):: TConc, TConcP, SConc
    Real(R8), Allocatable, Dimension(:,:,:):: DissolvedGasSediments
    Integer(8), Allocatable, Dimension(:) :: MFTBubbReleased, LastDiffVolume
    Integer(4), Allocatable, Dimension(:,:) :: BubblesLNumber, BubblesStatus
    Real(R8), Allocatable, Dimension(:,:) :: BubblesRadius, BubblesRiseV, BubblesCarried
    Real(R8), Allocatable, Dimension(:,:,:) :: BubblesGasConc, BRVoluAGas, BRRateAGas
    Real(R8), Allocatable, Dimension(:,:) :: BubblesReleaseAllValue, BRRateAGasNet 
    Real(R8), Allocatable, Dimension(:) :: BottomTurbulence
    Logical, Allocatable, Dimension(:) :: CrackOpen
    Logical, Allocatable, Dimension(:,:) :: FirstBubblesRelease, BubblesAtSurface
    
    
    Integer(4) :: CEMAFilN, NumConsolidRegns, CEMASedimentType
    Integer(4) :: CEMASNPOutFilN=2411, CEMATSR1OutFilN=2412, SegNumI, LayerNum
    Integer(4) :: CEMABtmLayFilN=2414, TempCntr1
    Integer(4) :: CEMASedFlxFilN1=2415, CEMASedFlxFilN2=2416, CEMASedFlxFilN3=2417, CEMALogFilN=2418, CEMASedFlxFilN4=2419
    Integer(4) :: CEMASedFlxFilN5=2420, CEMASedFlxFilN6=2421, CEMASedFlxFilN7=2422, CEMASedFlxFilN8=2423, CEMASedFlxFilN9=2424,CEMASedFlxFilN10=2441,CEMASedFlxFilN11=2442
    Integer(4) :: CEMASedFlxFilN12=3483,CEMASedFlxFilN13=3484,CEMASedFlxFilN14=3485,CEMASedFlxFilN15=3486,CEMASedFlxFilN16=3487,CEMASedFlxFilN17=3488,CEMASedFlxFilN18=3489
    Integer(4) :: CEMASedFlxFilN19=3490,CEMASedFlxFilN20=3491,CEMASedFlxFilN21=3492,CEMASedFlxFilN22=3493
    Integer(4) :: CEMASedFlxFilN23=3494,CEMASedFlxFilN24=3495,CEMASedFlxFilN25=3496,CEMASedFlxFilN26=3497,CEMASedFlxFilN27=3498,CEMASedFlxFilN28=3499
    Integer(4) :: CEMASedFlxFilN29=3500,CEMASedFlxFilN30=3501,CEMASedFlxFilN31=3502
    Integer(4) :: CEMAOutFilN1=2426, CEMAOutFilN2=2429, CEMAOutFilN3=2435, CEMAOutFilN4=2437,CEMAOutFilBub=2440
    Integer(4) :: CEMAOutFilN5=2438, CEMAOutFilN6=2439
    Integer(4) :: NumRegnsSedimentDiagenesis, NumRegnsSedimentBedComposition
    Integer(4) :: NumFFTActivePrds, FFTActPrd, NumBubRelArr, NumGas=4
    Real(R8) :: LayerAddThkFrac, BedElevationInit, BedPorosityInit, CEMASedimentDensity
    Real(R8) :: CEMAParticleSize, TotalPoreWatVolume, TotalSedimentsInBed, TotalPoreWatRemoved
    Real(R8) :: CEMASedimentSVelocity, CEMAPWpH
    Real(R8) :: NH4_NH3_Eqb_Const, HS_H2S_Eqb_Const,VolumeIncreasedConsolid
    Real(R8) :: GasConst_R=0.0821    !L.atm/mol/K
    Real(R8) :: HenryConst_NH3, HenryConst_CH4, HenryConst_H2S, HenryConst_CO2
    Real(R8) :: InitFFTLayerConc, FFTLayerSettVel
    Real(R8) :: GasDiff_Sed, CalibParam_R1, YoungModulus, CritStressIF
    Real(R8) :: BubbRelScale, CrackCloseFraction, MaxBubbRad, BubbRelFraction, BubbAccFraction
    Real(R8) :: BubbRelFractionAtm, BubbWatGasExchRate
    Real(R8) :: CEMATurbulenceScaling
    Real(R8) :: IceThicknessChange  ! cb 2/5/13
    Real(R8) :: TAUCRPOM, crshields, spgrav_POM, dia_POM,  GasReleaseCH4  !,GasReleaseCO2   ! SW 10/10/2017   ! SW 10/19/2017
    REAL(R8) :: NXTSEDIAG,SEDIAGFREQ   ! SW 5/25/2017  
    
    Logical CEMARelatedCode, IncludeBedConsolidation, IncludeCEMASedDiagenesis, IncludeFFTLayer, FFTActive, FirstTimeInFFTCode
    Logical IncludeIron, IncludeManganese, IncludeDynamicpH, IncludeAlkalinity, SD_global                                         ! cb 5/22/15
    Logical CEMASedimentProcessesInc, WriteBESnp, WritePWSnp, WriteCEMAMFTSedFlx, CEMA_POM_Resuspension
    Logical FirstTimeinCEMAMFTSedDiag, MoveFFTLayerDown
    Logical LimBubbSize, UseReleaseFraction, FirstTimeInBubbles, ApplyBubbTurb
    LOGICAL sediment_diagenesis, cao_method, Bubbles_Calculation  ! SW 10/20/2017
    
    Real(R8), Allocatable, Dimension(:) :: SDRegnH2S_T, SDRegnNH3_T, SDRegnCH4_T, SDRegnNO3_T        
    Real(R8), Allocatable, Dimension(:) :: SDRegnTIC_T, SDRegnALK_T, SDRegnPO4_T
    Real(R8), Allocatable, Dimension(:) :: SDRegnFe2_T, SDRegnFeOOH_T, SDRegnMn2_T, SDRegnMnO2_T
    Real(R8), Allocatable, Dimension(:) :: SDRegnT_T, SDRegnpH
    Real(R8), Allocatable, Dimension(:) :: SDPFLUX, SDNH4FLUX,SDNO3FLUX
        
    End Module CEMAVars
    MODULE ENVIRPMOD
CHARACTER*3, SAVE, ALLOCATABLE, DIMENSION (:) :: CC_E,CD_E
CHARACTER*3, SAVE :: VEL_VPR,TEMP_VPR,SELECTIVEC,DEPTH_VPR
REAL,ALLOCATABLE, SAVE, DIMENSION (:) :: C_CNT,CD_CNT,C_TOT,CD_TOT,T_CLASS,V_CLASS,C_SUM,CD_SUM,C_AVG,CD_AVG,D_CLASS
REAL,ALLOCATABLE, SAVE, DIMENSION (:,:) :: C_CLASS,CD_CLASS,CONC_C,CONC_CD
REAL, ALLOCATABLE, SAVE, DIMENSION (:) :: C_INT,C_TOP,CD_INT,CD_TOP,CN_E,CDN_E
INTEGER NUMCLASS,NAC_E,NACD_E,JJ,JACD,ISTART,IEND
    END MODULE ENVIRPMOD
    MODULE ALGAE_TOX
    INTEGER, PARAMETER                                 :: NUMATOXINS=4
    INTEGER                                            :: NATS, NATE, ATOXIN_DEBUG_FN=2501
    LOGICAL                                            :: ALGAE_TOXIN_FILE, ALGAE_TOXINS 
    REAL, ALLOCATABLE, DIMENSION (:,:)                 :: CTP,CTB
    REAL, DIMENSION (:)                                :: CTREL(NUMATOXINS),CTD(NUMATOXINS)
    CHARACTER(2)                                       :: ATOX, ATOX_DEBUG
    END MODULE ALGAE_TOX

    

!***********************************************************************************************************************************
!**                                          C E - Q U A L - W 2  P R E P R O C E S S O R                                         **
!***********************************************************************************************************************************

PROGRAM W2_PRE
  USE IFPORT, only: chdir,GETDRIVEDIRQQ,FILE$CURDRIVE   
  USE IFQWIN
  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_NAN    ! FOR NAN CHECK
  !USE DFLIB ! SW 11/9/2017 USED IFPORT INSTEAD
  USE BUILDVERSION   ! SW 3/2015
  USE CONTROL; USE GRIDC; USE KINETICS; USE PREC; USE FILE_ERROR; use macrophytec; use zooplanktonc; USE HABITAT; USE CEMAVars
  USE GRID
  USE EPIPHYTON; USE ENVIRPMOD; USE TDG_SYSTDG; USE ALGAE_TOX
! Variable declarations

  TYPE (WINDOWCONFIG) WC
  TYPE (QWINFO)       WINFO
  CHARACTER(360):: MODDIR, DIRC   ! SW 11/9/2017 If you are likely to encounter a long path, allocate a buffer of size $MAXPATH ($MAXPATH = 260). USE IFPORT CHARACTER($MAXPATH) dir dir = FILE$CURDRIVE length = GETDRIVEDIRQQ(dir)
  REAL          :: JDAY,   JDAYO, BCONSTRICTION   ! SW 8/5/2018
  REAL          :: WAIT_TIME, TIME_BUFFER,CO2ATMPPM  ! SW 9/19/019
  LOGICAL       :: CONSTRICTION
  LOGICAL   :: LAKE_RIVER_CONTOURC       ! SW 2/27/2020
  CHARACTER(2) :: LAKE_RIVER_CONTOUR_ON
  CHARACTER(80):: FILE_LAKE_CONTOUR_T(20),FILE_LAKE_CONTOUR_DO(20),FILE_RIVER_CONTOUR_T(20),FILE_RIVER_CONTOUR_DO(20)
  INTEGER      :: NUM_LAKE_CONTOUR,NUM_RIVER_CONTOUR,LAKE_CONTOUR_SEG(20),RIVER_CONTOUR_BR1(20),RIVER_CONTOUR_BR2(20)
  INTEGER      :: JW_RIVER_CONTOUR(20), JW_LAKE_CONTOUR(20),LAKE_CONTOUR_FORMAT,RIVER_CONTOUR_FORMAT
  REAL         :: LAKE_CONTOUR_START(20),LAKE_CONTOUR_FREQ(20),RIVER_CONTOUR_START(20),RIVER_CONTOUR_FREQ(20),NXT_RIVER_CONTOUR(20), NXT_LAKE_CONTOUR(20)
  INTEGER       :: ISEG, ICON,LAKE_RIVER_CONTOUR=9540
  INTEGER       :: nbodcs, nbodce, nbodps, nbodpe, nbodns, nbodne, ibod, jcb
  INTEGER       :: CON,    WRN,    NPROC, NPOLY   !   NBR, NPROC   SW 4/6/15   ERR,    ! 7/14/16   NPOLY=1 THEN LINEAR SYSTEM - SET LAT/LONG AT SEG2 FIRST WB, ELSE DOWNSTREAM DAM - SET LAT/LONG AT DAM
  INTEGER       :: YEAR,   CUF, EUHF,TUHF,CUHF,EDHF,TDHF,CDHF,INUM,IMWB
  INTEGER, DIMENSION(:) :: INUMMWB(20)
  INTEGER       ::  NTR,    NSP,    NWD,    NPU,    NGT,    NPI    ! NWB,   
  INTEGER       :: ICST, IOK=0
  INTEGER       :: DOC_DER,POC_DER,TOC_DER,DON_DER,PON_DER,TON_DER,TKN_DER,TN_DER,NH3_DER,DOP_DER,POP_DER,TOP_DER,TP_DER,APR_DER
  INTEGER       :: CHLA_DER,ATOT_DER,O2DG_DER,TDG_DER,TURB_DER,TOTSS_DER,TISS_DER,CBODU_DER,PH_DER,CO2_DER,HCO3_DER,CO3_DER,SECCHI_DER
  INTEGER(I2)   :: X=320,  Y=320,  ONE=1
  CHARACTER(1)  :: ESC, ICHAR1,ICHAR2, CHAR1,INFORMAT1
  CHARACTER(2)  :: DEG,IIN,INFORMAT2,MWB
  CHARACTER(2), DIMENSION(:) :: MWB2(20)
  CHARACTER(3)  :: BRA,    TRA 
  CHARACTER(5)  :: FMTA
  CHARACTER(7)  :: BLANK
  CHARACTER(8)  :: AID,     CCC,    LIMC, CLOSEC,SELECTC,HABTATC,ENVIRPC,AERATEC,INITUWL,tcyearly,tspltcntr,tspltc,dltintr,GT2CHAR,dynsel,tcelevcon,elcont,aid1
  CHARACTER(8)  :: tsdynsel,elcontspl,tsshare,ORGCC, SED_DIAG
  CHARACTER(8)  :: RSIC,   RSOC,   TSRC,   WDOC, WLC,  FLOWBALC, NPBALC, CO2ATM_REGRESSION
  CHARACTER(11) :: FMTI
  CHARACTER(13) :: FMTF
  CHARACTER(18) :: FMTFI,  FMTIF,  FMT2FI, FMTF2I
  CHARACTER(21) :: FMT2IF
  CHARACTER(20) :: CCONSTIT(150)
  CHARACTER(389) :: WINTITLE                                                                                            !SW 10/15/00
  CHARACTER(72) :: BTHFN1, METFN1, QOTFN1
  CHARACTER(72) :: QINFN1, QTRFN1, QDTFN1, PREFN1, EUHFN1, EDHFN1
  CHARACTER(72) :: TINFN1, TTRFN1, TDTFN1, TPRFN1, TUHFN1, TDHFN1
  CHARACTER(72) :: CINFN1, CTRFN1, CDTFN1, CPRFN1, CUHFN1, CDHFN1
  CHARACTER(72) :: SNPFN1, PRFFN1, VPLFN1, CPLFN1, SPRFN1, FLXFN1
  CHARACTER(72) :: QWDFN,  RSIFN,  QGTFN,  TSRFN,  WDOFN
  CHARACTER(72) :: SHDFN,  WSCFN, blank_check                                                                                       !SW 04/03/02
  CHARACTER(72) :: CONFN='w2_con.npt',INIFN='pre.opt',WRNFN='pre.wrn',ERRFN='pre.err', conaer
  CHARACTER(240):: TITLE240, title80
  CHARACTER(240), DIMENSION(:) :: FILEDIRMWB(20)
  CHARACTER(80) :: HEADER(3),SEGNUM
  LOGICAL       :: TRIBUTARIES,      WITHDRAWALS, MULTIPLEWB, ALGAE_MIGRATION_FILE
  LOGICAL       :: PIPES,            GATES,            SPILLWAYS,   PUMPS
  LOGICAL       :: DISSOLVED_SOLIDS, ALGAE,            CBO_DEMAND,  GENERIC_CONST
  LOGICAL       :: LABILE_DOM,       REFRACTORY_DOM,   LABILE_POM,  REFRACTORY_POM
  LOGICAL       :: PHOSPHORUS,       AMMONIUM,         NITRATE,     IRON
  LOGICAL       :: DISSOLVED_OXYGEN, SUSPENDED_SOLIDS, TOT_IC
  LOGICAL       :: ACTIVE_SEGMENT,   RESTART_IN
  LOGICAL       :: STATUS,           RESULT
  INTEGER(4)    :: length,istatus,tcnelev,tciseg,tspltjb,numtsplt,tempn,jstsplt(11),priority(11)  ! SW 3.72
  REAL          :: depth(11), minfrac(11),minhead(11),maxhead(11),maxflow(11)                     ! SW 3.72
  INTEGER(4), ALLOCATABLE, DIMENSION(:)     :: TCJB,TCJS
  CHARACTER(8), ALLOCATABLE,DIMENSION(:)    :: TCNTR

  real :: tcelev(20),tempcrit(50),XX
  real :: nxqin1,nxqin2,qinnx,qind,qino,qratio,qin1,qin2,jday1,jday2                                                                  ! cb 9/7/10
  REAL :: VOLG1,  QOTAVW1,  QINAVW1, tspltfreq, tsconv, NPBALF, WLF, FLOWBALF
  REAL, ALLOCATABLE, DIMENSION(:) ::  SAGR1, CVGR1, HGR1, BGR1
  real, allocatable, dimension(:) ::   sromax,sromin,tairmax,tairmin,tdewmax,tdewmin,windmax,windmin,phimin,phimax,cloudmax,TDEWAV,cloudmin,sroav,tairav,phiav,windav,cloudav
  INTEGER, ALLOCATABLE, DIMENSION(:) ::  NCCGR1                                              ! SW 11/1/13

  Logical SkipLoop, file_exists, DIAG_EXISTS
  Character(256) MessageTemp
  CHARACTER(2) :: PARTON,ALINE,DXTHEORYC
  REAL :: OUTFREQP
  INTEGER :: ILINEAR,HTSTBOT,HTSTSIDE,IDEBUG,NFISHSEG,NFISHPCEL

! Allocation declarations

  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: NACATD
  REAL,          ALLOCATABLE, DIMENSION(:)     :: ATMDEP_P,ATMDEP_N
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: ATM_DEPOSITIONC,ATM_DEPOSITION_INTERPOLATION
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:,:)   :: C_ATM_DEPOSITION
  INTEGER,       ALLOCATABLE, DIMENSION(:,:)   :: ATMDCN
  CHARACTER(72), ALLOCATABLE, DIMENSION(:)     :: ATMDEPFN  
  REAL(R8),      ALLOCATABLE, DIMENSION(:,:)   :: TSSUH1, TSSDH1, ATM_DEP_LOADING,ATM_DEP_AVG,ATM_DEP_MIN,ATM_DEP_MAX
  REAL,          ALLOCATABLE, DIMENSION(:)     :: ESP,    A1SP,   B1SP,   A2SP,   B2SP,   ETUSP,  EBUSP,  ETDSP,  EBDSP, TINMAX, TINMIN,TRTMAX, TRTMIN
  REAL,          ALLOCATABLE, DIMENSION(:)     :: AGASSP, BGASSP, CGASSP
  REAL,          ALLOCATABLE, DIMENSION(:)     :: EUPI,   EDPI,   WPI,    DLXPI,  FPI,    FMINPI, ETUPI,  EBUPI,  ETDPI,  EBDPI
  REAL,          ALLOCATABLE, DIMENSION(:)     :: EPU,    STRTPU, ENDPU,  EONPU,  EOFFPU, QPU,    ETPU,   EBPU
  REAL,          ALLOCATABLE, DIMENSION(:)     :: EGT,    ETUGT,  EBUGT,  ETDGT,  EBDGT, EGT2
  REAL,          ALLOCATABLE, DIMENSION(:)     :: GTA1,   GTB1,   GTA2,   GTB2
  REAL,          ALLOCATABLE, DIMENSION(:)     :: A1GT,   B1GT,   G1GT,   A2GT,   B2GT,   G2GT,    BGATE
  REAL,          ALLOCATABLE, DIMENSION(:)     :: AGASGT, BGASGT, CGASGT
  REAL,          ALLOCATABLE, DIMENSION(:)     :: ETTR,   EBTR,   EWD, TDH, TUH
  REAL,          ALLOCATABLE, DIMENSION(:)     :: ICEI,   ALBEDO, HWI,    BETAI,  GAMMAI, ICEMIN, ICET2,  EXH2O,  BETA,   EXOM, Z0    ! SW 3.6
  REAL,          ALLOCATABLE, DIMENSION(:)     :: EXSS,   SSS,    COLQ10, COLDK,  DX,     CBHE, TAUCR         !TC 08/20/03 SW 1/16/04
  REAL,          ALLOCATABLE, DIMENSION(:)     :: AFW,    BFW,    CFW,    WINDH,  TSED,   FI,     TSEDF, WSC
  REAL,          ALLOCATABLE, DIMENSION(:)     :: SEDCI,  SEDDK,  SODT1,  SODT2,  SODK1,  SODK2,  FSOD,   FSED, seds,sedbr  !SW 6/1/07
  REAL,          ALLOCATABLE, DIMENSION(:)     :: AX,     AZMAX,  T2I,    ELBOT,  SLOPE,SLOPEC   !   LAT,    LONG,  SW 4/6/15
  REAL,          ALLOCATABLE, DIMENSION(:)     :: ELWS,   SHADE,   FRIC  !   REAL,          ALLOCATABLE, DIMENSION(:)     :: DLX,    ELWS,   SHADE,  PHI0,   FRIC  SW 4/6/15
  REAL,          ALLOCATABLE, DIMENSION(:)     :: XBR,    XGR,    XI
  REAL,          ALLOCATABLE, DIMENSION(:)     :: QTRAVB, QTRMXB, QINAV,  QPRAV,  QDTAV,  QSTRAV
  REAL,          ALLOCATABLE, DIMENSION(:)     :: QSTR,   QSTRO
  REAL,          ALLOCATABLE, DIMENSION(:)     :: QINS,   QPRS,   QDTS,   QSTRS
  REAL,          ALLOCATABLE, DIMENSION(:)     :: QINMXB, QOTMXB, QINAVB, QOTAVB, QWDMXB, QWDAVB,  QDNAVB
  REAL,          ALLOCATABLE, DIMENSION(:)     :: QINMX,  QDTMX,  QSTRMX, QPRMX
  REAL,          ALLOCATABLE, DIMENSION(:)     :: QINAVW, QOTAVW, QINMXW, QOTMXW
  REAL,          ALLOCATABLE, DIMENSION(:)     :: RSOD,   DLTD,   TSRKTD
  REAL,          ALLOCATABLE, DIMENSION(:)     :: RSOF,   DLTF,   TSRKTF
  REAL,          ALLOCATABLE, DIMENSION(:)     :: DLTMAX, THETA
  REAL,          ALLOCATABLE, DIMENSION(:)     :: CIN,    CIN2, CPR                                                         !SW 01/07/01 4/11/2018

  REAL,          ALLOCATABLE, DIMENSION(:)     :: RCOEF1, RCOEF2, RCOEF3, RCOEF4,DGPO2,MINKL
  REAL,          ALLOCATABLE, DIMENSION(:)     :: SOD,ARODI,STRCKLR,BOUNDFR,AZE
  REAL,          ALLOCATABLE, DIMENSION(:)     :: ETSR
  REAL,          ALLOCATABLE, DIMENSION(:)     :: TSRD,   TSRF,   WDOD,   WDOF
  REAL,          ALLOCATABLE, DIMENSION(:)     :: TTLB,   TTRB,   CLLB,   CLRB,   SRLB1,  SRRB1,  SRLB2,  SRRB2,  ANG
  REAL,          ALLOCATABLE, DIMENSION(:)     :: SRFJDAY1, SRFJDAY2, SHADEIN,cinnload,cinpload,ctrnload,ctrpload
  REAL,          ALLOCATABLE, DIMENSION(:)     :: TDHMIN,TDHMAX,TUHMIN,TUHMAX
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: CUH,CDH,CUHMIN,CUHMAX,CUHAVG,CDHMIN,CDHMAX,CDHAVG
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: EPICI
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: TOPO                                                                 !SW 04/03/02
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: JDQIN,  JDQDT,  JDQTR,  JDQPR,  JDQWD
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: QIN,    QDT,    QTR,    QPR,    QWD
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: C2IWB   !   REAL,          ALLOCATABLE, DIMENSION(:,:)   :: B,      C2IWB  SW 4/6/15
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: SAGR,   SABR,   CVGR,   CVBR,   HGR,    HBR,    BGR,    BBR
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: ESTR,   WSTR
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: SNPD,   SCRD,   PRFD,   SPRD,   CPLD,   VPLD,   FLXD
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: SNPF,   SCRF,   PRFF,   SPRF,   CPLF,   VPLF,   FLXF
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: QTRS,   QTRAV,  QTRMX,  QTRMXS
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: QWDS,   QWDAV,  QWDMX
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: CINAVG, CINMIN, CINMAX, CTRAVG, CTRMIN, CTRMAX,cinload,ctrload ,cdtload,cdtavg,cdtmin,cdtmax                  !SW 01/07/01
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: TVP
  REAL,          ALLOCATABLE, DIMENSION(:,:,:) :: CVP,    EPIVP
  REAL(R8),      ALLOCATABLE, DIMENSION(:)     :: VOLB,   VOLG
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: BTH,    VPR,    LPR,    JJS
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: NSNP,   NSCR,   NSPR,   NPRF,   NVPL,   NFLX,   NCPL               !TC 06/11/02
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: NQIN,   NQDT,   NQTR,   NQPR,   NQOT
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: NISNP,  NIPRF,  NISPR,  KBMAX
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: IUPI,   IDPI,   KTUPI,  KBUPI,  KTDPI,  KBDPI
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: KTUSP,  KBUSP,  KTDSP,  KBDSP,  IUSP,   IDSP
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: KTUGT,  KBUGT,  KTDGT,  KBDGT,  IUGT,   IDGT
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: IUPU,   IDPU,   KTPU,   KBPU
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: EQSP,   EQGT,   NEQN,   ANEQN,  ENEQN
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: KTI
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: UQB,    DQB,    NSTR, NL     ! UHS,    DHS,    NPOINT   JBDN,   SW 4/8/2015
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: NTBR,   NWBR,   NUPIBR, NDPIBR, NUGTBR, NDGTBR
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: JBTR,   JBWD,   JBUPI,  JBDPI,  JBUGT,  JBDGT
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: IWD,    KTWD,   KBWD
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: IWR,    KTWR,   KBWR
  REAL,          ALLOCATABLE, DIMENSION(:)     :: EKTWR, EKBWR
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: nbodc, nbodn, nbodp
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: ITR,FBC,IBR
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: IWDO,   ITSR
  INTEGER,       ALLOCATABLE, DIMENSION(:,:)   :: INCN,   TRCN, incdt                                                         !SW 01/07/01
  INTEGER,       ALLOCATABLE, DIMENSION(:,:)   :: ISNP,   IPRF,   ISPR
  INTEGER,       ALLOCATABLE, DIMENSION(:,:)   :: NCCGR,  NCCBR
  INTEGER,       ALLOCATABLE, DIMENSION(:,:)   :: KTSTR,  KBSTR
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: EXC,    EXIC                                                         !SW 12/04/01
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: VISC,   CELC,   FETCHC, DYNSTRUC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: AZC,    AZSLC,  WTYPEC, BK,     SLICEC, REAERC, DTRC, TKECAL,gridcc
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: GASSPC, GASGTC, DYNPUMP, DYNPIPE
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: CAC,    CNAME2
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: PTRC,   PPUC,   PUGTC,  PDGTC,  PDSPC,  PUSPC,  PDPIC,  PUPIC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: ICEC,   RHEVC,  SEDC,   SDC,    PRNSC, DYNSEDK
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: SNPC,   SCRC,   SPRC,   VPLC,   CPLC,   PRFC,   FLXC, TCPL
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: QINIC,  TRIC,   DTRIC,  WDIC,   HDIC,   METIC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: SLTRC,  SLHTC,  FRICC, DLTADD
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: WINDC,  QINC,   QOUTC,  HEATC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: VBC,    MBC,    EBC,    PQC,    EVC,    PRC
  CHARACTER(72), ALLOCATABLE, DIMENSION(:)     :: BTHFN,  METFN,  VPRFN,  LPRFN,  QOTFN
  CHARACTER(72), ALLOCATABLE, DIMENSION(:)     :: QINFN,  QTRFN,  QDTFN,  PREFN,  EUHFN,  EDHFN
  CHARACTER(72), ALLOCATABLE, DIMENSION(:)     :: TINFN,  TTRFN,  TDTFN,  TPRFN,  TUHFN,  TDHFN
  CHARACTER(72), ALLOCATABLE, DIMENSION(:)     :: CINFN,  CTRFN,  CDTFN,  CPRFN,  CUHFN,  CDHFN
  CHARACTER(72), ALLOCATABLE, DIMENSION(:)     :: SNPFN,  PRFFN,  VPLFN,  CPLFN,  SPRFN,  FLXFN
  CHARACTER(72), ALLOCATABLE, DIMENSION(:)     :: EXTFN                                                                !TC 12/12/01
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: HPLTC,  CPLTC,  CDPLTC            ! SW 1/16/04
  CHARACTER(4),  ALLOCATABLE, DIMENSION(:)     :: CUNIT1                            ! SW 1/16/04
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: SEDRC                             ! SW 1/16/04
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: LATGTC, LATSPC, LATPIC, LATPUC, DYNGTC,GTIC

  CHARACTER(8),  ALLOCATABLE, DIMENSION(:,:)   :: SINKC,  STRIC,  CINBRC, CINTRC, CDTBRC, CPRBRC, CFWBC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:,:)   :: EPIPRC
  LOGICAL                                      :: CONSTITUENTS,     BRANCH_FOUND, FISH_PARTICLE_EXIST, CSVFORMAT, TECPLOTBR, ALGAE_MIGRATION                                      !SW 10/17/02
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: DIST_TRIB, ZERO_SLOPE, PUMP_DOWNSTREAM                                                !SW 10/17/02
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: SHIFT_DEMAND, SEDIMENT_CALC, PRECIPITATION
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: UH_EXTERNAL,  DH_EXTERNAL,   UQ_EXTERNAL, DQ_EXTERNAL, DAM_FLOW
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: UH_INTERNAL,  DH_INTERNAL,   UQ_INTERNAL, DQ_INTERNAL
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: OPEN_VPR,     OPEN_LPR, OPEN_LPRC
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: ISO_TEMP,     VERT_TEMP,     LONG_TEMP,   LONG_SED
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: INTERP_INFLOW                                                       ! cb 9/7/10
  LOGICAL,       ALLOCATABLE, DIMENSION(:,:)   :: ISO_CONC,     VERT_CONC,     LONG_CONC
  LOGICAL,       ALLOCATABLE, DIMENSION(:,:)   :: ISO_EPIPHYTON, VERT_EPIPHYTON, LONG_EPIPHYTON, EPIPHYTON_CALC, ISO_MACRO, VERT_MACRO, LONG_MACRO
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: ISO_SEDIMENT,   VERT_SEDIMENT,  LONG_SEDIMENT
  INTEGER,       ALLOCATABLE, DIMENSION(:)  :: IDUM
  REAL,          ALLOCATABLE, DIMENSION(:)  :: DDUM,FDUM
  INTEGER :: NDUM, JTEC, NTR1
  CHARACTER(8):: CDUM, CDNAME2
  CHARACTER(14)  :: KFNAME2
  LOGICAL :: DYNAMIC_SD



! Data declarations

  DATA BLANK /'       '/
  DATA DMO2  /2.04E-9/,  DMCO2 /1.63E-9/
  DATA CON   /21/, NPT /22/, WRN /371/, NEW /372/    !ERR /370/, 
  DATA NDC   /27/, NHY /15/, NFL /72/
  DATA IANGLES /18/         !SW 04/03/02
  
  DYNAMIC_SD=.FALSE.
  
! Intel Fortran specific code
OPEN(10,FILE='PREW2CodeCompilerVersion.opt',status='unknown')
write(10,*)'Preprocessor Compiler Version and Preprocessor Code Compile Date'
write(10,*)'INTEL_COMPILER_VERSION:',INTEL_COMPILER_VERSION
write(10,*)'INTEL_COMPILER_BUILD_DATE:',INTEL_COMPILER_BUILD_DATE
write(10,*)'PreW2 Version compile date:',BUILDTIME
close(10)

LENGTH=0
call get_command_argument(1,dirc,length,istatus)
dirc=TRIM(dirc)

if(length /= 0)then
    istatus=CHDIR(dirc)
    select case(istatus)
      case(2)  ! ENOENT
        write(WIN,*)'Error: Command argument not blank and the directory does not exist:',dirc
        write(ERR,*)'Error: Command argument not blank and the directory does not exist:',dirc
      case(20)   ! ENOTDIR
        write(WIN,*)'Error: Command argument not blank and is not a directory:', dirc
        write(ERR,*)'Error: Command argument not blank and is not a directory:', dirc
      case(0)    ! no error
    end select
endif

moddir = FILE$CURDRIVE              !  Get current directory
length = GETDRIVEDIRQQ(moddir)

  STAT = ABOUTBOXQQ ('CE-QUAL-W2 Preprocessor V4.5 - PSU Dept. Civil & Environmental Engineering S. Wells'C)   
  WRITE (WINTITLE,'(A29,A360)') 'CE-QUAL-W2 4.5 preprocessor. ',moddir                  !A260)') 'CE-QUAL-W2 4.1 preprocessor. ',moddir
  WINTITLE=TRIM(WINTITLE)
  OPEN  (WIN,FILE='user',TITLE=WINTITLE)
  WINFO%TYPE = QWIN$MAX
  RESULT     = SETWSIZEQQ   (WIN,WINFO)
  STATUS     = FLOODFILL    (X,Y,ONE)
  I          = SETBKCOLOR   (15)
  I          = SETTEXTCOLOR (ONE)
  I          = SETEXITQQ    (QWIN$EXITPERSIST)
!  I          = QWIN$SCROLLDOWN    ! Check this

! Variable initializations

  FMTA       = '(A)'
  DEG        =  CHAR(248)//'C'
  ESC        =  CHAR(027)
  UNIT       =  CON
  DELETE_ERR = .TRUE.
  DELETE_WRN = .TRUE.
  ICST=0
  EUHF=0;TUHF=0;CUHF=0;EDHF=0;TDHF=0;CDHF=0;XX=0.0
  sediment_diagenesis=.false.    ! initialization SW 7/24/2017

! Initialize input/output files

  OPEN (CON,FILE=CONFN,STATUS='OLD',IOSTAT=IERR)
    IF (IERR /= 0) THEN
    CLOSE(CON)
    WRITE (WIN,*)'Could not open w2_con.npt'
    CONFN='w2_con.csv'
    OPEN (CON,FILE=CONFN,STATUS='OLD',IOSTAT=I)
       IF (I /= 0) THEN
       WRITE (WIN,*)'Could not open w2_con.npt or w2_con.csv'
       STOP
       ENDIF
       WRITE (WIN,*) 'Found w2_con.csv'
  END IF
!IF (IERR /= 0) THEN
!    WRITE (WIN,*) 'Could not open control file [CONFN="w2_con.npt"]'
!    STOP
!  END IF
  OPEN (ERR,FILE=ERRFN,STATUS='UNKNOWN',IOSTAT=IERR)
  IF (IERR /= 0) THEN
    WRITE (WIN,*) 'Could not open preprocessor error output file [ERRFN="pre.err"]'
    STOP
  END IF
  OPEN (WRN,FILE=WRNFN,STATUS='UNKNOWN',IOSTAT=IERR)
  IF (IERR /= 0) THEN
    WRITE (WIN,*) 'Could not open preprocessor warning output file [WRNFN="pre.wrn"]'
    STOP
  END IF
  OPEN (INI,FILE=INIFN,STATUS='UNKNOWN',IOSTAT=IERR)
  IF (IERR /= 0) THEN
    WRITE (WIN,*) 'Could not open preprocessor output file [INIFN="pre.opt"]'
    STOP
  END IF

!***********************************************************************************************************************************
!*                                                     Task 1:  Read all inputs                                                   **
!***********************************************************************************************************************************

! Title cards

  WRITE (WIN,*) 'Control file'
  ALLOCATE (TITLE(11))
  WRITE (WIN,*) '  title cards'
  
  
IF(CONFN=='w2_con.npt')THEN

  READ (CON,'(//A8/(8X,A72))',ERR=400) AID, (TITLE(J),J=2,11)
  IF (AID /= 'TITLE C ')               GO TO 400

! Array dimensions

  WRITE (WIN,*) '  array dimensions'
  READ (CON,'(/A8/(8X,5I8,A8,I8))')       AID, NWB, NBR, IMX, KMX, NPROC, CLOSEC,NPOLY
  
  IF(NPOLY==1)DOWNSTREAMDAM = .FALSE.   ! This implies polygon coordinates are starting at US(JB=1) rather than at the dam - but currently this is disabled
  
  
  IF (AID /= 'GRID    ')               GO TO 400
  READ (CON,'(/A8/(8X,10I8))')         AID, NTR, NST, NIW, NWD,  NGT, NSP, NPI, NPU
  IF (AID /= 'IN/OUTFL')               GO TO 400
  READ (CON,'(/A8/(8X,10I8))')         AID, NGC, NSS, NAL, NEP,  NBOD, NMC, NZP
  IF (AID /= 'CONSTITU')               GO TO 400
  READ (CON,'(/A8/(8X,I8,8A8))')       AID, NOD,SELECTC,HABTATC,ENVIRPC,AERATEC,INITUWL,ORGCC,SED_DIAG
  IF (AID /= 'MISCELL ')               GO TO 400
  
ELSE
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)
  DO J=2,11
  READ (CON,*,ERR=400)TITLE(J)
  ENDDO
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400) AID
  IF (AID /= 'NWB     ') THEN
      AID1='NWB     ';ICST=3
      GO TO 400
  ENDIF
  
  READ (CON,*,ERR=400) NWB, NBR, IMX, KMX, NPROC, CLOSEC, NPOLY; CLOSEC=ADJUSTR(CLOSEC)                   !'(A,5I0,A)'  
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400) AID
  IF (AID /= 'NTR     ') THEN
      AID1='NTR     ';ICST=3
      GO TO 400
  ENDIF
    IF(NPOLY==1)DOWNSTREAMDAM = .FALSE.
  
  READ (CON,*,ERR=400) NTR, NST, NIW, NWD, NGT, NSP, NPI, NPU
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)  AID
  IF (AID /= 'NGC     ') THEN
    AID1='NGC     ';ICST=3;GO TO 400
  ENDIF
  
  READ (CON,*,ERR=400) NGC, NSS, NAL, NEP, NBOD, NMC, NZP  
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400) AID
  IF (AID /= 'NDAY    ') THEN
    AID1='NDAY    ';ICST=3;GO TO 400  
  ENDIF
  
  READ (CON,*,ERR=400) NOD,SELECTC,HABTATC,ENVIRPC,AERATEC,INITUWL,ORGCC,SED_DIAG           !'(I0,5(A))'   
  SELECTC=ADJUSTR(SELECTC);HABTATC=ADJUSTR(HABTATC);ENVIRPC=ADJUSTR(ENVIRPC);AERATEC=ADJUSTR(AERATEC);INITUWL=ADJUSTR(INITUWL)    
  SED_DIAG=ADJUSTR(SED_DIAG)
ENDIF

  IF(NTR==0)THEN
      NTR1=1
  ELSE
      NTR1=NTR
  ENDIF
! Constituent numbers
  CCONSTIT=''
  NTDS  = 1
  CCONSTIT(NTDS)='TDS or SALINITY'
  NGCS  = 2
  NGCE  = NGCS+NGC-1
  JJ=0
  DO J=NGCS,NGCE
  JJ=JJ+1
  WRITE(IIN,'(I2)')JJ
  CCONSTIT(J)='GENERIC#'//IIN
  ENDDO
  NSSS  = NGCE+1
  NSSE  = NSSS+NSS-1
  JJ=0
  DO J=NSSS,NSSE
  JJ=JJ+1
  WRITE(IIN,'(I2)')JJ
  CCONSTIT(J)='SS#'//IIN
  ENDDO
  NWAGE = NSSE+1;CCONSTIT(NWAGE)='WATERAGE'
  NBACT = NWAGE+1;CCONSTIT(NBACT)='BACT'
  NDGP  = NBACT+1;CCONSTIT(NDGP)='DGP'
  NN2   = NDGP+1;CCONSTIT(NN2)='N2'
  NH2S  = NN2+1;CCONSTIT(NH2S)='H2S'
  NCH4  = NH2S+1;CCONSTIT(NCH4)='CH4'
  NSO4  = NCH4+1;CCONSTIT(NSO4)='SO4'
  NFEII = NSO4+1;  CCONSTIT(NFEII)='FEII'
  NFEOOH= NFEII+1;  CCONSTIT(NFEOOH)='FEOOH'
  NMNII = NFEOOH+1;  CCONSTIT(NMNII)='MNII'
  NMNO2 = NMNII+1;  CCONSTIT(NMNO2)='MNO2'
  
  NPO4  = NMNO2+1
  CCONSTIT(NPO4)='PO4-P'
  NNH4  = NPO4+1
  CCONSTIT(NNH4)='NNH4-N'
  NNO3  = NNH4+1
  CCONSTIT(NNO3)='NO3-N'
  NDSI  = NNO3+1
  CCONSTIT(NDSI)='DSI-SI'
  NPSI  = NDSI+1
  CCONSTIT(NPSI)='PSI-SI'
  NLDOM = NPSI+1
  CCONSTIT(NLDOM)='LDOM'
  NRDOM = NLDOM+1
  CCONSTIT(NRDOM)='RDOM'
  NLPOM = NRDOM+1
  CCONSTIT(NLPOM)='LPOM'
  NRPOM = NLPOM+1
  CCONSTIT(NRPOM)='RPOM'
  NBODS = NRPOM+1
 ! NBODE = NBODS+NBOD-1
if(nbod.gt.0)then
 ALLOCATE (nbodc(nbod), nbodp(nbod), nbodn(nbod))
  ibod=nbods
  nbodcs=ibod
  do jcb=1,nbod
     nbodc(jcb)=ibod
     ibod=ibod+1
  end do
  nbodce=ibod-1
  JJ=0
  DO J=NBODS,NBODCE
  JJ=JJ+1
  WRITE(IIN,'(I2)')JJ
  CCONSTIT(J)='CBOD#'//IIN
  ENDDO
  nbodps=ibod
  do jcb=1,nbod
    nbodp(jcb)=ibod
    ibod=ibod+1
  end do
  nbodpe=ibod-1
   JJ=0
  DO J=NBODPS,NBODPE
  JJ=JJ+1
  WRITE(IIN,'(I2)')JJ
  CCONSTIT(J)='CBOD-P#'//IIN
  ENDDO
  nbodns=ibod
  do jcb=1,nbod
    nbodn(jcb)=ibod
    ibod=ibod+1
  end do
  nbodne=ibod-1
  JJ=0
  DO J=NBODNS,NBODNE
  JJ=JJ+1
  WRITE(IIN,'(I2)')JJ
  CCONSTIT(J)='CBOD-N#'//IIN
  ENDDO
  else
    nbodcs=1; nbodce=1; nbodps=1; nbodpe=1; nbodns=1; nbodne=1     ! SW 11/21/2016
  
end if
!    NBODE = NBODS+NBOD-1
  NBODE = NBODS+NBOD*3-1
  NAS   = NBODE+1
  NAE   = NAS+NAL-1
   JJ=0
  DO J=NAS,NAE
  JJ=JJ+1
  WRITE(IIN,'(I2)')JJ
  CCONSTIT(J)='ALGAE#'//IIN
  ENDDO
  NDO   = NAE+1
  CCONSTIT(NDO)='DO'
  NTIC  = NDO+1
  CCONSTIT(NTIC)='TIC-C'
  NALK  = NTIC+1
  CCONSTIT(NALK)='ALK as CaCO3'
  NZOOS = NALK+1
  NZOOE = NZOOS+NZP-1
   JJ=0
  DO J=NZOOS,NZOOE
  JJ=JJ+1
  WRITE(IIN,'(I2)')JJ
  CCONSTIT(J)='ZOOP#'//IIN
  ENDDO
  NLDOMP= NZOOE+1
  CCONSTIT(NLDOMP)='LDOM-P'
  NRDOMP= NLDOMP+1
  CCONSTIT(NRDOMP)='RDOM-P'
  NLPOMP= NRDOMP+1
  CCONSTIT(NLPOMP)='LPOM-P'
  NRPOMP= NLPOMP+1
  CCONSTIT(NRPOMP)='RPOM-P'
  NLDOMN= NRPOMP+1
  CCONSTIT(NLDOMN)='LDOM-N'
  NRDOMN= NLDOMN+1
  CCONSTIT(NRDOMN)='RDOM-N'
  NLPOMN= NRDOMN+1
  CCONSTIT(NLPOMN)='LPOM-N'
  NRPOMN= NLPOMN+1
  CCONSTIT(NRPOMN)='RPOM-N'
  NATS=NRPOMN+1
  NATE=NATS+3
  CCONSTIT(NATS)='MICROCYSTIN'
  CCONSTIT(NATS+1)='CYLINDROSPERMOPSIN'
  CCONSTIT(NATS+2)='ANATOXIN-A'
  CCONSTIT(NATS+3)='SAXITOXIN'
  
  NCT=NATE    !NRPOMN

! Constituent, tributary, and widthdrawal totals

  !NCT  = NRPOMN
  NTRT = NTR+NGT+NSP+NPI+NPU
  NWDT = NWD+NGT+NSP+NPI+NPU
  NEPT = MAX(NEP,1)     ! ,1
  NEPTT= MAX(NEP,5)
  Nmct = MAX(nmc,1)     ! ,1
  NZPt  = MAX(nzp,1)    ! ,1    ! SW 6/1/07
  NALT =MAX(NAL,5)
  NZPTT = MAX(NZPT,5)
  NMCTT=MAX(NMCT,5)
  IF(NST>5)THEN
      NSTT=NST    ! FIXED # OF ROWS IN W2_CON.CSV MIN IS 5 INCREASE IF NECESSARY
  ELSE
      NSTT=5
  ENDIF
  
! Allocation declarations
  ALLOCATE (ATMDCN(NCT,NWB), NACATD(NWB),ATM_DEP_AVG(NCT,NWB),ATM_DEP_MIN(NCT,NWB),ATM_DEP_MAX(NCT,NWB))
  ALLOCATE (C_ATM_DEPOSITION(NCT,NWB), ATM_DEPOSITIONC(NWB), ATM_DEP_LOADING(NCT,NWB),ATM_DEPOSITION_INTERPOLATION(NWB),ATMDEPFN(NWB),ATMDEP_P(NWB),ATMDEP_N(NWB))
  ALLOCATE (EXC(NWB),    EXIC(NWB))
  ALLOCATE (SSAC(NSS),   ALAC(NAL))
  ALLOCATE (CUNIT(NCT),  CNAME(NCT),  CNAME1(NCT), CNAME2(NCT), CDNAME(NDC), NADC(NWB), CDMULT(NDC), CMULT(NCT), CUNIT1(NCT), CUNIT2(NCT))    ! SW 1/16/04 8/25/05
  ALLOCATE (DMULT(NDC))
  ALLOCATE (HNAME(NHY),  HFMT(NHY))
  ALLOCATE (SLTRC(NWB),  THETA(NWB),  FRICC(NWB))
  ALLOCATE (WINDC(NWB),  QINC(NWB),   QOUTC(NWB),  HEATC(NWB),  SLHTC(NWB))
  ALLOCATE (QINIC(NBR),  DTRIC(NBR),  TRIC(NTR),   WDIC(NWD),   HDIC(NBR),   METIC(NWB))
  ALLOCATE (VBC(NWB),    EBC(NWB),    MBC(NWB),    PQC(NWB),    EVC(NWB),    PRC(NWB))
  ALLOCATE (ICEC(NWB),   SLICEC(NWB), ICEI(NWB),   ALBEDO(NWB), HWI(NWB),    BETAI(NWB),  GAMMAI(NWB), ICEMIN(NWB), ICET2(NWB))
  ALLOCATE (EXH2O(NWB),  BETA(NWB),   EXOM(NWB),   EXSS(NWB),   SROC(NWB), Z0(NWB))
  ALLOCATE (CBHE(NWB),   TSED(NWB),   FI(NWB),     TSEDF(NWB),  AFW(NWB),    BFW(NWB),    CFW(NWB),    WINDH(NWB),  RHEVC(NWB))
  ALLOCATE (SEDC(NWB),   PRNSC(NWB),  SEDCI(NWB),  SEDDK(NWB),  FSOD(NWB),   FSED(NWB), seds(nwb), sedbr(nwb),DYNSEDK(NWB))   ! SW 6/1/07
  ALLOCATE (SODT1(NWB),  SODT2(NWB),  SODK1(NWB),  SODK2(NWB))
  ALLOCATE (FETCHC(NWB), SDC(NWB),    DX(NWB))
  ALLOCATE (SNPC(NWB),   PRFC(NWB),   SPRC(NWB),   CPLC(NWB),   VPLC(NWB),   FLXC(NWB),   SCRC(NWB))
  ALLOCATE (NSNP(NWB),   NPRF(NWB),   NSPR(NWB),   NCPL(NWB),   NVPL(NWB),   NFLX(NWB),   NSCR(NWB), TCPL(NWB))
  ALLOCATE (NISNP(NWB),  NIPRF(NWB),  NISPR(NWB))
  ALLOCATE (SNPFN(NWB),  PRFFN(NWB),  SPRFN(NWB),  VPLFN(NWB),  CPLFN(NWB),  FLXFN(NWB))
  ALLOCATE (BTHFN(NWB),  METFN(NWB),  VPRFN(NWB),  LPRFN(NWB),  JJS(NST))
  ALLOCATE (QINFN(NBR),  TINFN(NBR),  CINFN(NBR),  CDHFN(NBR),  QDTFN(NBR),  TDTFN(NBR),  CDTFN(NBR),  PREFN(NBR),  TPRFN(NBR))
  ALLOCATE (CPRFN(NBR),  EUHFN(NBR),  TUHFN(NBR),  CUHFN(NBR),  EDHFN(NBR),  TDHFN(NBR),  QOTFN(NBR),TINMAX(NBR), TINMIN(NBR),TRTMAX(NTR),TRTMIN(NTR))   ! SW 3/28/13
  ALLOCATE (TDHMIN(NBR), TDHMAX(NBR), TUHMIN(NBR), TUHMAX(NBR))
  ALLOCATE (XF(IMX))
  ALLOCATE (TUH(KMX),TDH(KMX),PUMP_DOWNSTREAM(NPU))  ! SW 3/1/2020
  ALLOCATE (EXTFN(NWB))                                                                                                !TC 12/12/01
  ALLOCATE (BODP(NBOD),  BODN(NBOD),  BODC(NBOD))                                                                      !TC 01/15/02
  ALLOCATE (DTRC(NBR),   ALPHA(NBR))
  ALLOCATE (AX(NWB),     WTYPEC(NWB), GRIDCC(NWB),JBDN(NWB),   AZSLC(NWB),  AZMAX(NWB),  KBMAX(NWB),  AZC(NWB),FBC(NWB),AZE(NWB),ARODI(NWB),STRCKLR(NWB),BOUNDFR(NWB),TKECAL(NWB))
  ALLOCATE (VISC(NWB),   CELC(NWB), DLTADD(NWB))
  ALLOCATE (BTH(NWB),    VPR(NWB),    LPR(NWB))
  ALLOCATE (T2I(NWB),    KTWB(NWB))                                                                                    !TC 06/11/02
  ALLOCATE (QINAVW(NWB), QOTAVW(NWB), QINMXW(NWB), QOTMXW(NWB))
  ALLOCATE (LAT(NWB),    LONG(NWB),   ELBOT(NWB),  BS(NWB),     BE(NWB), YNORTHING(NWB), XEASTING(NWB), UTMZONE(NWB))  ! SW 4/6/2015
  ALLOCATE (XC(IMX),YC(IMX),Xnw(IMX),Ynw(IMX),xw(IMX),yw(IMX),Xsw(IMX),Ysw(IMX),Xse(IMX),Yse(IMX),Xne(IMX),Yne(IMX),xe(IMX),ye(IMX),XCN(IMX),YCN(IMX),XCS(IMX),YCS(IMX))
  ALLOCATE (IEPR(NWB),   KEPR(NWB))
  ALLOCATE (REAERC(NWB), NEQN(NWB))
  ALLOCATE (COLQ10(NWB), COLDK(NWB))
  ALLOCATE (LDOMDK(NWB), RDOMDK(NWB), LRDDK(NWB))
  ALLOCATE (LPOMDK(NWB), RPOMDK(NWB), LRPDK(NWB),  POMS(NWB),   APOM(NAL),   ANPR(NAL),   ANEQN(NAL))
  ALLOCATE (ORGP(NWB),   ORGN(NWB),   ORGC(NWB),   ORGSI(NWB))
  ALLOCATE (OMT1(NWB),   OMT2(NWB),   OMK1(NWB),   OMK2(NWB))
  ALLOCATE (KBOD(NBOD),  TBOD(NBOD),  RBOD(NBOD))                                                                      !TC 09/23/02
  ALLOCATE (PO4R(NWB),   PARTP(NWB))
  ALLOCATE (NH4R(NWB),   NH4DK(NWB))
  ALLOCATE (NH4T1(NWB),  NH4T2(NWB),  NH4K1(NWB),  NH4K2(NWB))
  ALLOCATE (NO3DK(NWB),  NO3S(NWB),FNO3SED(NWB))
  ALLOCATE (NO3T1(NWB),  NO3T2(NWB),  NO3K1(NWB),  NO3K2(NWB))
  ALLOCATE (DSIR(NWB),   PSIS(NWB),   PSIDK(NWB),  PARTSI(NWB))
  ALLOCATE (FER(NWB),    FES(NWB))
  ALLOCATE (CO2R(NWB),   DYNSTRUC(NBR))
  ALLOCATE (O2NH4(NWB),  O2OM(NWB),   O2AR(NAL),   O2AG(NAL))                                                          !TC 09/23/02
  ALLOCATE (CGQ10(NGC),  CG0DK(NGC),  CG1DK(NGC),  CGS(NGC),CGLDK(NGC),CGKLF(NGC),CGCS(NGC))
  ALLOCATE (sromax(nwb),sromin(nwb),tairmax(nwb),tairmin(nwb),tdewmax(nwb),tdewmin(nwb),windmax(nwb),windmin(nwb),phimin(nwb),phimax(nwb),cloudmax(nwb),cloudmin(nwb),tdewav(nwb),sroav(nwb),tairav(nwb),phiav(nwb),windav(nwb),cloudav(nwb))
  ALLOCATE (CAC(NCT))
  ALLOCATE (RCOEF1(NWB), RCOEF2(NWB), RCOEF3(NWB), RCOEF4(NWB),DGPO2(NWB),MINKL(NWB))
  ALLOCATE (CIN(NCT),    CIN2(NCT), CPR(NCT))                                                                                    !SW 01/07/01  4/11/2018
  ALLOCATE (TTLB(IMX),   TTRB(IMX),   CLLB(IMX),   CLRB(IMX),   SRLB1(IMX),  SRRB1(IMX),   SRLB2(IMX),  SRRB2(IMX))
  ALLOCATE (SRFJDAY1(IMX),SHADEIN(IMX),SRFJDAY2(IMX),TOPO(IMX,IANGLES),ANG(IANGLES))                                   !SW 04/03/02
  ALLOCATE (VOLB(NBR),   VOLG(NWB))
  ALLOCATE (XBR(NBR),    XGR(NWB),    XI(IMX))
  ALLOCATE (US(NBR),     DS(NBR),     CUS(NBR),    UHS(NBR),    DHS(NBR),    UQB(NBR),     DQB(NBR))
  ALLOCATE (NL(NBR),     NPOINT(NBR), SLOPE(NBR),SLOPEC(NBR))
  ALLOCATE (NWBR(NBR),   NTBR(NBR),   NSTR(NBR),   NUPIBR(NBR), NDPIBR(NBR), NUGTBR(NBR),  NDGTBR(NBR))
  ALLOCATE (QTRAVB(NBR), QTRMXB(NBR))
  ALLOCATE (QINAV(NBR),  QPRAV(NBR),  QDTAV(NBR),  QSTRAV(NBR))
  ALLOCATE (QINS(NBR),   QPRS(NBR),   QDTS(NBR),   QSTRS(NBR))
  ALLOCATE (QINMX(NBR),  QDTMX(NBR),  QSTRMX(NBR), QPRMX(NBR))
  ALLOCATE (QINMXB(NBR), QOTMXB(NBR))
  ALLOCATE (QINAVB(NBR), QOTAVB(NBR))
  ALLOCATE (QWDAVB(NBR), QWDMXB(NBR))
  ALLOCATE (QDNAVB(NBR))
  ALLOCATE (NQIN(NBR),   NQDT(NBR),   NQPR(NBR),   NQOT(NBR))
  ALLOCATE (Z(IMX),      DLX(IMX),    KB(IMX),     BK(IMX))
  ALLOCATE (ICETH(IMX),  KTI(IMX),    ELWS(IMX),   SHADE(IMX))
  ALLOCATE (SOD(IMX),    PHI0(IMX),   FRIC(IMX))
  ALLOCATE (CDN(NDC,NWB))
  ALLOCATE (CN(NCT),     INCN(NCT,NBR),   TRCN(NCT,NTR), incdt(nct,nbr))                                                               !SW 01/07/01
  ALLOCATE (SSS(NSS),    SEDRC(NSS),   TAUCR(NSS))                                                                      !TC 08/20/03 SW 1/16/04
  ALLOCATE (AG(NAL),     AR(NAL),     AE(NAL),     AM(NAL),     AS(NAL),     EXA(NAL),    ASAT(NAL),   ALGP(NAL),   ALGN(NAL))
  ALLOCATE (ALGC(NAL),   ALGSI(NAL),  ACHLA(NAL),  AHSP(NAL),   AHSN(NAL),   AHSSI(NAL))
  ALLOCATE (AT1(NAL),    AT2(NAL),    AT3(NAL),    AT4(NAL),    AK1(NAL),    AK2(NAL),    AK3(NAL),    AK4(NAL), AVERTM(NAL))
  ALLOCATE (EG(nept),     ER(nept),     EE(nept),     EM(nept),     EB(nept),     ESAT(nept),   EP(nept),     EN(nept))
  ALLOCATE (EC(nept),     ESI(nept),    ECHLA(nept),  EHSP(nept),   EHSN(nept),   EHSSI(nept),  EPOM(nept),   EHS(nept))
  ALLOCATE (ET1(nept),    ET2(nept),    ET3(nept),    ET4(nept),    EK1(nept),    EK2(nept),    EK3(nept),    EK4(nept))
  ALLOCATE (ENPR(nept),   ENEQN(nept),  O2ER(nept),   O2EG(nept),KG_H2O_CONSTANT(NWB))
  ALLOCATE (TSRKTD(NOD), TSRKTF(NOD), IWDO(NOD),   RSOD(NOD),   RSOF(NOD),   DLTD(NOD),   DLTF(NOD),   DLTMAX(NOD))
  ALLOCATE (IUPI(NPI),   IDPI(NPI),   EUPI(NPI),   EDPI(NPI),   WPI(NPI),    DLXPI(NPI),  FPI(NPI),    FMINPI(NPI), PUPIC(NPI), DYNPIPE(NPI))
  ALLOCATE (ETUPI(NPI),  EBUPI(NPI),  KTUPI(NPI),  KBUPI(NPI),  PDPIC(NPI),  ETDPI(NPI),  EBDPI(NPI),  KTDPI(NPI))
  ALLOCATE (KBDPI(NPI),  JBUPI(NPI),  JBDPI(NPI))
  ALLOCATE (PUSPC(NSP),  ETUSP(NSP),  EBUSP(NSP),  KTUSP(NSP),  KBUSP(NSP),  PDSPC(NSP), ETDSP(NSP), EBDSP(NSP))
  ALLOCATE (KTDSP(NSP),  KBDSP(NSP),  IUSP(NSP),   IDSP(NSP),   ESP(NSP),    A1SP(NSP),  B1SP(NSP),  A2SP(NSP))
  ALLOCATE (B2SP(NSP),   AGASSP(NSP), BGASSP(NSP), CGASSP(NSP), EQSP(NSP),   GASSPC(NSP))
  ALLOCATE (GTA1(NGT),   GTB1(NGT),   GTA2(NGT),   GTB2(NGT),   IUGT(NGT),   IDGT(NGT),   EGT(NGT), EGT2(NGT),  A1GT(NGT),   B1GT(NGT))
  ALLOCATE (G1GT(NGT),   A2GT(NGT),   B2GT(NGT),   G2GT(NGT),   BGATE(NGT),  PUGTC(NGT),  ETUGT(NGT),  EBUGT(NGT),  KTUGT(NGT))
  ALLOCATE (KBUGT(NGT),  PDGTC(NGT),  ETDGT(NGT),  EBDGT(NGT),  KTDGT(NGT),  KBDGT(NGT),  AGASGT(NGT), BGASGT(NGT), CGASGT(NGT))
  ALLOCATE (GASGTC(NGT), EQGT(NGT),   JBUGT(NGT),  JBDGT(NGT))
  ALLOCATE (LATGTC(NGT), LATSPC(NSP), LATPIC(NPI), LATPUC(NPU), DYNGTC(NGT),GTIC(NGT))
  ALLOCATE (IUPU(NPU),   IDPU(NPU),   EPU(NPU),    STRTPU(NPU), ENDPU(NPU),  EONPU(NPU),  EOFFPU(NPU), QPU(NPU),    PPUC(NPU))
  ALLOCATE (ETPU(NPU),   EBPU(NPU),   KTPU(NPU),   KBPU(NPU))
  ALLOCATE (IWR(NIW),    KTWR(NIW),   KBWR(NIW),EKTWR(NIW), EKBWR(NIW))
  ALLOCATE (IWD(NWDT),   JBWD(NWD),   EWD(NWDT),   KTWD(NWDT),  KBWD(NWDT))
  ALLOCATE (ITR(NTRT),   QTRFN(NTR),  TTRFN(NTR),  CTRFN(NTR),  ETTR(NTRT),  EBTR(NTRT),  PTRC(NTRT),IBR(NBR))
  ALLOCATE (JBTR(NTR),   NQTR(NTR),   NACTR(NTR),  NACIN(NBR),  NACPR(NWB),  NACDT(NBR))
  ALLOCATE (TSRD(NOD),   TSRF(NOD),   WDOD(NOD),   WDOF(NOD))
  ALLOCATE (T2(KMX,IMX),DYNPUMP(NPU))
  ALLOCATE (ITSR(IMX*KMX),   ETSR(IMX*KMX))
  ALLOCATE (ESTR(NSTT,NBR),   WSTR(NSTT,NBR),   SINKC(NSTT,NBR),  QSTR(NST),       QSTRO(NST),      STRIC(NSTT,NBR))
  ALLOCATE (KTSTR(NSTT,NBR),  KBSTR(NSTT,NBR))
  ALLOCATE (TDS(KMX,IMX),    COL(KMX,IMX),    NH4(KMX,IMX),    NO3(KMX,IMX),    PO4(KMX,IMX),    FE(KMX,IMX),     DSI(KMX,IMX))
  ALLOCATE (PSI(KMX,IMX),    LDOM(KMX,IMX),   RDOM(KMX,IMX),   LPOM(KMX,IMX),   RPOM(KMX,IMX),   O2(KMX,IMX),     TIC(KMX,IMX))
  ALLOCATE (ALK(KMX,IMX))
  ALLOCATE (TN(KMX,IMX),     TP(KMX,IMX))
  ALLOCATE (DON(KMX,IMX),    DOP(KMX,IMX),    DOC(KMX,IMX))
  ALLOCATE (PON(KMX,IMX),    POP(KMX,IMX),    POC(KMX,IMX))
  ALLOCATE (TON(KMX,IMX),    TOP(KMX,IMX),    TOC(KMX,IMX))
  ALLOCATE (CHLA(KMX,IMX),   ATOT(KMX,IMX))
  ALLOCATE (O2DG(KMX,IMX),   TOTSS(KMX,IMX),  TISS(KMX,IMX))
  ALLOCATE (CBODU(KMX,IMX),  TKN(KMX,IMX))
  ALLOCATE (CINBRC(NCT,NBR), CINTRC(NCT,NTR1), CDTBRC(NCT,NBR), CPRBRC(NCT,NBR), CPRWBC(NCT,NWB), CFWBC(NFL,NWB),  HPRWBC(NHY,NWB))
  ALLOCATE (QTRS(NTR,NBR),   QTRAV(NTR,NBR),  QTRMXS(NTR,NBR), QTRMX(NTR,NBR),  CTRAVG(NCT,NTR), CTRMIN(NCT,NTR), CTRMAX(NCT,NTR))
  ALLOCATE (WSC(IMX),    SNPD(NOD,NWB),   SNPF(NOD,NWB),   SCRD(NOD,NWB),   SCRF(NOD,NWB),   PRFD(NOD,NWB))
  ALLOCATE (QWDS(NWD,NBR),   QWDAV(NWD,NBR),  QWDMX(NWD,NBR))
  ALLOCATE (PRFF(NOD,NWB),   SPRD(NOD,NWB),   SPRF(NOD,NWB),   CPLD(NOD,NWB),   CPLF(NOD,NWB))
  ALLOCATE (VPLD(NOD,NWB),   VPLF(NOD,NWB),   FLXD(NOD,NWB),   FLXF(NOD,NWB))
  ALLOCATE (NCCBR(KMX,NBR),  CVBR(KMX,NBR),   SABR(KMX,NBR),   HBR(KMX,NBR),    BBR(KMX,NBR))
  ALLOCATE (SAGR(KMX,NWB),   NCCGR(KMX,NWB),  CVGR(KMX,NWB),   HGR(KMX,NWB),    BGR(KMX,NWB))
  ALLOCATE (TVP(KMX,NWB),    H(KMX,NWB))
  ALLOCATE (B(KMX,IMX),      CONV(KMX,IMX),   EL(KMX,IMX))
  ALLOCATE (CDWBC(NDC,NWB))
  ALLOCATE (C2IWB(NCT,NWB),  CINAVG(NCT,NBR), CINMIN(NCT,NBR), CINMAX(NCT,NBR),CDTAVG(NCT,NBR), CDTMIN(NCT,NBR), CDTMAX(NCT,NBR),cinload(nct,nbr),ctrload(nct,ntr),cdtload(nct,nbr),cinnload(nbr),cinpload(nbr),ctrnload(ntr),ctrpload(ntr))
  ALLOCATE (EPIC(NWB,NEPTT), EPICI(NWB,NEPTT),EPIPRC(NWB,NEPTT))                                                       ! cb 6/16/06
  ALLOCATE (IPR(IMX,NWB),    ISNP(IMX,NWB),   IPRF(IMX,NWB),   ISPR(IMX,NWB))
  ALLOCATE (CVP(KMX,NCT,NWB), EPIVP(KMX,IMX,NEPT),EPD(KMX,IMX,NEPT))
  ALLOCATE (C2(KMX,IMX,NCT), CD(KMX,IMX,NDC), CG(KMX,IMX,NGC), ALG(KMX,IMX,NAL), SS(KMX,IMX,NSS), CBOD(KMX,IMX,NBOD))
  ALLOCATE (SHIFT_DEMAND(NWB), SEDIMENT_CALC(NWB))
  ALLOCATE (OPEN_VPR(NWB),     OPEN_LPR(NWB),    OPEN_LPRC(NWB),  ISO_TEMP(NWB),    VERT_TEMP(NWB))
  ALLOCATE (UH_EXTERNAL(NBR),  DH_EXTERNAL(NBR),   UH_INTERNAL(NBR), DH_INTERNAL(NBR), DAM_FLOW(NBR))
  ALLOCATE (UQ_EXTERNAL(NBR),  DQ_EXTERNAL(NBR),   UQ_INTERNAL(NBR), DQ_INTERNAL(NBR))
  ALLOCATE (ISO_CONC(NCT,NWB), VERT_CONC(NCT,NWB), LONG_TEMP(NWB),   LONG_SED(NWB),    LONG_CONC(NCT,NWB))
  ALLOCATE (ICE_CALC(NWB),     DIST_TRIB(NBR),     PRECIPITATION(NWB))
  ALLOCATE (ZERO_SLOPE(NWB))                                                                                           !SW 10/17/02
  ALLOCATE (INTERP_INFLOW(NBR))                                                                                       ! cb 9/7/2010
  ALLOCATE (CMIN(NCT),   CMAX(NCT),   HYMIN(NHY),  HYMAX(NHY),  CDMIN(NDC),  CDMAX(NDC))                               ! SW 1/16/04
  ALLOCATE (FMTH(NHY),   HMULT(NHY),  FMTC(NCT),   FMTCD(NDC))                                                         ! SW 1/16/04
  ALLOCATE (CPLTC(NCT),  HPLTC(NHY),  CDPLTC(NDC))                                                                     ! SW 1/16/04
  ALLOCATE (zg(NZPt),zm(NZPt),zeff(NZPt),prefp(NZPt),zr(NZPt),zoomin(NZPt),zs2p(NZPt),ZS(NZPT),exz(NZPt),PREFZ(NZPtT,nzpt),o2zr(nzpt))
  ALLOCATE (zt1(NZPt),zt2(NZPt),zt3(NZPt),zt4(NZPt),zk1(NZPt),zk2(NZPt),zk3(NZPt),zk4(NZPt))
  ALLOCATE (zp(NZPt),zn(NZPt),zc(NZPt))
  allocate (prefa(nalT,nzpt))
  allocate (po4zr(kmx,imx),nh4zr(kmx,imx))
  allocate (zmu(kmx,imx,nzpt),tgraze(kmx,imx,nzpt),zrt(kmx,imx,nzpt),zmt(kmx,imx,nzpt)) ! MLM POINTERS:,zoo(kmx,imx,nzpt),zooss(kmx,imx,nzpt))
  allocate (zoorm(kmx,imx,nzpt),zoormr(kmx,imx,nzpt),zoormf(kmx,imx,nzpt))
  allocate (lpzooout(kmx,imx),lpzooin(kmx,imx),dozr(kmx,imx),ticzr(kmx,imx))
  allocate (agz(kmx,imx,nal,nzpt), zgz(kmx,imx,nzpt,nzpt),agzt(kmx,imx,nal)) !omnivorous zooplankton
  allocate (ORGPLD(kmx,imx), ORGPRD(kmx,imx), ORGPLP(kmx,imx), ORGPRP(kmx,imx), ORGNLD(kmx,imx), ORGNRD(kmx,imx), ORGNLP(kmx,imx))
  allocate (ORGNRP(kmx,imx))
  allocate (ldompmp(kmx,imx),ldomnmp(kmx,imx),lpompmp(kmx,imx),lpomnmp(kmx,imx),rpompmp(kmx,imx),rpomnmp(kmx,imx))
  allocate (lpzooinp(kmx,imx),lpzooinn(kmx,imx),lpzoooutp(kmx,imx),lpzoooutn(kmx,imx))
  allocate (SEDVPp(KMX,NWB),SEDVPc(KMX,NWB),SEDVPn(KMX,NWB))
  allocate (sedp(kmx,imx),sedn(kmx,imx),sedcar(kmx,imx),SED(KMX,IMX))
  allocate (sdkv(kmx,imx),seddktot(kmx,imx))
  allocate (sedcip(nwb),sedcin(nwb),sedcic(nwb),sedcis(nwb))
  ALLOCATE (cbods(NBOD), cbodns(kmx,imx), sedcb(kmx,imx), sedcbp(kmx,imx), sedcbn(kmx,imx), sedcbc(kmx,imx),nbodtot(kmx,imx),pbodtot(kmx,imx))
  allocate  (print_macrophyte(nwb,nmct), macrophyte_calc(nwb,nmct),macwbc(nwb,nmctT),conv2(kmx,kmx),mprwbc(nwb,nmctT))
  allocate  (mac(kmx,imx,nmct), macrc(kmx,kmx,imx,nmct),mact(kmx,kmx,imx), macrm(kmx,kmx,imx,nmct), macss(kmx,kmx,imx,nmct))
  allocate  (mgr(kmx,kmx,imx,nmct),mmr(kmx,imx,nmct), mrr(kmx,imx,nmct))
  allocate  (smacrc(kmx,kmx,imx,nmct), smacrm(kmx,kmx,imx,nmct))
  allocate  (smact(kmx,kmx,imx), smac(kmx,imx,nmct))
  allocate  (mt1(nmct),mt2(nmct),mt3(nmct),mt4(nmct),mk1(nmct),mk2(nmct),mk3(nmct),mk4(nmct),mg(nmct),mr(nmct),mm(nmct))
  allocate  (mbmp(nmct), mmax(nmct), cdstem(nmct),dmv(nmct),dwsa(nmct),anorm(nmct))
  allocate  (mp(nmct), mn(nmct), mc(nmct),psed(nmct),nsed(nmct),mhsp(nmct),mhsn(nmct),mhsc(nmct),msat(nmct),exm(nmct))
  allocate  (O2MG(nmct), O2MR(nmct),  LRPMAC(nmct),  MPOM(nmct))
  allocate  (kticol(imx),armac(imx),macwbci(nwb,nmctT))
  allocate  (macmbrs(nbr,nmct), macmbrt(nbr,nmct),ssmacmb(nbr,nmct))
  allocate  (cw(kmx,imx), bic(kmx,imx))
  allocate  (mactrmr(kmx,imx,nmct), mactrmf(kmx,imx,nmct),mactrm(kmx,imx,nmct), macrcvp(KMX,NWB,nmc),macrclp(KMX,imx,nmc))
  allocate  (mlfpr(kmx,kmx,imx,nmct))
  allocate  (mllim(kmx,kmx,imx,nmct), mplim(kmx,imx,nmct),mclim(kmx,imx,nmct),mnlim(kmx,imx,nmct))
  ALLOCATE  (GAMMAj(kmx,KMX,IMX))	
  allocate (por(kmx,imx),VOLKTi(imx),VOLi(Kmx,Imx),vstem(kmx,imx,nmct),vstemkt(imx,nmct),sarea(nmct))
  ALLOCATE (IWIND(NWB))
  allocate (cbodp(kmx,imx,nbod), cbodn(kmx,imx,nbod))
  ALLOCATE (ISO_EPIPHYTON(NWB,NEPt), VERT_EPIPHYTON(NWB,NEPt), LONG_EPIPHYTON(NWB,NEPt), EPIPHYTON_CALC(NWB,NEPt), ISO_MACRO(NWB,NMCT), VERT_MACRO(NWB,NMCT), LONG_MACRO(NWB,NMCT))         !TC 10/25/02
  ALLOCATE (ISO_SEDIMENT(NWB),    VERT_SEDIMENT(NWB),   LONG_SEDIMENT(NWB))
  ALLOCATE(SAGR1(KMX), CVGR1(KMX), NCCGR1(KMX), HGR1(KMX), BGR1(KMX))                                                      ! SW 11/1/13
  ALLOCATE (DDUM(NOD),FDUM(NOD), IDUM(IMX))
  ALLOCATE (CH4R(NWB), H2SR(NWB), FEIIR(NWB), MNIIR(NWB),SO4R(NWB))
  ALLOCATE (BACTQ10(NWB), BACT1DK(NWB), BACTLDK(NWB), BACTS(NWB))
  !ALLOCATE (A_DISG(NWB),B_DISG(NWB),C_DISG(NWB))
  ALLOCATE (CoeffA_Turb(NWB), CoeffB_Turb(NWB),SECC_PAR(NWB))
  ALLOCATE (CGR(NGC),SSCS(NSS))
  ALLOCATE (H2SQ10(NWB), H2S1DK(NWB), CH4Q10(NWB), CH41DK(NWB))
  ALLOCATE (KFE_OXID(NWB), KFE_RED(NWB), KFEOOH_HalfSat(NWB), FeSetVel(NWB))
  ALLOCATE (KMN_OXID(NWB), KMN_RED(NWB), KMNO2_HalfSat(NWB), MnSetVel(NWB))

  
! State variable pointers

  TDS  => C2(:,:,1);     CG   => C2(:,:,NGCS:NGCE); SS   => C2(:,:,NSSS:NSSE); PO4  => C2(:,:,NPO4);        NH4  => C2(:,:,NNH4)
  NO3  => C2(:,:,NNO3);  DSI  => C2(:,:,NDSI);      PSI  => C2(:,:,NPSI);      LDOM => C2(:,:,NLDOM)
  RDOM => C2(:,:,NRDOM); LPOM => C2(:,:,NLPOM);     RPOM => C2(:,:,NRPOM);  ALG  => C2(:,:,NAS:NAE)
  O2   => C2(:,:,NDO);   TIC  => C2(:,:,NTIC);      ALK  => C2(:,:,NALK);      zoo => C2(:,:,NZOOS:NZOOE)  ; LDOMP=> C2(:,:,NLDOMP)
  RDOMP=> C2(:,:,NRDOMP); LPOMP=> C2(:,:,NLPOMP); RPOMP=> C2(:,:,NRPOMP); LDOMN=> C2(:,:,NLDOMN);RDOMN=> C2(:,:,NRDOMN); LPOMN=> C2(:,:,NLPOMN); RPOMN=> C2(:,:,NRPOMN)
  CBOD => C2(:,:,nbodcs:nbodce); CBODp => C2(:,:,nbodps:nbodpe)  ; CBODn => C2(:,:,nbodns:nbodne)	

! Derived variable pointers

  !DOC  => CD(:,:,1);  POC  => CD(:,:,2); TOC => CD(:,:,3); DON => CD(:,:,4);  PON => CD(:,:,5);  TON  => CD(:,:,6)                             ! REVISED SW 2/18/16
  !TKN => CD(:,:,7); TN   => CD(:,:,8);  DOP  => CD(:,:,9); POP => CD(:,:,10); TOP => CD(:,:,11); TP  => CD(:,:,12); CHLA => CD(:,:,14)
  !ATOT => CD(:,:,15); O2DG => CD(:,:,16); TOTSS => CD(:,:,17); TISS => CD(:,:,18); CBODU => CD(:,:,19)
  DOC_DER=1;POC_DER=2;TOC_DER=3;DON_DER=4;PON_DER=5;TON_DER=6
  TKN_DER=7;TN_DER=8;NH3_DER=9;DOP_DER=10;POP_DER=11;TOP_DER=12
  TP_DER=13;APR_DER=14;CHLA_DER=15;ATOT_DER=16;O2DG_DER=17
  TDG_DER=18;TURB_DER=19; TOTSS_DER=20
  TISS_DER=21;CBODU_DER=22;PH_DER=23; CO2_DER=24; HCO3_DER=25; CO3_DER=26; SECCHI_DER=27

  DOC   => CD(:,:,DOC_DER);  POC  => CD(:,:,POC_DER);  TOC  => CD(:,:,TOC_DER);  DON  => CD(:,:,DON_DER)
  PON   => CD(:,:,PON_DER);  TON  => CD(:,:,TON_DER)
  TKN   => CD(:,:,TKN_DER);  TN   => CD(:,:,TN_DER);  NH3  => CD(:,:,NH3_DER)     
  DOP  => CD(:,:,DOP_DER);  POP  => CD(:,:,POP_DER); TOP   => CD(:,:,TOP_DER); TP   => CD(:,:,TP_DER)
  APR   => CD(:,:,APR_DER); CHLA => CD(:,:,CHLA_DER); ATOT => CD(:,:,ATOT_DER); O2DG => CD(:,:,O2DG_DER); 
  TDG   => CD(:,:,TDG_DER); TURB  => CD(:,:,TURB_DER)
  TOTSS => CD(:,:,TOTSS_DER); TISS => CD(:,:,TISS_DER)
  CBODU => CD(:,:,CBODU_DER); PH   => CD(:,:,PH_DER); CO2  => CD(:,:,CO2_DER)
  HCO3 => CD(:,:,HCO3_DER); CO3  => CD(:,:,CO3_DER); ; SECCHID => CD(:,:,SECCHI_DER)  


! Variable initializations

  DO JW=1,NWB
    KTWB(JW) = 2
    BTH(JW)  = NPT
    VPR(JW)  = NPT+1
    LPR(JW)  = NPT+2
    NPT      = NPT+3
  END DO

IF(CONFN=='w2_con.npt')THEN

! Time control cards

  WRITE (WIN,*) '  time control cards'
  READ (CON,'(/A8/8X,2F8.0,I8)',       ERR=400)  AID, TMSTRT, TMEND, YEAR
  IF (AID /= 'TIME CON')                         GO TO 400
  READ (CON,'(/A8/8X,I8,F8.0,A8)',        ERR=400)  AID, NDLT, DLTMIN,DLTINTR
  IF (AID /= 'DLT CON ')                          GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',       ERR=400)  AID, (DLTD(J),J=1,NDLT)
  IF (AID /= 'DLT DATE')                         GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',       ERR=400)  AID, (DLTMAX(J),J=1,NDLT)
  IF (AID /= 'DLT MAX ')                         GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',       ERR=400)  AID, (DLTF(J),J=1,NDLT)
  IF (AID /= 'DLT FRN ')                         GO TO 400
  READ (CON,'(/A8/(8X,3A8))',          ERR=400)  AID, (VISC(JW), CELC(JW),DLTADD(JW), JW=1,NWB)
  IF (AID /= 'DLT LIMI')                         GO TO 400

! Grid definition cards

  WRITE (WIN,*) '  grid definition cards'
  READ (CON,'(/A8/(8X,7I8,2F8.3))',    ERR=400)  AID, (US(JB),DS(JB),UHS(JB),DHS(JB),UQB(JB),DQB(JB),NL(JB),SLOPE(JB),SLOPEC(JB),JB=1,NBR)
  IF (AID /= 'BRANCH G')                        GO TO 400
  READ (CON,'(/A8/(8X,3F8.0,3I8))',   ERR=400)  AID, (LAT(JW),LONG(JW),ELBOT(JW),BS(JW),BE(JW),JBDN(JW),JW=1,NWB)
  IF (AID /= 'LOCATION')                        GO TO 400

! Initial condition cards

  WRITE (WIN,*) '  initial conditions cards'
  READ (CON,'(/A8/(8X,2F8.0,A8,A8))',    ERR=400)  AID, (T2I(JW),    ICEI(JW),   WTYPEC(JW), GRIDCC(JW), JW=1,NWB)
  IF (AID /= 'INIT CND')                        GO TO 400
  READ (CON,'(/A8/(8X,6A8))',         ERR=400)  AID, (VBC(JW),    EBC(JW),    MBC(JW),    PQC(JW),   EVC(JW),     PRC(JW),         &
                                                      JW=1,NWB)
  IF (AID /= 'CALCULAT')                        GO TO 400
  READ (CON,'(/A8/(8X,4A8))',         ERR=400)  AID, (WINDC(JW),  QINC(JW),   QOUTC(JW),  HEATC(JW), JW=1,NWB)
  IF (AID /= 'DEAD SEA')                        GO TO 400
  READ (CON,'(/A8/(8X,3A8))',         ERR=400)  AID, (QINIC(JB),  DTRIC(JB),  HDIC(JB),   JB=1,NBR)
  IF (AID /= 'INTERPOL')                        GO TO 400
  READ (CON,'(/A8/(8X,5A8,4F8.0))',   ERR=400)  AID, (SLHTC(JW),  SROC(JW),   RHEVC(JW),  METIC(JW), FETCHC(JW), AFW(JW),          &
                                                      BFW(JW),    CFW(JW),    WINDH(JW),  JW=1,NWB)
  IF (AID /= 'HEAT EXC')                        GO TO 400
  READ (CON,'(/A8/(8X,2A8,6F8.0))',   ERR=400)  AID, (ICEC(JW),   SLICEC(JW), ALBEDO(JW), HWI(JW),   BETAI(JW),  GAMMAI(JW),       &
                                                      ICEMIN(JW), ICET2(JW),  JW=1,NWB)
  IF (AID /= 'ICE COVE')                        GO TO 400
  READ (CON,'(/A8/(8X,A8,F8.0))',     ERR=400)  AID, (SLTRC(JW),  THETA(JW),  JW=1,NWB)
  IF (AID /= 'TRANSPOR')                        GO TO 400
  READ (CON,'(/A8/(8X,6F8.0,A8,F8.0))',    ERR=400)  AID, (AX(JW), DX(JW), CBHE(JW), TSED(JW), FI(JW), TSEDF(JW), FRICC(JW),Z0(JW), JW=1,NWB)
  IF (AID /= 'HYD COEF')                        GO TO 400

! Viscosity cards

  WRITE (WIN,*) '  viscosity/friction cards'
  READ (CON,'(/A8/(8X,2A8,F8.0,I8,4F8.0,A8))',    ERR=400)  AID, (AZC(JW), AZSLC(JW), AZMAX(JW),FBC(JW),AZE(JW),ARODI(JW),STRCKLR(JW),BOUNDFR(JW),TKECAL(JW),JW=1,NWB)
  IF (AID /= 'EDDY VIS')                        GO TO 400

! Inflow-outflow cards

  WRITE (WIN,*) '  inflow/outflow cards'
  READ (CON,'(/A8/(8X,I8,A8))',          ERR=400)  AID, (NSTR(JB), DYNSTRUC(JB),   JB=1,NBR)
  IF (AID /= 'N STRUC ')                        GO TO 400
  READ (CON,'(/A8/(8X,9A8))',         ERR=400)  AID, (STRIC(JS,1), JS=1,NSTR(1))
  IF (AID /= 'STR INT ')                        GO TO 400
  DO JB=2,NBR
    READ (CON,'(:8X,9A8)',            ERR=400) (STRIC(JS,JB),      JS=1,NSTR(JB))
  END DO
  READ (CON,'(/A8/(:8X,9I8))',        ERR=400)  AID, (KTSTR(JS,1), JS=1,NSTR(1))
  IF (AID /= 'STR TOP ')                        GO TO 400
  DO JB=2,NBR
    READ (CON,'(:8X,9I8)',            ERR=400) (KTSTR(JS,JB),      JS=1,NSTR(JB))
  END DO
  READ (CON,'(/A8/(:8X,9I8))',        ERR=400)  AID, (KBSTR(JS,1), JS=1,NSTR(1))
  IF (AID /= 'STR BOT ')                        GO TO 400
  DO JB=2,NBR
    READ (CON,'(:8X,9I8)',            ERR=400) (KBSTR(JS,JB),      JS=1,NSTR(JB))
  END DO
  READ (CON,'(/A8/(:8X,9A8))',        ERR=400)  AID, (SINKC(JS,1), JS=1,NSTR(1))
  IF (AID /= 'STR SINK')                        GO TO 400
  DO JB=2,NBR
    READ (CON,'(:8X,9A8)',            ERR=400) (SINKC(JS,JB),      JS=1,NSTR(JB))
  END DO
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (ESTR(JS,1),  JS=1,NSTR(1))
  IF (AID /= 'STR ELEV')                        GO TO 400
  DO JB=2,NBR
    READ (CON,'(:8X,9F8.0)',          ERR=400) (ESTR(JS,JB),       JS=1,NSTR(JB))
  END DO
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (WSTR(JS,1),  JS=1,NSTR(1))
  IF (AID /= 'STR WIDT')                        GO TO 400
  DO JB=2,NBR
    READ (CON,'(:8X,9F8.0)',          ERR=400) (WSTR(JS,JB),       JS=1,NSTR(JB))
  END DO
  READ (CON,'(/A8/(:8X,2I8,6F8.0,a8,a8))',   ERR=400)  AID, (IUPI(J),   IDPI(J),  EUPI(J),   EDPI(J),   WPI(J),    DLXPI(J),  FPI(J),     &
                                                      FMINPI(J), latpic(j), dynpipe(j), J=1,NPI)
  IF (AID /= 'PIPES   ')                        GO TO 400
  READ (CON,'(/A8/(:8X,A8,2F8.0,2I8))',ERR=400)  AID, (PUPIC(J),  ETUPI(J), EBUPI(J),  KTUPI(J),  KBUPI(J),  J=1,NPI)
  IF (AID /= 'PIPE UP ')                        GO TO 400
  READ (CON,'(/A8/(:8X,A8,2F8.0,2I8))',ERR=400)  AID, (PDPIC(J),  ETDPI(J), EBDPI(J),  KTDPI(J),  KBDPI(J),  J=1,NPI)
  IF (AID /= 'PIPE DOW')                        GO TO 400
  READ (CON,'(/A8/(:8X,2I8,5F8.0,a8))',  ERR=400)  AID, (IUSP(J),   IDSP(J),  ESP(J),    A1SP(J),   B1SP(J),   A2SP(J),   B2SP(J),    &
                                                      latSPc(j), J=1,NSP)
  IF (AID /= 'SPILLWAY')                        GO TO 400
  READ (CON,'(/A8/(:8X,A8,2F8.0,2I8))',ERR=400)  AID, (PUSPC(J),  ETUSP(J), EBUSP(J),  KTUSP(J),  KBUSP(J),  J=1,NSP)
  IF (AID /= 'SPILL UP')                        GO TO 400
  READ (CON,'(/A8/(:8X,A8,2F8.0,2I8))',ERR=400)  AID, (PDSPC(J),  ETDSP(J), EBDSP(J),  KTDSP(J),  KBDSP(J),  J=1,NSP)
  IF (AID /= 'SPILL DO')                        GO TO 400
  READ (CON,'(/A8/(:8X,A8,I8,3F8.0))', ERR=400)  AID, (GASSPC(J), EQSP(J),  AGASSP(J), BGASSP(J), CGASSP(J), J=1,NSP)
  IF (AID /= 'SPILL GA')                        GO TO 400
  READ (CON,'(/A8/(:8X,2I8,7F8.0,a8))',   ERR=400)  AID, (IUGT(J),   IDGT(J),  EGT(J),    A1GT(J),   B1GT(J),   G1GT(J),   A2GT(J),    &
                                                      B2GT(J),   G2GT(J), latgtc(j), J=1,NGT)   ! sw 7/14/06
  IF (AID /= 'GATES   ')                        GO TO 400
  READ (CON,'(/A8/(:8X,4F8.0,a8,a8))',       ERR=400)  AID, (GTA1(J),   GTB1(J),  GTA2(J),   GTB2(J),  dyngtc(j),gtic(j),        J=1,NGT)  ! sw 7/14/06
  IF (AID /= 'GATE WEI')                        GO TO 400
  READ (CON,'(/A8/(:8X,A8,2F8.0,2I8))',ERR=400)  AID, (PUGTC(J),  ETUGT(J), EBUGT(J),  KTUGT(J),  KBUGT(J),  J=1,NGT)
  IF (AID /= 'GATE UP ')                        GO TO 400
  READ (CON,'(/A8/(:8X,A8,2F8.0,2I8))',ERR=400)  AID, (PDGTC(J),  ETDGT(J), EBDGT(J),  KTDGT(J),  KBDGT(J),  J=1,NGT)
  IF (AID /= 'GATE DOW')                        GO TO 400
  READ (CON,'(/A8/(:8X,A8,I8,3F8.0))', ERR=400)  AID, (GASGTC(J), EQGT(J),  AGASGT(J), BGASGT(J), CGASGT(J), J=1,NGT)
  IF (AID /= 'GATE GAS')                        GO TO 400
  READ (CON,'(/A8/(:8X,2I8,6F8.0,2a8))',   ERR=400)  AID, (IUPU(J),   IDPU(J),  EPU(J),    STRTPU(J), ENDPU(J),  EONPU(J),  EOFFPU(J),  &
                                                      QPU(J),  latpuc(j), DYNPUMP(J),  J=1,NPU)    
  
  !PUMP_DOWNSTREAM=.FALSE.
  !  DO J=1,NPU
  !    IF(IDPU(J) < 0)THEN
  !        PUMP_DOWNSTREAM(J)=.TRUE.
  !        IDPU(J)=ABS(IDPU(J))
  !    ENDIF
  !  ENDDO     
  
  ! RA 8/14/07
  IF (AID /= 'PUMPS 1')                         GO TO 400
  READ (CON,'(/A8/(:8X,A8,2F8.0,2I8))',ERR=400)  AID, (PPUC(J),   ETPU(J),  EBPU(J),   KTPU(J),   KBPU(J),  J=1,NPU)
  IF (AID /= 'PUMPS 2')                         GO TO 400
  READ (CON,'(/A8/(:8X,9I8))',        ERR=400)  AID, (IWR(JW),  JW=1,NIW)
  IF (AID /= 'WEIR SEG')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',        ERR=400)  AID, (EKTWR(JW), JW=1,NIW)
  IF (AID /= 'WEIR TOP')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',        ERR=400)  AID, (EKBWR(JW), JW=1,NIW)
  IF (AID /= 'WEIR BOT')                        GO TO 400
  READ (CON,'(/A8/(:8X,9A8))',        ERR=400)  AID, (WDIC(JW), JW=1,NWD)
  IF (AID /= 'WD INT  ')                        GO TO 400
  READ (CON,'(/A8/(:8X,9I8))',        ERR=400)  AID, (IWD(JW),  JW=1,NWD)
  IF (AID /= 'WD SEG  ')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (EWD(JW),  JW=1,NWD)
  IF (AID /= 'WD ELEV ')                        GO TO 400
  READ (CON,'(/A8/(:8X,9I8))',        ERR=400)  AID, (KTWD(JW), JW=1,NWD)
  IF (AID /= 'WD TOP  ')                        GO TO 400
  READ (CON,'(/A8/(:8X,9I8))',        ERR=400)  AID, (KBWD(JW), JW=1,NWD)
  IF (AID /= 'WD BOT  ')                        GO TO 400
  READ (CON,'(/A8/(:8X,9A8))',        ERR=400)  AID, (PTRC(JT), JT=1,NTR)
  IF (AID /= 'TRIB PLA')                        GO TO 400
  READ (CON,'(/A8/(:8X,9A8))',        ERR=400)  AID, (TRIC(JW), JW=1,NTR)
  IF (AID /= 'TRIB INT')                        GO TO 400
  READ (CON,'(/A8/(:8X,9I8))',        ERR=400)  AID, (ITR(JT),  JT=1,NTR)
  IF (AID /= 'TRIB SEG')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (ETTR(JT), JT=1,NTR)
  IF (AID /= 'TRIB TOP')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (EBTR(JT), JT=1,NTR)
  IF (AID /= 'TRIB BOT')                        GO TO 400
  READ (CON,'(/A8/(8X,A8))',          ERR=400)  AID, (DTRC(JB), JB=1,NBR)
  IF (AID /= 'DST TRIB')                        GO TO 400
  
ELSE
  WRITE (WIN,*) '  time control TMSTR'
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400) AID
    IF (AID /= 'TMSTRT  ') THEN
     AID1='TMSTRT  ';  ICST=3;GO TO 400
    ENDIF
    
  !READ (CON, *, ERR=400) HEADER(1)
  !   DO J=1,15   ! CHECK FIRST 12 CHARACTERS TO SEE IF ':'
  !       IF(HEADER(1)(J:J)==':')THEN
  !           CHAR1='$'
  !           EXIT
  !       ENDIF
  !   ENDDO
  !   BACKSPACE(1010)
      
  READ (CON,*,ERR=400)  TMSTRT,   TMEND,    YEAR    
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400) AID
    IF (AID /= 'NDLT    ') THEN
    AID1='NDLT    '; ICST=3;GO TO 400
    ENDIF
    
  READ (CON,*,ERR=400)  NDLT,     DLTMIN, DLTINTR; DLTD=0.0; DLTINTR=ADJUSTR(DLTINTR)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)  (DLTD(J), J =1,NDLT)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)  (DLTMAX(J), J =1,NDLT)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400) AID
    IF (AID /= 'DLTF    ') THEN
    AID1='DLTF    ';  ICST=3;GO TO 400
    ENDIF
    
  READ (CON,*,ERR=400)  (DLTF(J),   J =1,NDLT)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)
  WRITE (WIN,*) '  viscosity/friction VISC'
  READ (CON,*,ERR=400)      (VISC(JW), JW=1,NWB)
  READ (CON,*,ERR=400)      (CELC(JW), JW=1,NWB)
  READ (CON,*,ERR=400)      (DLTADD(JW), JW=1,NWB)
  VISC=ADJUSTR(VISC);   CELC=ADJUSTR(CELC);   DLTADD=ADJUSTR(DLTADD)

! Grid definition cards
  READ (CON,*,ERR=400)
  WRITE (WIN,*) '  grid geometry US'
  READ (CON,*,ERR=400) AID
      IF (AID /= 'BR1     ')THEN
      AID1='BR1     ';  ICST=3;GO TO 400
      ENDIF

  READ (CON,*,ERR=400)  (US(JB),    JB=1,NBR)
  READ (CON,*,ERR=400)  (DS(JB),    JB=1,NBR)
  READ (CON,*,ERR=400)  (UHS(JB),   JB=1,NBR)
  READ (CON,*,ERR=400)  (DHS(JB),   JB=1,NBR)
  READ (CON,*,ERR=400)  (NL(JB),    JB=1,NBR)
  READ (CON,*,ERR=400)  (SLOPE(JB), JB=1,NBR)
  READ (CON,*,ERR=400)  (SLOPEC(JB),JB=1,NBR)
  READ (CON,*,ERR=400)
  WRITE (WIN,*) '  location LAT'
  READ (CON,*,ERR=400) AID
      IF (AID /= 'WB1     ') THEN
        AID1='WB1     '
        ICST=3;GO TO 400
      ENDIF

  READ (CON,*,ERR=400)   (LAT(JW),                     JW=1,NWB)
  READ (CON,*,ERR=400)   (LONG(JW),                    JW=1,NWB)
  READ (CON,*,ERR=400)   (ELBOT(JW),                   JW=1,NWB)
  READ (CON,*,ERR=400)   (BS(JW),                      JW=1,NWB)
  READ (CON,*,ERR=400)   (BE(JW),                      JW=1,NWB)
  READ (CON,*,ERR=400)   (JBDN(JW),                    JW=1,NWB)
  READ (CON,*,ERR=400)
  WRITE (WIN,*) '  initial temperature ice T2I'
  READ (CON,*,ERR=400)AID
    IF (AID /= 'WB1     ') THEN
         AID1='WB1     '
        ICST=3;GO TO 400
    ENDIF

! Initial condition cards

  READ (CON,*,ERR=400)      (T2I(JW),                              JW=1,NWB)
  READ (CON,*,ERR=400)      (ICEI(JW),                             JW=1,NWB)
  READ (CON,*,ERR=400)      (WTYPEC(JW),                           JW=1,NWB); WTYPEC=ADJUSTR(WTYPEC)
  READ (CON,*,ERR=400)      (GRIDCC(JW),                           JW=1,NWB); GRIDCC=ADJUSTR(GRIDCC)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)
  WRITE (WIN,*) '  volume heat balance VBC'

  READ (CON,*,ERR=400)            (VBC(JW),              JW=1,NWB); VBC=ADJUSTR(VBC)
  READ (CON,*,ERR=400)            (EBC(JW),              JW=1,NWB); EBC=ADJUSTR(EBC)
  READ (CON,*,ERR=400)            (MBC(JW),              JW=1,NWB); MBC=ADJUSTR(MBC)   
  READ (CON,*,ERR=400)            (PQC(JW),              JW=1,NWB); PQC=ADJUSTR(PQC)
  READ (CON,*,ERR=400)            (EVC(JW),              JW=1,NWB); EVC=ADJUSTR(EVC)
  READ (CON,*,ERR=400)            (PRC(JW),              JW=1,NWB); PRC=ADJUSTR(PRC)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)
  WRITE (WIN,*) '  dead sea WINDC'

  READ (CON,*,ERR=400)            (WINDC(JW),             JW=1,NWB); WINDC=ADJUSTR(WINDC)
  READ (CON,*,ERR=400)            (QINC(JW),              JW=1,NWB); QINC=ADJUSTR(QINC)
  READ (CON,*,ERR=400)            (QOUTC(JW),             JW=1,NWB); QOUTC=ADJUSTR(QOUTC)
  READ (CON,*,ERR=400)            (HEATC(JW),             JW=1,NWB); HEATC=ADJUSTR(HEATC)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)
  WRITE (WIN,*) '  interpolation QINIC'

  READ (CON,*,ERR=400)            (QINIC(JB),               JB=1,NBR); QINIC=ADJUSTR(QINIC)
  READ (CON,*,ERR=400)            (DTRIC(JB),               JB=1,NBR); DTRIC=ADJUSTR(DTRIC)
  READ (CON,*,ERR=400)            (HDIC(JB),                JB=1,NBR); HDIC=ADJUSTR(HDIC)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)
  WRITE (WIN,*) '  meteorological SLHTC '

  READ (CON,*,ERR=400)      (SLHTC(JW),                     JW=1,NWB); SLHTC=ADJUSTR(SLHTC)
  READ (CON,*,ERR=400)      (SROC(JW),                      JW=1,NWB); SROC=ADJUSTR(SROC)
  READ (CON,*,ERR=400)      (RHEVC(JW),                     JW=1,NWB); RHEVC=ADJUSTR(RHEVC)
  READ (CON,*,ERR=400)      (METIC(JW),                     JW=1,NWB); METIC=ADJUSTR(METIC)
  READ (CON,*,ERR=400)      (FETCHC(JW),                    JW=1,NWB); FETCHC=ADJUSTR(FETCHC)
  READ (CON,*,ERR=400)      (AFW(JW),                       JW=1,NWB)
  READ (CON,*,ERR=400)      (BFW(JW),                       JW=1,NWB)
  READ (CON,*,ERR=400)      (CFW(JW),                       JW=1,NWB)
  READ (CON,*,ERR=400)      (WINDH(JW),                     JW=1,NWB)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)
  WRITE (WIN,*) '  ice cover ICEC'

  READ (CON,*,ERR=400)      (ICEC(JW),                      JW=1,NWB); ICEC=ADJUSTR(ICEC)
  READ (CON,*,ERR=400)      (SLICEC(JW),                    JW=1,NWB); SLICEC=ADJUSTR(SLICEC)
  READ (CON,*,ERR=400)      (ALBEDO(JW),                    JW=1,NWB)
  READ (CON,*,ERR=400)      (HWI(JW),                       JW=1,NWB)
  READ (CON,*,ERR=400)      (BETAI(JW),                     JW=1,NWB)
  READ (CON,*,ERR=400)      (GAMMAI(JW),                    JW=1,NWB)
  READ (CON,*,ERR=400)      (ICEMIN(JW),                    JW=1,NWB)
  READ (CON,*,ERR=400)      (ICET2(JW),                     JW=1,NWB)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)

  READ (CON,*,ERR=400)        (SLTRC(JW),                   JW=1,NWB); SLTRC=ADJUSTR(SLTRC)
  READ (CON,*,ERR=400)        (THETA(JW),                   JW=1,NWB)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)

    WRITE (WIN,*) '  longitudinal viscosity AX'

  READ (CON,*,ERR=400)       (AX(JW),                       JW=1,NWB)
  READ (CON,*,ERR=400)       (DX(JW),                       JW=1,NWB)
  READ (CON,*,ERR=400)       (CBHE(JW),                     JW=1,NWB)
  READ (CON,*,ERR=400)       (TSED(JW),                     JW=1,NWB)
  READ (CON,*,ERR=400)       (FI(JW),                       JW=1,NWB)
  READ (CON,*,ERR=400)       (TSEDF(JW),                    JW=1,NWB)
  READ (CON,*,ERR=400)       (FRICC(JW),                    JW=1,NWB); FRICC=ADJUSTR(FRICC)
  READ (CON,*,ERR=400)       (Z0(JW),                       JW=1,NWB)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)

    WRITE (WIN,*) '  turbulence closure AZC'
  READ (CON,*,ERR=400)      (AZC(JW),    JW=1,NWB);AZC=ADJUSTR(AZC)          
  READ (CON,*,ERR=400)      (AZSLC(JW),  JW=1,NWB);AZSLC=ADJUSTR(AZSLC)          
  READ (CON,*,ERR=400)      (AZMAX(JW),  JW=1,NWB)          
  READ (CON,*,ERR=400)      (FBC(JW),  JW=1,NWB)          
  READ (CON,*,ERR=400)      (AZE(JW), JW=1,NWB)          
  READ (CON,*,ERR=400)      (ARODI(JW),  JW=1,NWB)          
  READ (CON,*,ERR=400)      (STRCKLR(JW), JW=1,NWB)          
  READ (CON,*,ERR=400)      (BOUNDFR(JW),JW=1,NWB)   
  READ (CON,*,ERR=400)      (TKECAL(JW),JW=1,NWB);TKECAL=ADJUSTR(TKECAL)          !FBC(JW),AZE(JW),ARODI(JW),STRCKLR(JW),BOUNDFR(JW),TKECAL(JW)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)
  WRITE (WIN,*) '  inflow/outflow NSTR'

  READ (CON,*,ERR=400)            (NSTR(JB),      JB=1,NBR)
  READ (CON,*,ERR=400)            (DYNSTRUC(JB),  JB=1,NBR); DYNSTRUC=adjustr(DYNSTRUC)
  
  DO JS=1,NSTT
    READ (CON,*,ERR=400)        (STRIC(JS,JB),  JB=1,NBR)
  END DO
  STRIC=adjustr(STRIC)
  DO JS=1,NSTT
    READ (CON,*,ERR=400)        (KTSTR(JS,JB),  JB=1,NBR)
  END DO
  DO JS=1,NSTT
    READ (CON,*,ERR=400)        (KBSTR(JS,JB),  JB=1,NBR)
  END DO
  DO JS=1,NSTT
    READ (CON,*,ERR=400)        (SINKC(JS,JB), JB=1,NBR)
  END DO
  SINKC=adjustr(SINKC)
  DO JS=1,NSTT
    READ (CON,*,ERR=400)        (ESTR(JS,JB),  JB=1,NBR)
  END DO
  DO JS=1,NSTT
    READ (CON,*,ERR=400)        (WSTR(JS,JB),  JB=1,NBR)
  END DO
  READ (CON,*,ERR=400)
    WRITE (WIN,*) '  pipes IUPI'
  READ (CON,*,ERR=400)AID
          IF (AID /= 'PIPE1   ') THEN
          AID1='PIPE1   '
          ICST=3;GO TO 400
      ENDIF

  READ (CON,*,ERR=400)  (IUPI(JP),  JP=1,NPI)
  READ (CON,*,ERR=400)  (IDPI(JP),  JP=1,NPI)
  READ (CON,*,ERR=400)  (EUPI(JP),  JP=1,NPI)
  READ (CON,*,ERR=400)  (EDPI(JP),  JP=1,NPI)
  READ (CON,*,ERR=400)  (WPI(JP),   JP=1,NPI)
  READ (CON,*,ERR=400)  (DLXPI(JP), JP=1,NPI)
  READ (CON,*,ERR=400)  (FPI(JP),   JP=1,NPI)
  READ (CON,*,ERR=400)  (FMINPI(JP),JP=1,NPI)
  READ (CON,*,ERR=400)  (LATPIC(JP),JP=1,NPI);LATPIC=adjustr(LATPIC)
  READ (CON,*,ERR=400)  (DYNPIPE(JP),JP=1,NPI);DYNPIPE=adjustr(DYNPIPE) 
  
  READ (CON,*,ERR=400)  (PUPIC(JP),  JP=1,NPI);PUPIC=adjustr(PUPIC)
  READ (CON,*,ERR=400)  (ETUPI(JP),  JP=1,NPI)
  READ (CON,*,ERR=400)  (EBUPI(JP),  JP=1,NPI)
  READ (CON,*,ERR=400)  (KTUPI(JP),  JP=1,NPI)
  READ (CON,*,ERR=400)  (KBUPI(JP),  JP=1,NPI)
  
  READ (CON,*,ERR=400)  (PDPIC(JP),  JP=1,NPI);PDPIC=adjustr(PDPIC)
  READ (CON,*,ERR=400)  (ETDPI(JP),  JP=1,NPI)
  READ (CON,*,ERR=400)  (EBDPI(JP),  JP=1,NPI)
  READ (CON,*,ERR=400)  (KTDPI(JP),  JP=1,NPI)
  READ (CON,*,ERR=400)  (KBDPI(JP),  JP=1,NPI)
  
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)
    WRITE (WIN,*) '  spillways IUSP'

  READ (CON,*,ERR=400) (IUSP(JS),    JS=1,NSP)
  READ (CON,*,ERR=400) (IDSP(JS),    JS=1,NSP)
  READ (CON,*,ERR=400) (ESP(JS),     JS=1,NSP)
  READ (CON,*,ERR=400) (A1SP(JS),    JS=1,NSP)
  READ (CON,*,ERR=400) (B1SP(JS),    JS=1,NSP)
  READ (CON,*,ERR=400) (A2SP(JS),    JS=1,NSP)
  READ (CON,*,ERR=400) (B2SP(JS),    JS=1,NSP)
  READ (CON,*,ERR=400) (LATSPC(JS),  JS=1,NSP);LATSPC=ADJUSTR(LATSPC)
  
  READ (CON,*,ERR=400) (PUSPC(JS),   JS=1,NSP);PUSPC=adjustr(PUSPC)
  READ (CON,*,ERR=400) (ETUSP(JS),   JS=1,NSP)
  READ (CON,*,ERR=400) (EBUSP(JS),   JS=1,NSP)
  READ (CON,*,ERR=400) (KTUSP(JS),   JS=1,NSP)
  READ (CON,*,ERR=400) (KBUSP(JS),   JS=1,NSP)
  
  READ (CON,*,ERR=400) (PDSPC(JS),   JS=1,NSP) ;PDSPC=ADJUSTR(PDSPC) 
  READ (CON,*,ERR=400) (ETDSP(JS),   JS=1,NSP)
  READ (CON,*,ERR=400) (EBDSP(JS),   JS=1,NSP)
  READ (CON,*,ERR=400) (KTDSP(JS),   JS=1,NSP)
  READ (CON,*,ERR=400) (KBDSP(JS),   JS=1,NSP)
  
  READ (CON,*,ERR=400) (GASSPC(JS),  JS=1,NSP);GASSPC=ADJUSTR(GASSPC)
  READ (CON,*,ERR=400) (EQSP(JS),    JS=1,NSP)
  READ (CON,*,ERR=400) (AGASSP(JS),  JS=1,NSP)
  READ (CON,*,ERR=400) (BGASSP(JS),  JS=1,NSP)
  READ (CON,*,ERR=400) (CGASSP(JS),  JS=1,NSP)
  
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)
    WRITE (WIN,*) '  gates IUGT'

  READ (CON,*,ERR=400) (IUGT(JG),   JG=1,NGT)
  READ (CON,*,ERR=400) (IDGT(JG),   JG=1,NGT)
  READ (CON,*,ERR=400) (EGT(JG),    JG=1,NGT)
  READ (CON,*,ERR=400) (A1GT(JG),   JG=1,NGT)
  READ (CON,*,ERR=400) (B1GT(JG),   JG=1,NGT)
  READ (CON,*,ERR=400) (G1GT(JG),   JG=1,NGT)
  READ (CON,*,ERR=400) (A2GT(JG),   JG=1,NGT)
  READ (CON,*,ERR=400) (B2GT(JG),   JG=1,NGT)
  READ (CON,*,ERR=400) (G2GT(JG),   JG=1,NGT)
  READ (CON,*,ERR=400) (LATGTC(JG), JG=1,NGT);LATGTC=ADJUSTR(LATGTC)
  
  READ (CON,*,ERR=400) (GTA1(JG),   JG=1,NGT)  
  READ (CON,*,ERR=400) (GTB1(JG),   JG=1,NGT)  
  READ (CON,*,ERR=400) (GTA2(JG),   JG=1,NGT)  
  READ (CON,*,ERR=400) (GTB2(JG),   JG=1,NGT)  
  READ (CON,*,ERR=400) (DYNGTC(JG), JG=1,NGT);DYNGTC=ADJUSTR(DYNGTC)  
  READ (CON,*,ERR=400) (GTIC(JG),   JG=1,NGT);GTIC=ADJUSTR(GTIC)  
  
  ETUGT=0.0;EBUGT=0.0;ETDGT=0.0;EBDGT=0.0
  
  READ (CON,*,ERR=400) (PUGTC(JG),  JG=1,NGT);PUGTC=ADJUSTR(PUGTC)
  READ (CON,*,ERR=400) (ETUGT(JG),  JG=1,NGT)
  READ (CON,*,ERR=400) (EBUGT(JG),  JG=1,NGT)
  READ (CON,*,ERR=400) (KTUGT(JG),  JG=1,NGT)
  READ (CON,*,ERR=400) (KBUGT(JG),  JG=1,NGT)
  READ (CON,*,ERR=400) (PDGTC(JG),  JG=1,NGT);PDGTC=ADJUSTR(PDGTC)
  READ (CON,*,ERR=400) (ETDGT(JG),  JG=1,NGT)
  READ (CON,*,ERR=400) (EBDGT(JG),  JG=1,NGT)
  READ (CON,*,ERR=400) (KTDGT(JG),  JG=1,NGT)
  READ (CON,*,ERR=400) (KBDGT(JG),  JG=1,NGT)
  
  READ (CON,*,ERR=400) (GASGTC(JG), JG=1,NGT);GASGTC=ADJUSTR(GASGTC)  
  READ (CON,*,ERR=400) (EQGT(JG),   JG=1,NGT)  
  READ (CON,*,ERR=400) (AGASGT(JG), JG=1,NGT)  
  READ (CON,*,ERR=400) (BGASGT(JG), JG=1,NGT)  
  READ (CON,*,ERR=400) (CGASGT(JG), JG=1,NGT)  
  
  READ (CON,*,ERR=400)
    WRITE (WIN,*) '  pumps IUPU'
  READ (CON,*,ERR=400)AID
          IF (AID /= 'Pump1   ') THEN
          AID1='Pump1   '
          ICST=3;GO TO 400
      ENDIF
  
  READ (CON,*,ERR=400) (IUPU(JP),   JP=1,NPU)
  READ (CON,*,ERR=400) (IDPU(JP),   JP=1,NPU)
  READ (CON,*,ERR=400) (EPU(JP),    JP=1,NPU)
  READ (CON,*,ERR=400) (STRTPU(JP), JP=1,NPU)
  READ (CON,*,ERR=400) (ENDPU(JP),  JP=1,NPU)
  READ (CON,*,ERR=400) (EONPU(JP),  JP=1,NPU)
  READ (CON,*,ERR=400) (EOFFPU(JP), JP=1,NPU)
  READ (CON,*,ERR=400) (QPU(JP),    JP=1,NPU)
  READ (CON,*,ERR=400) (LATPUC(JP), JP=1,NPU);LATPUC=adjustr(LATPUC)
  READ (CON,*,ERR=400) (DYNPUMP(JP),JP=1,NPU);DYNPUMP=ADJUSTR(DYNPUMP)
  READ (CON,*,ERR=400) (PPUC(JP),   JP=1,NPU);PPUC=ADJUSTR(PPUC)
  READ (CON,*,ERR=400) (ETPU(JP),   JP=1,NPU)
  READ (CON,*,ERR=400) (EBPU(JP),   JP=1,NPU)
  READ (CON,*,ERR=400) (KTPU(JP),   JP=1,NPU)
  READ (CON,*,ERR=400) (KBPU(JP),   JP=1,NPU)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)
  WRITE (WIN,*) '  internal weirs IWR'

  READ (CON,*,ERR=400)         (IWR(JW),    JW=1,NIW)
  READ (CON,*,ERR=400)        (EKTWR(JW),   JW=1,NIW)               
  READ (CON,*,ERR=400)        (EKBWR(JW),   JW=1,NIW)   
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)

    WRITE (WIN,*) '  withdrawals'

  READ (CON,*,ERR=400)         (WDIC(JW),   JW=1,NWD);WDIC=adjustr(WDIC)
  READ (CON,*,ERR=400)         (IWD(JW),    JW=1,NWD)
  READ (CON,*,ERR=400)         (EWD(JW),    JW=1,NWD)
  READ (CON,*,ERR=400)         (KTWD(JW),   JW=1,NWD)
  READ (CON,*,ERR=400)         (KBWD(JW),   JW=1,NWD); PTRC= '      '  ! SW 9/27/13 INITIALIZATION SINCE ALLOCATION IS TO NTRT
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)

    WRITE (WIN,*) '  tributaries'

  READ (CON,*,ERR=400)        (PTRC(JT),    JT=1,NTR);PTRC=ADJUSTR(PTRC)
  READ (CON,*,ERR=400)        (TRIC(JT),   JT=1,NTR);TRIC=ADJUSTR(TRIC)
  READ (CON,*,ERR=400)        (ITR(JT),    JT=1,NTR)
  READ (CON,*,ERR=400)        (ETTR(JT),  JT=1,NTR)
  READ (CON,*,ERR=400)        (EBTR(JT),  JT=1,NTR)
  READ (CON,*,ERR=400)        (QTRFN(JT),  JT=1,NTR)
  READ (CON,*,ERR=400)        (TTRFN(JT),  JT=1,NTR)
  READ (CON,*,ERR=400)        (CTRFN(JT),  JT=1,NTR)  
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)
  WRITE (WIN,*) '  distributed tributaries'

  READ (CON,*,ERR=400)      (DTRC(JB),   JB=1,NBR);DTRC=adjustr(DTRC)
  READ (CON,*,ERR=400)
  WRITE (WIN,*) '  output control HNAME'
  READ (CON,*,ERR=400)AID
          IF (AID /= 'HNAME   ') THEN
          AID1='HNAME   '
          ICST=3;GO TO 400
      ENDIF
ENDIF
  PUMP_DOWNSTREAM=.FALSE.
    DO J=1,NPU
      IF(IDPU(J) < 0)THEN
          PUMP_DOWNSTREAM(J)=.TRUE.
          IDPU(J)=ABS(IDPU(J))
      ENDIF
    ENDDO     
    
! Bathymetry definition CHECKS BEFORE PROCEEDING FURTHER
  WRITE (WIN,*) '  grid check: US, DS, BE, BS'
  AID='GridChek'
! Check bathymetry
  isum=0
  DO JW=1,NWB
    isum=isum+be(jw)-bs(jw)+1
  END DO
  if(isum /= nbr)then
  CALL ERRORS
  WRITE (ERR,'(A)') 'There is an error in BS and BE (or in NWB): The sum of BS-BE+1 for all waterbodies does not equal NBR. Check branch geometry and waterbody branch definition.'
  GO TO 405
  endif
  
  IF(NWB>1)THEN
      DO JW=1,NWB-1
      IF(BE(JW)>BE(JW+1))THEN
        CALL ERRORS
      WRITE (ERR,'(A,i3,A,I3,A,I3,A,I3,A)') 'Branch geometry and waterbody branch designation error. For WB:',jw ,' BE is branch:',be(jw),' but for WB:',jw+1,' BE is branch:',be(jw+1),'. BE(JW) must always be < BE(JW+1).'
      GO TO 405
      ENDIF
      IF(BS(JW)>BS(JW+1))THEN
        CALL ERRORS
      WRITE (ERR,'(A,i3,A,I3,A,I3,A,I3,A)') 'Branch geometry and waterbody branch designation error. For WB:',jw ,' BS is branch:',bS(jw),' but for WB:',jw+1,' BS is branch:',bS(jw+1),'. BS(JW) must always be < BS(JW+1).'
      GO TO 405
      ENDIF
      ENDDO
  ENDIF

  ! check to make sure that all US and DS have inactive segments on either side of them

  IF(US(1) /= 2)THEN
  CALL ERRORS
  WRITE(ERR,*)'Branch GEOM: US(1) must equal 2 since it is the first branch and the first segment is inactive'

  END IF
    DO JB=1,NBR-1
    DIFF=US(JB+1)-DS(JB)
    IF(DIFF /= 3)THEN
     CALL ERRORS
     WRITE(ERR,'(A,I3,A,I3)')'Branch GEOM: There must be inactive segments at the end and beginning of each branch. There are not enough or too many between branch:',jb,' and branch:',jb+1

    ENDIF
    END DO
  IF(DS(BE(NWB))+1 /= IMX)THEN
  CALL ERRORS
  WRITE(ERR,'(A)')'Branch GEOM: DS(NBR)+1 must equal IMX - check if DS of last branch is correct.'

  ENDIF

IF(CONFN=='w2_con.npt')THEN

! Output control cards (excluding constituents)

  WRITE (WIN,*) '  output control cards'
  READ (CON,'(/A8/(:8X,9A8))',        ERR=400)  AID, (HPRWBC(1,JW),JW=1,NWB)
  IF (AID /= 'HYD PRIN')                        GO TO 400
  DO JH=2,NHY
    READ (CON,'(:8X,9A8)',            ERR=400)  (HPRWBC(JH,JW), JW=1,NWB)
  END DO
  READ (CON,'(/A8/(:8X,A8,2I8))',     ERR=400)  AID, (SNPC(J), NSNP(J), NISNP(J), J=1,NWB)
  IF (AID /= 'SNP PRIN')                        GO TO 400

  !do j=1,nwb
  !if(nsnp(j) > nod)then
  !     CALL ERRORS
  !     WRITE (ERR,'(A,i3,A,i3,a,i3)') 'The # of output dates for Snapshot output file NSNP[',nsnp(j),'] is > NDAY[',nod,'] for waterbody:',j
  !     WRITE (WIN,*) 'Error - snapshot output - see pre.err'
  !     stop
  !endif
  !if(nisnp(j) > imx)then
  !     CALL ERRORS
  !     WRITE (ERR,'(A,i3,A,i3,a,i3)') 'The # segments for SNP output file NISNP[',nisnp(j),'] is > IMX[',imx,'], the max # of segments for waterbody:',j
  !     WRITE (WIN,*) 'Error - snapshot output - see pre.err'
  !     stop
  !endif
  !enddo

  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (SNPD(J,1),J=1,NSNP(1))
  IF (AID /= 'SNP DATE')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9F8.0)',          ERR=400) (SNPD(J,JW),     J=1,NSNP(JW))
  END DO
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (SNPF(J,1),J=1,NSNP(1))
  IF (AID /= 'SNP FREQ')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9F8.0)',          ERR=400) (SNPF(J,JW),     J=1,NSNP(JW))
  END DO
  READ (CON,'(/A8/(:8X,9I8))',        ERR=400)  AID, (ISNP(J,1),J=1,NISNP(1))
  IF (AID /= 'SNP SEG ')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9I8)',            ERR=400) (ISNP(J,JW),     J=1,NISNP(JW))
  END DO
  READ (CON,'(/A8/(8X,A8,I8))',       ERR=400)  AID, (SCRC(JW), NSCR(JW), JW=1,NWB)
  IF (AID /= 'SCR PRIN')                        GO TO 400
  READ (CON,'(/A8/(8X,9F8.0))',       ERR=400)  AID, (SCRD(J,1),J=1,NSCR(1))
  IF (AID /= 'SCR DATE')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9F8.0)',          ERR=400) (SCRD(J,JW),     J=1,NSCR(JW))
  END DO
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (SCRF(J,1),J=1,NSCR(1))
  IF (AID /= 'SCR FREQ')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9F8.0)',          ERR=400) (SCRF(J,JW),     J=1,NSCR(JW))
  END DO
  READ (CON,'(/A8/(8X,A8,2I8))',      ERR=400)  AID, (PRFC(JW), NPRF(JW), NIPRF(JW), JW=1,NWB)
  
  !do j=1,nwb
  !if(nprf(j) > nod)then
  !     CALL ERRORS
  !     WRITE (ERR,'(A,i3,A,i3,a,i3)') 'The # of output dates for Profile [PRF] output file NPRF[',nprf(j),'] is > NDAY[',nod,'] for waterbody:',j
  !     WRITE (WIN,*) 'Error - profile output - see pre.err'
  !     stop
  !endif
  !if(niprf(j) > imx)then
  !     CALL ERRORS
  !     WRITE (ERR,'(A,i3,A,i3,a,i3)') 'The # segments for Profile [PRF] output file NIPRF[',niprf(j),'] is > IMX[',imx,'], the max # of segments for waterbody:',j
  !     WRITE (WIN,*) 'Error - profile output - see pre.err'
  !     stop
  !endif
  !enddo
  
  
  
  
  
  IF (AID /= 'PRF PLOT')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (PRFD(J,1),J=1,NPRF(1))
  IF (AID /= 'PRF DATE')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9F8.0)',          ERR=400) (PRFD(J,JW),     J=1,NPRF(JW))
  END DO
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (PRFF(J,1),J=1,NPRF(1))
  IF (AID /= 'PRF FREQ')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9F8.0)',          ERR=400) (PRFF(J,JW),     J=1,NPRF(JW))
  END DO
  READ (CON,'(/A8/(:8X,9I8))',        ERR=400)  AID, (IPRF(J,1),J=1,NIPRF(1))
  IF (AID /= 'PRF SEG ')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9I8)',            ERR=400) (IPRF(J,JW),     J=1,NIPRF(JW))
  END DO
  READ (CON,'(/A8/(8X,A8,2I8))',      ERR=400)  AID, (SPRC(JW), NSPR(JW), NISPR(JW), JW=1,NWB)  
  
  IF (AID /= 'SPR PLOT')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (SPRD(J,1),J=1,NSPR(1))
  IF (AID /= 'SPR DATE')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9F8.0)',          ERR=400) (SPRD(J,JW),     J=1,NSPR(JW))
  END DO
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (SPRF(J,1),J=1,NSPR(1))
  IF (AID /= 'SPR FREQ')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9F8.0)',          ERR=400) (SPRF(J,JW),     J=1,NSPR(JW))
  END DO
  READ (CON,'(/A8/(:8X,9I8))',        ERR=400)  AID, (ISPR(J,1),J=1,NISPR(1))
  IF (AID /= 'SPR SEG ')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9I8)',            ERR=400) (ISPR(J,JW),     J=1,NISPR(JW))
  END DO
  READ (CON,'(/A8/(:8X,A8,I8))',      ERR=400)  AID, (VPLC(JW), NVPL(JW), JW=1,NWB)
  IF (AID /= 'VPL PLOT')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (VPLD(J,1),J=1,NVPL(1))
  IF (AID /= 'VPL DATE')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9F8.0)',          ERR=400) (VPLD(J,JW),     J=1,NVPL(JW))
  END DO
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (VPLF(J,1),J=1,NVPL(1))
  IF (AID /= 'VPL FREQ')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9F8.0)',          ERR=400) (VPLF(J,JW),     J=1,NVPL(JW))
  END DO
  READ (CON,'(/A8/(8X,A8,I8,A8))',       ERR=400)  AID, (CPLC(JW), NCPL(JW),TCPL(JW), JW=1,NWB)
  IF (AID /= 'CPL PLOT')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (CPLD(J,1),J=1,NCPL(1))
  IF (AID /= 'CPL DATE')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9F8.0)',          ERR=400) (CPLD(J,JW),     J=1,NCPL(JW))
  END DO
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (CPLF(J,1),J=1,NCPL(1))
  IF (AID /= 'CPL FREQ')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9F8.0)',          ERR=400) (CPLF(J,JW),     J=1,NCPL(JW))
  END DO
  READ (CON,'(/A8/(8X,A8,I8))',       ERR=400)  AID, (FLXC(JW), NFLX(JW), JW=1,NWB)
  IF (AID /= 'FLUXES  ')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (FLXD(J,1),J=1,NFLX(1))
  IF (AID /= 'FLX DATE')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9F8.0)',          ERR=400) (FLXD(J,JW),     J=1,NFLX(JW))
  END DO
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (FLXF(J,1),J=1,NFLX(1))
  IF (AID /= 'FLX FREQ')                        GO TO 400
  DO JW=2,NWB
    READ (CON,'(:8X,9F8.0)',          ERR=400) (FLXF(J,JW),     J=1,NFLX(JW))
  END DO
  READ (CON,'(/A8/(8X,A8,2I8))',      ERR=400)  AID, TSRC, NTSR, NIKTSR
  
  IF (AID /= 'TSR PLOT')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (TSRD(J),  J=1,NTSR)
  IF (AID /= 'TSR DATE')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (TSRF(J),  J=1,NTSR)
  IF (AID /= 'TSR FREQ')                        GO TO 400
  READ (CON,'(/A8/(:8X,9I8))',        ERR=400)  AID, (ITSR(J),  J=1,NIKTSR)
  IF (AID /= 'TSR SEG ')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (ETSR(J),  J=1,NIKTSR)
  IF (AID /= 'TSR LAYE' .AND. AID /= 'TSR ELEV')GO TO 400

  READ (CON,'(/A8/(:8X,A8,9F8.0))',      ERR=400)  AID, WLC,WLF
  IF (AID /= 'WLOUT   ')GO TO 400
  READ (CON,'(/A8/(:8X,A8,9F8.0))',      ERR=400)  AID, FLOWBALC,FLOWBALF
  IF (AID /= 'FLOWBAL ')GO TO 400
  READ (CON,'(/A8/(:8X,A8,9F8.0))',      ERR=400)  AID, NPBALC,NPBALF
  IF (AID /= 'NPBAL   ')GO TO 400  
  
  READ (CON,'(/A8/8X,A8,2I8)',        ERR=400)  AID, WDOC, NWDO, NIWDO
  IF (AID /= 'WITH OUT')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (WDOD(J),  J=1,NWDO)
  IF (AID /= 'WITH DAT')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (WDOF(J),  J=1,NWDO)
  IF (AID /= 'WITH FRE')                        GO TO 400
  READ (CON,'(/A8/(:8X,9I8))',        ERR=400)  AID, (IWDO(J),  J=1,NIWDO)
  IF (AID /= 'WITH SEG')                        GO TO 400
  READ (CON,'(/A8/8X,A8,I8,A8)',      ERR=400)  AID, RSOC, NRSO, RSIC
  IF (AID /= 'RESTART ')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (RSOD(J),  J=1,NRSO)
  IF (AID /= 'RSO DATE')                        GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (RSOF(J),  J=1,NRSO)
  IF (AID /= 'RSO FREQ')                        GO TO 400
ELSE
   ! Output control cards (excluding constituents)
  DO JH=1,NHY
    READ (CON,*,ERR=400)    HNAME(JH),  FMTH(JH),  HMULT(JH),(HPRWBC(JH,JW),JW=1,NWB);HPRWBC=ADJUSTR(HPRWBC)
  END DO
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)
    WRITE (WIN,*) '  SNP file control SNPC'

  READ (CON,*,ERR=400)        SNPC(1)
  READ (CON,*,ERR=400)        NSNP(1)
      SNPC(2:NWB)=SNPC(1);SNPC=ADJUSTR(SNPC)
      NSNP(2:NWB)=NSNP(1)
      
  NISNP=0     ! SW 3/31/2020
  DO JW=1,NWB
      DO JB=BS(JW),BE(JW)
          DO I=US(JB),DS(JB)
              NISNP(JW)=NISNP(JW)+1
              ISNP(NISNP(JW),JW)=I
          ENDDO
      ENDDO
  ENDDO
    READ (CON,*,ERR=400)      (SNPD(J,1),J=1,NSNP(1))
    READ (CON,*,ERR=400)      (SNPF(J,1),J=1,NSNP(1))
    DO J=1,NSNP(1)
    SNPD(J,2:NWB)=SNPD(J,1)
    SNPF(J,2:NWB)=SNPF(J,1)
    ENDDO

  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400) 
  READ (CON,*,ERR=400)         SCRC(1);  SCRC(1)=ADJUSTR(SCRC(1))    
  READ (CON,*,ERR=400)         NSCR(1)
  READ (CON,*,ERR=400)          (SCRD(J,1),J=1,NSCR(1))
  READ (CON,*,ERR=400)          (SCRF(J,1),J=1,NSCR(1))
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)  
  
  IF(NWB > 1)THEN
      SCRC(2:NWB)=SCRC(1)
      NSCR(2:NWB)=NSCR(1)
  DO J=1,NSCR(1)
      SCRD(J,2:NWB)=SCRD(J,1)
      SCRF(J,2:NWB)=SCRF(J,1)
  ENDDO
  ENDIF
   
  CDUM='     OFF'
  READ (CON,*,ERR=400)      CDUM;CDUM=ADJUSTR(CDUM)         
  READ (CON,*,ERR=400)      NDUM          
  READ (CON,*,ERR=400)      NIDUM          
  
  READ (CON,*,ERR=400)      (DDUM(J),J=1,NDUM)    
  READ (CON,*,ERR=400)      (FDUM(J),J=1,NDUM)    
  READ (CON,*,ERR=400)      (IDUM(J),J=1,NIDUM)   
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)  
  
  PRFC='     OFF'
  NPRF=0
  NIPRF=0
  IF(CDUM == '      ON')THEN                      
  DO J=1,NDUM
  PRFD(J,1:NWB)=DDUM(J)
  PRFF(J,1:NWB)=FDUM(J)
  ENDDO
  NPRF(1:NWB)=NDUM   
      DO JW=1,NWB 
          JJ=0
          DO J=1,NIDUM  
              IF(IDUM(J) >= US(BS(JW)) .AND. IDUM(J) <= DS(BE(JW)))THEN
                  JJ=JJ+1
                  IPRF(JJ,JW)=IDUM(J)
                  NIPRF(JW)=JJ
                  PRFC(JW)='      ON'
              ENDIF
          ENDDO
      ENDDO     
  ENDIF
   
  CDUM='     OFF'
  READ (CON,*,ERR=400)      CDUM;CDUM=ADJUSTR(CDUM)         
  READ (CON,*,ERR=400)      NDUM         
  READ (CON,*,ERR=400)      NIDUM          
  READ (CON,*,ERR=400)      (DDUM(J),J=1,NDUM)    
  READ (CON,*,ERR=400)      (FDUM(J),J=1,NDUM)    
  READ (CON,*,ERR=400)      (IDUM(J),J=1,NIDUM)   
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)  
  
  SPRC='     OFF'
  NSPR=0
  NISPR=0
  IF(CDUM == '      ON'.OR. CDUM=='     ONV')THEN               
  DO J=1,NDUM
  SPRD(J,1:NWB)=DDUM(J)
  SPRF(J,1:NWB)=FDUM(J)
  ENDDO
  NSPR(1:NWB)=NDUM  
      DO JW=1,NWB 
          JJ=0
          DO J=1,NIDUM 
              IF(IDUM(J) >= US(BS(JW)) .AND. IDUM(J) <= DS(BE(JW)))THEN
                  JJ=JJ+1
                  ISPR(JJ,JW)=IDUM(J)
                  NISPR(JW)=JJ
                  SPRC(JW)=CDUM
              ENDIF
          ENDDO
      ENDDO     
  ENDIF

  VPLC='     OFF'
  NVPL=0
  READ (CON,*,ERR=400)    VPLC(1);VPLC(1)=ADJUSTR(VPLC(1))
  READ (CON,*,ERR=400)    NVPL(1)
  READ (CON,*,ERR=400)   (VPLD(J,1), J=1,NVPL(1))
  READ (CON,*,ERR=400)   (VPLF(J,1), J=1,NVPL(1))
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)  
    WRITE (WIN,*) '  CPL output CPLC'

  READ (CON,*,ERR=400)      CPLC(1);CPLC(1)=ADJUSTR(CPLC(1))
  READ (CON,*,ERR=400)      NCPL(1)
  READ (CON,*,ERR=400)      TCPL(1);TCPL(1)=ADJUSTR(TCPL(1))
  
  CPLC(2:NWB)=CPLC(1)
  NCPL(2:NWB)=NCPL(1)
  TCPL(2:NWB)=TCPL(1)
  READ (CON,*,ERR=400)    (CPLD(J,1), J=1,NCPL(1))
  READ (CON,*,ERR=400)    (CPLF(J,1), J=1,NCPL(1))
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400) 
  DO J=1,NCPL(1)
    CPLD(J,2:NWB)=CPLD(J,1)
    CPLF(J,2:NWB)=CPLF(J,1)
  ENDDO    
    
  READ (CON,*,ERR=400)        FLXC(1);FLXC(1)=ADJUSTR(FLXC(1))
  READ (CON,*,ERR=400)        NFLX(1)
  FLXC(2:NWB)=FLXC(1)
  NFLX(2:NWB)=NFLX(1)

  READ (CON,*,ERR=400)       (FLXD(J,1), J=1,NFLX(1))
  READ (CON,*,ERR=400)       (FLXF(J,1), J=1,NFLX(1))
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400) 
  DO J=1,NFLX(1)
  FLXD(J,2:NWB)=FLXD(J,1)
  FLXF(J,2:NWB)=FLXF(J,1) 
  ENDDO
    WRITE (WIN,*) '  time series output'

  READ (CON,*,ERR=400)           TSRC;TSRC=ADJUSTR(TSRC)
  READ (CON,*,ERR=400)           NTSR
  READ (CON,*,ERR=400)          NIKTSR
  READ (CON,*,ERR=400)          TSRFN
  
  READ (CON,*,ERR=400)        (TSRD(J), J=1,NTSR)
  READ (CON,*,ERR=400)        (TSRF(J), J=1,NTSR)
  READ (CON,*,ERR=400)        (ITSR(J), J=1,NIKTSR)
  READ (CON,*,ERR=400)        (ETSR(J), J=1,NIKTSR)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)AID
          IF (AID /= 'WLEVEL') THEN
          AID1='WLEVEL'
          ICST=3;GO TO 400
          ENDIF  
  READ (CON,*,ERR=400)        WLC;WLC=ADJUSTR(WLC)
  READ (CON,*,ERR=400)        WLF
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)AID
          IF (AID /= 'FLOWBAL') THEN
          AID1='FLOWBAL'
          ICST=3;GO TO 400
          ENDIF  
  READ (CON,*,ERR=400)        FLOWBALC;FLOWBALC=ADJUSTR(FLOWBALC)
  READ (CON,*,ERR=400)        FLOWBALF
    READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)AID
          IF (AID /= 'NPBAL') THEN
          AID1='NPBAL'
          ICST=3;GO TO 400
          ENDIF  
  READ (CON,*,ERR=400)        NPBALC;NPBALC=ADJUSTR(NPBALC)
  READ (CON,*,ERR=400)        NPBALF

  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400) 
     WRITE (WIN,*) '  withdrawal output'

  READ (CON,*,ERR=400)           WDOC; WDOC=ADJUSTR(WDOC)   
  READ (CON,*,ERR=400)           NWDO 
  READ (CON,*,ERR=400)           NIWDO
  READ (CON,*,ERR=400)           WDOFN

!  ALLOCATE (IWDO(MAX(1,NIWDO)))
  READ (CON,*,ERR=400)       (WDOD(J), J=1,NWDO)
  READ (CON,*,ERR=400)       (WDOF(J), J=1,NWDO)
  READ (CON,*,ERR=400)       (IWDO(J), J=1,NIWDO)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)AID
          IF (AID /= 'RESTART') THEN
          AID1='RESTART'
          ICST=3;GO TO 400
          ENDIF  
  WRITE (WIN,*) '  restart RSOC'

  READ (CON,*,ERR=400)         RSOC;RSOC=ADJUSTR(RSOC)
  READ (CON,*,ERR=400)         NRSO
  READ (CON,*,ERR=400)         RSIC;RSIC=ADJUSTR(RSIC)
  READ (CON,*,ERR=400)         RSIFN

  RSOD=0.0 ! SW 9/27/13 INITIALIZE SINCE ALLOCATED AS NOD BUT ONLY NRSO USED
  READ (CON,*,ERR=400)           (RSOD(J), J=1,NRSO)
  READ (CON,*,ERR=400)           (RSOF(J), J=1,NRSO)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)AID
          IF (AID /= 'CCC     ') THEN
          AID1='CCC     '
          ICST=3;GO TO 400
          ENDIF    
ENDIF
    do j=1,nwb
  if(nsnp(j) > nod)then
       CALL ERRORS
       WRITE (ERR,'(A,i3,A,i3,a,i3)') 'The # of output dates for Snapshot output file NSNP[',nsnp(j),'] is > NDAY[',nod,'] for waterbody:',j
       WRITE (WIN,*) 'Error - snapshot output - see pre.err'
       stop
  endif
  if(nisnp(j) > imx)then
       CALL ERRORS
       WRITE (ERR,'(A,i3,A,i3,a,i3)') 'The # segments for SNP output file NISNP[',nisnp(j),'] is > IMX[',imx,'], the max # of segments for waterbody:',j
       WRITE (WIN,*) 'Error - snapshot output - see pre.err'
       stop
  endif
  enddo
do j=1,nwb
  if(nprf(j) > nod)then
       CALL ERRORS
       WRITE (ERR,'(A,i3,A,i3,a,i3)') 'The # of output dates for Profile [PRF] output file NPRF[',nprf(j),'] is > NDAY[',nod,'] for waterbody:',j
       WRITE (WIN,*) 'Error - profile output - see pre.err'
       stop
  endif
  if(niprf(j) > imx)then
       CALL ERRORS
       WRITE (ERR,'(A,i3,A,i3,a,i3)') 'The # segments for Profile [PRF] output file NIPRF[',niprf(j),'] is > IMX[',imx,'], the max # of segments for waterbody:',j
       WRITE (WIN,*) 'Error - profile output - see pre.err'
       stop
  endif
  enddo
  do j=1,nwb
  if(NSPR(j) > nod)then
       CALL ERRORS
       WRITE (ERR,'(A,i3,A,i3,a,i3)') 'The # of output dates for Spreadsheet [SPR] output file NSPR[',NSPR(j),'] is > NDAY[',nod,'] for waterbody:',j
       WRITE (WIN,*) 'Error - spreadsheet output - see pre.err'
       stop
  endif
  if(nispr(j) > imx)then
       CALL ERRORS
       WRITE (ERR,'(A,i3,A,i3,a,i3)') 'The # segments for Spreadsheet [SPR] output file NISPR[',niSPR(j),'] is > IMX[',imx,'], the max # of segments for waterbody:',j
       WRITE (WIN,*) 'Error - spreadsheet output - see pre.err'
       stop
  endif
  enddo
  if(NTSR > nod)then
       CALL ERRORS
       WRITE (ERR,'(A,i3,A,i3,a,i3)') 'The # of output dates for Time Series [TSR] output file NTSR[',NTSR,'] is > NDAY[',nod,'] for waterbody:',j
       WRITE (WIN,*) 'Error - time series output - see pre.err'
       stop
  endif
  if(NIKTSR > imx*kmx)then
       CALL ERRORS
       WRITE (ERR,'(A,i7,A,i7,a,i3)') 'The # segments for Time series [TSR] output file NIKTSR[',NIKTSR,'] is > IMX*KMX[',imx*kmx,'], for waterbody:',j
       WRITE (WIN,*) 'Error - time series output - see pre.err'
       stop
  endif

IF(CONFN=='w2_con.npt')THEN

! Constituent control cards

  WRITE (WIN,*) '  constituent control cards'
  READ (CON,'(/A8/8X,2A8,I8,F8.0,A8)',  ERR=400)  AID, CCC, LIMC, CUF,CO2ATMPPM,CO2ATM_REGRESSION
  IF (AID /= 'CST COMP')                        GO TO 400
  READ (CON,'(/A8/(:2A8))',           ERR=400)  AID, (CNAME2(JC),CAC(JC), JC=1,NCT)
  IF (AID /= 'CST ACTI')                        GO TO 400
  READ (CON,'(/A8/(:8X,9A8))',        ERR=400)  AID, (CDWBC(1,JW),  JW=1,NWB)
                 IF (AID /= 'CST DERI')THEN
                    ICST=1
                    GO TO 400
                 ENDIF
  DO JD=2,NDC
    READ (CON,'(:8X,9A8)',            ERR=400) (CDWBC(JD,JW),       JW=1,NWB)
  END DO
  READ (CON,'(/A8/(:8X,9A8))',        ERR=400)  AID, (CFWBC(1,JW),  JW=1,NWB)
  IF (AID /= 'CST FLUX')                        GO TO 400
  DO JF=2,NFL
    READ (CON,'(:8X,9A8)',            ERR=400) (CFWBC(JF,JW),       JW=1,NWB)
  END DO
  READ (CON,'(/A8/(:8X,9F8.0))',      ERR=400)  AID, (C2IWB(1,JW),  JW=1,NWB)
  IF (AID /= 'CST ICON')                        GO TO 400
  DO JC=2,NCT
    READ (CON,'(:8X,9F8.0)',          ERR=400) (C2IWB(JC,JW),       JW=1,NWB)
  END DO
  READ (CON,'(/A8/(:8X,9A8))',        ERR=400)  AID, (CPRWBC(1,JW), JW=1,NWB)
  IF (AID /= 'CST PRIN')                        GO TO 400
  DO JC=2,NCT
    READ (CON,'(:8X,9A8)',            ERR=400) (CPRWBC(JC,JW),      JW=1,NWB)
  END DO
  READ (CON,'(/A8/(:8X,9A8))',        ERR=400)  AID, (CINBRC(1,JB), JB=1,NBR)
  IF (AID /= 'CIN CON ')                        GO TO 400
  DO JC=2,NCT
    READ (CON,'(:8X,9A8)',            ERR=400) (CINBRC(JC,JB),      JB=1,NBR)
  END DO
  READ (CON,'(/A8/(:8X,9A8))',        ERR=400)  AID, (CINTRC(1,JT), JT=1,NTR)
  IF (AID /= 'CTR CON ')                        GO TO 400
  DO JC=2,NCT
    READ (CON,'(:8X,9A8)',            ERR=400) (CINTRC(JC,JT),      JT=1,NTR)
  END DO
  READ (CON,'(/A8/(:8X,9A8))',        ERR=400)  AID, (CDTBRC(1,JB), JB=1,NBR)
  IF (AID /= 'CDT CON ')                        GO TO 400
  DO JC=2,NCT
    READ (CON,'(:8X,9A8)',            ERR=400) (CDTBRC(JC,JB),      JB=1,NBR)
  END DO
  READ (CON,'(/A8/(:8X,9A8))',        ERR=400)  AID, (CPRBRC(1,JB), JB=1,NBR)
  IF (AID /= 'CPR CON ')                        GO TO 400
  DO JC=2,NCT
    READ (CON,'(:8X,9A8)',            ERR=400) (CPRBRC(JC,JB),      JB=1,NBR)
  END DO

! Kinetics coefficients

  WRITE (WIN,*) '  kinetic coefficient cards'
  READ (CON,'(/A8/(8X,4F8.0,2A8))',     ERR=400)  AID, (EXH2O(JW),  EXSS(JW),   EXOM(JW),   BETA(JW),   EXC(JW),    EXIC(JW),      &
                                                        JW=1,NWB)
  IF (AID /= 'EX COEF ')                          GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',        ERR=400)  AID, (EXA(JA),    JA=1,NAL)
  IF (AID /= 'ALG EX  ')                          GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',        ERR=400)  AID, (EXZ(Jz),    Jz=1,nzpt)
  IF (AID /= 'ZOO EX  ')                          GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',        ERR=400)  AID, (EXM(Jm),    Jm=1,nmct)
  IF (AID /= 'MACRO EX')                          GO TO 400
    WRITE (WIN,*) '      extinction cards'
  READ (CON,'(/A8/(8X,7F8.0))',         ERR=400)  AID, (CGQ10(JG),  CG0DK(JG),  CG1DK(JG),  CGS(JG),CGLDK(JG),CGKLF(JG),CGCS(JG),CGR(JG),    JG=1,NGC)
  IF (AID /= 'GENERIC ')                          GO TO 400
      WRITE (WIN,*) '      generic cards'
  READ (CON,'(/A8/(:8X,F8.0,A8,F8.0))',   ERR=400)  AID, (SSS(JS),    SEDRC(JS),   TAUCR(JS), SSCS(JS), JS=1,NSS)                  ! SW 1/16/04
  IF (AID /= 'S SOLIDS')                          GO TO 400
      WRITE (WIN,*) '      iss cards'
      
  AID='BACTERIA'
    WRITE (WIN,*) '  bacteria'
  READ (CON,'(//(8X,4F8.0))',         ERR=400)         (BACTQ10(JW), BACT1DK(JW),BACTS(JW),BACTLDK(JW),       JW=1,NWB)
  AID='H2S'
    WRITE (WIN,*) '  H2S'
  READ (CON,'(//(8X,4F8.0))',         ERR=400)         (H2SR(JW), H2SQ10(JW), H2S1DK(JW), SO4R(JW),              JW=1,NWB)
  AID='CH4'
    WRITE (WIN,*) '  CH4'
  READ (CON,'(//(8X,3F8.0))',         ERR=400)         (CH4R(JW), CH4Q10(JW), CH41DK(JW),                        JW=1,NWB)
  AID='FE'
    WRITE (WIN,*) '  Fe'
  READ (CON,'(//(8X,5F8.0))',         ERR=400)         (FEIIR(JW), KFE_OXID(JW), KFE_RED(JW), KFEOOH_HalfSat(JW), FeSetVel(JW), JW=1,NWB)
  AID='MN'
    WRITE (WIN,*) '  Mn'
  READ (CON,'(//(8X,5F8.0))',         ERR=400)         (MNIIR(JW), KMN_OXID(JW), KMN_RED(JW), KMNO2_HalfSat(JW),  MnSetVel(JW), JW=1,NWB)
      
  READ (CON,'(/A8/(8X,9F8.0))',         ERR=400)  AID, (AG(JA),     AR(JA),     AE(JA),     AM(JA),     AS(JA),     AHSP(JA),      &
                                                        AHSN(JA),   AHSSI(JA),  ASAT(JA),   JA=1,NAL)                  !TC 09/01/01
  IF (AID /= 'ALGAL RA')                          GO TO 400
  READ (CON,'(/A8/(8X,8F8.0))',         ERR=400)  AID, (AT1(JA),    AT2(JA),    AT3(JA),    AT4(JA),    AK1(JA),    AK2(JA),       &
                                                        AK3(JA),    AK4(JA),    JA=1,NAL)
  IF (AID /= 'ALGAL TE')                          GO TO 400
  READ (CON,'(/A8/(8X,6F8.0,I8,F8.0,A8))', ERR=400)  AID, (ALGP(JA),   ALGN(JA),   ALGC(JA),   ALGSI(JA),  ACHLA(JA),  APOM(JA),      &
                                                        ANEQN(JA),  ANPR(JA), AVERTM(JA),  JA=1,NAL)
  IF (AID /= 'ALG STOI')                          GO TO 400
      WRITE (WIN,*) '      algae cards'
  READ (CON,'(/A8/(:8X,9A8))',          ERR=400)  AID, (EPIC(JW,1),                                              JW=1,NWB)
  IF (AID /= 'EPIPHYTE')                          GO TO 400
  DO JE=2,nept
    READ (CON,'(8X,9A8)',               ERR=400) (EPIC(JW,JE),                                                   JW=1,NWB)
  END DO
  READ (CON,'(/A8/(:8X,9A8))',          ERR=400)  AID, (EPIPRC(JW,1),                                            JW=1,NWB)
  IF (AID /= 'EPI PRIN')                          GO TO 400
  DO JE=2,nept
    READ (CON,'(8X,9A8)',               ERR=400) (EPIPRC(JW,JE),                                                 JW=1,NWB)
  END DO
  READ (CON,'(/A8/(:8X,9F8.0))',        ERR=400)  AID, (EPICI(JW,1),                                             JW=1,NWB)
  IF (AID /= 'EPI INIT')                          GO TO 400
  DO JE=2,nept
    READ (CON,'(8X,9F8.0)',             ERR=400) (EPICI(JW,JE),                                                  JW=1,NWB)
  END DO
  READ (CON,'(/A8/(8X,8F8.0))',         ERR=400)  AID, (EG(JE),     ER(JE),     EE(JE),     EM(JE),     EB(JE),                    &
                                                        EHSP(JE),   EHSN(JE),   EHSSI(JE),                       JE=1,nept)
  IF (AID /= 'EPI RATE')                          GO TO 400
  READ (CON,'(/A8/(8X,2F8.0,I8,F8.0))', ERR=400)  AID, (ESAT(JE),   EHS(JE),    ENEQN(JE),  ENPR(JE),            JE=1,nept)
  IF (AID /= 'EPI HALF')                          GO TO 400
  READ (CON,'(/A8/(8X,8F8.0))',         ERR=400)  AID, (ET1(JE),    ET2(JE),    ET3(JE),    ET4(JE),    EK1(JE),                   &
                                                        EK2(JE),    EK3(JE),    EK4(JE),                         JE=1,nept)
  IF (AID /= 'EPI TEMP')                          GO TO 400
  READ (CON,'(/A8/(8X,6F8.0))',         ERR=400)  AID, (EP(JE),     EN(JE),     EC(JE),     ESI(JE),    ECHLA(JE),                 &
                                                        EPOM(JE),                                                JE=1,nept)
  IF (AID /= 'EPI STOI')                          GO TO 400
    WRITE (WIN,*) '      epiphyton cards'

! v3.5 start
  READ (CON,'(/A8/(8X,9F8.0))',         ERR=400)   AID, (zg(jz),zr(jz),zm(jz),zeff(jz),PREFP(jz),ZOOMIN(jz),ZS2P(jz),ZS(JZ),            Jz=1,nzpt)
  IF (AID /= 'ZOOP RAT')                          GO TO 400
   IF(nzpt /= 0)then                                                         !RA 12/12/2006
    READ (CON,'(/A8/(:8X,8F8.0))',         ERR=400)  AID,     (PREFA(ja,1),                                                           Ja=1,nal)
   ELSE                                                                     !RA 12/12/2006
    READ (CON,'(/A8/(:8X))',         ERR=400)  AID                          !RA 12/12/2006
   END IF
    IF (AID /= 'ZOOP ALG')                          GO TO 400
  do jz=2,nzpt
    READ (CON,'((8X,8F8.0))',         ERR=400)       (PREFA(ja,jz),                                                           Ja=1,nal)
  end do
    READ (CON,'(/A8/(:8X,8F8.0))',         ERR=400)   AID,    (PREFz(jjz,1),                                                          Jjz=1,nzpt)
  IF (AID /= 'ZOOP ZOO')                          GO TO 400
  do jz=2,nzpt
    READ (CON,'((8X,8F8.0))',         ERR=400)       (PREFz(jjz,jz),                                                          Jjz=1,nzpt)
  end do
  READ (CON,'(/A8/(8X,8F8.0))',         ERR=400)  AID,         (zT1(Jz),    zT2(Jz),    zT3(Jz),    zT4(Jz),    zK1(Jz),   zK2(Jz),                         &
                                       zK3(Jz),    zK4(Jz),                                                    Jz=1,nzpt)
  IF (AID /= 'ZOOP TEM')                          GO TO 400
  READ (CON,'(/A8/(8X,3F8.0))',         ERR=400)  AID,        (zP(Jz),     zN(Jz),     zC(Jz),                                         Jz=1,nzpt)      !MLM 7/11/06 3f8.0
  IF (AID /= 'ZOOP STO')                          GO TO 400
      WRITE (WIN,*) '      zooplankton cards'

  READ (CON,'(/A8/(:8X,9A8))',         ERR=400)  AID,      (macwbC(JW,1),                                                           JW=1,NWB)  ! mlm 7/11/06
    IF (AID /= 'MACROPHY')                          GO TO 400
  DO Jm=2,nmct
    READ (CON,'(8X,9A8)',         ERR=400)             (macwbC(JW,Jm),                                                          JW=1,NWB)
  END DO

  READ (CON,'(/A8/(:8X,9A8))',         ERR=400)  AID,           (mprwbC(JW,1),                                                           JW=1,NWB)
     IF (AID /= 'MAC PRIN')                          GO TO 400

  DO Jm=2,nmct
    READ (CON,'(8X,9A8)',         ERR=400)             (mprwbC(JW,Jm),                                                          JW=1,NWB)
  END DO


  READ (CON,'(/A8/(:8X,9F8.0))',         ERR=400)  AID,         (macwbCI(JW,1),                                                          JW=1,NWB)
  IF (AID /= 'MAC INI ')                          GO TO 400

  DO Jm=2,nmct
    READ (CON,'(8X,9F8.0)',         ERR=400)           (macwbcI(JW,Jm),                                                         JW=1,NWB)
  END DO

  READ (CON,'(/A8/(8X,9F8.0))',         ERR=400)  AID,         (mG(jm), mR(jm), mM(jm), msat(jm),mhsp(jm),mhsn(jm),mhsc(jm),                           &
                                          mpom(jm),lrpmac(jm),     jm=1,nmct)
    IF (AID /= 'MAC RATE')                          GO TO 400
  READ (CON,'(/A8/(8X,8F8.0))',         ERR=400)  AID,         (psed(jm), nsed(jm),                                                     jm=1,nmct)
    IF (AID /= 'MAC SED ')                          GO TO 400
  READ (CON,'(/A8/(8X,8F8.0))',         ERR=400)  AID,         (mbmp(jm), mmax(jm),                                                     jm=1,nmct)
    IF (AID /= 'MAC DIST')                          GO TO 400
  READ (CON,'(/A8/(8X,8F8.0))',         ERR=400)  AID,         (cdstem(jm), dmv(jm),dwsa(jm),anorm(jm),                                 jm=1,nmct)
    IF (AID /= 'MAC DRAG')                          GO TO 400
  READ (CON,'(/A8/(8X,8F8.0))',         ERR=400)  AID,         (mT1(Jm),    mT2(Jm),    mT3(Jm),    mT4(Jm),    mK1(Jm),   mK2(Jm),                         &
                                       mK3(Jm),    mK4(Jm),                                                    Jm=1,nmct)
    IF (AID /= 'MAC TEMP')                          GO TO 400
  READ (CON,'(/A8/(8X,6F8.0))',         ERR=400)  AID,         (mP(Jm),     mN(Jm),     mC(Jm),                                         Jm=1,nmct)
    IF (AID /= 'MAC STOI')                          GO TO 400
    WRITE (WIN,*) '      macrophyte cards'

  READ (CON,'(/A8/(8X,3F8.0))',         ERR=400)  AID, (LDOMDK(JW), RDOMDK(JW), LRDDK(JW),  JW=1,NWB)
  IF (AID /= 'DOM     ')                          GO TO 400
  READ (CON,'(/A8/(8X,4F8.0))',         ERR=400)  AID, (LPOMDK(JW), RPOMDK(JW), LRPDK(JW),  POMS(JW),   JW=1,NWB)
  IF (AID /= 'POM     ')                          GO TO 400
  READ (CON,'(/A8/(8X,4F8.0))',         ERR=400)  AID, (ORGP(JW),   ORGN(JW),   ORGC(JW),   ORGSI(JW), CoeffA_Turb(JW), CoeffB_Turb(JW),SECC_PAR(JW), JW=1,NWB)
  IF (AID /= 'OM STOIC')                          GO TO 400
  READ (CON,'(/A8/(8X,4F8.0))',         ERR=400)  AID, (OMT1(JW),   OMT2(JW),   OMK1(JW),   OMK2(JW),   JW=1,NWB)
  IF (AID /= 'OM RATE ')                          GO TO 400
  READ (CON,'(/A8/(8X,4F8.0))',         ERR=400)  AID, (KBOD(JB),   TBOD(JB),   RBOD(JB), cbods(jb),  JB=1,NBOD)
  IF (AID /= 'CBOD    ')                          GO TO 400
  READ (CON,'(/A8/(8X,3F8.0))',         ERR=400)  AID, (BODP(JB),   BODN(JB),   BODC(JB),   JB=1,NBOD)                   !TC 01/15/02
  IF (AID /= 'CBOD STO')                          GO TO 400
      WRITE (WIN,*) '      OM-BOD cards'

  READ (CON,'(/A8/(8X,2F8.0))',         ERR=400)  AID, (PO4R(JW),   PARTP(JW),  JW=1,NWB)
  IF (AID /= 'PHOSPHOR')                          GO TO 400
  READ (CON,'(/A8/(8X,2F8.0))',         ERR=400)  AID, (NH4R(JW),   NH4DK(JW),  JW=1,NWB)
  IF (AID /= 'AMMONIUM')                          GO TO 400
  READ (CON,'(/A8/(8X,4F8.0))',         ERR=400)  AID, (NH4T1(JW),  NH4T2(JW),  NH4K1(JW),  NH4K2(JW),  JW=1,NWB)
  IF (AID /= 'NH4 RATE')                          GO TO 400
  READ (CON,'(/A8/(8X,3F8.0))',          ERR=400)  AID, (NO3DK(JW),  NO3S(JW),FNO3SED(JW), JW=1,NWB)
  IF (AID /= 'NITRATE ')                          GO TO 400
  READ (CON,'(/A8/(8X,4F8.0))',         ERR=400)  AID, (NO3T1(JW),  NO3T2(JW),  NO3K1(JW),  NO3K2(JW),  JW=1,NWB)
  IF (AID /= 'NO3 RATE')                          GO TO 400
  READ (CON,'(/A8/(8X,4F8.0))',         ERR=400)  AID, (DSIR(JW),   PSIS(JW),   PSIDK(JW),  PARTSI(JW), JW=1,NWB)
  IF (AID /= 'SILICA  ')                          GO TO 400
      WRITE (WIN,*) '      N-P-Si cards'

  READ (CON,'(/A8/(8X,F8.0))',          ERR=400)  AID, (CO2R(JW),   JW=1,NWB)
  IF (AID /= 'SED CO2 ')                          GO TO 400
  READ (CON,'(/A8/(8X,2F8.0))',         ERR=400)  AID, (O2NH4(JW),  O2OM(JW),   JW=1,NWB)
  IF (AID /= 'STOICH 1')                          GO TO 400
  READ (CON,'(/A8/(8X,2F8.0))',         ERR=400)  AID, (O2AR(JA),   O2AG(JA),   JA=1,NAL)
  IF (AID /= 'STOICH 2')                          GO TO 400
  READ (CON,'(/A8/(8X,2F8.0))',         ERR=400)  AID, (O2ER(JE),   O2EG(JE),   JE=1,nept)
  IF (AID /= 'STOICH 3')                          GO TO 400
  READ (CON,'(/A8/(8X,1F8.0))',         ERR=400)  AID, (O2ZR(JZ),   JZ=1,nzpt)
  IF (AID /= 'STOICH 4')                          GO TO 400
  READ (CON,'(/A8/(8X,2F8.0))',         ERR=400)  AID, (O2MR(Jm),   O2MG(Jm),   Jm=1,nmct)
  IF (AID /= 'STOICH 5')                          GO TO 400

  READ (CON,'(/A8/(8X,F8.0))',          ERR=400)  AID, O2LIM
  IF (AID /= 'O2 LIMIT')                          GO TO 400
      WRITE (WIN,*) '      O2 cards'

  READ (CON,'(/A8/(8X,2A8,6F8.0,A8))',     ERR=400)  AID, (SEDC(JW),   PRNSC(JW),  SEDCI(JW),  SEDDK(JW),seds(jw),  FSOD(JW),                  &     ! CB 12/30/12
                                                        FSED(JW), sedbr(jw), DYNSEDK(JW), JW=1,NWB)                                         ! SW 6/1/07
  IF (AID /= 'SEDIMENT')                          GO TO 400
  READ (CON,'(/A8/(8X,4F8.0))',         ERR=400)  AID, (SODT1(JW),  SODT2(JW),  SODK1(JW),  SODK2(JW),  JW=1,NWB)
  IF (AID /= 'SOD RATE')                          GO TO 400
  READ (CON,'(/A8/(:8X,9F8.0))',        ERR=400)  AID, (SOD(I),     I=1,IMX)
  IF (AID /= 'S DEMAND')                          GO TO 400
      WRITE (WIN,*) '      sediment cards'

  READ (CON,'(/A8/(8X,A8,I8,4F8.2))',   ERR=400)  AID, (REAERC(JW), NEQN(JW),   RCOEF1(JW), RCOEF2(JW), RCOEF3(JW), RCOEF4(JW), DGPO2(JW), MINKL(JW),  &
                                                        JW=1,NWB)
  IF (AID /= 'REAERATI')                          GO TO 400

! Input filenames

  WRITE (WIN,*) '  input filename cards'
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, RSIFN
  IF (AID /= 'RSI FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, QWDFN
  IF (AID /= 'QWD FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, QGTFN
  IF (AID /= 'QGT FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, WSCFN                                                                       !TC 03/12/02                                                                                    !elo
  IF (AID /= 'WSC FILE')              GO TO 400                                                                        !TC 03/12/02
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, SHDFN                                                                       !SW 04/03/02                                                                            !elo
  IF (AID /= 'SHD FILE')              GO TO 400                                                                        !SW 04/03/02
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (BTHFN(JW),JW=1,NWB)
  IF (AID /= 'BTH FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (METFN(JW),JW=1,NWB)
  IF (AID /= 'MET FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (EXTFN(JW), JW=1,NWB)                                                       !SW 12/12/01
  IF (AID /= 'EXT FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (ATMDEPFN(JW), JW=1,NWB)                                                       !SW 12/12/01
  IF (AID /= 'ATD FILE')              GO TO 400
  
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (VPRFN(JW),JW=1,NWB)
  IF (AID /= 'VPR FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (LPRFN(JW),JW=1,NWB)
  IF (AID /= 'LPR FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (QINFN(JB),JB=1,NBR)
  IF (AID /= 'QIN FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (TINFN(JB),JB=1,NBR)
  IF (AID /= 'TIN FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (CINFN(JB),JB=1,NBR)
  IF (AID /= 'CIN FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (QOTFN(JB),JB=1,NBR)
  IF (AID /= 'QOT FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (QTRFN(JT),JT=1,NTR)
  IF (AID /= 'QTR FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (TTRFN(JT),JT=1,NTR)
  IF (AID /= 'TTR FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (CTRFN(JT),JT=1,NTR)
  IF (AID /= 'CTR FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (QDTFN(JB),JB=1,NBR)
  IF (AID /= 'QDT FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (TDTFN(JB),JB=1,NBR)
  IF (AID /= 'TDT FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (CDTFN(JB),JB=1,NBR)
  IF (AID /= 'CDT FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (PREFN(JB),JB=1,NBR)
  IF (AID /= 'PRE FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (TPRFN(JB),JB=1,NBR)
  IF (AID /= 'TPR FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (CPRFN(JB),JB=1,NBR)
  IF (AID /= 'CPR FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (EUHFN(JB),JB=1,NBR)
  IF (AID /= 'EUH FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (TUHFN(JB),JB=1,NBR)
  IF (AID /= 'TUH FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (CUHFN(JB),JB=1,NBR)
  IF (AID /= 'CUH FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (EDHFN(JB),JB=1,NBR)
  IF (AID /= 'EDH FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (TDHFN(JB),JB=1,NBR)
  IF (AID /= 'TDH FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (CDHFN(JB),JB=1,NBR)
  IF (AID /= 'CDH FILE')              GO TO 400

! Output filenames

  WRITE (WIN,*) '  output filename cards'
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (SNPFN(JW),JW=1,NWB)
  IF (AID /= 'SNP FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (PRFFN(JW),JW=1,NWB)
  IF (AID /= 'PRF FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (VPLFN(JW),JW=1,NWB)
  IF (AID /= 'VPL FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (CPLFN(JW),JW=1,NWB)
  IF (AID /= 'CPL FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (SPRFN(JW),JW=1,NWB)
  IF (AID /= 'SPR FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, (FLXFN(JW),JW=1,NWB)
  IF (AID /= 'FLX FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, TSRFN
  IF (AID /= 'TSR FILE')              GO TO 400
  READ (CON,'(/A8/(8X,A72))',ERR=400) AID, WDOFN
  IF (AID /= 'WDO FILE')              GO TO 400
  CLOSE (CON)
ELSE
    
    ! Constituent control cards
  WRITE (WIN,*) '  constituent control CCC'

  READ (CON,*,ERR=400)           CCC, LIMC, CUF,CO2ATMPPM,CO2ATM_REGRESSION;CCC=ADJUSTR(CCC);LIMC=ADJUSTR(LIMC);CO2ATM_REGRESSION=ADJUSTR(CO2ATM_REGRESSION) 
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400) 
  
  READ (CON,*,ERR=400)           (ATM_DEPOSITIONC(JW), JW=1,NWB);ATM_DEPOSITIONC=ADJUSTR(ATM_DEPOSITIONC)
  READ (CON,*,ERR=400)           (ATM_DEPOSITION_INTERPOLATION(JW), JW=1,NWB);ATM_DEPOSITION_INTERPOLATION=ADJUSTR(ATM_DEPOSITION_INTERPOLATION)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400) 
  
  DO JC=1,NCT
  READ (CON,*,ERR=400)CNAME2(JC),CNAME(JC),CAC(JC), FMTC(JC), CMULT(JC), (C2IWB(JC,JW), JW=1,NWB),(CPRWBC(JC,JW), JW=1,NWB),(C_ATM_DEPOSITION(JC,JW), JW=1,NWB), (CINBRC(JC,JB), JB=1,NBR),(CINTRC(JC,JT), JT=1,NTR1), (CDTBRC(JC,JB), JB=1,NBR), (CPRBRC(JC,JB), JB=1,NBR)  
  ENDDO
  CAC=ADJUSTR(CAC);CPRWBC=ADJUSTR(CPRWBC);CINBRC=ADJUSTR(CINBRC);CINTRC=ADJUSTR(CINTRC);CDTBRC=ADJUSTR(CDTBRC);CPRBRC=ADJUSTR(CPRBRC);C_ATM_DEPOSITION=ADJUSTR(C_ATM_DEPOSITION)
  READ (CON,*,ERR=400)
  WRITE (WIN,*) '  derived concentrations'
  READ (CON,*,ERR=400) AID
          IF (AID /= 'CDNAME2 ') THEN
          AID1='CDNAME2 '
          GO TO 400
      ENDIF

  DO JD=1,NDC
  READ (CON,*,ERR=400)  CDNAME2,CDNAME(JD),FMTCD(JD),CDMULT(JD),(CDWBC(JD,JW), JW=1,NWB)
  ENDDO
  CDWBC=ADJUSTR(CDWBC)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400) 

 WRITE (WIN,*) '  constituent fluxes'

  KFWBC   ='     '   ! SW 9/27/13 INITIALIZE ENTIRE ARRAY
  DO JF=1,NFL   ! Fix this later
    READ (CON,*,ERR=400)   KFNAME2,(CFWBC(JF,JW),  JW=1,NWB)
  END DO
  CFWBC=ADJUSTR(CFWBC)

! Kinetics coefficients
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400) 
  WRITE (WIN,*) '  extinction EXH2O'

  READ (CON,*,ERR=400)     (EXH2O(JW),  JW=1,NWB)
  READ (CON,*,ERR=400)     (EXSS(JW),   JW=1,NWB)
  READ (CON,*,ERR=400)     (EXOM(JW),   JW=1,NWB)
  READ (CON,*,ERR=400)     (BETA(JW),   JW=1,NWB)
  READ (CON,*,ERR=400)     (EXC(JW),    JW=1,NWB);EXC=ADJUSTR(EXC)
  READ (CON,*,ERR=400)     (EXIC(JW),   JW=1,NWB);EXIC=ADJUSTR(EXIC)
 
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400) 

  READ (CON,*,ERR=400)         (EXA(JA),  JA=1,NAL)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400) 

  READ (CON,*,ERR=400)         (EXZ(JZ),  JZ=1,NZPT)  
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400) 
  
  READ (CON,*,ERR=400)         (EXM(JM),   JM=1,NMCT)  
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400) 
  WRITE (WIN,*) '  generic constituents'

  READ (CON,*,ERR=400)         (CGQ10(JG),  JG=1,NGC)
  READ (CON,*,ERR=400)         (CG0DK(JG),  JG=1,NGC)
  READ (CON,*,ERR=400)         (CG1DK(JG),  JG=1,NGC) 
  READ (CON,*,ERR=400)         (CGS(JG),    JG=1,NGC) 
  READ (CON,*,ERR=400)         (CGLDK(JG),  JG=1,NGC) 
  READ (CON,*,ERR=400)         (CGKLF(JG),  JG=1,NGC) 
  READ (CON,*,ERR=400)         (CGCS(JG),   JG=1,NGC) 
  READ (CON,*,ERR=400)         (CGR(JG),    JG=1,NGC)
 
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)  
    WRITE (WIN,*) '  suspended solids'

  READ (CON,*,ERR=400)   (SSS(JS),    JS=1,NSS)     ! READ (CON,'(//(8X,F8.0,A8,2F8.0,I8))') (SSS(JS), SEDRC(JS),  TAUCR(JS),  SSFLOC(JS), FLOCEQN(JS),            JS=1,NSS) !SR 04/21/13
  READ (CON,*,ERR=400)   (SEDRC(JS),  JS=1,NSS)     ! READ (CON,'(//(8X,F8.0,A8,2F8.0,I8))') (SSS(JS), SEDRC(JS),  TAUCR(JS),  SSFLOC(JS), FLOCEQN(JS),            JS=1,NSS) !SR 04/21/13
  READ (CON,*,ERR=400)   (TAUCR(JS),  JS=1,NSS)     ! READ (CON,'(//(8X,F8.0,A8,2F8.0,I8))') (SSS(JS), SEDRC(JS),  TAUCR(JS),  SSFLOC(JS), FLOCEQN(JS),            JS=1,NSS) !SR 04/21/13
  READ (CON,*,ERR=400)   (SSCS(JS),   JS=1,NSS)     ! READ (CON,'(//(8X,F8.0,A8,2F8.0,I8))') (SSS(JS), SEDRC(JS),  TAUCR(JS),  SSFLOC(JS), FLOCEQN(JS),            JS=1,NSS) !SR 04/21/13
  SEDRC=ADJUSTR(SEDRC)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)  
      WRITE (WIN,*) '  bacteria'

AID='BACTERIA'
  READ (CON,*,ERR=400)     (BACTQ10(JW),  JW=1,NWB)
  READ (CON,*,ERR=400)     (BACT1DK(JW),  JW=1,NWB)
  READ (CON,*,ERR=400)     (BACTLDK(JW),  JW=1,NWB)
  READ (CON,*,ERR=400)     (BACTS(JW),    JW=1,NWB)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)
  
  !READ (CON,*)     (A_DISG(JW),  JW=1,NWB)
  !READ (CON,*)     (B_DISG(JW),  JW=1,NWB)
  !READ (CON,*)     (C_DISG(JW),  JW=1,NWB)
  !READ (CON,*)
  !READ (CON,*)
        WRITE (WIN,*) '  H2S'

  AID='H2S'
  READ (CON,*,ERR=400)     (H2SR(JW),    JW=1,NWB)
  READ (CON,*,ERR=400)     (H2SQ10(JW),  JW=1,NWB)
  READ (CON,*,ERR=400)     (H2S1DK(JW),  JW=1,NWB)
  READ (CON,*,ERR=400)     (SO4R(JW),    JW=1,NWB)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)
  WRITE (WIN,*) '  CH4'

  AID='CH4'
  READ (CON,*,ERR=400)     (CH4R(JW),    JW=1,NWB)
  READ (CON,*,ERR=400)     (CH4Q10(JW),  JW=1,NWB)
  READ (CON,*,ERR=400)     (CH41DK(JW),  JW=1,NWB)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)
  WRITE (WIN,*) '  Fe'
  AID='FE'
  READ (CON,*,ERR=400)     (FEIIR(JW),    JW=1,NWB)
  READ (CON,*,ERR=400)     (KFE_OXID(JW),  JW=1,NWB)
  READ (CON,*,ERR=400)     (KFE_RED(JW),  JW=1,NWB)
  READ (CON,*,ERR=400)     (KFEOOH_HalfSat(JW),    JW=1,NWB)
  READ (CON,*,ERR=400)     (FeSetVel(JW),  JW=1,NWB)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)
        WRITE (WIN,*) '  Mn'
  AID='MN'
  READ (CON,*,ERR=400)     (MNIIR(JW),    JW=1,NWB)
  READ (CON,*,ERR=400)     (KMN_OXID(JW),  JW=1,NWB)
  READ (CON,*,ERR=400)     (KMN_RED(JW),  JW=1,NWB)
  READ (CON,*,ERR=400)     (KMNO2_HalfSat(JW),    JW=1,NWB)
  READ (CON,*,ERR=400)     (MNSetVel(JW),  JW=1,NWB)  
  
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)  
    WRITE (WIN,*) '  algae'

  READ (CON,*,ERR=400) (AG(JA),          JA=1,NAL)
  READ (CON,*,ERR=400) (AR(JA),          JA=1,NAL)
  READ (CON,*,ERR=400) (AE(JA),          JA=1,NAL)
  READ (CON,*,ERR=400) (AM(JA),          JA=1,NAL)
  READ (CON,*,ERR=400) (AS(JA),          JA=1,NAL)
  READ (CON,*,ERR=400) (AHSP(JA),        JA=1,NAL)
  READ (CON,*,ERR=400) (AHSN(JA),        JA=1,NAL)
  READ (CON,*,ERR=400) (AHSSI(JA),       JA=1,NAL)
  READ (CON,*,ERR=400) (ASAT(JA),        JA=1,NAL)

  READ (CON,*,ERR=400) (AT1(JA),         JA=1,NAL)
  READ (CON,*,ERR=400) (AT2(JA),         JA=1,NAL)
  READ (CON,*,ERR=400) (AT3(JA),         JA=1,NAL)
  READ (CON,*,ERR=400) (AT4(JA),         JA=1,NAL)
  READ (CON,*,ERR=400) (AK1(JA),         JA=1,NAL)
  READ (CON,*,ERR=400) (AK2(JA),         JA=1,NAL)
  READ (CON,*,ERR=400) (AK3(JA),         JA=1,NAL)
  READ (CON,*,ERR=400) (AK4(JA),         JA=1,NAL)

  READ (CON,*,ERR=400) (ALGP(JA),          JA=1,NAL)
  READ (CON,*,ERR=400) (ALGN(JA),          JA=1,NAL)
  READ (CON,*,ERR=400) (ALGC(JA),          JA=1,NAL)
  READ (CON,*,ERR=400) (ALGSI(JA),         JA=1,NAL)
  READ (CON,*,ERR=400) (ACHLA(JA),       JA=1,NAL)
  READ (CON,*,ERR=400) (APOM(JA),        JA=1,NAL)
  READ (CON,*,ERR=400) (ANEQN(JA),       JA=1,NAL)
  READ (CON,*,ERR=400) (ANPR(JA),        JA=1,NAL)

  READ (CON,*,ERR=400) (O2AR(JA),        JA=1,NAL)
  READ (CON,*,ERR=400) (O2AG(JA),        JA=1,NAL)
  READ (CON,*,ERR=400) (AVERTM(JA),      JA=1,NAL)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)  
      WRITE (WIN,*) '  periphyton'

  DO JE=1,NEPTT
  READ (CON,*,ERR=400)         (EPIC(JW,JE),  JW=1,NWB)
  READ (CON,*,ERR=400)         (EPIPRC(JW,JE),JW=1,NWB)
  READ (CON,*,ERR=400)         (EPICI(JW,JE), JW=1,NWB)
  ENDDO
  EPIC=ADJUSTR(EPIC);EPIPRC=ADJUSTR(EPIPRC)

  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)  
  
  READ (CON,*,ERR=400) (EG(JE),           JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*,ERR=400) (ER(JE),           JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*,ERR=400) (EE(JE),           JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*,ERR=400) (EM(JE),           JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*,ERR=400) (EB(JE),           JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*,ERR=400) (EHSP(JE),         JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*,ERR=400) (EHSN(JE),         JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*,ERR=400) (EHSSI(JE),        JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13  
  
  READ (CON,*,ERR=400) (ESAT(JE),         JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*,ERR=400) (EHS(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*,ERR=400) (ENEQN(JE),        JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*,ERR=400) (ENPR(JE),         JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13  
  
  READ (CON,*,ERR=400) (ET1(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*,ERR=400) (ET2(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*,ERR=400) (ET3(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*,ERR=400) (ET4(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*,ERR=400) (EK1(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*,ERR=400) (EK2(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*,ERR=400) (EK3(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*,ERR=400) (EK4(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  
  READ (CON,*,ERR=400) (EP(JE),           JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*,ERR=400) (EN(JE),           JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*,ERR=400) (EC(JE),           JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*,ERR=400) (ESI(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*,ERR=400) (ECHLA(JE),        JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*,ERR=400) (EPOM(JE),         JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*,ERR=400) (O2ER(JE),         JE=1,NEPT)
  READ (CON,*,ERR=400) (O2EG(JE),         JE=1,NEPT)

  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)  
  WRITE (WIN,*) '  zooplankton'
   
  READ (CON,*,ERR=400)         (ZG(JZ),    JZ=1,NZPT)
  READ (CON,*,ERR=400)         (ZR(JZ),    JZ=1,NZPT)
  READ (CON,*,ERR=400)         (ZM(JZ),    JZ=1,NZPT)
  READ (CON,*,ERR=400)         (ZEFF(JZ),  JZ=1,NZPT)
  READ (CON,*,ERR=400)         (PREFP(JZ), JZ=1,NZPT)
  READ (CON,*,ERR=400)         (ZOOMIN(JZ),JZ=1,NZPT)
  READ (CON,*,ERR=400)         (ZS2P(JZ),  JZ=1,NZPT)
  READ (CON,*,ERR=400)         (ZS(JZ),    JZ=1,NZPT)
  READ (CON,*,ERR=400)         (ZT1(JZ),   JZ=1,NZPT)
  READ (CON,*,ERR=400)         (ZT2(JZ),   JZ=1,NZPT)
  READ (CON,*,ERR=400)         (ZT3(JZ),   JZ=1,NZPT)
  READ (CON,*,ERR=400)         (ZT4(JZ),   JZ=1,NZPT)
  READ (CON,*,ERR=400)         (ZK1(JZ),   JZ=1,NZPT)
  READ (CON,*,ERR=400)         (ZK2(JZ),   JZ=1,NZPT)
  READ (CON,*,ERR=400)         (ZK3(JZ),   JZ=1,NZPT)
  READ (CON,*,ERR=400)         (ZK4(JZ),   JZ=1,NZPT) 
  
  READ (CON,*,ERR=400)         (ZP(JZ),    JZ=1,NZPT)  
  READ (CON,*,ERR=400)         (ZN(JZ),    JZ=1,NZPT)
  READ (CON,*,ERR=400)         (ZC(JZ),    JZ=1,NZPT)

  READ (CON,*,ERR=400)         (O2ZR(JZ),  JZ=1,NZPT)

  DO JA=1,NALT
    READ (CON,*,ERR=400)    (PREFA(JA,JZ),     JZ=1,NZPT)
  END DO

  DO JZ=1,NZPTT
    READ (CON,*,ERR=400)    (PREFZ(JZ,JJZ),   JJZ=1,NZPT)       
  END DO
  
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)   
    WRITE (WIN,*) '  macrophytes'

  DO JM=1,NMCTT
    READ (CON,*,ERR=400)   (MACWBC(JW,JM),  JW=1,NWB)
  END DO
  DO JM=1,NMCTT
    READ (CON,*,ERR=400)   (MPRWBC(JW,JM),  JW=1,NWB)
  END DO
  DO JM=1,NMCTT
    READ (CON,*,ERR=400)   (MACWBCI(JW,JM),  JW=1,NWB)
  END DO
  MACWBC=ADJUSTR(MACWBC);MPRWBC=ADJUSTR(MPRWBC)
    
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)AID
      IF (AID /= 'MacGroup') THEN
          AID1='MacGroup'
          ICST=3;GO TO 400
      ENDIF

  
  READ (CON,*,ERR=400)         (MG(JM),     JM=1,NMCT)
  READ (CON,*,ERR=400)         (MR(JM),     JM=1,NMCT)
  READ (CON,*,ERR=400)         (MM(JM),     JM=1,NMCT)
  READ (CON,*,ERR=400)         (MSAT(JM),   JM=1,NMCT)
  READ (CON,*,ERR=400)         (MHSP(JM),   JM=1,NMCT)
  READ (CON,*,ERR=400)         (MHSN(JM),   JM=1,NMCT)
  READ (CON,*,ERR=400)         (MHSC(JM),   JM=1,NMCT)
  READ (CON,*,ERR=400)         (MPOM(JM),   JM=1,NMCT)
  READ (CON,*,ERR=400)         (LRPMAC(JM), JM=1,NMCT)
  
  READ (CON,*,ERR=400)         (PSED(JM),   JM=1,NMCT)  
  READ (CON,*,ERR=400)         (NSED(JM),   JM=1,NMCT)

  READ (CON,*,ERR=400)         (MBMP(JM),   JM=1,NMCT)
  READ (CON,*,ERR=400)         (MMAX(JM),   JM=1,NMCT)
  READ (CON,*,ERR=400)         (CDSTEM(JM), JM=1,NMCT)  !CB 6/29/06
  READ (CON,*,ERR=400)         (DMV(JM),    JM=1,NMCT)  !CB 6/29/06
  READ (CON,*,ERR=400)         (DWSA(JM),   JM=1,NMCT)  !CB 6/29/06
  READ (CON,*,ERR=400)         (ANORM(JM),  JM=1,NMCT)  !CB 6/29/06  
  
  READ (CON,*,ERR=400)         (MT1(JM),    JM=1,NMCT)
  READ (CON,*,ERR=400)         (MT2(JM),    JM=1,NMCT)
  READ (CON,*,ERR=400)         (MT3(JM),    JM=1,NMCT)
  READ (CON,*,ERR=400)         (MT4(JM),    JM=1,NMCT)
  READ (CON,*,ERR=400)         (MK1(JM),    JM=1,NMCT)
  READ (CON,*,ERR=400)         (MK2(JM),    JM=1,NMCT)
  READ (CON,*,ERR=400)         (MK3(JM),    JM=1,NMCT)
  READ (CON,*,ERR=400)         (MK4(JM),    JM=1,NMCT)
  
  READ (CON,*,ERR=400)         (MP(JM),     JM=1,NMCT)
  READ (CON,*,ERR=400)         (MN(JM),     JM=1,NMCT)
  READ (CON,*,ERR=400)         (MC(JM),     JM=1,NMCT)
 
  READ (CON,*,ERR=400)         (O2MR(JM),   JM=1,NMCT)
  READ (CON,*,ERR=400)         (O2MG(JM),   JM=1,NMCT)

  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)   
    WRITE (WIN,*) '  organic matter LDOMDK'

  READ (CON,*,ERR=400)         (LDOMDK(JW),  JW=1,NWB)
  READ (CON,*,ERR=400)         (RDOMDK(JW),  JW=1,NWB)
  READ (CON,*,ERR=400)         (LRDDK(JW),   JW=1,NWB) 

  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)   
   
  READ (CON,*,ERR=400)         (LPOMDK(JW),   JW=1,NWB)
  READ (CON,*,ERR=400)         (RPOMDK(JW),   JW=1,NWB)
  READ (CON,*,ERR=400)         (LRPDK(JW),    JW=1,NWB)
  READ (CON,*,ERR=400)         (POMS(JW),     JW=1,NWB)
  
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)   
  AID='ORGP'
  WRITE (WIN,*) '  organic matter stoich'

  READ (CON,*,ERR=400)         (ORGP(JW),     JW=1,NWB)
  READ (CON,*,ERR=400)         (ORGN(JW),     JW=1,NWB)
  READ (CON,*,ERR=400)         (ORGC(JW),     JW=1,NWB)
  READ (CON,*,ERR=400)         (ORGSI(JW),    JW=1,NWB)
  READ (CON,*,ERR=400)         (O2OM(JW),     JW=1,NWB)

  READ (CON,*,ERR=400)         (OMT1(JW),     JW=1,NWB)
  READ (CON,*,ERR=400)         (OMT2(JW),     JW=1,NWB)
  READ (CON,*,ERR=400)         (OMK1(JW),     JW=1,NWB)
  READ (CON,*,ERR=400)         (OMK2(JW),     JW=1,NWB)
  
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)   
  AID='Turbidity'
  WRITE (WIN,*) '  turbidity'
  READ (CON,*,ERR=400)         (CoeffA_Turb(JW),    JW=1,NWB)
  READ (CON,*,ERR=400)         (CoeffB_Turb(JW),    JW=1,NWB)
  READ (CON,*,ERR=400)         (SECC_PAR(JW),    JW=1,NWB)

  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)   
    WRITE (WIN,*) '  BOD'
  AID='BOD'

  READ (CON,*,ERR=400)         (KBOD(JB),     JB=1,NBOD)
  READ (CON,*,ERR=400)         (TBOD(JB),     JB=1,NBOD)
  READ (CON,*,ERR=400)         (RBOD(JB),     JB=1,NBOD)
  READ (CON,*,ERR=400)         (CBODS(JB),    JB=1,NBOD)
  READ (CON,*,ERR=400)         (BODP(JB),     JB=1,NBOD)
  READ (CON,*,ERR=400)         (BODN(JB),     JB=1,NBOD)
  READ (CON,*,ERR=400)         (BODC(JB),     JB=1,NBOD)
  
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)   
  WRITE (WIN,*) '  P, N, and Si'
  AID='P,N,AND Si'
  READ (CON,*,ERR=400)         (PO4R(JW),     JW=1,NWB)
  READ (CON,*,ERR=400)         (PARTP(JW),    JW=1,NWB)
  READ (CON,*,ERR=400)         (NH4R(JW),     JW=1,NWB)
  READ (CON,*,ERR=400)         (NH4DK(JW),    JW=1,NWB)
  READ (CON,*,ERR=400)         (NH4T1(JW),    JW=1,NWB)
  READ (CON,*,ERR=400)         (NH4T2(JW),    JW=1,NWB)
  READ (CON,*,ERR=400)         (NH4K1(JW),    JW=1,NWB)
  READ (CON,*,ERR=400)         (NH4K2(JW),    JW=1,NWB)
  READ (CON,*,ERR=400)         (KG_H2O_CONSTANT(JW), JW=1,NWB)
  READ (CON,*,ERR=400)         (O2NH4(JW),    JW=1,NWB)
  READ (CON,*,ERR=400)         (NO3DK(JW),    JW=1,NWB)
  READ (CON,*,ERR=400)         (NO3S(JW),     JW=1,NWB)
  READ (CON,*,ERR=400)         (FNO3SED(JW),  JW=1,NWB)
  READ (CON,*,ERR=400)         (NO3T1(JW),    JW=1,NWB)
  READ (CON,*,ERR=400)         (NO3T2(JW),    JW=1,NWB)
  READ (CON,*,ERR=400)         (NO3K1(JW),    JW=1,NWB)
  READ (CON,*,ERR=400)         (NO3K2(JW),    JW=1,NWB)
  READ (CON,*,ERR=400)         (DSIR(JW),     JW=1,NWB)
  READ (CON,*,ERR=400)         (PSIS(JW),     JW=1,NWB)
  READ (CON,*,ERR=400)         (PSIDK(JW),    JW=1,NWB)
  READ (CON,*,ERR=400)         (PARTSI(JW),   JW=1,NWB)
     
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)   
    WRITE (WIN,*) '  CO2R'
  AID='CO2R'
  READ (CON,*,ERR=400)          (CO2R(JW),     JW=1,NWB)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)   
      WRITE (WIN,*) '  O2LIM'
  AID='O2LIM'

  READ (CON,*,ERR=400)           O2LIM

  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)   
  WRITE (WIN,*) '   sediment O2 demand first order model'

  READ (CON,*,ERR=400)     (SEDC(JW),   JW=1,NWB); SEDC=adjustr(SEDC)  
  READ (CON,*,ERR=400)     (PRNSC(JW),  JW=1,NWB); PRNSC=adjustr(PRNSC) 
  READ (CON,*,ERR=400)     (SEDCI(JW),   JW=1,NWB)  
  READ (CON,*,ERR=400)     (SEDDK(JW),   JW=1,NWB)     
  READ (CON,*,ERR=400)     (SEDS(JW),    JW=1,NWB)  
  READ (CON,*,ERR=400)     (FSOD(JW),    JW=1,NWB)  
  READ (CON,*,ERR=400)     (FSED(JW),    JW=1,NWB) 
  READ (CON,*,ERR=400)     (SEDBR(JW),   JW=1,NWB)  
  READ (CON,*,ERR=400)     (DYNSEDK(JW), JW=1,NWB); DYNSEDK=adjustr(DYNSEDK)
  READ (CON,*,ERR=400)     (SODT1(JW),   JW=1,NWB)
  READ (CON,*,ERR=400)     (SODT2(JW),   JW=1,NWB)
  READ (CON,*,ERR=400)     (SODK1(JW),   JW=1,NWB)
  READ (CON,*,ERR=400)     (SODK2(JW),   JW=1,NWB)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)   
    WRITE (WIN,*) '  sediment O2 demand zero order'

  READ (CON,*,ERR=400)  (SOD(I),  I=1,IMX)
  
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)    
    WRITE (WIN,*) '  reaeration'

  READ (CON,*,ERR=400)   (REAERC(JW), JW=1,NWB); REAERC=adjustr(REAERC)
  READ (CON,*,ERR=400)   (NEQN(JW),   JW=1,NWB)
  READ (CON,*,ERR=400)   (RCOEF1(JW), JW=1,NWB)
  READ (CON,*,ERR=400)   (RCOEF2(JW), JW=1,NWB)
  READ (CON,*,ERR=400)   (RCOEF3(JW), JW=1,NWB)
  READ (CON,*,ERR=400)   (RCOEF4(JW), JW=1,NWB)
  READ (CON,*,ERR=400)   (DGPO2(JW), JW=1,NWB)
  READ (CON,*,ERR=400)   (MINKL(JW), JW=1,NWB)
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)    
  
! Input filenames
  WRITE (WIN,*) '  file names'

  READ (CON,*,ERR=400)  QWDFN
  READ (CON,*,ERR=400)  QGTFN
  READ (CON,*,ERR=400)  WSCFN
  READ (CON,*,ERR=400)  SHDFN
  READ (CON,*,ERR=400)  VPLFN(1)
  VPLFN(2:NWB)=VPLFN(1)
  
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)    
 
  READ (CON,*,ERR=400) (BTHFN(JW), JW=1,NWB)
  READ (CON,*,ERR=400) (METFN(JW), JW=1,NWB)
  READ (CON,*,ERR=400) (EXTFN(JW), JW=1,NWB)
  READ (CON,*,ERR=400) (ATMDEPFN(JW), JW=1,NWB)
  READ (CON,*,ERR=400) (VPRFN(JW), JW=1,NWB)
  READ (CON,*,ERR=400) (LPRFN(JW), JW=1,NWB)

  ! Output filenames

  READ (CON,*,ERR=400) (SNPFN(JW), JW=1,NWB)
  READ (CON,*,ERR=400) (PRFFN(JW), JW=1,NWB)
  READ (CON,*,ERR=400) (CPLFN(JW), JW=1,NWB)
  READ (CON,*,ERR=400) (SPRFN(JW), JW=1,NWB)
  READ (CON,*,ERR=400) (FLXFN(JW), JW=1,NWB)
  
  READ (CON,*,ERR=400)
  READ (CON,*,ERR=400)    
 
  READ (CON,*,ERR=400) (QINFN(JB), JB=1,NBR)
  READ (CON,*,ERR=400) (TINFN(JB), JB=1,NBR)
  READ (CON,*,ERR=400) (CINFN(JB), JB=1,NBR)
  READ (CON,*,ERR=400) (QOTFN(JB), JB=1,NBR)
  READ (CON,*,ERR=400) (QDTFN(JB), JB=1,NBR)
  READ (CON,*,ERR=400) (TDTFN(JB), JB=1,NBR)
  READ (CON,*,ERR=400) (CDTFN(JB), JB=1,NBR)
  READ (CON,*,ERR=400) (PREFN(JB), JB=1,NBR)
  READ (CON,*,ERR=400) (TPRFN(JB), JB=1,NBR)
  READ (CON,*,ERR=400) (CPRFN(JB), JB=1,NBR)
  READ (CON,*,ERR=400) (EUHFN(JB), JB=1,NBR)
  READ (CON,*,ERR=400) (TUHFN(JB), JB=1,NBR)
  READ (CON,*,ERR=400) (CUHFN(JB), JB=1,NBR)
  READ (CON,*,ERR=400) (EDHFN(JB), JB=1,NBR)
  READ (CON,*,ERR=400) (TDHFN(JB), JB=1,NBR)
  READ (CON,*,ERR=400) (CDHFN(JB), JB=1,NBR)
 
  CLOSE (CON) 
ENDIF


  DO JW=1,NWB
    WRITE (WIN,*) 'Bathymetry file:', adjustl(trim(bthfn(jw)))
    UNIT = BTH(JW)
    OPEN (BTH(JW),FILE=BTHFN(JW),STATUS='OLD',IOSTAT=IERR)
    IF (IERR == 0) THEN
      READ  (BTH(JW),'(a1)')char1
      if(char1=='$')then
      WRITE( WIN,*) '      file in csv format'
      WRITE (WIN,*) '  segment lengths'
      READ  (BTH(JW),*,ERR=400)
      READ  (BTH(JW),*,ERR=400) AID,(DLX(I),  I=US(BS(JW))-1,DS(BE(JW))+1)
      WRITE (WIN,*) 'Finished reading:', AID
      WRITE (WIN,*) '  water surface elevations'
      READ  (BTH(JW),*,ERR=400) AID,(ELWS(I), I=US(BS(JW))-1,DS(BE(JW))+1)
      WRITE (WIN,*) 'Finished reading:', AID
      WRITE (WIN,*) '  segment orientation'
      READ  (BTH(JW),*,ERR=400) AID,(PHI0(I), I=US(BS(JW))-1,DS(BE(JW))+1)
      WRITE (WIN,*) 'Finished reading:', AID
      WRITE (WIN,*) '  segment bottom friction'
      READ  (BTH(JW),*,ERR=400) AID,(FRIC(I), I=US(BS(JW))-1,DS(BE(JW))+1)
      WRITE (WIN,*) 'Finished reading:', AID
      WRITE (WIN,*) '  layer thickness and segment widths'
      READ  (BTH(JW),*)
      DO K=1,KMX
      READ  (BTH(JW),*,ERR=400) H(K,JW),(B(K,I),I=US(BS(JW))-1,DS(BE(JW))+1)
      END DO
      else
      WRITE (WIN,*) '  segment lengths'
      READ  (BTH(JW),'(/A8/(10F8.0))',ERR=400) AID,(DLX(I),  I=US(BS(JW))-1,DS(BE(JW))+1)
      WRITE (WIN,*) '  water surface elevations'
      READ  (BTH(JW),'(/A8/(10F8.0))',ERR=400) AID,(ELWS(I), I=US(BS(JW))-1,DS(BE(JW))+1)
      WRITE (WIN,*) '  segment orientation'
      READ  (BTH(JW),'(/A8/(10F8.0))',ERR=400) AID,(PHI0(I), I=US(BS(JW))-1,DS(BE(JW))+1)
      WRITE (WIN,*) '  segment bottom friction'
      READ  (BTH(JW),'(/A8/(10F8.0))',ERR=400) AID,(FRIC(I), I=US(BS(JW))-1,DS(BE(JW))+1)
      WRITE (WIN,*) '  layer thickness'
      READ  (BTH(JW),'(/A8/(10F8.0))',ERR=400) AID,(H(K,JW), K=1,KMX)
      WRITE (WIN,*) '  segment widths'
      DO I=US(BS(JW))-1,DS(BE(JW))+1
        READ (BTH(JW),'(/A8/(10F8.0))',ERR=400) AID,(B(K,I),  K=1,KMX)
      END DO
      endif
      CLOSE (BTH(JW))
    ELSE
      CALL ERRORS
      WRITE (ERR,FMTA) 'Could not open bathymetry file '//BTHFN(JW)
      WRITE (WIN,FMTA) 'Could not open bathymetry file '//BTHFN(JW)
      STOP
    END IF
  END DO

! Names/units text files
IF(CONFN=='w2_con.npt')THEN
  WRITE (WIN,*) 'Graph file: graph.npt'
  UNIT=NPT
  OPEN (NPT,FILE='graph.npt',STATUS='OLD',IOSTAT=IERR)
  IF (IERR == 0) THEN                                            ! SW 1/16/04 replaced section
      AID='graph.npt'
      ICST=2   ! FLAG FOR ERROR OF CONSTITUENT ORDER
  READ  (NPT,'(///(A43,1X,A9,3F8.0,A8))',ERR=400) (HNAME(J),  FMTH(J),  HMULT(J),  HYMIN(J), HYMAX(J), HPLTC(J), J=1,NHY)
  READ  (NPT,'(// (A43,1X,A9,3F8.0,A8))',ERR=400) (CNAME(J),  FMTC(J),  CMULT(J),  CMIN(J),  CMAX(J),  CPLTC(J), J=1,NCT)
  READ  (NPT,'(// (A43,1X,A9,3F8.0,A8))',END=400,ERR=400) (CDNAME(J), FMTCD(J), CDMULT(J), CDMIN(J), CDMAX(J), CDPLTC(J),J=1,NDC)
      ICST=0
  !CLOSE (NPT) 
  !DO JC=1,NCT
  !  L1         = SCAN (CNAME(JC),',')+2
  !  L2         = SCAN (CNAME(JC)(L1:43),'  ')+L1
  !  CUNIT(JC)  = CNAME(JC)(L1:L2)
  !  CNAME1(JC) = CNAME(JC)(1:L1-3)
  !  CUNIT1(JC) = CUNIT(JC)(1:1)
  !  CUNIT2(JC) = CUNIT(JC)
  !  IF (CUNIT(JC)(1:2) == 'mg') THEN
  !    CUNIT1(JC) = 'g'
  !    CUNIT2(JC) = 'g/m^3'
  !  END IF
  !  IF (CUNIT(JC)(1:2) /= 'g/' .AND. CUNIT(JC)(1:2) /= 'mg') CUNIT1(JC) = '  '
  !END DO
  !FMTH(1:NHY) = ADJUSTL (FMTH(1:NHY))
  ELSE
    CALL ERRORS
    WRITE (ERR,FMTA) 'Could not open graph.npt'
    WRITE (WIN,FMTA) 'Could not open graph.npt'
  END IF
  CLOSE(NPT)
ENDIF
  
  DO JC=1,NCT
    L1         = SCAN (CNAME(JC),',')+2
    L2         = SCAN (CNAME(JC)(L1:43),'  ')+L1
    CUNIT(JC)  = CNAME(JC)(L1:L2)
    CNAME1(JC) = CNAME(JC)(1:L1-3)
    CUNIT1(JC) = CUNIT(JC)(1:1)
    CUNIT2(JC) = CUNIT(JC)
    IF (CUNIT(JC)(1:2) == 'mg') THEN
      CUNIT1(JC) = 'g'
      CUNIT2(JC) = 'g/m^3'
    END IF
    IF (CUNIT(JC)(1:2) /= 'g/' .AND. CUNIT(JC)(1:2) /= 'mg') CUNIT1(JC) = '  '
  END DO
  FMTH(1:NHY) = ADJUSTL (FMTH(1:NHY))
  
  ! WSC FILE

  WRITE (WIN,*) 'Wind sheltering file:', adjustl(trim(wscfn))
  OPEN (NPT,FILE=WSCFN,STATUS='OLD',IOSTAT=IERR)
  IF (IERR == 0) THEN
      UNIT=NPT
      AID=WSCFN
  READ(NPT,'(A1)')ICHAR1
  IF(ICHAR1=='$')THEN
  WRITE(WIN,*) '  wsc file in csv format'
  READ (NPT,'(/)')
  READ (NPT,*,ERR=400,END=29000) SDAY,(WSC(I),I=1,IMX)
  ELSE
  READ (NPT,'(//10F8.0:/(8X,9F8.0))',ERR=400,END=29000) SDAY,(WSC(I),I=1,IMX)
  ENDIF
        IF (SDAY > TMSTRT) THEN
          CALL ERRORS
          WRITE (ERR,'(3(A,F0.3))') 'Starting date [',SDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//WSCFN
        ENDIF
        DO I=1,IMX
        IF(WSC(I) <= 0.0)THEN
        CALL ERRORS
        WRITE (ERR,'(A,F0.3,A,I4,A)') 'Julian date ',SDAY,': WSC AT SEG(I)=',I,' <= 0.0 in '//WSCFN
        ENDIF
        IF (WSC(I) > 2.0) THEN
          CALL WARNINGS
          WRITE (WRN,'(A,F0.3,A,I4,A)') 'Julian day ',SDAY,': WSC(I) AT SEG(I)=',I,' > 2.0 in '//WSCFN
        END IF
        IF (WSC(I) > 0.0 .and. wsc(i) < 0.5) THEN
          CALL WARNINGS
          WRITE (WRN,'(A,F0.3,A,I4,A)') 'Julian day ',SDAY,': WSC(I) AT SEG(I)=',I,' < 0.5 in '//WSCFN
        END IF
        ENDDO

  SDAYO=SDAY

!  DO J=1,100
28995 continue  ! cb 3/26/10
    IF(ICHAR1=='$')THEN
    READ (NPT,*,END=29000,ERR=400) SDAY,(WSC(I),I=1,IMX)
    ELSE
    READ (NPT,'(10F8.0:/(8X,9F8.0))',END=29000,ERR=400) SDAY,(WSC(I),I=1,IMX)
    ENDIF
          IF (SDAY <= SDAYO) THEN
          CALL ERRORS
          WRITE (ERR,'(3(A,F0.3))') 'Julian date ',SDAY,' <= previous date of ',SDAYO,' in '//WSCFN
          END IF
        DO I=1,IMX
        IF(WSC(I) <= 0.0)THEN
        CALL ERRORS
        WRITE (ERR,'(A,F0.3,A,I4,A)') 'Julian date ',SDAY,': WSC AT SEG(I)=',I,' <= 0.0 in '//WSCFN
        ENDIF
        IF (WSC(I) > 2.0) THEN
          CALL WARNINGS
          WRITE (WRN,'(A,F0.3,A,I4,A)') 'Julian day ',SDAY,': WSC(I) AT SEG(I)=',I,' > 2.0 in '//WSCFN
        END IF
        IF (WSC(I) > 0.0 .and. wsc(i) < 0.5) THEN
          CALL WARNINGS
          WRITE (WRN,'(A,F0.3,A,I4,A)') 'Julian day ',SDAY,': WSC(I) AT SEG(I)=',I,' < 0.5 in '//WSCFN
        END IF
        ENDDO
        SDAYO=SDAY
  ! ENDDO

  go to 28995   ! cb 3/26/10
  29000 IF (SDAY < TMEND) THEN
        CALL ERRORS
        WRITE (ERR,'(3(A,F0.3))') 'Ending time ',SDAY,' < ending simulation time [TMEND=',TMEND,'] in '//WSCFN
        WRITE(ERR,*) '[The above error might be caused by a blank line at the end of the file.]'
        END IF
        IF(I < IMX) THEN
        CALL ERRORS
        WRITE (ERR,'(3(A,F0.3))') 'Not enough records (IMX) for last day ',SDAY,' in '//WSCFN
        END IF
  CLOSE (NPT)

  ELSE
    CALL ERRORS
    WRITE (ERR,'(A)') 'Could not open wsc file:'//adjustl(trim(wscfn))
    WRITE (WIN,'(A)') 'Could not open wsc file:'//adjustl(trim(wscfn))
  END IF
  CLOSE (NPT)


!***********************************************************************************************************************************
!*                                                 Task 2:  Variable Initialization                                               **
!***********************************************************************************************************************************

  NWBR   = 0;   NTBR    = 0;   KB      = 0;   KBMAX   = 0;   JBU     = 0;   JBD     = 0;   NASS    = 0;   NAAL    = 0; KEPR    = 0
  NCCBR  = 0;   NCCGR   = 0;   CN      = 0;   INCN    = 0;   incdt   = 0;   TRCN    = 0;   NAC     = 0;   NADC    = 0;   NACIN   = 0; NACTR   = 0
  NACDT  = 0;   NACPR   = 0;   NQINMX  = 0;   NQTRMX  = 0;   NQDTMX  = 0;   NQPRMX  = 0;   NQWDMX  = 0;   NCCGR   = 0
  XI     = 0.0; XBR     = 0.0; XGR     = 0.0; T2      = 0.0; C2      = 0.0; CD      = 0.0; EL      = 0.0; VOLB    = 0.0
  VOLG   = 0.0; CVBR    = 0.0; CVGR    = 0.0; SABR    = 0.0; SAGR    = 0.0; QINMX   = 0.0; QDTMX   = 0.0; QTRMX   = 0.0
  QWDMX  = 0.0; QWDMXB  = 0.0; QPRMX   = 0.0; QSTRMX  = 0.0; QINS    = 0.0; QTRS    = 0.0; QDTS    = 0.0; QWDS    = 0.0
  QPRS   = 0.0; QSTRS   = 0.0; QINAVB  = 0.0; QINAV   = 0.0; QDTAV   = 0.0; QPRAV   = 0.0; QTRAVB  = 0.0; QOTAVB  = 0.0
  QSTRAV = 0.0; QWDAVB  = 0.0; QDNAVB  = 0.0; QTRAVB  = 0.0; QINAVW  = 0.0; QINAVB  = 0.0; QOTAVW  = 0.0; QINMXB  = 0.0
  QOTMXB = 0.0; QINMXW  = 0.0; QOTMXW  = 0.0; QTRMXB  = 0.0; CINMAX  = 0.0; CINAVG  = 0.0; SAGR    = 0.0; CVGR    = 0.0; cin=0.0;cin2=0.0
  HGR    = 0.0; BGR     = 0.0; CVGR    = 0.0; HBR     = 0.0; cinload=0.0; ctrload=0.0;cinpload=0.0;cinnload=0.0;ctrpload=0.0;ctrnload=0.0;cdtload=0.0
  cdtavg=0.0;cdtmin=100.0;cdtmax=0.0
  MAC = 0.0; EPD = 0.0
  CONV   = BLANK
  FMTA   = '(A)'
  FMTFI  = '(3(A,F0.3,A,I0))'
  FMTIF  = '(3(A,I0,A,F0.3))'
  FMTF   = '(3(A,F0.3))'
  FMTI   = '(4(A,I0))'
  FMT2FI = '(2(A,F0.3),A,I0)'
  FMTF2I = '(A,F0.2,2(A,I0))'
  FMT2IF = '(3(2(A,I0),A,F0.3))'
  TRTMAX=0.0
  TRTMIN=100.
  TINMAX=0.0
  TINMIN=100.
  NACATD=0;ATMDCN =0

! Initialize logical control variables

  SEDIMENT_CALC    = .FALSE.; DISSOLVED_SOLIDS = .FALSE.; GENERIC_CONST    = .FALSE.; SUSPENDED_SOLIDS = .FALSE.
  PHOSPHORUS       = .FALSE.; AMMONIUM         = .FALSE.; NITRATE          = .FALSE.; IRON             = .FALSE.
  LABILE_DOM       = .FALSE.; REFRACTORY_DOM   = .FALSE.; LABILE_POM       = .FALSE.; REFRACTORY_POM   = .FALSE.
  CBO_DEMAND       = .FALSE.; ALGAE            = .FALSE.; DISSOLVED_OXYGEN = .FALSE.; TOT_IC           = .FALSE.
  OPEN_VPR         = .FALSE.; OPEN_LPR         = .FALSE.; OPEN_LPRC         = .FALSE.
  WITHDRAWALS      =  NWD   >   0;   TRIBUTARIES =  NTR   >   0;   PIPES     =  NPI   >   0
  GATES            =  NGT   >   0;   SPILLWAYS   =  NSP   >   0;   PUMPS     =  NPU   >   0
  ISO_TEMP         =  T2I   >=  0.0; VERT_TEMP   =  T2I   == -1.0; LONG_TEMP =  T2I   == -2.0
  LONG_SED         =  SEDCI == -2.0
  OPEN_VPR         =  VERT_TEMP
  OPEN_LPR         =  LONG_TEMP
  PRECIPITATION    =  PRC  == '      ON'
  CONSTITUENTS     =  CCC  == '      ON'
  ICE_CALC         =  ICEC == '      ON'
  ICE_CALC         =  ICEC == '    ONWB'
  RESTART_IN       =  RSIC == '      ON'
  DO JW=1,NWB
    IF (LONG_SED(JW)) OPEN_LPRC(JW) = .TRUE.
    ISO_SEDIMENT(JW)     = SEDCI(JW)   >=  0   .AND. SEDC(JW)   == '      ON'
    VERT_SEDIMENT(JW)    = SEDCI(JW)   == -1.0 .AND. SEDC(JW)   == '      ON'
    LONG_SEDIMENT(JW)    = SEDCI(JW)   <  -1.0 .AND. SEDC(JW)   == '      ON'
    ISO_EPIPHYTON(JW,:)  = EPICI(JW,:) >=  0   .AND. EPIC(JW,:) == '      ON'
    VERT_EPIPHYTON(JW,:) = EPICI(JW,:) == -1.0 .AND. EPIC(JW,:) == '      ON'
    LONG_EPIPHYTON(JW,:) = EPICI(JW,:) <  -1.0 .AND. EPIC(JW,:) == '      ON'
          IF (any(VERT_EPIPHYTON(JW,:))) OPEN_VPR(JW) = .TRUE.
          IF (any(LONG_EPIPHYTON(JW,:))) OPEN_LPRC(JW) = .TRUE.
    ISO_MACRO(JW,:)  = MACWBCI(JW,:) >=  0   .AND. MACWBC(JW,:) == '      ON'
    VERT_MACRO(JW,:) = MACWBCI(JW,:) == -1.0 .AND. MACWBC(JW,:) == '      ON'
    LONG_MACRO(JW,:) = MACWBCI(JW,:) <  -1.0 .AND. MACWBC(JW,:) == '      ON'
          IF (any(VERT_MACRO(JW,:))) OPEN_VPR(JW) = .TRUE.
          IF (any(LONG_MACRO(JW,:))) OPEN_LPRC(JW) = .TRUE.
    DO JC=1,NCT
      ISO_CONC(JC,JW)  = C2IWB(JC,JW) >=  0.0
      VERT_CONC(JC,JW) = C2IWB(JC,JW) == -1.0  .AND.  CAC(JC) == '      ON'
      LONG_CONC(JC,JW) = C2IWB(JC,JW) <  -1.0  .AND.  CAC(JC) == '      ON'
      IF (VERT_CONC(JC,JW)) OPEN_VPR(JW) = .TRUE.
      IF (LONG_CONC(JC,JW)) OPEN_LPRC(JW) = .TRUE.
    END DO
  END DO
  DO JB=1,NBR
    DAM_FLOW(JB)    = UHS(JB) <  -1; UH_EXTERNAL(JB) = UHS(JB) == -1; DH_EXTERNAL(JB) = DHS(JB) == -1
    UH_INTERNAL(JB) = UHS(JB) >   0; DH_INTERNAL(JB) = DHS(JB) >   0; UQ_EXTERNAL(JB) = UHS(JB) ==  0
    DQ_EXTERNAL(JB) = DHS(JB) ==  0; DQ_INTERNAL(JB) = DQB(JB) >   0
    UQ_INTERNAL(JB) = UQB(JB) >   0 .AND. .NOT. DAM_FLOW(JB)
    DIST_TRIB(JB)   = DTRC(JB) == '      ON'
    UHS(JB)         = ABS(UHS(JB))
    DHS(JB)         = ABS(DHS(JB))
    INTERP_INFLOW      = QINIC  == '      ON'   ! cb 9/7/2010
  END DO

! Inflow/outflow branch locations

  IF (WITHDRAWALS) THEN
    JBWD = 0
    NWBR = 0
    DO JB=1,NBR
      DO JW=1,NWD
        IF (IWD(JW) >= US(JB) .AND. IWD(JW) <= DS(JB)) THEN
          JBWD(JW) = JB
          NWBR(JB) = NWBR(JB)+1
        END IF
      END DO
    END DO
  END IF
  IF (TRIBUTARIES) THEN
    JBTR = 0
    NTBR = 0
    DO JB=1,NBR
      DO JT=1,NTR
        IF (ITR(JT) >= US(JB) .AND. ITR(JT) <= DS(JB)) THEN
          JBTR(JT) = JB
          NTBR(JB) = NTBR(JB)+1
        END IF
      END DO
    END DO
  END IF
  IF (PIPES) THEN
    JBUPI  = 0
    NUPIBR = 0
    JBDPI  = 0
    NDPIBR = 0
    DO JB=1,NBR
      DO JP=1,NPI
        IF (IUPI(JP) >= US(JB) .AND. IUPI(JP) <= DS(JB)) THEN
          JBUPI(JP)  = JB
          NUPIBR(JB) = NUPIBR(JB)+1
        END IF
        IF (IDPI(JP) >= US(JB) .AND. IDPI(JP) <= DS(JB)) THEN
          JBDPI(JP)  = JB
          NDPIBR(JB) = NDPIBR(JB)+1
        END IF
      END DO
    END DO
  END IF
  IF (GATES) THEN
    JBUGT  = 0
    NUGTBR = 0
    JBDGT  = 0
    NDGTBR = 0
    DO JB=1,NBR
      DO JG=1,NGT
        IF (IUGT(JG) >= US(JB) .AND. IUGT(JG) <= DS(JB)) THEN
          JBUGT(JG)  = JB
          NUGTBR(JB) = NUGTBR(JB)+1
        END IF
        IF (IDGT(JG) >= US(JB) .AND. IDGT(JG) <= DS(JB)) THEN
          JBDGT(JG)  = JB
          NDGTBR(JB) = NDGTBR(JB)+1
        END IF
      END DO
    END DO
  END IF

! Active constituents

  IF (CAC(NTDS)  == '      ON') DISSOLVED_SOLIDS = .TRUE.
  IF (CAC(NPO4)  == '      ON') PHOSPHORUS       = .TRUE.
  IF (CAC(NNH4)  == '      ON') AMMONIUM         = .TRUE.
  IF (CAC(NNO3)  == '      ON') NITRATE          = .TRUE.
  IF (CAC(NLDOM) == '      ON') LABILE_DOM       = .TRUE.
  IF (CAC(NRDOM) == '      ON') REFRACTORY_DOM   = .TRUE.
  IF (CAC(NLPOM) == '      ON') LABILE_POM       = .TRUE.
  IF (CAC(NRPOM) == '      ON') REFRACTORY_POM   = .TRUE.
  IF (CAC(NDO)   == '      ON') DISSOLVED_OXYGEN = .TRUE.
  IF (CAC(NTIC)  == '      ON') TOT_IC           = .TRUE.
  DO JG=NGCS,NGCE
    IF (CAC(JG)  == '      ON') THEN
      GENERIC_CONST = .TRUE.
      EXIT
    END IF
  END DO
  DO JS=NSSS,NSSE
    IF (CAC(JS) == '      ON') THEN
      SUSPENDED_SOLIDS = .TRUE.
      EXIT
    END IF
  END DO
  DO JB=NBODS,NBODE
    IF (CAC(JB) == '      ON') THEN
      CBO_DEMAND = .TRUE.
      EXIT
    END IF
  END DO
  DO JA=NAS,NAE
    IF (CAC(JA) == '      ON') THEN
      ALGAE = .TRUE.
      EXIT
    END IF
  END DO
  DO JW=1,NWB
    IF (SEDC(JW) == '      ON') SEDIMENT_CALC(JW) = .TRUE.
  END DO
  ALGAE_MIGRATION=.FALSE.
  IF (CONSTITUENTS) THEN
    DO JC=1,NCT
      IF (CAC(JC) == '      ON') THEN
        NAC     = NAC+1
        CN(NAC) = JC
      END IF
      DO JB=1,NBR
        IF (CINBRC(JC,JB) == '      ON') THEN                                                                          !SW 01/07/01
          NACIN(JB)          = NACIN(JB)+1                                                                             !SW 01/07/01
          INCN(NACIN(JB),JB) = JC                                                                                      !SW 01/07/01
        END IF
        IF (CDTBRC(JC,JB) == '      ON')then
            NACDT(JB) = NACDT(JB)+1
            incdt(nacdt(jb),jb)=jc
        ENDIF
      END DO
      DO JT=1,NTR                                                                                                      !SW 01/07/01
        IF (CINTRC(JC,JT) == '      ON') THEN                                                                          !SW 01/07/01
          NACTR(JT)          = NACTR(JT)+1                                                                             !SW 01/07/01
          TRCN(NACTR(JT),JT) = JC                                                                                      !SW 01/07/01
        END IF
      END DO                                                                                                           !SW 01/07/01
    END DO
    DO JW=1,NWB
      DO JC=1,NCT
        IF (CPRBRC(JC,JW) == '      ON') NACPR(JW) = NACPR(JW)+1
      END DO
      DO JD=1,NDC
        IF (CDWBC(JD,JW) == '      ON') THEN
          NADC(JW)         = NADC(JW)+1
          CDN(NADC(JW),JW) = JD                                                                                        !TC 09/23/02
        END IF
      END DO
      IF(ATM_DEPOSITIONC(JW)=='      ON')THEN
        DO JC=1,NCT
        IF(C_ATM_DEPOSITION(JC,JW)=='      ON')THEN  
        NACATD(JW)          = NACATD(JW)+1
        ATMDCN(NACATD(JW),JW) = JC         
        ENDIF
        ENDDO
        ENDIF   
    END DO
    DO JS=1,NSS
      IF (CAC(NSSS) == '      ON') THEN
        NASS       = NASS+1
        SSAC(NASS) = JS
      END IF
    END DO
    DO JA=1,NAL
      IF (CAC(NAS) == '      ON') THEN
        NAAL       = NAAL+1
        ALAC(NAAL) = JA
      END IF
    END DO
    ALGAE_TOXINS = .FALSE.
    DO J=NATS,NATE
        IF(CAC(J) == '      ON')THEN
            ALGAE_TOXINS = .TRUE.
            EXIT
        ENDIF
    ENDDO

    DO JA=1,NAL
        IF(AVERTM(JA)=='      ON')THEN
            ALGAE_MIGRATION=.TRUE.
            EXIT
        ENDIF
    ENDDO
    
    
  END IF

  ALLOCATE (CUH(KMX,NAC),CDH(KMX,NAC),CUHMIN(NAC,NBR),CUHMAX(NAC,NBR),CUHAVG(NAC,NBR),CDHMIN(NAC,NBR),CDHMAX(NAC,NBR),CDHAVG(NAC,NBR))


! Grid linkage

  DO JB=1,NBR

  ! Check if UHS is negative

  IF(DAM_FLOW(JB))THEN
    ! Find branch that this is connected to
    DO JJB=1,NBR
        IF (UHS(JB) == DS(JJB))THEN
          IF(NSTR(JJB) == 0)THEN
           CALL ERRORS
           WRITE (ERR,FMTI) 'Error grid linkage:Since UHS(JB) is < 0, # STR must be > 0 for connected branch(JJB). NSTR(JJB)=0. JB:',JB,' JJB:',JJB
          ENDIF
        EXIT
        ENDIF
        IF(JJB == NBR)THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Error grid linkage:Since UHS is < 0, this MUST be the DS of another branch. JB:',JB
        ENDIF
    END DO

  END IF

    JCHECK = 0
    IF (UHS(JB) /= 0 .AND. .NOT. UH_EXTERNAL(JB) .AND. .NOT. DH_EXTERNAL(JB)) THEN                                     !TC 08/31/01
      DO JJB=1,NBR
        IF (ABS(UHS(JB)) >= US(JJB) .AND. ABS(UHS(JB)) <= DS(JJB)) EXIT
      END DO
      if(jjb.gt.nbr)then                                                                                               ! SW 6/1/07
            CALL ERRORS
            WRITE (ERR,FMTI) 'Error in grid linkage:check[US],[DS],[UHS],[DHS] and/or[JBDN]for waterbody ',JW,' JB=',JB
            WRITE (WIN,FMTI) 'Error in grid linkage:check[US],[DS],[UHS],[DHS] and/or[JBDN]for waterbody ',JW,' JB=',JB
            STOP
      endif
      IF (ABS(UHS(JB)) == DS(JJB)) THEN
        IF (DHS(JJB) == US(JB)) JCHECK = 1
        IF (UHS(JB) < 0) THEN
          IF (JCHECK == 1) THEN
            CALL ERRORS
            WRITE (ERR,FMTI) 'Error in grid linkage:check[US],[DS],[UHS],[DHS] and/or[JBDN]for waterbody ',JW,' JB=',JB
            WRITE (WIN,FMTI) 'Error in grid linkage:check[US],[DS],[UHS],[DHS] and/or[JBDN]for waterbody ',JW,' JB=',JB
            STOP
          END IF
        END IF
      END IF
    END IF
  END DO
  ZERO_SLOPE = .TRUE.
  DO JW=1,NWB
    DO JB=BS(JW),BE(JW)                                                                                                !SW 06/12/01
      IF (SLOPE(JB) /= 0.0) ZERO_SLOPE(JW) = .FALSE.                                                                   !SW 10/17/02
    END DO                                                                                                             !SW 06/12/01
  END DO                                                                                                               !SW 06/12/01

!***********************************************************************************************************************************
!*                                                Task 3: Input/Output File Checks                                                **
!***********************************************************************************************************************************

  ! DHS and UHS checks
  
  DO JB=1,NBR
      IF(DHS(JB) > 0 .and. .NOT.DH_EXTERNAL(JB))THEN
          DO JJB=1,NBR
              IF(JJB == JB)CYCLE
              IF(DHS(JB) >= US(JJB) .AND. DHS(JB) <= DS(JJB))THEN
                  GO TO 2345
              ENDIF
          ENDDO
      CALL ERRORS
      WRITE (ERR,*) 'DHS(JB) incorrectly specified for branch:', jb,' DHS(JB)=',DHS(JB)
      EXIT
2345  CONTINUE
      ENDIF
ENDDO
    DO JB=1,NBR
      IF(UHS(JB) >0 .AND. .NOT.UH_EXTERNAL(JB))THEN
          DO JJB=1,NBR
              IF(JJB == JB)CYCLE
              IF(UHS(JB) >= US(JJB) .AND. UHS(JB) <= DS(JJB))THEN
                  GO TO 2346
              ENDIF
          ENDDO
      CALL ERRORS
      WRITE (ERR,*) 'UHS(JB) incorrectly specified for branch:', jb,' UHS(JB)=',UHS(JB)
      EXIT
2346  CONTINUE
    ENDIF
ENDDO
  
  
  
! Layer elevations

  DO JW=1,NWB                                                                                                          !SW 10/17/02
    IF (ZERO_SLOPE(JW)) THEN
      DO I=US(BS(JW))-1,DS(BE(JW))+1
        EL(KMX,I) = ELBOT(JW)
        DO K=KMX-1,1,-1
          EL(K,I) = EL(K+1,I)+H(K,JW)
        END DO
      END DO
    ELSE
      NPOINT = 0                                                                                                       !TC 02/11/01
      NUP              = 0
      NNBR             = 1
      NCBP             = 0
      NINTERNAL        = 0
      JB               = JBDN(JW)
      NPOINT(JB)       = 1
      EL(KMX,DS(JB)+1) = ELBOT(JW)
      DO WHILE (NNBR <= (BE(JW)-BS(JW)+1))
        NCBP = NCBP+1
        IF (NCBP > NBR) THEN
          CALL ERRORS
          WRITE (ERR,FMTI) 'Model linkage not logical for waterbody ',JW
          EXIT
        END IF
        IF (NINTERNAL == 0) THEN
          IF (NUP == 0) THEN
            DO I=DS(JB),US(JB)-1,-1
              IF (I /= DS(JB)) THEN
                EL(KMX,I) = EL(KMX,I+1)+SLOPE(JB)*(DLX(I)+DLX(I+1))*0.5
              ELSE
                EL(KMX,I) = EL(KMX,I+1)                                                                                !bug?
              END IF
              DO K=KMX-1,1,-1
                EL(K,I) = EL(K+1,I)+H(K,JW)*COS(ATAN(SLOPE(JB)))
              END DO
            END DO
          ELSE
            DO I=US(JB),DS(JB)+1
              IF (I /= US(JB)) THEN
                EL(KMX,I) = EL(KMX,I-1)-SLOPE(JB)*(DLX(I)+DLX(I+1))*0.5
              ELSE
                EL(KMX,I) = EL(KMX,I-1)
              END IF
              DO K=KMX-1,1,-1
                EL(K,I) = EL(K+1,I)+H(K,JW)*COS(ATAN(SLOPE(JB)))
              END DO
            END DO
            NUP = 0
          END IF
          DO K=KMX,1,-1
            IF (.NOT.DAM_FLOW(JB)) THEN
              IF (UHS(JB) /= 0) THEN
                EL(K,US(JB)-1) = EL(K,US(JB))+SLOPE(JB)*DLX(US(JB))
              ELSE
                EL(K,US(JB)-1) = EL(K,US(JB))
              END IF
            ELSE
              EL(K,US(JB)-1) = EL(K,US(JB))
            END IF
            IF (DHS(JB) /= 0) THEN
              EL(K,DS(JB)+1) = EL(K,DS(JB))-SLOPE(JB)*DLX(DS(JB))
            ELSE
              EL(K,DS(JB)+1) = EL(K,DS(JB))
            END IF
          END DO
        ELSE
          DO K=KMX-1,1,-1
            EL(K,UHS(JJB)) = EL(K+1,UHS(JJB))+H(K,JW)*COS(ATAN(SLOPE(JB)))
          END DO
          DO I = UHS(JJB)+1, DS(JB)
            EL(KMX,I) = EL(KMX,I-1)-SLOPE(JB)*(DLX(I)+DLX(I-1))*0.5
            DO K=KMX-1,1,-1
              EL(K,I) = EL(K+1,I)+H(K,JW)*COS(ATAN(SLOPE(JB)))
            END DO
          END DO
          DO I=UHS(JJB)-1,US(JB),-1
            EL(KMX,I) = EL(KMX,I+1)+SLOPE(JB)*(DLX(I)+DLX(I+1))*0.5
            DO K=KMX-1,1,-1
              EL(K,I) = EL(K+1,I)+H(K,JW)*COS(ATAN(SLOPE(JB)))
            END DO
          END DO
          NINTERNAL = 0
        END IF
        IF (NNBR == (BE(JW)-BS(JW)+1)) EXIT
        NNBROLD = NNBR
        DO JB=BS(JW),BE(JW)
          IF (NPOINT(JB) /= 1) THEN
            DO JJB = BS(JW), BE(JW)
              IF (DHS(JB) >= US(JJB) .AND. DHS(JB) <= DS(JJB) .AND. NPOINT(JJB) == 1) THEN
                NPOINT(JB)       = 1
                NNBR             = NNBR+1
                EL(KMX,DS(JB)+1) = EL(KMX,DHS(JB))+SLOPE(JB)*(DLX(DS(JB))+DLX(DHS(JB)))*0.5
                EXIT
              END IF
              IF (UHS(JJB) == DS(JB) .AND. NPOINT(JJB) == 1) THEN
                NPOINT(JB)       = 1
                NNBR             = NNBR+1
                EL(KMX,DS(JB)+1) = EL(KMX,US(JJB))+SLOPE(JB)*(DLX(US(JJB))+DLX(DS(JB)))*0.5
                EXIT
              END IF
              IF (UHS(JJB) <= DS(JB) .AND. UHS(JJB) >= US(JB) .AND. NPOINT(JJB) == 1) THEN
                NNBR             = NNBR+1
                NINTERNAL        = 1
                NPOINT(JB)       = 1
                EL(KMX,UHS(JJB)) = EL(KMX,US(JJB))+SLOPE(JJB)*DLX(US(JJB))*0.5
                EXIT
              END IF
              IF (UHS(JB) <= DS(JJB) .AND. UHS(JB) >= US(JJB) .AND. NPOINT(JJB) == 1) THEN
                NNBR             = NNBR+1
                NPOINT(JB)       = 1
                EL(KMX,US(JB)-1) = EL(KMX,UHS(JB))-SLOPE(JB)*DLX(US(JB))*0.5
                IF (SLOPE(JB) /= 0.0) THEN
                  CALL WARNINGS
                  WRITE (WRN,FMTI) 'Branch ',JB,' has non-zero slope and is only'
                  WRITE (WRN,FMTI) 'connected through a UHS BC. This may not be'
                  WRITE (WRN,FMTI) 'a problem:if it is,change branch to a waterbody'                                   !TC 02/09/01
                END IF
                NUP = 1
                EXIT
              END IF
            END DO
            IF (NPOINT(JB) == 1) EXIT
          END IF
        END DO
        IF (NCBP > NBR .OR. NNBROLD == NNBR) THEN
          CALL ERRORS
          WRITE (ERR,FMTI) 'Error in grid linkage: check[US],[DS],[UHS],[DHS],[JBDN];each branch must be connected to another branch for waterbody ',JW
        END IF
      END DO
    END IF
  END DO

  XF=0.0
  DO JW=1,NWB
  DO JB=BS(JW),BE(JW)
      IF(JB==1)THEN
          XF(2)=DLX(2)*0.5
      ELSE
          XF(US(JB))=XF(DS(JB-1))+(DLX(DS(JB-1))+DLX(US(JB)))*0.5
      ENDIF
      DO I=US(JB)+1,DS(JB)
        XF(I)   = XF(I-1)+(DLX(I)+DLX(I-1))*0.5
      END DO
  ENDDO
  ENDDO
   
! CORRECTED SLOPE CHECK

DO JB=1,NBR
    IF(SLOPE(JB)<SLOPEC(JB))THEN   ! CORRECTED SLOPE CAN NOT BE GREATER THAN ACTUAL SLOPE
    CALL ERRORS
    WRITE(ERR,'(A,I3)')'SLOPE(JB) CANNOT BE LESS THAN SLOPEC(JB) FOR JB=',JB
    ENDIF
 ENDDO


! NPROC and closec

IF(NPROC < 0)THEN
 CALL ERRORS
 WRITE(ERR,*)'ERROR: NPROC < 0, MUST BE > 0'
ELSEIF(NPROC==0)THEN
 CALL ERRORS
 WRITE(ERR,*)'NPROC = 0, MUST BE > 0 [W2 CODE WILL SET NPROC=1 WHEN NPROC=0]'
ELSEIF(NPROC > 4)THEN
 CALL WARNINGS
 WRITE(WRN,*)'NPROC > 4. MAKE SURE THIS DOES NOT LEAD TO LONG RUN TIMES. TRY 2.'
ENDIF
IF(CLOSEC /= '     OFF' .AND. CLOSEC /= '      ON')THEN
 CALL ERRORS
 WRITE(ERR,*)'CLOSEC MUST BE "OFF" OR "ON". CLOSEC=',CLOSEC
ENDIF

! SELECTC CONTROL

IF(SELECTC /= '     OFF' .AND. SELECTC /= '      ON' .AND. SELECTC /= '    USGS')THEN
 CALL ERRORS
 WRITE(ERR,*)'SELECTC must be either "OFF" or "ON" or "USGS". SELECTC=',SELECTC
ENDIF

IF(SELECTC == '      ON' .or. SELECTC == '    USGS')THEN      ! CHECK FILE SELECTC
unit=1010
    OPEN (1010,FILE='w2_selective.npt',STATUS='UNKNOWN',IOSTAT=IERR)
  IF (IERR /= 0) THEN
    WRITE (WIN,*) 'Could not open w2_selective.npt input file'
    CALL ERRORS
    WRITE(ERR,*)'Could not open w2_selective.npt input file'
  ELSE
    ! check file
     AID='W2_SELECTIVE.NPT'
      IF(SELECTC=='      ON')THEN
          WRITE(WIN,*)'Reading w2_selective.npt'
      ELSE
          WRITE(WIN,*)'Reading w2_selective.npt - USGS format'
      ENDIF

      read(1010,'(a)',err=400)CHAR1
      read(1010,*,err=400)
      read(1010,*,err=400)
      IF(CHAR1=='$')GO TO 7770
      ! check for commas another sign that it is comma delimited
      read(1010,'(A)',ERR=400)HEADER(1)
     DO J=1,30
         IF(HEADER(1)(J:J)==',')THEN
             CHAR1='$'
             EXIT
         ENDIF
     ENDDO
     BACKSPACE(1010)

7770  CONTINUE
      if(char1=='$')then
      read(1010,*,err=400)aid1,tfrqtmp
      else
      read(1010,'(a8,f8.0)',err=400)tfrqtmp
      endif
      do j=1,2
      read(1010,*,err=400)
      end do
      if(char1=='$')then
      read(1010,*,err=400)aid1,aid1,numtempc,tcdfreq
      AID1=ADJUSTR(AID1)
      else
      read(1010,'(8x,a8,i8,f8.0)',err=400)aid1,numtempc,tcdfreq
      endif

      allocate(tcjb(numtempc),tcjs(numtempc),tcntr(numtempc))


      if(AID1 /= '      ON' .and. AID1 /= '     OFF')then
      call errors
      write(err,*)'w2_selective: DYNSTR1 CONTROL must be either ON or OFF. CONTROL=',aid1
      endif
      if(tcdfreq < 0.05)then
      call warnings
      write(wrn,*)'w2_selective: DYNSTR1 FREQ is < 0.05 days. This may be an unrealistically high frequency. FREQ=',tcdfreq
      endif

      do j=1,2
      read(1010,*)
      end do

    IF(AID1 == '     OFF')THEN
      do j=1,numtempc   ! skip lines
          read(1010,*,err=400)
      enddo
      do j=1,2*numtempc+4
          read(1010,*,err=400)
      enddo
    ELSE

      ncountc=0
      do j=1,numtempc
      if(char1=='$')then
      read(1010,*,err=400)I,tcntr(j),tcjb(j),tcjs(j),tcyearly,tctsrt,tctend,tctemp,tcnelev,(tcelev(n),n=1,tcnelev)
      tcntr(j)=ADJUSTR(TCNTR(J))
      tcyearly=ADJUSTR(tcyearly) 
      else
      read(1010,'(a8,a8,i8,i8,a8,f8.0,f8.0,f8.0,i8,10(f8.0))',err=400)aid1,tcntr(j),tcjb(j),tcjs(j),tcyearly,tctsrt,tctend,tctemp,tcnelev,(tcelev(n),n=1,tcnelev)
      endif

       if(tcntr(j) /='      ST' .and. tcntr(j) /='      WD' )then
        call errors
        write(err,*)'w2_selective: DYNSTR2 ST/WD must be either "ST" or "WD" for #',j,' ST/WD=',tcntr(j)
       endif
        if(tcntr(j) =='      ST')then
          if(NSTR(tcjb(j))==0)then
          call errors
          write(err,'(a,i3,a,i3)')'w2_selective: DYNSTR2 ST/WD is "ST" but no structures are defined. NSTR(JB)=',NSTR(tcjb(j)),' for JB=',tcjb(j)
          endif
        endif
          if(tcntr(j) =='      WD')then
          if(NWD==0)then
          call errors
          write(err,'(a,i3)')'w2_selective: DYNSTR2 ST/WD is "WD" but no withdrawals are defined. NWD=',NWD
          endif
       endif
       if(tcyearly /='      ON' .and. tcyearly/='     OFF' )then
        call errors
        write(err,'(a,i3,a,a)')'w2_selective: DYNSTR2 YEARLY must be either "ON" or "OFF" for #',j,' YEARLY=',tcyearly
       endif
       if(tcyearly /='      ON' .and. tctsrt > 366. )then
        call errors
        write(err,'(a,i3,a,a)')'w2_selective: DYNSTR2 YEARLY is "ON" but Start day is > JDAY 366. Either set YEARLY=OFF or set start and end dates with a JDAY less than 366.'
       endif
       if(tcjb(j) < 1 .or. tcjb(j) > nbr)then
        call errors
        write(err,'(a,i4,a)')'w2_selective: DYNSTR2 JB must be between 1 and # of branches. JB=',tcjb(j),' for #',j
       endif
       if(tcntr(j) == '      ST'.and. aid == '      ON')then
         if(tcjs(j) > NSTR(tcjb(j)))then
         call errors
         write(err,'(a,i3,a,i3,a,i3,a,i3)')'w2_selective: DYNSTR2 JS/NW must be <= # structures[NSTR(JB)=',nstr(tcjb(j)),']. JS/NW=',tcjs(j),' for #',j,' for JB=',tcjb(j)
         endif
       endif
        if(tcntr(j) == '      WD'.and. aid == '      ON')then
         if(tcjs(j) > NWD)then
         call errors
         write(err,'(a,i3,a,i3,a,i3)')'w2_selective: DYNSTR2 JS/NW must be <= # withdrawals[NWD=',nwd,']. JS/NW=',tcjs(j),' for #',j
         endif
       endif


       if(tctsrt < tmstrt  .and. tcyearly /= '      ON')then
         call errors
         write(err,'(a,f10.3,a,f10.3,a,i3)')'w2_selective: DYNSTR2 TSRT must be >= Time of model start [TMSTRT=',tmstrt,']. TSRT=',tctsrt,' for #',j
       endif
        if(tctsrt > tmend)then
         call errors
         write(err,'(a,f10.3,a,f10.3,a,i3)')'w2_selective: DYNSTR2 TSRT must be <= model end time [TMEND=',tmend,']. TSRT=',tctsrt,' for #',j
       endif
         if(tctend < tmstrt  .and. tcyearly /= '      ON')then
         call errors
         write(err,'(a,f10.3,a,f10.3,a,i3)')'w2_selective: DYNSTR2 TEND must be >= Time of model start [TMSTRT=',tmstrt,']. TEND=',tctend,' for #',j
       endif
        if(tctend > tmend)then
         call warnings   ! changed to a warning rather than an error
         write(wrn,'(a,f10.3,a,f10.3,a,i3)')'w2_selective: DYNSTR2 TEND must be <= model end time [TMEND=',tmend,']. TEND=',tctend,' for #',j
       endif


       if(tcnelev > 1)then
       do jj=1,tcnelev-1
        if(tcelev(jj) < tcelev(jj+1) .and. tctemp >= 0.0)then
        call ERRORS
        write(ERR,'(a,i2,a,f7.1,a,i2,a,f7.1,a,i3)')'w2_selective: DYNSTR2 ELEV(',jj,')=',tcelev(jj),' < ELEV(',jj+1,')=',tcelev(jj+1),' for #',j
        endif
        end do
       endif
      end do
      do j=1,2
      read(1010,*,err=400)
      end do
      do j=1,numtempc
        if(char1=='$')then
        read(1010,*,err=400)I,tciseg,tcklay,dynsel
        dynsel=ADJUSTR(dynsel) 
      else
        read(1010,'(a8,i8,f8.0,a8,f8.0,2i8)',err=400)aid1,tciseg,tcklay,dynsel   !aid1 is ignored - skips 8 spaces
      endif

        if(tciseg < 0)then
           IF (WDOC /= '      ON' .and. WDOC /= '   ONH' .and. WDOC /= '     ONS' )THEN
           CALL ERRORS
           WRITE (ERR,'(a,i3,a,i3)') 'w2_selective: Withdrawal outflow control [WDOC='//WDOC(6:8)//'] must be " ON (or OHN or ONS)" since TCISEG < 0 in w2_selective.npt: Tciseg=',tciseg
           END IF
           IF(NWDO < 1)THEN
           CALL ERRORS
           WRITE (ERR,'(a,i3,a,i3)') 'w2_selective: NWDO for Withdrawals [NWDO=',NWDO,'] must be > 1 since TCISEG < 0 in w2_selective.npt: Tciseg=',tciseg
           END IF
           IF(NIWDO < abs(tciseg))THEN
           CALL ERRORS
           WRITE (ERR,'(a,i3,a,i3)') 'w2_selective: NIWDO for Withdrawals [NIWDO=',INWDO,'] must be > ABS(TCISEG) since TCISEG<0 in w2_selective.npt: Tciseg=',tciseg
           END IF
           if(NWDO > 0)THEN
            IF(WDOF(1) >= TCDFREQ)THEN
            CALL ERRORS
            WRITE(ERR,'(A,F7.2,A,F7.2)')'w2_selective: Since TCISEG<0, WDOF >= TCDFREQ. This will lead to inaccuracies in checking temrpatures. Recommend WDOF<TCDFREQ. [TCDFREQ=',tcdfreq,' WDOF(1)=',wdof(1)
            ENDIF
           ENDIF
         endif
         IF (DYNSEL /= '      ON'  .and. DYNSEL /= '     OFF')THEN
           CALL ERRORS
           WRITE (ERR,'(a,i3,a,i3)') 'w2_selective: Dynamic selective withdrawal outflow control [DYNSEL='//DYNSEL//'] must be " ON" of "OFF" for control#:',j
         END IF
        IF (DYNSEL == '      ON')THEN
            WRITE (SEGNUM,'(I0)') J
            SEGNUM  = ADJUSTL(SEGNUM)
            L       = LEN_TRIM(SEGNUM)
            OPEN (5000,FILE='dynselective'//SEGNUM(1:L)//'.npt',STATUS='UNKNOWN',IOSTAT=IERR)
            IF (IERR /= 0) THEN
            WRITE (WIN,*) 'Could not open dynselective'//SEGNUM(1:L)//'.npt input file: DYNSEL=ON'
            CALL ERRORS
            WRITE(ERR,*)'Could not open dynselective'//SEGNUM(1:L)//'.npt input file: DYNSEL=ON'
            ENDIF
            CLOSE(5000)
        ENDIF
      end do

      do j=1,2
      read(1010,*,err=400)
      end do
      do j=1,numtempc
      if(char1=='$')then
      read(1010,*,err=400)AID1,tcelevcon
      tcelevcon=ADJUSTR(tcelevcon)
      else
      read(1010,'(a8,a8)',err=400)aid1,tcelevcon
      endif
         IF (tcelevcon /= '      ON'  .and. tcelevcon /= '     OFF')THEN
           CALL ERRORS
           WRITE (ERR,'(a,a,i3)') 'w2_selective: AUTO ELEVCONTROL='//tcelevcon//'] must be " ON" of "OFF" for control#:',j
         END IF
      end do

    ENDIF

    IF(SELECTC == '      ON')THEN
      do j=1,2
      read(1010,*)
      end do
      if(char1=='$')then
      read(1010,*,err=400)AID1,tspltc,numtsplt
      tspltc=ADJUSTR(tspltc)
      else
      read(1010,'(a8,a8,i8)',err=400)aid1,tspltc,numtsplt
      endif
      IF (tspltc /= '      ON'  .and. tspltc /= '     OFF')THEN
           CALL ERRORS
           WRITE (ERR,'(a)') 'w2_selective: SPLIT TEMP CNTR='//tspltc//'] must be " ON" of "OFF"'
         END IF
      do j=1,2
      read(1010,*)
      end do
      do j=1,numtsplt
      if(char1=='$')then
      read(1010,*,err=400)I,tspltcntr,tspltjb,tcyearly,tctsrt,tctend,tspltt,nouts,(jstsplt(n),n=1,nouts),elcont,dynsel
      tspltcntr=ADJUSTR(tspltcntr)
      tcyearly=ADJUSTR(tcyearly)
      elcont=ADJUSTR(elcont)
      dynsel=ADJUSTR(dynsel)
      else        
      read(1010,'(a8,a8,i8,a8,f8.0,f8.0,f8.0,i8,2i8,a8,a8)',err=400)aid1,tspltcntr,tspltjb,tcyearly,tctsrt,tctend,tspltt,nouts,(jstsplt(n),n=1,nouts),elcont,dynsel
      endif
        if(tcyearly /='      ON' .and. tctsrt > 366. )then
        call errors
        write(err,'(a,i3,a,a)')'w2_selective: SPLIT TEMP YEARLY is "ON" but Start day is > JDAY 366. Either set YEARLY=OFF or set start and end dates with a JDAY less than 366.'
       endif
IF (DYNSEL /= '      ON'  .and. DYNSEL /= '     OFF')THEN
           CALL ERRORS
           WRITE (ERR,'(a,i3,a,i3)') 'w2_selective: SPLIT TEMP dynamic temperature control [DYNSELSPLT='//DYNSEL//'] must be " ON" of "OFF" for split control#:',j
         END IF
if(tspltc == '      ON')then
        if(tspltcntr /='      ST' .and. tspltcntr /='      WD' )then
        call errors
        write(err,*)'w2_selective: SPLIT TEMP ST/WD must be either "ST" or "WD" for #',j,' ST/WD=',tspltcntr
        endif

        if(tspltcntr =='      ST')then
 !         if(NSTR(tcjb(j))==0)then
          if(NSTR(tspltjb)==0)then                ! cb 11/5/12
          call errors
          write(err,'(a,i3,a,i3)')'w2_selective: SPLIT TEMP ST/WD is "ST" but no structures are defined. NSTR(JB)=',NSTR(tspltjb),' for JB=',tspltjb
          endif
        endif
        if(tspltcntr =='      WD')then
          if(NWD==0)then
          call errors
          write(err,'(a,i3)')'w2_selective: SPLIT TEMP ST/WD is "WD" but no withdrawals are defined. NWD=,',NWD
          endif
       endif


       if(tcyearly /='      ON' .and. tcyearly/='     OFF' )then
        call errors
        write(err,'(a,i3,a,a8)')'w2_selective: SPLIT TEMP YEARLY must be either "ON" or "OFF" for #',j,' YEARLY=',tcyearly
       endif
      if(tctsrt < tmstrt  .and. tcyearly /= '      ON')then
         call errors
         write(err,'(a,f10.3,a,f10.3,a,i3)')'w2_selective: SPLIT TEMP TSRT must be >= Time of model start [TMSTRT=',tmstrt,']. TSRT=',tctsrt,' for #',j
       endif
        if(tctsrt > tmend)then
         call errors
         write(err,'(a,f10.3,a,f10.3,a,i3)')'w2_selective: SPLIT TEMP TSRT must be <= model end time [TMEND=',tmend,']. TSRT=',tctsrt,' for #',j
       endif
         if(tctend < tmstrt  .and. tcyearly /= '      ON')then
         call errors
         write(err,'(a,f10.3,a,f10.3,a,i3)')'w2_selective: SPLIT TEMP TEND must be >= Time of model start [TMSTRT=',tmstrt,']. TEND=',tctend,' for #',j
       endif
        if(tctend > tmend)then
         call errors
         write(err,'(a,f10.3,a,f10.3,a,i3)')'w2_selective: SPLIT TEMP TEND must be <= model end time [TMEND=',tmend,']. TEND=',tctend,' for #',j
        endif
         if(elcont /='      ON' .and. elcont /='     OFF' )then
           call errors
           write(err,'(a,i3,a,a8)')'w2_selective: SPLIT TEMP ELCONT must be either "ON" or "OFF" for #',j,' ELCONT=',elcont
         endif
         if(nouts > 2)then
             call errors
             write(err,'(a,i3,a,i8)')'w2_selective: SPLIT TEMP NOUTS must be 2 (not > 2) for #',j,' NOUTS=',nouts
         endif
       do n=1,nouts

       if(n <= nouts-1)then
         if(jstsplt(n) == jstsplt(n+1))then
         call errors
         write(err,'(a,i3,a,i3,a,i3,a,i3,a,i3)')'w2_selective: SPLIT TEMP JS1/NW1 and JS2/NW2 must be different outlets. JS1/NW1=',jstsplt(n),' JS2/NW2=',jstsplt(n+1),' for JS/NW#:',n,' and JS/NW#:',n+1,' for #',j
         endif
       endif

       if(tspltcntr == '      ST')then

       if(jstsplt(n) > nstr(tspltjb))then
         call errors
         write(err,'(a,i3,a,i3,a,i3)')'w2_selective: SPLIT TEMP JSn/NWn is > NSTR for that branch. JSn/NWn=',jstsplt(n),' NSTR(JB)=',nstr(tspltjb),' for #',j
       endif

       elseif(tspltcntr == '      WD')then

        if(jstsplt(n) > nwd)then
         call errors
         write(err,'(a,i3,a,i3,a,i3)')'w2_selective: SPLIT TEMP JSn/NWn is > NWD. JSn/NWn=',jstsplt(n),' NWD=',nwd,' for #',j
        endif

       endif
      enddo
 endif
      enddo
      do j=1,2
      read(1010,*)
      end do

      if(char1=='$')then
      read(1010,*,err=400)AID1,tempn
      else
      read(1010,'(a8,i8)',err=400)aid1,tempn
      endif

      do j=1,2
      read(1010,*)
      end do
      do j=1,tempn

      if(char1=='$')then
        read(1010,*,err=400)I,(tempcrit(jw),jw=1,nwb)
      else
        read(1010,'(a8,10f8.0)',err=400)aid1,(tempcrit(jw),jw=1,nwb)   ! Note max of 10 waterbodies
      endif

      enddo
    ! end check file

  ELSE

     do j=1,2
      read(1010,*)
      end do
      if(char1=='$')then
      read(1010,*,err=400)aid1,tspltc,numtsplt,tspltfreq, tsconv
      tspltc=adjustr(tspltc)
      else
      read(1010,'(a8,a8,i8,2f8.0)',err=400)aid1,tspltc,numtsplt,tspltfreq, tsconv
      endif
      IF (tspltc /= '      ON'  .and. tspltc /= '     OFF')THEN
           CALL ERRORS
           WRITE (ERR,'(a)') 'w2_selectiveUSGS: SPLIT TEMP CNTR='//tspltc//'] must be " ON" of "OFF"'
         END IF
      do j=1,2
      read(1010,*)
      end do
do j=1,numtsplt
      if(char1=='$')then
      read(1010,*,err=400)I,tspltcntr,tspltjb,tcyearly,tctsrt,tctend,tspltt,tsdynsel,elcontspl,nouts,tsshare
      tspltcntr=ADJUSTR(tspltcntr)
      tcyearly=ADJUSTR(tcyearly)
      elcontspl=ADJUSTR(elcont)
      tsdynsel=ADJUSTR(tsdynsel)
      tsshare=ADJUSTR(tsshare)
      else
      read(1010,'(a8,a8,i8,a8,f8.0,f8.0,f8.0,2a8,i8,a8)',err=400)aid1,tspltcntr,tspltjb,tcyearly,tctsrt,tctend,tspltt,tsdynsel,elcontspl,nouts,tsshare
      endif
if(tspltc == '      ON')then
        if(tspltcntr /='      ST' .and. tspltcntr /='      WD' )then
        call errors
        write(err,*)'w2_selectiveUSGS: SPLIT TEMP ST/WD must be either "ST" or "WD" for #',j,' ST/WD=',tspltcntr
        endif

        if(tspltcntr =='      ST')then
 !         if(NSTR(tcjb(j))==0)then
          if(NSTR(tspltjb)==0)then                ! cb 11/5/12
          call errors
          write(err,'(a,i3,a,i3)')'w2_selectiveUSGS: SPLIT TEMP ST/WD is "ST" but no structures are defined. NSTR(JB)=',NSTR(tspltjb),' for JB=',tspltjb
          endif
        endif
        if(tspltcntr =='      WD')then
          if(NWD==0)then
          call errors
          write(err,'(a,i3)')'w2_selectiveUSGS: SPLIT TEMP ST/WD is "WD" but no withdrawals are defined. NWD=,',NWD
          endif
       endif


       if(tcyearly /='      ON' .and. tcyearly/='     OFF' )then
        call errors
        write(err,'(a,i3,a,a8)')'w2_selectiveUSGS: SPLIT TEMP YEARLY must be either "ON" or "OFF" for #',j,' YEARLY=',tcyearly
       endif
      if(tctsrt < tmstrt  .and. tcyearly /= '      ON')then
         call errors
         write(err,'(a,f10.3,a,f10.3,a,i3)')'w2_selectiveUSGS: SPLIT TEMP TSRT must be >= Time of model start [TMSTRT=',tmstrt,']. TSRT=',tctsrt,' for #',j
       endif
        if(tctsrt > tmend)then
         call errors
         write(err,'(a,f10.3,a,f10.3,a,i3)')'w2_selectiveUSGS: SPLIT TEMP TSRT must be <= model end time [TMEND=',tmend,']. TSRT=',tctsrt,' for #',j
       endif
         if(tctend < tmstrt  .and. tcyearly /= '      ON')then
         call errors
         write(err,'(a,f10.3,a,f10.3,a,i3)')'w2_selectiveUSGS: SPLIT TEMP TEND must be >= Time of model start [TMSTRT=',tmstrt,']. TEND=',tctend,' for #',j
       endif
        if(tctend > tmend)then
         call errors
         write(err,'(a,f10.3,a,f10.3,a,i3)')'w2_selectiveUSGS: SPLIT TEMP TEND must be <= model end time [TMEND=',tmend,']. TEND=',tctend,' for #',j
        endif
         if(elcontspl /='      ON' .and. elcontspl /='     OFF' )then
           call errors
           write(err,'(a,i3,a,a8)')'w2_selectiveUSGS: SPLIT TEMP ELCONT must be either "ON" or "OFF" for #',j,' ELCONT=',elcontspl
         endif
        if(tsdynsel /='      ON' .and. tsdynsel /='     OFF' )then
           call errors
           write(err,'(a,i3,a,a8)')'w2_selectiveUSGS: SPLIT TEMP TSDYN must be either "ON" or "OFF" for #',j,' TSDYN=',tsdynsel
        endif

      if (nouts < 2) then
        CALL ERRORS
        write (err, '(A,I0)') 'w2_selectiveUSGS: Less than two outlets specified for blending group. NOUTS must be > 1. NOUTS=',nouts
      else if (nouts > 10) then
        CALL ERRORS
        write (err, '(A,I0)') 'w2_selectiveUSGS: More than ten outlets specified for blending group. NOUTS must be < 11. NOUTS=',nouts
      end if
  endif
  enddo

  read (1010,'(/)',err=400)
    do j=1,numtsplt
              if(char1=='$')then
  read (1010,*,err=400) AID1,(jstsplt(n),n=1,nouts)
                  ELSE
  read (1010,'(8x,10i8)',err=400) (jstsplt(n),n=1,nouts)
                  ENDIF
                  

    do n=1,nouts
       if(n <= nouts-1)then
         if(jstsplt(n) == jstsplt(n+1))then
         call errors
         write(err,'(a,i3,a,i3,a,i3,a,i3,a,i3)')'w2_selectiveUSGS: Duplicate split outlet numbers in group. JS1/NW1=',jstsplt(n),' JS2/NW2=',jstsplt(n+1),' for JS/NW#:',n,' and JS/NW#:',n+1,' for #',j
         endif
       endif
   enddo

    end do
  read (1010,'(/)',err=400)
    do j=1,numtsplt
                  if(char1=='$')then
        read (1010,*,err=400)AID1, (depth(n),n=1,nouts)
                  ELSE
  read (1010,'(8x,10f8.0)',err=400) (depth(n),n=1,nouts)
                  ENDIF
                  
    end do
  read (1010,'(/)',err=400)
    do j=1,numtsplt
                          if(char1=='$')then
                                read (1010,*,err=400) AID1,(minfrac(n),n=1,nouts)
                          ELSE
                                read (1010,'(8x,10f8.0)',err=400) (minfrac(n),n=1,nouts)
                          ENDIF    
    end do
  read (1010,'(/)',err=400)
    do j=1,numtsplt
                          if(char1=='$')then
                                read (1010,*,err=400) (priority(n),n=1,nouts)
                          ELSE
                                read (1010,'(8x,10i8)',err=400) (priority(n),n=1,nouts)
                          ENDIF    
  do n=1,nouts
  if (priority(n) < -1) then
    call errors
    write (err, '(A,I0,A,I0,A)') 'w2_selectiveUSGS: Priority input for outlet ', jstsplt(n), ' in group ', j, ' is less than -1.'
  end if
  enddo

    end do
  read (1010,'(/)',err=400)
    do j=1,numtsplt
                          IF(char1=='$')THEN
                                  read (1010,*,err=400) AID1,(minhead(n),n=1,nouts)
                          ELSE
                                  read (1010,'(8x,10f8.0)',err=400) (minhead(n),n=1,nouts)
                          ENDIF    
    end do
  read (1010,'(/)',err=400)
    do j=1,numtsplt
                          IF(char1=='$')THEN
                                  read (1010,*,err=400) AID1,(maxhead(n),n=1,nouts)
                          ELSE
                                  read (1010,'(8x,10f8.0)',err=400) (maxhead(n),n=1,nouts)
                          ENDIF          
    end do
  read (1010,'(/)',err=400)
    do j=1,numtsplt
                          IF(char1=='$')THEN
                                  read (1010,*,err=400) AID1,(maxflow(n),n=1,nouts)
                          ELSE
                                  read (1010,'(8x,10f8.0)',err=400) (maxflow(n),n=1,nouts)
                          ENDIF          
    end do

! NOTE: These other checks can only be implemented once we make the arrays above 2D rather than 1D SW 3.72 3/2015
!if (tspltcntr(j) .eq. tcntr(jj) .and. (tspltcntr(j) .eq. '      WD' .or. tspltjb(j) == tcjb(jj))) then
!              do n=1,nouts(j)
!                if (jstsplt(j,n) == tcjs(jj)) then
!                  write (w2err, '(A,I0,A)') 'ERROR-- Outlet number ',tcjs(jj),' used in tower and blending group at same time.'
!                  ERROR_OPEN = .TRUE.       ! will trigger the program to end when this subroutine is completed
!                end if
!              end do
!            end if
 !if (tsminhead(j,n) > 0.0 .and. tsmaxhead(j,n) > 0.0 .and. tsminhead(j,n) > tsmaxhead(j,n)) then
 !         write (wrn, '(A,I0,A,I0,A)') 'WARNING-- Minimum and maximum head constraints for outlet ', jstsplt(j,n), ' in group ',   &
 !                                       j, ' are such that the outlet cannot ever be used.'
 !         WARNING_OPEN = .TRUE.
 !       end if
 !       if (tsdepth(j,n) > 0.0 .and. tsminhead(j,n) > 0.0 .and. tsdepth(j,n) < tsminhead(j,n)) then
 !         write (wrn, '(A,I0,A,I0,A)') 'WARNING-- Depth of floating outlet ', jstsplt(j,n), ' in group ', j,                       &
 !             ' is shallower than the minimum head constraint.  To honor the head constraint, no flow is possible for that outlet.'
 !         WARNING_OPEN = .TRUE.
 !       end if
!if ((tstsrt(jj) >= tstsrt(j) .and. tstsrt(jj) <  tstend(j)) .or.                                                           &
!            (tstend(jj) >  tstsrt(j) .and. tstend(jj) <= tstend(j))) then
!          if (tspltcntr(j) .eq. tspltcntr(jj) .and. (tspltcntr(j) .eq. '      WD' .or. tspltjb(j) == tspltjb(jj))) then
!            do n=1,nouts(j)
!              do nj=1,nouts(jj)
!                if (jstsplt(j,n) == jstsplt(jj,nj)) then
!                  write (w2err, '(A,I0,A)') 'ERROR-- Split outlet number ', jstsplt(j,n), ' used in more than one group at a time.'
!                  ERROR_OPEN = .TRUE.       ! will trigger the program to end when this subroutine is completed
!                end if
!              end do
!            end do
!          end if
!        end if



      do j=1,2
      read(1010,*)
      end do

      if(char1=='$')then
      read(1010,*,err=400)AID1,tempn
      else
      read(1010,'(a8,i8)',err=400)aid1,tempn
      endif

      do j=1,2
      read(1010,*)
      end do
      do j=1,tempn

      if(char1=='$')then
        read(1010,*,err=400)I,(tempcrit(jw),jw=1,nwb)
      else
        read(1010,'(a8,10f8.0)',err=400)aid1,(tempcrit(jw),jw=1,nwb)   ! Note max of 10 waterbodies
      endif

      enddo
    ! end check file
  END IF
  ENDIF
  CLOSE(1010)
ENDIF

! HABTATC CONTROL
IF(HABTATC /= '     OFF' .AND. HABTATC /= '      ON')THEN
 CALL ERRORS
 WRITE(ERR,*)'HABTATC MUST BE "OFF" OR "ON". HABTATC=',HABTATC
ENDIF

IF(HABTATC == '      ON')THEN      ! CHECK FILE HABTATC
    unit=5000
    WRITE (WIN,*) 'Reading W2 Habitat File: w2_habitat.npt'       
   OPEN (5000,FILE='w2_habitat.npt',STATUS='UNKNOWN',IOSTAT=IERR)
  IF (IERR /= 0) THEN
    WRITE (WIN,*) 'Could not open w2_habitat.npt input file'
    CALL ERRORS
    WRITE(ERR,*)'Could not open w2_habitat.npt input file'
  ELSE
    ! check file
    AID='HABITAT'
     WRITE(WIN,*)'Reading w2_habitat.npt'
        read(5000,*,err=400)
        read(5000,*,err=400)
        read(5000,*,err=400)ifish,conhab
        IF(IFISH == 0)IFISH=1
        allocate (fishname(ifish),fishtemph(ifish),fishtempl(ifish),fishdo(ifish),habvol(ifish),phabvol(ifish))
        read(5000,*,err=400)
        do i=1,ifish
          read(5000,*,err=400)fishname(i),fishtempl(i),fishtemph(i),fishdo(i)
          IF(FISHTEMPH(I) > 35.0)THEN
          CALL WARNINGS
          WRITE(WRN,'(a,i3,a,f6.2)')'w2_habitat: Fish Temperature HIGH for habitat volume > 35C! for Fish#',i,' Temp=',fishtemph(i)
          endif
          IF(FISHTEMPL(I) > 30.0)THEN
          CALL WARNINGS
          WRITE(WRN,'(a,i3,a,f6.2)')'w2_habitat: Fish Temperature LOW for habitat volume > 30C! for Fish#',i,' Temp=',fishtempl(i)
          endif
          IF(FISHDO(I) < 1.5)THEN
          CALL WARNINGS
          WRITE(WRN,'(a,i3,a,f6.2)')'w2_habitat: Fish Dissolved Oxygen for habitat volume < 1.5 mg/l! for Fish#',i,' DO=',fishdo(i)
          endif
        enddo
        read(5000,*,err=400)
        read(5000,*,err=400)nseg,conavg     ! volume weighted averages of critical WQ parameters
        if(nseg>imx-2)then
        CALL ERRORS
        WRITE(ERR,FMTI)'w2_habitat: # of segments to compute volume averages [NSEG=',NSEG,'] is > IMX-2=',imx-2
        endif
        read(5000,*,err=400)
        allocate(isegvol(nseg),cdo(nseg),cpo4(nseg),cno3(nseg),cnh4(nseg),cchla(nseg),ctotp(nseg),cdos(nseg),cpo4s(nseg),cno3s(nseg),cnh4s(nseg),cchlas(nseg),ctotps(nseg),cgamma(nseg))
        allocate(ssedd(imx))
        read(5000,*,err=400)(isegvol(i),i=1,nseg)   ! check to make sure this i is a real segment****

        DO I=1,NSEG
        ACTIVE_SEGMENT = .FALSE.
        DO JB=1,NBR
          IF (ISEGVOL(I) >= US(JB) .AND. ISEGVOL(I) <= DS(JB)) ACTIVE_SEGMENT = .TRUE.
        END DO
        IF (.NOT.ACTIVE_SEGMENT) THEN
          CALL ERRORS
          WRITE (ERR,FMTI) 'w2_habitat: Segment for weighted averages [ISEGVOL=',ISEGVOL(I),'] is a boundary segment or not in the active grid for # ',I
        END IF
        ENDDO

        read(5000,*,err=400)
        read(5000,*,err=400)kseg,consurf     ! # of layers for surface average

        IF(KSEG > KMX-2)THEN
        CALL ERRORS
        WRITE(ERR,FMTI)'w2_habitat: # of vertical layers to average [KSEG=',KSEG,'] is greater than KMX-2=',kmx-2
        END IF

        read(5000,*,err=400)
        read(5000,*,err=400)consod
  END IF
  CLOSE(5000)
ENDIF

! ENVIRPC CHECK
    AID='ENVIRPRF'
IF(ENVIRPC /= '     OFF' .AND. ENVIRPC /= '      ON')THEN
 CALL ERRORS
 WRITE(ERR,*)'ENVIRPC MUST BE "OFF" OR "ON". ENVIRPC=',ENVIRPC
ENDIF

IF(ENVIRPC == '      ON')THEN      ! CHECK FILE ENVIRPC
     WRITE (WIN,*) 'Reading W2 Environmental Performance File: w2_envirprf.npt'  
unit=5000
   OPEN (5000,FILE='w2_envirprf.npt',STATUS='UNKNOWN',IOSTAT=IERR)
  IF (IERR /= 0) THEN
    WRITE (WIN,*) 'Could not open w2_envirprf.npt input file'
    CALL ERRORS
    WRITE(ERR,*)'Could not open w2_envirprf.npt input file'
  ELSE
    ! check file
      N=0
     ALLOCATE(cc_e(NCT),c_int(NCT),c_top(NCT),cd_e(NDC),cd_int(NDC),cd_top(NDC),c_avg(NCT),cd_avg(NDC),cn_e(NCT),cdn_e(NDC))
     CSVFORMAT=.FALSE.
     READ(UNIT,'(//A)')HEADER(1)
     DO J=1,30
         IF(HEADER(1)(J:J)==',')THEN
             CSVFORMAT=.TRUE.
             EXIT
         ENDIF
     ENDDO
     REWIND(UNIT)
      
      IF(CSVFORMAT)THEN
          READ(UNIT,*)
          READ(UNIT,*)
          READ (UNIT,*,ERR=400,END=5001) i,numclass,selectivec,sjday1,sjday2,istart1,iend1  ! NEED TO START AT 2ND FIELD FOR NUMCLASS
          SELECTIVEC=ADJUSTR(SELECTIVEC)
      ELSE
          READ (unit,1200,ERR=400,END=5001) i,numclass,selectivec,sjday1,sjday2,istart1,iend1
      ENDIF
      N=1
      IF(NUMCLASS>50)THEN
        CALL WARNINGS
        WRITE(WRN,*)'ENVIRPC: NUMCLASS > 50. That may be a lot of intervals. NUMCLASS=',NUMCLASS
      ENDIF
       IF(selectivec /= 'OFF' .AND. selectivec /= ' ON')THEN
        CALL ERRORS
        WRITE(ERR,*)'ENVIRPC: SELECTC MUST BE "OFF" OR "ON". SELECTC=',selectivec
       ENDIF
       IF(istart1>imx)THEN
        CALL ERRORS
        WRITE(ERR,'(a45,i4,a6,i4)')'ENVIRPC: ISTR1 MUST BE LESS THAN IMX. ISTART1=',ISTART1,'>IMX=',IMX
       ENDIF
        IF(iend1>imx)THEN
        CALL ERRORS
        WRITE(ERR,'(a45,i4,a6,i4)')'ENVIRPC: IEND1 MUST BE LESS THAN IMX. IEND1=',IEND1,'>IMX=',IMX
       ENDIF

       IF(istart1>iend1)THEN
        CALL ERRORS
        WRITE(ERR,'(a45,i4,a6,i4)')'ENVIRPC: ISTR1 MUST BE LESS THAN IEND1. ISTART1=',ISTART1,'>IEND1=',IEND1
       ENDIF
       IF(sjday1>tmend)THEN
        CALL ERRORS
        WRITE(ERR,*)'ENVIRPC: JDAY1 MUST BE LESS THAN TMEND (End of model run time). JDAY1=',SJDAY1
       ENDIF
        IF(sjday2>tmend)THEN
        CALL ERRORS
        WRITE(ERR,*)'ENVIRPC: JDAY2 MUST BE LESS THAN TMEND (End of model run time). JDAY1=',SJDAY2
       ENDIF
       IF(sjday1>sjday2)THEN
        CALL ERRORS
        WRITE(ERR,*)'ENVIRPC: JDAY1 > JDAY2. JDAY1 must be less than JDAY2. JDAY1=',SJDAY1
       ENDIF
      IF(CSVFORMAT)THEN
          READ(UNIT,*)
          READ(UNIT,*)
          Read (unit,*,ERR=400,END=5001) VEL_VPR, VEL_INT, VEL_TOP,TEMP_VPR,TEMP_INT,TEMP_TOP, depth_vpr,d_int,d_top
          VEL_VPR=ADJUSTR(VEL_VPR);TEMP_VPR=ADJUSTR(TEMP_VPR);DEPTH_VPR=ADJUSTR(DEPTH_VPR)
      ELSE
          Read (unit,1201,ERR=400,END=5001) VEL_VPR, VEL_INT, VEL_TOP,TEMP_VPR,TEMP_INT,TEMP_TOP, depth_vpr,d_int,d_top
      ENDIF
      N=2
       IF(VEL_VPR /= 'OFF' .AND. VEL_VPR /= ' ON')THEN
        CALL ERRORS
        WRITE(ERR,*)'ENVIRPC: VELOCITY CONTROL MUST BE "OFF" OR "ON". VELOCITY=',VEL_VPR
       ENDIF
      IF(TEMP_VPR /= 'OFF' .AND. TEMP_VPR /= ' ON')THEN
        CALL ERRORS
        WRITE(ERR,*)'ENVIRPC: TEMPC MUST BE "OFF" OR "ON". TEMPC=',TEMP_VPR
      ENDIF
       IF(DEPTH_VPR /= 'OFF' .AND. DEPTH_VPR /= ' ON')THEN
        CALL ERRORS
        WRITE(ERR,*)'ENVIRPC: DEPTHC MUST BE "OFF" OR "ON". DEPTHC=',DEPTH_VPR
       ENDIF
      IF(CSVFORMAT)THEN
          READ(UNIT,*)
          READ(UNIT,*)
          DO JC=1,NCT
            READ (unit,*,ERR=400,END=5001) AID1,CC_E(JC), C_INT(JC), C_TOP(JC)
            CC_E(JC)=ADJUSTR(CC_E(JC))
          ENDDO
          READ(UNIT,*)
          READ(UNIT,*)
          DO JD=1,NDC
            READ (unit,*,ERR=400,END=5001) AID1,CD_E(JD),CD_INT(JD), CD_TOP(JD)
            CD_E(JD)=ADJUSTR(CD_E(JD))
          ENDDO
          
      ELSE

      READ (unit,1050,ERR=400,END=5001) (CC_E(JC), C_INT(JC), C_TOP(JC), JC=1,NCT)
      N=3          ! ADD MORE FILE CHECKS HERE
      READ (unit,1050,ERR=400,END=5001) (CD_E(JD),CD_INT(JD), CD_TOP(JD), JD=1,NDC)
      ENDIF
      CLOSE(unit)
      
 1050 FORMAT(//(8X,(5X,A3,F8.0,F8.0)))
 1200 FORMAT(//i1,7X,I8,5x,a3,f8.0,f8.0,i8,i8)
 1201 format(//3(5x,a3,f8.3,f8.3))
      GO TO 5002
5001 CALL ERRORS
    WRITE(ERR,'(A125,I4)')'ENVIRPC: w2_envirprf.npt is not complete. End of file on reading. Of the 4 sections in the file, the last valid section was #',N
5002 CONTINUE      
      
  END IF
  
ENDIF
! CHECK ON/OFF FOR WL, FLOWBAL AND NPBAL
IF(WLC /= '     OFF' .AND. WLC /= '      ON')THEN
 CALL ERRORS
 WRITE(ERR,*)'WLC (Water Level output) MUST BE "OFF" OR "ON". WLC=',WLC
ENDIF
IF(WLC=='      ON')THEN
    IF(WLF<= 0.0)THEN
        CALL ERRORS
        WRITE(ERR,*)'WLF (Water Level output frequency) is <= 0.0. WLF=',WLF
    ENDIF
ENDIF
IF(FLOWBALC /= '     OFF' .AND. FLOWBALC /= '      ON')THEN
 CALL ERRORS
 WRITE(ERR,*)'FLOWBALC (Flow balance output) MUST BE "OFF" OR "ON". FLOWBALC=',FLOWBALC
ENDIF
IF(FLOWBALC == '      ON' .AND. VBC(1) /= '      ON')THEN
 CALL WARNINGS
 WRITE(WRN,'(A)')'FLOWBALC (Flow balance output) is ON but the Volume Balance is OFF. Hence no flow balance results will be output. Suggest turn Volume Balance (VBC) ON'
ENDIF

IF(FLOWBALC=='      ON')THEN
    IF(FLOWBALF<= 0.0)THEN
        CALL ERRORS
        WRITE(ERR,*)'FLOWBALF (Flow balance output frequency) is <= 0.0. FLOWBALF=',FLOWBALF
    ENDIF
ENDIF
IF(NPBALC /= '     OFF' .AND. NPBALC /= '      ON')THEN
 CALL ERRORS
 WRITE(ERR,*)'NPBALC (N&P mass balance output) MUST BE "OFF" OR "ON". NPBALC=',NPBALC
ENDIF
IF(NPBALC=='      ON')THEN
    IF(NPBALF<= 0.0)THEN
        CALL ERRORS
        WRITE(ERR,*)'NPBALF (N&P mass balance output frequency) is <= 0.0. NPBALF=',NPBALF
    ENDIF
ENDIF
IF(NPBALC=='      ON')THEN
    DO JW=1,NWB
    IF(MBC(JW) /='      ON')THEN
        CALL ERRORS
        WRITE(ERR,'(A,A,A,I5)')'NPBALF (N&P mass balance output) must also have MBC=ON for each water body. MBC(JW)=',MBC(JW),' for WB:',JW
    ENDIF
    ENDDO
ENDIF

! AERATEC CHECK

IF(AERATEC /= '     OFF' .AND. AERATEC /= '      ON')THEN
 CALL ERRORS
 WRITE(ERR,*)'AERATEC MUST BE "OFF" OR "ON". AERATEC=',AERATEC
ENDIF
CSVFORMAT=.FALSE.
IF(AERATEC == '      ON')THEN      ! CHECK FILE AERATEC
     WRITE (WIN,*) 'Reading W2 Aeration File: w2_aerate.npt'  
    unit=5000
       OPEN (5000,FILE='w2_aerate.npt',STATUS='UNKNOWN',IOSTAT=IERR)
  IF (IERR /= 0) THEN
    WRITE (WIN,*) 'Could not open w2_aerate.npt input file'
    CALL ERRORS
    WRITE(ERR,*)'Could not open w2_aerate.npt input file'
  ELSE
  ! check file
     open(5000,file='w2_aerate.npt',status='old')
   READ(5000,'(//A)')HEADER(1)
     DO J=1,30
         IF(HEADER(1)(J:J)==',')THEN
             CSVFORMAT=.TRUE.
             EXIT
         ENDIF
     ENDDO
    REWIND(5000)
    
     IF(CSVFORMAT)THEN
            READ (5000,*)
            READ (5000,*)
            READ (5000,*,err=400)NAER,conaer
         
         ELSE
     
           READ (5000,'(//i8,a16)',err=400)NAER,conaer
         ENDIF
         
     if(naer.eq.0)naer=1
     read(5000,1013)
1013 FORMAT(/)
      do i=1,naer
               IF(CSVFORMAT)THEN
               Read(5000,*,err=400)iaseg,ktopa,kbota,smass,atimon,atimoff,dzmult,dooff,doon,iprb,kprb
               ELSE
               Read(5000,'(i8,i8,i8,f8.0,f8.0,f8.0,3f8.0,2i8)',err=400)iaseg,ktopa,kbota,smass,atimon,atimoff,dzmult,dooff,doon,iprb,kprb
               ENDIF
               
        if(ktopa >= kmx  .or. ktopa <= 1)then
        CALL ERRORS
        WRITE(ERR,'(A,I4)')'w2_aerate: KTOPA >= KMX or <= 1. KTOPA is not in the active grid area for aerator #',I
        endif
        if(kbota >= kmx  .or. kbota <= 1)then
        CALL ERRORS
        WRITE(ERR,'(A,I4)')'w2_aerate: KBOPA >= KMX or <= 1. KBOPA is not in the active grid area for aerator #',I
        endif
        if(atimon >= TMEND  .or. atimon < TMSTRT)then
        CALL ERRORS
        WRITE(ERR,'(A,I4)')'w2_aerate: ATIMON is >= TMEND or < TMSTRT for aerator #',I
        endif
        if(atimon > atimoff)then
        CALL ERRORS
        WRITE(ERR,'(A,I4)')'w2_aerate: ATIMON is > ATIMOFF for aerator #',I
        endif
        if(atimoff > TMEND  .or. atimoff < TMSTRT)then
        CALL ERRORS
        WRITE(ERR,'(A,I4)')'w2_aerate: ATIMOFF is > TMEND or < TMSTRT for aerator #',I
        endif
        if(DZMULT < 1.0)then
        CALL ERRORS
        WRITE(ERR,'(A,I4)')'w2_aerate: DZMULT must be 1.0 or higher. DZMULT < 1.0 for aerator #',I
        endif
        if(dooff <  2.0)then
        CALL WARNINGS
        WRITE(WRN,'(A,I4)')'w2_aerate: [DOOFF] DO for turning OFF the aerator is < 2 mg/l for aerator #',I
        endif
        if(doon > 10.0)then
        CALL WARNINGS
        WRITE(WRN,'(A,I4)')'w2_aerate: [DOON] DO for turning ON the aerator is > 10 mg/l for aerator #',I
        endif
        if(kprb >= kmx  .or. kprb <= 1)then
        CALL ERRORS
        WRITE(ERR,'(A,I4)')'w2_aerate: KPRB >= KMX or <= 1. KPRB is not in the active grid area for aerator #',I
        endif
        ACTIVE_SEGMENT = .FALSE.
        DO JB=1,NBR
          IF (IASEG >= US(JB) .AND. IASEG <= DS(JB)) ACTIVE_SEGMENT = .TRUE.
        END DO
        IF (.NOT.ACTIVE_SEGMENT) THEN
          CALL ERRORS
          WRITE (ERR,FMTI) 'w2_aerate: Aeration segment [IASEG=',IASEG,'] is a boundary segment or not in the active grid for aerator # ',I
        END IF
        ACTIVE_SEGMENT = .FALSE.
        DO JB=1,NBR
          IF (IPRB >= US(JB) .AND. IPRB <= DS(JB)) ACTIVE_SEGMENT = .TRUE.
        END DO
        IF (.NOT.ACTIVE_SEGMENT) THEN
          CALL ERRORS
          WRITE (ERR,FMTI) 'w2_aerate: Aeration probe segment [IPRB=',IPRB,'] is a boundary segment or not in the active grid for aerator # ',I
        END IF

    ENDDO

!      end do

  END IF
  CLOSE(5000)
ENDIF

! INITUWL CHECK

IF(INITUWL /= '     OFF' .AND. INITUWL /= '      ON')THEN
 CALL ERRORS
 WRITE(ERR,*)'INITUWL MUST BE "OFF" OR "ON". INITUWL=',INITUWL
ENDIF

! Sediment Diagenesis Check
IF(SED_DIAG=='      ON')THEN

 INQUIRE(FILE="W2_CEMA_Input.npt", EXIST=DIAG_EXISTS)   ! file_exists will be TRUE if the file
    
 If(DIAG_EXISTS) Then	
     CALL ERRORS
     WRITE(ERR,*)'W2_CEMA_Input.npt file name must be changed to w2_diagenesis.npt for Version 4.2 and higher'
 endif
 
 INQUIRE(FILE="W2_diagenesis.npt", EXIST=DIAG_EXISTS)   ! file_exists will be TRUE if the file
  If(DIAG_EXISTS) Then	
    WRITE (WIN,*) 'Reading W2 Sediment Diagenesis File: w2_diagenesis.npt'  
    Open(5000, File = "W2_diagenesis.npt")
	
	!Read Header
	SkipLoop = .FALSE.
	Do While(.NOT. SkipLoop)
	  Read(5000,'(a)')MessageTemp
      If(index(MessageTemp, "$") == 0)SkipLoop = .TRUE.
    End Do

    BackSpace(5000)
    CEMAFilN=5000
    
    Read(CEMAFilN,*,err=6500)MessageTemp, SD_global
    IF(.NOT.SD_GLOBAL)GO TO 6501
    ! GROUP 2: FFT Layer
    Read(CEMAFilN,*,err=6500)MessageTemp, IncludeFFTLayer
    IF(IncludeFFTLayer) THEN
        FirstTimeInFFTCode = .TRUE.
        Read(CEMAFilN,*,err=6500)MessageTemp, NumFFTActivePrds
        Allocate(FFTActPrdSt(NumFFTActivePrds), FFTActPrdEn(NumFFTActivePrds))
        Read(CEMAFilN,*,err=6500)MessageTemp, (FFTActPrdSt(i), i = 1, NumFFTActivePrds) 
        Read(CEMAFilN,*,err=6500)MessageTemp, (FFTActPrdEn(i), i = 1, NumFFTActivePrds)
        Read(CEMAFilN,*,err=6500)MessageTemp, InitFFTLayerConc
        Read(CEMAFilN,*,err=6500)MessageTemp, FFTLayerSettVel
        Read(CEMAFilN,*,err=6500)MessageTemp, MoveFFTLayerDown
        
    IF (IncludeFFTLayer .and. NumFFTActivePrds <= 0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Sediment Diagenesis:Fine Fluids Tailing (FFT) Layer is ON but number of time periods ',        NumFFTActivePrds,' <= 0 '
    END IF

    do i=1,NumFFTActivePrds
      if(FFTActPrdEn(i) < FFTActPrdSt(i))then
        CALL ERRORS
        WRITE (ERR,FMTI) 'Sediment Diagenesis:Fine Fluids Tailing (FFT) Layer period end time ',        FFTActPrdEn(i),' < start time ',FFTActPrdSt(i)        
      end if
    end do

    IF (InitFFTLayerConc < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTF) 'Sediment Diagenesis:Initial tailings concentration in FFT (gm/m^3) ',        InitFFTLayerConc,' < 0.0 '
    END IF
    
    IF (InitFFTLayerConc == 0.0) THEN
      CALL WARNINGS
      WRITE (WRN,FMTF) 'Sediment Diagenesis:Initial tailings concentration in FFT (gm/m^3) ',        InitFFTLayerConc,' = 0.0 '
    END IF


    IF (FFTLayerSettVel < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTF) 'Sediment Diagenesis:Settling velocity of FFT to MFT (m/d) ',        FFTLayerSettVel,' < 0.0 '
    END IF
    
    IF (FFTLayerSettVel == 0.0) THEN
      CALL WARNINGS
      WRITE (WRN,FMTF) 'Sediment Diagenesis:Settling velocity of FFT to MFT (m/d) ',        FFTLayerSettVel,' = 0.0 '
    END IF


    ELSE
        DO JSKIP=1,6
            READ(CEMAFilN,'(A80)',err=6500)MessageTemp
        ENDDO
    END IF
    !
    ! GROUP 3: Bed Consolidation
    Read(CEMAFilN,*,err=6500)MessageTemp, IncludeBedConsolidation
    IF(IncludeBedConsolidation) THEN 
	    Read(CEMAFilN,*,err=6500)MessageTemp, LayerAddThkFrac
	    Read(CEMAFilN,*,err=6500)MessageTemp, NumConsolidRegns
	    Allocate(ConsolidationType(NumConsolidRegns),ConstConsolidRate(NumConsolidRegns))
	    Allocate(ConstPoreWtrRate(NumConsolidRegns),ConsolidRateTemp(NumConsolidRegns))
	    Allocate(ConsRegSegSt(NumConsolidRegns), ConsRegSegEn(NumConsolidRegns))
	    Read(CEMAFilN,*,err=6500)MessageTemp, (ConsRegSegSt(i), i = 1, NumConsolidRegns)
	    Read(CEMAFilN,*,err=6500)MessageTemp, (ConsRegSegEn(i), i = 1, NumConsolidRegns)
	    Read(CEMAFilN,*,err=6500)MessageTemp, (ConsolidationType(i), i = 1, NumConsolidRegns)
	    Read(CEMAFilN,*,err=6500)MessageTemp, (ConstConsolidRate(i), i = 1, NumConsolidRegns)
	    Read(CEMAFilN,*,err=6500)MessageTemp, ConsolidRateRegnFil
        Read(CEMAFilN,*,err=6500)MessageTemp, WriteBESnp
        Read(CEMAFilN,*,err=6500)MessageTemp, WritePWSnp
    if(IncludeBedConsolidation  .and. NumConsolidRegns == 0)then     ! CB 7/9/2018 changed .or. to .and.
        CALL ERRORS
        WRITE(ERR,*)'Sediment Diagenesis:Bed Consolidation is turned on but no regions specified'
    endif

    do i=1,NumConsolidRegns
      if(ConsRegSegSt(i) < 1 .or. ConsRegSegSt(i) > IMX)then
          CALL ERRORS
          WRITE(ERR,FMTI)'Sediment Diagenesis:Starting segment in bed consolidation region outside model grid: I=',ConsRegSegSt(i),' is less than 1 or greater than IMX=',imx
      endif
    
      if(ConsRegSegEn(i) < 1 .or. ConsRegSegEn(i) > IMX)then
          CALL ERRORS
          WRITE(ERR,FMTI)'Sediment Diagenesis:Ending segment in bed consolidation region outside model grid: I=',ConsRegSegEn(i),' is less than 1 or greater than IMX=',imx
      endif
      if(ConsRegSegSt(i) > ConsRegSegEn(i))then
          CALL ERRORS
          WRITE(ERR,FMTI)'Sediment Diagenesis:Starting segment in bed consolidation region  I=',ConsRegSegSt(i),' is greater than ending segment I=',ConsRegSegEn(i)
      endif

      if(ConsolidationType(i) /= 0 .and. ConsolidationType(i) /= 1)then
          CALL ERRORS
          WRITE(ERR,FMTI)'Sediment Diagenesis:Data type for bed consolidation for each region must be 0 or 1 (0: Constant, 1: Time varying)'
      endif

        IF (ConstConsolidRate(i) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Sediment Diagenesis:Bed Consolidation Rate ',        ConstConsolidRate(i),' < 0.0 for region ',    i
        ELSE IF (ConstConsolidRate(i) > 1.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Sediment Diagenesis:Bed Consolidation Rate ',        ConstConsolidRate(i),' > 1.0 m/d for region ',    i
        END IF
    end do 
    ELSE
        DO JSKIP=1,9
            READ(CEMAFilN,*,err=6500)
        ENDDO
    ENDIF
    ! GROUP 4: Sediment Diagenesis
    Read(CEMAFilN,*,err=6500)MessageTemp, IncludeCEMASedDiagenesis
    IF(IncludeCEMASedDiagenesis) THEN
        SEDIMENT_CALC = .FALSE.
    END IF
    Read(CEMAFilN,*,err=6500)MessageTemp, BedElevationInit
	Read(CEMAFilN,*,err=6500)MessageTemp, BedPorosityInit
    Read(CEMAFilN,*,err=6500)MessageTemp, CEMAParticleSize
    Read(CEMAFilN,*,err=6500)MessageTemp, CEMASedimentType
	Read(CEMAFilN,*,err=6500)MessageTemp, CEMASedimentDensity
	Read(CEMAFilN,*,err=6500)MessageTemp, CEMASedimentSVelocity
    Read(CEMAFilN,*,err=6500)MessageTemp, CEMASedimentProcessesInc
    
        IF (BedElevationInit <= 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Sediment Diagenesis:Initial sediment diagenesis bed thickness',bedelevationinit,' <= 0'
    END IF

    IF (BedPorosityInit < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Sediment Diagenesis:Initial sediment bed porosity ',        BedPorosityInit,' < 0.0 '
    ELSE IF (BedPorosityInit > 1.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Sediment Diagenesis:Initial sediment bed porosity ',        BedPorosityInit,' > 1.0  '
    END IF    

    IF (CEMAParticleSize <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Sediment Diagenesis:Sediment particle size ',        CEMAParticleSize,' <= 0.0 '
    END IF

    if(CEMASedimentType /= 1 .and. CEMASedimentType /= 2)then
          CALL ERRORS
          WRITE(ERR,FMTI)'Sediment Diagenesis:Sediment type must be 1: Cohesive or 2: Non-cohesive'
      endif

    IF (CEMASedimentDensity <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Sediment Diagenesis:Sediment bulk density ',        CEMAParticleSize,' <= 0.0 '
    END IF

    IF (CEMASedimentSVelocity < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Sediment Diagenesis:Sediment particle settling velocity ',        CEMASedimentSVelocity,' < 0.0 '
    END IF
    
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
    !
    !
    IF(IncludeCEMASedDiagenesis) THEN
        sediment_diagenesis=.true.
        Read(CEMAFilN,*,err=6500)MessageTemp, Bubbles_Calculation
        
            ! GROUP 5: Bubbles
    !IF (.NOT. IncludeCEMASedDiagenesis) Bubbles_Calculation = .FALSE.
    IF(Bubbles_Calculation) THEN    
        Read(CEMAFilN,*,err=6500)MessageTemp, GasDiff_Sed    ! in m^2/s
        Read(CEMAFilN,*,err=6500)MessageTemp, CalibParam_R1
        Read(CEMAFilN,*,err=6500)MessageTemp, YoungModulus
        Read(CEMAFilN,*,err=6500)MessageTemp, CritStressIF
        Read(CEMAFilN,*,err=6500)MessageTemp, BubbRelScale
        Read(CEMAFilN,*,err=6500)MessageTemp, CrackCloseFraction
        Read(CEMAFilN,*,err=6500)MessageTemp, LimBubbSize
        Read(CEMAFilN,*,err=6500)MessageTemp, MaxBubbRad
        Read(CEMAFilN,*,err=6500)MessageTemp, UseReleaseFraction
        Read(CEMAFilN,*,err=6500)MessageTemp, BubbRelFraction
        Read(CEMAFilN,*,err=6500)MessageTemp, BubbAccFraction
        Read(CEMAFilN,*,err=6500)MessageTemp, NumBubRelArr
        Read(CEMAFilN,*,err=6500)MessageTemp, BubbRelFractionAtm
        Read(CEMAFilN,*,err=6500)MessageTemp, BubbWatGasExchRate
        Read(CEMAFilN,*,err=6500)MessageTemp, ApplyBubbTurb
        Read(CEMAFilN,*,err=6500)MessageTemp, CEMATurbulenceScaling
        
      IF (  GasDiff_Sed <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTF) 'Sediment Diagenesis:Gas diffusion coefficient in sediment in m�/s is',  GasDiff_Sed,', which is <= 0.0'
      end if  

     IF ( CalibParam_R1 <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTF) 'Sediment Diagenesis:Calibration parameter R1 in m is',  CalibParam_R1,', which is <= 0.0'
    end if

    IF ( YoungModulus <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTF) 'Sediment Diagenesis:Youngs modulus E in N/m^2 is',  YoungModulus,', which is <= 0.0'
    end if

    IF ( CritStressIF <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTF) 'Sediment Diagenesis:Critical stress intensity factor for sediments K1c in N/m^3/2 is',CritStressIF,', which is <= 0.0'
    end if

    IF ( BubbRelScale <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTF) 'Sediment Diagenesis:Bubbles release scale is',BubbRelScale,', which is <= 0.0'
    end if

    IF ( CrackCloseFraction <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTF) 'Sediment Diagenesis:Fraction of critical pressure at which cracks close is',CrackCloseFraction,', which is <= 0.0'
    end if

    IF ( MaxBubbRad <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTF) 'Sediment Diagenesis:Maximum bubble radius in mm is',MaxBubbRad,', which is <= 0.0'
    end if

    IF ( BubbRelFraction <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTF) 'Sediment Diagenesis:Bubbles release fraction (sediments) is',BubbRelFraction,', which is <= 0.0'
    end if

    IF (BubbAccFraction <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTF) 'Sediment Diagenesis:Bubbles accumulation fraction (sediments) is',BubbAccFraction,', which is <= 0.0'
    end if

    IF (NumBubRelArr <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTF) 'Sediment Diagenesis:Number of bubbles release array is',NumBubRelArr,', which is <= 0.0'
    end if

    IF (BubbRelFractionAtm <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTF) 'Sediment Diagenesis:Bubbles release fraction (atmosphere) is',BubbRelFractionAtm,', which is <= 0.0'
    end if

    IF (BubbWatGasExchRate <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTF) 'Sediment Diagenesis:Bubbles-Water gas exchange rate (1/s) is',BubbWatGasExchRate,', which is <= 0.0'
    end if

    IF (CEMATurbulenceScaling <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTF) 'Sediment Diagenesis:Turbulence scaling factor for bubbles release (1/s) is',CEMATurbulenceScaling,', which is <= 0.0'
    end if

        
        
    ELSE
        DO JSKIP=1,16
            READ(CEMAFilN,*,err=6500)
        ENDDO
    END IF   
       
        Read(CEMAFilN,*,err=6500)MessageTemp, CEMA_POM_Resuspension
        
        IF(CEMA_POM_Resuspension) THEN
          Read(CEMAFilN,*,err=6500)MessageTemp, TAUCRPOM
          Read(CEMAFilN,*,err=6500)MessageTemp, crshields
          Read(CEMAFilN,*,err=6500)MessageTemp, cao_method
          Read(CEMAFilN,*,err=6500)MessageTemp, spgrav_POM
          Read(CEMAFilN,*,err=6500)MessageTemp, dia_POM
          
              IF (TAUCRPOM <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTF) 'Sediment Diagenesis:Critical shear stress for particulate organic matter resuspension dynes/cm^2 is',TAUCRPOM,', which is <= 0.0'
    end if
    IF (crshields <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTF) 'Sediment Diagenesis:Critical Shields parameter for POM is',crshields,', which is <= 0.0'
    end if
    IF (spgrav_POM <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTF) 'Sediment Diagenesis:Specific gravity of particulate organic matter is',spgrav_POM,', which is <= 0.0'
    end if
    IF (dia_POM <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTF) 'Sediment Diagenesis:Particulate organic matter particle diameter m is',dia_POM,', which is <= 0.0'
    end if

        ELSE
            DO JSKIP=1,5
            READ(CEMAFilN,*,err=6500)
            ENDDO
        END IF        
        
        Read(CEMAFilN,*,err=6500)MessageTemp, IncludeAlkalinity
        Read(CEMAFilN,*,err=6500)MessageTemp, IncludeIron        
        Read(CEMAFilN,*,err=6500)MessageTemp, IncludeManganese
        !
        IF(IncludeAlkalinity) IncludeDynamicpH = .TRUE.
        Read(CEMAFilN,*,err=6500)MessageTemp, NumRegnsSedimentBedComposition
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
        Read(CEMAFilN,*,err=6500)   ! skip line for header
        Read(CEMAFilN,*,err=6500)MessageTemp, (SedBedInitRegSegSt(i), i = 1, NumRegnsSedimentBedComposition)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SedBedInitRegSegEn(i), i = 1, NumRegnsSedimentBedComposition)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnT_T(i),   i = 1, NumRegnsSedimentBedComposition)
        IF(.NOT. IncludeDynamicpH) THEN
            Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnpH(i), i = 1, NumRegnsSedimentBedComposition) 
        ELSE
            Read(CEMAFilN,*,err=6500)
        ENDIF
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnPOC_T(i), i = 1, NumRegnsSedimentBedComposition)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnPON_T(i), i = 1, NumRegnsSedimentBedComposition)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnPOP_T(i), i = 1, NumRegnsSedimentBedComposition)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnSul_T(i), i = 1, NumRegnsSedimentBedComposition)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnNH3_T(i), i = 1, NumRegnsSedimentBedComposition)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnNO3_T(i), i = 1, NumRegnsSedimentBedComposition)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnPO4_T(i), i = 1, NumRegnsSedimentBedComposition)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnH2S_T(i), i = 1, NumRegnsSedimentBedComposition)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnCH4_T(i), i = 1, NumRegnsSedimentBedComposition)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnTIC_T(i), i = 1, NumRegnsSedimentBedComposition)
        
            do i=1,NumRegnsSedimentBedComposition
      if(SedBedInitRegSegSt(i) < 1 .or. SedBedInitRegSegSt(i) > IMX)then
          CALL ERRORS
          WRITE(ERR,FMTI)'Sediment Diagenesis:Starting segment for sediment diagenesis region outside model grid: I=',SedBedInitRegSegSt(i),' is less than 1 or greater than IMX=',imx
      endif
 
      if(SedBedInitRegSegEn(i) < 1 .or. SedBedInitRegSegEn(i) > IMX)then
          CALL ERRORS
          WRITE(ERR,FMTI)'Sediment Diagenesis:Ending segment for sediment diagenesis region outside model grid: I=',SedBedInitRegSegEn(i),' is less than 1 or greater than IMX=',imx
      endif
      if(SedBedInitRegSegSt(i) > SedBedInitRegSegEn(i))then
          CALL ERRORS
          WRITE(ERR,FMTI)'Sediment Diagenesis:Starting segment in sediment diagenesis region  I=',SedBedInitRegSegSt(i),' is greater than ending segment I=',SedBedInitRegSegEn(i)
      endif

      IF (SDRegnT_T(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Initial temperature of sediment region ',i, ' is ',SDRegnT_T(i),' and < 0.0 '
      END IF

      IF (SDRegnPOC_T(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Initial particulate organic carbon (total) concentration of sediment region ',i, ' is ',SDRegnPOC_T(i),' and < 0.0 '
      ELSE IF(SDRegnPOC_T(i) == 0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Initial particulate organic carbon (total) concentration of sediment region ',i, ' is 0'
      end if

      IF (SDRegnPON_T(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Initial particulate organic nitrogen (total) concentration of sediment region ',i, ' is ',SDRegnPON_T(i),' and < 0.0 '
      ELSE IF(SDRegnPON_T(i) == 0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Initial particulate organic nitrogen (total) concentration of sediment region ',i, ' is 0'
      end if

      IF (SDRegnPOP_T(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Initial particulate organic phosphorus (total) concentration of sediment region ',i, ' is ',SDRegnPOP_T(i),' and < 0.0 '
      ELSE IF(SDRegnPOP_T(i) == 0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Initial particulate organic phosphorus (total) concentration of sediment region ',i, ' is 0'
      end if

      IF (SDRegnSul_T(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Initial sulfate concentration of sediment region ',i, ' is ',SDRegnSul_T(i),' and < 0.0 '
      ELSE IF(SDRegnSul_T(i) == 0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Initial sulfate concentration of sediment region ',i, ' is 0'
      end if

      IF (SDRegnNH3_T(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Initial dissolved ammonia concentration of sediment region ',i, ' is ',SDRegnNH3_T(i),' and < 0.0 '
      ELSE IF(SDRegnNH3_T(i) == 0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Initial dissolved ammonia concentration of sediment region ',i, ' is 0'
      end if

      IF (SDRegnPO4_T(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Initial phosphate concentration of sediment region ',i, ' is ',SDRegnPO4_T(i),' and < 0.0 '
      ELSE IF(SDRegnPO4_T(i) == 0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Initial phosphate concentration of sediment region ',i, ' is 0'
      end if

      IF (SDRegnH2S_T(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Initial sulfide concentration of sediment region ',i, ' is ',SDRegnH2S_T(i),' and < 0.0 '
      ELSE IF(SDRegnH2S_T(i) == 0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Initial sulfide concentration of sediment region ',i, ' is 0'
      end if

      IF (SDRegnCH4_T(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Initial methane concentration of sediment region ',i, ' is ',SDRegnCH4_T(i),' and < 0.0 '
      ELSE IF(SDRegnCH4_T(i) == 0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Initial methane concentration of sediment region ',i, ' is 0'
      end if

      IF (SDRegnTIC_T(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Initial total inorganic carbon concentration of sediment region ',i, ' is ',SDRegnTIC_T(i),' and < 0.0 '
      ELSE IF(SDRegnTIC_T(i) == 0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Initial total inorganic carbon concentration of sediment region ',i, ' is 0'
      end if
    end do

        
        
        IF(IncludeAlkalinity)THEN
            Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnALK_T(i), i = 1, NumRegnsSedimentBedComposition)
              do i=1,NumRegnsSedimentBedComposition
              IF (SDRegnALK_T(i) < 0.0) THEN
                CALL ERRORS
                WRITE (ERR,FMTIF) 'Sediment Diagenesis:Initial alkalinity concentration of sediment region ',i, ' is ',SDRegnALK_T(i),' and < 0.0 '
              ELSE IF(SDRegnALK_T(i) == 0.0)then
                CALL warnings
                WRITE (wrn,FMTIF) 'Sediment Diagenesis:Initial alkalinity concentration of sediment region ',i, ' is 0'
              end if
             end do
        ELSE
            READ(CEMAFilN,*,err=6500)
        ENDIF
        !
        IF(IncludeIron) THEN
          Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnFe2_T(i),   i = 1, NumRegnsSedimentBedComposition)
          Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnFeOOH_T(i), i = 1, NumRegnsSedimentBedComposition)
    do i=1,NumRegnsSedimentBedComposition
      IF (SDRegnFe2_T(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Initial ferrous iron concentration of sediment region ',i, ' is ',SDRegnFe2_T(i),' and < 0.0 '
      ELSE IF(SDRegnFe2_T(i) == 0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Initial ferrous iron concentration of sediment region ',i, ' is 0'
      end if
      IF (SDRegnFeOOH_T(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Initial iron oxyhyroxide concentration of sediment region ',i, ' is ',SDRegnFeOOH_T(i),' and < 0.0 '
      ELSE IF(SDRegnFeOOH_T(i) == 0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Initial iron oxyhyroxide concentration of sediment region ',i, ' is 0'
      end if
    end do
          
        ELSE
            DO JSKIP=1,2
            READ(CEMAFilN,*,err=6500)
            ENDDO            
        END IF
        IF(IncludeManganese) THEN
          Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnMn2_T(i),  i = 1, NumRegnsSedimentBedComposition)
          Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnMnO2_T(i), i = 1, NumRegnsSedimentBedComposition)
    do i=1,NumRegnsSedimentBedComposition
      IF (SDRegnMn2_T(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Initial Mn(II) concentration of sediment region ',i, ' is ',SDRegnMn2_T(i),' and < 0.0 '
      ELSE IF(SDRegnMn2_T(i) == 0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Initial Mn(II) concentration of sediment region ',i, ' is 0'
      end if
    end do
    do i=1,NumRegnsSedimentBedComposition
      IF (SDRegnMnO2_T(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Initial manganese dioxide concentration of sediment region ',i, ' is ',SDRegnMnO2_T(i),' and < 0.0 '
      ELSE IF(SDRegnMnO2_T(i) == 0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Initial manganese dioxide concentration of sediment region ',i, ' is 0'
      end if
    end do          
        ELSE
            DO JSKIP=1,2
            READ(CEMAFilN,*,err=6500)
            ENDDO            
        END IF
        !
        Read(CEMAFilN,*,err=6500)MessageTemp, NumRegnsSedimentDiagenesis
        Read(CEMAFilN,*,err=6500)    ! SKIP LINE FOR HEADER
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
        Read(CEMAFilN,*,err=6500)MessageTemp, (SedBedDiaRCRegSegSt(i),     i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SedBedDiaRCRegSegEn(i),     i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnPOC_L_Fr(i),          i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnPOC_R_Fr(i),          i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnPON_L_Fr(i),          i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnPON_R_Fr(i),          i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnPOP_L_Fr(i),          i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnPOP_R_Fr(i),          i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnPW_DiffCoeff(i),      i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (PartMixVel(i),              i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (BurialVel(i),               i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (POCr(i),                    i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegn_CH4CompMethod(i),    i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnOx_Threshold(i),      i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnAe_NH3_NO3_L(i),      i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnAe_NH3_NO3_H(i),      i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnAe_NO3_N2_L(i),       i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnAe_NO3_N2_H(i),       i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnAn_NO3_N2(i),         i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnAe_CH4_CO2(i),        i = 1, NumRegnsSedimentDiagenesis)   !Eq. 10.35
        Read(CEMAFilN,*,err=6500)MessageTemp, (KsOxch(i),                  i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnAe_HS_NH4_Nit(i),     i = 1, NumRegnsSedimentDiagenesis)   !Eq. 3.3
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnAe_HS_O2_Nit(i),      i = 1, NumRegnsSedimentDiagenesis)   !Eq. 3.3
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegn_Theta_PW(i),         i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegn_Theta_PM(i),         i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegn_Theta_NH3_NO3(i),    i = 1, NumRegnsSedimentDiagenesis)   
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegn_Theta_NO3_N2(i),     i = 1, NumRegnsSedimentDiagenesis)  
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegn_Theta_CH4_CO2(i),    i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegn_Sulfate_CH4_H2S(i),  i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegnAe_H2S_SO4(i),        i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegn_Theta_H2S_SO4(i),    i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegn_NormConst_H2S_SO4(i),i = 1, NumRegnsSedimentDiagenesis)   !Eq. 9.6
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegn_MinRate_POC_Lab(i),  i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegn_MinRate_POC_Ref(i),  i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegn_MinRate_POC_Ine(i),  i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegn_MinRate_PON_Lab(i),  i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegn_MinRate_PON_Ref(i),  i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegn_MinRate_PON_Ine(i),  i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegn_MinRate_POP_Lab(i),  i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegn_MinRate_POP_Ref(i),  i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegn_MinRate_POP_Ine(i),  i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegn_Theta_POC_Lab(i),    i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegn_Theta_POC_Ref(i),    i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegn_Theta_POC_Ine(i),    i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegn_Theta_PON_Lab(i),    i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegn_Theta_PON_Ref(i),    i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegn_Theta_PON_Ine(i),    i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegn_Theta_POP_Lab(i),    i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegn_Theta_POP_Ref(i),    i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegn_Theta_POP_Ine(i),    i = 1, NumRegnsSedimentDiagenesis)    
        Read(CEMAFilN,*,err=6500)MessageTemp, (Kdp2(i),                    i = 1, NumRegnsSedimentDiagenesis) 
        Read(CEMAFilN,*,err=6500)MessageTemp, (delta_kpo41(i),             i = 1, NumRegnsSedimentDiagenesis) 
        Read(CEMAFilN,*,err=6500)MessageTemp, (DOcr(i),                    i = 1, NumRegnsSedimentDiagenesis) 
        Read(CEMAFilN,*,err=6500)MessageTemp, (KdNH31(i),                  i = 1, NumRegnsSedimentDiagenesis)    
        Read(CEMAFilN,*,err=6500)MessageTemp, (KdNH32(i),                  i = 1, NumRegnsSedimentDiagenesis) 
        Read(CEMAFilN,*,err=6500)MessageTemp, (KdH2S1(i),                  i = 1, NumRegnsSedimentDiagenesis) 
        Read(CEMAFilN,*,err=6500)MessageTemp, (KdH2S2(i),                  i = 1, NumRegnsSedimentDiagenesis) 
        Read(CEMAFilN,*,err=6500)MessageTemp, (SDRegn_POMResuspMethod(i),  i = 1, NumRegnsSedimentDiagenesis)

        Read(CEMAFilN,*,err=6500)MessageTemp, (KdFe1(i),               i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (KdFe2(i),               i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (KdMn1(i),               i = 1, NumRegnsSedimentDiagenesis)
        Read(CEMAFilN,*,err=6500)MessageTemp, (KdMn2(i),               i = 1, NumRegnsSedimentDiagenesis)
        
    IF (IncludeCEMASedDiagenesis .and. NumRegnsSedimentBedComposition <= 0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Sediment Diagenesis:Sediment diagenesis processes are turned on but the number of regions for different initial sediment concentrations ',        NumRegnsSedimentBedComposition,' <= 0 '
    END IF
    
    IF (IncludeCEMASedDiagenesis .and. NumRegnsSedimentDiagenesis <= 0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Sediment Diagenesis:Sediment diagenesis processes are turned on but the number of regions for for different diagenesis related rates ',        NumRegnsSedimentBedComposition,' <= 0 '
    END IF
    do i=1,NumRegnsSedimentDiagenesis
      if(SedBedDiaRCRegSegSt(i) < 1 .or. SedBedDiaRCRegSegSt(i) > IMX)then
          CALL ERRORS
          WRITE(ERR,FMTI)'Sediment Diagenesis:Starting segment for region of different diagenesis related rates outside model grid: I=',SedBedDiaRCRegSegSt(i),' is less than 1 or greater than IMX=',imx
      endif

      if(SedBedDiaRCRegSegEn(i) < 1 .or. SedBedDiaRCRegSegEn(i)> IMX)then
          CALL ERRORS
          WRITE(ERR,FMTI)'Sediment Diagenesis:Ending segment for region of different diagenesis related rates outside model grid: I=',SedBedDiaRCRegSegEn(i),' is less than 1 or greater than IMX=',imx
      endif
      if(SedBedDiaRCRegSegSt(i) > SedBedDiaRCRegSegEn(i))then
          CALL ERRORS
          WRITE(ERR,FMTI)'Sediment Diagenesis:Starting segment in region of different diagenesis related rates I=',SedBedDiaRCRegSegSt(i),' is greater than ending segment I=',SedBedDiaRCRegSegEn(i)
      endif
      IF (SDRegnPOC_L_Fr(i) < 0.0) THEN
          SDRegnPOC_L_Fr(1)=ABS(SDRegnPOC_L_Fr(1))
          DYNAMIC_SD=.TRUE.   ! DYNAMIC NUTRIENT ACCUMULATION BASED ON STOICHIOMETRY
      ENDIF
      
      IF (SDRegnPOC_L_Fr(i) >1.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Fraction of labile poc for diagenesis region ',i, ' is ',SDRegnPOC_L_Fr(i),' ,which is > 1.0'
      ELSE IF(SDRegnPOC_L_Fr(i) == 0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Fraction of labile poc for diagenesis region ',i, ' is 0'
      end if

      IF (SDRegnPOC_R_Fr(i) < 0.0 .or. SDRegnPOC_R_Fr(i) >1.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Fraction of refractory poc for diagenesis region ',i, ' is ',SDRegnPOC_R_Fr(i),', which is < 0.0 or > 1.0'
      ELSE IF(SDRegnPOC_R_Fr(i) == 0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Fraction of refractory poc for diagenesis region ',i, ' is 0'
      end if

      IF (SDRegnPON_L_Fr(i) < 0.0 .or. SDRegnPON_L_Fr(i) >1.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Fraction of labile pon for diagenesis region ',i, ' is ',SDRegnPON_L_Fr(i),', which is < 0.0 or > 1.0'
      ELSE IF(SDRegnPON_L_Fr(i) == 0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Fraction of labile pon for diagenesis region ',i, ' is 0'
      end if

      IF (SDRegnPON_R_Fr(i) < 0.0 .or. SDRegnPON_R_Fr(i) >1.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Fraction of refractory pon for diagenesis region ',i, ' is ',SDRegnPON_R_Fr(i),', which is < 0.0 or > 1.0'
      ELSE IF(SDRegnPON_R_Fr(i)== 0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Fraction of refractory pon for diagenesis region ',i, ' is 0'
      end if

      IF (SDRegnPOP_L_Fr(i) < 0.0 .or. SDRegnPOP_L_Fr(i) >1.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Fraction of labile pop for diagenesis region ',i, ' is ',SDRegnPOP_L_Fr(i),', which is < 0.0 or > 1.0'
      ELSE IF(SDRegnPOP_L_Fr(i) == 0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Fraction of labile pop for diagenesis region ',i, ' is 0'
      end if
 
      IF (SDRegnPOP_R_Fr(i) < 0.0 .or. SDRegnPOP_R_Fr(i) >1.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Fraction of refractory pop for diagenesis region ',i, ' is ',SDRegnPOP_R_Fr(i),', which is < 0.0 or > 1.0'
      ELSE IF(SDRegnPOP_R_Fr(i)== 0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Fraction of refractory pop for diagenesis region ',i, ' is 0'
      end if

      IF (SDRegnPW_DiffCoeff(i) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Pore water diffusion coefficient for diagenesis region ',i, ' is ',SDRegnPW_DiffCoeff(i),', which is <= 0.0'     
      end if

      IF (SDRegnOx_Threshold(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:DO Threshold for aerobic layer oxidation rates mgO2/l of region ',i, ' is ',SDRegnOx_Threshold(i),', which is < 0.00'
      ELSE IF(SDRegnOx_Threshold(i) > 5.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:DO Threshold for aerobic layer oxidation rates mgO2/l of region ',i, ' is greater than 5 mg/l'
      end if

      IF (SDRegnAe_NH3_NO3_L(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Nitrification rate in aerobic layer (NH3->NO3) at DO below threshold for region ',i, ' is ',SDRegnAe_NH3_NO3_L(i),', which is < 0.00'
      ELSE IF(SDRegnAe_NH3_NO3_L(i) > 0.5)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Nitrification rate in aerobic layer (NH3->NO3) at DO below threshold of region ',i, ' is greater than 0.5'
      end if

      IF (SDRegnAe_NH3_NO3_H(i) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Nitrification rate in aerobic layer (NH3->NO3) at DO above threshold for region ',i, ' is ',SDRegnAe_NH3_NO3_H(i),', which is <= 0.00'
      ELSE IF(SDRegnAe_NH3_NO3_H(i) > 0.5)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Nitrification rate in aerobic layer (NH3->NO3) at DO above threshold of region ',i, ' is greater than 0.5'
      end if

      IF (SDRegnAe_NO3_N2_L(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Denitrification rate in aerobic layer (NO3->N2) at DO below threshold for region ',i, ' is ',SDRegnAe_NO3_N2_L(i),', which is < 0.00'
      ELSE IF(SDRegnAe_NO3_N2_L(i) > 0.5)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Denitrification rate in aerobic layer (NO3->N2) at DO below threshold for region ',i, ' is greater than 0.5'
      end if
  
      IF (SDRegnAe_NO3_N2_H(i) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Denitrification rate in aerobic layer (NO3->N2) at DO above threshold for region ',i, ' is ',SDRegnAe_NO3_N2_H(i),', which is <= 0.00'
      ELSE IF(SDRegnAe_NO3_N2_H(i) > 0.5)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Denitrification rate in aerobic layer (NO3->N2) at DO above threshold for region ',i, ' is greater than 0.5'
      end if
 
      IF (SDRegnAn_NO3_N2(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Denitrification rate in anerobic layer (NO3->N2) for region ',i, ' is ',SDRegnAn_NO3_N2(i),', which is < 0.00'
      ELSE IF(SDRegnAn_NO3_N2(i) > 0.5)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Denitrification rate in anerobic layer (NO3->N2) for region ',i, ' is greater than 0.5'
      end if
   
      IF (SDRegnAe_CH4_CO2(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Methane oxidation rate in aerobic layer for region ',i, ' is ',SDRegnAe_CH4_CO2(i),', which is < 0.0'
      ELSE IF(SDRegnAe_CH4_CO2(i) > 3.0 .or. SDRegnAe_CH4_CO2(i) == 0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Methane oxidation rate in aerobic layer for region ',i, ' is greater than 3.0 or equal to zero'
      end if

      IF (SDRegnAe_HS_NH4_Nit(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Nitrification half-saturation constant for NH4N in aerobic layer mgN/l for region ',i, ' is ',SDRegnAe_HS_NH4_Nit(i),', which is < 0.0'
      ELSE IF(SDRegnAe_HS_NH4_Nit(i) > 3.0 .or. SDRegnAe_HS_NH4_Nit(i) == 0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Nitrification half-saturation constant for NH4N in aerobic layer mgN/l for region ',i, ' is greater than 3.0'
      end if
 
      IF (SDRegnAe_HS_O2_Nit(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Nitrification half-saturation constant for O2 in aerobic layer mgO2/l for region ',i, ' is ',SDRegnAe_HS_O2_Nit(i),', which is < 0.0'
      ELSE IF(SDRegnAe_HS_O2_Nit(i) > 3.0 .or. SDRegnAe_HS_O2_Nit(i) == 0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Nitrification half-saturation constant for O2 in aerobic layer mgO2/l for region ',i, ' is greater than 3.0 or equal to zero'
      end if

      IF (SDRegn_Sulfate_CH4_H2S(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Sulfate concentration above which sulfide over methane is produced mgS/l for region ',i, ' is ',SDRegn_Sulfate_CH4_H2S(i),', which is < 0.0'
      End if

      IF (SDRegnAe_H2S_SO4(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Sulfide oxidation rate in aerobic layer m/d for region ',i, ' is ',SDRegnAe_H2S_SO4(i),', which is < 0.0'
      ELSE IF(SDRegnAe_H2S_SO4(i) > 3.0  .or. SDRegnAe_H2S_SO4(i) ==0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Sulfide oxidation rate in aerobic layer m/d for region ',i, ' is greater than 3.0 or equal to zero'
      end if

      IF (SDRegn_MinRate_PON_Lab(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Mineralization rate for labile PON 1/d for region ',i, ' is ',SDRegn_MinRate_PON_Lab(i),', which is < 0.0'
      ELSE IF(SDRegn_MinRate_PON_Lab(i) > 0.3 .or. SDRegn_MinRate_PON_Lab(i)==0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Mineralization rate for labile PON 1/d for region ',i, ' is greater than 0.3 or equal to zero'
      end if
 
      IF (SDRegn_MinRate_PON_Ref(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Mineralization rate for refractory PON 1/d for region ',i, ' is ',SDRegn_MinRate_PON_Ref(i),', which is < 0.0'
      ELSE IF(SDRegn_MinRate_PON_Ref(i) > 0.01 .or. SDRegn_MinRate_PON_Ref(i)==0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Mineralization rate for refractory PON 1/d for region ',i, ' is greater than 0.01 or equal to zero'
      end if
  
      IF (SDRegn_MinRate_PON_Ine(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Mineralization rate for inert/slow refractory PON 1/d for region ',i, ' is ',SDRegn_MinRate_PON_Ine(i),', which is < 0.0'
      ELSE IF(SDRegn_MinRate_PON_Ine(i) > 0.001 .or. SDRegn_MinRate_PON_Ine(i)==0.0 )then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Mineralization rate for inert/slow refractory PON 1/d for region ',i, ' is greater than 0.001 or equal to zero'
      end if
   
      IF (SDRegn_MinRate_POC_Lab(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Mineralization rate for labile POC 1/d for region ',i, ' is ',SDRegn_MinRate_POC_Lab(i),', which is < 0.0'
      ELSE IF(SDRegn_MinRate_POC_Lab(i) > 0.3 .or. SDRegn_MinRate_POC_Lab(i) == 0.0 )then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Mineralization rate for labile POC 1/d for region ',i, ' is greater than 0.3 or equal to zero'
      end if

      IF (SDRegn_MinRate_POC_Ref(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Mineralization rate for refractory POC 1/d for region ',i, ' is ',SDRegn_MinRate_POC_Ref(i),', which is < 0.0'
      ELSE IF(SDRegn_MinRate_POC_Ref(i) > 0.01 .or. SDRegn_MinRate_POC_Ref(i) == 0.00)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Mineralization rate for refractory POC 1/d for region ',i, ' is greater than 0.01 or equal to zero'
      end if

      IF (SDRegn_MinRate_POC_Ine(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Mineralization rate for inert/slow refractory PON 1/d for region ',i, ' is ',SDRegn_MinRate_POC_Ine(i),', which is < 0.0'
      ELSE IF(SDRegn_MinRate_POC_Ine(i) > 0.001 .or. SDRegn_MinRate_POC_Ine(i) == 0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Mineralization rate for inert/slow refractory PON 1/d for region ',i, ' is greater than 0.001 or equal to zero'
      end if

      IF (SDRegn_MinRate_POP_Lab(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Mineralization rate for labile POP 1/d for region ',i, ' is ',SDRegn_MinRate_POP_Lab(i),', which is < 0.0'
      ELSE IF(SDRegn_MinRate_POP_Lab(i) > 0.3 .or. SDRegn_MinRate_POP_Lab(i) == 0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Mineralization rate for labile POP 1/d for region ',i, ' is greater than 0.3 or equal to zero'
      end if

      IF (SDRegn_MinRate_POP_Ref(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Mineralization rate for refractory POP 1/d for region ',i, ' is ',SDRegn_MinRate_POP_Ref(i),', which is < 0.0'
      ELSE IF(SDRegn_MinRate_POP_Ref(i) > 0.01 .or. SDRegn_MinRate_POP_Ref(i) == 0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Mineralization rate for refractory POP 1/d for region ',i, ' is greater than 0.01 or equal to zero'
      end if

      IF (SDRegn_MinRate_POP_Ine(i) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Mineralization rate for inert/slow refractory POP 1/d for region ',i, ' is ',SDRegn_MinRate_POP_Ine(i),', which is < 0.0'
      ELSE IF(SDRegn_MinRate_POP_Ine(i) > 0.001 .or. SDRegn_MinRate_POP_Ine(i) == 0.0)then
        CALL warnings
        WRITE (wrn,FMTIF) 'Sediment Diagenesis:Mineralization rate for inert/slow refractory POP 1/d for region ',i, ' is greater than 0.001 or equal to zero'
      end if

      IF (SDRegn_Theta_PON_Lab(i) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Temperature coefficient for labile PON for region ',i, ' is ',SDRegn_Theta_PON_Lab(i),', which is <= 0.0'
      end if

      IF (SDRegn_Theta_PON_Ref(i) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Temperature coefficient refractory PON for region ',i, ' is ',SDRegn_Theta_PON_Ref(i),', which is <= 0.0'
      end if
 
      IF (SDRegn_Theta_PON_Ine(i) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Temperature coefficient inert PON for region ',i, ' is ',SDRegn_Theta_PON_Ine(i),', which is <= 0.0'
      end if
 
      IF (SDRegn_Theta_POC_Lab(i) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Temperature coefficient labile POC for region ',i, ' is ',SDRegn_Theta_POC_Lab(i),', which is <= 0.0'
      end if
 
      IF (SDRegn_Theta_POC_Ref(i) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Temperature coefficient refractory POC for region ',i, ' is ',SDRegn_Theta_POC_Ref(i),', which is <= 0.0'
      end if

      IF (SDRegn_Theta_POC_Ine(i) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Temperature coefficient inert POC for region ',i, ' is ',SDRegn_Theta_POC_Ine(i),', which is <= 0.0'
      end if
 
      IF (SDRegn_Theta_POP_Lab(i) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Temperature coefficient labile POP for region ',i, ' is ',SDRegn_Theta_POP_Lab(i),', which is <= 0.0'
      end if
   
      IF (SDRegn_Theta_POP_Ref(i) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Temperature coefficient refractory POP for region ',i, ' is ',SDRegn_Theta_POP_Ref(i),', which is <= 0.0'
      end if

      IF (SDRegn_Theta_POP_Ine(i) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Sediment Diagenesis:Temperature coefficient inert POP for region ',i, ' is ',SDRegn_Theta_POP_Ine(i),', which is <= 0.0'
      end if

      IF (Kdp2(i) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,'(a,f12.6,a,i3)') 'Sediment Diagenesis:Phosphorus Sorption Coefficient in Anaerobic Layer m^3/g is ',Kdp2(i),', which is <= 0.0 for region:', i
      end if   
 
      if(SDRegn_CH4CompMethod(i) /= 0 .and. SDRegn_CH4CompMethod(i) /= 1)then
          CALL ERRORS
          WRITE(ERR,FMTI)'Sediment Diagenesis:Methane production calculation method must be O or 1 for region ',i
      endif

    IF (PartMixVel(i) <= 0.0 .or. PartMixVel(i) > 5.0) THEN
      CALL ERRORS
      WRITE (ERR,'(a,f12.6,a,i3)') 'Sediment Diagenesis:Particle Mixing Velocity between Aerobic and Anaerobic Layers m/d is',PartMixVel(i),', which is <= 0.0 or > 5 m/d for region:', i
    end if
     IF (BurialVel(i) <= 0.0 .or. BurialVel(i) > 5.0) THEN
      CALL ERRORS
      WRITE (ERR,'(a,f12.6,a,i3)') 'Sediment Diagenesis:Sediment Diagenesis burial Velocity m/d is',BurialVel(i),', which is <= 0.0 or > 5.0 m/d for region:',i
     end if
 
      if(SDRegn_POMResuspMethod(i) /= 0 .and. SDRegn_POMResuspMethod(i) /= 1)then
          CALL ERRORS
          WRITE(ERR,FMTI)'Sediment Diagenesis:Algorithm for resuspension of particulate organic matter must be O or 1 for region ',i
      endif
    end do
    ! GROUP 6: Output
        Read(CEMAFilN,*,err=6500)MessageTemp, WriteCEMAMFTSedFlx
        READ(CEMAFilN,*,err=6500)MessageTemp, SEDIAGFREQ   ! FREQUENCY OF OUTPUT SW 5/25/2017
    IF (SEDIAGFREQ <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTF) 'Sediment Diagenesis:SEDIAGFREQ Frequency of output in days is',SEDIAGFREQ,', which is <= 0.0'
    end if

    END IF
    GO TO 6501
6500   CALL ERRORS
       WRITE(ERR,'(A,A)')'w2_diagenesis.npt input error. The last line Message was:',MessageTemp
       backspace(CEMAFilN)
       backspace(CEMAFilN)
       Read(CEMAFilN,'(a80)')title80
       write(err,*)'Last 2 lines of w2_diagenesis.npt were:'
       write(err,'(A80)')title80
       Read(CEMAFilN,'(a80)')title80
       write(err,'(A80)')title80
6501    close(CEMAFilN) 
 
  End If
ELSE
  INQUIRE(FILE="W2_diagenesis.npt", EXIST=DIAG_EXISTS)   ! file_exists will be TRUE if the file
  IF(DIAG_EXISTS) THEN	
      CALL WARNINGS
      WRITE(WRN,'(A)')'File w2_diagenesis.npt exists but SED_DIAG in control file is OFF, hence there will be no sediment diagenesis processes modeled.'
  ENDIF
ENDIF

 
! Grid geometry

  WRITE (WIN,*) 'Bathymetry file checks'
  WRITE (WIN,*) '  water surface elevations'
  DO JW=1,NWB
    ZMIN = -1000.0
    DO JB=BS(JW),BE(JW)
      ALPHA(JB) = ATAN(SLOPE(JB))
      IFLAG     = 0
      DO I=US(JB),DS(JB)
        KTI(I) = 2
        DO WHILE (EL(KTI(I),I) > ELWS(I) .AND. KTI(I) < KMX)
          KTI(I) = KTI(I)+1
          IF (KTI(I) == KMX .AND. EL(KTI(I),I)> ELWS(I)) THEN     !         IF (KTI(I) > KMX-1) THEN    SW 5/27/2017
            CALL ERRORS
            WRITE (ERR,FMTI) 'Initial water surface elevation below grid at segment ',I
          END IF
          ENDDO
        Z(I)     = (EL(KTI(I),I)-ELWS(I))/COS(ALPHA(JB))
        ZMIN     =  MAX(ZMIN,Z(I))
        ! KTI(I)   =  MAX(KTI(I)-1,2)   ! SW 5/27/17 MOVED FROM BELOW
        KTMAX    =  MAX(2,KTI(I))
        KTWB(JW) =  MAX(KTMAX,KTWB(JW))
        KTI(I)   =  MAX(KTI(I)-1,2)
        IF (KTWB(JW) >= KMX .AND. IFLAG == 0) THEN
          IFLAG =  1
          CALL ERRORS
          WRITE (ERR,FMTI) 'Water surface elevation too low in branch ',JB,' at segment ',I
        END IF
        DO K=2,KMX-1
          IF (B(K+1,I) == 0.0 .OR. K == KMX-1) THEN
            KB(I) = K
            EXIT
          END IF
        END DO
        KBMAX(JW) = MAX(KBMAX(JW),KB(I))
      END DO
      KB(US(JB)-1) = KB(US(JB))
      KB(DS(JB)+1) = KB(DS(JB))
    END DO
  END DO

  ! Constrictions                       ! SW 8/5/2018

   WRITE (WIN,*) 'Constriction File Check'
    INQUIRE(FILE='constriction.csv',EXIST=Constriction)    
    IF(CONSTRICTION)THEN
        CALL ERRORS
        WRITE (ERR,'(A,i4,A)') 'CONSTRICTION: File name must be changed to w2_constriction.csv(new file name for V4.3 and later) from constriction.csv(old file name)'
        CONSTRICTION=.FALSE.
    ENDIF
    INQUIRE(FILE='w2_constriction.csv',EXIST=Constriction)    

    IF(CONSTRICTION)THEN
        WRITE (WIN,*) '  Constriction ON: Reading constriction.csv'
        OPEN(NPT,FILE='w2_constriction.csv',STATUS='OLD')
        READ(NPT,*)
        READ(NPT,*)INFORMAT2,ICON    ! NUMBER OF CONSTRICTIONS
        IF(INFORMAT2 /= 'ON' .AND. INFORMAT2 /= 'OF')THEN
        CALL ERRORS
        WRITE (ERR,'(A,A,A)') 'CONSTRICTION: Reading file w2_constriction.csv: ON/OFF switch[',INFORMAT2,'] is not equal to ON of OFF'
        ENDIF
        
        
        IF (ICON >= IMX) THEN
        CALL ERRORS
        WRITE (ERR,'(A,i4,A)') 'CONSTRICTION: Reading file w2_constriction.csv: # of constrictions[',ICON,'] >= IMX '
        ENDIF
        IF (ICON <= 0) THEN
        CALL ERRORS
        WRITE (ERR,'(A,i4,A)') 'CONSTRICTION: Reading file w2_constriction.csv: # of constrictions[',ICON,'] <= 0 '
        ENDIF

        READ(NPT,*)
        DO J=1,ICON
            READ(NPT,*)ISEG,BCONSTRICTION
                    IF (BCONSTRICTION <= 0) THEN
                    CALL ERRORS
                    WRITE (ERR,'(A,F9.2,A,I4)') 'CONSTRICTION: Reading file w2_constriction.csv: CONSTRICTION WIDTH in m[',BCONSTRICTION,'] <= 0 for segment #:',ISEG
                    ENDIF
                    ACTIVE_SEGMENT = .FALSE.
                    DO JB=1,NBR
                      IF (ISEG >= US(JB) .AND. ISEG <= DS(JB)) ACTIVE_SEGMENT = .TRUE.
                    END DO
                    IF (.NOT.ACTIVE_SEGMENT) THEN
                      CALL ERRORS
                      WRITE (ERR,'(A,i4,A)') 'CONSTRICTION: Reading file w2_constriction.csv: Constriction segment [ISEG=',ISEG,'] is a boundary segment or not in the active grid'
                    END IF
        ENDDO
    ENDIF
    

    CLOSE(NPT)
    
! SysTDG input file 4.2
   WRITE (WIN,*) 'SYSTDG input file checks'
  SYSTDG=.FALSE.
  INQUIRE(FILE='w2_systdg.npt',EXIST=SYSTDG)   
  IF(SYSTDG)THEN
      POWNO = 0
      FLNO  = 0
      NBAY  = 0
      NRO   = 0
   WRITE (WIN,*) '  SYSTDG ON: Reading w2_systdg.npt'

    OPEN (NPT, FILE='w2_systdg.npt', STATUS='OLD')
    ! Check for csv file or old format
     CSVFORMAT=.FALSE.   
     READ (NPT,'(A)')HEADER(1)
     IF(HEADER(1)(1:1)=='$')CSVFORMAT=.TRUE.
     
     IF(.NOT.CSVFORMAT)THEN
     READ (NPT,'(A)')HEADER(1)    ! READ NEXT LINE - IF COMMAS IN FIRST FEW FIELDS IT IS IN CSV FORMAT
     DO J=1,7
         IF(HEADER(1)(J:J)==',')THEN
             CSVFORMAT=.TRUE.
             EXIT
         ENDIF
     ENDDO
     REWIND(NPT)
     ENDIF
      
      IF(CSVFORMAT)THEN
      AID = 'STDG:LN1'
      READ (NPT,*,ERR=400)
      READ (NPT,*,ERR=400)
      READ (NPT,*,ERR=400)
      DO I=1,10
      READ (NPT,*,ERR=400) AID1,TITLESYSTDG(I)
      ENDDO
      AID = 'STDG:LN2'
      READ (NPT,*,ERR=400)
      READ (NPT,*,ERR=400)
      READ (NPT,*,ERR=400)AID1,SYSTDGC, N2BNDC, DOBNDC, TDGTAC;SYSTDGC=ADJUSTR(SYSTDGC); N2BNDC=ADJUSTR(N2BNDC); DOBNDC=ADJUSTR(DOBNDC); TDGTAC=ADJUSTR(TDGTAC)
              
      AID = 'STDG:LN3'
      READ (NPT,*,ERR=400)
      READ (NPT,*,ERR=400)
      DO IG=1,NGT
      READ (NPT,*,ERR=400) AID1,GTTYP(IG), GTPC(IG)
      ENDDO
      GTTYP=ADJUSTR(GTTYP)

      DO IG = 1, NGT
         IF (GTTYP(IG)=='     POW') POWNO = POWNO + 1
         IF (GTTYP(IG)=='     FLD') FLNO  = FLNO  + 1
         IF (GTTYP(IG)=='      RO')  NRO   = NRO  + 1
         IF (GTTYP(IG)=='     SPB') NBAY  = NBAY  + 1
      END DO       
      AID = 'STDG:LN4'
      READ (NPT,*,ERR=400)
      READ (NPT,*,ERR=400)

      READ (NPT,*,ERR=400) AID1,FBE, TWCE, TWEMOD, TWE, TWETSC, TDGLOC, QSPILL, TDGSPMN;  TWETSC=ADJUSTR(TWETSC); TDGLOC=ADJUSTR(TDGLOC)
      AID = 'STDG:LN5'
      READ (NPT,*,ERR=400)
      READ (NPT,*,ERR=400)
      READ (NPT,*,ERR=400) AID1,TDGEQ, TDGP1, TDGP2, TDGP3, TDGP4, TDGP12, TDGP22, TDGP32, TDGP42   
      AID = 'STDG:LN6'
      IF (NRO>0) READ (NPT,*,ERR=400) AID1,TDGROEQ, ROP1, ROP2, ROP3, ROP4     
      AID = 'STDG:LN7'
      READ (NPT,*,ERR=400)
      READ (NPT,*,ERR=400)

      READ (NPT,*,ERR=400) AID1,TDGENTEQ, TDGE1, TDGE2, TDGE12, TDGE22   
      AID = 'STDG:LN8'
      READ (NPT,*,ERR=400)
      READ (NPT,*,ERR=400)

      READ (NPT,*,ERR=400) AID1, TWEFN                                                             
      CLOSE(NPT)       
    
    ELSE
      AID = 'STDG:LN1'
      READ (NPT, '(///(8X,A72))',ERR=400) (TITLESYSTDG(i), i=1,10)
      AID = 'STDG:LN2'
      READ (NPT, '(//8x,4A8)',ERR=400)SYSTDGC, N2BNDC, DOBNDC, TDGTAC
      AID = 'STDG:LN3'
      READ (NPT,'(//(:8X,A8,F8.2))',ERR=400) (GTTYP(IG), GTPC(IG), IG=1,NGT)   

      DO IG = 1, NGT
         IF (GTTYP(IG)=='     POW') POWNO = POWNO + 1
         IF (GTTYP(IG)=='     FLD') FLNO  = FLNO  + 1
         IF (GTTYP(IG)=='      RO')  NRO   = NRO  + 1
         IF (GTTYP(IG)=='     SPB') NBAY  = NBAY  + 1
      END DO                                                                                             

      AID = 'STDG:LN4'
      READ (NPT,'(//8X,2F8.3,I8,F8.3,2A8,2F8.3)',ERR=400) FBE, TWCE, TWEMOD, TWE, TWETSC, TDGLOC, QSPILL, TDGSPMN  
      AID = 'STDG:LN5'
      READ (NPT,'(//8X,I8,8F8.3)',ERR=400) TDGEQ, TDGP1, TDGP2, TDGP3, TDGP4, TDGP12, TDGP22, TDGP32, TDGP42   
      AID = 'STDG:LN6'
      IF (NRO>0) READ (NPT,'(8X,I8,4F8.5)',ERR=400) TDGROEQ, ROP1, ROP2, ROP3, ROP4     
      AID = 'STDG:LN7'
      READ (NPT,'(//8X,I8,4F8.3)',ERR=400) TDGENTEQ, TDGE1, TDGE2, TDGE12, TDGE22   
      AID = 'STDG:LN8'
      READ (NPT,'(//(8X,A72))',ERR=400)  TWEFN                                                             
      CLOSE(NPT)       
    ENDIF
        IF(SYSTDGC /= '     OFF' .AND. SYSTDGC /= '      ON')THEN
        CALL ERRORS
        WRITE (ERR,'(A,A8)') 'w2_systdg.npt file: SYSTDGC must be either OFF or ON. SYSTDGC=',SYSTDGC
        ENDIF
        IF(N2BNDC /= '     OFF' .AND. N2BNDC /= '      ON')THEN
        CALL ERRORS
        WRITE (ERR,'(A,A8)') 'w2_systdg.npt file: N2BNDC must be either OFF or ON. N2BNDC=',N2BNDC
        ENDIF
        IF(DOBNDC /= '     OFF' .AND. DOBNDC /= '      ON')THEN
        CALL ERRORS
        WRITE (ERR,'(A,A8)') 'w2_systdg.npt file: DOBNDC must be either OFF or ON. DOBNDC=',DOBNDC
        ENDIF
        IF(TDGTAC /= '     OFF' .AND. TDGTAC /= '      ON')THEN
        CALL ERRORS
        WRITE (ERR,'(A,A8)') 'w2_systdg.npt file: TDGTAC must be either OFF or ON. TDGTAC=',TDGTAC
        ENDIF

  ENDIF
  
  ! MULTIPLE WATERBODY input file 4.2
   WRITE (WIN,*) 'Multiple Waterbody Parallel simulation input file checks'
  MULTIPLEWB=.FALSE.    
  INQUIRE(FILE='multiple_WB.npt',EXIST=MULTIPLEWB)   
  IF(MULTIPLEWB)THEN
      CALL ERRORS
        WRITE (ERR,'(A,A2)') 'Rename multiple_WB.npt-->w2_multiple_WB.npt. File name convention was changed in Version 4.3. Add prefix: w2_.'
        SYSTDG=.FALSE.
  ENDIF
  
  INQUIRE(FILE='w2_multiple_WB.npt',EXIST=MULTIPLEWB)    
  IF(MULTIPLEWB)THEN
   WRITE (WIN,*) '  Multiple Waterbody Parallel simulation ON: Reading w2_multiple_WB.npt'

    OPEN (NPT, FILE='w2_multiple_WB.npt', STATUS='OLD')
     AID = 'MuWB:LN1'
     READ(NPT,*,ERR=400)
     AID = 'MuWB:LN2'
     READ(NPT,*,ERR=400)MWB
        IF(MWB /= 'OF' .AND. MWB /= 'ON')THEN
        CALL ERRORS
        WRITE (ERR,'(A,A2)') 'w2_multiple_WB.npt file: Multiple WB Wait must be either OFF or ON. Multiple WB wait=',MWB
        ENDIF    
     
     AID = 'MuWB:LN3'     
     READ(NPT,*,ERR=400)
     AID = 'MuWB:LN4'
     READ(NPT,*,ERR=400)IMWB
     AID = 'MuWB:LN5'
      IF(IMWB <= 0)THEN
        CALL ERRORS
        WRITE (ERR,'(A,A2)') 'w2_multiple_WB.npt file: Number of input types must be at least 1. Number of input types=',IMWB
      ENDIF   
     
     AID = 'MuWB:LN6'     
     READ(NPT,*,ERR=400)
     IF(IMWB <= 0)IMWB=1
     DO J=1,IMWB
     READ(NPT,*,ERR=400)MWB2(J), INUMMWB(J),FILEDIRMWB(J)
      IF(MWB2(J) /= 'TR' .AND. MWB2(J) /= 'BR')THEN
        CALL ERRORS
        WRITE (ERR,'(A,I2,A,A2)') 'w2_multiple_WB.npt file: Branch inflow or Trib inflow must use either TR or BR. Number of input:',J,' Input designator:',MWB2(J)
      ENDIF   

     ENDDO
     AID = 'MuWB:LN7'
     READ(NPT,*,ERR=400)
     AID = 'MuWB:LN8'
     READ(NPT,*,ERR=400)TIME_BUFFER
        IF(TIME_BUFFER>10.)THEN
        CALL WARNINGS
        WRITE (WRN,'(A,A8)') 'w2_multiple_WB.npt file: Time buffer > 10 days. This may be too long a buffer. Time_Buffer(days)=',time_buffer
        ENDIF    
     
          AID = 'MuWB:LN8'
     READ(NPT,*,ERR=400)
          AID = 'MuWB:LN9'
     READ(NPT,*,ERR=400)WAIT_TIME         
     
        IF(WAIT_TIME>120.)THEN
        CALL WARNINGS
        WRITE (WRN,'(A,A8)') 'w2_multiple_WB.npt file: Wait time between runs > 120 s. This may be too long to wait. Wait_time(s)=',wait_time
        ENDIF    

     CLOSE(NPT)                                                    
  ENDIF
! ATMOSPHERIC DEPOSITION CHECKS
     UNIT=NPT
   WRITE (WIN,*) 'Atmospheric Deposition checks'

    IF(CONSTITUENTS)THEN
        ATM_DEP_LOADING=0.0
        DO JW=1,NWB
        IF(ATM_DEPOSITIONC(JW)=='      ON')THEN
            ATM_DEP_AVG(:,JW)=0.0
            ATM_DEP_MIN(:,JW)=0.0
            ATM_DEP_MAX(:,JW)=0.0
            NN=0
          WRITE (WIN,*) '  Atmospheric Deposition ON: Reading Atm Dep input file for WB:',JW
         IF(ATM_DEPOSITION_INTERPOLATION(JW) /= '      ON'.AND. ATM_DEPOSITION_INTERPOLATION(JW) /= '     OFF')THEN
          CALL ERRORS
          WRITE (ERR,'(A,I4,A,A)') 'ATM DEP:ATM_DEPOSITION_INTERPOLATION must be ON or OFF in control file FOR WB:',JW,'. ATM_DEPOSITION_INTERPOLATION(JW)=',ATM_DEPOSITION_INTERPOLATION(JW)
         ENDIF
         
          OPEN (NPT,FILE=ATMDEPFN(JW),STATUS='OLD')
          
          READ( NPT,'(A1)')INFORMAT1
          READ (NPT,'(/)')
          
          DO JAC=1,NACATD(JW)
          ATM_DEP_MIN(ATMDCN(JAC,JW),JW)=100000.
          ENDDO

          DO J=1,5000
         IF(INFORMAT1=='$')THEN
            READ (NPT,*,END=4509,ERR=400)  JDAY,(ATM_DEP_LOADING(ATMDCN(JAC,JW),JW),JAC=1,NACATD(JW))         
          ELSE
            READ (NPT,'(100F8.0)',END=4509,ERR=400) JDAY,(ATM_DEP_LOADING(ATMDCN(JAC,JW),JW),JAC=1,NACATD(JW)) 
          ENDIF
! ATM DEP CHECKS
        IF (J == 1 .AND. JDAY > TMSTRT) THEN
          CALL ERRORS
          WRITE (ERR,FMTF) 'ATM DEP:Starting date [',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//ATMDEPFN(JW)
        ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
          if(jday /= 0.0)then
          CALL ERRORS
          WRITE (ERR,FMTF) 'ATM DEP:Julian date ',JDAY,' <= previous date of ',JDAYO,' in '//ATMDEPFN(JW)
          else
          CALL ERRORS
          WRITE (ERR,FMTF) '[Note: This may be a result of end of file blank lines]:Julian date ',JDAY,' <= previous date of ',JDAYO,' in '//ATMDEPFN(JW)
          endif
        END IF
        NN=NN+1
        DO JAC=1,NACATD(JW)
        ATM_DEP_AVG(ATMDCN(JAC,JW),JW)=ATM_DEP_AVG(ATMDCN(JAC,JW),JW)+ATM_DEP_LOADING(ATMDCN(JAC,JW),JW)
        ATM_DEP_MIN(ATMDCN(JAC,JW),JW)=MIN(ATM_DEP_MIN(ATMDCN(JAC,JW),JW),ATM_DEP_LOADING(ATMDCN(JAC,JW),JW))
        ATM_DEP_MAX(ATMDCN(JAC,JW),JW)=MAX(ATM_DEP_MAX(ATMDCN(JAC,JW),JW),ATM_DEP_LOADING(ATMDCN(JAC,JW),JW))
        ENDDO
          
          
          JDAYO=JDAY
          ENDDO
          
4509      CLOSE(NPT)
        ATM_DEP_AVG(:,JW)=ATM_DEP_AVG(:,JW)/REAL(NN)
        
            !DO JB=BS(JW),BE(JW)
            !    DO I=CUS(JB),DS(JB)
            !        DO JAC=1,NACATD(JW)
            !                IF(ATMDCN(JAC,JW)==NPO4.OR.ATMDCN(JAC,JW)==NLPOMP.OR.ATMDCN(JAC,JW)==NRPOMP)THEN 
            !                ATMDEP_P(JW)=ATMDEP_P(JW)+ATM_DEP_LOADING(ATMDCN(JAC,JW),JW)*BI(KT,I)*DLX(I)*3.17098E-11*DLT/1000.    ! P MASS BALANCE IN KG - CUMULATIVE
            !                ELSEIF(ATMDCN(JAC,JW)==NNO3.OR.ATMDCN(JAC,JW)==NLPOMN.OR.ATMDCN(JAC,JW)==NRPOMN.OR.ATMDCN(JAC,JW)==NNH4)THEN 
            !                ATMDEP_N(JW)=ATMDEP_N(JW)+ATM_DEP_LOADING(ATMDCN(JAC,JW),JW)*BI(KT,I)*DLX(I)*3.17098E-11*DLT/1000.    ! N MASS BALANCE IN KG - CUMULATIVE
            !                ENDIF
            !        ENDDO
            !    ENDDO
            !ENDDO  
          
        END IF
        ENDDO
    ENDIF
  
  
!Particle tracking Version 4.1 
  FISH_PARTICLE_EXIST=.FALSE.
  INQUIRE(FILE='particle.csv',EXIST=FISH_PARTICLE_EXIST)   
  IF(FISH_PARTICLE_EXIST)THEN  
    CALL ERRORS
    WRITE (ERR,'(A)') 'Rename particle.csv-->w2_particle.csv. File name convention was changed in Version 4.3. Add prefix: w2_.'
    FISH_PARTICLE_EXIST=.FALSE.
  ENDIF
  
  INQUIRE(FILE='w2_particle.csv',EXIST=FISH_PARTICLE_EXIST)    
  IF(FISH_PARTICLE_EXIST)THEN
   UNIT=NPT
   WRITE (WIN,*) 'Particle Tracking file checks'
   WRITE (WIN,*) '  Particle Tracking ON: Reading w2_particle.csv'
 
    open(NPT,file='w2_particle.csv',status='old')
    READ(NPT,*,END=4510,ERR=400)
    READ(NPT,*,END=4510,ERR=400)
    READ(NPT,*,END=4510,ERR=400)PARTON,NFISHSEG,NFISHPCEL,ALINE,DXTHEORYC,OUTFREQP,ILINEAR,HTSTBOT,HTSTSIDE,IDEBUG   !,ALPHAX,ALPHAZ     !'(/10x,7x,a3,i10,i10,f10.0,(7x,a3),1f10.0,7X,A3)
    
    IF(PARTON /= 'OF' .AND. PARTON /= 'ON')THEN
    CALL ERRORS
    WRITE (ERR,*) 'w2_particle.csv file: PART CONTROL must be either OFF or ON. PARTON=',PARTON
    ENDIF
    
    IF(PARTON=='ON')THEN
        IF(NFISHSEG <= 0.0)THEN
        CALL ERRORS
        WRITE (ERR,'(A,I4)') 'w2_particle.csv file: PART CONTROL is ON. NFISHSEG must be > 0. NFISHSEG=',NFISHSEG
        ENDIF
        IF(NFISHPCEL <= 0.0)THEN
        CALL ERRORS
        WRITE (ERR,'(A,I4)') 'w2_particle.csv file: PART CONTROL is ON. NFISHPCEL must be > 0. NFISHPCEL=',NFISHPCEL
        ENDIF
        IF(ALINE /= 'OF' .AND. ALINE /= 'ON')THEN
        CALL ERRORS
        WRITE (ERR,'(A,A)') 'w2_Particle.csv file: PART CONTROL is ON. LINE MUST BE EITHER ON OR OFF. LINE=',ALINE
        ENDIF
        IF(DXTHEORYC /= 'OF' .AND. DXTHEORYC /= 'ON')THEN
        CALL ERRORS
        WRITE (ERR,'(A,A)') 'w2_Particle.csv file: PART CONTROL is ON. DXTHEORY MUST BE EITHER ON/OFF. DXTHEORY=',DXTHEORYC
        ENDIF
        IF(OUTFREQP <= 0.0)THEN
        CALL ERRORS
        WRITE (ERR,'(A,F8.3)') 'w2_Particle.csv file: PART CONTROL is ON. OUTFREQP must be > 0. OUTFREQP=',OUTFREQP
        ENDIF
        IF(ILINEAR /= 1 .AND. ILINEAR /= 0)THEN
        CALL ERRORS
        WRITE (ERR,'(A,I4)') 'w2_Particle.csv file: PART CONTROL is ON. ILINEAR must be 1 OR 0. ILINEAR=',ILINEAR
        ENDIF
        IF(HTSTBOT /= 1 .AND. HTSTBOT /= 0)THEN
        CALL ERRORS
        WRITE (ERR,'(A,I4)') 'w2_Particle.csv file: PART CONTROL is ON. HTSTBOT must be 1 OR 0. HTSTBOT=',HTSTBOT
        ENDIF
        IF(HTSTSIDE /= 1 .AND. HTSTSIDE /= 0)THEN
        CALL ERRORS
        WRITE (ERR,'(A,I4)') 'w2_Particle.csv file: PART CONTROL is ON. HTSTSIDE must be 1 OR 0. HTSTSIDE=',HTSTSIDE
        ENDIF
        IF(IDEBUG /= 1 .AND. IDEBUG /= 0)THEN
        CALL ERRORS
        WRITE (ERR,'(A,I4)') 'w2_Particle.csv file: PART CONTROL is ON. IDEBUG must be 1 OR 0. IDEBUG=',IDEBUG
        ENDIF

        READ(NPT,*)
        DO I=1,NFISHSEG
        READ(NPT,*,END=4510,ERR=400)J,IFISH,IFISHT,IFISHB,FXLOCINIT,FZLOCINIT,SEDVELINIT,DELAYDATEINIT
        JWW=1
        DO JW=1,NWB
            DO JB=BS(JW),BE(JW)
            IF(IFISH==US(JB)-1  .OR. IFISH==DS(JB)+1)THEN
                CALL ERRORS
                WRITE (ERR,'(A,I4,A,I4,A,I4)') 'Particle.csv file: PARTON is ON. IFISH is a bouundary inactive segment. It must be between US(JB) and DS(JB) for JB=',jb,' IFISH=',ifish,' Particle #:',J
            ELSEIF(IFISH >= US(JB)  .AND. IFISH <= DS(JB))THEN
                JWW=JW
                EXIT
            ENDIF
            ENDDO
        ENDDO
        IF(FXLOCINIT > DLX(IFISH) .OR. FXLOCINIT < 0.0)THEN
            CALL ERRORS
            WRITE (ERR,'(A,F8.3,A,I4)') 'Particle.csv file: PARTON is ON. FXLOC(M) must be >= 0.0 or < DLX. FXLOC=',FXLOCINIT,' Particle #:',J
        ENDIF
        IF(FZLOCINIT > H(IFISHT,JWW) .OR. FZLOCINIT < 0.0)THEN
            CALL ERRORS
            WRITE (ERR,'(A,F8.3,A,I4)') 'Particle.csv file: PARTON is ON. FZLOC(M) must be >= 0.0 or < DZ. FZLOC=',FZLOCINIT,' Particle #:',J
        ENDIF
        IF(IFISHB > KB(IFISH))THEN
            CALL ERRORS
            WRITE (ERR,'(A,I4,A,I4,A,I4)') 'Particle.csv file: PARTON is ON. IFISHB must be <= KB[KB=',KB(IFISH),']. IFISHB=',IFISHB,' Particle #:',J
        ENDIF
        IF(IFISHB < KTWB(JWW))THEN
            CALL ERRORS
            WRITE (ERR,'(A,I4,A,I4,A,I4)') 'Particle.csv file: PARTON is ON. IFISHB must be >= KT[KT=',KTWB(JWW),']. IFISHB=',IFISHB,' Particle #:',J
        ENDIF
        IF(DELAYDATEINIT > tmend)THEN
        CALL ERRORS
        WRITE(ERR,'(A,F10.3,A,F10.3,A,I4)')'Particle.csv file: PARTON is ON. DELAYDATEINIT must be < End Model Run Time[TMEND=',tmend,']. DELAYDATEINIT=',DELAYDATEINIT,' Particle #:',J
        ENDIF
        IF(SEDVELINIT > 25.)THEN
        CALL WARNINGS
        WRITE(WRN,'(A,F8.3,A,I4)')'Particle.csv file: PARTON is ON. SEDVELINIT > 25 m/day. SEDVELINIT=',SEDVELINIT,' Particle #:',J
       ENDIF

        ENDDO
        READ(NPT,*,END=4510,ERR=400)  
        READ(NPT,*,END=4510,ERR=400)NUMCLASS
        
        IF(NUMCLASS == 0)THEN
        CALL WARNINGS
        WRITE(WRN,'(A)')'Particle.csv file: NUMCLASS = 0. Hence, there will be no histogram output for each particle.'
        ELSEIF(NUMCLASS < 0)THEN
             CALL ERRORS
             WRITE(ERR,'(A,F10.3,A,F10.3,A,I4)')'Particle.csv file: NUMCLASS is < 0. NUMCLASS must be >= 0.'
       ENDIF
        
        READ(NPT,*,END=4510,ERR=400)
        READ(NPT,*,END=4510,ERR=400)VEL_INT,VEL_TOP,TEMP_INT,TEMP_TOP,D_INT,D_TOP
        
        IF(NUMCLASS > 0 .AND. VEL_INT <= 0.0)THEN
        CALL WARNINGS
        WRITE(WRN,'(A)')'Particle.csv file: Velocity Interval (VEL_INT) <= 0.0. This must be > 0.0 to define velocity histogram bins.'
       ENDIF
        IF(NUMCLASS > 0 .AND. TEMP_INT <= 0.0)THEN
        CALL WARNINGS
        WRITE(WRN,'(A)')'Particle.csv file: Temperature Interval (TEMP_INT) <= 0.0. This must be > 0.0 to define temperature histogram bins.'
       ENDIF
        IF(NUMCLASS > 0 .AND. D_INT <= 0.0)THEN
        CALL WARNINGS
        WRITE(WRN,'(A)')'Particle.csv file: Depth Interval (D_INT) <= 0.0. This must be > 0.0 to define depth histogram bins.'
       ENDIF

        
    ENDIF 
    CLOSE(NPT)
 
  ELSE
  WRITE (WIN,*) 'Particle Tracking file checks'
  WRITE (WIN,*) '  Particle Tracking OFF'
    ENDIF
    GO TO 4520
4510 CONTINUE 
     CALL ERRORS
     WRITE (ERR,*) 'Particle.csv file read error'
4520 CONTINUE
! CONVERT LAT/LONG TO UTM PROJECTION NAD83
 ! CALL UTMS
! CONVERT UTM POINT AT DOWNSTREAM POINT TO SEGMENT POLYGONS AND SEGMENT CENTERS IN UTM COORDINATES
 ! CALL POLYGONS 
  
! Boundary cells

  WRITE (WIN,*) '  widths'
  DO JB=1,NBR
    DO K=1,KMX
      IF (B(K,US(JB)-1) /= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMT2IF) 'Upstream boundary segment width [B(',K,',',US(JB)-1,')=',B(K,US(JB)-1),'] /= 0 for branch ',JB
      END IF
      IF (B(K,DS(JB)+1) /= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMT2IF) 'Downstream boundary segment width [B(',K,',',DS(JB)+1,')=',B(K,DS(JB)+1),'] /= 0 for branch ',JB
      END IF
    END DO
    DO I=US(JB),DS(JB)
      DO K=2,KB(I)
        IF (B(K+1,I) > B(K,I) .AND. B(K+1,I) /= 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMT2IF) 'Cell width [B(',K+1,',',I,')=',B(K+1,I),'] > [B(',K,',',I,')=',B(K,I)
        END IF
        IF (B(K,I) < 5.0 .AND. B(K,I) > 0.0) THEN                                                                      !SW 01/23/01
          CALL WARNINGS
          WRITE (WRN,FMT2IF) 'Cell width [B(',K,',',I,')=',B(K,I),'] < 5m which can cause stability problems'
        END IF
      END DO
      IF (B(1,I) /= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTIF) 'Surface boundary layer cell width [B(1,',I,')=',B(1,I),'] /= 0'
      ELSE IF (B(KMX,I) /= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMT2IF) 'Bottom boundary layer cell width [B(',KMX,',',I,')=',B(KMX,I),'] /= 0'
      END IF
    END DO
  END DO

! Boundary widths

  DO JB=1,NBR
    DO K=2,KMX
      B(K,US(JB)-1) = B(K,US(JB))
      B(K,DS(JB)+1) = B(K,DS(JB))
    END DO
  END DO

! Width change

  DO JB=1,NBR
    DO I=US(JB),DS(JB)-1
      DO K=2,KB(I)
        IF (B(K,I+1) /= 0.0 .AND. B(K,I) /= 0.0) THEN
          IF (B(K,I)/B(K,I+1) >= 10.0 .OR. B(K,I+1)/B(K,I) >= 10.0) THEN
            CALL WARNINGS
            WRITE (WRN,FMT2IF) 'Cell width [B(',K,',',I,')=',B(K,I),'] < or > 10X width [B(',K,',',I+1,')=',B(K,I+1),']'

          END IF
        END IF
      END DO
    END DO
  END DO

! Layer thickness

  WRITE (WIN,*) '  layer thickness'
  DO JW=1,NWB
    DO K=1,KMX
      IF (H(K,JW) <= 0) THEN
        CALL   ERRORS
        WRITE (ERR,FMTF2I) 'Layer thickness [H =',H(K,JW),'] <= 0.0 m for layer ',K,' in waterbody ',JW
      ELSE IF (H(K,JW) < 0.1) THEN
        CALL   WARNINGS
        WRITE (WRN,FMTF2I) 'Layer thickness [H =',H(K,JW),'] <  0.1 m for layer ',K,' in waterbody ',JW
      ELSE IF (H(K,JW) > 3.0) THEN
        CALL   WARNINGS
        WRITE (WRN,FMTF2I) 'Layer thickness [H =',H(K,JW),'] >  3.0 m for layer ',K,' in waterbody ',JW
      END IF
    END DO
  END DO

! Meteorologic data

  sromax=0.0
  sromin=1000.
  tairmax=0.0
  tairmin=1000.
  tdewmax=0.0
  tdewmin=1000.
  windmax=0.0
  windmin=1000.
  phimin=1000.
  phimax=0.0
  cloudmax=0.0
  cloudmin=10.
  sroav=0.0
  tairav=0.0
  phiav=0.0
  windav=0.0
  phiav=0.0
  cloudav=0.0
  tdewav=0.0
  
  WRITE (WIN,*) 'Input files'
  WRITE (WIN,*) '  meteorology'
  DO JW=1,NWB
    write(win,*) '    ',metfn(jw)
    OPEN (UNIT=NPT,FILE=METFN(JW),STATUS='OLD',IOSTAT=IERR)
    IF (IERR == 0) THEN
        AID=METFN(JW)
         READ(NPT,'(A1)')ICHAR1
         IF(ICHAR1=='$')WRITE(WIN,*) '      met file in csv format'
        READ (NPT,'(/)',END=401)
   !   READ (NPT,'(//)')
        
        ! CHECK IF ERROR -CSV FORMAT CHECK
        IF(ICHAR1 /= '$')THEN
        READ(NPT,'(A72)')TITLE(1)
            DO JJ=1,50    ! CHECK ONLY FIRST 50 CHARACTERS
                IF(TITLE(1)(JJ:JJ) == ",")THEN
                    CALL ERRORS
                    WRITE (ERR,'(A175)') 'Meteorological file: File not designated as a CSV file (no "$" as first character) but commas found between input values for '//METFN(JW)
                    EXIT
                ENDIF
            ENDDO
        BACKSPACE(NPT)
        ENDIF
        
      J     = 1
      JDAYO = 0.0
      nmet=0
      
      DO WHILE (.TRUE.)
        IF(SROC(JW) == '      ON')then
        IF(ICHAR1=='$')THEN
        READ (NPT,*,END=110,ERR=400) JDAY, TAIR, TDEW, WIND, PHI, CLOUD, SRO
        ELSE
        READ (NPT,'(10F8.0/8X,9F8.0)',END=110,ERR=400) JDAY, TAIR, TDEW, WIND, PHI, CLOUD, SRO
        ENDIF
        nmet=nmet+1
        sromax(jw)=max(sro,sromax(jw))
        sromin(jw)=min(sro,sromin(jw))
        sroav(jw)=sroav(jw)+sro
        tairmax(jw)=max(tair,tairmax(jw))
        tairmin(jw)=min(tair,tairmin(jw))
        tairav(jw)=tairav(jw)+tair
        phimax(jw)=max(phi,phimax(jw))
        phimin(jw)=min(phi,phimin(jw))
        phiav(jw)=phiav(jw)+phi
        cloudmax(jw)=max(cloud,cloudmax(jw))
        cloudmin(jw)=min(cloud,cloudmin(jw))
        cloudav(jw)=cloudav(jw)+cloud
        windmax(jw)=max(wind,windmax(jw))
        windmin(jw)=min(wind,windmin(jw))
        windav(jw)=windav(jw)+wind
        tdewmax(jw)=max(tdew,tdewmax(jw))
        tdewmin(jw)=min(tdew,tdewmin(jw))
        tdewav(jw)=tdewav(jw)+tdew         
        
        IF (SRO < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTF) 'Solar Radiation [SRO=',SRO,'] < 0 W/m2 on day ',JDAY,' in '//METFN(JW)
        END IF
        IF (SRO > 1500.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTF) 'Solar Radiation [SRO=',SRO,'] > 1500 W/m2 (Solar Radiation at top of atmosphere == 1360 W/m2) on day ',JDAY,' in '//METFN(JW)
        ELSE IF (SRO > 1360.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTF) 'Solar Radiation [SRO=',SRO,'] > 1360 W/m2 (Solar Radiation at top of atmosphere == 1360 W/m2) on day ',JDAY,' in '//METFN(JW)
        END IF
        else
        IF(ICHAR1=='$')THEN
        READ (NPT,*,END=110,ERR=400) JDAY, TAIR, TDEW, WIND, PHI, CLOUD
        ELSE
        READ (NPT,'(10F8.0/8X,9F8.0)',END=110,ERR=400) JDAY, TAIR, TDEW, WIND, PHI, CLOUD
        ENDIF
        nmet=nmet+1
        tairmax(jw)=max(tair,tairmax(jw))
        tairmin(jw)=min(tair,tairmin(jw))
        tairav(jw)=tairav(jw)+tair
        phimax(jw)=max(phi,phimax(jw))
        phimin(jw)=min(phi,phimin(jw))
        phiav(jw)=phiav(jw)+phi
        cloudmax(jw)=max(cloud,cloudmax(jw))
        cloudmin(jw)=min(cloud,cloudmin(jw))
        cloudav(jw)=cloudav(jw)+cloud
        windmax(jw)=max(wind,windmax(jw))
        windmin(jw)=min(wind,windmin(jw))
        windav(jw)=windav(jw)+wind
        tdewav(jw)=tdewav(jw)+tdew
        tdewmin(jw)=min(tdew,tdewmin(jw))
        tdewmax(jw)=max(tdew,tdewmax(jw))
        endif
        IF (J == 1 .AND. JDAY > TMSTRT) THEN
          CALL ERRORS
          WRITE (ERR,FMTF) 'Starting date [',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//METFN(JW)
        ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
          if(jday /= 0.0)then
          CALL ERRORS
          WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date of ',JDAYO,' in '//METFN(JW)
          else
          CALL ERRORS
          WRITE (ERR,FMTF) '[Note: This may be a result of end of file blank lines]:Julian date ',JDAY,' <= previous date of ',JDAYO,' in '//METFN(JW)
          endif
        END IF
        IF(IEEE_IS_NAN(JDAY))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'JDAY is NaN on line of data ',J,' in '//METFN(JW)
        ENDIF
        IF(IEEE_IS_NAN(TAIR))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'TAIR is NaN on Julian date ',JDAY,' in '//METFN(JW)
        ENDIF
        IF(IEEE_IS_NAN(TDEW))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'TDEW is NaN on Julian date ',JDAY,' in '//METFN(JW)
        ENDIF
        IF(IEEE_IS_NAN(WIND))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'WIND is NaN on Julian date ',JDAY,' in '//METFN(JW)
        ENDIF        
        IF(IEEE_IS_NAN(PHI))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'PHI is NaN on Julian date ',JDAY,' in '//METFN(JW)
        ENDIF
        IF(IEEE_IS_NAN(CLOUD))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'CLOUD is NaN on Julian date ',JDAY,' in '//METFN(JW)
        ENDIF        
        
        IF (TAIR > 50.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTF) 'Air temperature [TAIR=',TAIR,'] > 50 deg on day ',JDAY,' in '//METFN(JW)
        END IF
        IF (TDEW > 50.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTF) 'Dew point temperature [TDEW=',TDEW,'] > 50 deg on day ',JDAY,' in '//METFN(JW)
        END IF
        IF (WIND > 20.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTF) 'Wind speed [WIND= ',WIND,'] > 20 m/s on day ',JDAY,' in '//METFN(JW)
        END IF
        IF (WIND < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTF) 'Wind speed [WIND=',WIND,'] < 0 m/s on day ',JDAY,' in '//METFN(JW)
        END IF
        IF (CLOUD > 10.0 .OR. CLOUD < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTF) 'Cloud cover [CLOUD=',CLOUD,'] < 0 or > 10 on day ',JDAY,' in '//METFN(JW)
        END IF
        IF (TAIR == 0.0 .AND. TDEW == 0.0 .AND. ABS(TAIR-TAIRO) > 5.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTF) 'Possible missing value on day ',JDAY,' in '//METFN(JW)
        END IF
        if(jday-jdayo >= 0.5 .and. SROC(JW) == '      ON' .AND. J /= 1)then
            CALL ERRORS
            WRITE (ERR,'(a,f10.3,a,f9.3,a,a)') 'MET Error: Inappropriate interpolation of solar radiation over time steps longer than a few hours. Difference between Jday values is',(JDAY-jdayo),' at JDAY:',JDAY,' in '//METFN(JW)
        elseif(jday-jdayo >= 1.0 .AND. J /= 1)then
            CALL WARNINGS
            WRITE (WRN,'(a,f10.3,a,f9.3,a,a)') 'MET Warning: Interpolation of met data over longer than a few hours is not recommended. Difference between Jday values is',(JDAY-jdayo),' at JDAY:',JDAY,' in '//METFN(JW)
        END IF
        J     = J+1
        JDAYO = JDAY
        TAIRO = TAIR
      END DO
 110  CONTINUE
       IF(SROC(JW) == '      ON')then
           sroav(jw)=sroav(jw)/nmet
       else
           sromin(jw)=0.0
           sromax(jw)=0.0
       endif
       
       tairav(jw)=tairav(jw)/nmet
       tdewav(jw)=tdewav(jw)/nmet
       windav(jw)=windav(jw)/nmet
       phiav(jw)=phiav(jw)/nmet
       cloudav(jw)=cloudav(jw)/nmet   
       
        if(cloudav(jw) < 1.0)then
            CALL ERRORS
            WRITE (ERR,'(a,f7.3,a,i3,a,a)') 'MET Error: Possible inappropriate cloud cover since average [CloudCoverAverage:',cloudav(jw),'] < 1.0 [Note that CLOUD should be between 0-10] for WB:',JW,' in '//METFN(JW)
        elseif(cloudav(jw) < 3.0)then
            CALL WARNINGS
            WRITE (WRN,'(a,f7.3,a,i3,a,a)') 'MET Warning: Check cloud cover since average [CloudCoverAverage:',cloudav(jw),'] < 3.0 [Note that CLOUD should be between 0-10] for WB:',JW,' in '//METFN(JW)
        END IF
        
        if(tairav(jw) < tdewav(jw))then
            CALL ERRORS
            WRITE (ERR,'(a,f7.2,a,f7.2,a,i3,a,a)') 'MET Error: Tair-average is less than Tdew-average. [Tair-av(C):',tairav(jw),'< Tdew-av(C):',tdewav(jw),'] for WB:',JW,' in '//METFN(JW)
        endif 
       
      
      IF (JDAY < TMEND) THEN
        CALL ERRORS
        WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//METFN(JW)
      END IF
      CLOSE (NPT)
    ELSE
      CALL ERRORS
      WRITE (ERR,FMTI) 'Could not open waterbody ',JW,' meteorological file '//METFN(JW)
    END IF
  END DO
  DO JW=1,NWB
  if(SROC(JW) == '      ON' .and. sromax(jw) < 150.0)then
      CALL ERRORS
      WRITE (ERR,*) 'SROC is ON but the maximum solar radiation is less than 150 W/m2. Check SRO in meteorological file:', METFN(JW)
  endif
   if(SROC(JW) == '      ON' .and. sromax(jw) < 400.0)then
      CALL WARNINGS
      WRITE (WRN,*) 'SROC is ON but the maximum solar radiation is less than 400 W/m2. Check SRO in meteorological file:', METFN(JW)
  endif
  ENDDO


! Number of records for inflow files

  DO JB=1,NBR
    OPEN (UNIT=NPT,FILE=QINFN(JB),STATUS='OLD',IOSTAT=IERR)
    IF (IERR == 0) THEN
      DO J=1,50000000
        READ (NPT,*,END=130)
      END DO
130   CONTINUE
      NQIN(JB) = J-4
      NQINMX   = MAX(NQINMX,NQIN(JB))+5
      CLOSE (NPT)
    END IF
    OPEN (UNIT=NPT,FILE=QDTFN(JB),STATUS='OLD',IOSTAT=IERR)
    IF (IERR == 0) THEN
      DO J=1,50000000
        READ (NPT,*,END=131)
      END DO
131   CONTINUE
      NQDT(JB) = J-4
      NQDTMX   = MAX(NQDTMX,NQDT(JB))+5
      CLOSE (NPT)
    END IF
    OPEN (UNIT=NPT,FILE=PREFN(JB),STATUS='OLD',IOSTAT=IERR)
    IF (IERR == 0) THEN
      DO J=1,50000000
        READ (NPT,*,END=132)
      END DO
132   CONTINUE
      NQPR(JB) = J-4
      NQPRMX   = MAX(NQPRMX,NQPR(JB))+5
      CLOSE (NPT)
    END IF
  END DO
  DO JT=1,NTR
    OPEN (UNIT=NPT,FILE=QTRFN(JT),STATUS='OLD',IOSTAT=IERR)
    IF (IERR == 0) THEN
      DO J=1,50000000
        READ (NPT,*,END=133)
      END DO
133   CONTINUE
      NQTR(JT) = J-3
      NQTRMX   = MAX(NQTRMX,NQTR(JT))+5
      CLOSE (NPT)
    END IF
  END DO
  OPEN (UNIT=NPT,FILE=QWDFN,STATUS='OLD',IOSTAT=IERR)
  IF (IERR == 0) THEN
    DO J=1,50000000
      READ (NPT,*,END=134)
    END DO
134 CONTINUE
    NQWD    = J-4
    NQWDMX  = MAX(NQWDMX,NQWD)+5
    CLOSE (NPT)
  END IF
  ALLOCATE (JDQIN(NQINMX,NBR),QIN(NQINMX,NBR))
  ALLOCATE (JDQDT(NQDTMX,NBR),QDT(NQDTMX,NBR))
  ALLOCATE (JDQPR(NQPRMX,NBR),QPR(NQPRMX,NBR))
  ALLOCATE (JDQTR(NQTRMX,NTR),QTR(NQTRMX,NTR))
  ALLOCATE (JDQWD(NQWDMX,NWD),QWD(NQWDMX,NWD))
  QIN = 0.0; QDT = 0.0; QPR = 0.0; QTR = 0.0; JDQIN = 0.0; JDQDT = 0.0; JDQPR = 0.0; JDQTR = 0.0

! Time varying data

  DO JW=1,NWB
    DO JB=BS(JW),BE(JW)
      WRITE (BRA,'(I0)') JB
      IF (UQ_EXTERNAL(JB)) THEN

!****** Inflows

        WRITE (WIN,*) '  inflow'
        OPEN (UNIT=NPT,FILE=QINFN(JB),STATUS='OLD',IOSTAT=IERR)
        IF (IERR == 0) THEN
            AID=QINFN(JB)
       !   READ  (NPT,'(A)') HEADER

          J       = 1
          JDAYO   = 0.0
          QINO    = 0.0

        READ(NPT,'(A1)')ICHAR1
        IF(ICHAR1=='$')WRITE(WIN,*) '      QIN file in csv format'
        READ (NPT,'(/)',END=401)
        
        ! CHECK IF ERROR -CSV FORMAT CHECK
        IF(ICHAR1 /= '$')THEN
        READ(NPT,'(A72)')TITLE(1)
            DO JJ=1,50    ! CHECK ONLY FIRST 50 CHARACTERS
                IF(TITLE(1)(JJ:JJ) == ",")THEN
                    CALL ERRORS
                    WRITE (ERR,'(A175)') 'QIN file: File not designated as a CSV file (no "$" as first character) but commas found between input values for '//QINFN(JB)
                    EXIT
                ENDIF
            ENDDO
        BACKSPACE(NPT)
        ENDIF       
        
          DO WHILE (.TRUE.)
                IF(ICHAR1=='$')THEN
            READ (NPT,*,END=135,ERR=400) JDQIN(J,JB), QIN(J,JB)
                  ELSE
            READ (NPT,'(10F8.0/8X,9F8.0)',END=135,ERR=400) JDQIN(J,JB), QIN(J,JB)
                  ENDIF
                  
        IF(IEEE_IS_NAN(JDQIN(J,JB)))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'JDAY is NaN on line of data  ',J,' in '//QINFN(JB)
        ENDIF
        IF(IEEE_IS_NAN(QIN(J,JB)))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'QIN is NaN on Julian date ',JDQIN(J,JB),' in '//QINFN(JB)
        ENDIF
             

            IF (QIN(J,JB)<0.0) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'QIN cannot be negative:',QIN(J,JB),' at date [JDAY=',JDQIN(J,JB),'] in '//QINFN(JB)
            END IF

            IF (JDAYO >= TMSTRT .AND. JDAYO <= TMEND) THEN
              QINMX(JB) = MAX(QIN(J,JB),QINMX(JB))
              IF (JDQIN(J,JB) > TMEND) THEN
                QINS(JB) = QINS(JB)+QINO*(TMEND-JDAYO)
              ELSE
                QINS(JB) = QINS(JB)+QINO*(JDQIN(J,JB)-JDAYO)
              END IF
            END IF
            IF ((J == 1) .AND. (JDQIN(J,JB) > TMSTRT)) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Starting date [',JDQIN(J,JB),'] > simulation start date [TMSTRT=',TMSTRT,'] in '//QINFN(JB)
            ELSE IF (JDQIN(J,JB) <= JDAYO .AND. J /= 1) THEN
                if(jdqin(j,jb) /= 0.0)then
                CALL ERRORS
              WRITE (ERR,FMTF) 'Julian date ',JDQIN(J,JB),' <= previous date ',JDAYO,' in '//QINFN(JB)
                else
                CALL ERRORS
              WRITE (ERR,FMTF) '[Note: This could be a result of blank lines at the end of the file]: Julian date ',JDQIN(J,JB),' <= previous date ',JDAYO,' in '//QINFN(JB)
                endif
            END IF
            QINO  = QIN(J,JB)
            JDAYO = JDQIN(J,JB)
            J     = J+1
          END DO
 135      CONTINUE
          IF (JDQIN(J-1,JB) < TMEND) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Ending time ',JDQIN(J-1,JB),' < ending simulation time [TMEND=',TMEND,'] in '//QINFN(JB)
          END IF
! start of check to make sure flow /= 0 when intitial water surface and flow turned on   cb 9/7/2010
          IF(INITUWL == '      ON' .and.  slope(jb) /= 0.0)THEN
            rewind(npt)

            READ(NPT,'(A1)')ICHAR1
            READ (NPT,'(/)')

            IF(ICHAR1=='$')THEN
            READ (npt, *,end=137) NXQIN2,QINNX
            ELSE
            READ (npt, '(2F8.0)',end=137) NXQIN2,QINNX
            ENDIF
            QIND   = QINNX
            QINO   = QINNX
            IF(ICHAR1=='$')THEN
            READ (npt, *,end=137) NXQIN1,QINNX
                ELSE
            READ (npt, '(2F8.0)',end=137) NXQIN1,QINNX
            ENDIF
            jday=tmstrt
            DO WHILE (JDAY >= NXQIN1)
              QIND   = QINNX
              QINO   = QINNX
              NXQIN2 = NXQIN1
              IF(ICHAR1=='$')THEN
              READ (npt,*,end=137) NXQIN1,QINNX
              ELSE
              READ (npt,'(2F8.0)',end=137) NXQIN1,QINNX
              ENDIF
            END DO
            IF (INTERP_INFLOW(JB)) THEN
              QRATIO = (NXQIN1-JDAY)/(NXQIN1-NXQIN2)
              QIND                      = (1.0-QRATIO)*QINNX                  +QRATIO*QINO
            end if
            if(qind == 0.0)then
              CALL ERRORS
              WRITE(ERR,'(4(A,I0,A))')'Initial water level and velocity computation (INITUWL) is ON but flow rate in upstream external flow boundary condition of branch ',jb,' is zero at start of simulation'
            end if
          end if
! start of check to make sure flow /= 0 when intitial water surface and flow turned on
 137      CLOSE (NPT)
        ELSE
          DO J=1,IMWB   ! CHECK FOR MULTIPLE WATER BODIES AND INPUT FILE IS NOT READY YET
              IF(MWB2(J) == 'BR')THEN
                  IF(INUMMWB(J)==JB)THEN
                      IOK=1
                      EXIT
                  ENDIF
              ENDIF
          ENDDO
          IF(IOK==0)THEN  
          CALL ERRORS
          WRITE (ERR,FMTI) 'Could not open branch ',JB,' inflow file '//QINFN(JB)
          ELSE
             CALL WARNINGS 
              WRITE (WRN,FMTI) 'Multiple Waterbodies: Note the following FLOW inflow file for branch ',JB,' is not yet available: '//QINFN(JB)
           ENDIF
          IOK=0
        END IF

!****** Inflow temperatures

        WRITE (WIN,*) '  inflow temperature'
        OPEN (UNIT=NPT,FILE=TINFN(JB),STATUS='OLD',IOSTAT=IERR)
        IF (IERR == 0) THEN
            AID=TINFN(JB)
       !   READ  (NPT,'(A)') HEADER
          J       = 1
          JDAYO   = 0.0

          READ(NPT,'(A1)')ICHAR1
          IF(ICHAR1=='$')WRITE(WIN,*) '      TIN file in csv format'
          READ (NPT,'(/)',END=401,ERR=401)
          
        ! CHECK IF ERROR -CSV FORMAT CHECK
        IF(ICHAR1 /= '$')THEN
        READ(NPT,'(A72)')TITLE(1)
            DO JJ=1,50    ! CHECK ONLY FIRST 50 CHARACTERS
                IF(TITLE(1)(JJ:JJ) == ",")THEN
                    CALL ERRORS
                    WRITE (ERR,'(A175)') 'TIN file: File not designated as a CSV file (no "$" as first character) but commas found between input values for '//TINFN(JB)
                    EXIT
                ENDIF
            ENDDO
        BACKSPACE(NPT)
        ENDIF       

          DO WHILE (.TRUE.)
            IF(ICHAR1=='$')THEN
            READ (NPT,*,END=140,ERR=400) JDAY, TIN
                  ELSE
            READ (NPT,'(10F8.0/8X,9F8.0)',END=140,ERR=400) JDAY, TIN
                  ENDIF
                  
        IF(IEEE_IS_NAN(JDAY))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'JDAY is NaN on line of data ',J,' in '//TINFN(JB)
        ENDIF
        IF(IEEE_IS_NAN(TIN))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'TIN is NaN on Julian date ',JDAY,' in '//TINFN(JB)
        ENDIF
                  
                  
                  
            IF (J == 1 .AND. JDAY > TMSTRT) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Starting date [',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in ' //TINFN(JB)
            ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
                if(jday /= 0.0)then
                CALL ERRORS
              WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//TINFN(JB)
                else
                CALL ERRORS
              WRITE (ERR,FMTF) '[Note: This may be a result of blank lines at end of the file]: Julian date ',JDAY,' <= previous date ',JDAYO,' in '//TINFN(JB)
                endif
            END IF
            TINMAX(JB)=MAX(TIN,TINMAX(JB))
            TINMIN(JB)=MIN(TIN,TINMIN(JB))
            IF(TIN > 40. .AND. TIN < 50.)THEN
                CALL WARNINGS
                WRITE(WRN,'(A,F8.2,A,F8.2,A,I4,A)')'Inflow temperature[TIN=',tin,']>40C on JDAY=',JDAY,' for branch:',JB,' in file:'//TINFN(JB)
            ENDIF
            IF(TIN > 50.)THEN
                CALL ERRORS
                WRITE(ERR,'(A,F8.2,A,F8.2,A,I4,A)')'Inflow temperature[TIN=',tin,']>50C on JDAY=',JDAY,' for branch:',JB,' in file:'//TINFN(JB)
            ENDIF
            IF(TIN < -5.)THEN
                    IF(TIN==-99.0)THEN
                CALL WARNINGS
                WRITE(WRN,'(A,F8.1,A,F8.2,A,I4,A)')'Inflow temperature[TIN=',tin,'] may be from model output with a ZERO inflow rate on JDAY=',JDAY,' for branch:',JB,' in file:'//TINFN(JB)
                    ELSE
                CALL ERRORS
                WRITE(ERR,'(A,F10.1,A,F8.2,A,I4,A)')'Inflow temperature[TIN=',tin,']<-5C on JDAY=',JDAY,' for branch:',JB,' in file:'//TINFN(JB)
                    ENDIF
            ENDIF
            JDAYO = JDAY
            J     = J+1
          END DO
 140      CONTINUE
          IF (JDAY < TMEND) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//TINFN(JB)
          END IF
          CLOSE (NPT)
        ELSE
            
           DO J=1,IMWB   ! CHECK FOR MULTIPLE WATER BODIES AND INPUT FILE IS NOT READY YET
              IF(MWB2(J) == 'BR')THEN
                  IF(INUMMWB(J)==JB)THEN
                      IOK=1
                      EXIT
                  ENDIF
              ENDIF
          ENDDO
          IF(IOK==0)THEN  
          CALL ERRORS
          WRITE (ERR,FMTI) 'Could not open branch ',JB,' inflow temperature file '//TINFN(JB)
          ELSE
             CALL WARNINGS 
              WRITE (WRN,FMTI) 'Multiple Waterbodies: Note the following TEMPERATURE inflow file for branch ',JB,' is not yet available: '//TINFN(JB)
           ENDIF
          IOK=0
        END IF

!****** Inflow concentrations

        IF (CONSTITUENTS .AND. NACIN(JB) > 0) THEN
          WRITE (WIN,*) '  inflow concentrations'
          OPEN (UNIT=NPT,FILE=CINFN(JB),IOSTAT=IERR,STATUS='OLD')
          WRITE(WIN,*) '     cin filename:', adjustl(trim(CINFN(jb)))
          jday1=0.0
          jday2=0.0
          IF (IERR == 0) THEN
              AID=CINFN(JB)
              READ(NPT,'(A1)')ICHAR1
              IF(ICHAR1=='$')THEN
              WRITE(WIN,*) '   cin file in csv format'
              ENDIF
              READ (NPT,'(/)',END=401)
              
        ! CHECK IF ERROR -CSV FORMAT CHECK
        IF(ICHAR1 /= '$')THEN
        READ(NPT,'(A72)')TITLE(1)
            DO JJ=1,50    ! CHECK ONLY FIRST 50 CHARACTERS
                IF(TITLE(1)(JJ:JJ) == ",")THEN
                    CALL ERRORS
                    WRITE (ERR,'(A175)') 'CIN file: File not designated as a CSV file (no "$" as first character) but commas found between input values for '//CINFN(JB)
                    EXIT
                ENDIF
            ENDDO
        BACKSPACE(NPT)
        ENDIF            
              
            ! READ (NPT,'(//)')
            J                               = 1
            JDAYO                           = 0.0
            CINMIN(INCN(1:NACIN(JB),JB),JB) = 1.0E10
            DO WHILE (.TRUE.)
              IF(ICHAR1=='$')THEN
              READ (NPT,*,END=150,ERR=400) JDAY, (CIN(INCN(JC,JB)),JC=1,NACIN(JB))
              ELSE
              READ (NPT,'(F8.0,1000(F8.0))',END=150,ERR=400) JDAY, (CIN(INCN(JC,JB)),JC=1,NACIN(JB))
              ENDIF
              
         IF(IEEE_IS_NAN(JDAY))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'JDAY is NaN on line of data ',J,' in '//CINFN(JB)
         ENDIF
        DO JC=1,NACIN(JB)
         IF(IEEE_IS_NAN(CIN(INCN(JC,JB))))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'CIN is NaN on Julian date ',JDAY,' in '//CINFN(JB)
        ENDIF
        ENDDO
                  
              
              
              
              IF (J == 1) TSTART = JDAY
              CINMIN(INCN(1:NACIN(JB),JB),JB) = MIN(CIN(INCN(1:NACIN(JB),JB)),CINMIN(INCN(1:NACIN(JB),JB),JB))
              CINMAX(INCN(1:NACIN(JB),JB),JB) = MAX(CIN(INCN(1:NACIN(JB),JB)),CINMAX(INCN(1:NACIN(JB),JB),JB))
              IF (J > 1) CINAVG(INCN(1:NACIN(JB),JB),JB) =  CINAVG(INCN(1:NACIN(JB),JB),JB)+(CIN(INCN(1:NACIN(JB),JB))             &
                                                            +CIN2(INCN(1:NACIN(JB),JB)))*0.5*(JDAY-TDAY2)

! compute mass loading from inflows
             if(j==1)then
             UNIT=NPT+1
             open(npt+1,file=qinfn(jb),status='old')
             AID=QINFN(JB)
           !  READ  (NPT+1,'(A)') HEADER
             QINO    = 0.0
             jj=1
             READ(NPT+1,'(A1)')ICHAR2
             READ (NPT+1,'(/)')
              IF(ICHAR2=='$')THEN
              READ (NPT+1,*,ERR=402) jday2,qin2    !JDQIN(JJ,JB), QIN(JJ,JB)
              ELSE
              READ (NPT+1,'(10F8.0/8X,9F8.0)',ERR=402) jday2,qin2    !JDQIN(JJ,JB), QIN(JJ,JB)
              ENDIF
              jday1=jday2
              qin1=qin2
             endif
              !  do while(jdqin(jj,jb) < jday)
                AID=QINFN(JB)
                do while(jday2 < jday)
                jj=jj+1
                IF(ICHAR2=='$')THEN
                READ (NPT+1,*,end=151) jday2,qin2 !JDQIN(JJ,JB), QIN(JJ,JB)
                    ELSE
                READ (NPT+1,'(10F8.0/8X,9F8.0)',end=151) jday2,qin2 !JDQIN(JJ,JB), QIN(JJ,JB)
                ENDIF
                if(jday2<jday)then
                 jday1=jday2
                 qin1=qin2
                endif
                enddo
              !  if(jday == jdqin(jj,jb).or.jj==1)then
                if(jday == jday2.or.jj==1)then
              !  qq=qin(jj,jb)
                qq=qin2
                else
              !  qq=qin(jj,jb)+((qin(jj,jb)-qin(jj-1,jb))/(jdqin(jj,jb)-jdqin(jj-1,jb)))*(jdqin(jj,jb)-jday)
                 qq=qin2+(qin2-qin1)/(jday2-jday1)*(jday2-jday)
                endif
           !     jday1=jday2
           !     qin1=qin2
                if(j>1)cinload(INCN(1:NACIN(JB),JB),jb)=cinload(INCN(1:NACIN(JB),JB),jb)+qq*(CIN(INCN(1:NACIN(JB),JB))+CIN2(INCN(1:NACIN(JB),JB)))*0.5*(JDAY-TDAY2)*86400./1000.  ! units of kg
                ! end of inflow loading

151           CIN2(INCN(1:NACIN(JB),JB)) = CIN(INCN(1:NACIN(JB),JB))
              TDAY2                      = JDAY
              IF (J == 1 .AND. JDAY > TMSTRT) THEN
                CALL ERRORS
                WRITE (ERR,FMTF) 'Starting date [',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//CINFN(JB)
              ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
                if(jday /= 0.0)then
                  CALL ERRORS
                WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//CINFN(JB)
                else
                CALL ERRORS
                WRITE (ERR,FMTF) '[This may be a result of blank lines at the end of the file]: Julian date ',JDAY,' <= previous date ',JDAYO,' in '//CINFN(JB)
                endif
              END IF
              JDAYO = JDAY
              J     = J+1
           END DO
 150        CONTINUE
            IF (JDAY < TMEND) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//CINFN(JB)
            END IF
            CLOSE (NPT)
            close(npt+1)
            cinload(INCN(1:NACIN(JB),JB),jb)=cinload(INCN(1:NACIN(JB),JB),jb)/(jday-tstart)
            CINAVG(INCN(1:NACIN(JB),JB),JB) = CINAVG(INCN(1:NACIN(JB),JB),JB)/(JDAY-TSTART)
          ELSE
            DO J=1,IMWB   ! CHECK FOR MULTIPLE WATER BODIES AND INPUT FILE IS NOT READY YET
              IF(MWB2(J) == 'BR')THEN
                  IF(INUMMWB(J)==JB)THEN
                      IOK=1
                      EXIT
                  ENDIF
              ENDIF
          ENDDO
          IF(IOK==0)THEN  
          CALL ERRORS
            WRITE (ERR,FMTI) 'Could not open branch ',JB,' inflow concentration file '//CINFN(JB)
          ELSE
             CALL WARNINGS 
              WRITE (WRN,FMTI) 'Multiple Waterbodies: Note the following CONCENTRATION inflow file for branch ',JB,' is not yet available: '//CINFN(JB)
           ENDIF
          IOK=0
          END IF
        END IF
      END IF

!**** Distributed tributaries

      IF (DIST_TRIB(JB)) THEN

!****** Inflows

        WRITE (WIN,*) '  distributed tributary inflow'
        OPEN (UNIT=NPT,FILE=QDTFN(JB),STATUS='OLD',IOSTAT=IERR)
        IF (IERR == 0) THEN
            AID=QDTFN(JB)
        WRITE(WIN,*)'       filename:', adjustl(trim(qdtfn(jb)))
       !   READ (NPT,'(//)')
          J     = 1
          QDTO  = 0.0
          JDAYO = 0.0
        READ(NPT,'(A1)')ICHAR1
        IF(ICHAR1=='$')WRITE(WIN,*) '      QDTR file in csv format'
        READ (NPT,'(/)',END=401)
        
        ! CHECK IF ERROR -CSV FORMAT CHECK
        IF(ICHAR1 /= '$')THEN
        READ(NPT,'(A72)')TITLE(1)
            DO JJ=1,50    ! CHECK ONLY FIRST 50 CHARACTERS
                IF(TITLE(1)(JJ:JJ) == ",")THEN
                    CALL ERRORS
                    WRITE (ERR,'(A175)') 'QDT file: File not designated as a CSV file (no "$" as first character) but commas found between input values for '//QDTFN(JB)
                    EXIT
                ENDIF
            ENDDO
        BACKSPACE(NPT)
        ENDIF               
        
          DO WHILE (.TRUE.)
              IF(ICHAR1=='$')THEN
               READ (NPT,*,END=160,ERR=400) JDQDT(J,JB), QDT(J,JB)
               ELSE
               READ (NPT,'(10F8.0/8X,9F8.0)',END=160,ERR=400) JDQDT(J,JB), QDT(J,JB)
               ENDIF
               
        IF(J>1)THEN
        IF(IEEE_IS_NAN(JDQDT(J,JB)))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'JDAY for Distributed Tributary is NaN on JDAY after ',JDQDT(J-1,JB),' in '//QDTFN(JB)
        ENDIF
        ELSE
            IF(IEEE_IS_NAN(JDQDT(J,JB)))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'JDAY for Distributed Tributary is NaN on first line of data, data line:',J,' in '//QDTFN(JB)
        ENDIF
        ENDIF
        
        IF(IEEE_IS_NAN(QDT(J,JB)))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'QDT for Distributed Tributary is NaN on Julian date ',JDQDT(J,JB),' in '//QDTFN(JB)
        ENDIF
               
               
               
               
            IF (JDAYO >= TMSTRT .AND. JDQDT(J,JB) <= TMEND) THEN
              QDTMX(JB) = MAX(QDT(J,JB),QDTMX(JB))
              IF (JDQDT(J,JB) > TMEND) THEN
                QDTS(JB) = QDTS(JB)+QDTO*(TMEND-JDAYO)
              ELSE
                QDTS(JB) = QDTS(JB)+QDTO*(JDQDT(J,JB)-JDAYO)
              END IF
            END IF
            IF (J == 1 .AND. JDQDT(J,JB) > TMSTRT) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Starting date [',JDQDT(J,JB),'] > simulation start date [TMSTRT=',TMSTRT,'] in '//QDTFN(JB)
            ELSE IF (JDQDT(J,JB) <= JDAYO .AND. J /= 1) THEN
              if(jdqdt(j,jb) /= 0.0)then
                CALL ERRORS
              WRITE (ERR,FMTF) 'Julian date ',JDQDT(J,JB),' <= previous date ',JDAYO,' in '//QDTFN(JB)
              else
                  CALL ERRORS
              WRITE (ERR,FMTF) '[This may be a result of blank lines at the end of the file]: Julian date ',JDQDT(J,JB),' <= previous date ',JDAYO,' in '//QDTFN(JB)
              endif
            END IF
            JDAYO = JDQDT(J,JB)
            QDTO  = QDT(J,JB)
            J     = J+1
          END DO
 160      CONTINUE
          IF (JDQDT(J-1,JB) < TMEND) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Ending time ',JDQDT(J,JB),' < ending simulation time [TMEND=',TMEND,'] in '//QDTFN(JB)
          END IF
          CLOSE (NPT)
        ELSE
          CALL ERRORS
          WRITE (ERR,FMTI) 'Could not open branch ',JB,' distributed tributary inflow file '//QDTFN(JB)
        END IF

!****** Inflow temperatures

        WRITE (WIN,*) '  distributed tributary inflow temperature'
        OPEN (UNIT=NPT,FILE=TDTFN(JB),STATUS='OLD',IOSTAT=IERR)
        IF (IERR == 0) THEN
            AID=TDTFN(JB)
      !    READ (NPT,'(//)')
          J     = 1
          JDAYO = 0.0
          READ(NPT,'(A1)')ICHAR1
          IF(ICHAR1=='$')WRITE(WIN,*) '      QDTR file in csv format'
          READ (NPT,'(/)',END=401)
          
        ! CHECK IF ERROR -CSV FORMAT CHECK
        IF(ICHAR1 /= '$')THEN
        READ(NPT,'(A72)')TITLE(1)
            DO JJ=1,50    ! CHECK ONLY FIRST 50 CHARACTERS
                IF(TITLE(1)(JJ:JJ) == ",")THEN
                    CALL ERRORS
                    WRITE (ERR,'(A175)') 'TDT file: File not designated as a CSV file (no "$" as first character) but commas found between input values for '//TDTFN(JB)
                    EXIT
                ENDIF
            ENDDO
        BACKSPACE(NPT)
        ENDIF        
          
          DO WHILE (.TRUE.)
              IF(ICHAR1=='$')THEN
                READ (NPT,*,END=170,ERR=400) JDAY, TDT
                ELSE
                READ (NPT,'(10F8.0/8X,9F8.0)',END=170,ERR=400) JDAY, TDT
                ENDIF

        IF(IEEE_IS_NAN(JDAY))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'JDAY for Distrbuted Tributary is NaN on line of data ',J,' in '//TDTFN(JB)
        ENDIF
        IF(IEEE_IS_NAN(TDT))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'TDT for Distributed Tributary is NaN on Julian date ',JDAY,' in '//TDTFN(JB)
        ENDIF
                
                
                
                
                IF (J == 1 .AND. JDAY > TMSTRT) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Starting date [',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//TDTFN(JB)
            ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
              if(jday /= 0.0)then
                CALL ERRORS
              WRITE (ERR,FMTF) 'Julian date ',JDAY,' < previous date =',JDAYO,'] in '//TDTFN(JB)
              else
              CALL ERRORS
              WRITE (ERR,FMTF) '[This may be a result of blank lines at end of file]: Julian date ',JDAY,' < previous date=',JDAYO,'] in '//TDTFN(JB)
              end if
            END IF
            J     = J+1
            JDAYO = JDAY
          END DO
 170      CONTINUE
          IF (JDAY < TMEND) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//TDTFN(JB)
          END IF
          CLOSE (NPT)
        ELSE
          CALL ERRORS
          WRITE (ERR,FMTI) 'Could not open branch ',JB,' distributed tributary inflow temperature file '//TDTFN(JB)    !TC 09/30/03
        END IF

!****** Inflow concentrations

        IF (CONSTITUENTS .AND. NACDT(JB) > 0) THEN
        cin=0.0
        cin2=0.0
          WRITE (WIN,*) '  distributed tributary concentrations'
          OPEN (UNIT=NPT,FILE=CDTFN(JB),IOSTAT=IERR,STATUS='OLD')
          IF (IERR == 0) THEN
              AID=CDTFN(JB)
              WRITE(WIN,*)'       filename:', adjustl(trim(cdtfn(jb)))
              READ(NPT,'(A1)')ICHAR1
              IF(ICHAR1=='$')THEN
              WRITE(WIN,*) '      cdt file in csv format'
              ENDIF
              READ (NPT,'(/)',END=401)
              
        ! CHECK IF ERROR -CSV FORMAT CHECK
        IF(ICHAR1 /= '$')THEN
        READ(NPT,'(A72)')TITLE(1)
            DO JJ=1,50    ! CHECK ONLY FIRST 50 CHARACTERS
                IF(TITLE(1)(JJ:JJ) == ",")THEN
                    CALL ERRORS
                    WRITE (ERR,'(A175)') 'CDT file: File not designated as a CSV file (no "$" as first character) but commas found between input values for '//CDTFN(JB)
                    EXIT
                ENDIF
            ENDDO
        BACKSPACE(NPT)
        ENDIF       
              
            !READ (NPT,'(//)')
            J     = 1
            JDAYO = 0.0
            DO WHILE (.TRUE.)
            IF(ICHAR1=='$')THEN
              READ (NPT,*,END=180,ERR=400) JDAY, (CIN(incdt(jc,jb)),JC=1,NACDT(JB))
                    ELSE
              READ (NPT,'(F8.0,100(F8.0))',END=180,ERR=400) JDAY, (CIN(incdt(jc,jb)),JC=1,NACDT(JB))
                    ENDIF
                    
        IF(IEEE_IS_NAN(JDAY))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'JDAY for Distributed Tributary is NaN on line of data  ',J,' in '//CDTFN(JB)
         ENDIF
        DO JC=1,NACDT(JB)
         IF(IEEE_IS_NAN(CIN(incdt(jc,jb))))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'CIN for Distributed Tributary is NaN on Julian date ',JDAY,' in '//CDTFN(JB)
        ENDIF
        ENDDO             
                    
                    
                    IF (J == 1) TSTART = JDAY
              CDTMIN(INCDT(1:NACDT(JB),JB),JB) = MIN(CIN(INCDT(1:NACDT(JB),JB)),CDTMIN(INCDT(1:NACDT(JB),JB),JB))
              CDTMAX(INCDT(1:NACDT(JB),JB),JB) = MAX(CIN(INCDT(1:NACDT(JB),JB)),CDTMAX(INCDT(1:NACDT(JB),JB),JB))
              IF (J > 1) CDTAVG(INCDT(1:NACDT(JB),JB),JB) =  CDTAVG(INCDT(1:NACDT(JB),JB),JB)+(CIN(INCDT(1:NACDT(JB),JB))             &
                                                            +CIN2(INCDT(1:NACDT(JB),JB)))*0.5*(JDAY-TDAY2)
              !CIN2(INCDT(1:NACDT(JB),JB)) = CIN(INCDT(1:NACDT(JB),JB))
              !TDAY2                      = JDAY
              
              IF (J == 1 .AND. JDAY > TMSTRT) THEN
                CALL ERRORS
                WRITE (ERR,FMTF) 'Starting date [',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//CDTFN(JB)
              ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
                 if(jday /= 0.0)then
                  CALL ERRORS
                WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//CDTFN(JB)
                 else
                CALL ERRORS
                WRITE (ERR,FMTF) '[This may be a result of blank lines at end of file]: Julian date ',JDAY,' <= previous date ',JDAYO,' in '//CDTFN(JB)
                 endif
              END IF

                 ! compute mass loading from dist trib inflows

             if(j==1)then
             open(npt+1,FILE=QDTFN(JB),STATUS='OLD')
             AID=QDTFN(JB)
             UNIT=NPT+1
       !      READ  (NPT+1,'(A)') HEADER
             QINO    = 0.0
             jj=1

             READ(NPT+1,'(A1)')ICHAR2
             READ (NPT+1,'(/)')
             IF(ICHAR2=='$')THEN
             READ (NPT+1,*,ERR=402)  JDQDT(J,JB), QDT(J,JB)
                ELSE
             READ (NPT+1,'(10F8.0/8X,9F8.0)',ERR=402)  JDQDT(J,JB), QDT(J,JB)
             ENDIF
             endif
             AID=QDTFN(JB)
                do while(JDQDT(jj,JB) < jday)
                jj=jj+1
                IF(ICHAR2=='$')THEN
                READ (NPT+1,*,end=180,ERR=402) JDQDT(JJ,JB), QDT(JJ,JB)
                    ELSE
                READ (NPT+1,'(10F8.0/8X,9F8.0)',end=180,ERR=402) JDQDT(JJ,JB), QDT(JJ,JB)
                ENDIF
                enddo
                if(jday == jdqdt(jj,jb).or.jj==1)then
                qq=qdt(jj,jb)
                else
                qq=qdt(jj,jb)+((qdt(jj,jb)-qdt(jj-1,jb))/(jdqdt(jj,jb)-jdqdt(jj-1,jb)))*(jdqdt(jj,jb)-jday)
                endif
                if(qq<0.0)qq=0.0
                if(j>1)cdtload(incdt(1:NACdt(Jb),Jb),jb)=cdtload(incdt(1:NACdt(Jb),Jb),jb)+qq*(Cin(incdt(1:NACdt(Jb),Jb))+CIN2(incdt(1:NACdt(Jb),Jb)))*0.5*(JDAY-TDAY2)*86400./1000.  ! units of kg
                ! end of inflow loading

              JDAYO = JDAY
              CIN2(incdt(1:NACdt(Jb),Jb)) = CIN(incdt(1:NACdt(Jb),Jb))
              TDAY2=JDAY
              J     = J+1
            END DO
 180        CONTINUE
            IF (JDAY < TMEND) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//CDTFN(JB)
            END IF
            CLOSE (NPT)
            close(npt+1)
            cdtload(incdt(1:NACdt(Jb),Jb),jb)=cdtload(incdt(1:NACdt(Jb),Jb),jb)/(jday-tstart)                        ! CB 5/10/10 Change units to kg/day
            CDTAVG(INCDT(1:NACDT(JB),JB),JB) = CDTAVG(INCDT(1:NACDT(JB),JB),JB)/(JDAY-TSTART)
          ELSE
            CALL WARNINGS
            WRITE (ERR,FMTI) 'Could not open branch ',JB,' distributed tributary inflow concentration file '//CDTFN(JB)
          END IF
        END IF
      END IF

!**** Precipitation

      IF (PRECIPITATION(JW)) THEN

!****** Inflows

        WRITE (WIN,*) '  precipitation'
        OPEN (UNIT=NPT,FILE=PREFN(JB),STATUS='OLD')
        IF (IERR == 0) THEN
            AID=PREFN(JB)
          !READ (NPT,'(//)')
          
          J     = 1
          JDAYO = 0.0
          READ(NPT,'(A1)')ICHAR1
          IF(ICHAR1=='$')WRITE(WIN,*) '      PREQ file in csv format'
          READ (NPT,'(/)',END=401)
          
        ! CHECK IF ERROR -CSV FORMAT CHECK
        IF(ICHAR1 /= '$')THEN
        READ(NPT,'(A72)')TITLE(1)
            DO JJ=1,50    ! CHECK ONLY FIRST 50 CHARACTERS
                IF(TITLE(1)(JJ:JJ) == ",")THEN
                    CALL ERRORS
                    WRITE (ERR,'(A175)') 'PRECIPQ file: File not designated as a CSV file (no "$" as first character) but commas found between input values for '//PREFN(JB)
                    EXIT
                ENDIF
            ENDDO
        BACKSPACE(NPT)
        ENDIF       
          
          
          DO WHILE (.TRUE.)
              
            IF(ICHAR1=='$')THEN
            READ (NPT,*,END=190,ERR=400) JDQPR(J,JB), QPR(J,JB)
             ELSE
            READ (NPT,'(10F8.0/8X,9F8.0)',END=190,ERR=400) JDQPR(J,JB), QPR(J,JB)
             ENDIF
             
        IF(J>1)THEN
        IF(IEEE_IS_NAN(JDQPR(J,JB)))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'JDAY for Precipitation is NaN on JDAY after ',JDQPR(J-1,JB),' in '//PREFN(JB)
        ENDIF
        ELSE
            IF(IEEE_IS_NAN(JDQPR(J,JB)))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'JDAY for Precipitation is NaN on first line of data, data line:',J,' in '//PREFN(JB)
        ENDIF
        ENDIF
        
        IF(IEEE_IS_NAN(QPR(J,JB)))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'QPR for Precipitation is NaN on Julian date ',JDQPR(J,JB),' in '//PREFN(JB)
        ENDIF
              
            !READ (NPT,'(10F8.0/8X,9F8.0)',END=190,ERR=400) JDQPR(J,JB), QPR(J,JB)
            IF (JDAYO >= TMSTRT .AND. JDQPR(J,JB) <= TMEND) THEN
              QPRMX(JB) = MAX(QPR(J,JB),QPRMX(JB))
              IF (JDQPR(J,JB) > TMEND) THEN
                QPRS(JB) = QPRS(JB)+QPR(J,JB)*(TMEND-JDAYO)
              ELSE
              QPRS(JB) = QPRS(JB)+QPR(J,JB)*(JDQPR(J,JB)-JDAYO)
              END IF
            END IF
            IF (J == 1 .AND. JDQPR(J,JB) > TMSTRT) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Starting date [',JDQPR(J,JB),'] > simulation start date [TMSTRT=',TMSTRT,'] in '//PREFN(JB)
            ELSE IF (JDQPR(J,JB) <= JDAYO .AND. J /= 1) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Julian date ',JDQPR(J,JB),' <= previous date ',JDAYO,' in '//PREFN(JB)
            END IF
            JDAYO = JDQPR(J,JB)
            J     = J+1
          END DO
 190      CONTINUE
          IF (JDQPR(J-1,JB) < TMEND) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Ending time ',JDQPR(J,JB),' < ending simulation time [TMEND=',TMEND,'] in '//PREFN(JB)
          END IF
        ELSE
          CALL ERRORS
          WRITE (ERR,FMTI) 'Could not open branch ',JB,' precipitation file '//PREFN(JB)
        END IF
        CLOSE (NPT)

!****** Inflow temperatures

        WRITE (WIN,*) '  precipitation temperature'
        OPEN (UNIT=NPT,FILE=TPRFN(JB),STATUS='OLD',IOSTAT=IERR)
        IF (IERR == 0) THEN
            AID=TPRFN(JB)
          !READ (NPT,'(//)')
          READ(NPT,'(A1)')ICHAR1
          IF(ICHAR1=='$')WRITE(WIN,*) '      PRET file in csv format'
          READ (NPT,'(/)',END=401)
            
        ! CHECK IF ERROR -CSV FORMAT CHECK
        IF(ICHAR1 /= '$')THEN
        READ(NPT,'(A72)')TITLE(1)
            DO JJ=1,50    ! CHECK ONLY FIRST 50 CHARACTERS
                IF(TITLE(1)(JJ:JJ) == ",")THEN
                    CALL ERRORS
                    WRITE (ERR,'(A175)') 'PRECIPT file: File not designated as a CSV file (no "$" as first character) but commas found between input values for '//TPRFN(JB)
                    EXIT
                ENDIF
            ENDDO
        BACKSPACE(NPT)
        ENDIF       
          
          
          J     = 1
          JDAYO = 0.0
          DO WHILE (.TRUE.)
            if(ICHAR1=='$')then
            READ (NPT,*,END=200,ERR=400) JDAY, TPR
                  else
            READ (NPT,'(10F8.0/8X,9F8.0)',END=200,ERR=400) JDAY, TPR
            endif
                  
         IF(IEEE_IS_NAN(JDAY))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'JDAY for Precipitation is NaN on line of data  ',J,' in '//TPRFN(JB)
        ENDIF
        IF(IEEE_IS_NAN(TPR))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'TPR for Precipitation is NaN on Julian date ',JDAY,' in '//TPRFN(JB)
        ENDIF      
                  
                  
            IF (J == 1 .AND. JDAY > TMSTRT) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Starting date [',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//TPRFN(JB)
            ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Julian date',JDAY,' <= previous date ',JDAYO,' in '//TPRFN(JB)
            END IF
            J     = J+1
            JDAYO = JDAY
          END DO
 200      CONTINUE
          IF (JDAY < TMEND) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//TPRFN(JB)
          END IF
          CLOSE (NPT)
        ELSE
          CALL ERRORS
          WRITE (ERR,FMTI) 'Could not open branch ',JB,' precipitation temperature file '//TPRFN(JB)
        END IF

!****** Inflow concentrations

        IF (CONSTITUENTS .AND. NACPR(JW) > 0) THEN
          WRITE (WIN,*) '  precipitation concentrations'
          OPEN (UNIT=NPT,FILE=CPRFN(JB),IOSTAT=IERR,STATUS='OLD')
          IF (IERR == 0) THEN
              AID=CPRFN(JB)
              READ(NPT,'(A1)')ICHAR1
              IF(ICHAR1=='$')THEN
              WRITE(WIN,*) '  cpr file in csv format'
              ENDIF
              READ (NPT,'(/)',END=401)
              
         ! CHECK IF ERROR -CSV FORMAT CHECK
        IF(ICHAR1 /= '$')THEN
        READ(NPT,'(A72)')TITLE(1)
            DO JJ=1,50    ! CHECK ONLY FIRST 50 CHARACTERS
                IF(TITLE(1)(JJ:JJ) == ",")THEN
                    CALL ERRORS
                    WRITE (ERR,'(A175)') 'PRECIPC file: File not designated as a CSV file (no "$" as first character) but commas found between input values for '//CPRFN(JB)
                    EXIT
                ENDIF
            ENDDO
        BACKSPACE(NPT)
        ENDIF       
              
              
            !READ (NPT,'(//)')
            J     = 1
            JDAYO = 0.0
            DO WHILE (.TRUE.)
                IF(CHAR1=='$')THEN
                READ (NPT,*,END=210,ERR=400) JDAY, (CPR(JC),JC=1,NACPR(JW))
                    ELSE
              READ (NPT,'(F8.0,1000(F8.0))',END=210,ERR=400) JDAY, (CPR(JC),JC=1,NACPR(JW))
                    ENDIF
                    
         IF(IEEE_IS_NAN(JDAY))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'JDAY for Precipitation is NaN on LINE of data ',J,' in '//CPRFN(JB)
         ENDIF
        DO JC=1,NACPR(JW)
         IF(IEEE_IS_NAN(CPR(JC)))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'CPR for Precipitation is NaN on Julian date ',JDAY,' in '//CPRFN(JB)
        ENDIF
        ENDDO       
                    
                    
              IF (J == 1 .AND. JDAY > TMSTRT) THEN
                CALL ERRORS
                WRITE (ERR,FMTF) 'Starting date [',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//CPRFN(JB)
              ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
                CALL ERRORS
                WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//CPRFN(JB)
              END IF
              JDAYO = JDAY
              J     = J+1
            END DO
 210        CONTINUE
            IF (JDAY < TMEND) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//CPRFN(JB)
            END IF
            CLOSE (NPT)
          ELSE
            CALL WARNINGS
          WRITE (ERR,FMTI) 'Could not open branch ',JB,' precipitation concentration file '//CPRFN(JB)
          END IF
        END IF
      END IF

!**** Upstream head

      IF (UH_EXTERNAL(JB)) THEN

!****** Elevations

        IFLAGEU = 0
        WRITE (WIN,*) '  upstream head elevations'
        OPEN (UNIT=NPT,FILE=EUHFN(JB),STATUS='OLD',IOSTAT=IERR)
        IF (IERR == 0) THEN
            AID=EUHFN(JB)
            READ( NPT,'(A1)')INFORMAT1
             IF(INFORMAT1=='$')THEN
              EUHF=1
             ENDIF
          !READ (NPT,'(//)')
          READ(NPT,'(/)')
          J     = 1
          JDAYO = 0.0
          DO WHILE (.TRUE.)
            !READ (NPT,'(10F8.0/8X,9F8.0)',END=220,ERR=400) JDAY, EUH
              IF(EUHF>0)THEN
                READ (NPT,*,END=220,ERR=400)  JDAY,EUH
                ELSE
                READ (NPT,'(10F8.0/8X,9F8.0)',END=220,ERR=400)JDAY,EUH
                ENDIF
            IF (JDAY >= TMSTRT .AND. IFLAGEU == 0) THEN
              ELDIFF = ABS(EUH-ELWS(US(JB)))
              IF (ELDIFF > 0.1) THEN
                CALL WARNINGS
                WRITE (WRN,FMTI) 'Boundary elevation and downstream water surface elevation difference > 0.1m for branch ',JB
              END IF
              IFLAGEU = 1
            END IF
            IF (J == 1 .AND. JDAY > TMSTRT) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Starting date [',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//EUHFN(JB)
            ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//EUHFN(JB)
            END IF
            IF(J>1)THEN
                !IF(ABS(EUHO-EUH)/(JDAY-JDAYO) > 20.)THEN
                !    CALL ERRORS
                !    WRITE (ERR,FMTF) 'Upstream Head BC: Water level change > 20 m/day on Julian date ',JDAY,' compared to previous date ',JDAYO,' in '//EDHFN(JB)
                IF(ABS(EDHO-EDH)/(JDAY-JDAYO) > 15.)THEN
                    CALL WARNINGS
                    WRITE (WRN,FMTF) 'Upstream Head BC: Water level change > 15 m/day on Julian date ',JDAY,' compared to previous date ',JDAYO,' in '//EDHFN(JB) 
                ENDIF
            ENDIF
            
            EUHO=EUH
            JDAYO = JDAY
            J     = J+1
          END DO
 220      CONTINUE
          IF (JDAY < TMEND) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//EUHFN(JB)
          END IF
          CLOSE (NPT)
        ELSE
          CALL ERRORS
          WRITE (ERR,FMTI) 'Could not open branch ',JB,' external upstream head file '//EUHFN(JB)
        END IF

!****** Temperatures

        WRITE (WIN,*) '  upstream head temperatures'
        OPEN (UNIT=NPT,FILE=TUHFN(JB),STATUS='OLD',IOSTAT=IERR)
        IF (IERR == 0) THEN
            AID=TUHFN(JB)
            READ( NPT,'(A2)')INFORMAT2
            IF(INFORMAT2=='$T')THEN
              TUHF=1
            ELSEIF(INFORMAT2(1:1)=='$')THEN
              TUHF=2
            ENDIF
          !READ (NPT,'(//)')
            READ (NPT,'(/)')
          J     = J+1
          JDAYO = 0.0
          TUHMAX(JB)=0.0
          TUHMIN(JB)=100.
          DO WHILE (.TRUE.)
            !READ (NPT,'(10F8.0/8X,9F8.0)',END=230,ERR=400) JDAY, (TUH,K=2,KB(US(JB)))
          IF(TUHF==1)THEN
          READ (NPT,*,END=230,ERR=400)  JDAY,XX
          TUH(2:KB(US(JB)))=XX
          ELSEIF(TUHF==2)THEN
          READ (NPT,*,END=230,ERR=400)  JDAY, (TUH(K),K=2,KB(US(JB)))
          ELSE
          READ (NPT,'(10F8.0:/(8X,9F8.0))',END=230,ERR=400) JDAY, (TUH(K),K=2,KB(US(JB)))
          ENDIF      
          
          TUHMAX(JB)=MAX(TUH(2),TUHMAX(JB))
          TUHMIN(JB)=MIN(TUH(2),TUHMIN(JB))
          
            IF (J == 1 .AND. JDAY > TMSTRT) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Starting date [',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//TUHFN(JB)
            ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//TUHFN(JB)
            END IF
            DO K=2,KB(US(JB))-1
                IF(TUH(K) < TUH(K+1))THEN
                    CALL WARNINGS
                    WRITE (WRN,'(A,F8.3,A,I3,A,A)') 'Upstream head temperature: On Julian date ',JDAY,' temperature inversion at K=',K,' in '//TUHFN(JB)
                ENDIF
            ENDDO
            J     = J+1
            JDAYO = JDAY
          END DO
 230      CONTINUE
          IF (JDAY < TMEND) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//TUHFN(JB)
         END IF
          CLOSE (NPT)
        ELSE
          CALL ERRORS
          WRITE (ERR,FMTI) 'Could not open branch ',JB,' external upstream head temperature file '//TUHFN(JB)
        END IF

!****** Constituent concentrations

        IF (CONSTITUENTS .AND. NAC > 0) THEN
          WRITE (WIN,*) '  upstream head concentrations'
          OPEN (UNIT=NPT,FILE=CUHFN(JB),STATUS='OLD',IOSTAT=IERR)
          IF (IERR == 0) THEN
              AID=CUHFN(JB)
 
          READ( NPT,'(A2)')INFORMAT2
          IF(INFORMAT2=='$T')THEN
              CUHF=1
          ELSEIF(INFORMAT2(1:1)=='$')THEN
              CUHF=2
          ENDIF
          READ (NPT,'(/)')
          
            J     = 1
            JDAYO = 0.0
            CUHMIN(:,JB)=100000.0
            CUHMAX(:,JB)=0.0
            CUHAVG(:,JB)=0.0
            DO WHILE (.TRUE.)
                IF(CUHF==1)THEN
                READ (NPT,*,END=240,ERR=400)  JDAY,(CUH(2,JAC),JAC=1,NAC)
                    !DO JAC=1,NAC
                    !CUH=XX(CN(JAC)) 
                    !ENDDO
                ELSE
                    DO JAC=1,NAC          
                    IF(CUHF==2)THEN
                    READ (NPT,*,END=240,ERR=400)  JDAY,(CUH(K,JAC),K=2,KB(US(JB)))
                    ELSE
                    READ (NPT,'(10F8.0:/(8X,9F8.0))',END=240,ERR=400) JDAY,(CUH(K,JAC),K=2,KB(US(JB)))
                    ENDIF
                    END DO
                END IF
          
            DO JAC=1,NAC
                CUHMIN(JAC,JB)=MIN(CUH(2,JAC),CUHMIN(JAC,JB))
                CUHMAX(JAC,JB)=MAX(CUH(2,JAC),CUHMAX(JAC,JB))
                CUHAVG(JAC,JB)=CUHAVG(JAC,JB)+CUH(2,JAC)
            ENDDO
            
                
            !READ (NPT,'(//)')
            !J     = 1
            !JDAYO = 0.0
            !DO WHILE (.TRUE.)
            !  DO JC=1,NAC
            !    READ (NPT,*,END=240,ERR=400) JDAY, (CUH,K=2,KB(US(JB)))
                IF (J == 1 .AND. JDAY > TMSTRT) THEN
                  CALL ERRORS
                  WRITE (ERR,FMTF) 'Starting date [',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//CUHFN(JB)
                ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
                  CALL ERRORS
                  WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//CUHFN(JB)
                END IF
              !END DO
              JDAYO = JDAY
              J     = J+1
            END DO
240         CONTINUE
            DO JAC=1,NAC
            CUHAVG(JAC,JB)=CUHAVG(JAC,JB)/REAL(J-1)
            ENDDO
            IF (JDAY < TMEND) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//CUHFN(JB)
           END IF
            CLOSE (NPT)
          ELSE
            CALL ERRORS
            WRITE (ERR,FMTI) 'Could not open branch ',JB,' external upstream head concentration file '//CUHFN(JB)
          END IF
        END IF
      END IF

!**** Downstream head

      IF (DH_EXTERNAL(JB)) THEN

!****** Elevations

        IFLAGED = 0
        WRITE (WIN,*) '  downstream head elevations'
        OPEN (UNIT=NPT,FILE=EDHFN(JB),STATUS='OLD',IOSTAT=IERR)
        IF (IERR == 0) THEN
            AID=EDHFN(JB)
             READ( NPT,'(A1)')INFORMAT1
             IF(INFORMAT1=='$')THEN
              EDHF=1
             ENDIF
          !READ (NPT,'(//)')
          READ(NPT,'(/)')
         J     = 1
          JDAYO = 0.0
          DO WHILE (.TRUE.)
           ! READ (NPT,'(10F8.0/8X,9F8.0)',END=250,ERR=400) JDAY, EDH
                IF(EDHF>0)THEN
                READ (NPT,*,END=250,ERR=400)  JDAY,EDH
                ELSE
                READ (NPT,'(10F8.0/8X,9F8.0)',END=250,ERR=400)JDAY,EDH
                ENDIF
            IF (JDAY >= TMSTRT .AND. IFLAGED == 0) THEN
              ELDIFF = ABS(EDH-ELWS(DS(JB)))
              IF (ELDIFF > 0.1) THEN
                CALL WARNINGS
                WRITE (WRN,FMTI) 'Boundary elevation and downstream water surface elevation difference > 0.1m for branch ',JB
              END IF
              IFLAGED = 1
            END IF
            IF (J == 1 .AND. JDAY > TMSTRT) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Starting date [',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//EDHFN(JB)
            ELSE IF (JDAY < JDAYO .AND. J /= 1) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//EDHFN(JB)
            END IF
            IF(J>1)THEN
                !IF(ABS(EDHO-EDH)/(JDAY-JDAYO) > 20.)THEN
                !    CALL ERRORS
                !    WRITE (ERR,FMTF) 'Downstream Head BC: Water level change > 20 m/day on Julian date ',JDAY,' compared to previous date ',JDAYO,' in '//EDHFN(JB)
                IF(ABS(EDHO-EDH)/(JDAY-JDAYO) > 15.)THEN
                    CALL WARNINGS
                    WRITE (WRN,FMTF) 'Downstream Head BC: Water level change > 15 m/day on Julian date ',JDAY,' compared to previous date ',JDAYO,' in '//EDHFN(JB) 
                ENDIF
            ENDIF
            
            EDHO=EDH
            JDAYO = JDAY
            J     = J+1
          END DO
 250      CONTINUE
          IF (JDAY < TMEND) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//EDHFN(JB)
          END IF
          CLOSE (NPT)
        ELSE
          CALL ERRORS
          WRITE (ERR,FMTI) 'Could not open branch ',JB,' external downstream head file '//EDHFN(JB)
        END IF

!****** Temperatures

        WRITE (WIN,*) '  downstream head temperatures'
        OPEN (UNIT=NPT,FILE=TDHFN(JB),STATUS='OLD',IOSTAT=IERR)
        IF (IERR == 0) THEN
            AID=TDHFN(JB)
            READ( NPT,'(A2)')INFORMAT2
            IF(INFORMAT2=='$T')THEN
              TDHF=1
            ELSEIF(INFORMAT2(1:1)=='$')THEN
              TDHF=2
            ENDIF
          !READ (NPT,'(//)')
            READ (NPT,'(/)')
          J     = 1
          JDAYO = 0.0
          TDHMAX(JB)=0.0
          TDHMIN(JB)=100.
          DO WHILE (.TRUE.)
            !READ (NPT,'(10F8.0:/(8X,9F8.0))',END=260,ERR=400) JDAY, (TDH,K=2,KB(DS(JB)))
          IF(TDHF==1)THEN
          READ (NPT,*,END=260,ERR=400)  JDAY,XX
          TDH(2:KB(DS(JB)))=XX
          ELSEIF(TDHF==2)THEN
          READ (NPT,*,END=260,ERR=400)  JDAY, (TDH(K),K=2,KB(DS(JB)))
          ELSE
          READ (NPT,'(10F8.0:/(8X,9F8.0))',END=260,ERR=400) JDAY, (TDH(K),K=2,KB(DS(JB)))
          ENDIF    
          
          TDHMAX(JB)=MAX(TDH(2),TDHMAX(JB))
          TDHMIN(JB)=MIN(TDH(2),TDHMIN(JB))
          
            IF (J == 1 .AND. JDAY > TMSTRT) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Starting date [',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//TDHFN(JB)
            ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//TDHFN(JB)
            END IF
            DO K=2,KB(DS(JB))-1
                IF(TDH(K) < TDH(K+1))THEN
                    CALL WARNINGS
                    WRITE (WRN,'(A,F8.3,A,I3,A,A)') 'Downstream head temperature: On Julian date ',JDAY,' temperature inversion at K=',K,' in '//TDHFN(JB)
                ENDIF
            ENDDO
            
            J     = J+1
            JDAYO = JDAY
          END DO
 260      CONTINUE
          IF (JDAY < TMEND) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//TDHFN(JB)
          END IF
          CLOSE (NPT)
        ELSE
          CALL ERRORS
          WRITE (ERR,FMTI) 'Could not open branch ',JB,' external downstream head temperature file '//TDHFN(JB)
        END IF

!****** Constituent concentrations

        IF (CONSTITUENTS .AND. NAC > 0) THEN
          WRITE (WIN,*) '  downstream head concentrations'
          OPEN (UNIT=NPT,FILE=CDHFN(JB),STATUS='OLD',IOSTAT=IERR)
          IF (IERR == 0) THEN
              AID=CDHFN(JB)
            !READ (NPT,'(A)') (AID1,I=1,3)
          READ( NPT,'(A2)')INFORMAT2
          IF(INFORMAT2=='$T')THEN
              CDHF=1
          ELSEIF(INFORMAT2(1:1)=='$')THEN
              CDHF=2
          ENDIF
          READ (NPT,'(/)')
          
            J     = 1
            JDAYO = 0.0
            CDHMIN(:,JB)=100000.0
            CDHMAX(:,JB)=0.0
            CDHAVG(:,JB)=0.0
            DO WHILE (.TRUE.)
                IF(CDHF==1)THEN
                
                READ(NPT,'(A10)',END=270)TITLE(1)      ! CHECK FOR JUST BLANK SPACES INTERPRETED AS DATA
                IF(TITLE(1)(1:10) == ',,,,,,,,,,')THEN
                    call errors
                    WRITE(ERR,'(A,I4)')'Concentration Downstream Head BC has extra blank lines that should be deleted. Branch #:',jb
                    EXIT
                ELSE
                    BACKSPACE(NPT)
                ENDIF
                
                READ (NPT,*,END=270,ERR=400)  JDAY,(CDH(2,JAC),JAC=1,NAC)
                    !DO JAC=1,NAC
                    !CDH=XX(CN(JAC)) 
                    !ENDDO
                ELSE
                    DO JAC=1,NAC          
                    IF(CDHF==2)THEN
                    READ (NPT,*,END=270,ERR=400)  JDAY,(CDH(K,JAC),K=2,KB(DS(JB)))
                    ELSE
                    READ (NPT,'(10F8.0:/(8X,9F8.0))',END=270,ERR=400) JDAY,(CDH(K,JAC),K=2,KB(DS(JB)))
                    ENDIF
                    END DO
                END IF
            DO JAC=1,NAC
                CDHMIN(JAC,JB)=MIN(CDH(2,JAC),CDHMIN(JAC,JB))
                CDHMAX(JAC,JB)=MAX(CDH(2,JAC),CDHMAX(JAC,JB))
                CDHAVG(JAC,JB)=CDHAVG(JAC,JB)+CDH(2,JAC)
            ENDDO
            !J     = 1
            !JDAYO = 0.0
            !DO WHILE (.TRUE.)
              !DO JC=1,NAC
                !READ (NPT,'(10F8.0/(:8X,9F8.0))',END=270,ERR=400) JDAY, (CDH,K=2,KB(DS(JB)))
                IF (J == 1 .AND. JDAY > TMSTRT) THEN
                  CALL ERRORS
                  WRITE (ERR,FMTF) 'Starting date [',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] OR Constituent is missing (Must have same active constituents as in control file) in '//CDHFN(JB)
                ELSE IF (JDAY < JDAYO .AND. J /= 1) THEN
                  CALL ERRORS
                  WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//CDHFN(JB)
                END IF
              !END DO
              JDAYO = JDAY
              J     = J+1
            END DO
270         CONTINUE
            DO JAC=1,NAC
            CDHAVG(JAC,JB)=CDHAVG(JAC,JB)/REAL(J-1)
            ENDDO
            IF (JDAY < TMEND) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] OR Constituent is missing (Must have same active constituents as in control file) in '//CDHFN(JB)
            END IF
            CLOSE (NPT)
          ELSE
            CALL ERRORS
            WRITE (ERR,FMTI) 'Could not open branch ',JB,' external downstream head concentration file '//CDHFN(JB)
          END IF
        END IF
      END IF
    END DO
  END DO

! Tributaries

  IF (TRIBUTARIES) THEN
    DO JT=1,NTR
      WRITE (TRA,'(I0)') JT

!**** Inflows

      WRITE (WIN,*) '  tributary inflows'
      OPEN (UNIT=NPT,FILE=QTRFN(JT),STATUS='OLD',IOSTAT=IERR)
      IF (IERR == 0) THEN
          AID=QTRFN(JT)
     !   READ  (NPT,'(A)') HEADER

        READ(NPT,'(A1)')ICHAR1
        IF(ICHAR1=='$')WRITE(WIN,*) '      QTR file in csv format'
        READ (NPT,'(/)')
        
        ! CHECK IF ERROR -CSV FORMAT CHECK
        IF(ICHAR1 /= '$')THEN
        READ(NPT,'(A72)')TITLE(1)
            DO JJ=1,50    ! CHECK ONLY FIRST 50 CHARACTERS
                IF(TITLE(1)(JJ:JJ) == ",")THEN
                    CALL ERRORS
                    WRITE (ERR,'(A175)') 'QTR file: File not designated as a CSV file (no "$" as first character) but commas found between input values for '//QTRFN(JT)
                    EXIT
                ENDIF
            ENDDO
        BACKSPACE(NPT)
        ENDIF       
        
        J       =  1
        JDAYO   =  0.0
        DO WHILE (.TRUE.)
          IF(ICHAR1=='$')THEN
          READ (NPT,*,END=280,ERR=400) JDQTR(J,JT), QTR(J,JT)
          ELSE
          READ (NPT,'(10F8.0/8X,9F8.0)',END=280,ERR=400) JDQTR(J,JT), QTR(J,JT)
          ENDIF
          
        IF(IEEE_IS_NAN(JDQTR(J,JT)))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'JDAY for Tributary is NaN on line of data ',J,' in '//QTRFN(JT)
        ENDIF
        IF(IEEE_IS_NAN(QTR(J,JT)))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'QTR for Tributary is NaN on Julian date ',JDQTR(J,JT),' in '//QTRFN(JT)
        ENDIF
          
          
          IF (JDAYO >= TMSTRT .AND. JDQTR(J,JT) <= TMEND) THEN
            IF (JBTR(JT) /= 0) QTRMX(JT,JBTR(JT)) = MAX(QTR(J,JT),QTRMX(JT,JBTR(JT)))
            IF (JDQTR(J,JT) > TMEND) THEN
              QTRS(JT,JBTR(JT)) = QTRS(JT,JBTR(JT))+QTR(J,JT)*(TMEND-JDAYO)
            ELSE
              QTRS(JT,JBTR(JT)) = QTRS(JT,JBTR(JT))+QTR(J,JT)*(JDQTR(J,JT)-JDAYO)
            END IF
          END IF
          IF (J == 1 .AND. JDQTR(J,JT) > TMSTRT) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Starting date [',JDQTR(J,JT),'] > simulation start date [TMSTRT=',TMSTRT,'] in '//QTRFN(JT)
          ELSE IF (JDQTR(J,JT) <= JDAYO .AND. J /= 1) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Julian date ',JDQTR(J,JT),' <= previous date ',JDAYO,' in '//QTRFN(JT)
          END IF
          JDAYO = JDQTR(J,JT)
          J     = J+1
        END DO
 280    CONTINUE
        IF (JDQTR(J-1,JT) < TMEND) THEN
          CALL ERRORS
          WRITE (ERR,FMTF) 'Ending time ',JDQTR(J-1,JT),' < ending simulation time [TMEND=',TMEND,'] in '//QTRFN(JT)
        END IF
        CLOSE (NPT)
      ELSE
          DO J=1,IMWB   ! CHECK FOR MULTIPLE WATER BODIES AND INPUT FILE IS NOT READY YET
              IF(MWB2(J) == 'TR')THEN
                  IF(INUMMWB(J)==JT)THEN
                      IOK=1
                      EXIT
                  ENDIF
              ENDIF
          ENDDO
          IF(IOK==0)THEN  
          CALL ERRORS
          WRITE (ERR,FMTI) 'Could not open tributary ',JT,' inflow file '//QTRFN(JT)
          ELSE
             CALL WARNINGS 
              WRITE (WRN,FMTI) 'Multiple Waterbodies: Note the following FLOW inflow file for TRIBUTARY ',JT,' is not yet available: '//QTRFN(JT)
           ENDIF
          IOK=0
      END IF

!**** Inflow temperatures

      WRITE (WIN,*) '  tributary inflow temperatures'
      OPEN (UNIT=NPT,FILE=TTRFN(JT),STATUS='OLD',IOSTAT=IERR)
      IF (IERR == 0) THEN
          AID=TTRFN(JT)
     !   READ  (NPT,'(A)') HEADER
        J       = 1
        JDAYO   = 0.0

        READ(NPT,'(A1)')ICHAR1
        IF(ICHAR1=='$')WRITE(WIN,*) '      TTR file in csv format'
        READ (NPT,'(/)')
        
        ! CHECK IF ERROR -CSV FORMAT CHECK
        IF(ICHAR1 /= '$')THEN
        READ(NPT,'(A72)')TITLE(1)
            DO JJ=1,50    ! CHECK ONLY FIRST 50 CHARACTERS
                IF(TITLE(1)(JJ:JJ) == ",")THEN
                    CALL ERRORS
                    WRITE (ERR,'(A175)') 'TTR file: File not designated as a CSV file (no "$" as first character) but commas found between input values for '//TTRFN(JT)
                    EXIT
                ENDIF
            ENDDO
        BACKSPACE(NPT)
        ENDIF       

        DO WHILE (.TRUE.)
           IF(ICHAR1=='$')THEN
           READ (NPT,*,END=290,ERR=400) JDAY, TTR
                 ELSE
           READ (NPT,'(10F8.0/8X,9F8.0)',END=290,ERR=400) JDAY, TTR
                 ENDIF
                 
            IF(IEEE_IS_NAN(JDAY))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'JDAY for Tributary is NaN on line of data ',J,' in '//TTRFN(JT)
        ENDIF
        IF(IEEE_IS_NAN(TTR))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'TTR for Tributary is NaN on Julian date ',JDAY,' in '//TTRFN(JT)
        ENDIF
                 
                 
                 
          IF (J == 1 .AND. JDAY > TMSTRT) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Starting date [',JDAY,' > simulation start date [TMSTRT=',TMSTRT,'] in '//TTRFN(JT)
          ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//TTRFN(JT)
          END IF
            TRTMAX(JT)=MAX(TTR,TRTMAX(JT))
            TRTMIN(JT)=MIN(TTR,TRTMIN(JT))
            IF(TTR > 40. .AND. TTR < 50.)THEN
                CALL WARNINGS
                WRITE(WRN,'(A,F12.1,A,F8.2,A,I4,A)')'Tributary temperature[TTR=',ttr,']>40C on JDAY=',JDAY,' for tributary:',JT,' in file:'//TTRFN(JT)
            ENDIF
            IF(TTR > 50.)THEN
                CALL ERRORS
                WRITE(ERR,'(A,F12.1,A,F8.2,A,I4,A)')'Tributary temperature[TTR=',ttr,']>50C on JDAY=',JDAY,' for tributary:',JT,' in file:'//TTRFN(JT)
            ENDIF
            IF(TTR < -5.)THEN
                CALL ERRORS
                WRITE(ERR,'(A,F12.1,A,F8.2,A,I4,A)')'Tributary temperature[TTR=',ttr,']<-5C on JDAY=',JDAY,' for tributary:',JT,' in file:'//TTRFN(JT)
            ENDIF

          JDAYO = JDAY
          J     = J+1
        END DO
 290    CONTINUE
        IF (JDAY < TMEND) THEN
          CALL ERRORS
          WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//TTRFN(JT)
        END IF
        CLOSE (NPT)
      ELSE
          DO J=1,IMWB   ! CHECK FOR MULTIPLE WATER BODIES AND INPUT FILE IS NOT READY YET
              IF(MWB2(J) == 'TR')THEN
                  IF(INUMMWB(J)==JT)THEN
                      IOK=1
                      EXIT
                  ENDIF
              ENDIF
          ENDDO
          IF(IOK==0)THEN  
          CALL ERRORS
        WRITE (ERR,FMTI) 'Could not open tributary ',JT,' inflow temperature file '//TTRFN(JT)
          ELSE
             CALL WARNINGS 
              WRITE (WRN,FMTI) 'Multiple Waterbodies: Note the following TEMPERATURE inflow file for TRIBUTARY ',JT,' is not yet available: '//TTRFN(JT)
           ENDIF
          IOK=0
      END IF

!**** Inflow concentrations
      cin=0.0
      cin2=0.0
      CTRMIN(TRCN(1:NACTR(JT),JT),JT) = 1000.0
      CTRMAX(TRCN(1:NACTR(JT),JT),JT) = 0.0
      CTRAVG(TRCN(1:NACTR(JT),JT),JT) = 0.0
      IF (CONSTITUENTS .AND. NACTR(JT) > 0) THEN
        WRITE (WIN,*) '  tributary inflow concentrations'
        OPEN (UNIT=NPT,FILE=CTRFN(JT),IOSTAT=IERR,STATUS='OLD')
        WRITE(WIN,*) '       filename:', adjustl(trim(ctrfn(jt)))
        IF (IERR == 0) THEN
            AID=CTRFN(JT)
        READ(NPT,'(A1)')ICHAR1
              IF(ICHAR1=='$')THEN
              WRITE(WIN,*) '       ctr file in csv format'
              ENDIF
              READ (NPT,'(/)')
              
        ! CHECK IF ERROR -CSV FORMAT CHECK
        IF(ICHAR1 /= '$')THEN
        READ(NPT,'(A72)')TITLE(1)
            DO JJ=1,50    ! CHECK ONLY FIRST 50 CHARACTERS
                IF(TITLE(1)(JJ:JJ) == ",")THEN
                    CALL ERRORS
                    WRITE (ERR,'(A175)') 'CTR file: File not designated as a CSV file (no "$" as first character) but commas found between input values for '//CTRFN(JT)
                    EXIT
                ENDIF
            ENDDO
        BACKSPACE(NPT)
        ENDIF       
              
              
         ! READ (NPT,'(//)')
          J     = 1
          JDAYO = 0.0
          DO WHILE (.TRUE.)
            IF(ICHAR1=='$')THEN
                                
            READ (NPT,*,END=300,ERR=400) JDAY, (CIN(TRCN(JC,JT)),JC=1,NACTR(JT))
                ELSE
            READ (NPT,'(F8.0,1000(F8.0))',END=300,ERR=400) JDAY, (CIN(TRCN(JC,JT)),JC=1,NACTR(JT))
            ENDIF
                
         IF(IEEE_IS_NAN(JDAY))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'JDAY for Tributary is NaN on line of data ',J,' in '//CTRFN(JT)
         ENDIF
        DO JC=1,NACTR(JT)
         IF(IEEE_IS_NAN(CIN(TRCN(JC,JT))))THEN
          CALL ERRORS
          WRITE (ERR,FMTF)'CTR for Tributary is NaN on Julian date ',JDAY,' in '//CTRFN(JT)
        ENDIF
        ENDDO                    
                
                
            IF (J == 1) TSTART = JDAY
            CTRMIN(TRCN(1:NACTR(JT),JT),JT) = MIN(CIN(TRCN(1:NACTR(JT),JT)),CTRMIN(TRCN(1:NACTR(JT),JT),JT))
            CTRMAX(TRCN(1:NACTR(JT),JT),JT) = MAX(CIN(TRCN(1:NACTR(JT),JT)),CTRMAX(TRCN(1:NACTR(JT),JT),JT))

    ! compute mass loading from trib inflows
    !   WRITE(WIN,*)'       compute inflow loading'
             if(j==1)then
             open(npt+1,FILE=QTRFN(JT),STATUS='OLD')
             AID=QTRFN(JT);UNIT=NPT+1
           !  READ  (NPT+1,'(A)') HEADER
             QINO    = 0.0
             jj=1
             READ(NPT+1,'(A1)')ICHAR2
             READ (NPT+1,'(/)')
             IF(ICHAR2=='$')THEN
             READ (NPT+1,*,ERR=402) JDQTR(J,JT), QTR(J,JT)
             ELSE
             READ (NPT+1,'(10F8.0/8X,9F8.0)',ERR=402) JDQTR(J,JT), QTR(J,JT)
             ENDIF
             endif
                AID=QTRFN(JT)
                do while(jdqtr(jj,jt) < jday)
                jj=jj+1
                IF(ICHAR2=='$')THEN
                READ (NPT+1,*,end=299,ERR=402) JDQtr(JJ,Jt), Qtr(JJ,Jt)
                    ELSE
                READ (NPT+1,'(10F8.0/8X,9F8.0)',end=299,ERR=402) JDQtr(JJ,Jt), Qtr(JJ,Jt)
                ENDIF
                enddo
                if(jday == jdqtr(jj,jt).or.jj==1)then
                qq=qtr(jj,jt)
                else
                qq=qtr(jj,jt)+((qtr(jj,jt)-qtr(jj-1,jt))/(jdqtr(jj,jt)-jdqtr(jj-1,jt)))*(jdqtr(jj,jt)-jday)
                endif
                if(qq<0.0)qq=0.0
                if(j>1)ctrload(trcn(1:NACtr(Jt),Jt),jt)=ctrload(trcn(1:NACtr(Jt),Jt),jt)+qq*(Cin(TRCN(1:NACTR(JT),JT))+CIN2(TRCN(1:NACTR(JT),JT)))*0.5*(JDAY-TDAY2)*86400./1000.  ! units of kg
                ! end of inflow loading

299         IF (J > 1) CTRAVG(TRCN(1:NACTR(JT),JT),JT) =  CTRAVG(TRCN(1:NACTR(JT),JT),JT)+(CIN(TRCN(1:NACTR(JT),JT))               &
                                                          +CIN2(TRCN(1:NACTR(JT),JT)))*0.5*(JDAY-TDAY2)
            CIN2(TRCN(1:NACTR(JT),JT)) = CIN(TRCN(1:NACTR(JT),JT))
            TDAY2                      = JDAY
            IF (J == 1 .AND. JDAY > TMSTRT) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Starting date [',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//CTRFN(JT)
            ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//CTRFN(JT)
            END IF
            JDAYO = JDAY
            J     = J+1
          END DO
 300      CONTINUE
          IF (JDAY < TMEND) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//CTRFN(JT)
          END IF
          CLOSE (NPT)
          close(npt+1)
          ctrload(trcn(1:NACtr(Jt),Jt),jt)=ctrload(trcn(1:NACtr(Jt),Jt),jt)/(JDAY-TSTART)                              !CB 5/11/10 convert to units of kg/day
          CTRAVG(TRCN(1:NACTR(JT),JT),JT) = CTRAVG(TRCN(1:NACTR(JT),JT),JT)/(JDAY-TSTART)                              !SW 01/07/01
        ELSE
          DO J=1,IMWB   ! CHECK FOR MULTIPLE WATER BODIES AND INPUT FILE IS NOT READY YET
              IF(MWB2(J) == 'TR')THEN
                  IF(INUMMWB(J)==JT)THEN
                      IOK=1
                      EXIT
                  ENDIF
              ENDIF
          ENDDO
          IF(IOK==0)THEN  
          CALL ERRORS
          WRITE (ERR,FMTI) 'Could not open tributary ',JT, ' inflow concentration file '//CTRFN(JT)
          ELSE
             CALL WARNINGS 
              WRITE (WRN,FMTI) 'Multiple Waterbodies: Note the following CONCENTRATION inflow file for TRIBUTARY ',JT,' is not yet available: '//CTRFN(JT)
           ENDIF
          IOK=0 
        END IF
      END IF
    END DO
  END IF

! Outflows

  WRITE (WIN,*) '  outflows'
  DO JW=1,NWB
    DO JB=BS(JW),BE(JW)
      WRITE (BRA,'(I0)') JW
      IF (DQ_EXTERNAL(JB)) THEN
        IF (NSTR(JB) /= 0) THEN
          OPEN (UNIT=NPT,FILE=QOTFN(JB),STATUS='OLD',IOSTAT=IERR)
          IF (IERR == 0) THEN
              AID=QOTFN(JB)
!            READ  (NPT,'(A)') HEADER
             READ(NPT,'(A1)')ICHAR1
             WRITE(WIN,*) '     for branch:',JB,'  file:', adjustl(trim(qotfn(jb)))
             IF(ICHAR1=='$')WRITE(WIN,*)'     csv format'
             READ(NPT,'(/)')
             
        ! CHECK IF ERROR -CSV FORMAT CHECK
        IF(ICHAR1 /= '$')THEN
        READ(NPT,'(A72)')TITLE(1)
            DO JJ=1,50    ! CHECK ONLY FIRST 50 CHARACTERS
                IF(TITLE(1)(JJ:JJ) == ",")THEN
                    CALL ERRORS
                    WRITE (ERR,'(A175)') 'QOT file: File not designated as a CSV file (no "$" as first character) but commas found between input values for '//QOTFN(JB)
                    EXIT
                ENDIF
            ENDDO
        BACKSPACE(NPT)
        ENDIF       
             
            J       = 1
            JDAYO   = 0.0
            QSTRO   = 0.0
            DO WHILE (.TRUE.)
              IF(ICHAR1=='$')THEN
              READ (NPT,*,END=310,ERR=400) JDAY, (QSTR(JS),JS=1,NSTR(JB))
              ELSE
              READ (NPT,'(10F8.0/(9X,9F8.0))',END=310,ERR=400) JDAY, (QSTR(JS),JS=1,NSTR(JB))
              END IF

              backspace(npt)
              read(npt,'(8x,a)')blank_check

            DO JS=1,NSTR(JB)
              if(nstr(jb) < 10)then                    ! this only really works for the old format and if there are less than 10 withdrawals/branch
              if(blank_check(((js-1)*8)+1:((js-1)*8)+8)=='        ' .and. ichar1/='$')then
              CALL ERRORS
              WRITE (ERR,'(a,i3,a,i3,a,f10.3,a,a)')'QSTR is blank for JS=',js,' for branch JB=',jb,' at date [JDAY=',JDAY,'] in '//QOTFN(JB)
              endif
              endif

            IF (QSTR(JS)<0.0) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'QSTR cannot be negative:',QSTR(JS),' at date [JDAY=',JDAY,'] in '//QOTFN(JB)
            END IF
            ENDDO

              IF (JDAYO >= TMSTRT .AND. JDAYO <= TMEND) THEN
                DO JS=1,NSTR(JB)
                  QSTRMX(JB) = MAX(QSTR(JS),QSTRMX(JB))
                  IF (JDAY > TMEND) THEN
                    QSTRS(JB) = QSTRS(JB)+QSTRO(JS)*(TMEND-JDAYO)
                  ELSE
                    QSTRS(JB) = QSTRS(JB)+QSTRO(JS)*(JDAY-JDAYO)
                  END IF
                END DO
              END IF
              IF (J == 1 .AND. JDAY > TMSTRT) THEN
                CALL ERRORS
                WRITE (ERR,FMTF) 'Starting date [',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//QOTFN(JB)
              ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
                CALL ERRORS
                WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//QOTFN(JB)
              END IF
              JDAYO = JDAY
              QSTRO = QSTR
              J     = J+1
            END DO
 310        CONTINUE
            IF (JDAY < TMEND) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//QOTFN(JB)
            END IF
            CLOSE (NPT)
          ELSE
            CALL ERRORS
            WRITE (ERR,FMTI) 'Could not open branch ',JB,' outflow file '//QOTFN(JB)
          END IF
        END IF
      END IF

      IF(DYNSTRUC(JB) /= '      ON'  .AND. DYNSTRUC(JB) /= '     OFF')THEN
         CALL ERRORS
         WRITE (ERR,'(A,A8,A,I3)') 'Structures DYNELEV must be ON or OFF. DYNELEV=',DYNSTRUC(JB),' for branch:',JB
      ENDIF
      IF(DYNSTRUC(JB) == '      ON')THEN
               WRITE (SEGNUM,'(I0)') JB
               SEGNUM = ADJUSTL(SEGNUM)
               L      = LEN_TRIM(SEGNUM)
               OPEN (UNIT=NPT,FILE='dynselev'//SEGNUM(1:L)//'.npt',STATUS='OLD',IOSTAT=IERR)
          IF (IERR == 0) THEN
              AID1='dynselev'//SEGNUM(1:L)//'.npt'
               READ (NPT,'(A1)')INFORMAT
               AID='#JS line2 IN'//AID1
               READ(NPT,*,ERR=400)NJS
               IF (NJS>NST) THEN
               CALL ERRORS
               WRITE (ERR,'(A,I5,A)') 'DYNSELEV: NJS: ',NJS,' < # OF STRUCTURES [NST=',NST,'] in dynselev'//SEGNUM(1:L)//'.npt'
               END IF
               IF (NJS <= 0) THEN
               CALL ERRORS
               WRITE (ERR,'(A,I3,A)') 'DYNSELEV: NJS: ',NJS,' < =0 in dynselev'//SEGNUM(1:L)//'.npt'
               END IF
               AID='JJS line3 IN'//AID1
               DO J=1,NJS
               READ(NPT,*,ERR=400)JJS(J)
               IF(JJS(J)<=0 .OR. JJS(J) > NST)THEN
               CALL ERRORS
               WRITE (ERR,'(A,I3,A,I3,A,I3,A,I3)') 'DYNSELEV: JJS: ',JJS(J),' < =0 or > NST[',NST,'] in dynselev'//SEGNUM(1:L)//'.npt for J:',J,' and JB:',JB
               WRITE (WIN,'(A,I3,A,I3,A,I3,A,I3)') 'CRITICAL ERROR: DYNSELEV: JJS: ',JJS(J),' < =0 or > NST[',NST,'] in dynselev'//SEGNUM(1:L)//'.npt for J:',J,' and JB:',JB
               ENDIF
               ENDDO
               AID=' JDAY  ESTR IN'//AID1
               IF(INFORMAT == '$')THEN
               READ(NPT,*)
               READ (NPT,*,ERR=400)JDAY,(ESTR(JJS(J),JB), J=1,NJS)
               ELSE
               READ (NPT,'(/10F8.0:/(8X,9F8.0))',ERR=400)JDAY,(ESTR(JJS(J),JB), J=1,NJS)
               ENDIF
               ID=DS(JB)
               IF(JDAY > TMSTRT)THEN
                   CALL ERRORS
                   WRITE (ERR,'(A,f8.3,A,f8.3,A)') 'DYNESTRUC: Starting date:',JDAY,' > simulation start date:',TMSTRT,' in dynselev'//SEGNUM(1:L)//'.npt'
               ENDIF

               DO J=1,NJS
               IF (ESTR(JJS(J),JB) < EL(KB(ID+1)+1,ID)) THEN
                 CALL ERRORS
                 WRITE (ERR,'(A,F8.2,A,I3,A,F8.2,A,I3,A,F8.3)') 'DYNESTRUC:Selective withdrawal elevation [ESTR=',ESTR(JJS(J),JB),'] < the bottom active cell elevation [EL(',     &
                               KB(ID+1),')=',EL(KB(ID)+1,ID),'] for structure ',JJS(J),' at jday:',JDAY
               END IF
               IF (ESTR(JJS(J),JB) > EL(2,ID)) THEN
                 CALL ERRORS
                 WRITE (ERR,'(A,F8.2,A,I3,A,F8.2,A,I3,A,F8.3)') 'DYNESTRUC:Selective withdrawal elevation [ESTR=',ESTR(JJS(J),JB),'] > grid top layer elevation [EL(2)=',EL(2,ID),&
                               '] for structure ',JJS(J),' at jday:',JDAY
               END IF
               ENDDO
               JDAYO=JDAY
               DO WHILE(.TRUE.)
               IF(INFORMAT == '$')THEN
               READ (NPT,*,END=725,ERR=400)JDAY,(ESTR(JJS(J),JB), J=1,NJS)
               ELSE
               READ (NPT,'(10F8.0:/(8X,9F8.0))',END=725,ERR=400)JDAY,(ESTR(JJS(J),JB), J=1,NJS)
               ENDIF
               IF (JDAY <= JDAYO) THEN
               CALL ERRORS
               WRITE (ERR,'(A,F12.3,A,F12.3,A)') 'DYNESTRUC:Julian date ',JDAY,' <= previous date ',JDAYO,' in dynselev'//SEGNUM(1:L)//'.npt'
               END IF
               JDAYO = JDAY
                  DO J=1,NJS
                  IF (ESTR(JJS(J),JB) < EL(KB(ID+1)+1,ID)) THEN
                  CALL ERRORS
                  WRITE (ERR,'(A,F8.2,A,I3,A,F8.2,A,I3,A,F12.3)') 'DYNESTRUC:Selective withdrawal elevation [ESTR=',ESTR(JJS(J),JB),'] < the bottom active cell elevation [EL(',     &
                               KB(ID+1),')=',EL(KB(ID)+1,ID),'] for structure ',JJS(J),' at jday:',JDAY
                  END IF
                  IF (ESTR(JJS(J),JB) > EL(2,ID)) THEN
                  CALL ERRORS
                  WRITE (ERR,'(A,F8.2,A,I3,A,F8.2,A,I3,A,F12.3)') 'DYNESTRUC:Selective withdrawal elevation [ESTR=',ESTR(JJS(J),JB),'] > grid top layer elevation [EL(2)=',EL(2,ID),&
                               '] for structure ',JJS(J),' at jday:',JDAY
               END IF
               ENDDO
               ENDDO

725            CLOSE(NPT)
               IF(JDAY < TMEND)THEN
                   CALL ERRORS
                   WRITE (ERR,'(A,f12.3,A,f12.3,A)') 'DYNESTRUC: Ending date:',JDAY,' < simulation end date:',TMEND,' in dynselev'//SEGNUM(1:L)//'.npt'
               ENDIF

           ELSE
            CALL ERRORS
            WRITE (ERR,FMTI) 'Could not open branch ',JB,' dynamic elevation structure file: dynselev'//SEGNUM(1:L)//'.npt'
           END IF

        ENDIF

    END DO    !JB
  END DO      !JW

! Withdrawals

  WRITE (WIN,*) '  withdrawals'
  IF (WITHDRAWALS) THEN
    OPEN (UNIT=NPT,FILE=QWDFN,STATUS='OLD',IOSTAT=IERR)
    IF (IERR == 0) THEN
        AID=QWDFN
      READ(NPT,'(A1)')ICHAR1
      IF(ICHAR1=='$')WRITE(WIN,*) '     csv format for withdrawals'
      READ (NPT,'(/)')
      
        ! CHECK IF ERROR -CSV FORMAT CHECK
        IF(ICHAR1 /= '$')THEN
        READ(NPT,'(A72)')TITLE(1)
            DO JJ=1,50    ! CHECK ONLY FIRST 50 CHARACTERS
                IF(TITLE(1)(JJ:JJ) == ",")THEN
                    CALL ERRORS
                    WRITE (ERR,'(A175)') 'QWD file: File not designated as a CSV file (no "$" as first character) but commas found between input values for '//QWDFN
                    EXIT
                ENDIF
            ENDDO
        BACKSPACE(NPT)
        ENDIF       
      
      J     = 1
      JDAYO = 0.0
      DO WHILE (.TRUE.)
        IF(ICHAR1=='$')THEN
          READ (NPT,*,END=320,ERR=400) JDAY, (QWD(J,JWD),JWD=1,NWD)
          ELSE
          READ (NPT,'(10F8.0/8X,9F8.0)',END=320,ERR=400) JDAY, (QWD(J,JWD),JWD=1,NWD)
        ENDIF
        IF (JDAYO >= TMSTRT .AND. JDAY <= TMEND) THEN
          JWD = 1
          DO JB=1,NBR
            DO JWB=1,NWBR(JB)
              QWDMX(JWB,JB) = MAX(QWD(J,JWD),QWDMX(JWB,JB))
              IF (JDAY > TMEND) THEN
                QWDS(JWB,JB) = QWDS(JWB,JB)+QWD(J,JWD)*(TMEND-JDAYO)
              ELSE
                QWDS(JWB,JB) = QWDS(JWB,JB)+QWD(J,JWD)*(JDAY-JDAYO)
              END IF
              JWD = JWD+1
            END DO
          END DO
        END IF
        IF (J == 1 .AND. JDAY > TMSTRT) THEN
          CALL ERRORS
          WRITE (ERR,FMTF) 'Starting date [',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//QWDFN
        ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
          CALL ERRORS
          WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//QWDFN
        END IF
        JDAYO = JDAY
        J     = J+1
      END DO
 320  CONTINUE
      IF (JDAY < TMEND) THEN
        CALL ERRORS
        WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//QWDFN
      END IF
      CLOSE (NPT)
    ELSE
      CALL ERRORS
      WRITE (ERR,FMTA) 'Could not open withdrawal file or early end of file'//QWDFN
    END IF
  END IF

! Gates

  IF (GATES) THEN
    WRITE (WIN,*) '  gates'
    OPEN (UNIT=NPT,FILE=QGTFN,IOSTAT=IERR,STATUS='OLD')
    WRITE(WIN,*) '       reading gate file:'//qgtfn
    IF (IERR == 0) THEN
        AID=QGTFN
        
    READ(NPT,'(A1)')ICHAR1
      IF(ICHAR1=='$')WRITE(WIN,*) '     csv format for gate file'
             
    !READ(NPT,*)

    READ(NPT,'(A8)',ERR=400)GT2CHAR
      IF(GT2CHAR == 'EGT2ELEV')THEN
      BACKSPACE(NPT)
       IF(ICHAR1=='$')THEN
        READ(NPT,*,ERR=121)GT2CHAR,(EGT2(JG),JG=1,NGT)   
        ELSE
        READ(NPT,'(8X,1000F8.0)',ERR=121)(EGT2(JG),JG=1,NGT)   
      ENDIF

      DO JG=1,NGT
      IF(EGT2(JG) < 0.0)THEN
         CALL ERRORS
         WRITE (ERR,*) 'EGT2<0.0 FOR GATE #:',JG,' IN FILE:', QGTFN
      ENDIF
      IF(EGT2(JG) > EL(2,IUGT(JG)))THEN
         CALL ERRORS
         WRITE (ERR,*) 'EGT2 > TOP OF GRID FOR GATE #:',JG,' IN FILE:', QGTFN
      ENDIF
      ENDDO
      GO TO 122
      121 CALL ERRORS
      WRITE(ERR,'(a,a25)')'ERROR READING 2ND LINE IN QGT FILE WITH EGT2ELEV (see User Manual GATE FILE) in file:',QGTFN
      122 CONTINUE
      ENDIF
      READ(NPT,*)
      
     ! CHECK IF ERROR -CSV FORMAT CHECK
        IF(ICHAR1 /= '$')THEN
        READ(NPT,'(A72)')TITLE(1)
            DO JJ=1,50    ! CHECK ONLY FIRST 50 CHARACTERS
                IF(TITLE(1)(JJ:JJ) == ",")THEN
                    CALL ERRORS
                    WRITE (ERR,'(A175)') 'GATE file: File not designated as a CSV file (no "$" as first character) but commas found between input values for '//QGTFN
                    EXIT
                ENDIF
            ENDDO
        BACKSPACE(NPT)
        ENDIF              
      
!      READ (NPT,'(//)')
      J     = 1
      JDAYO = 0.0
      DO WHILE (.TRUE.)
          IF(ICHAR1=='$')THEN
                READ (NPT,*,END=120,ERR=400) JDAY, (BGATE(JG),JG=1,NGT)
          ELSE
                READ (NPT,'(1000F8.0)',END=120,ERR=400) JDAY, (BGATE(JG),JG=1,NGT)
          ENDIF
          
        IF (J == 1 .AND. JDAY > TMSTRT) THEN
          CALL ERRORS
          WRITE (ERR,FMTF) 'Starting date [',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//QGTFN
        ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
          CALL ERRORS
          WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//QGTFN
        END IF
        J     = J+1
        JDAYO = JDAY
      END DO
 120  CONTINUE
      IF (JDAY < TMEND) THEN
        CALL ERRORS
        WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//QGTFN
      END IF
      CLOSE (NPT)
    ELSE
      CALL ERRORS
      WRITE (ERR,FMTA) 'Could not open gate filename or early end of file for '//QGTFN
    END IF
  END IF

! Shading

  WRITE (WIN,*) '  shading'
  OPEN  (NPT,FILE=SHDFN,IOSTAT=IERR,STATUS='OLD')                                                                      !TC 10/22/02
  WRITE(WIN,*)'       reading SHADE file:'//shdfn
  IF (IERR == 0) THEN
      AID=SHDFN

      READ(NPT,'(A1)',ERR=400)ICHAR1
      IF(ICHAR1=='$')WRITE(WIN,*) '     csv format for SHADE file'
      READ (NPT,'(/)')
      
     ! CHECK IF ERROR -CSV FORMAT CHECK
        IF(ICHAR1 /= '$')THEN
        READ(NPT,'(A72)')TITLE(1)
            DO JJ=1,50    ! CHECK ONLY FIRST 50 CHARACTERS
                IF(TITLE(1)(JJ:JJ) == ",")THEN
                    CALL ERRORS
                    WRITE (ERR,'(A175)') 'SHADE file: File not designated as a CSV file (no "$" as first character) but commas found between input values for '//SHDFN
                    EXIT
                ENDIF
            ENDDO
        BACKSPACE(NPT)
        ENDIF       
      
      
  !  READ  (NPT,'(//)')
    SHADEIN=1.0   ! SW 7/15/14 FOR CSV FILE BLANKS ARE NOT REGISTERED AS A NUMBER, HENCE NEED TO SET ARRAY
    DO I=1,IMX
      IF(ICHAR1=='$')THEN
      READ (NPT,*,ERR=400,END=4701) J,SHADEIN(I)
            ELSE
      READ (NPT,'(8X,F8.0)',ERR=400, END=4701) SHADEIN(I)
      ENDIF
      IF (SHADEIN(I) == 0.0) THEN
        CALL   WARNINGS
        WRITE (WRN,FMTI) '100% shade for segment ',I
      ELSE IF (SHADEIN(I) > 1.0) THEN
        CALL   ERRORS
        WRITE (ERR,FMTI) 'shade > 1.0 for segment ',I
      END IF
      IF (SHADEIN(I) < 0.0) THEN
        BACKSPACE (NPT)
        IF(ICHAR1=='$')THEN
        READ (NPT,*,ERR=400) J,SHADEIN(I),TTLB(I),TTRB(I),CLLB(I),CLRB(I),SRLB1(I),SRLB2(I),SRRB1(I),SRRB2(I),(TOPO(I,JJ),JJ=1,IANGLES),     &      ! SW 7/14/14
                                  SRFJDAY1(I),SRFJDAY2(I)
         ELSE
        READ (NPT,'(16X,28F8.0)',ERR=400) TTLB(I),TTRB(I),CLLB(I),CLRB(I),SRLB1(I),SRLB2(I),SRRB1(I),SRRB2(I),(TOPO(I,J),J=1,IANGLES),     &
                                  SRFJDAY1(I),SRFJDAY2(I)
        ENDIF
        IF (TTLB(I) < EL(KB(I),I)) THEN
          CALL   ERRORS
          WRITE (ERR,FMTI) 'Left bank vegetative shade cover elevation less than bottom channel elevation for segment ',I
        END IF
        IF (TTRB(I) < EL(KB(I),I)) THEN
          CALL   ERRORS
          WRITE (ERR,FMTI) 'Right bank vegetative shade cover elevation less than bottom channel elevation for segment ',I
        END IF
        IF (CLLB(I) <= 0.0) THEN
          CALL   ERRORS
          WRITE (ERR,FMTI) 'Left bank distance to vegetative shade cover <= 0 m for segment ',I
        END IF
        IF (CLRB(I) <= 0.0) THEN
          CALL   ERRORS
          WRITE (ERR,FMTI) 'Right bank distance to vegetative shade cover <= 0 m for segment ',I
        END IF
        IF (SRLB1(I) < 0.0 .OR. SRLB1(I) > 1.0) THEN
          CALL   ERRORS
          WRITE (ERR,FMTI) 'Shading reduction factor LB1 < 0 or > 1 for segment ',I
        END IF
        IF (SRLB2(I) < 0.0 .OR. SRLB2(I) > 1.0) THEN
          CALL   ERRORS
          WRITE (ERR,FMTI) 'Shading reduction factor LB2 < 0 or > 1 for segment ',I
        END IF
        IF (SRRB1(I) < 0.0 .OR. SRRB1(I) > 1.0) THEN
          CALL   ERRORS
          WRITE (ERR,FMTI) 'Shading reduction factor RB1 < 0 or > 1 for segment ',I
        END IF
        IF (SRRB2(I) < 0.0 .OR. SRRB2(I) > 1.0) THEN
          CALL   ERRORS
          WRITE (ERR,FMTI) 'Shading reduction factor RB2 < 0 or > 1 for segment ',I
        END IF
        DO JJ=1,IANGLES
          IF(TOPO(I,JJ) < 0 .OR. TOPO(I,JJ) > 3.1419/2.0) THEN
            CALL   ERRORS
            WRITE (ERR,FMTI) 'Shading topographic angle < 0 or > 3.1419/2 for segment ',I
          END IF
        END DO
        IF (SRFJDAY1(I) > SRFJDAY2(I)) THEN
          CALL   ERRORS
          WRITE (ERR,FMTI) 'Shade reduction factor JDAY1 (leaf out) > JDAY2 (leaf off) for segment ',I
        END IF                                                                                                         !SW 04/03/02
      END IF                                                                                                           !SW 04/03/02
    END DO                                                                                                             !SW 04/03/02
  ELSE                                                                                                                 !TC 10/22/02
    CALL ERRORS                                                                                                        !TC 10/22/02
    WRITE (ERR,FMTA) 'Could not open shade filename '//SHDFN                                                           !TC 10/22/02
  END IF                                                                                                               !TC 10/22/02
  GO TO 4700
4701   CALL ERRORS                                                                                                        !TC 10/22/02
    WRITE (ERR,'(a,a,a,i3)') 'Shade file missing row(s) or shade segment information in file:',trim(SHDFN),', while reading information for segment#',I                                                           !TC 10/22/02
4700 CLOSE (NPT)                                                                                                          !SW 04/03/02

! Flux Output
     DO JW=1,NWB
         IF(FLXC(JW) /= '      ON' .AND. FLXC(JW) /= '     OFF')THEN
            CALL ERRORS
            WRITE (ERR,'(a,a8,a)') 'Flux Output Control: [FLXC=',FLXC(JW),'] /= " ON" or "OFF"  '
         ENDIF
     ENDDO
     
   ! HMULT, CMULT, CDMULT
     
  DO J=1,NHY
      IF(HMULT(J) <= 0.0)THEN
          CALL ERRORS
          WRITE (ERR,'(a,a8,a)') 'Hydraulic Multiplier: HMULT is <= 0.0 [HMULT=',HMULT(J),'] for hydraulic variable #:',J
      ENDIF
  ENDDO

   DO J=1,NCT
      IF(CMULT(J) <= 0.0)THEN
          CALL ERRORS
          WRITE (ERR,'(a,F12.5,a,I3)') 'Concentration Multiplier: CMULT is <= 0.0 [CMULT=',CMULT(J),'] for constituent #:',J
      ELSEIF(CMULT(J) > 1.0)THEN
          CALL WARNINGS
          WRITE (WRN,'(a,F12.5,a,I3)') 'Concentration Multiplier:CMULT is > 1.0. All output concentrations are multiplied by CMULT. CMULT=',CMULT(J),' for constituent #:',J
      ENDIF
   ENDDO
   
    DO J=1,NDC
      IF(CDMULT(J) <= 0.0)THEN
          CALL ERRORS
          WRITE (ERR,'(a,F12.5,a,I3)') 'Concentration Derived Multiplier: CDMULT is <= 0.0 [CDMULT=',CDMULT(J),'] for derived constituent #:',J
        ELSEIF(CDMULT(J) > 1.0)THEN
          CALL WARNINGS
          WRITE (WRN,'(a,F12.5,a,I3)') 'Concentration Derived Multiplier: CDMULT is > 1.0. All output derived concentrations are multiplied by CDMULT. CDMULT=',CDMULT(J),' for derived constituent #:',J
      ENDIF
   ENDDO 
     
! Restart file

  WRITE (WIN,*) '  restart'
  IF (RESTART_IN) THEN
    INQUIRE(FILE=rsifn, EXIST=file_exists)   ! file_exists will be TRUE if the file
    
	If(.NOT. file_exists) Then
      CALL ERRORS
      WRITE (ERR,'(a,a,a)') 'Restart control RSIC=ON but restart file [RSIFN=',trim(rsifn),'] does not exist'
	End If !OPEN (NPT,FILE=RSIFN,STATUS='OLD',IOSTAT=IERR)
    !IF (IERR /= 0) THEN
    !  CALL ERRORS
    !  WRITE (ERR,FMTA) 'Could not open restart file '//RSIFN
    !ELSE
    !  CLOSE (NPT)
    !END IF
  END IF

 
   
  
! Water balance

  JTR = 0; JWD = 0                                                                                                     !TC 07/10/03
  DO JW=1,NWB
    DO JB=BS(JW),BE(JW)
      JTR = 0; JWD = 0                                                                                                  
      DO JT=1,NTBR(JB)
        JTR           = JTR+1                                                                                          !TC 07/10/03
        QTRMXB(JB)    = MAX(QTRMXB(JB),QTRMX(JTR,JB))
        QINMXB(JB)    = MAX(QINMXB(JB),QTRMXB(JB))
        QTRAV(JTR,JB) = QTRS(JTR,JB)/(TMEND-TMSTRT)
        QTRAVB(JB)    = QTRAVB(JB)+QTRAV(JTR,JB)
      END DO
      DO J=1,NWBR(JB)
        JWD           = JWD+1                                                                                          !TC 07/10/03
        QWDMXB(JB)    = MAX(QWDMXB(JB),QWDMX(JWD,JB))
        QOTMXB(JB)    = MAX(QOTMXB(JB),QWDMXB(JB))
        QWDAV(JWD,JB) = QWDS(JWD,JB)/(TMEND-TMSTRT)
        QWDAVB(JB)    = QWDAVB(JB)+QWDAV(JWD,JB)
      END DO
      QINAV(JB)  = QINS(JB) /(TMEND-TMSTRT)
      QPRAV(JB)  = QPRS(JB) /(TMEND-TMSTRT)
      QDTAV(JB)  = QDTS(JB) /(TMEND-TMSTRT)
      QSTRAV(JB) = QSTRS(JB)/(TMEND-TMSTRT)
      QINAVB(JB) = QINAV(JB) +QDTAV(JB)+QPRAV(JB)+QTRAVB(JB)
      QOTAVB(JB) = QSTRAV(JB)+QWDAVB(JB)
      QOTAVW(JW) = QOTAVW(JW)+QOTAVB(JB)
      QINAVW(JW) = QINAVW(JW)+QINAVB(JB)
      QINMXB(JB) = MAX(QINMX(JB), QDTMX(JB))
      QOTMXB(JB) = MAX(QSTRMX(JB),QWDMXB(JB))
      QOTMXW(JW) = MAX(QOTMXW(JW),QOTMXB(JB))
      QINMXW(JW) = MAX(QINMXW(JW),QINMXB(JB))
    END DO
  END DO
!  DO JW=1,NWB
!    DO JB=BS(JW),BE(JW)
!      IF (UQ_INTERNAL(JB)) THEN
!        QINAVB(JB) = QINAVB(UQB(JB))+QDTAV(JB)+QPRAV(JB)+QTRAVB(JB)
!        QINMXB(JB) = MAX(QINMXB(UQB(JB)),QINMXB(JB))
!        QINMXW(JW) = MAX(QINMXB(UQB(JB)),QINMXB(JB))
!        QINAVW(JW) = QINAVB(UQB(JB))
!      END IF
!      IF (DAM_FLOW(JB)) THEN
!        QINAVB(JB) = QOTAVB(UQB(JB))
!        QINMXB(JB) = MAX(QOTMXB(UQB(JB)),QINMX(JB))
!        QINMXW(JW) = MAX(QINMXW(JW),QINMXB(JB))
!        QINAVW(JW) = QSTRAV(UQB(JB))+QDTAV(JB)+QPRAV(JB)+QTRAVB(JB)
!      END IF
!    END DO
!  END DO

!***********************************************************************************************************************************
!*                                                      Task 4: Control File Inputs                                               **
!***********************************************************************************************************************************

  WRITE (WIN,*) 'Control file'
  WRITE (WIN,*) '  time controls'
  IF (TMEND < TMSTRT) THEN
    CALL ERRORS
    WRITE (ERR,FMTF) 'Starting time [TMSTRT=',TMSTRT,'] < ending time [TMEND=',TMEND,']'
  END IF

! Autostepping controls

  IF (NDLT <= 0) THEN
    CALL ERRORS
    WRITE (ERR,FMTI) 'Number of timestep intervals [NDLT=',NDLT,'] <= 0'
  END IF
  IF (DLTD(1) > TMSTRT) THEN
    CALL ERRORS
    WRITE (ERR,FMTF) 'Timestep date [DLTD(1)=',DLTD(1),'] > starting time [TMSTRT=',TMSTRT,']'
  END IF
  DO J=1,NDLT-1
    IF (DLTD(J+1) < DLTD(J)) THEN
      CALL ERRORS
      WRITE (ERR,FMTIF) 'Timestep date [DLTD(',J,')=',DLTD(J),'] < timestep date [DLTD(',J+1,')=',DLTD(J+1),']'
    END IF
  END DO
  DO J=1,NDLT
    IF (DLTD(J) < 1.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTIF) 'Timestep date [DLTD(',J,     ')=',DLTD(J),  '] <  1.0'
    END IF
    IF (DLTMAX(J) <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTIF) 'Maximum timestep [DLTMAX(',J,')=',DLTMAX(J),'] <= 0.0'
    END IF
    IF (DLTF(J) <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTIF) 'Timestep fraction [DLTF(',J,  ')=',DLTF(J), '] <= 0.0'
    ELSE IF (DLTF(J) >= 1.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTIF) 'Timestep fraction [DLTF(',J,  ')=',DLTF(J), '] >= 1.0'
    END IF
  END DO
  IF (DLTINTR /= '      ON' .AND. DLTINTR /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,'(a,a8,a)') 'DLTINTR: DLT interpolation control [DLTINTR=',DLTINTR,'] /= " ON" or "OFF"  '
    END IF
  IF (DLTMIN <= 0.0) THEN
    CALL ERRORS
    WRITE (ERR,FMTF) 'Minimum timestep [DLTMIN=',DLTMIN,'] <= 0.0'
  END IF
  IF (DLTMIN > 0.0 .AND. DLTMIN < 0.1) THEN
    CALL WARNINGS
    WRITE (WRN,FMTF) 'Minimum timestep [DLTMIN=',DLTMIN,'] <  0.1'
  END IF
  DO JW=1,NWB
    IF (VISC(JW) /= '      ON' .AND. VISC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Vertical eddy viscosity limitation control [VISC='//VISC(JW)(6:8)//'] must be " ON" or "OFF" for waterbody ',JW
    END IF
    IF (CELC(JW) /= '      ON' .AND. CELC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Internal gravity wave limitation control   [CELC='//CELC(JW)(6:8)//'] must be " ON" or "OFF" for waterbody ',JW
    END IF
        IF (DLTADD(JW) /= '      ON' .AND. DLTADD(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Additional DLT criterion DLTADD [DLTADD='//DLTADD(JW)(6:8)//'] must be " ON" or "OFF" for waterbody ',JW
    END IF

  END DO

! Channel slope and minimum number of layers for active segment

  WRITE (WIN,*) '  channel slope'
  DO JB=1,NBR
    IF (SLOPE(JB) > 0.1) THEN
      CALL WARNINGS
      WRITE (WRN,FMTFI) 'Channel slope [SLOPE=',SLOPE(JB),'] > 0.1 for branch ',JB
    ELSE IF (SLOPE(JB) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Channel slope [SLOPE=',SLOPE(JB),'] < 0.0 for branch ',JB
    END IF
    IF (NL(JB) > 3) THEN
      CALL WARNINGS
      WRITE (WRN,FMTI) 'Minimum number of layers for segment to be active [NL=',NL(JB),'] > 3'
    END IF
    IF (NL(JB) <= 0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Minimum number of layers for segment to be active [NL=',NL(JB),'] <= 0'
    END IF
  END DO

! Branch linkage

  DO JW=1,NWB
    BRANCH_FOUND = .FALSE.
    DO JB=BS(JW),BE(JW)
      IF (JBDN(JW) == JB) THEN
        BRANCH_FOUND = .TRUE.
        EXIT
      END IF
    END DO
    IF (.NOT. BRANCH_FOUND) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Downstream branch [JBDN=',JBDN(JW),'] is not located in waterbody ',JW
    END IF
  END DO

! Temperature intitial condition

  WRITE (WIN,*) '  initial conditions'
  DO JW=1,NWB
    IF (T2I(JW) < -2.0) THEN
      CALL WARNINGS
      WRITE (WRN,FMTFI) 'Initial temperature [T2I=',T2I(JW),'] <= 0.0 for waterbody ',JW
    END IF
  END DO

! Ice initial condition

  DO JW=1,NWB
    IF (ICE_CALC(JW)) THEN
      IF (ICEI(JW) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Ice initial thickness [ICEI=',ICEI(JW),'] <= 0.0 for waterbody ',JW
      ELSE IF (ICEI(JW) > 0.1) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Ice initial thickness [ICEI=',ICEI(JW),'] >  0.1 for waterbody ',JW
      END IF
      IF (ICEI(JW) > 0.0 .AND. ICEI(JW) < ICEMIN(JW)) THEN
        CALL WARNINGS
        WRITE (WRN,FMT2FI) 'Ice initial thickness [ICEI=',ICEI(JW),'] < minimum ice thickness [ICEMIN=',ICEMIN(JW), '] for '//     &
                           'waterbody ',JW
      END IF
    END IF
  END DO

! Water body type

  WRITE (WIN,*) '  waterbody type'
  DO JW=1,NWB
    IF (WTYPEC(JW) == '    SALT') THEN
      IF (.NOT.CONSTITUENTS) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Since water type is SALT, constituents [CCC] must be set to "ON" for waterbody ',JW
      ELSE IF (.NOT.DISSOLVED_SOLIDS) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Since water type is SALT, CAC(4) must be set to "ON" for waterbody ',JW
      END IF
    ELSE IF (WTYPEC(JW) /= '   FRESH' .AND. WTYPEC(JW) /= '    SALT') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Waterbody type [WTYPEC(',JW,')='//WTYPEC(JW)//'] /= "   FRESH" or "    SALT" for waterbody ',JW
    END IF
  END DO

  WRITE (WIN,*) '  waterbody grid'
  DO JW=1,NWB
    IF (GRIDCC(JW) /= '    RECT' .AND. GRIDCC(JW) /= '   TRAP') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Waterbody grid [GRIDC(',JW,')='//GRIDCC(JW)//'] /= "   RECT" or "    TRAP" for waterbody ',JW
    END IF
  END DO

! Calculations

  WRITE (WIN,*) '  calculations'
  DO JW=1,NWB
    IF (VBC(JW) /= '      ON' .AND. VBC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Volume balance control   [VBC='//VBC(JW)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
    IF (MBC(JW) /= '      ON' .AND. MBC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Mass balance control     [MBC='//MBC(JW)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
    IF (EBC(JW) /= '      ON' .AND. EBC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Energy balance control   [EBC='//EBC(JW)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
    IF (PQC(JW) /= '      ON' .AND. PQC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Inflow placement control [PQC='//PQC(JW)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
    IF (EVC(JW) /= '      ON' .AND. EVC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Evaporation control      [EVC='//EVC(JW)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
    IF (PRC(JW) /= '      ON' .AND. PRC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Precipitation control    [PRC='//PRC(JW)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
  END DO

! Interpolation controls

  WRITE (WIN,*) '  interpolation'
  DO JW=1,NWB
    IF (METIC(JW) /= '      ON' .AND. METIC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Meteorology interpolation control [METIC='//METIC(JW)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
  END DO
  DO JB=1,NBR
    IF (QINIC(JB) /= '      ON' .AND. QINIC(JB) /= '     OFF' .AND. UQ_EXTERNAL(JB)) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Inflow interpolation control [QINIC='//QINIC(JB)//'] /= " ON" or "OFF" for branch ',JB
    END IF
    IF (HDIC(JB) /= '      ON' .AND. HDIC(JB) /= '     OFF' .AND. (UH_EXTERNAL(JB) .OR. DH_EXTERNAL(JB))) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Head boundary interpolation control [HDIC='//HDIC(JB)//'] /= " ON" or "OFF" for branch ',JB
    END IF
    IF (DTRIC(JB) /= '      ON' .AND. DTRIC(JB) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Distributed tributary inflow interpolation control [DTRIC='//DTRIC(JB)//'] /= " ON" or "OFF" for branch ',JB
    END IF
    DO JS=1,NSTR(JB)
      IF (DQ_EXTERNAL(JB)) THEN
        IF (STRIC(JS,JB) /= '      ON' .AND. STRIC(JS,JB) /= '     OFF') THEN
          CALL ERRORS
          WRITE (ERR,FMTI) 'Outflow interpolation control [QOUTIC='//STRIC(JS,JB)//'] /= " ON" or "OFF" for structure ',JS,        &
                           ' in branch ',JB
        END IF
      END IF
    END DO
  END DO
  DO JT=1,NTR
    IF (TRIC(JT) /= '      ON' .AND. TRIC(JT) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Tributary inflow interpolation control [TRIC='//TRIC(JT)//'] /= " ON" or "OFF" for tributary ',JT
    END IF
  END DO
  DO JWD=1,NWD
    IF (WDIC(JWD) /= '      ON' .AND. WDIC(JWD) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Withdrawal interpolation control [WDIC='//WDIC(JWD)//'] /= " ON" or "OFF" for withdrawal ',JWD
    END IF
  END DO

! Dead sea

  WRITE (WIN,*) '  dead sea'
  DO JW=1,NWB
    IF (WINDC(JW) /= '      ON' .AND. WINDC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Wind control [WINDC='//WINDC(JW)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
    IF (QINC(JW) /= '      ON' .AND. QINC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Inflow control [QINC='//QINC(JW)//'] /= " ON" or "OFF" for waterbody ', JW
    END IF
    IF (QOUTC(JW) /= '      ON' .AND. QOUTC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Outflow control [QOUTC='//QOUTC(JW)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
    IF (HEATC(JW) /= '      ON' .AND. HEATC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Surface heat exchange control [HEATC='//HEATC(JW)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
  END DO

!** Surface heat exchange

  WRITE (WIN,*) '  surface heat exchange'
  DO JW=1,NWB
    IF (SLHTC(JW) /= '    TERM' .AND. SLHTC(JW) /= '      ET') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Heat exchange solution control [SLHTC='//SLHTC(JW)//'] /= "    TERM" or "      ET" for waterbody ',JW
    END IF
    IF (SLHTC(JW) == '      ET') THEN
      CALL WARNINGS
      WRITE (WRN,FMTI) 'Heat exchange solution control [SLHTC='//SLHTC(JW)//'] is not recommended. TERM should be more accurate. Waterbody:',JW
    END IF
! Radiation and evaporation

    IF (SROC(JW) /= '      ON' .AND. SROC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Solar radiation input control [SROC='//SROC(JW)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
    IF (AFW(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Wind speed function coefficient [AFW=',AFW(JW),'] < 0 for waterbody ',JW
    END IF
    IF (BFW(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Wind speed function coefficient [BFW=',BFW(JW),'] < 0 for waterbody ',JW
    END IF
    IF (CFW(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Wind speed function coefficient [CFW=',CFW(JW),'] < 0 for waterbody ',JW
    END IF
    IF (WINDH(JW) <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Wind speed measurement height [WINDH=',WINDH(JW),'] <= 0m for waterbody ',JW
    END IF
    IF (WINDH(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Wind speed measurement height [WINDH=',WINDH(JW),'] < 0m for waterbody ',JW
    END IF
    IF (RHEVC(JW) /= '      ON' .AND. RHEVC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Ryan-Harleman evaporation control [RHEVC='//RHEVC(JW)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
    IF (FETCHC(JW) /= '      ON' .AND. FETCHC(JW) /= '     OFF') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Fetch computation control         [FETCHC='//FETCHC(JW)//'] /= " ON" or "OFF" for waterbody ',JW
    END IF
  END DO

! Ice cover

  WRITE (WIN,*) '  ice cover'
  DO JW=1,NWB
    IF (ICEC(JW) /= '      ON' .AND. ICEC(JW) /= '     OFF' .AND. ICEC(JW) /= '    ONWB') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Ice control [ICEC='//ICEC(JW)//'] /= "ON" or "OFF" or "ONWB" for waterbody ',JW
    END IF
    IF (SLICEC(JW) /= '  DETAIL' .AND. SLICEC(JW) /= '  SIMPLE') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Ice solution control [SLICEC='//SLICEC(JW)//'] /= "  DETAIL" or "  SIMPLE" for waterbody ',JW
    END IF
    IF (ALBEDO(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Ice albedo [ALBEDO=',ALBEDO(JW),'] < 0.0 for waterbody ',JW
    END IF
    IF (HWI(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Water/ice heat exchange [HWI=',HWI(JW),'] < 0.0 for waterbody ',JW
    END IF
    IF (BETAI(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Ice solar radiation absorption [BETAI=',BETAI(JW),'] < 0.0 for waterbody ',JW
    END IF
    IF (BETAI(JW) > 0.9) THEN
      CALL WARNINGS
      WRITE (WRN,FMTFI) 'Ice solar radiation absorption [BETAI=',BETAI(JW),'] > 0.9 for waterbody ',JW
    END IF
    IF (GAMMAI(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Ice solar radiation extinction [GAMMAI=',GAMMAI(JW),'] < 0.0 for waterbody ',JW
    END IF
    IF (ICEMIN(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Minimum ice thickness for ice formation [ICEMIN=',ICEMIN(JW),'] < 0.0 for waterbody ',JW
    END IF
    IF (ICEMIN(JW) > 0.1) THEN
      CALL WARNINGS
      WRITE (WRN,FMTFI) 'Minimum ice thickness for ice formation [ICEMIN=',ICEMIN(JW),'] > 0.1 for waterbody ',JW
    END IF
    IF (ICET2(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Ice formation temperature [ICET2=',ICET2(JW),'] < 0.0 for waterbody ',JW
    END IF
    IF (ICET2(JW) > 4.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Ice formation temperature [ICET2=',ICET2(JW),'] > 4.0 for waterbody ',JW
    END IF
  END DO

! Transport solution

  WRITE (WIN,*) '  transport solution'
  DO JW=1,NWB
    IF (SLTRC(JW) /= 'QUICKEST' .AND. SLTRC(JW) /= '  UPWIND' .AND. SLTRC(JW) /= 'ULTIMATE') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Transport solution [SLTRC='//SLTRC(JW)//'] /= "QUICKEST" or "  UPWIND" or "ULTIMATE" for waterbody ',JW
    END IF
    IF (THETA(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Vertical advection time weighting [THETA=',THETA(JW),'] < 0.0 for waterbody ',JW
    ELSE IF (THETA(JW) > 1.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Vertical advection time weighting [THETA=',THETA(JW),'] > 1.0 for waterbody ',JW
    END IF
  END DO

! Hydraulic coefficients

  WRITE (WIN,*) '  hydraulic coefficients'
  DO JW=1,NWB
      
    IF (AX(JW) < -10.0) THEN
      CALL WARNINGS
      WRITE (WRN,FMTFI) 'Longitudinal eddy viscosity [AX=',AX(JW),'] < -10.0 for waterbody ',JW
    ENDIF
    
    IF(AX(JW) > 0.0)THEN
    IF (AX(JW) < 0.01 .OR. AX(JW) > 100.0) THEN
      CALL WARNINGS
      WRITE (WRN,FMTFI) 'Longitudinal eddy viscosity [AX=',AX(JW),'] < 0.01 or > 100 for waterbody ',JW
    END IF
    ENDIF
    
    IF (DX(JW) < -10.0) THEN
      CALL WARNINGS
      WRITE (WRN,FMTFI) 'Longitudinal eddy diffusivity [DX=',DX(JW),'] < -10.0 for waterbody ',JW
    ENDIF
    
    IF(DX(JW) > 0.0)THEN
    IF (DX(JW) < 0.01 .OR. DX(JW) > 100.0) THEN
      CALL WARNINGS
      WRITE (WRN,FMTFI) 'Longitudinal eddy diffusivity [DX=',DX(JW),'] < 0.01 or > 100 for waterbody ',JW
    END IF
    ENDIF
    
      DXMIN=1.E6
      DO JB=BS(JW),BE(JW)
      DO I=US(JB),DS(JB)
      DXMIN=MIN(DLX(I),DXMIN)
      ENDDO
      ENDDO
      DXTHEORY=5.84e-4*DXMIN**1.1
      IF (DX(JW) > DXTHEORY*10.) THEN
      CALL WARNINGS
      WRITE (WRN,"('Longitudinal eddy diffusivity [DX=',f8.3,'] 10x > DX-theory=',f8.4,' for waterbody ',i4,':Model may be unstable')")dx(jw),dxtheory,JW
      ENDIF
      IF (AX(JW) > DXTHEORY*10.) THEN
      CALL WARNINGS
      WRITE (WRN,"('Longitudinal eddy viscosity [AX=',f8.3,'] 10x > AX-theory=',f8.4,' for waterbody ',i4,':Model may be unstable')")ax(jw),dxtheory,JW
      ENDIF
    IF (CBHE(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,'(A,E8.1,A,I0)') 'Coefficient of bottom heat exchange [CBHE=',CBHE(JW),'] < 0.0 for waterbody ',JW
    END IF
    IF (CBHE(JW) > 3.0) THEN
      CALL WARNINGS
      WRITE (WRN,'(A,E8.1,A,I0)') 'Coefficient of bottom heat exchange [CBHE=',CBHE(JW),'] > 3.0 for waterbody ',JW    !TC 09/30/03
    END IF
    IF (TSED(JW) > 30 .OR. TSED(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Sediment temperature [TSED=',TSED(JW),'] < 0 or > 30 for waterbody ',JW
    END IF
    IF (FI(JW) > 0.05) THEN
      CALL WARNINGS
      WRITE (WRN,FMTFI) 'Sediment temperature adjustment [FI=',FI(JW),'] > 0.05 for waterbody ',JW
    END IF
    IF (FI(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Sediment temperature adjustment [FI=',FI(JW),'] < 0 for waterbody ',JW
    END IF
    IF (TSEDF(JW) < 0.00 .OR. TSEDF(JW) > 1.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Sediment temperature factor [TSEDF=',TSEDF(JW),'] < 0 or > 1 for waterbody ',JW
    END IF
  END DO

! Vertical eddy viscosity formulation

  WRITE (WIN,*) '  vertical eddy viscosity'
  DO JW=1,NWB
    !NOTE ADD CHECKS FOR FBC,E,ARODI,STRCKLR,BOUNDFR,TKECAL !!

    IF (AZC(JW) /= '    NICK' .AND. AZC(JW) /= '      W2' .AND. AZC(JW) /= '     W2N' .AND. AZC(JW) /= '     RNG' .AND.            &
        AZC(JW) /= '   PARAB' .and. azc(jw) /= '     TKE' .and. azc(jw) /= '    TKE1') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Vertical eddy viscosity formulation [AZC='//AZC(JW)//'] /= "      W2", "     W2N", "    NICK",'//          &
                       ' "     RNG", or "   PARAB" or "     TKE" or "    TKE1" for waterbody ',JW
    END IF
    IF (AZC(JW) == '    TKE1')THEN
        IF(FBC(JW) < 0 .OR. FBC(JW) >3)THEN
         CALL ERRORS
         WRITE (ERR,FMTI) 'FBC for the TKE1 model [FBC=',FBC(JW),'] must be either 1, 2, or 3'//          &
                       ' for waterbody ',JW
        ENDIF
        IF(AZE(JW) /= 9.535)THEN
         CALL WARNINGS
         WRITE (WRN,FMTFI) 'E for the TKE1 model [E=',AZE(JW),'] is set at a value other than the default of 9.535'//          &
                       ' for waterbody ',JW
        ENDIF
        IF(AZE(JW) < 0.0 .or. aze(jw) > 30.0)THEN
         CALL ERRORS
         WRITE (ERR,FMTFI) 'E for the TKE1 model [E=',AZE(JW),'] is set at a value < 0 or > 30;default is 9.535'//          &
                       ' for waterbody ',JW
        ENDIF
        IF(ARODI(JW) < 0.01 .OR. ARODI(JW) > 0.7)THEN
         CALL ERRORS
         WRITE (ERR,FMTFI) 'ARODI for the TKE1 model [ARODI=',ARODI(JW),'] is set at a value < 0.01 or > 0.7; default is 0.43 for FBC=1'//          &
                       ' and 0.07 for FBC=2 (ARODI not used for FBC=3) for waterbody ',JW
        END IF
        IF(STRCKLR(JW) < 0.0 .OR. STRCKLR(JW) > 50.0)THEN
         CALL ERRORS
         WRITE (ERR,FMTFI) 'STRCKLR for the TKE1 model [STRCKLR=',STRCKLR(JW),'] is set at a value < 0.0 or > 50.; default is 24.0'//          &
                       ' for waterbody ',JW
        END IF
        IF(BOUNDFR(JW) < 0.0 .OR. BOUNDFR(JW) > 20.0)THEN
         CALL ERRORS
         WRITE (ERR,FMTFI) 'BOUNDFR for the TKE1 model [BOUNDFR=',BOUNDFR(JW),'] is set at a value < 0.0 or > 20.; default is 10.0'//          &
                       ' for waterbody ',JW
        END IF
        IF (tkecal(jw) /= '     IMP'.and. tkecal(jw) /= '     EXP') THEN
            CALL ERRORS
            WRITE (ERR,FMTI) 'TKE1 calculation technique [TKECAL='//TKECAL(JW)//'] /= "     IMP" or "     EXP" for waterbody ',JW
        END IF

    ENDIF


    IF (AZSLC(JW) /= '     IMP' .AND. AZSLC(JW) /= '     EXP') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Vertical eddy viscosity numerical solution scheme [AZSLC='//AZSLC(JW)//'] /= "IMP" or "EXP" for '//        &
                       'waterbody ',JW
    END IF
    IF (AZMAX(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Maximum vertical eddy viscosity [AZMAX]=',AZMAX(JW),'] < 0.0 for waterbody ',JW
    END IF
    IF (AZMAX(JW) > 1.0) THEN
      CALL WARNINGS
      WRITE (WRN,FMTFI) 'Maximum vertical eddy viscosity [AZMAX]=',AZMAX(JW),'] > 1.0 for waterbody ',JW
    END IF
    DO JB=BS(JW),BE(JW)
      IF (SLOPE(JB) == 0.0) THEN
        IF (AZC(JW) /= '     W2N' .AND. AZC(JW) /= '      W2' .AND. AZC(JW) /= '     TKE' .AND. AZC(JW) /= '    TKE1' .AND. (.NOT.UH_EXTERNAL(JB) .AND. .NOT.DH_EXTERNAL(JB))) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Vertical eddy viscosity formulation [AZC='//AZC(JW)//'] should probably be "      W2" or '//           &
                           '"     W2N" or "     TKE" or "   TKE1" for waterbody ',JW
        END IF
        IF (AZMAX(JW) > 0.01 .AND. AZSLC(JW) == '     EXP') THEN
          CALL WARNINGS
          WRITE (WRN,'(A,1E10.3,A,I0)') 'Explicit vertical eddy viscosity [AZC="     EXP"] and a maximum vertical eddy '//         &
                                        'viscosity of [AZMAX=',AZMAX(JW),'] will result in low timesteps for waterbody ',JW
        END IF
      END IF
    END DO
  END DO

! Bottom friction type

  WRITE (WIN,*) '  bottom friction'
  DO JW=1,NWB
    IF (FRICC(JW) /= '   CHEZY' .AND. FRICC(JW) /= '    MANN') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Bottom friction type [FRICC='//FRICC(JW)//'] /= "CHEZY" or "MANN" for waterbody ',JW
    END IF
    DO JB=BS(JW),BE(JW)
      DO I=US(JB),DS(JB)
        IF (FRICC(JW) == '   CHEZY' .AND. (FRIC(I) < 30.0 .OR. FRIC(I) > 100.0)) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Chezy coefficient [FRIC=',FRIC(I),'] < 20 or > 100 for segment ',I
        ELSE IF (FRICC(JW) == '    MANN' .AND. (FRIC(I) < 0.001 .OR. FRIC(I) > 0.20)) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Mannings N [FRIC=',FRIC(I),'] < 0.001 or > 0.2 for segment ',I
        END IF
        IF (FRICC(JW) == '   CHEZY' .AND. (FRIC(I) < 10.0 .OR. FRIC(I) > 120.0)) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Chezy coefficient [FRIC=',FRIC(I),'] < 10 or > 120 for segment ',I
        ELSE IF (FRICC(JW) == '    MANN' .AND. (FRIC(I) < 0.00001 .OR. FRIC(I) > 1.50)) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Mannings N [FRIC=',FRIC(I),'] < 0.00001 or > 1.5 for segment ',I
        END IF
      END DO
    END DO
  END DO

  ! WIND ROUGHNESS HEIGHT Z0
  WRITE (WIN,*) '  wind roughness z0'
  DO JW=1,NWB
    IF (Z0(JW) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,'("Wind roughness height Z0 < 0.0 m",1X," Z0=",F10.5," FOR JW=",I4)') Z0(JW),JW
    ELSEIF (Z0(JW) > 1.0)THEN
      CALL WARNINGS
      WRITE (WRN,'("Wind roughness height Z0 > 1.0 m",1X," Z0=",F10.5," FOR JW=",I4)') Z0(JW),JW
    ELSEIF(Z0(JW) < 1.0E-7)THEN
      CALL ERRORS
      WRITE (ERR,'("Wind roughness height Z0 < 1.E-7 m",1X," Z0=",F10.5," FOR JW=",I4)') Z0(JW),JW
    END IF
  ENDDO



! Outlet structures

  WRITE (WIN,*) '  inflow/outflow'
  WRITE (WIN,*) '    outlet structures'
  DO JW=1,NWB
    DO JB=BS(JW),BE(JW)
      IF (DQ_EXTERNAL(JB)) THEN
        ID = DS(JB)
        IF (NSTR(JB) == 0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'No structures defined [NSTR=0] for downstream flow boundary in branch ',JB
        END IF
        DO JS=1,NSTR(JB)
          KSTR = 2
          DO WHILE (EL(KSTR,DS(JB)) >= ESTR(JS,JB) .AND. KSTR < KB(ID)+1)
            KSTR = KSTR+1
          END DO
          KSTR = KSTR-1
          IF (KTSTR(JS,JB) < 2) THEN
            CALL ERRORS
            WRITE (ERR,FMTI) 'Structure ',JS,' selective withdrawal top layer [KTSTR=',KTSTR(JS,JB),'] < 2 for branch ',JB
          END IF
          IF (KTSTR(JS,JB) > KBSTR(JS,JB)) THEN
            CALL ERRORS
            WRITE (ERR,FMTI) 'Structure ',JS,' selective withdrawal top layer [KTSTR=',KTSTR(JS,JB),'] > bottom selective '//      &
                             'withdrawal layer [KBSTR=',KBSTR(JS,JB),'] in branch ',JB
          END IF
          IF (KTSTR(JS,JB) > KSTR) THEN
            CALL ERRORS
            WRITE (ERR,FMTI) 'Structure ',JS,' selective withdrawal top layer [KTSTR=',KTSTR(JS,JB),'] > centerline elevation '//  &
                             '[KSTR=',KSTR,'] in branch ',JB
          END IF
          IF (KBSTR(JS,JB) > KB(ID)) THEN
            CALL ERRORS
            WRITE (ERR,FMTI) 'Structure ',JS,' selective withdrawal bottom layer [KBSTR=',KBSTR(JS,JB),'] > bottom active layer '//&
                             '[KB=',KB(ID),'] in branch ',JB
          END IF
          IF (KBSTR(JS,JB) < KB(ID)) THEN
            CALL WARNINGS
            WRITE (WRN,FMTI) 'Structure ',JS,' selective withdrawal bottom layer [KBSTR=',KBSTR(JS,JB),'] < bottom active layer '//&
                             '[KB=',KB(ID),'] in branch ',JB
          END IF
          IF (KBSTR(JS,JB) < KSTR) THEN
            CALL ERRORS
            WRITE (ERR,FMTI) 'Structure ',JS,' Selective withdrawal bottom layer [KBSTR=',KBSTR(JS,JB),'] < centerline layer'//    &
                             ' elevation [KSTR=',KSTR,'] in branch ',JB
          END IF
          IF (SINKC(JS,JB) /= '    LINE' .AND. SINKC(JS,JB) /= '   POINT') THEN
            CALL ERRORS
            WRITE (ERR,FMTI) 'Structure ',JS,' sink type [SINKC='//SINKC(JS,JB)//'] /= "    LINE" or "   POINT" in branch ',JB
          END IF
          IF (SINKC(JS,JB) == '    LINE') THEN
            IF (WSTR(JS,JB) <= 0.0) THEN
              CALL ERRORS
              WRITE (ERR,FMTIF) 'Structure ',JS,' width [WSTR=',WSTR(JS,JB),'] <= 0 in branch ',JB
            END IF
            IF (WSTR(JS,JB) > B(KSTR,ID)) THEN
              CALL ERRORS
              WRITE (ERR,'(A,I0,A,2(F0.3,2(A,I0,A,I0,A,F0.2,A,I0)))') 'Structure ',JS,' width [WSTR=',WSTR(JS,JB),'] > cell width [B(',KSTR, ',',ID,&
                                                     ')=',B(KSTR,ID),'] in branch ',JB
            END IF
          END IF
          IF (ESTR(JS,JB) < EL(KB(ID+1)+1,ID)) THEN
            CALL ERRORS
            WRITE (ERR,FMTFI) 'Selective withdrawal elevation [ESTR=',ESTR(JS,JB),'] < the bottom active cell elevation [EL(',     &
                               KB(ID+1),')=',EL(KB(ID)+1,ID),'] for structure ',JS
          END IF
          IF (ESTR(JS,JB) > EL(2,ID)) THEN
            CALL ERRORS
            WRITE (ERR,FMT2FI) 'Selective withdrawal elevation [ESTR=',ESTR(JS,JB),'] > grid top layer elevation [EL(2)=',EL(2,ID),&
                               '] for structure ',JS
          END IF
        END DO
      ELSE IF (.NOT. DQ_EXTERNAL(JB) .AND. NSTR(JB) > 0) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Number of outlet structures [NSTR=',NSTR(JB),'] > 0 for branch ',JB,', but this branch has an internal'//&
                         ' head or flow downstream boundary'
      END IF
    END DO
  END DO

! Pipes

  WRITE (WIN,*) '    pipes'
  DO JP=1,NPI
    IF (JBUPI(JP) == 0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Upstream segment [IUPI=',IUPI(JP),'] is located in an inactive boundary segment for pipe ',JP
    END IF
    IF (JBUPI(JP) == JBDPI(JP)) THEN
      CALL WARNINGS
      WRITE (WRN,FMTI) 'Upstream [IUPI=',IUPI(JP),'] and downstream [IDPI=',IDPI(JP),'] segments are in same branch for pipe ',JP
    END IF
    IF (IUPI(JP) > 0 .AND. IUPI(JP) < IMX) THEN
      IF (EUPI(JP) > EL(2,IUPI(JP))) THEN
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Upstream elevation [EUPI=',EUPI(JP),'] > top elevation of the grid [EL=',EL(1,IUPI(JP)),'] for pipe ',JP
      END IF
      IF (EUPI(JP) < EL(KB(IUPI(JP))+1,IUPI(JP))) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Upstream elevation [EUPI=',EUPI(JP),'] < elevation of the bottom active cell for pipe ',JP
      END IF
      IF (PUPIC(JP) == ' SPECIFY') THEN
        IF (ETUPI(JP) > EL(2,IUPI(JP))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement top elevation [ETUPI=',ETUPI(JP),'] > top elevation of the grid for pipe ',JP
        END IF
        IF (ETUPI(JP) < EL(KB(IUPI(JP))+1,IUPI(JP))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement top elevation [ETUPI=',ETUPI(JP),'] < elevation of the bottom active cell'//&
                            ' for pipe ',JP
        END IF
        IF (KBUPI(JP) > KB(IUPI(JP))) THEN
          CALL ERRORS
          WRITE (ERR,FMTI) 'Upstream selective withdrawal bottom layer [KBUPI=',KBUPI(JP),'] > segment bottom active layer [KB=',  &
                            KB(IUPI(JP)),'] for pipe ',JP
        END IF
        IF (KBUPI(JP) < KB(IUPI(JP))) THEN    ! cb 9/5/14
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Upstream selective withdrawal bottom layer [KBUPI=',KBUPI(JP),'] < segment bottom active layer [KB=',  &
                            KB(IUPI(JP)),'] for pipe ',JP
        END IF
      END IF
    ELSE
      CALL ERRORS
      WRITE (ERR,FMTI) 'Upstream segment [IUPI=',IUPI(JP),'] is outside of computational grid for pipe ',JP
    END IF
    IF (WPI(JP) <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Pipe diameter [WPI=',WPI(JP),'] <= 0 for pipe ',JP
    END IF
    IF (WPI(JP) > 5.0) THEN
      CALL WARNINGS
      WRITE (WRN,FMTFI) 'Pipe diameter [WPI=',WPI(JP),'] > 5 for pipe ',JP
    END IF
    IF (DLXPI(JP) <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Pipe length [DLXPI=',DLXPI(JP),'] <= 0 for pipe ',JP
    END IF
    IF (DLXPI(JP) > 1000.0) THEN
      CALL WARNINGS
      WRITE (WRN,FMTFI) 'Pipe length [DLXPI=',DLXPI(JP),'] > 1000 m for pipe ',JP
    END IF
    IF (FPI(JP) <= 0.0 .OR. FPI(JP) > 1.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Friction factor [FPI=',FPI(JP),'] < 0 or > 1 for pipe ',JP
    END IF
    IF (FMINPI(JP) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Minor friction [FMINPI=',FMINPI(JP),'] <  0 for pipe ',JP
    END IF
    IF (FMINPI(JP) > 10.0) THEN
      CALL WARNINGS
      WRITE (WRN,FMTFI) 'Minor friction [FMINPI=',FMINPI(JP),'] > 10 for pipe ',JP
    END IF
    IF (PUPIC(JP) /= ' SPECIFY' .AND. PUPIC(JP) /= '   DISTR' .AND. PUPIC(JP) /= ' DENSITY') THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Upstream inflow placement [PUPIC='//PUPIC(JP)//'] /= "SPECIFY" or "DISTR" or "DENSITY" for pipe ',JP
    END IF
    IF (KTUPI(JP) <= 1) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Selective withdrawal top layer [KTUPI=',KTUPI(JP),'] < 2 for pipe ',JP
    END IF
    IF (KTUPI(JP) > KBUPI(JP)) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Upstream selective withdrawal top layer [KTUPI=',KTUPI(JP),'] > bottom selective withdrawal layer [KBUPI=',&
                        KBUPI(JP),'] for pipe ',JP
    END IF
    IF (KTUPI(JP) == KBUPI(JP)) THEN          ! cb 9/5/14
      CALL WARNINGS
      WRITE (WRN,FMTI) 'Upstream selective withdrawal top layer [KTUPI=',KTUPI(JP),'] = bottom selective withdrawal layer [KBUPI=',&
                        KBUPI(JP),'] for pipe ',JP
    END IF
    IF (IDPI(JP) /= 0) THEN
      IF (IDPI(JP) == US(JBDPI(JP)-1) .OR. IDPI(JP) == DS(JBDPI(JP))+1 .OR. IDPI(JP) < 2 .OR. IDPI(JP) >= IMX) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Downstream segment [IDPI=',IDPI(JP),'] is located in an inactive boundary segment for pipe ',JP
      END IF
      IF (KTDPI(JP) < 2) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Downstream selective withdrawal top layer [KTDPI=',KTDPI(JP),'] < 2 for pipe ',JP
      END IF
      IF (KTDPI(JP) > KBDPI(JP)) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Downstream selective withdrawal top layer [KTDPI=',KTDPI(JP),'] > bottom selective withdrawal layer '//  &
                         '[KBDPI=',KBDPI(JP),'] for pipe ',JP
      END IF
      IF (KBDPI(JP) > KB(IDPI(JP))) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Downstream selective withdrawal bottom layer [KBDPI=',KBDPI(JP),'] > bottom active layer [KBDPI=',       &
                                 KB(IDPI(JP)),'] for pipe ',JP
      END IF
      IF (KBDPI(JP) < KB(IDPI(JP))) THEN      ! cb 9/5/14
        CALL WARNINGS
        WRITE (WRN,FMTI) 'Downstream selective withdrawal bottom layer [KBDPI=',KBDPI(JP),'] < bottom active layer [KBDPI=',       &
                                 KB(IDPI(JP)),'] for pipe ',JP
      END IF
      IF (PDPIC(JP) /= ' SPECIFY' .AND. PDPIC(JP) /= '   DISTR' .AND. PDPIC(JP) /= ' DENSITY') THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Downstream inflow placement [PDPIC='//PDPIC(JP)//'] /= "SPECIFY" or "DISTR" or "DENSITY" for pipe ',JP
      END IF
      IF (PDPIC(JP) == ' SPECIFY') THEN
        IF (ETDPI(JP) > EL(2,IDPI(JP))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement top elevation [ETDPI=',ETDPI(JP),'] > grid top elevation for pipe ',JP
        END IF
        IF (ETDPI(JP) < EL(KB(IDPI(JP))+1,IDPI(JP))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement top elevation [ETDPI=',ETDPI(JP),'] < elevation of the bottom active '//  &
                            'cell for pipe ',JP
        END IF
      END IF
    END IF
    IF ((LATPIC(JP) /= '    DOWN') .AND. (LATPIC(JP) /= '     LAT')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Pipe lateral control [LATPIC='//LATPIC(JP)(4:8)//'] must be either "DOWN" or "LAT" for pipe ',JP
    END IF
     IF ((DYNPIPE(JP) /= '      ON') .AND. (DYNPIPE(JP) /= '     OFF')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Pipe DYNAMIC PIPE control [DYNPIPE='//DYNPIPE(JP)(4:8)//'] must be either "ON" or "OFF" for pipe ',JP
    END IF
    IF(LATPIC(JP) == '    DOWN') THEN
      IF(IUPI(JP) /= DS(JBUPI(JP)))THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Since PIPE is specified as LATPIC='//LATPIC(JP)(4:8)//', it must be at end of a branch or change this to LAT for pipe ',JP
      ENDIF
    END IF
  END DO

! Spillways/weirs

  WRITE (WIN,*) '    spillways'
  DO JS=1,NSP
    DO JB=1,NBR
      IF (IUSP(JS) == (US(JB)-1) .OR. IUSP(JS) == (DS(JB)+1) .OR. IUSP(JS) < 2 .OR. IUSP(JS) > IMX) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Upstream segment [IUSP=',IUSP(JS),'] is located in a boundary segment for spillway ',JS
      END IF
      IF ((IDSP(JS) == US(JB)-1 .OR. IDSP(JS) == DS(JB)+1 .OR. IDSP(JS) < 2 .OR. IDSP(JS) >= IMX) .AND. IDSP(JS) /= 0) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Downstream segment [IDSP=',IDSP(JS),'] is located in a boundary segment for spillway ',JS
      END IF
      DO JBB=1,NBR                                                                                                     !SW 05/16/02
        IF (IUSP(JS) >= US(JBB) .AND. IUSP(JS) <= DS(JBB)) JBU = JBB                                                   !SW 05/16/02
        IF (IDSP(JS) >= US(JBB) .AND. IDSP(JS) <= DS(JBB)) JBD = JBB                                                   !SW 05/16/02
      END DO                                                                                                           !SW 05/16/02
      IF (JBD == JBU .AND. JBU /= 0 .AND. IDSP(JS) /= 0) THEN                                                          !SW 05/16/02
        CALL ERRORS
        WRITE (ERR,FMTI) 'Upstream [IUSP=',IUSP(JS),'] and downstream [IDSP=',IDSP(JS),'] are in same branch for spillway ',JS
      END IF
    END DO
    IF (ESP(JS) > EL(2,IUSP(JS))) THEN
      CALL ERRORS
      WRITE (ERR,FMT2FI) 'Upstream centerline elevation [ESP=',ESP(JS),'] > the top elevation of the grid [EL=',EL(2,IUSP(JS)),    &
                         '] for spillway ',JS
    END IF
    IF (IUSP(JS) > 0 .AND. IUSP(JS) < IMX) THEN
      IF (ESP(JS) < EL(KB(IUSP(JS))+1,IUSP(JS)) .and. LATSPC(JS) == '     LAT') THEN                                                                  !TC 09/30/03
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Upstream centerline elevation [ESP=',ESP(JS),'] < the bottom active cell elevation [EL=',              &
                            EL(KB(IUSP(JS)),IUSP(JS)),'] for spillway ',JS                                             !SW 10/16/02
      END IF
      IF (ESP(JS) < EL(KB(IUSP(JS))+1,IUSP(JS))-slope(jbu)*dlx(iusp(js))*0.5   .and. LATSPC(JS) == '    DOWN') THEN                                                                  !TC 09/30/03
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Upstream centerline elevation [ESP=',ESP(JS),'] < the bottom active cell elevation computed at DOWN edge of cell [EL at DOWN edge=',              &
                            (EL(KB(IUSP(JS))+1,IUSP(JS))-slope(jbu)*dlx(iusp(js))*0.5),'] for spillway ',JS                                             !SW 10/16/02
      END IF
    END IF
    IF (KTUSP(JS) < 2) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Selective withdrawal top layer [KTUSP=',KTUSP(JS),'] < 2 for spillway',JS
    END IF
    IF (KTUSP(JS) > KBUSP(JS)) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Upstream selective withdrawal top layer [KTUSP=',KTUSP(JS),'] > bottom selective withdrawal layer [KBUSP=',&
                        KBUSP(JS),'] for spillway ',JS
    END IF
    IF (KBUSP(JS) > KB(IUSP(JS))) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Upstream selective withdrawal bottom layer [KBUSP=',KBUSP(JS),'] > segment bottom active layer [KB=',     &
                        KB(IUSP(JS)),'] for spillway ',JS
    END IF
    IF (KBUSP(JS) < KB(IUSP(JS))) THEN                  ! cb 9/5/14
      CALL WARNINGS
      WRITE (WRN,FMTI) 'Upstream selective withdrawal bottom layer [KBUSP=',KBUSP(JS),'] < segment bottom active layer [KB=',     &
                        KB(IUSP(JS)),'] for spillway ',JS
    END IF
    IF (IDSP(JS) /= 0) THEN
      IF (KTDSP(JS) < 2) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Selective withdrawal top layer [KTDSP=',KTDSP(JS),'] < 2 for spillway ',JS
      END IF
      IF (KTDSP(JS) > KBDSP(JS)) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Downstream selective withdrawal top layer [KTDSP=',KTDSP(JS),'] > bottom selective withdrawal layer '//  &
                         '[KBUSP=',KBDSP(JS),'] for spillway ',JS
      END IF
      IF (KBDSP(JS) > KB(IDSP(JS))) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Downstream selective withdrawal bottom layer [KBDSP=',KBDSP(JS),'] > segment bottom active layer [KB=', &
                          KB(IDSP(JS)),'] for spillway ',JS
      END IF
      IF (KBDSP(JS) < KB(IDSP(JS))) THEN
        CALL WARNINGS
        WRITE (WRN,FMTI) 'Downstream selective withdrawal bottom layer [KBDSP=',KBDSP(JS),'] < segment bottom active layer [KB=', &
                          KB(IDSP(JS)),'] for spillway ',JS
      END IF
      IF (PUSPC(JS) /= ' SPECIFY' .AND. PUSPC(JS) /= '   DISTR' .AND. PUSPC(JS) /= ' DENSITY') THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Upstream inflow placement [PUSPC='//PUSPC(JS)//'] /= "SPECIFY" or "DISTR" or "DENSITY" for spillway ',JS
      END IF
      IF (PDSPC(JS) /= ' SPECIFY' .AND. PDSPC(JS) /= '   DISTR' .AND. PDSPC(JS) /= ' DENSITY') THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Downstream inflow placement [PDSPC='//PDSPC(JS)//'] /= "SPECIFY" or "DISTR" or "DENSITY" for spillway ',JS
      END IF
      IF (PUSPC(JS) == ' SPECIFY') THEN
        IF (ETUSP(JS) > EL(2,IUSP(JS))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement top elevation [ETUSP=',ETUSP(JS),'] > grid top elevationfor spillway ',JS
        END IF
        IF (ETUSP(JS) < EL(KB(IUSP(JS)+1),IUSP(JS))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement top elevation [ETUSP=',ETUSP(JS),'] < elevation of the bottom active '//    &
                            'cell for spillway ',JS
        END IF
        IF (ETUSP(JS) < ESP(JS)) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement top elevation [ETUSP=',ETUSP(JS),'] < spillway ',JS,' elevation'
        END IF
        IF (EBUSP(JS) > EL(2,IUSP(JS))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement bottom elevation [EBUSP=',EBUSP(JS),'] > spillway ',JS,' top elevation'
        END IF
        IF (EBUSP(JS) < EL(KB(IUSP(JS)+1),IUSP(JS))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement bottom elevation [EBUSP=',EBUSP(JS),'] < elevation of the bottom active '// &
                            'cell for spillway ',JS
        END IF
        IF (EBUSP(JS) > ESP(JS)) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement bottom elevation [EBUSP=',EBUSP(JS),'] > elevation of the spillway for '//  &
                            'spillway ',JS
        END IF
        IF (ETUSP(JS) < EBUSP(JS)) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement top elevation [ETUSP=',ETUSP(JS),'] < elevation of the bottom inflow '//    &
                            'placement for spillway ',JS
        END IF
      END IF
      IF (PDSPC(JS) == ' SPECIFY') THEN
        IF (ETDSP(JS) > EL(2,IDSP(JS))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement top elevation [ETDSP=',ETDSP(JS),'] > grid top elevation for spillway ',JS
        END IF
        IF (ETDSP(JS) < EL(KB(IDSP(JS))+1,IDSP(JS))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement top elevation [ETDSP=',ETDSP(JS),'] < bottom active cell elevation for '//&
                            'spillway ',JS
        END IF
        IF (ETDSP(JS) < ESP(JS)) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement top elevation [ETDSP=',ETDSP(JS),'] < spillway elevation for spillway ',JS
        END IF
        IF (EBDSP(JS) > EL(2,IDSP(JS))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement bottom elevation [EBDSP=',EBDSP(JS),'] > spillway top elevation for '//   &
                            'spillway ',JS
        END IF
        IF (EBDSP(JS) < EL(KB(IDSP(JS))+1,IDSP(JS))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement bottom elevation [EBDSP=',EBDSP(JS),'] < bottom active cell elevation '// &
                            'for spillway ',JS
        END IF
        IF (EBDSP(JS) > ESP(JS)) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement bottom elevation [EBDSP=',EBDSP(JS),'] > spillway elevation for '//       &
                            'spillway ',JS
        END IF
        IF (ETDSP(JS) < EBDSP(JS)) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement top elevation [ETDSP=',ETDSP(JS),'] < bottom inflow placement elevation'//&
                            ' for spillway ',JS
        END IF
      END IF
    END IF
    IF ((LATSPC(JS) /= '    DOWN') .AND. (LATSPC(JS) /= '     LAT')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Spillway lateral control [LATSPC='//LATSPC(JS)(4:8)//'] must be either "DOWN" or "LAT" for spillway ',JS
    END IF
    IF (LATSPC(JS) == '    DOWN') THEN     ! SW 3/26/10
      iflag=0
      DO JB=1,NBR
        IF (IUSP(JS) == DS(JB)) THEN
          iflag=1
          exit
        END IF
      END DO
        IF(IFLAG == 0)THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Since SPILLWAY is specified as LATSPC='//LATSPC(JS)(4:8)//', it must be at end of a branch or change this to LAT for spillway ',JS
        ENDIF
    END IF

!** Coefficients

    IF (A1SP(JS) <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Spillway coefficient a1 [A1SP] <= 0 for spillway ',JS
    END IF
    IF (A2SP(JS) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Spillway coefficient a2 [A2SP] < 0 for spillway ',JS
    END IF
    IF (B1SP(JS) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Spillway coefficient b1 [B1SP] < 0 for spillway ',JS
    END IF
    IF (B2SP(JS) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Spillway coefficient b2 [B2SP] < 0 for spillway ',JS
    END IF
    IF (B1SP(JS) == 0.0) THEN
      CALL WARNINGS
      WRITE (WRN,FMTI) 'Spillway coefficient b1 [B1SP] = 0 for spillway ',JS
    END IF
    IF (B1SP(JS) > 2.0) THEN
      CALL WARNINGS
      WRITE (WRN,FMTI) 'Spillway coefficient b1 [B1SP] > 2 for spillway ',JS
    END IF

!** Total dissolved gas equations

    IF ((GASSPC(JS) /= '      ON') .AND. (GASSPC(JS) /= '     OFF')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Spillway gas control [GASSPC='//GASSPC(JS)(6:8)//'] must be either " ON" or "OFF" for spillway ',JS
    END IF
 !   IF (IDSP(JS) /= 0 .AND. GASSPC(JS) == '      ON') THEN
     IF (GASSPC(JS) == '      ON') THEN    ! SW 7/2012
      IF (EQSP(JS) == 1) THEN
        IF (AGASSP(JS) < 0.0 .OR. AGASSP(JS) > 10.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Spillway gas coefficient a [AGASSP] < 0 or > 10 for spillway ',JS
        END IF
        IF (BGASSP(JS) < 100.0 .OR. BGASSP(JS) > 150.) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Spillway gas coefficient b [BGASSP] < 100 or > 150 for spillway ',JS
        END IF
      ELSE IF (EQSP(JS) == 2) THEN
        IF (AGASSP(JS) < 100.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Spillway gas coefficient a [AGASSP] < 100 for spillway ',JS
        END IF
        IF (BGASSP(JS) > 0.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Spillway gas coefficient b [BGASSP] > 0 for spillway ',JS
        END IF
        IF (CGASSP(JS) > 0.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Spillway gas coefficient c [CGASSP] > 0 for spillway ',JS
        END IF
      ELSE IF (EQSP(JS) == 3) THEN
        IF (AGASSP(JS) < 0.65 .OR. AGASSP(JS) > 1.85) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Spillway gas coefficient a [AGASSP] < 0.65 or > 1.85 for spillway ',JS
        END IF
        IF (BGASSP(JS) > 1.0 .OR. BGASSP(JS) < 0.05) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Spillway gas coefficient b [BGASSP] > 1 or < 0.05 for spillway ',JS
        END IF
        IF (CGASSP(JS) < 0.0 .OR. CGASSP(JS) > 10.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Spillway gas coefficient c [CGASSP] < 0 or > 10 for spillway ',JS
        END IF
    ELSE IF (EQSP(JT) == 4) THEN
        IF (AGASSP(JS) < 0.5 .OR. AGASSP(JS) > 1.5) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Spillway gas coefficient a [AGASSP] < 0.5 or > 1.5 for spillway ',JS
        END IF
        IF (BGASSP(JS) /= 1.0 .AND. BGASSP(JS) /= 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTI) 'Spillway gas coefficient b [BGASSP] must be 1 or 0 for EQN 4 for spillway ',JS
        END IF
        IF (CGASSP(JS) /= 1.0 .AND. CGASSP(JS) /= 0.0) THEN
          CALL ERRORS
          WRITE (WRN,FMTI) 'Spillway gas coefficient c [CGASSP] must be 1 or 0 for EQN 4 for spillway ',JS
        END IF
        
        IF(BGASSP(JS)==1)THEN
         WRITE (SEGNUM,'(I0)') JS
         SEGNUM = ADJUSTL(SEGNUM)  
         L = LEN_TRIM(SEGNUM) 
         TITLE240='w2_sp'//SEGNUM(1:L)//'DO.csv'
         UNIT=2450
         OPEN (UNIT, FILE=TRIM(TITLE240),STATUS='OLD')
         READ(UNIT,*,END=1350,ERR=400);READ(UNIT,*,END=1350,ERR=400);READ(UNIT,*,END=1350,ERR=400)   ! SKIP 3 LINES
         
         DO J=1,10000
         READ(UNIT,*,END=1350,ERR=400)JDAY,AGASSP(JS)
              IF (J == 1 .AND. JDAY > TMSTRT) THEN
                CALL ERRORS
                WRITE (ERR,FMTF) 'Starting date [',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//TRIM(TITLE240)
              ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
                CALL ERRORS
                WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//TRIM(TITLE240)
              END IF
              JDAYO = JDAY
              IF(AGASSP(JS) < 0.5  .OR. AGASSP(JS) > 1.5)THEN
                  CALL WARNINGS
                  WRITE (ERR,'(A,F10.3,A)') 'Spillway gas coefficient AGASSP < 0.5 or > 1.5 for JDAY:',JDAY,' in '//TRIM(TITLE240)
              ENDIF
         ENDDO
 1350    CONTINUE   
                IF (JDAY < TMEND) THEN
                CALL ERRORS
                WRITE (ERR,FMTF) 'Ending date [',JDAY,'] < simulation ending date [TMSTRT=',TMEND,'] in '//TRIM(TITLE240)
                ENDIF
                CLOSE(2450)
        ENDIF        
      ELSE
        CALL ERRORS
        WRITE (ERR,FMTI) 'Spillway gas equation number [EQSP]=',EQSP(JS),'] must be between 1 and 4 for spillway ',JS
      END IF
    END IF
  END DO

! Gates

  WRITE (WIN,*) '    gates'
  DO JT=1,NGT
    DO JB=1,NBR
      IF (IUGT(JT) == (US(JB)-1) .OR. IUGT(JT) == (DS(JB)+1) .OR. IUGT(JT) < 2 .OR. IUGT(JT) > IMX) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Upstream segment [IUGT=',IUGT(JT),'] is located in a boundary segment for gate ',JT
      END IF
      IF ((IDGT(JT) == US(JB)-1 .OR. IDGT(JT) == DS(JB)+1 .OR. IDGT(JT) < 2 .OR. IDGT(JT) >= IMX) .AND. IDGT(JT) /= 0) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Downstream segment [IDGT=',IDGT(JT),'] is located in a boundary segment for gate ',JT
      END IF
      DO JBB=1,NBR                                                                                                     !SW 05/16/02
        IF (IUGT(JT) >= US(JBB) .AND. IUGT(JT) <= DS(JBB)) JBU = JBB                                                   !SW 05/16/02
        IF (IDGT(JT) >= US(JBB) .AND. IDGT(JT) <= DS(JBB)) JBD = JBB                                                   !SW 05/16/02
      END DO                                                                                                           !SW 05/16/02
      IF (JBD == JBU .AND. JBU /= 0 .AND. IDGT(JT) /= 0) THEN                                                          !SW 05/16/02
        CALL ERRORS
        WRITE (ERR,FMTI) 'Upstream [IUGT=',IUGT(JT),'] and downstream [IDGT=',IDGT(JT),'] are in same branch for gate ',JT
      END IF
    END DO
    IF (EGT(JT) > EL(2,IUGT(JT))) THEN
      CALL ERRORS
      WRITE (ERR,FMT2FI) 'Upstream centerline elevation [EGT=',EGT(JT),'] > the top elevation of the grid [EL=',EL(2,IUGT(JT)),    &
                         '] for gate ',JT
    END IF
    IF (IUGT(JT) > 0 .AND. IUGT(JT) < IMX) THEN
      IF (EGT(JT) < EL(KB(IUGT(JT))+1,IUGT(JT))   .and. LATGTC(JT) == '     LAT') THEN                                                                  !TC 09/30/03
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Upstream centerline elevation [EGT=',EGT(JT),'] < the bottom active cell elevation [EL=',              &
                            EL(KB(IUGT(JT)),IUGT(JT)),'] for gate ',JT                                             !SW 10/16/02
      END IF
      IF (EGT(JT) < EL(KB(IUGT(JT))+1,IUGT(JT))-slope(jbu)*dlx(iugt(jt))*0.5   .and. LATGTC(JT) == '    DOWN') THEN                                                                  !TC 09/30/03
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Upstream centerline elevation [EGT=',EGT(JT),'] < the bottom active cell elevation computed at DOWN edge of cell [EL at DOWN edge=',              &
                            (EL(KB(IUGT(JT))+1,IUGT(JT))-slope(jbu)*dlx(iugt(jt))*0.5),'] for gate ',JT                                             !SW 10/16/02
      END IF
    END IF
    IF (KTUGT(JT) < 2) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Selective withdrawal top layer [KTUGT=',KTUGT(JT),'] < 2 for gate',JT
    END IF
    IF (KTUGT(JT) > KBUGT(JT)) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Upstream selective withdrawal top layer [KTUGT=',KTUGT(JT),'] > bottom selective withdrawal layer [KBUGT=',&
                        KBUGT(JT),'] for gate ',JT
    END IF
    IF (KBUGT(JT) > KB(IUGT(JT))) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Upstream selective withdrawal bottom layer [KBUGT=',KBUGT(JT),'] > segment bottom active layer [KB=',     &
                        KB(IUGT(JT)),'] for gate ',JT
    END IF
    IF (KBUGT(JT) < KB(IUGT(JT))) THEN                          ! cb 9/5/14
      CALL Warnings
      WRITE (WRN,FMTI) 'Upstream selective withdrawal bottom layer [KBUGT=',KBUGT(JT),'] < segment bottom active layer [KB=',     &
                        KB(IUGT(JT)),'] for gate ',JT
    END IF
    IF (IDGT(JT) /= 0) THEN
      IF (KTDGT(JT) < 2) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Selective withdrawal top layer [KTDGT=',KTDGT(JT),'] < 2 for gate ',JT
      END IF
      IF (KTDGT(JT) > KBDGT(JT)) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Downstream selective withdrawal top layer [KTDGT=',KTDGT(JT),'] > bottom selective withdrawal layer '//  &
                         '[KBUGT=',KBDGT(JT),'] for gate ',JT
      END IF
      IF (KBDGT(JT) > KB(IDGT(JT))) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Downstream selective withdrawal bottom layer [KBDGT=',KBDGT(JT),'] > segment bottom active layer [KB=', &
                          KB(IDGT(JT)),'] for gate ',JT
      END IF
      IF (KBDGT(JT) < KB(IDGT(JT))) THEN                          ! cb 9/5/14
        CALL Warnings
        WRITE (WRN,FMTI) 'Downstream selective withdrawal bottom layer [KBDGT=',KBDGT(JT),'] < segment bottom active layer [KB=', &
                          KB(IDGT(JT)),'] for gate ',JT
      END IF
      IF (PUGTC(JT) /= ' SPECIFY' .AND. PUGTC(JT) /= '   DISTR' .AND. PUGTC(JT) /= ' DENSITY') THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Upstream inflow placement [PUGTC='//PUGTC(JT)//'] /= "SPECIFY" or "DISTR" or "DENSITY" for gate ',JT
      END IF
      IF (PDGTC(JT) /= ' SPECIFY' .AND. PDGTC(JT) /= '   DISTR' .AND. PDGTC(JT) /= ' DENSITY') THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Downstream inflow placement [PDGTC='//PDGTC(JT)//'] /= "SPECIFY" or "DISTR" or "DENSITY" for gate ',JT
      END IF
      IF (PUGTC(JT) == ' SPECIFY') THEN
        IF (ETUGT(JT) > EL(2,IUGT(JT))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement top elevation [ETUGT=',ETUGT(JT),'] > grid top elevation for gate ',JT
        END IF
        IF (ETUGT(JT) < EL(KB(IUGT(JT))+1,IUGT(JT))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement top elevation [ETUGT=',ETUGT(JT),'] < elevation of the bottom active '//    &
                            'cell for gate  ',JT
        END IF
        IF (ETUGT(JT) < EGT(JT)) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement top elevation [ETUGT=',ETUGT(JT),'] < gate ',JT,' elevation'
        END IF
        IF (EBUGT(JT) > EL(2,IUGT(JT))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement bottom elevation [EBUGT=',EBUGT(JT),'] > gate ',JT,' top elevation'
        END IF
        IF (EBUGT(JT) < EL(KB(IUGT(JT)+1),IUGT(JT))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement bottom elevation [EBUGT=',EBUGT(JT),'] < elevation of the bottom active '// &
                            'cell for gate ',JT
        END IF
        IF (EBUGT(JT) > EGT(JT)) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement bottom elevation [EBUGT=',EBUGT(JT),'] > elevation of the gate for '//  &
                            'gate ',JT
        END IF
        IF (ETUGT(JT) < EBUGT(JT)) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Upstream inflow placement top elevation [ETUGT=',ETUGT(JT),'] < elevation of the bottom inflow '//    &
                            'placement for gate ',JT
        END IF
      END IF
      IF (PDGTC(JT) == ' SPECIFY') THEN
        IF (ETDGT(JT) > EL(2,IDGT(JT))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement top elevation [ETDGT=',ETDGT(JT),'] > grid top elevation for gate ',JT
        END IF
        IF (ETDGT(JT) < EL(KB(IDGT(JT))+1,IDGT(JT))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement top elevation [ETDGT=',ETDGT(JT),'] < bottom active cell elevation for '//&
                            'gate ',JT
        END IF
        IF (ETDGT(JT) < EGT(JT)) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement top elevation [ETDGT=',ETDGT(JT),'] < gate elevation for gate ',JT
        END IF
        IF (EBDGT(JT) > EL(2,IDGT(JT))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement bottom elevation [EBDGT=',EBDGT(JT),'] > gate top elevation for '//   &
                            'gate ',JT
        END IF
        IF (EBDGT(JT) < EL(KB(IDGT(JT))+1,IDGT(JT))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement bottom elevation [EBDGT=',EBDGT(JT),'] < bottom active cell elevation '// &
                            'for gate ',JT
        END IF
        IF (EBDGT(JT) > EGT(JT)) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement bottom elevation [EBDGT=',EBDGT(JT),'] > gate elevation for '//       &
                            'gate ',JT
        END IF
        IF (ETDGT(JT) < EBDGT(JT)) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Downstream inflow placement top elevation [ETDGT=',ETDGT(JT),'] < bottom inflow placement elevation'//&
                            ' for gate ',JT
        END IF
      END IF
    END IF


    IF ((LATGTC(JT) /= '    DOWN') .AND. (LATGTC(JT) /= '     LAT')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Gate lateral control [LATGTC='//LATGTC(JT)(4:8)//'] must be either "DOWN" or "LAT" for gate ',JT
    END IF
    IF ((DYNGTC(JT) /= '    FLOW') .AND. (DYNGTC(JT) /= '     ZGT').AND. (DYNGTC(JT) /= '       B')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Gate dynamic variable control [DYNGTC='//DYNGTC(JT)(4:8)//'] must be either "FLOW" or "ZGT" or "B" for gate ',JT
    END IF
    IF ((GTIC(JT) /= '     OFF') .AND. (GTIC(JT) /= '      ON')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Gate interpolation control [GTIC='//GTIC(JT)(5:8)//'] must be either "ON" or "OFF" for gate ',JT
    END IF
    IF (LATGTC(JT) == '    DOWN') THEN                   ! SW 3/26/10
      IF(IUGT(JT) /= DS(JBUGT(JT)))THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Since GATE is specified as LATGTC='//LATGTC(JT)(4:8)//', it must be at end of a branch or change this to LAT for gate ',JT
      ENDIF
    END IF


!** Coefficients

    IF (A1GT(JT) <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Gate coefficient a1 [A1GT] <= 0 for gate ',JT
    END IF
    IF (A2GT(JT) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Gate coefficient a2 [A2GT] < 0 for gate ',JT
    END IF
    IF (B1GT(JT) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Gate coefficient b1 [B1GT] < 0 for gate ',JT
    END IF
    IF (B2GT(JT) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Gate coefficient b2 [B2GT] < 0 for gate ',JT
    END IF
    IF (B1GT(JT) == 0.0) THEN
      CALL WARNINGS
      WRITE (WRN,FMTI) 'Gate coefficient b1 [B1GT] = 0 for gate ',JT
    END IF
    IF (B1GT(JT) > 2.0) THEN
      CALL WARNINGS
      WRITE (WRN,FMTI) 'Gate coefficient b1 [B1GT] > 2 for gate ',JT
    END IF
       IF (G1GT(JT) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Gate coefficient g1 [G1GT] < 0 for gate ',JT
    END IF
       IF (G2GT(JT) < 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Gate coefficient g2 [G2GT] < 0 for gate ',JT
    END IF
     IF (GTA1(JT) <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Gate coefficient gta1 [GTA1] <= 0 for gate ',JT
    END IF
     IF (GTA2(JT) <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Gate coefficient gta2 [GTA2] <= 0 for gate ',JT
    END IF
      IF (GTB1(JT) <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Gate coefficient gtb1 [GTB1] <= 0 for gate ',JT
    END IF
     IF (GTB2(JT) <= 0.0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Gate coefficient gtb2 [GTB2] <= 0 for gate ',JT
    END IF

!** Total dissolved gas equations

    IF ((GASGTC(JT) /= '      ON') .AND. (GASGTC(JT) /= '     OFF')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Gate gas control [GASGTC='//GASGTC(JT)(6:8)//'] must be either " ON" or "OFF" for gate ',JT
    END IF
 !   IF (IDGT(JT) /= 0 .AND. GASGTC(JT) == '      ON') THEN
     IF (GASGTC(JT) == '      ON') THEN                ! SW 7/2012
      IF (EQGT(JT) == 1) THEN
        IF (AGASGT(JT) < 0.0 .OR. AGASGT(JT) > 10.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Gate gas coefficient a [AGASGT] < 0 or > 10 for gate ',JT
        END IF
        IF (BGASGT(JT) < 100.0 .OR. BGASGT(JT) > 150.) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Gate gas coefficient b [BGASGT] < 100 or > 150 for gate ',JT
        END IF
      ELSE IF (EQGT(JT) == 2) THEN
        IF (AGASGT(JT) < 100.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Gate gas coefficient a [AGASGT] < 100 for gate ',JT
        END IF
        IF (BGASGT(JT) > 0.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Gate gas coefficient b [BGASGT] > 0 for gate ',JT
        END IF
        IF (CGASGT(JT) > 0.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Gate gas coefficient c [CGASGT] > 0 for gate ',JT
        END IF
      ELSE IF (EQGT(JT) == 3) THEN
        IF (AGASGT(JT) < 0.65 .OR. AGASGT(JT) > 1.85) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Gate gas coefficient a [AGASGT] < 0.65 or > 1.85 for gate ',JT
        END IF
        IF (BGASGT(JT) > 1.0 .OR. BGASGT(JT) < 0.05) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Gate gas coefficient b [BGASGT] > 1 or < 0.05 for gate ',JT
        END IF
        IF (CGASGT(JT) < 0.0 .OR. CGASGT(JT) > 10.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Gate gas coefficient c [CGASGT] < 0 or > 10 for gate ',JT
        END IF
        ELSE IF (EQGT(JT) == 4) THEN
        IF (AGASGT(JT) < 0.5 .OR. AGASGT(JT) > 1.5) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Gate gas coefficient a [AGASGT] < 0.5 or > 1.5 for gate ',JT
        END IF
        IF (BGASGT(JT) /= 1.0 .AND. BGASGT(JT) /= 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTI) 'Gate gas coefficient b [BGASGT] must be 1 or 0 for EQN 4 for gate ',JT
        END IF
        IF (CGASGT(JT) /= 1.0 .AND. CGASGT(JT) /= 0.0) THEN
          CALL ERRORS
          WRITE (WRN,FMTI) 'Gate gas coefficient c [CGASGT] must be 1 or 0 for EQN 4 for gate ',JT
        END IF
        
        IF(BGASGT(JT)==1)THEN
         WRITE (SEGNUM,'(I0)') JT
         SEGNUM = ADJUSTL(SEGNUM)  
         L = LEN_TRIM(SEGNUM) 
         TITLE240='w2_gt'//SEGNUM(1:L)//'DO.csv'
         UNIT=2450
         OPEN (UNIT, FILE=TRIM(TITLE240),STATUS='OLD')
         READ(UNIT,*,END=1360,ERR=400);READ(UNIT,*,END=1360,ERR=400);READ(UNIT,*,END=1360,ERR=400)   ! SKIP 3 LINES
         
         DO J=1,10000
         READ(UNIT,*,END=1360,ERR=400)JDAY,AGASGT(JT)
              IF (J == 1 .AND. JDAY > TMSTRT) THEN
                CALL ERRORS
                WRITE (ERR,FMTF) 'Starting date [',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//TRIM(TITLE240)
              ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
                CALL ERRORS
                WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//TRIM(TITLE240)
              END IF
              JDAYO = JDAY
              IF(AGASGT(JT) < 0.5  .OR. AGASGT(JT) > 1.5)THEN
                  CALL WARNINGS
                  WRITE (ERR,'(A,F10.3,A)') 'Gate gas coefficient AGASGT for EQN=4 < 0.5 or > 1.5 for JDAY:',JDAY,' in '//TRIM(TITLE240)
              ENDIF
         ENDDO
 1360    CONTINUE   
                IF (JDAY < TMEND) THEN
                CALL ERRORS
                WRITE (ERR,FMTF) 'Ending date [',JDAY,'] < simulation ending date [TMSTRT=',TMEND,'] in '//TRIM(TITLE240)
                ENDIF
        CLOSE(2450)
        ENDIF
        
        
        
        
        
        
        
      ELSE
        CALL ERRORS
        WRITE (ERR,FMTI) 'Gate gas equation number [EQGT]=',EQGT(JT),'] must be between 1 and 4 for gate ',JT
      END IF
    END IF
  END DO




! Pumps

  WRITE (WIN,*) '    pumps'
  DO JP=1,NPU

  IF(DYNPUMP(JP) /= '      ON'  .AND. DYNPUMP(JP) /= '     OFF')THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Pump dynamic control [DYNPUMP='//DYNPUMP(JP)(4:8)//'] must be either "ON" or "OFF" for pump ',JP
  END IF

  IF(DYNPUMP(JP) == '      ON')THEN
   WRITE (SEGNUM,'(I0)') JP
     SEGNUM = ADJUSTL(SEGNUM)
     L      = LEN_TRIM(SEGNUM)
     NPT=NPT+1
     OPEN (NPT,FILE='dynpump'//SEGNUM(1:L)//'.npt',STATUS='OLD',IOSTAT=I)
     IF (I /= 0) THEN
     WRITE (WIN,*) 'Could not open DYNPUMP FILE="dynpumpX.npt" for pump#:',JP
     WRITE (ERR,*) 'Could not open DYNPUMP FILE="dynpumpX.npt" for pump#:',JP
     STOP
     ENDIF

     !READ(NPT,*)
     !READ(NPT,*)
     !READ(NPT,*)
     !DO J=1,1000
     !READ(NPT,'(4F8.0)',END=7200)JDAY,EPU1,EON,EOFF,QP
     READ(NPT,'(A1)')ICHAR1
     READ(NPT,*)
     READ(NPT,*)
     DO J=1,15000
            IF(ICHAR1=='$')THEN
             READ(NPT,*,END=7200)JDAY,EPU1,EON,EOFF,QP
             ELSE
             READ(NPT,'(4F8.0)',END=7200)JDAY,EPU1,EON,EOFF,QP
             ENDIF 


      IF (J == 1 .AND. JDAY > TMSTRT) THEN
          CALL ERRORS
          WRITE (ERR,FMTF) 'DYNPUMP:Starting date [',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//'dynpump'//SEGNUM(1:L)//'.npt'
        ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
          CALL ERRORS
          WRITE (ERR,FMTF) 'DYNPUMP:Julian date ',JDAY,' <= previous date of ',JDAYO,' in '//'dynpump'//SEGNUM(1:L)//'.npt'
      END IF
      JDAYO=JDAY


     IF (EPU1 > EL(2,IUPU(JP))) THEN
      CALL ERRORS
      WRITE (ERR,FMTF2I) 'DYNPUMP:Upstream pump elevation [EPU=',EPU1,'] > layer 2 for pump ',JP,' FOR JD:',INT(JDAY)
    END IF
    IF (EPU1 < EL(KB(IUPU(JP))+1,IUPU(JP))) THEN
      CALL ERRORS
      WRITE (ERR,FMTF2I) 'DYNPUMP:Upstream pump elevation [EPU=',EPU1,'] < bottom active layer [KB=',KB(IUPU(JP)+1),'] for pump ',JP
    END IF
    IF (EON > EL(2,IUPU(JP))) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'DYNPUMP:Upstream elevation for turning pump on [EONPU=',EON,'] > layer 2 for pump ',JP
    END IF
    IF (EON < EL(KB(IUPU(JP))+1,IUPU(JP))) THEN
      CALL ERRORS
      WRITE (ERR,FMTF2I) 'DYNPUMP:Upstream elevation for turning pump on [EONPU=',EON,'] < bottom active layer [KB=',KB(IUPU(JP)),   &
                         '] for pump ',JP
    END IF
    IF (EON <= EOFF) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'DYNPUMP:Upstream elevation for turning pump on [EONPU=',EON,'] <= the elevation for turning the pump off '// &
                        'for pump ',JP
    END IF
    IF (EOFF > EL(2,IUPU(JP))) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'DYNPUMP:Upstream elevation for turning pump off [EOFFPU=',EOFF,'] > grid top elevation for pump ',JP
    END IF
    IF (EOFF < EL(KB(IUPU(JP))+1,IUPU(JP))) THEN
      CALL ERRORS
      WRITE (ERR,FMTF2I) 'DYNPUMP:Upstream elevation for turning pump off [EOFFPU=',EOFF,'] < bottom active layer [KB=',KB(IUPU(JP)),&
                         '] for pump ',JP
    END IF


     ENDDO
  ENDIF

7200 DO JB=1,NBR
      IF (IUPU(JP) == US(JB)-1 .OR. IUPU(JP) == DS(JB)+1 .OR. IUPU(JP) <= 0 .OR. IUPU(JP) >= IMX) THEN                 ! RA 8/14/07
        CALL ERRORS
        WRITE (ERR,FMTI) 'Upstream pump segment [IUPU=',IUPU(JP),'] is located in a boundary segment for pump ',JP
      END IF
      IF (IDPU(JP) == US(JB)-1 .OR. IDPU(JP) == DS(JB)+1 .OR. IDPU(JP) < 0 .OR. IDPU(JP) >= IMX) THEN                  !TC 09/30/03
        CALL ERRORS
        WRITE (ERR,FMTI) 'Downstream pump segment [IDPU=',IDPU(JP),'] is located in a boundary segment for pump ',JP
      END IF
    END DO
    IF (EPU(JP) > EL(2,IUPU(JP))) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Upstream pump elevation [EPU=',EPU(JP),'] > layer 2 for pump ',JP
    END IF
    IF (EPU(JP) < EL(KB(IUPU(JP))+1,IUPU(JP))) THEN
      CALL ERRORS
      WRITE (ERR,FMTF2I) 'Upstream pump elevation [EPU=',EPU(JP),'] < bottom active layer [KB=',KB(IUPU(JP)+1),'] for pump ',JP
    END IF
    IF (EONPU(JP) > EL(2,IUPU(JP))) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Upstream elevation for turning pump on [EONPU=',EONPU(JP),'] > layer 2 for pump ',JP
    END IF
    IF (EONPU(JP) < EL(KB(IUPU(JP))+1,IUPU(JP))) THEN
      CALL ERRORS
      WRITE (ERR,FMTF2I) 'Upstream elevation for turning pump on [EONPU=',EONPU(JP),'] < bottom active layer [KB=',KB(IUPU(JP)),   &
                         '] for pump ',JP
    END IF
    IF (EONPU(JP) <= EOFFPU(JP)) THEN
      IF(.NOT.PUMP_DOWNSTREAM(JP))THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Upstream elevation for turning pump on [EONPU=',EONPU(JP),'] < the elevation for turning the pump off '// &
                        'for pump ',JP
      ENDIF
    END IF
    IF (EOFFPU(JP) > EL(2,IUPU(JP))) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Upstream elevation for turning pump off [EOFFPU=',EOFFPU(JP),'] < grid top elevation for pump ',JP
    END IF
    IF (EOFFPU(JP) < EL(KB(IUPU(JP))+1,IUPU(JP))) THEN
      CALL ERRORS
      WRITE (ERR,FMTF2I) 'Upstream elevation for turning pump off [EOFFPU=',EOFFPU(JP),'] > bottom active layer [KB=',KB(IUPU(JP)),&
                         '] for pump ',JP
    END IF
    IF (KTPU(JP) < 2) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Top selective withdrawal layer [KTPU=',KTPU(JP),'] < layer 2 for pump ',JP
    END IF
    !IF (KBPU(JP) > KB(IUPU(JP))+1) THEN
    IF (KBPU(JP) > KB(IUPU(JP))) THEN                    ! cb 9/5/14
      CALL ERRORS
      WRITE (ERR,FMTI) 'Bottom selective withdrawal layer [KBPU=',KBPU(JP),'] > bottom active layer [KB=',KB(IUPU(JP)),'] for pump ',JP
    END IF
    IF (KBPU(JP) < KB(IUPU(JP))) THEN                    ! cb 9/5/14
      CALL WARNINGS
      WRITE (WRN,FMTI) 'Bottom selective withdrawal layer [KBPU=',KBPU(JP),'] < bottom active layer [KB=',KB(IUPU(JP)),'] for pump ',JP
    END IF
    IF (KTPU(JP) > KBPU(JP)) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Top selective withdrawal layer [KTPU=',KTPU(JP),'] > bottom selective withdrawal layer [KBPU=',KBPU(JP),       &
                       '] for pump ',JP
    END IF
    IF (EPU(JP) > EONPU(JP)) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Centerline elevation [EPU=',EPU(JP),'] > on elevation for pump ',JP
    END IF
    IF (EPU(JP) > EOFFPU(JP)) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Centerline elevation [EPU=',EPU(JP),'] > off elevation for pump ',JP
    END IF
    IF (PPUC(JP) == ' SPECIFY' .AND. IDPU(JP) /= 0) THEN                                                               !TC 09/30/03
      IF (PPUC(JP) /= ' SPECIFY' .AND. PPUC(JP) /= '   DISTR' .AND. PPUC(JP) /= ' DENSITY') THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Inflow placement [PPPUC='//PPUC(JP)//'] /= " SPECIFY" or "   DISTR" or " DENSITY" for pump ',JP
      END IF
      IF (ETPU(JP) > EL(2,IDPU(JP))) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Downstream inflow placement top elevation [ETPU=',ETPU(JP),'] > grid top elevation for pump ',JP
      END IF
      IF (ETPU(JP) < EL(KB(IDPU(JP))+1,IDPU(JP))) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Downstream inflow placement top elevation [ETPU=',ETPU(JP),'] < bottom active cell elevation for pump ',&
                           JP
      END IF
      IF (EBPU(JP) > EL(2,IDPU(JP))) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Downstream inflow placement bottom elevation [EBPU=',EBPU(JP),'] > grid top elevation for pump ',JP
      END IF
      IF (EBPU(JP) < EL(KB(IDPU(JP))+1,IDPU(JP))) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Downstream inflow placement bottom elevation [EBPU=',EBPU(JP),'] < bottom active cell elevation for '// &
                          'pump ',JP
      END IF
      !IF (EBPU(JP) > EPU(JP)) THEN    ! Deleted SW 10/18/2017
      !  CALL ERRORS
      !  WRITE (ERR,FMTFI) 'Downstream inflow placement bottom elevation [EBPU=',EBPU(JP),'] > elevation of pump ',JP
      !END IF
      IF (ETPU(JP) < EBPU(JP)) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Downstream inflow placement top elevation [ETPU=',ETPU(JP),'] < inflow placement bottom elevation for'//&
                          ' pump ',JP
      END IF
    END IF
    IF ((LATPUC(JP) /= '    DOWN') .AND. (LATPUC(JP) /= '     LAT')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Pump lateral control [LATPUC='//LATPUC(JP)(4:8)//'] must be either "DOWN" or "LAT" for pump ',JP
    END IF

     IF (LATPUC(JP) == '    DOWN') THEN     ! SW 3/26/10
      iflag=0
      DO JB=1,NBR
        IF (IUPU(JP) == DS(JB)) THEN
          iflag=1
          exit
        END IF
      END DO
        IF(IFLAG == 0)THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Since PUMP is specified as LATPUC='//LATPUC(JP)(4:8)//', it must be at end of a branch or change this to LAT for pump ',JP
        ENDIF
    END IF


  END DO

! Internal weir

  WRITE (WIN,*) '    internal weir'
  DO JWR=1,NIW
    IFLAG = 0
    DO JW=1,NWB
      DO JB=BS(JW),BE(JW)
        IF (IWR(JWR) >= US(JB) .AND. IWR(JWR) <= DS(JB)) THEN
          KT    = KTWB(JW)
          IFLAG = 1
          EXIT
        END IF
      END DO
      IF (IFLAG == 1) EXIT
    END DO
    
    IF(EKTWR(JWR) > 0.1)THEN     ! FOR OLDER VERSIONS WHERE ktwr /= 0 IF KTWR = 0 THEN IT IS A FLOATING WEIR
        
        KTWR(JWR)=INT(EKTWR(JWR))
        KBWR(JWR)=INT(EKBWR(JWR))
        
    IF (KTWR(JWR) <= KT) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Internal weir surface layer [KTWR=',KTWR(JWR),'] < water surface layer [KT=',KT,'] for weir ',JWR
    END IF
    IF (KTWR(JWR) > KBWR(JWR)) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Internal weir top layer [KTWR=',KTWR(JWR),'] > weir bottom layer [KBWR=',KBWR(JWR), '] for weir ',JWR
    END IF
    IF (KBWR(JWR) > KB(IWR(JWR))+1) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Internal weir bottom layer [KBWR=',KBWR(JWR),'] > bottom active layer [KB=',KB(IWR(JWR)),'] for weir ',JWR
    END IF
    
    ELSE
        IF (EKTWR(JWR) /= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,'(A,F8.2,A,I3,A,I3)') 'Internal weir surface layer [EKTWR=',EKTWR(JWR),'] is not ==0.0 or equal to a layer value for weir ',JWR
       END IF
        IF (EKBWR(JWR) > 0.0) THEN
        CALL ERRORS
        WRITE (ERR,'(A,F8.2,A,I3,A,I3)') 'Internal weir bottom layer [EKBWR=',EKBWR(JWR),'] is not < 0. It must be if EKTWR is equal to 0.0 indicating floating skimmer weir for weir ',JWR
        END IF
         IF (abs(EKBWR(JWR)) < 1.0) THEN
        CALL warnings
        WRITE (WRN,'(A,F8.2,A,I3,A,I3)') 'Internal weir bottom layer [EKBWR=',EKBWR(JWR),'] is less than 1 m in depth for this floating curtain weir #:',JWR
       END IF
        
    ENDIF
    
    
    DO JB=1,NBR
      IF (IWR(JWR) == US(JB)-1) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Internal weir segment [IWR=',IWR(JWR),'] is an upstream boundary segment in branch ',JB
      ELSE IF (IWR(JWR) == DS(JB)+1) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Internal weir segment [IWR=',IWR(JWR),'] is a downstream boundary segment in branch ',JB
      ELSE IF (IWR(JWR) == DS(JB)) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Internal weir segment [IWR=',IWR(JWR),'] is a downstream segment in branch ',JB
      END IF
    END DO
    
  END DO

! Lateral withdrawals

  WRITE (WIN,*) '    lateral withdrawals'
  DO JW=1,NWD
    IF (EWD(JW) > EL(2,IWD(JW))) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Withdrawal elevation [EWD=',EWD(JW),'] > grid top elevation for withdrawal ',JW
    ELSE IF (EWD(JW) < EL(KB(IWD(JW))+1,IWD(JW))) THEN
      CALL ERRORS
      WRITE (ERR,FMT2FI) 'Withdrawal elevation [EWD=',EWD(JW),'] < bottom active layer elevation [EL=',EL(KB(IWD(JW))+1,IWD(JW)),  &
                          '] for withdrawal ',JW
    END IF
    DO JB=1,NBR
      IF (IWD(JW) == US(JB)-1 .OR. (JB == 1 .AND. IWD(JW) == 0)) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Withdrawal segment [IWD=',IWD(JW),'] is an upstream boundary segment in branch ',JB
      ELSE IF (IWD(JW) == DS(JB)+1) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Withdrawal segment [IWD=',IWD(JW),'] is a downstream boundary segment in branch ',JB
      END IF
    END DO
    IF (KTWD(JW) < 2) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Top selective withdrawal layer [KTWD=',KTWD(JW),'] < layer 2 for withdrawal ',JW
    END IF
    !IF (KBWD(JW) > KB(IWD(JW))+1) THEN
    IF (KBWD(JW) > KB(IWD(JW))) THEN              ! cb 9/5/14
      CALL ERRORS
      WRITE (ERR,FMTI) 'Bottom selective withdrawal layer [KBWD=',KBWD(JW),'] > bottom active layer [KB=',KB(IWD(JW)),             &
                       '] for withdrawal ',JW
    END IF
    IF (KBWD(JW) < KB(IWD(JW))) THEN              ! cb 9/5/14
      CALL WARNINGS
      WRITE (WRN,FMTI) 'Bottom selective withdrawal layer [KBWD=',KBWD(JW),'] < bottom active layer [KB=',KB(IWD(JW)),             &
                       '] for withdrawal ',JW
    END IF
    IF (KTWD(JW) > KBWD(JW)) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Top selective withdrawal layer [KTWD=',KTWD(JW),'] > bottom selective withdrawal layer [KBWD=',KBWD(JW),   &
                       '] for withdrawal ', JW
    END IF
  END DO

! Tributaries

  WRITE (WIN,*) '    tributaries'
  IF (TRIBUTARIES) THEN
    DO JT=1,NTR
      IF (PTRC(JT) /= ' SPECIFY' .AND. PTRC(JT) /= '   DISTR' .AND. PTRC(JT) /= ' DENSITY') THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Inflow placement [PTRC='//PTRC(JT)//'] /= " SPECIFY" or "   DISTR" or " DENSITY" for tributary ',JP
      END IF
      ACTIVE_SEGMENT = .FALSE.
      DO JB=1,NBR
        IF (ITR(JT) >= US(JB) .AND. ITR(JT) <= DS(JB)) ACTIVE_SEGMENT = .TRUE.
      END DO
      IF (.NOT.ACTIVE_SEGMENT) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Tributary segment [ITR(',JT,')=',ITR(JT),'] is a boundary segment or not in the active grid for tributary ',JT
      END IF
      IF (PTRC(JT) == ' SPECIFY') THEN
        IF (ETTR(JT) > EL(2,ITR(JT))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Inflow placement top elevation [ETTR=',ETTR(JT),'] > grid top elevation for tributary ',JT
        END IF
        IF (ETTR(JT) < EL(KB(ITR(JT))+1,ITR(JT))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Inflow placement top elevation [ETTR=',ETTR(JT),'] < bottom active cell elevation for tributary ',JT
        END IF
        IF (EBTR(JT) > EL(2,ITR(JT))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Inflow placement bottom elevation [EBTR=',EBTR(JT),'] > grid top elevation for tributary ',JT
        END IF
        IF (EBTR(JT) < EL(KB(ITR(JT))+1,ITR(JT))) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Inflow placement bottom elevation [EBTR=',EBTR(JT),'] < bottom active cell elevation for tributary ',JT
        END IF
        IF (ETTR(JT) < EBTR(JT)) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Inflow placement top elevation [ETTR=',ETTR(JT),'] < bottom elevation inflow placement for '//        &
                            'tributary ',JT
        END IF
      END IF
    END DO
  END IF

! Distributed tributaries

  WRITE (WIN,*) '    distributed tributaries'
  DO JB=1,NBR
    IF ((DTRC(JB) /= '      ON') .AND. (DTRC(JB) /= '     OFF')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Distributed tributary control [DTRC='//DTRC(JB)(6:8)//'] must be either " ON" or "OFF" for branch ',JB
    END IF
  END DO

! Hydrodynamic output

  WRITE (WIN,*) 'Output controls'
  WRITE (WIN,*) '  hydrodynamic variables'
  DO JW=1,NWB
    DO JH=1,NHY
      IF ((HPRWBC(JH,JW) /= '      ON') .AND. (HPRWBC(JH,JW) /= '     OFF')) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Hydrodynamic output control [HPRWBC='//HPRWBC(JH,JW)(6:8)//'] must be either " ON" or "OFF" for '//      &
                         'variable ',JH,' in waterbody ',JW
      END IF
    END DO
  END DO

! Snapshot output

  WRITE (WIN,*) '  snapshot'
  DO JW=1,NWB
    IF ((SNPC(JW) /= '      ON') .AND. (SNPC(JW) /= '     OFF')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Snapshot output: Snapshot control [SNPC='//SNPC(JW)(6:8)//'] must be either " ON" or "OFF" for waterbody ',JW
    END IF
    IF (SNPC(JW) == '      ON') THEN
      IF (NSNP(JW) <= 0) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Snapshot output: Number of snapshot dates [NSNP=',NSNP(JW),'] <= 0 for waterbody ',JW
      END IF
      DO J=1,NSNP(JW)-1
        IF (SNPD(J+1,JW) <= SNPD(J,JW)) THEN
          CALL ERRORS
          WRITE (ERR,FMT2FI) 'Snapshot output: Snapshot date [SNPD=',SNPD(J+1,JW),'] <= previous date [SNPD=',SNPD(J,JW),'] in waterbody ',JW
        END IF
      END DO
      DO J=1,NSNP(JW)
        IF (SNPF(J,JW) <= 0.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMT2FI) 'Snapshot output: Snapshot frequency [SNPF=',SNPF(J,JW),'] <= 0 for date [SNPD=',SNPD(J,JW),'] in waterbody ',JW
        END IF
      END DO
      IF (NISNP(JW) <= 0) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Snapshot output: Number of snapshot segment [NISNP=',NISNP(JW),'] s <= 0 for waterbody ',JW
      END IF
      DO J=1,NISNP(JW)
      ACTIVE_SEGMENT = .FALSE.
        DO JB=1,NBR
          IF (ISNP(J,JW) >= US(JB) .AND. ISNP(J,JW) <= DS(JB)) ACTIVE_SEGMENT = .TRUE.
        END DO
        IF (.NOT.ACTIVE_SEGMENT) THEN
          CALL ERRORS
          WRITE (ERR,FMTI) 'Snapshot output: Snapshot segment [ISNP=',ISNP(J,JW),'] is a boundary segment or not in the active grid in waterbody ',JW
        END IF
      END DO
    END IF
  END DO

! Screen output

  WRITE (WIN,*) '  screen'
  DO JW=1,NWB
    IF ((SCRC(JW) /= '      ON') .AND. (SCRC(JW) /= '     OFF')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Screen output: Screen control [SCRC='//SCRC(JW)(6:8)//'] must be either " ON" or "OFF" for waterbody ',JW
    END IF
    IF (SCRC(JW) == '      ON') THEN
      IF (NSCR(JW) <= 0) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Screen output: Number of screen dates [NSCR=',NSCR(JW),'] <= 0 for waterbody ',JW
      END IF
      DO J=1,NSCR(JW)-1
        IF (SCRD(J+1,JW) <= SCRD(J,JW)) THEN
          CALL ERRORS
          WRITE (ERR,FMT2FI) 'Screen output: Screen date [SCRD=',SCRD(J+1,JW),'] <= previous date [SCRD=',SCRD(J,JW),'] in waterbody ',JW
        END IF
      END DO
      DO J=1,NSCR(JW)
        IF (SCRF(J,JW) <= 0.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Screen output: Screen frequency [SCRF=',SCRF(J,JW),'] <= 0 for waterbody ',JW
        END IF
      END DO
    END IF
  END DO

! Time series output

  WRITE (WIN,*) '  time series'
  IF ((TSRC /= '      ON') .AND. (TSRC /= '     OFF')) THEN
    CALL ERRORS
    WRITE (ERR,FMTA) 'Time series output: Time series control [TSRC='//TSRC(6:8)//'] must be either " ON" or "OFF"'
  END IF
  IF (TSRC == '      ON') THEN
    IF (NTSR <= 0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Time series output: Number of time series dates [NTSR=',NTSR,' <= 0. Since TSRC=ON, it must be >=1.'
    END IF
    IF (NIKTSR <= 0) THEN
      CALL ERRORS
      WRITE (WRN,FMTI) 'Time series output: Number of time series segments [NITSR=',NIKTSR,' <= 0. Since TSRC=ON, it must be >=1.'
    END IF
    DO J=1,NTSR-1
      IF (TSRD(J+1) <= TSRD(J)) THEN
        CALL ERRORS
        WRITE (ERR,FMTF) 'Time series output: Time series date [TSRD=',TSRD(J+1),'] <= previous date [TSRD=',TSRD(J),']'
      END IF
    END DO
    DO J=1,NTSR
      IF (TSRF(J) <= 0.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTF) 'Time series output: Time series frequency [TSRF=',TSRF(J),'] <= 0 for date [TSRD=',TSRD(J),']'
      END IF
    END DO
    DO JB=1,NBR                 ! SW 8/12/05
    DO J=1,NIKTSR                 ! SW 7/5/05
      IF (ITSR(J) == (US(JB)-1) .OR. ITSR(J) == (DS(JB)+1) .OR. ITSR(J) < 2 .OR. ITSR(J) > IMX) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Time series output: Time series segment[ITSR=',ITSR(J),'] is located at a boundary segment or outside the grid for NITSR ',J
      END IF
    END DO
    END DO
  END IF

! Spreadsheet output

  WRITE (WIN,*) '  spreadsheet'
  DO JW=1,NWB
    IF ((SPRC(JW) /= '      ON') .AND. (SPRC(JW) /= '     OFF') .AND.SPRC(JW) /= '     ONV' ) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Spreadsheet output: Spreadsheet control [SPRC='//SPRC(JW)(6:8)//'] must be either " ON" or "ONV" or "OFF" for waterbody ',JW
    END IF
    IF (SPRC(JW) == '      ON' .OR. SPRC(JW) == '     ONV') THEN
      IF (NSPR(JW) <= 0) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Spreadsheet output: Number of spreadsheet dates [NSPR=',NSPR(JW),']  <= 0 for waterbody ',JW
      END IF
      DO J=1,NSPR(JW)-1
        IF (SPRD(J+1,JW) <= SPRD(J,JW)) THEN
          CALL ERRORS
          WRITE (ERR,FMT2FI) 'Spreadsheet output: Spreadsheet date [SPRD=',SPRD(J+1,JW),'] <= previous date [SPRD=', SPRD(J,JW),'] in waterbody ',JW
        END IF
      END DO
      DO J=1,NSPR(JW)
        IF (SPRF(J,JW) <= 0.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMT2FI) 'Spreadsheet output: Spreadsheet frequency [SPRF=',SPRF(J,JW),'] <= 0 for date [SPRD=',SPRD(J,JW),'] in waterbody ',JW
        END IF
      END DO
      IF (NISPR(JW) <= 0) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Spreadsheet output: Number of spreadsheet segments [NISPR=',NISPR(JW),'] <= 0 for waterbody ',JW
      END IF
      DO J=1,NISPR(JW)
      ACTIVE_SEGMENT = .FALSE.
        DO JB=BS(JW),BE(JW)     ! 1,NBR
          IF (ISPR(J,JW) >= US(JB) .AND. ISPR(J,JW) <= DS(JB)) ACTIVE_SEGMENT = .TRUE.
        END DO
        IF (.NOT.ACTIVE_SEGMENT) THEN
          CALL ERRORS
          WRITE (ERR,FMTI) 'Spreadsheet output: Spreadsheet segment [ISPR=',ISPR(J,JW),'] is a boundary segment or not in the active grid in waterbody ',JW
        END IF
      END DO
    END IF
  END DO

! Profile output

  WRITE (WIN,*) '  profile plots'
  DO JW=1,NWB
    IF ((PRFC(JW) /= '      ON') .AND. (PRFC(JW) /= '     OFF')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Profile output: Profile control [PRFC='//PRFC(JW)(6:8)//'] must be either " ON" or "OFF" for waterbody ',JW
    END IF
    IF (PRFC(JW) == '      ON') THEN
      IF (NPRF(JW) <= 0) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Profile output: Number of profile dates [NPRF=',NPRF(JW),'] <= 0 for waterbody ',JW
      END IF
      DO J=1,NPRF(JW)-1
        IF (PRFD(J+1,JW) <= PRFD(J,JW)) THEN
          CALL ERRORS
        WRITE (ERR,FMT2FI) 'Profile output: Profile date [PRFD=',PRFD(J+1,JW),'] <= previous date [PRFD=',PRFD(J,JW),'] in waterbody ',JW
        END IF
      END DO
      DO J=1,NPRF(JW)
        IF (PRFF(J,JW) <= 0.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMT2FI) 'Profile output: Profile frequency [PRFF=',PRFF(J,JW),'] <= 0 for date [PRFD=',PRFD(J,JW),'] in waterbody ',JW
        END IF
      END DO
      IF (NIPRF(JW) <= 0) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Profile output: Number of profile segments [NIPRF=',NIPRF(JW),'] <= 0 for waterbody ',JW
      END IF
      DO J=1,NIPRF(JW)
    IF(iprf(1,1) /= -1)then  ! SW 4/1/16
      ACTIVE_SEGMENT = .FALSE.
        DO JB=BS(JW),BE(JW)     ! 1,NBR
          IF (IPRF(J,JW) >= US(JB) .AND. IPRF(J,JW) <= DS(JB)) ACTIVE_SEGMENT = .TRUE.
        END DO
        IF (.NOT.ACTIVE_SEGMENT) THEN
          CALL ERRORS
          WRITE (ERR,FMTI) 'Profile output: Profile segment [IPRF=',IPRF(J,JW),'] is a boundary segment or not in the active grid in waterbody ',JW
        END IF
    ENDIF
    
      END DO
    END IF
  END DO

! Vector output

  WRITE (WIN,*) '  vector plots'
  DO JW=1,NWB
    IF ((VPLC(JW) /= '      ON') .AND. (VPLC(JW) /= '     OFF')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'VPL output: Vector plot control [VPLC='//VPLC(JW)(6:8),'] must be either " ON" or "OFF" for waterbody ',JW
    END IF
    IF (VPLC(JW) == '      ON') THEN
      IF (NVPL(JW) <= 0) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'VPL output: VPL is ON, but number of vector plot dates [NVPL=',NVPL(JW),'] <= 0 for waterbody ',JW
      END IF
      DO J=1,NVPL(JW)-1
        IF (VPLD(J+1,JW) <= VPLD(J,JW)) THEN
          CALL ERRORS
          WRITE (ERR,FMT2FI) 'VPL output: Vector plot date [VPLD=',VPLD(J,JW),'] <= previous date [VPLD=',VPLD(J-1,JW),'] in waterbody ',JW
        END IF
      END DO
      DO J=1,NVPL(JW)
        IF (VPLF(J,JW) <= 0.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMT2FI) 'Vector output: Vector plot frequency [VPLF=',VPLF(J,JW),'] <= 0 for date [VPLD=',VPLD(J,JW),'] in waterbody ',JW
        END IF
      END DO
    END IF
  END DO

! Contour output

  WRITE (WIN,*) '  contour plots'
  DO JW=1,NWB
    IF ((CPLC(JW) /= '      ON') .AND. (CPLC(JW) /= '     OFF')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Contour output: Contour plot control [VPLC='//VPLC(JW)(6:8)//'] must be either " ON" or "OFF" for waterbody ',JW
    END IF
    IF ((TCPL(JW) /= '      ON') .AND. (TCPL(JW) /= '     OFF')) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Contour output: Contour plot control [TECPLOT='//TCPL(JW)(6:8)//'] must be either " ON" or "OFF" for waterbody ',JW
    END IF
    IF (CPLC(JW) == '      ON') THEN
      IF (NCPL(JW) <= 0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTI) 'Contour output: Number of contour plot dates [NCPL=',NCPL(JW),'] <= 0 for waterbody ',JW
      END IF
      DO J=1,NCPL(JW)-1
        IF (CPLD(J+1,JW) <= CPLD(J,JW)) THEN
          CALL ERRORS
          WRITE (ERR,FMT2FI) 'Contour output: Contour plot date [VPLD=',VPLD(J,JW),'] <= previous date [VPLD=',VPLD(J-1,JW),'] in waterbody ',JW
        END IF
      END DO
      DO J=1,NCPL(JW)
        IF (CPLF(J,JW) <= 0.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMT2FI) 'Contour output: Contour plot frequency [VPLF=',VPLF(J,JW),'] <= 0 for date [VPLD=',VPLD(J,JW),'] in waterbody ',JW
        END IF
      END DO
    END IF
  END DO
  
  TECPLOTBR=.FALSE.
  DO JW=1,NWB
  IF(TCPL(JW)=='      ON')THEN         ! CHECK FOR W2_TECPLOTBR.CSV FILE
    INQUIRE(FILE='w2_tecplotbr.csv',EXIST=TECPLOTBR) 
    WRITE (WIN,*) '  tecplot contour plots'
      IF(TECPLOTBR)THEN
      WRITE (WIN,*) '   reading w2_tecplotbr.csv'    
      OPEN(CON,FILE='w2_tecplotbr.csv',STATUS='OLD')
        READ(CON,*)
        READ(CON,*)JTEC    ! NUMBER OF BRANCHES TO PLOT, MUST BE LESS THAN NBR
        READ(CON,*)
        READ(CON,*)(IBR(J),J=1,JTEC)
         DO J=1,JTEC
          IF(IBR(J) < 1 .OR. IBR(J) > NBR)THEN
          CALL ERRORS
          WRITE (ERR,'(A,I4,A,I4)') 'W2_tecplotbr.csv: Branch number for plotting must be between 1 and NBR. Branch #:',IBR(J),' for #:',J
          ENDIF
         ENDDO
        CLOSE(CON)
      ENDIF
  ENDIF
  ENDDO
  

! Withdrawal output

  WRITE (WIN,*) '  withdrawal'
  IF ((WDOC /= '      ON') .AND. (WDOC /= '     OFF') .and. (WDOC /= '     ONH') .and. (WDOC /= '     ONS')) THEN
    CALL ERRORS
    WRITE (ERR,FMTA) 'Withdrawal outflow control [WDOC='//WDOC(6:8)//'] must be either ON or ONH or ONS or OFF'
  END IF
  IF ((WDOC == '      ON') .or. (WDOC == '     ONS') .or. (WDOC == '     ONH')) THEN
        IF(NIWDO <= 0)THEN
        CALL ERRORS
        WRITE(ERR,*)'WITHDRAWAL OUTPUT:NIWDO <=0. It must be 1 or more since WDOC=ON or ONH or ONS.'
        endif
        IF(NWDO <= 0)THEN
        CALL ERRORS
        WRITE(ERR,*)'WITHDRAWAL OUTPUT:NWDO <=0. It must be 1 or more since WDOC=ON or ONH or ONS.'
        endif
    DO J=1,NIWDO
      ACTIVE_SEGMENT = .FALSE.
      DO JB=1,NBR
        IF (IWDO(J) >= US(JB) .AND. IWDO(J) <= DS(JB)) ACTIVE_SEGMENT = .TRUE.
      END DO
      IF (.NOT.ACTIVE_SEGMENT) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Withdrawal out segment [IWDO=',IWDO(J),'] is not an active segment'
      END IF
    END DO
    DO J=1,NWDO
      IF (WDOF(J) <= 0.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTF) 'Withdrawal out frequency [WDOF=',WDOF(J),'] <= 0'
      END IF
    END DO
  END IF

! Restart

  WRITE (WIN,*) '  restart'
  IF ((RSOC /= '      ON') .AND. (RSOC /= '     OFF')) THEN
    CALL ERRORS
    WRITE (ERR,FMTA) 'Restart output control [RSOC='//RSOC(6:8)//'] must be either " ON" or "OFF"'
  END IF
  IF (RSOC == '      ON') THEN
    IF (NRSO <= 0) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Number of restart output dates [NRSO=',NRSO,'] <= 0'
    END IF
    DO J=1,NRSO-1
      IF (RSOD(J+1) <= RSOD(J)) THEN
        CALL ERRORS
        WRITE (ERR,FMTF) 'Restart output date [RSOD=',RSOD(J+1),'] <= previous date [RSOD=',RSOD(J),']'
      END IF
    END DO
    DO J=1,NRSO
      IF (RSOF(J) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTF) 'Restart output frequency [RSOF=',RSOF(J),'] <= 0'
      END IF
    END DO
  END IF

! Constituent controls

  WRITE (WIN,*) 'Constituents'
  WRITE (WIN,*) '  active constituents'
  IF ((CCC /= '      ON') .AND. (CCC /= '     OFF')) THEN
    CALL ERRORS
    WRITE (ERR,FMTA) 'Constituent computation control [CCC='//CCC(6:8)//'] must be either " ON" or "OFF"'
  END IF
  IF (CONSTITUENTS) THEN
    DO JW=1,NWB
      IF ((SEDC(JW) /= '      ON') .AND. (SEDC(JW) /= '     OFF')) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Sediment computation control [SEDC='//SEDC(JW)(6:8)//'] must be either " ON" or "OFF" for waterbody ',JW
      END IF
      IF ((PRNSC(JW) /= '      ON') .AND. (PRNSC(JW) /= '     OFF')) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Sediment print control [PRNSC='//PRNSC(JW)(6:8)//'] must be either " ON" or "OFF" for waterbody ',JW
      END IF
    END DO
    IF (CUF > 12) THEN
      CALL WARNINGS
      WRITE (WRN,FMTI) 'Kinetic update frequency [CUF=',CUF,'] > 12'
    END IF
    IF (CUF < 1) THEN                                                      ! SW 5/19/10
      CALL ERRORS
      WRITE (ERR,FMTI) 'Kinetic update frequency [CUF=',CUF,'] < 1. CUF must be >= 1.'
    END IF

!** Active constituents

    DO JAC=1,NAC
      JC = CN(JAC)
      IF ((CAC(JC) /= '      ON') .AND. (CAC(JC) /= '     OFF')) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Active constituent computation control [CAC='//CAC(JC)(6:8)//'] must be either " ON" or "OFF" '//        &
                         'for constituent ',JC
      END IF
    END DO

!** pH calculations

    DO JW=1,NWB
      IF (CDWBC(PH_DER,JW) == '      ON') THEN
        IF (CAC(NTIC) == '     OFF' .OR. CAC(NALK) == '     OFF') THEN
          CALL ERRORS
          WRITE (ERR,FMTA) 'pH computations are =" ON", so TIC and alkalanity must also be turned " ON"'
        END IF
      END IF
    END DO
    
    
 !** TURBIDITY calculations

    DO JW=1,NWB
      IF (CDWBC(TURB_DER,JW) == '      ON') THEN
        IF (CDWBC(TOTSS_DER,JW) == '     OFF') THEN
          CALL ERRORS
          WRITE (ERR,FMTA) 'Turbidity computations are =" ON", so TSS must also be turned " ON" as a derived constituent'
        END IF
      END IF
    END DO


!** Initial constituent concentrations

    WRITE (WIN,*) '  initial concentrations'
    DO JW=1,NWB
      DO JAC=1,NAC
        JC = CN(JAC)
        IF (C2IWB(JC,JW) < -2.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTI) 'Initial constituent concentration < -2 for constituent ',JC,' in waterbody ',JW
        END IF
        IF (CAC(JC) == '      ON') THEN
          IF (C2IWB(JC,JW) == 0.0) THEN
            NAME = TRIM(CNAME1(JC))
        !    IF (NAME /= 'Residence time' .AND. NAME /= 'Water age') THEN
            IF (NAME /= 'WaterAge' .AND. NAME /= 'Residence Time') THEN     ! SW 7/15/14
              CALL WARNINGS
              WRITE (WRN,FMTI) 'Initial constituent concentration = 0 for constituent ',JC,' in waterbody ',JW
            END IF
          END IF
        END IF
      END DO
    END DO

!** Constituent output

    WRITE (WIN,*) '  constituent output'
    DO JW=1,NWB
      DO JC=1,NCT     !JAC=1,NAC
       ! JC = CN(JAC)
        IF (CPRWBC(JC,JW) /= '      ON' .AND. CPRWBC(JC,JW) /= '     OFF') THEN
          CALL ERRORS
          WRITE (ERR,FMTI) 'Constituent print control [CPRWBC='//CPRWBC(JC,JW)(6:8)//'] must be either " ON" or "OFF" for '//      &
                           'constituent ',JC
        END IF
        IF (CPRWBC(JC,JW) == '      ON') THEN
          IF (CAC(JC) /= '      ON') THEN
            CALL WARNINGS
            WRITE (WRN,FMTI) 'Constituent print control [CPRWBC='//CPRWBC(JC,JW)(6:8)//'] is " ON" but active constituent'//       &
                             ' control [CAC='//CAC(JC)(6:8)//'] is not " ON" for constituent ',JC
          END IF
        END IF
      END DO
    END DO

!** Inflow constituent concentrations

    WRITE (WIN,*) '  branch inflow concentrations'
    DO JB=1,NBR
      DO JC=1,NCT        !JAC=1,NAC           !JC = CN(JAC)
        IF ((CINBRC(JC,JB) /= '      ON') .AND. (CINBRC(JC,JB) /= '     OFF')) THEN
          CALL ERRORS
          WRITE (ERR,FMTI) 'Inflow active constituent control [CINBRC='//CINBRC(JC,JB)(6:8)//'] must be either " ON" or "OFF" '//  &
                           'for constituent ',JC,' for branch ',JB
        END IF
        IF (CAC(JC) == '      ON') THEN
          IF (CINBRC(JC,JB) /= '      ON') THEN
            NAME = TRIM(CNAME1(JC))
            IF (NAME /= 'Residence time' .AND. NAME /= 'Water age') THEN
              CALL WARNINGS
              WRITE (WRN,FMTI) 'Constituent ',JC,' is " ON", but inflow constituent concentration is "OFF" for branch ',JB
            END IF
          END IF
        END IF
      END DO
    END DO

!** Tributary constituent concentrations

    WRITE (WIN,*) '  tributary inflow concentrations'
    IF (TRIBUTARIES) THEN
      DO JT=1,NTR
        DO JC=1,NCT     !JAC=1,NAC        !JC = CN(JAC)
          IF ((CINTRC(JC,JT) /= '      ON') .AND. (CINTRC(JC,JT) /= '     OFF')) THEN
            CALL ERRORS
            WRITE (ERR,FMTI) 'Inflow file constituent control [CINTRC='//CINTRC(JC,JT)(6:8)//'] must be either " ON" or "OFF" '//&
                             'for constituent ',JC,' for tributary ',JT
          END IF
          IF (CAC(JC) == '      ON') THEN
            IF (CINTRC(JC,JT) /= '      ON') THEN
              NAME = TRIM(CNAME1(JC))
              IF (NAME /= 'WaterAge' .AND. NAME /= 'Water Age') THEN
                CALL WARNINGS
                WRITE (WRN,FMTI) 'Constituent ',JC,' is " ON", but inflow constituent concentration is "OFF" for tributary ',JT
              END IF
            END IF
          END IF
        END DO
      END DO
    END IF

!** Distributed tributary concentrations

    WRITE (WIN,*) '  distributed tributary inflow concentrations'
    DO JB=1,NBR
      IF (DIST_TRIB(JB)) THEN
        DO jc=1,nct     !JAC=1,NAC
         ! JC = CN(JAC)
          IF ((CDTBRC(JC,JB) /= '      ON') .AND. (CDTBRC(JC,JB) /= '     OFF')) THEN
            CALL ERRORS
            WRITE (ERR,FMTI) 'Distributed tributary inflow file constituent control [CDTBRC='//CDTBRC(JC,JB)(6:8)//'] must be '//&
                             'either " ON" or "OFF" for constituent ',JC,' for branch ',JB
          END IF
          IF (CAC(JC) == '      ON') THEN
            IF (CDTBRC(JC,JB) /= '      ON') THEN
              CALL WARNINGS
              WRITE (WRN,FMTI) 'Constituent ',JC,' is " ON", but distributed tributary inflow constituent concentration '//        &
                               'is "OFF" for branch ',JB
            END IF
          END IF
        END DO
      END IF
    END DO

!** Precipitation concentrations

    WRITE (WIN,*) '  precipitation concentrations'
    DO JW=1,NWB
      IF (PRECIPITATION(JW)) THEN
        DO JB=BS(JW),BE(JW)
          DO jc=1,NCT     !JAC=1,NAC
            !JC = CN(JAC)
            IF ((CPRBRC(JC,JB) /= '      ON') .AND. (CPRBRC(JC,JB) /= '     OFF')) THEN
              CALL ERRORS
              WRITE (ERR,FMTI) 'Precipitation active constituent control [CINBRC='//CINBRC(JC,JB)(6:8)//'] must be either " ON" '//&
                               'or "OFF" for constituent ',JC,' for branch ',JB
            END IF
            IF (CAC(JC) == '      ON') THEN
              IF (CPRBRC(JC,JB) /= '      ON') THEN
                CALL WARNINGS
                WRITE (WRN,FMTI) 'Constituent ',JC,' is " ON", but precipitation constituent concentration is "OFF" for  branch ',JB
              END IF
            END IF
          END DO
        END DO
      END IF
    END DO
  ELSEIF (CUF < 1) THEN                                                      ! SW 5/19/10
      CALL ERRORS
      WRITE (ERR,FMTI) 'Kinetic update frequency [CUF=',CUF,'] < 1. CUF must be >= 1 even when CONSTITUENTS=OFF.'
  END IF

! Light extinction

  WRITE (WIN,*) '  kinetic coefficients'
  WRITE (WIN,*) '    light extinction'
  DO JW=1,NWB
    IF (EXH2O(JW) > 0.95) THEN
      CALL WARNINGS
      WRITE (WRN,FMTI) 'Pure water light extinction > 0.95 for waterbody ',JW
    ELSE IF (EXH2O(JW) < 0.25) THEN
      CALL ERRORS
      WRITE (ERR,FMTI) 'Pure water light extinction < 0.25 for waterbody ',JW
    END IF
    IF (BETA(JW) > 0.75) THEN
      CALL WARNINGS
      WRITE (WRN,FMTFI) 'Short wave solar radiation absorption [BETA=',BETA(JW),'] > 0.75 for waterbody ',JW
    ELSE IF (BETA(JW) <= 0.3 .AND. BETA(JW) >= 0.1) THEN
      CALL WARNINGS
      WRITE (WRN,FMTFI) 'Short wave solar radiation absorption [BETA=',BETA(JW),'] < 0.3 for waterbody ', JW
    ELSE IF (BETA(JW) < 0.1) THEN
      CALL ERRORS
      WRITE (ERR,FMTFI) 'Short wave solar radiation absorption [BETA=',BETA(JW),'] < 0.1 for waterbody ', JW
    END IF
    IF ((EXC(JW) /= '      ON') .AND. (EXC(JW)   /= '     OFF')) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Light Extinction File Read Control [EXC='//EXC(JW)//'] must be either " ON" '//&
                               'or "OFF" for waterbody ',JW
    END IF
        IF ((EXIC(JW) /= '      ON') .AND. (EXIC(JW)   /= '     OFF')) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Light Extinction Interpolation Control [EXIC='//EXIC(JW)//'] must be either " ON" '//&
                               'or "OFF" for waterbody ',JW
    END IF
    IF (CONSTITUENTS) THEN
      IF (SUSPENDED_SOLIDS) THEN
        IF (EXSS(JW) > 0.1) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Suspended solids light extinction [EXSS=',EXSS(JW),'] > 0.1 for waterbody ', JW
        ELSE IF (EXSS(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Suspended solids light extinction [EXSS=',EXSS(JW),'] < 0 for waterbody ',   JW
        END IF
      END IF
      IF (LABILE_DOM .OR. REFRACTORY_DOM .OR. LABILE_POM .OR. REFRACTORY_POM) THEN
        IF (EXOM(JW) > 0.85) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Organic matter light extinction [EXOM=',EXOM(JW),'] > 0.85 for waterbody ',  JW
        ELSE IF (EXOM(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Organic matter light extinction [EXOM=',EXOM(JW),'] < 0 for waterbody ',     JW
        END IF
      END IF
    END IF
    IF(EXC(JW) == '      ON')THEN

        ! Light Extinction File

  WRITE (WIN,*) '  light extinction file'
  OPEN  (NPT,FILE=EXTFN(JW),IOSTAT=IERR,STATUS='OLD')                                                                      !TC 10/22/02
  WRITE(WIN,*)'       reading light ext file:'//extfn(jw)
  IF (IERR == 0) THEN
      AID=EXTFN(JW)

      READ(NPT,'(A1)',ERR=400)ICHAR1
      IF(ICHAR1=='$')WRITE(WIN,*) '     csv format for EXT file'
      READ (NPT,'(/)')
      
     ! CHECK IF ERROR -CSV FORMAT CHECK
        IF(ICHAR1 /= '$')THEN
        READ(NPT,'(A72)')TITLE(1)
            DO JJ=1,50    ! CHECK ONLY FIRST 50 CHARACTERS
                IF(TITLE(1)(JJ:JJ) == ",")THEN
                    CALL ERRORS
                    WRITE (ERR,'(A175)') 'Light Extinction file: File not designated as a CSV file (no "$" as first character) but commas found between input values for '//EXTFN(JW)
                    EXIT
                ENDIF
            ENDDO
        BACKSPACE(NPT)
        ENDIF       
        
        J       = 1
        JDAYO   = 0.0
        
          DO WHILE (.TRUE.)
            IF(ICHAR1=='$')THEN
            READ (NPT,*,END=335,ERR=400) JDAY,EXH2O(JW)
                  ELSE
            READ (NPT,'(10F8.0)',END=335,ERR=400) JDAY,EXH2O(JW)
            ENDIF

            IF ((J == 1) .AND. (JDAY > TMSTRT)) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'Starting date [',JDAY,'] > simulation start date [TMSTRT=',TMSTRT,'] in '//EXTFN(JW)
            ELSE IF (JDAY <= JDAYO .AND. J /= 1) THEN
                if(JDAY /= 0.0)then
                CALL ERRORS
              WRITE (ERR,FMTF) 'Julian date ',JDAY,' <= previous date ',JDAYO,' in '//EXTFN(JW)
                else
                CALL ERRORS
              WRITE (ERR,FMTF) '[**Note: These errors could be a result of blank lines at the end of the file**]: Julian date ',JDAY,' <= previous date ',JDAYO,' in '//EXTFN(JW)
                endif
            END IF

            IF (EXH2O(JW) < 0.15) THEN
              CALL ERRORS
              WRITE (ERR,FMTF) 'EXH2O(JW) cannot usually be less than 0.15 m-1:',EXH2O(JW),' at date [JDAY=',JDAY,'] in '//EXTFN(JW)
            ELSEIF(EXH2O(JW) < 0.25) THEN
              CALL WARNINGS
              WRITE (WRN,FMTF) 'EXH2O(JW) is less than 0.25 m-1:',EXH2O(JW),' at date [JDAY=',JDAY,'] in '//EXTFN(JW) 
            END IF

            
            JDAYO = JDAY
            J     = J+1
          END DO
335      CONTINUE
          IF (JDAY < TMEND) THEN
            CALL ERRORS
            WRITE (ERR,FMTF) 'Ending time ',JDAY,' < ending simulation time [TMEND=',TMEND,'] in '//EXTFN(JW)
          END IF
      
  ELSE                                                                                                                 !TC 10/22/02
    CALL ERRORS                                                                                                        !TC 10/22/02
    WRITE (ERR,FMTA) 'Could not open light extinction filename '//EXTFN(JW)                                                         !TC 10/22/02
  END IF                                                                                                               !TC 10/22/02
  CLOSE (NPT)                                                                                                          !SW 04/03/02
        
        
    ENDIF
  END DO
  IF (ALGAE) THEN
    DO JA=1,NAAL
      IF (EXA(ALAC(JA)) > 0.8) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal light extinction [EXA=',EXA(ALAC(JA)),'] > 0.8 for algal group ',JW
      ELSE IF (EXA(ALAC(JA)) < 0.1) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal light extinction [EXA=',EXA(ALAC(JA)),'] < 0.1 for algal group ',JW
      END IF
    END DO
  END IF

! Constituents

  IF (CONSTITUENTS) THEN

!** Generic constituent

    WRITE (WIN,*) '    generic'
    IF (GENERIC_CONST) THEN
      DO JG=1,NGC
        IF (CGQ10(JG) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Temperature multiplier [CGQ10]=',CGQ10(JG),'] < 0.0 for generic constituent ', JG
        ELSE IF (CGQ10(JG) > 1.5) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Temperature multiplier [CGQ10]=',CGQ10(JG),'] > 1.5 for generic constituent ', JG
        END IF
        IF (CG0DK(JG) < 0.0 .AND. CG0DK(JG) /= -1.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) '0-order decay rate [CG0DK]=',    CG0DK(JG),'] < 0.0 for generic constituent ', JG
        ELSE IF (CG0DK(JG) > 10.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) '0-order decay rate [CG0DK]=',    CG0DK(JG),'] > 10.0 for generic constituent ',JG
        END IF
        IF (CGS(JG) < 0.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Settling rate [CGS]=',           CGS(JG),'] < 0.0. Hence, it will be bouyant. Generic constituent: ',   JG
        ELSE IF (CGS(JG) > 10.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Settling rate [CGS]=',           CGS(JG),'] > 10.0 for generic constituent ',  JG
        END IF
        IF (CGLDK(JG) < 0.0) THEN   
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Photodegradation parameter [CGLDK]=',           CGLDK(JG),'] < 0.0 for generic constituent ',   JG
        ELSE IF (CGLDK(JG) > 1.0E-4) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Photodegradation parameter [CGLDK]=',           CGLDK(JG),'] > 1.0E-4 m2/J for generic constituent ',  JG
        END IF
        IF (CGKLF(JG) < 0.0) THEN   
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Volatilization reaeration rate fraction of O2 rate [CGKLF]=',           CGKLF(JG),'] < 0.0 for generic constituent ',   JG
        ELSE IF (CGKLF(JG) > 2.5) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Volatilization reaeration rate fraction of O2 rate [CGKLF]=',           CGKLF(JG),'] > 2.5 for generic constituent ',  JG
        END IF
        IF (CGCS(JG) == -1.00) THEN   !CGLDK(JG),CGKLF(JG),CGCS(JG)
          WRITE(WIN,*)'TDG IS ON'
          IF(CGKLF(JG) /= 1.034)THEN
                CALL WARNINGS
                WRITE (WRN,FMTFI) 'N2 volatilization reaeration rate fraction of O2 rate [CGKLF]=',           CGKLF(JG),'] /= 1.034 for generic constituent ',  JG
          ENDIF
        ENDIF
        IF (CGCS(JG) > 100.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Saturation concentration [CGCS]=',           CGCS(JG),'] > 100.0 for generic constituent ',   JG
        ELSE IF (CGCS(JG) > 50.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Saturation concentration [CGCS]=',           CGCS(JG),'] > 50.0 for generic constituent ',   JG
        END IF
      END DO
    END IF

!** Suspended solids

    WRITE (WIN,*) '    suspended solids'
    IF (SUSPENDED_SOLIDS) THEN
      DO JS=1,NASS
        IF (SSS(SSAC(JS)) > 10.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Suspended solids settling rate [SSS=',SSS(JS),'] > 10.0 for suspended solids group ',JS
        ELSE IF (SSS(SSAC(JS)) == 0.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Suspended solids settling rate [SSS=',SSS(JS),']  = 0.0 for suspended solids group ',JS
        ELSE IF (SSS(SSAC(JS)) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Suspended solids settling rate [SSS=',SSS(JS),']  < 0.0 for suspended solids group ',JS
        END IF
        IF (SEDRC(JS) /= '      ON' .AND. SEDRC(JS) /= '     OFF') THEN                        ! SW 1/16/04
          CALL ERRORS
          WRITE (ERR,FMTI) 'Wind sediment resuspension control   [SEDRC='//SEDRC(JS)//'] /= " ON" or "OFF" for SS group ',JS
        END IF

	  END DO
    END IF

!** Algae
    WRITE (WIN,*) '    algae'
    DO J=1,NAAL
      JA = ALAC(J)
      IF (AG(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal growth rate [AG=',        AG(JA),'] < 0.0 for algal group ',    JA
      ELSE IF (AG(JA) < 0.1) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal growth rate [AG=',        AG(JA),'] < 0.1 for algal group ',    JA
      ELSE IF (AG(JA) > 5.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal growth rate [AG=',        AG(JA),'] > 5.0 for algal group ',    JA
      END IF
      IF (AM(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal mortality rate [AM=',     AM(JA),'] < 0.0 for algal group ',    JA
      ELSE IF (AM(JA) < 0.001) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal mortality rate [AM=',     AM(JA),'] < 0.001 for algal group ',  JA
      ELSE IF (AM(JA) > 1.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal mortality rate [AM=',     AM(JA),'] > 1.0 for algal group ',    JA
      END IF
      IF (AE(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal excretion rate [AE=',     AE(JA),'] < 0.0 for algal group ',    JA
      ELSE IF (AE(JA) < 0.001) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal excretion rate [AE=',     AE(JA),'] < 0.001 for algal group ',  JA
      ELSE IF (AE(JA) > 1.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal excretion rate [AE=',     AE(JA),'] < 1.0 for algal group ',    JA
      END IF
      IF (AR(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal excretion rate [AE=',     AE(JA),'] < 0.0 for algal group ',    JA
      ELSE IF (AR(JA) < 0.001) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal excretion rate [AE=',     AE(JA),'] < 0.001 for algal group ',  JA
      ELSE IF (AR(JA) > 1.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal excretion rate [AE=',     AE(JA),'] > 1.0 for algal group ',    JA
      END IF
      IF (AS(JA) > 5.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal settling rate [AS=',      AS(JA),'] > 5.0 for algal group ',    JA
      ELSE IF (AS(JA) < 0.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal settling rate [AS=',      AS(JA),'] < 0.0 for algal group ',    JA
      END IF
      IF (ASAT(JA) < 20.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal light saturation [ASAT=', ASAT(JA),'] < 20.0 for algal group ', JA
      ELSE IF (ASAT(JA) > 300.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal light saturation [ASAT=', ASAT(JA),'] > 300.0 for algal group ',JA
      END IF

!**** Algal temperature rate constants

      IF (AT1(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal minimum temperature [AT1=',AT1(JA),'] <  0.0 for algal group ',JA
      END IF
      IF (AT2(JA) > 40.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal maximum temperature [AT2=',AT2(JA),'] > 40.0 for algal group ',JA
      ELSE IF (AT2(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal maximum temperature [AT2=',AT2(JA),'] <  0.0 for algal group ',JA
      END IF
      IF (AT3(JA) > 40.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal maximum temperature [AT3=',AT3(JA),'] > 40.0 for algal group ',JA
      ELSE IF (AT3(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal maximum temperature [AT3=',AT3(JA),'] <  0.0 for algal group ',JA
      END IF
      IF (AT4(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal minimum temperature [AT4=',AT4(JA),'] <  0.0 for algal group ',JA
      END IF
      IF (AK1(JA) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal minimum temperature rate multiplier [AK1=',AK1(JA),'] <= 0.0 for algal group ',JA
      END IF
      IF (AK2(JA) > 1.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal maximum temperature rate multiplier [AK2=',AK2(JA),'] > 1.0 for algal group ',JA
      ELSE IF (AK2(JA) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal maximum temperature rate multiplier [AK2=',AK2(JA),'] <= 0.0 for algal group ',JA
      END IF
      IF (AK3(JA) > 1.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal maximum temperature rate multiplier [AK3=',AK3(JA),'] > 1.0 for algal group ',JA
      ELSE IF (AK3(JA) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal maximum temperature rate multiplier [AK3=',AK3(JA),'] <= 0.0 for algal group ',JA
      END IF
      IF (AK4(JA) > 1) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal minimum temperature rate multiplier [AK4=',AK4(JA),'] > 1.0 for algal group ',JA
      ELSE IF (AK4(JA) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal minimum temperature rate multiplier [AK4=',AK4(JA),'] <= 0.0 for algal group ',JA
      END IF
      IF (AT1(JA) >= AT2(JA)) THEN
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Algal temperature rate multiplier [AT1=',AT1(JA),'] > [AT2=',AT2(JA),'] for algal group ',JA
      END IF
      IF (AT2(JA) >= AT3(JA)) THEN
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Algal temperature rate multiplier [AT2=',AT2(JA),'] > [AT3=',AT3(JA),'] for algal group ',JA
      END IF
      IF (AT3(JA) >= AT4(JA)) THEN
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Algal temperature rate multiplier [AT3=',AT3(JA),'] > [AT4=',AT4(JA),'] for algal group ',JA
      END IF

!**** Algal stoichiometry

      IF (ALGP(JA) > 0.1) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Phosphorus to algal stoichiometry [ALGP=',ALGP(JA),'] > 0.1 for algal group ',   JA
      ELSE IF (ALGP(JA) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Phosphorus to algal stoichiometry [ALGP=',ALGP(JA),'] <= 0.0 for algal group ',   JA
      END IF
      IF (ALGN(JA) > 0.3) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Nitrogen to algal stoichiometry [ALGN=',  ALGN(JA),'] > 0.3 for algal group ',   JA
      ELSE IF (ALGN(JA) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Nitrogen to algal stoichiometry [ALGN=',  ALGN(JA),'] <= 0.0 for algal group ',   JA
      END IF
      IF (ALGC(JA) > 0.6) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Carbon to algal stoichiometry [ALGC=',    ALGC(JA),'] > 0.6 for algal group ',   JA
      ELSE IF (ALGC(JA) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Carbon to algal stoichiometry [ALGC=',    ALGC(JA),'] <= 0.0 for algal group ',   JA
      END IF
      IF (ALGSI(JA) > 0.5) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Silica to algal stoichiometry [ALGSI=',   ALGSI(JA),'] > 0.5 for algal group ',  JA
      ELSE IF (ALGSI(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Silica to algal stoichiometry [ALGSI=',   ALGSI(JA),'] < 0.0 for algal group ',  JA
      END IF
      IF (ACHLA(JA) > 0.8) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal to chlorophyll a ratio mg/ug [ACHLA=',    ACHLA(JA),'] > 0.8 for algal group ',JA
      ELSE IF (ACHLA(JA) < 0.01) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal to chlorophyll a ratio mg/ug [ACHLA=',    ACHLA(JA),'] < 0.01 for algal group ', JA
      ELSE IF (ACHLA(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal to chlorophyll a ratio mg/ug [ACHLA=',    ACHLA(JA),'] < 0.0 for algal group ',  JA
      END IF
      IF (APOM(JA) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Algal mortality to POM fraction [APOM=',  APOM(JA),'] < 0.0 for algal group ',   JA
      ELSE IF (APOM(JA) < 0.5) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Algal mortality to POM fraction [APOM='  ,APOM(JA),'] < 0.5 for algal group ',   JA
      END IF
      IF (ANEQN(JA) /= 1  .AND. ANEQN(JA) /= 2) THEN
        CALL ERRORS
        WRITE (ERR,'(A,I2,A,I2)') 'Algal N preference factor equation [ANEQN=',aneqn(ja),'] must be 1 or 2 for algal group:',   JA
      ENDIF
       IF(ANEQN(JA) == 2)THEN
       IF (ANPR(JA) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,'(A,f10.3,A,I2)') 'Algal N preference half saturation constant [ANPR=',anpr(ja),'] must be > 0.0 for algal group:',   JA
       ENDIF
       ENDIF

    END DO

        WRITE (WIN,*) '    algae toxins'
        IF(ALGAE_TOXINS)THEN
            INQUIRE(FILE='w2_Algae_Toxin.csv',EXIST=ALGAE_TOXIN_FILE)    
            IF(ALGAE_TOXIN_FILE)THEN
                 UNIT=2450
                 OPEN(2450,FILE='w2_Algae_Toxin.csv',STATUS='OLD')
                  AID='AlgToxin'
                  READ(2450,*,ERR=400)  ! SKIP HEADER
                  READ(2450,*,ERR=400)ATOX,ATOX_DEBUG     ! '(A2)'
                  IF(ATOX /= 'ON')THEN
                      CALL WARNINGS
                      WRITE(WRN,'(a,a2)')'Algae Toxins are active and w2_Algae_Toxins.csv exists, but global ON/OFF switch is not ON. ATOX=',atox
                  ENDIF
                  
                  
                  IF(ATOX == 'ON')THEN
                  ALLOCATE(CTP(NUMATOXINS,NAL),CTB(NUMATOXINS,NAL))
                  READ(2450,*,ERR=400)
                  READ(2450,*,ERR=400)(CTP(1,JA),JA=1,NAL)
                  READ(2450,*,ERR=400)(CTB(1,JA),JA=1,NAL)              
                  
                  READ(2450,*,ERR=400) CTREL(1)      !(CTL(1,JA),JA=1,NAL)
                  READ(2450,*,ERR=400) CTD(1)     !(CTDE(1,JA),JA=1,NAL)
                  
                  IF(CAC(NATS)=='      ON')THEN
                      DO JA=1,NAL
                          IF(CTP(1,JA) > 1.0 .OR. CTP(1,JA) < 0.0)THEN
                              CALL ERRORS
                              WRITE(ERR,'(A,F8.2,A,I2)')'Algae Toxins: In file w2_Algae_Toxin.csv for MC, CTP(1,JA)=',ctp(1,ja),' is > 1.0 or < 0.0 for algae group:',ja
                          ENDIF
                          IF(CTB(1,JA) >= 1.0 .OR. CTB(1,JA) < 0.0)THEN
                              CALL ERRORS
                              WRITE(ERR,'(A,F8.2,A,I2)')'Algae Toxins: In file w2_Algae_Toxin.csv for MC, CTB(1,JA)=',ctb(1,ja),' is >= 1.0 or < 0.0 for algae group:',ja
                          ENDIF
                      ENDDO
                  ENDIF

                  READ(2450,*,ERR=400)
                  READ(2450,*,ERR=400)(CTP(2,JA),JA=1,NAL)
                  READ(2450,*,ERR=400)(CTB(2,JA),JA=1,NAL)
                  READ(2450,*,ERR=400) CTREL(2)     !(CTL(2,JA),JA=1,NAL)
                  READ(2450,*,ERR=400) CTD(2)    !(CTDE(2,JA),JA=1,NAL)
                  
                  IF(CAC(NATS+1)=='      ON')THEN
                      DO JA=1,NAL
                          IF(CTP(2,JA) > 1.0 .OR. CTP(2,JA) < 0.0)THEN
                              CALL ERRORS
                              WRITE(ERR,'(A,F8.2,A,I2)')'Algae Toxins: In file w2_Algae_Toxin.csv for CYN, CTP(2,JA)=',ctp(2,ja),' is > 1.0 or < 0.0 for algae group:',ja
                          ENDIF
                          IF(CTB(2,JA) >= 1.0 .OR. CTB(2,JA) < 0.0)THEN
                              CALL ERRORS
                              WRITE(ERR,'(A,F8.2,A,I2)')'Algae Toxins: In file w2_Algae_Toxin.csv for CYN, CTB(2,JA)=',ctb(2,ja),' is >= 1.0 or < 0.0 for algae group:',ja
                          ENDIF
                      ENDDO
                  ENDIF
                  
                  READ(2450,*,ERR=400)
                  READ(2450,*,ERR=400)(CTP(3,JA),JA=1,NAL)
                  READ(2450,*,ERR=400)(CTB(3,JA),JA=1,NAL)
                  READ(2450,*,ERR=400) CTREL(3)     !(CTL(3,JA),JA=1,NAL)
                  READ(2450,*,ERR=400) CTD(3)    !(CTDE(3,JA),JA=1,NAL)
                  IF(CAC(NATS+2)=='      ON')THEN
                      DO JA=1,NAL
                          IF(CTP(3,JA) > 1.0 .OR. CTP(3,JA) < 0.0)THEN
                              CALL ERRORS
                              WRITE(ERR,'(A,F8.2,A,I2)')'Algae Toxins: In file w2_Algae_Toxin.csv for ATX, CTP(3,JA)=',ctp(3,ja),' is > 1.0 or < 0.0 for algae group:',ja
                          ENDIF
                          IF(CTB(3,JA) >= 1.0 .OR. CTB(3,JA) < 0.0)THEN
                              CALL ERRORS
                              WRITE(ERR,'(A,F8.2,A,I2)')'Algae Toxins: In file w2_Algae_Toxin.csv for ATX, CTB(3,JA)=',ctb(3,ja),' is >= 1.0 or < 0.0 for algae group:',ja
                          ENDIF
                      ENDDO
                  ENDIF
                  
                  READ(2450,*,ERR=400)
                  READ(2450,*,ERR=400)(CTP(4,JA),JA=1,NAL)
                  READ(2450,*,ERR=400)(CTB(4,JA),JA=1,NAL)
                  READ(2450,*,ERR=400) CTREL(4)     !(CTL(4,JA),JA=1,NAL)
                  READ(2450,*,ERR=400) CTD(4)    !(CTDE(4,JA),JA=1,NAL)
                  IF(CAC(NATS+3)=='      ON')THEN
                      DO JA=1,NAL
                          IF(CTP(4,JA) > 1.0 .OR. CTP(4,JA) < 0.0)THEN
                              CALL ERRORS
                              WRITE(ERR,'(A,F8.2,A,I2)')'Algae Toxins: In file w2_Algae_Toxin.csv for STX, CTP(4,JA)=',ctp(4,ja),' is > 1.0 or < 0.0 for algae group:',ja
                          ENDIF
                          IF(CTB(4,JA) >= 1.0 .OR. CTB(4,JA) < 0.0)THEN
                              CALL ERRORS
                              WRITE(ERR,'(A,F8.2,A,I2)')'Algae Toxins: In file w2_Algae_Toxin.csv for STX, CTB(4,JA)=',ctb(4,ja),' is >= 1.0 or < 0.0 for algae group:',ja
                          ENDIF
                      ENDDO
                  ENDIF       
                  CLOSE(2450)
            ELSE
                CALL ERRORS
                WRITE (ERR,'(A)') 'Algae Toxins are active but the input file w2_Algae_Toxin.csv does not exist.'
            ENDIF
            ENDIF
        ENDIF
        WRITE (WIN,*) '    algae vertical migration'
        IF(ALGAE_MIGRATION)THEN
            INQUIRE(FILE='w2_AlgaeMigration.csv',EXIST=ALGAE_MIGRATION_FILE)    
            IF(.NOT.ALGAE_MIGRATION_FILE)THEN
                 CALL ERRORS
                 WRITE (ERR,'(A)') 'AVERTM is ON for an algae group, but the file, w2_AlgaeMigration.csv, does not exist.'
            ELSE
                AID='w2_AlgaeMigration.csv';UNIT=2450
                      OPEN(2450,FILE='w2_AlgaeMigration.csv',STATUS='OLD')
      READ(2450,*,ERR=400)  ! SKIP HEADER
      READ(2450,*,ERR=400)MIGRATION     ! '(A2)'
      IF(MIGRATION /= 'ON' .OR. MIGRATION /='OF')THEN
                CALL ERRORS
                WRITE (ERR,'(A,A2)') 'w2_AlgaeMigration.csv: MIGRATION must be ON or OF[F]. [MIGRATION]=',MIGRATION    
      ENDIF
      IF(MIGRATION /= 'ON')GO TO 1001
      READ(2450,*,ERR=400)
      READ(2450,*,ERR=400)NMIG,ALGMIGRATION_DEBUG
      IF(ALGMIGRATION_DEBUG /= 1 .OR. ALGMIGRATION_DEBUG/=0)THEN
                CALL ERRORS
                WRITE (ERR,'(A,I3)') 'w2_AlgaeMigration.csv: ALGMIGRATION_DEBUG must be 1 or 1. [ALGMIGRATION_DEBUG]=',ALGMIGRATION_DEBUG    
      ENDIF
      
      READ(2450,*,ERR=400)
      ALLOCATE(MIGRATE_GROUP(NMIG),MIGRATE_MODEL(NMIG),AMP(NMIG),PHASE(NMIG),C_COEFF_EXT(NMIG),RAD(NMIG),MIND(NMIG),MAXD(NMIG),DENSI(NMIG),DENBI(NMIG),T_DEC(NMIG),TS_DEC(NMIG),C_DENINC(NMIG),C_DENDEC(NMIG),&
          DEPTH_LIM_ONOFF(NMIG),DEPTH_LIM(NMIG),LOSS_FRAC(NMIG),I_C(NMIG),C_DENINC_1(NMIG),C_DENINC_2(NMIG),C_DENDEC_1(NMIG),C_DENDEC_2(NMIG),DENP_MINS(NMIG),DENP_MINB(NMIG),DEN_COR(MIGI),&
          NMINT(NMIG),MIGON(NMIG,12),MIGOFF(NMIG,12),DEPTH_CALC_ONOFF(NMIG),EXP_DEPTH(NMIG))
      DO I=1,NMIG
        READ(2450,*,ERR=400)MIGRATE_GROUP(I)
        IF(MIGRATE_GROUP(I)==0 .OR. MIGRATE_GROUP(I)>NAL)THEN
                CALL ERRORS
                WRITE (ERR,'(A,I2,a,i2)') 'w2_AlgaeMigration.csv: MIGRATE_GROUP(I) must be >0 AND <# of algal groups. [MIGRATE_GROUP(I)]=',MIGRATE_GROUP(I),' for group I=',I    
        ENDIF

        READ(2450,*,ERR=400)
        READ(2450,*,ERR=400)NMINT(I)
        READ(2450,*,ERR=400)
        READ(2450,*,ERR=400)MIGON(I,1:NMINT(I))
        READ(2450,*,ERR=400)
        READ(2450,*,ERR=400)MIGOFF(I,1:NMINT(I))
        READ(2450,*,ERR=400)
        READ(2450,*,ERR=400)MIGRATE_MODEL(I)
        IF(MIGRATE_MODEL(I) <=0 .OR. MIGRATE_MODEL(I)>4)THEN
                CALL ERRORS
                WRITE (ERR,'(A,I2,a,i2)') 'w2_AlgaeMigration.csv: MIGRATE_MODEL(I) must be between 1 AND 4. [MIGRATE_MODEL(I)]=',MIGRATE_MODEL(I),' for group I=',I    
        ENDIF
        READ(2450,*,ERR=400)
        IF(MIGRATE_MODEL(I)==1 .OR. MIGRATE_MODEL(I)==2)THEN
            READ(2450,*,ERR=400)AMP(I),PHASE(I),C_COEFF_EXT(I),DEPTH_CALC_ONOFF(I),EXP_DEPTH(I),DEPTH_LIM_ONOFF(I),DEPTH_LIM(I),LOSS_FRAC(I)
            READ(2450,*,ERR=400)
            READ(2450,*,ERR=400)
            READ(2450,*,ERR=400)
            READ(2450,*,ERR=400)
            READ(2450,*,ERR=400)
        ELSEIF(MIGRATE_MODEL(I)==3)THEN
            READ(2450,*,ERR=400)
            READ(2450,*,ERR=400)
            READ(2450,*,ERR=400)RAD(I),MIND(I),MAXD(I),DENSI(I),DENBI(I),T_DEC(I),TS_DEC(I),C_DENINC(I),C_DENDEC(I),DEPTH_LIM_ONOFF(I),DEPTH_LIM(I),LOSS_FRAC(I)
            READ(2450,*,ERR=400)
            READ(2450,*,ERR=400)
            READ(2450,*,ERR=400)
        ELSE
            READ(2450,*,ERR=400)
            READ(2450,*,ERR=400)            
            READ(2450,*,ERR=400)
            READ(2450,*,ERR=400)            
            READ(2450,*,ERR=400)RAD(I),MIND(I),MAXD(I),DENSI(I),DENBI(I),I_C(I),C_DENINC_1(I),C_DENINC_2(I),C_DENDEC_1(I),C_DENDEC_2(I),DENP_MINS(I),DENP_MINB(I),DEN_COR(I),DEPTH_LIM_ONOFF(I),&
                DEPTH_LIM(I),LOSS_FRAC(I)
            READ(2450,*,ERR=400)            
        ENDIF
      ENDDO

1001   CLOSE(2450)

            ENDIF
                    
        ELSE
             INQUIRE(FILE='w2_AlgaeMigration.csv',EXIST=ALGAE_MIGRATION_FILE) 
             IF(ALGAE_MIGRATION)THEN
                 CALL WARNINGS
                 WRITE (WRN,'(A)') 'AVERTM is OFF for all algae groups, but the file, w2_AlgaeMigration.csv, exists but is not implemented in the model.'
             ENDIF
        ENDIF
        


!** Epiphyton

    WRITE (WIN,*) '    epiphyton'
    DO JE=1,NEP
      IF (EG(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton growth rate [EG=',        EG(JE),'] < 0.0 for epiphyton group ',    JE
      ELSE IF (EG(JE) < 0.1) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Epiphyton growth rate [EG=',        EG(JE),'] < 0.1 for epiphyton group ',    JE
      ELSE IF (EG(JE) > 5.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Epiphyton growth rate [EG=',        EG(JE),'] > 5.0 for epiphyton group ',    JE
      END IF
      IF (EM(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton mortality rate [EM=',     EM(JE),'] < 0.0 for epiphyton group ',    JE
      ELSE IF (EM(JE) < 0.001) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Epiphyton mortality rate [EM=',     EM(JE),'] < 0.001 for epiphyton group ',  JE
      ELSE IF (EM(JE) > 1.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Epiphyton mortality rate [EM=',     EM(JE),'] > 1.0 for epiphyton group ',    JE
      END IF
      IF (EE(JE) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton excretion rate [EE=',     EE(JE),'] < 0.0 for epiphyton group ',    JE
      ELSE IF (EE(JE) < 0.001) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Epiphyton excretion rate [EE=',     EE(JE),'] < 0.001 for epiphyton group ',  JE
      ELSE IF (EE(JE) > 1.0) THEN
         CALL WARNINGS
         WRITE (WRN,FMTFI) 'Epiphyton excretion rate [EE=',    EE(JE),'] < 1.0 for epiphyton group ',    JE
      END IF
      IF (ER(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton respiration rate [ER=',   ER(JE),'] < 0.0 for epiphyton group ',    JE
      ELSE IF (ER(JE) < 0.001) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Epiphyton respiration rate [ER=',   ER(JE),'] < 0.001 for epiphyton group ',  JE
      ELSE IF (ER(JE) > 1.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Epiphyton respiration rate [ER=',   ER(JE),'] > 1.0 for epiphyton group ',    JE
      END IF
      IF (EB(JE) > 0.01) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Epiphyton burial rate [ES=',        EB(JE),'] > 0.01 for epiphyton group ',   JE
      ELSE IF (EB(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton burial rate [EB=',        EB(JE),'] < 0.0 for epiphyton group ',    JE
      ELSE IF (EB(JE) == 0.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Epiphyton burial rate [EB=',        EB(JE),'] = 0.0 for epiphyton group ',    JE

      END IF
      IF (ESAT(JE) < 20.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton light saturation [ESAT=', ESAT(JE),'] < 20.0 for epiphyton group ', JE
      ELSE IF (ESAT(JE) > 300.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton light saturation [ESAT=', ESAT(JE),'] > 300.0 for epiphyton group ',JE
      END IF

!**** Epiphyton temperature rate constants

      IF (ET1(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton minimum temperature [ET1=',ET1(JE),'] <  0.0 for epiphyton group ',JE
      END IF
      IF (ET2(JE) > 40.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Epiphyton maximum temperature [ET2=',ET2(JE),'] > 40.0 for epiphyton group ',JE
      ELSE IF (ET2(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton maximum temperature [ET2=',ET2(JE),'] <  0.0 for epiphyton group ',JE
      END IF
      IF (ET3(JE) > 40.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Epiphyton maximum temperature [ET3=',ET3(JE),'] > 40.0 for epiphyton group ',JE
      ELSE IF (ET3(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton maximum temperature [ET3=',ET3(JE),'] <  0.0 for epiphyton group ',JE
      END IF
      IF (ET4(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton minimum temperature [ET4=',ET4(JE),'] <  0.0 for epiphyton group ',JE
      END IF
      IF (EK1(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton minimum temperature rate multiplier [EK1=',EK1(JE),'] < 0.0 for epiphyton group ',JE
      END IF
      IF (EK2(JE) > 1.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton maximum temperature rate multiplier [EK2=',EK2(JE),'] > 1.0 for epiphyton group ',JE
      ELSE IF (EK2(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton maximum temperature rate multiplier [EK2=',EK2(JE),'] < 0.0 for epiphyton group ',JE
      END IF
      IF (EK3(JE) > 1.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton maximum temperature rate multiplier [EK3=',EK3(JE),'] > 1.0 for epiphyton group ',JE
      ELSE IF (EK3(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton maximum temperature rate multiplier [EK3=',EK3(JE),'] < 0.0 for epiphyton group ',JE
      END IF
      IF (EK4(JE) > 1) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton minimum temperature rate multiplier [EK4=',EK4(JE),'] > 1.0 for epiphyton group ',JE
      ELSE IF (EK4(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton minimum temperature rate multiplier [EK4=',EK4(JE),'] < 0.0 for epiphyton group ',JE
      END IF
      IF (ET1(JE) >= ET2(JE)) THEN
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Epiphyton temperature rate multiplier [ET1=',ET1(JE),'] > [ET2=',ET2(JE),'] for epiphyton group ',JE
      END IF
      IF (ET2(JE) >= ET3(JE)) THEN
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Epiphyton temperature rate multiplier [ET2=',ET2(JE),'] > [ET3=',ET3(JE),'] for epiphyton group ',JE
      END IF
      IF (ET3(JE) >= ET4(JE)) THEN
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Epiphyton temperature rate multiplier [ET3=',ET3(JE),'] > [ET4=',ET4(JE),'] for epiphyton group ',JE
      END IF

!**** Epiphyton stoichiometry

      IF (EP(JE) > 0.1) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Phosphorus to epiphyton stoichiometry [EP=', EP(JE),'] > 0.1 for epiphyton group ',     JE
      ELSE IF (EP(JE) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Phosphorus to epiphyton stoichiometry [EP=', EP(JE),'] <= 0.0 for epiphyton group ',     JE
      END IF
      IF (EN(JE) > 0.3) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Nitrogen to epiphyton stoichiometry [EN=',   EN(JE),'] > 0.3 for epiphyton group ',     JE
      ELSE IF (EN(JE) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Nitrogen to epiphyton stoichiometry [EN=',   EN(JE),'] <= 0.0 for epiphyton group ',     JE
      END IF
      IF (EC(JE) > 0.6) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Carbon to epiphyton stoichiometry [EC=',     EC(JE),'] > 0.6 for epiphyton group ',     JE
      ELSE IF (EC(JE) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Carbon to epiphyton stoichiometry [EC=',     EC(JE),'] <= 0.0 for epiphyton group ',     JE
      END IF
      IF (ESI(JE) > 0.5) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Silica to epiphyton stoichiometry [ESI=',    ESI(JE),'] > 0.5 for epiphyton group ',    JE
      ELSE IF (ESI(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Silica to epiphyton stoichiometry [ESI=',    ESI(JE),'] < 0.0 for epiphyton group ',    JE
      END IF
      IF (ECHLA(JE) > 0.8) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton to chlorophyll a ratio mg/ug [ECHLA=',   ECHLA(JE),'] > 0.8 for epiphyton group ',JE
      ELSE IF (ECHLA(JE) < 0.01) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Epiphyton to chlorophyll a ratio mg/ug [ECHLA=',   ECHLA(JE),'] < 0.01 for epiphyton group ', JE
      ELSE IF (ECHLA(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton to chlorophyll a ratio mg/ug [ECHLA=',   ECHLA(JE),'] < 0.0 for epiphyton group ',  JE
      END IF
      IF (EPOM(JE) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Epiphyton mortality to POM fraction [EPOM=', EPOM(JE),'] < 0.0 for epiphyton group ',   JE
      ELSE IF (EPOM(JE) < 0.5) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Epiphyton mortality to POM fraction [EPOM=' ,EPOM(JE),'] < 0.5 for epiphyton group ',   JE
      END IF

      do jw=1,nwb
        IF (EPIC(JW,JE) /= '     OFF' .and. epic(JW,JE) /= '      ON') THEN
        CALL ERRORS
        WRITE (ERR,'(3(A,A,A,I0,a,i0))') 'Epiphyton ON/OFF control [EPIC=', EPIC(JW,JE),'] must be ON or OFF for epiphyton group ',   JE,' and JW=',JW
        endif
        IF (EPIPRC(JW,JE) /= '     OFF' .and. epiPRc(JW,JE) /= '      ON') THEN
        CALL ERRORS
        WRITE (ERR,'(3(A,A,A,I0,a,i0))') 'Epiphyton ON/OFF print control [EPIPRC=', EPIPRC(JW,JE),'] must be ON or OFF for epiphyton group ',   JE,' and JW=',JW
        endif
      enddo

    END DO

!** Macrophytes ******************

    WRITE (WIN,*) '    macropytes'

    IF(NMC == 0)THEN
    DO jw=1,nwb
        IF (macwbc(JW,1) == '      ON') THEN
        CALL ERRORS
        WRITE (ERR,'(3(A,A,A,I0,a,i0))') 'Macrophyte ON/OFF control [MACWBC=', MACWBC(JW,1),'] is ON even though NMC[# of macrophyte groups] is set to zero for waterbody=',JW
        endif
    ENDDO



    ENDIF


    DO JM=1,NMC

      do jw=1,nwb
        IF (MACWBC(JW,JM) /= '     OFF' .and. macwbc(JW,Jm) /= '      ON') THEN
        CALL ERRORS
        WRITE (ERR,'(3(A,A,A,I0,a,i0))') 'Macrophyte ON/OFF control [MACWBC=', MACWBC(JW,JM),'] must be ON or OFF for macrophyte group ',   JM,' and JW=',JW
        endif
        IF (MPRWBC(JW,JM) /= '     OFF' .and. mprwbc(JW,Jm) /= '      ON') THEN
        CALL ERRORS
        WRITE (ERR,'(3(A,A,A,I0,a,i0))') 'Macrophyte ON/OFF print control [MPRWBC=', MPRWBC(JW,JM),'] must be ON or OFF for macrophyte group ',   JM,' and JW=',JW
        endif
      enddo

      IF (MG(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte growth rate [MG=',        MG(JM),'] < 0.0 for macrophyte group ',    JM
      ELSE IF (MG(JM) < 0.01) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte growth rate [MG=',        MG(JM),'] < 0.01 for macrophyte group ',    JM
      ELSE IF (MG(JM) > 2.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte growth rate [MG=',        MG(JM),'] > 2.0 for macrophyte group ',    JM
      END IF
      IF (MM(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte mortality rate [MM=',     MM(JM),'] < 0.0 for macrophyte group ',    JM
      ELSE IF (MM(JM) < 0.001) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte mortality rate [MM=',     MM(JM),'] < 0.001 for macrophyte group ',  JM
      ELSE IF (MM(JM) > 1.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte mortality rate [MM=',     MM(JM),'] > 1.0 for macrophyte group ',    JM
      END IF
      IF (MR(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte respiration rate [MR=',   MR(JM),'] < 0.0 for macrophyte group ',    JM
      ELSE IF (MR(JM) < 0.001) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte respiration rate [MR=',   MR(JM),'] < 0.001 for macrophyte group ',  JM
      ELSE IF (MR(JM) > 1.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte respiration rate [MR=',   MR(JM),'] > 1.0 for macrophyte group ',    JM
      END IF

      IF (MSAT(JM) < 15.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte light saturation [MSAT=', MSAT(JM),'] < 15.0 for macrophyte group ', JM
      ELSE IF (MSAT(JM) > 300.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte light saturation [MSAT=', MSAT(JM),'] > 300.0 for macrophyte group ',JM
      END IF


      IF (PSED(JM) < 0.0 .or. PSED(JM) > 1.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte sediment P uptake fraction [PSED=', PSED(JM),'] < 0 or > 1 for macrophyte group ', JM
      END IF
      IF (NSED(JM) < 0.0 .or. NSED(JM) > 1.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte sediment N uptake fraction [NSED=', NSED(JM),'] < 0 or > 1 for macrophyte group ', JM
      END IF
      IF (MBMP(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte threshold conc to move upward [MBMP=', MBMP(JM),'] < 0 for macrophyte group ', JM
       ELSE IF (MBMP(JM) > 150.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte threshold conc to move upward [MBMP=', MBMP(JM),'] > 150 for macrophyte group ', JM
      END IF
      IF (MMAX(JM) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte maximum concentration [MMAX=', MMAX(JM),'] < or = 0 for macrophyte group ', JM
       ELSE IF (MMAX(JM) > 1500.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte maximum concentration [MMAX=', MMAX(JM),'] > 1500 for macrophyte group ', JM
      END IF

      IF (dmv(JM) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte dry weight to wet volume ratio [DMV=', DMV(JM),'] < or = 0 for macrophyte group ', JM
       ELSE IF (dmv(JM) > 1.0e7) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte dry weight to wet volume ratio [DMV=', DMV(JM),'] > 1E7 for macrophyte group ', JM
      END IF
      IF (cdstem(JM) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte CD coefficient [CDDRAG=', cdstem(JM),'] < or = 0 for macrophyte group ', JM
       ELSE IF (cdstem(JM) < 1.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte CD coefficient [CDDRAG=', cdstem(JM),'] < 1 for macrophyte group ', JM
       ELSE IF (cdstem(JM) > 3.5) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte CD coefficient [CDDRAG=', cdstem(JM),'] > 3.5 for macrophyte group ', JM
      END IF

        IF (dwsa(JM) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte dry weight to surface area [DWSA=', DWSA(JM),'] < or = 0 for macrophyte group ', JM
       ELSE IF (dwsa(JM) > 100.) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte dry weight to surface area [DWSA=', DWSA(JM),'] > 100. for macrophyte group ', JM
      END IF

      IF (anorm(JM) <= 0.0 .or. anorm(jm) > 1.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte area fraction normal to flow [ANORM=', ANORM(JM),'] < or = 0 or > 1.0 for macrophyte group ', JM
      END IF

!**** Macrophyte temperature rate constants

      IF (MT1(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte minimum temperature [MT1=',MT1(JM),'] <  0.0 for macrophyte group ',JM
      END IF
      IF (MT2(JM) > 40.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte maximum temperature [MT2=',MT2(JM),'] > 40.0 for macrophyte group ',JM
      ELSE IF (MT2(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte maximum temperature [MT2=',MT2(JM),'] <  0.0 for macrophyte group ',JM
      END IF
      IF (MT3(JM) > 40.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte maximum temperature [MT3=',MT3(JM),'] > 40.0 for macrophyte group ',JM
      ELSE IF (MT3(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte maximum temperature [MT3=',MT3(JM),'] <  0.0 for macrophyte group ',JM
      END IF
      IF (MT4(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte minimum temperature [MT4=',MT4(JM),'] <  0.0 for macrophyte group ',JM
      END IF
      IF (MK1(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte minimum temperature rate multiplier [MK1=',MK1(JM),'] < 0.0 for macrophyte group ',JM
      END IF
      IF (MK2(JM) > 1.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte maximum temperature rate multiplier [MK2=',MK2(JM),'] > 1.0 for macrophyte group ',JM
      ELSE IF (MK2(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte maximum temperature rate multiplier [MK2=',MK2(JM),'] < 0.0 for macrophyte group ',JM
      END IF
      IF (MK3(JM) > 1.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte maximum temperature rate multiplier [MK3=',MK3(JM),'] > 1.0 for macrophyte group ',JM
      ELSE IF (MK3(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte maximum temperature rate multiplier [MK3=',MK3(JM),'] < 0.0 for macrophyte group ',JM
      END IF
      IF (MK4(JM) > 1) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte minimum temperature rate multiplier [MK4=',MK4(JM),'] > 1.0 for macrophyte group ',JM
      ELSE IF (MK4(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte minimum temperature rate multiplier [MK4=',MK4(JM),'] < 0.0 for macrophyte group ',JM
      END IF
      IF (MT1(JM) >= MT2(JM)) THEN
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Macrophyte temperature rate multiplier [MT1=',MT1(JM),'] > [MT2=',MT2(JM),'] for macrophyte group ',JM
      END IF
      IF (MT2(JM) >= MT3(JM)) THEN
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Macrophyte temperature rate multiplier [MT2=',MT2(JM),'] > [MT3=',MT3(JM),'] for macrophyte group ',JM
      END IF
      IF (MT3(JM) >= MT4(JM)) THEN
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Macrophyte temperature rate multiplier [MT3=',MT3(JM),'] > [MT4=',MT4(JM),'] for macrophyte group ',JM
      END IF

!**** Macrophyte stoichiometry

      IF (MP(JM) > 0.1) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Phosphorus to macrophyte stoichiometry [MP=', MP(JM),'] > 0.1 for macrophyte group ',     JM
      ELSE IF (MP(JM) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Phosphorus to macrophyte stoichiometry [MP=', MP(JM),'] <= 0.0 for macrophyte group ',     JM
      END IF
      IF (MN(JM) > 0.3) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Nitrogen to macrophyte stoichiometry [MN=',   MN(JM),'] > 0.3 for macrophyte group ',     JM
      ELSE IF (MN(JM) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Nitrogen to macrophyte stoichiometry [MN=',   MN(JM),'] <= 0.0 for macrophyte group ',     JM
      END IF
      IF (MC(JM) > 0.6) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Carbon to macrophyte stoichiometry [MC=',     MC(JM),'] > 0.6 for macrophyte group ',     JM
      ELSE IF (MC(JM) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Carbon to macrophyte stoichiometry [MC=',     MC(JM),'] <= 0.0 for macrophyte group ',     JM
      END IF

      IF (MPOM(JM) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Macrophyte mortality to POM fraction [MPOM=', MPOM(JM),'] < 0.0 for macrophyte group ',   JM
      ELSE IF (MPOM(JM) < 0.5) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Macrophyte mortality to POM fraction [MPOM=' ,MPOM(JM),'] < 0.5 for macrophyte group ',   JM
      END IF
    END DO

!**Zooplankton **************
   IF(NZP>0)THEN 
    do jz=1,nzp
    IF(ZG(JZ) <= 0.0)THEN
      CALL ERRORS
        WRITE (ERR,FMTFI) 'Zooplankton growth rate [ZG=', ZG(JZ),'] <= 0.0 for zooplankton group ',   JZ  
    ENDIF
    IF(ZG(JZ) > 3.0)THEN
      CALL WARNINGS
        WRITE (WRN,FMTFI) 'Zooplankton growth rate [ZG=', ZG(JZ),'] > 3.0 for zooplankton group ',   JZ  
    ENDIF
    IF(ZR(JZ) <= 0.0)THEN
      CALL ERRORS
        WRITE (ERR,FMTFI) 'Zooplankton respiration rate [ZR=', ZR(JZ),'] <= 0.0 for zooplankton group ',   JZ  
    ENDIF
    IF(ZR(JZ) > 1.0)THEN
      CALL WARNINGS
        WRITE (WRN,FMTFI) 'Zooplankton respiration rate [ZR=', ZR(JZ),'] > 1.0 for zooplankton group ',   JZ  
    ENDIF
    IF(ZM(JZ) <= 0.0)THEN
      CALL ERRORS
        WRITE (ERR,FMTFI) 'Zooplankton mortality rate [ZM=', ZM(JZ),'] <= 0.0 for zooplankton group ',   JZ  
    ENDIF
    IF(ZM(JZ) > 0.5)THEN
      CALL WARNINGS
        WRITE (WRN,FMTFI) 'Zooplankton mortality rate [ZM=', ZM(JZ),'] > 0.5 for zooplankton group ',   JZ  
    ENDIF
    
    IF (ZP(JZ) > 0.1) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Phosphorus to zooplankton stoichiometry [ZP=', ZP(JZ),'] > 0.1 for zooplankton group ',     JZ
      ELSE IF (ZP(JZ) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Phosphorus to zooplankton stoichiometry [ZP=', ZP(JZ),'] <= 0.0 for zooplankton group ',     JZ
      END IF
      IF (ZN(JZ) > 0.3) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Nitrogen to zooplankton stoichiometry [ZN=',   ZN(JZ),'] > 0.3 for zooplankton group ',     JZ
      ELSE IF (ZN(JZ) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Nitrogen to zooplankton stoichiometry [ZN=',   ZN(JZ),'] <= 0.0 for zooplankton group ',     JZ
      END IF
      IF (ZC(JZ) > 0.6) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Carbon to zooplankton stoichiometry [ZC=',     ZC(JZ),'] > 0.6 for zooplankton group ',     JZ
      ELSE IF (ZC(JZ) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Carbon to zooplankton stoichiometry [ZC=',     ZC(JZ),'] <= 0.0 for zooplankton group ',     JZ
      END IF    
      
      ! zooplankton temperature rates
        IF (ZT1(JZ) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'zooplankton minimum temperature [ZT1=',ZT1(JZ),'] <  0.0 for zooplankton group ',JZ
      END IF
      IF (ZT2(JZ) > 40.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'zooplankton maximum temperature [ZT2=',ZT2(JZ),'] > 40.0 for zooplankton group ',JZ
      ELSE IF (ZT2(JZ) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'zooplankton maximum temperature [ZT2=',ZT2(JZ),'] <  0.0 for zooplankton group ',JZ
      END IF
      IF (ZT3(JZ) > 40.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'zooplankton maximum temperature [ZT3=',ZT3(JZ),'] > 40.0 for zooplankton group ',JZ
      ELSE IF (ZT3(JZ) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'zooplankton maximum temperature [ZT3=',ZT3(JZ),'] <  0.0 for zooplankton group ',JZ
      END IF
      IF (ZT4(JZ) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'zooplankton minimum temperature [ZT4=',ZT4(JZ),'] <  0.0 for zooplankton group ',JZ
      END IF
      IF (ZK1(JZ) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'zooplankton minimum temperature rate multiplier [ZK1=',ZK1(JZ),'] <= 0.0 for zooplankton group ',JZ
      END IF
      IF (ZK2(JZ) > 1.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'zooplankton maximum temperature rate multiplier [ZK2=',ZK2(JZ),'] > 1.0 for zooplankton group ',JZ
      ELSE IF (ZK2(JZ) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Zooplankton maximum temperature rate multiplier [ZK2=',ZK2(JZ),'] <= 0.0 for zooplankton group ',JZ
      END IF
      IF (ZK3(JZ) > 1.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Zooplankton maximum temperature rate multiplier [ZK3=',ZK3(JZ),'] > 1.0 for zooplankton group ',JZ
      ELSE IF (ZK3(JZ) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Zooplankton maximum temperature rate multiplier [ZK3=',ZK3(JZ),'] <= 0.0 for zooplankton group ',JZ
      END IF
      IF (ZK4(JZ) > 1) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Zooplankton minimum temperature rate multiplier [ZK4=',ZK4(JZ),'] > 1.0 for zooplankton group ',JZ
      ELSE IF (ZK4(JZ) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Zooplankton minimum temperature rate multiplier [ZK4=',ZK4(JZ),'] <= 0.0 for zooplankton group ',JZ
      END IF
      IF (ZT1(JZ) >= ZT2(JZ)) THEN
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Zooplankton temperature rate multiplier [ZT1=',ZT1(JZ),'] > [ZT2=',ZT2(JZ),'] for zooplankton group ',JZ
      END IF
      IF (ZT2(JZ) >= ZT3(JZ)) THEN
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Zooplankton temperature rate multiplier [ZT2=',ZT2(JZ),'] > [ZT3=',ZT3(JZ),'] for zooplankton group ',JZ
      END IF
      IF (ZT3(JZ) >= ZT4(JZ)) THEN
        CALL ERRORS
        WRITE (ERR,FMT2FI) 'Zooplankton temperature rate multiplier [ZT3=',ZT3(JZ),'] > [ZT4=',ZT4(JZ),'] for zooplankton group ',JZ
      END IF     
    enddo
   ENDIF
   

!** Dissolved organic matter

    WRITE (WIN,*) '    dissolved organic matter'
    IF (LABILE_DOM) THEN
      DO JW=1,NWB
        IF (LDOMDK(JW) > 0.3) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Labile dissolved organic matter decay rate [LDOMDK=',LDOMDK(JW),'] > 0.3 for waterbody ',JW
        ELSE IF (LDOMDK(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Labile dissolved organic matter decay rate [LDOMDK=',LDOMDK(JW),'] < 0.0 for waterbody ',JW
        END IF
        IF (LRDDK(JW) > 0.01) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Labile to refractory DOM decay rate [LRDDK=',LRDDK(JW),'] > 0.01 for waterbody ',JW
        ELSE IF (LRDDK(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Labile to refractory DOM decay rate [LRDDK=',LRDDK(JW),'] <  0.0 for waterbody ',JW
        END IF
      END DO
    END IF
    IF (REFRACTORY_DOM) THEN
      DO JW=1,NWB
        IF (RDOMDK(JW) > 0.01) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Refractory dissolved organic matter decay rate [RDOMDK=',RDOMDK(JW),'] > 0.01 for waterbody ',JW
        ELSE IF (RDOMDK(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Refractory dissolved organic matter decay rate [RDOMDK=',RDOMDK(JW),'] < 0.0 for waterbody ',JW
        END IF
      END DO
    END IF

!** Particulate organic matter

    WRITE (WIN,*) '    particulate organic matter'
    IF (LABILE_POM) THEN
      DO JW=1,NWB
        IF (LPOMDK(JW) > 0.5) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Labile particulate organic matter decay rate [LPOMDK=',LPOMDK(JW),'] > 0.5 for waterbody ',JW
        ELSE IF (LPOMDK(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Labile particulate organic matter decay rate [LPOMDK=',LPOMDK(JW),'] < 0.0 for waterbody ',JW
        END IF
        IF (POMS(JW) > 1.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Particulate organic matter settling rate [POMS=',POMS(JW),'] > 1.0 for waterbody ',JW
        ELSE IF (POMS(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Particulate organic matter settling rate [POMS=',POMS(JW),'] < 0.0 for waterbody ',JW
        END IF
      END DO
    END IF
    IF (REFRACTORY_POM) THEN
      DO JW=1,NWB
        IF (RPOMDK(JW) > 0.5) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Refractory particulate organic matter decay rate [RPOMDK=',RPOMDK(JW),'] > 0.5 for waterbody ',JW
        ELSE IF (RPOMDK(JW) < 0.0) THEN
          WRITE (WRN,FMTFI) 'Refractory particulate organic matter decay rate [RPOMDK=',RPOMDK(JW),'] < 0.0 for waterbody ',JW
          CALL ERRORS
        END IF
      END DO
    END IF

!** Organic matter temperature rate multipliers

    IF (LABILE_DOM .OR. REFRACTORY_DOM .OR. LABILE_POM .OR. REFRACTORY_POM) THEN
      DO JW=1,NWB
        IF (OMT1(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Organic matter decay minimum temperature [OMT1=',OMT1(JW),'] < 0.0 for waterbody ', JW
        ELSE IF (OMT1(JW) > 15) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Organic matter decay minimum temperature [OMT1=',OMT1(JW),'] > 15.0 for waterbody ',JW
        END IF
        IF (OMT2(JW) > 40.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Organic matter decay maximum temperature [OMT2=',OMT2(JW),'] > 40.0 for waterbody ',JW
        ELSE IF (OMT2(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Organic matter decay maximum temperature [OMT2=',OMT2(JW),'] < 0.0 for waterbody ', JW
        END IF
        IF (OMK1(JW) <= 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Organic matter decay minimum temperature rate multiplier [OMK1=',OMK1(JW),'] <= 0.0 for waterbody ',JW
        ELSE IF (OMK1(JW) > 0.5) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Organic matter decay minimum temperature rate multiplier [OMK1=',OMK1(JW),'] > 0.5 for waterbody ',JW
        END IF
        IF (OMK2(JW) > 1.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Organic matter decay maximum temperature rate multiplier [OMK2=',OMK2(JW),'] > 1.0 for waterbody ',JW
        ELSE IF (OMK2(JW) <= 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Organic matter decay maximum temperature rate multiplier [OMK2=',OMK2(JW),'] <= 0.0 for waterbody ',JW
        END IF
      END DO
    END IF

!** Sediments

    WRITE (WIN,*) '    sediments'
    DO JW=1,NWB
      IF (SEDIMENT_CALC(JW)) THEN
        IF (SEDDK(JW) > 0.3) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Sediment decay rate [SEDDK=',SEDDK(JW),'] > 0.3 for waterbody ',JW
        ELSE IF (SEDDK(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Sediment decay rate [SEDDK=',SEDDK(JW),'] < 0.0 for waterbody ',JW
        END IF
        IF (SEDs(JW) > 1.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Sediment focusing velocity or settling rate [SEDS=',SEDS(JW),'] > 1.0 m/d for waterbody ',JW
        ELSE IF (SEDs(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Sediment focusing velocity or settling rate [SEDS=',SEDS(JW),'] < 0.0 for waterbody ',JW
        END IF

        IF (SEDbr(JW) > 1.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Sediment burial rate [SEDBR=',SEDBR(JW),'] > 1.0 d^-1 for waterbody ',JW
        ELSE IF (SEDbr(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Sediment burial rate [SEDBR=',SEDBR(JW),'] < 0.0 for waterbody ',JW
        END IF

        IF (DYNSEDK(JW) /= '     OFF' .AND. DYNSEDK(JW) /= '      ON')THEN
          CALL ERRORS
          WRITE (ERR,"('Dynamic calculation of sediment decay control [DYNSEDK=',A8,'] must be ON or OFF for waterbody ',I4)")DYNSEDK(JW),JW
        END IF
      END IF

        IF (SODT1(JW) <= 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Sediment decay minimum temperature [SODT1=',SODT1(JW),'] <= 0 C for waterbody ', JW
        ELSE IF (SODT1(JW) > 10.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Sediment decay minimum temperature [SODT1=',SODT1(JW),'] > 10 C for waterbody ',JW
        END IF
        IF (SODT2(JW) > 40.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Sediment decay maximum temperature [SODT2=',SODT2(JW),'] > 40 C for waterbody ',JW
        ELSE IF (SODT2(JW) <= 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Sediment decay maximum temperature [SODT2=',SODT2(JW),'] <= 0 C for waterbody ', JW
        END IF
        IF (SODK1(JW) <= 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Sediment decay minimum temperature rate multiplier [SODK1=',SODK1(JW),'] <= 0.0 for waterbody ',JW
        ELSE IF (SODK1(JW) > 0.5) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Sediment decay minimum temperature rate multiplier [SODK1=',SODK1(JW),'] > 0.5 for waterbody ',JW
        END IF
        IF (SODK2(JW) > 1.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Sediment decay maximum temperature rate multiplier [SODK1=',SODK1(JW),'] > 1.0 for waterbody ',JW
        ELSE IF (SODK2(JW) <= 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Sediment decay maximum temperature rate multiplier [SODK1=',SODK1(JW),'] <= 0.0 for waterbody ',JW
        END IF
        !IF(FSOD(JW)>0.0 .AND. SEDIMENT_DIAGENESIS)THEN    ! SW 7/24/2017 CHECK ON SOD=0 IF SED DIAG IS ON
        !    DO JB=BS(JW),BE(JW)
        !        DO I=US(JB),DS(JB)
        !            IF(SOD(I)>0.0)THEN
        !                CALL ERRORS
        !                WRITE (ERR,FMTFI) 'Sediment diagenesis is ON. Zero order SOD must be set to zero (or just set FSOD=0). SOD=',SOD(I),' for segment ',I
        !            ENDIF
        !        ENDDO
        !    ENDDO
        !ENDIF
        
    END DO
    

!** Carbonaceous biochemical oxygen demand

    WRITE (WIN,*) '    CBOD'
    IF (CBO_DEMAND) THEN
      DO JB=1,NBOD
        IF (KBOD(JB) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'CBOD decay rate [KBOD=',             KBOD(JB),'] < 0.0 for CBOD group ',JB
        END IF
        IF (TBOD(JB) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'CBOD temperature coefficient [TBOD=',TBOD(JB),'] < 0.0 for CBOD group ',JB
        END IF
        IF (RBOD(JB) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'CBOD5/ultimate CBOD ratio [RBOD=',   RBOD(JB),'] < 0.0 for CBOD group ',JB
        END IF
        IF (CBODS(JB) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'CBOD settling rate [CBODS=',   CBODS(JB),'] < 0.0 for CBOD group ',JB
        END IF
      END DO
    END IF

!** Phosphorus

    WRITE (WIN,*) '    phosphorus'
    IF (PHOSPHORUS) THEN
      DO JW=1,NWB
        IF (PO4R(JW) > 0.05) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Phosphorus sediment release rate [PO4R=',PO4R(JW),'] > 0.05 for waterbody ',JW
        ELSE IF (PO4R(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Phosphorus sediment release rate [PO4R=',PO4R(JW),'] < 0.0 for waterbody ', JW
        END IF
        IF (PARTP(JW) > 2.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Phosphorus/suspended solids partitioning coefficient [PARTP=',PARTP(JW),'] > 2.0 for waterbody ',JW
        ELSE IF (PARTP(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Phosphorus/suspended solids partitioning coefficient [PARTP=',PARTP(JW),'] < 0.0 for waterbody ',JW
        END IF
      END DO
      DO J=1,NAAL
        JA = ALAC(J)
        IF (AHSP(JA) > 0.05) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Algal phosphorus half-saturation concentration [AHSP=',AHSP(JA),'] > 0.05 for algal group ',JA
        ELSE IF (AHSP(JA) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Algal phosphorus half-saturation concentration [AHSP=',AHSP(JA),'] < 0.00 for algal group ',JA
        END IF
      END DO
    END IF

!** Ammonium

    WRITE (WIN,*) '    ammonium'
    IF (AMMONIUM) THEN
      DO JW=1,NWB
        IF (NH4R(JW) > 0.5) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Ammonium sediment release rate [NH4R=',NH4R(JW),'] > 0.5 for waterbody ',JW
        ELSE IF (NH4R(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Ammonium sediment release rate [NH4R=',NH4R(JW),'] < 0.0 for waterbody ',JW
        END IF
        IF (NH4DK(JW) > 0.5) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Nitrification rate [NH4DK=',NH4DK(JW),'] > 0.5 for waterbody ',JW
        ELSE IF (NH4DK(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Nitrification rate [NH4DK=',NH4DK(JW),'] < 0.0 for waterbody ',JW
        END IF
      END DO
      DO J=1,NAAL
        JA = ALAC(JA)
        IF (AHSN(JA) > 0.05) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Nitrogen half-saturation concentration for algal uptake [AHSN=',AHSN(JA),'] > 0.05 for algae group ',JA
        ELSE IF (AHSN(JA) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Nitrogen half-saturation concentration for algal uptake [AHSN=',AHSN(JA),'] < 0.00 for algae group ',JA
        ELSEIF (AHSN(JA) < 0.005)THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Nitrogen half-saturation concentration for algal uptake [AHSN=',AHSN(JA),'] < 0.005 for algae group ',JA
        END IF
      END DO

!**** Nitrification temperature rate multipliers

      DO JW=1,NWB
        IF (NH4T1(JW) <= 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Nitrification minimum temperature [NH4T1=',NH4T1(JW),'] <= 0.0C for waterbody ',JW
        ELSE IF (NH4T1(JW) > 10.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Nitrification minimum temperature [NH4T1=',NH4T1(JW),'] > 10.0C for waterbody ',JW
        END IF
        IF (NH4T2(JW) > 40.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Nitrification maximum temperature [NH4T2=',NH4T2(JW),'] > 40.0C for waterbody ',JW
        ELSE IF (NH4T2(JW) < 10.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Nitrification maximum temperature [NH4T2=',NH4T2(JW),'] < 10.0C for waterbody ',JW
        END IF
        IF (NH4K1(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Nitrification minimum temperature rate multiplier [NH4K2=',NH4K1(JW),'] < 0.0 for waterbody ',JW
        ELSE IF (NH4K1(JW) > 1.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Nitrification minimum temperature rate multiplier [NH4K2=',NH4K1(JW),'] > 1.0 for waterbody ',JW
        END IF
        IF (NH4K2(JW) > 1.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Nitrification maximum temperature rate multiplier [NH4K2=',NH4K2(JW),'] > 1.0 for waterbody ',JW
        ELSE IF (NH4K2(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Nitrification maximum temperature rate multiplier [NH4K2=',NH4K2(JW),'] < 0.0 for waterbody ',JW
        END IF
      END DO
    END IF

!** Nitrate

    WRITE (WIN,*) '    nitrate'
    DO JW=1,NWB
      IF (NITRATE) THEN
        IF (NO3DK(JW) > 0.5) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Denitrification rate [NO3DK=',NO3DK(JW),'] > 0.5 for waterbody ',JW
        ELSE IF (NO3DK(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Denitrification rate [NO3DK=',NO3DK(JW),'] < 0.0 for waterbody ',JW
        END IF
        IF (NO3S(JW) > 0.5) THEN
          CALL WARNINGS
          WRITE (WRN,"('Rate of diffusion of NO3-N into bottom muds [NO3S=',f8.4,'] > 0.5 m/d for waterbody ',i4)")NO3S(JW),JW
        ELSE IF (NO3S(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,"('Rate of diffusion of NO3-N into bottom muds [NO3S=',f8.4,']  < 0.0 for waterbody ',i4)")NO3S(JW),JW
        END IF
        IF (FNO3SED(JW) > 0.0.and. NO3S(JW)>0.0) THEN
          CALL WARNINGS
          WRITE (WRN,"('NO3 is being diffused into sediments and incorporated into sediment organic matter for waterbody ',i4)")JW
          WRITE (WRN,"('  [This may have unintended consequences for OM:SED-N ratio and unexpected high release of NH4 from sediments as SED-N increases]')")
        ENDIF

        IF (FNO3SED(JW) > 1.0) THEN
          CALL ERRORS
          WRITE (ERR,"('Fraction[0-1] of NO3 diffused into sediments incorporated into OM [FNO3SED=',f8.4,']  > 1.0 for waterbody ',i4)")FNO3SED(JW),JW
        ELSE IF (FNO3SED(JW)> 0.5) THEN
          CALL WARNINGS
          WRITE (WRN,"('Fraction of NO3 diffused into sediments incorporated into OM [FNO3SED=',f8.4,'] > 0.5 for waterbody ',i4)")FNO3SED(JW),JW
        ELSE IF (FNO3SED(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,"('Fraction[0-1] of NO3 diffused into sediments incorporated into OM [FNO3SED=',f8.4,']  < 0.0 for waterbody ',i4)")FNO3SED(JW),JW
        END IF

        IF (NO3T1(JW) <= 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Denitrification minimum temperature [NO3T1=',NO3T1(JW),'] <= 0.0 C for waterbody ',JW
        END IF
        IF (NO3T2(JW) > 40.0) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Denitrification maximum temperature [NO3T2=',NO3T2(JW),'] > 40.0C for waterbody ',JW
        ELSE IF (NO3T2(JW) <= 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Denitrification maximum temperature [NO3T2=',NO3T2(JW),'] <= 0.0C for waterbody ', JW
        END IF
        IF (NO3K1(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Denitrification minimum temperature rate multiplier [NO3K1=',NO3K1(JW),'] < 0.0 for waterbody ',JW
        END IF
        IF (NO3K2(JW) > 1.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Denitrification maximum temperature rate multiplier [NO3K1=',NO3K1(JW),'] > 1.0 for waterbody ',JW
        ELSE IF (NO3K2(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Denitrification maximum temperature rate multiplier [NO3K1=',NO3K1(JW),'] < 0.0 for waterbody ',JW
        END IF
      END IF
    END DO

!** Carbon dioxide sediment release

    WRITE (WIN,*) '    carbon dioxide'
    DO JW=1,NWB
      IF (TOT_IC) THEN
        IF (CO2R(JW) > 1.4) THEN
          CALL WARNINGS
          WRITE (WRN,FMTFI) 'Carbon dioxide sediment release rate [CO2R=',CO2R(JW),'] > 1.4 for waterbody ',JW
        ELSE IF (CO2R(JW) < 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTFI) 'Carbon dioxide sediment release rate [CO2R=',CO2R(JW),'] < 0.0 for waterbody ',JW
        END IF
      END IF
    END DO

!!** Iron
!
!    WRITE (WIN,*) '    iron'
!    DO JW=1,NWB
!      IF (IRON) THEN
!        IF (FER(JW) > 2.0) THEN
!          CALL WARNINGS
!          WRITE (WRN,FMTFI) 'Iron sediment release rate [FER=',FER(JW),'] > 2.0 for waterbody ',JW
!        ELSE IF (FER(JW) < 0.0) THEN
!          CALL ERRORS
!          WRITE (ERR,FMTFI) 'Iron sediment release rate [FER=',FER(JW),'] < 0.0 for waterbody ',JW
!        END IF
!        IF (FES(JW) > 10.0) THEN
!          CALL WARNINGS
!          WRITE (WRN,FMTFI) 'Iron settling rate [FES=',       FES(JW),'] > 10.0 for waterbody ',JW
!        ELSE IF (FES(JW) < 0.0) THEN
!          CALL ERRORS
!          WRITE (ERR,FMTFI) 'Iron settling rate [FES=',       FES(JW),'] < 0.0 for waterbody ', JW
!        END IF
!      END IF
!    END DO

!** Stoichiometry

    WRITE (WIN,*) '    stoichiometry'
    DO JW=1,NWB
      IF (O2NH4(JW) /= 4.57) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Oxygen to ammonium stoichiometry [O2NH4=',       O2NH4(JW),'] /= 4.57 for waterbody ', JW
      END IF
      IF (O2OM(JW) /= 1.4) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Oxygen to organic matter stoichiometry [O2OM=',  O2OM(JW),'] /= 1.4 for waterbody ',   JW
      END IF
      IF (ORGP(JW) /= 0.005) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Phosphorus/organic matter stoichiometry [ORGP=',  ORGP(JW),'] /= 0.005 for waterbody ',JW
      END IF
      IF (ORGN(JW) /= 0.08) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Nitrogen/organic matter stoichiometry [ORGN=',    ORGN(JW),'] /= 0.08 for waterbody ', JW
      END IF
      IF (ORGC(JW) /= 0.45) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Carbon/organic matter stoichiometry [ORGC=',      ORGC(JW),'] /= 0.45 for waterbody ', JW
      END IF
      
      IF (O2NH4(JW) < 1.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Oxygen to ammonium stoichiometry [O2NH4=',       O2NH4(JW),'] < 1.0 [DEFAULT=4.2-4.57] for waterbody ', JW
      END IF
      IF (O2OM(JW) < 0.25) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Oxygen to organic matter stoichiometry [O2OM=',  O2OM(JW),'] < 0.25 [DEFAULT=1.4] for waterbody ',   JW
      END IF
      IF (ORGP(JW) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Phosphorus/organic matter stoichiometry [ORGP=',  ORGP(JW),'] <= 0.0 for waterbody ',JW
      END IF
      IF (ORGN(JW) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Nitrogen/organic matter stoichiometry [ORGN=',    ORGN(JW),'] <= 0.0 for waterbody ', JW
      END IF
      IF (ORGC(JW) <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Carbon/organic matter stoichiometry [ORGC=',      ORGC(JW),'] <= 0.0 for waterbody ', JW
      END IF
    END DO
    DO JA=1,NAL
      IF (O2AR(JA) /= 1.1) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Oxygen to algal respiration stoichiometry [O2AR=',O2AR(JA),'] /= 1.1 for algal group', JA
      END IF
    IF (O2AR(JA) < 0.25) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Oxygen to algal respiration stoichiometry [O2AR=',O2AR(JA),'] < 0.25 [DEFAULT=1.1] for algal group', JA
      END IF
      IF (O2AG(JA) /= 1.4) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Oxygen to algal production stoichiometry [O2AG=', O2AG(JA),'] /= 1.4 for algal group', JA
      END IF
    IF (O2AG(JA) < 0.25) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Oxygen to algal production stoichiometry [O2AG=', O2AG(JA),'] < 0.25 [DEFAULT=1.4] for algal group', JA
      END IF
    END DO
    DO JE=1,NEP
     IF (O2ER(JE) < 0.2) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Oxygen to epiphyton respiration stoichiometry [O2ER=',O2ER(JE),'] < 0.2 [DEFAULT=1.1] for epiphyton group', JE
      END IF
      IF (O2ER(JE) /= 1.1) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Oxygen to epiphyton respiration stoichiometry [O2ER=',O2ER(JE),'] /= 1.1 for epiphyton group', JE
      END IF
      IF (O2EG(JE) /= 1.4) THEN                                                                                        !CB 04/08/03
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Oxygen to epiphyton production stoichiometry [O2EG=', O2EG(JE),'] /= 1.4 for epiphyton group', JE
      END IF
      IF (O2EG(JE) < 0.2) THEN                                                                                        !CB 04/08/03
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Oxygen to epiphyton production stoichiometry [O2EG=', O2EG(JE),'] < 0.2 [DEFAULT=1.4] for epiphyton group', JE
      END IF
    END DO

!** Oxygen limit

    WRITE (WIN,*) '    DO half saturation coeff KDO'
    IF (DISSOLVED_OXYGEN) THEN
      IF (O2LIM > 2.0) THEN
        CALL WARNINGS
        WRITE (WRN,"('DO half-saturation coeff for anoxia [KDO=',f8.3,'] > 2.0 mg/l')")o2lim
      ELSE IF (O2LIM <= 0.0) THEN
        CALL ERRORS
        WRITE (ERR,"('DO half-saturation coeff for anoxia [KDO=',f8.3,'] <= 0.0')")o2lim
      END IF
    END IF


! Sediment oxygen demand

  WRITE (WIN,*) '    sediment oxygen demand'
  IF (DISSOLVED_OXYGEN) THEN
    DO I=1,IMX
      IF (SOD(I) > 5.0) THEN
        CALL WARNINGS
        WRITE (WRN,FMTFI) 'Sediment oxygen demand [SOD=',SOD(I),'] > 5.0 for segment ',I
      ELSE IF (SOD(I) < 0.0) THEN
        CALL ERRORS
        WRITE (ERR,FMTFI) 'Sediment oxygen demand [SOD=',SOD(I),'] < 0.0 for segment ',I
      END IF
    END DO
  END IF

! Reaeration formulae

  WRITE (WIN,*) '    reaeration'
  DO JW=1,NWB
    IF (REAERC(JW) == '   RIVER') THEN
      DO JB=BS(JW),BE(JW)
        IF (SLOPE(JB) == 0) THEN
          IF (NEQN(JW) == 3 .OR. NEQN(JW) == 7 .OR. NEQN(JW) == 8 .OR. NEQN(JW) == 9) THEN
            CALL ERRORS
            WRITE (ERR,FMTI) 'Reaeration equation # [NEQN=',NEQN(JW),'] is inappropriate for waterbody ',JW
          END IF
        END IF
      END DO
      IF (NEQN(JW) < 0 .OR. NEQN(JW) > 10) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Reaeration equation # [NEQN=',NEQN(JW),'] is inappropriate for waterbody ',JW
      END IF
      IF (NEQN(JW) == 9) THEN
        IF (RCOEF1(JW) == 0.0 .AND. RCOEF2(JW) == 0.0 .AND. RCOEF3(JW) == 0.0 .AND. RCOEF4(JW) == 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTI) 'Reaeration equation # [NEQN=',NEQN(JW),'] is inappropriate for waterbody ',JW
        END IF
      END IF
    ELSE IF (REAERC(JW) == '    LAKE') THEN
      IF (NEQN(JW) < 1 .OR. NEQN(JW) > 14) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Reaeration equation # [NEQN=',NEQN(JW),'] is inappropriate for waterbody ',JW
      END IF
      IF (NEQN(JW) == 14) THEN
        IF (RCOEF1(JW) == 0.0 .AND. RCOEF2(JW) == 0.0 .AND. RCOEF3(JW) == 0.0 .AND. RCOEF4(JW) == 0.0) THEN
          CALL ERRORS
          WRITE (ERR,FMTI) 'Reaeration equation # [NEQN=',NEQN(JW),'] is inappropriate for waterbody ',JW
        END IF
      END IF
    ELSE IF (REAERC(JW) == ' ESTUARY') THEN
      IF (NEQN(JW) < 0 .OR. NEQN(JW) > 2) THEN
        CALL ERRORS
        WRITE (ERR,FMTI) 'Reaeration equation # [NEQN=',NEQN(JW),'] is inappropriate for waterbody ',JW
      END IF
    ELSE
      CALL ERRORS
      WRITE (ERR,*) JW, CONFN
      WRITE (ERR,FMTI) 'Reaeration type [REAERC='//REAERC(JW)//'] /= "RIVER", "LAKE" or "ESTUARY" for waterbody ',JW
    END IF
    IF(MINKL(JW) > 0.8)THEN
        CALL WARNINGS
        IF (REAERC(JW) == '    LAKE') THEN
        WRITE (WRN,FMTFI) 'Minimum reaeration coefficient is  [MINKL=',MINKL(JW),'] > 0.8 m/d for waterbody ',JW
        ELSE
        WRITE (WRN,FMTFI) 'Minimum reaeration coefficient is  [MINKL=',MINKL(JW),'] > 0.8 1/d for waterbody ',JW
        ENDIF
    ENDIF
    IF(CAC(NDGP)=='      ON')THEN
        IF(DGPO2(JW) > 1.04)THEN
           CALL WARNINGS
           WRITE (WRN,FMTFI) 'Ratio of dissolved gas to O2 is  [DGPO2=',DGPO2(JW),'] > 1.04 for waterbody ',JW
        ENDIF
        IF(DGPO2(JW) < 1.01)THEN
           CALL WARNINGS
           WRITE (WRN,FMTFI) 'Ratio of dissolved gas to O2 is  [DGPO2=',DGPO2(JW),'] < 1.01 for waterbody ',JW
        ENDIF
    ENDIF 
    
  END DO
  
  ENDIF
  

  ! Filenames

  WRITE (WIN,*) 'Filenames'
  DO JW=1,NWB
    BTHFN1 = BTHFN(JW); METFN1 = METFN(JW); SNPFN1 = SNPFN(JW); PRFFN1 = PRFFN(JW); SPRFN1 = SPRFN(JW); VPLFN1 = VPLFN(JW)
    CPLFN1 = CPLFN(JW); FLXFN1 = FLXFN(JW)
    DO JJW=1,NWB
      IF (JJW /= JW .AND. BTHFN1 == BTHFN(JJW)) CALL DUPLICATE_FILENAME (ERR,BTHFN1)
      IF (JJW /= JW .AND. METFN1 == METFN(JJW)) CALL DUPLICATE_FILENAME (ERR,METFN1)
      IF (JJW /= JW .AND. SNPFN1 == SNPFN(JJW)) CALL DUPLICATE_FILENAME (ERR,SNPFN1)
      IF (JJW /= JW .AND. PRFFN1 == PRFFN(JJW)) CALL DUPLICATE_FILENAME (ERR,PRFFN1)
      IF (JJW /= JW .AND. SPRFN1 == SPRFN(JJW)) CALL DUPLICATE_FILENAME (ERR,SPRFN1)
     ! IF (JJW /= JW .AND. VPLFN1 == VPLFN(JJW)) CALL DUPLICATE_FILENAME (ERR,VPLFN1)   ! not needed since only 1 anyway
      IF (JJW /= JW .AND. CPLFN1 == CPLFN(JJW)) CALL DUPLICATE_FILENAME (ERR,CPLFN1)
      IF (JJW /= JW .AND. FLXFN1 == FLXFN(JJW)) CALL DUPLICATE_FILENAME (ERR,FLXFN1)
      IF (PRECIPITATION(JW)) THEN
        IF (JJW /= JW .AND. PREFN1 == PREFN(JJW)) CALL DUPLICATE_FILENAME (ERR,PREFN1)
        IF (JJW /= JW .AND. TPRFN1 == TPRFN(JJW)) CALL DUPLICATE_FILENAME (ERR,TPRFN1)
        IF (CONSTITUENTS) THEN
          IF (JJW /= JW .AND. CPRFN1 == CPRFN(JJW)) CALL DUPLICATE_FILENAME (ERR,CPRFN1)
        END IF
      END IF
    END DO
  END DO
  DO JB=1,NBR
    QINFN1 = QINFN(JB); TINFN1 = TINFN(JB); CINFN1 = CINFN(JB); QOTFN1 = QOTFN(JB); QDTFN1 = QDTFN(JB); TDTFN1 = TDTFN(JB)
    CDTFN1 = CDTFN(JB); PREFN1 = PREFN(JB); TPRFN1 = TPRFN(JB); CPRFN1 = CPRFN(JB); EUHFN1 = EUHFN(JB); TUHFN1 = TUHFN(JB)
    CUHFN1 = CUHFN(JB); EDHFN1 = EDHFN(JB); TDHFN1 = TDHFN(JB); CDHFN1 = CDHFN(JB)
    DO JJB=1,NBR
      IF (JJB /= JB .AND. QINFN1 == QINFN(JJB) .AND. UQ_EXTERNAL(JB))     CALL DUPLICATE_FILENAME (ERR,QINFN1)
      IF (JJB /= JB .AND. TINFN1 == TINFN(JJB) .AND. UQ_EXTERNAL(JB))     CALL DUPLICATE_FILENAME (ERR,TINFN1)
      IF (JJB /= JB .AND. QOTFN1 == QOTFN(JJB) .AND. DQ_EXTERNAL(JB))     CALL DUPLICATE_FILENAME (ERR,QOTFN1)
      IF (JJB /= JB .AND. QDTFN1 == QDTFN(JJB) .AND. DIST_TRIB(JB))       CALL DUPLICATE_FILENAME (ERR,QDTFN1)
      IF (JJB /= JB .AND. TDTFN1 == TDTFN(JJB) .AND. DIST_TRIB(JB))       CALL DUPLICATE_FILENAME (ERR,TDTFN1)
      IF (JJB /= JB .AND. EUHFN1 == EUHFN(JJB) .AND. UH_EXTERNAL(JB))     CALL DUPLICATE_FILENAME (ERR,EUHFN1)
      IF (JJB /= JB .AND. TUHFN1 == TUHFN(JJB) .AND. UH_EXTERNAL(JB))     CALL DUPLICATE_FILENAME (ERR,TUHFN1)
      IF (JJB /= JB .AND. EDHFN1 == EDHFN(JJB) .AND. DH_EXTERNAL(JB))     CALL DUPLICATE_FILENAME (ERR,EDHFN1)
      IF (JJB /= JB .AND. TDHFN1 == TDHFN(JJB) .AND. DH_EXTERNAL(JB))     CALL DUPLICATE_FILENAME (ERR,TDHFN1)
      IF (CONSTITUENTS) THEN
        IF (JJB /= JB .AND. CINFN1 == CINFN(JJB) .AND. UQ_EXTERNAL(JB))   CALL DUPLICATE_FILENAME (ERR,CINFN1)
        IF (JJB /= JB .AND. CDTFN1 == CDTFN(JJB) .AND. DIST_TRIB(JB))     CALL DUPLICATE_FILENAME (ERR,CDTFN1)
!        IF (JJB /= JB .AND. CPRFN1 == CPRFN(JJB) .AND. PRECIPITATION(JB)) CALL DUPLICATE_FILENAME (ERR,CPRFN1)   ! SW 1/16/04
        IF (JJB /= JB .AND. CUHFN1 == CUHFN(JJB) .AND. UH_EXTERNAL(JB))   CALL DUPLICATE_FILENAME (ERR,CUHFN1)
        IF (JJB /= JB .AND. CDHFN1 == CDHFN(JJB) .AND. DH_EXTERNAL(JB))   CALL DUPLICATE_FILENAME (ERR,CDHFN1)
      END IF
    END DO
  END DO
  DO JT=1,NTR
    QTRFN1 = QTRFN(JT); TTRFN1 = TTRFN(JT); CTRFN1 = CTRFN(JT)
    DO JJT=1,NTR
      IF (JJT /= JT .AND. QTRFN1 == QTRFN(JJT)) CALL DUPLICATE_FILENAME (ERR,QTRFN1)
      IF (JJT /= JT .AND. TTRFN1 == TTRFN(JJT)) CALL DUPLICATE_FILENAME (ERR,TTRFN1)
      IF (JJT /= JT .AND. CTRFN1 == CTRFN(JJT) .AND. CONSTITUENTS) CALL DUPLICATE_FILENAME (ERR,CTRFN1)
    END DO
  END DO

!***********************************************************************************************************************************
!*                                            Task 5:  Calculate Volume-Area-Elevation Table                                      **
!***********************************************************************************************************************************

  WRITE (WIN,*) 'Volume-area-elevation'
  DO JW=1,NWB
    VOLG(JW) = 0.0
    DO JB=BS(JW),BE(JW)
      IU  = US(JB)
      ID  = DS(JB)
      IF (.NOT.DAM_FLOW(JB)) THEN
        IF (UHS(JB) /= 0 .AND. UHS(JB) > 1) THEN
          DIFF = ABS(ELWS(IU)-ELWS(UHS(JB)))
          IF (DIFF > 1.5) THEN
            CALL WARNINGS
            WRITE (WRN,FMTF2I) 'Water surface elevation difference of ',DIFF,' m between segment ',IU,' and segment ',UHS(JB)
          END IF
        END IF
      END IF
      IF (DHS(JB) /= 0 .AND. DHS(JB) > 1) THEN
        DIFF = ABS(ELWS(ID)-ELWS(DHS(JB)))
        IF (DIFF > 1.5) THEN
          CALL WARNINGS
          WRITE (WRN,FMTF2I) 'Water surface elevation difference of ',DIFF,' m between segment ',ID,' and segment ', DHS(JB)
        END IF
      END IF
      DO I=ID,IU,-1
        XI(I)   = XI(I+1)+(DLX(I)+DLX(I+1))*0.5
        XBR(JB) = XBR(JB)+DLX(I)
      END DO
      XGR(JW) = XGR(JW)+XBR(JB)

!**** Water surface and bottom layers

      DO I=IU-1,ID+1
        IF (ELWS(I) <= EL(KB(I)+1,I)) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Water surface elevation is below bottom elevation at segment ',I
        ELSE IF (ELWS(I) <=  EL(KB(I)+1,I)+0.4*H(KB(I),JW)) THEN
          CALL WARNINGS
          WRITE (WRN,FMTI) 'Water surface elevation is close to bottom elevation at segment ',I
        END IF
        IF (I < ID+1) THEN
          DIFF = ABS(ELWS(I)-ELWS(I+1))
          IF (DIFF > 1.5) THEN
            CALL WARNINGS
            WRITE (WRN,FMTF2I) 'Water surface elevation difference of ',DIFF,'m between segment ',I,'and segment ',I+1
          END IF
        END IF
      END DO

!**** Upstream active segment and single layer

      IUT = IU
      DO I=IU,ID
        IF (KB(I)-KTWB(JW) < NL(JB)-1) IUT = I+1
      END DO
      IF (IUT > DS(JB)-1) THEN
        CALL ERRORS
        WRITE (ERR,'(A,i3,A,i3,a,i3,a,6(i3,1x))') 'No active segments in branch [JB=',JB,'] of waterbody [JW=',jw,'] KTWB=',ktwb(jw),' KB for last 6 segments of Branch:',(kb(i),i=id-5,id)
      END IF
      CUS(JB) = IUT
      IF (CUS(JB) /= US(JB)) THEN
        DO JBB=1,NBR
          IF (DHS(JBB) == US(JB)) THEN
            CALL ERRORS
            WRITE (ERR,'(A,i3,A,i3,a,i3,a,i3,a,i3,a)') 'Water level too low in branch [JB=',JB,'] of waterbody [JW=',jw,'] where [CUS(JB)=',CUS(JB),'] connected to branch [JB=',JBB,'] by downstream head. CUS must equal US [US(JB)=',us(jb),'] since used as a boundary head.'
          END IF
          IF (UHS(JBB) == US(JB)) THEN
            CALL ERRORS
            WRITE (ERR,'(A,i3,A,i3,a,i3,a,i3,a,i3,a)') 'Water level too low in branch [JB=',JB,'] of waterbody [JW=',jw,'] where [CUS(JB)=',CUS(JB),'] connected to branch [JB=',JBB,'] by upstream head. CUS must equal US [US(JB)=',us(jb),'] since used as a boundary head.'
          END IF
        END DO
      END IF

!**** Branch and grid total volume

      DO I=IU,ID
        VOLB(JB) = VOLB(JB)+DLX(I)*B(KTWB(JW),I)*(H(KTWB(JW),JW)-Z(I))
        VOLG(JW) = VOLG(JW)+DLX(I)*B(KTWB(JW),I)*(H(KTWB(JW),JW)-Z(I))
        DO K=KTWB(JW)+1,KB(I)
          VOLB(JB) = VOLB(JB)+DLX(I)*B(K,I)*H(K,JW)
          VOLG(JW) = VOLG(JW)+DLX(I)*B(K,I)*H(K,JW)
        END DO
      END DO

!**** Branch and grid area and volume by layer

      DO K=KMX-1,2,-1
        NCCBR(K,JB) = NCCBR(K+1,JB)
        CVBR(K,JB)  = CVBR(K+1,JB)
        DO I=IU,ID
          IF (K <= KB(I)) THEN
            SABR(K,JB)  = SABR(K,JB)+B(K,I)        *DLX(I)
            CVBR(K,JB)  = CVBR(K,JB)+B(K,I)*H(K,JW)*DLX(I)
            NCCBR(K,JB) = NCCBR(K,JB)+1
          END IF
        END DO
        SAGR(K,JW)  = SAGR(K,JW) +SABR(K,JB)
        NCCGR(K,JW) = NCCGR(K,JW)+NCCBR(K,JB)
        CVGR(K,JW)  = CVGR(K,JW) +CVBR(K,JB)
      END DO
    END DO
  END DO

! Average depths and widths

  DO K=KMX-1,2,-1
    DO JW=1,NWB
      IF (SAGR(K,JW) /= 0.0) HGR(K,JW) = CVGR(K,JW)/SAGR(K,JW)
      BGR(K,JW) = SAGR(K,JW)/XGR(JW)
      DO JB=BS(JW),BE(JW)
        IF (SABR(K,JB) /= 0.0) HBR(K,JB) = CVBR(K,JB)/SABR(K,JB)
        BBR(K,JB) = SABR(K,JB)/XBR(JB)
      END DO
    END DO
  END DO

! Beginning and ending segment and bottom layer for snapshots




  DO JW=1,NWB
    IEPR(JW) = NISNP(JW)

    DO I=1,IEPR(JW)
        ACTIVE_SEGMENT = .FALSE.
                IF (ISNP(I,JW) >= US(1) .AND. ISNP(I,JW) <= DS(NBR)) ACTIVE_SEGMENT = .TRUE.
                IF(.NOT.ACTIVE_SEGMENT)THEN
                CALL ERRORS
                WRITE (ERR,'(a,i4,a,i4)') 'Snapshot segment error: ISNP(IEPR,JW) is not set to an active segment. ISNP=', ISNP(I,JW), ' FOR JW=',JW
                IPR(I,JW) = (US(1)+DS(NBR))/2     ! SET TO FICTICIOUS LOCATION SO ARRAY LIMITS NOT VIOLATED IN THE NEXT LINES
                WRITE(WIN,*)'Error in ISNP alters SNP.OPT initial condition output - see pre.err'

                ELSE
                IPR(I,JW) = ISNP(I,JW)

                ENDIF
    END DO
    DO I=1,IEPR(JW)
      KEPR(JW) = MAX(KB(IPR(I,JW)),KEPR(JW))
    END DO

  END DO

!***********************************************************************************************************************************
!*                                           Task 6:  Initial Conditions for Simulation                                           **
!***********************************************************************************************************************************

! Temperature and constituents

  WRITE (WIN,*) 'Initial conditions'
  DO JW=1,NWB
    IF (OPEN_VPR(JW)) THEN
      WRITE (WIN,*) '        vpr file input:', adjustl(trim(vprfn(jw)))
      OPEN (VPR(JW),FILE=VPRFN(JW),STATUS='OLD',IOSTAT=IERR)
      IF (IERR == 0) THEN
      UNIT=VPR(JW)
      READ(VPR(JW),'(A1)')ICHAR1
      IF(ICHAR1=='$')THEN
        READ (VPR(JW),'(/)')
        IF (VERT_TEMP(JW)) READ (VPR(JW),*,ERR=400)        AID, (TVP(K,JW),   K=KTWB(JW),KBMAX(JW))
        WRITE (WIN,*) '        vpr file input: FINISHED TEMPERATURE next Constituents'
        IF (CONSTITUENTS) THEN
          DO JC=1,NCT
            IF (VERT_CONC(JC,JW)) READ (VPR(JW),*,ERR=400) AID, (CVP(K,JC,JW),K=KTWB(JW),KBMAX(JW))
          END DO
        END IF
        WRITE (WIN,*) '        vpr file input: FINISHED CONSTITUENTS next Epiphyton and Macrophytes'
        DO J=1,NEP   ! EPIPHYTON
        IF(VERT_EPIPHYTON(JW,J))READ (VPR(JW),*,ERR=400) AID, (EPIVP(K,JW,J),K=KTWB(JW),KBMAX(JW))
        ENDDO
        DO J=1,NMC  ! MACROPHYTES
        IF(VERT_MACRO(JW,J))READ (VPR(JW),*,ERR=400) AID, (MACRCVP(K,JW,J),K=KTWB(JW),KBMAX(JW))
        ENDDO
      ELSE
        IF (VERT_TEMP(JW)) READ (VPR(JW),'(/A8/(8X,9F8.0))',ERR=400)        AID, (TVP(K,JW),   K=KTWB(JW),KBMAX(JW))
               WRITE (WIN,*) '        vpr file input: FINISHED TEMPERATURE next Constituents'
        IF (CONSTITUENTS) THEN
          DO JC=1,NCT
            IF (VERT_CONC(JC,JW)) READ (VPR(JW),'(/A8/(8X,9F8.0))',ERR=400) AID, (CVP(K,JC,JW),K=KTWB(JW),KBMAX(JW))
          END DO
        END IF
               WRITE (WIN,*) '        vpr file input: FINISHED CONSTITUENTS next Epiphyton and Macrophytes'
        DO J=1,NEP   ! EPIPHYTON
        IF(VERT_EPIPHYTON(JW,J))READ (VPR(JW),'(/A8/(8X,9F8.0))',ERR=400) AID, (EPIVP(K,JW,J),K=KTWB(JW),KBMAX(JW))
        ENDDO
        DO J=1,NMC  ! MACROPHYTES
        IF(VERT_MACRO(JW,J))READ (VPR(JW),'(/A8/(8X,9F8.0))',ERR=400) AID, (MACRCVP(K,JW,J),K=KTWB(JW),KBMAX(JW))
        ENDDO
        
      ENDIF
      ELSE
        CALL ERRORS
        WRITE (ERR,FMTA) 'Vertical initial conditions for temperature or constituents or epiphyton or macrophytes, but could not open '//VPRFN(JW)
      END IF
    END IF

!** Longitudinal/vertical initial profiles

    IF (OPEN_LPR(JW)) THEN
      OPEN (LPR(JW),FILE=LPRFN(JW),STATUS='OLD',IOSTAT=IERR)
      WRITE (WIN,*) '        lpr file input:', adjustl(trim(lprfn(jw)))
      IF (IERR == 0) THEN
        UNIT=LPR(JW)
        READ (LPR(JW),'(A1)')ICHAR1
        IF(ICHAR1=='$')READ( LPR(JW),*)AID
! Read LPR File - later in code - see below
        
      ELSE
        CALL ERRORS
        WRITE (ERR,FMTA) 'Longitudinal/vertical initial conditions for temperature or constituents or epihpyton or macrophytes, but could not open '//LPRFN(JW)
        WRITE(BRA,'(I0)') JW
        NAME = 'lpr_wb'//TRIM(ADJUSTL(BRA))//'.npt'
        OPEN  (LPR(JW),FILE=NAME,STATUS='UNKNOWN')
        WRITE (LPR(JW),FMTA) 'Longitudinal/vertical initial constituents concentrations'
        IF (LONG_TEMP(JW)) THEN
          DO JB=BS(JW),BE(JW)
            DO I=CUS(JB),DS(JB)
              WRITE(LPR(JW),'(/10A,3X,A,I0/(:8X,9F8.3))') 'TEMP    ',('     T2I',J=1,9),'Segment ',I,(0.0,K=KTWB(JW),KB(I))
            END DO
          END DO
        END IF
     IF(CONSTITUENTS)THEN
        DO JC=1,NCT
          IF (LONG_CONC(JC,JW)) THEN
            DO JB=BS(JW),BE(JW)
              DO I=CUS(JB),DS(JB)
                WRITE(LPR(JW),'(/10A,3X,A,1X,I0/(:8X,9F8.3))') CNAME2(JC),('     C2I',J=1,9),'Segment',I,(0.0,K=KTWB(JW),KB(I))
              END DO
            END DO
          END IF
        END DO
        IF (LONG_SED(JW)) THEN
          DO JB=BS(JW),BE(JW)
            DO I=CUS(JB),DS(JB)
              WRITE(LPR(JW),'(/10A,3X,A,I0/(:8X,9F8.3))') 'SED     ',('   SEDCI',J=1,9),'Segment ',I,(0.0,K=KTWB(JW),KB(I))
            END DO
          END DO
        END IF
        END IF
      END IF
    END IF

    IF (.NOT.OPEN_LPR(JW).AND.CONSTITUENTS.AND.OPEN_LPRC(JW)) THEN
      OPEN (LPR(JW),FILE=LPRFN(JW),STATUS='OLD',IOSTAT=IERR)
      WRITE (WIN,*) '        lpr file input:', adjustl(trim(lprfn(jw)))
      IF (IERR == 0) THEN
          UNIT=LPR(JW)
          READ (LPR(JW),'(A1)')ICHAR1
        IF(ICHAR1=='$')READ( LPR(JW),*)AID
      ELSE
        CALL ERRORS
        WRITE (ERR,FMTA) 'Longitudinal/vertical initial conditions for constituents, but could not open '//LPRFN(JW)
        WRITE(BRA,'(I0)') JW
        NAME = 'lpr_wb'//TRIM(ADJUSTL(BRA))//'.npt'
        OPEN  (LPR(JW),FILE=NAME,STATUS='UNKNOWN')
        WRITE (LPR(JW),FMTA) 'Longitudinal/vertical initial constituents concentrations'
        IF (LONG_TEMP(JW)) THEN
          DO JB=BS(JW),BE(JW)
            DO I=CUS(JB),DS(JB)
              WRITE(LPR(JW),'(/10A,3X,A,I0/(:8X,9F8.3))') 'TEMP    ',('     T2I',J=1,9),'Segment ',I,(0.0,K=KTWB(JW),KB(I))
            END DO
          END DO
        END IF
        DO JC=1,NCT
          IF (LONG_CONC(JC,JW)) THEN
            DO JB=BS(JW),BE(JW)
              DO I=CUS(JB),DS(JB)
                WRITE(LPR(JW),'(/10A,3X,A,1X,I0/(:8X,9F8.3))') CNAME2(JC),('     C2I',J=1,9),'Segment',I,(0.0,K=KTWB(JW),KB(I))
              END DO
            END DO
          END IF
        END DO
        IF (LONG_SED(JW)) THEN
          DO JB=BS(JW),BE(JW)
            DO I=CUS(JB),DS(JB)
              WRITE(LPR(JW),'(/10A,3X,A,I0/(:8X,9F8.3))') 'SED     ',('   SEDCI',J=1,9),'Segment ',I,(0.0,K=KTWB(JW),KB(I))
            END DO
          END DO
        END IF
      END IF
    END IF

    DO JB=BS(JW),BE(JW)

!**** Ice thickness

      ICETH(CUS(JB)-1:DS(JB)+1) = ICEI(JW)

!**** Temperature

    IF(LONG_TEMP(JW).AND.ICHAR1=='$')READ (LPR(JW),*,ERR=400)
      DO I=CUS(JB),DS(JB)
        IF (LONG_TEMP(JW))THEN
            IF(ICHAR1=='$')THEN
            READ (LPR(JW),*,ERR=400) AID,(T2(K,I),K=KTWB(JW),KB(I))
            ELSE
            READ (LPR(JW),'(/A8/(8X,9F8.0))',IOSTAT=IERR,ERR=400) AID, (T2(K,I),K=KTWB(JW),KB(I))
            ENDIF
        ENDIF
        IF (IERR == 0) THEN
          DO K=KTWB(JW),KB(I)
            IF (ISO_TEMP(JW))  T2(K,I) = T2I(JW)
            IF (VERT_TEMP(JW)) T2(K,I) = TVP(K,JW)
          END DO
        END IF
      END DO
    END DO
        IF (LONG_TEMP(JW))  WRITE (WIN,*) '        lpr file input: FINISHED TEMPERATURE next Constituents'

!**** Constituents


      DO JC=1,NAC
          IF (LONG_CONC(CN(JC),JW).AND.ICHAR1=='$')READ (LPR(JW),*,ERR=400)AID
       DO JB=BS(JW),BE(JW)
        DO I=CUS(JB),DS(JB)
          IF (LONG_CONC(CN(JC),JW)) THEN
              IF(ICHAR1=='$')THEN
              READ (LPR(JW),*,ERR=400) AID,(C2(K,I,CN(JC)),K=KTWB(JW),KB(I))
              ELSE
              READ (LPR(JW),'(/A8/(8X,9F8.0))',IOSTAT=IERR,ERR=400) AID,(C2(K,I,CN(JC)),K=KTWB(JW),KB(I))
              ENDIF
          ENDIF
          IF (IERR == 0) THEN
            DO K=KTWB(JW),KB(I)
              IF (ISO_CONC(CN(JC),JW))  C2(K,I,CN(JC)) = C2IWB(CN(JC),JW)
              IF (VERT_CONC(CN(JC),JW)) C2(K,I,CN(JC)) = CVP(K,CN(JC),JW)
            END DO
          END IF
        END DO
      END DO
      END DO
      !**** Epiphyton


    DO JE=1,NEP
        IF (EPIPHYTON_CALC(JW,JE)) THEN
        IF (LONG_EPIPHYTON(JW,JE).AND.ICHAR1=='$')READ (LPR(JW),*,ERR=400)
           DO JB=BS(JW),BE(JW)
           DO I=CUS(JB),DS(JB)
              IF (LONG_EPIPHYTON(JW,JE))THEN
                  IF(ICHAR1=='$')THEN
                    READ (LPR(JW),*,ERR=400) AID,(EPD(K,I,JE),K=KTWB(JW),KB(I))
                      ELSE
                    READ (LPR(JW),'(/A8/(8X,9F8.0))',IOSTAT=IERR,ERR=400) AID,(EPD(K,I,JE),K=KTWB(JW),KB(I))
                ENDIF
              ENDIF
              IF (ISO_EPIPHYTON(JW,JE))  EPD(:,I,JE) = EPICI(JW,JE)
              IF (VERT_EPIPHYTON(JW,JE)) EPD(:,I,JE) = EPIVP(:,JW,JE)
            END DO
           END DO
        END IF
    END DO
    
        DO JE=1,NMC
        IF (MACWBC(JW,JE) == '      ON') THEN
        IF (LONG_MACRO(JW,JE).AND.ICHAR1=='$')READ (LPR(JW),*,ERR=400)
           DO JB=BS(JW),BE(JW)
           DO I=CUS(JB),DS(JB)
              IF (LONG_MACRO(JW,JE))THEN
                  IF(ICHAR1=='$')THEN
                    READ (LPR(JW),*,ERR=400) AID,(MAC(K,I,JE),K=KTWB(JW),KB(I))
                      ELSE
                    READ (LPR(JW),'(/A8/(8X,9F8.0))',IOSTAT=IERR,ERR=400) AID,(MAC(K,I,JE),K=KTWB(JW),KB(I))
                   ENDIF
              ENDIF
              IF (ISO_MACRO(JW,JE))  MAC(:,I,JE) = MACWBCI(JW,JE)
              IF (VERT_MACRO(JW,JE)) MAC(:,I,JE) = MACRCVP(:,JW,JE)
            END DO
           END DO
        END IF
      END DO

!**** Sediments

      DO JB=BS(JW),BE(JW)
        IF (SEDIMENT_CALC(JW)) THEN
            IF(LONG_SEDIMENT(JW).AND.JB==BS(JW))READ (LPR(JW),*,ERR=400)
          DO I=CUS(JB),DS(JB)
            IF (LONG_SEDIMENT(JW))THEN
                IF(ICHAR1=='$')THEN
                    READ (LPR(JW),*,ERR=400)AID, (SED(K,I),K=KTWB(JW),KB(I))      ! SW 6/11/14
                    ELSE
                    READ (LPR(JW),'(/A8/(8X,9F8.0))',IOSTAT=IERR,ERR=400) AID,(SED(K,I),K=KTWB(JW),KB(I))
                ENDIF
            ENDIF
          END DO
        END IF
      END DO

  END DO    ! JW
  
  DO JW=1,NWB
        IF(CDWBC(PH_DER,JW)/='      ON')THEN
            IF(TOT_IC)THEN
            CALL ERRORS
            WRITE(ERR,'(A,I4)')'Since TIC is ON as an active state variable, you must turn on pH as a derived variable for waterbody #:',JW
            ENDIF
        ENDIF
        IF(CDWBC(PH_DER,JW) =='      ON')THEN
            IF(.not.TOT_IC)THEN
            CALL ERRORS
            WRITE(ERR,'(A,I4)')'Since pH is ON as an active derived variable, you must turn on TIC as an active state variable for waterbody #:',JW
            ENDIF
        ENDIF

  ENDDO
  
INQUIRE(FILE='w2_lake_river_contour.csv',EXIST=LAKE_RIVER_CONTOURC)    
IF(LAKE_RIVER_CONTOURC)THEN
    UNIT=LAKE_RIVER_CONTOUR
    AID='w2_lake_river_contour.csv'
  OPEN (LAKE_RIVER_CONTOUR,FILE='w2_lake_river_contour.csv',STATUS='OLD',IOSTAT=I)
  READ(LAKE_RIVER_CONTOUR,*,ERR=400)
  READ(LAKE_RIVER_CONTOUR,*,ERR=400)LAKE_RIVER_CONTOUR_ON
  READ(LAKE_RIVER_CONTOUR,*,ERR=400)
  READ(LAKE_RIVER_CONTOUR,*,ERR=400)NUM_LAKE_CONTOUR,LAKE_CONTOUR_FORMAT
  READ(LAKE_RIVER_CONTOUR,*,ERR=400)
  DO JJ=1,NUM_LAKE_CONTOUR
  READ(LAKE_RIVER_CONTOUR,*,ERR=400)LAKE_CONTOUR_SEG(JJ),LAKE_CONTOUR_START(JJ),LAKE_CONTOUR_FREQ(JJ)
  ENDDO
  READ(LAKE_RIVER_CONTOUR,*,ERR=400)
  READ(LAKE_RIVER_CONTOUR,*,ERR=400)NUM_RIVER_CONTOUR,RIVER_CONTOUR_FORMAT
  READ(LAKE_RIVER_CONTOUR,*,ERR=400)
  DO JJ=1,NUM_RIVER_CONTOUR
  READ(LAKE_RIVER_CONTOUR,*,ERR=400)RIVER_CONTOUR_BR1(JJ),RIVER_CONTOUR_BR2(JJ),RIVER_CONTOUR_START(JJ),RIVER_CONTOUR_FREQ(JJ)
      DO JW=1,NWB
      IF(BS(JW)>= RIVER_CONTOUR_BR1(JJ) .AND. BE(JW) <= RIVER_CONTOUR_BR1(JJ))THEN
             IF(BS(JW)>= RIVER_CONTOUR_BR2(JJ) .AND. BE(JW) <= RIVER_CONTOUR_BR2(JJ))THEN
                 EXIT
             ELSE
                 CALL ERRORS
                 WRITE(ERR,'(A,I2,A,I2,A,I2,A)')'River Temperature/DO Contours Error: All branches must be within 1 waterbody. Branch 1[',RIVER_CONTOUR_BR1(JJ),'] is in waterbody [',JW,'] but Branch2[',RIVER_CONTOUR_BR2(JJ),'] is not in same waterbody.'
                 EXIT
             ENDIF
      ENDIF
      
    ENDDO

  ENDDO
  CLOSE(LAKE_RIVER_CONTOUR)
  ENDIF


!***********************************************************************************************************************************
!*                                                     Task 7:  Output section                                                    **
!***********************************************************************************************************************************

! Initial input and conditions

  WRITE (WIN,*) 'Preprocessor output'
  WRITE (INI,'(1X,A/)') 'CE-QUAL-W2 preprocessor - V4.5'
  WRITE (INI,11130) (TITLE(J),J=2,11)                                                                                  !TC 07/23/02
  WRITE (INI,11140) 'Time Control', TMSTRT, TMEND, YEAR
  WRITE (INI,11150)  NDLT, DLTMIN
  WRITE (INI,11160) (DLTD(J),J=1,NDLT)
  WRITE (INI,11170) (DLTMAX(J),J=1,NDLT)
  WRITE (INI,11180) (DLTF(J),J=1,NDLT)
  WRITE (INI,11190) 'Timestep limitation', (JW,VISC(JW)(6:8), CELC(JW)(6:8), DLTADD(JW)(6:8), JW=1,NWB)
  WRITE (INI,'(/1X,A)') 'Initial Conditions'
  DO JW=1,NWB
    WRITE (INI,FMTI) '   Waterbody ',JW
    IF (ISO_TEMP(JW))  WRITE (INI,'(5X,A25,F5.1,A3)') 'Temperature      [T2I] = ',T2I(JW),'C'
    IF (VERT_TEMP(JW)) WRITE (INI,'(5X,A52)')         'Temperature      [T2I] = Downstream vertical profile'
    WRITE (INI,FMTA) '     Water type    [WTYPEC] = '//WTYPEC(JW)(4:8)//' water'
    WRITE (INI,FMTF) '     Ice thickness   [ICEI] = ',ICEI(JW),' m'
  END DO
  WRITE (INI,11250) 'Calculations', (JW,EVC(JW)(6:8),PRC(JW)(6:8),VBC(JW)(6:8),EBC(JW)(6:8),MBC(JW)(6:8),PQC(JW)(6:8),             &
                                     WINDC(JW)(6:8), QINC(JW)(6:8),QOUTC(JW)(6:8),HEATC(JW)(6:8),SLHTC(JW)(5:8),JW=1,NWB)
  DO JW=1,NWB
    WRITE (INI,11262) JW,SROC(JW)(6:8)
    WRITE (INI,11263) AFW(JW)
    WRITE (INI,11264) BFW(JW)
    WRITE (INI,11265) CFW(JW)
    WRITE (INI,11266) WINDH(JW)
    WRITE (INI,11267) RHEVC(JW)(6:8)
  END DO
  WRITE (INI,11270)  (JB,QINIC(JB)(6:8), DTRIC(JB)(6:8), HDIC(JB)(6:8), JB=1,NBR)
  WRITE (INI,11271)  (JW,METIC(JW)(6:8), JW =1,NWB)
  IF (TRIBUTARIES) WRITE (INI,11272)  (JT,TRIC(JT)(6:8),  JT =1,NTR)
  DO JB=1,NBR
    IF (NSTR(JB) /= 0) WRITE (INI,11273)  JB,(JS,STRIC(JS,JB)(6:8),JS=1,NSTR(JB))
  END DO
  WRITE (INI,11274) (JWD,WDIC(JWD)(6:8),JWD=1,NWD)
  WRITE (INI,FMTA) ' Meteorological Parameters'
  DO JW=1,NWB
    WRITE (INI,11300)  JW, LAT(JW), LONG(JW)
  END DO
  
  WRITE(INI,'(//"Meteorological Data Input Summary",/,"    Parameter   Waterbody   Average Value    Maximum     Minimum",/)')
  DO JW=1,NWB
      write(INI,'(4x,"TAIR(C)",8X,I3,7X,F10.3,3X,F10.3,2X,F10.3)')JW,tairav(jw),tairmax(jw),tairmin(jw)
      write(INI,'(4x,"TDEW(C)",8X,I3,7X,F10.3,3X,F10.3,2X,F10.3)')JW,tdewav(jw),tdewmax(jw),tdewmin(jw)
      write(INI,'(4x,"WIND(m/s)",6X,I3,7X,F10.3,3X,F10.3,2X,F10.3)')JW,windav(jw),windmax(jw),windmin(jw)
      write(INI,'(4x,"PHI(rad)",7X,I3,7X,F10.3,3X,F10.3,2X,F10.3)')JW,phiav(jw),phimax(jw),phimin(jw)
      write(INI,'(4x,"CLOUD(0-10)",4X,I3,7X,F10.3,3X,F10.3,2X,F10.3)')JW,cloudav(jw),cloudmax(jw),cloudmin(jw)
      write(INI,'(4x,"SRO(W/m2)",6X,I3,7X,F10.3,3X,F10.3,2X,F10.3,/)')JW,sroav(jw),sromax(jw),sromin(jw)
  ENDDO
  WRITE(INI,'(/)')
  
  IEGR = 1
  WRITE (INI,11330)
  DO WHILE (IEGR < IMX-1)
    IBGR = IEGR+1
    IEGR = IEGR+19
    IF (IEGR > IMX-1) IEGR = IMX-1
    WRITE (INI,11340) (I,I=IBGR,IEGR)
    WRITE (INI,11350) (PHI0(I),I=IBGR,IEGR)
  END DO
  WRITE (INI,11360) (JW,SLTRC(JW),THETA(JW),JW=1,NWB)
  WRITE (INI,11365)
  DO JW=1,NWB
    WRITE (INI,11400) JW
    WRITE (INI,11370) AX(JW),DX(JW),TSED(JW),CBHE(JW)   !,selectc
  END DO
  WRITE (INI,11375)
  DO JW=1,NWB
    WRITE (INI,11380) JW,ICEC(JW)(6:8),SLICEC(JW),ALBEDO(JW),HWI(JW),BETAI(JW),GAMMAI(JW)
  END DO
  WRITE (INI,11385)
  
  IF(WLC=='      ON')WRITE(INI,'(2X,"Water Level Output File: ON",T40,"Output file: wl.csv",T65,"Output frequency(days):",f10.4)')WLF
  
  IF(FLOWBALC == '      ON' .AND. VBC(1) == '      ON')WRITE(INI,'(2X,"Flow Balance Output File: ON",T40,"Output file: flowbal.csv",T65,"Output frequency(days):",f10.4)')FLOWBALF
  IF(FLOWBALC==  '      ON' .AND. VBC(1) /= '      ON')WRITE(INI,'(2X,"Flow Balance Output File is ON but VBC=OFF. No output file will be produced unless Volume Balance is ON.")') 
  IF(NPBALC=='      ON')WRITE(INI,'(2X,"N & P Mass Balance Output File: ON",T40,"Output file: massbal.csv",T65,"Output frequency(days):",f10.4)')NPBALF
  IF(LAKE_RIVER_CONTOUR_ON=='ON')THEN
      WRITE(INI,'(2x,"Lake River Contouor Output is ON and w2_lake_river_contour.csv has been read in."/)')
  ELSE
      WRITE(INI,'(2x,"Lake River Contouor Output is NOT ON but w2_lake_river_contour.csv has been read in."/)')
  ENDIF
  
  DO JW=1,NWB
    WRITE (INI,11390) JW,(HNAME(J),HPRWBC(J,JW)(6:8),J=1,NHY)
  END DO

  DO JW=1,NWB
    WRITE (INI,11400)  JW
    WRITE (INI,11410)  SNPC(JW)(6:8)
    IF (SNPC(JW)(6:8) == ' ON') WRITE (INI,11420)  NSNP(JW)
    IF (SNPC(JW)(6:8) == ' ON') WRITE (INI,11430) (SNPD(J,JW),J=1,NSNP(JW))
    IF (SNPC(JW)(6:8) == ' ON') WRITE (INI,11440) (SNPF(J,JW),J=1,NSNP(JW))
    WRITE (INI,11450)  SCRC(JW)(6:8)
    IF (SCRC(JW)(6:8) == ' ON') WRITE (INI,11455)  NSCR(JW)
    IF (SCRC(JW)(6:8) == ' ON') WRITE (INI,11460) (SCRD(J,JW),J=1,NSCR(JW))
    IF (SCRC(JW)(6:8) == ' ON') WRITE (INI,11470) (SCRF(J,JW),J=1,NSCR(JW))
    WRITE (INI,11520)  FLXC(JW)(6:8)
    IF (FLXC(JW)(6:8) == ' ON') WRITE (INI,11530)  NFLX(JW)
    IF (FLXC(JW)(6:8) == ' ON') WRITE (INI,11540) (FLXD(J,JW),J=1,NFLX(JW))
    IF (FLXC(JW)(6:8) == ' ON') WRITE (INI,11550) (FLXF(J,JW),J=1,NFLX(JW))
    WRITE (INI,11560)  VPLC(JW)(6:8)
    IF (VPLC(JW)(6:8) == ' ON') WRITE (INI,11570)  NVPL(JW)
    IF (VPLC(JW)(6:8) == ' ON') WRITE (INI,11580) (VPLD(J,JW),J=1,NVPL(JW))
    IF (VPLC(JW)(6:8) == ' ON') WRITE (INI,11590) (VPLF(J,JW),J=1,NVPL(JW))
    WRITE (INI,11600)  PRFC(JW)(6:8)
    IF (PRFC(JW)(6:8) == ' ON') WRITE (INI,11610)  NPRF(JW)
    IF (PRFC(JW)(6:8) == ' ON') WRITE (INI,11620)  NIPRF(JW)
    IF (PRFC(JW)(6:8) == ' ON') WRITE (INI,11630) (IPRF(I,JW),I=1,NIPRF(JW))
    IF (PRFC(JW)(6:8) == ' ON') WRITE (INI,11640) (PRFD(J,JW),J=1,NPRF(JW))
    IF (PRFC(JW)(6:8) == ' ON') WRITE (INI,11650) (PRFF(J,JW),J=1,NPRF(JW))
    WRITE (INI,11660)  SPRC(JW)(6:8)
    IF (SPRC(JW)(6:8) == ' ON') WRITE (INI,11670)  NSPR(JW)
    IF (SPRC(JW)(6:8) == ' ON') WRITE (INI,11680)  NISPR(JW)
    IF (SPRC(JW)(6:8) == ' ON') WRITE (INI,11690) (ISPR(I,JW),I=1,NISPR(JW))
    IF (SPRC(JW)(6:8) == ' ON') WRITE (INI,11700) (SPRD(J,JW),J=1,NSPR(JW))
    IF (SPRC(JW)(6:8) == ' ON') WRITE (INI,11710) (SPRF(J,JW),J=1,NSPR(JW))
    IF (SPRC(JW)(6:8) == 'ONV') WRITE (INI,11670)  NSPR(JW)
    IF (SPRC(JW)(6:8) == 'ONV') WRITE (INI,11680)  NISPR(JW)
    IF (SPRC(JW)(6:8) == 'ONV') WRITE (INI,11690) (ISPR(I,JW),I=1,NISPR(JW))
    IF (SPRC(JW)(6:8) == 'ONV') WRITE (INI,11700) (SPRD(J,JW),J=1,NSPR(JW))
    IF (SPRC(JW)(6:8) == 'ONV') WRITE (INI,11710) (SPRF(J,JW),J=1,NSPR(JW))

    WRITE (INI,11720)  CPLC(JW)(6:8)
    IF (CPLC(JW)(6:8) == ' ON') WRITE (INI,11730)  NCPL(JW)
    IF (CPLC(JW)(6:8) == ' ON') WRITE (INI,11740) (CPLD(J,JW),J=1,NCPL(JW))
    IF (CPLC(JW)(6:8) == ' ON') WRITE (INI,11750) (CPLF(J,JW),J=1,NCPL(JW))
  END DO
  WRITE (INI,11480)  TSRC(6:8)
  IF (TSRC(6:8) == ' ON') WRITE (INI,11490)  NTSR
  IF (TSRC(6:8) == ' ON') WRITE (INI,11500) (TSRD(J),J=1,NTSR)
  IF (TSRC(6:8) == ' ON') WRITE (INI,11510) (TSRF(J),J=1,NTSR)
  WRITE (INI,11760)  RSOC(6:8), RSIC(6:8)
  IF (RSOC(6:8) == ' ON') WRITE (INI,11770)  NRSO
  IF (RSOC(6:8) == ' ON') WRITE (INI,11780) (RSOD(J),J =1,NRSO)
  IF (RSOC(6:8) == ' ON') WRITE (INI,11790) (RSOF(J),J =1,NRSO)
  WRITE (INI,11800) (JB,NSTR(JB),JB=1,NBR)
  DO JB=1,NBR
    IF (NSTR(JB) > 0)  THEN
      WRITE (INI,11805)  JB
      WRITE (INI,11810) (JS,SINKC(JS,JB),WSTR(JS,JB),ESTR(JS,JB),KBSTR(JS,JB),JS=1,NSTR(JB))
    END IF
  END DO
  WRITE (INI,11840)  NWD, (IWD(JW),JW=1,NWD)
  WRITE (INI,11860)  NTR, (ITR(JT),JT=1,NTR)
  IF (TRIBUTARIES) THEN
    WRITE (INI,11870) (PTRC(JT),     JT=1,NTR)
    WRITE (INI,11880) (ETTR(JT),     JT=1,NTR)
    WRITE (INI,11890) (EBTR(JT),     JT=1,NTR)
  END IF
  IF (SPILLWAYS) THEN
    WRITE (INI,11821)  NSP,(ESP(JS),JS=1,NSP)
    WRITE (INI,11822) (IUSP(JS),    JS=1,NSP)
    WRITE (INI,11823) (IDSP(JS),    JS=1,NSP)
    WRITE (INI,11824) (PUSPC(JS),   JS=1,NSP)
    WRITE (INI,11825) (PDSPC(JS),   JS=1,NSP)
    WRITE (INI,11826) (ETUSP(JS),   JS=1,NSP)
    WRITE (INI,11827) (ETDSP(JS),   JS=1,NSP)
    WRITE (INI,11828) (EBUSP(JS),   JS=1,NSP)
    WRITE (INI,11829) (EBDSP(JS),   JS=1,NSP)
    WRITE (INI,11830) (A1SP(JS),B1SP(JS),A2SP(JS),B2SP(JS),JS=1,NSP)
    WRITE (INI,11865) (GASSPC(JS),EQSP(JS),AGASSP(JS),BGASSP(JS),CGASSP(JS),JS=1,NSP)
    DO JS=1,NSP
    IF(GASSPC(JS) =='      ON' .AND. EQSP(JS)==4 .AND. BGASSP(JS)==1.0)THEN
         WRITE (SEGNUM,'(I0)') JS
         SEGNUM = ADJUSTL(SEGNUM)  
         L = LEN_TRIM(SEGNUM) 
         TITLE240='w2_sp'//SEGNUM(1:L)//'DO.csv'
        WRITE(INI,'(5X,A,I2,A,A)')'For spillway [',JS,'], since EQSP=4, an external file is read in for spillway gas fraction:',trim(title240)
    ENDIF
    ENDDO

  END IF
  IF (GATES) THEN
    WRITE (INI,11891)  NGT,(EGT(JG),JG=1,NGT)
    WRITE (INI,11892) (IUGT(JG),    JG=1,NGT)
    WRITE (INI,11893) (IDGT(JG),    JG=1,NGT)
    WRITE (INI,11894) (PUGTC(JG),   JG=1,NGT)
    WRITE (INI,11895) (PDGTC(JG),   JG=1,NGT)
    WRITE (INI,11896) (ETUGT(JG),   JG=1,NGT)
    WRITE (INI,11897) (ETDGT(JG),   JG=1,NGT)
    WRITE (INI,11898) (EBUGT(JG),   JG=1,NGT)
    WRITE (INI,11899) (EBDGT(JG),   JG=1,NGT)
    WRITE (INI,11820) (A1GT(JG),B1GT(JG),G1GT(JG),A2GT(JG),B2GT(JG),G2GT(JG),JG=1,NGT)
    WRITE (INI,11866) (GASGTC(JT),EQGT(JT),AGASGT(JT),BGASGT(JT),CGASGT(JT),JT=1,NGT)
    DO JT=1,NGT
    IF(GASGTC(JT)=='      ON' .AND. EQGT(JT)==4 .AND. BGASGT(JT)==1.0)THEN
         WRITE (SEGNUM,'(I0)') JT
         SEGNUM = ADJUSTL(SEGNUM)  
         L = LEN_TRIM(SEGNUM) 
         TITLE240='w2_gt'//SEGNUM(1:L)//'DO.csv'
        WRITE(INI,'(5X,A,I2,A,A)')'For gate [',JT,'], since EQGT=4, an external file is read in for gate gas fraction:',trim(title240)
    ENDIF
    ENDDO
    
  END IF
  IF (PIPES) THEN
    WRITE (INI,11831)  NPI,(EUPI(JP),JP=1,NPI)
    WRITE (INI,11841) (EDPI(JP),     JP=1,NPI)
    WRITE (INI,11842) (WPI(JP),      JP=1,NPI)
    WRITE (INI,11843) (DLXPI(JP),    JP=1,NPI)
    WRITE (INI,11844) (FPI(JP),      JP=1,NPI)
    WRITE (INI,11845) (FMINPI(JP),   JP=1,NPI)
    WRITE (INI,11832) (IUPI(JP),     JP=1,NPI)
    WRITE (INI,11833) (IDPI(JP),     JP=1,NPI)
    WRITE (INI,11834) (PUPIC(JP),    JP=1,NPI)
    WRITE (INI,11835) (PDPIC(JP),    JP=1,NPI)
    WRITE (INI,11836) (ETUPI(JP),    JP=1,NPI)
    WRITE (INI,11837) (ETDPI(JP),    JP=1,NPI)
    WRITE (INI,11838) (EBUPI(JP),    JP=1,NPI)
    WRITE (INI,11839) (EBDPI(JP),    JP=1,NPI)
  END IF
  IF (PUMPS) THEN
    WRITE (INI,11851)  NPU,(EPU(JP),JP=1,NPU)
    WRITE (INI,11852) (IUPU(JP),    JP=1,NPU)
    WRITE (INI,11853) (IDPU(JP),    JP=1,NPU)
    WRITE (INI,11854) (STRTPU(JP),  JP=1,NPU)
    WRITE (INI,11855) (ENDPU(JP),   JP=1,NPU)
    WRITE (INI,11856) (EONPU(JP),   JP=1,NPU)
    WRITE (INI,11857) (EOFFPU(JP),  JP=1,NPU)
    WRITE (INI,11858) (QPU(JP),     JP=1,NPU)
    WRITE (INI,11859) (PPUC(JP),    JP=1,NPU)
    WRITE (INI,11861) (ETPU(JP),    JP=1,NPU)
    WRITE (INI,11862) (EBPU(JP),    JP=1,NPU)
    WRITE (INI,11863) (KTPU(JP),    JP=1,NPU)
    WRITE (INI,11864) (KBPU(JP),    JP=1,NPU)
  END IF
  WRITE (INI,11900)
  WRITE (INI,11910) (JB,DTRC(JB)(6:8),JB=1,NBR)
  WRITE (INI,11930) CONFN, RSIFN, QWDFN
  IF(SELECTC=='      ON'.or. SELECTC == '    USGS')WRITE(INI,'(3X,"w2_selective.npt read in and used during simulation.",T65,"SELECTC=",A)')ADJUSTR(SELECTC)
  IF(HABTATC=='      ON')WRITE(INI,'(3X,"w2_habitat.npt read in and used during simulation.",T65,"# of fish classes:",i5)')IFISH
  IF(ENVIRPC=='      ON')WRITE(INI,'(3X,"w2_envirprf.npt read in and used during simulation.",T65,"# of bins:",i5)')NUMCLASS
  IF(AERATEC=='      ON')WRITE(INI,'(3X,"w2_aerate.npt read in and used during simulation.",T65,"# of aerators:",i5)')NAER
  IF(INITUWL=='      ON')WRITE(INI,'(A)')'   INITUWL=ON. W2 will calculate initial water surface profile and initial velocity regime.'
  IF(CONSTRICTION)WRITE(INI,'(3X,"w2_constriction.csv read in and used during simulation.",T65,"# of constrictions:",i5)')ICON
  IF(FISH_PARTICLE_EXIST.AND.PARTON == 'ON')WRITE(INI,'(3X,"w2_particle.csv read in and used during simulation.",T65,"# of particle releases:",i5)')NFISHSEG
  IF(FISH_PARTICLE_EXIST.AND.PARTON /= 'ON')WRITE(INI,'(3X,"w2_particle.csv read in but NOT used during simulation.")')
  IF(MULTIPLEWB.and.MWB/='ON')WRITE(INI,'(3X,"w2_multiple_WB.npt read in but NOT used during simulation.")')
  IF(MULTIPLEWB.and.MWB=='ON')THEN
    WRITE(INI,'(3X,"w2_multiple_WB.npt read in and used during simulation.",T65,"# of inputs:",i3)')IMWB
      DO J=1,IMWB
          WRITE(INI,'(T12,"Multiple WB input #:",i3,T40,"Type of input:",a2,t65," Directory:",a)')J,MWB2(J),FILEDIRMWB(J)
      ENDDO
  ENDIF
  
  IF(SYSTDG .and. SYSTDGC=='      ON')WRITE(INI,'(A)')'   w2_systdg.npt read in for gate TDG impacts and the model is used for TDG calculations.'
  IF(SYSTDG .and. SYSTDGC=='     OFF')WRITE(INI,'(A)')'   w2_systdg.npt read in for gate TDG impacts BUT NOT used for TDG calculations since SYSTDGC=OFF.'
  IF(TECPLOTBR)WRITE(INI,'(A,T65,A,<JTEC>(I3,","))')'   w2_tecplotbr.csv file is read in and TECPLOT output is ON.','ACTIVE BRANCHES:',(IBR(I),I=1,JTEC)
  IF(SED_DIAG=='      ON'.AND.DIAG_EXISTS)THEN
        IF(sd_global)WRITE(INI,'(A)')'   w2_diagenesis.npt file is read in and Global Diagenesis is ON:'
        IF(.Not.sd_global)WRITE(INI,'(A)')'   w2_diagenesis.npt file is read in and Global Diagenesis is OFF:'
        IF(IncludeFFTLayer)then
            WRITE(INI,'(T20,A)')'FFT layer is ON'
        else
            WRITE(INI,'(T20,A)')'FFT layer is OFF'
        endif
        if(IncludeBedConsolidation)then
            WRITE(INI,'(T20,A)')'Bed Consolidation is ON'
        else
            WRITE(INI,'(T20,A)')'Bed Consolidation is OFF'
        endif
        if(IncludeCEMASedDiagenesis)then
        WRITE(INI,'(T20,A)')'Sediment Diagenesis Processes are ON'
        else
        WRITE(INI,'(T20,A)')'Sediment Diagenesis Processes are OFF'
        endif
        if(Bubbles_Calculation)WRITE(INI,'(T20,A)')'Sediment Diagenesis Bubbles are ON'
        if(CEMA_POM_Resuspension)WRITE(INI,'(T20,A)')'Sediment POM resuspension processes are ON'
        if(IncludeAlkalinity)WRITE(INI,'(T20,A)')'Sediment pH and Alkalinity processes are ON'
        if(IncludeIron)WRITE(INI,'(T20,A)')'Sediment Fe processes are ON'
        if(IncludeManganese)WRITE(INI,'(T20,A)')'Sediment Mn processes are ON'
        if(.not.Bubbles_Calculation)WRITE(INI,'(T20,A)')'Sediment Diagenesis Bubbles are OFF'
        if(.not.CEMA_POM_Resuspension)WRITE(INI,'(T20,A)')'Sediment POM resuspension processes are OFF'
        if(.not.IncludeAlkalinity)WRITE(INI,'(T20,A)')'Sediment pH and Alkalinity processes are OFF'
        if(.not.IncludeIron)WRITE(INI,'(T20,A)')'Sediment Fe processes are OFF'
        if(.not.IncludeManganese)WRITE(INI,'(T20,A)')'Sediment Mn processes are OFF'
        if(DYNAMIC_SD)WRITE(INI,'(T20,A)')'Sediment dynamic accumulation of labile to refractory organics is ON since SDRegnPOC_L_Fr<0. Hence, the ratio is not fixed over time.'
  endif
  

  WRITE(INI,*)
  DO JW=1,NWB
    WRITE (INI,11940) JW, BTHFN(JW), METFN(JW), VPRFN(JW), LPRFN(JW)
    IF(ATM_DEPOSITIONC(JW)=='      ON')THEN
    WRITE (INI,'(5X,"Atmospheric Deposition           = ",A72)') ATMDEPFN(JW)
    ENDIF
  END DO
  DO JB=1,NBR
    WRITE (INI,11950) JB, QINFN(JB), TINFN(JB), CINFN(JB), QOTFN(JB), QDTFN(JB), TDTFN(JB), CDTFN(JB), PREFN(JB), TPRFN(JB),       &
                          CPRFN(JB), EUHFN(JB), TUHFN(JB), CUHFN(JB), EDHFN(JB), TDHFN(JB), CDHFN(JB)
  END DO
  WRITE (INI,11970) (JT,QTRFN(JT),TTRFN(JT),CTRFN(JT),JT=1,NTR)
  WRITE (INI,11980)  ERRFN, WRNFN, TSRFN, WDOFN
  DO JW=1,NWB
    WRITE (INI,11990) JW, SNPFN(JW), FLXFN(JW), PRFFN(JW), VPLFN(JW), CPLFN(JW)  !put in shade, wsc, graph, etc. "to do"
  END DO
  IF (CONSTITUENTS) THEN
    WRITE (INI,12000) CCC(6:8), LIMC(6:8), CUF
    IF(CO2ATM_REGRESSION=='      ON')THEN
        IF(YEAR<1980)THEN
          CO2ATM  = 0.000041392*YEAR**3 - 0.231409975*YEAR**2 + 430.804190829*YEAR - 266735.857433224  
        ELSE
          CO2ATM  = 0.015903*YEAR**2 - 61.799598*YEAR + 60357.055057
        ENDIF
        IF(YEAR<1920)THEN
            CALL ERRORS
            WRITE (ERR,'(A)') 'CO2 atmospheric concentrations outside range of regression. YEAR must be > 1920. Suggest entering a specific value of CO2ATM in ppm rather than use the regression.'
        ENDIF
        
        WRITE (INI,'(A,I4,A,F10.2,A)')'   CO2 atmospheric concentration is computed internally. For the start year of simulation [YEAR=',YEAR,'], CO2 atmosphere:',CO2ATM,' ppm CO2'
    ELSE
        WRITE (INI,'(A,F10.2,A)') '   CO2 atmospheric concentration was specified by the user as',CO2ATMPPM,' ppm CO2'
        IF(CO2ATMPPM > 450. .OR. CO2ATMPPM < 300.)THEN
            CALL WARNINGS
            WRITE(WRN,'(A,F10.2)') 'CO2 atmospheric value entered by the user is either > 450 ppm or < 300 ppm. The CO2 atm concentration in ppm is ',CO2ATMPPM
        ENDIF
        
    ENDIF
        
    DO JW=1,NWB
      WRITE (INI,12030) JW
      IF(CDWBC(NH3_DER,JW)=='      ON')THEN
        WRITE(INI,'(A,I3,T65,A,F10.2)')'     NH3 derived constituent active and gas loss modeled.'       !JW=',JW,'KG_H2O_CONSTANT=',KG_H2O_CONSTANT(JW)
        IF(KG_H2O_CONSTANT(JW) /= 168.)THEN
            CALL WARNINGS
            WRITE (WRN,'(A,f10.2,A,I4)')'NH3 gas loss: KG_H2O_CONSTANT is not equal to recommended value of 168. The value used is',kg_h2o_constant(jw),' for WB:',jw  
        ENDIF
        IF(CDWBC(PH_DER,JW)/='      ON')THEN
            CALL ERRORS
            WRITE(ERR,'(A,I4)')'NH3 derived variable is ON but pH is not ON. pH must be on to compute NH3 for waterbody #:',JW
        ENDIF
        
      ENDIF
      WRITE (INI,12020)
      WRITE (INI,12040) (CNAME1(JC),CAC(JC)(6:8),C2IWB(JC,JW),CPRWBC(JC,JW)(6:8),JC=1,NCT),'Sediments          ',& ! (CNAME1(JC),CAC(JC)(6:8),C2IWB(JC,JW),CFWBC(JC,JW)(6:8),CPRWBC(JC,JW)(6:8),JC=1,NCT),'Sediments          ',&
                         SEDC(JW)(6:8), SEDCI(JW), PRNSC(JW)(6:8)                                            !SEDC(JW)(6:8), SEDCI(JW),'OFF', PRNSC(JW)(6:8)
     ! WRITE OUT IF FLUXES ARE on/off  WRITE (INI,12041) (CFWBC(JC,JW)(6:8),JC=1,NFL) - add name of flux variable to write out - change kfname2 to an array - low priority
      DO JB=BS(JW),BE(JW)
        WRITE (INI,12045) JB,(CNAME1(JC),CINBRC(JC,JB)(6:8),CDTBRC(JC,JB)(6:8),CPRBRC(JC,JB)(6:8),JC=1,NCT)
      END DO
      WRITE (INI,12050)
      WRITE (INI,12060)
      DO JD=1,NDC
        WRITE (INI,12070) CDNAME(JD),CDWBC(JD,JW)(6:8)
      END DO
    END DO
    DO JT=1,NTR
      WRITE (INI,12046) JT,(CNAME1(JC),CINTRC(JC,JT)(6:8),JC=1,NCT)
    END DO
    DO JW=1,NWB
      WRITE (INI,12030)  JW
      WRITE (INI,12080) 'Constituent Rates'
      WRITE (INI,12100) (CNAME1(JG+1), CGQ10(JG), CG0DK(JG), CG1DK(JG), CGS(JG),JG=1,NGC)
      WRITE (INI,12090) (SSS(NS),NS=1,NSS)
      WRITE (INI,12110)  LDOMDK(JW), LRDDK(JW)
      WRITE (INI,12120)  RDOMDK(JW)
      WRITE (INI,12130)  LPOMDK(JW), LRPDK(JW), POMS(JW)
      WRITE (INI,12150)  RPOMDK(JW)
      DO JA=1,NAL
        WRITE (INI,12160) JA, AG(JA), AM(JA), AE(JA), AR(JA), AS(JA), ALGP(JA), ALGN(JA), ALGC(JA), ALGSI(JA), ACHLA(JA), APOM(JA), AVERTM(JA)
      END DO
        DO JA=1,NEPT
        WRITE (INI,12162) JA, EG(JA), EM(JA), EE(JA), ER(JA), EB(JA), EP(JA), EN(JA), EC(JA), ESI(JA), ECHLA(JA), EPOM(JA)
      END DO

      DO JZ=1,NZP
        WRITE (INI,12161) JZ, ZG(JZ), ZM(JZ), ZR(JZ), ZS(JZ), ZP(JZ), ZN(JZ), ZC(JZ)
      END DO

      WRITE (INI,12180)  PO4R(JW)
      WRITE (INI,12190)  NH4DK(JW), NH4R(JW),KG_H2O_CONSTANT(JW)
      WRITE (INI,12200)  NO3DK(JW),NO3S(JW),FNO3SED(JW)
      WRITE (INI,12210)  DSIR(JW), PSIDK(JW), PSIS(JW)
      WRITE (INI,12220)  SEDDK(JW)
!      WRITE (INI,12230)  FES(JW), FER(JW)
      WRITE (INI,12240) (SOD(I),I=US(BS(JW)),DS(BE(JW)))
      WRITE (INI,12250)  FSOD(JW), FSED(JW)
      WRITE (INI,12255)  SEDS(jw), SEDBR(jw)
    END DO
    WRITE (INI,12330) 'Upper Temperature Bounds'
    DO JW=1,NWB
      WRITE (INI,'(A,I2)')'   Waterbody #',JW
      WRITE (INI,12290)  NH4T1(JW), NH4T2(JW)
      WRITE (INI,12300)  NO3T1(JW), NO3T2(JW)
      WRITE (INI,12310)  OMT1(JW),  OMT2(JW)
      WRITE (INI,12320)  SODT1(JW), SODT2(JW)
    END DO
    DO JA=1,NAL
      WRITE (INI,12280) JA, AT1(JA), AT2(JA)
    END DO
    WRITE (INI,12270) 'Lower Temperature Bounds'
    DO JA=1,NAL
      WRITE (INI,12340) JA, AT3(JA), AT4(JA)
    END DO
      WRITE (INI,12355) 'Algae Oxygen',(JA,O2AR(JA),O2AG(JA),JA=1,NAL)
      WRITE (INI,12356) 'Periphyton Oxygen',(JA,O2ER(JA),O2EG(JA),JA=1,NEPT) 
      WRITE (INI,12380) 'Half Saturation', (JA,AHSP(JA), AHSN(JA),  AHSSI(JA),JA=1,NAL)
      WRITE (INI,12400) 'Attentuation', (JA,EXA(JA), JA=1,NAL)
      WRITE (INI,12410) 'Light Saturation',(JA,ASAT(JA),JA=1,NAL)
      WRITE (INI,12420) 'Diffusion',                DMO2,      DMCO2
      
      DO JW=1,NWB
      WRITE (INI,'(A,I2)')'   Waterbody #',JW
      WRITE (INI,12350) 'Stoichiometric Equivalence', O2NH4(JW), O2OM(JW)
      WRITE (INI,12360)     ORGC(JW),ORGP(JW), ORGN(JW), ORGSI(JW)
      WRITE (INI,12390) 'Light',               BETA(JW), EXH2O(JW), EXSS(JW), EXOM(JW)
      WRITE (INI,12430) 'Partitioning Coefficients',PARTP(JW), PARTSI(JW)
      WRITE (INI,12440) 'Miscellaneous Constants',  O2LIM,     CO2R(JW)
      WRITE (INI,12391) 'Bacteria',bact1dk(jw),bactq10(jw),bactldk(jw),bacts(jw)
      WRITE (INI,12392) 'H2S',H2SR(JW),H2SQ10(JW),H2S1DK(JW),SO4R(JW)
      WRITE (INI,12393) 'CH4',CH4R(JW),CH4Q10(JW),CH41DK(JW)
      WRITE (INI,12394) 'FE',FEIIR(JW),KFE_OXID(JW),KFE_RED(JW),KFEOOH_HALFSAT(JW),FESETVEL(JW)
      WRITE (INI,12395) 'MN',MNIIR(JW),KMN_OXID(JW),KMN_RED(JW),KMNO2_HALFSAT(JW),MNSETVEL(JW)
      WRITE (INI,12396) 'Reaeration',REAERC(JW),NEQN(JW),DGPO2(JW),MINKL(JW)
    END DO
    IF (NBOD > 0) WRITE (INI,12370)  (JB, KBOD(JB), TBOD(JB), RBOD(JB), JB=1,NBOD)                                     !TC 06/18/02


!** Inflow constituent statistics                                                                                      !SW 01/07/01

    WRITE (INI,'(/A)') '  Inflow Constituent Statistics'
    DO JW=1,NWB
      DO JB=BS(JW),BE(JW)
        IF (UQ_EXTERNAL(JB)) THEN
          WRITE (INI,FMTI) '    Branch ',JB
          WRITE (INI,FMTA)    '       Constituent name        Average     Maximum     Minimum ApproxLoading(kg/d)'
          DO JC=1,NACIN(JB)
            WRITE (INI,'(7X,A,3(2X,F10.3),2x,e12.5)') CNAME1(INCN(JC,JB))(1:19), CINAVG(INCN(JC,JB),JB), CINMAX(INCN(JC,JB),JB),            &
                                             CINMIN(INCN(JC,JB),JB),cinload(INCN(JC,JB),JB)
          END DO
        END IF
      END DO
    END DO

!** Distributary Tributary constituent statistics           - need to include average and max/min

    DO JB=1,NBR
      if(dist_trib(jb))then
      WRITE (INI,FMTI) '    Distributed Tributary ',JB
      WRITE (INI,FMTA)    '       Constituent name        Average     Maximum     Minimum ApproxLoading(kg/d)'
      DO JC=1,NACDT(JB)
        WRITE (INI,'(7X,A,3(2X,F10.3),2x,e12.5)') CNAME1(INCDT(JC,JB))(1:19),CDTAVG(INCDT(JC,JB),JB), CDTMAX(INCDT(JC,JB),JB),            &
                                             CDTMIN(INCDT(JC,JB),JB),cdtload(incdt(jc,jb),jb)            !CTRAVG(TRCN(JC,JT),JT), CTRMAX(TRCN(JC,JT),JT),CTRMIN(TRCN(JC,JT),JT)
      END DO
      endif
    END DO
    
    !** External Head constituent statistics          

    DO JB=1,NBR
      if(uh_external(jb))then
      WRITE (INI,FMTI) '    External Head Upstream ',JB
      WRITE (INI,FMTA)    '       Constituent name             Average            Maximum        Minimum '
      DO JC=1,NAC
        WRITE (INI,'(7X,A,3(2X,F16.3),2x,e12.5)') CNAME1(CN(JC))(1:19),CUHAVG(JC,JB), CUHMAX(JC,JB),  CUHMIN(JC,JB)            
      END DO
      endif
    END DO
    
    !** External Head constituent statistics          

    DO JB=1,NBR
      if(DH_external(jb))then
      WRITE (INI,FMTI) '    External Head Downstream ',JB
      WRITE (INI,FMTA)    '       Constituent name             Average            Maximum        Minimum '
      DO JC=1,NAC
        WRITE (INI,'(7X,A,3(2X,F16.3),2x,e12.5)') CNAME1(CN(JC))(1:19),CDHAVG(JC,JB), CDHMAX(JC,JB),  CDHMIN(JC,JB)            
      END DO
      endif
    END DO

!** Tributary constituent statistics                                                                                   !SW 01/07/01

    DO JT=1,NTR
      WRITE (INI,FMTI) '     Tributary ',JT
      WRITE (INI,FMTA)    '       Constituent name        Average     Maximum     Minimum ApproxLoading(kg/d)'
      DO JC=1,NACTR(JT)
        WRITE (INI,'(7X,A,3(2X,F10.3), 2x,e12.5)') CNAME1(TRCN(JC,JT))(1:19), CTRAVG(TRCN(JC,JT),JT), CTRMAX(TRCN(JC,JT),JT),                &
                                         CTRMIN(TRCN(JC,JT),JT),ctrload(trCN(JC,Jt),Jt)
      END DO
    END DO


!** Overall Approx Loading for all inflows, tribs, dist tribs
      WRITE (INI,FMTA) '    Overall approx loading for all inflows, tributaries, and distributed tributaries'
      WRITE (INI,FMTA)    '       Constituent name      ApproxLoading(kg/d)'
DO J=1,NCT
  xx=0.0
    do jb=1,nbr
      do jc=1,nacin(jb)
        IF(incn(jc,jb) == j)then
        xx=xx+cinload(incn(jc,jb),jb)
        endif
      enddo
      do jc=1,nacdt(jb)
        IF(incdt(jc,jb) == j)then
        xx=xx+cdtload(incdt(jc,jb),jb)
        endif
      enddo
    enddo
    do jt=1,ntr
        do jc=1,nactr(jt)
        IF(trcn(jc,jt) == j)then
        xx=xx+ctrload(trcn(jc,jt),jt)
        endif
        enddo
    enddo
   if(xx>0.0)WRITE (INI,'(7X,A,2x,e12.5)') CNAME1(j)(1:19),xx
enddo

WRITE(INI,*)
WRITE(INI,*)'Atmospheric Loading in kg/km2/year'
DO JW=1,NWB
IF(ATM_DEPOSITIONC(JW)=='      ON')THEN
WRITE(INI,'(4x,A)')'WQ STATE VARIABLE        WATERBODY        AVERAGE         MINIMUM           MAXIMUM'
DO J=1,NCT
    WRITE(INI,'(7X,A,3X,I5,10X,E12.5,5X,E12.5,5X,E12.5)')CNAME1(j)(1:19),JW,ATM_DEP_AVG(J,JW),ATM_DEP_MIN(J,JW),ATM_DEP_MAX(J,JW)
ENDDO
ELSE
    WRITE(INI,'(A,I5)')'NO ATMOSPHERIC DEPOSITION FOR WATERBODY:',JW
ENDIF

ENDDO

DO JW=1,NWB
   IF(ATM_DEPOSITIONC(JW)=='      ON')THEN
     WRITE(INI,'(A)')'Approximate Atmospheric Deposition Loading in kg/day'
     EXIT
   ENDIF
ENDDO

DO J=1,NCT
   DO JW=1,NWB
   IF(ATM_DEPOSITIONC(JW)=='      ON')THEN
       WRITE(INI,'(7X,A,2x,e12.5)') CNAME1(j)(1:19),ATM_DEP_AVG(J,JW)*SAGR(KTWB(JW),JW)/365./1000./1000.
   ENDIF
   ENDDO
ENDDO
  ENDIF
  
  WRITE (INI,'(/1X,A)') 'Water Balance Summary'
  DO JW=1,NWB
    WRITE (INI,12460) JW, QINAVW(JW), QINMXW(JW), QOTAVW(JW), QOTMXW(JW)
    DO JB=BS(JW),BE(JW)
      WRITE (INI,12470) JB, QINAVB(JB),QINMXB(JB),QINAV(JB), QINMX(JB), QTRAVB(JB),QTRMXB(JB),QDTAV(JB),QDTMX(JB),QPRAV(JB),       &
                            QPRMX(JB), QSTRAV(JB),QSTRMX(JB),QWDAVB(JB),QWDMXB(JB)
    END DO
  END DO

  !** Inflow TEMPERATURE STATS                                                                                      !SW 01/07/01

    WRITE (INI,'(/A)') '  Branch Inflow Temperature Min/Max'
    WRITE (INI,FMTA)    '  Branch(JB)    Maximum Temp(C)     Minimum Temp(C) '
    DO JW=1,NWB
      DO JB=BS(JW),BE(JW)
        IF (UQ_EXTERNAL(JB)) THEN

          WRITE (INI,'(7X,I4,6X,F10.3,9x,F10.3)') JB,TINMAX(JB),TINMIN(JB)
        END IF
      END DO
    END DO

     !** Tributary TEMPERATURE STATS                                                                                      !SW 01/07/01

   IF(NTR>0) WRITE (INI,'(/A)') '  Tributary Inflow Temperature Min/Max'
    DO JT=1,NTR
          WRITE (INI,FMTA)    '  Tributary(JT)    Maximum Temp(C)     Minimum Temp(C) '
          WRITE (INI,'(7X,I4,6X,F10.3,9x,F10.3)') JT,TRTMAX(JT),TRTMIN(JT)
    END DO
    
       !** EXTERNAL HEAD TEMPERATURE STATS                                                                                      !SW 01/07/01

   DO JB=1,NBR
       IF(UH_EXTERNAL(JB))THEN
           WRITE (INI,'(/A,i4)') '  Upstream Head Boundary Temperature Min/Max for Branch:',jb

          WRITE (INI,FMTA)    'Maximum Temp(C)     Minimum Temp(C) '
          WRITE (INI,'(F10.3,9x,F10.3)') tuhmax(jb),tuhmin(jb)
       ENDIF
       IF(DH_EXTERNAL(JB))THEN
           WRITE (INI,'(/A,i4)') '  Downstream Head Boundary Temperature Min/Max for Branch:',jb

          WRITE (INI,FMTA)    'Maximum Temp(C)     Minimum Temp(C) '
          WRITE (INI,'(F10.3,9x,F10.3)') tdhmax(jb),tdhmin(jb)
       ENDIF
    ENDDO



  WRITE (INI,12480) 'Geometry', IMX, KMX, NBR
  DO JW=1,NWB
    KBGR = 0
    KEGR = 0
    WRITE (INI,12490) JW, US(BS(JW))-1, DS(BE(JW))+1, BS(JW), BE(JW), ELBOT(JW), KTWB(JW)
    DO WHILE (KEGR < KMX)
      KBGR = KEGR+1
      KEGR = MIN(KEGR+19,KMX)
      WRITE (INI,12500) (K,K=KBGR,KEGR)
      WRITE (INI,12510) (H(K,JW),K=KBGR,KEGR)
    END DO
    DO JB=BS(JW),BE(JW)
      WRITE (INI,12520) JB, US(JB), DS(JB), UHS(JB), DHS(JB), NL(JB)
    END DO
  END DO

! Volume-area-elevation table

  DO JW=1,NWB
    IF (NBR > 1) THEN
      DO JB=BS(JW),BE(JW)
        WRITE (INI,12600) VOLB(JB)
        WRITE (INI,12620) JB
        WRITE (INI,12640)
        DO K=2,KMX-1
          IF (K == KTWB(JW)) THEN
            WRITE (INI,12650) K, EL(K,DS(JBDN(JW))), SABR(K,JB)/1.E6, CVBR(K,JB)/1.E6, NCCBR(K,JB), HBR(K,JB), BBR(K,JB)
          ELSE
            WRITE (INI,12655) K, EL(K,DS(JBDN(JW))), SABR(K,JB)/1.E6, CVBR(K,JB)/1.E6, NCCBR(K,JB), HBR(K,JB), BBR(K,JB)
          END IF
        END DO
      END DO
    END IF
    WRITE (INI,12670) JW, VOLG(JW)
    WRITE (INI,12690) JW
    WRITE (INI,12630)
    DO K=2,KMX-1
      IF (QOTAVW(JW) > 0.0) THEN
        IF (K == KTWB(JW)) THEN
          WRITE (INI,12650) K, EL(K,DS(JBDN(JW))), SAGR(K,JW)/1.E6, CVGR(K,JW)/1.E6, NCCGR(K,JW), HGR(K,JW), BGR(K,JW), CVGR(K,JW) &
                            /(QOTAVW(JW)*86400.0)
        ELSE
          WRITE (INI,12655) K, EL(K,DS(JBDN(JW))), SAGR(K,JW)/1.E6, CVGR(K,JW)/1.E6, NCCGR(K,JW), HGR(K,JW), BGR(K,JW), CVGR(K,JW) &
                            /(QOTAVW(JW)*86400.0)
        END IF
      ELSE IF (QINAVW(JW) > 0.0) THEN
        IF (K == KTWB(JW)) THEN
          WRITE (INI,12650) K, EL(K,DS(JBDN(JW))), SAGR(K,JW)/1.E6, CVGR(K,JW)/1.E6, NCCGR(K,JW), HGR(K,JW), BGR(K,JW), CVGR(K,JW) &
                            /(QINAVW(JW)*86400.0)
        ELSE
          WRITE (INI,12655) K, EL(K,DS(JBDN(JW))), SAGR(K,JW)/1.E6, CVGR(K,JW)/1.E6, NCCGR(K,JW), HGR(K,JW), BGR(K,JW), CVGR(K,JW) &
                            /(QINAVW(JW)*86400.0)
        END IF
      ELSE
        IF (K == KTWB(JW)) THEN
          WRITE (INI,12650) K, EL(K,DS(JBDN(JW))), SAGR(K,JW)/1.E6, CVGR(K,JW)/1.E6, NCCGR(K,JW), HGR(K,JW), BGR(K,JW)
        ELSE
          WRITE (INI,12655) K, EL(K,DS(JBDN(JW))), SAGR(K,JW)/1.E6, CVGR(K,JW)/1.E6, NCCGR(K,JW), HGR(K,JW), BGR(K,JW)
        END IF
      END IF
    END DO
  END DO



    VOLG1=0.0
    QOTAVW1=0.0
    QINAVW1=0.0
    SAGR1=0.0
    CVGR1=0.0
    NCCGR1=0
    HGR1=0.0
    BGR1=0.0
  DO JW=1,NWB
    VOLG1=VOLG1+VOLG(JW)
    QOTAVW1=QOTAVW1+QOTAVW(JW)
    QINAVW1=QINAVW1+QINAVW(JW)
    DO K=2,KMX-1
        SAGR1(K)=SAGR1(K)+SAGR(K,JW)
        CVGR1(K)=CVGR1(K)+CVGR(K,JW)
        NCCGR1(K)=NCCGR1(K)+NCCGR(K,JW)
        HGR1(K)=HGR1(K)+HGR(K,JW)
        BGR1(K)=BGR1(K)+BGR(K,JW)
    ENDDO
  ENDDO
    WRITE (INI,*)
    WRITE (INI,'(A140)')'The following grid summary for all waterbodies is only valid for all waterbodies having the same EBOT. Otherwise ignore this summary.'
    WRITE (INI,12671) VOLG1
    WRITE (INI,12691)
    WRITE (INI,12630)
    DO K=2,KMX-1
      IF (QOTAVW1 > 0.0) THEN

          WRITE (INI,12655) K, EL(K,DS(JBDN(1))), SAGR1(K)/1.E6, CVGR1(K)/1.E6, NCCGR1(K), HGR1(K), BGR1(K), CVGR1(K) &
                            /(QOTAVW1*86400.0)

      ELSE IF (QINAVW1 > 0.0) THEN                     ! SW 6/11/14

          WRITE (INI,12655) K, EL(K,DS(JBDN(1))), SAGR1(K)/1.E6, CVGR1(K)/1.E6, NCCGR1(K), HGR1(K), BGR1(K), CVGR1(K) &
                            /(QINAVW1*86400.0)

      ELSE

          WRITE (INI,12655) K, EL(K,DS(JBDN(1))), SAGR1(K)/1.E6, CVGR1(K)/1.E6, NCCGR1(K), HGR1(K), BGR1(K)

      END IF
    END DO
    WRITE(INI,*)


!** Bathymetry
  DO JW=1,NWB
    IEGR = US(BS(JW))-1
    WRITE (INI,12530) JW
    DO WHILE (IEGR < DS(BE(JW)))
      IBGR = IEGR+1
      IEGR = IEGR+18
      IF (IEGR > DS(BE(JW))) IEGR = DS(BE(JW))
      WRITE (INI,'(/7X,21I10)') (I,I=IBGR,IEGR)
      WRITE (INI,'(7X,21I10)')  (INT(XI(I)),I=IBGR,IEGR)
      DO K=1,KMX
        DO I=IBGR,IEGR
          WRITE (BK(I),'(F8.0)') B(K,I)
        END DO
        IF (K == KTWB(JW)) THEN
          WRITE (INI,'(1X,I3,1X,"KT",1X,25A10)') K,(BK(I),I=IBGR,IEGR)
        ELSE
          WRITE (INI,'(1X,I3,4X,25A10)')         K,(BK(I),I=IBGR,IEGR)
        END IF
      END DO
    END DO
    CALL PRINT_INITIAL_GRID (KTWB(JW))
  END DO

  CLOSE (INI,ERR=10)

! Close error/warning files

10 CONTINUE
  IF (DELETE_WRN) THEN
    CLOSE (WRN,STATUS='DELETE')
    WRITE (WIN,'(/A)') 'Number of warnings = 0'
  ELSE
    CLOSE (WRN)
    WRITE (WIN,'(/A,I0,A)') 'Number of warnings = ',NWRN,' - see the file pre.wrn for more information'
  END IF
  IF (DELETE_ERR) THEN
    CLOSE (ERR,STATUS='DELETE')
    WRITE (WIN,FMTA)  'Number of errors   = 0'
  ELSE
    CLOSE (ERR)
    WRITE (WIN,FMTI)  'Number of errors   = ',NERR,' - see the file pre.err for more information'
  END IF
  WRITE (WIN,FMTA)
  WRITE (WIN,FMTA) 'Normal termination'
  GO TO 99999

! Error message for input cards
  
405 CONTINUE
    WRITE(ERR,'(A)')'Severe error in branch geometry or waterbody-branch specification. Pre-processing stopped.'
    WRITE(WIN,'(A)')'Severe error in branch geometry or waterbody-branch specification. Pre-processing stopped. See pre.err for information.'
    CLOSE(ERR);STOP
  
402 CONTINUE
    CALL INPUT_ERROR (UNIT,ERR,AID,'Either illegal value or incorrect card somewhere in the following cards')
    CLOSE(ERR);STOP
    
  
401 CONTINUE
    WRITE(ERR,'(A,A)')'Read Error in file:',AID
    WRITE(ERR,'(A)')'Check to make sure that the file is in Windows file format <CR><LF> and not Unix file format <CR>.'
    CALL INPUT_ERROR (UNIT,ERR,AID,'Either illegal value or incorrect card somewhere in the following cards')
    CLOSE(ERR)
    STOP

400 CONTINUE
    IF(ICST==2)THEN
    CALL INPUT_ERROR (UNIT,ERR,'GRAPH.NPT','Either illegal value or incorrect card somewhere in the following cards in GRAPH.NPT')
    ELSEIF(ICST==3)THEN
    CALL INPUT_ERROR2 (UNIT,ERR,AID,AID1,'Header value was not expected-error in file row order or the header value was incorrect   ')   
    ELSE
    CALL INPUT_ERROR (UNIT,ERR,AID,'Either illegal value or incorrect card somewhere in the following cards')
    ENDIF
  IF(ICST == 1)THEN
  WRITE(ERR,*)'Probable error in active constituent list - in w2_con.npt and graph.npt. The list of constituents should be:'
  WRITE(ERR,*)'#     PROPER_CONST_ORDER      CONST_IN_CONTROL_FILE'
  DO JC=1,NCT
  WRITE(ERR,'(I3,5X,A,4X,A)')JC,CCONSTIT(JC),CNAME2(JC)
  ENDDO
  ENDIF

  IF(ICST==2)THEN
      WRITE(ERR,'(A)')'Probable error in order of active constituents - in graph.npt. The list of constituents should be:'
      WRITE(ERR,*)'#     PROPER_CONST_ORDER      CONST_IN_W2_CON.NPT   CONST_IN_GRAPH.NPT'
      DO JC=1,NCT
      WRITE(ERR,'(I3,5X,A,4X,A,16X,A)')JC,CCONSTIT(JC),CNAME2(JC),CNAME(JC)
  ENDDO
  ENDIF

  CLOSE(ERR)
  STOP
!**** Output formats

11130 FORMAT (1X,A)
11140 FORMAT (/1X,A/                                                                                                               &
              3X,'Starting time (Julian day) [TMSTRT] =',F8.2/                                                                     &
              3X,'Ending time (Julian day)    [TMEND] =',F8.2/                                                                     &
              3X,'Year                         [YEAR] =',I8)
11150 FORMAT (3X,'# Timestep intervals         [NDLT] =',I8/                                                                       &
              3X,'Minimum timestep (sec)     [DLTMIN] =',F8.1)
11160 FORMAT (3X,'Timestep day (Julian day)    [DLTD] =',F8.1)
11170 FORMAT (3X,'Maximum timestep (sec)     [DLTMAX] =',7F8.1,' sec')
11180 FORMAT (3X,'Fraction of timestep         [DLTF] =',7F8.2)
11190 FORMAT (3X,A/                                                                                                                &
             (5X,'Waterbody ',I0/                                                                                                  &
              5X,'  Vertical eddy viscosity    [VISC] = ',A/                                                                       &
              5X,'  Internal gravity wave      [CELC] = ',A/                                                                       &
              5X,'  Additional stability     [DLTADD] = ',A))
11250 FORMAT (/1X,A/                                                                                                               &
             (3X,'Waterbody ',I0/                                                                                                  &
              5X,'Evaporation     [EVC] = ',A/                                                                                     &
              5X,'Precipitation   [PRC] = ',A/                                                                                     &
              5X,'Volume balance  [VBC] = ',A/                                                                                     &
              5X,'Energy balance  [EBC] = ',A/                                                                                     &
              5X,'Mass balance    [MBC] = ',A/                                                                                     &
              5X,'Place inflows   [PQC] = ',A/                                                                                     &
              5X,'Wind          [WINDC] = ',A/                                                                                     &
              5X,'Inflow         [QINC] = ',A/                                                                                     &
              5X,'Outflow       [QOUTC] = ',A/                                                                                     &
              5X,'Heat exchange [HEATC] = ',A/                                                                                     &
              5X,'Heat exchange [SLHTC] = ',A))
11262 FORMAT (3X,'Waterbody ',I0/                                                                                                  &
              5X,'read radiation               [SROC] =',7(5X,A):/                                                                 &
             (42X,7(5X,A)))
11263 FORMAT (5X,'wind function coefficient a   [AFW] =',7F8.2:/                                                                   &
             (42X,7(F8.2)))
11264 FORMAT (5X,'wind function coefficient b   [BFW] =',7F8.2:/                                                                   &
             (42X,7(F8.2)))
11265 FORMAT (5X,'wind function coefficient c   [CFW] =',7F8.2:/                                                                   &
             (42X,7(F8.2)))
11266 FORMAT (5X,'wind height                 [WINDH] =',7F8.2:/                                                                   &
             (42X,7(F8.2)))
11267 FORMAT (5X,'Ryan-Harleman evaporation   [RHEVC] =',7(5X,A):/                                                                 &
             (42X,7(5X,A)))
11270 FORMAT (/1X,'Input Interpolations'/                                                                                          &
             (3X,'Branch ',I0/                                                                                                     &
              5X,'Inflow                [QINIC] = ',A/                                                                             &
              5X,'Distributed tributary [DTRIC] = ',A/                                                                             &
              5X,'Head boundary          [HDIC] = ',A))
11271 FORMAT (3X,'Waterbody ',I0/                                                                                                  &
              5X,'Meteorology           [METIC] = ',A)
11272 FORMAT (:3X,'Tributary ',I0,'              [TRIC] = ',A)
11273 FORMAT (3X,'Branch ',I0/                                                                                                     &
             (5X,'Structure ',I0,'           [STRIC] = ',A))
11274 FORMAT (:5X,'Withdrawal ',I0,'            [WDIC] = ',A)
11300 FORMAT (3X,'Waterbody ',I0/                                                                                                  &
              5X,'Latitude                  [LAT] =',F8.2/                                                                         &
              5X,'Longitude                [LONG] =',F8.2)
11330 FORMAT (3X,'Axis orientation')
11340 FORMAT (5X,'Segment #    ',(1X,21I5))
11350 FORMAT (5X,'[PHI0] (rads)',(1X,21F5.2))
11360 FORMAT (/1X,'Transport Solution'/                                                                                            &
             (3X,'Waterbody ',I0/                                                                                                  &
              5X,'Transport [SLTRC] = ',A/                                                                                         &
              5X,'Theta     [THETA] =',F9.2))
11365 FORMAT (/1X,'Hydraulic coefficients')
11370 FORMAT (3X,'Longitudinal eddy viscosity           [AX] =',F9.2,' m^2/sec'/                                                   &
              3X,'Longitudinal eddy diffusivity         [DX] =',F9.2,' m^2/sec'/                                                   &
              3X,'Sediment temperature                [TSED] =',F9.2,1X,'oC'/                                                         &
              3X,'Coefficient of bottom heat exchange [CBHE] =',1PE9.2,1X,'W m^2/oC')
11375 FORMAT (/1X,'Ice cover')
11380 FORMAT (3X,'Waterbody ',I0/                                                                                                  &
              5X,'Ice calculations          [ICEC] = ',5X,A/                                                                       &
              5X,'Solution                [SLICEC] = ',A/                                                                          &
              5X,'Albedo                  [ALBEDO] = ',F8.2/                                                                       &
              5X,'Ice-water heat exchange    [HWI] = ',F8.2/                                                                       &
              5X,'Light absorption         [BETAI] = ',F8.2/                                                                       &
              5X,'Light decay             [GAMMAI] = ',F8.2)
11385 FORMAT (/1X,'Output Control')
11390 FORMAT (3X,'Waterbody ',I0/                                                                                                  &
             (5X,A41,' = ',A))
11400 FORMAT (3X,'Waterbody ',I0)
11410 FORMAT (5X,'Snapshot         [SNPC] = ',A)
11420 FORMAT (7X,'Number of time intervals [NSNP] =',I7)
11430 FORMAT (7X,'Date  (Julian day)       [SNPD] =',11F7.2)
11440 FORMAT (7X,'Frequency  (days)        [SNPF] =',11F7.2)
11450 FORMAT (5X,'Screen           [SCRC] = ',A)
11455 FORMAT (7X,'Number of time intervals [NSCR] =',I7)
11460 FORMAT (7X,'Date  (Julian day)       [SCRD] =',11F7.2)
11470 FORMAT (7X,'Frequency  (days)        [SCRF] =',11F7.2)
11480 FORMAT (5X,'Time series      [TSRC] = ',A)
11490 FORMAT (7X,'Number of time intervals [NTSR] =',I7)
11500 FORMAT (7X,'Date  (Julian day)       [TSRD] =',11F7.2)
11510 FORMAT (7X,'Frequency  (days)        [TSRF] =',11F7.2)
11520 FORMAT (5X,'Fluxes           [FLXC] = ',A)
11530 FORMAT (7X,'Number of time intervals [NFLX] =',I7)
11540 FORMAT (7X,'Date  (Julian day)       [FLXD] =',11F7.2)
11550 FORMAT (7X,'Frequency  (days)        [FLXF] =',11F7.2)
11560 FORMAT (5X,'Vector plot      [VPLC] = ',A)
11570 FORMAT (7X,'Number of time intervals [NVPL] =',I7)
11580 FORMAT (7X,'Date  (Julian day)       [VPLD] =',11F7.2)
11590 FORMAT (7X,'Frequency  (days)        [VPLF] =',11F7.2)
11600 FORMAT (5X,'Profile plot     [PRFC] = ',A)
11610 FORMAT (7X,'Number of time intervals [NPRF] =',I7)
11620 FORMAT (7X,'Number of stations      [NIPRF] =',I7)
11630 FORMAT (7X,'Segment location         [IPRF] =',11I7)
11640 FORMAT (7X,'Date  (Julian day)       [PRFD] =',11F7.2)
11650 FORMAT (7X,'Frequency  (days)        [PRFF] =',11F7.2)
11660 FORMAT (5X,'Spreadsheet plot [SPRC] = ',A)
11670 FORMAT (7X,'Number of time intervals [NSPR] =',I7)
11680 FORMAT (7X,'Number of stations      [NISPR] =',I7)
11690 FORMAT (7X,'Segment location         [ISPR] =',11I7)
11700 FORMAT (7X,'Date  (Julian day)       [SPRD] =',11F7.2)
11710 FORMAT (7X,'Frequency  (days)        [SPRF] =',11F7.2)
11720 FORMAT (5X,'Contour plot     [CPLC] = ',A)
11730 FORMAT (7X,'Number of time intervals [NCPL] =',I7)
11740 FORMAT (7X,'Date (Julian day)        [CPLD] =',11F7.2)
11750 FORMAT (7X,'Frequency  (days)        [CPLF] =',11F7.2)
11760 FORMAT (5X,'Restart out      [RSOC] = ',A/                                                                                   &
              5X,'Restart in       [RSIC] = ',A)
11770 FORMAT (7X,'Number of time intervals [NRSO] =',I7)
11780 FORMAT (7X,'Date (Julian day)        [RSOD] =',11F7.2)
11790 FORMAT (7X,'Frequency  (days)        [RSOF] =',11F7.2)
11800 FORMAT (/'Inflow/Outflow'/                                                                                                   &
              3X,'Selective Withdrawal'/                                                                                           &
              5X,'Branch ',3X,'# of structures [NSTR]':/                                                                           &
             (7X,I2,14X,I3))
11805 FORMAT (5X,'Branch ',I0)
11810 FORMAT (:7X,'Structure',3X,'Type',3X,'Width (m)',3X,'Elevation (m)',3X,'Bottom Layer'/                                       &
             (10X,I2,4X,A8,1X,F8.1,5X,F8.1,9X,I5))
11820 FORMAT (5X,'Gate coefficients'/                                                                                              &
              7X,'[A1GT]',4X,'[B1GT]',4X,'[G1GT]',4X,'[A2GT]',4X,'[B2GT]',4X,'[G2GT]'/                                             &
              :(2X,6(2X,F8.2)))
11821 FORMAT (3X,'Number of spillways/external weirs [NSP] = ',I0:/                                                                &
              5X,'centerline elevation             [ESP] = ',7F8.2,                                                                &
             (:/T26,7F8.2))
11822 FORMAT (5X,'upstream segment number         [IUSP] = ',7I8,                                                                  &
             (:/T26,40I3))
11823 FORMAT (5X,'downstream segment number       [IDSP] = ',7I8,                                                                  &
             (:/T26,40I3))
11824 FORMAT (5X,'upstream inflow placement      [PUSPC] = ',7A8,                                                                  &
             (:T33,7A8))
11825 FORMAT (5X,'downstream inflow placement    [PDSPC] = ',7A8,                                                                  &
             (:T33,7A8))
11826 FORMAT (5X,'upstream top elevation         [EUTSP] = ',7F8.2,                                                                &
             (:T33,7F8.2))
11827 FORMAT (5X,'downstream top elevation       [EDTSP] = ',7F8.2,                                                                &
             (:T33,7F8.2))
11828 FORMAT (5X,'upstream bottom elevation      [EUBSP] = ',7F8.2,                                                                &
             (:T33,7F8.2))
11829 FORMAT (5X,'downstream bottom elevation    [EDBSP] = ',7F8.2,                                                                &
             (:T33,7F8.2))
11830 FORMAT (5X,'Spillway/external weir coefficients'/                                                                            &
              7X,'[A1SP]',4X,'[B1SP]',4X,'[A2SP]',4X,'[B2SP]'/                                                                     &
              :(2X,4(2X,F8.2)))
             
11831 FORMAT (3X,'Number of pipes [NPI] = ',I0:/                                                                                   &
              5X,'upstream centerline elevation   [EUPI] = ',7F8.2,                                                                &
             (:/T26,7F8.2))
11841 FORMAT (5X,'downstream centerline elevation [EDPI] = ',7F8.2,                                                                &
             (:/T26,7F8.2))
11842 FORMAT (5X,'diameter                         [WPI] = ',7F8.2,                                                                &
             (:/T26,7F8.2))
11843 FORMAT (5X,'length                         [DLXPI] = ',7F8.2,                                                                &
             (:/T26,7F8.2))
11844 FORMAT (5X,'Mannings n                       [FPI] = ',7F8.2,                                                                &
             (:/T26,7F8.2))
11845 FORMAT (5X,'minor friction losses         [FMINPI] = ',7F8.2,                                                                &
             (:/T26,7F8.2))
11832 FORMAT (5X,'upstream segment number         [IUPI] = ',7I8,                                                                  &
             (:/T26,40I3))
11833 FORMAT (5X,'downstream segment number       [IDPI] = ',7I8,                                                                  &
             (:/T26,40I3))
11834 FORMAT (5X,'upstream inflow placement      [PUPIC] = ',7A8,                                                                  &
             (:T33,7A8))
11835 FORMAT (5X,'downstream inflow placement    [PDPIC] = ',7A8,                                                                  &
             (:T33,7A8))
11836 FORMAT (5X,'upstream top elevation         [EUTPI] = ',7F8.2,                                                                &
             (:T33,7F8.2))
11837 FORMAT (5X,'downstream top elevation       [EDTPI] = ',7F8.2,                                                                &
             (:T33,7F8.2))
11838 FORMAT (5X,'upstream bottom elevation      [EUBPI] = ',7F8.2,                                                                &
             (:T33,7F8.2))
11839 FORMAT (5X,'downstream bottom elevation    [EDBPI] = ',7F8.2,                                                                &
             (:T33,7F8.2))
11840 FORMAT (3X,'Number of withdrawals [NWD] = ',I0:/                                                                             &
              5X,'segment number [IWD] = ',10I3,                                                                                   &
             (:/T26,10I3))
11851 FORMAT (3X,'Number of pumps [NPU] = ',I0:/                                                                                   &
              5X,'centerline elevation             [EPU] = ',7F8.2,                                                                &
             (:/T26,7F8.2))
11852 FORMAT (5X,'upstream segment number         [IUPU] = ',7I8,                                                                  &
             (:/T26,7I8))
11853 FORMAT (5X,'downstream segment number       [IDPU] = ',7I8,                                                                  &
             (:/T26,7I8))
11854 FORMAT (5X,'starting time                 [STRTPU] = ',7F8.2,                                                                &
             (:T33,7F8.2))
11855 FORMAT (5X,'ending time                    [ENDPU] = ',7F8.2,                                                                &
             (:T33,7F8.2))
11856 FORMAT (5X,'starting elevation             [EONPU] = ',7F8.2,                                                                &
             (:T33,7F8.2))
11857 FORMAT (5X,'ending elevaton               [EOFFPU] = ',7F8.2,                                                                &
             (:T33,7F8.2))
11858 FORMAT (5X,'pump flowrate                    [QPU] = ',7F8.2,                                                                 &
             (:T33,7F8.2))
11859 FORMAT (5X,'inflow placement                [PPUC] = ',7A8,                                                                  &
             (:T33,7A8))
11861 FORMAT (5X,'top elevation                   [ETPU] = ',7F8.2,                                                                &
             (:T33,7F8.2))
11862 FORMAT (5X,'bottom elevation                [EBPU] = ',7F8.2,                                                                &
             (:T33,7F8.2))
11863 FORMAT (5X,'top layer                       [KTPU] = ',7I8,                                                                  &
             (:T33,7I8))
11864 FORMAT (5X,'bottom layer                    [KBPU] = ',7I8,                                                                  &
             (:T33,7I8))
11860 FORMAT (3X,'Number of tributaries [NTR] = ',I0:/                                                                             &
              5X,'segment number     [ITR] = ',7I8,                                                                                &
             (:/T33,7I8))
             
11865 FORMAT (5X,'Spillway/external gas coefficients'/                                                                            &
              8X,'[GASSPC]',2X,'[EQSP]',3x,'[AGASSP]',2X,'[BGASSP]',3X,'[CGASSP]'/                                                                     &
              :(2X,5X,A8,5x,i2,2x,3(2X,F8.2)))
11866 FORMAT (5X,'Gate/external gas coefficients'/                                                                            &
              8X,'[GASGTC]',2X,'[EQGT]',3x,'[AGASGT]',2X,'[BGASGT]',3X,'[CGASGT]'/                                                                     &
              :(2X,5X,A8,5x,i2,2x,3(2X,F8.2)))

11870 FORMAT (:5X,'Inflow placement  [PTRC] = ',7A8,                                                                               &
             (/:T33,7A8))
11880 FORMAT (:5X,'Top elevation     [ETTR] = ',7F8.2,                                                                             &
             (:/T33,7F8.2))
11890 FORMAT (:5X,'Bottom elevation  [EBTR] = ',7F8.2,                                                                             &
             (:/T33,7F8.2))
11891 FORMAT (3X,'Number of gates [NGT] =',I3:/                                                                                    &
              5X,'centerline elevation             [EGT] = ',7F8.2,                                                                &
             (:/T26,7F8.2))
11892 FORMAT (:5X,'upstream segment number         [IUGT] = ',7I8,                                                                 &
             (:/T26,40I3))
11893 FORMAT (:5X,'downstream segment number       [IDGT] = ',7I8,                                                                 &
             (:/T26,40I3))
11894 FORMAT (:5X,'upstream inflow placement      [PUGTC] = ',7A8,                                                                 &
             (:T33,7A8))
11895 FORMAT (:5X,'downstream inflow placement    [PDGTC] = ',7A8,                                                                 &
             (:T33,7A8))
11896 FORMAT (:5X,'upstream top elevation         [EUTGT] = ',7F8.2,                                                               &
             (:T33,7F8.2))
11897 FORMAT (:5X,'downstream top elevation       [EDTGT] = ',7F8.2,                                                               &
             (:T33,7F8.2))
11898 FORMAT (:5X,'upstream bottom elevation      [EUBGT] = ',7F8.2,                                                               &
             (:T33,7F8.2))
11899 FORMAT (:5X,'downstream bottom elevation    [EDBGT] = ',7F8.2,                                                               &
             (:T33,7F8.2))
11900 FORMAT (3X,'Distributed tributaries [DTRC]')
11910 FORMAT (5X,'Branch ',I0,' = ',A3)

11930 FORMAT (/1X,'Input Filenames'/3X,'Control    = ',A72/3X,'Restart    = ',A72/3X,'Withdrawal = ',A72)
11940 FORMAT (3X,'Waterbody ',I0/                                                                                                  &
              5X,'Bathymetry           = ',A72/                                                                                    &
              5X,'Meteorology          = ',A72/                                                                                    &
              5X,'Vertical profile     = ',A72/                                                                                    &
              5X,'Longitudinal profile = ',A72)
11950 FORMAT (3X,'Branch ',I0/                                                                                                     &
              5X,'Inflow                               = ',A72/                                                                    &
              5X,'Inflow temperature                   = ',A72/                                                                    &
              5X,'Inflow concentrations                = ',A72/                                                                    &
              5X,'Outflow                              = ',A72/                                                                    &
              5X,'Distributed tributary inflows        = ',A72/                                                                    &
              5X,'Distributed tributary temperatures   = ',A72/                                                                    &
              5X,'Distributed tributary concentrations = ',A72/                                                                    &
              5X,'Precipitation                        = ',A72/                                                                    &
              5X,'Precipitation temperatures           = ',A72/                                                                    &
              5X,'Precipitation concentrations         = ',A72/                                                                    &
              5X,'Upstream head                        = ',A72/                                                                    &
              5X,'Upstream head temperatures           = ',A72/                                                                    &
              5X,'Upstream head concentrations         = ',A72/                                                                    &
              5X,'Downstream head                      = ',A72/                                                                    &
              5X,'Downstream head temperatures         = ',A72/                                                                    &
              5X,'Downstream head concentrations       = ',A72)
11970 FORMAT (:(3X,'Tributary ',I0/                                                                                                &
              5X,'Inflow               = ',A72/                                                                                    &
              5X,'Inflow temperature   = ',A72/                                                                                    &
              5X,'Inflow concentration = ',A72)/)
11980 FORMAT (1X,'Output Filenames'/                                                                                               &
              3X,'Error       = ',A72/                                                                                             &
              3X,'Warning     = ',A72/                                                                                             &
              3X,'Time series = ',A72/                                                                                             &
              3X,'Withdrawal  = ',A72)
11990 FORMAT (3X,'Waterbody ',I0/                                                                                                  &
              5X,'Snapshot     = ',A72/                                                                                            &
              5X,'Fluxes       = ',A72/                                                                                            &
              5X,'Profile      = ',A72/                                                                                            &
              5X,'Vector plot  = ',A72/                                                                                            &
              5X,'Contour plot = ',A72)
12000 FORMAT (/1X,'Constituents [CCC] = ',A3/                                                                                      &
              3X,'Algal limiting nutrient  [LIMC] = ',A3/                                                                          &
              3X,'Kinetics update frequency [CUF] = ',I3)
12020          FORMAT (5X,'State Variables'/                                                                                                &
             ! 7X,'Constituent',T29,'Computation',T42,'Initial Conc',T57,'Fluxes',T67,'Printout'/                                   &
             ! 9X,'[CNAME]',T32,'[CAC]',T42,'[C2IWB,g/m^3]',T57,'[CFWBC]',T67,'[CPRWBC]')
              7X,'Constituent',T29,'Computation',T42,'Initial Conc',T57,'Printout'/                                   &
              9X,'[CNAME]',T32,'[CAC]',T42,'[C2IWB,g/m^3]',T57,'[CPRWBC]')
12030 FORMAT (/3X,'Waterbody ',I0)
12040 FORMAT (7X,A19,T33,A3,T42,F10.3,T59,A3)   !,T69,A3)
12045 FORMAT (/3X,'Branch ',I0/                                                                                                    &
              5X,'State Variables'/                                                                                                &
              7X,'Constituent',T32,'Inflow',T42,'Distributed trib',T62,'Precipitation'/                                            &
              9X,'[CNAME]',T31,'[CINBRC]',T46,'[CDTBRC]',T64,'[CPRBRC]'/                                                           &
             (7X,A19,T33,A3,T48,A3,T67,A3))
12046 FORMAT (/3X,'Tributary ',I0/                                                                                                 &
              5X,'State Variables'/                                                                                                &
              7X,'Constituent',T32,'Inflow'/                                                                                       &
              9X,'[CNAME]',T31,'[CINTRC]'/                                                                                         &
             (7X,A19,T33,A3))
12050 FORMAT (/5X,'Derived Variables')
12060 FORMAT (19X,'Constituent',T46,'Computation'/                                                                                 &
              20X,'[CDNAME]',T48,'[CDWBC]')
12070 FORMAT (7X,A,T50,1000A3)
12080 FORMAT (5X,A/                                                                                                                &
              7X,'Constituent',13X,'Rate/Coefficient')
12090 FORMAT (7X,'Suspended solids',T27,'Settling',T50,'[SSS] =',F6.3,' m/day')
12100 FORMAT (7X,A21,T27,'Temperature mult',T48,'[CGQ10] =',F6.3/                                                                  &
              T27,'0-Order Decay        [CG0DK] =',F6.3,'g/m^3/day'/                                                               &
              T27,'1-Order Decay        [CG1DK] =',F6.3,'/day'/                                                                    &
              T27,'Settling               [CGS] =',F6.3,' m/day')
12110 FORMAT (7X,'Labile DOM',T27,'Decay',T47,'[LDOMDK] =',F6.3,' /day'/                                                           &
              T27,'to refractory        [LRDDK] =',F6.3,' /day')
12120 FORMAT (7X,'Refractory DOM',T27,'Decay',T47,'[RDOMDK] =',F6.3,' /day')
12130 FORMAT (7X,'Labile POM',T27,'Decay',T47,'[LPOMDK] =',F6.3,' /day'/                                                           &
              T27,'to refractory        [LRPDK] =',F6.3,' /day'/                                                                   &
              T27,'Settling',T49,'[POMS] =',F6.3,' m/day')
12150 FORMAT (7X,'Refractory POM',T27,'Decay',T47,'[RPOMDK] =',F6.3,' /day')
12160 FORMAT (7X,'Algal group',I3,T27,'Growth',T51,'[AG] =',F6.3,' /day'/                                                          &
              T27,'Mortality               [AM] =',F6.3,' /day'/                                                                   &
              T27,'Excretion               [AE] =',F6.3,' /day'/                                                                   &
              T27,'Respiration             [AR] =',F6.3,' /day'/                                                                   &
              T27,'Settling                [AS] =',F6.3,' m/day'/                                                                  &
              T27,'Org-P                 [ALGP] =',F6.3/                                                                           &
              T27,'Org-N                 [ALGN] =',F6.3/                                                                           &
              T27,'Org-C                 [ALGC] =',F6.3/                                                                           &
              T27,'Org-Si               [ALGSI] =',F6.3/                                                                           &
              T27,'Algae/chl a ratio    [ACHLA] =',F6.2,' mg/ug'/                                                                  &
              T27,'Fraction algae to POM [APOM] =',F6.2/                                                                           &       
              T27,'Vertical Migration  [AVERTM] =',A8) 
12161         FORMAT (7X,'Zooplankton group',I0,T27,'Growth',T51,'[ZG] =',F6.3,' /day'/                                                    &
              T27,'Mortality               [ZM] =',F6.3,' /day'/                                                                   &
              T27,'Respiration             [ZR] =',F6.3,' /day'/                                                                   &
              T27,'Settling                [ZS] =',F6.3,' m/day'/                                                                  &
              T27,'Org-P                   [ZP] =',F6.3/                                                                           &
              T27,'Org-N                   [ZN] =',F6.3/                                                                           &
              T27,'Org-C                   [ZC] =',F6.3)    
12162 FORMAT (7X,'Periphyton group',I0,T27,'Growth',T51,'[EG] =',F6.3,' /day'/                                                          &
              T27,'Mortality               [EM] =',F6.3,' /day'/                                                                   &
              T27,'Excretion               [EE] =',F6.3,' /day'/                                                                   &
              T27,'Respiration             [ER] =',F6.3,' /day'/                                                                   &
              T27,'Burial                  [EB] =',F6.3,' m/day'/                                                                  &
              T27,'Org-P                   [EP] =',F6.3/                                                                           &
              T27,'Org-N                   [EN] =',F6.3/                                                                           &
              T27,'Org-C                   [EC] =',F6.3/                                                                           &
              T27,'Org-Si                 [ESI] =',F6.3/                                                                           &
              T27,'Periphyton/chla ratio[ECHLA] =',F6.2,' mg/ug'/                                                                  &
              T27,'Fraction Periphyton to POM [EPOM] =',F6.2/)      ! MLM 7/11/06
12180 FORMAT (7X,'Phosphorous',T27,'Release',T46,'[PO4R] =',F6.3,' g/m^2/day')
12190 FORMAT (7X,'Ammonium',T27,'Decay',T45,'[NH4DK] =',F6.3,' /day'/                                                              &
              T27,'Release',T46,'[NH4R] =',F6.3,' g/m^2/day'/T27,'NH3 gas loss',T40,'[KG_H2O_CONSTANT]=',F7.1)
12200 FORMAT (7X,'Nitrate-Nitrite',T27,'Decay',T45,'[NO3DK] =',F6.3,' /day'/   &
              T27,'Diffusion to sed',T45,'[NO3S] =',F6.3,' m/day'/                           &
              T27,'Frac to OM',T45,'[FNO3SED] =',F6.3,' [-]')
12210 FORMAT (7X,'Silica',T27,'Decay',T45,'[PSIDK] =',F6.3,' /day'/                                                                &
              T27,'Release',T46,'[DSIR] =',F6.3,' g/m^2/day'/                                                                      &
              T27,'Settling',T46,'[PSIS] =',F6.3,' m/day')
12220 FORMAT (7X,'Sediment',T27,'Decay',T45,'[SEDDK] =',F6.3,' /day')
12240 FORMAT (7X,'Oxygen',T27,'Sediment demand     [SOD] = ',13F5.1:/                                                              &
             (T55,13F5.1))
12250 FORMAT (T27,'SOD fraction',T46,'[FSOD] = ',F5.1/                                                                             &
              T27,'Sediment fraction',T46,'[FSED] = ',F5.1)
12255 FORMAT (T27,'Sediment focusing ',T46,'[SEDS] = ',F7.3,' m/day'/                                                                             &
              T27,'Sediment burial ',T46,'[SEDBR] = ',F7.3,' /day')
12270 FORMAT (5X,A/                                                                                                                &
              9X,'Constituent',T24,'Rate',T37,'Lower',T50,'Max Lower')
12280 FORMAT (7X,'Algal group ',I0,T23,'Growth',T34,'[AT1] =',F5.1,T49,'[AT2] =',F5.1)
12290 FORMAT (7X,'Ammonium',T23,'Decay',T32,'[NH4T1] =',F5.1,T47,'[NH4T2] =',F5.1)
12300 FORMAT (7X,'Nitrate',T23,'Decay',T32,'[NO3T1] =',F5.1,T47,'[NO3T2] =',F5.1)
12310 FORMAT (7X,'Organic',T23,'Decay',T33,'[OMT1] =',F5.1,T48,'[OMT2] =',F5.1)
12320 FORMAT (7X,'Sediment',T23,'Decay',T32,'[SEDT1] =',F5.1,T47,'[SEDT2] =',F5.1)
12330 FORMAT (5X,A/                                                                                                                &
              9X,'Constituent',T24,'Rate',T37,'Upper',T50,'Max Upper')
12340 FORMAT (7X,'Algal group ',I0,T23,'Growth',T34,'[AT3] =',F5.1,T49,'[AT4] =',F5.1)
12350 FORMAT (5X,A/                                                                                                                &
              7X,'Oxygen'/                                                                                                         &
              9X,'Ammonium       [O2NH4] =',F5.2/                                                                                  &
              9X,'Organic matter  [O2OM] =',F5.2/)
12355 FORMAT (5X,A/(7X,'Algal group ',I0/9X,'Respiration     [O2AR] =',F5.2/                                                       &
              9X,'Algal growth    [O2AG] =',F5.2))
12356 FORMAT (5X,A/(7X,'Periphyton group ',I0/      9X,'Respiration      [O2ER] =',F5.2/                                                       &
              9X,'Periphyton growth [O2EG] =',F5.2))
12360 FORMAT (7X,'Organic Matter'/                                                                                                 &
              9X,'Carbon      [BIOC] =',F6.3/                                                                                      &
              9X,'Phosphorous [BIOP] =',F6.3/                                                                                      &
              9X,'Nitrogen    [BION] =',F6.3/                                                                                      &
              9X,'Silica     [BIOSI] =',F6.3)
12370 FORMAT (:5X,'CBOD ',I0/                                                                                                      &
              7X,'Decay rate                   [KBOD] =',F6.3,' /day'/                                                             &
              7X,'Temperature adjustment       [TBOD] =',F6.3/                                                                     &
              7X,'Ultimate CBOD to CBOD5 ratio [RBOD] =',F6.3)
12380 FORMAT (5X,A15/                                                                                                              &
             (7X,'Algal group ',I0/                                                                                                &
              9X,'Phosphorous [AHSP] =',F6.3,' g/m^3'/                                                                             &
              9X,'Nitrogen    [AHSN] =',F6.3,' g/m^3'/                                                                             &
              9X,'Silica     [AHSSI] =',F6.3,' g/m^3'))
12390 FORMAT (5X,A/                                                                                                                &
              9X,'Surface layer      [BETA] =',F5.2/                                                                               &
              9X,'Water             [EXH2O] =',F5.2,' /m'/                                                                         &
              9X,'Inorganic solids   [EXSS] =',F5.2,' /m'/                                                                         &
              9X,'Organic solids     [EXOM] =',F5.2,' /m')
12391 FORMAT (5X,A/                                                                                                                &        
              (9X,'Bacteria decay            [BACT1DK] =',F5.3 ' /day'/                                                                               &
              9X,'Theta Temp Rate Arrhenius [BACTQ10] =',F5.3/                                                                         &
              9X,'Photo degradation         [BACTLDK] =',F5.3,' m^2/J'/                                                                         &
              9X,'Settling rate               [BACTS] =',F5.3,' m/day'))
12392 FORMAT (5X,A/                                                                                                                &        
              (9X,'H2S sediment release         [H2SR] =',F5.3 '[-]'/                                                                               &
              9X,'Theta Temp Rate Arrhenius  [H2SQ10] =',F5.3/                                                                         &
              9X,'H2S decay rate             [H2S1DK] =',F5.3,' /day'/                                                                         &
              9X,'SO4 sediment release         [SO4R] =',F5.3,'[-]'))
12393 FORMAT (5X,A/                                                                                                                &        
              (9X,'CH4 sediment release         [CH4R] =',F5.3 '[-]'/                                                                               &
              9X,'Theta Temp Rate Arrhenius  [CH4Q10] =',F5.3/                                                                         &
              9X,'CH4 decay rate             [CH41DK] =',F5.3,' /day'))
12394 FORMAT (5X,A/                                                                                                                &        
              (9X,'FE2 sediment release        [FEIIR] =',F5.3 '[-]'/                                                                               &
              9X,'Rate oxidation           [KFE_OXID] =',F5.3,' /day'/                                                                         &
              9X,'Rate reduction            [KFE_RED] =',F5.3,' /day'/                                                                         &
              9X,'Half-saturation Fe3[KFEOOH_HALFSAT] =',F5.3,/                                                                         &
              9X,'FE Settling velocity     [FESETVEL] =',F5.3,' m/day'))
12395 FORMAT (5X,A/                                                                                                                &        
              (9X,'MN2 sediment release        [MNIIR] =',F5.3 '[-]'/                                                                               &
              9X,'Rate oxidation           [KMN_OXID] =',F5.3,' /day'/                                                                         &
              9X,'Rate reduction            [KMN_RED] =',F5.3,' /day'/                                                                         &
              9X,'Half-saturation Mn  [KMNO2_HALFSAT] =',F5.3,/                                                                         &
              9X,'MN Settling velocity     [MNSETVEL] =',F5.3,' m/day'))
12396 FORMAT (5X,A/                                                                                                           &        
              (9X,'Reaeration water body      [REAERC] =',A /                                                                               &
              9X,'Equation #                   [NEQN] =',I2/                                                                         &
              9X,'Dissolved gas pressure/O2   [DGPO2] =',F5.3,/                                                                         &
              9X,'Minimum reaeration          [MINKL] =',F5.3, ' m/day or 1/day'/))             
               
12400 FORMAT (5X,A/(9X,'Algal group ',I0,'       [EXA] =',F5.2,' /m'))
12410 FORMAT (5X,A/7X,'Saturation Intensity'/                                                                                           &
             (9X,'Algal group ',I0,' [ASAT] =',F6.1,' W/m^2'))
12420 FORMAT (5X,A/                                                                                                                &
              7X,'Oxygen          [DMO2] =',1PE10.3,' m^2/g'/                                                                      &
              7X,'Carbon dioxide [DMCO2] =',1PE10.3,' m^2/g')
12430 FORMAT (5X,A/                                                                                                                &
              7X,'Phosphorous [PARTP] =',F6.3,' m^3/g'/                                                                            &
              7X,'Silica     [PARTSI] =',F6.3,' m^3/g')
12440 FORMAT (5X,A/                                                                                                                &
              7X,'Half-saturation coeff anoxia          [KDO] =',F6.2,' g/m^3'/                                                    &
              7X,'CO2 sediment release fraction of SOD [CO2R] =',F6.2,' [-]')
12460 FORMAT (3X,'Waterbody ',I0/                                                                                                  &
              T7,'total inflows',T25,'total outflows'/                                                                             &
              T5,'average  maximum',T24,'average  maximum'/                                                                        &
              T4,F8.2,T13,F8.2,T23,F8.2,T32,F8.2/)
12470 FORMAT (5X,'Branch ',I0/                                                                                                     &
              7X,'Inflows'/                                                                                                        &
              14X,'total'/                                                                                                         &
              9X,'average  maximum'/                                                                                               &
              9X,F8.2,1X,F8.2/                                                                                                     &
              13X,'upstream',9X,'tributaries',5X,'distributed tributaries',5X,'precipitation'/                                     &
              9X,'average  maximum',3X,'average  maximum',6X,'average  maximum',6X,'average  maximum'/                             &
              8X,F8.2,1X,F8.2,2X,F8.2,1X,F8.2,5X,F8.2,1X,F8.2,5X,F8.2,1X,F8.2/                                                     &
              7X,'Outflows'/                                                                                                       &
              14X,'outlets',9X,'withdrawals'/                                                                                      &
              9X,'average  maximum',3X,'average  maximum'/                                                                         &
              8X,F8.2,1X,F8.2,2X,F8.2,1X,F8.2)
12480 FORMAT (/A8/                                                                                                                 &
              3X,'Overall Grid'/                                                                                                   &
              5X,'Total'/7X,'segments [IMX] = ',I0/                                                                                &
              7X,'layers   [KMX] = ',I0/                                                                                           &
              7X,'branches [NBR] = ',I0)
12490 FORMAT (/3X,'Waterbody ',I0/                                                                                                 &
              5X,'Segments                 = ',I0,'-',I0/                                                                          &
              5X,'Branches                 = ',I0,'-',I0/                                                                          &
              5X,'Bottom elevation [ELBOT] =',F9.2,' m'/                                                                           &
              5X,'Surface layer       [KT] =',I9/                                                                                  &
              5X,'Vertical spacing     [H]')
12500 FORMAT (7X,'Layer',6X,14I5)
12510 FORMAT (7X,'Height (m) ',14F5.1)
12520 FORMAT (5X,'Branch ',I0/                                                                                                     &
              7X,'Upstream segment',T31,'[US] = ',I0,T44,'Downstream ','segment       [DS] = ',I0/                                 &
              7X,'Upstream head segment [UHS] = ',I0,T44,'Downstream ','head segment [DHS] = ',I0/                                 &
              7X,'Minimum # of layers    [NL] = ',I0)
12530 FORMAT (8X,'Waterbody ',I0,' Bathymetry [B], m')
12600 FORMAT (/1X,'Initial Branch Volume [VOLB] = ',F0.1,' m^3'/)
12620 FORMAT (T34,'Branch ',I0,' Volume-Area-Elevation Table'/T34,'Note: Elevation is at top of layer'/)
12630 FORMAT (1X,'Layer',T12,'Elevation',T27,'Area',T39,'Volume',T48,'Active Cells',T62,'Average depth',T77,'Average width',T92,   &
                 'Residence time'/                                                                                                 &
              T15,'(m)',T24,'(1.0E6 m^2)',T37,'(1.0E6 m^3)',T66,'(m)',T82,'(m)',T96,'(days)'/)
12640 FORMAT (1X,'Layer',T12,'Elevation',T27,'Area',T39,'Volume',T48,'Active Cells',T62,'Average depth',T77,'Average width',/      &
              T15,'(m)',T24,'(1.0E6 m^2)',T37,'(1.0E6 m^3)',T66,'(m)',T82,'(m)'/)
12650 FORMAT (1X,I3,' KT',T12,F7.2,T24,F8.3,T36,F12.3,T50,I5,T62,F8.1,T78,F9.2,T93,F12.2)
12655 FORMAT (1X,I3,      T12,F7.2,T24,F8.3,T36,F12.3,T50,I5,T62,F8.1,T78,F9.2,T93,F12.2)
12670 FORMAT (/1X,'Waterbody ',I0,' Initial Volume [VOLG] = ',F0.1,' m^3'/)
12671 FORMAT (/1X,'**ALL WATERBODIES**  Initial Volume [VOLG] = ',F0.1,' m^3'/)
12690 FORMAT (T34,'Waterbody ',I0,' Volume-Area-Elevation Table'/T34,'Note: Elevation is at top of layer'/)
12691 FORMAT (T34,'All Waterbodies   Volume-Area-Elevation Table'/T34,'Note: Elevation is at top of layer'/)
99999 CONTINUE
      IF(NERR>0)THEN
          WINFO%TYPE = QWIN$SIZECURR
          RESULT     = SETWSIZEQQ   (WIN,WINFO)
          
           STAT = ABOUTBOXQQ ('CE-QUAL-W2 Preprocessor V4.3 - PSU Dept. Civil & Environmental Engineering S. Wells'C)   
         !  WRITE (WINTITLE,'(A29,A360)') 'CE-QUAL-W2 4.3 preprocessor. ',moddir 
           ! STAT = ABOUTBOXQQ ('CE-QUAL-W2 Preprocessor V4.3: PRE.ERR'C)   
            WRITE (WINTITLE,'(A14,a360)') 'PRE.ERR file: ',moddir
            WINTITLE=TRIM(WINTITLE)
     OPEN  (10,FILE='user',TITLE=WINTITLE)
 
  STATUS     = FLOODFILL    (X,Y,ONE)
  I          = SETBKCOLOR   (15)
  I          = SETTEXTCOLOR (ONE)
  I          = SETEXITQQ    (QWIN$EXITPERSIST)    
     OPEN(20,FILE='PRE.ERR',STATUS='OLD')
     DO J=1,1000
     READ(20,'(A)',END=9977)TITLE240
     WRITE(10,'(A)')TRIM(ADJUSTL(TITLE240))
     ENDDO
9977 CONTINUE 
     RESULT = CLICKMENUQQ (QWIN$TILE)    ! TILE THE MULTIPLE WINDOWS ONCE DONE
     CLOSE(20)
      ENDIF
          
END PROGRAM W2_PRE

!***********************************************************************************************************************************
!*                                         S U B R O U T I N E   P R I N T  I N I T I A L  G R I D                                **
!***********************************************************************************************************************************

SUBROUTINE PRINT_INITIAL_GRID (KT)
  USE CONTROL; USE GRIDC; USE EPIPHYTON; USE MACROPHYTEC

! Variable initialization

  CONV = '          '

! Water surface variables

  WRITE (INI,'(//5X,A18/)')   'Water Surface [Z], m'
  DO I=1,IEPR(JW)
    WRITE (CONV(1,I),'(F10.4)') Z(IPR(I,JW))
  END DO
  WRITE (INI,'(2X,1000I10)')  (IPR(I,JW),I=1,IEPR(JW))
  WRITE (INI,'(3X,1000A10/)') (CONV(1,I),I=1,IEPR(JW))
  IF (ICE_CALC(JW)) WRITE (INI,'(/2X,A24/3X,1000F10.4/)') 'Ice Thickness (m)', (ICETH(IPR(I,JW)),I=1,IEPR(JW))

! Temperatures

  IF (HPRWBC(4,JW) == '      ON') THEN
    DO I=1,IEPR(JW)
      DO K=KT,KB(IPR(I,JW))
        WRITE (CONV(K,I),'(F10.3)') T2(K,IPR(I,JW))
      END DO
    END DO
    WRITE (INI,'(/2X,A35/)')    HNAME(4)
    WRITE (INI,'(2X,1000I10)') (IPR(I,JW),I=1,IEPR(JW))
    DO K=KT,KEPR(JW)
      WRITE (INI,'(1X,I3,1000A10)') K,(CONV(K,I),I=1,IEPR(JW))
    END DO
  END IF

! Constituent concentrations

  DO J=1,NAC
    IF (CPRWBC(CN(J),JW) == '      ON') THEN
      DO I=1,IEPR(JW)
        DO K=KT,KB(IPR(I,JW))
          WRITE (CONV(K,I),'(F10.3)') C2(K,IPR(I,JW),CN(J))*CMULT(CN(J))
        END DO
      END DO
      WRITE (INI,'(/5X,A37/)')    CNAME(CN(J))
      WRITE (INI,'(2X,1000I10)') (IPR(I,JW),I=1,IEPR(JW))
      DO K=KT,KEPR(JW)
        WRITE (INI,'(1X,I3,1000A10)') K,(CONV(K,I),I=1,IEPR(JW))
      END DO
    END IF
  END DO
  
  ! EPIPHYTON/PERIPHYTON concentrations

  DO J=1,NEP
    IF (EPIC(JW,J) == '      ON') THEN
      DO I=1,IEPR(JW)
        DO K=KT,KB(IPR(I,JW))
          WRITE (CONV(K,I),'(F10.3)') EPD(K,IPR(I,JW),J)
        END DO
      END DO
      WRITE (INI,'(/5X,A37/)')    'Epiphyton/periphyton g/m^2'
      WRITE (INI,'(2X,1000I10)') (IPR(I,JW),I=1,IEPR(JW))
      DO K=KT,KEPR(JW)
        WRITE (INI,'(1X,I3,1000A10)') K,(CONV(K,I),I=1,IEPR(JW))
      END DO
    END IF
  END DO
  
   ! MACROPHYTE concentrations

  DO J=1,NMC
    IF (MACWBC(JW,J) == '      ON') THEN
      DO I=1,IEPR(JW)
        DO K=KT,KB(IPR(I,JW))
          WRITE (CONV(K,I),'(F10.3)') mac(k,IPR(I,JW),j)
        END DO
      END DO
      WRITE (INI,'(/5X,A37/)')    'Macrophyte concentration g/m^3'
      WRITE (INI,'(2X,1000I10)') (IPR(I,JW),I=1,IEPR(JW))
      DO K=KT,KEPR(JW)
        WRITE (INI,'(1X,I3,1000A10)') K,(CONV(K,I),I=1,IEPR(JW))
      END DO
    END IF
  END DO

! Derived constituent concentrations

  CALL DERIVED_CONSTITUENTS
  DO J=1,NADC(JW)
    IF (CDWBC(CDN(J,JW),JW) == '      ON') THEN
      DO I=1,IEPR(JW)
        DO K=KT,KB(IPR(I,JW))
          WRITE (CONV(K,I),'(F10.3)') CD(K,IPR(I,JW),CDN(J,JW))*cDMULT(CDN(J,JW))     ! RA 7/11/08
        END DO
      END DO
      WRITE (INI,'(/5X,A43/)')    CDNAME(CDN(J,JW))
      WRITE (INI,'(2X,1000I10)') (IPR(I,JW),I=1,IEPR(JW))
      DO K=KT,KEPR(JW)
        WRITE (INI,'(1X,I3,1000A10)') K,(CONV(K,I),I=1,IEPR(JW))
      END DO
    END IF
  END DO
END SUBROUTINE PRINT_INITIAL_GRID

!***********************************************************************************************************************************
!*                                  S U B R O U T I N E    D E R I V E D  C O N S T I T U E N T S                                 **
!***********************************************************************************************************************************

SUBROUTINE DERIVED_CONSTITUENTS
  USE CONTROL; USE KINETICS; USE GRIDC
  REAL, ALLOCATABLE, DIMENSION(:,:) :: DOM, POM
  ALLOCATE (DOM(KMX,IMX), POM(KMX,IMX))

  IU    =  CUS(BS(JW))
  ID    =  DS(BE(JW))
  TISS  =  0.0
  ATOT  =  0.0
  TISS  =  0.0
  TOTSS =  0.0
  CHLA  =  0.0
  NBODtot  =  0.0
  PBODtot  =  0.0
  PALT  = (1.0-((EL(KTWB(JW),(IU+ID)/2)-Z((IU+ID)/2)*COS(ALPHA(BE(JW))))/1000.0)/44.3)**5.25
  DO I=IU,ID
    DO K=KTWB(JW),KB(I)
      DO JS=1,NASS
        TISS(K,I) = TISS(K,I)+SS(K,I,SSAC(JS))
      END DO
      DO N=1,NBOD
        NBODtot(K,I)=NBODtot(K,I)+CBODN(K,I,N)
        PBODtot(K,I)=PBODtot(K,I)+CBODP(K,I,N)
!        CBOD(K,I)=CBOD(K,I)+BOD(K,I,N) ! *** NEED TO FIX***
      ENDDO
    END DO
    DO K=KTWB(JW),KB(I)
      DO JA=1,NAAL
        ATOT(K,I) = ATOT(K,I)+ALG(K,I,ALAC(JA))
      END DO
      DOM(K,I) = LDOM(K,I)+RDOM(K,I)
      POM(K,I) = LPOM(K,I)+RPOM(K,I)
      DOC(K,I) = DOM(K,I)*ORGC(JW)
      POC(K,I) = POM(K,I)*ORGC(JW)
      DO JA=1,NAAL
        POC(K,I) = POC(K,I)+ALG(K,I,ALAC(JA))*ALGC(ALAC(JA))
      END DO
      TOC(K,I)   = DOC(K,I)+POC(K,I)   ! ***NEED TO ADD CBOD C TO THIS***
      DOP(K,I)   = LDOMP(K,I)+RDOMP(K,I)    !DOM(K,I)*ORGP(JW)   SW 2/18/16
      DON(K,I)   = LDOMN(K,I)+RDOMN(K,I)    !DOM(K,I)*ORGN(JW)
      POP(K,I)   = LPOMP(K,I)+RPOMP(K,I)    !POM(K,I)*ORGP(JW)
      PON(K,I)   = LPOMN(K,I)+RPOMN(K,I)    !POM(K,I)*ORGN(JW)
      TOP(K,I)   = DOP(K,I)+POP(K,I)+PBODtot(K,I)
      TON(K,I)   = DON(K,I)+PON(K,I)+NBODtot(K,I)
      TKN(K,I)   = TON(K,I)+NH4(K,I)
      CBODU(K,I) = O2OM(JW)*(DOM(K,I)+POM(K,I)+ATOT(K,I))
      TPSS       = 0.0
      DO JS=1,NASS
        TPSS = TPSS+SS(K,I,SSAC(JS))*PARTP(JW)
      END DO
      TPALG=0.0
      TNALG=0.0
      DO JA=1,NAAL                                ! SW 2/18/16
          TPALG=TPALG+ALG(K,I,JA)*ALGP(JA)
          TNALG=TNALG+ALG(K,I,JA)*ALGN(JA)
      ENDDO
      
      TP(K,I)   =  TOP(K,I)+PO4(K,I)+TPSS+TPALG
      TN(K,I)   =  TON(K,I)+NH4(K,I)+NO3(K,I)+TNALG
      SATO      =  EXP(7.7117-1.31403*(LOG(T2(K,I)+45.93)))*PALT
      O2DG(K,I) = (O2(K,I)/SATO)*100.0
      DO JA=1,NAAL
        TOTSS(K,I) = TOTSS(K,I)+ALG(K,I,JA)
        IF (ACHLA(ALAC(JA)) /= 0.0) CHLA(K,I)  = CHLA(K,I)+ALG(K,I,ALAC(JA))/ACHLA(ALAC(JA))
      END DO
      TOTSS(K,I) = TOTSS(K,I)+TISS(K,I)+POM(K,I)
    END DO
  END DO
END SUBROUTINE DERIVED_CONSTITUENTS

!***********************************************************************************************************************************
!*                                            S U B R O U T I N E    W A R N I N G S                                              **
!***********************************************************************************************************************************

SUBROUTINE WARNINGS
  USE CONTROL
  NWRN       =  NWRN+1
  DELETE_WRN = .FALSE.
END SUBROUTINE WARNINGS

!***********************************************************************************************************************************
!*                                              S U B R O U T I N E    E R R O R S                                                **
!***********************************************************************************************************************************

SUBROUTINE ERRORS
  USE CONTROL
  NERR       =  NERR+1
  DELETE_ERR = .FALSE.
END SUBROUTINE ERRORS

!***********************************************************************************************************************************
!*                                           S U B R O U T I N E    I N P U T  E R R O R                                          **
!***********************************************************************************************************************************

SUBROUTINE INPUT_ERROR (UNIT,ERR,AID,MESSAGE)
  INTEGER       :: UNIT,    ERR,      WIN=6
  CHARACTER(8)  :: AID
  CHARACTER(70) :: MESSAGE, ERROR(5)
  ERROR = ' '
  DO J=1,5
    BACKSPACE (UNIT)
  END DO
  DO J=1,5
    READ (UNIT,'(A70)',ERR=1,END=1) ERROR(J)
  END DO
1 CONTINUE
  WRITE (ERR,'(A/(A))') MESSAGE,ERROR
  WRITE (WIN,'(A/(A))') MESSAGE,ERROR
  WRITE (ERR,'(/A/)')  'Last value of card header was '//AID
  WRITE (WIN,'(/A/)')  'Last value of card header was '//AID
  WRITE(WIN,*) 'SEE PRE.ERR'
!  STOP
 END SUBROUTINE INPUT_ERROR

 SUBROUTINE INPUT_ERROR2 (UNIT,ERR,AID,AID1,MESSAGE)
  INTEGER       :: UNIT,    ERR,      WIN=6
  CHARACTER(8)  :: AID,AID1
  CHARACTER(90) :: MESSAGE, ERROR(7)
  ERROR = ' '
  DO J=1,7
    BACKSPACE (UNIT)
  END DO
  DO J=1,7
    READ (UNIT,'(A90)',ERR=1,END=1) ERROR(J)
  END DO
1 CONTINUE
  WRITE (ERR,'(A/(A))') MESSAGE,ERROR
  WRITE (WIN,'(A/(A))') MESSAGE,ERROR
  WRITE (ERR,'(/A/)')  'Read value of header was '//AID//' but the header should be '//AID1
  WRITE (WIN,'(/A/)')  'Read value of header was '//AID//' but the header should be '//AID1
  WRITE(WIN,*) 'SEE PRE.ERR'
 END SUBROUTINE INPUT_ERROR2