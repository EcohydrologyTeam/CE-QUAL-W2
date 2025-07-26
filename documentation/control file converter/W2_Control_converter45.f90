!c control_file_read_write.f90

!***********************************************************************************************************************************
!**                                                      Module Declaration                                                       **
!***********************************************************************************************************************************
MODULE PREC
  INTEGER, PARAMETER :: I2=SELECTED_INT_KIND (3)
  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(15)
END MODULE PREC
MODULE RSTART
  USE PREC
  REAL                                               :: DLTS,   CURMAX, DLTFF, DLTMAXX    ! SW 7/13/2010
  INTEGER                                            :: RSODP,  DLTDP,  TSRDP,  WDODP,  CUF,    RSO=31
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: SNPDP,  VPLDP,  CPLDP,  PRFDP,  SCRDP,  SPRDP,  FLXDP, NSPRF
  REAL                                               :: NXTMRS, NXTMWD, NXTMTS
  REAL,              ALLOCATABLE, DIMENSION(:)       :: NXTMSN, NXTMPR, NXTMSP, NXTMCP, NXTMVP, NXTMSC, NXTMFL
  REAL(R8),          ALLOCATABLE, DIMENSION(:)       :: SBKT,   ELTMF
  REAL(R8),          ALLOCATABLE, DIMENSION(:,:)     :: TSSUH2, TSSDH2, SAVH2,  SAVHR,  SU,     SW,     SAZ
  REAL(R8),          ALLOCATABLE, DIMENSION(:,:,:)   :: CSSUH2, CSSDH2
  REAL(R8)                                           :: ELTM
  REAL(R8),          ALLOCATABLE, DIMENSION(:)       :: VOLIN,  VOLOUT, VOLUH,  VOLDH,  VOLPR,  VOLTRB, VOLDT, VOLWD,  VOLEV
  REAL(R8),          ALLOCATABLE, DIMENSION(:)       :: VOLSBR, VOLTBR, VOLSR,  VOLTR
  REAL(R8),          ALLOCATABLE, DIMENSION(:)       :: TSSEV,  TSSPR,  TSSTR,  TSSDT,  TSSWD,  TSSUH,  TSSDH, TSSIN,  TSSOUT
  REAL(R8),          ALLOCATABLE, DIMENSION(:)       :: TSSS,   TSSB,   TSSICE
  REAL(R8),          ALLOCATABLE, DIMENSION(:)       :: ESBR,   ETBR,   EBRI,   SZ
  REAL(R8),          ALLOCATABLE, DIMENSION(:,:)     :: CMBRT
END MODULE RSTART
MODULE GLOBAL
  USE PREC
  REAL,   PARAMETER                                  :: DAY=86400.0,  NONZERO=1.0E-20, REFL=0.94, FRAZDZ=0.14, DZMIN=1.4E-7
  REAL,   PARAMETER                                  :: AZMIN=1.4E-6, DZMAX=1.0E3,     RHOW=1000.0
  REAL                                               :: DLT,    DLTMIN, DLTTVD, XX
  REAL                                               :: BETABR, START,  HMAX2,  CURRENT
  REAL(R8),   POINTER,            DIMENSION(:,:)     :: U,      W,      T2,     AZ,     RHO,    ST,     SB
  REAL(R8),   POINTER,            DIMENSION(:,:)     :: DLTLIM, VSH,    ADMX,   DM,     ADMZ,   HDG,    HPG,    GRAV
  REAL(R8),   TARGET,ALLOCATABLE, DIMENSION(:,:)     :: T1,     TSS
  REAL(R8),   TARGET,ALLOCATABLE, DIMENSION(:,:,:)   :: C1,     C2,     C1S,    CSSB,   CSSK
  REAL,       TARGET,ALLOCATABLE, DIMENSION(:,:,:)   :: KF,     CD
  REAL(R8),   TARGET,ALLOCATABLE, DIMENSION(:,:,:)   :: HYD
  REAL,   TARGET,    ALLOCATABLE, DIMENSION(:,:,:,:) :: AF,     EF
  REAL,              ALLOCATABLE, DIMENSION(:)       :: ICETH,  ELKT,   HMULT,  CMULT,  CDMULT, WIND2,  AZMAX,  PALT, Z0
  REAL(R8),          ALLOCATABLE, DIMENSION(:,:)     :: QSS,    VOLUH2, VOLDH2, QUH1,   QDH1,   UXBR,   UYBR,   VOL
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: ALLIM,  APLIM,  ANLIM,  ASLIM,  KFS
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: ELLIM,  EPLIM,  ENLIM,  ESLIM
  INTEGER                                            :: W2ERR,  WRN
  INTEGER                                            :: IMX,    KMX,    NBR,    NTR,    NWD,    NWB,    NCT,    NBOD
  INTEGER                                            :: NST,    NSP,    NGT,    NPI,    NPU,    NWDO,   NIKTSR, NUNIT
  INTEGER                                            :: JW,     JB,     JC,     IU,     ID,     KT,     I,      JJB
  INTEGER                                            :: NOD,    NDC,    NAL,    NSS,    NHY,    NFL,    NEP,    NEPT, NDC2
  INTEGER                                            :: NZP,    nzpt, JZ,     NZOOS,  NZOOE,  nmc,   nmct  ! number of zooplankton groups, CONSTIUENT NUMBER FOR ZOOPLANKTON, START AND END
  INTEGER, POINTER,               DIMENSION(:)       :: SNP,    PRF,    VPL,    CPL,    SPR,    FLX,    FLX2
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: BS,     BE,     US,     CUS,    DS,     JBDN
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: KB,     KTI,    SKTI,   KTWB,   KBMIN,  CDHS
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: UHS,    DHS,    UQB,    DQB
  INTEGER, TARGET,   ALLOCATABLE, DIMENSION(:,:)     :: OPT
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: nbodc, nbodn, nbodp        ! cb 6/6/10
  LOGICAL,           ALLOCATABLE, DIMENSION(:)       :: ICE,    ICE_CALC,LAYERCHANGE
  CHARACTER(10)                                      :: CCTIME
  CHARACTER(12)                                      :: CDATE
  CHARACTER(72)                                      :: RSIFN
  CHARACTER(180)                                     :: moddir                     ! current working directory
  REAL(R8),     SAVE, ALLOCATABLE, DIMENSION(:,:)    :: RATZ,   CURZ1,  CURZ2,  CURZ3    ! SW 5/15/06
  real                                               :: g,pi
  REAL(R8)                                           :: DENSITY
  DATA                                        NDC /27/, NHY /15/, NFL /73/, NDC2 /23/      ! cb 6/6/10
  DATA                                        G /9.81/, PI/3.14159265359/
  DATA                                        WRN /32/, W2ERR /33/
  EXTERNAL DENSITY
END MODULE GLOBAL
MODULE GEOMC
  USE PREC
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: JBUH,   JBDH,   JWUH,   JWDH
  REAL(R8),          ALLOCATABLE, DIMENSION(:)       :: ALPHA,  SINA,   COSA,   SLOPE,  BKT,    DLX,    DLXR, SLOPEC, SINAC
  REAL(R8),          ALLOCATABLE, DIMENSION(:,:)     :: H,      H1,     H2,     BH1,    BH2,    BHR1,    BHR2,   AVHR
  REAL(R8),          ALLOCATABLE, DIMENSION(:,:)     :: B,      BI,     BB,     BH,     BHR,    BR,      EL,     AVH1,  AVH2, bnew ! SW 1/23/06
  REAL(R8),          ALLOCATABLE, DIMENSION(:,:)     :: DEPTHB, DEPTHM, FETCHU, FETCHD
  REAL(R8),          ALLOCATABLE, DIMENSION(:)       :: Z, ELWS
END MODULE GEOMC
MODULE NAMESC
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: LNAME
  CHARACTER(6),      ALLOCATABLE, DIMENSION(:)       :: CUNIT,  CUNIT2
  CHARACTER(8),      ALLOCATABLE, DIMENSION(:)       :: CNAME2, CDNAME2
  CHARACTER(9),      ALLOCATABLE, DIMENSION(:)       :: FMTH,   FMTC,   FMTCD
  CHARACTER(19),     ALLOCATABLE, DIMENSION(:)       :: CNAME1
  CHARACTER(43),     ALLOCATABLE, DIMENSION(:)       :: CNAME,  CNAME3, CDNAME, CDNAME3, HNAME
  CHARACTER(72),     ALLOCATABLE, DIMENSION(:)       :: TITLE
  CHARACTER(10),     ALLOCATABLE, DIMENSION(:,:)     :: CONV
  CHARACTER(8) :: NAME1(11)
  CHARACTER(43):: NAME2(11)
  DATA NAME1 /'WaterAge','Bacteria','DGP','N2','H2S','CH4','SO4','FEII','FEOOH','MnII','MnO2'/
  DATA NAME2 /'Age,  days','Bacteria,  col/100ml','Dissolved Gas Pressure,atm','N2 dissolved gas, mg/l','H2S, dissolved gas, mg/l','CH4 dissolved gas, mg/l','SO4 dissolved, mg/l','Reduced FE(II), mg/l','Oxidized FeOOH, mg/l','Reduced Mn(II), mg/l','Oxidized MnO2, mg/l'/


END MODULE NAMESC
MODULE STRUCTURES
  REAL                                               :: DIA,    FMAN,   CLEN,   CLOSS,  UPIE,   DNIE
  REAL,              ALLOCATABLE, DIMENSION(:)       :: QOLD,   QOLDS,  VMAX,   DTP,    DTPS
  REAL,              ALLOCATABLE, DIMENSION(:)       :: EGT,    A1GT,   B1GT,   G1GT,   A2GT,   B2GT,   G2GT, EGT2
  REAL,              ALLOCATABLE, DIMENSION(:)       :: QGT,    GTA1,   GTB1,   GTA2,   GTB2,   BGT
  REAL,              ALLOCATABLE, DIMENSION(:)       :: QSP,    A1SP,   B1SP,   A2SP,   B2SP,   ESP
  REAL,              ALLOCATABLE, DIMENSION(:)       :: EUPI,   EDPI,   WPI,    DLXPI,  FPI,    FMINPI, QPI, BP
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: YS,     VS,     YSS,    VSS,    YST,    VST,    YSTS,   VSTS
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: IUPI,   IDPI,   JWUPI,  JWDPI,  JBDPI,  JBUPI
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: IUSP,   IDSP,   JWUSP,  JWDSP,  JBUSP,  JBDSP
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: IUGT,   IDGT,   JWUGT,  JWDGT,  JBUGT,  JBDGT
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: IWR,    KTWR,   KBWR
  LOGICAL,           ALLOCATABLE, DIMENSION(:)       :: LATERAL_SPILLWAY, LATERAL_PIPE, LATERAL_GATE, LATERAL_PUMP, BEGIN, WLFLAG
  CHARACTER(8),      ALLOCATABLE, DIMENSION(:)       :: LATGTC, LATSPC, LATPIC, LATPUC, DYNGTC, DYNPIPE, DYNPUMP                         ! SW 5/10/10
  CHARACTER(8)                                       :: GT2CHAR 
  REAL,          ALLOCATABLE, DIMENSION(:)     :: EPU,    STRTPU, ENDPU,  EONPU,  EOFFPU, QPU
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: IUPU,   IDPU,   KTPU,   KBPU,   JWUPU,  JWDPU,  JBUPU,  JBDPU
  real :: THR, OMEGA, EPS2
  integer :: NN, NNPIPE, NC
  REAL,              ALLOCATABLE, DIMENSION(:)       :: EGTo,bgto       ! cb/8/13/ 2010
  CHARACTER(8),      ALLOCATABLE, DIMENSION(:)       :: gtic            ! cb/8/13/ 2010
  DATA                                             THR/0.01/, OMEGA/0.8/, EPS2/0.0001/
  DATA                                          NN/19/ ,   NNPIPE /19/, NC/7/
END MODULE STRUCTURES
MODULE TRANS
  USE PREC
  REAL,              ALLOCATABLE, DIMENSION(:)       :: THETA
  REAL(R8),POINTER,               DIMENSION(:,:)     :: COLD,   CNEW,   SSB,    SSK
  REAL(R8),          ALLOCATABLE, DIMENSION(:,:)     :: DX,     DZ,     DZQ
  REAL(R8),          ALLOCATABLE, DIMENSION(:,:)     :: ADX,    ADZ,    AT,     VT,     CT,     DT
END MODULE TRANS
MODULE SURFHE
  REAL                                               :: RHOWCP, PHISET
  REAL,              ALLOCATABLE, DIMENSION(:)       :: ET,     CSHE,   LAT,    LONGIT, SHADE,  RB,     RE,     RC
  REAL,              ALLOCATABLE, DIMENSION(:)       :: WIND,   WINDH,  WSC,    AFW,    BFW,    CFW,    PHI0
  LOGICAL,           ALLOCATABLE, DIMENSION(:)       :: RH_EVAP
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: IWIND  !MLM 08/12/05
END MODULE SURFHE
MODULE TVDC
  REAL,              ALLOCATABLE, DIMENSION(:)       :: QIN,    QTR,    QDTR,   PR,     ELUH,   ELDH,   QWD,    QSUM
  REAL,              ALLOCATABLE, DIMENSION(:)       :: TIN,    TTR,    TDTR,   TPR,    TOUT,   TWDO,   TIND,   QIND
  REAL,              ALLOCATABLE, DIMENSION(:)       :: TAIR,   TDEW,   CLOUD,  PHI,    SRON
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: CUH,    CDH
  INTEGER                                            :: NAC,    NOPEN
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: NACPR,  NACIN,  NACDT,  NACTR,  NACD,   CN
  INTEGER,           ALLOCATABLE, DIMENSION(:,:)     :: TRCN,   INCN,   DTCN,   PRCN
  LOGICAL                                            :: CONSTITUENTS
  CHARACTER(72)                                      :: QGTFN,  QWDFN,  WSCFN,  SHDFN
  CHARACTER(72),     ALLOCATABLE, DIMENSION(:)       :: METFN,  QOTFN,  QINFN,  TINFN,  CINFN,  QTRFN,  TTRFN,  CTRFN,  QDTFN
  CHARACTER(72),     ALLOCATABLE, DIMENSION(:)       :: TDTFN,  CDTFN,  PREFN,  TPRFN,  CPRFN,  EUHFN,  TUHFN,  CUHFN,  EDHFN
  CHARACTER(72),     ALLOCATABLE, DIMENSION(:)       :: EXTFN,  CDHFN,  TDHFN
END MODULE TVDC
MODULE KINETIC
  USE PREC
  REAL                                               :: kdo, CO2PPM                        
  REAL(R8),    POINTER,               DIMENSION(:,:)     :: TDS,    COL,    NH4,    NO3,    PO4,    FE,     DSI,    PSI,    LDOM
  REAL(R8),    POINTER,               DIMENSION(:,:)     :: RDOM,   LPOM,   RPOM,   O2,     TIC,    ALK
  REAL(R8),    POINTER,               DIMENSION(:,:)     :: COLSS,  NH4SS,  NO3SS,  PO4SS,  FESS,   DSISS,  PSISS,  LDOMSS
  REAL(R8),    POINTER,               DIMENSION(:,:)     :: RDOMSS, LPOMSS, RPOMSS, DOSS,   TICSS,  CASS
  REAL,    POINTER,               DIMENSION(:,:)     :: PH,     CO2,    HCO3,   CO3
  REAL,    POINTER,               DIMENSION(:,:)     :: TN,     TP,     TKN
  REAL,    POINTER,               DIMENSION(:,:)     :: DON,    DOP,    DOC
  REAL,    POINTER,               DIMENSION(:,:)     :: PON,    POP,    POC
  REAL,    POINTER,               DIMENSION(:,:)     :: TON,    TOP,    TOC
  REAL,    POINTER,               DIMENSION(:,:)     :: APR,    CHLA,   ATOT
  REAL,    POINTER,               DIMENSION(:,:)     :: O2DG
  REAL,    POINTER,               DIMENSION(:,:)     :: SSSI,   SSSO,   TISS,   TOTSS
  REAL,    POINTER,               DIMENSION(:,:)     :: PO4AR,  PO4AG,  PO4AP,  PO4SD,  PO4SR,  PO4NS,  PO4POM, PO4DOM, PO4OM
  REAL,    POINTER,               DIMENSION(:,:)     :: PO4ER,  PO4EG,  PO4EP,  TICEP,  DOEP,   DOER
  REAL,    POINTER,               DIMENSION(:,:)     :: NH4ER,  NH4EG,  NH4EP,  NO3EG,  DSIEG,  LDOMEP, LPOMEP
  REAL,    POINTER,               DIMENSION(:,:)     :: NH4AR,  NH4AG,  NH4AP,  NH4SD,  NH4SR,  NH4D,   NH4POM, NH4DOM, NH4OM
  REAL,    POINTER,               DIMENSION(:,:)     :: NO3AG,  NO3D,   NO3SED
  REAL,    POINTER,               DIMENSION(:,:)     :: DSIAG,  DSID,   DSISD,  DSISR,  DSIS
  REAL,    POINTER,               DIMENSION(:,:)     :: PSIAM,  PSID,   PSINS
  REAL,    POINTER,               DIMENSION(:,:)     :: FENS,   FESR
  REAL,    POINTER,               DIMENSION(:,:)     :: LDOMAP, LDOMD,  LRDOMD, RDOMD
  REAL,    POINTER,               DIMENSION(:,:)     :: LPOMAP, LPOMD,  LRPOMD, RPOMD,  LPOMNS, RPOMNS
  REAL,    POINTER,               DIMENSION(:,:)     :: DOAP,   DOAR,   DODOM,  DOPOM,  DOOM,   DONIT
  REAL,    POINTER,               DIMENSION(:,:)     :: DOSED,  DOSOD,  DOBOD,  DOAE
  REAL,    POINTER,               DIMENSION(:,:)     :: CBODU,  CBODDK, TICAP
  REAL,    POINTER,               DIMENSION(:,:)     :: SEDD,   SODD,   SEDAS,  SEDOMS, SEDNS
  REAL(R8),POINTER,               DIMENSION(:,:,:)   :: SS,     ALG,    CBOD,   CG
  REAL(R8),POINTER,               DIMENSION(:,:,:)   :: SSSS,   ASS,    CBODSS, CGSS
  REAL,    POINTER,               DIMENSION(:,:,:)   :: AGR,    ARR,    AER,    AMR,    ASR
  REAL,    POINTER,               DIMENSION(:,:,:)   :: EGR,    ERR,    EER,    EMR,    EBR
  REAL(R8),POINTER,               DIMENSION(:,:)     :: LDOMP,  RDOMP,  LPOMP,  RPOMP,  LDOMN,  RDOMN,  LPOMN,  RPOMN
  REAL(R8),POINTER,               DIMENSION(:,:)     :: LDOMPSS,  RDOMPSS, LPOMPSS, RPOMPSS, LDOMNSS, RDOMNSS
  REAL(R8),POINTER,               DIMENSION(:,:)     :: LPOMNSS,  RPOMNSS
  REAL,    POINTER,               DIMENSION(:,:)     :: LDOMPAP,  LDOMPEP, LPOMPAP, LPOMPNS, RPOMPNS
  REAL,    POINTER,               DIMENSION(:,:)     :: LDOMNAP,  LDOMNEP, LPOMNAP, LPOMNNS, RPOMNNS
  REAL,    POINTER,               DIMENSION(:,:)     :: SEDDP,    SEDASP,  SEDOMSP, SEDNSP,  LPOMEPP
  REAL,    POINTER,               DIMENSION(:,:)     :: SEDDN,    SEDASN,  SEDOMSN, SEDNSN,  LPOMEPN, SEDNO3
  REAL,    POINTER,               DIMENSION(:,:)     :: SEDDC,    SEDASC,  SEDOMSC, SEDNSC,  LPOMEPC
  REAL,    POINTER,               DIMENSION(:,:)     :: CBODNS,   SEDCB,   SEDCBP,  SEDCBN,  SEDCBC
  REAL,    POINTER,               DIMENSION(:,:)     :: sedbr,    sedbrp,  sedbrc,  sedbrn, co2reaer        !cb 11/30/06
  REAL(R8),POINTER,               DIMENSION(:,:,:)   :: cbodp,    cbodn       ! cb 6/6/10
  REAL(R8),POINTER,               DIMENSION(:,:,:)   :: cbodpss,    cbodnss       ! cb 6/6/10
  REAL,    POINTER,               DIMENSION(:,:)     :: CBODNSp,  CBODNSn          ! cb 6/6/10
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: EPM,    EPD,    EPC
  REAL,              ALLOCATABLE, DIMENSION(:)       :: CGQ10,  CG0DK,  CG1DK,  CGS, CGLDK, CGKLF, CGCS  !LCJ 2/26/15 SW 10/16/15
  REAL,              ALLOCATABLE, DIMENSION(:)       :: SOD,    SDK,    LPOMDK, RPOMDK, LDOMDK, RDOMDK, LRDDK,  LRPDK
  REAL,              ALLOCATABLE, DIMENSION(:)       :: SSS,    TAUCR,  POMS,   FES, seds, sedb  !cb 11/27/06
  REAL,              ALLOCATABLE, DIMENSION(:)       :: AG,     AR,     AE,     AM,     AS,     AHSN,   AHSP,   AHSSI,  ASAT
  REAL,              ALLOCATABLE, DIMENSION(:)       :: AP,     AN,     AC,     ASI,    ACHLA,  APOM,   ANPR
  REAL,              ALLOCATABLE, DIMENSION(:)       :: EG,     ER,     EE,     EM,     EB
  REAL,              ALLOCATABLE, DIMENSION(:)       :: EHSN,   EHSP,   EHSSI,  ESAT,   EHS,    ENPR
  REAL,              ALLOCATABLE, DIMENSION(:)       :: EP,     EN,     EC,     ESI,    ECHLA,  EPOM
  REAL,              ALLOCATABLE, DIMENSION(:)       :: BETA,   EXH2O,  EXSS,   EXOM,   EXA
  REAL,              ALLOCATABLE, DIMENSION(:)       :: DSIR,   PSIS,   PSIDK,  PARTSI
  REAL,              ALLOCATABLE, DIMENSION(:)       :: ORGP,   ORGN,   ORGC,   ORGSI
  REAL,              ALLOCATABLE, DIMENSION(:)       :: BODP,   BODN,   BODC
  REAL,              ALLOCATABLE, DIMENSION(:)       :: PO4R,   PARTP
  REAL,              ALLOCATABLE, DIMENSION(:)       :: NH4DK,  NH4R,   NO3DK,  NO3S, FNO3SED
  REAL,              ALLOCATABLE, DIMENSION(:)       :: O2AG,   O2AR,   O2OM,   O2NH4
  REAL,              ALLOCATABLE, DIMENSION(:)       :: O2EG,   O2ER
  REAL,              ALLOCATABLE, DIMENSION(:)       :: CO2R,   FER
  REAL,              ALLOCATABLE, DIMENSION(:)       :: KBOD,   TBOD,   RBOD
  REAL,              ALLOCATABLE, DIMENSION(:)       :: CAQ10,  CADK,   CAS
  REAL,              ALLOCATABLE, DIMENSION(:)       :: OMT1,   OMT2,   SODT1,  SODT2,  NH4T1,  NH4T2,  NO3T1,  NO3T2
  REAL,              ALLOCATABLE, DIMENSION(:)       :: OMK1,   OMK2,   SODK1,  SODK2,  NH4K1,  NH4K2,  NO3K1,  NO3K2
  REAL,              ALLOCATABLE, DIMENSION(:)       :: AT1,    AT2,    AT3,    AT4
  REAL,              ALLOCATABLE, DIMENSION(:)       :: AK1,    AK2,    AK3,    AK4
  REAL,              ALLOCATABLE, DIMENSION(:)       :: ET1,    ET2,    ET3,    ET4
  REAL,              ALLOCATABLE, DIMENSION(:)       :: EK1,    EK2,    EK3,    EK4
  REAL,              ALLOCATABLE, DIMENSION(:)       :: REAER,  WIND10, CZ,     QC,     QERR
  REAL,              ALLOCATABLE, DIMENSION(:)       :: RCOEF1, RCOEF2, RCOEF3, RCOEF4
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: DO1,    DO2,    DO3,    GAMMA
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: SED,    FPSS,   FPFE
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: CBODD
  REAL,              ALLOCATABLE, DIMENSION(:)       :: CBODS
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: ORGPLD,  ORGPRD,   ORGPLP,    ORGPRP,  ORGNLD,  ORGNRD, ORGNLP, ORGNRP
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: LDOMPMP, LDOMNMP,  LPOMPMP,   LPOMNMP, RPOMPMP, RPOMNMP
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: LPZOOINP,LPZOOINN, LPZOOOUTP, LPZOOOUTN
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: SEDC,    SEDN, SEDP
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: SEDVPC,  SEDVPP, SEDVPN
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: SDKV,    SEDDKTOT
  INTEGER                                            :: nldomp,nrdomp,nlpomp,nrpomp,nldomn,nrdomn,nlpomn,nrpomn
  INTEGER,           ALLOCATABLE, DIMENSION(:)       :: NAF,    NEQN,   ANEQN,  ENEQN
  INTEGER,           ALLOCATABLE, DIMENSION(:,:)     :: KFCN
  LOGICAL,           ALLOCATABLE, DIMENSION(:)       :: SEDIMENT_RESUSPENSION
  CHARACTER(8),      ALLOCATABLE, DIMENSION(:)       :: CAC,    REAERC
  CHARACTER(10),     ALLOCATABLE, DIMENSION(:,:)     :: LFPR
  CHARACTER(8)                                       :: CO2YRLY
END MODULE KINETIC
MODULE SELWC
  REAL,              ALLOCATABLE, DIMENSION(:)   :: EWD,    VNORM,  QNEW, tavgw
  REAL,              ALLOCATABLE, DIMENSION(:,:) :: QSTR,   QSW,    ESTR,   WSTR, TAVG            ! SW Selective 7/30/09
  REAL,              ALLOCATABLE, DIMENSION(:,:) :: CAVGW, CDAVGW
  REAL,              ALLOCATABLE, DIMENSION(:,:,:):: CAVG, CDAVG
  INTEGER,           ALLOCATABLE, DIMENSION(:)   :: NSTR,   NOUT,   KTWD,   KBWD,   KTW,   KBW
  INTEGER,           ALLOCATABLE, DIMENSION(:,:) :: KTSW,   KBSW,   KOUT
  CHARACTER(8),      ALLOCATABLE, DIMENSION(:) :: DYNSTRUC
END MODULE SELWC
MODULE GDAYC
  REAL                                           :: DAYM,   EQTNEW
  INTEGER                                        :: JDAYG,  M,      YEAR,   GDAY
  LOGICAL                                        :: LEAP_YEAR
  CHARACTER(9)                                   :: MONTH
END MODULE GDAYC
MODULE SCREENC
  USE PREC
  REAL                                           :: JDAY,   DLTS1,  JDMIN,  MINDLT, DLTAV,  ELTMJD
  REAL(R8),          ALLOCATABLE, DIMENSION(:)   :: ZMIN,   CMIN,   CMAX,   HYMIN,  HYMAX,  CDMIN,  CDMAX
  INTEGER                                        :: ILOC,   KLOC,   IMIN,   KMIN,   NIT,    NV,     JTT,     JWW
  INTEGER,           ALLOCATABLE, DIMENSION(:)   :: IZMIN
  CHARACTER(8),      ALLOCATABLE, DIMENSION(:)   :: ACPRC,  AHPRC,  ACDPRC
END MODULE SCREENC
MODULE TDGAS
  REAL,              ALLOCATABLE, DIMENSION(:)   :: AGASSP, BGASSP, CGASSP, AGASGT, BGASGT, CGASGT
  INTEGER,           ALLOCATABLE, DIMENSION(:)   :: EQSP,   EQGT
END MODULE TDGAS
MODULE LOGICC
  LOGICAL                                        :: SUSP_SOLIDS,        OXYGEN_DEMAND,    UPDATE_GRAPH,     INITIALIZE_GRAPH
  LOGICAL                                        :: WITHDRAWALS,        TRIBUTARIES,      GATES, PIPES
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: NO_WIND,            NO_INFLOW,        NO_OUTFLOW,       NO_HEAT
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: UPWIND,             ULTIMATE,         FRESH_WATER,      SALT_WATER
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: LIMITING_DLT,       TERM_BY_TERM,     MANNINGS_N,       PH_CALC
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: ONE_LAYER,          DIST_TRIBS,       PRECIPITATION
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: PRINT_SEDIMENT,     LIMITING_FACTOR,  READ_EXTINCTION,  READ_RADIATION
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: UH_INTERNAL,        DH_INTERNAL,      UH_EXTERNAL,      DH_EXTERNAL
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: UQ_INTERNAL,        DQ_INTERNAL,      UQ_EXTERNAL,      DQ_EXTERNAL
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: UP_FLOW,            DN_FLOW,          INTERNAL_FLOW
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: DAM_INFLOW,         DAM_OUTFLOW                                    !TC 08/03/04
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: INTERP_METEOROLOGY, INTERP_INFLOW,    INTERP_DTRIBS,    INTERP_TRIBS
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: INTERP_WITHDRAWAL,  INTERP_HEAD,      INTERP_EXTINCTION
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: VISCOSITY_LIMIT,    CELERITY_LIMIT,   IMPLICIT_AZ,      TRAPEZOIDAL !SW 07/16/04
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: HYDRO_PLOT,         CONSTITUENT_PLOT, DERIVED_PLOT
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: INTERP_gate     ! cb 8/13/2010
  LOGICAL,           ALLOCATABLE, DIMENSION(:,:) :: PRINT_DERIVED,      PRINT_HYDRO,      PRINT_CONST,      PRINT_EPIPHYTON
  LOGICAL,           ALLOCATABLE, DIMENSION(:,:) :: POINT_SINK,         INTERNAL_WEIR,    INTERP_OUTFLOW
END MODULE LOGICC
MODULE SHADEC
  integer, PARAMETER :: IANG=18
  REAL,PARAMETER                                 :: GAMA=(3.1415926*2.)/REAL(IANG)                         ! SW 10/17/05
  REAL,                           DIMENSION(IANG):: ANG                                                    ! SW 10/17/05
  REAL,              ALLOCATABLE, DIMENSION(:)   :: A00,    DECL,   HH,     TTLB,   TTRB,   CLLB,   CLRB   ! SW 10/17/05
  REAL,              ALLOCATABLE, DIMENSION(:)   :: SRLB1,  SRRB1,  SRLB2,  SRRB2,  SRFJD1, SRFJD2, SHADEI
  REAL,              ALLOCATABLE, DIMENSION(:,:) :: TOPO
  LOGICAL,           ALLOCATABLE, DIMENSION(:)   :: DYNAMIC_SHADE
  DATA ANG  /0.00000, 0.34907, 0.69813, 1.04720, 1.39626, 1.74533, 2.09440, 2.44346, &
            2.79253, 3.14159, 3.49066, 3.83972, 4.18879, 4.53786, 4.88692, 5.23599, 5.58505, 5.93412/      ! SW 10/17/05
END MODULE SHADEC
MODULE EDDY
USE PREC
  CHARACTER(8),      ALLOCATABLE, DIMENSION(:)      :: AZC,IMPTKE
  REAL,              ALLOCATABLE, DIMENSION(:)      :: WSHY,   FRIC
  REAL,              ALLOCATABLE, DIMENSION(:,:)    :: FRICBR, DECAY
  REAL(R8),          ALLOCATABLE, DIMENSION (:,:,:) :: TKE
  REAL(R8),          ALLOCATABLE, DIMENSION (:,:)   :: AZT, DZT
  REAL,              ALLOCATABLE, DIMENSION(:)      :: USTARBTKE, E
  REAL,              ALLOCATABLE, DIMENSION(:)      :: EROUGH, ARODI, TKELATPRDCONST, STRICK
  INTEGER,           ALLOCATABLE, DIMENSION(:)      :: FIRSTI, LASTI, WALLPNT, TKEBC
  LOGICAL,           ALLOCATABLE, DIMENSION(:)      :: STRICKON, TKELATPRD
END MODULE EDDY
MODULE MACROPHYTEC
  REAL,    POINTER,               DIMENSION(:,:)     :: NH4MR,  NH4MG,  LDOMMAC, RPOMMAC, LPOMMAC, DOMP, DOMR, TICMC
  REAL,    POINTER,               DIMENSION(:,:)     :: PO4MR,  PO4MG
  REAL,              ALLOCATABLE, DIMENSION(:)       :: MG,     MR,     MM, MMAX,   MBMP
  REAL,              ALLOCATABLE, DIMENSION(:)       :: MT1,    MT2,    MT3,    MT4,    MK1,    MK2,    MK3,    MK4
  REAL,              ALLOCATABLE, DIMENSION(:)       :: MP,     MN,     MC
  REAL,              ALLOCATABLE, DIMENSION(:)       :: PSED,   NSED,   MHSP,   MHSN,   MHSC,   msat,   exm
  REAL,              ALLOCATABLE, DIMENSION(:)       :: CDdrag, dwv,    dwsa,  anorm    ! cb 6/29/06
  REAL,              ALLOCATABLE, DIMENSION(:)       :: ARMAC
  REAL,              ALLOCATABLE, DIMENSION(:)       :: O2MG,   O2MR,   LRPMAC,  MPOM
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: MACMBRS,MACMBRT,SSMACMB
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: CW,     BIC, macwbci
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: MACTRMR,MACTRMF,MACTRM
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: MMR,    MRR
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: MAC,    MACT
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: MPLIM,  MNLIM, MCLIM
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: SMAC,   SMACT
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: GAMMAJ
  REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: MGR
  REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: MACRC,  MACRM
  REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: MLLIM
  REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: MACSS
  REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: SMACRC, SMACRM
  LOGICAL,           ALLOCATABLE, DIMENSION(:)       :: KTICOL
  LOGICAL,           ALLOCATABLE, DIMENSION(:,:)     :: PRINT_MACROPHYTE, MACROPHYTE_CALC
  LOGICAL                                            :: MACROPHYTE_ON
  CHARACTER(3),      ALLOCATABLE, DIMENSION(:,:)     :: mprwbc, macwbc
  CHARACTER(10),      ALLOCATABLE, DIMENSION(:,:)    :: CONV2
  CHARACTER(10),     ALLOCATABLE, DIMENSION(:,:,:,:) :: MLFPR
!  DATA                                                  SAVOLRAT /9000.0/, DEN /6.0E4/   !cb 6/30/06
END MODULE MACROPHYTEC
MODULE POROSITYC
    REAL,              ALLOCATABLE, DIMENSION(:)     :: SAREA, VOLKTI
    REAL,              ALLOCATABLE, DIMENSION(:,:)   :: POR,   VOLI,   VSTEMKT
    REAL,              ALLOCATABLE, DIMENSION(:,:,:) :: VSTEM
    LOGICAL,       ALLOCATABLE, DIMENSION(:)         :: HEAD_FLOW
    LOGICAL,       ALLOCATABLE, DIMENSION(:)         :: UP_HEAD
END MODULE POROSITYC
MODULE ZOOPLANKTONC
  USE PREC
  LOGICAL                                            :: ZOOPLANKTON_CALC
  REAL,              ALLOCATABLE, DIMENSION(:)       :: zg,zm,zeff,PREFP,zr,ZOOMIN,ZS2P,EXZ
  REAL,              ALLOCATABLE, DIMENSION(:)       :: Zt1,Zt2,Zt3,Zt4,Zk1,Zk2,Zk3,Zk4
  REAL,              ALLOCATABLE, DIMENSION(:)       :: ZP,ZN,ZC,o2zr
    REAL,              ALLOCATABLE, DIMENSION(:,:)   :: PREFA, PREFZ ! OMNIVOROUS ZOOPLANKTON
  REAL,              ALLOCATABLE, DIMENSION(:,:)     :: po4zr,NH4ZR,DOZR,TICZR,LPZOOOUT,LPZOOIN
  REAL(R8),POINTER,               DIMENSION(:,:,:)   :: ZOO, ZOOSS
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: ZMU,TGRAZE,ZRT,ZMT
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: ZOORM,ZOORMR,ZOORMF
  REAL,              ALLOCATABLE, DIMENSION(:,:,:)   :: agzt
  REAL,              ALLOCATABLE, DIMENSION(:,:,:,:) :: AGZ, ZGZ ! OMNIVOROUS ZOOPLANKTON
END MODULE ZOOPLANKTONC
module initialvelocity
  LOGICAL                                            :: init_vel, once_through
  REAL,          ALLOCATABLE, DIMENSION(:)           :: qssi,elwss,uavg
  REAL,          ALLOCATABLE, DIMENSION(:,:)         :: bsave
  LOGICAL,           ALLOCATABLE, DIMENSION(:)       :: loop_branch
end module initialvelocity
MODULE ENVIRPMOD
character*3, save, allocatable, dimension (:) :: CC_E,CD_E
character*3, save :: vel_vpr,temp_vpr
real, save :: temp_int,temp_top,vel_int,vel_top,sumvolt,timlast,dltt,v_cnt,v_sum,v_tot,volgl,t_crit,v_crit,c_crit,cd_crit,temp_c,v_avg,vel_c
real, save :: t_cnt,t_sum,t_tot,t_avg
real,allocatable, save, dimension (:) :: c_cnt,cd_cnt,c_tot,cd_tot,t_class,v_class,c_sum,cd_sum,c_avg,cd_avg
real,allocatable, save, dimension (:,:) :: c_class,cd_class,conc_c,conc_cd
real, allocatable, save, dimension (:) :: c_int,c_top,cd_int,cd_top,cn_e,cdn_e
integer CONE,numclass,iopenfish,nac_e,nacd_e,jj,jacd
data CONE/1500/
END MODULE ENVIRPMOD
Module MAIN
USE PREC
! Variable declaration
  CHARACTER(3), ALLOCATABLE, DIMENSION(:) :: DLTADD
  INTEGER       :: J,NIW,NGC,NGCS,NTDS,NCCS,NGCE,NSSS,NSSE,NPO4,NNH4
  INTEGER       :: NNO3,NDSI,NPSI,NFE,NLDOM,NRDOM,NLPOM,NRPOM,NBODS
  INTEGER       :: NBODE, NAS, NAE, NDO, NTIC, NALK, NTRT, NWDT
  INTEGER       :: NDT, JS, JP, JG,JT, JH, NTSR, NIWDO, NRSO,JD
  INTEGER       :: nbodcs, nbodce, nbodps, nbodpe, nbodns, nbodne, ibod, jcb       ! cb 6/6/10
  INTEGER       :: JF,JA,JM,JE,JJZ,K,L3,L1,L2,NTAC,NDSP,NTACMX,NTACMN,JFILE
  INTEGER       :: KBP,JWR,JJJB,JDAYNX,NIT1,JWD,L,IUT,IDT,KTIP
  INTEGER       :: INCRIS,IE,II,NDLT,NRS,INCR,IS,JAC
  REAL          :: JDAYTS, JDAY1, TMSTRT, TMEND,HMAX, DLXMAX,CELRTY
  REAL(R8)      :: DLMR, TICE                        ! SW 4/19/10
  REAL          :: TAU1,TAU2, ELTMS, EPI,HMIN,DLXMIN, RHOICP
  REAL          :: NXTVD,TTIME, ZB,WWT,DFC,GC2,HRAD,EFFRIC,UDR,UDL,AB
  REAL          :: DEPKTI,COLB,COLDEP,SSTOT,RHOIN,VQIN,VQINI
  REAL          :: QINFR, ELT,RHOIRL1,V1,BHSUM,BHRSUM,WT1,WT2
  REAL          :: ICETHU, ICETH1, ICETH2, ICE_TOL,DEL,HICE            ! SW 4/19/10
  REAL          :: DLTCAL,HEATEX,SROOUT,SROSED,SROIN,SRONET,TFLUX,HIA
  REAL          :: TAIRV,EA,ES,DTV
  REAL          :: T2R4
  INTEGER       :: CON,    RSI,    GRF,  NDG=16             ! cb 1/26/09
  integer       :: vsf,    sif 
  LOGICAL       :: ADD_LAYER,      SUB_LAYER,          WARNING_OPEN,    ERROR_OPEN,      VOLUME_WARNING, SURFACE_WARNING
  LOGICAL       :: END_RUN,        BRANCH_FOUND,       NEW_PAGE,        UPDATE_KINETICS, UPDATE_RATES
  LOGICAL       :: WEIR_CALC,      DERIVED_CALC,       RESTART_IN,      RESTART_OUT
  LOGICAL       :: SPILLWAY,       PUMPS
  LOGICAL       :: TIME_SERIES,    DOWNSTREAM_OUTFLOW, ICE_COMPUTATION            !, WINTER    ! SW/RC 4/28/11 eliminate WINTER
  CHARACTER(1)  :: ESC
  CHARACTER(2)  :: DEG
  CHARACTER(3)  :: GDCH
  CHARACTER(8)  :: RSOC,   RSIC,   CCC,   LIMC,   WDOC,   TSRC,   EXT, SELECTC, CLOSEC, HABTATC,ENVIRPC, AERATEC, inituwl, DLTINTER      ! SW 7/31/09; 8/24/09
  CHARACTER(8)  ::SYSTDG,   N2BND,   DOBND, CDUM  
  CHARACTER(10) :: BLANK,  BLANK1, sedch,   sedpch,   sednch,   sedcch 
  CHARACTER(72) :: WDOFN,  RSOFN,  TSRFN, SEGNUM, LINE, SEGNUM2
  LOGICAL       :: RETLOG

! Allocatable array declarations

  REAL,          ALLOCATABLE, DIMENSION(:)     :: ETUGT,  EBUGT,  ETDGT,  EBDGT,DDUM,FDUM
  REAL,          ALLOCATABLE, DIMENSION(:)     :: ETUSP,  EBUSP,  ETDSP,  EBDSP
  REAL,          ALLOCATABLE, DIMENSION(:)     :: ETUPI,  EBUPI,  ETDPI,  EBDPI,  ETPU,   EBPU,   TSEDF
  REAL,          ALLOCATABLE, DIMENSION(:)     :: CSUM,   CDSUM,  X1
  REAL,          ALLOCATABLE, DIMENSION(:)     :: RSOD,   RSOF,   DLTD,   DLTF,   DLTMAX, QWDO
  REAL,          ALLOCATABLE, DIMENSION(:)     :: ICETHI, ALBEDO, HWI,    BETAI,  GAMMAI, ICEMIN, ICET2,  CBHE,   TSED
  REAL,          ALLOCATABLE, DIMENSION(:)     :: FI,     SEDCI,  FSOD,   FSED,   AX,     RANLW,    T2I,    ELBOT,  DXI
  REAL,          ALLOCATABLE, DIMENSION(:)     :: QINT,   QOUTT
  REAL,          ALLOCATABLE, DIMENSION(:)     :: WSHX,   SROSH,  EV
  REAL,          ALLOCATABLE, DIMENSION(:)     :: QDT,    QPR,    ICESW,  RS,     RN
  REAL,          ALLOCATABLE, DIMENSION(:)     :: XBR,    QPRBR,  EVBR,   TPB
  REAL(R8),      ALLOCATABLE, DIMENSION(:)     :: DLXRHO, Q,      QSSUM
  REAL,          ALLOCATABLE, DIMENSION(:)     :: ELTRT,  ELTRB
  REAL,          ALLOCATABLE, DIMENSION(:)     :: TSRD,   TSRF,   WDOD,   WDOF
  REAL,          ALLOCATABLE, DIMENSION(:)     :: QOAVR,  QIMXR,  QOMXR,  QTAVB,  QTMXB
  REAL,          ALLOCATABLE, DIMENSION(:)     :: FETCH,  ETSR
  REAL(R8),      ALLOCATABLE, DIMENSION(:)     :: QINSUM, TINSUM
  REAL,          ALLOCATABLE, DIMENSION(:)     :: CDTOT
  REAL,          ALLOCATABLE, DIMENSION(:)     :: SEDCIp, sedcin, sedcic, sedcis   !v3.5
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: ESTRT,  WSTRT,  CINSUM
  REAL(R8),      ALLOCATABLE, DIMENSION(:,:)   :: P,      HSEG,   QTOT
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: CPB,    COUT,   CWDO,   CDWDO, KFJW
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: C2I,    EPICI
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: QTRF
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: SNPD,   SCRD,   PRFD,   SPRD,   CPLD,   VPLD,   FLXD
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: SNPF,   SCRF,   PRFF,   SPRF,   CPLF,   VPLF,   FLXF
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: TVP,    SEDVP,  QINF
  REAL(R8),      ALLOCATABLE, DIMENSION(:,:)   :: TSSUH1, TSSDH1
  REAL(R8),      ALLOCATABLE, DIMENSION(:,:,:) :: CSSUH1, CSSDH1
  REAL,          ALLOCATABLE, DIMENSION(:,:,:) :: EPIVP,  CVP
  REAL(R8),      ALLOCATABLE, DIMENSION(:)     :: VOLB
  REAL(R8),      ALLOCATABLE, DIMENSION(:)     :: DLVOL,  VOLG
  REAL(R8),      ALLOCATABLE, DIMENSION(:)     :: A,      C,      D,      F,      V,      BTA,    GMA,    BHRHO
  REAL(R8),      ALLOCATABLE, DIMENSION(:)     :: DLVR,   ESR,    ETR
  REAL(R8),      ALLOCATABLE, DIMENSION(:,:)   :: CMBRS
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: KTUGT,  KBUGT,  KTDGT,  KBDGT, IDUM
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: KTUSP,  KBUSP,  KTDSP,  KBDSP
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: KTUPI,  KBUPI,  KTDPI,  KBDPI
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: NSNP,   NSCR,   NSPR,   NVPL,   NFLX,   NCPL,   BTH
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: VPR,    LPR,    NIPRF,  NISPR,  NPRF
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: NISNP
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: NBL,    KBMAX,  KBI
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: KBR,    IBPR
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: TSR
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: NPOINT, NL,     KTQIN,  KBQIN, ilayer    ! SW 1/23/06
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: ITR,    KTTR,   KBTR,   JBTR
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: IWD,    KWD,    JBWD
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: IWDO,   ITSR
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: ILAT,   JBDAM,  JSS
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: icpl                                      ! cb 1/26/09
  INTEGER,       ALLOCATABLE, DIMENSION(:,:)   :: KTSWT,  KBSWT
  INTEGER,       ALLOCATABLE, DIMENSION(:,:)   :: IPRF,   ISPR,   ISNP,   BL,     WDO,    CDN, WDO2
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: ALLOW_ICE,           PUMPON,        FETCH_CALC,   ICE_IN  !     RC/SW 4/28/11
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: DN_HEAD,        HEAD_BOUNDARY
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: PLACE_QIN,      PLACE_QTR,      SPECIFY_QTR
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: OPEN_VPR,       OPEN_LPR
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: ISO_TEMP,       VERT_TEMP,      LONG_TEMP,     VERT_PROFILE,  LONG_PROFILE
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: SEDIMENT_CALC,  DETAILED_ICE,   IMPLICIT_VISC, SNAPSHOT,      PROFILE
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: VECTOR,         CONTOUR,        SPREADSHEET,   SCREEN_OUTPUT
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: FLUX,           EVAPORATION,    ZERO_SLOPE
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: ISO_SEDIMENT,   VERT_SEDIMENT,  LONG_SEDIMENT
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: VOLUME_BALANCE, ENERGY_BALANCE, MASS_BALANCE, BOD_CALC, ALG_CALC
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: bod_calcp,bod_calcn                                                  ! cb 5/19/2011
  LOGICAL,       ALLOCATABLE, DIMENSION(:,:)   :: ISO_EPIPHYTON,  VERT_EPIPHYTON, LONG_EPIPHYTON, EPIPHYTON_CALC
  LOGICAL,       ALLOCATABLE, DIMENSION(:,:)   :: ISO_CONC,       VERT_CONC,      LONG_CONC,      TDG_SPILLWAY,   TDG_GATE
  CHARACTER(4),  ALLOCATABLE, DIMENSION(:)     :: CUNIT1
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: SEG,    SEDRC,  TECPLOT
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: HPLTC,  CPLTC,  CDPLTC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: EXC,    EXIC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: GASGTC, GASSPC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: CWDOC,  CDWDOC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: ICEC,   SEDCc,  SEDPRC, SNPC,   SCRC,   SPRC,   PRFC,DYNSEDK
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: RHEVC,  VPLC,   CPLC,   AZSLC,  FETCHC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: DTRC,   SROC,   KFAC,   CDAC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: INCAC,  TRCAC,  DTCAC,  PRCAC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: WTYPEC, GRIDC                                                        !SW 07/16/04
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: PUSPC,  PDSPC,  PUGTC,  PDGTC,  PDPIC,  PUPIC,  PPUC,   TRC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: SLICEC, FLXC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: VBC,    MBC,    EBC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: PQC,    EVC,    PRC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: QINC,   QOUTC,  WINDC,  HEATC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: VISC,   CELC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: SLTRC,  SLHTC,  FRICC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: QINIC,  TRIC,   DTRIC,  WDIC,   HDIC,   METIC, KFNAME2
  CHARACTER(10), ALLOCATABLE, DIMENSION(:)     :: C2CH,   CDCH,   EPCH,   macch, KFCH
  CHARACTER(45), ALLOCATABLE, DIMENSION(:)     :: KFNAME
  CHARACTER(72), ALLOCATABLE, DIMENSION(:)     :: SNPFN,  PRFFN,  VPLFN,  CPLFN,  SPRFN,  FLXFN,  FLXFN2, BTHFN,  VPRFN,  LPRFN
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:,:)   :: SINKC,  SINKCT
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:,:)   :: CPRBRC, CDTBRC, CPRWBC, CINBRC, CTRTRC, HPRWBC, STRIC,  CDWBC,  KFWBC
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:,:)   :: EPIC,   EPIPRC
  CHARACTER(10), ALLOCATABLE, DIMENSION(:,:)   :: CONV1
  CHARACTER(72)                                :: CONFN='w2_con.npt'
  CHARACTER(72)                                :: TEXT
  integer nproc

! Data declarations
  Real  :: RK1, RL1, RIMT, RHOA, RHOI, VTOL, CP, thrkti
  DATA RK1   /2.12/,         RL1    /333507.0/, RIMT /0.0/, RHOA /1.25/, RHOI /916.0/, VTOL /1.0E3/, CP /4186.0/
  DATA ICE_TOL /0.005/
  DATA BLANK /'          '/, BLANK1 /'    -99.00'/
  DATA CON   /10/,  RSI /11/
  data thrkti /0.10/  

END Module Main

Module Selective1 
 REAL                                          :: nxtstr, nxttcd, nxtsplit,tcdfreq,tfrqtmp
  CHARACTER(8)                                 :: tempc,tspltc
  CHARACTER(8), ALLOCATABLE, DIMENSION(:)      :: tcelevcon,tcyearly,tcntr,tspltcntr,monctr,tsyearly,DYNSEL, ELCONTSPL
  INTEGER                                      :: numtempc,numtsplt, tempn,nstt    
  INTEGER, ALLOCATABLE, DIMENSION(:)           :: tcnelev,tcjb,tcjs,tciseg,tspltjb,nouts,kstrsplt, jbmon, jsmon, ncountcw,SELD
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: tcelev, tempcrit,qstrfrac
  REAL,          ALLOCATABLE, DIMENSION(:)     :: tctemp,tctend,tctsrt,tcklay,tspltt,volm,qwdfrac,tstend,tstsrt,NXSEL,TEMP2
  INTEGER, ALLOCATABLE, DIMENSION(:,:)         :: jstsplt, ncountc
   REAL,          ALLOCATABLE, DIMENSION(:,:)  :: volmc 

End Module Selective1

module habitat
  integer :: ifish, n, nseg,iopenfish,kkmax,kseg,jjw
  real    :: voltot,O2CORR,DOSAT
  character*80, allocatable, dimension(:) :: fishname
  character*80 :: conhab,conavg,consurf,consod
  character*300 :: habline1, habline2,habline3,habline4,habline5,habline6,habline7
  real, allocatable, dimension(:) :: fishtempl,fishtemph,fishdo,habvol,phabvol,cdo,cpo4,cno3,cnh4,cchla,ctotp,cdos,cpo4s,cno3s,cnh4s,cchlas,ctotps,cgamma,ssedd
  integer, allocatable, dimension (:) :: isegvol
end module habitat

MODULE extra
  CHARACTER(1),  ALLOCATABLE, DIMENSION(:)     :: bthtype, vprtype
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: hyd_prin,cst_icon,cst_prin
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: cin_con, ctr_con, cdt_con, cpr_con
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: gen_name,alg_name, ss_name
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: zoo_name,epi_name, mac_name, bod_name
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: br_name,wb_name, tr_name
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: pipe_name,spill_name, gate_name, pump_name  
  CHARACTER(72)                                :: CONFNwrite, selfnwrite, habfnwrite
  CHARACTER(72),  ALLOCATABLE, DIMENSION(:)     :: bthline1,bthline2,vprline1,vprline2
  CHARACTER(72),  ALLOCATABLE, DIMENSION(:)    :: bthfnwrite, vprfnwrite
  CHARACTER(43),     ALLOCATABLE, DIMENSION(:)       :: CNAME4, cdunit
  logical vertprf, bthchng,vprchng, w2select, w2hab
  INTEGER :: NSTT,NEPTT
end module extra

! CE-QUAL-W2 computations
PROGRAM control_file_read_write
  USE MAIN
  USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC
  USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS
  USE RSTART
  use macrophytec; use porosityc; use zooplanktonc  !v3.5
  use extra
  use ifport, only:chdir  
  IMPLICIT NONE
  INTEGER(4) length,istatus,IERR
  character*512 dirc  
  logical S_OPEN, S_EXIST,h_open, h_exist,changes_exist
  integer s_number,h_number, LL2, LL4
   
  w2select=.false.
  w2hab=.false.
  
call get_command_argument(1,dirc,length,istatus)  
dirc=TRIM(dirc)

if(length /= 0)then
    istatus=CHDIR(dirc)
    select case(istatus)
      case(2)  ! ENOENT
        write(*,*)'The directory does not exist:',dirc
      case(20)   ! ENOTDIR
        write(*,*)'This is not a directory:', dirc
      case(0)    ! no error
    end select  
endif

     vertprf=.true.
     bthchng=.true.
     vprchng=.true.
! reading input files
    OPEN (CON,FILE=CONFN,STATUS='OLD',IOSTAT=IERR)
    IF (IERR /= 0) THEN
    CLOSE(CON)
    CONFN='w2_con.csv'
    OPEN (CON,FILE=CONFN,STATUS='OLD',IOSTAT=I)
       IF (I /= 0) THEN
           OPEN(200,FILE='CONVERTER_ERRORS.TXT',STATUS='UNKNOWN')
           WRITE(200,*)'Neither w2_con.npt nor w2_con.csv found'
       STOP
       ENDIF
       CALL READ_CONTROL_FILECSV
  ELSE 
     call read_control_fileNPT            !also reads bathymetry and vertical profile file(s)
  ENDIF
  

! writing files with changes
     CONFNwrite='w2_con45.csv'
     call control_file_csv1
     call control_file_csv2
     CALL WRITE_CSV
     call WRITE_CSV2
      
  STOP
    END PROGRAM control_file_read_write
!***************************************
    SUBROUTINE READ_CONTROL_FILECSV
    USE MAIN
    USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  use macrophytec; use porosityc; use zooplanktonc  
  use initialvelocity
  use extra  
  IMPLICIT NONE
  
  character*1 char1,ichar1
  character*8 AID
  character*72 linetest
  CHARACTER(8)  :: IBLANK,DMULT
  character(2)   :: ichar2
  INTEGER :: NPT, NNBP, NCBP,NINTERNAL,NUP,KTMAX,NTR1,JJ,NALT,NZPTT,NMCTT,NDUM,NIDUM
  
  ALLOCATE (TITLE(11))
  READ (CON,*)
  READ (CON,*)
  READ (CON,*)
  DO J=1,10
  READ (CON,*)TITLE(J)
  ENDDO
  READ (CON,*)
  READ (CON,*) 
  
  READ (CON,*) NWB, NBR, IMX, KMX, NPROC, CLOSEC; CLOSEC=ADJUSTR(CLOSEC)                   !'(A,5I0,A)'  
  READ (CON,*)
  READ (CON,*)
  
  READ (CON,*) NTR, NST, NIW, NWD, NGT, NSP, NPI, NPU
  READ (CON,*)
  READ (CON,*)  
  
  READ (CON,*) NGC, NSS, NAL, NEP, NBOD, NMC, NZP  
  READ (CON,*)
  READ (CON,*) 
  
  READ (CON,*) NOD,SELECTC,HABTATC,ENVIRPC,AERATEC,INITUWL           !'(I0,5(A))'   
  SELECTC=ADJUSTR(SELECTC);HABTATC=ADJUSTR(HABTATC);ENVIRPC=ADJUSTR(ENVIRPC);AERATEC=ADJUSTR(AERATEC);INITUWL=ADJUSTR(INITUWL)    

  IF(NTR==0)THEN
      NTR1=1
  ELSE
      NTR1=NTR
  ENDIF
! Constituent numbers
  NTDS  = 1
  NGCS  = 2
  NGCE  = NGCS+NGC-1
 
  NSSS  = NGCE+1
  NSSE  = NSSS+NSS-1
 
  NPO4  = NSSE+1
  NNH4  = NPO4+1
  NNO3  = NNH4+1
  NDSI  = NNO3+1
  NPSI  = NDSI+1
  NFE   = NPSI+1
  NLDOM = NFE+1
  NRDOM = NLDOM+1
  NLPOM = NRDOM+1
  NRPOM = NLPOM+1
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

  nbodps=ibod
  do jcb=1,nbod
    nbodp(jcb)=ibod
    ibod=ibod+1
  end do
  nbodpe=ibod-1

  nbodns=ibod
  do jcb=1,nbod
    nbodn(jcb)=ibod
    ibod=ibod+1
  end do
  nbodne=ibod-1
 
  else
    nbodcs=1; nbodce=1; nbodps=1; nbodpe=1; nbodns=1; nbodne=1     ! SW 11/21/2016
  
end if
!    NBODE = NBODS+NBOD-1
  NBODE = NBODS+NBOD*3-1
  NAS   = NBODE+1
  NAE   = NAS+NAL-1

  NDO   = NAE+1
  NTIC  = NDO+1
  NALK  = NTIC+1
  NZOOS = NALK+1
  NZOOE = NZOOS+NZP-1
 
  NLDOMP= NZOOE+1
  NRDOMP= NLDOMP+1
  NLPOMP= NRDOMP+1
  NRPOMP= NLPOMP+1
  NLDOMN= NRPOMP+1
  NRDOMN= NLDOMN+1
  NLPOMN= NRDOMN+1
  NRPOMN= NLPOMN+1

! Constituent, tributary, and widthdrawal totals

  NCT  = NRPOMN
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
    
  ALLOCATE (EXC(NWB),    EXIC(NWB))
  ALLOCATE (CUNIT(NCT),  CNAME(NCT),  CNAME1(NCT), CNAME2(NCT), CDNAME(NDC), CDMULT(NDC), CMULT(NCT), CUNIT1(NCT), CUNIT2(NCT))    ! SW 1/16/04 8/25/05
  ALLOCATE (HNAME(NHY),CDNAME2(NDC))
  ALLOCATE (SLTRC(NWB),  THETA(NWB),  FRICC(NWB))
  ALLOCATE (WINDC(NWB),  QINC(NWB),   QOUTC(NWB),  HEATC(NWB),  SLHTC(NWB))
  ALLOCATE (QINIC(NBR),  DTRIC(NBR),  TRIC(NTR),   WDIC(NWD),   HDIC(NBR),   METIC(NWB))
  ALLOCATE (VBC(NWB),    EBC(NWB),    MBC(NWB),    PQC(NWB),    EVC(NWB),    PRC(NWB))
  ALLOCATE (ICEC(NWB),   SLICEC(NWB), ICETHI(NWB),   ALBEDO(NWB), HWI(NWB),    BETAI(NWB),  GAMMAI(NWB), ICEMIN(NWB), ICET2(NWB))
  ALLOCATE (EXH2O(NWB),  BETA(NWB),   EXOM(NWB),   EXSS(NWB),   SROC(NWB), Z0(NWB))
  ALLOCATE (CBHE(NWB),   TSED(NWB),   FI(NWB),     TSEDF(NWB),  AFW(NWB),    BFW(NWB),    CFW(NWB),    WINDH(NWB),  RHEVC(NWB))
  ALLOCATE (SEDCC(NWB),   SEDPRC(NWB),  SEDCI(NWB),  SDK(NWB),  FSOD(NWB),   FSED(NWB), seds(nwb), sedb(nwb),DYNSEDK(NWB))   ! SW 6/1/07
  ALLOCATE (SODT1(NWB),  SODT2(NWB),  SODK1(NWB),  SODK2(NWB))
  ALLOCATE (FETCHC(NWB))
  ALLOCATE (SNPC(NWB),   PRFC(NWB),   SPRC(NWB),   CPLC(NWB),   VPLC(NWB),   FLXC(NWB),   SCRC(NWB))
  ALLOCATE (NSNP(NWB),   NPRF(NWB),   NSPR(NWB),   NCPL(NWB),   NVPL(NWB),   NFLX(NWB),   NSCR(NWB), TECPLOT(NWB))
  ALLOCATE (NISNP(NWB),  NIPRF(NWB),  NISPR(NWB))
  ALLOCATE (SNPFN(NWB),  PRFFN(NWB),  SPRFN(NWB),  VPLFN(NWB),  CPLFN(NWB),  FLXFN(NWB))
  ALLOCATE (BTHFN(NWB),  METFN(NWB),  VPRFN(NWB),  LPRFN(NWB))
  ALLOCATE (QINFN(NBR),  TINFN(NBR),  CINFN(NBR),  CDHFN(NBR),  QDTFN(NBR),  TDTFN(NBR),  CDTFN(NBR),  PREFN(NBR),  TPRFN(NBR))
  ALLOCATE (CPRFN(NBR),  EUHFN(NBR),  TUHFN(NBR),  CUHFN(NBR),  EDHFN(NBR),  TDHFN(NBR),  QOTFN(NBR))   ! SW 3/28/13
  ALLOCATE (EXTFN(NWB))                                                                                                !TC 12/12/01
  ALLOCATE (BODP(NBOD),  BODN(NBOD),  BODC(NBOD))                                                                      !TC 01/15/02
  ALLOCATE (DTRC(NBR),   ALPHA(NBR))
  ALLOCATE (AX(NWB),     WTYPEC(NWB), GRIDC(NWB),JBDN(NWB),   AZSLC(NWB),  AZMAX(NWB),  KBMAX(NWB),  AZC(NWB),TKEBC(NWB),EROUGH(NWB),ARODI(NWB),STRICK(NWB),TKELATPRDCONST(NWB),IMPTKE(NWB))
  ALLOCATE (VISC(NWB),   CELC(NWB), DLTADD(NWB))
  ALLOCATE (BTH(NWB),    VPR(NWB),    LPR(NWB))
  ALLOCATE (T2I(NWB),    KTWB(NWB),DXI(NWB))                                                                                    !TC 06/11/02

  ALLOCATE (LAT(NWB),    LONGIT(NWB),   ELBOT(NWB),  BS(NWB),     BE(NWB))  ! SW 4/6/2015
  !ALLOCATE (XC(IMX),YC(IMX),Xnw(IMX),Ynw(IMX),xw(IMX),yw(IMX),Xsw(IMX),Ysw(IMX),Xse(IMX),Yse(IMX),Xne(IMX),Yne(IMX),xe(IMX),ye(IMX),XCN(IMX),YCN(IMX),XCS(IMX),YCS(IMX))
  ALLOCATE (REAERC(NWB), NEQN(NWB))
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
  ALLOCATE (CAC(NCT))
  ALLOCATE (RCOEF1(NWB), RCOEF2(NWB), RCOEF3(NWB), RCOEF4(NWB))
  ALLOCATE (TTLB(IMX),   TTRB(IMX),   CLLB(IMX),   CLRB(IMX),   SRLB1(IMX),  SRRB1(IMX),   SRLB2(IMX),  SRRB2(IMX))
  ALLOCATE (VOLB(NBR),   VOLG(NWB))
  ALLOCATE (XBR(NBR))
  ALLOCATE (US(NBR),     DS(NBR),     CUS(NBR),    UHS(NBR),    DHS(NBR),    UQB(NBR),     DQB(NBR))
  ALLOCATE (NL(NBR),     NPOINT(NBR), SLOPE(NBR),SLOPEC(NBR))
  ALLOCATE (NSTR(NBR))
 
  ALLOCATE (Z(IMX),      DLX(IMX),    KB(IMX))
  ALLOCATE (ICETH(IMX),  KTI(IMX),    ELWS(IMX),   SHADE(IMX))
  ALLOCATE (SOD(IMX),    PHI0(IMX),   FRIC(IMX))
  ALLOCATE (CDN(NDC,NWB))
  ALLOCATE (CN(NCT),     INCN(NCT,NBR),   TRCN(NCT,NTR))                                                               !SW 01/07/01
  ALLOCATE (SSS(NSS),    SEDRC(NSS),   TAUCR(NSS))                                                                      !TC 08/20/03 SW 1/16/04
  ALLOCATE (AG(NAL),     AR(NAL),     AE(NAL),     AM(NAL),     AS(NAL),     EXA(NAL),    ASAT(NAL),   AP(NAL),   AN(NAL))
  ALLOCATE (AC(NAL),   ASI(NAL),  ACHLA(NAL),  AHSP(NAL),   AHSN(NAL),   AHSSI(NAL))
  ALLOCATE (AT1(NAL),    AT2(NAL),    AT3(NAL),    AT4(NAL),    AK1(NAL),    AK2(NAL),    AK3(NAL),    AK4(NAL))
  ALLOCATE (EG(nept),     ER(nept),     EE(nept),     EM(nept),     EB(nept),     ESAT(nept),   EP(nept),     EN(nept))
  ALLOCATE (EC(nept),     ESI(nept),    ECHLA(nept),  EHSP(nept),   EHSN(nept),   EHSSI(nept),  EPOM(nept),   EHS(nept))
  ALLOCATE (ET1(nept),    ET2(nept),    ET3(nept),    ET4(nept),    EK1(nept),    EK2(nept),    EK3(nept),    EK4(nept))
  ALLOCATE (ENPR(nept),   ENEQN(nept),  O2ER(nept),   O2EG(nept))
  ALLOCATE (IWDO(NOD),   RSOD(NOD),   RSOF(NOD),   DLTD(NOD),   DLTF(NOD),   DLTMAX(NOD))
  ALLOCATE (IUPI(NPI),   IDPI(NPI),   EUPI(NPI),   EDPI(NPI),   WPI(NPI),    DLXPI(NPI),  FPI(NPI),    FMINPI(NPI), PUPIC(NPI), DYNPIPE(NPI))
  ALLOCATE (ETUPI(NPI),  EBUPI(NPI),  KTUPI(NPI),  KBUPI(NPI),  PDPIC(NPI),  ETDPI(NPI),  EBDPI(NPI),  KTDPI(NPI))
  ALLOCATE (KBDPI(NPI),  JBUPI(NPI),  JBDPI(NPI))
  ALLOCATE (PUSPC(NSP),  ETUSP(NSP),  EBUSP(NSP),  KTUSP(NSP),  KBUSP(NSP),  PDSPC(NSP), ETDSP(NSP), EBDSP(NSP))
  ALLOCATE (KTDSP(NSP),  KBDSP(NSP),  IUSP(NSP),   IDSP(NSP),   ESP(NSP),    A1SP(NSP),  B1SP(NSP),  A2SP(NSP))
  ALLOCATE (B2SP(NSP),   AGASSP(NSP), BGASSP(NSP), CGASSP(NSP), EQSP(NSP),   GASSPC(NSP))
  ALLOCATE (GTA1(NGT),   GTB1(NGT),   GTA2(NGT),   GTB2(NGT),   IUGT(NGT),   IDGT(NGT),   EGT(NGT), EGT2(NGT),  A1GT(NGT),   B1GT(NGT))
  ALLOCATE (G1GT(NGT),   A2GT(NGT),   B2GT(NGT),   G2GT(NGT),   PUGTC(NGT),  ETUGT(NGT),  EBUGT(NGT),  KTUGT(NGT))
  ALLOCATE (KBUGT(NGT),  PDGTC(NGT),  ETDGT(NGT),  EBDGT(NGT),  KTDGT(NGT),  KBDGT(NGT),  AGASGT(NGT), BGASGT(NGT), CGASGT(NGT))
  ALLOCATE (GASGTC(NGT), EQGT(NGT),   JBUGT(NGT),  JBDGT(NGT))
  ALLOCATE (LATGTC(NGT), LATSPC(NSP), LATPIC(NPI), LATPUC(NPU), DYNGTC(NGT),GTIC(NGT))
  ALLOCATE (IUPU(NPU),   IDPU(NPU),   EPU(NPU),    STRTPU(NPU), ENDPU(NPU),  EONPU(NPU),  EOFFPU(NPU), QPU(NPU),    PPUC(NPU))
  ALLOCATE (ETPU(NPU),   EBPU(NPU),   KTPU(NPU),   KBPU(NPU))
  ALLOCATE (IWR(NIW),    KTWR(NIW),   KBWR(NIW))
  ALLOCATE (IWD(NWDT),   JBWD(NWD),   EWD(NWDT),   KTWD(NWDT),  KBWD(NWDT))
  ALLOCATE (ITR(NTRT),   QTRFN(NTR),  TTRFN(NTR),  CTRFN(NTR),  ELTRT(NTRT),  ELTRB(NTRT),  TRC(NTRT))
  
  ALLOCATE (JBTR(NTR),    NACTR(NTR),  NACIN(NBR),  NACPR(NWB),  NACDT(NBR))
  ALLOCATE (TSRD(NOD),   TSRF(NOD),   WDOD(NOD),   WDOF(NOD))
  ALLOCATE (T2(KMX,IMX),DYNPUMP(NPU))
  ALLOCATE (ITSR(IMX*KMX),   ETSR(IMX*KMX))
  ALLOCATE (ESTRT(NSTT,NBR),   WSTRT(NSTT,NBR),   SINKCT(NSTT,NBR),  STRIC(NSTT,NBR))
  ALLOCATE (KTSWT(NSTT,NBR),  KBSWT(NSTT,NBR))
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
  ALLOCATE (CINBRC(NCT,NBR), CTRTRC(NCT,NTR1), CDTBRC(NCT,NBR), CPRBRC(NCT,NBR), CPRWBC(NCT,NWB), KFWBC(NFL,NWB),  HPRWBC(NHY,NWB))
  ALLOCATE (WSC(IMX),    SNPD(NOD,NWB),   SNPF(NOD,NWB),   SCRD(NOD,NWB),   SCRF(NOD,NWB),   PRFD(NOD,NWB))
  ALLOCATE (PRFF(NOD,NWB),   SPRD(NOD,NWB),   SPRF(NOD,NWB),   CPLD(NOD,NWB),   CPLF(NOD,NWB))
  ALLOCATE (VPLD(NOD,NWB),   VPLF(NOD,NWB),   FLXD(NOD,NWB),   FLXF(NOD,NWB))

  ALLOCATE (TVP(KMX,NWB),    H(KMX,NWB))
  ALLOCATE (B(KMX,IMX),      CONV(KMX,IMX),   EL(KMX,IMX))
  ALLOCATE (CDWBC(NDC,NWB))
  ALLOCATE (C2I(NCT,NWB))
  ALLOCATE (EPIC(NWB,NEPTT), EPICI(NWB,NEPTT),EPIPRC(NWB,NEPTT))                                                       ! cb 6/16/06
  ALLOCATE (ISNP(IMX,NWB),   IPRF(IMX,NWB),   ISPR(IMX,NWB))
  ALLOCATE (CVP(KMX,NCT,NWB), EPIVP(KMX,IMX,NEPT),EPD(KMX,IMX,NEPT))
  ALLOCATE (C2(KMX,IMX,NCT), CD(KMX,IMX,NDC), CG(KMX,IMX,NGC), ALG(KMX,IMX,NAL), SS(KMX,IMX,NSS), CBOD(KMX,IMX,NBOD))

  ALLOCATE (OPEN_VPR(NWB),     OPEN_LPR(NWB),     ISO_TEMP(NWB),    VERT_TEMP(NWB))
  ALLOCATE (UH_EXTERNAL(NBR),  DH_EXTERNAL(NBR),   UH_INTERNAL(NBR), DH_INTERNAL(NBR))
  ALLOCATE (UQ_EXTERNAL(NBR),  DQ_EXTERNAL(NBR),   UQ_INTERNAL(NBR), DQ_INTERNAL(NBR))
  ALLOCATE (ISO_CONC(NCT,NWB), VERT_CONC(NCT,NWB), LONG_TEMP(NWB),   LONG_CONC(NCT,NWB))
  ALLOCATE (ICE_CALC(NWB),     PRECIPITATION(NWB))
  ALLOCATE (ZERO_SLOPE(NWB))                                                                                           !SW 10/17/02
  ALLOCATE (INTERP_INFLOW(NBR))                                                                                       ! cb 9/7/2010
  ALLOCATE (CMIN(NCT),   CMAX(NCT),   HYMIN(NHY),  HYMAX(NHY),  CDMIN(NDC),  CDMAX(NDC))                               ! SW 1/16/04
  ALLOCATE (FMTH(NHY),   HMULT(NHY),  FMTC(NCT),   FMTCD(NDC))                                                         ! SW 1/16/04
  ALLOCATE (CPLTC(NCT),  HPLTC(NHY),  CDPLTC(NDC))                                                                     ! SW 1/16/04
  ALLOCATE (zg(NZPt),zm(NZPt),zeff(NZPt),prefp(NZPt),zr(NZPt),zoomin(NZPt),zs2p(NZPt),exz(NZPt),PREFZ(NZPtT,nzpt),o2zr(nzpt))
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
  allocate (sedp(kmx,imx),sedn(kmx,imx),SED(KMX,IMX))
  allocate (sdkv(kmx,imx),seddktot(kmx,imx))
  allocate (sedcip(nwb),sedcin(nwb),sedcic(nwb),sedcis(nwb))
  ALLOCATE (cbods(NBOD), cbodns(kmx,imx), sedcb(kmx,imx), sedcbp(kmx,imx), sedcbn(kmx,imx), sedcbc(kmx,imx))
  allocate  (print_macrophyte(nwb,nmct), macrophyte_calc(nwb,nmct),macwbc(nwb,nmctT),conv2(kmx,kmx),mprwbc(nwb,nmctT))
  allocate  (mac(kmx,imx,nmct), macrc(kmx,kmx,imx,nmct),mact(kmx,kmx,imx), macrm(kmx,kmx,imx,nmct), macss(kmx,kmx,imx,nmct))
  allocate  (mgr(kmx,kmx,imx,nmct),mmr(kmx,imx,nmct), mrr(kmx,imx,nmct))
  allocate  (smacrc(kmx,kmx,imx,nmct), smacrm(kmx,kmx,imx,nmct))
  allocate  (smact(kmx,kmx,imx), smac(kmx,imx,nmct))
  allocate  (mt1(nmct),mt2(nmct),mt3(nmct),mt4(nmct),mk1(nmct),mk2(nmct),mk3(nmct),mk4(nmct),mg(nmct),mr(nmct),mm(nmct))
  allocate  (mbmp(nmct), mmax(nmct), cdDRAG(nmct),dWv(nmct),dwsa(nmct),anorm(nmct))
  
  allocate  (mp(nmct), mn(nmct), mc(nmct),psed(nmct),nsed(nmct),mhsp(nmct),mhsn(nmct),mhsc(nmct),msat(nmct),exm(nmct))
  allocate  (O2MG(nmct), O2MR(nmct),  LRPMAC(nmct),  MPOM(nmct))
  allocate  (kticol(imx),armac(imx),macwbci(nwb,nmctT))
  allocate  (macmbrs(nbr,nmct), macmbrt(nbr,nmct),ssmacmb(nbr,nmct))
  allocate  (cw(kmx,imx), bic(kmx,imx))
  allocate  (mactrmr(kmx,imx,nmct), mactrmf(kmx,imx,nmct),mactrm(kmx,imx,nmct))
  allocate  (mlfpr(kmx,kmx,imx,nmct))
  allocate  (mllim(kmx,kmx,imx,nmct), mplim(kmx,imx,nmct),mclim(kmx,imx,nmct),mnlim(kmx,imx,nmct))
  ALLOCATE  (GAMMAj(kmx,KMX,IMX))	
  allocate (por(kmx,imx),VOLKTi(imx),VOLi(Kmx,Imx),vstem(kmx,imx,nmct),vstemkt(imx,nmct),sarea(nmct))
  ALLOCATE (IWIND(NWB))
  allocate (cbodp(kmx,imx,nbod), cbodn(kmx,imx,nbod))
  ALLOCATE (ISO_EPIPHYTON(NWB,NEPt), VERT_EPIPHYTON(NWB,NEPt), LONG_EPIPHYTON(NWB,NEPt), EPIPHYTON_CALC(NWB,NEPt))         !TC 10/25/02
  ALLOCATE (ISO_SEDIMENT(NWB),    VERT_SEDIMENT(NWB),   LONG_SEDIMENT(NWB))
  ALLOCATE (DDUM(NOD),FDUM(NOD), IDUM(IMX),KFNAME(NFL), KFNAME2(NFL))

  READ (CON,*)
  READ (CON,*) 
    
  READ (CON,*)  TMSTRT,   TMEND,    YEAR    
  READ (CON,*)
  READ (CON,*) 
    
  READ (CON,*)  NDLT,     DLTMIN, DLTINTER; DLTD=0.0; DLTINTER=ADJUSTR(DLTINTER)
  READ (CON,*)
  READ (CON,*)
  READ (CON,*)  (DLTD(J), J =1,NDLT)
  READ (CON,*)
  READ (CON,*)
  READ (CON,*)  (DLTMAX(J), J =1,NDLT)
  READ (CON,*)
  READ (CON,*) 
    
  READ (CON,*)  (DLTF(J),   J =1,NDLT)
  READ (CON,*)
  READ (CON,*)
  READ (CON,*)      (VISC(JW), JW=1,NWB)
  READ (CON,*)      (CELC(JW), JW=1,NWB)
  READ (CON,*)      (DLTADD(JW), JW=1,NWB)
  VISC=ADJUSTR(VISC);   CELC=ADJUSTR(CELC);   DLTADD=ADJUSTR(DLTADD)

! Grid definition cards
  READ (CON,*)
  READ (CON,*) 

  READ (CON,*)  (US(JB),    JB=1,NBR)
  READ (CON,*)  (DS(JB),    JB=1,NBR)
  READ (CON,*)  (UHS(JB),   JB=1,NBR)
  READ (CON,*)  (DHS(JB),   JB=1,NBR)
  READ (CON,*)  (NL(JB),    JB=1,NBR)
  READ (CON,*)  (SLOPE(JB), JB=1,NBR)
  READ (CON,*)  (SLOPEC(JB),JB=1,NBR)
  READ (CON,*)
  READ (CON,*) 

  READ (CON,*)   (LAT(JW),                     JW=1,NWB)
  READ (CON,*)   (LONGIT(JW),                    JW=1,NWB)
  READ (CON,*)   (ELBOT(JW),                   JW=1,NWB)
  READ (CON,*)   (BS(JW),                      JW=1,NWB)
  READ (CON,*)   (BE(JW),                      JW=1,NWB)
  READ (CON,*)   (JBDN(JW),                    JW=1,NWB)
  READ (CON,*)
  READ (CON,*)

! Initial condition cards

  READ (CON,*)      (T2I(JW),                              JW=1,NWB)
  READ (CON,*)      (ICETHI(JW),                             JW=1,NWB)
  READ (CON,*)      (WTYPEC(JW),                           JW=1,NWB); WTYPEC=ADJUSTR(WTYPEC)
  READ (CON,*)      (GRIDC(JW),                           JW=1,NWB); GRIDC=ADJUSTR(GRIDC)
  READ (CON,*)
  READ (CON,*)

  READ (CON,*)            (VBC(JW),              JW=1,NWB); VBC=ADJUSTR(VBC)
  READ (CON,*)            (EBC(JW),              JW=1,NWB); EBC=ADJUSTR(EBC)
  READ (CON,*)            (MBC(JW),              JW=1,NWB); MBC=ADJUSTR(MBC)   
  READ (CON,*)            (PQC(JW),              JW=1,NWB); PQC=ADJUSTR(PQC)
  READ (CON,*)            (EVC(JW),              JW=1,NWB); EVC=ADJUSTR(EVC)
  READ (CON,*)            (PRC(JW),              JW=1,NWB); PRC=ADJUSTR(PRC)
  READ (CON,*)
  READ (CON,*)

  READ (CON,*)            (WINDC(JW),             JW=1,NWB); WINDC=ADJUSTR(WINDC)
  READ (CON,*)            (QINC(JW),              JW=1,NWB); QINC=ADJUSTR(QINC)
  READ (CON,*)            (QOUTC(JW),             JW=1,NWB); QOUTC=ADJUSTR(QOUTC)
  READ (CON,*)            (HEATC(JW),             JW=1,NWB); HEATC=ADJUSTR(HEATC)
  READ (CON,*)
  READ (CON,*)

  READ (CON,*)            (QINIC(JB),               JB=1,NBR); QINIC=ADJUSTR(QINIC)
  READ (CON,*)            (DTRIC(JB),               JB=1,NBR); DTRIC=ADJUSTR(DTRIC)
  READ (CON,*)            (HDIC(JB),                JB=1,NBR); HDIC=ADJUSTR(HDIC)
  READ (CON,*)
  READ (CON,*)

  READ (CON,*)      (SLHTC(JW),                     JW=1,NWB); SLHTC=ADJUSTR(SLHTC)
  READ (CON,*)      (SROC(JW),                      JW=1,NWB); SROC=ADJUSTR(SROC)
  READ (CON,*)      (RHEVC(JW),                     JW=1,NWB); RHEVC=ADJUSTR(RHEVC)
  READ (CON,*)      (METIC(JW),                     JW=1,NWB); METIC=ADJUSTR(METIC)
  READ (CON,*)      (FETCHC(JW),                    JW=1,NWB); FETCHC=ADJUSTR(FETCHC)
  READ (CON,*)      (AFW(JW),                       JW=1,NWB)
  READ (CON,*)      (BFW(JW),                       JW=1,NWB)
  READ (CON,*)      (CFW(JW),                       JW=1,NWB)
  READ (CON,*)      (WINDH(JW),                     JW=1,NWB)
  READ (CON,*)
  READ (CON,*)

  READ (CON,*)      (ICEC(JW),                      JW=1,NWB); ICEC=ADJUSTR(ICEC)
  READ (CON,*)      (SLICEC(JW),                    JW=1,NWB); SLICEC=ADJUSTR(SLICEC)
  READ (CON,*)      (ALBEDO(JW),                    JW=1,NWB)
  READ (CON,*)      (HWI(JW),                       JW=1,NWB)
  READ (CON,*)      (BETAI(JW),                     JW=1,NWB)
  READ (CON,*)      (GAMMAI(JW),                    JW=1,NWB)
  READ (CON,*)      (ICEMIN(JW),                    JW=1,NWB)
  READ (CON,*)      (ICET2(JW),                     JW=1,NWB)
  READ (CON,*)
  READ (CON,*)

  READ (CON,*)        (SLTRC(JW),                   JW=1,NWB); SLTRC=ADJUSTR(SLTRC)
  READ (CON,*)        (THETA(JW),                   JW=1,NWB)
  READ (CON,*)
  READ (CON,*)

  READ (CON,*)       (AX(JW),                       JW=1,NWB)
  READ (CON,*)       (DXI(JW),                       JW=1,NWB)
  READ (CON,*)       (CBHE(JW),                     JW=1,NWB)
  READ (CON,*)       (TSED(JW),                     JW=1,NWB)
  READ (CON,*)       (FI(JW),                       JW=1,NWB)
  READ (CON,*)       (TSEDF(JW),                    JW=1,NWB)
  READ (CON,*)       (FRICC(JW),                    JW=1,NWB); FRICC=ADJUSTR(FRICC)
  READ (CON,*)       (Z0(JW),                       JW=1,NWB)
  READ (CON,*)
  READ (CON,*)

  READ (CON,*)      (AZC(JW),    JW=1,NWB);AZC=ADJUSTR(AZC)          
  READ (CON,*)      (AZSLC(JW),  JW=1,NWB);AZSLC=ADJUSTR(AZSLC)          
  READ (CON,*)      (AZMAX(JW),  JW=1,NWB)          
  READ (CON,*)      (TKEBC(JW),  JW=1,NWB)          
  READ (CON,*)      (EROUGH(JW), JW=1,NWB)          
  READ (CON,*)      (ARODI(JW),  JW=1,NWB)          
  READ (CON,*)      (STRICK(JW), JW=1,NWB)          
  READ (CON,*)      (TKELATPRDCONST(JW),JW=1,NWB)   
  READ (CON,*)      (IMPTKE(JW),JW=1,NWB);IMPTKE=ADJUSTR(IMPTKE)          !FBC(JW),AZE(JW),ARODI(JW),STRCKLR(JW),BOUNDFR(JW),TKECAL(JW)
  READ (CON,*)
  READ (CON,*)

  READ (CON,*)            (NSTR(JB),      JB=1,NBR)
  READ (CON,*)            (DYNSTRUC(JB),  JB=1,NBR); DYNSTRUC=adjustr(DYNSTRUC)
  
  DO JS=1,NSTT
    READ (CON,*)        (STRIC(JS,JB),  JB=1,NBR)
  END DO
  STRIC=adjustr(STRIC)
  DO JS=1,NSTT
    READ (CON,*)        (KTSWT(JS,JB),  JB=1,NBR)
  END DO
  DO JS=1,NSTT
    READ (CON,*)        (KBSWT(JS,JB),  JB=1,NBR)
  END DO
  SINKCT=''
  DO JS=1,NSTT
    READ (CON,*)        (SINKCT(JS,JB), JB=1,NBR)
  END DO
  SINKCT=adjustr(SINKCT)
  DO JS=1,NSTT
    READ (CON,*)        (ESTRT(JS,JB),  JB=1,NBR)
  END DO
  DO JS=1,NSTT
    READ (CON,*)        (WSTRT(JS,JB),  JB=1,NBR)
  END DO
  READ (CON,*)
  READ (CON,*)

  READ (CON,*)  (IUPI(JP),  JP=1,NPI)
  READ (CON,*)  (IDPI(JP),  JP=1,NPI)
  READ (CON,*)  (EUPI(JP),  JP=1,NPI)
  READ (CON,*)  (EDPI(JP),  JP=1,NPI)
  READ (CON,*)  (WPI(JP),   JP=1,NPI)
  READ (CON,*)  (DLXPI(JP), JP=1,NPI)
  READ (CON,*)  (FPI(JP),   JP=1,NPI)
  READ (CON,*)  (FMINPI(JP),JP=1,NPI)
  READ (CON,*)  (LATPIC(JP),JP=1,NPI);LATPIC=adjustr(LATPIC)
  READ (CON,*)  (DYNPIPE(JP),JP=1,NPI);DYNPIPE=adjustr(DYNPIPE) 
  
  READ (CON,*)  (PUPIC(JP),  JP=1,NPI);PUPIC=adjustr(PUPIC)
  READ (CON,*)  (ETUPI(JP),  JP=1,NPI)
  READ (CON,*)  (EBUPI(JP),  JP=1,NPI)
  READ (CON,*)  (KTUPI(JP),  JP=1,NPI)
  READ (CON,*)  (KBUPI(JP),  JP=1,NPI)
  
  READ (CON,*)  (PDPIC(JP),  JP=1,NPI);PDPIC=adjustr(PDPIC)
  READ (CON,*)  (ETDPI(JP),  JP=1,NPI)
  READ (CON,*)  (EBDPI(JP),  JP=1,NPI)
  READ (CON,*)  (KTDPI(JP),  JP=1,NPI)
  READ (CON,*)  (KBDPI(JP),  JP=1,NPI)
  
  READ (CON,*)
  READ (CON,*)

  READ (CON,*) (IUSP(JS),    JS=1,NSP)
  READ (CON,*) (IDSP(JS),    JS=1,NSP)
  READ (CON,*) (ESP(JS),     JS=1,NSP)
  READ (CON,*) (A1SP(JS),    JS=1,NSP)
  READ (CON,*) (B1SP(JS),    JS=1,NSP)
  READ (CON,*) (A2SP(JS),    JS=1,NSP)
  READ (CON,*) (B2SP(JS),    JS=1,NSP)
  READ (CON,*) (LATSPC(JS),  JS=1,NSP);LATSPC=ADJUSTR(LATSPC)
  
  READ (CON,*) (PUSPC(JS),   JS=1,NSP);PUSPC=adjustr(PUSPC)
  READ (CON,*) (ETUSP(JS),   JS=1,NSP)
  READ (CON,*) (EBUSP(JS),   JS=1,NSP)
  READ (CON,*) (KTUSP(JS),   JS=1,NSP)
  READ (CON,*) (KBUSP(JS),   JS=1,NSP)
  
  READ (CON,*) (PDSPC(JS),   JS=1,NSP) ;PDSPC=ADJUSTR(PDSPC) 
  READ (CON,*) (ETDSP(JS),   JS=1,NSP)
  READ (CON,*) (EBDSP(JS),   JS=1,NSP)
  READ (CON,*) (KTDSP(JS),   JS=1,NSP)
  READ (CON,*) (KBDSP(JS),   JS=1,NSP)
  
  READ (CON,*) (GASSPC(JS),  JS=1,NSP);GASSPC=ADJUSTR(GASSPC)
  READ (CON,*) (EQSP(JS),    JS=1,NSP)
  READ (CON,*) (AGASSP(JS),  JS=1,NSP)
  READ (CON,*) (BGASSP(JS),  JS=1,NSP)
  READ (CON,*) (CGASSP(JS),  JS=1,NSP)
  
  READ (CON,*)
  READ (CON,*)

  READ (CON,*) (IUGT(JG),   JG=1,NGT)
  READ (CON,*) (IDGT(JG),   JG=1,NGT)
  READ (CON,*) (EGT(JG),    JG=1,NGT)
  READ (CON,*) (A1GT(JG),   JG=1,NGT)
  READ (CON,*) (B1GT(JG),   JG=1,NGT)
  READ (CON,*) (G1GT(JG),   JG=1,NGT)
  READ (CON,*) (A2GT(JG),   JG=1,NGT)
  READ (CON,*) (B2GT(JG),   JG=1,NGT)
  READ (CON,*) (G2GT(JG),   JG=1,NGT)
  READ (CON,*) (LATGTC(JG), JG=1,NGT);LATGTC=ADJUSTR(LATGTC)
  
  READ (CON,*) (GTA1(JG),   JG=1,NGT)  
  READ (CON,*) (GTB1(JG),   JG=1,NGT)  
  READ (CON,*) (GTA2(JG),   JG=1,NGT)  
  READ (CON,*) (GTB2(JG),   JG=1,NGT)  
  READ (CON,*) (DYNGTC(JG), JG=1,NGT);DYNGTC=ADJUSTR(DYNGTC)  
  READ (CON,*) (GTIC(JG),   JG=1,NGT);GTIC=ADJUSTR(GTIC)  
  
  READ (CON,*) (PUGTC(JG),  JG=1,NGT);PUGTC=ADJUSTR(PUGTC)
  READ (CON,*) (ETUGT(JG),  JG=1,NGT)
  READ (CON,*) (EBUGT(JG),  JG=1,NGT)
  READ (CON,*) (KTUGT(JG),  JG=1,NGT)
  READ (CON,*) (KBUGT(JG),  JG=1,NGT)
  READ (CON,*) (PDGTC(JG),  JG=1,NGT);PDGTC=ADJUSTR(PDGTC)
  READ (CON,*) (ETDGT(JG),  JG=1,NGT)
  READ (CON,*) (EBDGT(JG),  JG=1,NGT)
  READ (CON,*) (KTDGT(JG),  JG=1,NGT)
  READ (CON,*) (KBDGT(JG),  JG=1,NGT)
  
  READ (CON,*) (GASGTC(JG), JG=1,NGT);GASGTC=ADJUSTR(GASGTC)  
  READ (CON,*) (EQGT(JG),   JG=1,NGT)  
  READ (CON,*) (AGASGT(JG), JG=1,NGT)  
  READ (CON,*) (BGASGT(JG), JG=1,NGT)  
  READ (CON,*) (CGASGT(JG), JG=1,NGT)  
  
  READ (CON,*)
  READ (CON,*)

  
  READ (CON,*) (IUPU(JP),   JP=1,NPU)
  READ (CON,*) (IDPU(JP),   JP=1,NPU)
  READ (CON,*) (EPU(JP),    JP=1,NPU)
  READ (CON,*) (STRTPU(JP), JP=1,NPU)
  READ (CON,*) (ENDPU(JP),  JP=1,NPU)
  READ (CON,*) (EONPU(JP),  JP=1,NPU)
  READ (CON,*) (EOFFPU(JP), JP=1,NPU)
  READ (CON,*) (QPU(JP),    JP=1,NPU)
  READ (CON,*) (LATPUC(JP), JP=1,NPU);LATPUC=adjustr(LATPUC)
  READ (CON,*) (DYNPUMP(JP),JP=1,NPU);DYNPUMP=ADJUSTR(DYNPUMP)
  READ (CON,*) (PPUC(JP),   JP=1,NPU);PPUC=ADJUSTR(PPUC)
  READ (CON,*) (ETPU(JP),   JP=1,NPU)
  READ (CON,*) (EBPU(JP),   JP=1,NPU)
  READ (CON,*) (KTPU(JP),   JP=1,NPU)
  READ (CON,*) (KBPU(JP),   JP=1,NPU)
  READ (CON,*)
  READ (CON,*)

  READ (CON,*)         (IWR(JW),    JW=1,NIW)
  READ (CON,*)        (KTWR(JW),   JW=1,NIW)               
  READ (CON,*)        (KBWR(JW),   JW=1,NIW)   
  READ (CON,*)
  READ (CON,*)

  READ (CON,*)         (WDIC(JW),   JW=1,NWD);WDIC=adjustr(WDIC)
  READ (CON,*)         (IWD(JW),    JW=1,NWD)
  READ (CON,*)         (EWD(JW),    JW=1,NWD)
  READ (CON,*)         (KTWD(JW),   JW=1,NWD)
  READ (CON,*)         (KBWD(JW),   JW=1,NWD); TRC= '      '  ! SW 9/27/13 INITIALIZATION SINCE ALLOCATION IS TO NTRT
  READ (CON,*)
  READ (CON,*)

  READ (CON,*)        (TRC(JT),    JT=1,NTR);TRC=ADJUSTR(TRC)
  READ (CON,*)        (TRIC(JT),   JT=1,NTR);TRIC=ADJUSTR(TRIC)
  READ (CON,*)        (ITR(JT),    JT=1,NTR)
  READ (CON,*)        (ELTRT(JT),  JT=1,NTR)
  READ (CON,*)        (ELTRB(JT),  JT=1,NTR)
  READ (CON,*)        (QTRFN(JT),  JT=1,NTR)
  READ (CON,*)        (TTRFN(JT),  JT=1,NTR)
  READ (CON,*)        (CTRFN(JT),  JT=1,NTR)  
  READ (CON,*)
  READ (CON,*)
  
  READ (CON,*)      (DTRC(JB),   JB=1,NBR);DTRC=adjustr(DTRC)
  READ (CON,*)
  READ (CON,*)

   ! Output control cards (excluding constituents)
  DO JH=1,NHY
    READ (CON,*)    HNAME(JH),  FMTH(JH),  HMULT(JH),(HPRWBC(JH,JW),JW=1,NWB);HPRWBC=ADJUSTR(HPRWBC)
  END DO
  READ (CON,*)
  READ (CON,*)

  READ (CON,*)        SNPC(1)
  READ (CON,*)        NSNP(1)
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
    READ (CON,*)      (SNPD(J,1),J=1,NSNP(1))
    READ (CON,*)      (SNPF(J,1),J=1,NSNP(1))
    DO J=1,NSNP(1)
    SNPD(J,2:NWB)=SNPD(J,1)
    SNPF(J,2:NWB)=SNPF(J,1)
    ENDDO

  READ (CON,*)
  READ (CON,*) 
  READ (CON,*)         SCRC(1);  SCRC(1)=ADJUSTR(SCRC(1))    
  READ (CON,*)         NSCR(1)
  READ (CON,*)          (SCRD(J,1),J=1,NSCR(1))
  READ (CON,*)          (SCRF(J,1),J=1,NSCR(1))
  READ (CON,*)
  READ (CON,*)  
  
  IF(NWB > 1)THEN
      SCRC(2:NWB)=SCRC(1)
      NSCR(2:NWB)=NSCR(1)
  DO J=1,NSCR(1)
      SCRD(J,2:NWB)=SCRD(J,1)
      SCRF(J,2:NWB)=SCRF(J,1)
  ENDDO
  ENDIF
   
  CDUM='     OFF'
  READ (CON,*)      CDUM;CDUM=ADJUSTR(CDUM)         
  READ (CON,*)      NDUM       ! PRFC(1), NPRF(1), NIPRF(1)    
  READ (CON,*)      NIDUM       ! PRFC(1), NPRF(1), NIPRF(1)    
  
  READ (CON,*)      (DDUM(J),J=1,NDUM)    !(PRFD(J,1),J=1,NPRF(1))
  READ (CON,*)      (FDUM(J),J=1,NDUM)    !(PRFF(J,1),J=1,NPRF(1))
  READ (CON,*)      (IDUM(J),J=1,NIDUM)   !(IPRF(J,1),J=1,NIPRF(1))
  READ (CON,*)
  READ (CON,*)  
  
  PRFC='     OFF'
  NPRF=0
  NIPRF=0
  IF(CDUM == '      ON')THEN                      !NWB>1 .AND. 
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
  READ (CON,*)      CDUM;CDUM=ADJUSTR(CDUM)         
  READ (CON,*)      NDUM         
  READ (CON,*)      NIDUM          
  READ (CON,*)      (DDUM(J),J=1,NDUM)    
  READ (CON,*)      (FDUM(J),J=1,NDUM)    
  READ (CON,*)      (IDUM(J),J=1,NIDUM)   
  READ (CON,*)
  READ (CON,*)  
  
  SPRC='     OFF'
  NSPR=0
  NISPR=0
  IF(CDUM == '      ON')THEN               !NWB>1 .AND. 
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
                  SPRC(JW)='      ON'
              ENDIF
          ENDDO
      ENDDO     
  ENDIF

  VPLC='     OFF'
  NVPL=0
  READ (CON,*)    VPLC(1);VPLC(1)=ADJUSTR(VPLC(1))
  READ (CON,*)    NVPL(1)
  READ (CON,*)   (VPLD(J,1), J=1,NVPL(1))
  READ (CON,*)   (VPLF(J,1), J=1,NVPL(1))
  READ (CON,*)
  READ (CON,*)  

  READ (CON,*)      CPLC(1);CPLC(1)=ADJUSTR(CPLC(1))
  READ (CON,*)      NCPL(1)
  READ (CON,*)      TECPLOT(1);TECPLOT(1)=ADJUSTR(TECPLOT(1))
  
  CPLC(2:NWB)=CPLC(1)
  NCPL(2:NWB)=NCPL(1)
  TECPLOT(2:NWB)=TECPLOT(1)
  READ (CON,*)    (CPLD(J,1), J=1,NCPL(1))
  READ (CON,*)    (CPLF(J,1), J=1,NCPL(1))
  READ (CON,*)
  READ (CON,*) 
  DO J=1,NCPL(1)
    CPLD(J,2:NWB)=CPLD(J,1)
    CPLF(J,2:NWB)=CPLF(J,1)
  ENDDO    
    
  READ (CON,*)        FLXC(1);FLXC(1)=ADJUSTR(FLXC(1))
  READ (CON,*)        NFLX(1)
  FLXC(2:NWB)=FLXC(1)
  NFLX(2:NWB)=NFLX(1)

  READ (CON,*)       (FLXD(J,1), J=1,NFLX(1))
  READ (CON,*)       (FLXF(J,1), J=1,NFLX(1))
  READ (CON,*)
  READ (CON,*) 
  DO J=1,NFLX(1)
  FLXD(J,2:NWB)=FLXD(J,1)
  FLXF(J,2:NWB)=FLXF(J,1) 
  ENDDO

  READ (CON,*)           TSRC;TSRC=ADJUSTR(TSRC)
  READ (CON,*)           NTSR
  READ (CON,*)          NIKTSR
  READ (CON,*)          TSRFN
  
  READ (CON,*)        (TSRD(J), J=1,NTSR)
  READ (CON,*)        (TSRF(J), J=1,NTSR)
  READ (CON,*)        (ITSR(J), J=1,NIKTSR)
  READ (CON,*)        (ETSR(J), J=1,NIKTSR)
  READ (CON,*)
  READ (CON,*) 

  READ (CON,*)           WDOC; WDOC=ADJUSTR(WDOC)   
  READ (CON,*)           NWDO 
  READ (CON,*)           NIWDO
  READ (CON,*)           WDOFN

!  ALLOCATE (IWDO(MAX(1,NIWDO)))
  READ (CON,*)       (WDOD(J), J=1,NWDO)
  READ (CON,*)       (WDOF(J), J=1,NWDO)
  READ (CON,*)       (IWDO(J), J=1,NIWDO)
  READ (CON,*)
  READ (CON,*) 

  READ (CON,*)         RSOC;RSOC=ADJUSTR(RSOC)
  READ (CON,*)         NRSO
  READ (CON,*)         RSIC;RSIC=ADJUSTR(RSIC)
  READ (CON,*)         RSIFN

  RSOD=0.0 ! SW 9/27/13 INITIALIZE SINCE ALLOCATED AS NOD BUT ONLY NRSO USED
  READ (CON,*)           (RSOD(J), J=1,NRSO)
  READ (CON,*)           (RSOF(J), J=1,NRSO)
  READ (CON,*)
  READ (CON,*)    

    
    ! Constituent control cards

  READ (CON,*)           CCC, LIMC, CUF,CO2PPM,CO2YRLY;CCC=ADJUSTR(CCC);LIMC=ADJUSTR(LIMC);CO2YRLY=ADJUSTR(CO2YRLY) 
  READ (CON,*)
  READ (CON,*) 
  DO JC=1,NCT
  READ (CON,*)CNAME2(JC),CNAME(JC),CAC(JC), FMTC(JC), CMULT(JC), (C2I(JC,JW), JW=1,NWB),(CPRWBC(JC,JW), JW=1,NWB), (CINBRC(JC,JB), JB=1,NBR),(CTRTRC(JC,JT), JT=1,NTR1), (CDTBRC(JC,JB), JB=1,NBR), (CPRBRC(JC,JB), JB=1,NBR)  
  ENDDO
  CAC=ADJUSTR(CAC);CPRWBC=ADJUSTR(CPRWBC);CINBRC=ADJUSTR(CINBRC);CTRTRC=ADJUSTR(CTRTRC);CDTBRC=ADJUSTR(CDTBRC);CPRBRC=ADJUSTR(CPRBRC)
  READ (CON,*)
  READ (CON,*) 

  DO JD=1,NDC2
  READ (CON,*)  CDNAME2(JD),CDNAME(JD),FMTCD(JD),CDMULT(JD),(CDWBC(JD,JW), JW=1,NWB)
  ENDDO
  CDWBC=ADJUSTR(CDWBC)
  READ (CON,*)
  READ (CON,*) 

  KFWBC   ='     '   ! SW 9/27/13 INITIALIZE ENTIRE ARRAY
  DO JF=1,73   ! Fix this later
    READ (CON,*)   KFNAME2(JF),(KFWBC(JF,JW),  JW=1,NWB)
  END DO
  KFWBC=ADJUSTR(KFWBC)

! Kinetics coefficients
  READ (CON,*)
  READ (CON,*) 

  READ (CON,*)     (EXH2O(JW),  JW=1,NWB)
  READ (CON,*)     (EXSS(JW),   JW=1,NWB)
  READ (CON,*)     (EXOM(JW),   JW=1,NWB)
  READ (CON,*)     (BETA(JW),   JW=1,NWB)
  READ (CON,*)     (EXC(JW),    JW=1,NWB);EXC=ADJUSTR(EXC)
  READ (CON,*)     (EXIC(JW),   JW=1,NWB);EXIC=ADJUSTR(EXIC)
 
  READ (CON,*)
  READ (CON,*) 

  READ (CON,*)         (EXA(JA),  JA=1,NAL)
  READ (CON,*)
  READ (CON,*) 

  READ (CON,*)         (EXZ(JZ),  JZ=1,NZPT)  
  READ (CON,*)
  READ (CON,*) 
  
  READ (CON,*)         (EXM(JM),   JM=1,NMCT)  
  READ (CON,*)
  READ (CON,*) 

  READ (CON,*)         (CGQ10(JG),  JG=1,NGC)
  READ (CON,*)         (CG0DK(JG),  JG=1,NGC)
  READ (CON,*)         (CG1DK(JG),  JG=1,NGC) 
  READ (CON,*)         (CGS(JG),    JG=1,NGC) 
  READ (CON,*)         (CGLDK(JG),  JG=1,NGC) 
  READ (CON,*)         (CGKLF(JG),  JG=1,NGC) 
  READ (CON,*)         (CGCS(JG),   JG=1,NGC) 
 
  READ (CON,*)
  READ (CON,*)  

  READ (CON,*)   (SSS(JS),    JS=1,NSS)     ! READ (CON,'(//(8X,F8.0,A8,2F8.0,I8))') (SSS(JS), SEDRC(JS),  TAUCR(JS),  SSFLOC(JS), FLOCEQN(JS),            JS=1,NSS) !SR 04/21/13
  READ (CON,*)   (SEDRC(JS),  JS=1,NSS)     ! READ (CON,'(//(8X,F8.0,A8,2F8.0,I8))') (SSS(JS), SEDRC(JS),  TAUCR(JS),  SSFLOC(JS), FLOCEQN(JS),            JS=1,NSS) !SR 04/21/13
  READ (CON,*)   (TAUCR(JS),  JS=1,NSS)     ! READ (CON,'(//(8X,F8.0,A8,2F8.0,I8))') (SSS(JS), SEDRC(JS),  TAUCR(JS),  SSFLOC(JS), FLOCEQN(JS),            JS=1,NSS) !SR 04/21/13
  SEDRC=ADJUSTR(SEDRC)
  READ (CON,*)
  READ (CON,*)  

  READ (CON,*) (AG(JA),          JA=1,NAL)
  READ (CON,*) (AR(JA),          JA=1,NAL)
  READ (CON,*) (AE(JA),          JA=1,NAL)
  READ (CON,*) (AM(JA),          JA=1,NAL)
  READ (CON,*) (AS(JA),          JA=1,NAL)
  READ (CON,*) (AHSP(JA),        JA=1,NAL)
  READ (CON,*) (AHSN(JA),        JA=1,NAL)
  READ (CON,*) (AHSSI(JA),       JA=1,NAL)
  READ (CON,*) (ASAT(JA),        JA=1,NAL)

  READ (CON,*) (AT1(JA),         JA=1,NAL)
  READ (CON,*) (AT2(JA),         JA=1,NAL)
  READ (CON,*) (AT3(JA),         JA=1,NAL)
  READ (CON,*) (AT4(JA),         JA=1,NAL)
  READ (CON,*) (AK1(JA),         JA=1,NAL)
  READ (CON,*) (AK2(JA),         JA=1,NAL)
  READ (CON,*) (AK3(JA),         JA=1,NAL)
  READ (CON,*) (AK4(JA),         JA=1,NAL)

  READ (CON,*) (AP(JA),          JA=1,NAL)
  READ (CON,*) (AN(JA),          JA=1,NAL)
  READ (CON,*) (AC(JA),          JA=1,NAL)
  READ (CON,*) (ASI(JA),         JA=1,NAL)
  READ (CON,*) (ACHLA(JA),       JA=1,NAL)
  READ (CON,*) (APOM(JA),        JA=1,NAL)
  READ (CON,*) (ANEQN(JA),       JA=1,NAL)
  READ (CON,*) (ANPR(JA),        JA=1,NAL)

  READ (CON,*) (O2AR(JA),        JA=1,NAL)
  READ (CON,*) (O2AG(JA),        JA=1,NAL)
  READ (CON,*)
  READ (CON,*)  

  DO JE=1,NEPTT
  READ (CON,*)         (EPIC(JW,JE),  JW=1,NWB)
  READ (CON,*)         (EPIPRC(JW,JE),JW=1,NWB)
  READ (CON,*)         (EPICI(JW,JE), JW=1,NWB)
  ENDDO
  EPIC=ADJUSTR(EPIC);EPIPRC=ADJUSTR(EPIPRC)

  READ (CON,*)
  READ (CON,*)  
  
  READ (CON,*) (EG(JE),           JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*) (ER(JE),           JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*) (EE(JE),           JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*) (EM(JE),           JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*) (EB(JE),           JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*) (EHSP(JE),         JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*) (EHSN(JE),         JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*) (EHSSI(JE),        JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13  
  
  READ (CON,*) (ESAT(JE),         JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*) (EHS(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*) (ENEQN(JE),        JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*) (ENPR(JE),         JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13  
  
  READ (CON,*) (ET1(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*) (ET2(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*) (ET3(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*) (ET4(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*) (EK1(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*) (EK2(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*) (EK3(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*) (EK4(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  
  READ (CON,*) (EP(JE),           JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*) (EN(JE),           JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*) (EC(JE),           JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*) (ESI(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*) (ECHLA(JE),        JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*) (EPOM(JE),         JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  READ (CON,*) (O2ER(JE),         JE=1,NEPT)
  READ (CON,*) (O2EG(JE),         JE=1,NEPT)

  READ (CON,*)
  READ (CON,*)  
   
  READ (CON,*)         (ZG(JZ),    JZ=1,NZPT)
  READ (CON,*)         (ZR(JZ),    JZ=1,NZPT)
  READ (CON,*)         (ZM(JZ),    JZ=1,NZPT)
  READ (CON,*)         (ZEFF(JZ),  JZ=1,NZPT)
  READ (CON,*)         (PREFP(JZ), JZ=1,NZPT)
  READ (CON,*)         (ZOOMIN(JZ),JZ=1,NZPT)
  READ (CON,*)         (ZS2P(JZ),  JZ=1,NZPT)
  READ (CON,*)         (ZT1(JZ),   JZ=1,NZPT)
  READ (CON,*)         (ZT2(JZ),   JZ=1,NZPT)
  READ (CON,*)         (ZT3(JZ),   JZ=1,NZPT)
  READ (CON,*)         (ZT4(JZ),   JZ=1,NZPT)
  READ (CON,*)         (ZK1(JZ),   JZ=1,NZPT)
  READ (CON,*)         (ZK2(JZ),   JZ=1,NZPT)
  READ (CON,*)         (ZK3(JZ),   JZ=1,NZPT)
  READ (CON,*)         (ZK4(JZ),   JZ=1,NZPT) 
  
  READ (CON,*)         (ZP(JZ),    JZ=1,NZPT)  
  READ (CON,*)         (ZN(JZ),    JZ=1,NZPT)
  READ (CON,*)         (ZC(JZ),    JZ=1,NZPT)

  READ (CON,*)         (O2ZR(JZ),  JZ=1,NZPT)

  DO JA=1,NALT
    READ (CON,*)    (PREFA(JA,JZ),     JZ=1,NZPT)
  END DO

  DO JZ=1,NZPTT
    READ (CON,*)    (PREFZ(JZ,JJZ),   JJZ=1,NZPT)       
  END DO
  
  READ (CON,*)
  READ (CON,*)   

  DO JM=1,NMCTT
    READ (CON,*)   (MACWBC(JW,JM),  JW=1,NWB)
  END DO
  DO JM=1,NMCTT
    READ (CON,*)   (MPRWBC(JW,JM),  JW=1,NWB)
  END DO
  DO JM=1,NMCTT
    READ (CON,*)   (MACWBCI(JW,JM),  JW=1,NWB)
  END DO
  MACWBC=ADJUSTR(MACWBC);MPRWBC=ADJUSTR(MPRWBC)
    
  READ (CON,*)
  READ (CON,*)

  
  READ (CON,*)         (MG(JM),     JM=1,NMCT)
  READ (CON,*)         (MR(JM),     JM=1,NMCT)
  READ (CON,*)         (MM(JM),     JM=1,NMCT)
  READ (CON,*)         (MSAT(JM),   JM=1,NMCT)
  READ (CON,*)         (MHSP(JM),   JM=1,NMCT)
  READ (CON,*)         (MHSN(JM),   JM=1,NMCT)
  READ (CON,*)         (MHSC(JM),   JM=1,NMCT)
  READ (CON,*)         (MPOM(JM),   JM=1,NMCT)
  READ (CON,*)         (LRPMAC(JM), JM=1,NMCT)
  
  READ (CON,*)         (PSED(JM),   JM=1,NMCT)  
  READ (CON,*)         (NSED(JM),   JM=1,NMCT)

  READ (CON,*)         (MBMP(JM),   JM=1,NMCT)
  READ (CON,*)         (MMAX(JM),   JM=1,NMCT)
  READ (CON,*)         (CDDRAG(JM), JM=1,NMCT)  !CB 6/29/06
  READ (CON,*)         (DWV(JM),    JM=1,NMCT)  !CB 6/29/06
  READ (CON,*)         (DWSA(JM),   JM=1,NMCT)  !CB 6/29/06
  READ (CON,*)         (ANORM(JM),  JM=1,NMCT)  !CB 6/29/06  

  READ (CON,*)         (MT1(JM),    JM=1,NMCT)
  READ (CON,*)         (MT2(JM),    JM=1,NMCT)
  READ (CON,*)         (MT3(JM),    JM=1,NMCT)
  READ (CON,*)         (MT4(JM),    JM=1,NMCT)
  READ (CON,*)         (MK1(JM),    JM=1,NMCT)
  READ (CON,*)         (MK2(JM),    JM=1,NMCT)
  READ (CON,*)         (MK3(JM),    JM=1,NMCT)
  READ (CON,*)         (MK4(JM),    JM=1,NMCT)
  
  READ (CON,*)         (MP(JM),     JM=1,NMCT)
  READ (CON,*)         (MN(JM),     JM=1,NMCT)
  READ (CON,*)         (MC(JM),     JM=1,NMCT)
 
  READ (CON,*)         (O2MR(JM),   JM=1,NMCT)
  READ (CON,*)         (O2MG(JM),   JM=1,NMCT)

  READ (CON,*)
  READ (CON,*)   

  READ (CON,*)         (LDOMDK(JW),  JW=1,NWB)
  READ (CON,*)         (RDOMDK(JW),  JW=1,NWB)
  READ (CON,*)         (LRDDK(JW),   JW=1,NWB) 

  READ (CON,*)
  READ (CON,*)   
   
  READ (CON,*)         (LPOMDK(JW),   JW=1,NWB)
  READ (CON,*)         (RPOMDK(JW),   JW=1,NWB)
  READ (CON,*)         (LRPDK(JW),    JW=1,NWB)
  READ (CON,*)         (POMS(JW),     JW=1,NWB)
  
  READ (CON,*)
  READ (CON,*)   
  
  READ (CON,*)         (ORGP(JW),     JW=1,NWB)
  READ (CON,*)         (ORGN(JW),     JW=1,NWB)
  READ (CON,*)         (ORGC(JW),     JW=1,NWB)
  READ (CON,*)         (ORGSI(JW),    JW=1,NWB)
  READ (CON,*)         (O2OM(JW),     JW=1,NWB)

  READ (CON,*)         (OMT1(JW),     JW=1,NWB)
  READ (CON,*)         (OMT2(JW),     JW=1,NWB)
  READ (CON,*)         (OMK1(JW),     JW=1,NWB)
  READ (CON,*)         (OMK2(JW),     JW=1,NWB)

  READ (CON,*)
  READ (CON,*)   

  READ (CON,*)         (KBOD(JB),     JB=1,NBOD)
  READ (CON,*)         (TBOD(JB),     JB=1,NBOD)
  READ (CON,*)         (RBOD(JB),     JB=1,NBOD)
  READ (CON,*)         (CBODS(JB),    JB=1,NBOD)
  READ (CON,*)         (BODP(JB),     JB=1,NBOD)
  READ (CON,*)         (BODN(JB),     JB=1,NBOD)
  READ (CON,*)         (BODC(JB),     JB=1,NBOD)
  
  READ (CON,*)
  READ (CON,*)   

  READ (CON,*)         (PO4R(JW),     JW=1,NWB)
  READ (CON,*)         (PARTP(JW),    JW=1,NWB)
  READ (CON,*)         (NH4R(JW),     JW=1,NWB)
  READ (CON,*)         (NH4DK(JW),    JW=1,NWB)
  READ (CON,*)         (NH4T1(JW),    JW=1,NWB)
  READ (CON,*)         (NH4T2(JW),    JW=1,NWB)
  READ (CON,*)         (NH4K1(JW),    JW=1,NWB)
  READ (CON,*)         (NH4K2(JW),    JW=1,NWB)
  READ (CON,*)         (O2NH4(JW),    JW=1,NWB)
  READ (CON,*)         (NO3DK(JW),    JW=1,NWB)
  READ (CON,*)         (NO3S(JW),     JW=1,NWB)
  READ (CON,*)         (FNO3SED(JW),  JW=1,NWB)
  READ (CON,*)         (NO3T1(JW),    JW=1,NWB)
  READ (CON,*)         (NO3T2(JW),    JW=1,NWB)
  READ (CON,*)         (NO3K1(JW),    JW=1,NWB)
  READ (CON,*)         (NO3K2(JW),    JW=1,NWB)
  READ (CON,*)         (DSIR(JW),     JW=1,NWB)
  READ (CON,*)         (PSIS(JW),     JW=1,NWB)
  READ (CON,*)         (PSIDK(JW),    JW=1,NWB)
  READ (CON,*)         (PARTSI(JW),   JW=1,NWB)
    
  READ (CON,*)
  READ (CON,*)   

  READ (CON,*)         (FER(JW),       JW=1,NWB)
  READ (CON,*)         (FES(JW),       JW=1,NWB)
     
  READ (CON,*)
  READ (CON,*)   

  READ (CON,*)          (CO2R(JW),     JW=1,NWB)
  READ (CON,*)
  READ (CON,*)   
  
  READ (CON,*)          KDO

  READ (CON,*)
  READ (CON,*)   

  READ (CON,*)     (SEDCC(JW),   JW=1,NWB); SEDCC=adjustr(SEDCC)  
  READ (CON,*)     (SEDPRC(JW),  JW=1,NWB); SEDPRC=adjustr(SEDPRC) 
  READ (CON,*)     (SEDCI(JW),   JW=1,NWB)  
  READ (CON,*)     (SDK(JW),   JW=1,NWB)     
  READ (CON,*)     (SEDS(JW),    JW=1,NWB)  
  READ (CON,*)     (FSOD(JW),    JW=1,NWB)  
  READ (CON,*)     (FSED(JW),    JW=1,NWB) 
  READ (CON,*)     (SEDB(JW),   JW=1,NWB)  
  READ (CON,*)     (DYNSEDK(JW), JW=1,NWB); DYNSEDK=adjustr(DYNSEDK)
  READ (CON,*)     (SODT1(JW),   JW=1,NWB)
  READ (CON,*)     (SODT2(JW),   JW=1,NWB)
  READ (CON,*)     (SODK1(JW),   JW=1,NWB)
  READ (CON,*)     (SODK2(JW),   JW=1,NWB)
  READ (CON,*)
  READ (CON,*)   

  READ (CON,*)  (SOD(I),  I=1,IMX)
  
  READ (CON,*)
  READ (CON,*)    

  READ (CON,*)   (REAERC(JW), JW=1,NWB); REAERC=adjustr(REAERC)
  READ (CON,*)   (NEQN(JW),   JW=1,NWB)
  READ (CON,*)   (RCOEF1(JW), JW=1,NWB)
  READ (CON,*)   (RCOEF2(JW), JW=1,NWB)
  READ (CON,*)   (RCOEF3(JW), JW=1,NWB)
  READ (CON,*)   (RCOEF4(JW), JW=1,NWB)
  READ (CON,*)
  READ (CON,*)    
  
! Input filenames

  READ (CON,*)  QWDFN
  READ (CON,*)  QGTFN
  READ (CON,*)  WSCFN
  READ (CON,*)  SHDFN
  READ (CON,*)  VPLFN(1)
  VPLFN(2:NWB)=VPLFN(1)
  
  READ (CON,*)
  READ (CON,*)    
 
  READ (CON,*) (BTHFN(JW), JW=1,NWB)
  READ (CON,*) (METFN(JW), JW=1,NWB)
  READ (CON,*) (EXTFN(JW), JW=1,NWB)
  READ (CON,*) (VPRFN(JW), JW=1,NWB)
  READ (CON,*) (LPRFN(JW), JW=1,NWB)

  ! Output filenames

  READ (CON,*) (SNPFN(JW), JW=1,NWB)
  READ (CON,*) (PRFFN(JW), JW=1,NWB)
  READ (CON,*) (CPLFN(JW), JW=1,NWB)
  READ (CON,*) (SPRFN(JW), JW=1,NWB)
  READ (CON,*) (FLXFN(JW), JW=1,NWB)
  
  READ (CON,*)
  READ (CON,*)    
 
  READ (CON,*) (QINFN(JB), JB=1,NBR)
  READ (CON,*) (TINFN(JB), JB=1,NBR)
  READ (CON,*) (CINFN(JB), JB=1,NBR)
  READ (CON,*) (QOTFN(JB), JB=1,NBR)
  READ (CON,*) (QDTFN(JB), JB=1,NBR)
  READ (CON,*) (TDTFN(JB), JB=1,NBR)
  READ (CON,*) (CDTFN(JB), JB=1,NBR)
  READ (CON,*) (PREFN(JB), JB=1,NBR)
  READ (CON,*) (TPRFN(JB), JB=1,NBR)
  READ (CON,*) (CPRFN(JB), JB=1,NBR)
  READ (CON,*) (EUHFN(JB), JB=1,NBR)
  READ (CON,*) (TUHFN(JB), JB=1,NBR)
  READ (CON,*) (CUHFN(JB), JB=1,NBR)
  READ (CON,*) (EDHFN(JB), JB=1,NBR)
  READ (CON,*) (TDHFN(JB), JB=1,NBR)
  READ (CON,*) (CDHFN(JB), JB=1,NBR)
 
  CLOSE (CON) 
RETURN

    END SUBROUTINE READ_CONTROL_FILECSV


!**************************************************************
  subroutine read_control_fileNPT
	USE MAIN
    USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  use macrophytec; use porosityc; use zooplanktonc  
  use initialvelocity
  use extra  
  IMPLICIT NONE
  
  character*1 char1,ichar1
  character*8 AID
  character*72 linetest
  CHARACTER(8)  :: IBLANK
  character(2)   :: ichar2
  INTEGER :: NPT, NNBP, NCBP,NINTERNAL,NUP,KTMAX

  !open(CON,file=confn,status='old')

! Title and array dimensions

  ALLOCATE (TITLE(11))
  READ (CON,'(///(8X,A72))') (TITLE(J),J=1,10)
  READ (CON,'(//8X,5I8,2A8)') NWB, NBR, IMX, KMX, NPROC, CLOSEC                     ! SW 7/31/09
  READ (CON,'(//8X,8I8)')     NTR, NST, NIW, NWD, NGT, NSP, NPI, NPU
  READ (CON,'(//8X,7I8,a8)')  NGC, NSS, NAL, NEP, NBOD, nmc, nzp  
  READ (CON,'(//8X,I8,8A8)')  NOD,SELECTC,HABTATC,ENVIRPC,AERATEC,inituwl,SYSTDG,N2BND,DOBND

  if(NPROC == 0)NPROC=1                                                                 ! SW 7/31/09
  !call omp_set_num_threads(NPROC)   ! set # of processors to NPROC  Moved to INPUT subroutine  TOGGLE FOR DEBUG
  if(SELECTC=='        ')then
     SELECTC='     OFF'
  endif
  
! Constituent numbers

  NTDS  = 1
  NGCS  = 2
  NGCE  = NGCS+NGC-1
  NSSS  = NGCE+1
  NSSE  = NSSS+NSS-1
  NPO4  = NSSE+1
  NNH4  = NPO4+1
  NNO3  = NNH4+1
  NDSI  = NNO3+1
  NPSI  = NDSI+1
  NFE   = NPSI+1
  NLDOM = NFE+1
  NRDOM = NLDOM+1
  NLPOM = NRDOM+1
  NRPOM = NLPOM+1
  NBODS = NRPOM+1
if(nbod.gt.0)then    ! variable stoichiometry for CBOD    ! cb 6/6/10
 ALLOCATE (nbodc(nbod), nbodp(nbod), nbodn(nbod))
  ibod=nbods
  nbodcs=ibod  
  do jcb=1,nbod
     nbodc(jcb)=ibod
     ibod=ibod+1    
  end do
  nbodce=ibod-1
  nbodps=ibod
  do jcb=1,nbod
    nbodp(jcb)=ibod
    ibod=ibod+1
  end do
  nbodpe=ibod-1
  nbodns=ibod
  do jcb=1,nbod
    nbodn(jcb)=ibod
    ibod=ibod+1
  end do
  nbodne=ibod-1
  ELSE
    NBODNS=1;NBODNE=1;NBODPS=1;NBODPE=1;NBODCS=1;NBODCE=1
end if
!    NBODE = NBODS+NBOD-1 
  NBODE = NBODS+NBOD*3-1
  NAS   = NBODE+1
  NAE   = NAS+NAL-1
  NDO   = NAE+1
  NTIC  = NDO+1
  NALK  = NTIC+1

! v3.5 start
  NZOOS = NALK + 1
  NZOOE = NZOOS + NZP - 1
  NLDOMP=NZOOE+1
  NRDOMP=nldomp+1
  NLPOMP=nrdomp+1
  NRPOMP=nlpomp+1
  NLDOMn=nrpomp+1
  NRDOMn=nldomn+1
  NLPOMn=nrdomn+1
  NRPOMn=nlpomn+1
  nct=nrpomn
! v3.5 end

! Constituent, tributary, and widthdrawal totals

  NTRT = MAX(NTR+NGT+NSP+NPI+NPU,1)
  NWDT = MAX(NWD+NGT+NSP+NPI+NPU,1)
  NEPT = MAX(NEP,1)
  Nmct = MAX(nmc,1)
  nzpt=max(nzp,1)

  ALLOCATE (CDAC(NDC), X1(IMX), TECPLOT(NWB))
  ALLOCATE (WSC(IMX),    KBI(IMX))
  ALLOCATE (VBC(NWB),    EBC(NWB),    MBC(NWB),    PQC(NWB),    EVC(NWB),    PRC(NWB))
  ALLOCATE (WINDC(NWB),  QINC(NWB),   QOUTC(NWB),  HEATC(NWB),  SLHTC(NWB))
  ALLOCATE (QINIC(NBR),  DTRIC(NBR),  TRIC(NTR),   WDIC(NWD),   HDIC(NBR),   METIC(NWB))
  ALLOCATE (EXC(NWB),    EXIC(NWB))
  ALLOCATE (SLTRC(NWB),  THETA(NWB),  FRICC(NWB),  NAF(NWB),    ELTMF(NWB), Z0(NWB))
  ALLOCATE (ZMIN(NWB),   IZMIN(NWB))
  ALLOCATE (C2CH(NCT),   CDCH(NDC),   EPCH(NEPT),  macch(nmct), KFCH(NFL))  !v3.5
  ALLOCATE (CPLTC(NCT),  HPLTC(NHY),  CDPLTC(NDC))
  ALLOCATE (CMIN(NCT),   CMAX(NCT),   HYMIN(NHY),  HYMAX(NHY),  CDMIN(NDC),  CDMAX(NDC))
  ALLOCATE (JBDAM(NBR),  ILAT(NWDT))
  ALLOCATE (QINSUM(NBR), TINSUM(NBR), TIND(NBR),   JSS(NBR),    QIND(NBR))
  ALLOCATE (QOLD(NPI),   DTP(NPI),    DTPS(NPI),   QOLDS(NPI))
  ALLOCATE (LATGTC(NGT), LATSPC(NSP), LATPIC(NPI), DYNPIPE(NPI),DYNPUMP(NPU),LATPUC(NPU), DYNGTC(NGT))
  ALLOCATE (gtic(NGT),BGTo(NGT),   EGTo(NGT) )   ! cb 8/13/2010
  ALLOCATE (INTERP_gate(ngt))                     ! cb 8/13/2010  
  ALLOCATE (OPT(NWB,7))
  ALLOCATE (CDWBC(NDC,NWB),     KFWBC(NFL,NWB),        CPRWBC(NCT,NWB),    CINBRC(NCT,NBR),     CTRTRC(NCT,NTRT))
  ALLOCATE (CDTBRC(NCT,NBR),    CPRBRC(NCT,NBR))
  ALLOCATE (YSS(NNPIPE,NPI),    VSS(NNPIPE,NPI),       YS(NNPIPE,NPI),     VS(NNPIPE,NPI),      VSTS(NNPIPE,NPI))
  ALLOCATE (YSTS(NNPIPE,NPI),   YST(NNPIPE,NPI),       VST(NNPIPE,NPI))
  ALLOCATE (CBODD(KMX,IMX,NBOD))
  ALLOCATE (ALLIM(KMX,IMX,NAL), APLIM(KMX,IMX,NAL),    ANLIM(KMX,IMX,NAL), ASLIM(KMX,IMX,NAL))
  ALLOCATE (ELLIM(KMX,IMX,NEP), EPLIM(KMX,IMX,NEP),    ENLIM(KMX,IMX,NEP), ESLIM(KMX,IMX,NEP))
  ALLOCATE (CSSK(KMX,IMX,NCT),  C1(KMX,IMX,NCT),       C2(KMX,IMX,NCT),    CD(KMX,IMX,NDC),     KF(KMX,IMX,NFL))
  ALLOCATE (KFS(KMX,IMX,NFL),   AF(KMX,IMX,NAL,5),     EF(KMX,IMX,NEP,5),  HYD(KMX,IMX,NHY), KFJW(NWB,NFL))
  ALLOCATE (TKE(KMX,IMX,3), AZT(KMX,IMX),DZT(KMX,IMX))
  ALLOCATE (USTARBTKE(IMX),E(IMX),EROUGH(NWB), ARODI(NWB), STRICK(NWB), TKELATPRDCONST(NWB))
  ALLOCATE (FIRSTI(NWB), LASTI(NWB), TKELATPRD(NWB), STRICKON(NWB), WALLPNT(NWB), IMPTKE(NWB), TKEBC(NWB))
  ALLOCATE (HYDRO_PLOT(NHY),    CONSTITUENT_PLOT(NCT), DERIVED_PLOT(NDC))
  ALLOCATE (ZERO_SLOPE(NWB),    DYNAMIC_SHADE(IMX))
  ALLOCATE (AZSLC(NWB))
  ALLOCATE (NSPRF(NWB))
  ALLOCATE (KBMAX(NWB),  ELKT(NWB),   WIND2(IMX))
  ALLOCATE (VISC(NWB),   CELC(NWB),   REAERC(NWB))
  ALLOCATE (QOAVR(NWB),  QIMXR(NWB),  QOMXR(NWB))
  ALLOCATE (LAT(NWB),    LONGIT(NWB), ELBOT(NWB))
  ALLOCATE (BTH(NWB),    VPR(NWB),    LPR(NWB))
  ALLOCATE (NISNP(NWB),  NIPRF(NWB),  NISPR(NWB))
  ALLOCATE (icpl(NWB))                                 ! cb 1/26/09
  ALLOCATE (A00(NWB),    HH(NWB),     DECL(NWB))
  ALLOCATE (T2I(NWB),    KTWB(NWB),   KBR(NWB),    IBPR(NWB))
  ALLOCATE (DLVR(NWB),   ESR(NWB),    ETR(NWB),    NBL(NWB))
  ALLOCATE (LPRFN(NWB),  EXTFN(NWB),  BTHFN(NWB),  METFN(NWB),  VPRFN(NWB))
  ALLOCATE (SNPFN(NWB),  PRFFN(NWB),  SPRFN(NWB),  CPLFN(NWB),  VPLFN(NWB),  FLXFN(NWB),FLXFN2(NWB))
  ALLOCATE (AFW(NWB),    BFW(NWB),    CFW(NWB),    WINDH(NWB),  RHEVC(NWB),  FETCHC(NWB))
  ALLOCATE (SDK(NWB),    FSOD(NWB),   FSED(NWB),   SEDCI(NWB),  SEDCc(NWB),   SEDPRC(NWB), seds(nwb), sedb(nwb), DYNSEDK(NWB))  !cb 11/28/06
  ALLOCATE (ICEC(NWB),   SLICEC(NWB), ICETHI(NWB), ALBEDO(NWB), HWI(NWB),    BETAI(NWB),  GAMMAI(NWB), ICEMIN(NWB), ICET2(NWB))
  ALLOCATE (EXH2O(NWB),  BETA(NWB),   EXOM(NWB),   EXSS(NWB),   DXI(NWB),    CBHE(NWB),   TSED(NWB),   TSEDF(NWB),  FI(NWB))
  ALLOCATE (AX(NWB),     WTYPEC(NWB), JBDN(NWB),   AZC(NWB),    AZMAX(NWB),  QINT(NWB),   QOUTT(NWB),  GRIDC(NWB))     !SW 07/14/04
  ALLOCATE (TAIR(NWB),   TDEW(NWB),   WIND(NWB),   PHI(NWB),    CLOUD(NWB),  CSHE(IMX),   SRON(NWB),   RANLW(NWB))
  ALLOCATE (SNPC(NWB),   SCRC(NWB),   PRFC(NWB),   SPRC(NWB),   CPLC(NWB),   VPLC(NWB),   FLXC(NWB))
  ALLOCATE (NXTMSN(NWB), NXTMSC(NWB), NXTMPR(NWB), NXTMSP(NWB), NXTMCP(NWB), NXTMVP(NWB), NXTMFL(NWB))
  ALLOCATE (SNPDP(NWB),  SCRDP(NWB),  PRFDP(NWB),  SPRDP(NWB),  CPLDP(NWB),  VPLDP(NWB),  FLXDP(NWB))
  ALLOCATE (NSNP(NWB),   NSCR(NWB),   NPRF(NWB),   NSPR(NWB),   NCPL(NWB),   NVPL(NWB),   NFLX(NWB))
  ALLOCATE (NEQN(NWB),   PO4R(NWB),   PARTP(NWB))
  ALLOCATE (NH4DK(NWB),  NH4R(NWB))
  ALLOCATE (NO3DK(NWB),  NO3S(NWB), FNO3SED(NWB))
  ALLOCATE (FER(NWB),    FES(NWB))
  ALLOCATE (CO2R(NWB),   SROC(NWB))
  ALLOCATE (O2ER(NEPT),  O2EG(NEPT))
  ALLOCATE (CAQ10(NWB),  CADK(NWB),   CAS(NWB))
  ALLOCATE (BODP(NBOD),  BODN(NBOD),  BODC(NBOD))
  ALLOCATE (KBOD(NBOD),  TBOD(NBOD),  RBOD(NBOD))
  ALLOCATE (LDOMDK(NWB), RDOMDK(NWB), LRDDK(NWB))
  ALLOCATE (OMT1(NWB),   OMT2(NWB),   OMK1(NWB),   OMK2(NWB))
  ALLOCATE (LPOMDK(NWB), RPOMDK(NWB), LRPDK(NWB),  POMS(NWB))
  ALLOCATE (ORGP(NWB),   ORGN(NWB),   ORGC(NWB),   ORGSI(NWB))
  ALLOCATE (RCOEF1(NWB), RCOEF2(NWB), RCOEF3(NWB), RCOEF4(NWB))
  ALLOCATE (NH4T1(NWB),  NH4T2(NWB),  NH4K1(NWB),  NH4K2(NWB))
  ALLOCATE (NO3T1(NWB),  NO3T2(NWB),  NO3K1(NWB),  NO3K2(NWB))
  ALLOCATE (DSIR(NWB),   PSIS(NWB),   PSIDK(NWB),  PARTSI(NWB))
  ALLOCATE (SODT1(NWB),  SODT2(NWB),  SODK1(NWB),  SODK2(NWB))
  ALLOCATE (O2NH4(NWB),  O2OM(NWB))
  ALLOCATE (O2AR(NAL),   O2AG(NAL))
  ALLOCATE (CGQ10(NGC),  CG0DK(NGC),  CG1DK(NGC),  CGS(NGC), CGLDK(NGC),CGKLF(NGC),CGCS(NGC))
  ALLOCATE (CUNIT(NCT),  CUNIT1(NCT), CUNIT2(NCT))
  ALLOCATE (CAC(NCT),    INCAC(NCT),  TRCAC(NCT),  DTCAC(NCT),  PRCAC(NCT))
  ALLOCATE (CNAME(NCT),  CNAME1(NCT), CNAME2(NCT), CNAME3(NCT), CMULT(NCT),  CSUM(NCT))
  ALLOCATE (CN(NCT))
  ALLOCATE (SSS(NSS),    TAUCR(NSS),  SEDRC(NSS))
  ALLOCATE (CDSUM(NDC))
  ALLOCATE (DTRC(NBR))
  ALLOCATE (NSTR(NBR),   XBR(NBR),DYNSTRUC(NBR))
  ALLOCATE (QTAVB(NBR),  QTMXB(NBR))
  ALLOCATE (BS(NWB),     BE(NWB),     JBUH(NBR),   JBDH(NBR),   JWUH(NBR),   JWDH(NBR))
  ALLOCATE (TSSS(NBR),   TSSB(NBR),   TSSICE(NBR))
  ALLOCATE (ESBR(NBR),   ETBR(NBR),   EBRI(NBR))
  ALLOCATE (QIN(NBR),    PR(NBR),     QPRBR(NBR),  QDTR(NBR),   EVBR(NBR))
  ALLOCATE (TIN(NBR),    TOUT(NBR),   TPR(NBR),    TDTR(NBR),   TPB(NBR))
  ALLOCATE (NACPR(NBR),  NACIN(NBR),  NACDT(NBR),  NACTR(NTR),  NACD(NWB))
  ALLOCATE (QSUM(NBR),   NOUT(NBR),   KTQIN(NBR),  KBQIN(NBR),  ELUH(NBR),   ELDH(NBR))
  ALLOCATE (NL(NBR),     NPOINT(NBR), SLOPE(NBR),  SLOPEC(NBR), ALPHA(NBR),  COSA(NBR),   SINA(NBR),   SINAC(NBR), ilayer(imx))   
  ALLOCATE (CPRFN(NBR),  EUHFN(NBR),  TUHFN(NBR),  CUHFN(NBR),  EDHFN(NBR),  TDHFN(NBR),  QOTFN(NBR),  PREFN(NBR))
  ALLOCATE (QINFN(NBR),  TINFN(NBR),  CINFN(NBR),  CDHFN(NBR),  QDTFN(NBR),  TDTFN(NBR),  CDTFN(NBR),  TPRFN(NBR))
  ALLOCATE (VOLWD(NBR),  VOLSBR(NBR), VOLTBR(NBR), DLVOL(NBR),  VOLG(NWB),   VOLSR(NWB),  VOLTR(NWB),  VOLEV(NBR))
  ALLOCATE (VOLB(NBR),   VOLPR(NBR),  VOLTRB(NBR), VOLDT(NBR),  VOLUH(NBR),  VOLDH(NBR),  VOLIN(NBR),  VOLOUT(NBR))
  ALLOCATE (US(NBR),     DS(NBR),     CUS(NBR),    UHS(NBR),    DHS(NBR),    UQB(NBR),    DQB(NBR),    CDHS(NBR))
  ALLOCATE (TSSEV(NBR),  TSSPR(NBR),  TSSTR(NBR),  TSSDT(NBR),  TSSWD(NBR),  TSSUH(NBR),  TSSDH(NBR),  TSSIN(NBR),  TSSOUT(NBR))
  ALLOCATE (ET(IMX),     RS(IMX),     RN(IMX),     RB(IMX),     RC(IMX),     RE(IMX),     SHADE(IMX))
  ALLOCATE (DLTMAX(NOD), QWDO(IMX),   TWDO(IMX))                                                                        ! SW 1/24/05
  ALLOCATE (SOD(IMX),    ELWS(IMX),   BKT(IMX),    REAER(IMX))
  ALLOCATE (ICETH(IMX),  ICE(IMX),    ICESW(IMX))
  ALLOCATE (Q(IMX),      QC(IMX),     QERR(IMX),   QSSUM(IMX))
  ALLOCATE (KTI(IMX),    SKTI(IMX),   SROSH(IMX),  SEG(IMX),    DLXRHO(IMX))
  ALLOCATE (DLX(IMX),    DLXR(IMX))
  ALLOCATE (A(IMX),      C(IMX),      D(IMX),      F(IMX),      V(IMX),      BTA(IMX),    GMA(IMX))
  ALLOCATE (KBMIN(IMX),  EV(IMX),     QDT(IMX),    QPR(IMX),    SBKT(IMX),   BHRHO(IMX))
  ALLOCATE (SZ(IMX),     WSHX(IMX),   WSHY(IMX),   WIND10(IMX), CZ(IMX),     FETCH(IMX),  PHI0(IMX),   FRIC(IMX))
  ALLOCATE (Z(IMX),      KB(IMX),     PALT(IMX))
  ALLOCATE (VNORM(KMX))
  ALLOCATE (ANPR(NAL),   ANEQN(NAL),  APOM(NAL))
  ALLOCATE (AC(NAL),     ASI(NAL),    ACHLA(NAL),  AHSP(NAL),   AHSN(NAL),   AHSSI(NAL))
  ALLOCATE (AT1(NAL),    AT2(NAL),    AT3(NAL),    AT4(NAL),    AK1(NAL),    AK2(NAL),    AK3(NAL),    AK4(NAL))
  ALLOCATE (AG(NAL),     AR(NAL),     AE(NAL),     AM(NAL),     AS(NAL),     EXA(NAL),    ASAT(NAL),   AP(NAL),   AN(NAL))
  ALLOCATE (ENPR(NEPT),  ENEQN(NEPT))
  ALLOCATE (EG(NEPT),    ER(NEPT),    EE(NEPT),    EM(NEPT),    EB(NEPT),    ESAT(NEPT),  EP(NEPT),    EN(NEPT))
  ALLOCATE (EC(NEPT),    ESI(NEPT),   ECHLA(NEPT), EHSP(NEPT),  EHSN(NEPT),  EHSSI(NEPT), EPOM(NEPT),  EHS(NEPT))
  ALLOCATE (ET1(NEPT),   ET2(NEPT),   ET3(NEPT),   ET4(NEPT),   EK1(NEPT),   EK2(NEPT),   EK3(NEPT),   EK4(NEPT))
  ALLOCATE (HNAME(NHY),  FMTH(NHY),   HMULT(NHY),  FMTC(NCT),   FMTCD(NDC))
  ALLOCATE (KFAC(NFL),   KFNAME(NFL), KFNAME2(NFL),KFCN(NFL,NWB))
  ALLOCATE (C2I(NCT,NWB),    TRCN(NCT,NTR))
  ALLOCATE (CDN(NDC,NWB),    CDNAME(NDC),     CDNAME2(NDC),    CDNAME3(NDC),    CDMULT(NDC))
  ALLOCATE (CMBRS(NCT,NBR),  CMBRT(NCT,NBR),  INCN(NCT,NBR),   DTCN(NCT,NBR),   PRCN(NCT,NBR))
  ALLOCATE (FETCHU(IMX,NBR), FETCHD(IMX,NBR))
  ALLOCATE (IPRF(IMX,NWB),   ISNP(IMX,NWB),   ISPR(IMX,NWB),   BL(IMX,NWB))
  ALLOCATE (H1(KMX,IMX),     H2(KMX,IMX),     BH1(KMX,IMX),    BH2(KMX,IMX),    BHR1(KMX,IMX),   BHR2(KMX,IMX),   QTOT(KMX,IMX))
  ALLOCATE (SAVH2(KMX,IMX),  AVH1(KMX,IMX),   AVH2(KMX,IMX),   AVHR(KMX,IMX),   SAVHR(KMX,IMX))
  ALLOCATE (LFPR(KMX,IMX),   BI(KMX,IMX), bnew(kmx,imx))        ! SW 1/23/06
  ALLOCATE (ADX(KMX,IMX),    ADZ(KMX,IMX),    DO1(KMX,IMX),    DO2(KMX,IMX),    DO3(KMX,IMX),    SED(KMX,IMX))
  ALLOCATE (B(KMX,IMX),      CONV(KMX,IMX),   CONV1(KMX,IMX),  EL(KMX,IMX),     DZ(KMX,IMX),     DZQ(KMX,IMX),    DX(KMX,IMX))
  ALLOCATE (P(KMX,IMX),      SU(KMX,IMX),     SW(KMX,IMX),     SAZ(KMX,IMX),    T1(KMX,IMX),     TSS(KMX,IMX),    QSS(KMX,IMX))
  ALLOCATE (BB(KMX,IMX),     BR(KMX,IMX),     BH(KMX,IMX),     BHR(KMX,IMX),    VOL(KMX,IMX),    HSEG(KMX,IMX),   DECAY(KMX,IMX))
  ALLOCATE (DEPTHB(KMX,IMX), DEPTHM(KMX,IMX), FPSS(KMX,IMX),   FPFE(KMX,IMX),   FRICBR(KMX,IMX), UXBR(KMX,IMX),   UYBR(KMX,IMX))
  ALLOCATE (QUH1(KMX,NBR),   QDH1(KMX,NBR),   VOLUH2(KMX,NBR), VOLDH2(KMX,NBR))
  ALLOCATE (TSSUH1(KMX,NBR), TSSUH2(KMX,NBR), TSSDH1(KMX,NBR), TSSDH2(KMX,NBR))
  ALLOCATE (TVP(KMX,NWB),    SEDVP(KMX,NWB),  H(KMX,NWB))
  ALLOCATE (QINF(KMX,NBR),   KOUT(KMX,NBR))
  ALLOCATE (CT(KMX,IMX),     AT(KMX,IMX),     VT(KMX,IMX),     DT(KMX,IMX),     GAMMA(KMX,IMX))
  ALLOCATE (CWDO(NCT,NOD),   CDWDO(NDC,NOD),  CWDOC(NCT),      CDWDOC(NDC),     CDTOT(NDC))
  ALLOCATE (RSOD(NOD),       RSOF(NOD),       DLTD(NOD),       DLTF(NOD))
  ALLOCATE (TSRD(NOD),       TSRF(NOD),       WDOD(NOD),       WDOF(NOD))
  ALLOCATE (SNPD(NOD,NWB),   SNPF(NOD,NWB),   SPRD(NOD,NWB),   SPRF(NOD,NWB))
  ALLOCATE (SCRD(NOD,NWB),   SCRF(NOD,NWB),   PRFD(NOD,NWB),   PRFF(NOD,NWB))
  ALLOCATE (CPLD(NOD,NWB),   CPLF(NOD,NWB),   VPLD(NOD,NWB),   VPLF(NOD,NWB),   FLXD(NOD,NWB),   FLXF(NOD,NWB))
  ALLOCATE (EPIC(NWB,NEPT),  EPICI(NWB,NEPT), EPIPRC(NWB,NEPT))
  ALLOCATE (EPIVP(KMX,NWB,NEP))
  ALLOCATE (CUH(KMX,NCT,NBR),     CDH(KMX,NCT,NBR))
  ALLOCATE (EPM(KMX,IMX,NEPT),    EPD(KMX,IMX,NEPT),    EPC(KMX,IMX,NEPT))
  ALLOCATE (C1S(KMX,IMX,NCT),     CSSB(KMX,IMX,NCT),    CVP(KMX,NCT,NWB))
  ALLOCATE (CSSUH1(KMX,NCT,NBR),  CSSUH2(KMX,NCT,NBR),  CSSDH2(KMX,NCT,NBR), CSSDH1(KMX,NCT,NBR))
  ALLOCATE (OPEN_VPR(NWB),        OPEN_LPR(NWB))
  ALLOCATE (READ_EXTINCTION(NWB), READ_RADIATION(NWB))
  ALLOCATE (DIST_TRIBS(NBR),      LIMITING_FACTOR(NAL))
  ALLOCATE (UPWIND(NWB),          ULTIMATE(NWB))
  NSTT=MAX(NST,5)
  ALLOCATE (STRIC(NSTT,NBR), ESTRT(NSTT,NBR), WSTRT(NSTT,NBR), KTSWT(NSTT,NBR), KBSWT(NSTT,NBR),SINKCT(NSTT,NBR))
  allocate (estr(nstT,nbr))
  ALLOCATE (FRESH_WATER(NWB),     SALT_WATER(NWB),      TRAPEZOIDAL(NWB))                                              !SW 07/16/04
  ALLOCATE (UH_EXTERNAL(NBR),     DH_EXTERNAL(NBR),     UH_INTERNAL(NBR),    DH_INTERNAL(NBR))
  ALLOCATE (UQ_EXTERNAL(NBR),     DQ_EXTERNAL(NBR),     UQ_INTERNAL(NBR),    DQ_INTERNAL(NBR))
  ALLOCATE (UP_FLOW(NBR),         DN_FLOW(NBR),         UP_HEAD(NBR),        DN_HEAD(NBR))
  ALLOCATE (INTERNAL_FLOW(NBR),   DAM_INFLOW(NBR),      DAM_OUTFLOW(NBR),    HEAD_FLOW(NBR),      HEAD_BOUNDARY(NWB))  !TC 08/03/04
  ALLOCATE (ISO_CONC(NCT,NWB),    VERT_CONC(NCT,NWB),   LONG_CONC(NCT,NWB))
  ALLOCATE (ISO_SEDIMENT(NWB),    VERT_SEDIMENT(NWB),   LONG_SEDIMENT(NWB))
  ALLOCATE (VISCOSITY_LIMIT(NWB), CELERITY_LIMIT(NWB),  IMPLICIT_AZ(NWB))
  ALLOCATE (FETCH_CALC(NWB),      ONE_LAYER(IMX),       IMPLICIT_VISC(NWB))
  ALLOCATE (LIMITING_DLT(NWB),    TERM_BY_TERM(NWB),    MANNINGS_N(NWB))
  ALLOCATE (PLACE_QIN(NWB),       PLACE_QTR(NTRT),      SPECIFY_QTR(NTRT))
  ALLOCATE (PRINT_CONST(NCT,NWB), PRINT_HYDRO(NHY,NWB), PRINT_SEDIMENT(NWB))
  ALLOCATE (VOLUME_BALANCE(NWB),  ENERGY_BALANCE(NWB),  MASS_BALANCE(NWB))
  ALLOCATE (DETAILED_ICE(NWB),    ICE_CALC(NWB),               ALLOW_ICE(IMX))           !   ICE_IN(NBR),    RC/SW 4/28/11
  ALLOCATE (EVAPORATION(NWB),     PRECIPITATION(NWB),   RH_EVAP(NWB),         PH_CALC(NWB))
  ALLOCATE (NO_INFLOW(NWB),       NO_OUTFLOW(NWB),      NO_HEAT(NWB),         NO_WIND(NWB))
  ALLOCATE (ISO_TEMP(NWB),        VERT_TEMP(NWB),       LONG_TEMP(NWB),       VERT_PROFILE(NWB),  LONG_PROFILE(NWB))
  ALLOCATE (SNAPSHOT(NWB),        PROFILE(NWB),         VECTOR(NWB),          CONTOUR(NWB),       SPREADSHEET(NWB))
  ALLOCATE (SCREEN_OUTPUT(NWB),   FLUX(NWB))
  ALLOCATE (PRINT_DERIVED(NDC,NWB),  PRINT_EPIPHYTON(NWB,NEPT))
  ALLOCATE (SEDIMENT_CALC(NWB),      EPIPHYTON_CALC(NWB,NEPT), SEDIMENT_RESUSPENSION(NSS),BOD_CALC(NBOD),ALG_CALC(NAL))
  ALLOCATE (BOD_CALCp(NBOD), BOD_CALCn(NBOD))                                                ! cb 5/19/2011
  ALLOCATE (TDG_SPILLWAY(NWDT,NSP),  TDG_GATE(NWDT,NGT),       INTERNAL_WEIR(KMX,IMX))
  ALLOCATE (ISO_EPIPHYTON(NWB,NEPT), VERT_EPIPHYTON(NWB,NEPT), LONG_EPIPHYTON(NWB,NEPT))
  ALLOCATE (LATERAL_SPILLWAY(NSP),   LATERAL_GATE(NGT),        LATERAL_PUMP(NPU),        LATERAL_PIPE(NPI))
  ALLOCATE (INTERP_HEAD(NBR),        INTERP_WITHDRAWAL(NWD),   INTERP_EXTINCTION(NWB),   INTERP_DTRIBS(NBR))
  ALLOCATE (INTERP_OUTFLOW(NST,NBR), INTERP_INFLOW(NBR),       INTERP_METEOROLOGY(NWB),  INTERP_TRIBS(NTR))
  ALLOCATE (LNAME(NCT+NHY+NDC))
  ALLOCATE (IWR(NIW),    KTWR(NIW),   KBWR(NIW))
  ALLOCATE (JWUSP(NSP),  JWDSP(NSP),  QSP(NSP))
  ALLOCATE (KTWD(NWDT),  KBWD(NWDT),  JBWD(NWDT))
  ALLOCATE (GTA1(NGT),   GTB1(NGT),   GTA2(NGT),   GTB2(NGT))
  ALLOCATE (BGT(NGT),    IUGT(NGT),   IDGT(NGT),   EGT(NGT),    EGT2(NGT))
  ALLOCATE (QTR(NTRT),   TTR(NTRT),   KTTR(NTRT),  KBTR(NTRT))
  ALLOCATE (AGASGT(NGT), BGASGT(NGT), CGASGT(NGT), GASGTC(NGT))
  ALLOCATE (PUGTC(NGT),  ETUGT(NGT),  EBUGT(NGT),  KTUGT(NGT),  KBUGT(NGT))
  ALLOCATE (PDGTC(NGT),  ETDGT(NGT),  EBDGT(NGT),  KTDGT(NGT),  KBDGT(NGT))
  ALLOCATE (A1GT(NGT),   B1GT(NGT),   G1GT(NGT),   A2GT(NGT),   B2GT(NGT),   G2GT(NGT))
  ALLOCATE (EQGT(NGT),   JBUGT(NGT),  JBDGT(NGT),  JWUGT(NGT),  JWDGT(NGT),  QGT(NGT))
  ALLOCATE (JBUPI(NPI),  JBDPI(NPI),  JWUPI(NPI),  JWDPI(NPI),  QPI(NPI), BP(NPI))                              ! SW 5/10/10
  ALLOCATE (IUPI(NPI),   IDPI(NPI),   EUPI(NPI),   EDPI(NPI),   WPI(NPI),    DLXPI(NPI),  FPI(NPI),    FMINPI(NPI), PUPIC(NPI))
  ALLOCATE (ETUPI(NPI),  EBUPI(NPI),  KTUPI(NPI),  KBUPI(NPI),  PDPIC(NPI),  ETDPI(NPI),  EBDPI(NPI),  KTDPI(NPI),  KBDPI(NPI))
  ALLOCATE (PUSPC(NSP),  ETUSP(NSP),  EBUSP(NSP),  KTUSP(NSP),  KBUSP(NSP),  PDSPC(NSP),  ETDSP(NSP),  EBDSP(NSP))
  ALLOCATE (KTDSP(NSP),  KBDSP(NSP),  IUSP(NSP),   IDSP(NSP),   ESP(NSP),    A1SP(NSP),   B1SP(NSP),   A2SP(NSP))
  ALLOCATE (B2SP(NSP),   AGASSP(NSP), BGASSP(NSP), CGASSP(NSP), EQSP(NSP),   GASSPC(NSP), JBUSP(NSP),  JBDSP(NSP))
  ALLOCATE (IUPU(NPU),   IDPU(NPU),   EPU(NPU),    STRTPU(NPU), ENDPU(NPU),  EONPU(NPU),  EOFFPU(NPU), QPU(NPU),   PPUC(NPU))
  ALLOCATE (ETPU(NPU),   EBPU(NPU),   KTPU(NPU),   KBPU(NPU),   JWUPU(NPU),  JWDPU(NPU),  JBUPU(NPU),  JBDPU(NPU), PUMPON(NPU))
  ALLOCATE (IWD(NWDT),   KWD(NWDT),   QWD(NWDT),   EWD(NWDT),   KTW(NWDT),   KBW(NWDT))
  ALLOCATE (ITR(NTRT),   QTRFN(NTR),  TTRFN(NTR),  CTRFN(NTR),  ELTRT(NTRT), ELTRB(NTRT), TRC(NTRT),   JBTR(NTRT), QTRF(KMX,NTRT))
  ALLOCATE (TTLB(IMX),   TTRB(IMX),   CLLB(IMX),   CLRB(IMX))
  ALLOCATE (SRLB1(IMX),  SRRB1(IMX),  SRLB2(IMX),  SRRB2(IMX),  SRFJD1(IMX), SHADEI(IMX), SRFJD2(IMX))
  ALLOCATE (TOPO(IMX,IANG))                                                                                        ! SW 10/17/05
  ALLOCATE (QSW(KMX,NWDT),   HPRWBC(NHY,NWB))
  ALLOCATE (RATZ(KMX,NWB),   CURZ1(KMX,NWB),  CURZ2(KMX,NWB),   CURZ3(KMX,NWB))   ! SW 5/15/06

  ALLOCATE (zg(NZPt),zm(NZPt),zeff(NZPt),prefp(NZPt),zr(NZPt),zoomin(NZPt),zs2p(NZPt),exz(NZPt),PREFZ(NZPt,nzpt))
  ALLOCATE (zt1(NZPt),zt2(NZPt),zt3(NZPt),zt4(NZPt),zk1(NZPt),zk2(NZPt),zk3(NZPt),zk4(NZPt),o2zr(nzpt))
  ALLOCATE (zp(NZPt),zn(NZPt),zc(NZPt))
  allocate (prefa(nal,nzpt))
  allocate (po4zr(kmx,imx),nh4zr(kmx,imx))
  allocate (zmu(kmx,imx,nzp),tgraze(kmx,imx,nzp),zrt(kmx,imx,nzp),zmt(kmx,imx,nzp)) ! MLM POINTERS:,zoo(kmx,imx,NZP),zooss(kmx,imx,NZP))
  allocate (zoorm(kmx,imx,nzp),zoormr(kmx,imx,nzp),zoormf(kmx,imx,nzp))
  allocate (lpzooout(kmx,imx),lpzooin(kmx,imx),dozr(kmx,imx),ticzr(kmx,imx))
  allocate (agz(kmx,imx,nal,nzp), zgz(kmx,imx,nzp,nzp),agzt(kmx,imx,nal)) !omnivorous zooplankton
  allocate (ORGPLD(kmx,imx), ORGPRD(kmx,imx), ORGPLP(kmx,imx), ORGPRP(kmx,imx), ORGNLD(kmx,imx), ORGNRD(kmx,imx), ORGNLP(kmx,imx))
  allocate (ORGNRP(kmx,imx))
  allocate (ldompmp(kmx,imx),ldomnmp(kmx,imx),lpompmp(kmx,imx),lpomnmp(kmx,imx),rpompmp(kmx,imx),rpomnmp(kmx,imx))
  allocate (lpzooinp(kmx,imx),lpzooinn(kmx,imx),lpzoooutp(kmx,imx),lpzoooutn(kmx,imx))
  allocate (SEDVPp(KMX,NWB),SEDVPc(KMX,NWB),SEDVPn(KMX,NWB))
  allocate (sedp(kmx,imx),sedn(kmx,imx),sedc(kmx,imx))
  allocate (sdkv(kmx,imx),seddktot(kmx,imx))
  allocate (sedcip(nwb),sedcin(nwb),sedcic(nwb),sedcis(nwb))
  ALLOCATE (cbods(NBOD), cbodns(kmx,imx), sedcb(kmx,imx), sedcbp(kmx,imx), sedcbn(kmx,imx), sedcbc(kmx,imx))
  allocate  (print_macrophyte(nwb,nmct), macrophyte_calc(nwb,nmct),macwbc(nwb,nmct),conv2(kmx,kmx),mprwbc(nwb,nmct))
  allocate  (mac(kmx,imx,nmct), macrc(kmx,kmx,imx,nmct),mact(kmx,kmx,imx), macrm(kmx,kmx,imx,nmct), macss(kmx,kmx,imx,nmct))
  allocate  (mgr(kmx,kmx,imx,nmct),mmr(kmx,imx,nmct), mrr(kmx,imx,nmct))
  allocate  (smacrc(kmx,kmx,imx,nmct), smacrm(kmx,kmx,imx,nmct))
  allocate  (smact(kmx,kmx,imx), smac(kmx,imx,nmct))
  allocate  (mt1(nmct),mt2(nmct),mt3(nmct),mt4(nmct),mk1(nmct),mk2(nmct),mk3(nmct),mk4(nmct),mg(nmct),mr(nmct),mm(nmct))
  allocate  (mbmp(nmct), mmax(nmct), cddrag(nmct), dwv(nmct), dwsa(nmct), anorm(nmct))
  allocate  (mp(nmct), mn(nmct), mc(nmct),psed(nmct),nsed(nmct),mhsp(nmct),mhsn(nmct),mhsc(nmct),msat(nmct),exm(nmct))
  allocate  (O2MG(nmct), O2MR(nmct),  LRPMAC(nmct),  MPOM(nmct))
  allocate  (kticol(imx),armac(imx),macwbci(nwb,nmct))
  allocate  (macmbrs(nbr,nmct), macmbrt(nbr,nmct),ssmacmb(nbr,nmct))
  allocate  (cw(kmx,imx), bic(kmx,imx))
  allocate  (mactrmr(kmx,imx,nmct), mactrmf(kmx,imx,nmct),mactrm(kmx,imx,nmct))
  allocate  (mlfpr(kmx,kmx,imx,nmct))
  allocate  (mllim(kmx,kmx,imx,nmct), mplim(kmx,imx,nmct),mclim(kmx,imx,nmct),mnlim(kmx,imx,nmct))
  ALLOCATE  (GAMMAj(kmx,KMX,IMX))	
  allocate (por(kmx,imx),VOLKTi(imx),VOLi(Kmx,Imx),vstem(kmx,imx,nmct),vstemkt(imx,nmct),sarea(nmct))
  ALLOCATE (IWIND(NWB))
  ALLOCATE (LAYERCHANGE(NWB))
  ALLOCATE(DLTADD(NWB))
  allocate (cbodp(kmx,imx,nbod), cbodn(kmx,imx,nbod))     
  
  ! start specific allocations for control_file_read_write.f90 program
  allocate (hyd_prin(nhy),cst_icon(nct),cst_prin(nct))
  allocate  (cin_con(nct), ctr_con(nct), cdt_con(nct), cpr_con(nct))
  allocate  (gen_name(NGC),alg_name(NAL), ss_name(nss))
  allocate  (zoo_name(nzpt),epi_name(Nept),mac_name(nmct),bod_name(nbod))
  allocate  (tr_name(ntr),wb_name(nwb),br_name(nbr))
  allocate  (pipe_name(npi),spill_name(nsp), gate_name(ngt), pump_name(npu))
   ALLOCATE (BTHFNwrite(NWB), bthtype(nwb), bthline1(nwb),bthline2(nwb),vprfnwrite(nwb),vprtype(nwb),vprline1(nwb),vprline2(nwb))  
  zoo_name=blank;epi_name=blank; mac_name=blank; bod_name=blank  
! end specific allocations for control_file_read_write.f90 program


! Zero variables

  ITR  = 0;   JBTR = 0;   KTTR = 0;   KBTR = 0;   QTR  = 0.0; TTR  = 0.0; QTRF = 0.0; SNPD  = 0.0; TSRD  = 0.0
  PRFD = 0.0; SPRD = 0.0; CPLD = 0.0; VPLD = 0.0; SCRD = 0.0; FLXD = 0.0; WDOD = 0.0; RSOD = 0.0; ELTRB = 0.0; ELTRT = 0.0

! Input file unit numbers

  NUNIT = 40
  DO JW=1,NWB
    BTH(JW) = NUNIT
    VPR(JW) = NUNIT+1
    LPR(JW) = NUNIT+2
    NUNIT   = NUNIT+3
  END DO
  GRF = NUNIT; NUNIT = NUNIT+1

! Time control cards

  READ (CON,'(//8X,2F8.0,I8)')         TMSTRT,   TMEND,    YEAR
  READ (CON,'(//8X,I8,F8.0,a8)')         NDLT,     DLTMIN, DLTINTER
  READ (CON,'(//(:8X,9F8.0))')        (DLTD(J),            J =1,NDLT)
  READ (CON,'(//(:8X,9F8.0))')        (DLTMAX(J),          J =1,NDLT)
  READ (CON,'(//(:8X,9F8.0))')        (DLTF(J),            J =1,NDLT)
  READ (CON,'(//(8X,2A8))')           (VISC(JW), CELC(JW), JW=1,NWB)

! Grid definition cards

  READ (CON,'(//(8X,7I8,F8.0,F8.0))') (US(JB),  DS(JB),     UHS(JB),   DHS(JB), UQB(JB), DQB(JB),  NL(JB), SLOPE(JB),SLOPEC(JB), JB=1,NBR)
  READ (CON,'(//(8X,3F8.0,3I8))')     (LAT(JW), LONGIT(JW), ELBOT(JW), BS(JW),  BE(JW),  JBDN(JW),                    JW=1,NWB)

! Initial condition cards

  READ (CON,'(//(8X,2F8.0,2A8))')     (T2I(JW),    ICETHI(JW),  WTYPEC(JW),  GRIDC(JW),                               JW=1,NWB)
  READ (CON,'(//(8X,6A8))')           (VBC(JW),    EBC(JW),     MBC(JW),     PQC(JW),   EVC(JW),   PRC(JW),           JW=1,NWB)
  READ (CON,'(//(8X,4A8))')           (WINDC(JW),  QINC(JW),    QOUTC(JW),   HEATC(JW),                               JW=1,NWB)
  READ (CON,'(//(8X,3A8))')           (QINIC(JB),  DTRIC(JB),   HDIC(JB),                                             JB=1,NBR)
  READ (CON,'(//(8X,5A8,4F8.0))')     (SLHTC(JW),  SROC(JW),    RHEVC(JW),   METIC(JW), FETCHC(JW), AFW(JW),                       &
                                       BFW(JW),    CFW(JW),     WINDH(JW),                                            JW=1,NWB)
  READ (CON,'(//(8X,2A8,6F8.0))')     (ICEC(JW),   SLICEC(JW),  ALBEDO(JW),  HWI(JW),   BETAI(JW),  GAMMAI(JW),                    &
                                       ICEMIN(JW), ICET2(JW),                                                         JW=1,NWB)
  READ (CON,'(//(8X,A8,F8.0))')       (SLTRC(JW),  THETA(JW),                                                         JW=1,NWB)
  READ (CON,'(//(8X,6F8.0,A8,F8.0))')      (AX(JW),     DXI(JW),     CBHE(JW),    TSED(JW),  FI(JW),     TSEDF(JW),                     &
                                       FRICC(JW), Z0(JW),                                                                    JW=1,NWB)
  READ (CON,'(//(8X,2A8,F8.0,I8,F8.0,F8.0,F8.0,F8.0,A8))')     (AZC(JW),    AZSLC(JW),   AZMAX(JW),   TKEBC(JW),EROUGH(JW),       &
                                       ARODI(JW),STRICK(JW),TKELATPRDCONST(JW),IMPTKE(JW),JW=1,NWB)          !,PHISET(JW

  DO JW=1,NWB
  IF(Z0(JW) <= 0.0)Z0(JW)=0.001      ! SW 11/28/07
   DO JB=BS(JW),BE(JW)
    DO I=US(JB),DS(JB)
    E(I)=EROUGH(JW)
    ENDDO
   ENDDO
  ENDDO

! Inflow-outflow cards
  STRIC='     OFF'
  KTSWT=0
  KBSWT=0
  SINKCT='        '
  ESTRT=0.0
  WSTRT=0.0

  READ (CON,'(//(8X,I8,A8))')            (NSTR(JB), DYNSTRUC(JB),     JB=1,NBR)
  READ (CON,'(/)')
  DO JB=1,NBR
    READ (CON,'(:8X,9A8)')            (STRIC(JS,JB),  JS=1,NSTR(JB))
  END DO
  READ (CON,'(/)')
  DO JB=1,NBR
    READ (CON,'(:8X,9I8)')            (KTSWT(JS,JB), JS=1,NSTR(JB))
  END DO
  READ (CON,'(/)')
  DO JB=1,NBR
    READ (CON,'(:8X,9I8)')            (KBSWT(JS,JB), JS=1,NSTR(JB))
  END DO
  READ (CON,'(/)')
  DO JB=1,NBR
    READ (CON,'(:8X,9A8)')            (SINKCT(JS,JB),JS=1,NSTR(JB))
  END DO
  READ (CON,'(/)')
  DO JB=1,NBR
    READ (CON,'(:8X,9F8.0)')          (ESTRT(JS,JB), JS=1,NSTR(JB))
  END DO
  READ (CON,'(/)')
  DO JB=1,NBR
    READ (CON,'(:8X,9F8.0)')          (WSTRT(JS,JB), JS=1,NSTR(JB))
  END DO
  READ (CON,'(//(:a8,2I8,6F8.0,A8,A8))') (pipe_name(jp),IUPI(JP),   IDPI(JP),   EUPI(JP),   EDPI(JP),    WPI(JP),                                   &
                                       DLXPI(JP),  FPI(JP),    FMINPI(JP), LATPIC(JP),  DYNPIPE(JP),JP=1,NPI)
  READ (CON,'(//(:8X,A8,2F8.0,2I8))') (PUPIC(JP),  ETUPI(JP),  EBUPI(JP),  KTUPI(JP),   KBUPI(JP),  JP=1,NPI)
  READ (CON,'(//(:8X,A8,2F8.0,2I8))') (PDPIC(JP),  ETDPI(JP),  EBDPI(JP),  KTDPI(JP),   KBDPI(JP),  JP=1,NPI)
  READ (CON,'(//(:a8,2I8,5F8.0,A8))') (spill_name(js),IUSP(JS),   IDSP(JS),   ESP(JS),    A1SP(JS),    B1SP(JS),                                  &
                                       A2SP(JS),   B2SP(JS),   LATSPC(JS),                          JS=1,NSP)
  READ (CON,'(//(:8X,A8,2F8.0,2I8))') (PUSPC(JS),  ETUSP(JS),  EBUSP(JS),  KTUSP(JS),   KBUSP(JS),  JS=1,NSP)
  READ (CON,'(//(:8X,A8,2F8.0,2I8))') (PDSPC(JS),  ETDSP(JS),  EBDSP(JS),  KTDSP(JS),   KBDSP(JS),  JS=1,NSP)
  READ (CON,'(//(:8X,A8,I8,3F8.0))')  (GASSPC(JS), EQSP(JS),   AGASSP(JS), BGASSP(JS),  CGASSP(JS), JS=1,NSP)
  READ (CON,'(//(:a8,2I8,7F8.0,A8))') (gate_name(jg),IUGT(JG),   IDGT(JG),   EGT(JG),    A1GT(JG),    B1GT(JG),                                  &
                                       G1GT(JG),   A2GT(JG),   B2GT(JG),   G2GT(JG),    LATGTC(JG), JG=1,NGT)
  READ (CON,'(//(:8X,4F8.0,2A8))')     (GTA1(JG),   GTB1(JG),   GTA2(JG),   GTB2(JG),    DYNGTC(JG),gtic(jg), JG=1,NGT)  ! cb 8/13/2010
  READ (CON,'(//(:8X,A8,2F8.0,2I8))') (PUGTC(JG),  ETUGT(JG),  EBUGT(JG),  KTUGT(JG),   KBUGT(JG),  JG=1,NGT)
  READ (CON,'(//(:8X,A8,2F8.0,2I8))') (PDGTC(JG),  ETDGT(JG),  EBDGT(JG),  KTDGT(JG),   KBDGT(JG),  JG=1,NGT)
  READ (CON,'(//(:8X,A8,I8,3F8.0))')  (GASGTC(JG), EQGT(JG),   AGASGT(JG), BGASGT(JG),  CGASGT(JG), JG=1,NGT)
  READ (CON,'(//(:a8,2I8,6F8.0,2A8))') (pump_name(jp),IUPU(JP),   IDPU(JP),   EPU(JP),    STRTPU(JP),  ENDPU(JP),                                 &
                                       EONPU(JP),  EOFFPU(JP), QPU(JP),    LATPUC(JP),  DYNPUMP(JP),      JP=1,NPU)
  READ (CON,'(//(:8X,A8,2F8.0,2I8))') (PPUC(JP),   ETPU(JP),   EBPU(JP),   KTPU(JP),    KBPU(JP),   JP=1,NPU)
  READ (CON,'(//(:8X,9I8))')          (IWR(JW),    JW=1,NIW)
  READ (CON,'(//(:8X,9I8))')          (KTWR(JW),   JW=1,NIW)
  READ (CON,'(//(:8X,9I8))')          (KBWR(JW),   JW=1,NIW)
  READ (CON,'(//(:8X,9A8))')          (WDIC(JW),   JW=1,NWD)
  READ (CON,'(//(:8X,9I8))')          (IWD(JW),    JW=1,NWD)
  READ (CON,'(//(:8X,9F8.0))')        (EWD(JW),    JW=1,NWD)
  READ (CON,'(//(:8X,9I8))')          (KTWD(JW),   JW=1,NWD)
  READ (CON,'(//(:8X,9I8))')          (KBWD(JW),   JW=1,NWD)
  READ (CON,'(//(:8X,9A8))')          (TRC(JT),    JT=1,NTR)
  READ (CON,'(//(:8X,9A8))')          (TRIC(JT),   JT=1,NTR)
  READ (CON,'(//(:8X,9I8))')          (ITR(JT),    JT=1,NTR)
  READ (CON,'(//(:8X,9F8.0))')        (ELTRT(JT),  JT=1,NTR)
  READ (CON,'(//(:8X,9F8.0))')        (ELTRB(JT),  JT=1,NTR)
  READ (CON,'(//(8X,A8))')            (DTRC(JB),   JB=1,NBR)

! Output control cards (excluding constituents)

  READ (CON,'(/)')
  DO JH=1,NHY
    READ (CON,'(10a8:/(:8x,9a8))')            hyd_prin(jh),(HPRWBC(JH,JW),JW=1,NWB)
  END DO
  READ (CON,'(//(8X,A8,2I8))')        (SNPC(JW), NSNP(JW), NISNP(JW), JW=1,NWB)
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (SNPD(J,JW),J=1,NSNP(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (SNPF(J,JW),J=1,NSNP(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9I8)')            (ISNP(I,JW),I=1,NISNP(JW))
  END DO
  READ (CON,'(//(8X,A8,I8))')         (SCRC(JW), NSCR(JW), JW=1,NWB)
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (SCRD(J,JW),J=1,NSCR(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (SCRF(J,JW),J=1,NSCR(JW))
  END DO
  READ (CON,'(//(8X,A8,2I8))')        (PRFC(JW), NPRF(JW), NIPRF(JW), JW=1,NWB)
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (PRFD(J,JW),J=1,NPRF(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (PRFF(J,JW),J=1,NPRF(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9I8)')            (IPRF(J,JW),J=1,NIPRF(JW))
  END DO
  READ (CON,'(//(8X,A8,2I8))')        (SPRC(JW), NSPR(JW), NISPR(JW), JW=1,NWB)
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (SPRD(J,JW),J=1,NSPR(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (SPRF(J,JW),J=1,NSPR(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9I8)')            (ISPR(J,JW), J=1,NISPR(JW))
  END DO
  READ (CON,'(//(8X,A8,I8))')         (VPLC(JW),  NVPL(JW),  JW=1,NWB)
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (VPLD(J,JW), J=1,NVPL(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (VPLF(J,JW), J=1,NVPL(JW))
  END DO
  READ (CON,'(//(8X,A8,I8,A8))')      (CPLC(JW),   NCPL(JW), TECPLOT(JW),JW=1,NWB)
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (CPLD(J,JW), J=1,NCPL(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (CPLF(J,JW), J=1,NCPL(JW))
  END DO
  READ (CON,'(//(8X,A8,I8))')         (FLXC(JW),   NFLX(JW), JW=1,NWB)
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (FLXD(J,JW), J=1,NFLX(JW))
  END DO
  READ (CON,'(/)')
  DO JW=1,NWB
    READ (CON,'(:8X,9F8.0)')          (FLXF(J,JW), J=1,NFLX(JW))
  END DO
  READ (CON,'(//8X,A8,2I8)')           TSRC,    NTSR,    NIKTSR; ALLOCATE (ITSR(MAX(1,NIKTSR)), ETSR(MAX(1,NIKTSR)))
  READ (CON,'(//(:8X,9F8.0))')        (TSRD(J), J=1,NTSR)
  READ (CON,'(//(:8X,9F8.0))')        (TSRF(J), J=1,NTSR)
  READ (CON,'(//(:8X,9I8))')          (ITSR(J), J=1,NIKTSR)
  READ (CON,'(//(:8X,9F8.0))')        (ETSR(J), J=1,NIKTSR)
  READ (CON,'(//8X,A8,2I8)')           WDOC,    NWDO,    NIWDO;  ALLOCATE (IWDO(MAX(1,NIWDO)))
  READ (CON,'(//(:8X,9F8.0))')        (WDOD(J), J=1,NWDO)
  READ (CON,'(//(:8X,9F8.0))')        (WDOF(J), J=1,NWDO)
  READ (CON,'(//(8X,9I8))')           (IWDO(J), J=1,NIWDO)
  READ (CON,'(//8X,A8,I8,A8)')         RSOC,    NRSO,    RSIC
  READ (CON,'(//(:8X,9F8.0))')        (RSOD(J), J=1,NRSO)
  READ (CON,'(//(:8X,9F8.0))')        (RSOF(J), J=1,NRSO)

! Constituent control cards

  READ (CON,'(//8X,2A8,I8,F8.0,A8)')           CCC, LIMC, CUF,CO2PPM,CO2YRLY
  IF(CO2YRLY /='      ON' .AND.  CO2YRLY /='     OFF')CO2YRLY ='      ON'
  READ (CON,'(//(2A8))')              (CNAME2(JC),  CAC(JC),      JC=1,NCT)
  READ (CON,'(/)')
  
  DO JD=1,NDC2
    READ (CON,'(A8,9A8:/(:8X,9A8))')           CDNAME2(JD),(CDWBC(JD,JW), JW=1,NWB)
  END DO
  READ (CON,'(/)')
  KFNAME2 ='     '   ! SW 9/27/13 INITIALIZE ENTIRE ARRAY
  KFWBC   ='     '   ! SW 9/27/13 INITIALIZE ENTIRE ARRAY
!  DO JF=1,NFL
  do jf=1,73   ! Fix this later
    if(nwb < 10)READ (CON,'(A8,(:9A8))')         KFNAME2(JF),(KFWBC(JF,JW),  JW=1,NWB)
    if(nwb >= 10)READ (CON,'(A8,9A8:/(:8X,9A8))')         KFNAME2(JF),(KFWBC(JF,JW),  JW=1,NWB)          !cb 9/13/12  sw2/18/13  Foramt 6/16/13 8/13/13
  END DO
  !kfname2(121) = 'CO2GASX'

  READ (CON,'(/)')
  DO JC=1,NCT
    READ (CON,'(a8,9F8.0:/(:8X,9F8.0))')          cst_icon(jc),(C2I(JC,JW),    JW=1,NWB)
  END DO
  READ (CON,'(/)')
  DO JC=1,NCT
    READ (CON,'(a8,9A8:/(:8X,9A8))')            cst_prin(jc),(CPRWBC(JC,JW), JW=1,NWB)
  END DO
  READ (CON,'(/)')
  DO JC=1,NCT
    READ (CON,'(a8,9A8:/(:8X,9A8))')            cin_con(jc),(CINBRC(JC,JB), JB=1,NBR)
  END DO
  READ (CON,'(/)')
  DO JC=1,NCT
    READ (CON,'(a8,9A8:/(:8X,9A8))')            ctr_con(jc),(CTRTRC(JC,JT), JT=1,NTR)
  END DO
  IF(NTR==0)THEN
      CTRTRC(:,1)='     OFF'
  ENDIF
  
  READ (CON,'(/)')
  DO JC=1,NCT
    READ (CON,'(a8,9A8:/(:8X,9A8))')            cdt_con(jc),(CDTBRC(JC,JB), JB=1,NBR)
  END DO
  READ (CON,'(/)')
  DO JC=1,NCT
    READ (CON,'(a8,9A8:/(:8X,9A8))')            cpr_con(jc),(CPRBRC(JC,JB), JB=1,NBR)
  END DO

! Kinetics coefficients

  READ (CON,'(//(8X,4F8.0,2A8))')     (EXH2O(JW),  EXSS(JW),   EXOM(JW),   BETA(JW),   EXC(JW),   EXIC(JW),    JW=1,NWB)
  READ (CON,'(//(8X,9F8.0))')         (EXA(JA),                                                                JA=1,NAL)
  READ (CON,'(//(8X,9F8.0))')         (EXZ(Jz),                                                                Jz=1,Nzpt)  !v3.5
  READ (CON,'(//(8X,9F8.0))')         (EXM(Jm),                                                                Jm=1,nmct)  !v3.5
  !READ (CON,'(//(a8,4F8.0))')         (gen_name(Jg),CGQ10(JG),  CG0DK(JG),  CG1DK(JG),  CGS(JG),                            JG=1,NGC)
  READ (CON,'(//(a8,7F8.0))')         (gen_name(Jg),CGQ10(JG),  CG0DK(JG),  CG1DK(JG),  CGS(JG), CGLDK(JG),CGKLF(JG),CGCS(JG),           JG=1,NGC) !LCJ 2/26/15
  READ (CON,'(//(a8,F8.0,A,F8.0))')   (ss_name(js),SSS(JS),    SEDRC(JS),  TAUCR(JS),                                      JS=1,NSS)
  READ (CON,'(//(a8,9F8.0))')         (alg_name(ja),AG(JA),     AR(JA),     AE(JA),     AM(JA),     AS(JA),                                     &
                                       AHSP(JA),   AHSN(JA),   AHSSI(JA),  ASAT(JA),                           JA=1,NAL)
  READ (CON,'(//(8X,8F8.0))')         (AT1(JA),    AT2(JA),    AT3(JA),    AT4(JA),    AK1(JA),   AK2(JA),                         &
                                       AK3(JA),    AK4(JA),                                                    JA=1,NAL)
  READ (CON,'(//(8X,6F8.0,I8,F8.0))') (AP(JA),     AN(JA),     AC(JA),     ASI(JA),    ACHLA(JA), APOM(JA),                        &
                                       ANEQN(JA),  ANPR(JA),   JA=1,NAL)
  READ (CON,'(//(8X,9A8))')           (EPIC(JW,1),                                                             JW=1,NWB)
  DO JE=2,NEPT
    READ (CON,'(8X,9A8)')             (EPIC(JW,JE),                                                            JW=1,NWB)
  END DO
  READ (CON,'(//(8X,9A8))')           (EPIPRC(JW,1),                                                           JW=1,NWB)
  DO JE=2,NEPT
    READ (CON,'(8X,9A8)')             (EPIPRC(JW,JE),                                                          JW=1,NWB)
  END DO
  READ (CON,'(//(8X,9F8.0))')         (EPICI(JW,1),                                                            JW=1,NWB)
  DO JE=2,NEPT
    READ (CON,'(8X,9F8.0)')           (EPICI(JW,JE),                                                           JW=1,NWB)
  END DO
  READ (CON,'(//(a8,8F8.0))')         (epi_name(je),EG(JE),     ER(JE),     EE(JE),     EM(JE),     EB(JE),    EHSP(JE),                        &
                                       EHSN(JE),   EHSSI(JE),                                                  JE=1,NEP)
  READ (CON,'(//(8X,2F8.0,I8,F8.0))') (ESAT(JE),   EHS(JE),    ENEQN(JE),  ENPR(JE),                           JE=1,NEP)
  READ (CON,'(//(8X,8F8.0))')         (ET1(JE),    ET2(JE),    ET3(JE),    ET4(JE),    EK1(JE),   EK2(JE),                         &
                                       EK3(JE),    EK4(JE),                                                    JE=1,NEP)
  READ (CON,'(//(8X,6F8.0))')         (EP(JE),     EN(JE),     EC(JE),     ESI(JE),    ECHLA(JE), EPOM(JE),    JE=1,NEP)
! v3.5 start
  READ (CON,'(//(a8,7F8.0))')         (zoo_name(jz),zg(jz),zr(jz),zm(jz),zeff(jz),PREFP(jz),ZOOMIN(jz),ZS2P(jz),            Jz=1,Nzpt)

  READ (CON,'(//(8X,8F8.0))')         (PREFA(ja,1),                                                            Ja=1,nal)          ! MM 7/13/06
  do jz=2,nzpt
    READ (CON,'((8X,8F8.0))')       (PREFA(ja,jz),                                                           Ja=1,nal)
  end do
  READ (CON,'(//(8X,8F8.0))')       (PREFz(jjz,1),                                                          Jjz=1,nzpt)
  do jz=2,nzpt
    READ (CON,'((8X,8F8.0))')       (PREFz(jjz,jz),                                                          Jjz=1,nzpt)           ! MM 7/13/06
  end do
  READ (CON,'(//(8X,8F8.0))')         (zT1(Jz),    zT2(Jz),    zT3(Jz),    zT4(Jz),    zK1(Jz),   zK2(Jz),                         &
                                       zK3(Jz),    zK4(Jz),                                                    Jz=1,Nzpt)
  READ (CON,'(//(8X,3F8.0))')         (zP(Jz),     zN(Jz),     zC(Jz),                                         Jz=1,Nzpt)
  READ (CON,'(//(8X,9A8))')           (macwbC(JW,1),                                                           JW=1,NWB)
  DO Jm=2,nmct
    READ (CON,'(8X,9A8)')             (macwbC(JW,Jm),                                                          JW=1,NWB)
  END DO
  READ (CON,'(//(8X,9A8))')           (mprwbC(JW,1),                                                           JW=1,NWB)
  DO Jm=2,nmct
    READ (CON,'(8X,9A8)')             (mprwbC(JW,Jm),                                                          JW=1,NWB)
  END DO
  READ (CON,'(//(8X,9F8.0))')         (macwbCI(JW,1),                                                          JW=1,NWB)
  DO Jm=2,nmct
    READ (CON,'(8X,9F8.0)')           (macwbcI(JW,Jm),                                                         JW=1,NWB)
  END DO
  READ (CON,'(//(a8,9F8.0))')         (mac_name(jm),mG(jm), mR(jm), mM(jm), msat(jm),mhsp(jm),mhsn(jm),mhsc(jm),                           &
                                          mpom(jm),lrpmac(jm),     jm=1,nmct)
  READ (CON,'(//(8X,2F8.0))')         (psed(jm), nsed(jm),                                                     jm=1,nmct)
  READ (CON,'(//(8X,2F8.0))')         (mbmp(jm), mmax(jm),                                                     jm=1,nmct)
  READ (CON,'(//(8X,4F8.0))')         (cddrag(jm),dwv(jm),dwsa(jm),anorm(jm),                                 jm=1,nmct)  !cb 6/29/06
  READ (CON,'(//(8X,8F8.0))')         (mT1(Jm),    mT2(Jm),    mT3(Jm),    mT4(Jm),    mK1(Jm),   mK2(Jm),                         &
                                       mK3(Jm),    mK4(Jm),                                                    Jm=1,nmct)
  READ (CON,'(//(8X,3F8.0))')         (mP(Jm),     mN(Jm),     mC(Jm),                                         Jm=1,nmct)
  READ (CON,'(//(8X,3F8.0))')         (LDOMDK(JW), RDOMDK(JW), LRDDK(JW),                                      JW=1,NWB)
  READ (CON,'(//(8X,4F8.0))')         (LPOMDK(JW), RPOMDK(JW), LRPDK(JW),  POMS(JW),                           JW=1,NWB)
  READ (CON,'(//(8X,4F8.0))')         (ORGP(JW),   ORGN(JW),   ORGC(JW),   ORGSI(JW),                          JW=1,NWB)
  READ (CON,'(//(8X,4F8.0))')         (OMT1(JW),   OMT2(JW),   OMK1(JW),   OMK2(JW),                           JW=1,NWB)
  READ (CON,'(//(8X,4F8.0))')         (KBOD(JB),   TBOD(JB),   RBOD(JB), cbods(jb),                           JB=1,NBOD)  !v3.5
  READ (CON,'(//(8X,3F8.0))')         (BODP(JB),   BODN(JB),   BODC(JB),                                       JB=1,NBOD)
  READ (CON,'(//(8X,2F8.0))')         (PO4R(JW),   PARTP(JW),                                                  JW=1,NWB)
  READ (CON,'(//(8X,2F8.0))')         (NH4R(JW),   NH4DK(JW),                                                  JW=1,NWB)
  READ (CON,'(//(8X,4F8.0))')         (NH4T1(JW),  NH4T2(JW),  NH4K1(JW),  NH4K2(JW),                          JW=1,NWB)
  READ (CON,'(//(8X,3F8.0))')         (NO3DK(JW),  NO3S(JW),   FNO3SED(JW),                                    JW=1,NWB)
  READ (CON,'(//(8X,4F8.0))')         (NO3T1(JW),  NO3T2(JW),  NO3K1(JW),  NO3K2(JW),                          JW=1,NWB)
  READ (CON,'(//(8X,4F8.0))')         (DSIR(JW),   PSIS(JW),   PSIDK(JW),  PARTSI(JW),                         JW=1,NWB)
  READ (CON,'(//(8X,2F8.0))')         (FER(JW),    FES(JW),                                                    JW=1,NWB)
  READ (CON,'(//(8X,F8.0))')          (CO2R(JW),                                                               JW=1,NWB)
  READ (CON,'(//(8X,2F8.0))')         (O2NH4(JW),  O2OM(JW),                                                   JW=1,NWB)
  READ (CON,'(//(8X,2F8.0))')         (O2AR(JA),   O2AG(JA),                                                   JA=1,NAL)
  READ (CON,'(//(8X,2F8.0))')         (O2ER(JE),   O2EG(JE),                                                   JE=1,NEPT)
  READ (CON,'(//(8X,F8.0))')          (O2zR(Jz),                                                               Jz=1,Nzpt)
  READ (CON,'(//(8X,2F8.0))')         (O2mR(Jm),   O2mG(jm),                                                   Jm=1,nmct)
  READ (CON,'(//(8X,F8.0))')           KDO
  READ (CON,'(//(8X,2A8,6F8.0,A8))')     (SEDCc(JW),   SEDPRC(JW), SEDCI(JW),  SDK(JW), seds(jw),   FSOD(JW),   FSED(JW), sedb(jw),DYNSEDK(JW),   JW=1,NWB)  ! cb 11/28/06
  READ (CON,'(//(8X,4F8.0))')         (SODT1(JW),  SODT2(JW),  SODK1(JW),  SODK2(JW),                          JW=1,NWB)
  READ (CON,'(//(8X,9F8.0))')         (SOD(I),                                                                  I=1,IMX)
  READ (CON,'(//(8X,A8,I8,4F8.2))')   (REAERC(JW), NEQN(JW),   RCOEF1(JW), RCOEF2(JW), RCOEF3(JW), RCOEF4(JW), JW=1,NWB)

! Input filenames

  READ (CON,'(//(8X,A72))')  RSIFN
  READ (CON,'(//(8X,A72))')  QWDFN
  READ (CON,'(//(8X,A72))')  QGTFN
  READ (CON,'(//(8X,A72))')  WSCFN
  READ (CON,'(//(8X,A72))')  SHDFN
  READ (CON,'(//(8X,A72))') (BTHFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (METFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (EXTFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (VPRFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (LPRFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (QINFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (TINFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (CINFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (QOTFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (QTRFN(JT), JT=1,NTR)
  READ (CON,'(//(8X,A72))') (TTRFN(JT), JT=1,NTR)
  READ (CON,'(//(8X,A72))') (CTRFN(JT), JT=1,NTR)
  READ (CON,'(//(8X,A72))') (QDTFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (TDTFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (CDTFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (PREFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (TPRFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (CPRFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (EUHFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (TUHFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (CUHFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (EDHFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (TDHFN(JB), JB=1,NBR)
  READ (CON,'(//(8X,A72))') (CDHFN(JB), JB=1,NBR)

! Output filenames

  READ (CON,'(//(8X,A72))') (SNPFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (PRFFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (VPLFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (CPLFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (SPRFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))') (FLXFN(JW), JW=1,NWB)
  READ (CON,'(//(8X,A72))')  TSRFN
  READ (CON,'(//(8X,A72))')  WDOFN
  CLOSE (CON)

! Bathymetry file

  DO JW=1,NWB
    OPEN (BTH(JW),FILE=BTHFN(JW),STATUS='OLD')
	READ  (BTH(JW),'(a1)')char1                                 ! New Bathymetry format option SW 6/22/09
      IF(CHAR1=='$')THEN
      READ  (BTH(JW),*)
      READ  (BTH(JW),*) AID,(DLX(I),  I=US(BS(JW))-1,DS(BE(JW))+1)
      READ  (BTH(JW),*) AID,(ELWS(I), I=US(BS(JW))-1,DS(BE(JW))+1)
      READ  (BTH(JW),*) AID,(PHI0(I), I=US(BS(JW))-1,DS(BE(JW))+1)
      READ  (BTH(JW),*) AID,(FRIC(I), I=US(BS(JW))-1,DS(BE(JW))+1)
      READ  (BTH(JW),*)
      DO K=1,KMX
      READ  (BTH(JW),*) H(K,JW),(B(K,I),I=US(BS(JW))-1,DS(BE(JW))+1)
      END DO
      DO I=US(BS(JW))-1,DS(BE(JW))+1
      H2(:,I) = H(:,JW)
      END DO
      ELSE	
    rewind(bth(jw))
    read(bth(jw),'(a)')bthline1
    read(bth(jw),'(a)')bthline2
    read(bth(jw),*)
    READ (BTH(JW),'((10F8.0))') (DLX(I),  I=US(BS(JW))-1,DS(BE(JW))+1)
    READ (BTH(JW),'(//(10F8.0))') (ELWS(I), I=US(BS(JW))-1,DS(BE(JW))+1)
    READ (BTH(JW),'(//(10F8.0))') (PHI0(I), I=US(BS(JW))-1,DS(BE(JW))+1)
    READ (BTH(JW),'(//(10F8.0))') (FRIC(I), I=US(BS(JW))-1,DS(BE(JW))+1)
    READ (BTH(JW),'(//(10F8.0))') (H(K,JW), K=1,KMX)
    DO I=US(BS(JW))-1,DS(BE(JW))+1
      READ (BTH(JW),'(//(10F8.0))') (B(K,I), K=1,KMX)
      H2(:,I) = H(:,JW)
    END DO
       if(jw<10)then
          write(ichar1,'(i1)')jw
          bthfn(jw)='bth'//ichar1//'.csv'
       else
          write(ichar2,'(i2)')jw
          bthfn(jw)='bth'//ichar2//'.csv'
       endif

    OPEN (11,FILE=BTHFN(jw),STATUS='UNKNOWN')
    
	WRITE  (11,'(a1)')'$'                                ! New Bathymetry format option SW 6/22/09
    
      Write(11,'(",",1000(i3,","))')(I,I=us(bs(jw))-1,ds(be(jw))+1 )
      Write(11,'("DLX,",1000(f10.3,","))') (DLX(I), I=us(bs(jw))-1,ds(be(jw))+1)
      Write(11,'("ELWS,",1000(f10.3,","))') (ELWS(I), I=us(bs(jw))-1,ds(be(jw))+1)
      Write(11,'("PHI0,",1000(f10.3,","))') (PHI0(I), I=us(bs(jw))-1,ds(be(jw))+1)
      Write(11,'("FRIC,",1000(f10.3,","))') (FRIC(I), I=us(bs(jw))-1,ds(be(jw))+1)
      Write(11,'("LAYERH,",<IMX>(","),"K,ELEV")')
      DO K=1,KMX
      WRITE (11,'(f8.3,",",<imx>(f10.3,","),i3)') H(K,jw),(B(K,I),I=us(bs(jw))-1,ds(be(jw))+1),K
      END DO
    close(11)
    
      endif
    CLOSE (BTH(JW))
  END DO
  H1 = H2
  BI = B

  ALLOCATE(BSAVE(KMX,IMX))
  BSAVE=0.0
  BSAVE = B
  
  npt=11
  OPEN (NPT,FILE='graph.npt',STATUS='OLD')                                           ! SW 1/16/04 replaced section
  READ  (NPT,'(///(A43,1X,A9,3F8.0,A8))') (HNAME(J),  FMTH(J),  HMULT(J),  HYMIN(J), HYMAX(J), HPLTC(J), J=1,NHY)
  READ  (NPT,'(// (A43,1X,A9,3F8.0,A8))') (CNAME(J),  FMTC(J),  CMULT(J),  CMIN(J),  CMAX(J),  CPLTC(J), J=1,NCT)
  READ  (NPT,'(// (A43,1X,A9,3F8.0,A8))') (CDNAME(J), FMTCD(J), CDMULT(J), CDMIN(J), CDMAX(J), CDPLTC(J),J=1,NDC2)
  CLOSE(NPT)

  
  DO JC=1,NCT
    L1         = SCAN (CNAME(JC),',')+2
    L2         = SCAN (CNAME(JC)(L1:43),'  ')+L1
    CUNIT(JC)  = CNAME(JC)(L1:L2)
    CNAME1(JC) = CNAME(JC)(1:L1-3)
    CUNIT1(JC) = CUNIT(JC)(1:1)
    CUNIT2(JC) = CUNIT(JC)
    IF (CUNIT(JC)(1:2) == 'mg') THEN
      CUNIT1(JC) = 'g'
      CUNIT2(JC) = 'gm^3'
    END IF
    IF (CUNIT(JC)(1:2) /= 'g' .AND. CUNIT(JC)(1:2) /= 'mg') CUNIT1(JC) = '  '
  END DO
  
  
  ! Initialize logical control variables
  CONSTITUENTS =  CCC  == '      ON'
  VERT_PROFILE = .FALSE.  
  DO JW=1,NWB    
    VERT_TEMP(JW)        = T2I(JW)     == -1    
    IF(CONSTITUENTS)THEN                     ! CB 12/04/08    
    VERT_SEDIMENT(JW)    = SEDCI(JW)   == -1.0 .AND. SEDCC(JW)   == '      ON'    
    VERT_EPIPHYTON(JW,:) = EPICI(JW,:) == -1.0 .AND. EPIC(JW,:) == '      ON'
    DO JC=1,NCT      
      VERT_CONC(JC,JW) = C2I(JC,JW) == -1.0 .AND. CAC(JC) == '      ON'      
      IF (VERT_CONC(JC,JW)) VERT_PROFILE(JW) = .TRUE.      
    END DO
    IF (VERT_SEDIMENT(JW))         VERT_PROFILE(JW) = .TRUE.    
    IF (ANY(VERT_EPIPHYTON(JW,:))) VERT_PROFILE(JW) = .TRUE.    
    END IF                          ! cb 12/04/08
    IF (VERT_TEMP(JW))             VERT_PROFILE(JW) = .TRUE.    
  END DO
  PH_CALC               = CONSTITUENTS .AND. CDWBC(20,:) == '      ON'
  IZMIN  = 0;   KTWB   = 2;   KMIN   = 1;   IMIN   = 1; KBMAX  = 0
  
  ! Convert slope to angle alpha in radians

  ALPHA = ATAN(SLOPE)
  SINA  = SIN(ALPHA)
  SINAC = SIN(ATAN(SLOPEC))
  COSA  = COS(ALPHA)
    
  TIME_SERIES           = TSRC        == '      ON'
  if(WDOC== '      ON' .or. WDOC   == '     ONH' .or. WDOC        == '     ONS')DOWNSTREAM_OUTFLOW = .TRUE.  ! cb 4/11/18
  SEDIMENT_CALC         = CONSTITUENTS .AND. SEDCc        == '      ON'
  ICE_CALC              = ICEC        == '      ON'
  ICE_COMPUTATION       = ANY(ICE_CALC)
  DERIVED_CALC          = CONSTITUENTS .AND. ANY(CDWBC   == '      ON')
  VECTOR             = VPLC   == '      ON'
  
  NAC    = 0;   NTAC   = 0;   NACD   = 0;   NACIN  = 0;   NACTR  = 0;   NACDT  = 0;   NACPR  = 0; NAF    = 0
  IF (CONSTITUENTS) THEN
    DO JC=1,NCT
      IF (CAC(JC) == '      ON') THEN
        NAC     = NAC+1
        CN(NAC) = JC
      END IF
      DO JB=1,NBR
        IF (CINBRC(JC,JB) == '      ON') THEN
          NACIN(JB)       = NACIN(JB)+1
          INCN(NACIN(JB),JB) = JC
        END IF
        IF (CDTBRC(JC,JB) == '      ON') THEN
          NACDT(JB)       = NACDT(JB)+1
          DTCN(NACDT(JB),JB) = JC
        END IF
        IF (CPRBRC(JC,JB) == '      ON') THEN
          NACPR(JB)       = NACPR(JB)+1
          PRCN(NACPR(JB),JB) = JC
        END IF
      END DO
      DO JT=1,NTR
        IF (CTRTRC(JC,JT) == '      ON') THEN
          NACTR(JT)          = NACTR(JT)+1
          TRCN(NACTR(JT),JT) = JC
        END IF
      END DO
    END DO
    DO JW=1,NWB
      DO JD=1,NDC
        IF (CDWBC(JD,JW) == '      ON') THEN
          NACD(JW)         = NACD(JW)+1
          CDN(NACD(JW),JW) = JD
        END IF
      END DO
      DO JF=1,NFL
        IF (KFWBC(JF,JW) == '      ON') THEN
          NAF(JW)          = NAF(JW)+1
          KFCN(NAF(JW),JW) = JF
        elseif(ph_calc(jw).and. jf==121)then
          NAF(JW)          = NAF(JW)+1
          KFCN(NAF(JW),JW) = JF       
        END IF
      END DO
    END DO
  END IF
  
  return
END SUBROUTINE


!**************************************************************
 SUBROUTINE CONTROL_FILE_CSV1
  USE MAIN
  USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  use macrophytec; use porosityc; use zooplanktonc;   use initialvelocity  ;  use extra, ONLY:CONFNwrite
  IMPLICIT NONE
  
  INTEGER :: N    !,NNDC


  open(CON,file=CONFNwrite,status='unknown')

! Title and array dimensions
  
  write (CON,'(A,",")')'CE-QUAL-W2 Version,,4.5'
  write (CON,'(A,",")')'Control File version,,4.5,w2_con45.csv'
  write (CON,'(A,",")')'Title comments: next 10 lines'
  DO J=1,10
  write (CON,'(A,A,A)')'"',TRIM(TITLE(J)),'",'
  ENDDO
  write (CON,*)
  write (CON,'(A,",")')'NWB, NBR, IMX, KMX, NPROC, CLOSEC'
  write (CON,'(5(I8,","),A,",")') NWB, NBR, IMX, KMX, NPROC, ADJUSTL(CLOSEC)  
  write (CON,*)
  write (CON,'(A,",")')'NTR, NST, NIW, NWD, NGT, NSP, NPI, NPU'
  write (CON,'(8(I4,","),",")') NTR, NST, NIW, NWD, NGT, NSP, NPI, NPU
  write (CON,*)
  write (CON,'(A,",")')'NGC, NSS, NAL, NEP, NBOD, NMC, NZP'
  write (CON,'(7(I4,","),",")') NGC, NSS, NAL, NEP, NBOD, NMC, NZP  
  write (CON,*)
  write (CON,'(A,",")')'NDAY,SELECTC,HABTATC,ENVIRPC,AERATEC,INITUWL'
  write (CON,'(I3,",",5(A,","))') NOD,ADJUSTL(SELECTC),ADJUSTL(HABTATC),ADJUSTL(ENVIRPC),ADJUSTL(AERATEC),ADJUSTL(INITUWL)           !'(I0,5(A))'   


! Time control cards

 
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'TMSTRT, TMEND,YEAR'
  WRITE (CON,'(2(F14.4,","),I4,",")')  TMSTRT,   TMEND,    YEAR    
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'NDLT,DLTMIN, DLTINTER'
  WRITE (CON,'((I8,","),F12.4,",",A,",")')  NDLT,     DLTMIN, ADJUSTL(DLTINTER)
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'DLTD,DLTD,DLTD,DLTD,DLTD,DLTD,DLTD,DLTD,DLTD,DLTD,'
  WRITE (CON,'(<NDLT>(F12.3,","),",")')  (DLTD(J), J =1,NDLT)
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'DLTMAX,DLTMAX,DLTMAX,DLTMAX,DLTMAX,DLTMAX,DLTMAX,DLTMAX,DLTMAX,'
  WRITE (CON,'(<NDLT>(F12.3,","),",")')  (DLTMAX(J), J =1,NDLT)
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'DLTF,DLTF,DLTF,DLTF,DLTF,DLTF,DLTF,DLTF,DLTF,DLTF,'
  WRITE (CON,'(<NDLT>(F12.3,","),",")')  (DLTF(J),   J =1,NDLT)
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10,'
  WRITE (CON,'(<NWB>(A,","))')(ADJUSTL(VISC(JW)), JW=1,NWB)
  WRITE (CON,'(<NWB>(A,","))')(ADJUSTL(CELC(JW)), JW=1,NWB)
  DLTADD='OFF'
  WRITE (CON,'(<NWB>(A,","))')(ADJUSTL(DLTADD(JW)), JW=1,NWB)


! Grid definition cards

  WRITE (CON,*)
  WRITE (CON,'(A,",")')'BR1,BR2,BR3,BR4,BR5,BR6,BR7,BR8,BR9,BR10,'
  WRITE (CON,'(<NBR>(I8,","))')  (US(JB),    JB=1,NBR)
  WRITE (CON,'(<NBR>(I8,","))')  (DS(JB),    JB=1,NBR)
  WRITE (CON,'(<NBR>(I8,","))')  (UHS(JB),   JB=1,NBR)
  WRITE (CON,'(<NBR>(I8,","))')  (DHS(JB),   JB=1,NBR)
  WRITE (CON,'(<NBR>(I8,","))')  (NL(JB),    JB=1,NBR)
  WRITE (CON,'(<NBR>(F12.5,","))')  (SLOPE(JB), JB=1,NBR)
  WRITE (CON,'(<NBR>(F12.5,","))')  (SLOPEC(JB),JB=1,NBR)
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10,'


  WRITE (CON,'(<NWB>(F12.5,","))')   (LAT(JW),                     JW=1,NWB)
  WRITE (CON,'(<NWB>(F12.5,","))')   (LONGIT(JW),                  JW=1,NWB)
  WRITE (CON,'(<NWB>(F12.5,","))')   (ELBOT(JW),                   JW=1,NWB)
  WRITE (CON,'(<NWB>(I8,","))')   (BS(JW),                      JW=1,NWB)
  WRITE (CON,'(<NWB>(I8,","))')   (BE(JW),                      JW=1,NWB)
  WRITE (CON,'(<NWB>(I8,","))')   (JBDN(JW),                    JW=1,NWB)
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10,'
  
  
! Initial condition cards

  WRITE (CON,'(<NWB>(F10.4,","))')      (T2I(JW),                              JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.4,","))')      (ICETHI(JW),                           JW=1,NWB)
  WRITE (CON,'(<NWB>(A,","))')      (ADJUSTL(WTYPEC(JW)),  JW=1,NWB)
  WRITE (CON,'(<NWB>(A,","))')      (ADJUSTL(GRIDC(JW)),   JW=1,NWB)
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10,'

  WRITE (CON,'(<NWB>(A,","))')            (ADJUSTL(VBC(JW)),              JW=1,NWB)
  WRITE (CON,'(<NWB>(A,","))')            (ADJUSTL(EBC(JW)),              JW=1,NWB)
  WRITE (CON,'(<NWB>(A,","))')            (ADJUSTL(MBC(JW)),              JW=1,NWB)  
  WRITE (CON,'(<NWB>(A,","))')            (ADJUSTL(PQC(JW)),              JW=1,NWB)
  WRITE (CON,'(<NWB>(A,","))')            (ADJUSTL(EVC(JW)),              JW=1,NWB)
  WRITE (CON,'(<NWB>(A,","))')            (ADJUSTL(PRC(JW)),              JW=1,NWB)
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10,'

  WRITE (CON,'(<NWB>(A,","))')            (ADJUSTL(WINDC(JW)),             JW=1,NWB)
  WRITE (CON,'(<NWB>(A,","))')            (ADJUSTL(QINC(JW)),              JW=1,NWB)
  WRITE (CON,'(<NWB>(A,","))')            (ADJUSTL(QOUTC(JW)),             JW=1,NWB)
  WRITE (CON,'(<NWB>(A,","))')            (ADJUSTL(HEATC(JW)),             JW=1,NWB)
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'BR1,BR2,BR3,BR4,BR5,BR6,BR7,BR8,BR9,BR10,'

  WRITE (CON,'(<NBR>(A,","))')            (ADJUSTL(QINIC(JB)),               JB=1,NBR)
  WRITE (CON,'(<NBR>(A,","))')            (ADJUSTL(DTRIC(JB)),               JB=1,NBR)
  WRITE (CON,'(<NBR>(A,","))')            (ADJUSTL(HDIC(JB)),                JB=1,NBR)
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10,'

  WRITE (CON,'(<NWB>(A,","))')      (ADJUSTL(SLHTC(JW)),                     JW=1,NWB)
  WRITE (CON,'(<NWB>(A,","))')      (ADJUSTL(SROC(JW)),                      JW=1,NWB)
  WRITE (CON,'(<NWB>(A,","))')      (ADJUSTL(RHEVC(JW)),                     JW=1,NWB)
  WRITE (CON,'(<NWB>(A,","))')      (ADJUSTL(METIC(JW)),                     JW=1,NWB)
  WRITE (CON,'(<NWB>(A,","))')      (ADJUSTL(FETCHC(JW)),                    JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.4,","))')      (AFW(JW),                       JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.4,","))')      (BFW(JW),                       JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.4,","))')      (CFW(JW),                       JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.5,","))')      (WINDH(JW),                     JW=1,NWB)
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10,'

  WRITE (CON,'(<NWB>(A,","))')      (ADJUSTL(ICEC(JW)),                      JW=1,NWB)
  WRITE (CON,'(<NWB>(A,","))')      (ADJUSTL(SLICEC(JW)),                    JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.4,","))')      (ALBEDO(JW),                    JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.4,","))')      (HWI(JW),                       JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.4,","))')      (BETAI(JW),                     JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.4,","))')      (GAMMAI(JW),                    JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.4,","))')      (ICEMIN(JW),                    JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.4,","))')      (ICET2(JW),                     JW=1,NWB)
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10,'

  WRITE (CON,'(<NWB>(A,","))')        (ADJUSTL(SLTRC(JW)),                   JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.4,","))')        (THETA(JW),                   JW=1,NWB)
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10,'

  WRITE (CON,'(<NWB>(F10.4,","))')       (AX(JW),                       JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.4,","))')       (DXI(JW),                      JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.4,","))')       (CBHE(JW),                     JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.4,","))')       (TSED(JW),                     JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.4,","))')       (FI(JW),                       JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.4,","))')       (TSEDF(JW),                    JW=1,NWB)
  WRITE (CON,'(<NWB>(A,","))')       (ADJUSTL(FRICC(JW)),                    JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.4,","))')       (Z0(JW),                       JW=1,NWB)
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10,'

  WRITE (CON,'(<NWB>(A,","))')      (ADJUSTL(AZC(JW)),    JW=1,NWB)        
  WRITE (CON,'(<NWB>(A,","))')      (ADJUSTL(AZSLC(JW)),  JW=1,NWB)        
  WRITE (CON,'(<NWB>(F10.4,","))')      (AZMAX(JW),  JW=1,NWB)          
  WRITE (CON,'(<NWB>(I5,","))')      (TKEBC(JW),  JW=1,NWB)          
  WRITE (CON,'(<NWB>(F10.4,","))')      (EROUGH(JW), JW=1,NWB)          
  WRITE (CON,'(<NWB>(F10.4,","))')      (ARODI(JW),  JW=1,NWB)          
  WRITE (CON,'(<NWB>(F10.4,","))')      (STRICK(JW), JW=1,NWB)          
  WRITE (CON,'(<NWB>(F10.4,","))')      (TKELATPRDCONST(JW),JW=1,NWB)   
  WRITE (CON,'(<NWB>(A,","))')      (ADJUSTL(IMPTKE(JW)),JW=1,NWB)        

! Inflow-outflow cards
  
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'BR1,BR2,BR3,BR4,BR5,BR6,BR7,BR8,BR9,BR10,'
  WRITE (CON,'(<NBR>(I4,","))')           (NSTR(JB),      JB=1,NBR)
  WRITE (CON,'(<NBR>(A,","))')            (ADJUSTL(DYNSTRUC(JB)),  JB=1,NBR)
  
  DO JS=1,NST
    WRITE (CON,'(<NBR>(A,","))')        (ADJUSTL(STRIC(JS,JB)),  JB=1,NBR)
  END DO
  IF(NST<5)THEN
      DO JS=NST+1,5
          WRITE(CON,*)',,,,,,,'
      ENDDO
  ENDIF
  
  DO JS=1,NST
    WRITE (CON,'(<NBR>(I4,","))')        (KTSWT(JS,JB),  JB=1,NBR)
  END DO
    IF(NST<5)THEN
      DO JS=NST+1,5
          WRITE(CON,*)',,,,,,,'
      ENDDO
  ENDIF
  DO JS=1,NST
    WRITE (CON,'(<NBR>(I4,","))')        (KBSWT(JS,JB),  JB=1,NBR)
  END DO
    IF(NST<5)THEN
      DO JS=NST+1,5
          WRITE(CON,*)',,,,,,,'
      ENDDO
  ENDIF
  DO JS=1,NST
    WRITE (CON,'(<NBR>(A,","))')        (ADJUSTL(SINKCT(JS,JB)), JB=1,NBR)
  END DO
    IF(NST<5)THEN
      DO JS=NST+1,5
          WRITE(CON,*)',,,,,,,'
      ENDDO
  ENDIF
  DO JS=1,NST
    WRITE (CON,'(<NBR>(F10.3,","))')        (ESTRT(JS,JB),  JB=1,NBR)
  END DO
    IF(NST<5)THEN
      DO JS=NST+1,5
          WRITE(CON,*)',,,,,,,'
      ENDDO
  ENDIF
  DO JS=1,NST
    WRITE (CON,'(<NBR>(F10.3,","))')        (WSTRT(JS,JB),  JB=1,NBR)
  END DO
    IF(NST<5)THEN
      DO JS=NST+1,5
          WRITE(CON,*)',,,,,,,'
      ENDDO
  ENDIF
  WRITE (CON,*)
  
return
    end SUBROUTINE CONTROL_FILE_CSV1
    !**************************************************************
 SUBROUTINE CONTROL_FILE_CSV2
  USE MAIN
  USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  use macrophytec; use porosityc; use zooplanktonc;   use initialvelocity  ;  use extra, ONLY:CONFNwrite
  IMPLICIT NONE
  
  INTEGER :: N !,NNDC

  WRITE (CON,'(A,",")')'PIPE1,PIPE2,PIPE3,PIPE4,PIPE5,PIPE6,PIPE7,PIPE8,PIPE9,PIPE10,'
  WRITE (CON,'(<NPI>(I5,","))')  (IUPI(JP),  JP=1,NPI)
  WRITE (CON,'(<NPI>(I5,","))')  (IDPI(JP),  JP=1,NPI)
  WRITE (CON,'(<NPI>(F10.4,","))')  (EUPI(JP),  JP=1,NPI)
  WRITE (CON,'(<NPI>(F10.4,","))')  (EDPI(JP),  JP=1,NPI)
  WRITE (CON,'(<NPI>(F10.4,","))')  (WPI(JP),   JP=1,NPI)
  WRITE (CON,'(<NPI>(F10.4,","))')  (DLXPI(JP), JP=1,NPI)
  WRITE (CON,'(<NPI>(F8.4,","))')  (FPI(JP),   JP=1,NPI)
  WRITE (CON,'(<NPI>(F8.4,","))')  (FMINPI(JP),JP=1,NPI)
  WRITE (CON,'(<NPI>(A,","))')  (ADJUSTL(LATPIC(JP)),JP=1,NPI)
  WRITE (CON,'(<NPI>(A,","))')  (ADJUSTL(DYNPIPE(JP)),JP=1,NPI) 
  
  WRITE (CON,'(<NPI>(A,","))')  (ADJUSTL(PUPIC(JP)),  JP=1,NPI)
  WRITE (CON,'(<NPI>(F8.3,","))')  (ETUPI(JP),  JP=1,NPI)
  WRITE (CON,'(<NPI>(F8.3,","))')  (EBUPI(JP),  JP=1,NPI)
  WRITE (CON,'(<NPI>(I5,","))')  (KTUPI(JP),  JP=1,NPI)
  WRITE (CON,'(<NPI>(I5,","))')  (KBUPI(JP),  JP=1,NPI)
  
  WRITE (CON,'(<NPI>(A,","))')  (ADJUSTL(PDPIC(JP)),  JP=1,NPI)
  WRITE (CON,'(<NPI>(F8.4,","))')  (ETDPI(JP),  JP=1,NPI)
  WRITE (CON,'(<NPI>(F8.4,","))')  (EBDPI(JP),  JP=1,NPI)
  WRITE (CON,'(<NPI>(I5,","))')  (KTDPI(JP),  JP=1,NPI)
  WRITE (CON,'(<NPI>(I5,","))')  (KBDPI(JP),  JP=1,NPI)
  
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'SP1,SP2,SP3,SP4,SP5,SP6,SP7,SP8,SP9,SP10,'
  WRITE (CON,'(<NSP>(I5,","))') (IUSP(JS),    JS=1,NSP)
  WRITE (CON,'(<NSP>(I5,","))') (IDSP(JS),    JS=1,NSP)
  WRITE (CON,'(<NSP>(F8.3,","))') (ESP(JS),     JS=1,NSP)
  WRITE (CON,'(<NSP>(F8.3,","))') (A1SP(JS),    JS=1,NSP)
  WRITE (CON,'(<NSP>(F8.3,","))') (B1SP(JS),    JS=1,NSP)
  WRITE (CON,'(<NSP>(F8.3,","))') (A2SP(JS),    JS=1,NSP)
  WRITE (CON,'(<NSP>(F8.3,","))') (B2SP(JS),    JS=1,NSP)
  WRITE (CON,'(<NSP>(A,","))') (ADJUSTL(LATSPC(JS)),  JS=1,NSP)
  
  WRITE (CON,'(<NSP>(A,","))') (ADJUSTL(PUSPC(JS)),   JS=1,NSP)
  WRITE (CON,'(<NSP>(F8.3,","))') (ETUSP(JS),   JS=1,NSP)
  WRITE (CON,'(<NSP>(F8.3,","))') (EBUSP(JS),   JS=1,NSP)
  WRITE (CON,'(<NSP>(I5,","))') (KTUSP(JS),   JS=1,NSP)
  WRITE (CON,'(<NSP>(I5,","))') (KBUSP(JS),   JS=1,NSP)
  
  WRITE (CON,'(<NSP>(A,","))') (ADJUSTL(PDSPC(JS)),   JS=1,NSP) 
  WRITE (CON,'(<NSP>(F8.3,","))') (ETDSP(JS),   JS=1,NSP)
  WRITE (CON,'(<NSP>(F8.3,","))') (EBDSP(JS),   JS=1,NSP)
  WRITE (CON,'(<NSP>(I5,","))') (KTDSP(JS),   JS=1,NSP)
  WRITE (CON,'(<NSP>(I5,","))') (KBDSP(JS),   JS=1,NSP)
  
  WRITE (CON,'(<NSP>(A,","))') (ADJUSTL(GASSPC(JS)),  JS=1,NSP)
  WRITE (CON,'(<NSP>(I5,","))') (EQSP(JS),    JS=1,NSP)
  WRITE (CON,'(<NSP>(F8.3,","))') (AGASSP(JS),  JS=1,NSP)
  WRITE (CON,'(<NSP>(F8.3,","))') (BGASSP(JS),  JS=1,NSP)
  WRITE (CON,'(<NSP>(F8.3,","))') (CGASSP(JS),  JS=1,NSP)
  
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'GATE1,GATE2,GTE3,GATE4,GATE5,GATE6,GATE7,GATE8,GATE9,GATE10,'
  WRITE (CON,'(<NGT>(I5,","))') (IUGT(JG),   JG=1,NGT)
  WRITE (CON,'(<NGT>(I5,","))') (IDGT(JG),   JG=1,NGT)
  WRITE (CON,'(<NGT>(F8.3,","))') (EGT(JG),    JG=1,NGT)
  WRITE (CON,'(<NGT>(F8.3,","))') (A1GT(JG),   JG=1,NGT)
  WRITE (CON,'(<NGT>(F8.3,","))') (B1GT(JG),   JG=1,NGT)
  WRITE (CON,'(<NGT>(F8.3,","))') (G1GT(JG),   JG=1,NGT)
  WRITE (CON,'(<NGT>(F8.3,","))') (A2GT(JG),   JG=1,NGT)
  WRITE (CON,'(<NGT>(F8.3,","))') (B2GT(JG),   JG=1,NGT)
  WRITE (CON,'(<NGT>(F8.3,","))') (G2GT(JG),   JG=1,NGT)
  WRITE (CON,'(<NGT>(A,","))') (ADJUSTL(LATGTC(JG)), JG=1,NGT)
  
  WRITE (CON,'(<NGT>(F8.3,","))') (GTA1(JG),   JG=1,NGT)  
  WRITE (CON,'(<NGT>(F8.3,","))') (GTB1(JG),   JG=1,NGT)  
  WRITE (CON,'(<NGT>(F8.3,","))') (GTA2(JG),   JG=1,NGT)  
  WRITE (CON,'(<NGT>(F8.3,","))') (GTB2(JG),   JG=1,NGT)  
  WRITE (CON,'(<NGT>(A,","))') (ADJUSTL(DYNGTC(JG)), JG=1,NGT)
  WRITE (CON,'(<NGT>(A,","))') (ADJUSTL(GTIC(JG)),   JG=1,NGT)
  WRITE (CON,'(<NGT>(A,","))') (ADJUSTL(PUGTC(JG)),  JG=1,NGT)
  WRITE (CON,'(<NGT>(F8.3,","))') (ETUGT(JG),  JG=1,NGT)
  WRITE (CON,'(<NGT>(F8.3,","))') (EBUGT(JG),  JG=1,NGT)
  WRITE (CON,'(<NGT>(I5,","))') (KTUGT(JG),  JG=1,NGT)
  WRITE (CON,'(<NGT>(I5,","))') (KBUGT(JG),  JG=1,NGT)
  WRITE (CON,'(<NGT>(A,","))') (ADJUSTL(PDGTC(JG)),  JG=1,NGT)
  WRITE (CON,'(<NGT>(F8.3,","))') (ETDGT(JG),  JG=1,NGT)
  WRITE (CON,'(<NGT>(F8.3,","))') (EBDGT(JG),  JG=1,NGT)
  WRITE (CON,'(<NGT>(I5,","))') (KTDGT(JG),  JG=1,NGT)
  WRITE (CON,'(<NGT>(I5,","))') (KBDGT(JG),  JG=1,NGT)
  
  WRITE (CON,'(<NGT>(A,","))') (ADJUSTL(GASGTC(JG)), JG=1,NGT)
  WRITE (CON,'(<NGT>(I5,","))') (EQGT(JG),   JG=1,NGT)  
  WRITE (CON,'(<NGT>(F8.3,","))') (AGASGT(JG), JG=1,NGT)  
  WRITE (CON,'(<NGT>(F8.3,","))') (BGASGT(JG), JG=1,NGT)  
  WRITE (CON,'(<NGT>(F8.3,","))') (CGASGT(JG), JG=1,NGT)  
  
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'Pump1,Pump2,Pump3,Pump4,Pump5,Pump6,Pump7,Pump8,Pump9,Pump10,'
  
  WRITE (CON,'(<NPU>(I5,","))') (IUPU(JP),   JP=1,NPU)
  WRITE (CON,'(<NPU>(I5,","))') (IDPU(JP),   JP=1,NPU)
  WRITE (CON,'(<NPU>(F8.3,","))') (EPU(JP),    JP=1,NPU)
  WRITE (CON,'(<NPU>(F8.3,","))') (STRTPU(JP), JP=1,NPU)
  WRITE (CON,'(<NPU>(F8.3,","))') (ENDPU(JP),  JP=1,NPU)
  WRITE (CON,'(<NPU>(F8.3,","))') (EONPU(JP),  JP=1,NPU)
  WRITE (CON,'(<NPU>(F8.3,","))') (EOFFPU(JP), JP=1,NPU)
  WRITE (CON,'(<NPU>(F8.3,","))') (QPU(JP),    JP=1,NPU)
  WRITE (CON,'(<NPU>(A,","))') (ADJUSTL(LATPUC(JP)), JP=1,NPU)
  WRITE (CON,'(<NPU>(A,","))') (ADJUSTL(DYNPUMP(JP)),JP=1,NPU)

! Output control cards (excluding constituents)

  WRITE (CON,'(<NPU>(A,","))') (ADJUSTL(PPUC(JP)),   JP=1,NPU)
  WRITE (CON,'(<NPU>(F8.3,","))') (ETPU(JP),   JP=1,NPU)
  WRITE (CON,'(<NPU>(F8.3,","))') (EBPU(JP),   JP=1,NPU)
  WRITE (CON,'(<NPU>(I5,","))') (KTPU(JP),   JP=1,NPU)
  WRITE (CON,'(<NPU>(I5,","))') (KBPU(JP),   JP=1,NPU)
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'IW1,IW2,IW3,IW4,IW5,IW6,IW7,IW8,IW9,IW10,'

  WRITE (CON,'(<NIW>(I5,","))')        (IWR(JW),    JW=1,NIW)
  WRITE (CON,'(<NIW>(I5,","))')        (KTWR(JW),   JW=1,NIW)               
  WRITE (CON,'(<NIW>(I5,","))')        (KBWR(JW),   JW=1,NIW)   
  WRITE (CON,*)
  
return
    end SUBROUTINE CONTROL_FILE_CSV2
!***************************************    
    
    !**************************************************************
 SUBROUTINE WRITE_CSV
  USE MAIN
  USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  use macrophytec; use porosityc; use zooplanktonc;   use initialvelocity  ;  use extra, ONLY:CONFNwrite
  IMPLICIT NONE
  
  INTEGER :: N,NTFORMAT,NNEW,JJ  !NNDC,
  CHARACTER(43) :: TITLE43

  WRITE (CON,'(A,",")')'WD1,WD2,WD3,WD4,WD5,WD6,WD7,WD8,WD9,WD10,'
  
  WRITE (CON,'(<NWD>(A,","))')         (ADJUSTL(WDIC(JW)),   JW=1,NWD)
  WRITE (CON,'(<NWD>(I5,","))')         (IWD(JW),    JW=1,NWD)
  WRITE (CON,'(<NWD>(F10.3,","))')      (EWD(JW),    JW=1,NWD)
  WRITE (CON,'(<NWD>(I5,","))')         (KTWD(JW),   JW=1,NWD)
  WRITE (CON,'(<NWD>(I5,","))')         (KBWD(JW),   JW=1,NWD)  
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'TR1,TR2,TR3,TR4,TR5,TR6,TR7,TR8,TR9,TR10,'

  WRITE (CON,'(<NTR>(A,","))')        (ADJUSTL(TRC(JT)),   JT=1,NTR)
  WRITE (CON,'(<NTR>(A,","))')        (ADJUSTL(TRIC(JT)),   JT=1,NTR)
  WRITE (CON,'(<NTR>(I5,","))')       (ITR(JT),    JT=1,NTR)
  WRITE (CON,'(<NTR>(F10.3,","))')    (ELTRT(JT),  JT=1,NTR)
  WRITE (CON,'(<NTR>(F10.3,","))')    (ELTRB(JT),  JT=1,NTR)
  WRITE (CON,'(<NTR>(A,A,A,","))')        (('"',QTRFN(JT),'"'),  JT=1,NTR)
  WRITE (CON,'(<NTR>(A,A,A,","))')        (('"',TTRFN(JT),'"'),  JT=1,NTR)
  WRITE (CON,'(<NTR>(A,A,A,","))')        (('"',CTRFN(JT),'"'),  JT=1,NTR)  
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'BR1,BR2,BR3,BR4,BR5,BR6,BR7,BR8,BR9,BR10,'
  
  WRITE (CON,'(<NBR>(A,","))')      (ADJUSTL(DTRC(JB)),   JB=1,NBR)
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'HNAME,FMTH,HMULT,HPRWBC1,HPRWBC2,HPRWBC3,HPRWBC4,HPRWBC5,HPRWBC6,HPRWBC7,HPRWBC8,HPRWBC9,HPRWBC10,'

  DO JH=1,NHY
    WRITE (CON,'(A,A,A,",",A,",",F10.3,",",<NWB>(A,","))')   '"', HNAME(JH),'"',  FMTH(JH),  HMULT(JH),(ADJUSTL(HPRWBC(JH,JW)),JW=1,NWB)    !'(""",A,"",",""",A,"",",F10.3,",",<NWB>(A,","))'
  END DO
  WRITE (CON,*)
  WRITE (CON,*)'SNP'
  
  WRITE (CON,'(A,",")')  ADJUSTL(SNPC(1))
  WRITE (CON,*)  NSNP(1)
      
  WRITE (CON,'(<NSNP(1)>(F12.4,","))')  (SNPD(J,1),J=1,NSNP(1))
  WRITE (CON,'(<NSNP(1)>(F12.4,","))')  (SNPF(J,1),J=1,NSNP(1))

  WRITE (CON,*)
  WRITE (CON,'(A,",")')'SCR'
  WRITE (CON,'(A,",")')         ADJUSTL(SCRC(1))    
  WRITE (CON,*)    NSCR(1)
  WRITE (CON,'(<NSCR(1)>(F12.4,","))')   (SCRD(J,1),J=1,NSCR(1))
  WRITE (CON,'(<NSCR(1)>(F12.4,","))')   (SCRF(J,1),J=1,NSCR(1))
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'PRFC'  
   
  WRITE (CON,'(A,",")')    ADJUSTL(PRFC(1))
  N=NIPRF(1)
  DO JW=2,NWB
  DO J=1,NIPRF(JW)
      N=N+1
      IPRF(N,1)=IPRF(J,JW)      
  ENDDO
  ENDDO
  NIPRF(1)=N
  
  WRITE (CON,*)  NPRF(1)       
  WRITE (CON,*)  NIPRF(1)     
  
  WRITE (CON,'(<NPRF(1)>(F12.4,","))')   (PRFD(J,1),J=1,NPRF(1))
  WRITE (CON,'(<NPRF(1)>(F12.4,","))')   (PRFF(J,1),J=1,NPRF(1))
  WRITE (CON,'(<NIPRF(1)>(I5,","))')   (IPRF(J,1),J=1,NIPRF(1))
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'SPR'  
  
  N=NISPR(1)
  DO JW=2,NWB
  DO J=1,NISPR(JW)
      N=N+1
      ISPR(N,1)=ISPR(J,JW)      
  ENDDO
  ENDDO
  NISPR(1)=N
   
  WRITE (CON,'(A,",")')      ADJUSTL(SPRC(1))         
  WRITE (CON,*)      NSPR(1)         
  WRITE (CON,*)      NISPR(1)          
  WRITE (CON,'(<NSPR(1)>(F12.3,","))')      (SPRD(J,1),J=1,NSPR(1))    
  WRITE (CON,'(<NSPR(1)>(F12.3,","))')      (SPRF(J,1),J=1,NSPR(1))    
  WRITE (CON,'(<NISPR(1)>(I5,","))')        (ISPR(J,1),J=1,NISPR(1))   
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'W2L'  
  
  WRITE (CON,'(A,",")')    ADJUSTL(VPLC(1))
  WRITE (CON,*)    NVPL(1)
  WRITE (CON,'(<NVPL(1)>(F12.3,","))')   (VPLD(J,1), J=1,NVPL(1))
  WRITE (CON,'(<NVPL(1)>(F12.3,","))')   (VPLF(J,1), J=1,NVPL(1))
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'CPL'  
  
  WRITE (CON,'(A,",")')      ADJUSTL(CPLC(1))
  WRITE (CON,*)      NCPL(1)
  WRITE (CON,'(A,",")')      ADJUSTL(TECPLOT(1))

  WRITE (CON,'(<NCPL(1)>(F12.3,","))')    (CPLD(J,1), J=1,NCPL(1))
  WRITE (CON,'(<NCPL(1)>(F12.3,","))')    (CPLF(J,1), J=1,NCPL(1))
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'FLUX' 
    
  WRITE (CON,'(A,",")') ADJUSTL(FLXC(1))
  WRITE (CON,*)        NFLX(1)

  WRITE (CON,'(<NFLX(1)>(F12.3,","))')       (FLXD(J,1), J=1,NFLX(1))
  WRITE (CON,'(<NFLX(1)>(F12.3,","))')       (FLXF(J,1), J=1,NFLX(1))
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'TSR' 
  
  WRITE (CON,'(A,",")')   ADJUSTL(TSRC)
  WRITE (CON,*)           NTSR
  WRITE (CON,*)          NIKTSR
  WRITE (CON,'(A,",")')  TSRFN
  
  WRITE (CON,'(<NTSR>(F12.3,","))')        (TSRD(J), J=1,NTSR)
  WRITE (CON,'(<NTSR>(F12.3,","))')        (TSRF(J), J=1,NTSR)
  WRITE (CON,'(<NIKTSR>(I5,","))')      (ITSR(J), J=1,NIKTSR)
  WRITE (CON,'(<NIKTSR>(F12.3,","))')   (ETSR(J), J=1,NIKTSR)
  WRITE (CON,*)
  
  
  WRITE (CON,'(A,",")')'WLEVEL' 
  WRITE (CON,'(A)')'OFF'
  WRITE (CON,*)'14.0'
  WRITE (CON,*)
  
  WRITE (CON,'(A,",")')'FLOWBAL' 
  WRITE (CON,'(A)')'ON'
  WRITE (CON,*)'14.0'
  WRITE (CON,*)
  
  WRITE (CON,'(A,",")')'NPBAL' 
  WRITE (CON,'(A)')'OFF'
  WRITE (CON,*)'14.0'
  WRITE (CON,*)
  
  WRITE (CON,'(A,",")')'WDO' 
   
  WRITE (CON,'(A,",")')           ADJUSTL(WDOC)   
  WRITE (CON,*)           NWDO 
  WRITE (CON,*)           NIWDO
  WRITE (CON,'(A,",")')           WDOFN

  WRITE (CON,'(<NWDO>(F12.3,","))')       (WDOD(J), J=1,NWDO)
  WRITE (CON,'(<NWDO>(F12.3,","))')       (WDOF(J), J=1,NWDO)
  WRITE (CON,'(<NIWDO>(I5,","))')         (IWDO(J), J=1,NIWDO)
  WRITE (CON,*)
  WRITE (CON,'(A,",")') 'RESTART'

  WRITE (CON,'(A,",")')        ADJUSTL(RSOC)
  WRITE (CON,*)         NRSO
  WRITE (CON,'(A,",")')         ADJUSTL(RSIC)
  WRITE (CON,'(A,A,A,",")')   '"',RSIFN,'"'

  WRITE (CON,'(<NRSO>(F12.3,","))')           (RSOD(J), J=1,NRSO)
  WRITE (CON,'(<NRSO>(F12.3,","))')           (RSOF(J), J=1,NRSO)
  WRITE (CON,*)
  WRITE (CON,'(A,",")') 'CCC,LIMC,CUF,CO2PPM,CO2YRLY'
  
  ! Constituent control cards

  WRITE (CON,'(A,",",A,",",I4,",",F10.3,",",A,",")')          ADJUSTL(CCC), ADJUSTL(LIMC), CUF,CO2PPM,ADJUSTL(CO2YRLY)
  WRITE (CON,*)
  
  WRITE (CON,'(A)')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10,'
  WRITE (CON,'(<NWB>(A,","))') ('OFF', JW=1,NWB)   !ATM_DEPOSITIONC(JW)
  WRITE (CON,'(<NWB>(A,","))') ('ON', JW=1,NWB)     !ATM_DEPOSITION_INTERPOLATION(JW)

  
  WRITE (CON,'(A)')',,,,,Repeat columns as necessary - EVEN IF NO TRIBUTARIES INCLUDE A COLUMN FOR CTRTRC1 (It can be blanks)'
    IF(NTR==0)THEN
        NTFORMAT=1
    ELSE
        NTFORMAT=NTR
    ENDIF  
    WRITE (CON,'(A,A,A,A,A,<NWB>(A),<NWB>(A),<NWB>(A),<NBR>(A),<NTFORMAT>(A),<NBR>(A),<NBR>(A))')'CNAME2-short name,','CNAME-long name,','CAC-active,','FMTC-Fortran output,', 'CMULT-output multiplier,', ('C2I(JW)-initial conc,', JW=1,NWB),('CPRWBC(JW)-print,', JW=1,NWB), ('C_ATM_DEPOSITION(JW),', JW=1,NWB),('CINBRC(JB)-inflow c,', JB=1,NBR),('CTRTRC(JT)-trib c,', JT=1,NTFORMAT), ('CDTBRC(JB)-distributed c,', JB=1,NBR), ('CPRBRC(JB)-precip c,', JB=1,NBR)  

  NNEW=1+NGC+NSS  
  DO JC=1,NCT
      
  IF(cname2(jc)=='FE      ' .or. cname2(jc)=='      FE')cycle
  WRITE (CON,'(A,",",A,A,A,",",A,",",A,",",F14.5,",",<NWB>(F12.5,","),<NWB>(A,","),<NWB>(A,","),<NBR>(A,","),<NTFORMAT>(A,","),<NBR>(A,","),<NBR>(A,","))')CNAME2(JC),'"',CNAME(JC),'"',ADJUSTL(CAC(JC)), FMTC(JC), CMULT(JC), (C2I(JC,JW), JW=1,NWB),(ADJUSTL(CPRWBC(JC,JW)), JW=1,NWB),('OFF',JW=1,NWB), (ADJUSTL(CINBRC(JC,JJB)), JJB=1,NBR),(ADJUSTL(CTRTRC(JC,JT)), JT=1,NTFORMAT), (ADJUSTL(CDTBRC(JC,JJB)), JJB=1,NBR), (ADJUSTL(CPRBRC(JC,JJB)), JJB=1,NBR)  
  ! ADD NEW STATE VARIABLES

  IF(JC==NNEW)THEN
  DO JJ=1,11
  WRITE (CON,'(A,",",A,A,A,",",A,",",A,",",F14.5,",",<NWB>(F12.5,","),<NWB>(A,","),<NWB>(A,","),<NBR>(A,","),<NTFORMAT>(A,","),<NBR>(A,","),<NBR>(A,","))')NAME1(JJ),'"',NAME2(JJ),'"','OFF','(F10.3)',1.0, (0.0, JW=1,NWB),('OFF', JW=1,NWB),('OFF',JW=1,NWB), ('OFF', JJB=1,NBR),('OFF', JT=1,NTFORMAT), ('OFF', JJB=1,NBR), ('OFF', JJB=1,NBR) 
  ENDDO
  ENDIF
  
  ENDDO

  WRITE (CON,*)
  WRITE (CON,'(A,",")')'CDNAME2,CDNAME,FMTCD,CDMULT,CDWBC1,CDWBC2,CDWBC3,CDWBC4,CDWBC5,CDWBC6,CDWBC7,CDWBC8,CDWBC9'
  
  DO J=23,17,-1
     CDNAME(J+3)=CDNAME(J);CDNAME2(J+3)=CDNAME2(J);FMTCD(J+3)=FMTCD(J);CDWBC(J+3,:)=CDWBC(J,:);CDMULT(J+3)=CDMULT(J)
  ENDDO
  CDNAME2(18)='TDG';CDNAME(18)='Total dissolved gas, % ';CDNAME2(19)='Turbidity';CDNAME(19)='Turbidity, NTU';CDMULT(19)=1.0;CDMULT(18)=1.0;CDWBC(18,:)='    OFF';CDWBC(19,:)='    OFF'
  DO J=16,9,-1
      CDNAME(J+1)=CDNAME(J);CDNAME2(J+1)=CDNAME2(J);FMTCD(J+1)=FMTCD(J);CDWBC(J+1,:)=CDWBC(J,:);CDMULT(J+1)=CDMULT(J)
  ENDDO
  CDNAME2(9)='NH3';CDNAME(9)='Unionized ammonia, g/m3 as N';CDWBC(19,:)='    OFF';CDMULT(9)=1.0;FMTCD(9)=FMTCD(8)
  CDNAME2(27)='Secchi';CDNAME(27)='Secchi disk depth, m';CDWBC(27,:)='    OFF';CDMULT(27)=1.0;FMTCD(27)=FMTCD(26)
  
  DO JD=1,NDC
  WRITE (CON,'(A,",",A,A,A,",",A,",",F14.5,",",<NWB>(A,","))')  CDNAME2(JD),'"',CDNAME(JD),'"',FMTCD(JD),CDMULT(JD),(ADJUSTL(CDWBC(JD,JW)), JW=1,NWB)
  ENDDO
  WRITE (CON,*)
  WRITE (CON,'(A,",")') 'KFNAME2,CFWBC1,CFWBC2,CFWBC3,CFWBC4,CFWBC5,CFWBC6,CFWBC7,CFWBC8'

  
  DO J=42,73
      KFNAME2(J-1)=KFNAME2(J);KFWBC(J-1,:)=KFWBC(J,:)
  ENDDO
  DO J=39,27,-1
      KFNAME2(J+1)=KFNAME2(J);KFWBC(J+1,:)=KFWBC(J,:)
  ENDDO
 KFNAME2(27)='NH3GAS';KFWBC(J,:)='     OFF'
  
  DO JF=1,72
    WRITE (CON,'(A,",",<NWB>(A,","))')   KFNAME2(JF),(ADJUSTL(KFWBC(JF,JW)),  JW=1,NWB)
  END DO

! Kinetics coefficients
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10'
 

  WRITE (CON,'(<NWB>(F10.3,","))')     (EXH2O(JW),  JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')     (EXSS(JW),   JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')     (EXOM(JW),   JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')     (BETA(JW),   JW=1,NWB)
  WRITE (CON,'(<NWB>(A,","))')     (ADJUSTL(EXC(JW)),    JW=1,NWB)
  WRITE (CON,'(<NWB>(A,","))')     (ADJUSTL(EXIC(JW)),   JW=1,NWB)
 
  WRITE (CON,*)
  WRITE (CON,'(A,",")') 'EXA1,EXA2,EXA3,EXA4,EXA5,EXA6,EXA7,EXA8,EXA9,EXA10'

  WRITE (CON,'(<NAL>(F10.3,","))')         (EXA(JA),  JA=1,NAL)
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'EXZ1,EXZ2,EXZ3,EXZ4,EXZ5,EXZ6,EXZ7,EXZ8,EXZ9,EXZ10' 

  WRITE (CON,'(<NZPT>(F10.3,","))')         (EXZ(JZ),  JZ=1,NZPT)  
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'EXM1,EXM2,EXM3,EXM4,EXM5,EXM6,EXM7' 
  
  WRITE (CON,'(<NMCT>(F10.3,","))')         (EXM(JM),   JM=1,NMCT)  
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'CG1,CG2,CG3,CG4,CG5,CG6,CG7,CG8,CG9,CG10,CG11,CG12' 

  WRITE (CON,'(<NGC>(F10.3,","))')         (CGQ10(JG),  JG=1,NGC)
  WRITE (CON,'(<NGC>(F10.3,","))')         (CG0DK(JG),  JG=1,NGC)
  WRITE (CON,'(<NGC>(F10.3,","))')         (CG1DK(JG),  JG=1,NGC) 
  WRITE (CON,'(<NGC>(F10.3,","))')         (CGS(JG),    JG=1,NGC) 
  WRITE (CON,'(<NGC>(F10.3,","))')         (CGLDK(JG),  JG=1,NGC) 
  WRITE (CON,'(<NGC>(F10.3,","))')         (CGKLF(JG),  JG=1,NGC) 
  WRITE (CON,'(<NGC>(F10.3,","))')         (CGCS(JG),   JG=1,NGC) 
  WRITE (CON,'(<NGC>(F10.3,","))')         (0.0,   JG=1,NGC)    ! CGR(JG)
 
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'SS1,SS2,SS3,SS4,SS5,SS6,SS7,SS8'  
  
  WRITE (CON,'(<NSS>(F10.3,","))')   (SSS(JS),    JS=1,NSS)     ! WRITE (CON,'(//(8X,F8.0,A8,2F8.0,I8))') (SSS(JS), SEDRC(JS),  TAUCR(JS),  SSFLOC(JS), FLOCEQN(JS),            JS=1,NSS) !SR 04/21/13
  WRITE (CON,'(<NSS>(A,","))')      (ADJUSTL(SEDRC(JS)),  JS=1,NSS)     ! WRITE (CON,'(//(8X,F8.0,A8,2F8.0,I8))') (SSS(JS), SEDRC(JS),  TAUCR(JS),  SSFLOC(JS), FLOCEQN(JS),            JS=1,NSS) !SR 04/21/13
  WRITE (CON,'(<NSS>(F10.3,","))')   (TAUCR(JS),  JS=1,NSS)     ! WRITE (CON,'(//(8X,F8.0,A8,2F8.0,I8))') (SSS(JS), SEDRC(JS),  TAUCR(JS),  SSFLOC(JS), FLOCEQN(JS),            JS=1,NSS) !SR 04/21/13
  WRITE (CON,'(<NSS>(F10.3,","))')   (0.0,  JS=1,NSS)   ! SSCS IF =-1 THEN THIS IS MATURE FINE TAILINGS
  WRITE (CON,*)
    
  WRITE (CON,'(A,",")')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10'
  WRITE (CON,'(<NWB>(F10.3,","))')     (1.04,  JW=1,NWB)    ! BACTQ10(JW)
    XX=0.001
  WRITE (CON,'(<NWB>(F10.3,","))')      (XX,  JW=1,NWB)    ! BACT1DK(JW)
  WRITE (CON,'(<NWB>(F10.3,","))')      (XX,  JW=1,NWB) ! (BACTLDK(JW),  JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')      (XX,  JW=1,NWB) ! (BACTS(JW),    JW=1,NWB)
  WRITE (CON,*)
  
  WRITE (CON,'(A,",")')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10'
  WRITE (CON,'(<NWB>(F10.3,","))')      (XX,  JW=1,NWB) !(H2SR(JW),    JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')      (1.04,  JW=1,NWB)  !(H2SQ10(JW),  JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')      (XX,  JW=1,NWB)  !(H2S1DK(JW),  JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')      (XX,  JW=1,NWB)  !(SO4R(JW),    JW=1,NWB)
  WRITE (CON,*)

  WRITE (CON,'(A,",")')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10'
  WRITE (CON,'(<NWB>(F10.3,","))')     (XX,  JW=1,NWB)   !(CH4R(JW),    JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')     (1.04,  JW=1,NWB)   !(CH4Q10(JW),  JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')     (0.1,  JW=1,NWB)   !(CH41DK(JW),  JW=1,NWB)
  WRITE (CON,*)
 
  WRITE (CON,'(A,",")')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10'
  WRITE (CON,'(<NWB>(F10.3,","))')     (XX,  JW=1,NWB)   !(FEIIR(JW),    JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')     (1.04,  JW=1,NWB)   !(KFE_OXID(JW),  JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')     (XX,  JW=1,NWB)  ! (KFE_RED(JW),  JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')     (XX,  JW=1,NWB)  ! (KFEOOH_HalfSat(JW),    JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')     (XX,  JW=1,NWB)  ! (FeSetVel(JW),  JW=1,NWB)
  WRITE (CON,*)

  WRITE (CON,'(A,",")')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10'
  WRITE (CON,'(<NWB>(F10.3,","))')    (XX,  JW=1,NWB)   ! (MNIIR(JW),    JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')    (1.04,  JW=1,NWB)   ! (KMN_OXID(JW),  JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')    (XX,  JW=1,NWB)   ! (KMN_RED(JW),  JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')    (XX,  JW=1,NWB)   ! (KMNO2_HalfSat(JW),    JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')    (XX,  JW=1,NWB)  ! (MNSetVel(JW),  JW=1,NWB)  
  
  
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'ALG1,ALG2,ALG3,ALG4,ALG5,ALG6,ALG7,ALG8'  
  
  WRITE (CON,'(<NAL>(F10.3,","))') (AG(JA),          JA=1,NAL)
  WRITE (CON,'(<NAL>(F10.3,","))') (AR(JA),          JA=1,NAL)
  WRITE (CON,'(<NAL>(F10.3,","))') (AE(JA),          JA=1,NAL)
  WRITE (CON,'(<NAL>(F10.3,","))') (AM(JA),          JA=1,NAL)
  WRITE (CON,'(<NAL>(F10.3,","))') (AS(JA),          JA=1,NAL)
  WRITE (CON,'(<NAL>(F10.4,","))') (AHSP(JA),        JA=1,NAL)
  WRITE (CON,'(<NAL>(F10.4,","))') (AHSN(JA),        JA=1,NAL)
  WRITE (CON,'(<NAL>(F10.4,","))') (AHSSI(JA),       JA=1,NAL)
  WRITE (CON,'(<NAL>(F10.3,","))') (ASAT(JA),        JA=1,NAL)

  WRITE (CON,'(<NAL>(F10.3,","))') (AT1(JA),         JA=1,NAL)
  WRITE (CON,'(<NAL>(F10.3,","))') (AT2(JA),         JA=1,NAL)
  WRITE (CON,'(<NAL>(F10.3,","))') (AT3(JA),         JA=1,NAL)
  WRITE (CON,'(<NAL>(F10.3,","))') (AT4(JA),         JA=1,NAL)
  WRITE (CON,'(<NAL>(F10.3,","))') (AK1(JA),         JA=1,NAL)
  WRITE (CON,'(<NAL>(F10.3,","))') (AK2(JA),         JA=1,NAL)
  WRITE (CON,'(<NAL>(F10.3,","))') (AK3(JA),         JA=1,NAL)
  WRITE (CON,'(<NAL>(F10.3,","))') (AK4(JA),         JA=1,NAL)

  WRITE (CON,'(<NAL>(F10.3,","))') (AP(JA),          JA=1,NAL)
  WRITE (CON,'(<NAL>(F10.3,","))') (AN(JA),          JA=1,NAL)
  WRITE (CON,'(<NAL>(F10.3,","))') (AC(JA),          JA=1,NAL)
  WRITE (CON,'(<NAL>(F10.3,","))') (ASI(JA),         JA=1,NAL)
  WRITE (CON,'(<NAL>(F10.3,","))') (ACHLA(JA),       JA=1,NAL)
  WRITE (CON,'(<NAL>(F10.3,","))') (APOM(JA),        JA=1,NAL)
  WRITE (CON,'(<NAL>(I3,","))') (ANEQN(JA),       JA=1,NAL)
  WRITE (CON,'(<NAL>(F10.3,","))') (ANPR(JA),        JA=1,NAL)

  WRITE (CON,'(<NAL>(F10.3,","))') (O2AR(JA),        JA=1,NAL)
  WRITE (CON,'(<NAL>(F10.3,","))') (O2AG(JA),        JA=1,NAL)
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10'  
   
  DO JE=1,NEPT
  WRITE (CON,'(<NWB>(A,","))')         (adjustl(EPIC(JW,JE)),  JW=1,NWB)
  WRITE (CON,'(<NWB>(A,","))')         (adjustl(EPIPRC(JW,JE)),JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.4,","))')   (EPICI(JW,JE), JW=1,NWB)
  ENDDO
  IF(NEPT<5)THEN
  DO JE=NEPT+1,5
  WRITE (CON,*) ',,,,,,,'  
  WRITE (CON,*) ',,,,,,,'  
  WRITE (CON,*) ',,,,,,,'  
  ENDDO
  ENDIF
      WRITE (CON,*)
  WRITE (CON,'(A,",")')'EP1,EP2,EP3,EP4,EP5,EP6,EP7,EP8'  
  
  WRITE (CON,'(<NEPT>(F10.3,","))') (EG(JE),           JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  WRITE (CON,'(<NEPT>(F10.3,","))') (ER(JE),           JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  WRITE (CON,'(<NEPT>(F10.3,","))') (EE(JE),           JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  WRITE (CON,'(<NEPT>(F10.3,","))') (EM(JE),           JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  WRITE (CON,'(<NEPT>(F10.3,","))') (EB(JE),           JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  WRITE (CON,'(<NEPT>(F10.3,","))') (EHSP(JE),         JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  WRITE (CON,'(<NEPT>(F10.3,","))') (EHSN(JE),         JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  WRITE (CON,'(<NEPT>(F10.3,","))') (EHSSI(JE),        JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13  
  
  WRITE (CON,'(<NEPT>(F10.3,","))') (ESAT(JE),         JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  WRITE (CON,'(<NEPT>(F10.3,","))') (EHS(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  WRITE (CON,'(<NEPT>(I5,","))')    (ENEQN(JE),        JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  WRITE (CON,'(<NEPT>(F10.3,","))') (ENPR(JE),         JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13  
  
  WRITE (CON,'(<NEPT>(F10.3,","))') (ET1(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  WRITE (CON,'(<NEPT>(F10.3,","))') (ET2(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  WRITE (CON,'(<NEPT>(F10.3,","))') (ET3(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  WRITE (CON,'(<NEPT>(F10.3,","))') (ET4(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  WRITE (CON,'(<NEPT>(F10.3,","))') (EK1(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  WRITE (CON,'(<NEPT>(F10.3,","))') (EK2(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  WRITE (CON,'(<NEPT>(F10.3,","))') (EK3(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  WRITE (CON,'(<NEPT>(F10.3,","))') (EK4(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  
  WRITE (CON,'(<NEPT>(F10.3,","))') (EP(JE),           JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  WRITE (CON,'(<NEPT>(F10.3,","))') (EN(JE),           JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  WRITE (CON,'(<NEPT>(F10.3,","))') (EC(JE),           JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  WRITE (CON,'(<NEPT>(F10.3,","))') (ESI(JE),          JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  WRITE (CON,'(<NEPT>(F10.3,","))') (ECHLA(JE),        JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  WRITE (CON,'(<NEPT>(F10.3,","))') (EPOM(JE),         JE=1,NEPT)        !JE=1,NEP)  SW 9/27/13
  WRITE (CON,'(<NEPT>(F10.3,","))') (O2ER(JE),         JE=1,NEPT)
  WRITE (CON,'(<NEPT>(F10.3,","))') (O2EG(JE),         JE=1,NEPT)

  WRITE (CON,*)
return
    end SUBROUTINE WRITE_CSV
    
     SUBROUTINE WRITE_CSV2
  USE MAIN
  USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  use macrophytec; use porosityc; use zooplanktonc;   use initialvelocity  ;  use extra, ONLY:CONFNwrite
  IMPLICIT NONE
  
  INTEGER :: N   

  WRITE (CON,'(A,",")')'Zoo1,Zoo2,Zoo3,Zoo4,Zoo5,Zoo6,Zoo7'
  
  WRITE (CON,'(<NZPT>(F10.3,","))')         (ZG(JZ),    JZ=1,NZPT)
  WRITE (CON,'(<NZPT>(F10.3,","))')         (ZR(JZ),    JZ=1,NZPT)
  WRITE (CON,'(<NZPT>(F10.3,","))')         (ZM(JZ),    JZ=1,NZPT)
  WRITE (CON,'(<NZPT>(F10.3,","))')         (ZEFF(JZ),  JZ=1,NZPT)
  WRITE (CON,'(<NZPT>(F10.3,","))')         (PREFP(JZ), JZ=1,NZPT)
  WRITE (CON,'(<NZPT>(F10.3,","))')         (ZOOMIN(JZ),JZ=1,NZPT)
  WRITE (CON,'(<NZPT>(F10.3,","))')         (ZS2P(JZ),  JZ=1,NZPT)
  WRITE (CON,'(<NZPT>(F10.3,","))')         (0.0,  JZ=1,NZPT)     ! ZS SETTLING RATE
  WRITE (CON,'(<NZPT>(F10.3,","))')         (ZT1(JZ),   JZ=1,NZPT)
  WRITE (CON,'(<NZPT>(F10.3,","))')         (ZT2(JZ),   JZ=1,NZPT)
  WRITE (CON,'(<NZPT>(F10.3,","))')         (ZT3(JZ),   JZ=1,NZPT)
  WRITE (CON,'(<NZPT>(F10.3,","))')         (ZT4(JZ),   JZ=1,NZPT)
  WRITE (CON,'(<NZPT>(F10.3,","))')         (ZK1(JZ),   JZ=1,NZPT)
  WRITE (CON,'(<NZPT>(F10.3,","))')         (ZK2(JZ),   JZ=1,NZPT)
  WRITE (CON,'(<NZPT>(F10.3,","))')         (ZK3(JZ),   JZ=1,NZPT)
  WRITE (CON,'(<NZPT>(F10.3,","))')         (ZK4(JZ),   JZ=1,NZPT) 
  
  WRITE (CON,'(<NZPT>(F10.3,","))')         (ZP(JZ),    JZ=1,NZPT)  
  WRITE (CON,'(<NZPT>(F10.3,","))')         (ZN(JZ),    JZ=1,NZPT)
  WRITE (CON,'(<NZPT>(F10.3,","))')         (ZC(JZ),    JZ=1,NZPT)

  WRITE (CON,'(<NZPT>(F10.3,","))')         (O2ZR(JZ),  JZ=1,NZPT)

  DO JA=1,NAL
    WRITE (CON,'(<NZPT>(F10.3,","))')    (PREFA(JA,JZ),     JZ=1,NZPT)
  END DO
  IF(NAL<5)THEN
      DO J=NAL+1,5
          WRITE (CON,*) ',,,,,,,'
      ENDDO
  ENDIF

  DO JZ=1,NZPT
    WRITE (CON,'(<NZPT>(F10.3,","))')    (PREFZ(JZ,JJZ),   JJZ=1,NZPT)       
  END DO
  IF(NZPT<5)THEN
      DO J=NZPT+1,5
          WRITE (CON,*) ',,,,,,,'
      ENDDO
  ENDIF
  
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10'   
  
  DO JM=1,NMCT
    WRITE (CON,'(<NWB>(A,","))')   (ADJUSTL(MACWBC(JW,JM)),  JW=1,NWB)
  END DO
  IF(NMCT<5)THEN
      DO J=NMCT+1,5
          WRITE (CON,*) ',,,,,,,'
      ENDDO
  ENDIF
  
  DO JM=1,NMCT
    WRITE (CON,'(<NWB>(A,","))')   (ADJUSTL(MPRWBC(JW,JM)),  JW=1,NWB)
  END DO
    IF(NMCT<5)THEN
      DO J=NMCT+1,5
          WRITE (CON,*) ',,,,,,,'
      ENDDO
  ENDIF

  DO JM=1,NMCT
    WRITE (CON,'(<NWB>(F10.3,","))')   (MACWBCI(JW,JM),  JW=1,NWB)
  END DO
    IF(NMCT<5)THEN
      DO J=NMCT+1,5
          WRITE (CON,*) ',,,,,,,'
      ENDDO
  ENDIF

  WRITE (CON,*)
  WRITE (CON,'(A,",")')'MacGroup1,MacGroup2,MacGroup3,MacGroup4,MacGroup5,MacGroup6,MacGroup7,MacGroup8,'   
  
  WRITE (CON,'(<NMCT>(F10.3,","))')         (MG(JM),     JM=1,NMCT)
  WRITE (CON,'(<NMCT>(F10.3,","))')         (MR(JM),     JM=1,NMCT)
  WRITE (CON,'(<NMCT>(F10.3,","))')         (MM(JM),     JM=1,NMCT)
  WRITE (CON,'(<NMCT>(F10.3,","))')         (MSAT(JM),   JM=1,NMCT)
  WRITE (CON,'(<NMCT>(F10.3,","))')         (MHSP(JM),   JM=1,NMCT)
  WRITE (CON,'(<NMCT>(F10.3,","))')         (MHSN(JM),   JM=1,NMCT)
  WRITE (CON,'(<NMCT>(F10.3,","))')         (MHSC(JM),   JM=1,NMCT)
  WRITE (CON,'(<NMCT>(F10.3,","))')         (MPOM(JM),   JM=1,NMCT)
  WRITE (CON,'(<NMCT>(F10.3,","))')         (LRPMAC(JM), JM=1,NMCT)
  
  WRITE (CON,'(<NMCT>(F10.3,","))')         (PSED(JM),   JM=1,NMCT)  
  WRITE (CON,'(<NMCT>(F10.3,","))')         (NSED(JM),   JM=1,NMCT)

  WRITE (CON,'(<NMCT>(F10.3,","))')         (MBMP(JM),   JM=1,NMCT)
  WRITE (CON,'(<NMCT>(F10.3,","))')         (MMAX(JM),   JM=1,NMCT)
  WRITE (CON,'(<NMCT>(F10.3,","))')         (CDDRAG(JM), JM=1,NMCT)  !CB 6/29/06
  WRITE (CON,'(<NMCT>(F10.3,","))')         (DWV(JM),    JM=1,NMCT)  !CB 6/29/06
  WRITE (CON,'(<NMCT>(F10.3,","))')         (DWSA(JM),   JM=1,NMCT)  !CB 6/29/06
  WRITE (CON,'(<NMCT>(F10.3,","))')         (ANORM(JM),  JM=1,NMCT)  !CB 6/29/06  
  
  WRITE (CON,'(<NMCT>(F10.3,","))')         (MT1(JM),    JM=1,NMCT)
  WRITE (CON,'(<NMCT>(F10.3,","))')         (MT2(JM),    JM=1,NMCT)
  WRITE (CON,'(<NMCT>(F10.3,","))')         (MT3(JM),    JM=1,NMCT)
  WRITE (CON,'(<NMCT>(F10.3,","))')         (MT4(JM),    JM=1,NMCT)
  WRITE (CON,'(<NMCT>(F10.3,","))')         (MK1(JM),    JM=1,NMCT)
  WRITE (CON,'(<NMCT>(F10.3,","))')         (MK2(JM),    JM=1,NMCT)
  WRITE (CON,'(<NMCT>(F10.3,","))')         (MK3(JM),    JM=1,NMCT)
  WRITE (CON,'(<NMCT>(F10.3,","))')         (MK4(JM),    JM=1,NMCT)
  
  WRITE (CON,'(<NMCT>(F10.3,","))')         (MP(JM),     JM=1,NMCT)
  WRITE (CON,'(<NMCT>(F10.3,","))')         (MN(JM),     JM=1,NMCT)
  WRITE (CON,'(<NMCT>(F10.3,","))')         (MC(JM),     JM=1,NMCT)
 
  WRITE (CON,'(<NMCT>(F10.3,","))')         (O2MR(JM),   JM=1,NMCT)
  WRITE (CON,'(<NMCT>(F10.3,","))')         (O2MG(JM),   JM=1,NMCT)

  WRITE (CON,*)
  WRITE (CON,'(A,",")')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10'   
  
  WRITE (CON,'(<NWB>(F10.3,","))')         (LDOMDK(JW),  JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (RDOMDK(JW),  JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (LRDDK(JW),   JW=1,NWB) 

  WRITE (CON,*)
  WRITE (CON,'(A,",")')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10'   
   
  WRITE (CON,'(<NWB>(F10.3,","))')         (LPOMDK(JW),   JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (RPOMDK(JW),   JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (LRPDK(JW),    JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (POMS(JW),     JW=1,NWB)
  
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10'   
  
  WRITE (CON,'(<NWB>(F10.3,","))')         (ORGP(JW),     JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (ORGN(JW),     JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (ORGC(JW),     JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (ORGSI(JW),    JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (O2OM(JW),     JW=1,NWB)

  WRITE (CON,'(<NWB>(F10.3,","))')         (OMT1(JW),     JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (OMT2(JW),     JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (OMK1(JW),     JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (OMK2(JW),     JW=1,NWB)
  
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10'
  WRITE (CON,'(<NWB>(F10.3,","))')         (1.1,     JW=1,NWB)  ! COEFFA_TURB
  WRITE (CON,'(<NWB>(F10.3,","))')         (0.05,     JW=1,NWB) ! COEFFB_TURB
  WRITE (CON,'(<NWB>(F10.3,","))')         (1.5,     JW=1,NWB)  ! SECC_PAR
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'BOD1,BOD2,BOD3,BOD4,BOD5,BOD6,BOD7,BOD8,BOD9,BOD10'   
  
  WRITE (CON,'(<NBOD>(F10.3,","))')         (KBOD(JB),     JB=1,NBOD)
  WRITE (CON,'(<NBOD>(F10.3,","))')         (TBOD(JB),     JB=1,NBOD)
  WRITE (CON,'(<NBOD>(F10.3,","))')         (RBOD(JB),     JB=1,NBOD)
  WRITE (CON,'(<NBOD>(F10.3,","))')         (CBODS(JB),    JB=1,NBOD)
  WRITE (CON,'(<NBOD>(F10.3,","))')         (BODP(JB),     JB=1,NBOD)
  WRITE (CON,'(<NBOD>(F10.3,","))')         (BODN(JB),     JB=1,NBOD)
  WRITE (CON,'(<NBOD>(F10.3,","))')         (BODC(JB),     JB=1,NBOD)
  
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10'   

  WRITE (CON,'(<NWB>(F10.3,","))')         (PO4R(JW),     JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (PARTP(JW),    JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (NH4R(JW),     JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (NH4DK(JW),    JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (NH4T1(JW),    JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (NH4T2(JW),    JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (NH4K1(JW),    JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (NH4K2(JW),    JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (168.,    JW=1,NWB)    ! KG_H2O FOR NH3 GAS TRANSFER
  WRITE (CON,'(<NWB>(F10.3,","))')         (O2NH4(JW),    JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (NO3DK(JW),    JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (NO3S(JW),     JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (FNO3SED(JW),  JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (NO3T1(JW),    JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (NO3T2(JW),    JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (NO3K1(JW),    JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (NO3K2(JW),    JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (DSIR(JW),     JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (PSIS(JW),     JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (PSIDK(JW),    JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')         (PARTSI(JW),   JW=1,NWB)
    
  !WRITE (CON,*)
  !WRITE (CON,'(A,",")')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10'   
  !
  !WRITE (CON,'(<NWB>(F10.3,","))')         (FER(JW),       JW=1,NWB)
  !WRITE (CON,'(<NWB>(F10.3,","))')         (FES(JW),       JW=1,NWB)
     
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10'   
  
  WRITE (CON,'(<NWB>(F10.3,","))')          (CO2R(JW),     JW=1,NWB)
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'O2LIMIT'   
  
  WRITE (CON,'(F10.3,",")')  KDO

  WRITE (CON,*)
  WRITE (CON,'(A,",")')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10'   

  WRITE (CON,'(<NWB>(A,","))')     (ADJUSTL(SEDCC(JW)),   JW=1,NWB)  
  WRITE (CON,'(<NWB>(A,","))')     (ADJUSTL(SEDPRC(JW)),  JW=1,NWB) 
  WRITE (CON,'(<NWB>(F10.3,","))')     (SEDCI(JW),   JW=1,NWB)  
  WRITE (CON,'(<NWB>(F10.3,","))')     (SDK(JW),     JW=1,NWB)     
  WRITE (CON,'(<NWB>(F10.3,","))')     (SEDS(JW),    JW=1,NWB)  
  WRITE (CON,'(<NWB>(F10.3,","))')     (FSOD(JW),    JW=1,NWB)  
  WRITE (CON,'(<NWB>(F10.3,","))')     (FSED(JW),    JW=1,NWB) 
  WRITE (CON,'(<NWB>(F10.3,","))')     (SEDB(JW),    JW=1,NWB)  
  WRITE (CON,'(<NWB>(A,","))')     (ADJUSTL(DYNSEDK(JW)), JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')     (SODT1(JW),   JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')     (SODT2(JW),   JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')     (SODK1(JW),   JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')     (SODK2(JW),   JW=1,NWB)
  WRITE (CON,*)
  WRITE (CON,'(<IMX>(I5,","))')(I,I=1,IMX)   
  
  WRITE (CON,'(<IMX>(F10.4,","))')  (SOD(I),  I=1,IMX)
  
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10'    
  
  WRITE (CON,'(<NWB>(A,","))')   (ADJUSTL(REAERC(JW)), JW=1,NWB)
  WRITE (CON,'(<NWB>(I3,","))')   (NEQN(JW),   JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')   (RCOEF1(JW), JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')   (RCOEF2(JW), JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')   (RCOEF3(JW), JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')   (RCOEF4(JW), JW=1,NWB)
  WRITE (CON,'(<NWB>(F10.3,","))')   (1.0, JW=1,NWB)   ! RCOEF5
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'FILE NAMES'    
  
! Input filenames

  WRITE (CON,'(A,A,A,",")')  '"',QWDFN,'"'
  WRITE (CON,'(A,A,A,",")')   '"',QGTFN,'"'
  WRITE (CON,'(A,A,A,",")')   '"',WSCFN,'"'
  WRITE (CON,'(A,A,A,",")')   '"',SHDFN,'"'
  WRITE (CON,'(A,A,A,",")')   '"',VPLFN(1),'"'
  
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'WB1,WB2,WB3,WB4,WB5,WB6,WB7,WB8,WB9,WB10'    
 
  WRITE (CON,'(<NWB>(A,A,A,","))') ( ('"',BTHFN(JW),'"'), JW=1,NWB)
  WRITE (CON,'(<NWB>(A,A,A,","))') ( ('"',METFN(JW),'"'), JW=1,NWB)
  WRITE (CON,'(<NWB>(A,A,A,","))') ( ('"',EXTFN(JW),'"'), JW=1,NWB)
  WRITE (CON,'(<NWB>(A,A,A,","))') ( ('"','ATM_DEP.CSV','"'), JW=1,NWB)
 
  WRITE (CON,'(<NWB>(A,A,A,","))') ( ('"',VPRFN(JW),'"'), JW=1,NWB)
  WRITE (CON,'(<NWB>(A,A,A,","))') ( ('"',LPRFN(JW),'"'), JW=1,NWB)

  ! Output filenames

  WRITE (CON,'(<NWB>(A,A,A,","))') ( ('"',SNPFN(JW),'"'), JW=1,NWB)
  WRITE (CON,'(<NWB>(A,A,A,","))') ( ('"',PRFFN(JW),'"'), JW=1,NWB)
  WRITE (CON,'(<NWB>(A,A,A,","))') ( ('"',CPLFN(JW),'"'), JW=1,NWB)
  WRITE (CON,'(<NWB>(A,A,A,","))') ( ('"',SPRFN(JW),'"'), JW=1,NWB)
  WRITE (CON,'(<NWB>(A,A,A,","))') ( ('"',FLXFN(JW),'"'), JW=1,NWB)
  
  WRITE (CON,*)
  WRITE (CON,'(A,",")')'BR1,BR2,BR3,BR4,BR5,BR6,BR7,BR8,BR9,BR10'    
 
  WRITE (CON,'(<NBR>(A,A,A,","))') ( ('"',QINFN(JB),'"'), JB=1,NBR)
  WRITE (CON,'(<NBR>(A,A,A,","))') ( ('"',TINFN(JB),'"'), JB=1,NBR)
  WRITE (CON,'(<NBR>(A,A,A,","))') ( ('"',CINFN(JB),'"'), JB=1,NBR)
  WRITE (CON,'(<NBR>(A,A,A,","))') ( ('"',QOTFN(JB),'"'), JB=1,NBR)
  WRITE (CON,'(<NBR>(A,A,A,","))') ( ('"',QDTFN(JB),'"'), JB=1,NBR)
  WRITE (CON,'(<NBR>(A,A,A,","))') ( ('"',TDTFN(JB),'"'), JB=1,NBR)
  WRITE (CON,'(<NBR>(A,A,A,","))') ( ('"',CDTFN(JB),'"'), JB=1,NBR)
  WRITE (CON,'(<NBR>(A,A,A,","))') ( ('"',PREFN(JB),'"'), JB=1,NBR)
  WRITE (CON,'(<NBR>(A,A,A,","))') ( ('"',TPRFN(JB),'"'), JB=1,NBR)
  WRITE (CON,'(<NBR>(A,A,A,","))') ( ('"',CPRFN(JB),'"'), JB=1,NBR)
  WRITE (CON,'(<NBR>(A,A,A,","))') ( ('"',EUHFN(JB),'"'), JB=1,NBR)
  WRITE (CON,'(<NBR>(A,A,A,","))') ( ('"',TUHFN(JB),'"'), JB=1,NBR)
  WRITE (CON,'(<NBR>(A,A,A,","))') ( ('"',CUHFN(JB),'"'), JB=1,NBR)
  WRITE (CON,'(<NBR>(A,A,A,","))') ( ('"',EDHFN(JB),'"'), JB=1,NBR)
  WRITE (CON,'(<NBR>(A,A,A,","))') ( ('"',TDHFN(JB),'"'), JB=1,NBR)
  WRITE (CON,'(<NBR>(A,A,A,","))') ( ('"',CDHFN(JB),'"'), JB=1,NBR)

  CLOSE (CON)
  
return
end SUBROUTINE WRITE_CSV2