**  upper columbia slough model - full model
** w2us.for model scott wells/psu/ce (503) 725-4276
* includes zooplankton (see ce-qual-r1), numerous additions
* all alterations are in lower case
************************************************************************
**                                                                    **
**                          CE-QUAL-W2                                **
**                                                                    **
**                       A Two-dimensional,                           **
**                Hydrodynamic - Water Quality Model                  **
**                                                                    **
**                         Version 2.02                               **
**                          March, 1993                               **
**                                                                    **
**                   Revisions by: Thomas M. Cole                     **
**                   Water Quality Modeling Group                     **
**                   U.S. Army Corps of Engineers                     **
**                   Waterways Experiment Station                     **
**                   Vicksburg, Mississippi 39180                     **
**                   Phone number: (601) 634-3283                     **
**                                                                    **
************************************************************************


************************************************************************
**             Task 1: Common Block Data Initialization               **
************************************************************************

      BLOCK DATA COMMON_DATA
c        PARAMETER (NCP=21)
        PARAMETER (NCP=22)
        INTEGER    SNP
        CHARACTER  CUNIT*6,    CNAME*16,   HNAME*24
        DIMENSION  CUNIT(NCP), CNAME(NCP), HNAME(3)
        COMMON /NAMESC/ HNAME, CNAME, CUNIT
        COMMON /SNPUC/  SNP
        DATA HNAME /'Horizontal velocity [U]',
     .              'Vertical velocity   [W]',
     .              'Temperature        [T2]'/
        DATA CNAME /'Tracer         ', 'Suspended solids',
     .              'Coliform',        'Dissolved solids',
     .              'Labile DOM',      'Refractory DOM',
     .              'Algae',           'Detritus',
     .              'Phosphorous',     'Ammonia',
     .              'Nitrate-Nitrite', 'Dissolved oxygen',
     .              'Sediment',        'Inorganic carbon',
     .              'Alkalinity',      'pH',
     .              'Carbon dioxide',  'Bicarbonate',
     .              'Carbonate',       'Iron',
     .              'CBOD',            'Zooplankton'/
        DATA CUNIT /8*'g/m^3', 3*'mg/m^3', 4*'g/m^3', ' ',
     .              3*'g/m^3',   'mg/m^3', 2*'g/m^3'/
        DATA SNP   /6/
      END

************************************************************************
**                  Task 2: Program Initialization                    **
************************************************************************

      PROGRAM  CE_QUAL_W2
      INCLUDE  'w2.inc'

***** Type declarations

      REAL    MINDLT
      REAL    LAT,    LONG
      REAL    JDAY,   JDMIN
      REAL    IT2,    IDX,     IICETH
      REAL    NH3T1,  NH3T2,   NO3T1,  NO3T2,  NO3K1,  NO3K2,  NH3K1,
     .        NH3K2,  LABDK,   LRFDK,  NH3DK,  NO3DK,  NH3D,   NO3D,
     .        NH3RM,  NO3RM,   NH3REL, KBOD
      REAL    NXTVD,  NXTMSN,  NXTMTS, NXTMPR, NXTMCP, NXTMVP, NXTMRS
      real    nxtep,nxgsmcdd
      REAL    ICETH,  ICESW,   ICETHU, ICETH1, ICETH2, ICEMIN, ICET2
      INTEGER CON,    BTH,     RSI,    RSO,    TSR,    PRF,    VPL,
     .        CPL,    VPR,     SNP,    ERR,    WRN
      INTEGER TSRDP,  SNPDP,   VPLDP,  CPLDP,  PRFDP,  RSODP,  WSCDP,
     .        DLTDP
      INTEGER UHS,    DHS,     US,     DS,     CUS
      INTEGER CN,     TRCN,    DTCN,   PRCN,   UHCN,   DHCN
      INTEGER FREQUK, YEAR,    GDAY,   SKTT
      LOGICAL SNAPSHOT,        PROFILE,        TIME_SERIES,
     .        VECTOR,          CONTOUR,        RESTART_IN,
     .        RESTART_OUT,     test,           test2, test3
      LOGICAL TRIBUTARIES,     DIST_TRIBS,     WITHDRAWALS,
     .        EVAPORATION,     PRECIPITATION, pumpon, notflowaug,
     .        pumpon4,pumpinc,pumpdec
      LOGICAL UP_FLOW,         DN_FLOW,        UP_HEAD,
     .        DN_HEAD,         UH_INTERNAL,    DH_INTERNAL,
     .        UH_EXTERNAL,     DH_EXTERNAL,    HEAD_BOUNDARY
      LOGICAL ISO_TEMP,        VERT_TEMP,      LONG_TEMP,
     .        ISO_CONC,        VERT_CONC,      LONG_CONC
      LOGICAL ADD_LAYER,       SUB_LAYER,      LONG_FORM,
     .        SHORT_FORM
      LOGICAL OPEN_FILES,      OPEN_VPR,       OPEN_LPR
      LOGICAL FRESH_WATER,     SALT_WATER,     SUSP_SOLIDS,
     .        CONSTITUENTS,    LIMITING_FACTOR
      LOGICAL NO_WIND,         NO_INFLOW,      NO_OUTFLOW,
     .        NO_HEAT
      LOGICAL WINTER,          ICE_IN,         ALLOW_ICE,
     .        ICE,             ICE_CALC,       DETAILED_ICE,
     .        TERM_BY_TERM
      LOGICAL INTERP_INFLOW,   INTERP_OUTFLOW, INTERP_MET,
     .        INTERP_DTRIBS,   INTERP_TRIBS,   INTERP_HEAD,
     .        INTERP_WITHDRWL, INTERPOLATE
      LOGICAL VOL_BALANCE,     PLACE_QIN,      END_RUN,
     .        BRANCH_FOUND,    UPWIND,         POINT_SINK,
     .        SHIFT_DEMAND,    LEAP_YEAR,      SEL_WITHDRAWAL,
     .        TRANSPORT,       MASS_BALANCE,   VIOLATION,
     .        UPDATE_KINETICS, PLACE_QTR
      LOGICAL WARNING_OPEN,    VOLUME_WARNING
      CHARACTER*1  EXT1
      CHARACTER*2  EXT2
      CHARACTER*3  CPRN,  HPRN,  GDCH,  EXT3
      CHARACTER*3  ACC,   INACC, TRACC, DTACC, PRACC
      CHARACTER*3  SNPC,  VPLC,  TSRC,  PRFC,  CPLC,  RSOC, RSIC
      CHARACTER*3  INFIC, OUTIC, TRIC,  DTRIC, WDIC,  HDIC, METIC
      CHARACTER*3  CCC,   PQC,   VBC,   EVC,   ICEC,  PRC,  LIMC,
     .             DTRC,  QINC,  QOUTC, WINDC, HEATC, SWC,  SDC,
     .             MBC,   PQTC
      CHARACTER*4  EXT4,  LFAC
      CHARACTER*5  FORM,  WTYPE, EXT5
      CHARACTER*6  CUNIT
      CHARACTER*7  BLANK, CONV,  LFPR
      CHARACTER*8  SLTR,  SINK,  SLICE, SLHEAT, stest, stest2, stest3
      CHARACTER*9  MONTH
      CHARACTER*12 RSOFN
      CHARACTER*13 CONVD1
      CHARACTER*16 CNAME
      CHARACTER*24 HNAME
      CHARACTER*72 TITLE, DUMMY
      CHARACTER*72 RSIFN, VPRFN, TSRFN, PRFFN, VPLFN, CPLFN, METFN,
     .             QINFN, TINFN, CINFN, QOTFN, QWDFN, QTRFN, TTRFN,
     .             CTRFN, QDTFN, TDTFN, CDTFN, PREFN, TPRFN, CPRFN,
     .             EUHFN, TUHFN, CUHFN, EDHFN, TDHFN, CDHFN, SNPFN,
     .             LPRFN, BTHFN
      DOUBLE PRECISION VOLSG,  VOLTG, VOLSBR, VOLTBR, VOLEV, VOLPR,
     .                 VOLTR,  VOLDT, VOLWD,  VOLUH,  VOLDH, VOLIN,
     .                 VOLOUT, DLVOL
      DOUBLE PRECISION ICMBR,  CMBR,  DLMBR,  CSSMBR
      DOUBLE PRECISION BHRHO,  GMA,   BTA,    A,      C,     D,
     .                 V,      F
      DOUBLE PRECISION Z1,     SZ1,   Z2,     ZMIN, el1, el2
      DOUBLE PRECISION BTAT,   GMAT,  AT,     CT,     DT,    VT
      DOUBLE PRECISION TADL,   TADV,  CADL,   CADV
      DOUBLE PRECISION AD1L,   AD2L,  AD3L,   AD1V,   AD2V,  AD3V
      DOUBLE PRECISION DX1,    DX2,   DX3,    ALFA
      DOUBLE PRECISION T1L,    T2L,   T3L,    C1L,    C2L,   C3L
      DOUBLE PRECISION DPT

***** Dimension statements

      dimension sman(imp),talgae(kmp,imp),zmu(kmp,imp),zrt(kmc,imc),
     .          zoormr(kmc,imc),zoormf(kmc,imc),zmt(kmc,imc),volb(nbp)
      DIMENSION ADMX(KMP,IMP), U(KMP,IMP),   SU(KMP,IMP),   W(KMP,IMP),
     .          SW(KMP,IMP),   AZ(KMP,IMP),  ADMZ(KMP,IMP), DZ(KMP,IMP),
     .          DX(KMP,IMP),   DM(KMP,IMP),  RHO(KMP,IMP),  P(KMP,IMP),
     .          HPG(KMP,IMP),  SB(KMP,IMP),  ST(KMP,IMP)
c      DIMENSION T1(KMP,IMP),   T2(KMP,IMP),  TSS(KMP,IMP),  QSS(KMP,IMP)
       DIMENSION  TSS(KMP,IMP),  QSS(KMP,IMP)
      DIMENSION B(KMP,IMP),    BB(KMP,IMP),  BR(KMP,IMP),   BH(KMP,IMP),
     .          BHR(KMP,IMP),  CONV(KMP,17), HSEG(KMP,IMP)
      DIMENSION HEIGHT(IMP),   HKT1(IMP),    HKT2(IMP),     AVHKT(IMP),
     .          BKT(IMP),      BHKT1(IMP),   BHKT2(IMP),    BHRKT1(IMP),
     .          BHRKT2(IMP),   DLX(IMP),     banksh(imp)
      DIMENSION KB(IMP),       KBMIN(IMP)
      DIMENSION BHRHO(IMP),    A(IMP),       C(IMP),        D(IMP),
     .          F(IMP),        V(IMP),       BTA(IMP),      GMA(IMP),
     .          Z1(IMP),       SZ1(IMP),     Z2(IMP)
      DIMENSION EV(IMP),       QDT(IMP),     AZMXKT(IMP),   QPR(IMP),
     .          IPRF(IMP),     RN(IMP),      INIT(IMP)
      DIMENSION ICETH(IMP),    ICE(IMP),     ICESW(IMP)
      DIMENSION SOD(IMP),      SODS(IMP)
      DIMENSION KTT(IMP),      SKTT(IMP)
      DIMENSION Q(IMP),        QC(IMP),      QSSUM(IMP)
      DIMENSION H(KMP),        AVH(KMP),     EL(KMP),       TVP(KMP),
     .          AZMAX(KMP),    qc1(imp),     qc2(imp)
      DIMENSION QWD(NWP),      TWD(NWP),     IWD(NWP),      KWD(NWP),
     .          JBWD(NWP)
      DIMENSION QTR(NTP),      TTR(NTP),     JBTR(NTP),     ITR(NTP),
     .          KTR(NTP),      QTRFN(NTP),   TTRFN(NTP),    CTRFN(NTP)
      DIMENSION CIC(NCP),      AVCOUT(NCP),  CPRN(NCP),     CNAME(NCP),
     .          CUNIT(NCP),    ACC(NCP),     CN(NCP),       INCN(NCP),
     .          TRCN(NCP),     DTCN(NCP),    PRCN(NCP),     UHCN(NCP),
     .          DHCN(NCP),     INACC(NCP),   TRACC(NCP),    DTACC(NCP),
     .          PRACC(NCP)
      dimension cicbr(ncp,nbp)
      DIMENSION KTOPSW(NBP),   KBOTSW(NBP)
      DIMENSION DTRC(NBP),     SWC(NBP)
      DIMENSION QSUM(NBP),     PHI0(NBP),    NOUT(NBP)
c      DIMENSION QSUM(NBP), PHI0(imp), NOUT(NBP) later for phi0 problem is
c adjustment of fetch length when each cell has its own orientation
      DIMENSION ELUH(NBP),     ELDH(NBP),    JBUH(NBP),     JBDH(NBP)
      DIMENSION US(NBP),       CUS(NBP),     DS(NBP),       UHS(NBP),
     .          DHS(NBP)
      DIMENSION QIN(NBP),      PR(NBP),      QPRBR(NBP),    QDTR(NBP),
     .          EVBR(NBP)
      DIMENSION TIN(NBP),      TPR(NBP),     TDTR(NBP)
      DIMENSION VOLEV(NBP),    VOLPR(NBP),   VOLTR(NBP),    VOLDT(NBP),
     .          VOLUH(NBP),    VOLDH(NBP),   VOLIN(NBP),    VOLOUT(NBP),
     .          VOLWD(NBP),    VOLSBR(NBP),  VOLTBR(NBP),   DLVOL(NBP)
      DIMENSION QINFN(NBP),    TINFN(NBP),   CINFN(NBP),    QDTFN(NBP),
     .          TDTFN(NBP),    CDTFN(NBP),   PREFN(NBP),    TPRFN(NBP),
     .          CPRFN(NBP),    EUHFN(NBP),   TUHFN(NBP),    CUHFN(NBP),
     .          EDHFN(NBP),    TDHFN(NBP),   CDHFN(NBP),    QOTFN(NBP)
      DIMENSION UP_FLOW(NBP),      DN_FLOW(NBP),     UP_HEAD(NBP),
     .          DN_HEAD(NBP),      UH_EXTERNAL(NBP), DH_EXTERNAL(NBP),
     .          UH_INTERNAL(NBP),  DH_INTERNAL(NBP), DIST_TRIBS(NBP)
      DIMENSION CIN(NCP,NBP),      CDTR(NCP,NBP),    CPR(NCP,NBP),
     .          CUH(KMC,NCP,NBP),  CDH(KMC,NCP,NBP), CTR(NCP,NTP),
     .          CWD(NWP,NCP)
      DIMENSION CMBR(NCP,NBP),     DLMBR(NCP,NBP),   ICMBR(NCP,NBP),
     .          CSSMBR(NCP,NBP)
      DIMENSION AGR(KMC,IMC),      ARR(KMC,IMC),     AMR(KMC,IMC),
     .          AER(KMC,IMC),      A1(KMC,IMC),      A2(KMC,IMC),
     .          A3(KMC,IMC),
     .          DETD(KMC,IMC),     SEDD(KMC,IMC),    ORGD(KMC,IMC),
     .          SO2D(KMC,IMC),     NH3D(KMC,IMC),    NO3D(KMC,IMC),
     .          CBODD(KMC,IMC),    OMRM(KMC,IMC),    NH3RM(KMC,IMC),
     .          NO3RM(KMC,IMC),    AGRMR(KMC,IMC),   AGRMF(KMC,IMC),
     .          SETIN(KMC,IMC),    SETOUT(KMC,IMC),  LFAC(KMC,IMC),
     .          LFPR(KMC,IMC)
      DIMENSION QUH1(KMP,NBP),     QUH2(KMP,NBP),    QDH1(KMP,NBP),
     .          QDH2(KMP,NBP),     TUH(KMP,NBP),     TDH(KMP,NBP),
     .          TSSUH1(KMP,NBP),   TSSUH2(KMP,NBP),  TSSDH1(KMP,NBP),
     .          TSSDH2(KMP,NBP),   AKBR(KMP,NBP),    QINF(KMP,NBP),
     .          QOUT(KMP,NBP),     KOUT(KMP,NBP),    QTRF(KMP,NTP)
      DIMENSION CONVD1(KMP,11)
      DIMENSION FETCHU(IMP,NBP),   FETCHD(IMP,NBP)
      DIMENSION ISO_CONC(NCP),     VERT_CONC(NCP),   LONG_CONC(NCP)
      DIMENSION C1S(KMC,IMC,NCP),  CSSB(KMC,IMC,NCP)
      DIMENSION CVP(KMC,NCP)
      DIMENSION SF2L(IMP),         SF3L(IMP,2),      SF4L(IMP,2),
     .          SF5L(IMP,2),       SF6L(IMP,2),      SF7L(IMP,2),
     .          SF8L(IMP,2),       SF9L(IMP,2),      SF10L(IMP,2),
     .          SF11L(IMP,2),      SF12L(IMP,2),     SF13L(IMP,2),
     .          SF14L(IMP,2)
      DIMENSION SF2V(KMP),         SF3V(KMP,2),      SF4V(KMP,2),
     .          SF5V(KMP,2),       SF6V(KMP,2),      SF7V(KMP,2),
     .          SF8V(KMP,2),       SF9V(KMP,2),      SF10V(KMP,2),
     .          SF11V(KMP,2)
      DIMENSION DX1(KMP,IMP),      DX2(KMP,IMP),     DX3(KMP,IMP),
     .          AD1L(KMP,IMP),     AD2L(KMP,IMP),    AD3L(KMP,IMP),
     .          AD1V(KMP,IMP),     AD2V(KMP,IMP),    AD3V(KMP,IMP),
     .          TADL(KMP,IMP),     TADV(KMP,IMP)
      DIMENSION CT(KMP,IMP),       AT(KMP,IMP),      BTAT(KMP,IMP),
     .          VT(KMP),           DT(KMP),          GMAT(KMP)
      DIMENSION NSTR(NBP),         QSTR(NSP,NBP),    ESTR(NSP,NBP),
     .          WSTR(NSP,NBP),     KBSW(NBP),        SINK(NSP,NBP)
      DIMENSION CADL(KMC,IMC,NCP),    CADV(KMC,IMC,NCP)
      DIMENSION CSSUH1(KMP,NCP,NBP),  CSSUH2(KMP,NCP,NBP),
     .          CSSDH1(KMP,NCP,NBP),  CSSDH2(KMP,NCP,NBP)
      DIMENSION SNPD(NDP), VPLD(NDP), PRFD(NDP), CPLD(NDP), RSOD(NDP),
     .          TSRD(NDP), DLTD(NDP), WSCD(NDP)
      DIMENSION SNPF(NDP), VPLF(NDP), PRFF(NDP), CPLF(NDP), RSOF(NDP),
     .          TSRF(NDP), DLTF(NDP)
      DIMENSION WSC(NDP),  DLTMAX(NDP)
      DIMENSION IPR(17),   IPRLF(17), IPRSF(11), TITLE(5),  HNAME(3),
     .          HPRN(3)
      DIMENSION SEL_WITHDRAWAL(NBP), POINT_SINK(NSP,NBP),   ICE_IN(NBP),
     .          ALLOW_ICE(IMP),      TRANSPORT(NCP)
      dimension gw1(nbp),gw2(nbp),cctr(22),ccdo(imp),ccbod(imp)
      double precision T1(KMP,IMP),   T2(KMP,IMP), ttt, tice

***** Common declarations

      COMMON /GLOBLC/ JB,    JC,     IU,     ID,    KT,    DLT,  KB
      COMMON /GEOMHC/ EL,    H,      HKT1,   HKT2
      COMMON /GEOMBC/ B,     BKT,    BH,     BHKT1, BHKT2, BHRKT1
      COMMON /NAMESC/ HNAME, CNAME,  CUNIT
      COMMON /SNPUC/  SNP
      COMMON /TEMPC/  T1,    T2
      COMMON /ICEC/   ICE,   ICETH,  ICE_CALC
      COMMON /HYDRC1/ U,     W,      AZ,     RHO
      COMMON /HYDRC2/ Z1,    Z2
      COMMON /PRNTC1/ IPR,   IBPR,   IEPR,   KEPR
      COMMON /PRNTC2/ TITLE, CPRN,   HPRN,   CONV,  CONVD1
      COMMON /PRNTC3/ CUS,   DS
      COMMON /GRDLGC/ LONG_FORM,      SHORT_FORM,     LIMITING_FACTOR
      COMMON /DNSPHC/ FRESH_WATER,    SALT_WATER,     SUSP_SOLIDS
      COMMON /GRTVDC/ CONSTITUENTS,   CN,             NAC
      COMMON /INTERC/ INTERP_INFLOW,  INTERP_OUTFLOW, INTERP_MET,
     .                INTERP_DTRIBS,  INTERP_TRIBS,   INTERP_HEAD,
     .                INTERP_WITHDRWL
      COMMON /TVDDSC/ NO_WIND,        NO_INFLOW,      NO_OUTFLOW,
     .                NO_HEAT
      COMMON /TVDLC1/ PRECIPITATION,  WITHDRAWALS,    TRIBUTARIES,
     .                DIST_TRIBS
      COMMON /TVDLC2/ UP_FLOW,        DN_FLOW,        UH_INTERNAL,
     .                UH_EXTERNAL,    DH_INTERNAL,    DH_EXTERNAL
      COMMON /TVDLC3/ OPEN_FILES,     TERM_BY_TERM
      COMMON /TVDMTC/ TAIR,   TDEW,   CLOUD,  PHI,    ET,     CSHE,
     .                SRO,    LAT,    LONG,   banksh
      COMMON /TVDFNC/ METFN,  QWDFN,  QOTFN,  QINFN,  TINFN,  CINFN,
     .                QTRFN,  TTRFN,  CTRFN,  QDTFN,  TDTFN,  CDTFN,
     .                PREFN,  TPRFN,  CPRFN,  EUHFN,  TUHFN,  CUHFN,
     .                EDHFN,  TDHFN,  CDHFN
      COMMON /TVDQC/  QIN,    QTR,    QDTR,   PR,     ELUH,   ELDH,
     .                QOUT,   QWD
      COMMON /TVDTC/  TIN,    TTR,    TDTR,   TPR,    TUH,    TDH
      COMMON /TVDCC1/ CIN,    CTR,    CDTR,   CPR,    CUH,    CDH
      COMMON /TVDCC2/ INCN,   TRCN,   DTCN,   PRCN,   UHCN,   DHCN
      COMMON /TVDCC3/ NACIN,  NACTR,  NACPR,  NACDT
      COMMON /RTMLTC/ OMT1,   OMT2,   NH3T1,  NH3T2,  NO3T1,  NO3T2,
     .                AGT1,   AGT2,   AGT3,   AGT4,   OMK1,   OMK2,
     .                NH3K1,  NH3K2,  NO3K1,  NO3K2,  AGK1,   AGK2,
     .                AGK3,   AGK4,   zoot1,  zoot2,  zoot3,  zoot4,
     .                zook1,  zook2,  zook3,  zook4
      COMMON /DKORGC/ DETD,   ORGD
      COMMON /DKSEDC/ SEDD,   SO2D,   SOD
      COMMON /DKNITC/ NH3D,   NO3D
      COMMON /DKBODC/ CBODD
      COMMON /GBLRTC/ OMRM,   NH3RM,  NO3RM,  AGRMR,AGRMF,zoormr,zoormf
      COMMON /GLBLCC/ PALT,   ALGDET, O2LIM,  WIND, wscdp, wsc
      COMMON /DKMLTC/ A1,     A2, A3

      COMMON /CLFRMC/ COLQ10, COLDK
      COMMON /ORGDKC/ SEDDK,  DETDK,  LABDK,  REFDK,  LRFDK
      COMMON /SETLC1/ SETIN,  SETOUT
      COMMON /SETLC2/ SSETL,  DSETL,  ASETL,  FESETL
      COMMON /PHYTGC/ AGR,    ARR,    AMR,    AER
      COMMON /PHYTC1/ AEXCR,  AMORT,  AGROW,  ARESP,  ASATUR, AHSN,
     .                AHSP
      COMMON /PHYTC2/ BETA,   EXH2O,  EXINOR, EXORG
      COMMON /PHOSPC/ PO4REL, BIOP,   PARTP, obiop
      COMMON /NITROC/ BION,   PARTN,  NH3DK,  NH3REL, NO3DK, obion
      COMMON /OXYGNC/ O2ORG,  O2ALG,  O2NH3,  O2RESP, o2xfact
      COMMON /CARBNC/ CO2REL, BIOC, obioc
      COMMON /IRONC/  FEREL
      COMMON /CBODC/  KBOD,   TBOD,   RBOD
      COMMON /LFACC/  LFAC,   LFPR
      COMMON /SELWC/  NSTR,   QSTR,   ESTR,   WSTR,   KBSW,   KTOPSW,
     .                KBOTSW, NOUT,   KOUT
      COMMON /GDAYC1/ GDAY,   DAYM,   JDAYG,  LEAP_YEAR
      COMMON /GDAYC2/ MONTH
      COMMON /TVDSWC/ SEL_WITHDRAWAL, POINT_SINK
      common /zoopc/ zmax,zmort,zeffic,pref1,pref2,zresp,zoomin,zs2p,
     .               zmu,talgae,zrt,zmt
      common/prison/aqbox,aqppipe
      common/pumpreg/qpstor
      common /massbal/ c1s,dlx,itr,iwd
     .                   ,up_head,dn_head,jbwd,kwd
     .                   ,qpr,qinf,quh1,
     .                    qdh1,end_run,jbtr,
     .                    qtrf,qdt
        COMMON /SCRNC2/ NWD,    NTR


***** Data statements

      DATA RHOA /1.25/,   RHOI   /1000.0/,   RHOICE /916.0/,
     .     RK1  /2.12/,   RL1    /333507.0/, RIMT   /0.0/,
     .     G    /9.8/,    VTOL   /100.0/,    BLANK  /'       '/,
     .     CP   /4186.0/, FRAZDZ /0.14/
c      DATA CON /10/, BTH /11/, RSI /12/, VPR /13/, LPR /14/
      DATA CON /1/, BTH /2/, RSI /3/, VPR /4/, LPR /5/
C      DATA RSO /21/, TSR /22/, PRF /23/, VPL /24/, CPL /25/, WRN /26/,
C     .     ERR /27/
      DATA RSO /6/, TSR /7/, PRF /8/, VPL /9/, CPL /10/, WRN /11/,
     .     ERR /12/
      NB = NBP

************************************************************************
**                     Task 3: Input Section                          **
************************************************************************

***** Open control file

      OPEN (CON,FILE=CONFN,STATUS='OLD')

***** Title cards

      READ (CON,*)
      READ (CON,1000) (TITLE(J),J=1,4)

***** Time control cards

      READ (CON,1010)  TMSTRT,TMEND,YEAR
      READ (CON,1020)  NDT,DLTMIN
      READ (CON,1030) (DLTD(J),J=1,NDT)
      READ (CON,1030) (DLTMAX(J),J=1,NDT)
      READ (CON,1030) (DLTF(J),J=1,NDT)

***** Grid definition cards

c      READ (CON,1020)  KT,DATUM
      READ (CON,1020)  KT,DATUM,addkt,subkt
      READ (CON,1030)  elpmpon,elpmpof,qpump,elpon4,elpof4,qpump4
      READ (CON,1040) (US(JB),DS(JB),UHS(JB),DHS(JB),PHI0(JB),JB=1,NBP)
      READ (CON,1030)  LAT,LONG

***** Initial condition cards

      READ (CON,1050)  IT2,IICETH,WTYPE
      READ (CON,1060)  VBC,MBC,PQC,PQTC,EVC,PRC
      READ (CON,1060)  INFIC,TRIC,DTRIC,HDIC,OUTIC,WDIC,METIC
      READ (CON,1060)  WINDC,QINC,QOUTC,HEATC
      READ (CON,1070)  ICEC,SLICE,SLHEAT,ALBEDO,HWI,BETAI,GAMMAI,ICEMIN,
     .                 ICET2
      READ (CON,1080)  SLTR,THETA,stest,stest2,stest3
      READ (CON,1090)  NWSC
      READ (CON,1030) (WSCD(J),J=1,NWSC)
      READ (CON,1030) (WSC(J),J=1,NWSC)
c      READ (CON,1030)  AX,IDX,AZMIN,DZMIN,DZMAX,CHEZY
      READ (CON,1030)  AX,IDX,AZMIN,DZMIN,DZMAX,CHEZY,cbhe,tsed

***** Inflow-outflow cards

      READ (CON,1060) (SWC(JB),JB=1,NBP)
      READ (CON,1090) (NSTR(JB),JB=1,NBP)
      READ (CON,1090) (KBSW(JB),JB=1,NBP)
      READ (CON,1100) (SINK(JS,1),JS=1,NSTR(1))
      DO JB=2,NB
        READ (CON,1110) (SINK(JS,JB),JS=1,NSTR(JB))
      END DO
      READ (CON,1030) (ESTR(JS,1),JS=1,NSTR(1))
      DO JB=2,NB
        READ (CON,1120) (ESTR(JS,JB),JS=1,NSTR(JB))
      END DO
      READ (CON,1030) (WSTR(JS,1),JS=1,NSTR(1))
      DO JB=2,NB
        READ (CON,1120) (WSTR(JS,JB),JS=1,NSTR(JB))
      END DO
      READ (CON,1090) (NOUT(JB),JB=1,NBP)
      READ (CON,1090) (KOUT(JO,1),JO=1,NOUT(1))
      DO JB=2,NB
        READ (CON,1130) (KOUT(JO,JB),JO=1,NOUT(JB))
      END DO
      READ (CON,1090)  NWD
      READ (CON,1090) (IWD(JW),JW=1,NWD)
      READ (CON,1090) (KWD(JW),JW=1,NWD)
c read in irrigation amounts from control file
      read(con,1030)qwd1,qwd2,qwd3,qwd4,qwd5,qwd6,qwd7,qwd8,qwd9,
     . qwd10,qwd11,qwd13
      READ (CON,1090)  NTR
      READ (CON,1090) (ITR(JT),JT=1,NTP)
      READ (CON,1060) (DTRC(JB),JB=1,NBP)
c read in fraction of qflow that goes to each distributed tributary
c      read(con,1030)qdtr11,qdtr12,qdtr15,qdtr16,qdtr25,
c     . qdtr27,qdtr30,qdtr31
c reading groundwater coefficient 1 and 2 - chris
      read(con,1030)(gw1(jb),jb=1,nbp)
       read(con,1030)(gw2(jb),jb=1,nbp)

***** Output control cards (excluding constituents)

      READ (CON,1140)  FORM,(HPRN(J),J=1,3)
      READ (CON,1090) (IPRSF(I),I=1,11)
      READ (CON,1090) (IPRLF(I),I=1,17)
      READ (CON,1150)  SNPC,NSNP
      READ (CON,1030) (SNPD(J),J=1,NSNP)
      READ (CON,1030) (SNPF(J),J=1,NSNP)
      READ (CON,1160)  PRFC,NPRF,NIPRF
      READ (CON,1030) (PRFD(J),J=1,NPRF)
      READ (CON,1030) (PRFF(J),J=1,NPRF)
      READ (CON,1090) (IPRF(J),J=1,NIPRF)
      READ (CON,1150)  TSRC,NTSR
      READ (CON,1030) (TSRD(J),J=1,NTSR)
      READ (CON,1030) (TSRF(J),J=1,NTSR)
      READ (CON,1150)  VPLC,NVPL
      READ (CON,1030) (VPLD(J),J=1,NVPL)
      READ (CON,1030) (VPLF(J),J=1,NVPL)
      READ (CON,1150)  CPLC,NCPL
      READ (CON,1030) (CPLD(J),J=1,NCPL)
      READ (CON,1030) (CPLF(J),J=1,NCPL)
      READ (CON,1150)  RSOC,NRSO,RSIC
      READ (CON,1030) (RSOD(J),J=1,NRSO)
      READ (CON,1030) (RSOF(J),J=1,NRSO)

***** Constituent control cards

      READ (CON,1170)  CCC,LIMC,SDC,FREQUK
      READ (CON,1060) (ACC(JC),JC=1,NCP)
c      READ (CON,1030) (CIC(JC),JC=1,NCP)
      READ (CON,1030) (CICbr(JC,1),JC=1,NCP)
      DO JB=2,NB
        READ (CON,1030) (cicbr(jc,JB),Jc=1,ncp)
      END DO
      READ (CON,1060) (CPRN(JC),JC=1,NCP)
      READ (CON,1060) (INACC(JC),JC=1,NCP)
      READ (CON,1060) (TRACC(JC),JC=1,NCP)
      READ (CON,1060) (DTACC(JC),JC=1,NCP)
      READ (CON,1060) (PRACC(JC),JC=1,NCP)

***** Kinetics coefficients

      READ (CON,1030)  EXH2O,EXINOR,EXORG,BETA
      READ (CON,1030)  COLQ10,COLDK
      READ (CON,1030)  SSETL
      READ (CON,1030)  AGROW,AMORT,AEXCR,ARESP,ASETL,ASATUR,ALGDET
      READ (CON,1030)  AGT1,AGT2,AGT3,AGT4,AGK1,AGK2,AGK3,AGK4
      READ (CON,1030)  LABDK,LRFDK,REFDK
      READ (CON,1030)  DETDK,DSETL
      READ (CON,1030)  OMT1,OMT2,OMK1,OMK2
      READ (CON,1030)  SEDDK
      READ (CON,1030) (SOD(I),I=1,IMP)
      READ (CON,1030)  KBOD,TBOD,RBOD
      READ (CON,1030)  PO4REL,PARTP,AHSP
      READ (CON,1030)  NH3REL,NH3DK,PARTN,AHSN
      READ (CON,1030)  NH3T1,NH3T2,NH3K1,NH3K2
      READ (CON,1030)  NO3DK
      READ (CON,1030)  NO3T1,NO3T2,NO3K1,NO3K2
      READ (CON,1030)  CO2REL
      READ (CON,1030)  FEREL,FESETL
c zooplankton variables
      read(con,1030)zmax,zmort,zeffic,pref1,pref2,zresp,zoomin,
     .              zs2p
      read(con,1030)zoot1,zoot2,zoot3,zoot4,zook1,zook2,zook3,zook4
      READ (CON,1030)  O2NH3,O2ORG,O2RESP,O2ALG,BIOP,BION,BIOC
      READ (CON,1030)  O2LIM,o2xfact,obiop,obion,obioc

***** Input filenames

      READ (CON,1000)  BTHFN
      READ (CON,1000)  VPRFN
      READ (CON,1000)  LPRFN
      READ (CON,1000)  RSIFN
      READ (CON,1000)  METFN
      READ (CON,1000)  QWDFN
      READ (CON,1000) (QINFN(JB),JB=1,NBP)
      READ (CON,1000) (TINFN(JB),JB=1,NBP)
      READ (CON,1000) (CINFN(JB),JB=1,NBP)
      READ (CON,1000) (QOTFN(JB),JB=1,NBP)
      READ (CON,1000) (QTRFN(JT),JT=1,NTP)
      READ (CON,1000) (TTRFN(JT),JT=1,NTP)
      READ (CON,1000) (CTRFN(JT),JT=1,NTP)
      READ (CON,1000) (QDTFN(JB),JB=1,NBP)
      READ (CON,1000) (TDTFN(JB),JB=1,NBP)
      READ (CON,1000) (CDTFN(JB),JB=1,NBP)
      READ (CON,1000) (PREFN(JB),JB=1,NBP)
      READ (CON,1000) (TPRFN(JB),JB=1,NBP)
      READ (CON,1000) (CPRFN(JB),JB=1,NBP)
      READ (CON,1000) (EUHFN(JB),JB=1,NBP)
      READ (CON,1000) (TUHFN(JB),JB=1,NBP)
      READ (CON,1000) (CUHFN(JB),JB=1,NBP)
      READ (CON,1000) (EDHFN(JB),JB=1,NBP)
      READ (CON,1000) (TDHFN(JB),JB=1,NBP)
      READ (CON,1000) (CDHFN(JB),JB=1,NBP)

***** Output filenames

      READ (CON,1000)  SNPFN
      READ (CON,1000)  TSRFN
      READ (CON,1000)  PRFFN
      READ (CON,1000)  VPLFN
      READ (CON,1000)  CPLFN

***** Bathymetry file

      OPEN (BTH,FILE=BTHFN,STATUS='OLD')
      READ (BTH,*)
      READ (BTH,1180) (DLX(I),I=1,IMP)
      READ (BTH,1180) (Z2(I),I=1,IMP)
      READ (BTH,1180) (H(K),K=1,KMP)
c read in manning's n and phi0 for each cell
      read(bth,1030)(sman(i),i=1,imp)
c bank shading factor - reduce shortwave inputs to each cell by shading
c factor is between 0 (full shading) and 1 (no shading)
      read(bth,1030)(banksh(i),i=1,imp)
c      read(bth,1030)(phi0(i),i=1,imp) ...first fix problem with variable fetch
c
      DO I=1,IMP
        READ (BTH,1180) (B(K,I),K=1,KMP)
      END DO

***** Initialize logical control variables

      CONSTITUENTS = CCC.EQ.' ON'
      RESTART_IN   = RSIC.EQ.' ON'
      ISO_TEMP     = INT(IT2).GE.0
      VERT_TEMP    = INT(IT2).EQ.-1
      LONG_TEMP    = INT(IT2).LT.-1
      OPEN_VPR     = VERT_TEMP
      OPEN_LPR     = LONG_TEMP
c     don't need cause assigning initial constituent concentrations
c     by branch
c      DO JC=1,NCP
c        ISO_CONC(JC)  = INT(CIC(JC)).GE.0
c        VERT_CONC(JC) = INT(CIC(JC)).EQ.-1
c        LONG_CONC(JC) = INT(CIC(JC)).LT.-1
c        IF (VERT_CONC(JC)) OPEN_VPR = .TRUE.
c        IF (LONG_CONC(JC)) OPEN_LPR = .TRUE.
c      END DO

***** Temperature and constituents

      IF (OPEN_VPR) THEN
        OPEN (VPR,FILE=VPRFN,STATUS='OLD')
        READ (VPR,*)
        IF (VERT_TEMP) READ (VPR,1030) (TVP(K),K=KT,KMP)
        IF (CONSTITUENTS) THEN
          DO JC=1,NCP
            IF (VERT_CONC(JC)) READ (VPR,1030) (CVP(K,JC),K=KT,KMP)
          END DO
        ENDIF
      ENDIF
      IF (OPEN_LPR) THEN
        OPEN (LPR,FILE=LPRFN,STATUS='OLD')
        READ (LPR,*)
      END IF

***** Restart data

      IF (RESTART_IN) THEN
        OPEN (RSI,FILE=RSIFN,STATUS='OLD',FORM='UNFORMATTED')
        READ (RSI) KT,     NIT,    NV,     KMIN,   IMIN
        READ (RSI) DLTDP,  SNPDP,  TSRDP,  VPLDP,  PRFDP,  CPLDP,
     .             RSODP,  WSCDP
        READ (RSI) JDAY,   JDAYG,  YEAR,   ELTM,   DLT,    AVDLT,
     .             MINDLT, JDMIN,  CURMAX
        READ (RSI) NXTMSN, NXTMTS, NXTMPR, NXTMCP, NXTMVP, NXTMRS
        READ (RSI) VOLIN,  VOLOUT, VOLUH,  VOLDH,  VOLPR,  VOLTR,
     .             VOLDT,  VOLWD,  VOLEV,  VOLSBR, ICMBR,  CSSMBR
        READ (RSI) TSSUH2, TSSDH2, CSSUH2, CSSDH2, QUH2,   QDH2
        READ (RSI) Z1,     Z2,     SZ1,    SKTT,   ICE,    ICETH
        READ (RSI) U,      W,      SU,     SW,     AZ,     T1,
     .             T2,     C1,     C2
      ENDIF

***** Close files

      CLOSE (CON)
      CLOSE (BTH)
      IF (OPEN_VPR)   CLOSE (VPR)
      IF (RESTART_IN) CLOSE (RSI)

***** Input FORMATs

 1000 FORMAT(//(8X,A72))
 1010 FORMAT(//(8X,2F8.0,I8))
 1020 FORMAT(//8X,I8,8F8.0)
 1030 FORMAT(//(8X,9F8.0))
 1040 FORMAT(//(8X,4I8,F8.0))
 1050 FORMAT(//8X,2F8.0,3X,A5,6(5X,A3))
 1060 FORMAT(//(8X,9(5X,A3)))
 1070 FORMAT(//13X,A3,2A8,6F8.0)
 1080 FORMAT(//8X,A8,F8.0,3A8)
 1090 FORMAT(//(8X,9I8))
 1100 FORMAT(//(8X,9A8))
 1110 FORMAT(:8X,9A8)
 1120 FORMAT(:8X,9F8.0)
 1130 FORMAT(:8X,9I8)
 1140 FORMAT(//11X,A5,8(5X,A3))
 1150 FORMAT(//13X,A3,I8,5X,A3)
 1160 FORMAT(//13X,A3,8I8)
 1170 FORMAT(//8X,3(5X,A3),I8)
 1180 FORMAT(//(10F8.0))

************************************************************************
**                  Task 4: Initialization Section                    **
************************************************************************

************************************************************************
**                     Task 4.1: Zero Variables                       **
************************************************************************

      DO I=1,17
        DO K=1,KMP
          CONV(K,I) = BLANK
        END DO
      END DO
      DO I=1,11
        DO K=1,KMP
          CONVD1(K,I) = '             '
        END DO
      END DO
      DO I=1,IMP
        DO K=1,KMP
          LFPR(K,I) = '       '
        END DO
      END DO
      DO I=1,IMP
        INIT(I)  = 0
        ICESW(I) = 1.0
        DO JB=1,NBP
          FETCHU(I,JB) = 0.0
          FETCHD(I,JB) = 0.0
        END DO
        IF (.NOT.RESTART_IN) THEN
          Z1(I) = 0.0
          DO K=1,KMP
            U(K,I) = 0.0
            W(K,I) = 0.0
          END DO
        END IF
        DO K=1,KMP
          TSS(K,I) = 0.0
          QSS(K,I) = 0.0
          IF (CONSTITUENTS) THEN
            DO JC=1,NCP
              IF (.NOT.RESTART_IN) THEN
                CSSB(K,I,JC) = 0.0
                CSSK(K,I,JC) = 0.0
              END IF
            END DO
          END IF
        END DO
      END DO
      DO JB=1,NBP
        IF (.NOT.RESTART_IN) THEN
          VOLEV(JB)  = 0.0
          VOLPR(JB)  = 0.0
          VOLTR(JB)  = 0.0
          VOLDT(JB)  = 0.0
          VOLWD(JB)  = 0.0
          VOLUH(JB)  = 0.0
          VOLDH(JB)  = 0.0
          VOLIN(JB)  = 0.0
          VOLOUT(JB) = 0.0
          VOLSBR(JB) = 0.0
        END IF
        DO K=1,KMP
          AKBR(K,JB) = 0.0
        END DO
      END DO
      DLMR  = 0.0
      DLVR  = 0.0
      KEPR  = 0
      NACIN = 0
      NACTR = 0
      NACDT = 0
      NACPR = 0

************************************************************************
**            Task 4.2: Initialize Miscellaneous Variables            **
************************************************************************

***** Logical controls

      OPEN_FILES      = .TRUE.
      VOLUME_WARNING  = .TRUE.
      UPDATE_KINETICS = .TRUE.
      END_RUN         = .FALSE.
      VIOLATION       = .FALSE.
      WARNING_OPEN    = .FALSE.
      TRIBUTARIES     = NTR.GT.0
      WITHDRAWALS     = NWD.GT.0
      PLACE_QIN       = PQC.EQ.' ON'
      EVAPORATION     = EVC.EQ.' ON'
      VOL_BALANCE     = VBC.EQ.' ON'
      MASS_BALANCE    = MBC.EQ.' ON'
      SHIFT_DEMAND    = SDC.EQ.' ON'
      PRECIPITATION   = PRC.EQ.' ON'
      SNAPSHOT        = SNPC.EQ.' ON'
      CONTOUR         = CPLC.EQ.' ON'
      VECTOR          = VPLC.EQ.' ON'
      PROFILE         = PRFC.EQ.' ON'
      ICE_CALC        = ICEC.EQ.' ON'
      PLACE_QTR       = PQTC.EQ.' ON'
      TIME_SERIES     = TSRC.EQ.' ON'
      RESTART_OUT     = RSOC.EQ.' ON'
      INTERP_TRIBS    = TRIC.EQ.' ON'
      INTERP_HEAD     = HDIC.EQ.' ON'
      INTERP_WITHDRWL = WDIC.EQ.' ON'
      INTERP_MET      = METIC.EQ.' ON'
      INTERP_INFLOW   = INFIC.EQ.' ON'
      INTERP_DTRIBS   = DTRIC.EQ.' ON'
      INTERP_OUTFLOW  = OUTIC.EQ.' ON'
      NO_INFLOW       = QINC.EQ.'OFF'
      NO_HEAT         = HEATC.EQ.'OFF'
      NO_WIND         = WINDC.EQ.'OFF'
      NO_OUTFLOW      = QOUTC.EQ.'OFF'
      LONG_FORM       = FORM.EQ.' LONG'
      SHORT_FORM      = FORM.NE.' LONG'
      UPWIND          = SLTR.EQ.'  UPWIND'
      test            = stest.eq.'      ON'
      test2            = stest2.eq.'      ON'
      test3            = stest3.eq.'      ON'
      TERM_BY_TERM    = SLHEAT.NE.'      ET'
      DETAILED_ICE    = SLICE.EQ.'  DETAIL'
     .                  .OR.(ICE_CALC.AND.DETAILED_ICE)
      SUSP_SOLIDS     = CONSTITUENTS.AND.ACC(2).EQ.' ON'
      LIMITING_FACTOR = CONSTITUENTS.AND.ACC(7).EQ.' ON'
     .                  .AND.LIMC.EQ.' ON'
      FRESH_WATER     = CONSTITUENTS.AND.ACC(4).EQ.' ON'
     .                  .AND.WTYPE.EQ.'FRESH'
      SALT_WATER      = CONSTITUENTS.AND.ACC(4).EQ.' ON'
     .                  .AND.WTYPE.EQ.' SALT'
      INTERPOLATE     = INTERP_INFLOW.OR.INTERP_TRIBS.OR.INTERP_DTRIBS
     .                  .OR.INTERP_OUTFLOW.OR.INTERP_WITHDRWL
     .                  .OR.INTERP_HEAD.OR.INTERP_MET
      LEAP_YEAR       = MOD(YEAR,4).EQ.0
      IF (RESTART_IN) THEN
        WINTER = .FALSE.
        IF (JDAY.GT.300.0.OR.JDAY.LT.40.0) WINTER = .TRUE.
      ELSE
        WINTER = .FALSE.
        IF (TMSTRT.GT.300.0.OR.TMSTRT.LT.40.0) WINTER = .TRUE.
      END IF
      DO JC=1,NCP
        TRANSPORT(JC) = .TRUE.
        IF (JC.EQ.13.OR.(JC.GT.15.AND.JC.LT.20)) TRANSPORT(JC) = .FALSE.
      END DO
      DO JB=1,NBP
        UP_FLOW(JB)        = UHS(JB).EQ.0
        DN_FLOW(JB)        = DHS(JB).EQ.0
        UP_HEAD(JB)        = UHS(JB).NE.0
        DN_HEAD(JB)        = DHS(JB).NE.0
        UH_INTERNAL(JB)    = UHS(JB).GT.0
        DH_INTERNAL(JB)    = DHS(JB).GT.0
        UH_EXTERNAL(JB)    = UHS(JB).EQ.-1
        DH_EXTERNAL(JB)    = DHS(JB).EQ.-1
        DIST_TRIBS(JB)     = DTRC(JB).EQ.' ON'
        SEL_WITHDRAWAL(JB) = SWC(JB).EQ.' ON'
        DO JS=1,NSTR(JB)
          POINT_SINK(JS,JB) = SINK(JS,JB).EQ.'   POINT'
        END DO
        IF (UH_EXTERNAL(JB).OR.DH_EXTERNAL(JB)) HEAD_BOUNDARY = .TRUE.
      END DO

***** Convert rates from per-day to per-second

      KBOD   = KBOD/86400.0
      ASETL  = ASETL/86400.0
      SSETL  = SSETL/86400.0
      DSETL  = DSETL/86400.0
      LABDK  = LABDK/86400.0
      REFDK  = REFDK/86400.0
      LRFDK  = LRFDK/86400.0
      NH3DK  = NH3DK/86400.0
      NO3DK  = NO3DK/86400.0
      DETDK  = DETDK/86400.0
      COLDK  = COLDK/86400.0
      SEDDK  = SEDDK/86400.0
      AEXCR  = AEXCR/86400.0
      AMORT  = AMORT/86400.0
      ARESP  = ARESP/86400.0
      AGROW  = AGROW/86400.0
      FESETL = FESETL/86400.0
      zmax=zmax/86400.
      zresp=zresp/86400.
      zmort=zmort/86400.
      DO I=1,IMP
        SOD(I)  = SOD(I)/86400.0
        SODS(I) = SOD(I)
      END DO

***** Time and printout control variables

      IF (.NOT.RESTART_IN) THEN
        JDAY   = TMSTRT
        JDAYG  = JDAY
        ELTM   = TMSTRT*86400.0
        DLT    = DLTMAX(1)
        DLTS   = DLT
        MINDLT = DLT
        NIT    = 0
        NV     = 0
        DLTDP  = 1
        WSCDP  = 1
        SNPDP  = 1
        TSRDP  = 1
        VPLDP  = 1
        PRFDP  = 1
        CPLDP  = 1
        RSODP  = 1
        NXTMSN = SNPD(1)
        NXTMTS = TSRD(1)
        NXTMPR = PRFD(1)
        NXTMCP = CPLD(1)
        NXTMVP = VPLD(1)
        NXTMRS = RSOD(1)
c    time counter for environmental performance criterion
        NXTEP=Tmstrt+0.02083
      END IF
      DO J=NWSC+1,NDP
        WSCD(J) = TMEND+1.0
      END DO
      DO J=NDT+1,NDP
        DLTD(J) = TMEND+1.0
      END DO
      DO J=NSNP+1,NDP
        SNPD(J) = TMEND+1.0
      END DO
      DO  J=NTSR+1,NDP
        TSRD(J) = TMEND+1.0
      END DO
      DO J=NPRF+1,NDP
        PRFD(J) = TMEND+1.0
      END DO
      DO J=NVPL+1,NDP
        VPLD(J) = TMEND+1.0
      END DO
      DO J=NCPL+1,NDP
        CPLD(J) = TMEND+1.0
      END DO
      DO J=NRSO+1,NDP
        RSOD(J) = TMEND+1.0
      END DO
      CURMAX = DLTMAX(DLTDP)/DLTF(DLTDP)
      NXTVD  = JDAY
      JDAYNX = JDAYG+1
      NSPRF  = 0
      boddonx=45.0

***** Active constituents

      IF (CONSTITUENTS) THEN
        DO JC=1,NCP
          IF (ACC(JC).EQ.' ON') THEN
            NAC     = NAC+1
            CN(NAC) = JC
          END IF
          IF (INACC(JC).EQ.' ON') THEN
            NACIN       = NACIN+1
            INCN(NACIN) = JC
          END IF
          IF (TRACC(JC).EQ.' ON') THEN
            NACTR       = NACTR+1
            TRCN(NACTR) = JC
          END IF
          IF (DTACC(JC).EQ.' ON') THEN
            NACDT       = NACDT+1
            DTCN(NACDT) = JC
          END IF
          IF (PRACC(JC).EQ.' ON') THEN
            NACPR       = NACPR+1
            PRCN(NACPR) = JC
          END IF
        END DO
      END IF
      TITLE(5) = ' '

************************************************************************
**                 Task 4.3: Initialize Geometry                      **
************************************************************************

***** Geometry branch by branch

      DO JB=1,NBP
        IU = US(JB)
        ID = DS(JB)

******* Bottom layers

        DO I=IU,ID
          K = 2
          DO WHILE (B(K,I).GT.0.0)
            KB(I) = K
            K     = K+1
          END DO
        END DO
        KB(IU-1) = KB(IU)
        KB(ID+1) = KB(ID)
      END DO
      DO JB=1,NBP
        IU = US(JB)
        ID = DS(JB)

******* Upstream active cell

        I = IU
        DO WHILE (KB(I)-KT.LT.1.AND.I.LT.ID)
          I = I+1
        END DO
        IUT     = I
        CUS(JB) = IUT

******* Boundary bottom layers

        IF (UH_EXTERNAL(JB)) KB(IUT-1) = KB(IUT)
        IF (DH_EXTERNAL(JB)) KB(ID+1)  = KB(ID)
        IF (UH_INTERNAL(JB)) KB(IUT-1) = MIN(KB(UHS(JB)),KB(IUT))
        IF (DH_INTERNAL(JB)) KB(ID+1)  = MIN(KB(DHS(JB)),KB(ID))

******* Minimum bottom layers

        DO I=IU-1,ID
          KBMIN(I) = MIN(KB(I),KB(I+1))
        END DO

******* Boundary widths

        DO I=IU,ID
          B(1,I) = B(2,I)
          DO K=KB(I)+1,KMP
            B(K,I) = B(KB(I),I)
          END DO
        END DO
        DO K=1,KB(IU)
          B(K,IU-1) = B(K,IU)
          IF (UH_INTERNAL(JB)) B(K,IU-1) = B(K,UHS(JB))
        END DO
        DO K=1,KB(ID)
          B(K,ID+1) = B(K,ID)
          IF (DH_INTERNAL(JB)) B(K,ID+1) = B(K,DHS(JB))
        END DO

******* Initial water surface and derived geometry

        DO I=IU-1,ID+1
          IF (.NOT.RESTART_IN) THEN
            Z1(I)  = Z2(I)
            SZ1(I) = Z2(I)
          END IF
          ZMIN      = MAX(-1000.0,SNGL(Z2(I)))
          HKT2(I)   = H(KT)-Z2(I)
          AVHKT(I)  = (HKT2(I)+H(KT+1))*0.5
          BHKT2(I)  = 0.0
          HEIGHT(I) = 0.0
          K         = KT-1
          KTT(I)    = KT
          IF (.NOT.RESTART_IN) SKTT(I)  = KTT(I)
          DO WHILE (HEIGHT(I).LT.-Z2(I))
            HEIGHT(I) = HEIGHT(I)+H(K)
            KTT(I)    = MAX(K,2)
            SKTT(I)   = KTT(I)
            K         = K-1
          END DO
          HEIGHT(I) = HEIGHT(I)-H(K+1)
          DO K=KT,KTT(I)+1,-1
            BHKT2(I) = BHKT2(I)+B(K,I)*H(K)
          END DO
          IF (Z2(I).GE.0.0) THEN
            BHKT2(I) = B(KT,I)*HKT2(I)
          ELSE
            BHKT2(I) = BHKT2(I)-(B(KTT(I),I)*(HEIGHT(I)+Z2(I)))
          END IF
          BKT(I) = BHKT2(I)/HKT2(I)
          DO K=1,KMP-1
            BH(K,I) = B(K,I)*H(K)
            BB(K,I) = (B(K,I)+B(K+1,I))*0.5
          END DO
          BH(KB(I)+1,I) = BH(KB(I),I)
        END DO
        IDT = ID+1
        IF (JB.EQ.NBP) IDT = ID
        DO I=IU-1,IDT
          BHRKT2(I) = (BHKT2(I)+BHKT2(I+1))*0.5
          DO K=1,KMP-1
            BR(K,I)  = (B(K,I)+B(K,I+1))*0.5
            BHR(K,I) = (BH(K,I)+BH(K,I+1))*0.5
          END DO
        END DO
        DO K=1,KMP-1
          AVH(K) = (H(K)+H(K+1))*0.5
        END DO

******* Branch numbers corresponding to tributaries, withdrawals, and head

        IF (TRIBUTARIES) THEN
          DO JT=1,NTP
            IF (ITR(JT).GE.US(JB).AND.ITR(JT).LE.DS(JB)) THEN
              JBTR(JT) = JB
            END IF
          END DO
        ENDIF
        IF (WITHDRAWALS) THEN
          DO JW=1,NWP
            IF (IWD(JW).GE.US(JB).AND.IWD(JW).LE.DS(JB)) THEN
              JBWD(JW) = JB
            END IF
          END DO
        END IF
        IF (UH_INTERNAL(JB)) THEN
          JBUH(JB)     = 0
          BRANCH_FOUND = .FALSE.
          DO WHILE (.NOT.BRANCH_FOUND)
            JBUH(JB) = JBUH(JB)+1
            DO I=CUS(JBUH(JB)),DS(JBUH(JB))
              IF (I.EQ.UHS(JB)) BRANCH_FOUND = .TRUE.
            END DO
          END DO
        END IF
        IF (DH_INTERNAL(JB)) THEN
          JBDH(JB)     = 0
          BRANCH_FOUND = .FALSE.
          DO WHILE (.NOT.BRANCH_FOUND)
            JBDH(JB) = JBDH(JB)+1
            DO I=CUS(JBDH(JB)),DS(JBDH(JB))
              IF (I.EQ.DHS(JB)) BRANCH_FOUND = .TRUE.
            END DO
          END DO
        END IF

******* Branch layer area

        DO K=KMP-1,2,-1
          DO I=IUT,ID
            IF (K.LE.KB(I)) AKBR(K,JB) = AKBR(K,JB)+B(K,I)*DLX(I)
          END DO
        END DO

******* Wind fetch lengths

        DO I=IUT,ID
          FETCHD(I,JB) = FETCHD(I-1,JB)+DLX(I)
        END DO
        DO I=ID,IUT,-1
          FETCHU(I,JB) = FETCHU(I+1,JB)+DLX(I)
        END DO
      END DO

***** Boundary bottom counters and segment heights

      DO I=1,IMP
        B(KB(I)+1,I) = B(KB(I),I)
        DO K=KB(I),2,-1
          HSEG(K,I) = HSEG(K+1,I)+H(K)
        END DO
      END DO

***** Layer elevations

      EL(KMP) = DATUM
      DO K=KMP-1,1,-1
        EL(K) = EL(K+1)+H(K)
      END DO

***** Beginning and ending segment and bottom layer for snapshots

      IF (LONG_FORM) THEN
        IEPR = MIN(IMP,17)
        DO I=1,IEPR
          IPR(I) = IPRLF(I)
        END DO
      ELSE
        IEPR = MIN(IMP,11)
        DO I=1,IEPR
          IPR(I) = IPRSF(I)
        END DO
      END IF
      DO I=1,IEPR
        KEPR = MAX(KB(IPR(I)),KEPR)
      END DO
      IBPR = 1
      DO WHILE (CUS(1).GT.IPR(IBPR))
        IBPR = IBPR+1
      END DO

***** Interpolation multipliers

      DO I=2,IMP-1
        DO K=2,KB(I)

********* Positive flows

          DLXT = DLX(I-1)
          IF (K.GT.KB(I-1)) DLXT = DLX(I)
          DLXMIN     = MIN(DLX(I+1),DLX(I))
          SF2L(I)    = (DLX(I+1)+DLX(I))/2.0
          SF3L(I,1)  = DLX(I)/(DLX(I)+DLX(I+1))
          SF4L(I,1)  = DLX(I)**2
          SF5L(I,1)  = DLX(I+1)/(DLX(I)+DLX(I+1))
          SF6L(I,1)  = 0.25*(DLXT+2.0*DLX(I)+DLX(I+1))*(DLXT+DLX(I))
          SF7L(I,1)  = -0.25*(DLX(I)+DLX(I+1))*(DLXT+DLX(I))
          SF8L(I,1)  = 0.25*(DLX(I)+DLX(I+1))*(DLXT+2.0*DLX(I)+DLX(I+1))
          SF9L(I,1)  = 0.5*(DLX(I)-DLX(I+1))*DLXMIN
          SF10L(I,1) = 0.5*(DLXT+2.0*DLX(I)-DLX(I+1))*DLXMIN
          SF11L(I,1) = 0.5*(DLXT+3.0*DLX(I))*DLXMIN
          SF12L(I,1) = SF9L(I,1)/SF6L(I,1)/SF2L(I)
          SF13L(I,1) = SF10L(I,1)/SF7L(I,1)/SF2L(I)
          SF14L(I,1) = SF11L(I,1)/SF8L(I,1)/SF2L(I)
          HTOP       = H(K-1)
          HMID       = H(K)
          HBOT       = H(K+1)
          HMIN       = MIN(HBOT,HMID)
          SF2V(K)    = (HBOT+HMID)/2.0
          SF3V(K,1)  = HMID**2
          SF4V(K,1)  = HMID/(HMID+HBOT)
          SF5V(K,1)  = HBOT/(HMID+HBOT)
          SF6V(K,1)  = 0.25*(HTOP+2.0*HMID+HBOT)*(HTOP+HMID)
          SF7V(K,1)  = -0.25*(HMID+HBOT)*(HTOP+HMID)
          SF8V(K,1)  = 0.25*(HMID+HBOT)*(HTOP+2.0*HMID+HBOT)
          SF9V(K,1)  = 0.5*(HMID-HBOT)*HMIN
          SF10V(K,1) = 0.5*(HTOP+2.0*HMID-HBOT)*HMIN
          SF11V(K,1) = 0.5*(HTOP+3.0*HMID)*HMIN

********* Negative flows

          IF (I.LT.IMP-1) THEN
            DLXT = DLX(I+2)
            IF (K.GT.KB(I+2)) DLXT = DLX(I+1)
            DLXMIN     = MIN(DLX(I),DLX(I+1))
            SF2L(I)    = (DLX(I+1)+DLX(I))/2.0
            SF3L(I,2)  = DLX(I+1)/(DLX(I)+DLX(I+1))
            SF4L(I,2)  = DLX(I+1)**2
            SF5L(I,2)  = DLX(I)/(DLX(I)+DLX(I+1))
            SF6L(I,2)  = 0.25*(DLX(I)+2.0*DLX(I+1)+DLXT)*(DLX(I)
     .                   +DLX(I+1))
            SF7L(I,2)  = -0.25*(DLX(I+1)+DLXT)*(DLX(I)+DLX(I+1))
            SF8L(I,2)  = 0.25*(DLX(I)+2.0*DLX(I+1)+DLXT)*(DLX(I+1)+DLXT)
            SF9L(I,2)  = -0.5*(3.0*DLX(I+1)+DLXT)*DLXMIN
            SF10L(I,2) = 0.5*(DLX(I)-2.0*DLX(I+1)-DLXT)*DLXMIN
            SF11L(I,2) = 0.5*(DLX(I)-DLX(I+1))*DLXMIN
            SF12L(I,2) = SF9L(I,2)/SF6L(I,2)/SF2L(I)
            SF13L(I,2) = SF10L(I,2)/SF7L(I,2)/SF2L(I)
            SF14L(I,2) = SF11L(I,2)/SF8L(I,2)/SF2L(I)
          END IF
          HTOP = H(K)
          HMID = H(K+1)
          IF (K.LT.KB(I)) THEN
            HBOT = H(K+2)
            IF (K.EQ.KB(I)-1) HBOT = H(K+1)
            HMIN       = MIN(HTOP,HMID)
            SF2V(K)    = (HMID+HTOP)/2.0
            SF3V(K,2)  = HMID**2
            SF4V(K,2)  = HMID/(HTOP+HMID)
            SF5V(K,2)  = HTOP/(HTOP+HMID)
            SF6V(K,2)  = 0.25*(HTOP+2.0*HMID+HBOT)*(HTOP+HMID)
            SF7V(K,2)  = -0.25*(HMID+HBOT)*(HTOP+HMID)
            SF8V(K,2)  = 0.25*(HTOP+2.0*HMID+HBOT)*(HMID+HBOT)
            SF9V(K,2)  = -0.5*(3.0*HMID+HBOT)*HMIN
            SF10V(K,2) = 0.5*(HTOP-2.0*HMID-HBOT)*HMIN
            SF11V(K,2) = 0.5*(HTOP-HMID)*HMIN
          END IF
        END DO
      END DO

************************************************************************
**            Task 4.4: Initial Conditions for Simulation             **
************************************************************************

      DO JB=1,NBP
        IU = CUS(JB)
        ID = DS(JB)
        IF (.NOT.RESTART_IN) THEN

********* Temperature

          DO I=IU,ID
            IF (LONG_TEMP) READ (LPR,1030) (T1(K,I),K=KT,KB(I))
            DO K=KT,KB(I)
              IF (ISO_TEMP)  T1(K,I) = IT2
              IF (VERT_TEMP) T1(K,I) = TVP(K)
              T2(K,I) = T1(K,I)
            END DO
          END DO

********* Constituents

          DO JC=1,NAC
            DO I=IU,ID
c     assigning initial constituent concentrations
c     by branch
c              IF (LONG_CONC(CN(JC))) THEN
c                READ (LPR,1030) (C2(K,I,CN(JC)),K=KT,KB(I))
c              END IF
              DO K=KT,KB(I)
c                IF (ISO_CONC(CN(JC)))  C2(K,I,CN(JC)) = CIC(CN(JC))
c                IF (VERT_CONC(CN(JC))) C2(K,I,CN(JC)) = CVP(K,CN(JC))
                c2(k,i,cn(jc)) = cicbr(cn(jc),jb)
                C1(K,I,CN(JC)) = C2(K,I,CN(JC))
              END DO
            END DO
          END DO

********* Constituent mass

          DO JC=1,NAC
            JAC           = CN(JC)
            ICMBR(JAC,JB) = 0.0
            DO I=IU,ID
              ICMBR(JAC,JB) = ICMBR(JAC,JB)+C2(KT,I,JAC)*DLX(I)
     .                        *BHKT2(I)
              DO K=KT+1,KB(I)
                ICMBR(JAC,JB) = ICMBR(JAC,JB)+C2(K,I,JAC)*DLX(I)
     .                          *BH(K,I)
              END DO
            END DO
          END DO

********* Ice cover

          DO I=IU,ID
            ICETH(I) = IICETH
            ICE(I)   = ICETH(I).GT.0.0
          END DO

********* Vertical eddy viscosity

          IUT = IU
          IDT = ID-1
          IF (UP_HEAD(JB)) IUT = IU-1
          IF (DN_HEAD(JB)) IDT = ID
          DO I=IUT,IDT
            DO K=KT,KBMIN(I)-1
              AZ(K,I) = AZMIN
            END DO
          END DO
        END IF

******* Saved concentrations

        DO JC=1,NAC
          DO I=IU,ID
            DO K=KT,KB(I)
              C1S(K,I,CN(JC)) = C1(K,I,CN(JC))
            END DO
          END DO
        END DO

******* Density

        DO I=IU,ID
          DO K=KT,KB(I)
            ttt=T2(K,I)
            RHO(K,I) = DENSITY (ttt,SS(K,I),DISS(K,I))
          END DO
        END DO

******* Horizontal diffusion

        DO I=IU,ID-1
          DO K=KT,KBMIN(I)
            DX(K,I) = IDX
          END DO
        END DO
      END DO
      IF (OPEN_LPR) CLOSE (LPR)
      CALL GREGORIAN_DATE (YEAR)

************************************************************************
**                    Task 5: Initialize Output                       **
************************************************************************

      IF (SNAPSHOT)    OPEN (SNP,FILE=SNPFN,STATUS='UNKNOWN')
      IF (TIME_SERIES) OPEN (TSR,FILE=TSRFN,STATUS='UNKNOWN')
      IF (VECTOR)      OPEN (VPL,FILE=VPLFN,STATUS='UNKNOWN')
      IF (PROFILE)     OPEN (PRF,FILE=PRFFN,STATUS='UNKNOWN')
      IF (CONTOUR)     OPEN (CPL,FILE=CPLFN,STATUS='UNKNOWN')

***** Plot files

      IF (.NOT.RESTART_IN) THEN
        IF (PROFILE) THEN
          WRITE (PRF,2500)  TITLE
          WRITE (PRF,2510) (CPRN(JC),JC=1,NCP),' ON'
          WRITE (PRF,2520) (CNAME(JC),CUNIT(JC),JC=1,NCP),
     .                     'Temperature     ','deg C  '
          WRITE (PRF,2530)  CONSTITUENTS,NAC+1,(CN(JC),JC=1,NAC),22
          WRITE (PRF,2540)  PRFDP,KT,(KB(IPRF(I)),I=1,NIPRF)
          WRITE (PRF,2550)  H
          DO JC=1,NAC
            IF (CPRN(CN(JC)).EQ.' ON') THEN
              DO JPRF=1,NIPRF
                I   = IPRF(JPRF)
                NRS = KB(I)-KT+1
                WRITE (PRF,2560) CN(JC),NRS,(C2(K,I,CN(JC)),K=KT,KB(I))
              END DO
            END IF
          END DO
          DO JPRF=1,NIPRF
            I   = IPRF(JPRF)
            NRS = KB(I)-KT+1
            WRITE (PRF,2560) 22,NRS,(T2(K,I),K=KT,KB(I))
          END DO
        END IF
        IF (TIME_SERIES) THEN
          WRITE (TSR,5000) TITLE
          WRITE (TSR,5010) NAC,(CN(JC),JC=1,NAC)
          WRITE (TSR,5020) (CNAME(CN(JC)),CUNIT(CN(JC)),JC=1,NAC)
        END IF
        IF (CONTOUR) THEN
          WRITE (CPL,*) TITLE
          WRITE (CPL,*) DLX
          WRITE (CPL,*) DATUM,KMP,NBP,H,KB,DS,(CN(JC),JC=1,NAC)
        END IF
        IF (VECTOR) THEN
          WRITE (VPL,*) TITLE
          WRITE (VPL,*) H,KB,US,DS,DLX
        END IF
      END IF

************************************************************************
**                    Task 6: Begin Simulation                        **
************************************************************************
c        open(111,file='gsmcdd95.dat',status='old')
c        nxgsmcdd=0.0
        iflag1=0
        iflag2=0
        pumpon=.false.
        pumpon4=.false.
c        pumpinc=.false.
c        pumpdec=.false.
        qtr61=0.0
        notflowaug=.false.
      DO WHILE (.NOT.END_RUN)

c change boundary condition at fairview lake to an open system on day 293
c only valid for year 1992

      if(jday.ge.293.344.and.iflag1.eq.0)then
      iflag1=1
      dn_flow(8)=.false.
      dn_head(8)=.true.
      dh_external(8)=.true.
      head_boundary=.true.
      id=ds(8)
      kb(id+1)=kb(id)
            KBMIN(id) = MIN(KB(id),KB(id+1))
            DO K=KT,KBMIN(id)-1
              AZ(K,id) = AZMIN
            END DO
      end if

        IF (JDAY.GE.NXTVD)  CALL TIME_VARYING_DATA  (JDAY,NXTVD)
        IF (INTERPOLATE)   CALL INTERPOLATE_INPUTS (JDAY)
c scan water levels through the upper slough branches 4,5,6,7,8
c if any are above 10.50 ft, turn off pumping from columbia river, flow
c augmentation
        if(notflowaug)then
c if water level at cell 72 is at 9.9 ft turn it back on again
          if(((el(kt)-z2(72))*3.2808).le.9.9)then
          notflowaug=.false.
          qtr(61)=qtr61
          end if
        end if
        do jb=4,8
          iu=cus(jb)
          id=ds(jb)
          do i=iu,id
          if(((el(kt)-z2(i))*3.2808).ge.10.50)then
          if(notflowaug)then
          go to 1532
          else
          qtr61=qtr(61)
          end if
1532      notflowaug=.true.
          qtr(61)=0.0
          go to 1533
          end if
          end do
        end do
1533    continue

c      IF (JDAY.GE.NXgsmcdd) THEN
c21045   IF (JDAY.GE.NXgsmcdd) THEN
c          READ (111,'(f8.0)') NXgsmcdd
c          GO TO 21045
c        END IF
c        BACKSPACE 111
c        BACKSPACE 111
c        READ (111,'(2f8.0)') DUMMY2,gstemp
c      END IF
c      gsmcdd=gstemp*0.3048



c computing pumping outflow to lower slough at mcdd1 pump station
        if((el(kt)-z2(214)).gt.elpmpon)then
        pumpon=.true.
        qwd(29)=qpump
        end if
        if(pumpon)then
        if((el(kt)-z2(214)).lt.elpmpof)then
        pumpon=.false.
        qwd(29)=0.0
        end if
        end if

c regulating pumping outflow to lower slough at mcdd1 pump station
c        qwd(29)=qpstor
c	maxqpump=8.494
c        if((el(kt)-z2(214)).gt.(gsmcdd+0.3*0.3048))then
c          pumpinc=.true.
c        end if
c        if((el(kt)-z2(214)).lt.gsmcdd)then
c            pumpinc=.false.
c        end if
c        if(pumpinc)then
c          if((el(kt)-z2(214)).gt.(gsmcdd+0.3*0.3048))then	
c           qwd(29)=qpstor*1.2
c          end if
c          if((el(kt)-z2(214)).gt.(gsmcdd+0.6*0.3048))then
c           qwd(29)=qpstor*2.0
c          end if
c regulating the maximum pumping rate = 300 cfs
c		  if(qwd(29).gt.maxqpump) qwd(29)=8.494
c          if(qwd(29).eq.0.0) qwd(29)=qpump
c        end if
c
c        if((el(kt)-z2(214)).lt.(gsmcdd-0.3*0.3048))then
c          pumpdec=.true.
c        end if
c        if((el(kt)-z2(214)).gt.gsmcdd)then
c          pumpdec=.false.
c        end if
c        if((el(kt)-z2(214))/0.3048.lt.4.5)then
c          pumpdec=.true.
c        end if
c        if(pumpdec)then
c          qwd(29)=qpstor*0.5
c          if((el(kt)-z2(214)).lt.(gsmcdd-0.5*0.3048))then
c            qwd(29)=qpstor*0.2
c          end if
c          if((el(kt)-z2(214))/0.3048.lt.4.5)then
c            qwd(29)=0.0
c          end if
c        end if

c compute pumping outflow to Columbia at mcdd4 pump station
        if((el(kt)-z2(221)).gt.elpon4)then
        pumpon4=.true.
        qwd(40)=qpump4
        end if
        if(pumpon4)then
        if((el(kt)-z2(221)).lt.elpof4)then
        pumpon4=.false.
        qwd(40)=0.0
        end if
        end if

        if(nit.eq.0)then
        open(158,file='test.opt',status='new')
        end if

        if((nit/10)*10.eq.nit)then
          write(158,2601)jday,qwd(40),3.2808*(el(kt)-z2(221)),
     .      pumpon4
 2601     format(f9.4,2f8.2,L9)
        end if


c regulating pumping outflow to Columbia at mcdd4 pump station
c        if((el(kt)-z2(221)).gt.(gsmcdd+0.5*0.3048))then
c          pumpon4=.true.
c        end if
c        if((el(kt)-z2(221)).lt.(gsmcdd-0.3*0.3048) )then
c            pumpon4=.false.
c        end if
c        if((el(kt)-z2(221))/0.3048.lt.4.5)then
c            pumpon4=.false.
c        end if
c        if(pumpon4)then
c          qwd(40)=qpump4
c        else
c          qwd(40)=0.0
c        end if

c set irrigation withdrawals based on time of day/season
        qwd(30)=0.0
        qwd(31)=0.0
        qwd(32)=0.0
        qwd(33)=0.0
        qwd(34)=0.0
        qwd(35)=0.0
        qwd(36)=0.0
        qwd(37)=0.0
        qwd(38)=0.0
        qwd(39)=0.0
c following line is for code debugging only, set test='OFF' for production runs
      if(test)go to 1549
        if(jday.le.304.and.jday.ge.121)then
           today=jday-int(jday)
               if(today.ge.0.0.and.today.lt.0.17)then
                 qwd(30)=qwd1
                 qwd(31)=qwd2
                 qwd(32)=qwd3+qwd7
                 qwd(33)=qwd4
                 qwd(34)=qwd5
                 qwd(35)=qwd6
                 qwd(36)=qwd8
                 qwd(37)=qwd9+qwd10
                 qwd(38)=qwd11
                 qwd(39)=qwd13
               end if
               if(today.ge.0.17.and.today.lt.0.25)then
                 qwd(30)=qwd1
                 qwd(31)=qwd2
                 qwd(32)=qwd7
                 qwd(33)=qwd4
                 qwd(35)=qwd6
                 qwd(37)=qwd9
                 qwd(38)=qwd11
                 qwd(39)=qwd13
               end if
               if(today.ge.0.25.and.today.lt.0.27)then
                 qwd(30)=qwd1
                 qwd(32)=qwd7
                 qwd(33)=qwd4
                 qwd(35)=qwd6
                 qwd(38)=qwd11
                 qwd(39)=qwd13
               end if
               if(today.ge.0.27.and.today.lt.0.33)then
                 qwd(32)=qwd7
                 qwd(33)=qwd4
                 qwd(35)=qwd6
                 qwd(38)=qwd11
                 qwd(39)=qwd13
               end if
               if(today.ge.0.33.and.today.lt.0.42)then
                 qwd(38)=qwd11
               end if
         end if
c
c compute groundwater inflow distributed among the branches based on
c branch surface area
c        elmcdd1=(el(kt)-z2(2))*3.2808
c compute average GW elev at cells 119, 116, 154, 163, and 209
c        elavg=el(kt)-(z2(119)+z2(116)+z2(154)+z2(163)+z2(209))/5.
c        elavg=elavg*3.2808
c        qgw=(34514.4*elmcdd1**(-3.018))/35.313
c         qgw=(111.158-6.512*elavg)/35.313
c        qdtr(31)=gw1(31)*(el(kt)-z2(209))+gw2(31)
c        qdtr(30)=gw1(30)*(el(kt)-z2(204))+gw2(30)
c        qdtr(27)=gw1(27)*(el(kt)-z2(189))+gw2(27)
c        qdtr(25)=gw1(25)*(el(kt)-z2(178))+gw2(25)
c        qdtr(15)=gw1(15)*(el(kt)-z2(127))+gw2(15)
c        qdtr(16)=gw1(16)*(el(kt)-z2(132))+gw2(16)
c        qdtr(12)=gw1(12)*(el(kt)-z2(109))+gw2(12)
c        qdtr(11)=gw1(11)*(el(kt)-z2(104))+gw2(11)
c        qdtr(17)=gw1(17)*(el(kt)-z2(137))+gw2(17)
c calculating groundwater inflow of branches - chris
        do ig=1,nbp
          if(dist_tribs(ig))then
           wsurf=el(kt)-z2(us(ig))
           qdtr(ig)=gw1(ig)*wsurf+gw2(ig)
           if(qdtr(ig).lt.0.0)qdtr(ig)=0.0
          end if
        end do
c compute flow between branches
c flow from mcdd1 gravity structure
1549    mup=1
        mdn=100
        el1=el(kt)-z2(cus(mup))
        call qculv(el1,el2,qflow,jday,mup,mdn)
        qin(1)=0.0
        qout(1,8)=0.0
        qout(1,31)=0.0
        qout(1,12)=0.0
        qout(2,8)=0.0
        qout(2,31)=0.0
        qout(2,12)=0.0
        qwd(1)=qflow
c loop through the number of culverts #=38 but only 26 between branches
10010   CONTINUE
        do j=1,26
       if(j.eq.1)then
        mup=2
        mdn=1
        mupc=29
        mdnc=26
        mtr=1
        mwd=2
        mqin=2
        mout=1
        go to 1523
       end if
       if(j.eq.2)then
        mup=3
        mdn=2
        mupc=46
        mdnc=43
        mtr=2
        mwd=3
        mqin=3
        mout=2
        go to 1523
       end if
       if(j.eq.3)then
        mup=4
        mdn=3
        mupc=56
        mdnc=53
        mtr=3
        mwd=4
        mqin=4
        mout=3
        go to 1523
       end if
       if(j.eq.4)then
        mup=5
        mdn=4
        mupc=61
        mdnc=58
        mtr=4
        mwd=5
        mqin=5
        mout=4
        go to 1523
       end if
       if(j.eq.5)then
        mup=6
        mdn=5
        mupc=67
        mdnc=64
        mtr=5
        mwd=6
        mqin=6
        mout=5
        go to 1523
       end if
       if(j.eq.6)then
        mup=7
        mdn=6
        mupc=81
        mdnc=78
        mtr=6
        mwd=7
        mqin=7
        mout=6
        go to 1523
       end if
       if(j.eq.7)then
        mup=8
        mdn=7
        mupc=89
        mdnc=86
        mtr=7
        mwd=8
        mqin=8
        mout=7
        go to 1523
       end if
       if(j.eq.8)then
        mup=10
        mdn=9
        mupc=99
        mdnc=96
        mtr=9
        mwd=9
        mqin=10
        mout=9
        go to 1523
       end if
       if(j.eq.9)then
        mup=11
        mdn=10
        mupc=104
        mdnc=101
        mtr=10
        mwd=10
        mqin=11
        mout=10
        go to 1523
       end if
       if(j.eq.10)then
        mup=12
        mdn=11
        mupc=109
        mdnc=106
        mtr=11
        mwd=11
        mqin=12
        mout=11
        go to 1523
       end if
       if(j.eq.11)then
        mup=14
        mdn=13
        mupc=119
        mdnc=116
        mtr=13
        mwd=12
        mqin=14
        mout=13
        go to 1523
       end if
       if(j.eq.12)then
        mup=15
        mdn=14
        mupc=127
        mdnc=120
        mtr=15
        mwd=14
        mqin=15
        go to 1523
       end if
       if(j.eq.13)then
        mup=16
        mdn=15
        mupc=132
        mdnc=129
        mtr=16
        mwd=15
        mqin=16
        mout=15
        go to 1523
       end if
       if(j.eq.14)then
        mup=17
        mdn=14
        mupc=137
        mdnc=124
        mtr=14
        mwd=16
        mqin=17
        mout=14
        go to 1523
       end if
       if(j.eq.15)then
        mup=18
        mdn=17
        mupc=143
        mdnc=140
        mtr=18
        mwd=17
        mqin=18
        mout=17
        go to 1523
       end if
       if(j.eq.16)then
        mup=19
        mdn=18
        mupc=148
        mdnc=145
        mtr=19
        mwd=18
        mqin=19
        mout=18
        go to 1523
       end if
       if(j.eq.17)then
        mup=20
        mdn=19
        mupc=153
        mdnc=150
        mtr=20
        mwd=19
        mqin=20
        mout=19
        go to 1523
       end if
       if(j.eq.18)then
        mup=21
        mdn=20
        mupc=158
        mdnc=155
        mtr=21
        mwd=20
        mqin=21
        mout=20
        go to 1523
       end if
       if(j.eq.19)then
        mup=22
        mdn=21
        mupc=163
        mdnc=160
        mtr=22
        mwd=21
        mqin=22
        mout=21
        go to 1523
       end if
       if(j.eq.20)then
        mup=23
        mdn=22
        mupc=168
        mdnc=165
        mtr=23
        mwd=22
        mqin=23
        mout=22
        go to 1523
       end if
       if(j.eq.21)then
        mup=24
        mdn=23
        mupc=173
        mdnc=170
        mtr=24
        mwd=23
        mqin=24
        mout=23
        go to 1523
       end if
       if(j.eq.22)then
        mup=26
        mdn=24
        mupc=184
        mdnc=175
        mtr=25
        mwd=24
        mqin=26
        mout=24
        go to 1523
       end if
       if(j.eq.23)then
        mup=27
        mdn=26
        mupc=189
        mdnc=186
        mtr=27
        mwd=25
        mqin=27
        mout=26
        go to 1523
       end if
       if(j.eq.24)then
        mup=29
        mdn=28
        mupc=199
        mdnc=196
        mtr=28
        mwd=26
        mqin=29
        mout=28
        go to 1523
       end if
       if(j.eq.25)then
        mup=30
        mdn=29
        mupc=204
        mdnc=201
        mtr=29
        mwd=27
        mqin=30
        mout=29
        go to 1523
       end if
       if(j.eq.26)then
        mup=31
        mdn=30
        mupc=209
        mdnc=206
        mtr=30
        mwd=28
        mqin=31
        mout=30
        go to 1523
       end if
c calibration run no flow across cross-levee
1523    if(j.eq.3.and.jday.ge.278.788.and.jday.lt.284.292)then
        qflow=0.0
        go to 1529
        end if
c1523    continue

        el1=el(kt)-z2(mupc)
        el2=el(kt)-z2(mdnc)

        call qculv(el1,el2,qflow,jday,mup,mdn)
        if(qflow.gt.0.0.and.qin(mqin).ne.0.0)then
c flow oscillation, set q to 0.0, lags everything by 1 dt
        qflow=0.0
        end if
        if(qflow.lt.0.0.and.qwd(mwd).ne.0.0)then
        qflow=0.0
        end if
c
1529          if(qflow.gt.0)then
c calculate maximum flow allowed in time step
c must have a smooth transition in q from time step to time step
                ell1=dlx(mupc)*b(ktt(mupc),mupc)
                ell2=dlx(mdnc)*b(ktt(mdnc),mdnc)
                dzz2=qflow*dlt/(ell2)
                dzz1=qflow*dlt/(ell1)
                if((dzz2+el2).gt.(el1-dzz1))then
                equilz=((el2*ell2)+(el1*ell1))/
     .           (ell1+ell2)
                 dzz=0.95*(el1-equilz)
                 qflow=(ell1*dzz/dlt)
                end if

            qwd(mwd)=qflow
            if(mqin.eq.15)qwd(mwd-1)=0
            qtr(mtr)=qflow
            if(mqin.ne.15)qout(1,mout)=0.0
            if(mqin.ne.15)qout(2,mout)=0.0
            qin(mqin)=0.0
          else
                if(qflow.ne.0.0)then
                ell1=dlx(mupc)*b(ktt(mupc),mupc)
                ell2=dlx(mdnc)*b(ktt(mdnc),mdnc)
                dzz2=qflow*dlt/(ell2)
                dzz1=qflow*dlt/(ell1)
                if((dzz1+el1).gt.(el2-dzz2))then
                equilz=((el2*ell2)+(el1*ell1))/
     .           (ell1+ell2)
                 dzz=0.95*(equilz-el1)
                 qflow=-(ell1*dzz/dlt)
                end if
                end if

            qwd(mwd)=0.0
            qtr(mtr)=0.0
            if(mqin.eq.15)qwd(mwd-1)=-qflow
              if(mqin.ne.15)then
                 if(kt.le.5)then
c note that this outlet is distributed around cells 5(=1) and 6(=2)
                    qout(1,mout)=-qflow/2.
                    qout(2,mout)=-qflow/2.
                 else
                    qout(1,mout)=0.0
                    qout(2,mout)=-qflow
                 end if
               qin(mqin)=-qflow
               end if
           end if
        ttr(mtr)=t2(kt,mupc)
        tin(mqin)=t2(kt,mdnc)
                  if(constituents)then
                  do jc=1,nac
                  ctr(cn(jc),mtr)= c2(kt,mupc,cn(jc))
                  cin(cn(jc),mqin)= c2(kt,mdnc,cn(jc))
                  end do
                  end if
        end do

************************************************************************
**            Task 6.1: Update Hydrodynamic Sources/Sinks             **
************************************************************************

******* Entry point for timestep violation

        DO JB=1,NBP
          IF (SEL_WITHDRAWAL(JB)) THEN
            KTOPSW(JB) = KMP
            KBOTSW(JB) = 2
            CALL SELECTIVE_WITHDRAWAL
          END IF
        END DO
c10010   CONTINUE
        DO JB=1,NBP
          IU = CUS(JB)
          ID = DS(JB)
          IF (EVAPORATION) THEN
            EVBR(JB) = 0.0
            FW       = 9.2+0.46*WIND*WIND
            DO I=IU,ID
              TM    = (T2(KT,I)+TDEW)*0.5
              VPTG  = 0.35+0.015*TM+0.0012*TM*TM
              EV(I) = VPTG*(T2(KT,I)-TDEW)*FW*B(KT,I)*DLX(I)/2.45E9
              IF (EV(I).LT.0.0.OR.ICE(I)) EV(I) = 0.0
              QSS(KT,I) = QSS(KT,I)-EV(I)
              EVBR(JB)  = EVBR(JB)+EV(I)
            END DO
          END IF
          IF (PRECIPITATION) THEN
            QPRBR(JB) = 0.0
            DO I=IU,ID
              QPR(I)    = PR(JB)*B(KTT(I),I)*DLX(I)
              QPRBR(JB) = QPRBR(JB)+QPR(I)
              QSS(KT,I) = QSS(KT,I)+QPR(I)
            END DO
          END IF
          IF (TRIBUTARIES) THEN
            DO JT=1,NTP
              IF (JB.EQ.JBTR(JT)) THEN
                I = ITR(JT)
                IF (I.LT.IU) I = IU
                IF (PLACE_QTR) THEN
                  K       = KT
                  KTR(JT) = K
                  ttt=TTR(JT)
                  RHOIN   = DENSITY (ttt,CTR(2,JT),CTR(4,JT))
                  DO WHILE (RHOIN.GT.RHO(K,I).AND.KTR(JT).LT.KB(I))
                    K       = K+1
                    KTR(JT) = K
                  END DO
                ELSE
                  BHSUM       = BHKT2(I)
                  QTRF(KT,JT) = 0.0
                  DO K=KT+1,KB(I)
                    BHSUM      = BHSUM+BH(K,I)
                    QTRF(K,JT) = 0.0
                  END DO
                  QTRF(KT,JT) = BHKT2(I)/BHSUM
                  DO K=KT+1,KB(I)
                    QTRF(K,JT) = BH(K,I)/BHSUM
                  END DO
                END IF
                DO K=KT,KB(I)
                  QSS(K,I) = QSS(K,I)+QTR(JT)*QTRF(K,JT)
                END DO
              END IF
            END DO
          END IF
          IF (DIST_TRIBS(JB)) THEN
            DO I=IU,ID
              QDT(I)    = QDTR(JB)*B(KT,I)*DLX(I)/AKBR(KT,JB)
              QSS(KT,I) = QSS(KT,I)+QDT(I)
            END DO
          END IF
          IF (WITHDRAWALS) THEN
          if(kt.le.5)qf=4.
          if(kt.eq.6)qf=3.
          if(kt.eq.7)qf=2.
            DO JW=1,NWP
              IF (JB.EQ.JBWD(JW)) THEN
                IF (KWD(JW).GE.KT) THEN
                  QSS(KWD(JW),IWD(JW))=QSS(KWD(JW),IWD(JW))-QWD(JW)/qf
c note that the kwd for all these is cell 8, below this outflow
c is spread among cells 5, 6, and 7 also
                  if(5.ge.kt)QSS(5,IWD(JW))=QSS(5,IWD(JW))-QWD(JW)/qf
                  if(6.ge.kt)QSS(6,IWD(JW))=QSS(6,IWD(JW))-QWD(JW)/qf
                  if(7.ge.kt)QSS(7,IWD(JW))=QSS(7,IWD(JW))-QWD(JW)/qf
                END IF
              END IF
            END DO
          END IF
          IF (UH_INTERNAL(JB)) THEN
            DO K=KT,KB(IU-1)
              QSS(K,UHS(JB)) = QSS(K,UHS(JB))-QUH2(K,JB)/DLT
            END DO
          END IF
          IF (DH_INTERNAL(JB)) THEN
            DO K=KT,KB(ID+1)
              QSS(K,DHS(JB)) = QSS(K,DHS(JB))+QDH2(K,JB)/DLT
            END DO
          END IF
        END DO

************************************************************************
**               Task 6.2: Hydrodynamic Calculations                  **
************************************************************************

        DO JB=1,NBP
          IU = CUS(JB)
          ID = DS(JB)

************************************************************************
**                    Task 6.2.1: Momentum Terms                      **
************************************************************************

********* Set boundary concentrations, temperatures, and densities

          IF (.NOT.VIOLATION) THEN
            IUT = IU
            IDT = ID
            IF (UP_FLOW(JB)) THEN
              DO K=KT,KB(IU)
                DO JC=1,NAC
                  C1(K,IU-1,CN(JC))  = CIN(CN(JC),JB)
                  C1S(K,IU-1,CN(JC)) = CIN(CN(JC),JB)
                END DO
                T1(K,IU-1) = TIN(JB)
                T2(K,IU-1) = TIN(JB)
c below are new lines of code
               if(qin(jb).eq.0.0)then
              T1(K,IU-1) = T1(K,IU)
              T2(K,IU-1) = T1(K,IU)
               end if

              END DO
            END IF
            IF (DN_FLOW(JB)) THEN
              DO K=KT,KB(IU)
                DO JC=1,NAC
                  C1(K,ID+1,CN(JC))  = C1S(K,ID,CN(JC))
                  C1S(K,ID+1,CN(JC)) = C1S(K,ID,CN(JC))
                END DO
                T1(K,ID+1) = T2(K,ID)
                T2(K,ID+1) = T2(K,ID)
              END DO
            END IF
            IF (UP_HEAD(JB)) THEN
              IUT = IU-1
              IF (UH_INTERNAL(JB)) THEN
                DO K=KT,KB(IU-1)
                  DO JC=1,NAC
                    C1(K,IUT,CN(JC)) = C2(K,UHS(JB),CN(JC))
                    C2(K,IUT,CN(JC)) = C2(K,UHS(JB),CN(JC))
                  END DO
                  T1(K,IUT)  = T2(K,UHS(JB))
                  T2(K,IUT)  = T2(K,UHS(JB))
                  RHO(K,IUT) = RHO(K,UHS(JB))
                END DO
              ELSE IF (UH_EXTERNAL(JB)) THEN
                DO K=KT,KB(IU-1)
                  DO JC=1,NAC
                    C1(K,IUT,CN(JC)) = CUH(K,CN(JC),JB)
                    C2(K,IUT,CN(JC)) = CUH(K,CN(JC),JB)
                  END DO
                  T1(K,IUT)  = TUH(K,JB)
                  T2(K,IUT)  = TUH(K,JB)
                  ttt=T2(K,IUT)
                  RHO(K,IUT) = DENSITY (ttt,SS(K,IUT),DISS(K,IUT))
                END DO
              END IF
            END IF
            IF (DN_HEAD(JB)) THEN
              IDT = ID+1
              IF (DH_INTERNAL(JB)) THEN
                DO K=KT,KB(ID+1)
                  DO JC=1,NAC
                    C1(K,IDT,CN(JC)) = C2(K,DHS(JB),CN(JC))
                    C2(K,IDT,CN(JC)) = C2(K,DHS(JB),CN(JC))
                  END DO
                  T1(K,IDT)  = T2(K,DHS(JB))
                  T2(K,IDT)  = T2(K,DHS(JB))
                  RHO(K,IDT) = RHO(K,DHS(JB))
                END DO
              ELSE IF (DH_EXTERNAL(JB)) THEN
                DO K=KT,KB(ID+1)
                  DO JC=1,NAC
                    C1(K,IDT,CN(JC)) = CDH(K,CN(JC),JB)
                    C2(K,IDT,CN(JC)) = CDH(K,CN(JC),JB)
                  END DO
                  T1(K,IDT)  = TDH(K,JB)
                  T2(K,IDT)  = TDH(K,JB)
                  ttt=T2(K,IDT)
                  RHO(K,IDT) = DENSITY (ttt,SS(K,IDT),DISS(K,IDT))
                END DO
              END IF
            END IF

*********** Densities

            DO I=IU,ID
              DO K=KT,KB(I)
                ttt=T2(K,I)
                RHO(K,I) = DENSITY (ttt,SS(K,I),DISS(K,I))
              END DO
            END DO

*********** Density pressures

            DO I=IUT,IDT
              if(test3)then
              P(KT,I) = RHO(KT,I)*G*HKT2(I)
              else
              P(KT,I) = RHO(KT,I)*G*H(KT)
              end if
              DO K=KT+1,KB(I)
                P(K,I) = P(K-1,I)+RHO(K,I)*G*H(K)
              END DO
            END DO

*********** Horizontal density gradients

            DO I=IUT,IDT-1
              DLXRHO    = 1.0/((DLX(I)+DLX(I+1))*RHOI)
              HPG(KT,I) = DLXRHO*(BKT(I)+BKT(I+1))*0.5*(HKT2(I+1)
     .                    *P(KT,I+1)-HKT2(I)*P(KT,I))
              DO K=KT+1,KBMIN(I)
                HPG(K,I) = DLXRHO*BHR(K,I)*((P(K-1,I+1)-P(K-1,I))
     .                     +(P(K,I+1)-P(K,I)))
              END DO
            END DO
          END IF

********* Vertical eddy viscosity

c done in tvds          WINDS = WSC(WSCDP)*WIND
          CZ    = 0.0
          IF (WIND.GE.1.0)  CZ = 0.0005*SQRT(WIND)
          IF (WIND.GE.15.0) CZ = 0.0026
          SSC    = RHOA*CZ*WIND**2/RHOI
          SSCCOS = SSC*COS(PHI-PHI0(JB))
          SSCSIN = SSC*ABS(SIN(PHI-PHI0(JB)))
          DO I=IUT,IDT-1
            FETCH = FETCHD(I,JB)
            IF (COS(PHI-PHI0(JB)).LT.0.0) FETCH = FETCHU(I,JB)
            WWT = 0.0
            IF (WIND.NE.0.0) THEN
              WWT = 6.95E-2*(FETCH**0.233)*ABS(WIND)**0.534
            END IF
            DFC        = -8.0*3.14159**2/(G*WWT*WWT+1.0E-20)
            DEPTH      = HKT2(I)
            EXPDF      = EXP(MAX(DFC*DEPTH,-20.0))
            SSCCOS     = ICESW(I)*SSCCOS
            SSCSIN     = ICESW(I)*SSCSIN
c br must be for the avearge upper cell
c            ST(KT,I)   = SSCCOS*(BR(KT-1,I)+BR(KT,I))*0.5
             brkt=(BKT(I)+BKT(I+1))*0.5
             bbkt=(b(ktt(i),i)+b(ktt(i+1),i+1))*0.5
             ST(KT,I)   = SSCCOS*bbkt
C STILL SOME DISCUSSION ON ABOVE LINE, tom thinks it should be
C            ST(KT,I)=SSCOS*(BR(KTT(I)-1,I)+BR(KT,I))*0.5
c next line is OK
            ST(KT+1,I) = SSCCOS*EXPDF*(BR(KT,I)+BR(KT+1,I))*0.5
C            ST(KT+1,I) = SSCCOS*EXPDF*(brkt+BR(KT+1,I))*0.5
            SHEARS     = ((U(KT+1,I)-U(KT,I))/AVHKT(I))**2
            if(az(kt,i).eq.0.0)az(kt,i)=azmin
            AZ0        = 0.1*AVHKT(I)**2*SQRT(SHEARS+(SSCSIN
     .                   *EXPDF/AZ(KT,I))**2)+AZMIN
            AZMXKT(I) = 0.1*HKT2(I)*HKT2(I)/DLT
            RIAZ0      = LOG(AZ0/AZMXKT(I))/1.5
            BUOY       = (RHO(KT+1,I)-RHO(KT,I)+RHO(KT+1,I+1)
     .                   -RHO(KT,I+1))/(2.0*AVHKT(I))
            RI         = G*BUOY/(RHOI*SHEARS+1.0E-10)
            RIAZ1      = MAX(RI,RIAZ0)
            RIAZ1      = MIN(RIAZ1,10.0)
            EXPRAZ     = EXP(-1.5*RIAZ1)
            AZ(KT,I)   = MAX(AZMIN,AZ0*EXPRAZ+AZMIN*(1.0-EXPRAZ))
            RI         = MIN(RI, 10.0)
            RI         = MAX(RI,-10.0)
            EXPRDZ     = EXP(-1.5*RI)
            DZ(KT,I)   = MAX(DZMIN,FRAZDZ*(AZ0*EXPRDZ+DZMIN
     .                   *(1.0-EXPRDZ)))
            KBT        = KBMIN(I)
            DO K=KT+2,KBT
              DEPTH      = DEPTH+H(K-1)
              AZMAX(K-1) = 0.1*H(K)*H(K)/DLT
              EXPDF      = EXP(MAX(DFC*DEPTH,-20.0))
              ST(K,I)    = SSCCOS*EXPDF*(BR(K-1,I)+BR(K,I))*0.5
              SHEARS     = ((U(K,I)-U(K-1,I))/AVH(K-1))**2
              AZ0        = 0.1*AVH(K-1)**2*SQRT(SHEARS+(SSCSIN
     .                     *EXPDF/AZ(K-1,I))**2)+AZMIN
              RIAZ0      = LOG(AZ0/AZMAX(K-1))/1.5
              BUOY       = (RHO(K,I)-RHO(K-1,I)+RHO(K,I+1)
     .                      -RHO(K-1,I+1))/(2.0*AVH(K-1))
              RI         = G*BUOY/(RHOI*SHEARS+0.000001)
              RIAZ1      = MAX(RI,RIAZ0)
              RIAZ1      = MIN(RIAZ1,10.0)
              EXPRAZ     = EXP(-1.5*RIAZ1)
              AZ(K-1,I)  = MAX(AZMIN,AZ0*EXPRAZ+AZMIN*(1.0-EXPRAZ))
              RI         = MIN(RI, 10.0)
              RI         = MAX(RI,-10.0)
              EXPRDZ     = EXP(-1.5*RI)
              DZ(K-1,I)  = MAX(DZMIN,FRAZDZ*(AZ0*EXPRDZ+DZMIN
     .                     *(1.0-EXPRDZ)))
            END DO
            SB(KBT,I) = SSCCOS*EXPDF*(BR(KBT-1,I)+BR(KBT,I))*0.5
          END DO
          DO K=KT,KB(IDT)-1
            DZ(K,IDT) = DZ(K,IDT-1)
          END DO

********* Check for density inversions

          DO I=IUT,IDT
            DO K=KT,KB(I)-1
              IF (RHO(K,I).GT.RHO(K+1,I)) DZ(K,I) = DZMAX
            END DO
          END DO

********* Shear stresses

          DO I=IUT,IDT-1
            KBT = KBMIN(I)
            ST(KT+1,I)=ST(KT+1,I)+AZ(KT,I)*(BR(KT,I)+BR(KT+1,I))*0.5
     .                  *(U(KT,I)-U(KT+1,I))/AVHKT(I)
            DO K=KT+2,KBT
              ST(K,I) = ST(K,I)+AZ(K-1,I)*(BR(K-1,I)+BR(K,I))*0.5
     .                  *(U(K-1,I)-U(K,I))/AVH(K-1)
            END DO
            DO K=KT,KBT-1
c Manning's approach instead of Chezy
c Chezy coef=(Rh^1/6)/n
c              SB(K,I) = ST(K+1,I)+G/(CHEZY*CHEZY)*(BR(K,I)-BR(K+1,I))
c     .                  *U(K,I)*ABS(U(K,I))
c            END DO
              if(k.eq.kt)then
c              brkt=(bkt(i)+bkt(i+1))*0.5
        rh=bhrkt2(i)/(hkt2(i)+hkt2(i+1)+br(ktt(i),i)-br(k+1,i))
        SB(K,I)=ST(K+1,I)+G*sman(i)*sman(i)/(rh**0.333)
     .                  *(br(ktt(i),i)-BR(K+1,I))
     .                  *U(K,I)*ABS(U(K,I))
              end if
              if(k.ne.kt)then
              rh=bhr(k,i)/(2.*h(k)+br(k,i)-br(k+1,i))
        SB(K,I)=ST(K+1,I)+G*sman(i)*sman(i)/(rh**0.333)
     .                  *(BR(K,I)-BR(K+1,I))
     .                  *U(K,I)*ABS(U(K,I))
              end if
            END DO
        rh=bhr(kbt,i)/(2.*h(kbt)+br(kbt,i))
        SB(KBT,I) = SB(KBT,I)+G*(sman(i)*sman(i))/(rh**0.333)
     .                 *BR(KBT,I)*U(KBT,I)
     .                 *ABS(U(KBT,I))
          END DO
c            SB(KBT,I) = SB(KBT,I)+G/(CHEZY*CHEZY)*BR(KBT,I)*U(KBT,I)
c     .                  *ABS(U(KBT,I))
c          END DO

********* Horizontal momentum

          IF (.NOT.VIOLATION) THEN
            DO I=IU,ID-1
              AR         = (1.0+SIGN(1.0,(U(KT,I)+U(KT,I+1))*0.5))*0.5
              AL         = (1.0+SIGN(1.0,(U(KT,I)+U(KT,I-1))*0.5))*0.5
              ADMX(KT,I) = BHKT2(I+1)*(U(KT,I+1)+U(KT,I))*0.5
     .                     *(AR*U(KT,I)+(1.0-AR)*U(KT,I+1))/DLX(I+1)
     .                     -BHKT2(I)*(U(KT,I)+U(KT,I-1))*0.5
     .                     *(AL*U(KT,I-1)+(1.0-AL)*U(KT,I))/DLX(I)
              DM(KT,I)   = AX/DLX(I)*(2.0/(DLX(I+1)+DLX(I))*(BHRKT2(I+1)
     .                     *U(KT,I+1)-BHRKT2(I)*U(KT,I))-2.0/(DLX(I-1)
     .                     +DLX(I))*(BHRKT2(I)*U(KT,I)-BHRKT2(I-1)
     .                     *U(KT,I-1)))
              DO K=KT+1,KBMIN(I)
                AR        = (1.0+SIGN(1.0,(U(K,I)+U(K,I+1))*0.5))*0.5
                AL        = (1.0+SIGN(1.0,(U(K,I)+U(K,I-1))*0.5))*0.5
                ADMX(K,I) = BH(K,I+1)*(U(K,I+1)+U(K,I))*0.5*(AR*U(K,I)
     .                      +(1.0-AR)*U(K,I+1))/DLX(I+1)-BH(K,I)*(U(K,I)
     .                      +U(K,I-1))*0.5*(AL*U(K,I-1)+(1.0-AL)*U(K,I))
     .                      /DLX(I)
                DM(K,I)   = AX/DLX(I)*(2.0/(DLX(I+1)+DLX(I))*(BHR(K,I+1)
     .                      *U(K,I+1)-BHR(K,I)*U(K,I))-2.0/(DLX(I-1)
     .                      +DLX(I))*(BHR(K,I)*U(K,I)-BHR(K,I-1)
     .                      *U(K,I-1)))
              END DO
            END DO

*********** Vertical momentum

            DO I=IU,ID-1
              DO K=KT,KB(I)-1
                AB        = (1.0+SIGN(1.0,(W(K,I+1)+W(K,I))*0.5))*0.5
                ADMZ(K,I) = (BR(K,I)+BR(K+1,I))*0.5*(W(K,I+1)+W(K,I))
     .                      *0.5*(AB*U(K,I)+(1.0-AB)*U(K+1,I))
c              if(k.eq.kt)then
c              bbb=(bkt(i)+bkt(i+1))*0.5
c              else
c              bbb=br(k,i)
c              end if
c              ADMZ(K,I) = (bbb+BR(K+1,I))*0.5*(W(K,I+1)+W(K,I))*0.5
c     .                    *(AB*U(K,I)+(1.0-AB)*U(K+1,I))
              END DO
            END DO
          END IF

************************************************************************
**         Task 6.2.2: Compute the Water Surface Elevation            **
************************************************************************

********* Tridiagonal coefficients (D and F)

          DO I=IU,ID-1
            BHRHO(I) = BHKT2(I+1)/RHO(KT,I+1)+BHKT2(I)/RHO(KT,I)
            DO K=KT+1,KBMIN(I)
              BHRHO(I) = BHRHO(I)+(BH(K,I+1)/RHO(K,I+1)+BH(K,I)
     .                   /RHO(K,I))
            END DO
            D(I) = U(KT,I)*BHRKT2(I)-U(KT,I-1)*BHRKT2(I-1)-QSS(KT,I)
            F(I) = -SB(KT,I)+ST(KT,I)-ADMX(KT,I)+DM(KT,I)-HPG(KT,I)
            DO K=KT+1,KB(I)
              D(I) = D(I)+(U(K,I)*BHR(K,I)-U(K,I-1)*BHR(K,I-1)
     .               -QSS(K,I))
              F(I) = F(I)+(-SB(K,I)+ST(K,I)-ADMX(K,I)+DM(K,I)-HPG(K,I))
            END DO
          END DO
          D(IU) = U(KT,IU)*BHRKT2(IU)-QSS(KT,IU)
          DO K=KT+1,KB(IU)
            D(IU) = D(IU)+(U(K,IU)*BHR(K,IU)-QSS(K,IU))
          END DO

********* Boundary tridiagonal coefficients (D and F)

          IF (UP_FLOW(JB)) D(IU) = D(IU)-QIN(JB)
          IF (DN_FLOW(JB)) THEN
            D(ID) = -U(KT,ID-1)*BHRKT2(ID-1)-QSS(KT,ID)
            DO K=KT+1,KB(ID)
              D(ID) = D(ID)-(U(K,ID-1)*BHR(K,ID-1)+QSS(K,ID))
            END DO
            DO JO=1,NOUT(JB)
              IF (KOUT(JO,JB).GE.KT) THEN
                D(ID) = D(ID)+QOUT(JO,JB)
              END IF
            END DO
          END IF
          IF (UP_HEAD(JB)) THEN
            BHRHO(IU-1) = BHKT2(IU)/RHO(KT,IU)+BHKT2(IU-1)/RHO(KT,IU-1)
            DO K=KT+1,KB(IU-1)
              BHRHO(IU-1) = BHRHO(IU-1)+(BH(K,IU)/RHO(K,IU)+BH(K,IU-1)
     .                      /RHO(K,IU-1))
            END DO
            D(IU)   = D(IU)-U(KT,IU-1)*BHRKT2(IU-1)
            F(IU-1) = -SB(KT,IU-1)+ST(KT,IU-1)-HPG(KT,IU-1)
            DO K=KT+1,KB(IU)
              D(IU)   = D(IU)-U(K,IU-1)*BHR(K,IU-1)
              F(IU-1) = F(IU-1)-(SB(K,IU-1)-ST(K,IU-1)+HPG(K,IU-1))
            END DO
          END IF
          IF (DN_HEAD(JB)) THEN
            BHRHO(ID) = BHKT2(ID+1)/RHO(KT,ID+1)+BHKT2(ID)/RHO(KT,ID)
            DO K=KT+1,KB(ID+1)
              BHRHO(ID) = BHRHO(ID)+(BH(K,ID+1)/RHO(K,ID+1)+BH(K,ID)
     .                    /RHO(K,ID))
            END DO
            D(ID) = U(KT,ID)*BHRKT2(ID)-U(KT,ID-1)*BHRKT2(ID-1)
     .              -QSS(KT,ID)
            F(ID) = -SB(KT,ID)+ST(KT,ID)-HPG(KT,ID)
            DO K=KT+1,KB(ID)
              D(ID) = D(ID)+(U(K,ID)*BHR(K,ID)-U(K,ID-1)*BHR(K,ID-1)
     .                -QSS(K,ID))
              F(ID) = F(ID)+(-SB(K,ID)+ST(K,ID)-HPG(K,ID))
            END DO
          END IF
        END DO
        DO JB=1,NBP
          IU = CUS(JB)
          ID = DS(JB)

********* Set boundary surface elevations

          IF (UH_INTERNAL(JB)) Z1(IU-1) = Z1(UHS(JB))
          IF (UH_EXTERNAL(JB)) Z1(IU-1) = EL(KT)-ELUH(JB)
          IF (DH_INTERNAL(JB)) Z1(ID+1) = Z1(DHS(JB))
          IF (DH_EXTERNAL(JB))then
          Z1(ID+1) = EL(KT)-ELDH(JB)
c this code was added to smooth out the sudden transition of opening up the
c pipes from Fairview Lake into the Upper Slough
          diff=z1(id+1)-z1(id)
          if(abs(diff).gt.0.05)z1(id+1)=z1(id)+0.05*diff/abs(diff)
          if((el(kt)-z1(id+1)).gt.(el(2)-0.10))z1(id+1)=el(kt)-
     .         (el(2)-0.10)
          end if

********* Tridiagonal coefficients (A, C, D, and V)

          DO I=IU,ID
            A(I) = -RHO(KT,I-1)*G*DLT**2*BHRHO(I-1)*1.0/(DLX(I)
     .             +DLX(I-1))
            C(I) = -RHO(KT,I+1)*G*DLT**2*BHRHO(I)*1.0/(DLX(I)+DLX(I+1))
            V(I) = RHO(KT,I)*G*DLT**2*(BHRHO(I)*1.0/(DLX(I)+DLX(I+1))
     .             +BHRHO(I-1)*1.0/(DLX(I)+DLX(I-1)))+DLX(I)*B(KTT(I),I)
            D(I) = DLT*(D(I)+DLT*(F(I)-F(I-1)))+DLX(I)*B(KTT(I),I)*Z2(I)
          END DO
          IF (UP_HEAD(JB)) D(IU) = D(IU)-A(IU)*Z1(IU-1)
          IF (DN_HEAD(JB)) D(ID) = D(ID)-C(ID)*Z1(ID+1)

********* Thomas algorithm for water surface elevation solution

          BTA(IU) = V(IU)
          GMA(IU) = D(IU)/BTA(IU)
          DO I=IU+1,ID
            BTA(I) = V(I)-A(I)*C(I-1)/BTA(I-1)
            GMA(I) = (D(I)-A(I)*GMA(I-1))/BTA(I)
          END DO
          Z1(ID) = GMA(ID)
          DO K=1,ID-IU
            I     = ID-K
            Z1(I) = GMA(I)-C(I)*Z1(I+1)/BTA(I)
          END DO
          IF (UP_FLOW(JB)) Z1(IU-1) = Z1(IU)
          IF (DN_FLOW(JB)) Z1(ID+1) = Z1(ID)

********* Update surface layer and related variables

          DO I=IU-1,ID+1
            K         = KT-1
            KTT(I)    = KT
            HEIGHT(I) = 0.0
            DO WHILE (HEIGHT(I).LT.-Z1(I))
              HEIGHT(I) = HEIGHT(I)+H(K)
              KTT(I)    = MAX(K,2)
              K         = MAX(K-1,1)
            END DO
            HEIGHT(I) = MAX(HEIGHT(I)-H(K+1),0.0)
            RATIO     = 1.0
            IF (KTT(I).NE.SKTT(I)) RATIO = B(SKTT(I),I)/B(KTT(I),I)
            IF (KTT(I).EQ.KT.OR.SKTT(I).EQ.KT) THEN
              Z1(I) = Z1(I)*RATIO
            ELSE
              IF (KTT(I).GT.SKTT(I)) THEN
                Z1(I) = (HEIGHT(I)+Z1(I)+H(KTT(I)))*RATIO-(HEIGHT(I)
     .                  +H(KTT(I)))
              ELSE IF (KTT(I).LT.SKTT(I)) THEN
                Z1(I) = (HEIGHT(I)+Z1(I))*RATIO-HEIGHT(I)
              END IF
            END IF
            HKT1(I)  = H(KT)-Z1(I)
            AVHKT(I) = (HKT1(I)+H(KT+1))*0.5
            IF (KTT(I).EQ.KT) THEN
              BHKT1(I) = B(KT,I)*HKT1(I)
            ELSE
              BHKT1(I) = 0.0
              DO K=KT,KTT(I)+1,-1
                BHKT1(I) = BHKT1(I)+BH(K,I)
              END DO
              BHKT1(I) = BHKT1(I)-B(KTT(I),I)*(HEIGHT(I)+Z1(I))
            END IF
            BKT(I) = BHKT1(I)/HKT1(I)
          END DO
          DO I=IU-1,ID
            BHRKT1(I) = (BHKT1(I)+BHKT1(I+1))*0.5
          END DO
          BHRKT1(ID+1) = BHKT1(ID+1)
          DLVOL(JB)    = 0.0
          DO I=IU,ID
            DLVOL(JB) = DLVOL(JB)+(BHKT1(I)-BHKT2(I))*DLX(I)
          END DO

************************************************************************
**           Task 6.2.3: Compute Longitudinal Velocities              **
************************************************************************

          IUT = IU
          IDT = ID
          IF (UP_HEAD(JB)) IUT = IU-1
          IF (DN_HEAD(JB)) IDT = ID+1

********* Pressures

          DO I=IUT,IDT
            P(KT,I) = RHO(KT,I)*G*HKT1(I)
            DO K=KT+1,KB(I)
              P(K,I) = P(K-1,I)+RHO(K,I)*G*H(K)
            END DO
          END DO

********* Horizontal pressure gradients

          DO I=IUT,IDT-1
            DLXRHO    = 1.0/((DLX(I)+DLX(I+1))*RHOI)
            HPG(KT,I) = DLXRHO*(BKT(I)+BKT(I+1))*0.5*(HKT2(I+1)
     .                  *P(KT,I+1)-HKT2(I)*P(KT,I))
            DO K=KT+1,KBMIN(I)
              HPG(K,I) = DLXRHO*BHR(K,I)*((P(K-1,I+1)-P(K-1,I))
     .                   +(P(K,I+1)-P(K,I)))
            END DO
          END DO

********* Boundary horizontal velocities

          IF (UP_FLOW(JB)) THEN
            BHSUM       = BHKT1(IU)
            QINF(KT,JB) = 0.0
            DO K=KT+1,KB(IU)
              BHSUM      = BHSUM+BH(K,IU)
              QINF(K,JB) = 0.0
            END DO

*********** Place inflow according to density

            IF (PLACE_QIN) THEN
              KQIN  = KT
              ttt=TIN(JB)
              RHOIN = DENSITY (ttt,CIN(2,JB),CIN(4,JB))
              DO WHILE (RHOIN.GT.RHO(KQIN,IU).AND.KQIN.LT.KB(IU))
                KQIN = KQIN+1
              END DO
              TOTQIN = QIN(JB)*DLT
              VICELL = BH(KQIN,IU)*DLX(IU)
              IF (KQIN.EQ.KT) VICELL = BHKT1(IU)*DLX(IU)
              QINF(KQIN,JB) = 1.0
              IF (TOTQIN.GE.0.9*VICELL) THEN
                QINF(KQIN,JB) = 0.9
                IF (KQIN.EQ.KT) QINF(KT+1,JB) = 0.1
                IF (KQIN.NE.KT) QINF(KT-1,JB) = 0.1
              END IF

*********** Place inflows according to weighted cell volumes

            ELSE
              QINF(KT,JB) = BHKT1(IU)/BHSUM
              DO K=KT+1,KB(IU)
                QINF(K,JB) = BH(K,IU)/BHSUM
              END DO
            END IF
            U(KT,IU-1) = QINF(KT,JB)*QIN(JB)/BHRKT1(IU-1)
            DO K=KT+1,KB(IU)
              U(K,IU-1) = QINF(K,JB)*QIN(JB)/BHR(K,IU-1)
            END DO
          END IF
          IF (DN_FLOW(JB)) THEN
            DO JO=1,NOUT(JB)
              K = KOUT(JO,JB)
              IF (K.GE.KT) THEN
                BHRT = BHR(K,ID)
                IF (K.EQ.KT) BHRT = BHRKT1(ID)
                U(K,ID) = QOUT(JO,JB)/BHRT
              END IF
            END DO
          END IF
          IF (UP_HEAD(JB)) THEN
            U(KT,IU-1) = (BHRKT2(IU-1)*U(KT,IU-1)+DLT*(-SB(KT,IU-1)
     .                   +ST(KT,IU-1)-HPG(KT,IU-1)))/BHRKT1(IU-1)
            DO K=KT+1,KB(IU-1)
              U(K,IU-1) = (BHR(K,IU-1)*U(K,IU-1)+DLT*(-SB(K,IU-1)
     .                    +ST(K,IU-1)-HPG(K,IU-1)))/BHR(K,IU-1)
            END DO
          END IF
          IF (DN_HEAD(JB)) THEN
            U(KT,ID) = (BHRKT2(ID)*U(KT,ID)+DLT*(-SB(KT,ID)+ST(KT,ID)
     .                 -HPG(KT,ID)))/BHRKT1(ID)
            DO K=KT+1,KB(ID+1)
              U(K,ID) = (BHR(K,ID)*U(K,ID)+DLT*(-SB(K,ID)+ST(K,ID)
     .                  -HPG(K,ID)))/BHR(K,ID)
            END DO
          END IF

********* Horizontal velocities

          DO I=IU,ID-1
            U(KT,I) = (BHRKT2(I)*U(KT,I)+DLT*(-SB(KT,I)+ST(KT,I)
     .                -ADMZ(KT,I)+DM(KT,I)-ADMX(KT,I)-HPG(KT,I)))
     .                /BHRKT1(I)
            DO K=KT+1,KBMIN(I)
              U(K,I) = U(K,I)+DLT/BHR(K,I)*(-SB(K,I)+ST(K,I)-ADMZ(K,I)
     .                 +ADMZ(K-1,I)-ADMX(K,I)+DM(K,I)-HPG(K,I))
            END DO
          END DO

********* Corrected horizontal velocities



c set test2=ON for skipping flow correction
          if(test2)go to 1790
c compute qc from both directions and use the one with the lowest mean error
c          if(z2(iu).lt.z2(id))then

          Q(IU-1) = U(KT,IU-1)*BHRKT1(IU-1)
          DO K=KB(IU),KT+1,-1
            Q(IU-1) = Q(IU-1)+U(K,IU-1)*BHR(K,IU-1)
          END DO
          QC1(IU-1) = Q(IU-1)
          DO I=IU,ID
            IF (I.EQ.ID.AND.DN_FLOW(JB)) THEN
              Q(I)     = 0.0
              QSSUM(I) = 0.0
              do k=kt,kb(id)
                  BHRT = BHR(K,I)
                  IF (K.EQ.KT) BHRT = BHRKT1(I)
                  Q(I)     = Q(I)+U(K,I)*BHRT
                  QSSUM(I) = QSSUM(I)+QSS(K,I)
              end do
            ELSE
              Q(I)     = U(KT,I)*BHRKT1(I)
              QSSUM(I) = QSS(KT,I)
              DO K=KBMIN(I),KT+1,-1
                Q(I)     = Q(I)+U(K,I)*BHR(K,I)
               QSSUM(I) = QSSUM(I)+QSS(K,I)
              END DO
            END IF
            QC1(I) = QC1(I-1)+DLX(I)/DLT*(BHKT2(I)-BHKT1(I))+QSSUM(I)
          END DO

          QC2(ID) = Q(ID)
          iut=iu
          if(up_head(jb))iut=iu-1
          DO I=id-1,iut,-1
            QC2(I) = QC2(I+1)-DLX(I+1)/DLT*
     .            (BHKT2(I+1)-BHKT1(I+1))-QSSUM(I+1)

          END DO

c compute which qc has lowest error
           err1=0.0
           err2=0.0

           if(up_head(jb))then
           do i=iu-1,id
           qc(i)=qc2(i)
           end do
           go to 1678
           end if

           if(dn_head(jb))then
           do i=iu-1,id
           qc(i)=qc1(i)
           end do
           go to 1678
           end if

c for this model if qout.ne.0.0 it would apply to k=1 and 2
           if(qout(1,jb).eq.0.0)then
           do jw=1,nwp
           if(iwd(jw).eq.iu)go to 1679
           end do
           go to 1680
1679       if(qwd(jw).ne.0.0)then
           do i=iu-1,id
           qc(i)=qc2(i)
           end do
           go to 1678
           end if
           end if
1680       continue

           do i=iu,id
           if(q(i).ne.0.0)then
           err1=abs((q(i)-qc1(i))/q(i))+err1
           err2=abs((q(i)-qc2(i))/q(i))+err2
           end if
           end do

           if(err1.lt.err2)then
             do i=iu-1,id
             qc(i)=qc1(i)
             end do
           else
             do i=iu-1,id
             qc(i)=qc2(i)
             end do
          end if
1678        iut=iu
            if(up_head(jb))iut=iu-1
          do i=iut,id
              BHRSUM = 0.0
              BHRSUM   = BHRKT1(I)
              DO K=KBMIN(I),KT+1,-1
                BHRSUM   = BHRSUM+BHR(K,I)
              END DO
            IF (I.EQ.ID.AND.DN_FLOW(JB)) THEN
                  do k=kb(i),kt,-1
                  BHRT = BHR(K,I)
                  IF (K.EQ.KT) BHRT = BHRKT1(I)
                  IF (U(K,I).NE.0.0) BHRSUM = BHRSUM+BHRT
                  end do
              DO JO=1,NOUT(JB)
                K = KOUT(JO,JB)
                IF (K.GE.KT) THEN
                  IF (U(K,I).NE.0.0) U(K,I) = U(K,I)+(QC(I)-Q(I))
     .                                        /BHRSUM
                END IF
              END DO
            ELSE
              DO K=KBMIN(I),KT,-1
                U(K,I) = U(K,I)+(QC(I)-Q(I))/BHRSUM
              END DO
              end if
          end do

1790      continue

********* Head boundary flows

          IF (UP_HEAD(JB)) THEN
            QUH1(KT,JB) = U(KT,IU-1)*BHRKT1(IU-1)
            DO K=KT+1,KB(IU-1)
              QUH1(K,JB) = U(K,IU-1)*BHR(K,IU-1)
            END DO
          END IF
          IF (DN_HEAD(JB)) THEN
            QDH1(KT,JB) = U(KT,ID)*BHRKT1(ID)
            DO K=KT+1,KB(ID+1)
              QDH1(K,JB) = U(K,ID)*BHR(K,ID)
            END DO
          END IF

c              DO JO=1,NOUT(JB)
c                K = KOUT(JO,JB)
c                  QSS(K,ID)=QSS(K,ID)-QOUT(JO,JB)
c              END DO

************************************************************************
**            Task 6.2.4: Compute Vertical Velocities                 **
************************************************************************

          DO I=IU,ID
            DO K=KB(I)-1,KT,-1
              WT1    = W(K+1,I)*BB(K+1,I)
              WT2    = (BHR(K+1,I)*U(K+1,I)-BHR(K+1,I-1)*U(K+1,I-1)
     .                 -QSS(K+1,I))/DLX(I)
              W(K,I) = (WT1+WT2)/BB(K,I)
            END DO
            WT1       = W(KT,I)*BB(KT,I)
            WT2       = (BHRKT1(I)*U(KT,I)-BHRKT1(I-1)*U(KT,I-1)
     .                  -QSS(KT,I))/DLX(I)
            W(KT-1,I) = (WT1+WT2)/BB(KT-1,I)
          END DO
        END DO

************************************************************************
**               Task 6.2.5: Determine Next Timestep                  **
************************************************************************

        DO JB=1,NBP
          DO I=CUS(JB),DS(JB)
            DEPTH = HSEG(KT,I)-Z1(I)
            IF (DEPTH.LT.0.0) THEN
              IF (.NOT.WARNING_OPEN) THEN
                OPEN (WRN,FILE='wrn.opt',STATUS='UNKNOWN')
                WARNING_OPEN = .TRUE.
              END IF
              WRITE (WRN,6020) JDAY,I,Z1(I),HSEG(KT,I)-Z1(I)
              DEPTH = HSEG(KT,I)
            END IF
            TAU1   = 2.0*MAX(AX,IDX)/(DLX(I)*DLX(I))
            TAU2   = 2.0*AZ(KT,I)/(HKT1(I)*HKT1(I))
            CELRTY = SQRT((ABS(RHO(KB(I),I)-RHO(KT,I)))/1000.0*G
     .               *DEPTH*0.5)
            QTOT   = (ABS(U(KT,I))*BHKT1(I)+ABS(U(KT,I-1))*BHKT1(I-1)
     .               +ABS(W(KT,I))*BB(KT,I)*DLX(I))*0.5+BB(KT,I)
     .               *DLX(I)*ABS(HKT2(I)-HKT1(I))/DLT
            DLTCAL = 1.0/((QTOT/BHKT1(I)+CELRTY)/DLX(I)+TAU1+TAU2)
            IF (DLTCAL.LT.CURMAX) THEN
              KLOC   = KT
              ILOC   = I
              CURMAX = INT(DLTCAL)
              IF (CURMAX.LT.MINDLT) THEN
                KMIN = KT
                IMIN = I
              END IF
            END IF
            DO K=KT+1,KB(I)
              TAU2   = 2.0*AZ(K,I)/(H(K)*H(K))
              QTOT   = ABS(U(K,I))*BHR(K,I)+ABS(U(K,I-1))*BHR(K,I-1)
     .                 +(ABS(W(K,I))*BB(K,I)+ABS(W(K-1,I))*BB(K-1,I))
     .                 *DLX(I)*0.5
              DLTCAL = 1.0/((QTOT/BH(K,I)+CELRTY)/DLX(I)+TAU1+TAU2)
              IF (DLTCAL.LT.CURMAX) THEN
                KLOC   = K
                ILOC   = I
                CURMAX = INT(DLTCAL)
                IF (CURMAX.LT.MINDLT) THEN
                  KMIN = K
                  IMIN = I
                END IF
              END IF
            END DO
          END DO
        END DO

******* Restore timestep dependent variables and restart calculations

        VIOLATION = .FALSE.
        IF (CURMAX.LT.DLT.AND.DLT.GT.DLTMIN) THEN
          DLT = INT(MAX(DLTMIN,DLTF(DLTDP)*CURMAX))
          IF (DLT.LE.DLTMIN) THEN
            IF (.NOT.WARNING_OPEN) THEN
              OPEN (WRN,FILE='wrn.opt',STATUS='UNKNOWN')
              WARNING_OPEN = .TRUE.
            END IF
            WRITE (WRN,6000) JDAY,DLT
          END IF
          CURMAX = DLTMAX(DLTDP)/DLTF(DLTDP)
          IF (DLT.LT.MINDLT) THEN
            MINDLT = DLT
            JDMIN  = JDAY+DLT/86400.0
          END IF
          DO JB=1,NBP
            DO I=CUS(JB)-1,DS(JB)+1
              Z1(I)  = SZ1(I)
              KTT(I) = SKTT(I)
              DO K=KT,KB(I)
                U(K,I)   = SU(K,I)
                W(K,I)   = SW(K,I)
                QSS(K,I) = 0.0
              END DO
            END DO
          END DO
          NV        = NV+1
          VIOLATION = .FALSE.
          GO TO 10010
        END IF

************************************************************************
**     Task 6.3: Compute Temporal Balance Terms and Temperatures      **
************************************************************************

c set bank shading back to almost no shading = = 0.97 when winter comes,
c define loss of deciduous trees in FALL as occurring in 1992 after JD 294
c ==October 20
c
        if(iflag2.eq.0)then
        if(jday.le.91.0.or.jday.ge.294.0)then
        iflag2=1
          do i=1,imp
          if(banksh(i).lt.0.97)banksh(i)=0.97
          end do
        end if
        end if

        IF (.NOT.NO_HEAT) THEN
          IF (TERM_BY_TERM) CALL HEAT_EXCHANGE (JDAY)
          RS = SRO*RHOI*CP
          IF (TERM_BY_TERM) CALL RADIATION (RS,CLOUD,TAIR,RSN,RAN)
        END IF
        DO JB=1,NBP
          IU = CUS(JB)
          ID = DS(JB)

********* Surface heat exchange

          IF (.NOT.NO_HEAT) THEN
            DO I=IU,ID
              IF (.NOT.ICE(I)) THEN
                GAMMA     = EXH2O+EXINOR*SS(KT,I)+EXORG*(ALGAE(KT,I)
     .                      +DETRIT(KT,I)+zoop(kt,i))
                DEPTH     = HKT1(I)
                TSS(KT,I) = TSS(KT,I)-(1.0-BETA)*SRO*EXP(-GAMMA*DEPTH)
     .                      *B(KT,I)*DLX(I)
                IF (TERM_BY_TERM) THEN
                  CALL SURFACE_TERMS (TAIR,TDEW,WIND,T2(KT,I),RB,RE,RC)
                  RN(I)     = (RSN*banksh(i)+RAN-RB-RE-RC)
                  TSS(KT,I) = TSS(KT,I)+RN(I)/RHOI/CP*B(KTT(I),I)*DLX(I)
                ELSE
                  TSS(KT,I) = TSS(KT,I)+CSHE*(ET-T2(KT,I))*B(KTT(I),I)
     .                        *DLX(I)
                END IF
                DO K=KT+1,KB(I)
                  GAMMA    = EXH2O+EXINOR*SS(K,I)+EXORG*(ALGAE(K,I)
     .                       +DETRIT(K,I)+zoop(kt,i))
                  TSS(K,I) = TSS(K,I)+(1.0-BETA)*SRO*(EXP(-GAMMA*DEPTH)
     .                       -EXP(-GAMMA*(DEPTH+H(K))))*B(K,I)*DLX(I)
                  DEPTH    = DEPTH+H(K)
                END DO
              END IF

************* Bottom                                                    !10/9/93
c  added to upper slough model 9/15/95, chris
                                                                        !10/9/93
              TSS(KT,I) = TSS(KT,I)+CBHE*(TSED-T2(KT,I))*(B(KT,I)       !10/9/93
     .                    -B(KT+1,I))*DLX(I)                            !10/9/93
              DO K=KT+1,KB(I)-1                                         !10/9/93
                TSS(K,I) = TSS(K,I)+CBHE*(TSED-T2(K,I))*(B(K,I)         !10/9/93
     .                     -B(K+1,I))*DLX(I)                            !10/9/93
              END DO                                                    !10/9/93
              TSS(KB(I),I) = TSS(KB(I),I)+CBHE*(TSED-T2(KB(I),I))       !10/9/93
     .                       *(B(KB(I),I))*DLX(I)                       !10/9/93

            END DO

            IF (ICE_CALC) THEN
              HIA = 0.2367*CSHE/5.65E-8
              DO I=IU,ID
                ALLOW_ICE(I) = .TRUE.
                DO K=KT,KB(I)
                  IF (T2(K,I).GT.ICET2) THEN
                    ALLOW_ICE(I) = .FALSE.
                    GO TO 10020
                  END IF
                END DO
              END DO
10020         CONTINUE
              ICE_IN(JB) = .TRUE.
              DO I=IU,ID
                IF (ICETH(I).LT.ICEMIN) THEN
                  ICE_IN(JB) = .FALSE.
                  GO TO 10030
                END IF
              END DO
10030         CONTINUE
              DO I=IU,ID
                IF (DETAILED_ICE) THEN
                  IF (T2(KT,I).LT.0.0) THEN
                    IF (INIT(I).EQ.0) THEN
                      ICETH2 = -T2(KT,I)*RHO(KT,I)*CP*HKT2(I)/(RHOICE
     .                         *RL1)
                      IF (ICETH2.LT.0.005) THEN
                        ICETH2  = 0.0
                        INIT(I) = 0
                      ELSE
                        TSS(KT,I) = TSS(KT,I)-T2(KT,I)*RHO(KT,I)*CP
     .                              *HKT2(I)*B(KTT(I),I)/4.186E6/DLT
     .                              *DLX(I)
                        INIT(I)   = 1
                      END IF
                    END IF
                  END IF

***************** Ice present

                  IF (ICE(I)) THEN
                    ICOUNT = 0
                    TICE   = TAIR
10040               CONTINUE
                    ICOUNT = ICOUNT+1
                    IF (ICOUNT.GT.2000) GO TO 10050
                    CALL SURFACE_TERMS (TAIR,TDEW,WIND,TICE,RB,RE,RC)
                    RN(I) = RS*(1.0-ALBEDO)*BETAI+RAN-RB-RE-RC

******************* Heat balance at air-ice interface

                    DEL = RN(I)+RK1*(RIMT-TICE)/ICETH(I)
                    IF (ABS(DEL).GT.1.0) THEN
                      TICE = TICE+DEL/500.0
                      GO TO 10040
                    END IF
10050               CONTINUE

******************* Solar radiation attenuation

                    SRAEI     = SRO*(1.0-ALBEDO)*(1.0-BETAI)*EXP(-GAMMAI
     .                          *ICETH(I))*B(KTT(I),I)
                    TSS(KT,I) = TSS(KT,I)+SRAEI*DLX(I)
                    IF (TICE.GT.0.0) THEN
                      HICE   = RHOICE*CP*0.5*TICE/2.0*ICETH(I)
     .                         *B(KTT(I),I)/4.186E6/DLT
                      ICETHU = -DLT*HICE/B(KT,I)*4.186E6/(RHOICE*RL1)
                      TICE   = 0.0
                    END IF

******************* Ice growth at ice-water interface

                    IF (TICE.LT.0.0) THEN
                      ICETH1 = DLT*(RK1*(RIMT-TICE)/ICETH(I))/(RHOICE
     .                         *RL1)
                    END IF

******************* Ice melt at ice-water interface

                    IF (T2(KT,I).GT.0.0) THEN
                      ICETH2    = -DLT*HWI*(T2(KT,I)-RIMT)/(RHOICE*RL1)
                      TSS(KT,I) = TSS(KT,I)+2.392E-7*HWI*(RIMT-T2(KT,I))
     .                            *B(KTT(I),I)*DLX(I)
                    END IF
                  END IF

***************** Total ice thickness

                  ICETH(I) = ICETH(I)+ICETHU+ICETH1+ICETH2
                  IF (ICETH(I).LT.0.005) ICETH(I) = 0.0
                  IF (WINTER.AND.(.NOT.ICE_IN(JB))) THEN
                    IF (.NOT.ALLOW_ICE(I)) ICETH(I) = 0.0
                  END IF
                  ICE(I)   = ICETH(I).GT.0.0
                  ICESW(I) = 1.0
                  INIT(I)  = 0
                  IF (ICE(I)) THEN
                    ICESW(I) = 0.0
                    INIT(I)  = 1
                  END IF
                  ICETHU = 0.0
                  ICETH1 = 0.0
                  ICETH2 = 0.0
                  IF (ICETH(I).LT.0.005.AND.ICETH(I).GT.0.0) THEN
                    ICETH(I) = 0.005
                  END IF
                ELSE
                  HIA      = 0.2367*CSHE/5.65E-8
                  ICETH(I) = ICETH(I)+DLT*((RIMT-ET)/(ICETH(I)/RK1+1.0
     .                       /HIA)-(T2(KT,I)-RIMT))/(RHOICE*RL1)
                  ICETH(I) = MAX(ICETH(I),0.0)
                  ICE(I)   = ICETH(I).GT.0.0
                  IF (ICE(I)) THEN
                    TSS(KT,I) = TSS(KT,I)+2.392E-7*(RIMT-T2(KT,I))
     .                          *B(KTT(I),I)*DLX(I)
                  END IF
                END IF
              END DO
            END IF
          END IF

********* External heat sources/sinks and total inflow/outflow

          IF (EVAPORATION) THEN
            DO I=IU,ID
              VOLEV(JB) = VOLEV(JB)-EV(I)*DLT
              TSS(KT,I) = TSS(KT,I)-EV(I)*T2(KT,I)
            END DO
          END IF
          IF (PRECIPITATION) THEN
            DO I=IU,ID
              TSS(KT,I) = TSS(KT,I)+TPR(JB)*QPR(I)
              VOLPR(JB) = VOLPR(JB)+QPR(I)*DLT
            END DO
          END IF
          IF (TRIBUTARIES) THEN
            DO JT=1,NTP
              IF (JB.EQ.JBTR(JT)) THEN
                I = ITR(JT)
                IF (I.LT.CUS(JB)) I = CUS(JB)
                IF (PLACE_QTR) THEN
                   TSS(KTR(JT),I)=TSS(KTR(JT),I)+TTR(JT)*QTR(JT)
                ELSE
                DO K=KT,KB(I)
                  TSS(K,I) = TSS(K,I)+TTR(JT)*QTR(JT)*QTRF(K,JT)
                END DO
                END IF
                VOLTR(JB) = VOLTR(JB)+QTR(JT)*DLT
              END IF
            END DO
          ENDIF
          IF (DIST_TRIBS(JB)) THEN
            DO I=IU,ID
              TSS(KT,I) = TSS(KT,I)+TDTR(JB)*QDT(I)
              VOLDT(JB) = VOLDT(JB)+QDT(I)*DLT
            END DO
          END IF
          IF (WITHDRAWALS) THEN
          if(kt.le.5)qf=4.
          if(kt.eq.6)qf=3.
          if(kt.eq.7)qf=2.
            DO JW=1,NWP
              IF (JB.EQ.JBWD(JW)) THEN
                IF (KWD(JW).GE.KT) THEN
                  TSS(KWD(JW),IWD(JW))=TSS(KWD(JW),IWD(JW))-
     .             T2(KWD(JW),IWD(JW))*QWD(JW)/qf
c note that the kwd for all these is cell 8, below this outflow
c is spread among cells 5, 6, and 7 also
                  if(5.ge.kt)TSS(5,IWD(JW))=TSS(5,IWD(JW))-
     .             T2(5,IWD(JW))*QWD(JW)/qf
                  if(6.ge.kt)TSS(6,IWD(JW))=TSS(6,IWD(JW))-
     .             T2(6,IWD(JW))*QWD(JW)/qf
                  if(7.ge.kt)TSS(7,IWD(JW))=TSS(7,IWD(JW))-
     .             T2(7,IWD(JW))*QWD(JW)/qf
                  VOLWD(JB) = VOLWD(JB)-QWD(JW)*DLT
                END IF
              END IF
            END DO
          ENDIF
          IF (UP_FLOW(JB)) THEN
            VOLIN(JB) = VOLIN(JB)+QIN(JB)*DLT
            DO K=KT,KB(IU)
              TSS(K,IU) = TSS(K,IU)+QINF(K,JB)*QIN(JB)*TIN(JB)
            END DO
          END IF
          IF (DN_FLOW(JB)) THEN
            DO JO=1,NOUT(JB)
              K = KOUT(JO,JB)
              IF (K.GE.KT) THEN
c                TSS(K,ID)  = TSS(K,ID)-QOUT(JO,JB)*T2(K,ID)
         if(k.eq.kt)tss(k,id)=tss(k,id)-u(k,id)*bhrkt2(id)*t2(k,id)
         if(k.ne.kt)TSS(K,ID)=TSS(K,ID)-u(k,id)*bhr(k,id)*t2(k,id)
                VOLOUT(JB) = VOLOUT(JB)-QOUT(JO,JB)*DLT
              END IF
            END DO
          END IF
          IF (UP_HEAD(JB)) THEN
            IUT = IU
            IF (QUH1(KT,JB).GE.0.0) IUT = IU-1
            TSSUH1(KT,JB) = T2(KT,IUT)*QUH1(KT,JB)
            TSS(KT,IU)    = TSS(KT,IU)+TSSUH1(KT,JB)
            VOLUH(JB)     = VOLUH(JB)+QUH1(KT,JB)*DLT
            DO K=KT+1,KB(IU)
              IUT = IU
              IF (QUH1(K,JB).GE.0.0) IUT = IU-1
              TSSUH1(K,JB) = T2(K,IUT)*QUH1(K,JB)
              TSS(K,IU)    = TSS(K,IU)+TSSUH1(K,JB)
              VOLUH(JB)    = VOLUH(JB)+QUH1(K,JB)*DLT
            END DO
          END IF
          IF (UH_INTERNAL(JB)) THEN
            DO K=KT,KB(IU-1)
              TSS(K,UHS(JB))  = TSS(K,UHS(JB))-TSSUH2(K,JB)/DLT
              VOLUH(JBUH(JB)) = VOLUH(JBUH(JB))-QUH2(K,JB)
            END DO
          END IF
          IF (DN_HEAD(JB)) THEN
            IUT = ID+1
            IF (QDH1(KT,JB).GE.0.0) IUT = ID
            TSSDH1(KT,JB) = T2(KT,IUT)*QDH1(KT,JB)
            TSS(KT,ID)    = TSS(KT,ID)-TSSDH1(KT,JB)
            VOLDH(JB)     = VOLDH(JB)-QDH1(KT,JB)*DLT
            DO K=KT+1,KB(ID+1)
              IUT = ID+1
              IF (QDH1(K,JB).GE.0.0) IUT = ID
              TSSDH1(K,JB) = T2(K,IUT)*QDH1(K,JB)
              TSS(K,ID)    = TSS(K,ID)-TSSDH1(K,JB)
              VOLDH(JB)    = VOLDH(JB)-QDH1(K,JB)*DLT
            END DO
          END IF
          IF (DH_INTERNAL(JB)) THEN
            DO K=KT,KB(ID+1)
              TSS(K,DHS(JB))  = TSS(K,DHS(JB))+TSSDH2(K,JB)/DLT
              VOLDH(JBDH(JB)) = VOLDH(JBDH(JB))+QDH2(K,JB)
            END DO
          END IF

********* Horizontal advection and diffusion multipliers

          DO I=IU,ID-1
            DO K=KT,KB(I)
              COUR = U(K,I)*DLT/DLX(I)
              IF (U(K,I).GE.0.0) THEN
                T1L = T2(K,I-1)
                T2L = T2(K,I)
                T3L = T2(K,I+1)
                IF (K.GT.KB(I-1).OR.I.EQ.IU) T1L = T2(K,I)
                IF (UPWIND) THEN
                  DX1(K,I)  = 0.0
                  DX2(K,I)  = -DX(K,I)/SF2L(I)
                  DX3(K,I)  = DX(K,I)/SF2L(I)
                  AD1L(K,I) = 0.0
                  AD2L(K,I) = 1.0
                  AD3L(K,I) = 0.0
                ELSE
                  DX1(K,I)  = DX(K,I)*SF12L(I,1)
                  DX2(K,I)  = DX(K,I)*SF13L(I,1)
                  DX3(K,I)  = DX(K,I)*SF14L(I,1)
                  ALFA      = 2.0*(DX(K,I)*DLT/(SF2L(I)*SF2L(I))
     .                        -(1.0-COUR*COUR)/6.0)*SF4L(I,1)
                  AD1L(K,I) = (ALFA-COUR*SF9L(I,1)/2.0)/SF6L(I,1)
                  AD2L(K,I) = SF5L(I,1)+(ALFA-COUR*SF10L(I,1)/2.0)
     .                        /SF7L(I,1)
                  AD3L(K,I) = SF3L(I,1)+(ALFA-COUR*SF11L(I,1)/2.0)
     .                        /SF8L(I,1)
                END IF
              ELSE
                T1L = T2(K,I)
                T2L = T2(K,I+1)
                T3L = T2(K,I+2)
                IF (K.GT.KB(I+2).OR.I.EQ.ID-1) T3L = T2(K,I+1)
                IF (UPWIND) THEN
                  DX1(K,I)  = -DX(K,I)/SF2L(I)
                  DX2(K,I)  = DX(K,I)/SF2L(I)
                  DX3(K,I)  = 0.0
                  AD1L(K,I) = 0.0
                  AD2L(K,I) = 1.0
                  AD3L(K,I) = 0.0
                ELSE
                  DX1(K,I)  = DX(K,I)*SF12L(I,2)
                  DX2(K,I)  = DX(K,I)*SF13L(I,2)
                  DX3(K,I)  = DX(K,I)*SF14L(I,2)
                  ALFA      = 2.0*(DX(K,I)*DLT/(SF2L(I)*SF2L(I))
     .                        -(1.0-COUR*COUR)/6.0)*SF4L(I,2)
                  AD1L(K,I) = SF3L(I,2)+(ALFA-COUR*SF9L(I,2)/2.0)
     .                        /SF6L(I,2)
                  AD2L(K,I) = SF5L(I,2)+(ALFA-COUR*SF10L(I,2)/2.0)
     .                        /SF7L(I,2)
                  AD3L(K,I) = (ALFA-COUR*SF11L(I,2)/2.0)/SF8L(I,2)
                END IF
              END IF
              TADL(K,I) =  (DX1(K,I)-U(K,I)*AD1L(K,I))*T1L
     .                    +(DX2(K,I)-U(K,I)*AD2L(K,I))*T2L
     .                    +(DX3(K,I)-U(K,I)*AD3L(K,I))*T3L
            END DO
          END DO

********* Vertical advection multipliers

          DO I=IU,ID
            DO K=KT,KB(I)-1
              IF (W(K,I).GE.0.0) THEN
                T1V = T2(K-1,I)
                T2V = T2(K,I)
                T3V = T2(K+1,I)
                IF (K.LE.KT+1) THEN
                  T1V  = T2(KT,I)
                  HTOP = HKT2(I)
                  HBOT = H(K+1)
                  HMID = H(K)
                  IF (K.EQ.KT) HMID = HKT2(I)
                  HMIN       = MIN(HBOT,HMID)
                  SF2V(K)    = (HBOT+HMID)/2.0
                  SF3V(K,1)  = HMID**2
                  SF4V(K,1)  = HMID/(HMID+HBOT)
                  SF5V(K,1)  = HBOT/(HMID+HBOT)
                  SF6V(K,1)  = 0.25*(HTOP+2.0*HMID+HBOT)*(HTOP+HMID)
                  SF7V(K,1)  = -0.25*(HMID+HBOT)*(HTOP+HMID)
                  SF8V(K,1)  = 0.25*(HMID+HBOT)*(HTOP+2.0*HMID+HBOT)
                  SF9V(K,1)  = 0.5*(HMID-HBOT)*HMIN
                  SF10V(K,1) = 0.5*(HTOP+2.0*HMID-HBOT)*HMIN
                  SF11V(K,1) = 0.5*(HTOP+3.0*HMID)*HMIN
                END IF
                IF (UPWIND) THEN
                  AD1V(K,I) = 0.0
                  AD2V(K,I) = 1.0
                  AD3V(K,I) = 0.0
                ELSE
                  COUR      = W(K,I)*DLT/SF2V(K)
                  ALFA      = 2.0*(DZ(K,I)*DLT/(SF2V(K)*SF2V(K))
     .                        -(1.0-COUR*COUR)/6.0)*SF3V(K,1)
                  AD1V(K,I) = (ALFA-COUR*SF9V(K,1)/2.0)/SF6V(K,1)
                  AD2V(K,I) = SF5V(K,1)+(ALFA-COUR*SF10V(K,1)/2.0)
     .                        /SF7V(K,1)
                  AD3V(K,I) = SF4V(K,1)+(ALFA-COUR*SF11V(K,1)/2.0)
     .                        /SF8V(K,1)
                END IF
              ELSE
                T1V = T2(K,I)
                T2V = T2(K+1,I)
                T3V = T2(K+2,I)
                IF (K.EQ.KB(I)-1) T3V = T2(K+1,I)
                IF (K.EQ.KT) THEN
                  HTOP       = HKT2(I)
                  HMID       = H(KT+1)
                  HBOT       = H(KT+2)
                  HMIN       = MIN(HTOP,HMID)
                  SF2V(K)    = (HMID+HTOP)/2.0
                  SF3V(K,2)  = HMID**2
                  SF4V(K,2)  = HMID/(HTOP+HMID)
                  SF5V(K,2)  = HTOP/(HTOP+HMID)
                  SF6V(K,2)  = 0.25*(HTOP+2.0*HMID+HBOT)*(HTOP+HMID)
                  SF7V(K,2)  = -0.25*(HMID+HBOT)*(HTOP+HMID)
                  SF8V(K,2)  = 0.25*(HTOP+2.0*HMID+HBOT)*(HMID+HBOT)
                  SF9V(K,2)  = -0.5*(3.0*HMID+HBOT)*HMIN
                  SF10V(K,2) = 0.5*(HTOP-2.0*HMID-HBOT)*HMIN
                  SF11V(K,2) = 0.5*(HTOP-HMID)*HMIN
                END IF
                IF (UPWIND) THEN
                  AD1V(K,I) = 0.0
                  AD2V(K,I) = 1.0
                  AD3V(K,I) = 0.0
                ELSE
                  COUR      = W(K,I)*DLT/SF2V(K)
                  ALFA      = 2.0*(DZ(K,I)*DLT/(SF2V(K)*SF2V(K))
     .                        -(1.0-COUR*COUR)/6.0)*SF3V(K,2)
                  AD1V(K,I) = SF4V(K,2)+(ALFA-COUR*SF9V(K,2)/2.0)
     .                        /SF6V(K,2)
                  AD2V(K,I) = SF5V(K,2)+(ALFA-COUR*SF10V(K,2)/2.0)
     .                        /SF7V(K,2)
                  AD3V(K,I) = (ALFA-COUR*SF11V(K,2)/2.0)/SF8V(K,2)
                END IF
              END IF
              TADV(K,I) = -W(K,I)*(AD1V(K,I)*T1V+AD2V(K,I)*T2V
     .                    +AD3V(K,I)*T3V)
            END DO
          END DO
        END DO

******* Transport heat

        DO JB=1,NBP
          IU = CUS(JB)
          ID = DS(JB)
          DO I=IU,ID
            T1(KT,I) = (T2(KT,I)*BHKT2(I)/DLT+(TADL(KT,I)*BHRKT1(I)
     .                 -TADL(KT,I-1)*BHRKT1(I-1))/DLX(I)+(1.0-THETA)
     .                 *TADV(KT,I)*BB(KT,I)+TSS(KT,I)/DLX(I))*DLT
     .                 /BHKT1(I)
            DO K=KT+1,KB(I)
              T1(K,I) = (T2(K,I)*BH(K,I)/DLT+(TADL(K,I)*BHR(K,I)
     .                  -TADL(K,I-1)*BHR(K,I-1))/DLX(I)+(1.0-THETA)
     .                  *(TADV(K,I)*BB(K,I)-TADV(K-1,I)*BB(K-1,I))
     .                  +TSS(K,I)/DLX(I))*DLT/BH(K,I)
            END DO

*********** Vertical advection and implicit diffusion

            K       = KT
            AT(K,I) = 0.0
            CT(K,I) = DLT/BHKT1(I)*(BB(K,I)*(THETA*0.5*W(K,I)-DZ(K,I)
     .                /AVHKT(I)))
            VT(K)   = 1.0+DLT/BHKT1(I)*(BB(K,I)*(DZ(K,I)/AVHKT(I)+THETA
     .                *0.5*W(K,I)))
            DT(K)   = T1(K,I)
            K       = KT+1
            AT(K,I) = -DLT/BH(K,I)*(BB(K-1,I)*(DZ(K-1,I)/AVHKT(I)+THETA
     .                *0.5*W(K-1,I)))
            CT(K,I) = DLT/BH(K,I)*(BB(K,I)*(THETA*0.5*W(K,I)-DZ(K,I)
     .                /AVH(K)))
            VT(K)   = 1.0+DLT/BH(K,I)*(BB(K,I)*(DZ(K,I)/AVH(K)+THETA
     .                *0.5*W(K,I))+BB(K-1,I)*(DZ(K-1,I)/AVHKT(I)
     .                -THETA*0.5*W(K-1,I)))
            DT(K)   = T1(K,I)
            DO K=KT+2,KB(I)-1
              AT(K,I) = -DLT/BH(K,I)*(BB(K-1,I)*(DZ(K-1,I)/AVH(K-1)
     .                  +THETA*0.5*W(K-1,I)))
              CT(K,I) = DLT/BH(K,I)*(BB(K,I)*(THETA*0.5*W(K,I)-DZ(K,I)
     .                  /AVH(K)))
              VT(K)   = 1.0+DLT/BH(K,I)*(BB(K,I)*(DZ(K,I)/AVH(K)+THETA
     .                  *0.5*W(K,I))+BB(K-1,I)*(DZ(K-1,I)/AVH(K-1)
     .                  -THETA*0.5*W(K-1,I)))
              DT(K)   = T1(K,I)
            END DO
            K = KB(I)
            IF (KB(I)-KT.GT.1) THEN
              AT(K,I) = -DLT/BH(K,I)*(BB(K-1,I)*(DZ(K-1,I)/AVH(K-1)
     .                  +THETA*0.5*W(K-1,I)))
              CT(K,I) = 0.0
              VT(K)   = 1.0+DLT/BH(K,I)*(BB(K-1,I)*(DZ(K-1,I)/AVH(K-1)
     .                  -THETA*0.5*W(K-1,I)))
              DT(K)   = T1(K,I)
            ELSE
              AT(K,I) = -DLT/BH(K,I)*(BB(K-1,I)*(DZ(K-1,I)/AVHKT(I)
     .                  +THETA*0.5*W(K-1,I)))
              CT(K,I) = 0.0
              VT(K)   = 1.0+DLT/BH(K,I)*(BB(K-1,I)*(DZ(K-1,I)/AVHKT(I)
     .                  -THETA*0.5*W(K-1,I)))
              DT(K)   = T1(K,I)
            END IF

*********** Tridiagonal solution

            BTAT(KT,I) = VT(KT)
            DO K=KT+1,KB(I)
              BTAT(K,I) = VT(K)-AT(K,I)/BTAT(K-1,I)*CT(K-1,I)
            END DO
            GMAT(KT) = DT(KT)
            DO K=KT+1,KB(I)
              GMAT(K) = DT(K)-AT(K,I)/BTAT(K-1,I)*GMAT(K-1)
            END DO
            T1(KB(I),I) = GMAT(KB(I))/BTAT(KB(I),I)
            DO K=KB(I)-1,KT,-1
              T1(K,I) = (GMAT(K)-CT(K,I)*T1(K+1,I))/BTAT(K,I)
            END DO
          END DO
        END DO

************************************************************************
**                 Task 6.4:  Update Constituents                     **
************************************************************************

        IF (CONSTITUENTS) THEN
          PALT = (1.-((EL(KT)-Z2(DS(1)))/1000.0)/44.3)**5.25
          DO JB=1,NBP
            IU = CUS(JB)
            ID = DS(JB)

*********** Internal sources/sinks

            IF (UPDATE_KINETICS) THEN
              CALL RATE_MULTIPLIERS
              CALL DECAY_CONSTANTS
              DO JAC=1,NAC
                JC = CN(JAC)
                IF (JC.EQ.2)  CALL SUSPENDED_SOLIDS
                IF (JC.EQ.3)  CALL COLIFORM
                IF (JC.EQ.5)  CALL LABILE_DOM
                IF (JC.EQ.6)  CALL REFRACTORY_DOM
                IF (JC.EQ.7)  CALL PHYTOPLANKTON
                IF (JC.EQ.8)  CALL DETRITUS
                IF (JC.EQ.9)  CALL PHOSPHOROUS
                IF (JC.EQ.10) CALL AMMONIA
                IF (JC.EQ.11) CALL NITRATE
                IF (JC.EQ.12) CALL DISSOLVED_OXYGEN
                IF (JC.EQ.13) CALL SEDIMENT
                IF (JC.EQ.14) CALL INORGANIC_CARBON
                IF (JC.EQ.16) CALL PH_CO2
                IF (JC.EQ.20) CALL IRON
                if(jc.eq.21)call biochemical_o2_demand
                if(jc.eq.22)call zooplankton
              END DO
            END IF
            DO JAC=1,NAC
              JC = CN(JAC)

************* External sources/sinks

              IF (TRIBUTARIES) THEN
                DO JT=1,NTP
                  IF (JB.EQ.JBTR(JT)) THEN
                    I = ITR(JT)
                    IF(I.LT.CUS(JB)) I = CUS(JB)
                    IF(PLACE_QTR) THEN
                    CSSB(KTR(JT),I,JC) = CSSB(KTR(JT),I,JC)+CTR(JC,JT)
     .                     *QTR(JT)
                    ELSE
                    DO K=KT,KB(I)
                      CSSB(K,I,JC) = CSSB(K,I,JC)+CTR(JC,JT)*QTR(JT)
     .                               *QTRF(K,JT)
                    END DO
                    END IF
                  END IF
                END DO
              END IF
              IF (DIST_TRIBS(JB)) THEN
                DO I=IU,ID
                  CSSB(KT,I,JC) = CSSB(KT,I,JC)+CDTR(JC,JB)*QDT(I)
                END DO
              END IF
          IF (WITHDRAWALS) THEN
          if(kt.le.5)qf=4.
          if(kt.eq.6)qf=3.
          if(kt.eq.7)qf=2.
            DO JW=1,NWP
              IF (JB.EQ.JBWD(JW)) THEN
                IF (KWD(JW).GE.KT) THEN
                  CSSB(KWD(JW),IWD(JW),JC)=CSSB(KWD(JW),IWD(JW),JC)
     .             -C2(KWD(JW),IWD(JW),JC)*QWD(JW)/qf
c note that the kwd for all these is cell 8, below this outflow
c is spread among cells 5, 6, and 7 also
                  if(5.ge.kt)CSSB(5,IWD(JW),JC)=CSSB(5,IWD(JW),JC)
     .             -C2(5,IWD(JW),JC)*QWD(JW)/qf
                  if(6.ge.kt)CSSB(6,IWD(JW),JC)=CSSB(6,IWD(JW),JC)
     .             -C2(6,IWD(JW),JC)*QWD(JW)/qf
                  if(7.ge.kt)CSSB(7,IWD(JW),JC)=CSSB(7,IWD(JW),JC)
     .             -C2(7,IWD(JW),JC)*QWD(JW)/qf
                END IF
              END IF
            END DO
          ENDIF
              IF (PRECIPITATION) THEN
                DO I=IU,ID
                  CSSB(KT,I,JC) = CSSB(KT,I,JC)+CPR(JC,JB)*QPR(I)
                END DO
              END IF
              IF (UP_FLOW(JB)) THEN
                DO K=KT,KB(IU)
                  CSSB(K,IU,JC) = CSSB(K,IU,JC)+QINF(K,JB)*QIN(JB)
     .                            *CIN(JC,JB)
                END DO
              END IF
              IF (DN_FLOW(JB)) THEN
                DO JO=1,NOUT(JB)
                  K = KOUT(JO,JB)
                  IF (K.GE.KT) THEN
       if(k.eq.kt)cssb(k,id,jc)=cssb(k,id,jc)
     .                         -u(k,id)*bhrkt2(id)*c2(k,id,jc)
       if(k.ne.kt)cssb(K,ID,jc)=cssb(K,ID,jc)
     .                         -u(k,id)*bhr(k,id)*c2(k,id,jc)
c                    CSSB(K,ID,JC) = CSSB(K,ID,JC)-QOUT(JO,JB)
c     .                              *C2(K,ID,JC)
                  END IF
                END DO
              END IF
              IF (UP_HEAD(JB)) THEN
                IUT = IU
                IF (QUH1(KT,JB).GE.0.0) IUT = IU-1
                CSSUH1(KT,JC,JB) = C2(KT,IUT,JC)*QUH1(KT,JB)
                CSSB(KT,IU,JC)   = CSSB(KT,IU,JC)+CSSUH1(KT,JC,JB)
                DO K=KT+1,KB(IU)
                  IUT = IU
                  IF (QUH1(K,JB).GE.0.0) IUT = IU-1
                  CSSUH1(K,JC,JB) = C2(K,IUT,JC)*QUH1(K,JB)
                  CSSB(K,IU,JC)   = CSSB(K,IU,JC)+CSSUH1(K,JC,JB)
                END DO
                IF (UH_INTERNAL(JB)) THEN
                  I = UHS(JB)
                  DO K=KT,KB(IU)
                    CSSB(K,I,JC) = CSSB(K,I,JC)-CSSUH2(K,JC,JB)/DLT
                  END DO
                END IF
              END IF
              IF (DN_HEAD(JB)) THEN
                IUT = ID+1
                IF (QDH1(KT,JB).GE.0.0) IUT = ID
                CSSDH1(KT,JC,JB) = C2(KT,IUT,JC)*QDH1(KT,JB)
                CSSB(KT,ID,JC)   = CSSB(KT,ID,JC)-CSSDH1(KT,JC,JB)
                DO K=KT+1,KB(ID+1)
                  IUT = ID+1
                  IF (QDH1(K,JB).GE.0.0) IUT = ID
                  CSSDH1(K,JC,JB) = C2(K,IUT,JC)*QDH1(K,JB)
                  CSSB(K,ID,JC)   = CSSB(K,ID,JC)-CSSDH1(K,JC,JB)
                END DO
                IF (DH_INTERNAL(JB)) THEN
                  I = DHS(JB)
                  DO K=KT,KB(ID+1)
                    CSSB(K,I,JC) = CSSB(K,I,JC)+CSSDH2(K,JC,JB)/DLT
                  END DO
                END IF
              END IF
            END DO
          END DO

********* Transport constituents

          DO JB=1,NBP
            IU = CUS(JB)
            ID = DS(JB)
            DO JAC=1,NAC
              JC = CN(JAC)
              IF (TRANSPORT(JC)) THEN

*************** Horizontal advection and diffusion terms

                DO I=IU,ID-1
                  DO K=KT,KB(I)
                    IF (U(K,I).GE.0.0) THEN
                      C1L = C1S(K,I-1,JC)
                      C2L = C1S(K,I,JC)
                      C3L = C1S(K,I+1,JC)
                      IF (K.GT.KB(I-1).OR.I.EQ.IU) C1L = C1S(K,I,JC)
                    ELSE
                      C1L = C1S(K,I,JC)
                      C2L = C1S(K,I+1,JC)
                      C3L = C1S(K,I+2,JC)
                      IF (K.GT.KB(I+2).OR.I.EQ.ID-1) C3L = C1S(K,I,JC)
                    END IF
                    CADL(K,I,JC) =  (DX1(K,I)-U(K,I)*AD1L(K,I))*C1L
     .                             +(DX2(K,I)-U(K,I)*AD2L(K,I))*C2L
     .                             +(DX3(K,I)-U(K,I)*AD3L(K,I))*C3L
                  END DO
                END DO

*************** Vertical advection terms

                DO I=IU,ID
                  DO K=KT,KB(I)-1
                    IF (W(K,I).GE.0.0) THEN
                      C1V = C1S(K-1,I,JC)
                      C2V = C1S(K,I,JC)
                      C3V = C1S(K+1,I,JC)
                      IF (K.LE.KT+1) C1V = C1S(KT,I,JC)
                    ELSE
                      C1V = C1S(K,I,JC)
                      C2V = C1S(K+1,I,JC)
                      C3V = C1S(K+2,I,JC)
                      IF (K.EQ.KB(I)-1) C3V = C1S(K+1,I,JC)
                    END IF
                    CADV(K,I,JC) = -W(K,I)*(AD1V(K,I)*C1V+AD2V(K,I)*C2V
     .                             +AD3V(K,I)*C3V)
                  END DO
                END DO
              END IF
            END DO
          END DO

********* Transport constituents

          DO JB=1,NBP
            IU = CUS(JB)
            ID = DS(JB)
            DO JAC=1,NAC
              JC = CN(JAC)
              IF (TRANSPORT(JC)) THEN
                DO I=IU,ID
                  C1(KT,I,JC) = ((C1S(KT,I,JC)*BHKT2(I)/DLT
     .                          +(CADL(KT,I,JC)*BHRKT1(I)
     .                          -CADL(KT,I-1,JC)*BHRKT1(I-1))/DLX(I)
     .                          +(1.0-THETA)*CADV(KT,I,JC)*BB(KT,I))
     .                          /BHKT1(I)+CSSB(KT,I,JC)/(DLX(I)
     .                          *BHKT1(I))+CSSK(KT,I,JC))*DLT
                  DO K=KT+1,KB(I)
                    C1(K,I,JC) = ((C1S(K,I,JC)*BH(K,I)/DLT
     .                           +(CADL(K,I,JC)*BHR(K,I)-CADL(K,I-1,JC)
     .                           *BHR(K,I-1))/DLX(I)+(1.0-THETA)
     .                           *(CADV(K,I,JC)*BB(K,I)-CADV(K-1,I,JC)
     .                           *BB(K-1,I)))/BH(K,I)+CSSB(K,I,JC)
     .                           /(DLX(I)*BH(K,I))+CSSK(K,I,JC))*DLT
                  END DO

***************** Crank-Nicholson vertical advection and implicit diffusion

                  DT(KB(I)) = C1(KB(I),I,JC)
                  DO K=KT,KB(I)-1
                    DT(K) = C1(K,I,JC)
                  END DO
                  GMAT(KT) = DT(KT)
                  DO K=KT+1,KB(I)
                    GMAT(K) = DT(K)-AT(K,I)/BTAT(K-1,I)*GMAT(K-1)
                  END DO
                  C1(KB(I),I,JC) = GMAT(KB(I))/BTAT(KB(I),I)
                  DO K=KB(I),KT,-1
                    C1(K,I,JC) = (GMAT(K)-CT(K,I)*C1(K+1,I,JC))
     .                           /BTAT(K,I)
c            minimum value of zooplankton=zoomin
             if(jc.eq.22.and.c1(k,i,jc).lt.zoomin)c1(k,i,jc)=zoomin
                  END DO
                END DO
              END IF
            END DO
          END DO
        END IF

************************************************************************
**        Task 7: Layer - Segment Additions and Subtractions          **
************************************************************************

******* Determine water surface minimum height

        ZMIN = -1000.0
        DO JB=1,NBP
          IU = CUS(JB)
          ID = DS(JB)
          DO I=IU,ID
            ZMIN = MAX(ZMIN,Z1(I))
          END DO
        END DO
c        ADD_LAYER = ZMIN.LT.-0.80*H(KT-1).AND.KT.NE.2
c        SUB_LAYER = ZMIN.GT.0.60*H(KT)
        ADD_LAYER = ZMIN.LT.-addkt*H(KT-1).AND.KT.NE.2
        SUB_LAYER = ZMIN.GT.subkt*H(KT)
        DO WHILE (ADD_LAYER)
          IF (SNAPSHOT) WRITE (SNP,4000) KT-1,JDAY

********* Initialize variables in new layer

          KT = KT-1
          DO JB=1,NBP
            IU = CUS(JB)
            ID = DS(JB)
            DO I=IU-1,ID+1
              Z1(I)     = H(KT)+Z1(I)
              z2(i)     = h(kt)+z2(i)
              HKT1(I)   = H(KT)-Z1(I)
              HEIGHT(I) = HEIGHT(I)-H(KT+1)
              BHKT1(I)  = BHKT1(I)-BH(KT+1,I)
              BKT(I)    = BHKT1(I)/HKT1(I)
              U(KT,I)   = U(KT+1,I)
              W(KT,I)   = W(KT+1,I)
              T1(KT,I)  = T1(KT+1,I)
              DO JC=1,NAC
                C1(KT,I,CN(JC))   = C1(KT+1,I,CN(JC))
                CSSK(KT,I,CN(JC)) = CSSK(KT+1,I,CN(JC))
              END DO
            END DO
            DO I=IU-1,ID
              AZMXKT(I) = 0.1*HKT1(I)*HKT1(I)/DLT
              BHRKT1(I) = (BHKT1(I)+BHKT1(I+1))*0.5
            END DO
            IF (UP_HEAD(JB)) THEN
              BHSUM           = BHRKT1(IU-1)+BHR(KT+1,IU-1)
              QUH1(KT,JB)     = QUH1(KT+1,JB)*BHRKT1(IU-1)/BHSUM
              QUH1(KT+1,JB)   = QUH1(KT+1,JB)*BHR(KT+1,IU-1)/BHSUM
              TSSUH1(KT,JB)   = TSSUH1(KT+1,JB)*BHRKT1(IU-1)/BHSUM
              TSSUH1(KT+1,JB) = TSSUH1(KT+1,JB)*BHR(KT+1,IU-1)/BHSUM
              DO JC=1,NAC
                CSSUH1(KT,CN(JC),JB)   = CSSUH1(KT+1,CN(JC),JB)
     .                                   *BHRKT1(IU-1)/BHSUM
                CSSUH1(KT+1,CN(JC),JB) = CSSUH1(KT+1,CN(JC),JB)
     .                                   *BHR(KT+1,IU-1)/BHSUM
              END DO
            END IF
            IF (DN_HEAD(JB)) THEN
              BHSUM           = BHRKT1(ID)+BHR(KT+1,ID)
              QDH1(KT,JB)     = QDH1(KT+1,JB)*BHRKT1(ID)/BHSUM
              QDH1(KT+1,JB)   = QDH1(KT+1,JB)*BHR(KT+1,ID)/BHSUM
              TSSDH1(KT,JB)   = TSSDH1(KT+1,JB)*BHRKT1(ID)/BHSUM
              TSSDH1(KT+1,JB) = TSSDH1(KT+1,JB)*BHR(KT+1,ID)/BHSUM
              DO JC=1,NAC
                CSSDH1(KT,CN(JC),JB)   = CSSDH1(KT+1,CN(JC),JB)
     .                                   *BHRKT1(ID)/BHSUM
                CSSDH1(KT+1,CN(JC),JB) = CSSDH1(KT+1,CN(JC),JB)
     .                                   *BHR(KT+1,ID)/BHSUM
              END DO
            END IF
            DO I=IU,ID-1
              DX(KT,I) = IDX
            END DO
            IDT = ID-1
            IF (DN_HEAD(JB)) IDT = ID
            if(up_head(jb)) iut=iu-1
c            DO I=IU,IDT
            DO I=iut,IDT
              AZ(KT,I) = AZMIN
            END DO

*********** Upstream active cell

            I = US(JB)
            DO WHILE (KB(I)-KT.LT.1.AND.I.LT.DS(JB))
              I = I+1
            END DO
            IUT = I

*********** Add segments

            IF (IUT.NE.IU) THEN
              IF (SNAPSHOT) WRITE (SNP,4010) IU,IUT,JDAY
              DO I=IUT-1,IU-1
                Z1(I)     = Z1(IU)
                z2(i)     = z2(iu)
                KTT(I)    = KTT(IU)
                HKT1(I)   = H(KT)-Z1(IU)
                BHKT1(I)  = 0.0
                HEIGHT(I) = HEIGHT(IU)
                DO K=KT,KTT(I)+1,-1
                  BHKT1(I) = BHKT1(I)+BH(K,I)
                END DO
                BHKT1(I)  = BHKT1(I)-B(KTT(I),I)*(HEIGHT(IU)+Z1(IU))
                BKT(I)    = BHKT1(I)/HKT1(I)
                AZMXKT(I) = 0.1*HKT1(I)*HKT1(I)/DLT
              END DO
              DO I=IUT-1,IU-1
                BHRKT1(I) = (BHKT1(I)+BHKT1(I+1))*0.5
              END DO
              DO I=IUT,IU-1
                ICE(I)   = ICE(IU)
                ICETH(I) = ICETH(IU)
                DO K=KT,KB(I)
                  DX(K,I) = IDX
                  T1(K,I) = T1(K,IU)
                  DO JC=1,NAC
                    BHT = BH(K,I)
                    IF (K.EQ.KT) BHT = BHKT1(I)
                    C1(K,I,CN(JC))   = C1(K,IU,CN(JC))
                    ICMBR(CN(JC),JB) = ICMBR(CN(JC),JB)
     .                                 +C1(K,IU,CN(JC))*DLX(I)*BHT
                  END DO
                END DO
                DO K=KT,KB(I)-1
                  AZ(K,I) = AZ(K,IU)
                END DO
              END DO
              DO K=KB(IUT),KB(IU)
                U(K,IU-1) = 0.0
              END DO
              IF (CONSTITUENTS.AND.SHIFT_DEMAND) THEN
                IDIFF = IUT-US(JB)
                DO I=IUT,ID
                  SOD(I) = SODS(I-IDIFF)
                END DO
              END IF
              IU      = IUT
              CUS(JB) = IU
              IF (UH_EXTERNAL(JB)) KB(IU-1) = KB(IU)
              IF (UH_INTERNAL(JB)) KB(IU-1) = MIN(KB(UHS(JB)),KB(IU))
              IF (JB.EQ.1) THEN
                IBPR = 1
                DO WHILE (IU.GT.IPR(IBPR))
                  IBPR = IBPR+1
                END DO
              END IF
            END IF
          END DO
          ZMIN = -1000.0
          DO JB=1,NBP
            DO I=CUS(JB),DS(JB)
              ZMIN = MAX(ZMIN,Z1(I))
            END DO
          END DO
          ADD_LAYER = ZMIN.LT.-0.80*H(KT-1).AND.KT.NE.2
        END DO
        DO WHILE (SUB_LAYER)
          IF (SNAPSHOT) WRITE (SNP,4020) KT,JDAY

********* Initialize variables in subtracted layer

          KT = KT+1
          DO JB=1,NBP
            IU = CUS(JB)
            ID = DS(JB)
            DO I=IU-1,ID+1
              A1(KT-1,I) = 0.0
            END DO

*********** Initialize variables in new layer

            DO I=IU-1,ID+1
              Z1(I)     = Z1(I)-H(KT-1)
              z2(i)     = z2(i)-h(kt-1)
              HKT1(I)   = H(KT)-Z1(I)
              BHKT1(I)  = BHKT1(I)+BH(KT,I)
              BKT(I)    = BHKT1(I)/HKT1(I)
              HEIGHT(I) = HEIGHT(I)+H(KT)
              AZMXKT(I) = 0.1*HKT1(I)*HKT1(I)/DLT
              T1(KT,I)  = (T1(KT-1,I)*(BHKT1(I)-BH(KT,I))
     .                   +T1(KT,I)*BH(KT,I))/BHKT1(I)
              DO JC=1,NAC
                JAC              = CN(JC)
                C1(KT,I,JAC)     = (C1(KT-1,I,JAC)*(BHKT1(I)-BH(KT,I))
     .                             +C1(KT,I,JAC)*BH(KT,I))/BHKT1(I)
                CSSB(KT,I,JAC)   = CSSB(KT-1,I,JAC)+CSSB(KT,I,JAC)
                CSSK(KT,I,JAC)   = (CSSK(KT-1,I,JAC)*(BHKT1(I)-BH(KT,I))
     .                             +CSSK(KT,I,JAC)*BH(KT,I))/BHKT1(I)
                CSSB(KT-1,I,JAC) = 0.0
              END DO
            END DO
            DO I=IU-1,ID
              BHRKT1(I) = (BHKT1(I)+BHKT1(I+1))*0.5
            END DO
            IF (UP_HEAD(JB)) THEN
              QUH1(KT,JB)   = QUH1(KT,JB)+QUH1(KT-1,JB)
              TSSUH1(KT,JB) = TSSUH1(KT-1,JB)+TSSUH1(KT,JB)
              DO JC=1,NAC
                CSSUH1(KT,CN(JC),JB) = CSSUH1(KT-1,CN(JC),JB)
     .                                 +CSSUH1(KT,CN(JC),JB)
              END DO
            END IF
            IF (DN_HEAD(JB)) THEN
              QDH1(KT,JB)   = QDH1(KT,JB)+QDH1(KT-1,JB)
              TSSDH1(KT,JB) = TSSDH1(KT-1,JB)+TSSDH1(KT,JB)
              DO JC=1,NAC
                CSSDH1(KT,CN(JC),JB) = CSSDH1(KT-1,CN(JC),JB)
     .                                 +CSSDH1(KT,CN(JC),JB)
              END DO
            END IF

*********** Upstream active cell

            I = US(JB)
            DO WHILE (KB(I)-KT.LT.1.AND.I.LT.DS(JB))
              I = I+1
            END DO
            IUT = I
            IF (IUT.GT.DS(JB)-2) THEN
              OPEN  (ERR,FILE='err.opt',STATUS='UNKNOWN')
              WRITE (ERR,7000) JB,JDAY,KT
              STOP
            END IF

*********** Subtract segments

            IF (IUT.NE.IU) THEN
              IF (SNAPSHOT) WRITE (SNP,4030) IU,IUT,JDAY
              DO I=IU-1,IUT-1
                F(I)     = 0.0
                Z1(I)    = 0.0
                z2(i)    = 0.0
                ICETH(I) = 0.0
                BHRHO(I) = 0.0
                ICE(I)   = .FALSE.
                DO K=KT,KB(I)
                  AZ(K,I)   = 0.0
                  DX(K,I)   = 0.0
                  U(K,I)    = 0.0
                  T1(K,I)   = 0.0
                  TADL(K,I) = 0.0
                  QSS(K,I)  = 0.0
                  DO JC=1,NAC
                    JAC = CN(JC)
                    BHT = BH(K,I)
                    IF (K.EQ.KT) BHT = BHKT1(I)
                    ICMBR(JAC,JB) = ICMBR(JAC,JB)-C1(K,I,JAC)*DLX(I)*BHT
                    C1(K,I,JAC)   = 0.0
                    C2(K,I,JAC)   = 0.0
                    CADL(K,I,JAC) = 0.0
                  END DO
                END DO
              END DO
              IF (CONSTITUENTS.AND.SHIFT_DEMAND) THEN
                IDIFF = IUT-US(JB)
                DO I=IUT,ID
                  SOD(I) = SODS(I-IDIFF)
                END DO
              END IF
              CUS(JB)      = IUT
              IU           = IUT
              Z1(IU-1)     = Z1(IU)
              z2(iu-1)     = z2(iu)
              KTT(IU-1)    = KTT(IU)
              HKT1(IU-1)   = HKT1(IU)
              BHKT1(IU-1)  = 0.0
              HEIGHT(IU-1) = HEIGHT(IU)
              DO K=KT,KTT(IU-1),-1
                BHKT1(IU-1) = BHKT1(IU-1)+B(K,IU-1)*H(K)
              END DO
              BHKT1(IU-1)  = BHKT1(IU-1)-B(KTT(IU-1),IU-1)
     .                       *(HEIGHT(IU-1)+Z1(IU-1))
              BKT(IU-1)    = BHKT1(IU-1)/HKT1(IU-1)
              BHRKT1(IU-1) = (BHKT1(IU-1)+BHKT1(IU))*0.5
              IF (UH_EXTERNAL(JB)) KB(IU-1) = KB(IU)
              IF (UH_INTERNAL(JB)) KB(IU-1) = MIN(KB(UHS(JB)),KB(IU))
              IF (JB.EQ.1) THEN
                IBPR = 1
                DO WHILE (IU.GT.IPR(IBPR))
                  IBPR = IBPR+1
                END DO
              END IF
            END IF
          END DO
          ZMIN = -1000.0
          DO JB=1,NBP
            DO I=CUS(JB),DS(JB)
              ZMIN = MAX(ZMIN,Z1(I))
            END DO
          END DO
          SUB_LAYER = ZMIN.GT.0.60*H(KT)
        END DO

************************************************************************
**                         Task 8: Balances                           **
************************************************************************

        IF (VOL_BALANCE) THEN
          VOLSG = 0.0
          VOLTG = 0.0
          DO JB=1,NBP
            VOLSBR(JB) = VOLSBR(JB)+DLVOL(JB)
            VOLTBR(JB) = VOLEV(JB)+VOLPR(JB)+VOLTR(JB)+VOLDT(JB)
     .                   +VOLWD(JB)+VOLUH(JB)+VOLDH(JB)+VOLIN(JB)
     .                   +VOLOUT(JB)
            IF (ABS(VOLSBR(JB)-VOLTBR(JB)).GT.VTOL) THEN
              IF (VOLUME_WARNING) THEN
                IF (.NOT.WARNING_OPEN) THEN
                  OPEN (WRN,FILE='wrn.opt',STATUS='UNKNOWN')
                  WARNING_OPEN = .TRUE.
                END IF
                WRITE (WRN,6010) JDAY,VOLSBR(JB),VOLTBR(JB),VOLSBR(JB)
     .                           -VOLTBR(JB)
                VOLUME_WARNING = .FALSE.
              END IF
            END IF
            VOLSG = VOLSG+VOLSBR(JB)
            VOLTG = VOLTG+VOLTBR(JB)
          END DO
        END IF
        IF (MASS_BALANCE) THEN
          DO JB=1,NBP
            DO JC=1,NAC
              JAC  = CN(JC)
              IF (TRANSPORT(JAC)) THEN
                CMBR(JAC,JB) = 0.0
                DO I=CUS(JB),DS(JB)
                  CMBR(JAC,JB)   = CMBR(JAC,JB)+C1(KT,I,JAC)*DLX(I)
     .                             *BHKT1(I)
                  CSSMBR(JAC,JB) = CSSMBR(JAC,JB)+(CSSB(KT,I,JAC)
     .                             +CSSK(KT,I,JAC)*BHKT1(I)*DLX(I))*DLT
                  DO K=KT+1,KB(I)
                    CMBR(JAC,JB)   = CMBR(JAC,JB)+C1(K,I,JAC)*DLX(I)
     .                               *BH(K,I)
                    CSSMBR(JAC,JB) = CSSMBR(JAC,JB)+(CSSB(K,I,JAC)
     .                               +CSSK(K,I,JAC)*BH(K,I)*DLX(I))*DLT
                  END DO
                END DO
              END IF
              DLMBR(JAC,JB) = ICMBR(JAC,JB)+CSSMBR(JAC,JB)
            END DO
          END DO
        END IF

************************************************************************
**            Task 9: Update Variables for Next Timestep              **
************************************************************************


        DO JB=1,NBP
          IU = CUS(JB)
          ID = DS(JB)
          DO I=IU-1,ID+1
            DPT       = Z1(I)
            Z1(I)     = Z2(I)
            Z2(I)     = DPT
            SZ1(I)    = Z1(I)
            SKTT(I)   = KTT(I)
            HKT2(I)   = HKT1(I)
            AVHKT(I)  = (HKT2(I)+H(KT+1))*0.5
            BHKT2(I)  = BHKT1(I)
            BHRKT2(I) = BHRKT1(I)
            IF (HKT2(I).LT.0.0) THEN
              OPEN  (ERR,FILE='err.opt',STATUS='UNKNOWN')
              WRITE (ERR,7010) JDAY,I,HKT2(I)
              STOP
            END IF
            DO K=KT,KB(I)
              T2(K,I)  = T1(K,I)
              SU(K,I)  = U(K,I)
              SW(K,I)  = W(K,I)
              QSS(K,I) = 0.0
              TSS(K,I) = 0.0
            END DO
          END DO
          IF (UH_INTERNAL(JB)) THEN
            DO K=KT,KB(IU-1)
              QUH2(K,JB)   = QUH1(K,JB)*DLT
              TSSUH2(K,JB) = TSSUH1(K,JB)*DLT
            END DO
          END IF
          IF (DH_INTERNAL(JB)) THEN
            DO K=KT,KB(ID+1)
              QDH2(K,JB)   = QDH1(K,JB)*DLT
              TSSDH2(K,JB) = TSSDH1(K,JB)*DLT
            END DO
          END IF
          DO JC=1,NAC
            DO I=IU-1,ID+1
              DO K=KT,KB(I)
                C1S(K,I,CN(JC))  = C1(K,I,CN(JC))
c tom says that a negative ensures mass conservation, so recommend
c keeping the negative quantity (if so) intact
                C2(K,I,CN(JC))   = MAX(C1(K,I,CN(JC)),0.0)
                CSSB(K,I,CN(JC)) = 0.0
              END DO
            END DO
            IF (UH_INTERNAL(JB)) THEN
              DO K=KT,KB(IU-1)
                CSSUH2(K,CN(JC),JB) = CSSUH1(K,CN(JC),JB)*DLT
              END DO
            END IF
            IF (DH_INTERNAL(JB)) THEN
              DO K=KT,KB(ID+1)
                CSSDH2(K,CN(JC),JB) = CSSDH1(K,CN(JC),JB)*DLT
              END DO
            END IF
          END DO
        END DO
        ELKT = EL(KT)-Z2(DS(1))

******* Time related variables

c open custom output files
        if(nit.eq.0)then
        open(193,file='uswl.opt',status='new')
        open(194,file='usct.opt',status='new')
        open(191,file='flows.opt',status='new')
        end if
c write out data every 200 time steps
       if((nit/100)*100.eq.nit)then
       write(193,3500)jday,3.2808*(el(kt)-z2(2)),3.2808*
     . (el(kt)-z2(10)),3.2808*(el(kt)-z2(119)),3.2808*
     . (el(kt)-z2(116)),3.2808*(el(kt)-z2(22)),3.2808*
     . (el(kt)-z2(140)),3.2808*(el(kt)-z2(29)),3.2808*
     . (el(kt)-z2(96)),3.2808*(el(kt)-z2(36)),3.2808*
     . (el(kt)-z2(209)),3.2808*(el(kt)-z2(43)),3.2808*
     . (el(kt)-z2(51)),3.2808*(el(kt)-z2(186)),3.2808*
     . (el(kt)-z2(199)),3.2808*
     . (el(kt)-z2(178)),3.2808*(el(kt)-z2(53)),
     . 3.2808*(el(kt)-z2(56)),3.2808*(el(kt)-z2(58)),
     . 3.2808*(el(kt)-z2(61)),3.2808*(el(kt)-z2(64)),
     . 3.2808*(el(kt)-z2(67)),3.2808*(el(kt)-z2(72)),
     . 3.2808*(el(kt)-z2(221))
3500   format(f7.3,1x,24(f5.2,1x))
       write(194,9504)jday,u(kt,4),u(kt,20),u(kt,34),
     . u(kt,51),u(kt,70),u(kt,83)
       write(194,9504)jday,t2(kt,2),
     . ((t2(kt,119)+t2(kt+1,119))/2.),((t2(kt,22)+t2(kt+1,22))/2.),
     . ((t2(kt,37)+t2(kt+1,37))/2.),((t2(kt,50)+t2(kt+1,50))/2.),
     . ((t2(kt,221)+t2(kt+1,221))/2.),((t2(kt,91)+t2(kt+1,91))/2.),
     . ((t2(kt,96)+t2(kt+1,96))/2.),((t2(kt,140)+t2(kt+1,140))/2.),
     . ((t2(kt,180)+t2(kt+1,180))/2.),((t2(kt,196)+t2(kt+1,196))/2.)
       do jc=1,nac
       jx=cn(jc)
       write(194,9504)jday,((c2(kt,2,jx)+c2(kt+1,2,jx))/2.),
     . ((c2(kt,119,jx)+c2(kt+1,119,jx))/2.),
     . ((c2(kt,22,jx)+c2(kt+1,22,jx))/2.),
     . ((c2(kt,37,jx)+c2(kt+1,37,jx))/2.),
     . ((c2(kt,50,jx)+c2(kt+1,50,jx))/2.),
     . ((c2(kt,221,jx)+c2(kt+1,221,jx))/2.),
     . ((c2(kt,91,jx)+c2(kt+1,91,jx))/2.),
     . ((c2(kt,96,jx)+c2(kt+1,96,jx))/2.),
     . ((c2(kt,140,jx)+c2(kt+1,140,jx))/2.),
     . ((c2(kt,180,jx)+c2(kt+1,180,jx))/2.),
     . ((c2(kt,196,jx)+c2(kt+1,196,jx))/2.)
       end do
9504   format(f8.3,11(2x,e15.4))
        end if
c output file for mcdd1 output for lower slough model
        if(nit.eq.0)then
        open(195,file='qtmcdd1.opt',status='new')
        open(196,file='cmcdd1.opt',status='new')
        write(195,9500)
        write(196,9501)
        end if
9500    format(1x,'qtmcdd1 file for flow and temp from mcdd1'/1x/4x,
     .  'jday', 7x,'q',4x,'temp')
9503    FORMAT(F8.3,F8.3,F8.3)
9501    FORMAT(1X,'cmcdd1 inflow conc file '/1X/
     .    3X,'JDAY',1X,'TRACER',1X,'SSOLID',1X,'CLFORM',1X,
     .   'DSOLID',2X,'BOD-L',2X,'BOD-R',2X,'ALGAE',1X,'DETRIT',
     .   4X,'PO4',4X,'NH4',4X,'NO3',5X,'O2',2X,'SDMNT',1X,'CARBON',
     .   4X,'ALK',5X,'PH',5X,'FE',3x,'ZOOP')
9502     FORMAT(F7.2,F7.2,F7.2,F7.0,F7.2,F7.2,F7.1,F7.2,F7.2,3(F7.3),
     .    F7.2,f7.3,F7.2,F7.2,F7.2,f7.2,f7.3)
        if((nit/200)*200.eq.nit)then
        write(195,9503)jday,(qwd(1)+qwd(29)),
     .   ((t2(kt,2)+t2(kt+1,2))/2.)
        write(196,9502)jday,((c2(kt,2,1)+c2(kt+1,2,1))/2.),
     .  ((c2(kt,2,2)+c2(kt+1,2,2))/2.),((c2(kt,2,3)+c2(kt+1,2,3))/2.),
     .  ((c2(kt,2,4)+c2(kt+1,2,4))/2.),((c2(kt,2,5)+c2(kt+1,2,5))/2.),
     .  ((c2(kt,2,6)+c2(kt+1,2,6))/2.),((c2(kt,2,7)+c2(kt+1,2,7))/2.),
     .  ((c2(kt,2,8)+c2(kt+1,2,8))/2.),((c2(kt,2,9)+c2(kt+1,2,9))/2.),
     .  ((c2(kt,2,10)+c2(kt+1,2,10))/2.),
     .  ((c2(kt,2,11)+c2(kt+1,2,11))/2.),
     .  ((c2(kt,2,12)+c2(kt+1,2,12))/2.),
     .  ((c2(kt,2,13)+c2(kt+1,2,13))/2.),
     .  ((c2(kt,2,14)+c2(kt+1,2,14))/2.),
     .  ((c2(kt,2,15)+c2(kt+1,2,15))/2.),
     .  ((c2(kt,2,16)+c2(kt+1,2,16))/2.),
     .  ((c2(kt,2,20)+c2(kt+1,2,20))/2.),
     .  ((c2(kt,2,22)+c2(kt+1,2,22))/2.)
        end if

c  writing out z1 - chris
c        if(jday.gt.232.0.and.jday.lt.232.5)then
          elevc=3.2808*(el(kt)-z2(2))
          if(elevc.gt.7.80.and.elevc.lt.7.85)then
            open(99,file='z1.opt',status='unknown')
            write(99,*)kt
            write(99,'(10f8.3)')(z2(iii), iii=1,imp)
            close(99)
          end if
c        end if

c  flow data write-out - chris
        qmcdd=(qwd(1)+qwd(29))/(0.3048**3)
         qmcdd4=qwd(40)/(0.3048**3)  
        if((nit/20)*20.eq.nit)then
          qpris=-qtr(28)/(0.3048**3)
          if(qpris.eq.0.0)then
            qpris=0.0
            jb=28
            DO JO=1,NOUT(JB)
              IF (KOUT(JO,JB).GE.KT) qpris = qpris+
     .           QOUT(JO,JB)/(0.3048**3)
            END DO
          end if

          q82nd=-qtr(1)/(0.3048**3)
          if(q82nd.eq.0.0)then
            q82nd=0.0
            jb=1
            DO JO=1,NOUT(JB)
              IF (KOUT(JO,JB).GE.KT) q82nd = q82nd+
     .          QOUT(JO,JB)/(0.3048**3)
            END DO
          end if

          q92nd=-qtr(24)/(0.3048**3)
          if(q92nd.eq.0.0)then
            q92nd=0.0
            jb=23
            DO JO=1,NOUT(JB)
              IF (KOUT(JO,JB).GE.KT) q92nd = q92nd+
     .          QOUT(JO,JB)/(0.3048**3)
            END DO
          end if

          q47th=-qtr(13)/(0.3048**3)
          if(q47th.eq.0.0)then
            q47th=0.0
            jb=13
            DO JO=1,NOUT(JB)
              IF (KOUT(JO,JB).GE.KT) q47th = q47th+
     .          QOUT(JO,JB)/(0.3048**3)
            END DO
          end if


          q33rd=-qtr(9)/(0.3048**3)
          if(q33rd.eq.0.0)then
            q33rd=0.0
            jb=9
            DO JO=1,NOUT(JB)
              IF (KOUT(JO,JB).GE.KT) q33rd = q33rd+
     .                        QOUT(JO,JB)/(0.3048**3)
            END DO
          end if

          qdike=-qtr(3)/(0.3048**3)
          if(qdike.eq.0.0)then
            qdike=0.0
            jb=3
            DO JO=1,NOUT(JB)
              IF (KOUT(JO,JB).GE.KT)qdike = qdike+
     .                  QOUT(JO,JB)/(0.3048**3)
            END DO
          end if
        write(191,'(f8.3,10f8.2)')jday,qmcdd,qpris,q82nd,q92nd,
     .     q47th,q33rd,aqbox,aqppipe,qdike,qmcdd4
        end if


c  special tracer output files - chris
      open(175,file='trace1.opt',status='new')
      open(176,file='trace2.opt',status='new')
      open(177,file='trace3.opt',status='new')
c write out data every 200 time steps
       if((nit/100)*100.eq.nit)then
       
       jx=cn(1)
       do itt=1,22
         if(itt.eq.1)icl=34
         if(itt.eq.2)icl=25
         if(itt.eq.3)icl=17
         if(itt.eq.4)icl=16
         if(itt.eq.5)icl=15
         if(itt.eq.6)icl=14
         if(itt.eq.7)icl=4
         if(itt.eq.8)icl=3
         if(itt.eq.9)icl=115
         if(itt.eq.10)icl=85
         if(itt.eq.11)icl=83
         if(itt.eq.12)icl=2
         if(itt.eq.13)icl=81
         if(itt.eq.14)icl=219
         if(itt.eq.15)icl=57
         if(itt.eq.16)icl=211
         if(itt.eq.17)icl=50
         if(itt.eq.18)icl=210
         if(itt.eq.19)icl=205
         if(itt.eq.20)icl=190
         if(itt.eq.21)icl=159
         if(itt.eq.22)icl=122
         
         nlay=0
         cctr(itt)=0.0
         do k=kt,kb(icl)
           cctr(itt)=cctr(itt)+c2(k,icl,jx)
           nlay=nlay+1
         end do
         cctr(itt)=cctr(itt)/real(nlay)
       end do

       write(175,9561)jday,(cctr(ii),ii=1,8)
       write(176,9561)jday,(cctr(ii),ii=9,16) 
       write(177,9562)jday,(cctr(ii),ii=17,22)

9561   format(f8.3,8(2x,e15.4))
 9562  format(f8.3,6(2x,e15.4))
        end if

c  special write-out to track bod and do plume created by 
c airport deicing
      open(185,file='boddo.opt',status='new')

c write out data every 200 time steps
       if(jday.ge.45.0.and.jday.le.50.0)then
       
       if(jday.ge.boddonx)then
       boddonx=boddonx+0.1
       do icl=1,imp
         
         jx=cn(12)
         nlay=0
         ccdo(icl)=0.0
         do k=kt,kb(icl)
           ccdo(icl)=ccdo(icl)+c2(k,icl,jx)
           nlay=nlay+1
         end do
         ccdo(icl)=ccdo(icl)/real(nlay)

         jx=cn(5)
         nlay=0
         ccbod(icl)=0.0
         do k=kt,kb(icl)
           ccbod(icl)=ccbod(icl)+c2(k,icl,jx)
           nlay=nlay+1
         end do
         ccbod(icl)=ccbod(icl)/real(nlay)

       end do
        write(185,9561)jday
       write(185,9563)(ccdo(ii),ii=1,imp)
       write(185,9561)jday
       write(185,9563)(ccbod(ii),ii=1,imp)

       end if

       end if

9563   format(e15.4,8(2x,e15.4))


        NIT     = NIT+1
        ELTM    = ELTM+DLT
        JDAY    = ELTM/86400.0
        ELTMJD  = JDAY-TMSTRT
        END_RUN = JDAY.GE.TMEND
        DLTS    = DLT
        DLT     = INT(MAX(DLTMIN,DLTF(DLTDP)*CURMAX))
        DLT     = MIN(DLT,DLTS+1.0*DLTS)
        AVDLT   = (ELTM-TMSTRT*86400.0)/NIT
        IF (JDAY.GE.WSCD(WSCDP+1)) WSCDP = WSCDP+1
        IF (JDAY.GE.DLTD(DLTDP+1)) DLTDP = DLTDP+1
        IF (DLT.GT.DLTMAX(DLTDP))  DLT   = DLTMAX(DLTDP)
        IF (DLTS.LT.MINDLT) THEN
          MINDLT = DLTS
          JDMIN  = JDAY
        END IF
        CURMAX = DLTMAX(DLTDP)/DLTF(DLTDP)
        IF (INT(JDAY).EQ.JDAYNX) THEN
          JDAYG  = JDAYG+1
          JDAYNX = JDAYNX+1
        END IF
        IF (JDAYG.GT.300) WINTER = .TRUE.
        IF (JDAYG.LT.40)  WINTER = .FALSE.
        WRITE (GDCH,'(I3)') GDAY
        CALL GREGORIAN_DATE (YEAR)
        update_kinetics = .false.
        IF (MOD(NIT,FREQUK).EQ.0) UPDATE_KINETICS = .TRUE.

************************************************************************
**                      Task 10: Output Results                       **
************************************************************************
c  calling environmental performance criterion subroutine
        IF(JDAY.GE.NXTEP)THEN
          NXTEP=NXTEP+0.02083
          do jb=1,nbp
            iu=cus(jb)
            id=ds(jb)
*******     Branch volumes
            VOLB(JB)=0.0
              DO I=IU,ID
               VOLB(JB) = VOLB(JB)+DLX(I)*BHKT2(I)
                 DO K=KT+1,KB(I)
                 VOLB(JB) = VOLB(JB)+DLX(I)*B(K,I)*H(K)
                 END DO
              END DO
      if(jb.eq.4.or.jb.eq.5.or.jb.eq.6.or.jb.eq.7.or.jb.eq.8.or.
     .   jb.eq.33)then
c reach 3 envir performance criterion
         CALL ENVIRP2(tmstrt,jday,end_run,volb(jb),dlx)
         else
c reach 2 envir perfromance criterion
         CALL ENVIRP(tmstrt,jday,end_run,volb(jb),dlx)
      end if

          end do
        call mass_load
        END IF

******* Snapshots

        IF (SNAPSHOT) THEN
          IF (JDAY.GE.NXTMSN.OR.JDAY.GE.SNPD(SNPDP+1)) THEN
            IF (JDAY.GE.SNPD(SNPDP+1)) THEN
              SNPDP  = SNPDP+1
              NXTMSN = SNPD(SNPDP)
            END IF
            NXTMSN = NXTMSN+SNPF(SNPDP)
            DO JB=1,NBP
              QSUM(JB) = 0.0
              DO JO=1,NOUT(JB)
                IF (KOUT(JO,JB).GE.KT) QSUM(JB) = QSUM(JB)+QOUT(JO,JB)
              END DO
            END DO
            WRITE (SNP,3000) (TITLE(I),I=1,5)
            WRITE (SNP,3010) 'Time Parameters',MONTH,GDAY,YEAR,
     .                       INT(JDAY),(JDAY-INT(JDAY))*24.0,
     .                       INT(ELTMJD),(ELTMJD-INT(ELTMJD))*24.0,
     .                       INT(DLTS),KLOC,ILOC,INT(MINDLT),INT(JDMIN),
     .                       (JDMIN-INT(JDMIN))*24.0,KMIN,IMIN,
     .                       INT(AVDLT),NIT,NV
            IF(TERM_BY_TERM) THEN
             WRITE (SNP,3020) 'Meteorlogical Parameters',TAIR,TDEW,
     .                        WIND*1.0/WSC(WSCDP),
     .                        PHI,WSC(WSCDP),CLOUD,ET,CSHE,SRO
            ELSE
             WRITE (SNP,3025) 'Meteorlogical Parameters',ET,TDEW,
     .        WIND*1.0/WSC(WSCDP),
     .                        PHI,WSC(WSCDP),CLOUD,CSHE,SRO
            END IF
            WRITE (SNP,3030) 'Inflows','Upstream inflows'
            WRITE (SNP,3040) (JB,QIN(JB),TIN(JB),JB=1,NBP)
            DO JB=1,NBP
              IF (DIST_TRIBS(JB)) THEN
                WRITE (SNP,3050)
                WRITE (SNP,3060) JB,QDTR(JB),TDTR(JB)
              END IF
            END DO
            WRITE (SNP,3070) (ITR(JT),JT=1,NTR)
            IF (PLACE_QTR) THEN
              WRITE (SNP,3080) (KTR(JT),JT=1,NTR)
            ELSE
              WRITE (SNP,3085) (KT,KB(ITR(JT)),JT=1,NTR)
            END IF
            WRITE (SNP,3090) (QTR(JT),JT=1,NTR)
            WRITE (SNP,3100) (TTR(JT),JT=1,NTR)
            DO JB=1,NBP
              IF (DN_FLOW(JB)) THEN
                WRITE (SNP,3110) QSUM(JB),(KOUT(JO,JB),JO=1,NOUT(JB))
                WRITE (SNP,3120) (QOUT(JO,JB),JO=1,NOUT(JB))
              END IF
            END DO
            WRITE (SNP,3130) (IWD(JW),JW=1,NWD)
            WRITE (SNP,3140) (KWD(JW),JW=1,NWD)
            WRITE (SNP,3150) (QWD(JW),JW=1,NWD)
            IF (CONSTITUENTS) THEN
              WRITE (SNP,3160) 'Constituent Inflow Concentrations'
              DO JB=1,NBP
                WRITE (SNP,3170) JB,(CNAME(CN(JC)),CIN(CN(JC),JB),
     .                           JC=1,NAC)
              END DO
              DO JT=1,NTR
                WRITE (SNP,3180) JT,(CNAME(CN(JC)),CTR(CN(JC),JT),
     .                           JC=1,NAC)
              END DO
              DO JB=1,NBP
                IF (DIST_TRIBS(JB)) THEN
                  WRITE (SNP,3190) JB,(CNAME(CN(JC)),
     .                             CDTR(CN(JC),JB),JC=1,NAC)
                END IF
              END DO
            END IF
            IF (EVAPORATION.OR.PRECIPITATION) WRITE(SNP,3200)
     .        'Surface Calculations'
            IF (EVAPORATION)   WRITE (SNP,3210) (JB,EVBR(JB),JB=1,NBP)
            IF (PRECIPITATION) WRITE (SNP,3220) (JB,QPRBR(JB),JB=1,NBP)
            IF (HEAD_BOUNDARY) THEN
              WRITE (SNP,3230)
              DO JB=1,NBP
                IF (UH_EXTERNAL(JB)) WRITE (SNP,3240) JB,ELUH(JB)
                IF (DH_EXTERNAL(JB)) WRITE (SNP,3250) JB,ELDH(JB)
              END DO
            END IF
            IF (VOL_BALANCE) THEN
              WRITE (SNP,3260)
              DO JB=1,NBP
                IF (VOLTBR(JB).NE.0.0) DLVR = (VOLTBR(JB)-VOLSBR(JB))
     .                                        /VOLTBR(JB)
                WRITE (SNP,3270) JB,VOLTBR(JB),VOLSBR(JB),VOLTBR(JB)
     .                           -VOLSBR(JB),DLVR*100.0
              END DO
              IF (VOLTG.NE.0.0) DLVR = (VOLTG-VOLSG)/VOLTG
              WRITE (SNP,3280) VOLTG,VOLSG,VOLTG-VOLSG,DLVR*100.0
            END IF
            IF (MASS_BALANCE) THEN
              WRITE (SNP,3290)
              DO JB=1,NBP
                WRITE (SNP,3300) JB
                DO JC=1,NAC
                  JAC = CN(JC)
                  IF (TRANSPORT(JAC)) THEN
                    IF (CMBR(JAC,JB).NE.0.0) THEN
                      DLMR = (CMBR(JAC,JB)-DLMBR(JAC,JB))/CMBR(JAC,JB)
                    END IF
                    WRITE (SNP,3310) CNAME(JAC),CMBR(JAC,JB)/1000.0,
     .                              (CMBR(JAC,JB)-DLMBR(JAC,JB))/1000.0,
     .                               DLMR*100.0
                  END IF
                END DO
              END DO
            END IF
            WRITE (SNP,3320) 'Geometry',KT,ELKT,(JB,CUS(JB),JB=1,NBP)
            CALL PRINT_GRID (JDAY,GDAY,MONTH,YEAR)
          END IF
        END IF

******* Vertical profiles

        IF (PROFILE) THEN
          IF (JDAY.GE.NXTMPR.OR.JDAY.GE.PRFD(PRFDP+1)) THEN
            IF (JDAY.GE.PRFD(PRFDP+1)) THEN
              PRFDP  = PRFDP+1
              NXTMPR = PRFD(PRFDP)
            END IF
            NXTMPR = NXTMPR+PRFF(PRFDP)
            NSPRF  = NSPRF+1
            WRITE (PRF,2570) JDAY,MONTH,GDAY,YEAR,KT,SNGL(Z1(DS(1))),
     .                       NSPRF
            DO JC=1,NAC
              IF (CPRN(CN(JC)).EQ.' ON') THEN
                DO JPRF=1,NIPRF
                  I   = IPRF(JPRF)
                  NRS = KB(I)-KT+1
                  WRITE (PRF,2560) CN(JC),NRS,(C2(K,I,CN(JC)),
     .                             K=KT,KB(I))
                END DO
              END IF
            END DO
            DO JPRF=1,NIPRF
              I   = IPRF(JPRF)
              NRS = KB(I)-KT+1
              WRITE (PRF,2560) 22,NRS,(T2(K,I),K=KT,KB(I))
            END DO
          END IF
        END IF

******* Velocity vectors

        IF (VECTOR) THEN
          IF (JDAY.GE.NXTMVP.OR.JDAY.GE.VPLD(VPLDP+1)) THEN
            IF (JDAY.GE.VPLD(VPLDP+1)) THEN
              VPLDP  = VPLDP+1
              NXTMVP = VPLD(VPLDP)
            END IF
            NXTMVP = NXTMVP+VPLF(VPLDP)
            WRITE (VPL,*) JDAY,MONTH//GDCH//', ',YEAR,KT,US
            WRITE (VPL,*) Z2
            WRITE (VPL,*) U
            WRITE (VPL,*) W
          END IF
        END IF

******* Contours

        IF (CONTOUR) THEN
          IF (JDAY.GE.NXTMCP.OR.JDAY.GE.CPLD(CPLDP+1)) THEN
            IF (JDAY.GE.CPLD(CPLDP+1)) THEN
              CPLDP  = CPLDP+1
              NXTMCP = CPLD(CPLDP)
            END IF
            NXTMCP = NXTMCP+CPLF(CPLDP)
c            WRITE (CPL,*) JDAY,MONTH//GDCH//', ',YEAR,KT,US
            WRITE (CPL,*) JDAY,YEAR,KT,US,z2
            DO JB=1,NBP
              DO I=US(JB),DS(JB)
                WRITE (CPL,8000) (T2(K,I),K=KT,KB(I))
                DO JC=1,NAC
                  IF (CPRN(CN(JC)).EQ.' ON') THEN
                    WRITE (CPL,8000) (C2(K,I,CN(JC)),K=KT,KB(I))
                  END IF
                END DO
              END DO
            END DO
          END IF
        END IF

******* Time series

        IF (TIME_SERIES) THEN
          IF (JDAY.GE.NXTMTS.OR.JDAY.GE.TSRD(TSRDP+1)) THEN
            IF (JDAY.GE.TSRD(TSRDP+1)) THEN
              TSRDP  = TSRDP+1
              NXTMTS = TSRD(TSRDP)
            END IF
            NXTMTS = NXTMTS+TSRF(TSRDP)
            DO JB=1,NBP
              IF (DN_FLOW(JB)) THEN
                QSUM(JB) = 0.0
                AVTOUT   = 0.0
                DO JC=1,NAC
                  AVCOUT(CN(JC)) = 0.0
                END DO
                DO JO=1,NOUT(JB)
                  K = KOUT(JO,JB)
                  IF (K.GE.KT) THEN
                    QSUM(JB) = QSUM(JB)+QOUT(JO,JB)
                    AVTOUT   = AVTOUT+T2(K,DS(JB))*QOUT(JO,JB)
                    DO JC=1,NAC
                      AVCOUT(CN(JC)) = AVCOUT(CN(JC))+QOUT(JO,JB)
     .                                 *C2(K,DS(JB),CN(JC))
                    END DO
                  END IF
                END DO
                IF (QSUM(JB).GT.0.0) THEN
                  AVTOUT = AVTOUT/QSUM(JB)
                  DO JC=1,NAC
                    AVCOUT(CN(JC)) = AVCOUT(CN(JC))/QSUM(JB)
                  END DO
                END IF
              END IF
            END DO
            IF (WITHDRAWALS) THEN
              DO JW=1,NWP
                TWD(JW) = T2(KWD(JW),IWD(JW))
                DO JC=1,NAC
                  CWD(JW,CN(JC)) = C2(KWD(JW),IWD(JW),CN(JC))
                END DO
              END DO
            END IF
            WRITE (TSR,5030) JDAY,DLT,ELKT,AVTOUT,(AVCOUT(CN(JC)),
     .                       JC=1,NAC)
          END IF
        END IF

******* Restart

        IF (RESTART_OUT) THEN
          IF (JDAY.GE.NXTMRS.OR.JDAY.GE.RSOD(RSODP+1)) THEN
            IF (JDAY.GE.RSOD(RSODP+1)) THEN
              RSODP  = RSODP+1
              NXTMRS = RSOD(RSODP)
            END IF
            NXTMRS = NXTMRS+RSOF(RSODP)
            IDAY   = INT(JDAY)
            IF (IDAY.LT.10) THEN
              WRITE (EXT1,'(I1)') IDAY
              RSOFN = 'rso'//EXT1//'.opt'
            ELSE IF (IDAY.LT.100) THEN
              WRITE (EXT2,'(I2)') IDAY
              RSOFN = 'rso'//EXT2//'.opt'
            ELSE IF (IDAY.LT.1000) THEN
              WRITE (EXT3,'(I3)') IDAY
              RSOFN = 'rso'//EXT3//'.opt'
            ELSE IF (IDAY.LT.10000) THEN
              WRITE (EXT4,'(I4)') IDAY
              RSOFN = 'rso'//EXT4//'.opt'
            ELSE
              WRITE (EXT5,'(I5)') IDAY
              RSOFN = 'rso'//EXT5//'.opt'
            END IF
            OPEN  (RSO,FILE=RSOFN,STATUS='UNKNOWN',FORM='UNFORMATTED')
            WRITE (RSO) KT,     NIT,    NV,     KMIN,   IMIN
            WRITE (RSO) DLTDP,  SNPDP,  TSRDP,  VPLDP,  PRFDP,  CPLDP,
     .                  RSODP,  WSCDP
            WRITE (RSO) JDAY,   JDAYG,  YEAR,   ELTM,   DLT,    AVDLT,
     .                  MINDLT, JDMIN,  CURMAX
            WRITE (RSO) NXTMSN, NXTMTS, NXTMPR, NXTMCP, NXTMVP, NXTMRS
            WRITE (RSO) VOLIN,  VOLOUT, VOLUH,  VOLDH,  VOLPR,  VOLTR,
     .                  VOLDT,  VOLWD,  VOLEV,  VOLSBR, ICMBR,  CSSMBR
            WRITE (RSO) TSSUH2, TSSDH2, CSSUH2, CSSDH2, QUH2,   QDH2
            WRITE (RSO) Z1,     Z2,     SZ1,    SKTT,   ICE,    ICETH
            WRITE (RSO) U,      W,      SU,     SW,     AZ,     T1,
     .                  T2,     C1,     C2
            CLOSE (RSO)
          END IF
        END IF
      END DO

************************************************************************
**                     Task 11: End Simulation                        **
************************************************************************
      CALL ENVIRP(tmstrt,jday,end_run,volb(1),dlx)
      CALL ENVIRP2(tmstrt,jday,end_run,volb(1),dlx)
      call mass_load
      IF (SNAPSHOT)     CLOSE (SNP)
      IF (TIME_SERIES)  CLOSE (TSR)
      IF (VECTOR)       CLOSE (VPL)
      IF (PROFILE)      CLOSE (PRF)
      IF (CONTOUR)      CLOSE (CPL)
      IF (WARNING_OPEN) CLOSE (WRN)

***** Profile FORMATs

 2500 FORMAT(A72)
 2510 FORMAT(20(1X,A3))
 2520 FORMAT(3(1X,A16,1X,A6))
 2530 FORMAT(L2,(:/20I4))
 2540 FORMAT(20I4)
 2550 FORMAT(10F8.2)
 2560 FORMAT(2I4/(8(1PE10.2E2)))
 2570 FORMAT(F8.3,A9,I3,', ',2I4,F8.4,I4)

***** Snapshot FORMATs

 3000 FORMAT('1',5(A72/1X))
 3010 FORMAT(1X,A15/
     .       '+',4('_'),1X,10('_')//
     .       3X,'Gregorian date      [GDAY] =',A16,I3,',',I5/
     .       3X,'Julian date         [JDAY] =',I8,' days',F6.2,
     .         ' hours'/
     .       3X,'Elapsed time      [ELTMJD] =',I8,' days',F6.2,
     .         ' hours'/
     .       3X,'Timestep             [DLT] =',I8,' sec'/
     .       5X,'at location  [KLOC,ILOC] = (',I2,',',I3,')'/
     .       3X,'Minimum timestep  [MINDLT] =',I8,' sec '/
     .       5X,'at Julian day    [JDMIN] =',I8,' days',F6.2,' hours'/
     .       5X,'at location  [KMIN,IMIN] = (',I2,',',I3,')'/
     .       3X,'Average timestep   [AVDLT] =',I8,' sec'/
     .       3X,'Number of iterations [NIT] =',I8/
     .       3X,'Number of violations  [NV] =',I8/)
 3020 FORMAT(1X,A24/
     .       '+',13('_'),1X,10('_')//
     .       3x,'Input'/
     .       5X,'Air temperature       [TAIR] =',F9.2,' deg C'/
     .       5X,'Dewpoint temperature  [TDEW] =',F9.2,' deg C'/
     .       5X,'Wind speed            [WIND] =',F9.2,' m/sec'/
     .       5X,'Wind direction         [PHI] =',F9.2,' rad'/
     .       5X,'Wind sheltering        [WSC] =',F9.2/
     .       5X,'Cloud cover          [CLOUD] =',F9.2/
     .       3X,'Calculated'/
     .       5X,'Equilibrium temperature [ET] =',F9.2,' deg C'/
     .       5X,'Surface heat exchange [CSHE] =',1PE9.2,' m/sec'/
     .       5X,'Short wave radiation   [SRO] =',1PE9.2,' c-m/sec'/)
 3025 FORMAT(1X,A24/
     .       '+',13('_'),1X,10('_')//
     .       3x,'Input'/
     .       5X,'Equilibrium temperature [ET] =',F9.2,' deg C'/
     .       5X,'Dewpoint temperature  [TDEW] =',F9.2,' deg C'/
     .       5X,'Wind speed            [WIND] =',F9.2,' m/sec'/
     .       5X,'Wind direction         [PHI] =',F9.2,' rad'/
     .       5X,'Wind sheltering        [WSC] =',F9.2/
     .       5X,'Cloud cover          [CLOUD] =',F9.2/
     .       5X,'Surface heat exchange [CSHE] =',1PE9.2,' m/sec'/
     .       5X,'Short wave radiation   [SRO] =',1PE9.2,' c-m/sec'/)
 3030 FORMAT(1X,A7/
     .       '+',7('_')//
     .       3X,A16)
 3040 FORMAT(5X,'Branch',I3/
     .       7X,'Inflow       [QIN] =',F8.2,' m^3/sec'/
     .       7X,'Temperature  [TIN] =',F8.2,' deg C')
 3050 FORMAT(/3X,'Distributed Tributaries'/
     .       '+',2X,11('_'),1X,11('_'))
 3060 FORMAT(5X,'Branch',I3/
     .       7X,'Inflow      [QDTR] =',F8.2,' m^3/sec'/
     .       7X,'Temperature [TDTR] =',F8.2,' deg C')
 3070 FORMAT('+'://3X,'Tributaries'/
     .       '+',2X,11('_')/
     .       5X,'Segment     [ITR] =',8I8/(T26,8I8))
 3080 FORMAT('+':/5X,'Layer       [KTR] =',8I8/
     .      (T26,8I8))
 3085 FORMAT('+':/5X,'Layer       [KTR] =',8(I5,'-',I2)/
     .      (T26,8(I5,'-',I2)))
 3090 FORMAT('+':/5X,'Inflow      [QTR] =',8F8.1/
     .      (T26,8F8.1))
 3100 FORMAT('+':/5X,'Temperature [TTR] =',8F8.1/
     .      (T26,8F8.1))
 3110 FORMAT('+'://
     .       1X,'Total outflow =',F8.2,' m^3/s'//
     .       3X,'Outlets'/
     .       '+',2X,7('_')/
     .       (5X,'Layer             [KOUT] =',9I8))
 3120 FORMAT('+':/
     .       (5X,'Outflow (m^3/sec) [QOUT] =',9F8.2))
 3130 FORMAT('+'://
     .       3X,'Withdrawals'/
     .       '+',11('_')/
     .       5X,'Segment           [IWD] =',9I7,(:/29X,9I7))
 3140 FORMAT(:5X,'Layer             [KWD] =',9I7,(:/29X,9I7))
 3150 FORMAT(:5X,'Outflow (m^3/sec) [QWD] =',9F7.2,(:/29X,9F7.2))
 3160 FORMAT('1',A33/
     .       '+',11('_'),1X,6('_'),1X,14('_')/)
 3170 FORMAT(3X,'Branch',I3/
     .      (5X,A16,T24,'=',F8.3,' g/m^3'))
 3180 FORMAT('+':/
     .       3X,'Tributary',I3/
     .      (5X,A16,T24,'=',F8.3,' g/m^3'))
 3190 FORMAT('+':/
     .       3X,'Distributed tributary',I3/
     .      (5X,A16,T24,'=',F8.3,' g/m^3'))
 3200 FORMAT(1X,A20/
     .       '+',7('_'),1X,12('_')/)
 3210 FORMAT(3X,'Evaporation'/
     .       '+',2X,11('_')/
     .      (:/5X,'Branch',I3,' =',F8.2))
 3220 FORMAT(3X,'Precipitation'/
     .       '+',2X,13('_')//
     .      (5X,'Branch',I3,' =',F8.2/))
 3230 FORMAT(/1X,'External head boundary elevations'/
     .       '+',8('_'),1X,4('_'),1X,8('_'),1X,10('_')/)
 3240 FORMAT(3X,'Branch',I3/
     .       5X,'Upstream elevation   [ELUH] =',F8.3,' m')
 3250 FORMAT(3X,'Branch',I3/
     .       5X,'Downstream elevation [ELDH] =',F8.3,' m')
 3260 FORMAT(/1X,'Water Balance'/
     .       '+',5('_'),1X,7('_')/)
 3270 FORMAT(3X,'Branch',I3/
     .       5X,'Spatial change  [VOLSBR] = ',1PE15.8E2,' m^3'/
     .       5X,'Temporal change [VOLTBR] = ',1PE15.8E2,' m^3'/
     .       5X,'Volume error             = ',1PE15.8E2,' m^3'/
     .       5X,'Percent error            = ',1PE15.8E2,' %')
 3280 FORMAT(3X,'Grid'/
     .       5X,'Spatial change   [VOLSG] = ',1PE15.8E2,' m^3'/
     .       5X,'Temporal change  [VOLTG] = ',1PE15.8E2,' m^3'/
     .       5X,'Volume error             = ',1PE15.8E2,' m^3'/
     .       5X,'Percent error            = ',1PE15.8E2,' %')
 3290 FORMAT(/1X,'Mass Balance'/
     .       '+',4('_'),1X,7('_')/)
 3300 FORMAT(3X,'Branch',I3)
 3310 FORMAT(5X,A16/
     .       7X,'Total mass   [CMBR] = ',1PE15.8E2,' kg'/
     .       7X,'Mass error          = ',1PE15.8E2,' kg'/
     .       7X,'Percent error       = ',1PE15.8E2,' %')
 3320 FORMAT(/1X,A8/
     .       '+',8('_')//
     .       3X,'Surface layer [KT] =',I8/
     .       3X,'Elevation   [ELKT] =',F8.2,' m'//
     .       3X,'Current upstream segment  [IU]'/
     .       '+',2X,7('_'),1X,8('_'),1X,7('_')//
     .       (5X,'Branch',I3,' =',I3))

***** Run time information FORMATs

 4000 FORMAT(/1X,20('*'),' Add layer',I3,' at Julian day =',F9.3,
     .         1X,20('*'))
 4010 FORMAT(/1X,11('*'),' Add from segments ',I3,' up to ',I3,
     .         ' at Julian day =',F9.3,1X,10('*'))
 4020 FORMAT(/1X,18('*'),' Subtract layer',I3,' at Julian day =',F9.3,
     .         1X,17('*'))
 4030 FORMAT(/1X,8('*'),' Subtract from segments ',I3,' up to ',I3,
     .         ' at Julian day =',F9.3,1X,8('*'))
 4040 FORMAT(//18('*'),2X,A40,2X,18('*')/'1'/'1'/'1')

***** Time series FORMAT's

 5000 FORMAT(A72)
 5010 FORMAT(30I3)
 5020 FORMAT(4(A16,1X,A6))
 5030 FORMAT(F10.3,2F10.2,(6(1PE10.3E2):/30X))

***** Run time warning FORMATs

 6000 FORMAT('Computational warning at Julian day = ',F9.2/
     .       '  timestep        = ',F4.0,' sec')
 6010 FORMAT('Computational warning at Julian day = ',F9.2/
     .       '  spatial change  =',1PE15.8E2,' m^3'/
     .       '  temporal change =',1PE15.8E2,' m^3'/
     .       '  volume error    =',1PE15.8E2,' m^3')
 6020 FORMAT('**SEVERE** computational warning at Julian day = ',F9.2/
     .       '  at segment ',I3/
     .       '    water surface elevation [Z1] =',F10.3,' m'/
     .       '    layer thickness              =',F10.3,' m')

***** Run time error FORMATs

 7000 FORMAT('Fatal error'/
     .       '  Insufficient segments in branch',I3/
     .       '  Julian day = ',F9.2,'  water surface layer =',I3)
 7010 FORMAT('Fatal error'/
     .       '  Negative surface layer height'/
     .       '  Julian day = ',F9.2,'  segment =',I3,'  height = ',
     .       F8.2,' m')

***** Contour plot FORMATs

 8000 FORMAT(13(1PE10.3E2))

      STOP
      END

************************************************************************
**                F U N C T I O N   D E N S I T Y                     **
************************************************************************

      FUNCTION DENSITY (T,SS,DS)
      LOGICAL          FRESH_WATER, SALT_WATER, SUSP_SOLIDS
      double precision t
      COMMON  /DNSPHC/ FRESH_WATER, SALT_WATER, SUSP_SOLIDS

      DENSITY = ((((6.536332E-9*T-1.120083E-6)*T+1.001685E-4)*T
     .          -9.09529E-3)*T+6.793952E-2)*T+0.842594
      IF (SUSP_SOLIDS) DENSITY = DENSITY+6.2E-4*SS
      IF (FRESH_WATER) DENSITY = DENSITY+DS*((4.99E-8*T-3.87E-6)*T
     .                           +8.221E-4)
      IF (SALT_WATER)  DENSITY = DENSITY+DS*((((5.3875E-9*T-8.2467E-7)
     .                           *T+7.6438E-5)*T-4.0899E-3)*T+0.824493)
     .                           +((-1.6546E-6*T+1.0227E-4)
     .                           *T-5.72466E-3)*DS**1.5+4.8314E-4*DS*DS
      DENSITY = DENSITY+999.0
      END

************************************************************************
**        S U B R O U T I N E    G R E G O R I A N   D A T E          **
************************************************************************

      SUBROUTINE GREGORIAN_DATE (YEAR)

***** Variable declarations

      INTEGER   YEAR, GDAY
      LOGICAL   LEAP_YEAR
      CHARACTER MONTH*9

***** Common declarations

      COMMON /GDAYC1/ GDAY, DAYM, JDAYG, LEAP_YEAR
      COMMON /GDAYC2/ MONTH

***** Determine if new year (regular or leap) and increment year

      IF (.NOT.LEAP_YEAR.AND.JDAYG.EQ.365) THEN
        JDAYG     = JDAYG-365
        YEAR      = YEAR+1
        LEAP_YEAR = MOD(YEAR,4).EQ.0
      ELSE IF (JDAYG.EQ.366) THEN
        JDAYG     = JDAYG-366
        YEAR      = YEAR+1
        LEAP_YEAR = .FALSE.
      END IF

***** Determine month and day of year

      IF (LEAP_YEAR) THEN
        IF (JDAYG.GE.0.AND.JDAYG.LT.31) THEN
          GDAY  = JDAYG+1
          DAYM  = 31.0
          MONTH = '  January'
        ELSE IF (JDAYG.GE.31.AND.JDAYG.LT.60) THEN
          GDAY  = JDAYG-30
          DAYM  = 29.0
          MONTH = ' February'
        ELSE IF (JDAYG.GE.60.AND.JDAYG.LT.91) THEN
          GDAY  = JDAYG-59
          DAYM  = 31.0
          MONTH = '    March'
        ELSE IF (JDAYG.GE.91.AND.JDAYG.LT.121) THEN
          GDAY  = JDAYG-90
          DAYM  = 30.0
          MONTH = '    April'
        ELSE IF (JDAYG.GE.121.AND.JDAYG.LT.152) THEN
          GDAY  = JDAYG-120
          DAYM  = 31.0
          MONTH = '      May'
        ELSE IF (JDAYG.GE.152.AND.JDAYG.LT.182) THEN
          GDAY  = JDAYG-151
          DAYM  = 30.0
          MONTH = '     June'
        ELSE IF (JDAYG.GE.182.AND.JDAYG.LT.213) THEN
          GDAY  = JDAYG-181
          DAYM  = 31.0
          MONTH = '     July'
        ELSE IF (JDAYG.GE.213.AND.JDAYG.LT.244) THEN
          GDAY  = JDAYG-212
          DAYM  = 31.0
          MONTH = '   August'
        ELSE IF (JDAYG.GE.244.AND.JDAYG.LT.274) THEN
          GDAY  = JDAYG-243
          DAYM  = 30.0
          MONTH = 'September'
        ELSE IF (JDAYG.GE.274.AND.JDAYG.LT.305) THEN
          GDAY  = JDAYG-273
          DAYM  = 31.0
          MONTH = '  October'
        ELSE IF (JDAYG.GE.305.AND.JDAYG.LT.335) THEN
          GDAY  = JDAYG-304
          DAYM  = 30.0
          MONTH = ' November'
        ELSE IF (JDAYG.GE.335.AND.JDAYG.LT.366) THEN
          GDAY  = JDAYG-334
          DAYM  = 31.0
          MONTH = ' December'
        END IF
      ELSE
        IF (JDAYG.GE.0.AND.JDAYG.LT.31) THEN
          GDAY  = JDAYG+1
          DAYM  = 31.0
          MONTH = '  January'
        ELSE IF (JDAYG.GE.31.AND.JDAYG.LT.59) THEN
          GDAY  = JDAYG-30
          DAYM  = 29.0
          MONTH = ' February'
        ELSE IF (JDAYG.GE.59.AND.JDAYG.LT.90) THEN
          GDAY  = JDAYG-58
          DAYM  = 31.0
          MONTH = '    March'
        ELSE IF (JDAYG.GE.90.AND.JDAYG.LT.120) THEN
          GDAY  = JDAYG-89
          DAYM  = 30.0
          MONTH = '    April'
        ELSE IF (JDAYG.GE.120.AND.JDAYG.LT.151) THEN
          GDAY  = JDAYG-119
          DAYM  = 31.0
          MONTH = '      May'
        ELSE IF (JDAYG.GE.151.AND.JDAYG.LT.181) THEN
          GDAY  = JDAYG-150
          DAYM  = 30.0
          MONTH = '     June'
        ELSE IF (JDAYG.GE.181.AND.JDAYG.LT.212) THEN
          GDAY  = JDAYG-180
          DAYM  = 31.0
          MONTH = '     July'
        ELSE IF (JDAYG.GE.212.AND.JDAYG.LT.243) THEN
          GDAY  = JDAYG-211
          DAYM  = 31.0
          MONTH = '   August'
        ELSE IF (JDAYG.GE.243.AND.JDAYG.LT.273) THEN
          GDAY  = JDAYG-242
          DAYM  = 30.0
          MONTH = 'September'
        ELSE IF (JDAYG.GE.273.AND.JDAYG.LT.304) THEN
          GDAY  = JDAYG-272
          DAYM  = 31.0
          MONTH = '  October'
        ELSE IF (JDAYG.GE.304.AND.JDAYG.LT.334) THEN
          GDAY  = JDAYG-303
          DAYM  = 30.0
          MONTH = ' November'
        ELSE IF (JDAYG.GE.334.AND.JDAYG.LT.365) THEN
          GDAY  = JDAYG-333
          DAYM  = 31.0
          MONTH = ' December'
        END IF
      END IF
      RETURN
      END

************************************************************************
**      S U B R O U T I N E   T I M E   V A R Y I N G   D A T A       **
************************************************************************

      SUBROUTINE TIME_VARYING_DATA (JDAY,NXTVD)
      INCLUDE   'w2.inc'
      SAVE
      INTEGER  OTQ, TRQ, WDQ, TRT, TRC, DTQ, DTT, DTC, PRE, PRT, PRC,
     .         UHE, UHT, UHC, DHE, DHT, DHC
      INTEGER  CN,     INCN,   TRCN,   DTCN,   PRCN,   UHCN,   DHCN,
     .         UNIT, wscdp
      REAL     NXMET1, NXQOT1, NXQIN1, NXTIN1, NXCIN1, NXQWD1, NXQTR1,
     .         NXTTR1, NXCTR1, NXQDT1, NXTDT1, NXCDT1, NXPR1,  NXTPR1,
     .         NXCPR1, NXEUH1, NXTUH1, NXCUH1, NXEDH1, NXTDH1, NXCDH1,
     .         NXTVD,  JDAY
      REAL     NXMET2, NXQOT2, NXQIN2, NXTIN2, NXCIN2, NXQWD2, NXQTR2,
     .         NXTTR2, NXCTR2, NXQDT2, NXTDT2, NXCDT2, NXPR2,  NXTPR2,
     .         NXCPR2, NXEUH2, NXTUH2, NXCUH2, NXEDH2, NXTDH2, NXCDH2
      LOGICAL  OPEN_FILES,    CONSTITUENTS, WITHDRAWALS, TRIBUTARIES,
     .         PRECIPITATION, DIST_TRIBS,   UH_EXTERNAL, DH_EXTERNAL,
     .         NO_WIND,       NO_INFLOW,    NO_OUTFLOW,  NO_HEAT,
     .         UP_FLOW,       DN_FLOW,      POINT_SINK,
     .         SEL_WITHDRAWAL,TERM_BY_TERM
      CHARACTER*72 METFN, QINFN, TINFN, CINFN, QOTFN, QWDFN, QTRFN,
     .             TTRFN, CTRFN, QDTFN, TDTFN, CDTFN, PREFN, TPRFN,
     .             CPRFN, EUHFN, TUHFN, CUHFN, EDHFN, TDHFN, CDHFN
      DIMENSION TRQ(NTP),   TRT(NTP),   TRC(NTP)
      DIMENSION INQ(NBP),   INT(NBP),   INC(NBP),   DTQ(NBP),
     .          DTT(NBP),   DTC(NBP),   PRE(NBP),   PRT(NBP),
     .          PRC(NBP),   UHE(NBP),   UHT(NBP),   UHC(NBP),
     .          DHE(NBP),   DHT(NBP),   DHC(NBP),   OTQ(NBP)

***** Common declarations

      COMMON /GLBLCC/ PALT,   ALGDET, O2LIM,  WIND, wscdp, wsc(ndp)
      COMMON /TVDMTC/ TAIR,   TDEW,   CLOUD,   PHI,   ET,   CSHE,
     .                 SRO,   LAT,    LONG,    banksh(imp)
      COMMON /INTN1C/ TAIRNX, TDEWNX, CLOUDNX, PHINX, ETNX, CSHENX,
     .                SRONX,  WINDNX
      COMMON /INTO1C/ TAIRO,  TDEWO,  CLOUDO,  PHIO,  ETO,  CSHEO,
     .                SROO,   WINDO
      COMMON /GLOBLC/ JB,     JC,     IU,      ID,    KT,   DLT,
     .                KB(IMP)
      COMMON /INTN2C/ CTRNX(NCP,NTP),     CINNX(NCP,NBP),
     .                QOUTNX(KMP,NBP),    CDTRNX(NCP,NBP),
     .                CPRNX(NCP,NBP),     TUHNX(KMP,NBP),
     .                TDHNX(KMP,NBP),     QSTRNX(NSP,NBP),
     .                CUHNX(KMP,NCP,NBP), CDHNX(KMP,NCP,NBP),
     .                QDTRNX(NBP),        TDTRNX(NBP),
     .                PRNX(NBP),          TPRNX(NBP),
     .                ELUHNX(NBP),        ELDHNX(NBP),
     .                QWDNX(NWP),         QTRNX(NTP),
     .                TTRNX(NTP),         QINNX(NBP),
     .                TINNX(NBP)
      COMMON /INTO2C/ CTRO(NCP,NTP),      CINO(NCP,NBP),
     .                QOUTO(KMP,NBP),     CDTRO(NCP,NBP),
     .                TUHO(KMP,NBP),      TDHO(KMP,NBP),
     .                QSTRO(NSP,NBP),     CUHO(KMP,NCP,NBP),
     .                CDHO(KMP,NCP,NBP),  QDTRO(NBP),
     .                TDTRO(NBP),         ELUHO(NBP),
     .                ELDHO(NBP),         QWDO(NWP),
     .                QTRO(NTP),          TTRO(NTP),
     .                QINO(NBP),          TINO(NBP)
      COMMON /INTN3C/ NXQTR1(NTP),     NXTTR1(NTP),     NXCTR1(NTP),
     .                NXQIN1(NBP),     NXTIN1(NBP),     NXCIN1(NBP),
     .                NXQDT1(NBP),     NXTDT1(NBP),     NXCDT1(NBP),
     .                NXPR1(NBP),      NXTPR1(NBP),     NXCPR1(NBP),
     .                NXEUH1(NBP),     NXTUH1(NBP),     NXCUH1(NBP),
     .                NXEDH1(NBP),     NXTDH1(NBP),     NXCDH1(NBP),
     .                NXQOT1(NBP),     NXMET1,          NXQWD1
      COMMON /INTN4C/ NXQTR2(NTP),     NXTTR2(NTP),     NXCTR2(NTP),
     .                NXQIN2(NBP),     NXTIN2(NBP),     NXCIN2(NBP),
     .                NXQDT2(NBP),     NXTDT2(NBP),     NXCDT2(NBP),
     .                NXPR2(NBP),      NXTPR2(NBP),     NXCPR2(NBP),
     .                NXEUH2(NBP),     NXTUH2(NBP),     NXCUH2(NBP),
     .                NXEDH2(NBP),     NXTDH2(NBP),     NXCDH2(NBP),
     .                NXQOT2(NBP),     NXMET2,          NXQWD2
      COMMON /TVDFNC/ METFN,           QWDFN,           QOTFN(NBP),
     .                QINFN(NBP),      TINFN(NBP),      CINFN(NBP),
     .                QTRFN(NTP),      TTRFN(NTP),      CTRFN(NTP),
     .                QDTFN(NBP),      TDTFN(NBP),      CDTFN(NBP),
     .                PREFN(NBP),      TPRFN(NBP),      CPRFN(NBP),
     .                EUHFN(NBP),      TUHFN(NBP),      CUHFN(NBP),
     .                EDHFN(NBP),      TDHFN(NBP),      CDHFN(NBP)
      COMMON /TVDQC/  QIN(NBP),        QTR(NTP),        QDTR(NBP),
     .                PR(NBP),         ELUH(NBP),       ELDH(NBP),
     .                QOUT(KMP,NBP),   QWD(NWP)
      COMMON /TVDTC/  TIN(NBP),        TTR(NTP),        TDTR(NBP),
     .                TPR(NBP),        TUH(KMP,NBP),    TDH(KMP,NBP)
      COMMON /TVDCC1/ CIN(NCP,NBP),    CTR(NCP,NTP),    CDTR(NCP,NBP),
     .                CPR(NCP,NBP),    CUH(KMC,NCP,NBP),CDH(KMC,NCP,NBP)
      COMMON /TVDCC2/ INCN(NCP),       TRCN(NCP),       DTCN(NCP),
     .                PRCN(NCP),       UHCN(NCP),       DHCN(NCP)
      COMMON /TVDCC3/ NACIN, NACTR,    NACPR, NACDT
      COMMON /TVDDSC/ NO_WIND,         NO_INFLOW,       NO_OUTFLOW,
     .                NO_HEAT
      COMMON /TVDLC1/ PRECIPITATION,   WITHDRAWALS,     TRIBUTARIES,
     .                DIST_TRIBS(NBP)
      COMMON /TVDLC3/ OPEN_FILES,      TERM_BY_TERM
      COMMON /GRTVDC/ CONSTITUENTS,    CN(NCP),         NAC
      COMMON /SELWC/  NSTR(NBP),       QSTR(NSP,NBP),   ESTR(NSP,NBP),
     .                WSTR(NSP,NBP),   KBSW(NBP),       KTOPSW(NBP),
     .                KBOTSW(NBP),     NOUT(NBP),       KOUT(KMP,NBP)
      COMMON /TVDLC2/ UP_FLOW(NBP),        DN_FLOW(NBP),
     .                UH_INTERNAL(NBP),    UH_EXTERNAL(NBP),
     .                DH_INTERNAL(NBP),    DH_EXTERNAL(NBP)
      COMMON /TVDSWC/ SEL_WITHDRAWAL(NBP), POINT_SINK(NSP,NBP)
      common /pumpreg/qpstor
c      DATA UNIT /30/
      DATA UNIT /13/

***** Open input files

      IF (OPEN_FILES) THEN
        iflag1=0
        MET  = UNIT
        UNIT = UNIT+1
        OPEN (MET,FILE=METFN,STATUS='OLD')
        READ (MET,1000)
        IF (WITHDRAWALS) THEN
          WDQ  = UNIT
          UNIT = UNIT+1
          OPEN (WDQ,FILE=QWDFN,STATUS='OLD')
          READ (WDQ,1000)
        END IF
        IF (TRIBUTARIES) THEN
c          DO JT=1,NTP
            do i=1,4
            if(i.eq.1)jt=8
            if(i.eq.2)jt=12
            if(i.eq.3)jt=17
            if(i.eq.4)jt=26
            TRQ(jt) = UNIT
            UNIT    = UNIT+1
            OPEN (TRQ(JT),FILE=QTRFN(JT),STATUS='OLD')
            READ (TRQ(JT),1000)
C            READ (TRT(JT),1000)
            IF (CONSTITUENTS) THEN
              TRC(JT) = UNIT
              UNIT    = UNIT+1
              OPEN (TRC(JT),FILE=CTRFN(JT),STATUS='OLD')
              READ (TRC(JT),1000)
            END IF
            end do
c
          DO JT=31,NTP
            TRQ(JT) = UNIT
            UNIT    = UNIT+1
c            TRT(JT) = UNIT
c            UNIT    = UNIT+1
            OPEN (TRQ(JT),FILE=QTRFN(JT),STATUS='OLD')
c            OPEN (TRT(JT),FILE=TTRFN(JT),STATUS='OLD')
            READ (TRQ(JT),1000)
c            READ (TRT(JT),1000)
            IF (CONSTITUENTS) THEN
              TRC(JT) = UNIT
              UNIT    = UNIT+1
              OPEN (TRC(JT),FILE=CTRFN(JT),STATUS='OLD')
              READ (TRC(JT),1000)
            END IF
          END DO
        END IF
          DO JB=1,NBP
c          IF (UP_FLOW(JB)) THEN
c            INQ(JB) = UNIT
c            UNIT    = UNIT+1
c            INT(JB) = UNIT
c            UNIT    = UNIT+1
c            OPEN (INQ(JB),FILE=QINFN(JB),STATUS='OLD')
c            OPEN (INT(JB),FILE=TINFN(JB),STATUS='OLD')
c            READ (INQ(JB),1000)
c            READ (INT(JB),1000)
c            IF (CONSTITUENTS) THEN
c              INC(JB) = UNIT
c              UNIT    = UNIT+1
c              OPEN (INC(JB),FILE=CINFN(JB),STATUS='OLD')
c              READ (INC(JB),1000)
c            END IF
c          END IF
          IF (DN_FLOW(JB).and.jb.eq.33) THEN
            OTQ(JB) = UNIT
            UNIT    = UNIT+1
            OPEN (OTQ(JB),FILE=QOTFN(JB),STATUS='OLD')
            READ (OTQ(JB),1000)
          END IF
          IF (PRECIPITATION) THEN
            PRE(JB) = UNIT
            UNIT    = UNIT+1
            PRT(JB) = UNIT
            UNIT    = UNIT+1
            OPEN (PRE(JB),FILE=PREFN(JB),STATUS='OLD')
            OPEN (PRT(JB),FILE=TPRFN(JB),STATUS='OLD')
            READ (PRE(JB),1000)
            READ (PRT(JB),1000)
            IF (CONSTITUENTS) THEN
              PRC(JB) = UNIT
              UNIT    = UNIT+1
              OPEN (PRC(JB),FILE=CPRFN(JB),STATUS='OLD')
              READ (PRC(JB),1000)
            END IF
          END IF
c this file is only opened for one gw input to the system all other inputs are
c the same...read in only for branch 11 and no other branches
          IF (DIST_TRIBS(JB).and.jb.eq.11) THEN
c             DTQ(JB) = UNIT
c            UNIT    = UNIT+1
             DTT(JB) = UNIT
             UNIT    = UNIT+1
c            OPEN (DTQ(JB),FILE=QDTFN(JB),STATUS='OLD')
             OPEN (DTT(JB),FILE=TDTFN(JB),STATUS='OLD')
c            READ (DTQ(JB),1000)
             READ (DTT(JB),1000)
            IF (CONSTITUENTS) THEN
              DTC(JB) = UNIT
              UNIT    = UNIT+1
              OPEN (DTC(JB),FILE=CDTFN(JB),STATUS='OLD')
              READ (DTC(JB),1000)
            END IF
          END IF
          IF (UH_EXTERNAL(JB)) THEN
            UHE(JB) = UNIT
            UNIT    = UNIT+1
            UHT(JB) = UNIT
            UNIT    = UNIT+1
            OPEN (UHE(JB),FILE=EUHFN(JB),STATUS='OLD')
            OPEN (UHT(JB),FILE=TUHFN(JB),STATUS='OLD')
            READ (UHE(JB),1000)
            READ (UHT(JB),1000)
            IF (CONSTITUENTS) THEN
              UHC(JB) = UNIT
              UNIT    = UNIT+1
              OPEN (UHC(JB),FILE=CUHFN(JB),STATUS='OLD')
              READ (UHC(JB),1000)
            END IF
          END IF
          IF (DH_EXTERNAL(JB)) THEN
            DHE(JB) = UNIT
            UNIT    = UNIT+1
            DHT(JB) = UNIT
            UNIT    = UNIT+1
            OPEN (DHE(JB),FILE=EDHFN(JB),STATUS='OLD')
            OPEN (DHT(JB),FILE=TDHFN(JB),STATUS='OLD')
            READ (DHE(JB),1000)
            READ (DHT(JB),1000)
            IF (CONSTITUENTS) THEN
              DHC(JB) = UNIT
              UNIT    = UNIT+1
              OPEN (DHC(JB),FILE=CDHFN(JB),STATUS='OLD')
              READ (DHC(JB),1000)
            END IF
          END IF
        END DO
        OPEN_FILES = .FALSE.
      END IF

c open downstream head bc files if ready to

       if(dh_external(8))then
         if(iflag1.eq.0)then
         iflag1=1
            JB=8
            DHE(JB) = UNIT
            UNIT    = UNIT+1
            DHT(JB) = UNIT
            UNIT    = UNIT+1
            OPEN (DHE(JB),FILE=EDHFN(JB),STATUS='OLD')
            OPEN (DHT(JB),FILE=TDHFN(JB),STATUS='OLD')
            READ (DHE(JB),1000)
            READ (DHT(JB),1000)
            IF (CONSTITUENTS) THEN
              DHC(JB) = UNIT
              UNIT    = UNIT+1
              OPEN (DHC(JB),FILE=CDHFN(JB),STATUS='OLD')
              READ (DHC(JB),1000)
            END IF
         end if
       end if

***** Meteorlogical data

      IF (JDAY.GE.NXMET1) THEN
        DO WHILE (JDAY.GE.NXMET1)
          NXMET2 = NXMET1
          TDEW   = TDEWNX
          WIND   = WINDNX
          PHI    = PHINX
          TDEWO  = TDEWNX
          WINDO  = WINDNX
          PHIO   = PHINX
          IF (TERM_BY_TERM) THEN
            TAIR   = TAIRNX
            CLOUD  = CLOUDNX
            TAIRO  = TAIRNX
            CLOUDO = CLOUDNX
            READ (MET,1010) NXMET1,TAIRNX,TDEWNX,WINDNX,PHINX,CLOUDNX
          ELSE
            ET    = ETNX
            CSHE  = CSHENX
            SRO   = SRONX
            ETO   = ETNX
            CSHEO = CSHENX
            SROO  = SRONX
            READ (MET,1010) NXMET1,ETNX,TDEWNX,WINDNX,PHINX,CSHENX,SRONX
          END IF
          windnx=windnx*wsc(wscdp)
        END DO
      END IF
      NXTVD = MIN(NXTVD,NXMET1)

***** Withdrawals

      IF (WITHDRAWALS) THEN
        IF (JDAY.GE.NXQWD1) THEN
          DO WHILE (JDAY.GE.NXQWD1)
            NXQWD2 = NXQWD1
c            DO JW=1,NWP
             jw=29
              QWD(JW)  = QWDNX(JW)
              if(jw.eq.29)qpstor=qwd(29)               
              QWDO(JW) = QWDNX(JW)
c            END DO
c            READ (WDQ,1020) NXQWD1,(QWDNX(JW),JW=1,NWP)
            READ (WDQ,1020) NXQWD1,QWDNX(JW)
          END DO
        END IF
        NXTVD = MIN(NXTVD,NXQWD1)
      END IF

***** Tributaries

      IF (TRIBUTARIES) THEN
        do i=1,4
            if(i.eq.1)jt=8
            if(i.eq.2)jt=12
            if(i.eq.3)jt=17
            if(i.eq.4)jt=26
          IF (JDAY.GE.NXQTR1(JT)) THEN
            DO WHILE (JDAY.GE.NXQTR1(JT))
              NXQTR2(JT) = NXQTR1(JT)
              QTR(JT)    = QTRNX(JT)
              QTRO(JT)   = QTRNX(JT)
              TTR(JT)    = TTRNX(JT)
              TTRO(JT)   = TTRNX(JT)
              READ (TRQ(JT),1020) NXQTR1(JT),QTRNX(JT),ttrnx(jt)
            END DO
          END IF
          NXTVD = MIN(NXTVD,NXQTR1(JT))

********* Inflow constituent concentrations

          IF (CONSTITUENTS) THEN
            IF (JDAY.GE.NXCTR1(JT)) THEN
              DO WHILE (JDAY.GE.NXCTR1(JT))
                NXCTR2(JT) = NXCTR1(JT)
                DO JC=1,NACTR
                  CTR(TRCN(JC),JT)  = CTRNX(TRCN(JC),JT)
                  CTRO(TRCN(JC),JT) = CTRNX(TRCN(JC),JT)
                END DO
c                READ (TRC(JT),1020) NXCTR1(JT),(CTRNX(TRCN(JC),JT),
                READ (TRC(JT),1021) NXCTR1(JT),(CTRNX(TRCN(JC),JT),
     .                              JC=1,NACTR)
              END DO
            END IF
            NXTVD = MIN(NXTVD,NXCTR1(JT))
          END IF
        end do
        DO JT=31,NTP

********* Inflow

          IF (JDAY.GE.NXQTR1(JT)) THEN
            DO WHILE (JDAY.GE.NXQTR1(JT))
              NXQTR2(JT) = NXQTR1(JT)
              QTR(JT)    = QTRNX(JT)
              QTRO(JT)   = QTRNX(JT)
              TTR(JT)    = TTRNX(JT)
              TTRO(JT)   = TTRNX(JT)
              READ (TRQ(JT),1020) NXQTR1(JT),QTRNX(JT),ttrnx(jt)
            END DO
          END IF
          NXTVD = MIN(NXTVD,NXQTR1(JT))

********* Inflow temperatures

c          IF (JDAY.GE.NXTTR1(JT)) THEN
c            DO WHILE (JDAY.GE.NXTTR1(JT))
c              NXTTR2(JT) = NXTTR1(JT)
c              TTR(JT)    = TTRNX(JT)
c              TTRO(JT)   = TTRNX(JT)
c              READ (TRT(JT),1020) NXTTR1(JT),TTRNX(JT)
c            END DO
c          END IF
c          NXTVD = MIN(NXTVD,NXTTR1(JT))

********* Inflow constituent concentrations

          IF (CONSTITUENTS) THEN
            IF (JDAY.GE.NXCTR1(JT)) THEN
              DO WHILE (JDAY.GE.NXCTR1(JT))
                NXCTR2(JT) = NXCTR1(JT)
                DO JC=1,NACTR
                  CTR(TRCN(JC),JT)  = CTRNX(TRCN(JC),JT)
                  CTRO(TRCN(JC),JT) = CTRNX(TRCN(JC),JT)
                END DO
c                READ (TRC(JT),1020) NXCTR1(JT),(CTRNX(TRCN(JC),JT),
                READ (TRC(JT),1021) NXCTR1(JT),(CTRNX(TRCN(JC),JT),
     .                              JC=1,NACTR)
              END DO
            END IF
            NXTVD = MIN(NXTVD,NXCTR1(JT))
          END IF
        END DO
      END IF

***** Branch related inputs

      DO JB=1,NBP

******* Inflow

c        IF (UP_FLOW(JB)) THEN
c          IF (JDAY.GE.NXQIN1(JB)) THEN
c            DO WHILE (JDAY.GE.NXQIN1(JB))
c              NXQIN2(JB) = NXQIN1(JB)
c              QIN(JB)    = QINNX(JB)
c              QINO(JB)   = QINNX(JB)
c              READ (INQ(JB),1020) NXQIN1(JB),QINNX(JB)
c            END DO
c          END IF
c          NXTVD = MIN(NXTVD,NXQIN1(JB))
c
********* Inflow temperature

c          IF (JDAY.GE.NXTIN1(JB)) THEN
c            DO WHILE (JDAY.GE.NXTIN1(JB))
c              NXTIN2(JB) = NXTIN1(JB)
c              TIN(JB)    = TINNX(JB)
c              TINO(JB)   = TINNX(JB)
c              READ (INT(JB),1020) NXTIN1(JB),TINNX(JB)
c            END DO
c          END IF
c          NXTVD = MIN(NXTVD,NXTIN1(JB))
c
********* Inflow constituent concentrations

c          IF (CONSTITUENTS) THEN
c            IF (JDAY.GE.NXCIN1(JB)) THEN
c              DO WHILE (JDAY.GE.NXCIN1(JB))
c                NXCIN2(JB) = NXCIN1(JB)
c                DO JC=1,NACIN
c                  CIN(INCN(JC),JB)  = CINNX(INCN(JC),JB)
c                  CINO(INCN(JC),JB) = CINNX(INCN(JC),JB)
c                END DO
cc                READ (INC(JB),1030) NXCIN1(JB),(CINNX(INCN(JC),JB),
c                READ (INC(JB),1021) NXCIN1(JB),(CINNX(INCN(JC),JB),
c     .                              JC=1,NACIN)
c              END DO
c            END IF
c            NXTVD = MIN(NXTVD,NXCIN1(JB))
c          END IF
c        END IF

******* Outflow

        IF (DN_FLOW(JB).and.jb.eq.33) THEN
          IF (JDAY.GE.NXQOT1(JB)) THEN
            DO WHILE (JDAY.GE.NXQOT1(JB))
              NXQOT2(JB) = NXQOT1(JB)
              IF (SEL_WITHDRAWAL(JB)) THEN
                DO JS=1,NSTR(JB)
                  QSTR(JS,JB)  = QSTRNX(JS,JB)
                  QSTRO(JS,JB) = QSTRNX(JS,JB)
                END DO
                READ (OTQ(JB),1020) NXQOT1(JB),(QSTRNX(JS,JB),
     .                              JS=1,NSTR(JB))
              ELSE
                DO JO=1,NOUT(JB)
                  QOUT(JO,JB)  = QOUTNX(JO,JB)
                  QOUTO(JO,JB) = QOUTNX(JO,JB)
                END DO
                READ (OTQ(JB),1020) NXQOT1(JB),(QOUTNX(JO,JB),
     .                              JO=1,NOUT(JB))
              END IF
            END DO
          END IF
          NXTVD = MIN(NXTVD,NXQOT1(JB))
        END IF

******* Distributed tributaries

        IF (DIST_TRIBS(JB).and.jb.eq.11) THEN

********* Inflow

c          IF (JDAY.GE.NXQDT1(JB)) THEN
c            DO WHILE (JDAY.GE.NXQDT1(JB))
c              NXQDT2(JB) = NXQDT1(JB)
c              QDTR(JB)   = QDTRNX(JB)
c              QDTRO(JB)  = QDTRNX(JB)
c              READ (DTQ(JB),1020) NXQDT1(JB),QDTRNX(JB)
c            END DO
c          END IF
c          NXTVD = MIN(NXTVD,NXQDT1(JB))

********* Temperature

          IF (JDAY.GE.NXTDT1(JB)) THEN
            DO WHILE (JDAY.GE.NXTDT1(JB))
              NXTDT2(JB) = NXTDT1(JB)
              TDTR(JB)   = TDTRNX(JB)
c set other distributed temps for other jb's
              do ii=1,nbp
               IF (DIST_TRIBS(ii))then
                TDTR(ii)   = TDTRNX(JB)
               end if
              end do
c              TDTR(16)   = TDTRNX(JB)
c              TDTR(15)   = TDTRNX(JB)
c              TDTR(25)   = TDTRNX(JB)
c              TDTR(27)   = TDTRNX(JB)
c              TDTR(30)   = TDTRNX(JB)
c              TDTR(31)   = TDTRNX(JB)
c
              TDTRO(JB)  = TDTRNX(JB)
              READ (DTT(JB),1020) NXTDT1(JB),TDTRNX(JB)
            END DO
          END IF
          NXTVD = MIN(NXTVD,NXTDT1(JB))

********* Constituent concentrations

          IF (CONSTITUENTS) THEN
            IF (JDAY.GE.NXCDT1(JB)) THEN
              DO WHILE (JDAY.GE.NXCDT1(JB))
                NXCDT2(JB) = NXCDT1(JB)
                DO JC=1,NACDT
                  CDTR(DTCN(JC),JB)  = CDTRNX(DTCN(JC),JB)
                   do ii=1,nbp
                    IF (DIST_TRIBS(ii))then
                     CDTR(DTCN(JC),ii)  = CDTRNX(DTCN(JC),JB)
                    end if
                   end do
c set other gw dst loadings
c                  CDTR(DTCN(JC),12)  = CDTRNX(DTCN(JC),JB)
c                  CDTR(DTCN(JC),16)  = CDTRNX(DTCN(JC),JB)
c                  CDTR(DTCN(JC),15)  = CDTRNX(DTCN(JC),JB)
c                  CDTR(DTCN(JC),25)  = CDTRNX(DTCN(JC),JB)
c                  CDTR(DTCN(JC),27)  = CDTRNX(DTCN(JC),JB)
c                  CDTR(DTCN(JC),30)  = CDTRNX(DTCN(JC),JB)
c                  CDTR(DTCN(JC),31)  = CDTRNX(DTCN(JC),JB)
c
                  CDTRO(DTCN(JC),JB) = CDTRNX(DTCN(JC),JB)
                END DO
c                READ (DTC(JB),1020) NXCDT1(JB),(CDTRNX(DTCN(JC),JB),
                READ (DTC(JB),1021) NXCDT1(JB),(CDTRNX(DTCN(JC),JB),
     .                              JC=1,NACDT)
              END DO
            END IF
            NXTVD = MIN(NXTVD,NXCDT1(JB))
          END IF
        END IF

******* Precipitation

        IF (PRECIPITATION) THEN
          IF (JDAY.GE.NXPR1(JB)) THEN
           DO WHILE (JDAY.GE.NXPR1(JB))
              NXPR2(JB) = NXPR1(JB)
              PR(JB)    = PRNX(JB)
              READ (PRE(JB),1020) NXPR1(JB),PRNX(JB)
            END DO
          END IF
          NXTVD = MIN(NXTVD,NXPR1(JB))

********* Temperature

          IF (JDAY.GE.NXTPR1(JB)) THEN
            DO WHILE (JDAY.GE.NXTPR1(JB))
              NXTPR2(JB) = NXTPR1(JB)
              TPR(JB)    = TPRNX(JB)
              READ (PRT(JB),1020) NXTPR1(JB),TPRNX(JB)
            END DO
          END IF
          NXTVD = MIN(NXTVD,NXTPR1(JB))

********* Constituent concentrations

          IF (CONSTITUENTS) THEN
            IF (JDAY.GE.NXCPR1(JB)) THEN
              DO WHILE (JDAY.GE.NXCPR1(JB))
                NXCPR2(JB) = NXCPR1(JB)
                DO JC=1,NACPR
                  CPR(PRCN(JC),JB) = CPRNX(PRCN(JC),JB)
                END DO
c                READ (PRC(JB),1020) NXCPR1(JB),(CPRNX(PRCN(JC),JB),
                READ (PRC(JB),1021) NXCPR1(JB),(CPRNX(PRCN(JC),JB),
     .                              JC=1,NACPR)
              END DO
            END IF
            NXTVD = MIN(NXTVD,NXCPR1(JB))
          END IF
        END IF

******* Upstream head conditions

        IF (UH_EXTERNAL(JB)) THEN

********* Elevations

          IF (JDAY.GE.NXEUH1(JB)) THEN
            DO WHILE (JDAY.GE.NXEUH1(JB))
              NXEUH2(JB) = NXEUH1(JB)
              ELUH(JB)   = ELUHNX(JB)
              ELUHO(JB)  = ELUHNX(JB)
              READ (UHE(JB),1020) NXEUH1(JB),ELUHNX(JB)
            END DO
          END IF
          NXTVD = MIN(NXTVD,NXEUH1(JB))

********* Temperatures

          IF (JDAY.GE.NXTUH1(JB)) THEN
            DO WHILE (JDAY.GE.NXTUH1(JB))
              NXTUH2(JB) = NXTUH1(JB)
              DO K=2,KMP-1
                TUH(K,JB)  = TUHNX(K,JB)
                TUHO(K,JB) = TUHNX(K,JB)
              END DO
              READ (UHT(JB),1020) NXTUH1(JB),(TUHNX(K,JB),K=2,KMP-1)
            END DO
          END IF
          NXTVD = MIN(NXTVD,NXTUH1(JB))

********* Constituent concentrations

          IF (CONSTITUENTS) THEN
            IF (JDAY.GE.NXCUH1(JB)) THEN
              DO WHILE (JDAY.GE.NXCUH1(JB))
                NXCUH2(JB) = NXCUH1(JB)
c                DO JC=1,NAC
                 DO JC=1,nacin
                  DO K=2,KMP-1
                    CUH(K,CN(JC),JB)  = CUHNX(K,CN(JC),JB)
                    CUHO(K,CN(JC),JB) = CUHNX(K,CN(JC),JB)
                  END DO
                END DO
c                DO JC=1,NAC
                DO JC=1,nacin
                  READ (UHC(JB),1020) NXCUH1(JB),(CUHNX(K,CN(JC),JB),
     .                                K=2,KMP-1)
                END DO
              END DO
            END IF
            NXTVD = MIN(NXTVD,NXCUH1(JB))
          END IF
        END IF

******* Downstream head

        IF (DH_EXTERNAL(JB)) THEN

********* Elevation

          IF (JDAY.GE.NXEDH1(JB)) THEN
            DO WHILE (JDAY.GE.NXEDH1(JB))
              NXEDH2(JB) = NXEDH1(JB)
              ELDH(JB)   = ELDHNX(JB)
              ELDHO(JB)  = ELDHNX(JB)
              READ (DHE(JB),1020) NXEDH1(JB),ELDHNX(JB)
            END DO
          END IF
          NXTVD = MIN(NXTVD,NXEDH1(JB))

********* Temperature

          IF (JDAY.GE.NXTDH1(JB)) THEN
            DO WHILE (JDAY.GE.NXTDH1(JB))
              NXTDH2(JB) = NXTDH1(JB)
              DO K=2,KMP-1
                TDH(K,JB)  = TDHNX(K,JB)
                TDHO(K,JB) = TDHNX(K,JB)
              END DO
              READ (DHT(JB),1020) NXTDH1(JB),(TDHNX(K,JB),K=2,KMP-1)
            END DO
          END IF
          NXTVD = MIN(NXTVD,NXTDH1(JB))

********* Constituents

          IF (CONSTITUENTS) THEN
            IF (JDAY.GE.NXCDH1(JB)) THEN
              DO WHILE (JDAY.GE.NXCDH1(JB))
                NXCDH2(JB) = NXCDH1(JB)
c                DO JC=1,NAC
                DO JC=1,nacin
                  DO K=2,KMP-1
                    CDH(K,CN(JC),JB)  = CDHNX(K,CN(JC),JB)
                    CDHO(K,CN(JC),JB) = CDHNX(K,CN(JC),JB)
                  END DO
                END DO
c                DO JC=1,NAC
                DO JC=1,nacin
                  READ (DHC(JB),1020) NXCDH1(JB),(CDHNX(K,CN(JC),JB),
     .                                K=2,KMP-1)
                END DO
              END DO
            END IF
            NXTVD = MIN(NXTVD,NXCDH1(JB))
          END IF
        END IF
      END DO

***** Dead sea case

      IF (NO_WIND) THEN
        WIND   = 0.0
        WINDO  = 0.0
        WINDNX = 0.0
      END IF
      IF (NO_INFLOW) THEN
        DO JB=1,NBP
          QIN(JB)    = 0.0
          QINO(JB)   = 0.0
          QINNX(JB)  = 0.0
          QDTRO(JB)  = 0.0
          QDTRNX(JB) = 0.0
          PR(JB)     = 0.0
          PRNX(JB)   = 0.0
        END DO
c        DO JT=1,NTP
c          QTR(JT)   = 0.0
c          QTRO(JT)  = 0.0
c          QTRNX(JT) = 0.0
c        END DO
      END IF
      IF (NO_OUTFLOW) THEN
        DO JB=1,NBP
          IF (SEL_WITHDRAWAL(JB)) THEN
            DO JS=1,NSTR(JB)
              QSTR(JS,JB)   = 0.0
              QSTRO(JS,JB)  = 0.0
              QSTRNX(JS,JB) = 0.0
            END DO
          ELSE
            DO JO=1,NOUT(JB)
              QOUT(JO,JB)   = 0.0
              QOUTO(JO,JB)  = 0.0
              QOUTNX(JO,JB) = 0.0
            END DO
          ENDIF
        END DO
c        DO JW=1,NWP
c          QWD(JW)   = 0.0
c          QWDO(JW)  = 0.0
c          QWDNX(JW) = 0.0
c        END DO
      END IF
      IF (NO_HEAT) THEN
        CSHE   = 0.0
        CSHEO  = 0.0
        CSHENX = 0.0
        SRO    = 0.0
        SROO   = 0.0
        SRONX  = 0.0
      END IF

***** Input FORMATs

 1000 FORMAT(//)
 1010 FORMAT(10F8.0)
 1020 FORMAT(10F8.0/(8X,9F8.0))
 1021 format(18f7.0)
 1030 FORMAT(18F8.0)
      END

***********************************************************************
**            S U B R O U T I N E   I N T E R P O L A T E            **
***********************************************************************

      SUBROUTINE INTERPOLATE_INPUTS (JDAY)
      INCLUDE  'w2.inc'
      REAL     JDAY
      INTEGER  CN,     INCN,   TRCN,   DTCN,   PRCN,   UHCN,   DHCN,
     .         wscdp
      REAL     NXMET1, NXQOT1, NXQIN1, NXTIN1, NXCIN1, NXQWD1, NXQTR1,
     .         NXTTR1, NXCTR1, NXQDT1, NXTDT1, NXCDT1, NXPR1,  NXTPR1,
     .         NXCPR1, NXEUH1, NXTUH1, NXCUH1, NXEDH1, NXTDH1, NXCDH1
      REAL     NXMET2, NXQOT2, NXQIN2, NXTIN2, NXCIN2, NXQWD2, NXQTR2,
     .         NXTTR2, NXCTR2, NXQDT2, NXTDT2, NXCDT2, NXPR2,  NXTPR2,
     .         NXCPR2, NXEUH2, NXTUH2, NXCUH2, NXEDH2, NXTDH2, NXCDH2
      LOGICAL  OPEN_FILES,    CONSTITUENTS,    WITHDRAWALS,
     .         TRIBUTARIES,   PRECIPITATION,   DIST_TRIBS,
     .         UH_EXTERNAL,   DH_EXTERNAL,     UP_FLOW,
     .         DN_FLOW,       POINT_SINK,      SEL_WITHDRAWAL,
     .         TERM_BY_TERM
      LOGICAL  INTERP_INFLOW, INTERP_OUTFLOW,  INTERP_MET,
     .         INTERP_DTRIBS, INTERP_TRIBS,    INTERP_HEAD,
     .         INTERP_WITHDRWL
      COMMON /GLBLCC/ PALT,   ALGDET, O2LIM,   WIND, wscdp, wsc(ndp)
      COMMON /TVDMTC/ TAIR,   TDEW,   CLOUD,   PHI,   ET,   CSHE,
     .                SRO,    LAT,    LONG,    banksh(imp)
      COMMON /INTN1C/ TAIRNX, TDEWNX, CLOUDNX, PHINX, ETNX, CSHENX,
     .                SRONX,  WINDNX
      COMMON /INTO1C/ TAIRO,  TDEWO,  CLOUDO,  PHIO,  ETO,  CSHEO,
     .                SROO,   WINDO
      COMMON /INTN2C/ CTRNX(NCP,NTP),     CINNX(NCP,NBP),
     .                QOUTNX(KMP,NBP),    CDTRNX(NCP,NBP),
     .                CPRNX(NCP,NBP),     TUHNX(KMP,NBP),
     .                TDHNX(KMP,NBP),     QSTRNX(NSP,NBP),
     .                CUHNX(KMP,NCP,NBP), CDHNX(KMP,NCP,NBP),
     .                QDTRNX(NBP),        TDTRNX(NBP),
     .                PRNX(NBP),          TPRNX(NBP),
     .                ELUHNX(NBP),        ELDHNX(NBP),
     .                QWDNX(NWP),         QTRNX(NTP),
     .                TTRNX(NTP),         QINNX(NBP),
     .                TINNX(NBP)
      COMMON /INTO2C/ CTRO(NCP,NTP),      CINO(NCP,NBP),
     .                QOUTO(KMP,NBP),     CDTRO(NCP,NBP),
     .                TUHO(KMP,NBP),      TDHO(KMP,NBP),
     .                QSTRO(NSP,NBP),     CUHO(KMP,NCP,NBP),
     .                CDHO(KMP,NCP,NBP),  QDTRO(NBP),
     .                TDTRO(NBP),         ELUHO(NBP),
     .                ELDHO(NBP),         QWDO(NWP),
     .                QTRO(NTP),          TTRO(NTP),
     .                QINO(NBP),          TINO(NBP)
      COMMON /INTN3C/ NXQTR1(NTP),     NXTTR1(NTP),     NXCTR1(NTP),
     .                NXQIN1(NBP),     NXTIN1(NBP),     NXCIN1(NBP),
     .                NXQDT1(NBP),     NXTDT1(NBP),     NXCDT1(NBP),
     .                NXPR1(NBP),      NXTPR1(NBP),     NXCPR1(NBP),
     .                NXEUH1(NBP),     NXTUH1(NBP),     NXCUH1(NBP),
     .                NXEDH1(NBP),     NXTDH1(NBP),     NXCDH1(NBP),
     .                NXQOT1(NBP),     NXMET1,          NXQWD1
      COMMON /INTN4C/ NXQTR2(NTP),     NXTTR2(NTP),     NXCTR2(NTP),
     .                NXQIN2(NBP),     NXTIN2(NBP),     NXCIN2(NBP),
     .                NXQDT2(NBP),     NXTDT2(NBP),     NXCDT2(NBP),
     .                NXPR2(NBP),      NXTPR2(NBP),     NXCPR2(NBP),
     .                NXEUH2(NBP),     NXTUH2(NBP),     NXCUH2(NBP),
     .                NXEDH2(NBP),     NXTDH2(NBP),     NXCDH2(NBP),
     .                NXQOT2(NBP),     NXMET2,          NXQWD2
      COMMON /TVDQC/  QIN(NBP),        QTR(NTP),        QDTR(NBP),
     .                PR(NBP),         ELUH(NBP),       ELDH(NBP),
     .                QOUT(KMP,NBP),   QWD(NWP)
      COMMON /TVDTC/  TIN(NBP),        TTR(NTP),        TDTR(NBP),
     .                TPR(NBP),        TUH(KMP,NBP),    TDH(KMP,NBP)
      COMMON /TVDCC1/ CIN(NCP,NBP),    CTR(NCP,NTP),    CDTR(NCP,NBP),
     .                CPR(NCP,NBP),    CUH(KMC,NCP,NBP),CDH(KMC,NCP,NBP)
      COMMON /TVDCC2/ INCN(NCP),       TRCN(NCP),       DTCN(NCP),
     .                PRCN(NCP),       UHCN(NCP),       DHCN(NCP)
      COMMON /TVDCC3/ NACIN, NACTR,    NACPR, NACDT
      COMMON /TVDLC1/ PRECIPITATION,   WITHDRAWALS,     TRIBUTARIES,
     .                DIST_TRIBS(NBP)
      COMMON /TVDLC3/ OPEN_FILES,      TERM_BY_TERM
      COMMON /SELWC/  NSTR(NBP),       QSTR(NSP,NBP),   ESTR(NSP,NBP),
     .                WSTR(NSP,NBP),   KBSW(NBP),       KTOPSW(NBP),
     .                KBOTSW(NBP),     NOUT(NBP),       KOUT(KMP,NBP)
      COMMON /GRTVDC/ CONSTITUENTS,    CN(NCP),         NAC
      COMMON /INTERC/ INTERP_INFLOW,   INTERP_OUTFLOW,  INTERP_MET,
     .                INTERP_DTRIBS,   INTERP_TRIBS,    INTERP_HEAD,
     .                INTERP_WITHDRWL
      COMMON /TVDLC2/ UP_FLOW(NBP),        DN_FLOW(NBP),
     .                UH_INTERNAL(NBP),    UH_EXTERNAL(NBP),
     .                DH_INTERNAL(NBP),    DH_EXTERNAL(NBP)
      COMMON /TVDSWC/ SEL_WITHDRAWAL(NBP), POINT_SINK(NSP,NBP)

***** Meteorlogical data

      IF (INTERP_MET) THEN
        RATIO = (NXMET1-JDAY)/(NXMET1-NXMET2)
        TDEW  = (1.0-RATIO)*TDEWNX+RATIO*TDEWO
        WIND  = (1.0-RATIO)*WINDNX+RATIO*WINDO
        PHI   = (1.0-RATIO)*PHINX+RATIO*PHIO
        IF (TERM_BY_TERM) THEN
          TAIR  = (1.0-RATIO)*TAIRNX+RATIO*TAIRO
          CLOUD = (1.0-RATIO)*CLOUDNX+RATIO*CLOUDO
        ELSE
          ET   = (1.0-RATIO)*ETNX+RATIO*ETO
          CSHE = (1.0-RATIO)*CSHENX+RATIO*CSHEO
          SRO  = (1.0-RATIO)*SRONX+RATIO*SROO
        END IF
      END IF

***** Withdrawals

      IF (WITHDRAWALS) THEN
        IF (INTERP_WITHDRWL) THEN
          QRATIO = (NXQWD1-JDAY)/(NXQWD1-NXQWD2)
c          DO JW=1,NWP
           jw=29
            QWD(JW) = (1.0-QRATIO)*QWDNX(JW)+QRATIO*QWDO(JW)
c          END DO
        END IF
      END IF

***** Tributaries

      IF (TRIBUTARIES) THEN
        IF (INTERP_TRIBS) THEN
          DO JT=1,NTP
            QRATIO = (NXQTR1(JT)-JDAY)/(NXQTR1(JT)-NXQTR2(JT))
            TRATIO = (NXTTR1(JT)-JDAY)/(NXTTR1(JT)-NXTTR2(JT))
            IF (CONSTITUENTS) THEN
              CRATIO = (NXCTR1(JT)-JDAY)/(NXCTR1(JT)-NXCTR2(JT))
            END IF
            QTR(JT) = (1.0-QRATIO)*QTRNX(JT)+QRATIO*QTRO(JT)
            TTR(JT) = (1.0-TRATIO)*TTRNX(JT)+TRATIO*TTRO(JT)
            DO JC=1,NACTR
              CTR(TRCN(JC),JT) = (1.0-CRATIO)*CTRNX(TRCN(JC),JT)+CRATIO
     .                           *CTRO(TRCN(JC),JT)
            END DO
          END DO
        END IF
      END IF

***** Branch related inputs

      DO JB=1,NBP

******* Inflow

c        IF (UP_FLOW(JB)) THEN
c          IF (INTERP_INFLOW) THEN
c            QRATIO = (NXQIN1(JB)-JDAY)/(NXQIN1(JB)-NXQIN2(JB))
c            TRATIO = (NXTIN1(JB)-JDAY)/(NXTIN1(JB)-NXTIN2(JB))
c            IF (CONSTITUENTS) THEN
c              CRATIO = (NXCIN1(JB)-JDAY)/(NXCIN1(JB)-NXCIN2(JB))
c            END IF
c            QIN(JB) = (1.0-QRATIO)*QINNX(JB)+QRATIO*QINO(JB)
c            TIN(JB) = (1.0-TRATIO)*TINNX(JB)+TRATIO*TINO(JB)
c            DO JC=1,NACIN
c              CIN(INCN(JC),JB) = (1.0-CRATIO)*CINNX(INCN(JC),JB)+CRATIO
c     .                           *CINO(INCN(JC),JB)
c            END DO
c          END IF
c        END IF

******* Outflow

        IF (DN_FLOW(JB).and.jb.eq.33) THEN
          IF (INTERP_OUTFLOW) THEN
            QRATIO = (NXQOT1(JB)-JDAY)/(NXQOT1(JB)-NXQOT2(JB))
            IF (SEL_WITHDRAWAL(JB)) THEN
              DO JS=1,NSTR(JB)
                QSTR(JS,JB) = (1.0-QRATIO)*QSTRNX(JS,JB)
     .                        +QRATIO*QSTRO(JS,JB)
              END DO
            ELSE
              DO JO=1,NOUT(JB)
                QOUT(JO,JB) = (1.0-QRATIO)*QOUTNX(JO,JB)
     .                        +QRATIO*QOUTO(JO,JB)
              END DO
            END IF
          END IF
        END IF

******* Distributed tributaries

        IF (DIST_TRIBS(JB).and.jb.eq.11) THEN
          IF (INTERP_DTRIBS) THEN
            
c            QRATIO = (NXQDT1(JB)-JDAY)/(NXQDT1(JB)-NXQDT2(JB))
            TRATIO = (NXTDT1(JB)-JDAY)/(NXTDT1(JB)-NXTDT2(JB))
            IF (CONSTITUENTS) THEN
              CRATIO = (NXCDT1(JB)-JDAY)/(NXCDT1(JB)-NXCDT2(JB))
            END IF
c            QDTR(JB) = (1.0-QRATIO)*QDTRNX(JB)+QRATIO*QDTRO(JB)
             TDTR(JB) = (1.0-TRATIO)*TDTRNX(JB)+TRATIO*TDTRO(JB)
c set other branch inputs
             do ii=1,nbp
              IF (DIST_TRIBS(ii))then
                TDTR(ii) = (1.0-TRATIO)*TDTRNX(jb)+TRATIO*TDTRO(jb)
              end if
             end do
c            TDTR(12) = (1.0-TRATIO)*TDTRNX(JB)+TRATIO*TDTRO(JB)
c            TDTR(16) = (1.0-TRATIO)*TDTRNX(JB)+TRATIO*TDTRO(JB)
c            TDTR(15) = (1.0-TRATIO)*TDTRNX(JB)+TRATIO*TDTRO(JB)
c            TDTR(25) = (1.0-TRATIO)*TDTRNX(JB)+TRATIO*TDTRO(JB)
c            TDTR(27) = (1.0-TRATIO)*TDTRNX(JB)+TRATIO*TDTRO(JB)
c            TDTR(30) = (1.0-TRATIO)*TDTRNX(JB)+TRATIO*TDTRO(JB)
c            TDTR(31) = (1.0-TRATIO)*TDTRNX(JB)+TRATIO*TDTRO(JB)
c
            DO JC=1,NACDT
              CDTR(DTCN(JC),JB) = (1.0-CRATIO)*CDTRNX(DTCN(JC),JB)
     .                            +CRATIO*CDTRO(DTCN(JC),JB)
c set other branch inputs
              do ii=1,nbp
              IF (DIST_TRIBS(ii))then
               CDTR(DTCN(JC),ii) = (1.0-CRATIO)*CDTRNX(DTCN(JC),jb)
     .                 +CRATIO*CDTRO(DTCN(JC),JB) 
              end if
             end do
c              CDTR(DTCN(JC),12) = (1.0-CRATIO)*CDTRNX(DTCN(JC),JB)
c     .                            +CRATIO*CDTRO(DTCN(JC),JB)
c              CDTR(DTCN(JC),16) = (1.0-CRATIO)*CDTRNX(DTCN(JC),JB)
c     .                            +CRATIO*CDTRO(DTCN(JC),JB)
c              CDTR(DTCN(JC),15) = (1.0-CRATIO)*CDTRNX(DTCN(JC),JB)
c     .                            +CRATIO*CDTRO(DTCN(JC),JB)
c              CDTR(DTCN(JC),25) = (1.0-CRATIO)*CDTRNX(DTCN(JC),JB)
c     .                            +CRATIO*CDTRO(DTCN(JC),JB)
c              CDTR(DTCN(JC),27) = (1.0-CRATIO)*CDTRNX(DTCN(JC),JB)
c     .                            +CRATIO*CDTRO(DTCN(JC),JB)
c              CDTR(DTCN(JC),30) = (1.0-CRATIO)*CDTRNX(DTCN(JC),JB)
c     .                            +CRATIO*CDTRO(DTCN(JC),JB)
c              CDTR(DTCN(JC),31) = (1.0-CRATIO)*CDTRNX(DTCN(JC),JB)
c     .                            +CRATIO*CDTRO(DTCN(JC),JB)
c
            END DO
          END IF
        END IF

******* Upstream head conditions

        IF (UH_EXTERNAL(JB)) THEN
          IF (INTERP_HEAD) THEN
            HRATIO = (NXEUH1(JB)-JDAY)/(NXEUH1(JB)-NXEUH2(JB))
            TRATIO = (NXTUH1(JB)-JDAY)/(NXTUH1(JB)-NXTUH2(JB))
            IF (CONSTITUENTS) THEN
              CRATIO = (NXCUH1(JB)-JDAY)/(NXCUH1(JB)-NXCUH2(JB))
            END IF
            ELUH(JB) = (1.0-HRATIO)*ELUHNX(JB)+HRATIO*ELUHO(JB)
            DO K=2,KMP-1
              TUH(K,JB) = (1.0-TRATIO)*TUHNX(K,JB)+TRATIO*TUHO(K,JB)
            END DO
c            DO JC=1,NAC
            DO JC=1,nacin
              DO K=2,KMP-1
                CUH(K,CN(JC),JB) = (1.0-CRATIO)*CUHNX(K,CN(JC),JB)
     .                             +CRATIO*CUHO(K,CN(JC),JB)
              END DO
            END DO
          END IF
        END IF

******* Downstream head

        IF (DH_EXTERNAL(JB)) THEN
          IF (INTERP_HEAD) THEN
            HRATIO = (NXEDH1(JB)-JDAY)/(NXEDH1(JB)-NXEDH2(JB))
            TRATIO = (NXTDH1(JB)-JDAY)/(NXTDH1(JB)-NXTDH2(JB))
            IF (CONSTITUENTS) THEN
              CRATIO = (NXCDH1(JB)-JDAY)/(NXCDH1(JB)-NXCDH2(JB))
            END IF
            ELDH(JB) = (1.0-HRATIO)*ELDHNX(JB)+HRATIO*ELDHO(JB)
            DO K=2,KMP-1
              TDH(K,JB) = (1.0-TRATIO)*TDHNX(K,JB)+TRATIO*TDHO(K,JB)
            END DO
c            DO JC=1,NAC
            DO JC=1,nacin
              DO K=2,KMP-1
                CDH(K,CN(JC),JB) = (1.0-CRATIO)*CDHNX(K,CN(JC),JB)
     .                             +CRATIO*CDHO(K,CN(JC),JB)
              END DO
            END DO
          END IF
        END IF
      END DO
      END

************************************************************************
**          S U B R O U T I N E   H E A T   E X C H A N G E           **
************************************************************************

      SUBROUTINE HEAT_EXCHANGE (JDAY)

        INCLUDE 'w2.inc'

******* Variable declarations

        REAL      JDAY, LAT, LONG, LONG0
        integer   wscdp
        LOGICAL   LEAP_YEAR
        CHARACTER MONTH*9

******* Dimension declarations

        DIMENSION EQT(12)

******* Common declarations

        COMMON /GLBLCC/ PALT, ALGDET, O2LIM, WIND, wscdp, wsc(ndp)
        COMMON /GDAYC1/ GDAY, DAYM,   JDAYG, LEAP_YEAR
        COMMON /GDAYC2/ MONTH
        COMMON /TVDMTC/ TAIR, TDEW,   CLOUD, PHI, ET, CSHE, SRO,
     .                  LAT,  LONG,   banksh(imp)

******* Data declarations

        DATA EQT   /-0.13, -0.23, -0.16, -0.02, 0.06, 0.00, -0.09,
     .              -0.08,  0.06,  0.22,  0.25, 0.10/
        DATA LONG0 /90.0/

******* Convert W2 units to subroutine units

        WIND1 = WIND*2.23714
        TDEW1 = TDEW*9.0/5.0+32.0
        TAIR1 = TAIR*9.0/5.0+32.0

******* Determine month

        IF (MONTH.EQ.'  January') M = 1
        IF (MONTH.EQ.' February') M = 2
        IF (MONTH.EQ.'    March') M = 3
        IF (MONTH.EQ.'    April') M = 4
        IF (MONTH.EQ.'      May') M = 5
        IF (MONTH.EQ.'     June') M = 6
        IF (MONTH.EQ.'     July') M = 7
        IF (MONTH.EQ.'   August') M = 8
        IF (MONTH.EQ.'September') M = 9
        IF (MONTH.EQ.'  October') M = 10
        IF (MONTH.EQ.' November') M = 11
        IF (MONTH.EQ.' December') M = 12

******* Compute solar radiation

        D     = 0.409280*COS(0.017214*(172.0-JDAYG))
        X     = (JDAY-JDAYG)*24.0
        H     = 0.261799*(X-(LONG-LONG0)*0.066667+EQT(M)-12.0)
        SINAL = SIN(LAT*.017453)*SIN(D)+COS(LAT*.017453)*COS(D)*COS(H)
        AL    = ASIN(SINAL)
        A0    = 57.2985*AL
        SRO   = 2.044*A0+0.1296*A0**2-0.001941*A0**3+7.591E-6*A0**4
c        SRO   = (1.0-0.0081*CLOUD*CLOUD)*SRO*24.0
        SRO   = (1.0-0.0065*CLOUD*CLOUD)*SRO*24.0
        IF (A0.LT.0.0) SRO = 0.0

******* Compute equilibrium temperature and heat exchange coefficient

        ET    = TDEW1
        FW    = 70.0+(0.7*WIND1*WIND1)
        HA    = 3.1872E-08*(TAIR1+460.0)*(TAIR1+460.0)*(TAIR1+460.0)
     .          *(TAIR1+460.0)
        TSTAR = (ET+TDEW1)/2.0
        BETA  = 0.255-(0.0085*TSTAR)+(0.000204*TSTAR*TSTAR)
        CSHE  = 15.7+(0.26+BETA)*FW
        ETP   = (SRO+HA-1801.0)/CSHE+(CSHE-15.7)*(0.26*TAIR1+BETA*TDEW1)
     .          /(CSHE*(0.26+BETA))
        DO J=1,50
          IF (ABS(ETP-ET).LT.0.05) GO TO 10000
          ET    = ETP
          TSTAR = (ET+TDEW1)/2.0
          BETA  = 0.255-(0.0085*TSTAR)+(0.000204*TSTAR**2)
          CSHE  = 15.7+(0.26+BETA)*FW
          ETP   = (SRO+HA-1801.0)/CSHE+(CSHE-15.7)*(0.26*TAIR1+BETA
     .            *TDEW1)/(CSHE*(0.26+BETA))
        END DO
10000   CONTINUE

******* Convert to W2 units

        ET   = (ET-32.0)*5.0/9.0
        SRO  = 3.14E-8*SRO
        CSHE = 5.65E-8*CSHE
      END

************************************************************************
**             S U B R O U T I N E   R A D I A T I O N                **
************************************************************************

      SUBROUTINE RADIATION (SRO,CLOUD,TAIR,RSN,RAN)

******* Data declarations

        DATA SBC /2.0411E-7/

******* Net solar radiation

        RSN = 0.94*SRO

******* Net atmospheric radiation

        T2K = 273.2+TAIR
        FAC = 1.0+0.0017*CLOUD**2
        RAN = 1000.0/3600.0*9.37E-6*SBC*T2K**6*FAC*(1.0-0.03)

      RETURN
      END

************************************************************************
**          S U B R O U T I N E   S U R F A C E   T E R M S           **
************************************************************************

      SUBROUTINE SURFACE_TERMS (TAIR,TDEW,WIND,TS,RB,RE,RC)
        double precision ts

******* Data declarations

        DATA SBC /2.0411E-7/, ETAW /0.96/

******* Partial water vapor pressure of the air

        IF (TDEW.GT.0.0) THEN
          A = 7.5
          B = 237.3
          C = 0.6609
        ELSE
          A = 9.5
          B = 265.5
          C = 0.6609
        END IF
        EA = EXP(2.3026*(A*TDEW/(TDEW+B)+C))

******* Partial water vapor pressure at the water surface temperature

        IF (TS.GT.0.0) THEN
          A = 7.5
          B = 237.3
          C = 0.6609
        ELSE
          A = 9.5
          B = 265.5
          C = 0.6609
        END IF
        ES = EXP(2.3026*(A*TS/(TS+B)+C))
        DE = ES-EA

******* Longwave radiation

        RB = 1000.0/3600.0*ETAW*SBC*(TS+273.2)**4

******* Evaporation

        W2A = WIND*3600.0/(0.3048*5280.0)
        FWA = 70.0+0.7*W2A*W2A
        FW  = 0.1313*FWA
        RE  = FW*DE

******* Conduction

        RC = FW*0.47*(TS-TAIR)

      RETURN
      END

************************************************************************
**            S U B R O U T I N E   P R I N T   G R I D               **
************************************************************************

      SUBROUTINE PRINT_GRID (JDAY,GDAY,MONTH,YEAR)
      INCLUDE   'w2.inc'
      REAL       JDAY, ICETH
      INTEGER    SNP,  CN,  YEAR, GDAY
      LOGICAL    CONSTITUENTS, ICE, ICE_CALC, LONG_FORM, SHORT_FORM,
     .           LIMITING_FACTOR
      CHARACTER  CPRN*3,   HPRN*3,  LFAC*4,  CONV*7,    HUNIT*6,
     .           CUNIT*6,  LFPR*7,  MONTH*9, CONVD1*13, CNAME*16,
     .           HNAME*24, TITLE*72
      DOUBLE PRECISION Z1, Z2,t1,t2
      DIMENSION HUNIT(3)
      COMMON /SNPUC/  SNP
      COMMON /GLOBLC/ JB, JC, IU, ID, KT, DLT, KB(IMP)
      COMMON /TEMPC/  T1(KMP,IMP),    T2(KMP,IMP)
      COMMON /HYDRC1/ U(KMP,IMP),     W(KMP,IMP),    AZ(KMP,IMP),
     .                RHO(KMP,IMP)
      COMMON /HYDRC2/ Z1(IMP),        Z2(IMP)
      COMMON /PRNTC1/ IPR(17),        IBPR,          IEPR,      KEPR
      COMMON /PRNTC2/ TITLE(5),       CPRN(NCP),     HPRN(3),
     .                CONV(KMP,17),   CONVD1(KMP,11)
      COMMON /LFACC/  LFAC(KMC,IMC),  LFPR(KMC,IMC)
      COMMON /DKSEDC/ SEDD(KMC,IMC),  SO2D(KMC,IMC), SOD(IMP)
      COMMON /GRDLGC/ LONG_FORM,      SHORT_FORM,    LIMITING_FACTOR
      COMMON /GRTVDC/ CONSTITUENTS,   CN(NCP),       NAC
      COMMON /CBODC/  KBOD,           TBOD,          RBOD
      COMMON /NAMESC/ HNAME(3),       CNAME(NCP),    CUNIT(NCP)
      COMMON /ICEC/   ICE(IMP),       ICETH(IMP),    ICE_CALC
      DATA   HUNIT(1) /'m/sec '/, HUNIT(2) /'mm/sec'/,
     .       HUNIT(3) /'deg C '/

***** Water surface and ice cover

      WRITE (SNP,3000) 'Water Surface [Z1] (m)',(IPR(I),I=IBPR,IEPR)
      DO I=IBPR,IEPR
        WRITE (CONV(1,I),'(F7.4)') Z2(IPR(I))
      END DO
      WRITE (SNP,3010) (CONV(1,I),I=IBPR,IEPR)
      IF (ICE_CALC) THEN
        DO I=IBPR,IEPR
          WRITE (CONV(2,I),'(F7.4)') ICETH(IPR(I))
        END DO
        WRITE (SNP,3020) 'Ice Thickness (m)',(CONV(2,I),I=IBPR,IEPR)
      END IF
      DO I=IBPR,IEPR
        WRITE (CONV(1,I),'(F7.2)') SOD(IPR(I))*86400.0
      END DO
      IF (CONSTITUENTS) THEN
        WRITE (SNP,3030) 'Sediment Oxygen Demand [SOD] (g/m^2/day)',
     .                    (CONV(1,I),I=IBPR,IEPR)
      END IF

***** Velocities, temperatures, and vertical eddy viscosities

      DO J=1,3
        IF (HPRN(J).EQ.' ON') THEN
          DO I=IBPR,IEPR
            DO K=KT,KB(IPR(I))
              IF (J.EQ.1) WRITE (CONV(K,I),'(F7.4)') U(K,IPR(I))
              IF (J.EQ.2) WRITE (CONV(K,I),'(F7.4)') W(K,IPR(I))*1000.0
              IF (J.EQ.3) WRITE (CONV(K,I),'(F7.2)') T2(K,IPR(I))
            END DO
          END DO
          WRITE (SNP,3040) (TITLE(I),I=1,5)
          WRITE (SNP,3050)  MONTH,GDAY,YEAR,INT(JDAY),(JDAY-INT(JDAY))
     .                      *24.0
          IF (LONG_FORM)    WRITE (SNP,3060) HNAME(J),HUNIT(J)
          IF (SHORT_FORM)   WRITE (SNP,3070) HNAME(J),HUNIT(J)
          WRITE (SNP,3080) (IPR(I),I=IBPR,IEPR)
          DO K=KT,KEPR
            WRITE (SNP,3090) K,(CONV(K,I),I=IBPR,IEPR)
          END DO
        END IF
      END DO

***** Constituent concentrations

      DO J=1,NAC
        IF (CPRN(CN(J)).EQ.' ON') THEN
          MULT = 1.0
          IF (CN(J).GE.9.AND.CN(J).LE.11) MULT = 1000.0
          DO I=IBPR,IEPR
            DO K=KT,KB(IPR(I))
              WRITE (CONV(K,I),'(F7.2)') C2(K,IPR(I),CN(J))*MULT
            END DO
          END DO
          WRITE (SNP,3040) (TITLE(I),I=1,5)
          WRITE (SNP,3050)  MONTH,GDAY,YEAR,INT(JDAY),(JDAY-INT(JDAY))
     .                      *24.0
          IF (LONG_FORM)   WRITE (SNP,3120)  CNAME(CN(J)),CUNIT(CN(J))
          IF (SHORT_FORM)  WRITE (SNP,3130)  CNAME(CN(J)),CUNIT(CN(J))
          WRITE (SNP,3080) (IPR(I),I=IBPR,IEPR)
          DO K=KT,KEPR
            WRITE (SNP,3090) K,(CONV(K,I),I=IBPR,IEPR)
          END DO
        END IF
      END DO
      IF (CONSTITUENTS.AND.LIMITING_FACTOR) THEN
        WRITE (SNP,3040) (TITLE(I),I=1,5)
        WRITE (SNP,3050)  MONTH,GDAY,YEAR,INT(JDAY),(JDAY-INT(JDAY))*24.
        IF (LONG_FORM)    WRITE (SNP,3100)
        IF (SHORT_FORM)   WRITE (SNP,3110)
        WRITE (SNP,3080) (IPR(I),I=IBPR,IEPR)
        DO K=KT,KEPR
          WRITE (SNP,3090) K,(LFPR(K,IPR(I)),I=IBPR,IEPR)
        END DO
      END IF
      RETURN

***** Snapshot FORMATs

 3000 FORMAT(/3X,A22/
     .       '+',2X,5('_'),1X,7('_')//
     .       3X,17I7)
 3010 FORMAT(3X,17A7/)
 3020 FORMAT(/3X,A17/
     .       '+',2X,3('_'),1X,9('_')//
     .       3X,17A7)
 3030 FORMAT(3X,A40/
     .       '+',2X,8('_'),1X,6('_'),1X,6('_')//
     .       3X,17A7/)
 3035 FORMAT(3X,A44/
     .       '+',2X,7('_'),1X,6('_'),1X,6('_'),1X,6('_')//
     .       3X,17A7/)
 3040 FORMAT('1',5(A72/1X))
 3050 FORMAT(1X,A9,I3,', ',I4,'       Julian Date =',I6,' days',F6.2,
     .       ' hours'/)
 3060 FORMAT(47X,A24,1X,A6/)
 3070 FORMAT(22X,A24,1X,A6/)
 3080 FORMAT(2X,17I7)
 3090 FORMAT(1X,I2,17A7)
 3100 FORMAT(45X,'Limiting Factor'/)
 3110 FORMAT(20X,'Limiting Factor'/)
 3120 FORMAT(45X,A16,1X,A6/)
 3130 FORMAT(20X,A16,1X,A6/)
      END

************************************************************************
**   S U B R O U T I N E   S E L E C T I V E   W I T H D R A W A L    **
************************************************************************

      SUBROUTINE SELECTIVE_WITHDRAWAL
      INCLUDE 'w2.inc'

***** Type declarations

      LOGICAL SEL_WITHDRAWAL, POINT_SINK

***** Dimension declarations

      DIMENSION VNORM(KMP)

***** Common declarations

      COMMON /GLOBLC/ JB, JC, IU, ID, KT, DLT, KB(IMP)
      COMMON /GEOMHC/ EL(KMP),       H(KMP),        HKT1(IMP),
     .                HKT2(IMP)
      COMMON /TVDQC/  QIN(NBP),      QTR(NTP),      QDTR(NBP),
     .                PR(NBP),       ELUH(NBP),     ELDH(NBP),
     .                QOUT(KMP,NBP), QWD(NWP)
      COMMON /SELWC/  NSTR(NBP),     QSTR(NSP,NBP), ESTR(NSP,NBP),
     .                WSTR(NSP,NBP), KBSW(NBP),     KTOPSW(NBP),
     .                KBOTSW(NBP),   NOUT(NBP),     KOUT(KMP,NBP)
      COMMON /HYDRC1/ U(KMP,IMP),    W(KMP,IMP),    AZ(KMP,IMP),
     .                RHO(KMP,IMP)
      COMMON /TVDSWC/ SEL_WITHDRAWAL(NBP), POINT_SINK(NSP,NBP)

***** Data declarations

      DATA G /9.81/

***** Initialize variables

      DO K=1,KMP
        KOUT(K,JB) = 0
        QOUT(K,JB) = 0.0
        VNORM(K)   = 0.0
      END DO

***** Calculations

      DO JS=1,NSTR(JB)
        DO K=KT,KB(ID)
          IF (EL(K).GE.ESTR(JS,JB)) KSTR = K
        END DO
        KTOP = KT
        KBOT = KBSW(JB)

******* Boundary interference

        RATIO = (ESTR(JS,JB)-EL(KBOT))/(EL(KT)-EL(KBOT))
        COEF  = 1.0
        IF (RATIO.LT.0.10.OR.RATIO.GT.0.90) COEF = 2.0

******* Average density frequency of region above the structure

        DO K=KSTR-1,KT,-1
          DTOP = EL(K)-ESTR(JS,JB)
          DFT  = SQRT((ABS(RHO(K,ID)-RHO(KSTR,ID)))/DTOP
     .           /RHO(KSTR,ID)*G)
          DFT  = MAX(DFT,1.0E-10)

********* Half-height of withdrawal zone

          IF (POINT_SINK(JS,JB)) THEN
            ZTOP = (COEF*QSTR(JS,JB)/DFT)**.333333
          ELSE
            ZTOP = SQRT(2.0*COEF*QSTR(JS,JB)/WSTR(JS,JB)/DFT)
          END IF
          IF (DTOP.GE.ZTOP) THEN
            KTOP = K
            GO TO 10000
          END IF
        END DO
10000   CONTINUE

******* Upper withdrawal layer limit

        ELTOP = ESTR(JS,JB)+ZTOP
        IF (ELTOP.LT.EL(KT)) THEN
          RDT = ABS(RHO(KSTR,ID)-RHO(KTOP,ID))
        ELSE
          RDT = ABS(RHO(KSTR,ID)-RHO(KT,ID))*ZTOP/(EL(KT)
     .          -ESTR(JS,JB))
        END IF
        RDT = MAX(RDT,1.0E-10)

******* Average density frequency of region below the structure

        DO K=KSTR+1,KBSW(JB)
          DBOT = ESTR(JS,JB)-EL(K)
          DFB  = SQRT((ABS(RHO(K,ID)-RHO(KSTR,ID)))/DBOT
     .           /RHO(KSTR,ID)*G)
          DFB  = MAX(DFB,1.0E-10)

********* Withdrawal zone half-height

          IF (POINT_SINK(JS,JB)) THEN
            ZBOT = (COEF*QSTR(JS,JB)/DFB)**0.333333
          ELSE
            ZBOT = SQRT(2.0*COEF*QSTR(JS,JB)/WSTR(JS,JB)/DFB)
          END IF
          IF (DBOT.GE.ZBOT) THEN
            KBOT = K
            GO TO 10010
          END IF
        END DO
10010   CONTINUE

******* Reference density

        ELBOT = ESTR(JS,JB)-ZBOT
        IF (ELBOT.GT.EL(KBSW(JB)+1)) THEN
          RDB = ABS(RHO(KSTR,ID)-RHO(KBOT,ID))
        ELSE
          RDB = ABS(RHO(KSTR,ID)-RHO(KBSW(JB),ID))*ZBOT
     .          /(ESTR(JS,JB)-EL(KBSW(JB)+1))
        END IF
        RDB = MAX(RDB,1.0E-10)

******* Velocity profile

        VT     = 0.0
        DENDIF = MAX(RDT,RDB)
        DO K=KTOP,KBOT
          JO        = K-KTOP+1
          VNORM(JO) = 1.0-((RHO(K,ID)-RHO(KSTR,ID))/DENDIF)**2
          VNORM(JO) = MAX(VNORM(JO),0.0)
          VT        = VT+VNORM(JO)
        END DO

******* Layer flow rate

        IF (KBOT.EQ.KMP) KBOT = KMP-1
        DO JO=1,KBOT-KTOP+1
          QOUT(JO,JB) = VNORM(JO)/VT*QSTR(JS,JB)+QOUT(JO,JB)
        END DO

******* Outflow limits for all structures

        KTOPSW(JB) = MIN(KTOP,KTOPSW(JB))
        KBOTSW(JB) = MAX(KBOT,KBOTSW(JB))
      END DO

***** Number of outflows and layer locations for all structures combined

      NOUT(JB) = KBOTSW(JB)-KTOPSW(JB)+1
      DO K=KTOPSW(JB),KBOTSW(JB)
        KOUT(K-KTOPSW(JB)+1,JB) = K
      END DO

***** Reset boundary horizontal velocities for layers with no flow

      DO K=KT,KTOPSW(JB)-1
        U(K,ID) = 0.0
      END DO
      DO K=KBOTSW(JB)+1,KMP-1
        U(K,ID) = 0.0
      END DO
      DO JO=1,NOUT(JB)
        IF (QOUT(JO,JB).EQ.0.0) U(KOUT(JO,JB),ID) = 0.0
      END DO
      RETURN
      END

************************************************************************
**       S U B R O U T I N E   R A T E   M U L T I P L I E R S        **
************************************************************************

      SUBROUTINE RATE_MULTIPLIERS
      INCLUDE   'w2.inc'
      REAL    LAM1,  LAM2,  NH3K1, NH3K2, NO3K1, NO3K2, NH3T1, NH3T2,
     .        NO3T1, NO3T2, NH3RM, NO3RM, K1,     K2,   K3,    K4
      double precision t1,t2
      COMMON /GLOBLC/ JB, JC, IU, ID, KT, DLT, KB(IMP)
      COMMON /GBLRTC/ OMRM(KMC,IMC),  NH3RM(KMC,IMC), NO3RM(KMC,IMC),
     .                AGRMR(KMC,IMC), AGRMF(KMC,IMC),
     .                zoormr(kmc,imc), zoormf(kmc,imc)
      COMMON /RTMLTC/ OMT1,  OMT2,  NH3T1, NH3T2, NO3T1, NO3T2, AGT1,
     .                AGT2,  AGT3,  AGT4,  OMK1,  OMK2,  NH3K1, NH3K2,
     .                NO3K1, NO3K2, AGK1,  AGK2,  AGK3,  AGK4,
     .                zoot1,  zoot2,  zoot3,  zoot4, zook1,  zook2,
     .                zook3,  zook4
      COMMON /TEMPC/  T1(KMP,IMP),  T2(KMP,IMP)

***** Rising and falling temperature rate functions

      FR(TT,TT1,TT2,K1,K2) = K1*EXP(LOG(K2*(1.0-K1)/(K1*(1.0-K2)))
     .                       /(TT2-TT1)*(TT-TT1))
      FF(TT,TT3,TT4,K3,K4) = K4*EXP(LOG(K3*(1.0-K4)/(K4*(1.0-K3)))
     .                       /(TT4-TT3)*(TT4-TT))
      DO I=IU,ID
        DO K=KT,KB(I)
          ttt=T1(K,I)
          LAM1       = FR(TTT,NH3T1,NH3T2,NH3K1,NH3K2)
          NH3RM(K,I) = LAM1/(1.0+LAM1-NH3K1)
          LAM1       = FR(TTT,NO3T1,NO3T2,NO3K1,NO3K2)
          NO3RM(K,I) = LAM1/(1.0+LAM1-NO3K1)
          LAM1       = FR(TTT,OMT1,OMT2,OMK1,OMK2)
          OMRM(K,I)  = LAM1/(1.0+LAM1-OMK1)
          LAM1       = FR(TTT,AGT1,AGT2,AGK1,AGK2)
          LAM2       = FF(TTT,AGT3,AGT4,AGK3,AGK4)
          AGRMR(K,I) = LAM1/(1.0+LAM1-AGK1)
          AGRMF(K,I) = LAM2/(1.0+LAM2-AGK4)
c for zooplankton
          LAM1       = FR(TTT,zoot1,zoot2,zook1,zook2)
          LAM2       = FF(TTT,zoot3,zoot4,zook3,zook4)
          zoormr(k,i)= lam1/(1.+lam1-zook1)
          zoormf(k,i)= lam2/(1.+lam2-zook4)
        END DO
      END DO
      END

************************************************************************
**       S U B R O U T I N E   D E C A Y   C O N S T A N T S          **
************************************************************************

      SUBROUTINE DECAY_CONSTANTS
      INCLUDE   'w2.inc'
      REAL    LABDK,  NH3DK, NO3DK, NH3D, NO3D, NH3RM, NO3RM, KBOD
      integer wscdp
      double precision t1,t2
      COMMON /GLOBLC/ JB, JC, IU, ID, KT, DLT, KB(IMP)
      COMMON /DKORGC/ DETD(KMC,IMC),  ORGD(KMC,IMC)
      COMMON /DKSEDC/ SEDD(KMC,IMC),  SO2D(KMC,IMC),  SOD(IMP)
      COMMON /DKNITC/ NH3D(KMC,IMC),  NO3D(KMC,IMC)
      COMMON /DKBODC/ CBODD(KMC,IMC)
      COMMON /DKMLTC/ A1(KMC,IMC),    A2(KMC,IMC), A3(KMC,IMC)
      COMMON /GLBLCC/ PALT, ALGDET,   O2LIM, WIND, wscdp, wsc(ndp)
      COMMON /GBLRTC/ OMRM(KMC,IMC),  NH3RM(KMC,IMC), NO3RM(KMC,IMC),
     .                AGRMR(KMC,IMC), AGRMF(KMC,IMC),
     .                zoormr(kmc,imc), zoormf(kmc,imc)
      COMMON /GEOMHC/ EL(KMP),        H(KMP),         HKT1(IMP),
     .                HKT2(IMP)
      COMMON /GEOMBC/ B(KMP,IMP),     BKT(IMP),       BH(KMP,IMP),
     .                BHKT1(IMP),     BHKT2(IMP),     BHRKT1(IMP)
      COMMON /SETLC1/ SETIN(KMC,IMC), SETOUT(KMC,IMC)
      COMMON /TEMPC/  T1(KMP,IMP),    T2(KMP,IMP)
      COMMON /SETLC2/ SSETL, DSETL, ASETL, FESETL
      COMMON /ORGDKC/ SEDDK, DETDK, LABDK, REFDK,  LRFDK
      COMMON /NITROC/ BION,  PARTN, NH3DK, NH3REL, NO3DK, obion
      COMMON /CBODC/  KBOD,  TBOD,  RBOD

      DO I=IU,ID
        DO K=KT,KB(I)
          A1(K,I)     = (1.0+SIGN(1.0,DO(K,I)-O2LIM))*0.5
          A2(K,I)     = (1.0+SIGN(1.0,O2LIM-DO(K,I)))*0.5
          A3(K,I)     = (1.0+SIGN(1.0,DO(K,I)))*0.5
          ORGD(K,I)   = OMRM(K,I)*(LABDK*LABDOM(K,I)+REFDK*REFDOM(K,I))
     .                  *A3(K,I)
          SO2D(K,I)   = SOD(I)/BH(K,I)*OMRM(K,I)*(B(K,I)-B(K+1,I))
          NH3D(K,I)   = NH3DK*NH3RM(K,I)*NH3(K,I)*A1(K,I)
          NO3D(K,I)   = NO3DK*NO3RM(K,I)*NO3(K,I)*A2(K,I)
          CBODD(K,I)  = KBOD*TBOD**(T1(K,I)-20.0)*A3(K,I)
          DETD(K,I)   = DETDK*OMRM(K,I)*DETRIT(K,I)*A3(K,I)
          SEDD(K,I)   = SEDDK*OMRM(K,I)*SEDMNT(K,I)*A3(K,I)
          SETIN(K,I)  = (SSETL*SS(K-1,I)+FESETL*FE(K-1,I))/H(K)
     .                  *A1(K,I)
          SETOUT(K,I) = (SSETL*SS(K,I)+FESETL*FE(K,I))/H(K)*A1(K,I)
        END DO
        SETOUT(KT,I)  = SETOUT(KT,I)*H(K)/HKT2(I)*A1(KT,I)
        SO2D(KT,I)    = SOD(I)/BHKT2(I)*OMRM(KT,I)*(BKT(I)-B(KT+1,I))
        SO2D(KB(I),I) = SOD(I)/BH(KB(I),I)*OMRM(KB(I),I)*B(KB(I),I)
      END DO
      END

************************************************************************
**      S U B R O U T I N E   S U S P E N D E D   S O L I D S         **
************************************************************************

      SUBROUTINE SUSPENDED_SOLIDS
      INCLUDE   'w2.inc'
      COMMON /GLOBLC/ JB, JC,  IU, ID, KT, DLT, KB(IMP)
      COMMON /GEOMHC/ EL(KMP), H(KMP), HKT1(IMP), HKT2(IMP)
      COMMON /SETLC2/ SSETL,   DSETL,  ASETL, FESETL

      DO I=IU,ID
        SSSS(KT,I) = -SSETL*SS(KT,I)/HKT2(I)
        DO K=KT+1,KB(I)
          SSSS(K,I) = SSETL*(SS(K-1,I)-SS(K,I))/H(K)
        END DO
      END DO
      END

************************************************************************
**               S U B R O U T I N E   C O L I F O R M                **
************************************************************************

      SUBROUTINE COLIFORM
      INCLUDE   'w2.inc'
      double precision t1,t2
      COMMON /CLFRMC/ COLQ10, COLDK
      COMMON /GLOBLC/ JB, JC, IU, ID, KT, DLT, KB(IMP)
      COMMON /TEMPC/  T1(KMP,IMP), T2(KMP,IMP)

      DO I=IU,ID
        DO K=KT,KB(I)
          COLSS(K,I) = -COLDK*COLQ10**(T1(K,I)-20.0)*colfrm(k,i)
        END DO
      END DO
      END

************************************************************************
**            S U B R O U T I N E   L A B I L E   D O M               **
************************************************************************

      SUBROUTINE LABILE_DOM
      INCLUDE   'w2.inc'
      REAL    LABDK, LRFDK
      integer  wscdp
      COMMON /GLOBLC/ JB, JC, IU, ID, KT, DLT, KB(IMP)
      COMMON /GLBLCC/ PALT, ALGDET, O2LIM, WIND, wscdp, wsc(ndp)
      COMMON /GBLRTC/ OMRM(KMC,IMC),  NH3RM(KMC,IMC), NO3RM(KMC,IMC),
     .                AGRMR(KMC,IMC), AGRMF(KMC,IMC),
     .                zoormr(kmc,imc), zoormf(kmc,imc)
      COMMON /DKMLTC/ A1(KMC,IMC),    A2(KMC,IMC), A3(KMC,IMC)
      COMMON /PHYTGC/ AGR(KMC,IMC),   ARR(KMC,IMC),   AMR(KMC,IMC),
     .                AER(KMC,IMC)
      COMMON /ORGDKC/ SEDDK, DETDK, LABDK,  REFDK, LRFDK

      DO I=IU,ID
        DO K=KT,KB(I)
          DECAY       = OMRM(K,I)*A3(K,I)*(LABDK+LRFDK)*LABDOM(K,I)
          ALGP        = (AER(K,I)+(1.0-ALGDET)*AMR(K,I))*ALGAE(K,I)
          LDOMSS(K,I) = ALGP-DECAY
        END DO
      END DO
      END

************************************************************************
**         S U B R O U T I N E   R E F R A C T O R Y   D O M          **
************************************************************************

      SUBROUTINE REFRACTORY_DOM
      INCLUDE   'w2.inc'
      REAL    LRFDK
      COMMON /GLOBLC/ JB, JC, IU, ID, KT, DLT, KB(IMP)
      COMMON /GBLRTC/ OMRM(KMC,IMC),  NH3RM(KMC,IMC), NO3RM(KMC,IMC),
     .                AGRMR(KMC,IMC), AGRMF(KMC,IMC),
     .                zoormr(kmc,imc), zoormf(kmc,imc)
      COMMON /DKMLTC/ A1(KMC,IMC),    A2(KMC,IMC), A3(KMC,IMC)
      COMMON /ORGDKC/ SEDDK, DETDK, LABDK, REFDK, LRFDK

      DO I=IU,ID
        DO K=KT,KB(I)
          RDOMSS(K,I) = OMRM(K,I)*(LRFDK*LABDOM(K,I)-REFDK*REFDOM(K,I)
     .                  *A3(K,I))
        END DO
      END DO
      END

************************************************************************
**          S U B R O U T I N E   P H Y T O P L A N K T O N           **
************************************************************************

      SUBROUTINE PHYTOPLANKTON
      INCLUDE   'w2.inc'
      REAL      LAM1, LAM2, LTCOEF, LLIM, NLIM, LIMIT, NETSET
      CHARACTER LF*3, LFAC*4, LFPR*7
      COMMON /PHYTGC/ AGR(KMC,IMC),   ARR(KMC,IMC), AMR(KMC,IMC),
     .                AER(KMC,IMC)
      COMMON /LFACC/  LFAC(KMC,IMC),  LFPR(KMC,IMC)
      COMMON /DKMLTC/ A1(KMC,IMC),    A2(KMC,IMC), A3(KMC,IMC)
      COMMON /GBLRTC/ OMRM(KMC,IMC),  NH3RM(KMC,IMC), NO3RM(KMC,IMC),
     .                AGRMR(KMC,IMC), AGRMF(KMC,IMC),
     .                zoormr(kmc,imc), zoormf(kmc,imc)
      COMMON /GEOMHC/ EL(KMP),        H(KMP),         HKT1(IMP),
     .                HKT2(IMP)
      COMMON /PHYTC1/ AEXCR,  AMORT,  AGROW,  ARESP, ASATUR, AHSN, AHSP
      COMMON /PHYTC2/ BETA,   EXH2O,  EXINOR, EXORG
      COMMON /SETLC2/ SSETL,  DSETL,  ASETL,  FESETL
      COMMON /GLOBLC/ JB,     JC,     IU,    ID,    KT,     DLT,
     .                KB(IMP)
      COMMON /TVDMTC/ TAIR,   TDEW,   CLOUD,  PHI,   ET,     CSHE,
     .                SRO,    LAT,    LONG,   banksh(imp)
      COMMON /PHOSPC/ PO4REL, BIOP,   PARTP, obiop
      COMMON /NITROC/ BION,   PARTN,  NH3DK,  NH3REL, NO3DK, obion
      common /zoopc/ zmax,zmort,zeffic,pref1,pref2,zresp,zoomin,zs2p,
     .               zmu(kmc,imc),talgae(kmc,imc),zrt(kmc,imc),
     .               zmt(kmc,imc)
c
      LTCOEF = (1.0-BETA)*SRO*4.186E6/ASATUR
      DO I=IU,ID

******* Limiting factor

        GAMMA = EXH2O+EXINOR*SS(KT,I)+EXORG*(ALGAE(KT,I)+DETRIT(KT,I)
     .          +zoop(kt,i))
        DEPTH = HKT2(I)
        LAM1  = LTCOEF*banksh(i)
        LAM2  = LTCOEF*banksh(i)*EXP(-GAMMA*DEPTH)
        LLIM  = 2.718282/(GAMMA*DEPTH)*(EXP(-LAM2)-EXP(-LAM1))
        PLIM  = PO4(KT,I)/(PO4(KT,I)+AHSP)
        NLIM  = (NH3(KT,I)+NO3(KT,I))/(NH3(KT,I)+NO3(KT,I)+AHSN)
        LIMIT = MIN(PLIM,NLIM,LLIM)
        IF (LIMIT.EQ.PLIM) THEN
          WRITE (LFAC(KT,I),'(F4.3)') PLIM
          LF         = ' P '
          LFPR(KT,I) = LF//LFAC(KT,I)
        ELSE IF (LIMIT.EQ.NLIM) THEN
          WRITE (LFAC(KT,I),'(F4.3)') NLIM
          LF         = ' N '
          LFPR(KT,I) = LF//LFAC(KT,I)
        ELSE IF (LIMIT.EQ.LLIM) THEN
          WRITE (LFAC(KT,I),'(F4.3)') LLIM
          LF         = ' L '
          LFPR(KT,I) = LF//LFAC(KT,I)
        END IF

******* Sources/sinks

        ARR(KT,I)   = AGRMR(KT,I)*AGRMF(KT,I)*ARESP*A3(KT,I)
        AMR(KT,I)   = (AGRMR(KT,I)+1.0-AGRMF(KT,I))*AMORT
        AGR(KT,I)   = AGRMR(KT,I)*AGRMF(KT,I)*AGROW*LIMIT
        AGR(KT,I)   = MIN(AGR(KT,I),PO4(KT,I)/(BIOP*DLT*ALGAE(KT,I)
     .                +1.0E-12),(NH3(KT,I)+NO3(KT,I))/(BION*DLT
     .                *ALGAE(KT,I)+1.0E-12))
        AER(KT,I)   = MIN((1.0-LLIM)*AEXCR,AGR(KT,I))
        GROWTH      = (AGR(KT,I)-ARR(KT,I)-AER(KT,I)-AMR(KT,I))
     .                *ALGAE(KT,I)
        NETSET      = -ASETL*ALGAE(KT,I)/HKT2(I)
c zooplankton grazing
        talgae(kt,i)  = algae(kt,i)*pref1+
     .                  pref2*detrit(kt,i)
        zmu(kt,i)     = (zoormr(kt,i)*zoormf(kt,i))*zmax*
     .                  ((talgae(kt,i)-zoomin)/(talgae(kt,i)+zs2p))
        if(zmu(kt,i).lt.0..or.zoomin.gt.talgae(kt,i))then
        zmu(kt,i)=0.0
        zoopp=0.0
        go to 55
        end if
        if(do(kt,i).lt.2.0)zmu(kt,i)=0.0
        zoopp         = zmu(kt,i)*zoop(kt,i)*pref1*(algae(kt,i)/
     .                  talgae(kt,i))
55      ALGSS(KT,I) = GROWTH+NETSET-zoopp
        DO K=KT+1,KB(I)

********* Limiting factor

          GAMMA = EXH2O+EXINOR*SS(KT,I)+EXORG*(ALGAE(K,I)+DETRIT(K,I)
     .            +zoop(k,i))
          LAM1  = LTCOEF*EXP(-GAMMA*DEPTH)
          LAM2  = LTCOEF*EXP(-GAMMA*(DEPTH+H(K)))
          LLIM  = 2.718282*(EXP(-LAM2)-EXP(-LAM1))/(GAMMA*H(K))
          DEPTH = DEPTH+H(K)
          PLIM  = PO4(K,I)/(PO4(K,I)+AHSP)
          NLIM  = (NH3(K,I)+NO3(K,I))/(NH3(K,I)+NO3(K,I)+AHSN)
          LIMIT = MIN(PLIM,NLIM,LLIM)
          IF (LIMIT.EQ.PLIM) THEN
            WRITE (LFAC(K,I),'(F4.3)') PLIM
            LF        = ' P '
            LFPR(K,I) = LF//LFAC(K,I)
          ELSE IF (LIMIT.EQ.NLIM) THEN
            WRITE (LFAC(K,I),'(F4.3)') NLIM
            LF        = ' N '
            LFPR(K,I) = LF//LFAC(K,I)
          ELSE IF (LIMIT.EQ.LLIM) THEN
            WRITE (LFAC(K,I),'(F4.3)') LLIM
            LF        = ' L '
            LFPR(K,I) = LF//LFAC(K,I)
          END IF

********* Sources/sinks

          ARR(K,I)   = AGRMR(K,I)*AGRMF(K,I)*ARESP*A3(K,I)
          AMR(K,I)   = (AGRMR(K,I)+1.0-AGRMF(K,I))*AMORT
          AGR(K,I)   = AGRMR(K,I)*AGRMF(K,I)*AGROW*LIMIT
          AGR(K,I)   = MIN(AGR(K,I),PO4(K,I)/(BIOP*DLT*ALGAE(K,I)
     .                 +1.0E-16),(NH3(K,I)+NO3(K,I))/(BION*DLT
     .                 *ALGAE(K,I)+1.0E-16))
          AER(K,I)   = MIN((1.0-LLIM)*AEXCR,AGR(K,I))
          GROWTH     = (AGR(K,I)-ARR(K,I)-AER(K,I)-AMR(K,I))*ALGAE(K,I)
          NETSET     = ASETL*(ALGAE(K-1,I)-ALGAE(K,I))/H(K)
c zooplankton grazing
        talgae(k,i)   = algae(k,i)*pref1+
     .                  pref2*detrit(k,i)
        zmu(k,i)      = (zoormr(k,i)*zoormf(k,i))*zmax*
     .                 ((talgae(k,i)-zoomin)/(talgae(k,i)+zs2p))
        if(zmu(k,i).lt.0..or.zoomin.gt.talgae(k,i))then
        zmu(k,i)=0.0
        zoopp=0.0
        go to 56
        end if
        if(do(k,i).lt.2.0)zmu(k,i)=0.0
        zoopp         = zmu(k,i)*zoop(k,i)*pref1*(algae(k,i)/
     .                  talgae(k,i))
56          ALGSS(K,I) = GROWTH+NETSET-zoopp
        END DO
      END DO
      END

************************************************************************
**               S U B R O U T I N E   D E T R I T U S                **
************************************************************************

      SUBROUTINE DETRITUS
      INCLUDE   'w2.inc'
      REAL    NETSET
      integer wscdp
      COMMON /GLOBLC/ JB, JC, IU, ID, KT, DLT, KB(IMP)
      COMMON /GLBLCC/ PALT,  ALGDET, O2LIM, WIND, wscdp, wsc(ndp)
      COMMON /SETLC2/ SSETL, DSETL,  ASETL, FESETL
      COMMON /GEOMHC/ EL(KMP),       H(KMP),       HKT1(IMP),
     .                HKT2(IMP)
      COMMON /DKORGC/ DETD(KMC,IMC), ORGD(KMC,IMC)
      COMMON /PHYTGC/ AGR(KMC,IMC),  ARR(KMC,IMC), AMR(KMC,IMC),
     .                AER(KMC,IMC)
      common /zoopc/ zmax,zmort,zeffic,pref1,pref2,zresp,zoomin,zs2p,
     .               zmu(kmc,imc),talgae(kmc,imc),zrt(kmc,imc),
     .               zmt(kmc,imc)
      COMMON /GBLRTC/ OMRM(KMC,IMC),  NH3RM(KMC,IMC), NO3RM(KMC,IMC),
     .                AGRMR(KMC,IMC), AGRMF(KMC,IMC),
     .                zoormr(kmc,imc), zoormf(kmc,imc)

      DO I=IU,ID
        ALGP        = ALGDET*AMR(KT,I)*ALGAE(KT,I)
        NETSET      = -DSETL*DETRIT(KT,I)/HKT2(I)
        if(talgae(kt,i).ne.0.0)then
        detzooin      = zoop(kt,i)*pref2*zmu(kt,i)*detrit(kt,i)/
     .                  talgae(kt,i)
        else
        detzooin=0.0
        end if
c determination of zmt and zrt for zooplankton
        tmpfac=1.-zoormf(kt,i)
        if(tmpfac.lt.0.02)tmpfac=0.02
        fmf=1.0
        if(do(kt,i).lt.2.)fmf=2.0
        zmt(kt,i)=fmf*tmpfac*zmort
        zrt(kt,i)=zoormr(kt,i)*zresp
        if(zoop(kt,i).lt.0.001)then
                       zmt(kt,i)=0.0
                       zrt(kt,i)=0.0
        end if
c
        detzooout    = zoop(kt,i)*(zmt(kt,i)+(zmu(kt,i)-
     .                 (zmu(kt,i)*zeffic)))
c
        DETSS(KT,I) = ALGP-DETD(KT,I)+NETSET+detzooout-detzooin
        DO K=KT+1,KB(I)
          ALGP       = ALGDET*AMR(K,I)*ALGAE(K,I)
          NETSET     = DSETL*(DETRIT(K-1,I)-DETRIT(K,I))/H(K)
          if(talgae(k,i).ne.0.0)then
          detzooin      = zoop(k,i)*pref2*zmu(k,i)*detrit(k,i)/
     .                  talgae(k,i)
          else
          detzooin = 0.0
          end if
c determination of zmt and zrt for zooplankton
        tmpfac=1.-zoormf(k,i)
        if(tmpfac.lt.0.02)tmpfac=0.02
        fmf=1.0
        if(do(k,i).lt.2.)fmf=2.0
        zmt(k,i)=fmf*tmpfac*zmort
        zrt(k,i)=zoormr(k,i)*zresp
        if(zoop(k,i).lt.0.001)then
                       zmt(k,i)=0.0
                       zrt(k,i)=0.0
        end if
c
        detzooout    = zoop(k,i)*(zmt(k,i)+
     .                 (zmu(k,i)-(zmu(k,i)*zeffic)))
c
        DETSS(K,I) = ALGP-DETD(K,I)+NETSET+detzooout-detzooin
        END DO
      END DO
      END

************************************************************************
**            S U B R O U T I N E   P H O S P H O R O U S             **
************************************************************************

      SUBROUTINE PHOSPHOROUS
      INCLUDE   'w2.inc'
      COMMON /GLOBLC/ JB, JC, IU, ID, KT, DLT, KB(IMP)
      COMMON /SETLC1/ SETIN(KMC,IMC), SETOUT(KMC,IMC)
      COMMON /DKORGC/ DETD(KMC,IMC),  ORGD(KMC,IMC)
      COMMON /DKSEDC/ SEDD(KMC,IMC),  SO2D(KMC,IMC), SOD(IMP)
      COMMON /DKMLTC/ A1(KMC,IMC),    A2(KMC,IMC), A3(KMC,IMC)
      COMMON /PHYTGC/ AGR(KMC,IMC),   ARR(KMC,IMC), AMR(KMC,IMC),
     .                AER(KMC,IMC)
      COMMON /PHOSPC/ PO4REL, BIOP, PARTP, obiop
      common /zoopc/ zmax,zmort,zeffic,pref1,pref2,zresp,zoomin,zs2p,
     .               zmu(kmc,imc),talgae(kmc,imc),zrt(kmc,imc),
     .               zmt(kmc,imc)

      DO I=IU,ID
        DO K=KT,KB(I)
          ALGP       = (ARR(K,I)-AGR(K,I))*ALGAE(K,I)
c zooplankton
          zoopp=zrt(k,i)*zoop(k,i)
c          PO4SS(K,I) = BIOP*(ALGP+DETD(K,I)+ORGD(K,I)+SEDD(K,I)
c       PO4SS(K,I)=BIOP*(zoopp+ALGP+DETD(K,I)+ORGD(K,I)+SEDD(K,I))
c     .                 +PO4REL*SO2D(K,I)*A2(K,I)+PARTP
c     .                 *(SETIN(K,I)*PO4(K-1,I)-SETOUT(K,I)*PO4(K,I))
      PO4SS(K,I)=BIOP*(zoopp+ALGP+DETD(K,I))+
     .                  obiop*(ORGD(K,I)+SEDD(K,I))
     .                 +PO4REL*SO2D(K,I)*A2(K,I)+PARTP
     .                 *(SETIN(K,I)*PO4(K-1,I)-SETOUT(K,I)*PO4(K,I))
        END DO
      END DO
      END

************************************************************************
**               S U B R O U T I N E   A M M O N I A                  **
************************************************************************

      SUBROUTINE AMMONIA
      INCLUDE   'w2.inc'
      REAL    NH3REL, NH3D, NO3D,  NH3DK, NO3DK
      COMMON /NITROC/ BION, PARTN, NH3DK, NH3REL, NO3DK, obion
      COMMON /GLOBLC/ JB, JC, IU, ID, KT, DLT, KB(IMP)
      COMMON /DKMLTC/ A1(KMC,IMC),    A2(KMC,IMC), A3(KMC,IMC)
      COMMON /SETLC1/ SETIN(KMC,IMC), SETOUT(KMC,IMC)
      COMMON /DKORGC/ DETD(KMC,IMC),  ORGD(KMC,IMC)
      COMMON /DKSEDC/ SEDD(KMC,IMC),  SO2D(KMC,IMC), SOD(IMP)
      COMMON /DKNITC/ NH3D(KMC,IMC),  NO3D(KMC,IMC)
      COMMON /PHYTGC/ AGR(KMC,IMC),   ARR(KMC,IMC),  AMR(KMC,IMC),
     .                AER(KMC,IMC)
      common /zoopc/ zmax,zmort,zeffic,pref1,pref2,zresp,zoomin,zs2p,
     .               zmu(kmc,imc),talgae(kmc,imc),zrt(kmc,imc),
     .               zmt(kmc,imc)

      DO I=IU,ID
        DO K=KT,KB(I)
          ALGP       = (ARR(K,I)-AGR(K,I)*NH3(K,I)/(NH3(K,I)
     .                 +NO3(K,I)+1.0E-20))*ALGAE(K,I)
c zooplankton
          zoopp=zrt(k,i)*zoop(k,i)
c          NH3SS(K,I) = BION*(ALGP+DETD(K,I)+ORGD(K,I)+SEDD(K,I)
       NH3SS(K,I)=BION*(zoopp+ALGP+DETD(K,I))
     .                 +obion*(ORGD(K,I)+SEDD(K,I))
     .                 +NH3REL*SO2D(K,I)*A2(K,I)+NO3D(K,I)
     .                 -NH3D(K,I)+PARTN*(SETIN(K,I)*NH3(K-1,I)
     .                 -SETOUT(K,I)*NH3(K,I))
        END DO
      END DO
      END

************************************************************************
**                S U B R O U T I N E   N I T R A T E                 **
************************************************************************

      SUBROUTINE NITRATE
      INCLUDE   'w2.inc'
      REAL    NH3D, NO3D, NH3DK, NO3DK
      COMMON /NITROC/ BION, PARTN, NH3DK, NH3REL, NO3DK, obion
      COMMON /GLOBLC/ JB, JC, IU, ID, KT, DLT, KB(IMP)
      COMMON /DKNITC/ NH3D(KMC,IMC), NO3D(KMC,IMC)
      COMMON /PHYTGC/ AGR(KMC,IMC),  ARR(KMC,IMC),  AMR(KMC,IMC),
     .                AER(KMC,IMC)

      DO I=IU,ID
        DO K=KT,KB(I)
          ALGC       = BION*(1.0-NH3(K,I)/(NH3(K,I)+NO3(K,I)+1.0E-20))
     .                 *AGR(K,I)*ALGAE(K,I)
          NO3SS(K,I) = NH3D(K,I)-NO3D(K,I)-ALGC
        END DO
      END DO
      END

************************************************************************
**       S U B R O U T I N E   D I S S O L V E D   O X Y G E N        **
************************************************************************

      SUBROUTINE DISSOLVED_OXYGEN
      INCLUDE   'w2.inc'
      REAL    NH3D, KBOD, ICETH
      LOGICAL ICE, ICE_CALC
      integer wscdp
      double precision t1,t2
      COMMON /GLOBLC/ JB, JC, IU, ID, KT, DLT, KB(IMP)
      COMMON /DKMLTC/ A1(KMC,IMC),    A2(KMC,IMC), A3(KMC,IMC)
      COMMON /GLBLCC/ PALT,   ALGDET, O2LIM, WIND, wscdp, wsc(ndp)
      COMMON /OXYGNC/ O2ORG,  O2ALG,  O2NH3, O2RESP, o2xfact
      COMMON /CBODC/  KBOD,   TBOD,   RBOD
      COMMON /GEOMHC/ EL(KMP),        H(KMP),        HKT1(IMP),
     .                HKT2(IMP)
      COMMON /ICEC/   ICE(IMP),       ICETH(IMP),    ICE_CALC
      COMMON /TEMPC/  T1(KMP,IMP),    T2(KMP,IMP)
      COMMON /DKORGC/ DETD(KMC,IMC),  ORGD(KMC,IMC)
      COMMON /DKSEDC/ SEDD(KMC,IMC),  SO2D(KMC,IMC), SOD(IMP)
      COMMON /DKNITC/ NH3D(KMC,IMC),  NO3D(KMC,IMC)
      COMMON /PHYTGC/ AGR(KMC,IMC),   ARR(KMC,IMC),  AMR(KMC,IMC),
     .                AER(KMC,IMC)
      COMMON /DKBODC/ CBODD(KMC,IMC)
      common /zoopc/ zmax,zmort,zeffic,pref1,pref2,zresp,zoomin,zs2p,
     .               zmu(kmc,imc),talgae(kmc,imc),zrt(kmc,imc),
     .               zmt(kmc,imc)

      SATO(X) = EXP(7.7117-1.31403*(LOG(X+45.93)))*PALT
      O2EX=o2xfact*2.04E-9/((200.0-60.0*SQRT(MIN(WIND,11.0)))*1.E-6)
      DO I=IU,ID
        DOSS(KT,I) = 0.0
        DO K=KT,KB(I)
          ALGP      = (O2ALG*AGR(K,I)-O2RESP*ARR(K,I))*ALGAE(K,I)
          zoopp=o2resp*zrt(k,i)*zoop(k,i)
          DOSS(K,I) = ALGP-O2NH3*NH3D(K,I)-O2ORG*(DETD(K,I)+SEDD(K,I))
     .                -SO2D(K,I)*A3(K,I)-O2ORG*ORGD(K,I)-CBODD(K,I)
     .                *CBOD(K,I)*RBOD-zoopp
        END DO
        ttt=T1(KT,I)
        SATDO = SATO(ttt)
        IF (.NOT.ICE(I)) DOSS(KT,I) = DOSS(KT,I)+(SATDO-DO(KT,I))
     .                                *O2EX/HKT2(I)
      END DO
      END

************************************************************************
**               S U B R O U T I N E   S E D I M E N T                **
************************************************************************

      SUBROUTINE SEDIMENT
      INCLUDE   'w2.inc'
      COMMON /GLOBLC/ JB, JC, IU, ID, KT, DLT, KB(IMP)
      COMMON /GEOMHC/ EL(KMP),    H(KMP),      HKT1(IMP),   HKT2(IMP)
      COMMON /GEOMBC/ B(KMP,IMP), BKT(IMP),    BH(KMP,IMP), BHKT1(IMP),
     .                BHKT2(IMP), BHRKT1(IMP)
      COMMON /SETLC2/ SSETL, DSETL,  ASETL, FESETL
      COMMON /DKSEDC/ SEDD(KMC,IMC), SO2D(KMC,IMC), SOD(IMP)

      DO I=IU,ID
        SETTLE       = (ASETL*ALGAE(KT,I)+DSETL*DETRIT(KT,I))*DLT
     .                 /HKT2(I)*(1.0-B(KT+1,I)/BKT(I))
        SEDMNT(KT,I) = MAX(SEDMNT(KT,I)+SETTLE-SEDD(KT,I)*DLT,0.0)
        DO K=KT+1,KB(I)-1
          SETTLE      = (ASETL*ALGAE(K,I)+DSETL*DETRIT(K,I))*DLT/H(K)
     .                  *(1.0-B(K+1,I)/B(K,I))
          SEDMNT(K,I) = MAX(SEDMNT(K,I)+SETTLE-SEDD(K,I)*DLT,0.0)
        END DO
        SETTLE          = (ASETL*ALGAE(KB(I),I)+DSETL*DETRIT(KB(I),I))
     .                    *DLT/H(K)
        SEDMNT(KB(I),I) = MAX(SEDMNT(KB(I),I)+SETTLE-SEDD(KB(I),I)*DLT,
     .                        0.0)
      END DO
      END

************************************************************************
**       S U B R O U T I N E   I N O R G A N I C   C A R B O N        **
************************************************************************

      SUBROUTINE INORGANIC_CARBON
      INCLUDE   'w2.inc'
      LOGICAL ICE, ICE_CALC
      integer wscdp
      double precision t1,t2
      COMMON /CARBNC/ CO2REL, BIOC, obioc
      COMMON /GLOBLC/ JB, JC, IU, ID, KT, DLT, KB(IMP)
      COMMON /GLBLCC/ PALT, ALGDET, O2LIM, WIND, wscdp, wsc(ndp)
      COMMON /GEOMHC/ EL(KMP),       H(KMP),        HKT1(IMP),
     .                HKT2(IMP)
      COMMON /ICEC/   ICE(IMP),      ICETH(IMP),    ICE_CALC
      COMMON /TEMPC/  T1(KMP,IMP),   T2(KMP,IMP)
      COMMON /DKORGC/ DETD(KMC,IMC), ORGD(KMC,IMC)
      COMMON /DKMLTC/ A1(KMC,IMC),    A2(KMC,IMC), A3(KMC,IMC)
      COMMON /DKSEDC/ SEDD(KMC,IMC), SO2D(KMC,IMC), SOD(IMP)
      COMMON /PHYTGC/ AGR(KMC,IMC),  ARR(KMC,IMC),  AMR(KMC,IMC),
     .                AER(KMC,IMC)
      common /zoopc/ zmax,zmort,zeffic,pref1,pref2,zresp,zoomin,zs2p,
     .               zmu(kmc,imc),talgae(kmc,imc),zrt(kmc,imc),
     .               zmt(kmc,imc)
      COMMON /OXYGNC/ O2ORG,  O2ALG,  O2NH3, O2RESP, o2xfact

      CO2EX = o2xfact*1.63E-9/
     .      ((200.0-60.0*SQRT(MIN(WIND,11.0)))*1.0E-6)
      DO I=IU,ID
        CSS(KT,I) = 0.0
        DO K=KT,KB(I)
          ALGP     = (ARR(K,I)-AGR(K,I))*ALGAE(K,I)
          zoopp=zrt(k,i)*zoop(k,i)
        CSS(K,I) = BIOC*(ALGP+zoopp+detd(k,i))+
     .             obioc*(ORGD(K,I)+SEDD(K,I))
     .               +CO2REL*SO2D(K,I)*A3(K,I)
        END DO
        IF (.NOT.ICE(I)) CSS(KT,I) = CSS(KT,I)+CO2EX*(0.286
     .                               *EXP(-0.0314*(T1(KT,I))*PALT)
     .                               -co2(KT,I))/HKT2(I)

      END DO
      END

************************************************************************
**                S U B R O U T I N E   P H   C O 2                   **
************************************************************************

      SUBROUTINE PH_CO2
      INCLUDE   'w2.inc'
      REAL    KW, K1, K2, INCR
      LOGICAL FRESH_WATER, SALT_WATER, SUSP_SOLIDS
      double precision t1,t2
      COMMON /GLOBLC/ JB, JC, IU,  ID, KT, DLT, KB(IMP)
      COMMON /TEMPC/  T1(KMP,IMP), T2(KMP,IMP)
      COMMON /DNSPHC/ FRESH_WATER, SALT_WATER,  SUSP_SOLIDS

C NOTE INPUT UNITS OF ALKAL IS MG/L AS CACO3, INPUT UNITS OF CARBON
C IS MG/L AS C
      DO I=IU,ID
        DO K=KT,KB(I)
          CAR1 = CARBON(K,I)/12000.0
          ALK1 = ALKAL(K,I)/5.0E+04
          TEM2 = T1(K,I)+273.15

********* Ionic strength

          IF (FRESH_WATER) S2 = 2.5E-05*DISS(K,I)
          IF (SALT_WATER)  S2 = 0.00147+0.019885*DISS(K,I)+0.000038
     .                          *DISS(K,I)*DISS(K,I)

********* Debye-Huckel terms and activity coefficients

c the following line was commented out in my version
C          IF (ABS(S1-S2).GE.0.1) THEN
            SQRS2  = SQRT(S2)
            DH1    = -.5085*SQRS2/(1.0+1.3124*SQRS2)+4.745694E-03
     .               +4.160762E-02*S2-9.284843E-03*S2*S2
            DH2    = -2.0340*SQRS2/(1.0+1.4765*SQRS2)+1.205665E-02
     .               +9.715745E-02*S2-2.067746E-02*S2*S2
            H2CO31 = 10.0**(0.0755*S2)
            HCO31  = 10.0**DH1
            CO31   = 10.0**DH2
            OH     = HCO31
c the following 2 lines were commented out in my version
C            S1     = S2
C          END IF

********* Temperature adjustment

c the following line was commented out in my version
C          IF (ABS(TEM1-TEM2).GE.1.0) THEN
            KW   = 10.0**(-5242.39/TEM2+35.3944-0.00835*TEM2-11.8261
     .             *LOG10(TEM2))/OH
            K1   = 10.0**(-3404.71/TEM2+14.8435-0.032786*TEM2)*H2CO31
     .             /HCO31
            K2   = 10.0**(-2902.39/TEM2+6.498-0.02379*TEM2)*HCO31/CO31
c the following 2 lines were commented out
C            TEM1 = TEM2
C          END IF

********* pH evaluation

          PH1 = -PH(K,I)-2.1
          IF (PH(K,I).LE.0.0) PH1 = -14.0
          INCR = 10.0
          DO J=1,3
            FF   = 1.0
            INCR = INCR/10.0
            ITER = 0
            DO WHILE (FF.GT.0.0.AND.ITER.LT.12)
              PH1    = PH1+INCR
              HION   = 10.0**PH1
              BICAR1 = CAR1*K1/(K1+K1*K2+HION)
              FF     = BICAR1*(HION+2.0*K2)/HION+KW/HION-ALK1-HION
              ITER   = ITER+1
            END DO
            PH1 = PH1-INCR
          END DO

********* pH, carbon dioxide, bicarbonate, and carbonate concentrations

          PH(K,I)   = -PH1
          CO2(K,I)  = CARBON(K,I)/(1.0+K1/HION+K1*K2/(HION*HION))
          HCO3(K,I) = CARBON(K,I)/(1.0+HION/K1+K2/HION)
          CO3(K,I)  = CARBON(K,I)/((HION*HION)/(K1*K2)+HION/K2+1.0)
c reassign variable names
         c1(k,i,16)=ph(k,i)
         c1(k,i,17)=co2(k,i)
         c1(k,i,18)=hco3(k,i)
         c1(k,i,19)=co3(k,i)
        END DO
      END DO
      END

************************************************************************
**                  S U B R O U T I N E    I R O N                    **
************************************************************************

      SUBROUTINE IRON
        INCLUDE   'w2.inc'
        REAL    NETSET
        COMMON /IRONC/  FEREL
        COMMON /GLOBLC/ JB, JC, IU, ID, KT, DLT, KB(IMP)
        COMMON /SETLC2/ SSETL,  DSETL,  ASETL,   FESETL
        COMMON /GEOMHC/ EL(KMP),       H(KMP),        HKT1(IMP),
     .                  HKT2(IMP)
        COMMON /DKSEDC/ SEDD(KMC,IMC), SO2D(KMC,IMC), SOD(IMP)
        COMMON /DKMLTC/ A1(KMC,IMC),   A2(KMC,IMC), A3(KMC,IMC)

        DO I=IU,ID
          NETSET     = -FESETL*FE(KT,I)*A1(KT,I)/HKT2(I)
          SEDREL     = FEREL*SO2D(KT,I)*A2(KT,I)
          FESS(KT,I) = SEDREL+NETSET
          DO K=KT+1,KB(I)
            NETSET    = FESETL*(FE(K-1,I)*A1(K-1,I)-FE(K,I)
     .                  *A1(K,I))/H(K)
            SEDREL    = FEREL*SO2D(K,I)*A2(K,I)
            FESS(K,I) = SEDREL+NETSET
          END DO
        END DO
      END

************************************************************************
**  S U B R O U T I N E   B I O C H E M I C A L   O 2   D E M A N D   **
************************************************************************

      SUBROUTINE BIOCHEMICAL_O2_DEMAND
        INCLUDE   'w2.inc'
        COMMON /GLOBLC/ JB, JC, IU, ID, KT, DLT, KB(IMP)
        COMMON /DKBODC/ CBODD(KMC,IMC)

        DO I=IU,ID
          DO K=KT,KB(I)
            CBODSS(K,I) = -CBODD(K,I)*CBOD(K,I)
          END DO
        END DO
      END
c**********************************************************
c**           SUBROUTINE ZOOPLANKTON                     **
C**********************************************************
      subroutine zooplankton
      include 'w2.inc'
      COMMON /GLOBLC/ JB, JC, IU, ID, KT, DLT, KB(IMP)
      common /zoopc/ zmax,zmort,zeffic,pref1,pref2,zresp,zoomin,zs2p,
     .               zmu(kmc,imc),talgae(kmc,imc),zrt(kmc,imc),
     .               zmt(kmc,imc)
c note that zmu, talgae were calculated in subroutine algae
c and that zmt and zrt were calculated in subroutine detritus
          DO 10010 I=IU,ID
          DO 10000 K=KT,KB(I)
            CSSK(K,I,JC) = (zmu(k,i)*zeffic-zrt(k,i)-zmt(k,i))
     .                     *zoop(k,i)
10000     CONTINUE
10010   CONTINUE
      END
***********************************************************************
***********************************************************************

       subroutine qculv(el1,el2,aqout,aday,mup,mdn)

c     This subroutine is to be included in the CE-QUAL-W2 Army Corps of
c     Engineers model to calculate flow through culverts on the upper slough.

c     variables:
c       el1   - upper reach surface elevation
c       el2   - lower reach surface elevation
c       mup   - upper reach branch number
c       mdn   - lower reach branch number
c       dia   - diameter of culvert
c      clen   - length of culver
c      upie   - upper reach invert elevation
c      dnie   - lower reach invert elevation
c        ec   - contraction coefficient for entrance
c      fman   - Manning's friction coefficient
c      closs  - minor loss coefficient
c     width   - width of weir or width of box culvert
c       hdi   - height of drop inlet
c        wc   - weir coefficient
c        cm   - Manning's unit coefficient (cm=1.486, English; cm=1.0, metric)
c      qout   - culvert flow
c         g   - gravitational constant
c       eps   - threshold

      parameter(pi=3.14159265359)
      implicit double precision(b-h,o-z)
      dimension adia(50),aclen(50),aupie(50),adnie(50),aec(50),
     .          afman(50),acloss(50),awidth(50),ahdi(50),awc(50)
      common/culvert/dia,fman,fall,dist,ec,hdn,closs,g,cm,hdif,upie,
     .               dnie,wc,width,hdi
      common/prison/aqbox,aqppipe

      g=9.81
      day=aday
      cm=1.0
      tqout=0.0

      if(nit.eq.1)go to 900
      nit=1
      xday=0.
      dday=0.
      open(99,file='culvert.npt',STATUS='OLD')
      open(98,file='mcdd1.npt',status='old')
      open(97,file='wl_ls92.npt',status='old')
      read(97,688)
      read(97,688)
      read(97,688)
      read(98,688)
      read(98,688)
      read(98,688)
      read(99,688)
      read(99,688)
      read(99,688)
      read(99,788,end=900)(dummy,adia(i),aclen(i),aupie(i),adnie(i),
     .       aec(i),afman(i),acloss(i),awidth(i),ahdi(i),awc(i),i=1,39)
      close(99)
688   format(1X)
788   format(11f8.3)
900   continue

      if(mup.eq.2.and.mdn.eq.1)then
        jc=1
        i1=1
      end if
      if(mup.eq.3.and.mdn.eq.2)then
        jc=1
        i1=2
      end if
      if(mup.eq.4.and.mdn.eq.3)then
        jc=2
        i1=3
        i2=4
      end if
      if(mup.eq.5.and.mdn.eq.4)then
        jc=2
        i1=5
        i2=6
      end if
      if(mup.eq.6.and.mdn.eq.5)then
        jc=2
        i1=7
        i2=8
      end if
      if(mup.eq.7.and.mdn.eq.6)then
        jc=2
        i1=9
        i2=10
      end if
      if(mup.eq.8.and.mdn.eq.7)then
        jc=3
        i1=11
        i2=12
        i3=13
      end if
      if(mup.eq.10.and.mdn.eq.9)then
        jc=1
        i1=14
      end if
      if(mup.eq.11.and.mdn.eq.10)then
        jc=1
        i1=15
      end if
      if(mup.eq.12.and.mdn.eq.11)then
        jc=1
        i1=16
      end if
      if(mup.eq.14.and.mdn.eq.13)then
        jc=3
        i1=17
        i2=18
        i3=19
      end if
      if(mup.eq.15.and.mdn.eq.14)then
        jc=1
        i1=20
      end if
      if(mup.eq.16.and.mdn.eq.15)then
        jc=1
        i1=21
      end if
      if(mup.eq.17.and.mdn.eq.14)then
        jc=1
        i1=22
      end if
      if(mup.eq.18.and.mdn.eq.17)then
        jc=1
        i1=23
      end if
      if(mup.eq.19.and.mdn.eq.18)then
        jc=1
        i1=24
      end if
      if(mup.eq.20.and.mdn.eq.19)then
        jc=1
        i1=25
      end if
      if(mup.eq.21.and.mdn.eq.20)then
        jc=1
        i1=26
      end if
      if(mup.eq.22.and.mdn.eq.21)then
        jc=1
        i1=27
      end if
      if(mup.eq.23.and.mdn.eq.22)then
        jc=2
        i1=28
        i2=29
      end if
      if(mup.eq.24.and.mdn.eq.23)then
        jc=2
        i1=30
        i2=31
      end if
      if(mup.eq.26.and.mdn.eq.24)then
        jc=1
        i1=32
      end if
      if(mup.eq.27.and.mdn.eq.26)then
        jc=1
        i1=33
      end if
      if(mup.eq.29.and.mdn.eq.28)then
        jc=2
        i1=34
        i2=39
      end if
      if(mup.eq.30.and.mdn.eq.29)then
        jc=2
        i1=35
        i2=36
      end if
      if(mup.eq.31.and.mdn.eq.30)then
        jc=2
        i1=37
        i2=38
      end if
      if(mup.eq.1.and.mdn.eq.100)then
        jc=1
        i1=40
      end if

      do 525 i=1,jc
        if(i.eq.1)ic=i1
        if(i.eq.2)ic=i2
        if(i.eq.3)ic=i3

c   reading in input data for mcdd#1
        if(ic.eq.40)then
           IF (DAY.GE.XDAY) THEN
545          IF (DAY.GE.XDAY) THEN
               READ (98,791)XDAY
               GO TO 545
             END IF
             BACKSPACE 98
             BACKSPACE 98
      read(98,289)dummy,adia(ic),aclen(ic),aupie(ic),adnie(ic),
     .       aec(ic),afman(ic),acloss(ic),awidth(ic),ahdi(ic),awc(ic)
           END IF
c   reading lower slough water levels
           IF (DAY.GE.dDAY) THEN
345          IF (DAY.GE.dDAY) THEN
               READ (97,791)dDAY
               GO TO 345
             END IF
             BACKSPACE 97
             BACKSPACE 97
             read(97,589)dummy,ells
           END IF
589        FORMAT(2F8.3)
289        format(11f8.3)
791        format(f8.3)
      end if

        dia=adia(ic)
        clen=aclen(ic)
        upie=aupie(ic)
        dnie=adnie(ic)
        ec=aec(ic)
        fman=afman(ic)
        closs=acloss(ic)
        width=awidth(ic)
        hdi=ahdi(ic)
        wc=awc(ic)
        if(ic.eq.40)el2=ells
        day=aday

      hdif=abs(el1-el2)
      fall=upie-dnie
      dist=clen
      hie=max(upie,dnie)

c     test if either upstream or downstream surface elevations exceed
c     the highest (upstream) I.E., then q=0.
      eps=1.e-9
      if((hie+eps).ge.el1.and.(hie+eps).ge.el2)then
        qout=0.
        go to 450
      end if

c     when upstream and downstream surface elevations are approximately
c     equal, q=0.

      if(abs(hdif).lt.eps)then
        qout=0.
        go to 450
      end if

c     special algorithm for flow at riser and mcdd#1

      if(ic.eq.40.or.ic.eq.21)then
        fall=(upie-hdi)-dnie
        qout=riser(el1,el2,ic)
        go to 450
      end if

c     special algorithm for flow through box culvert

      if(ic.eq.34)then
        qout=boxculv(el1,el2)
        aqbox=qout/(0.3048**3)
        go to 450
      end if

c     decision sequence for circular culverts

      if(el1.gt.el2)then
        head=el1-upie
        hdn=el2-dnie
        if(upie.ge.dnie)then
          go to 225
        else
          go to 325
        end if
      end if

      if(el2.gt.el1)then
        head=el2-dnie
        hdn=el1-upie
        fall=-fall
        if(upie.ge.dnie)then
          go to 325
        else
          go to 225
        end if
      end if
c     if direction of flow is in the downhill slope direction

225     if(head.ge.dia.and.hdn.ge.dia)then
          qout=type4(head)
          go to 450
        end if
        if(head.le.dia.and.hdn.lt.dia)then
          if(hdn.gt.0.0)then
            qout=type7(head)
            go to 450
          else
            hdn=0.0
            if(el1.gt.el2)then
              hdif=el1-dnie
            else
              hdif=el2-upie
            end if
            qout=type7(head)
            go to 450
          end if
        end if


c     if outlet is not submerged and inlet head is > 1.2*diameter,
c     flow is either type 5 or type 6

        if(head.gt.dia.and.hdn.lt.dia)then
          qt5=type5(head)
          sf=10.3*(fman*qt5)**2/(cm**2*dia**5.33)
          cl= (((qt5/(3.14*dia**2/4.))**2/(2.*g))*(1.-1./ec**2)
     .         + dia - 0.74*dia)/(fall/dist - sf)
          if(cl.le.dist)then
            qout=qt5
            go to 450
          else
            qout=type6(head)
            go to 450
          end if
        end if

c     when the direction of flow is in the uphill slope direction

325     if(head.ge.dia.and.hdn.ge.dia)then
          qout=type4(head)
          go to 450
        else if(head.ge.dia.and.(hdn+abs(fall)).gt.dia)then
c error next line
c          dist=(hdn+abs(fall)-dnie-dia)/abs(fall/clen)
           theta=acos(abs(fall)/clen)
           dist=clen-(dia-hdn)/(tan((pi/2.)-theta))
          qout=type4(head)
          go to 450
        else if(head.ge.dia)then
           qout=type6(head)
          go to 450
        else if(head.lt.dia.and.hdn.lt.0.0)then
          hdn=0.0
          if(el1.gt.el2)then
            hdif=el1-dnie
          else
            hdif=el2-upie
          end if
          qout=type7(head)
          go to 450
        else
          qout=type7(head)
          go to 450
        end if

450   continue
      if(ic.eq.39)aqppipe=qout/(0.3048**3)
      tqout=tqout+qout
525   continue
      qout=tqout
      if(ic.eq.40)qout=5.*qout
      if(el2.gt.el1)qout=-qout
      aqout=qout
      return

      end

c     this function computes the flow of a culvert that has a steep
c     slope and is inlet controlled and has a unsubmerged entrance

      function type7(head)
      implicit double precision(b-h,o-z)
      common/culvert/dia,fman,fall,dist,ec,hdn,closs,g,cm,hdif,upie,
     .               dnie,wc,width,hdi

      depth=head
      bac = barea(depth,dia)

      type7 = ec *bac * sqrt((2.*g)*(hdif))

      end

c     this function computes the flow of a culvert that has a steep
c     slope and is inlet controlled and has a unsubmerged entrance

      function type1(head)
      implicit double precision(b-h,o-z)
      external t1func

      x1=1.e-4
      x2=10.
      tol=0.001

      type1 = zbrent2(t1func,x1,x2,tol,head,1)

      end

c     when the function T1FUNC equals zero, the flow in the culvert
c     is type1 flow

      function t1func(flow,head)
      implicit double precision(b-h,o-z)
      common/culvert/dia,fman,fall,dist,ec,hdn,closs,g,cm,hdif,upie,
     .               dnie,wc,width,hdi

      cd = depthcrit(flow)
      bac = barea(cd,dia)

      t1func = flow**2 - (ec * bac)**2 * ((2.*g)*(head-cd))

      end

c     this function computes the flow of a culvert that has a mild
c     slope and is outlet controlled with an unsubmerged entrance

      function type2(head)
      implicit double precision(b-h,o-z)
      external t2func
      common/culvert/dia,fman,fall,dist,ec,hdn,closs,g,cm,hdif,upie,
     .               dnie,wc,width,hdi

      x1=0.1e-3
      x2=0.33*dia*cm*sqrt(hdif/dist)/fman
      tol=0.00001


      type2 = zbrent2(t2func,x1,x2,tol,head,2)

      end

c     when the function T2FUNC equals zero, the flow in the culvert
c     is type2 flow

      function t2func(flow,head)
      implicit double precision(b-h,o-z)
      common/culvert/dia,fman,fall,dist,ec,hdn,closs,g,cm,hdif,upie,
     .               dnie,wc,width,hdi

      cd = depthcrit(flow)
      bac = barea(cd,dia)
      fl = floss(flow,2)
      if(fl.gt.hdif)fl=hdif
c      fl =fall

      t2func = flow**2. - (ec * bac)**2. * ((2.*g)*(head+fall-fl-cd))

      end

c     this function computes the flow of a culvert that has a mild
c     slope and is outlet controlled with an unsubmerged entrance
c     and whose tailwater height is greater than the critical depth

      function type3(head)
      implicit double precision(b-h,o-z)
      common/culvert/dia,fman,fall,dist,ec,hdn,closs,g,cm,hdif,upie,
     .               dnie,wc,width,hdi
      external t3func

      x1=0.1e-3
      x2=0.33*dia*cm*sqrt(hdif/dist)/fman
      tol=0.00001


      type3 = zbrent2(t3func,x1,x2,tol,head,3)

      end

c     when the function T3FUNC equals zero, the flow in the culvert
c     is type 3 flow

      function t3func(flow,head)
      implicit double precision(b-h,o-z)
      common/culvert/dia,fman,fall,dist,ec,hdn,closs,g,cm,hdif,upie,
     .               dnie,wc,width,hdi

      ba3 = barea(hdn,dia)
      fl = floss(flow,3)
      if(fl.gt.hdif)fl=hdif
c      fl=fall

      t3func = flow**2 - (ec * ba3)**2 * ((2.*g)*(head+fall
     .         -fl-hdn))
      end

c     this function computes the flow of a culvert that has a
c     submerged entrance and whose tailwater height is greater than
c     its diameter

      function type4(head)
      implicit double precision(b-h,o-z)
      common/culvert/dia,fman,fall,dist,ec,hdn,closs,g,cm,hdif,upie,
     .               dnie,wc,width,hdi

      parameter(pi=3.14159265359)

      bao = pi*dia**2/4.
      ro = bao/(pi*dia)
      fl = (2.*g*dist*fman**2)/(cm**2*ro**(4./3.))

c      type4=ec * bao * sqrt(2*g*(head+fall-hdn)/(1.+closs+fl))
      type4=ec * bao * sqrt(2*g*hdif/(1.+closs+fl))

      end

c     this function computes the flow of a culvert that flows full and
c     has a submerged entrance and a tailwater height which is less
c     than its diameter

      function type5(head)
      implicit double precision(b-h,o-z)
      common/culvert/dia,fman,fall,dist,ec,hdn,closs,g,cm,hdif,upie,
     .               dnie,wc,width,hdi

      parameter(pi=3.14159265359)

      bao = pi*dia**2/4.
      ro = bao/(pi*dia)
      fl = (2.*g*dist*fman**2)/(cm**2*ro**(4./3.))

      type5=ec * bao * sqrt(2*g*(head+fall-dia)/(1.+closs+fl))

      end

c     this function computes the flow of a culvert that flows part full and
c     has a submerged entrance and a tailwater height which is less
c     than its diameter

      function type6(head)
      implicit double precision(b-h,o-z)
      common/culvert/dia,fman,fall,dist,ec,hdn,closs,g,cm,hdif,upie,
     .               dnie,wc,width,hdi
      parameter(pi=3.14159265359)

      bao = pi*dia**2/4.

      type6 = bao*ec*sqrt(2.*g*head)

      end

c     the function DEPTHCRIT calculates the critical depth using
c     the root finding function zbrent1

      function depthcrit(flow)
      implicit double precision(b-h,o-z)
      common/culvert/dia,fman,fall,dist,ec,hdn,closs,g,cm,hdif,upie,
     .               dnie,wc,width,hdi
      external cdfunc

      x1=dia/50.
      x2=dia
      tol=.001

      depthcrit=zbrent1(cdfunc,x1,x2,tol,flow,7)


      end


c     When the function cdfunc equals zero, the depth is at critical
c     depth

      function cdfunc(depth,flow)
      implicit double precision(b-h,o-z)
      common/culvert/dia,fman,fall,dist,ec,hdn,closs,g,cm,hdif,upie,
     .               dnie,wc,width,hdi



      cdfunc=(flow**2.*twidth(depth,dia))/(barea(depth,dia)**3.*g)-1.

      end


c     the function DEPTHNORM calculates the normal depth using
c     the root finding function zbrent1

      function depthnorm(flow)
      implicit double precision(b-h,o-z)
      common/culvert/dia,fman,fall,dist,ec,hdn,closs,g,cm,hdif,upie,
     .               dnie,wc,width,hdi
      external dnfunc


       x1=dia/50.
       x2=dia*.94
       tol=.001

        depthnorm = zbrent1(dnfunc,x1,x2,tol,flow,8)

      end

c   when the function DNFUNC equals zero, the depth of the pipe
c   is the normal depth

      function dnfunc(depth,flow)
      implicit double precision(b-h,o-z)
      common/culvert/dia,fman,fall,dist,ec,hdn,closs,g,cm,hdif,upie,
     .               dnie,wc,width,hdi

      dnfunc = flow-cm*(barea(depth,dia)**(5./3.))*sqrt(hdif/dist)/
     .                 (fman*(wetper(depth,dia))**(2./3.))

      end


c     this function calculates the friction head loss in the culvert

      function floss(flow,num)
      implicit double precision(b-h,o-z)
      common/culvert/dia,fman,fall,dist,ec,hdn,closs,g,cm,hdif,upie,
     .               dnie,wc,width,hdi

      if(num.eq.3)depth=hdn
      if(num.eq.2)depth=depthcrit(flow)

      floss = (flow * fman)**2 * wetper(depth,dia)**(4./3.) * dist/
     .        (cm**2 * barea(depth,dia)**(10./3.))

      end

c   this function calculates the cross-sectional area occupied
c   by water in a circular pipe given the water depth and pipe diameter

      function barea(depth,dia)
      implicit double precision(b-h,o-z)

      parameter(pi=3.14159265359)
c note that the argument to asin has to between 1 and -1 only
c hence if depth>dia then just calculate area as pipe area
        if(depth.lt.dia)then
        barea= (depth-dia/2.)*sqrt(depth*dia-depth**2.)+(dia**2./4.)*
     .         asin((2./dia)*(depth-dia/2.)) + (pi*dia**2.)/8.
        else
        barea=(pi*dia**2.)/4.
        end if

      end

c     this function calculates the width of the liquid cross-section
c     at the liquid surface within a circular  pipe

      function twidth(depth,dia)
      implicit double precision(b-h,o-z)

      if(depth.lt.dia)then
      twidth = 2.*sqrt((dia*depth)-depth**2.)
      else
      twidth=0.0
      end if

      end

c   this function calculates the wetted perimeter
c   in a circular pipe given the water depth and pipe diameter

      function wetper(depth,dia)
      implicit double precision(b-h,o-z)
      parameter(pi=3.14159265359)
c asin arg only valid from 1 to -1
c hence note changes
        if(depth.lt.dia)then
        wetper= dia * (asin((2./dia)*(depth-dia/2.)) + pi/2.)
        else
        wetper=pi*dia
        end if

      end


c     The function zbrent1 finds the root of the function FUNC known to
c     lie between x1 and x2. The root is returned as zbrent1 when its
c     accuracy is known within TOL.

c      From:
c     "NUMERICAL RECIPES, THE ART OF SCIENTIFIC COMPUTING", Press, et. al.,
c      Cambridge University Press, pages 243-254.

      function zbrent1(func,x1,x2,tol,barg,num)
      implicit double precision(b-h,o-z)
      external func


      parameter(factor=0.1,ntry=50,itmax=100,eps=3.e-8)

c     root bracketing algorithm used to ensure x1 and x2  bracket the root

      if(x1.eq.x2)pause 'you have to guess an initial range'
      f1 = func(x1,barg)
      f2 = func(x2,barg)

c     special bracketing for normal and critical depth functions

      if(num.ge.7.and.f1.le.0.)then
        do 35 i=1,40
          x1=x1/10.
          f1 = func(x1,barg)
          if(f1.gt.0.)go to 79
35        continue
          if(num.eq.7)pause 'could not make cdfunc(x1) > zero'
          if(num.eq.8)pause 'could not make ndfunc(x1) > 0'
      end if
79        continue


      do 33 j=1,ntry
        if(f1*f2.lt.0.)go to 44
        if(abs(f1).lt.abs(f2))then
          x1=x1+factor*(x1-x2)
          f1 = func(x1,barg)
        else
          x2=x2+factor*(x2-x1)
          f2 = func(x2,barg)
        end if
33    continue
      pause 'could not bracket root'
44    continue

c     Brent's root finding algorithm

      ba=x1
      b=x2
      fa = func(ba,barg)
      fb = func(b,barg)

      if(fb*fa.gt.0.)pause 'root must be bracketed for zbrent1'
      fc=fb
      do 11 iter=1,itmax
        if(fb*fc.gt.0.)then
          c=ba
          fc=fa
          d=b-ba
          e=d
        end if
        if(abs(fc).lt.abs(fb))then
          ba=b
          b=c
          c=ba
          fa=fb
          fb=fc
          fc=fa
        end if
c     convergence check
        tol1=2.*eps*abs(b)+0.5*tol
        xm=.5*(c-b)
        if(abs(xm).le.tol1.or.fb.eq.0.)then
         zbrent1=b
         go to 22
        end if
        if(abs(e).ge.tol1.and.abs(fa).gt.abs(fb))then
          s=fb/fa
          if(ba.eq.c)then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-ba)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          end if
          if(p.gt.0.)q=-q
          p=abs(p)
          if(2.*p.lt.min(3.*xm*q-abs(tol1*q),abs(e*q)))then
            e=d
            d=p/q
          else
            d=xm
            e=d
          end if
        else
          d=xm
          e=d
        end if
        ba=b
        fa=fb
        if(abs(d).gt.tol1)then
          b=b+d
        else
          b=b+sign(tol1,xm)
        end if
        fb= func(b,barg)
11    continue
      pause 'zbrent1 exceeding maximum number of iterations'
      zbrent1=b
22    continue

      end



c     The function zbrent finds the root of the function FUNC known to
c     lie between x1 and x2. The root is returned as ZBRENT when its
c     accuracy is known within TOL.

c      From:
c     "NUMERICAL RECIPES, THE ART OF SCIENTIFIC COMPUTING", Press, et. al.,
c      Cambridge University Press, pages 243-254.

      function zbrent2(func,x1,x2,tol,barg,num)
      implicit double precision(b-h,o-z)
      external func


      parameter(ntry=50,itmax=100,eps=3.e-8)

c     root bracketing algorithm used to ensure x1 and x2  bracket the root

      factor = 2.00

      if(x1.eq.x2)pause 'you have to guess an initial range'
      f1 = func(x1,barg)
      f2 = func(x2,barg)

c     setting lower bracket

      if(num.le.5.and.f1.gt.0.)then
        do 25 i=1,100
          x1=x1/10.
          f1 = func(x1,barg)
           if(f1.lt.0.)go to 69
25        continue
          pause 'could not make t#func(x1) less than zero'
      end if
69        continue

c     setting upper bracket

      if(num.gt.4)factor=2.
      if(num.le.5.and.f2.lt.0.)then
        do 45 i=1,100
          x2=x2*factor
          f2 = func(x2,barg)
           if(f2.gt.0.)go to 89
45        continue
          pause 'could not make typef(x2) greater than zero'
      end if
89        continue



c     Brent's root finding algorithm

      ba=x1
      b=x2
      fa = func(ba,barg)
      fb = func(b,barg)

      if(fb*fa.gt.0.)pause 'root must be bracketed for zbrent2'
      fc=fb
      do 11 iter=1,itmax
        if(fb*fc.gt.0.)then
          c=ba
          fc=fa
          d=b-ba
          e=d
        end if
        if(abs(fc).lt.abs(fb))then
          ba=b
          b=c
          c=ba
          fa=fb
          fb=fc
          fc=fa
        end if
c     convergence check
        tol1=2.*eps*abs(b)+0.5*tol
        xm=.5*(c-b)
        if(abs(xm).le.tol1.or.fb.eq.0.)then
         zbrent2=b
         go to 22
        end if
        if(abs(e).ge.tol1.and.abs(fa).gt.abs(fb))then
          s=fb/fa
          if(ba.eq.c)then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-ba)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          end if
          if(p.gt.0.)q=-q
          p=abs(p)
          if(2.*p.lt.min(3.*xm*q-abs(tol1*q),abs(e*q)))then
            e=d
            d=p/q
          else
            d=xm
            e=d
          end if
        else
          d=xm
          e=d
        end if
        ba=b
        fa=fb
        if(abs(d).gt.tol1)then
          b=b+d
        else
          b=b+sign(tol1,xm)
        end if
        fb= func(b,barg)
11    continue
      pause 'zbrent2 exceeding maximum number of iterations'
      zbrent2=b
22    continue

      end

c     algorithm for calculating flow at risers and mcdd#1

      function riser(el1,el2,ic)
      implicit double precision(b-h,o-z)
      common/culvert/dia,fman,fall,dist,ec,hdn,closs,g,cm,hdif,upie,
     .               dnie,wc,width,hdi

      parameter(pi=3.14159265359)

      fall=(upie-hdi)-dnie
      if(ic.eq.40)width=width/5.

c     if water level of lower reach exceeds that of the upper reach

       if(el2.gt.el1)then
         if(ic.eq.40)then
           riser=0.0
           go to 550
         end if
         head=el2-dnie
         hdn=el1-(upie-hdi)
         riser=type4(head)
         go to 550
       end if

c     if surface level of lower slough exceeds the upper slough weir
c     crest elevation but is still less than the surface level of the
c     upper slough
      if(el2.ge.upie.and.el1.gt.el2)then
        head=el1-(upie-hdi)
        hdn=el2-dnie
        riser=type4(head)
        go to 550
      end if

      head=el1-upie
      hdn=el2-dnie

c     if flow is into lower slough and downstream surface elevation is less
c     than crest of the upstream weir

      qweir=0.385*wc*sqrt(2.*g)*width*head**(3./2.)

      bao = pi*dia**2/4.
      ro = bao/(pi*dia)
      fl = (dist*(fman*qweir)**2)/((cm*bao)**2*ro**(4./3.))

c     downstream head is greater than pipe diameter

      if(hdn.ge.dia)then
        droph=(qweir/(ec*bao))**2*(1.+closs)/(2.*g)-fall+hdn+fl
        if(droph.lt.hdi)then
          riser=qweir
          go to 550
        else
          hd=head+hdi
          riser=type4(hd)
        end if
      end if

c     downstream head is less than pipe diameter

      if(hdn.lt.dia)then
        droph=(qweir/(ec*bao))**2*(1.+closs)/(2.*g)-fall+hdn+fl
        if(droph.lt.hdi)then
          riser=qweir
         else
          hd=head+hdi
          riser=type5(hd)
        end if
      end if

550   continue
      end

c    algorithm for calculating flow in box culvert

      function boxculv(el1,el2)
      implicit double precision(b-h,o-z)
      common/culvert/dia,fman,fall,dist,ec,hdn,closs,g,cm,hdif,upie,
     .               dnie,wc,width,hdi

      if(el1.gt.el2)then
        head=el1-upie
        hdn=el2-dnie
      else
        head=el2-dnie
        hdn=el1-upie
      end if

      bar=head*width
      per=2.*head+width

      boxculv=sqrt((head+abs(fall)-hdn)/((fman**2*per**(4./3.)*dist)/
     .        (cm**2*bar**(10./3.)) + 1./(bar**2*2.*g)))

      end

*******************************************************************
**           S U B R O U T I N E   E N V I R P
*******************************************************************

      subroutine envirp(tmstrt,jday,end_run,volg,dlx)
      include 'w2.inc'
      real phcnt,docnt,algaecnt,colcnt,velcnt,phcnt2
      real jday
      logical end_run
      COMMON /GLOBLC/ JB, JC, IU, ID, KT, DLT, KB(IMP)
      common/hydrc1/ u(kmp,imp),w(kmp,imp),az(kmp,imp),rho(kmp,imp)
      COMMON /GEOMHC/ EL(KMP),  H(KMP),  HKT1(IMP), HKT2(IMP)
      COMMON /GEOMBC/ B(KMP,IMP),  BKT(IMP),  BH(KMP,IMP),
     .      BHKT1(IMP), BHKT2(IMP), BHRKT1(IMP)
      parameter(numclass=20)
      dimension phclass(numclass),doclass(numclass),
     .          algaeclass(numclass),calgae(kmp,imp),
     .          colclass(numclass),velclass(numclass),dlx(imp),
     .          phclass2(numclass)
       data  timinv /0.02083/

c  initializing variables at first call
      if(jday.lt.tmstrt+1.1*timinv.and.ilast.ne.9)then
        dltt=timinv*86400.
        do i=1,numclass
          phclass(i)=0.0
          phclass2(i)=0.0
          doclass(i)=0.0
          algaeclass(i)=0.0
          velclass(i)=0.0
          colclass(i)=0.0
        end do
        phsum=0.0
        phsum2=0.0
        algaesum=0.0
        dosum=0.0
        velsum=0.0
        colsum=0.0
        ilast=9

c in this implementation these are added to during the simulation SW 1/95

      phcnt=0
      phcnt2=0
      docnt=0
      algaecnt=0
      colcnt=0
      velcnt=0
      phtot=0.
      algaetot=0.
      dotot=0.
      coltot=0.0
      veltot=0.0
      sumvolt=0.0
      open(192,file='envirp.npt',status='old')
      read(192,'(1x)')
      read(192,'(1x)')
      read(192,'(2f8.3)')phcrit,phincr
      read(192,'(1x)')
      read(192,'(2f8.3)')phcrit2,phincr2
      read(192,'(1x)')
      read(192,'(2f8.3)')docrit,doincr
      read(192,'(1x)')
      read(192,'(2f8.3)')algaecrit,algaeincr
      read(192,'(1x)')
      read(192,'(2f8.3)')velcrit,velincr
      read(192,'(1x)')
      read(192,'(2f8.3)')colcrit,colincr
      read(192,'(1x)')
      read(192,'(f8.3)')chlaconv
      close(192)

      end if

      if(end_run)go to 650

c start loop for succeeding calls to subroutine

      do i=iu,id
        do k=kt,kb(i)
              IF(K.EQ.KT)VOL=DLX(I)*BHKT2(I)
              IF(K.NE.KT)VOL=DLX(I)*B(K,I)*H(K)

c   checking and grouping ph violations - upper limit
          if(ph(k,i).ge.phcrit)then
            phtot=phtot+ph(k,i)*VOL*dltt
            phcnt=phcnt+VOL*dltt
            phlow=phcrit
c  determines which histogram class each violation falls within
            do jj=1,numclass
              if(ph(k,i).ge.phlow.and.ph(k,i).le.phlow+phincr)then
c  instead of summing violations in each histo class, i just add another
c  delta t to that class total (the result is the same)
                phclass(jj)=phclass(jj)+dltt*VOL
                go to 200
              else
                phlow=phlow+phincr
              end if
c  groups all violations with values exceeding second-to-last class
c   into last class
              if(jj.eq.numclass.and.ph(k,i).gt.phlow+phincr)then
                phclass(jj)=phclass(jj)+dltt*VOL
              end if
            end do
          end if
200       continue

c   checking and grouping ph violations -lower limit
          if(ph(k,i).lt.phcrit2)then
            phtot2=phtot2+ph(k,i)*VOL*dltt
            phcnt2=phcnt2+VOL*dltt
            phhigh2=phcrit2
            do jj=1,numclass
              if(ph(k,i).lt.phhigh2.and.ph(k,i).ge.phhigh2-phincr2)then
                phclass2(jj)=phclass2(jj)+dltt*VOL
                go to 202
              else
                phhigh2=phhigh2-phincr2
              end if
              if(jj.eq.numclass.and.ph(k,i).lt.phhigh2-phincr2)then
                phclass2(jj)=phclass2(jj)+dltt*VOL
              end if
            end do
          end if
202       continue

c  likewise for do
          if(do(k,i).lt.docrit)then
             dotot=dotot+do(k,i)*VOL*dltt
             docnt=docnt+VOL*dltt
             dohigh=docrit
             do jj=1,numclass
               if(do(k,i).lt.dohigh.and.do(k,i).ge.dohigh-doincr)then
                 doclass(jj)=doclass(jj)+dltt*VOL
                 go to 300
               else
                 dohigh=dohigh-doincr
               end if
               if(jj.eq.numclass.and.do(k,i).lt.dohigh-doincr)then
                 doclass(jj)=doclass(jj)+dltt*VOL
               end if
            end do
          end if
300       continue

c  likewise for algae
          calgae(k,i)=algae(k,i)*chlaconv
          if(calgae(k,i).ge.algaecrit)then
            algaetot=algaetot+calgae(k,i)*VOL*dltt
            algaecnt=algaecnt+VOL*dltt
            algaelow=algaecrit
            do jj=1,numclass
              if(calgae(k,i).ge.algaelow.and.calgae(k,i).le.algaelow+
     .           algaeincr)then
                algaeclass(jj)=algaeclass(jj)+dltt*VOL
                go to 400
              else
                algaelow=algaelow+algaeincr
              end if
              if(jj.eq.numclass.and.calgae(k,i).gt.algaelow
     .           +algaeincr)then
                algaeclass(jj)=algaeclass(jj)+dltt*VOL
              end if
            end do
          end if
400       continue

c   checking and grouping velocities
          if(abs(u(k,i)).ge.velcrit)then
            veltot=veltot+abs(u(k,i))*VOL*dltt
            velcnt=velcnt+VOL*dltt
            vellow=velcrit
c  determines which histogram class each violation falls within
            do jj=1,numclass
              if(abs(u(k,i)).ge.vellow.and.abs(u(k,i)).le.vellow+
     .         velincr)then
c  instead of summing violations in each histo class, i just add another
c  delta t to that class total (the result is the same)
                velclass(jj)=velclass(jj)+dltt*VOL
                go to 500
              else
                vellow=vellow+velincr
              end if
c  groups all violations with values exceeding second-to-last class
c   into last class
              if(jj.eq.numclass.and.abs(u(k,i)).gt.vellow+velincr)then
                velclass(jj)=velclass(jj)+dltt*VOL
              end if
            end do
          end if
500       continue

c   checking and grouping coliform concentrations
          if(colfrm(k,i).ge.colcrit)then
            coltot=coltot+colfrm(k,i)*VOL*dltt
            colcnt=colcnt+VOL*dltt
            collow=colcrit
c  determines which histogram class each violation falls within
            do jj=1,numclass
              if(colfrm(k,i).ge.collow.and.colfrm(k,i).le.
     .             collow*colincr)then
c  instead of summing violations in each histo class, i just add another
c  delta t to that class total (the result is the same)
                colclass(jj)=colclass(jj)+dltt*VOL
                go to 600
              else
                collow=collow*colincr
              end if
c  groups all violations with values exceeding second-to-last class
c   into last class
              if(jj.eq.numclass.and.colfrm(k,i).gt.collow*colincr)then
                colclass(jj)=colclass(jj)+dltt*VOL
              end if
            end do
          end if
600       continue

        end do
       end do

c sum of volg*dltt for volume fraction calculation

       sumvolt=sumvolt+VOLG*dltt


650   continue



      if(end_run)then
c  calculating average violation concentration and writing to file
        if(phcnt.gt.0.0)then
          avgph=phtot/phcnt
        else
          avgph=0.0
        end if
        if(phcnt2.gt.0.0)then
          avgph2=phtot2/phcnt2
        else
          avgph2=0.0
        end if
        if(docnt.gt.0.0)then
          avgdo=dotot/docnt
        else
          avgdo=0.0
        end if
        if(algaecnt.gt.0.0)then
          avgalgae=algaetot/algaecnt
        else
          avgalgae=0.0
        end if
        if(velcnt.gt.0.0)then
          avgvel=veltot/velcnt
        else
          avgvel=0.0
        end if
        if(colcnt.gt.0.0)then
          avgcol=coltot/colcnt
        else
          avgcol=0.0
        end if

        open(192,file='viol2.opt',status='new')
c  writing number representing violations for each class and also number
c  that is minimum value of that class
        phint=phcrit
        phint2=phcrit2-phincr2*real(numclass)
        doint=docrit-doincr*real(numclass)
        algaeint=algaecrit
        velint=velcrit
        colint=colcrit
        do i=1,numclass
        write(192,125)phint,phclass(i)/sumvolt,
     .  phint2,phclass2(numclass+1-i)/sumvolt, doint,
     .  doclass(numclass+1-i)/sumvolt, algaeint,algaeclass(i)/sumvolt,
     .  velint,velclass(i)/sumvolt, colint,colclass(i)/sumvolt
125     format(5(f6.2,3x,e12.4,3x),e12.4,3x,e12.4)
          phsum=phsum+phclass(i)/sumvolt
          phsum2=phsum2+phclass2(numclass+1-i)/sumvolt  
          dosum=dosum+doclass(numclass+1-i)/sumvolt
          colsum=colsum+colclass(i)/sumvolt
          velsum=velsum+velclass(i)/sumvolt
          algaesum=algaesum+algaeclass(i)/sumvolt
          phint=phint+phincr
          phint2=phint2+phincr2
          doint=doint+doincr
          algaeint=algaeint+algaeincr
          velint=velint+velincr
          colint=colint*colincr
        end do
        write(192,'(1x)')
        write(192,126)phsum,phsum2,dosum,algaesum,velsum,colsum
126     format(' 0 ',e12.4,' 0 ',e12.4, ' 0 ',e12.4,' 0 ',e12.4,' 0 ',
     .   e12.4,' 0 ',e12.4)
        write(192,'(1x)')
        write(192,126)avgph,avgph2,avgdo,avgalgae,avgvel,avgcol
        close(192)
      end if

      return

      end
*******************************************************************
**           S U B R O U T I N E   E N V I R P 2
*******************************************************************

      subroutine envirp2(tmstrt,jday,end_run,volg,dlx)
      include 'w2.inc'
      real phcnt,docnt,algaecnt,colcnt,velcnt,phcnt2
      real jday
      logical end_run
      COMMON /GLOBLC/ JB, JC, IU, ID, KT, DLT, KB(IMP)
      common/hydrc1/ u(kmp,imp),w(kmp,imp),az(kmp,imp),rho(kmp,imp)
      COMMON /GEOMHC/ EL(KMP),  H(KMP),  HKT1(IMP), HKT2(IMP)
      COMMON /GEOMBC/ B(KMP,IMP),  BKT(IMP),  BH(KMP,IMP),
     .      BHKT1(IMP), BHKT2(IMP), BHRKT1(IMP)
      parameter(numclass=20)
      dimension phclass(numclass),doclass(numclass),
     .          algaeclass(numclass),calgae(kmp,imp),
     .          colclass(numclass),velclass(numclass),dlx(imp),
     .          phclass2(numclass)
       data  timinv /0.02083/

c  initializing variables at first call
      if(jday.lt.tmstrt+1.1*timinv.and.ilast.ne.9)then
        dltt=timinv*86400.
        do i=1,numclass
          phclass(i)=0.0
          phclass2(i)=0.0
          doclass(i)=0.0
          algaeclass(i)=0.0
          velclass(i)=0.0
          colclass(i)=0.0
        end do
        phsum=0.0
        phsum2=0.0
        algaesum=0.0
        dosum=0.0
        velsum=0.0
        colsum=0.0
        ilast=9

c in this implementation these are added to during the simulation SW 1/95

      phcnt=0
      phcnt2=0
      docnt=0
      algaecnt=0
      colcnt=0
      velcnt=0
      phtot=0.
      algaetot=0.
      dotot=0.
      coltot=0.0
      veltot=0.0
      sumvolt=0.0
      open(192,file='envirp.npt',status='old')
      read(192,'(1x)')
      read(192,'(1x)')
      read(192,'(2f8.3)')phcrit,phincr
      read(192,'(1x)')
      read(192,'(2f8.3)')phcrit2,phincr2
      read(192,'(1x)')
      read(192,'(2f8.3)')docrit,doincr
      read(192,'(1x)')
      read(192,'(2f8.3)')algaecrit,algaeincr
      read(192,'(1x)')
      read(192,'(2f8.3)')velcrit,velincr
      read(192,'(1x)')
      read(192,'(2f8.3)')colcrit,colincr
      read(192,'(1x)')
      read(192,'(f8.3)')chlaconv
      close(192)

      end if

      if(end_run)go to 650

c start loop for succeeding calls to subroutine

      do i=iu,id
        do k=kt,kb(i)
              IF(K.EQ.KT)VOL=DLX(I)*BHKT2(I)
              IF(K.NE.KT)VOL=DLX(I)*B(K,I)*H(K)

c   checking and grouping ph violations - upper limit
          if(ph(k,i).ge.phcrit)then
            phtot=phtot+ph(k,i)*VOL*dltt
            phcnt=phcnt+VOL*dltt
            phlow=phcrit
c  determines which histogram class each violation falls within
            do jj=1,numclass
              if(ph(k,i).ge.phlow.and.ph(k,i).le.phlow+phincr)then
c  instead of summing violations in each histo class, i just add another
c  delta t to that class total (the result is the same)
                phclass(jj)=phclass(jj)+dltt*VOL
                go to 200
              else
                phlow=phlow+phincr
              end if
c  groups all violations with values exceeding second-to-last class
c   into last class
              if(jj.eq.numclass.and.ph(k,i).gt.phlow+phincr)then
                phclass(jj)=phclass(jj)+dltt*VOL
              end if
            end do
          end if
200       continue

c   checking and grouping ph violations -lower limit
          if(ph(k,i).lt.phcrit2)then
            phtot2=phtot2+ph(k,i)*VOL*dltt
            phcnt2=phcnt2+VOL*dltt
            phhigh2=phcrit2
            do jj=1,numclass
              if(ph(k,i).lt.phhigh2.and.ph(k,i).ge.phhigh2-phincr2)then
                phclass2(jj)=phclass2(jj)+dltt*VOL
                go to 202
              else
                phhigh2=phhigh2-phincr2
              end if
              if(jj.eq.numclass.and.ph(k,i).lt.phhigh2-phincr2)then
                phclass2(jj)=phclass2(jj)+dltt*VOL
              end if
            end do
          end if
202       continue

c  likewise for do
          if(do(k,i).lt.docrit)then
             dotot=dotot+do(k,i)*VOL*dltt
             docnt=docnt+VOL*dltt
             dohigh=docrit
             do jj=1,numclass
               if(do(k,i).lt.dohigh.and.do(k,i).ge.dohigh-doincr)then
                 doclass(jj)=doclass(jj)+dltt*VOL
                 go to 300
               else
                 dohigh=dohigh-doincr
               end if
               if(jj.eq.numclass.and.do(k,i).lt.dohigh-doincr)then
                 doclass(jj)=doclass(jj)+dltt*VOL
               end if
            end do
          end if
300       continue

c  likewise for algae
          calgae(k,i)=algae(k,i)*chlaconv
          if(calgae(k,i).ge.algaecrit)then
            algaetot=algaetot+calgae(k,i)*VOL*dltt
            algaecnt=algaecnt+VOL*dltt
            algaelow=algaecrit
            do jj=1,numclass
              if(calgae(k,i).ge.algaelow.and.calgae(k,i).le.algaelow+
     .           algaeincr)then
                algaeclass(jj)=algaeclass(jj)+dltt*VOL
                go to 400
              else
                algaelow=algaelow+algaeincr
              end if
              if(jj.eq.numclass.and.calgae(k,i).gt.algaelow
     .           +algaeincr)then
                algaeclass(jj)=algaeclass(jj)+dltt*VOL
              end if
            end do
          end if
400       continue

c   checking and grouping velocities
          if(abs(u(k,i)).ge.velcrit)then
            veltot=veltot+abs(u(k,i))*VOL*dltt
            velcnt=velcnt+VOL*dltt
            vellow=velcrit
c  determines which histogram class each violation falls within
            do jj=1,numclass
              if(abs(u(k,i)).ge.vellow.and.abs(u(k,i)).le.vellow+
     .         velincr)then
c  instead of summing violations in each histo class, i just add another
c  delta t to that class total (the result is the same)
                velclass(jj)=velclass(jj)+dltt*VOL
                go to 500
              else
                vellow=vellow+velincr
              end if
c  groups all violations with values exceeding second-to-last class
c   into last class
              if(jj.eq.numclass.and.abs(u(k,i)).gt.vellow+velincr)then
                velclass(jj)=velclass(jj)+dltt*VOL
              end if
            end do
          end if
500       continue

c   checking and grouping coliform concentrations
          if(colfrm(k,i).ge.colcrit)then
            coltot=coltot+colfrm(k,i)*VOL*dltt
            colcnt=colcnt+VOL*dltt
            collow=colcrit
c  determines which histogram class each violation falls within
            do jj=1,numclass
              if(colfrm(k,i).ge.collow.and.colfrm(k,i).le.
     .             collow*colincr)then
c  instead of summing violations in each histo class, i just add another
c  delta t to that class total (the result is the same)
                colclass(jj)=colclass(jj)+dltt*VOL
                go to 600
              else
                collow=collow*colincr
              end if
c  groups all violations with values exceeding second-to-last class
c   into last class
              if(jj.eq.numclass.and.colfrm(k,i).gt.collow*colincr)then
                colclass(jj)=colclass(jj)+dltt*VOL
              end if
            end do
          end if
600       continue

        end do
       end do

c sum of volg*dltt for volume fraction calculation

       sumvolt=sumvolt+VOLG*dltt


650   continue



      if(end_run)then
c  calculating average violation concentration and writing to file
        if(phcnt.gt.0.0)then
          avgph=phtot/phcnt
        else
          avgph=0.0
        end if
        if(phcnt2.gt.0.0)then
          avgph2=phtot2/phcnt2
        else
          avgph2=0.0
        end if
        if(docnt.gt.0.0)then
          avgdo=dotot/docnt
        else
          avgdo=0.0
        end if
        if(algaecnt.gt.0.0)then
          avgalgae=algaetot/algaecnt
        else
          avgalgae=0.0
        end if
        if(velcnt.gt.0.0)then
          avgvel=veltot/velcnt
        else
          avgvel=0.0
        end if
        if(colcnt.gt.0.0)then
          avgcol=coltot/colcnt
        else
          avgcol=0.0
        end if

        open(192,file='viol3.opt',status='new')
c  writing number representing violations for each class and also number
c  that is minimum value of that class
        phint=phcrit
        phint2=phcrit2-phincr2*real(numclass)
        doint=docrit-doincr*real(numclass)
        algaeint=algaecrit
        velint=velcrit
        colint=colcrit
        do i=1,numclass
        write(192,125)phint,phclass(i)/sumvolt,
     .  phint2,phclass2(numclass+1-i)/sumvolt, doint,
     .  doclass(numclass+1-i)/sumvolt, algaeint,algaeclass(i)/sumvolt,
     .  velint,velclass(i)/sumvolt, colint,colclass(i)/sumvolt
125     format(5(f6.2,3x,e12.4,3x),e12.4,3x,e12.4)
          phsum=phsum+phclass(i)/sumvolt
          phsum2=phsum2+phclass2(numclass+1-i)/sumvolt  
          dosum=dosum+doclass(numclass+1-i)/sumvolt
          colsum=colsum+colclass(i)/sumvolt
          velsum=velsum+velclass(i)/sumvolt
          algaesum=algaesum+algaeclass(i)/sumvolt
          phint=phint+phincr
          phint2=phint2+phincr2
          doint=doint+doincr
          algaeint=algaeint+algaeincr
          velint=velint+velincr
          colint=colint*colincr
        end do
        write(192,'(1x)')
        write(192,126)phsum,phsum2,dosum,algaesum,velsum,colsum
126     format(' 0 ',e12.4,' 0 ',e12.4, ' 0 ',e12.4,' 0 ',e12.4,' 0 ',
     .   e12.4,' 0 ',e12.4)
        write(192,'(1x)')
        write(192,126)avgph,avgph2,avgdo,avgalgae,avgvel,avgcol
        close(192)
      end if

      return

      end

  
c***************************************
        subroutine mass_load
        INCLUDE 'w2.inc'
        SAVE
****** Variable declaration

        INTEGER  OTQ,    TRQ,    WDQ,    TRT,    TRC,    DTQ,    DTT,
     .           DTC,    PRE,    PRT,    PRC,    UHE,    UHT,    UHC,
     .           DHE,    DHT,    DHC
        INTEGER  CN,     INCN,   TRCN,   DTCN,   PRCN,   UHCN,   DHCN
        INTEGER  UNIT,   WSCDP, cus, ds
        REAL NH3REL, NH3D, NO3D, NH3DK, NO3DK
        REAL JDAY
        LOGICAL  OPEN_FILES,     CONSTITUENTS, WITHDRAWALS, TRIBUTARIES,
     .           PRECIPITATION,  DIST_TRIBS,   UH_EXTERNAL, DH_EXTERNAL,
     .           UP_FLOW,        DN_FLOW,      POINT_SINK, END_RUN,
     .           UP_HEAD, DN_HEAD, UH_INTERNAL, DH_INTERNAL

******* Dimension declaration

        DIMENSION TRQ(NTP), TRT(NTP), TRC(NTP)
        DIMENSION INQ(NBP), INT(NBP), INC(NBP), DTQ(NBP), DTT(NBP),
     .            DTC(NBP), PRE(NBP), PRT(NBP), PRC(NBP), UHE(NBP),
     .            UHT(NBP), UHC(NBP), DHE(NBP), DHT(NBP), DHC(NBP),
     .            OTQ(NBP),wdq(nwp)                                     !4/5/95

******* Common declaration

        COMMON /GEOMBC/ B(KMP,IMP),     BKT(IMP),       BH(KMP,IMP),
     .                  BHKT1(IMP),     BHKT2(IMP),     BHRKT1(IMP)
        COMMON /SCRNC2/ NWD,    NTR
        COMMON /PRNTC3/ CUS(nbp),   DS(nbp)

        COMMON /TVDMTC/ TAIR,   TDEW,   CLOUD,  PHI,   ET,     CSHE,
     .                SRO,    LAT,    LONG,   banksh(imp)
        COMMON /GLOBLC/ JB, JC, IU, ID, KT, DLT, KB(IMP)
        COMMON /GLBLCC/ PALT,   ALGDET,  O2LIM,   WIND,  WSCDP, WSC(NDP)
        COMMON /TVDQC/  QIN(NBP),      QTR(NTP),      QDTR(NBP),
     .                PR(NBP),       ELUH(NBP),     ELDH(NBP),
     .                QOUT(KMP,NBP), QWD(NWP)
        COMMON /TVDTC/  TIN(NBP),         TTR(NTP),      TDTR(NBP),
     .                  TPR(NBP),         TUH(KMP,NBP),  TDH(KMP,NBP)
        COMMON /TVDCC1/ CIN(NCP,NBP),     CTR(NCP,NTP),
     .                  CDTR(NCP,NBP),    CPR(NCP,NBP),
     .                  CUH(KMC,NCP,NBP), CDH(KMC,NCP,NBP)
        COMMON /OXYGNC/ O2ORG,  O2ALG,  O2NH3, O2RESP, o2xfact
        COMMON /TVDCC2/ INCN(NCP),        TRCN(NCP),     DTCN(NCP),
     .                  PRCN(NCP),        UHCN(NCP),     DHCN(NCP)
        COMMON /TVDCC3/ NACIN, NACTR,     NACPR, NACDT
        COMMON /TVDLC1/ PRECIPITATION,    WITHDRAWALS,   TRIBUTARIES,
     .                  DIST_TRIBS(NBP)
        COMMON /GRTVDC/ CONSTITUENTS,     CN(NCP),       NAC
       COMMON /SELWC/  NSTR(NBP),     QSTR(NSP,NBP), ESTR(NSP,NBP),
     .                WSTR(NSP,NBP), KBSW(NBP),     KTOPSW(NBP),
     .                KBOTSW(NBP),   NOUT(NBP),     KOUT(KMP,NBP)
        COMMON /TVDLC2/ UP_FLOW(NBP),        DN_FLOW(NBP),
     .                  UH_INTERNAL(NBP),    UH_EXTERNAL(NBP),
     .                  DH_INTERNAL(NBP),    DH_EXTERNAL(NBP)

        COMMON /NITROC/ BION, PARTN, NH3DK, NH3REL, NO3DK, obion
        COMMON /DKMLTC/ A1(KMC,IMC),    A2(KMC,IMC),   A3(KMC,IMC)
        COMMON /DKSEDC/ SEDD(KMC,IMC),  SO2D(KMC,IMC), SOD(IMP)
        COMMON /PHOSPC/ PO4REL, BIOP, PARTP, obiop
        common /massbal/ c1s(kmc,imc,ncp),dlx(imp),itr(ntp),iwd(nwp)
     .                   ,up_head(nbp),dn_head(nbp),jbwd(nwp),kwd(nwp)
     .                   ,qpr(imp),qinf(kmp,nbp),quh1(kmp,nbp),
     .                    qdh1(kmp,nbp),end_run,jbtr(ntp),
     .                    qtrf(kmp,ntp),qdt(imp)

        real mbtrib,mbqwd,mbpre,mbufl,mbdfl,mbuh,mbuhi,mbdh,mbdhi,
     .       mbdtrib
        dimension mbtrib(ntp,ncp),mbqwd(nwp,ncp),mbpre(nbp,ncp),
     .    mbufl(nbp,ncp),mbdtrib(nbp,ncp),
     .    mbdfl(nbp,ncp),mbuh(nbp,ncp),mbuhi(nbp,ncp),mbdh(nbp,ncp),
     .    mbdhi(nbp,ncp)

c note this conversion converts the time step of 30 min into fraction of a
c day and then converts days to seconds and divides by 1000 to convert
c Q*C*dt in m3/s*mg/l*day to kg
         dtm=(0.5/24.)*86.4
c conversion for sod units since sod units are in g/m2/day and so2d is
c in sod*area/vol==g/m3/day  then *dtsod==kg/m3 then multiply by cell vol
c to get kg
         dtsod=(0.5/24.)*1000.
         if(.not.end_run)then

         do jb=1,nbp
         iu=cus(jb)
         id=ds(jb)

         DO 100 JAC=1,NAC
              JC = CN(JAC)
c define which parameters to do mass balance on
c BOD=5,PO4=9,NH3=10,NO3=11
         if(jc.ne.5.and.jc.ne.9.and.jc.ne.10.and.jc.ne.11)go to 100

************* External sources/sinks

              IF (TRIBUTARIES) THEN
                DO JT=1,NTR                                             !4/10/94
                  IF (JB.EQ.JBTR(JT)) THEN
                    I = ITR(JT)
                    IF (I.LT.CUS(JB)) I = CUS(JB)
                    DO K=KT,KB(I)
                     MBTRIB(JT,JC) = MBTRIB(JT,JC)+CTR(JC,JT)*QTR(JT)
     .                               *QTRF(K,JT)*dtm
                    END DO
                  END IF
                END DO
              END IF
              IF (DIST_TRIBS(JB)) THEN
                DO I=IU,ID
                MBDTRIB(JB,JC) = MBDTRIB(JB,JC)+CDTR(JC,JB)*QDT(I)*dtm
                END DO
              END IF
              IF (WITHDRAWALS) THEN
                DO JW=1,NWD                                             !4/10/94
                  IF (JB.EQ.JBWD(JW)) THEN
                    I = IWD(JW)
                    K = KWD(JW)
                    IF (K.LT.KT) K = KT                                 !12/9/94
c                  MBQWD(JW,JC) = MBQWD(JW,JC)-C1S(K,I,JC)*QWD(JW)*dtm  !12/9/94
                 vol=bhkt1(i)+bh(kt+1,i)+bh(kt+2,i)                      !4/5/95
                 qfrac1=bhkt1(i)/vol                                     !4/5/95
                 qfrac2=bh(kt+1,i)/vol                                   !4/5/95
                 qfrac3=bh(kt+2,i)/vol                                   !4/5/95
                MBQWD(JW,JC) = MBQWD(JW,JC)-C1S(kt,I,JC)*
     .                            QWD(JW)*QFRAC1*dtm                         !4/5/95
                MBQWD(JW,JC) = MBQWD(JW,JC)-C1S(kt+1,I,JC)*
     .                            QWD(JW)*QFRAC2*dtm                         !4/5/95
                MBQWD(JW,JC) = MBQWD(JW,JC)-C1S(kt+2,I,JC)*
     .                            QWD(JW)*QFRAC3*dtm                         !4/5/95
                  END IF
                END DO
              END IF
              IF (PRECIPITATION) THEN
                DO I=IU,ID
                  MBPRE(JB,JC) = MBPRE(JB,JC)+CPR(JC,JB)*QPR(I)*dtm
                END DO
              END IF
              IF (UP_FLOW(JB)) THEN
                DO K=KT,KB(IU)
                  MBUFL(JB,JC) = MBUFL(JB,JC)+QINF(K,JB)*QIN(JB)
     .                            *CIN(JC,JB)*dtm
                END DO
              END IF
              IF (DN_FLOW(JB)) THEN
c                DO K=KT,KB(ID)                                          !3/01/94
                 do jo=1,nout(jb)
                 k=kout(jo,jb)
                  MBDFL(JB,JC) = MBDFL(JB,JC)-QOUT(jo,JB)                !3/01/94
     .                            *C1S(K,ID,JC)*dtm                         !3/01/94
                END DO
              END IF
              IF (UP_HEAD(JB)) THEN
                IUT = IU
                IF (QUH1(KT,JB).GE.0.0) IUT = IU-1
                CSSU = C1S(KT,IUT,JC)*QUH1(KT,JB)
                MBUH(JB,JC)   = MBUH(JB,JC)+CSSU*dtm
                DO K=KT+1,KB(IU)
                  IUT = IU
                  IF (QUH1(K,JB).GE.0.0) IUT = IU-1
                  CSSU = C1S(K,IUT,JC)*QUH1(K,JB)
                  MBUH(JB,JC)   = MBUH(JB,JC)+CSSU*dtm
                END DO
c                IF (UH_INTERNAL(JB)) THEN
c                  I = UHS(JB)
c                  DO K=KT,KB(IU)
c                    MBUHI(JB,JC) = MBUHI(JB,JC)-CSSUH2(K,JC,JB)*dtm/DLT
c                  END DO
c                END IF
              END IF
              IF (DN_HEAD(JB)) THEN
                IDT = ID+1
                IF (QDH1(KT,JB).GE.0.0) IDT = ID
                CSSD = C1S(KT,IDT,JC)*QDH1(KT,JB)
                MBDH(JB,JC)   = MBDH(JB,JC)-CSSD*dtm
                DO K=KT+1,KB(ID+1)
                  IDT = ID+1
                  IF (QDH1(K,JB).GE.0.0) IDT = ID
                  CSSD = C1S(K,IDT,JC)*QDH1(K,JB)
                  MBDH(JB,JC)   = MBDH(JB,JC)-CSSD*dtm
                END DO
c                IF (DH_INTERNAL(JB)) THEN
c                  I = DHS(JB)
c                  DO K=KT,KB(ID+1)
c                    MBDHI(JB,JC) = MBDHI(JB,JC)+CSSDH2(K,JC,JB)*dtm/DLT
c                  END DO
c                END IF
              END IF
100         continue

c calculate internal loading of BOD, NH3, and PO4 through SOD only
c and Sediment decay compartments

c NH4,PO4,BOD
        DO I=IU,ID
          DO K=KT,KB(I)
          if(k.eq.kt)then
            vol=bhkt1(i)*dlx(i)
          else
            vol=bh(k,i)*dlx(i)
          end if
          sodnh3 = sodnh3+(NH3REL*SO2D(K,I)*A2(K,I))*vol*dtsod
          sodpo4 = sodpo4+(PO4REL*SO2D(K,I)*A2(K,I))*vol*dtsod
          sodbod = sodbod+(SO2D(K,I)*A3(K,I))*vol*dtsod
          sednh3 = sednh3+(BION*SEDD(K,I))*vol*dtsod
          sedpo4 = sedpo4+(BIOP*SEDD(K,I))*vol*dtsod
          sedbod = sedbod+(O2ORG*SEDD(K,I))*vol*dtsod
          END DO
        END DO

c end of do loop for each branch
        END DO
        ELSE
C END_RUN AND WRITE OUT RESULTS
        open(1,file='massbal.opt',status='new')
	write(1,300)
        if(tributaries)then
         do jt=1,ntr
         write(1,400)jt,itr(jt),mbtrib(jt,5),mbtrib(jt,9),
     .                          mbtrib(jt,10),mbtrib(jt,11)
         end do
        end if
        do jb=1,nbp
        if(dist_tribs(jb))then
         write(1,301)jb,mbdtrib(jb,5),mbdtrib(jb,9),
     .                  mbdtrib(jb,10),mbdtrib(jb,11)
        end if
        end do
   	  write(1,308)sodbod,sodpo4,sodnh3
          write(1,508)sedbod,sedpo4,sednh3

        do jb=1,nbp
        if(up_flow(jb))then
         write(1,304)jb,mbufl(jb,5),mbufl(jb,9),
     .                  mbufl(jb,10),mbufl(jb,11)
        end if
        end do
        do jb=1,nbp
        if(precipitation)then
         write(1,303)jb,mbpre(jb,5),mbpre(jb,9),
     .                  mbpre(jb,10),mbpre(jb,11)
        end if
        end do
        do jb=1,nbp
        if(dn_flow(jb))then
         write(1,305)jb,mbdfl(jb,5),mbdfl(jb,9),
     .                  mbdfl(jb,10),mbdfl(jb,11)
        end if
        end do
         if(withdrawals)then
         do jw=1,nwd
         write(1,302)jw,iwd(jw),mbqwd(jw,5),mbqwd(jw,9),
     .                  mbqwd(jw,10),mbqwd(jw,11)
         end do
        end if
        do jb=1,nbp
        if(up_head(jb))then
         write(1,306)jb,mbuh(jb,5),mbuh(jb,9),
     .                  mbuh(jb,10),mbuh(jb,11)
        end if
        end do
        do jb=1,nbp
        if(dn_head(jb))then
         write(1,307)jb,mbdh(jb,5),mbdh(jb,9),
     .                  mbdh(jb,10),mbdh(jb,11)
        end if
        end do

308	  format('"SOD Sources"',3(e12.3,2x))
508	  format('"SED Sources"',3(e12.3,2x))
       END IF
c format statements
300    format('"Type of Source/Sink"',1x,'"soluble BOD-kg"',
     .    1x,'"soluble PO4-P, kg"',1x,'"soluble NH4-N, kg"',1x,
     .    '"soluble NO3-N, kg"')
400    format('"Trib',i3,'/',i3,'"',4(e12.3,2x))
302    format('"QWD ',i3,'/',i3,'"',4(e12.3,2x))
301    format('"QDTRIB',i3,'"',4(e12.3,2x))
303    format('"Precip',i3,'"',4(e12.3,2x))
304    format('"Inflow',i3,'"',4(e12.3,2x))
305    format('"Outflow',i3,'"',4(e12.3,2x))
306    format('"Ext Up Head',i3,'"',4(e12.3,2x))
307    format('"Ext Down Head',i3,'"',4(e12.3,2x))
       return
       end





