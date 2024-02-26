SUBROUTINE ENDSIMULATION

USE MAIN
USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  USE MACROPHYTEC; USE POROSITYC; USE ZOOPLANKTONC; USE INITIALVELOCITY; USE BIOENERGETICS; USE TRIDIAG_V
  USE modSYSTDG, ONLY: DEALLOCATE_SYSTDG                      !  systdg
  ! CEMA testing start
  use CEMAVars
  ! CEMA testing end
  IMPLICIT NONE
  EXTERNAL RESTART_OUTPUT
  INTEGER IFILE 
!***********************************************************************************************************************************
!*                                                    Task 3: End Simulation                                                      **
!***********************************************************************************************************************************

 ! CEMA testing start
   !real jcs
   !!jcs=jcinzz*2.67/1000.0  ! converting C flux to DO flux, assuming 2.67 gO/gC
   !jcs=SD_jctest*2.67  ! converting C flux to DO flux, assuming 2.67 gO/gC
   ! !write(1081,'(5g12.5)')xjnh4,Jcinzz/1000.0,Jcs,MFTSedFlxVars(2,26)
   ! write(1081,'(5g12.5)')xjnh4,SD_Jctest,Jcs,MFTSedFlxVars(2,26)
 ! CEMA testing end

  CALL DATE_AND_TIME (CDATE,CCTIME)
  IF (.NOT. ERROR_OPEN) TEXT = 'Normal termination at '//CCTIME(1:2)//':'//CCTIME(3:4)//':'//CCTIME(5:6)//' on '//CDATE(5:6)//'/'     &
                                                       //CDATE(7:8)//'/'//CDATE(3:4)
  TEXT = ADJUSTL (TRIM(TEXT))
  CALL CPU_TIME (CURRENT)
  DO JW=1,NWB
    IF (SNAPSHOT(JW)) THEN
      WRITE (SNP(JW),'(/A/)')            ADJUSTL(TRIM(TEXT))
      WRITE (SNP(JW),'(A)')             'Runtime statistics'
      WRITE (SNP(JW),'(2(A,I0))')       '  Grid                 = ', IMX,' x ',KMX
      WRITE (SNP(JW),'(A,I0)')          '  Maximum active cells = ', NTACMX,'  Minimum active cells = ',NTACMN
      WRITE (SNP(JW),'(3(A,F0.1))')     '  Segment lengths, m   = ', DLXMIN,'-',DLXMAX
      WRITE (SNP(JW),'(3(A,F0.1))')     '  Layer heights, m     = ', HMIN,  '-',HMAX
      WRITE (SNP(JW),'(A)')             '  Timestep'
      WRITE (SNP(JW),'(A,I0)')          '    Total iterations   = ', NIT
      WRITE (SNP(JW),'(A,I0)')          '    # of violations    = ', NV
      WRITE (SNP(JW),'(A,F0.2)')        '    % violations       = ', FLOAT(NV)/FLOAT(NIT)*100.0
      WRITE (SNP(JW),'(A,I0,A)')        '    Average timestep   = ', INT(DLTAV),' sec'
      WRITE (SNP(JW),'(A,I0,A,F0.2,A)') '  Simulation time      = ', INT(ELTMJD),' days ',(ELTMJD-INT(ELTMJD))*24.0,' hours'
      WRITE (SNP(JW),'(A,F0.2,A)')      '  Total CPU runtime    = ',(CURRENT-START)/60.0,' min'
      CLOSE (SNP(JW))
    END IF
    !IF (VECTOR(JW))      CLOSE (VPL(JW))
    IF (PROFILE(JW))     CLOSE (PRF(JW))
    IF (SPREADSHEET(JW)) CLOSE (SPR(JW))
    IF (CONTOUR(JW))     CLOSE (CPL(JW))
  END DO
  
  ! *** DSI W2_TOOL LINKAGE
  IF(VECTOR(1))CLOSE(VPL(1))
  
  IF (TIME_SERIES) THEN
    DO J=1,NIKTSR
      CLOSE (TSR(J))
    END DO
    close(WLFN)   ! WL output file  ! SW 9/25/13
  END IF
  IF (WARNING_OPEN) THEN
    CLOSE (WRN)
  ELSE
    CLOSE (WRN,STATUS='DELETE')
  END IF
  IF (ERROR_OPEN) THEN
    CLOSE (W2ERR)
  ELSE
    CLOSE (W2ERR,STATUS='DELETE')
  END IF
  DO J=40,NOPEN
   CLOSE (J)
  END DO

      IF (FLOWBALC=='      ON') THEN        
      CLOSE(FLOWBFN)   ! flowbal file
      ENDIF

      IF (NPBALC=='      ON') THEN    
      CLOSE(MASSBFN)   ! MASS BALANCE file
      ENDIF
    
  IF(SELECTC == '      ON')then          ! SW 9/25/13 New Section on closing files
  ifile=1949
  do jb=1,nbr
      if(nstr(jb) > 0)then
          ifile=ifile+1
          close(ifile)
      endif
  enddo
  if(nwd > 0)then
      ifile=ifile+1
      close(ifile)
  endif
  do jw=1,nwb   ! sw 4/20/15
      ifile=ifile+1
      close(ifile)
  enddo

  endif
  
      IF (DOWNSTREAM_OUTFLOW) THEN  
      JFILE=0  
      DO JWD=1,NIWDO  
        CLOSE(WDO(JWD,1))
        CLOSE(WDO(JWD,2))
        IF (CONSTITUENTS) THEN  
          CLOSE (WDO(JWD,3))
        END IF  
        IF (DERIVED_CALC) THEN  
          CLOSE(WDO(JWD,4))
        END IF  
          
        ! Determine the # of withdrawals at the WITH SEG  
        DO JB=1,NBR  ! structures
          IF(IWDO(JWD)==DS(JB) .AND. NSTR(JB) /= 0)THEN  
            DO JS=1,NSTR(JB)  
              JFILE=JFILE+1  
              CLOSE(WDO2(JFILE,1))  
              CLOSE(WDO2(JFILE,2))
              IF (CONSTITUENTS) THEN  
                CLOSE(WDO2(JFILE,3))
              ENDIF   
              IF (DERIVED_CALC) THEN  
                CLOSE(WDO2(JFILE,4))
              ENDIF  
            ENDDO  
          ENDIF  
        ENDDO  
          
        DO JS=1,NWD  ! withdrawals
          IF(IWDO(JWD) == IWD(JS))THEN  
            JFILE=JFILE+1  
            CLOSE(WDO2(JFILE,1))
            CLOSE(WDO2(JFILE,2))
            IF (CONSTITUENTS) THEN  
                CLOSE(WDO2(JFILE,3)) 
            ENDIF  
            IF (DERIVED_CALC) THEN  
             CLOSE(WDO2(JFILE,4))
            ENDIF  
          ENDIF  
        ENDDO  
          
        DO JS=1,NSP  ! spillways
          IF(IWDO(JWD) == IUSP(JS))THEN  
            JFILE=JFILE+1  
            CLOSE(WDO2(JFILE,1))
            CLOSE(WDO2(JFILE,2))
            IF (CONSTITUENTS) THEN  
                CLOSE(WDO2(JFILE,3))
            ENDIF  
            IF (DERIVED_CALC) THEN  
               CLOSE(WDO2(JFILE,4))
            ENDIF  
          ENDIF  
        ENDDO  
          
        DO JS=1,NPU  ! pumps
          IF(IWDO(JWD) == IUPU(JS))THEN  
            JFILE=JFILE+1  
            CLOSE(WDO2(JFILE,1))
            CLOSE(WDO2(JFILE,2))
            IF (CONSTITUENTS) THEN  
                CLOSE(WDO2(JFILE,3))
            ENDIF  
            IF (DERIVED_CALC) THEN  
                CLOSE(WDO2(JFILE,4))
            ENDIF  
          ENDIF  
        ENDDO  
           
        DO JS=1,NPI  ! pipes
          IF(IWDO(JWD) == IUPI(JS))THEN  
            JFILE=JFILE+1  
            CLOSE(WDO2(JFILE,1))
            CLOSE(WDO2(JFILE,2))
            IF (CONSTITUENTS) THEN  
                CLOSE(WDO2(JFILE,3))
            ENDIF  
            IF (DERIVED_CALC) THEN  
                CLOSE(WDO2(JFILE,4))
            ENDIF  
          ENDIF  
        ENDDO  
          
        DO JS=1,NGT  ! gates
          IF(IWDO(JWD) == IUGT(JS))THEN  
            JFILE=JFILE+1  
            CLOSE(WDO2(JFILE,1))
            CLOSE(WDO2(JFILE,2))
            IF (CONSTITUENTS) THEN  
                CLOSE(WDO2(JFILE,3))
            ENDIF  
            IF (DERIVED_CALC) THEN  
                CLOSE(WDO2(JFILE,4))
            ENDIF  
          ENDIF  
        ENDDO  
          
      END DO  
    END IF  
  
 
  
  !OPEN(W2ERR,FILE='W2Errordump.opt',status='unknown')
  !WRITE(w2err,*)'JDAY',jday,'SZ',sz,'Z',z,'H2KT',h2(kt,1:imx),'H1KT',h1(kt,1:imx),'BHR1',bhr1(kt,1:imx),'BHR2',bhr2(kt,1:imx),'WSE',elws,'Q',q,'QC',qc,'QERR',qerr,'T1',t1(kt,1:imx),'T2',t2(kt,1:imx),'SUKT',su(kt,1:imx),&
  !                        'UKT',u(kt,1:imx),'QIN',qin,'QTR',qtr,'QWD',qwd
IF (ERROR_OPEN) THEN                                           ! modified to be more organized and comma-delimited  !SR 12/26/2019
    OPEN (W2ERR,FILE='W2Errordump.csv',status='unknown')         ! changed to csv                                     !SR 12/26/2019
    WRITE (W2ERR,*) 'JDAY = ', JDAY                                                                                   !SR 12/26/2019
    WRITE (W2ERR,'(A,1000(",",F0.6))') 'QIN:', (QIN(J),J=1,NBR)                                                       !SR 12/26/2019
    WRITE (W2ERR,'(A,1000(",",F0.6))') 'QTR:', (QTR(J),J=1,NTRT)                                                      !SR 12/26/2019
    WRITE (W2ERR,'(A,1000(",",F0.6))') 'QDT:', (QDTR(J),J=1,NBR)                                                      !SR 12/26/2019
    WRITE (W2ERR,'(A,1000(",",F0.6))') 'QWD:', (QWD(J),J=1,NWDT)                                                      !SR 12/26/2019
    WRITE (W2ERR,'(/A)') 'SEG,BRANCH,KT,WSE,SZ,Z,Q,QC,QERR,H2KT,H1KT,BHR1,BHR2,T1,T2,SUKT,UKT'                        !SR 12/26/2019
    DO JW=1,NWB                                                                                                       !SR 12/26/2019
      KT = KTWB(JW)                                                                                                   !SR 12/26/2019
      DO JB=BS(JW),BE(JW)                                                                                             !SR 12/26/2019
        DO I=US(JB)-1,DS(JB)+1                                                                                        !SR 12/26/2019
          WRITE (W2ERR,'(I0,",",I0,",",I0,14(",",F0.6))') I, JB, KT, ELWS(I), SZ(I), Z(I), Q(I), QC(I), QERR(I),                   &
                 H2(KT,I), H1(KT,I), BHR1(KT,I), BHR2(KT,I), T1(KT,I), T2(KT,I), SU(KT,I), U(KT,I)                    !SR 12/26/2019
        END DO                                                                                                        !SR 12/26/2019
      END DO                                                                                                          !SR 12/26/2019
    END DO                                                                                                            !SR 12/26/2019
    CLOSE(W2ERR)
  END IF

  DEALLOCATE (LAYERCHANGE,TECPLOT,X1, HAB)
  DEALLOCATE (TSR,    WDO,    WDO2, ETSR,   IWDO,   ITSR,   TITLE,  CDAC,   WSC,    ESTR,   WSTR,   QSTR,   KTSW,   KBSW,   SINKC, JBTSR, CONSTRICTION, BCONSTRICTION)   ! SW 7/24/2018  8/5/2018
  DEALLOCATE (EBC,    MBC,    PQC,    EVC,    PRC,    WINDC,  QINC,   QOUTC,  HEATC,  SLHTC,  QINIC,  DTRIC,  TRIC,   WDIC)
  DEALLOCATE (EXC,    EXIC,   VBC,    METIC,  SLTRC,  THETA,  FRICC,  NAF,    ELTMF,  ZMIN,   IZMIN,  C2CH,   CDCH,   EPCH, KFCH, APCH, ANCH, ALCH)
  DEALLOCATE (CPLTC,  HPLTC,  CMIN,   CMAX,   HYMIN,  HYMAX,  CDMIN,  CDMAX,  JBDAM,  ILAT,   CDPLTC, QINSUM, TINSUM, TIND)
  DEALLOCATE (QOLD,   DTP,    DTPS,   QOLDS,  QIND,   JSS,    HDIC,   QNEW,   YSS,    VSS,    YS,     VS,     VSTS,   NSPRF)
  DEALLOCATE (LATGTC, LATSPC, LATPIC, DYNPIPE,LATPUC, DYNGTC, OPT,    CIND,   CINSUM, CDWBC,  KFWBC,  CPRWBC, CINBRC, CTRTRC, CDTBRC, DYNPUMP)   ! SW 5/10/10
  DEALLOCATE (YSTS,   YST,    VST,    ALLIM,  APLIM,  ANLIM,  ASLIM,  ELLIM,  EPLIM,  ENLIM,  ESLIM,  CSSK,   C1,     C2, Z0)
  DEALLOCATE (KFS,    AF,     EF,     HYD,    KF,     AZSLC,  STRIC,  CPRBRC, CD,     KBMAX,  ELKT,   WIND2,  VISC,   CELC, DLTADD)
  DEALLOCATE (QOAVR,  QIMXR,  QOMXR,  REAERC, LAT,    LONGIT, ELBOT,  BTH,    VPR,    LPR,    NISNP,  NIPRF,  NISPR,  DECL)
  DEALLOCATE (A00,    HH,     T2I,    KTWB,   KBR,    IBPR,   DLVR,   ESR,    ETR,    NBL,    LPRFN,  EXTFN,  BTHFN,  METFN)
  DEALLOCATE (SNPFN,  PRFFN,  SPRFN,  CPLFN,  VPLFN,  FLXFN,  FLXFN2, VPRFN,  AFW,    BFW,    CFW,    WINDH,  RHEVC,  FETCHC, JBDN, SPRVFN)   ! SW 9/28/2018
  DEALLOCATE (KBI, MACCH, GRIDC, GMA, BTA, QTOT, SEDCIP, SEDCIN, SEDCIC,SEDCIS)            ! SW 9/27/2007
  DEALLOCATE (SDK,    FSOD,   FSED,   SEDCI,  SEDCC,   SEDPRC, ICEC,   SLICEC, ICETHI, ALBEDO, HWI,    BETAI,  GAMMAI,ICEMIN)
  DEALLOCATE (SEDS,   SEDB)    !CB 11/28/06
  DEALLOCATE (EXH2O,  BETA,   EXOM,   EXSS,   DXI,    CBHE,   TSED,   TSEDF,  FI,     ICET2,  AZC,    AZMAX)     ! QINT,   QOUTT
  DEALLOCATE (AX,     WTYPEC, TAIR,   TDEW,   WIND,   PHI,    CLOUD,  CSHE,   SRON,   RANLW,    RB,     RC,     RE,     SHADE)
  DEALLOCATE (ET,     RS,     RN,     SNPC,   SCRC,   PRFC,   SPRC,   CPLC,   VPLC,   FLXC,   NXTMCP, NXTMVP, NXTMFL, GAMMA, F_NH3)
  DEALLOCATE (NXTMSN, NXTMSC, NXTMPR, NXTMSP, SNPDP,  SCRDP,  PRFDP,  SPRDP,  CPLDP,  VPLDP,  FLXDP,  NCPL,   NVPL,   NFLX)
  DEALLOCATE (NSNP,   NSCR,   NPRF,   NSPR,   NEQN,   PO4R,   PARTP,  NH4DK,  NH4R,   NO3DK,  NO3S,   CDSUM)
  DEALLOCATE (SSCS, CH4R, H2SR, FEIIR, MNIIR,SO4R)
  DEALLOCATE (BACTQ10, BACT1DK, BACTLDK, BACTS)
  DEALLOCATE (CoeffA_Turb, CoeffB_Turb,SECC_PAR)
  DEALLOCATE (H2SQ10, H2S1DK, CH4Q10, CH41DK)
  DEALLOCATE (KFE_OXID, KFE_RED, KFEOOH_HalfSat, FeSetVel)
  DEALLOCATE (KMN_OXID, KMN_RED, KMNO2_HalfSat, MnSetVel,DGPO2,MINKL)
  DEALLOCATE (LPOMHK, RPOMHK)
  DEALLOCATE (LDOMPDK, LRDOMPDK, RDOMPDK, LDOMNDK, LRDOMNDK, RDOMNDK, LDOMCDK, LRDOMCDK, RDOMCDK)
  DEALLOCATE (LPOMPDK, LRPOMPDK, RPOMPDK, LPOMNDK, LRPOMNDK, RPOMNDK, LPOMCDK, LRPOMCDK, RPOMCDK)
  DEALLOCATE (LDOMCMP, LPOMCMP, RPOMCMP)
  DEALLOCATE (LPZOOINC,LPZOOOUTC)
  DEALLOCATE (LDOP, RDOP, LPOP, RPOP, LDON, RDON, LPON, RPON)
  DEALLOCATE (LDOC, RDOC, LPOC, RPOC)
  DEALLOCATE (PSIEM, SEDEB)                                                                  
  DEALLOCATE (LPOMPEP, LPOMNEP, LPOMCEP)
  DEALLOCATE (LPOMHD,  RPOMHD)
  DEALLOCATE (LDOMCAP, LDOMCEP, LPOMCAP, LPOMCNS, RPOMCNS)
  DEALLOCATE (LDOMPD,  LRDOMPD, RDOMPD,  LPOMPD,  LRPOMPD, RPOMPD,  LPOMPHD, RPOMPHD) 
  DEALLOCATE (LDOMND,  LRDOMND, RDOMND,  LPOMND,  LRPOMND, RPOMND,  LPOMNHD, RPOMNHD) 
  DEALLOCATE (LDOMCD,  LRDOMCD, RDOMCD,  LPOMCD,  LRPOMCD, RPOMCD,  LPOMCHD, RPOMCHD)

  DEALLOCATE (CO2R,   SROC,   O2ER,   O2EG,   CAQ10,  CADK,   CAS,    BODP,   BODN,   BODC,   KBOD,   TBOD,   RBOD,   DTRC)
  DEALLOCATE (LDOMDK, RDOMDK, LRDDK,  OMT1,   OMT2,   OMK1,   OMK2,   LPOMDK, RPOMDK, LRPDK,  POMS,   ORGP,   ORGN,   ORGC)
  DEALLOCATE (RCOEF1, RCOEF2, RCOEF3, RCOEF4, ORGSI,  NH4T1,  NH4T2,  NH4K1,  NH4K2,  NO3T1,  NO3T2,  NO3K1,  NO3K2,  NSTR)
  DEALLOCATE (DSIR,   PSIS,   PSIDK,  PARTSI, SODT1,  SODT2,  SODK1,  SODK2,  O2NH4,  O2OM,   O2AR,   O2AG,   CG1DK,  CGS)
  DEALLOCATE (CGQ10,  CG0DK,  CGLDK, CGKLF,CGCS,CGR,CUNIT,  CUNIT1, CUNIT2, CAC,    INCAC,  TRCAC,  DTCAC,  PRCAC,  CNAME,  CNAME1, CNAME2, CMULT) !LCJ 2/26/15
  DEALLOCATE (CN,     INCN,   DTCN,   PRCN,   CSUM,   DLTMAX, QWDO,   TWDO,   SSS,    SEDRC,  TAUCR,  XBR, FNO3SED, DYNSTRUC)
!  DEALLOCATE (SSFLOC, FLOCEQN)                                                 
  !DEALLOCATE (SEDCC1,SEDCC2, ICEQSS,SDK1,sdk2,SEDCI1,SEDCI2,SEDPRC1,SEDPRC2,SEDVP1,SEDVP2,SED1,SED2) 
  !DEALLOCATE (SEDCC1,SEDCC2, ICEQSS,SDK1,sdk2,SEDCI1,SEDCI2,SEDPRC1,SEDPRC2,SEDVP1,SEDVP2,SED1,SED2,fsedc1,fsedc2,pbiom,nbiom,cbiom)   ! Amaila, cb 6/7/17
  DEALLOCATE (SEDCC1,SEDCC2, ICEQSS,SDK1,sdk2,SEDCI1,SEDCI2,SEDPRC1,SEDPRC2,SEDVP1,SEDVP2,SED1,SED2,fsedc1,fsedc2,pbiom,nbiom,cbiom,sed1ic,sed2ic,sdfirstadd)   ! cb 9/3/17
  DEALLOCATE (ISO_SEDIMENT1, VERT_SEDIMENT1,LONG_SEDIMENT1,ISO_SEDIMENT2,VERT_SEDIMENT2,LONG_SEDIMENT2,PRINT_SEDIMENT1,PRINT_SEDIMENT2) 
  DEALLOCATE (SEDIMENT_CALC1,SEDIMENT_CALC2)
  DEALLOCATE (QTAVB,  QTMXB,  BS,     BE,     JBUH,   JBDH,   TSSS,   TSSB,   TSSICE, ESBR,   ETBR,   EBRI,   QDTR,   EVBR)
  DEALLOCATE (QIN,    PR,     QPRBR,  TIN,    TOUT,   TPR,    TDTR,   TPB,    NACPR,  NACIN,  NACDT,  NACTR,  NACD,   ELDH)
  DEALLOCATE (QSUM,   NOUT,   KTQIN,  KBQIN,  ELUH,   NL,     NPOINT, SLOPE,  SLOPEC, ALPHA,  COSA,   SINA,   SINAC,  TDHFN,  QOTFN,  PREFN)
  DEALLOCATE (CPRFN,  EUHFN,  TUHFN,  CUHFN,  EDHFN,  QINFN,  TINFN,  CINFN,  CDHFN,  QDTFN,  TDTFN,  CDTFN,  TPRFN,  VOLEV)
  DEALLOCATE (VOLWD,  VOLSBR, VOLTBR, DLVOL,  VOLG,   VOLSR,  VOLTR,  VOLB,   VOLPR,  VOLTRB, VOLDT,  VOLUH,  VOLDH,  VOLIN, VOLICE, ICEBANK)
  DEALLOCATE (US,     DS,     CUS,    UHS,    DHS,    UQB,    DQB,    CDHS,   VOLOUT, TSSWD,  TSSUH,  TSSDH,  TSSIN,  TSSOUT)
  DEALLOCATE (TSSEV,  TSSPR,  TSSTR,  TSSDT,  SOD,    ELWS,   BKT,    REAER,  ICETH,  ICE,    ICESW,  Q,      QC,     QERR)
  DEALLOCATE (KTI,    SROSH,  SEG,    DLXRHO, QSSUM,  DLX,    DLXR,   QUH1,   QDH1,   BI,     JWUH,   JWDH)
  DEALLOCATE (A,      C,      D,      F,      V,      SKTI,   KBMIN,  EV,     QDT,    QPR,    SBKT,   BHRHO)
  DEALLOCATE (SZ,     WSHX,   WSHY,   WIND10, CZ,     FETCH,  PHI0,   FRIC,   ADZ,    HMULT,  FMTC,   FMTCD,  CNAME3, CDNAME3)
  DEALLOCATE (Z,      KB,     VNORM,  ANPR,   ANEQN,  APOM,   ACHLA,  AHSP,   AHSN,   AHSSI)
  DEALLOCATE (AC,     ASI,    AT1,    AT2,    AT3,    AT4,    AK1,    AK2,    AK3,    AK4,    EXA,    ASAT,   AP,     AN, AVERTM)
  DEALLOCATE (AG,     AR,     AE,     AM,     AS,     ENPR,   ENEQN,  EG,     ER,     EE,     EM,     EB,     ESAT,   EP)
  DEALLOCATE (EC,     ESI,    ECHLA,  EHSP,   EHSN,   EHSSI,  EPOM,   EHS,    EN,     ET4,    EK1,    EK2,    EK3,    EK4)
  DEALLOCATE (ET1,    ET2,    ET3,    HNAME,  FMTH,    KFAC,  KFNAME, KFNAME2,KFCN,   C2I,    TRCN,   CDN,    CDNAME, CDNAME2,CDMULT)
  DEALLOCATE (CMBRS,  CMBRT,  FETCHU, FETCHD, IPRF,   ISNP,   ISPR,   BL,     LFPR,   DO3,    SED,    TKE,    PALT)
  DEALLOCATE (ADX,    DO1,    DO2,    B,      CONV,   CONV1,  EL,     DZ,     DZQ,    DX,     SAZ,    T1,TSS,QSS,BNEW, ILAYER)   ! SW 1/23/06
  DEALLOCATE (P,      SU,     SW,     BB,     BR,     BH,     BHR,    VOL,    HSEG,   DECAY,  FPFE,   FRICBR, UXBR,   UYBR)
  DEALLOCATE (DEPTHB, DEPTHM, FPSS,   TUH,    TDH,    TSSUH1, TSSUH2, TSSDH1, TSSDH2, SEDVP,  H,      EPC)
  DEALLOCATE (TVP,    QINF,   QOUT,   KOUT,   VOLUH2, VOLDH2, CWDO,   CDWDO,  CWDOC,  CDWDOC, CDTOT,  CPR,    CPB,    COUT)
  DEALLOCATE (CIN,    CDTR,   RSOD,   RSOF,   DLTD,   DLTF,   TSRD,   TSRF,   WDOD,   WDOF,   SNPD,   SNPF,   SPRD,   SPRF)
  DEALLOCATE (SCRD,   SCRF,   PRFD,   PRFF,   CPLD,   CPLF,   VPLD,   VPLF,   FLXD,   FLXF,   EPIC,   EPICI,  EPIPRC, EPIVP)
  DEALLOCATE (CUH,    CDH,    EPM,    EPD,    C1S,    CSSB,   CVP,    CSSUH1, CSSUH2, CSSDH2, CSSDH1, LNAME,  IWR,    KTWR, EKTWR, EKBWR)
  DEALLOCATE (JWUSP,  JWDSP,  QSP,    KBWR,   KTWD,   KBWD,   JBWD,   GTA1,   GTB1,   GTA2,   GTB2,   BGT,    IUGT,   IDGT)
  DEALLOCATE (QTR,    TTR,    KTTR,   KBTR,   EGT,    EGT2,   AGASGT, BGASGT, CGASGT, GASGTC, PUGTC,  ETUGT,  EBUGT,  KTUGT,  KBUGT)
  DEALLOCATE (PDGTC,  ETDGT,  EBDGT,  KTDGT,  KBDGT,  A1GT,   B1GT,   G1GT,   A2GT,   B2GT,   G2GT,   JWUGT,  JWDGT,  QGT)
  DEALLOCATE (EQGT,   JBUGT,  JBDGT,  JBUPI,  JBDPI,  JWUPI,  JWDPI,  QPI,    IUPI,   IDPI,   EUPI,   EDPI,   WPI,    DLXPI, BP)   ! SW 5/5/10
  DEALLOCATE (ETUPI,  EBUPI,  KTUPI,  KBUPI,  PDPIC,  ETDPI,  EBDPI,  KTDPI,  KBDPI,  FPI,    FMINPI, PUPIC,  ETDSP,  EBDSP)
  DEALLOCATE (PUSPC,  ETUSP,  EBUSP,  KTUSP,  KBUSP,  PDSPC,  KTDSP,  KBDSP,  IUSP,   IDSP,   ESP,    A1SP,   B1SP,   A2SP)
  DEALLOCATE (B2SP,   AGASSP, BGASSP, CGASSP, EQSP,   GASSPC, JBUSP,  JBDSP,  STRTPU, ENDPU,  EONPU,  EOFFPU, QPU,    PPUC)
  DEALLOCATE (IUPU,   IDPU,   EPU,    ETPU,   EBPU,   KTPU,   KBPU,   JWUPU,  JWDPU,  JBUPU,  JBDPU,  PUMPON, KTW,    KBW, PUMP_DOWNSTREAM)
  DEALLOCATE (IWD,    KWD,    QWD,    EWD,    ITR,    QTRFN,  TTRFN,  CTRFN,  ELTRT,  ELTRB,  TRC,    JBTR,   QTRF,   CLRB)
  DEALLOCATE (TTLB,   TTRB,   CLLB,   SRLB1,  SRRB1,  SRLB2,  SRRB2,  SRFJD1, SHADEI, SRFJD2, TOPO,   QSW,    CTR)    ! SW 10/17/05
  DEALLOCATE (H1,     H2,     BH1,    BH2,    BHR1,   BHR2,   AVH1,   AVH2,   SAVH2,  AVHR,   SAVHR,  CBODD)
  DEALLOCATE (POINT_SINK,         HPRWBC,   READ_EXTINCTION, READ_RADIATION)
  DEALLOCATE (DIST_TRIBS,     UPWIND,               ULTIMATE,           FRESH_WATER,      SALT_WATER,      LIMITING_FACTOR)
  DEALLOCATE (UH_EXTERNAL,    DH_EXTERNAL,          UH_INTERNAL,        DH_INTERNAL,      UQ_INTERNAL,     DQ_INTERNAL)
  DEALLOCATE (UQ_EXTERNAL,    DQ_EXTERNAL,          UP_FLOW,            DN_FLOW,          UP_HEAD,         DN_HEAD)
  DEALLOCATE (INTERNAL_FLOW,  DAM_INFLOW,           DAM_OUTFLOW,        HEAD_FLOW,        HEAD_BOUNDARY)      !TC 08/03/04
  DEALLOCATE (ISO_CONC,             VERT_CONC,          LONG_CONC,        VERT_SEDIMENT,   LONG_SEDIMENT)
  DEALLOCATE (ISO_SEDIMENT,   VISCOSITY_LIMIT,      CELERITY_LIMIT,     IMPLICIT_AZ,      ONE_LAYER,       IMPLICIT_VISC)
  DEALLOCATE (FETCH_CALC,     LIMITING_DLT,         TERM_BY_TERM,       MANNINGS_N,       PLACE_QTR,       SPECIFY_QTR)
  DEALLOCATE (PLACE_QIN,      PRINT_CONST,          PRINT_HYDRO,        PRINT_SEDIMENT,   ENERGY_BALANCE,  MASS_BALANCE)
  DEALLOCATE (VOLUME_BALANCE, DETAILED_ICE,         ICE_CALC,                ALLOW_ICE,       PH_CALC, BR_INACTIVE)     ! ICE_IN,       RC/SW 4/28/11
  DEALLOCATE (BOD_CALCP, BOD_CALCN) 
  DEALLOCATE (EVAPORATION,    PRECIPITATION,        RH_EVAP,            NO_INFLOW,        NO_OUTFLOW,      NO_HEAT, BR_NOTECPLOT)   ! SW 8/27/2019
  DEALLOCATE (ISO_TEMP,       VERT_TEMP,            LONG_TEMP,          VERT_PROFILE,     LONG_PROFILE,    NO_WIND)
  DEALLOCATE (SNAPSHOT,       PROFILE,              VECTOR,             CONTOUR,          SPREADSHEET,     INTERNAL_WEIR)
  DEALLOCATE (SCREEN_OUTPUT,  FLUX,                 DYNAMIC_SHADE,      TRAPEZOIDAL, BOD_CALC, ALG_CALC)
  DEALLOCATE (SEDIMENT_CALC,  EPIPHYTON_CALC,       PRINT_DERIVED,      PRINT_EPIPHYTON,  TDG_SPILLWAY,    TDG_GATE, DYNSEDK)
  DEALLOCATE (ISO_EPIPHYTON,  VERT_EPIPHYTON,       LONG_EPIPHYTON,     LATERAL_SPILLWAY, LATERAL_GATE,    LATERAL_PUMP)
  DEALLOCATE (iso_macrophyte,  vert_macrophyte,       long_macrophyte, macrcvp,   macrclp)  ! cb 8/21/15
  DEALLOCATE (INTERP_HEAD,    INTERP_WITHDRAWAL,    INTERP_EXTINCTION,  INTERP_DTRIBS,    LATERAL_PIPE,    INTERP_TRIBS)
  DEALLOCATE (INTERP_OUTFLOW, INTERP_INFLOW,        INTERP_METEOROLOGY, ZERO_SLOPE)
  DEALLOCATE (SEDIMENT_RESUSPENSION, ACTIVE_RULE_W2SELECTIVE)   !HYDRO_PLOT, CONSTITUENT_PLOT, DERIVED_PLOT,        
  DEALLOCATE (ORGPLD, ORGPRD, ORGPLP, ORGPRP, ORGNLD, ORGNRD, ORGNLP)
  DEALLOCATE (ICPL,TAVG,TAVGW,CAVG,CAVGW,CDAVG,CDAVGW) 
  DEALLOCATE (ORGNRP,KG_H2O_CONSTANT)
  DEALLOCATE  (PRINT_MACROPHYTE, MACROPHYTE_CALC,MACWBC,CONV2)
  DEALLOCATE  (MAC, MACRC,MACT, MACRM, MACSS)
  DEALLOCATE  (MGR,MMR, MRR)
  DEALLOCATE  (SMACRC, SMACRM)
  DEALLOCATE  (SMACT, SMAC)
  DEALLOCATE  (MT1,MT2,MT3,MT4,MK1,MK2,MK3,MK4,MG,MR,MM)
  DEALLOCATE  (MP, MN, MC,PSED,NSED,MHSP,MHSN,MHSC,MSAT)
  DEALLOCATE  (CDDRAG,KTICOL,ARMAC,MACWBCI, ANORM, DWV, DWSA)
  DEALLOCATE  (MBMP,MMAX,MPOM,LRPMAC,O2MR,O2MG)
  DEALLOCATE  (MACMBRS, MACMBRT,SSMACMB)
  DEALLOCATE  (CW, BIC)
  DEALLOCATE  (MACTRMR, MACTRMF,MACTRM)
  DEALLOCATE  (MLFPR)
  DEALLOCATE  (MLLIM, MPLIM,MCLIM,MNLIM)
  DEALLOCATE  (GAMMAJ)
  DEALLOCATE (POR,VOLKTI,VOLI,VSTEM,VSTEMKT,SAREA)
  DEALLOCATE (IWIND) ! MLM 08/12/05
  DEALLOCATE (ZG,ZM,ZEFF,PREFP,ZR,ZOOMIN,ZS2P,EXZ,ZT1,ZT2,ZT3,ZT4,ZK1,ZK2)
  DEALLOCATE (LDOMPMP,LDOMNMP,LPOMPMP,LPOMNMP,RPOMPMP,RPOMNMP,O2ZR) ! MLM 06/10/06
  DEALLOCATE (MPRWBC)                                               ! MLM 06/10/06
  DEALLOCATE (EXM)                                                  ! MLM 06/10/06
  DEALLOCATE (USTARBTKE,E,EROUGH, ARODI, STRICK, TKELATPRDCONST,AZT,DZT)
  DEALLOCATE(FIRSTI, LASTI, TKELATPRD, STRICKON, WALLPNT, IMPTKE, TKEBC)
  DEALLOCATE (ZK3,ZK4,ZP,ZN,ZC,PREFA,ZMU,TGRAZE,ZRT,ZMT,ZOORM,ZOORMR,ZOORMF,ZSR, ZS) ! POINTERS ,ZOO,ZOOSS,   SW 1/28/2019
  DEALLOCATE (LPZOOOUT,LPZOOIN,PO4ZR,NH4ZR,DOZR,TICZR,AGZ,AGZT)
  DEALLOCATE (GTIC,BGTO,EGTO)   ! CB 8/13/2010
  DEALLOCATE (INTERP_GATE)                     ! CB 8/13/2010  
  DEALLOCATE (ZGZ,PREFZ) !OMNIVOROUS ZOOPLANKTON
  DEALLOCATE (LPZOOINP,LPZOOINN,LPZOOOUTP,LPZOOOUTN)
  DEALLOCATE (SEDC, SEDN, SEDP,SEDNINFLUX, SEDPINFLUX, PFLUXIN,NFLUXIN)
  DEALLOCATE (SEDVPC, SEDVPP, SEDVPN)
  DEALLOCATE (SDKV,SEDDKTOT)
  DEALLOCATE (CBODS,KFJW)
  DEALLOCATE(BSAVE, GMA1,BTA1,ATMDEP_P,ATMDEP_N)
  DEALLOCATE(C_ATM_DEPOSITION, ATM_DEPOSITION,ATM_DEPOSITIONC,ATM_DEP_LOADING,ATM_DEPOSITION_INTERPOLATION,ATMDEPFN,ATMDCN,NACATD)
  DEALLOCATE(TN_SEDSOD_NH4,TP_SEDSOD_PO4,TPOUT,TPTRIB,TPDTRIB,TPWD,TPPR,TPIN,TNOUT,TNTRIB,TNDTRIB,TNWD,TNPR,TNIN,NH3GASLOSS)   ! TP_SEDBURIAL,TN_SEDBURIAL,
  IF(NBOD > 0)DEALLOCATE(NBODC,NBODN,NBODP)
  !IF(MWB_EXIST)THEN
  DEALLOCATE (WAIT_FOR_TRIB_INPUT,  WAIT_FOR_BRANCH_INPUT, TR_FILEDIR, BR_FILEDIR)                                      !SR 11/26/19
  !ENDIF
  IF (WAIT_FOR_INFLOW_RESULTS) THEN                                                                                     !SR 11/26/19
    DEALLOCATE (WAIT_TYPE, WAIT_INDEX, FILEDIR)                                                                         !SR 11/26/19
    CLOSE (9911)                                                                                                        !SR 11/26/19
  END IF

IF(FISHBIO)THEN
    DEALLOCATE (BIOEXPFN,WEIGHTNUM,C2ZOO,VOLROOS,C2W,IBIO)
    DEALLOCATE (BIOD, BIOF,BIODP)
ENDIF

  CALL DEALLOCATE_TIME_VARYING_DATA
  CALL DEALLOCATE_TRANSPORT
  IF(CONSTITUENTS)CALL DEALLOCATE_KINETICS
  CALL DEALLOCATE_WATERBODY
  CALL DEALLOCATE_PIPE_FLOW
  CALL DEALLOCATE_OPEN_CHANNEL
  IF(CONSTITUENTS .AND. AERATEC  == '      ON')CALL DEALLOCATE_AERATE
  IF(SELECTC == '      ON')CALL DEALLOCATE_SELECTIVE
  IF(SELECTC == '    USGS')CALL DEALLOCATE_SELECTIVEUSGS
  DEALLOCATE(WBSEG, PALT_JW, ELWS_INI, GTTYP, GTPC)                                               ! systdg 
  IF (SYSTDG) CALL DEALLOCATE_SYSTDG                                                              ! systdg
  IF (TDGTA) CLOSE (targetfnno)                                                                   ! systdg TDGtarget
  IF(TDGTA) CALL DEALLOCATE_TDGtarget                                                             ! systdg TDGtarget
  IF(CONSTITUENTS .AND. PHBUFF_EXIST)THEN
      DEALLOCATE(SDENI,PKI,PKSD)
    IF (OM_BUFFERING) THEN
     IF (OMTYPE == '    DIST') THEN
       DEALLOCATE (SDEN,PK,FRACT)
     ELSE
       DEALLOCATE (SDEN,PK)
     ENDIF
    ENDIF
  ENDIF
  If(CEMARelatedCode) Then         
        CALL DEALLOCATE_CEMA
	End If
    
  RETURN
  
  END SUBROUTINE ENDSIMULATION
