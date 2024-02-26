!***********************************************************************
!***********************************************************************
!***********************************************************************
!***                                                                ****
!***        N U M E R I C A L   F I S H   S U R R O G A T E         ****
!***                                                                ****
!***********************************************************************
!***********************************************************************
!***********************************************************************
! Andy Goodwin F77 Version 1/2001
! S. Wells F90 Version + Enhancements for Random Water Movement + Bug Fixes + etc 1/2001
    ! MANY CHANGES SW 5/1/2015
! Module Definitions

Module Fishy
      Integer, PARAMETER :: FPARA=16,FFP=4,WQP=2    ! The number of fish parameters. The number of Flow Field Parameters evaluated &
                                     !   the number of Water Quality Parameters evaluated at nodes for use in the NFS
      Integer NFISH,NUMGILLNETS,NUMACOUSTICS,NDT,JJ, JBP, KLAST,NL,NR
      DOUBLE PRECISION VALUE,GRADXVEL,GRADZVEL,GRADTEMP,GRADDO
      DOUBLE PRECISION XGRADHV,ZGRADHV,XGRADVV,ZGRADVV
      DOUBLE PRECISION XGRADDO,ZGRADDO,XGRADTP,ZGRADTP
      DOUBLE PRECISION FXVEL,FZVEL,FISHTEMP,FISHDO,TEMPORARYGRAD
      DOUBLE PRECISION FXHVCOUNT,FXVVCOUNT,FXDOCOUNT,FXTPCOUNT,FXRDCOUNT,TOTFXWGT
      DOUBLE PRECISION FZHVCOUNT,FZVVCOUNT,FZDOCOUNT,FZTPCOUNT,FZRDCOUNT,TOTFZWGT
      DOUBLE PRECISION FX00COUNT,FZ00COUNT,nfsfreq,RUNDIFF
      REAL ::   wmax,zmin,zmax  ! SW 2/16/01
      REAL            FXLOC,FYLOC,FZLOC,XREFL,YREFL,ZBOTREFL,ZSURREFL
      REAL            UFISH,VFISH,WFISH,FYVEL,SEDVEL,ALPHAX,ALPHAZ   
      REAL            OUTFREQ,JDAYDIFF,TSETP2MULT 
      REAL            FROMKTBOT,RRR,RMAT,RDX,RDY,RDZ,FYLOCTMP
      REAL            FSIZE,FAGE,LJDAY
      REAL            OUTFREQJAN,OUTFREQFEB,OUTFREQMAR,OUTFREQAPR
      REAL            OUTFREQMAY,OUTFREQJUN,OUTFREQJUL,OUTFREQAUG
      REAL            OUTFREQSEP,OUTFREQOCT,OUTFREQNOV,OUTFREQDEC
      REAL            HR,MILTIME,ASPRATIO,TEMPOPTD,TEMPOPTN
      REAL            FXSENSPH,BXSENSPH,UZSENSPH,DZSENSPH
      REAL            FBDYSEARCH,BBDYSEARCH,UBDYSEARCH,DBDYSEARCH
      REAL            XLOK,ZLOK,DIFFROMOPT,TEMPOPT
      REAL            MAXUFISH,MAXWFISH,MXXSPDL,MXZSPDL,DOTHRES
      REAL            FXGRADHV,FXGRADVV,FXGRADDO,FXGRADTP,FRDX,FRDY
      REAL            FZGRADHV,FZGRADVV,FZGRADDO,FZGRADTP,FRDZ
      REAL            HVXWEIGT,VVXWEIGT,DOXWEIGT,TPXWEIGT,RDXWEIGT,RDYWEIGT
      REAL            HVZWEIGT,VVZWEIGT,DOZWEIGT,TPZWEIGT,RDZWEIGT
      REAL            XURGENCY,YURGENCY,ZURGENCY
      REAL            CURRTFISHX,CURRTFISHZ,OTHERFISHX,OTHERFISHZ,FNEARX,FNEARZ
      REAL            XSCHZONEINN,ZSCHZONEINN,XSCHZONEMID,ZSCHZONEMID
      REAL            XSCHZONEOUT,ZSCHZONEOUT,XINZONE,ZINZONE,XMDZONE,ZMDZONE
      REAL            XOTZONE,ZOTZONE,EPSILONRD,LRUNDAY
      REAL            TOPTEMP,BOTTEMP,BOTK,MIDKTEMP,MIDKT,LFTTEMP,RGTTEMP
      REAL            MIDITEMP,MIDIT,TOPDO,BOTDO,MIDKDO,MIDKD,LFTDO,RGTDO
      REAL            MIDIDO,MIDID,TOPHVEL,BOTHVEL,BOTKK,MIDKHVEL,MIDKH
      REAL            LFTHVEL,RGTHVEL,MIDIHVEL,MIDIH,TOPVVEL,BOTVVEL
      REAL            MIDKVVEL,MIDKV,LFTVVEL,RGTVVEL,MIDIVVEL,MIDIV,VVELCAP
      REAL            SKYNIGHT,SKYDAWN,SKYDAY,SKYDUSK,COUNTFORTRUC
      REAL            XSEPCOMFORT,ZSEPCOMFORT,XDISPCOEFF,ZDISPCOEFF
      REAL            MAXREACTXHV,MINREACTXHV,DISTREACTXHV,MAXREACTZHV,MINREACTZHV,DISTREACTZHV
      REAL            MAXREACTXVV,MINREACTXVV,DISTREACTXVV,MAXREACTZVV,MINREACTZVV,DISTREACTZVV
      REAL            MAXREACTXDO,MINREACTXDO,DISTREACTXDO,MAXREACTZDO,MINREACTZDO,DISTREACTZDO
      REAL            MAXREACTXTP,MINREACTXTP,DISTREACTXTP,MAXREACTZTP,MINREACTZTP,DISTREACTZTP
      REAL            XURGTRAK,ZURGTRAK
      REAL            XHVTRUNC,XHVTRUCSUM,ZHVTRUNC,ZHVTRUCSUM,XVVTRUNC,XVVTRUCSUM,ZVVTRUNC,ZVVTRUCSUM
      REAL            XDOTRUNC,XDOTRUCSUM,ZDOTRUNC,ZDOTRUCSUM,XTPTRUNC,XTPTRUCSUM,ZTPTRUNC,ZTPTRUCSUM
      REAL            FXREACTVEL,FXREACTTMP,FXREACTDO,FXREACTRND
      REAL            FZREACTVEL,FZREACTTMP,FZREACTDO,FZREACTRND
      REAL            FXHVSPAN,FXVVSPAN,FXDOSPAN,FXTPSPAN,FXRDSPAN,TOTFXSPAN
      REAL            FZHVSPAN,FZVVSPAN,FZDOSPAN,FZTPSPAN,FZRDSPAN,TOTFZSPAN
      REAL            OLDFXLOC,OLDFYLOC,OLDFZLOC,DEPTHINT
      REAL            HADEPTHINT
      REAL            TEMPTHRES,TSTEP1,TSTEP1MULT,TSTEP2,TSTEP2MULT,DOTHRES2,DOSTEP1MULT,DELAYDATE

      INTEGER         FCOUNT,ISWITCH,FIMP,FKMP,FNBP,SURFCALC,BOTVVVEL
      INTEGER         KTWBF
      INTEGER         RD,SEED,TAG,FJR
      INTEGER         UNBP,DNBP,FN,WBRUN,NZONES
      INTEGER         NGRPFH,FI,FK,DIR,VARIABLE,LOCATE
      INTEGER         ILOK,KLOK,FSCHL,OTHFK,OTHFI,VARYTEMP,VARYDO,FKMPTEMP
      INTEGER         VARYHVEL,VARYVVEL,FATPT
      INTEGER         FSNAG,NET,SND
      integer         nfishpcel,nfishseg,ncollector    ! SW 1/16/01

      CHARACTER*2     AMPM
      CHARACTER*4     TXTNULLWQ,TXTNULLFF,TXTSTIMRULE,TXTVELORULE,&
                      TXTTEMPRULE,TXTDORULE,TXTRANDRULE,TXTPASSTRAN,&
                      TXTSCHOOLG,TXTMULTIRULE,txtparticle
      DIMENSION       RMAT(4)
      DIMENSION       FXVEL(5),FZVEL(5),FISHTEMP(5),FISHDO(5)
      DIMENSION       GRADXVEL(4),GRADZVEL(4),GRADTEMP(4),GRADDO(4)

      LOGICAL         DEBUG,LINEAR,WBSKIP,PREVENTBRCHSWITCH,SHOWSKY
      LOGICAL         STIMULIRULES,PASSIVETRANSPORT,RANDOMIZATION
      LOGICAL         VELOCITYRULES,TEMPRULES,DORULES,SCHOOLING
      LOGICAL         NULLFIELDWQ,WQNULLLINEAR,NULLFIELDFF,FFNULLLINEAR
      LOGICAL         MULTIPLERESPONSE,GILLNETSAMPLING,ACOUSTICSAMPLING,particle,line
      character*3     fishon,collector,partmod         ! SW 2/16/01

      real, allocatable, dimension (:,:) :: gillnets, HAOPERAT,OLDHALOC
      real, allocatable, dimension (:,:,:) :: flowfield, wqfield,sndcatch,nodes,lastfield
      integer, allocatable, dimension (:,:,:) :: netcatch
      integer, allocatable, dimension (:,:) :: rimpbr,limpbr
      real, allocatable, dimension (:,:) :: fishes
      integer, allocatable, dimension  (:) :: SNAGCOUNT,SNDCOUNT,nbrf
      integer, allocatable, dimension (:,:,:) :: corners
      integer, allocatable, dimension (:,:) :: ndinfo
      integer, allocatable, dimension  (:) :: ifish,ifisht,ifishb  ! SW 1/16/01
      integer, allocatable, dimension  (:) :: icoll,icollt,icollb  ! SW 2/16/01
      real, allocatable, dimension (:,:,:) :: BRCHFISH
      
      INTEGER :: DIAGFN=8001, DATADEBUGFN=8002, BARCHRTXFN=8003,BARCHRTZFN=8004,GILFN=8005,HYACFN=8006,FINALFN=8007

End Module Fishy


!  This Subroutine Calls the Following Subroutines:
!      FIMPBR,GRIDPLOT,WHATJR,SPLINE,RANDOM,VGILLNETS,ACOUSTICS,
!      FINDNEWBR,TAG124578,TAG369,FISHPLOT,INTERCONST,INTERFLOWF


      SUBROUTINE FISH

     USE SURFHE; Use Fishy; Use GDAYC;  Use SCREENC; Use GEOMC; USE GLOBAL
     IMPLICIT NONE
     
     REAL :: DZ
     INTEGER :: JF,KK,N

  IF (NIT.EQ.0) THEN             ! The very first time the Subroutine FISH is called

!***** Open 'Numerical Fish Surrogate' files
      Call Read_Fish_Data
      if(fishon.eq.'OFF')go to 97
      OPEN (DIAGFN,FILE='DIAGNOSTICS.OUT',STATUS='UNKNOWN')                         !FISH
      OPEN (DATADEBUGFN,FILE='DATADEBUG.OUT',STATUS='UNKNOWN')                           !FISH
      OPEN (BARCHRTXFN,FILE='BarChtX.dat',STATUS='UNKNOWN')                             !FISH
      OPEN (BARCHRTZFN,FILE='BarChtZ.dat',STATUS='UNKNOWN')
      open (FINALFN,file='finalfish.dat',status='unknown')                             !FISH

! Allocate arrays
! Read input data files

  Allocate(FISHES(NFISH,FPARA),WQFIELD(KMX,IMX,WQP),FLOWFIELD(KMX,IMX,FFP),NETCATCH(NUMGILLNETS,NFISH,5),&
           SNDCATCH(NUMACOUSTICS,3*NFISH,2),RIMPBR(IMX,2),LIMPBR(IMX,2), NODES(KMX,IMX,2),&
           NDINFO(IMX*KMX+IMX+KMX-1,4),CORNERS(NBR,IMX*KMX+IMX+KMX-1,4),LASTFIELD(KMX,IMX,2),&
           BRCHFISH(NBR,NFISH,FPARA),NBRF(NBR))
  flowfield=0.0
  wqfield=0.0
  corners=0
  ndinfo=0
  nodes=0.0
  netcatch=0
  sndcatch=0.0
  lastfield=0.0
  brchfish=0.0
  !rimpJB=0
  !limpJB=0
  nbrf=0
  fishes=0.0
        CALL FIMPBR                  ! Establishes data set of what segments are in each Branch
        SEED           = 92          ! Seed for the Random Number Generator


!------------------------------------------ FISH DATA INITIALIZATION ------------------------------------------

!Multiple Fish
        NGRPFH = 0                                       ! NGRPFH  = The number of fish being initiated
        DO jj=1,nfishseg
!        DO 140 FI=112,135,1                              ! Fish will be placed at each one of these Segment nodes
          DO FK=ifisht(jj),ifishb(jj)                 !KTWB(2),(KB(FI)-2),1   ! Fish will be placed at each one of these Layer nodes
            DO FATPT=1,nfishpcel                     ! The Number of Fish that will be placed at each drop location
              NGRPFH = NGRPFH + 1                        ! NGRPFH  = The number of fish being initiated
              IF (NGRPFH.GT.NFISH) THEN                  ! CHECK
                WRITE(DATADEBUGFN,*) 'ERROR: Must increase NFISH in fish.npt (Overestimate if necessary)'
                WRITE(DATADEBUGFN,9450) NGRPFH,NFISH
 9450           FORMAT(' NGRPFH =',I8,'     NFISH=',I8)
                STOP
              END IF
              FISHES(NGRPFH,1) = ifish(jj)    !FI        ! FIMP    = Initial segment IMP where fish is released
              FISHES(NGRPFH,3) = FK                      ! FKMP    = Initial layer KMP where fish is released
            end do
           end do
!  140   CONTINUE
        end do
        IF (NGRPFH.NE.NFISH) THEN                        ! CHECK
          WRITE(DATADEBUGFN,*) 'ERROR: NGRPFH in Subroutine FISH not equal to NFISH in fish.npt'
          WRITE(DATADEBUGFN,9450) NGRPFH,NFISH
          STOP
        END IF
!_____________________________________________________
!Data for All Fish
        klast=fishes(1,3)
        DO FN=1,NFISH
         IF(.NOT.LINE)THEN
          FISHES(FN,4)  = FZLOC   !0.0   ! FZLOC   = Location of fish within layer KMP from top side
         ELSE
           if(fn.eq.1.or.fishes(fn,3).ne.real(klast))then
            FISHES(FN,4)  = 0.           ! FZLOC   = assumed to be zero for LINE
            fimp=fishes(fn,1)
            call whatjr
            dz=h(fishes(fn,3),fjr)/nfishpcel
            klast=fishes(fn,3)
           else
            fishes(fn,4)=fishes(fn-1,4)+dz          
           endif 
         ENDIF
          FISHES(FN,2)  = FXLOC   !0.0   ! FXLOC   = Location of fish within segment IMP from upstream side
          FISHES(FN,5)  = B(INT(FISHES(FN,3)),INT(FISHES(FN,1)))*.5
                                     ! FYLOC   = Lateral fish release location (from left bank in plan view)
          call findbranch(fishes(fn,1))
          FISHES(FN,6)  = real(FNBP)   !6     ! FNBP    = Branch where fish is released
          FISHES(FN,7)  = FSIZE   !0.178 ! FSIZE   = Size (i.e., length) of fish in meters (1 inch = 0.0254 meters)
          FISHES(FN,8)  = FAGE   !2.0   ! FAGE    = Age of the fish
          FISHES(FN,9)  = UFISH   !0.0   ! UFISH   = Initial longitudinal velocity of the fish relative to water
          FISHES(FN,10) = VFISH   !0.0   ! VFISH   = Initial lateral velocity of the fish relative to water
          FISHES(FN,11) = WFISH   !0.0  ! WFISH   = Initial vertical velocity of the fish relative to water
          FISHES(FN,12) = 0          ! ISWITCH = Is fish still in system: Yes=0  No=1
          FISHES(FN,13) = 0          ! FSNAG   = Gillnet # fish is snagged in; = 0 if fish not in a gillnet
        END DO
      ELSE
  END IF

!*********************************************************************************************
!******************************  NO MORE INPUT BELOW THIS POINT  *****************************
!*********************************************************************************************

  if(fishon.eq.'OFF')go to 97

!Initialize Parameters Used to Evaluate Implementation of the Stimuli-Response Rules,
!  Prepare Text for Output Purposes, and Output Input Data to DIAGNOSTICS.OUT File

!Set Counting and Averaging Variables (used for the current timestep only) to Zero

      FXHVSPAN   = 0                             ! Counts # of fish for which the Horizontal Velocity is Dominant in the X-direction
      FXVVSPAN   = 0                             ! Counts # of fish for which the Vertical Velocity is Dominant in the X-direction
      FXDOSPAN   = 0                             ! Counts # of fish for which the Dissolved Oxygen is Dominant in the X-direction
      FXTPSPAN   = 0                             ! Counts # of fish for which the Temperature is Dominant in the X-direction
      FXRDSPAN   = 0                             ! Counts # of fish for which the Randomization is Dominant in the X-direction
      TOTFXSPAN  = 0                             ! Sums Weights Used in making X-Movement Decisions
      FZHVSPAN   = 0                             ! Counts # of fish for which the Horizontal Velocity is Dominant in the Z-direction
      FZVVSPAN   = 0                             ! Counts # of fish for which the Vertical Velocity is Dominant in the Z-direction
      FZDOSPAN   = 0                             ! Counts # of fish for which the Dissolved Oxygen is Dominant in the Z-direction
      FZTPSPAN   = 0                             ! Counts # of fish for which the Temperature is Dominant in the Z-direction
      FZRDSPAN   = 0                             ! Counts # of fish for which the Randomization is Dominant in the Z-direction
      TOTFZSPAN  = 0                             ! Sums Weights Used in making Z-Movement Decisions

      IF (NIT.EQ.0) THEN

!Set Counting and Averaging Variables to Zero

        FX00COUNT     = 0                        ! Counts # of times No Gradient is Dominant in the X-direction
        FXHVCOUNT     = 0                        ! Counts # of times Horizontal Velocity is Dominant in the X-direction
        FXVVCOUNT     = 0                        ! Counts # of times Vertical Velocity is Dominant in the X-direction
        FXDOCOUNT     = 0                        ! Counts # of times Dissolved Oxygen is Dominant in the X-direction
        FXTPCOUNT     = 0                        ! Counts # of times Temperature is Dominant in the X-direction
        FXRDCOUNT     = 0                        ! Counts # of times Randomization is Dominant in the X-direction
        TOTFXWGT      = 0                        ! Sums All Weights Used in making X-Movement Decisions
        FZ00COUNT     = 0                        ! Counts # of times No Gradient is Dominant in the Z-direction
        FZHVCOUNT     = 0                        ! Counts # of times Horizontal Velocity is Dominant in the Z-direction
        FZVVCOUNT     = 0                        ! Counts # of times Vertical Velocity is Dominant in the Z-direction
        FZDOCOUNT     = 0                        ! Counts # of times Dissolved Oxygen is Dominant in the Z-direction
        FZTPCOUNT     = 0                        ! Counts # of times Temperature is Dominant in the Z-direction
        FZRDCOUNT     = 0                        ! Counts # of times Randomization is Dominant in the Z-direction
        TOTFZWGT      = 0                        ! Sums All Weights Used in making Z-Movement Decisions
        XHVTRUNC      = 0                        ! Counts # of times Scaled Horz Vel Gradient in X-dir must be truncated to 1.0
        XHVTRUCSUM    = 0                        ! Sums the Values of all Scaled Horz Vel Gradients that were Truncated in X-dir
        ZHVTRUNC      = 0                        ! Counts # of times Scaled Horz Vel Gradient in Z-dir must be truncated to 1.0
        ZHVTRUCSUM    = 0                        ! Sums the Values of all Scaled Horz Vel Gradients that were Truncated in Z-dir
        XVVTRUNC      = 0                        ! Counts # of times Scaled Vert Vel Gradient in X-dir must be truncated to 1.0
        XVVTRUCSUM    = 0                        ! Sums the Values of all Scaled Vert Vel Gradients that were Truncated in X-dir
        ZVVTRUNC      = 0                        ! Counts # of times Scaled Vert Vel Gradient in Z-dir must be truncated to 1.0
        ZVVTRUCSUM    = 0                        ! Sums the Values of all Scaled Vert Vel Gradients that were Truncated in Z-dir
        XDOTRUNC      = 0                        ! Counts # of times Scaled Diss Oxyg Gradient in X-dir must be truncated to 1.0
        XDOTRUCSUM    = 0                        ! Sums the Values of all Scaled Diss Oxyg Gradients that were Truncated in X-dir
        ZDOTRUNC      = 0                        ! Counts # of times Scaled Diss Oxyg Gradient in Z-dir must be truncated to 1.0
        ZDOTRUCSUM    = 0                        ! Sums the Values of all Scaled Diss Oxyg Gradients that were Truncated in Z-dir
        XTPTRUNC      = 0                        ! Counts # of times Scaled Temperature Gradient in X-dir must be truncated to 1.0
        XTPTRUCSUM    = 0                        ! Sums the Values of all Scaled Temperature Gradients that were Truncated in X-dir
        ZTPTRUNC      = 0                        ! Counts # of times Scaled Temperature Gradient in Z-dir must be truncated to 1.0
        ZTPTRUCSUM    = 0                        ! Sums the Values of all Scaled Temperature Gradients that were Truncated in Z-dir
        COUNTFORTRUC  = 0                        ! Counts total # of times Stimuli-Response Rules are implemented

!Preparing Text for TecPlot Animation and Output Files

        IF (NULLFIELDWQ.EQV..TRUE.) THEN
          TXTNULLWQ = '  ON'
        ELSE
          TXTNULLWQ = ' OFF'
        END IF
        IF (NULLFIELDFF.EQV..TRUE.) THEN
          TXTNULLFF = '  ON'
        ELSE
          TXTNULLFF = ' OFF'
        END IF
        IF (MULTIPLERESPONSE.EQV..TRUE.) THEN
          TXTMULTIRULE = '  ON'
        ELSE
          TXTMULTIRULE = ' OFF'
        END IF
        IF (VELOCITYRULES.EQV..TRUE.) THEN
          TXTVELORULE = '  ON'
        ELSE
          TXTVELORULE = ' OFF'
        END IF
        IF (TEMPRULES.EQV..TRUE.) THEN
          TXTTEMPRULE = '  ON'
        ELSE
          TXTTEMPRULE = ' OFF'
        END IF
        IF (DORULES.EQV..TRUE.) THEN
          TXTDORULE = '  ON'
        ELSE
          TXTDORULE = ' OFF'
        END IF
        IF (RANDOMIZATION.EQV..TRUE.) THEN
          TXTRANDRULE = '  ON'
        ELSE
          TXTRANDRULE = ' OFF'
        END IF
        IF ((VELOCITYRULES.EQV..FALSE.).AND.&
            (TEMPRULES.EQV..FALSE.)    .AND.&
            (DORULES.EQV..FALSE.)      .AND.&
            (RANDOMIZATION.EQV..FALSE.)) THEN      ! This Additional Test is Used to Prevent Confusion in Output and Also to Prevent
          STIMULIRULES = .FALSE.                   !   Division by Zero when Calculating Weighted Averages of Each Factors' Influence
        END IF
        IF (STIMULIRULES.EQV..TRUE.) THEN
          TXTSTIMRULE = '  ON'
        ELSE
          TXTSTIMRULE = ' OFF'
          TXTVELORULE = ' OFF'
          TXTTEMPRULE = ' OFF'
          TXTDORULE   = ' OFF'
          TXTRANDRULE = ' OFF'
        END IF
        IF (PASSIVETRANSPORT.EQV..TRUE.) THEN
          TXTPASSTRAN = '  ON'
        ELSE
          TXTPASSTRAN = ' OFF'
        END IF
        IF (SCHOOLING.EQV..TRUE.) THEN
          TXTSCHOOLG = '  ON'
        ELSE
          TXTSCHOOLG = ' OFF'
        END IF
        if(particle)then
        txtparticle='  ON'
        else
        txtparticle=' OFF'
        end if

!Output Selected Input Data to the DIAGNOSTICS.OUT File

        WRITE(DIAGFN,*) ' Numerical Fish Surrogate Diagnostic Output'
        write(DIAGFN,*) '   PARTICLE ONLY TRANSPORT:',txtparticle
        WRITE(DIAGFN,*) '   Null Water Quality Field :',TXTNULLWQ
        WRITE(DIAGFN,*) '   Null Hydrodynamic Field  :',TXTNULLFF
        WRITE(DIAGFN,*) '   Passive Transport      :',TXTPASSTRAN
        WRITE(DIAGFN,*) '   Schooling Behavior     :',TXTSCHOOLG
        WRITE(DIAGFN,*) '   Multi-Response Behavior:',TXTMULTIRULE
        WRITE(DIAGFN,*) '   Stimuli-Response Rules :',TXTSTIMRULE
        WRITE(DIAGFN,*) '     - Velocity Rules    ',TXTVELORULE
        WRITE(DIAGFN,*) '     - Temperature Rules ',TXTTEMPRULE
        WRITE(DIAGFN,*) '     - Diss Oxygen Rules ',TXTDORULE
        WRITE(DIAGFN,*) '     - Randomization     ',TXTRANDRULE
        WRITE(DIAGFN,*) '   '
        WRITE(DIAGFN,*) 'Numerical Fish Surrogate Timestep (JDAY)'
        WRITE(DIAGFN,*) '  NFSFREQ  = ',DLT/(real(NDT)*86400.)
        WRITE(DIAGFN,*) '   '
        WRITE(DIAGFN,*) 'Weights of Influence Parameters (Range: 0 --> 1)'
        WRITE(DIAGFN,*) '  Horizontal Velocity Parameters'
        WRITE(DIAGFN,*) '    HVXWEIGT = ',HVXWEIGT
        WRITE(DIAGFN,*) '    VVXWEIGT = ',VVXWEIGT
        WRITE(DIAGFN,*) '  Vertical Velocity Parameters'
        WRITE(DIAGFN,*) '    HVZWEIGT = ',HVZWEIGT
        WRITE(DIAGFN,*) '    VVZWEIGT = ',VVZWEIGT
        WRITE(DIAGFN,*) '  Temperature Parameters'
        WRITE(DIAGFN,*) '    TPXWEIGT = ',TPXWEIGT
        WRITE(DIAGFN,*) '    TPZWEIGT = ',TPZWEIGT
        WRITE(DIAGFN,*) '    TEMPTHRES     = ',TEMPTHRES
        WRITE(DIAGFN,*) '    TSTEP1        = ',TSTEP1
        WRITE(DIAGFN,*) '    TSTEP1MULT    = ',TSTEP1MULT
        WRITE(DIAGFN,*) '    TSTEP2        = ',TSTEP2
        WRITE(DIAGFN,*) '    TSTEP2MULT    = ',TSTEP2MULT
        WRITE(DIAGFN,*) '  Dissolved Oxygen Parameters'
        WRITE(DIAGFN,*) '    DOXWEIGT = ',DOXWEIGT
        WRITE(DIAGFN,*) '    DOZWEIGT = ',DOZWEIGT
        WRITE(DIAGFN,*) '    DOTHRES2      = ',DOTHRES2
        WRITE(DIAGFN,*) '    DOSTEP1MULT   = ',DOSTEP1MULT
        WRITE(DIAGFN,*) '  Random Displacement Parameters'
        WRITE(DIAGFN,*) '    RDXWEIGT  = ',RDXWEIGT
        WRITE(DIAGFN,*) '    RDYWEIGT  = ',RDYWEIGT
        WRITE(DIAGFN,*) '    RDZWEIGT  = ',RDZWEIGT
        WRITE(DIAGFN,*) '    EPSILONRD = ',EPSILONRD
        WRITE(DIAGFN,*) '   '
        WRITE(DIAGFN,*) 'Values Used to Scale the Influence Parameters&
      Gradients'
        WRITE(DIAGFN,*) '  Horizontal Velocity'
        WRITE(DIAGFN,*) '    X-direction'
        WRITE(DIAGFN,*) '      MAXREACTXHV   = ',MAXREACTXHV
        WRITE(DIAGFN,*) '      MINREACTXHV   = ',MINREACTXHV
        WRITE(DIAGFN,*) '      DISTREACTXHV  = ',DISTREACTXHV
        WRITE(DIAGFN,*) '    Z-direction'
        WRITE(DIAGFN,*) '      MAXREACTZHV   = ',MAXREACTZHV
        WRITE(DIAGFN,*) '      MINREACTZHV   = ',MINREACTZHV
        WRITE(DIAGFN,*) '      DISTREACTZHV  = ',DISTREACTZHV
        WRITE(DIAGFN,*) '  Vertical Velocity'
        WRITE(DIAGFN,*) '    X-direction'
        WRITE(DIAGFN,*) '      MAXREACTXVV   = ',MAXREACTXVV
        WRITE(DIAGFN,*) '      MINREACTXVV   = ',MINREACTXVV
        WRITE(DIAGFN,*) '      DISTREACTXVV  = ',DISTREACTXVV
        WRITE(DIAGFN,*) '    Z-direction'
        WRITE(DIAGFN,*) '      MAXREACTZVV   = ',MAXREACTZVV
        WRITE(DIAGFN,*) '      MINREACTZVV   = ',MINREACTZVV
        WRITE(DIAGFN,*) '      DISTREACTZVV  = ',DISTREACTZVV
        WRITE(DIAGFN,*) '  Temperature'
        WRITE(DIAGFN,*) '    X-direction'
        WRITE(DIAGFN,*) '      MAXREACTXTP   = ',MAXREACTXTP
        WRITE(DIAGFN,*) '      MINREACTXTP   = ',MINREACTXTP
        WRITE(DIAGFN,*) '      DISTREACTXTP  = ',DISTREACTXTP
        WRITE(DIAGFN,*) '    Z-direction'
        WRITE(DIAGFN,*) '      MAXREACTZTP   = ',MAXREACTZTP
        WRITE(DIAGFN,*) '      MINREACTZTP   = ',MINREACTZTP
        WRITE(DIAGFN,*) '      DISTREACTZTP  = ',DISTREACTZTP
        WRITE(DIAGFN,*) '  Dissolved Oxygen'
        WRITE(DIAGFN,*) '    X-direction'
        WRITE(DIAGFN,*) '      MAXREACTXDO   = ',MAXREACTXDO
        WRITE(DIAGFN,*) '      MINREACTXDO   = ',MINREACTXDO
        WRITE(DIAGFN,*) '      DISTREACTXDO  = ',DISTREACTXDO
        WRITE(DIAGFN,*) '    Z-direction'
        WRITE(DIAGFN,*) '      MAXREACTZDO   = ',MAXREACTZDO
        WRITE(DIAGFN,*) '      MINREACTZDO   = ',MINREACTZDO
        WRITE(DIAGFN,*) '      DISTREACTZDO  = ',DISTREACTZDO
        WRITE(DIAGFN,*) '   '
        WRITE(DIAGFN,*) 'SCHOOLING AND DISPERSION Parameters'
        WRITE(DIAGFN,*) '  XINZONE    = ',XINZONE
        WRITE(DIAGFN,*) '  ZINZONE    = ',ZINZONE
        WRITE(DIAGFN,*) '  XMDZONE    = ',XMDZONE
        WRITE(DIAGFN,*) '  ZMDZONE    = ',ZMDZONE
        WRITE(DIAGFN,*) '  XOTZONE    = ',XOTZONE
        WRITE(DIAGFN,*) '  ZOTZONE    = ',ZOTZONE
        WRITE(DIAGFN,*) '  XDISPCOEFF = ',XDISPCOEFF
        WRITE(DIAGFN,*) '  ZDISPCOEFF = ',ZDISPCOEFF
        WRITE(DIAGFN,*) '   '
        WRITE(DIAGFN,*) 'SENSORY SPHERE Parameters'
        WRITE(DIAGFN,*) '  FBDYSEARCH = ',FBDYSEARCH
        WRITE(DIAGFN,*) '  BBDYSEARCH = ',BBDYSEARCH
        WRITE(DIAGFN,*) '  UBDYSEARCH = ',UBDYSEARCH
        WRITE(DIAGFN,*) '  DBDYSEARCH = ',DBDYSEARCH
        WRITE(DIAGFN,*) '   '
        WRITE(DIAGFN,*) 'MISC Parameters '
        WRITE(DIAGFN,*) '  TEMPOPTD  = ',TEMPOPTD
        WRITE(DIAGFN,*) '  TEMPOPTN  = ',TEMPOPTN
        WRITE(DIAGFN,*) '  MXXSPDL  = ',MXXSPDL
        WRITE(DIAGFN,*) '  MXZSPDL  = ',MXZSPDL
        WRITE(DIAGFN,*) '  DOTHRES  = ',DOTHRES
        WRITE(DIAGFN,*) '  XREFL    = ',XREFL
        WRITE(DIAGFN,*) '  YREFL    = ',YREFL
        WRITE(DIAGFN,*) '  ZBOTREFL = ',ZBOTREFL
        WRITE(DIAGFN,*) '  ZSURREFL = ',ZSURREFL
        write(DIAGFN,*)' FISHES:FN,FISHES(FN,I)'
! initial fish output
    do jf=1,nfish
    write(DIAGFN,'(i7,1x,<fpara>(f10.3,1x))')jf,(fishes(jf,i),i=1,fpara)
    end do
    !close(DIAGFN)
! End fish output section

      END IF


!Specification of the Time Stepping used in Running the NFS (W2 may, at times, run at timesteps of down to 2 minutes!)

! COMPUTE DYNAMIC NFSFREQ SW 1/15/01

      nfsfreq=dlt/(real(ndt)*86400.)    ! SW

      IF (NIT.EQ.0) LRUNDAY = JDAY                            ! LRUNDAY is used in determining NFS run frequency
      RUNDIFF = JDAY - LRUNDAY                                ! RUNDIFF is used in determining NFS run frequency
      IF ((RUNDIFF.LT.NFSFREQ).AND.(NIT.NE.0)) THEN           ! Skip the Numerical Fish Surrogate Subroutine for this timestep
        GOTO 97
      ELSE                                                    ! Run the Numerical Fish Surrogate Subroutine for this timestep
      END IF


!Specification of when and how often to output information/data

      IF (NIT.EQ.0) LJDAY = JDAY                              ! LJDAY is used for TecPlot output frequency purposes
      JDAYDIFF = JDAY - LJDAY                                 ! JDAYDIFF is used for TecPlot output frequency purposes
          
      IF (IMON.EQ.1)   OUTFREQ = OUTFREQJAN
      IF (IMON.EQ.2)   OUTFREQ = OUTFREQFEB
      IF (IMON.EQ.3)   OUTFREQ = OUTFREQMAR
      IF (IMON.EQ.4)   OUTFREQ = OUTFREQAPR
      IF (IMON.EQ.5)   OUTFREQ = OUTFREQMAY
      IF (IMON.EQ.6)   OUTFREQ = OUTFREQJUN
      IF (IMON.EQ.7)   OUTFREQ = OUTFREQJUL
      IF (IMON.EQ.8)   OUTFREQ = OUTFREQAUG
      IF (IMON.EQ.9)   OUTFREQ = OUTFREQSEP
      IF (IMON.EQ.10)  OUTFREQ = OUTFREQOCT
      IF (IMON.EQ.11)  OUTFREQ = OUTFREQNOV
      IF (IMON.EQ.12)  OUTFREQ = OUTFREQDEC
      
     ! IF (JDAY.GT.366) THEN
     !   WRITE(DATADEBUGFN,*) 'ERROR: Must modify OUTFREQ calculations in Subroutine FISH'
     !   STOP
     ! ELSE
     ! END IF
      IF ((JDAYDIFF.GE.OUTFREQ).AND.(OUTFREQ.NE.99)) THEN     ! Set FCOUNT = 1, so information will be outputted
        FCOUNT = 1                                            ! Information outputted iff FCOUNT = 1
        LJDAY = JDAY
      ELSE
        FCOUNT = 0                                            ! Set FCOUNT = 0, so information won't be outputted
      END IF


!Calculate Time of Day/Night

      !IF (((JDAY-REAL(JDAYG)).LT.0).OR.&                       ! CHECK: To make sure JDAY remainder is always between
      !    ((JDAY-REAL(JDAYG)).GE.1)) THEN                     !   0 and 1
      !  WRITE(DATADEBUGFN,*) 'ERROR: (JDAY-JDAYG) Calculation Incorrect'
      !  STOP
      !ELSE
      !END IF
      MILTIME =  (JDAY-REAL(INT(JDAY)))*24 !(JDAY-REAL(JDAYG))*24                         ! Military Time - MILTIME used for NFS Calculations
      HR = MILTIME                                            ! Convert MILTIME to hour of day - HR used for TecPlot
      AMPM = 'am'                                             ! Is HR am or pm?
      IF (HR.GE.12) AMPM = 'pm'                               ! Is HR am or pm?
      IF (HR.GE.13) HR = HR - 12                              ! Value of HR after 12 noon
      IF (HR.LT.1) HR = HR + 12                               ! Value of HR if time is between 12 midnight and 1am


!Determining Optimum Temperature for Fish Species, Given the Time of Day
IF(.not.particle)then
      IF ((MILTIME.GT.SKYNIGHT).OR.&
          (MILTIME.LE.SKYDAWN)) THEN                          ! During the Night
        TEMPOPT = TEMPOPTN
      ELSE IF ((MILTIME.GT.SKYDAWN).AND.&
               (MILTIME.LE.SKYDAY)) THEN                      ! During the Morning
        TEMPOPT = (TEMPOPTD-TEMPOPTN)/(SKYDAY-SKYDAWN)*&
                  (MILTIME-SKYDAY) + TEMPOPTD
      ELSE IF ((MILTIME.GT.SKYDAY).AND.&
               (MILTIME.LE.SKYDUSK)) THEN                     ! During the Day
        TEMPOPT = TEMPOPTD
      ELSE IF ((MILTIME.GT.SKYDUSK).AND.&
               (MILTIME.LE.SKYNIGHT)) THEN                    ! During the Evening
        TEMPOPT = (TEMPOPTN-TEMPOPTD)/(SKYNIGHT-SKYDUSK)*&
                  (MILTIME-SKYDUSK) + TEMPOPTD
      ELSE                                                    ! CHECK
        WRITE(DATADEBUGFN,*) 'ERROR: MILTIME not Corresponding with Intervals Set for Determining Optimum Temperature for Fish'
        STOP
      END IF
ENDIF

!Determining What Virtual Gillnets are Active

      DO 98 NET=1,NUMGILLNETS
        IF ((JDAY.GE.(GILLNETS(NET,1)+GILLNETS(NET,3)/24)).AND.&
            (JDAY.LE.(GILLNETS(NET,2)+GILLNETS(NET,4)/24))) THEN
          GILLNETS(NET,11) = 1                                ! Turn Gillnet ON
        ELSE
          GILLNETS(NET,11) = 0                                ! Turn Gillnet OFF
        END IF
   98 CONTINUE


!Determining What Virtual Hydroacoustic Surveys are Active

      DO 17 SND=1,NUMACOUSTICS                                ! NUMACOUSTICS = Total # of Hydroacoustic Samples
        IF ((JDAY.GE.(HAOPERAT(SND,1)+HAOPERAT(SND,3)/24)).AND.&
            (JDAY.LE.(HAOPERAT(SND,2)+HAOPERAT(SND,4)/24))) THEN
          HAOPERAT(SND,13) = 1                                ! Turn Hydroacoustic Survey ON
          OLDHALOC(SND,1) = HAOPERAT(SND,5)                   ! Previous Location of Upstream End of HA Beam
          OLDHALOC(SND,2) = HAOPERAT(SND,6)                   ! Previous Location of Downstream End of HA Beam
          HAOPERAT(SND,5) = HAOPERAT(SND,5) + (JDAY-LRUNDAY)*&
                           (HAOPERAT(SND,7)*3600*24)          ! Updated Location of Upstream End of HA Beam
          HAOPERAT(SND,6) = HAOPERAT(SND,6) + (JDAY-LRUNDAY)*&
                           (HAOPERAT(SND,7)*3600*24)          ! Updated Location of Downstream End of HA Beam
        ELSE
          HAOPERAT(SND,13) = 0                                ! Turn Hydroacoustic Survey OFF
        END IF
   17 CONTINUE


!Call the Subroutine GRIDPLOT to set up the FE Grid for TecPlot and to Calc node information for Traffic Rules

      IF ((FCOUNT.EQ.1).OR.(NIT.EQ.0)) CALL GRIDPLOT          ! This prepares the FE Grid for output display and
                                                              !   calculates the Flow and WQ information at nodes


!Load Fish Information and Begin NFS Logic

      IF (JDAY .LT. DELAYDATE) GOTO 28

      CALL INTERCONST                                         ! Subroutine to interpolate Constituent values
      CALL INTERFLOWF                                         ! Subroutine to interpolate Flow Field values
  if(debug)then   ! SW 2/01/01
    write(DIAGFN,*)'JDAY=',JDAY,'NDT,NFSFREQ:',NDT,NFSFREQ,' DELAYDATE:',delaydate
  end if          ! SW 2/01/01
DO N=1,NDT   ! LOOP OVER NDT Time periods
      DO 131 FN=1,NFISH              ! This is the only loop where statements are not indented

      FIMP =    INT(FISHES(FN,1))    ! Segment IMP where fish is located
      FXLOC =       FISHES(FN,2)     ! Location of fish within segment IMP from upstream side
      FKMP =    INT(FISHES(FN,3))    ! Layer KMP where fish is located
      FZLOC =       FISHES(FN,4)     ! Location of fish within layer KMP from top side
      FYLOC =       FISHES(FN,5)     ! Location of fish laterally from left bank (in plan view)
      FNBP =    INT(FISHES(FN,6))    ! Branch where fish is located
      FSIZE =       FISHES(FN,7)     ! Size (i.e., length) of fish in meters
      FAGE =        FISHES(FN,8)     ! Age of the fish
      UFISH =       FISHES(FN,9)     ! Longitudinal velocity of fish relative to water
      VFISH =       FISHES(FN,10)    ! Lateral velocity of fish relative to water
      WFISH =       FISHES(FN,11)    ! Vertical velocity of fish relative to water
      ISWITCH = INT(FISHES(FN,12))   ! Is fish still in system?: Yes=0  No=1
      FSNAG =   INT(FISHES(FN,13))   ! Gillnet # fish is snagged in; = 0 if fish not in a gillnet

      IF (ISWITCH.EQ.0) THEN         ! If fish has already left the system, then skip to the next fish

!Determine Horizontal & Vertical Velocity, Temperature, and DO Gradients in Fish Sensory Sphere:
!  How far to seek out gradients (i.e. edges of the fish sensory sphere given timestep used)?

        CALL RANDOM(SEED,RRR)                                           ! Subroutine which calculates a random number

        FXSENSPH = (FBDYSEARCH * FSIZE) * RRR                           ! FXSENSPH = Leading edge of the sensory sphere
        BXSENSPH = (BBDYSEARCH * FSIZE) * RRR                           ! BXSENSPH = Trailing edge of the sensory sphere
        UZSENSPH = (UBDYSEARCH * FSIZE) * RRR                           ! UZSENSPH = Upper edge of the sensory sphere
        DZSENSPH = (DBDYSEARCH * FSIZE) * RRR                           ! DZSENSPH = Lower edge of the sensory sphere

        CALL WHATJR                                                     ! What waterbody,JR, is the fish in?
        KTWBF = KTWB(FJR)                                                 ! Water Surface Layer of waterbody fish is located in

!  Find Flow and Water Quality Constituent Values at places of interest:

        DIR = 1                                               ! DIR = Direction Fish is swimming: 1=downstream, 2=upstream
        IF (UFISH.LT.0) DIR = -1                              ! UFISH is (+) when Fish is swimming downstream
        DO 142 VARIABLE=1,4                                   ! There are 4 Variables: 1=Horz Vel, 2=Vert Vel, 3=Temp, 4=Diss Oxy
          DO 143 LOCATE=1,5                                   ! There are 5 Locations around Fish of interest-See Subroutine SPLINE
            ILOK = FIMP                                       ! Segment I where Variable Value is desired
            XLOK = FXLOC                                      ! Location within Segment I where Variable Value is desired
            KLOK = FKMP                                       ! Layer K where Variable Value is desired
            ZLOK = FZLOC                                      ! Location within Layer K where Variable Value is desired
            IF ((LOCATE.EQ.1).AND.(DIR.EQ.1)) XLOK=FXLOC+FXSENSPH       ! Location forward of Fish when Fish swimming DOWNstream
            IF  (LOCATE.EQ.2)                 ZLOK=FZLOC+DZSENSPH       ! Location below Fish
            IF ((LOCATE.EQ.3).AND.(DIR.EQ.1)) XLOK=FXLOC-BXSENSPH       ! Location behind Fish when Fish swimming DOWNstream
            IF  (LOCATE.EQ.4)                 ZLOK=FZLOC-UZSENSPH       ! Location above Fish
            IF ((LOCATE.EQ.1).AND.(DIR.EQ.-1)) XLOK=FXLOC-FXSENSPH      ! Location forward of Fish when Fish swimming UPstream
            IF ((LOCATE.EQ.3).AND.(DIR.EQ.-1)) XLOK=FXLOC+BXSENSPH      ! Location behind Fish when Fish swimming UPstream
            CALL SPLINE                                                 ! Bilinear Spline Interpolation Subroutine
            IF (VARIABLE.EQ.1) FXVEL(LOCATE)     = VALUE                ! Horiz Vel at Location of Interest: (+) is Downstream !!!
            IF (VARIABLE.EQ.2) FZVEL(LOCATE)     = VALUE                ! Vert Vel at Location of Interest:  (+) is Downward !!!
            IF (VARIABLE.EQ.3) FISHTEMP(LOCATE)  = VALUE                ! Temperature at Location of Interest
            IF (VARIABLE.EQ.4) FISHDO(LOCATE)    = VALUE                ! Dissolved Oxygen at Location of Interest
  143     CONTINUE
  142   CONTINUE


! Determine lateral velocity based on lateral withdrawals or connected branches  SW 2/01/01
     call lateral_velocity
!        FYVEL = 0.0                                                     ! Lateral Velocity


!  Calculate (Linear) Gradients within Fish Sensory Sphere:
IF(.not.particle)then
!    Horizontal Velocity
        GRADXVEL(1) = (ABS(FXVEL(1))-ABS(FXVEL(5)))/FXSENSPH            ! Location=1 ==> Forward of Fish
        GRADXVEL(2) = (ABS(FXVEL(2))-ABS(FXVEL(5)))/DZSENSPH            ! Location=2 ==> Below Fish
        GRADXVEL(3) = (ABS(FXVEL(3))-ABS(FXVEL(5)))/BXSENSPH            ! Location=3 ==> Behind Fish
        GRADXVEL(4) = (ABS(FXVEL(4))-ABS(FXVEL(5)))/UZSENSPH            ! Location=4 ==> Above Fish

!    Vertical Velocity
        GRADZVEL(1) = (ABS(FZVEL(1))-ABS(FZVEL(5)))/FXSENSPH            ! Location=1 ==> Forward of Fish
        GRADZVEL(2) = (ABS(FZVEL(2))-ABS(FZVEL(5)))/DZSENSPH            ! Location=2 ==> Below Fish
        GRADZVEL(3) = (ABS(FZVEL(3))-ABS(FZVEL(5)))/BXSENSPH            ! Location=3 ==> Behind Fish
        GRADZVEL(4) = (ABS(FZVEL(4))-ABS(FZVEL(5)))/UZSENSPH            ! Location=4 ==> Above Fish

!    Temperature
        GRADTEMP(1) = (FISHTEMP(1)-FISHTEMP(5))/FXSENSPH                ! Location=1 ==> Forward of Fish
        GRADTEMP(2) = (FISHTEMP(2)-FISHTEMP(5))/DZSENSPH                ! Location=2 ==> Below Fish
        GRADTEMP(3) = (FISHTEMP(3)-FISHTEMP(5))/BXSENSPH                ! Location=3 ==> Behind Fish
        GRADTEMP(4) = (FISHTEMP(4)-FISHTEMP(5))/UZSENSPH                ! Location=4 ==> Above Fish

!    Dissolved Oxygen
        GRADDO(1)   = (FISHDO(1)-FISHDO(5))/FXSENSPH                    ! Location=1 ==> Forward of Fish
        GRADDO(2)   = (FISHDO(2)-FISHDO(5))/DZSENSPH                    ! Location=2 ==> Below Fish
        GRADDO(3)   = (FISHDO(3)-FISHDO(5))/BXSENSPH                    ! Location=3 ==> Behind Fish
        GRADDO(4)   = (FISHDO(4)-FISHDO(5))/UZSENSPH                    ! Location=4 ==> Above Fish


!  Switch Gradients for Fish Facing Upstream (All Notation Below Assumes Fish Facing Downstream)

        IF (DIR.EQ.-1) THEN
          TEMPORARYGRAD = GRADXVEL(1)
          GRADXVEL(1)   = GRADXVEL(3)
          GRADXVEL(3)   = TEMPORARYGRAD
          TEMPORARYGRAD = GRADZVEL(1)
          GRADZVEL(1)   = GRADZVEL(3)
          GRADZVEL(3)   = TEMPORARYGRAD
          TEMPORARYGRAD = GRADTEMP(1)
          GRADTEMP(1)   = GRADTEMP(3)
          GRADTEMP(3)   = TEMPORARYGRAD
          TEMPORARYGRAD = GRADDO(1)
          GRADDO(1)     = GRADDO(3)
          GRADDO(3)     = TEMPORARYGRAD
        END IF


!**********************************************
!  F I S H   B E H A V I O R A L   L O G I C  *
!**********************************************

!Evaluate Gradients and Select Maximum Gradient for Each Factor (Horizontal Vel, Vertical Vel, DO, and Temp)
!  NOTE: Gradients are Defined Relative to Conditions at Current Fish Location

!Horizontal Velocity Gradients                                     ! (-) Sign used to archive what direction fish will want to move
!!  X-direction
!        IF (GRADXVEL(1).GT.GRADXVEL(3)) XGRADHV =  GRADXVEL(1)     ! Influence of X-direct'l Horz Vel Grad proportional to GRADXVEL(1)
!        IF (GRADXVEL(1).LT.GRADXVEL(3)) XGRADHV = -GRADXVEL(3)     ! Influence of X-direct'l Horz Vel Grad proportional to GRADXVEL(3)
!        IF (GRADXVEL(1).EQ.GRADXVEL(3))                XGRADHV = 0 ! Influence of X-direct'l Horiz Veloc Gradient is zero
!        IF ((GRADXVEL(1).LE.0).AND.(GRADXVEL(3).LE.0)) XGRADHV = 0 ! Influence of X-direct'l Horiz Veloc Gradient is zero
!  X-direction
        IF (GRADXVEL(1).GT.GRADXVEL(3)) XGRADHV =  ABS(FXVEL(1))   ! Influence of X-direct'l Horz Vel Grad proportional to FXVEL(5)
        IF (GRADXVEL(1).LT.GRADXVEL(3)) XGRADHV = -ABS(FXVEL(3))   ! Influence of X-direct'l Horz Vel Grad proportional to FXVEL(5)
        IF (GRADXVEL(1).EQ.GRADXVEL(3))                XGRADHV = 0 ! Influence of X-direct'l Horiz Veloc Gradient is zero
        IF ((GRADXVEL(1).LE.0).AND.(GRADXVEL(3).LE.0)) XGRADHV = 0 ! Influence of X-direct'l Horiz Veloc Gradient is zero

!  Z-direction
        IF (GRADXVEL(4).GT.GRADXVEL(2)) ZGRADHV = -GRADXVEL(4)     ! Influence of Z-direct'l Horz Vel Grad proportional to GRADXVEL(4)
        IF (GRADXVEL(4).LT.GRADXVEL(2)) ZGRADHV =  GRADXVEL(2)     ! Influence of Z-direct'l Horz Vel Grad proportional to GRADXVEL(2)
        IF (GRADXVEL(4).EQ.GRADXVEL(2))                ZGRADHV = 0 ! Influence of Z-direct'l Horiz Veloc Gradient is zero
        IF ((GRADXVEL(4).LE.0).AND.(GRADXVEL(2).LE.0)) ZGRADHV = 0 ! Influence of Z-direct'l Horiz Veloc Gradient is zero

!Vertical Velocity Gradients                                       ! (-) Sign used to archive what direction fish will want to move
!  X-direction
        IF (GRADZVEL(1).GT.GRADZVEL(3)) XGRADVV =  GRADZVEL(1)     ! Influence of X-direct'l Vert Vel Grad proportional to GRADZVEL(1)
        IF (GRADZVEL(1).LT.GRADZVEL(3)) XGRADVV = -GRADZVEL(3)     ! Influence of X-direct'l Vert Vel Grad proportional to GRADZVEL(3)
        IF (GRADZVEL(1).EQ.GRADZVEL(3))                XGRADVV = 0 ! Influence of X-direct'l Vertical Veloc Gradient is zero
        IF ((GRADZVEL(1).LE.0).AND.(GRADZVEL(3).LE.0)) XGRADVV = 0 ! Influence of X-direct'l Vertical Veloc Gradient is zero

!  Z-direction
        IF (GRADZVEL(4).GT.GRADZVEL(2)) ZGRADVV = -GRADZVEL(4)     ! Influence of Z-direct'l Vert Vel Grad proportional to GRADZVEL(4)
        IF (GRADZVEL(4).LT.GRADZVEL(2)) ZGRADVV =  GRADZVEL(2)     ! Influence of Z-direct'l Vert Vel Grad proportional to GRADZVEL(2)
        IF (GRADZVEL(4).EQ.GRADZVEL(2))                ZGRADVV = 0 ! Influence of Z-direct'l Vertical Veloc Gradient is zero
        IF ((GRADZVEL(4).LE.0).AND.(GRADZVEL(2).LE.0)) ZGRADVV = 0 ! Influence of Z-direct'l Vertical Veloc Gradient is zero

!Dissolved Oxygen Gradients                                        ! (-) Sign used to archive what direction fish will want to move
!  X-direction
        IF (GRADDO(1).GT.GRADDO(3)) XGRADDO =  GRADDO(1)           ! Influence of X-direct'l DO Grad proportional to GRADDO(1)
        IF (GRADDO(1).LT.GRADDO(3)) XGRADDO = -GRADDO(3)           ! Influence of X-direct'l DO Grad proportional to GRADDO(3)
        IF (GRADDO(1).EQ.GRADDO(3))                XGRADDO = 0     ! Influence of X-direct'l DO Gradient is zero
        IF ((GRADDO(1).LE.0).AND.(GRADDO(3).LE.0)) XGRADDO = 0     ! Influence of X-direct'l DO Gradient is zero

!  Z-direction
        IF (GRADDO(4).GT.GRADDO(2)) ZGRADDO = -GRADDO(4)           ! Influence of Z-direct'l DO Grad proportional to GRADDO(4)
        IF (GRADDO(4).LT.GRADDO(2)) ZGRADDO =  GRADDO(2)           ! Influence of Z-direct'l DO Grad proportional to GRADDO(2)
        IF (GRADDO(4).EQ.GRADDO(2))                ZGRADDO = 0     ! Influence of Z-direct'l DO Gradient is zero
        IF ((GRADDO(4).LE.0).AND.(GRADDO(2).LE.0)) ZGRADDO = 0     ! Influence of Z-direct'l DO Gradient is zero

!Temperature Gradients (Based on optimum preferred Temperature)      ! (-) Sign used to archive what direction fish will want to move
        DIFFROMOPT = FISHTEMP(5) - TEMPOPT                           !DIFFROMOPT = Temperature difference from Optimum Temperature
        IF (DIFFROMOPT.LT.0) THEN                                    !Temperature at Fish location is colder than optimum
!  X-direction
          IF (GRADTEMP(1).GT.GRADTEMP(3)) XGRADTP =  GRADTEMP(1)     ! Influence of X-direct'l Temp Grad proportional to GRADTEMP(1)
          IF (GRADTEMP(1).LT.GRADTEMP(3)) XGRADTP = -GRADTEMP(3)     ! Influence of X-direct'l Temp Grad proportional to GRADTEMP(3)
          IF ((GRADTEMP(1).LE.0).AND.(GRADTEMP(3).LE.0)) XGRADTP = 0 ! No X-direct'l Temp Grad infl. if colder forward and behind fish

!  Z-direction
          IF (GRADTEMP(4).GT.GRADTEMP(2)) ZGRADTP = -GRADTEMP(4)     ! Influence of Z-direct'l Temp Grad proportional to GRADTEMP(2)
          IF (GRADTEMP(4).LT.GRADTEMP(2)) ZGRADTP =  GRADTEMP(2)     ! Influence of Z-direct'l Temp Grad proportional to GRADTEMP(4)
          IF ((GRADTEMP(4).LE.0).AND.(GRADTEMP(2).LE.0)) ZGRADTP = 0 ! No Z-direct'l Temp Grad infl. if colder above and below fish
        ELSE IF (DIFFROMOPT.GT.0) THEN                               !Temperature at Fish location is warmer than optimum
!  X-direction
          IF (GRADTEMP(1).GT.GRADTEMP(3)) XGRADTP =  GRADTEMP(3)     ! Influence of X-direct'l Temp Grad proportional to GRADTEMP(3)
          IF (GRADTEMP(1).LT.GRADTEMP(3)) XGRADTP = -GRADTEMP(1)     ! Influence of X-direct'l Temp Grad proportional to GRADTEMP(1)
          IF ((GRADTEMP(1).GE.0).AND.(GRADTEMP(3).GE.0)) XGRADTP = 0 ! No X-direct'l Temp Grad infl. if warmer forward and behind fish

!  Z-direction
          IF (GRADTEMP(4).GT.GRADTEMP(2)) ZGRADTP = -GRADTEMP(2)     ! Influence of Z-direct'l Temp Grad proportional to GRADTEMP(4)
          IF (GRADTEMP(4).LT.GRADTEMP(2)) ZGRADTP =  GRADTEMP(4)     ! Influence of Z-direct'l Temp Grad proportional to GRADTEMP(2)
          IF ((GRADTEMP(4).GE.0).AND.(GRADTEMP(2).GE.0)) ZGRADTP = 0 ! No Z-direct'l Temp Grad infl. if warmer above and below fish
        ELSE                                                         ! Temperature at Fish location is optimum
          XGRADTP = 0                                                ! Influence of X-direct'l Temperature Gradient is zero
          ZGRADTP = 0                                                ! Influence of Z-direct'l Temperature Gradient is zero
        END IF
        IF (GRADTEMP(1).EQ.GRADTEMP(3))   XGRADTP = 0                ! Influence of X-direct'l Temperature Gradient is zero
        IF (GRADTEMP(4).EQ.GRADTEMP(2))   ZGRADTP = 0                ! Influence of Z-direct'l Temperature Gradient is zero


!Gradient Scaling Calculations

!Scaling of Horizontal Velocity Gradients (by a 'Practical' Gradient Likely to Induce the Maximum Fish Response)
!  X-direction
        XGRADHV = XGRADHV / ABS((MAXREACTXHV-MINREACTXHV)/DISTREACTXHV)
        IF (ABS(XGRADHV).GT.1) THEN
          XHVTRUNC = XHVTRUNC + 1                                    ! Tally # of times Scaled ABS(XGRADHV) must be Truncated to 1.0
          XHVTRUCSUM = XHVTRUCSUM + ABS(XGRADHV)                     ! Sum all ABS(XGRADHV) which must be Truncated, for later analysis
          IF (XGRADHV.GT.0) XGRADHV = 1                              ! Truncate XGRADHV to  1.0
          IF (XGRADHV.LT.0) XGRADHV = -1                             ! Truncate XGRADHV to -1.0
        END IF

!  Z-direction
        ZGRADHV = ZGRADHV / ABS((MAXREACTZHV-MINREACTZHV)/DISTREACTZHV)
        IF (ABS(ZGRADHV).GT.1) THEN
          ZHVTRUNC = ZHVTRUNC + 1                                    ! Tally # of times Scaled ABS(ZGRADHV) must be Truncated to 1.0
          ZHVTRUCSUM = ZHVTRUCSUM + ABS(ZGRADHV)                     ! Sum all ABS(ZGRADHV) which must be Truncated, for later analysis
          IF (ZGRADHV.GT.0) ZGRADHV = 1                              ! Truncate ZGRADHV to  1.0
          IF (ZGRADHV.LT.0) ZGRADHV = -1                             ! Truncate ZGRADHV to -1.0
        END IF


!Scaling of Vertical Velocity Gradients (by a 'Practical' Gradient Likely to Induce the Maximum Fish Response)
!  X-direction
        XGRADVV = XGRADVV / ABS((MAXREACTXVV-MINREACTXVV)/DISTREACTXVV)
        IF (ABS(XGRADVV).GT.1) THEN
          XVVTRUNC = XVVTRUNC + 1                                    ! Tally # of times Scaled ABS(XGRADVV) must be Truncated to 1.0
          XVVTRUCSUM = XVVTRUCSUM + ABS(XGRADVV)                     ! Sum all ABS(XGRADVV) which must be Truncated, for later analysis
          IF (XGRADVV.GT.0) XGRADVV = 1                              ! Truncate XGRADVV to  1.0
          IF (XGRADVV.LT.0) XGRADVV = -1                             ! Truncate XGRADVV to -1.0
        END IF

!  Z-direction
        ZGRADVV = ZGRADVV / ABS((MAXREACTZVV-MINREACTZVV)/DISTREACTZVV)
        IF (ABS(ZGRADVV).GT.1) THEN
          ZVVTRUNC = ZVVTRUNC + 1                                    ! Tally # of times Scaled ABS(ZGRADVV) must be Truncated to 1.0
          ZVVTRUCSUM = ZVVTRUCSUM + ABS(ZGRADVV)                     ! Sum all ABS(ZGRADVV) which must be Truncated, for later analysis
          IF (ZGRADVV.GT.0) ZGRADVV = 1                              ! Truncate ZGRADVV to  1.0
          IF (ZGRADVV.LT.0) ZGRADVV = -1                             ! Truncate ZGRADVV to -1.0
        END IF


!Scaling of Dissolved Oxygen Gradients (by a 'Practical' Gradient Likely to Induce the Maximum Fish Response)
!  X-direction
        XGRADDO = XGRADDO / ABS((MAXREACTXDO-MINREACTXDO)/DISTREACTXDO)
        IF (ABS(XGRADDO).GT.1) THEN
          XDOTRUNC = XDOTRUNC + 1                                    ! Tally # of times Scaled ABS(XGRADDO) must be Truncated to 1.0
          XDOTRUCSUM = XDOTRUCSUM + ABS(XGRADDO)                     ! Sum all ABS(XGRADDO) which must be Truncated, for later analysis
          IF (XGRADDO.GT.0) XGRADDO = 1                              ! Truncate XGRADDO to  1.0
          IF (XGRADDO.LT.0) XGRADDO = -1                             ! Truncate XGRADDO to -1.0
        END IF

!  Z-direction
        ZGRADDO = ZGRADDO / ABS((MAXREACTZDO-MINREACTZDO)/DISTREACTZDO)
        IF (ABS(ZGRADDO).GT.1) THEN
          ZDOTRUNC = ZDOTRUNC + 1                                    ! Tally # of times Scaled ABS(ZGRADDO) must be Truncated to 1.0
          ZDOTRUCSUM = ZDOTRUCSUM + ABS(ZGRADDO)                     ! Sum all ABS(ZGRADDO) which must be Truncated, for later analysis
          IF (ZGRADDO.GT.0) ZGRADDO = 1                              ! Truncate ZGRADDO to  1.0
          IF (ZGRADDO.LT.0) ZGRADDO = -1                             ! Truncate ZGRADDO to -1.0
        END IF


!Scaling of Temperature Gradients (by a 'Practical' Gradient Likely to Induce the Maximum Fish Response)
!  X-direction
        XGRADTP = XGRADTP / ABS((MAXREACTXTP-MINREACTXTP)/DISTREACTXTP)
        IF (ABS(XGRADTP).GT.1) THEN
          XTPTRUNC = XTPTRUNC + 1                                    ! Tally # of times Scaled ABS(XGRADTP) must be Truncated to 1.0
          XTPTRUCSUM = XTPTRUCSUM + ABS(XGRADTP)                     ! Sum all ABS(XGRADTP) which must be Truncated, for later analysis
          IF (XGRADTP.GT.0) XGRADTP = 1                              ! Truncate XGRADTP to  1.0
          IF (XGRADTP.LT.0) XGRADTP = -1                             ! Truncate XGRADTP to -1.0
        END IF

!  Z-direction
        ZGRADTP = ZGRADTP / ABS((MAXREACTZTP-MINREACTZTP)/DISTREACTZTP)
        IF (ABS(ZGRADTP).GT.1) THEN
          ZTPTRUNC = ZTPTRUNC + 1                                    ! Tally # of times Scaled ABS(ZGRADTP) must be Truncated to 1.0
          ZTPTRUCSUM = ZTPTRUCSUM + ABS(ZGRADTP)                     ! Sum all ABS(ZGRADTP) which must be Truncated, for later analysis
          IF (ZGRADTP.GT.0) ZGRADTP = 1                              ! Truncate ZGRADTP to  1.0
          IF (ZGRADTP.LT.0) ZGRADTP = -1                             ! Truncate ZGRADTP to -1.0
        END IF

        COUNTFORTRUC = COUNTFORTRUC + 1                              ! Tally # of times the above tests are implemented


!Modify Selected Gradients into Stepwise Functions

        IF (FISHTEMP(5) .GT. TEMPTHRES) THEN
          IF (XGRADTP .LT. 0) XGRADTP = -1
          IF (XGRADTP .GT. 0) XGRADTP =  1
          IF (ZGRADTP .LT. 0) ZGRADTP = -1
          IF (ZGRADTP .GT. 0) ZGRADTP =  1
        ELSE IF (ABS(DIFFROMOPT) .LT. TSTEP1) THEN
          XGRADTP = TSTEP1MULT * XGRADTP
          ZGRADTP = TSTEP1MULT * ZGRADTP
        ELSE IF (ABS(DIFFROMOPT) .LT. TSTEP2) THEN
          XGRADTP = TSTEP2MULT * XGRADTP
          ZGRADTP = TSTEP2MULT * ZGRADTP
        END IF
        IF (FISHDO(5) .GT. DOTHRES2) THEN
          XGRADDO = DOSTEP1MULT * XGRADDO
          ZGRADDO = DOSTEP1MULT * ZGRADDO
        END IF


!Calculation and Scaling of Random Displacement Values

        DO 104 RD=1,4                                                 ! Calculate 4 random numbers
          CALL RANDOM(SEED,RRR)                                       ! Subroutine which calculates a random number
          RMAT(RD) = RRR
  104   CONTINUE

        RDX = (2*RMAT(1)-1) * SQRT(ABS(UFISH)+EPSILONRD)              ! Calc Random Parameter in X-direction (Longitudinal)
        RDY = (2*RMAT(2)-1) * SQRT(ABS(VFISH)+EPSILONRD)              ! Calc Random Parameter in Y-direction (Lateral)
        RDZ = (2*RMAT(3)-1) * SQRT(ABS(WFISH)+EPSILONRD)              ! Calc Random Parameter in Z-direction (Vertical)

        RDX = RDX / SQRT(FSIZE*MXXSPDL+EPSILONRD)                     ! Scale Random Parameter to limit range to -1 <--> 1

! Check this - why only +1 and -1 ??? Scott Wells
        IF (RDX .LT. 0) THEN
          RDX = -1
        ELSE
          RDX = 1
        END IF
        RDY = RDY / SQRT(     1       +EPSILONRD)                     ! Scale Random Parameter to limit range to -1 <--> 1
        RDZ = RDZ / SQRT(FSIZE*MXZSPDL+EPSILONRD)                     ! Scale Random Parameter to limit range to -1 <--> 1


!Logic to Turn Different Numerical Fish Surrogate Options 'ON' or 'OFF'

        IF (VELOCITYRULES.EQV..FALSE.) THEN    ! If .TRUE., velocity stimuli-response rules contribute to the movement of fish
          HVXWEIGT = 0.0                       ! If .FALSE., velocity stimuli-response rules do not contribute to the movement of fish
          VVXWEIGT = 0.0
          HVZWEIGT = 0.0
          VVZWEIGT = 0.0
        END IF

        IF (TEMPRULES.EQV..FALSE.) THEN        ! If .TRUE., temp stimuli-response rules contribute to the movement of fish
          TPXWEIGT = 0.0                       ! If .FALSE., temp stimuli-response rules do not contribute to the movement of fish
          TPZWEIGT = 0.0
        END IF

        IF (DORULES.EQV..FALSE.) THEN          ! If .TRUE., DO stimuli-response rules contribute to the movement of fish
          DOXWEIGT = 0.0                       ! If .FALSE., DO stimuli-response rules do not contribute to the movement of fish
          DOZWEIGT = 0.0
        END IF

        IF (RANDOMIZATION.EQV..FALSE.) THEN    ! If .TRUE., random displacement contributes to the movement of fish
          RDXWEIGT = 0.0                       ! If .FALSE., random displacement does not contribute to the movement of fish
          RDYWEIGT = 0.0
          RDZWEIGT = 0.0
        END IF

        IF (STIMULIRULES.EQV..FALSE.) THEN ! If .TRUE., stimuli-response rules (excluding passive transport) contribute to fish movement
          HVXWEIGT = 0.0                   ! If .FALSE., stimuli-resp. rules (excl. passive transp) do not contribute to fish movement
          VVXWEIGT = 0.0
          DOXWEIGT = 0.0
          TPXWEIGT = 0.0
          RDXWEIGT = 0.0
          RDYWEIGT = 0.0
          HVZWEIGT = 0.0
          VVZWEIGT = 0.0
          DOZWEIGT = 0.0
          TPZWEIGT = 0.0
          RDZWEIGT = 0.0
        END IF

        IF (PASSIVETRANSPORT.EQV..FALSE.) THEN ! If .TRUE., passive transport contributes to the movement of fish
          FXVEL(5) = 0.0                       ! If .FALSE., passive transport does not contribute to the movement of fish
          FYVEL    = 0.0
          FZVEL(5) = 0.0
        END IF


!Calculation of Values Which Will Influence Fish Behavior/Movement

!  X-direction
        FXGRADHV = XGRADHV * HVXWEIGT    ! FXGRADHV  = Horiz Vel Gradient Value in X-direction the Fish uses to determine movement
                                         !  HVXWEIGT = Horiz Vel Weight for X-direction (Range: 0 --> 1)
        FXGRADVV = XGRADVV * VVXWEIGT    ! FXGRADVV  = Vert Vel Gradient Value in X-direction the Fish uses to determine movement
                                         !  VVXWEIGT = Vert Vel Weight for X-direction (Range: 0 --> 1)
        FXGRADDO = XGRADDO * DOXWEIGT    ! FXGRADDO  = DO Gradient Value in X-direction the Fish uses to determine movement
                                         !  DOXWEIGT = DO Weight for X-direction (Range: 0 --> 1)
        FXGRADTP = XGRADTP * TPXWEIGT    ! FXGRADTP  = Temperature Gradient Value in X-direction the Fish uses to determine movement
                                         !  TPXWEIGT = Temperature Weight for X-direction (Range: 0 --> 1)
        FRDX     =   RDX   * RDXWEIGT    ! FRDX      = X-directional Random Parameter that contributes to Fish movement
                                         !  RDXWEIGT = Random Movement Weight for X-direction (Range: 0 --> 1)

!  Y-direction
        FRDY     =   RDY   * RDYWEIGT    ! FRDY      = Y-directional Random Parameter that contributes to Fish movement
                                         !  RDYWEIGT = Random Movement Weight for Y-direction (Range: 0 --> 1)

!  Z-direction
        FZGRADHV = ZGRADHV * HVZWEIGT    ! FZGRADHV  = HorzVel Grad Value in Z-dir Fish uses to determine movement (Shear Zone Response)
                                         !  HVZWEIGT = Horiz Vel Weight for Z-direction (Range: 0 --> 1)
        FZGRADVV = ZGRADVV * VVZWEIGT    ! FZGRADVV  = Vert Vel Gradient Value in Z-direction the Fish uses to determine movement
                                         !  VVZWEIGT = Vert Vel Weight for Z-direction (Range: 0 --> 1)
        FZGRADDO = ZGRADDO * DOZWEIGT    ! FZGRADDO  = DO Gradient Value in Z-direction the Fish uses to determine movement
                                         !  DOZWEIGT = DO Weight for Z-direction (Range: 0 --> 1)
        FZGRADTP = ZGRADTP * TPZWEIGT    ! FZGRADTP  = Temperature Gradient Value in Z-direction the Fish uses to determine movement
                                         !  TPZWEIGT = Temperature Weight for Z-direction (Range: 0 --> 1)
        FRDZ     =   RDZ   * RDZWEIGT    ! FRDZ      = Z-directional Random Parameter that contributes to Fish movement
                                         !  RDZWEIGT = Random Movement Weight for Z-direction (Range: 0 --> 1)


!Calculate 'URGENCY' of Fish to Move Towards Better Environmental Conditions (-1.0 .LE. URGENCY .LE. 1.0)

        IF ((MULTIPLERESPONSE).AND.(STIMULIRULES.EQV..TRUE.)) THEN ! URGENCY Calculated as a Weighted Average of Selected
!  X-direction                                                     !   Scaled Gradients
          XURGENCY = (FXGRADHV+FXGRADVV+FXGRADDO+FXGRADTP+FRDX)/&     ! Weighted Average of Selected Scaled Gradients in X-direction
                     (HVXWEIGT+VVXWEIGT+DOXWEIGT+TPXWEIGT+RDXWEIGT)
!    Tallying the '% Influence' Over All Timesteps
          FXHVCOUNT = FXHVCOUNT + FXGRADHV                           ! Tally the weight of FXGRADHV in XURGENCY
          FXVVCOUNT = FXVVCOUNT + FXGRADVV                           ! Tally the weight of FXGRADVV in XURGENCY
          FXDOCOUNT = FXDOCOUNT + FXGRADDO                           ! Tally the weight of FXGRADDO in XURGENCY
          FXTPCOUNT = FXTPCOUNT + FXGRADTP                           ! Tally the weight of FXGRADTP in XURGENCY
          FXRDCOUNT = FXRDCOUNT + FRDX                               ! Tally the weight of FRDX in XURGENCY
          TOTFXWGT  = TOTFXWGT +&
                     (HVXWEIGT+VVXWEIGT+DOXWEIGT+TPXWEIGT+RDXWEIGT)  ! Tally the Total of Weights
!    Tallying the '% Influence' for Current Timestep Only
          FXHVSPAN  = FXHVSPAN + FXGRADHV                            ! Tally the weight of FXGRADHV in XURGENCY
          FXVVSPAN  = FXVVSPAN + FXGRADVV                            ! Tally the weight of FXGRADVV in XURGENCY
          FXDOSPAN  = FXDOSPAN + FXGRADDO                            ! Tally the weight of FXGRADDO in XURGENCY
          FXTPSPAN  = FXTPSPAN + FXGRADTP                            ! Tally the weight of FXGRADTP in XURGENCY
          FXRDSPAN  = FXRDSPAN + FRDX                                ! Tally the weight of FRDX in XURGENCY
          TOTFXSPAN = TOTFXSPAN +&
                     (HVXWEIGT+VVXWEIGT+DOXWEIGT+TPXWEIGT+RDXWEIGT)  ! Tally the Total Weight

!  Y-direction
          YURGENCY = FRDY

!  Z-direction
          ZURGENCY = (FZGRADHV+FZGRADVV+FZGRADDO+FZGRADTP+FRDZ)/&     ! Weighted Average of Selected Scaled Gradients in Z-direction
                     (HVZWEIGT+VVZWEIGT+DOZWEIGT+TPZWEIGT+RDZWEIGT)
!    Tallying the '% Influence' Over All Timesteps
          FZHVCOUNT = FZHVCOUNT + FZGRADHV                           ! Tally the weight of FZGRADHV in ZURGENCY
          FZVVCOUNT = FZVVCOUNT + FZGRADVV                           ! Tally the weight of FZGRADVV in ZURGENCY
          FZDOCOUNT = FZDOCOUNT + FZGRADDO                           ! Tally the weight of FZGRADDO in ZURGENCY
          FZTPCOUNT = FZTPCOUNT + FZGRADTP                           ! Tally the weight of FZGRADTP in ZURGENCY
          FZRDCOUNT = FZRDCOUNT + FRDZ                               ! Tally the weight of FRDZ in ZURGENCY
          TOTFZWGT  = TOTFZWGT +&
                     (HVZWEIGT+VVZWEIGT+DOZWEIGT+TPZWEIGT+RDZWEIGT)  ! Tally the Total of Weights
!    Tallying the '% Influence' for Current Timestep Only
          FZHVSPAN  = FZHVSPAN + FZGRADHV                            ! Tally the weight of FZGRADHV in ZURGENCY
          FZVVSPAN  = FZVVSPAN + FZGRADVV                            ! Tally the weight of FZGRADVV in ZURGENCY
          FZDOSPAN  = FZDOSPAN + FZGRADDO                            ! Tally the weight of FZGRADDO in ZURGENCY
          FZTPSPAN  = FZTPSPAN + FZGRADTP                            ! Tally the weight of FZGRADTP in ZURGENCY
          FZRDSPAN  = FZRDSPAN + FRDZ                                ! Tally the weight of FRDZ in ZURGENCY
          TOTFZSPAN = TOTFZSPAN +&
                     (HVZWEIGT+VVZWEIGT+DOZWEIGT+TPZWEIGT+RDZWEIGT)  ! Tally the Total Weight

        ELSE                                                       ! URGENCY Calculated as the Maximum of All
!  X-direction                                                     !   Scaled Gradients
          XURGENCY = 0                                               ! Set X-directional URGENCY to 0
          XURGTRAK = 0                                               ! Set URGENCY Tracker to 0 (will be used for tallying purposes)
          IF (ABS(FXGRADHV).GT.ABS(XURGENCY)) THEN                   ! If Scaled Gradient FXGRADHV is Greater Than 0:
            XURGENCY = FXGRADHV                                      !   Set X-directional URGENCY to FXGRADHV, and
            XURGTRAK = 1                                             !   Set URGENCY Tracker to 1 (will be used for tallying purposes)
          END IF
          IF (ABS(FXGRADVV).GT.ABS(XURGENCY)) THEN                   ! If Scaled Gradient FXGRADVV is Greater Than XURGENCY:
            XURGENCY = FXGRADVV                                      !   Set X-directional URGENCY to FXGRADVV, and
            XURGTRAK = 2                                             !   Set URGENCY Tracker to 2 (will be used for tallying purposes)
          END IF
          IF (ABS(FXGRADDO).GT.ABS(XURGENCY)) THEN                   ! If Scaled Gradient FXGRADDO is Greater Than XURGENCY:
            XURGENCY = FXGRADDO                                      !   Set X-directional URGENCY to FXGRADDO, and
            XURGTRAK = 3                                             !   Set URGENCY Tracker to 3 (will be used for tallying purposes)
          END IF
          IF (ABS(FXGRADTP).GT.ABS(XURGENCY)) THEN                   ! If Scaled Gradient FXGRADTP is Greater Than XURGENCY:
            XURGENCY = FXGRADTP                                      !   Set X-directional URGENCY to FXGRADTP, and
            XURGTRAK = 4                                             !   Set URGENCY Tracker to 4 (will be used for tallying purposes)
          END IF
          IF (ABS(FRDX)    .GT.ABS(XURGENCY)) THEN                   ! If Scaled Gradient FRDX is Greater Than XURGENCY:
            XURGENCY = FRDX                                          !   Set X-directional URGENCY to FRDX, and
            XURGTRAK = 5                                             !   Set URGENCY Tracker to 5 (will be used for tallying purposes)
          END IF
!    Tallying the '% Influence' Over All Timesteps
          IF (XURGTRAK.EQ.0) FX00COUNT = FX00COUNT + 1               ! Tally the number of times XURGENCY = 0
          IF (XURGTRAK.EQ.1) FXHVCOUNT = FXHVCOUNT + 1               ! Tally the number of times FXGRADHV is used as XURGENCY
          IF (XURGTRAK.EQ.2) FXVVCOUNT = FXVVCOUNT + 1               ! Tally the number of times FXGRADVV is used as XURGENCY
          IF (XURGTRAK.EQ.3) FXDOCOUNT = FXDOCOUNT + 1               ! Tally the number of times FXGRADDO is used as XURGENCY
          IF (XURGTRAK.EQ.4) FXTPCOUNT = FXTPCOUNT + 1               ! Tally the number of times FXGRADTP is used as XURGENCY
          IF (XURGTRAK.EQ.5) FXRDCOUNT = FXRDCOUNT + 1               ! Tally the number of times FRDX is used as XURGENCY
          TOTFXWGT  = TOTFXWGT + 1                                   ! Tally the number of times decision immediately above is made
!    Tallying the '% Influence' for Current Timestep Only
          IF (XURGTRAK.EQ.1) FXHVSPAN  = FXHVSPAN + 1                ! Tally the weight of FXGRADHV in XURGENCY
          IF (XURGTRAK.EQ.2) FXVVSPAN  = FXVVSPAN + 1                ! Tally the weight of FXGRADVV in XURGENCY
          IF (XURGTRAK.EQ.3) FXDOSPAN  = FXDOSPAN + 1                ! Tally the weight of FXGRADDO in XURGENCY
          IF (XURGTRAK.EQ.4) FXTPSPAN  = FXTPSPAN + 1                ! Tally the weight of FXGRADTP in XURGENCY
          IF (XURGTRAK.EQ.5) FXRDSPAN  = FXRDSPAN + 1                ! Tally the weight of FRDX in XURGENCY
          TOTFXSPAN = TOTFXSPAN + 1                                  ! Tally the Total Weight

!  Y-direction
          YURGENCY = FRDY

!  Z-direction
          ZURGENCY = 0                                               ! Set Z-directional URGENCY to 0
          ZURGTRAK = 0                                               ! Set URGENCY Tracker to 0 (will be used for tallying purposes)
          IF (ABS(FZGRADHV).GT.ABS(ZURGENCY)) THEN                   ! If Scaled Gradient FZGRADHV is Greater Than 0:
            ZURGENCY = FZGRADHV                                      !   Set Z-directional URGENCY to FZGRADHV, and
            ZURGTRAK = 1                                             !   Set URGENCY Tracker to 1 (will be used for tallying purposes)
          END IF
          IF (ABS(FZGRADVV).GT.ABS(ZURGENCY)) THEN                   ! If Scaled Gradient FZGRADVV is Greater Than ZURGENCY:
            ZURGENCY = FZGRADVV                                      !   Set Z-directional URGENCY to FZGRADVV, and
            ZURGTRAK = 2                                             !   Set URGENCY Tracker to 2 (will be used for tallying purposes)
          END IF
          IF (ABS(FZGRADDO).GT.ABS(ZURGENCY)) THEN                   ! If Scaled Gradient FXGRADDO is Greater Than ZURGENCY:
            ZURGENCY = FZGRADDO                                      !   Set Z-directional URGENCY to FXGRADDO, and
            ZURGTRAK = 3                                             !   Set URGENCY Tracker to 3 (will be used for tallying purposes)
          END IF
          IF (ABS(FZGRADTP).GT.ABS(ZURGENCY)) THEN                   ! If Scaled Gradient FXGRADTP is Greater Than ZURGENCY:
            ZURGENCY = FZGRADTP                                      !   Set Z-directional URGENCY to FXGRADTP, and
            ZURGTRAK = 4                                             !   Set URGENCY Tracker to 4 (will be used for tallying purposes)
          END IF
          IF (ABS(FRDZ)    .GT.ABS(ZURGENCY)) THEN                   ! If Scaled Gradient FRDZ is Greater Than ZURGENCY:
            ZURGENCY = FRDZ                                          !   Set Z-directional URGENCY to FRDZ, and
            ZURGTRAK = 5                                             !   Set URGENCY Tracker to 5 (will be used for tallying purposes)
          END IF
!    Tallying the '% Influence' Over All Timesteps
          IF (ZURGTRAK.EQ.0) FZ00COUNT = FZ00COUNT + 1               ! Tally the number of times ZURGENCY = 0
          IF (ZURGTRAK.EQ.1) FZHVCOUNT = FZHVCOUNT + 1               ! Tally the number of times FZGRADHV is used as ZURGENCY
          IF (ZURGTRAK.EQ.2) FZVVCOUNT = FZVVCOUNT + 1               ! Tally the number of times FZGRADVV is used as ZURGENCY
          IF (ZURGTRAK.EQ.3) FZDOCOUNT = FZDOCOUNT + 1               ! Tally the number of times FZGRADDO is used as ZURGENCY
          IF (ZURGTRAK.EQ.4) FZTPCOUNT = FZTPCOUNT + 1               ! Tally the number of times FZGRADTP is used as ZURGENCY
          IF (ZURGTRAK.EQ.5) FZRDCOUNT = FZRDCOUNT + 1               ! Tally the number of times FRDZ is used as ZURGENCY
          TOTFZWGT  = TOTFZWGT + 1                                   ! Tally the number of times decision immediately above is made
!    Tallying the '% Influence' for Current Timestep Only
          IF (ZURGTRAK.EQ.1) FZHVSPAN  = FZHVSPAN + 1                ! Tally the weight of FZGRADHV in ZURGENCY
          IF (ZURGTRAK.EQ.2) FZVVSPAN  = FZVVSPAN + 1                ! Tally the weight of FZGRADVV in ZURGENCY
          IF (ZURGTRAK.EQ.3) FZDOSPAN  = FZDOSPAN + 1                ! Tally the weight of FZGRADDO in ZURGENCY
          IF (ZURGTRAK.EQ.4) FZTPSPAN  = FZTPSPAN + 1                ! Tally the weight of FZGRADTP in ZURGENCY
          IF (ZURGTRAK.EQ.5) FZRDSPAN  = FZRDSPAN + 1                ! Tally the weight of FRDZ in ZURGENCY
          TOTFZSPAN = TOTFZSPAN + 1                                  ! Tally the Total Weight
        END IF


!CHECK Whether the URGENCY Values Violate their Restricted Range of: -1.0 <--> 1.0

        IF (XURGENCY.LT.-1) THEN                        ! Output useful info if XURGENCY .LT. -1.0 so...(-1.0 .LE. XURGENCY .LE. 1.0)
          WRITE(DATADEBUGFN,9460) XURGENCY,FXGRADHV,FXGRADVV,FXGRADDO,FXGRADTP,FRDX
 9460     FORMAT('ERROR: XURGENCY LT -1; XURGENCY=',F9.4,' FXGRADHV='&
                 ,F9.4,' FXGRADVV=',F9.4,' FXGRADDO=',F9.4,' FXGRADTP='&
                 ,F9.4,' FRDX=',F9.4)
!          STOP
        ELSE IF (XURGENCY.GT.1) THEN                    ! Output useful info if XURGENCY .GT. 1.0 so...(-1.0 .LE. XURGENCY .LE. 1.0)
          WRITE(DATADEBUGFN,9470) XURGENCY,FXGRADHV,FXGRADVV,FXGRADDO,FXGRADTP,FRDX
 9470     FORMAT('ERROR: XURGENCY GT 1; XURGENCY=',F9.4,' FXGRADHV='&
                 ,F9.4,' FXGRADVV=',F9.4,' FXGRADDO=',F9.4,' FXGRADTP='&
                 ,F9.4,' FRDX=',F9.4)
!          STOP
        ELSE
        END IF

        IF (ZURGENCY.LT.-1) THEN                        ! Output useful info if ZURGENCY .LT. -1.0 so...(-1.0 .LE. ZURGENCY .LE. 1.0)
          WRITE(DATADEBUGFN,9480) ZURGENCY,FZGRADHV,FZGRADVV,FZGRADDO,FZGRADTP,FRDZ
 9480     FORMAT('ERROR: ZURGENCY LT -1; ZURGENCY=',F9.4,' FZGRADHV='&
                 ,F9.4,' FZGRADVV=',F9.4,' FZGRADDO=',F9.4,' FZGRADTP='&
                 ,F9.4,' FRDZ=',F9.4)
!          STOP
        ELSE IF (ZURGENCY.GT.1) THEN                    ! Output useful info if ZURGENCY .GT. 1.0 so...(-1.0 .LE. ZURGENCY .LE. 1.0)
          WRITE(DATADEBUGFN,9490) ZURGENCY,FZGRADHV,FZGRADVV,FZGRADDO,FZGRADTP,FRDZ
 9490     FORMAT('ERROR: ZURGENCY GT 1; ZURGENCY=',F9.4,' FZGRADHV='&
                 ,F9.4,' FZGRADVV=',F9.4,' FZGRADDO=',F9.4,' FZGRADTP='&
                 ,F9.4,' FRDZ=',F9.4)
!          STOP
        ELSE
        END IF


!Estimate Maximum Fish Cruising Speed - Adjust for Low DO and High Temperature
!  Maximum fish cruising speed rules should be different for day and night
!  Adjust Maximum Fish Cruising Speed for Temperature
        MAXUFISH = FSIZE*MXXSPDL                                        ! MXXSPDL = # of Fish Lengths covered in X-direction
                                                                        !           per second at Overall Maximum Fish Velocity
        MAXUFISH = MAXUFISH*(1-ABS(FISHTEMP(5)-TEMPOPT)/TEMPOPT)        ! PENALTY FUNCTION can be replaced by more elaborate function

        IF (MAXUFISH.LT.0) MAXUFISH = 0                                 ! Prevents horiz cruise spd from going (-) when temperature
                                                                        !   is greater than 30 deg C or less than 0 deg C

        MAXWFISH = FSIZE*MXZSPDL                                        ! MXZSPDL = # of Fish Lengths covered in Z-direction
                                                                        !           per second at Overall Maximum Fish Velocity
        MAXWFISH = MAXWFISH*(1-ABS(FISHTEMP(5)-TEMPOPT)/TEMPOPT)        ! PENALTY FUNCTION can be replaced by more elaborate function

        IF (MAXWFISH.LT.0) MAXWFISH = 0                                 ! Prevents vert cruise spd from going (-) when temperature
                                                                        !   is greater than 30 deg C or less than 0 deg C

!  Adjust Maximum Fish Cruising Speed for Dissolved Oxygen              ! DOTHRES = DO Threshold below which Fish cruising speed
        IF (FISHDO(5).LT.DOTHRES) THEN                                  !           starts to suffer
          MAXUFISH = MAXUFISH*(1-(DOTHRES-FISHDO(5))/DOTHRES)           ! PENALTY FUNCTION can be replaced by more elaborate function
          MAXWFISH = MAXWFISH*(1-(DOTHRES-FISHDO(5))/DOTHRES)           ! PENALTY FUNCTION can be replaced by more elaborate function
          IF (FISHDO(5).LT.0) THEN                                      ! CHECK
            WRITE(DATADEBUGFN,*) 'ERROR: Dissolved Oxygen (FISHDO(5)) LT 0.0'
!            STOP
          ELSE
          END IF
        ELSE
        END IF


!Calculate New Fish Location

        UFISH = MAXUFISH * XURGENCY                                    ! The velocity at which the Fish swims depends on the maximum
        VFISH =    1     * YURGENCY                                    !   speed possible AND the urgency to move in that direction
        WFISH = MAXWFISH * ZURGENCY                                    !   (-1.0 .LE. URGENCY .LE. 1.0)

        OLDFXLOC = FXLOC
        OLDFYLOC = FYLOC
        OLDFZLOC = FZLOC

!        FXLOC = FXLOC + (FXVEL(5) + UFISH) * (JDAY-LRUNDAY)*24*3600    ! New updated Fish X-Location
!        FYLOC = FYLOC + (FYVEL    + VFISH) * (JDAY-LRUNDAY)*24*3600    ! New updated Fish Y-Location
!        FZLOC = FZLOC + (FZVEL(5) + WFISH) * (JDAY-LRUNDAY)*24*3600    ! New updated Fish Z-Location
        FXLOC = FXLOC + (FXVEL(5) + UFISH) * (nfsfreq)*24*3600    ! New updated Fish X-Location
        FYLOC = FYLOC + (FYVEL    + VFISH) * (nfsfreq)*24*3600    ! New updated Fish Y-Location
        FZLOC = FZLOC + (FZVEL(5) + WFISH) * (nfsfreq)*24*3600    ! New updated Fish Z-Location
ELSE
CALL PART_TRANSPORT


ENDIF

!Virtual Sampling only on last iteration of NDT
!  Virtual Gillnet Sampling

        IF (GILLNETSAMPLING.and.n.eq.ndt) THEN
          CALL VGILLNETS
        END IF

!  Virtual Hydroacoustics Sampling

        IF (ACOUSTICSAMPLING.and.n.eq.ndt) THEN
          CALL ACOUSTICS
        END IF


!Schooling and Dispersion Process
IF(.NOT.PARTICLE)THEN
        IF (SCHOOLING) THEN                                       ! Schooling and Dispersion Processes are activated
          DO 144 FSCHL=1,FN                                       ! FN = Fish that have been updated during this timestep so far
            IF (INT(FISHES(FSCHL,6)).EQ.FNBP) THEN                ! FISHES(FSCHL,6) = Branch where other fish are located
              CURRTFISHX = NODES(FKMP,FIMP,1) + FXLOC             ! CURRTFISHX = X Distance (within Branch) where current fish is located
              CURRTFISHZ = NODES(FKMP,FIMP,2) - FZLOC             ! CURRTFISHZ = Z Distance (within Branch) where current fish is located
              OTHFK = INT(FISHES(FSCHL,3))                        ! FISHES(FSCHL,3) = Layer KMP where other fish is located
              OTHFI = INT(FISHES(FSCHL,1))                        ! FISHES(FSCHL,1) = Segment IMP where other fish is located
              OTHERFISHX = NODES(OTHFK,OTHFI,1) + FISHES(FSCHL,2) ! OTHERFISHX = X Distance (within Branch) where other fish is located
              OTHERFISHZ = NODES(OTHFK,OTHFI,2) - FISHES(FSCHL,4) ! OTHERFISHZ = Z Distance (within Branch) where other fish is located
              FNEARX = OTHERFISHX - CURRTFISHX                    ! FNEARX = X-distance from current fish under review to other fish
              FNEARZ = ABS(OTHERFISHZ) - ABS(CURRTFISHZ)          ! FNEARZ = Z-distance from current fish under review to other fish
              XSCHZONEINN = FSIZE * XINZONE                       ! XSCHZONEINN = Edge of Inner Zone of Repulsion (X-direction)
              ZSCHZONEINN = FSIZE * ZINZONE                       ! ZSCHZONEINN = Edge of Inner Zone of Repulsion (Z-direction)
              XSCHZONEMID = FSIZE * XMDZONE                       ! XSCHZONEMID = Edge of Middle Zone of Orientation (X-direction)
              ZSCHZONEMID = FSIZE * ZMDZONE                       ! ZSCHZONEMID = Edge of Middle Zone of Orientation (Z-direction)
              XSCHZONEOUT = FSIZE * XOTZONE                       ! XSCHZONEOUT = Edge of Outer Zone of Attraction (X-direction)
              ZSCHZONEOUT = FSIZE * ZOTZONE                       ! ZSCHZONEOUT = Edge of Outer Zone of Attraction (Z-direction)
              IF ((MILTIME.GT.SKYDAY).AND.&
                  (MILTIME.LE.SKYDUSK)) THEN                      ! During Daytime Hours...
!  This following schooling logic comes from Huth and Wissel (1992) as in "Animal Groups in Three Dimensions" by Parrish & Hamner
                IF ((ABS(FNEARX).LE.XSCHZONEINN).AND.&
                    (ABS(FNEARZ).LE.ZSCHZONEINN)) THEN            ! Other fish is 'too close for comfort'
                  IF (FNEARX.GT.0) FXLOC = FXLOC - FSIZE/4        ! Fish retreats backward by 1/4 fish length (FSIZE/4)
                  IF (FNEARX.LT.0) FXLOC = FXLOC + FSIZE/4        ! Fish retreats forward by 1/4 fish length (FSIZE/4)
                  IF (FNEARZ.GT.0) FZLOC = FZLOC - FSIZE/4        ! Fish retreats upward by 1/4 fish length (FSIZE/4)
                  IF (FNEARZ.LT.0) FZLOC = FZLOC + FSIZE/4        ! Fish retreats downward by 1/4 fish length (FSIZE/4)
                ELSE IF (((ABS(FNEARX).GT.XSCHZONEINN).AND.&
                          (ABS(FNEARX).LE.XSCHZONEMID)).AND.&
                         ((ABS(FNEARZ).GT.ZSCHZONEINN).AND.&
                          (ABS(FNEARZ).LE.ZSCHZONEMID))) THEN     ! Other fish at a comfortable distance apart within the school
                            ! Do Nothing - The fish is comfortable at current location
                ELSE IF (((ABS(FNEARX).GT.XSCHZONEMID).AND.&
                          (ABS(FNEARX).LE.XSCHZONEOUT)).AND.&
                         ((ABS(FNEARZ).GT.ZSCHZONEMID).AND.&
                          (ABS(FNEARZ).LE.ZSCHZONEOUT))) THEN     ! Fish close enough to be attracted to one another
                  IF (FNEARX.GT.0) FXLOC = FXLOC + FNEARX - FSIZE ! Fish closes gap between fish to one fish length (FSIZE)
                  IF (FNEARX.LT.0) FXLOC = FXLOC - FNEARX + FSIZE ! Fish closes gap between fish to one fish length (FSIZE)
                  IF (FNEARZ.GT.0) FZLOC = FZLOC + FNEARZ - FSIZE ! Fish closes gap between fish to one fish length (FSIZE)
                  IF (FNEARZ.LT.0) FZLOC = FZLOC - FNEARZ + FSIZE ! Fish closes gap between fish to one fish length (FSIZE)
                ELSE
                END IF
              ELSE IF ((MILTIME.GT.SKYNIGHT).OR.&
                       (MILTIME.LE.SKYDAWN)) THEN                 ! During Nighttime Hours...
                IF ((ABS(FNEARX).LE.XSCHZONEOUT).AND.&
                    (ABS(FNEARZ).LE.ZSCHZONEOUT)) THEN            ! Fish too close while feeding and will actively disperse
!  This original dispersion logic comes from Thompson et al. (1974) as in "Animal Groups in Three Dimensions" by Parrish & Hamner
                  XSEPCOMFORT = XSCHZONEOUT - ABS(FNEARX)                  ! Separation distance for comfort
                  IF (FNEARX.GT.0) FXLOC = FXLOC -&
                      XSEPCOMFORT*EXP(XDISPCOEFF*XSEPCOMFORT)/2            ! Fish increases gap between fish (the "/2" is arbitrary)
                  IF (FNEARX.LE.0) FXLOC = FXLOC +&
                      XSEPCOMFORT*EXP(XDISPCOEFF*XSEPCOMFORT)/2            ! Fish increases gap between fish (the "/2" is arbitrary)
                  ZSEPCOMFORT = ZSCHZONEOUT - ABS(FNEARZ)                  ! Separation distance for comfort
                  IF (FNEARZ.GT.0) FZLOC = FZLOC -&
                      ZSEPCOMFORT*EXP(ZDISPCOEFF*ZSEPCOMFORT)/2            ! Fish increases gap between fish (the "/2" is arbitrary)
                  IF (FNEARZ.LE.0) FZLOC = FZLOC +&
                      ZSEPCOMFORT*EXP(ZDISPCOEFF*ZSEPCOMFORT)/2            ! Fish increases gap between fish (the "/2" is arbitrary)
!  This original dispersion logic comes from Huth and Wissel (1992) as in "Animal Groups in Three Dimensions" by Parrish & Hamner
!                           IF (FNEARX.GT.0) FXLOC = FXLOC - XSCHZONEOUT   ! Fish increases gap between fish by XSCHZONEOUT
!                           IF (FNEARX.LE.0) FXLOC = FXLOC + XSCHZONEOUT   ! Fish increases gap between fish by XSCHZONEOUT
!                           IF (FNEARZ.GT.0) FZLOC = FZLOC - ZSCHZONEMID   ! Fish increases gap between fish by ZSCHZONEMID
!                           IF (FNEARZ.LE.0) FZLOC = FZLOC + ZSCHZONEMID   ! Fish increases gap between fish by ZSCHZONEMID
                ELSE
                END IF
              ELSE                                                ! During Morning(Dawn) and Evening(Dusk)...
                ! Do Nothing, Fish will move individually without schooling or dispersion
              END IF
            ELSE                                                  ! Other fish not located in same branch as current fish
            END IF
  144     CONTINUE
        ELSE
        END IF
ENDIF

!***************************************************************
!  C H E C K   F O R   B O U N D A R Y   V I O L A T I O N S   *
!***************************************************************

! Debug check SW 2/01/02
!FIMPTMP=FIMP
!FKMPTMP=FKMP
!FXLOCTMP=FXLOC
!FZLOCTMP=FZLOC
!FYLOCTMP=FYLOC
!if(abs(fzloc).gt.2*h(fkmp,fjr).or.abs(fxloc).gt.2*dlx(fimp))then
!write(DIAGFN,*)'Jday:',jday,' FZLOC>2*DZ or FXLOC>2DLX'
!write(DIAGFN,*)'fn,fxloc,fyloc,fzloc,fimp,fkmp'
!write(DIAGFN,*)fn,fxloc,fyloc,fzloc,fimp,fkmp
!write(DIAGFN,*)'fxvel(5),fzvel(5)',fxvel(5),fzvel(5)
!end if
!Check for Boundary Violations: Lateral Direction Check (1 of 2)

  if(fyvel.eq.0.0)then   ! SW 2/01/01 logic for removing particles laterally also
        IF (FYLOC.LT.0) THEN    ! reflect particles
           FYLOC = B(FKMP,FIMP)*YREFL 
        ELSE IF (FYLOC.GT.B(FKMP,FIMP)) THEN
           FYLOC = B(FKMP,FIMP)*(1-YREFL) 
        ELSE
        END IF
  else
        IF (FYLOC.LT.0.and.fyvel.lt.0.0) THEN  ! Particle is reflected since only remove from RHS Inflow has pushed it to one side of bank
           FYLOC = B(FKMP,FIMP)*YREFL 
        ELSE IF (FYLOC.GT.B(FKMP,FIMP).and.fyvel.gt.0.0) THEN   
           ISWITCH = 1
           ! SW 2/01/01 Track timing of fish movement from system
           fishes(fn,14)=JDAY               ! Time fish left system
           fishes(fn,15)=JDAY-DELAYDATE     ! Detention time of fish in system
           fishes(fn,16)=1.0   ! Signifies the lateral removal
           ! End Section SW 2/01/01
           GOTO 20
           ! Particle is removed
        ELSEIF (FYLOC.GT.B(FKMP,FIMP).and.fyvel.lt.0.0) THEN   ! reflect - rarely occurs if ever
           FYLOC = B(FKMP,FIMP)*(1-YREFL) 
        ELSEIF (FYLOC.LT.0.and.fyvel.gt.0.0) THEN  ! reflect - rarely occurs if ever
           FYLOC = B(FKMP,FIMP)*YREFL 
        END IF
  end if

!Check for Boundary Violations: Horizontal Direction

   10   IF (FXLOC.GT.DLX(FIMP)) THEN
          IF (DEBUG) WRITE(DATADEBUGFN,*) 'FXLOC GT DLX'
          TAG = 0
          IF (LIMPBR(FIMP+1,1).GT.0) TAG = 1                            ! Branch joining left bank at new segment
          IF (RIMPBR(FIMP+1,1).GT.0) TAG = 2                            ! Branch joining right bank at new segment
          IF ((RIMPBR(FIMP+1,1).GT.0).AND.(LIMPBR(FIMP+1,1).GT.0))TAG = 3  ! Branches joining at both banks
          IF (PREVENTBRCHSWITCH) TAG = 0                         ! Prevents fish from moving upstream into another branch
          IF (FIMP.EQ.DS(DNBP)) THEN  
            if(fxvel(5).gt.0.0)then   ! fish leaves the system if there is a structure outflow
              IF (DEBUG) WRITE(DATADEBUGFN,*) 'FIMP+1 GT DS(DNBP)'
              WRITE(DATADEBUGFN,9110) JDAY,FNBP
9110          FORMAT(/,'Fish left the system on day ',F9.3,'  thru Branch ',I6,' at the downstream end')
              ISWITCH = 1
              ! SW 2/01/01 Track timing of fish movement from system
              fishes(fn,14)=JDAY               ! Time fish left system
              fishes(fn,15)=JDAY-DELAYDATE     ! Detention time of fish in system
              fishes(fn,16)=2.     ! Signifies longitudinal removal
              ! End Section SW 2/01/01
              GOTO 20
            else   ! fish is "reflected"  ! SW 2/01/01
              FXLOC = DLX(FIMP)*(1-XREFL)
            end if
          !!*** This code only allows fish movement downstream when there is pumpback operation - need to adjust for new V3 model ! SW 1/14/01 ***
          !ELSE IF (FNBP.gt.0.and.FIMP.EQ.DS(FNBP)) THEN           ! SW Check for JBG       ! Logic for moving fish to downstream waterbody
          !    IF (DEBUG) WRITE(DATADEBUGFN,*) 'FIMP+1 GT DS(FNBP)'
          !    FXLOC = FXLOC - DLX(FIMP)
          !    FIMP = US(JBP)
          !    FNBP = JBP
          !    CALL WHATJR
          !    KTWBF = KTWB(FJR)                                      ! Water Surface Layer of waterbody fish is located in
          ELSE IF ((FIMP+1).GT.DS(FNBP)) THEN                    ! Logic for moving fish to downstream branch
              IF (DEBUG) WRITE(DATADEBUGFN,*) 'FIMP+1 GT DS(FNBP)'        !   within a given waterbody,JR.
              IF (DHS(FNBP).GT.0) THEN
                FYLOCTMP = FYLOC
                 IF (LIMPBR(DHS(FNBP),2).EQ.FNBP) THEN              ! Moving downstream to left bank of new segmt
                   FYLOC = FXLOC - DLX(FIMP)                        ! Lateral direction becomes the longitudinal
                   FXLOC = (1-FYLOCTMP/B(FKMP,FIMP))*DLX(DHS(FNBP)) !   direction in the new branch
                 ELSE IF (RIMPBR(DHS(FNBP),2).EQ.FNBP) THEN         ! Moving downstream to right bank of new segmt
                   FYLOC = B(FKMP,DHS(FNBP)) - (FXLOC - DLX(FIMP))  ! Lateral direction becomes the longitudinal
                   FXLOC = FYLOCTMP/B(FKMP,FIMP) * DLX(DHS(FNBP))   !   direction in the new branch
                 ELSE
                   WRITE(FINALFN,*) 'ERROR: Unable to determine whether moving downstream from right or left bank to new branch'
                   write(FINALFN,*)'JDAY=',jday,' FN,FIMP,FKMP:',fn,fimp,fkmp
                   do jf=1,3     !fn
                   write(FINALFN,'(i7,1x,<fpara>(f10.3,1x))')jf,(fishes(jf,i),i=1,fpara)
                   end do
                   STOP
                 END IF
                FIMP = DHS(FNBP)
                IF (DEBUG) WRITE(DATADEBUGFN,*) 'CALLING FINDNEWBR'
                CALL FINDNEWBR
                IF (DEBUG) WRITE(DATADEBUGFN,*) 'FINISHED FINDNEWBR CALL'
!            ELSE   ! SW not a correct error message
!              WRITE(*,*) 'ERROR: Inappropriate DHS(FNBP)'
!              WRITE(*,9160) FNBP, DHS(FNBP)
! 9160         FORMAT(' Branch = ',I6,'   DHS(FNBP) = ',I6)
!              STOP 
             END IF
          ELSE IF ((TAG.EQ.1).OR.(TAG.EQ.2)) THEN
            IF (DEBUG) WRITE(DATADEBUGFN,*) 'CALLING TAG124578:FXLOC GT DLX'
            CALL TAG124578
            IF (DEBUG) WRITE(DATADEBUGFN,*) 'FINISHED TAG124578 CALL'
          ELSE IF (TAG.EQ.3) THEN
            IF (DEBUG) WRITE(DATADEBUGFN,*) 'CALLING TAG369:FXLOC GT DLX'
            CALL TAG369
            IF (DEBUG) WRITE(DATADEBUGFN,*) 'FINISHED TAG369 CALL'
          ELSE
            IF (DEBUG) WRITE(DATADEBUGFN,*) 'ORDINARY FXLOC GT DLX'
            FXLOC = FXLOC - DLX(FIMP)
            FIMP = FIMP + 1                                      ! Fish advances one segment downstream
          END IF
          IF (FXLOC.GT.DLX(FIMP).or.FXLOC.LT.0.)then
             write(FINALFN,*)'Possible Error - looping GOTO10 - Set DNBP to Proper Branch'
             write(FINALFN,*)'JDAY=',jday,' FN,FIMP,FKMP:',fn,fimp,fkmp
             do jf=1,3         !fn
             write(FINALFN,'(i7,1x,<fpara>(f10.3,1x))')jf,(fishes(jf,i),i=1,fpara)
             end do 
             GOTO 10                        ! More thought needs to go into this overshooting logic
          ENDIF
        ELSE IF (FXLOC.LT.0) THEN
          IF (DEBUG) WRITE(DATADEBUGFN,*) 'FXLOC LT 0'
          TAG = 0
!          IF ((FIMP.EQ.US(JBP)).OR.(FIMP.EQ.US(UNBP))) THEN      ! The UpStream end of Branch JBP or UNBP
          IF ((FIMP.EQ.US(UNBP))) THEN      ! The UpStream end of Branch JBP or UNBP ! SW 2/01/01 removed JBP condition
             IF (DEBUG) WRITE(DATADEBUGFN,*) 'FIMP EQ US(JBP or UNBP)'
!            WRITE(DATADEBUGFN,9150) JDAY,FNBP,FKMP
! 9150       FORMAT(/,'The fish left the system on day ',F9.3,    ! UpStream end of Branch JBP or UNBP where fish are
!     .               '  thru FNBP ',I6,' at the upstream end',   !   collected and tallied when moving upstream
!     .             /,' via Layer FKMP ',I6)
!            ISWITCH = 1
!            GOTO 20
             FXLOC = DLX(FIMP)*XREFL   
          ELSE IF (FIMP.EQ.CUS(FNBP)) THEN
             IF (DEBUG) WRITE(DIAGFN,*) 'FIMP EQ CUS(FNBP)'
!             if(fxloc.lt.(dlx(fimp)*xrefl))FXLOC = DLX(FIMP)*XREFL   
             FXLOC = DLX(FIMP)*XREFL  
          ELSE
             IF (LIMPBR(FIMP-1,1).GT.0) TAG = 4
             IF (RIMPBR(FIMP-1,1).GT.0) TAG = 5
             IF ((RIMPBR(FIMP-1,1).GT.0).AND.(LIMPBR(FIMP-1,1).GT.0)) THEN 
             TAG = 6
          ELSE
          END IF
            IF (PREVENTBRCHSWITCH) TAG = 0                       ! Prevents fish from moving upstream into another branch
          END IF
          IF ((TAG.EQ.4).OR.(TAG.EQ.5)) THEN
            IF (DEBUG) WRITE(DATADEBUGFN,*) 'CALLING TAG124578:FXLOC LT 0'
            CALL TAG124578
            IF (DEBUG) WRITE(DATADEBUGFN,*) 'FINISHED TAG124578 CALL'
!            write(DIAGFN,*)'Calling TAG124578 A'  ! Debug SW
          ELSE IF (TAG.EQ.6) THEN
            IF (DEBUG) WRITE(DATADEBUGFN,*) 'CALLING TAG369:FXLOC LT 0'
!            write(DIAGFN,*)'Calling TAG369 A'  ! Debug SW
            CALL TAG369
            IF (DEBUG) WRITE(DATADEBUGFN,*) 'FINISHED TAG369 CALL'
          ELSE
            IF (DEBUG) WRITE(DATADEBUGFN,*) 'ORDINARY FXLOC LT 0'
            FIMP = FIMP - 1                                    ! There are still segments upstream for the fish to
            FXLOC = DLX(FIMP) + FXLOC                          !    move to
          END IF
! These are problems that need correcting !!!     SW 2/1/01
          IF (FXLOC.GT.DLX(FIMP)) GOTO 10 
          IF (FXLOC.LT.0) GOTO 10
!        ELSEIF(FXLOC.EQ.0.0.and.cus(fnbp).eq.fimp.and.fxvel(5).eq.0.)then
!             FXLOC = DLX(FIMP)*XREFL    ! REflect partilce off boundary 
        ELSE
          TAG = 0
          IF (LIMPBR(FIMP,1).GT.0) TAG = 7                          ! Branch joining left bank at current segment
          IF (RIMPBR(FIMP,1).GT.0) TAG = 8                          ! Branch joining right bank at current segment
          IF ((RIMPBR(FIMP,1).GT.0).AND.(LIMPBR(FIMP,1).GT.0)) THEN ! Branches joining at both banks
            TAG = 9
          ELSE
          END IF
          IF (PREVENTBRCHSWITCH) TAG = 0                       ! Prevents fish from moving upstream into another branch
          IF ((TAG.EQ.7).OR.(TAG.EQ.8)) THEN
            IF (DEBUG) WRITE(DATADEBUGFN,*) 'CALLING TAG124578:&
                                      AT INTERSECTING FIMP'
            CALL TAG124578
            IF (DEBUG) WRITE(DATADEBUGFN,*) 'FINISHED TAG124578 CALL'
!             write(DIAGFN,*)'Calling TAG124578 B'  ! Debug SW

          ELSE IF (TAG.EQ.9) THEN
            IF (DEBUG) WRITE(DATADEBUGFN,*) 'CALLING TAG369:&
                                      AT INTERSECTING FIMP'
            CALL TAG369
!            write(DIAGFN,*)'Calling TAG369 B'  ! Debug SW
            IF (DEBUG) WRITE(DATADEBUGFN,*) 'FINISHED TAG369 CALL'
          ELSE
            IF (DEBUG) WRITE(DATADEBUGFN,*) ' NO SPECIAL RULES IMPLEMENTED'
          END IF
        END IF
        IF (DEBUG) WRITE(DATADEBUGFN,*) 'HORIZONTAL PLANE CALC OK'


! Debug

!IF(FIMP.eq.cus(fnbp).and.fxloc.lt.10..and.jday.gt.1.0)then
!write(DIAGFN,*)'Jday:',jday,' FXLOC<10 and FIMP=CUS'
!write(DIAGFN,*)'fn,fxloc,fyloc,fzloc,fimp,fkmp'
!write(DIAGFN,*)fn,fxloc,fyloc,fzloc,fimp,fkmp
!write(DIAGFN,*)'fxvel(5),fzvel(5),u(fk,fi):',fxvel(5),fzvel(5),u(fkmp,fimp)
!write(DIAGFN,*)'CUS(FNBP),KT(US),u(fkmp+1,fimp),u(fkmp+2,fimp):',cus(fnbp),KTWBF,u(fkmp+1,fimp),u(fkmp+2,fimp)
!write(DIAGFN,*)'u(fkmp,fimp+1),u(fkmp+1,fimp+1):',u(fkmp,fimp+1),u(fkmp+1,fimp+1)
!do i=cus(fnbp),cus(fnbp)+3
!do k=KTWBF,kb(i)
!write(DIAGFN,*)k,i,u(k,i)
!end do
!end do
!end if


!Check for Boundary Violations: Vertical Direction Check (1 of 2)

        KTWBF = KTWB(FJR)                            ! Water Surface Layer of waterbody fish is located in

   70   IF ((FKMP.Le.KTI(FIMP)).OR.&   ! SW 2/01/01 Change so that do a check if in "air" or not
            (FKMP.GT.KB(FIMP))) THEN
          IF (DEBUG) WRITE(DATADEBUGFN,*) 'B(FKMP,FIMP) EQ 0 CHECK&
                     ACTIVATED PRECEEDING FKMP UPDATE'
          IF (FKMP.GT.KB(FIMP)) THEN               !  leaving fish hanging in air?   --  NO
            IF (DEBUG) WRITE(DATADEBUGFN,*) 'FKMP GT KB'
            FKMP = FKMP - 1
            FZLOC = H(FKMP,FJR)*(1-ZBOTREFL)           !    Logic for rebounding off of the bottom
            GOTO 70
!          ELSEif(fkmp.lt.kti(fimp))then                                     !  leaving fish hanging in air?   --  YES
!            IF (DEBUG) WRITE(DATADEBUGFN,*) 'FKMP LE KTWBF'
!            FKMP = KTI(FIMP)                       ! May want to modify this since now any fish in KTI or above
!            FZLOCTEMP = H(KTWBF,FJR)-Z(FIMP)            !    will be placed at the water surface - see ELSE Statemnt
!            write(DIAGFN,*)'Layer elevation check A SW, jday,fn',jday,fn   ! SW Debug
!            IF (KTI(FIMP).LT.KTWBF) THEN
!              DO 60 KK = KTWBF,(FKMP+1),-1
!                FZLOCTEMP = FZLOCTEMP - H(KK,FJR)
!   60         CONTINUE
!              FZLOC = H(KTI(FIMP),FJR) - FZLOCTEMP
!            ELSE IF (KTI(FIMP).EQ.KTWBF) THEN
!              FZLOC = -Z(FIMP)
!            ELSE
!              WRITE(*,*) 'ERROR: Moving fish back into water from air'
!            END IF
          elseif(fkmp.le.kti(fimp))then
          fkmp=kti(fimp)
! check and make sure the fish/particles are below the water surface   ! SW 2/01/01
               if(kti(fimp).eq.KTWBF)then
                  if(fzloc.lt.z(fimp))then
                   fzloc=z(fimp)+h(fkmp,fjr)*(1-zsurrefl)
!                   write(DIAGFN,*)'Layer elevation check B SW, jday,fn',jday,fn   ! SW Debug
                  end if
               else
                   if(fzloc.le.h(kti(fimp),fjr)+z(fimp))then
                     if(z(fimp).gt.0.0)write(*,*)'Debug Error in vertical reflection:FN,I,K',fn,fimp,fkmp
                   fzloc=(h(kti(fimp),fjr)+z(fimp))+h(fkmp,fjr)*(1-zsurrefl)
!                   write(DIAGFN,*)'Layer elevation check C SW, jday,fn',jday,fn   ! SW Debug
                   end if
               end if
          END IF
        ELSE                                           ! Moving downstream result in fish below the bottom? -- NO
        END IF

!Check for Boundary Violations: Vertical Direction Check (2 of 2)

   30   IF (FZLOC.GT.H(FKMP,FJR)) THEN                     ! Did fish move down below current layer?  -- YES
          IF (DEBUG) WRITE(DATADEBUGFN,*) 'FZLOC GT H'
          IF (FKMP+1.GT.KB(FIMP)) THEN                 ! Did fish move below bottom layer? -- YES
            FZLOC = H(FKMP,FJR)*(1-ZBOTREFL)               !    Logic for rebounding off of the bottom
          ELSE                                         ! Did fish move below bottom layer? -- NO
            FZLOC = FZLOC - H(FKMP,FJR)
            FKMP = FKMP + 1
            IF (FZLOC.GT.H(FKMP,FJR)) GOTO 30              ! Did fish move down more than one layer - must improve this.
          END IF
        ELSE IF (FZLOC.LT.0) THEN                      ! Did fish move up?  --  YES
          IF (DEBUG) WRITE(DATADEBUGFN,*) 'FZLOC LT 0'
          IF (FKMP.EQ.(KTWBF+1)) THEN                   ! Is fish just below surface layer KTWBF?  --  YES
            IF (DEBUG) WRITE(DATADEBUGFN,*) 'FKMP EQ KTWBF+1'
            IF (ABS(FZLOC).GT.(H(KTWBF,FJR)-Z(FIMP))) THEN  ! Did fish move above water surface?  --  YES
              FROMKTBOT = (H(KTWBF,FJR)-Z(FIMP))*(1-ZSURREFL)
            ELSE                                       ! Did fish move above water surface?  --  NO
              FROMKTBOT = ABS(FZLOC)
            END IF
            SURFCALC = 1
          ELSE IF (FKMP.EQ.KTWBF) THEN
            IF (DEBUG) WRITE(DATADEBUGFN,*) 'FKMP EQ KTWBF'
            IF ((ABS(FZLOC)+H(KTWBF,FJR)).GT.(H(KTWBF,FJR)-Z(FIMP))) THEN
              FROMKTBOT = (H(KTWBF,FJR)-Z(FIMP))*(1-ZSURREFL)
            ELSE
              FROMKTBOT = ABS(FZLOC) + H(KTWBF,FJR)
            END IF
            SURFCALC = 1
          ELSE IF (FKMP.LE.(KTWBF-1)) THEN
            IF (DEBUG) WRITE(DATADEBUGFN,*) 'FKMP LE KTWBF-1'
            FROMKTBOT = 0.0
            DO 40 KK = KTWBF,(FKMP+1),-1
              FROMKTBOT = FROMKTBOT + H(KK,FJR)
   40       CONTINUE
            FROMKTBOT = FROMKTBOT + ABS(FZLOC)
            IF (FROMKTBOT.GT.(H(KTWBF,FJR)-Z(FIMP))) THEN
              FROMKTBOT = (H(KTWBF,FJR)-Z(FIMP))*(1-ZSURREFL)
            ELSE
              FROMKTBOT = FROMKTBOT
            END IF
            SURFCALC = 1
          ELSE
            IF (DEBUG) WRITE(DATADEBUGFN,*) 'ORDINARY FZLOC LT 0'
            FKMP = FKMP - 1
            FZLOC = H(FKMP,FJR) + FZLOC
            IF (FZLOC.LT.0) GOTO 30                    ! Did fish move up more than one layer - must improve this.
            SURFCALC = 0 
          END IF
          IF (SURFCALC.EQ.1) THEN
            FKMPTEMP = KTWBF
   50       IF (FROMKTBOT.GT.H(FKMPTEMP,FJR)) THEN
              FROMKTBOT = FROMKTBOT - H(FKMPTEMP,FJR)
              FKMPTEMP = FKMPTEMP - 1
              GOTO 50
            ELSE
              FZLOC = H(FKMPTEMP,FJR) - FROMKTBOT
              FKMP = FKMPTEMP
            END IF
          ELSE
          END IF
        ELSE
        END IF
        IF (FKMP.LT.KTI(FIMP)) THEN
          WRITE(DATADEBUGFN,9350) FKMP,FZLOC,KTI(FIMP),(H(KTWBF,FJR)-Z(FIMP))
 9350     FORMAT('FKMP.LT.KTI(FIMP) ==> FKMP=',I6,' FZLOC=',F10.3&
                 ,' KTI(FIMP)=',I6,' and (H(KTWBF)-Z(FIMP))=',F10.3)
        ELSE
        END IF
        IF (DEBUG) WRITE(DATADEBUGFN,*) 'VERTICAL PLANE CALC OK'

!Check for Boundary Violations: Lateral Direction Check (2 of 2)

        IF (FYLOC.LT.0) THEN
          FYLOC = B(FKMP,FIMP)*YREFL 
        ELSE IF (FYLOC.GT.B(FKMP,FIMP)) THEN
          FYLOC = B(FKMP,FIMP)*(1-YREFL) 
        ELSE
        END IF

      END IF

   20 CONTINUE

! Debug output  SW 2/01/01

!if(abs(fimp-fimptmp).gt.1.or.abs(fkmp-fkmptmp).gt.1)then
!write(DIAGFN,*)'jday:',jday
!write(DIAGFN,*)'fn,fxloc,fyloc,fzloc,fkmp,fimp'
!write(DIAGFN,*)fn,fxloctmp,fyloctmp,fzloctmp,fkmptmp,fimptmp
!write(DIAGFN,*)'Updated values:fn,fxloc,fyloc,fzloc,fkmp,fimp'
!write(DIAGFN,*)fn,fxloc,fyloc,fzloc,fkmp,fimp
!endif

! Collector Check  SW 2/16/01
          if(collector.eq.' ON')then
            do i=1,ncollector
            if(fimp.eq.icoll(i).and.fkmp.ge.icollt(i).and.fkmp.le.icollb(i))then
              ! SW 2/01/01 Track timing of fish movement from system
              fishes(fn,14)=JDAY               ! Time fish left system
              fishes(fn,15)=JDAY-DELAYDATE     ! Detention time of fish in system
              fishes(fn,16)=3.     ! Signifies collector removal=3
              ! End Section SW 2/01/01
              exit
            endif
            end do
          end if
! End Collector Check

      call findnewbr
      FISHES(FN,1)  = REAL(FIMP)      ! Segment IMP where fish is located
      FISHES(FN,2)  =      FXLOC      ! Location of fish within segment IMP from upstream side
      FISHES(FN,3)  = REAL(FKMP)      ! Layer KMP where fish is located
      FISHES(FN,4)  =      FZLOC      ! Location of fish within layer KMP from top side
      FISHES(FN,5)  =      FYLOC      ! Location of fish laterally from left bank (in plan view)
      FISHES(FN,6)  = REAL(FNBP)      ! Branch where fish is located
      FISHES(FN,7)  =      FSIZE      ! Size (i.e., length) of fish in meters
      FISHES(FN,8)  =      FAGE+nfsfreq    ! Age of the fish  ! SW 2/01/01
      FISHES(FN,9)  =      UFISH      ! Longitudinal velocity of fish relative to water
      FISHES(FN,10) =      VFISH      ! Lateral velocity of fish relative to water
      FISHES(FN,11) =      WFISH      ! Vertical velocity of fish relative to water
      FISHES(FN,12) = REAL(ISWITCH)   ! Is fish still in system?: Yes=0  No=1
      FISHES(FN,13) = REAL(FSNAG)     ! Gillnet # fish is snagged in; = 0 if fish not in a gillnet

  131 CONTINUE
END DO   ! END OF NDT LOOP

   28 CONTINUE

      LRUNDAY = JDAY


!*******************************************************
!  O U T P U T   A N D   E N D   S U B R O U T I N E   *
!*******************************************************

      IF ((FCOUNT.EQ.1).OR.(NIT.EQ.0)) THEN
        CALL FISHPLOT                                                   ! This plots the location of virtual fish

        IF (JDAY .LT. DELAYDATE) THEN                                   ! This prevents division by zero when fish
          TOTFXSPAN = 1                                                 !   are not moving.
          TOTFZSPAN = 1
        END IF

!Output Information to the Bar Chart File (X-Direction)
        IF(TOTFXSPAN /= 0)THEN
        FXREACTVEL = ((FXHVSPAN+FXVVSPAN)/TOTFXSPAN)*100
        FXREACTTMP = (FXTPSPAN/TOTFXSPAN)*100
        FXREACTDO  = (FXDOSPAN/TOTFXSPAN)*100
        FXREACTRND = (FXRDSPAN/TOTFXSPAN)*100
        ELSE
        FXREACTVEL = 0
        FXREACTTMP = 0
        FXREACTDO  = 0
        FXREACTRND = 0
        ENDIF

        IF (NIT.EQ.0) THEN
          WRITE(BARCHRTXFN,9552) 
 9552     FORMAT('TITLE = "Bar Chart: % of X-Directional Movement due to&
       Each Influence Factor"')
          WRITE(BARCHRTXFN,9553) 
 9553     FORMAT('VARIABLES = "% X-dir Movement", "Infl Factor"')
        ELSE
          WRITE(BARCHRTXFN,*) ' '
        END IF
        WRITE(BARCHRTXFN,9529) JDAY
 9529   FORMAT('ZONE T="JDAY ',F9.2,'", I=4, F=POINT')
        WRITE(BARCHRTXFN,9530) FXREACTVEL, 4
        WRITE(BARCHRTXFN,9530) FXREACTTMP, 3
        WRITE(BARCHRTXFN,9530) FXREACTDO, 2
        WRITE(BARCHRTXFN,9530) FXREACTRND, 1
 9530   FORMAT(F5.1,'   ',I2)
        IF (NIT.EQ.0) THEN
          WRITE(BARCHRTXFN,9542) 
 9542     FORMAT('TEXT X=51.5, Y=5.1, F=HELV-BOLD, CS=GRID, HU=FRAME,&
      AN=MIDCENTER, C=BLACK, H=11, T="X-Movement"')
          WRITE(BARCHRTXFN,9543) 
 9543     FORMAT('TEXT X=-5, Y=4, F=HELV-BOLD, CS=GRID, HU=FRAME,&
      AN=MIDRIGHT, C=CUSTOM7, H=10, T="VEL"')
          WRITE(BARCHRTXFN,9544) 
 9544     FORMAT('TEXT X=-5, Y=3, F=HELV-BOLD, CS=GRID, HU=FRAME,&
      AN=MIDRIGHT, C=CUSTOM7, H=10, T="TMP"')
          WRITE(BARCHRTXFN,9545) 
 9545     FORMAT('TEXT X=-5, Y=2, F=HELV-BOLD, CS=GRID, HU=FRAME,&
      AN=MIDRIGHT, C=CUSTOM7, H=10, T="DO"')
          WRITE(BARCHRTXFN,9546) 
 9546     FORMAT('TEXT X=-5, Y=1, F=HELV-BOLD, CS=GRID, HU=FRAME,&
      AN=MIDRIGHT, C=CUSTOM7, H=10, T="RND"')
        END IF
        WRITE(BARCHRTXFN,9531) FXREACTVEL,NZONES
 9531   FORMAT('TEXT X=105, Y=4, F=HELV-BOLD, CS=GRID, HU=FRAME,&
      AN=MIDLEFT, C=CUSTOM7, H=10, T="',F5.1,'", ZN=',I6)
        WRITE(BARCHRTXFN,9532) FXREACTTMP,NZONES
 9532   FORMAT('TEXT X=105, Y=3, F=HELV-BOLD, CS=GRID, HU=FRAME,&
      AN=MIDLEFT, C=CUSTOM7, H=10, T="',F5.1,'", ZN=',I6)
        WRITE(BARCHRTXFN,9533) FXREACTDO,NZONES
 9533   FORMAT('TEXT X=105, Y=2, F=HELV-BOLD, CS=GRID, HU=FRAME,&
      AN=MIDLEFT, C=CUSTOM7, H=10, T="',F5.1,'", ZN=',I6)
        WRITE(BARCHRTXFN,9534) FXREACTRND,NZONES
 9534   FORMAT('TEXT X=105, Y=1, F=HELV-BOLD, CS=GRID, HU=FRAME,&
      AN=MIDLEFT, C=CUSTOM7, H=10, T="',F5.1,'", ZN=',I6)


!Output Information to the Bar Chart File (Z-Direction)
        IF(TOTFZSPAN /= 0)THEN    ! SW 5/1/15
        FZREACTVEL = ((FZHVSPAN+FZVVSPAN)/TOTFZSPAN)*100
        FZREACTTMP = (FZTPSPAN/TOTFZSPAN)*100
        FZREACTDO  = (FZDOSPAN/TOTFZSPAN)*100
        FZREACTRND = (FZRDSPAN/TOTFZSPAN)*100
        ELSE
        FZREACTVEL = 0
        FZREACTTMP = 0
        FZREACTDO  = 0
        FZREACTRND = 0
            
        ENDIF
        

        IF (NIT.EQ.0) THEN
          WRITE(BARCHRTZFN,9554) 
 9554     FORMAT('TITLE = "Bar Chart: % of Z-Directional Movement due to&
       Each Influence Factor"')
          WRITE(BARCHRTZFN,9555) 
 9555     FORMAT('VARIABLES = "Infl Factor", "% Z-dir Movement"')
        ELSE
          WRITE(BARCHRTZFN,*) ' '
        END IF
        WRITE(BARCHRTZFN,9529) JDAY
        WRITE(BARCHRTZFN,9535) 1, FZREACTVEL
        WRITE(BARCHRTZFN,9535) 2, FZREACTTMP
        WRITE(BARCHRTZFN,9535) 3, FZREACTDO
        WRITE(BARCHRTZFN,9535) 4, FZREACTRND
 9535   FORMAT(I2,'   ',F5.1)
        IF (NIT.EQ.0) THEN
          WRITE(BARCHRTZFN,9547) 
 9547     FORMAT('TEXT X=2.5, Y=130, F=HELV-BOLD, CS=GRID, HU=FRAME,&
      AN=CENTER, C=BLACK, H=11, T="Z-Movement"')
          WRITE(BARCHRTZFN,9548) 
 9548     FORMAT('TEXT X=1, Y=-5, F=HELV-BOLD, CS=GRID, HU=FRAME,&
      AN=HEADCENTER, C=CUSTOM7, H=10, T="VEL"')
          WRITE(BARCHRTZFN,9549) 
 9549     FORMAT('TEXT X=2, Y=-5, F=HELV-BOLD, CS=GRID, HU=FRAME,&
      AN=HEADCENTER, C=CUSTOM7, H=10, T="TMP"')
          WRITE(BARCHRTZFN,9550) 
 9550     FORMAT('TEXT X=3, Y=-5, F=HELV-BOLD, CS=GRID, HU=FRAME,&
      AN=HEADCENTER, C=CUSTOM7, H=10, T="DO"')
          WRITE(BARCHRTZFN,9551) 
 9551     FORMAT('TEXT X=4, Y=-5, F=HELV-BOLD, CS=GRID, HU=FRAME,&
      AN=HEADCENTER, C=CUSTOM7, H=10, T="RND"')
        END IF
        WRITE(BARCHRTZFN,9536) FZREACTVEL,NZONES
 9536   FORMAT('TEXT X=1, Y=105, F=HELV-BOLD, CS=GRID, HU=FRAME,&
      AN=CENTER, C=CUSTOM7, H=10, T="',F5.1,'", ZN=',I6)
        WRITE(BARCHRTZFN,9537) FZREACTTMP,NZONES
 9537   FORMAT('TEXT X=2, Y=105, F=HELV-BOLD, CS=GRID, HU=FRAME,&
      AN=CENTER, C=CUSTOM7, H=10, T="',F5.1,'", ZN=',I6)
        WRITE(BARCHRTZFN,9538) FZREACTDO,NZONES
 9538   FORMAT('TEXT X=3, Y=105, F=HELV-BOLD, CS=GRID, HU=FRAME,&
      AN=CENTER, C=CUSTOM7, H=10, T="',F5.1,'", ZN=',I6)
        WRITE(BARCHRTZFN,9539) FZREACTRND,NZONES
 9539   FORMAT('TEXT X=4, Y=105, F=HELV-BOLD, CS=GRID, HU=FRAME,&
      AN=CENTER, C=CUSTOM7, H=10, T="',F5.1,'", ZN=',I6)
      END IF

      IF (DEBUG) WRITE(DATADEBUGFN,*) 'SUBROUTINE FISH COMPLETED'
   

   97 CONTINUE

      return
      END

!***********************************************************************
!***********************************************************************
!*                                                                    **
!*           R A N D O M   N U M B E R   G E N E R A T O R            **
!*                                                                    **
!***********************************************************************
!***********************************************************************
!  This Subroutine Calls No Other Subroutines


      SUBROUTINE RANDOM(SEED,RNDX)             ! Random Number Generator Subroutine
      IMPLICIT NONE                            !   obtained from "STRUCTURED FORTRAN
                                               !   77 for Engineers and Scientists, 5th
      INTEGER   SEED                           !   Edition - Author: Delores M. Etter"
      REAL      RNDX                           !   The random # is between 0.0 and 1.0
      

      SEED = 2045*SEED + 1
      SEED = SEED - (SEED/1048576)*1048576
      RNDX = REAL(SEED+1)/1048577.0

      RETURN
      END


!***********************************************************************
!***********************************************************************
!*                                                                    **
!* S E G M E N T S   W I T H   C O N N E C T I N G   B R A N C H E S  **
!*                                                                    **
!***********************************************************************
!***********************************************************************
!  This Subroutine Calls No Other Subroutines

      SUBROUTINE FIMPBR

      Use FISHY; Use GEOMC; USE GLOBAL;USE SURFHE
      
      IMPLICIT NONE
      REAL    :: BRANGLE

! *** ONLY HERE IS PHI0 used SW 1/14/01, DHS is used here also

      DO 105 I=1,IMX                                !For all segments in the system
        NR = 0                                      ! Number of branches connecting to the right bank
        NL = 0                                      ! Number of branches connecting to the left bank
        RIMPBR(I,1) = 0                             ! Right bank joining segment
        RIMPBR(I,2) = 0                             ! Right bank joining branch
        LIMPBR(I,1) = 0                             ! Left bank joining segment
        LIMPBR(I,2) = 0                             ! Left bank joining branch
        DO 100 JB=1,NBR                             !Examine all branches in the system
        IF(BR_INACTIVE(JB))CYCLE  ! SW 6/12/2017
          IF (I.EQ.DHS(JB)) THEN                    ! If find a connecting branch
            IF (PHI0(DS(JB)).GT.PHI0(I)) THEN
              BRANGLE = PHI0(DS(JB))-PHI0(I)        ! BRANGLE = Angle of Incoming Branch relative to
            ELSE                                    !           Branch I
              BRANGLE = PHI0(DS(JB))+(2*(22/7)-PHI0(I))
            END IF
            IF ((BRANGLE.GT.0).AND.(BRANGLE.LT.(22/7))) THEN  ! From BRANGLE, one can determine if the
              NR = NR + 1                                     !   Incoming Branch is on the right or
              RIMPBR(I,1) = DS(JB)                            !   left of Branch I
              RIMPBR(I,2) = JB
            ELSE
              NL = NL + 1
              LIMPBR(I,1) = DS(JB) 
              LIMPBR(I,2) = JB
            END IF
          ELSE
          END IF
  100   CONTINUE
        IF ((NR.GT.1).OR.(NL.GT.1)) THEN                        ! CHECK: Can only have one Branch join a
          WRITE(*,*) 'ERROR: More than 1 branch joining segment !   segment from each side&
                      from the same side'
          WRITE(*,9140) I
 9140     FORMAT('ERROR occurred at segment ',I6)
          STOP
        ELSE
        END IF
  105 CONTINUE

      RETURN
      END


!***********************************************************************
!***********************************************************************
!*                                                                    **
!*   F I N D   C U R R E N T   B R A N C H   O F   L O C A T I O N    **
!*                                                                    **
!***********************************************************************
!***********************************************************************
!  This Subroutine Calls No Other Subroutines


      SUBROUTINE FINDNEWBR
      
    Use Fishy; Use GEOMC; USE GLOBAL
    
    IMPLICIT NONE
    integer trib

      TRIB = 1                                              ! TRIB = Tributary (or Branch)
   80 IF ((FIMP.GE.US(TRIB)).AND.(FIMP.LE.DS(TRIB))) THEN
        FNBP = TRIB
      ELSE
        TRIB = TRIB + 1
        IF (TRIB.GT.NBR) THEN                                 ! CHECK: TRIB can't be more
          WRITE(*,*) 'ERROR: Finding new FNBP in Subroutine'  !        than NBR # BRANCHES
          STOP
        ELSE
          GOTO 80
        END IF
      END IF

      RETURN
      END


!***********************************************************************
!***********************************************************************
!*                                                                    **
!*        F I N D   W A T E R B O D Y   O F   L O C A T I O N         **
!*                                                                    **
!***********************************************************************
!***********************************************************************
!  This Subroutine Calls No Other Subroutines

      SUBROUTINE WHATJR
      
      Use FISHY
      USE GLOBAL
      IMPLICIT NONE
      
      INTEGER  WBDY

      WBDY = 0                                              ! WBDY = Water Body
   90 WBDY = WBDY + 1
      IF (WBDY.GT.NWB) THEN                                 ! CHECK: WBDY can't be more
        WRITE(*,*) 'ERROR: Finding current JR in Subroutine'!        than NWB
        STOP
      ELSE
      END IF
      DO 92 JB=BS(WBDY),BE(WBDY)                            ! Branches within a particular
        IF ((FIMP.GE.US(JB)).AND.(FIMP.LE.DS(JB))) THEN     !   Water Body
          FJR = WBDY
          GOTO 91
        ELSE
        END IF
   92 CONTINUE
      GOTO 90

   91 CONTINUE

      RETURN
      END


!***********************************************************************
!***********************************************************************
!*                                                                    **
!*   I N F L U E N C E   O F   O N E   J O I N I N G   B R A N C H    **
!*                                                                    **
!***********************************************************************
!***********************************************************************
!  This Subroutine Calls the Following Subroutines:
!      FINDNEWBR


      SUBROUTINE TAG124578
      Use FISHY;USE GLOBAL; Use GEOMC
      IMPLICIT NONE
      
      REAL       USBRFLOW,FLOWHERE,CBID,CBIDL,CBIDR
      REAL       WIDTH,TOTBRAREA,AVEBRVEL

      INTEGER    USBRIMP,FKMPTMP
      INTEGER    FIMPTMP

      IF (TAG.EQ.1) USBRIMP = LIMPBR(FIMP+1,1)     ! USBRIMP = UpStream left BRanch IMP
      IF (TAG.EQ.2) USBRIMP = RIMPBR(FIMP+1,1)     ! USBRIMP = UpStream right BRanch IMP
      IF (TAG.EQ.4) USBRIMP = LIMPBR(FIMP-1,1)     ! USBRIMP = UpStream left BRanch IMP
      IF (TAG.EQ.5) USBRIMP = RIMPBR(FIMP-1,1)     ! USBRIMP = UpStream right BRanch IMP
      IF (TAG.EQ.7) USBRIMP = LIMPBR(FIMP,1)       ! USBRIMP = UpStream left BRanch IMP
      IF (TAG.EQ.8) USBRIMP = RIMPBR(FIMP,1)       ! USBRIMP = UpStream right BRanch IMP
      
      IF (DEBUG) THEN
        WRITE(DATADEBUGFN,9170) TAG,FIMP,USBRIMP
 9170   FORMAT(' TAG=',I4,' INITIAL FIMP=',I6,' USBRIMP=',I6)
      ELSE
      END IF
      FKMPTMP = KTI(USBRIMP)
      USBRFLOW = 0.0                               ! USBRFLOW = UpStream (incoming) BRanch FLOW
      TOTBRAREA = 0.0                              ! X-Sectional Area of incoming branch
  102 IF (FKMPTMP.LE.KB(USBRIMP)) THEN
        USBRFLOW =  USBRFLOW + U(FKMPTMP,USBRIMP)*(B(FKMPTMP,USBRIMP)*H(FKMPTMP,FJR))
        TOTBRAREA = TOTBRAREA + B(FKMPTMP,USBRIMP)*H(FKMPTMP,FJR)
        FKMPTMP = FKMPTMP + 1
        GOTO 102
      ELSE
      END IF
      IF ((TAG.EQ.1).OR.(TAG.EQ.2)) THEN
        FKMPTMP = KTI(FIMP+1)
        FIMPTMP = FIMP+1
      ELSE IF ((TAG.EQ.4).OR.(TAG.EQ.5)) THEN
        FKMPTMP = KTI(FIMP-1)
        FIMPTMP = FIMP-1
      ELSE IF ((TAG.EQ.7).OR.(TAG.EQ.8)) THEN
        FKMPTMP = KTI(FIMP)
        FIMPTMP = FIMP
      ELSE
      END IF
      FLOWHERE = 0.0                                              ! FLOWHERE = FLOW in current branch
  101 IF (FKMPTMP.LE.KB(FIMPTMP)) THEN
        FLOWHERE = FLOWHERE + ABS(U(FKMPTMP,FIMPTMP))*(B(FKMPTMP,FIMPTMP)*H(FKMPTMP,FJR))
        FKMPTMP = FKMPTMP + 1
        GOTO 101
      ELSE
      END IF
      FKMPTMP = FKMP                                              ! FKMPTMP=where fish is. FIMPTMP=where fish will be.
  109 IF((FKMPTMP.LT.KTI(FIMPTMP)).OR.(FKMPTMP.GT.KB(FIMPTMP)))THEN   ! Check to make sure there is an active cell
        IF (DEBUG) WRITE(DATADEBUGFN,*) ' B EQ 0 CHECK ACTIVATED (1)'    !   where calculating CBID
        IF (FKMPTMP.LT.KTI(FIMPTMP)) THEN
          FKMPTMP = FKMPTMP + 1
        ELSE IF (FKMPTMP.GT.KB(FIMPTMP)) THEN
          FKMPTMP = FKMPTMP - 1
        ELSE
          WRITE(*,*) 'ERROR: Problem finding where to calculate CBID'
        END IF
        GOTO 109
      ELSE
      END IF
      CBID = ABS(USBRFLOW)/(ABS(USBRFLOW)+FLOWHERE)*&   ! CBID = Combining Branch Influence Distance
             B(FKMPTMP,FIMPTMP)                        !        into current branch
      AVEBRVEL = USBRFLOW/TOTBRAREA                    ! AVEBRVEL = Average BRanch VELocity, = Flow/X-sect Area
      IF (DEBUG) THEN
        WRITE(DATADEBUGFN,9180) B(FKMPTMP,FIMPTMP),CBID,AVEBRVEL,FYLOC
 9180   FORMAT(' B=',F9.2,' CBID=',F9.2,' AVEBRVEL=',E9.3,' INITIAL FYLOC=',F9.3)
      ELSE
      END IF
      IF ((TAG.EQ.1).OR.(TAG.EQ.4).OR.(TAG.EQ.7)) THEN                 !Incoming left branch
        CBIDL = CBID                                                   ! CBID from the left bank
        IF ((FYLOC.GT.0).AND.(FYLOC.LT.CBIDL)) THEN
          FYLOC = FYLOC + (AVEBRVEL+FYVEL+VFISH)*DLT + RDY             ! New updated Fish Y-LOCation
          IF (DEBUG) WRITE(DATADEBUGFN,*) ' FYLOC INFLUENCED BY LEFT BRANCH'
        ELSE
        END IF
      ELSE IF ((TAG.EQ.2).OR.(TAG.EQ.5).OR.(TAG.EQ.8)) THEN            !Incoming right branch
        CBIDR = B(FKMPTMP,FIMPTMP) - CBID                              ! CBID from the right bank
        IF ((FYLOC.GT.CBIDR).AND.(FYLOC.LT.B(FKMPTMP,FIMPTMP))) THEN
          FYLOC = FYLOC - (AVEBRVEL-FYVEL-VFISH)*DLT + RDY             ! New updated Fish Y-LOCation
          IF (DEBUG) WRITE(DATADEBUGFN,*) ' FYLOC INFLUENCED BY RIGHT BRANCH'
        ELSE
        END IF
      ELSE
      END IF
      IF (DEBUG) WRITE(DATADEBUGFN,9190) FYLOC
 9190 FORMAT(' THE NEW FYLOC=',F9.3)

! Calculating New Fish Location: FIMP, FXLOC, FYLOC, and FNBP (and FKMP, if necessary) 

      WIDTH = B(FKMPTMP,FIMPTMP)
      IF ((FYLOC.LT.0).AND.((TAG.EQ.1).OR.(TAG.EQ.4).OR.(TAG.EQ.7))) THEN          ! ((CASE 1))
          IF (DEBUG) WRITE(DATADEBUGFN,*) ' ((CASE 1)) ACTIVATED'
          FYLOCTMP = FYLOC
          FIMP = USBRIMP
  110     IF ((FKMPTMP.LT.KTI(FIMP)).OR.&
              (FKMPTMP.GT.KB(FIMP))) THEN                           ! Check to make sure there is an active cell
            IF (DEBUG) WRITE(DATADEBUGFN,*) ' B EQ 0 CHECK ACTIVATED (2)'  !   where calculating new fish location
            IF (FKMPTMP.LT.KTI(FIMP)) THEN
              FKMPTMP = FKMPTMP + 1
            ELSE IF (FKMPTMP.GT.KB(FIMP)) THEN
              FKMPTMP = FKMPTMP - 1
            ELSE
              WRITE(*,*) 'ERROR: Problem finding upstream KMP for fish in left upstream branch'
              STOP
            END IF
            GOTO 110
          ELSE
          END IF
          FYLOC = B(FKMPTMP,FIMP)/2                                 ! When fish moves into upstream branch I
          CALL FINDNEWBR                                            !   arbitrarily set its new FYLOC as the middle
          FXLOC = DLX(FIMP) + FYLOCTMP                              ! A negative FXLOC will be caught by the Subroutine FISH
      ELSE IF (FYLOC.LT.0) THEN                                     ! ((CASE 2))
          IF (DEBUG) WRITE(DATADEBUGFN,*) ' ((CASE 2)) ACTIVATED'
          FYLOC = B(FKMPTMP,FIMPTMP)*YREFL                          ! Logic for rebounding off of left bank
          IF (TAG.EQ.2) THEN                                        ! Moving downstream to FIMP+1
            FXLOC = FXLOC - DLX(FIMP)
            FIMP = FIMPTMP
          ELSE IF (TAG.EQ.5) THEN                                   ! Moving upstream to FIMP-1
            FIMP = FIMPTMP
            FXLOC = DLX(FIMP) + FXLOC
          ELSE IF (TAG.EQ.8) THEN
            ! NOTHING
          ELSE
            WRITE(*,*) 'ERROR: Unable to calculate fish movement off of left bank'
            STOP
          END IF
      ELSE IF ((FYLOC.GT.WIDTH).AND.((TAG.EQ.2).OR.&
                                    (TAG.EQ.5).OR.(TAG.EQ.8))) THEN ! ((CASE 3))
          IF (DEBUG) WRITE(DATADEBUGFN,*) ' ((CASE 3)) ACTIVATED'
          FYLOCTMP = FYLOC
          FIMP = USBRIMP
  111     IF ((FKMPTMP.LT.KTI(FIMP)).OR.&
              (FKMPTMP.GT.KB(FIMP))) THEN                           ! Check to make sure there is an active cell
            IF (DEBUG) WRITE(DATADEBUGFN,*) ' B EQ 0 CHECK ACTIVATED (3)'  !   where calculating new fish location
            IF (FKMPTMP.LT.KTI(FIMP)) THEN
              FKMPTMP = FKMPTMP + 1
            ELSE IF (FKMPTMP.GT.KB(FIMP)) THEN
              FKMPTMP = FKMPTMP - 1
            ELSE
              WRITE(*,*) 'ERROR: Problem finding upstream KMP for fish in right upstream branch'
              STOP
            END IF
            GOTO 111
          ELSE
          END IF
          FYLOC = B(FKMPTMP,FIMP)/2                                 ! When fish moves into upstream branch I
          CALL FINDNEWBR                                            !   arbitrarily set its new FYLOC as the middle
          FXLOC = DLX(FIMP)-(FYLOCTMP-WIDTH)                        ! A negative FXLOC will be caught by the Subroutine FISH
      ELSE IF (FYLOC.GT.WIDTH) THEN                                 ! ((CASE 4))
          IF (DEBUG) WRITE(DATADEBUGFN,*) ' ((CASE 4)) ACTIVATED'
          FYLOC = B(FKMPTMP,FIMPTMP)*(1-YREFL)                      ! Logic for rebounding off of right bank
          IF (TAG.EQ.1) THEN                                        ! Moving downstream to FIMP+1
            FXLOC = FXLOC - DLX(FIMP)
            FIMP = FIMPTMP
          ELSE IF (TAG.EQ.4) THEN                                   ! Moving upstream to FIMP-1
            FIMP = FIMPTMP
            FXLOC = DLX(FIMP) + FXLOC
          ELSE IF (TAG.EQ.7) THEN
            ! NOTHING
          ELSE
            WRITE(*,*) 'ERROR: Unable to calculate fish movement off of right bank'
            STOP
          END IF
      ELSE                                                          ! ((CASE 5))
          IF (DEBUG) WRITE(DATADEBUGFN,*) ' ((CASE 5)) ACTIVATED'
          IF ((TAG.EQ.1).OR.(TAG.EQ.2)) THEN                        ! Moving downstream to FIMP+1
            FXLOC = FXLOC - DLX(FIMP)
            FIMP = FIMPTMP
          ELSE IF ((TAG.EQ.4).OR.(TAG.EQ.5)) THEN                   ! Moving upstream to FIMP-1
            FIMP = FIMPTMP
            FXLOC = DLX(FIMP) + FXLOC
          ELSE IF ((TAG.EQ.7).OR.(TAG.EQ.8)) THEN
            ! NOTHING
          ELSE
            WRITE(*,*) 'ERROR: Unable to calculate fish FXLOC in TAG124578'
            STOP
          END IF
      END IF
      IF (DEBUG) THEN
        WRITE(DATADEBUGFN,9200) FNBP,FIMP,FXLOC,FYLOC
 9200   FORMAT(' NEW: FNBP=',I6,' FIMP=',I6,' FXLOC=',F9.2,' FYLOC=',F9.3)
      ELSE
      END IF

      RETURN
      END


!***********************************************************************
!***********************************************************************
!*                                                                    **
!*  I N F L U E N C E   O F   T W O   J O I N I N G   B R A N C H E S **
!*                                                                    **
!***********************************************************************
!***********************************************************************
!  This Subroutine Calls the Following Subroutines:
!      FINDNEWBR


      SUBROUTINE TAG369

      
      Use FISHY
      USE GLOBAL
      Use GEOMC
      
      
      IMPLICIT NONE

      REAL       USBRFLOW,FLOWHERE,CBID,CBIDR
      REAL       WIDTH,TOTBRAREA,AVEBRVEL

      INTEGER    USBRIMP,FKMPTMP,JBR
      INTEGER    FIMPTMP

      DIMENSION  USBRIMP(2),FKMPTMP(2),USBRFLOW(2),TOTBRAREA(2),CBID(2),AVEBRVEL(2)

      IF (TAG.EQ.3) THEN
        USBRIMP(1) = LIMPBR(FIMP+1,1)              ! USBRIMP(1) = UpStream left BRanch IMP
        USBRIMP(2) = RIMPBR(FIMP+1,1)              ! USBRIMP(2) = UpStream right BRanch IMP
      ELSE IF (TAG.EQ.6) THEN
        USBRIMP(1) = LIMPBR(FIMP-1,1)              ! USBRIMP(1) = UpStream left BRanch IMP
        USBRIMP(2) = RIMPBR(FIMP-1,1)              ! USBRIMP(2) = UpStream right BRanch IMP
      ELSE IF (TAG.EQ.9) THEN
        USBRIMP(1) = LIMPBR(FIMP,1)                ! USBRIMP(1) = UpStream left BRanch IMP
        USBRIMP(2) = RIMPBR(FIMP,1)                ! USBRIMP(2) = UpStream right BRanch IMP
      ELSE
      END IF
      FKMPTMP(1) = KTI(USBRIMP(1))                 ! Surface Layer of incoming left branch
      FKMPTMP(2) = KTI(USBRIMP(2))                 ! Surface Layer of incoming right branch
      DO 106 JBR=1,2                                ! BR: =1 incoming left branch, =2 incoming right branch
        USBRFLOW(JBR) = 0.0                         ! USBRFLOW = UpStream (incoming) BRanch FLOW
        TOTBRAREA(JBR) = 0.0                        ! X-Sectional Area of incoming branch
  103   IF (FKMPTMP(JBR).LE.KB(USBRIMP(JBR))) THEN
          USBRFLOW(JBR) =  USBRFLOW(JBR) + U(FKMPTMP(JBR),USBRIMP(JBR))*&
                         (B(FKMPTMP(JBR),USBRIMP(JBR))*H(FKMPTMP(JBR),FJR))
          TOTBRAREA(JBR) = TOTBRAREA(JBR) +&
                          B(FKMPTMP(JBR),USBRIMP(JBR))*H(FKMPTMP(JBR),FJR)
          FKMPTMP(JBR) = FKMPTMP(JBR) + 1
          GOTO 103
        ELSE
        END IF
  106 CONTINUE
      IF (TAG.EQ.3) THEN
        FKMPTEMP = KTI(FIMP+1)
        FIMPTMP = FIMP+1
      ELSE IF (TAG.EQ.6) THEN
        FKMPTEMP = KTI(FIMP-1)
        FIMPTMP = FIMP-1
      ELSE IF (TAG.EQ.9) THEN
        FKMPTEMP = KTI(FIMP)
        FIMPTMP = FIMP
      ELSE
      END IF
      FLOWHERE = 0.0                                      ! FLOWHERE = FLOW in current branch
  107 IF (FKMPTEMP.LE.KB(FIMPTMP)) THEN
        FLOWHERE = FLOWHERE + ABS(U(FKMPTEMP,FIMPTMP))*&
                  (B(FKMPTEMP,FIMPTMP)*H(FKMPTEMP,FJR))
        FKMPTEMP = FKMPTEMP + 1
        GOTO 107
      ELSE
      END IF
      FKMPTEMP = FKMP                                     ! FKMPTEMP=where fish is. FIMPTMP=where fish will be.
  108 IF ((FKMPTEMP.LT.KTI(FIMPTMP)).OR.&
          (FKMPTEMP.GT.KB(FIMPTMP))) THEN                 ! Check to make sure there is an active cell
        IF (FKMPTEMP.LT.KTI(FIMPTMP)) THEN                !   where calculating CBID
          FKMPTEMP = FKMPTEMP + 1
        ELSE IF (FKMPTEMP.GT.KB(FIMPTMP)) THEN
          FKMPTEMP = FKMPTEMP - 1
        ELSE
          WRITE(*,*) 'ERROR: Problem finding where to&
                             calculate CBID with 2&
                             incoming branches'
        END IF
        GOTO 108
      ELSE
      END IF
      DO 112 JBR=1,2
        CBID(JBR) = ABS(USBRFLOW(JBR))/((ABS(USBRFLOW(1))+&
                   ABS(USBRFLOW(2))+FLOWHERE)*&
                   B(FKMPTEMP,FIMPTMP))                   ! CBID = Combining Branch Influence Distance
        AVEBRVEL(JBR) = USBRFLOW(JBR)/TOTBRAREA(JBR)         ! AVEBRVEL = Average BRanch VELocity, = Flow/X-sect Area
  112 CONTINUE
      WIDTH = B(FKMPTEMP,FIMPTMP)
      IF ((FYLOC.GT.0).AND.(FYLOC.LT.CBID(1))) THEN
        FYLOC = FYLOC + (AVEBRVEL(1)+FYVEL+&
                         VFISH)*DLT + RDY                 ! New updated Fish Y-LOCation
      ELSE
        CBIDR = WIDTH - CBID(2)                           ! CBIDR = CBID from the right bank
        IF ((FYLOC.GT.CBIDR).AND.(FYLOC.LT.WIDTH)) THEN
          FYLOC = FYLOC - (AVEBRVEL(2)-FYVEL-&
                           VFISH)*DLT + RDY               ! New updated Fish Y-LOCation
        ELSE
        END IF
      END IF

! Calculating New Fish Location: FIMP, FXLOC, FYLOC, and FNBP (and FKMP, if necessary) 

      IF (FYLOC.LT.0) THEN
          FYLOCTMP = FYLOC
          FIMP = USBRIMP(1)
  113     IF ((FKMPTEMP.LT.KTI(FIMP)).OR.&
              (FKMPTEMP.GT.KB(FIMP))) THEN                ! Check to make sure there is an active cell
            IF (FKMPTEMP.LT.KTI(FIMP)) THEN               !   where calculating new fish location
              FKMPTEMP = FKMPTEMP + 1
            ELSE IF (FKMPTEMP.GT.KB(FIMP)) THEN
              FKMPTEMP = FKMPTEMP - 1
            ELSE
              WRITE(*,*) 'ERROR: Problem finding upstream&
                          KMP for fish in left upstream&
                          branch with 2 incoming branches'
              STOP
            END IF
            GOTO 113
          ELSE
          END IF
          FYLOC = B(FKMPTEMP,FIMP)/2                      ! When fish moves into upstream branch I
          CALL FINDNEWBR                                  !   arbitrarily set its new FYLOC as the middle
          FXLOC = DLX(FIMP) + FYLOCTMP                    ! A negative FXLOC will be caught by the Subroutine FISH
      ELSE IF (FYLOC.GT.WIDTH) THEN
          FYLOCTMP = FYLOC
          FIMP = USBRIMP(2)
  114     IF ((FKMPTEMP.LT.KTI(FIMP)).OR.&
              (FKMPTEMP.GT.KB(FIMP))) THEN                ! Check to make sure there is an active cell
            IF (FKMPTEMP.LT.KTI(FIMP)) THEN               !   where calculating new fish location
              FKMPTEMP = FKMPTEMP + 1
            ELSE IF (FKMPTEMP.GT.KB(FIMP)) THEN
              FKMPTEMP = FKMPTEMP - 1
            ELSE
              WRITE(*,*) 'ERROR: Problem finding upstream&
                          KMP for fish in right upstream&
                          branch with 2 incoming branches'
              STOP
            END IF
            GOTO 114
          ELSE
          END IF
          FYLOC = B(FKMPTEMP,FIMP)/2                      ! When fish moves into upstream branch I
          CALL FINDNEWBR                                  !   arbitrarily set its new FYLOC as the middle
          FXLOC = DLX(FIMP)-(FYLOCTMP-WIDTH)              ! A negative FXLOC will be caught by the Subroutine FISH
      ELSE
          IF (TAG.EQ.3) THEN                              ! Moving downstream to FIMP+1
            FXLOC = FXLOC - DLX(FIMP)
            FIMP = FIMPTMP
          ELSE IF (TAG.EQ.6) THEN                         ! Moving upstream to FIMP-1
            FIMP = FIMPTMP
            FXLOC = DLX(FIMP) + FXLOC
          ELSE IF (TAG.EQ.9) THEN
            ! NOTHING
          ELSE
            WRITE(*,*) 'ERROR: Unable to calculate fish&
                        FXLOC in TAG369'
            STOP
          END IF
      END IF

      RETURN
      END


!***********************************************************************
!***********************************************************************
!*                                                                    **
!*   P R E P A R E   S Y S T E M   G R I D   F O R   D I S P L A Y    **
!*                                                                    **
!***********************************************************************
!***********************************************************************
! Number of TecPlot Zones       = # of timesteps            (Animation by Zone)
!           TecPlot Limit     ==> 32,700 per Data Set
! Number of TecPlot Data Points = # of nodes                (for Grid Display)
!                                 # of fish                 (for Fish Display)
!           TecPlot Limit     ==> 2 billion per Variable
! Number of TecPlot Variables   = # of flow & const. param  (for Grid Display)
!                                 # of fish parameters      (for Fish Display)
!           TecPlot Limit     ==> 32,700 per Data Set
! Number of TecPlot Data Sets   = 2*(# of branches)         (for Grid & Fish Displays)
!           TecPlot Limit     ==> 128

!  This Subroutine Calls the Following Subroutines:
!      INTERCONST,INTERFLOWF


      SUBROUTINE GRIDPLOT    ! This Subroutine preps Grid, Flow, & WQ info for TecPlot

      Use FISHY
      USE GLOBAL
      Use GEOMC
      Use SCREENC
      
      IMPLICIT NONE

      INTEGER     WB,NN,NNLAST,WATER
      INTEGER     BRCHNN,KK,II,ND,ELE,EE,K
      INTEGER     TOTNN
      INTEGER     JOININGTRIB
      integer, allocatable, save, DIMENSION (:) :: TOTELE,TOTBRNN
      logical, SAVE  ::  FEGRID
      REAL :: XDIST,ZDIST,VERTFLOW

      CHARACTER*1  DTYP1
      CHARACTER*2  DTYP2
      CHARACTER*6  FILE_PREFIX
      CHARACTER*4  FILE_SUFFIX
      CHARACTER*12 FILENAME
      DATA         FILE_PREFIX/ 'Branch'/
      DATA         FILE_SUFFIX/ '.dat'/

      IF (NIT.EQ.0)then
      FEGRID = .FALSE.                  ! The first pass is used to set up the FE grid
      allocate(totele(NBR),totbrnn(NBR))    ! SW 1/14/01
      end if

   93 IF (FEGRID) CALL INTERCONST                     ! Subroutine to interpolate Constituent values
      IF (FEGRID) CALL INTERFLOWF                     ! Subroutine to interpolate Flow Field values
      NN = 0                                          ! NN = Global Node # (UpperLeft Corner of cell(K,I))
      DO 115 WB=1,NWB                                 ! WB = Water Body # / NWB = Total # of Water Bodies
        DO 116 JB=BS(WB),BE(WB)                     ! JB = Branch #
          IF (FEGRID.EQV..FALSE.) GOTO 94             ! Skip to setting up FE grid if NIT = 0
          IF (JB.LT.10) THEN                        ! Creating the output file names based on Branch #
            WRITE(DTYP1,9210) JB
            FILENAME=FILE_PREFIX//DTYP1//FILE_SUFFIX
 9210       FORMAT(I1)
          ELSE IF (JB.LE.64) THEN                   ! TecPlot is limited to 128 Data Sets (i.e. Data Files)
            WRITE(DTYP2,9220) JB                    !   That leaves: 64 Data Sets for Grid Display and
            FILENAME=FILE_PREFIX//DTYP2//FILE_SUFFIX  !                64 Data Sets for Fish Display
 9220       FORMAT(I2)
          ELSE                                        ! CHECK
            WRITE(*,*) 'ERROR: Exceeded Maximum Number of TecPlot Data Sets ==> # Data Sets = # of Branches'
            STOP
          END IF
          IF (NIT.EQ.0) THEN
            OPEN(20000+JB,FILE=FILENAME,STATUS='UNKNOWN')       ! Opening the Output Files Based on Branch #
            WRITE(20000+JB,9240) JB                           ! Creating Output File Header for TecPlot
 9240       FORMAT('TITLE = "GRID NODE INFO for Branch',I6,'"')
            WRITE(20000+JB,9270)                                ! Creating Output File Header for TecPlot
 9270       FORMAT('VARIABLES = "X", "Z", "WATER", "HorizVel", "VertVel"&
      , "Temp", "DO"')
            WRITE(20000+JB,9250) JDAY,TOTBRNN(JB),TOTELE(JB)! Creating Output File Header for TecPlot
 9250       FORMAT('ZONE T="JDAY ',F9.2,'", N=',I8,', E=',I8,&
      ', F=FEPOINT, ET=QUADRILATERAL')
          ELSE                                                    ! Subsequent Passes thru Subroutine
            WRITE(20000+JB,*) ' '
            WRITE(20000+JB,9280) JDAY,TOTBRNN(JB),TOTELE(JB)
 9280       FORMAT('ZONE T="JDAY ',F9.2,'", N=',I8,', E=',I8,&
      ', F=FEPOINT, ET=QUADRILATERAL, D=(FECONNECT)')
          END IF
   94     ZDIST = H(2,WB)                                            ! Begin Calculations to set up FE grid
          NNLAST = NN
	if(nnlast.eq.0)nnlast=1    ! sw 1/9/01
          BRCHNN = 0                                ! BRCHNN = Branch Node # (each Branch starts off with '0')
          DO 117 K=2,KMX                                          ! K = Layer #
         ! DO K=2,KMX
            ZDIST = ZDIST - H(K,WB)                                  ! ZDIST = Depth to Node in Branch JB
            XDIST = 0                                             ! XDIST = X Distance to Node in Branch JB
             DO 118  I=US(JB),DS(JB)+1         ! I = Segment #!  DO I=CUS(JB),DS(JB)+1             ! DO 118  I=US(JB),DS(JB)+1         ! I = Segment #
              !DO K=KTWB(WB),KB(I)   ! SW 5/2015
              IF ((K-1.LE.KB(I)).OR.(K-1.LE.KB(I-1))) THEN        ! All Nodes at or above Water Body Bottom
                NN = NN + 1                         ! NN = Global Node # (UpperLeft Corner of cell(K,I))
                BRCHNN = BRCHNN + 1                 ! BRCHNN = Branch Node # (each Branch starts off with '0')
                IF (FEGRID.EQV..FALSE.) THEN        ! This info is stored for later use
                  NDINFO(NN,1) = NN
                  NDINFO(NN,2) = BRCHNN
                  NDINFO(NN,3) = K
                  NDINFO(NN,4) = I
                  NODES(K,I,1) = XDIST
                  NODES(K,I,2) = ZDIST
                  TOTNN = NN
              !  ELSE
                END IF
                IF (FEGRID) THEN
                  IF (K.LT.KTWB(WB)) THEN                               ! WATER used for value-blanking cells above 9952 surface in TecPlot
                    IF (SHOWSKY) THEN                                  ! Show Colors of Day Above Water Surface
                      IF ((MILTIME.GT.SKYNIGHT).OR.&
                          (MILTIME.LE.SKYDAWN)) THEN                   ! During the Night, the Sky is
                        WATER = -100                                   !   BLACK
                      ELSE IF ((MILTIME.GT.SKYDAWN).AND.&
                               (MILTIME.LE.SKYDAY)) THEN               ! During the Morning, the Sky is
                        WATER = -75                                    !   YELLOW
                      ELSE IF ((MILTIME.GT.SKYDAY).AND.&
                               (MILTIME.LE.SKYDUSK)) THEN              ! During the Day, the Sky is
                        WATER = -50                                    !   BLUE
                      ELSE IF ((MILTIME.GT.SKYDUSK).AND.&
                               (MILTIME.LE.SKYNIGHT)) THEN             ! During the Evening, the Sky is
                        WATER = -25                                    !   ORANGE
                      ELSE                                             ! CHECK
                        WRITE(DATADEBUGFN,*) 'ERROR: MILTIME not Corresponding with Intervals Set for Determining Sky Color for TecPlot'
                        STOP
                      END IF
                    ELSE
                      WATER = -1                                       ! WATER = -1 ==> Cells above water surface are WHITE
                    END IF
                    WQFIELD(K,I,1) = WQFIELD(KTWB(WB),I,1) ! This is done to prevent odd-looking contours in TecPlot near
                    WQFIELD(K,I,2) = WQFIELD(KTWB(WB),I,2) !   the water surface since WQ does not exist above water surface
                  ELSE
                    WATER = 1                                          ! WATER = 1  ==> cell at or below current water surface
                  END IF
                  IF (WBSKIP.AND.(WB.NE.WBRUN)) GOTO 136               ! Option:Skip output for Water Bodies, except WBRUN
                  IF (ABS(FLOWFIELD(K,I,2)).GE.VVELCAP) THEN           ! This keeps relatively large Vert Vel vectors in TecPlot
                    IF (FLOWFIELD(K,I,2).GE.0) THEN                    !   from blacking out the screen since these vectors
                      VERTFLOW = VVELCAP                               !   would be so big
                    ELSE
                      VERTFLOW = -VVELCAP
                    END IF
                  ELSE
                    VERTFLOW = FLOWFIELD(K,I,2)
                  END IF
!                  JOININGTRIB = 0
!                  IF (RIMPBR(I,1).GT.0) JOININGTRIB = 1                ! Branch joining right bank at Segment I
!                  IF (LIMPBR(I,1).GT.0) JOININGTRIB = 1                ! Branch joining left bank at Segment I
!                  WRITE(20000+JB,9260) XDIST,ZDIST,WATER,&            ! Outputting node information
                  if(k==ktwb(wb))then
                      zdist=elws(i)
                  else
                      zdist=elws(i)-depthm(k,i)
                  endif
                  
                WRITE(20000+JB,9260) XDIST,zdist,WATER,&            ! Outputting node information  SW output elevation in m
                        FLOWFIELD(K,I,1),-VERTFLOW/ASPRATIO,&           !   to a file for TecPlot to display
                        WQFIELD(K,I,1),WQFIELD(K,I,2)
!     .                  B(K,I),JOININGTRIB
 9260             FORMAT(F15.1,F10.2,I7,E15.2,E15.2,F8.2,F8.2)
!     .                   F10.2,I5)
                                                 ! XDIST                 = X Distance to Node in Branch JB
                                                 ! ZDIST                 = Depth to Node in Branch JB
                                                 ! NN                    = Global Node # 
                                                 ! BRCHNN                = Branch Node #
                                                 ! K                     = Layer #
                                                 ! I                     = Segment #
                                                 ! JB                  = Branch #
                                                 ! WATER                 = Water in cell? [-1=No , 1=Yes]
                                                 ! FLOWFIELD(K,I,1)      = Velocity: Horizontal (m/s)
                                                 ! FLOWFIELD(K,I,2)                  Vertical   (m/s):W(K,I) is (+) Downward
                                                 ! FLOWFIELD(K,I,3)      = Acceleration: Horizontal (m/s^2)
                                                 ! FLOWFIELD(K,I,4)                      Vertical   (m/s^2)
                                                 ! WQFIELD(K,I,1)        = Temperature              (deg C)
                                                 ! WQFIELD(K,I,2)        = Dissolved Oxygen         (g/m^3)
                                                 ! B(K,I)                = Width of Cell (K,I)
                                                 ! JOININGTRIB           = Incoming Tributary at Segment I? [0=No , 1=Yes]

  136           END IF

              END IF
              XDIST = XDIST + DLX(I)
    !ENDDO
    !ENDDO
    
  118       CONTINUE
  117     CONTINUE
          IF (FEGRID) THEN                                 ! If FEGRID = .TRUE., then FE connectivity already set up
            IF (NIT.EQ.0) THEN                             ! Once FE connectivity is established (below), it still
              WRITE(20000+JB,*) ' '                      !   must be outputted the first time node info is outputted
              IF (WBSKIP.AND.(WB.NE.WBRUN)) GOTO 137       ! Option:Skip output for all Water Bodies, except WBRUN
              DO 122 EE=1,TOTELE(JB)
                WRITE(20000+JB,9290) CORNERS(JB,EE,1),& ! Writing connectivity information to output file
                                       CORNERS(JB,EE,2),&
                                       CORNERS(JB,EE,3),&
                                       CORNERS(JB,EE,4)
 9290           FORMAT(I8,I8,I8,I8)
  122         CONTINUE
            ELSE
  137       END IF
          ELSE                                             ! Determining the Node Connectivities
            ELE = 0                                        ! ELE = Element #
            DO 119 KK=2,KMX                                ! KK = Layer #
              DO 120 II=US(JB),DS(JB)                  ! II = Segment #
                IF (KK.LE.KB(II)) THEN                     ! Any Layer above Water Body Bottom
                  ELE = ELE + 1
                  DO 121 ND=NNLAST,NN                      ! NNLAST prevents having to scan the entire matrix for
                    IF ((NDINFO(ND,3).EQ.KK).AND.&          !   the desired information needed to determine connectivity
                        (NDINFO(ND,4).EQ.II)) THEN
                      CORNERS(JB,ELE,1) = NDINFO(ND,2)
                    ELSE IF ((NDINFO(ND,3).EQ.KK).AND.&
                             (NDINFO(ND,4).EQ.II+1)) THEN
                      CORNERS(JB,ELE,2) = NDINFO(ND,2)
                    ELSE IF ((NDINFO(ND,3).EQ.KK+1).AND.&
                             (NDINFO(ND,4).EQ.II+1)) THEN
                      CORNERS(JB,ELE,3) = NDINFO(ND,2)
                    ELSE IF ((NDINFO(ND,3).EQ.KK+1).AND.&
                             (NDINFO(ND,4).EQ.II)) THEN
                      CORNERS(JB,ELE,4) = NDINFO(ND,2)
                    ELSE
                    END IF
  121             CONTINUE
                ELSE
                END IF
  120         CONTINUE
  119       CONTINUE
            TOTELE(JB) = ELE                             ! TOTELE(JB) = Total # of elements in Branch JB
            TOTBRNN(JB) = BRCHNN                         ! TOTBRNN(JB) = Total # of nodes in Branch JB
          END IF
  116   CONTINUE
  115 CONTINUE

      IF (FEGRID) THEN
      ELSE
        FEGRID = .TRUE.                                    ! Set FEGRID = .TRUE. once connectivity has been
        GOTO 93                                            !   determined
      END IF

      RETURN
      END


!***********************************************************************
!***********************************************************************
!*                                                                    **
!* I N T E R P O L A T   A L G O R I T H M - C O N S T I T U E N T S  **
!*                                                                    **
!***********************************************************************
!***********************************************************************
!    1                   2                                        3    *
!     O-------------------O----------------------------------------O   *
!     |                   |                                        |   *
!     |                   |                                        |   *
!     |                   |                                        |   *
!     |         +         |                    +                   |   *
!     |                   |                                        |   *
!     |                   |                                        |   *
!    4|                  5|                                       6|   *
!     O-------------------O----------------------------------------O   *
!     |                   |                                        |   *
!     |                   |                                        |   *
!     |         +         |                    +                   |   *
!     |                   |                                        |   *
!    7|                  8|                                       9|   *
!     O-------------------O----------------------------------------O   *
!                                                                      *
!     + = Constituent Information Given in W2                          *
!     O = Node where Information is Needed for TecPlot Purposes        *
!                                                                      *
!***********************************************************************
!  This Subroutine Calls No Other Subroutines

      SUBROUTINE INTERCONST

      Use FISHY; USE GLOBAL; Use GEOMC; Use SCREENC; USE KINETIC
      
      IMPLICIT NONE

      REAL        DV,DW,DX,DY
      REAL        H1T,H2T,H3T,H4T,H1DO,H2DO,H3DO,H4DO
      REAL        XRANGE,YRANGE
      REAL        TOPK,LFTI,RGTI
      REAL        CT,BT,AT
      REAL        CV,BD,AD

      INTEGER     WB,NN,K,KK,II
      INTEGER     TOTNN

      SAVE        TOPK

      NN = 0                                                       ! NN = Global Node # (UpperLeft Corner of cell(K,I))
      DO 127 WB=1,NWB                                              ! WB = Water Body # / NWB = Total # of Water Bodies
        DO 128 JB=BS(WB),BE(WB)                                    ! JB = Branch #
          DO 129 K=2,KMX                                           ! K = Layer #
            DO 130 I=US(JB),(DS(JB)+1)                         ! I = Segment #  changed from us to cus sw 5/2015
              IF ((K-1.LE.KB(I)).OR.(K-1.LE.KB(I-1))) THEN         ! All Nodes at or above Water Body Bottom
                NN = NN + 1                                        ! PRE-CHECK CALCULATION
                KK = NDINFO(NN,3)                                  ! PRE-CHECK CALCULATION
                II = NDINFO(NN,4)                                  ! PRE-CHECK CALCULATION
                IF ((I.NE.II).OR.(K.NE.KK)) THEN                   ! CHECK
                  
                  WRITE(DATADEBUGFN,*) 'ERROR: Nodes not matching up for Constituent Interpolation'
                  write(DATADEBUGFN,*)'JDAY=',jday
                  write(DATADEBUGFN,*)'WB=',wb,' JR=', JB
                  write(DATADEBUGFN,*)'NN=',nn,' kk=',kk,' ii=',ii,' k=',k,' i=',i
                  STOP
                ELSE
                END IF

! LINEAR INTERPOLATION for Boundary Nodes; those analogous to Node 5 use Bilinear Splines below
                IF (K.LT.KTWB(WB)) THEN                                  !Above Water Surface
                  WQFIELD(K,I,1) = 0.0
                  WQFIELD(K,I,2) = 0.0
                ELSE IF ((K.EQ.KTWB(WB)).AND.(I.EQ.US(JB))) THEN       !Analogous to Node 1 above
                  WQFIELD(K,I,1) = T2(K,I)
                  WQFIELD(K,I,2) = O2(K,I)
                ELSE IF ((K.EQ.KTWB(WB)).AND.(I.EQ.DS(JB)+1)) THEN     !Analogous to Node 3 above
                  WQFIELD(K,I,1) = T2(K,I-1)
                  WQFIELD(K,I,2) = O2(K,I-1)
                ELSE IF (K.EQ.KTWB(WB)) THEN                             !Analogous to Node 2 above
                  XRANGE = .5*DLX(I-1)+.5*DLX(I)
                  WQFIELD(K,I,1) = (T2(K,I-1)*(XRANGE-.5*DLX(I-1))+ T2(K,I)*(XRANGE-.5*DLX(I)))/XRANGE
                  WQFIELD(K,I,2) = (O2(K,I-1)*(XRANGE-.5*DLX(I-1))+O2(K,I)*(XRANGE-.5*DLX(I)))/XRANGE
                ELSE IF ((K.EQ.KB(I)+1).AND.(I.EQ.US(JB))) THEN       !Analogous to Node 7 above
                  WQFIELD(K,I,1) = T2(K-1,I)
                  WQFIELD(K,I,2) = O2(K-1,I)
                ELSE IF ((K.EQ.KB(I)+1).AND.(I.EQ.DS(JB)+1)) THEN     !Analogous to Node 9 above
                  WQFIELD(K,I,1) = T2(K-1,I-1)
                  WQFIELD(K,I,2) = O2(K-1,I-1)
                ELSE IF (K.EQ.KB(I)+1) THEN                             !Analogous to Node 8 above
                  XRANGE = .5*DLX(I-1)+.5*DLX(I)
                  WQFIELD(K,I,1) = (T2(K-1,I-1)*(XRANGE-.5*DLX(I-1))+&
                                   T2(K-1,I)*(XRANGE-.5*DLX(I)))/&
                                   XRANGE
                  WQFIELD(K,I,2) = (O2(K-1,I-1)*(XRANGE-.5*DLX(I-1))+&
                                   O2(K-1,I)*(XRANGE-.5*DLX(I)))/&
                                   XRANGE
                ELSE IF (I.EQ.US(JB)) THEN                            !Analogous to Node 4 above
                  YRANGE = .5*H(K-1,WB)+.5*H(K,WB)
                  WQFIELD(K,I,1) = (T2(K-1,I)*(YRANGE-.5*H(K-1,WB))+&
                                    T2(K,I)*(YRANGE-.5*H(K,WB)))/YRANGE
                  WQFIELD(K,I,2) = (O2(K-1,I)*(YRANGE-.5*H(K-1,WB))+&
                                    O2(K,I)*(YRANGE-.5*H(K,WB)))/YRANGE
                ELSE IF (I.EQ.DS(JB)+1) THEN                          !Analogous to Node 6 above
                  YRANGE = .5*H(K-1,WB)+.5*H(K,WB)
                  WQFIELD(K,I,1) = (T2(K-1,I-1)*(YRANGE-.5*H(K-1,WB))+&
                                   T2(K,I-1)*(YRANGE-.5*H(K,WB)))/&
                                   YRANGE
                  WQFIELD(K,I,2) = (O2(K-1,I-1)*(YRANGE-.5*H(K-1,WB))+&
                                   O2(K,I-1)*(YRANGE-.5*H(K,WB)))/&
                                   YRANGE
                ELSE                                               ! Analogous to Node 5 above (i.e., in the middle)
! Begin Bilinear Spline Interpolation Algorithm (Algorithm obtained from "TWO DIMENSIONAL SPLINE INTERPOLATION
!                                                ALGORITHMS - Author: Helmuth Spath")
                  DV = 0.0-(-0.5*DLX(I-1))           ! 0.0 = X-coord of the desired interpolated location
                  DW = 0.0-(-0.5*H(K-1,WB))             ! 0.0 = Y-coord of the desired interpolated location
                  DX = 0.5*DLX(I)-(-0.5*DLX(I-1))
                  DY = 0.5*H(K,WB)-(-0.5*H(K-1,WB))
                  H1T = T2(K-1,I-1)                  ! H1T,H2T,H3T,H4T = Value of Temperature at 4 cell centers as
                  H2T = T2(K-1,I)                    !   established in W2
                  H3T = T2(K,I-1)
                  H4T = T2(K,I)
                  H1DO = O2(K-1,I-1)                 ! H1DO,H2DO,H3DO,H4DO = Value of Diss. Oxyg. at 4 cell centers as
                  H2DO = O2(K-1,I)                   !   established in W2
                  H3DO = O2(K,I-1)
                  H4DO = O2(K,I)
                  WQFIELD(K,I,1) = H1T+(H2T-H1T)/DX*DV+(H3T-H1T)  &                ! Temperature at grid node
                          /DY*DW+(H4T-H3T-H2T+H1T)/(DX*DY)*DV*DW
                  WQFIELD(K,I,2) = H1DO+(H2DO-H1DO)/DX*DV+(H3DO-H1DO) &            ! Dissolved Oxygen at grid node
                          /DY*DW+(H4DO-H3DO-H2DO+H1DO)/(DX*DY)*DV*DW
! End Bilinear Spline Interpolation Algorithm
                END IF

!NULL WATER QUALITY FIELD (if so desired)
                IF (NULLFIELDWQ) THEN                              ! Establish a Null Water Quality Field
                  IF (NIT.EQ.0) TOPK = REAL(KTWB(WB))                ! Null Field established for WB
                  LFTI = REAL(US(JB))
                  RGTI = REAL(DS(JB)+1)
                  IF (VARYTEMP.EQ.1) THEN                          !Vary Temperature Vertically
                    IF (WQNULLLINEAR) THEN                         !  Use Linear Variation in Null Temperature Values
                      WQFIELD(K,I,1) = (TOPTEMP-BOTTEMP)/&
                                       (TOPK-BOTK)*(REAL(K)-TOPK)&
                                       + TOPTEMP
                    ELSE                                           !  Use Parabolic Variation in Null Temperature Values
                      CT = TOPTEMP
                      BT = (MIDKTEMP-TOPTEMP)/(MIDKT-TOPK)
                      AT = ((BOTTEMP-MIDKTEMP)/(BOTK-MIDKT) - BT)/&
                           (BOTK-TOPK)
                      WQFIELD(K,I,1) = AT*(REAL(K))**2+BT*REAL(K)+CT
                    END IF
                  ELSE IF (VARYTEMP.EQ.2) THEN                     !Vary Temperature Horizontally
                    IF (WQNULLLINEAR) THEN                         !  Use Linear Variation in Null Temperature Values
                      WQFIELD(K,I,1) = (LFTTEMP-RGTTEMP)/&
                                       (LFTI-RGTI)*(REAL(I)-LFTI)&
                                       + LFTTEMP
                    ELSE                                           !  Use Parabolic Variation in Null Temperature Values
                      CT = LFTTEMP
                      BT = (MIDITEMP-LFTTEMP)/(MIDIT-LFTI)
                      AT = ((RGTTEMP-MIDITEMP)/(RGTI-MIDIT) - BT)/&
                           (RGTI-LFTI)
                      WQFIELD(K,I,1) = AT*(REAL(I))**2+BT*REAL(I)+CT
                    END IF
                  ELSE
                  END IF
                  IF (VARYDO.EQ.1) THEN                            !Vary Diss. Oxyg. Values Vertically
                    IF (WQNULLLINEAR) THEN                         !  Use Linear Variation in Null Diss. Oxyg. Values
                      WQFIELD(K,I,2) = (TOPDO-BOTDO)/&
                                       (TOPK-BOTK)*(REAL(K)-TOPK)&
                                       + TOPDO
                    ELSE                                           !  Use Parabolic Variation in Null Diss. Oxyg. Values
                      CV = TOPDO
                      BD = (MIDKDO-TOPDO)/(MIDKD-TOPK)
                      AD = ((BOTDO-MIDKDO)/(BOTK-MIDKD) - BD)/&
                           (BOTK-TOPK)
                      WQFIELD(K,I,2) = AD*(REAL(K))**2+BD*REAL(K)+CV
                    END IF
                  ELSE IF (VARYDO.EQ.2) THEN                       !Vary Diss. Oxyg. Values Horizontally
                    IF (WQNULLLINEAR) THEN                         !  Use Linear Variation in Null Diss. Oxyg. Values
                      WQFIELD(K,I,2) = (LFTDO-RGTDO)/&
                                       (LFTI-RGTI)*(REAL(I)-LFTI)&
                                       + LFTDO
                    ELSE                                           !  Use Parabolic Variation in Null Diss. Oxyg. Values
                      CV = LFTDO
                      BD = (MIDIDO-LFTDO)/(MIDID-LFTI)
                      AD = ((RGTDO-MIDIDO)/(RGTI-MIDID) - BD)/&
                           (RGTI-LFTI)
                      WQFIELD(K,I,2) = AD*(REAL(I))**2+BD*REAL(I)+CV
                    END IF
                  ELSE
                  END IF
                ELSE
                END IF

              ELSE
              END IF
  130       CONTINUE
  129     CONTINUE
  128   CONTINUE
  127 CONTINUE

      RETURN
      END


!***********************************************************************
!***********************************************************************
!*                                                                    **
!*  I N T E R P O L A T    A L G O R I T H M  -  F L O W   F I E L D  **
!*                                                                    **
!***********************************************************************
!***********************************************************************
!    1                   2                                        3    *
!     O--------[VV]-------O------------------[VV]------------------O   *
!     |                   |                                        |   *
!     |                   |                                        |   *
!     |                   |                                        |   *
!    [HV]                [HV]                                     [HV] *
!     |                   |                                        |   *
!     |                   |                                        |   *
!    4|                  5|                                       6|   *
!     O--------[VV]-------O------------------[VV]------------------O   *
!     |                   |                                        |   *
!     |                   |                                        |   *
!    [HV]                [HV]                                     [HV] *
!     |                   |                                        |   *
!    7|                  8|                                       9|   *
!     O--------[VV]-------O------------------[VV]------------------O   *
!                                                                      *
!     [HV] = Horizontal Velocity Given                                 *
!     [VV] = Vertical Velocity Given                                   *
!      O   = Node where Information is Needed for TecPlot Purposes     *
!                                                                      *
!***********************************************************************

!  This Subroutine Calls No Other Subroutines


      SUBROUTINE INTERFLOWF

      Use FISHY; USE GLOBAL; Use GEOMC; Use SCREENC
      IMPLICIT NONE
      
      REAL        UVERT,HVERT,WHORZ,DLXHO
      REAL        XV,FDDU,FDDW,LASTJDAY
      REAL        XRANGE,YRANGE
      REAL        TOPK,LFTI,RGTI
      REAL        CT,BT,AT,BD,CV,AD

      INTEGER     WB,NN,KK,II
      INTEGER     TOTNN,N,JJZ,J,K

      PARAMETER(N=3)                            ! The order of the Newton Interpolating Polynomial

      SAVE        LASTJDAY,TOPK

      DIMENSION   UVERT(N+1),HVERT(N+1),WHORZ(N+1),DLXHO(N+1),FDDU(N+1,N+1),FDDW(N+1,N+1)
                  
      NN = 0                                                       ! NN = Global Node # (UpperLeft Corner of cell(K,I))
      DO 123 WB=1,NWB                                              ! WB = Water Body # / NWB = Total # of Water Bodies
        DO 124 JB=BS(WB),BE(WB)                                  ! JB = Branch #
          DO 125 K=2,KMX                                           ! K = Layer #
            DO 126 I=US(JB),(DS(JB)+1)                         ! I = Segment #  SW 5/2015
              IF ((K-1.LE.KB(I)).OR.(K-1.LE.KB(I-1))) THEN         ! All Nodes at or above Water Body Bottom
                NN = NN + 1                                        ! PRE-CHECK CALCULATION
                KK = NDINFO(NN,3)                                  ! PRE-CHECK CALCULATION
                II = NDINFO(NN,4)                                  ! PRE-CHECK CALCULATION
                IF ((I.NE.II).OR.(K.NE.KK)) THEN                   ! CHECK
                  WRITE(DATADEBUGFN,*) 'ERROR: Nodes not matching up for Flowfield Interpolation'
                  write(DATADEBUGFN,*)'WB=',wb,' JB=', JB
                  write(DATADEBUGFN,*)'NN=',nn,' kk=',kk,' ii=',ii,' k=',k,' i=',i

                  STOP
                ELSE
                END IF

!LINEAR INTERPOLATION of the Flow Field (i.e. Velocities)
                IF (LINEAR) THEN              ! LINEAR INTERPOLATION
                  IF (K.LT.KTWB(WB)) THEN                           !Above Water Surface
                    FLOWFIELD(K,I,1) = 0.0                         !   Horiz. Vel.
                    FLOWFIELD(K,I,2) = 0.0                         !   Vert.  Vel.
                  ELSE IF (K.EQ.KTWB(WB)) THEN                      !At Water Surface
                    FLOWFIELD(K,I,1) = U(K,I-1)                    !   Horiz. Vel.
                    FLOWFIELD(K,I,2) = 0.0                         !   Vert.  Vel.
                  ELSE IF (K.EQ.KB(I)+1) THEN                      !Water Body Bottom
                    FLOWFIELD(K,I,1) = U(K-1,I-1)
                    FLOWFIELD(K,I,2) = 0.0
                  ELSE IF (I.EQ.US(JB)) THEN                     !Upstream-most Nodes
                    YRANGE = .5*H(K-1,WB)+.5*H(K,WB)
                    FLOWFIELD(K,I,1) = (U(K-1,I-1)*(YRANGE-&
                      .5*H(K-1,WB)) + U(K,I-1)*(YRANGE-.5*H(K,WB)))&
                      /YRANGE
                    FLOWFIELD(K,I,2) = W(K-1,I)
                  ELSE IF (I.EQ.DS(JB)+1) THEN                   !Downstream-most Nodes
                    YRANGE = .5*H(K-1,WB)+.5*H(K,WB)
                    FLOWFIELD(K,I,1) = (U(K-1,I-1)*(YRANGE-&
                      .5*H(K-1,WB)) + U(K,I-1)*(YRANGE-.5*H(K,WB)))&
                      /YRANGE
                    FLOWFIELD(K,I,2) = W(K-1,I-1)
                  ELSE                                             !Analogous to Node 5 above (i.e., in the middle)
                    YRANGE = .5*H(K-1,WB)+.5*H(K,WB)
                    FLOWFIELD(K,I,1) = (U(K-1,I-1)*(YRANGE-&
                      .5*H(K-1,WB)) + U(K,I-1)*(YRANGE-.5*H(K,WB)))&
                      /YRANGE
                    XRANGE = .5*DLX(I-1)+.5*DLX(I)
                    FLOWFIELD(K,I,2) = ((W(K-1,I-1)*(XRANGE-&
                      .5*DLX(I-1)) + W(K-1,I)*(XRANGE-.5*DLX(I)))&
                      /XRANGE)
                  END IF
                ELSE

!3rd ORDER NEWTON INTERPOLATING POLYNOMIAL used for interpolation of the Flow Field (i.e. Velocities)
!                                                 1
!                   O-----[VV]----O------[VV]------O-----[VV]----O------[VV]------O
!                   |             |                |             |                |
!                   |             |                |             |                |
!                  [HV]          [HV]             [HV]          [HV]             [HV]
!                   |             |                |             |                |
!                   |             |               2|             |                |
!                   O-----[VV]----O------[VV]------O-----[VV]----O------[VV]------O
!                   |             |                |             |                |
!                   |             |                |             |                |
!                  [HV]          [HV]             [HV]          [HV]             [HV]
!                   |             |                |             |                |
!                  1|            2|                |            3|               4|
!                   O-----[VV]----O------[VV]------X-----[VV]----O------[VV]------O
!                   |             |                |             |                |
!                   |             |                |             |                |
!                  [HV]          [HV]             [HV]  (K,I)   [HV]             [HV]
!                   |             |                |             |                |
!                   |             |               3|             |                |
!                   O-----[VV]----O------[VV]------O-----[VV]----O------[VV]------O
!                   |             |                |             |                |
!                   |             |                |             |                |
!                  [HV]          [HV]             [HV]          [HV]             [HV]
!                   |             |                |             |                |
!                   |             |               4|             |                |
!                   O-----[VV]----O------[VV]------O-----[VV]----O------[VV]------O
!
                  XV = 0.0                          ! Location (X ABOVE) of Interpolated Velocity (both U & W)
                  IF (K-1.LT.KTWB(WB)) THEN          !Both Top Nodes are above Water Surface
                    UVERT(1) = U(K,I-1)             ! This Preserves the Horizontal Velocity at the Water Surface
                    UVERT(2) = U(K,I-1)             ! This Preserves the Horizontal Velocity at the Water Surface
                    HVERT(1) = -H(K,WB)/10000          ! Arbitrary value
                    HVERT(2) = -H(K,WB)/10001          ! Arbitrary value
                  ELSE IF (K-2.LT.KTWB(WB)) THEN     !Only Top-most Node above Water Surface
                    UVERT(1) = U(K-1,I-1)           ! This Preserves the Horizontal Velocity at the Water Surface
                    UVERT(2) = U(K-1,I-1)
                    HVERT(1) = -H(K-1,WB)              ! Arbitrary value
                    HVERT(2) = -0.5*H(K-1,WB)
                  ELSE                              !Both Top Nodes exist
                    UVERT(1) = U(K-2,I-1)
                    UVERT(2) = U(K-1,I-1)
                    HVERT(1) = -H(K-1,WB)-0.5*H(K-2,WB)
                    HVERT(2) = -0.5*H(K-1,WB)
                  END IF
                  IF (K.GT.KB(I)) THEN              !Both Bottom Nodes are below Water Body Bottom
                    UVERT(3) = 0                    ! Helps satisfy no-slip condition on bottom (i.e., U=0)
                    UVERT(4) = 0                    ! Helps satisfy no-slip condition on bottom (i.e., U=0)
                    HVERT(3) = H(K-1,WB)/10001         ! Arbitrary value
                    HVERT(4) = H(K-1,WB)/10000         ! Arbitrary value
                  ELSE IF (K+1.GT.KB(I)) THEN       !Only Bottom-most Node below Water Body Bottom
                    UVERT(3) = U(K,I-1)
                    UVERT(4) = 0                    ! Helps satisfy no-slip condition on bottom (i.e., U=0)
                    HVERT(3) = 0.5*H(K,WB)
                    HVERT(4) = H(K,WB)                 ! Helps satisfy no-slip condition on bottom (i.e., U=0)
                  ELSE                              !Both Bottom Nodes exist
                    UVERT(3) = U(K,I-1)
                    UVERT(4) = U(K+1,I-1)
                    HVERT(3) = 0.5*H(K,WB)
                    HVERT(4) = H(K,WB)+0.5*H(K+1,WB)
                  END IF
                  IF (I-1.LT.US(JB)) THEN         !Both Left Nodes are Upstream of Water Body
                    WHORZ(1) = 0                    ! Helps satisfy no-slip condition at Upstream Boundary
                    WHORZ(2) = 0                    ! Helps satisfy no-slip condition at Upstream Boundary
                    DLXHO(1) = -DLX(I)/10000        ! Arbitrary value
                    DLXHO(2) = -DLX(I)/10001        ! Arbitrary value
                  ELSE IF (I-2.LT.US(JB)) THEN    !Only Left-most Node Upstream of Water Body
                    WHORZ(1) = 0                    ! Helps satisfy no-slip condition at Upstream Boundary
                    WHORZ(2) = W(K-1,I-1)
                    DLXHO(1) = -DLX(I-1)            ! Arbitrary value
                    DLXHO(2) = -0.5*DLX(I-1)
                  ELSE                              !Both Upstream Nodes exist
                    WHORZ(1) = W(K-1,I-2)
                    WHORZ(2) = W(K-1,I-1)
                    DLXHO(1) = -DLX(I-1)-0.5*DLX(I-2)
                    DLXHO(2) = -0.5*DLX(I-1)
                  END IF
                  IF (I.GT.DS(JB)) THEN           !Both Right Nodes are Downstream of Water Body 
                    WHORZ(3) = 0                    ! Helps satisfy no-slip condition at Downstream Boundary
                    WHORZ(4) = 0                    ! Helps satisfy no-slip condition at Downstream Boundary
                    DLXHO(3) = DLX(I-1)/10001       ! Arbitrary value
                    DLXHO(4) = DLX(I-1)/10000       ! Arbitrary value
                  ELSE IF (I+1.GT.DS(JB)) THEN    !Only Right-most Node Downstream of Water Body
                    WHORZ(3) = W(K-1,I)
                    WHORZ(4) = 0                    ! Helps satisfy no-slip condition at Downstream Boundary
                    DLXHO(3) = 0.5*DLX(I)
                    DLXHO(4) = DLX(I)               ! Arbitrary value
                  ELSE                              !Both Downstream Nodes exist
                    WHORZ(3) = W(K-1,I)
                    WHORZ(4) = W(K-1,I+1)
                    DLXHO(3) = 0.5*DLX(I)
                    DLXHO(4) = DLX(I)+0.5*DLX(I+1)
                  END IF
! Begin Newton Interpolating Polynomial Algorithm(Algorithm obtained from "NUMERICAL METHODS for Engineers -
!                                                 Authors: Steven Chapra & Raymond Canale")
                  DO JJZ=1,N+1
                    FDDU(JJZ,1) = UVERT(JJZ)
                    FDDW(JJZ,1) = WHORZ(JJZ)
                  END DO
                  DO J=2,N+1
                    DO JJZ=1,(N+1-J)+1
                      FDDU(JJZ,J)=(FDDU(JJZ+1,J-1)-FDDU(JJZ,J-1))/(HVERT(JJZ+J-1)-HVERT(JJZ))
                      FDDW(JJZ,J)=(FDDW(JJZ+1,J-1)-FDDW(JJZ,J-1))/(DLXHO(JJZ+J-1)-DLXHO(JJZ))
                    END DO
                  END DO
                  IF (K.LT.KTWB(WB)) THEN                           !Above Water Surface
                    FLOWFIELD(K,I,1) = 0.0                         ! Horiz. Vel.
                    FLOWFIELD(K,I,2) = 0.0                         ! Vert.  Vel. - Original Vert Vel from W2 are in millimeters
                  ELSE
                    FLOWFIELD(K,I,1) = FDDU(1,1)+FDDU(1,2)*&
                         (XV-HVERT(1))+FDDU(1,3)*(XV-HVERT(1))*&
                         (XV-HVERT(2))+FDDU(1,4)*(XV-HVERT(1))*&
                         (XV-HVERT(2))*(XV-HVERT(3))
                    FLOWFIELD(K,I,2) = (FDDW(1,1)+FDDW(1,2)*&
                         (XV-DLXHO(1))+FDDW(1,3)*(XV-DLXHO(1))*&
                         (XV-DLXHO(2))+FDDW(1,4)*(XV-DLXHO(1))*&
                         (XV-DLXHO(2))*(XV-DLXHO(3)))              ! Original Vertical Vel from W2 are in millimeters
                  END IF
! End Newton Interpolating Polynomial Algorithm
                END IF

!NULL FLOW FIELD (if so desired)
                BOTK = BOTKK
                IF (NULLFIELDFF) THEN                              ! Establish a Null Flow Field
                  IF (NIT.EQ.0) TOPK = REAL(KTWB(WB))                ! Null Field established for WB
                  LFTI = REAL(US(JB))
                  RGTI = REAL(DS(JB)+1)
                  IF (VARYHVEL.EQ.1) THEN                          !Vary Horizontal Velocity Vertically
                    IF (FFNULLLINEAR) THEN                         !  Use Linear Variation in Null Horizontal Velocity Values
                      FLOWFIELD(K,I,1) = (TOPHVEL-BOTHVEL)/(TOPK-BOTK)*(REAL(K)-TOPK)+ TOPHVEL
                    ELSE                                           !  Use Parabolic Variation in Null Horizontal Velocity Values
                      CT = TOPHVEL
                      BT = (MIDKHVEL-TOPHVEL)/(MIDKH-TOPK)
                      AT = ((BOTHVEL-MIDKHVEL)/(BOTK-MIDKH) - BT)/(BOTK-TOPK)
                      FLOWFIELD(K,I,1)=AT*(REAL(K))**2+BT*REAL(K)+CT
                    END IF
                  ELSE IF (VARYHVEL.EQ.2) THEN                     !Vary Horizontal Velocity Horizontally
                    IF (FFNULLLINEAR) THEN                         !  Use Linear Variation in Null Horizontal Velocity Values
                      FLOWFIELD(K,I,1) = (LFTHVEL-RGTHVEL)/&
                                         (LFTI-RGTI)*(REAL(I)-LFTI)&
                                         + LFTHVEL
                    ELSE                                           !  Use Parabolic Variation in Null Horizontal Velocity Values
                      CT = LFTHVEL
                      BT = (MIDIHVEL-LFTHVEL)/(MIDIH-LFTI)
                      AT = ((RGTHVEL-MIDIHVEL)/(RGTI-MIDIH) - BT)/&
                           (RGTI-LFTI)
                      FLOWFIELD(K,I,1)=AT*(REAL(I))**2+BT*REAL(I)+CT
                    END IF
                  ELSE
                  END IF
                  IF (VARYVVEL.EQ.1) THEN                          !Vary Vertical Velocity Values Vertically
                    IF (FFNULLLINEAR) THEN                         !  Use Linear Variation in Null Vertical Velocity Values
                      FLOWFIELD(K,I,2) = (TOPVVEL-BOTVVEL)/&
                                         (TOPK-BOTK)*(REAL(K)-TOPK)&
                                         + TOPVVEL
                    ELSE                                           !  Use Parabolic Variation in Null Vertical Velocity Values
                      CV = TOPVVEL
                      BD = (MIDKVVEL-TOPVVEL)/(MIDKV-TOPK)
                      AD = ((BOTVVEL-MIDKVVEL)/(BOTK-MIDKV) - BD)/&
                           (BOTK-TOPK)
                      FLOWFIELD(K,I,2)=AD*(REAL(K))**2+BD*REAL(K)+CV
                    END IF
                  ELSE IF (VARYVVEL.EQ.2) THEN                     !Vary Vertical Velocity Values Horizontally
                    IF (FFNULLLINEAR) THEN                         !  Use Linear Variation in Null Vertical Velocity Values
                      FLOWFIELD(K,I,2) = (LFTVVEL-RGTVVEL)/&
                                         (LFTI-RGTI)*(REAL(I)-LFTI)&
                                         + LFTVVEL
                    ELSE                                           !  Use Parabolic Variation in Null Vertical Velocity Values
                      CV = LFTVVEL
                      BD = (MIDIVVEL-LFTVVEL)/(MIDIV-LFTI)
                      AD = ((RGTVVEL-MIDIVVEL)/(RGTI-MIDIV) - BD)/&
                           (RGTI-LFTI)
                      FLOWFIELD(K,I,2)=AD*(REAL(I))**2+BD*REAL(I)+CV
                    END IF
                  ELSE
                  END IF
                ELSE
                END IF

! Begin Calculation of Accelerations ==> Time Derivative of Velocity (Backward Differencing)
                IF (NIT.NE.0) THEN
                  IF ((JDAY-LASTJDAY).EQ.0) THEN
                    FLOWFIELD(K,I,3) = 0.0
                    FLOWFIELD(K,I,4) = 0.0
                  ELSE
                    FLOWFIELD(K,I,3) = (FLOWFIELD(K,I,1)-&
                         LASTFIELD(K,I,1))/((JDAY-LASTJDAY)*24*3600)
                    FLOWFIELD(K,I,4) = (FLOWFIELD(K,I,2)-&
                         LASTFIELD(K,I,2))/((JDAY-LASTJDAY)*24*3600)
                  END IF
                  LASTFIELD(K,I,1) = FLOWFIELD(K,I,1)
                  LASTFIELD(K,I,2) = FLOWFIELD(K,I,2)
                  LASTJDAY = JDAY
                ELSE                              ! No Time Derivative if Time=0 (i.e., NIT=0)
                  FLOWFIELD(K,I,3) = 0
                  FLOWFIELD(K,I,4) = 0
                  LASTFIELD(K,I,1) = FLOWFIELD(K,I,1)
                  LASTFIELD(K,I,2) = FLOWFIELD(K,I,2)
                  LASTJDAY = JDAY
                END IF
! End Calculation of Accelerations

              ELSE
              END IF
  126       CONTINUE
  125     CONTINUE
  124   CONTINUE
  123 CONTINUE

      RETURN
      END


!***********************************************************************
!***********************************************************************
!*                                                                    **
!*         P R E P A R E   F I S H   F O R   D I S P L A Y            **
!*                                                                    **
!***********************************************************************
!***********************************************************************
! Number of TecPlot Zones       = # of timesteps            (Animation by Zone)
!           TecPlot Limit     ==> 32,700 per Data Set
! Number of TecPlot Data Points = # of nodes                (for Grid Display)
!                                 # of fish                 (for Fish Display)
!           TecPlot Limit     ==> 2 billion per Variable
! Number of TecPlot Variables   = # of flow & const. param  (for Grid Display)
!                                 # of fish parameters      (for Fish Display)
!           TecPlot Limit     ==> 32,700 per Data Set
! Number of TecPlot Data Sets   = 2*(# of branches)         (for Grid & Fish Displays)
!           TecPlot Limit     ==> 128

!  This Subroutine Calls No Other Subroutines


      SUBROUTINE FISHPLOT                ! This Subroutine preps fish information for TecPlot

      
      Use FISHY; Use GEOMC; USE GLOBAL; Use GDAYC; Use SCREENC
      
      IMPLICIT NONE
            
      REAL        XLOC,ZLOC,YLOC,FSZ,FAG,UF,VF,WF
      REAL        NETXPUT,NETZTOP,NETZBOT,SNDXUPS
      REAL        SNDXDWN,SNDZTOP,WBBOTTOM

      INTEGER     WB,FISHNBP
      INTEGER     BRF,FBRCH,WATER,K
      INTEGER     HAIMP,NETIMP

      CHARACTER    FTYP1*1,FTYP2*2
      CHARACTER*4  FILE_PREFIX
      CHARACTER*4  FILE_SUFFIX
      CHARACTER*10 FILENAME
      DATA         FILE_PREFIX/ 'Fish'/
      DATA         FILE_SUFFIX/ '.dat'/

      IF (NIT.EQ.0) NZONES = 0                        ! This # is inputted into a TecPlot Macro for animation
      NZONES = NZONES + 1                             ! This # is inputted into a TecPlot Macro for animation
      DO 132 WB=1,NWB                                 ! WB = Water Body # / NWB = Total # of Water Bodies
        DO 133 JB=BS(WB),BE(WB)                     ! JB = Branch #
          NBRF(JB) = 0                              ! NBRF(JB) = # of Fish in Branch JB
          DO 134 FN=1,NFISH                           ! NFISH = Total # of Fish in System
            IF (INT(FISHES(FN,12)).EQ.1) GOTO 134     ! If Fish left system do not output its info to TecPlot
            FISHNBP = INT(FISHES(FN,6))               ! FISHNBP = Branch where current fish resides
            IF (FISHNBP.EQ.JB) THEN
              NBRF(JB) = NBRF(JB) + 1                     ! NBRF    = # of Fish in Branch JB
              BRCHFISH(JB,NBRF(JB),1) =  FISHES(FN,1)     ! FIMP    = Segment IMP where fish is located
              BRCHFISH(JB,NBRF(JB),2) =  FISHES(FN,2)     ! FXLOC   = Fish location within segment IMP
              BRCHFISH(JB,NBRF(JB),3) =  FISHES(FN,3)     ! FKMP    = Layer KMP where fish is located
              BRCHFISH(JB,NBRF(JB),4) =  FISHES(FN,4)     ! FZLOC   = Fish location within layer KMP
              BRCHFISH(JB,NBRF(JB),5) =  FISHES(FN,5)     ! FYLOC   = Lateral fish location
              BRCHFISH(JB,NBRF(JB),6) =  FISHES(FN,6)     ! FNBP    = Branch where fish is located
              BRCHFISH(JB,NBRF(JB),7) =  FISHES(FN,7)     ! FSIZE   = Size (i.e., length) of fish in meters
              BRCHFISH(JB,NBRF(JB),8) =  FISHES(FN,8)     ! FAGE    = Age of the fish
              BRCHFISH(JB,NBRF(JB),9) =  FISHES(FN,9)     ! UFISH   = Longitudinal velocity of the fish
              BRCHFISH(JB,NBRF(JB),10) = FISHES(FN,10)    ! VFISH   = Lateral velocity of the fish
              BRCHFISH(JB,NBRF(JB),11) = FISHES(FN,11)    ! WFISH   = Vertical velocity of the fish
            ELSE
            if(JB.eq.1)then
            write(DATADEBUGFN,*)'Lost particle JB=1:JDAY,FN,fimp,fkmp,fishes(fn,6),fishnbp:',jday,fn,fimp,fkmp,fishes(fn,6),fishnbp  ! DEBUG SW 
            end if
            END IF
  134     CONTINUE
          IF (JB.LT.10) THEN                        ! Creating the output file names based on Branch #
            WRITE(FTYP1,9120) JB
            FILENAME=FILE_PREFIX//FTYP1//FILE_SUFFIX
 9120       FORMAT(I1)
          ELSE IF (JB.LE.64) THEN                   ! TecPlot is limited to 128 Data Sets (i.e. Data Files)
            WRITE(FTYP2,9130) JB                    !   That leaves: 64 Data Sets for Grid Display and
            FILENAME=FILE_PREFIX//FTYP2//FILE_SUFFIX  !                64 Data Sets for Fish Display
 9130       FORMAT(I2)
          ELSE                                        ! CHECK
            WRITE(DATADEBUGFN,*) 'ERROR: Exceeded Maximum Number of TecPlot Data Sets ==># Data Sets = # of Branches'
            STOP
          END IF
          IF (NIT.EQ.0) THEN
            OPEN(30000+JB,FILE=FILENAME,STATUS='UNKNOWN')             ! Opening the Output Files Based on Branch #
            WRITE(30000+JB,9300) JB                                 ! Creating Output File Header for TecPlot
 9300       FORMAT('TITLE = "FISH INFO for Branch',I6,'"')
            WRITE(30000+JB,9310)                                      ! Creating Output File Header for TecPlot
 9310       FORMAT('VARIABLES = "X", "Z", "FSIZE", "FAGE", "WATER"')
          ELSE
            WRITE(30000+JB,*) ' '
          END IF
          WRITE(30000+JB,9320) JDAY,(1+NBRF(JB))                    ! Creating Output File Header for TecPlot
 9320     FORMAT('ZONE T="JDAY ',F9.2,'", I=',I7,', F=POINT')
          WRITE(30000+JB,9330) 0.01,0.01,0.0,0.0,-1         ! This is a dumby fish so TecPlot can at least
                                                              !   open the file if no fish are in JB
          DO 135 BRF=1,NBRF(JB)                             ! NBRF = # of Fish in Branch JB
            K = INT(BRCHFISH(JB,BRF,3))                     ! K = Layer # where Fish is located
            I = INT(BRCHFISH(JB,BRF,1))                     ! I = Segment # where Fish is located
            XLOC = NODES(K,I,1) + BRCHFISH(JB,BRF,2)        ! XLOC = X Distance (within Branch) where Fish is located
            IF (K.LT.KTWB(WB)) THEN
         !     ZLOC = NODES(KTWB(WB),I,2) - 0                   ! ZLOC(Adjusted): Due to inability of TecPlot to show
            zloc=elws(i)   ! SW 5/2015 CHANGE TO ELEVATION RATHER THAN DEPTH
            ELSE                                              !                 actual Water Surface
         !     ZLOC = NODES(K,I,2) - BRCHFISH(JB,BRF,4)      ! ZLOC = Z Distance where Fish is located
            zloc=el(k,i)- BRCHFISH(JB,BRF,4)   ! SW 5/2015
            END IF
            YLOC = BRCHFISH(JB,BRF,5)                       ! YLOC = Y Distance where Fish is located
            FBRCH = INT(BRCHFISH(JB,BRF,6))                 ! FBRCH = Branch # where Fish is located
            FSZ = BRCHFISH(JB,BRF,7)                        ! FSZ = Fish Size
            FAG = BRCHFISH(JB,BRF,8)                        ! FAG = Fish Age
            WATER = 1                                         ! WATER = 1 (for value-blanking purposes in TecPlot)
            UF = BRCHFISH(JB,BRF,9)                         ! UF = Fish Horizontal Velocity (m/s)
            VF = BRCHFISH(JB,BRF,10)                        ! VF = Fish Lateral Velocity (m/s)
            WF = BRCHFISH(JB,BRF,11)                        ! WF = Fish Vertical Velocity (m/s)
            IF (WBSKIP.AND.(WB.NE.WBRUN)) GOTO 135            ! Option:Skip output for all Water Bodies, except WBRUN
            WRITE(30000+JB,9330) XLOC,ZLOC,FSZ,FAG,WATER
 9330       FORMAT(F15.2,F10.2,F9.2,F9.2,I4)
  135     CONTINUE
          WRITE(30000+JB,9440) MONTH,GDAY,YEAR,INT(HR),AMPM,&
                                 JDAY,NZONES                  ! Output info used as dynamic text in TecPlot Animation
 9440     FORMAT('TEXT X=47, Y=90, F=HELV-BOLD, HU=FRAME, AN=MIDCENTER,&
       C=RED, H=3, T="',A10,I3,',',I5,4X,I3,A2,'    (JDAY',F8.3,')", ZN=',I6)
          WRITE(30000+JB,9518) NBRF(JB),NZONES            ! Output info used as dynamic text in TecPlot Animation
 9518     FORMAT('TEXT X=30.0, Y=21.0, F=HELV-BOLD, HU=FRAME,&
      AN=MIDRIGHT, C=BLACK, H=2.1, T="',I6,'", ZN=',I6)


!  Display Virtual Gillnets

          IF (GILLNETSAMPLING) THEN
            DO 15 NET=1,NUMGILLNETS                                     ! NUMGILLNETS = Total # of Gillnets
              IF ((GILLNETS(NET,11).EQ.1).AND. &                         ! Is Gillnet Active? [0=No , 1=Yes]
                  (GILLNETS(NET,12).EQ.JB)) THEN                      ! Is Gillnet Located in Branch "JB"?
                NETXPUT = GILLNETS(NET,5)
                NETZTOP = GILLNETS(NET,8)*-1 +  &                        ! Depth Below Water Surface (m) to Top of Gillnet
                          NODES(KTWBF,CUS(GILLNETS(NET,12)),2)           !   Any Segment I in Branch JB would suffice for NODES(K,I,2)
                NETZBOT = GILLNETS(NET,9)*-1 +   &                       ! Depth Below Water Surface (m) to Bottom of Gillnet
                          NODES(KTWBF,CUS(GILLNETS(NET,12)),2)           !   Any Segment I in Branch JB would suffice for NODES(K,I,2)

!    Need to Determine Where the Water Body Bottom is for TecPlot Display of Gillnet
                NETIMP = US(GILLNETS(NET,12))
   26           IF ((NETXPUT.GE.NODES(KTWBF,NETIMP,1)).AND.  &            ! Any Layer K in Branch JB would suffice for NODES(K,I,1)
                    (NETXPUT.LT.NODES(KTWBF,NETIMP+1,1))) THEN
                  WBBOTTOM = NODES(KB(NETIMP)+1,NETIMP,2)               ! Depth (m) of Water Body Bottom in Segment NETIMP
                ELSE
                  NETIMP = NETIMP + 1
                  IF (NETIMP.GT.DS(GILLNETS(NET,12))+1) THEN            ! CHECK TO SEE IF Gillnet was placed OUT-OF-BOUNDS
                    WRITE(*,9578) NET
 9578               FORMAT('Gillnet #',I4,' exceeded Upstream     or&
      Downstream Boundaries of the Branch')
                    GOTO 27
                  END IF
                  GOTO 26
                END IF
                IF (NETZBOT.LT.WBBOTTOM) NETZBOT = WBBOTTOM
!    Finished Determining Where the Bottom of the Gillnet is

                WRITE(30000+JB,9509) NZONES
 9509           FORMAT('GEOMETRY X=0, Y=0, CS=GRID, C=BLACK,&
      L=DASHDOTDOT, PL=0.3, LT=1.0, T=LINE, F=POINT, ZN=',I6,&
      ', MFC="SCOPE = LOCAL"')
                WRITE(30000+JB,9557) NETXPUT,NETZTOP,NETXPUT,NETZBOT
 9557           FORMAT('1',/,'2',/,F15.1,'  ',F10.2,/,F15.1,'  ',F10.2)
              END IF
   27         CONTINUE
   15       CONTINUE
          END IF


!  Display Virtual Hydroacoustics Sampling 'Box'

          IF (ACOUSTICSAMPLING) THEN
            DO 22 SND=1,NUMACOUSTICS                                    ! NUMACOUSTICS = Total # of Hydroacoustic Surveys
              IF (HAOPERAT(SND,13).EQ.1) THEN                           ! Is Hydroacoustic Survey Active? [0=No , 1=Yes]
                SNDXUPS = HAOPERAT(SND,5)                               ! X-Location (m) of Upstream End of HA Sampling 'Box'
                SNDXDWN = HAOPERAT(SND,6)                               ! X-Location (m) of Downstream End of HA Sampling 'Box'
                SNDZTOP = HAOPERAT(SND,10)*-1 +   &                         ! Depth Below Water Surface (m) to Top of HA Sampling 'Box'
                          NODES(KTWBF,CUS(HAOPERAT(SND,14)),2)           !   Any Segment I in Branch JB would suffice for NODES(K,I,2)

!    Need to Determine Where the Water Body Bottom is for TecPlot Display of HA Beam
                HAIMP = US(HAOPERAT(SND,14))
   24           IF ((SNDXUPS.GE.NODES(KTWBF,HAIMP,1)).AND. &              ! Any Layer K in Branch JB would suffice for NODES(K,I,1)
                    (SNDXUPS.LT.NODES(KTWBF,HAIMP+1,1))) THEN
                  WBBOTTOM = NODES(KB(HAIMP)+1,HAIMP,2)                 ! Depth (m) of Water Body Bottom in Segment HAIMP
                ELSE
                  HAIMP = HAIMP + 1
                  IF (HAIMP.GT.DS(HAOPERAT(SND,14))+1) THEN             ! CHECK TO SEE IF HA Survey WENT OUT-OF-BOUNDS
                    WRITE(*,9577) SND
 9577               FORMAT('HA Survey #',I4,' exceeded Upstream     or&
      Downstream Boundaries of the Branch')
                    GOTO 25
                  END IF
                  GOTO 24
                END IF
!    Finished Determining Where the Bottom of the HA Beam is

                WRITE(30000+JB,9567) SNDXUPS,WBBOTTOM,NZONES
 9567           FORMAT('GEOMETRY X=',F15.1,', Y=',F10.2,', CS=GRID,&
      C=BLACK, L=SOLID, LT=1.0, T=RECTANGLE, ZN=',I6,&
      ', MFC="SCOPE = LOCAL"')
                WRITE(30000+JB,9566) (SNDXDWN-SNDXUPS),&
                                       (SNDZTOP-WBBOTTOM)
 9566           FORMAT(F15.1,'  ',F10.2)
              END IF
   25         CONTINUE
   22       CONTINUE
          END IF


!  The Following is Used for Creating and Positioning Static Text in the TecPlot Animation

          IF (NIT.EQ.0) THEN
            WRITE(30000+JB,9541) 
 9541       FORMAT('TEXT X=25.0, Y=25.0, F=HELV-BOLD, HU=FRAME,&
      AN=MIDCENTER, C=BLACK, H=2.1, T="# of Fish/Particles"')
!            WRITE(30000+JB,9519) TXTSTIMRULE
! 9519       FORMAT('TEXT X=11.2, Y=47.2, F=COURIER, HU=FRAME,&
!      AN=MIDLEFT, C=BLACK, H=2.1, T="Stimuli Rules  :',A4,'"')
!            WRITE(30000+JB,9520) TXTMULTIRULE
! 9520       FORMAT('TEXT X=11.2, Y=44.9, F=COURIER, HU=FRAME,&
!      AN=MIDLEFT, C=BLACK, H=2.1, T="(Multi-Response)',A4,'"')
!            WRITE(30000+JB,9521) TXTVELORULE
! 9521       FORMAT('TEXT X=11.2, Y=42.6, F=COURIER, HU=FRAME,&
!      AN=MIDLEFT, C=BLACK, H=2.1, T="- Velocity      ',A4,'"')
!            WRITE(30000+JB,9522) TXTTEMPRULE
! 9522       FORMAT('TEXT X=11.2, Y=40.3, F=COURIER, HU=FRAME,&
!      AN=MIDLEFT, C=BLACK, H=2.1, T="- Temperature   ',A4,'"')
!            WRITE(30000+JB,9523) TXTDORULE
! 9523       FORMAT('TEXT X=11.2, Y=38.0, F=COURIER, HU=FRAME,&
!      AN=MIDLEFT, C=BLACK, H=2.1, T="- Diss Oxygen   ',A4,'"')
!            WRITE(30000+JB,9524) TXTRANDRULE
! 9524       FORMAT('TEXT X=11.2, Y=35.7, F=COURIER, HU=FRAME,&
!      AN=MIDLEFT, C=BLACK, H=2.1, T="- Randomization ',A4,'"')
!            WRITE(30000+JB,9525) TXTNULLWQ
! 9525       FORMAT('TEXT X=31.8, Y=42.6, F=COURIER, HU=FRAME,&
!      AN=MIDLEFT, C=BLACK, H=2.1, T="WQ  Null Field:',A4,'"')
!            WRITE(30000+JB,9526) TXTNULLFF
! 9526       FORMAT('TEXT X=31.8, Y=40.3, F=COURIER, HU=FRAME,&
!      AN=MIDLEFT, C=BLACK, H=2.1, T="HYD Null Field:',A4,'"')
!            WRITE(30000+JB,9527) TXTPASSTRAN
! 9527       FORMAT('TEXT X=31.8, Y=38.0, F=COURIER, HU=FRAME,&
!      AN=MIDLEFT, C=BLACK, H=2.1, T="Passive Transp:',A4,'"')
!            WRITE(30000+JB,9528) TXTSCHOOLG
! 9528       FORMAT('TEXT X=31.8, Y=35.7, F=COURIER, HU=FRAME,&
!      AN=MIDLEFT, C=BLACK, H=2.1, T="Schooling Behv:',A4,'"')
!already commented out            WRITE(30000+JB,9540) INFLSPAN     ! This variable is not defined ! SW 1/9/01
!      WRITE(30000+JB,9540)               ! SW 1/9/01
! 9540       FORMAT('TEXT X=38, Y=31, F=HELV-BOLD, HU=FRAME,&
!      AN=MIDCENTER, C=BLACK, H=2.2, T="% Influence for Timestep"') 
          END IF
  133   CONTINUE
  132 CONTINUE

      RETURN
      END


!***********************************************************************
!***********************************************************************
!*                                                                    **
!*     B I L I N E A R   S P L I N E   I N T E R P O L A T I O N      **
!*                                                                    **
!***********************************************************************
!***********************************************************************
!                                                                      *
!   Locations:1,2,3,4,5            4 (Above Fish)                      *
!    (LOCATE)                      ^                                   *
!                               UZSENSPH                               *
!                                  |                                   *
!                                _____Fish                             *
!  (Behind Fish)             |\_/     \_     (Forward of Fish)         *
!       3  <--- BXSENSPH --- | _   5   _> ----- FXSENSPH ----->  1     *
!                            |/ \_____/                                *
!                                  |                                   *
!                               DZSENSPH                               *
!                                  V                                   *
!                                  2 (Below Fish)                      *
!                                                                      *
!     Variables are relative to the fish, and the NFS compensates      *
!                 if fish x-axis orientation changes                   *
!                                                                      *
!***********************************************************************

!  This Subroutine Calls No Other Subroutines


      SUBROUTINE SPLINE             ! Bilinear Spline Algorithm obtained from
                                    !   "TWO DIMENSIONAL SPLINE INTERPOLATION
      Use FISHY
      USE GLOBAL                             !    ALGORITHMS - Author: Helmuth Spath"
      Use GEOMC
      
      IMPLICIT NONE
     
      REAL      XH1,XH2,H3,H4,DV,DW,DX,DY

 !Check for Boundary Violations: Horizontal Plane

   95 IF (XLOK.GT.DLX(ILOK)) THEN                             ! Must look into the next segment downstream
        IF ((ILOK+1).GT.DS(FNBP)) THEN                        ! Already at most downstream segment in branch
          XLOK = DLX(ILOK)
        ELSE IF (KLOK.GT.KB(ILOK+1)) THEN
          XLOK = DLX(ILOK)
        ELSE                                                  ! More segments exist downstream of fish location
          XLOK = XLOK - DLX(ILOK)
          ILOK = ILOK + 1                                     ! Looking into next segment downstream
        END IF
        IF (XLOK.GT.DLX(ILOK)) GOTO 95
      ELSE IF (XLOK.LT.0) THEN                                ! Must look into the next segment upstream
        IF (ILOK.EQ.CUS(FNBP)) THEN                           ! Already at most upstream segment in branch
          XLOK = 0.0
        ELSE IF (KLOK.GT.KB(ILOK-1)) THEN
          XLOK = 0.0
        ELSE                                                  ! More segments exist upstream of current segment
          ILOK = ILOK - 1                                     ! Looking into next segment upstream
          XLOK = DLX(ILOK) + XLOK
        END IF
        IF (XLOK.LT.0) GOTO 95
      ELSE
      END IF

!Check for Boundary Violations: Vertical Plane

   96 IF (ZLOK.GT.H(KLOK,FJR)) THEN                               ! Must look into the next layer below
        IF (KLOK.EQ.KB(ILOK)) THEN                            ! Already at water body bottom
          ZLOK = H(KLOK,FJR)
        ELSE                                                  ! More layers exist below fish location
          ZLOK = ZLOK - H(KLOK,FJR)
          KLOK = KLOK + 1                                     ! Looking into next layer below
        END IF
        IF (ZLOK.GT.H(KLOK,FJR)) GOTO 96
      ELSE IF (ZLOK.LT.0) THEN                                ! Must look into the next layer above
        IF (KLOK.EQ.KTWBF.or.klok.eq.1) THEN       ! SW 1/17/01        ! Already at water surface layer
          ZLOK = 0.0
        ELSE                                                  ! More layers exist above fish location
          KLOK = KLOK - 1                                     ! Looking into next layer above
          ZLOK = H(KLOK,FJR) + ZLOK
        END IF
        IF (ZLOK.LT.0) GOTO 96
      ELSE
      END IF

!Begin Bilinear Spline Interpolation Algorithm

      DV = XLOK-0.0                       ! XLOK = X-coord of interpolated location; 0.0 = X-coord of left nodes
      DW = ZLOK-0.0                       ! ZLOK = Y-coord of interpolated location; 0.0 = Y-coord of upper nodes
      DX = DLX(ILOK)                      ! DX = X distance between left nodes and right nodes
      DY = H(KLOK,FJR)                        ! DY = Y distance between upper nodes and lower nodes

      IF (VARIABLE.EQ.1) THEN             ! Horizontal Velocity Calculations
        XH1 = FLOWFIELD(KLOK,ILOK,1)       ! UpperLeft Node
        XH2 = FLOWFIELD(KLOK,ILOK+1,1)     ! UpperRight Node
        H3 = FLOWFIELD(KLOK+1,ILOK,1)     ! LowerLeft Node
        H4 = FLOWFIELD(KLOK+1,ILOK+1,1)   ! LowerRight Node
      ELSE IF (VARIABLE.EQ.2) THEN        ! Vertical Velocity Calculations
        XH1 = FLOWFIELD(KLOK,ILOK,2)       ! XH1,XH2,H3,H4 = Value of Vertical Velocity at 4 cell nodes
        XH2 = FLOWFIELD(KLOK,ILOK+1,2)
        H3 = FLOWFIELD(KLOK+1,ILOK,2)
        H4 = FLOWFIELD(KLOK+1,ILOK+1,2)
      ELSE IF (VARIABLE.EQ.3) THEN        ! Temperature Calculations
        XH1 = WQFIELD(KLOK,ILOK,1)         ! XH1,XH2,H3,H4 = Value of Temperature at 4 cell nodes
        XH2 = WQFIELD(KLOK,ILOK+1,1)
        H3 = WQFIELD(KLOK+1,ILOK,1)
        H4 = WQFIELD(KLOK+1,ILOK+1,1)
      ELSE IF (VARIABLE.EQ.4) THEN        ! Dissolved Oxygen Calculations
        XH1 = WQFIELD(KLOK,ILOK,2)         ! XH1,XH2,H3,H4 = Value of Dissolved Oxygen at 4 cell nodes
        XH2 = WQFIELD(KLOK,ILOK+1,2)
        H3 = WQFIELD(KLOK+1,ILOK,2)
        H4 = WQFIELD(KLOK+1,ILOK+1,2)
      ELSE
      END IF

      VALUE = XH1+(XH2-XH1)/DX*DV+(H3-XH1)/DY*DW+(H4-H3-XH2+XH1)/(DX*DY)*DV*DW

!End Bilinear Spline Interpolation Algorithm

      RETURN
      END


!***********************************************************************
!***********************************************************************
!*                                                                    **
!*       V I R T U A L   S A M P L I N G   -   G I L L N E T S        **
!*                                                                    **
!***********************************************************************
!***********************************************************************
!  This Subroutine Calls the Following Subroutines:
!      SPLINE

      SUBROUTINE VGILLNETS      
      USE GLOBAL; Use FISHY; USE GEOMC
      IMPLICIT NONE

      DOUBLE PRECISION TEMPVALUE,DOVALUE

      REAL            OLDXDIST,OLDYDISTFROMCENTER,OLDZDIST
      REAL            NEWXDIST,NEWYDISTFROMCENTER,NEWZDIST
      REAL            NETXPUT,NETYLFT,NETYRGT,NETZTOP,NETZBOT,FRACTIONX
      REAL            YDISTATNET,ZDISTATNET,XLOCATNET,ZLOCATNET

      INTEGER         NNUM

      OLDXDIST           = NODES(FKMP,FIMP,1) + OLDFXLOC              ! Previous Distance Downstream from Dam (m) to Location of Fish
      OLDYDISTFROMCENTER = OLDFYLOC - B(FKMP,FIMP)/2                  ! Previous Width Distance (m) from Lake Center
      OLDZDIST           = (NODES(FKMP,FIMP,2)-NODES(KTWBF,FIMP,2))&
                            - OLDFZLOC                                ! Previous Depth (m) of Fish
      NEWXDIST           = NODES(FKMP,FIMP,1) + FXLOC                 ! New Distance Downstream from Dam (m) to Location of Fish
      NEWYDISTFROMCENTER = FYLOC - B(FKMP,FIMP)/2                     ! New Width Distance (m) from Lake Center
      NEWZDIST           = (NODES(FKMP,FIMP,2)-NODES(KTWBF,FIMP,2))&
                            - FZLOC                                   ! New Depth (m) of Fish

      DO 99 NNUM=1,NUMGILLNETS                                        ! NUMGILLNETS = Total Number of Gillnets
        IF ((GILLNETS(NNUM,11).EQ.0).AND.(FSNAG.EQ.NNUM)) THEN        ! Release all Fish caught in a gillnet no longer active
          FSNAG = 0
        ELSE IF ((GILLNETS(NNUM,11).EQ.1).AND.(FSNAG.EQ.NNUM)) THEN   ! All Fish caught in a gillnet do not move
          FXLOC = OLDFXLOC
          FYLOC = OLDFYLOC
          FZLOC = OLDFZLOC
        ELSE IF ((GILLNETS(NNUM,11).EQ.1).AND.(FSNAG.EQ.0)) THEN      ! Test to see if the free Fish will get snagged in a gillnet
          NETXPUT = GILLNETS(NNUM,5)                                  ! X-Location (m) Downstream Where Gillnet is Placed
          NETYLFT = GILLNETS(NNUM,6)                                  ! Meters from Center of Lake Gillnet will extend to Left Bank
          NETYRGT = GILLNETS(NNUM,7)                                  ! Meters from Center of Lake Gillnet will extend to Right Bank
          NETZTOP = GILLNETS(NNUM,8)*-1                               ! Depth Below Water Surface (m) to Top of Gillnet: (+) Number
          NETZBOT = GILLNETS(NNUM,9)*-1                               ! Depth Below Water Surface (m) to Bottom of Gillnet: (+) Number
          FRACTIONX = 0
          IF ((OLDXDIST.LE.NETXPUT).AND.(NEWXDIST.GE.NETXPUT)) THEN   ! Fish May be Caught in Gillnet Approaching from Upstream
            FRACTIONX = (NETXPUT-OLDXDIST)/(NEWXDIST-OLDXDIST)
          ELSE IF ((OLDXDIST.GE.NETXPUT).AND.&
                   (NEWXDIST.LE.NETXPUT)) THEN                        ! Fish May be Caught in Gillnet Approaching from Downstream
            FRACTIONX = (OLDXDIST-NETXPUT)/(OLDXDIST-NEWXDIST)
          END IF
          IF (FRACTIONX.NE.0) THEN                                    ! Fish May be Caught in Gillnet
            YDISTATNET =  OLDYDISTFROMCENTER +&
                         (NEWYDISTFROMCENTER-OLDYDISTFROMCENTER)&
                         *FRACTIONX                                   ! Location of Fish in Width Dimension when it crosses Gillnet
            ZDISTATNET = OLDZDIST + (NEWZDIST-OLDZDIST)*FRACTIONX     ! Depth of Fish when it crosses Gillnet
            IF (((YDISTATNET.GE.NETYLFT).AND.&
                 (YDISTATNET.LE.NETYRGT)).AND.&
                ((ZDISTATNET.LE.NETZTOP).AND.&
                 (ZDISTATNET.GE.NETZBOT))) THEN                       ! Fish has just been caught in Gillnet
              FSNAG = NNUM                                            ! Fish caught in Gillnet # NNUM
              SNAGCOUNT(NNUM) = SNAGCOUNT(NNUM) + 1                   ! Tally # of Fish caught in Gillnet # NNUM
              FXLOC = NETXPUT - NODES(FKMP,FIMP,1)                    ! Location of Fish in Gillnet
              FYLOC = YDISTATNET + B(FKMP,FIMP)/2                     ! Location of Fish in Gillnet
              FZLOC = ZDISTATNET -&
                     (NODES(FKMP,FIMP,2)-NODES(KTWBF,FIMP,2))          ! Location of Fish in Gillnet
              XLOCATNET =      NETXPUT   - NODES(FKMP,FIMP,1)
              ZLOCATNET = -1*(ZDISTATNET - NODES(FKMP,FIMP,2))
              ILOK = FIMP                                             ! Segment I where Variable Value is desired
              XLOK = XLOCATNET                                        ! Location within Segment I where Variable Value is desired
              KLOK = FKMP                                             ! Layer K where Variable Value is desired
              ZLOK = ZLOCATNET                                        ! Location within Layer K where Variable Value is desired
              NETCATCH(NNUM,SNAGCOUNT(NNUM),1) = SNAGCOUNT(NNUM)      ! Tally # of Fish Caught
              NETCATCH(NNUM,SNAGCOUNT(NNUM),2) = NETXPUT              ! X-Location (m) of Gillnet Downstream from RBR Dam
              NETCATCH(NNUM,SNAGCOUNT(NNUM),3) = ZDISTATNET           ! Depth (m) at Which the Fish was Caught
              Variable=3 ! SW 1/15/01
              CALL SPLINE                                             ! VARIABLE: 1=Horz Vel, 2=Vert Vel, 3=Temp, 4=Diss Oxy
              NETCATCH(NNUM,SNAGCOUNT(NNUM),4) = VALUE            ! Temperature Where the Fish was Caught
              Variable=4 ! SW 1/15/01
              CALL SPLINE                                             ! VARIABLE: 1=Horz Vel, 2=Vert Vel, 3=Temp, 4=Diss Oxy
              NETCATCH(NNUM,SNAGCOUNT(NNUM),5) = VALUE              ! Dissolved Oxygen Where the Fish was Caught
            END IF
          END IF
        END IF
   99 CONTINUE

      RETURN
      END


!***********************************************************************
!***********************************************************************
!*                                                                    **
!*  V I R T U A L   S A M P L I N G  -  H Y D R O A C O U S T I C S   **
!*                                                                    **
!***********************************************************************
!***********************************************************************
!  This Subroutine Calls No Other Subroutines


      SUBROUTINE ACOUSTICS
    Use FISHY; USE GEOMC
    IMPLICIT NONE
    
      REAL            OLDXDIST,OLDYDISTFROMCENTER,OLDZDIST,NEWXDIST
      REAL            NEWYDISTFROMCENTER,NEWZDIST
      REAL            SNDYLFT,SNDYRGT,SNDZTOP,SNDZBOT
      REAL            TAGDEPTH,YDIST

      INTEGER         ANUM
                                                                      ! (3*NFISH) = Only a Guess; Impossible to Know How
                                                                      !   Many Fish will be Detected during an HA Survey


      OLDXDIST           = NODES(FKMP,FIMP,1) + OLDFXLOC              ! Previous Distance Downstream from Dam (m) to Location of Fish
      OLDYDISTFROMCENTER = OLDFYLOC - B(FKMP,FIMP)/2                  ! Previous Width Distance (m) from Lake Center
      OLDZDIST           = (NODES(FKMP,FIMP,2)-NODES(KTWBF,FIMP,2))&
                            - OLDFZLOC                                ! Previous Depth (m) of Fish
      NEWXDIST           = NODES(FKMP,FIMP,1) + FXLOC                 ! New Distance Downstream from Dam (m) to Location of Fish
      NEWYDISTFROMCENTER = FYLOC - B(FKMP,FIMP)/2                     ! New Width Distance (m) from Lake Center
      NEWZDIST           = (NODES(FKMP,FIMP,2)-NODES(KTWBF,FIMP,2))&
                            - FZLOC                                   ! New Depth (m) of Fish

      DO 23 ANUM=1,NUMACOUSTICS                                       ! NUMACOUSTICS = Total # of Hydroacoustic Surveys
        IF (HAOPERAT(ANUM,13).EQ.1) THEN                              ! Hydroacoustic Survey is Active
          SNDYLFT = HAOPERAT(ANUM,8)                                  ! Meters from Center of Lake HA Beam will extend to Left Bank
          SNDYRGT = HAOPERAT(ANUM,9)                                  ! Meters from Center of Lake HA Beam will extend to Right Bank
          SNDZTOP = HAOPERAT(ANUM,10)*-1                              ! Depth Below Water Surface (m) to Top of HA Beam: (+) Number
          SNDZBOT = HAOPERAT(ANUM,11)*-1                              ! Depth Below Water Surface (m) to Bottom of HA Beam: (+) Number
          TAGDEPTH = 0
          IF (((OLDXDIST.LE.OLDHALOC(ANUM,1)).AND.&
               (NEWXDIST.GE.HAOPERAT(ANUM,5))).OR.&                    ! Fish May Have Crossed HA Beam Approaching from Upstream
              ((OLDXDIST.GE.OLDHALOC(ANUM,2)).AND.&
               (NEWXDIST.LE.HAOPERAT(ANUM,6)))) THEN                  ! Fish May Have Crossed HA Beam Approaching from Downstream
            TAGDEPTH = (NEWZDIST+OLDZDIST) / 2                        ! Approx. Depth of Fish when it crosses HA Beam
            YDIST = (NEWYDISTFROMCENTER+OLDYDISTFROMCENTER) / 2       ! Approx. Location of Fish in Width Dimension when crossing HA Beam
          END IF
          IF ((TAGDEPTH.NE.0).OR.&
              ((NEWXDIST.GE.HAOPERAT(ANUM,5)).AND.&
               (NEWXDIST.LE.HAOPERAT(ANUM,6)))) THEN
            IF (TAGDEPTH.EQ.0) THEN
              TAGDEPTH = NEWZDIST
              YDIST = NEWYDISTFROMCENTER
            END IF
            IF (((YDIST.GE.SNDYLFT).AND.(YDIST.LE.SNDYRGT)).AND.&
                ((TAGDEPTH.LE.SNDZTOP).AND.(TAGDEPTH.GE.SNDZBOT))) THEN ! Fish is Detected by Hydroacoustic Beam
              SNDCOUNT(ANUM) = SNDCOUNT(ANUM) + 1                       ! Tally # of Fish Detected by Hydroacoustic Beam
              SNDCATCH(ANUM,SNDCOUNT(ANUM),1) = SNDCOUNT(ANUM)          ! Tally # of Fish Detected by Hydroacoustic Beam
              SNDCATCH(ANUM,SNDCOUNT(ANUM),2) = TAGDEPTH                ! Depth (m) at Which Fish is Detected by Hydroacoustic Beam
            END IF
          END IF
        END IF
   23 CONTINUE

      RETURN
      END
!*****************************************
Subroutine fishoutput
!*****************************************
! Numerical Fish Surrogate Output
!*****************************************

Use FISHY; USE GLOBAL
IMPLICIT NONE

Real ::  MAXTEMP, MAXDO, MINTEMP, MINDO, DINTCOUNT, AVEDEPTH,AVEINTDO,AVEINTTEMP,NETDROPDAY
REAL ::  NETPULLDAY, DATEND,DATSTRT,HASPEED,HRSEND,HRSSTRT,NETDROPHR,NETX,NETYLFT,NETYRGT,NETZBOT,NETZTOP
REAL ::  SNDYLFT, SNDYRGT,SNDZTOP,SNDZBOT, NETPULLHR
INTEGER :: CATCH,DINT,DINTLAST,CATCHDEPTH,DOTALLY,TEMPTALLY, FIJ, JF, TAGGED, FVAR2,FVAR3, TAGGEDDEPTH

if(fishon.eq.' ON')then

!Output to Screen After Simulation Finishes                                           !FISH

      WRITE(*,*) '  '
      WRITE(*,*) '  '
      WRITE(*,*) '  '
      WRITE(*,*) 'RESULTS: How Often Factors were Dominant.'
      WRITE(*,*) '  X-direction:'
      WRITE(*,9500) (FXHVCOUNT/TOTFXWGT)*100
      WRITE(*,9501) (FXVVCOUNT/TOTFXWGT)*100
      WRITE(*,9503) (FXTPCOUNT/TOTFXWGT)*100
      WRITE(*,9502) (FXDOCOUNT/TOTFXWGT)*100
      WRITE(*,9504) (FXRDCOUNT/TOTFXWGT)*100
      WRITE(*,9556) (FX00COUNT/TOTFXWGT)*100
      WRITE(*,*) '  '
      WRITE(*,*) '  Z-direction:'
      WRITE(*,9500) (FZHVCOUNT/TOTFZWGT)*100
      WRITE(*,9501) (FZVVCOUNT/TOTFZWGT)*100
      WRITE(*,9503) (FZTPCOUNT/TOTFZWGT)*100
      WRITE(*,9502) (FZDOCOUNT/TOTFZWGT)*100
      WRITE(*,9504) (FZRDCOUNT/TOTFZWGT)*100
      WRITE(*,9556) (FZ00COUNT/TOTFZWGT)*100
      WRITE(*,*) '  '
      WRITE(*,*) '  '
      WRITE(*,*) 'RESULTS: % of the Time the Following Scaled Gradients&
      were Truncated to 1.0:'
      WRITE(*,9510) (XHVTRUNC/COUNTFORTRUC)*100, XHVTRUCSUM/XHVTRUNC
      WRITE(*,9512) (XVVTRUNC/COUNTFORTRUC)*100, XVVTRUCSUM/XVVTRUNC
      WRITE(*,9516) (XTPTRUNC/COUNTFORTRUC)*100, XTPTRUCSUM/XTPTRUNC
      WRITE(*,9514) (XDOTRUNC/COUNTFORTRUC)*100, XDOTRUCSUM/XDOTRUNC
      WRITE(*,9511) (ZHVTRUNC/COUNTFORTRUC)*100, ZHVTRUCSUM/ZHVTRUNC
      WRITE(*,9512) (ZVVTRUNC/COUNTFORTRUC)*100, ZVVTRUCSUM/ZVVTRUNC
      WRITE(*,9516) (ZTPTRUNC/COUNTFORTRUC)*100, ZTPTRUCSUM/ZTPTRUNC
      WRITE(*,9514) (ZDOTRUNC/COUNTFORTRUC)*100, ZDOTRUCSUM/ZDOTRUNC
      WRITE(*,*) '  '
      WRITE(*,*) '  '
      WRITE(*,9340) NZONES
      WRITE(*,*) '  '
      WRITE(*,*) '  '
      WRITE(*,*) '  '

 9500 FORMAT('    Horz Velocity Dominant --> ',F6.1,'%')
 9501 FORMAT('    Vert Velocity Dominant --> ',F6.1,'%')
 9502 FORMAT('    Diss Oxygen   Dominant --> ',F6.1,'%')
 9503 FORMAT('    Temperature   Dominant --> ',F6.1,'%')
 9504 FORMAT('    Randomization Dominant --> ',F6.1,'%')
 9556 FORMAT('    Nothing       Dominant --> ',F6.1,'%')
 9510 FORMAT('   X-direction:',/&
            ,'     Horz Velocity Gradient ----------------->',F6.1,'%',/&
            ,'       ABS(Average Exceedence Value) = ',F6.3)
 9511 FORMAT('   Z-direction:',/&
            ,'     Horz Velocity Gradient ----------------->',F6.1,'%',/&
            ,'       ABS(Average Exceedence Value) = ',F6.3)
 9512 FORMAT('     Vert Velocity Gradient ----------------->',F6.1,'%',/&
            ,'       ABS(Average Exceedence Value) = ',F6.3)
 9514 FORMAT('     Dissolved Oxygen Gradient -------------->',F6.1,'%',/&
            ,'       ABS(Average Exceedence Value) = ',F6.3,/)
 9516 FORMAT('     Temperature Gradient ------------------->',F6.1,'%',/&
            ,'       ABS(Average Exceedence Value) = ',F6.3)
 9340 FORMAT('Enter NZONES =',I6,&
             ' into TecPlot Macro "Animate_Temp_Flood.mcr"')


!Output to File After Simulation Finishes (Exact Same as Screen Output)             !FISH

      WRITE(DIAGFN,*) '  '
      WRITE(DIAGFN,*) '  '
      WRITE(DIAGFN,*) 'RESULTS: How Often Factors were Dominant.'
      WRITE(DIAGFN,*) '  X-direction:'
      WRITE(DIAGFN,9500) (FXHVCOUNT/TOTFXWGT)*100
      WRITE(DIAGFN,9501) (FXVVCOUNT/TOTFXWGT)*100
      WRITE(DIAGFN,9503) (FXTPCOUNT/TOTFXWGT)*100
      WRITE(DIAGFN,9502) (FXDOCOUNT/TOTFXWGT)*100
      WRITE(DIAGFN,9504) (FXRDCOUNT/TOTFXWGT)*100
      WRITE(DIAGFN,9556) (FX00COUNT/TOTFXWGT)*100
      WRITE(DIAGFN,*) '  '
      WRITE(DIAGFN,*) '  Z-direction:'
      WRITE(DIAGFN,9500) (FZHVCOUNT/TOTFZWGT)*100
      WRITE(DIAGFN,9501) (FZVVCOUNT/TOTFZWGT)*100
      WRITE(DIAGFN,9503) (FZTPCOUNT/TOTFZWGT)*100
      WRITE(DIAGFN,9502) (FZDOCOUNT/TOTFZWGT)*100
      WRITE(DIAGFN,9504) (FZRDCOUNT/TOTFZWGT)*100
      WRITE(DIAGFN,9556) (FZ00COUNT/TOTFZWGT)*100
      WRITE(DIAGFN,*) '  '
      WRITE(DIAGFN,*) '  '
      WRITE(DIAGFN,*) 'RESULTS: % of the Time the Following Scaled&
      Gradients were Truncated to 1.0:'
      WRITE(DIAGFN,9510) (XHVTRUNC/COUNTFORTRUC)*100, XHVTRUCSUM/XHVTRUNC
      WRITE(DIAGFN,9512) (XVVTRUNC/COUNTFORTRUC)*100, XVVTRUCSUM/XVVTRUNC
      WRITE(DIAGFN,9516) (XTPTRUNC/COUNTFORTRUC)*100, XTPTRUCSUM/XTPTRUNC
      WRITE(DIAGFN,9514) (XDOTRUNC/COUNTFORTRUC)*100, XDOTRUCSUM/XDOTRUNC
      WRITE(DIAGFN,9511) (ZHVTRUNC/COUNTFORTRUC)*100, ZHVTRUCSUM/ZHVTRUNC
      WRITE(DIAGFN,9512) (ZVVTRUNC/COUNTFORTRUC)*100, ZVVTRUCSUM/ZVVTRUNC
      WRITE(DIAGFN,9516) (ZTPTRUNC/COUNTFORTRUC)*100, ZTPTRUCSUM/ZTPTRUNC
      WRITE(DIAGFN,9514) (ZDOTRUNC/COUNTFORTRUC)*100, ZDOTRUCSUM/ZDOTRUNC
      WRITE(DIAGFN,*) '  '
      WRITE(DIAGFN,*) '  '
      WRITE(DIAGFN,9340) NZONES


!Manipulate Virtual Gillnet Sampling Data for Output                                !FISH

      OPEN (GILFN,FILE='VGILLNETS.OUT',STATUS='UNKNOWN')
      WRITE(GILFN,9563) 
 9563 FORMAT('TITLE = "Gillnet Results"')
      WRITE(GILFN,9508) 
 9508 FORMAT('VARIABLES = "Ave Depth", "Num of Fish", "Ave Temp",&
      "Ave DO", "Max Temp", "Min Temp", "Max DO", "Min DO"')

      IF (GILLNETSAMPLING) THEN                               ! If GILLNETSAMPLING = .TRUE., Virtual Gillnets were Active
      DO 11 NET=1,NUMGILLNETS                                 ! NUMGILLNETS = Total # of Gillnets
        NETDROPDAY = GILLNETS(NET,1)                          ! JDay (INTEGER) to START Gillnet Sampling
        NETPULLDAY = GILLNETS(NET,2)                          ! JDay (INTEGER) to END Gillnet Sampling
        NETDROPHR  = GILLNETS(NET,3)                          ! Hour (Military Time) to START Gillnet Sampling
        NETPULLHR  = GILLNETS(NET,4)                          ! Hour (Military Time) to END Gillnet Sampling
        NETX       = GILLNETS(NET,5)                          ! X-Location (m) Downstream from RBR Dam Where Gillnet is Placed
        NETYLFT    = GILLNETS(NET,6)                          ! Meters from Center of Lake Gillnet will extend to Left Bank: (-) Number
        NETYRGT    = GILLNETS(NET,7)                          ! Meters from Center of Lake Gillnet will extend to Right Bank: (+) Number
        NETZTOP    = GILLNETS(NET,8)*-1                       ! Depth Below Water Surface (m) to Top of Gillnet
        NETZBOT    = GILLNETS(NET,9)*-1                       ! Depth Below Water Surface (m) to Bottom of Gillnet
        IF (NET.NE.1) WRITE(GILFN,*) '  '
        WRITE(GILFN,9505) NET,INT(-NETZBOT/DEPTHINT)+1
 9505   FORMAT('ZONE T="Gillnet #',I4,'", I=',I4,', F=POINT')
        DINTLAST = 0
        DO DINT=-DEPTHINT,(NETZBOT-DEPTHINT),-DEPTHINT
          DINTCOUNT = 0.0                                       ! Used to Tally the Number of Fish in Each Depth Interval
          TEMPTALLY = 0                                       ! Used to Sum the Temperatures where Each Fish is Caught in Gillnet
          DOTALLY   = 0                                       ! Used to Sum the Dissolved Oxygens where Each Fish is Caught in Gillnet
          MAXTEMP   = 0                                       ! Used to Track the Max Temperature in Each Depth Interval
          MINTEMP   = 1E6                                     ! Used to Track the Min Temperature in Each Depth Interval
          MAXDO     = 0                                       ! Used to Track the Max Dissolved Oxygen in Each Depth Interval
          MINDO     = 1E6                                     ! Used to Track the Min Dissolved Oxygen in Each Depth Interval
          
          DO CATCH=1,SNAGCOUNT(NET)                        ! SNAGCOUNT(NET) = The Total # of Fish Caught in Each Gillnet
            CATCHDEPTH = NETCATCH(NET,CATCH,3)                ! Depth Below Water Surface (m) of Fish Caught in Gillnet
            IF ((CATCHDEPTH.LE.DINTLAST).AND.(CATCHDEPTH.GT.DINT)) THEN
              DINTCOUNT = DINTCOUNT + 1.0
              TEMPTALLY = TEMPTALLY + NETCATCH(NET,CATCH,4)
              DOTALLY   = DOTALLY   + NETCATCH(NET,CATCH,5)
              IF (NETCATCH(NET,CATCH,4).GT.MAXTEMP) &
                  MAXTEMP = NETCATCH(NET,CATCH,4)
              IF (NETCATCH(NET,CATCH,4).LT.MINTEMP) &
                  MINTEMP = NETCATCH(NET,CATCH,4)
              IF (NETCATCH(NET,CATCH,5).GT.MAXDO) &
                  MAXDO   = NETCATCH(NET,CATCH,5)
              IF (NETCATCH(NET,CATCH,5).LT.MINDO) &
                  MINDO   = NETCATCH(NET,CATCH,5)
            END IF
          ENDDO
          
          AVEDEPTH   = DINTLAST + (DINT-DINTLAST)/2
          IF (DINTCOUNT.GT.0) THEN
            AVEINTTEMP = TEMPTALLY / DINTCOUNT
            AVEINTDO   = DOTALLY   / DINTCOUNT
          ELSE
            AVEINTTEMP = -99
            AVEINTDO   = -99
            MAXTEMP    = -99
            MINTEMP    = -99
            MAXDO      = -99
            MINDO      = -99
          END IF
          WRITE(GILFN,9506) AVEDEPTH,INT(DINTCOUNT),AVEINTTEMP,AVEINTDO,&
                           MAXTEMP,MINTEMP,MAXDO,MINDO
 9506     FORMAT(F6.2,4X,I7,6(4X,F6.2))
          DINTLAST = DINT
       ENDDO
        WRITE(GILFN,9564) NET,NET
 9564   FORMAT('TEXT X=49, Y=97, F=HELV-BOLD, CS=FRAME, HU=FRAME,&
      AN=MIDCENTER, C=CUSTOM7, H=3, T="Gillnet #',I4,'", ZN=',I4)
        WRITE(GILFN,9574) NETX,NET
 9574   FORMAT('TEXT X=49, Y=94.5, F=HELV-BOLD, CS=FRAME, HU=FRAME,&
      AN=MIDCENTER, C=BLACK, H=2, T="Downstream Location =',F8.1,&
      ' (m)", ZN=',I4)
        WRITE(GILFN,9507) INT(NETDROPDAY),INT(NETPULLDAY),&
                         NETDROPHR,NETPULLHR,NET
 9507   FORMAT('TEXT X=49, Y=92, F=HELV-BOLD, CS=FRAME, HU=FRAME,&
      AN=MIDCENTER, C=BLACK, H=2, T="JDAY:',I4,' -',I4,'      TIME:',&
      F5.1,' hrs -',F5.1,' hrs", ZN=',I4)
        WRITE(GILFN,9572) SNAGCOUNT(NET),&
                         (REAL(SNAGCOUNT(NET))/REAL(NFISH))*100,NET
 9572   FORMAT('TEXT X=49, Y=89.5, F=HELV-BOLD, CS=FRAME, HU=FRAME,&
      AN=MIDCENTER, C=BLACK, H=2, T="# Fish Caught =',I7,&
      '      % of Total Population =',F7.2,'", ZN=',I4)
        WRITE(GILFN,9575) NETZTOP,NETZBOT,NET
 9575   FORMAT('TEXT X=49, Y=87, F=HELV-BOLD, CS=FRAME, HU=FRAME,&
      AN=MIDCENTER, C=BLACK, H=2, T="Top of Gillnet = ',F6.2,&
      ' (m)      Bottom of Gillnet = ',F6.2,' (m)", ZN=',I4)
        WRITE(GILFN,9576) (NETYRGT-NETYLFT),NET
 9576   FORMAT('TEXT X=49, Y=84.5, F=HELV-BOLD, CS=FRAME, HU=FRAME,&
      AN=MIDCENTER, C=BLACK, H=2, T="Width of Gillnet =',F8.2,&
      ' (m)", ZN=',I4)
   11 CONTINUE
      END IF


!Manipulate Virtual Hydroacoustics Sampling Results for Output                      !FISH

      OPEN (HYACFN,FILE='VACOUSTICS.OUT',STATUS='UNKNOWN')
      WRITE(HYACFN,9565) 
 9565 FORMAT('TITLE = "Hydroacoustic Survey Results"')
      WRITE(HYACFN,9558) 
 9558 FORMAT('VARIABLES = "Ave Depth", "Num of Fish"')

      IF (ACOUSTICSAMPLING) THEN
      DO 18 SND=1,NUMACOUSTICS                              ! NUMACOUSTICS = Total # of Hydroacoustic Surveys
        DATSTRT = HAOPERAT(SND,1)                           ! JDay (INTEGER) to START Hydroacoustic Survey
        DATEND  = HAOPERAT(SND,2)                           ! JDay (INTEGER) to END Hydroacoustic Survey
        HRSSTRT = HAOPERAT(SND,3)                           ! Hour (Military Time) to START Hydroacoustic Survey
        HRSEND  = HAOPERAT(SND,4)                           ! Hour (Military Time) to END Hydroacoustic Survey
        HASPEED = HAOPERAT(SND,7)                           ! Speed (m/s) of Hydroacoustic Sampling 'Box'
        SNDYLFT = HAOPERAT(SND,8)                           ! Meters from Center of Lake Hydroacoustic Sampling will extend to Left Bank
        SNDYRGT = HAOPERAT(SND,9)                           ! Meters from Center of Lake Hydroacoustic Sampling will extend to Right Bank
        SNDZTOP = HAOPERAT(SND,10)*-1                       ! Depth Below Water Surface (m) to Top of Hydroacoustic Sampling Zone
        SNDZBOT = HAOPERAT(SND,11)*-1                       ! Depth Below Water Surface (m) to Bottom of Hydroacoustic Sampling Zone
        IF (SND.NE.1) WRITE(HYACFN,*) '  '
        WRITE(HYACFN,9559) SND,INT(-SNDZBOT/HADEPTHINT)+1
 9559   FORMAT('ZONE T="HA Survey #',I4,'", I=',I4,', F=POINT')
        DINTLAST = 0
        DO 19 DINT=-HADEPTHINT,(SNDZBOT-HADEPTHINT),-HADEPTHINT
          DINTCOUNT = 0.0                                   ! Used to Tally the Number of Fish in Each Depth Interval
          DO 21 TAGGED=1,SNDCOUNT(SND)                      ! SNDCOUNT(SND) = Total # of Fish Detected in Each HA Survey
            TAGGEDDEPTH = SNDCATCH(SND,TAGGED,2)
            IF ((TAGGEDDEPTH.LE.DINTLAST).AND.(TAGGEDDEPTH.GT.DINT)) THEN
              DINTCOUNT = DINTCOUNT + 1.0                   ! Tally the # of Fish in Each Depth Interval
            END IF
   21     CONTINUE
          AVEDEPTH = DINTLAST + (DINT-DINTLAST)/2
          WRITE(HYACFN,9560) AVEDEPTH,INT(DINTCOUNT)
 9560     FORMAT(F6.2,4X,I8)
          DINTLAST = DINT
   19   CONTINUE
        WRITE(HYACFN,9561) SND,SND
 9561   FORMAT('TEXT X=49, Y=97, F=HELV-BOLD, CS=FRAME, HU=FRAME,&
      AN=MIDCENTER, C=CUSTOM7, H=3, T="Hydroacoustic Survey #',I4,&
      '", ZN=',I4)
        WRITE(HYACFN,9568) (HAOPERAT(SND,5)+HAOPERAT(SND,6)) / 2,SND
 9568   FORMAT('TEXT X=49, Y=94, F=HELV-BOLD, CS=FRAME, HU=FRAME,&
      AN=MIDCENTER, C=BLACK, H=2, T="HA Survey Box: Final Location =',&
      F8.1,' (m)", ZN=',I4)
        WRITE(HYACFN,9569) (HAOPERAT(SND,6)-HAOPERAT(SND,5)),SND
 9569   FORMAT('TEXT X=53.87, Y=91.5, F=HELV-BOLD, CS=FRAME, HU=FRAME,&
      AN=MIDCENTER, C=BLACK, H=2, T="Length =',F8.1,' (m)", ZN=',I4)
        WRITE(HYACFN,9570) HAOPERAT(SND,7),SND
 9570   FORMAT('TEXT X=49, Y=89, F=HELV-BOLD, CS=FRAME, HU=FRAME,&
      AN=MIDCENTER, C=BLACK, H=2, T="Speed of HA Survey Box = ',&
      F6.2,' (m/s)", ZN=',I4)
        WRITE(HYACFN,9562) INT(DATSTRT),INT(DATEND),HRSSTRT,HRSEND,SND
 9562   FORMAT('TEXT X=49, Y=86.5, F=HELV-BOLD, CS=FRAME, HU=FRAME,&
      AN=MIDCENTER, C=BLACK, H=2, T="JDAY:',I4,' -',I4,&
      '      TIME:',F5.1,' hrs -',F5.1,' hrs", ZN=',I4)
        WRITE(HYACFN,9573) SNDCOUNT(SND),&
                         (REAL(SNDCOUNT(SND))/REAL(NFISH))*100,SND
 9573   FORMAT('TEXT X=49, Y=84, F=HELV-BOLD, CS=FRAME, HU=FRAME,&
      AN=MIDCENTER, C=BLACK, H=2, T="# Fish Detected =',I7,&
      '      % of Total Population =',F7.2,'", ZN=',I4)
   18 CONTINUE
      END IF


!Close All the Output Files Used for the Numerical Fish Surrogate

      CLOSE(DIAGFN)                                                                   !FISH
      CLOSE(DATADEBUGFN)                                                                   !FISH
      CLOSE(BARCHRTXFN)                                                                   !FISH
      CLOSE(BARCHRTZFN)                                                                   !FISH
      DO FIJ=1,NBR                                                                  !FISH
        FVAR2 = 20000+FIJ
        FVAR3 = 30000+FIJ
        CLOSE(FVAR2)
        CLOSE(FVAR3)
      END DO
      CLOSE(GILFN)                                                                   !FISH
      CLOSE(HYACFN) 

! final fish output
    do jf=1,3   !nfish
    write(FINALFN,'(i7,1x,<fpara>(f10.3,1x))')jf,(fishes(jf,i),i=1,fpara)
    end do
    close(FINALFN)
! End fish output section

end if                                                                  !FISH

      Return                                                               !FISH
End Subroutine FishOutput

! ************************ READ FISH INPUT DATA FILES

Subroutine Read_Fish_Data

Use Fishy
IMPLICIT NONE

character*3 Char(10)
INTEGER :: NG,I,J,NA

    open(DIAGFN,file='fish.npt',status='old')

    READ(DIAGFN,'(/)')

    READ(DIAGFN,'(/10x,7x,a3,i10,i10,f10.0,2(7x,a3),3f10.0)')FISHON,NFISHPCEL,NFISHSEG,DELAYDATE,char(1),char(2),SEDVEL,ALPHAX,ALPHAZ
    If(char(1).eq.' ON')then
    PARTICLE=.true.
    else
    PARTICLE=.false.
    end if
    If(char(2).eq.' ON')then
    LINE=.true.
    else
    LINE=.false.
    end if

! LINE: If LINE, distribute NFISHCEL not at a point but linearly along the line defined by IFISH and IFISHT and IFISHB; IF LINE, FZLOC=0
! ALPHAX/ALPHAZ: Turbulent Disp Coeff multipliers for random water movement for x and z directions, DISPX=ALPHAX*DX, DISPZ=ALPHAZ*DZ
! PARTICLE = ON only "dumb" particle transport with random fluid motion + Sedimentation Vel (SEDVEL)
! NFISHPCEL: # of fish added per cel
! NFISHSEG: # of segements to add fish to
! DELAYDATE:  Delay Fish Behavioral Movement Until After this Date.

  if(fishon.ne.'OFF')then

    if(nfishseg.eq.0)nfishseg=1
    Allocate(ifish(nfishseg),ifisht(nfishseg),ifishb(nfishseg))

    READ(DIAGFN,'(//(10X,9I10))')(IFISH(I),I=1,NFISHSEG)
    READ (DIAGFN,'(//(10X,9I10))')(IFISHT(I),I=1,NFISHSEG)
    READ (DIAGFN,'(//(10X,9I10))')(IFISHB(I),I=1,NFISHSEG)

! IFISH: SEG#s where fish are added
! IFISHT:Uppermost top Layer where Fish are added in Segment IFISH
! IFISHB: Bottommost Layer where Fish are added in Segment IFISH

    NFISH=0
    do i=1,nfishseg
    NFISH=NFISH+(ifishb(i)-ifisht(i)+1)*nfishpcel
    end do

    IF(PARTICLE)THEN
        READ(DIAGFN,1001)FXLOC,FZLOC,OUTFREQP
    ELSE
        
    READ(DIAGFN,1001)FXLOC,FZLOC,FSIZE,FAGE,UFISH,VFISH,WFISH
    
    ENDIF
    

!          ! FXLOC   = Location of fish within segment IMP from upstream side
!          ! FZLOC   = Location of fish within layer KMP from top side
!          ! FNBP    = Branch where fish is released
!          ! FSIZE   = Size (i.e., length) of fish in meters (1 inch = 0.0254 meters)
!          ! FAGE    = Age of the fish in years internally changed to days
           fage=fage*365.
!          ! UFISH   = Initial longitudinal velocity of the fish relative to water
!          ! VFISH   = Initial lateral velocity of the fish relative to water
!          ! WFISH   = Initial vertical velocity of the fish relative to water

    READ(DIAGFN,1000)NDT    ! SW 
    READ(DIAGFN,1002)Char(1)
    !If(char(1).eq.' ON')then
    !MULTIRESPONSE=.true.
    !else
    !MULTIRESPONSE=.false.
    !end if
    READ(DIAGFN,1001)HVXWEIGT,VVXWEIGT
    READ(DIAGFN,1001)HVZWEIGT,VVZWEIGT
    READ(DIAGFN, 1001)TPXWEIGT,TPZWEIGT,TEMPTHRES,TSTEP1,TSTEP1MULT,TSTEP2,TSETP2MULT 
    READ(DIAGFN,1001)  DOXWEIGT,  DOZWEIGT,  DOTHRES2, DOSTEP1MULT  
    READ(DIAGFN,1001) RDXWEIGT,  RDYWEIGT,  RDZWEIGT, EPSILONRD  
!      NFSFREQ  = 0.006250            ! The number (or fraction thereof) of JDAYs between successive runs of the NFS module
! REPLACED WITH NDT
!______________________________________________________
!  SENSORY WEIGHT and SCALING Parameters
!      MULTIPLERESPONSE = .FALSE.     ! If .TRUE., URGENCY is Calculated as a Weighted Average of the
                                     !            Selected Scaled Horz Vel, Vert Vel, DO, and/or Temp Gradients
                                     ! If .FALSE., URGENCY is Calculated as the Maximum of the following:
                                     !             Horz Vel, Vert Vel, DO, and/or Temp Gradients
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Horizontal Movement Parameters
!      HVXWEIGT = 1.0                 ! HVXWEIGT = Horiz Vel Weight for X-direction         (Range: 0 --> 1)
!      VVXWEIGT = 0.3                 ! VVXWEIGT = Vert Vel Weight for X-direction          (Range: 0 --> 1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Vertical Movement Parameters
!      HVZWEIGT = 0.2                 ! HVZWEIGT = Horiz Vel Weight for Z-direction         (Range: 0 --> 1)
!      VVZWEIGT = 0.1                 ! VVZWEIGT = Vert Vel Weight for Z-direction          (Range: 0 --> 1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Temperature Parameters
!      TPXWEIGT = 0.90                ! TPXWEIGT = Temperature Weight for X-direction       (Range: 0 --> 1)
!      TPZWEIGT = 0.75                ! TPZWEIGT = Temperature Weight for Z-direction       (Range: 0 --> 1)
!      TEMPTHRES     = 26.0           ! Temperature Threshold, Above Which Fish Swims with Maximum Urgency Towards Cooler Water
!      TSTEP1        = 1.5            ! ABS(Temperature Difference) from Optimum Temperature where Urgency is Multiplied by TSTEP1MULT
!      TSTEP1MULT    = 0.0            ! Urgency is Multiplied by TSTEP1MULT if the Fish is Within TSTEP1 of the Optimum Temperature
!      TSTEP2        = 3.5            ! ABS(Temperature Difference) from Optimum Temperature where Urgency is Multiplied by TSTEP2MULT
!      TSTEP2MULT    = 0.5            ! Urgency is Multiplied by TSTEP2MULT if the Fish is Within TSTEP2 of the Optimum Temperature
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Dissolved Oxygen Parameters
!      DOXWEIGT = 0.76                ! DOXWEIGT = Dissolved Oxygen Weight for X-direction  (Range: 0 --> 1)
!      DOZWEIGT = 0.12                ! DOZWEIGT = Dissolved Oxygen Weight for Z-direction  (Range: 0 --> 1)
!      DOTHRES2      = 6.5            ! Dissolved Oxygen Threshold, Above Which Fish Swims with a Decremented Urgency
!      DOSTEP1MULT   = 0.1            ! Urgency is Multiplied by DOSTEP1MULT if DO at Fish Location is Greater Than DOTHRES2
!- - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Random Displacement Parameters
!      RDXWEIGT  = 0.999              ! Random Movement Weight for X-direction              (Range: 0 --> 1)
!      RDYWEIGT  = 0.5                ! Random Movement Weight for Y-direction              (Range: 0 --> 1)
!      RDZWEIGT  = 0.47               ! Random Movement Weight for Z-direction              (Range: 0 --> 1)
!      EPSILONRD = 1E-8               ! Random Parameter Epsilon - Used only to keep RDX, RDY, and/or RDZ from getting stuck on 0
                                     !   when UFISH, VFISH, or WFISH = 0, respectively.

    READ(DIAGFN,1001) MAXREACTXHV, MINREACTXHV, DISTREACTXHV, MAXREACTZHV, MINREACTZHV, DISTREACTZHV
    READ(DIAGFN,1001) MAXREACTXVV, MINREACTXVV, DISTREACTXVV, MAXREACTZVV, MINREACTZVV, DISTREACTZVV
    READ(DIAGFN,1001) MAXREACTXTP, MINREACTXTP, DISTREACTXTP, MAXREACTZTP, MINREACTZTP, DISTREACTZTP
    READ(DIAGFN,1001) MAXREACTXDO, MINREACTXDO, DISTREACTXDO, MAXREACTZDO, MINREACTZDO, DISTREACTZDO

!  The following values are used to scale the gradients each fish detects within its sensory sphere.  The
!   parameters immediately below should be valued according to the gradient likely to induce the maximum fish
!   response, with higher magnitude gradients not inducing any more of a response from the species of fish
!   under study.  The following parameters define a representative example of the practical maximum gradient
!   likely to induce the greatest fish response/reaction/movement.  ALL NUMBERS MUST BE POSITIVE !!!
!   Horizontal Velocity
!!    X-direction
!      MAXREACTXHV   = 1E-2           ! MAXREACTXHV  = For example gradient, the greatest Horiz Velocity value (m/s) in X-direction
!      MINREACTXHV   = 1E-3           ! MINREACTXHV  = For example gradient, the lowest Horiz Velocity value (m/s) in X-direction
!      DISTREACTXHV  = 250            ! DISTREACTXHV = For example gradient, the distance between MAXREACTXHV and MINREACTXHV (m)
!    X-direction
!      MAXREACTXHV   = 0.035          ! MAXREACTXHV  = For example gradient, the greatest Horiz Velocity value (m/s) in X-direction
!      MINREACTXHV   = 0              ! MINREACTXHV  = For example gradient, the lowest Horiz Velocity value (m/s) in X-direction
!      DISTREACTXHV  = 1              ! DISTREACTXHV = For example gradient, the distance between MAXREACTXHV and MINREACTXHV (m)
!    Z-direction
!      MAXREACTZHV   = 1E-2           ! MAXREACTZHV  = For example gradient, the greatest Horiz Velocity value (m/s) in Z-direction
!      MINREACTZHV   = 1E-3           ! MINREACTZHV  = For example gradient, the lowest Horiz Velocity value (m/s) in Z-direction
!      DISTREACTZHV  = 0.135          ! DISTREACTZHV = For example gradient, the distance between MAXREACTZHV and MINREACTZHV (m)
!  -   -   -   -   -   -   -   -   -   -   -   -   -
!   Vertical Velocity
!    X-direction
!      MAXREACTXVV   = 1E-2           ! MAXREACTXVV  = For example gradient, the greatest Vert Velocity value (m/s) in X-direction
!      MINREACTXVV   = 1E-3           ! MINREACTXVV  = For example gradient, the lowest Vert Velocity value (m/s) in X-direction
!      DISTREACTXVV  = 3000           ! DISTREACTXVV = For example gradient, the distance between MAXREACTXVV and MINREACTXVV (m)
!    Z-direction
!      MAXREACTZVV   = 1E-2           ! MAXREACTZVV  = For example gradient, the greatest Vert Velocity value (m/s) in Z-direction
!      MINREACTZVV   = 1E-3           ! MINREACTZVV  = For example gradient, the lowest Vert Velocity value (m/s) in Z-direction
!      DISTREACTZVV  = 110            ! DISTREACTZVV = For example gradient, the distance between MAXREACTZVV and MINREACTZVV (m)
!  -   -   -   -   -   -   -   -   -   -   -   -   -
!   Temperature
!    X-direction
!      MAXREACTXTP   = 30             ! MAXREACTXTP  = For example gradient, the greatest Temp value (deg C) in X-direction
!      MINREACTXTP   = 15             ! MINREACTXTP  = For example gradient, the lowest Temp value (deg C) in X-direction
!      DISTREACTXTP  = 1700           ! DISTREACTXTP = For example gradient, the distance between MAXREACTXTP and MINREACTXTP (m)
!    Z-direction
!      MAXREACTZTP   = 30             ! MAXREACTZTP  = For example gradient, the greatest Temp value (deg C) in Z-direction
!      MINREACTZTP   = 15             ! MINREACTZTP  = For example gradient, the lowest Temp value (deg C) in Z-direction
!      DISTREACTZTP  = 3.80           ! DISTREACTZTP = For example gradient, the distance between MAXREACTZTP and MINREACTZTP (m)
!  -   -   -   -   -   -   -   -   -   -   -   -   -
!   Dissolved Oxygen
!    X-direction
!      MAXREACTXDO   = 8              ! MAXREACTXDO  = For example gradient, the greatest Diss Oxyg value (mg/L) in X-direction
!      MINREACTXDO   = 4              ! MINREACTXDO  = For example gradient, the lowest Diss Oxyg value (mg/L) in X-direction
!      DISTREACTXDO  = 7000           ! DISTREACTXDO = For example gradient, the distance between MAXREACTXDO and MINREACTXDO (m)
!    Z-direction
!      MAXREACTZDO   = 8              ! MAXREACTZDO  = For example gradient, the greatest Diss Oxyg value (mg/L) in Z-direction
!      MINREACTZDO   = 4              ! MINREACTZDO  = For example gradient, the lowest Diss Oxyg value (mg/L) in Z-direction
!      DISTREACTZDO  = 3.80           ! DISTREACTZDO = For example gradient, the distance between MAXREACTZDO and MINREACTZDO (m)
!

    READ(DIAGFN,1001) XINZONE, ZINZONE, XMDZONE, ZMDZONE, XOTZONE, ZOTZONE, XDISPCOEFF, ZDISPCOEFF
    READ(DIAGFN,1001) FBDYSEARCH, BBDYSEARCH, UBDYSEARCH, DBDYSEARCH

!  SCHOOLING AND DISPERSION Parameters
!      XINZONE    = 0.25              ! Number of body lengths (FSIZE) to edge of Inner Zone of Repulsion (X-direction)
!      ZINZONE    = 0.25              ! Number of body lengths (FSIZE) to edge of Inner Zone of Repulsion (Z-direction)
!      XMDZONE    = 2.0               ! Number of body lengths (FSIZE) to edge of Middle Zone of Orientation (X-direction)
!      ZMDZONE    = 2.0               ! Number of body lengths (FSIZE) to edge of Middle Zone of Orientation (Z-direction)
!      XOTZONE    = 100.0             ! Number of body lengths (FSIZE) to edge of Outerzone of Attraction (X-direction)
!      ZOTZONE    = 20.0              ! Number of body lengths (FSIZE) to edge of Outerzone of Attraction (Z-direction)
!      XDISPCOEFF = 0.001             ! Scaling Constant used in the X-directional Dispersion (Threshold-Tendency Model) Formulation
!      ZDISPCOEFF = 0.001             ! Scaling Constant used in the Z-directional Dispersion (Threshold-Tendency Model) Formulation
!______________________________________________________
!  SENSORY SPHERE Parameters
!      FBDYSEARCH = 7500              ! The # of body lengths (FSIZE) in front of the fish the fish will likely search
!                                     !   during each timestep (approx 7 min) to determine conditions/gradients
!      BBDYSEARCH = 7500              ! The # of body lengths (FSIZE) behind the fish the fish will likely search
!                                     !   during each timestep (approx 7 min) to determine conditions/gradients
!      UBDYSEARCH = 100               ! The # of body lengths (FSIZE) above the fish the fish will likely search
!                                     !   during each timestep (approx 7 min) to determine conditions/gradients
!      DBDYSEARCH = 100               ! The # of body lengths (FSIZE) below the fish the fish will likely search
!                                     !   during each timestep (approx 7 min) to determine conditions/gradients
!___________
    READ(DIAGFN,1001) TEMPOPTD, TEMPOPTN, MXXSPDL,  MXZSPDL, DOTHRES, XREFL, YREFL, ZBOTREFL, ZSURREFL

!  MISC Parameters
!     TEMPOPTD  = 17.0               ! Optimum Temperature for Fish Species during Day
!     TEMPOPTN  = 22.5               ! Optimum Temperature for Fish Species at Night
!     MXXSPDL  = 12.0                ! # of Fish Lengths covered in X-direction per second at Max Fish Urgency (perfect conditions)
!                                    !   Note:  RDXWEIGT * MXXSPDL = Approximate Fish Cruising Speed (Fish Lengths/sec)
!     MXZSPDL  = 0.13                ! # of Fish Lengths covered in Z-direction per second at Max Fish Urgency (perfect conditions)
!     DOTHRES  = 6.0                 ! DO Threshold below which Maximum Fish Cruising Speed Starts to Suffer
!     XREFL    = 0.1                 ! Reflect this percentage of segment length
!                                    !   when fish encounters horizontal obstacle
!     YREFL    = 0.1                 ! Reflect this % of cell width when encounter stream bank
!     ZBOTREFL = 2.0                 ! Reflect this percentage of bottom layer height
!                                    !   when fish encounters bottom
!     ZSURREFL = 0.5                 ! Reflect this percentage of surface layer height
!                                     !   when fish encounters water surface

    READ(DIAGFN,'(//i10,7x,a3,f10.0)')NUMGILLNETS,char(1), DEPTHINT
    If(char(1).eq.' ON')then
    GILLNETSAMPLING = .TRUE.
    else
    GILLNETSAMPLING = .FALSE.
    end if

    Allocate(gillnets(NUMGILLNETS,12),SNAGCOUNT(NUMGILLNETS))
    gillnets=0.0
    snagcount=0.0

    READ(DIAGFN,'(/)')
    do ng=1,numgillnets
    READ(DIAGFN,'(14f10.0)')(gillnets(ng,j),j=1,9),gillnets(ng,12)
    gillnets(ng,10)=real(ng)
    end do

!        GILLNETSAMPLING = .TRUE.            ! If GILLNETSAMPLING = .TRUE., Virtual Gillnet Sampling is Allowed
!        DEPTHINT        = 1.0               ! Depth Interval (m) Used to Display Gillnet Results; (+) Number
!Gillnet #6
!        GILLNETS(6,1)  = 235              ! JDay (INTEGER) to START Gillnet Sampling
!        GILLNETS(6,2)  = 236              ! JDay (INTEGER) to END Gillnet Sampling
!        GILLNETS(6,3)  = 15.0             ! Hour (Military Time) to START Gillnet Sampling
!        GILLNETS(6,4)  = 9.0              ! Hour (Military Time) to END Gillnet Sampling
!        GILLNETS(6,5)  = 30380.0          ! X-Location (m) Downstream from RBR Dam Where Gillnet is Placed
!        GILLNETS(6,6)  = -3000.0          ! Width (meters from Center of Lake) Gillnet will extend towards Left Bank: (-) Number
!        GILLNETS(6,7)  = 3000.0           ! Width (meters from Center of Lake) Gillnet will extend towards Right Bank: (+) Number
!        GILLNETS(6,8)  = 0.0              ! Depth Below Water Surface (m) to Top of Gillnet: (+) Number
!        GILLNETS(6,9)  = 73.0             ! Depth Below Water Surface (m) to Bottom of Gillnet: (+) Number
!        GILLNETS(6,10) = 6                ! Gillnet #
!        GILLNETS(NET,11) = 0              ! Gillnet Active? [0=No , 1=Yes]
!        GILLNETS(NET,12) = 6              ! Branch # Where Gillnet is Placed
!
    READ(DIAGFN,'(//i10,7x,a3,f10.0)')NUMACOUSTICS,char(1), HADEPTHINT
    if(char(1).eq.' ON')then
    ACOUSTICSAMPLING = .TRUE. 
    else
    ACOUSTICSAMPLING = .FALSE. 
    end if

    Allocate(HAOPERAT(NUMACOUSTICS,14),OLDHALOC(NUMACOUSTICS,2),SNDCOUNT(NUMACOUSTICS))
    haoperat=0.0
    oldhaloc=0.0
    sndcount=0.0

    READ(DIAGFN,'(/)')
    do na=1,NUMACOUSTICS
    READ(DIAGFN,'(14f10.0)')(haoperat(na,j),j=1,11),haoperat(na,14)
    haoperat(na,12)=real(na)
    end do


!        ACOUSTICSAMPLING = .TRUE.           ! If ACOUSTICSAMPLING = .TRUE., Virtual Hydroacoustic Sampling is Allowed
!        HADEPTHINT       = 1.0              ! Depth Interval (m) Used to Display Hydroacoustic Sampling Results; (+) Number
!Hydroacoustic Sample #11
!        HAOPERAT(11,1)  = 227               ! JDay (INTEGER) to START Hydroacoustic Sampling Operation
!        HAOPERAT(11,2)  = 227               ! JDay (INTEGER) to END Hydroacoustic Sampling Operation
!        HAOPERAT(11,3)  = 21.5              ! Hour (Military Time) to START Hydroacoustic Sampling Operation
!        HAOPERAT(11,4)  = 22.12             ! Hour (Military Time) to END Hydroacoustic Sampling Operation
!        HAOPERAT(11,5)  = 59998.7           ! X-Location (m) of Upstream End of Hydroacoustic Sampling 'Box' at Start of Sampling
!        HAOPERAT(11,6)  = 60000.0           ! X-Location (m) of Downstream End of Hydroacoustic Sampling 'Box' at Start of Sampling
!        HAOPERAT(11,7)  = -2.235            ! Speed (m/s) of Hydroacoustic  Sampling 'Box' [ (-)=Moves Upstream , (+)=Moves Downstream ]
!        HAOPERAT(11,8)  = -2000.0           ! Width (meters from Center of Lake) that HA Sampling extends to Left Bank: (-) Number
!        HAOPERAT(11,9)  = 2000.0            ! Width (meters from Center of Lake) that HA Sampling extends to Right Bank: (+) Number
!        HAOPERAT(11,10) = 2.0               ! Depth Below Water Surface (m) to Top of Hydroacoustic Sampling Zone: (+) Number
!        HAOPERAT(11,11) = 60.0              ! Depth Below Water Surface (m) to Bottom of Hydroacoustic Sampling Zone; (+) Number
!                                            !   Overestimate Depth in order to Cover the Entire Water Column
!        HAOPERAT(11,12) = 11                ! Hydroacoustic Sampling Operation #
!        HAOPERAT(SND,13) = 0              ! Hydroacoustic Sampling Active? [0=No , 1=Yes]
!        HAOPERAT(SND,14) = 6              ! Branch # Where Hydroacoustic Sampling is Done

    READ(DIAGFN,1001)OUTFREQJAN,OUTFREQFEB,OUTFREQMAR,OUTFREQAPR,OUTFREQMAY,OUTFREQJUN,OUTFREQJUL,OUTFREQAUG,OUTFREQSEP,OUTFREQOCT,OUTFREQNOV,OUTFREQDEC

! Output as close to this fraction or every number of JDAYs as possible
!   OR enter 99 to skip output for a particular month altogether


    READ(DIAGFN,'(//7x,a3,7x,a3,2i10)')char(1),char(2),VARYTEMP, VARYDO
    if(char(1).eq.' ON')then
    NULLFIELDWQ = .TRUE. 
    else
    NULLFIELDWQ = .FALSE. 
    end if
    if(char(2).eq.' ON')then
    WQNULLLINEAR = .TRUE. 
    else
    WQNULLLINEAR = .FALSE. 
    end if

    READ(DIAGFN,'(//2f10.0,i10,f10.0,i10)') TOPTEMP, MIDKTEMP,  MIDKT,  BOTTEMP,   BOTK
    READ(DIAGFN,'(//2f10.0,i10,f10.0)') LFTTEMP, MIDITEMP,  MIDIT,  RGTTEMP 
    READ(DIAGFN,'(//2f10.0,i10,f10.0,i10)')TOPDO, MIDKDO,MIDKD, BOTDO   
    READ(DIAGFN,'(//2f10.0,i10,f10.0)') LFTDO,    MIDIDO,      MIDID,     RGTDO 

    READ(DIAGFN,'(//f10.0,7x,a3,7x,a3,2i10)')VVELCAP,char(1),char(2),VARYHVEL,VARYVVEL

    if(char(1).eq.' ON')then
    NULLFIELDFF = .TRUE. 
    else
    NULLFIELDFF = .FALSE. 
    end if
    if(char(2).eq.' ON')then
    FFNULLLINEAR = .TRUE. 
    else
    FFNULLLINEAR = .FALSE. 
    end if

    READ(DIAGFN,'(//2f10.0,i10,f10.0,i10)')TOPHVEL,  MIDKHVEL, MIDKH, BOTHVEL, BOTKK  
    READ(DIAGFN,'(//2f10.0,i10,f10.0)')  LFTHVEL,  MIDIHVEL,     MIDIH,   RGTHVEL  
    READ(DIAGFN,'(//2f10.0,i10,f10.0,i10)')TOPVVEL,  MIDKVVEL, MIDKV, BOTVVVEL
    READ(DIAGFN,'(//2f10.0,i10,f10.0)')  LFTVVEL,  MIDIVVEL,     MIDIV,   RGTVVEL  


    READ(DIAGFN,1002)(char(i),i=1,7)
    if(char(1).eq.' ON')then
    STIMULIRULES = .TRUE. 
    else
    STIMULIRULES = .FALSE. 
    end if
    if(char(2).eq.' ON')then
    VELOCITYRULES = .TRUE. 
    else
    VELOCITYRULES = .FALSE. 
    end if
    if(char(3).eq.' ON')then
    TEMPRULES = .TRUE. 
    else   
    TEMPRULES = .FALSE. 
    end if
    if(char(4).eq.' ON')then
    DORULES = .TRUE. 
    else
    DORULES = .FALSE. 
    end if
    if(char(5).eq.' ON')then
    RANDOMIZATION = .TRUE. 
    else
    RANDOMIZATION = .FALSE. 
    end if
    if(char(6).eq.' ON')then
    PASSIVETRANSPORT = .TRUE. 
    else
    PASSIVETRANSPORT = .FALSE. 
    end if
    if(char(7).eq.' ON')then
    SCHOOLING = .TRUE. 
    else
    SCHOOLING = .FALSE. 
    end if
    
    READ(DIAGFN,'(//2(7x,a3),i10,2(7x,a3))') char(1),char(2),WBRUN,char(3),char(4)
    if(char(1).eq.' ON')then
    DEBUG = .TRUE. 
    else
    DEBUG  = .FALSE. 
    end if
    if(char(2).eq.' ON')then
    WBSKIP= .TRUE. 
    else
    WBSKIP = .FALSE. 
    end if
    if(char(3).eq.' ON')then
    LINEAR= .TRUE. 
    else
    LINEAR = .FALSE. 
    end if
    if(char(4).eq.' ON')then
    PREVENTBRCHSWITCH = .TRUE. 
    else
    PREVENTBRCHSWITCH = .FALSE. 
    end if

    READ(DIAGFN,'(//f10.0,7x,a3,4(f10.0))')ASPRATIO,char(1), SKYNIGHT, SKYDAWN, SKYDAY, SKYDUSK  
    if(char(1).eq.' ON')then
    SHOWSKY = .TRUE. 
    else
    SHOWSKY = .FALSE. 
    end if

    READ(DIAGFN,'(//2i10,7x,a3,i10)')UNBP, DNBP, collector, ncollector

    if(ncollector.eq.0)ncollector=1
    Allocate(icoll(ncollector),icollt(ncollector),icollb(ncollector))

    READ(DIAGFN,'(//(10X,9I10))')(Icoll(I),I=1,Ncollector)
    READ (DIAGFN,'(//(10X,9I10))')(IcollT(I),I=1,Ncollector)
    READ (DIAGFN,'(//(10X,9I10))')(IcollB(I),I=1,Ncollector)

!Collector: ON/OFF turn on a collector - if any fish gets in the collector it is captured and removed from system: flag is "3"
!NCOLLECTOR: # of segments where there are collectors
!ICOLL: seg # of collector
!ICOLLT, ICOLLB: top and bottom of collector in layer #

! Particle Behavior Modifocation
    READ(DIAGFN,'(//7x,a3,3f10.0)')PartMod, WMAX, ZMIN,ZMAX

! PartMod: ON/OFF modify particle behavior ?
! WMAX: maximum vertical velocity in m/s (negaive is upward)
! ZMIN: minimum depth to turn on behavior maodification
! ZMAX: maximum depth to to attain 100% of WMAX (note exponential transition between ZMIN and ZMAX)
   ENDIF
   
 end if

!      NULLFIELDWQ      = .FALSE.     ! If NULLFIELDWQ = .TRUE., actual WQ conditions are replaced with desired (Null) conditions
!                                     ! If NULLFIELDWQ = .FALSE., actual Water Quality conditions are modeled
!      WQNULLLINEAR     = .TRUE.      ! If WQNULLLINEAR = .TRUE., use Linear Variation in desired (Null) WQ values
!                                     ! If WQNULLLINEAR = .FALSE., use Parabolic Variation in desired (Null) WQ values
!
!      VARYTEMP         = 1           ! VARYTEMP = 1 Null Temp Values vary vertically; = 2 Null Temp Values vary horizontally
!      !-Relevant if VARYTEMP = 1     !          = 0 No Null Temp Values are implemented, while still allowing Null DO Values
!        TOPTEMP        = 24          ! TOPTEMP = Desired Temperature at the Surface Layer (KTWBF or KTWB(WB))
!        MIDKTEMP       = 15          ! MIDKTEMP = Desired Temp at selected middle Layer (MIDKT) if using Parabolic Variation
!        MIDKT          = 20          ! MIDKT = Layer # where the Temperature = MIDKTEMP if using Parabolic Variation
!        BOTTEMP        = 5           ! BOTTEMP = Desired Temperature at the absolute Water Body Bottom (BOTK)
!        BOTK           = 50          ! BOTK = Layer # of the absolute Water Body Bottom
!      !-Relevant if VARYTEMP = 2
!        LFTTEMP        = 30          ! LFTTEMP = Desired Temperature at the Upstream End of Branch
!        MIDITEMP       = 20          ! MIDITEMP = Desired Temp at selected middle Segment (MIDIT) if using Parabolic Variation
!        MIDIT          = 130         ! MIDIT = Segment # where the Temperature = MIDITEMP if using Parabolic Variation
!        RGTTEMP        = 5           ! RGTTEMP = Desired Temperature at the Downstream End of Branch
!
!      VARYDO           = 1           ! VARYDO = 1 Null DO Values vary vertically; = 2 Null DO Values vary horizontally
!      !-Relevant if VARYDO = 1       !        = 0 No Null DO Values are implemented, while still allowing Null Temp Values
!        TOPDO          = 10          ! TOPDO = Desired Dissolved Oxygen at the Surface Layer (KTWBF or KTWB(WB))
!        MIDKDO         = 6           ! MIDKDO = Desired DO at selected middle Layer (MIDKD) if using Parabolic Variation
!        MIDKD          = 25          ! MIDKD = Layer # where the Dissolved Oxygen = MIDKDO if using Parabolic Variation
!        BOTDO          = 2           ! BOTDO = Desired Dissolved Oxygen at the absolute Water Body Bottom (BOTK)
!      !-Relevant if VARYDO = 2
!        LFTDO          = 8           ! LFTDO = Desired Dissolved Oxygen at the Upstream End of Branch
!        MIDIDO         = 6           ! MIDIDO = Desired DO at selected middle Segment (MIDID) if using Parabolic Variation
!        MIDID          = 135         ! MIDID = Segment # where the Dissolved Oxygen = MIDIDO if using Parabolic Variation
!        RGTDO          = 4           ! RGTDO = Desired Dissolved Oxygen at the Downstream End of Branch
!______________________________________________________
!      VVELCAP          = 4E-4        ! Any Vert Vel that exceeds this value will have its TECPLOT VECTOR truncated back to this value
!                                     !   Suggested Value = 4E-4
!      NULLFIELDFF      = .FALSE.     ! If NULLFIELDFF = .TRUE., actual Flow conditions are replaced w/desired (Null) conditions
!                                     ! If NULLFIELDFF = .FALSE., actual Flow conditions are modeled
!      FFNULLLINEAR     = .TRUE.      ! If FFNULLLINEAR = .TRUE., use Linear Variation in desired (Null) Flow values
!                                     ! If FFNULLLINEAR = .FALSE., use Parabolic Variation in desired (Null) Flow values
!
!      VARYHVEL         = 2           ! VARYHVEL = 1 Null Horz Vel Values vary vertically; = 2 Null Horz Vel Values vary horizontally
!      !-Relevant if VARYHVEL = 1     !          = 0 No Null Horz Vel Values are implemented, while still allowing Null Vert Vel Values
!        TOPHVEL        = 3E-4        ! TOPHVEL = Desired Horz Vel at the Surface Layer (KTWBF or KTWB(WB))
!        MIDKHVEL       = 1E-4        ! MIDKHVEL = Desired Horz Vel at selected middle Layer (MIDKH) if using Parabolic Variation
!        MIDKH          = 20          ! MIDKH = Layer # where the Horz Vel = MIDKHVEL if using Parabolic Variation
!        BOTHVEL        = 6E-6        ! BOTHVEL = Desired Horz Vel at the absolute Water Body Bottom (BOTKK)
!        BOTKK          = 50          ! BOTKK = Layer # of the absolute Water Body Bottom
!      !-Relevant if VARYHVEL = 2
!        LFTHVEL        = 5E-1        ! LFTHVEL = Desired Horz Vel at the Upstream End of Branch
!        MIDIHVEL       = 8E-5        ! MIDIHVEL = Desired Horz Vel at selected middle Segment (MIDIH) if using Parabolic Variation
!        MIDIH          = 130         ! MIDIH = Segment # where the Horz Vel = MIDIHVEL if using Parabolic Variation
!        RGTHVEL        = 1E-3        ! RGTHVEL = Desired Horz Vel at the Downstream End of Branch
!
!      VARYVVEL         = 1           ! VARYVVEL = 1 Null Vert Vel Values vary vertically; = 2 Null Vert Vel Values vary horizontally
!      !-Relevant if VARYVVEL = 1     !          = 0 No Null Vert Vel Values are implemented, while still allowing Null Horz Vel Values
!        TOPVVEL        = 1E-9        ! TOPVVEL = Desired Vert Vel at the Surface Layer (KTWBF or KTWB(WB))
!        MIDKVVEL       = 2E-5        ! MIDKVVEL = Desired Vert Vel at selected middle Layer (MIDKV) if using Parabolic Variation
!        MIDKV          = 25          ! MIDKV = Layer # where the Vert Vel = MIDKVVEL if using Parabolic Variation
!        BOTVVEL        = 1E-9        ! BOTVVEL = Desired Vert Vel at the absolute Water Body Bottom (BOTKK)
!      !-Relevant if VARYVVEL = 2
!!        LFTVVEL        = 3E-4        ! LFTVVEL = Desired Vert Vel at the Upstream End of Branch
!        MIDIVVEL       = 2E-5        ! MIDIVVEL = Desired Vert Vel at selected middle Segment (MIDIV) if using Parabolic Variation
!        MIDIV          = 135         ! MIDIV = Segment # where the Vert Vel = MIDIVVEL if using Parabolic Variation
!        RGTVVEL        = 1E-5        ! RGTVVEL = Desired Vert Vel at the Downstream End of Branch
!______________________________________________________
!      STIMULIRULES     = .FALSE.      ! If .TRUE., stimuli-response rules (excluding passive transport) contribute to fish movement
!           ! TRUE                         ! If .FALSE., stimuli-resp. rules (excl. passive transp) DO NOT contribute to fish movement
!      VELOCITYRULES    = .FALSE.      ! If .TRUE., velocity stimuli-response rules contribute to the movement of fish
!           ! TRUE                         ! If .FALSE., velocity stimuli-response rules DO NOT contribute to the movement of fish
!      TEMPRULES        = .FALSE.      ! If .TRUE., temperature stimuli-response rules contribute to the movement of fish
!           ! TRUE                          ! If .FALSE., temperature stimuli-response rules DO NOT contribute to the movement of fish
!      DORULES          = .FALSE.      ! If .TRUE., dissolved oxygen stimuli-response rules contribute to the movement of fish
!           ! TRUE                          ! If .FALSE., dissolved oxygen stimuli-response rules DO NOT contribute to the movement of fish
!      RANDOMIZATION    = .TRUE.      ! If .TRUE., random displacement terms (i.e. RDX,RDY,RDZ) are calculated
!                                     ! If .FALSE., random displacement terms (i.e. RDX,RDY,RDZ) are equal to 0.0
!______________________________________________________
!      PASSIVETRANSPORT = .TRUE.      ! If .TRUE., passive transport contributes to the movement of fish
!                                     ! If .FALSE., passive transport DOES NOT contribute to the movement of fish
!      SCHOOLING        = .FALSE.     ! If .TRUE., fish school forming function is activated
!                                     ! If .FALSE., fish school forming function is DEACTIVATED
!______________________________________________________
!      DEBUG            = .FALSE.     ! If .TRUE., user informed when important steps are completed successfully
!                                     ! If .FALSE., no extra output generated that may assist in debugging Numerical Fish Surrogate
!      WBSKIP           = .TRUE.      ! If .TRUE., output files will be generated for Water Body "WBRUN" only
!                                     ! If .FALSE., output files will be generated for all Water Bodies
!      WBRUN            = 2           ! If WBSKIP = .TRUE., output files will be generated for this Water Body only
!      LINEAR           = .FALSE.     ! If .TRUE., velocity interpolation scheme ==> linear interpolation
!                                     ! If .FALSE., velocity interpolation scheme ==> 3rd Order Newton Interpolating Polynomial
!      PREVENTBRCHSWITCH= .TRUE.      ! If .TRUE., prevents fish from moving upstream (only) into another branch
!                                     ! If .FALSE., allows fish to move upstream and/or downstream into other branches
!______________________________________________________
!      ASPRATIO         = 1.0         ! Used for scaling Vertical Velocity output for TecPlot purposes only
!      SHOWSKY          = .TRUE.      ! If .TRUE., cells above the water surface will be colored according to the time of day
!                                     ! If .FALSE., cells above the water surface will be blank (i.e. WHITE)
!      SKYNIGHT         = 20.5        ! SKYNIGHT = The time of day when night begins; after which sky is BLACK (military time)
!      SKYDAWN          = 6.5         ! SKYDAWN = The time of day when morning begins; after which sky is YELLOW (military time)
!      SKYDAY           = 9.5         ! SKYDAY = The time of day when 'day' begins; after which sky is BLUE (military time)
!      SKYDUSK          = 17          ! SKYDUSK = The time of day when evening begins; after which sky is ORANGE (military time)
!______________________________________________________
!      UNBP             = 6           ! The Branch, NBP, at who's upstream end fish are collected
!                                     !   Fish are already collected and tallied when moving upstream
!                                     !   into the dam (i.e., US(JBP))
!      DNBP             = 6           ! The Branch, NBP, at who's downstream end fish are collected
!                                     !   This must be the Branch which has no Branch downstream of it

    close(DIAGFN)

1000 Format(//14I10)
1001 Format(//14f10.0)
1002 Format(//14(7x,a3))


End Subroutine Read_Fish_Data

Subroutine findbranch(seg)

     Use Fishy; USE GLOBAL
     IMPLICIT NONE
     
    real seg
    integer iseg,j

     iseg=int(seg)
      do j=1,NBR
      if(iseg.ge.us(j).and.iseg.le.ds(j))exit
      end do
      fnbp=j
      return

End Subroutine findbranch

Subroutine Part_transport
! Compute passive particle transport including sedimentation
    
    Use Fishy
    Use TRANS
    USE GLOBAL
    Use GEOMC
    IMPLICIT NONE
    
    REAL :: DXMIN, DZMIN1,DZMAX1, COSTHETA, SINTHETA, DX1, DX2, DXAVG, SK, RZ, RX, DZ1, DZ2, DZAVG, r1, r2
    REAL :: DISPX, DISPZ, VEL,WPART,XAREA
    Data DXMIN,DZMIN1,DZMAX1 /1.0,0.2,5.0/   ! SW 2/01/01   ****change this since dzmin and dzmax are set in main file

        FXLOC = FXLOC + (FXVEL(5)) * (nfsfreq)*24.*3600.    ! New updated Part X-Location
        FYLOC = FYLOC + (FYVEL) * (nfsfreq)*24.*3600.    ! New updated Part Y-Location
        FZLOC = FZLOC + (FZVEL(5) + SEDVEL) * (nfsfreq)*24.*3600.    ! New updated Part Z-Location

if(partmod.eq.' ON')then
  if(depthm(fkmp,fimp).ge.zmin.and.depthm(fkmp,fimp).le.zmax)then
  sk=-log(0.001)/(zmin-zmax)
  WPART=WMAX*(1.0-exp(-sk*(depthm(fkmp,fimp)-zmax)))
  FZLOC=FZLOC+WPART*nfsfreq*24.*3600.
  elseif(depthm(fkmp,fimp).gt.zmax)then
  FZLOC=FZLOC+WMAX*nfsfreq*24.*3600.
  end if
end if

        vel=sqrt(fxvel(5)*fxvel(5)+fzvel(5)*fzvel(5))

    ! Compute "Dispersion Processes"
        if(vel.gt.1.0e-09)then
        costheta=fxvel(5)/vel
        sintheta=fzvel(5)/vel
!        Dispx=alphax*abs(fxvel(5))
!        Dispz=alphaz*abs(fzvel(5))

! average to segment centers - should be interpolated though   ! SW 2/01/01
       DX1=DX(fkmp,fimp)
       DX2=DX(fkmp,fimp-1)
       if(dx1.le.0.0)dx1=dxmin
       if(dx2.le.0.0)dx2=dxmin
       DXAVG=0.5*(dx1+dx2)

       DZ1=DZ(fkmp,fimp)
       if(fkmp.le.KTWBF)then
       DZ2=DZ(fkmp,fimp)
       else
       DZ2=DZ(fkmp-1,fimp)
       endif
! Note during density inversions DZ is set to DZMAX1=1000. This leads to incredible
! variations in RZ - constrain to DZ=5.0
       if(dz1.gt.10.)dz1=dzmax1
       if(dz2.gt.10.)dz2=dzmax1

       if(dz1.le.0.0)dz1=dzmin1
       if(dz2.le.0.0)dz2=dzmin1
       DZAVG=0.5*(dz1+dz2)

!if(dzavg.gt.10)then
!write(DIAGFN,*)'DZAVG>10 JDAY:',jday,'dz1,dz2,fkmp,fimp,dz(fkmp,fimp),dz(fkmp-1,fimp)'
!write(DIAGFN,*)dz1,dz2,fkmp,fimp,dz(fkmp,fimp),dz(fkmp-1,fimp)
!end if
!if(dxavg.gt.20)then
!write(DIAGFN,*)'DXAVG>10 JDAY:',jday,'dx1,dx2,fkmp,fimp,dx(fkmp,fimp),dx(fkmp-1,fimp)'
!write(DIAGFN,*)dx1,dx2,fkmp,fimp,dx(fkmp,fimp),dx(fkmp,fimp-1)
!end if


        Dispx=DXAVG*alphax            ! Note these should be interpolated rather than using nearest cell #
        Dispz=DZAVG*alphaz

! COMPUTE DZ and DX by interpolating

        call random(seed,r1)
        call random(seed,r2)

!        r1=gasdev(seed)
!        r2=gasdev(seed)

        r1=(r1-0.5)*2.
        r2=(r2-0.5)*2.

        RX=SQRT(6.0*Dispx*nfsfreq)*(r1*costheta-r2*sintheta)
        RZ=SQRT(6.0*Dispz*nfsfreq)*(r1*sintheta+r2*costheta)

! constain random component to segment length and cell layer height
        if(abs(rz).gt.h(KTWBF,fjr))rz=sign(1.0,rz)*h(KTWBF,fjr)   !.gt.h(fkmp,KTWBF))rz=sign(1.0,rz)*h(fkmp,KTWBF)    SW 7/1/2017
        if(abs(rx).gt.dlx(fimp))rx=sign(1.0,rx)*dlx(fimp)

!if(rz.gt.2..or.rx.gt.1000.)then
!write(DIAGFN,*)'rx,rz,dispz,nfsfreq,r1,sintheta,r2,costheta,dzavg,alphaz,dxavg,alphax,dz1,dz2,dx1,dx2,vel,fxvel(5),fzvel(5)'
!write(DIAGFN,*)rx,rz,dispz,nfsfreq,r1,sintheta,r2,costheta,dzavg,alphaz,dxavg,alphax,dz1,dz2,dx1,dx2,vel,fxvel(5),fzvel(5)
!end if

        !RX=SQRT(24.0*Dispx*nfsfreq)*(r1-0.5)
        !RZ=SQRT(24.0*Dispz*nfsfreq)*(r2-0.5)

        else
        RX=0.0
        RZ=0.0

        end if

        FXLOC = FXLOC + RX   ! New updated Part X-Location
        FYLOC = FYLOC + RX   ! New updated Part Y-Location
        FZLOC = FZLOC + RZ    ! New updated Part Z-Location

!if(rx.eq.0.0)then
!write(DIAGFN,*)'RX=0.0:rx,rz,dispx,r1,sintheta,r2,costheta,vel,fxvel(5)'  ! debug
!write(DIAGFN,*)rx,rz,dispx,r1,sintheta,r2,costheta,vel,fxvel(5)  ! debug
!end if

    return

End Subroutine Part_transport

    Function Gasdev(idum3)

! From S. Li, PSU, PArticle Transport Algorithm
      data iset /0/
      if (iset.eq.0)then
 1       v1=2.*ran(idum3)-1
	 v2=2.*ran(idum3)-1
	 r=v1**2+v2**2
	 if(r.ge.1)goto 1
	 fac=sqrt(-2.*log(r)/r)
	 gset=v1*fac
	 gasdev=v2*fac
	 iset=1
      else
	 gasdev=gset
	 iset=0
      endif
      return

    End Function Gasdev

!**********************************
   Subroutine lateral_velocity
! SW 2/01/02

   Use FISHY
   USE GLOBAL
   Use GEOMC
   
   IMPLICIT NONE
   REAL :: XAREA
   
! Concept all QSS from each cell will be treated as a lateral withdrawal - each withdrawal will
! be assigned a RHS or LHS looking downstream location; lateral velocity origin is the
! segment/cell center. Velocities to the RHS are + and those to the LHS are -
!  


! Assume at first that all withdrawals and inputs are on RHS (hence inflow would generate - and outflow + velocities)


    if(fkmp.le.KTWB(fjr))then
    xarea=h1(KTWB(FJR),fimp)*dlx(fimp)
    fyvel=-qss(KTWB(fjr),fimp)/xarea
    else
    xarea=h1(fkmp,fIMP)*dlx(fimp)
    fyvel=-qss(fkmp,fimp)/xarea
    end if

    return
    end
!************************************