!***********************************************************************
!***********************************************************************
!***********************************************************************
!***                                                                ****
!***        N U M E R I C A L   PARTICLE/FISH   S U R R O G A T E   ****
!***                                                                ****
!***********************************************************************
!***********************************************************************
!***********************************************************************
! Andy Goodwin F77 Version 1/2001
! Scott Wells F90 Version + Enhancements for Random Water Movement + Bug Fixes/Link to W2V3.1  1/2001
! Scott Wells Updates to latest version V3.7 5/1/2015 
! Scott Wells Updates and particle tracking V4.1 8/1/2017 8/2018
! Module Definitions

Module Fishy
      Integer, PARAMETER :: FPARA=16,FFP=4,WQP=2    ! The number of fish parameters. The number of Flow Field Parameters evaluated &
                                     !   the number of Water Quality Parameters evaluated at nodes for use in the NFS
      Integer NFISH,NUMGILLNETS,NUMACOUSTICS,JJ, JBP, KLAST,NL,NR
      DOUBLE PRECISION VALUE,GRADXVEL,GRADZVEL,GRADTEMP,GRADDO
      DOUBLE PRECISION XGRADHV,ZGRADHV,XGRADVV,ZGRADVV
      DOUBLE PRECISION XGRADDO,ZGRADDO,XGRADTP,ZGRADTP
      DOUBLE PRECISION FXVEL,FZVEL,FISHTEMP,FISHDO,TEMPORARYGRAD
      DOUBLE PRECISION FXHVCOUNT,FXVVCOUNT,FXDOCOUNT,FXTPCOUNT,FXRDCOUNT,TOTFXWGT
      DOUBLE PRECISION FZHVCOUNT,FZVVCOUNT,FZDOCOUNT,FZTPCOUNT,FZRDCOUNT,TOTFZWGT
      DOUBLE PRECISION FX00COUNT,FZ00COUNT,nfsfreq,RUNDIFF
      REAL ::   wmax,zmin,zmax  ! SW 2/16/01
      REAL            FXLOC,FYLOC,FZLOC,XREFL,YREFL,ZBOTREFL,ZSURREFL
      REAL            UFISH,VFISH,WFISH,FYVEL,ALPHAX,ALPHAZ   
      REAL            OUTFREQ,JDAYDIFF,TSETP2MULT 
      REAL            FROMKTBOT,RRR,RMAT,RDX,RDY,RDZ,FYLOCTMP
      REAL            FSIZE,FAGE,LJDAY
      REAL            OUTFREQJAN,OUTFREQFEB,OUTFREQMAR,OUTFREQAPR
      REAL            OUTFREQMAY,OUTFREQJUN,OUTFREQJUL,OUTFREQAUG
      REAL            OUTFREQSEP,OUTFREQOCT,OUTFREQNOV,OUTFREQDEC,OUTFREQP   ! SW 7/1/2017
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
      REAL            MIDIDO,MIDID,TOPHVEL,BOTHVEL,MIDKHVEL,MIDKH
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
      REAL            TEMPTHRES,TSTEP1,TSTEP1MULT,TSTEP2,TSTEP2MULT,DOTHRES2,DOSTEP1MULT       !,DELAYDATE

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

      LOGICAL         DEBUG,LINEAR,WBSKIP,SHOWSKY
      LOGICAL         particle,line,HITSTICKBOTTOM,HITSTICKSIDE,HIST_V,HIST_T,HIST_D  ! SW 4/30/2018
      REAL         :: VEL_INT,VEL_TOP,TEMP_INT,TEMP_TOP,D_INT,D_TOP
      INTEGER      :: NUMCLASS, NMONITOR  
      character*3     PARTON,collector,partmod,DXTHEORY        ! SW 2/16/01 7/25/2017 4/30/2018

      real, allocatable, dimension (:,:,:) :: flowfield, wqfield,sndcatch,nodes,lastfield
      integer, allocatable, dimension (:,:,:) :: netcatch
      integer, allocatable, dimension (:,:) :: rimpbr,limpbr
      real, allocatable, dimension (:,:) :: fishes
      integer, allocatable, dimension  (:) :: SNAGCOUNT,SNDCOUNT,nbrf
      integer, allocatable, dimension (:,:,:) :: corners
      integer, allocatable, dimension (:,:) :: ndinfo,MONITORONOFF
      integer, allocatable, dimension  (:) :: ifish,ifisht,ifishb,GROUP,FNEXT  ! SW 1/16/01
      integer, allocatable, dimension  (:) :: icoll,icollt,icollb,IMONITOR  ! SW 2/16/01
      REAL, ALLOCATABLE, DIMENSION (:) :: DELAYDATE,fxloci,fzloci,sedvel,XSHIFTMONITOR
      real, allocatable, dimension (:,:,:) :: BRCHFISH
      
      REAL, ALLOCATABLE, DIMENSION (:)   :: v_tot,v_cnt, d_tot, d_cnt, t_tot, t_cnt, d_avg, t_avg, v_avg, d_sum, t_sum, v_sum, sumvolt  ! SW 4/30/2018
      REAL, ALLOCATABLE, DIMENSION (:,:) :: v_class, t_class, d_class, TMONITOR
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: NEXT_BRANCH
      
      INTEGER :: DIAGFN=8001, DATADEBUGFN=8002, BARCHRTXFN=8003,BARCHRTZFN=8004,FINALFN=8005,INITIALFN=8006, NMONITORS

End Module Fishy


!  This Subroutine Calls the Following Subroutines:
!      FIMPBR,GRIDPLOT,WHATJR,SPLINE,RANDOM
!      FINDNEWBR,FISHPLOT,INTERCONST,INTERFLOWF


      SUBROUTINE FISH

     USE SURFHE; Use Fishy; Use GDAYC;  Use SCREENC; Use GEOMC; USE GLOBAL;   Use MAIN, only: IWD, FISH_PARTICLE_EXIST, KBI 
     IMPLICIT NONE
     
     REAL :: DZ
     INTEGER :: JF,KK,N, IW, GROUPLAST

  IF (NIT.EQ.0) THEN             ! The very first time the Subroutine FISH is called

!***** Open 'Numerical Fish Surrogate' files
      Call Read_Fish_Data
      IF(.NOT.FISH_PARTICLE_EXIST)RETURN  !STOP ALL PROCESSING    
      OPEN (DIAGFN,FILE='DIAGNOSTICS.OUT',STATUS='UNKNOWN')                         !FISH
      OPEN (DATADEBUGFN,FILE='DATADEBUG.OUT',STATUS='UNKNOWN')                           !FISH
      OPEN (FINALFN,file='finalparticle.csv',status='unknown')                             !FISH
      OPEN (INITIALFN,file='initialparticle.csv',status='unknown')

! Allocate arrays
! Read input data files

  Allocate(FISHES(NFISH,FPARA),WQFIELD(KMX,IMX,WQP),FLOWFIELD(KMX,IMX,FFP),NETCATCH(NUMGILLNETS,NFISH,5),&
           SNDCATCH(NUMACOUSTICS,3*NFISH,2),RIMPBR(IMX,2),LIMPBR(IMX,2), NODES(KMX,IMX,2),&
           NDINFO(IMX*KMX+IMX+KMX-1,4),CORNERS(NBR,IMX*KMX+IMX+KMX-1,4),LASTFIELD(KMX,IMX,2),&
           BRCHFISH(NBR,NFISH,FPARA),NBRF(NBR))
  ALLOCATE(NEXT_BRANCH(NBR),FNEXT(NBR))
  flowfield=0.0
  wqfield=0.0
  corners=0
  ndinfo=0
  nodes=0.0
  netcatch=0
  sndcatch=0.0
  lastfield=0.0
  brchfish=0.0
  nbrf=0
  fishes=0.0
  wbskip=.false.
  showsky=.false.
        CALL FIMPBR                  ! Establishes data set of what segments are in each Branch
        SEED           = 92          ! Seed for the Random Number Generator


!------------------------------------------ PARTICLE DATA INITIALIZATION ------------------------------------------

!Multiple Particles
        NGRPFH = 0                                   ! NGRPFH  = The number of particles/fish being initiated
        DO jj=1,nfishseg
          DO FK=ifisht(jj),ifishb(jj)                 !KTWB(2),(KB(FI)-2),1   ! Fish will be placed at each one of these Layer nodes
            DO FATPT=1,nfishpcel                     ! The Number of Fish that will be placed at each drop location
              NGRPFH = NGRPFH + 1                        ! NGRPFH  = The number of particles/fish being initiated
              FISHES(NGRPFH,1) = ifish(jj)    !FI        ! FIMP    = Initial segment IMP where particle/fish is released
              FISHES(NGRPFH,3) = FK                      ! FKMP    = Initial layer KMP where particle/fish is released
            end do
           end do
        end do
!_____________________________________________________
!Data for All Fish
        GROUPLAST=1
        N=NFISHPCEL
        DO FN=1,NFISH
         IF(.NOT.LINE)THEN
          FISHES(FN,4)  = FZLOCI(FN)   !0.0   ! FZLOC   = Location of fish within layer KMP from top side
         ELSE
           if(N==NFISHPCEL .OR. GROUPLAST /= GROUP(FN))then
            FISHES(FN,4)  = 0.           ! FZLOC   = assumed to be zero for LINE
            fimp=fishes(fn,1)
            call whatjr
            dz=h(fishes(fn,3),fjr)/nfishpcel
            GROUPLAST=GROUP(FN)
            N=1
           else
            fishes(fn,4)=fishes(fn-1,4)+dz
            N=N+1
           endif 
         ENDIF
          FISHES(FN,2)  = FXLOCI(FN)   !0.0   ! FXLOC   = Location of fish within segment IMP from upstream side
          FISHES(FN,5)  = B(INT(FISHES(FN,3)),INT(FISHES(FN,1)))*0.5
                                     ! FYLOC   = Lateral fish release location (from left bank in plan view)
          !                     Looking down on a segment 
          !                                downstream
          !        y=0                        y=B/2                       y=B
          !         +---------------------------+--------------------------+
          !         |                                                      |
          !         |                                                      |
          !  L      |                           X                          |       R     ! location of particle in lateral
          !         |                                                      |          
          !         |                                                      |          
          !         |                                                      |
          !         +----------------------------+-------------------------+
          !                                 upstream
          
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

        CALL NEXTBRANCH    ! SW 8/2018 ARE BRANCHES CONNECTED BY SPILLWAYS/GATES/PUMPS TO TRANSFER PARTICLES TO OTHER PARTS OF THE SYSTEM
  END IF

!*********************************************************************************************
!******************************  NO MORE INPUT BELOW THIS POINT  *****************************
!*********************************************************************************************


!Initialize Parameters Used to Evaluate Implementation of the Stimuli-Response Rules,
!  Prepare Text for Output Purposes, and Output Input Data to DIAGNOSTICS.OUT File


      IF (NIT.EQ.0) THEN
        VVELCAP          = 4.E-4    ! Any Vert Vel that exceeds this value will have its TECPLOT VECTOR truncated back to this value
                                   !   Suggested Value = 4E-4
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


        if(particle)then
        txtparticle='  ON'
        else
        txtparticle=' OFF'
        end if


! initial fish output
    write(INITIALFN,'(a291)')'Part#,Seg#,XLocationwithinSegmentfromUpstreamSide(m),Layer#,VerticalDistfromTop(m),LateralDistfromLeftBank,Branch#,ParticleInModel(=0),JDAYleftsystem,DetentionTime(days),RemovalMechanism,SedVelocity(m/d),DateStart'
    do jf=1,nfish     
    write(INITIALFN,'(i7,",",1x,10(f10.3,","),f8.4,",",f12.4)')jf,fishes(jf,1),fishes(jf,2),fishes(jf,3),fishes(jf,4),fishes(jf,5),fishes(jf,6),fishes(jf,12),fishes(jf,14),fishes(jf,15),fishes(jf,16),sedvel(jf),delaydate(jf)
    end do

    close(INITIALFN)

      END IF


!Specification of the Time Stepping used in Running the NFS (W2 may, at times, run at timesteps of down to 2 minutes!)

! COMPUTE DYNAMIC NFSFREQ SW 1/15/01

      nfsfreq=dlt/86400.    

      IF (NIT.EQ.0) LRUNDAY = JDAY                            ! LRUNDAY is used in determining NFS run frequency
      RUNDIFF = JDAY - LRUNDAY                                ! RUNDIFF is used in determining NFS run frequency


!Specification of when and how often to output information/data

      IF (NIT.EQ.0) LJDAY = JDAY                              ! LJDAY is used for TecPlot output frequency purposes
      JDAYDIFF = JDAY - LJDAY                                 ! JDAYDIFF is used for TecPlot output frequency purposes
      OUTFREQ=OUTFREQP
      
      IF (JDAYDIFF.GE.OUTFREQ) THEN     ! Set FCOUNT = 1, so information will be outputted
        FCOUNT = 1                                            ! Information outputted iff FCOUNT = 1
        LJDAY = JDAY
      ELSE
        FCOUNT = 0                                            ! Set FCOUNT = 0, so information won't be outputted
      END IF


!Calculate Time of Day/Night

      MILTIME =  (JDAY-REAL(INT(JDAY)))*24.       !(JDAY-REAL(JDAYG))*24                         ! Military Time - MILTIME used for NFS Calculations
      HR = MILTIME                                            ! Convert MILTIME to hour of day - HR used for TecPlot
      AMPM = 'am'                                             ! Is HR am or pm?
      IF (HR.GE.12) AMPM = 'pm'                               ! Is HR am or pm?
      IF (HR.GE.13) HR = HR - 12                              ! Value of HR after 12 noon
      IF (HR.LT.1) HR = HR + 12                               ! Value of HR if time is between 12 midnight and 1am


!Call the Subroutine GRIDPLOT to set up the FE Grid for TecPlot and to Calc node information for Traffic Rules

      IF ((FCOUNT.EQ.1).OR.(NIT.EQ.0)) CALL GRIDPLOT          ! This prepares the FE Grid for output display and
                                                              !   calculates the Flow and WQ information at nodes

!Load Fish Information and Begin NFS Logic

      !IF (JDAY .LT. DELAYDATE) GOTO 28
      
      CALL INTERCONST                                         ! Subroutine to interpolate Constituent values
      CALL INTERFLOWF                                         ! Subroutine to interpolate Flow Field values
 
      DO 131 FN=1,NFISH              ! This is the only loop where statements are not indented
      IF (JDAY .LT. DELAYDATE(FN)) CYCLE
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
        ! SW 8/2018 Correct for layer subtraction
        IF(KTWBF > FKMP)THEN
            IF(KB(FIMP)==KBI(FIMP))THEN
            FKMP=KTWBF
            FISHES(FN,3)=FKMP    ! Layer KMP where fish is located
            FZLOC=-(H(FKMP-1,FJR)-FZLOC)
            FISHES(FN,4)=FZLOC
            ENDIF
        ENDIF
        
            
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
CALL PART_TRANSPORT



!***************************************************************
!  C H E C K   F O R   B O U N D A R Y   V I O L A T I O N S   *
!***************************************************************


!Check for Boundary Violations: Lateral Direction Check (1 of 2)

  if(fyvel.eq.0.0)then   ! SW 2/01/01 logic for removing particles laterally also
        IF (FYLOC.LT.0) THEN    ! reflect particles
            IF(HITSTICKSIDE)THEN
                ISWITCH = 1
                fishes(fn,14)=JDAY               ! Time fish left system
                fishes(fn,15)=JDAY-DELAYDATE(FN)     ! Detention time of fish in system
                fishes(fn,12)=1.0                ! FISH LEAVES SYSTEM
                fishes(fn,16)=6.0                ! Particle leaves by hitting side wall and sticking-LHS
            ELSE
           FYLOC = B(FKMP,FIMP)*YREFL 
           ENDIF
        ELSE IF (FYLOC.GT.B(FKMP,FIMP)) THEN
            IF(HITSTICKSIDE)THEN
                ISWITCH = 1
                fishes(fn,14)=JDAY               ! Time fish left system
                fishes(fn,15)=JDAY-DELAYDATE(FN)     ! Detention time of fish in system
                fishes(fn,12)=1.0                ! FISH LEAVES SYSTEM
                fishes(fn,16)=7.0                ! Particle leaves by hitting side wall and sticking-RHS
            ELSE
           FYLOC = B(FKMP,FIMP)*(1-YREFL) 
           ENDIF
        END IF
  else
        IF (FYLOC.LT.0.and.fyvel.lt.0.0.AND.LIMPBR(FIMP,1)==0) THEN  ! 
            ! check for a withdrawal
            DO IW=1,NWD
              IF(FIMP==IWD(IW))THEN    ! remove particle through withdrawal
                ISWITCH = 1
                fishes(fn,14)=JDAY               ! Time fish left system
                fishes(fn,15)=JDAY-DELAYDATE(FN)     ! Detention time of fish in system
                fishes(fn,12)=1.0                ! FISH LEAVES SYSTEM
                fishes(fn,16)=1.0                ! Signifies the lateral removal by withdrawal
                WRITE(DIAGFN,*) 'Withdrawal: Particle Leaves System on JDAY:',jday, 'Data:',fishes(fn,:)
                GOTO 20
              ENDIF     
            ENDDO     
            IF(HITSTICKSIDE)THEN
                 ISWITCH = 1
                fishes(fn,14)=JDAY               ! Time fish left system
                fishes(fn,15)=JDAY-DELAYDATE(FN)     ! Detention time of fish in system
                fishes(fn,12)=1.0                ! FISH LEAVES SYSTEM
                fishes(fn,16)=6.0                ! Particle leaves by hitting side wall and sticking-LHS
                WRITE(DIAGFN,*) 'Hit Stick LHSide: Particle sticks to side wall on JDAY:',jday, 'Data:',fishes(fn,:)
            ELSE
            FYLOC = B(FKMP,FIMP)*YREFL     ! reflect
            ENDIF
        ELSE IF(FYLOC.GT.B(FKMP,FIMP).and.fyvel.gt.0.0.AND.RIMPBR(FIMP,1)>0) THEN   
            WRITE(DIAGFN,*) 'LateralRIGHT: Particle move to new branch on JDAY:',jday,'FYVEL=',fyvel,'OLD I:',FIMP,'NEW I:', RIMPBR(FIMP,1),' OLD branch:',FNBP,' NEW branch:',RIMPBR(FIMP,2)
            FIMP=RIMPBR(FIMP,1)
            FYLOC= B(FKMP,FIMP)*0.5   ! PLACE IN MIDDLE LATERALLY
            FXLOC= DLX(FIMP)*0.9      ! PLACE IN NEAR END OF SEGMENT
           
           !ISWITCH = 1
           !! SW 2/01/01 Track timing of fish movement from system
           !fishes(fn,14)=JDAY               ! Time fish left system
           !fishes(fn,15)=JDAY-DELAYDATE     ! Detention time of fish in system
           !fishes(fn,16)=1.0                ! Signifies the lateral removal
           !fishes(fn,12)=1.0                ! FISH LEAVES SYSTEM
           !WRITE(DIAGFN,*) 'LateralRIGHT: Particle Leaves System on JDAY:',jday,'FYVEL=',fyvel, 'Data:',fishes(fn,:)
           ! End Section SW 2/01/01
           !GOTO 20
           ! Particle is removed       ***need to transfer to another branch - it is only removed****
        ELSE IF(FYLOC.LT.0.0.and.fyvel.Lt.0.0.AND.LIMPBR(FIMP,1)>0) THEN  
            WRITE(DIAGFN,*) 'LateralLEFT: Particle move to new branch on JDAY:',jday,'FYVEL=',fyvel,'OLD I:',FIMP,'NEW I:', LIMPBR(FIMP,1),' OLD branch:',FNBP,' NEW branch:',LIMPBR(FIMP,2)
            FIMP=LIMPBR(FIMP,1)
            FYLOC= B(FKMP,FIMP)*0.5   ! PLACE IN MIDDLE LATERALLY
            FXLOC= DLX(FIMP)*0.9      ! PLACE NEAR END OF SEGMENT                 
           !ISWITCH = 1
           !! SW 2/01/01 Track timing of fish movement from system
           !fishes(fn,14)=JDAY               ! Time fish left system
           !fishes(fn,15)=JDAY-DELAYDATE     ! Detention time of fish in system
           !fishes(fn,16)=1.0                ! Signifies the lateral removal
           !fishes(fn,12)=1.0                ! FISH LEAVES SYSTEM
           !WRITE(DIAGFN,*) 'LateralLEFT: Particle Leaves System on JDAY:',jday,'FYVEL=',fyvel, 'Data:',fishes(fn,:)
           ! End Section SW 2/01/01
           !GOTO 20
  
        ELSEIF (FYLOC.GT.B(FKMP,FIMP).and.RIMPBR(FIMP,1)==0.0) THEN   ! reflect 
            IF(HITSTICKSIDE)THEN
                ISWITCH = 1
                fishes(fn,14)=JDAY               ! Time fish left system
                fishes(fn,15)=JDAY-DELAYDATE(FN)     ! Detention time of fish in system
                fishes(fn,12)=1.0                ! FISH LEAVES SYSTEM
                fishes(fn,16)=7.0                ! Particle leaves by hitting side wall and sticking-RHS
                WRITE(DIAGFN,*) 'Hit Stick RHSide: Particle sticks to side wall on JDAY:',jday, 'Data:',fishes(fn,:)
            ELSE
                
           FYLOC = B(FKMP,FIMP)*(1-YREFL) 
            ENDIF
            
        ELSEIF (FYLOC.LT.0.0 .and.LIMPBR(FIMP,1)==0) THEN  ! reflect 
            IF(HITSTICKSIDE)THEN
                ISWITCH = 1
                fishes(fn,14)=JDAY               ! Time fish left system
                fishes(fn,15)=JDAY-DELAYDATE(FN)     ! Detention time of fish in system
                fishes(fn,12)=1.0                ! FISH LEAVES SYSTEM
                fishes(fn,16)=6.0                ! Particle leaves by hitting side wall and sticking-LHS
                WRITE(DIAGFN,*) 'Hit Stick LHSide: Particle sticks to side wall on JDAY:',jday, 'Data:',fishes(fn,:)
            ELSE
           FYLOC = B(FKMP,FIMP)*YREFL 
           ENDIF
        END IF
  end if

!Check for Boundary Violations: Horizontal Direction

   10   IF (FXLOC.GT.DLX(FIMP)) THEN
 
           IF(FIMP < DS(FNBP))THEN
                FIMP = FIMP + 1                                   
                FXLOC = FXLOC-DLX(FIMP) 
           ELSEIF(FXVEL(5) > 0.0 .AND. DHS(FNBP)==0)THEN
              ! CHECK IF A BRANCH DOWNSTREAM TO PASS THE PARTICLE TO
               IF(NEXT_BRANCH(FNBP))THEN
                WRITE(DIAGFN,*) 'At end of branch: Particle moves to next branch on JDAY:',jday,'OLD I:',FIMP,'NEW I:',FNEXT(FNBP),'Data:',fishes(fn,:)
                FIMP=FNEXT(FNBP)
                FYLOC= B(FKMP,FIMP)*0.5   ! PLACE IN MIDDLE LATERALLY
                FXLOC= DLX(FIMP)*0.25      ! PLACE IN FIRST QUARTER OF SEGMENT     **** FIX THIS ****
                ! KEEP THE SAME VERTICAL LAYER AS PREVIOUS PARTICLE  *** FIX THIS ****
               ELSE
               ! LOSE PARTCLE THROUGH DAM
                ISWITCH = 1
                fishes(fn,14)=JDAY               ! Time fish left system
                fishes(fn,15)=JDAY-DELAYDATE(FN)     ! Detention time of fish in system
                fishes(fn,12)=1.0                ! FISH LEAVES SYSTEM
                fishes(fn,16)=5.0  ! Particle leaves at downstream structure/dam/hydraulic structure
                WRITE(DIAGFN,*) 'At Dam: Particle Leaves System on JDAY:',jday, 'Data:',fishes(fn,:)
                ENDIF
           GOTO 20
            ELSEIF(FXVEL(5) > 0.0 .AND. DHS(FNBP)== -1)THEN
               ! LOSE PARTCLE to external head BC
                ISWITCH = 1
                fishes(fn,14)=JDAY               ! Time fish left system
                fishes(fn,15)=JDAY-DELAYDATE(FN)     ! Detention time of fish in system
                fishes(fn,12)=1.0                ! FISH LEAVES SYSTEM
                fishes(fn,16)=4.0  ! Particle leaves at external head BC Downstream
                WRITE(DIAGFN,*) 'At External Downstream Head Boundary: Particle Leaves System on JDAY:',jday, 'Data:',fishes(fn,:)
           GOTO 20
           ELSEIF(FXVEL(5) > 0.0 .AND. DHS(FNBP) > 0)THEN
                WRITE(DIAGFN,*) 'DHS: Particle move to new branch on JDAY:',jday,'FXVEL=',fxvel(5),'OLD I:',FIMP,'NEW I:', DHS(FNBP),' OLD branch:',FNBP
                FIMP=DHS(FNBP)
                FYLOC= B(FKMP,FIMP)*0.5   ! PLACE IN MIDDLE LATERALLY
                FXLOC= DLX(FIMP)*0.5      ! PLACE IN MIDDLE OF SEGMENT                  
               ! MOVE PARTICLE TO NEW BRANCH
                !ISWITCH = 1
                !fishes(fn,14)=JDAY               ! Time fish left system
                !fishes(fn,15)=JDAY-DELAYDATE     ! Detention time of fish in system
                !fishes(fn,12)=1.0                ! FISH LEAVES SYSTEM
                !WRITE(DIAGFN,*) 'At Dam: Particle Leaves System on JDAY:',jday, 'Data:',fishes(fn,:)
           ELSE
               ! REFLECT OFF DAM
                FXLOC=DLX(FIMP)*(1.0-XREFL)
           ENDIF
            
        ELSE IF (FXLOC.LT.0) THEN

          IF (FIMP.EQ.CUS(FNBP) .AND. UHS(FNBP)==0) THEN      ! The UpStream end of Branch JBP or UNBP
              FXLOC = DLX(FIMP)*XREFL    ! REflect partilce off UPSTREAM boundary 
              
              ELSEIF(FIMP.EQ.CUS(FNBP) .AND. UHS(FNBP)==-1 .AND. FXVEL(5) < 0.0) THEN 
                            ! LOSE PARTCLE to external head BC
                ISWITCH = 1
                fishes(fn,14)=JDAY               ! Time fish left system
                fishes(fn,15)=JDAY-DELAYDATE(FN)     ! Detention time of fish in system
                fishes(fn,12)=1.0                ! FISH LEAVES SYSTEM
                fishes(fn,16)=3.0  ! Particle leaves at external head BC Upstream
                WRITE(DIAGFN,*) 'At External Upstream Head Boundary: Particle Leaves System on JDAY:',jday, 'Data:',fishes(fn,:)
           GOTO 20
              ELSEIF(FIMP.EQ.CUS(FNBP) .AND. UHS(FNBP)==-1 .AND. FXVEL(5) >= 0.0) THEN 
             FXLOC = DLX(FIMP)*XREFL    ! REflect partilce off UPSTREAM boundary 
              ELSEIF(FIMP.EQ.CUS(FNBP) .AND. UHS(FNBP)>0 .AND. FXVEL(5) >= 0.0) THEN     
                    FXLOC = DLX(FIMP)*XREFL    ! REflect partilce off UPSTREAM boundary   
              ELSEIF(FIMP.EQ.CUS(FNBP) .AND. UHS(FNBP)>0 .AND. FXVEL(5) < 0.0) THEN     
                WRITE(DIAGFN,*) 'UHS: Particle move to new branch on JDAY:',jday,'FXVEL=',fxvel(5),'OLD I:',FIMP,'NEW I:', DHS(FNBP),' OLD branch:',FNBP
                FIMP=UHS(FNBP)
                FYLOC= B(FKMP,FIMP)*0.5   ! PLACE IN MIDDLE LATERALLY
                FXLOC= DLX(FIMP)*0.5      ! PLACE IN MIDDLE OF SEGMENT                           
              ELSEIF(FIMP.EQ.CUS(FNBP) .AND. UHS(FNBP)<-1 ) THEN       ! at a DAM
                    FXLOC = DLX(FIMP)*XREFL    ! REflect partilce off UPSTREAM boundary                                     
              ELSEIF(FIMP>CUS(FNBP))THEN
            FIMP = FIMP - 1                                    ! There are still segments upstream for the fish to
            FXLOC = DLX(FIMP) + FXLOC                          !    move to
              ELSE
            WRITE(DIAGFN,*)'***NO RESOLUTION OF FXLOC<0*** JDAY,I,K,BR#,FN#,FYLOC,FXLOC:',JDAY,FIMP,FKMP,FNBP,FN,FYLOC,FXLOC
          END IF

        !ELSEIF(FXLOC.EQ.0.0.and.cus(fnbp).eq.fimp.and.fxvel(5).eq.0.)then
        !     FXLOC = DLX(FIMP)*XREFL    ! REflect partilce off UPSTREAM boundary 
        !     
        !ELSEIF(FXLOC == DLX(FIMP) .AND. DS(FNBP)==FIMP   .AND.  FXVEL(5) <= 0.0)THEN
        !     FXLOC = DLX(FIMP)*(1.0-XREFL)  ! REFLECT OFF DOWNSTREAM BOUNDARY
        !ELSE

        END IF


!Check for Boundary Violations: Vertical Direction Check (1 of 2)

        KTWBF = KTWB(FJR)                            ! Water Surface Layer of waterbody fish is located in
! This first section below checks and corrects fish vertically above or below the bottom
   70   IF ((FKMP.Le.KTI(FIMP)).OR.&   ! SW 2/01/01 Change so that do a check if in "air" or not
            (FKMP.GT.KBI(FIMP))) THEN                            ! SW 8/2018 USE KBI (INITIAL KB) RATHER THAN ALTERED KB IN CASE OF LOWERING KB IN RIVER SLOT
          IF (FKMP.GT.KBI(FIMP)) THEN                   !  FKMP is below KB 
            FKMP = FKMP - 1
            FZLOC = H(FKMP,FJR)*(1-ZBOTREFL)           !    Logic for rebounding off of the bottom
            GOTO 70
          elseif(fkmp.le.kti(fimp))then               ! Don't leave fish in the air - move to surface layer
          fkmp=kti(fimp)                 ! check and make sure the fish/particles are below the water surface   ! SW 2/01/01
               if(kti(fimp).eq.KTWBF)then
                  if(fzloc.lt.z(fimp))then
                   fzloc=z(fimp)+h(fkmp,fjr)*(1-zsurrefl)
                  end if
               else
                   if(fzloc.le.h(kti(fimp),fjr)+z(fimp))then
                     if(z(fimp).gt.0.0)write(DATADEBUGFN,*)'Debug Error in vertical reflection:FN,I,K,JDAY',fn,fimp,fkmp,jday
                   fzloc=(h(kti(fimp),fjr)+z(fimp))+h(fkmp,fjr)*(1-zsurrefl)
                   end if
               end if
          END IF                                        
        END IF

!Check for Boundary Violations: Vertical Direction Check (2 of 2)

   30   IF (FZLOC >= H(FKMP,FJR)) THEN                     ! Did fish move down below current layer?  -- YES
          IF (DEBUG) WRITE(DATADEBUGFN,*) 'FZLOC GT H'
          IF (FKMP+1.GT.KBI(FIMP)) THEN                 ! Did fish move below bottom layer? -- YES   USE KBI RATHER THAN KB SW 8/2018
                IF(HITSTICKBOTTOM)THEN
                  ! PARTICLE LOST FROM SYSTEM
                                              ! LOSE PARTCLE 
                ISWITCH = 1
                fishes(fn,14)=JDAY               ! Time fish left system
                fishes(fn,15)=JDAY-DELAYDATE(FN)     ! Detention time of fish in system
                fishes(fn,12)=1.0                ! FISH LEAVES SYSTEM
                fishes(fn,16)=2.0  ! Particle hits bottom and sticks
                WRITE(DIAGFN,*) 'Particle Hit and Stick Bottom: Particle Leaves System on JDAY:',jday, 'Data:',fishes(fn,:)
                  ELSE
                    FZLOC = H(FKMP,FJR)*(1-ZBOTREFL)               !    Logic for rebounding off of the bottom
                  ENDIF
                  
          ELSE                                         ! Did fish move below bottom layer? -- NO
            FZLOC = FZLOC - H(FKMP,FJR)
            FKMP = FKMP + 1
            IF (FZLOC.GT.H(FKMP,FJR)) GOTO 30              ! Did fish move down more than one layer 
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
        END IF

!Check for Boundary Violations: Lateral Direction Check (2 of 2)

        IF (FYLOC.LT.0) THEN
          FYLOC = B(FKMP,FIMP)*YREFL 
        ELSE IF (FYLOC.GT.B(FKMP,FIMP)) THEN
          FYLOC = B(FKMP,FIMP)*(1-YREFL) 
        ELSE
        END IF

      END IF

   20 CONTINUE

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
      
      ! CHECK MONITORING STATIONS
      DO N=1,NMONITORS
          IF(MONITORONOFF(FN,N)==0)THEN
              IF(IMONITOR(N)==FIMP)THEN
                  IF(XSHIFTMONITOR(N) >= FXLOC)THEN
                      MONITORONOFF(FN,N)=1
                      TMONITOR(FN,N)=JDAY
                  ENDIF
              ENDIF
          ENDIF
      ENDDO
      

      IF(ISWITCH==0)CALL HISTOGRAM
      
  131 CONTINUE   !end of NFISH loop

   28 CONTINUE

      LRUNDAY = JDAY


!*******************************************************
!  O U T P U T   A N D   E N D   S U B R O U T I N E   *
!*******************************************************

      IF ((FCOUNT.EQ.1).OR.(NIT.EQ.0)) THEN
        CALL FISHPLOT                                                   ! This plots the location of virtual fish
      ENDIF

      return
    END
    
    SUBROUTINE NEXTBRANCH   ! SW 8/7/2016
    USE Fishy; Use GLOBAL; USE STRUCTURES
    INTEGER :: JS, JG, JP
    NEXT_BRANCH=.FALSE.      ! Set all branch connections to FALSE
    DO JB=1,NBR
      IF(DHS(JB) == 0)THEN
          DO JS=1,NSP
              IF(IUSP(JS)==DS(JB).AND.IDSP(JS).NE.0)THEN
                 FNEXT(JB)=IDSP(JS)
                 NEXT_BRANCH(JB)=.TRUE.
                 EXIT
              ENDIF
          ENDDO
          DO JG=1,NGT
               IF(IUGT(JG)==DS(JB).AND.IDGT(JG).NE.0)THEN
                 FNEXT(JB)=IDGT(JG)
                 NEXT_BRANCH(JB)=.TRUE.
                 EXIT
              ENDIF
          ENDDO
          DO JP=1,NPU
                 IF(IUPU(JP)==DS(JB).AND.IDPU(JP).NE.0)THEN
                 FNEXT(JB)=IDPU(JP)
                 NEXT_BRANCH(JB)=.TRUE.
                 EXIT
                ENDIF
          ENDDO
      ENDIF
    ENDDO
    
      
      RETURN
    
    
    
    END SUBROUTINE NEXTBRANCH

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
       ! IF(BR_INACTIVE(JB))CYCLE  ! SW 6/12/2017  No Need to do this for an inactive branch since this is called only once to set up geometery
          IF (I.EQ.DHS(JB)) THEN                    ! If find a connecting branch
            IF (PHI0(DS(JB)).GT.PHI0(I)) THEN
              BRANGLE = PHI0(DS(JB))-PHI0(I)        ! BRANGLE = Angle of Incoming Branch relative to
            ELSE                                    !           Branch I
              BRANGLE = PHI0(DS(JB))+(2.*(22./7.)-PHI0(I))
            END IF
            IF ((BRANGLE.GT.0).AND.(BRANGLE.LT.(22./7.))) THEN  ! From BRANGLE, one can determine if the
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
        IF ((NR.GT.1).OR.(NL.GT.1)) THEN                        ! CHECK: Can only have one Branch join a segment from each side from the same side
          WRITE(DATADEBUGFN,*) 'ERROR: More than 1 branch joining segment'
          WRITE(DATADEBUGFN,9140) I
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
      
    Use Fishy; Use GEOMC; USE GLOBAL;Use SCREENC
    
    IMPLICIT NONE
    integer trib

   81 TRIB = 1                                              ! TRIB = Tributary (or Branch)
   80 IF ((FIMP.GE.US(TRIB)).AND.(FIMP.LE.DS(TRIB))) THEN
        FNBP = TRIB
      ELSE
        TRIB = TRIB + 1
        IF (TRIB.GT.NBR) THEN                                 ! CHECK: TRIB can't be more
          WRITE(DATADEBUGFN,*) 'ERROR: Finding new FNBP in Subroutine'  !        than NBR # BRANCHES
          WRITE(DATADEBUGFN,*) 'ERROR: FIMP=',FIMP,' JDAY=',JDAY
          WRITE(DATADEBUGFN,*) 'FIMP=FIMP-1'
          FIMP=FIMP-1
          GO TO 81
          !STOP
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
        WRITE(DATADEBUGFN,*) 'ERROR: Finding current WB in Subroutine'!        than NWB
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
      Use SCREENC; USE MAIN, ONLY: KBI
      
      IMPLICIT NONE

      INTEGER     WB,NN,NNLAST,WATER
      INTEGER     BRCHNN,KK,II,ND,ELE,EE,K
      INTEGER     TOTNN
      INTEGER     JOININGTRIB
      integer, allocatable, save, DIMENSION (:) :: TOTELE,TOTBRNN
      logical, SAVE  ::  FEGRID
      REAL :: XDIST,ZDIST,VERTFLOW

      CHARACTER*1  :: DTYP1
      CHARACTER*2  :: DTYP2
      CHARACTER*6  :: FILE_PREFIX='Branch'
      CHARACTER*4  :: FILE_SUFFIX='.dat'
      CHARACTER*12 FILENAME

      IF (NIT.EQ.0)then
      FEGRID = .FALSE.                  ! The first pass is used to set up the FE grid
      allocate(totele(NBR),totbrnn(NBR))    ! SW 1/14/01
      end if

   93 IF (FEGRID) CALL INTERCONST                     ! Subroutine to interpolate Constituent values
      IF (FEGRID) CALL INTERFLOWF                     ! Subroutine to interpolate Flow Field values
      NN = 0                                          ! NN = Global Node # (UpperLeft Corner of cell(K,I))
      !DO 115 WB=1,NWB                                 ! WB = Water Body # / NWB = Total # of Water Bodies
      DO WB=1,NWB 
        !DO 116 JB=BS(WB),BE(WB)                     ! JB = Branch #
          DO JB=BS(WB),BE(WB) 
          IF (FEGRID.EQV..FALSE.) GOTO 94             ! Skip to setting up FE grid if NIT = 0

          IF (NIT.EQ.0) THEN
            IF (JB.LT.10) THEN                        ! Creating the output file names based on Branch #
            WRITE(DTYP1,9210) JB
            FILENAME=FILE_PREFIX//DTYP1//FILE_SUFFIX
 9210       FORMAT(I1)
            ELSE IF (JB.LE.64) THEN                   ! TecPlot is limited to 128 Data Sets (i.e. Data Files)
            WRITE(DTYP2,9220) JB                    !   That leaves: 64 Data Sets for Grid Display and
            FILENAME=FILE_PREFIX//DTYP2//FILE_SUFFIX  !                64 Data Sets for Fish Display
 9220       FORMAT(I2)
            ELSE                                        ! CHECK
            WRITE(DATADEBUGFN,*) 'ERROR: Exceeded Maximum Number of TecPlot Data Sets ==> # Data Sets = # of Branches'
            STOP
            END IF          
              
            OPEN(20000+JB,FILE=FILENAME,STATUS='UNKNOWN')       ! Opening the Output Files Based on Branch #
            WRITE(20000+JB,9240) JB                           ! Creating Output File Header for TecPlot
 9240       FORMAT('TITLE = "GRID NODE INFO for Branch',I6,'"')
            WRITE(20000+JB,9270)                                ! Creating Output File Header for TecPlot
 9270       FORMAT('VARIABLES = "X", "Z", "WATER", "HorizVel", "VertVel"&
      , "Temp", "DO", "KLayer", "ISegment"')
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
         ! DO 117 K=2,KMX                                          ! K = Layer #
         DO K=2,KMX
            ZDIST = ZDIST - H(K,WB)                                  ! ZDIST = Depth to Node in Branch JB
            XDIST = 0                                             ! XDIST = X Distance to Node in Branch JB
             DO I=US(JB),DS(JB)+1          ! I = Segment #!  DO I=CUS(JB),DS(JB)+1             ! DO 118  I=US(JB),DS(JB)+1         ! I = Segment #
              !DO K=KTWB(WB),KB(I)   ! SW 5/2015
              IF ((K-1.LE.KBI(I)).OR.(K-1.LE.KBI(I-1))) THEN        ! All Nodes at or above Water Body Bottom  USE KBI INIITIAL KB SW 8/2018
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
                        FLOWFIELD(K,I,1),-VERTFLOW,&                !   to a file for TecPlot to display
                        WQFIELD(K,I,1),WQFIELD(K,I,2),K,I

 9260             FORMAT(F15.1,F10.2,I7,E15.2,E15.2,F8.2,F8.2,I5,I5)
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
        ENDDO
    ENDDO
    
 ! 118       CONTINUE        
  !117     CONTINUE
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
                IF (KK.LE.KBI(II)) THEN                     ! Any Layer above Water Body Bottom   USE KBI RATHER THAN KB SW 8/2018
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
  !116   CONTINUE
  !115 CONTINUE
          ENDDO
      ENDDO
      

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

      Use FISHY; USE GLOBAL; Use GEOMC; Use SCREENC; USE KINETIC; USE MAIN, ONLY: KBI
      
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
      !DO 127 WB=1,NWB                                              ! WB = Water Body # / NWB = Total # of Water Bodies
      DO WB=1,NWB
          DO JB=BS(WB),BE(WB)
              DO K=2,KMX
                  DO I=US(JB),DS(JB)+1
        !DO 128 JB=BS(WB),BE(WB)                                    ! JB = Branch #
        !  DO 129 K=2,KMX                                           ! K = Layer #
        !    DO 130 I=US(JB),(DS(JB)+1)                         ! I = Segment #  changed from us to cus sw 5/2015
              IF ((K-1.LE.KBI(I)).OR.(K-1.LE.KBI(I-1))) THEN         ! All Nodes at or above Water Body Bottom,  KBI IS INITIAL KB  !((K-1.LE.KB(I)).OR.(K-1.LE.KB(I-1))) THEN 
                NN = NN + 1                                        ! PRE-CHECK CALCULATION
                KK = NDINFO(NN,3)                                  ! PRE-CHECK CALCULATION
                II = NDINFO(NN,4)                                  ! PRE-CHECK CALCULATION
                IF ((I.NE.II).OR.(K.NE.KK)) THEN                   ! CHECK
                  
                  WRITE(DATADEBUGFN,*) 'ERROR: Nodes not matching up for Constituent Interpolation'
                  write(DATADEBUGFN,*)'JDAY=',jday
                  write(DATADEBUGFN,*)'Waterbody=',wb,' Branch JB=', JB
                  write(DATADEBUGFN,*)'NN=',nn,' kk=',kk,' ii=',ii,' k=',k,' i=',i
                  !STOP
                !ELSE
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
                ELSE IF ((K.EQ.KBI(I)+1).AND.(I.EQ.US(JB))) THEN       !Analogous to Node 7 above                   K.EQ.KB(I)+1).AND.(I.EQ.US(JB))) 
                  WQFIELD(K,I,1) = T2(K-1,I)
                  WQFIELD(K,I,2) = O2(K-1,I)
                ELSE IF ((K.EQ.KBI(I)+1).AND.(I.EQ.DS(JB)+1)) THEN     !Analogous to Node 9 above                    K.EQ.KB(I)+1).AND.(I.EQ.DS(JB)+1))
                  WQFIELD(K,I,1) = T2(K-1,I-1)
                  WQFIELD(K,I,2) = O2(K-1,I-1)
                ELSE IF (K.EQ.KBI(I)+1) THEN                             !Analogous to Node 8 above                (K.EQ.KB(I)+1)
                  IF(KB(I) /= KBI(I))THEN
                    WQFIELD(K,I,1) = T2(KTWB(WB),I-1)   ! SW 8/20/2018 IN CASE OF LOWERING OF RIVER WATER LEVEL BECAUSE OF A "NOTCH"
                    WQFIELD(K,I,2) = O2(KTWB(WB),I-1)
                  ELSE
                  XRANGE = .5*DLX(I-1)+.5*DLX(I)
                  WQFIELD(K,I,1) = (T2(K-1,I-1)*(XRANGE-.5*DLX(I-1))+&
                                   T2(K-1,I)*(XRANGE-.5*DLX(I)))/XRANGE
                  WQFIELD(K,I,2) = (O2(K-1,I-1)*(XRANGE-.5*DLX(I-1))+&
                                   O2(K-1,I)*(XRANGE-.5*DLX(I)))/XRANGE
                  ENDIF
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


              !ELSE
              END IF
  !130       CONTINUE
  !129     CONTINUE
  !128   CONTINUE
  !127 CONTINUE
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    

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

      Use FISHY; USE GLOBAL; Use GEOMC; Use SCREENC; USE MAIN, ONLY: KBI
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
      !DO 123 WB=1,NWB                                              ! WB = Water Body # / NWB = Total # of Water Bodies
      DO WB=1,NWB
        !DO 124 JB=BS(WB),BE(WB)                                  ! JB = Branch #
          DO JB=BS(WB),BE(WB)
          !DO 125 K=2,KMX                                           ! K = Layer #
              DO K=2,KMX
            !DO 126 I=US(JB),(DS(JB)+1)                         ! I = Segment #  SW 5/2015
               DO I=US(JB),DS(JB)+1
              IF ((K-1.LE.KBI(I)).OR.(K-1.LE.KBI(I-1))) THEN         ! All Nodes at or above Water Body Bottom     ((K-1.LE.KB(I)).OR.(K-1.LE.KB(I-1))) 
                NN = NN + 1                                        ! PRE-CHECK CALCULATION
                KK = NDINFO(NN,3)                                  ! PRE-CHECK CALCULATION
                II = NDINFO(NN,4)                                  ! PRE-CHECK CALCULATION
                IF ((I.NE.II).OR.(K.NE.KK)) THEN                   ! CHECK
                  WRITE(DATADEBUGFN,*) 'ERROR: Nodes not matching up for Flowfield Interpolation'
                  write(DATADEBUGFN,*)'WB=',wb,' JB=', JB
                  write(DATADEBUGFN,*)'NN=',nn,' kk=',kk,' ii=',ii,' k=',k,' i=',i

                  !STOP
                !ELSE
                END IF

!LINEAR INTERPOLATION of the Flow Field (i.e. Velocities)
                IF (LINEAR) THEN              ! LINEAR INTERPOLATION
                  IF (K.LT.KTWB(WB)) THEN                           !Above Water Surface
                    IF(K < KTI(I))THEN   ! SW 12/12/2018
                    FLOWFIELD(K,I,1) = U(KTWB(WB),I-1)         ! SW 12/26/2018 was =0 before      !   Horiz. Vel.   IN AIR should still be in the upper layer - nothing goes up into the air - corrected soon by BC
                    FLOWFIELD(K,I,2) = 0.0                         !   Vert.  Vel.
                    ELSE 
                    FLOWFIELD(K,I,1) = U(KTWB(WB),I-1)             !   Horiz. Vel.  NOT IN AIR  SW 12/12/2018
                    FLOWFIELD(K,I,2) = 0.0                         !   Vert.  Vel.  
                    ENDIF
                  ELSE IF (K.EQ.KTWB(WB)) THEN                      !At Water Surface
                    FLOWFIELD(K,I,1) = U(K,I-1)                    !   Horiz. Vel.
                    FLOWFIELD(K,I,2) = 0.0                         !   Vert.  Vel.
                  ELSE IF (K.EQ.KBI(I)+1) THEN                      !Water Body Bottom    (K.EQ.KB(I)+1)
                    IF(KB(I) /= KBI(I))THEN
                        FLOWFIELD(K,I,1) = U(KTWB(WB),I-1)    ! IN CASE WATER LEVELS ARE LOWERED IN A RIVER BELOW KB SW 8/20/2018
                        FLOWFIELD(K,I,2) = 0.0
                    ELSE
                        FLOWFIELD(K,I,1) = U(K-1,I-1)
                        FLOWFIELD(K,I,2) = 0.0
                    ENDIF
                    
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
                    HVERT(1) = -H(K,WB)/10000.          ! Arbitrary value
                    HVERT(2) = -H(K,WB)/10001.          ! Arbitrary value
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
                  IF (K.GT.KB(I)) THEN              !Both Bottom Nodes are below Water Body Bottom   12/12/2018 Change back to KB USE KBI SW 8/2018
                    UVERT(3) = 0.0                    ! Helps satisfy no-slip condition on bottom (i.e., U=0)
                    UVERT(4) = 0.0                    ! Helps satisfy no-slip condition on bottom (i.e., U=0)
                    HVERT(3) = H(K-1,WB)/10001.         ! Arbitrary value
                    HVERT(4) = H(K-1,WB)/10000.         ! Arbitrary value
                  ELSE IF (K+1.GT.KB(I)) THEN       !Only Bottom-most Node below Water Body Bottom   12/12/2018 Change back to KB  USE KBI SW 8/2018
                    UVERT(3) = U(K,I-1)
                    UVERT(4) = 0.0                    ! Helps satisfy no-slip condition on bottom (i.e., U=0)
                    HVERT(3) = 0.5*H(K,WB)
                    HVERT(4) = H(K,WB)                 ! Helps satisfy no-slip condition on bottom (i.e., U=0)
                  ELSE                              !Both Bottom Nodes exist
                    UVERT(3) = U(K,I-1)
                    UVERT(4) = U(K+1,I-1)
                    HVERT(3) = 0.5*H(K,WB)
                    HVERT(4) = H(K,WB)+0.5*H(K+1,WB)
                  END IF
                  IF (I-1.LT.US(JB)) THEN         !Both Left Nodes are Upstream of Water Body
                    WHORZ(1) = 0.0                    ! Helps satisfy no-slip condition at Upstream Boundary
                    WHORZ(2) = 0.0                    ! Helps satisfy no-slip condition at Upstream Boundary
                    DLXHO(1) = -DLX(I)/10000.        ! Arbitrary value
                    DLXHO(2) = -DLX(I)/10001.        ! Arbitrary value
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


! Begin Calculation of Accelerations ==> Time Derivative of Velocity (Backward Differencing)
                IF (NIT.NE.0) THEN
                  IF ((JDAY-LASTJDAY).EQ.0) THEN
                    FLOWFIELD(K,I,3) = 0.0
                    FLOWFIELD(K,I,4) = 0.0
                  ELSE
                    FLOWFIELD(K,I,3) = (FLOWFIELD(K,I,1)-LASTFIELD(K,I,1))/((JDAY-LASTJDAY)*86400.)    ! 24.*3600.
                    FLOWFIELD(K,I,4) = (FLOWFIELD(K,I,2)-LASTFIELD(K,I,2))/((JDAY-LASTJDAY)*86400.)
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

              !ELSE
              END IF
  !126       CONTINUE
  !125     CONTINUE
  !124   CONTINUE
  !123 CONTINUE
               ENDDO
              ENDDO
          ENDDO
      ENDDO
      

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
      CHARACTER*4 :: FILE_PREFIX='Part'
      CHARACTER*4 :: FILE_SUFFIX='.dat'
      CHARACTER*10 FILENAME

      IF (NIT.EQ.0) NZONES = 0                        ! This # is inputted into a TecPlot Macro for animation
      NZONES = NZONES + 1                             ! This # is inputted into a TecPlot Macro for animation
      DO 132 WB=1,NWB                                 ! WB = Water Body # / NWB = Total # of Water Bodies
        DO 133 JB=BS(WB),BE(WB)                     ! JB = Branch #
          NBRF(JB) = 0                              ! NBRF(JB) = # of Fish in Branch JB
          DO 134 FN=1,NFISH                           ! NFISH = Total # of Fish in System
            IF (INT(FISHES(FN,12)).EQ.1 .or. DELAYDATE(FN)> JDAY) GOTO 134     ! If Fish left system do not output its info to TecPlot
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
            END IF
  134     CONTINUE

          IF (NIT.EQ.0) THEN
              
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
              
            OPEN(30000+JB,FILE=FILENAME,STATUS='UNKNOWN')             ! Opening the Output Files Based on Branch #
            WRITE(30000+JB,9300) JB                                 ! Creating Output File Header for TecPlot
 9300       FORMAT('TITLE = "PARTICLE INFO for Branch',I6,'"')
            WRITE(30000+JB,9310)                                      ! Creating Output File Header for TecPlot
 9310       FORMAT('VARIABLES = "X", "Z", "FSIZE", "FAGE", "WATER", "ISEG", "K"')
          ELSE
            WRITE(30000+JB,*) ' '
          END IF
          WRITE(30000+JB,9320) JDAY,(1+NBRF(JB))                    ! Creating Output File Header for TecPlot
 9320     FORMAT('ZONE T="JDAY ',F9.2,'", I=',I7,', F=POINT')
          WRITE(30000+JB,9330) 0.01,0.01,0.0,0.0,-1,0,0         ! This is a dumby fish so TecPlot can at least
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
            WRITE(30000+JB,9330) XLOC,ZLOC,FSZ,FAG,WATER,I,K
 9330       FORMAT(F15.2,F10.2,F9.2,F9.2,I4,1X,I3,1X,I3)
  135     CONTINUE
          WRITE(30000+JB,9440) MONTH,GDAY,YEAR,INT(HR),AMPM,&
                                 JDAY,NZONES                  ! Output info used as dynamic text in TecPlot Animation
 9440     FORMAT('TEXT X=47, Y=90, F=HELV-BOLD, HU=FRAME, AN=MIDCENTER,&
       C=RED, H=2.5, T="',A10,I3,',',I5,4X,I3,A2,'    (JDAY',F8.3,')", ZN=',I6)
          WRITE(30000+JB,9518) NBRF(JB),NZONES            ! Output info used as dynamic text in TecPlot Animation
9518      FORMAT('TEXT X=30.0, Y=21.0, F=HELV-BOLD, HU=FRAME,&
      AN=MIDRIGHT, C=BLACK, H=2.1, T="',I6,'", ZN=',I6)

          !  The Following is Used for Creating and Positioning Static Text in the TecPlot Animation

          IF (NIT.EQ.0) THEN
            WRITE(30000+JB,9541) 
 9541       FORMAT('TEXT X=25.0, Y=25.0, F=HELV-BOLD, HU=FRAME,&
      AN=MIDCENTER, C=BLACK, H=2.1, T="# of Particles"')
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
      Use GEOMC; USE MAIN, ONLY: KBI
      
      IMPLICIT NONE
     
      REAL      XH1,XH2,H3,H4,DV,DW,DX,DY

 !Check for Boundary Violations: Horizontal Plane

   95 IF (XLOK.GT.DLX(ILOK)) THEN                             ! Must look into the next segment downstream
        IF ((ILOK+1).GT.DS(FNBP)) THEN                        ! Already at most downstream segment in branch
          XLOK = DLX(ILOK)
        ELSE IF (KLOK.GT.KBI(ILOK+1)) THEN                    ! SW USE KBI 8/2018
          XLOK = DLX(ILOK)
        ELSE                                                  ! More segments exist downstream of fish location
          XLOK = XLOK - DLX(ILOK)
          ILOK = ILOK + 1                                     ! Looking into next segment downstream
        END IF
        IF (XLOK.GT.DLX(ILOK)) GOTO 95
      ELSE IF (XLOK.LT.0) THEN                                ! Must look into the next segment upstream
        IF (ILOK.EQ.CUS(FNBP)) THEN                           ! Already at most upstream segment in branch
          XLOK = 0.0
        ELSE IF (KLOK.GT.KBI(ILOK-1)) THEN                    ! SW 8/2018
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
        IF (KLOK.EQ.KBI(ILOK)) THEN                            ! Already at water body bottom    SW 8/2018
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
INTEGER :: CATCH,DINT,DINTLAST,CATCHDEPTH,DOTALLY,TEMPTALLY, FIJ, JF, TAGGED, FVAR2,FVAR3, TAGGEDDEPTH, II


! final fish output
   IF(NMONITORS<10)THEN
       write(FINALFN,'(A214,<NMONITORS>(A11,I1,","))')'Part#,Seg#,XLocationwithinSegmentfromUpstreamSide(m),Layer#,VerticalDistfromTop(m),LateralDistfromLeftBank,Branch#,ParticleInModel(=0),JDAYleftsystem,DetentionTime(days),RemovalMechanism,SedVelocity(m/d),DateStart,',('MonitorDate',I,I=1,NMONITORS) !Monitor2Date, Monitor3Date'
   ELSE
       write(FINALFN,'(A214,<NMONITORS>(A11,I2,","))')'Part#,Seg#,XLocationwithinSegmentfromUpstreamSide(m),Layer#,VerticalDistfromTop(m),LateralDistfromLeftBank,Branch#,ParticleInModel(=0),JDAYleftsystem,DetentionTime(days),RemovalMechanism,SedVelocity(m/d),DateStart,',('MonitorDate',I,I=1,NMONITORS) !Monitor2Date, Monitor3Date'
   ENDIF
   
    do jf=1,nfish     
    write(FINALFN,'(i7,",",1x,10(f10.3,","),f8.4,",",f12.4,",",100(f10.3,","))')jf,fishes(jf,1),fishes(jf,2),fishes(jf,3),fishes(jf,4),fishes(jf,5),fishes(jf,6),fishes(jf,12),fishes(jf,14),fishes(jf,15),fishes(jf,16),sedvel(jf),delaydate(jf),(tmonitor(jf,ii),ii=1,nmonitors)
    end do

    
    CALL HISTOGRAM_OUTPUT
    close(FINALFN)
! End fish output section
                                                             !FISH

      Return                                                               !FISH
End Subroutine FishOutput

! ************************ READ FISH INPUT DATA FILES

Subroutine Read_Fish_Data

Use Fishy; Use SCREENC, ONLY:JDAY; USE MAIN, ONLY: FISH_PARTICLE_EXIST 
IMPLICIT NONE

REAL, ALLOCATABLE, DIMENSION(:) :: FXLOCINIT,FZLOCINIT,SEDVELINIT,DELAYDATEINIT
character*3 ALINE
INTEGER :: NG,I,J,NA, IDEBUG, NFISHLAST,ILINEAR,HTSTBOT,HTSTSIDE

            HIST_V=.FALSE.
            HIST_T=.FALSE.
            HIST_D=.FALSE.

    open(DIAGFN,file='w2_particle.csv',status='old')
    READ(DIAGFN,*)
    READ(DIAGFN,*)
    READ(DIAGFN,*)PARTON,NFISHSEG,NFISHPCEL,ALINE,DXTHEORY,OUTFREQP,ILINEAR,HTSTBOT,HTSTSIDE,IDEBUG   !,ALPHAX,ALPHAZ     !'(/10x,7x,a3,i10,i10,f10.0,(7x,a3),1f10.0,7X,A3)
    If(PARTON.eq.'ON')then
    PARTICLE=.true.
    else
    PARTICLE=.false.
    FISH_PARTICLE_EXIST=.FALSE.  ! STOPS ALL FURTHER 
    RETURN
    end if
    If(ALINE.eq.'ON')then
    LINE=.true.
    else
    LINE=.false.
    end if
    if(outfreqp==0.0)outfreqp=1.0   ! error trapping

! LINE: If LINE, distribute NFISHCEL not at a point but linearly along the line defined by IFISH and IFISHT and IFISHB; IF LINE, FZLOC=0
! PARTICLE = ON only "dumb" particle transport with random fluid motion + Sedimentation Vel (SEDVEL)
! NFISHPCEL: # of particles added per cel
! NFISHSEG: # of segements to add particles to
! DATE:  Date particles are deposited
    if(nfishseg.eq.0)nfishseg=1
    Allocate(ifish(nfishseg),ifisht(nfishseg),ifishb(nfishseg),delaydateinit(nfishseg),fxlocinit(nfishseg),fzlocinit(nfishseg),sedvelinit(nfishseg)) 
        READ(DIAGFN,*)
        DO I=1,NFISHSEG
        READ(DIAGFN,*)NA,IFISH(I),IFISHT(I),IFISHB(I),FXLOCINIT(I),FZLOCINIT(I),SEDVELINIT(I),DELAYDATEINIT(I)
        
        IF(DELAYDATEINIT(I) < JDAY)DELAYDATEINIT(I)=JDAY
        ENDDO
   
    !READ(DIAGFN,'(//(10X,9I10))')(IFISH(I),I=1,NFISHSEG)
    !READ (DIAGFN,'(//(10X,9I10))')(IFISHT(I),I=1,NFISHSEG)
    !READ (DIAGFN,'(//(10X,9I10))')(IFISHB(I),I=1,NFISHSEG)
    ! IFISH: SEG#s where fish are added
    ! IFISHT:Uppermost top Layer where Fish are added in Segment IFISH
    ! IFISHB: Bottommost Layer where Fish are added in Segment IFISH
  NFISH=0
  do i=1,nfishseg
    NFISH=NFISH+(ifishb(i)-ifisht(i)+1)*nfishpcel
  end do
  
  ALLOCATE (GROUP(NFISH),FXLOCI(NFISH),FZLOCI(NFISH),DELAYDATE(NFISH),SEDVEL(NFISH))
  
    NFISHLAST=1
    do i=1,nfishseg
        DO J=NFISHLAST,(NFISHLAST-1)+(ifishb(i)-ifisht(i)+1)*nfishpcel
        GROUP(J)=I
        FZLOCI(J)=FZLOCINIT(I)
        FXLOCI(J)=FXLOCINIT(I)
        DELAYDATE(J)=DELAYDATEINIT(I)
        SEDVEL(J)=SEDVELINIT(I)/86400.    ! convert m/d to m/s
        ENDDO
    NFISHLAST=J
    end do
   
        READ(DIAGFN,*,END=200)  ! TEST CODE IN CASE SOMEONE DOES NOT HAVE THESE PARTS OF THE CONTROL FILE
        READ(DIAGFN,*)NUMCLASS
        IF(NUMCLASS == 0)GO TO 200
            ALLOCATE(SUMVOLT(NFISH))
            sumvolt=0.0
        READ(DIAGFN,*)
        READ(DIAGFN,*)VEL_INT,VEL_TOP,TEMP_INT,TEMP_TOP,D_INT,D_TOP
        IF(VEL_INT /= 0.0)THEN
            HIST_V=.TRUE.
            allocate (v_class(NFISH,numclass))
            ALLOCATE(V_TOT(NFISH),V_CNT(NFISH),V_SUM(NFISH),V_AVG(NFISH))   
            v_cnt=0.0
            v_class=0.0
            v_tot=0.0
            V_SUM=0.0
            V_AVG=0.0
        ENDIF
        
        IF(TEMP_INT /= 0.0)THEN
            HIST_T=.TRUE.
             ALLOCATE (T_TOT(NFISH),T_CNT(NFISH),T_SUM(NFISH),T_AVG(NFISH))
             allocate (t_class(NFISH,numclass))
            t_cnt=0.0
            t_class=0.0
            t_tot=0.0
            T_AVG=0.0
            T_SUM=0.0
        ENDIF
        
        IF(D_INT /= 0.0)THEN
            HIST_D=.TRUE.  
             ALLOCATE (D_TOT(NFISH),D_CNT(NFISH),D_SUM(NFISH),D_AVG(NFISH))
            allocate(d_class(NFISH,numclass))
           d_cnt=0.0
           d_class=0.0
           d_tot=0.0
           D_AVG=0.0
           D_SUM=0.0
        ENDIF   
        
        READ(DIAGFN,*)              ! SW 11/1/2018
        READ(DIAGFN,*)NMONITORS
        READ(DIAGFN,*)              ! SW 11/1/2018
        ALLOCATE(IMONITOR(NMONITORS),TMONITOR(NFISH,NMONITORS),XSHIFTMONITOR(NMONITORS),MONITORONOFF(NFISH,NMONITORS))
        MONITORONOFF=0   ! IF SET TO ZERO IT IS ACTIVE - WHEN PARTICLE SHOWS UP IT GETS SWITCHED TO 1
        TMONITOR=0.0
        READ(DIAGFN,*)(IMONITOR(I),I=1,NMONITORS)              ! SW 11/1/2018   
        
        READ(DIAGFN,*)  
        READ(DIAGFN,*)(XSHIFTMONITOR(I),I=1,NMONITORS)              ! SW 11/1/2018 
!        READ(DIAGFN,'(//3F10.0,I10)')FXLOC,FZLOC,OUTFREQP,IDEBUG

!          ! FXLOC   = Location of fish within segment IMP from upstream side
!          ! FZLOC   = Location of fish within layer KMP from top side
!          ! FNBP    = Branch where fish is released

200     XREFL    = 0.1                 ! Reflect this percentage of segment length
                                       !   when fish encounters horizontal obstacle
        YREFL    = 0.1                 ! Reflect this % of cell width when encounter stream bank
        ZBOTREFL = 0.5  !SW 12/26/2018 2.0                 ! Reflect this percentage of bottom layer height
                                       !   when fish encounters bottom
        ZSURREFL = 0.5                 ! Reflect this percentage of surface layer height
                                       !   when fish encounters water surface
        
        
        DEBUG=.FALSE.
        LINEAR=.FALSE.
        HITSTICKBOTTOM=.FALSE.
        HITSTICKSIDE=.FALSE.
        IF(HTSTBOT==1)HITSTICKBOTTOM=.TRUE.
        IF(HTSTSIDE==1)HITSTICKSIDE=.TRUE.
        IF(ILINEAR==1)LINEAR=.TRUE.
        IF(IDEBUG==1)DEBUG=.TRUE.
  
    close(DIAGFN)
    RETURN
    
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
    
    REAL :: DXMIN=0.0005, DZMIN1=0.0001,DZMAX1=0.5, COSTHETA, SINTHETA, DX1, DX2, DXAVG, SK, RZ, RX, DZ1, DZ2, DZAVG, r1, r2
    REAL :: DISPX, DISPZ, VEL,WPART,XAREA
    ! DXMIN,DZMIN1,DZMAX1 ! SW 2/01/01   in units of m2/s

        FXLOC = FXLOC + (FXVEL(5)) * (nfsfreq)*24.*3600.    ! New updated Part X-Location
        FYLOC = FYLOC + (FYVEL) * (nfsfreq)*24.*3600.    ! New updated Part Y-Location
        FZLOC = FZLOC + (FZVEL(5) + SEDVEL(FN)) * (nfsfreq)*24.*3600.    ! New updated Part Z-Location

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
! variations in RZ - constrain to DZ=0.5
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

        IF(DXTHEORY==' ON')THEN
            DISPX=5.84e-4*DLX(FIMP)**1.1
            DISPZ=DZAVG
        ELSE
            Dispx=DXAVG  !*alphax            ! Note these should be interpolated rather than using nearest cell #
            Dispz=DZAVG  !*alphaz
        ENDIF

! COMPUTE DZ and DX by interpolating

        call random(seed,r1)
        call random(seed,r2)

!        r1=gasdev(seed)
!        r2=gasdev(seed)

        r1=(r1-0.5)*2.
        r2=(r2-0.5)*2.

        RX=SQRT(6.0*Dispx*nfsfreq*86400.)*(r1*costheta-r2*sintheta)     ! DX is m2/s    nfsfreq is days
        RZ=SQRT(6.0*Dispz*nfsfreq*86400.)*(r1*sintheta+r2*costheta)
        
        !write(9500,*)(r1*costheta-r2*sintheta),(r1*sintheta+r2*costheta)

! constrain random component to segment length and cell layer height
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
! SW 2/01/02   7/31/2017

   Use FISHY
   USE GLOBAL
   Use GEOMC
   Use MAIN, only: EV, QPR, IWD
   
   IMPLICIT NONE
   REAL :: XAREA
   LOGICAL :: WITH_AT_PARTICLE_LOCATION
   INTEGER :: IW
   
! Concept all QSS from each cell will be treated as a lateral withdrawal - each withdrawal will
! be assigned a RHS or LHS looking downstream location; lateral velocity origin is the
! segment/cell center. Velocities to the RHS are + and those to the LHS are -
!  IF RHS LOOKING DOWNSTREAM THEN VELOCITY IS POSITIVE, RIMPBR(I,1)>0
!  IF LHS LOOKING DOWNSTREAM THEN VELOCITY IS NEGATIVE, LIMPBR(I,1)>0   
   ! TRIBS AND INFLOWS IN QSS ARE ASSIGNED POSITIVE VALUES, WITHDRAWALS ARE NEGATIVE
   ! 
    WITH_AT_PARTICLE_LOCATION=.FALSE.
! CHECK FOR WITHDRAWALS AT FIMP
      DO IW=1,NWD
          IF(FIMP==IWD(IW))THEN
              WITH_AT_PARTICLE_LOCATION=.TRUE.
              EXIT
          ENDIF     
      ENDDO

      IF(RIMPBR(FIMP,1) == 0 .AND. LIMPBR(FIMP,1) == 0 .AND. .NOT.WITH_AT_PARTICLE_LOCATION)THEN
        FYVEL=0.0
      ELSE
          
        if(fkmp.le.KTWB(fjr))then   
            xarea=h1(KTWB(FJR),fimp)*dlx(fimp)
            fyvel=(qss(KTWB(fjr),fimp)+ev(fimp)-QPR(FIMP))/xarea        ! ADD EVAPORATION AND REMOVE PRECIP BACK TO QSS...be careful about evaporation which is also a -QSS flow

        IF(RIMPBR(FIMP,1) > 0)THEN
            IF(BR_INACTIVE(RIMPBR(FIMP,2)))THEN
                FYVEL=0.0
            ELSE
                fyvel=-fyvel  
            ENDIF     
        ELSEIF(LIMPBR(FIMP,1) > 0)THEN
            IF(BR_INACTIVE(LIMPBR(FIMP,2)))THEN
            FYVEL=0.0
            ELSE
            fyvel=fyvel   
            ENDIF
        ENDIF

        else
        xarea=h1(fkmp,fIMP)*dlx(fimp)
        fyvel=qss(fkmp,fimp)/xarea
        
        IF(RIMPBR(FIMP,1) > 0)THEN
            IF(.NOT.BR_INACTIVE(RIMPBR(FIMP,2)) .AND. FKMP <= KB(DS(RIMPBR(FIMP,2))))THEN
            fyvel=-fyvel  
            ELSE
            FYVEL=0.0
            ENDIF        
        ELSEIF(LIMPBR(FIMP,1) > 0)THEN
            IF(.NOT.BR_INACTIVE(LIMPBR(FIMP,2)) .AND. FKMP <= KB(DS(LIMPBR(FIMP,2))))THEN
            fyvel=fyvel   
            ELSE
            FYVEL=0.0
            ENDIF
        ENDIF
        
        end if
    ENDIF
    
    return
    end
!************************************   
   SUBROUTINE HISTOGRAM
   USE Fishy; USE GLOBAL; USE GEOMC; USE MAIN, ONLY: KBI
   REAL :: DepthParticle
    
    N=FN   ! N IS THE FISH NUMBER - THIS ROUTINE IS CALLED FOR EACH FISH FN  
   ! Velocity
    K=FKMP
    I=FIMP
    DepthParticle=ELWS(I)-EL(K,I)+FZLOC
    
    IF(HIST_V.OR.HIST_T.OR.HIST_D)sumvolt(N)=sumvolt(N)+dlt   ! ONLY COMPUTED ONCE
     
        IF(HIST_V)THEN
            v_tot(N)=v_tot(N)+u(k,i)*dlt
            v_cnt(N)=v_cnt(N)+dlt
            v_crit=vel_top
            if(u(k,i).ge.vel_top)v_class(N,1)=v_class(N,1)+dlt
            do jj=2,numclass
              if(u(k,i).lt.v_crit.and.u(k,i).ge.v_crit-vel_int)then
                v_class(N,jj)=v_class(N,jj)+dlt
                go to 210
              else
                v_crit=v_crit-vel_int
              end if
              if(jj.eq.numclass.and.u(k,i).lt.v_crit+vel_int)then
                v_class(N,jj)=v_class(N,jj)+dlt
              end if
            end do
        ENDIF
        
210       continue

             if(HIST_T)then

            t_tot(N)=t_tot(N)+t2(k,i)*dlt
            t_cnt(N)=t_cnt(N)+dlt
            t_crit=temp_top
            if(t2(k,i).ge.temp_top)t_class(N,1)=t_class(N,1)+dlt
            do jj=2,numclass
              if(t2(k,i).lt.t_crit.and.t2(k,i).ge.t_crit-temp_int)then
                t_class(N,jj)=t_class(N,jj)+dlt
                go to 200
              else
                t_crit=t_crit-temp_int
              end if
              if(jj.eq.numclass.and.t2(k,i).lt.t_crit+temp_int)then
                t_class(N,jj)=t_class(N,jj)+dlt
              end if
            end do
          end if
200       continue    
       
! Depth
           IF(HIST_D)THEN
            d_tot(N)=d_tot(N)+DepthParticle*dlt           ! SW changed all KB(I to KBI(I    8/2018              d_tot(N)=d_tot(N)+depthb(kbi(i),i)*dlt
            d_cnt(N)=d_cnt(N)+dlt
            d_crit=d_top
            if(DepthParticle.ge.d_top)d_class(N,1)=d_class(N,1)+dlt                   ! if(depthb(kbi(i),i).ge.d_top)d_class(N,1)=d_class(N,1)+dlt
            do jj=2,numclass
              if(DepthParticle.lt.d_crit.and.DepthParticle.ge.d_crit-d_int)then    ! if(depthb(kbi(i),i).lt.d_crit.and.depthb(kbi(i),i).ge.d_crit-d_int)then
                d_class(N,jj)=d_class(N,jj)+dlt
                go to 300
              else
                d_crit=d_crit-d_int
              end if
              if(jj.eq.numclass.and.DepthParticle.lt.d_crit+d_int)then                !if(jj.eq.numclass.and.depthb(kbi(i),i).lt.d_crit+d_int)then
                d_class(N,jj)=d_class(N,jj)+dlt
              end if
            end do
          end if
300       continue

          RETURN
          
   
    END SUBROUTINE HISTOGRAM
    SUBROUTINE HISTOGRAM_OUTPUT
        USE Fishy
        if(HIST_T)then
          DO N=1,NFISH
          if(t_cnt(N).gt.0.0)then
          t_avg(N)=t_tot(N)/t_cnt(N)
          else
          t_avg(N)=0.0
          end if
          ENDDO

        open(CONE,file='envrprf_t_particle.csv',status='unknown')
        write(CONE,*)'Temperature interval,Fraction of time,Number of Particles:',nfish,','
        write(CONE,'("Interval,",<NFISH>("Particle",i4,","))')(I,I=1,NFISH)   !"Interval,",<NFISH>("Particle",i4,",")
        temp_c=temp_top
          do i=1,numclass
          write(CONE,125)temp_c,(t_class(N,i)/sumvolt(N),N=1,NFISH)
          temp_c=temp_c-temp_int
          DO N=1,NFISH
          t_sum(N)=t_sum(N)+t_class(N,i)/sumvolt(N)
          ENDDO
          end do
        write(CONE,'(1x)')
        write(CONE,'(" Sum of fractions, ",<nfish>(e12.4,","))')(t_sum(N),N=1,NFISH)
        write(CONE,'(1x)')
        write(CONE,'(" Average, ",<nfish>(e12.4,","))')(t_avg(N), N=1,NFISH)
        close(CONE)
        end if

        if(HIST_V)then
            DO N=1,NFISH
          if(v_cnt(N).gt.0.0)then
          v_avg(N)=v_tot(N)/v_cnt(N)
          else
          v_avg(N)=0.0
          end if
          ENDDO
        open(CONE,file='envrprf_v_particle.csv',status='unknown')
        write(CONE,*)'Velocity interval,Fraction of volume,Number of Particles:',nfish,','
        write(CONE,'("Interval,",<NFISH>("Particle",i4,","))')(I,I=1,NFISH)
        vel_c=vel_top
          do i=1,numclass
          write(CONE,125)vel_c,(v_class(N,i)/sumvolt(N),N=1,NFISH)
          vel_c=vel_c-vel_int
          DO N=1,NFISH
          v_sum(N)=v_sum(N)+v_class(N,i)/sumvolt(N)
          ENDDO
          end do
        write(CONE,'(1x)')
        write(CONE,'(" Sum of fractions, ",<nfish>(e12.4,","))')(v_sum(N),N=1,NFISH)
        write(CONE,'(1x)')
        write(CONE,'(" Average, ",<nfish>(e12.4,","))')(v_avg(N),N=1,NFISH)
        close(CONE)
        end if
        
        
        if(HIST_D)then
          DO N=1,NFISH
          if(d_cnt(N).gt.0.0)then
          d_avg(N)=d_tot(N)/d_cnt(N)
          else
          d_avg(N)=0.0
          end if
          END DO
        open(CONE,file='envrprf_depth_particle.csv',status='unknown')
        write(CONE,*)'Depth interval,Fraction of time,Number of Particles:',nfish,','
        write(CONE,'("Interval,",<NFISH>("Particle",i4,","))')(I,I=1,NFISH)
        d_c=d_top
          do i=1,numclass
          write(CONE,125)d_c,(d_class(N,i)/d_cnt(N),n=1,nfish)
          d_c=d_c-d_int
          DO N=1,NFISH
          d_sum(N)=d_sum(N)+d_class(N,i)/d_cnt(N)
          ENDDO
          end do
        write(CONE,'(1x)')
        write(CONE,'(" Sum of fractions, ",<nfish>(e12.4,","))')(d_sum(N), N=1,NFISH)
        write(CONE,'(1x)')
        write(CONE,'(" Average, ",<nfish>(e12.4,","))')(d_avg(N), N=1,NFISH)
        close(CONE)
        end if
        

       
125       format((f8.2,',',<NFISH>(e12.4,',')))

    RETURN
    END SUBROUTINE HISTOGRAM_OUTPUT
   