! CE-QUAL-W2 computations
INTEGER(4) FUNCTION CE_QUAL_W2 (DLG)

! IVF/CVF specific code
  USE DFLOGM; USE MSCLIB; USE DFWIN, RENAMED => DLT;

 !DEC$ATTRIBUTES STDCALL   :: ce_qual_w2
 !DEC$ATTRIBUTES REFERENCE :: Dlg
  USE IFPORT                 ! to get current working directory
  USE MAIN
  USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  USE MACROPHYTEC; USE POROSITYC; USE ZOOPLANKTONC  
  USE modSYSTDG, ONLY: INPUT_SYSTDG                              ! systdg  
  Use CEMAVars
  USE CEMASedimentDiagenesis, only: C2SF,InitCond_SedFlux
  USE INITIALVELOCITY; USE ENVIRPMOD
  USE BIOENERGETICS
  USE BUILDVERSION
  IMPLICIT NONE
 ! include "omp_lib.h"      ! OPENMP directive to adjust the # of processors TOGGLE FOR DEBUG

  EXTERNAL RESTART_OUTPUT
  TYPE (DIALOG) :: DLG
  INTEGER       :: RESULT         !, RESULT1, IRESULT   ! SW 2/2019
  CHARACTER(240):: MODDIR1  
  REAL          :: DEPTH
  !INTEGER                                   :: N_WAITS, NWAIT                                                           !SR 11/26/19
  !INTEGER,        ALLOCATABLE, DIMENSION(:) :: WAIT_INDEX                                                               !SR 11/26/19
  !CHARACTER(2),   ALLOCATABLE, DIMENSION(:) :: WAIT_TYPE                                                                !SR 11/26/19
  !CHARACTER(240), ALLOCATABLE, DIMENSION(:) :: FILEDIR
  LOGICAL       :: CSVFORMAT
  CHARACTER(30) :: CHAR30
  CHARACTER(8)  :: CHAR8

!***********************************************************************************************************************************
!**                                                       Task 1: Inputs                                                          **
!***********************************************************************************************************************************

INTEGER(4) length,istatus
character*255 dirc
!  call omp_set_num_threads(4)   ! set # of processors to NPROC  Moved to INPUT subroutine

IF(END_RUN.or.ERROR_OPEN)STOP    ! SW 6/26/15 3/18/16 Added code to prevent a thread from reinitializing output files as dialog box is closing...intermittant error Updated 8/23/2017

CALL GET_COMMAND_ARGUMENT(1,DIRC,LENGTH,ISTATUS)  
DIRC=TRIM(DIRC)

! IF(ISTATUS.NE.0)WRITE(*,*)'GET_COMMAND_ARGUMENT FAILED: STATUS=',ISTATUS

IF(LENGTH /= 0)THEN
    ISTATUS=CHDIR(DIRC)
    SELECT CASE(ISTATUS)
      CASE(2)  ! ENOENT
        WRITE(W2ERR,*)'The directory does not exist:',DIRC
        WRITE(W2ERR,*)'Run stopped'
      CASE(20)   ! ENOTDIR
        WRITE(W2ERR,*)'This is not a directory:', DIRC
        WRITE(W2ERR,*)'Run stopped'
      CASE(0)    ! NO ERROR
    END SELECT  
ENDIF

MODDIR = FILE$CURDRIVE              !  GET CURRENT DIRECTORY
LENGTH = GETDRIVEDIRQQ(MODDIR)

OPEN(CON,FILE='W2CodeCompilerVersion.opt',status='unknown')
write(CON,'(A,F5.2)')' CE-QUAL-W2 Version #:',W2VER
write(CON,*)'Compiler Version and Code Compile Date'
write(CON,*)'INTEL_COMPILER_VERSION:',INTEL_COMPILER_VERSION
write(CON,*)'INTEL_COMPILER_BUILD_DATE:',INTEL_COMPILER_BUILD_DATE
write(CON,*)'CE-QUAL-W2 Version compile date:',BUILDTIME
close(CON)


! open ancillary control file for fish bioenergetics !mlm bioexp
  BIOEXP = .FALSE.              ! INITIALIZE LOGICAL VARIABLE THIS IS READ IN THE CONTROL FILE
  FISHBIO= .FALSE.
  INQUIRE(FILE='W2_con_anc.npt',EXIST=FISHBIO)     ! SW 5/26/15
  IF(FISHBIO)THEN
  
  OPEN(FISHBIOFN,FILE='W2_con_anc.npt',status='old')
   DO II = 1,16
    READ(FISHBIOFN,'(A8)') BIOC ! DUMMY VARIABLE AT THIS POINT FIX THIS
   END DO
  ENDIF
  
  FISH_PARTICLE_EXIST=.FALSE.
  INQUIRE(FILE='w2_particle.csv',EXIST=FISH_PARTICLE_EXIST)    ! SW 4/30/15
  
! Open control file
  IOPENFISH=0
  OPEN (CON,FILE=CONFN,STATUS='OLD',IOSTAT=I)
  IF (I /= 0) THEN
    CLOSE(CON)
    TEXT = 'Could not open w2_con.npt'
    CONFN='w2_con.csv'
    OPEN (CON,FILE=CONFN,STATUS='OLD',IOSTAT=I)
       IF (I /= 0) THEN
       TEXT = 'Could not open w2_con.npt and w2_con.csv'
       GO TO 240
       ENDIF
       TEXT = 'Found w2_con.csv'
  END IF

CALL INPUT

!SP CEMA
Call CEMA_W2_Input
!End SP CEMA

! Open Error File
  OPEN (W2ERR,FILE='w2.err',STATUS='UNKNOWN')

  ! DYN PIPE ADJUSTMENT CODE
   INQUIRE(FILE='w2_dynpipe_adjust.npt',EXIST=DYNPIPEADJUST)     ! SW 5/26/15
   DYNPAD='OF'  ! IT IS ONLY CHECKED IF IT IS ='ON' 2 CHARACTER VARIABLE
   DYNPAD_PIPE=1
  IF(DYNPIPEADJUST)THEN
  OPEN(CON,FILE='w2_dynpipe_adjust.npt',status='old')
    READ(CON,*) ! SKIP LINE
    READ(CON,*) ! SKIP LINE
    READ(CON,*)DYNPAD
    READ(CON,*)! SKIP LINE
    READ(CON,*)DYNPAD_PIPE
    READ(CON,*) ! SKIP LINE
    READ(CON,*)DYNPAD_SEG
    READ(CON,*) ! SKIP LINE
    READ(CON,*)DYNPAD_WL
    READ(CON,*) ! SKIP LINE
    READ(CON,*)DYNPAD_MAXRATE
    READ(CON,*) ! SKIP LINE
    READ(CON,*)DYNPAD_PERCENTCHANGE
    OPEN(DYNPIPELOG,FILE='dynpipe_adjustment_log.csv',status='unknown')
    WRITE(DYNPIPELOG,'(A)')'JDAY,BP(DYNPAD_PIPE),Z(DYNPAD_SEG),SZ(DYNPAD_SEG),DLT,(Z(DYNPAD_SEG)-SZ(DYNPAD_SEG))/DLT'
    CLOSE(CON)  
  ENDIF
  
  ! END DYN PIPE ADJUSTMENT CODE SW 2/18/2020
  
  
! Read multiple seperate waterbody file  2/9/2019 SW
WAIT_FOR_INFLOW_RESULTS=.FALSE.
MWB_EXIST= .FALSE.   ! USING OLD FILE NAME
INQUIRE(FILE='w2_multiple_WB.npt',EXIST=MWB_EXIST)     ! SW 5/26/15
IF(MWB_EXIST)THEN
OPEN(CON,FILE='w2_multiple_WB.npt',STATUS='OLD')
READ(CON,*)
READ(CON,*)DEG   ! USING OLD CHARACTER NAME
IF(DEG=='ON')THEN
!--------------
!new variables:
!--------------
!  N_WAITS    -- integer specifying the number of input file sets to await
!  NWAIT      -- integer counter index for do loops
!  WAIT_TYPE  -- character array holding the type of input file we're awaiting ('BR' or 'TR')
!  WAIT_INDEX -- integer array holding the branch or tributary index for a set of files we're awaiting
!  FILEDIR    -- character array to hold the directory names of the awaited files
    
    WAIT_FOR_INFLOW_RESULTS=.TRUE.
    READ (CON,*)                                                                                                        !SR 11/26/19
    READ (CON,*) N_WAITS                                                                                                !SR 11/26/19
    IF (N_WAITS.GT.0) THEN                                                                                              !SR 11/26/19
      ALLOCATE (WAIT_TYPE(MAX(1,N_WAITS)), WAIT_INDEX(MAX(1,N_WAITS)), FILEDIR(MAX(1,N_WAITS)))                         !SR 11/26/19
    ELSE                                                                                                                !SR 11/26/19
      WRITE (W2ERR,'(A)') 'ERROR-- Number of inputs to wait for must be greater than zero in multiple_WB.npt file.'             !SR 11/26/19
      STOP                                                                                                              !SR 11/26/19
    END IF                                                                                                              !SR 11/26/19
    READ (CON,*)                                                                                                        !SR 11/26/19
    DO NWAIT=1,N_WAITS                                                                                                  !SR 11/26/19
      READ (CON,*) WAIT_TYPE(NWAIT), WAIT_INDEX(NWAIT), FILEDIR(NWAIT)                                                  !SR 11/26/19
    END DO                                                                                                              !SR 11/26/19
    READ (CON,*)
    READ (CON,*) TIME_BUFFER
    IF (TIME_BUFFER < 0) THEN                                                                                           !SR 11/26/19
      WRITE (W2ERR,'(A)') 'ERROR-- Buffer time (days of data to await) in multiple_WB.npt cannot be less than zero.'            !SR 11/26/19
      STOP                                                                                                              !SR 11/26/19
    END IF                                                                                                              !SR 11/26/19
    READ (CON,*)
    READ (CON,*) WAIT_TIME
    IF (WAIT_TIME <= 0) THEN                                                                                            !SR 11/26/19
      WRITE (W2ERR,'(A)') 'ERROR-- Delay time in multiple_WB.npt must be greater than zero.'                                    !SR 11/26/19
      STOP                                                                                                              !SR 11/26/19
    END IF                                                                                                              !SR 11/26/19
    OPEN (9911,FILE='WaitForRunLog.opt',STATUS='unknown')         !moved: no need to create file if not used            !SR 11/26/19
  END IF
  CLOSE (CON)
END IF
!-----------------------------------------------------------------------------------------------------
!Setting up some variables so that the program knows which inputs to wait for... just before CALL INIT
!-----------------------------------------------------------------------------------------------------
  ! Set up some variables so that the program knows which inputs to wait for and where to find them                     !SR 11/26/19
  WAIT_FOR_TRIB_INPUT   = .FALSE.                                                                                       !SR 11/26/19
  WAIT_FOR_BRANCH_INPUT = .FALSE.                                                                                       !SR 11/26/19
  IF (WAIT_FOR_INFLOW_RESULTS) THEN
    DO NWAIT=1,N_WAITS                                                                                                  !SR 11/26/19
      IF (WAIT_TYPE(NWAIT) == 'BR') THEN                                   ! BR for branch input                        !SR 11/26/19
        WAIT_FOR_BRANCH_INPUT(WAIT_INDEX(NWAIT)) = .TRUE.                                                               !SR 11/26/19
        BR_FILEDIR(WAIT_INDEX(NWAIT))            = FILEDIR(NWAIT)                                                       !SR 11/26/19
      ELSEIF (WAIT_TYPE(NWAIT) == 'TR') THEN                               ! TR for tributary input                     !SR 11/26/19
        WAIT_FOR_TRIB_INPUT(WAIT_INDEX(NWAIT))   = .TRUE.                                                               !SR 11/26/19
        TR_FILEDIR(WAIT_INDEX(NWAIT))            = FILEDIR(NWAIT)                                                       !SR 11/26/19
      END IF                                                                                                            !SR 11/26/19
    END DO                                                                                                              !SR 11/26/19
  END IF


  RESTART_IN   =  RSIC == '      ON'.OR. RESTART_PUSHED
! Restart data
  IF (RESTART_PUSHED) RSIFN = 'rso.opt'

  JDAY = TMSTRT
  IF (RESTART_IN) THEN
    RSI=10   ! SW 5/26/15
    VERT_PROFILE = .FALSE.
    LONG_PROFILE = .FALSE.
    OPEN  (RSI,FILE=RSIFN,FORM='UNFORMATTED',STATUS='OLD')
    READ  (RSI) NIT,    NV,     KMIN,   IMIN,   NSPRF,  CMBRT,  ZMIN,   IZMIN,  START,  CURRENT
    READ  (RSI) DLTDP,  SNPDP,  TSRDP,  VPLDP,  PRFDP,  CPLDP,  SPRDP,  RSODP,  SCRDP,  FLXDP,  WDODP
    READ  (RSI) JDAY,   ELTM,   ELTMF,  DLT,    DLTAV,  DLTS,   MINDLT, JDMIN,  CURMAX
    READ  (RSI) NXTMSN, NXTMTS, NXTMPR, NXTMCP, NXTMVP, NXTMRS, NXTMSC, NXTMSP, NXTMFL, NXTMWD, NXWL, NXFLOWBAL, NXNPBAL,NXTMWD_SEC
    READ  (RSI) VOLIN,  VOLOUT, VOLUH,  VOLDH,  VOLPR,  VOLTRB, VOLDT,  VOLWD,  VOLEV,  VOLSBR, VOLTR, VOLSR,VOLICE,ICEBANK
    READ  (RSI) TSSEV,  TSSPR,  TSSTR,  TSSDT,  TSSWD,  TSSIN,  TSSOUT, TSSS,   TSSB,   TSSICE
    READ  (RSI) TSSUH,  TSSDH,  TSSUH2, TSSDH2, CSSUH2, CSSDH2, VOLUH2, VOLDH2, QUH1
    READ  (RSI) ESBR,   ETBR,   EBRI
    READ  (RSI) Z,      SZ,     ELWS,   SAVH2,  SAVHR,  H2
    READ  (RSI) KTWB,   KTI,    SKTI,   SBKT
    READ  (RSI) ICE,    ICETH,  CUF,    QSUM
    READ  (RSI) U,      W,      SU,     SW,     AZ,     SAZ,    DLTLIM
    READ  (RSI) T1,     T2,     C1,     C2,     C1S,    SED,    KFS,    CSSK
    READ  (RSI) EPD,    EPM
    READ  (RSI) MACMBRT,MACRC,  SMACRC, MAC,    SMAC,   MACRM,  MACSS
    READ  (RSI) SEDC, SEDN, SEDP, ZOO, CD  ! mlm 10/06
    READ  (RSI) SDKV                       ! cb 11/30/06
    READ  (RSI) TKE                        ! sw 10/4/07
    READ  (RSI) BR_INACTIVE,WARNING_OPEN                ! SW 8/1/2018
    if(envirpc == '      ON')THEN
    
      allocate(cc_e(NCT),c_int(NCT),c_top(NCT),cd_e(NDC),cd_int(NDC),cd_top(NDC),c_avg(NCT),cd_avg(NDC),cn_e(NCT),cdn_e(NDC))
      cc_e='   '
      c_int=0.0
      c_top=0.0
      cd_e='   '
      cd_int=0.0
      cd_top=0.0
      c_avg=0.0
      cd_avg=0.0
      cn_e=0.0
      cdn_e=0.0
      NAC_E=0
      NACD_E=0
      OPEN(CONE,file='w2_envirprf.npt',status='old')
      
     CSVFORMAT=.FALSE.
     READ(CONE,'(//A)')CHAR30
     DO J=1,30
         IF(CHAR30(J:J)==',')THEN
             CSVFORMAT=.TRUE.
             EXIT
         ENDIF
     ENDDO
     REWIND(CONE)

      IF(CSVFORMAT)THEN
        READ(CONE,*)
        READ(CONE,*)
        READ (CONE,*) I_SEGINT,numclass,selectivec,sjday1,sjday2,istart(1),iend(1),(istart(I),iend(I),I=2,I_SEGINT)
        SELECTIVEC=ADJUSTR(SELECTIVEC)
        IF(I_SEGINT==0)I_SEGINT=1
        READ(CONE,*)
        READ(CONE,*)
        Read (CONE,*) VEL_VPR,VEL_VPR, VEL_INT, VEL_TOP,TEMP_VPR,TEMP_INT,TEMP_TOP, depth_vpr,d_int,d_top
        VEL_VPR=ADJUSTR(VEL_VPR);TEMP_VPR=ADJUSTR(TEMP_VPR);DEPTH_VPR=ADJUSTR(DEPTH_VPR)
        READ(CONE,*)
        READ(CONE,*)
        DO JC=1,NCT
        READ (CONE,*) CHAR8, CC_E(JC), C_INT(JC), C_TOP(JC)
        CC_E(JC)=ADJUSTR(CC_E(JC))
        ENDDO
        READ(CONE,*)
        READ(CONE,*)
        DO JD=1,NDC
        READ (CONE,*) CHAR8, CD_E(JD),CD_INT(JD), CD_TOP(JD)
        CD_E(JD)=ADJUSTR(CD_E(JD))        
        ENDDO
      ELSE    
      READ (CONE,'(//I1,7X,I8,5x,a3,f8.0,f8.0,9(i8,i8))') I_SEGINT,numclass,selectivec,sjday1,sjday2,istart(1),iend(1),(istart(I),iend(I),I=2,I_SEGINT)
      IF(I_SEGINT==0)I_SEGINT=1
      Read (CONE,'(//8x,3(5x,a3,f8.3,f8.3))') VEL_VPR, VEL_INT, VEL_TOP,TEMP_VPR,TEMP_INT,TEMP_TOP, depth_vpr,d_int,d_top
      READ (CONE,'(//(8X,(5X,A3,F8.0,F8.0)))') (CC_E(JC), C_INT(JC), C_TOP(JC), JC=1,NCT)
      READ (CONE,'(//(8X,(5X,A3,F8.0,F8.0)))') (CD_E(JD),CD_INT(JD), CD_TOP(JD), JD=1,NDC)
      ENDIF
      CLOSE(CONE)

          DO JC=1,NCT
          IF (CC_E(JC).EQ.' ON') THEN
            NAC_E     = NAC_E+1
            CN_E(NAC_E) = JC
          END IF
          End DO
          DO JD=1,NDC
          IF (CD_E(JD).EQ.' ON') THEN
            NACD_E     = NACD_E+1
            CDN_E(NACD_E) = JD
          END IF
          End DO

    allocate (c_cnt(I_SEGINT,NCT),cd_cnt(I_SEGINT,NDC),c_class(I_SEGINT,NCT,numclass),cd_class(I_SEGINT,NDC,numclass),c_tot(I_SEGINT,NCT),cd_tot(I_SEGINT,NDC),t_class(I_SEGINT,numclass),v_class(I_SEGINT,numclass),c_sum(NCT),cd_sum(NDC))
    allocate (conc_c(NCT,numclass),conc_cd(NDC,numclass))
    allocate(d_class(I_SEGINT,numclass))
    ALLOCATE (D_TOT(I_SEGINT),D_CNT(I_SEGINT),T_TOT(I_SEGINT),T_CNT(I_SEGINT))
    ALLOCATE(V_TOT(I_SEGINT),V_CNT(I_SEGINT),VOLGL(I_SEGINT),SUMVOLT(I_SEGINT))
   
    READ(RSI)T_CLASS,V_CLASS,C_CLASS,CD_CLASS,T_TOT,T_CNT,SUMVOLT,V_CNT,V_TOT,C_TOT,C_CNT,CD_TOT,CD_CNT
    
    ENDIF
    IF(NPI > 0)READ(RSI)YS,VS,VST,YST,DTP,QOLD
    READ(RSI)TPOUT,TPTRIB,TPDTRIB,TPWD,TPPR,TPIN,TP_SEDSOD_PO4,PFLUXIN,TNOUT,TNTRIB,TNDTRIB,TNWD,TNPR,TNIN,TN_SEDSOD_NH4,NFLUXIN,ATMDEP_P,ATMDEP_N,NH3GASLOSS     ! TP_SEDBURIAL,TN_SEDBURIAL,
    ! SEDIMENT DIAGENESIS RESTART SW 2/2019
    IF(IncludeCEMASedDiagenesis)THEN
    READ(RSI)C2SF,CellArea,BedPorosity
    IF(Bubbles_Calculation) THEN
    READ(RSI)TConc,SConc,CrackOpen,BubbleRelWB,MFTBubbReleased,GasReleaseCH4
    END IF   
    ENDIF
    
    CLOSE (RSI)
  END IF
  CE_QUAL_W2 =  1
  
  ! Open warning file

  IF(.NOT.WARNING_OPEN)THEN
      OPEN (WRN,FILE='w2.wrn',STATUS='UNKNOWN')
  ELSE
      OPEN (WRN,FILE='w2.wrn',POSITION='APPEND')
      WRITE(WRN,*)'***RESTART*** APPENDING ON JDAY',JDAY
  ENDIF
  
CALL INIT


! determining initial horizontal velocities and water levels
    once_through=.true.
    IF(inituwl == '      ON')init_vel=.true.    
    if(.not. restart_in)then       
      if(init_vel)then
        allocate (qssi(imx),loop_branch(nbr),elwss(imx),uavg(imx))
        elwss=elws
        call initial_water_level
        b=bsave
        call initgeom
        call initial_u_velocity       
        open(NUNIT,file='init_wl_u_check.dat',status='unknown')
        write(NUNIT,'("       i elws_calc    qssi       u   depth elws_init")')
        DO JW=1,NWB        
          DO JB=BS(JW),BE(JW)
            IU = CUS(JB)
            ID = DS(JB)
            do i=iu,id
              depth=elws(i)-el(kbi(i)+1,i)
              write(NUNIT,'(i8,f8.3,f8.2,f8.3,2f8.2)')i,elws(i),qssi(i),u(kt,i),depth,elwss(i)
            end do
          end do
        end do
        close(NUNIT)
        deallocate (qssi,loop_branch,elwss)
      end if    
    end if

  IF (.NOT. RESTART_IN) THEN
    LINE    = CCTIME(1:2)//':'//CCTIME(3:4)//':'//CCTIME(5:6)
    RESULT  = DLGSET (DLG,STARTING_TIME,TRIM(LINE))                                                   !Display starting time
    RESULT  = DLGSET (DLG,STATUS,'Executing')                                                         !Display execution status
    CURRENT = 0.0
  ELSE
    CALL CPU_TIME (CURRENT)
  END IF

  CALL OUTPUTINIT
    IF (RESTART_IN) THEN
    DO JW=1,NWB
      IF (SCREEN_OUTPUT(JW))CALL SCREEN_UPDATE(DLG)
    ENDDO
    ENDIF

  IF (.NOT. RESTART_IN) CALL CPU_TIME (START)

  if (macrophyte_on.and.constituents) call porosity  
  IF(SELECTC == '      ON')CALL SELECTIVEINIT   ! new subroutine for selecting water temperature target
  IF(SELECTC == '    USGS')CALL SELECTIVEINITUSGS   ! new subroutine for selecting water temperature target
  IF (TDGTA) CALL InitTDGtarget                  ! tdgtarget - initial
  IF(AERATEC == '      ON' .and. oxygen_demand)CALL AERATE
    If(CEMARelatedCode .and. IncludeBedConsolidation)Call SetupCEMASedimentModel
    If(IncludeFFTLayer)Call CEMAFFTLayerCode
    !If(CEMARelatedCode .and. IncludeCEMASedDiagenesis)Call CEMASedimentDiagenesis  
 
!***********************************************************************************************************************************
!**                                                   Task 2: Calculations                                                        **
!***********************************************************************************************************************************
  DO WHILE (.NOT. END_RUN.AND. .NOT. STOP_PUSHED)    
    IF (JDAY >= NXTVD) CALL READ_INPUT_DATA (NXTVD)
    CALL INTERPOLATE_INPUTS
    DLTTVD = (NXTVD-JDAY)*DAY
    DLT    =  MIN(DLT,DLTTVD+1.0)
    DLTS1  =  DLT
    IF (DLT <= DLTTVD+0.999) THEN
      DLTS = DLT
    ELSE
      KLOC = 1
      ILOC = 1
    END IF

    ! update wind at 2m for evaopration and evaoprative heat flux computations  ! SW 5/21/15
   DO JW=1,NWB
    DO I=CUS(BS(JW)),DS(BE(JW))
      WIND2(I) = WIND(JW)*WSC(I)*DLOG(2.0D0/Z0(JW))/DLOG(WINDH(JW)/Z0(JW))    
    END DO
   ENDDO
    
210 continue   ! timestep violation entry point
 IF(SELECTC == '      ON')CALL SELECTIVE   ! new subroutine for selecting water temperature target
 IF(SELECTC == '    USGS')CALL SELECTIVEUSGS   ! new subroutine for selecting water temperature target
CALL HYDROINOUT

!SP CEMA
  !if(sediment_diagenesis)then
    If(CEMARelatedCode .and. IncludeBedConsolidation)Call ComputeCEMARelatedSourceSinks
  !  If(CEMARelatedCode .and. IncludeCEMASedDiagenesis)Call ComputeCEMADiagenesisSourceSinks  SW 6/27/2017
  !end if
!End SP CEMA

!***********************************************************************************************************************************
!**                                           Task 2.2: Hydrodynamic calculations                                                 **
!***********************************************************************************************************************************
!!$OMP PARALLEL DO default(Private)   !(KT,IU,ID,IUT,IDT,K,I,JB,ZB,WWT,DFC,JJB,BETABR,GC2,HRAD,UDR,UDL,AB)
    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        IF(BR_INACTIVE(JB))CYCLE    ! SW 6/12/17
        IU = CUS(JB)
        ID = DS(JB)

!***********************************************************************************************************************************
!**                                Task 2.2.1: Boundary concentrations, temperatures, and densities                               **
!***********************************************************************************************************************************

        IUT = IU
        IDT = ID
        IF (UP_FLOW(JB)) THEN
          IF (.NOT. INTERNAL_FLOW(JB)) THEN
            DO K=KT,KB(IU)
                IF (QIND(JB)+QINSUM(JB).GT.0.0) THEN
                  TIN(JB)               = (TINSUM(JB)               *QINSUM(JB)+TIND(JB)          *QIND(JB))/(QIND(JB)+QINSUM(JB))
                  CIN(CN(1:NAC),JB)     =  MAX((CINSUM(CN(1:NAC),JB)*QINSUM(JB)+CIND(CN(1:NAC),JB)*QIND(JB))/(QIND(JB)+QINSUM(JB)),&
                                                0.0)
                  T1(K,IU-1)            =  TIN(JB)
                  T2(K,IU-1)            =  TIN(JB)
                  C1S(K,IU-1,CN(1:NAC)) =  CIN(CN(1:NAC),JB)
                  QIN(JB)               =  QIND(JB)+QINSUM(JB)
                ELSE
                  QIN(JB)               =  0.0
                  TIN(JB)               =  TIND(JB)
                  T1(K,IU-1)            =  TIND(JB)
                  T2(K,IU-1)            =  TIND(JB)
                  C1S(K,IU-1,CN(1:NAC)) =  CIND(CN(1:NAC),JB)
                END IF
            END DO
          ELSE IF (.NOT. DAM_INFLOW(JB)) THEN                                                                          !TC 08/03/04
            IF (JBUH(JB) >= BS(JW) .AND. JBUH(JB) <= BE(JW)) THEN
              TIN(JB)           = T1(KT,UHS(JB))
              CIN(CN(1:NAC),JB) = MAX(C1S(KT,UHS(JB),CN(1:NAC)),0.0)
              DO K=KT,KB(IU)    !CONCURRENT(K=KT:KB(IU))   !FORALL(K=KT:KB(IU))           !DO K=KT,KB(IU)
                T1(K,IU-1)            = T1(K,UHS(JB))
                T2(K,IU-1)            = T1(K,UHS(JB))
                C1S(K,IU-1,CN(1:NAC)) = C1S(K,UHS(JB),CN(1:NAC))
                C1(K,IU-1,CN(1:NAC))  = C1S(K,UHS(JB),CN(1:NAC))
                C2(K,IU-1,CN(1:NAC))  = C1S(K,UHS(JB),CN(1:NAC))
              END DO                      
            ELSE
              CALL UPSTREAM_WATERBODY
              TIN(JB)           = T1(KT,IU-1)
              CIN(CN(1:NAC),JB) = MAX(C1(KT,IU-1,CN(1:NAC)),0.0)
            END IF
          ELSE
            TIN(JB)           = TINSUM(JB)
            QIN(JB)           = QINSUM(JB)
            CIN(CN(1:NAC),JB) = MAX(CINSUM(CN(1:NAC),JB),0.0)
            DO K=KT,KB(ID)
              T1(K,IU-1)            = TIN(JB)
              T2(K,IU-1)            = TIN(JB)
              C1S(K,IU-1,CN(1:NAC)) = CIN(CN(1:NAC),JB)
            END DO
           END IF
        END IF
        IF (DN_FLOW(JB)) THEN
          DO K=KT,KB(ID)
              T1(K,ID+1)            = T2(K,ID)
              T2(K,ID+1)            = T2(K,ID)
              C1S(K,ID+1,CN(1:NAC)) = C1S(K,ID,CN(1:NAC))
          END DO
        END IF
        IF (UP_HEAD(JB)) THEN
          IUT = IU-1
          IF (UH_INTERNAL(JB)) THEN
            IF (JBUH(JB) >= BS(JW) .AND. JBUH(JB) <= BE(JW)) THEN
              DO K=KT,KB(IUT)
                RHO(K,IUT)           = RHO(K,UHS(JB))
                T1(K,IUT)            = T2(K,UHS(JB))
                T2(K,IUT)            = T2(K,UHS(JB))
                C1S(K,IUT,CN(1:NAC)) = C1S(K,UHS(JB),CN(1:NAC))
                C1(K,IUT,CN(1:NAC))  = C1S(K,UHS(JB),CN(1:NAC))
                C2(K,IUT,CN(1:NAC))  = C1S(K,UHS(JB),CN(1:NAC))
              END DO
            ELSE
              CALL UPSTREAM_WATERBODY
            END IF
            DO K=KT,KB(IUT)
              RHO(K,IUT) = DENSITY(T2(K,IUT),DMAX1(TDS(K,IUT),0.0D0),DMAX1(TISS(K,IUT),0.0D0))
            END DO
          ELSE IF (UH_EXTERNAL(JB)) THEN
            DO K=KT,KB(IUT)
              RHO(K,IUT)           = DENSITY(TUH(K,JB),DMAX1(TDS(K,IUT),0.0D0),DMAX1(TISS(K,IUT),0.0D0))
              T1(K,IUT)            = TUH(K,JB)
              T2(K,IUT)            = TUH(K,JB)
              C1S(K,IUT,CN(1:NAC)) = CUH(K,CN(1:NAC),JB)
              C1(K,IUT,CN(1:NAC))  = CUH(K,CN(1:NAC),JB)
              C2(K,IUT,CN(1:NAC))  = CUH(K,CN(1:NAC),JB)
            END DO
          END IF
        END IF
        IF (DN_HEAD(JB)) THEN
          IDT = ID+1
          IF (DH_INTERNAL(JB)) THEN
            IF (JBDH(JB) >= BS(JW) .AND. JBDH(JB) <= BE(JW)) THEN
             DO K=KT,KB(IDT)
                RHO(K,IDT)           = RHO(K,DHS(JB))
                T1(K,IDT)            = T2(K,DHS(JB))
                T2(K,IDT)            = T2(K,DHS(JB))
                C1S(K,IDT,CN(1:NAC)) = C1S(K,DHS(JB),CN(1:NAC))
                C1(K,IDT,CN(1:NAC))  = C1S(K,DHS(JB),CN(1:NAC))
                C2(K,IDT,CN(1:NAC))  = C1S(K,DHS(JB),CN(1:NAC))
              END DO
            ELSE
              CALL DOWNSTREAM_WATERBODY
            END IF
            DO K=KT,KB(ID)
              RHO(K,IDT) = DENSITY(T2(K,IDT),DMAX1(TDS(K,IDT),0.0D0),DMAX1(TISS(K,IDT),0.0D0))
            END DO
          ELSE IF (DH_EXTERNAL(JB)) THEN
            DO K=KT,KB(IDT)
              RHO(K,IDT)           = DENSITY(TDH(K,JB),DMAX1(TDS(K,IDT),0.0D0),DMAX1(TISS(K,IDT),0.0D0))
              T1(K,IDT)            = TDH(K,JB)
              T2(K,IDT)            = TDH(K,JB)
              C1S(K,IDT,CN(1:NAC)) = CDH(K,CN(1:NAC),JB)
              C1(K,IDT,CN(1:NAC))  = CDH(K,CN(1:NAC),JB)
              C2(K,IDT,CN(1:NAC))  = CDH(K,CN(1:NAC),JB)
            END DO
          END IF
        END IF

!***********************************************************************************************************************************
!**                                                 Task 2.2.2: Momentum terms                                                    **
!***********************************************************************************************************************************

!****** Density pressures

        DO I=IUT,IDT
          DO K=KT,KB(I)
            P(K,I) = P(K-1,I)+RHO(K,I)*G*H(K,JW)*COSA(JB)
          END DO
        END DO

!****** Horizontal density gradients

        DO I=IUT,IDT-1
          HDG(KT,I) = DLXRHO(I)*(BKT(I)+BKT(I+1))*0.5D0*H(KT,JW)*(P(KT,I+1)-P(KT,I))
          DO K=KT+1,KBMIN(I)
            HDG(K,I) = DLXRHO(I)*BHR2(K,I)*((P(K-1,I+1)-P(K-1,I))+(P(K,I+1)-P(K,I)))
          END DO
        END DO

!****** Adjusted wind speed and surface wind shear drag coefficient

        DO I=IU-1,ID+1
          WIND10(I) = WIND(JW)*WSC(I)*DLOG(10.0D0/Z0(JW))/DLOG(WINDH(JW)/Z0(JW))     ! older  version z0=0.01                      ! SW 11/28/07
          FETCH(I)  = FETCHD(I,JB)
          IF (COS(PHI(JW)-PHI0(I)) < 0.0) FETCH(I) = FETCHU(I,JB)
          IF (FETCH(I) <= 0.0) FETCH(I) = DLX(I)
          IF (FETCH_CALC(JW)) THEN
            ZB        = 0.8D0*DLOG(FETCH(I)*0.5D0)-1.0718D0
            WIND10(I) = WIND10(I)*(5.0D0*ZB+4.6052D0)/(3.0D0*ZB+9.2103D0)
          END IF
          
          IF(WIND10(I) >= 15.0)THEN                     ! SW 1/19/2008
          CZ(I) = 0.0026D0
          ELSEIF(WIND10(I) >= 4.0)THEN
          CZ(I) = 0.0005D0*DSQRT(WIND10(I)) 
          ELSEIF(WIND10(I) >= 0.5)THEN
          CZ(I)= 0.0044D0*WIND10(I)**(-1.15D0)
          ELSE
          CZ(I)= 0.01D0
          ENDIF
          
  !        CZ(I) = 0.0
  !        IF (WIND10(I) >= 1.0)  CZ(I) = 0.0005*SQRT(WIND10(I))
  !        IF (WIND10(I) >= 4.0) CZ(I) = 0.0005*SQRT(WIND10(I))          
  !        IF (WIND10(I) >= 15.0) CZ(I) = 0.0026
        END DO

!****** Longitudinal and lateral surface wind shear and exponential decay

        DO I=IUT,IDT-1
          !WSHX(I) = CZ(I)*WIND10(I)**2*RHOA/RHOW*    DCOS(PHI(JW)-PHI0(I))* ICESW(I)
          !WSHY(I) = CZ(I)*WIND10(I)**2*RHOA/RHOW*DABS(DSIN(PHI(JW)-PHI0(I)))*ICESW(I)
          WSHX(I) = CZ(I)*WIND10(I)*WIND10(I)*RHOA/RHOW*    DCOS(PHI(JW)-PHI0(I))* ICESW(I)    ! SW 4/20/16 SPEED
          WSHY(I) = CZ(I)*WIND10(I)*WIND10(I)*RHOA/RHOW*DABS(DSIN(PHI(JW)-PHI0(I)))*ICESW(I)
          WWT     = 0.0
          IF (WIND10(I) /= 0.0) WWT = 6.95D-2*(FETCH(I)**0.233D0)*WIND10(I)**0.534D0
          DFC = -8.0D0*PI*PI/(G*WWT*WWT+NONZERO)
          DO K=KT,KBMIN(I)
            DECAY(K,I) = DEXP(DMAX1(DFC*DEPTHB(K,I),-30.0D0))
          END DO

!******** Branch inflow lateral shear and friction

          DO JJB=1,NBR
            IF(BR_INACTIVE(JJB))CYCLE  ! SW 6/12/2017
            IF (I == UHS(JJB) .AND. .NOT. INTERNAL_FLOW(JJB)) THEN
              BETABR = (PHI0(I)-PHI0(US(JJB)))
              IF (JJB >= BS(JW) .AND. JJB <= BE(JW)) THEN
                DO K=KT,KBMIN(I)
                  IF (U(K,US(JJB)) < 0.0) THEN
                    UXBR(K,I) = UXBR(K,I)+ABS(U(K,US(JJB)))*DCOS(BETABR)     *VOLUH2(K,JJB)/(DLT*DLX(I))
                    UYBR(K,I) = UYBR(K,I)              +ABS(DSIN(BETABR))*ABS(VOLUH2(K,JJB))/DLT
                  END IF
                END DO
              ELSE
                CALL UPSTREAM_BRANCH
              END IF
            END IF
            IF (I == DHS(JJB)) THEN
              BETABR = (PHI0(I)-PHI0(DS(JJB)))
              IF (I == US(JB) .AND. UHS(JB) /= DS(JJB)) THEN
                IF (JJB >= BS(JW) .AND. JJB <= BE(JW)) THEN
                  DO K=KT,KBMIN(I)
                    IF (U(K,DS(JJB)) >= 0.0) THEN
                      UXBR(K,I) = UXBR(K,I)+U(K,DS(JJB))*   DCOS(BETABR) *VOLDH2(K,JJB)/(DLT*DLX(I))
                      UYBR(K,I) = UYBR(K,I)            +ABS(DSIN(BETABR))*VOLDH2(K,JJB)/DLT
                    END IF
                  END DO
                ELSE
                  CALL DOWNSTREAM_BRANCH
                END IF
              ELSE IF (I /= US(JB)) THEN
                IF (JJB >= BS(JW) .AND. JJB <= BE(JW)) THEN
                  DO K=KT,KBMIN(I)
                    IF (U(K,DS(JJB)) >= 0.0) THEN
                      UXBR(K,I) = UXBR(K,I)+U(K,DS(JJB))*   DCOS(BETABR) *VOLDH2(K,JJB)/(DLT*DLX(I))
                      UYBR(K,I) = UYBR(K,I)            +ABS(DSIN(BETABR))*VOLDH2(K,JJB)/DLT
                    END IF
                  END DO
                ELSE
                  CALL DOWNSTREAM_BRANCH
                END IF
              END IF
            END IF
          END DO
          DO K=KT,KBMIN(I)
            FRICBR(K,I) = (FI(JW)/8.0D0)*RHO(K,I)*(UYBR(K,I)/(DLX(I)*H2(K,I)))**2
          END DO
        END DO

!****** Vertical eddy viscosities/diffusivities
        FIRSTI(JW) = IUT
		LASTI(JW) = IDT
        DO I=IUT,IDT-1
          CALL CALCULATE_AZ
          !SP CEMA
          if(sediment_diagenesis)then
            If(CEMARelatedCode .and. IncludeCEMASedDiagenesis .and. ApplyBubbTurb)Call CEMABubblesTurbulence
          end if
          !End SP CEMA
          IF (KBMIN(I) <= KT+1 .AND. KB(I) > KBMIN(I)) THEN
            AZ(KBMIN(I),I) = AZMIN
            DZ(KBMIN(I),I) = DZMIN
          END IF
        END DO
        IF (AZC(JW) == '     TKE'.OR.AZC(JW) == '    TKE1') THEN
          AZT(:,IDT-1)  = AZ(:,IDT-1)
          DO I=IUT,IDT-2
            DO K=KT,KBMIN(I)
              AZT(K,I)  = 0.5*(AZ(K,I)+AZ(K,I+1))
            END DO
          AZ(KBMIN(I),I) = AZMIN              !SG 10/4/07 SW 10/4/07
          END DO
          AZ(KT:KMX-1,IUT:IDT-1)=AZT(KT:KMX-1,IUT:IDT-1)
        END IF
        DO JWR=1,NIW
        IF (WEIR_CALC) AZ(KTWR(JWR)-1:KBWR(JWR),IWR(1:NIW)) = 0.0
        END DO

!****** Average eddy diffusivities

        IF(AZC(JW) /= '     TKE'.AND.AZC(JW) /= '    TKE1')THEN
        DZ(KT:KB(IDT)-1,IDT) = DZT(KT:KB(IDT)-1,IDT-1)    ! DZT is only used for non-TKE algorithms
        ELSE
        DZ(KT:KB(IDT)-1,IDT) = DZ(KT:KB(IDT)-1,IDT-1)
        ENDIF
        DO I=IUT,IDT-1
          DO K=KT,KB(I)-1
            IF (K >= KBMIN(I)) THEN
              IF (KB(I-1) >= KB(I) .AND. I /= IUT) THEN
                DZ(K,I) = DZ(K,I-1)
              ELSE
                DZ(K,I) = DZMIN
              END IF
            ELSE
              IF(AZC(JW) /= '     TKE'.AND.AZC(JW) /= '    TKE1')THEN
                 IF(I == IUT)THEN                             ! SW 10/20/07
                    DZ(K,I)=DZT(K,I)
                 ELSE
                    DZ(K,I) = (DZT(K,I)+DZT(K,I-1))*0.5D0        ! SW 10/20/07  (DZT(K,I)+DZT(K+1,I))*0.5 ! FOR NON-TKE ALGORITHMS, AVERAGE DZ FROM EDGES TO CELL CENTERS
                 ENDIF
              ENDIF
            END IF
          END DO
        END DO

! Hypolimnetic aeration

        IF(AERATEC == '      ON' .and. oxygen_demand)CALL DZAERATE

!****** Density inversions

        DO I=IUT,IDT
          DO K=KT,KB(I)-1
            DZQ(K,I) = MIN(1.0D-2,DZ(K,I))                                    !MIN(1.0E-4,DZ(K,I)) No reason to limit DZ in rivers/estuaries-used in ULTIMATE scheme
            IF (RHO(K,I) > RHO(K+1,I)) DZ(K,I) = DZMAX
          END DO
        END DO

!****** Wind, velocity, and bottom shear stresses @ top and bottom of cell

        SB(:,IUT:IDT-1) = 0.0
        DO I=IUT,IDT-1
          ST(KT,I) = WSHX(I)*BR(KTI(I),I)
          DO K=KT+1,KBMIN(I)
            ST(K,I) = WSHX(I)*DECAY(K-1,I)*BR(K,I)
            IF (.NOT. IMPLICIT_VISC(JW)) ST(K,I) = ST(K,I)+AZ(K-1,I)*(BR(K-1,I)+BR(K,I))*0.5D0*(U(K-1,I)-U(K,I))/((AVH2(K-1,I)       &
                                                   +AVH2(K-1,I+1))*0.5D0)
          END DO
          GC2 = 0.0
          IF (FRIC(I) /= 0.0) GC2 = G/(FRIC(I)*FRIC(I))

          HRAD=BHR2(KT,I)/(BR(KTI(I),I)-BR(KT+1,I)+2.0D0*AVHR(KT,I))
          IF(MACROPHYTE_ON.AND.MANNINGS_N(JW))THEN
            CALL MACROPHYTE_FRICTION(HRAD,FRIC(I),EFFRIC,KT,I)
            GC2=G*EFFRIC*EFFRIC/HRAD**0.33333333D0
          ELSE IF(.NOT.MACROPHYTE_ON.AND.MANNINGS_N(JW))THEN
            GC2=G*FRIC(I)*FRIC(I)/HRAD**0.33333333D0
          END IF
          IF (ONE_LAYER(I)) THEN
            SB(KT,I) = ST(KT+1,I)+GC2*(BR(KTI(I),I)+2.0D0*AVHR(KT,I))*U(KT,I)*DABS(U(KT,I))
          ELSE
            SB(KT,I) = GC2*(BR(KTI(I),I)-BR(KT+1,I)+2.0D0*AVHR(KT,I))*U(KT,I)*DABS(U(KT,I))
            DO K=KT+1,KBMIN(I)-1
              HRAD=(BHR2(K,I)/(BR(K,I)-BR(K+1,I)+2.0D0*H(K,JW)))
              IF(MACROPHYTE_ON.AND.MANNINGS_N(JW))THEN
                CALL MACROPHYTE_FRICTION(HRAD,FRIC(I),EFFRIC,K,I)
                GC2=G*EFFRIC*EFFRIC/HRAD**0.33333333D0
              ELSE IF(.NOT.MACROPHYTE_ON.AND.MANNINGS_N(JW))THEN
                GC2=G*FRIC(I)*FRIC(I)/HRAD**0.33333333D0
              END IF
              SB(K,I) = GC2*(BR(K,I)-BR(K+1,I)+2.0D0*H(K,JW))*U(K,I)*DABS(U(K,I))
            END DO
            IF (KT /= KBMIN(I)) THEN
              HRAD=(BHR2(KBMIN(I),I)/(BR(KBMIN(I),I)+2.0D0*H(KBMIN(I),JW)))
              IF(MACROPHYTE_ON.AND.MANNINGS_N(JW))THEN
                CALL MACROPHYTE_FRICTION(HRAD,FRIC(I),EFFRIC,KBMIN(I),I)
                GC2=G*EFFRIC*EFFRIC/HRAD**0.33333333D0
              ELSE IF(.NOT.MACROPHYTE_ON.AND.MANNINGS_N(JW))THEN
                GC2=G*FRIC(I)*FRIC(I)/HRAD**0.33333333D0
              END IF

              IF (KBMIN(I) /= KB(I)) THEN
                SB(KBMIN(I),I) = GC2*(BR(KBMIN(I),I)-BR(KBMIN(I)+1,I)+2.0D0*H2(K,I))*U(KBMIN(I),I)*DABS(U(KBMIN(I),I))
              ELSE
                SB(KBMIN(I),I) = GC2*(BR(KBMIN(I),I)+2.0D0*H2(K,I))*U(KBMIN(I),I)*DABS(U(KBMIN(I),I))
              END IF
            END IF
          END IF
          DO K=KT,KBMIN(I)-1
            SB(K,I) = SB(K,I)+ST(K+1,I)
          END DO
          SB(KBMIN(I),I) = SB(KBMIN(I),I)+WSHX(I)*DECAY(KBMIN(I),I)*(BR(KBMIN(I)-1,I)+BR(KBMIN(I),I))*0.5D0
        END DO

!****** Horizontal advection of momentum

        DO I=IU,ID-1
          DO K=KT,KBMIN(I)
            UDR       = (1.0D0+DSIGN(1.0D0,(U(K,I)+U(K,I+1))*0.5D0))*0.5D0
            UDL       = (1.0D0+DSIGN(1.0D0,(U(K,I)+U(K,I-1))*0.5D0))*0.5D0
            ADMX(K,I) = (BH2(K,I+1)*(U(K,I+1)+U(K,I))*0.5D0*(UDR*U(K,I)+(1.0-UDR)*U(K,I+1))-BH2(K,I)*(U(K,I)+U(K,I-1))               &
                        *0.5D0*(UDL*U(K,I-1)+(1.0D0-UDL)*U(K,I)))/DLXR(I)
          END DO
        END DO

!****** Horizontal dispersion of momentum

        DO I=IU,ID-1
          DO K=KT,KBMIN(I)
            IF(AX(JW) >= 0.0)THEN
            DM(K,I) = AX(JW)*(BH2(K,I+1)*(U(K,I+1)-U(K,I))/DLX(I+1)-BH2(K,I)*(U(K,I)-U(K,I-1))/DLX(I))/DLXR(I)
            ELSE
            DM(K,I) = ABS(U(K,I))*ABS(AX(JW))*H(K,JW)*(BH2(K,I+1)*(U(K,I+1)-U(K,I))/DLX(I+1)-BH2(K,I)*(U(K,I)-U(K,I-1))/DLX(I))/DLXR(I)     ! SW 8/2/2017 SCALE AX WITH U, FOR EXAMPLE AX=0.1U
            ENDIF
          END DO
        END DO

!****** Vertical advection of momentum

        DO I=IU,ID-1
          DO K=KT,KB(I)-1
            AB        = (1.0D0+DSIGN(1.0D0,(W(K,I+1)+W(K,I))*0.5D0))*0.5D0
            ADMZ(K,I) = (BR(K,I)+BR(K+1,I))*0.5D0*(W(K,I+1)+W(K,I))*0.5D0*(AB*U(K,I)+(1.0-AB)*U(K+1,I))
          END DO
        END DO

!****** Gravity force due to channel slope

        DO I=IU-1,ID
          GRAV(KT,I) = AVHR(KT,I)*(BKT(I)+BKT(I+1))*0.5D0*G*SINAC(JB)                                                
          DO K=KT+1,KB(I)                                                                                              
            GRAV(K,I) = BHR2(K,I)*G*SINAC(JB)
          END DO
        END DO

        IF(ICEC(JW)  == '    ONWB')THEN
        DO I = IU,ID            ! water loss due to ice formation or water gain due to ice melting
          VolIce(jb)=VolIce(jb)+iceqss(i)*dlt
          QSS(KT,I) = QSS(KT,I) + IceQSS(I)
          IceQSS(I) = 0.0d00
        END DO
        END IF


!***********************************************************************************************************************************
!**                                            Task 2.2.3: Water surface elevation                                                **
!***********************************************************************************************************************************

!****** Tridiagonal coefficients

        BHRHO(IU-1:ID+1) = 0.0D0; D(IU-1:ID+1) = 0.0D0; F(IU-1:ID+1) = 0.0D0
        DO I=IU,ID-1
          DO K=KT,KBMIN(I)
            BHRHO(I) = BHRHO(I)+(BH2(K,I+1)/RHO(K,I+1)+BH2(K,I)/RHO(K,I))
          END DO
          DO K=KT,KB(I)
            D(I) = D(I)+(U(K,I)*BHR2(K,I)-U(K,I-1)*BHR2(K,I-1)-QSS(K,I)+(UXBR(K,I)-UXBR(K,I-1))*DLT)
            F(I) = F(I)+(-SB(K,I)+ST(K,I)-ADMX(K,I)+DM(K,I)-HDG(K,I)+GRAV(K,I))
          END DO
        END DO

!****** Boundary tridiagonal coefficients

        D(IU) = 0.0D0
        DO K=KT,KB(IU)
          D(IU) = D(IU)+(U(K,IU)*BHR2(K,IU)-QSS(K,IU))+UXBR(K,IU)*DLT
        END DO
        IF (DN_FLOW(JB)) THEN
          DO K=KT,KB(ID)
            D(ID) = D(ID)-U(K,ID-1)*BHR2(K,ID-1)-QSS(K,ID)+(UXBR(K,ID)-UXBR(K,ID-1))*DLT+QOUT(K,JB)
          END DO
        END IF
        IF (UP_HEAD(JB)) THEN
          DO K=KT,KBMIN(IU-1)
            BHRHO(IU-1) = BHRHO(IU-1)+(BH2(K,IU)/RHO(K,IU)+BH2(K,IU-1)/RHO(K,IU-1))
          END DO
          DO K=KT,KB(IU)
            D(IU)   = D(IU)-U(K,IU-1)*BHR2(K,IU-1)
            F(IU-1) = F(IU-1)-(SB(K,IU-1)-ST(K,IU-1)+HDG(K,IU-1)-GRAV(K,IU-1))
          END DO
        END IF
        IF (DN_HEAD(JB)) THEN
          DO K=KT,KBMIN(ID)
            BHRHO(ID) = BHRHO(ID)+(BH2(K,ID+1)/RHO(K,ID+1)+BH2(K,ID)/RHO(K,ID))
          END DO
          DO K=KT,KB(ID)
            D(ID) = D(ID)+(U(K,ID)*BHR2(K,ID)-U(K,ID-1)*BHR2(K,ID-1)-QSS(K,ID))+(UXBR(K,ID)-UXBR(K,ID-1))*DLT
            F(ID) = F(ID)+(-SB(K,ID)+ST(K,ID)-HDG(K,ID)+GRAV(K,ID))
          END DO
        END IF
      END DO
    END DO
 !!$OMP END PARALLEL DO
    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
      IF(BR_INACTIVE(JB))CYCLE
        IU = CUS(JB)
        ID = DS(JB)
        IF (INTERNAL_FLOW(JB) .AND. .NOT. DAM_INFLOW(JB)) THEN                                                         !TC 08/03/04
          QIN(JB) = 0.0D0
          DO K=KTWB(JWUH(JB)),KB(UHS(JB))
            QIN(JB) = QIN(JB)+U(K,UHS(JB))*BHR2(K,UHS(JB))
          END DO
        END IF
        IF (UP_FLOW(JB)) D(IU) = D(IU)-QIN(JB)

!****** Boundary surface elevations

        IF (UH_INTERNAL(JB)) THEN
          Z(IU-1)    = ((-EL(KTWB(JWUH(JB)),UHS(JB))+Z(UHS(JB))*COSA(JBUH(JB)))+EL(KT,IU-1)+SINA(JB)*DLXR(IU-1))/COSA(JB)
          ELWS(IU-1) = EL(KT,IU-1)-Z(IU-1)*COSA(JB)
          KTI(IU-1)  = 2
          DO WHILE (EL(KTI(IU-1),IU-1) > ELWS(IU-1))
            KTI(IU-1) = KTI(IU-1)+1
          END DO
          KTI(IU-1) = MAX(KTI(IU-1)-1,2)
        END IF
        IF (UH_EXTERNAL(JB)) Z(IU-1) = (EL(KT,IU-1)-(ELUH(JB)+SINA(JB)*DLX(IU)*0.5D0))/COSA(JB)
        IF (DH_INTERNAL(JB)) THEN
          Z(ID+1)    = ((-EL(KTWB(JWDH(JB)),DHS(JB))+Z(DHS(JB))*COSA(JBDH(JB)))+EL(KT,ID+1))/COSA(JB)
          ELWS(ID+1) = EL(KT,ID+1)-Z(ID+1)*COSA(JB)
          KTI(ID+1)  = 2
          DO WHILE (EL(KTI(ID+1),ID+1) > ELWS(ID+1))
            KTI(ID+1) = KTI(ID+1)+1
          END DO
          KTI(ID+1) = MAX(KTI(ID+1)-1,2)
          IF (KTI(ID+1) >= KB(ID)) THEN
            Z(ID+1)    = Z(ID)-SLOPE(JB)*DLX(ID)/2.0D0
            ELWS(ID+1) = EL(KT,ID+1)-Z(ID+1)*COSA(JB)
            KTI(ID+1)  = 2
            DO WHILE (EL(KTI(ID+1),ID+1) > ELWS(ID+1))
              KTI(ID+1) = KTI(ID+1)+1
            END DO
            KTI(ID+1) = MAX(KTI(ID+1)-1,2)
          END IF
        END IF
        IF (DH_EXTERNAL(JB)) Z(ID+1) = (EL(KT,ID+1)-(ELDH(JB)-SINA(JB)*DLX(ID)*0.5D0))/COSA(JB)

!****** Implicit water surface elevation solution

        DO I=IU,ID   !CONCURRENT(I=IU:ID)   !FORALL(I=IU:ID)                       !DO I=IU,ID
          A(I) = -RHO(KT,I-1)*G*COSA(JB)*DLT*DLT* BHRHO(I-1)*0.5D0/DLXR(I-1)
          C(I) = -RHO(KT,I+1)*G*COSA(JB)*DLT*DLT* BHRHO(I)  *0.5D0/DLXR(I)
          V(I) =  RHO(KT,I)  *G*COSA(JB)*DLT*DLT*(BHRHO(I)  *0.5D0/DLXR(I)+BHRHO(I-1)*0.5D0/DLXR(I-1))+DLX(I)*BI(KT,I)
          D(I) =  DLT*(D(I)+DLT*(F(I)-F(I-1)))+DLX(I)*BI(KT,I)*Z(I)
        END DO                   
        IF (UP_HEAD(JB)) D(IU) = D(IU)-A(IU)*Z(IU-1)
        IF (DN_HEAD(JB)) D(ID) = D(ID)-C(ID)*Z(ID+1)
        BTA(IU) = V(IU)
        GMA(IU) = D(IU)
        DO I=IU+1,ID
          BTA(I) = V(I)-A(I)/BTA(I-1)*C(I-1)
          GMA(I) = D(I)-A(I)/BTA(I-1)*GMA(I-1)
        END DO
        Z(ID) = GMA(ID)/BTA(ID)
        DO I=ID-1,IU,-1
          Z(I) = (GMA(I)-C(I)*Z(I+1))/BTA(I)
        END DO

!****** Boundary water surface elevations

        IF (UP_FLOW(JB) .AND. .NOT. HEAD_FLOW(JB)) Z(IU-1) = Z(IU)
        IF (UP_FLOW(JB) .AND.       HEAD_FLOW(JB)) Z(IU-1) = (-EL(KTWB(JWUH(JB)),UHS(JB))+Z(UHS(JB))*COSA(JBUH(JB))+EL(KT,IU-1)  &
                                                             +SINA(JBUH(JB))*DLXR(IU-1))/COSA(JBUH(JB))
        IF (DN_FLOW(JB))                           Z(ID+1) = Z(ID)

!****** Updated surface layer and geometry

        IF (.NOT. TRAPEZOIDAL(JW)) THEN                                                                                !SW 07/16/04
          DO I=IU-1,ID+1
            IF (EL(KT,I)-Z(I)*COSA(JB) > EL(KTI(I),I)) THEN
              DO WHILE ( EL(KT,I)-Z(I)*COSA(JB) > EL(KTI(I),I) .AND. KTI(I) /= 2)
                Z(I)   = (EL(KT,I)-EL(KTI(I),I)-(EL(KT,I)-EL(KTI(I),I)-Z(I)*COSA(JB))*(B(KTI(I),I)/B(KTI(I)-1,I)))/COSA(JB)

                IF(MACROPHYTE_ON)THEN
                  KTIP=KTI(I)
!C  KEEPING TRACK IF COLUMN KTI HAS MACROPHYTES
                  IF(KTIP.GT.2)KTICOL(I)=.FALSE.
                END IF

                KTI(I) =  MAX(KTI(I)-1,2)
              END DO
            ELSE IF (EL(KT,I)-Z(I)*COSA(JB) < EL(KTI(I)+1,I)) THEN
              DO WHILE (EL(KT,I)-Z(I)*COSA(JB) < EL(KTI(I)+1,I) .AND. KTI(I) < KB(I))                   ! sw 7/18/11
                Z(I)   = (EL(KT,I)-EL(KTI(I)+1,I)-(EL(KT,I)-EL(KTI(I)+1,I)-Z(I)*COSA(JB))*(B(KTI(I),I)/B(KTI(I)+1,I)))/COSA(JB)
                KTI(I) =  KTI(I)+1
                IF(MACROPHYTE_ON)KTICOL(I)=.TRUE.  
                IF (KTI(I) >= KB(I)) EXIT
              END DO
            END IF
            BI(KT:KB(I),I) =  B(KT:KB(I),I)
            BI(KT,I)       =  B(KTI(I),I)
            H1(KT,I)       =  H(KT,JW)-Z(I)
            AVH1(KT,I)     = (H1(KT,I)+H1(KT+1,I))*0.5D0
            IF (KT == KTI(I) .OR. KTI(I) >= KB(I)) THEN
              BH1(KT,I) = B(KT,I)*H1(KT,I)
            ELSE
              BH1(KT,I) = BI(KT,I)*(EL(KT,I)-Z(I)*COSA(JB)-EL(KTI(I)+1,I))/COSA(JB)
            END IF
            DO K=KTI(I)+1,KT
              BH1(KT,I) = BH1(KT,I)+BNEW(K,I)*H(K,JW) !BNEW(K,I)*H(K,JW)   ! SW 1/23/06
            END DO
            BKT(I)    = BH1(KT,I)/H1(KT,I)
            IF(KBI(I) < KB(I))BKT(I)=BH1(KT,I)/(H1(KT,I)-(EL(KBI(I)+1,I)-EL(KB(I)+1,I)))    ! SW 1/23/06
            VOL(KT,I) = BH1(KT,I)*DLX(I)
          END DO
          DO I=IU-1,ID
            AVHR(KT,I) = H1(KT,I)  +(H1(KT,I+1) -H1(KT,I))/(0.5D0*(DLX(I)+DLX(I+1)))*0.5D0*DLX(I)                          !SW 07/29/04  (H1(KT,I+1) +H1(KT,I))*0.5   
            IF(KBI(I) < KB(I))AVHR(KT,I)=(H1(KT,I)-(EL(KBI(I)+1,I)-EL(KB(I)+1,I)))  &
               +(H1(KT,I+1)-(EL(KBI(I)+1,I+1)-EL(KB(I)+1,I+1)) -H1(KT,I)+(EL(KBI(I)+1,I)&
               -EL(KB(I)+1,I)))/(0.5D0*(DLX(I)+DLX(I+1)))*0.5D0*DLX(I)        ! SW 1/23/06
            BHR1(KT,I) =  BH1(KT,I)+(BH1(KT,I+1)-BH1(KT,I))/(0.5D0*(DLX(I)+DLX(I+1)))*0.5D0*DLX(I)                          !SW 07/29/04 (BH1(KT,I+1)+BH1(KT,I))*0.5 
            IF(CONSTRICTION(KT,I))THEN    ! SW 6/26/2018
              IF(BHR1(KT,I) > BCONSTRICTION(I)*H1(KT,I))BHR1(KT,I)= BCONSTRICTION(I)*H1(KT,I)
            ENDIF
          END DO
          AVHR(KT,ID+1) = H1(KT,ID+1)
          BHR1(KT,ID+1) = BH1(KT,ID+1)
          DLVOL(JB)        = 0.0
        ELSE                                                                                                           !SW 07/16/04
          DO I=IU-1,ID+1
            BI(KT:KB(I),I) =  B(KT:KB(I),I)
            CALL GRID_AREA2
            H1(KT,I)   =  H(KT,JW)-Z(I)
            AVH1(KT,I) = (H1(KT,I)+H1(KT+1,I))*0.5
            CALL GRID_AREA1 (EL(KT,I)-Z(I),EL(KT+1,I),BH1(KT,I),BI(KT,I))
            BKT(I)    = BH1(KT,I)/H1(KT,I)
            if(kbi(i) < kb(i))bkt(i)=bh1(kt,i)/(h1(kt,i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))    ! SW 1/23/06
            VOL(KT,I) = BH1(KT,I)*DLX(I)
          END DO
          DO I=IU-1,ID
            AVHR(KT,I) = H1(KT,I)  +(H1(KT,I+1) -H1(KT,I))/(0.5D0*(DLX(I)+DLX(I+1)))*0.5D0*DLX(I)                          !SW 07/29/04
            if(kbi(i) < kb(i))avhr(kt,i)=(h1(kt,i)-(el(kbi(i)+1,i)-el(kb(i)+1,i))) &
               +(H1(KT,I+1)-(el(kbi(i)+1,i+1)-el(kb(i)+1,i+1)) -H1(KT,I)+(el(kbi(i)+1,i)&
               -el(kb(i)+1,i)))/(0.5D0*(DLX(I)+DLX(I+1)))*0.5D0*DLX(I)                                                     ! SW 1/23/06
            BHR1(KT,I) = BH1(KT,I)+(BH1(KT,I+1)-BH1(KT,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                          !SW 07/29/04
          END DO
          AVHR(KT,ID+1) = H1(KT,ID+1)
          BHR1(KT,ID+1) = BH1(KT,ID+1)
          DLVOL(JB)     = 0.0
        END IF
        ELWS(CUS(JB):DS(JB)+1) = EL(KT,CUS(JB):DS(JB)+1)-Z(CUS(JB):DS(JB)+1)*COSA(JB)
        DO I=IU,ID
          DLVOL(JB) = DLVOL(JB)+(BH1(KT,I)-BH2(KT,I))*DLX(I)
          IF (KT == 2 .AND. H1(KT,I) > H(2,JW) .AND. .NOT. SURFACE_WARNING) THEN
            WRITE (WRN,'(A,I0,A,F0.3)') 'Water surface is above the top of layer 2 in segment ',I,' at day ',JDAY
            WARNING_OPEN    = .TRUE.
            SURFACE_WARNING = .TRUE.
          END IF
        END DO

        IF(MACROPHYTE_ON)THEN
!C  IF DEPTH IN KTI LAYER BECOMES GREATER THAN THRESHOLD, SETTING
!C      MACROPHYTE CONC. IN KTI COLUMN TO INITIAL CONC.
          DO I=IU,ID
            DEPKTI=ELWS(I)-EL(KTI(I)+1,I)

!******* MACROPHYTES, SETTING CONC. OF MACROPHYTES IN NEW COLUMNS TO
!********* INITIAL CONCENTRATION IF COLUMN DEPTH IS GREATER THAN 'THRKTI'
            IF(.NOT.KTICOL(I).AND.DEPKTI.GE.THRKTI)THEN
              KTICOL(I)=.TRUE.
              JT=KTI(I)
              MACT(JT,KT,I)=0.0
              DO M=1,NMC
                !MACRC(JT,KT,I,M)=MACWBCI(JW,M)
                IF (ISO_macrophyte(JW,m))  macrc(jt,kt,I,m) = macwbci(JW,m)     ! cb 3/7/16
                IF (VERT_macrophyte(JW,m)) macrc(jt,kt,I,m) = 0.1
                IF (long_macrophyte(JW,m)) macrc(jt,kt,I,m) = 0.1
                COLB=EL(KTI(I)+1,I)
                COLDEP=ELWS(I)-COLB
                !MACRM(JT,KT,I,M)=MACWBCI(JW,M)*COLDEP*CW(JT,I)*DLX(I)
                MACRM(JT,KT,I,M)=macrc(jt,kt,I,m)*COLDEP*CW(JT,I)*DLX(I)         ! cb 3/17/16                 
                MACT(JT,KT,I)=MACT(JT,KT,I)+MACWBCI(JW,M)
                MACMBRT(JB,M) = MACMBRT(JB,M)+MACRM(JT,KT,I,M)
              END DO
            END IF

!****** MACROPHYTES, WHEN COLUMN DEPTH IS LESS THAN 'THRKTI', ZEROING OUT CONC.
            IF(KTICOL(I).AND.DEPKTI.LT.THRKTI)THEN
              KTICOL(I)=.FALSE.
              JT=KTI(I)
              MACT(JT,KT,I)=0.0
              DO M=1,NMC
                MACMBRT(JB,M) = MACMBRT(JB,M)-MACRM(JT,KT,I,M)
                MACRC(JT,KT,I,M)=0.0
                MACRM(JT,KT,I,M)=0.0
              END DO
            END IF
          END DO
        END IF

!***********************************************************************************************************************************
!**                                             Task 2.2.4: Longitudinal velocities                                               **
!***********************************************************************************************************************************

        IUT = IU
        IDT = ID
        IF (UP_HEAD(JB)) IUT = IU-1
        IF (DN_HEAD(JB)) IDT = ID+1

!****** Pressures

        DO I=IUT,IDT
          DO K=KT,KB(I)
            P(K,I) = P(K-1,I)+RHO(K,I)*G*H1(K,I)*COSA(JB)
          END DO
        END DO

!****** Horizontal pressure gradients

        DO I=IUT,IDT-1
          HPG(KT,I) = DLXRHO(I)*(BKT(I)+BKT(I+1))*0.5D0*(H1(KT,I+1)*P(KT,I+1)-H1(KT,I)*P(KT,I))
          DO K=KT+1,KBMIN(I)
            HPG(K,I) = DLXRHO(I)*BHR2(K,I)*((P(K-1,I+1)-P(K-1,I))+(P(K,I+1)-P(K,I)))
          END DO
        END DO

!****** Boundary horizontal velocities

        IF (UP_FLOW(JB)) THEN
          IF (.NOT. HEAD_FLOW(JB)) THEN
            QINF(:,JB) = 0.0
            IF (PLACE_QIN(JW)) THEN

!************ Inflow layer

              K     = KT
              SSTOT = 0.0
              DO JC=NSSS,NSSE
                SSTOT = SSTOT+CIN(JC,JB)
              END DO
              RHOIN = DENSITY(TIN(JB),DMAX1(CIN(1,JB),0.0D0),DMAX1(SSTOT,0.0D0))
              DO WHILE (RHOIN > RHO(K,IU) .AND. K < KB(IU))
                K = K+1
              END DO
              KTQIN(JB) = K
              KBQIN(JB) = K

!************ Layer inflows

              VQIN  =  QIN(JB)*DLT
              VQINI =  VQIN
              QINFR =  1.0
              INCR  = -1
              DO WHILE (QINFR > 0.0D0)
                V1 = VOL(K,IU)
                IF (K <= KB(IU)) THEN
                  IF (VQIN > 0.5D0*V1) THEN
                    QINF(K,JB) = 0.5D0*V1/VQINI
                    QINFR      = QINFR-QINF(K,JB)
                    VQIN       = VQIN-QINF(K,JB)*VQINI
                    IF (K == KT) THEN
                      K    = KBQIN(JB)
                      INCR = 1
                    END IF
                  ELSE
                    QINF(K,JB) = QINFR
                    QINFR      = 0.0D0
                  END IF
                  IF (INCR < 0) KTQIN(JB) = K
                  IF (INCR > 0) KBQIN(JB) = MIN(KB(IU),K)
                  K = K+INCR
                ELSE
                  QINF(KT,JB) = QINF(KT,JB)+QINFR
                  QINFR       = 0.0D0
                END IF
              END DO
            ELSE
              KTQIN(JB) = KT
              KBQIN(JB) = KB(IU)
              BHSUM     = 0.0D0
              DO K=KT,KB(IU)
                BHSUM = BHSUM+BH1(K,IU)
              END DO
              DO K=KT,KB(IU)
                QINF(K,JB) = BH1(K,IU)/BHSUM
              END DO
            END IF
            DO K=KT,KB(IU)
              U(K,IU-1) = QINF(K,JB)*QIN(JB)/BHR1(K,IU-1)
            END DO
          ELSE
            KTQIN(JB) = KT
            KBQIN(JB) = KB(IU)
            IF (JBUH(JB) <= BE(JW) .AND. JBUH(JB) >= BS(JW)) THEN
              DO K=KT,KB(IU)
                U(K,IU-1) = U(K,UHS(JB))*BHR1(K,UHS(JB))/BHR1(K,IU-1)
              END DO
            ELSE
              CALL UPSTREAM_VELOCITY
            END IF
          END IF
        END IF
        IF (DN_FLOW(JB)) THEN
          DO K=KT,KB(ID)
            U(K,ID) = QOUT(K,JB)/BHR1(K,ID)
          END DO
        END IF
        IF (UP_HEAD(JB)) THEN
          DO K=KT,KB(IU-1)
            U(K,IU-1) = (BHR2(K,IU-1)*U(K,IU-1)+DLT*(-SB(K,IU-1)+ST(K,IU-1)-HPG(K,IU-1)+GRAV(K,IU-1)))/BHR1(K,IU-1)
          END DO
        END IF
        IF (DN_HEAD(JB)) THEN
          DO K=KT,KB(ID+1)
            U(K,ID) = (BHR2(K,ID)*U(K,ID)+DLT*(-SB(K,ID)+ST(K,ID)-HPG(K,ID)+GRAV(K,ID)))/BHR1(K,ID)
          END DO
        END IF

!****** Horizontal velocities

        DO I=IU,ID-1
          DO K=KT,KBMIN(I)
            U(K,I) = (BHR2(K,I)*U(K,I))/BHR1(K,I)+(DLT*(-SB(K,I)+ST(K,I)-ADMZ(K,I)+ADMZ(K-1,I)-ADMX(K,I)+DM(K,I)-HPG(K,I)+GRAV(K,I)&
                     +UXBR(K,I)/H2(K,I)))/BHR1(K,I)
            IF (INTERNAL_WEIR(K,I)) U(K,I) = 0.0D0
          END DO
        END DO

!****** Implicit vertical eddy viscosity

        IF (IMPLICIT_VISC(JW)) THEN
        !  AT = 0.0D0; CT = 0.0D0; VT = 0.0D0; DT = 0.0D0
        DO I=IUT,IDT-1                ! SW CODE SPEEDUP
            DO K=KT,KBMIN(I) 
            AT(K,I) = 0.0D0; CT(K,I) = 0.0D0; VT(K,I) = 0.0D0; DT(K,I) = 0.0D0
            ENDDO
        ENDDO
          DO I=IUT,IDT-1
            DO K=KT,KBMIN(I)            !KB(I)  SW 10/7/07
              AT(K,I) = -DLT/BHR1(K,I)*AZ(K-1,I)*(BHR1(K-1,I)/AVHR(K-1,I)+BR(K,I))  /(AVH1(K-1,I)+AVH1(K-1,I+1))
              CT(K,I) = -DLT/BHR1(K,I)*AZ(K,I)  *(BHR1(K,I)  /AVHR(K,I)  +BR(K+1,I))/(AVH1(K,I)  +AVH1(K,I+1))
              VT(K,I) =  1.0D0-AT(K,I)-CT(K,I)
              DT(K,I) =  U(K,I)
            END DO
            CALL TRIDIAG(AT(:,I),VT(:,I),CT(:,I),DT(:,I),KT,KBMIN(I),KMX,U(:,I))
          END DO
        END IF

!****** Corrected horizontal velocities

        IF (UP_HEAD(JB)) THEN
          IS    =  ID
          IE    =  IU-1
          INCR  = -1
          Q(IS) =  0.0D0
          DO K=KT,KB(ID)
            Q(IS) = Q(IS)+U(K,IS)*BHR1(K,IS)
          END DO
          QSSUM(IS) = 0.0D0
          DO K=KT,KB(IS)
            QSSUM(IS) = QSSUM(IS)+QSS(K,IS)
          END DO
        ELSE
          IS   = IU-1
          IE   = ID
          INCR = 1
          IF (DN_FLOW(JB)) IE = ID-1
          Q(IS) = 0.0D0
          DO K=KT,KB(IU)
            Q(IS) = Q(IS)+U(K,IS)*BHR1(K,IS)
          END DO
        END IF
        QC(IS) = Q(IS)
        DO I=IS+INCR,IE,INCR
          QSSUM(I) = 0.0D0
          DO K=KT,KB(I)
            QSSUM(I) = QSSUM(I)+QSS(K,I)
          END DO
          BHRSUM = 0.0D0
          Q(I)   = 0.0D0
          DO K=KT,KBMIN(I)
            IF (.NOT. INTERNAL_WEIR(K,I)) THEN
              BHRSUM = BHRSUM+BHR1(K,I)
              Q(I)   = Q(I)+U(K,I)*BHR1(K,I)
            END IF
          END DO
          IF (UP_HEAD(JB)) THEN
            QC(I) = QC(I+1)+(BH1(KT,I+1)-BH2(KT,I+1))*DLX(I+1)/DLT-QSSUM(I+1)
          ELSE
            QC(I) = QC(I-1)-(BH1(KT,I)  -BH2(KT,I))  *DLX(I)  /DLT+QSSUM(I)
          END IF
          DO K=KT,KBMIN(I)
            IF (INTERNAL_WEIR(K,I)) THEN
              U(K,I) = 0.0D0
            ELSE
              U(K,I) =  U(K,I)+(QC(I)-Q(I))/BHRSUM
              IF (Q(I) /= 0.0) QERR(I) = (Q(I)-QC(I))/Q(I)*100.0
            END IF
          END DO
        END DO

!****** Head boundary flows

        IF (UP_HEAD(JB)) QUH1(KT:KB(IU-1),JB) = U(KT:KB(IU-1),IU-1)*BHR1(KT:KB(IU-1),IU-1)
        IF (DN_HEAD(JB)) QDH1(KT:KB(ID+1),JB) = U(KT:KB(ID+1),ID)  *BHR1(KT:KB(ID+1),ID)

!***********************************************************************************************************************************
!**                                              Task 2.2.5: Vertical velocities                                                  **
!***********************************************************************************************************************************

        DO I=IU,ID
          DO K=KB(I)-1,KT,-1
            WT1    =  W(K+1,I)*BB(K+1,I)
            WT2    = (BHR(K+1,I)*U(K+1,I)-BHR(K+1,I-1)*U(K+1,I-1)-QSS(K+1,I))/DLX(I)
            W(K,I) = (WT1+WT2)/BB(K,I)
          END DO
        END DO
      END DO
    END DO

!***********************************************************************************************************************************
!**                                                  Task 2.2.6: Autostepping                                                     **
!***********************************************************************************************************************************

    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        DO I=CUS(JB),DS(JB)

        IF(DLTADD(JW)=='      ON'.and.ABS(H1(KT,I)-H2(KT,I))/H2(KT,I) > 0.35)THEN
            WRITE (WRN,'(A,F0.3,A,I0,A,F0.3/3(A,F0.3),a,i10,a)') 'Computational warning |h1-h2|/h2>0.35 on Julian day = ',JDAY,' at segment ',I,' timestep DLT= ',DLT,&
                                                  '   Water surface deviation [Z,m] = ',Z(I),' H1 layer thickness(m) = ',H1(KT,I),' H2 layer thickenss(m)=',h2(kt,i),' Iteration[NIT]=',nit,' DLT reduced'
            WARNING_OPEN = .TRUE.
            IF(H1(KT,I)>0.0)THEN
            CURMAX=DLT*0.5
            ELSE
            CURMAX=DLTMIN
            ENDIF
            GO TO 220

          ELSEIF (H1(KT,I) < 0.0) THEN
            WRITE (WRN,'(A,F0.3,A,I0/4(A,F0.3))') 'Computational warning at Julian day = ',JDAY,' at segment ',I,'timestep = ',DLT,&
                                                  ' water surface deviation [Z] = ',Z(I),' m  layer thickness = ',H1(KT,I),' m'
            WARNING_OPEN = .TRUE.
            IF (DLT > DLTMIN) THEN
              WRITE (WRN,'(A,I0/2(A,F0.3),A,I0)') 'Negative surface layer thickness in segment ',I,'  time step reduced to ',  &
                                                   DLTMIN,' s on day ',JDAY,' at iteration ',NIT
              WARNING_OPEN = .TRUE.
              CURMAX       =  DLTMIN
              GO TO 220
            ELSE
              WRITE (W2ERR,'(A,F0.3/A,I0)') 'Unstable water surface elevation on day ',JDAY,'negative surface layer thickness '//  &
                                            'using minimum timestep at iteration ',NIT
              WRITE(W2ERR,*)'Branch #:',jb,' in Waterbody:',jw,' Surface layer KT:',ktwb(jw)
              WRITE (W2ERR,'(A)') 'Segment, Surface layer thickness, m, Flow m3/s, U(KT,I) m/s, ELWS, m'
              DO II=MAX(CUS(JB),I-3),MIN(DS(JB),I+3)
                WRITE (W2ERR,'(T4,I3,T19,F10.2,t37,f10.2,1x,f10.2,2x,f10.2)') II,H1(KT,II),QC(II),U(KT,II),ELWS(II)                           ! SW 7/13/10
              END DO
              TEXT = 'Runtime error - see w2.err'
              ERROR_OPEN = .TRUE.
              GO TO 230
            END IF
          END IF
        END DO
        DO I=CUS(JB),DS(JB)
           IF (VISCOSITY_LIMIT(JW))THEN
              IF(AX(JW) >= 0.0)TAU1   = 2.0*AX(JW)/(DLX(I)*DLX(I))
           ENDIF   
          IF (CELERITY_LIMIT(JW))  CELRTY = SQRT((ABS(RHO(KB(I),I)-RHO(KT,I)))/1000.0*G*DEPTHB(KBI(I),I)*0.5)               ! SW 1/23/06
          DO K=KT,KB(I)
            IF (VISCOSITY_LIMIT(JW) .AND. .NOT. IMPLICIT_VISC(JW)) TAU2 = 2.0*AZ(K,I)/(H1(K,JW)*H1(K,JW))
            QTOT(K,I) = (ABS(U(K,I))*BHR1(K,I)+ABS(U(K,I-1))*BHR1(K,I-1)+(ABS(W(K,I))*BB(K,I)+ABS(W(K-1,I))*BB(K-1,I))*DLX(I)      &
                        +DLX(I)*ABS(BH2(K,I)-BH1(K,I))/DLT+ABS(QSS(K,I)))*0.5
              IF (VISCOSITY_LIMIT(JW).AND.AX(JW)<0.0)THEN
              TAU1   = 2.0*ABS(U(K,I))*ABS(AX(JW))*H(K,JW) /(DLX(I)*DLX(I))
              ENDIF  
            DLTCAL    = 1.0/((QTOT(K,I)/BH1(K,I)+CELRTY)/DLX(I)+TAU1+TAU2+NONZERO)
            IF (DLTCAL < CURMAX) THEN
              KLOC   = K
              ILOC   = I
              CURMAX = DLTCAL
              IF (DLTFF*CURMAX < MINDLT) THEN
                KMIN = K
                IMIN = I
              END IF
            END IF
          END DO
        END DO
      END DO
    END DO

!** Restore timestep dependent variables and restart calculations

220 CONTINUE
    IF (CURMAX < DLT .AND. DLT > DLTMIN) THEN
      DLT = DLTFF*CURMAX
      IF (DLT <= DLTMIN) THEN
        WRITE (WRN,'(A,F0.3/A,F0.3,A)') 'Computational warning at Julian day = ',JDAY,' timestep = ',DLT,' sec: DLT<DLTMIN set DLT=DLTMIN'
        WARNING_OPEN = .TRUE.
        DLT          =  DLTMIN
      END IF
      NV        = NV+1
      Z         = SZ
      U         = SU
      W         = SW
      AZ        = SAZ
      AVH2      = SAVH2
      AVHR      = SAVHR
      KTI       = SKTI
      BKT       = SBKT
      QSS       = 0.0
      !SP CEMA
      !if(sediment_diagenesis)then
      !  If(CEMARelatedCode .and. IncludeBedConsolidation)TSS       = 0.0  ! SW 7/27/2017
      !end if
      !End SP CEMA
      SB        = 0.0
      DLTS      = DLT

        do jw=1,nwb                                                                            ! SW 8/25/05
        do jb=bs(jw),be(jw)
        do i=us(jb)-1,ds(jb)+1
            VOL(KTWB(JW),I) = BH2(KTWB(JW),I)*DLX(I)
            BI(KTWB(JW),I) = B(KTI(I),I)
        end do
        end do
        end do


      CURMAX    = DLTMAXX/DLTFF
      IF (PIPES) THEN
        YS   = YSS
        VS   = VSS
        VST  = VSTS
        YST  = YSTS
        DTP  = DTPS
        QOLD = QOLDS
      END IF

!********** Macrophytes
      DO JW=1,NWB
        DO M=1,NMC
          IF (MACROPHYTE_CALC(JW,M)) THEN
            KT = KTWB(JW)
              DO JB=BS(JW),BE(JW)
                DO I=CUS(JB),DS(JB)
                  DO K=KT,KB(I)
                    MAC(K,I,M)=SMAC(K,I,M)
                    IF(KTICOL(I))THEN
                      JT=KTI(I)
                    ELSE
                      JT=KTI(I)+1
                    END IF
                    JE=KB(I)
                    DO J=JT,JE
                      MACRC(J,K,I,M)=SMACRC(J,K,I,M)
                      MACRM(J,K,I,M)=SMACRM(J,K,I,M)
                    END DO
                  END DO
                END DO
             END DO
          END IF
        END DO
      END DO

      GO TO 210
    END IF
    DLTLIM(KMIN,IMIN) = DLTLIM(KMIN,IMIN)+1.0

!** Layer bottom and middle depths

    DO JW=1,NWB
      DO JB=BS(JW),BE(JW)
        DO I=CUS(JB)-1,DS(JB)
          DEPTHB(KTWB(JW),I) = H1(KTWB(JW),I)
          DEPTHM(KTWB(JW),I) = H1(KTWB(JW),I)*0.5
             if(kbi(i) < kb(i)  .and. (el(kbi(i)+1,i)-el(kb(i)+1,i)) <  h1(ktwb(jw),i))then   ! SW 7/22/10 if h1 < elev diff this means depth is below the bottom - if we ignore that the run will continue but if dpethb is negative it will bomb in computing DECAY
             depthb(ktwb(jw),i)=(h1(ktwb(jw),i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))    ! SW 1/23/06
             depthm(ktwb(jw),i)=(h1(ktwb(jw),i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))*0.5   
             endif
          DO K=KTWB(JW)+1,KMX
            DEPTHB(K,I) = DEPTHB(K-1,I)+ H1(K,I)
            DEPTHM(K,I) = DEPTHM(K-1,I)+(H1(K-1,I)+H1(K,I))*0.5
          END DO
        END DO
      END DO
    END DO

! CHECK FOR DYNAMIC PIPE ADJUSTMENT SW 2/18/2020
      IF(NPI>0)THEN     ! SW 12/20/2020
        IF(DYNPAD=='ON'.AND.dynpipe(DYNPAD_PIPE) == '      ON')THEN    ! CHECK FOR WL VIOLATIONS
        IF(ELWS(DYNPAD_SEG)<DYNPAD_WL)THEN
            IF(DYNPAD_MAXRATE < (Z(DYNPAD_SEG)-SZ(DYNPAD_SEG))/DLT)THEN
                BP(DYNPAD_PIPE)=BP(DYNPAD_PIPE)*(1.+DYNPAD_PERCENTCHANGE/100.)
                 IF(BP(DYNPAD_PIPE)<=0.0)BP(DYNPAD_PIPE)=DYNPAD_PERCENTCHANGE/100.
                 IF(BP(DYNPAD_PIPE)>=1.0)BP(DYNPAD_PIPE)=1.0
                 WRITE(DYNPIPELOG,'(F9.3,",",4(F10.4,","),e12.4,",")')JDAY,BP(DYNPAD_PIPE),Z(DYNPAD_SEG),SZ(DYNPAD_SEG),DLT,(Z(DYNPAD_SEG)-SZ(DYNPAD_SEG))/DLT
            ENDIF
        ENDIF
        ENDIF
      ENDIF
! END DYN PIPE ADJUSTMENT
        
CALL temperature

IF (CONSTITUENTS) CALL wqconstituents

IF(FISH_PARTICLE_EXIST)CALL FISH ! SW 4/30/15

!if(sediment_diagenesis)then
  If(CEMARelatedCode .and. IncludeBedConsolidation)Call CEMAUpdateVerticalLayering
!end if

CALL LAYERADDSUB
if(error_open)go to 230

CALL BALANCES

CALL UPDATE

if(restart_in)then
  if(iopenfish==0)nxtmts=jday
endif

IF (JDAY.GE.NXTMTS.OR.JDAY.GE.TSRD(TSRDP+1).or.nit==1) THEN       ! OUTPUT AT FREQUENCY OF TSR FILES
IF(HABTATC  == '      ON')CALL FISHHABITAT(iopenfish)                                                  ! OUTPUT AT FREQUENCY OF TSR FILES
IF(AERATEC  == '      ON' .and. oxygen_demand)CALL AERATEOUTPUT
IF(RESTART_IN)IOPENFISH=1
IF(ENVIRPC  == '      ON')CALL ENVIRP
iopenfish=1
ENDIF                                                             ! OUTPUT AT FREQUENCY OF TSR FILES


CALL OUTPUTA
!**** Screen output
DO JW=1,NWB
      IF (SCREEN_OUTPUT(JW)) THEN
        IF (JDAY >= NXTMSC(JW) .OR. JDAY >= SCRD(SCRDP(JW)+1,JW)) THEN
          IF (JDAY >= SCRD(SCRDP(JW)+1,JW)) THEN
            SCRDP(JW)  = SCRDP(JW)+1
            NXTMSC(JW) = SCRD(SCRDP(JW),JW)
          END IF
          KT         = KTWB(JW)
          NXTMSC(JW) = NXTMSC(JW)+SCRF(SCRDP(JW),JW)
          CALL SCREEN_UPDATE (DLG)
          CALL DATE_AND_TIME (CDATE,CCTIME)
 !         DO JH=1,NHY
 !           IF (HYDRO_PLOT(JH))       CALL GRAPH_UPDATE (JH,HYD(:,:,JH),     HNAME(JH), HYMIN(JH),1.0,       LNAME(JH))
 !         END DO
 !         DO JC=1,NCT
 !           IF (CONSTITUENT_PLOT(JC)) CALL GRAPH_UPDATE (JH+JC,C2(:,:,JC),   CNAME(JC), CMIN(JC), CMULT(JC), LNAME(JH+JC))
 !         END DO
 !         DO JD=1,NDC
 !           IF (DERIVED_PLOT(JD))     CALL GRAPH_UPDATE (JH+JC+JD,CD(:,:,JD),CDNAME(JD),CDMIN(JD),CDMULT(JD),LNAME(JH+JC+JD))
 !         END DO
        END IF
      END IF
 END DO

END DO    ! END OF MAIN DO WHILE LOOP

230 CONTINUE
  IF (STOP_PUSHED) THEN
    TEXT  = 'Execution stopped at '//CCTIME(1:2)//':'//CCTIME(3:4)//':'//CCTIME(5:6)//' on '//CDATE(5:6)//'/'//CDATE(7:8)//'/'        &
                                   //CDATE(3:4)
    CALL RESTART_OUTPUT ('rso.opt')
  END IF

  IF(END_RUN .and. RESTART_OUT)CALL RESTART_OUTPUT ('rso.opt')  ! cb 4/9/15 writing restart output at end of simulation if RSOC='ON'
IOPENFISH=3
IF(ENVIRPC  == "      ON")CALL ENVIRP
IF(HABTATC  == "      ON")call fishhabitat(iopenfish)                                         ! FINAL OUTPUT FOR ENVIR PERFORMANCE
CALL ENDSIMULATION
  !IF (WAIT_FOR_INFLOW_RESULTS) THEN                                                                                     !SR 11/26/19
  !  DEALLOCATE (WAIT_TYPE, WAIT_INDEX, FILEDIR)                                                                         !SR 11/26/19
  !  CLOSE (9911)                                                                                                        !SR 11/26/19
  !END IF
! FISH OUTPUT SW 4/30/15  *************
        IF(FISH_PARTICLE_EXIST)call fishoutput  

240 CONTINUE
!  CALL DEALLOCATE_GRAPH
  
  IF(CLOSEC=='      ON' .AND. END_RUN)THEN
  CALL EXITDIALOG(DLG,TEXT)
  ELSE
  CALL STOP_W2 (DLG,TEXT)
  ENDIF
  RETURN
END FUNCTION CE_QUAL_W2
