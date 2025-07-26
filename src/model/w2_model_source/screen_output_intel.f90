!***********************************************************************************************************************************
!**                                                 F U N C T I O N   W I N M A I N                                               **
!***********************************************************************************************************************************

INTEGER*4 FUNCTION WINMAIN (HINSTANCE,HPREVINSTANCE,LPSZCMDLINE,NCMDSHOW)
  !DEC$ IF DEFINED(_X86_)
  !DEC$ ATTRIBUTES STDCALL, ALIAS : '_WinMain@16' :: WINMAIN
  !DEC$ ELSE
  !DEC$ ATTRIBUTES STDCALL, ALIAS: 'WinMain':: WinMain
  !DEC$ ENDIF
  USE DFLIB; USE DFWIN, RENAMED => DLT; USE DFLOGM
  INTEGER(4) :: HINSTANCE, HPREVINSTANCE, LPSZCMDLINE, NCMDSHOW

  CALL W2_DIALOG
  WINMAIN = 0
  RETURN
END FUNCTION WINMAIN

!***********************************************************************************************************************************
!**                                                P R O G R A M   W 2   D I A L O G                                              **
!***********************************************************************************************************************************

subroutine W2_DIALOG
  USE DFLOGM; USE MSCLIB
  INTEGER       :: RESULT
  LOGICAL       :: RSO_EXISTS=.FALSE., RESTARTED, RESULTLOG
  TYPE (DIALOG) :: DLG
  EXTERNAL RUN_W2

  RESTARTED = .FALSE.                                                                               
  OPEN (1,FILE='rso.opt',STATUS='OLD',IOSTAT=RESULT); CLOSE (1)                                   !Does restart file exist?
  RSO_EXISTS = RESULT == 0                                                                        !Restart file exists
  RESULTLOG = DLGINIT   (OUTPUT_DIALOG,DLG)                                                          !Initialize dialog box
  RESULTLOG = DLGSET    (DLG,STATUS,'Pending execution')                                             !Display execution status
  RESULTLOG = DLGSET    (DLG,RESTART,        RSO_EXISTS,DLG_ENABLE)                                  !Enable  'Restart' button
  RESULTLOG = DLGSET    (DLG,RUN,           .FALSE.,    DLG_ENABLE)
  RESULTLOG = DLGSET    (DLG,STOP_EXECUTION,.FALSE.,    DLG_ENABLE)                                  !Disable 'Stop'    button
  RESULTLOG = DLGSETSUB (DLG,OUTPUT_DIALOG,RUN_W2)
  RESULTLOG = DLGSETSUB (DLG,RUN,            RUN_W2)                                                 !Set 'Run'     callback
  RESULTLOG = DLGSETSUB (DLG,STOP_EXECUTION, RUN_W2)                                                 !Set 'Stop'    callback
  RESULTLOG = DLGSETSUB (DLG,RESTART,        RUN_W2)                                                 !Set 'Restart' callback
  RESULTLOG = DLGSETSUB (DLG,HIGHEST,        RUN_W2)                                                 !Set 'Highest' callback
  RESULTLOG = DLGSETSUB (DLG,HIGH,           RUN_W2)                                                 !Set 'High'    callback
  RESULTLOG = DLGSETSUB (DLG,NORMAL,         RUN_W2)                                                 !Set 'Normal'  callback
  RESULTLOG = DLGSETSUB (DLG,LOW,            RUN_W2)                                                 !Set 'Low'     callback
  RESULTLOG = DLGSETSUB (DLG,LOWEST,         RUN_W2)                                                 !Set 'Lowest'  callback
  RESULTLOG = DLGSETSUB (DLG,IDLE,           RUN_W2)                                                 !Set 'Idle'    callback
  RESULTLOG = DLGSETSUB (DLG,IDOK,RUN_W2)
  RESULTLOG = DLGMODAL  (DLG)                                                                        !Show dialog box
  CALL DLGUNINIT     (DLG)                                                                        !Close dialog box
END subroutine W2_DIALOG

!***********************************************************************************************************************************
!*                                                  S U B R O U T I N E   R U N   W 2                                             **
!***********************************************************************************************************************************

SUBROUTINE RUN_W2 (DLG,CONTROL_NAME,ACTION)
  USE DFWIN, RENAMED => DLT; USE DFLOGM; USE MSCLIB; USE GLOBAL, ONLY: CDATE, CCTIME;  USE MAIN, ONLY: END_RUN               !Rename DLT in DFWIN
  IMPLICIT NONE
  INTEGER                             :: ACTION, IDTHREAD, CONTROL_NAME, I
  LOGICAL                             :: RESULT
  CHARACTER(8)                        :: TIME
  CHARACTER(72)                       :: TEXT
  CHARACTER(1000)                     :: TEXT1
  TYPE(DIALOG)                        :: DLG
  TYPE(T_SECURITY_ATTRIBUTES),POINTER :: NULL_SA
  INTERFACE
    INTEGER(4) FUNCTION CE_QUAL_W2 (H)
    !DEC$ATTRIBUTES STDCALL :: CE_QUAL_W2
      INTEGER H
    END FUNCTION
  END INTERFACE

  SELECT CASE (CONTROL_NAME)
    CASE (RUN)
      CALL BLANK_DIALOG  (DLG)
      CALL DATE_AND_TIME (CDATE,CCTIME)
      RESTART_PUSHED = .FALSE.
      STOP_PUSHED    = .FALSE.
      HTHREAD        =  CREATETHREAD (NULL_SA,0,LOC(CE_QUAL_W2),LOC(DLG),0,LOC(IDTHREAD))          !Start W2 in a new thread
      TIME           =  CCTIME(1:2)//':'//CCTIME(3:4)//':'//CCTIME(5:6)
      RESULT         =  DLGSET (DLG,RUN,           .FALSE.,DLG_ENABLE)                             !Disable 'Run'     button
      RESULT         =  DLGSET (DLG,CLOSE,         .FALSE.,DLG_ENABLE)                             !Disable 'Run'     button
      RESULT         =  DLGSET (DLG,RESTART,       .FALSE.,DLG_ENABLE)                             !Disable 'Restart' button
      RESULT         =  DLGSET (DLG,STOP_EXECUTION,.TRUE., DLG_ENABLE)                             !Enable  'Stop'    button
      RESULT         =  DLGSET (DLG,STARTING_TIME,TIME)                                            !Display starting time
      RESULT         =  DLGSET (DLG,STATUS,'Executing')                                            !Display execution status
    CASE (OUTPUT_DIALOG)
      CALL BLANK_DIALOG  (DLG)
      CALL DATE_AND_TIME (CDATE,CCTIME)
      RESTART_PUSHED = .FALSE.
      STOP_PUSHED    = .FALSE.
      HTHREAD        =  CREATETHREAD (NULL_SA,0,LOC(CE_QUAL_W2),LOC(DLG),0,LOC(IDTHREAD))          !Start W2 in a new thread
      TIME           =  CCTIME(1:2)//':'//CCTIME(3:4)//':'//CCTIME(5:6)
      RESULT         =  DLGSET (DLG,RUN,           .FALSE.,DLG_ENABLE)                             !Disable 'Run'     button
      RESULT         =  DLGSET (DLG,CLOSE,         .FALSE.,DLG_ENABLE)                             !Disable 'Run'     button
      RESULT         =  DLGSET (DLG,RESTART,       .FALSE.,DLG_ENABLE)                             !Disable 'Restart' button
      RESULT         =  DLGSET (DLG,STOP_EXECUTION,.TRUE., DLG_ENABLE)                             !Enable  'Stop'    button
      RESULT         =  DLGSET (DLG,STARTING_TIME,TIME)                                            !Display starting time
      RESULT         =  DLGSET (DLG,STATUS,'Executing')                                            !Display execution status
    CASE (STOP_EXECUTION)
      STOP_PUSHED = .TRUE.
      END_RUN     = .TRUE.
      RESULT = DLGSET      (DLG,STOP_EXECUTION,.FALSE.,DLG_ENABLE)                                 !Disable 'Stop'    button
      RESULT = DLGSET      (DLG,RUN,           .TRUE., DLG_ENABLE)                                 !Enable  'Run'     button
      RESULT = CLOSEHANDLE (HTHREAD)                                                               !Close thread handle
    CASE (RESTART)
      STOP_PUSHED    = .FALSE.
      RESTART_PUSHED = .TRUE.
      END_RUN        = .FALSE.
      HTHREAD        =  CREATETHREAD (NULL_SA,0,LOC(CE_QUAL_W2),LOC(DLG),0,LOC(IDTHREAD))          !Start W2 in a new thread
      RESULT         =  DLGSET       (DLG,ENDING_TIME,   ' ')
      RESULT         =  DLGSET       (DLG,RUN,           .FALSE.,DLG_ENABLE)                       !Disable 'Run'     button
      RESULT         =  DLGSET       (DLG,CLOSE,         .FALSE.,DLG_ENABLE)                       !Disable 'Run'     button
      RESULT         =  DLGSET       (DLG,RESTART,       .FALSE.,DLG_ENABLE)                       !Disable 'Restart' button
      RESULT         =  DLGSET       (DLG,STOP_EXECUTION,.TRUE., DLG_ENABLE)                       !Enable  'Stop'    button
      RESULT         =  DLGSET       (DLG,STATUS,'Executing')                                      !Display execution status
    CASE (HIGHEST)
      CALL ENABLE (DLG)                                                                            !Enable  priority buttons
      I      = SETTHREADPRIORITY (HTHREAD,THREAD_PRIORITY_HIGHEST)                                 !Set highest priority
      RESULT = DLGSET            (DLG,HIGHEST,.FALSE.,DLG_ENABLE)                                  !Disable 'Highest' button
    CASE (HIGH)
      CALL ENABLE (DLG)                                                                            !Enable  priority buttons
      I      = SETTHREADPRIORITY (HTHREAD,THREAD_PRIORITY_ABOVE_NORMAL)                            !Set high    priority
      RESULT = DLGSET            (DLG,HIGH,   .FALSE.,DLG_ENABLE)                                  !Disable 'High'   button
    CASE (NORMAL)
      CALL ENABLE (DLG)                                                                            !Enable  priority buttons
      I      = SETTHREADPRIORITY (HTHREAD,THREAD_PRIORITY_NORMAL)                                  !Set normal  priority
      RESULT = DLGSET            (DLG,NORMAL, .FALSE.,DLG_ENABLE)                                  !Disable 'Normal' button
    CASE (LOW)
      CALL ENABLE (DLG)                                                                            !Enable  priority buttons
      I      = SETTHREADPRIORITY (HTHREAD,THREAD_PRIORITY_BELOW_NORMAL)                            !Set low     priority
      RESULT = DLGSET            (DLG,LOW,    .FALSE.,DLG_ENABLE)                                  !Disable 'Low'    button
    CASE (LOWEST)
      CALL ENABLE (DLG)                                                                            !Enable  priority buttons
      I      = SETTHREADPRIORITY (HTHREAD,THREAD_PRIORITY_LOWEST)                                  !Set lowest  priority
      RESULT = DLGSET            (DLG,LOWEST, .FALSE.,DLG_ENABLE)                                  !Disable 'Lowest' button
    CASE (IDLE)
      CALL ENABLE (DLG)                                                                            !Enable  priority buttons
      I      = SETTHREADPRIORITY (HTHREAD,THREAD_PRIORITY_IDLE)                                    !Set idle priority
      RESULT = DLGSET            (DLG,IDLE,   .FALSE.,DLG_ENABLE)                                  !Disable 'Idle'   button
  END SELECT
RETURN

ENTRY EXITDIALOG(DLG,TEXT)
 CALL DLGEXIT(DLG)
RETURN

ENTRY STOP_W2 (DLG,TEXT)
  CALL DATE_AND_TIME (CDATE,CCTIME)
  TEXT1  = CCTIME(1:2)//':'//CCTIME(3:4)//':'//CCTIME(5:6)
  RESULT = DLGSET (DLG,ENDING_TIME,TEXT1)                                                          !Display ending time
  RESULT = DLGSET (DLG,STATUS,TEXT)                                                                !Execution status
  RESULT = DLGSET (DLG,STOP_EXECUTION,.FALSE.,DLG_ENABLE)                                          !Disable 'Stop'    button
  RESULT = DLGSET (DLG,RUN,           .TRUE., DLG_ENABLE)                                          !Enable  'Run'     button
  RESULT = DLGSET (DLG,CLOSE,         .TRUE., DLG_ENABLE)                                          !Enable  'Close'    button
  IF (STOP_PUSHED) RESULT = DLGSET (DLG,RESTART,.TRUE.,DLG_ENABLE)                                 !Enable  'Restart' button
END SUBROUTINE RUN_W2

SUBROUTINE ENABLE (DLG)
  USE DFWIN, RENAMED => DLT; USE DFLOGM; USE MSCLIB                                                !Rename DLT in DFWIN
  LOGICAL :: RESULT                                                                                ! 9/27/07 SW
  TYPE(DIALOG) :: DLG

  RESULT = DLGSET (DLG,IDLE,   .TRUE.,DLG_ENABLE)                                                  !Enable 'Idle'    button
  RESULT = DLGSET (DLG,HIGHEST,.TRUE.,DLG_ENABLE)                                                  !Enable 'Highest' button
  RESULT = DLGSET (DLG,HIGH,   .TRUE.,DLG_ENABLE)                                                  !Enable 'High'    button
  RESULT = DLGSET (DLG,NORMAL, .TRUE.,DLG_ENABLE)                                                  !Enable 'Normal'  button
  RESULT = DLGSET (DLG,LOW,    .TRUE.,DLG_ENABLE)                                                  !Enable 'Low'     button
  RESULT = DLGSET (DLG,LOWEST, .TRUE.,DLG_ENABLE)                                                  !Enable 'Lowest'  button
END SUBROUTINE ENABLE

!***********************************************************************************************************************************
!*                                                E N T R Y   S C R E E N  U P D A T E                                            **
!***********************************************************************************************************************************

SUBROUTINE SCREEN_UPDATE (DLG)
  USE DFLOGM; USE DFLIB; USE MSCLIB; USE GEOMC; USE GLOBAL; USE GDAYC; USE SCREENC; USE SURFHE; USE TVDC; USE LOGICC; USE NAMESC
  USE STRUCTURES; USE MAIN, ONLY:TMSTRT, TMEND
  IMPLICIT NONE
  CHARACTER(8)    :: TIME
  CHARACTER(3000) :: TEXT1
  TYPE(DIALOG)    :: DLG
  INTEGER         :: IPROG
  LOGICAL         :: RESULT
  TIME = CCTIME(1:2)//':'//CCTIME(3:4)//':'//CCTIME(5:6)
  CALL DATE_AND_TIME (CDATE,CCTIME)
  CALL CPU_TIME (CURRENT)
  WRITE (TEXT1,'(A," ",I0,", ",I0)')  MONTH,GDAY,YEAR;            RESULT = DLGSET (DLG,GREGORIAN_DAY,                    TEXT1)
  WRITE (TEXT1,'(I0)')                INT(JDAY);                  RESULT = DLGSET (DLG,JULIAN_DAY,                       TEXT1)
  WRITE (TEXT1,'(F0.2)')             (JDAY-INT(JDAY))*24.0;       RESULT = DLGSET (DLG,JULIAN_HOUR,                      TEXT1)
  WRITE (TEXT1,'(I0)')                INT(ELTMJD);                RESULT = DLGSET (DLG,ELAPSED_DAY,                      TEXT1)
  WRITE (TEXT1,'(F0.2)')             (ELTMJD-INT(ELTMJD))*24.0;   RESULT = DLGSET (DLG,ELAPSED_HOUR,                     TEXT1)
  WRITE (TEXT1,'(I0)')                INT(DLTS1);                 RESULT = DLGSET (DLG,TIMESTEP,                         TEXT1)
  WRITE (TEXT1,'("(",I0,",",I0,")")') KLOC,ILOC;                  RESULT = DLGSET (DLG,TIMESTEP_LOCATION,                TEXT1)
  WRITE (TEXT1,'(I0)')                INT(MINDLT);                RESULT = DLGSET (DLG,MIN_TIMESTEP,                     TEXT1)
  WRITE (TEXT1,'("(",I0,",",I0,")")') KMIN,IMIN;                  RESULT = DLGSET (DLG,MIN_TIMESTEP_LOCATION,            TEXT1)
  WRITE (TEXT1,'(I0)')                INT(JDMIN);                 RESULT = DLGSET (DLG,MIN_TIMESTEP_DAY,                 TEXT1)
  WRITE (TEXT1,'(F0.2)')             (JDMIN-INT(JDMIN))*24.0;     RESULT = DLGSET (DLG,MIN_TIMESTEP_HOUR,                TEXT1)
  WRITE (TEXT1,'(I0)')                INT(DLTAV);                 RESULT = DLGSET (DLG,AVERAGE_TIMESTEP,                 TEXT1)
  WRITE (TEXT1,'(I0)')                NIT;                        RESULT = DLGSET (DLG,ITERATIONS,                       TEXT1)
  WRITE (TEXT1,'(I0)')                NV;                         RESULT = DLGSET (DLG,TIMESTEP_VIOLATIONS,              TEXT1)
  WRITE (TEXT1,'(F0.2)')              FLOAT(NV)/FLOAT(NIT)*100.0; RESULT = DLGSET (DLG,PERCENT,                          TEXT1)
  WRITE (TEXT1,'(F0.2)')              TAIR(JW);                   RESULT = DLGSET (DLG,AIR_TEMPERATURE,                  TEXT1)
  WRITE (TEXT1,'(F0.2)')              TDEW(JW);                   RESULT = DLGSET (DLG,DEW_POINT_TEMPERATURE,            TEXT1)
  WRITE (TEXT1,'(F0.2)')              WIND(JW);                   RESULT = DLGSET (DLG,WIND_SPEED,                       TEXT1)
  WRITE (TEXT1,'(F0.2)')              PHI(JW);                    RESULT = DLGSET (DLG,WIND_DIRECTION,                   TEXT1)
  WRITE (TEXT1,'(F0.2)')              CLOUD(JW);                  RESULT = DLGSET (DLG,CLOUD_COVER,                      TEXT1)
  WRITE (TEXT1,'(F0.1)')              ET(DS(JBDN(JW)));           RESULT = DLGSET (DLG,EQUILIBRIUM_TEMP,                 TEXT1)
  WRITE (TEXT1,'(F0.1)')              CSHE(DS(JBDN(JW)))*RHOWCP;  RESULT = DLGSET (DLG,SURFACE_HEAT_EXCHANGE,            TEXT1)
  WRITE (TEXT1,'(F0.1)')              SRON(JW);                   RESULT = DLGSET (DLG,SOLAR_RADIATION,                  TEXT1)
  WRITE (TEXT1,'(A)')                 TIME;                       RESULT = DLGSET (DLG,CURRENT_TIME,                     TEXT1)
  WRITE (TEXT1,'(I0)')                KTWB(JW);                   RESULT = DLGSET (DLG,SURFACE_LAYER,                    TEXT1)
  WRITE (TEXT1,'(F0.2)')              ELKT(JW);                   RESULT = DLGSET (DLG,SURFACE_ELEVATION,                TEXT1)
  WRITE (TEXT1,'(F0.2)')              ZMIN(JW);                   RESULT = DLGSET (DLG,MIN_DEVIATION,                    TEXT1)
  WRITE (TEXT1,'(I0)')                IZMIN(JW);                  RESULT = DLGSET (DLG,MIN_DEVIATION_SEGMENT,            TEXT1)
  WRITE (TEXT1,'(*(F0.2,2X))')        QIN;                        RESULT = DLGSET (DLG,BRANCH_INFLOW,                    TEXT1)
  WRITE (TEXT1,'(*(F0.2,2X))')        TIN;                        RESULT = DLGSET (DLG,BRANCH_INFLOW_TEMPERATURE,        TEXT1)
  WRITE (TEXT1,'(*(F0.2,2X))')        QDTR;                       RESULT = DLGSET (DLG,DIST_TRIBUTARY_INFLOW,            TEXT1)
  WRITE (TEXT1,'(*(F0.2,2X))')        TDTR;                       RESULT = DLGSET (DLG,DIST_TRIBUTARY_INFLOW_TEMPERATURE,TEXT1)
  WRITE (TEXT1,'(*(F0.3,2X))')        QTR;                        RESULT = DLGSET (DLG,TRIBUTARY_INFLOW,                 TEXT1)
  WRITE (TEXT1,'(*(F0.2,2X))')        TTR;                        RESULT = DLGSET (DLG,TRIBUTARY_INFLOW_TEMPERATURE,     TEXT1)
  WRITE (TEXT1,'(*(F0.2,2X))')        QSUM;                       RESULT = DLGSET (DLG,OUTFLOW,                          TEXT1)
  WRITE (TEXT1,'(*(F0.2,2X))')        QWD;                        RESULT = DLGSET (DLG,WITHDRAWAL,                       TEXT1)
  WRITE (TEXT1,'(*(F0.2,2X))')        QPU;                        RESULT = DLGSET (DLG,PUMPFLOW,                         TEXT1)
  WRITE (TEXT1,'(*(F0.2,2X))')        QSP;                        RESULT = DLGSET (DLG,SPILLWAYFLOW,                     TEXT1)
  WRITE (TEXT1,'(*(F0.2,2X))')        QPI;                        RESULT = DLGSET (DLG,PIPEFLOW,                         TEXT1)
  WRITE (TEXT1,'(*(F0.2,2X))')        QGT;                        RESULT = DLGSET (DLG,GATEFLOW,                         TEXT1)
  WRITE (TEXT1,'(A180)')           MODDIR;                        RESULT = DLGSET (DLG,MODELDIRECTORY,                 TEXT1) 

  WRITE (TEXT1,'(F0.2)')             (CURRENT)/60.0;              RESULT = DLGSET (DLG,CPU_TIMES,                        TEXT1)
  IF (MINDLT >= 1.0) WRITE (TEXT1,'(I0)') INT(MINDLT);            IF (MINDLT < 1.0) WRITE (TEXT1,'(F0.3)') MINDLT
  RESULT = DLGSET (DLG,MIN_TIMESTEP,                     TEXT1)
  
  IPROG=INT(((JDAY-TMSTRT)/(TMEND-TMSTRT))*100)   ! range is 0 to 100
  RESULT = DLGSET (DLG,PROGRESSBAR,IPROG,DLG_POSITION)
RETURN

ENTRY BLANK_DIALOG (DLG)
  RESULT = DLGSET (DLG,GREGORIAN_DAY,         ' ');               RESULT = DLGSET (DLG,JULIAN_DAY,                       ' ')
  RESULT = DLGSET (DLG,JULIAN_HOUR,           ' ');               RESULT = DLGSET (DLG,ELAPSED_DAY,                      ' ')
  RESULT = DLGSET (DLG,ELAPSED_HOUR,          ' ');               RESULT = DLGSET (DLG,TIMESTEP,                         ' ')
  RESULT = DLGSET (DLG,TIMESTEP_LOCATION,     ' ');               RESULT = DLGSET (DLG,MIN_TIMESTEP,                     ' ')
  RESULT = DLGSET (DLG,MIN_TIMESTEP_LOCATION, ' ');               RESULT = DLGSET (DLG,MIN_TIMESTEP_DAY,                 ' ')
  RESULT = DLGSET (DLG,MIN_TIMESTEP_HOUR,     ' ');               RESULT = DLGSET (DLG,AVERAGE_TIMESTEP,                 ' ')
  RESULT = DLGSET (DLG,ITERATIONS,            ' ');               RESULT = DLGSET (DLG,TIMESTEP_VIOLATIONS,              ' ')
  RESULT = DLGSET (DLG,PERCENT,               ' ');               RESULT = DLGSET (DLG,AIR_TEMPERATURE,                  ' ')
  RESULT = DLGSET (DLG,DEW_POINT_TEMPERATURE, ' ');               RESULT = DLGSET (DLG,WIND_SPEED,                       ' ')
  RESULT = DLGSET (DLG,WIND_DIRECTION,        ' ');               RESULT = DLGSET (DLG,CLOUD_COVER,                      ' ')
  RESULT = DLGSET (DLG,EQUILIBRIUM_TEMP,      ' ');               RESULT = DLGSET (DLG,SURFACE_HEAT_EXCHANGE,            ' ')
  RESULT = DLGSET (DLG,SOLAR_RADIATION,       ' ');               RESULT = DLGSET (DLG,CURRENT_TIME,                     ' ')
  RESULT = DLGSET (DLG,SURFACE_LAYER,         ' ');               RESULT = DLGSET (DLG,SURFACE_ELEVATION,                ' ')
  RESULT = DLGSET (DLG,MIN_DEVIATION,         ' ');               RESULT = DLGSET (DLG,MIN_DEVIATION_SEGMENT,            ' ')
  RESULT = DLGSET (DLG,BRANCH_INFLOW,         ' ');               RESULT = DLGSET (DLG,BRANCH_INFLOW_TEMPERATURE,        ' ')
  RESULT = DLGSET (DLG,DIST_TRIBUTARY_INFLOW, ' ');               RESULT = DLGSET (DLG,DIST_TRIBUTARY_INFLOW_TEMPERATURE,' ')
  RESULT = DLGSET (DLG,TRIBUTARY_INFLOW,      ' ');               RESULT = DLGSET (DLG,TRIBUTARY_INFLOW_TEMPERATURE,     ' ')
  RESULT = DLGSET (DLG,OUTFLOW,               ' ');               RESULT = DLGSET (DLG,ENDING_TIME,                      ' ')
  RESULT = DLGSET (DLG,WITHDRAWAL,            ' ');               RESULT = DLGSET (DLG,PIPEFLOW,                        ' ')
  RESULT = DLGSET (DLG,SPILLWAYFLOW,         ' ');                RESULT = DLGSET (DLG,GATEFLOW,                        ' ')
  RESULT = DLGSET (DLG,PUMPFLOW,             ' ');                RESULT = DLGSET (DLG,MODELDIRECTORY,             ' ')
  RESULT = DLGSET (DLG,PROGRESSBAR,0,DLG_POSITION)
END SUBROUTINE SCREEN_UPDATE

