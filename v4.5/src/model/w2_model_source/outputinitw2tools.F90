SUBROUTINE OUTPUTINIT  
    
  USE MAIN  
  USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY  
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART  
  USE MACROPHYTEC; USE POROSITYC; USE ZOOPLANKTONC;  USE BIOENERGETICS; USE CEMAVars, ONLY: SEDIMENT_DIAGENESIS
  USE ALGAE_TOXINS
  IMPLICIT NONE
  EXTERNAL RESTART_OUTPUT  
    
  REAL    DIST  
  REAL(4) :: ELC, DLX_OLD
  INTEGER*4,ALLOCATABLE,DIMENSION(:,:) :: ICOMP  
    
  INTEGER JN,IW, JO, IC, JJ 
  INTEGER*4 IFLAG  
  CHARACTER(60) :: TITLEWITH2  
  CHARACTER(100) :: TITLEWITH  
  CHARACTER(3) :: ICHAR3
  INTEGER, ALLOCATABLE, DIMENSION(:) :: IBR
    
  ALLOCATE(ICOMP(KMX,IMX))  
    
  !***********************************************************************************************************************************  
  !*                                                           Task 1.5: Outputs                                                    **  
  !***********************************************************************************************************************************  
    
      ALLOCATE(IBR(NBR))
      INQUIRE(FILE='w2_tecplotbr.csv',EXIST=BR_NOTECPLOT(1))     
    IF(BR_NOTECPLOT(1))THEN
        OPEN(CON,FILE='w2_tecplotbr.csv',STATUS='OLD')
        READ(CON,*)
        READ(CON,*)IC    ! NUMBER OF BRANCHES TO PLOT, MUST BE LESS THAN NBR
        READ(CON,*)
        READ(CON,*)(IBR(J),J=1,IC)
          DO J=1,IC
            BR_NOTECPLOT(IBR(J))=.FALSE.
          ENDDO
        CLOSE(CON)
    ELSE
        BR_NOTECPLOT=.FALSE.
    ENDIF
    DEALLOCATE(IBR)
  
  
  ! Open output files  
    
  IF (RESTART_IN) THEN  

      DO JW=1,NWB  ! INITIALIZE OUTPUT FOR TECPLOT FOR RESTART  9/4/2019
         IF(JW == 1)DIST=0.0  
      DO JB=BS(JW),BE(JW)  
        X1(US(JB))=DIST+DLX(US(JB))/2.  
        DO I=US(JB)+1,DS(JB)  
          DIST=DIST+(DLX(I)+DLX(I-1))/2.0  
          X1(I)=DIST  
        END DO  
        DIST=DIST+DLX(DS(JB))  
        X1(DS(JB)+1)=DIST  
      ENDDO  
    END DO        
      
    DO JW=1,NWB  
      IF (SNAPSHOT(JW))    OPEN (SNP(JW),FILE=SNPFN(JW),POSITION='APPEND')  
      !IF (VECTOR(JW))      OPEN (VPL(JW),FILE=VPLFN(JW),POSITION='APPEND')   *** DSI  
      IF (SPREADSHEET(JW)) OPEN (SPR(JW),FILE=SPRFN(JW),POSITION='APPEND')  
      IF(SPRC(JW) == '     ONV')OPEN (SPRV(JW),FILE=SPRVFN(JW),POSITION='APPEND')   ! SW 9/28/2018
      IF (CONTOUR(JW))     OPEN (CPL(JW),FILE=CPLFN(JW),POSITION='APPEND')  
      IF (PROFILE(JW))     OPEN (PRF(JW),FILE=PRFFN(JW),POSITION='APPEND')  
      IF (FLUX(JW))THEN  
        OPEN (FLX(JW),FILE=FLXFN(JW),POSITION='APPEND')  
        JDAY1=0.0  
        REWIND (FLX(JW))  
        DO WHILE (JDAY1<JDAY)  
          READ (FLX(JW),'(A72)',END=12) LINE  
                IF (LINE(1:8) == 'New date') THEN
            BACKSPACE(FLX(JW))  
            READ (FLX(JW),'(8X,F10.0)',END=12) JDAY1  
          ENDIF  
        ENDDO  
        DO J=1,14  
          BACKSPACE (FLX(JW))  
        END DO  
          
12      CONTINUE  
        WRITE (SEGNUM,'(I0)') JW  
        SEGNUM = ADJUSTL(SEGNUM)  
        L      = LEN_TRIM(SEGNUM)  
              OPEN (FLX2(JW), FILE='kflux_wb'//SEGNUM(1:L)//'.csv',POSITION='APPEND')
        JDAY1=0.0  
        REWIND (FLX2(JW))  
        READ   (FLX2(JW),'(/)',END=13)  
        DO WHILE (JDAY1 < JDAY)  
          READ (FLX2(JW),'(F10.0)',END=13) JDAY1  
        END DO  
        BACKSPACE (FLX2(JW))  
13      JDAY1=0.0  
      ENDIF  
      IF (SNAPSHOT(JW)) THEN  
        REWIND (SNP(JW))  
        DO WHILE (.TRUE.)  
          READ (SNP(JW),'(A72)',END=100) LINE  
          IF (LINE(26:28) == 'NIT') THEN  
            BACKSPACE SNP(JW)  
            READ (SNP(JW),'(31X,I10)',END=100) NIT1  
            IF (NIT1 > NIT) THEN  
              DO J=1,24  
                BACKSPACE (SNP(JW))  
              END DO  
              EXIT  
            END IF  
          END IF  
        END DO  
      END IF  
100   CONTINUE  
      IF (SPREADSHEET(JW)) THEN  
        REWIND (SPR(JW))  
        READ (SPR(JW),*)  
        DO WHILE (JDAY1 < JDAY)  
          READ (SPR(JW),*,END=101) LINE(1:38),JDAY1        !'(A,F10.0)'
        END DO  
        BACKSPACE (SPR(JW))  
101     CONTINUE  
        JDAY1 = 0.0  
        IF(SPRC(JW)=='     ONV')THEN
        REWIND (SPRV(JW))  
        READ (SPRV(JW),*)  
        DO WHILE (JDAY1 < JDAY)  
          READ (SPRV(JW),*,END=104) LINE(1:38),JDAY1        !'(A,F10.0)'
        END DO  
        BACKSPACE (SPRV(JW))  
104     CONTINUE  
        JDAY1 = 0.0   
        ENDIF
      END IF  
      IF (PROFILE(JW).and. iprf(1,1) /= -1) THEN    ! SW 4/1/2016
        REWIND (PRF(JW))  
        READ   (PRF(JW),'(A)')        (LINE,J=1,11)  
        READ   (PRF(JW),'(8I8)')       I  
        READ   (PRF(JW),'(10I8)')     (I,J=1,NIPRF(JW))  
        READ   (PRF(JW),'(20(1X,A))')  LINE (1:8), (LINE (1:3), JC=1,NCT),(LINE (1:3), JD=1,NDC)  
        READ   (PRF(JW),'(2A)')        LINE (1:26),(LINE (1:26),JC=1,NCT),(LINE (1:43),JD=1,NDC)  
        DO WHILE (JDAY1 < JDAY)  
          READ (PRF(JW),'(A72)',END=102) LINE  
          L1 = 0  
          L1 = SCAN(LINE,',')  
          IF (L1 /= 0) THEN  
            BACKSPACE (PRF(JW))  
            READ (PRF(JW),'(F8.0)',END=102) JDAY1  
          END IF  
        END DO  
        BACKSPACE (PRF(JW))  
        JDAY1 = 0.0  
      END IF  
102   CONTINUE  
      JDAY1 = 0.0  
      IF (CONTOUR(JW)) THEN  
        IF(TECPLOT(JW) /= '      ON')THEN 
        REWIND (CPL(JW))  
        DO WHILE (JDAY1 < JDAY)  
          READ (CPL(JW),'(A72)',END=103) LINE  
          IF (LINE(1:8) == 'New date') THEN
            BACKSPACE (CPL(JW))  
            READ (CPL(JW),'(A,F12.4)',END=103) LINE(1:9),JDAY1  
          END IF  
        END DO  
        BACKSPACE (CPL(JW))  
        JDAY1 = 0.0 
        ELSE
        REWIND (CPL(JW))  
        DO WHILE (JDAY1 < JDAY)  
          READ (CPL(JW),'(A72)',END=103) LINE  
          IF (LINE(1:8) == 'ZONE T="') THEN
            BACKSPACE (CPL(JW))  
            READ (CPL(JW),'(A,F9.0)',END=103) LINE(1:8),JDAY1  
          END IF  
        END DO  
        BACKSPACE (CPL(JW))  
        JDAY1 = 0.0                        
        ENDIF
      END IF  
103   CONTINUE  
    END DO  

    IF (VECTOR(1)) THEN  
    OPEN (VPL(1),FILE=VPLFN(1),STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='BINARY',POSITION='APPEND')  
    ENDIF
    
    IF (DOWNSTREAM_OUTFLOW) THEN  
      JFILE=0  
      L1 = SCAN(WDOFN,'.',BACK=.TRUE.)   ! SW 8/22/14 CHECKS FROM RIGHTHAND SIDE NOT LEFTHANDSIDE  
      DO JWD=1,NIWDO  
        WRITE (SEGNUM,'(I0)') IWDO(JWD)  
        SEGNUM = ADJUSTL(SEGNUM)  
        L      = LEN_TRIM(SEGNUM)  
        OPEN   (WDO(JWD,1),FILE='qwo_'//SEGNUM(1:L)//WDOFN(L1:L1+4),POSITION='APPEND')  ! '.opt' SW 4/14/2017
        REWIND (WDO(JWD,1))  
        READ   (WDO(JWD,1),'(//)',END=106)  
        DO WHILE (JDAY1 < JDAY)  
          READ (WDO(JWD,1),*,END=106) JDAY1                    !'(F8.0)'
        END DO  
        BACKSPACE (WDO(JWD,1))  
106     JDAY1 = 0.0  
        OPEN   (WDO(JWD,2),FILE='two_'//SEGNUM(1:L)//WDOFN(L1:L1+4),POSITION='APPEND')    ! '.opt' SW 4/14/2017 repeated for all .opt below for Downstream outflow
        REWIND (WDO(JWD,2))  
        READ   (WDO(JWD,2),'(//)',END=107)  
        DO WHILE (JDAY1 < JDAY)  
          READ (WDO(JWD,2),*,END=107) JDAY1                    !'(F8.0)'
        END DO  
        BACKSPACE (WDO(JWD,2))  
107     JDAY1=0.0  
        IF (CONSTITUENTS) THEN  
          OPEN   (WDO(JWD,3),FILE='cwo_'//SEGNUM(1:L)//WDOFN(L1:L1+4),POSITION='APPEND')
          REWIND (WDO(JWD,3))  
          READ   (WDO(JWD,3),'(//)',END=108)  
          DO WHILE (JDAY1 < JDAY)  
            READ (WDO(JWD,3),*,END=108) JDAY1                   ! '(F8.0)'
          END DO  
          BACKSPACE (WDO(JWD,3))  
108       CONTINUE  
          JDAY1=0.0  
        END IF  
        IF (DERIVED_CALC) THEN  
          OPEN   (WDO(JWD,4),FILE='dwo_'//SEGNUM(1:L)//WDOFN(L1:L1+4),POSITION='APPEND')
          REWIND (WDO(JWD,4))  
          READ   (WDO(JWD,4),'(//)')  
          DO WHILE (JDAY1 < JDAY)  
            READ (WDO(JWD,4),*,END=109) JDAY1                    ! '(F8.0)'
          END DO  
          BACKSPACE (WDO(JWD,4))  
109       CONTINUE  
          JDAY1 = 0.0  
        END IF  
          
        ! Determine the # of withdrawals at the WITH SEG  
        DO JB=1,NBR  ! structures
          IF(IWDO(JWD)==DS(JB) .AND. NSTR(JB) /= 0)THEN  
            DO JS=1,NSTR(JB)  
              JFILE=JFILE+1  
              WRITE(SEGNUM2,'(I0)')IWDO(JWD)  
              SEGNUM2 = ADJUSTL(SEGNUM2)  
              L2      = LEN_TRIM(SEGNUM2)  
              TITLEWITH2=' STR WITHDRAWAL AT SEG'//SEGNUM2(1:L2)  
              WRITE(SEGNUM,'(F10.0)')ESTR(JS,JB)  
              SEGNUM = ADJUSTL(SEGNUM)  
              L      = LEN_TRIM(SEGNUM)  
              TITLEWITH=ADJUSTL(TITLEWITH2)//' ELEV OF WITHDRAWAL CENTERLINE:'//SEGNUM(1:L)  
              WRITE(SEGNUM,'(I0)')JS  
              SEGNUM = ADJUSTL(SEGNUM)  
              L      = LEN_TRIM(SEGNUM)  
              WDO2(JFILE,1) = NUNIT; NUNIT = NUNIT+1  
            OPEN(WDO2(JFILE,1),FILE='qwo_str'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),POSITION='APPEND')
              REWIND(WDO2(JFILE,1))  
              READ   (WDO2(JFILE,1),'(//)')  
              DO WHILE (JDAY1 < JDAY)  
                READ (WDO2(JFILE,1),*,END=110) JDAY1           ! '(F8.0)'
              END DO  
              BACKSPACE (WDO2(JFILE,1))  
110           CONTINUE  
              JDAY1 = 0.0  
                
              WDO2(JFILE,2) = NUNIT; NUNIT = NUNIT+1  
            OPEN(WDO2(JFILE,2),FILE='two_str'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),POSITION='APPEND')
              REWIND(WDO2(JFILE,2))  
              READ   (WDO2(JFILE,2),'(//)')  
              DO WHILE (JDAY1 < JDAY)  
                READ (WDO2(JFILE,2),*,END=111) JDAY1          ! '(F8.0)'
              END DO  
              BACKSPACE (WDO2(JFILE,2))  
111           CONTINUE  
              JDAY1 = 0.0  
              IF (CONSTITUENTS) THEN  
                WDO2(JFILE,3) = NUNIT; NUNIT = NUNIT+1  
                OPEN(WDO2(JFILE,3),FILE='cwo_str'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),POSITION='APPEND')
                REWIND(WDO2(JFILE,3))  
                READ   (WDO2(JFILE,3),'(//)')  
                DO WHILE (JDAY1 < JDAY)  
                  READ (WDO2(JFILE,3),*,END=112) JDAY1        ! '(F8.0)'
                END DO  
                BACKSPACE (WDO2(JFILE,3))  
112             CONTINUE  
                JDAY1 = 0.0  
              ENDIF  
                
              IF (DERIVED_CALC) THEN  
                WDO2(JFILE,4) = NUNIT; NUNIT = NUNIT+1  
                OPEN(WDO2(JFILE,4),FILE='dwo_str'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),POSITION='APPEND')
                REWIND(WDO2(JFILE,4))  
                READ   (WDO2(JFILE,4),'(//)')  
                DO WHILE (JDAY1 < JDAY)  
                  READ (WDO2(JFILE,4),*,END=113) JDAY1          !'(F8.0)'
                END DO  
                BACKSPACE (WDO2(JFILE,4))  
113             CONTINUE  
                JDAY1 = 0.0  
              ENDIF  
            ENDDO  
          ENDIF  
        ENDDO  
          
        DO JS=1,NWD  ! withdrawals
          IF(IWDO(JWD) == IWD(JS))THEN  
            JFILE=JFILE+1  
            WRITE(SEGNUM2,'(I0)')IWDO(JWD)  
            SEGNUM2 = ADJUSTL(SEGNUM2)  
            L2      = LEN_TRIM(SEGNUM2)  
            TITLEWITH2=' WITHDRAWAL AT SEG'//SEGNUM2(1:L2)  
            WRITE(SEGNUM,'(F10.0)')EWD(JS)  
            SEGNUM = ADJUSTL(SEGNUM)  
            L      = LEN_TRIM(SEGNUM)  
            TITLEWITH=ADJUSTL(TITLEWITH2)//' ELEV OF WITHDRAWAL CENTERLINE:'//SEGNUM(1:L)  
            WRITE(SEGNUM,'(I0)')JS  
            SEGNUM = ADJUSTL(SEGNUM)  
            L      = LEN_TRIM(SEGNUM)  
            WDO2(JFILE,1) = NUNIT; NUNIT = NUNIT+1  
            OPEN(WDO2(JFILE,1),FILE='qwo_wd'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),POSITION='APPEND')
            REWIND(WDO2(JFILE,1))  
            READ   (WDO2(JFILE,1),'(//)')  
            DO WHILE (JDAY1 < JDAY)  
              READ (WDO2(JFILE,1),*,END=114) JDAY1           !'(F8.0)'
            END DO  
            BACKSPACE (WDO2(JFILE,1))  
114         CONTINUE  
            JDAY1 = 0.0  
            WDO2(JFILE,2) = NUNIT; NUNIT = NUNIT+1  
            OPEN(WDO2(JFILE,2),FILE='two_wd'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),POSITION='APPEND')
            REWIND(WDO2(JFILE,2))  
            READ   (WDO2(JFILE,2),'(//)')  
            DO WHILE (JDAY1 < JDAY)  
              READ (WDO2(JFILE,2),*,END=115) JDAY1       !'(F8.0)'
            END DO  
            BACKSPACE (WDO2(JFILE,2))  
115         CONTINUE  
            JDAY1 = 0.0  
            IF (CONSTITUENTS) THEN  
              WDO2(JFILE,3) = NUNIT; NUNIT = NUNIT+1  
                OPEN(WDO2(JFILE,3),FILE='cwo_wd'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),POSITION='APPEND')
              REWIND(WDO2(JFILE,3))  
              READ   (WDO2(JFILE,3),'(//)')  
              DO WHILE (JDAY1 < JDAY)  
                READ (WDO2(JFILE,3),*,END=116) JDAY1      !'(F8.0)'
              END DO  
              BACKSPACE (WDO2(JFILE,3))  
116           CONTINUE  
              JDAY1 = 0.0  
            ENDIF  
            IF (DERIVED_CALC) THEN  
              WDO2(JFILE,4) = NUNIT; NUNIT = NUNIT+1  
                OPEN(WDO2(JFILE,4),FILE='dwo_wd'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),POSITION='APPEND')
              REWIND(WDO2(JFILE,4))  
              READ   (WDO2(JFILE,4),'(//)')  
              DO WHILE (JDAY1 < JDAY)  
                READ (WDO2(JFILE,4),*,END=117) JDAY1        !'(F8.0)'
              END DO  
              BACKSPACE (WDO2(JFILE,4))  
117           CONTINUE  
              JDAY1 = 0.0  
            ENDIF  
          ENDIF  
        ENDDO  
          
        DO JS=1,NSP  ! spillways
          IF(IWDO(JWD) == IUSP(JS))THEN  
            JFILE=JFILE+1  
            WRITE(SEGNUM2,'(I0)')IWDO(JWD)  
            SEGNUM2 = ADJUSTL(SEGNUM2)  
            L2      = LEN_TRIM(SEGNUM2)  
            TITLEWITH2=' SPILLWAY WITHDRAWAL AT SEG'//SEGNUM2(1:L2)  
            WRITE(SEGNUM,'(F10.0)')ESP(JS)  
            SEGNUM = ADJUSTL(SEGNUM)  
            L      = LEN_TRIM(SEGNUM)  
            TITLEWITH=ADJUSTL(TITLEWITH2)//' ELEV OF WITHDRAWAL CENTERLINE:'//SEGNUM(1:L)  
            WRITE(SEGNUM,'(I0)')JS  
            SEGNUM = ADJUSTL(SEGNUM)  
            L      = LEN_TRIM(SEGNUM)  
            WDO2(JFILE,1) = NUNIT; NUNIT = NUNIT+1  
            OPEN(WDO2(JFILE,1),FILE='qwo_sp'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),POSITION='APPEND')
            REWIND(WDO2(JFILE,1))  
            READ   (WDO2(JFILE,1),'(//)')  
            DO WHILE (JDAY1 < JDAY)  
              READ (WDO2(JFILE,1),*,END=118) JDAY1      !'(F8.0)'
            END DO  
            BACKSPACE (WDO2(JFILE,1))  
118         CONTINUE  
            JDAY1 = 0.0  
            WDO2(JFILE,2) = NUNIT; NUNIT = NUNIT+1  
            OPEN(WDO2(JFILE,2),FILE='two_sp'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),POSITION='APPEND')
            REWIND(WDO2(JFILE,2))  
            READ   (WDO2(JFILE,2),'(//)')  
            DO WHILE (JDAY1 < JDAY)  
              READ (WDO2(JFILE,2),*,END=119) JDAY1      !'(F8.0)'
            END DO  
            BACKSPACE (WDO2(JFILE,2))  
119         CONTINUE  
            JDAY1 = 0.0  
            IF (CONSTITUENTS) THEN  
              WDO2(JFILE,3) = NUNIT; NUNIT = NUNIT+1  
                OPEN(WDO2(JFILE,3),FILE='cwo_sp'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),POSITION='APPEND')
              REWIND(WDO2(JFILE,3))  
              READ   (WDO2(JFILE,3),'(//)')  
              DO WHILE (JDAY1 < JDAY)  
                READ (WDO2(JFILE,3),*,END=120) JDAY1      !'(F8.0)'
              END DO  
              BACKSPACE (WDO2(JFILE,3))  
120           CONTINUE  
              JDAY1 = 0.0  
            ENDIF  
            IF (DERIVED_CALC) THEN  
              WDO2(JFILE,4) = NUNIT; NUNIT = NUNIT+1  
                OPEN(WDO2(JFILE,4),FILE='dwo_sp'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),POSITION='APPEND')
              REWIND(WDO2(JFILE,4))  
              READ   (WDO2(JFILE,4),'(//)')  
              DO WHILE (JDAY1 < JDAY)  
                READ (WDO2(JFILE,4),*,END=121) JDAY1          !'(F8.0)'
              END DO  
              BACKSPACE (WDO2(JFILE,4))  
121           CONTINUE  
              JDAY1 = 0.0  
            ENDIF  
          ENDIF  
        ENDDO  
          
          
        DO JS=1,NPU  ! pumps
          IF(IWDO(JWD) == IUPU(JS))THEN  
            JFILE=JFILE+1  
            WRITE(SEGNUM2,'(I0)')IWDO(JWD)  
            SEGNUM2 = ADJUSTL(SEGNUM2)  
            L2      = LEN_TRIM(SEGNUM2)  
            TITLEWITH2='PUMP WITHDRAWAL AT SEG'//SEGNUM2(1:L2)  
            WRITE(SEGNUM,'(F10.0)')EPU(JS)  
            SEGNUM = ADJUSTL(SEGNUM)  
            L      = LEN_TRIM(SEGNUM)  
            TITLEWITH=ADJUSTL(TITLEWITH2)//' ELEV OF WITHDRAWAL CENTERLINE:'//SEGNUM(1:L)  
            WRITE(SEGNUM,'(I0)')JS  
            SEGNUM = ADJUSTL(SEGNUM)  
            L      = LEN_TRIM(SEGNUM)  
            WDO2(JFILE,1) = NUNIT; NUNIT = NUNIT+1  
            OPEN(WDO2(JFILE,1),FILE='qwo_pmp'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),POSITION='APPEND')
            REWIND(WDO2(JFILE,1))  
            READ   (WDO2(JFILE,1),'(//)')  
            DO WHILE (JDAY1 < JDAY)  
              READ (WDO2(JFILE,1),*,END=122) JDAY1        !'(F8.0)'
            END DO  
            BACKSPACE (WDO2(JFILE,1))  
122         CONTINUE  
            JDAY1 = 0.0  
            WDO2(JFILE,2) = NUNIT; NUNIT = NUNIT+1  
            OPEN(WDO2(JFILE,2),FILE='two_pmp'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),POSITION='APPEND')
            REWIND(WDO2(JFILE,2))  
            READ   (WDO2(JFILE,2),'(//)')  
            DO WHILE (JDAY1 < JDAY)  
              READ (WDO2(JFILE,2),*,END=123) JDAY1         !'(F8.0)'
            END DO  
            BACKSPACE (WDO2(JFILE,2))  
123         CONTINUE  
            JDAY1 = 0.0  
            IF (CONSTITUENTS) THEN  
              WDO2(JFILE,3) = NUNIT; NUNIT = NUNIT+1  
                OPEN(WDO2(JFILE,3),FILE='cwo_pmp'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),POSITION='APPEND')
              REWIND(WDO2(JFILE,3))  
              READ   (WDO2(JFILE,3),'(//)')  
              DO WHILE (JDAY1 < JDAY)  
                READ (WDO2(JFILE,3),*,END=124) JDAY1      !'(F8.0)'
              END DO  
              BACKSPACE (WDO2(JFILE,3))  
124           CONTINUE  
              JDAY1 = 0.0  
            ENDIF  
            IF (DERIVED_CALC) THEN  
              WDO2(JFILE,4) = NUNIT; NUNIT = NUNIT+1  
                OPEN(WDO2(JFILE,4),FILE='dwo_pmp'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),POSITION='APPEND')
              REWIND(WDO2(JFILE,4))  
              READ   (WDO2(JFILE,4),'(//)')  
              DO WHILE (JDAY1 < JDAY)  
                READ (WDO2(JFILE,4),*,END=125) JDAY1     
              END DO  
              BACKSPACE (WDO2(JFILE,4))  
125           CONTINUE  
              JDAY1 = 0.0  
            ENDIF  
          ENDIF  
        ENDDO  
          
          
        DO JS=1,NPI  ! pipes
          IF(IWDO(JWD) == IUPI(JS))THEN  
            JFILE=JFILE+1  
            WRITE(SEGNUM2,'(I0)')IWDO(JWD)  
            SEGNUM2 = ADJUSTL(SEGNUM2)  
            L2      = LEN_TRIM(SEGNUM2)  
            TITLEWITH2='PIPE WITHDRAWAL AT SEG'//SEGNUM2(1:L2)  
            WRITE(SEGNUM,'(F10.0)')EUPI(JS)  
            SEGNUM = ADJUSTL(SEGNUM)  
            L      = LEN_TRIM(SEGNUM)  
            TITLEWITH=ADJUSTL(TITLEWITH2)//' ELEV OF WITHDRAWAL CENTERLINE:'//SEGNUM(1:L)  
            WRITE(SEGNUM,'(I0)')JS  
            SEGNUM = ADJUSTL(SEGNUM)  
            L      = LEN_TRIM(SEGNUM)  
            WDO2(JFILE,1) = NUNIT; NUNIT = NUNIT+1  
            OPEN(WDO2(JFILE,1),FILE='qwo_pipe'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),POSITION='APPEND')
            REWIND(WDO2(JFILE,1))  
            READ   (WDO2(JFILE,1),'(//)')  
            DO WHILE (JDAY1 < JDAY)  
              READ (WDO2(JFILE,1),*,END=126) JDAY1  
            END DO  
            BACKSPACE (WDO2(JFILE,1))  
126         CONTINUE  
            JDAY1 = 0.0  
            WDO2(JFILE,2) = NUNIT; NUNIT = NUNIT+1  
            OPEN(WDO2(JFILE,2),FILE='two_pipe'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),POSITION='APPEND')
            REWIND(WDO2(JFILE,2))  
            READ   (WDO2(JFILE,2),'(//)')  
            DO WHILE (JDAY1 < JDAY)  
              READ (WDO2(JFILE,2),*,END=127) JDAY1  
            END DO  
            BACKSPACE (WDO2(JFILE,2))  
127         CONTINUE  
            JDAY1 = 0.0  
            IF (CONSTITUENTS) THEN  
              WDO2(JFILE,3) = NUNIT; NUNIT = NUNIT+1  
                OPEN(WDO2(JFILE,3),FILE='cwo_pipe'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),POSITION='APPEND')
              REWIND(WDO2(JFILE,3))  
              READ   (WDO2(JFILE,3),'(//)')  
              DO WHILE (JDAY1 < JDAY)  
                READ (WDO2(JFILE,3),*,END=128) JDAY1  
              END DO  
              BACKSPACE (WDO2(JFILE,3))  
128           CONTINUE  
              JDAY1 = 0.0  
            ENDIF  
            IF (DERIVED_CALC) THEN  
              WDO2(JFILE,4) = NUNIT; NUNIT = NUNIT+1  
                OPEN(WDO2(JFILE,4),FILE='dwo_pipe'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),POSITION='APPEND')
              REWIND(WDO2(JFILE,4))  
              READ   (WDO2(JFILE,4),'(//)')  
              DO WHILE (JDAY1 < JDAY)  
                READ (WDO2(JFILE,4),*,END=129) JDAY1  
              END DO  
              BACKSPACE (WDO2(JFILE,4))  
129           CONTINUE  
              JDAY1 = 0.0  
            ENDIF  
          ENDIF  
        ENDDO  
          
        DO JS=1,NGT  ! gates
          IF(IWDO(JWD) == IUGT(JS))THEN  
            JFILE=JFILE+1  
            WRITE(SEGNUM2,'(I0)')IWDO(JWD)  
            SEGNUM2 = ADJUSTL(SEGNUM2)  
            L2      = LEN_TRIM(SEGNUM2)  
            TITLEWITH2='GATE WITHDRAWAL AT SEG'//SEGNUM2(1:L2)  
            WRITE(SEGNUM,'(F10.0)')EGT(JS)  
            SEGNUM = ADJUSTL(SEGNUM)  
            L      = LEN_TRIM(SEGNUM)  
            TITLEWITH=ADJUSTL(TITLEWITH2)//' ELEV OF WITHDRAWAL CENTERLINE:'//SEGNUM(1:L)  
            WRITE(SEGNUM,'(I0)')JS  
            SEGNUM = ADJUSTL(SEGNUM)  
            L      = LEN_TRIM(SEGNUM)  
            WDO2(JFILE,1) = NUNIT; NUNIT = NUNIT+1  
            OPEN(WDO2(JFILE,1),FILE='qwo_gate'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),POSITION='APPEND')
            REWIND(WDO2(JFILE,1))  
            READ   (WDO2(JFILE,1),'(//)')  
            DO WHILE (JDAY1 < JDAY)  
              READ (WDO2(JFILE,1),*,END=130) JDAY1  
            END DO  
            BACKSPACE (WDO2(JFILE,1))  
130         CONTINUE  
            JDAY1 = 0.0  
            WDO2(JFILE,2) = NUNIT; NUNIT = NUNIT+1  
            OPEN(WDO2(JFILE,2),FILE='two_gate'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),POSITION='APPEND')
            REWIND(WDO2(JFILE,2))  
            READ   (WDO2(JFILE,2),'(//)')  
            DO WHILE (JDAY1 < JDAY)  
              READ (WDO2(JFILE,2),*,END=131) JDAY1  
            END DO  
            BACKSPACE (WDO2(JFILE,2))  
131         CONTINUE  
            JDAY1 = 0.0  
            IF (CONSTITUENTS) THEN  
              WDO2(JFILE,3) = NUNIT; NUNIT = NUNIT+1  
                OPEN(WDO2(JFILE,3),FILE='cwo_gate'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),POSITION='APPEND')
              REWIND(WDO2(JFILE,3))  
              READ   (WDO2(JFILE,3),'(//)')  
              DO WHILE (JDAY1 < JDAY)  
                READ (WDO2(JFILE,3),*,END=132) JDAY1  
              END DO  
              BACKSPACE (WDO2(JFILE,3))  
132           CONTINUE  
              JDAY1 = 0.0  
            ENDIF  
            IF (DERIVED_CALC) THEN  
              WDO2(JFILE,4) = NUNIT; NUNIT = NUNIT+1  
                OPEN(WDO2(JFILE,4),FILE='dwo_gate'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),POSITION='APPEND')
              REWIND(WDO2(JFILE,4))  
              READ   (WDO2(JFILE,4),'(//)')  
              DO WHILE (JDAY1 < JDAY)  
                READ (WDO2(JFILE,4),*,END=133) JDAY1  
              END DO  
              BACKSPACE (WDO2(JFILE,4))  
133           CONTINUE  
              JDAY1 = 0.0  
            ENDIF  
          ENDIF  
        ENDDO  
          
      END DO  
    END IF  
    ! BIOENERGETICS mlm  
    IF (BIOEXP) THEN  
      JDAY1=0.0  
      L1 = SCAN(BIOFN,'.')  
      DO J=1,NIBIO  
        WRITE (SEGNUM,'(I0)') J  
        SEGNUM = ADJUSTL(SEGNUM)  
        L      = LEN_TRIM(SEGNUM)  
        BIOFN  = BIOFN(1:L1-1)//'_'//SEGNUM(1:L)//'.OPT'  
        OPEN  (BIOEXPFN(J),FILE=BIOFN,POSITION='APPEND')  
        REWIND (BIOEXPFN(J))  
        READ(bioexpfn(j),'(/)',end=1111)  
        DO WHILE (JDAY1<JDAY)  
          read (bioexpfn(j),'(F8.0)',end = 1111)JDAY1  
        END DO  
        BACKSPACE(BIOEXPFN(J))  
1111    CONTINUE  
        JDAY1=0.0  
      END DO  
        
      L1 = SCAN(WEIGHTFN,'.')  
      DO J=1,NIBIO  
        WRITE (SEGNUM,'(I0)') J  
        SEGNUM = ADJUSTL(SEGNUM)  
        L      = LEN_TRIM(SEGNUM)  
        WEIGHTFN  = WEIGHTFN(1:L1-1)//'_'//SEGNUM(1:L)//'.OPT'  
        OPEN  (WEIGHTNUM(J),FILE=WEIGHTFN,POSITION='APPEND')  
        REWIND (WEIGHTNUM(J))  
        READ(weightnum(j),'(/)',end=1112)  
        DO WHILE (JDAY1<JDAY)  
          read (weightnum(j),'(F10.0)',end = 1112)JDAY1  
        END DO  
        BACKSPACE(WEIGHTNUM(J))  
1112    CONTINUE  
        JDAY1=0.0  
      END DO  
    END IF  

    IF (TIME_SERIES) THEN  
      L1 = SCAN(TSRFN1,'.',BACK=.TRUE.)   ! SW 8/22/14 CHECKS FROM RIGHTHAND SIDE NOT LEFTHANDSIDE  
      DO J=1,NIKTSR  
        WRITE (SEGNUM,'(I0)') ITSR(J)  
        SEGNUM = ADJUSTL(SEGNUM)  
        L      = LEN_TRIM(SEGNUM)  
        WRITE (SEGNUM2,'(I0)')J  
        SEGNUM2 = ADJUSTL(SEGNUM2)  
        L2      = LEN_TRIM(SEGNUM2)  
        TSRFN  = TSRFN1(1:L1-1)//'_'//SEGNUM2(1:L2)//'_seg'//SEGNUM(1:L)//'.'//TSRFN1(L1+1:L1+4)                  !'.opt'
        OPEN   (TSR(J),FILE=TSRFN,POSITION='APPEND')  
        REWIND (TSR(J))  
       ! READ   (TSR(J),'(A72)',END=140)   (LINE,I=1,11)  
        READ   (TSR(J),'(/F10.3)',END=140) JDAYTS  
        DO WHILE (JDAYTS < JDAY)  
          READ (TSR(J),'(F10.0)',END=140) JDAYTS  
        END DO  
        BACKSPACE (TSR(J))  
140     CONTINUE  
      END DO  
    END IF  
    DO JW=1,NWB  
      IF (FLOWBALC=='      ON') THEN  
        OPEN(FLOWBFN,FILE='flowbal.csv',POSITION='APPEND')  
        JDAY1=0.0  
        REWIND (FLOWBFN)  
        READ   (FLOWBFN,'(/)',END=141)  
        DO WHILE (JDAY1 < JDAY)  
          READ (FLOWBFN,*,END=141) JDAY1  
        END DO  
        BACKSPACE (FLOWBFN)  
141     JDAY1=0.0  
        EXIT  
      ENDIF  
    ENDDO
    DO JW=1,NWB  
        IF (NPBALC=='      ON') THEN
        OPEN(MASSBFN,FILE='massbal.csv',POSITION='APPEND')  
        JDAY1=0.0  
        REWIND (MASSBFN)  
        READ   (MASSBFN,'(/)',END=143)  
        DO WHILE (JDAY1 < JDAY)  
          READ (MASSBFN,*,END=143) JDAY1  
        END DO  
        BACKSPACE (MASSBFN)  
143     JDAY1=0.0  
        EXIT  
      ENDIF  
    ENDDO  
    IF(WLC=='      ON')THEN  
      OPEN(WLFN,FILE='wl.csv',POSITION='APPEND')  
      JDAY1=0.0  
      REWIND (WLFN)  
      READ   (WLFN,'(//)',END=142)  
      DO WHILE (JDAY1 < JDAY)  
        READ (WLFN,'(F10.0)',END=142) JDAY1  
      END DO  
      BACKSPACE (WLFN)  
142   JDAY1=0.0  
    ENDIF  
  ELSE  
    ! *** Cold Start Output File Initialization  
      
    DO JW=1,NWB  ! moved so that x1 get initialized even if tecplot (and contour) OFF for first water body but ON for others   ! cb 10/10/10
      !c calculating longitudinal distance of segments  
      IF(JW == 1)DIST=0.0  
      DO JB=BS(JW),BE(JW)  
        X1(US(JB))=DIST+DLX(US(JB))/2.  
        DO I=US(JB)+1,DS(JB)  
          DIST=DIST+(DLX(I)+DLX(I-1))/2.0  
          X1(I)=DIST  
        END DO  
        DIST=DIST+DLX(DS(JB))  
        X1(DS(JB)+1)=DIST  
      ENDDO  
    END DO  
      
  INQUIRE(FILE='w2_lake_river_contour.csv',EXIST=LAKE_RIVER_CONTOURC)    
  
IF(LAKE_RIVER_CONTOURC)THEN
  OPEN (LAKE_RIVER_CONTOUR,FILE='w2_lake_river_contour.csv',STATUS='OLD',IOSTAT=I)
  READ(LAKE_RIVER_CONTOUR,*)
  READ(LAKE_RIVER_CONTOUR,*)LAKE_RIVER_CONTOUR_ON
  READ(LAKE_RIVER_CONTOUR,*)
  READ(LAKE_RIVER_CONTOUR,*)NUM_LAKE_CONTOUR,LAKE_CONTOUR_FORMAT
  READ(LAKE_RIVER_CONTOUR,*)
  DO JJ=1,NUM_LAKE_CONTOUR
  READ(LAKE_RIVER_CONTOUR,*)LAKE_CONTOUR_SEG(JJ),LAKE_CONTOUR_START(JJ),LAKE_CONTOUR_FREQ(JJ)
  IF(JDAY >= LAKE_CONTOUR_START(JJ))THEN
      NXT_LAKE_CONTOUR(JJ)=JDAY
      ELSE
      NXT_LAKE_CONTOUR(JJ)=LAKE_CONTOUR_START(JJ)
  ENDIF
  
    DO JW=1,NWB
      IF(US(BS(JW))<= LAKE_CONTOUR_SEG(JJ) .AND. DS(BE(JW)) >= LAKE_CONTOUR_SEG(JJ))THEN
      JW_LAKE_CONTOUR(JJ)=JW
      EXIT
      ENDIF
    ENDDO
    ENDDO
  READ(LAKE_RIVER_CONTOUR,*)
  READ(LAKE_RIVER_CONTOUR,*)NUM_RIVER_CONTOUR,RIVER_CONTOUR_FORMAT
  READ(LAKE_RIVER_CONTOUR,*)
  DO JJ=1,NUM_RIVER_CONTOUR
  READ(LAKE_RIVER_CONTOUR,*)RIVER_CONTOUR_BR1(JJ),RIVER_CONTOUR_BR2(JJ),RIVER_CONTOUR_START(JJ),RIVER_CONTOUR_FREQ(JJ)
  IF(JDAY >= RIVER_CONTOUR_START(JJ))THEN
      NXT_RIVER_CONTOUR(JJ)=JDAY
      ELSE
      NXT_RIVER_CONTOUR(JJ)=RIVER_CONTOUR_START(JJ)
  ENDIF
    DO JW=1,NWB
      IF(BS(JW)<= RIVER_CONTOUR_BR1(JJ) .AND. BE(JW) >= RIVER_CONTOUR_BR1(JJ))THEN
      JW_RIVER_CONTOUR(JJ)=JW
      EXIT
      ENDIF
    ENDDO
  ENDDO
  CLOSE(LAKE_RIVER_CONTOUR)
  IF(LAKE_RIVER_CONTOUR_ON=='ON')THEN
          DO JJ=1,NUM_LAKE_CONTOUR
          WRITE(ICHAR3,'(I3)')LAKE_CONTOUR_SEG(JJ)
          FILE_LAKE_CONTOUR_T(JJ)='LakeContour_T_Seg'//ADJUSTL(trim(ICHAR3))//'.csv'
          OPEN(LAKE_RIVER_CONTOUR+JJ-1,FILE=FILE_LAKE_CONTOUR_T(JJ),status='unknown')
          IF(LAKE_CONTOUR_FORMAT == 1)THEN
          WRITE(LAKE_RIVER_CONTOUR+JJ-1,*)'JDAY,ELEVATION(m),TEMPERATURE(C)'
          ELSE
          WRITE(LAKE_RIVER_CONTOUR+JJ-1,'("TIME,",*(F8.2,","))')((EL(K,LAKE_CONTOUR_SEG(JJ))+EL(K,LAKE_CONTOUR_SEG(JJ)))*0.5,K=2,KB(LAKE_CONTOUR_SEG(JJ)))
          ENDIF
              IF(OXYGEN_DEMAND)THEN
              WRITE(ICHAR3,'(I3)')LAKE_CONTOUR_SEG(JJ)
              FILE_LAKE_CONTOUR_DO(JJ)='LakeContour_DO_Seg'//ADJUSTL(trim(ICHAR3))//'.csv'
              OPEN(LAKE_RIVER_CONTOUR+10+JJ-1,FILE=FILE_LAKE_CONTOUR_DO(JJ),status='unknown')
              IF(LAKE_CONTOUR_FORMAT == 1)THEN
              WRITE(LAKE_RIVER_CONTOUR+10+JJ-1,*)'JDAY,ELEVATION(m),DisslvedOxygen(mg/l)'
              ELSE
              WRITE(LAKE_RIVER_CONTOUR+10+JJ-1,'("TIME,",*(F8.2,","))')((EL(K,LAKE_CONTOUR_SEG(JJ))+EL(K,LAKE_CONTOUR_SEG(JJ)))*0.5,K=2,KB(LAKE_CONTOUR_SEG(JJ)))
              ENDIF
              ENDIF
          ENDDO
          
          DO JJ=1,NUM_RIVER_CONTOUR
          WRITE(ICHAR3,'(I3)')RIVER_CONTOUR_BR1(JJ)
          FILE_RIVER_CONTOUR_T(JJ)='RiverContour_T_Br'//ADJUSTL(trim(ICHAR3))//'.csv'
          OPEN(LAKE_RIVER_CONTOUR+20+JJ-1,FILE=FILE_RIVER_CONTOUR_T(JJ),status='unknown')
          IF(RIVER_CONTOUR_FORMAT==1)THEN
          WRITE(LAKE_RIVER_CONTOUR+20+JJ-1,*)'JDAY,DISTANCE(M),TEMPERATURE(C)'
          ELSE
          WRITE(LAKE_RIVER_CONTOUR+20+JJ-1,'("TIME,",*(F8.2,","))')((X1(I),I=US(JB),DS(JB)),JB=RIVER_CONTOUR_BR1(JJ),RIVER_CONTOUR_BR2(JJ))
          ENDIF
              IF(OXYGEN_DEMAND)THEN
              WRITE(ICHAR3,'(I3)')LAKE_CONTOUR_SEG(JJ)
              FILE_RIVER_CONTOUR_DO(JJ)='RiverContour_DO_Br'//ADJUSTL(trim(ICHAR3))//'.csv'
              OPEN(LAKE_RIVER_CONTOUR+30+JJ-1,FILE=FILE_RIVER_CONTOUR_DO(JJ),status='unknown')
              IF(RIVER_CONTOUR_FORMAT == 1)THEN
              WRITE(LAKE_RIVER_CONTOUR+30+JJ-1,*)'JDAY,ELEVATION(M),DisslvedOxygen(mg/l)'
              ELSE
              WRITE(LAKE_RIVER_CONTOUR+30+JJ-1,'("TIME,",*(F8.2,","))')((EL(K,LAKE_CONTOUR_SEG(JJ))+EL(K,LAKE_CONTOUR_SEG(JJ)))*0.5,K=2,KB(LAKE_CONTOUR_SEG(JJ)))
              ENDIF
              ENDIF
          ENDDO
    ENDIF  
ENDIF

    DO JW=1,NWB  
      IF (SNAPSHOT(JW))    OPEN (SNP(JW),FILE=SNPFN(JW),STATUS='UNKNOWN')  
      !IF (VECTOR(JW))      OPEN (VPL(JW),FILE=VPLFN(JW),STATUS='UNKNOWN')  
      IF (PROFILE(JW))     OPEN (PRF(JW),FILE=PRFFN(JW),STATUS='UNKNOWN')  
      IF (SPREADSHEET(JW)) OPEN (SPR(JW),FILE=SPRFN(JW),STATUS='UNKNOWN')  
      IF(SPRC(JW)=='     ONV') OPEN (SPRV(JW),FILE=SPRVFN(JW),STATUS='UNKNOWN')     ! SW 9/28/2018
      IF (CONTOUR(JW))     OPEN (CPL(JW),FILE=CPLFN(JW),STATUS='UNKNOWN')  
      IF (FLUX(JW))THEN  
        OPEN (FLX(JW),FILE=FLXFN(JW),STATUS='UNKNOWN')  
        WRITE (SEGNUM,'(I0)') JW  
        SEGNUM = ADJUSTL(SEGNUM)  
        L      = LEN_TRIM(SEGNUM)  
              OPEN (FLX2(JW), FILE='kflux_wb'//SEGNUM(1:L)//'.csv',STATUS='UNKNOWN')    ! SW 7/1/2019
        WRITE(FLX2(JW), '("JDAY,  ELTM,",*(A,","))')(KFNAME2(KFCN(JF,JW)),JF=1,NAF(JW))  
      ENDIF  
        
      !**** Output files  
        
      IF (PROFILE(JW).and. iprf(1,1) /= -1) THEN    ! SW 4/1/2016
        TTIME = TMSTRT  
        DO WHILE (TTIME <= TMEND)  
          NDSP = NDSP+1  
          TTIME = TTIME+PRFF(PRFDP(JW),JW)  
          IF (TTIME >= PRFD(PRFDP(JW)+1,JW)) PRFDP(JW) = PRFDP(JW)+1  
        END DO  
        PRFDP(JW) = 1  
        WRITE (PRF(JW),'(A)')         TITLE  
        WRITE (PRF(JW),'(8I8,L2)')    KMX,NIPRF(JW),NDSP,NCT,NDC,NAC+NACD(JW)+1,PRFDP(JW),KTWB(JW),CONSTITUENTS  
        WRITE (PRF(JW),'(10I8)')      IPRF(1:NIPRF(JW),JW)  
        WRITE (PRF(JW),'(20(1X,A))') ' ON',CPRWBC(:,JW)(6:8),CDWBC(:,JW)(6:8)  
        WRITE (PRF(JW),'(2A)')       'Temperature,oC                            ',ADJUSTL(CNAME),ADJUSTL(CDNAME)  
        WRITE (PRF(JW),'(20I4)')      1,CN(1:NAC)+1,CDN(1:NACD(JW),JW)+NCT+1  
        WRITE (PRF(JW),'(10F8.0)')    1.0,CMULT,CDMULT  
        WRITE (PRF(JW),'(20I4)')      (KB(IPRF(I,JW)),I=1,NIPRF(JW))        ! KB(IPRF(1:NIPRF(JW),JW))  
        WRITE (PRF(JW),'(10F8.2)')    H  
        DO JP=1,NIPRF(JW)  
          NRS = KB(IPRF(JP,JW))-KTWB(JW)+1  
          WRITE (PRF(JW),'(A8,I4/(8(F10.2)))') 'TEMP    ',NRS,(T2(K,IPRF(JP,JW)),K=KTWB(JW),KB(IPRF(JP,JW)))  
        END DO  
        DO JC=1,NAC  
          IF (PRINT_CONST(CN(JC),JW)) THEN  
            DO JP=1,NIPRF(JW)  
              NRS = KB(IPRF(JP,JW))-KTWB(JW)+1  
              WRITE (PRF(JW),'(A,I4/(8(E13.6,1x)))') ADJUSTL(CNAME2(CN(JC))),NRS,(C2(K,IPRF(JP,JW),CN(JC))*CMULT(CN(JC)),K=KTWB(JW),  &                    ! CB 1/25/05  
              KB(IPRF(JP,JW)))  
            END DO  
          END IF  
        END DO  
        DO JD=1,NACD(JW)  
          DO JP=1,NIPRF(JW)  
            NRS = KB(IPRF(JP,JW))-KTWB(JW)+1  
            WRITE (PRF(JW),'(A,I4/(8(E13.6,1x)))') ADJUSTL(CDNAME2(CDN(JD,JW))),NRS,(CD(K,IPRF(JP,JW),CDN(JD,JW))*CDMULT(CDN(JD,JW)), &                     ! CB 1/25/05  
            K=KTWB(JW),KB(IPRF(JP,JW)))  
          END DO  
        END DO  
      END IF  
      IF (SPREADSHEET(JW)) THEN  
        DO J=1,NISPR(JW)  
          WRITE (SEGNUM,'(I0)') ISPR(J,JW)  
          SEGNUM = ADJUSTL(SEGNUM)  
          L      = LEN_TRIM(SEGNUM)  
          SEG(J) = 'Seg_'//SEGNUM(1:L)
        END DO  
        IF(SPRC(JW)   == '     ONV')THEN
            WRITE (SPRV(JW),'(A,A,*(A7,","))') 'Constituent,','Julian_day,',(SEG(J),J=1,NISPR(JW))     ! SW 9/28/2018
        ENDIF
            WRITE (SPR(JW),'(A,A,A,*("Elevation,",A7,","))') 'Constituent,','Julian_day,','Depth,',(SEG(J),J=1,NISPR(JW))  
        
      END IF  
      IF (CONTOUR(JW)) THEN  
        IF(TECPLOT(JW) /= '      ON')THEN  
          WRITE (CPL(JW),'(A)')           TITLE  
          WRITE (CPL(JW),'(8(I8,2X))')    NBR  
          WRITE (CPL(JW),'(8(I8,2X))')    IMX,KMX  
          DO JB=BS(JW),BE(JW)  
            WRITE (CPL(JW),'(9(I8,2X))')  US(JB),DS(JB)  
            WRITE (CPL(JW),'(9(I8,2X))')  KB(US(JB):DS(JB))  
          END DO  
          WRITE (CPL(JW),'(8(E13.6,2X))') DLX  
          WRITE (CPL(JW),'(8(E13.6,2X))') H  
          WRITE (CPL(JW),'(8(I8,2X))')    NAC  
          WRITE (CPL(JW),'(A)')           (CNAME1(CN(JN)),JN=1,NAC)       ! SW 3/1/2017
        ELSE  
          !c calculating longitudinal distance of segments  
          !           IF(JW == 1)DIST=0.0  
          !           do jb=BS(JW),BE(JW)  
          !               x1(us(jb))=dist+dlx(us(jb))/2.  
          !               DO I=US(JB)+1,DS(JB)  
          !                   DIST=DIST+(DLX(I)+dlx(i-1))/2.0  
          !                   X1(I)=DIST  
          !               END DO  
          !               DIST=DIST+DLX(DS(JB))  
          !               X1(DS(JB)+1)=DIST  
          !           ENDDO  
          WRITE (CPL(JW), *)'TITLE="CE-QUAL-W2"'  
          IF(HABTATC  == '      ON')THEN
              IF(NAC==0)THEN
                  WRITE (CPL(JW),19231)
                  ELSEIF(NACD(JW)==0)THEN
                  WRITE (CPL(JW),19232)(CNAME2(CN(JN)),JN=1,NAC)       !WRITE (CPL(JW),19233)(CNAME2(CN(JN)),JN=1,NAC)    SW 1/17/17
                  ELSEIF(NACD(JW)/=0)THEN
                  WRITE (CPL(JW),19233)(CNAME2(CN(JN)),JN=1,NAC),(CDNAME2(CDN(JD,JW)),JD=1,NACD(JW)) 
              ENDIF
          ELSE
               IF(NAC==0)THEN
                  WRITE (CPL(JW),19230)
                  ELSEIF(NACD(JW)==0)THEN
                  WRITE (CPL(JW),19234)(CNAME2(CN(JN)),JN=1,NAC)    !    1/17/17          !WRITE (CPL(JW),19233)(CNAME2(CN(JN)),JN=1,NAC)    SW 1/17/17
                  ELSEIF(NACD(JW)/=0)THEN
                  WRITE (CPL(JW),19235)(CNAME2(CN(JN)),JN=1,NAC),(CDNAME2(CDN(JD,JW)),JD=1,NACD(JW))    !    1/17/17      !WRITE (CPL(JW),19234)(CNAME2(CN(JN)),JN=1,NAC)  SW 9/28/13  
              ENDIF                                                    
          ENDIF
        19230 FORMAT('VARIABLES="Distance, m","Elevation, m","U(m/s)","W(m/s)","T(C)","RHO" ')
        19231 FORMAT('VARIABLES="Distance, m","Elevation, m","U(m/s)","W(m/s)","T(C)","RHO", "HABITAT" ')
        19232 FORMAT('VARIABLES="Distance, m","Elevation, m","U(m/s)","W(m/s)","T(C)","RHO", "HABITAT" ',<NAC>(',"',A8,'"'))
        19233 FORMAT('VARIABLES="Distance, m","Elevation, m","U(m/s)","W(m/s)","T(C)","RHO", "HABITAT" ',<NAC>(',"',A8,'"'),<NACD(JW)>(',"',A8,'"'))
        19234 FORMAT('VARIABLES="Distance, m","Elevation, m","U(m/s)","W(m/s)","T(C)","RHO" ',<NAC>(',"',A8,'"'))  ! SW 9/28/13
        19235 FORMAT('VARIABLES="Distance, m","Elevation, m","U(m/s)","W(m/s)","T(C)","RHO" ',<NAC>(',"',A8,'"'),<NACD(JW)>(',"',A8,'"'))  ! SW 9/28/13
        ENDIF  
      END IF  
        
      !IF (VECTOR(JW)) THEN     *** DSI  
      !  WRITE (VPL(JW),*)  TITLE  
      !  WRITE (VPL(JW),*)  H,KB,US,DS,DLX  
      !END IF  
        
    END DO  
      
    ! ***  
    IF (TIME_SERIES) THEN  
      L1 = SCAN(TSRFN1,'.',BACK=.TRUE.)    ! SW 8/22/14
      DO J=1,NIKTSR  
        WRITE (SEGNUM,'(I0)') ITSR(J)  
        SEGNUM = ADJUSTL(SEGNUM)  
        L      = LEN_TRIM(SEGNUM)  
        WRITE (SEGNUM2,'(I0)')J  
        SEGNUM2 = ADJUSTL(SEGNUM2)  
        L2      = LEN_TRIM(SEGNUM2)  
        TSRFN  = TSRFN1(1:L1-1)//'_'//SEGNUM2(1:L2)//'_seg'//SEGNUM(1:L)//'.'//TSRFN1(L1+1:L1+4)                                  !'.opt'
        OPEN  (TSR(J),FILE=TSRFN,STATUS='UNKNOWN')  
        !WRITE (TSR(J),'(A)') (TITLE(I),I=1,11)  
        I = ITSR(J)  ! SR 5/10/05
        DO JW=1,NWB  
          IF (I >= US(BS(JW)) .AND. I <= DS(BE(JW))) EXIT  
        END DO  

        IF (ICE_COMPUTATION) THEN  
          IF(SEDIMENT_CALC(JW))THEN  
            WRITE (TSR(J),'(*(A,","))') 'JDAY','DLT(s)','ELWS(m)','T2(C)','U(ms-1)','Q(m3s-1)','SRON(Wm-2)','EXT(m-1)',   &  
            'DEPTH(m)','WIDTH(m)','SHADE','ICETH(m)',                         &
            'Tvolavg(C)','NetRad(Wm-2)','SWSolar(Wm-2)','LWRad(Wm-2)','BackRad(Wm-2)','EvapF(Wm-2)','ConducF(Wm-2)','ReaerationCoeff(day-1)',   &
            (CNAME2(CN(JC)),JC=1,NAC),                                           &  
            ('     EPI',JE=1,NEP),('     MAC',JM=1,NMC),'     SED(Organic matter in sediments g/m3)','    SEDP(OrgP in sediments gP/m3)','    SEDN(OrgN in sediments gN/m3)','    SEDC(OrgCSediments gC/m3)',   &  
            (CDNAME2(CDN(JD,JW)),JD=1,NACD(JW)),(KFNAME2(KFCN(JF,JW)),JF=1,NAF(JW)),('PLIM_'//ADJUSTL(CNAME2(NAS+JA-1)),JA=1,NAL),('NLIM_'//ADJUSTL(CNAME2(NAS+JA-1)),JA=1,NAL),('LLIM_'//ADJUSTL(CNAME2(NAS+JA-1)),JA=1,NAL) 
          ELSE  
            WRITE (TSR(J),'(*(A,","))') 'JDAY','DLT(s)','ELWS(m)','T2(C)','U(ms-1)','Q(m3s-1)','SRON(Wm-2)','EXT(m-1)',   &  
            'DEPTH(m)','WIDTH(m)','SHADE','ICETH(m)','Tvolavg(C)','NetRad(Wm-2)','SWSolar(Wm-2)','LWRad(Wm-2)','BackRad(Wm-2)','EvapF(Wm-2)','ConducF(Wm-2)','ReaerationCoeff(day-1)',  (CNAME2(CN(JC)),JC=1,NAC),                     &  
            ('     EPI',JE=1,NEP),('     MAC',JM=1,NMC),(CDNAME2(CDN(JD,JW)),JD=1,NACD(JW)),(KFNAME2(KFCN(JF,JW)),JF=1,NAF(JW)),('PLIM_'//ADJUSTL(CNAME2(NAS+JA-1)),JA=1,NAL),('NLIM_'//ADJUSTL(CNAME2(NAS+JA-1)),JA=1,NAL),('LLIM_'//ADJUSTL(CNAME2(NAS+JA-1)),JA=1,NAL)     ! SW 10/20/15  
          END IF  
        ELSE  
          IF(SEDIMENT_CALC(JW))THEN  !mlm 7/25/06
            WRITE (TSR(J),'(*(A,","))')'JDAY','DLT(s)','ELWS(m)','T2(C)','U(ms-1)','Q(m3s-1)','SRON(Wm-2)','EXT(m-1)',   &  
            'DEPTH(m)','WIDTH(m)','SHADE','Tvolavg(C)','NetRad(Wm-2)','SWSolar(Wm-2)','LWRad(Wm-2)','BackRad(Wm-2)','EvapF(Wm-2)','ConducF(Wm-2)','ReaerationCoeff(day-1)',   &
            (CNAME2(CN(JC)),JC=1,NAC),                     &  
            ('     EPI',JE=1,NEP),('     MAC',JM=1,NMC),'     SED(Organic matter in sediments g/m3)','    SEDP(OrgP in sediments gP/m3)','    SEDN(OrgN in sediments gN/m3)','    SEDC(OrgCSediments gC/m3)',   &  
            (CDNAME2(CDN(JD,JW)),JD=1,NACD(JW)),(KFNAME2(KFCN(JF,JW)),JF=1,NAF(JW)),('PLIM_'//ADJUSTL(CNAME2(NAS+JA-1)),JA=1,NAL),('NLIM_'//ADJUSTL(CNAME2(NAS+JA-1)),JA=1,NAL),('LLIM_'//ADJUSTL(CNAME2(NAS+JA-1)),JA=1,NAL)   
          ELSE  
            WRITE (TSR(J),'(*(A,","))') 'JDAY','DLT(s)','ELWS(m)','T2(C)','U(ms-1)','Q(m3s-1)','SRON(Wm-2)','EXT(m-1)',   &  
             'DEPTH(m)','WIDTH(m)','SHADE','Tvolavg(C)','NetRad(Wm-2)','SWSolar(Wm-2)','LWRad(Wm-2)','BackRad(Wm-2)','EvapF(Wm-2)','ConducF(Wm-2)','ReaerationCoeff(day-1)',   &
            (CNAME2(CN(JC)),JC=1,NAC),                     &  
            ('     EPI',JE=1,NEP),('     MAC',JM=1,NMC),(CDNAME2(CDN(JD,JW)),JD=1,NACD(JW)),(KFNAME2(KFCN(JF,JW)),JF=1,NAF(JW)),('PLIM_'//ADJUSTL(CNAME2(NAS+JA-1)),JA=1,NAL),('NLIM_'//ADJUSTL(CNAME2(NAS+JA-1)),JA=1,NAL),('LLIM_'//ADJUSTL(CNAME2(NAS+JA-1)),JA=1,NAL)  
          END IF  
        END IF  
          
      END DO  
    END IF  
    IF(ALGAE_TOXIN)THEN
      IF(ATOX_DEBUG=='ON')THEN
          OPEN(ATOXIN_DEBUG_FN,FILE='algae_toxin_debug.csv',STATUS='UNKNOWN')
          WRITE(ATOXIN_DEBUG_FN,'(A,<NUMATOXINS>("EX_TOXIN_",I1,","),<NUMATOXINS>("IN_TOXIN_",I1,","),<NUMATOXINS>("CTESS_TOXIN_",I1,","),<NAL>("ALG_",I1,","))')'JDAY,K,I,',(J,J=1,numatoxins),(J,J=1,numatoxins),(J,J=1,numatoxins),(J,J=1,NAL)
      ENDIF
    ENDIF
    ! BIOENERGETICS mlm  
    IF (BIOEXP) THEN  
      L1 = SCAN(BIOFN,'.')  
      DO J=1,NIBIO  
        WRITE (SEGNUM,'(I0)') J  
        SEGNUM = ADJUSTL(SEGNUM)  
        L      = LEN_TRIM(SEGNUM)  
        BIOFN  = BIOFN(1:L1-1)//'_'//SEGNUM(1:L)//'.opt'  
        OPEN  (BIOEXPFN(J),FILE=BIOFN,STATUS='UNKNOWN')  
        WRITE (bioexpfn(J),'(12A8)') (adjustr(bhead(ii)),ii=1,11)  
      END DO  
      L1 = SCAN(WEIGHTFN,'.')  
      DO J=1,NIBIO  
        WRITE (SEGNUM,'(I0)') J  
        SEGNUM = ADJUSTL(SEGNUM)  
        L      = LEN_TRIM(SEGNUM)  
        WEIGHTFN  = WEIGHTFN(1:L1-1)//'_'//SEGNUM(1:L)//'.opt'  
        OPEN  (WEIGHTNUM(J),FILE=WEIGHTFN,STATUS='UNKNOWN')  
        WRITE (weightnum(J),'(50(A10,","))') adjustr('jday'),(adjustr(cname2(cn(jc))),jc=1,nac),adjustr('pH'),adjustr('TP')  
      END DO  
        
    END IF  
    IF (DOWNSTREAM_OUTFLOW) THEN  
      L1 = SCAN(WDOFN,'.',BACK=.TRUE.)   ! SW 8/22/14 CHECKS FROM RIGHTHAND SIDE NOT LEFTHANDSIDE  ! SW 8/3/2018  
      DO JWD=1,NIWDO  
        WRITE (SEGNUM,'(I0)') IWDO(JWD)  
        SEGNUM = ADJUSTL(SEGNUM)  
        L      = LEN_TRIM(SEGNUM)  
        OPEN  (WDO(JWD,1),FILE='qwo_'//SEGNUM(1:L)//WDOFN(L1:L1+4),STATUS='UNKNOWN')                   !      OPEN  (WDO(JWD,1),FILE='qwo_'//SEGNUM(1:L)//'.opt',STATUS='UNKNOWN')
        OPEN  (WDO(JWD,2),FILE='two_'//SEGNUM(1:L)//WDOFN(L1:L1+4),STATUS='UNKNOWN')                   !'.opt',STATUS='UNKNOWN')
        WRITE (WDO(JWD,1),'(A,I0/A/A)') '$Flow file for segment ',       IWDO(JWD),'To the right of the sum of flows are individual flows starting with QWD then QSTR','JDAY,QWD(m3s-1),'  
        WRITE (WDO(JWD,2),'(A,I0/A/A)') '$Temperature file for segment ',IWDO(JWD),'To the right of the sum of temperatures are individual temperatures starting with QWD then QSTR','JDAY,T(C),'  
        DO JW=1,NWB  
          IF (IWDO(JWD) >= US(BS(JW)) .AND. IWDO(JWD) <= DS(BE(JW))) EXIT  
        END DO  
        IF (CONSTITUENTS) THEN  
          OPEN  (WDO(JWD,3),FILE='cwo_'//SEGNUM(1:L)//WDOFN(L1:L1+4),STATUS='UNKNOWN')                 !'.opt',STATUS='UNKNOWN')
          WRITE (WDO(JWD,3),'(A,I0//(*(A,",")))') '$Concentration file for segment ',IWDO(JWD),'JDAY', (CNAME2(CN(J)),J=1,NAC)                !CNAME2(CN(1:NAC))  
        END IF  
        IF (DERIVED_CALC) THEN  
          OPEN  (WDO(JWD,4),FILE='dwo_'//SEGNUM(1:L)//WDOFN(L1:L1+4),STATUS='UNKNOWN')                 !'.opt',STATUS='UNKNOWN')
          WRITE (WDO(JWD,4),'(A,I0//(*(A,",")))') 'Derived constituent file for segment ',IWDO(JWD),'JDAY', (CDNAME2(CDN(J,JW)),J=1,NACD(JW))                 !CDNAME2(CDN(1:NACD(JW),JW))  
        END IF  
      END DO  
    END IF  
      
    IF (DOWNSTREAM_OUTFLOW) THEN  
      JFILE=0  
      DO JWD=1,NIWDO  
          
        ! Determine the # of withdrawals at the WITH SEG  
        DO JB=1,NBR  ! structures
          IF(IWDO(JWD)==DS(JB) .AND. NSTR(JB) /= 0)THEN  
            DO JS=1,NSTR(JB)  
              JFILE=JFILE+1  
              WRITE(SEGNUM2,'(I0)')IWDO(JWD)  
              SEGNUM2 = ADJUSTL(SEGNUM2)  
              L2      = LEN_TRIM(SEGNUM2)  
              TITLEWITH2='$STR WITHDRAWAL AT SEG'//SEGNUM2(1:L2)  
              WRITE(SEGNUM,'(F10.0)')ESTR(JS,JB)  
              SEGNUM = ADJUSTL(SEGNUM)  
              L      = LEN_TRIM(SEGNUM)  
              TITLEWITH=ADJUSTL(TITLEWITH2)//' ELEV OF WITHDRAWAL CENTERLINE:'//SEGNUM(1:L)  
              WRITE(SEGNUM,'(I0)')JS  
              SEGNUM = ADJUSTL(SEGNUM)  
              L      = LEN_TRIM(SEGNUM)  
              WDO2(JFILE,1) = NUNIT; NUNIT = NUNIT+1  
            OPEN(WDO2(JFILE,1),FILE='qwo_str'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),STATUS='UNKNOWN')      !'.opt',STATUS='UNKNOWN')
              WRITE (WDO2(JFILE,1),'(A//A)') TITLEWITH,'JDAY,QWD(m3s-1)'  
              WDO2(JFILE,2) = NUNIT; NUNIT = NUNIT+1  
            OPEN(WDO2(JFILE,2),FILE='two_str'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),STATUS='UNKNOWN')      !'.opt',STATUS='UNKNOWN')
              WRITE (WDO2(JFILE,2),'(A//A)') TITLEWITH,'JDAY,T(C)'  
              IF (CONSTITUENTS) THEN  
                WDO2(JFILE,3) = NUNIT; NUNIT = NUNIT+1  
                OPEN(WDO2(JFILE,3),FILE='cwo_str'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),STATUS='UNKNOWN')   !'.opt',STATUS='UNKNOWN')
                WRITE (WDO2(JFILE,3),'(A//(*(A,",")))') TITLEWITH,'JDAY',(CNAME2(CN(J)),J=1,NAC)     !CNAME2(CN(1:NAC))  
              ENDIF  
              IF (DERIVED_CALC) THEN  
                WDO2(JFILE,4) = NUNIT; NUNIT = NUNIT+1  
                OPEN(WDO2(JFILE,4),FILE='dwo_str'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),STATUS='UNKNOWN')     !'.opt',STATUS='UNKNOWN')
                DO JW=1,NWB  
                  IF (IWDO(JWD) >= US(BS(JW)) .AND. IWDO(JWD) <= DS(BE(JW))) EXIT  
                END DO  
                WRITE (WDO2(JFILE,4),'(A//(*(A,",")))') TITLEWITH,'JDAY',(CDNAME2(CDN(J,JW)),J=1,NACD(JW))  !CDNAME2(CDN(1:NACD(JW),JW))  
              ENDIF  
            ENDDO  
          ENDIF  
        ENDDO  
          
        DO JS=1,NWD  ! withdrawals
          IF(IWDO(JWD) == IWD(JS))THEN  
            JFILE=JFILE+1  
            WRITE(SEGNUM2,'(I0)')IWDO(JWD)  
            SEGNUM2 = ADJUSTL(SEGNUM2)  
            L2      = LEN_TRIM(SEGNUM2)  
            TITLEWITH2='$WITHDRAWAL AT SEG'//SEGNUM2(1:L2)  
            WRITE(SEGNUM,'(F10.0)')EWD(JS)  
            SEGNUM = ADJUSTL(SEGNUM)  
            L      = LEN_TRIM(SEGNUM)  
            TITLEWITH=ADJUSTL(TITLEWITH2)//' ELEV OF WITHDRAWAL CENTERLINE:'//SEGNUM(1:L)  
            WRITE(SEGNUM,'(I0)')JS  
            SEGNUM = ADJUSTL(SEGNUM)  
            L      = LEN_TRIM(SEGNUM)  
            WDO2(JFILE,1) = NUNIT; NUNIT = NUNIT+1  
            OPEN(WDO2(JFILE,1),FILE='qwo_wd'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),STATUS='UNKNOWN')        !'.opt',STATUS='UNKNOWN')
            WRITE (WDO2(JFILE,1),'(A//A)') TITLEWITH,'JDAY,QWD(m3s-1)'  
            WDO2(JFILE,2) = NUNIT; NUNIT = NUNIT+1  
            OPEN(WDO2(JFILE,2),FILE='two_wd'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),STATUS='UNKNOWN')         !'.opt',STATUS='UNKNOWN')
            WRITE (WDO2(JFILE,2),'(A//A)') TITLEWITH,'JDAY,T(C)'  
            IF (CONSTITUENTS) THEN  
              WDO2(JFILE,3) = NUNIT; NUNIT = NUNIT+1  
                OPEN(WDO2(JFILE,3),FILE='cwo_wd'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),STATUS='UNKNOWN')     !'.opt',STATUS='UNKNOWN')
              WRITE (WDO2(JFILE,3),'(A//(*(A,",")))') TITLEWITH,'JDAY',(CNAME2(CN(J)),J=1,NAC)     !CNAME2(CN(1:NAC))  
            ENDIF  
            IF (DERIVED_CALC) THEN  
              WDO2(JFILE,4) = NUNIT; NUNIT = NUNIT+1  
                OPEN(WDO2(JFILE,4),FILE='dwo_wd'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),STATUS='UNKNOWN')     !'.opt',STATUS='UNKNOWN')
              DO JW=1,NWB  
                IF (IWDO(JWD) >= US(BS(JW)) .AND. IWDO(JWD) <= DS(BE(JW))) EXIT  
              END DO  
              WRITE (WDO2(JFILE,4),'(A//(*(A,",")))') TITLEWITH,'JDAY',(CDNAME2(CDN(J,JW)),J=1,NACD(JW))   !CDNAME2(CDN(1:NACD(JW),JW))  
            ENDIF  
          ENDIF  
        ENDDO  
          
        DO JS=1,NSP  ! spillways
          IF(IWDO(JWD) == IUSP(JS))THEN  
            JFILE=JFILE+1  
            WRITE(SEGNUM2,'(I0)')IWDO(JWD)  
            SEGNUM2 = ADJUSTL(SEGNUM2)  
            L2      = LEN_TRIM(SEGNUM2)  
            TITLEWITH2='$SPILLWAY WITHDRAWAL AT SEG'//SEGNUM2(1:L2)  
            WRITE(SEGNUM,'(F10.0)')ESP(JS)  
            SEGNUM = ADJUSTL(SEGNUM)  
            L      = LEN_TRIM(SEGNUM)  
            TITLEWITH=ADJUSTL(TITLEWITH2)//' ELEV OF WITHDRAWAL CENTERLINE:'//SEGNUM(1:L)  
            WRITE(SEGNUM,'(I0)')JS  
            SEGNUM = ADJUSTL(SEGNUM)  
            L      = LEN_TRIM(SEGNUM)  
            WDO2(JFILE,1) = NUNIT; NUNIT = NUNIT+1  
            OPEN(WDO2(JFILE,1),FILE='qwo_sp'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),STATUS='UNKNOWN')      !'.opt',STATUS='UNKNOWN')
            WRITE (WDO2(JFILE,1),'(A//A)') TITLEWITH,'JDAY,QWD(m3s-1)'  
            WDO2(JFILE,2) = NUNIT; NUNIT = NUNIT+1  
            OPEN(WDO2(JFILE,2),FILE='two_sp'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),STATUS='UNKNOWN')      !'.opt',STATUS='UNKNOWN')
            WRITE (WDO2(JFILE,2),'(A//A)') TITLEWITH,'JDAY,T(C)'  
            IF (CONSTITUENTS) THEN  
              WDO2(JFILE,3) = NUNIT; NUNIT = NUNIT+1  
                OPEN(WDO2(JFILE,3),FILE='cwo_sp'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),STATUS='UNKNOWN')     !'.opt',STATUS='UNKNOWN')
              WRITE (WDO2(JFILE,3),'(A//(*(A,",")))') TITLEWITH,'JDAY',(CNAME2(CN(J)),J=1,NAC)     !CNAME2(CN(1:NAC))  
            ENDIF  
            IF (DERIVED_CALC) THEN  
              WDO2(JFILE,4) = NUNIT; NUNIT = NUNIT+1  
                OPEN(WDO2(JFILE,4),FILE='dwo_sp'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),STATUS='UNKNOWN')     !'.opt',STATUS='UNKNOWN')
              DO JW=1,NWB  
                IF (IWDO(JWD) >= US(BS(JW)) .AND. IWDO(JWD) <= DS(BE(JW))) EXIT  
              END DO  
              WRITE (WDO2(JFILE,4),'(A//(*(A,",")))') TITLEWITH,'JDAY',(CDNAME2(CDN(J,JW)),J=1,NACD(JW))    !CDNAME2(CDN(1:NACD(JW),JW))  
            ENDIF  
          ENDIF  
        ENDDO  
          
          
        DO JS=1,NPU  ! pumps
          IF(IWDO(JWD) == IUPU(JS))THEN  
            JFILE=JFILE+1  
            WRITE(SEGNUM2,'(I0)')IWDO(JWD)  
            SEGNUM2 = ADJUSTL(SEGNUM2)  
            L2      = LEN_TRIM(SEGNUM2)  
            TITLEWITH2='$PUMP WITHDRAWAL AT SEG'//SEGNUM2(1:L2)  
            WRITE(SEGNUM,'(F10.0)')EPU(JS)  
            SEGNUM = ADJUSTL(SEGNUM)  
            L      = LEN_TRIM(SEGNUM)  
            TITLEWITH=ADJUSTL(TITLEWITH2)//' ELEV OF WITHDRAWAL CENTERLINE:'//SEGNUM(1:L)  
            WRITE(SEGNUM,'(I0)')JS  
            SEGNUM = ADJUSTL(SEGNUM)  
            L      = LEN_TRIM(SEGNUM)  
            WDO2(JFILE,1) = NUNIT; NUNIT = NUNIT+1  
            OPEN(WDO2(JFILE,1),FILE='qwo_pmp'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),STATUS='UNKNOWN')     !'.opt',STATUS='UNKNOWN')
            WRITE (WDO2(JFILE,1),'(A//A)') TITLEWITH,'JDAY,QWD(m3s-1)'  
            WDO2(JFILE,2) = NUNIT; NUNIT = NUNIT+1  
            OPEN(WDO2(JFILE,2),FILE='two_pmp'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),STATUS='UNKNOWN')     !'.opt',STATUS='UNKNOWN')
            WRITE (WDO2(JFILE,2),'(A//A)') TITLEWITH,'JDAY,T(C)'  
            IF (CONSTITUENTS) THEN  
              WDO2(JFILE,3) = NUNIT; NUNIT = NUNIT+1  
                OPEN(WDO2(JFILE,3),FILE='cwo_pmp'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),STATUS='UNKNOWN')   !'.opt',STATUS='UNKNOWN')
              WRITE (WDO2(JFILE,3),'(A//(*(A,",")))') TITLEWITH,'JDAY',(CNAME2(CN(J)),J=1,NAC)      !CNAME2(CN(1:NAC))  
            ENDIF  
            IF (DERIVED_CALC) THEN  
              WDO2(JFILE,4) = NUNIT; NUNIT = NUNIT+1  
                OPEN(WDO2(JFILE,4),FILE='dwo_pmp'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),STATUS='UNKNOWN')    !'.opt',STATUS='UNKNOWN')
              DO JW=1,NWB  
                IF (IWDO(JWD) >= US(BS(JW)) .AND. IWDO(JWD) <= DS(BE(JW))) EXIT  
              END DO  
              WRITE (WDO2(JFILE,4),'(A//(*(A,",")))') TITLEWITH,'JDAY',(CDNAME2(CDN(J,JW)),J=1,NACD(JW))     !CDNAME2(CDN(1:NACD(JW),JW))  
            ENDIF  
          ENDIF  
        ENDDO  
          
          
        DO JS=1,NPI  ! pipes
          IF(IWDO(JWD) == IUPI(JS))THEN  
            JFILE=JFILE+1  
            WRITE(SEGNUM2,'(I0)')IWDO(JWD)  
            SEGNUM2 = ADJUSTL(SEGNUM2)  
            L2      = LEN_TRIM(SEGNUM2)  
            TITLEWITH2='$PIPE WITHDRAWAL AT SEG'//SEGNUM2(1:L2)  
            WRITE(SEGNUM,'(F10.0)')EUPI(JS)  
            SEGNUM = ADJUSTL(SEGNUM)  
            L      = LEN_TRIM(SEGNUM)  
            TITLEWITH=ADJUSTL(TITLEWITH2)//' ELEV OF WITHDRAWAL CENTERLINE:'//SEGNUM(1:L)  
            WRITE(SEGNUM,'(I0)')JS  
            SEGNUM = ADJUSTL(SEGNUM)  
            L      = LEN_TRIM(SEGNUM)  
            WDO2(JFILE,1) = NUNIT; NUNIT = NUNIT+1  
            OPEN(WDO2(JFILE,1),FILE='qwo_pipe'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),STATUS='UNKNOWN')     !'.opt',STATUS='UNKNOWN')
            WRITE (WDO2(JFILE,1),'(A//A)') TITLEWITH,'JDAY,QWD(m3s-1)'  
            WDO2(JFILE,2) = NUNIT; NUNIT = NUNIT+1  
            OPEN(WDO2(JFILE,2),FILE='two_pipe'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),STATUS='UNKNOWN')     !'.opt',STATUS='UNKNOWN')
            WRITE (WDO2(JFILE,2),'(A//A)') TITLEWITH,'JDAY,T(C)'  
            IF (CONSTITUENTS) THEN  
              WDO2(JFILE,3) = NUNIT; NUNIT = NUNIT+1  
                OPEN(WDO2(JFILE,3),FILE='cwo_pipe'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),STATUS='UNKNOWN')    !'.opt',STATUS='UNKNOWN')
              WRITE (WDO2(JFILE,3),'(A//(*(A,",")))') TITLEWITH,'JDAY',(CNAME2(CN(J)),J=1,NAC)      !CNAME2(CN(1:NAC))  
            ENDIF  
            IF (DERIVED_CALC) THEN  
              WDO2(JFILE,4) = NUNIT; NUNIT = NUNIT+1  
                OPEN(WDO2(JFILE,4),FILE='dwo_pipe'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),STATUS='UNKNOWN')    !'.opt',STATUS='UNKNOWN')
              DO JW=1,NWB  
                IF (IWDO(JWD) >= US(BS(JW)) .AND. IWDO(JWD) <= DS(BE(JW))) EXIT  
              END DO  
              WRITE (WDO2(JFILE,4),'(A//(*(A,",")))') TITLEWITH,'JDAY',(CDNAME2(CDN(J,JW)),J=1,NACD(JW))      !CDNAME2(CDN(1:NACD(JW),JW))  
            ENDIF  
          ENDIF  
        ENDDO  
          
        DO JS=1,NGT  ! gates
          IF(IWDO(JWD) == IUGT(JS))THEN  
            JFILE=JFILE+1  
            WRITE(SEGNUM2,'(I0)')IWDO(JWD)  
            SEGNUM2 = ADJUSTL(SEGNUM2)  
            L2      = LEN_TRIM(SEGNUM2)  
            TITLEWITH2='$GATE WITHDRAWAL AT SEG'//SEGNUM2(1:L2)  
            WRITE(SEGNUM,'(F10.0)')EGT(JS)  
            SEGNUM = ADJUSTL(SEGNUM)  
            L      = LEN_TRIM(SEGNUM)  
            TITLEWITH=ADJUSTL(TITLEWITH2)//' ELEV OF WITHDRAWAL CENTERLINE:'//SEGNUM(1:L)  
            WRITE(SEGNUM,'(I0)')JS  
            SEGNUM = ADJUSTL(SEGNUM)  
            L      = LEN_TRIM(SEGNUM)  
            WDO2(JFILE,1) = NUNIT; NUNIT = NUNIT+1  
            OPEN(WDO2(JFILE,1),FILE='qwo_gate'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),STATUS='UNKNOWN')     !'.opt',STATUS='UNKNOWN')
            WRITE (WDO2(JFILE,1),'(A//A)') TITLEWITH,'JDAY,QWD(m3s-1)'  
            WDO2(JFILE,2) = NUNIT; NUNIT = NUNIT+1  
            OPEN(WDO2(JFILE,2),FILE='two_gate'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),STATUS='UNKNOWN')     !'.opt',STATUS='UNKNOWN')
            WRITE (WDO2(JFILE,2),'(A//A)') TITLEWITH,'JDAY,T(C)'  
            IF (CONSTITUENTS) THEN  
              WDO2(JFILE,3) = NUNIT; NUNIT = NUNIT+1  
                OPEN(WDO2(JFILE,3),FILE='cwo_gate'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),STATUS='UNKNOWN')    !'.opt',STATUS='UNKNOWN')
              WRITE (WDO2(JFILE,3),'(A//(*(A,",")))') TITLEWITH,'JDAY',  (CNAME2(CN(J)),J=1,NAC)      !CNAME2(CN(1:NAC))  
            ENDIF  
            IF (DERIVED_CALC) THEN  
              WDO2(JFILE,4) = NUNIT; NUNIT = NUNIT+1  
                OPEN(WDO2(JFILE,4),FILE='dwo_gate'//SEGNUM(1:L)//'_seg'//SEGNUM2(1:L2)//WDOFN(L1:L1+4),STATUS='UNKNOWN')     !'.opt',STATUS='UNKNOWN')
              DO JW=1,NWB  
                IF (IWDO(JWD) >= US(BS(JW)) .AND. IWDO(JWD) <= DS(BE(JW))) EXIT  
              END DO  
              WRITE (WDO2(JFILE,4),'(A//(*(A,",")))') TITLEWITH,'JDAY',(CDNAME2(CDN(J,JW)),J=1,NACD(JW))       !CDNAME2(CDN(1:NACD(JW),JW))  
            ENDIF  
          ENDIF  
        ENDDO  
          
      END DO  
    END IF  
      
     ! OUTPUT FLOW LOADING AND POLLUTANT LOADING DEBUGGING INFORMATION OUTPUT FREQUENCY AT OUTPUT OF CPL
    DO JW=1,NWB  
      IF (FLOWBALC=='      ON') THEN  
        OPEN(FLOWBFN,FILE='flowbal.csv',STATUS='UNKNOWN')  
        IF(VOLUME_BALANCE(JW))THEN
            WRITE(FLOWBFN,'(a102)')'JDAY,WB,VOLIN(m3),VOLPR(m3),VOLOUT(m3),VOLWD(m3),VOLEV(m3),VOLDT(m3),VOLTRB(m3),VOLICE(m3),%VOLerror,'
        ELSE
            WRITE(FLOWBFN,'(a93)')'JDAY,WB,VOLIN(m3),VOLPR(m3),VOLOUT(m3),VOLWD(m3),VOLEV(m3),VOLDT(m3),VOLTRB(m3),VOLICE(m3),'
        ENDIF
        EXIT  
      ENDIF  
    ENDDO  
        DO JW=1,NWB  
        IF (NPBALC=='      ON') THEN
        OPEN(MASSBFN,FILE='massbal.csv',STATUS='UNKNOWN')  
        IF(SEDIMENT_DIAGENESIS)THEN
            WRITE(MASSBFN,'(A574)')'JDAY,WB,TP-Waterbody(kg),TP-Sediment(kg),TP-Plants(kg),OutflowTP(kg),TributaryTP(kg),DistributedTributaryTP(kg),WithdrawalTP(kg),PrecipitationTP(kg),InflowTP(kg),AtmosphericDepositionP(kg),SED+SOD_PRelease(kg),PFluxtoSediments(kg),SedimentDiagenesisPFlux(kg),TN-Waterbody(kg),TN-Sediment(kg),TN-Plants(kg),OutflowTN(kg),TributaryTN(kg),DistributedTributaryTN(kg),WithdrawalTN(kg),PrecipitationTN(kg),InflowTN(kg),AtmosphericDepositionN(kg),NH3GasLoss(kg),SED+SOD_NRelease(kg),NFluxtoSediments(kg),SedimentDiagenesisNH4Flux(kg),SedimentDiagenesisNO3Flux(kg)'   ! 'JDAY,WB,TP-Waterbody(kg),TP-Sediment(kg),TP-Plants(kg),OutflowTP(kg),TributaryTP(kg),DistributedTributaryTP(kg),WithdrawalTP(kg),PrecipitationTP(kg),InflowTP(kg),BurialP_SED(kg),SED+SOD_PRelease(kg),TN-Waterbody(kg),TN-Sediment(kg),TN-Plants(kg),OutflowTN(kg),TributaryTN(kg),DistributedTributaryTN(kg),WithdrawalTN(kg),PrecipitationTN(kg),InflowTN(kg),BurialN_SED(kg),SED+SOD_NRelease(kg)'
            ELSE
            WRITE(MASSBFN,'(A467)')'JDAY,WB,TP-Waterbody(kg),TP-Sediment(kg),TP-Plants(kg),OutflowTP(kg),TributaryTP(kg),DistributedTributaryTP(kg),WithdrawalTP(kg),PrecipitationTP(kg),InflowTP(kg),AtmosphericDepositionP(kg),SED+SOD_PRelease(kg),PFluxtoSediments(kg),TN-Waterbody(kg),TN-Sediment(kg),TN-Plants(kg),OutflowTN(kg),TributaryTN(kg),DistributedTributaryTN(kg),WithdrawalTN(kg),PrecipitationTN(kg),InflowTN(kg),AtmosphericDepositionN(kg),NH3GasLoss(kg),SED+SOD_NRelease(kg),NFluxtoSediments(kg)'   ! 'JDAY,WB,TP-Waterbody(kg),TP-Sediment(kg),TP-Plants(kg),OutflowTP(kg),TributaryTP(kg),DistributedTributaryTP(kg),WithdrawalTP(kg),PrecipitationTP(kg),InflowTP(kg),BurialP_SED(kg),SED+SOD_PRelease(kg),TN-Waterbody(kg),TN-Sediment(kg),TN-Plants(kg),OutflowTN(kg),TributaryTN(kg),DistributedTributaryTN(kg),WithdrawalTN(kg),PrecipitationTN(kg),InflowTN(kg),BurialN_SED(kg),SED+SOD_NRelease(kg)'
        ENDIF
            
        EXIT  
      ENDIF  
    ENDDO  
      
    IF(WLC=='      ON')THEN  
      OPEN(WLFN,FILE='wl.csv',STATUS='UNKNOWN')  
      write(WLFN,'("JDAY,",*("SEG",i4,","))')((i,i=us(jb),ds(jb)),jb=1,nbr)         
    ENDIF  
      
    !**** DSI W2 Linkage File (W2L) (Supercedes Old Velocity vectors)  
    IF (VECTOR(1)) THEN  
      ! *** Apply the same linkage settings for all waterbodies  
      OPEN (VPL(1),FILE=VPLFN(1),STATUS='UNKNOWN',ACCESS='SEQUENTIAL',FORM='BINARY')  
        
      ! *** W2 Version  
      WRITE(VPL(1)) W2VER  
        
      WRITE(VPL(1)) TITLE  
        
      WRITE(VPL(1)) INT4(NWB), INT4(NBR), INT4(IMX), INT4(KMX), INT4(NCT), INT4(NAC)  
        
      ! *** Flag the output file if using Outlet time series  
      IF (TIME_SERIES) THEN  
        IFLAG =  -1  
      ELSE  
        IFLAG =  0  
      END IF  
      WRITE(VPL(1)) IFLAG  
        
      ! *** MODEL CONFIGURATION  
      DO JW=1,NWB  
        WRITE(VPL(1)) INT4(BS(JW)),INT4(BE(JW))  
      END DO  
        
      WRITE(VPL(1)) (REAL(DLX(I),4) ,I=1,IMX)  
      WRITE(VPL(1)) (REAL(PHI0(I),4),I=1,IMX)  
      WRITE(VPL(1)) (INT4(US(K)) ,K=1,NBR)  
      WRITE(VPL(1)) (INT4(DS(K)) ,K=1,NBR)  
      WRITE(VPL(1)) (INT4(UHS(K)),K=1,NBR)  
      WRITE(VPL(1)) (INT4(DHS(K)),K=1,NBR)  
      WRITE(VPL(1)) ((REAL(H(K,JW),4) ,K=1,KMX),JW=1,NWB)  
        
      ! *** INITIALIZE ALL CELLS AS INACTIVE (SET FLAG = 0)  
      ICOMP = 0  
        
      ! *** Set up the array to handle cell types : Outflows  
      DO JB = 1, NBR  
        IF (NSTR(JB) > 0 ) THEN  
          DO JS=1,NSTR(JB)  
            DO K=KMX-1,2,-1  
              IF (ESTR(JS,JB).LE.EL(K,DS(JB))) THEN  
                ICOMP(K,DS(JB)) = 4  
              END IF  
            END DO  
          END DO  
    !    ELSE      ! SW NEEDS FIXED SINCE WE DO NOT HAVE NOUT AND KOUT DEFINED 9/28/13
    !     DO JO = 1, NOUT(JB)  
    !        ICOMP(KOUT(JO,JB),DS(JB)) = 3  
    !      END DO  
        END IF  
      END DO  
        
      ! *** Set up the array to handle cell types : Withdrawals  
      DO IW=1,NWD
        !DO K=KTW(IW),KBW(IW),-1  
        DO K=KMX,2,-1  
          IF (EWD(IW) <= EL(K,IWD(IW))) THEN  
            ICOMP(K-1,IWD(IW)) = 2
            EXIT
          END IF  
        END DO  
      END DO  
        
      ! *** Set up the array to handle cell types : Tributaries  
      IF (NTR.GT.0) THEN  
        DO JT=1,NTR  
          DO K=2,KB(ITR(JT))  
            IF (ICOMP(K,ITR(JT)) == 0) THEN  
              ICOMP(K,ITR(JT)) = 5  
            ENDIF  
          END DO
        END DO  
      ENDIF  
        
      ! *** FLAG ACTIVE CELLS  
      ! *** (IF BATHYMETRY(B)>0, FLAG CELL WITH VALUE OF 1, IF IT HAS  
      ! *** NOT BEEN SET ABOVE TO A VALUE OF 2,3, OR 4 )  
      DO JB=1,NBR  
        IU = US(JB)  
        ID = DS(JB)  
        DO I=IU,ID  
          DO K=2,KB(I)  
            IF((B(K,I) > 0.).AND.(ICOMP(K,I) == 0)) ICOMP(K,I) = 1  
          END DO  
        END DO  
      END DO  
        
      ! *** NOW WRITE SEGMENT WIDTHS (B)  
      DO I = 1, IMX  
        DO K = 1,KMX  
          IF(ICOMP(K,I) > 0)THEN  
            WRITE(VPL(1)) REAL(B(K,I),4)  
          ELSE  
            WRITE(VPL(1)) 0.0_4  
          ENDIF  
        ENDDO  
      ENDDO  
        
      ! *** Output X -  
      ! *** THIS SEQUENCE SWEEPS ACROSS COLUMNS WITH BLOCK CENTERED  
      ! *** X COORDS WORKING UPSTREAM FOR EACH BRANCH  
      DO K = 1,KMX  
        DO JB = 1, NBR  
          IU = US(JB)-1  
          ID = DS(JB)+1  
          DO I = ID,IU,-1  
            IF (I.EQ.ID) THEN  
              DLX_OLD = 0.5*DLX(I)  
            ELSE  
              DLX_OLD = DLX_OLD + 0.5*DLX(I) + 0.5*DLX(I+1)  
            END IF  
            WRITE(VPL(1)) REAL(DLX_OLD,4)  
          END DO  
        END DO  
      END DO  
        
      ! *** Output Y  
      ! *** THIS SEQUENCE SWEEPS ACROSS COLUMNS WITH BLOCK CENTERED  
      ! *** ELEVATIONS WORKING FROM THE BOTTOM LAYER TO THE TOP  
      DO K = 1,KMX  
        DO JB = 1, NBR  
          ! *** BOTTOM INACTIVE CELL  
          !ELC = ELBOT(JW)-0.5*H(KMX,JW)  
          !WRITE(VPL(1)) ELC  
            
          IU = US(JB)-1  
          ID = DS(JB)+1  
          DO I = ID,IU,-1  
            WRITE(VPL(1)) REAL(EL(K,I),4)  
          END DO  
        END DO  
      END DO  
        
      WRITE(VPL(1)) ((ICOMP(K,I),K=1,KMX),I=1,IMX)  
        
      ! *** CONSTITUENT NUMBER LIST  
      WRITE (VPL(1)) (INT4(CN(JC)),JC=1,NAC)  
        
      ! *** CONSTITUENT NAME  
      WRITE (VPL(1)) (CNAME1(CN(JC)),JC=1,NAC)  
        
      ! *** CONSTITUENT UNITS  
      WRITE (VPL(1)) (CUNIT2(CN(JC)),JC=1,NAC)  
        
      IF (TIME_SERIES) THEN  
        !OPEN(ITSRNAG,FILE='ITSRNAG.TMP',STATUS='SCRATCH', ACCESS='SEQUENTIAL',FORM='BINARY')    PMC - FIX LATER
      END IF  
        
    ENDIF  
      
  END IF  
    
  RETURN  
END SUBROUTINE OUTPUTINIT  

