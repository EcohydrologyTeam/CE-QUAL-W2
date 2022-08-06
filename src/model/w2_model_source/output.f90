
!***********************************************************************************************************************************
!**                                           S U B R O U T I N E   O U T P U T                                                   **
!***********************************************************************************************************************************

SUBROUTINE OUTPUT (JDAY,IUPR,IDPR,KBR,ISNP,BL,NBL)
  USE GLOBAL; USE GDAYC;  USE GEOMC;  USE KINETIC; USE TVDC; USE NAMESC; USE LOGICC
  USE MACROPHYTEC
  Use CEMAVars  
  IMPLICIT NONE

! Type declaration

  REAL                        :: JDAY, LIMIT
  INTEGER                     :: NBL,IUPR,IDPR,JH,L,K,J,NLINES,KBR,JAC,JD,JE,JT,JJ,JA
  INTEGER, DIMENSION(IMX)     :: BL
  INTEGER, DIMENSION(IMX,NWB) :: ISNP
  LOGICAL                     :: NEW_PAGE
  CHARACTER(8)                :: LFAC
  CHARACTER(10)               :: BLANK='          '

! Variable initialization

  NEW_PAGE = .TRUE.

! Blank inactive cells

  NBL  = 1
  JB   = 1
  IUPR = 1
  DO I=1,IDPR-1
    IF (CUS(JB) > ISNP(I,JW)) THEN
      BL(NBL) = I
      NBL     = NBL+1
      IF (JB == 1) IUPR = I+1
    END IF
    IF (ISNP(I+1,JW) > DS(JB)) JB = JB+1
  END DO
  NBL = NBL-1

! Water surface elevation, water surface deviation, ice cover, and sediment oxygen demand

  CONV(1,:) = BLANK
  DO I=IUPR,IDPR
    DO JJB=1,NBR
      IF (ISNP(I,JW) >= US(JJB)-1 .AND. ISNP(I,JW) <= DS(JJB)+1) EXIT
    END DO
    WRITE (CONV(1,I),'(F10.3)') EL(KTWB(JW),ISNP(I,JW))-Z(ISNP(I,JW))*COSA(JJB)
  END DO
  WRITE (SNP(JW),'(/A//2X,1000I10)') '          Water Surface, m',(ISNP(I,JW),I=IUPR,IDPR)
  WRITE (SNP(JW),'(2X,1000A10/)') (CONV(1,I),I=IUPR,IDPR)
  DO I=IUPR,IDPR
    WRITE (CONV(1,I),'(F10.4)') SNGL(Z(ISNP(I,JW)))
  END DO
  WRITE (SNP(JW),'(/A//2X,1000I10)') '          Water Surface Deviation (positive downwards), m',(ISNP(I,JW),I=IUPR,IDPR)
  WRITE (SNP(JW),'(2X,1000A10/)') (CONV(1,I),I=IUPR,IDPR)
  IF (ICE_CALC(JW)) THEN
    DO I=IUPR,IDPR
      WRITE (CONV(1,I),'(F10.3)') ICETH(ISNP(I,JW))
    END DO
    WRITE (SNP(JW),'(/A//3X,1000A10)') '          Ice Thickness, m',(CONV(1,I),I=IUPR,IDPR)
  END IF
  IF (CONSTITUENTS) THEN
    DO I=IUPR,IDPR
      WRITE (CONV(1,I),'(F10.3)') SOD(ISNP(I,JW))*DAY
    END DO
    IF (OXYGEN_DEMAND) THEN
      WRITE (SNP(JW),'(/A//3X,1000A10/)') '          Sediment Oxygen Demand, g/m^2/day',(CONV(1,I),I=IUPR,IDPR)
    END IF
  END IF

! Hydrodynamic variables and temperature

  CONV = BLANK
  DO JH=1,NHY
    L = LEN_TRIM(FMTH(JH))
    IF (PRINT_HYDRO(JH,JW)) THEN
      DO I=IUPR,IDPR
        IF (JH == 1) THEN
          DO K=KTWB(JW),KB(ISNP(I,JW))
            WRITE (CONV(K,I),FMTH(JH)) INT(HYD(K,ISNP(I,JW),JH))
          END DO
        ELSE IF (JH > 6) THEN
          DO K=KTWB(JW),KB(ISNP(I,JW))
            WRITE (CONV(K,I),FMTH(JH)(1:L)) HYD(K,ISNP(I,JW),JH)*DLT
          END DO
        ELSE
          DO K=KTWB(JW),KB(ISNP(I,JW))
            WRITE (CONV(K,I),FMTH(JH)(1:L)) HYD(K,ISNP(I,JW),JH)*HMULT(JH)
          END DO
        END IF
      END DO
      IF (NEW_PAGE) THEN
        WRITE (SNP(JW),'("1",11(A/1X))') (TITLE(J),J=1,11)
        NLINES = KMX-KTWB(JW)+14
      END IF
      NLINES   = NLINES+KMX-KTWB(JW)+11
      NEW_PAGE = NLINES > 72
      WRITE (SNP(JW),'(/1X,3(A,1X,I0),A,F0.2,A/)') MONTH,GDAY,',',YEAR,'    Julian day = ',INT(JDAY),' days ',                     &
                                                  (JDAY-INT(JDAY))*24.0,' hours   '//HNAME(JH)
      WRITE (SNP(JW),'(1X,A,1000I10)') 'Layer  Depth',(ISNP(I,JW),I=IUPR,IDPR)
      DO K=KTWB(JW),KBR
        WRITE (SNP(JW),'(1X,I4,F8.2,1000A10)') K,DEPTHM(K,DS(BS(JW))),(CONV(K,I),I=IUPR,IDPR)
      END DO
    END IF
  END DO

! Constituent concentrations

  IF (CONSTITUENTS) THEN
    DO JAC=1,NAC
      JC = CN(JAC)
      L  = LEN_TRIM(FMTC(JC))
      IF (PRINT_CONST(JC,JW)) THEN
        DO I=IUPR,IDPR
          DO K=KTWB(JW),KB(ISNP(I,JW))
            WRITE (CONV(K,I),FMTC(JC)(1:L)) C2(K,ISNP(I,JW),JC)*CMULT(JC)
          END DO
        END DO
        IF (NEW_PAGE) THEN
          WRITE (SNP(JW),'("1",11(A/1X))') (TITLE(J),J=1,11)
          NLINES = KMX-KTWB(JW)+14
        END IF
        NLINES   = NLINES+KMX-KTWB(JW)+11
        NEW_PAGE = NLINES > 72
        WRITE (SNP(JW),'(/1X,3(A,1X,I0),A,F0.2,A/)') MONTH,GDAY,',',YEAR,'    Julian Date ',INT(JDAY),' days ',(JDAY-INT(JDAY))    &
                                                     *24.0,' hours   '//CNAME(JC)
        WRITE (SNP(JW),'(1X,A,1000I10)') 'Layer  Depth',(ISNP(I,JW),I=IUPR,IDPR)
        DO K=KTWB(JW),KBR
          WRITE (SNP(JW),'(1X,I4,F8.2,1000A10)') K,DEPTHM(K,DS(BS(JW))),(CONV(K,I),I=IUPR,IDPR)
        END DO
      END IF
    END DO

!** Derived constituent concentrations

    DO JD=1,NDC
      L = LEN_TRIM(FMTCD(JD))
      IF (PRINT_DERIVED(JD,JW)) THEN
        DO I=IUPR,IDPR
          DO K=KTWB(JW),KB(ISNP(I,JW))
            WRITE (CONV(K,I),FMTCD(JD)(1:L)) CD(K,ISNP(I,JW),JD)*CDMULT(JD)
          END DO
        END DO
        IF (NEW_PAGE) THEN
          WRITE (SNP(JW),'("1",11(A/1X))') (TITLE(J),J=1,11)
          NLINES = KMX-KTWB(JW)+14
        END IF
        NLINES   = NLINES+KMX-KTWB(JW)+11
        NEW_PAGE = NLINES > 72
        WRITE (SNP(JW),'(/1X,3(A,1X,I0),A,F0.2,A/)') MONTH,GDAY,',',YEAR,'    Julian Date ',INT(JDAY),' days ',(JDAY-INT(JDAY))    &
                                                     *24.0,' hours    '//CDNAME(JD)
        WRITE (SNP(JW),'(1X,A,1000I10)') 'Layer  Depth',(ISNP(I,JW),I=IUPR,IDPR)
        DO K=KTWB(JW),KBR
          WRITE (SNP(JW),'(1X,I4,F8.2,1000A10)') K,DEPTHM(K,DS(BS(JW))),(CONV(K,I),I=IUPR,IDPR)
        END DO
      END IF
    END DO

!** Sediment

    IF (PRINT_SEDIMENT(JW)) THEN
      DO I=IUPR,IDPR
        DO K=KTWB(JW),KB(ISNP(I,JW))
          WRITE (CONV(K,I),'(F10.2)') SED(K,ISNP(I,JW))
        END DO
      END DO
      IF (NEW_PAGE) THEN
        WRITE (SNP(JW),'("1",11(A/1X))') (TITLE(J),J=1,11)
        NLINES = KMX-KTWB(JW)+14
      END IF
      NLINES   = NLINES+KMX-KTWB(JW)+11
      NEW_PAGE = NLINES > 72
      WRITE (SNP(JW),'(/1X,3(A,1X,I0),A,F0.2,A/)') MONTH,GDAY,',',YEAR,'    Julian Date ',INT(JDAY),' days ',(JDAY-INT(JDAY))*24.0,&
                                                 ' hours     Organic sediments, g/m^3'
      WRITE (SNP(JW),'(1X,A,1000I10)') 'Layer  Depth',(ISNP(I,JW),I=IUPR,IDPR)
      DO K=KTWB(JW),KBR
        WRITE (SNP(JW),'(1X,I4,F8.2,1000A10)') K,DEPTHM(K,DS(BS(JW))),(CONV(K,I),I=IUPR,IDPR)
      END DO
    END IF

    IF (PRINT_SEDIMENT(JW)) THEN
      DO I=IUPR,IDPR
        DO K=KTWB(JW),KB(ISNP(I,JW))
          WRITE (CONV(K,I),'(F10.2)') SEDp(K,ISNP(I,JW))
        END DO
      END DO
      IF (NEW_PAGE) THEN
        WRITE (SNP(JW),'("1",11(A/1X))') (TITLE(J),J=1,11)
        NLINES = KMX-KTWB(JW)+14
      END IF
      NLINES   = NLINES+KMX-KTWB(JW)+11
      NEW_PAGE = NLINES > 72
      WRITE (SNP(JW),'(/1X,3(A,1X,I0),A,F0.2,A/)') MONTH,GDAY,',',YEAR,'    Julian Date ',INT(JDAY),' days ',(JDAY-INT(JDAY))*24.0,&
                                                 ' hours     Organic phosphorus sediments, g/m^3'
      WRITE (SNP(JW),'(1X,A,1000I10)') 'Layer  Depth',(ISNP(I,JW),I=IUPR,IDPR)
      DO K=KTWB(JW),KBR
        WRITE (SNP(JW),'(1X,I4,F8.2,1000A10)') K,DEPTHM(K,DS(BS(JW))),(CONV(K,I),I=IUPR,IDPR)
      END DO
    END IF

    IF (PRINT_SEDIMENT(JW)) THEN
      DO I=IUPR,IDPR
        DO K=KTWB(JW),KB(ISNP(I,JW))
          WRITE (CONV(K,I),'(F10.2)') SEDn(K,ISNP(I,JW))
        END DO
      END DO
      IF (NEW_PAGE) THEN
        WRITE (SNP(JW),'("1",11(A/1X))') (TITLE(J),J=1,11)
        NLINES = KMX-KTWB(JW)+14
      END IF
      NLINES   = NLINES+KMX-KTWB(JW)+11
      NEW_PAGE = NLINES > 72
      WRITE (SNP(JW),'(/1X,3(A,1X,I0),A,F0.2,A/)') MONTH,GDAY,',',YEAR,'    Julian Date ',INT(JDAY),' days ',(JDAY-INT(JDAY))*24.0,&
                                                 ' hours     Organic nitrogen sediments, g/m^3'
      WRITE (SNP(JW),'(1X,A,1000I10)') 'Layer  Depth',(ISNP(I,JW),I=IUPR,IDPR)
      DO K=KTWB(JW),KBR
        WRITE (SNP(JW),'(1X,I4,F8.2,1000A10)') K,DEPTHM(K,DS(BS(JW))),(CONV(K,I),I=IUPR,IDPR)
      END DO
    END IF

    IF (PRINT_SEDIMENT(JW)) THEN
      DO I=IUPR,IDPR
        DO K=KTWB(JW),KB(ISNP(I,JW))
          WRITE (CONV(K,I),'(F10.2)') SEDc(K,ISNP(I,JW))
        END DO
      END DO
      IF (NEW_PAGE) THEN
        WRITE (SNP(JW),'("1",11(A/1X))') (TITLE(J),J=1,11)
        NLINES = KMX-KTWB(JW)+14
      END IF
      NLINES   = NLINES+KMX-KTWB(JW)+11
      NEW_PAGE = NLINES > 72
      WRITE (SNP(JW),'(/1X,3(A,1X,I0),A,F0.2,A/)') MONTH,GDAY,',',YEAR,'    Julian Date ',INT(JDAY),' days ',(JDAY-INT(JDAY))*24.0,&
                                                 ' hours     Organic carbon sediments, g/m^3'
      WRITE (SNP(JW),'(1X,A,1000I10)') 'Layer  Depth',(ISNP(I,JW),I=IUPR,IDPR)
      DO K=KTWB(JW),KBR
        WRITE (SNP(JW),'(1X,I4,F8.2,1000A10)') K,DEPTHM(K,DS(BS(JW))),(CONV(K,I),I=IUPR,IDPR)
      END DO
    END IF

! Amaila Start
    IF (PRINT_SEDIMENT1(JW)) THEN
      DO I=IUPR,IDPR
        DO K=KTWB(JW),KB(ISNP(I,JW))
          WRITE (CONV(K,I),'(F10.2)') SED1(K,ISNP(I,JW))
        END DO
      END DO
      IF (NEW_PAGE) THEN
        WRITE (SNP(JW),'("1",11(A/1X))') (TITLE(J),J=1,11)
        NLINES = KMX-KTWB(JW)+14
      END IF
      NLINES   = NLINES+KMX-KTWB(JW)+11
      NEW_PAGE = NLINES > 72
      WRITE (SNP(JW),'(/1X,3(A,1X,I0),A,F0.2,A/)') MONTH,GDAY,',',YEAR,'    Julian Date ',INT(JDAY),' days ',(JDAY-INT(JDAY))*24.0,&
                                                 ' hours     Tree Organic Matter-Labile, g/m^3'
      WRITE (SNP(JW),'(1X,A,1000I10)') 'Layer  Depth',(ISNP(I,JW),I=IUPR,IDPR)
      DO K=KTWB(JW),KBR
        WRITE (SNP(JW),'(1X,I4,F8.2,1000A10)') K,DEPTHM(K,DS(BS(JW))),(CONV(K,I),I=IUPR,IDPR)
      END DO
    END IF

    IF (PRINT_SEDIMENT2(JW)) THEN
      DO I=IUPR,IDPR
        DO K=KTWB(JW),KB(ISNP(I,JW))
          WRITE (CONV(K,I),'(F10.2)') SED2(K,ISNP(I,JW))
        END DO
      END DO
      IF (NEW_PAGE) THEN
        WRITE (SNP(JW),'("1",11(A/1X))') (TITLE(J),J=1,11)
        NLINES = KMX-KTWB(JW)+14
      END IF
      NLINES   = NLINES+KMX-KTWB(JW)+11
      NEW_PAGE = NLINES > 72
      WRITE (SNP(JW),'(/1X,3(A,1X,I0),A,F0.2,A/)') MONTH,GDAY,',',YEAR,'    Julian Date ',INT(JDAY),' days ',(JDAY-INT(JDAY))*24.0,&
                                                 ' hours     Tree Organic Matter-Refactory, g/m^3'
      WRITE (SNP(JW),'(1X,A,1000I10)') 'Layer  Depth',(ISNP(I,JW),I=IUPR,IDPR)
      DO K=KTWB(JW),KBR
        WRITE (SNP(JW),'(1X,I4,F8.2,1000A10)') K,DEPTHM(K,DS(BS(JW))),(CONV(K,I),I=IUPR,IDPR)
      END DO
    END IF

! Amaila End

!** Epiphyton

    DO JE=1,NEP
      IF (PRINT_EPIPHYTON(JW,JE)) THEN
        DO I=IUPR,IDPR
          DO K=KTWB(JW),KB(ISNP(I,JW))
            WRITE (CONV(K,I),'(F10.2)') EPD(K,ISNP(I,JW),JE)
          END DO
        END DO
        IF (NEW_PAGE) THEN
          WRITE (SNP(JW),'("1",11(A/1X))') (TITLE(J),J=1,11)
          NLINES = KMX-KTWB(JW)+14
        END IF
        NLINES   = NLINES+KMX-KTWB(JW)+11
        NEW_PAGE = NLINES > 72
        WRITE (SNP(JW),'(/1X,3(A,1X,I0),A,F0.2,A/)') MONTH,GDAY,',',YEAR,'    Julian Date ',INT(JDAY),' days ',                    &
                                                    (JDAY-INT(JDAY))*24.0,' hours     Epiphyton, g/m^2'
        WRITE (SNP(JW),'(1X,A,1000I10)') 'Layer  Depth',(ISNP(I,JW),I=IUPR,IDPR)
        DO K=KTWB(JW),KBR
          WRITE (SNP(JW),'(1X,I4,F8.2,1000A10)') K,DEPTHM(K,DS(BS(JW))),(ADJUSTR(CONV(K,I)),I=IUPR,IDPR)
        END DO
      END IF
    END DO

!********* macrophytes
    do L=1,nmc
      CONV = BLANK
      IF (PRINT_macrophyte(jw,L)) THEN
        DO I=IUPR,IDPR
          DO K=ktwb(jw),KB(ISNP(I,jw))
            WRITE (CONV(K,I),'(F10.2)') mac(K,ISNP(I,jw),L)
          END DO
        END DO
        IF (NEW_PAGE) THEN
          WRITE (SNP(JW),'("1",11(A/1X))') (TITLE(J),J=1,11)
          NLINES = KMX-KTWB(JW)+14
        END IF
        NLINES   = NLINES+KMX-KTWB(JW)+11
        NEW_PAGE = NLINES > 72
        WRITE (SNP(jw),'(/1X,3(A,1X,I0),A,F0.2,2A,i0,A/)')  MONTH,GDAY,',',YEAR,'    Julian Date ',INT(JDAY),' days ',(JDAY-INT(JDAY))    &    ! cb 12/20/11
                                                     *24.0,' hours   ','  Macrophyte Group #',L,' g/m^3'
        WRITE (SNP(JW),'(1X,A,1000I10)') 'Layer  Depth',(ISNP(I,JW),I=IUPR,IDPR)
        DO K=ktwb(jw),KBR
          WRITE (SNP(JW),'(1X,I4,F8.2,1000A10)') K,DEPTHM(K,DS(BS(JW))),(CONV(K,I),I=IUPR,IDPR)
        END DO
        do i=iupr,idpr
          conv2=blank
          DO K=ktwb(jw),KB(ISNP(I,jw))
            if(k.eq.ktwb(jw))then
              jt=kti(isnp(i,jw))
            else
              jt=k
            end if
            je=kb(isnp(i,jw))
            DO jj=jt,je
              WRITE (CONV2(K,jj),'(F10.2)') macrc(jj,K,ISNP(I,jw),L)
            END DO
          END DO
          IF (NEW_PAGE) THEN
            WRITE (SNP(JW),'("1",11(A/1X))') (TITLE(J),J=1,11)
            NLINES = KMX-KTWB(JW)+14
          END IF
          NLINES   = NLINES+KMX-KTWB(JW)+11
          NEW_PAGE = NLINES > 72
!          WRITE (SNP(jw),'(/1X,3(A,1X,I0),A,F0.2,A/)')  MONTH,GDAY,',',YEAR,'    Julian Date ',INT(JDAY),' days ',(JDAY-INT(JDAY))    &
!                                                     *24.0,' hours   ','  Macrophyte Columns g/m^3'
          WRITE (SNP(JW),'(/1X,3(A,1X,I0),A,F0.2,2A,i0,A/)')  MONTH,GDAY,',',YEAR,'    Julian Date ',INT(JDAY),' days ',(JDAY-INT(JDAY))    &  ! cb 12/20/11
                                                     *24.0,' hours   ','  Macrophyte Group #',L,' Columns g/m^3'
          WRITE (SNP(JW),3052)ISNP(I,JW)
3052      FORMAT(7X,"SEGMENT   ",I8)
          JT=KTI(ISNP(I,JW))
          JE=KB(ISNP(I,JW))
          WRITE (SNP(JW),'(A,200I10)') ' LAYER  DEPTH',(JJ,JJ=JT,JE)
          DO K=KTWB(JW),KB(ISNP(I,JW))
            WRITE (SNP(JW),'(1X,I4,F8.2,1000A10)') K,DEPTHM(K,DS(BS(JW))),(CONV2(K,JJ),JJ=JT,JE)
          END DO
          WRITE (SNP(JW),3053)L,ISNP(I,JW)
3053      FORMAT(7X,"MACROPHYTE GROUP #",I0," LIMITATION SEGMENT   ",I8)   ! CB 12/20/11
          WRITE (SNP(JW),'(A,200I10)')' LAYER  DEPTH', (JJ,JJ=JT,JE)
          MLFPR=BLANK
          DO K=KTWB(JW),KB(ISNP(I,JW))
            IF(K.EQ.KTWB(JW))THEN
              JT=KTI(ISNP(I,JW))
            ELSE
              JT=K
            END IF
            DO JJ=JT,JE
              LIMIT = MIN(MPLIM(K,ISNP(I,JW),L),MNLIM(K,ISNP(I,JW),L),MCLIM(K,ISNP(I,JW),L),MLLIM(JJ,K,ISNP(I,JW),L))
              IF (LIMIT == MPLIM(K,ISNP(I,JW),L)) THEN
                WRITE (LFAC,'(F8.4)') MPLIM(K,ISNP(I,JW),L)
                MLFPR(JJ,K,ISNP(I,JW),L) = ' P'//LFAC
              ELSE IF (LIMIT == MNLIM(K,I,L)) THEN
                WRITE (LFAC,'(F8.4)') MNLIM(K,ISNP(I,JW),L)
                MLFPR(JJ,K,ISNP(I,JW),L) = ' N'//LFAC
              ELSE IF (LIMIT == MCLIM(K,I,L)) THEN
                WRITE (LFAC,'(F8.4)') MCLIM(K,ISNP(I,JW),L)
                MLFPR(JJ,K,ISNP(I,JW),L) = ' C'//LFAC
              ELSE IF (LIMIT == MLLIM(JJ,K,ISNP(I,JW),L)) THEN
                WRITE (LFAC,'(F8.4)') MLLIM(JJ,K,ISNP(I,JW),L)
                MLFPR(JJ,K,ISNP(I,JW),L) = ' L'//LFAC
              END IF
            END DO
          END DO
          JT=KTI(ISNP(I,JW))
          JE=KB(ISNP(I,JW))
          DO K=KTWB(JW),KB(ISNP(I,JW))
            WRITE (SNP(JW),'(1X,I4,F8.2,1000A10)') K,DEPTHM(K,DS(BS(JW))),(MLFPR(JJ,K,ISNP(I,JW),L),JJ=JT,JE)
          END DO
        END DO
      END IF
    END DO

!** Algal nutrient limitations

    DO JA=1,NAL
      IF (LIMITING_FACTOR(JA)) THEN
        IF (NEW_PAGE) THEN
          WRITE (SNP(JW),'("1",11(A/1X))') (TITLE(J),J=1,11)
          NLINES = KMX-KTWB(JW)+14
        END IF
        NLINES   = NLINES+KMX-KTWB(JW)+11
        NEW_PAGE = NLINES > 72
        DO I=IUPR,IDPR     !mlm 6/30/06
          DO K=KTWB(JW),KB(ISNP(I,JW))   !mlm  6/30/06
            LIMIT = MIN(APLIM(K,ISNP(I,JW),JA),ANLIM(K,ISNP(I,JW),JA),ASLIM(K,ISNP(I,JW),JA),ALLIM(K,ISNP(I,JW),JA))
            IF (LIMIT == APLIM(K,ISNP(I,JW),JA)) THEN
              WRITE (LFAC,'(F8.4)') APLIM(K,ISNP(I,JW),JA)
              LFPR(K,ISNP(I,JW)) = ' P'//LFAC
            ELSE IF (LIMIT == ANLIM(K,ISNP(I,JW),JA)) THEN
              WRITE (LFAC,'(F8.4)') ANLIM(K,ISNP(I,JW),JA)
              LFPR(K,ISNP(I,JW)) = ' N'//LFAC
            ELSE IF (LIMIT == ASLIM(K,ISNP(I,JW),JA)) THEN
              WRITE (LFAC,'(F8.4)') ASLIM(K,ISNP(I,JW),JA)
              LFPR(K,ISNP(I,JW)) = ' S'//LFAC
            ELSE IF (LIMIT == ALLIM(K,ISNP(I,JW),JA)) THEN
              WRITE (LFAC,'(F8.4)') ALLIM(K,ISNP(I,JW),JA)
              LFPR(K,ISNP(I,JW)) = ' L'//LFAC
            END IF
          END DO
        END DO
        WRITE (SNP(JW),'(/1X,3(A,1X,I0),A,F0.2,A,I0,A/)') MONTH,GDAY,',',YEAR,'    Julian Date',INT(JDAY),' days ',                &
                                                         (JDAY-INT(JDAY))*24.0,' hours    Algal group ',JA,' limiting factor'
        WRITE (SNP(JW),'(1X,A,1000I10)') 'Layer  Depth',(ISNP(I,JW),I=IUPR,IDPR)
        DO K=KTWB(JW),KBR
          WRITE (SNP(JW),'(1X,I4,F8.2,1000A)') K,DEPTHM(K,DS(BS(JW))),(LFPR(K,ISNP(I,JW)),I=IUPR,IDPR)
        END DO
      END IF
    END DO

!** Epiphyton nutrient limitations

    DO JE=1,NEP
      IF (PRINT_EPIPHYTON(JW,JE)) THEN
        IF (NEW_PAGE) THEN
          WRITE (SNP(JW),'("1",11(A/1X))') (TITLE(J),J=1,11)
          NLINES = KMX-KTWB(JW)+14
        END IF
        NLINES   = NLINES+KMX-KTWB(JW)+11
        NEW_PAGE = NLINES > 72
        DO I=IUPR,IDPR                                            !mlm   6/30/2006
          DO K=KTWB(JW),KB(ISNP(I,JW))                            !mlm  6/30/2006
            LIMIT = MIN(EPLIM(K,ISNP(I,JW),JE),ENLIM(K,ISNP(I,JW),JE),ESLIM(K,ISNP(I,JW),JE),ELLIM(K,ISNP(I,JW),JE))
            IF (LIMIT == EPLIM(K,ISNP(I,JW),JE)) THEN
              WRITE (LFAC,'(F8.4)') EPLIM(K,ISNP(I,JW),JE)
              LFPR(K,ISNP(I,JW)) = ' P'//LFAC
            ELSE IF (LIMIT == ENLIM(K,ISNP(I,JW),JE)) THEN
              WRITE (LFAC,'(F8.4)') ENLIM(K,ISNP(I,JW),JE)
              LFPR(K,ISNP(I,JW)) = ' N'//LFAC
            ELSE IF (LIMIT == ESLIM(K,ISNP(I,JW),JE)) THEN
              WRITE (LFAC,'(F8.4)') ESLIM(K,ISNP(I,JW),JE)
              LFPR(K,ISNP(I,JW)) = ' S'//LFAC
            ELSE IF (LIMIT == ELLIM(K,ISNP(I,JW),JE)) THEN
              WRITE (LFAC,'(F8.4)') ELLIM(K,ISNP(I,JW),JE)
              LFPR(K,ISNP(I,JW)) = ' L'//LFAC
            END IF
          END DO
        END DO
        WRITE (SNP(JW),'(/1X,3(A,1X,I0),A,F0.2,A,I0,A/)') MONTH,GDAY,',',YEAR,'    Julian Date',INT(JDAY),' days ',                &
                                                         (JDAY-INT(JDAY))*24.0,' hours    Epiphyton group ',JE,' limiting factor'
        WRITE (SNP(JW),'(1X,A,1000I10)') 'Layer  Depth',(ISNP(I,JW),I=IUPR,IDPR)
        DO K=KTWB(JW),KBR
          WRITE (SNP(JW),'(1X,I4,F8.2,1000A10)') K,DEPTHM(K,DS(BS(JW))),(LFPR(K,ISNP(I,JW)),I=IUPR,IDPR)
        END DO
      END IF
    END DO
  END IF
END SUBROUTINE OUTPUT
