!***********************************************************************************************************************************
!**                                                  T D G   T A R G E T    I N I T                                               **
!***********************************************************************************************************************************
! Input and output files changed into the csv format, 7/2021
Subroutine InitTDGtarget
  Use Selective1TDGtarget; USE MAIN; USE modSYSTDG, ONLY: POWNO, FLNO, NBAY, BEGNO, ENDNO, POWGTNO, FLGTNO, NRO, TDGLOC, GTNAME
  USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  !
  IMPLICIT NONE
  INTEGER       :: it, n, ig
  CHARACTER(8)  :: AID1
  CHARACTER(72) :: TGAFN
  
  targetfnno = 8888
  open (targetfnno+1, file='w2_TDGtarget.csv', status='old')
  !read (targetfnno+1, '(///(8X,A72))') (TITLETDGTARGET(it), it=1,10)
  !read (targetfnno+1,'(//8x,2f8.3)') tsfreq, tsconv
  !read (targetfnno+1,'(//8x,a8,3f8.3,a8,i8,a8,i8)')tsyearly, tstsrt, tstend, tstarget, tsdynsel, tsiteration, dyupdate, dygroup   
  read (targetfnno+1,*) 
  read (targetfnno+1,*) 
  read (targetfnno+1,*)
  do it = 1,10
    read (targetfnno+1,*) TITLETDGTARGET(it)
  end do
  read (targetfnno+1,*) 
  read (targetfnno+1,*)
  read (targetfnno+1,*) AID1, tsfreq, tsconv     
  read (targetfnno+1,*) 
  read (targetfnno+1,*)
  read (targetfnno+1,*) AID1, tsyearly, tstsrt, tstend, tstarget, tsdynsel, tsiteration, dyupdate, dygroup 
  
  tsyearly=ADJUSTR(tsyearly); tsdynsel=ADJUSTR(tsdynsel); dyupdate=ADJUSTR(dyupdate)
  if (tsyearly=='     OFF' .AND. tstsrt<TMSTRT) tstsrt = TMSTRT
  NXTSPLIT  = TMSTRT
  NXTSPLIT2 = TMSTRT
  if (tsyearly=='     OFF') then
    DAYTEST = jday
  else
    DAYTEST = real(jdayg) + jday - int(jday)
  end if
  IF (NXTSPLIT>TMSTRT)  NXTSPLIT  = TMSTRT
  IF (DAYTEST<=tstsrt)  DAYTEST   = tstsrt
  IF (tstsrt>TMSTRT)    NXTSPLIT  = tstsrt
  IF (NXTSPLIT2>TMSTRT) NXTSPLIT2 = TMSTRT
  IF (tstsrt>TMSTRT)    NXTSPLIT2 = tstsrt
  NGSP  = NBAY+NRO
  NGPH  = POWNO
  NGFL  = FLNO 
  NOUTS = NGT
  NGSPPH = NGSP+NGPH
  allocate(SPGTNO(NGSP), SPPRIOR(NGSP), SPMINFRAC(NGSP))
  if (NGPH>0) allocate(PHGTNO(NGPH), PHMAXFLOW(NGPH))
  it = 0
  do ig = 1,NGT
    if (GTNAME(ig)) then
      it = it+1
      SPGTNO(it) = ig
    end if
  end do
  
  if (NGPH>0) then
    it = 0
    do ig = 1,NGT
      if (GTTYP(ig)=='     POW') then
        it = it+1
        PHGTNO(it) = ig
      end if
    end do
  end if
  
  !read (targetfnno+1,'(//8x,<NGSP>i8)') (SPPRIOR(n),n=1,NGSP)  
  !read (targetfnno+1,'(//8x,<NGSP>f8.3)') (SPMINFRAC(n),n=1,NGSP)
  !if (NGPH>0) read (targetfnno+1,'(//8x,<NGPH>f8.3)') (PHMAXFLOW(n),n=1,NGPH) 
  read (targetfnno+1,*) 
  read (targetfnno+1,*)
  read (targetfnno+1,*) AID1, (SPPRIOR(n),n=1,NGSP) 
  read (targetfnno+1,*)
  read (targetfnno+1,*)
  read (targetfnno+1,*) AID1, (SPMINFRAC(n),n=1,NGSP) 
  if (NGPH>0) then
    read (targetfnno+1,*)
    read (targetfnno+1,*)
    read (targetfnno+1,*) AID1, (PHMAXFLOW(n),n=1,NGPH)
  end if

  do n=1,NGSP
    if (SPMINFRAC(n)>1.0) SPMINFRAC(n) = 1.0    ! remove unrealistic input value
  end do
  do n=1,NGPH    
    if (PHMAXFLOW(n)<0.0) PHMAXFLOW(n) = 0.0    ! remove unrealistic input value
  end do
  if (tsconv<0.001) tsconv = 0.001  ! constrain the convergence criterion to be >= 0.001 and <= 5.0
  if (tsconv> 5.0)  tsconv = 5.0
  
  ! OPEN DYNAMIC TDG TARGET FILES    
  if (tsdynsel=='      ON') then
    !read (targetfnno+1,'(//(8X,A72))') TGAFN
    read (targetfnno+1,*)
    read (targetfnno+1,*)
    read (targetfnno+1,*) AID1, TGAFN
    open (targetfnno+3, file=TGAFN, status='old')
    read (targetfnno+3,*)
    read (targetfnno+3,*)
    read (targetfnno+3,*)
    read (targetfnno+3,*) nxtjday, tstarget2         
    tstarget = tstarget2                                         
    read (targetfnno+3,*) nxtjday, tstarget2              
  end if
  close(targetfnno+1)
  
  !!  Initial output file  
  open (targetfnno, FILE='TDGTarget_output.csv', status='unknown')
  !if (NGT>0) write (targetfnno,'(3A,<NGT>(A,i2))')'     JDAY','         TDG' ,'      SUM Q   ',('      Q',n, n=1,NGT)
  if (NGT>0) write (targetfnno,'("     JDAY,", " C,", "       TDG,", "    SUM Q,", <NGT>(A,i2,","))') ('      Q',n, n=1,NGT)
  open (targetfnno+2, FILE='TDGTarget_warning.opt', action="READWRITE", status='unknown')
  return
End subroutine InitTDGtarget

!***********************************************************************************************************************************
!**                                                 T D G    T A R G E T   S U B R O U T I N E                                    **
!***********************************************************************************************************************************
Subroutine TDGtarget
  Use Selective1TDGtarget; USE modSYSTDG, ONLY : TDG_TDG, TDGLOC, SYSTDG_TDG; USE MAIN, ONLY:targetfnno,WARNING_OPEN
  USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART; 
  ! 
  IMPLICIT NONE
  !
  REAL,    ALLOCATABLE, DIMENSION(:)  :: QGTSAVE
  LOGICAL, ALLOCATABLE, DIMENSION(:)  :: SP_ACTIVE, PH_ACTIVE
  INTEGER                             :: ig, ii, ITERATION, prior_top, priortop_n
  REAL                                :: Q_ALL, Q_SP, Q_PH, QSP_AVL, QPH_AVL, Q_TEMP
  REAL                                :: SUM_SP_FRAC, SUM_PH_MAXFLOW
  REAL                                :: Q_MAX, Q_MIN, Q_CUT, Q_CUTTED, Q_LEFT, QPH_ADDED
  REAL                                :: SUM_TOP_FLOW, SUM_TOP_MINFRAC
  REAL                                :: SUM_QGT, SUM_QGT2
  REAL                                :: MINV=0.0000000001
  INTEGER, ALLOCATABLE, DIMENSION(:)  :: priortop_spno
  CHARACTER(1)                        :: CO

  ALLOCATE(QGTSAVE(NGT), SP_ACTIVE(NGSP))
  IF (NGPH>0) ALLOCATE(PH_ACTIVE(NGPH))
  QGTSAVE = QGT   
  SUM_QGT = 0.0
  DO ig = 1,NGT
    SUM_QGT = SUM_QGT + QGT(ig)
  END DO
  SP_ACTIVE =.FALSE.
  PH_ACTIVE =.FALSE.
  !  
  SUM_SP_FRAC = 0.0
  DO ig = 1,NGSP
    IF (QGT(SPGTNO(ig))>0.0) SUM_SP_FRAC = SUM_SP_FRAC + SPMINFRAC(ig)
  END DO
  IF (SUM_SP_FRAC>1.0) THEN   ! IF SUM SP MIN FRAC IS BIGGER THAN 1.0 THEN ADJUST THE SPMINFRAC BY PERCENTAGE
    DO ig = 1,NGSP
      IF (QGT(SPGTNO(ig))>0.0) SPMINFRAC(ig) = SPMINFRAC(ig) * (1.0/SUM_SP_FRAC)
    END DO
    SUM_SP_FRAC = 1.0
  END IF
  !  
  SUM_PH_MAXFLOW = 0.0
  DO ig = 1,NGPH
    IF (QGT(PHGTNO(ig))>0.0) SUM_PH_MAXFLOW = SUM_PH_MAXFLOW + PHMAXFLOW(ig)
  END DO  
  !  
  Q_SP = 0.0
  DO ig = 1,NGSP
    Q_SP = Q_SP + QGT(SPGTNO(ig))
  END DO
  
  IF (Q_SP>0.0) THEN     
    IF (tsyearly=='     OFF') THEN
      DAYTEST = JDAY
    ELSE
      DAYTEST = real(JDAYG) + JDAY - int(JDAY)
    END IF
    
    CALL SYSTDG_TDG           ! CALCULATE CURRENT TDG
    IF (tsdynsel=='      ON') THEN
      IF (DAYTEST>tstsrt) tstarget = tstarget2 
      DO WHILE (JDAY>=nxtjday)
        READ (targetfnno+3,*) nxtjday, tstarget2 
        tstarget = tstarget2 
      END DO
    END IF
    
    CO = ' '
    IF (DAYTEST>=tstsrt .AND. DAYTEST<=tstend .AND. TDG_TDG>tstarget+tsconv) THEN
      ! INITIAL VARIABLES FOR ITERATIONS
      CO = 'R'
      Q_SP  = 0.0
      Q_PH  = 0.0
      Q_ALL = 0.0
      DO ig = 1,NGT
        Q_ALL = Q_ALL + QGT(ig)
      END DO
      DO ig = 1,NGSP
        Q_SP = Q_SP + QGT(SPGTNO(ig))
        IF (QGT(SPGTNO(ig))>0.0 .AND. (QGT(SPGTNO(ig))>Q_ALL*SPMINFRAC(SPGTNO(ig))+MINV)) SP_ACTIVE(ig) = .TRUE.
      END DO
      DO ig = 1,NGPH
        Q_PH = Q_PH + QGT(PHGTNO(ig))
        IF (QGT(PHGTNO(ig))>0.0 .AND. (QGT(SPGTNO(ig))+MINV<PHMAXFLOW(ig))) PH_ACTIVE(ig) = .TRUE.
      END DO
          
      QSP_AVL = Q_SP - Q_ALL*SUM_SP_FRAC
      IF (SUM_PH_MAXFLOW>0.0) QPH_AVL = SUM_PH_MAXFLOW - Q_PH
      IF (QPH_AVL<0.0) QPH_AVL = 0.0
      Q_MAX = Q_SP
      IF (TDGLOC=='     REL' .AND. SUM_PH_MAXFLOW/=0.0 .AND. QPH_AVL>0.0) Q_MAX = MIN(Q_SP, QPH_AVL)   ! IF RELEASE TDG, FLOW CUT TO POWERHOUSE MUST LESS THAN MAX FLOW
      Q_MIN = MAX(0.0, Q_ALL*SUM_SP_FRAC)
      IF (TDGLOC=='     REL' .AND. SUM_PH_MAXFLOW/=0.0 ) Q_MIN = MAX(Q_ALL*SUM_SP_FRAC, Q_SP-QPH_AVL, 0.0) 
      Q_CUT = 0.5*(Q_MAX+Q_MIN)
      ITERATION = 1
      IF (JDAY>=NXTSPLIT2 .and. dyupdate=='      ON') CALL Dy_Priority            ! DYNAMIC PRIORITY UPDATE
      
      DO WHILE (ABS(TDG_TDG-tstarget)>tsconv .AND. ITERATION<=tsiteration .AND. Q_CUT>0.0)
        Q_CUTTED = 0.0 
        Q_LEFT = Q_SP-Q_CUTTED
        DO ig = 1, NGSP
          IF (QGT(SPGTNO(ig))>0.0 .AND. (QGT(SPGTNO(ig))>Q_ALL*SPMINFRAC(SPGTNO(ig))+MINV)) SP_ACTIVE(ig) = .TRUE.
        END DO
        DO ig = 1,NGPH
          IF (QGT(PHGTNO(ig))>0.0 .AND. (QGT(SPGTNO(ig))+MINV<PHMAXFLOW(ig))) PH_ACTIVE(ig) = .TRUE.
        END DO
        Q_TEMP = 0.0
        
        DO WHILE (Q_CUTTED+MINV<Q_CUT .AND. (Q_LEFT>Q_ALL*SUM_SP_FRAC+MINV))  ! CUT FLOW TO SP UNTILL Q_CUTTED = Q_CUT
          ! UPDATE PRIOR TOP
          prior_top = -999                  ! highest prior
          DO ig = 1,NGSP
            IF (prior_top==-999 .OR. SPPRIOR(ig)<prior_top ) THEN
              IF (SP_ACTIVE(ig) .AND. SPPRIOR(ig)>0) prior_top = SPPRIOR(ig)
            END IF
          END DO
          priortop_n = 0
          DO ig = 1,ngsp                   ! initial highest prior out NO
            IF (SPPRIOR(ig)==prior_top .AND. SP_ACTIVE(ig)) priortop_n = priortop_n + 1
          END DO
          IF (priortop_n>0) ALLOCATE(priortop_spno(priortop_n))
          ii=0
          SUM_TOP_FLOW = 0.0
          SUM_TOP_MINFRAC = 0.0
          DO ig = 1,NGSP                 ! initial highest prior GTNO, SUM FLOW AND SUM MINFRAC
            IF (SPPRIOR(ig)==prior_top .AND. SP_ACTIVE(ig)) THEN
              ii = ii+1
              priortop_spno(ii) = ig
              SUM_TOP_FLOW = SUM_TOP_FLOW + QGT(SPGTNO(ig))
              SUM_TOP_MINFRAC = SUM_TOP_MINFRAC + SPMINFRAC(ig)
            END IF
          END DO
          
          IF ((SUM_TOP_FLOW-Q_ALL*SUM_TOP_MINFRAC)<=Q_CUT-Q_CUTTED) THEN  ! top prior flow is not enough to cut all Q_CUT
            DO ig = 1,priortop_n
              Q_CUTTED = Q_CUTTED + QGT(SPGTNO(priortop_spno(ig))) - Q_ALL*SPMINFRAC(priortop_spno(ig))
              QGT(SPGTNO(priortop_spno(ig))) = Q_ALL*SPMINFRAC(priortop_spno(ig))
              SP_ACTIVE(priortop_spno(ig)) = .FALSE.
            END DO
            Q_LEFT = Q_SP-Q_CUTTED
          ELSE 
            Q_TEMP = Q_CUT-Q_CUTTED       ! top prior available flow if more than Q_CUT-Q_CUTTED, so cut it by flow percentage
            DO ig = 1,priortop_n
              Q_CUTTED = Q_CUTTED + Q_TEMP*(QGT(SPGTNO(priortop_spno(ig))) - Q_ALL*SPMINFRAC(priortop_spno(ig)))/(SUM_TOP_FLOW-Q_ALL*SUM_TOP_MINFRAC)
              QGT(SPGTNO(priortop_spno(ig))) = QGT(SPGTNO(priortop_spno(ig))) - Q_TEMP*(QGT(SPGTNO(priortop_spno(ig))) &
                                             - Q_ALL*SPMINFRAC(priortop_spno(ig)))/(SUM_TOP_FLOW-Q_ALL*SUM_TOP_MINFRAC)
            END DO
            Q_LEFT = Q_SP-Q_CUTTED
          END IF
          
          DO ig = 1,priortop_n
            IF (QGT(SPGTNO(priortop_spno(ig)))>(Q_ALL*SPMINFRAC(priortop_spno(ig))+MINV)) SP_ACTIVE(priortop_spno(ig)) = .TRUE.
          END DO
          IF (priortop_n>0) DEALLOCATE(priortop_spno)
          ! THIS PRIOR TOP ENDED
        END DO  ! END DO WHILE LOOP FOR Q_CUT
        
        IF (Q_LEFT>Q_SP-Q_CUT) THEN       ! KEEP FLOW BALANCE
          DO ig = 1,NGSP
            IF (QGT(SPGTNO(ig))>(Q_ALL*SPMINFRAC(ig)+(Q_LEFT+Q_CUT-Q_SP)) .AND. SP_ACTIVE(ig)) THEN
              QGT(SPGTNO(ig)) = QGT(SPGTNO(ig))-(Q_LEFT+Q_CUT-Q_SP)
              Q_LEFT = Q_SP - Q_CUT
              Q_CUTTED = Q_CUT
              EXIT
            END IF
          END DO
        END IF 
        
        ! CUT FLOW TO POWERHOUSE
        QPH_ADDED = 0.0
        IF (SUM_PH_MAXFLOW==0.0 .AND. Q_CUTTED>=Q_CUT) THEN   ! NO MAX FLOW LIMIT FOR ALL POWERHOUSE, FLOW ADDED BY PERCENTAGE        
          DO ig = 1,NGPH
            IF (QGT(PHGTNO(ig))>0.0) THEN
              QPH_ADDED = QPH_ADDED + Q_CUTTED*QGT(PHGTNO(ig))/Q_PH
              QGT(PHGTNO(ig)) = QGT(PHGTNO(ig)) + Q_CUTTED*QGT(PHGTNO(ig))/Q_PH
            END IF
          END DO    
        ELSE IF (SUM_PH_MAXFLOW/=0.0 .AND. Q_CUTTED>=Q_CUT) THEN
          DO ig =1, NGPH
            IF (PHMAXFLOW(ig)/=0.0 .AND. QGT(PHGTNO(ig))>0.0) THEN
              QPH_ADDED = QPH_ADDED + Q_CUTTED*(PHMAXFLOW(ig)-QGT(PHGTNO(ig)))/QPH_AVL
              QGT(PHGTNO(ig)) = QGT(PHGTNO(ig)) + Q_CUTTED*(PHMAXFLOW(ig)-QGT(PHGTNO(ig)))/QPH_AVL                         
            ELSE IF (QGT(PHGTNO(ig))>0.0) THEN
              QGT(PHGTNO(ig)) = QGT(PHGTNO(ig))+Q_CUTTED
              QPH_ADDED = Q_CUTTED
            END IF
          END DO
        END IF
        
        IF (QPH_ADDED<Q_CUTTED) THEN      ! FLOW BALANCE
          DO ig = 1,NGPH
            IF (QGT(PHGTNO(ig))>0.0 .AND. (QGT(PHGTNO(ig))+Q_CUTTED-QPH_ADDED<PHMAXFLOW(ig))) THEN
              QGT(PHGTNO(ig)) = QGT(PHGTNO(ig))+Q_CUTTED-QPH_ADDED
              QPH_ADDED = Q_CUTTED
              EXIT
            END IF
          END DO
        END IF
        
        CALL SYSTDG_TDG
        IF (ABS(TDG_TDG-tstarget)<=tsconv) THEN
          EXIT
        END IF
        IF (TDG_TDG-tstarget>tsconv) THEN
          Q_MIN = Q_CUT
          QGT = QGTSAVE                 ! NOT CONVERGENT
          CALL SYSTDG_TDG
        ELSE IF (tstarget-TDG_TDG>tsconv) THEN
          Q_MAX = Q_CUT
          QGT = QGTSAVE                 ! NOT CONVERGENT
          CALL SYSTDG_TDG
        END IF
        Q_CUT = 0.5*(Q_MAX+Q_MIN)
        ITERATION = ITERATION + 1
        !
      END DO                            ! END DO WHILE FOR DICHONOMY FLOW CUT
      
      IF (ITERATION==1 .AND. Q_SP==0.0) WRITE (targetfnno+2, '(A,F12.3)') 'SPILL FLOW IS ZERO ON JDAY', JDAY
      IF (ITERATION==1 .AND. TDG_TDG < tstarget+tsconv) WRITE (targetfnno+2,'(A,F12.3)') 'TDG IS LOWER THAN TARGET, NO ITERATION CALCULATION NEEDED ON JDAY', JDAY
      IF (TDG_TDG-tstarget>tsconv .AND. ITERATION>tsiteration .AND. Q_SP>0.0) THEN
        WRITE (targetfnno+2,'(A,F12.3)') 'THE ITERATION LIMIT IS EXCEEDED, AND TDG HAS STILL NOT CONVERGED TO TARGET ON JDAY', JDAY
        WARNING_OPEN = .TRUE.
        QGT = QGTSAVE                    ! UNDO THE FLOW CUT
        CALL SYSTDG_TDG
      END IF
      !
    END IF
    
    IF (JDAY>=NXTSPLIT) THEN
      IF ((TDG_TDG-tstarget)>tsconv) THEN
        CO = 'U'
      ELSE
        IF (CO/='R') CO = ' '
      END IF
      SUM_QGT2 = 0.0
      DO ig = 1,NGT
        SUM_QGT2 = SUM_QGT2+QGT(ig)
      END DO
      !WRITE (targetfnno, '(A, F10.3, 2A, F10.3, A, F9.3, A, <NGT>(F9.3))')' ',JDAY,'  ', CO,TDG_TDG,'  ',SUM_QGT2,'  ',(QGT(ig), ig = 1, NGT)
      WRITE (targetfnno, '(F10.3, ",", A, ",", F10.3, ",", F9.3, ",", <NGT>(F9.3,","))') JDAY, CO, TDG_TDG, SUM_QGT2, (QGT(ig), ig = 1, NGT)
      NXTSPLIT = NXTSPLIT + tsfreq
    END IF
    IF (JDAY>=NXTSPLIT2) NXTSPLIT2 = NXTSPLIT2 + tsfreq
    !
  END IF
  DEALLOCATE(QGTSAVE, SP_ACTIVE)
  IF (NGPH>0) DEALLOCATE(PH_ACTIVE)  
End Subroutine TDGtarget

Subroutine DEALLOCATE_TDGtarget
  USE Selective1TDGtarget; USE MAIN, ONLY: targetfnno
  IMPLICIT NONE
  !
  close(targetfnno)  
  close(targetfnno+2)
  close(targetfnno+3)
  IF (tsdynsel=='      ON') close(targetfnno+3)
  deallocate(SPGTNO, SPPRIOR, SPMINFRAC)
  if (NGPH>0) deallocate(PHGTNO,PHMAXFLOW)
End Subroutine DEALLOCATE_TDGtarget

Subroutine Dy_Priority
  USE Selective1TDGtarget; USE STRUCTURES, ONLY: QGT
  IMPLICIT NONE
  !  
  INTEGER                  :: ib, IBB, ig, TOP_SPNO, TOP_PRIOR, COUNT, CONTU, NGROUP, NLEFT
  INTEGER, DIMENSION(NGSP) :: PRIOR_SP
  LOGICAL, DIMENSION(NGSP) :: PRIOR_UPDATED
  REAL,    ALLOCATABLE, DIMENSION(:)  :: QSP_TEMP
  INTEGER, ALLOCATABLE, DIMENSION(:)  :: GTNO
  TOP_PRIOR = 1
  PRIOR_SP = NGSP+1
  PRIOR_UPDATED = .FALSE.
  
  COUNT = 0
  DO ib = 1,NGSP
    IF (QGT(SPGTNO(ib))>0.0) THEN
      COUNT = COUNT + 1
    ELSE
      PRIOR_UPDATED(ib) = .TRUE.   
    END IF
  END DO
  IF (COUNT>0) ALLOCATE (GTNO(COUNT))
  IF (COUNT>0) ALLOCATE(QSP_TEMP(COUNT))
  GTNO = 0
  COUNT = 0
  DO ib = 1,NGSP
    IF (QGT(SPGTNO(ib))>0.0) THEN
      COUNT = COUNT + 1
      GTNO(COUNT) = ib
      QSP_TEMP(COUNT) = QGT(SPGTNO(ib))
    END IF
  END DO
  IF (COUNT==1) PRIOR_SP(GTNO(1)) = 1
  IF (COUNT>1) THEN
    CALL BUBBLE_SORT(QSP_TEMP, GTNO, COUNT)
    PRIOR_SP(GTNO(1)) = 1
    DO ib = 1, COUNT-1
      IF (QSP_TEMP(ib)/=QSP_TEMP(ib+1)) THEN
        TOP_PRIOR = TOP_PRIOR + 1
        PRIOR_SP(GTNO(ib+1)) = TOP_PRIOR
      ELSE
        PRIOR_SP(GTNO(ib+1)) = TOP_PRIOR
      END IF
    END DO
  END IF
  
  IF (dygroup>1) THEN         ! GROUP THE PRIORITY BY DYGROUP
    NGROUP = INT(COUNT/dygroup)
    NLEFT = MOD(COUNT, dygroup)
    DO ib = 1,NGROUP
      DO ig = 1,dygroup
        PRIOR_SP(GTNO((ib-1)*dygroup+ig)) = ib
      END DO
    END DO
    IF (NLEFT>0) THEN
      DO ig = 1,NLEFT
        PRIOR_SP(GTNO(NGROUP*dygroup+ig)) = NGROUP + 1
      END DO
    END IF
  END IF
  IF (COUNT>0) DEALLOCATE(GTNO,QSP_TEMP)
  SPPRIOR = PRIOR_SP
  RETURN  
End Subroutine Dy_Priority

Subroutine BUBBLE_SORT (QSP_TEMP, GTNO, COUNT)
  IMPLICIT NONE
  INTEGER                   :: COUNT
  REAL, DIMENSION(COUNT)    :: QSP_TEMP
  INTEGER, DIMENSION(COUNT) :: GTNO
  INTEGER                   :: I, J
  REAL                      :: TEMP, TEMP_NO
  
  DO I = COUNT-1, 1, -1
    DO J = 1,I
      IF (QSP_TEMP(J)<QSP_TEMP(J+1)) THEN
        TEMP = QSP_TEMP(J)
        QSP_TEMP(J) = QSP_TEMP(J+1)
        QSP_TEMP(J+1) = TEMP
            
        TEMP_NO   = GTNO(J)
        GTNO(J)   = GTNO(J+1)
        GTNO(J+1) = TEMP_NO
      END IF
    END DO
  END DO
  RETURN
End Subroutine BUBBLE_SORT