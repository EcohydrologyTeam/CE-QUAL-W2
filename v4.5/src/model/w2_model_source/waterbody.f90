
!***********************************************************************************************************************************
!**                                            S U B R O U T I N E    W A T E R B O D Y                                           **
!***********************************************************************************************************************************

SUBROUTINE WATERBODY
  USE GLOBAL; USE GEOMC; USE TVDC; USE LOGICC; USE PREC

! Type declarations
  IMPLICIT NONE
  
  REAL, SAVE, ALLOCATABLE, DIMENSION(:)   :: ELL,    ELR,    CL
  REAL, SAVE, ALLOCATABLE, DIMENSION(:,:) :: QU,     QD
  REAL(R8)                                :: C(KMX,IMX), SS(KMX,IMX)
  REAL                                    :: ELW, EL1
  REAL                                    :: Q1,HT,FRAC,BRTOT,T1L,T2L,U1,B1,B2
  INTEGER                                 :: IUT,IDT,JJW,K,KL,KR

! Allocation declarations

  ALLOCATE (ELL(KMX), ELR(KMX), CL(NCT), QU(KMX,IMX), QD(KMX,IMX))

! Variable initialization

  ELL = 0.0; ELR = 0.0; CL = 0.0; QU = 0.0; QD = 0.0

! Debug variable
!  ncount=0
! End debug

RETURN

!***********************************************************************************************************************************
!**                                               U P S T R E A M   V E L O C I T Y                                               **
!***********************************************************************************************************************************

ENTRY UPSTREAM_VELOCITY
  DO JJB=1,NBR
    IF (UHS(JB) >= US(JJB) .AND. UHS(JB) <= DS(JJB)) EXIT
  END DO
  DO JJW=1,NWB
    IF (JJB >=BS(JJW) .AND. JJB <= BE(JJW)) EXIT
  END DO
  DO K=KTWB(JJW),KB(UHS(JB))+1
    ELL(K) = EL(K,UHS(JB))-SINA(JJB)*DLX(UHS(JB))*0.5
  END DO
  DO K=KT,KB(IU)+1
    ELR(K) = EL(K,IU)+SINA(JB)*DLX(IU)*0.5
  END DO
  ELW = ELWS(UHS(JB))                     !EL(KTWB(JJW),UHS(JB))-Z(UHS(JB))*COSA(JJB)
  EL1  = ELW-SINA(JJB)*DLX(UHS(JB))*0.5
!   IF(SLOPE(JB) /= 0.0)THEN
!   EL1  = ELW-(ELWS(UHS(JB)+1)-ELW)/(0.5*(DLX(UHS(JB)+DLX(JB)+1))*DLX(UHS(JB))*0.5   ! SW 7/17/09
!   ELSE
!   EL1=ELW
!   ENDIF
  ELW = EL1
  KL   = KTWB(JJW)+1
  DO K=KT+1,KB(IU)+1
    IF (ELR(K) >= ELL(KL)) THEN
      IF (KL == KTWB(JJW)+1) THEN
        Q1 = U(KL-1,UHS(JB))*BHR1(KTWB(JJW),UHS(JB))
        IF (KL == KB(UHS(JB))+1 .AND. ELL(KL) < ELR(KB(IU)+1)) THEN
          HT = ELW-ELR(KB(IU)+1)
        ELSE
          HT = H1(KTWB(JJW),UHS(JB))
        END IF
      ELSE
        Q1 = U(KL-1,UHS(JB))*BHR1(KL-1,UHS(JB))
        HT = H1(KL-1,UHS(JB))
      END IF
      IF (K == KT+1) THEN
        EL1         = EL(KTWB(JJW),UHS(JB))-Z(UHS(JB))*COSA(JJB)-SINA(JJB)*DLX(UHS(JB))*0.5
        U(K-1,IU-1) = Q1*((EL1-ELR(K))/HT)/BHR1(KT,IU-1)
      ELSE
        U(K-1,IU-1) = Q1*((EL1-ELR(K))/HT)/BHR1(K-1,IU-1)
      END IF
      EL1 = ELR(K)
      IF (ELR(K) == ELL(KL)) KL = KL+1
    ELSE
      Q1 = 0.0
      DO WHILE (ELR(K) <= ELL(KL))
        IF (KL == KTWB(JJW)+1 .AND. K == KT+1) THEN
          Q1 = Q1+U(KL-1,UHS(JB))*BHR1(KTWB(JJW),UHS(JB))
        ELSE IF (KL == KTWB(JJW)+1) THEN
          Q1 = Q1+U(KL-1,UHS(JB))*BHR1(KTWB(JJW),UHS(JB))*(EL1-ELL(KL))/H1(KTWB(JJW),UHS(JB))
        ELSE
          FRAC = (EL1-ELL(KL))/H1(KL-1,UHS(JB))
          Q1   =  Q1+U(KL-1,UHS(JB))*BHR1(KL-1,UHS(JB))*FRAC
        END IF
        EL1 = ELL(KL)
        KL  = KL+1
        IF (KL > KB(UHS(JB))) EXIT
      END DO
      IF (K == KT+1) THEN
        BRTOT = BHR1(KT,IU-1)
      ELSE
        BRTOT = BHR1(K-1,IU-1)
      END IF
      FRAC = 0.0
!      IF (KL < KMX) THEN                             ! SW 6/29/06
        HT   =  H1(KL-1,UHS(JB))
        FRAC = (EL1-ELR(K))/HT
        IF (KB(UHS(JB)) >= KB(IU-1) .AND. K > KB(IU-1)) FRAC = (EL1-ELL(KL))/HT     ! SW 6/29/06
        IF (KL >= KMX .and. ELR(K) < ELL(KL)) FRAC=1.0          ! SW 6/29/06
        Q1   =  Q1+U(KL-1,UHS(JB))*BHR1(KL-1,UHS(JB))*FRAC
!      ELSE                                           ! SW 6/29/06
!        Q1 = Q1+U(KL-1,UHS(JB))*BHR1(KL-1,UHS(JB))   ! SW 6/29/06
!      END IF                                         ! SW 6/29/06
      U(K-1,IU-1) = Q1/BRTOT
      IF (KL > KB(UHS(JB))) THEN
        IF (FRAC < 1.0 .AND. FRAC /= 0.0) THEN
          IF (K == KB(IU)+1) THEN
            U(K-1,IU-1) = (Q1+U(KL-1,UHS(JB))*BHR1(KL-1,UHS(JB))*(1.0-FRAC))/BRTOT
          ELSE
            U(K,IU-1) = U(KL-1,UHS(JB))*BHR1(KL-1,UHS(JB))*(1.0-FRAC)/BHR(K,IU-1)
          END IF
        END IF
        GO TO 100
      END IF
      EL1 = ELR(K)
    END IF
  END DO
100 CONTINUE

! Debug
!    QL=0.0; QR=0.0
!  do k=ktwb(jjw),kb(uhs(jb))
!    QL=QL+u(k,uhs(jb))*BHR1(k,uhs(jb))
!  enddo
!  do k=kt,kb(iu-1)
!    QR=QR+u(k,iu-1)*BHR1(k,iu-1)
!  enddo
!  if((QR-QL)/QL > 0.02)then
!     if(ncount.eq.0)open(1299,file='debug_out.txt',status='unknown')
!     ncount=ncount+1
!     write(1299,*)'***QL=',QL,' QR=',QR
!     write(1299,*)'UHS(JB)=',uhs(jb),'  iu-1=',iu-1,' ELWS=', elw
!     write(1299,*)'LEFT: K      U     BHR1    ELEV      H1     BHR    ELL'
!     do k=ktwb(jjw),kb(uhs(jb))
!       write(1299,'(i8,f8.3,f8.3,5f8.3)')k,u(k,uhs(jb)),BHR1(k,uhs(jb)),el(k,uhs(jb)),h1(k,uhs(jb)),bhr(k,uhs(jb)),ell(k)
!     enddo
!     write(1299,*)'RIGHT: K     U     BHR1    ELEV      H1     BHR    ELR'
!     do k=kt,kb(iu-1)
!       write(1299,'(i8,f8.3,f8.3,5f8.3)')k,u(k,iu-1),BHR1(k,iu-1),el(k,iu-1),h1(k,iu-1),bhr(k,iu-1),elr(k)
!     enddo
!  end if
! End debug

RETURN

!***********************************************************************************************************************************
!**                                              U P S T R E A M   W A T E R B O D Y                                              **
!***********************************************************************************************************************************

ENTRY UPSTREAM_WATERBODY
  DO JJB=1,NBR
    IF (UHS(JB) >= US(JJB) .AND. UHS(JB) <= DS(JJB)) EXIT
  END DO
  DO JJW=1,NWB
    IF (JJB >=BS(JJW) .AND. JJB <= BE(JJW)) EXIT
  END DO
  DO K=KTWB(JJW),KB(UHS(JB))+1
    ELL(K) = EL(K,UHS(JB))-SINA(JJB)*DLX(UHS(JB))*0.5
  END DO
  DO K=KT,KB(IU)+1
    ELR(K) = EL(K,IU)+SINA(JB)*DLX(IU)*0.5
  END DO
  ELW  = ELWS(UHS(JB))                                !EL(KTWB(JJW),UHS(JB))-Z(UHS(JB))*COSA(JJB)
  EL1  = ELW-SINA(JJB)*DLX(UHS(JB))*0.5
 ! EL1  = ELW-(ELWS(UHS(JB)+1)-ELW)/(0.5*(DLX(UHS(JB)+DLX(JB)+1))*DLX(UHS(JB))*0.5   ! SW 7/17/09
  KL   = KTWB(JJW)+1
  DO K=KT+1,KB(IU)+1
    IF (ELR(K) >= ELL(KL)) THEN
      T1(K-1,IU-1)            = T1(KL-1,UHS(JB))
      T2(K-1,IU-1)            = T2(KL-1,UHS(JB))
      C1S(K-1,IU-1,CN(1:NAC)) = C1S(KL-1,UHS(JB),CN(1:NAC))
      C1(K-1,IU-1,CN(1:NAC))  = C1S(KL-1,UHS(JB),CN(1:NAC))
      C2(K-1,IU-1,CN(1:NAC))  = C1S(KL-1,UHS(JB),CN(1:NAC))
      EL1                     = ELR(K)
      IF (ELR(K) == ELL(KL)) KL = KL+1
    ELSE
      BRTOT = 0.0
      CL    = 0.0
      T1L   = 0.0
      T2L   = 0.0
      DO WHILE (ELR(K) <= ELL(KL))
        IF (KL == KTWB(JJW)+1 .AND. K == KT+1) THEN
          B1 = BH2(KTWB(JJW),UHS(JB))
        ELSE
          B1 = B(KL-1,UHS(JB))*(EL1-ELL(KL))
        END IF
        BRTOT         = BRTOT+B1
        T1L           = T1L+B1*T1(KL-1,UHS(JB))
        T2L           = T2L+B1*T2(KL-1,UHS(JB))
        CL(CN(1:NAC)) = CL(CN(1:NAC))+B1*C1S(KL-1,UHS(JB),CN(1:NAC))
        EL1 = ELL(KL)
        KL  = KL+1
        IF (KL > KB(UHS(JB))+1) EXIT
      END DO
      IF (KL <= KB(UHS(JB))+1) THEN
        B1    = B(KL-1,UHS(JB))*(EL1-ELR(K))
        BRTOT = BRTOT+B1
        IF (BRTOT > 0.0) THEN
          T1(K-1,IU-1)            = (T1L+B1*T1(KL-1,UHS(JB)))/BRTOT
          T2(K-1,IU-1)            = (T2L+B1*T2(KL-1,UHS(JB)))/BRTOT
          C1S(K-1,IU-1,CN(1:NAC)) = (CL(CN(1:NAC))+B1*C1S(KL-1,UHS(JB),CN(1:NAC)))/BRTOT
          C1(K-1,IU-1,CN(1:NAC))  =  C1S(K-1,IU-1,CN(1:NAC))
          C2(K-1,IU-1,CN(1:NAC))  =  C1S(K-1,IU-1,CN(1:NAC))
        ELSE
          T1(K-1,IU-1)            = T1(KL-1,UHS(JB))
          T2(K-1,IU-1)            = T2(KL-1,UHS(JB))
          C1S(K-1,IU-1,CN(1:NAC)) = C1S(KL-1,UHS(JB),CN(1:NAC))
          C1(K-1,IU-1,CN(1:NAC))  = C1S(K-1,IU-1,CN(1:NAC))
          C2(K-1,IU-1,CN(1:NAC))  = C1S(K-1,IU-1,CN(1:NAC))
        END IF
      ELSE
        IF (BRTOT > 0.0) THEN
          T1(K-1,IU-1)            = T1L          /BRTOT
          T2(K-1,IU-1)            = T2L          /BRTOT
          C1S(K-1,IU-1,CN(1:NAC)) = CL(CN(1:NAC))/BRTOT
          C1(K-1,IU-1,CN(1:NAC))  = C1S(K-1,IU-1,CN(1:NAC))
          C2(K-1,IU-1,CN(1:NAC))  = C1S(K-1,IU-1,CN(1:NAC))
        ELSE
          T1(K-1,IU-1)            = T1(KL-1,UHS(JB))
          T2(K-1,IU-1)            = T2(KL-1,UHS(JB))
          C1S(K-1,IU-1,CN(1:NAC)) = C1S(KL-1,UHS(JB),CN(1:NAC))
          C1(K-1,IU-1,CN(1:NAC))  = C1S(K-1,IU-1,CN(1:NAC))
          C2(K-1,IU-1,CN(1:NAC))  = C1S(K-1,IU-1,CN(1:NAC))
        END IF
        EXIT
      END IF
      EL1 = ELR(K)
    END IF
  END DO
RETURN

!***********************************************************************************************************************************
!**                                            D O W N S T R E A M   W A T E R B O D Y                                            **
!***********************************************************************************************************************************

ENTRY DOWNSTREAM_WATERBODY
  DO JJB=1,NBR
    IF (CDHS(JB) >= CUS(JJB) .AND. CDHS(JB) <= DS(JJB)) EXIT
  END DO
  DO JJW=1,NWB
    IF (JJB >= BS(JJW) .AND. JJB <= BE(JJW)) EXIT
  END DO
  KR = KTWB(JJW)+1
  DO K=KT,KB(ID)+1
    ELL(K) = EL(K,ID)-SINA(JB)*DLX(ID)*0.5
  END DO
  DO K=KTWB(JJW),KB(CDHS(JB))+1
    ELR(K) = EL(K,CDHS(JB))+SINA(JJB)*DLX(CDHS(JB))*0.5
  END DO
  ELW = ELWS(ID)                                             !EL(KTWB(JW),ID)-Z(ID)*COSA(JB)
  EL1  = ELW-SINA(JB)*DLX(ID)*0.5
!   EL1  = ELW+(ELW-ELWS(ID-1))/(0.5*(DLX(ID)+DLX(ID-1))*DLX(ID)*0.5   ! SW 7/17/09
  IDT  = ID+1
  DO K=KT+1,KB(ID)+1
    IF (ELL(K) >= ELR(KR)) THEN
      T1(K-1,IDT)            = T1(KR-1,CDHS(JB))
      T2(K-1,IDT)            = T2(KR-1,CDHS(JB))
      C1S(K-1,IDT,CN(1:NAC)) = C1S(KR-1,CDHS(JB),CN(1:NAC))
      C1(K-1,IDT,CN(1:NAC))  = C1S(KR-1,CDHS(JB),CN(1:NAC))
      C2(K-1,IDT,CN(1:NAC))  = C1S(KR-1,CDHS(JB),CN(1:NAC))
      EL1                    = ELL(K)
      IF (ELL(K) == ELR(KR)) KR = KR+1
      IF (KR > KB(CDHS(JB))+1) EXIT
    ELSE
      BRTOT = 0.0
      CL    = 0.0
      T1L   = 0.0
      T2L   = 0.0
      DO WHILE (ELL(K) <= ELR(KR))
        IF (KR == KTWB(JJW)+1 .AND. K == KT+1) THEN
          B1    = BH2(KTWB(JJW),CDHS(JB))
          BRTOT = BRTOT+B1
        ELSE
          B1    = B(KR-1,CDHS(JB))*(EL1-ELR(KR))
          BRTOT = BRTOT+B1
        END IF
        T1L           = T1L+B1*T1(KR-1,CDHS(JB))
        T2L           = T2L+B1*T2(KR-1,CDHS(JB))
        CL(CN(1:NAC)) = CL(CN(1:NAC))+B1*C1S(KR-1,CDHS(JB),CN(1:NAC))
        EL1           = ELR(KR)
        KR            = KR+1
        IF (KR > KB(CDHS(JB))+1) EXIT
      END DO
      IF (KR <= KB(CDHS(JB)+1)) THEN
        B1 = B(KR-1,CDHS(JB))*(EL1-ELL(K))
      ELSE
        B1 = 0.0
      END IF
      BRTOT = BRTOT+B1
      IF (BRTOT == 0.0) EXIT
      T1(K-1,IDT)            = (T1L+B1*T1(KR-1,CDHS(JB)))/BRTOT
      T2(K-1,IDT)            = (T2L+B1*T2(KR-1,CDHS(JB)))/BRTOT
      C1S(K-1,IDT,CN(1:NAC)) = (CL(CN(1:NAC))+B1*C1S(KR-1,CDHS(JB),CN(1:NAC)))/BRTOT
      C1(K-1,IDT,CN(1:NAC))  =  C1S(K-1,IDT,CN(1:NAC))
      C2(K-1,IDT,CN(1:NAC))  =  C1S(K-1,IDT,CN(1:NAC))
      EL1                    =  ELL(K)
    END IF
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                 U P S T R E A M   B R A N C H                                                 **
!***********************************************************************************************************************************

ENTRY UPSTREAM_BRANCH
  DO JJW=1,NWB
    IF (JJB >=BS(JJW) .AND. JJB <= BE(JJW)) EXIT
  END DO
  DO K=KT,KB(I)+1
    ELL(K) = EL(K,I)
  END DO
  DO K=KTWB(JWUH(JB)),KB(CUS(JJB))+1
    ELR(K) = EL(K,CUS(JJB))+SINA(JJB)*DLX(CUS(JJB))*0.5
  END DO
  ELW  = ELWS(CUS(JJB))                               !EL(KTWB(JWUH(JB)),CUS(JJB))-Z(CUS(JJB))*COSA(JJB)
  EL1  = ELW+SINA(JJB)*DLX(CUS(JJB))*0.5
!   EL1  = ELW-(ELWS(CUS(JJB)+1)-ELW)/(0.5*(DLX(CUS(JJB)+DLX(CUS(JJB)+1))*DLX(CUS(JJB))*0.5   ! SW 7/17/09
  KR   = KTWB(JWUH(JB))+1
  DO K=KT+1,KB(I)+1
    IF (ELL(K) >= ELR(KR)) THEN
      Q1 = VOLUH2(KR-1,JJB)/DLT
      B1 = BHR(KR-1,CUS(JJB)-1)
      IF (KR == KTWB(JWUH(JB))+1) B1 = BHR2(KT,CUS(JJB)-1)
      U1 = U(KR-1,CUS(JJB)-1)*B1
      HT = H2(KR-1,JWUH(JB))
      IF (KR == KTWB(JWUH(JB))+1) HT = H2(KT,CUS(JJB))
      Q1 = Q1*((EL1-ELL(K))/HT)
      B2 = BHR(K-1,I)
      IF (K == KTWB(JW)+1) B2 = BHR2(KT,I)
      UXBR(K-1,I) = UXBR(K-1,I)+(ABS((U1/B2)*COS(BETABR)*Q1)/DLX(I))
      UYBR(K-1,I) = UYBR(K-1,I)+ABS(Q1*SIN(BETABR))
      EL1         = ELL(K)
      IF (ELL(K) == ELR(KR)) KR = KR+1
    ELSE
      U1    = 0.0
      Q1    = 0.0
      BRTOT = 0.0
      DO WHILE (ELL(K) <= ELR(KR))
        IF (KR /= KTWB(JWUH(JB))+1) THEN
          FRAC = (EL1-ELR(KR))/(H(KR-1,JWUH(JB)))
          B1   =  BHR(KR-1,CUS(JJB)-1)
        ELSE
          FRAC = (EL1-ELR(KR))/H2(KT,CUS(JJB))
          B1   =  BHR2(KT,CUS(JJB)-1)
        END IF
        U1  = U1+U(KR-1,CUS(JJB)-1)*B1*FRAC
        Q1  = Q1+VOLUH2(KR-1,JJB)/DLT*FRAC
        EL1 = ELR(KR)
        KR  = KR+1
        IF (KR > KB(CUS(JJB)+1)) EXIT
      END DO
      IF (K == KTWB(JW)+1) THEN
        B2 = BHR2(KT,I)
      ELSE
        B2 = BHR(K-1,I)
      END IF
      IF (H(KR-1,JWUH(JB)) /= 0.0) THEN
        IF (KR-1 == KTWB(JWUH(JB))) THEN
          HT = H2(KT,CUS(JJB))
        ELSE
          HT = H2(KR-1,JWUH(JB))
        END IF
        FRAC        = (EL1-ELL(KR-1))/HT
        Q1          =  Q1+FRAC*VOLUH2(KR-1,JJB)/DLT
        UXBR(K-1,I) =  UXBR(K-1,I)+(ABS((U1/B2)*COS(BETABR)*Q1)/DLX(I))
        UYBR(K-1,I) =  UYBR(K-1,I)+ABS(Q1*SIN(BETABR))
      END IF
      IF (KR > KB(CUS(JJB)+1)) EXIT
      EL1 = ELL(K)
    END IF
  END DO
RETURN

!***********************************************************************************************************************************
!**                                               D O W N S T R E A M   B R A N C H                                               **
!***********************************************************************************************************************************

ENTRY DOWNSTREAM_BRANCH
  DO JJW=1,NWB
    IF (JJB >= BS(JJW) .AND. JJB <= BE(JJW)) EXIT
  END DO
  DO K=KT,KB(I)+1
    ELR(K) = EL(K,I)
  END DO
  DO K=KTWB(JJW),KB(DS(JJB))+1
    ELL(K) = EL(K,DS(JJB))+SINA(JJB)*DLX(DS(JJB))*0.5
  END DO
  ELW = ELWS(DS(JJB))                                      !EL(KTWB(JJW),DS(JJB))-Z(DS(JJB))*COSA(JJB)
  EL1  = ELW-SINA(JJB)*DLX(DS(JJB))*0.5
!  IF(SLOPE(JJB) /= 0.0)THEN
!  EL1  = ELW+(ELW-ELWS(DS(JJB)-1))/(0.5*(DLX(DS(JJB)+DLX(DS(JJB)+1))*DLX(DS(JJB))*0.5   ! SW 7/17/09
!  ELSE
!  EL1=ELW
!  ENDIF
  KL   = KTWB(JJW)+1
  DO K=KT+1,KB(I)+1
    IF (ELR(K) >= ELL(KL)) THEN
      Q1 = VOLDH2(KL-1,JJB)/DLT
      B1 = BHR(KL-1,DS(JJB)-1)
      IF (KL == KTWB(JJW)+1) B1 = BHR2(KT,DS(JJB))
      U1 = U(KL-1,DS(JJB)-1)*B1
      HT = H2(KL-1,JJW)
      IF (KL == KTWB(JJW)+1) HT = H2(KT,DS(JJB))
      Q1 = Q1*((EL1-ELR(K))/HT)
      B2 = BHR(K-1,I)
      IF (K == KTWB(JW)+1) B2 = BHR2(KT,I)
      UXBR(K-1,I) = UXBR(K-1,I)+(ABS((U1/B2)*COS(BETABR)*Q1)/DLX(I))
      UYBR(K-1,I) = UYBR(K-1,I)+ABS(Q1*SIN(BETABR))
      EL1 = ELR(K)
      IF (ELR(K) == ELL(KL)) KL = KL+1
    ELSE
      U1    = 0.0
      Q1    = 0.0
      BRTOT = 0.0
      DO WHILE (ELR(K) <= ELL(KL))
        IF (KL /= KTWB(JJW)+1) THEN
          FRAC = (EL1-ELL(KL))/H(KL-1,JJW)
          B1   = BHR(KL-1,DS(JJB)-1)
        ELSE
          FRAC = (EL1-ELL(KL))/H2(KT,DS(JJB))
          B1   = BHR2(KT,DS(JJB))
        END IF
        U1  = U1+U(KL-1,DS(JJB))*B1*FRAC
        Q1  = Q1+(VOLDH2(KL-1,JJB)/DLT)*FRAC
        EL1 = ELL(KL)
        KL  = KL+1
        IF (KL > KB(DS(JJB)+1)) EXIT
      END DO
      IF (K == KTWB(JW)+1) THEN
        B2 = BHR2(KT,I)
      ELSE
        B2 = BHR(K-1,I)
      END IF
      IF (H(KL-1,JJW) /= 0.0) THEN
        IF (KL-1 == KTWB(JJW)) THEN
          HT = H2(KT,DS(JJB))
        ELSE
          HT = H2(KL-1,JJW)
        END IF
        FRAC        = (EL1-ELL(KL-1))/HT
        Q1          =  Q1+FRAC*(VOLDH2(KL-1,JJB)/DLT)
        UXBR(K-1,I) =  UXBR(K-1,I)+(ABS((U1/B2)*COS(BETABR)*Q1)/DLX(I))
        UYBR(K-1,I) =  UYBR(K-1,I)+ABS(Q1*SIN(BETABR))
      END IF
      IF (KL > KB(DS(JJB)+1)) EXIT
      EL1 = ELR(K)
    END IF
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                   U P S T R E A M   F L O W                                                   **
!***********************************************************************************************************************************

ENTRY UPSTREAM_FLOW
  DO JJB=1,NBR
    IF (UHS(JB) >= US(JJB) .AND. UHS(JB) <= DS(JJB)) EXIT
  END DO
  DO JJW=1,NWB
    IF (JJB >=BS(JJW) .AND. JJB <= BE(JJW)) EXIT
  END DO
  DO K=KTWB(JWUH(JB)), KB(UHS(JB))+1
    ELL(K) = EL(K,UHS(JB))
  END DO
  DO K=KT,KB(IU)+1
    ELR(K) = EL(K,IU)+SINA(JB)*DLX(IU)*0.5
  END DO
  ELW = ELWS(IU)                                  !EL(KTWB(JW),IU)-Z(IU)*COSA(JB)
  EL1  = ELW+SINA(JB)*DLX(IU)*0.5
!  IF(SLOPE(JB) /= 0.0)THEN
!  EL1  = ELW-(ELWS(IU+1)-ELW)/(0.5*(DLX(IU)+DLX(IU+1))*DLX(IU)*0.5   ! SW 7/17/09
!  ELSE
!  EL1=ELW
!  ENDIF
  KR   = KT+1
  DO K=KTWB(JWUH(JB))+1,KB(UHS(JB))+1
    IF (ELL(K) >= ELR(KR)) THEN
      Q1 = VOLUH2(KR-1,JB)/DLT
      HT = H2(KR-1,JW)
      IF (KR == KTWB(JW)+1) THEN
        QU(K-1,UHS(JB))  = Q1
        QSS(K-1,UHS(JB)) = QSS(K-1,UHS(JB))-Q1
      ELSE
        QU(K-1,UHS(JB))  = Q1*((EL1-ELL(K))/HT)
        QSS(K-1,UHS(JB)) = QSS(K-1,UHS(JB))-QU(K-1,UHS(JB))
      END IF
      EL1 = ELL(K)
      IF (ELL(K) == ELR(KR)) KR = KR+1
    ELSE
      Q1 = 0.0
      DO WHILE (ELL(K) <= ELR(KR))
        IF (KR /= KTWB(JW)+1) THEN
          FRAC = (EL1-ELR(KR))/(H(KR-1,JW))
        ELSE
          FRAC = (EL1-ELR(KR))/(H2(KT,IU))
        END IF
        Q1  = Q1+VOLUH2(KR-1,JB)/DLT*FRAC
        EL1 = ELR(KR)
        KR  = KR+1
        IF (KR > KB(IU)+1) EXIT
      END DO
      IF (H(KR-1,JW) /= 0.0) THEN
        FRAC = (EL1-ELL(K))/H(KR-1,JW)
        QU(K-1,UHS(JB)) = Q1+FRAC*VOLUH2(KR-1,JB)/DLT
      ELSE
        QU(K-1,UHS(JB)) = Q1
      END IF
      QSS(K-1,UHS(JB)) = QSS(K-1,UHS(JB))-QU(K-1,UHS(JB))
      IF (KR > KB(IU)+1) EXIT
      EL1 = ELL(K)
    END IF
  END DO
RETURN

!***********************************************************************************************************************************
!**                                                 D O W N S T R E A M   F L O W                                                 **
!***********************************************************************************************************************************

ENTRY DOWNSTREAM_FLOW
  DO JJB=1,NBR
    IF (CDHS(JB) >= US(JJB) .AND. CDHS(JB) <= DS(JJB)) EXIT
  END DO
  DO JJW=1,NWB
    IF (JJB >=BS(JJW) .AND. JJB <= BE(JJW)) EXIT
  END DO
  DO K=KTWB(JJW),KB(CDHS(JB))+1
    ELR(K) = EL(K,CDHS(JB))
  END DO
  DO K=KT,KB(ID)+1
    ELL(K) = EL(K,ID)-SINA(JB)*DLX(ID)*0.5
  END DO
  ELW  = ELWS(ID)                            !EL(KTWB(JW),ID)-Z(ID)*COSA(JB)
  EL1  = ELW-SINA(JB)*DLX(ID)*0.5
!  IF(SLOPE(JB) /= 0.0)THEN
!  EL1  = ELW+(ELW-ELWS(ID-1))/(0.5*(DLX(ID)+DLX(ID-1))*DLX(ID)*0.5   ! SW 7/17/09
!  ELSE
!  EL1=ELW
!  ENDIF
  KL   = KT+1
  DO K=KTWB(JJW)+1,KB(CDHS(JB))+1
    IF (ELR(K) >= ELL(KL)) THEN
      Q1 = VOLDH2(KL-1,JB)/DLT
      HT = H2(KL-1,JW)
      IF (KL == KTWB(JW)+1) HT = H2(KT,ID)
      QD(K-1,CDHS(JB))  = Q1*((EL1-ELR(K))/HT)
      QSS(K-1,CDHS(JB)) = QSS(K-1,CDHS(JB))+QD(K-1,DHS(JB))
      EL1               = ELR(K)
      IF (ELR(K) == ELL(KL)) KL = KL+1
    ELSE
      Q1 = 0.0
      DO WHILE (ELR(K) <= ELL(KL))
        IF (KL /= KTWB(JW)+1) THEN
          FRAC = (EL1-ELL(KL))/(H(KL-1,JW))
        ELSE
          FRAC = (EL1-ELL(KL))/(H2(KT,ID))
        END IF
        Q1  = Q1+(VOLDH2(KL-1,JB)/DLT)*FRAC
        EL1 = ELL(KL)
        KL  = KL+1
        IF (KL > KB(ID)+1) EXIT
      END DO
      IF (H(KL-1,JW) /= 0.0) THEN
        FRAC             = (EL1-ELR(K))/H(KL-1,JW)
        QD(K-1,CDHS(JB)) =  Q1+FRAC*(VOLDH2(KL-1,JB)/DLT)
      ELSE
        QD(K-1,CDHS(JB)) = Q1
      END IF
      QSS(K-1,CDHS(JB)) = QSS(K-1,CDHS(JB))+QD(K-1,CDHS(JB))
      IF (KL > KB(ID)+1) EXIT
      EL1 = ELR(K)
    END IF
  END DO
RETURN

!***********************************************************************************************************************************
!**                                            U P S T R E A M   C O N S T I T U E N T                                            **
!***********************************************************************************************************************************

ENTRY UPSTREAM_CONSTITUENT(C,SS)
  DO JJB=1,NBR
    IF (UHS(JB) >= US(JJB) .AND. UHS(JB) <= DS(JJB)) EXIT
  END DO
  DO JJW=1,NWB
    IF (JJB >=BS(JJW) .AND. JJB <= BE(JJW)) EXIT
  END DO
  DO K=KTWB(JWUH(JB)),KB(UHS(JB))+1
    ELL(K) = EL(K,UHS(JB))
  END DO
  DO K=KT,KB(IU)+1
    ELR(K) = EL(K,IU)+SINA(JB)*DLX(IU)*0.5
  END DO
  ELW  = ELWS(IU)                            !EL(KTWB(JW),IU)-Z(IU)*COSA(JB)
  EL1  = ELW+SINA(JB)*DLX(IU)*0.5
!  IF(SLOPE(JB) /= 0.0)THEN
!  EL1  = ELW-(ELWS(IU+1)-ELW)/(0.5*(DLX(IU)+DLX(IU+1))*DLX(IU)*0.5   ! SW 7/17/09
!  ELSE
!  EL1=ELW
!  ENDIF
  KR   = KT+1
  DO K=KTWB(JWUH(JB))+1,KB(UHS(JB))+1
    IUT = IU
    IF (QU(K-1,UHS(JB)) >= 0.0) IUT = IU-1
    IF (ELL(K) >= ELR(KR)) THEN
      T1L             = C(KR-1,IUT)
      SS(K-1,UHS(JB)) = SS(K-1,UHS(JB))-T1L*QU(K-1,UHS(JB))
      EL1             = ELL(K)
      IF (ELL(K) == ELR(KR)) KR = KR+1
    ELSE
      T1L   = 0.0
      BRTOT = 0.0
      DO WHILE (ELL(K) <= ELR(KR))
        IF (KR == KT+1 .AND. K == KTWB(JWUH(JB))+1) THEN
          B1    = BH2(KT,IU)
          BRTOT = BRTOT+B1
        ELSE
          B1    = B(KR-1,IU)*(EL1-ELR(KR))
          BRTOT = BRTOT+B1
        END IF
        IUT = IU
        IF (QU(K-1,UHS(JB)) >= 0.0) IUT = IU-1
        T1L = T1L+B1*C(KR-1,IUT)
        EL1 = ELR(KR)
        KR  = KR+1
        IF (KR > KB(IU)+1) EXIT
      END DO
      IUT = IU
      IF (QU(K-1,UHS(JB)) >= 0.0) IUT = IU-1
      B1              =  B(KR-1,IU)*(EL1-ELL(K))
      BRTOT           =  BRTOT+B1
      T1L             = (T1L+B1*C(KR-1,IUT))/BRTOT
      SS(K-1,UHS(JB)) = TSS(K-1,UHS(JB))-T1L*QU(K-1,UHS(JB))
      IF (KR > KB(IU)+1) EXIT
      EL1 = ELL(K)
    END IF
  END DO
RETURN

!***********************************************************************************************************************************
!**                                          D O W N S T R E A M   C O N S T I T U E N T                                          **
!***********************************************************************************************************************************

ENTRY DOWNSTREAM_CONSTITUENT (C,SS)
  DO JJB=1,NBR
    IF (CDHS(JB) >= US(JJB) .AND. CDHS(JB) <= DS(JJB)) EXIT
  END DO
  DO JJW=1,NWB
    IF (JJB >=BS(JJW) .AND. JJB <= BE(JJW)) EXIT
  END DO
  DO K=KTWB(JJW),KB(CDHS(JB))+1
    ELR(K) = EL(K,CDHS(JB))
  END DO
  DO K=KT,KB(ID)+1
    ELL(K) = EL(K,ID)+SINA(JB)*DLX(ID)*0.5
  END DO
  ELW  = ELWS(ID)                                     !EL(KTWB(JW),ID)-Z(ID)*COSA(JB)
  EL1  = ELW+SINA(JB)*DLX(ID)*0.5
!  IF(SLOPE(JB) /= 0.0)THEN
!  EL1  = ELW+(ELW-ELWS(ID-1))/(0.5*(DLX(ID)+DLX(ID-1))*DLX(ID)*0.5   ! SW 7/17/09
!  ELSE
!  EL1=ELW
!  ENDIF
  KL   = KT+1
  DO K=KTWB(JJW)+1,KB(CDHS(JB))+1
    IDT = ID+1
    IF (QD(K-1,CDHS(JB)) >= 0.0) IDT = ID
    IF (ELR(K) >= ELL(KL)) THEN
      T1L              = C(KL-1,IDT)
      SS(K-1,CDHS(JB)) = SS(K-1,CDHS(JB))+T1L*QD(K-1,CDHS(JB))
      EL1              = ELR(K)
      IF (ELR(K) == ELL(KL)) KL = KL+1
    ELSE
      T1L   = 0.0
      BRTOT = 0.0
      DO WHILE (ELR(K) <= ELL(KL))
        IF (KL == KTWB(JW)+1 .AND. K == KTWB(JJW)+1) THEN
          B1    = BH2(KT,ID)
          BRTOT = BRTOT+B1
        ELSE
          B1    = B(KL-1,ID)*(EL1-ELL(KL))
          BRTOT = BRTOT+B1
        END IF
        T1L = T1L+B1*C(KL-1,IDT)
        EL1 = ELL(KL)
        KL  = KL+1
        IF (KL > KB(ID)) THEN
          B1               =  B(KL-1,ID)*(EL1-ELR(K))
          BRTOT            =  BRTOT+B1
          T1L              = (T1L+B1*C(KL-1,IDT))/BRTOT
          SS(K-1,CDHS(JB)) =  SS(K-1,CDHS(JB))+T1L*QD(K-1,CDHS(JB))
          GO TO 200
        END IF
      END DO
      B1               =  B(KL-1,ID)*(EL1-ELR(K))
      BRTOT            =  BRTOT+B1
      T1L              = (T1L+B1*C(KL-1,IDT))/BRTOT
      SS(K-1,CDHS(JB)) =  SS(K-1,CDHS(JB))+T1L*QD(K-1,CDHS(JB))
      EL1              =  ELL(KL)
    END IF
  END DO
200 CONTINUE
RETURN
ENTRY DEALLOCATE_WATERBODY
  DEALLOCATE (ELL, ELR, CL, QU, QD)
RETURN
END SUBROUTINE WATERBODY
