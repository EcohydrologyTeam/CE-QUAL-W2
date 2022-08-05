
!***********************************************************************************************************************************
!**                                            S U B R O U T I N E   W I T H D R A W A L                                          **
!***********************************************************************************************************************************

SUBROUTINE WITHDRAWAL
  USE GLOBAL; USE GEOMC; USE TVDC; USE SELWC; USE LOGICC
  USE MAIN, ONLY: DERIVED_CALC,CDN,TDGON,JSG,NNSG,NDO,JWD,EA,SYSTDG,NN2,NDGP,O2DG_DER,TDG_DER
  USE modSYSTDG, ONLY: GTNAME, UPDATE_TDGC, SYSTDG_TDG;    ! systdg
  USE SCREENC, ONLY: JDAY
  USE KINETIC, ONLY: CAC
  IMPLICIT NONE
  REAL :: HSWT,HSWB,ELR,WSEL,ELSTR,COEF,RATIO,HT,RHOFT,DLRHOT,HB,RHOFB,DLRHOB,VSUM,DLRHOMAX,HWDT,HWDB,ELWD,TEMPEST,ESTRTEST,QSUMJS
  REAL :: FRACV,QSUMWD
  REAL(R8)  :: dosat, n2sat   ! cb 11/7/17
  INTEGER :: K,JS,KSTR,KTOP,KBOT,KWD,JJWD     
RETURN

!***********************************************************************************************************************************
!**                                             D O W N S T R E A M   W I T H D R A W A L                                         **
!***********************************************************************************************************************************

ENTRY DOWNSTREAM_WITHDRAWAL (JS)

! Variable initialization

  HSWT = 0.0; HSWB = 0.0; VNORM = 0.0; QNEW = 0.0

! Water surface elevation

  ELR  = SINA(JB)*DLX(ID)*0.5
  WSEL = ELWS(ID)-ELR                   !EL(KT,ID)-Z(ID)*COSA(JB)

! Structure layer

  DO K=KT,KB(ID)
    IF (EL(K,ID)-ELR < ESTR(JS,JB)) EXIT
  END DO
  KSTR = MAX(K-1,KT)
  KSTR = MIN(KSTR,KB(ID))

! Initial withdrawal limits

  KTOP = MAX(KTSW(JS,JB),KT)
  IF (KSTR < KTOP) KTOP = KSTR
  KBOT = MIN(KBSW(JS,JB),KB(ID))
  IF (KBOT <= KT .AND. KBOT /= KB(ID)) KBOT = KT+1
  IF (KBOT > KB(ID)) KBOT = KB(ID)
  ELSTR = ESTR(JS,JB)
  IF (ESTR(JS,JB) <= EL(KB(ID)+1,ID+1)-ELR) THEN
    KSTR  = KB(ID)
    ELSTR = EL(KB(ID),ID)-ELR
  END IF
  IF (ESTR(JS,JB) > EL(KT,ID)-ELR) ELSTR = WSEL
  IF (KBSW(JS,JB) < KSTR) THEN
    KSTR  = KT
    ELSTR = WSEL
  END IF

! Boundary interference

  COEF = 1.0
  IF (WSEL-(EL(KBOT,ID)-ELR) /= 0.0) THEN   ! SR 11/2021
    RATIO = (ELSTR-(EL(KBOT,ID)-ELR))/(WSEL-(EL(KBOT,ID)-ELR))
    IF (RATIO < 0.1 .OR. RATIO > 0.9) COEF = 2.0
  END IF

! Withdrawal zone above structure

  DO K=KSTR-1,KTOP,-1

!** Density frequency

    HT    = (EL(K,ID)-ELR)-ELSTR
    RHOFT = MAX(SQRT((ABS(RHO(K,ID)-RHO(KSTR,ID)))/(HT*RHO(KSTR,ID)+NONZERO)*G),NONZERO)

!** Thickness

    IF (POINT_SINK(JS,JB)) THEN
      HSWT = (COEF*QSTR(JS,JB)/RHOFT)**0.333333
    ELSE
      HSWT = SQRT(2.0*COEF*QSTR(JS,JB)/(WSTR(JS,JB)*RHOFT))
    END IF
    IF (HT >= HSWT) THEN
      KTOP = K; EXIT
    END IF
  END DO

! Reference density

  IF ((ELSTR+HSWT) < WSEL) THEN
    DLRHOT = ABS(RHO(KSTR,ID)-RHO(KTOP,ID))
  ELSE IF (WSEL == ELSTR) THEN
    DLRHOT = NONZERO
  ELSE
    DLRHOT = ABS(RHO(KSTR,ID)-RHO(KT,ID))*HSWT/(WSEL-ELSTR)
  END IF
  DLRHOT = MAX(DLRHOT,NONZERO)

! Withdrawal zone below structure

  DO K=KSTR+1,KBOT

!** Density frequency

    HB    = ELSTR-(EL(K,ID)-ELR)
    RHOFB = MAX(SQRT((ABS(RHO(K,ID)-RHO(KSTR,ID)))/(HB*RHO(KSTR,ID)+NONZERO)*G),NONZERO)

!** Thickness

    IF (POINT_SINK(JS,JB)) THEN
      HSWB = (COEF*QSTR(JS,JB)/RHOFB)**0.333333
    ELSE
      HSWB = SQRT(2.0*COEF*QSTR(JS,JB)/(WSTR(JS,JB)*RHOFB))
    END IF
    IF (HB >= HSWB) THEN
      KBOT = K; EXIT
    END IF
  END DO

! Reference density

  IF ((ELSTR-HSWB) > EL(KBOT+1,ID)) THEN
    DLRHOB = ABS(RHO(KSTR,ID)-RHO(KBOT,ID))
  ELSE IF ((EL(KBOT+1,ID)-ELR) == ELSTR) THEN                                                                          !SR 03/24/13
    DLRHOB = NONZERO                                                                                                   !SR 03/24/13
  ELSE
    DLRHOB = ABS(RHO(KSTR,ID)-RHO(KBOT,ID))*HSWB/(ELSTR-(EL(KBOT+1,ID)-ELR))
  END IF
  DLRHOB = MAX(DLRHOB,NONZERO)

! Velocity profile

  VSUM     = 0.0
!  DLRHOMAX = MAX(DLRHOT,DLRHOB,1.0E-10)                      ! GH 1/31/08
  DO K=KTOP,KBOT
!    VNORM(K) = ABS(1.0-((RHO(K,ID)-RHO(KSTR,ID))/DLRHOMAX)**2)*BHR2(K,ID)
 	   IF(K.GT.KSTR)THEN
       DLRHOMAX = MAX(DLRHOB,1.0E-10)                          !GH 1/31/08
       ELSE
       DLRHOMAX = MAX(DLRHOT,1.0E-10)                          !GH 1/31/08
       ENDIF
     VNORM(K) = 1.0-((RHO(K,ID)-RHO(KSTR,ID))/DLRHOMAX)**2
 	 IF(VNORM(K).GT.1.0) VNORM(K)=1.0                         !GH 1/31/08
	 IF(VNORM(K).LT.0.0) VNORM(K)=0.0                         !GH 1/31/08
	 VNORM(K)=VNORM(K)*BHR2(K,ID)
     VSUM     = VSUM+VNORM(K)
  END DO

! OUTFLOWS
  QSUMJS=0.0                                                  ! SW 7/30/09
  TAVG(JS,JB)=0.0                                                    ! CB 5/12/10
  IF(CONSTITUENTS)CAVG(JS,JB,CN(1:NAC))=0.0
  IF(DERIVED_CALC)CDAVG(JS,JB,CDN(1:NACD(JW),JW))=0.0
  DO K=KTOP,KBOT
    QNEW(K)    = (VNORM(K)/VSUM)*QSTR(JS,JB)
    QOUT(K,JB) =  QOUT(K,JB)+QNEW(K)
    TAVG(JS,JB)=TAVG(JS,JB)+QNEW(K)*T2(K,ID)                  ! SW 7/30/09
    IF(CONSTITUENTS)CAVG(JS,JB,CN(1:NAC))=CAVG(JS,JB,CN(1:NAC))+QNEW(K)*C2(K,ID,CN(1:NAC))  
    IF(DERIVED_CALC)CDAVG(JS,JB,CDN(1:NACD(JW),JW))=CDAVG(JS,JB,CDN(1:NACD(JW),JW))+QNEW(K)*CD(K,ID,CDN(1:NACD(JW),JW))
    QSUMJS=QSUMJS+QNEW(K)
  END DO
IF(QSUMJS.GT.0.0)THEN
  TAVG(JS,JB)=TAVG(JS,JB)/QSUMJS
  IF(CONSTITUENTS)then                    ! cb 1/16/13
    CAVG(JS,JB,CN(1:NAC))=CAVG(JS,JB,CN(1:NAC))/QSUMJS
    if(tdgon)then
     IF (nnsg==1) THEN                                ! nnsg==1 is gate flow
           !
           ! systedg
           IF (SYSTDG) THEN                          
               IF(GTNAME(jsg)) THEN                   
                  CALL  UPDATE_TDGC (0,palt(id),jsg,tavg(js,jb),cavg(js,jb,NDO))    
                  CALL  UPDATE_TDGC (1,palt(id),jsg,tavg(js,jb),cavg(js,jb,NN2))  
                  CALL  UPDATE_TDGC (2,palt(id),jsg,tavg(js,jb),cavg(js,jb,NDGP))
               ELSE
                  call total_dissolved_gas (0,palt(id),nnsg,jsg,tavg(js,jb),cavg(js,jb,NDO))    
                  call total_dissolved_gas (1,palt(id),nnsg,jsg,tavg(js,jb),cavg(js,jb,NN2))     ! n2 GAS
                  call total_dissolved_gas (2,palt(id),nnsg,jsg,tavg(js,jb),cavg(js,jb,NDGP)) 
               END IF
           ELSE
               call total_dissolved_gas (0,palt(id),nnsg,jsg,tavg(js,jb),cavg(js,jb,NDO))    
               call total_dissolved_gas (1,palt(id),nnsg,jsg,tavg(js,jb),cavg(js,jb,NN2))     ! n2 GAS
               call total_dissolved_gas (2,palt(id),nnsg,jsg,tavg(js,jb),cavg(js,jb,NDGP))
           END IF
           !
      ELSE
      call total_dissolved_gas (0,palt(id),nnsg,jsg,tavg(js,jb),cavg(js,jb,NDO))    
      call total_dissolved_gas (1,palt(id),nnsg,jsg,tavg(js,jb),cavg(js,jb,NN2))     ! n2 GAS
      call total_dissolved_gas (2,palt(id),nnsg,jsg,tavg(js,jb),cavg(js,jb,NDGP))
      END IF
    end if
  end if
  IF(DERIVED_CALC)then                    ! cb 1/16/13
    CDAVG(JS,JB,CDN(1:NACD(JW),JW))=CDAVG(JS,JB,CDN(1:NACD(JW),JW))/QSUMJS
    !if(tdgon)then                  ! cb 11/6/17
      !cdavg(js,jb,16)  = (cavg(js,jb,ndo)/exp(7.7117-1.31403*(log(tavg(js,jb)+45.93)))*palt(id))*100.0 
      dosat=exp(7.7117-1.31403*(log(tavg(js,jb)+45.93)))*palt(id)
      cdavg(js,jb,O2DG_DER)=(cavg(js,jb,ndo)/dosat)*100.0 
      !IF(ngctdg /= 0)THEN
      If(CAC(NN2)== '      ON') THEN
          EA = DEXP(2.3026D0*(7.5D0*TDEW(JW)/(TDEW(JW)+237.3D0)+0.6609D0))*0.001316   ! in mm Hg   0.0098692atm=7.5006151mmHg  
          !cdavg(js,jb,NDC)  = (cavg(js,jb,NGN2)/(1.5568D06*0.79*(PALT(ID)-EA)*(1.8816D-5 - 4.116D-7 * Tavg(js,jb) + 4.6D-9 * Tavg(js,jb)**2)))*100.0    ! SW 10/27/15      
          n2sat=1.5568D06*0.79*(PALT(ID)-EA)*(1.8816D-5 - 4.116D-7 * Tavg(js,jb) + 4.6D-9 * Tavg(js,jb)**2)
          cdavg(js,jb,TDG_DER)  = 100.*(0.79*(cavg(js,jb,NN2)/n2sat) + 0.21*(cavg(js,jb,ndo)/dosat))
      ELSE IF(CAC(NDGP)== '      ON') THEN
          cdavg(js,jb,TDG_DER)  = cavg(js,jb,NDGP)/palt(id)*100.0
      END IF
      !ENDIF
    !end if
  end if
ELSE
  TAVG(JS,JB)=-99.0
  IF(CONSTITUENTS)CAVG(JS,JB,CN(1:NAC))=-99.0
  IF(DERIVED_CALC)CDAVG(JS,JB,CDN(1:NACD(JW),JW))=-99.0
END IF

! Inactive layers and total outflow

  IF (JS == NST) THEN
    WHERE (QOUT(:,JB) == 0.0) U(:,ID) = 0.0
  END IF
RETURN
!***********************************************************************************************************************************
!**                                             D O W N S T R E A M   W I T H D R A W A L  ESTIMATE                               **
!***********************************************************************************************************************************

ENTRY DOWNSTREAM_WITHDRAWAL_ESTIMATE(JS,TEMPEST,ESTRTEST)

! VARIABLE INITIALIZATION

  HSWT = 0.0; HSWB = 0.0; VNORM = 0.0; QNEW = 0.0

! Water surface elevation

  ELR  = SINA(JB)*DLX(ID)*0.5
  WSEL = EL(KT,ID)-Z(ID)*COSA(JB)-ELR

! Structure layer

  DO K=KT,KB(ID)
    IF (EL(K,ID)-ELR < estrtest) EXIT
  END DO
  KSTR = MAX(K-1,KT)
  KSTR = MIN(KSTR,KB(ID))

! Initial withdrawal limits

  KTOP = MAX(KTSW(JS,JB),KT)
  IF (KSTR < KTOP) KTOP = KSTR
  KBOT = MIN(KBSW(JS,JB),KB(ID))
  IF (KBOT <= KT .AND. KBOT /= KB(ID)) KBOT = KT+1
  IF (KBOT > KB(ID)) KBOT = KB(ID)                                                                                     !SW 06/03/02
  ELSTR = ESTRTEST
  IF (ESTRTEST <= EL(KB(ID)+1,ID+1)-ELR) THEN                                                                       !SW 10/17/01
    KSTR  = KB(ID)
    ELSTR = EL(KB(ID),ID)-ELR                                                                                          !SW 10/17/01
  END IF
  IF (ESTRTEST > EL(KT,ID)-ELR) ELSTR = WSEL
  IF (KBSW(JS,JB) < KSTR) THEN
    KSTR  = KT
    ELSTR = WSEL                                                                                                       !SW 10/05/00
  END IF

! Boundary interference

  COEF = 1.0
  IF (WSEL-(EL(KBOT,ID)-ELR) /= 0.0) THEN    ! SR 11/2021
    RATIO = (ELSTR-(EL(KBOT,ID)-ELR))/(WSEL-(EL(KBOT,ID)-ELR))                                                         !SW 10/17/01
    IF (RATIO < 0.1 .OR. RATIO > 0.9) COEF = 2.0
  END IF

! Withdrawal zone above structure

  DO K=KSTR-1,KTOP,-1

!** Density frequency

    HT    = (EL(K,ID)-ELR)-ELSTR
    RHOFT = MAX(SQRT((ABS(RHO(K,ID)-RHO(KSTR,ID)))/(HT*RHO(KSTR,ID)+NONZERO)*G),NONZERO)

!** Thickness

    IF (POINT_SINK(JS,JB)) THEN
      HSWT = (COEF*QSTR(JS,JB)/RHOFT)**0.333333
    ELSE
      HSWT = SQRT(2.0*COEF*QSTR(JS,JB)/(WSTR(JS,JB)*RHOFT))
    END IF
    IF (HT >= HSWT) THEN
      KTOP = K; EXIT
    END IF
  END DO

! Reference density

  IF ((ELSTR+HSWT) < WSEL) THEN
    DLRHOT = ABS(RHO(KSTR,ID)-RHO(KTOP,ID))
  ELSE IF (WSEL == ELSTR) THEN
    DLRHOT = NONZERO
  ELSE
    DLRHOT = ABS(RHO(KSTR,ID)-RHO(KT,ID))*HSWT/(WSEL-ELSTR)
  END IF
  DLRHOT = MAX(DLRHOT,NONZERO)

! Withdrawal zone below structure

  DO K=KSTR+1,KBOT

!** Density frequency

    HB    = ELSTR-(EL(K,ID)-ELR)                                                                                       !SW 10/17/01
    RHOFB = MAX(SQRT((ABS(RHO(K,ID)-RHO(KSTR,ID)))/(HB*RHO(KSTR,ID)+NONZERO)*G),NONZERO)

!** Thickness

    IF (POINT_SINK(JS,JB)) THEN
      HSWB = (COEF*QSTR(JS,JB)/RHOFB)**0.333333
    ELSE
      HSWB = SQRT(2.0*COEF*QSTR(JS,JB)/(WSTR(JS,JB)*RHOFB))
    END IF
    IF (HB >= HSWB) THEN
      KBOT = K; EXIT
    END IF
  END DO

! Reference density

  IF ((ELSTR-HSWB) > EL(KBOT+1,ID)) THEN
    DLRHOB = ABS(RHO(KSTR,ID)-RHO(KBOT,ID))
  ELSE IF ((EL(KBOT+1,ID)-ELR) == ELSTR) THEN                                                                          !SR 03/24/13
    DLRHOB = NONZERO                                                                                                   !SR 03/24/13
  ELSE
    DLRHOB = ABS(RHO(KSTR,ID)-RHO(KBOT,ID))*HSWB/(ELSTR-(EL(KBOT+1,ID)-ELR))                                           !SW 10/17/01
  END IF
  DLRHOB = MAX(DLRHOB,NONZERO)

! Velocity profile

  VSUM     = 0.0
!  DLRHOMAX = MAX(DLRHOT,DLRHOB,1.0E-10)
  DO K=KTOP,KBOT
  IF(K.GT.KSTR)THEN
       DLRHOMAX = MAX(DLRHOB,1.0E-10)                          !GH 1/31/08
       ELSE
       DLRHOMAX = MAX(DLRHOT,1.0E-10)                          !GH 1/31/08
       ENDIF
       VNORM(K) = 1.0-((RHO(K,ID)-RHO(KSTR,ID))/DLRHOMAX)**2
     IF(VNORM(K).GT.1.0) VNORM(K)=1.0                          !GH 1/31/08
	 IF(VNORM(K).LT.0.0) VNORM(K)=0.0                          !GH 1/31/08
	 VNORM(K)=VNORM(K)*BHR2(K,ID)
     VSUM= VSUM+VNORM(K)
  END DO

! Outflows

  tempest=0.0
  DO K=KTOP,KBOT
    tempest=tempest+t2(k,id)*(VNORM(K)/VSUM)*QSTR(JS,JB)
  END DO

  if(qstr(js,jb).gt.0.0)tempest=tempest/qstr(js,jb)

RETURN
!***********************************************************************************************************************************
!**                                                L A T E R A L   W I T H D R A W A L                                            **
!***********************************************************************************************************************************

ENTRY LATERAL_WITHDRAWAL

! Variable initialization

  VNORM = 0.0; QSW(:,JWD) = 0.0; HWDT = 0.0; HWDB = 0.0

! Structure layer

  K = KT
  DO K=KT,KB(I)
    IF (EL(K,I) < EWD(JWD)) EXIT
  END DO
  KWD = MAX(K-1,KT)
  KWD = MIN(KWD,KB(I))

! Initial withdrawal limits

  KTOP = MAX(KTWD(JWD),KT)
  IF (KWD < KTOP) KTOP = KWD
  KBOT = MIN(KBWD(JWD),KB(I))
  IF (KBOT <= KT .AND. KB(I) /= KBOT) KBOT = KT+1
  IF (KBOT > KB(I)) KBOT = KB(I)
  ELWD = EWD(JWD)
  IF (EWD(JWD) <= EL(KB(I)+1,I)) THEN
    KWD  = KB(I)
    ELWD = EL(KB(I),I)
  END IF
  IF (EWD(JWD) > EL(KT,I)) ELWD = EL(KT,I)
  IF (KBWD(JWD) < KWD) THEN
    KWD  = KT
    ELWD = EL(KT,I)
  END IF

! Boundary interference

  COEF = 1.0
  IF (KT /= KBOT) THEN
    RATIO = (ELWD-EL(KBOT,I))/(EL(KT,I)-EL(KBOT,I))
    IF (RATIO < 0.1 .OR. RATIO > 0.9) COEF = 2.0
  END IF

! Withdrawal zone above structure

  DO K=KWD-1,KTOP,-1

!** Density frequency

    HT    = EL(K,I)-ELWD
    RHOFT = MAX(SQRT((ABS(RHO(K,I)-RHO(KWD,I)))/(HT*RHO(KWD,I)+NONZERO)*G),NONZERO)

!** Thickness

    HWDT = (COEF*QWD(JWD)/RHOFT)**0.333333
    IF (HT >= HWDT) THEN
      KTOP = K; EXIT
    END IF
  END DO

! Reference density

  IF ((ELWD+HWDT) < EL(KT,I)) THEN
    DLRHOT = ABS(RHO(KWD,I)-RHO(KTOP,I))
  ELSE IF (EL(KT,I) == ELWD) THEN
    DLRHOT = NONZERO
  ELSE
    DLRHOT = ABS(RHO(KWD,I)-RHO(KT,I))*HWDT/(EL(KT,I)-ELWD)
  END IF
  DLRHOT = MAX(DLRHOT,NONZERO)

! Withdrawal zone below structure

  DO K=KWD+1,KBOT

!** Density frequency

    HB    = ELWD-EL(K,I)
    RHOFB = MAX(SQRT((ABS(RHO(K,I)-RHO(KWD,I)))/(HB*RHO(KWD,I)+NONZERO)*G),NONZERO)

!** Thickness

    HWDB = (COEF*QWD(JWD)/RHOFB)**0.333333
    IF (HB >= HWDB) THEN
      KBOT = K; EXIT
    END IF
  END DO

! Reference density

  IF ((ELWD-HWDB) > EL(KBOT+1,I)) THEN
    DLRHOB = ABS(RHO(KWD,I)-RHO(KBOT,I))
  ELSE IF (EL(KBOT+1,I) == ELWD) THEN                                                                                  !SR 03/24/13
    DLRHOB = NONZERO                                                                                                   !SR 03/24/13
  ELSE
    DLRHOB = ABS(RHO(KWD,I)-RHO(KBOT,I))*HWDB/(ELWD-EL(KBOT+1,I))
  END IF
  DLRHOB = MAX(DLRHOB,NONZERO)

! Velocity profile

  VSUM     = 0.0
!  DLRHOMAX = MAX(DLRHOT,DLRHOB,1.0E-10)                                                                             ! SW 1/24/05
  DO K=KTOP,KBOT
!    VNORM(K) = ABS(1.0-((RHO(K,I)-RHO(KWD,I))/DLRHOMAX)**2)*BHR2(K,I)
 	   IF(K.GT.KWD)THEN
       DLRHOMAX = MAX(DLRHOB,1.0E-10)                          !GH 1/31/08
       ELSE
       DLRHOMAX = MAX(DLRHOT,1.0E-10)                          !GH 1/31/08
       ENDIF
     VNORM(K) = 1.0-((RHO(K,I)-RHO(KWD,I))/DLRHOMAX)**2
 	 IF(VNORM(K).GT.1.0) VNORM(K)=1.0                         !GH 1/31/08
	 IF(VNORM(K).LT.0.0) VNORM(K)=0.0                         !GH 1/31/08
	 VNORM(K)=VNORM(K)*BHR2(K,I)
     VSUM     = VSUM+VNORM(K)
  END DO

! Outflows
  QSUMWD=0.0                                                  ! SW 7/30/09
  TAVGW(JWD)=0.0
  IF(CONSTITUENTS)CAVGW(JWD,CN(1:NAC))=0.0
  IF(DERIVED_CALC)CDAVGW(JWD,CDN(1:NACD(JW),JW))=0.0

  DO K=KTOP,KBOT
    FRACV=(VNORM(K)/VSUM)
    QSW(K,JWD) = QSW(K,JWD)+FRACV*QWD(JWD)
    TAVGW(JWD)=TAVGW(JWD)+FRACV*QWD(JWD)*T2(K,I)                  ! SW 7/30/09
    IF(CONSTITUENTS)CAVGW(JWD,CN(1:NAC))=CAVGW(JWD,CN(1:NAC))+FRACV*QWD(JWD)*C2(K,I,CN(1:NAC))  
    IF(DERIVED_CALC)CDAVGW(JWD,CDN(1:NACD(JW),JW))=CDAVGW(JWD,CDN(1:NACD(JW),JW))+FRACV*QWD(JWD)*CD(K,I,CDN(1:NACD(JW),JW))
    QSUMWD=QSUMWD+FRACV*QWD(JWD)
  END DO
  ! Debug
  !if(qwd(jwd)>0.0 .and. qsumwd <= 0.0)then
  !    write(9575,'(A,f8.3,1x,i5,1x,i5,1x,i5,f8.4,1x,f8.4,1x,f10.2)')'JDAY, ktop, kbot, kwd, qwd, qsumwd, elwd:',JDAY, ktop, kbot, kwd, qwd, qsumwd, elwd
  !endif
  ! Debug
  IF(QSUMWD.GT.0.0)THEN
    TAVGW(JWD)=TAVGW(JWD)/QSUMWD               ! SW 7/30/09
    IF(CONSTITUENTS)then                       ! cb 1/16/13
      CAVGW(JWD,CN(1:NAC))=CAVGW(JWD,CN(1:NAC))/QSUMWD  
      if(tdgon)then
        IF (nnsg==1) THEN                                      ! systdg        nnsg==1 is gate flow
            IF (SYSTDG) THEN                                   ! systdg 
                IF (GTNAME(jsg)) THEN                          ! systdg 
                   CALL  UPDATE_TDGC (0,palt(i),jsg,tavgw(jwd),cavgw(jwd,NDO))         ! systdg 
                   CALL  UPDATE_TDGC(1,palt(i),jsg,tavgw(jwd),cavgw(jwd,NN2))          ! systdg 
                   CALL  UPDATE_TDGC(2,palt(i),jsg,tavgw(jwd),cavgw(jwd,NDGP))
                ELSE
                   call total_dissolved_gas (0,palt(i),nnsg,jsg,tavgw(jwd),cavgw(jwd,NDO))  
                   call total_dissolved_gas (1,palt(i),nnsg,jsg,tavgw(jwd),cavgw(jwd,NN2))     ! n2 GAS
                   call total_dissolved_gas (2,palt(i),nnsg,jsg,tavgw(jwd),cavgw(jwd,NDGP)) 
                END IF
            ELSE
        call total_dissolved_gas (0,palt(i),nnsg,jsg,tavgw(jwd),cavgw(jwd,NDO))  
        call total_dissolved_gas (1,palt(i),nnsg,jsg,tavgw(jwd),cavgw(jwd,NN2))     ! n2 GAS
        call total_dissolved_gas (2,palt(i),nnsg,jsg,tavgw(jwd),cavgw(jwd,NDGP))    
            END IF   
        ELSE
            call total_dissolved_gas (0,palt(i),nnsg,jsg,tavgw(jwd),cavgw(jwd,NDO))  
            call total_dissolved_gas (1,palt(i),nnsg,jsg,tavgw(jwd),cavgw(jwd,NN2))     ! n2 GAS
            call total_dissolved_gas (2,palt(i),nnsg,jsg,tavgw(jwd),cavgw(jwd,NDGP)) 
        END IF
      end if
    end if
    IF(DERIVED_CALC)then
      CDAVGW(JWD,CDN(1:NACD(JW),JW))=CDAVGW(JWD,CDN(1:NACD(JW),JW))/QSUMWD
      !if(tdgon)then                ! cb 11/6/17
        !cdavgw(jwd,O2DG_DER)  = (cavgw(jwd,ndo)/exp(7.7117-1.31403*(log(tavgw(jwd)+45.93)))*palt(i))*100.0       
        dosat=exp(7.7117-1.31403*(log(tavgw(jwd)+45.93)))*palt(i)
        cdavgw(jwd,O2DG_DER)=(cavgw(jwd,ndo)/dosat)*100.0 
        If(CAC(NN2)== '      ON') THEN
          EA = DEXP(2.3026D0*(7.5D0*TDEW(JW)/(TDEW(JW)+237.3D0)+0.6609D0))*0.001316   ! in mm Hg   0.0098692atm=7.5006151mmHg  
          !cdavgw(jwd,NDC)  = (cavgw(jwd,NGN2)/(1.5568D06*0.79*(PALT(I)-EA)*(1.8816D-5 - 4.116D-7 * Tavgw(jwd) + 4.6D-9 * Tavgw(jwd)**2)))*100.0    ! SW 10/27/15      
          n2sat=1.5568D06*0.79*(PALT(I)-EA)*(1.8816D-5 - 4.116D-7 * Tavgw(jwd) + 4.6D-9 * Tavgw(jwd)**2)
          cdavgw(jwd,TDG_DER)  = 100.*(0.79*(cavgw(jwd,NN2)/n2sat) + 0.21*(cavgw(jwd,ndo)/dosat))
        ELSE IF(CAC(NDGP)== '      ON') THEN
          cdavgw(jwd,TDG_DER)  =  cavgw(jwd,NDGP)/ palt(i) * 100.0
        END IF  
    end if
  ELSE
    TAVGW(JWD)=-99.0
    IF(CONSTITUENTS)CAVGW(JWD,CN(1:NAC))=-99.0 
    IF(DERIVED_CALC)CDAVGW(JWD,CDN(1:NACD(JW),JW))=-99.0
  ENDIF
  KTW(JWD) = KTOP
  KBW(JWD) = KBOT
  RETURN
!***********************************************************************************************************************************
!**                                                L A T E R A L   W I T H D R A W A L ESTIMATE                                   **
!***********************************************************************************************************************************

  ENTRY LATERAL_WITHDRAWAL_ESTIMATE (JJWD,TEMPEST,ESTRTEST)

! VARIABLE INITIALIZATION

  VNORM = 0.0; QSW(:,JJWD) = 0.0; HWDT = 0.0; HWDB = 0.0

! Structure layer

  K = KT
  DO K=KT,KB(I)
    IF (EL(K,I) < estrtest) EXIT
  END DO
  KWD = MAX(K-1,KT)
  KWD = MIN(KWD,KB(I))

! Initial withdrawal limits

  KTOP = MAX(KTWD(JJWD),KT)
  IF (KWD < KTOP) KTOP = KWD
  KBOT = MIN(KBWD(JJWD),KB(I))
  IF (KBOT <= KT .AND. KB(I) /= KBOT) KBOT = KT+1
  IF (KBOT > KB(I)) KBOT = KB(I)
  ELWD = ESTRTEST
  IF (ESTRTEST <= EL(KB(I)+1,I)) THEN
    KWD  = KB(I)
    ELWD = EL(KB(I),I)
  END IF
  IF (ESTRTEST > EL(KT,I)) ELWD = EL(KT,I)
  IF (KBWD(JJWD) < KWD) THEN
    KWD  = KT
    ELWD = EL(KT,I)
  END IF

! Boundary interference

  COEF = 1.0
  IF (KT /= KBOT) THEN
    RATIO = (ELWD-EL(KBOT,I))/(EL(KT,I)-EL(KBOT,I))
    IF (RATIO < 0.1 .OR. RATIO > 0.9) COEF = 2.0
  END IF

! Withdrawal zone above structure

  DO K=KWD-1,KTOP,-1

!** Density frequency

    HT    = EL(K,I)-ELWD
    RHOFT = MAX(SQRT((ABS(RHO(K,I)-RHO(KWD,I)))/(HT*RHO(KWD,I)+NONZERO)*G),NONZERO)

!** Thickness

    HWDT = (COEF*QWD(JJWD)/RHOFT)**0.333333
    IF (HT >= HWDT) THEN
      KTOP = K; EXIT
    END IF
  END DO

! Reference density

  IF ((ELWD+HWDT) < EL(KT,I)) THEN
    DLRHOT = ABS(RHO(KWD,I)-RHO(KTOP,I))
  ELSE IF (EL(KT,I) == ELWD) THEN
    DLRHOT = NONZERO
  ELSE
    DLRHOT = ABS(RHO(KWD,I)-RHO(KT,I))*HWDT/(EL(KT,I)-ELWD)
  END IF
  DLRHOT = MAX(DLRHOT,NONZERO)

! Withdrawal zone below structure

  DO K=KWD+1,KBOT

!** Density frequency

    HB    = ELWD-EL(K,I)
    RHOFB = MAX(SQRT((ABS(RHO(K,I)-RHO(KWD,I)))/(HB*RHO(KWD,I)+NONZERO)*G),NONZERO)

!** Thickness

    HWDB = (COEF*QWD(JJWD)/RHOFB)**0.333333
    IF (HB >= HWDB) THEN
      KBOT = K; EXIT
    END IF
  END DO

! Reference density

  IF ((ELWD-HWDB) > EL(KBOT+1,I)) THEN
    DLRHOB = ABS(RHO(KWD,I)-RHO(KBOT,I))
  ELSE IF (EL(KBOT+1,I) == ELWD) THEN                                                                                  !SR 03/24/13
    DLRHOB = NONZERO                                                                                                   !SR 03/24/13
  ELSE
    DLRHOB = ABS(RHO(KWD,I)-RHO(KBOT,I))*HWDB/(ELWD-EL(KBOT+1,I))
  END IF
  DLRHOB = MAX(DLRHOB,NONZERO)

! Velocity profile

  VSUM     = 0.0
!  DLRHOMAX = MAX(DLRHOT,DLRHOB,1.0E-10)                                                                             ! SW 1/24/05
  DO K=KTOP,KBOT
!    VNORM(K) = ABS(1.0-((RHO(K,I)-RHO(KWD,I))/DLRHOMAX)**2)*BHR2(K,I)
 	   IF(K.GT.KWD)THEN
       DLRHOMAX = MAX(DLRHOB,1.0E-10)                          !GH 1/31/08
       ELSE
       DLRHOMAX = MAX(DLRHOT,1.0E-10)                          !GH 1/31/08
       ENDIF
     VNORM(K) = 1.0-((RHO(K,I)-RHO(KWD,I))/DLRHOMAX)**2
 	 IF(VNORM(K).GT.1.0) VNORM(K)=1.0                          !GH 1/31/08
	 IF(VNORM(K).LT.0.0) VNORM(K)=0.0                          !GH 1/31/08
	 VNORM(K)=VNORM(K)*BHR2(K,I)
     VSUM = VSUM+VNORM(K)
  END DO

! Outflows

  DO K=KTOP,KBOT
    tempest=tempest+t2(k,i)*(VNORM(K)/VSUM)*QWD(JJWD)
  END DO
  if(qwd(Jjwd).gt.0.0)tempest=tempest/qwd(Jjwd)
  KTW(JJWD) = KTOP
  KBW(JJWD) = KBOT
  return


END SUBROUTINE WITHDRAWAL



MODULE SELECTIVE1 
 REAL                                          :: NXTSTR, NXTTCD, NXTSPLIT,TCDFREQ,TFRQTMP
  CHARACTER(8)                                 :: TEMPC,TSPLTC
  CHARACTER(8), ALLOCATABLE, DIMENSION(:)      :: TCELEVCON,TCYEARLY,TCNTR,TSPLTCNTR,MONCTR,TSYEARLY,DYNSEL,ELCONTSPL,DYNSELSPLT
  INTEGER                                      :: NUMTEMPC,NUMTSPLT, TEMPN        
  INTEGER, ALLOCATABLE, DIMENSION(:)           :: TCNELEV,TCJB,TCJS,TCISEG,TSPLTJB,NOUTS,KSTRSPLT, JBMON, JSMON, NCOUNTCW,SELD
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: TCELEV, TEMPCRIT,QSTRFRAC
  REAL,          ALLOCATABLE, DIMENSION(:)     :: TCTEMP,TCTEND,TCTSRT,TCKLAY,TSPLTT,VOLM,QWDFRAC,TSTEND,TSTSRT,NXSEL,TEMP2,NXSELSPLT,TEMP3
  INTEGER, ALLOCATABLE, DIMENSION(:,:)         :: JSTSPLT, NCOUNTC, JSTSPLTT
  REAL,          ALLOCATABLE, DIMENSION(:,:)  :: VOLMC 
  LOGICAL, ALLOCATABLE, DIMENSION(:)          :: DYNSF,DYNSPF
  REAL, ALLOCATABLE, DIMENSION(:)             :: MINWL   ! Minimum water level above centerline of outlet if TCELEVCON is ON
END MODULE SELECTIVE1

SUBROUTINE SELECTIVEINIT

USE SELECTIVE1;   USE MAIN
  USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART

  IMPLICIT NONE
  
  INTEGER N, IFILE
  REAL DAYTEST
  CHARACTER(1) :: INFORMAT,CHAR1
  CHARACTER(8) :: CHAR8
  CHARACTER(30):: CHAR30

!**                                                   Task 2: Calculations                                                        **
!***********************************************************************************************************************************
      
      IFILE=1949
      TAVG=0.0
      TAVGW=0.0
      DO JB=1,NBR
        IF(NSTR(JB) > 0)THEN
        IFILE=IFILE+1
        WRITE (SEGNUM,'(I0)') JB
        SEGNUM = ADJUSTL(SEGNUM)
        L      = LEN_TRIM(SEGNUM)

        IF(RESTART_IN)THEN
                OPEN  (IFILE,FILE='str_br'//segnum(1:l)//'.csv',POSITION='APPEND')
                JDAY1=0.0
                REWIND (IFILE)
                READ   (IFILE,'(//)',END=13)
                DO WHILE (JDAY1 < JDAY)
                READ (IFILE,*,END=13) JDAY1            !READ (IFILE,'(F10.0)',END=13) JDAY1
                END DO
                BACKSPACE (IFILE)
                13     JDAY1=0.0    
        ELSE
                OPEN  (IFILE,FILE='str_br'//segnum(1:l)//'.csv',status='unknown')
                WRITE(IFILE,*)'Branch:,',jb,', # of structures:,',nstr(jb),', outlet temperatures'
                WRITE(IFILE,'("      JDAY,",<nstr(jb)>(6x,"T(C),"),<nstr(jb)>(3x,"Q(m3/s),"),<nstr(jb)>(4x,"ELEVCL,"))')
        ENDIF
        ENDIF
      END DO
 
      IF(NWD > 0)THEN
       IFILE=IFILE+1  
       IF(RESTART_IN)THEN
                OPEN  (IFILE,FILE='wd_out.opt',POSITION='APPEND')
                JDAY1=0.0
                REWIND (IFILE)
                READ   (IFILE,'(/)',END=14)
                DO WHILE (JDAY1 < JDAY)
                READ (IFILE,'(F10.0)',END=14) JDAY1
                END DO
                BACKSPACE (IFILE)
                14    JDAY1=0.0   
       ELSE
        OPEN  (IFILE,FILE='wd_out.opt',STATUS='unknown')
        WRITE(IFILE,*)'Withdrawals: # of withdrawals:',nwd,' outlet temperatures'
        WRITE(IFILE,'("      JDAY",<nwd>(6x,"T(C)"),<nwd>(3x,"Q(m3/s)"),<nwd>(4x,"ELEVCL"))')
       ENDIF
      end if
 
      
      OPEN(NUNIT,FILE='w2_selective.npt',STATUS='old')
      
      read(NUNIT,'(a)')CHAR1
      read(NUNIT,*)
      read(NUNIT,*)
      IF(CHAR1=='$')GO TO 7770
      ! check for commas another sign that it is comma delimited
      read(NUNIT,'(A)')CHAR30
     DO J=1,30
         IF(CHAR30(J:J)==',')THEN
             CHAR1='$'
             EXIT
         ENDIF
     ENDDO
     BACKSPACE(NUNIT)

7770  CONTINUE
      IF(CHAR1=='$')THEN
      READ(NUNIT,*)CHAR8,TFRQTMP
      READ(NUNIT,*)
      READ(NUNIT,*)
      READ(NUNIT,*)CHAR8,TEMPC,NUMTEMPC,TCDFREQ; TEMPC=ADJUSTR(TEMPC)
      ELSE
      READ(NUNIT,'(8X,F8.0)')TFRQTMP
      READ(NUNIT,'(//8X,A8,I8,F8.0)')TEMPC,NUMTEMPC,TCDFREQ
      ENDIF
      NXTSTR=TMSTRT
      NXTTCD=TMSTRT
      NXTSPLIT=TMSTRT
      
  ALLOCATE (TCNELEV(NUMTEMPC),TCJB(NUMTEMPC),TCJS(NUMTEMPC), TCELEV(NUMTEMPC,11),TCTEMP(NUMTEMPC),TCTEND(NUMTEMPC),TCTSRT(NUMTEMPC),NCOUNTC(NST,NBR),TCISEG(NUMTEMPC),TCKLAY(NUMTEMPC),TCELEVCON(NUMTEMPC)) 
  ALLOCATE (TCYEARLY(NUMTEMPC), JBMON(NUMTEMPC),JSMON(NUMTEMPC),TCNTR(NUMTEMPC)) 
  ALLOCATE (VOLM(NWB),MONCTR(NUMTEMPC),NCOUNTCW(NWD),QWDFRAC(NWD),QSTRFRAC(NST,NBR),DYNSEL(NUMTEMPC),SELD(NUMTEMPC),NXSEL(NUMTEMPC),TEMP2(NUMTEMPC))       
  ALLOCATE (DYNSF(NUMTEMPC),MINWL(NUMTEMPC))    
  DYNSF=.FALSE.;MINWL=0.0
  
      DO J=1,2
      READ(NUNIT,*)
      END DO
      NCOUNTC=0
      DO J=1,NUMTEMPC
            IF(CHAR1=='$')THEN
              READ(NUNIT,*)CHAR8,TCNTR(J),TCJB(J),TCJS(J),TCYEARLY(J),TCTSRT(J),TCTEND(J),TCTEMP(J),TCNELEV(J),(TCELEV(J,N),N=1,TCNELEV(J))
              TCNTR(J)=ADJUSTR(TCNTR(J))
              TCYEARLY=ADJUSTR(TCYEARLY) 
            ELSE
              READ(NUNIT,'(8X,A8,I8,I8,A8,F8.0,F8.0,F8.0,I8,10(F8.0))')TCNTR(J),TCJB(J),TCJS(J),TCYEARLY(J),TCTSRT(J),TCTEND(J),TCTEMP(J),TCNELEV(J),(TCELEV(J,N),N=1,TCNELEV(J))
            ENDIF
            
        IF(TCNTR(J)=='      ST')THEN      
        TCELEV(J,TCNELEV(J)+1)=ESTR(TCJS(J),TCJB(J))   ! ALWAYS PUT THE ORIGINAL ELEVATION AS THE LAST ELEVATION
        ELSE
        TCELEV(J,TCNELEV(J)+1)=EWD(TCJS(J))   ! ALWAYS PUT THE ORIGINAL ELEVATION AS THE LAST ELEVATION
        ENDIF
      END DO
      DO J=1,2
      READ(NUNIT,*)
      END DO
      DO J=1,NUMTEMPC
                      IF(CHAR1=='$')THEN
                                READ(NUNIT,*)CHAR8,TCISEG(J),TCKLAY(J),DYNSEL(J);DYNSEL(J)=ADJUSTR(DYNSEL(J)) 
                      ELSE 
                                READ(NUNIT,'(8X,I8,F8.0,A8)')TCISEG(J),TCKLAY(J),DYNSEL(J) 
                      ENDIF
                      
      END DO
      DO J=1,2
      READ(NUNIT,*)
      END DO
      DO J=1,NUMTEMPC
                      IF(CHAR1=='$')THEN
                                      READ(NUNIT,*)CHAR8,TCELEVCON(J),MINWL(J);tcelevcon(J)=ADJUSTR(tcelevcon(J)) 
                      ELSE 
                                      READ(NUNIT,'(8X,A8,F8.0)')TCELEVCON(J),MINWL(J) 
                      ENDIF
      END DO
      DO J=1,2
      READ(NUNIT,*)
      END DO
                      IF(CHAR1=='$')THEN
                                    READ(NUNIT,*)CHAR8,TSPLTC,NUMTSPLT;TSPLTC=ADJUSTR(TSPLTC) 
                      ELSE 
                                    READ(NUNIT,'(8X,A8,I8)')TSPLTC,NUMTSPLT
                      ENDIF

      
      ALLOCATE(TSYEARLY(NUMTSPLT),TSTSRT(NUMTSPLT),TSTEND(NUMTSPLT),TSPLTJB(NUMTSPLT),TSPLTT(NUMTSPLT),NOUTS(NUMTSPLT),   &
               JSTSPLT(NUMTSPLT,10),KSTRSPLT(NUMTSPLT),TSPLTCNTR(NUMTSPLT))
      ALLOCATE(JSTSPLTT(NUMTSPLT,10),ELCONTSPL(NUMTSPLT),DYNSELSPLT(NUMTSPLT),DYNSPF(NUMTSPLT),NXSELSPLT(NUMTSPLT),TEMP3(NUMTSPLT))
      DYNSPF=.FALSE.
      
      DO J=1,2
      READ(NUNIT,*)
      END DO
      DO J=1,NUMTSPLT
              IF(CHAR1=='$')THEN
                    READ(NUNIT,*)CHAR8,TSPLTCNTR(J),TSPLTJB(J),TSYEARLY(J),TSTSRT(J),TSTEND(J),TSPLTT(J),NOUTS(J),(JSTSPLTT(J,N),N=1,2),ELCONTSPL(J),DYNSELSPLT(J)
                    TSPLTCNTR(J)=ADJUSTR(TSPLTCNTR(J))
                    TSYEARLY(J)=ADJUSTR(TSYEARLY(J))
                    ELCONTSPL(J)=ADJUSTR(ELCONTSPL(J))
                    DYNSELSPLT(J)=ADJUSTR(DYNSELSPLT(J))
              ELSE 
                    READ(NUNIT,'(8X,A8,I8,A8,F8.0,F8.0,F8.0,I8,2I8,A8,A8)')TSPLTCNTR(J),TSPLTJB(J),TSYEARLY(J),TSTSRT(J),TSTEND(J),TSPLTT(J),NOUTS(J),(JSTSPLTT(J,N),N=1,2),ELCONTSPL(J),DYNSELSPLT(J)
              ENDIF
      NOUTS(J)=2                ! NUMBER OF OUTLETS FOR EACH SPLIT FLOW PERIOD LIMITED TO 2
      !IF(NOUTS(J).GT.2)WRITE(*,*)'TCD NOUTS > 2 - ONLY FIRST 2 WILL BE USED'
      ENDDO
      JSTSPLT=JSTSPLTT                                                                                             ! CB 10/14/11 START
      DO J=1,NUMTSPLT  !REODERING OUTLETS SO THAT HIGHEST ELEVATION STRUCTURE ON TOP (ASSUMING 2 SPLIT OUTLETS) 
!        IF(TCNTR(J) == '      ST')THEN
        IF(TSPLTCNTR(J) == '      ST')THEN                                                                        ! cb 11/11/12
          IF(ESTR(JSTSPLTT(J,1),TSPLTJB(J)) < ESTR(JSTSPLTT(J,2),TSPLTJB(J)))THEN                               
            JSTSPLT(J,1)=JSTSPLTT(J,2)                                                                          
            JSTSPLT(J,2)=JSTSPLTT(J,1)                                                                          
          END IF                                                                                                
!        ELSE IF(TCNTR(J) == '      WD')THEN
        ELSE IF(TSPLTCNTR(J) == '      WD')THEN                                                                        ! cb 11/11/12
          IF(EWD(JSTSPLTT(J,1)) < EWD(JSTSPLTT(J,2)))THEN                                    
            JSTSPLT(J,1)=JSTSPLTT(J,2)                                                                          
            JSTSPLT(J,2)=JSTSPLTT(J,1)                                                                          
          END IF                                                                                                
        END IF
      END DO                                                                                                       ! CB 10/14/11 END
      DO J=1,2
      READ(NUNIT,*)
      END DO
               IF(CHAR1=='$')THEN
                    READ(NUNIT,*)CHAR8,TEMPN
              ELSE 
                    READ(NUNIT,'(8X,I8)')TEMPN
              ENDIF
      DO J=1,2
      READ(NUNIT,*)
      END DO
      ALLOCATE(TEMPCRIT(NWB,TEMPN),VOLMC(NWB,TEMPN))
      DO J=1,TEMPN
        IF(CHAR1=='$')THEN
        READ(NUNIT,*)CHAR8,(TEMPCRIT(JW,J),JW=1,NWB)   ! NOTE MAX OF 100 WATERBODIES   sw 4/20/15  
        ELSE
        READ(NUNIT,'(8X,100F8.0)')(TEMPCRIT(JW,J),JW=1,NWB)   ! NOTE MAX OF 100 WATERBODIES   sw 4/20/15
        ENDIF
      END DO
      CLOSE(NUNIT)

      
      DO JW=1,NWB
        IFILE=IFILE+1
        WRITE (SEGNUM,'(I0)') JW
        SEGNUM = ADJUSTL(SEGNUM)
        L      = LEN_TRIM(SEGNUM)
         IF(RESTART_IN)THEN
                OPEN  (IFILE,FILE='VOLUME_WB'//SEGNUM(1:L)//'.OPT',POSITION='APPEND')
                JDAY1=0.0
                REWIND (IFILE)
                READ   (IFILE,'(/)',END=15)
                DO WHILE (JDAY1 < JDAY)
                READ (IFILE,'(F10.0)',END=15) JDAY1
                END DO
                BACKSPACE (IFILE)
                15    JDAY1=0.0   
       ELSE
        OPEN  (IFILE,FILE='VOLUME_WB'//SEGNUM(1:L)//'.OPT',STATUS='UNKNOWN')
        WRITE(IFILE,4315)
       ENDIF
      ENDDO
      
4315  FORMAT("JDAY    VOLUME    ",<TEMPN>("VOLCRIT      "))


! INITIALIZING STRUCTURE ELEVATION IF STRUCTURE
IF(TEMPC=='      ON')THEN     
  DO JW=1,NWB
   DO JB=BS(JW),BE(JW)
    DO JS=1,NST
     DO J=1,NUMTEMPC        
       IF(TCJB(J) == JB .AND. TCJS(J) == JS .AND. TCNTR(J) == '      ST')THEN
           IF(TCYEARLY(J) == '     OFF')THEN
             DAYTEST=JDAY
           ELSE
             DAYTEST=REAL(JDAYG)+JDAY-INT(JDAY)
           END IF
           IF(DAYTEST >= TCTSRT(J) .AND. DAYTEST < TCTEND(J))THEN               
               ! MAKING SURE THAT STRUCTURE IS BELOW WATER SURFACE
             DO NN=1,TCNELEV(J)
               IF(TCELEV(J,NN) < ELWS(DS(JB)))THEN
                 NCOUNTC(JS,JB)=NN
                 ESTR(JS,JB)=TCELEV(J,NCOUNTC(JS,JB))
                 EXIT
               END IF                 
             END DO
		   END IF
	   END IF
	 END DO
	END DO
   END DO
  END DO
   
   
   ! INITIALIZING STRUCTURE ELEVATION IF WITHDRAWAL

  DO JWD=1,NWD
   
     DO J=1,NUMTEMPC        
       IF(TCJS(J) == JWD .AND. TCNTR(J) == '      WD')THEN
           IF(TCYEARLY(J) == '     OFF')THEN
             DAYTEST=JDAY
           ELSE
             DAYTEST=REAL(JDAYG)+JDAY-INT(JDAY)
           END IF
           IF(DAYTEST >= TCTSRT(J) .AND. DAYTEST < TCTEND(J))THEN               
               ! MAKING SURE THAT STRUCTURE IS BELOW WATER SURFACE
             DO NN=1,TCNELEV(J)
               IF(TCELEV(J,NN) < ELWS(IWD(JWD)))THEN
                 NCOUNTCW(JWD)=NN
                 EWD(JWD)=TCELEV(J,NCOUNTCW(JWD))
                 EXIT
               END IF                 
             END DO
		   END IF
	   END IF
	 END DO

  END DO

  
  ! OPEN DYNAMIC SELECTIVE WITHDRAWAL FILES
  
  DO J=1,numtempc
     if(DYNSEL(J) == '      ON')then
     WRITE (SEGNUM,'(I0)') J     
     SEGNUM = ADJUSTL(SEGNUM)
     L      = LEN_TRIM(SEGNUM)   
     SELD(J) = 1009+J  
     OPEN (SELD(J),FILE='dynselective'//SEGNUM(1:L)//'.npt',STATUS='OLD') 
     
      READ(SELD(J),'(A1)')INFORMAT    ! SW 8/28/2019
      IF(INFORMAT=='$')DYNSF(J)=.TRUE.
    
        IF(DYNSF(J))THEN
        READ (SELD(J),'(/)')
        READ (SELD(J),*) NXSEL(J),TEMP2(J)
        tctemp(J)=TEMP2(J)
        READ (SELD(J),*) NXSEL(J),TEMP2(J)       
        
        ELSE
        READ (SELD(J),'(//1000F8.0)') NXSEL(J),TEMP2(J)
        tctemp(J)=TEMP2(J)
        READ (SELD(J),'(1000F8.0)') NXSEL(J),TEMP2(J)
        ENDIF
     
      !READ (SELD(J),'(///1000F8.0)') NXSEL(J),TEMP2(J)
      !  tctemp(J)=TEMP2(J)
      !READ (SELD(J),'(1000F8.0)') NXSEL(J),TEMP2(J)
     END IF
  ENDDO 
END IF

IF(TSPLTC == '      ON')THEN
    DO J=1,NUMTSPLT
     IF(DYNSELSPLT(J) == '      ON')then
     WRITE (SEGNUM,'(I0)') J     
     SEGNUM = ADJUSTL(SEGNUM)
     L      = LEN_TRIM(SEGNUM)   
     SELD(J) = 1059+J  
     OPEN (SELD(J),FILE='dynselective'//SEGNUM(1:L)//'_splt.npt',STATUS='OLD') 
     
      READ(SELD(J),'(A1)')INFORMAT    ! SW 8/28/2019
      IF(INFORMAT=='$')DYNSPF(J)=.TRUE.
    
        IF(DYNSPF(J))THEN
        READ (SELD(J),'(/)')
        READ (SELD(J),*) NXSELSPLT(J),TEMP3(J)
        TSPLTT(J)=TEMP3(J)
        READ (SELD(J),*) NXSELSPLT(J),TEMP3(J)       
        
        ELSE
        READ (SELD(J),'(//1000F8.0)') NXSELSPLT(J),TEMP3(J)
        TSPLTT(J)=TEMP2(J)
        READ (SELD(J),'(1000F8.0)') NXSELSPLT(J),TEMP3(J)
        ENDIF
     
     END IF
    ENDDO 
    OPEN(2900,FILE='Split_Temp_Debug.csv',status='unknown')
        WRITE(2900,*)'JDAY,TSPLTT,TTOP,TBOT,QALL,QSTR1,QSTR2,ESTR1,ESTR2'
ENDIF

  
 RETURN

END SUBROUTINE SELECTIVEINIT

SUBROUTINE SELECTIVE
 USE SELECTIVE1
  USE MAIN
  USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  
  IMPLICIT NONE
  !** Timestep violation entry point  210 CONTINUE                
  INTEGER JJ, JJW, KK, KS, IFILE, KSTR
  REAL DAYTEST, ELR, QALL, TCOMP, TEMPBOT, TEMPEST, TEMPTOP, TMOD, WSEL

  IF(TSPLTC=='      ON')THEN    
    DO J=1,NUMTSPLT
      IF(TSYEARLY(J) == '     OFF')THEN
        DAYTEST=JDAY
      ELSE
        DAYTEST=REAL(JDAYG)+JDAY-INT(JDAY)
      END IF
      IF(NXTSPLIT > TSTSRT(J) .AND. DAYTEST <= TSTSRT(J))THEN
        NXTSPLIT=TSTSRT(J)
      END IF


  IF(DYNSELSPLT(J) == '      ON')THEN
    SELD(J)=1059+J
     DO WHILE (JDAY >= NXSELSPLT(J))
        TSPLTT(J)=TEMP3(J)
                IF(DYNSPF(J))THEN
                    READ (SELD(J),*) NXSELSPLT(J),TEMP3(J)
                ELSE
                    READ (SELD(J),'(1000F8.0)') NXSELSPLT(J),TEMP3(J)
                ENDIF   
    END DO
   ENDIF
   ENDDO
END IF

 IF(TSPLTC=='      ON'.AND.JDAY.GE.NXTSPLIT)THEN  
 
  DO J=1,NUMTSPLT
        IF(TSYEARLY(J) == '     OFF')THEN
            DAYTEST=JDAY
          ELSE
        DAYTEST=REAL(JDAYG)+JDAY-INT(JDAY)
        END IF
   IF(DAYTEST >= TSTSRT(J) .AND. DAYTEST < TSTEND(J))THEN 
    ! DO STRUCTURES FIRST
    DO JW=1,NWB
        DO JB=BS(JW),BE(JW)
            IF(TSPLTJB(J) == JB .AND. TSPLTCNTR(J) == '      ST')THEN
                QALL=0.0
                DO JJ=1,NOUTS(J)
                QALL=QALL+QSTR(JSTSPLT(J,JJ),TSPLTJB(J))   ! SUM UP ALL THE FLOWS
                ELR  = SINA(JB)*DLX(DS(JB))*0.5
                    DO K=KTWB(JW),KB(DS(JB))
                    IF (EL(K,DS(JB))-ELR < ESTR(JSTSPLT(J,JJ),TSPLTJB(J))) EXIT                                                                               !SW 10/17/01
                    END DO
                KSTR = K-1
                KSTRSPLT(JJ) = MIN(KSTR,KB(DS(JB)))
                ENDDO               
              DO JJ=1,NOUTS(J)               ! cb 11/11/12 dividing total flow between outlets for temperature test - if no flow there is no temperature test
                  QSTR(JSTSPLT(J,JJ),TSPLTJB(J)) = qall/real(nouts(j))
              ENDDO               
              ID=DS(JB)
              ELR  = SINA(JB)*DLX(ID)*0.5          ! CB 10/14/11
              WSEL = ELWS(ID)-ELR                  ! CB 10/14/11
              kt=ktwb(jw)      ! cb 07/24/19
              CALL DOWNSTREAM_WITHDRAWAL_ESTIMATE(JSTSPLT(J,1),TEMPTOP,ESTR(JSTSPLT(J,1),TSPLTJB(J)))
              CALL DOWNSTREAM_WITHDRAWAL_ESTIMATE(JSTSPLT(J,2),TEMPBOT,ESTR(JSTSPLT(J,2),TSPLTJB(J)))
             IF(ESTR(JSTSPLT(J,1),TSPLTJB(J)) > WSEL .AND. ELCONTSPL(J) =='     OFF') THEN   ! NO FLOWS THROUG THIS OUTLET IF WSEL BELOW LEVEL OF OUTLET  ! CB 10/14/11
               QSTR(JSTSPLT(J,1),TSPLTJB(J))=0.0
               QSTRFRAC(JSTSPLT(J,1),TSPLTJB(J))=0.0
              
            ELSE IF(TEMPTOP > TSPLTT(J)  .AND.  TEMPBOT > TSPLTT(J) ) THEN   ! NO FLOWS THROUG THIS OUTLET IF T1 AND T2 > TCRITERIA
               QSTR(JSTSPLT(J,1),TSPLTJB(J))=0.0
               QSTRFRAC(JSTSPLT(J,1),TSPLTJB(J))=0.0

              !ELSEIF(T2(KSTRSPLT(1),DS(JB)) < TSPLTT(J)) THEN   ! ALL FLOWS FROM TOP IF TCRITERIA < TOUTLET
              ELSEIF(TEMPTOP < TSPLTT(J)) THEN   ! ALL FLOWS FROM TOP IF TCRITERIA < TOUTLET
               QSTR(JSTSPLT(J,1),TSPLTJB(J))=QALL
               QSTRFRAC(JSTSPLT(J,1),TSPLTJB(J))=1.0

              ELSE
                !QSTR(JSTSPLT(J,1),TSPLTJB(J))=QALL*(TSPLTT(J)-T2(KSTRSPLT(2),DS(JB)))/(T2(KSTRSPLT(1),DS(JB))-T2(KSTRSPLT(2),DS(JB)))
                IF(ABS(TEMPTOP-TEMPBOT) < 0.0001)THEN
                  QSTR(JSTSPLT(J,1),TSPLTJB(J))=QALL
                  QSTRFRAC(JSTSPLT(J,1),TSPLTJB(J))=1.0
                ELSE
                  QSTR(JSTSPLT(J,1),TSPLTJB(J))=QALL*(TSPLTT(J)-TEMPBOT)/(TEMPTOP-TEMPBOT)
                  QSTRFRAC(JSTSPLT(J,1),TSPLTJB(J))=QSTR(JSTSPLT(J,1),TSPLTJB(J))/QALL
                END IF
              ENDIF

              QSTR(JSTSPLT(J,2),TSPLTJB(J))=QALL-QSTR(JSTSPLT(J,1),TSPLTJB(J))
              QSTRFRAC(JSTSPLT(J,2),TSPLTJB(J))=QSTR(JSTSPLT(J,2),TSPLTJB(J))/QALL
              write(2900,'(9(f12.4,","))')jday,tspltt(j),TEMPTOP,TEMPBOT,QALL,QSTR(JSTSPLT(J,1),TSPLTJB(J)),QSTR(JSTSPLT(J,2),TSPLTJB(J)),ESTR(JSTSPLT(J,1),TSPLTJB(J)),ESTR(JSTSPLT(J,2),TSPLTJB(J))
             EXIT
             END IF
         END DO
       END DO
    ! DO WITHDRAWALS NEXT
      DO JWD=1,NWD
            IF(TSPLTCNTR(J) == '      WD')THEN
                QALL=0.0
               DO JJB=1,NBR
                 IF (IWD(JWD) >= US(JJB) .AND. IWD(JWD) <= DS(JJB)) EXIT
               END DO
               DO JJW=1,NWB
                 IF (JJB >= BS(JJW) .AND. JJB <= BE(JJW)) EXIT
               END DO
               DO JJ=1,NOUTS(J)
                QALL=QALL+QWD(JSTSPLT(J,JJ))   ! SUM UP ALL THE FLOWS
                ELR  = SINA(JJB)*DLX(IWD(JWD))*0.5
                    DO K=KTWB(JJW),KB(IWD(JWD))
                    IF (EL(K,IWD(JWD))-ELR < EWD(JSTSPLT(J,JJ))) EXIT                                                                               !SW 10/17/01
                    END DO
                KSTR = K-1
                KSTRSPLT(JJ) = MIN(KSTR,KB(IWD(JWD)))
               ENDDO
               JJ=1               ! ASSIGN FLOW TO FIRST OUTLET               
               WSEL = ELWS(IWD(JWD))-ELR                  ! CB 10/14/11
               I=IWD(JWD)
               kt=ktwb(jjw)      ! cb 07/24/19
               CALL LATERAL_WITHDRAWAL_ESTIMATE(JSTSPLT(J,1),TEMPTOP,EWD(JSTSPLT(J,1)))
               CALL LATERAL_WITHDRAWAL_ESTIMATE(JSTSPLT(J,2),TEMPBOT,EWD(JSTSPLT(J,2)))              
              IF(EWD(JSTSPLT(J,1)) > WSEL .AND. TCELEVCON(J) =='     OFF') THEN
                QWD(JSTSPLT(J,1))=0.0
                QWDFRAC(JSTSPLT(J,1))=0.0             
             ELSE IF(TEMPTOP > TSPLTT(J)  .AND.  TEMPBOT > TSPLTT(J) ) THEN   ! NO FLOWS THROUG THIS OUTLET IF T1 AND T2 > TCRITERIA
               QWD(JSTSPLT(J,1))=0.0
               QWDFRAC(JSTSPLT(J,1))=0.0

              ELSEIF(TEMPTOP < TSPLTT(J)) THEN   ! ALL FLOWS FROM TOP IF TCRITERIA < TOUTLET
               QWD(JSTSPLT(J,1))=QALL
               QWDFRAC(JSTSPLT(J,1))=1.0

              ELSE
                !QWD(JSTSPLT(J,1))=QALL*(TSPLTT(J)-T2(KSTRSPLT(2),IWD(JWD)))/(T2(KSTRSPLT(1),IWD(JWD))-T2(KSTRSPLT(2),IWD(JWD)))
                IF(ABS(TEMPTOP-TEMPBOT) < 0.0001)THEN
                  QWD(JSTSPLT(J,1))=QALL
                  QWDFRAC(JSTSPLT(J,1))=1.0
                ELSE
                  QWD(JSTSPLT(J,1))=QALL*(TSPLTT(J)-TEMPBOT)/(TEMPTOP-TEMPBOT)
                  QWDFRAC(JSTSPLT(J,1))=QWD(JSTSPLT(J,1))/QALL
                END IF
              ENDIF

              QWD(JSTSPLT(J,2))=QALL-QWD(JSTSPLT(J,1))
              QWDFRAC(JSTSPLT(J,2))=QWD(JSTSPLT(J,2))/QALL
             EXIT
             END IF
       END DO
     ENDIF
     ENDDO
  
   NXTSPLIT=NXTSPLIT+TCDFREQ
  END IF
  IF(TSPLTC=='      ON')THEN

    DO J=1,NUMTSPLT
    IF(TSYEARLY(J) == '     OFF')THEN
            DAYTEST=JDAY
          ELSE
        DAYTEST=REAL(JDAYG)+JDAY-INT(JDAY)
        END IF
    IF(DAYTEST >= TSTSRT(J) .AND. DAYTEST < TSTEND(J))THEN 
    ! DO STRUCTURES FIRST
      DO JW=1,NWB
        DO JB=BS(JW),BE(JW)
            IF(TSPLTJB(J) == JB .AND. TSPLTCNTR(J) == '      ST')THEN
                QALL=0.0
                DO JJ=1,NOUTS(J)
                QALL=QALL+QSTR(JSTSPLT(J,JJ),TSPLTJB(J))   ! SUM UP ALL THE FLOWS
                ELR  = SINA(JB)*DLX(DS(JB))*0.5
                    DO K=KTWB(JW),KB(DS(JB))
                    IF (EL(K,DS(JB))-ELR < ESTR(JSTSPLT(J,JJ),TSPLTJB(J))) EXIT                                                                               !SW 10/17/01
                    END DO
                KSTR = K-1
                KSTRSPLT(JJ) = MIN(KSTR,KB(DS(JB)))
                ENDDO               
              QSTR(JSTSPLT(J,1),TSPLTJB(J))=QSTRFRAC(JSTSPLT(J,1),TSPLTJB(J))*QALL
              QSTR(JSTSPLT(J,2),TSPLTJB(J))=QSTRFRAC(JSTSPLT(J,2),TSPLTJB(J))*QALL
             EXIT
             END IF
        END DO
      END DO
    ! DO WITHDRAWALS NEXT
      DO JWD=1,NWD
            IF(TSPLTCNTR(J) == '      WD')THEN
                QALL=0.0
                DO JJB=1,NBR
                  IF (IWD(JWD) >= US(JJB) .AND. IWD(JWD) <= DS(JJB)) EXIT
                END DO
                DO JJW=1,NWB
                  IF (JJB >= BS(JJW) .AND. JJB <= BE(JJW)) EXIT
                END DO
               DO JJ=1,NOUTS(J)
                QALL=QALL+QWD(JSTSPLT(J,JJ))   ! SUM UP ALL THE FLOWS
                ELR  = SINA(JJB)*DLX(IWD(JWD))*0.5
                    DO K=KTWB(JJW),KB(IWD(JWD))
                    IF (EL(K,IWD(JWD))-ELR < EWD(JSTSPLT(J,JJ))) EXIT                                                                               !SW 10/17/01
                    END DO
                KSTR = K-1
                KSTRSPLT(JJ) = MIN(KSTR,KB(IWD(JWD)))
               ENDDO               
               QWD(JSTSPLT(J,1))=  QWDFRAC(JSTSPLT(J,1))*QALL
               QWD(JSTSPLT(J,2))=  QWDFRAC(JSTSPLT(J,2))*QALL
             EXIT
             END IF
      END DO
      ENDIF
    ENDDO
   
  ENDIF   
 

      IF (JDAY.GE.NXTSTR) THEN
        NXTSTR = NXTSTR+TFRQTMP   
        IFILE=1949
        DO JB=1,NBR
            IF(NSTR(JB) > 0)THEN
            IFILE=IFILE+1
            WRITE (IFILE,'(F10.4,",",<NSTR(JB)>(F10.2,","),<NSTR(JB)>(F10.2,","),<NSTR(JB)>(F10.2,","))') JDAY,(TAVG(I,JB),I=1,NSTR(JB)),(QSTR(I,JB),I=1,NSTR(JB)),(ESTR(I,JB),I=1,NSTR(JB))
            END IF
         ENDDO          
          IF(NWD > 0)THEN
            IFILE=IFILE+1
            WRITE (IFILE,'(F10.4,<NWD>F10.2,<NWD>F10.2,<NWD>F10.2)') JDAY,(TAVGW(I),I=1,NWD),(QWD(I),I=1,NWD),(EWD(I),I=1,NWD)
          END IF
         ! TEMPERATURE CONTROL LOGIC 

         ! COMPUTING RESERVOIR VOLUME AND VOLUME BELOW 'TEMPCRIT'        
        VOLMC=0.0
        VOLM=0.0
        DO JW=1,NWB
         KT = KTWB(JW)
           DO JB=BS(JW),BE(JW)           
             DO I=CUS(JB),DS(JB)
               VOLM(JW) = VOLM(JW) +BH2(KT,I)*DLX(I)               
               DO K=KT+1,KB(I)
                 VOLM(JW) = VOLM(JW)+BH(K,I)*DLX(I)               
               END DO
               DO KK=1,TEMPN                                         
                 IF(T2(KT,I).LE.TEMPCRIT(JW,KK))VOLMC(JW,KK) = VOLMC(JW,KK)+BH2(KT,I)*DLX(I)                                                 
                 DO K=KT+1,KB(I)                 
                   IF(T2(K,I).LE.TEMPCRIT(JW,KK))VOLMC(JW,KK) = VOLMC(JW,KK)+BH(K,I)*DLX(I)
                 END DO
               END DO               
             END DO         
           END DO
     
         IFILE=IFILE+1
         WRITE(IFILE,5315)JDAY,VOLM(JW),(VOLMC(JW,KK), KK=1,TEMPN)
5315     FORMAT(F8.2,100(G12.4,G12.4))
       ENDDO
         
      ENDIF



IF(TEMPC=='      ON'.AND.JDAY.GE.NXTTCD)THEN  

! IF DYNAMIC SELECTIVE CHANGE TEMPERATURE 

 DO J=1,NUMTEMPC
  IF(DYNSEL(J) == '      ON')THEN
    SELD(J)=1009+J
     DO WHILE (JDAY >= NXSEL(J))
        TCTEMP(J)=TEMP2(J)
                IF(DYNSF(J))THEN
                    READ (SELD(J),*) NXSEL(J),TEMP2(J)
                ELSE
                    READ (SELD(J),'(1000F8.0)') NXSEL(J),TEMP2(J)
                ENDIF   
    END DO
   ENDIF
  ENDDO

  
! STRUCTURES  
  
  DO JW=1,NWB
   DO JB=BS(JW),BE(JW)
    DO JS=1,NST
     DO J=1,NUMTEMPC      

          
      IF(TCJB(J) == JB .AND. TCJS(J) == JS .AND.  TCNTR(J) == '      ST')THEN
          IF(TCISEG(J).EQ.0)THEN
            TCOMP=TAVG(TCJS(J),TCJB(J))   !CB 9/8/06   TAVG(JSMON(J),JBMON(J))
          ELSEIF(TCISEG(J) < 0)THEN
            TCOMP=TWDO(ABS(TCISEG(J)))      ! SW 11/26/10       
          ELSE

! CHECKING TO SEE IF THE MONITORING SEGMENT TCISEG IS IN THE SAME BRANCH AND WATER BODY AS THE STRUCTURE
            DO JJB=1,NBR
              IF (TCISEG(J) >= US(JJB) .AND. TCISEG(J) <= DS(JJB)) EXIT
            END DO
            DO JJW=1,NWB
              IF (JJB >= BS(JJW) .AND. JJB <= BE(JJW)) EXIT
            END DO

            IF (TCKLAY(J)< 0) THEN                                                                                       
              K = INT(ABS(TCKLAY(J)))
            ELSE
              DO K=KTWB(JJW),KB(TCISEG(J))
                IF (DEPTHB(K,TCISEG(J)) > TCKLAY(J)) EXIT                                                                      
              END DO
              K = MIN(K,KB(TCISEG(J)))                                                                                         
            END IF
            TCOMP=T2(K,TCISEG(J))
          ENDIF
          IF(TCYEARLY(J) == '     OFF')THEN
            DAYTEST=JDAY
          ELSE
            DAYTEST=REAL(JDAYG)+JDAY-INT(JDAY)
          END IF          
          IF(DAYTEST >= TCTSRT(J) .AND. DAYTEST < TCTEND(J))THEN
               ACTIVE_RULE_W2SELECTIVE(JS,JB)=.TRUE.
               IF(TCOMP > TCTEMP(J) .AND. TCNELEV(J) > NCOUNTC(JS,JB))THEN
               ! MAKING SURE THAT THE NEXT LOWER STRUCTURE FOR A PARTICULAR 'J' IS FOUND
                 DO NN=NCOUNTC(JS,JB)+1,TCNELEV(J)                 
                   IF(TCELEV(J,NN) < ESTR(JS,JB))THEN
                      NCOUNTC(JS,JB)=NN
                      ESTR(JS,JB)=TCELEV(J,NCOUNTC(JS,JB))
                      EXIT
                   END IF                 
                 END DO                                               
               ELSEIF(TCOMP < TCTEMP(J) .AND.  NCOUNTC(JS,JB).GT. 1)THEN
                 ! TO PREVENT THIS HAPPENING AT EACH TIME IT CHECKS IT AND HENCE OSCIALLTING BACK AND FORTH - CHECK THE TEMP AT THE UPPER OUTLET ALSO
                 IF(TCISEG(J) > 0)THEN  
                    IF(JB.EQ.JJB)THEN
                      DO KS=KTWB(JW),KB(DS(JB))
                        IF (DEPTHB(KS,TCISEG(J)) > TCELEV(J,NCOUNTC(JS,JB)-1)) EXIT                                                                          !TC 01/03/02
                      END DO
                      KS = MIN(KS,KB(TCISEG(J)))
                      TMOD= T2(KS,DS(JB))
                    ELSE
                      TMOD=T2(K,TCISEG(J))                      
                    END IF
                    IF(TMOD < TCTEMP(J) .AND. TCELEV(J,NCOUNTC(JS,JB)-1) < ELWS(DS(JB)))THEN                      
                      ! MAKING SURE THAT THE NEXT UPPER STRUCTURE FOR A PARTICULAR 'J' IS FOUND
                      DO NN=NCOUNTC(JS,JB)-1,1,-1
                        IF(TCELEV(J,NN) > ESTR(JS,JB))THEN
                          NCOUNTC(JS,JB)=NN
                          ESTR(JS,JB)=TCELEV(J,NCOUNTC(JS,JB))
                          EXIT
                        END IF                 
                      END DO
                    ENDIF
                 END IF  ! CB 9/8/06
                 IF(TCISEG(J).EQ.0)THEN
! CALCULATE THE ESTIMATED OUTFLOW TEMPERATURE AT HIGHER PORTS WHEN TCOMP<TCTEMP(J), AND MOVE UP IF HIGHER PORT STILL MEETS TO CRITERIA - THIS DOESN'T HAPPEN WHEN TCISEG < 0
                   DO NN=1,NCOUNTC(JS,JB)-1
                     ID=DS(JB)
                     KT=KTWB(JW)                     
                     CALL DOWNSTREAM_WITHDRAWAL_ESTIMATE(JS,TEMPEST,TCELEV(J,NN))
                     IF(TEMPEST < TCTEMP(J) .AND. TCELEV(J,NN) < ELWS(DS(JB)))THEN
                       NCOUNTC(JS,JB)=NN
                       ESTR(JS,JB)=TCELEV(J,NCOUNTC(JS,JB))
                       EXIT
                     END IF
                   END DO
                 END IF
               ENDIF
               IF(TCELEVCON(J) =='      ON' .AND. TCNELEV(J) > NCOUNTC(JS,JB).AND. ESTR(JS,JB) > (ELWS(DS(JB))-MINWL(J)))THEN  
                 NCOUNTC(JS,JB)=NCOUNTC(JS,JB)+1
                 ESTR(JS,JB)=TCELEV(J,NCOUNTC(JS,JB))
               END IF
          ELSE
              ACTIVE_RULE_W2SELECTIVE(JS,JB)=.FALSE.
          ENDIF          
      ENDIF            
     END DO
    END DO
   END DO
  ENDDO
  
 ! Withdrawals 
  

  DO JWD=1,NWD
     DO J=1,NUMTEMPC
      IF(TCJS(J) == JWD .AND.  TCNTR(J) == '      WD')THEN
          IF(TCISEG(J).EQ.0)THEN
!           TCOMP=TOUT(JB)
            TCOMP=TAVGW(TCJS(J))   !CB 9/8/06   TAVGW(JSMON(J))
          ELSEIF(TCISEG(J) < 0)THEN  
            TCOMP=TWDO(ABS(TCISEG(J)))
          ELSE

! CHECKING TO SEE IF THE MONITORING SEGMENT TCISEG IS IN THE SAME BRANCH AND WATER BODY AS THE WITHDRAWAL
            DO JJB=1,NBR
              IF (TCISEG(J) >= US(JJB) .AND. TCISEG(J) <= DS(JJB)) EXIT
            END DO
            DO JJW=1,NWB
              IF (JJB >= BS(JJW) .AND. JJB <= BE(JJW)) EXIT
            END DO

            IF (TCKLAY(J)< 0) THEN                                                                                       
              K = INT(ABS(TCKLAY(J)))
            ELSE
              DO K=KTWB(JJW),KB(TCISEG(J))
                IF (DEPTHB(K,TCISEG(J)) > TCKLAY(J)) EXIT                                                                      
              END DO
              K = MIN(K,KB(TCISEG(J)))                                                                                         
            END IF
            TCOMP=T2(K,TCISEG(J))
          ENDIF
          IF(TCYEARLY(J) == '     OFF')THEN
            DAYTEST=JDAY
          ELSE
            DAYTEST=REAL(JDAYG)+JDAY-INT(JDAY)
          END IF
          IF(DAYTEST >= TCTSRT(J) .AND. DAYTEST < TCTEND(J))THEN
               IF(TCOMP > TCTEMP(J) .AND. TCNELEV(J) > NCOUNTCW(JWD))THEN
               ! MAKING SURE THAT THE NEXT LOWER STRUCTURE FOR A PARTICULAR 'J' IS FOUND
                 DO NN=NCOUNTCW(JWD)+1,TCNELEV(J)                 
                   IF(TCELEV(J,NN) < EWD(JWD))THEN
                      NCOUNTCW(JWD)=NN
                      EWD(JWD)=TCELEV(J,NCOUNTCW(JWD))
                      EXIT
                   END IF                 
                 END DO                                               
               ELSEIF(TCOMP < TCTEMP(J) .AND.  NCOUNTCW(JWD).GT. 1)THEN
                 ! TO PREVENT THIS HAPPENING AT EACH TIME IT CHECKS IT AND HENCE OSCIALLTING BACK AND FORTH - CHECK THE TEMP AT THE UPPER OUTLET ALSO
                 IF(TCISEG(J) >  0)THEN  
                      TMOD=T2(K,TCISEG(J))                      
                    IF(TMOD < TCTEMP(J) .AND. TCELEV(J,NCOUNTCW(JWD)-1) < ELWS(IWD(JWD)))THEN                      
                      ! MAKING SURE THAT THE NEXT UPPER STRUCTURE FOR A PARTICULAR 'J' IS FOUND
                      DO NN=NCOUNTCW(JWD)-1,1,-1
                        IF(TCELEV(J,NN) > EWD(JWD))THEN
                          NCOUNTCW(JWD)=NN
                          EWD(JWD)=TCELEV(J,NCOUNTCW(JWD))
                          EXIT
                        END IF                 
                      END DO
                    ENDIF
                 END IF  ! CB 9/8/06
                 IF(TCISEG(J) == 0)THEN
! CALCULATE THE ESTIMATED OUTFLOW TEMPERATURE AT HIGHER PORTS WHEN TCOMP<TCTEMP(J), AND MOVE UP IF HIGHER PORT STILL MEETS TO CRITERIA
                   I         = MAX(CUS(JBWD(JWD)),IWD(JWD))
                   DO NN=1,NCOUNTCW(JWD)-1                     
                     CALL LATERAL_WITHDRAWAL_ESTIMATE(JWD,TEMPEST,TCELEV(J,NN))
                     IF(TEMPEST < TCTEMP(J) .AND. TCELEV(J,NN) < ELWS(IWD(JWD)))THEN
                       NCOUNTCW(JWD)=NN
                       EWD(JWD)=TCELEV(J,NCOUNTCW(JWD))
                       EXIT
                     END IF
                   END DO
                 END IF
               ENDIF
               IF(TCELEVCON(J) =='      ON' .AND. TCNELEV(J) > NCOUNTCW(JWD).AND. EWD(JWD) > ELWS(IWD(JWD)))THEN  
                 NCOUNTCW(JWD)=NCOUNTCW(JWD)+1
                 EWD(JWD)=TCELEV(J,NCOUNTCW(JWD))
               END IF
          ENDIF          
      ENDIF            
     END DO
  ENDDO
  
  NXTTCD = NXTTCD+TCDFREQ    
ENDIF  
  
RETURN
ENTRY DEALLOCATE_SELECTIVE
CLOSE(2900)
  DEALLOCATE (TCNELEV,TCJB,TCJS, TCELEV,TCTEMP,TCTEND,TCTSRT,NCOUNTC,TCISEG,TCKLAY,TCELEVCON,ELCONTSPL) 
  DEALLOCATE (TSPLTJB,TSPLTT,NOUTS,JSTSPLT,KSTRSPLT,TCYEARLY, JBMON,JSMON,TCNTR,TSPLTCNTR,JSTSPLTT) 
  DEALLOCATE (VOLM,MONCTR,NCOUNTCW,QWDFRAC,QSTRFRAC,MINWL)    
  DEALLOCATE(TEMPCRIT,VOLMC,DYNSEL,SELD,NXSEL,TEMP2,TSYEARLY,TSTEND,TSTSRT,DYNSF,DYNSPF,DYNSELSPLT,TEMP3,NXSELSPLT)
RETURN   

END SUBROUTINE SELECTIVE
!***********************************************************************************************************************************
!**                                                S E L E C T I V E   I N I T                                                    **
!***********************************************************************************************************************************

Module Selective1USGS
  INTEGER                                      :: numtempc, numtsplt, tempn, ng1, ng2
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: tcnelev, tcjb, tcjs, tciseg, kstrsplt, ncountcw, seld
  INTEGER,       ALLOCATABLE, DIMENSION(:)     :: tspltjb, tsseld, nouts, nout0, nout1, nout2
  INTEGER,       ALLOCATABLE, DIMENSION(:,:)   :: jstsplt, ncountc, tsprior
  REAL                                         :: nxtstr, nxttcd, nxtsplit, tcdfreq, tfrqtmp, tspltfreq
  REAL                                         :: sum_minfrac1, sum_minfrac2, qfrac1, qfrac2, tsconv
  REAL,          ALLOCATABLE, DIMENSION(:)     :: tctemp, tctend, tctsrt, tcklay, tspltt, volm, qwdfrac, tstend, tstsrt, nxsel,temp2
  REAL,          ALLOCATABLE, DIMENSION(:)     :: nxtssel, tstemp2, ewdsav, minfrac1, minfrac2, maxfrac1, maxfrac2, splt2t, splt2e
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: tcelev, tempcrit, qstrfrac, volmc
  REAL,          ALLOCATABLE, DIMENSION(:,:)   :: tsdepth, tsminfrac, tsminhead, tsmaxhead, tsmaxflow, estrsav
  CHARACTER(8)                                 :: tempc, tspltc
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: tcelevcon, tcyearly, tcntr, dynsel
  CHARACTER(8),  ALLOCATABLE, DIMENSION(:)     :: tspltcntr, tsyearly, elcontspl, tsdynsel
  CHARACTER(5),  ALLOCATABLE, DIMENSION(:,:)   :: tstype
  LOGICAL,       ALLOCATABLE, DIMENSION(:)     :: wd_active, share_flow
  LOGICAL,       ALLOCATABLE, DIMENSION(:,:)   :: no_flow, str_active
  LOGICAL, ALLOCATABLE, DIMENSION(:)           :: DYNSF

End Module Selective1USGS


Subroutine SelectiveInitUSGS

  Use Selective1USGS; USE MAIN
  USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  IMPLICIT NONE
  
  integer      :: ifile, nj, N, JJ
  REAL         :: DAYTEST
  character(8) :: tsshare
  CHARACTER(1) INFORMAT

  ifile=1949
  tavg=0.0
  tavgw=0.0
  do jb=1,nbr
    if(nstr(jb) > 0)then
      ifile=ifile+1
      write (segnum,'(i0)') jb
      segnum = adjustl(segnum)
      l      = len_trim(segnum)

      IF(RESTART_IN)THEN
        OPEN  (ifile,file='str_br'//segnum(1:l)//'.csv',POSITION='APPEND')
        JDAY1=0.0
        REWIND (IFILE)
        READ   (IFILE,'(//)',END=13)
        DO WHILE (JDAY1 < JDAY)
          READ (IFILE,*,END=13) JDAY1                !'(F10.0)'
        END DO
        BACKSPACE (IFILE)
13      JDAY1=0.0
      ELSE
        open  (ifile,file='str_br'//segnum(1:l)//'.csv',status='unknown')
        write (ifile,*)'Branch:,',jb,', # of structures:,',nstr(jb),', outlet temperatures'
        write (ifile,'("      JDAY,",<nstr(jb)>(6x,"T(C),"),<nstr(jb)>(3x,"Q(m3/s),"),<nstr(jb)>(4x,"ELEVCL,"))')
      ENDIF
    endif
  end do

  if(nwd > 0)then
    ifile=ifile+1
    IF(RESTART_IN)THEN
      OPEN  (ifile,file='wd_out.opt',POSITION='APPEND')
      JDAY1=0.0
      REWIND (IFILE)
      READ   (IFILE,'(/)',END=14)
      DO WHILE (JDAY1 < JDAY)
        READ (IFILE,'(F10.0)',END=14) JDAY1
      END DO
      BACKSPACE (IFILE)
14    JDAY1=0.0
    ELSE
      open  (ifile,file='wd_out.opt',status='unknown')
      write (ifile,*)'Withdrawals: # of withdrawals:',nwd,' outlet temperatures'
      write (ifile,'("      JDAY",<nwd>(6x,"T(C)"),<nwd>(3x,"Q(m3/s)"),<nwd>(4x,"ELEVCL"))')
    ENDIF
  end if

  open (NUNIT,file='w2_selective.npt',status='old')
  read (NUNIT,'(///8x,f8.0)') tfrqtmp
  read (NUNIT,'(//8x,a8,i8,f8.0)') tempc, numtempc, tcdfreq
  nxtstr   = tmstrt
  nxttcd   = tmstrt
  nxtsplit = tmstrt

  allocate (tcnelev(numtempc),tcjb(numtempc),tcjs(numtempc),tcelev(numtempc,11),tctemp(numtempc),tctend(numtempc),tctsrt(numtempc))
  allocate (ncountc(nst,nbr),tciseg(numtempc),tcklay(numtempc),tcelevcon(numtempc))
  Allocate (tcyearly(numtempc), tcntr(numtempc))
  allocate (volm(nwb),ncountcw(nwd),qwdfrac(nwd),qstrfrac(nst,nbr),DYNSEL(numtempc),SELD(numtempc),NXSEL(numtempc),TEMP2(numtempc))
  ALLOCATE(DYNSF(NUMTEMPC))
  DYNSF=.FALSE.
  
  ncountc=0
  read (NUNIT,'(/)')
  do j=1,numtempc
    read(NUNIT,'(8x,a8,i8,i8,a8,f8.0,f8.0,f8.0,i8,10(f8.0))')tcntr(j),tcjb(j),tcjs(j),tcyearly(j),tctsrt(j),tctend(j),tctemp(j),tcnelev(j),(tcelev(j,n),n=1,tcnelev(j))
    if(tcntr(j)=='      ST')then
      tcelev(j,tcnelev(j)+1)=ESTR(tcjs(j),tcjb(j))   ! always put the original elevation as the last elevation
    else
      tcelev(j,tcnelev(j)+1)=EWD(tcjs(j))   ! always put the original elevation as the last elevation
    endif
  end do
  read (NUNIT,'(/)')
  do j=1,numtempc
    read (NUNIT,'(8x,i8,f8.0,A8)') tciseg(j), tcklay(j), DYNSEL(J)
  end do
  read (NUNIT,'(/)')
  do j=1,numtempc
    read (NUNIT,'(8x,a8)') tcelevcon(j)
  end do
  read (NUNIT,'(//8x,a8,i8,2f8.0)') tspltc, numtsplt, tspltfreq, tsconv

  allocate (tsyearly(numtsplt), tstsrt(numtsplt), tstend(numtsplt), tspltjb(numtsplt), tspltt(numtsplt), nouts(numtsplt))
  allocate (jstsplt(numtsplt,10), kstrsplt(10), tspltcntr(numtsplt), elcontspl(numtsplt))
  allocate (tsdepth(numtsplt,10), tstype(numtsplt,10), tsminfrac(numtsplt,10), tsprior(numtsplt,10))
  allocate (tsminhead(numtsplt,10), tsmaxhead(numtsplt,10), tsmaxflow(numtsplt,10), no_flow(numtsplt,10), share_flow(numtsplt))
  allocate (tsdynsel(numtsplt), tsseld(numtsplt), nxtssel(numtsplt), tstemp2(numtsplt))
  allocate (nout0(10), nout1(10), nout2(10), minfrac1(10), maxfrac1(10), minfrac2(10), maxfrac2(10), splt2t(10), splt2e(10))
  allocate (ewdsav(nwd), wd_active(nwd), estrsav(nst,nbr), str_active(nst,nbr))

  read (NUNIT,'(/)')
  do j=1,numtsplt
    read (NUNIT,'(8x,a8,i8,a8,3f8.0,2a8,i8,a8)') tspltcntr(j), tspltjb(j), tsyearly(j), tstsrt(j), tstend(j),                       &
                                                tspltt(j), tsdynsel(j), elcontspl(j), nouts(j), tsshare
    if (tspltc == '      ON') then
      if (nouts(j) < 2) then
        write (w2err, '(A,I0)') 'ERROR-- Less than two outlets specified for blending group ',j
        ERROR_OPEN = .TRUE.     ! will trigger the program to end when this subroutine is completed
        return
      else if (nouts(j) > 10) then
        write (w2err, '(A,I0)') 'ERROR-- More than ten outlets specified for blending group ',j
        ERROR_OPEN = .TRUE.
        return
      end if
    end if
    share_flow(j) = tsshare == '      ON'
  end do
  read (NUNIT,'(/)')
  do j=1,numtsplt
    read (NUNIT,'(8x,10i8)') (jstsplt(j,n),n=1,nouts(j))
  end do
  read (NUNIT,'(/)')
  do j=1,numtsplt
    read (NUNIT,'(8x,10f8.0)') (tsdepth(j,n),n=1,nouts(j))
  end do
  read (NUNIT,'(/)')
  do j=1,numtsplt
    read (NUNIT,'(8x,10f8.0)') (tsminfrac(j,n),n=1,nouts(j))
  end do
  read (NUNIT,'(/)')
  do j=1,numtsplt
    read (NUNIT,'(8x,10i8)') (tsprior(j,n),n=1,nouts(j))
  end do
  read (NUNIT,'(/)')
  do j=1,numtsplt
    read (NUNIT,'(8x,10f8.0)') (tsminhead(j,n),n=1,nouts(j))
  end do
  read (NUNIT,'(/)')
  do j=1,numtsplt
    read (NUNIT,'(8x,10f8.0)') (tsmaxhead(j,n),n=1,nouts(j))
  end do
  read (NUNIT,'(/)')
  do j=1,numtsplt
    read (NUNIT,'(8x,10f8.0)') (tsmaxflow(j,n),n=1,nouts(j))
  end do

  estrsav = estr    ! Save the original structure elevations
  ewdsav  = ewd     ! Save the original withdrawal elevations
  do j=1,numtsplt
    do n=1,nouts(j)
      tstype(j,n) = "FIXED"
      if (tsdepth(j,n)   > 0.0) tstype(j,n) = "FLOAT"
      if (tsminfrac(j,n) > 1.0) tsminfrac(j,n) = 1.0    ! remove unrealistic input value
      if (tsminhead(j,n) < 0.0) tsminhead(j,n) = 0.0    ! remove unrealistic input value
      if (tsmaxhead(j,n) < 0.0) tsmaxhead(j,n) = 0.0    ! remove unrealistic input value
      if (tsmaxflow(j,n) < 0.0) tsmaxflow(j,n) = 0.0    ! remove unrealistic input value
    end do
  end do
  if (tsconv <= 0.0) tsconv = 0.005   ! constrain the convergence criterion to be > 0.0 and <= 0.1
  if (tsconv >  0.1) tsconv = 0.1

  read (NUNIT,'(//8x,i8)') tempn
  allocate (tempcrit(nwb,tempn),volmc(nwb,tempn))
  read (NUNIT,'(/)')
  do j=1,tempn
    read (NUNIT,'(8x,10f8.0)') (tempcrit(jw,j), jw=1,nwb)   ! Note max of 10 waterbodies
  end do
  close (NUNIT)

  do jw=1,nwb
    ifile=ifile+1
    write (segnum,'(i0)') jw
    segnum = adjustl(segnum)
    l      = len_trim(segnum)
    IF(RESTART_IN)THEN
      OPEN  (ifile,file='Volume_wb'//segnum(1:l)//'.opt',POSITION='APPEND')
      JDAY1=0.0
      REWIND (IFILE)
      READ   (IFILE,'(/)',END=15)
      DO WHILE (JDAY1 < JDAY)
        READ (IFILE,'(F10.0)',END=15) JDAY1
      END DO
      BACKSPACE (IFILE)
15    JDAY1=0.0
    ELSE
      open (ifile,file='Volume_wb'//segnum(1:l)//'.opt',status='unknown')
      write(ifile,4315)
    END IF
  end do

4315  format("jday    Volume    ",<tempn>("Volcrit      "))

  if (tempc == '      ON') then
    do j=1,numtempc
      if (tcyearly(j) == '     OFF') then
        daytest=jday
      else
        daytest=real(jdayg)+jday-int(jday)
      end if
      if (daytest >= tctsrt(j) .and. daytest < tctend(j)) then

      ! initializing structure elevation
        if (tcntr(j) == '      ST') then
          jb = tcjb(j)                                         ! set branch index
          js = tcjs(j)                                         ! set structure index
          do nj=1,tcnelev(j)                                   ! making sure that structure is below water surface
            if (tcelev(j,nj) < elws(ds(jb))) then
              ncountc(js,jb)=nj
            ! ESTR(js,jb)=tcelev(j,ncountc(js,jb))    ! don't alter the elevation at this point.  Set it later.
              exit
            end if
          end do

      ! initializing withdrawal elevation
        else if (tcntr(j) == '      WD') then
          jwd = tcjs(j)                                        ! set withdrawal index
          do nj=1,tcnelev(j)                                   ! making sure that structure is below water surface
            if (tcelev(j,nj) < elws(iwd(jwd))) then
              ncountcw(jwd)=nj
            ! EWD(jwd)=tcelev(j,ncountcw(jwd))        ! don't alter the elevation at this point.  Set it later.
              exit
            end if
          end do
        end if
      end if

    ! Open dynamic selective withdrawal files
      if (DYNSEL(J) == '      ON') then
        WRITE (SEGNUM,'(I0)') J
        SEGNUM  = ADJUSTL(SEGNUM)
        L       = LEN_TRIM(SEGNUM)
        SELD(J) = 1009+J
        OPEN (SELD(J),FILE='dynselective'//SEGNUM(1:L)//'.npt',STATUS='OLD')
        
          READ(SELD(J),'(A1)')INFORMAT    ! SW 8/28/2019
        IF(INFORMAT=='$')DYNSF(J)=.TRUE.
    
        IF(DYNSF(J))THEN
        READ (SELD(J),'(/)')
        READ (SELD(J),*) NXSEL(J),TEMP2(J)
        tctemp(J)=TEMP2(J)
        READ (SELD(J),*) NXSEL(J),TEMP2(J)       
        
        ELSE
        READ (SELD(J),'(//1000F8.0)') NXSEL(J),TEMP2(J)
        tctemp(J)=TEMP2(J)
        READ (SELD(J),'(1000F8.0)') NXSEL(J),TEMP2(J)
        ENDIF       
        
        !READ (SELD(J),'(///1000F8.0)') NXSEL(J),TEMP2(J)
        !tctemp(J)=TEMP2(J)
        !READ (SELD(J),'(1000F8.0)') NXSEL(J),TEMP2(J)
      end if
    end do
  end if

! Open dynamic temperature target files for blending outlets
  if (tspltc == '      ON') then
    do j=1,numtsplt
      if (tsdynsel(j) == '      ON') then
        write (segnum,'(i0)') j
        segnum    = adjustl(segnum)
        L         = len_trim(segnum)
        tsseld(j) = 1009+numtempc+j
        open (tsseld(j),file='dynsplit_selective'//segnum(1:L)//'.npt',status='old')
        read (tsseld(j),'(///2F8.0)') nxtssel(j), tstemp2(j)
        tspltt(j) = tstemp2(j)
        read (tsseld(j),'(2F8.0)') nxtssel(j), tstemp2(j)
      end if
    end do
  end if

! Test to see if the user specified inconsistent inputs. If so, stop with an error message.
  if (tspltc == '      ON') then
    do j=1,numtsplt
      do n=1,nouts(j)-1
        do nj=n+1,nouts(j)
          if (jstsplt(j,n) == jstsplt(j,nj)) then
            write (w2err, '(A,I0)') 'w2_selective.npt USGS ERROR-- Duplicate split outlet numbers in group ', j
            ERROR_OPEN = .TRUE.     ! will trigger the program to end when this subroutine is completed
          end if
        end do
      end do
    end do
    do j=1,numtsplt-1
      do jj=j+1,numtsplt
        if ((tstsrt(jj) >= tstsrt(j) .and. tstsrt(jj) <  tstend(j)) .or.                                                           &
            (tstend(jj) >  tstsrt(j) .and. tstend(jj) <= tstend(j))) then
          if (tspltcntr(j) .eq. tspltcntr(jj) .and. (tspltcntr(j) .eq. '      WD' .or. tspltjb(j) == tspltjb(jj))) then
            do n=1,nouts(j)
              do nj=1,nouts(jj)
                if (jstsplt(j,n) == jstsplt(jj,nj)) then
                  write (w2err, '(A,I0,A)') 'w2_selective.npt USGS ERROR-- Split outlet number ', jstsplt(j,n), ' used in more than one group at a time.'
                  ERROR_OPEN = .TRUE.       ! will trigger the program to end when this subroutine is completed
                end if
              end do
            end do
          end if
        end if
      end do
    end do
    do j=1,numtsplt
      do n=1,nouts(j)
        if (tsprior(j,n) < -1) then
          write (w2err, '(A,I0,A,I0,A)') 'w2_selective.npt USGS ERROR-- Priority input for outlet ', jstsplt(j,n), ' in group ', j, ' is less than -1.'
          ERROR_OPEN = .TRUE.            ! will trigger the program to end when this subroutine is completed
        end if
        if (tsminhead(j,n) > 0.0 .and. tsmaxhead(j,n) > 0.0 .and. tsminhead(j,n) > tsmaxhead(j,n)) then
          write (wrn, '(A,I0,A,I0,A)') 'w2_selective.npt USGS WARNING-- Minimum and maximum head constraints for outlet ', jstsplt(j,n), ' in group ',   &
                                        j, ' are such that the outlet cannot ever be used.'
          WARNING_OPEN = .TRUE.
        end if
        if (tsdepth(j,n) > 0.0 .and. tsminhead(j,n) > 0.0 .and. tsdepth(j,n) < tsminhead(j,n)) then
          write (wrn, '(A,I0,A,I0,A)') 'w2_selective.npt USGS WARNING-- Depth of floating outlet ', jstsplt(j,n), ' in group ', j,                       &
              ' is shallower than the minimum head constraint.  To honor the head constraint, no flow is possible for that outlet.'
          WARNING_OPEN = .TRUE.
        end if
      end do
    end do

    if (tempc == '      ON') then
      do j=1,numtsplt
        do jj=1,numtempc
          if ((tctsrt(jj) >= tstsrt(j) .and. tctsrt(jj) <  tstend(j)) .or.                                                         &
              (tctend(jj) >  tstsrt(j) .and. tctend(jj) <= tstend(j))) then
            if (tspltcntr(j) .eq. tcntr(jj) .and. (tspltcntr(j) .eq. '      WD' .or. tspltjb(j) == tcjb(jj))) then
              do n=1,nouts(j)
                if (jstsplt(j,n) == tcjs(jj)) then
                  write (w2err, '(A,I0,A)') 'w2_selective.npt USGS ERROR-- Outlet number ',tcjs(jj),' used in tower and blending group at same time.'
                  ERROR_OPEN = .TRUE.       ! will trigger the program to end when this subroutine is completed
                end if
              end do
            end if
          end if
        end do
      end do
    end if
  end if

  if (tempc == '      ON') then
    do j=1,numtempc-1
      do jj=j+1,numtempc
        if ((tctsrt(jj) >= tctsrt(j) .and. tctsrt(jj) <  tctend(j)) .or.                                                           &
            (tctend(jj) >  tctsrt(j) .and. tctend(jj) <= tctend(j))) then
          if (tcntr(j) .eq. tcntr(jj) .and. (tcntr(j) .eq. '      WD' .or. tcjb(j) == tcjb(jj))) then
            if (tcjs(j) == tcjs(jj)) then
              write (w2err, '(A,I0,A)') 'w2_selective.npt USGS ERROR-- Tower outlet number ', tcjs(j), ' used more than once for overlapping dates.'
              ERROR_OPEN = .TRUE.       ! will trigger the program to end when this subroutine is completed
            end if
          end if
        end if
      end do
    end do
  end if

return

End subroutine SelectiveInitUSGS


!***********************************************************************************************************************************
!**                                                   S E L E C T I V E                                                           **
!***********************************************************************************************************************************

Subroutine SelectiveUSGS
  Use Selective1USGS; USE MAIN
  USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  IMPLICIT NONE

  integer :: jj, jst, n, nj, num_noflow, ng0, prior1, prior2, ng1max, ng1min, num_left
  integer :: j2hi, j2lo, j2max, j2min, j2pref
  real    :: qall, elr, wsel, q_notblended, sum_minfrac0, sum_maxfrac1, sum_maxfrac2, sumfrac
  real    :: maxelev, minelev, blendfrac, excess_frac, addfrac, maxtemp, mintemp
  real    :: lastfrac, lastfrac2, ttarg, sumtemp, etemp, etemp1, etemp2, sumelev, elev1, elev2
  
  INTEGER JJW, KK, KS, IFILE, KSTR
  REAL DAYTEST, TCOMP, TEMPEST, TMOD

! qstr = qstrsav   ! xxx not sure how to do this yet -- need to reset QSTR when control periods expire
! qwd  = qwdsav    ! xxx not sure how to do this yet -- need to reset QWD when control periods expire

  str_active = .FALSE.
  wd_active  = .FALSE.
  j2pref     = 1

! Update some date variables and determine which outlets are being actively blended or adjusted.
  if (tspltc == '      ON') then
    do j=1,numtsplt
      if (tsyearly(j) == '     OFF') then
        daytest = jday
      else
        daytest = real(jdayg) + jday - int(jday)
      end if
      if (nxtsplit > tstsrt(j) .and. daytest <= tstsrt(j)) then
        nxtsplit = tstsrt(j)
      end if
      if (daytest >= tstsrt(j) .and. daytest < tstend(j)) then
        do jj=1,nouts(j)
          if (tspltcntr(j) == '      ST') then
            str_active(jstsplt(j,jj),tspltjb(j)) = .TRUE.
          else if (tspltcntr(j) == '      WD') then
            wd_active(jstsplt(j,jj)) = .TRUE.
          end if
        end do
      end if
    end do
  end if
  if (tempc == '      ON') then
    do j=1,numtempc
      if (tcyearly(j) == '     OFF') then
        daytest = jday
      else
        daytest = real(jdayg) + jday - int(jday)
      end if
      if (daytest >= tctsrt(j) .and. daytest < tctend(j)) then
        if (tcntr(j) == '      ST') then
          str_active(tcjs(j),tcjb(j)) = .TRUE.
        else if (tcntr(j) == '      WD') then
          wd_active(tcjs(j)) = .TRUE.
        end if
      end if
    end do
  end if

! Reset elevations of outlets back to original values outside of control periods
  do jst=1,nst
    do jb=1,nbr
      if (.not. str_active(jst,jb)) estr(jst,jb) = estrsav(jst,jb)
    end do
  end do
  do jwd=1,nwd
    if (.not. wd_active(jwd)) ewd(jwd) = ewdsav(jwd)
  end do


! Check to see if it's time to update temperature targets and flow fractions for blended groups.
  if (tspltc=='      ON' .and. jday .ge. nxtsplit) then

  ! Update the temperature targets
    do j=1,numtsplt
      if (tsdynsel(j) == '      ON') then
        do while (jday >= nxtssel(j))
          tspltt(j) = tstemp2(j)
          read (tsseld(j),'(2F8.0)') nxtssel(j), tstemp2(j)
        end do
      end if
    end do

    do j=1,numtsplt
      qall    = 0.0                                                          ! sum up all the flows
      sumfrac = 0.0                                                          ! sum of flow fraction multipliers
      do jj=1,nouts(j)
        if (tspltcntr(j) == '      ST') then
          qall    = qall    + qstr(jstsplt(j,jj),tspltjb(j))
          sumfrac = sumfrac + qstrfrac(jstsplt(j,jj),tspltjb(j))
        else if (tspltcntr(j) == '      WD') then
          qall    = qall    + qwd(jstsplt(j,jj))
          sumfrac = sumfrac + qwdfrac(jstsplt(j,jj))
        end if
      end do
      if (tsyearly(j) == '     OFF') then
        daytest = jday
      else
        daytest = real(jdayg) + jday - int(jday)
      end if

    ! Do blending calculations if date is in correct window or the window was just entered.
    ! If flows are zero and flow fractions are already initialized, then leave everything alone.
      if (daytest >= tstsrt(j) .and. daytest < tstend(j) .and.                                                                     &
          (qall > 0.0 .or. sumfrac < 0.001 .or. daytest < tstsrt(j) + tspltfreq)) then

      ! Do structures first
        if (tspltcntr(j) == '      ST') then
          jb = tspltjb(j)                                                    ! set branch index
          do jw=1,nwb                                                        ! set waterbody index
            if (jb >= bs(jw) .and. jb <= be(jw)) exit
          end do
          no_flow(j,:) = .FALSE.
          num_noflow   = 0
          do jj=1,nouts(j)
            jst          = jstsplt(j,jj)
            elr          = sina(jb) * dlx(ds(jb)) * 0.5
            wsel         = elws(ds(jb)) - elr                                ! compute water-surface elevation
            estr(jst,jb) = estrsav(jst,jb)                                   ! reset outlet elevation to original
            if (tstype(j,jj) .eq. "FLOAT") then
              estr(jst,jb) = wsel - tsdepth(j,jj)
            else if (estr(jst,jb) > wsel) then
              if (elcontspl(j) == '     OFF') then
                no_flow(j,jj) = .TRUE.                                       ! no flow-- high and dry
                num_noflow    = num_noflow + 1
              else
                estr(jst,jb) = wsel                                          ! poor man's floating outlet
              end if
            end if
            if (.not. no_flow(j,jj) .and. tsminhead(j,jj) > 0.0 .and. wsel - estr(jst,jb) < tsminhead(j,jj)) then
              no_flow(j,jj) = .TRUE.                                         ! minimum head criterion not met -- no flow
              num_noflow    = num_noflow + 1
            end if
            if (.not. no_flow(j,jj) .and. tsmaxhead(j,jj) > 0.0 .and. wsel - estr(jst,jb) > tsmaxhead(j,jj)) then
              no_flow(j,jj) = .TRUE.                                         ! maximum head criterion exceeded -- no flow
              num_noflow    = num_noflow + 1
            end if
            do k=ktwb(jw),kb(ds(jb))
              if (el(k,ds(jb))-elr < estr(jst,jb)) exit
            end do
            kstrsplt(jj)     = min(k-1,kb(ds(jb)))
            qstrfrac(jst,jb) = 0.0                                           ! initialize flow fractions
          end do

        ! Use priority inputs to determine which outlets to use
          prior1 = -999
          do jj=1,nouts(j)
            if (.not. no_flow(j,jj) .and. tsprior(j,jj) >= 0) then
              if (prior1 == -999 .or. tsprior(j,jj) < prior1) prior1 = tsprior(j,jj)
            end if
          end do
          prior2 = -999
          do jj=1,nouts(j)
            if (.not. no_flow(j,jj) .and. tsprior(j,jj) >= 0 .and. tsprior(j,jj) > prior1) then
              if (prior2 == -999 .or. tsprior(j,jj) < prior2) prior2 = tsprior(j,jj)
            end if
          end do

        ! Outlets with a priority of -1 get used, but are not blended
          ng0 = 0
          q_notblended = 0.0
          do jj=1,nouts(j)
            jst = jstsplt(j,jj)
            if (.not. no_flow(j,jj) .and. tsprior(j,jj) == -1) then
              ng0 = ng0 + 1
              nout0(ng0) = jj
              if (qstr(jst,jb) > tsmaxflow(j,jj) .and. tsmaxflow(j,jj) > 0.0) then
                q_notblended = q_notblended + tsmaxflow(j,jj)
                qstrfrac(jst,jb) = tsmaxflow(j,jj) / qall
              else if (qall > 0.0) then
                q_notblended = q_notblended + qstr(jst,jb)
                qstrfrac(jst,jb) = qstr(jst,jb) / qall
              end if
            end if
          end do
          sum_minfrac0 = 0.0
          if (qall > 0.0) sum_minfrac0 = q_notblended / qall

        ! Outlets with priority 1 and 2 may be used and blended.
          ng1 = 0
          ng2 = 0
          sum_minfrac1 = 0.0
          sum_minfrac2 = 0.0
          sum_maxfrac1 = 0.0
          sum_maxfrac2 = 0.0
          do jj=1,nouts(j)
            if (.not. no_flow(j,jj)) then
              if (tsprior(j,jj) == prior1) then
                ng1 = ng1 + 1
                nout1(ng1) = jj
                maxfrac1(ng1) = 1.0
                if (qall > 0.0 .and. tsmaxflow(j,jj) > 0.0) maxfrac1(ng1) = min(1.0, tsmaxflow(j,jj)/qall)
                minfrac1(ng1) = tsminfrac(j,jj)
                if (tsminfrac(j,jj) < 0.0) then
                  minfrac1(ng1) = 0.0
                  if (qall > 0.0) minfrac1(ng1) = min(1.0, abs(tsminfrac(j,jj))/qall)
                end if
                if (minfrac1(ng1) > maxfrac1(ng1)) minfrac1(ng1) = maxfrac1(ng1)
                sum_minfrac1 = sum_minfrac1 + minfrac1(ng1)
                sum_maxfrac1 = sum_maxfrac1 + maxfrac1(ng1)

              else if (tsprior(j,jj) == prior2) then
                ng2 = ng2 + 1
                nout2(ng2) = jj
                maxfrac2(ng2) = 1.0
                if (qall > 0.0 .and. tsmaxflow(j,jj) > 0.0) maxfrac2(ng2) = min(1.0, tsmaxflow(j,jj)/qall)
                minfrac2(ng2) = tsminfrac(j,jj)
                if (tsminfrac(j,jj) < 0.0) then
                  minfrac2(ng2) = 0.0
                  if (qall > 0.0) minfrac2(ng2) = min(1.0, abs(tsminfrac(j,jj))/qall)
                end if
                if (minfrac2(ng2) > maxfrac2(ng2)) minfrac2(ng2) = maxfrac2(ng2)
                sum_minfrac2 = sum_minfrac2 + minfrac2(ng2)
                sum_maxfrac2 = sum_maxfrac2 + maxfrac2(ng2)
              end if
            end if
          end do

        ! If minimum flows are overspecified, then the priority 2 minimum flow fractions are decreased.
          if (ng2 > 0 .and. sum_minfrac0 + sum_minfrac1 + sum_minfrac2 > 1.0) then
            if (sum_minfrac0 + sum_minfrac1 >= 1.0) then
              ng2 = 0
              sum_minfrac2 = 0.0
            else
              do n=1,ng2
                minfrac2(n) = minfrac2(n) * (1.0 - sum_minfrac0 - sum_minfrac1) / sum_minfrac2
              end do
              sum_minfrac2 = 1.0 - sum_minfrac0 - sum_minfrac1
            end if
          end if

        ! If minimum flows are still overspecified, then the priority 1 minimum flow fractions are decreased.
          if (ng1 > 0 .and. sum_minfrac0 + sum_minfrac1 > 1.0) then
            if (sum_minfrac0 >= 1.0) then
              ng1 = 0
              sum_minfrac1 = 0.0
            else
              do n=1,ng1
                minfrac1(n) = minfrac1(n) * (1.0 - sum_minfrac0) / sum_minfrac1
              end do
              sum_minfrac1 = 1.0 - sum_minfrac0
            end if
          end if

        ! If group 1 has 3 or more outlets and group 2 has no outlets, then redistribute priorities based on elevation.
        ! Keep the highest and lowest elevation outlets in group 1, and put other active outlets into nonblended category
        ! with their minimum flows.  If ties in elevation exist, go with the first in the list.
          if (ng1 > 2 .and. ng2 == 0) then
            ng1max  = 1
            ng1min  = 1
            jst     = jstsplt(j,nout1(1))
            maxelev = estr(jst,jb)
            minelev = estr(jst,jb)
            do n=2,ng1
              jst = jstsplt(j,nout1(n))
              if (estr(jst,jb) > maxelev) then
                maxelev = estr(jst,jb)
                ng1max  = n
              else if (estr(jst,jb) < minelev) then
                minelev = estr(jst,jb)
                ng1min  = n
              end if
            end do
            blendfrac = 1.0 - sum_minfrac0 - sum_minfrac1 + minfrac1(ng1max) + minfrac1(ng1min)
            if (maxfrac1(ng1max) + maxfrac1(ng1min) < blendfrac) then
              if (sum_maxfrac1 < 1.0 - sum_minfrac0) then
                write (wrn,'(A,I0,A,F0.3)') 'Warning-- Maximum flows for outlets exceeded for group ', j, ' at day ', jday
                WARNING_OPEN = .TRUE.
                do n=1,ng1
                  if (n .ne. ng1max .and. n .ne. ng1min) minfrac1(n) = maxfrac1(n)
                end do
              else
                excess_frac = blendfrac - maxfrac1(ng1max) - maxfrac1(ng1min)
                num_left = ng1 - 2
                do nj=1,ng1                                        ! iterative process to redistribute excess flows
                  if (num_left > 0 .and. excess_frac > 0.0) then
                    addfrac = excess_frac / num_left
                    do n=1,ng1
                      if (n .ne. ng1max .and. n .ne. ng1min .and. maxfrac1(n) - minfrac1(n) > 0.00001) then
                        if (minfrac1(n) + addfrac > maxfrac1(n)) then
                          num_left    = num_left - 1
                          excess_frac = excess_frac - (maxfrac1(n) - minfrac1(n))
                          minfrac1(n) = maxfrac1(n)
                        else
                          excess_frac = excess_frac - addfrac
                          minfrac1(n) = minfrac1(n) + addfrac
                        end if
                      end if
                    end do
                  end if
                end do
              end if
            end if
            do n=1,ng1                                       ! assign the other priority 1 outlets to nonblended status
              if (n .ne. ng1max .and. n .ne. ng1min) then
                ng0              = ng0 + 1
                nout0(ng0)       = nout1(n)
                jst              = jstsplt(j,nout1(n))
                sum_minfrac0     = sum_minfrac0 + minfrac1(n)
                q_notblended     = q_notblended + qall * minfrac1(n)
                qstrfrac(jst,jb) = minfrac1(n)
              end if
            end do
            ng1 = 1                                ! rearrange outlets-- one in each priority group, but same priority
            ng2 = 1
            nout2(1)     = nout1(ng1min)
            minfrac2(1)  = minfrac1(ng1min)
            maxfrac2(1)  = maxfrac1(ng1min)
            sum_minfrac2 = minfrac1(ng1min)
            sum_maxfrac2 = maxfrac1(ng1min)
            nout1(1)     = nout1(ng1max)
            minfrac1(1)  = minfrac1(ng1max)
            maxfrac1(1)  = maxfrac1(ng1max)
            sum_minfrac1 = minfrac1(ng1max)
            sum_maxfrac1 = maxfrac1(ng1max)
            prior2       = prior1
          end if

        ! If only two blended outlets, ensure that they are in separate groups.
          if (ng1 == 2 .and. ng2 == 0) then
            ng1 = 1
            ng2 = 1
            nout2(1)     = nout1(2)
            minfrac2(1)  = minfrac1(2)
            maxfrac2(1)  = maxfrac1(2)
            sum_minfrac2 = minfrac1(2)
            sum_maxfrac2 = maxfrac1(2)
            sum_minfrac1 = minfrac1(1)
            sum_maxfrac1 = maxfrac1(1)
            prior2       = prior1
          end if


        ! Begin the blending decisions.
        ! No usable outlets.  All flow fractions remain at zero.
          if (nouts(j) == num_noflow) then
            write (wrn,'(A,I0,A,F0.3)') 'Warning-- All outlets dry or unusable for group ', j, ' at day ', jday
            WARNING_OPEN = .TRUE.

        ! Only nonblended outlets.
          else if (nouts(j) == ng0) then
            write (wrn,'(A,I0,A,F0.3)') 'Warning-- Only nonblended outlets present in group ', j, ' at day ', jday
            WARNING_OPEN = .TRUE.

        ! Only one blended outlet.  It gets all of the blended flow, but must not exceed its maximum flow criterion.
          else if (ng1 + ng2 == 1) then
            jst = jstsplt(j,nout1(1))
            qstrfrac(jst,jb) = 1.0 - sum_minfrac0
            if (qall - q_notblended > tsmaxflow(j,nout1(1)) .and. tsmaxflow(j,nout1(1)) > 0.0) then
              qstrfrac(jst,jb) = tsmaxflow(j,nout1(1)) / qall
              write (wrn,'(A,A,I0,A,I0,A,F0.3)') 'Warning-- Total release flow rate decreased to comply with maximum flow ',       &
                                                 'criterion for structure ', jst, ' in group ', j, ' at day ', jday
              WARNING_OPEN = .TRUE.
            end if

        ! Minimum flows comprise entire release.  No blending calculations required.
          else if (abs(1.0 - sum_minfrac0 - sum_minfrac1 - sum_minfrac2) <= 0.000001) then
            do n=1,ng1
              jst = jstsplt(j,nout1(n))
              qstrfrac(jst,jb) = minfrac1(n)
            end do
            do n=1,ng2
              jst = jstsplt(j,nout2(n))
              qstrfrac(jst,jb) = minfrac2(n)
            end do

        ! More than one usable outlet, and blending among priority 1 outlet(s) and priority 2 outlet(s) required.
          else
            id = ds(jb)                                                      ! needed for downstream_withdrawal_estimate
            kt = ktwb(jw)                                                    ! needed for downstream_withdrawal_estimate

          ! Warn the user if maximum flow criteria are likely to decrease the specified outflows.
            if (sum_minfrac0 + sum_maxfrac1 + sum_maxfrac2 < 1.0) then
              write (wrn,'(A,A,I0,A,F0.3)') 'Warning-- Total release flow rate may be decreased to comply with maximum flow ',     &
                                            'criteria for structures in group ', j, ' at day ', jday
              WARNING_OPEN = .TRUE.
            end if

          ! Set the initial release fractions.  Ensure that maximum flows are not exceeded.
            qfrac1 = sum_minfrac1 + 0.5 * (1.0 - sum_minfrac0 - sum_minfrac1 - sum_minfrac2)
            qfrac1 = min(sum_maxfrac1, qfrac1)
            qfrac2 = 1.0 - sum_minfrac0 - qfrac1
            if (qfrac2 > sum_maxfrac2) then
              excess_frac = qfrac2 - sum_maxfrac2
              if (qfrac1 + excess_frac <= sum_maxfrac1) then
                qfrac1 = qfrac1 + excess_frac
              else
                qfrac1 = sum_maxfrac1
              end if
              qfrac2 = sum_maxfrac2
            end if
            j2pref = 1
            call Set_Flow_Fracs(j, jb, j2pref)                               ! set flow fractions; redistribute if maxfrac exceeded

          ! If priority 2 outlets are not sharing flows, identify the two with the highest and lowest elevations
            if (.not. share_flow(j) .and. ng2 > 1) then
              j2hi    = 1
              j2lo    = 1
              jst     = jstsplt(j,nout2(1))
              maxelev = estr(jst,jb)
              minelev = estr(jst,jb)
              do n=2,ng2
                jst = jstsplt(j,nout2(n))
                if (estr(jst,jb) > maxelev) then
                  maxelev = estr(jst,jb)
                  j2hi    = n
                else if (estr(jst,jb) < minelev) then
                  minelev = estr(jst,jb)
                  j2lo    = n
                end if
              end do
            end if

          ! Get weighted blend of all nonblended release temperatures
            ttarg = tspltt(j)
            if (sum_minfrac0 > 0.0) then
              sumtemp = 0.0
              do n=1,ng0
                jst = jstsplt(j,nout0(n))
                qstr(jst,jb) = qall * qstrfrac(jst,jb)
                if (qstr(jst,jb) > 0.0) then
                  call downstream_withdrawal_estimate(jst,etemp,estr(jst,jb))
                else
                  etemp = t2(kstrsplt(nout0(n)),ds(jb))                      ! Use temperature at outlet elevation if no flow
                end if
                sumtemp = sumtemp + qstrfrac(jst,jb) * etemp
              end do
              etemp = sumtemp / sum_minfrac0
              ttarg = (ttarg - sum_minfrac0 * etemp) / (1.0 - sum_minfrac0)  ! New temperature target for blended releases
            end if

          ! Need an iterative approach because released T depends on Q
            lastfrac = qfrac1
            do jj=1,8                                                        ! Maximum of eight iterations
              lastfrac2 = lastfrac
              lastfrac  = qfrac1

              sumtemp = 0.0
              sumelev = 0.0
              do n=1,ng1                                                     ! Get weighted temp and elevation for group 1
                jst = jstsplt(j,nout1(n))
                qstr(jst,jb) = qall * qstrfrac(jst,jb)
                if (qstr(jst,jb) > 0.0) then
                  call downstream_withdrawal_estimate(jst,etemp,estr(jst,jb))
                else
                  etemp = t2(kstrsplt(nout1(n)),ds(jb))                      ! Use temperature at outlet elevation if no flow
                end if
                if (qfrac1 > 0.0) then
                  sumtemp = sumtemp + qstrfrac(jst,jb) * etemp
                  sumelev = sumelev + qstrfrac(jst,jb) * estr(jst,jb)
                else
                  sumtemp = sumtemp + etemp
                  sumelev = sumelev + estr(jst,jb)
                end if
              end do
              if (qfrac1 > 0.0) then
                etemp1 = sumtemp / qfrac1                                    ! Weighted temperature from group 1 outlets
                elev1  = sumelev / qfrac1                                    ! Weighted elevation of group 1 outlets
              else
                etemp1 = sumtemp / ng1
                elev1  = sumelev / ng1
              end if

              if (share_flow(j) .or. ng2 < 2) then                           ! Get weighted temp and elevation for group 2
                sumtemp = 0.0                                                ! ...when flows are shared among outlets
                sumelev = 0.0
                do n=1,ng2
                  jst = jstsplt(j,nout2(n))
                  qstr(jst,jb) = qall * qstrfrac(jst,jb)
                  if (qstr(jst,jb) > 0.0) then
                    call downstream_withdrawal_estimate(jst,etemp,estr(jst,jb))
                  else
                    etemp = t2(kstrsplt(nout2(n)),ds(jb))                    ! Use temperature at outlet elevation if no flow
                  end if
                  if (qfrac2 > 0.0) then
                    sumtemp = sumtemp + qstrfrac(jst,jb) * etemp
                    sumelev = sumelev + qstrfrac(jst,jb) * estr(jst,jb)
                  else
                    sumtemp = sumtemp + etemp
                    sumelev = sumelev + estr(jst,jb)
                  end if
                end do
                if (qfrac2 > 0.0) then
                  etemp2 = sumtemp / qfrac2                                  ! Weighted temperature from group 2 outlets
                  elev2  = sumelev / qfrac2                                  ! Weighted elevation of group 2 outlets
                else
                  etemp2 = sumtemp / ng2
                  elev2  = sumelev / ng2
                end if

              else                                                           ! ...and when flows are not shared
                if (qfrac2 == 0.0) then
                  do n=1,ng2
                    jst       = jstsplt(j,nout2(n))
                    splt2t(n) = t2(kstrsplt(nout2(n)),ds(jb))
                    splt2e(n) = estr(jst,jb)
                  end do
                else
                  do nj=1,ng2                                                ! Find the temperatures produced in group 2
                    sumtemp = 0.0                                            ! by testing when each outlet is preferred
                    sumelev = 0.0
                    call Set_Flow_Fracs2(j, jb, nj)
                    do n=1,ng2
                      jst = jstsplt(j,nout2(n))
                      qstr(jst,jb) = qall * qstrfrac(jst,jb)
                      if (qstr(jst,jb) > 0.0) then
                        call downstream_withdrawal_estimate(jst,etemp,estr(jst,jb))
                      else
                        etemp = t2(kstrsplt(nout2(n)),ds(jb))
                      end if
                      sumtemp = sumtemp + qstrfrac(jst,jb) * etemp
                      sumelev = sumelev + qstrfrac(jst,jb) * estr(jst,jb)
                    end do
                    splt2t(nj) = sumtemp / qfrac2
                    splt2e(nj) = sumelev / qfrac2
                  end do
                end if
                j2max   = 1
                j2min   = 1
                maxtemp = splt2t(1)
                mintemp = splt2t(1)
                do n=2,ng2
                  if (splt2t(n) > maxtemp) then
                    maxtemp = splt2t(n)
                    j2max   = n
                  else if (splt2t(n) < mintemp) then
                    mintemp = splt2t(n)
                    j2min   = n
                  end if
                end do
                if (ttarg < etemp1 - 0.001) then                             ! need a colder temp from group 2
                  if (maxtemp - mintemp > 0.001) then
                    etemp2 = splt2t(j2min)                                   ! preferred outlet is the coldest one
                    elev2  = splt2e(j2min)
                    j2pref = j2min
                  else
                    etemp2 = splt2t(j2lo)                                    ! preferred outlet is the lowest one
                    elev2  = splt2e(j2lo)
                    j2pref = j2lo
                  end if
                else if (ttarg > etemp1 + 0.001) then                        ! need a warmer temp from group 2
                  if (maxtemp - mintemp > 0.001) then
                    etemp2 = splt2t(j2max)                                   ! preferred outlet is the warmest one
                    elev2  = splt2e(j2max)
                    j2pref = j2max
                  else
                    etemp2 = splt2t(j2hi)                                    ! preferred outlet is the highest one
                    elev2  = splt2e(j2hi)
                    j2pref = j2hi
                  end if
                else
                  etemp2 = splt2t(1)                                         ! if temp is close to target, choose first outlet
                  elev2  = splt2e(1)
                  j2pref = 1
                end if
              end if

            ! Target temperature is less than either outlet temperature.
              if (ttarg < etemp1 .and. ttarg < etemp2) then
                qfrac1 = sum_minfrac1                                        ! default for if/then cases
                if (abs(etemp1 - etemp2) < 0.001) then
                  if (prior1 == prior2) then                                 ! If each outlet has the same priority level, then...
                    if (elev1 <= elev2) then                                 ! Choose lower outlet if both have same temperature.
                      qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2
                    end if
                  else if (prior1 < prior2) then
                    qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2               ! Choose higher priority outlet if temps are the same.
                  end if
                else if (etemp1 < etemp2) then                               ! If temps are different, choose the one closer
                  qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2                 ! to target temperature
                end if

            ! Target temperature is greater than either outlet temperature.
              else if (ttarg > etemp1 .and. ttarg > etemp2) then
                qfrac1 = sum_minfrac1                                        ! default for if/then cases
                if (abs(etemp1 - etemp2) < 0.001) then
                  if (prior1 == prior2) then                                 ! If each outlet has the same priority level, then...
                    if (elev1 >= elev2) then                                 ! Choose upper outlet if both have same temperature.
                      qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2
                    end if
                  else if (prior1 < prior2) then
                    qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2               ! Choose higher priority outlet if temps are the same.
                  end if
                else if (etemp1 > etemp2) then                               ! If temps are different, choose the one closer
                  qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2                 ! to target temperature
                end if

            ! Target temperature is essentially the same as the two outlet temperatures.
              else if (abs(etemp1 - etemp2) < 0.001) then
                qfrac1 = sum_minfrac1                                        ! default for if/then cases
                if (prior1 == prior2) then                                   ! If each outlet has the same priority level, then...
                  qfrac1 = sum_minfrac1 + 0.5 * (1.0 - sum_minfrac0 - sum_minfrac1 - sum_minfrac2)    ! Split the flow equally.
                else if (prior1 < prior2) then
                  qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2                 ! Choose higher priority outlet if temps are the same.
                end if

            ! Target temperature is between the two outlet temperatures.
              else
                qfrac1 = (1.0 - sum_minfrac0) * abs((ttarg-etemp2)/(etemp1-etemp2+NONZERO))
                qfrac1 = max(sum_minfrac1, qfrac1)
                qfrac1 = min(1.0 - sum_minfrac0 - sum_minfrac2, qfrac1)
              end if
              qfrac1 = min(sum_maxfrac1, qfrac1)
              qfrac2 = 1.0 - sum_minfrac0 - qfrac1
              if (qfrac2 > sum_maxfrac2) then
                excess_frac = qfrac2 - sum_maxfrac2
                if (qfrac1 + excess_frac <= sum_maxfrac1) then
                  qfrac1 = qfrac1 + excess_frac
                else
                  qfrac1 = sum_maxfrac1
                end if
                qfrac2 = sum_maxfrac2
              end if

            ! Set flow fractions for individual outlets and redistribute flows if maximum flow fractions exceeded.
              call Set_Flow_Fracs(j, jb, j2pref)

            ! Exit the loop if the latest flow fraction calculation agrees with the previous one.
            ! Exit if no flow, because no iteration requried in that case.
              if (abs(lastfrac - qfrac1) < tsconv .or. qall == 0.0) exit
            end do

          ! Check to see if iterative solution did not converge.
            if (abs(lastfrac - qfrac1) >= tsconv .and. qall > 0.0) then
              write (wrn,'(A,F0.3,3(A,F0.4))') 'Flow fraction calculations not converging at day ', jday,                          &
                                               '  Current: ', qfrac1, ' Last: ', lastfrac, ' Next-to-last: ', lastfrac2
              WARNING_OPEN = .TRUE.

            ! Check to see if the iterative solution is unstable.  If so, use priorities to assign releases.
            ! Criteria:  change is at least 0.1 and most recent change is in opposite direction from previous change.
              if (abs(lastfrac - qfrac1) >= 0.1 .and. (qfrac1-lastfrac)*(lastfrac-lastfrac2) < 0.0) then
                qfrac1 = sum_minfrac1                                        ! default for if/then cases
                if (prior1 < prior2) then                                    ! group 1 is higher priority
                  qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2
                else                                                         ! else, fulfill minima and split the rest
                  qfrac1 = sum_minfrac1 + 0.5 * (1.0 - sum_minfrac0 - sum_minfrac1 - sum_minfrac2)
                end if
                qfrac1 = min(sum_maxfrac1, qfrac1)
                qfrac2 = 1.0 - sum_minfrac0 - qfrac1
                if (qfrac2 > sum_maxfrac2) then
                  excess_frac = qfrac2 - sum_maxfrac2
                  if (qfrac1 + excess_frac <= sum_maxfrac1) then
                    qfrac1 = qfrac1 + excess_frac
                  else
                    qfrac1 = sum_maxfrac1
                  end if
                  qfrac2 = sum_maxfrac2
                end if
                call Set_Flow_Fracs(j, jb, j2pref)                           ! set flow fractions; redistribute if maxfrac exceeded
              end if
            end if
          end if

        ! Set final flows to go with the flow fractions.  May not be necessary, but do it anyway.
          do jj=1,nouts(j)
            qstr(jstsplt(j,jj),jb) = qall * qstrfrac(jstsplt(j,jj),jb)
          end do


      ! Do Withdrawals next
        else if (tspltcntr(j) == '      WD') then
          jwd = jstsplt(j,1)                                                 ! assume withdrawals are from same branch and waterbody
          do jb=1,nbr
            if (iwd(jwd) >= us(jb) .and. iwd(jwd) <= ds(jb)) exit
          end do
          do jw=1,nwb
            if (jb >= bs(jw) .and. jb <= be(jw)) exit
          end do
          no_flow(j,:) = .FALSE.
          num_noflow   = 0
          do jj=1,nouts(j)
            jwd      = jstsplt(j,jj)
            elr      = sina(jb) * dlx(iwd(jwd)) * 0.5
            wsel     = elws(iwd(jwd)) - elr                                  ! compute water-surface elevation
            ewd(jwd) = ewdsav(jwd)                                           ! reset outlet elevation to original
            if (tstype(j,jj) .eq. "FLOAT") then
              ewd(jwd) = wsel - tsdepth(j,jj)
            else if (ewd(jwd) > wsel) then
              if (elcontspl(j) == '     OFF') then
                no_flow(j,jj) = .TRUE.                                       ! no flow-- high and dry
                num_noflow    = num_noflow + 1
              else
                ewd(jwd) = wsel                                              ! poor man's floating outlet
              end if
            end if
            if (.not. no_flow(j,jj) .and. tsminhead(j,jj) > 0.0 .and. wsel - ewd(jwd) < tsminhead(j,jj)) then
              no_flow(j,jj) = .TRUE.                                         ! minimum head criterion not met -- no flow
              num_noflow    = num_noflow + 1
            end if
            if (.not. no_flow(j,jj) .and. tsmaxhead(j,jj) > 0.0 .and. wsel - ewd(jwd) > tsmaxhead(j,jj)) then
              no_flow(j,jj) = .TRUE.                                         ! maximum head criterion exceeded -- no flow
              num_noflow    = num_noflow + 1
            end if
            do k=ktwb(jw),kb(iwd(jwd))
              if (el(k,iwd(jwd))-elr < ewd(jwd)) exit
            end do
            kstrsplt(jj) = min(k-1,kb(iwd(jwd)))
            qwdfrac(jwd) = 0.0                                               ! initialize flow fractions
          end do

        ! Use priority inputs to determine which outlets to use
          prior1 = -999
          do jj=1,nouts(j)
            if (.not. no_flow(j,jj) .and. tsprior(j,jj) >= 0) then
              if (prior1 == -999 .or. tsprior(j,jj) < prior1) prior1 = tsprior(j,jj)
            end if
          end do
          prior2 = -999
          do jj=1,nouts(j)
            if (.not. no_flow(j,jj) .and. tsprior(j,jj) >= 0 .and. tsprior(j,jj) > prior1) then
              if (prior2 == -999 .or. tsprior(j,jj) < prior2) prior2 = tsprior(j,jj)
            end if
          end do

        ! Outlets with a priority of -1 get used, but are not blended
          ng0 = 0
          q_notblended = 0.0
          do jj=1,nouts(j)
            jwd = jstsplt(j,jj)
            if (.not. no_flow(j,jj) .and. tsprior(j,jj) == -1) then
              ng0 = ng0 + 1
              nout0(ng0) = jj
              if (qwd(jwd) > tsmaxflow(j,jj) .and. tsmaxflow(j,jj) > 0.0) then
                q_notblended = q_notblended + tsmaxflow(j,jj)
                qwdfrac(jwd) = tsmaxflow(j,jj) / qall
              else if (qall > 0.0) then
                q_notblended = q_notblended + qwd(jwd)
                qwdfrac(jwd) = qwd(jwd) / qall
              end if
            end if
          end do
          sum_minfrac0 = 0.0
          if (qall > 0.0) sum_minfrac0 = q_notblended / qall

        ! Outlets with priority 1 and 2 may be used and blended.
          ng1 = 0
          ng2 = 0
          sum_minfrac1 = 0.0
          sum_minfrac2 = 0.0
          sum_maxfrac1 = 0.0
          sum_maxfrac2 = 0.0
          do jj=1,nouts(j)
            if (.not. no_flow(j,jj)) then
              if (tsprior(j,jj) == prior1) then
                ng1 = ng1 + 1
                nout1(ng1) = jj
                maxfrac1(ng1) = 1.0
                if (qall > 0.0 .and. tsmaxflow(j,jj) > 0.0) maxfrac1(ng1) = min(1.0, tsmaxflow(j,jj)/qall)
                minfrac1(ng1) = tsminfrac(j,jj)
                if (tsminfrac(j,jj) < 0.0) then
                  minfrac1(ng1) = 0.0
                  if (qall > 0.0) minfrac1(ng1) = min(1.0, abs(tsminfrac(j,jj))/qall)
                end if
                if (minfrac1(ng1) > maxfrac1(ng1)) minfrac1(ng1) = maxfrac1(ng1)
                sum_minfrac1 = sum_minfrac1 + minfrac1(ng1)
                sum_maxfrac1 = sum_maxfrac1 + maxfrac1(ng1)

              else if (tsprior(j,jj) == prior2) then
                ng2 = ng2 + 1
                nout2(ng2) = jj
                maxfrac2(ng2) = 1.0
                if (qall > 0.0 .and. tsmaxflow(j,jj) > 0.0) maxfrac2(ng2) = min(1.0, tsmaxflow(j,jj)/qall)
                minfrac2(ng2) = tsminfrac(j,jj)
                if (tsminfrac(j,jj) < 0.0) then
                  minfrac2(ng2) = 0.0
                  if (qall > 0.0) minfrac2(ng2) = min(1.0, abs(tsminfrac(j,jj))/qall)
                end if
                if (minfrac2(ng2) > maxfrac2(ng2)) minfrac2(ng2) = maxfrac2(ng2)
                sum_minfrac2 = sum_minfrac2 + minfrac2(ng2)
                sum_maxfrac2 = sum_maxfrac2 + maxfrac2(ng2)
              end if
            end if
          end do

        ! If minimum flows are overspecified, then the priority 2 minimum flow fractions are decreased.
          if (ng2 > 0 .and. sum_minfrac0 + sum_minfrac1 + sum_minfrac2 > 1.0) then
            if (sum_minfrac0 + sum_minfrac1 >= 1.0) then
              ng2 = 0
              sum_minfrac2 = 0.0
            else
              do n=1,ng2
                minfrac2(n) = minfrac2(n) * (1.0 - sum_minfrac0 - sum_minfrac1) / sum_minfrac2
              end do
              sum_minfrac2 = 1.0 - sum_minfrac0 - sum_minfrac1
            end if
          end if

        ! If minimum flows are still overspecified, then the priority 1 minimum flow fractions are decreased.
          if (ng1 > 0 .and. sum_minfrac0 + sum_minfrac1 > 1.0) then
            if (sum_minfrac0 >= 1.0) then
              ng1 = 0
              sum_minfrac1 = 0.0
            else
              do n=1,ng1
                minfrac1(n) = minfrac1(n) * (1.0 - sum_minfrac0) / sum_minfrac1
              end do
              sum_minfrac1 = 1.0 - sum_minfrac0
            end if
          end if

        ! If group 1 has 3 or more outlets and group 2 has no outlets, then redistribute priorities based on elevation.
        ! Keep the highest and lowest elevation outlets in group 1, and put other active outlets into nonblended category
        ! with their minimum flows.  If ties in elevation exist, go with the first in the list.
          if (ng1 > 2 .and. ng2 == 0) then
            ng1max  = 1
            ng1min  = 1
            jwd     = jstsplt(j,nout1(1))
            maxelev = ewd(jwd)
            minelev = ewd(jwd)
            do n=2,ng1
              jwd = jstsplt(j,nout1(n))
              if (ewd(jwd) > maxelev) then
                maxelev = ewd(jwd)
                ng1max  = n
              else if (ewd(jwd) < minelev) then
                minelev = ewd(jwd)
                ng1min  = n
              end if
            end do
            blendfrac = 1.0 - sum_minfrac0 - sum_minfrac1 + minfrac1(ng1max) + minfrac1(ng1min)
            if (maxfrac1(ng1max) + maxfrac1(ng1min) < blendfrac) then
              if (sum_maxfrac1 < 1.0 - sum_minfrac0) then
                write (wrn,'(A,I0,A,F0.3)') 'Warning-- Maximum flows for outlets exceeded for group ', j, ' at day ', jday
                WARNING_OPEN = .TRUE.
                do n=1,ng1
                  if (n .ne. ng1max .and. n .ne. ng1min) minfrac1(n) = maxfrac1(n)
                end do
              else
                excess_frac = blendfrac - maxfrac1(ng1max) - maxfrac1(ng1min)
                num_left = ng1 - 2
                do nj=1,ng1                                        ! iterative process to redistribute excess flows
                  if (num_left > 0 .and. excess_frac > 0.0) then
                    addfrac = excess_frac / num_left
                    do n=1,ng1
                      if (n .ne. ng1max .and. n .ne. ng1min .and. maxfrac1(n) - minfrac1(n) > 0.00001) then
                        if (minfrac1(n) + addfrac > maxfrac1(n)) then
                          num_left    = num_left - 1
                          excess_frac = excess_frac - (maxfrac1(n) - minfrac1(n))
                          minfrac1(n) = maxfrac1(n)
                        else
                          excess_frac = excess_frac - addfrac
                          minfrac1(n) = minfrac1(n) + addfrac
                        end if
                      end if
                    end do
                  end if
                end do
              end if
            end if
            do n=1,ng1                                       ! assign the other priority 1 outlets to nonblended status
              if (n .ne. ng1max .and. n .ne. ng1min) then
                ng0          = ng0 + 1
                nout0(ng0)   = nout1(n)
                jwd          = jstsplt(j,nout1(n))
                sum_minfrac0 = sum_minfrac0 + minfrac1(n)
                q_notblended = q_notblended + qall * minfrac1(n)
                qwdfrac(jwd) = minfrac1(n)
              end if
            end do
            ng1 = 1                                ! rearrange outlets-- one in each priority group, but same priority
            ng2 = 1
            nout2(1)     = nout1(ng1min)
            minfrac2(1)  = minfrac1(ng1min)
            maxfrac2(1)  = maxfrac1(ng1min)
            sum_minfrac2 = minfrac1(ng1min)
            sum_maxfrac2 = maxfrac1(ng1min)
            nout1(1)     = nout1(ng1max)
            minfrac1(1)  = minfrac1(ng1max)
            maxfrac1(1)  = maxfrac1(ng1max)
            sum_minfrac1 = minfrac1(ng1max)
            sum_maxfrac1 = maxfrac1(ng1max)
            prior2       = prior1
          end if

        ! If only two blended outlets, ensure that they are in separate groups.
          if (ng1 == 2 .and. ng2 == 0) then
            ng1 = 1
            ng2 = 1
            nout2(1)     = nout1(2)
            minfrac2(1)  = minfrac1(2)
            maxfrac2(1)  = maxfrac1(2)
            sum_minfrac2 = minfrac1(2)
            sum_maxfrac2 = maxfrac1(2)
            sum_minfrac1 = minfrac1(1)
            sum_maxfrac1 = maxfrac1(1)
            prior2       = prior1
          end if


        ! Begin the blending decisions.
        ! No usable outlets.  All flow fractions remain at zero.
          if (nouts(j) == num_noflow) then
            write (wrn,'(A,I0,A,F0.3)') 'Warning-- All outlets dry or unusable for group ', j, ' at day ', jday
            WARNING_OPEN = .TRUE.

        ! Only nonblended outlets.
          else if (nouts(j) == ng0) then
            write (wrn,'(A,I0,A,F0.3)') 'Warning-- Only nonblended outlets present in group ', j, ' at day ', jday
            WARNING_OPEN = .TRUE.

        ! Only one blended outlet.  It gets all of the blended flow, but must not exceed its maximum flow criterion.
          else if (ng1 + ng2 == 1) then
            jwd = jstsplt(j,nout1(1))
            qwdfrac(jwd) = 1.0 - sum_minfrac0
            if (qall - q_notblended > tsmaxflow(j,nout1(1)) .and. tsmaxflow(j,nout1(1)) > 0.0) then
              qwdfrac(jwd) = tsmaxflow(j,nout1(1)) / qall
              write (wrn,'(A,A,I0,A,I0,A,F0.3)') 'Warning-- Total release flow rate decreased to comply with maximum flow ',       &
                                                 'criterion for withdrawal ', jwd, ' in group ', j, ' at day ', jday
              WARNING_OPEN = .TRUE.
            end if

        ! Minimum flows comprise entire release.  No blending calculations required.
          else if (abs(1.0 - sum_minfrac0 - sum_minfrac1 - sum_minfrac2) <= 0.000001) then
            do n=1,ng1
              jwd = jstsplt(j,nout1(n))
              qwdfrac(jwd) = minfrac1(n)
            end do
            do n=1,ng2
              jwd = jstsplt(j,nout2(n))
              qwdfrac(jwd) = minfrac2(n)
            end do

        ! More than one usable outlet, and blending among priority 1 outlet(s) and priority 2 outlet(s) required.
          else
            I  = iwd(jstsplt(j,nout1(1)))                                    ! needed for lateral_withdrawal_estimate
            kt = ktwb(jw)                                                    ! needed for lateral_withdrawal_estimate

          ! Warn the user if maximum flow criteria are likely to decrease the specified outflows.
            if (sum_minfrac0 + sum_maxfrac1 + sum_maxfrac2 < 1.0) then
              write (wrn,'(A,A,I0,A,F0.3)') 'Warning-- Total release flow rate may be decreased to comply with maximum flow ',     &
                                            'criteria for withdrawals in group ', j, ' at day ', jday
              WARNING_OPEN = .TRUE.
            end if

          ! Set the initial release fractions.  Ensure that maximum flows are not exceeded.
            qfrac1 = sum_minfrac1 + 0.5 * (1.0 - sum_minfrac0 - sum_minfrac1 - sum_minfrac2)
            qfrac1 = min(sum_maxfrac1, qfrac1)
            qfrac2 = 1.0 - sum_minfrac0 - qfrac1
            if (qfrac2 > sum_maxfrac2) then
              excess_frac = qfrac2 - sum_maxfrac2
              if (qfrac1 + excess_frac <= sum_maxfrac1) then
                qfrac1 = qfrac1 + excess_frac
              else
                qfrac1 = sum_maxfrac1
              end if
              qfrac2 = sum_maxfrac2
            end if
            j2pref = 1
            call Set_Flow_Fracs(j, jb, j2pref)                               ! set flow fractions; redistribute if maxfrac exceeded

          ! If priority 2 outlets are not sharing flows, identify the two with the highest and lowest elevations
            if (.not. share_flow(j) .and. ng2 > 1) then
              j2hi    = 1
              j2lo    = 1
              jwd     = jstsplt(j,nout2(1))
              maxelev = ewd(jwd)
              minelev = ewd(jwd)
              do n=2,ng2
                jwd = jstsplt(j,nout2(n))
                if (ewd(jwd) > maxelev) then
                  maxelev = ewd(jwd)
                  j2hi    = n
                else if (ewd(jwd) < minelev) then
                  minelev = ewd(jwd)
                  j2lo    = n
                end if
              end do
            end if

          ! Get weighted blend of all nonblended release temperatures
            ttarg = tspltt(j)
            if (sum_minfrac0 > 0.0) then
              sumtemp = 0.0
              do n=1,ng0
                jwd = jstsplt(j,nout0(n))
                qwd(jwd) = qall * qwdfrac(jwd)
                if (qwd(jwd) > 0.0) then
                  call lateral_withdrawal_estimate(jwd,etemp,ewd(jwd))       ! Get an estimate of the temperature of the outflow
                else
                  etemp = t2(kstrsplt(nout0(n)),iwd(jwd))                    ! Use temperature at outlet elevation if no flow
                end if
                sumtemp = sumtemp + qwdfrac(jwd) * etemp
              end do
              etemp = sumtemp / sum_minfrac0
              ttarg = (ttarg - sum_minfrac0 * etemp) / (1.0 - sum_minfrac0)  ! New temperature target for blended releases
            end if

          ! Need an iterative approach because released T depends on Q
            lastfrac = qfrac1
            do jj=1,8                                                        ! Maximum of eight iterations
              lastfrac2 = lastfrac
              lastfrac  = qfrac1

              sumtemp = 0.0
              sumelev = 0.0
              do n=1,ng1                                                     ! Get weighted temp and elevation for group 1
                jwd = jstsplt(j,nout1(n))
                qwd(jwd) = qall * qwdfrac(jwd)
                if (qwd(jwd) > 0.0) then
                  call lateral_withdrawal_estimate(jwd,etemp,ewd(jwd))
                else
                  etemp = t2(kstrsplt(nout1(n)),iwd(jwd))                    ! Use temperature at outlet elevation if no flow
                end if
                if (qfrac1 > 0.0) then
                  sumtemp = sumtemp + qwdfrac(jwd) * etemp
                  sumelev = sumelev + qwdfrac(jwd) * ewd(jwd)
                else
                  sumtemp = sumtemp + etemp
                  sumelev = sumelev + ewd(jwd)
                end if
              end do
              if (qfrac1 > 0.0) then
                etemp1 = sumtemp / qfrac1                                    ! Weighted temperature from group 1 outlets
                elev1  = sumelev / qfrac1                                    ! Weighted elevation of group 1 outlets
              else
                etemp1 = sumtemp / ng1
                elev1  = sumelev / ng1
              end if

              if (share_flow(j) .or. ng2 < 2) then                           ! Get weighted temp and elevation for group 2
                sumtemp = 0.0                                                ! ...when flows are shared among outlets
                sumelev = 0.0
                do n=1,ng2
                  jwd = jstsplt(j,nout2(n))
                  qwd(jwd) = qall * qwdfrac(jwd)
                  if (qwd(jwd) > 0.0) then
                    call lateral_withdrawal_estimate(jwd,etemp,ewd(jwd))
                  else
                    etemp = t2(kstrsplt(nout2(n)),iwd(jwd))                  ! Use temperature at outlet elevation if no flow
                  end if
                  if (qfrac2 > 0.0) then
                    sumtemp = sumtemp + qwdfrac(jwd) * etemp
                    sumelev = sumelev + qwdfrac(jwd) * ewd(jwd)
                  else
                    sumtemp = sumtemp + etemp
                    sumelev = sumelev + ewd(jwd)
                  end if
                end do
                if (qfrac2 > 0.0) then
                  etemp2 = sumtemp / qfrac2                                  ! Weighted temperature from group 2 outlets
                  elev2  = sumelev / qfrac2                                  ! Weighted elevation of group 2 outlets
                else
                  etemp2 = sumtemp / ng2
                  elev2  = sumelev / ng2
                end if

              else                                                           ! ...and when flows are not shared
                if (qfrac2 == 0.0) then
                  do n=1,ng2
                    jwd       = jstsplt(j,nout2(n))
                    splt2t(n) = t2(kstrsplt(nout2(n)),iwd(jwd))
                    splt2e(n) = ewd(jwd)
                  end do
                else
                  do nj=1,ng2                                                ! Find the temperatures produced in group 2
                    sumtemp = 0.0                                            ! by testing when each outlet is preferred
                    sumelev = 0.0
                    call Set_Flow_Fracs2(j, jb, nj)
                    do n=1,ng2
                      jwd = jstsplt(j,nout2(n))
                      qwd(jwd) = qall * qwdfrac(jwd)
                      if (qwd(jwd) > 0.0) then
                        call lateral_withdrawal_estimate(jwd,etemp,ewd(jwd))
                      else
                        etemp = t2(kstrsplt(nout2(n)),iwd(jwd))
                      end if
                      sumtemp = sumtemp + qwdfrac(jwd) * etemp
                      sumelev = sumelev + qwdfrac(jwd) * ewd(jwd)
                    end do
                    splt2t(nj) = sumtemp / qfrac2
                    splt2e(nj) = sumelev / qfrac2
                  end do
                end if
                j2max   = 1
                j2min   = 1
                maxtemp = splt2t(1)
                mintemp = splt2t(1)
                do n=2,ng2
                  if (splt2t(n) > maxtemp) then
                    maxtemp = splt2t(n)
                    j2max   = n
                  else if (splt2t(n) < mintemp) then
                    mintemp = splt2t(n)
                    j2min   = n
                  end if
                end do
                if (ttarg < etemp1 - 0.001) then                             ! need a colder temp from group 2
                  if (maxtemp - mintemp > 0.001) then
                    etemp2 = splt2t(j2min)                                   ! preferred outlet is the coldest one
                    elev2  = splt2e(j2min)
                    j2pref = j2min
                  else
                    etemp2 = splt2t(j2lo)                                    ! preferred outlet is the lowest one
                    elev2  = splt2e(j2lo)
                    j2pref = j2lo
                  end if
                else if (ttarg > etemp1 + 0.001) then                        ! need a warmer temp from group 2
                  if (maxtemp - mintemp > 0.001) then
                    etemp2 = splt2t(j2max)                                   ! preferred outlet is the warmest one
                    elev2  = splt2e(j2max)
                    j2pref = j2max
                  else
                    etemp2 = splt2t(j2hi)                                    ! preferred outlet is the highest one
                    elev2  = splt2e(j2hi)
                    j2pref = j2hi
                  end if
                else
                  etemp2 = splt2t(1)                                         ! if temp is close to target, choose first outlet
                  elev2  = splt2e(1)
                  j2pref = 1
                end if
              end if

            ! Target temperature is less than either outlet temperature.
              if (ttarg < etemp1 .and. ttarg < etemp2) then
                qfrac1 = sum_minfrac1                                        ! default for if/then cases
                if (abs(etemp1 - etemp2) < 0.001) then
                  if (prior1 == prior2) then                                 ! If each outlet has the same priority level, then...
                    if (elev1 <= elev2) then                                 ! Choose lower outlet if both have same temperature.
                      qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2
                    end if
                  else if (prior1 < prior2) then
                    qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2               ! Choose higher priority outlet if temps are the same.
                  end if
                else if (etemp1 < etemp2) then                               ! If temps are different, choose the one closer
                  qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2                 ! to target temperature
                end if

            ! Target temperature is greater than either outlet temperature.
              else if (ttarg > etemp1 .and. ttarg > etemp2) then
                qfrac1 = sum_minfrac1                                        ! default for if/then cases
                if (abs(etemp1 - etemp2) < 0.001) then
                  if (prior1 == prior2) then                                 ! If each outlet has the same priority level, then...
                    if (elev1 >= elev2) then                                 ! Choose upper outlet if both have same temperature.
                      qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2
                    end if
                  else if (prior1 < prior2) then
                    qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2               ! Choose higher priority outlet if temps are the same.
                  end if
                else if (etemp1 > etemp2) then                               ! If temps are different, choose the one closer
                  qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2                 ! to target temperature
                end if

            ! Target temperature is essentially the same as the two outlet temperatures.
              else if (abs(etemp1 - etemp2) < 0.001) then
                qfrac1 = sum_minfrac1                                        ! default for if/then cases
                if (prior1 == prior2) then                                   ! If each outlet has the same priority level, then...
                  qfrac1 = sum_minfrac1 + 0.5 * (1.0 - sum_minfrac0 - sum_minfrac1 - sum_minfrac2)    ! Split the flow equally.
                else if (prior1 < prior2) then
                  qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2                 ! Choose higher priority outlet if temps are the same.
                end if

            ! Target temperature is between the two outlet temperatures.
              else
                qfrac1 = (1.0 - sum_minfrac0) * abs((ttarg-etemp2)/(etemp1-etemp2+NONZERO))
                qfrac1 = max(sum_minfrac1, qfrac1)
                qfrac1 = min(1.0 - sum_minfrac0 - sum_minfrac2, qfrac1)
              end if
              qfrac1 = min(sum_maxfrac1, qfrac1)
              qfrac2 = 1.0 - sum_minfrac0 - qfrac1
              if (qfrac2 > sum_maxfrac2) then
                excess_frac = qfrac2 - sum_maxfrac2
                if (qfrac1 + excess_frac <= sum_maxfrac1) then
                  qfrac1 = qfrac1 + excess_frac
                else
                  qfrac1 = sum_maxfrac1
                end if
                qfrac2 = sum_maxfrac2
              end if

            ! Set flow fractions for individual outlets and redistribute flows if maximum flow fractions exceeded.
              call Set_Flow_Fracs(j, jb, j2pref)

            ! Exit the loop if the latest flow fraction calculation agrees with the previous one.
            ! Exit if no flow, because no iteration requried in that case.
              if (abs(lastfrac - qfrac1) < tsconv .or. qall == 0.0) exit
            end do

          ! Check to see if iterative solution did not converge.
            if (abs(lastfrac - qfrac1) >= tsconv .and. qall > 0.0) then
              write (wrn,'(A,F0.3,3(A,F0.4))') 'Flow fraction calculations not converging at day ', jday,                          &
                                               '  Current: ', qfrac1, ' Last: ', lastfrac, ' Next-to-last: ', lastfrac2
              WARNING_OPEN = .TRUE.

            ! Check to see if the iterative solution is unstable.  If so, use priorities to assign releases.
            ! Criteria:  change is at least 0.1 and most recent change is in opposite direction from previous change.
              if (abs(lastfrac - qfrac1) >= 0.1 .and. (qfrac1-lastfrac)*(lastfrac-lastfrac2) < 0.0) then
                qfrac1 = sum_minfrac1                                        ! default for if/then cases
                if (prior1 < prior2) then                                    ! group 1 is higher priority
                  qfrac1 = 1.0 - sum_minfrac0 - sum_minfrac2
                else                                                         ! else, fulfill minima and split the rest
                  qfrac1 = sum_minfrac1 + 0.5 * (1.0 - sum_minfrac0 - sum_minfrac1 - sum_minfrac2)
                end if
                qfrac1 = min(sum_maxfrac1, qfrac1)
                qfrac2 = 1.0 - sum_minfrac0 - qfrac1
                if (qfrac2 > sum_maxfrac2) then
                  excess_frac = qfrac2 - sum_maxfrac2
                  if (qfrac1 + excess_frac <= sum_maxfrac1) then
                    qfrac1 = qfrac1 + excess_frac
                  else
                    qfrac1 = sum_maxfrac1
                  end if
                  qfrac2 = sum_maxfrac2
                end if
                call Set_Flow_Fracs(j, jb, j2pref)                           ! set flow fractions; redistribute if maxfrac exceeded
              end if
            end if
          end if

        ! Set final flows to go with the flow fractions.  May not be necessary, but do it anyway.
          do jj=1,nouts(j)
            qwd(jstsplt(j,jj)) = qall * qwdfrac(jstsplt(j,jj))
          end do
        end if
      end if
    end do

    nxtsplit = nxtsplit + tspltfreq
  end if

! Use the flow fractions to set flows in blended groups.
  if (tspltc=='      ON') then
    do j=1,numtsplt
      if (tsyearly(j) == '     OFF') then
        daytest = jday
      else
        daytest = real(jdayg) + jday - int(jday)
      end if
      if (daytest >= tstsrt(j) .and. daytest < tstend(j)) then
        qall = 0.0

      ! Do structures first
        if (tspltcntr(j) == '      ST') then
          do jj=1,nouts(j)
            qall = qall + qstr(jstsplt(j,jj),tspltjb(j))                            ! sum up all the flows
          end do
          do jj=1,nouts(j)                                                          ! set the flows and honor the maximum flow
            jst = jstsplt(j,jj)
            qstr(jst,tspltjb(j)) = qstrfrac(jst,tspltjb(j)) * qall
            if (tsmaxflow(j,jj) > 0.0 .and. qstr(jst,tspltjb(j)) > tsmaxflow(j,jj)) qstr(jst,tspltjb(j)) = tsmaxflow(j,jj)
          end do

      ! Do Withdrawals next
        else if (tspltcntr(j) == '      WD') then
          do jj=1,nouts(j)
            qall = qall + qwd(jstsplt(j,jj))                                        ! sum up all the flows
          end do
          do jj=1,nouts(j)                                                          ! set the flows and honor the maximum flow
            jwd = jstsplt(j,jj)
            qwd(jwd) = qwdfrac(jwd) * qall
            if (tsmaxflow(j,jj) > 0.0 .and. qwd(jwd) > tsmaxflow(j,jj)) qwd(jwd) = tsmaxflow(j,jj)
          end do
        end if
      end if
    end do
  end if

! Output some results.
  if (jday.ge.nxtstr) then
    nxtstr = nxtstr+tfrqtmp
    ifile=1949
    do jb=1,nbr
      if (nstr(jb) > 0) then
        ifile=ifile+1
        write (ifile,'(f10.4,",",<nstr(jb)>(f10.2,","),<nstr(jb)>(f10.2,","),<nstr(jb)>(f10.2,","))') jday,(tavg(i,jb),i=1,nstr(jb)),(qstr(i,jb),i=1,nstr(jb)),(estr(i,jb),i=1,nstr(jb))                   ! SW 8/28/2019
      end if
    end do
    if (nwd > 0) then
      ifile=ifile+1
      write (ifile,'(f10.4,<nwd>f10.2,<nwd>f10.2,<nwd>f10.2)') jday,(tavgw(i),i=1,nwd),(qwd(i),i=1,nwd),(ewd(i),i=1,nwd)
    end if

  ! computing reservoir volume and volume below 'tempcrit'
    volmc=0.0
    volm=0.0
    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        DO I=cus(jb),ds(jb)
          volm(jw) = volm(jw) +BH2(KT,I)*DLX(I)
          DO K=kt+1,kb(i)
            volm(jw) = volm(jw)+BH(K,I)*DLX(I)
          END DO
          do kk=1,tempn
            if(t2(kt,i).le.tempcrit(jw,kk))volmc(jw,kk) = volmc(jw,kk)+BH2(KT,I)*DLX(I)
            DO K=kt+1,kb(i)
              if(t2(k,i).le.tempcrit(jw,kk))volmc(jw,kk) = volmc(jw,kk)+BH(K,I)*DLX(I)
            END DO
          end do
        end do
      end do

      ifile=ifile+1
      write(ifile,5315)jday,volm(jw),(volmc(jw,kk), kk=1,tempn)
5315  format(f8.2,100(g12.4,g12.4))
    end do
  end if

! Check elevations and status of temperature control towers.
  if (tempc == '      ON' .and. jday .ge. nxttcd) then

  ! Update the temperature targets
    DO J=1,NUMTEMPC
      IF (DYNSEL(J) == '      ON') THEN
        SELD(J)=1009+J
        DO WHILE (JDAY >= NXSEL(J))
          tctemp(j) = TEMP2(J)
             IF(DYNSF(J))THEN
                    READ (SELD(J),*) NXSEL(J),TEMP2(J)
                ELSE
                    READ (SELD(J),'(1000F8.0)') NXSEL(J),TEMP2(J)
                ENDIF   
          !READ (SELD(J),'(1000F8.0)') NXSEL(J),TEMP2(J)
        END DO
      END IF
    END DO

    do j=1,numtempc

    ! Structures
      if (tcntr(j) == '      ST') then
        js = tcjs(j)                                 ! set structure index
        jb = tcjb(j)                                 ! set branch index
        DO JW=1,NWB                                  ! set waterbody index
          IF (JB >= BS(JW) .AND. JB <= BE(JW)) EXIT
        END DO

        if (tciseg(j) .eq. 0) then
          tcomp = tavg(js,jb)               !cb 9/8/06
        else if (tciseg(j) < 0) then
          tcomp = twdo(abs(tciseg(j)))      ! sw 11/26/10
        else

        ! Check to see if the monitoring segment tciseg is in the same branch and waterbody as the structure
          DO JJB=1,NBR
            IF (tciseg(j) >= US(JJB) .AND. tciseg(j) <= DS(JJB)) exit
          end do
          DO JJW=1,NWB
            IF (JjB >= BS(JJW) .AND. JJB <= BE(JJW)) EXIT
          END DO

          IF (tcklay(j) < 0) THEN
            K = INT(ABS(tcklay(j)))
          ELSE
            DO K=KTWB(JJW),KB(tciseg(j))
              IF (DEPTHB(K,tciseg(j)) > tcklay(j)) EXIT
            END DO
            K = MIN(K,KB(tciseg(j)))
          END IF
          tcomp = t2(k,tciseg(j))
        end if
        if (tcyearly(j) == '     OFF') then
          daytest = jday
        else
          daytest = real(jdayg) + jday - int(jday)
        end if
        if (daytest >= tctsrt(j) .and. daytest < tctend(j)) then
          estr(js,jb) = tcelev(j,ncountc(js,jb))                  ! initialize the structure elevation

          if (tcomp > tctemp(j) .and. tcnelev(j) > ncountc(js,jb)) then
          ! making sure that the next lower structure for a particular 'j' is found
            do nj=ncountc(js,jb)+1,tcnelev(j)
              if (tcelev(j,nj) < estr(js,jb)) then
                ncountc(js,jb) = nj
                estr(js,jb)    = tcelev(j,ncountc(js,jb))
                exit
              end if
            end do

          else if (tcomp < tctemp(j) .and. ncountc(js,jb) .gt. 1) then
          ! to prevent this happening at each time it checks it and hence oscillating back and forth - check the temp at the upper outlet also
            if (tciseg(j) > 0) then
              if (jb .eq. jjb) then
                wsel = elws(ds(jb)) - sina(jb) * dlx(ds(jb)) * 0.5                        ! compute water-surface elevation !SR 03/24/13
                do ks=ktwb(jw),kb(ds(jb))
                ! if (depthb(ks,tciseg(j)) > tcelev(j,ncountc(js,jb)-1)) exit             !??can't be right-- SR 03/24/13
                  if (wsel - depthb(ks,tciseg(j)) < tcelev(j,ncountc(js,jb)-1)) exit      !SR 03/24/13
                end do
                ks   = min(ks,kb(tciseg(j)))
                tmod = t2(ks,ds(jb))
              else
                tmod = t2(k,tciseg(j))
              end if
              if (tmod < tctemp(j) .and. tcelev(j,ncountc(js,jb)-1) < elws(ds(jb))) then
              ! making sure that the next upper structure for a particular 'j' is found
                do nj=ncountc(js,jb)-1,1,-1
                  if (tcelev(j,nj) > estr(js,jb)) then
                    ncountc(js,jb) = nj
                    estr(js,jb)    = tcelev(j,ncountc(js,jb))
                    exit
                  end if
                end do
              end if

            else if (tciseg(j) .eq. 0) then
            ! calculate the estimated outflow temperature at higher ports when tcomp < tctemp(j),
            ! and move up if higher port still meets criteria - this doesn't happen when tciseg < 0
              do nj=1,ncountc(js,jb)-1
                id = ds(jb)
                kt = ktwb(jw)
                call downstream_withdrawal_estimate(js,tempest,tcelev(j,nj))
                if (tempest < tctemp(j) .and. tcelev(j,nj) < elws(ds(jb))) then
                  ncountc(js,jb) = nj
                  estr(js,jb)    = tcelev(j,ncountc(js,jb))
                  exit
                end if
              end do
            end if
          end if
          if (tcelevcon(j) == '      ON' .and. tcnelev(j) > ncountc(js,jb) .and. estr(js,jb) > elws(ds(jb))) then
            ncountc(js,jb) = ncountc(js,jb)+1
            estr(js,jb)    = tcelev(j,ncountc(js,jb))
          end if
        end if

    ! Withdrawals
      else if (tcntr(j) == '      WD') then
        jwd = tcjs(j)
        if (tciseg(j) .eq. 0) then
        ! tcomp = tout(jb)
          tcomp = tavgw(tcjs(j))   !cb 9/8/06
        else if (tciseg(j) < 0) then
          tcomp = twdo(abs(tciseg(j)))
        else

        ! checking to see if the monitoring segment tciseg is in the same branch and water body as the withdrawal
          DO JJB=1,NBR
            IF (tciseg(j) >= US(JJB) .AND. tciseg(j) <= DS(JJB)) exit
          end do
          DO JJW=1,NWB
            IF (JjB >= BS(JJW) .AND. JJB <= BE(JJW)) EXIT
          END DO

          IF (tcklay(j) < 0) THEN
            K = INT(ABS(tcklay(j)))
          ELSE
            DO K=KTWB(JJW),KB(tciseg(j))
              IF (DEPTHB(K,tciseg(j)) > tcklay(j)) EXIT
            END DO
            K = MIN(K,KB(tciseg(j)))
          END IF
          tcomp = t2(k,tciseg(j))
        end if
        if (tcyearly(j) == '     OFF') then
          daytest = jday
        else
          daytest = real(jdayg) + jday - int(jday)
        end if
        if (daytest >= tctsrt(j) .and. daytest < tctend(j)) then
          ewd(jwd) = tcelev(j,ncountcw(jwd))                  ! initialize the withdrawal elevation

          if (tcomp > tctemp(j) .and. tcnelev(j) > ncountcw(jwd)) then
          ! making sure that the next lower structure for a particular 'j' is found
            do nj=ncountcw(jwd)+1,tcnelev(j)
              if (tcelev(j,nj) < ewd(jwd)) then
                ncountcw(jwd) = nj
                ewd(jwd)      = tcelev(j,ncountcw(jwd))
                exit
              end if
            end do

          else if (tcomp < tctemp(j) .and. ncountcw(jwd) .gt. 1) then
          ! to prevent this happening at each time it checks it and hence oscillating back and forth - check the temp at the upper outlet also
            if (tciseg(j) > 0) then
              tmod = t2(k,tciseg(j))
              if (tmod < tctemp(j) .and. tcelev(j,ncountcw(jwd)-1) < elws(iwd(jwd))) then
              ! making sure that the next upper structure for a particular 'j' is found
                do nj=ncountcw(jwd)-1,1,-1
                  if (tcelev(j,nj) > ewd(jwd)) then
                    ncountcw(jwd) = nj
                    ewd(jwd)      = tcelev(j,ncountcw(jwd))
                    exit
                  end if
                end do
              end if

            else if (tciseg(j) == 0) then
            ! calculate estimated outflow temperature at higher ports when tcomp < tctemp(j),
            ! and move up if higher port still meets criteria
              I = MAX(CUS(JBWD(JWD)),IWD(JWD))
              DO JJB=1,NBR
                IF (I >= US(JJB) .AND. I <= DS(JJB)) exit
              end do
              DO JJW=1,NWB
                IF (JJB >= BS(JJW) .AND. JJB <= BE(JJW)) EXIT
              END DO
              kt = ktwb(jjw)
              do nj=1,ncountcw(jwd)-1
                call lateral_withdrawal_estimate(jwd,tempest,tcelev(j,nj))
                if (tempest < tctemp(j) .and. tcelev(j,nj) < elws(iwd(jwd))) then
                  ncountcw(jwd) = nj
                  ewd(jwd)      = tcelev(j,ncountcw(jwd))
                  exit
                end if
              end do
            end if
          end if
          if (tcelevcon(j) == '      ON' .and. tcnelev(j) > ncountcw(jwd) .and. ewd(jwd) > elws(iwd(jwd))) then
            ncountcw(jwd) = ncountcw(jwd)+1
            ewd(jwd)      = tcelev(j,ncountcw(jwd))
          end if
        end if
      end if
    end do

    nxttcd = nxttcd + tcdfreq
  end if
return
		  
		  
ENTRY DEALLOCATE_SELECTIVEUSGS
  DEAllocate (tcnelev,tcjb,tcjs, tcelev,tctemp,tctend,tctsrt,ncountc,tciseg,tcklay,tcelevcon,elcontspl)
  DEAllocate (tspltjb,tspltt,nouts,jstsplt,kstrsplt,tcyearly, tcntr,tspltcntr)
  DEallocate (volm,ncountcw,qwdfrac,qstrfrac)
  DEallocate (tempcrit,volmc,DYNSEL,SELD,NXSEL,TEMP2,TSYEARLY,TSTEND,TSTSRT)
  deallocate (tsdepth, tstype, tsminfrac, tsprior, tsminhead, tsmaxhead, tsmaxflow, no_flow)
  deallocate (tsdynsel, tsseld, nxtssel, tstemp2, ewdsav, estrsav, share_flow, wd_active, str_active)
  deallocate (nout0, nout1, nout2, minfrac1, minfrac2, maxfrac1, maxfrac2, splt2t, splt2e,DYNSF)
RETURN

End Subroutine SelectiveUSGS


!***********************************************************************************************************************************
!**                                                S E T _ F L O W _ F R A C T I O N S                                            **
!***********************************************************************************************************************************
                                                                              ! Entire routine added/modified by S. Rounds, 06/26/13
Subroutine Set_Flow_Fractions
  Use Selective1USGS                                                          ! This routine sets the flow fractions, and then
  IMPLICIT NONE                                                               ! redistributes flow to other outlets in the same
  integer :: n, nj, nexcess, j, jb, jst, jwd, j2pref                          ! group when one or more outlets exceed their maximum
  real    :: excess_frac, addfrac                                             ! flow rates.  Excess flow that cannot be accommodated
Return                                                                        ! within each group will be discarded.


Entry Set_Flow_Fracs(j, jb, j2pref)
  excess_frac = 0.0                                                           ! Excess flow above maximum release rates
  nexcess     = 0                                                             ! Number of group 1 outlets exceeding maximum rates

! Set and rebalance release fractions for group 1 structures
  if (tspltcntr(j) == '      ST') then
    do n=1,ng1                                                                ! Find maxed-out outlets in group 1; set flows to max
      jst = jstsplt(j,nout1(n))
      qstrfrac(jst,jb) = minfrac1(n) + (qfrac1 - sum_minfrac1) / ng1
      if (qstrfrac(jst,jb) > maxfrac1(n)) then
        nexcess = nexcess + 1
        excess_frac = excess_frac + qstrfrac(jst,jb) - maxfrac1(n)
        qstrfrac(jst,jb) = maxfrac1(n)
      end if
    end do
    if (excess_frac > 0.0 .and. ng1 - nexcess > 0) then                       ! Redistribute excess flow to other outlets in group
      do nj=1,ng1                                                             ! Iterative process, in case others get maxed-out
        if (ng1 == nexcess .or. excess_frac <= 0.00001) exit
        addfrac = excess_frac / (ng1 - nexcess)
        do n=1,ng1
          jst = jstsplt(j,nout1(n))
          if (maxfrac1(n) - qstrfrac(jst,jb) > 0.00001) then
            if (qstrfrac(jst,jb) + addfrac > maxfrac1(n)) then
              nexcess = nexcess + 1
              excess_frac = excess_frac - (maxfrac1(n) - qstrfrac(jst,jb))
              qstrfrac(jst,jb) = maxfrac1(n)
            else
              excess_frac = excess_frac - addfrac
              qstrfrac(jst,jb) = qstrfrac(jst,jb) + addfrac
            end if
          end if
        end do
      end do
    end if

! Set and rebalance release fractions for group 1 withdrawals
  else
    do n=1,ng1                                                                ! Find maxed-out outlets in group 1; set flows to max
      jwd = jstsplt(j,nout1(n))
      qwdfrac(jwd) = minfrac1(n) + (qfrac1 - sum_minfrac1) / ng1
      if (qwdfrac(jwd) > maxfrac1(n)) then
        nexcess = nexcess + 1
        excess_frac = excess_frac + qwdfrac(jwd) - maxfrac1(n)
        qwdfrac(jwd) = maxfrac1(n)
      end if
    end do
    if (excess_frac > 0.0 .and. ng1 - nexcess > 0) then                       ! Redistribute excess flow to other outlets in group
      do nj=1,ng1                                                             ! Iterative process, in case others get maxed-out
        if (ng1 == nexcess .or. excess_frac <= 0.00001) exit
        addfrac = excess_frac / (ng1 - nexcess)
        do n=1,ng1
          jwd = jstsplt(j,nout1(n))
          if (maxfrac1(n) - qwdfrac(jwd) > 0.00001) then
            if (qwdfrac(jwd) + addfrac > maxfrac1(n)) then
              nexcess = nexcess + 1
              excess_frac = excess_frac - (maxfrac1(n) - qwdfrac(jwd))
              qwdfrac(jwd) = maxfrac1(n)
            else
              excess_frac = excess_frac - addfrac
              qwdfrac(jwd) = qwdfrac(jwd) + addfrac
            end if
          end if
        end do
      end do
    end if
  end if

Entry Set_Flow_Fracs2(j, jb, j2pref)                                          ! Separate entry just for group 2 outlets

  if (j2pref == 0) j2pref = 1                                                 ! Preferred outlet number, if not sharing
  excess_frac = 0.0                                                           ! Excess flow above maximum release rates
  nexcess     = 0                                                             ! Number of group 2 outlets exceeding maximum rates

! Set and rebalance release fractions for group 2 structures
  if (tspltcntr(j) == '      ST') then
    do n=1,ng2                                                                ! Find maxed-out outlets in group 2; set flows to max
      jst = jstsplt(j,nout2(n))
      if (.not. share_flow(j) .and. ng2 > 1) then                             ! Direct flow to preferred outlet if not shared
        if (n == j2pref) then
          qstrfrac(jst,jb) = qfrac2 - sum_minfrac2 + minfrac2(n)
        else
          qstrfrac(jst,jb) = minfrac2(n)
        end if
      else
        qstrfrac(jst,jb) = minfrac2(n) + (qfrac2 - sum_minfrac2) / ng2
      end if
      if (qstrfrac(jst,jb) > maxfrac2(n)) then
        nexcess = nexcess + 1
        excess_frac = excess_frac + qstrfrac(jst,jb) - maxfrac2(n)
        qstrfrac(jst,jb) = maxfrac2(n)
      end if
    end do
    if (excess_frac > 0.0 .and. ng2 - nexcess > 0) then                       ! Redistribute excess flow to other outlets in group
      do nj=1,ng2                                                             ! Iterative process, in case others get maxed-out
        if (ng2 == nexcess .or. excess_frac <= 0.00001) exit
        addfrac = excess_frac / (ng2 - nexcess)
        do n=1,ng2
          jst = jstsplt(j,nout2(n))
          if (maxfrac2(n) - qstrfrac(jst,jb) > 0.00001) then
            if (qstrfrac(jst,jb) + addfrac > maxfrac2(n)) then
              nexcess = nexcess + 1
              excess_frac = excess_frac - (maxfrac2(n) - qstrfrac(jst,jb))
              qstrfrac(jst,jb) = maxfrac2(n)
            else
              excess_frac = excess_frac - addfrac
              qstrfrac(jst,jb) = qstrfrac(jst,jb) + addfrac
            end if
          end if
        end do
      end do
    end if

! Set and rebalance release fractions for group 2 withdrawals
  else if (tspltcntr(j) == '      WD') then
    do n=1,ng2                                                                ! Find maxed-out outlets in group 2; set flows to max
      jwd = jstsplt(j,nout2(n))
      if (.not. share_flow(j) .and. ng2 > 1) then                             ! Direct flow to preferred outlet if not shared
        if (n == j2pref) then
          qwdfrac(jwd) = qfrac2 - sum_minfrac2 + minfrac2(n)
        else
          qwdfrac(jwd) = minfrac2(n)
        end if
      else
        qwdfrac(jwd) = minfrac2(n) + (qfrac2 - sum_minfrac2) / ng2
      end if
      if (qwdfrac(jwd) > maxfrac2(n)) then
        nexcess = nexcess + 1
        excess_frac = excess_frac + qwdfrac(jwd) - maxfrac2(n)
        qwdfrac(jwd) = maxfrac2(n)
      end if
    end do
    if (excess_frac > 0.0 .and. ng2 - nexcess > 0) then                       ! Redistribute excess flow to other outlets in group
      do nj=1,ng2                                                             ! Iterative process, in case others get maxed-out
        if (ng2 == nexcess .or. excess_frac <= 0.00001) exit
        addfrac = excess_frac / (ng2 - nexcess)
        do n=1,ng2
          jwd = jstsplt(j,nout2(n))
          if (maxfrac2(n) - qwdfrac(jwd) > 0.00001) then
            if (qwdfrac(jwd) + addfrac > maxfrac2(n)) then
              nexcess = nexcess + 1
              excess_frac = excess_frac - (maxfrac2(n) - qwdfrac(jwd))
              qwdfrac(jwd) = maxfrac2(n)
            else
              excess_frac = excess_frac - addfrac
              qwdfrac(jwd) = qwdfrac(jwd) + addfrac
            end if
          end if
        end do
      end do
    end if
  end if
Return
End Subroutine Set_Flow_Fractions