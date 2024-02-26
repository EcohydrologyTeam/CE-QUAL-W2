
!***********************************************************************************************************************************
!**                                        S U B R O U T I N E   C A L C U L A T E  A Z                                           **
!***********************************************************************************************************************************
MODULE AZ_LOCAL
  USE GEOMC; USE GLOBAL; USE TRANS; USE PREC; USE EDDY; USE KINETIC; use macrophytec; USE LOGICC; USE MAIN, ONLY:WARNING_OPEN
  REAL(R8), DIMENSION (2) :: SIG
  REAL(R8)            :: TKEMIN1=1.25D-7,TKEMIN2=1.0D-9,RIAZ0,AZ0,RIAZ1,EXPAZ,DEPTHL,DEPTHR,ZDLR,DEPTH,USTAR,ZD
  REAL(R8)            :: SLM,VISCK,VISCF,HRAD,EFFRIC,GC2,RI,BUOY,USTARBKT,WALLFUNCTION,UDR,UDL,CONVKX,CONVEX,CONVKZ,CONVEZ,ATA,AB
  INTEGER             :: K,KL,J,IUT,IDT,JJ
  REAL(R8)            :: BOUK,PRDK,PRHE,PRHK,UNST,UNSE,USTARB
  EXTERNAL WALLFUNCTION
END MODULE AZ_LOCAL

  SUBROUTINE CALCULATE_AZ
 ! USE GEOMC; USE GLOBAL; USE TRANS; USE PREC; USE EDDY; USE KINETIC; use macrophytec; USE LOGICC; USE MAIN, ONLY:WARNING_OPEN
  USE AZ_LOCAL
  IMPLICIT NONE
  SAVE
  SIG(1)=1.0D0; SIG(2)=1.3D0

  IF     (AZC(JW) == '     TKE') THEN
    CALL CALCULATE_TKE
  ELSEIF (AZC(JW) == '    TKE1') THEN
    CALL CALCULATE_TKE1
  ELSE
    DO K=KT,KBMIN(I)-1
    VSH(K,I) = ((U(K+1,I)-U(K,I))/((AVH2(K,I)+AVH2(K,I+1))*0.5))**2
    END DO

    DO K=KT,KBMIN(I)-1
      CALL CALCULATE_AZ0
      BUOY    = (RHO(K+1,I)-RHO(K,I)+RHO(K+1,I+1)-RHO(K,I+1))/(AVH2(K,I)+AVH2(K,I+1))            ! 2.0*AVH2(K,I)
      RIAZ0   =  DLOG(AZ0/AZMAX(JW))*0.666666667D0                                                  ! /1.5
      RI      =  G*BUOY/(RHOW*VSH(K,I)+NONZERO)
      RIAZ1   =  DMAX1(RI,RIAZ0)
      RIAZ1   =  DMIN1(RIAZ1,10.0D0)
      EXPAZ   =  DEXP(-1.5D0*RIAZ1)
      AZ(K,I) = DMAX1(AZMIN,AZ0*EXPAZ+AZMIN*(1.0D0-EXPAZ))              ! AZ computed at lower edge of cell
      DZT(K,I) = DMAX1(DZMIN,FRAZDZ*(AZ0*EXPAZ+DZMIN*(1.0-EXPAZ)))     ! DZ computed at lower edge of cell - later averaged to cell center lower edge
    END DO
  END IF
RETURN
END SUBROUTINE CALCULATE_AZ
    
SUBROUTINE CALCULATE_AZ0
  USE AZ_LOCAL
  IMPLICIT NONE
  SAVE

  IF (K == KT) THEN
    DEPTHL = (ELWS(I)-EL(KB(I),I)    +H2(KB(I),I)  *COSA(JB))/COSA(JB)      !(EL(KT,I)  -Z(I)  *COSA(JB)-EL(KB(I),I)    +H2(KB(I),I)  *COSA(JB))/COSA(JB)
    DEPTHR = (ELWS(I+1)-EL(KB(I+1),I+1)+H2(KB(I+1),I)*COSA(JB))/COSA(JB)    !(EL(KT,I+1)-Z(I+1)*COSA(JB)-EL(KB(I+1),I+1)+H2(KB(I+1),I)*COSA(JB))/COSA(JB)
    ZDLR   =  DEPTHL-H2(KT,I)+DEPTHR-H2(KT,I+1)
    IF (AZC(JW) == '     RNG' .OR. AZC(JW) == '   PARAB') THEN
      DEPTH = (DEPTHR+DEPTHL)*0.5
      USTAR =  SQRT(G*DEPTH*SLOPEC(JB))
      IF (SLOPEC(JB) == 0.0) THEN
        USTAR = 0.0
        DO KL=KT,KBMIN(I)
          USTAR = USTAR+SQRT(AZ(KL-1,I)*SQRT(VSH(KL-1,I))/RHO(KL,I))
        END DO
        USTAR = USTAR/(KBMIN(I)-KT+1)
      END IF
    END IF
  ELSE
    ZDLR = (EL(K-1,I)-EL(KB(I),I)+H2(KB(I),I)*COSA(JB))/COSA(JB)+(EL(K-1,I+1)-EL(KB(I),I+1)+H2(KB(I+1),I)*COSA(JB))/COSA(JB)
  END IF
  ZD = ZDLR/(DEPTHL+DEPTHR)
  IF (AZC(JW) == '    NICK') THEN
    SLM = (DEPTH*(0.14-0.08*(1.0-ZD)**2-0.06*(1.0-ZD)**4))**2
    AZ0 =  MAX(AZMIN,SLM*SQRT(VSH(K,I)))
  ELSE IF (AZC(JW) == '     RNG') THEN
    VISCK = DEXP((T2(K,I)+495.691)/(-37.3877))
    IF (T2(K,I) > 30.0)  VISCK = DEXP((T2(K,I)+782.190)/(-57.7600))
    VISCF = MAX(0.0,0.08477*((ZDLR*0.5*USTAR/VISCK)**3)*((1.0-ZDLR*0.5/DEPTH)**3)-100.0)
    VISCF = (1.0+VISCF)**0.33333333333333
    AZ0   = MAX(AZMIN,VISCK*VISCF)
  ELSE IF (AZC(JW) == '   PARAB') THEN
    AZ0 = MAX(AZMIN,0.41*USTAR*ZDLR*0.5*(1.0-ZD))
  ELSE
    SLM = HMAX2
    IF (AZC(JW) == '     W2N') SLM = ((DEPTHR+DEPTHL)*0.5*(0.14-0.08*(1.0-ZD)**2-0.06*(1.0-ZD)**4))**2
    AZ0 = 0.4*SLM*SQRT(VSH(K,I)+((FRICBR(K,I)+WSHY(I)*DECAY(K,I))/(AZ(K,I)+NONZERO))**2)+AZMIN
  END IF
RETURN
END SUBROUTINE CALCULATE_AZ0

SUBROUTINE CALCULATE_TKE
  USE AZ_LOCAL
  IMPLICIT NONE
  SAVE
  USTAR = SQRT(1.25*CZ(I)*WIND10(I)**2/RHO(KT,I))
  IF (MANNINGS_N(JW)) THEN
   HRAD = BH1(KT,I)/(B(KTI(I),I)-B(KT+1,I)+2.*AVH1(KT,I))      ! HRAD = BHR1(KT,I)/(BR(KTI(I),I)-BR(KT+1,I)+2.*AVH1(KT,I))  SW 10/5/07  These calculations are at the segment centers and vertical center of a layer
    if(macrophyte_on.and.mannings_n(jw))then
      call macrophyte_friction(hrad,fric(i),effric,kt,i)
      gc2=g*effric*effric/hrad**0.33333333
    else if(.not.macrophyte_on.and.mannings_n(jw))then
      gc2=g*fric(i)*fric(i)/hrad**0.33333333
    end if
  ELSE
    GC2 = 0.0
    IF (FRIC(I) /= 0.0) GC2 = G/(FRIC(I)*FRIC(I))
  END IF
  USTARBKT     = SQRT(GC2)*ABS(0.5*(U(KT,I)+U(KT,I-1)))                                                 ! SG 10/4/07
  TKE(KT,I,1) = (3.33*(USTAR*USTAR+USTARBKT*USTARBKT))*(BH2(KT,I)/BH1(KT,I))                            ! SG 10/4/07
  TKE(KT,I,2) = (USTAR*USTAR*USTAR+USTARBKT*USTARBKT*USTARBKT)*5.0/H1(KT,I)*(BH2(KT,I)/BH1(KT,I))       ! SG 10/4/07
  DO K=KT+1,KB(I)-1
    BOUK = MAX(AZ(K,I)*G*(RHO(K+1,I)-RHO(K,I))/(H(K,JW)*RHOW),0.0)
    PRDK = AZ(K,I)*(0.5*(U(K,I)+U(K,I-1)-U(K+1,I)-U(K+1,I-1))/(H(K,JW)*0.5 + H(K+1,JW)*0.5))**2.0       ! SG 10/4/07
    PRHE = 10.0*GC2**1.25*ABS(0.5*(U(K,I)+U(K,I-1)))**4.0/(0.5*B(K,I))**2.0
    IF (MANNINGS_N(JW)) THEN
      HRAD = BH(K,I)/(B(K,I)-B(K+1,I)+2.0*H(K,JW))          ! HRAD = BHR(K,I)/(BR(K,I)-BR(K+1,I)+2.0*H(K,JW))  SW 10/5/07
      if(macrophyte_on.and.mannings_n(jw))then
        call macrophyte_friction(hrad,fric(i),effric,k,i)
        gc2=g*effric*effric/hrad**0.33333333
      else if (.not.macrophyte_on.and.mannings_n(jw))then
        gc2=g*fric(i)*fric(i)/hrad**0.33333333
      end if
    END IF
    PRHK       = GC2/(0.5*B(K,I))*ABS(0.5*(U(K,I)+U(K,I-1)))**3.0
    UNST       = PRDK-TKE(K,I,2)
    UNSE       = 1.44*TKE(K,I,2)/TKE(K,I,1)*PRDK-1.92*(TKE(K,I,2)/TKE(K,I,1)*TKE(K,I,2))
    TKE(K,I,1) = TKE(K,I,1)+DLT*(UNST+PRHK-BOUK)
    TKE(K,I,2) = TKE(K,I,2)+DLT*(UNSE+PRHE)
  END DO
  USTARB         = SQRT(GC2)*ABS(0.5*(U(KB(I),I)+U(KB(I),I-1)))
  TKE(KB(I),I,1) = 0.5*(3.33*USTARB*USTARB+TKE(KB(I),I,1))
  TKE(KB(I),I,2) = 0.5*(USTARB*USTARB*USTARB*5.0/H(KB(I),JW)+TKE(KB(I),I,2))

  DO J=1,2                                                    ! SG 10/4/07 Series of bug fixes for TKE
	  K=KT			
	  AT(K,I) = 0.0
	  CT(K,I) = 0.0
	  VT(K,I) = 1.0
	  DT(K,I) = TKE(K,I,J)	
	  DO K = KT+1,KB(I)-1
      AT(K,I) = -DLT/BH1(K,I)*BB(K-1,I)/SIG(J)*AZ(K-1,I)/AVH1(K-1,I)
      CT(K,I) = -DLT/BH1(K,I)*BB(K,I)/SIG(J)*AZ(K,I)/AVH1(K,I)
      VT(K,I) =  1.0-AT(K,I)-CT(K,I)
      DT(K,I) =  TKE(K,I,J)
    END DO
	  K = KB(I)
	  AT(K,I) = 0.0
	  CT(K,I) = 0.0
	  VT(K,I) = 1.0
	  DT(K,I) = TKE(K,I,J)
	  CALL TRIDIAG(AT(:,I),VT(:,I),CT(:,I),DT(:,I),KT,KB(I),KMX,TKE(:,I,J))
  END DO

  DO K=KT,KB(I)
    TKE(K,I,1) = MAX(TKE(K,I,1),TKEMIN1)
    TKE(K,I,2) = MAX(TKE(K,I,2),TKEMIN2)
    AZT(K,I)     = 0.09*TKE(K,I,1)*TKE(K,I,1)/TKE(K,I,2)
  END DO

  DO K=KT,KB(I)-1
    AZ(K,I) = 0.5*(AZT(K,I)+AZT(K+1,I))
    AZ(K,I) = MAX(AZMIN,AZ(K,I))
    AZ(K,I) = MIN(AZMAX(JW),AZ(K,I))
    DZ(K,I) = MAX(DZMIN,FRAZDZ*AZ(K,I))            ! No need to average DZ further since defined at cell, center bottom. AZ needs to be averaged to RHS of cell.
  END DO
  AZ(KB(I),I) = AZMIN
  AZT(KB(I),I)= AZMIN
  RETURN
    END SUBROUTINE CALCULATE_TKE
    
SUBROUTINE CALCULATE_TKE1
  USE AZ_LOCAL
  IMPLICIT NONE
  SAVE
! Subroutine based on work of Gould(2006) RODI WITHOUT WIND 			
! SW 6-28-04 DEADSEA  TKE
! BASED ON SUBROUTINE FIRST DEVELOPED BY CHAPMAN AND COLE (1995)
! BUG FIXES, CORRECTIONS, ADDITIONS SCOTT WELLS (2001)

! BOTTOM LAYER BC
  IF(ABS((U(KB(I),I)+U(KB(I),I-1))*.5)>NONZERO)THEN
	USTARB = WALLFUNCTION()
! ERROR TRAPPING
    IF(USTARB == 0.0)THEN
!	        IF (MANNINGS_N(JW)) THEN
!	        HRAD = BH1(KT,I)/(B(KTI(I),I)-B(KT+1,I)+2.*AVH1(KT,I))
!            if(macrophyte_on.and.mannings_n(jw))then
!                call macrophyte_friction(hrad,fric(i),effric,kt,i)
!                gc2=g*effric*effric/hrad**0.33333333
!            else if(.not.macrophyte_on.and.mannings_n(jw))then
!                gc2=g*fric(i)*fric(i)/hrad**0.33333333
!            end if
!            ELSE
!            GC2 = 0.0
!            IF (FRIC(I) /= 0.0) GC2 = G/(FRIC(I)*FRIC(I))
!            END IF
!	    USTARB = SQRT(GC2)*ABS(0.5*(U(KB(I),I)+U(KB(I),I-1)))
! END ERROR TRAPPING	
	ENDIF
	TKE(KB(I),I,1) = (3.33*USTARB*USTARB)
	TKE(KB(I),I,2) = (USTARB*USTARB*USTARB)/(0.41*(H(KB(I),JW))*0.5)                      !/2.0
    USTARBTKE(I) = USTARB
	IF(STRICKON(JW))THEN
	  CALL CALCFRIC()
	ENDIF
  ELSE
	USTARB = 0
	TKE(KB(I),I,1) = TKEMIN1
	TKE(KB(I),I,2) = TKEMIN2
    USTARBTKE(I) = USTARB
  ENDIF

!INLET SET TO MOLECULAR CONDITIONS
  IUT = FIRSTI(JW)
  IF(I==IUT)THEN
	DO K= KT, KB(I)
	  TKE(K,IUT-1,2) = TKEMIN2
	  TKE(K,IUT-1,1) = TKEMIN1
	ENDDO
  ENDIF
!OUTLET SET TO EQUAL TO I=I_OUTLET - 1
  IDT = LASTI(JW)
  IF(I==IDT-1)THEN
    DO K= KT, KB(I)
	  TKE(K,IDT,2) = TKE(K,IDT-1,2)
	  TKE(K,IDT,1) = TKE(K,IDT-1,1)
	ENDDO
  ENDIF

!CALCULATE HORIZONTAL PORTION OF THE SPLIT K-E EQUATIONS FOR INTERIOR LAYERS
  DO K = KT+1,KB(I)-1		
	IF (MANNINGS_N(JW)) THEN
	   HRAD = BH1(KT,I)/(B(KTI(I),I)-B(KT+1,I)+2.D0*AVH1(KT,I))                                !HRAD = BHR(K,I)/(BR(K,I)-BR(K+1,I)+2.0*H(K,JW))
	   GC2  = G*FRIC(I)*FRIC(I)/HRAD**0.333
	ELSE
           GC2 = 0.0
           IF (FRIC(I) /= 0.0) GC2 = G/(FRIC(I)*FRIC(I))
	END IF

	BOUK = DMAX1(AZ(K,I)*G*(RHO(K+1,I)-RHO(K,I))/(H(K,JW)*RHOW),0.0D0)
	PRDK = AZ(K,I)*(0.5D0*(U(K,I)+U(K,I-1)-U(K+1,I)-U(K+1,I-1))/(H(K,JW)*0.5D0 + H(K+1,JW)*0.5D0))**2.0                 !/2.0  /2.0
	IF(.NOT.TKELATPRD(JW))THEN
	  PRHE = 0.0
	  PRHK = 0.0
	ELSE
	  PRHE = TKELATPRDCONST(JW)*GC2**1.25*DABS(0.5D0*(U(K,I)+U(K,I-1)))**4.0/(0.5D0*B(K,I))**2.0
      PRHK = GC2/(0.5D0*B(K,I))*DABS(0.5D0*(U(K,I)+U(K,I-1)))**3.0
    ENDIF
	UNST	   = PRDK-TKE(K,I,2)
	UNSE	   = 1.44D0*TKE(K,I,2)/TKE(K,I,1)*PRDK-1.92D0*(TKE(K,I,2)/TKE(K,I,1)*TKE(K,I,2))

	UDR		= (1.0D0+DSIGN(1.0D0,U(K,I)))*0.5D0
	UDL		= (1.0D0+DSIGN(1.0D0,U(K,I-1)))*0.5D0
	CONVKX 	= DLT/(DLX(I)*BH1(K,I))*((U(K,I)*(UDR*TKE(K,I,1)+(1.0-UDR)*TKE(K,I+1,1))*BR(K,I)*H2(K,I)-U(K,I-1)*(UDL*TKE(K,I-1,1)+(1.0-UDL)*TKE(K,I,1))*BR(K,I-1)*H2(K,I-1)))
	CONVEX 	= DLT/(DLX(I)*BH1(K,I))*((U(K,I)*(UDR*TKE(K,I,2)+(1.0-UDR)*TKE(K,I+1,2))*BR(K,I)*H2(K,I)-U(K,I-1)*(UDL*TKE(K,I-1,2)+(1.0-UDL)*TKE(K,I,2))*BR(K,I-1)*H2(K,I-1)))
	IF(IMPTKE(JW) == '     IMP')THEN
	  TKE(K,I,1) = TKE(K,I,1) - CONVKX + DLT*(UNST+PRHK-BOUK)		
	  TKE(K,I,2) = TKE(K,I,2) - CONVEX + DLT*(UNSE+PRHE)
	ELSE
	  ATA		  = (1.D00+DSIGN(1.0D0,W(K-1,I)))*0.5D0
	  AB		  = (1.0D0+DSIGN(1.0D0,W(K,I)))*0.5D0
	
	  CONVKZ	  = DLT/(H2(K,I)*BH1(K,I))*((AB*TKE(K,I,1)+(1.0-AB)*TKE(K+1,I,1))*W(K,I)*BB(K,I)*AVH2(K,I)-(ATA*(TKE(K-1,I,1)+(1.0-ATA)*TKE(K,I,1))*W(K-1,I)*BB(K-1,I)*AVH2(K-1,I)))
      CONVEZ	  = DLT/(H2(K,I)*BH1(K,I))*((AB*TKE(K,I,2)+(1.0-AB)*TKE(K+1,I,2))*W(K,I)*BB(K,I)*AVH2(K,I)-(ATA*(TKE(K-1,I,2)+(1.0-ATA)*TKE(K,I,2))*W(K-1,I)*BB(K-1,I)*AVH2(K-1,I)))

	  TKE(K,I,1) = TKE(K,I,1) - CONVKX - CONVKZ + DLT*(UNST+PRHK-BOUK)		
	  TKE(K,I,2) = TKE(K,I,2) - CONVEX - CONVEZ + DLT*(UNSE+PRHE)
	ENDIF
  END DO

!CALCULATE TOP LAYER BOUNDARY CONDITION AND SOLVE THE K-E EQUATION USING TRIDIAGONAL SOLVER
500 IF(TKEBC(JW)==1)THEN
!RODI NOWIND
    DEPTHL = (ELWS(I)-EL(KB(I),I)	  +H2(KB(I),I)	*COSA(JB))/COSA(JB)                                !EL(KT,I)  -Z(I)  *COSA(JB)
    DEPTHR = (ELWS(I+1)-EL(KB(I+1),I+1)+H2(KB(I+1),I)*COSA(JB))/COSA(JB)                                !EL(KT,I+1)-Z(I+1)*COSA(JB)
    DEPTH = (DEPTHR+DEPTHL)*0.5
	DO J=1,2
	  K=KT	
	  IF(J==1)THEN		
		AT(K,I) = 0.0
		CT(K,I) =-1.0
		VT(K,I) = 1.0
		DT(K,I) = 0.0
	  ELSE
		TKE(KT,I,2) = TKE(K,I,1)**(1.5)/(ARODI(JW)*DEPTH)*BH2(KT,I)/BH1(KT,I)                                     !3.0/2.0
		AT(K,I) = 0.0
		CT(K,I) = 0.0
		VT(K,I) = 1.0
		DT(K,I) = TKE(KT,I,2)
	  ENDIF
	  DO K = KT+1,KB(I)-1
		 IF(IMPTKE(JW)=='     IMP')THEN
		   AT(K,I) =  -DLT/BH1(K,I)*(BB(K-1,I)/SIG(J)*AZ(K-1,I)/AVH1(K-1,I) + 0.5*W(K-1,I)*BB(K-1,I))
		   CT(K,I) =  -DLT/BH1(K,I)*(BB(K,I)/SIG(J)*AZ(K,I)/AVH1(K,I) - 0.5*W(K,I)*BB(K,I))
		   VT(K,I) =  1 +&
				      DLT/BH1(K,I)*(BB(K,I)/SIG(J)*AZ(K,I)/AVH1(K,I) +&
				      BB(K-1,I)/SIG(J)*AZ(K-1,I)/AVH1(K-1,I) +&
					  0.5*(W(K,I)*BB(K,I)-W(K-1,I)*BB(K-1,I)))
		   DT(K,I) =  TKE(K,I,J)
		 ELSE
		   AT(K,I) = -DLT/BH1(K,I)*BB(K-1,I)/SIG(J)*AZ(K-1,I)/AVH1(K-1,I)
		   CT(K,I) = -DLT/BH1(K,I)*BB(K,I)/SIG(J)*AZ(K,I)/AVH1(K,I)
		   VT(K,I) =  1.0-AT(K,I)-CT(K,I)
		   DT(K,I) =  TKE(K,I,J)
		 ENDIF
	  END DO
	  K = KB(I)
	  AT(K,I) = 0.0
	  CT(K,I) = 0.0
	  VT(K,I) = 1.0
	  DT(K,I) = TKE(K,I,J)
	  CALL TRIDIAG(AT(:,I),VT(:,I),CT(:,I),DT(:,I),KT,KB(I),KMX,TKE(:,I,J))
	END DO
  ELSEIF(TKEBC(JW)==2)THEN
!RODI WIND
    DEPTHL = (EL(KT,I)  -Z(I)  *COSA(JB)-EL(KB(I),I)	  +H2(KB(I),I)	*COSA(JB))/COSA(JB)
    DEPTHR = (EL(KT,I+1)-Z(I+1)*COSA(JB)-EL(KB(I+1),I+1)+H2(KB(I+1),I)*COSA(JB))/COSA(JB)
    DEPTH = (DEPTHR+DEPTHL)*0.5
	DO JJ=1,3
	  IF(JJ==1)THEN
		J=1
		DO K = KT+1,KB(I)-1
		  TKE(K,I,3) = TKE(K,I,1)
		ENDDO
	  ELSEIF(JJ==2)THEN
		J=1
	  ELSEIF(JJ==3)THEN
		J=2
	  ENDIF
	  ! WIND SHEAR AT SURFACE
	  USTAR = SQRT(1.25*CZ(I)*WIND10(I)**2/RHO(KT,I))	
	  IF(JJ==1)THEN 	
		AT(KT,I) = 0.0
		CT(KT,I) =-1.0
		VT(KT,I) = 1.0
		DT(KT,I) = 0.0
	  ELSEIF(JJ==2)THEN
		IF(TKE(KT,I,1)*0.3<(USTAR*USTAR))THEN
		  TKE(KT,I,1) = (USTAR*USTAR)*3.33333*BH2(KT,I)/BH1(KT,I)
		  AT(KT,I) = 0.0
		  CT(KT,I) = 0.0
		  VT(KT,I) = 1.0
		  DT(KT,I) = TKE(KT,I,J)
		  DO K = KT+1,KB(I)-1
			TKE(K,I,1) = TKE(K,I,3)
		  ENDDO
		ELSE
		  AT(KT,I) = 0.0
		  CT(KT,I) =-1.0
		  VT(KT,I) = 1.0
		  DT(KT,I) = 0.0		
		ENDIF
	  ENDIF
	  IF(JJ==3)THEN
		TKE(KT,I,2) = (TKE(KT,I,1)*0.3)**(1.5)/(0.41*(H2(KT,I)*0.5 + ARODI(JW)*DEPTH*   &
		              (1-(USTAR*USTAR)/(TKE(KT,I,1)*0.3))))*BH2(KT,I)/BH1(KT,I)		   !3.0/2.0    /2.0
		AT(KT,I) = 0.0
		CT(KT,I) = 0.0
		VT(KT,I) = 1.0
		DT(KT,I) = TKE(KT,I,2)
	  ENDIF
	  DO K = KT+1,KB(I)-1
		 IF(IMPTKE(JW)=='     IMP')THEN
		   AT(K,I) =  -DLT/BH1(K,I)*(BB(K-1,I)/SIG(J)*AZ(K-1,I)/AVH1(K-1,I) + 0.5*W(K-1,I)*BB(K-1,I))
		   CT(K,I) =  -DLT/BH1(K,I)*(BB(K,I)/SIG(J)*AZ(K,I)/AVH1(K,I) - 0.5*W(K,I)*BB(K,I))
		   VT(K,I) =  1 +&
				      DLT/BH1(K,I)*(BB(K,I)/SIG(J)*AZ(K,I)/AVH1(K,I) +&
				      BB(K-1,I)/SIG(J)*AZ(K-1,I)/AVH1(K-1,I) +&
					  0.5*(W(K,I)*BB(K,I)-W(K-1,I)*BB(K-1,I)))
		   DT(K,I) =  TKE(K,I,J)
		 ELSE
		   AT(K,I) = -DLT/BH1(K,I)*BB(K-1,I)/SIG(J)*AZ(K-1,I)/AVH1(K-1,I)
		   CT(K,I) = -DLT/BH1(K,I)*BB(K,I)/SIG(J)*AZ(K,I)/AVH1(K,I)
		   VT(K,I) =  1.0-AT(K,I)-CT(K,I)
		   DT(K,I) =  TKE(K,I,J)
		 ENDIF
	  END DO
	  K = KB(I)
	  AT(K,I) = 0.0
	  CT(K,I) = 0.0
	  VT(K,I) = 1.0
	  DT(K,I) = TKE(K,I,J)
	  CALL TRIDIAG(AT(:,I),VT(:,I),CT(:,I),DT(:,I),KT,KB(I),KMX,TKE(:,I,J))
	END DO
  ELSEIF(TKEBC(JW)==3)THEN
!CEQUALW2 WIND
! TOP LAYER BC
  USTAR = SQRT(1.25*CZ(I)*WIND10(I)**2/RHO(KT,I))
  IF (MANNINGS_N(JW)) THEN
	HRAD = BH1(KT,I)/(B(KTI(I),I)-B(KT+1,I)+2.*AVH1(KT,I))                   !HRAD = BHR1(KT,I)/(BR(KTI(I),I)-BR(KT+1,I)+2.*AVH1(KT,I))
	GC2  = G*FRIC(I)*FRIC(I)/HRAD**0.333
  ELSE
	GC2 = 0.0
	IF (FRIC(I) /= 0.0) GC2 = G/(FRIC(I)*FRIC(I))
  END IF
  USTARBKT	  = SQRT(GC2)*ABS(0.5*(U(KT,I)+U(KT,I-1)))
  TKE(KT,I,1) = (3.33*(USTAR*USTAR+USTARBKT*USTARBKT))*(BH2(KT,I)/BH1(KT,I))
  TKE(KT,I,2) = (USTAR*USTAR*USTAR+USTARBKT*USTARBKT*USTARBKT)*5.0/H1(KT,I)*(BH2(KT,I)/BH1(KT,I))
	DO J=1,2
	  K=KT			
	  AT(K,I) = 0.0
	  CT(K,I) = 0.0
	  VT(K,I) = 1.0
	  DT(K,I) = TKE(K,I,J)	
	  DO K = KT+1,KB(I)-1
		 IF(IMPTKE(JW)=='     IMP')THEN
		   AT(K,I) =  -DLT/BH1(K,I)*(BB(K-1,I)/SIG(J)*AZ(K-1,I)/AVH1(K-1,I) + 0.5*W(K-1,I)*BB(K-1,I))
		   CT(K,I) =  -DLT/BH1(K,I)*(BB(K,I)/SIG(J)*AZ(K,I)/AVH1(K,I) - 0.5*W(K,I)*BB(K,I))
		   VT(K,I) =  1 +&
				      DLT/BH1(K,I)*(BB(K,I)/SIG(J)*AZ(K,I)/AVH1(K,I) +&
				      BB(K-1,I)/SIG(J)*AZ(K-1,I)/AVH1(K-1,I) +&
					  0.5*(W(K,I)*BB(K,I)-W(K-1,I)*BB(K-1,I)))
		   DT(K,I) =  TKE(K,I,J)
		 ELSE
		   AT(K,I) = -DLT/BH1(K,I)*BB(K-1,I)/SIG(J)*AZ(K-1,I)/AVH1(K-1,I)
		   CT(K,I) = -DLT/BH1(K,I)*BB(K,I)/SIG(J)*AZ(K,I)/AVH1(K,I)
		   VT(K,I) =  1.0-AT(K,I)-CT(K,I)
		   DT(K,I) =  TKE(K,I,J)
		 ENDIF
	  END DO
	  K = KB(I)
	  AT(K,I) = 0.0
	  CT(K,I) = 0.0
	  VT(K,I) = 1.0
	  DT(K,I) = TKE(K,I,J)
	  CALL TRIDIAG(AT(:,I),VT(:,I),CT(:,I),DT(:,I),KT,KB(I),KMX,TKE(:,I,J))
	END DO
  ELSE
    WARNING_OPEN    = .TRUE.
    WRITE(WRN,*)'UNKNOWN TKEBC OPTION CHECK W2_CON.NPT'
    WRITE(WRN,*)'Setting TKE calculation to TKEBC(JW)=3'
    TKEBC(JW)=3
    GOTO 500
  ENDIF

  DO K=KT,KB(I)
    TKE(K,I,1) = MAX(TKE(K,I,1),TKEMIN1)
    TKE(K,I,2) = MAX(TKE(K,I,2),TKEMIN2)
    AZT(K,I)     = 0.09*TKE(K,I,1)*TKE(K,I,1)/TKE(K,I,2)
  END DO

  DO K=KT,KB(I)-1
	AZ(K,I) = 0.5*(AZT(K,I)+AZT(K+1,I))
	AZ(K,I) = MAX(AZMIN,AZ(K,I))
	AZ(K,I) = MIN(AZMAX(JW),AZ(K,I))
	DZ(K,I) = MAX(DZMIN,FRAZDZ*AZ(K,I))      ! No need to average further, DZ is now defined at cell bottom, AZ will need to be averaged to the RHS of cell
  END DO
  AZ(KB(I),I) = AZMIN
  AZT(KB(I),I)= AZMIN
RETURN
END SUBROUTINE CALCULATE_TKE1

REAL(R8) FUNCTION WALLFUNCTION()
USE GLOBAL;USE GEOMC;USE EDDY;USE LOGICC;USE MAIN, ONLY:WARNING_OPEN; USE SCREENC, ONLY:JDAY
IMPLICIT NONE
EXTERNAL SEMILOG, RTBIS
REAL(R8) :: KS, UBOUND,LBOUND,VISCK,SEMILOG,RTBIS,UST,PERIMETER,AREA,BS1                                  !,USTARB
REAL(R8) :: V
INTEGER :: K,MULTIPLIER
    VISCK = DEXP(-(T2(KB(I),I)+495.691)*0.026746764)                                                     !/(-37.3877)
    IF (T2(KB(I),I) > 30.0)  VISCK = DEXP(-(T2(KB(I),I)+782.190)*0.01731301939)                          !/(-57.7600)
  	UBOUND = ABS((U(KB(I),I)+U(KB(I),I-1))*.5)
    IF(SEMILOG(UBOUND,ABS((U(KB(I),I)+U(KB(I),I-1))*.5),(H(KB(I),JW))/2.0,VISCK)>0) THEN               !FRIC(I),
	  MULTIPLIER = 2
	  DO WHILE(SEMILOG(UBOUND,ABS((U(KB(I),I)+U(KB(I),I-1))*.5),(H(KB(I),JW))/2.0,VISCK)>0)            !,FRIC(I)
  		UBOUND = ABS((U(KB(I),I)+U(KB(I),I-1))*.5)*REAL(MULTIPLIER)
		MULTIPLIER = MULTIPLIER + 1
		IF(MULTIPLIER>30)THEN
		  WARNING_OPEN    = .TRUE.
		  WRITE(WRN,*)'WALLFUNCTION: UPPER BOUND NOT FOUND IN TKE1 ROUTINE ON JDAY:',JDAY
		  WRITE(WRN,*)'Setting TKE calculation to TKE from TKE1'
		  AZC(JW)='     TKE'
		  WALLFUNCTION=0.0
		  RETURN
		ENDIF
	  ENDDO
	ENDIF
	LBOUND = VISCK/(E(I)*(H(KB(I),JW))/2.0D0)+NONZERO
	MULTIPLIER = 1
	DO WHILE(SEMILOG(LBOUND,DABS((U(KB(I),I)+U(KB(I),I-1))*.5D0),(H(KB(I),JW))/2.0D0,VISCK)<0)
 		LBOUND = VISCK/(E(I)*(H(KB(I),JW))/2.0)+10**MULTIPLIER*NONZERO
		MULTIPLIER = MULTIPLIER + 1
		IF(MULTIPLIER>30)THEN
		  WARNING_OPEN    = .TRUE.
		  WRITE(WRN,*)'WALLFUNCTION: LOWER BOUND NOT FOUND IN TKE1 ROUTINE ON JDAY:',JDAY
		  WRITE(WRN,*)'Setting TKE calculation to TKE from TKE1'
		  AZC(JW)='     TKE'
		  WALLFUNCTION=0.0
		  RETURN
		ENDIF
	ENDDO    	
    WALLFUNCTION = RTBIS(LBOUND,UBOUND,1D-15,DABS((U(KB(I),I)+U(KB(I),I-1))*.5D0),(H(KB(I),JW))/2.0D0,VISCK)           !,FRIC(I)
    IF(STRICKON(JW))THEN
	  V = VISCK
	  UST = WALLFUNCTION
	  IF (MANNINGS_N(JW)) THEN
        KS = ((STRICK(JW)*FRIC(I))**6)
	  ELSE
        PERIMETER = B(KB(I),I)+2.0D0*H(KB(I),JW)
        DO K=KT,KB(I)-1
          PERIMETER = PERIMETER + (B(K,I)-B(K+1,I)+2.0D0*H(K,JW))
        END DO
        AREA = 0
        DO K=KT,KB(I)
          AREA = AREA + BH1(K,I)
        END DO
	    KS = ((STRICK(JW)*1.0D0/FRIC(I)*(AREA/PERIMETER)**(0.16666666667D0))**6)        !1.0/6.0
	  ENDIF
      IF(UST*KS/V>1.0)THEN
        BS1 = (5.5D0+2.5D0*DLOG(UST*KS/V))*DEXP(-0.217D0*DLOG(UST*KS/V)**2)+8.5D0*(1-DEXP(-0.217D0*DLOG(UST*KS/V)**2))
        E(I) = DEXP((.41*BS1)/(UST*KS/V))
      ELSE
        E(I) = 9.535D0
      ENDIF	
	ENDIF
ENDFUNCTION WALLFUNCTION

FUNCTION SEMILOG(UST,URES,Y,V)
USE EDDY;USE GLOBAL
IMPLICIT NONE
REAL(R8), INTENT(IN) :: UST,URES,Y,V
REAL(R8) :: SEMILOG
  SEMILOG = 0.41*URES/DLOG(E(I)*Y*UST/V) - UST
END FUNCTION SEMILOG

FUNCTION RTBIS(X1,X2,XACC,U,Y,V)
USE GLOBAL, ONLY:WRN,JW; USE EDDY, ONLY:AZC; USE SCREENC, ONLY:JDAY; USE PREC
IMPLICIT NONE
REAL(R8), INTENT(IN) :: X1,X2,XACC,U,Y,V
REAL(R8) :: RTBIS,SEMILOG
EXTERNAL SEMILOG
INTEGER, PARAMETER :: MAXIT=70                                                    !80
! USING BISECTION, FIND THE ROOT OF A FUNCTION FUNC KNOWN TO LIE BETWEEN X1 AND X2. THE
! ROOT, RETURNED AS RTBIS, WILL BE REFINED UNTIL ITS ACCURACY IS ±XACC.
! PARAMETER: MAXIT IS THE MAXIMUM ALLOWED NUMBER OF BISECTIONS.
INTEGER  :: J
REAL(R8) :: DX,F,FMID,XMID
  FMID = SEMILOG(X2,U,Y,V)
  F = SEMILOG(X1,U,Y,V)
  IF (F*FMID >= 0.0) THEN
          WRITE(WRN,*)'RTBIS: ROOT MUST BE BRACKETED IN TKE1 ROUTINE ON JDAY:',JDAY
		  WRITE(WRN,*)'Setting TKE calculation to TKE from TKE1'
		  RTBIS=X1
		  AZC(JW)='     TKE'
          RETURN
  ENDIF
  IF (F < 0.0) THEN !ORIENT THE SEARCH SO THAT F>0 LIES AT X+DX.
    RTBIS = X1
    DX = X2-X1
  ELSE
    RTBIS = X2
    DX = X1-X2
  END IF
  DO J=1,MAXIT !BISECTION LOOP.
    DX = DX*0.5
    XMID = RTBIS+DX
    FMID = SEMILOG(XMID,U,Y,V)
    IF (FMID <= 0.0) RTBIS=XMID
    IF (ABS(DX) < XACC .OR. FMID == 0.0) THEN
	  RETURN
	ENDIF
  END DO
  WRITE(WRN,*)'RTBIS: TOO MANY BISECTIONS IN TKE1 ROUTINE ON JDAY:',JDAY
  WRITE(WRN,*)'Setting TKE calculation to TKE from TKE1'
  AZC(JW)='     TKE'
END FUNCTION RTBIS

SUBROUTINE CALCFRIC()
USE GLOBAL;USE GEOMC;USE EDDY;USE KINETIC;USE LOGICC
IMPLICIT NONE
REAL(R8)    :: DEPTHL,DEPTHR,DEPTH,PERIMETER,AREA,UAVG
INTEGER :: K
  DEPTHL = (ELWS(I)-EL(KB(I),I)    +H2(KB(I),I)  *COSA(JB))/COSA(JB)                                              !EL(KT,I)  -Z(I)  *COSA(JB)
  DEPTHR = (ELWS(I+1)-EL(KB(I+1),I+1)+H2(KB(I+1),I)*COSA(JB))/COSA(JB)                                                     !EL(KT,I+1)-Z(I+1)*COSA(JB)
  DEPTH = (DEPTHR+DEPTHL)*0.5
  PERIMETER = 0.0
  PERIMETER = BR(KB(I),I)+2.0*H(KB(I),JW)
  DO K=KT,KB(I)-1
    PERIMETER = PERIMETER + (B(K,I)-B(K+1,I)+2.0*H(K,JW))
  END DO
  AREA = 0.0
  DO K=KT,KB(I)
    AREA = AREA + BH1(K,I)
  END DO
  UAVG = 0.5*(QC(I)+QC(I-1))/AREA
  IF(ABS(UAVG)>NONZERO)THEN
    IF (MANNINGS_N(JW)) THEN
	  FRIC(I) = (AREA/PERIMETER)**(0.16666666667)/SQRT(9.81*(UAVG/USTARBTKE(I))**2)                      !1.0/6.0   (9.81/(USTARBTKE(I)/UAVG)**2)
	ELSE
	  FRIC(I) = SQRT(9.81*(UAVG/USTARBTKE(I))**2)
	ENDIF
  ENDIF
END SUBROUTINE



