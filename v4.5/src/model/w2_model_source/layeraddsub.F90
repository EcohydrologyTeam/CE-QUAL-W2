SUBROUTINE LAYERADDSUB
USE MAIN
USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART
  USE MACROPHYTEC; USE POROSITYC; USE ZOOPLANKTONC
 IMPLICIT NONE

INTEGER :: KTMAX,JBIZ,KKB,I_BR_NUM
REAL(R8):: TMAC
REAL(R8):: W1,W2,W3, DUMMY

!***********************************************************************************************************************************
!**                                       Task 2.5: Layer - Segment Additions and Subtractions                                    **
!***********************************************************************************************************************************

!** Water surface minimum thickness

    DO JW=1,NWB
      KT       =  KTWB(JW)
      ZMIN(JW) = -1000.0
      KTMAX    =  2                                                                                                 ! SR 10/17/05
      DO JB=BS(JW),BE(JW)
         IF(BR_INACTIVE(JB))THEN 
          IF(DS(JB)-US(JB)+1>=3)THEN
              I_BR_NUM=DS(JB)-2
          ELSE
              I_BR_NUM=DS(JB)-1    ! FOR BRANCHES WITH LESS THAN 3 SEGMENTS
          ENDIF
          
          IF(CUS(JBDH(JB))<=DHS(JB) .AND. ELWS(DHS(JB))>EL(KB(I_BR_NUM),I_BR_NUM))THEN     ! ***
              BR_INACTIVE(JB)=.FALSE.
               IF (SNAPSHOT(JW)) WRITE (SNP(JW),'(/1X,13("*"),1X,A,I0,A,F0.3,A,I0,1X,A,I0,13("*"))') '   Branch Active: ',jb,&
                                                        ' at Julian day = ',JDAY,'   NIT = ',NIT 
               WRITE(WRN,'(1X,13("*"),1X,A,I0,A,F0.3,A,I0)') '   Branch Active: ',jb,' at Julian day = ',JDAY,'   NIT = ',NIT       ! SW 7/24/2018
               IF (SNAPSHOT(JW)) WRITE (SNP(JW),'(/17X,2(A,I0))') ' Add segments ',DS(JB)-1,' through ',DS(JB)
               WRITE (WRN,'(/17X,2(A,I0))') ' Add segments ',DS(JB)-1,' through ',DS(JB)
               CUS(JB)=DS(JB)-1
                        DO I=DS(JB)-1,DS(JB)                 
                          Z(I)         =  Z(DHS(JB))
                          KTI(I)       =  KTI(DHS(JB))
                          H1(KT+1,I)   =  H(KT+1,JW)
                          AVH1(KT+1,I) = (H1(KT+1,I)+H1(KT+2,I))*0.5
                          AVHR(KT+1,I) =  H1(KT+1,I)+(H1(KT+1,I+1)-H1(KT+1,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                  !SW 07/29/04
                          IF (.NOT. TRAPEZOIDAL(JW)) THEN
                            BH1(KT+1,I) =  B(KT+1,I)*H(KT+1,JW)
                            H1(KT,I)    =  H(KT,JW)-Z(I)
                            BI(KT:KB(I),I) =  B(KT:KB(I),I)                                                                        ! SW 4/18/07
                            BI(KT,I)    =  B(KTI(I),I)
                            AVH1(KT,I)  = (H1(KT,I)+H1(KT+1,I))*0.5
                            AVHR(KT,I)  =  H1(KT,I)+(H1(KT,I+1)-H1(KT,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                       !SW 07/29/04
                            BH1(KT,I)   =  BI(KT,I)*(EL(KT,I)-Z(I)*COSA(JB)-EL(KTI(I)+1,I))/COSA(JB)
                            IF (KTI(I) >= KB(I)) BH1(KT,I) = B(KT,I)*H1(KT,I)
                            DO K=KTI(I)+1,KT
                              BH1(KT,I) = BH1(KT,I)+BH(K,I)
                            END DO
                          ELSE
                            CALL GRID_AREA1 (EL(KT,I)-Z(I),EL(KT+1,I),BH1(KT,I),BI(KT,I))                                                    !SW 08/03/04
                            BH1(KT+1,I) =  0.25*H(KT+1,JW)*(BB(KT,I)+2.*B(KT+1,I)+BB(KT+1,I))
                            H1(KT,I)    =  H(KT,JW)-Z(I)
                            AVH1(KT,I)  = (H1(KT,I)+H1(KT+1,I))*0.5
                            AVHR(KT,I)  =  H1(KT,I)+(H1(KT,I+1)-H1(KT,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                       !SW 07/29/04
                          END IF
                          BKT(I)      = BH1(KT,I)/H1(KT,I)
                          DEPTHB(KT,I) = H1(KT,I)
                          DEPTHM(KT,I) = H1(KT,I)*0.5
                          DO K=KT+1,KB(I)
                            DEPTHB(K,I) = DEPTHB(K-1,I)+ H1(K,I)
                            DEPTHM(K,I) = DEPTHM(K-1,I)+(H1(K-1,I)+H1(K,I))*0.5
                          END DO
                        END DO
                        DO I=DS(JB)-1,DS(JB)
                          BHR1(KT+1,I) = BH1(KT+1,I)+(BH1(KT+1,I+1)-BH1(KT+1,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                !SW 07/29/04
                          BHR1(KT,I)   = BH1(KT,I)  +(BH1(KT,I+1)  -BH1(KT,I))  /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                !SW 07/29/04
                                      IF(CONSTRICTION(KT,I))THEN    ! SW 6/26/2018
                                        IF(BHR1(KT,I) > BCONSTRICTION(I)*H1(KT,I))BHR1(KT,I)= BCONSTRICTION(I)*H1(KT,I)
                                        IF(BHR1(KT+1,I) > BCONSTRICTION(I)*H1(KT+1,I))BHR1(KT+1,I)= BCONSTRICTION(I)*H1(KT+1,I)
                                      ENDIF
                        END DO
                        DO I=DS(JB)-1,DS(JB)
                          WIND2(I) = WIND2(DHS(JB))
                          IF (DYNAMIC_SHADE(I)) CALL SHADING
                          DO K=KT,KB(I)
                             U(K,I)  = 0.0
                            SDKV(K,I)=SDK(JW)                            
                            IF(DXI(JW) >= 0.0)THEN
                                 DX(K,I) = DXI(JW)
                            ELSE
                                 DX(K,I) = ABS(U(K,I))*ABS(DXI(JW))*H(K,JW)   ! SW 8/2/2017
                            ENDIF
                            IF (INTERNAL_WEIR(K,I)) THEN
                              DX(K,I) = 0.0
                              U(K,I)  = 0.0
                            END IF
                            T1(K,I)           = T1(K,DHS(JB))
                            T2(K,I)           = T1(K,DHS(JB))
                            SU(K,I)           = 0.0
                            C1(K,I,CN(1:NAC)) = C1(K,DHS(JB),CN(1:NAC))
                            C2(K,I,CN(1:NAC)) = C1(K,DHS(JB),CN(1:NAC))
                            DO JE=1,NEP
                              EPD(K,I,JE) = 0.01
                              EPC(K,I,JE) = 0.01/H1(K,I)
                            END DO
                            CMBRT(CN(1:NAC),JB) = CMBRT(CN(1:NAC),JB)+C1(K,DHS(JB),CN(1:NAC))*DLX(I)*BH1(K,I)
                            EBRI(JB)            = EBRI(JB)           +T1(K,DHS(JB))          *DLX(I)*BH1(K,I)
                          END DO
                          DO K=KT,KB(I)-1
                            AZ(K,I)  = AZ(K,DHS(JB))
		                    TKE(K,I,1) = TKE(K,DHS(JB),1) !sg 10/4/07
                            TKE(K,I,2) = TKE(K,DHS(JB),2) !sg 10/4/07
                            SAZ(K,I) = AZ(K,DHS(JB))
                            IF (INTERNAL_WEIR(K,I)) THEN
                              AZ(K,I)  = 0.0
	                          TKE(K,I,1) = 0.0 !sg 10/4/07
                              TKE(K,I,2) = 0.0 !sg 10/4/07
                              SAZ(K,I) = 0.0
                            END IF
                          END DO
                        ENDDO
          ENDIF
          ENDIF
        IF(.NOT.BR_INACTIVE(JB))THEN
        DO I=CUS(JB),DS(JB)
          IF(KB(I) > KTMAX) KTMAX = KB(I)                                                                           ! SR 10/17/05
          IF (Z(I) > ZMIN(JW)) THEN
            IZMIN(JW) = I
            JBIZ      = JB
          END IF
          ZMIN(JW) = MAX(ZMIN(JW),Z(I))
        END DO
        ENDIF
      END DO
      ADD_LAYER = ZMIN(JW) < -0.85*H(KT-1,JW) .AND. KT /= 2
      SUB_LAYER = ZMIN(JW) >  0.60*H(KT,JW)   .AND. KT < KTMAX                                                       ! SR 10/17/05
      IF (KTWB(JW) == KMX-1 .AND. SLOPE(JBIZ) > 0.0 .AND. SUB_LAYER .AND. ONE_LAYER(IZMIN(JW))) THEN
        IF (ZMIN(JW) > 0.99*H(KT,JW)) THEN
            WRITE (WRN,'(A,I0,2(A,F0.3))') 'Low water in segment ',IZMIN(JW),&
                                              ' water surface deviation'//' = ',ZMIN(JW),' at day ',JDAY
            WARNING_OPEN = .TRUE.
        ENDIF
        SUB_LAYER    = .FALSE.
      END IF

      IF(ADD_LAYER == .TRUE. .OR. SUB_LAYER == .TRUE.)THEN
      LAYERCHANGE(JW)=.TRUE.
      ELSE
      LAYERCHANGE(JW)=.FALSE.
      ENDIF

!**** Add layers

      DO WHILE (ADD_LAYER)
        IF (SNAPSHOT(JW)) WRITE (SNP(JW),'(/1X,13("*"),1X,A,I0,A,F0.3,A,I0,1X,A,I0,13("*"))') '   Add layer ',KT-1,&
                                                        ' at Julian day = ',JDAY,'   NIT = ',NIT,' IZMIN =',IZMIN(JW)   ! SW 1/23/06
        WARNING_OPEN = .TRUE.
        WRITE (WRN,'(/1X,13("*"),1X,A,I0,A,F0.3,A,I0,1X,A,I0,13("*"))') '   Add layer ',KT-1,&
                                                        ' at Julian day = ',JDAY,'   NIT = ',NIT,' IZMIN =',IZMIN(JW)   ! SW 1/23/06

!****** Variable initialization

        KTWB(JW) = KTWB(JW)-1
        KT       = KTWB(JW)
        ilayer = 0
        
  ! RECOMPUTE INTERNAL WEIR FOR FLOATING WEIR      
    IF (WEIR_CALC) THEN   !  SW 3/16/18
    DO JWR=1,NIW
     IF(IWR(JWR) >= US(BS(JW)) .AND. IWR(JWR) <= DS(BE(JW)))THEN
        IF (EKTWR(JWR) == 0.0) THEN  
            KTWR(JWR)=KTWB(JW)
        ELSE  
          KTWR(JWR) = INT(EKTWR(JWR))  
        END IF 
        IF (EKBWR(JWR) <= 0.0) THEN  
            DO K=KTWR(JWR),KB(IWR(JWR))  
            IF (DEPTHB(K,IWR(JWR)) > ABS(EKBWR(JWR))) THEN
                KBWR(JWR)=K
                EXIT  
            ENDIF
            END DO   
        ELSE  
          KBWR(JWR) = INT(EKBWR(JWR))  
        END IF  
        
      DO K=2,KMX-1
        IF ((K >= KTWR(JWR) .AND. K <= KBWR(JWR))) INTERNAL_WEIR(K,IWR(JWR)) = .TRUE.
      END DO
     ENDIF     
    END DO    
  END IF
        
        
        DO JB=BS(JW),BE(JW)
        IF(BR_INACTIVE(JB))CYCLE    ! SW 6/12/2017
          IU = CUS(JB)
          ID = DS(JB)
          DO I=IU-1,ID+1
            Z(I)          =  H(KT,JW)+Z(I)
            H1(KT,I)      =  H(KT,JW)-Z(I)
            H1(KT+1,I)    =  H(KT+1,JW)
            H2(KT+1,I)    =  H(KT+1,JW)
            AVH1(KT,I)    = (H1(KT,I)  +H1(KT+1,I))*0.5D0
            AVH1(KT+1,I)  = (H1(KT+1,I)+H1(KT+2,I))*0.5D0
            IF (.NOT. TRAPEZOIDAL(JW)) THEN
              BH1(KT,I)   = BH1(KT+1,I)-Bnew(KT+1,I)*H1(KT+1,I)                              ! SW 1/23/06
              BH1(KT+1,I) = BNEW(KT+1,I)*H1(KT+1,I)                                          ! SW 1/23/06
            ELSE
              CALL GRID_AREA1(EL(KT,I)-Z(I),EL(KT+1,I),BH1(KT,I),DUMMY)                                                          !SW 08/03/04
              BH1(KT+1,I) = 0.25D0*H1(KT+1,JW)*(BB(KT,I)+2.D0*B(KT+1,I)+BB(KT+1,I))
            ENDIF
            VOL(KT,I)     = BH1(KT,I)  *DLX(I)
            VOL(KT+1,I)   = BH1(KT+1,I)*DLX(I)
            BKT(I)        = BH1(KT,I)/H1(KT,I)
            DEPTHB(KT,I)  = H1(KT,I)
            DEPTHM(KT,I)  = H1(KT,I)*0.5
            BI(KT:KB(I),I) =  B(KT:KB(I),I)   ! SW 8/26/05
            BI(KT,I)       =  B(KTI(I),I)
            T1(KT,I)      = T1(KT+1,I)
            SDKV(KT,I)    = SDKV(KT+1,I)      ! SW 1/18/08
            SED(KT,I)     = SED(KT+1,I)       ! SW 1/18/08
            SEDN(KT,I)    = SEDN(KT+1,I)      ! SW 1/18/08
            SEDP(KT,I)    = SEDP(KT+1,I)      ! SW 1/18/08
            SEDC(KT,I)    = SEDC(KT+1,I)      ! SW 1/18/08
!            RHO(KT,I)     = DENSITY(T1(KT,I),MAX(TDS(KT,I),0.0),MAX(TISS(KT,I),0.0))    ! SR 5/15/06
            if(sdfirstadd(kt,i))then
            sed1(kt,i)=sed1ic(kt,i)    ! cb 6/17/17
            sed2(kt,i)=sed2ic(kt,i)    ! cb 6/17/17
              sdfirstadd(kt,i)=.false.
            end if
            DO K=KT+1,KMX
              DEPTHB(K,I) = DEPTHB(K-1,I)+ H1(K,I)
              DEPTHM(K,I) = DEPTHM(K-1,I)+(H1(K-1,I)+H1(K,I))*0.5D0
            END DO
            C1(KT,I,CN(1:NAC))             = C1(KT+1,I,CN(1:NAC))
            CSSK(KT,I,CN(1:NAC))           = CSSK(KT+1,I,CN(1:NAC))
            KF(KT,I,KFCN(1:NAF(JW),JW))    = KF(KT+1,I,KFCN(1:NAF(JW),JW))
            KFS(KT,I,KFCN(1:NAF(JW),JW))   = KFS(KT+1,I,KFCN(1:NAF(JW),JW))                 !KF(KT+1,I,KFCN(1:NAF(JW),JW))   CODE ERROR FIX SW 10/24/2017
            KF(KT+1,I,KFCN(1:NAF(JW),JW))  = 0.0
            KFS(KT+1,I,KFCN(1:NAF(JW),JW)) = 0.0
            IF(KT >= KBI(I))THEN         ! CB 5/24/06
            ADX(KT+1,I)=0.0             ! CB 5/15/06
            C1(KT+1,I,CN(1:NAC))=0.0    ! CB 5/15/06
            CSSK(KT+1,I,CN(1:NAC))=0.0  ! CB 5/15/06
            ENDIF                       ! CB 5/15/06
            DO JE=1,NEP
              if(kt < kbi(i))then    ! CB 4/28/06
              EPD(KT,I,JE)   = EPD(KT+1,I,JE)
              EPM(KT,I,JE)   = EPD(KT,I,JE)*(Bi(KT,I)-B(KT+1,I)+2.0*H1(KT,I))*DLX(I)             ! SR 5/15/06
              EPM(KT+1,I,JE) = EPM(KT+1,I,JE)-EPM(KT,I,JE)
              EPC(KT,I,JE)   = EPM(KT,I,JE)/VOL(KT,I)
              EPC(KT+1,I,JE) = EPM(KT+1,I,JE)/VOL(KT+1,I)
              else
              EPD(KT,I,JE)   = EPD(KT+1,I,JE) ! SW 5/15/06
              EPM(KT,I,JE)   = EPM(KT+1,I,JE)
              EPC(KT,I,JE) =   EPC(KT+1,I,JE)
              EPD(KT+1,I,JE) = 0.0
              EPM(KT+1,I,JE) = 0.0
              EPC(KT+1,I,JE) = 0.0
              END IF                ! CB 4/28/06
            END DO
          END DO

!********macrophytes...
          DO I=IU,ID
            JT=KTI(I)
            JE=KB(I)
            DO J=JT,JE
              IF(J.LT.KT)THEN
                COLB=EL(J+1,I)
              ELSE
                COLB=EL(KT+1,I)
              END IF
         !     COLDEP=ELWS(I)-COLB
               coldep=EL(KT,i)-Z(i)*COSA(JB)-colb  ! cb 3/7/16
         !     MACT(J,KT,I)=MACT(J,KT+1,I)
              IF(MACROPHYTE_ON)MACT(J,KT,I)=MACT(J,KT+1,I)    ! SW 9/28/13
              DO M=1,NMC
                IF(MACROPHYTE_CALC(JW,M))THEN
                  MACRC(J,KT,I,M)=MACRC(J,KT+1,I,M)
                  MACRM(J,KT,I,M)=MACRC(J,KT,I,M)*CW(J,I)*COLDEP*DLX(I)
                END IF
              END DO
            END DO

            JT=KT+1
            JE=KB(I)
            DO J=JT,JE
              DO M=1,NMC
                IF(MACROPHYTE_CALC(JW,M))THEN
                  MACRM(J,KT+1,I,M)=MACRC(J,KT+1,I,M)*CW(J,I)*H(KT+1,JW)*DLX(I)
                END IF
              END DO
            END DO
          END DO
          DO I=IU-1,ID
            AVHR(KT+1,I) =   H1(KT+1,I) +(H1(KT+1,I+1) -H1(KT+1,I)) /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                  !SW 07/29/04 (H1(KT+1,I+1) +H1(KT+1,I))*0.5
            AVHR(KT,I)   =   H1(KT,I)   +(H1(KT,I+1)   -H1(KT,I))   /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                  !SW 07/29/04 (H1(KT,I+1)   +H1(KT,I))*0.5
            BHR1(KT,I)   =   BH1(KT,I)  +(BH1(KT,I+1)  -BH1(KT,I))  /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                  !SW 07/29/04  (BH1(KT,I+1)  +BH1(KT,I))*0.5
            BHR1(KT+1,I) =   BH1(KT+1,I)+(BH1(KT+1,I+1)-BH1(KT+1,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                  !SW 07/29/04  (BH1(KT+1,I+1)+BH1(KT+1,I))*0.5  
                             IF(CONSTRICTION(KT,I))THEN    ! SW 6/26/2018 Valid for all K
                                        IF(BHR1(KT,I) > BCONSTRICTION(I)*H1(KT,I))BHR1(KT,I)= BCONSTRICTION(I)*H1(KT,I)
                                        IF(BHR1(KT+1,I) > BCONSTRICTION(I)*H1(KT+1,I))BHR1(KT+1,I)= BCONSTRICTION(I)*H1(KT+1,I)
                             ENDIF
            U(KT,I)      = U(KT+1,I)
          END DO
          DO I=IU,ID
            IF (ONE_LAYER(I)) THEN
              W(KT,I) = 0.0
            ELSE
              W1      =  W(KT+1,I)*BB(KT+1,I)
              W2      = (BHR1(KT+1,I)*U(KT+1,I)-BHR1(KT+1,I-1)*U(KT+1,I-1))/DLX(I)
              W3      = (-QSS(KT+1,I)*BH1(KT+1,I)/(BH1(KT+1,I)+BH1(KT,I)))/DLX(I)
              W(KT,I) = (W1+W2+W3)/BB(KT,I)
            END IF
          END DO
          IF (UP_HEAD(JB)) THEN
            BHSUM                     = BHR1(KT,IU-1)            +BHR1(KT+1,IU-1)
            QUH1(KT,JB)               = QUH1(KT+1,JB)            *BHR1(KT,IU-1)  /BHSUM
            QUH1(KT+1,JB)             = QUH1(KT+1,JB)            *BHR1(KT+1,IU-1)/BHSUM
            TSSUH1(KT,JB)             = TSSUH1(KT+1,JB)          *BHR1(KT,IU-1)  /BHSUM
            TSSUH1(KT+1,JB)           = TSSUH1(KT+1,JB)          *BHR1(KT+1,IU-1)/BHSUM
            CSSUH1(KT,CN(1:NAC),JB)   = CSSUH1(KT+1,CN(1:NAC),JB)*BHR1(KT,IU-1)  /BHSUM
            CSSUH1(KT+1,CN(1:NAC),JB) = CSSUH1(KT+1,CN(1:NAC),JB)*BHR1(KT+1,IU-1)/BHSUM
          END IF
          IF (DN_HEAD(JB)) THEN
            BHSUM                     = BHR1(KT,ID)              +BHR1(KT+1,ID)
            QDH1(KT,JB)               = QDH1(KT+1,JB)            *BHR1(KT,ID)    /BHSUM
            QDH1(KT+1,JB)             = QDH1(KT+1,JB)            *BHR1(KT+1,ID)  /BHSUM
            TSSDH1(KT,JB)             = TSSDH1(KT+1,JB)          *BHR1(KT,ID)    /BHSUM
            TSSDH1(KT+1,JB)           = TSSDH1(KT+1,JB)          *BHR1(KT+1,ID)  /BHSUM
            CSSDH1(KT,CN(1:NAC),JB)   = CSSDH1(KT+1,CN(1:NAC),JB)*BHR1(KT,ID)    /BHSUM
            CSSDH1(KT+1,CN(1:NAC),JB) = CSSDH1(KT+1,CN(1:NAC),JB)*BHR1(KT+1,ID)  /BHSUM
          END IF
          DO I=IU,ID-1
                IF(DXI(JW)>=0.0)THEN
                        DX(KT,I) = DXI(JW)
                ELSE
                        DX(KT,I) = ABS(U(KT,I))*ABS(DXI(JW))*H(K,JW)    ! SW 8/2/2017
                ENDIF
                IF (INTERNAL_WEIR(KT,I)) DX(KT,I) = 0.0
          END DO
          IUT = IU
          IDT = ID-1
          IF (UP_HEAD(JB)) IUT = IU-1
          IF (DN_HEAD(JB)) IDT = ID
          DO I=IUT,IDT
            AZ(KT,I)  = AZMIN
            TKE(KT,I,1) = 1.25E-7     !sg 10/4/07
            TKE(KT,I,2) = 1.0E-9      !sg 10/4/07
            SAZ(KT,I) = AZMIN
            IF (INTERNAL_WEIR(KT,I)) THEN
              AZ(KT,I)  = 0.0
              TKE(KT,I,1) = 0.0       !sg  10/4/07
              TKE(KT,I,2) = 0.0       !sg  10/4/07
              SAZ(KT,I) = 0.0
            END IF
          END DO
          !IF (CONSTITUENTS) THEN     ! SW 6/26/2019 No need to call twice - just do at end of BRANCH loop
          !  CALL TEMPERATURE_RATES
          !  CALL KINETIC_RATES
          !END IF

!******** Upstream active segment

          IUT = US(JB)
          IF (SLOPE(JB) == 0.0) THEN
            DO I=US(JB),DS(JB)
              IF (KB(I)-KT < NL(JB)-1) IUT = I+1
            END DO
          ELSE
            DO I=US(JB)-1,DS(JB)+1
              IF (KB(I) > KBI(I)) THEN
                BNEW(KB(I),I)  = B(KB(I),I)                                                    ! SW 1/23/06   ! SW 3/2/05
                DX(KB(I),I) = 0.0
                KB(I)       = KB(I)-1
                ILAYER(I) = 1
                U(KB(I)+1,I)=0.0                                                               ! SW 1/23/06
                WRITE (WRN,'(2(A,I8),A,F0.3)') 'Raising bottom layer at segment ',I,' at iteration ',NIT,' at Julian day ',JDAY
                WARNING_OPEN = .TRUE.
              END IF
            END DO                   ! SW 1/23/06
            DO I=US(JB)-1,DS(JB)+1   ! SW 1/23/06
!              IF (KB(I)-KT < NL(JB)-1) IUT = I+1    ! SW 1/23/06
!                IF (I /= DS(JB)+1) KBMIN(I)   = MIN(KB(I),KB(I+1))                    ! SW 1/23/06                                 ! SW 3/2/05
                IF (I /= US(JB)-1) KBMIN(I-1) = MIN(KB(I-1),KB(I))                     ! SW 1/23/06
                iF(KBI(I) < KB(I))THEN
                BKT(I)=BH1(KT,I)/(H1(KT,I)-(EL(KBI(I)+1,I)-EL(KB(I)+1,I)))    ! SW 1/23/06
                DEPTHB(KTWB(JW),I)=(H1(KTWB(JW),I)-(EL(KBI(I)+1,I)-EL(KB(I)+1,I)))    ! SW 1/23/06
                DEPTHM(KTWB(JW),I)=(H1(KTWB(JW),I)-(EL(KBI(I)+1,I)-EL(KB(I)+1,I)))*0.5    ! SW 1/23/06
                IF(I<=DS(JB))THEN  ! SW 8/6/2018
                    AVHR(KT,I)=(H1(KT,I)-(EL(KBI(I)+1,I)-EL(KB(I)+1,I)))   +(H1(KT,I+1)-(EL(KBI(I)+1,I+1)-EL(KB(I)+1,I+1))&
                 -H1(KT,I)+(EL(KBI(I)+1,I)-EL(KB(I)+1,I)))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)   
                ELSE
                    AVHR(KT,I)=AVHR(KT,I-1)
                ENDIF
                ! SW 1/23/06
                end if
            ENDDO
        DO I=US(JB)-1,DS(JB)+1   ! SW 1/23/06
         DO K=KBMIN(I)+1,KB(I)
         U(K,I)=0.0
         END DO
        END DO    ! SW 11/9/07

        DO I=US(JB),DS(JB)     ! SW 11/9/07
         IF(ILAYER(I).EQ.1.AND.ILAYER(I+1).EQ.0)THEN
         BHRSUM=0.0
         Q(I)=0.0
          DO K=KT,KBMIN(I)
            IF (.NOT. INTERNAL_WEIR(K,I)) THEN
              BHRSUM = BHRSUM+BHR1(K,I)
              Q(I)   = Q(I)+U(K,I)*BHR1(K,I)
            END IF
          END DO
          DO K=KT,KBMIN(I)
            IF (INTERNAL_WEIR(K,I)) THEN
              U(K,I) = 0.0
            ELSE
              U(K,I) =  U(K,I)+(QC(I)-Q(I))/BHRSUM
            END IF
          END DO
          ELSEIF(ILAYER(I).EQ.1.AND.ILAYER(I-1).EQ.0)THEN
          BHRSUM=0.0
          Q(I-1)=0.0
          DO K=KT,KBMIN(I-1)
            IF (.NOT. INTERNAL_WEIR(K,I-1)) THEN
              BHRSUM = BHRSUM+BHR1(K,I-1)
              Q(I-1)   = Q(I-1)+U(K,I-1)*BHR1(K,I-1)
            END IF
          END DO
          DO K=KT,KBMIN(I-1)
            IF (INTERNAL_WEIR(K,I-1)) THEN
              U(K,I-1) = 0.0
            ELSE
              U(K,I-1) =  U(K,I-1)+(QC(I-1)-Q(I-1))/BHRSUM
            END IF
          END DO
          endif

            END DO
          END IF

!******** Segment addition

          IF (IUT /= IU) THEN
            IF (SNAPSHOT(JW)) WRITE (SNP(JW),'(/17X,2(A,I0))') ' Add segments ',IUT,' through ',IU-1
            WARNING_OPEN = .TRUE.
            WRITE (WRN,'(/17X,2(A,I0))') ' Add segments ',IUT,' through ',IU-1
            DO I=IUT-1,IU-1
              Z(I)         =  Z(IU)
              KTI(I)       =  KTI(IU)
              H1(KT+1,I)   =  H(KT+1,JW)
              AVH1(KT+1,I) = (H1(KT+1,I)+H1(KT+2,I))*0.5
              AVHR(KT+1,I) =  H1(KT+1,I)+(H1(KT+1,I+1)-H1(KT+1,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                  !SW 07/29/04
              IF (.NOT. TRAPEZOIDAL(JW)) THEN
                BH1(KT+1,I) =  B(KT+1,I)*H(KT+1,JW)
                H1(KT,I)    =  H(KT,JW)-Z(I)
                BI(KT:KB(I),I) =  B(KT:KB(I),I)                                                                        ! SW 4/18/07
                BI(KT,I)    =  B(KTI(I),I)
                AVH1(KT,I)  = (H1(KT,I)+H1(KT+1,I))*0.5
                AVHR(KT,I)  =  H1(KT,I)+(H1(KT,I+1)-H1(KT,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                       !SW 07/29/04
                BH1(KT,I)   =  BI(KT,I)*(EL(KT,I)-Z(I)*COSA(JB)-EL(KTI(I)+1,I))/COSA(JB)
                IF (KTI(I) >= KB(I)) BH1(KT,I) = B(KT,I)*H1(KT,I)
                DO K=KTI(I)+1,KT
                  BH1(KT,I) = BH1(KT,I)+BH(K,I)
                END DO
              ELSE
                CALL GRID_AREA1 (EL(KT,I)-Z(I),EL(KT+1,I),BH1(KT,I),BI(KT,I))                                                    !SW 08/03/04
                BH1(KT+1,I) =  0.25*H(KT+1,JW)*(BB(KT,I)+2.*B(KT+1,I)+BB(KT+1,I))
                H1(KT,I)    =  H(KT,JW)-Z(I)
                AVH1(KT,I)  = (H1(KT,I)+H1(KT+1,I))*0.5
                AVHR(KT,I)  =  H1(KT,I)+(H1(KT,I+1)-H1(KT,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                       !SW 07/29/04
              END IF
              BKT(I)      = BH1(KT,I)/H1(KT,I)
              DEPTHB(KT,I) = H1(KT,I)
              DEPTHM(KT,I) = H1(KT,I)*0.5
              DO K=KT+1,KB(I)
                DEPTHB(K,I) = DEPTHB(K-1,I)+ H1(K,I)
                DEPTHM(K,I) = DEPTHM(K-1,I)+(H1(K-1,I)+H1(K,I))*0.5
              END DO
            END DO
            DO I=IUT-1,IU-1
              BHR1(KT+1,I) = BH1(KT+1,I)+(BH1(KT+1,I+1)-BH1(KT+1,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                !SW 07/29/04
              BHR1(KT,I)   = BH1(KT,I)  +(BH1(KT,I+1)  -BH1(KT,I))  /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                !SW 07/29/04
                             IF(CONSTRICTION(KT,I))THEN    ! SW 6/26/2018 Valid for all K
                                        IF(BHR1(KT,I) > BCONSTRICTION(I)*H1(KT,I))BHR1(KT,I)= BCONSTRICTION(I)*H1(KT,I)
                                        IF(BHR1(KT+1,I) > BCONSTRICTION(I)*H1(KT+1,I))BHR1(KT+1,I)= BCONSTRICTION(I)*H1(KT+1,I)
                             ENDIF
            END DO
            DO I=IUT,IU-1
              !ICE(I)   = ICE(IU)   ! SW 9/29/15
              !ICETH(I) = ICETH(IU)
              WIND2(I) = WIND2(IU)
              IF (DYNAMIC_SHADE(I)) CALL SHADING
              DO K=KT,KBMIN(I)              !KB(I)   SW 12/18/2018
                U(K,I)  = U(K,IU)
                            IF(DXI(JW) >= 0.0)THEN
                                 DX(K,I) = DXI(JW)
                            ELSE
                                 DX(K,I) = ABS(U(K,I))*ABS(DXI(JW))*H(K,JW)    ! SW 8/2/2017
                            ENDIF
                SDKV(K,I)=SDK(JW)     ! SW 1/18/08
                IF (INTERNAL_WEIR(K,I)) THEN
                  DX(K,I) = 0.0
                  U(K,I)  = 0.0
                END IF
              ENDDO
              
            DO K=KT,KB(I)  ! SW 12/18/2018
                T1(K,I)           = T1(K,IU)
                T2(K,I)           = T1(K,IU)
                SU(K,I)           = U(K,IU)
                C1(K,I,CN(1:NAC)) = C1(K,IU,CN(1:NAC))
                C2(K,I,CN(1:NAC)) = C1(K,IU,CN(1:NAC))
                DO JE=1,NEP
                  EPD(K,I,JE) = 0.01
                  EPC(K,I,JE) = 0.01/H1(K,I)
                END DO
                CMBRT(CN(1:NAC),JB) = CMBRT(CN(1:NAC),JB)+C1(K,IU,CN(1:NAC))*DLX(I)*BH1(K,I)
                EBRI(JB)            = EBRI(JB)           +T1(K,IU)          *DLX(I)*BH1(K,I)
              END DO
              DO K=KT,KB(I)-1
                AZ(K,I)  = AZ(K,IU)
		        TKE(K,I,1) = TKE(K,IU,1) !sg 10/4/07
                TKE(K,I,2) = TKE(K,IU,2) !sg 10/4/07
                SAZ(K,I) = AZ(K,IU)
                IF (INTERNAL_WEIR(K,I)) THEN
                  AZ(K,I)  = 0.0
	              TKE(K,I,1) = 0.0 !sg 10/4/07
                  TKE(K,I,2) = 0.0 !sg 10/4/07
                  SAZ(K,I) = 0.0
                END IF
              END DO
            !END DO
!*********macrophytes
              DO M=1,NMC
                IF(MACROPHYTE_CALC(JW,M))THEN
                  JT=KTI(I)
                  JE=KB(I)
                  DO J=JT,JE
                    IF(J.LT.KT)THEN
                      COLB=EL(J+1,I)
                    ELSE
                      COLB=EL(KT+1,I)
                    END IF
                    !COLDEP=ELWS(I)-COLB
                    coldep=EL(KT,i)-Z(i)*COSA(JB)-colb  ! cb 3/7/16
                    !MACRC(J,KT,I,M)=MACWBCI(JW,M)
                    IF (ISO_macrophyte(JW,m))  macrc(j,kt,I,m) = macwbci(JW,m)     ! cb 3/7/16
                    IF (VERT_macrophyte(JW,m)) macrc(j,kt,I,m) = 0.1
                    IF (long_macrophyte(JW,m)) macrc(j,kt,I,m) = 0.1
                    MACRM(J,KT,I,M)=MACRC(J,KT,I,M)*CW(J,I)*COLDEP*DLX(I)
                    MACMBRT(JB,M) = MACMBRT(JB,M)+MACRM(J,KT,I,M)
                  END DO
                  DO K=KT+1,KB(I)
                    JT=K
                    JE=KB(I)
                    DO J=JT,JE
                      !MACRC(J,K,I,M)=MACWBCI(JW,M)                      
                      IF (ISO_macrophyte(JW,m))  macrc(j,k,I,m) = macwbci(JW,m)     ! cb 3/7/16
                      IF (VERT_macrophyte(JW,m)) macrc(j,k,I,m) = 0.1
                      IF (long_macrophyte(JW,m)) macrc(j,k,I,m) = 0.1
                      MACRM(J,K,I,M)=MACRC(J,K,I,M)*CW(J,I)*H2(K,I)*DLX(I)
                      MACMBRT(JB,M) = MACMBRT(JB,M)+MACRM(J,K,I,M)
                    END DO
                  END DO
                END IF
              END DO
            END DO  ! cb 3/7/16  moved enddo to include macrophytes
            U(KB(IUT):KB(IU),IU-1)  = 0.0
            SU(KB(IUT):KB(IU),IU-1) = 0.0
            ADX(KB(IUT):KB(IU),IU)  = 0.0
            IU                      = IUT
            CUS(JB)                 = IU
            IF (UH_EXTERNAL(JB)) KB(IU-1) = KB(IU)
            IF (UH_INTERNAL(JB)) THEN
              IF (JBUH(JB) >= BS(JW) .AND. JBUH(JB) <= BE(JW)) THEN
                KB(IU-1) = MIN(KB(UHS(JB)),KB(IU))
              ELSE
                DO KKB=KT,KMX
                  IF (EL(KKB,IU) <= EL(KB(UHS(JB)),UHS(JB))) EXIT
                END DO
                KB(IU-1) = MIN(KKB,KB(IU))
              END IF
            END IF
            IF (UP_HEAD(JB)) THEN
              AZ(KT:KB(IU-1)-1,IU-1)  = AZMIN
              TKE(KT:KB(IU-1)-1,IU-1,1) = 1.25E-7     !SG 10/4/07
              TKE(KT:KB(IU-1)-1,IU-1,2) = 1.0E-9      !SG 10/4/07
              SAZ(KT:KB(IU-1)-1,IU-1) = AZMIN
            END IF
          END IF
          IF (CONSTITUENTS) THEN
            CALL TEMPERATURE_RATES
            CALL KINETIC_RATES
          END IF

!******** Total active cells and single layers

          DO I=IU,ID
            NTAC         = NTAC+1
            ONE_LAYER(I) = KTWB(JW) == KB(I)
          END DO
          NTACMX = MAX(NTAC,NTACMX)
        END DO
        CALL INTERPOLATION_MULTIPLIERS

!****** Additional layers

        ZMIN(JW) = -1000.0
        DO JB=BS(JW),BE(JW)
        IF(BR_INACTIVE(JB))CYCLE    ! SW 6/12/2017
          DO I=CUS(JB),DS(JB)
            ZMIN(JW) = MAX(ZMIN(JW),Z(I))
          END DO
        END DO
        ADD_LAYER = ZMIN(JW) < -0.80*H(KT-1,JW) .AND. KT /= 2
      END DO

!**** Subtract layers

      DO WHILE (SUB_LAYER)
        IF (SNAPSHOT(JW)) WRITE (SNP(JW),'(/1X,13("*"),1X,A,I0,A,F0.3,A,I0,1X,A,I0,1x,13("*"))') 'Subtract layer ',&
        KT,' at Julian day = ', JDAY,' NIT = ',NIT,' IZMIN =',IZMIN(JW)      ! SW 1/23/06   
        WRITE (WRN,'(/1X,13("*"),1X,A,I0,A,F0.3,A,I0,1X,A,I0,1x,13("*"))') 'Subtract layer ',KT,' at Julian day = ', JDAY,' NIT = ',NIT,' IZMIN =',IZMIN(JW)   

!****** Variable initialization

        KTWB(JW) = KTWB(JW)+1
        KT       = KTWB(JW)
        ILAYER=0       ! SW 1/23/06  11/7/07
        
! RECOMPUTE INTERNAL WEIR FOR FLOATING WEIR      
    IF (WEIR_CALC) THEN   !  SW 3/16/18
    DO JWR=1,NIW
     IF(IWR(JWR) >= US(BS(JW)) .AND. IWR(JWR) <= DS(BE(JW)))THEN
        IF (EKTWR(JWR) == 0.0) THEN  
            KTWR(JWR)=KTWB(JW)
        ELSE  
          KTWR(JWR) = INT(EKTWR(JWR))  
        END IF 
        IF (EKBWR(JWR) <= 0.0) THEN  
            DO K=KTWR(JWR),KB(IWR(JWR))  
            IF (DEPTHB(K,IWR(JWR)) > ABS(EKBWR(JWR))) THEN
                KBWR(JWR)=K
                EXIT  
                ENDIF
            END DO   
        ELSE  
          KBWR(JWR) = INT(EKBWR(JWR))  
        END IF  
        
      DO K=2,KMX-1
        IF ((K >= KTWR(JWR) .AND. K <= KBWR(JWR))) INTERNAL_WEIR(K,IWR(JWR)) = .TRUE.
      END DO
     ENDIF     
    END DO    
  END IF
        
        
        
        DO JB=BS(JW),BE(JW)
        IF(BR_INACTIVE(JB))CYCLE    ! SW 6/12/2017
          IU = CUS(JB)
          ID = DS(JB)
          IF (CONSTITUENTS) DO1(KT-1,IU-1:ID+1) = 0.0
          DO I=IU-1,ID+1
            Z(I)         =  Z(I)-H(KT-1,JW)
            H1(KT-1,I)   =  H(KT-1,JW)
            H1(KT,I)     =  H(KT,JW)-Z(I)
            BI(KT,I)     =  B(KTI(I),I)
            BI(KT-1,I)   =  B(KT-1,I)
            AVH1(KT-1,I) = (H(KT-1,JW)+H(KT,JW))*0.5
            AVH1(KT,I)   = (H1(KT,I)+H1(KT+1,I))*0.5
            IF(.NOT. TRAPEZOIDAL(JW))THEN
              IF(KB(I) >= KT)THEN                            ! SW 1/23/06
              BH1(KT,I)   =  BH1(KT-1,I)+BH1(KT,I)           ! SW 1/23/06
              BH1(KT-1,I) =  B(KT-1,I)*H(KT-1,JW)
              ELSE  ! SW 1/23/06
              BH1(KT,I)   =  BH1(KT-1,I)         ! SW 1/23/06
              END IF  ! SW 1/23/06
            ELSE
              CALL GRID_AREA1 (EL(KT,I)-Z(I),EL(KT+1,I),BH1(KT,I),DUMMY)                                                         !SW 08/03/04
              BH1(KT-1,I) = 0.25*H1(KT-1,JW)*(BB(KT-2,I)+2.*B(KT-1,I)+BB(KT-1,I))
            ENDIF
            VOL(KT,I)   =  BH1(KT,I)*DLX(I)
            VOL(KT-1,I) =  BH1(KT-1,I)*DLX(I)
            BKT(I)      =  BH1(KT,I)/H1(KT,I)
            if(kb(i) >= kt)then                             ! SW 1/23/06
              U(KT,I)  = (U(KT-1,I)*BHR1(KT-1,I)+U(KT,I)*BHR(KT,I))/(BHR1(KT-1,I)+BHR(KT,I))
              T1(KT,I) = (T1(KT-1,I)*(BH1(KT,I)-BH(KT,I))+T1(KT,I)*BH(KT,I))/BH1(KT,I)
            ELSE
!              EBRI(JB) = EBRI(JB)-T1(KT,I)*VOL(KT,I)   1/23/06
              u(kt,i) = u(kt-1,i)    ! SW 1/23/06
              t1(kt,i) = t1(kt-1,i)  ! SW 1/23/06
            END IF
            IF(KB(I) >= KT)THEN       ! SW 1/23/06
            C1(KT,I,CN(1:NAC))             = (C1(KT-1,I,CN(1:NAC))  *(BH1(KT,I)-BH(KT,I))+C1(KT,I,CN(1:NAC))  *BH(KT,I))/BH1(KT,I)
            CSSK(KT,I,CN(1:NAC))           = (CSSK(KT-1,I,CN(1:NAC))*(BH1(KT,I)-BH(KT,I))+CSSK(KT,I,CN(1:NAC))*BH(KT,I))/BH1(KT,I)
            ELSE
            C1(KT,I,CN(1:NAC))             = C1(KT-1,I,CN(1:NAC))    ! SW 1/23/06
            CSSK(KT,I,CN(1:NAC))           = CSSK(KT-1,I,CN(1:NAC))   ! SW 1/23/06
            ENDIF   ! SW 1/23/06
            CSSB(KT,I,CN(1:NAC))           =  CSSB(KT-1,I,CN(1:NAC))+CSSB(KT,I,CN(1:NAC))
            !KF(KT,I,KFCN(1:NAF(JW),JW))    =  KF(KT-1,I,KFCN(1:NAF(JW),JW))
            !KFS(KT,I,KFCN(1:NAF(JW),JW))   =  KFS(KT-1,I,KFCN(1:NAF(JW),JW))
            KF(KT,I,KFCN(1:NAF(JW),JW))    =  (KF(KT-1,I,KFCN(1:NAF(JW),JW))*VOL(KT-1,I)+  KF(KT,I,KFCN(1:NAF(JW),JW))*VOL(KT,I))/(VOL(KT-1,I)+VOL(KT,I))   ! SW Fix suggested by Taylor Adams Hydros 25Feb2021  ! KF is in units of g/m3/s
            KFS(KT,I,KFCN(1:NAF(JW),JW))   =  KFS(KT-1,I,KFCN(1:NAF(JW),JW))+KFS(KT,I,KFCN(1:NAF(JW),JW))          ! SW Fix suggested by Taylor Adams Hydros 25Feb2021 ! KFS is in units of g KFS=KF*VOL*DT          C1(KT-1,I,CN(1:NAC))           =  0.0
            C1(KT-1,I,CN(1:NAC))           =  0.0
            C2(KT-1,I,CN(1:NAC))           =  0.0
            CSSB(KT-1,I,CN(1:NAC))         =  0.0
            CSSK(KT-1,I,CN(1:NAC))         =  0.0
            KF(KT-1,I,KFCN(1:NAF(JW),JW))  =  0.0
            KFS(KT-1,I,KFCN(1:NAF(JW),JW)) =  0.0
            DO JE=1,NEP
              IF(KT <= KBI(I))THEN    ! CB 4/28/06
              EPM(KT,I,JE)   = EPM(KT-1,I,JE)+EPM(KT,I,JE)
              EPD(KT,I,JE)   = EPM(KT,I,JE)/((BI(KT,I)-BI(KT+1,I)+2.0*H1(KT,I))*DLX(I))
              EPC(KT,I,JE)   = EPM(KT,I,JE)/VOL(KT,I)
              EPM(KT-1,I,JE) = 0.0
              EPD(KT-1,I,JE) = 0.0
              EPC(KT-1,I,JE) = 0.0
              ELSE                   ! SW 5/15/06
              EPM(KT,I,JE)   = EPM(KT-1,I,JE)
              EPD(KT,I,JE)   = EPD(KT-1,I,JE)
              EPC(KT,I,JE)   = EPC(KT-1,I,JE)
              EPM(KT-1,I,JE) = 0.0
              EPD(KT-1,I,JE) = 0.0
              EPC(KT-1,I,JE) = 0.0
              ENDIF                   ! CB 4/28/06
            END DO
          END DO
          DO I=IU,ID
            DO M=1,NMC
              IF(MACROPHYTE_CALC(JW,M))THEN
                IF(KTICOL(I))THEN
                  JT=KTI(I)
                ELSE
                  JT=KTI(I)+1
                END IF
                JE=KB(I)
                DO J=JT,JE
                  IF(J.LT.KT)THEN
                    COLB=EL(J+1,I)
                  ELSE
                    COLB=EL(KT+1,I)
                  END IF
                  !COLDEP=ELWS(I)-COLB
                   coldep=EL(KT,i)-Z(i)*COSA(JB)-colb  ! cb 3/7/16
                  IF(J.LT.KT)THEN
                    MACRM(J,KT,I,M)=MACRM(J,KT-1,I,M)
                  ELSE
                    MACRM(J,KT,I,M)=MACRM(J,KT-1,I,M)+MACRM(J,KT,I,M)
                  END IF
                  IF(CW(J,I).GT.0.0)THEN
                    MACRC(J,KT,I,M)=MACRM(J,KT,I,M)/(CW(J,I)*COLDEP*DLX(I))
                  ELSE
                    MACRC(J,KT,I,M)=0.0
                  END IF
                  MACRM(J,KT-1,I,M)=0.0
                  MACRC(J,KT-1,I,M)=0.0

                END DO
              END IF
            END DO
            JT=KTI(I)
            JE=KB(I)
            DO J=JT,JE
              MACT(J,KT,I)=0.0
              MACT(J,KT-1,I)=0.0
            END DO
            DO M=1,NMC
              IF(MACROPHYTE_CALC(JW,M))THEN
                 DO J=JT,JE
                   MACT(J,KT,I)=MACRC(J,KT,I,M)+MACT(J,KT,I)
                 END DO
              END IF
            END DO
            DO M=1,NMC
              TMAC=0.0
              IF(MACROPHYTE_CALC(JW,M))THEN
                JT=KTI(I)
                JE=KB(I)
                DO J=JT,JE
                  TMAC=TMAC+MACRM(J,KT,I,M)
                END DO
              END IF
              MAC(KT,I,M)=TMAC/(BH1(KT,I)*DLX(I))
            END DO
            DO M=1,NMC
              IF(MACROPHYTE_CALC(JW,M))THEN
                MAC(KT-1,I,M)=0.0
              END IF
            END DO
          END DO

          DO I=IU-1,ID
            AVHR(KT-1,I) =  H1(KT-1,I)+( H1(KT-1,I+1)- H1(KT-1,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                  !SW 07/29/04
            AVHR(KT,I)   =  H1(KT,I)  +( H1(KT,I+1)  - H1(KT,I))  /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                  !SW 07/29/04
            BHR1(KT,I)   = BH1(KT,I)  +(BH1(KT,I+1)  -BH1(KT,I))  /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                  !SW 07/29/04
            BHR1(KT-1,I) = BH1(KT-1,I)+(BH1(KT-1,I+1)-BH1(KT-1,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                  !SW 07/29/04
                             IF(CONSTRICTION(KT,I))THEN    ! SW 6/26/2018 Valid for all K
                                        IF(BHR1(KT,I) > BCONSTRICTION(I)*H1(KT,I))BHR1(KT,I)= BCONSTRICTION(I)*H1(KT,I)
                                        IF(BHR1(KT-1,I) > BCONSTRICTION(I)*H1(KT-1,I))BHR1(KT-1,I)= BCONSTRICTION(I)*H1(KT-1,I)
                             ENDIF            
          END DO
          U(KT-1,IU-1:ID+1)     = 0.0
          W(KT-1,IU-1:ID+1)     = 0.0
          SW(KT-1,IU-1:ID+1)    = 0.0                                                                                  !TC 3/9/05
          P(KT-1,IU-1:ID+1)     = 0.0
          AZ(KT-1,IU-1:ID+1)    = 0.0
          TKE(KT-1,IU-1:ID+1,1) = 0.0   !sg 10/4/07
          TKE(KT-1,IU-1:ID+1,2) = 0.0   !sg 10/4/07
          DZ(KT-1,IU-1:ID+1)    = 0.0
          ADMZ(KT-1,IU-1:ID+1)  = 0.0
          ADZ(KT-1,IU-1:ID+1)   = 0.0
          DECAY(KT-1,IU-1:ID+1) = 0.0
          IF (UP_HEAD(JB)) THEN
            QUH1(KT,JB)             = QUH1(KT,JB)+QUH1(KT-1,JB)
            TSSUH1(KT,JB)           = TSSUH1(KT-1,JB)          +TSSUH1(KT,JB)
            CSSUH1(KT,CN(1:NAC),JB) = CSSUH1(KT-1,CN(1:NAC),JB)+CSSUH1(KT,CN(1:NAC),JB)
          END IF
          IF (DN_HEAD(JB)) THEN
            QDH1(KT,JB)             = QDH1(KT,JB)+QDH1(KT-1,JB)
            TSSDH1(KT,JB)           = TSSDH1(KT-1,JB)          +TSSDH1(KT,JB)
            CSSDH1(KT,CN(1:NAC),JB) = CSSDH1(KT-1,CN(1:NAC),JB)+CSSDH1(KT,CN(1:NAC),JB)
          END IF

!******** Upstream active segment

          IUT = US(JB)
          IF (SLOPE(JB) /= 0.0) THEN
            DO I=US(JB)-1,DS(JB)+1
              IF (KB(I) < KT ) THEN                                                                                  ! SR 10/17/05
                KB(I)                 = KT
                Bnew(KB(I),I)         = 0.000001   ! sw 1/23/06
                IF(DXI(JW) >= 0.0)THEN
                     DX(KB(I),I) = DXI(JW)
                ELSE
                     DX(KB(I),I) = ABS(U(KB(I),I))*ABS(DXI(JW))*H(K,JW)    ! SW 8/2/2017
                ENDIF
                ilayer(i)=1
                T1(KB(I),I)           = T1(KT,I)                   !    SW 5/15/06    T1(KB(I)-1,I)
                C1(KB(I),I,CN(1:NAC)) = C1(KT,I,CN(1:NAC))         !    SW 5/15/06    C1(KB(I)-1,I,CN(1:NAC))
                
                IF (SEDIMENT_CALC(JW))THEN      ! SW 5/26/2022
                SED(KB(I),I)=SED(KT-1,I);SEDC(KB(I),I)=SEDC(KT-1,I);SEDN(KB(I),I)=SEDN(KT-1,I);SEDP(KB(I),I)=SEDP(KT-1,I)
                ENDIF              
                
                WRITE (WRN,'(2(A,I8),A,F0.3,A,F0.3)') 'Lowering bottom segment ',I,' at iteration ',NIT,' at Julian day ',&
                                                       JDAY,' Z(I)=',Z(I)
                WARNING_OPEN = .TRUE.
              END IF
            ENDDO                    ! SW 1/23/06
            DO I=US(JB)-1,DS(JB)+1   ! SW 1/23/06
!               IF (I /= DS(JB)+1) KBMIN(I)   = MIN(KB(I),KB(I+1))  ! SW 1/23/06
                IF (I /= US(JB)-1) KBMIN(I-1) = MIN(KB(I-1),KB(I))  ! SW 1/23/06
                IF(KBI(I) < KB(I))THEN
                BKT(I)=BH1(KT,I)/(H1(KT,I)-(EL(KBI(I)+1,I)-EL(KB(I)+1,I)))    ! SW 1/23/06
                DEPTHB(KTWB(JW),I)=(H1(KTWB(JW),I)-(EL(KBI(I)+1,I)-EL(KB(I)+1,I)))    ! SW 1/23/06
                DEPTHM(KTWB(JW),I)=(H1(KTWB(JW),I)-(EL(KBI(I)+1,I)-EL(KB(I)+1,I)))*0.5    ! SW 1/23/06
                AVHR(KT,I)=(H1(KT,I)-(EL(KBI(I)+1,I)-EL(KB(I)+1,I)))   +(H1(KT,I+1)-(EL(KBI(I)+1,I+1)-EL(KB(I)+1,I+1))&
                                           -H1(KT,I)+(EL(KBI(I)+1,I)-EL(KB(I)+1,I)))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                                                                           ! SW 1/23/06
                END IF
            ENDDO

        DO I=US(JB),DS(JB)   ! SW 1/23/06   11/13/07   US(JB)-1,DS(JB)+1

         IF(ILAYER(I).EQ.1.AND.ILAYER(I+1).EQ.0)THEN  ! SW 1/23/06
         BHRSUM=0.0
         Q(I)=0.0
          DO K=KT,KBMIN(I)
            IF (.NOT. INTERNAL_WEIR(K,I)) THEN
              BHRSUM = BHRSUM+BHR1(K,I)
              Q(I)   = Q(I)+U(K,I)*BHR1(K,I)
            END IF
          END DO
          DO K=KT,KBMIN(I)
            IF (INTERNAL_WEIR(K,I)) THEN
              U(K,I) = 0.0
            ELSE
              U(K,I) =  U(K,I)+(QC(I)-Q(I))/BHRSUM
            END IF
          END DO
          ELSEIF(ILAYER(I).EQ.1.AND.ILAYER(I-1).EQ.0)THEN
          BHRSUM=0.0
          Q(I-1)=0.0
          DO K=KT,KBMIN(I-1)
            IF (.NOT. INTERNAL_WEIR(K,I-1)) THEN
              BHRSUM = BHRSUM+BHR1(K,I-1)
              Q(I-1)   = Q(I-1)+U(K,I-1)*BHR1(K,I-1)
            END IF
          END DO
          DO K=KT,KBMIN(I-1)
            IF (INTERNAL_WEIR(K,I-1)) THEN
              U(K,I-1) = 0.0
            ELSE
              U(K,I-1) =  U(K,I-1)+(QC(I-1)-Q(I-1))/BHRSUM
            END IF
          END DO
          ENDIF    ! SW 1/23/06

            END DO
          END IF
          DO I=US(JB),DS(JB)
            IF (KB(I)-KT < NL(JB)-1) IUT = I+1
            ONE_LAYER(I) = KTWB(JW) == KB(I)
          END DO
          IF (IUT > DS(JB)) THEN
            IF(JB==1)THEN
            WRITE (W2ERR,'(A,I0/A,F0.2,2(A,I0))') 'Fatal error - insufficient segments in branch ',JB,'Julian day = ',JDAY,      &    ! SEE NEW CODE BELOW
                                                  ' at iteration ',NIT,' with water surface layer = ',KT
            WRITE (W2ERR,'(2(A,I0))')             'Minimum water surface located at segment ',IZMIN(JW),' with bottom layer at ',&
                                                   KB(IZMIN(JW))
            TEXT = 'Runtime error - see w2.err'
            ERROR_OPEN = .TRUE.
            RETURN
            ELSE
            BR_INACTIVE(JB)=.TRUE.  ! SW 6/12/2017
            IF (SNAPSHOT(JW)) WRITE (SNP(JW),'(/1X,13("*"),1X,A,I0,A,F0.3,A,I0,1X,A,I0,13("*"))') '   Branch Inactive: ',jb,&
                                                        ' at Julian day = ',JDAY,'   NIT = ',NIT 
            WARNING_OPEN = .TRUE.
            WRITE (WRN,'(/1X,13("*"),1X,A,I0,A,F0.3,A,I0,1X,A,I0,13("*"))') '   Branch Inactive: ',jb,' at Julian day = ',JDAY,'   NIT = ',NIT    ! SW 11/16/2018
            DO I=IU,DS(JB)   ! SW 12/17/2018
              DO K=KT-1,KB(I)       ! SW 12/18/2018  KT,KB(I)
                EBRI(JB)            = EBRI(JB)-T1(K,I)*BH1(K,I)*DLX(I)    ! VOL(K,I)   SW 12/18/2018
                CMBRT(CN(1:NAC),JB) = CMBRT(CN(1:NAC),JB)-C1(K,I,CN(1:NAC))*BH1(K,I)*DLX(I)   !VOL(K,I)              !+(CSSB(K,I,CN(1:NAC))+CSSK(K,I,CN(1:NAC))*VOL(K,I))*DLT
              END DO
            END DO
            CYCLE
            ENDIF
            
       ! Go to 230
          END IF

!******** Segment subtraction

          IF (IUT /= IU) THEN
            IF (SNAPSHOT(JW)) WRITE (SNP(JW),'(/17X,A,I0,A,I0)') ' Subtract segments ',IU,' through ',IUT-1
            WARNING_OPEN = .TRUE.
            WRITE(WRN,'(/17X,A,I0,A,I0)') ' Subtract segments ',IU,' through ',IUT-1
            DO I=IU,IUT-1
              DO K=KT-1,KB(I)     ! SW 12/18/2018 KT,KB(I)
                EBRI(JB)            = EBRI(JB)-T1(K,I)*BH1(K,I)*DLX(I)      !*VOL(K,I)    SW 12/18/2018
                CMBRT(CN(1:NAC),JB) = CMBRT(CN(1:NAC),JB)-C1(K,I,CN(1:NAC))*BH1(K,I)*DLX(I)    !*VOL(K,I)        !+(CSSB(K,I,CN(1:NAC))+CSSK(K,I,CN(1:NAC))*VOL(K,I))*DLT
              END DO
            END DO

            DO I=IU,IUT-1
              DO M=1,NMC
                IF(MACROPHYTE_CALC(JW,M))THEN
                  JT=KTI(I)
                  JE=KB(I)
                  DO J=JT,JE
                    IF(J.LT.KT)THEN
                      COLB=EL(J+1,I)
                    ELSE
                      COLB=EL(KT+1,I)
                    END IF
                    !COLDEP=ELWS(I)-COLB
                     coldep=EL(KT,i)-Z(i)*COSA(JB)-colb  ! cb 3/7/16
!                    MACMBRT(JB,M) = MACMBRT(JB,M)-MACRM(J,KT,I,M)+(MACSS(J,KT,I,M)*COLDEP*CW(J,I)*DLX(I))*DLT
                    MACMBRT(JB,M) = MACMBRT(JB,M)-MACRM(J,KT,I,M)
                  END DO
                  DO K=KT+1,KB(I)
                    JT=K
                    JE=KB(I)
                    DO J=JT,JE
!                      MACMBRT(JB,M) = MACMBRT(JB,M)-MACRM(J,K,I,M)+(MACSS(J,K,I,M)*H2(K,I)*CW(J,I)*DMX(I))*DLT
                      MACMBRT(JB,M) = MACMBRT(JB,M)-MACRM(J,K,I,M)
                    END DO
                  END DO
                END IF
              END DO
            END DO
            
            F(IU-1:IUT-1)     =  0.0
            Z(IU-1:IUT-1)     =  0.0
            !ICETH(IU-1:IUT-1) =  0.0    ! SW 9/29/15
            BHRHO(IU-1:IUT-1) =  0.0
            !ICE(IU-1:IUT-1)   = .FALSE.  ! SW 9/29/15
            DO K=KT,KB(IUT)
              ADX(K,IU-1:IUT-1)            = 0.0
              DX(K,IU-1:IUT-1)             = 0.0
              AZ(K,IU-1:IUT-1)             = 0.0
              TKE(K,IU-1:IUT-1,1)          = 0.0 !SG  10/4/07
              TKE(K,IU-1:IUT-1,2)          = 0.0 !SG  10/4/07
              SAZ(K,IU-1:IUT-1)            = 0.0
              U(K,IU-1:IUT-1)              = 0.0
              SU(K,IU-1:IUT-1)             = 0.0
              T1(K,IU-1:IUT-1)             = 0.0
              TSS(K,IU-1:IUT-1)            = 0.0
              QSS(K,IU-1:IUT-1)            = 0.0
              C1(K,IU-1:IUT-1,CN(1:NAC))   = 0.0
              C2(K,IU-1:IUT-1,CN(1:NAC))   = 0.0
              C1S(K,IU-1:IUT-1,CN(1:NAC))  = 0.0
              CSSB(K,IU-1:IUT-1,CN(1:NAC)) = 0.0
              CSSK(K,IU-1:IUT-1,CN(1:NAC)) = 0.0
            END DO
            
            DO M=1,NMC
              IF (MACROPHYTE_CALC(JW,M)) THEN
                MAC(K,I,M)=0.0
                MACT(J,K,I)=0.0
              END IF
            END DO
            JT=KTI(I)
            JE=KB(I)
            DO J=JT,JE
              DO M=1,NMC
                IF(MACROPHYTE_CALC(JW,M))THEN
                  MACRC(J,K,I,M) = 0.0
                END IF
              END DO
            END DO

            IU           =  IUT
            CUS(JB)      =  IU
            Z(IU-1)      = (EL(KT,IU-1)-(EL(KT,IU)-Z(IU)*COSA(JB)))/COSA(JB)
            SZ(IU-1)     =  Z(IU)
            KTI(IU-1)    =  KTI(IU)
            IF (.NOT. TRAPEZOIDAL(JW)) THEN
              BI(KT,IU-1)  = B(KTI(IU-1),I)
              H1(KT,IU-1)  = H(KT,JW)-Z(IU-1)
              BH1(KT,IU-1) = Bnew(KTI(IU-1),IU-1)*(EL(KT,IU-1)-EL(KTI(IU-1)+1,IU-1)-Z(IU-1)*COSA(JB))/COSA(JB)   ! sw 1/23/06  Bnew(KTI(IU-1),IU-1)*(EL(KT,IU-1)-EL(KTI(IU-1)+1,IU-1)-Z(IU-1)*COSA(JB))/COSA(JB)     ! SR 10/17/05
              IF (KT >= KB(IU-1)) BH1(KT,IU-1) = Bnew(KT,IU-1)*H1(KT,IU-1)   ! sw 1/23/06
              DO K=KTI(IU-1)+1,KT
                BH1(KT,IU-1) = BH1(KT,IU-1)+BH1(K,IU-1)
              END DO
            ELSE
              CALL GRID_AREA1 (EL(KT,I)-Z(I),EL(KT+1,IU-1),BH1(KT,IU-1),BI(KT,IU-1))                                             !SW 08/03/04
              BH1(KT,I) = 0.25*H(KT,JW)*(BB(KT-1,I)+2.*B(KT,I)+BB(KT,I))
              H1(KT,I)  = H(KT,JW)-Z(I)
            END IF
            BKT(IU-1)     =  BH1(KT,IU-1)/H1(KT,IU-1)
            BHR1(KT,IU-1) =  BH1(KT,IU-1)+(BH1(KT,IU)-BH1(KT,IU-1))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                 !SW 07/29/04
                            IF(CONSTRICTION(KT,IU-1))THEN    ! SW 6/26/2018 Valid for all K
                                        IF(BHR1(KT,IU-1) > BCONSTRICTION(IU-1)*H1(KT,IU-1))BHR1(KT,IU-1)= BCONSTRICTION(IU-1)*H1(KT,IU-1)
                             ENDIF
            IF (UH_EXTERNAL(JB)) KB(IU-1) = KB(IU)
            IF (UH_INTERNAL(JB)) THEN
              IF (JBUH(JB) >= BS(JW) .AND. JBUH(JB) <= BE(JW)) THEN
                KB(IU-1) = MIN(KB(UHS(JB)),KB(IU))
              ELSE
                DO KKB = KT, KMX
                  IF (EL(KKB,IU) <= EL(KB(UHS(JB)),UHS(JB))) EXIT
                END DO
                KB(IU-1) = MIN(KKB,KB(IU))
              END IF
            END IF
          END IF
          IF (CONSTITUENTS) THEN      ! SW 5/15/06
            CALL TEMPERATURE_RATES
            CALL KINETIC_RATES
          END IF


!******** Total active cells

          DO I=IU,ID
            NTAC = NTAC-1
          END DO
          NTACMN = MIN(NTAC,NTACMN)
        END DO
        CALL INTERPOLATION_MULTIPLIERS

!****** Additional layer subtractions

        ZMIN(JW) = -1000.0
        DO JB=BS(JW),BE(JW)
        IF(BR_INACTIVE(JB))CYCLE    ! SW 6/12/2017
          DO I=CUS(JB),DS(JB)
            ZMIN(JW) = MAX(ZMIN(JW),Z(I))
          END DO
        END DO
        SUB_LAYER = ZMIN(JW) > 0.60*H(KT,JW) .AND. KT < KTMAX                                                         ! SR 10/17/05
      END DO
    END DO

!** Temporary downstream head segment

    DO JB=1,NBR
    IF(BR_INACTIVE(JB))CYCLE    ! SW 6/12/2017
      IF (DHS(JB) > 0) THEN
        DO JJB=1,NBR
          IF (DHS(JB) >= US(JJB) .AND. DHS(JB) <= DS(JJB)) EXIT
        END DO
        IF (CUS(JJB) > DHS(JB)) CDHS(JB) = CUS(JJB)
      END IF
    END DO
RETURN
END SUBROUTINE LAYERADDSUB
