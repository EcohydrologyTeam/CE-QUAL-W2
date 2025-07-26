!************************************************************************
!**               S U B R O U T I N E    porosity                      **
!************************************************************************

      SUBROUTINE POROSITY

      USE GEOMC;USE GLOBAL;USE MACROPHYTEC;USE POROSITYC;USE LOGICC
      USE SCREENC
      IMPLICIT NONE

      INTEGER :: K,M,KUP,IEXIT,K1,KDN,IUT
      REAL    :: B11,VSTOT,ELR,ELL,ELL1,ELL2,ELR2,EL1,EL2,ERR1,ERR2

    IF(NIT.EQ.0)THEN

 1040 FORMAT((8X,I8,3F8.0))

      DO JW=1,NWB
        KT = KTWB(JW)
        DO JB=BS(JW),BE(JW)
          IU   = CUS(JB)
          ID   = DS(JB)
          DO I=IU,ID
            DO K=2,KB(I)
              VOLI(K,I) = BH(K,I)*DLX(I)
            END DO
            VOLI(KT,I)    = BH2(KT,I)*DLX(I)
          END DO
        END DO
      END DO

    END IF

    DO JB=1,NBR
      COSA(JB)=COS(ALPHA(JB))
    END DO

!C  CALCULATING # OF MACROPHYTE STEMS IN EACH CELL

    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        IU   = CUS(JB)
        ID   = DS(JB)
        DO I=IU,ID
!          HKTI  = H(KT,JW)-Z(I)    REPLACED BY H1(KT,I)
          IF(KT.EQ.KTI(I))THEN
            VOLKTI(I)=H1(KT,I)*BIC(KT,I)*DLX(I)
          ELSE
            VOLKTI(I) = BIC(KTI(I),I)*(EL(KT,I)-EL(KTI(I)+1,I)-Z(I)*COSA(JB))/COSA(JB)*DLX(I)
          END IF
          DO K=KTI(I)+1,KT
            VOLKTI(I) = VOLKTI(I)+VOLI(K,I)
          END DO

          DO M=1,NMC
            VSTEMKT(I,M)=(MAC(KT,I,M)*VOLKTI(I))/DWV(M)    !CB 6/29/06
          END DO

          DO K=KT+1,KB(I)
            DO M=1,NMC
              VSTEM(K,I,M)=(MAC(K,I,M)*VOLI(K,I))/DWV(M)   !CB 6/29/06
            END DO
          END DO
        END DO
      END DO
    END DO

    POR=1.0
    DO JW=1,NWB
      KT = KTWB(JW)
      DO JB=BS(JW),BE(JW)
        IU = CUS(JB)
        ID = DS(JB)

        DO I=IU,ID
          DO K=KT,KB(I)
            IF(K.EQ.KT)THEN
              VSTOT=0.0
              DO M=1,NMC
                  VSTOT=VSTOT+VSTEMKT(I,M)
              END DO
              POR(KT,I)=(VOLKTI(I)-VSTOT)/VOLKTI(I)
            ELSE
              VSTOT=0.0
              DO M=1,NMC
                VSTOT=VSTOT+VSTEM(K,I,M)
              END DO
              POR(K,I)=(VOLI(K,I)-VSTOT)/VOLI(K,I)
            END IF
          END DO
        END DO

        DO I=IU,ID
          DO K=KTI(I),KB(I)
            IF(K.LE.KT)THEN
              B(K,I)=POR(KT,I)*BIC(K,I)
            ELSE
              B(K,I)=POR(K,I)*BIC(K,I)
            END IF

          END DO
        END DO

      END DO
    END DO



! BOUNDARY WIDTHS

  DO JW=1,NWB
    KT = KTWB(JW)
    DO JB=BS(JW),BE(JW)
      IU = US(JB)
      ID = DS(JB)
      DO I=IU-1,ID+1
        B(1,I) = B(2,I)
        DO K=KB(I)+1,KMX
          B(K,I) = B(KB(I),I)
        END DO
      END DO
    END DO
  END DO
  DO JW=1,NWB
    KT = KTWB(JW)
    DO JB=BS(JW),BE(JW)
      IU    = US(JB)
      ID    = DS(JB)
      IEXIT = 0
      DO K=1,KMX-1
        B(K,IU-1) = B(K,IU)
        IF (UH_INTERNAL(JB) .OR. HEAD_FLOW(JB)) THEN
          IF (JBUH(JB) >= BS(JW) .AND. JBUH(JB) <= BE(JW)) THEN
            B(K,IU-1) = B(K,UHS(JB))
          ELSE
            ELR = EL(K,IU)+SINA(JB)*DLX(IU)*0.5
            ELL = EL(2,UHS(JB))-SINA(JBUH(JB))*DLX(UHS(JB))*0.5
            IF (ELR >= ELL) THEN
              B(K,IU-1) = B(2,UHS(JB))
            ELSE
              DO KUP=2,KMX-1
                ELL1 = EL(KUP,UHS(JB))-SINA(JBUH(JB))*DLX(UHS(JB))*0.5
                ELL2 = EL(KUP+1,UHS(JB))-SINA(JBUH(JB))*DLX(UHS(JB))*0.5
                IF (ELL1 > ELR .AND. ELL2 <= ELR) THEN
                  IF (KUP > KB(UHS(JB))) THEN
                    KB(IU-1)    = K-1
                    KBMIN(IU-1) = MIN(KB(IU),KB(IU-1))
                    IEXIT       = 1
                    EXIT
                  END IF
                  ELR2 = EL(K+1,IU)+SINA(JB)*DLX(IU)*0.5
                  IF (ELR2 >= ELL2) THEN
                    B(K,IU-1) = B(KUP,UHS(JB)); EXIT
                  ELSE
                    K1 = KUP+1
                    IF (K1 > KMX) EXIT
                    B11 = 0.0
                    EL1 = ELR
                    EL2 = EL(K1,UHS(JB))-SINA(JBUH(JB))*DLX(IU)*0.5
                    DO WHILE (ELR2 <= EL2)
                      B11 = B11+(EL1-EL2)*B(K1-1,UHS(JB))
                      EL1 = EL2
                      K1  = K1+1
                      IF (K1 >= KMX+1 .OR. EL2 == ELR2) EXIT
                      EL2 = EL(K1,UHS(JB))-SINA(JBUH(JB))*DLX(UHS(JB))*0.5
                      IF (EL2 <= ELR2) EL2 = ELR2
                    END DO
                    B(K,IU-1) = B11/H(K,JW); EXIT
                  END IF
                END IF
              END DO
              IF (EL(KMX,UHS(JB)) > EL(K,IU)) B(K,IU-1) = B(K-1,IU-1)
              IF (B(K,IU-1) == 0.0) B(K,IU-1) = B(K-1,IU-1)
              IF (IEXIT == 1) EXIT
            END IF
          END IF
        END IF
      END DO
      IEXIT = 0
      DO K=1,KMX-1
        B(K,ID+1) = B(K,ID)
        IF (DH_INTERNAL(JB)) THEN
          IF (JBDH(JB) >= BS(JW) .AND. JBDH(JB) <= BE(JW)) THEN
            B(K,ID+1) = B(K,DHS(JB))
          ELSE
            ELL = EL(K,ID)-SINA(JB)*DLX(ID)*0.5
            ELR = EL(2,DHS(JB))+SINA(JBDH(JB))*DLX(DHS(JB))*0.5
            IF (ELL >= ELR) THEN
              B(K,ID+1) = B(2,DHS(JB))
            ELSE
              DO KDN=2,KMX-1
                ERR1 = EL(KDN,DHS(JB))  +SINA(JBDH(JB))*DLX(DHS(JB))*0.5
                ERR2 = EL(KDN+1,DHS(JB))+SINA(JBDH(JB))*DLX(DHS(JB))*0.5
                IF (ERR1 >= ELL .AND. ERR2 < ELL) THEN
                  IF (KDN > KB(DHS(JB))) THEN
                    KB(ID+1)  = K-1
                    KBMIN(ID) = MIN(KB(ID),KB(ID+1))
                    IEXIT     = 1
                    EXIT
                  END IF
                  ELL2 = EL(K+1,ID)-SINA(JB)*DLX(ID)*0.5
                  IF (ELL2 >= ERR2) THEN
                    B(K,ID+1) = B(KDN,DHS(JB)); EXIT
                  ELSE
                    K1  = KDN+1
                    IF (K1 > KMX) EXIT
                    B11 = 0.0
                    EL2 = ELL
                    EL1 = EL(K1,DHS(JB))+SINA(JBDH(JB))*DLX(DHS(JB))*0.5
                    DO WHILE (ELL2 <= EL1)
                      B11 = B11+(EL2-EL1)*B(K1-1,DHS(JB))
                      EL2 = EL1
                      K1  = K1+1
                      IF (K1 >= KMX+1 .OR. EL1 == ELL2) EXIT
                      EL1 = EL(K1,DHS(JB))+SINA(JBDH(JB))*DLX(DHS(JB))*0.5
                      IF (EL1 <= ELL2) EL1 = ELL2
                    END DO
                    B(K,ID+1) = B11/H(K,JW); EXIT
                  END IF
                END IF
              END DO
              IF (EL(KMX,DHS(JB)) > EL(K,ID)) B(K,ID+1) = B(K-1,ID+1)
              IF (B(K,ID+1) == 0.0) B(K,ID+1) = B(K-1,ID+1)
              IF (IEXIT == 1) EXIT
            END IF
          END IF
        END IF
      END DO

!**** AREAS AND BOTTOM WIDTHS

      IF (.NOT. TRAPEZOIDAL(JW)) THEN                                                                                  !SW 07/16/04
        DO I=IU-1,ID+1
          DO K=1,KMX-1
            BH2(K,I) = B(K,I)*H(K,JW)
            BH(K,I)  = B(K,I)*H(K,JW)
            BB(K,I)  = B(K,I)-(B(K,I)-B(K+1,I))/(0.5*(H(K,JW)+H(K+1,JW)))*H(K,JW)*0.5                                  !SW 08/02/04
          END DO
          BH(KMX,I) = BH(KMX-1,I)
        END DO
!****** DERIVED GEOMETRY

        DO I=IU-1,ID+1
          BH2(KT,I) = B(KTI(I),I)*(EL(KT,I)-EL(KTI(I)+1,I)-Z(I)*COSA(JB))/COSA(JB)
          IF (KT == KTI(I)) BH2(KT,I) = H2(KT,I)*B(KT,I)
          DO K=KTI(I)+1,KT
            BH2(KT,I) = BH2(KT,I)+BH(K,I)
          END DO
          BKT(I)   = BH2(KT,I)/H2(KT,I)
          BI(KT,I) = B(KTI(I),I)
        END DO
      ELSE                                                                                                             !SW 07/16/04
        DO I=IU-1,ID+1
          DO K=1,KMX-1
            BB(K,I)  = B(K,I)-(B(K,I)-B(K+1,I))/(0.5*(H(K,JW)+H(K+1,JW)))*H(K,JW)*0.5
          END DO
          BB(KB(I),I) = B(KB(I),I)*0.5
          BH2(1,I)    = B(1,I)*H(1,JW)
          BH(1,I)     = BH2(1,I)
          DO K=2,KMX-1
            BH2(K,I) = 0.25*H(K,JW)*(BB(K-1,I)+2.*B(K,I)+BB(K,I))
            BH(K,I)  = BH2(K,I)
          END DO
          BH(KMX,I) = BH(KMX-1,I)
        END DO
        DO I=IU-1,ID+1
          CALL GRID_AREA1(EL(KT,I)-Z(I),EL(KT+1,I),BH2(KT,I),BI(KT,I))
          BKT(I) = BH2(KT,I)/H2(KT,I)
        END DO
      END IF
      DO I=IU-1,ID
        DO K=1,KMX-1
          AVH2(K,I) = (H2(K,I) +H2(K+1,I)) *0.5
          AVHR(K,I) =  H2(K,I)+(H2(K,I+1)-H2(K,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                                  !SW 07/29/04
        END DO
        AVH2(KMX,I) = H2(KMX,I)
        DO K=1,KMX
          BR(K,I)   = B(K,I)  +(B(K,I+1)  -B(K,I))  /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                                !SW 07/29/04
          BHR(K,I)  = BH(K,I) +(BH(K,I+1) -BH(K,I)) /(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                                !SW 07/29/04
          BHR2(K,I) = BH2(K,I)+(BH2(K,I+1)-BH2(K,I))/(0.5*(DLX(I)+DLX(I+1)))*0.5*DLX(I)                                !SW 07/29/04
           IF(CONSTRICTION(K,I))THEN    ! SW 6/26/2018
              IF(BR(K,I) > BCONSTRICTION(I))          BR(K,I)  = BCONSTRICTION(I)
              IF(BHR(K,I) > BCONSTRICTION(I)*H(K,JW)) BHR(K,I) = BCONSTRICTION(I)*H(K,JW)
              IF(BHR2(K,I) > BCONSTRICTION(I)*H(K,JW))BHR2(K,I)= BCONSTRICTION(I)*H(K,JW)
          ENDIF

        END DO
      END DO
      DO K=1,KMX-1
        AVH2(K,ID+1) = (H2(K,ID+1)+H2(K+1,ID+1))*0.5
        BR(K,ID+1)   =   B(K,ID+1)
        BHR(K,ID+1)  =   BH(K,ID+1)
      END DO
      AVH2(KMX,ID+1) = H2(KMX,ID+1)
      AVHR(KT,ID+1)  = H2(KT,ID+1)
      BHR2(KT,ID+1)  = BH2(KT,ID+1)
      IUT = IU
      IF (UP_HEAD(JB)) IUT = IU-1
      DO I=IUT,ID
        DO K=1,KMX-1
          VOL(K,I) = B(K,I)*H2(K,I)*DLX(I)
        END DO
        VOL(KT,I)    = BH2(KT,I)*DLX(I)
        DEPTHB(KT,I) = H2(KT,I)
        DEPTHM(KT,I) = H2(KT,I)*0.5
        DO K=KT+1,KMX
          DEPTHB(K,I) = DEPTHB(K-1,I)+ H2(K,I)
          DEPTHM(K,I) = DEPTHM(K-1,I)+(H2(K-1,I)+H2(K,I))*0.5
        END DO
      END DO
    END DO
  END DO

10 CONTINUE

    RETURN
    END SUBROUTINE POROSITY

!************************************************************************
!**               S U B R O U T I N E    MACROPHYTE_FRICTION           **
!************************************************************************

    SUBROUTINE MACROPHYTE_FRICTION(HRAD,BEDFR,EFFRIC,K,II)

    USE GEOMC;USE GLOBAL;USE MACROPHYTEC;USE POROSITYC
    IMPLICIT NONE
      INTEGER :: K,M,II
      REAL(R8) :: BEDFR, EFFRIC, HRAD
      REAL     :: SAVOLRAT,XSAREA,TSAREA,ARTOT,SCTOT,CDAVG,FRIN

  DO M=1,NMC
    SAVOLRAT=DWV(M)/DWSA(M)     !CB 6/29/2006
    IF(K.EQ.KT)THEN
!      SAREA(M)=VSTEMKT(II,M)*SAVOLRAT/PI
      SAREA(M)=VSTEMKT(II,M)*SAVOLRAT*ANORM(M)     !CB 6/29/2006
    ELSE
!      SAREA(M)=VSTEM(K,II,M)*SAVOLRAT/PI
      SAREA(M)=VSTEM(K,II,M)*SAVOLRAT*ANORM(M)     !CB 6/29/2006
    END IF
  END DO
  XSAREA=BH2(K,II)

  TSAREA=0.0
  ARTOT=0.0
  SCTOT=0.0
  DO M=1,NMC
    ARTOT=ARTOT+SAREA(M)
    SCTOT=SCTOT+CDDRAG(M)*SAREA(M)
    TSAREA=TSAREA+SAREA(M)
  END DO

  IF(ARTOT.GT.0.0)THEN
    CDAVG=SCTOT/ARTOT
    FRIN=CDAVG*TSAREA*HRAD**(4./3.)/(2.0*G*XSAREA*DLX(II)*BEDFR**2)
    EFFRIC=BEDFR*SQRT(1.0+FRIN)
  ELSE
    EFFRIC=BEDFR
  END IF

  RETURN
  END SUBROUTINE MACROPHYTE_FRICTION
