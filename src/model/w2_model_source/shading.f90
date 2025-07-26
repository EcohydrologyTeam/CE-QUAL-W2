
!***********************************************************************************************************************************
!**                                                S U B R O U T I N E   S H A D I N G                                            **
!***********************************************************************************************************************************

SUBROUTINE SHADING
  USE SHADEC; USE GLOBAL; USE GDAYC; USE SURFHE; USE GEOMC; USE SCREENC; USE LOGICC
  IMPLICIT NONE
  CHARACTER(1) :: BANK
  REAL         :: LOCAL,STANDARD,HOUR,TAUD,SINAL,A02,AZ00,A0,AX,ANG1,ANG2,TOPOANG,SFACT,HT,CLINE,SRED,STLEN,EDGE,EDAZ,SN,AZT
  INTEGER      :: IDAY,J

! Calculate solar altitude, declination, and local hour angle when short-wave solar radiation is provided as input

  IF (READ_RADIATION(JW)) THEN
    LOCAL    =  LONGIT(JW)
    STANDARD =  15.0*INT(LONGIT(JW)/15.0)
    HOUR     = (JDAY-INT(JDAY))*24.0
    IDAY     =  JDAY-((INT(JDAY/365))*365)
    !IDAY     =  IDAY+INT(INT(JDAY/365)/4)    SR 12/2018 IDAY FIX
    IDAY     =  IDAY-INT(INT(JDAY/365)/4)
    TAUD     = (2*PI*(IDAY-1))/365
    EQTNEW   =  0.170*SIN(4*PI*(IDAY-80)/373)-0.129*SIN(2*PI*(IDAY-8)/355)
    HH(JW)   =  0.261799*(HOUR-(LOCAL-STANDARD)*0.0666667+EQTNEW-12.0)
    DECL(JW) =  0.006918-0.399912*COS(TAUD)+0.070257*SIN(TAUD)-0.006758*COS(2*TAUD)+0.000907*SIN(2*TAUD)-0.002697*COS(3*TAUD)      &
                +0.001480*SIN(3*TAUD)
    SINAL    =  SIN(LAT(JW)*.0174533)*SIN(DECL(JW))+COS(LAT(JW)*.0174533)*COS(DECL(JW))*COS(HH(JW))
    A00(JW)  =  57.2957795*ASIN(SINAL)
  END IF

! If the sun is below the horizon, set SHADE(I) to 0

  IF (A00(JW) < 0.0) THEN
    SHADE(I) = 0.0
  ELSE

!** Calculate solar azimuth angle

    A02 = A00(JW)/57.2957795
    AX  = (SIN(DECL(JW))*COS(LAT(JW)*0.017453)-COS(DECL(JW))*COS(HH(JW))*SIN(LAT(JW)*0.017453))/COS(A02)
    IF (AX >  1.0) AX =  1.0
    IF (AX < -1.0) AX = -1.0
    AZT = ACOS(AX)
    IF (HH(JW) < 0.0) THEN
     AZ00 = AZT
    ELSE
     AZ00 = 2.0*PI-AZT
    END IF
    A0 = A02

!** Interpolate the topographic shade angle

    DO J=1,IANG-1
      IF (AZ00 > ANG(J) .AND. AZ00 <= ANG(J+1)) THEN
        ANG1    =  AZ00-ANG(J)
        ANG2    = (TOPO(I,J+1)-TOPO(I,J))/GAMA                 ! SW 10/17/05
        TOPOANG =  TOPO(I,J)+ANG2*ANG1
      END IF
    END DO
    IF (AZ00 > ANG(IANG) .AND. AZ00 <= 2*PI) THEN
      ANG1    =  AZ00-ANG(IANG)
      ANG2    = (TOPO(I,1)-TOPO(I,IANG))/GAMA                  ! SW 10/17/05
      TOPOANG =  TOPO(I,IANG)+ANG2*ANG1
    END IF

!** Complete topographic shading if solar altitude less than topo angle

    IF (A0 <= TOPOANG) THEN
      SFACT = 0.90
      GO TO 100
    END IF

!** No vegetative shading if azimuth angle is oriented parallel to stream

    IF (AZ00 == PHI0(I) .OR. AZ00 == PHI0(I)+PI .OR. AZ00+PI == PHI0(I)) THEN
      SFACT = 0.0
      GO TO 100
    END IF

!** Bank with the controlling vegetation

    IF (PHI0(I) > 0.0 .AND. PHI0(I) <= PI) THEN
      IF (AZ00 > PHI0(I)     .AND. AZ00 <= PHI0(I)+PI) BANK = 'L'
      IF (AZ00 > 0.0         .AND. AZ00 <= PHI0(I))    BANK = 'R'
      IF (AZ00 > PHI0(I)+PI  .AND. AZ00 <  2.0*PI)     BANK = 'R'
    ELSE IF (PHI0(I) > PI .AND. PHI0(I) <= 2.0*PI) THEN
      IF (AZ00 >= PHI0(I)    .AND. AZ00 < 2.0*PI)      BANK = 'L'
      IF (AZ00 >= 0.0        .AND. AZ00 < PHI0(I)-PI)  BANK = 'L'
      IF (AZ00 >= PHI0(I)-PI .AND. AZ00 < PHI0(I))     BANK = 'R'
    END IF

!** No topographic shading

    IF (BANK == 'L') THEN
      IF (TTLB(I) < ELWS(I)) THEN
        SFACT = 0.0
        GO TO 100
      ELSE
        HT    = TTLB(I)-ELWS(I)
        CLINE = CLLB(I)
        SRED  = SRLB2(I)
        IF (JDAYG > SRFJD1(I) .AND. JDAYG <= SRFJD2(I)) SRED = SRLB1(I)
      END IF
    ELSE
      IF (TTRB(I) < ELWS(I)) THEN
        SFACT = 0.0
        GO TO 100
      ELSE
        HT    = TTRB(I)-ELWS(I)
        CLINE = CLRB(I)
        SRED  = SRRB2(I)
        IF (JDAYG > SRFJD1(I) .AND. JDAYG <= SRFJD2(I)) SRED = SRRB1(I)
      END IF
    END IF
    STLEN = HT/TAN (A0)
    EDGE  = MAX (0.0,CLINE-BI(KT,I)/2.0)

!** Distance from vegetation to water edge on line parallel to azimuth

    EDAZ = EDGE/ABS(SIN(PHI0(I)-AZ00))
    IF (STLEN <= EDAZ) THEN
      SFACT = 0.0
      GO TO 100
    END IF

!** Distance shadow extends over water (perpendicular to segment orientation)

    SN    = MIN (HT*ABS (SIN (ABS (PHI0(I)-AZ00)))/TAN (A0)-EDGE,BI(KT,I))
    SFACT = SRED*SN/BI(KT,I)
100 CONTINUE
    SHADE(I) = MAX (0.0,1.0-SFACT)
    SHADE(I) = MIN(ABS(SHADEI(I)),SHADE(I))              ! SW 10/2/2017 Allows for fixed canopy cover over top of channel - only used if shade is less than shadei only valid for -0.99 and 0.0
  END IF
  RETURN
END SUBROUTINE SHADING
