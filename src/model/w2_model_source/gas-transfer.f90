
!***********************************************************************************************************************************
!**                                          S U B R O U T I N E   G A S   T R A N S F E R                                        **
!***********************************************************************************************************************************

SUBROUTINE GAS_TRANSFER
  USE GLOBAL; USE GEOMC; USE KINETIC
  IMPLICIT NONE
  REAL, PARAMETER :: THETA_REAERATION = 1.024, M_TO_FT = 3.2808
  REAL :: AREA,ADEPTH,UAVG,HDEPTH,S,USTAR,A,BCOEF,DMO2
  INTEGER :: K

  IF (REAERC(JW) == '   RIVER') THEN

!** Average depth in ft

    AREA = 0.0
    DO K=KT,KBMIN(I)
      AREA = AREA+BHR1(K,I)
    END DO
    ADEPTH = AREA/BR(KTI(I),I)*M_TO_FT

!** Average velocity in feet/second

    UAVG = ABS(QC(I))/AREA*M_TO_FT

!** Reaeration factor

    IF (NEQN(JW) == 0) THEN
      IF (ADEPTH <= 2.0) THEN
        REAER(I) = 21.64*UAVG**0.67/ADEPTH**1.85
      ELSE IF (UAVG <= 1.8) THEN
        REAER(I) = 12.96*SQRT(UAVG)/ADEPTH**1.5
      ELSE
        HDEPTH = -11.875*UAVG+23.375
        IF (HDEPTH >= ADEPTH) THEN
          REAER(I) = 12.96*SQRT(UAVG)/ADEPTH**1.5
        ELSE
          REAER(I) = 11.57*UAVG**0.969/ADEPTH**1.673
        END IF
      END IF
    ELSE IF (NEQN(JW) == 1) THEN                                                                              !O'connor-Dobbins
      REAER(I) = 12.96*SQRT(UAVG)/ADEPTH**1.5                                                                 ! units: day-1
    ELSE IF (NEQN(JW) == 2) THEN                                                                              !Churchill
      REAER(I) = 11.57*UAVG**0.969/ADEPTH**1.673                                                              ! units: day-1
    ELSE IF (NEQN(JW) == 3) THEN                                                                              !Tsivoglou
      S = SLOPEC(JB)*5280.0
      IF (ABS(QC(I))*35.5 >= 10.0) THEN
        REAER(I) = 0.88*S*UAVG
      ELSE
        REAER(I) = 1.8*S*UAVG
      END IF
    ELSE IF (NEQN(JW) == 4) THEN                                                                              !Owens
      REAER(I) = 21.64*UAVG**0.67/ADEPTH**1.85
    ELSE IF (NEQN(JW) == 5) THEN                                                                              !Thackston and Krenkel
      USTAR    = SQRT(ADEPTH*SLOPEC(JB)*32.2)                                                                  ! SR 5/10/05
      REAER(I) = 24.88*(1.0+SQRT(0.176*UAVG/SQRT(ADEPTH)))*USTAR/ADEPTH                                       ! SR 5/10/05
    ELSE IF (NEQN(JW) == 6) THEN                                                                              !Langbien and Durum
      REAER(I) = 7.60*UAVG/ADEPTH**1.33
    ELSE IF (NEQN(JW) == 7) THEN                                                                              !Melching and Flores
      UAVG = UAVG/M_TO_FT
      IF (QC(I) == 0.0) THEN
        REAER(I) = 0.0
      ELSE IF (ABS(QC(I)) < 0.556) THEN
        REAER(I) = 517.0*((UAVG*SLOPEC(JB))**0.524)*ABS(QC(I))**(-0.242)
      ELSE
        REAER(I) = 596.0*((UAVG*SLOPEC(JB))**0.528)*ABS(QC(I))**(-0.136)
      END IF
    ELSE IF (NEQN(JW) == 8) THEN                                                                              !Melching and Flores
      UAVG   = UAVG/M_TO_FT
      ADEPTH = ADEPTH/M_TO_FT
      IF (ABS(QC(I)) < 0.556) THEN
        REAER(I) = 88.0*((UAVG*SLOPEC(JB))**0.313)*ADEPTH**(-0.353)
      ELSE
        REAER(I) = 142.0*((UAVG*SLOPEC(JB))**0.333)*ADEPTH**(-0.66)*BI(KT,I)**(-0.243)
      END IF
    ELSE IF (NEQN(JW) == 9) THEN                                                                              !User defined SI units
      UAVG     = UAVG/M_TO_FT
      ADEPTH   = ADEPTH/M_TO_FT
      REAER(I) = RCOEF1(JW)*(UAVG**RCOEF2(JW))*(ADEPTH**RCOEF3(JW))*(SLOPEC(JB)**RCOEF4(JW))
    ELSE IF (NEQN(JW) == 10) THEN                                                                             ! Thackston and Krenkel - updated
      USTAR    = SQRT(ADEPTH*SLOPEC(JB)*32.2)                                                                  ! SR 5/10/05
      REAER(I) = 4.99*(1.0+9.0*(0.176*UAVG/SQRT(ADEPTH))**0.25)*USTAR/ADEPTH                                  ! SR 5/10/05
    END IF
    !REAER(I) = REAER(I)*ADEPTH/M_TO_FT - now eliminated since changed computation in DO water quality section
    IF(MINKL(JW)>REAER(I))REAER(I)=MINKL(JW)    ! minimum value units day-1
  ELSE IF (REAERC(JW) == '    LAKE') THEN
    IF (NEQN(JW) == 1) THEN                                                                                   !Broecker
      REAER(I) = 0.864*WIND10(I)
    ELSE IF (NEQN(JW) == 2) THEN
      IF (WIND10(I) <= 3.5) THEN                                                                              !Gelda
        A     = 0.2
        BCOEF = 1.0
      ELSE
        A     = 0.057
        BCOEF = 2.0
      END IF
      REAER(I) = A*WIND10(I)**BCOEF
    ELSE IF (NEQN(JW) == 3) THEN                                                                              !Banks & Herrera
      REAER(I) = (0.728*SQRT(WIND10(I))-0.317*WIND10(I)+0.0372*WIND10(I)**2)                                   !units of m/day
    ELSE IF (NEQN(JW) == 4) THEN                                                                              !Wanninkhof
      REAER(I) = 0.0986*WIND10(I)**1.64                                                                         !units of m/day
    ELSE IF (NEQN(JW) == 5) THEN                                                                              !Chen & Kanwisher
      DMO2     = 2.04E-9
      REAER(I) = DAY*DMO2/((200.0-60.0*SQRT(MIN(WIND10(I),11.0)))*1.E-6)
    ELSE IF (NEQN(JW) == 6) THEN                                                                              !Cole & Buchak
      REAER(I) = (0.5+0.05*WIND10(I)*WIND10(I))
    ELSE IF (NEQN(JW) == 7) THEN                                                                              !Banks
      IF (WIND10(I) <= 5.5) THEN
        REAER(I) = 0.362*SQRT(WIND10(I))
      ELSE
        REAER(I) = 0.0277*WIND10(I)**2
      END IF
    ELSE IF (NEQN(JW) == 8) THEN                                                                              !Smith
      REAER(I) = 0.64+0.128*WIND10(I)**2
    ELSE IF (NEQN(JW) == 9) THEN                                                                              !Liss
      IF (WIND10(I) <= 4.1) THEN
        REAER(I) = 0.156*WIND10(I)**0.63
      ELSE
        REAER(I) = 0.0269*WIND10(I)**1.9
      END IF
    ELSE IF (NEQN(JW) == 10) THEN                                                                             !Downing and Truesdale
      REAER(I) = 0.0276*WIND10(I)**2
    ELSE IF (NEQN(JW) == 11) THEN                                                                             !Kanwisher
      REAER(I) = 0.0432*WIND10(I)**2
    ELSE IF (NEQN(JW) == 12) THEN                                                                             !Yu, et al
      REAER(I) = 0.319*WIND10(I)
    ELSE IF (NEQN(JW) == 13) THEN                                                                             !Weiler
      IF (WIND10(I) <= 1.6) THEN
        REAER(I) = 0.398
      ELSE
        REAER(I) = 0.155*WIND10(I)**2
      END IF
    ELSE IF (NEQN(JW) == 14) THEN                                                                             !User defined
      REAER(I) = RCOEF1(JW)+RCOEF2(JW)*WIND10(I)**RCOEF3(JW)
    END IF
    IF(MINKL(JW)>REAER(I))REAER(I)=MINKL(JW)    ! minimum value units m day-1
    REAER(I) = REAER(I)*BI(KT,I)/BH2(KT,I)                                                         ! conversion from m/d to 1/d for all wind based equations
  ELSE IF (REAERC(JW) == ' ESTUARY') THEN
    AREA = 0.0
    DO K=KT,KBMIN(I)
      AREA = AREA+BHR1(K,I)
    END DO
 
    ADEPTH = AREA/BR(KTI(I),I)
    UAVG   = ABS(QC(I))/AREA

!** Reaeration factor

    IF (NEQN(JW) == 0) THEN
    ADEPTH = ADEPTH*M_TO_FT
    UAVG   = UAVG*M_TO_FT
      IF (ADEPTH <= 2.0) THEN
        REAER(I) = 21.64*UAVG**0.67/ADEPTH**1.85                                        ! units of day-1
      ELSE IF (UAVG <= 1.8) THEN
        REAER(I) = 12.96*SQRT(UAVG)/ADEPTH**1.5
      ELSE
        HDEPTH = -11.875*UAVG+23.375
        IF (HDEPTH >= ADEPTH) THEN
          REAER(I) = 12.96*SQRT(UAVG)/ADEPTH**1.5
        ELSE
          REAER(I) = 11.57*UAVG**0.969/ADEPTH**1.673
        END IF
      END IF
    ELSE IF (NEQN(JW) == 1) THEN                                                                           !Thomann and Fitzpatrick
      REAER(I) = (0.728*SQRT(WIND10(I))-0.317*WIND10(I)+0.0372*WIND10(I)**2)+3.93*SQRT(UAVG)/(ADEPTH)**0.5 ! units of m/day
      REAER(I) = REAER(I)/ADEPTH                                                                           ! units of 1/day
    ELSE IF (NEQN(JW) == 2) THEN                                                                           !User defined
      REAER(I) = RCOEF1(JW)*(UAVG**RCOEF2(JW))*(ADEPTH**RCOEF3(JW))+(0.5+RCOEF4(JW)*WIND10(I)*WIND10(I))   ! units of 1/day
    ENDIF
    IF(MINKL(JW)>REAER(I))REAER(I)=MINKL(JW)    ! minimum value units day-1

  END IF
  !IF(RCOEF1(JW)>REAER(I))REAER(I)=RCOEF1(JW)    ! minimum value units day-1
  REAER(I) = REAER(I)*THETA_REAERATION**(T1(KT,I)-20.0)      ! units of day-1
  REAER(I) = REAER(I)/DAY
END SUBROUTINE GAS_TRANSFER
