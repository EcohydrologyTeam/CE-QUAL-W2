
!***********************************************************************************************************************************
!**                                        S U B R O U T I N E   H E A T  E X C H A N G E                                         **
!***********************************************************************************************************************************

SUBROUTINE HEAT_EXCHANGE
  USE GLOBAL; USE GDAYC; USE SURFHE; USE TVDC; USE SHADEC; USE PREC
  IMPLICIT NONE

! Type declaration
  REAL     :: JDAY, LOCAL
  REAL(R8) :: TSUR
  REAL(R8) :: MPS_TO_MPH=2.23714D0,W_M2_TO_BTU_FT2_DAY=7.60796D0,FLUX_BR_TO_FLUX_SI=0.23659D0,BOWEN_CONSTANT=0.47D0
  REAL(R8) :: DEG_F,DEG_C,BTU_FT2_DAY_TO_W_M2=0.1314D0
  REAL(R8) :: STANDARD,HOUR,TAUD,SINAL,A0,TDEW_F,TAIR_F,SRO_BR,WIND_MPH,WIND2M,ACONV,BCONV
  REAL(R8) :: TSTAR,BETA,FW,RA,ETP,EA,ES,TAIRV,DTV,DTVL
  INTEGER  :: IDAY,J

RETURN

!***********************************************************************************************************************************
!**                                            S H O R T  W A V E  R A D I A T I O N                                              **
!***********************************************************************************************************************************

ENTRY SHORT_WAVE_RADIATION (JDAY)
  LOCAL    =  LONGIT(JW)
  STANDARD =  15.0*INT(LONGIT(JW)/15.0)
  HOUR     = (JDAY-INT(JDAY))*24.0
  IDAY     =  JDAY-((INT(JDAY/365))*365)
  !IDAY     =  IDAY+INT(INT(JDAY/365)/4)  ! SR 12/2018 IDAY FIX
  IDAY     =  IDAY-INT(INT(JDAY/365)/4)
  TAUD     = (2.*PI*(IDAY-1))/365.
  EQTNEW   =  0.170*SIN(4.*PI*(IDAY-80)/373.)-0.129*SIN(2.*PI*(IDAY-8)/355.)
  HH(JW)   =  0.261799*(HOUR-(LOCAL-STANDARD)*0.0666667+EQTNEW-12.0)
  DECL(JW) =  0.006918-0.399912*COS(TAUD)+0.070257*SIN(TAUD)-0.006758*COS(2.*TAUD)+0.000907*SIN(2.*TAUD)-0.002697*COS(3.*TAUD)        &
              +0.001480*SIN(3.*TAUD)
  SINAL    =  SIN(LAT(JW)*.0174533)*SIN(DECL(JW))+COS(LAT(JW)*.0174533)*COS(DECL(JW))*COS(HH(JW))
  A00(JW)  =  57.2957795*ASIN(SINAL)
  A0       =  A00(JW)
  IF (A0 > 0.0) THEN
    SRON(JW) = (1.0-0.0065*CLOUD(JW)**2)*24.0*(2.044*A0+0.1296*A0**2-1.941E-3*A0**3+7.591E-6*A0**4)*BTU_FT2_DAY_TO_W_M2
  ELSE
    SRON(JW) = 0.0
  END IF
RETURN

!***********************************************************************************************************************************
!**                                          E Q U I L I B R I U M  T E M P E R A T U R E                                         **
!***********************************************************************************************************************************

ENTRY EQUILIBRIUM_TEMPERATURE

! British units

  TDEW_F   = DEG_F(TDEW(JW))
  TAIR_F   = DEG_F(TAIR(JW))
  SRO_BR   = SRON(JW)*W_M2_TO_BTU_FT2_DAY*SHADE(I)
  !WIND_MPH = WIND(JW)*WSC(I)*MPS_TO_MPH
  !WIND2M   = WIND_MPH*DLOG(2.0D0/Z0(JW))/DLOG(WINDH(JW)/Z0(JW))+NONZERO     ! SW 11/28/07  old version z0=0.003
  WIND2M=WIND2(I)*MPS_TO_MPH    ! ALREADY COMPUTED IN w2 MAIN PROGRAM
  ACONV    = W_M2_TO_BTU_FT2_DAY
  IF (CFW(JW) == 1.0) BCONV = 3.401062
  IF (CFW(JW) == 2.0) BCONV = 1.520411
  ! SW Issues: CFW not determined for other values of CFW

! Equilibrium temperature and heat exchange coefficient

  ET(I)   =  TDEW_F
  TSTAR   = (ET(I)+TDEW_F)*0.5
  BETA    =  0.255-(8.5E-3*TSTAR)+(2.04E-4*TSTAR*TSTAR)
  FW      =  ACONV*AFW(JW)+BCONV*BFW(JW)*WIND2M**CFW(JW)
  CSHE(I) =  15.7+(0.26+BETA)*FW
  RA      =  3.1872E-08*(TAIR_F+459.67)**4
  ETP     = (SRO_BR+RA-1801.0)/CSHE(I)+(CSHE(I)-15.7)*(0.26*TAIR_F+BETA*TDEW_F)/(CSHE(I)*(0.26+BETA))
  J       =  0
  DO WHILE (ABS(ETP-ET(I)) > 0.05 .AND. J < 10)
    ET(I)   =  ETP
    TSTAR   = (ET(I)+TDEW_F)*0.5
    BETA    =  0.255-(8.5E-3*TSTAR)+(2.04E-4*TSTAR*TSTAR)
    CSHE(I) =  15.7+(0.26+BETA)*FW
    ETP     = (SRO_BR+RA-1801.0)/CSHE(I)+(CSHE(I)-15.7)*(0.26*TAIR_F+BETA*TDEW_F)/(CSHE(I)*(0.26+BETA))
    J       =  J+1
  END DO

! SI units

  ET(I)   = DEG_C(ET(I))
  CSHE(I) = CSHE(I)*FLUX_BR_TO_FLUX_SI/RHOWCP
RETURN

!***********************************************************************************************************************************
!**                                                   S U R F A C E   T E R M S                                                   **
!***********************************************************************************************************************************

ENTRY SURFACE_TERMS (TSUR)

! Partial water vapor pressure of air (mm hg)

!  EA = EXP(2.3026*(9.5*TDEW(JW)/(TDEW(JW)+265.5)+0.6609))         ! SW 6/10/2011
!  IF (TDEW(JW) > 0.0) EA = EXP(2.3026*(7.5*TDEW(JW)/(TDEW(JW)+237.3)+0.6609))
  EA = DEXP(2.3026D0*(7.5D0*TDEW(JW)/(TDEW(JW)+237.3D0)+0.6609D0))

! Partial water vapor pressure at the water surface

  IF(TSUR<0.0)THEN
  ES = DEXP(2.3026D0*(9.5D0*TSUR/(TSUR+265.5D0)+0.6609D0))
  ELSE
  ES = DEXP(2.3026D0*(7.5D0*TSUR/(TSUR+237.3D0)+0.6609D0))
  ENDIF

! Wind function

  IF (RH_EVAP(JW)) THEN
    TAIRV = (TAIR(JW)+273.0D0)/(1.0D0-0.378D0*EA/760.0D0)
    DTV   = (TSUR+273.0D0)/(1.0D0-0.378D0*ES/760.0D0)-TAIRV
    DTVL  =  0.0084D0*WIND2(I)**3
    IF (DTV < DTVL) DTV = DTVL
    FW = (3.59D0*DTV**0.3333D0+4.26D0*WIND2(I))
  ELSE
    FW = AFW(JW)+BFW(JW)*WIND2(I)**CFW(JW)
  END IF

! Evaporative flux

  RE(I) = FW*(ES-EA)
  !IF(RE(I) < 0.0)RE(I)=0.0     ! SW 6/22/2016  SHOULD WE USE THIS AS AN ANALOG FOR CONDENSATION WHEN LESS THAN 0? TVA(1972) SUGGESTS YOU CAN - SEE SECTION 4 AND 5 but Ryan, Harleman suggest setting RE=0 if less than zero since the condensation coefficients are unknown...

! Conductive flux

  RC(I) = FW*BOWEN_CONSTANT*(TSUR-TAIR(JW))

! Back radiation flux

  RB(I) = 5.51D-8*(TSUR+273.15D0)**4
  RETURN  
END SUBROUTINE HEAT_EXCHANGE

! Function declaration

   REAL(R8) FUNCTION DEG_F(X)
   USE PREC
   REAL(R8) :: X
   DEG_F =  X*1.8+32.0
   END FUNCTION DEG_F
   REAL(R8) FUNCTION DEG_C(X)
   USE PREC
   REAL(R8) :: X
   DEG_C = (X-32.0)*5.0/9.0
   END FUNCTION DEG_C

