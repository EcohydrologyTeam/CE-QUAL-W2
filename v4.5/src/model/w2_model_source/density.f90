!***********************************************************************************************************************************
!**                                            F U N C T I O N   D E N S I T Y                                                    **
!***********************************************************************************************************************************

REAL(R8) FUNCTION DENSITY (T,TDS,SS)
  use PREC
  USE LOGICC, ONLY: SUSP_SOLIDS, FRESH_WATER, SALT_WATER; USE GLOBAL, ONLY:JW
  IMPLICIT NONE
  REAL(R8) :: T,TDS,SS
                       DENSITY = ((((6.536332D-9*T-1.120083D-6)*T+1.001685D-4)*T-9.09529D-3)*T+6.793952D-2)*T+0.842594D0
  IF (SUSP_SOLIDS)     DENSITY = DENSITY+6.2D-4*SS
  IF (FRESH_WATER(JW)) DENSITY = DENSITY+TDS*((4.99D-8*T-3.87D-6)*T+8.221D-4)
  IF (SALT_WATER(JW))  DENSITY = DENSITY+TDS*((((5.3875D-9*T-8.2467D-7)*T+7.6438D-5)*T-4.0899D-3)*T+0.824493D0)                      &
                                 +((-1.6546D-6*T+1.0227D-4)*T-5.72466D-3)*TDS**1.5D0+4.8314D-4*TDS*TDS
  DENSITY = DENSITY+999.0D0
END FUNCTION DENSITY
