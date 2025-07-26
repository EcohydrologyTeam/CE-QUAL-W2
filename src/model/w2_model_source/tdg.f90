
!***********************************************************************************************************************************
!**                                       S U B R O U T I N E   T O T A L  D I S S O L V E D  G A S                               **
!***********************************************************************************************************************************

SUBROUTINE TOTAL_DISSOLVED_GAS (NSAT,P,NSG,N,T,C)
  USE TDGAS; USE STRUCTURES; USE GLOBAL; USE MAIN, ONLY: EA; USE TVDC, ONLY: TDEW
  IMPLICIT NONE
  INTEGER      :: N, NSG, NSAT
  REAL(R8)     :: T,P,C
  REAL         :: SAT,DB,DA,TDG

  IF(NSAT==0)THEN    ! DISSOLVED OXYGEN
      SAT = EXP(7.7117-1.31403*(LOG(T+45.93)))*P
  ELSE IF(NSAT==1)THEN                ! N2 GAS
  EA = DEXP(2.3026D0*(7.5D0*TDEW(JW)/(TDEW(JW)+237.3D0)+0.6609D0))*0.001316   ! in mm Hg   0.0098692atm=7.5006151mmHg     
   SAT=(1.5568D06*0.79*(P-EA)*(1.8816D-5 - 4.116D-7 * T + 4.6D-9 * T*T))   ! SW 10/27/15    4/20/16 SPEED 
  ELSE
    SAT = p    
  ENDIF
  
  IF (NSG == 0) THEN
    IF (EQSP(N) == 1) THEN
      TDG = AGASSP(N)*.035313*QSP(N)+BGASSP(N)
      IF (TDG > 145.0) TDG = 145.0
      C = SAT
      IF (TDG >= 100.0) C = TDG*SAT/100.0    
    ELSE IF (EQSP(N) == 2) THEN
      TDG = AGASSP(N)+BGASSP(N)*EXP(0.035313*QSP(N)*CGASSP(N))
      IF (TDG > 145.0)  TDG = 145.0    
      C = SAT
      IF (TDG >= 100.0) C = TDG*SAT/100.0
    ELSE IF (EQSP(N) == 3) THEN
      DA = SAT-C                                                                               ! MM 5/21/2009 DA: Deficit upstream
      DB = DA/(1.0+0.38*AGASSP(N)*BGASSP(N)*CGASSP(N)*(1.0-0.11*CGASSP(N))*(1.0+0.046*T))      ! DB: deficit downstream
      C  = SAT-DB
    ELSE IF (EQSP(N) == 4) THEN  
      C = AGASSP(N) * C
      IF(CGASSP(N)==1  .AND. C > SAT)C= SAT
    END IF
  ELSE IF (EQGT(N) == 1) THEN
    TDG = AGASGT(N)*0.035313*QGT(N)+BGASGT(N)
    IF (TDG > 145.0) TDG = 145.0
    C = SAT
    IF (TDG >= 100.0) C = TDG*SAT/100.0
  ELSE IF (EQGT(N) == 2) THEN
    TDG = AGASGT(N)+BGASGT(N)*EXP(.035313*QGT(N)*CGASGT(N))
    IF (TDG > 145.0) TDG = 145.0
    C = SAT
    IF (TDG >= 100.0) C = TDG*SAT/100.0
  ELSE IF (EQGT(N) == 3) THEN  
    DA = SAT-C                                                                               ! MM 5/21/2009 DA: Deficit upstream
    DB = DA/(1.0+0.38*AGASGT(N)*BGASGT(N)*CGASGT(N)*(1.0-0.11*CGASGT(N))*(1.0+0.046*T))      ! DB: deficit downstream
    C  = SAT-DB
  ELSE IF (EQGT(N) == 4) THEN  
      C = AGASGT(N) * C  
      IF(CGASGT(N)==1  .AND. C > SAT)C= SAT
  END IF
END SUBROUTINE TOTAL_DISSOLVED_GAS
