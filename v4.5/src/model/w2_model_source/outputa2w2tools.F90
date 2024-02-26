SUBROUTINE OUTPUTA  
    
  USE MAIN  
  USE GLOBAL;     USE NAMESC; USE GEOMC;  USE LOGICC; USE PREC;  USE SURFHE;  USE KINETIC; USE SHADEC; USE EDDY  
  USE STRUCTURES; USE TRANS;  USE TVDC;   USE SELWC;  USE GDAYC; USE SCREENC; USE TDGAS;   USE RSTART  
  USE MACROPHYTEC; USE POROSITYC; USE ZOOPLANKTONC;   USE BIOENERGETICS
  use CEMAVars
  USE CEMAOutputRoutines; USE ALGAE_TOXINS
  IMPLICIT NONE
  
  EXTERNAL RESTART_OUTPUT  
    
  INTEGER :: JAD,JAF,IFLAG,JWWD,NLINES,ITOT,JJ, NUMOUTLETS,JSSS(100),JJC, IC, KK
  REAL    :: TVOLAVG,QSUMM,CGASD,XDUM,QOUTLET(100),TOUTLET(100),VOLTOT, TBLANK  ! SW 2/28/2020  
  REAL(R8):: DLVBR,DLE,CGAS,TGATE,TSPILL  
    
  ! *** DSI W2_TOOL LINKAGE  
  REAL*4,SAVE,ALLOCATABLE,DIMENSION(:)::WSEL  
  REAL*4,SAVE,ALLOCATABLE,DIMENSION(:,:)::WDSI  
    
IF (VECTOR(1)) THEN  
  IF(.NOT.ALLOCATED(WSEL))THEN  
    ALLOCATE(WSEL(IMX))
    ALLOCATE(WDSI(KMX,IMX))  
  ENDIF  
ENDIF

    
  !***********************************************************************************************************************************  
  !*                                                    Task 2.8: Output Results                                                    **  
  !***********************************************************************************************************************************  
    
IF(LAKE_RIVER_CONTOUR_ON=='ON')THEN    ! SW 2/28/2020
          DO JJ=1,NUM_LAKE_CONTOUR
              IF(JDAY >= NXT_LAKE_CONTOUR(JJ))THEN
                  NXT_LAKE_CONTOUR(JJ)=NXT_LAKE_CONTOUR(JJ)+LAKE_CONTOUR_FREQ(JJ)
                  TBLANK=-99.0
                  IF(LAKE_CONTOUR_FORMAT == 1)THEN
                  WRITE(LAKE_RIVER_CONTOUR+JJ-1,'(F10.4,",",*(F8.2,","))')JDAY,ELWS(LAKE_CONTOUR_SEG(JJ))+0.1,TBLANK
                   WRITE(LAKE_RIVER_CONTOUR+JJ-1,'(F10.4,",",*(F8.2,","))')JDAY,ELWS(LAKE_CONTOUR_SEG(JJ)),T2(KTWB(JW_LAKE_CONTOUR(JJ)),LAKE_CONTOUR_SEG(JJ))
                  DO K=KTWB(JW_LAKE_CONTOUR(JJ))+1,KB(LAKE_CONTOUR_SEG(JJ))
                      WRITE(LAKE_RIVER_CONTOUR+JJ-1,'(F10.4,",",*(F8.2,","))')JDAY,(EL(K,LAKE_CONTOUR_SEG(JJ))+EL(K,LAKE_CONTOUR_SEG(JJ)))*0.5,T2(K,LAKE_CONTOUR_SEG(JJ))
                  ENDDO
                  ELSE
                  WRITE(LAKE_RIVER_CONTOUR+JJ-1,'(F10.4,",",*(F8.2,","))')JDAY,(TBLANK,K=2,KTWB(JW_LAKE_CONTOUR(JJ))-1),(T2(K,LAKE_CONTOUR_SEG(JJ)),K=KTWB(JW_LAKE_CONTOUR(JJ)),KB(LAKE_CONTOUR_SEG(JJ)))
                  ENDIF
                  IF(OXYGEN_DEMAND)THEN
                    IF(LAKE_CONTOUR_FORMAT == 1)THEN
                    WRITE(LAKE_RIVER_CONTOUR+10+JJ-1,'(F10.4,",",*(F8.2,","))')JDAY,ELWS(LAKE_CONTOUR_SEG(JJ))+0.1,TBLANK
                    WRITE(LAKE_RIVER_CONTOUR+10+JJ-1,'(F10.4,",",*(F8.2,","))')JDAY,ELWS(LAKE_CONTOUR_SEG(JJ)),C2(KTWB(JW_LAKE_CONTOUR(JJ)),LAKE_CONTOUR_SEG(JJ),NDO)
                    DO K=KTWB(JW_LAKE_CONTOUR(JJ))+1,KB(LAKE_CONTOUR_SEG(JJ))
                      WRITE(LAKE_RIVER_CONTOUR+10+JJ-1,'(F10.4,",",*(F8.2,","))')JDAY,(EL(K,LAKE_CONTOUR_SEG(JJ))+EL(K,LAKE_CONTOUR_SEG(JJ)))*0.5,C2(K,LAKE_CONTOUR_SEG(JJ),NDO)
                    ENDDO
                    ELSE
                    WRITE(LAKE_RIVER_CONTOUR+10+JJ-1,'(F10.4,",",*(F8.2,","))')JDAY,(TBLANK,K=2,KTWB(JW_LAKE_CONTOUR(JJ))-1),(C2(K,LAKE_CONTOUR_SEG(JJ),NDO),K=KTWB(JW_LAKE_CONTOUR(JJ)),KB(LAKE_CONTOUR_SEG(JJ)))
                    ENDIF                 
                  ENDIF
              ENDIF    
          ENDDO
          DO JJ=1,NUM_RIVER_CONTOUR
              IF(JDAY >= NXT_RIVER_CONTOUR(JJ))THEN
                  NXT_RIVER_CONTOUR(JJ)=NXT_RIVER_CONTOUR(JJ)+RIVER_CONTOUR_FREQ(JJ)
                  IF(RIVER_CONTOUR_FORMAT == 1)THEN
                  DO JB=RIVER_CONTOUR_BR1(JJ),RIVER_CONTOUR_BR2(JJ)
                  DO I=US(JB),DS(JB)
                      WRITE(LAKE_RIVER_CONTOUR+20+JJ-1,'(F10.4,",",*(F8.2,","))')JDAY,X1(I),T2(KTWB(JW_RIVER_CONTOUR(JJ)),I)
                  ENDDO
                  ENDDO
                  ELSE
                  WRITE(LAKE_RIVER_CONTOUR+20+JJ-1,'(F10.4,",",*(F8.2,","))')JDAY,((T2(KTWB(JW_RIVER_CONTOUR(JJ)),I),I=US(JB),DS(JB)),JB=RIVER_CONTOUR_BR1(JJ),RIVER_CONTOUR_BR2(JJ))
                  ENDIF
                    IF(OXYGEN_DEMAND)THEN
                    IF(RIVER_CONTOUR_FORMAT == 1)THEN
                    DO JB=RIVER_CONTOUR_BR1(JJ),RIVER_CONTOUR_BR2(JJ)
                    DO I=US(JB),DS(JB)
                      WRITE(LAKE_RIVER_CONTOUR+30+JJ-1,'(F10.4,",",*(F8.2,","))')JDAY,X1(I),C2(KTWB(JW_RIVER_CONTOUR(JJ)),I,NDO)
                    ENDDO
                    ENDDO
                    ELSE
                    WRITE(LAKE_RIVER_CONTOUR+30+JJ-1,'(F10.4,",",*(F8.2,","))')JDAY,((C2(KTWB(JW_RIVER_CONTOUR(JJ)),I,NDO),I=US(JB),DS(JB)),JB=RIVER_CONTOUR_BR1(JJ),RIVER_CONTOUR_BR2(JJ))
                    ENDIF             
                    ENDIF

              ENDIF
          ENDDO
    ENDIF  
    
  ! BIOEXP MLM  
  IF (BIOEXP) THEN  
    IF (JDAY.GE.NXTBIO) THEN   
      NXTBIO = NXTBIO+BIOF(1)  ! mlm3
      DO J=1,NIBIO  
          
        I = IBIO(J)  
        DO JW=1,NWB  
          IF (I >= US(BS(JW))-1 .AND. I <= DS(BE(JW))+1) EXIT  
        END DO  
        ! bioexp  
        DO K = KTWB(JW),KBI(I)  
          DO JC = NZOOS,NZOOE  
            C2ZOO(K,I,JC) =  C2(K,I,JC)*CMULT(JC)  
          END DO  
	WRITE(BIOEXPFN(J),'(F8.2,",",I8,",",3(F8.2,","),<NZOOE-NZOOS+1>(F8.3,","),I8,",",2(F8.2,","),A,",",I0,",",I0)') JDAY,IBIO(J),DEPTHM(K,I),T1(K,I),GAMMA(K,I),&          ! CB 1/6/17
			     (C2ZOO(K,I,JC), JC = NZOOS,NZOOE),K,BH(K,I),EL(K,I),MONTH,GDAY,YEAR
            
        END DO  
           ! VOLUME WEIGHTING OF ACTIVE CONSTITUENTS     !MLM 18.07.06
        KLIM = MIN(KB(I),KTWB(JW)+5)  
        DO K=KTWB(JW),KLIM  
          DO JJC = 1,NAC  
            C2W(I,JJC) = C2W(I,JJC)+C2(K,I,CN(JJC))*CMULT(CN(JJC))*BH(K,I)  
          END DO  
			     C2W(I,NAC+1) = C2W(I,NAC+1)+CD(K,I,PH_DER)*BH(K,I) !PH
			     C2W(I,NAC+2) = C2W(I,NAC+2)+CD(K,I,TP_DER)*BH(K,I) ! TOTAL PHOSPHOROUS
			 VOLROOS(I)=VOLROOS(I)+BH(K,I)                          ! VOLUME WEIGHTED
        END DO  
        DO JJC = 1,NAC+2  
          C2W(I,JJC) = C2W(I,JJC)/VOLROOS(I)  
        END DO  
          
        WRITE(WEIGHTNUM(J),'(F10.3,",",<NAC>(F9.3,","),2(F9.3,","))')JDAY,(C2W(I,JJC),JJC = 1,NAC+2)  
        VOLROOS = 0.0
          
      END DO  
    END IF  
  END IF  
  
  IF(WLC=='      ON')THEN
      IF(JDAY.GE.NXWL)THEN
        NXWL = NXWL+WLF  
        ! write out water level File  
      WRITE(WLFN,'(f10.3,",",*(f8.3,","))')jday,((elws(i),i=us(jb),ds(jb)),jb=1,nbr)      
      ENDIF
  ENDIF
  
  !** Time series  
    
  IF (TIME_SERIES) THEN  
    IF (JDAY.GE.NXTMTS.OR.JDAY.GE.TSRD(TSRDP+1)) THEN 
      IF (JDAY.GE.TSRD(TSRDP+1)) THEN  
        TSRDP  = TSRDP+1  
        NXTMTS = TSRD(TSRDP)  
      END IF  
      NXTMTS = NXTMTS+TSRF(TSRDP)  
                
      DO J=1,NIKTSR  
        I = ITSR(J)  
        ! find out if segment is inactive OR cell is inactive for fixed layer - do not write out tsr file  ! SW 7/24/2018
        IF(BR_INACTIVE(JBTSR(J)))CYCLE
        IF(CUS(JBTSR(J)) > I)CYCLE
        DO JW=1,NWB  
          IF (I >= US(BS(JW))-1 .AND. I <= DS(BE(JW))+1) EXIT  
        END DO          
        IF(ETSR(J) < 0)THEN                 ! SW 7/24/2018
            IF(INT(ABS(ETSR(J))) < KTWB(JW))CYCLE  
        ENDIF
        
        ! TEMP VOL WEIGHTED AVERAGE  SW 6/1/2015
            TVOLAVG=0.0     ! SW 6/1/2015
            VOLTOT=0.0
                DO K=KTWB(JW),KB(I)   
                VOLTOT=VOLTOT+VOL(K,I)
                TVOLAVG=TVOLAVG+T1(K,I)*VOL(K,I) 
                ENDDO
            IF(KB(I)>=KTWB(JW))tvolavg=tvolavg/voltot
            
        IF (ETSR(J) < 0) THEN  
          K = INT(ABS(ETSR(J)))  
        ELSE  
          DO K=KTWB(JW),KB(I)  
            IF (DEPTHB(K,I) > ETSR(J)) EXIT  
          END DO  
          IF (K > KB(I)) CYCLE  
        END IF  
        
    IF(K >= KTWB(JW))THEN      ! SW 4/4/2018 ADDED TO ELIMINATE THE LAST VALUE OF A VARIABLE BEING USED FOR CASE WHEN ESTR IS NEGATIVE
     
        DO JAC=1,NAC  
          L = LEN_TRIM(FMTC(CN(JAC)))  
          WRITE (C2CH(JAC),FMTC(CN(JAC))(1:L)) C2(K,I,CN(JAC))*CMULT(CN(JAC))  
        END DO  
        DO JF=1,NAF(JW)  
          !WRITE(KFCH(JF),'(E10.3)')KF(K,I,(KFCN(JF,JW)))*VOL(K,I)/1000./DAY
          WRITE(KFCH(JF),'(E10.3)')KF(K,I,(KFCN(JF,JW)))*VOL(K,I)/1000.*DAY   ! BG 12/09/18
        ENDDO  
        DO JAD=1,NACD(JW)  
          L = LEN_TRIM(FMTCD(CDN(JAD,JW)))  
          WRITE (CDCH(JAD),FMTCD(CDN(JAD,JW))(1:L)) CD(K,I,CDN(JAD,JW))*CDMULT(CDN(JAD,JW))  
        END DO  
        DO JE=1,NEP  
          WRITE (EPCH(JE),'(F10.3)') EPD(K,I,JE)                                                    ! SW 8/13/06  
        END DO  
        DO JA=1,NAL
          WRITE (APCH(JA),'(F10.3)') APLIM(K,I,JA)                                                    ! SW 8/13/06  
          WRITE (ANCH(JA),'(F10.3)') ANLIM(K,I,JA)                                                    ! SW 8/13/06  
          WRITE (ALCH(JA),'(F10.3)') ALLIM(K,I,JA)                                                    ! SW 8/13/06  
        END DO  
        DO JM=1,NMC  
          WRITE (macCH(Jm),'(F10.3)') mac(K,I,Jm)                                                   ! SW 8/13/06  
        END DO  
        IF(SEDIMENT_CALC(JW))THEN  
          write (sedch,'(F10.3)') sed(K,I)                                                          ! SW 8/13/06  
          write (sedpch,'(F10.3)') sedp(K,I)  
          write (sednch,'(F10.3)') sedn(K,I)  
          write (sedcch,'(F10.3)') sedc(K,I)  
        END IF  
        IF (ICE_COMPUTATION) THEN  
          IF(SEDIMENT_CALC(JW))THEN  
            WRITE (TSR(J),'(f10.3,",",19(F10.3,","),*(A,","))') JDAY,DLT,ELWS(I),T1(K,I),U(K,I),QC(I),SRON(JW)*1.06,GAMMA(K,I),DEPTHB(KB(I),I),    &     ! SW 8/13/06  
            BI(KTWB(JW),I),SHADE(I),ICETH(I),TVOLAVG,rn(i),rs(i),ranlw(jw),rb(i),re(i),rc(i),REAER(I)*86400.,(ADJUSTR(C2CH(JAC)),JAC=1,NAC),                      &  ! CB 7/26/07
            (ADJUSTR(EPCH(JE)),JE=1,NEP),(ADJUSTR(MACCH(JM)),JM=1,NMC),SEDCH,SEDPCH,SEDNCH,SEDCCH, &  
            (ADJUSTR(CDCH(JAD)),JAD=1,NACD(JW)),(ADJUSTR(KFCH(JF)),JF=1,NAF(JW)),(ADJUSTR(APCH(JA)),JA=1,NAL),(ADJUSTR(ANCH(JA)),JA=1,NAL),(ADJUSTR(ALCH(JA)),JA=1,NAL)    ! SW 10/20/15 
          ELSE  
            WRITE (TSR(J),'(f10.3,",",19(F10.3,","),*(A,","))') JDAY,DLT,ELWS(I),T1(K,I),U(K,I),QC(I),SRON(JW)*1.06,GAMMA(K,I),DEPTHB(KB(I),I),      &     ! SW 8/13/06  
            BI(KTWB(JW),I),SHADE(I),ICETH(I),TVOLAVG,rn(i),rs(i),ranlw(jw),rb(i),re(i),rc(i),REAER(I)*86400.,(ADJUSTR(C2CH(JAC)),JAC=1,NAC),                      &  ! CB 7/26/07
            (ADJUSTR(EPCH(JE)),JE=1,NEP),(ADJUSTR(MACCH(JM)),JM=1,NMC),                            &  
            (ADJUSTR(CDCH(JAD)),JAD=1,NACD(JW)),(ADJUSTR(KFCH(JF)),JF=1,NAF(JW)),(ADJUSTR(APCH(JA)),JA=1,NAL),(ADJUSTR(ANCH(JA)),JA=1,NAL),(ADJUSTR(ALCH(JA)),JA=1,NAL)    ! SW 10/20/15  
          END IF  
        ELSE  
          IF(SEDIMENT_CALC(JW))THEN  
            WRITE (TSR(J),'(f10.3,",",18(F10.3,","),*(A,","))') JDAY,DLT,ELWS(I),T1(K,I),U(K,I),QC(I),SRON(JW)*1.06,GAMMA(K,I),DEPTHB(KB(I),I),    &     ! SW 8/13/06  
            BI(KTWB(JW),I),SHADE(I),TVOLAVG,rn(i),rs(i),ranlw(jw),rb(i),re(i),rc(i),REAER(I)*86400.,(ADJUSTR(C2CH(JAC)),JAC=1,NAC),(ADJUSTR(EPCH(JE)),            &  ! CB 7/26/07
            JE=1,NEP),(ADJUSTR(MACCH(JM)),JM=1,NMC),SEDCH,SEDPCH,SEDNCH,SEDCCH,                   &  
            (ADJUSTR(CDCH(JAD)),JAD=1,NACD(JW)),(ADJUSTR(KFCH(JF)),JF=1,NAF(JW)),(ADJUSTR(APCH(JA)),JA=1,NAL),(ADJUSTR(ANCH(JA)),JA=1,NAL),(ADJUSTR(ALCH(JA)),JA=1,NAL)    ! SW 10/20/15  
          ELSE  
            WRITE (TSR(J),'(f10.3,",",18(F10.3,","),*(A,","))') JDAY,DLT,ELWS(I),T1(K,I),U(K,I),QC(I),SRON(JW)*1.06,GAMMA(K,I),DEPTHB(KB(I),I),      &      ! SW 8/13/06  
            BI(KTWB(JW),I),SHADE(I),TVOLAVG,rn(i),rs(i),ranlw(jw),rb(i),re(i),rc(i),REAER(I)*86400.,(ADJUSTR(C2CH(JAC)),JAC=1,NAC),(ADJUSTR(EPCH(JE)),            &  ! CB 7/26/07
            JE=1,NEP),(ADJUSTR(MACCH(JM)),JM=1,NMC),(ADJUSTR(CDCH(JAD)),JAD=1,NACD(JW)),    &  
            (ADJUSTR(KFCH(JF)),JF=1,NAF(JW)),(ADJUSTR(APCH(JA)),JA=1,NAL),(ADJUSTR(ANCH(JA)),JA=1,NAL),(ADJUSTR(ALCH(JA)),JA=1,NAL)    ! SW 10/20/15  
          END IF  
        END IF 
     ELSE      ! SW 4/4/2018
        IF (ICE_COMPUTATION) THEN  
          IF(SEDIMENT_CALC(JW))THEN  
            WRITE (TSR(J),'(f10.3,",",19(F10.3,","),*(F10.1,","))') JDAY,DLT,ELWS(I),-99.,-99.,QC(I),SRON(JW)*1.06,-99.,DEPTHB(KB(I),I),    &     ! SW 8/13/06  
            BI(KTWB(JW),I),SHADE(I),ICETH(I),TVOLAVG,rn(i),rs(i),ranlw(jw),rb(i),re(i),rc(i),REAER(I)*86400.,(-99.,JAC=1,NAC),                      &  ! CB 7/26/07
            (-99.,JE=1,NEP),(-99.,JM=1,NMC), SED(K,I),SEDP(K,I),SEDN(K,I),SEDC(K,I),               &                                !SEDCH,SEDPCH,SEDNCH,SEDCCH, &  
            (-99.,JAD=1,NACD(JW)),(-99.,JF=1,NAF(JW)),(-99.,JA=1,NAL),(-99.,JA=1,NAL),(-99.,JA=1,NAL)     ! SW 10/20/15 
          ELSE  
            WRITE (TSR(J),'(f10.3,",",19(F10.3,","),*(F10.1,","))') JDAY,DLT,ELWS(I),-99.,-99.,QC(I),SRON(JW)*1.06,-99.,DEPTHB(KB(I),I),      &     ! SW 8/13/06  
            BI(KTWB(JW),I),SHADE(I),ICETH(I),TVOLAVG,rn(i),rs(i),ranlw(jw),rb(i),re(i),rc(i),REAER(I)*86400.,(-99.,JAC=1,NAC),                      &  ! CB 7/26/07
            (-99.,JE=1,NEP),(-99.,JM=1,NMC),                            &  
            (-99.,JAD=1,NACD(JW)),(-99.,JF=1,NAF(JW)),(-99.,JA=1,NAL),(-99.,JA=1,NAL),(-99.,JA=1,NAL)     ! SW 10/20/15  
          END IF  
        ELSE  
          IF(SEDIMENT_CALC(JW))THEN  
            WRITE (TSR(J),'(f10.3,",",18(F10.3,","),*(F10.1,","))') JDAY,DLT,ELWS(I),-99.,-99.,QC(I),SRON(JW)*1.06,-99.,DEPTHB(KB(I),I),    &     ! SW 8/13/06  
            BI(KTWB(JW),I),SHADE(I),TVOLAVG,rn(i),rs(i),ranlw(jw),rb(i),re(i),rc(i),REAER(I)*86400.,(-99.,JAC=1,NAC),                      &  ! CB 7/26/07
            (-99.,JE=1,NEP),(-99.,JM=1,NMC), SED(K,I),SEDP(K,I),SEDN(K,I),SEDC(K,I),               &                                  !SEDCH,SEDPCH,SEDNCH,SEDCCH, &  
            (-99.,JAD=1,NACD(JW)),(-99.,JF=1,NAF(JW)),(-99.,JA=1,NAL),(-99.,JA=1,NAL),(-99.,JA=1,NAL)    ! SW 10/20/15 
          ELSE  
            WRITE (TSR(J),'(f10.3,",",18(F10.3,","),*(F10.1,","))') JDAY,DLT,ELWS(I),-99.,-99.,QC(I),SRON(JW)*1.06,-99.,DEPTHB(KB(I),I),      &     ! SW 8/13/06  
            BI(KTWB(JW),I),SHADE(I),TVOLAVG,rn(i),rs(i),ranlw(jw),rb(i),re(i),rc(i),REAER(I)*86400.,(-99.,JAC=1,NAC),                      &  ! CB 7/26/07
            (-99.,JE=1,NEP),(-99.,JM=1,NMC),                            &  
            (-99.,JAD=1,NACD(JW)),(-99.,JF=1,NAF(JW)),(-99.,JA=1,NAL),(-99.,JA=1,NAL),(-99.,JA=1,NAL)    ! SW 10/20/15  
          END IF  
        END IF        
     ENDIF   ! SW 4/4/2018
     
    IF(ALGAE_TOXIN)THEN
      IF(ATOX_DEBUG=='ON')THEN
          WRITE(ATOXIN_DEBUG_FN,'(F10.3,",",I4,",",I4,",",<NUMATOXINS>(E12.4,","),<NUMATOXINS>(E12.4,","),<NUMATOXINS>(E12.4,","),<NAL>(E12.4,","))')JDAY,K,I,(EX_TOXIN(K,I,JA),JA=1,numatoxins),(IN_TOXIN(K,I,JA),JA=1,numatoxins),(CTESS(K,I,JA),JA=1,numatoxins),(ALG(K,I,JA),JA=1,NAL)
      ENDIF
    ENDIF    
     
      END DO  
    END IF  
  END IF  
    
     ! SEDIMENT DIAGENESIS FREQUENCY OUTPUT      
          if(constituents.and.sediment_diagenesis.and.JDAY >= NXTSEDIAG)then 
          NXTSEDIAG=NXTSEDIAG+SEDIAGFREQ
          Call WriteCEMASedimentModelOutput
          Call WriteCEMASedimentFluxOutput
          end if           
   
  DO JW=1,NWB  
      
    !**** Inactive segments  
      
    JB       = BS(JW)  
    NBL(JW)  = 1  
    IBPR(JW) = 1  
    DO I=1,NISNP(JW)-1  
      IF (CUS(JB) > ISNP(I,JW)) THEN  
        BL(NBL(JW),JW) = I  
        NBL(JW)        = NBL(JW)+1  
        IBPR(JW)       = I+1  
      END IF  
      IF (ISNP(I+1,JW) > DS(JB)) JB = JB+1  
    END DO  
    NBL(JW) = NBL(JW)-1  
      
    !**** Snapshots  
      
    IF (SNAPSHOT(JW)) THEN  
      IF (JDAY >= NXTMSN(JW) .OR. JDAY >= SNPD(SNPDP(JW)+1,JW)) THEN  
        IF (JDAY >= SNPD(SNPDP(JW)+1,JW)) THEN  
          SNPDP(JW)  = SNPDP(JW)+1  
          NXTMSN(JW) = SNPD(SNPDP(JW),JW)  
        END IF  
        NXTMSN(JW) = NXTMSN(JW)+SNPF(SNPDP(JW),JW)  
        WRITE (SNP(JW),10490) W2VER,(TITLE(J),J=1,10)  
        WRITE (SNP(JW),10500) 'Time Parameters',MONTH,GDAY,YEAR,INT(JDAY),(JDAY-INT(JDAY))*24.0,INT(ELTMJD),                     &  
        (ELTMJD-INT(ELTMJD))*24.0,INT(DLTS1),KLOC,ILOC,INT(MINDLT),INT(JDMIN),(JDMIN-INT(JDMIN))*24.0,     &  
        KMIN,IMIN  
        IF (LIMITING_DLT(JW))  WRITE (SNP(JW),10510) KMIN,IMIN  
        WRITE (SNP(JW),10520)  INT(DLTAV),NIT,NV  
        WRITE (SNP(JW),10530) 'Meteorological Parameters'  
        WRITE (SNP(JW),10540)  TAIR(JW),DEG,TDEW(JW),DEG,PHI(JW),CLOUD(JW),ET(DS(1)),DEG,CSHE(DS(1)),SRON(JW),DEG  
        WRITE (SNP(JW),10550) 'Inflows','Upstream inflows'  
        DO JB=BS(JW),BE(JW)  
          IF (UP_FLOW(JB)) WRITE (SNP(JW),10560) JB,KTQIN(JB),KBQIN(JB),QIN(JB),TIN(JB),DEG  
        END DO  
        DO JB=BS(JW),BE(JW)  
          IF (DIST_TRIBS(JB)) THEN  
            WRITE (SNP(JW),10570)  
            WRITE (SNP(JW),10580) JB,QDTR(JB),TDTR(JB),DEG  
          END IF  
        END DO  
        IF (TRIBUTARIES) THEN  
          WRITE (SNP(JW),10590) (ITR(JT),          JT=1,JTT)  
          WRITE (SNP(JW),10600) (KTTR(JT),KBTR(JT),JT=1,JTT)  
          WRITE (SNP(JW),10610) (QTR(JT),          JT=1,JTT)  
          WRITE (SNP(JW),10620) (TTR(JT),          JT=1,JTT)  
        END IF  
        WRITE (SNP(JW),10630)  
        DO JB=BS(JW),BE(JW)  
          IF (DN_FLOW(JB)) THEN  
            WRITE (SNP(JW),10640)  JB,(QSTR(JS,JB),JS=1,JSS(JB))  
            WRITE (SNP(JW),10650)  QSUM(JB),(K,K=KTWB(JW),KB(DS(JB)))  
            WRITE (SNP(JW),10660) (QOUT(K,JB), K=KTWB(JW),KB(DS(JB)))  
            WRITE (SNP(JW),*)
            WRITE (SNP(JW),'(A)')'  LAYER   DEPTH(m)   T(C) DENSITY(kg/m3)  U(m/s)     Q(m3/s)'
            DO K=KTWB(JW),KB(DS(JB))
                WRITE (SNP(JW),'(I5,1X,F10.2,1X,F8.2,1X,F10.3,1X,F10.5,1X,F10.3)')K, DEPTHM(K,DS(JB)),T2(K,DS(JB)),RHO(K,DS(JB)),U(K,DS(JB)),QOUT(K,JB)
            ENDDO
          END IF  
        END DO  
        IF (WITHDRAWALS) THEN  
          DO JWD=1,JWW  
            WRITE (SNP(JW),10670) MAX(CUS(JBWD(JWD)),IWD(JWD)),QWD(JWD)  
            IF (QWD(JWD) /= 0.0) THEN  
              WRITE (SNP(JW),10680) (K,         K=KTW(JWD),KBW(JWD))  
              WRITE (SNP(JW),10690) (QSW(K,JWD),K=KTW(JWD),KBW(JWD))  
              
            WRITE (SNP(JW),*)
            WRITE (SNP(JW),'(A)')'  LAYER   DEPTH(m)   T(C) DENSITY(kg/m3)  Q(m3/s)'
            DO K=KTWB(JW),KB(IWD(JWD))
                WRITE (SNP(JW),'(I5,1X,F10.2,1X,F8.2,1X,F10.3,1X,F10.3)')K, DEPTHM(K,IWD(JWD)),T2(K,IWD(JWD)),RHO(K,IWD(JWD)),QSW(K,JWD)
            ENDDO             
              
            ELSE  
              WRITE (SNP(JW),10680)  
              WRITE (SNP(JW),10690) QWD(JWD)  
            END IF  
          END DO  
        END IF  
        IF (CONSTITUENTS) THEN  
          WRITE (SNP(JW),10700) 'Constituent Inflow Concentrations'  
          DO JB=BS(JW),BE(JW)  
            IF (UP_FLOW(JB) .AND. NACIN(JB) > 0)    WRITE (SNP(JW),10710) JB,(CNAME1(INCN(JC,JB))(1:18),CIN(INCN(JC,JB),JB),     &  
            CUNIT2(INCN(JC,JB)),JC=1,NACIN(JB))  
            IF (DIST_TRIBS(JB) .AND. NACDT(JB) > 0) WRITE (SNP(JW),10730) JB,(CNAME1(DTCN(JC,JB))(1:18),CDTR(DTCN(JC,JB),JB),    &  
            CUNIT2(DTCN(JC,JB)),JC=1,NACDT(JB))  
          END DO  
          DO JT=1,NTR  
            IF (NACTR(JT) > 0) WRITE (SNP(JW),10720) JT,(CNAME1(TRCN(JC,JT))(1:18),CTR(TRCN(JC,JT),JT),CUNIT2(TRCN(JC,JT)),      &  
            JC=1,NACTR(JT))  
          END DO  
        END IF  
        IF (EVAPORATION(JW) .OR. PRECIPITATION(JW)) WRITE (SNP(JW),10740)  
        IF (EVAPORATION(JW)) THEN  
          WRITE (SNP(JW),10750) (JB,EVBR(JB),JB=BS(JW),BE(JW))  ! SW 9/15/05
          WRITE (SNP(JW),10755) (JB,-VOLEV(JB),JB=BS(JW),BE(JW))  
        END IF  
        IF (PRECIPITATION(JW)) WRITE (SNP(JW),10760) (JB,PR(JB),JB=BS(JW),BE(JW))  
        IF (HEAD_BOUNDARY(JW)) THEN  
          WRITE (SNP(JW),10770)  
          DO JB=BS(JW),BE(JW)  
            IF (UH_EXTERNAL(JB)) WRITE (SNP(JW),10780) JB,ELUH(JB)  
            IF (DH_EXTERNAL(JB)) WRITE (SNP(JW),10790) JB,ELDH(JB)  
          END DO  
        END IF  
        IF (VOLUME_BALANCE(JW)) THEN  
          WRITE (SNP(JW),10800)  
          WRITE (SNP(JW),10810) JW,VOLSR(JW),VOLTR(JW),VOLTR(JW)-VOLSR(JW),DLVR(JW)  
          DO JB=BS(JW),BE(JW)  
            IF (VOLSBR(JB) /= 0.0) DLVBR = (VOLTBR(JB)-VOLSBR(JB))/VOLSBR(JB)  
            WRITE (SNP(JW),10820) JB,VOLSBR(JB),VOLTBR(JB),VOLTBR(JB)-VOLSBR(JB),DLVBR*100.0  
          END DO  
        END IF  
        IF (ENERGY_BALANCE(JW)) THEN  
          WRITE (SNP(JW),10830)  
          IF (ESR(JW) /= 0.0) DLE = (ESR(JW)-ETR(JW))/ESR(JW)  
          WRITE (SNP(JW),10840) JW,ESR(JW)*4.184E3,ETR(JW)*4.184E3,(ESR(JW)-ETR(JW))*4.184E3,DLE*100.0  
          DO JB=BS(JW),BE(JW)  
            WRITE (SNP(JW),10870) JB  
            IF (ESBR(JB) /= 0.0) DLE = (ESBR(JB)-ETBR(JB))/ESBR(JB)  
            WRITE (SNP(JW),10850) ESBR(JB)*4.184E3,ETBR(JB)*4.1843E3,(ESBR(JB)-ETBR(JB))*4.1843E3,DLE*100.0  
          END DO  
        END IF  
        IF (MASS_BALANCE(JW)) THEN  
          WRITE (SNP(JW),10860)  
          DO JB=BS(JW),BE(JW)  
            WRITE (SNP(JW),10870) JB  
            DO JC=1,NAC  
              IF (CMBRS(CN(JC),JB) /= 0.0) DLMR = (CMBRT(CN(JC),JB)-CMBRS(CN(JC),JB))/(CMBRS(CN(JC),JB)+NONZERO)*100.0  
              WRITE (SNP(JW),10880) CNAME1(CN(JC)),CMBRS(CN(JC),JB),CUNIT1(CN(JC)),CMBRT(CN(JC),JB),CUNIT1(CN(JC)),              &  
                                    (CMBRT(CN(JC),JB)-CMBRS(CN(JC),JB)),CUNIT1(CN(JC)),DLMR  
            END DO  
              
            DO M=1,NMC  
              IF(MACROPHYTE_CALC(JW,M))THEN  
                IF (MACMBRS(JB,M).NE.0.0) THEN  
                  DLMR = (MACMBRT(JB,M)-MACMBRS(JB,M))/(MACMBRS(JB,M)+NONZERO)  
                END IF  
                WRITE (SNP(JW),3312) M,MACMBRS(JB,M),MACMBRT(JB,M),(MACMBRT(JB,M)-MACMBRS(JB,M)),DLMR*100.0  
3312            FORMAT(5X,'Macrophyte spec ',i2,/7X,'Spatially integrated mass [MACMBRS] = ',1PE15.8E2,1X,'g ',/7X,   &  
                    'Temporally integrated mass [MACMBRT] = ',1PE15.8E2,1X,'g ',/7X,'Mass error                         = ',  &  
                    1PE15.8E2,1X,'g ',/7X,'Percent error                      = ',1PE15.8E2,' %')  
              END IF  
            END DO  
              
          END DO  
        END IF  
        WRITE (SNP(JW),10890) 'Geometry',KTWB(JW),ELKT(JW)  
        WRITE (SNP(JW),10900) (JB,CUS(JB),JB=BS(JW),BE(JW))  
        CALL OUTPUT (JDAY,IBPR(JW),NISNP(JW),KBR(JW),ISNP,BL(1,JW),NBL(JW))  

      END IF  
    END IF  
      
    !**** Vertical profiles  
      
    IF (PROFILE(JW)) THEN  
      IF (JDAY >= NXTMPR(JW) .OR. JDAY >= PRFD(PRFDP(JW)+1,JW)) THEN  
        IF (JDAY >= PRFD(PRFDP(JW)+1,JW)) THEN  
          PRFDP(JW)  = PRFDP(JW)+1  
          NXTMPR(JW) = PRFD(PRFDP(JW),JW)  
        END IF  
        NXTMPR(JW) = NXTMPR(JW)+PRFF(PRFDP(JW),JW)  
        NSPRF(JW)  = NSPRF(JW)+1  
        if(iprf(1,1) /= -1)then     ! SW 4/1/2016
        WRITE (PRF(JW),'(F8.3,1X,A3,I3,A,2I4,F8.4,I8)')JDAY,ADJUSTL(MONTH),GDAY,', ',YEAR,KTWB(JW),SNGL(Z(DS(BS(JW)))),NSPRF(JW)  
        DO JP=1,NIPRF(JW)  
          NRS = KB(IPRF(JP,JW))-KTWB(JW)+1  
          WRITE (PRF(JW),'(A8,I4/(8F10.2))') 'TEMP    ',NRS,(T2(K,IPRF(JP,JW)),K=KTWB(JW),KB(IPRF(JP,JW)))  
        END DO  
        DO JC=1,NAC  
          IF (PRINT_CONST(CN(JC),JW)) THEN  
            DO JP=1,NIPRF(JW)  
              NRS = KB(IPRF(JP,JW))-KTWB(JW)+1  
              WRITE (PRF(JW),'(A8,I4/(8(E13.6,2X)))') ADJUSTL(CNAME2(CN(JC))),NRS,(C2(K,IPRF(JP,JW),CN(JC))*CMULT(CN(JC)),       &     ! CB 1/24/05  
              K=KTWB(JW),KB(IPRF(JP,JW)))  
            END DO  
          END IF  
        END DO  
        IF (CONSTITUENTS) THEN  
          DO JD=1,NACD(JW)  
            DO JP=1,NIPRF(JW)  
              NRS = KB(IPRF(JP,JW))-KTWB(JW)+1  
              WRITE (PRF(JW),'(A8,I4/(8(E13.6,2X)))') ADJUSTL(CDNAME2(CDN(JD,JW))),NRS,(CD(K,IPRF(JP,JW),CDN(JD,JW))             &      ! CB 1/24/05  
              *CDMULT(CDN(JD,JW)),K=KTWB(JW), KB(IPRF(JP,JW)))  
            END DO  
          END DO  
        END IF
        
        elseif(jw == 1)then   ! write out individual files on these days
             WRITE (SEGNUM,'(F8.2)')JDAY            !   '(I0)'    INT(JDAY)
             SEGNUM = ADJUSTL(SEGNUM)
             L      = LEN_TRIM(SEGNUM)
             OPEN  (NUNIT,FILE='ProfLongJD'//SEGNUM(1:L)//'.csv',STATUS='UNKNOWN')
             IF(CONSTITUENTS)THEN
             WRITE(NUNIT,'(*(A,","))')'Seg#','ElevWaterSurf(m)','Q(m3/s)','SurfaceTemp(oC)','Depth(m)','Width(m)','VolWeighTemp(oC)' ,(CNAME2(CN(JC)),JC=1,NAC),(CDNAME2(CDN(JD,JW)),JD=1,NACD(JW))
             ELSE
             WRITE(NUNIT,'(*(A,","))')'Seg#','ElevWaterSurf(m)','Q(m3/s)','SurfaceTemp(oC)','Depth(m)','Width(m)','VolWeighTemp(oC)'   
             ENDIF
             
                 DO JJ=1,NWB
                     K=KTWB(JJ)
                     DO JB=BS(JJ),BE(JJ)
                         DO I=CUS(JB),DS(JB)
                             
                                     ! TEMP VOL WEIGHTED AVERAGE  SW 8/30/2018
                                        TVOLAVG=0.0     ! SW 8/30/2018
                                        VOLTOT=0.0
                                        DO KK=KTWB(JW),KB(I)   
                                            VOLTOT=VOLTOT+VOL(KK,I)
                                            TVOLAVG=TVOLAVG+T1(KK,I)*VOL(KK,I) 
                                        ENDDO
                                        IF(KB(I)>=KTWB(JW))tvolavg=tvolavg/voltot
    
                             IF(CONSTITUENTS)THEN
                             WRITE(NUNIT,'(I5,",",*(F12.4,","))')I,ELWS(I),QC(I),T2(K,I),DEPTHB(KB(I),I),B(KTI(I),I),TVOLAVG,(C2(K,I,CN(JAC))*CMULT(CN(JAC)),JAC=1,NAC),(CD(K,I,CDN(JD,JW))*CDMULT(CDN(JD,JW)),JD=1,NACD(JW))
                             ELSE
                            WRITE(NUNIT,'(I5,",",*(F12.3,","))')I,ELWS(I),QC(I),T2(K,I),DEPTHB(KB(I),I),B(KTI(I),I),TVOLAVG
                            ENDIF
                         ENDDO
                     ENDDO
                 ENDDO
                 CLOSE(NUNIT)
          endif     ! end of iprf /= -1
   
      END IF  
    END IF  
      
    !**** Spreadsheet  
      
    IF (SPREADSHEET(JW)) THEN  
      IF (JDAY >= NXTMSP(JW) .OR. JDAY >= SPRD(SPRDP(JW)+1,JW)) THEN  
        IF (JDAY >= SPRD(SPRDP(JW)+1,JW)) THEN  
          SPRDP(JW)  = SPRDP(JW)+1  
          NXTMSP(JW) = SPRD(SPRDP(JW),JW)  
        END IF  
        CONV1      = BLANK1  
        NXTMSP(JW) = NXTMSP(JW)+SPRF(SPRDP(JW),JW)  
        !DO J=1,NISPR(JW)  
        !  KBMAX(JW) = MAX(KB(ISPR(J,JW)),KBMAX(JW))  
        IF(SPRC(JW)   == '     ONV')THEN    ! SW 8/31/2018
        DO J=1,NISPR(JW)  
          !KBMAX(JW) = MAX(KB(ISPR(J,JW)),KBMAX(JW))
            K=KTWB(JW)
              ! TEMP VOL WEIGHTED AVERAGE  SW 8/30/2018
                                        TVOLAVG=0.0     ! SW 8/30/2018
                                        VOLTOT=0.0
                                        DO KK=KTWB(JW),KB(ISPR(J,JW))      ! cb 11/13/18
                                            VOLTOT=VOLTOT+VOL(KK,ISPR(J,JW))
                                            TVOLAVG=TVOLAVG+T1(KK,ISPR(J,JW))*VOL(KK,ISPR(J,JW)) 
                                        ENDDO
                                        IF(KB(ISPR(J,JW))>=KTWB(JW))THEN
                                            tvolavg=tvolavg/voltot
                                        ELSE
                                            TVOLAVG=-99.0
                                        ENDIF
                                        WRITE (CONV1(K,ISPR(J,JW)),'(F12.3)') TVOLAVG
        ENDDO          
        WRITE (SPRV(JW),'(A,",",F10.3,",",*(A,","))') 'Temperature',JDAY,(CONV1(K,ISPR(J,JW)),J=1,NISPR(JW))  
        
        DO JC=1,NAC  
          IF (PRINT_CONST(CN(JC),JW)) THEN  
            DO J=1,NISPR(JW)  
                                       TVOLAVG=0.0     ! SW 8/30/2018
                                        VOLTOT=0.0
                                        DO KK=KTWB(JW),KB(ISPR(J,JW))   ! cb 11/13/18
                                            VOLTOT=VOLTOT+VOL(KK,ISPR(J,JW))
                                            TVOLAVG=TVOLAVG+ C2(KK,ISPR(J,JW),CN(JC))*CMULT(CN(JC))*VOL(KK,ISPR(J,JW)) 
                                        ENDDO
                                        IF(KB(ISPR(J,JW))>=KTWB(JW))THEN
                                            tvolavg=tvolavg/voltot
                                        ELSE
                                            TVOLAVG=-99.0
                                        ENDIF
                WRITE (CONV1(K,ISPR(J,JW)),'(F12.4)') TVOLAVG                                                ! SW 8/31/18  
            END DO   
              WRITE (SPRV(JW),'(A38,",",F10.3,",",*(A,","))') CNAME3(CN(JC)),JDAY,(CONV1(K,ISPR(J,JW)),J=1,NISPR(JW))   
          END IF  
        END DO  
        IF (CONSTITUENTS) THEN  
          DO JD=1,NACD(JW)  
            IF (PRINT_DERIVED(CDN(JD,JW),JW)) THEN  
              DO J=1,NISPR(JW)  
                  TVOLAVG=0.0     ! SW 8/30/2018
                                        VOLTOT=0.0
                                        DO KK=KTWB(JW),KB(ISPR(J,JW))   ! cb 11/13/18
                                            VOLTOT=VOLTOT+VOL(KK,ISPR(J,JW))
                                            TVOLAVG=TVOLAVG+ CD(KK,ISPR(J,JW),CDN(JD,JW))*CDMULT(CDN(JD,JW))*VOL(KK,ISPR(J,JW)) 
                                        ENDDO
                                        IF(KB(ISPR(J,JW))>=KTWB(JW))THEN
                                            tvolavg=tvolavg/voltot
                                        ELSE
                                            TVOLAVG=-99.0
                                        ENDIF
                  WRITE (CONV1(K,ISPR(J,JW)),'(F12.4)') TVOLAVG                                     ! SW 8/31/18  
              END DO  
                WRITE (SPRV(JW),'(A38,",",F10.3,",",*(A,","))') CDNAME3(CDN(JD,JW)),JDAY,(CONV1(K,ISPR(J,JW)),J=1,NISPR(JW))  
            END IF  
            END DO  
            END IF    
          ENDIF         ! END SECTION ONV
        DO J=1,NISPR(JW)  
          KBMAX(JW) = MAX(KB(ISPR(J,JW)),KBMAX(JW))  
          DO K=KTWB(JW),KB(ISPR(J,JW))  
            WRITE (CONV1(K,J),'(F12.3)') T2(K,ISPR(J,JW))  
          END DO  

        END DO  
        DO K=KTWB(JW),KBMAX(JW)  
          WRITE (SPR(JW),'(A,",",2(F10.3,","),*(F10.3,",",A,","))') 'Temperature',JDAY,DEPTHM(K,DS(BS(JW))),                      &  
          (ELWS(ISPR(J,JW))-DEPTHM(K,ISPR(J,JW)),CONV1(K,J),J=1,NISPR(JW))  
        END DO  
        DO JC=1,NAC  
          IF (PRINT_CONST(CN(JC),JW)) THEN  
            DO J=1,NISPR(JW)  
              DO K=KTWB(JW),KB(ISPR(J,JW))  
                WRITE (CONV1(K,J),'(F12.4)') C2(K,ISPR(J,JW),CN(JC))*CMULT(CN(JC))                                                ! SW 8/13/06  
              END DO  
            END DO  
            DO K=KTWB(JW),KBMAX(JW)  
              WRITE (SPR(JW),'(A38,",",2(F10.3,","),*(F10.3,",",A,","))') CNAME3(CN(JC)),JDAY,DEPTHM(K,DS(BS(JW))),                            &  
              (ELWS(ISPR(J,JW))-DEPTHM(K,ISPR(J,JW)),CONV1(K,J),J=1,NISPR(JW))  
            END DO  
          END IF  
        END DO  
        IF (CONSTITUENTS) THEN  
          DO JD=1,NACD(JW)  
            IF (PRINT_DERIVED(CDN(JD,JW),JW)) THEN  
              DO J=1,NISPR(JW)  
                DO K=KTWB(JW),KB(ISPR(J,JW))  
                  WRITE (CONV1(K,J),'(F12.4)') CD(K,ISPR(J,JW),CDN(JD,JW))*CDMULT(CDN(JD,JW))                                     ! SW 8/13/06  
                END DO  
              END DO  
              DO K=KTWB(JW),KBMAX(JW)  
                WRITE (SPR(JW),'(A38,",",2(F10.3,","),*(F10.3,",",A,","))') CDNAME3(CDN(JD,JW)),JDAY,DEPTHM(K,DS(BS(JW))),                     &  
                (ELWS(ISPR(J,JW))-DEPTHM(K,ISPR(J,JW)),CONV1(K,J),J=1,NISPR(JW))  
              END DO  
            END IF  
          END DO  
        END IF 
      END IF  
    END IF  
      
    !**** Contours  
      
    IF (CONTOUR(JW)) THEN  
      IF (JDAY >= NXTMCP(JW) .OR. JDAY >= CPLD(CPLDP(JW)+1,JW)) THEN  
        IF (JDAY >= CPLD(CPLDP(JW)+1,JW)) THEN  
          CPLDP(JW)  = CPLDP(JW)+1  
          NXTMCP(JW) = CPLD(CPLDP(JW),JW)  
        END IF  
        NXTMCP(JW) = NXTMCP(JW)+CPLF(CPLDP(JW),JW)  
          
        IF(TECPLOT(JW) /= '      ON')THEN  
          WRITE (CPL(JW),'(A,F12.4,5X,A9,5X,I2,5X,I4)') 'New date ',JDAY,MONTH,GDAY,YEAR  
          WRITE (CPL(JW),'(9(I8,2X))')                   KTWB(JW)  
          WRITE (CPL(JW),'(9(E13.6,2X))')               (QTR(JT),JT=1,NTR)  
          WRITE (CPL(JW),'(9(E13.6,2X))')               (TTR(JT),JT=1,NTR)  
          DO JT=1,NTR  
            DO JAC=1,NACTR(JT)  
              IF (PRINT_CONST(TRCN(JAC,JT),JW)) WRITE (CPL(JW),'(9(E13.6,2X))') CTR(TRCN(JAC,JT),JT)  
            END DO  
          END DO  
          DO JB=BS(JW),BE(JW)  
            IF(BR_INACTIVE(JB))CYCLE  ! SW 7/1/2019
            WRITE (CPL(JW),'(9(I8,2X))')             CUS(JB)  
            WRITE (CPL(JW),'(9(E13.6,2X))')          QIN(JB), QSUM(JB)  
            DO I=CUS(JB),DS(JB)  
              WRITE (CPL(JW),'(A38/(9(E13.6,2X)))') 'BHR', (BHR1(K,I),K=KTWB(JW)+1,KB(I))  
            END DO  
            DO I=CUS(JB),DS(JB)  
              WRITE (CPL(JW),'(A38/(9(E13.6,2X)))') 'U',   (U(K,I),   K=KTWB(JW),KB(I))  
            END DO  
            WRITE (CPL(JW),'(A38/(9(E13.6,2X)))')   'QC',  (QC(I),    I=CUS(JB),DS(JB))  
            WRITE (CPL(JW),'(A38/(9(E13.6,2X)))')   'Z',   (Z(I),     I=CUS(JB),DS(JB))  
            WRITE (CPL(JW),'(A38/(9(I8,2X)))')   'KTI',   (kti(I),     I=CUS(JB),DS(JB))  ! v3.5  
            DO I=CUS(JB),DS(JB)  
              WRITE (CPL(JW),'(A38/(9(E13.6,2X)))') 'Temperature',(T2(K,I),K=KTWB(JW),KB(I))  
            END DO  
            DO JC=1,NAC  
              IF (PRINT_CONST(CN(JC),JW)) THEN  
                DO I=CUS(JB),DS(JB)  
                  WRITE (CPL(JW),'(A38/(9(E13.6,2X)))') CNAME(CN(JC)),(C2(K,I,CN(JC))*CMULT(CN(JC)),K=KTWB(JW),KB(I))  
                END DO  
              END IF  
            END DO  
            DO JE=1,NEP  
              DO I=CUS(JB),DS(JB)  
                IF (PRINT_EPIPHYTON(JW,JE)) WRITE (CPL(JW),'(A38/(9(E13.6,2X)))') 'Epiphyton',(EPD(K,I,JE),K=KTWB(JW),KB(I))  
              END DO  
            END DO  
            IF(PRINT_SEDIMENT(JW))THEN  
              DO I=CUS(JB),DS(JB)  
                WRITE (CPL(Jw),'(A38/(9(E13.6,2X)))')'Sediment',(seD(K,I),K=KTWB(JW),KB(I))  
              END DO  
              DO I=CUS(JB),DS(JB)  
                WRITE (CPL(Jw),'(A38/(9(E13.6,2X)))')'Sediment P',(seDp(K,I),K=KTWB(JW),KB(I))  
              END DO  
              DO I=CUS(JB),DS(JB)  
                WRITE (CPL(Jw),'(A38/(9(E13.6,2X)))')'Sediment N',(seDn(K,I),K=KTWB(JW),KB(I))  
              END DO  
              DO I=CUS(JB),DS(JB)  
                WRITE (CPL(Jw),'(A38/(9(E13.6,2X)))')'Sediment C',(seDc(K,I),K=KTWB(JW),KB(I))  
              END DO  
            END IF  
            DO M=1,NMC  
              IF (PRINT_MACROPHYTE(JW,M)) THEN  
                DO I=CUS(JB),DS(JB)  
                  WRITE (CPL(Jw),'(A38/(9(E13.6,2X)))')'Macrophytes',((macrc(j,K,I,m),j=kti(i),kb(i)), K=KTwb(Jw),KB(I))  
                END DO  
              END IF  
            END DO  
            IF (CONSTITUENTS) THEN  
              DO JD=1,NACD(JW)  
                IF (PRINT_DERIVED(CDN(JD,JW),JW)) THEN  
                  DO I=CUS(JB),DS(JB)  
                     WRITE (CPL(JW),'(A38/(9(F10.3,2X)))') CDNAME(CDN(JD,JW)),(CD(K,I,CDN(JD,JW))*CDMULT(CDN(JD,JW)),    K=KTWB(JW),KB(I))  ! cb 6/28/13
                  end do
                  !WRITE (CPL(JW),'(A38/(9(F10.3,2X)))') CDNAME(CDN(JD,JW)),((CD(K,I,CDN(JD,JW))*CDMULT(CDN(JD,JW)),             &        ! SW 8/12/06  
                  !K=KTWB(JW),KB(I)),I=CUS(JB),DS(JB))  ! CB 1/03/05
                END IF  
              END DO  
            END IF  
          END DO  
        ELSE  
          !         ICPL=ICPL+1  
          ICPL(JW)=ICPL(JW)+1  ! cb 1/26/09
          ITOT=0  
          !         do jb=1,nbr  
          DO JB=BS(JW),BE(JW)  ! cb 1/26/09
            IF(BR_INACTIVE(JB))CYCLE  ! SW 7/1/2019
            IF(BR_NOTECPLOT(JB))CYCLE  ! SW 8/27/2019
            ITOT=ITOT+DS(JB)-CUS(JB)+2  
          ENDDO  
          WRITE (CPL(JW),9864)JDAY,KMX-KTWB(JW)+2,ITOT  
9864      FORMAT('ZONE T="',f9.3,'"',' I=',I3,' J=',I3,' F=POINT')  
          DO JB=BS(JW),BE(JW) 
            IF(BR_INACTIVE(JB))CYCLE  ! SW 7/1/2019
            IF(BR_NOTECPLOT(JB))CYCLE  ! SW 8/27/2019
            DO I=CUS(JB), DS(JB)+1  
              K=KTWB(JW)  ! PRINT AN EXTRA LINE FOR THE SURFACE
              IF(I /= DS(JB)+1)THEN  
                IF(HABTATC  == '      ON')THEN
                  WRITE (CPL(JW),9999) X1(I),ELWS(I),U(K,I),-W(K,I),T1(K,I),RHO(K,I),HAB(K,I),(C2(K,I,CN(JC)),JC=1,NAC),(CD(K,I,CDN(JD,JW)),JD=1,NACD(JW))   ! SW 1/17/17
                ELSE
                  WRITE (CPL(JW),9999) X1(I),ELWS(I),U(K,I),-W(K,I),T1(K,I),RHO(K,I),(C2(K,I,CN(JC)),JC=1,NAC),(CD(K,I,CDN(JD,JW)),JD=1,NACD(JW))            ! SW 1/17/17
                ENDIF
                
              ELSE  
                XDUM=-99.0  
                IF(HABTATC  == '      ON')THEN                                                 ! SW 7/15/14
                WRITE (CPL(JW),9999) X1(I),ELWS(I),XDUM,XDUM,XDUM,XDUM,XDUM,(XDUM, JJ=1,NAC),(XDUM, JJ=1,NACD(JW))  
                ELSE                                                                           ! SW 7/15/14
                WRITE (CPL(JW),9999) X1(I),ELWS(I),XDUM,XDUM,XDUM,XDUM,(XDUM, JJ=1,NAC),(XDUM, JJ=1,NACD(JW))          ! SW 1/17/17 7/15/14
                ENDIF                                                                          ! SW 7/15/14
              ENDIF  
              DO K=KTWB(JW),KMX-1  
                IF(I /= DS(JB)+1 .AND. K <= KB(I))THEN  
                IF(HABTATC  == '      ON')THEN
                  WRITE (CPL(JW),9999) X1(I),ELWS(I)-DEPTHM(K,I),U(K,I),-W(K,I),T1(K,I),RHO(K,I),HAB(K,I),(C2(K,I,CN(JC)),JC=1,NAC),(CD(K,I,CDN(JD,JW)),JD=1,NACD(JW))   ! SW 1/17/17 
                ELSE
                  WRITE (CPL(JW),9999) X1(I),ELWS(I)-DEPTHM(K,I),U(K,I),-W(K,I),T1(K,I),RHO(K,I),(C2(K,I,CN(JC)),JC=1,NAC),(CD(K,I,CDN(JD,JW)),JD=1,NACD(JW))   ! SW 1/17/17 
                ENDIF
                
                  IF(K == KB(I))THEN  
                      IF(HABTATC  == '      ON')THEN
                    WRITE (CPL(JW),9999) X1(I),ELWS(I)-DEPTHB(K,I),U(K,I), -W(K,I),T1(K,I),RHO(K,I),HAB(K,I),(C2(K,I,CN(JC)),JC=1,NAC),(CD(K,I,CDN(JD,JW)),JD=1,NACD(JW))   ! SW 1/17/17   
                      ELSE
                    WRITE (CPL(JW),9999) X1(I),ELWS(I)-DEPTHB(K,I),U(K,I), -W(K,I),T1(K,I),RHO(K,I),(C2(K,I,CN(JC)),JC=1,NAC),(CD(K,I,CDN(JD,JW)),JD=1,NACD(JW))   ! SW 1/17/17   
                      ENDIF
                      
                  ENDIF  
                ELSE  
                  XDUM=-99.0  
                  IF(HABTATC  == '      ON')THEN
                  WRITE (CPL(JW),9999) X1(I),ELWS(I-1)-DEPTHM(K,I-1),XDUM,XDUM,XDUM,XDUM,XDUM,(XDUM, JJ=1,NAC),(XDUM, JJ=1,NACD(JW)) 
                  ELSE
                  WRITE (CPL(JW),9999) X1(I),ELWS(I-1)-DEPTHM(K,I-1),XDUM,XDUM,XDUM,XDUM,(XDUM, JJ=1,NAC),(XDUM, JJ=1,NACD(JW))      ! SW 7/15/14
                  ENDIF
                  
                  IF(K == KB(I))THEN  
                    IF(HABTATC  == '      ON')THEN
                    WRITE (CPL(JW),9999) X1(I),ELWS(I-1)-DEPTHB(K,I-1),XDUM,XDUM,XDUM,XDUM,XDUM,(XDUM, JJ=1,NAC),(XDUM, JJ=1,NACD(JW))    
                    ELSE
                    WRITE (CPL(JW),9999) X1(I),ELWS(I-1)-DEPTHB(K,I-1),XDUM,XDUM,XDUM,XDUM,(XDUM, JJ=1,NAC),(XDUM, JJ=1,NACD(JW))         ! SW 7/15/14
                    ENDIF
                    
                  ENDIF  
                ENDIF  
              END DO  
            END DO  
          END DO  
          WRITE(CPL(JW),9899)ICPL(JW),IMON,GDAY,YEAR  ! cb 1/26/09
9899      FORMAT('TEXT X=0.75, y=0.85, H=2.8,ZN=',i4,',',' C=BLACK,','T= "',i2,'/',i2,'/',i4,'"')  
          WRITE(CPL(JW),9863)ICPL(JW),JDAY  ! cb 1/26/09
9863      FORMAT('TEXT X=0.75, y=0.90, H=2.8,ZN=',i4,',',' C=BLACK,','T= "Julian Day ',f9.3,'"')  
9999      FORMAT (200(e13.6,1x))  
        ENDIF  
      END IF  
    END IF     
    
    !**** Fluxes   KF is the instantaneous flux in g/m3/s, KFS is the summed flux in g eventually divided by elapsed time between calls to FLUX output and converted to kg below  
    
    
    IF (FLUX(JW)) THEN  
      IF (JDAY >= NXTMFL(JW) .OR. JDAY >= FLXD(FLXDP(JW)+1,JW)) THEN  
        IF (JDAY >= FLXD(FLXDP(JW)+1,JW)) THEN  
          FLXDP(JW)  = FLXDP(JW)+1  
          NXTMFL(JW) = FLXD(FLXDP(JW),JW)  
        END IF          
        
        NLINES=0    ! SW 3/8/16
        NXTMFL(JW) = NXTMFL(JW)+FLXF(FLXDP(JW),JW)  
        CONV       = BLANK  
        DO JAF=1,NAF(JW)  
         IF(ELTMF(JW) > 0.0)THEN
          DO JB=BS(JW),BE(JW)  
            DO I=CUS(JB),DS(JB)  
              DO K=KTWB(JW),KB(I)  
                KFS(K,I,KFCN(JAF,JW)) = DAY*KFS(K,I,KFCN(JAF,JW))/(1000.*ELTMF(JW))  ! KFS IN G, 86400 S/D * G /ELAPSED TIME IN S/1000 G/KG == KG/D
                KFJW(JW,KFCN(JAF,JW)) = KFJW(JW,KFCN(JAF,JW))+KFS(K,I,KFCN(JAF,JW))  ! SUM UP FOR ENTIRE WATERBODY  
              END DO  
            END DO  
          END DO  
        ENDIF
          DO I=1,NISNP(JW)  
            DO K=KTWB(JW),KB(ISNP(I,JW))  
              WRITE (CONV(K,I),'(E10.3)') KFS(K,ISNP(I,JW),KFCN(JAF,JW))      ! KG/D
            END DO  
          END DO  
          IF (NEW_PAGE) THEN  
            WRITE (FLX(JW),'(/(A72))') (TITLE(J),J=1,11)  
            NLINES   =  KMX-KTWB(JW)+14  
            NEW_PAGE = .FALSE.  
          END IF  
          NLINES   = NLINES+KMX-KTWB(JW)+11  
          NEW_PAGE = NLINES > 72  
          WRITE (FLX(JW),'(/A,F10.3,X,3(A,I0),A,F0.2,A/)') 'New date ',JDAY,MONTH//' ',GDAY,', ',YEAR,'   Julian Date = ',       &  
          INT(JDAY),' days ',(JDAY-INT(JDAY))*24.0,                            &  
          ' hours           '//KFNAME(KFCN(JAF,JW))  
          WRITE (FLX(JW),'(3X,2000I10)')                  (ISNP(I,JW),I=1,NISNP(JW))  
          DO K=KTWB(JW),KBR(JW)  
            WRITE (FLX(JW),'(1X,I2,200A)') K,(CONV(K,I),I=1,NISNP(JW))  
          END DO  
        END DO  
        WRITE(FLX2(JW),'(F10.3,",",f8.3,",",*(E12.4,","))')JDAY,ELTMF(JW)/DAY,(KFJW(JW,KFCN(K,JW)),K=1,NAF(JW))  
        ELTMF(JW)                    = 0.0  
        KF(:,CUS(BS(JW)):DS(BE(JW)),KFCN(1:NAF(JW),JW))  = 0.0  
        KFS(:,CUS(BS(JW)):DS(BE(JW)),KFCN(1:NAF(JW),JW)) = 0.0  
        KFJW(JW,KFCN(1:NAF(JW),JW))  = 0.0  
      END IF  
    END IF  
      
  END DO  
    
  !** Downstream flow, temperature, and constituent files  
    
  IF (DOWNSTREAM_OUTFLOW) THEN  
    IF (JDAY >= NXTMWD .OR. JDAY >= WDOD(WDODP+1)) THEN  
      IF (JDAY >= WDOD(WDODP+1)) THEN  
        WDODP  = WDODP+1  
        NXTMWD = WDOD(WDODP)  
      END IF  
      IF(WDOC  =='      ON')THEN      
      NXTMWD = NXTMWD+WDOF(WDODP)  
      ELSEIF(WDOC  =='     ONS')THEN
      NXTMWD_SEC = NXTMWD_SEC+WDOF(WDODP)    ! cb 4/6/18 frequency test seconds
      NXTMWD = NXTMWD_SEC/86400.0
      ELSEIF(WDOC    =='     ONH')THEN
      NXTMWD_SEC = NXTMWD_SEC+WDOF(WDODP)*3600.    ! cb 4/6/18 frequency test hourly
      NXTMWD = NXTMWD_SEC/86400.0
      ENDIF

      JFILE=0  

      DO J=1,NIWDO  
        QWDO(J)    = 0.0  
        TWDO(J)    = 0.0  
        CWDO(:,J)  = 0.0  
        CDWDO(:,J) = 0.0  
        CDTOT      = 0.0  
        NUMOUTLETS = 0  
        JWD=0                      ! SW 5/17/13
          
        DO JW=1,NWB  
          DO JB=BS(JW),BE(JW)  
            IF (DS(JB) == IWDO(J)) THEN  
                
              DO JS=1,JSS(JB)  !NSTR(JB)
                NUMOUTLETS=NUMOUTLETS+1  
                IF(QSTR(JS,JB)==0.0)THEN
                    TAVG(JS,JB)=-99.0
                    CAVG(JS,JB,:)=-99.0
                    CDAVG(JS,JB,:)=-99.0
                ENDIF
                QOUTLET(NUMOUTLETS)=QSTR(JS,JB)  
                TOUTLET(NUMOUTLETS)=TAVG(JS,JB)  
              ENDDO  
                
              ! cb 1/16/13 removed old code
                
              ! OUTPUT INDIVIDUAL FILES  
              DO JS=1,NSTR(JB)  
                JFILE=JFILE+1  
                qwdo(j)           = qwdo(j)            +qstr(js,jb)                 ! cb 1/16/13  
                twdo(j)           = twdo(j)            +qstr(js,jb)*tavg(js,jb)
                WRITE(WDO2(JFILE,1),'(F10.3,",",F10.4)')JDAY,QSTR(JS,JB)  
                WRITE(WDO2(JFILE,2),'(F10.3,",",F8.2)')JDAY,TAVG(JS,JB)  
                IF (CONSTITUENTS)then
                  cwdo(cn(1:nac),j)= cwdo(cn(1:nac),j)+qstr(js,jb)*cavg(js,jb,cn(1:nac))                 ! cb 1/16/13
                  WRITE (WDO2(JFILE,3),'(F10.3,",",*(F10.4,","))') JDAY,(CAVG(JS,JB,CN(JC)),     JC=1,NAC)  
                end if
                IF (DERIVED_CALC)then
                  cdwdo(cdn(1:nacd(jw),jw),j) = cdwdo(cdn(1:nacd(jw),jw),j)+qstr(js,jb)*cdavg(js,jb,cdn(1:nacd(jw),jw))                 ! cb 1/16/13 
                  WRITE (WDO2(JFILE,4),'(F10.3,",",*(F10.4,","))') JDAY,(CDAVG(JS,JB,CDN(JD,JW)),JD=1,NACD(JW))  
                end if
              ENDDO  
              JSSS(JB)=NSTR(JB)  
            END IF  
          END DO  
        END DO  
          
        ! OUTPUT INDIVIDUAL FILES  
        ! Order: spillways NSP, pumps NPU, gates NGT, pipes NPI  
   !     JWD=0                        ! sw 5/17/13
        do jw=1,nwb                   ! cb 1/16/13
          if (iwdo(j) >= us(bs(jw)) .and. iwdo(j) <= ds(be(jw))) exit  
        end do
        DO JJ=1,NWD  
          JWD=JWD+1  
                IF(QWD(JWD)==0.0)THEN
                    TAVGW(JWD)=-99.0
                    CAVGW(JWD,:)=-99.0
                    CDAVGW(JWD,:)=-99.0
                ENDIF
          IF (IWD(JWD) == IWDO(J)) THEN  
            JFILE=JFILE+1
            qwdo(j)           = qwdo(j)            +qwd(jwd)                 ! cb 1/16/13  
            twdo(j)           = twdo(j)            +qwd(jwd)*tavgw(jwd)                 ! cb 1/16/13  
            WRITE(WDO2(JFILE,1),'(F10.3,",",F10.4)')JDAY,QWD(JWD)  
            WRITE(WDO2(JFILE,2),'(F10.3,",",F8.2)')JDAY,TAVGW(JWD)  
            IF (CONSTITUENTS)then
              cwdo(cn(1:nac),j)= cwdo(cn(1:nac),j)+qwd(jwd)*cavgw(jwd,cn(1:nac))                 ! cb 1/16/13  
              WRITE (WDO2(JFILE,3),'(F10.3,",",*(F10.4,","))') JDAY,(CAVGW(JWD,CN(JC)),     JC=1,NAC)  
            end if
            IF (DERIVED_CALC)THEN  
              !DO JW=1,NWB  
              !  IF (IWDO(J) >= US(BS(JW)) .AND. IWDO(J) <= DS(BE(JW))) EXIT  
              !END DO
              cdwdo(cdn(1:nacd(jw),jw),j) = cdwdo(cdn(1:nacd(jw),jw),j)+qwd(jwd)*cdavgw(jwd,cdn(1:nacd(jw),jw))                 ! cb 1/16/13  
              WRITE (WDO2(JFILE,4),'(F10.3,",",*(F10.4,","))') JDAY,(CDAVGW(JWD,CDN(JD,JW)),JD=1,NACD(JW))  
            ENDIF  
          ENDIF  
        ENDDO  
        DO JS=1,NSP  ! spillways
            IF(LATERAL_SPILLWAY(JS))THEN  
              JWD=JWD+1
            ELSE
              JSSS(JBUSP(JS))=JSSS(JBUSP(JS))+1
            ENDIF
            
          IF(IWDO(J) == IUSP(JS))THEN  
            JFILE=JFILE+1  
            WRITE(WDO2(JFILE,1),'(F10.3,",",F10.4)')JDAY,QSP(JS)  
            IF(LATERAL_SPILLWAY(JS))THEN  
            !  JWD=JWD+1
              qwdo(j)           = qwdo(j)            +qsp(js)                 ! cb 1/16/13  
              twdo(j)           = twdo(j)            +qsp(js)*tavgw(jwd)                 ! cb 1/16/13
                IF(QSP(JS)==0.0)THEN
                    TAVGW(JWD)=-99.0
                    CAVGW(JWD,:)=-99.0
                    CDAVGW(JWD,:)=-99.0
                ENDIF
              WRITE(WDO2(JFILE,2),'(F10.3,",",F8.2)')JDAY,TAVGW(JWD)  
              IF (CONSTITUENTS)then
                cwdo(cn(1:nac),j)= cwdo(cn(1:nac),j)+qsp(js)*cavgw(jwd,cn(1:nac))                 ! cb 1/16/13
                WRITE (WDO2(JFILE,3),'(F10.3,",",*(F10.4,","))') JDAY,(CAVGW(JWD,CN(JC)),     JC=1,NAC)  
              end if
              IF (DERIVED_CALC)then
                cdwdo(cdn(1:nacd(jw),jw),j) = cdwdo(cdn(1:nacd(jw),jw),j)+qsp(js)*cdavgw(jwd,cdn(1:nacd(jwusp(js)),jwusp(js)))        ! cb 1/16/13
                WRITE (WDO2(JFILE,4),'(F10.3,",",*(F10.4,","))') JDAY,(CDAVGW(JWD,CDN(JD,JWUSP(JS))),JD=1,NACD(JWUSP(JS)))  
              end if
            ELSE  
            !  JSSS(JBUSP(JS))=JSSS(JBUSP(JS))+1
                qwdo(j)           = qwdo(j)            +qsp(js)                 ! cb 1/16/13  
                twdo(j)           = twdo(j)            +qsp(js)*tavg(jsss(jbusp(js)),jbusp(js))    ! cb 1/16/13
                IF(QSP(JS)==0.0)THEN
                    TAVG(JSSS(JBUSP(JS)),JBUSP(JS))=-99.0
                    CAVG(JSSS(JBUSP(JS)),JBUSP(JS),:)=-99.0
                    CDAVG(JSSS(JBUSP(JS)),JBUSP(JS),:)=-99.0
                ENDIF
              WRITE(WDO2(JFILE,2),'(F10.3,",",F8.2)')JDAY,TAVG(JSSS(JBUSP(JS)),JBUSP(JS))  
              IF (CONSTITUENTS)then
                cwdo(cn(1:nac),j)= cwdo(cn(1:nac),j)+qsp(js)*cavg(jsss(jbusp(js)),jbusp(js),cn(1:nac))            ! cb 1/16/13
                WRITE (WDO2(JFILE,3),'(F10.3,",",*(F10.4,","))') JDAY,(CAVG(JSSS(JBUSP(JS)),JBUSP(JS),CN(JC)),     JC=1,NAC)  
              end if
              IF (DERIVED_CALC)then
                cdwdo(cdn(1:nacd(jw),jw),j) = cdwdo(cdn(1:nacd(jw),jw),j)+qsp(js)*cdavg(jsss(jbusp(js)),jbusp(js),cdn(1:nacd(jwusp(js)),jwusp(js)))    ! cb 1/16/13
                WRITE (WDO2(JFILE,4),'(F10.3,",",*(F10.4,","))') JDAY,(CDAVG(JSSS(JBUSP(JS)),JBUSP(JS),CDN(JD,JWUSP(JS))),JD=1,NACD(JWUSP(JS)))  
              end if
            ENDIF  
          ENDIF  
        ENDDO  
        DO JS=1,NPU  ! PUMP
             IF(LATERAL_PUMP(JS))THEN  
              JWD=JWD+1 
             ELSE
              JSSS(JBUPU(JS))=JSSS(JBUPU(JS))+1      
            ENDIF
          IF(IWDO(J) == IUPU(JS))THEN  
            JFILE=JFILE+1  
            IF(PUMPON(JS))THEN  
              WRITE(WDO2(JFILE,1),'(F10.3,",",F10.4)')JDAY,QPU(JS)  
            ELSE  
              WRITE(WDO2(JFILE,1),'(F10.3,",",F8.3)')JDAY,0.0  
            ENDIF  
            IF(LATERAL_PUMP(JS))THEN  
            !  JWD=JWD+1              
                IF(QPU(JS)==0.0)THEN
                    TAVGW(JWD)=-99.0
                    CAVGW(JWD,:)=-99.0
                    CDAVGW(JWD,:)=-99.0
                ENDIF
                if(pumpon(js))then                 ! cb 1/16/13
                  qwdo(j)           = qwdo(j)            +qpu(js)  
                  twdo(j)           = twdo(j)            +qpu(js)*tavgw(jwd)
                end if
              WRITE(WDO2(JFILE,2),'(F10.3,",",F8.2)')JDAY,TAVGW(JWD)  
              ! Debug
              !WRITE(WDO2(JFILE,2),'(F10.3,",",F8.2,",",f8.3,",",i5,",",i5)')JDAY,TAVGW(JWD),qpu(js),js,jwd      ! Debug
              IF (CONSTITUENTS)then
                if(pumpon(js))then                 ! cb 1/16/13
                  cwdo(cn(1:nac),j)= cwdo(cn(1:nac),j)+qpu(js)*cavgw(jwd,cn(1:nac))
                end if
                WRITE (WDO2(JFILE,3),'(F10.3,",",*(F10.4,","))') JDAY,(CAVGW(JWD,CN(JC)),     JC=1,NAC)  
              end if
              IF (DERIVED_CALC)then
                if(pumpon(js))then   ! cb 1/16/13
                  cdwdo(cdn(1:nacd(jw),jw),j) = cdwdo(cdn(1:nacd(jw),jw),j)+qpu(js)*cdavgw(jwd,cdn(1:nacd(jwupu(js)),jwupu(js))) 
                end if
                WRITE (WDO2(JFILE,4),'(F10.3,",",*(F10.4,","))') JDAY,(CDAVGW(JWD,CDN(JD,JWUPU(JS))),JD=1,NACD(JWUPU(JS)))  
              end if
            ELSE  
            !  JSSS(JBUPU(JS))=JSSS(JBUPU(JS))+1  
                IF(QPU(JS)==0.0)THEN
                    TAVG(JSSS(JBUPU(JS)),JBUPU(JS))=-99.0
                    CAVG(JSSS(JBUPU(JS)),JBUPU(JS),:)=-99.0
                    CDAVG(JSSS(JBUPU(JS)),JBUPU(JS),:)=-99.0                    
                ENDIF
              if(pumpon(js))then                 ! cb 1/16/13
                  qwdo(j)           = qwdo(j)            +qpu(js)  
                  twdo(j)           = twdo(j)            +qpu(js)*tavg(jsss(jbupu(js)),jbupu(js))
              end if
              WRITE(WDO2(JFILE,2),'(F10.3,",",F8.2)')JDAY,TAVG(JSSS(JBUPU(JS)),JBUPU(JS))  
              IF (CONSTITUENTS)then
                if(pumpon(js))then                ! cb 1/16/13
                  cwdo(cn(1:nac),j)= cwdo(cn(1:nac),j)+qpu(js)*cavg(jsss(jbupu(js)),jbupu(js),cn(1:nac)) 
                end if
                WRITE (WDO2(JFILE,3),'(F10.3,",",*(F10.4,","))') JDAY,(CAVG(JSSS(JBUPU(JS)),JBUPU(JS),CN(JC)),     JC=1,NAC)  
              end if
              IF (DERIVED_CALC)then
                if(pumpon(js))then   ! cb 1/16/13
                  cdwdo(cdn(1:nacd(jw),jw),j) = cdwdo(cdn(1:nacd(jw),jw),j)+qpu(js)*cdavg(jsss(jbupu(js)),jbupu(js),cdn(1:nacd(jwupu(js)),jwupu(js))) 
                end if
                WRITE (WDO2(JFILE,4),'(F10.3,",",*(F10.4,","))') JDAY,(CDAVG(JSSS(JBUPU(JS)),JBUPU(JS),CDN(JD,JWUPU(JS))),JD=1,NACD(JWUPU(JS)))  
              end if
            ENDIF  
          ENDIF  
        ENDDO  
        DO JS=1,NPI  ! pipes
             IF(LATERAL_PIPE(JS))THEN  
             JWD=JWD+1  
             ELSE
             JSSS(JBUPI(JS))=JSSS(JBUPI(JS))+1
             ENDIF
          IF(IWDO(J) == IUPI(JS))THEN  
            JFILE=JFILE+1  
            WRITE(WDO2(JFILE,1),'(F10.3,",",F10.4)')JDAY,QPI(JS)  
            IF(LATERAL_PIPE(JS))THEN  
           !   JWD=JWD+1  
                qwdo(j)           = qwdo(j)            +qpi(js)     ! cb 1/16/13
                twdo(j)           = twdo(j)            +qpi(js)*tavgw(jwd)
                IF(QPI(JS)==0.0)THEN
                    TAVGW(JWD)=-99.0
                    CAVGW(JWD,:)=-99.0
                    CDAVGW(JWD,:)=-99.0
                ENDIF
              WRITE(WDO2(JFILE,2),'(F10.3,",",F8.2)')JDAY,TAVGW(JWD)  
              IF (CONSTITUENTS)then
                cwdo(cn(1:nac),j)= cwdo(cn(1:nac),j)+qpi(js)*cavgw(jwd,cn(1:nac))     ! cb 1/16/13
                WRITE (WDO2(JFILE,3),'(F10.3,",",*(F10.4,","))') JDAY,(CAVGW(JWD,CN(JC)),     JC=1,NAC)  
              end if
              IF (DERIVED_CALC)then
                cdwdo(cdn(1:nacd(jw),jw),j) = cdwdo(cdn(1:nacd(jw),jw),j)+qpi(js)*cdavgw(jwd,cdn(1:nacd(jwupi(js)),jwupi(js)))     ! cb 1/16/13
                WRITE (WDO2(JFILE,4),'(F10.3,*(F10.4,","))') JDAY,(CDAVGW(JWD,CDN(JD,JWUPI(JS))),JD=1,NACD(JWUPI(JS)))  
              end if
            ELSE  
            !  JSSS(JBUPI(JS))=JSSS(JBUPI(JS))+1
                qwdo(j)           = qwdo(j)            +qpi(js)     ! cb 1/16/13
                twdo(j)           = twdo(j)            +qpi(js)*tavg(jsss(jbupi(js)),jbupi(js))
                IF(QPI(JS)==0.0)THEN
                    TAVG(JSSS(JBUPI(JS)),JBUPI(JS))=-99.0
                    CAVG(JSSS(JBUPI(JS)),JBUPI(JS),:)=-99.0
                    CDAVG(JSSS(JBUPI(JS)),JBUPI(JS),:)=-99.0
                ENDIF
              WRITE(WDO2(JFILE,2),'(F10.3,",",F8.2)')JDAY,TAVG(JSSS(JBUPI(JS)),JBUPI(JS))  
              IF (CONSTITUENTS)then
                cwdo(cn(1:nac),j)= cwdo(cn(1:nac),j)+qpi(js)*cavg(jsss(jbupi(js)),jbupi(js),cn(1:nac))     ! cb 1/16/13
                WRITE (WDO2(JFILE,3),'(F10.3,",",*(F10.4,","))') JDAY,(CAVG(JSSS(JBUPI(JS)),JBUPI(JS),CN(JC)),     JC=1,NAC)  
              end if
              IF (DERIVED_CALC)then
                cdwdo(cdn(1:nacd(jw),jw),j) = cdwdo(cdn(1:nacd(jw),jw),j)+qpi(js)*cdavg(jsss(jbupi(js)),jbupi(js),cdn(1:nacd(jwupi(js)),jwupi(js)))     ! cb 1/16/13
                WRITE (WDO2(JFILE,4),'(F10.3,",",*(F10.4,","))') JDAY,(CDAVG(JSSS(JBUPI(JS)),JBUPI(JS),CDN(JD,JWUPI(JS))),JD=1,NACD(JWUPI(JS)))  
              end if
            ENDIF  
          ENDIF  
        ENDDO  
        DO JS=1,NGT  ! gates
             IF(LATERAL_GATE(JS))THEN  
              JWD=JWD+1 
             ELSE
              JSSS(JBUGT(JS))=JSSS(JBUGT(JS))+1 
             ENDIF
             
          IF(IWDO(J) == IUGT(JS))THEN  
            JFILE=JFILE+1  
            WRITE(WDO2(JFILE,1),'(F10.3,",",F10.4)')JDAY,QGT(JS)  
            IF(LATERAL_GATE(JS))THEN  
          !    JWD=JWD+1  
              qwdo(j)           = qwdo(j)            +qgt(js)     ! cb 1/16/13
              twdo(j)           = twdo(j)            +qgt(js)*tavgw(jwd)
                IF(QGT(JS)==0.0)THEN
                    TAVGW(JWD)=-99.0
                    CAVGW(JWD,:)=-99.0
                    CDAVGW(JWD,:)=-99.0
                ENDIF
              WRITE(WDO2(JFILE,2),'(F10.3,",",F8.2)')JDAY,TAVGW(JWD)  
              IF (CONSTITUENTS)then
                cwdo(cn(1:nac),j)= cwdo(cn(1:nac),j)+qgt(js)*cavgw(jwd,cn(1:nac))     ! cb 1/16/13
                WRITE (WDO2(JFILE,3),'(F10.3,",",*(F10.4,","))') JDAY,(CAVGW(JWD,CN(JC)),     JC=1,NAC)  
              end if
              IF (DERIVED_CALC)then
                cdwdo(cdn(1:nacd(jw),jw),j) = cdwdo(cdn(1:nacd(jw),jw),j)+qgt(js)*cdavgw(jwd,cdn(1:nacd(jwugt(js)),jwugt(js)))     ! cb 1/16/13
                WRITE (WDO2(JFILE,4),'(F10.3,",",*(F10.4,","))') JDAY,(CDAVGW(JWD,CDN(JD,JWUGT(JS))),JD=1,NACD(JWUGT(JS)))  
              end if
            ELSE  
           !   JSSS(JBUGT(JS))=JSSS(JBUGT(JS))+1
               qwdo(j)           = qwdo(j)            +qgt(js)     ! cb 1/16/13
               twdo(j)           = twdo(j)            +qgt(js)*tavg(jsss(jbugt(js)),jbugt(js))
                IF(QGT(JS)==0.0)THEN
                    TAVG(JSSS(JBUGT(JS)),JBUGT(JS))=-99.0
                    CAVG(JSSS(JBUGT(JS)),JBUGT(JS),:)=-99.0
                    CDAVG(JSSS(JBUGT(JS)),JBUGT(JS),:)=-99.0
                ENDIF
              WRITE(WDO2(JFILE,2),'(F10.3,",",F8.2)')JDAY,TAVG(JSSS(JBUGT(JS)),JBUGT(JS))  
              IF (CONSTITUENTS)then
                cwdo(cn(1:nac),j)= cwdo(cn(1:nac),j)+qgt(js)*cavg(jsss(jbugt(js)),jbugt(js),cn(1:nac))     ! cb 1/16/13
                WRITE (WDO2(JFILE,3),'(F10.3,",",*(F10.4,","))') JDAY,(CAVG(JSSS(JBUGT(JS)),JBUGT(JS),CN(JC)),     JC=1,NAC)  
              end if
              IF (DERIVED_CALC)then
                cdwdo(cdn(1:nacd(jw),jw),j) = cdwdo(cdn(1:nacd(jw),jw),j)+qgt(js)*cdavg(jsss(jbugt(js)),jbugt(js),cdn(1:nacd(jwugt(js)),jwugt(js)))     ! cb 1/16/13
                WRITE (WDO2(JFILE,4),'(F10.3,",",*(F10.4,","))') JDAY,(CDAVG(JSSS(JBUGT(JS)),JBUGT(JS),CDN(JD,JWUGT(JS))),JD=1,NACD(JWUGT(JS)))  
              end if
            ENDIF  
          ENDIF  
        ENDDO  
        
        ! cb 1/16/13 deleted old withdrawal output code
          
        IF (QWDO(J) /= 0.0) TWDO(J) = TWDO(J)/QWDO(J)  
        DO JC=1,NAC  
          IF (QWDO(J) /= 0.0) CWDO(CN(JC),J) = CWDO(CN(JC),J)/QWDO(J)  
          WRITE (CWDOC(CN(JC)),'(F10.4)') CWDO(CN(JC),J)                        ! SW 9/23/13 Changed format from G8.3 to F8.3 to avoid format overflow
          CWDOC(CN(JC)) = ADJUSTR(CWDOC(CN(JC)))  
        END DO  
        DO JW=1,NWB  
          IF (IWDO(J) >= US(BS(JW)) .AND. IWDO(J) <= DS(BE(JW))) EXIT  
        END DO  
        DO JD=1,NACD(JW)  
          IF (QWDO(J) /= 0.0) CDWDO(CDN(JD,JW),J) = CDWDO(CDN(JD,JW),J)/QWDO(J)  
          WRITE (CDWDOC(CDN(JD,JW)),'(F10.4)') CDWDO(CDN(JD,JW),J)              ! SW 9/23/13 Changed format from G8.3 to F8.3 to avoid format overflow
          CDWDOC(CDN(JD,JW)) = ADJUSTR(CDWDOC(CDN(JD,JW)))  
        END DO  
2499      WRITE (WDO(J,1),'(F10.3,",",F9.3,",",8X,*(F9.3,","))',ERR=2500) JDAY, QWDO(J), (QOUTLET(I),I=1,NUMOUTLETS)     ! sw 3/2019 This code was necessary during multiple WB read/write possible error with one file copying and another writing at the same time
          GO TO 2501
2500      CALL SLEEPQQ(100)          
          WRITE(9911,'(A,f12.3,A,f12.3,A,f12.3)')'ERROR writing qwo file output on JDAY',JDAY,' retrying read after pausing 0.1 s'
          go to 2499
2501      CONTINUE
          WRITE (WDO(J,2),'(F10.3,",",F8.2,",",8X,*(F8.2,","))',ERR=2502) JDAY, TWDO(J), (TOUTLET(I),I=1,NUMOUTLETS)  
          GO TO 2503
2502      CALL SLEEPQQ(100)          
          WRITE(9911,'(A,f12.3,A,f12.3,A,f12.3)')'ERROR writing two file output on JDAY',JDAY,' retrying read after pausing 0.1 s'
          go to 2501
2503      CONTINUE
          IF (CONSTITUENTS) WRITE (WDO(J,3),'(F10.3,",",*(A10,","))',ERR=2504) JDAY,(CWDOC(CN(JC)),     JC=1,NAC)  
          GO TO 2505
2504      CALL SLEEPQQ(100)          
          WRITE(9911,'(A,f12.3,A,f12.3,A,f12.3)')'ERROR writing cwo file output on JDAY',JDAY,' retrying read after pausing 0.1 s'
          go to 2503
2505      CONTINUE
          IF (DERIVED_CALC) WRITE (WDO(J,4),'(F10.3,",",*(A10,","))') JDAY,(CDWDOC(CDN(JD,JW)),JD=1,NACD(JW))  
      END DO  
    END IF  
  END IF  
    
  !**** DSI W2 Linkage File (W2L) (Supercedes Old Velocity vectors)  
  IF (VECTOR(1)) THEN  
    ! *** Apply the same linkage settings for all waterbodies  
    IF (JDAY >= NXTMVP(1) .OR. JDAY >= VPLD(VPLDP(1)+1,1)) THEN  
      IF (JDAY >= VPLD(VPLDP(1)+1,1)) THEN  
        VPLDP(1)  = VPLDP(1)+1  
        NXTMVP(1) = VPLD(VPLDP(1),1)  
      END IF  
      NXTMVP(1) = NXTMVP(1)+VPLF(VPLDP(1),1)  
        
      ! *** Write the W2L Snapshot (DSI)
      WRITE(VPL(1)) REAL(JDAY,4)  
        
      ! *** Compute the elevation, with ELWS zeroed for segments < CUS  
      DO JW=1,NWB  
        DO JB=BS(JW),BE(JW)  
          DO I=US(JB)-1,DS(JB)+1  
            IF (I.LT.CUS(JB)) THEN  
              WSEL(I) = -9999  
            ELSE  
              WSEL(I) = ELWS(I)  
            END IF  
          END DO  
        END DO  
      END DO  
        
      WRITE(VPL(1)) ((WSEL(I)),I=1,IMX)  
      WRITE(VPL(1)) ((REAL(U(K,I),4),K=1,KMX),I=1,IMX)  
      WRITE(VPL(1)) ((REAL(W(K,I),4),K=1,KMX),I=1,IMX)  
        
      DO I = 1,IMX
        ! *** ABOVE ACTIVE LAYER 
        DO K = 1, KTI(I)-1  
          WDSI(K,I) = -9999.  
            END DO  
        ! *** ACTIVE LAYERS
        DO K = KTI(I),KB(I)
          WDSI(K,I) = T2(K,I)
          END DO  
        ! *** BELOW ACTIVE LAYER 
        DO K = KB(I)+1,KMX  
          WDSI(K,I) = -9999.  
        END DO  
      END DO  
      WRITE(VPL(1)) ((WDSI(K,I),K=1,KMX),I=1,IMX)  

      ! *** Constituent data  
      DO JC=1,NAC  
        IC=CN(JC)  
        DO I = 1,IMX
          ! *** ABOVE ACTIVE LAYER 
          DO K = 1, KTI(I)-1  
            WDSI(K,I) = -9999.  
          END DO  
          ! *** ACTIVE LAYERS
          DO K = KTI(I),KB(I)
            WDSI(K,I) = C2(K,I,IC)
          END DO  
          ! *** BELOW ACTIVE LAYER 
          DO K = KB(I)+1,KMX  
            WDSI(K,I) = -9999.  
          END DO  
        END DO  
        WRITE(VPL(1)) ((WDSI(K,I),K=1,KMX),I=1,IMX)  
      END DO  
        
    END IF  
  END IF  
    
  !** Restart  
    
  IF (RESTART_OUT) THEN  
    IF (JDAY >= NXTMRS .OR. JDAY >= RSOD(RSODP+1)) THEN  
      IF (JDAY >= RSOD(RSODP+1)) THEN  
        RSODP  = RSODP+1  
        NXTMRS = RSOD(RSODP)  
      END IF  
      NXTMRS = NXTMRS+RSOF(RSODP)  
      IF(RSOF(RSODP) >= 1.0)THEN  
        WRITE (EXT,'(I0)') INT(JDAY)  
      ELSE  
        WRITE(EXT,'(I0,"_",I2)')INT(JDAY),INT(100.*(JDAY-INT(JDAY)))   ! Allows for writing out file names with JDAY FREQ less than 1 day  
      ENDIF  
      EXT   = ADJUSTL(EXT)  
      L     = LEN_TRIM(EXT)  
      RSOFN = 'rso'//EXT(1:L)//'.opt'  
      CALL RESTART_OUTPUT (RSOFN)  
    END IF  
  END IF  
    
  ! Snapshot formats  
    
10490 FORMAT ('CE-QUAL-W2 VERSION',F4.2/                                                                                          &  
              (1X,A72))  
10500 FORMAT (/1X,A/                                                                                                               &  
              3X,'Gregorian date      [GDAY] =',A19,1X,I0,', ',I0/                                                                 &  
              3X,'Julian date         [JDAY] =',I10,' days',F6.2,' hours'/                                                         &  
              3X,'Elapsed time      [ELTMJD] =',I10,' days',F6.2,' hours'/                                                         &  
              3X,'Timestep             [DLT] =',I10,' sec'/                                                                        &  
              3X,'  at location  [KLOC,ILOC] = (',I0,',',I0,')'/                                                                   &  
              3X,'Minimum timestep  [MINDLT] =',I10,' sec '/                                                                       &  
              3X,'  at Julian day    [JDMIN] =',I10,' days',F6.2,' hours'/                                                         &  
              3X,'  at location  [KMIN,IMIN] = (',I0,',',I0,')')  
10510 FORMAT (3X,'Limiting timestep'/                                                                                              &  
              3X,'  at location  [KMIN,IMIN] = (',I0,',',I0,')')  
10520 FORMAT (3X,'Average timestep   [DLTAV] =',I10,' sec'/                                                                        &  
              3X,'Number of iterations [NIT] =',I10/                                                                               &  
              3X,'Number of violations  [NV] =',I10/)  
10530 FORMAT (1X,A)  
10540 FORMAT (3X,'Input'/                                                                                                          &  
              3X,'  Air temperature          [TAIR] =',F9.2,1X,A/                                                                  &  
              3X,'  Dewpoint temperature     [TDEW] =',F9.2,1X,A/                                                                  &  
              3X,'  Wind direction            [PHI] =',F9.2,' rad'/                                                                &  
              3X,'  Cloud cover             [CLOUD] =',F9.2/                                                                       &  
              3X,'  Calculated'/                                                                                                   &  
              5X,'  Equilibrium temperature    [ET] =',F9.2,1X,A/                                                                  &  
              5X,'  Surface heat exchange    [CSHE] =',E9.2,' m/sec'/                                                              &  
              5X,'  Net short wave radiation [SRON] =',E9.2,1X,A,' W/m^2'/)  
10550 FORMAT (1X,A/                                                                                                                &  
              3X,A)  
10560 FORMAT (5X,'Branch ',I0/                                                                                                     &  
              5X,'  Layer       [KQIN] = ',I0,'-',I0/                                                                              &  
              5X,'  Inflow       [QIN] =',F8.2,' m^3/sec'/                                                                         &  
              5X,'  Temperature  [TIN] =',F8.2,1X,A)  
10570 FORMAT (/3X,'Distributed Tributaries')  
10580 FORMAT (5X,'Branch ',I0/                                                                                                     &  
              5X,'  Inflow      [QDTR] =',F8.2,' m^3/sec'/                                                                         &  
              5X,'  Temperature [TDTR] =',F8.2,1X,A)  
10590 FORMAT (:/3X,'Tributaries'/                                                                                                  &  
              5X,'Segment     [ITR] =',11I8:/                                                                                      &  
              (T25,11I8))  
10600 FORMAT (:5X,'Layer      [KTWB] = ',11(I0,'-',I0,2X):/                                                                        &  
              (T25,11(I0,'-',I0)))  
10610 FORMAT (:5X,'Inflow      [QTR] =',11F8.2:/                                                                                   &  
              (T25,11F8.1))  
10620 FORMAT (:5X,'Temperature [TTR] =',11F8.2:/                                                                                   &  
              (T25,11F8.1))  
10630 FORMAT (/1X,'Outflows')  
10640 FORMAT (3X,'Structure outflows [QSTR]'/                                                                                      &  
              3X,'  Branch ',I0,' = ',11F8.2:/                                                                                     &  
              (T16,11F8.2))  
10650 FORMAT (:/3X,'Total outflow [QOUT] =',F8.2,' m^3/s'/                                                                         &  
              5X,'Outlets'/                                                                                                        &  
              5X,'  Layer             [KOUT] =',12I7:/                                                                             &  
              (33X,12I7))  
10660 FORMAT (:7X,'Outflow (m^3/sec) [QOUT] =',12F7.2:/                                                                            &  
              (33X,12F7.2))  
10670 FORMAT (:5X,'Withdrawals'/                                                                                                   &  
              5X,'  Segment            [IWD] =',I7/                                                                                &  
              5X,'  Outflow (m^3/sec)  [QWD] =',F7.2)  
10680 FORMAT (5X,'  Layer              [KWD] =',12I7/                                                                              &  
              (33X,12I7))  
10690 FORMAT (:5X,'  Outflow (m^3/sec)  [QSW] =',12F7.2/                                                                           &  
              (33X,12F7.2))  
10700 FORMAT (/'1',A)  
10710 FORMAT (3X,'Branch ',I0,' [CIN]'/                                                                                            &  
              (5X,A,T25,'=',F9.3,1X,A))  
10720 FORMAT (3X,'Tributary ',I0,' [CTR]'/                                                                                         &  
              (5X,A,T25,'=',F9.3,1X,A))  
10730 FORMAT (3X,'Distributed tributary ',I0,' [CDT]'/                                                                             &  
              (5X,A,T25,'=',F9.3,1X,A))  
10740 FORMAT (/'Surface calculations')  
10750 FORMAT (3X,'Evaporation rate [EV]'/                                                                                          &  
              (:3X,'  Branch ',I0,' = ',E10.3,' m^3/s'))  ! SW 9/15/05 4/21/10
10755 format (3x,'Cumulative evaporation [VOLEV]'/                                                                                 &  
              (:3X,'  Branch ',I0,' = ',F0.1,' m^3'))  
10760 FORMAT (3X,'Precipitation [PR]'/                                                                                             &  
              (3X,'  Branch ',I0,' = ',F8.6),' m/s')  
10770 FORMAT (/1X,'External head boundary elevations'/)  
10780 FORMAT (3X,'Branch ',I0/5X,'Upstream elevation   [ELUH] =',F8.3,' m')  
10790 FORMAT (3X,'Branch ',I0/5X,'Downstream elevation [ELDH] =',F8.3,' m')  
10800 FORMAT (/'Water Balance')  
10810 FORMAT (3X,'Waterbody ',I0/                                                                                                  &  
              3X,'  Spatial change  [VOLSR]  = ',E15.8,' m^3'/                                                                     &  
              3X,'  Temporal change [VOLTR]  = ',E15.8,' m^3'/                                                                     &  
              3X,'  Volume error             = ',E15.8,' m^3'/                                                                     &  
              3X,'  Percent error            = ',E15.8,' %')  
10820 FORMAT (3X,'Branch ',I0/                                                                                                     &  
              3X,'  Spatial change  [VOLSBR] = ',E15.8,' m^3'/                                                                     &  
              3X,'  Temporal change [VOLTBR] = ',E15.8,' m^3'/                                                                     &  
              3X,'  Volume error             = ',E15.8,' m^3'/                                                                     &  
              3X,'  Percent error            = ',E15.8,' %')  
10830 FORMAT (/1X,'Energy Balance')  
10840 FORMAT (3X,'Waterbody ',I0/                                                                                                  &  
              3X,'  Spatially integrated energy   [ESR] = ',E15.8,' kJ'/                                                           &  
              3X,'  Temporally integrated energy  [ETR] = ',E15.8,' kJ'/                                                           &  
              3X,'  Energy error                        = ',E15.8,' kJ'/                                                           &  
              3X,'  Percent error                       = ',E15.8,' %')  
10850 FORMAT (3X,'  Spatially integrated energy  [ESBR] = ',E15.8,' kJ'/                                                           &  
              3X,'  Temporally integrated energy [ETBR] = ',E15.8,' kJ'/                                                           &  
              3X,'  Energy error                        = ',E15.8,' kJ'/                                                           &  
              3X,'  Percent error                       = ',E15.8,' %')  
10860 FORMAT (/1X,'Mass Balance')  
10870 FORMAT (3X,'Branch ',I0)  
10880 FORMAT (5X,A/                                                                                                                &  
              5X,'  Spatially integrated mass  [CMBRS] = ',E15.8,1X,A/                                                             &  
              5X,'  Temporally integrated mass [CMBRT] = ',E15.8,1X,A/                                                             &  
              5X,'  Mass error                         = ',E15.8,1X,A/                                                             &  
              5X,'  Percent error                      = ',E15.8,' %')  
10890 FORMAT (/1X,A/                                                                                                               &  
              3X,'Surface layer [KT] = ',I0/                                                                                       &  
              3X,'Elevation   [ELKT] =',F10.3,' m')  
10900 FORMAT (/3X,'Current upstream segment [CUS]'/                                                                                &  
             (3X,'  Branch ',I0,' = ',I0))  
    
  RETURN  
  
END SUBROUTINE OUTPUTA  
