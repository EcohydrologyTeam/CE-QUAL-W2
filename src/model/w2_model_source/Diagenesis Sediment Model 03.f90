MODULE SEDMODEL_LOCAL
    Use MAIN 
    Use GLOBAL
    Use GEOMC
    Use SCREENC
    Use RSTART
    Use PREC
    Use EDDY
    Use LOGICC
    Use TVDC
    Use KINETIC
    Use TRANS
    Use CEMAVars
    Use CEMASedimentDiagenesis, only: SD_T1 
    ! Type declarations
    IMPLICIT NONE
    Integer(4) :: RegnNum
    Real(8)    :: VolumeofSedimentBed1, VolumeofSedimentBed2
    Real(8)    :: VolumeofPorewater1, VolumeofPorewater2
    Real(8)    :: VolumePresent, VolumeNew, VolumeUpdated
    Real(8)    :: SEChange, SurfaceAreaWB, VolumeIncreased
    Real(8)    :: CritShldPar, TCrit_E, TCrit_S, VelTemp, SScour
    Real(8)    :: TempVariable, TFlow_B, CEMAShieldsNumber
    Real(8)    :: Sgr, Dks, Uks, NuTemp, Dsks, Tsks, Usr, SedECoeff
    Real(8)    :: d1, d2, Depth, UctAvg, Wf, BetaD, Sum1, Sdfz
    Real(8)    :: rgha, rghapda, Caks, Srho, Capda, Yalinp, Value1
    Real(8)    :: SNetFlx, SSettle
    Real       :: Dummy
    Logical    :: BottomUpdated
    END MODULE
    
    Subroutine CEMASedimentModelW2
    Use SEDMODEL_LOCAL

    NMFT=0
    DO JG=1,NSS  
      IF(SSCS(JG)==-1.0)THEN
        NMFT=NSSS+JG-1        
        EXIT
      ENDIF
    ENDDO

    Entry SetupCEMASedimentModel
        
        !Call to initialize region consolidation data files
    Call InitializeBedConsolidationFiles(1210, ConsolidRateRegnFil)
        
        Do RegnNum = 1, NumConsolidRegns
            Do SegNumI = ConsRegSegSt(RegnNum), ConsRegSegEn(RegnNum)
                ConsolidRegnNum(SegNumI) = RegnNum
                bedconsolidationseg(segnumI)=.true.  ! cb 6/28/18
            End Do !SegNumI
        End Do !RegnNum 
                
        !Compute total porewater content
        TotalPoreWatVolume = 0.d00
        TotalSedimentsInBed = 0.d00
        Do JW=1, NWB
            KT = KTWB(JW)
            Do JB=BS(JW),BE(JW)
                IU = CUS(JB)
                ID = DS(JB)
                Do SegNumI = IU, ID
                    TotalPoreWatVolume = TotalPoreWatVolume + BedPorosity(SegNumI)*BedElevation(SegNumI)*B(KB(SegNumI),SegNumI)*DLX(SegNumI)
                    TotalSedimentsInBed = TotalSedimentsInBed + (1-BedPorosity(SegNumI))*BedElevation(SegNumI)*B(KB(SegNumI),SegNumI)*DLX(SegNumI)   
                   ! TotalPoreWatVolume = TotalPoreWatVolume + BedPorosity(SegNumI)*BedElevation(SegNumI)* sedcellwidth(SegNumI)*DLX(SegNumI)
                   ! TotalSedimentsInBed = TotalSedimentsInBed + (1-BedPorosity(SegNumI))*BedElevation(SegNumI)* sedcellwidth(SegNumI)*DLX(SegNumI)   
                End Do !SegNumI
            End Do !JB
        End Do !JW
              
    ! moved file opens and writes to CEMASedimentDiagenesis          
    Return


    Entry CEMASedimentModel
    
        PorewaterRelRate = 0.d00
        Do JW=1, NWB
            KT = KTWB(JW)
            Do JB=BS(JW),BE(JW)
                IU = CUS(JB)
                ID = DS(JB)
                Do SegNumI = IU, ID
                    If(CEMASSApplied(SegNumI))Then
                        CEMALayerAdded(SegNumI) = .FALSE. 
                        CEMACumPWRelease(SegNumI) = 0.d00
                        CEMASSApplied(SegNumI) = .FALSE.
                    End If
                End Do !SegNumO
            End Do !JB
        End Do !JW
        
        !Get Suspended sediment concentration
        Do JW=1, NWB
            KT = KTWB(JW)
            Do JB=BS(JW),BE(JW)
                IU = CUS(JB)
                ID = DS(JB)
                Do SegNumI = IU, ID
                    Do k = kt, KB(SegNumI)
                    
                        !CEMASedConc(SegNumI, K) = C1(K,SegNumI,NSSS)
                    CEMASedConc(SegNumI, K) = C1(K,SegNumI,NMFT)      ! cb 2/18/13
                    
                    End Do !k
                End Do !SegNumI
            End Do !JB
        End Do !JW
        
        !If(EndBedConsolidation)Return
        
        !Get bed elevation change rate
        Do JW=1, NWB
            KT = KTWB(JW)
            Do JB=BS(JW),BE(JW)
                
                IU = CUS(JB)
                ID = DS(JB)
                
                Do SegNumI = IU, ID
                    
                    If(EndBedConsolidation(SegNumI))Cycle
                    If(.not. BedConsolidationSeg(SegNumI))Cycle  ! cb 6/28/18
                    
                    RegnNum = ConsolidRegnNum(SegNumI)  
                    If(ConsolidationType(RegnNum) == 0)Then    !Constant value
                        BedConsolidRate(SegNumI) = ConstConsolidRate(RegnNum)/86400.d00 !m/d to m/s
                    End If
                    If(ConsolidationType(RegnNum) == 1)Then    !Time varying value
                    Call ReadBedConsolidationFiles(1210)
                    BedConsolidRate(SegNumI) = ConsolidRateTemp(RegnNum)/86400.d00 !m/d to m/s
                    End If
                
                End Do !SegNumI
                
            End Do !JB
        End Do !JW 
        
        !Compute net sediment bed elevation
        TotalPoreWatRemoved = 0.d00
        Do JW=1, NWB
            KT = KTWB(JW)
            Do JB=BS(JW),BE(JW)
                IU = CUS(JB)
                ID = DS(JB)
                Do SegNumI = IU, ID
                
                    If(EndBedConsolidation(SegNumI))Cycle
                    If(.not. BedConsolidationSeg(SegNumI))Cycle  ! cb 6/28/18
                    
                    VolumeofSedimentBed1 = BedElevation(SegNumI)*B(KB(SegNumI),SegNumI)*DLX(SegNumI)
                    !VolumeofSedimentBed1 = BedElevation(SegNumI)*B(KT,SegNumI)*DLX(SegNumI)
                    VolumeofPorewater1 = VolumeofSedimentBed1*BedPorosity(SegNumI)
                    BedElevation(SegNumI) = BedElevation(SegNumI) - BedConsolidRate(SegNumI)*dlt
                    BedElevationLayer(SegNumI) = BedElevationLayer(SegNumI) - BedConsolidRate(SegNumI)*dlt                    
                    VolumeofSedimentBed2 = BedElevation(SegNumI)*B(KB(SegNumI),SegNumI)*DLX(SegNumI)
                    !VolumeofSedimentBed2 = BedElevation(SegNumI)*B(KT,SegNumI)*DLX(SegNumI)
                    BedPorosity(SegNumI) = 1.d0 - VolumeofSedimentBed1/VolumeofSedimentBed2*(1-BedPorosity(SegNumI))
                    VolumeofPorewater2 = VolumeofSedimentBed2*BedPorosity(SegNumI)
                    TotalPoreWatRemoved = VolumeofPorewater1 - VolumeofPorewater2
                    PorewaterRelRate(SegNumI) = TotalPoreWatRemoved/dlt !m�/s
                    CEMACumPWRelease(SegNumI) = CEMACumPWRelease(SegNumI) + TotalPoreWatRemoved     !m�
                    CEMACumPWReleaseRate(SegNumI) = CEMACumPWRelease(SegNumI)/dlt !m�/s
                    
                    If(BedPorosity(SegNumI) < 0.d00)Then
                        TotalPoreWatRemoved = 0.d00
                        BedPorosity(SegNumI) = 0.d00
                        EndBedConsolidation(SegNumI) = .TRUE. 
                        
                        Write(W2ERR, '(a)')"Bed reached maximum compaction level. Porosity = 0"
                        Write(W2ERR, '(a,i4)')"No more bed compaction will occur at segment = ",SegNumI
                        ERROR_OPEN=.TRUE.    ! SR 7/27/2017
                    End If
                    
                End Do !SegNumI
            End Do !JB
        End Do !JW
        
        !Compute net sediment and porewater volume
        TotalPoreWatVolume = 0.d00
        TotalSedimentsInBed = 0.d00
        Do JW=1, NWB
            KT = KTWB(JW)
            Do JB=BS(JW),BE(JW)
                IU = CUS(JB)
                ID = DS(JB)
                Do SegNumI = IU, ID
                    
                    If(EndBedConsolidation(SegNumI))Cycle                                    
                    If(.not. BedConsolidationSeg(SegNumI))Cycle  ! cb 6/28/18
                    TotalPoreWatVolume = TotalPoreWatVolume + BedPorosity(SegNumI)*BedElevation(SegNumI)*B(KB(SegNumI),SegNumI)*DLX(SegNumI)   
                    TotalSedimentsInBed = TotalSedimentsInBed + (1-BedPorosity(SegNumI))*BedElevation(SegNumI)*B(KB(SegNumI),SegNumI)*DLX(SegNumI)   
                    !TotalPoreWatVolume = TotalPoreWatVolume + BedPorosity(SegNumI)*BedElevation(SegNumI)*B(KT,SegNumI)*DLX(SegNumI)   
                    !TotalSedimentsInBed = TotalSedimentsInBed + (1-BedPorosity(SegNumI))*BedElevation(SegNumI)*B(KT,SegNumI)*DLX(SegNumI)   
                                    
                End Do !SegNumI
            End Do !JB
        End Do !JW
        
        !Get Erosion/Deposition
        If(CEMASedimentProcessesInc)Call CEMASedimentProcesses
        
        !Update Suspended sediment concentration
        Do JW=1, NWB
            KT = KTWB(JW)
            Do JB=BS(JW),BE(JW)
                IU = CUS(JB)
                ID = DS(JB)
                Do SegNumI = IU, ID
                    Do k = kt, KB(SegNumI)
                    
                        !C1(K,SegNumI,NSSS) = CEMASedConc(SegNumI, K)
                    C1(K,SegNumI,NMFT) = CEMASedConc(SegNumI, K)     ! cb 2/18/13
                    
                    End Do !k
                End Do !SegNumI
            End Do !JB
        End Do !JW
        
    Return
       
    Entry CEMAUpdateVerticalLayering
    
        !Compute current volume
        VolumePresent = 0.d00
        Do JW=1, NWB
            KT = KTWB(JW)
            Do JB=BS(JW),BE(JW)
                IU = CUS(JB)
                ID = DS(JB)
                Do SegNumI = IU, ID
                    If(.not. BedConsolidationSeg(SegNumI))Cycle  ! cb 7/5/18
                    Do K = KT, KB(SegNumI)
                        If(K == KT)VolumePresent = VolumePresent + B(K,SegNumI)*DLX(SegNumI)*(H(K,JW)-Z(SegNumI))
                        If(K /= KT)VolumePresent = VolumePresent + B(K,SegNumI)*DLX(SegNumI)*H(K,JW)
                    End Do !K
                End Do !SegNumI
            End Do !JB
        End Do !JW   
        
        BottomUpdated = .FALSE.
        
        Do JW=1, NWB
            KT = KTWB(JW)
            Do JB=BS(JW),BE(JW)
                IU = CUS(JB)
                ID = DS(JB)
                Do SegNumI = IU, ID
                    
                    If(EndBedConsolidation(SegNumI))Cycle
                    
                    If(dabs(BedElevationLayer(SegNumI)) > LayerAddThkFrac*H(KB(SegNumI),JW))Then  !Add water layer
                        
                        BedElevationLayer(SegNumI) = BedElevationLayer(SegNumI) + H(KB(SegNumI),JW) !BedElevationLayer is negative due to bed consolidation
                        KB(SegNumI) = KB(SegNumI) + 1    !Add layer
                        
                        If(KB(SegNumI) > KMX-1)Then
                            KB(SegNumI) = KB(SegNumI) - 1
                            Write(W2ERR, '(a,i4)')"Maximum K layer reached for Segment = ",SegNumI
                            ERROR_OPEN=.TRUE.    ! SR 7/27/2017
                        Else
                            K = KB(SegNumI)
                            !Copy properties from old bottom layer
                            B(K,SegNumI)                =   B(K-1,SegNumI)
                            RHO(K,SegNumI)              =   RHO(K-1,SegNumI)
                            H1(K,SegNumI)               =   H(K,JW)
                            H2(K,SegNumI)               =   H(K,JW)
                            AVH1(K,SegNumI)             =   (H1(K,SegNumI)  + H1(K-1,SegNumI))*0.5
                            IF (.NOT. TRAPEZOIDAL(JW)) THEN
                                BH1(K,SegNumI)          =   BH1(K-1,SegNumI)-Bnew(K-1,SegNumI)*H1(K-1,SegNumI)                              ! SW 1/23/06
                            ELSE
                                CALL GRID_AREA1(EL(K,SegNumI)-Z(SegNumI),EL(K-1,SegNumI),BH1(K,SegNumI),DUMMY)                                                          !SW 08/03/04
                                BH1(K,SegNumI)          =   0.25*H1(K-1,JW)*(BB(K,SegNumI)+2.*B(K-1,SegNumI)+BB(K-1,SegNumI))
                            ENDIF
                            !VOL(K,SegNumI)              =   BH1(KT,SegNumI)  *DLX(SegNumI)
                            VOL(K,SegNumI)              =   BH1(K,SegNumI)  *DLX(SegNumI)
                            BKT(SegNumI)                =   BH1(KT,SegNumI)/H1(KT,SegNumI)
                            BI(KT:KB(SegNumI),SegNumI)  =   B(KT:KB(SegNumI),SegNumI)   ! SW 8/26/05
                            !T1(K,SegNumI)               =   T1(K-1,SegNumI)
                            !T2(K,SegNumI)               =   T2(K-1,SegNumI)
                            !C1(K,SegNumI,CN(1:NAC))     =   C1(K-1,SegNumI,CN(1:NAC))
                            !C2(K,SegNumI,CN(1:NAC))     =   C2(K-1,SegNumI,CN(1:NAC))
                            AZ(K-1,SegNumI)             =   AZ(K-2,SegNumI)
                            
                            !Copy Constituent properties from old bottom layer
                            !NH4(K,SegNumI)                =   NH4(K-1,SegNumI)
                            !NO3(K,SegNumI)                =   NO3(K-1,SegNumI)
                            !O2(K,SegNumI)                =   O2(K-1,SegNumI)
                            !C2(K,SegNumI,NDO)     =   C2(K-1,SegNumI,NDO)
                            
                            !SP added to test 04/01/2010
                            !DX(K,SegNumI)                =   DX(K-1,SegNumI)
                            !End SP added to test 04/01/2010
                            
                            BottomUpdated = .TRUE.
                            CEMALayerAdded(SegNumI) = .TRUE.
                            
                            Write(CEMALogFilN,*)"Layer added at Segment Number ",SegNumI, "on ", JDAY
                            Write(wrn,*)"Sed Diag:Layer added at Segment Number ",SegNumI, "on ", JDAY
                            
                        End If
                        
                    End If
                    
                End Do !SegNumI
            End Do !JB
        End Do !JW 
        
        !Compute new volume
        VolumeNew = 0.d00
        Do JW=1, NWB
            KT = KTWB(JW)
            Do JB=BS(JW),BE(JW)
                IU = CUS(JB)
                ID = DS(JB)
                Do SegNumI = IU, ID
                    If(.not. BedConsolidationSeg(SegNumI))Cycle  ! cb 7/5/18
                    Do K = KT, KB(SegNumI)
                        If(K == KT)VolumeNew = VolumeNew + B(K,SegNumI)*DLX(SegNumI)*(H(K,JW)-Z(SegNumI))
                        If(K /= KT)VolumeNew = VolumeNew + B(K,SegNumI)*DLX(SegNumI)*H(K,JW)
                    End Do !K
                End Do !SegNumI
            End Do !JB
        End Do !JW
        
        If(BottomUpdated)Then   !Update surface elevation to accomodate for layer addition
            
            VolumeIncreased = VolumeNew - VolumePresent
            
            !Calculate surface area
            SurfaceAreaWB = 0.d00
            Do JW=1, NWB
                KT = KTWB(JW)
                Do JB=BS(JW),BE(JW)
                    IU = CUS(JB)
                    ID = DS(JB)
                    Do SegNumI = IU, ID
                        If(.not. BedConsolidationSeg(SegNumI))Cycle  ! cb 7/5/18
                        SurfaceAreaWB = SurfaceAreaWB + B(KT,SegNumI)*DLX(SegNumI)
                    End Do !SegNumI
                End Do !JB
            End Do !JW
            
            SEChange = VolumeIncreased/SurfaceAreaWB
            Do JW=1, NWB
                KT = KTWB(JW)
                Do JB=BS(JW),BE(JW)
                    IU = CUS(JB)
                    ID = DS(JB)
                    Do SegNumI = IU, ID
                        If(.not. BedConsolidationSeg(SegNumI))Cycle  ! cb 7/5/18
                        z(SegNumI) = z(SegNumI) + SEChange
                    End Do !SegNumI
                End Do !JB
            End Do !JW
            
            !Update other variables
            Do JW=1, NWB
                KT = KTWB(JW)
                Do JB=BS(JW),BE(JW)
                    IU = CUS(JB)
                    ID = DS(JB)
                    Do SegNumI = IU, ID
                        If(.not. BedConsolidationSeg(SegNumI))Cycle  ! cb 7/5/18
                        SZ(SegNumI) = Z(SegNumI)
                        
                        Do K = KT, KB(SegNumI)
                            If(K == KT)H1(K,SegNumI)            =   H(K,JW)-Z(SegNumI)
                            If(K /= KT)H1(K,SegNumI)            =   H(K,JW)
                            BH1(K,SegNumI)                      =   B(K,SegNumI)*H1(K,SegNumI)
                            BH2(K,SegNumI)                      =   BH1(K,SegNumI)
                            BHR1(K,SegNumI)                     =   B(K,SegNumI)*H1(K,SegNumI)
                            BHR2(K,SegNumI)                     =   BHR1(K,SegNumI)
                        End Do !K
                    End Do !SegNumI
                End Do !JB
            End Do !JW
            
            !Move KB layer to padded cells
            Do JW=1,NWB
                KT = KTWB(JW)
                Do JB=BS(JW),BE(JW)
                    IU    = US(JB)
                    ID    = DS(JB)
                  Do SegNumI = IU, ID  ! cb 7/5/18
                    If(.not. BedConsolidationSeg(SegNumI))Cycle
                    if(SegNumI == IU)KB(IU-1) = KB(IU)
                    if(SegNumI == ID)KB(ID+1) = KB(ID)
                  end do
                End Do  !JB
            End Do !JW
            KBI = KB
            
            Do JW=1,NWB
                KT = KTWB(JW)
                Do JB=BS(JW),BE(JW)
                    IU    = US(JB)
                    ID    = DS(JB)
                    Do I=IU-1,ID   
                        If(.not. BedConsolidationSeg(SegNumI))Cycle  ! cb 7/5/18
                        KBMIN(I) =  MIN(KB(I),KB(I+1))
                    End Do  !I
                    KBMIN(ID+1) = KBMIN(ID)
                End Do  !JB
            End Do !JW
            
            !Compute updated volume
            VolumeUpdated = 0.d00
            Do JW=1, NWB
                KT = KTWB(JW)
                Do JB=BS(JW),BE(JW)
                    IU = CUS(JB)
                    ID = DS(JB)
                    Do SegNumI = IU, ID
                        If(.not. BedConsolidationSeg(SegNumI))Cycle  ! cb 7/5/18
                        Do K = KT, KB(SegNumI)
                            If(K == KT)VolumeUpdated = VolumeUpdated + B(K,SegNumI)*DLX(SegNumI)*(H(K,JW)-Z(SegNumI))
                            If(K /= KT)VolumeUpdated = VolumeUpdated + B(K,SegNumI)*DLX(SegNumI)*H(K,JW)
                        End Do !K
                    End Do !SegNumI
                End Do !JB
            End Do !JW
            
            DO JW=1,NWB
                DO JB=BS(JW),BE(JW)
                    DO I=CUS(JB)-1,DS(JB)
                        If(.not. BedConsolidationSeg(SegNumI))Cycle  ! cb 7/5/18
                        DEPTHB(KTWB(JW),I) = H1(KTWB(JW),I)
                        if(kbi(i) < kb(i))depthb(ktwb(jw),i)=(h1(ktwb(jw),i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))    ! SW 1/23/06
                        DEPTHM(KTWB(JW),I) = H1(KTWB(JW),I)*0.5
                        if(kbi(i) < kb(i))depthm(ktwb(jw),i)=(h1(ktwb(jw),i)-(el(kbi(i)+1,i)-el(kb(i)+1,i)))*0.5    ! SW 1/23/06
                        DO K=KTWB(JW)+1,KMX
                            DEPTHB(K,I) = DEPTHB(K-1,I)+ H1(K,I)
                            DEPTHM(K,I) = DEPTHM(K-1,I)+(H1(K-1,I)+H1(K,I))*0.5
                        END DO
                    END DO
                END DO
            END DO
            
            DO JW=1,NWB
                DO I=1,NISNP(JW)
                    If(.not. BedConsolidationSeg(SegNumI))Cycle  ! cb 7/5/18
                    KBR(JW) = MAX(KB(ISNP(I,JW)),KBR(JW))
                END DO
            END DO
            
            !SP added to test 04/01/2010
            DO JW=1,NWB
                CALL INTERPOLATION_MULTIPLIERS
        END DO  
        End If
        
        If(BottomUpdated .and. MoveFFTLayerDown)Call MoveFFTLayerConsolid
    
    Return
    
    Entry ComputeCEMARelatedSourceSinks
        
        !VOLCEMA = 0.d00
        Do JW=1, NWB
            KT = KTWB(JW)
            Do JB=BS(JW),BE(JW)
                IU = CUS(JB)
                ID = DS(JB)
                Do SegNumI = IU, ID
                
                    If(NumCEMAPWInst(SegNumI) > 5)Then
                        ApplyCEMAPWRelease(SegNumI) = .FALSE.
                        NumCEMAPWInst(SegNumI) = 0
                    End If
                    
                    If(CEMALayerAdded(SegNumI))Then
                        CEMACumPWToRelease(SegNumI) = CEMACumPWRelease(SegNumI)
                        CEMACumPWReleased(SegNumI) = CEMACumPWToRelease(SegNumI)/5.d0 !Release over 100 installments
                        ApplyCEMAPWRelease(SegNumI) = .TRUE.
                        CEMASSApplied(SegNumI) = .TRUE.
                    End If
                    
                    If(ApplyCEMAPWRelease(SegNumI))Then
                        NumCEMAPWInst(SegNumI) = NumCEMAPWInst(SegNumI) + 1
                        QSS(KB(SegNumI),SegNumI) = QSS(KB(SegNumI),SegNumI) + CEMACumPWReleased(SegNumI)/dlt !m�/s
                    TSS(KB(SegNumI),SegNumI) = TSS(KB(SegNumI),SegNumI) + CEMACumPWReleased(SegNumI)/dlt*SD_T1  !Cm�/s      ! cb 5/22/15
                        VOLCEMA(JB) = VOLCEMA(JB) + CEMACumPWRelease(SegNumI)
                    End If
                    
                End Do !SegNumI
            End Do !JB
        End Do !JW
        
    Return
    
    RETURN
    
    End Subroutine
    
    SUBROUTINE CEMASedimentProcesses
     USE SEDMODEL_LOCAL
        Do JW=1, NWB
            KT = KTWB(JW)
            Do JB=BS(JW),BE(JW)
                IU = CUS(JB)
                ID = DS(JB)
                Do SegNumI = IU, ID
                    If(EndBedConsolidation(SegNumI))Cycle
                    
                    !Calculate scouring
                    Call CEMACalculateScour
                    
                    !Calculate settling
                    If(IncludeFFTLayer .and. FFTActive)Then
                        SSettle = FFTLayerSettVel*CEMASedConc(SegNumI, KB(SegNumI)) 
                    End If
                    If(IncludeFFTLayer .and. .NOT. FFTActive)Then
                        SSettle = 0.d00
                    End If
                    If(.NOT. IncludeFFTLayer)Then
                        SSettle = CEMASedimentSVelocity*CEMASedConc(SegNumI, KB(SegNumI))
                    End If
                    
                    !Net Flux
                    SNetFlx = SScour - SSettle  !gm/m2/s
                    
                    !If(SNetFlx > 1)SNetFlx = 1
                    !If(SNetFlx < -1)SNetFlx = -1
                    
                    !New sediment concentration
                    CEMASedConc(SegNumI, KB(SegNumI)) = CEMASedConc(SegNumI, KB(SegNumI)) + SNetFlx/H(KB(SegNumI),JW)*dlt
                    
                BedElevation(SegNumI) = BedElevation(SegNumI) - SNetFlx*dlt/(CEMASedimentDensity*1000.d0)   !Density from kg/m?to gm/m?
                BedElevationLayer(SegNumI) = BedElevationLayer(SegNumI) - SNetFlx*dlt/(CEMASedimentDensity*1000.d0)   !Density from kg/m?to gm/m?
                    
                End Do !SegNumI
            End Do !JB
        End Do !JW
    
    Return
    END SUBROUTINE
    SUBROUTINE CEMACalculateScour
    USE SEDMODEL_LOCAL
        Sgr = CEMASedimentDensity/1000.d0
		Dks = CEMAParticleSize
		SScour = 0.d00
    
        !Calculate Erosion
        If(CEMASedimentType == 1)Then   !Cohesive Sediments
        	! Based on formulation in EFDC
			! Sediment Scour for Cohesive Sediments
			CritShldPar = 0.1
			
			! Estimate critical shear stress for erosion
			TCrit_E = 0.d0
			TCrit_S = 0.d0
			!Surface Erosion    !LB SW 3/2019 Mislabeled !Mass Erosion
			If(CEMASedimentDensity > 1065.)then	
				TCrit_S = (CEMASedimentDensity/1000. - 1.065)**0.2        ! SW 3/2019
				TCrit_S = 0.001*(0.883*TCrit_S + 0.05) !Bulk Density divided by 1000 to convert unit to gm/cm3
			End If		
			!Mass Erosion     !LB SW 3/2019 Mislabeled   !Surface Erosion
			If(CEMASedimentDensity > 1013.)then		
				TCrit_E = 0.001*(9.808*CEMASedimentDensity/1000.0 - 9.934) !Bulk Density divided by 1000 to convert unit to gm/cm3
			End If

			VelTemp = 0.5*(u(KB(SegNumI),SegNumI) + u(KB(SegNumI-1),SegNumI-1))
			VelTemp = sqrt(VelTemp**2)
			!Limit bottom velocity
			VelTemp = Min(VelTemp,0.01)
			IF(MANNINGS_N(JW))THEN
                HRAD = BH2(KB(SegNumI),SegNumI)/(BI(KB(SegNumI),SegNumI)+2.0*H2(KB(SegNumI),SegNumI))    !SW 3/2019
                VelTemp = sqrt(G)*VelTemp*FRIC(SegNumI)/HRAD**0.1667    !SW 3/2019
            ELSE
                VelTemp = sqrt(G)*VelTemp/FRIC(SegNumI)     !SW 3/2019
            ENDIF
            
    
			TFlow_B = RHO(KB(SegNumI),SegNumI)*(VelTemp**2)

			If(TFlow_B > TCrit_E)then              !If(TFlow_B > TCrit_S)then   ! Code fix LB SW 3/2019
				!Mass Erosion
				SScour = 0.62    ! g/m2/s
				!SScour = MAX(SScour,0.0)   ! Deleted code not necessary SW 3/2019
			Else
				If(TFlow_B > TCrit_S)then       !If(TFlow_B > TCrit_E)then      ! Code fix LB SW 3/2019
					!Surface Erosion 
					TempVariable = 0.23*dexp(0.198/(CEMASedimentDensity/1000.0 - 1.0023))
					TempVariable = (1/360.0)*(10**TempVariable)	!1/360 to convert from mg/hr-cm2 to gm/s-m2
					SScour = (TFlow_B - TCrit_S)/TCrit_S          !(TFlow_B - TCrit_E)/TCrit_E
					SScour = (SScour)**CritShldPar
					SScour = TempVariable*SScour
					SScour = dmax1(SScour,0.0)
				Else
					!No Erosion
					SScour = 0.d0
				End If
			End If   
        End If
        
        If(CEMASedimentType == 2)Then   !Non-Cohesive Sediments
            	NuTemp = 1.79e-6*dexp(-0.0266*0.5*(T1(KB(SegNumI),SegNumI)))
		        YalinP = ((CEMASedimentDensity - 1000.0)*G*CEMAParticleSize**3.0/(1000.0*NuTemp**2.0))**0.5
                IF(YALINP<=100.)THEN
                    Value1 = 0.041*(Log10(YalinP))**2.0 - 0.356*Log10(YalinP) - 0.977
                    CEMAShieldsNumber = (10.0)**Value1
                ELSEIF(YALINP>100. .AND. YALINP<=3000.)THEN
                    Value1 = 0.132*Log10(YalinP) - 1.804 
			        CEMAShieldsNumber = (10.0)**Value1
                ELSE
                    CEMAShieldsNumber = 0.045
                ENDIF                       
			Uks = Sqrt(CEMAShieldsNumber*(Sgr - 1.0)*G*Dks)	!Critical shear velocity (Shield's diagram)
			NuTemp = 1.79e-6*dexp(-0.0266*0.5*(T1(KB(SegNumI),SegNumI)))
			
			Sum1 = 0.d00
			kt = ktwb(jw)
			d1 = H(K,JW)-Z(SegNumI)
			Depth = 0.d00
			Do k = kt, KB(SegNumI)
				VelTemp = 0.5*(u(K,SegNumI) + u(K,SegNumI-1))
			    VelTemp = dsqrt(VelTemp**2)
				
				Sum1 = Sum1 + VelTemp*d1
				Depth = Depth + d1
				d1 = H(K+1,JW)
			End Do !k
			
			rgha = 0.01*Depth
			rghapda = rgha + 0.4*rgha
			
			UctAvg = Sum1/Depth				!Depth Averaged velocity Pg. 34
			!Limit bottom velocity
			UctAvg = Min(UctAvg,0.01)    ! This is a vertical average - probably not valid in a deep reservoir or lake but approximate order of magnitude
			If(Mannings_N(JW)) Then
                !HRAD = BH1(KT,I)/(B(KTI(I),I)-B(KT+1,I)+2.*AVH1(KT,I))
                HRAD = BH2(KB(SegNumI),SegNumI)/(BI(KB(SegNumI),SegNumI)+2.0*H2(KB(SegNumI),SegNumI))    !SW 3/2019
                !gc2=g*fric(i)*fric(i)/hrad**0.33333333
                gc2=g*fric(SegNumI)*fric(SegNumI)/hrad**0.33333333                  !SW 3/2019
            Else
                GC2 = 0.0
                !If(FRIC(I) /= 0.0) GC2 = G/(FRIC(I)*FRIC(I))
                If(FRIC(SegNumI) /= 0.0) GC2 = G/(FRIC(SegNumI)*FRIC(SegNumI))      !SW 3/2019
            End If
            Usr = Sqrt(GC2)*UctAvg 			
			!
			!Erosion process
			Dsks = Dks*((Sgr - 1.0)*G/NuTemp**2.0)**(1.0/3.0)		!Dimensionless particle diameter Pg. 34
			Tsks = (Usr**2.0 - Uks**2.0)/(Uks**2.0)				    !Transport stage parameter Pg. 34
			Wf = CEMASedimentSVelocity
			SRho = CEMASedimentDensity
			
			If(Wf/Usr >= 0.1 .and. Wf/Usr < 1.0) Then
				BetaD = 1.0 + 2.0*(CEMASedimentSVelocity/Usr)**2.0		!Equation 68 Page 39
			Else
				BetaD = 1.0
			End If
			SDfz = BetaD*az(KB(SegNumI)-1,SegNumI)		    !Vertical mass diffusion coefficient Equation 67 Page 39
			Caks = 0.0     !1.0e+20				! Near Bed Concentration Equation 54 Page 34	
			If(Tsks > 0.0) Caks = 0.015*(Dks/rgha)*(Tsks**1.5/Dsks**0.3)*SRho*1000.0 !gm/m^3
			!
			!Near Bed concentration extrapolated from suspended sediment Equation 31 Page 21
			!Linear interpolation Ref (
			d1 = H(KB(SegNumI),JW)
			d2 = H(KB(SegNumI)-1,JW)
			SedECoeff = (d1 + 0.5*d2 - rghapda)/(0.5*(d1 + d2)) 
			
			Capda = (1.0 - SedECoeff)*CEMASedConc(SegNumI, KB(SegNumI)) + SedECoeff*CEMASedConc(SegNumI, KB(SegNumI)-1)

			If(Capda < 0)Capda = 0.d0
			!
			!Erosion Equation 12 Page 9
			If(Caks > Capda) SScour = -BetaD*SDfz*((Capda - Caks)/(rghapda - rgha))

        End If
        
    Return
    END SUBROUTINE