  !8/2020: Add LayerNum
  !
  Subroutine GasBubblesFormation(Radius,DeltaT, Volume)
    
    Use GLOBAL
    Use SCREENC
    Use CEMAVars
    Use CEMASedimentDiagenesis, only: SD_T1         
    
    ! Type declarations
    IMPLICIT NONE
    
    Integer nGas, nRelArr, temp   ! nTry
    Real(8), Allocatable, Dimension(:) :: Ctot, CgB, C0B, C1B, Source0
    Real(8), Allocatable, Dimension(:) :: K, Mw, Henry
    Real(8) Volume, Porosity
    Real(8) DeltaT,Ro, CgT
    Real(8) RSI, P0, Nbubbles, NbubblesP, NbubbLost 
    Real(8) Pcrit, Pbubb, PbubbT
    Real(8) NetMass, DisMass, GasMass
    Real(8) C1T, C0T, CtT, Radius, Vbub, SourceT
    Real(8) Vbubbles, DiffVolume, BubSedT
    Real(8) Source
    Logical FoundOpenArray
        
    Allocate(Ctot(NumGas), CgB(NumGas), C0B(NumGas), C1B(NumGas))
    Allocate(Source0(NumGas), Henry(NumGas), K(NumGas), Mw(NumGas))
    
    Porosity    = BedPorosity(SegNumI)
    Henry(1)    = HenryConst_H2S        !L atm/M  H2S
    Henry(2)    = HenryConst_CH4        !L atm/M  CH4
    Henry(3)    = HenryConst_NH3        !L atm/M  NH3
    Henry(4)    = HenryConst_CO2        !L atm/M  CO2
    BubSedT     = 273.15 + SD_T1        !K           ! cb 5/22/15
    K(1)        = Henry(1)/GasConst_R/BubSedT
    K(2)        = Henry(2)/GasConst_R/BubSedT
    K(3)        = Henry(3)/GasConst_R/BubSedT
    K(4)        = Henry(4)/GasConst_R/BubSedT
    Ro          = 0.d0      !m
    RSI         = 8318.78   !l-N/m�/mol/K
    P0          = 9800.      !N/m�
    Mw(1)       = 36.        !H2S gm/mol
    Mw(2)       = 16.        !CH4 gm/mol
    Mw(3)       = 17.        !NH3 gm/mol
    Mw(4)       = 44.        !CO2 gm/mol
    If(CrackOpen(SegNumI))NbubbLost = MFTBubbReleased(SegNumI)
    
    Do nGas = 1, NumGas
        Ctot(nGas)    = TConc(nGas,LayerNum, SegNumI)
        Source0(nGas) = SConc(nGas,LayerNum, SegNumI)
    End Do !nGas
    
    If(FirstTimeInBubbles)Then
        
        CgT = 0.d0
        C1T = 0.d0
        C0T = 0.d0
        CtT = 0.d0
        Do nGas = 1, NumGas
            C0B(nGas) = Ctot(nGas)/(1+K(nGas))
            CgB(nGas) = C0B(nGas)*K(nGas)
            C1B(nGas) = Ctot(nGas)
            CgT = CgT + CgB(nGas)
            C1T = C1T + C1B(nGas)
            C0T = C0T + C0B(nGas)
            CtT = CtT + Ctot(nGas)
        End Do !nGas
        Radius = sqrt(2.0*Porosity*GasDiff_Sed*DeltaT*(C1T-C0T)/CgT + Ro**2)
        Vbub = (4.0/3.0)*3.1415927*Radius**3
        NetMass = CtT*Volume*Porosity
        DisMass = C0T*Volume*Porosity
        GasMass = NetMass - DisMass
        Nbubbles = GasMass/(Vbub*CgT)
        Pcrit = 1.32*(CritStressIF**6/(YoungModulus*Nbubbles*Vbub))**0.2 + P0
        PbubbT = 0.d0
        Do nGas = 1, NumGas
            Pbubb = CgB(nGas)*RSI*0.001*BubSedT/Mw(nGas)
            PbubbT = PbubbT + Pbubb
        End Do !nGas
        
    Else
        
        SourceT = 0.d0
        CgT = 0.d0
        C1T = 0.d0
        C0T = 0.d0
        CtT = 0.d0
        Do nGas = 1, NumGas
    
            Source = Source0(nGas)
            SourceT = SourceT + Source
            Ctot(nGas) = Ctot(nGas) + Source*DeltaT
            C0B(nGas) = Ctot(nGas)/(1+K(nGas))
            CgB(nGas) = C0B(nGas)*K(nGas)
            C1B(nGas) = C0B(nGas)
            CgT = CgT + CgB(nGas)
            C1T = C1T + C1B(nGas)
            C0T = C0T + C0B(nGas)
            CtT = CtT + Ctot(nGas)
            
        End Do !nGas
        
        Radius = Radius + Porosity*GasDiff_Sed/(Radius*CgT)*(SourceT*CalibParam_R1**2/(6*GasDiff_Sed)+(C1T-C0T))*DeltaT
            
        If(LimBubbSize)Then
            If(Radius > MaxBubbRad/1000.0)Radius = MaxBubbRad/1000.0
        End If
    
    End If
    
    Vbub = (4.0/3.0)*3.1415927*Radius**3
    NetMass = 0.d0
    DisMass = 0.d0
    CgT = 0.d0
    Do nGas = 1, NumGas
        NetMass = NetMass + Ctot(nGas)*Volume*Porosity
        DisMass = DisMass + C0B(nGas)*Volume*Porosity
        CgT = CgT + CgB(nGas)
    End Do !nGas
    GasMass = NetMass - DisMass
    Nbubbles = GasMass/(Vbub*CgT)
    
    Pcrit = 1.324*(CritStressIF**6/(YoungModulus*Nbubbles*Vbub))**0.2 + P0
    PbubbT = 0.d0
    Do nGas = 1, NumGas
        Pbubb = CgB(nGas)*RSI*0.001*BubSedT/Mw(nGas)
        PbubbT = PbubbT + Pbubb
    End Do !nGas
    
    PresBubbSed(SegNumI) = PbubbT
    PresCritSed(SegNumI) = Pcrit
        
    If(PbubbT < Pcrit*CrackCloseFraction)Then
        CrackOpen(SegNumI) = .FALSE.
        NbubbLost = 0
    End If
    
    If(LastDiffVolume(SegNumI) < 0)LastDiffVolume(SegNumI) = 0.d00
        
    If(PbubbT > Pcrit .and. .NOT. CrackOpen(SegNumI))Then
        
        CrackOpen(SegNumI) = .TRUE.
        Vbubbles = CritStressIF**6/(YoungModulus*((PbubbT-P0)/1.32)**5)
        DiffVolume = Vbub*Nbubbles - Vbubbles
        DiffVolume = DiffVolume*BubbRelScale
        If(UseReleaseFraction)Then
            Vbubbles = CritStressIF**6/(YoungModulus*((Pcrit*CrackCloseFraction)/1.32)**5)
            DiffVolume = BubbRelFraction*(Vbub*Nbubbles - Vbubbles)
        End If
        LastDiffVolume(SegNumI) = DiffVolume
        NbubbLost = (DiffVolume)/Vbub
        NbubblesP = Nbubbles
        Nbubbles = Nbubbles - NbubbLost
        Do nGas = 1, NumGas
            CgB(nGas) = CgB(nGas)*(Vbub*NbubblesP - DiffVolume)/(Vbub*NbubblesP)
            C0B(nGas) = CgB(nGas)/K(nGas)
            Ctot(nGas) = C0B(nGas)*(1+K(nGas))
        End Do !nGas
        
    End If
    
!!!!!!!!!!!!!!!!! debug    
!    CrackOpen(SegNumI)=.false.
!!!!!!!!!!!!!!!!!!!!!!! debug
    If(CrackOpen(SegNumI))Then
    
        NbubblesP = Nbubbles
        !Nbubbles = Nbubbles - NbubbLost  ! cb 2/21/13
        Do nGas = 1, NumGas
            CgB(nGas) = CgB(nGas)*(Vbub*NbubblesP - LastDiffVolume(SegNumI))/(Vbub*NbubblesP)
            C0B(nGas) = CgB(nGas)/K(nGas)
            Ctot(nGas) = C0B(nGas)*(1+K(nGas))
            
            TConcP(nGas,LayerNum,SegNumI) = Ctot(nGas)
             TConc(nGas,LayerNum,SegNumI) = Ctot(nGas)
    
        End Do !nGas
        
        temp = INT4(NbubbLost)
    
        FoundOpenArray = .FALSE.
        MFTBubbReleased(SegNumI) = KIDINT(NbubbLost)
        !nTry = 0
        Do nRelArr = 1, NumBubRelArr
            !nTry = nTry + 1
            If(BubblesStatus(SegNumI, nRelArr) == 0)Then
                FoundOpenArray = .TRUE.
                BubblesStatus(SegNumI, nRelArr) = 1
                BubblesCarried(SegNumI, nRelArr) = MFTBubbReleased(SegNumI)
                BubblesRadius(SegNumI, nRelArr) = Radius
                BubblesLNumber(SegNumI, nRelArr) = KB(SegNumI)
                Do nGas = 1, NumGas
                    BubblesGasConc(SegNumI, nRelArr, nGas) = CgB(nGas)
                End Do
                Exit
            End If
        End Do
        If(.NOT. FoundOpenArray)Then
            Write(CEMALogFilN,*)"Insufficient array size for bubbles release at JDAY = ", JDAY
            Write(w2err,*)"Insufficient array size for bubbles release at JDAY = ", JDAY
            Stop
        End if
    End If
    
    CgSed(SegNumI) = CgT
    C0Sed(SegNumI) = C0T
    CtSed(SegNumI) = CtT
     
    Return
  End Subroutine

  Subroutine CEMACalculateRiseVelocity

    Use MAIN 
    Use GLOBAL
    Use GEOMC
    Use CEMAVars
    IMPLICIT NONE
    
    Integer nGas, nRelArr
    Real(8) Rhog
    
    Do JW=1, NWB
        KT = KTWB(JW)
        Do JB=BS(JW),BE(JW)
            IU = CUS(JB)
            ID = DS(JB)
            Do SegNumI = IU, ID
                Do nRelArr = 1, NumBubRelArr
                    If(BubblesStatus(SegNumI, nRelArr) == 0)Then    
                        
                        Do nGas = 1, NumGas
                            BRVoluAGas(SegNumI, nRelArr, nGas) = 0.d00
                            BRRateAGas(SegNumI, nRelArr, nGas) = 0.d00
                            BRRateAGasNet(SegNumI, nGas) = 0.d00
                        End Do !nGas  
                        
                    End If     
                End Do                    
            End Do                    
        End Do                
    End Do            

    
    Do JW=1, NWB
        KT = KTWB(JW)
        Do JB=BS(JW),BE(JW)
            IU = CUS(JB)
            ID = DS(JB)
            Do SegNumI = IU, ID
                Do nRelArr = 1, NumBubRelArr
                    If(BubblesStatus(SegNumI, nRelArr) == 1)Then
                        
                        Rhog = 0.d0
                        Do nGas = 1, NumGas
                            Rhog = Rhog + BubblesGasConc(SegNumI, nRelArr, nGas)/1000.0          !kg/m�
                        End Do !nGas    
                        
                        Call CEMABubblesRiseVelocity(BubblesRadius(SegNumI, nRelArr), Rhog, BubblesRiseV(SegNumI, nRelArr))
                        
                    End If
                End Do                    
            End Do                    
        End Do                
    End Do                

    Return
  End Subroutine

  Subroutine CEMABubblesRiseVelocity(Radius, Rhog, RiseVelocity)

    IMPLICIT NONE

    Real(8) Radius, Rhog, RiseVelocity
    Real(8) Rhow, DynVisc
    Real(8) Sigma, Nd, W, Reynolds, M, Eo
    Real(8) J, H
    Real(8) Radius1, RiseVelocity1, RiseVelocity2

    !Calculate Rise Velocity
    Rhow = 1000             !kg/m3
    DynVisc = 0.001002      !kg/m/s
    Sigma = 0.0725          !N/m
    
    If(Radius*1000 <= 1)Then !<= 1 mm
        Nd = 4.*Rhow*(Rhow-Rhog)*9.8*Radius**3/(3.*DynVisc**2)
        W = dlog10(Nd)
        If(Nd <= 73.)Then
            Reynolds = Nd/24. - 1.7569d-4*Nd**2 + 6.9252d-7*Nd**3 - 2.3027d-10*Nd**4
        End If
        If(Nd > 73.0 .and. Nd <= 580.0)Then
            Reynolds = 10**(-1.7095 + 1.33438*W - 0.11591*W**2)
        End If
        If(Nd > 580.)Then
            Reynolds = 10**(-1.81391 + 1.34671*W - 0.12427*W**2 + 0.006344*W**3)
        End If
        RiseVelocity = Reynolds*DynVisc/(Rhow*Radius)
    End If
    
    If(Radius*1000.0 <= 15.0 .and. Radius*1000.0 > 1.0)Then !<= 15 mm
        M = 9.8*DynVisc**4*(Rhow-Rhog)/(Rhow**2*Sigma**3)
        Eo = 9.8*(Rhow-Rhog)*Radius**2/Sigma
        H = (4./3.)*Eo*M**(-0.149)*(DynVisc/DynVisc)**(-0.14)
        If(H < 59.3)Then
            J = 0.94*H**0.757
        Else
            J = 3.42*H**0.441
        End If
        RiseVelocity = DynVisc/(Rhow*Radius)*M**(-0.149)*(J-0.857)
    End If
    
    If(Radius*1000.0 > 15.0 .and. Radius*1000.0 <= 18.0)Then !<= 15 mm to 18 mm
        
        Radius1 = Radius
        Radius = 0.015
        M = 9.8*DynVisc**4*(Rhow-Rhog)/(Rhow**2*Sigma**3)
        Eo = 9.8*(Rhow-Rhog)*Radius**2/Sigma
        H = (4./3.)*Eo*M**(-0.149)*(DynVisc/DynVisc)**(-0.14)
        If(H < 59.3)Then
            J = 0.94*H**0.757
        Else
            J = 3.42*H**0.441
        End If
        RiseVelocity1 = DynVisc/(Rhow*Radius)*M**(-0.149)*(J-0.857)
        
        Radius = 0.018
        RiseVelocity2 = 0.711*sqrt(9.8*Radius*(Rhow-Rhog)/Rhow)
        
        Radius = Radius1
        
        RiseVelocity = ((Radius - 0.015)*RiseVelocity2 + RiseVelocity1*(0.018 - Radius))/(0.018-0.015)
        
    End If
    
    If(Radius*1000.0 > 18.0)Then !> 18 mm
        RiseVelocity = 0.711*sqrt(9.8*Radius*(Rhow-Rhog)/Rhow)
    End If

    Return
  End Subroutine


  Subroutine CEMABubblesTransport

    Use MAIN 
    Use GLOBAL
    Use GEOMC
    Use CEMAVars
    IMPLICIT NONE
    
    Integer BubbLayer, nRelArr
    Real(8) VLocationBubble, VDistTravBubble
    
    Do JW=1, NWB
        KT = KTWB(JW)
        Do JB=BS(JW),BE(JW)
            IU = CUS(JB)
            ID = DS(JB)
            Do SegNumI = IU, ID
                Do nRelArr = 1, NumBubRelArr
                    If(BubblesStatus(SegNumI, nRelArr) == 1 .and. .NOT. BubblesAtSurface(SegNumI, nRelArr))Then
                        
                        BubbLayer = BubblesLNumber(SegNumI, nRelArr)
                        VLocationBubble = 0.5*(el(BubbLayer+1, SegNumI) + el(BubbLayer, SegNumI))
                        VDistTravBubble = BubblesRiseV(SegNumI, nRelArr)*dlt
                        VLocationBubble = VLocationBubble + VDistTravBubble
                        
                        !Locate vertical location
                        BubblesLNumber(SegNumI, nRelArr) = KT
                        Do K = KT, KB(SegNumI)
                            If(VLocationBubble < el(K,SegNumI))BubblesLNumber(SegNumI, nRelArr) = K
                        End Do !K
                        
                        If(BubblesLNumber(SegNumI, nRelArr) == KT)Then
                            BubblesAtSurface(SegNumI, nRelArr) = .TRUE. 
                            FirstBubblesRelease(SegNumI, nRelArr) = .TRUE.
                            BubblesReleaseAllValue(SegNumI, nRelArr) = BubbRelFractionAtm*BubblesCarried(SegNumI, nRelArr)
                        End If
                        
                    End If
                End Do                    
            End Do                    
        End Do                
    End Do            

    Return
  End Subroutine

  Subroutine CEMABubblesRelease

    Use MAIN 
    Use GLOBAL
    Use GEOMC
    Use CEMAVars
    Use Screenc
    !IMPLICIT NONE    
    Integer nRelArr
    Real(8) TempBubblesRelVolume

    
    BRRateAGasNet = 0.d00
    
    Do JW=1, NWB
        KT = KTWB(JW)
        Do JB=BS(JW),BE(JW)
            IU = CUS(JB)
            ID = DS(JB)
            Do SegNumI = IU, ID
                Do nRelArr = 1, NumBubRelArr
                    If(BubblesStatus(SegNumI, nRelArr) == 1 .and. BubblesAtSurface(SegNumI, nRelArr) .and. .NOT. ICE(SegNumI))Then    
                        
                        Do nGas = 1, NumGas
                            !TempBubblesRelVolume = 4/3*3.14*BubblesRadius(SegNumI, nRelArr)**3
                            TempBubblesRelVolume = 4./3.*3.14*BubblesRadius(SegNumI, nRelArr)**3    ! SW 10/10/2017
                            BRVoluAGas(SegNumI, nRelArr, nGas) = BubblesReleaseAllValue(SegNumI, nRelArr)*TempBubblesRelVolume*BubblesGasConc(SegNumI, nRelArr, nGas)   !gm
                            BRRateAGas(SegNumI, nRelArr, nGas) = BRVoluAGas(SegNumI, nRelArr, nGas)/dlt !gm/s
                            BRRateAGasNet(SegNumI, nGas) = BRRateAGasNet(SegNumI, nGas) + BRRateAGas(SegNumI, nRelArr, nGas) !gm/s
                            BubbleRelWB(JW, nGas)= BubbleRelWB(JW, nGas)+DLT*BRRateAGasNet(SegNumI, nGas)/1000.   ! SW 7/1/2017 Convert from gm/s to kg
                        End Do !nGas  
                        
                        BubblesCarried(SegNumI, nRelArr) = BubblesCarried(SegNumI, nRelArr) - BubblesReleaseAllValue(SegNumI, nRelArr)
                        If(BubblesCarried(SegNumI, nRelArr) <= 0.d00)Then
                            BubblesReleaseAllValue(SegNumI, nRelArr) = 0.d00
                            BubblesCarried(SegNumI, nRelArr) = 0.d00
                            BubblesGasConc(SegNumI, nRelArr,:) = 0.d00
                            BubblesAtSurface(SegNumI, nRelArr) = .FALSE.
                            BubblesStatus(SegNumI, nRelArr) = 0
                        End If                           
                        
                    End If     
                End Do                    
            End Do                    
        End Do                
    End Do 
    
    Return
  End Subroutine


  Subroutine CEMABubbWatTransfer

    Use MAIN 
    Use GLOBAL
    Use GEOMC
    Use SCREENC
    Use KINETIC
    Use CEMAVars
    Use CEMASedimentDiagenesis, only: SD_T1
    IMPLICIT NONE
    Real(8), Allocatable, Dimension(:) :: KValue, Mw, Henry
    real(8) :: BubbDissSrcSnk, BubSedT, EqbDissConcentration
    Integer BubbLNumber, nRelArr,ngasconst, ngas  
    
    Allocate(Henry(NumGas), KValue(NumGas), Mw(NumGas))
    
    Henry(1)    = HenryConst_H2S        !L atm/M  H2S
    Henry(2)    = HenryConst_CH4        !L atm/M  CH4
    Henry(3)    = HenryConst_NH3        !L atm/M  NH3
    Henry(4)    = HenryConst_CO2        !L atm/M  CO2
    BubSedT     = 273.15 + SD_T1        !K
    Mw(1)       = 36.                    !H2S gm/mol
    Mw(2)       = 16.                    !CH4 gm/mol
    Mw(3)       = 17.                    !NH3 gm/mol
    Mw(4)       = 44.                    !CO2 gm/mol
    
    Do JW=1, NWB
        KT = KTWB(JW)
        Do JB=BS(JW),BE(JW)
            IU = CUS(JB)
            ID = DS(JB)
            Do SegNumI = IU, ID
                Do nRelArr = 1, NumBubRelArr
                    If(BubblesStatus(SegNumI, nRelArr) == 1)Then    
                        
                        Do nGas = 1, NumGas
!                        Do nGas = 1, NumGas-1 ! debug
                            
                            if(ngas == 1) ngasconst = NH2S   ! cb 2/18/13
                            if(ngas == 2) ngasconst = NCH4
                            if(ngas == 3) ngasconst = NSO4
                            if(ngas == 4) ngasconst = NTIC
                            
                            BubbLNumber = BubblesLNumber(SegNumI, nRelArr)
                            !KValue(1)        = Henry(1)/GasConst_R/(T1(BubbLNumber,SegNumI) + 273.15)
                            !KValue(2)        = Henry(2)/GasConst_R/(T1(BubbLNumber,SegNumI) + 273.15)
                            !KValue(3)        = Henry(3)/GasConst_R/(T1(BubbLNumber,SegNumI) + 273.15)
                            !KValue(4)        = Henry(4)/GasConst_R/(T1(BubbLNumber,SegNumI) + 273.15)
                            KValue(ngas)        = Henry(ngas)/GasConst_R/(T1(BubbLNumber,SegNumI) + 273.15)
                
                            EqbDissConcentration = BubblesGasConc(SegNumI, nRelArr, nGas)/KValue(nGas)
                            !BubbDissSrcSnk = BubbWatGasExchRate*(EqbDissConcentration - C1(BubbLNumber,SegNumI,1+nGas))    !g/m�/s
                            BubbDissSrcSnk = BubbWatGasExchRate*(EqbDissConcentration - C1(BubbLNumber,SegNumI,ngasconst))    !g/m�/s  cb 2/18/13
                            !CGSS(BubbLNumber,SegNumI,nGasconst) = CGSS(BubbLNumber,SegNumI,nGasconst) + BubbDissSrcSnk    !BubbDissSrcSnk > 0 Bubbles --> Water
                            C1(BubbLNumber,SegNumI,nGasconst) = C1(BubbLNumber,SegNumI,nGasconst) + BubbDissSrcSnk    !BubbDissSrcSnk > 0 Bubbles --> Water
                            BubblesGasConc(SegNumI, nRelArr, nGas) = BubblesGasConc(SegNumI, nRelArr, nGas) - BubbDissSrcSnk*dlt !BubbDissSrcSnk > 0 Bubbles --> Water
                            If(BubblesGasConc(SegNumI, nRelArr, nGas) < 0.d0)BubblesGasConc(SegNumI, nRelArr, nGas) = 0.d0
                            
                        End Do !nGas  
                        
                    End If     
                End Do                    
            End Do                    
        End Do                
    End Do
    
    Return
  End Subroutine

  Subroutine CEMABubblesTurbulence
    Use MAIN 
    Use GLOBAL
    Use GEOMC
    Use SCREENC
    Use KINETIC
    Use CEMAVars
    IMPLICIT NONE    
    Real(8) TempBubbDiam, TempRelVelocity
    Integer BubbCntr, nRelArr, BubbLNumber
    
    SegNumI = I
    Do k = KT, KBMIN(SegNumI)-1
        BubbCntr = 0
        TempBubbDiam = 0
        TempRelVelocity = 0
        Do nRelArr = 1, NumBubRelArr
            If(BubblesStatus(SegNumI, nRelArr) == 1)Then    
                BubbLNumber = BubblesLNumber(SegNumI, nRelArr)
                
                If(BubbLNumber == k)Then
                    BubbCntr = BubbCntr + 1
                    TempBubbDiam = TempBubbDiam + 2.0*BubblesRadius(SegNumI, nRelArr)
                    TempRelVelocity = TempRelVelocity + BubblesRiseV(SegNumI, nRelArr) + W(K-1,SegNumI) !Rise velocity is +ve upwards and W is +ve downwards
                End If
                
            End If
        End Do !nRelArr
        
        If(BubbCntr > 0)Then
            TempBubbDiam = TempBubbDiam/BubbCntr  
            TempRelVelocity = TempRelVelocity/BubbCntr
            If(k == KB(SegNumI))Then
                AZ(K,SegNumI) = AZ(K,SegNumI) + BottomTurbulence(SegNumI)
            Else
                !AZ(K,SegNumI) = AZ(K,SegNumI)  + TempBubbDiam/TempRelVelocity
                AZ(K,SegNumI) = AZ(K,SegNumI)  + TempBubbDiam*TempRelVelocity    ! cb 2/7/13
            End If
        End If
        
        AZ(KB(SegNumI)-1,SegNumI) = AZ(KB(SegNumI)-1,SegNumI) + BottomTurbulence(SegNumI)
        
    End Do                    
    
    Return
  End Subroutine

  Subroutine CEMABubblesReleaseTurbulence
    Use MAIN 
    Use GLOBAL
    Use GEOMC
    Use SCREENC
    Use KINETIC
    Use CEMAVars
    IMPLICIT NONE    
    Real(8) TempBubbDiam, TempRelVelocity
    Integer BubbCntr, nRelArr, BubbLNumber
    
    BottomTurbulence = 0.d00
    Do JW=1, NWB
        k = KB(SegNumI)
        Do JB=BS(JW),BE(JW)
            IU = CUS(JB)
            ID = DS(JB)
            Do SegNumI = IU, ID
    
                BubbCntr = 0
                TempBubbDiam = 0
                TempRelVelocity = 0
                Do nRelArr = 1, NumBubRelArr
                    If(BubblesStatus(SegNumI, nRelArr) == 1)Then    
                        BubbLNumber = BubblesLNumber(SegNumI, nRelArr)
                        
                        If(BubbLNumber == k)Then
                            BubbCntr = BubbCntr + 1
                            TempBubbDiam = TempBubbDiam + 2.0*BubblesRadius(SegNumI, nRelArr)
                            TempRelVelocity = TempRelVelocity + BubblesRiseV(SegNumI, nRelArr) + W(K-1,SegNumI) !Rise velocity is +ve upwards and W is +ve downwards
                        End If
                        
                    End If
                End Do !nRelArr
                
                If(BubbCntr > 0)Then
                    TempBubbDiam = TempBubbDiam/BubbCntr  
                    TempRelVelocity = TempRelVelocity/BubbCntr
                    BottomTurbulence(SegNumI) = CEMATurbulenceScaling*TempBubbDiam/TempRelVelocity
                End If
            
            End Do !SegNumI
        End Do  !JB       
    End Do !JW
    
    Do JW=1, NWB
        k = KB(SegNumI)
        Do JB=BS(JW),BE(JW)
            IU = CUS(JB)
            ID = DS(JB)
            Do SegNumI = IU, ID
                
                If(BottomTurbulence(SegNumI) > 0)Then
                    Continue
                End If    
                
            End Do
        End Do
    End Do
                       
  Return
  End Subroutine
