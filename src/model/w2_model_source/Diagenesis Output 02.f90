Module CEMAOutputRoutines
    Use MAIN 
    Use GLOBAL
    Use SCREENC
    Use GEOMC
    Use KINETIC
    USE NAMESC
    Use CEMAVars 
    USE CEMASedimentDiagenesis
    
    ! Type declarations
    IMPLICIT NONE

    Integer nGas
    
    contains
    !
    !
    Subroutine WriteCEMASedimentModelOutput
    
        !Bed Elevation
        If(WriteBESnp)Then
            Write(CEMASNPOutFilN,'("Bed elevation(m)")')
            Write(CEMASNPOutFilN,'("JDAY = ",f14.6)')JDAY
            Write(CEMASNPOutFilN,'(<IMX>i10)')(SegNumI, SegNumI = 1, IMX)
            Write(CEMASNPOutFilN,'(<IMX>f10.4)')(BedElevation(SegNumI), SegNumI = 1, IMX)
            
            Write(CEMASNPOutFilN,'("Bed elevation layer(m)")')
            Write(CEMASNPOutFilN,'("JDAY = ",f14.6)')JDAY
            Write(CEMASNPOutFilN,'(<IMX>i10)')(SegNumI, SegNumI = 1, IMX)
            Write(CEMASNPOutFilN,'(<IMX>f10.4)')(BedElevationLayer(SegNumI), SegNumI = 1, IMX) 
        End If
        
        !Bed porosity
        If(WritePWSnp)Then
            Write(CEMASNPOutFilN,'("Bed porosity (%)")')
            Write(CEMASNPOutFilN,'("JDAY = ",f8.2)')JDAY
            Write(CEMASNPOutFilN,'(<IMX>i10)')(SegNumI, SegNumI = 1, IMX)
            Write(CEMASNPOutFilN,'(<IMX>f10.4)')(BedPorosity(SegNumI)*100, SegNumI = 1, IMX)
            Write(CEMASNPOutFilN,'("Total Volume of Sediments = ",e14.6," m3")')TotalSedimentsInBed
            Write(CEMASNPOutFilN,'("Total Porewater Volume = ",e14.6," m3")')TotalPoreWatVolume
            Write(CEMASNPOutFilN,'("Total Porewater Removed = ",e14.6," m3")')TotalPoreWatRemoved
!            Write(CEMATSR1OutFilN,'(f14.6 , "," ,i5,",",i5, ",", f14.6, ",", f14.6, ",", f14.6, ",", f14.6, ",", f14.6, ",", f14.6, ",", f14.6, ",", f14.6)')JDAY, KTWB(1), KB(4), z(4), BedElevation(4), BedElevationLayer(4), PorewaterRelRate(4)*86400.d0, BedPorosity(4)*100, BedConsolidRate(4)*86400.d0, EL(KTWB(1),4)-Z(4), CEMACumPWRelease(4)
        End If
        
    End Subroutine
        
    Subroutine WriteCEMASedimentFluxOutput
        
        If(WriteCEMAMFTSedFlx)Then
            
            Write(CEMASedFlxFilN4,'("SOD(gO2/m2/d),",f12.4,",",<IMX>(f10.4,","),<IMX>(f10.4,","))')JDAY,(KFSF(KBI(SegNumI),SegNumI,4),  SegNumI = 1, IMX),(KFSFAV(SegNumI,1), SegNumI = 1, IMX)
            Write(CEMASedFlxFilN5,'("POCG1(g/m3),",f12.4,",",<IMX>(f12.4,","),<IMX>(f12.4,","))')  JDAY,(C2SF(KBI(SegNumI),SegNumI,19), SegNumI = 1, IMX),(KFSFAV(SegNumI,2), SegNumI = 1, IMX)
            Write(CEMASedFlxFilN6,'("POCG2(g/m3),",f12.4,",",<IMX>(f12.4,","),<IMX>(f12.4,","))')  JDAY,(C2SF(KBI(SegNumI),SegNumI,20), SegNumI = 1, IMX),(KFSFAV(SegNumI,3), SegNumI = 1, IMX)
            Write(CEMASedFlxFilN7,'("JC(gC/m2/d),",f12.4,",",<IMX>(f12.4,","),<IMX>(f12.4,","))')           JDAY,(KFSF(KBI(SegNumI),SegNumI,7),  SegNumI = 1, IMX),(KFSFAV(SegNumI,4), SegNumI = 1, IMX)
            Write(CEMASedFlxFilN8,'("JN(gN/m2/d),",f12.4,",",<IMX>(f12.4,","),<IMX>(f12.4,","))')           JDAY,(KFSF(KBI(SegNumI),SegNumI,8),  SegNumI = 1, IMX),(KFSFAV(SegNumI,5), SegNumI = 1, IMX)
            Write(CEMASedFlxFilN9,'("PONG1(g/m3),",f12.4,",",<IMX>(f12.4,","),<IMX>(f12.4,","))')  JDAY,(C2SF(KBI(SegNumI),SegNumI,22), SegNumI = 1, IMX),(KFSFAV(SegNumI,6), SegNumI = 1, IMX)
            Write(CEMASedFlxFilN10,'("PONG2(g/m3),",f12.4,",",<IMX>(f12.4,","),<IMX>(f12.4,","))') JDAY,(C2SF(KBI(SegNumI),SegNumI,23), SegNumI = 1, IMX),(KFSFAV(SegNumI,7), SegNumI = 1, IMX)
            
            Write(CEMASedFlxFilN11,'("JCH4(gO2/m2/d),",f12.4,",",<IMX>(f12.4,","),<IMX>(f12.4,","))')JDAY,(KFSF(KBI(SegNumI),SegNumI,10), SegNumI = 1, IMX),(KFSFAV(SegNumI,8),  SegNumI = 1, IMX)
            Write(CEMASedFlxFilN12,'("JNH4(gN/m2/d),",f12.4,",",<IMX>(f12.4,","),<IMX>(f12.4,","))')JDAY,(KFSF(KBI(SegNumI),SegNumI,14), SegNumI = 1, IMX),(KFSFAV(SegNumI,9),  SegNumI = 1, IMX)
            Write(CEMASedFlxFilN13,'("JNO3(gN/m2/d),",f12.4,",",<IMX>(f12.4,","),<IMX>(f12.4,","))')JDAY,(KFSF(KBI(SegNumI),SegNumI,15), SegNumI = 1, IMX),(KFSFAV(SegNumI,10), SegNumI = 1, IMX)
            Write(CEMASedFlxFilN14,'("JPO4(gP/m2/d), ,",f12.4,",",<IMX>(f12.4,","),<IMX>(f12.4,","))')JDAY,(KFSF(KBI(SegNumI),SegNumI,16), SegNumI = 1, IMX),(KFSFAV(SegNumI,11), SegNumI = 1, IMX)
            
            Write(CEMASedFlxFilN15,'("POPG1(g/m3),",f12.4,",",<IMX>(f12.4,","),<IMX>(f12.4,","))')JDAY,(C2SF(KBI(SegNumI),SegNumI,25), SegNumI = 1, IMX),(KFSFAV(SegNumI,12), SegNumI = 1, IMX)
            Write(CEMASedFlxFilN16,'("POPG2(g/m3),",f12.4,",",<IMX>(f12.4,","),<IMX>(f12.4,","))')JDAY,(C2SF(KBI(SegNumI),SegNumI,26), SegNumI = 1, IMX),(KFSFAV(SegNumI,13), SegNumI = 1, IMX)
            
            Write(CEMASedFlxFilN17,'("CSOD(gO2/m2/d),",f12.4,",",<IMX>(f10.4,","),<IMX>(f10.4,","))')JDAY,(KFSF(KBI(SegNumI),SegNumI,5), SegNumI = 1, IMX),(KFSFAV(SegNumI,14), SegNumI = 1, IMX)
            Write(CEMASedFlxFilN18,'("NSOD(gO2/m2/d),",f12.4,",",<IMX>(f10.4,","),<IMX>(f10.4,","))')JDAY,(KFSF(KBI(SegNumI),SegNumI,6), SegNumI = 1, IMX),(KFSFAV(SegNumI,15), SegNumI = 1, IMX)
            Write(CEMASedFlxFilN19,'("JP(gP/m2/d),",f12.4,",",<IMX>(f10.4,","),<IMX>(f10.4,","))')            JDAY,(KFSF(KBI(SegNumI),SegNumI,9), SegNumI = 1, IMX),(KFSFAV(SegNumI,16), SegNumI = 1, IMX)
            Write(CEMASedFlxFilN20,'("H1(m),",f12.4,",",<IMX>(f10.4,","))')JDAY,(SD_AerLayerThick(SegNumI), SegNumI = 1, IMX)
            Write(CEMASedFlxFilN21,'("Temperature_KB_Layer1(C),",f10.3,",",<IMX>(f10.2,","))')JDAY,(C2SF(KBI(SegNumI),SegNumI,15), SegNumI = 1, IMX)
            Write(CEMASedFlxFilN22,'("Temperature_KB_Layer2(C),",f10.3,",",<IMX>(f10.2,","))')JDAY,(C2SF(KBI(SegNumI),SegNumI,16), SegNumI = 1, IMX)
            
            Write(CEMASedFlxFilN23,'("NO3PorewaterAerobic_KB_Layer1(mgN/L),",f10.3,",",<IMX>(f10.2,","))')JDAY,(C2SF(KBI(SegNumI),SegNumI,3), SegNumI = 1, IMX)
            Write(CEMASedFlxFilN24,'("NO3PorewaterAnaerobic_KB_Layer2(mgN/L),",f10.3,",",<IMX>(f10.2,","))')JDAY,(C2SF(KBI(SegNumI),SegNumI,4), SegNumI = 1, IMX)
            Write(CEMASedFlxFilN25,'("NH3PorewaterAerobic_KB_Layer1(mgN/L),",f10.3,",",<IMX>(f10.2,","))')JDAY,(C2SF(KBI(SegNumI),SegNumI,1), SegNumI = 1, IMX)
            Write(CEMASedFlxFilN26,'("NH3PorewaterAnaerobic_KB_Layer2(mgN/L),",f10.3,",",<IMX>(f10.2,","))')JDAY,(C2SF(KBI(SegNumI),SegNumI,2), SegNumI = 1, IMX)
            Write(CEMASedFlxFilN27,'("PO4PorewaterAerobic_KB_Layer1(mgP/L),",f10.3,",",<IMX>(f10.2,","))')JDAY,(C2SF(KBI(SegNumI),SegNumI,5), SegNumI = 1, IMX)
            Write(CEMASedFlxFilN28,'("PO4PorewaterAnaerobic_KB_Layer2(mgP/L),",f10.3,",",<IMX>(f10.2,","))')JDAY,(C2SF(KBI(SegNumI),SegNumI,6), SegNumI = 1, IMX)
            Write(CEMASedFlxFilN29,'("SO4PorewaterAerobic_KB_Layer1(mgS/L),",f10.3,",",<IMX>(f10.2,","))')JDAY,(C2SF(KBI(SegNumI),SegNumI,9), SegNumI = 1, IMX)
            Write(CEMASedFlxFilN30,'("SO4PorewaterAnaerobic_KB_Layer2(mgS/L),",f10.3,",",<IMX>(f10.2,","))')JDAY,(C2SF(KBI(SegNumI),SegNumI,10), SegNumI = 1, IMX)

            Write(CEMASedFlxFilN31,'("FEIIPorewaterAerobic_KB_Layer1(mgFe/L),",f10.3,",",<IMX>(f10.2,","))')JDAY,(C2SF(KBI(SegNumI),SegNumI,30), SegNumI = 1, IMX)
            Write(CEMASedFlxFilN32,'("FEIIPorewaterAnaerobic_KB_Layer2(mgFe/L),",f10.3,",",<IMX>(f10.2,","))')JDAY,(C2SF(KBI(SegNumI),SegNumI,31), SegNumI = 1, IMX)
            Write(CEMASedFlxFilN33,'("MnIIPorewaterAerobic_KB_Layer1(mgMn/L),",f10.3,",",<IMX>(f10.2,","))')JDAY,(C2SF(KBI(SegNumI),SegNumI,34), SegNumI = 1, IMX)
            Write(CEMASedFlxFilN34,'("MnIIPorewaterAnaerobic_KB_Layer2(mgMn/L),",f10.3,",",<IMX>(f10.2,","))')JDAY,(C2SF(KBI(SegNumI),SegNumI,35), SegNumI = 1, IMX)
            Write(CEMASedFlxFilN35,'("CH4PorewaterAerobic_KB_Layer1(mgO2/L),",f10.3,",",<IMX>(f10.2,","))')JDAY,(C2SF(KBI(SegNumI),SegNumI,7), SegNumI = 1, IMX)
            Write(CEMASedFlxFilN36,'("CH4PorewaterAnaerobic_KB_Layer2(mgO2/L),",f10.3,",",<IMX>(f10.2,","))')JDAY,(C2SF(KBI(SegNumI),SegNumI,8), SegNumI = 1, IMX)

            IF(Bubbles_Calculation) Call CEMATempOutput
	        
        End If
    
    End Subroutine
        
    Subroutine CEMATempOutput
        !
        Implicit none
        Character(3) CrackStatus(IMX)
        !
        CrackStatus = " NO"
        Do SegNumI = 1, IMX
            If(CrackOpen(SegNumI))CrackStatus(SegNumI) = "YES"
        End Do !SegNumI
        
        DO SegNumI=1,IMX
            Write(CEMAOutFilN1,'(f12.4,",",i5,",",4(f12.4,","),a10,",")')JDAY,SegNumI,BubbleRadiusSed(SegNumI)*1000,CgSed(SegNumI),C0Sed(SegNumI),CtSed(SegNumI),CrackStatus(SegNumI)
        ENDDO
        
        !End Do !nGas
        DO SegNumI=1,IMX
            Write(CEMAOutFilN2,'(f12.4,",",i5,",",<IMX>(f12.4,","))')JDAY,SegNumI,(TConc(nGas,KB(SegNumI),SegNumI), nGas = 1, 4)
        ENDDO
        
        
        nGas=2 ! CH4 write
        !Write(CEMASedFlxFilN8,'(f8.2,<IMX>f10.2)')JDAY,(TConc(nGas,SegNumI), SegNumI = 1, IMX)
        Write(CEMASedFlxFilN8,'("DissolvedCH4(mg/l),",f10.3,",",<IMX>(f10.2,","))')JDAY,(DissolvedGasSediments(nGas,KB(SegNumI),SegNumI), SegNumI = 1, IMX)
        nGas=1 ! H2S write
        !Write(CEMASedFlxFilN9,'(f8.2,<IMX>f10.2)')JDAY,(TConc(nGas,SegNumI), SegNumI = 1, IMX)
        Write(CEMASedFlxFilN9,'("DissolvedH2S(mg/l),",f10.3,",",<IMX>(f10.2,","))')JDAY,(DissolvedGasSediments(nGas,KB(SegNumI),SegNumI), SegNumI = 1, IMX)

        !Do nGas = 1, 4
        !    If(nGas == 1)Write(CEMAOutFilN6,'("H2S Concentration (gm/m�")')
        !    If(nGas == 2)Write(CEMAOutFilN6,'("CH4  Concentration (gm/m�")')
        !    If(nGas == 3)Write(CEMAOutFilN6,'("NH3  Concentration (gm/m�")')
        !    If(nGas == 4)Write(CEMAOutFilN6,'("CO2  Concentration (gm/m�")')
        !    Write(CEMAOutFilN6,'(<IMX>f10.3)')(DissolvedGasSediments(nGas,SegNumI), SegNumI = 1, IMX)
        !End Do !nGas
        DO SegNumI=1,IMX
            Write(CEMAOutFilN6,'(f12.4,",",i5,","<IMX>(f12.4,","))')JDAY,SegNumI,(DissolvedGasSediments(nGas,KB(SegNumI),SegNumI), nGas = 1, 4)
        ENDDO
        
        !Write(CEMAOutFilN3,'("Gas Release to Atmosphere at JDAY = ",f8.2)')JDAY
        !Write(CEMAOutFilN3,'(<IMX>i10)')(SegNumI, SegNumI = 1, IMX)
        !Do nGas = 1, 4
        !    If(nGas == 1)Write(CEMAOutFilN3,'("H2S Release(gm/s)")')
        !    If(nGas == 2)Write(CEMAOutFilN3,'("CH4 Release (gm/s)")')
        !    If(nGas == 3)Write(CEMAOutFilN3,'("NH3 Release (gm/s)")')
        !    If(nGas == 4)Write(CEMAOutFilN3,'("CO2 Release (gm/s)")')
        !    Write(CEMAOutFilN3,'(<IMX>f10.3)')(BRRateAGasNet(SegNumI, nGas), SegNumI = 1, IMX)
        !End Do !nGas
        DO SegNumI=1,IMX
            Write(CEMAOutFilN3,'(f12.4,",",i5,","<IMX>(f12.4,","))')JDAY,SegNumI,(BRRateAGasNet(SegNumI, nGas), nGas = 1, 4)
        ENDDO
        
        DO JW=1,NWB
        Write(CEMAOutFilBub,'(F12.4,",",I5,",",6(E12.5,","))')JDAY,JW,(BubbleRelWB(JW,J),J=1,4), GasReleaseCH4/1000.    !,GasReleaseCO2/1000.   ! SW 7/1/2017 BubbleRelWB kg  GasReleaseCH4 kg C
        END DO
        
        !Write(CEMAOutFilN4,'(<5>f10.3)')JDAY,(BRRateAGasNet(5, nGas), nGas = 1, 4)
        !Write(227,'(f10.5, ",",f10.5)')jday, c1(kb(6),6,7)/1000.d0
        
    End Subroutine
        
End Module CEMAOutputRoutines
