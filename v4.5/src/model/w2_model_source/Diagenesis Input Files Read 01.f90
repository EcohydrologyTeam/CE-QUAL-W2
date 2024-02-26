Subroutine InitializeBedConsolidationFiles(TempFilNum, TempFilName)
    Use MAIN 
    Use GLOBAL
    Use SCREENC
    Use CEMAVars
    
    Implicit None
    
    Logical SkipLoop
    Integer(4) TempFilNum
    Character(256) MessageTemp, TempFilName
    
    !Open File	
	Open(TempFilNum, File = TempFilName(1:len_trim(TempFilName)-1))
	
	!Read Header
	SkipLoop = .FALSE.
	Do While(.NOT. SkipLoop)
		Read(TempFilNum,'(a)')MessageTemp
		If(index(MessageTemp, "$") == 0)SkipLoop = .TRUE.
	End Do
			
End Subroutine

Subroutine ReadBedConsolidationFiles(TempFilNum)
    Use MAIN 
    Use GLOBAL
    Use SCREENC
    Use CEMAVars
    
    Logical    :: SkipLoop
    Integer(4) :: TempFilNum
    Real(8) :: TimeJD1
    Real(8) :: TimeJD2
    Real(8) :: FactorInterp
	Real(8) :: ConsolidRateTemp11(NumConsolidRegns),ConsolidRateTemp1(NumConsolidRegns),ConsolidRateTemp2(NumConsolidRegns)																													
    
    !Read Data
	SkipLoop = .FALSE.
	Do While(.NOT. SkipLoop .or. EOF(TempFilNum))	
	  Read(TempFilNum,'(F8.0,<NumConsolidRegns>F8.0)')TimeJD1, (ConsolidRateTemp1(i),i=1,NumConsolidRegns)
	  Read(TempFilNum,'(F8.0,<NumConsolidRegns>F8.0)')TimeJD2, (ConsolidRateTemp2(i),i=1,NumConsolidRegns)

	    If(JDay >= TimeJD1 .and. JDay <= TimeJD2)then
		    SkipLoop = .TRUE.
	    End If
	    BackSpace(TempFilNum)
	End Do
	
	BackSpace(TempFilNum)
	
	FactorInterp = (JDay - TimeJD1)/(TimeJD2 - TimeJD1)

  DO i=1,NumConsolidRegns
    ConsolidRateTemp11(i) = ConsolidRateTemp1(i)*(1-FactorInterp) + ConsolidRateTemp2(i)*FactorInterp
  END DO
	
End Subroutine