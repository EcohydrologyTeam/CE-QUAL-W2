Module MetFileRegion    ! SW 12/13/2023

    Integer :: NMetFileRegions,IJUNK,NMet,FNMetRegion
    Integer, Allocatable, Dimension(:) :: MetRegWB,MetRegStart,MetRegEnd,WB_MetRegions,I_MetRegions
    Character*160, Allocatable,  Dimension(:) :: FNMetFileReg
    Character*2 :: MetRegOn
    Logical :: Met_Regions
    
    Contains
    
    Subroutine ReadMetRegions
    Use MAIN, ONLY: CON    
    Implicit None   
    Integer :: J

        
      Met_Regions=.FALSE.
      INQUIRE(FILE='w2_MetRegions.csv',EXIST=Met_Regions)    
     if(Met_regions)then
              FNMetRegion=CON
              OPEN(FNMetRegion,FILE='W2_MetRegions.csv',status='old')
              READ(FNMetRegion,*)
              READ(FNMetRegion,'(A2)')MetRegOn
          if(MetRegOn =='ON') then
                  READ(FNMetRegion,*)
                  READ(FNMetRegion,*)NMetFileRegions
                  ALLOCATE(MetRegWB(NMetFileRegions),MetRegStart(NMetFileRegions),MetRegEnd(NMetFileRegions),FNMetFileReg(NMetFileRegions))
                  READ(FNMetRegion,*)
              DO J=1,NMetFileRegions
                  READ(FNMetRegion,*)IJUNK,MetRegStart(J),MetRegEnd(J),FNMetFileReg(J)
              ENDDO
          else
                  Met_Regions=.FALSE.            
          endif
          CLOSE(CON)
     endif

    RETURN
    
    End Subroutine ReadMetRegions
    
    Subroutine MetRegionsWB
    Use GLOBAL
    Implicit None    
    !integer :: NMet

    Allocate(WB_MetRegions(NWB),I_MetRegions(IMX))
    MetRegWB=0
    WB_MetRegions=0
    I_MetRegions=0
    DO NMet=1,NMetFileRegions
    
    DO JW=1,NWB
            IF(MetRegStart(NMet) >= US(BS(JW))-1 .AND. MetRegEnd(NMet) <=  DS(BE(JW))+1)THEN
                 MetRegWB(NMet)=JW
                 WB_MetRegions(JW)=NMet
                 EXIT
            ENDIF
    ENDDO
    if(MetRegWB(NMet)==0)then
        DO JW=1,NWB
        if(MetRegStart(NMet) >= US(BS(JW))-1.and.MetRegStart(NMet) <= DS(BS(JW))-1 )THEN
                 MetRegWB(NMet)=JW
                 EXIT
        ENDIF
        ENDDO
    endif
    ENDDO
    
    DO JW=1,NWB
        if(WB_MetRegions(JW)==0)then
        DO NMet=1,NMetFileRegions
        IF(MetRegStart(NMet) <= US(BS(JW))-1 .AND. MetRegEnd(NMet) >=  DS(BE(JW))+1)then
            WB_MetRegions(JW)=NMet
            exit
        endif
        ENDDO
        endif
    ENDDO
    
    DO JW=1,NWB
    DO I=US(BS(JW))-1,DS(BE(JW))+1
        DO NMet=1,NMetFileRegions
        IF(MetRegStart(NMet) <= I .and. MetRegEnd(NMet) >= I)then
            I_MetRegions(I)=NMet
            exit
        ENDIF   
    ENDDO
    ENDDO
    ENDDO
    
    RETURN
    
    End Subroutine MetRegionsWB
    
    Subroutine EndMetRegion
    DEALLOCATE(MetRegWB,MetRegStart,MetRegEnd,FNMetFileReg,WB_MetRegions,I_MetRegions)
    RETURN
    
    End Subroutine EndMetRegion
    
End Module MetFileRegion
    