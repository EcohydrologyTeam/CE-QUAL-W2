  Subroutine CEMAFFTLayerCode

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
    Use CEMAVars
    
    IMPLICIT NONE
    
    Integer iTemp
    
    NMFT=0
    DO JG=1,NSS  
      IF(SSCS(JG)==-1.0)THEN
        NMFT=NSSS+JG-1       
        EXIT
      ENDIF
    ENDDO
    
    If(FirstTimeInFFTCode)Then
        FirstTimeInFFTCode = .FALSE.
        FFTActive = .TRUE.
        Do JW=1, NWB
            DO JB=BS(JW),BE(JW)
                DO I=CUS(JB),DS(JB)
                    !JAC = NSSS
                    JAC = NMFT      ! cb 2/18/13
                    K=KB(I)
                    C1(K,I,JAC)  = InitFFTLayerConc
                    C1S(K,I,JAC) = InitFFTLayerConc
                    C2(K,I,JAC)  = InitFFTLayerConc
                END DO
            END DO
        End Do
        Return
    End If   
    
    
    Do iTemp = 1, NumFFTActivePrds
    
        If(.NOT. FFTActive)Then
            If(JDAY > FFTActPrdSt(iTemp) .and. JDAY < FFTActPrdEn(iTemp))Then
            
                FFTActive = .TRUE.
                Do JW=1, NWB
                    DO JB=BS(JW),BE(JW)
                        DO I=CUS(JB),DS(JB)
                            !JAC = NSSS
                            JAC = NMFT  ! cb 2/18/13
                            K=KB(I)
                            C1(K,I,JAC)  = FFTLayConc(I)
                            C1S(K,I,JAC) = FFTLayConc(I)
                            C2(K,I,JAC)  = FFTLayConc(I)
                        END DO
                    END DO
                End Do  
                FFTActPrd = iTemp
                Exit  
            End If
        End If
        
    End Do
    
    If(FFTActive)Then
        iTemp = FFTActPrd
        If(JDAY > FFTActPrdEn(iTemp))Then
            FFTActive = .FALSE.
            Do JW=1, NWB
                DO JB=BS(JW),BE(JW)
                    DO I=CUS(JB),DS(JB)
                        !JAC = NSSS
                        JAC = NMFT  ! cb 2/18/13
                        K=KB(I)
                        FFTLayConc(I) = C1(K,I,JAC)
                        C1(K,I,JAC) = 0.d00
                        C1S(K,I,JAC) = 0.d00
                        C2(K,I,JAC) = 0.d00
                    END DO
                END DO
            End Do
        End If
    End If
    
  Return

  Entry MoveFFTLayerConsolid
    
    Do JW=1, NWB
        KT = KTWB(JW)
        Do JB=BS(JW),BE(JW)
            IU = CUS(JB)
            ID = DS(JB)
            Do SegNumI = IU, ID
                If(CEMALayerAdded(SegNumI))Then
                    !JAC = NSSS
                    jac = NMFT  ! cb 2/18/13
                    Do K = KB(SegNumI)-1, KT, -1
                        C1(K+1,SegNumI,JAC)  = C1(K,SegNumI,JAC)
                        C1S(K+1,SegNumI,JAC) = C1S(K,SegNumI,JAC)
                        C2(K+1,SegNumI,JAC)  = C2(K,SegNumI,JAC)
                    End Do !K
                    C1(KT,SegNumI,JAC)  = 0.d00
                    C1S(KT,SegNumI,JAC) = 0.d00
                    C2(KT,SegNumI,JAC)  = 0.d00
                End If
            End Do !SegNumI
        End Do !JB
    End Do !JW

    Return
  End Subroutine
