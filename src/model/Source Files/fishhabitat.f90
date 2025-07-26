SUBROUTINE FISHHABITAT(IOPENFISH)
USE GLOBAL;USE MAIN;USE SCREENC; USE KINETIC, ONLY:O2, CHLA, NO3,NH4,PO4, TP, GAMMA, SED, SATO; USE TVDC, ONLY: CONSTITUENTS; USE NAMESC, ONLY: CNAME2,CDNAME2; USE LOGICC
IMPLICIT NONE
INTEGER :: IFISH, N, NSEG,IOPENFISH,KKMAX,KSEG,JJW,JWFILE,JBFILE, X    !,JWFILE1,JBFILE1
REAL    :: VOLTOT,O2CORR,DOSAT,VOLDOALL
CHARACTER*80, ALLOCATABLE, DIMENSION(:) :: FISHNAME
CHARACTER*80 :: CONHAB,CONAVG,CONSURF,CONSOD
CHARACTER*2 :: DOVOL
REAL, ALLOCATABLE, DIMENSION(:)   :: FISHTEMPL,FISHTEMPH,FISHDO,HABVOL,PHABVOL,CDO,CPO4,CNO3,CNH4,CCHLA,CTOTP,CDOS,CPO4S,CNO3S,CNH4S,CCHLAS,CTOTPS,CGAMMA,SSEDD,VOLTOTBR,VOLTOTWB
REAL, ALLOCATABLE, DIMENSION(:,:) :: HABVOLBR,HABVOLWB,PHABVOLBR,PHABVOLWB
REAL, ALLOCATABLE, DIMENSION(:)   :: VOLDOBR,VOLDOWB
INTEGER, ALLOCATABLE, DIMENSION (:) :: ISEGVOL
SAVE

        !JWFILE1=9549
        !JBFILE1=9749
        
IF(IOPENFISH .NE. 3)THEN
! read input file

if(nit == 1.or.iopenfish==0)then
open(FISHHABFN,file='w2_habitat.npt',status='old')
! skip 1st 2 lines
read(FISHHABFN,*)
read(FISHHABFN,*)
read(FISHHABFN,*)ifish,conhab,DOVOL
conhab=TRIM(ADJUSTL(conhab))        !ZZ write fish_habitat output files into a folder specified here.
X = SCAN(conhab,'\',BACK=.TRUE.)      
allocate (fishname(ifish),fishtempl(ifish),fishtemph(ifish),fishdo(ifish),habvol(ifish),phabvol(ifish),habvolbr(nbr,ifish),habvolwb(nwb,ifish),phabvolbr(nbr,ifish),phabvolwb(nwb,ifish),voltotbr(nbr),voltotwb(nwb))
read(FISHHABFN,*)
do i=1,ifish
  read(FISHHABFN,*)fishname(i),fishtempl(i),fishtemph(i),fishdo(i)
enddo
read(FISHHABFN,*)
read(FISHHABFN,*)nseg,conavg     ! volume weighted averages of critical WQ 
conavg=TRIM(ADJUSTL(conavg))
read(FISHHABFN,*)
allocate(isegvol(nseg),cdo(nseg),cpo4(nseg),cno3(nseg),cnh4(nseg),cchla(nseg),ctotp(nseg),cdos(nseg),cpo4s(nseg),cno3s(nseg),cnh4s(nseg),cchlas(nseg),ctotps(nseg),cgamma(nseg))
allocate(ssedd(imx))
ALLOCATE(VOLDOBR(NBR),VOLDOWB(NWB))
read(FISHHABFN,*)(isegvol(i),i=1,nseg)
read(FISHHABFN,*)
read(FISHHABFN,*)kseg,consurf     ! # of layers for surface 
consurf=TRIM(ADJUSTL(consurf))
read(FISHHABFN,*)
read(FISHHABFN,*)consod
consod=TRIM(ADJUSTL(consod))
close(FISHHABFN)


if(restart_in)then
        open(FISHHABFN,file=conhab,POSITION='APPEND')
        JDAY1=0.0
        REWIND (FISHHABFN)
        READ   (FISHHABFN,'(///)')
        DO I=1,IFISH
        READ(FISHHABFN,*)
        ENDDO
        READ(FISHHABFN,'(//)',END=106)
        
        DO WHILE (JDAY1 < JDAY)
          READ (FISHHABFN,'(F10.0)',END=106) JDAY1
        END DO
        BACKSPACE (FISHHABFN)
        106     JDAY1=0.0
        jbfile=jbfile1;jwfile=jwfile1      
        do jw=1,nwb
        jwfile=jwfile+1
         WRITE (SEGNUM,'(I0)') JW
         SEGNUM = ADJUSTL(SEGNUM)
         L      = LEN_TRIM(SEGNUM)
         IF (X > 0) THEN
           OPEN(JWFILE,FILE=conhab(1:X)//'fish_habitat_wb'//SEGNUM(1:L)//'.csv',POSITION='APPEND')
         ELSE
           OPEN(JWFILE,FILE='fish_habitat_wb'//SEGNUM(1:L)//'.csv',POSITION='APPEND') 
         END IF
         JDAY1=0.0
         REWIND (JWFILE)
        READ   (JWFILE,'(///)')
        DO I=1,IFISH
        READ(JWFILE,*)
        ENDDO
        READ(JWFILE,'(//)',END=110)
        DO WHILE (JDAY1 < JDAY)
          READ (JWFILE,'(F10.0)',END=110) JDAY1
        END DO
        BACKSPACE (JWFILE)
        110     JDAY1=0.0
                
            do jb=bs(jw),be(jw)
            jbfile=jbfile+1
            WRITE (SEGNUM,'(I0)') JB
            SEGNUM = ADJUSTL(SEGNUM)
            L      = LEN_TRIM(SEGNUM)
            IF (X > 0) THEN
              OPEN(JBFILE,FILE=conhab(1:X)//'fish_habitat_br'//SEGNUM(1:L)//'.csv',POSITION='APPEND') 
            ELSE
              OPEN(JBFILE,FILE='fish_habitat_br'//SEGNUM(1:L)//'.csv',POSITION='APPEND')  
            END IF
            JDAY1=0.0
            REWIND (JBFILE)
            READ   (JBFILE,'(///)')
            DO I=1,IFISH
            READ(JBFILE,*)
            ENDDO
            READ(JBFILE,'(//)',END=111)
            DO WHILE (JDAY1 < JDAY)
            READ (JBFILE,'(F10.0)',END=111) JDAY1
            END DO
            BACKSPACE (JBFILE)
            111     JDAY1=0.0 
                        
            enddo
        enddo                 
                  
        if(oxygen_demand)then
            open(FISHHABFN+1,file=conavg,POSITION='APPEND')
            REWIND (FISHHABFN+1)
            READ   (FISHHABFN+1,'(//)')
            DO WHILE (JDAY1 < JDAY)
            READ (FISHHABFN+1,'(F10.0)',END=107) JDAY1
            END DO
            BACKSPACE (FISHHABFN+1)
            107     JDAY1=0.0
            open(FISHHABFN+2,file=consurf,POSITION='APPEND')
            REWIND (FISHHABFN+2)
            READ   (FISHHABFN+2,'(//)')
            DO WHILE (JDAY1 < JDAY)
            READ (FISHHABFN+2,'(F10.0)',END=108) JDAY1
            END DO
            BACKSPACE (FISHHABFN+2)
            108     JDAY1=0.0
            
                do jjw=1,nwb
                IF (SEDIMENT_CALC(JJW))then
                open(FISHHABFN+3,file=consod,POSITION='APPEND')
                REWIND (FISHHABFN+3)
                READ   (FISHHABFN+3,'(/)')
                DO WHILE (JDAY1 < JDAY)
                READ (FISHHABFN+3,'(F10.0)',END=109) JDAY1
                END DO
                BACKSPACE (FISHHABFN+3)
                109     JDAY1=0.0
                exit
                ENDIF
                enddo
         endif            
else

        open(FISHHABFN,file=conhab,status='unknown')
        write(FISHHABFN,*)'Fish habitat analysis: CE-QUAL-W2 model results'
        write(FISHHABFN,*)
        write(FISHHABFN,*)'Species, Temperature minimum, Temperature maximum, Dissolved oxygen minimum'
        do i=1,ifish
        write(FISHHABFN,"(a,',',t25,f8.2,',',f8.2,',',f8.2)")trim(fishname(i)),fishtempl(i),fishtemph(i),fishdo(i)
        enddo
        write(FISHHABFN,*)
        IF(DOVOL=='ON')THEN
        write(FISHHABFN,101)(trim(fishname(i)),trim(fishname(i)),i=1,ifish)
101     format('JDAY,',<ifish>('%VOL-',A,',','HAB-VOL(m3)-',A,','),'DO_VOL(m3)<=1mg/L',',')   
        ELSE
        write(FISHHABFN,100)(trim(fishname(i)),trim(fishname(i)),i=1,ifish)
100     format('JDAY,',<ifish>('%VOL-',A,',','HAB-VOL(m3)-',A,','))
        ENDIF
            
        jbfile=jbfile1;jwfile=jwfile1       
        do jw=1,nwb
        jwfile=jwfile+1
         WRITE (SEGNUM,'(I0)') JW
         SEGNUM = ADJUSTL(SEGNUM)
         L      = LEN_TRIM(SEGNUM)
         IF (X > 0) THEN
           OPEN(JWFILE,FILE=conhab(1:X)//'fish_habitat_wb'//SEGNUM(1:L)//'.csv',STATUS='UNKNOWN')
         ELSE
           OPEN(JWFILE,FILE='fish_habitat_wb'//SEGNUM(1:L)//'.csv',STATUS='UNKNOWN') 
         END IF
                write(JWFILE,*)'Fish habitat analysis: CE-QUAL-W2 model results'
                write(JWFILE,*)'FOR WATERBODY:',JW
                write(JWFILE,*)'Species, Temperature minimum, Temperature maximum, Dissolved oxygen minimum'
                do i=1,ifish
                write(JWFILE,"(a,',',t25,f8.2,',',f8.2,',',f8.2)")trim(fishname(i)),fishtempl(i),fishtemph(i),fishdo(i)
                enddo
                IF(DOVOL=='ON')THEN
                write(JWFILE,*)
                write(JWFILE,101)(trim(fishname(i)),trim(fishname(i)),i=1,ifish)   
                    ELSE
                write(JWFILE,*)
                write(JWFILE,100)(trim(fishname(i)),trim(fishname(i)),i=1,ifish)
                ENDIF
                
            do jb=bs(jw),be(jw)
            jbfile=jbfile+1
            WRITE (SEGNUM,'(I0)') JB
            SEGNUM = ADJUSTL(SEGNUM)
            L      = LEN_TRIM(SEGNUM)
            IF (X > 0) THEN
              OPEN(JBFILE,FILE=conhab(1:X)//'fish_habitat_br'//SEGNUM(1:L)//'.csv',STATUS='UNKNOWN')
            ELSE
              OPEN(JBFILE,FILE='fish_habitat_br'//SEGNUM(1:L)//'.csv',STATUS='UNKNOWN')  
            END IF
                write(JBFILE,*)'Fish habitat analysis: CE-QUAL-W2 model results'
                write(JBFILE,*)'FOR BRANCH:',JB
                write(JBFILE,*)'Species, Temperature minimum, Temperature maximum, Dissolved oxygen minimum'
                do i=1,ifish
                write(JBFILE,"(a,',',t25,f8.2,',',f8.2,',',f8.2)")trim(fishname(i)),fishtempl(i),fishtemph(i),fishdo(i)
                enddo
                IF(DOVOL=='ON')THEN
                write(JBFILE,*)
                write(JBFILE,101)(trim(fishname(i)),trim(fishname(i)),i=1,ifish)  
                    ELSE
                write(JBFILE,*)
                write(JBFILE,100)(trim(fishname(i)),trim(fishname(i)),i=1,ifish)
                ENDIF
            enddo
        enddo

  if(oxygen_demand)then
    !ZZ include additional rows in three output files
    open(FISHHABFN+1,file=conavg,status='unknown')
    write(FISHHABFN+1,*)'Fish habitat analysis: CE-QUAL-W2 model results'
    write(FISHHABFN+1,'(" FOR SEGMENTS:",*(1x,i4))') (isegvol(i),i=1,nseg)
    write(FISHHABFN+1,*)'Volume-weighted averages are computed at specified segments:'
    do i=1,ifish
    write(FISHHABFN+1,*)
    enddo
    !
    write(FISHHABFN+1,'(a,80(1x,i4))')'Volume weighted WQ parameters at segments:',(isegvol(i),i=1,nseg)
    write(FISHHABFN+1,103)(trim(cname2(NPO4)),isegvol(i),trim(cname2(NNH4)),isegvol(i),trim(cname2(NNO3)),isegvol(i),trim(cname2(NDO)),isegvol(i),trim(cdname2(12)),isegvol(i),trim(cdname2(14)),isegvol(i),i=1,nseg)   ! Chlor a and TP 
    103 format('JDAY,',<nseg>(6((A,'-',i3,','))))
    102 format('JDAY,',<nseg>(7((A,'-',i3,','))))  
    open(FISHHABFN+2,file=consurf,status='unknown')
    write(FISHHABFN+2,*)'Fish habitat analysis: CE-QUAL-W2 model results'
    write(FISHHABFN+2,'(" FOR SEGMENTS:",*(1x,i4))') (isegvol(i),i=1,nseg)
    write(FISHHABFN+2,*)'Volume weighted averages at specified number of layers of the surface:'
    do i=1,ifish
    write(FISHHABFN+2,*)
    enddo
    !
    write(FISHHABFN+2,'(a,i4,a,80(1x,i4))')'Surface (upper',kseg,' model layers) Volume weighted WQ parameters at segments:',(isegvol(i),i=1,nseg)
    write(FISHHABFN+2,102)(trim(cname2(NPO4)),isegvol(i),trim(cname2(NNH4)),isegvol(i),trim(cname2(NNO3)),isegvol(i),trim(cname2(NDO)),isegvol(i),trim(cdname2(12)),isegvol(i),trim(cdname2(14)),isegvol(i),'Gamma(m-1)',isegvol(i),i=1,nseg)

    do jjw=1,nwb
    IF (SEDIMENT_CALC(JJW))then
    open(FISHHABFN+3,file=consod,status='unknown')
    write(FISHHABFN+3,*)'Fish habitat analysis: CE-QUAL-W2 model results'
    write(FISHHABFN+3,*)'FOR ALL SEGMENTS'
    write(FISHHABFN+3,*)'Overall organic matter accumulation (g) at the bottom of each layer:'
    do i=1,ifish
    write(FISHHABFN+3,*)
    enddo
    write(FISHHABFN+3,*)
    !
    write(FISHHABFN+3,"('JDAY,',1000(i3,','))")(((i,i=us(jb),ds(jb)),jb=bs(jw),be(jw)),jw=1,nwb)
    exit
    ENDIF
    enddo
  endif

endif
endif

! compute total volume and habitat volume

habvol=0.0
voltot=0.0
HAB=100.0
voltotbr=0.0
voltotwb=0.0
habvolbr=0.0
habvolwb=0.0
VOLDOBR=0.0
VOLDOWB=0.0
VOLDOALL=0.0
do jw=1,nwb
    do jb=bs(jw),be(jw)
        do i=cus(jb),ds(jb)
            do k=ktwb(jw),kb(i)
                voltot=voltot+vol(k,i)
                voltotbr(jb)=voltotbr(jb)+vol(k,i)
                voltotwb(jw)=voltotwb(jw)+vol(k,i)
                IF(DOVOL=='ON')THEN
                    IF(O2(K,I) <= 1.0)THEN
                        VOLDOALL=VOLDOALL+VOL(K,I)
                        VOLDOBR(JB)=VOLDOBR(JB)+VOL(K,I)
                        VOLDOWB(JW)=VOLDOWB(JW)+VOL(K,I)
                    ENDIF   
                ENDIF
                    do ii=IFISH,1,-1
                    if(oxygen_demand)then
                        if(t2(k,i)<=fishtemph(ii).and.t2(k,i)>fishtempl(ii).and.o2(k,i)>=fishdo(ii))then
                            habvol(ii)=habvol(ii)+vol(k,i)
                            habvolbr(jb,ii)=habvolbr(jb,ii)+vol(k,i)
                            habvolwb(jw,ii)=habvolwb(jw,ii)+vol(k,i)
                            hab(k,i)=ii
                        endif
                    else
                        if(t2(k,i)<=fishtemph(ii).and.t2(k,i)>fishtempl(ii))then
                            habvol(ii)=habvol(ii)+vol(k,i)
                            habvolbr(jb,ii)=habvolbr(jb,ii)+vol(k,i)
                            habvolwb(jw,ii)=habvolwb(jw,ii)+vol(k,i)
                            hab(k,i)=ii
                        endif
                    endif
                    enddo
            enddo
        enddo
    enddo
enddo    

do ii=1,ifish
phabvol(ii)=habvol(ii)/voltot
    do jw=1,nwb
        phabvolwb(jw,ii)=habvolwb(jw,ii)/voltotwb(jw)   
        do jb=bs(jw),be(jw)
        phabvolbr(jb,ii)=habvolbr(jb,ii)/voltotbr(jb)   
        end do
    enddo
enddo

! write out results
IF(DOVOL=='ON')THEN
write(FISHHABFN,211)jday,(100.*phabvol(i),habvol(i),i=1,ifish),VOLDOALL    
ELSE
write(FISHHABFN,210)jday,(100.*phabvol(i),habvol(i),i=1,ifish)
ENDIF

210 format(f10.3,',',<ifish>(f8.2,',',e12.4,','))
211 format(f10.3,',',<ifish>(f8.2,',',e12.4,','),E12.4)
jbfile=jbfile1;jwfile=jwfile1 
do jw=1,nwb
    jwfile=jwfile+1
    IF(DOVOL=='ON')THEN
    write(jwfile,211)jday,(100.*phabvolwb(jw,i),habvolwb(jw,i),i=1,ifish),VOLDOWB(JW)   
        ELSE
    write(jwfile,210)jday,(100.*phabvolwb(jw,i),habvolwb(jw,i),i=1,ifish)
    ENDIF
    do jb=bs(jw),be(jw)
        jbfile=jbfile+1
        IF(DOVOL=='ON')THEN
        write(jbfile,211)jday,(100.*phabvolbr(jb,i),habvolbr(jb,i),i=1,ifish),VOLDOBR(JB)    
            ELSE
        write(jbfile,210)jday,(100.*phabvolbr(jb,i),habvolbr(jb,i),i=1,ifish)
        ENDIF
    enddo
enddo

if(oxygen_demand)then

    cno3=0.0
    cdo=0.0
    cpo4=0.0
    cchla=0.0
    cnh4=0.0
    cgamma=0.0
    ctotp=0.0

    do n=1,nseg

    i=isegvol(n)

    ! Find waterbody associated with this segment
    do jjw=1,nwb
        if(i >= us(bs(jjw)) .and. i <= ds(be(jjw)))exit
    enddo

    voltot=0.0
    kkmax=min(kseg,kb(i)-ktwb(jjw))      ! kseg is the # of layers
    if(kkmax < 0)cycle
        do k=ktwb(jjw),kb(i)   
        voltot=voltot+vol(k,i)
        cpo4(n)=cpo4(n)+po4(k,i)*vol(k,i)
        if(k <= ktwb(jjw)+kkmax)cgamma(n)=cgamma(n)+gamma(k,i)*vol(k,i)
        ! NOTE*********** No credit for superstauration - if DO > saturation, then set DO=100% saturation
        DOSAT=SATO(t2(k,i),0.d0,palt(i),SALT_WATER(jjw))   
            if(o2(k,i) > DOSAT )then 
            o2corr=DOSAT
            else
            o2corr=o2(k,i)
            endif
    
        cdo(n)=cdo(n)+o2corr*vol(k,i)
        cno3(n)=cno3(n)+no3(k,i)*vol(k,i)
        cchla(n)=cchla(n)+chla(k,i)*vol(k,i)
        ctotp(n)=ctotp(n)+tp(k,i)*vol(k,i)
        cnh4(n)=cnh4(n)+nh4(k,i)*vol(k,i)
    
        if(k == ktwb(jjw)+kkmax)then
        cdos(n)=cdo(n)/voltot
        cpo4s(n)=cpo4(n)/voltot
        cno3s(n)=cno3(n)/voltot
        cnh4s(n)=cnh4(n)/voltot
        cchlas(n)=cchla(n)/voltot
        ctotps(n)=ctotp(n)/voltot
        cgamma(n)=cgamma(n)/voltot   
        endif
    
    
        enddo
    cpo4(n)=cpo4(n)/voltot
    cdo(n)=cdo(n)/voltot
    cno3(n)=cno3(n)/voltot
    cnh4(n)=cnh4(n)/voltot
    cchla(n)=cchla(n)/voltot
    ctotp(n)=ctotp(n)/voltot
    enddo

    write(FISHHABFN+2,212)jday,(cpo4s(n),cnh4s(n),cno3s(n),cdos(n),ctotps(n),cchlas(n),cgamma(n),n=1,nseg)
    write(FISHHABFN+1,212)jday,(cpo4(n),cnh4(n),cno3(n),cdo(n),ctotp(n),cchla(n),n=1,nseg)
    212 format(f10.3,',',<nseg>(7(f10.4,',')))



! write out sed for each segment
ssedd=0.0
do jw=1,nwb
   IF (SEDIMENT_CALC(JW))then
    do jb=bs(jw),be(jw)
        do i=us(jb),ds(jb)
           if(ktwb(jw)<=kb(i))then
            do k=ktwb(jw),kb(i)
                ssedd(i)=ssedd(i)+sed(k,i)*vol(k,i)          
            end do 
           else
            ssedd(i)=-99.
           endif
        enddo
    enddo
  ENDIF
enddo

do jjw=1,nwb
IF (SEDIMENT_CALC(JJW))then
write(FISHHABFN+3,"(f10.3,',',1000(e13.4,','))")jday,(((ssedd(i),i=us(jb),ds(jb)),jb=bs(jw),be(jw)),jw=1,nwb)    ! OUTPUT IS IN GRAMS
exit
ENDIF
enddo

endif


else

deallocate(isegvol,cdo,cpo4,cno3,cnh4,cchla,ctotp,cdos,cpo4s,cno3s,cnh4s,cchlas,ctotps,cgamma,ssedd,fishname,fishtempl,fishtemph,fishdo,habvol,phabvol,habvolbr,habvolwb,phabvolbr,phabvolwb,voltotbr,voltotwb)
deallocate(VOLDOBR,VOLDOWB)
close(FISHHABFN)
close(FISHHABFN+1)
CLOSE(FISHHABFN+2)
CLOSE(FISHHABFN+3)
jbfile=jbfile1;jwfile=jwfile1 
    do jw=1,nwb
     jwfile=jwfile+1
    close(jwfile)
        do jb=bs(jw),be(jw)
        jbfile=jbfile+1
        close(jbfile)
        enddo
    enddo

endif


return
end subroutine fishhabitat