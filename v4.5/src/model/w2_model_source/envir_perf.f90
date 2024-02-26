!*******************************************************************
!**           S U B R O U T I N E   E N V I R P
!*******************************************************************
subroutine envirp

USE GLOBAL; USE MAIN;use NAMESC; use screenc, only:nit,jday;use tvdc, only: constituents; use rstart, only:eltm; use GEOMC, only:depthb
USE ENVIRPMOD; USE GDAYC
IMPLICIT NONE
CHARACTER(1) :: I_INT
CHARACTER(8) :: CHAR8
character(30):: CHAR30
INTEGER :: N
LOGICAL :: CSVFORMAT

save

!  initializing variables at first call
 if(NIT==1.or.iopenfish==0)then
      dltt=dlt
      allocate(cc_e(NCT),c_int(NCT),c_top(NCT),cd_e(NDC),cd_int(NDC),cd_top(NDC),c_avg(NCT),cd_avg(NDC),cn_e(NCT),cdn_e(NDC))
      cc_e='   '
      c_int=0.0
      c_top=0.0
      cd_e='   '
      cd_int=0.0
      cd_top=0.0
      c_avg=0.0
      cd_avg=0.0
      cn_e=0.0
      cdn_e=0.0
      NAC_E=0
      NACD_E=0
      CONE=NUNIT; NUNIT=NUNIT+1
      open(CONE,file='w2_envirprf.npt',status='old')
      
     CSVFORMAT=.FALSE.
     READ(CONE,'(//A)')CHAR30
     DO J=1,30
         IF(CHAR30(J:J)==',')THEN
             CSVFORMAT=.TRUE.
             EXIT
         ENDIF
     ENDDO
     REWIND(CONE)

      IF(CSVFORMAT)THEN
        READ(CONE,*)
        READ(CONE,*)
        READ (CONE,*) I_SEGINT,numclass,selectivec,sjday1,sjday2,istart(1),iend(1),(istart(I),iend(I),I=2,I_SEGINT)
        SELECTIVEC=ADJUSTR(SELECTIVEC)
        IF(I_SEGINT==0 .OR. SELECTIVEC=='OFF')I_SEGINT=1
        READ(CONE,*)
        READ(CONE,*)
        Read (CONE,*) VEL_VPR, VEL_INT, VEL_TOP,TEMP_VPR,TEMP_INT,TEMP_TOP, depth_vpr,d_int,d_top
        VEL_VPR=ADJUSTR(VEL_VPR);TEMP_VPR=ADJUSTR(TEMP_VPR);DEPTH_VPR=ADJUSTR(DEPTH_VPR)
        READ(CONE,*)
        READ(CONE,*)
        DO JC=1,NCT
        READ (CONE,*) CHAR8, CC_E(JC), C_INT(JC), C_TOP(JC)
        CC_E(JC)=ADJUSTR(CC_E(JC))
        ENDDO
        READ(CONE,*)
        READ(CONE,*)
        DO JD=1,NDC
        READ (CONE,*) CHAR8, CD_E(JD),CD_INT(JD), CD_TOP(JD)
        CD_E(JD)=ADJUSTR(CD_E(JD))        
        ENDDO
      ELSE
      READ (CONE,1200) I_SEGINT,numclass,selectivec,sjday1,sjday2,istart(1),iend(1),(istart(I),iend(I),I=2,I_SEGINT)
      IF(I_SEGINT==0 .OR. SELECTIVEC=='OFF')I_SEGINT=1
      Read (CONE,1201) VEL_VPR, VEL_INT, VEL_TOP,TEMP_VPR,TEMP_INT,TEMP_TOP, depth_vpr,d_int,d_top
      READ (CONE,1050) (CC_E(JC), C_INT(JC), C_TOP(JC), JC=1,NCT)
      READ (CONE,1050) (CD_E(JD),CD_INT(JD), CD_TOP(JD), JD=1,NDC)
      ENDIF
      CLOSE(CONE)

          DO JC=1,NCT
          IF (CC_E(JC).EQ.' ON') THEN
            NAC_E     = NAC_E+1
            CN_E(NAC_E) = JC
          END IF
          End DO
          DO JD=1,NDC
          IF (CD_E(JD).EQ.' ON') THEN
            NACD_E     = NACD_E+1
            CDN_E(NACD_E) = JD
          END IF
          End DO

 1050 FORMAT(//(8X,(5X,A3,F8.0,F8.0)))
 1200 FORMAT(//I1,7X,I8,5x,a3,f8.0,f8.0,9(i8,i8))        ! UP TO A MAX OF 9 INTERVALS
 1201 format(//8x,3(5x,a3,f8.3,f8.3))

    allocate (c_cnt(I_SEGINT,NCT),cd_cnt(I_SEGINT,NDC),c_class(I_SEGINT,NCT,numclass),cd_class(I_SEGINT,NDC,numclass),c_tot(I_SEGINT,NCT),cd_tot(I_SEGINT,NDC),t_class(I_SEGINT,numclass),v_class(I_SEGINT,numclass),c_sum(NCT),cd_sum(NDC))
    allocate (conc_c(NCT,numclass),conc_cd(NDC,numclass))
    allocate(d_class(I_SEGINT,numclass))
    ALLOCATE (D_TOT(I_SEGINT),D_CNT(I_SEGINT),T_TOT(I_SEGINT),T_CNT(I_SEGINT))
    ALLOCATE(V_TOT(I_SEGINT),V_CNT(I_SEGINT),VOLGL(I_SEGINT),SUMVOLT(I_SEGINT))
        cd_cnt=0.0
        c_cnt=0.0
        c_class=0.0
        cd_class=0.0
        c_tot=0.0
        cd_tot=0.0
        
        v_cnt=0.0
        v_class=0.0
        v_tot=0.0

        t_cnt=0.0
        t_class=0.0
        t_tot=0.0
        
        d_cnt=0.0
        d_class=0.0
        d_tot=0.0

      sumvolt=0.0
  else
      dltt=(jday-timlast)*86400.
  end if

  if(iopenfish.eq.3)go to 650     !iopenfish=3 is end of simulation deallocate arrays
  
  if(selectivec == ' ON')then
      if(jdayG < sjday1  .or. jdayG > sjday2)go to 650
  endif
  
DO N=1,I_SEGINT            ! LOOP OVER SEGMENT INTERVALS BASED ON INPUT DATA  
  
   volgL(N)=0.0
! start loop for succeeding calls to subroutine

do JW=1,NWB
  do jb=bs(JW),be(JW)
      do i=cus(jb),ds(jb)
        if(selectivec == ' ON')then
            !if(i < istart(N) .or. i > iend(N))then
            IF(i < istart(N))THEN 
                CYCLE    ! SW 2/16/2017   exit
            ELSEIF(i > iend(N))THEN
                EXIT
            endif
        endif
          
! Depth
           if(depth_vpr.eq.' ON')then
            d_tot(N)=d_tot(N)+depthb(kb(i),i)*dltt
            d_cnt(N)=d_cnt(N)+dltt
            d_crit=d_top
            if(depthb(kb(i),i).ge.d_top)d_class(N,1)=d_class(N,1)+dltt
            do jj=2,numclass
              if(depthb(kb(i),i).lt.d_crit.and.depthb(kb(i),i).ge.d_crit-d_int)then
                d_class(N,jj)=d_class(N,jj)+dltt
                go to 300
              else
                d_crit=d_crit-d_int
              end if
              if(jj.eq.numclass.and.depthb(kb(i),i).lt.d_crit+d_int)then
                d_class(N,jj)=d_class(N,jj)+dltt
              end if
            end do
          end if
300       continue

          
          
          
          
        do k=KTWB(JW),kb(i)
            volgL(N)=volgL(N)+vol(K,I)

! Temperature 

        if(temp_vpr.eq.' ON')then

            t_tot(N)=t_tot(N)+t2(k,i)*VOL(k,i)*dltt
            t_cnt(N)=t_cnt(N)+vol(k,i)*dltt
            t_crit=temp_top
            if(t2(k,i).ge.temp_top)t_class(N,1)=t_class(N,1)+dltt*vol(k,i)
            do jj=2,numclass
              if(t2(k,i).lt.t_crit.and.t2(k,i).ge.t_crit-temp_int)then
                t_class(N,jj)=t_class(N,jj)+dltt*vol(k,i)
                go to 200
              else
                t_crit=t_crit-temp_int
              end if
              if(jj.eq.numclass.and.t2(k,i).lt.t_crit+temp_int)then
                t_class(N,jj)=t_class(N,jj)+dltt*vol(k,i)
              end if
            end do
          end if
200       continue


! Velocity

        if(vel_vpr.eq.' ON')then

            v_tot(N)=v_tot(N)+u(k,i)*vol(k,i)*dltt
            v_cnt(N)=v_cnt(N)+vol(k,i)*dltt
            v_crit=vel_top
            if(u(k,i).ge.vel_top)v_class(N,1)=v_class(N,1)+dltt*vol(k,i)
            do jj=2,numclass
              if(u(k,i).lt.v_crit.and.u(k,i).ge.v_crit-vel_int)then
                v_class(N,jj)=v_class(N,jj)+dltt*vol(k,i)
                go to 210
              else
                v_crit=v_crit-vel_int
              end if
              if(jj.eq.numclass.and.u(k,i).lt.v_crit+vel_int)then
                v_class(N,jj)=v_class(N,jj)+dltt*vol(k,i)
              end if
            end do
          end if
210       continue


! Constituents
  IF(CONSTITUENTS)THEN
            do jc=1,nac_e
              jac=cn_e(jc)

            c_tot(N,jc)=c_tot(N,jc)+c2(k,i,jac)*cmult(jac)*vol(k,i)*dltt
            c_cnt(N,jc)=c_cnt(N,jc)+vol(k,i)*dltt
            c_crit=c_top(jac)
            if(c2(k,i,jac)*cmult(jac).ge.c_top(jac))c_class(N,jc,1)=c_class(N,jc,1)+dltt*vol(k,i)
            do jj=2,numclass
              if(c2(k,i,jac)*cmult(jac).lt.c_crit.and.c2(k,i,jac)*cmult(jac).ge.c_crit-c_int(jac))then
                c_class(N,jc,jj)=c_class(N,jc,jj)+dltt*vol(k,i)
                go to 220
              else
                c_crit=c_crit-c_int(jac)
              end if
              if(jj.eq.numclass.and.c2(k,i,jac)*cmult(jac).lt.c_crit+c_int(jac))then
                c_class(N,jc,jj)=c_class(N,jc,jj)+dltt*vol(k,i)
              end if
            end do
          
220       continue
    

            end do


! Derived Constituents

            do jc=1,nacd_e

            jacd=cdn_e(jc)

            cd_tot(N,jc)=cd_tot(N,jc)+cd(k,i,jacd)*CDMULT(jacd)*vol(k,i)*dltt
            cd_cnt(N,jc)=cd_cnt(N,jc)+vol(k,i)*dltt
            cd_crit=cd_top(jacd)
            if(cd(k,i,jacd)*CDMULT(jacd).ge.cd_top(jacd))cd_class(N,jc,1)=cd_class(N,jc,1)+dltt*vol(k,i)
            do jj=2,numclass
              if(cd(k,i,jacd)*CDMULT(jacd).lt.cd_crit.and.cd(k,i,jacd)*CDMULT(jacd).ge.cd_crit-cd_int(jacd))then
                cd_class(N,jc,jj)=cd_class(N,jc,jj)+dltt*vol(k,i)
                go to 240
              else
                cd_crit=cd_crit-cd_int(jacd)
              end if
              if(jj.eq.numclass.and.cd(k,i,jacd)*CDMULT(jacd).lt.cd_crit+cd_int(jacd))then
                cd_class(N,jc,jj)=cd_class(N,jc,jj)+dltt*vol(k,i)
              end if
            end do
  240       continue

            end do
  ENDIF            
       end do
     end do
  end do
end do

! sum of volgL*dltt for volume fraction calculation

   sumvolt(N)=sumvolt(N)+VOLGL(N)*dltt
   ENDDO  ! END OF INTERVAL FOR SEGMENTS

650   continue

      if(iopenfish == 3)then
!  calculating average violation concentration and writing to file
       DO N=1,I_SEGINT
        cd_sum=0.0
        c_sum=0.0
        T_SUM=0.0
        V_SUM=0.0
        WRITE (I_INT,'(I1)') N   
        
        if(temp_vpr.eq.' ON')then
          if(t_cnt(N).gt.0.0)then
          t_avg=t_tot(N)/t_cnt(N)
          else
          t_avg=0.0
          end if

        open(CONE,file='envrprf_t_'//I_INT//'.csv',status='unknown')
        write(CONE,*)'"Temperature interval,","Fraction of volume"'
        temp_c=temp_top
          do i=1,numclass
          write(CONE,125)temp_c,t_class(N,i)/sumvolt(N)
          temp_c=temp_c-temp_int
          t_sum=t_sum+t_class(N,i)/sumvolt(N)
          end do
        write(CONE,'(1x)')
        write(CONE,'(" Sum of fractions, ",e12.4)')t_sum
        write(CONE,'(1x)')
        write(CONE,'(" Average, ",e12.4)')t_avg
        close(CONE)
        end if

        if(vel_vpr.eq.' ON')then
          if(v_cnt(N).gt.0.0)then
          v_avg=v_tot(N)/v_cnt(N)
          else
          v_avg=0.0
          end if
        open(CONE,file='envrprf_v'//I_INT//'.csv',status='unknown')
        write(CONE,*)'"Velocity interval,","Fraction of volume"'
        vel_c=vel_top
          do i=1,numclass
          write(CONE,125)vel_c,v_class(N,i)/sumvolt(N)
          vel_c=vel_c-vel_int
          v_sum=v_sum+v_class(N,i)/sumvolt(N)
          end do
        write(CONE,'(1x)')
        write(CONE,'(" Sum of fractions, ",e12.4)')v_sum
        write(CONE,'(1x)')
        write(CONE,'(" Average, ",e12.4)')v_avg
        close(CONE)
        end if
        
        
        if(depth_vpr.eq.' ON')then
          if(d_cnt(N).gt.0.0)then
          d_avg=d_tot(N)/d_cnt(N)
          else
          d_avg=0.0
          end if
        open(CONE,file='envrprf_depth'//I_INT//'.csv',status='unknown')
        write(CONE,*)'"Depth interval,","Fraction of time"'
        d_c=d_top
          do i=1,numclass
          write(CONE,125)d_c,d_class(N,i)/d_cnt(N)
          d_c=d_c-d_int
          d_sum=d_sum+d_class(N,i)/d_cnt(N)
          end do
        write(CONE,'(1x)')
        write(CONE,'(" Sum of fractions, ",e12.4)')d_sum
        write(CONE,'(1x)')
        write(CONE,'(" Average, ",e12.4)')d_avg
        close(CONE)
        end if
        
if(nac_e > 0)then
       open(CONE,file='envrprf_c'//I_INT//'.csv',status='unknown')
       write(CONE,4000)(cname2(cn_e(jc)),jc=1,nac_e)
4000 format(<nac_e>(a8,'_interval, Fraction_of_volume, '))
       do jc=1,nac_e
        if(c_cnt(N,jc).gt.0.0)then
        c_avg(jc)=c_tot(N,jc)/c_cnt(N,jc)
        else
        c_avg(jc)=0.0
        end if
          do i=1,numclass
             if(i.eq.1)then
             conc_c(jc,i)=c_top(cn_e(jc))
             else
             conc_c(jc,i)=conc_c(jc,i-1)-c_int(cn_e(jc))
             end if
          c_sum(jc)=c_sum(jc)+c_class(N,jc,i)/sumvolt(N)
          end do
       end do
        do i=1,numclass
        write(CONE,126)(conc_c(jc,i),c_class(N,jc,i)/sumvolt(N),jc=1,nac_e)
        end do
        write(CONE,'(1x)')
        write(CONE,'(<nac_e>("Sum_of_fractions,",f9.4,","))')(c_sum(jc),jc=1,nac_e)     
        write(CONE,'(1x)')
        write(CONE,'(<nac_e>("Average,",e12.4,","))')(c_avg(jc),jc=1,nac_e)     
        close(CONE)
      
  if(nacd_e > 0)then     
       open(CONE,file='envrprf_cd'//I_INT//'.csv',status='unknown')
     write(CONE,4001)(cdname2(cdn_e(jc)),jc=1,nacd_e)
4001 format(<nacd_e>(a8,'_interval, Fraction_of_volume,'))
       do jc=1,nacd_e
        if(cd_cnt(N,jc).gt.0.0)then
        cd_avg(jc)=cd_tot(N,jc)/cd_cnt(N,jc)
        else
        cd_avg(jc)=0.0
        end if

          do i=1,numclass
              if(i.eq.1)then
              conc_cd(jc,i)=cd_top(cdn_e(jc))
              else
              conc_cd(jc,i)=conc_cd(jc,i-1)-cd_int(cdn_e(jc))
              end if
          cd_sum(jc)=cd_sum(jc)+cd_class(N,jc,i)/sumvolt(N)
          end do
       end do
        do i=1,numclass
        write(CONE,124)(conc_cd(jc,i),cd_class(N,jc,i)/sumvolt(N),jc=1,nacd_e)
        end do
        write(CONE,'(1x)')
        write(CONE,'(<nacd_e>("Sum_of_fractions,",f9.4,","))')(cd_sum(jc),jc=1,nacd_e)
        write(CONE,'(1x)')
        write(CONE,'(<nacd_e>("Average,",e12.4,","))')(cd_avg(jc),jc=1,nacd_e)
        close(CONE)
    endif
endif
ENDDO
       DEallocate(c_cnt,cd_cnt,c_class,cd_class,c_tot,cd_tot,t_class,v_class,c_sum,cd_sum)
       DEallocate(conc_c,conc_cd)
       DEallocate(cc_e,c_int,c_top,cd_e,cd_int,cd_top,c_avg,cd_avg,cn_e,cdn_e)
       DEALLOCATE(T_TOT,T_CNT,V_TOT,V_CNT,D_TOT,D_CNT,D_CLASS,VOLGL,SUMVOLT)
       
124       format(<nacd_e>(f10.4,',',e12.4,','))
125       format((f6.2,',',e12.4,','))
126       format(<nac_e>(f10.4,',',e12.4,','))
127       format(<nacd_e>(f10.4,',',e12.4,','))
128       format(<nacd_e>(" 0, ",e12.4,','))
129       format(<nac_e>(" 0, ",e12.4,','))

    end if

timlast=jday

      END SUBROUTINE ENVIRP

