! Hg Simulation Module
! Control file: W2_Hg.csv 
!
! 1) Equilibrium partitioning of HgII and MeHg on DOC, POM, ISS, and Algae are modeled. 
! State variables: total concentrations of Hg0, HgII, MeHg in the water column 
!                  total concentrations of HgII, MeHg in the active sediment layer
! 5//2022: Add partitioning of HgII and MeHg on multi-zooplankton groups
!          Update methylation/demethylation with the total organic C (POC and DOC) decomposition.   
!
! 2) Non-Equilibrium (kinetic) partitioning of HgII on POM, ISS will be included
! 
!=========================================================================================================================== 
MODULE HgModule
  USE PREC,     ONLY: R8
  USE MAIN,     ONLY: CONHG, HG_CALC, use_BedHg, ORGC_CALC, ALG_CALC, SEDIMENT_CALC, SED_DIAG, RESTART_IN, QTRF, QPR, QINF, ATM_DEP_LOADING, QDT, &
                      KBTR, IWD, ITR, JBTR, KTTR, JBWD, KWD, NHg0, NHgII, NMeHg
  USE GLOBAL,   ONLY: JW, IU, ID, KT, KB, IMX, KMX, DLT, DAY, T1, C2, NSS, NAL, NZP, ICE, KTWB, CUS, BS, BE, DS, VOL, NWB, C1S
  USE LOGICC,   ONLY: PH_CALC, WITHDRAWALS, TRIBUTARIES,UP_FLOW, DN_FLOW, PRECIPITATION,DIST_TRIBS, HEAD_FLOW
  USE KINETIC,  ONLY: SSS, SS, POMS, LPOM, RPOM, LPOMC, RPOMC, ORGC, LPOMCD, LRPOMCD, RPOMCD, LPOMCHD, RPOMCHD, &
                      LDOMCD, RDOMCD, DOC, PH, AS, ALG, BETA, GAMMA, REAER, KDO, O2, DO2, SODD, SEDC, SEDDC, &
                      Hg0, HgII, MeHg, Hg0SS, HgIISS, MeHgSS, ANOX, DO4
  USE SURFHE,   ONLY: SHADE
  USE TVDC,     ONLY: SRON, QIN,    QTR, QOUT, CIN,    CTR,    CDTR,   CPR, QWD, QIND, CIND
  USE GEOMC,    ONLY: BI, BH1, BH2, DEPTHB, H1, H2, DLX, B
  USE SCREENC,  ONLY: JDAY, JWW, JTT
  USE RSTART,   ONLY: CUF
  USE SELWC,    ONLY: QSW, KTW, KBW
  USE CEMASedimentDiagenesis, ONLY: C2SF, KFSF, SD_kdiaPOC,SD_ThtaPOC
  USE ZOOPLANKTONC, ONLY: ZOO, ZS, ZOOPLANKTON_CALC
  
  ! sediment and partitioning options
  !logical  :: use_BedHg   
  logical  :: use_vb
  logical  :: use_DOCSorbed(2:3)
  logical  :: use_POMSorbed(2:3)
  logical  :: use_AnySolidSorbed(2:3) 
  logical  :: use_AnyAlgaeSorbed(2:3) 
  logical  :: use_AnyZooSorbed(2:3)
  logical  :: use_NonDOCSorbed(2:3)
  logical  :: use_NonPOMSorbed(2:3)
  logical  :: use_NonAnySolidSorbed(2:3) 
  logical  :: use_NonAnyAlgaeSorbed(2:3) 
  logical  :: use_NonAnyZooSorbed(2:3)
  logical  :: use_NonEquilibrium(2:3)
  !
  ! equilibrium partitioning parameters of HgII and MeHg
  real(R8), allocatable, dimension(:,:,:)   :: Kap          ! Algae partition coefficient (L/kg)
  real(R8), allocatable, dimension(:,:,:)   :: Kzp          ! Zooplankton partition coefficient (L/kg)
  real(R8), allocatable, dimension(:,:)     :: Kdoc, Kdoc2  ! DOC partition coefficient (L/kg)
  real(R8), allocatable, dimension(:,:)     :: Kpom, Kpom2  ! POM partition coefficient (L/kg)
  real(R8), allocatable, dimension(:,:,:)   :: Kp, kp2      ! ISS partition coefficient (L/kg)
  !
  ! transformation parameters 
  real(R8), allocatable, dimension(:)       :: Hg00                         ! Air concentration of Hg0 (ng/L)
  real(R8), allocatable, dimension(:)       :: vv_20, vv_theta              ! User-defined volatilization velocity (m/day)
  real(R8), allocatable, dimension(:)       :: k_oxid1, k_oxid2, k_oxid3    ! Photo-oxidation rate from Hg0 to HgII
  real(R8), allocatable, dimension(:)       :: k_reduct, kdoc_reduct        ! Photoreduction rate from HgII to Hg0 (m^2/J)
  real(R8), allocatable, dimension(:)       :: k_decay,  kdoc_decay         ! Photodegradation rate from MeHg to Hg0 (m^2/J)
  !
  real(R8), allocatable, dimension(:)       :: k_meth, k_meth2              ! Methylation rate from HgII to MeHg (m^3/gC)  
  real(R8), allocatable, dimension(:)       :: kdoc_meth, kdoc_meth2                                  
  real(R8), allocatable, dimension(:)       :: k_demeth,  k_demeth2         ! Demethylation rate from MeHg to HgII (m^3/gC)    
  real(R8), allocatable, dimension(:)       :: kdoc_demeth, kdoc_demeth2     
  real(R8), allocatable, dimension(:)       :: KDOHg1, KDOHg2, KDOn         ! DO inhibition function coefficients for methylation/demethylation
  !
  real(R8), allocatable, dimension(:)       :: BedH        ! Active sediment layer thickness (m)
  real(R8), allocatable, dimension(:)       :: BedPor      ! Porosity                      
  real(R8), allocatable, dimension(:)       :: BedPs       ! Density of bed sediments
  real(R8), allocatable, dimension(:)       :: vb          ! Sediment burial velocity (m/d) 
  real(R8), allocatable, dimension(:,:)     :: vbs         ! Sediment burial velocity (m/s) internally
  real(R8), allocatable, dimension(:,:)     :: vm1,vm      ! Sediment-water mass transfer velocity (cm/d) vm1 is input vm is used in formulae
  !
  ! water column
  real(R8), allocatable, dimension(:,:,:)   :: Hgd         ! Dissolved Hg in water (ng/L)
  real(R8), allocatable, dimension(:,:,:)   :: Hgdoc       ! DOC sorbed Hg in water (ng/L)
  real(R8), allocatable, dimension(:,:,:)   :: Hgpom       ! POM sorbed Hg in water (ng/L)
  real(R8), allocatable, dimension(:,:,:,:) :: Hgp         ! ISS sorbed Hg in water (ng/L) 
  real(R8), allocatable, dimension(:,:,:,:) :: Hgap        ! Algae sorbed Hg in water (ng/L) 
  real(R8), allocatable, dimension(:,:,:,:) :: Hgzp        ! Zooplankton sorbed Hg in water (ng/L) 
  real(R8), allocatable, dimension(:,:,:)   :: Hgzpt       ! Zooplankton sorbed Hg in water (ng/L)  
  real(R8), allocatable, dimension(:,:,:)   :: Hgzpts      ! Zooplankton sorbed Hg in water based on the mass (ng/g)    
  real(R8), allocatable, dimension(:,:,:)   :: Hgpt        ! TSS sorbed Hg in water (ng/L)  
  real(R8), allocatable, dimension(:,:,:)   :: Hgpts       ! TSS sorbed Hg in water based on the mass (ng/g)
  !
  ! sediment layer
  real(R8) :: dHg2dt(2:3)
  real(R8), allocatable, dimension(:,:,:)   :: Hg2         ! Hg concentration in sediment (ng/L)
  real(R8), allocatable, dimension(:,:,:)   :: Hgd2        ! Dissolved Hg in pore water (ng/L)
  real(R8), allocatable, dimension(:,:,:)   :: Hgdoc2      ! DOC sorbed Hg in pore water (ng/L)
  real(R8), allocatable, dimension(:,:,:)   :: Hgpom2      ! POM sorbed Hg in sediment (ng/L)
  real(R8), allocatable, dimension(:,:,:,:) :: Hgp2        ! ISS sorbed Hg in sediment (ng/L) 
  real(R8), allocatable, dimension(:,:,:)   :: Hgpt2       ! TSS sorbed Hg in sediment (ng/L)
  real(R8), allocatable, dimension(:,:,:)   :: Hgpts2      ! TSS sorbed Hg in sediment based on the mass (ng/g) 
  real(R8), allocatable, dimension(:,:)     :: Hg2_Ini     ! Initial condition of Hg in sediment (ng/L)
  !
  ! derived variables and pathway fluxes
  real(R8), allocatable, dimension(:,:,:)   :: CON_HG, CON_HG2
  real(R8), allocatable, dimension(:,:,:)   :: KF_HG, KF_HG2
  real(R8), allocatable, dimension(:,:)     :: KF_HG_SUM,KF_HG2_SUM
  real(R8) :: JDAYSTARTHG,NXTMFLHG,HGFLUXINT
  character(20), allocatable, dimension(:)  :: CNAMEHG, CNAMEHG2
  character(30) :: KFNAMEHG(10), KFNAMEHG2(6)
  !
  real(R8) :: MW(3)    
  real(R8), allocatable, dimension(:,:) :: Hg0_Vol
  real(R8), allocatable, dimension(:,:) :: Hg0_Oxidation
  !
  real(R8), allocatable, dimension(:,:) :: HgII_Settling
  real(R8), allocatable, dimension(:,:) :: HgII_Reduction
  real(R8), allocatable, dimension(:,:) :: HgII_Burial
  real(R8), allocatable, dimension(:,:) :: HgII_Transfer   ! Diffusive flux water column to sediment
  real(R8), allocatable, dimension(:,:) :: HgII_Methylation, HgII2_Methylation 
  !
  real(R8), allocatable, dimension(:,:) :: MeHg_Settling
  real(R8), allocatable, dimension(:,:) :: MeHg_Decay  
  real(R8), allocatable, dimension(:,:) :: MeHg_Burial
  real(R8), allocatable, dimension(:,:) :: MeHg_Transfer   ! Diffusive flux water column to sediment
  real(R8), allocatable, dimension(:,:) :: MeHg_Demethylation, MeHg2_Demethylation 
  !
  integer  :: nRegion, nRegion_ini, FLUX_OUTPUT_FILE=8100, SREGION_INI          
  integer,  allocatable, dimension(:)     :: RegionS, RegionE, SREGION,EREGION              
  integer,  allocatable, dimension(:)     :: RegionS_ini, RegionE_ini         
  integer,  allocatable, dimension(:)     :: SegRegion             
  real(R8), allocatable, dimension(:)     :: dBedISSdt
  real(R8), allocatable, dimension(:,:)   :: BedISS_Ini,JPOC2, JPOC 
  real(R8), allocatable, dimension(:,:,:) :: BedISS
  !
  ! local variables  
  real(R8) :: ISS_Settling
  real(R8) :: POM, Bed_DOC, Bed_POM 
  real(R8) :: LIGHT                !, JPOC, JPOC2
  real(R8), allocatable, dimension(:,:) :: HgII_Bed, MeHg_Bed
  real(R8), allocatable, dimension(:,:) :: HgII2_Transfer, MeHg2_Transfer
  real(R8) :: KH_Hg, Hg_Value, xx
  real(R8),allocatable, dimension(:,:,:) :: Cd2, Cdoc2    !Cd2(2:3,1:KMX,1:IMX), Cdoc2(2:3,1:KMX,1:IMX)
  integer  :: I, K, JS, JA, JZ, r, i23, ia(2:3)
  logical  :: SkipLoop
  logical  :: IsWaterCell 
  character(256) :: Comments
  integer :: NSSmax,NALmax,NZPmax,idebug
  real(R8), allocatable, dimension(:) :: ATMDEP_HG0,ATMDEP_HG2,ATMDEP_MEHG
  real(R8) ::  THgOUT,THgIN,THgPR,THgWD,THgDTRIB,THgTRIB,MeHgOUT,MeHgIN,MeHgPR,MeHgWD,MeHgDTRIB,MeHgTRIB
  !real(R8) :: HgII_pom, HgII_ap, HgII_iss, HgII_zp, HgII_pomb, HgII_apb, HgII_issb, HgII_zpb  ! for debug
  !
  contains
  
  !===========================================================================================================================
  ! read Hg parameters
  subroutine HgModule_Input
    implicit none
      
	  Open(CONHG, file = "W2_Hg.csv", status='old')
	  SkipLoop = .FALSE.
	  Do While(.NOT. SkipLoop)
		  Read(CONHG,'(a)') Comments
		  If(index(Comments, "$") == 0) SkipLoop = .TRUE.
	  End Do
	  BackSpace(CONHG)
       
    ! global parameters
    read(CONHG,*) Comments,idebug
    
    read(CONHG,*) Comments, use_BedHg 
    read(CONHG,*) Comments, (ia(i23), i23=2,3)
    do i23 = 2,3 
      if (ia(i23) == 1) then
        use_DOCSorbed(i23) = .true.
      else
        use_NonDOCSorbed(i23) = .true.
      end if
    end do   
    !
    read(CONHG,*) Comments, (ia(i23), i23=2,3)
    do i23 = 2,3 
      if (ia(i23) == 1) then
        use_POMSorbed(i23) = .true.
      else
        use_NonPOMSorbed(i23) = .true.
      end if
    end do 
    !
    read(CONHG,*) Comments, (ia(i23), i23=2,3)
    do i23 = 2,3 
      if (ia(i23) == 1) then
        use_AnySolidSorbed(i23) = .true.
      else
        use_NonAnySolidSorbed(i23) = .true.
      end if
    end do
    !
    read(CONHG,*) Comments, (ia(i23), i23=2,3)
    do i23 = 2,3 
      if (ia(i23) == 1) then
        use_AnyAlgaeSorbed(i23) = .true.
      else
        use_NonAnyAlgaeSorbed(i23) = .true.
      end if
    end do   
    !
    read(CONHG,*) Comments, (ia(i23), i23=2,3)
    do i23 = 2,3 
      if (ia(i23) == 1) then
        use_AnyZooSorbed(i23) = .true.
      else
        use_NonAnyZooSorbed(i23) = .true.
      end if
    end do 
    !
    !do i23 = 2,3 
    !  if (use_NonDOCSorbed(i23) .or. use_NonPOMSorbed(i23) .or. use_NonAnySolidSorbed(i23) .or. &
    !    use_NonAnyAlgaeSorbed(i23) .or. use_NonAnyAlgaeSorbed(i23)) then
    !    use_NonEquilibrium(2:3) = .true.
    !  end if
    !end do
    !  
    read(CONHG,*) Comments, nRegion
    allocate(RegionS(nRegion), RegionE(nRegion))
    read(CONHG,*) Comments, (RegionS(r), r=1,nRegion)
    read(CONHG,*) Comments, (RegionE(r), r=1,nRegion)
    
    call InitializeHgModule
    
    ! water column
    ! Hg0
    read(CONHG,*) Comments
    read(CONHG,*) Comments, (Hg00(r),     r=1,nRegion)
    ! vv_20 < 0.0 --> O2 based rearation rate is used.
    read(CONHG,*) Comments, (vv_20(r),    r=1,nRegion)
    read(CONHG,*) Comments, (vv_theta(r), r=1,nRegion)
    read(CONHG,*) Comments, (k_oxid1(r),  r=1,nRegion)
    read(CONHG,*) Comments, (k_oxid2(r),  r=1,nRegion)
    read(CONHG,*) Comments, (k_oxid3(r),  r=1,nRegion)
    ! HgII
    read(CONHG,*) Comments,  MW(2)   
    read(CONHG,*) Comments, (Kdoc(2,r),   r=1,nRegion)  
    read(CONHG,*) Comments, (Kpom(2,r),   r=1,nRegion)
    ! partitioning coefficients for individual groups 
    do JS=1,NSSmax
    read(CONHG,*) Comments, (Kp(2,JS,r),r=1,nRegion)
    enddo
    do JA=1,NALmax
    read(CONHG,*) Comments, (Kap(2,JA,r),r=1,nRegion) 
    enddo
    do JZ=1,NZPmax
      read(CONHG,*) Comments, (Kzp(2,JZ,r),r=1,nRegion) 
    enddo
    !
    read(CONHG,*) Comments, (k_reduct(r),    r=1,nRegion)    
    read(CONHG,*) Comments, (kdoc_reduct(r), r=1,nRegion)     
    read(CONHG,*) Comments, (k_meth(r),      r=1,nRegion)
    read(CONHG,*) Comments, (kdoc_meth(r),   r=1,nRegion) 
    read(CONHG,*) Comments, (vm1(2,r),        r=1,nRegion) 
    ! MeHg
    read(CONHG,*) Comments, MW(3)  
    read(CONHG,*) Comments, (Kdoc(3,r),      r=1,nRegion)
    read(CONHG,*) Comments, (Kpom(3,r),      r=1,nRegion)
    ! 
    do JS=1,NSSmax
    read(CONHG,*) Comments, (Kp(3,JS,r),r=1,nRegion)
    enddo
    do JA=1,NALmax
    read(CONHG,*) Comments, (Kap(3,JA,r),r=1,nRegion)
    enddo
    do JZ=1,NZPmax
      read(CONHG,*) Comments, (Kzp(3,JZ,r),r=1,nRegion) 
    enddo
    !
    read(CONHG,*) Comments, (k_decay(r),     r=1,nRegion)
    read(CONHG,*) Comments, (kdoc_decay(r),  r=1,nRegion)
    read(CONHG,*) Comments, (k_demeth(r),    r=1,nRegion)
    read(CONHG,*) Comments, (kdoc_demeth(r), r=1,nRegion)
    read(CONHG,*) Comments, (vm1(3,r),        r=1,nRegion)
    !
    read(CONHG,*) Comments, (KDOHg1(r),      r=1,nRegion)
    read(CONHG,*) Comments, (KDOHg2(r),      r=1,nRegion)
    read(CONHG,*) Comments, (KDOn(r),        r=1,nRegion)
    !
    ! sediment layer
    read(CONHG,*) Comments
    if (use_BedHg) then
      read(CONHG,*) Comments, use_vb
      read(CONHG,*) Comments, (vb(r),         r=1,nRegion)
      vb=vb/DAY
      read(CONHG,*) Comments, (BedH(r),       r=1,nRegion)
      read(CONHG,*) Comments, (BedPor(r),     r=1,nRegion)  
      read(CONHG,*) Comments, (BedPs(r),      r=1,nRegion) 
      read(CONHG,*) Comments, (Kdoc2(2,r),    r=1,nRegion)      
      read(CONHG,*) Comments, (Kpom2(2,r),    r=1,nRegion)
      do js=1,nssmax
      read(CONHG,*) Comments, (Kp2(2,JS,r),r=1,nRegion) 
      enddo
      read(CONHG,*) Comments, (k_meth2(r),    r=1,nRegion)    
      read(CONHG,*) Comments, (kdoc_meth2(r), r=1,nRegion)  
      !
      read(CONHG,*) Comments, (Kdoc2(3,r),    r=1,nRegion)      
      read(CONHG,*) Comments, (Kpom2(3,r),    r=1,nRegion)  
      do js=1,nssmax
      read(CONHG,*) Comments, (Kp2(3,JS,r),r=1,nRegion) 
      enddo
      read(CONHG,*) Comments, (k_demeth2(r),  r=1,nRegion)    
      read(CONHG,*) Comments, (kdoc_demeth2(r), r=1,nRegion)  
        ELSE
      read(CONHG,*) Comments
      read(CONHG,*) Comments
      read(CONHG,*) Comments
      read(CONHG,*) Comments 
      read(CONHG,*) Comments
      read(CONHG,*) Comments     
      read(CONHG,*) Comments
          do js=1,nssmax
      read(CONHG,*) Comments
          enddo
      read(CONHG,*) Comments
      read(CONHG,*) Comments
      read(CONHG,*) Comments     
      read(CONHG,*) Comments 
          do js=1,nssmax
      read(CONHG,*) Comments
          enddo
      read(CONHG,*) Comments   
      read(CONHG,*) Comments     
    end if   
    !
    allocate(SegRegion(IMX))
    SegRegion = 1
    do r = 1,nRegion
      do I = RegionS(r),RegionE(r)
        SegRegion(I) = r
      end do
    end do
    
    
  end subroutine
  !===========================================================================================================================
  ! read initial conditions of HgII and MeHg and ISS in the sediment layer AND FLUX OUTPUT PARAMETERS
  subroutine HgModule_IC
    implicit none
    
    allocate(BedISS(KMX,IMX,NSSmax), dBedISSdt(NSSmax))
    
    if (.not. RESTART_IN) then 
      read(CONHG,*) Comments
      read(CONHG,*) Comments, nRegion_ini  
      allocate(RegionS_ini(nRegion_ini), RegionE_ini(nRegion_ini)) 
      allocate(Hg2_Ini(2:3, nRegion_ini))
      allocate(BedISS_Ini(NSSmax,nRegion_ini))
      !
      read(CONHG,*) Comments, (RegionS_ini(r), r=1,nRegion_ini)
      read(CONHG,*) Comments, (RegionE_ini(r), r=1,nRegion_ini)   
      read(CONHG,*) Comments, (Hg2_Ini(2,r), r=1,nRegion_ini)
      read(CONHG,*) Comments, (Hg2_Ini(3,r), r=1,nRegion_ini)
      do js=1,nssmax
      read(CONHG,*) Comments, (BedISS_Ini(JS,r), r=1,nRegion_ini)
      enddo
      
      do r = 1,nRegion_ini      
        do I = RegionS_ini(r),RegionE_ini(r)
          Hg2(2,1,I)  = Hg2_Ini(2,r)
          Hg2(3,1,I)  = Hg2_Ini(3,r)
          BedISS(1,I,:) = BedISS_Ini(:,r)
        end do    
      end do
      do K = 2,KMX
        Hg2(2,K,:)  = Hg2(2,1,:)
        Hg2(3,K,:)  = Hg2(3,1,:)
        BedISS(K,:,:) = BedISS(1,:,:)
      end do
      
      ! SEDIMENT FLUX OUTPUT REGIONS
      read(CONHG,*) Comments
      read(CONHG,*) Comments, SREGION_INI  
      allocate(SRegion(SRegion_ini),ERegion(SRegion_ini)) 
      read(CONHG,*) Comments, HGFLUXINT
      read(CONHG,*) Comments, (SREGION(JS),JS=1,SREGION_INI)
      read(CONHG,*) Comments, (EREGION(JS),JS=1,SREGION_INI)
    end if
    
    close(CONHG)   
  end subroutine  
  !===========================================================================================================================
  ! allocate and initialize Hg parameters and variables
  subroutine InitializeHgModule
    implicit none
    
    NSSmax=max(5,NSS)
    NALmax=max(5,NAL)
    NZPmax=max(5,NZP)
    
    ! Hg parameters
    if(allocated(Kap)) deallocate(Kap)
    allocate(Kap(2:3,NALmax,nRegion))
    Kap = 2.0E5
    !    
    if(allocated(Kzp)) deallocate(Kzp)
    allocate(Kzp(2:3,NZPmax,nRegion))
    Kzp = 2.0E5
    !
    if(allocated(Kdoc)) deallocate(Kdoc) 
    allocate(Kdoc(2:3,nRegion))
    Kdoc = 2.0E5
    if (use_BedHg) then
      if(allocated(Kdoc2)) deallocate(Kdoc2)
      allocate(Kdoc2(2:3,nRegion))
      Kdoc2 = 2.0E5
    end if
    !
    if(allocated(Kpom)) deallocate(Kpom) 
    allocate(Kpom(2:3,nRegion))
    Kpom = 2.0E5  
    if (use_BedHg) then
      if(allocated(Kpom2)) deallocate(Kpom2) 
      allocate(Kpom2(2:3,nRegion))
      Kpom2 = 2.0E5
    end if
    !
    if(allocated(Kp)) deallocate(Kp)
    allocate(Kp(2:3,NSSmax,nRegion))
    Kp = 2.0E5
    if (use_BedHg) then     
      if(allocated(Kp2)) deallocate(Kp2)
      allocate(Kp2(2:3,NSSmax,nRegion))
      Kp2 = 2.0E5      
    end if
    !
    if(allocated(vm1)) deallocate(vm1)
    allocate(vm1(2:3, nRegion))
    vm1 = 0.01 
    
    if(allocated(vm)) deallocate(vm)
    allocate(vm(2:3, nRegion));vm=0.01
    !
    call alloc_R8(nRegion, Hg00, 0.1d0)
    call alloc_R8(nRegion, vv_20, 0.144d0)
    call alloc_R8(nRegion, vv_theta, 1.024d0)
    call alloc_R8(nRegion, k_oxid1, 0.01d0)
    call alloc_R8(nRegion, k_oxid2, 0.0d0)
    call alloc_R8(nRegion, k_oxid3, 0.01d0)
    !  
    call alloc_R8(nRegion, k_reduct,    0.05d0)
    call alloc_R8(nRegion, kdoc_reduct, 0.0d0)
    call alloc_R8(nRegion, k_decay,    0.01d0)
    call alloc_R8(nRegion, kdoc_decay, 0.0d0)     
    call alloc_R8(nRegion, k_meth,    0.001d0)
    call alloc_R8(nRegion, kdoc_meth, 0.0d0)
    call alloc_R8(nRegion, k_demeth,    0.05d0)
    call alloc_R8(nRegion, kdoc_demeth, 0.0d0) 
    call alloc_R8(nRegion, KDOHg1,      1.0d0)
    call alloc_R8(nRegion, KDOHg2,      6.0d0)
    call alloc_R8(nRegion, KDOn,        3.5d0)
        
    if (use_BedHg) then
      call alloc_R8(nRegion, BedH,   0.1d0)
      call alloc_R8(nRegion, BedPor, 0.5d0)
      call alloc_R8(nRegion, BedPs,  2.7d0)
      call alloc_R8(nRegion, vb, 0.0025d0)
      call alloc_R8(nRegion, k_meth2,    0.001d0)
      call alloc_R8(nRegion, kdoc_meth2, 0.0d0)
      call alloc_R8(nRegion, k_demeth2,  0.02d0)
      call alloc_R8(nRegion, kdoc_demeth2, 0.0d0)         
    end if 
    !
    MW(1) = 200.59
    MW(2) = 271.52
    MW(3) = 230.66    
    !
    ! derived variables and pathway fuxes in the water column    
    allocate(Hgd(2:3, KMX,IMX), Hgdoc(2:3, KMX,IMX), Hgpom(2:3, KMX,IMX))
    allocate(Hgap(2:3,NAL, KMX,IMX), Hgp(2:3,NSS, KMX,IMX))
    allocate(Hgpt(2:3, KMX,IMX), Hgpts(2:3, KMX,IMX))
    ! Zooplankton
    allocate(Hgzp(2:3,NZP, KMX,IMX))
    allocate(Hgzpt(2:3, KMX,IMX), Hgzpts(2:3, KMX,IMX))    
    !
    allocate(CON_HG(11+2*(NSS+NAL+NZP),KMX,IMX))
    allocate(CNAMEHG(11+2*(NSS+NAL+NZP)), KF_HG(15, KMX,IMX),KF_HG_SUM(15,IMX))
    allocate(Hg0_Vol(KMX,IMX),Hg0_Oxidation(KMX,IMX),HgII_Reduction(KMX,IMX),HgII_Methylation(KMX,IMX),MeHg_Burial(KMX,IMX))
    allocate(HgII_Transfer(KMX,IMX),MeHg_Demethylation(KMX,IMX),MeHg_Decay(KMX,IMX),MeHg_Settling(KMX,IMX),MeHg_Transfer(KMX,IMX))
    
    Hg0_Oxidation=0.0;Hg0_Vol=0.0;HgII_Reduction=0.0;HgII_Methylation=0.0;HgII_Transfer=0.0;HgII_Settling=0.0
    MeHg_Demethylation=0.0;MeHg_Decay=0.0;MeHg_Settling=0.0;MeHg_Burial=0.0
    
    allocate(HgII2_Methylation(KMX,IMX),JPOC(KMX,IMX),JPOC2(KMX,IMX),HgII_Burial(KMX,IMX),MeHg2_Demethylation(KMX,IMX))
    HgII2_Methylation=0.0;JPOC=0.0;JPOC2=0.0;HgII_Burial=0.0;MeHg2_Demethylation=0.0
    KF_HG_SUM=0.0;KF_HG=0.0; CON_HG=0.0
    
    allocate(MeHg_Bed(KMX,IMX),HgII_Bed(KMX,IMX),HgII2_Transfer(KMX,IMX),MeHg2_Transfer(KMX,IMX),HgII_Settling(KMX,IMX))
    
    MeHg_Bed=0.0;HgII_Bed=0.0;HgII2_Transfer=0.0;MeHg2_Transfer=0.0;MeHg_Transfer=0.0
    
    allocate(ATMDEP_HG0(15),ATMDEP_HG2(15),ATMDEP_MEHG(15))
    
    allocate(Cd2(2:3,KMX,IMX), Cdoc2(2:3,KMX,IMX))
    Cd2=0.0;Cdoc2=0.0
    ATMDEP_HG0=0.0;ATMDEP_HG2=0.0;ATMDEP_MEHG=0.0
    THgOUT=0.0;THgIN=0.0;THgPR=0.0;THgWD=0.0;THgDTRIB=0.0;THgTRIB=0.0
    MeHgOUT=0.0;MeHgIN=0.0;MeHgPR=0.0;MeHgWD=0.0;MeHgDTRIB=0.0;MeHgTRIB=0.0
    !
    ! HgII
    CNAMEHG(1) = 'HgIId(ng/L)'
    CNAMEHG(2) = 'HgIIdoc(ng/L)'
    CNAMEHG(3) = 'HgIIpom(ng/L)'
    do JS = 1,NSS
      CNAMEHG(3+JS) = trim(addindex('HgIIp',JS)) //'(ng/L)'
    end do
    do JA = 1,NAL
      CNAMEHG((NSS+3)+JA) = trim(addindex('HgIIap',JA)) //'(ng/L)'
    end do
    ! Zooplankton
    do JZ = 1,NZP
      CNAMEHG((NSS+NAL+3)+JZ) = trim(addindex('HgIIzp',JZ)) //'(ng/L)'
    end do      
    CNAMEHG((NSS+NAL+NZP)+4) = 'HgIIpt(ng/L)'
    CNAMEHG((NSS+NAL+NZP)+5) = 'HgIIpts(ng/g)'  
    !
    ! MeHg
    CNAMEHG((NSS+NAL+NZP)+6) = 'MeHgd(ng/L)'
    CNAMEHG((NSS+NAL+NZP)+7) = 'MeHgdoc(ng/L)'
    CNAMEHG((NSS+NAL+NZP)+8) = 'MeHgpom(ng/L)' 
    do JS = 1,NSS
      CNAMEHG((NSS+NAL+NZP)+8+JS) = trim(addindex('MeHgp',JS)) //'(ng/L)'
    end do
    do JA = 1,NAL
      CNAMEHG(2*NSS+NAL+NZP+8+JA) = trim(addindex('MeHgap',JA)) //'(ng/L)'
    end do   
    ! Zooplankton
    do JZ = 1,NZP
      CNAMEHG(2*(NSS+NAL)+NZP+8+JZ) = trim(addindex('MeHgzp',JZ)) //'(ng/L)'
    end do     
    CNAMEHG(2*(NSS+NAL+NZP)+9) = 'MeHgpt(ng/L)'
    CNAMEHG(2*(NSS+NAL+NZP)+10) = 'MeHgpts(ng/g)'
    !
    CNAMEHG(2*(NSS+NAL+NZP)+11) = 'THg(ng/L)'
    CNAMEHG = ADJUSTR(CNAMEHG)
    !
    ! pathway fluxes
    KFNAMEHG(1) = 'Hg0 volatilization(ng/L/d)'
    KFNAMEHG(2) = 'Hg0 oxidation(ng/L/d)'
    KFNAMEHG(3) = 'HgII reduction(ng/L/d)'
    KFNAMEHG(4) = 'HgII methylation(ng/L/d)'
    KFNAMEHG(5) = 'HgII settling(ng/L/d)'
    KFNAMEHG(6) = 'HgII sediment transfer(ng/L/d)' 
    KFNAMEHG(7) = 'MeHg degradation(ng/L/d)'
    KFNAMEHG(8) = 'MeHg demethylation(ng/L/d)'
    KFNAMEHG(9) = 'MeHg settling(ng/L/d)'
    KFNAMEHG(10) = 'MeHg sediment transfer(ng/L/d)'
    KFNAMEHG = ADJUSTR(KFNAMEHG)
    !
    ! state variables, derived vaiables and pathway fluxes in the sediment layer
    if (use_BedHg) then
      allocate(Hg2(2:3,    KMX,IMX), Hgd2(2:3, KMX,IMX), Hgdoc2(2:3, KMX,IMX))
      allocate(Hgpom2(2:3, KMX,IMX), Hgp2(2:3,NSSmax, KMX,IMX))
      allocate(Hgpt2(2:3,  KMX,IMX), Hgpts2(2:3, KMX,IMX))
      allocate(CON_HG2(13+3*NSS, KMX,IMX))
      allocate(CNAMEHG2(13+3*NSS), KF_HG2(6, KMX,IMX),KF_HG2_SUM(6,IMX))  
      allocate(vbs(KMX,IMX))
      
      KF_HG2_SUM=0.0;CON_HG2=0.0;KF_HG2=0.0
      !
      CNAMEHG2(1) = 'BedHgII(ng/L)'
      CNAMEHG2(2) = 'BedHgIId(ng/L)' 
      CNAMEHG2(3) = 'BedHgIIdoc(ng/L)'
      CNAMEHG2(4) = 'BedHgIIpom(ng/L)'
      do JS = 1,NSS
        CNAMEHG2(4+JS) = trim(addindex('BedHgIIp',JS)) //'(ng/L)'
      end do
      CNAMEHG2(NSS+5) = 'BedHgIIpt(ng/L)'
      CNAMEHG2(NSS+6) = 'BedHgIIpts(ng/g)'  
      !
      CNAMEHG2(NSS+7) = 'BedMeHg(ng/L)'
      CNAMEHG2(NSS+8) = 'BedMeHgd(ng/L)'
      CNAMEHG2(NSS+9) = 'BedMeHgdoc(ng/L)'
      CNAMEHG2(NSS+10) = 'BedMeHgpom(ng/L)'
      do JS = 1,NSS
        CNAMEHG2(NSS+10+JS) = trim(addindex('BedMeHgp',JS)) //'(ng/L)'
      end do
      CNAMEHG2(2*NSS+11) = 'BedMeHgpt(ng/L)'
      CNAMEHG2(2*NSS+12) = 'BedMeHgpts(ng/g)'
      CNAMEHG2(2*NSS+13) = 'BedTHg(ng/L)'
      !
      do JS = 1,NSS
        CNAMEHG2(2*NSS+13+JS) = trim(addindex('BedISS',JS)) //'(mg/L)'
      end do
      CNAMEHG2 = ADJUSTR(CNAMEHG2)
      !
      KFNAMEHG2(1) = 'BedHgII methylation(ng/L/d)'
      KFNAMEHG2(2) = 'BedHgII burial(ng/L/d)'
      KFNAMEHG2(3) = 'BedHgII settling accumulation(ng/L/d)'      
      KFNAMEHG2(4) = 'BedMeHg demethylation(ng/L/d)'
      KFNAMEHG2(5) = 'BedMeHg burial(ng/L/d)'      
      KFNAMEHG2(6) = 'BedMeHg settling accumulation(ng/L/d)'
      KFNAMEHG2 = ADJUSTR(KFNAMEHG2)
    end if
  end subroutine
  !===========================================================================================================================
  ! compute sediment burial velocity according to the solid mass balance for the active sediment layer
  subroutine ComputeBedVb   
    real(R8) :: TISS(KMX,IMX), TISS2(KMX,IMX)
    
    DO I=IU,ID
      r = SegRegion(I)

    if(use_vb) then
        vbs(:,I) = vb(r)     ! m/s -- conversion from m/d to m/s made above
    else
        
      DO K=KT,KB(I)     
        TISS(K,I) = 0.0
        TISS2(K,I) = 0.0
        do JS = 1,NSS
          TISS(K,I) = TISS(K,I) + SSS(JS) * SS(K,I,JS)
          TISS2(K,I) = TISS2(K,I) + BedISS(K,I,JS)
        end do 
        !vbs(K,I) = TISS(K,I) / (1.0E6*TISS2(K,I)*(1.0-BedPor(r)))    ! m/d  Error in 1E6 since BedISS is in mg/l or g/m3 not ng/l also porosity not required  
        vbs(K,I) = TISS(K,I) / TISS2(K,I)                             !(1.0E6*BedPs(r)*(1.0-BedPor(r)))            ! m/s since BedISS=rhos*(1-porosity)

        !
        !if(use_vb) vbs(K,I) = vb(r)
      END DO
    endif
    END DO 
    
  end subroutine
  !===========================================================================================================================
  ! compute concentrations of ISS (1 to NSS) in the sediment layer
  !ZZ Note: sediment resuspension is not included.
  subroutine ComputeBedISS
    
    call ComputeBedVb 
    
    DO I=IU,ID
      r = SegRegion(I)     
      DO K=KT,KB(I)             
        do JS = 1,NSS
          !dBedISSdt(JS) = (SSS(JS)*SS(K,I,JS) - vb(r)/DAY*BedISS(K,I,JS)) / BedH(r)
          dBedISSdt(JS) = (SSS(JS)*SS(K,I,JS) - vbs(K,I)*BedISS(K,I,JS)) / BedH(r)        ! SSS in m/s  vbs in m/s
          BedISS(K,I,JS) = MAX((dBedISSdt(JS)*DLT + BedISS(K,I,JS)), 0.0)   
        end do 
      END DO
    END DO
  end subroutine
  !===========================================================================================================================
  ! call Hg kinetics subroutines  
  subroutine ComputeHgKinetics
    if(use_BedHg) call ComputeBedISS
    !
    call HgPartition
    call HgPathways
    call HgEquilibriumKinetics
    if(use_BedHg) call ComputeBedHg  
    call HgOutputs           
  end subroutine
  !===========================================================================================================================
  ! Hg partition
  subroutine HgPartition
    implicit none
    !
    integer :: j12
    
    DO I=IU,ID
      r = SegRegion(I)
      DO K=KT,KB(I) 
        !
        do i23 = 2, 3
          do j12 = 1, 2
            if (j12 == 1) then    
              IsWaterCell = .true.
            else
              IsWaterCell = .false.
            end if
            !
            if (i23 == 2) then
              Hg_Value = HgII(K,I)
            else if (i23 == 3) then
              Hg_Value = MeHg(K,I)
            end if
            !
            call EquilibriumPartition
            !if (use_NonEquilibrium(i23)) then
            !  call NonEquilibriumPartition
            !else  
            !  call EquilibriumPartition 
            !end if  
          end do      
        end do
        !
      END DO
    END DO
  end subroutine
  !===========================================================================================================================
  ! equilibrium partitioning of HgII and MeHg
  subroutine EquilibriumPartition
    implicit none
    !
    real(R8) :: Rd
          
    !ZZ Note: Sediment DOC is not modeled now.
    Bed_DOC = DOC(K,I)
    !
    IF (ORGC_CALC) THEN
      POM = (LPOMC(K,I)+RPOMC(K,I))/ORGC(JW)
    ELSE
      POM = LPOM(K,I)+RPOM(K,I)
    END IF
    ! 
    if (use_BedHg) then
      if (SEDIMENT_CALC(JW)) then
        ! First-order sediment compartment
        Bed_POM = SEDC(K,I) / ORGC(JW)        
      else if (SED_DIAG =='      ON') then
        Bed_POM = (C2SF(K,I,19)+C2SF(K,I,20)+C2SF(K,I,21)) / ORGC(JW)
      else
        Bed_POM = 0.0
      end if
    end if
    !
    ! water column
    if (IsWaterCell) then 
      Rd = 1.0
      if(use_DOCSorbed(i23))  Rd = Rd + Kdoc(i23,r)*DOC(K,I)/1.0E6
      if(use_POMSorbed(i23))  Rd = Rd + Kpom(i23,r)*POM/1.0E6
      if (use_AnySolidSorbed(i23)) then
        do JS = 1,NSS
          Rd = Rd + Kp(i23,JS,r)*SS(K,I,JS)/1.0E6
        end do
      end if
      if (use_AnyAlgaeSorbed(i23)) then
        do JA = 1,NAL
          if(ALG_CALC(JA)) Rd = Rd + Kap(i23,JA,r)*ALG(K,I,JA)/1.0E6
        end do
      end if
      ! Zooplankton
      if (use_AnyZooSorbed(i23)) then
        IF (ZOOPLANKTON_CALC) THEN
          do JZ = 1,NZP
            Rd = Rd + Kzp(i23,JZ,r)*ZOO(K,I,JZ)/1.0E6
          end do
        END IF
      end if      
      !
      Hgd(i23,K,I) = Hg_Value/Rd
      if (use_DOCSorbed(i23)) then
        Hgdoc(i23,K,I) = Kdoc(i23,r)*DOC(K,I)/1.0E6/Rd * Hg_Value
      else
        Hgdoc(i23,K,I) = 0.0
      end if
      if (use_POMSorbed(i23)) then
        Hgpom(i23,K,I) = Kpom(i23,r)*POM/1.0E6/Rd * Hg_Value
      else
        Hgpom(i23,K,I) = 0.0
      end if
      !
      do JS = 1,NSS
        if (use_AnySolidSorbed(i23)) then
          Hgp(i23,JS,K,I) = Kp(i23,JS,r)*SS(K,I,JS)/1.0E6/Rd * Hg_Value
        else
          Hgp(i23,JS,K,I) = 0.0
        end if
      end do       
      !
      if (use_AnyAlgaeSorbed(i23)) then
        do JA = 1,NAL
          if (ALG_CALC(JA)) then
            Hgap(i23,JA,K,I) = Kap(i23,JA,r)*ALG(K,I,JA)/1.0E6/Rd * Hg_Value
          else
            Hgap(i23,JA,K,I) = 0.0
          end if
        end do
      end if
      ! Zooplankton
      if (use_AnyZooSorbed(i23)) then
        IF (ZOOPLANKTON_CALC) THEN
          do JZ = 1,NZP
            Hgzp(i23,JZ,K,I) = Kzp(i23,JZ,r)*ZOO(K,I,JZ)/1.0E6/Rd * Hg_Value
          end do
        END IF
      end if
      !
    else if (use_BedHg) then
      Rd = BedPor(r) 
      if(use_DOCSorbed(i23)) Rd = Rd + Kdoc2(i23,r)*Bed_DOC/1.0E6*BedPor(r) 
      if(use_POMSorbed(i23)) Rd = Rd + Kpom2(i23,r)*Bed_POM/1.0E6
      if (use_AnySolidSorbed(i23)) then
        do JS = 1,NSS
          Rd = Rd + Kp2(i23,JS,r)*BedISS(K,I,JS)/1.0E6
        end do
      end if
      !
      Cd2(i23,K,I) = BedPor(r)/Rd * Hg2(i23,K,I)
      if (use_DOCSorbed(i23)) then
        Cdoc2(i23,K,I) = Kdoc2(i23,r)*Bed_DOC/1.0E6*BedPor(r)/Rd * Hg2(i23,K,I)
      else
        Cdoc2(i23,K,I) = 0.0
      end if
      if (use_POMSorbed(i23)) then
        Hgpom2(i23,K,I) = Kpom2(i23,r)*Bed_POM/1.0E6/Rd * Hg2(i23,K,I)
      else
        Hgpom2(i23,K,I) = 0.0
      end if
      !
      do JS = 1,NSS
        if (use_AnySolidSorbed(i23)) then
          Hgp2(i23,JS,K,I) = Kp2(i23,JS,r)*BedISS(K,I,JS)/1.0E6/Rd * Hg2(i23,K,I)
        else
          Hgp2(i23,JS,K,I) = 0.0
        end if
      end do          
    end if 
  end subroutine
  !===========================================================================================================================
  ! kinetic (non-equilibrium) partitioning of HgII 
  subroutine NonEquilibriumPartition  
  
  
  end subroutine
  !================================================================================================================================  
  ! Hg pathway fluxes
  subroutine HgPathways
    implicit none
    !
    real(R8) :: DGMratio, DGMratiot
    real(R8) :: BIBH2(KMX,IMX)
    real(R8) :: LAM1, LAM2
    
    DO I = IU,ID
      r = SegRegion(I)
      LIGHT = (1.0-BETA(JW))*SRON(JW)*SHADE(I)                
      LAM1  = LIGHT
      LAM2  = LIGHT
      DO K = KT,KB(I)
        LAM1  = LAM2
        LAM2  = LAM1*EXP(-GAMMA(K,I)*H2(K,I))
        LIGHT = LAM1*(1.-EXP(-GAMMA(K,I)*H2(K,I)))/(GAMMA(K,I)*H2(K,I))        
        !
        if (K == KB(I)) then
          BIBH2(K,I) = BI(K,I) / BH2(K,I)
        else
          BIBH2(K,I) = (BI(K,I)-BI(K+1,I)) / BH2(K,I)
        end if
        
        !JPOC = (LPOMCD(K,I)+LRPOMCD(K,I) + RPOMCD(K,I)+RPOMCHD(K,I))/DAY    ! gC/m3/s
        JPOC(K,I) = (LPOMCD(K,I)+LRPOMCD(K,I) + RPOMCD(K,I)+RPOMCHD(K,I) + (LDOMCD(K,I)+RDOMCD(K,I)))   !/DAY    ! add DOC flux
        JPOC(K,I) = JPOC(K,I)*ANOX(JW)*KDO(JW)/(KDO(JW)+O2(K,I))/DO4(K,I)   ! ANOXIC FRACTION OF TOTAL C TURNOVER RATE     ! gC/m3/s
        
        if (use_BedHg) then
          if (SEDIMENT_CALC(JW)) then
            ! First-order sediment compartment
            JPOC2(K,I) = SEDDC(K,I) / BedH(r)           ! gC/m3/s         
          else if (SED_DIAG =='      ON') then
            !JPOC2(K,I) = KFSF(K,I,1)/DAY / BedH(r)      ! gC/m3/s
            JPOC2(K,I) =  ((SD_kdiaPOC(1) * SD_ThtaPOC(1) ** (C2SF(K,I,16) - 20.))*C2SF(K,I,19) +  (SD_kdiaPOC(2) * SD_ThtaPOC(2) ** (C2SF(K,I,16) - 20.))*C2SF(K,I,20))/DAY      ! gC/m3/s
          else
            JPOC2(K,I) = 0.0
          end if
        end if       

        ! Hg0 volatilization
        ! A user-defined volatilization velocity or O2 based rearation rate is used.
        if (.NOT. ICE(I).AND.K == KT) then
          KH_Hg = 10.0**(-1078.0/(T1(KT,I)+273.15) - LOG10(T1(KT,I)+273.15) + 5.592)   ! Sanemasa(1975) equation
          if (vv_20(r)>0.0) then
            Hg0_Vol(KT,I) = vv_20(r)*(vv_theta(r)**(T1(KT,I)-20.0))/DAY * BI(KT,I)/BH2(KT,I) * (Hg0(KT,I)-Hg00(r)/KH_Hg)
          else
            Hg0_Vol(KT,I) = 0.632*REAER(I)  * (Hg0(KT,I)-Hg00(r)/KH_Hg)   ! REAER in units of 1/s  * BI(KT,I)/BH2(KT,I)
          end if
        else
          Hg0_Vol(KT,I) = 0.0  
        end if
        !
        ! HgII reduction 
        ! m^2/W/d --> m^2/W/s 
        HgII_Reduction(K,I)  =  Hgd(2,K,I) * k_reduct(r)/DAY*LIGHT
        if (use_DOCSorbed(2)) HgII_Reduction(K,I) = HgII_Reduction(K,I) + Hgdoc(2,K,I) * kdoc_reduct(r)/DAY*LIGHT
        !
        ! Hg0 oxidation
        IF (PH_CALC(JW)) THEN
          DGMratiot = k_oxid1(r) * DOC(K,I)**k_oxid2(r) * PH(K,I)**k_oxid3(r)
        ELSE  
          ! if PH is not calculated --> PH = 7.0
          DGMratiot = k_oxid1(r) * DOC(K,I)**k_oxid2(r) * 7.0**k_oxid3(r)
        END IF
        if ((Hgd(2,K,I)+Hgdoc(2,K,I)) > 0.0) then
          DGMratio  = Hg0(K,I) / (Hgd(2,K,I)+Hgdoc(2,K,I))
        else
          DGMratio  = 0.0
        end if
        if (DGMratiot > 0.0) then
          Hg0_Oxidation(K,I) = DGMratio/DGMratiot * HgII_Reduction(K,I)
        else
          Hg0_Oxidation(K,I) = 0.0
        end if
        ! 
        ! HgII settling and burial 
        if (use_POMSorbed(2)) then
          if (K == KT) then  
            HgII_Settling(K,I) = -POMS(JW) * Hgpom(2,K,I) * BI(K,I)/BH2(K,I)     ! ng/l/s
          else 
            HgII_Settling(K,I) = POMS(JW) * (Hgpom(2,K-1,I)-Hgpom(2,K,I)) * BI(K,I)/BH2(K,I)    !*BI(K,I)/BI(K-1,I)
          end if
          if (use_BedHg) then
            HgII_Bed(K,I) = Hgpom(2,K,I)  * POMS(JW)/BedH(r)  
            !HgII_Burial = Hgpom2(2,K,I) * vb(r)/DAY/BedH(r)  
            !HgII_Burial(K,I) = Hgpom2(2,K,I) * vbs(K,I)/DAY/BedH(r) ! SW code fix
            HgII_Burial(K,I) = Hgpom2(2,K,I) * vbs(K,I)/BedH(r) 

          end if    
        else
          HgII_Settling(K,I) = 0.0
          if (use_BedHg) then
            HgII_Bed(K,I) = 0.0
            HgII_Burial(K,I) = 0.0  
          end if
        end if
        !
        if (use_AnyAlgaeSorbed(2)) then
          do JA = 1,NAL
            if (ALG_CALC(JA)) then
              if (AS(JA) >= 0.0) then
                if (K == KT) then
                  HgII_Settling(K,I) = HgII_Settling(K,I) - AS(JA)*Hgap(2,JA,K,I) * BI(K,I)/BH2(K,I)  
                else
                  HgII_Settling(K,I) = HgII_Settling(K,I) + (AS(JA)*Hgap(2,JA,K-1,I)-AS(JA)*Hgap(2,JA,K,I)) * BI(K,I)/BH2(K,I)    !*BI(K,I)/BI(K-1,I)
                end if
                if (use_BedHg) then 
                  HgII_Bed(K,I) = HgII_Bed(K,I) + AS(JA)*Hgap(2,JA,K,I) / BedH(r) 
                end if
              else
                if (K == KB(I)) then
                  HgII_Settling(K,I) = HgII_Settling(K,I) + AS(JA)*Hgap(2,JA,K,I) * BI(K,I)/BH2(K,I)   
                else if (K == KT) then
                  HgII_Settling(K,I) = HgII_Settling(K,I) - AS(JA)*Hgap(2,JA,K+1,I) * BI(K+1,I)/BH2(K,I)
                else
                  HgII_Settling(K,I) = HgII_Settling(K,I) - AS(JA)*(Hgap(2,JA,K+1,I) * BI(K+1,I)/BH2(K,I) - Hgap(2,JA,K,I) * BI(K,I)/BH2(K,I)) 
                end if               
              end if          
            end if
          end do  
        end if
        ! Zooplankton
        if (use_AnyZooSorbed(2)) then
          IF (ZOOPLANKTON_CALC) THEN
            do JZ = 1,NZP
              if (ZS(JZ) >= 0.0) then
                if (K == KT) then
                  HgII_Settling(K,I) = HgII_Settling(K,I) - ZS(JZ)*Hgzp(2,JZ,K,I) * BI(K,I)/BH2(K,I)  
                else
                  HgII_Settling(K,I) = HgII_Settling(K,I) + ZS(JZ)*(Hgzp(2,JZ,K-1,I)-Hgzp(2,JZ,K,I)) * BI(K,I)/BH2(K,I)     !*BI(K,I)/BI(K-1,I)
                end if
                if (use_BedHg) then 
                  HgII_Bed(K,I) = HgII_Bed(K,I) + ZS(JZ)*Hgzp(2,JZ,K,I) / BedH(r)
                end if
              else
                if (K == KB(I)) then
                  HgII_Settling(K,I) = HgII_Settling(K,I) + ZS(JZ)*Hgzp(2,JZ,K,I) * BI(K,I)/BH2(K,I)
                else if (K == KT) then
                  HgII_Settling(K,I) = HgII_Settling(K,I) - ZS(JZ)*Hgzp(2,JZ,K+1,I) * BI(K+1,I)/BH2(K,I)  
                else
                  HgII_Settling(K,I) = HgII_Settling(K,I) - ZS(JZ)*(Hgzp(2,JZ,K+1,I) * BI(K+1,I)/BH2(K,I) - Hgzp(2,JZ,K,I) * BI(K,I)/BH2(K,I)) 
                end if               
              end if          
            end do  
          END IF
        end if        
        !
        if (use_AnySolidSorbed(2)) then
          do JS = 1,NSS
            if (K == KT) then   
              HgII_Settling(K,I) = HgII_Settling(K,I) - SSS(JS)*Hgp(2,JS,K,I) * BI(K,I)/BH2(K,I)  
            else
              HgII_Settling(K,I) = HgII_Settling(K,I) + SSS(JS)*(Hgp(2,JS,K-1,I)-Hgp(2,JS,K,I)) * BI(K,I)/BH2(K,I)     !*BI(K,I)/BI(K-1,I)
            end if
            if (use_BedHg) then 
              HgII_Bed(K,I) = HgII_Bed(K,I) + Hgp(2,JS,K,I) * SSS(JS)/BedH(r)
              !HgII_Burial = HgII_Burial + Hgp2(2,JS,K,I) * vb(r)/DAY/BedH(r) 
              !HgII_Burial(K,I) = HgII_Burial(K,I) + Hgp2(2,JS,K,I) * vbs(K,I)/DAY/BedH(r)  ! SW fix
              HgII_Burial(K,I) = HgII_Burial(K,I) + Hgp2(2,JS,K,I) * vbs(K,I)/BedH(r)   
            end if
          end do
        end if
        !
        ! HgII methylation under anoxic conditions
        ! m^3/gC * gC/m^3/s --> 1/s
        !xx = KDOHg(r)/(O2(K,I)+KDOHg(r)) 
        !xx = KDOHg1(r)/(KDOHg1(r)+(KDOn(r)*O2(K,I))**KDOHg2(r))    ! D. Hutchinson    
        !HgII_Methylation(K,I) = Hgd(2,K,I) * xx*k_meth(r)*JPOC(K,I)          ! ng/l/s
        HgII_Methylation(K,I) = Hgd(2,K,I) * k_meth(r)*JPOC(K,I)          ! ng/l/s
        !if(use_DOCSorbed(2)) HgII_Methylation(K,I) = HgII_Methylation(K,I) + Hgdoc(2,K,I) * xx*kdoc_meth(r)*JPOC(K,I) 
        if(use_DOCSorbed(2)) HgII_Methylation(K,I) = HgII_Methylation(K,I) + Hgdoc(2,K,I) * kdoc_meth(r)*JPOC(K,I) 

        if (use_BedHg) then
          HgII2_Methylation(K,I) = Cd2(2,K,I) * k_meth2(r) * JPOC2(K,I) 
          if(use_DOCSorbed(2)) HgII2_Methylation(K,I) = HgII2_Methylation(K,I) + Cdoc2(2,K,I) * kdoc_meth2(r)*JPOC2(K,I) 
        end if
        !
        ! MeHg settling and burial 
        if (use_POMSorbed(3)) then
          if (K == KT) then  
            MeHg_Settling(K,I) = -POMS(JW) * Hgpom(3,K,I) * BI(K,I)/BH2(K,I)
          else 
            MeHg_Settling(K,I) = POMS(JW) * (Hgpom(3,K-1,I)-Hgpom(3,K,I)) * BI(K,I)/BH2(K,I)      !*BI(K,I)/BI(K-1,I)
          end if
          if (use_BedHg) then
            MeHg_Bed(K,I) = Hgpom(3,K,I)  * POMS(JW)/BedH(r)
            !MeHg_Burial = Hgpom2(3,K,I) * vb(r)/DAY/BedH(r)
            !MeHg_Burial(K,I) = Hgpom2(3,K,I) * vbs(K,I)/DAY/BedH(r)   !SW fix
            MeHg_Burial(K,I) = Hgpom2(3,K,I) * vbs(K,I)/BedH(r)
          end if    
          !
        else
          MeHg_Settling(K,I) = 0.0
          if (use_BedHg) then
            MeHg_Bed(K,I) = 0.0
            MeHg_Burial(K,I) = 0.0
          end if
        end if  
        !
        if (use_AnyAlgaeSorbed(3)) then
          do JA = 1,NAL
            if (ALG_CALC(JA)) then
              if (AS(JA) >= 0.0) then
                if (K == KT) then
                  MeHg_Settling(K,I) = MeHg_Settling(K,I) - AS(JA)*Hgap(3,JA,K,I) * BI(K,I)/BH2(K,I)
                else
                  MeHg_Settling(K,I) = MeHg_Settling(K,I) + AS(JA)*(Hgap(3,JA,K-1,I)-Hgap(3,JA,K,I)) * BI(K,I)/BH2(K,I)      !*BI(K,I)/BI(K-1,I)
                end if
                if (use_BedHg) then 
                  MeHg_Bed(K,I) = MeHg_Bed(K,I) + AS(JA)*Hgap(3,JA,K,I) / BedH(r)
                end if
              else
                if (K == KB(I)) then
                  MeHg_Settling(K,I) = MeHg_Settling(K,I) + AS(JA)*Hgap(3,JA,K,I) * BI(K,I)/BH2(K,I)
                else if (K == KT) then
                  MeHg_Settling(K,I) = MeHg_Settling(K,I) - AS(JA)*Hgap(3,JA,K+1,I) * BI(K+1,I)/BH2(K,I)   
                else
                  MeHg_Settling(K,I) = MeHg_Settling(K,I) - AS(JA)*(Hgap(3,JA,K+1,I) * BI(K+1,I)/BH2(K,I) - Hgap(3,JA,K,I) * BI(K,I)/BH2(K,I)) 
                end if               
              end if  
            end if
          end do 
        end if 
        ! Zooplankton
        if (use_AnyZooSorbed(3)) then
          IF (ZOOPLANKTON_CALC) THEN
            do JZ = 1,NZP
              if (ZS(JZ) >= 0.0) then
                if (K == KT) then
                  MeHg_Settling(K,I) = MeHg_Settling(K,I) - ZS(JZ)*Hgzp(3,JZ,K,I) * BI(K,I)/BH2(K,I)
                else
                  MeHg_Settling(K,I) = MeHg_Settling(K,I) + ZS(JZ)*(Hgzp(3,JZ,K-1,I)-Hgzp(3,JZ,K,I)) * BI(K,I)/BH2(K,I)     !*BI(K,I)/BI(K-1,I)
                end if
                if (use_BedHg) then 
                  MeHg_Bed = MeHg_Bed + ZS(JZ)*Hgzp(3,JZ,K,I) / BedH(r)
                end if
              else
                if (K == KB(I)) then
                  MeHg_Settling(K,I) = MeHg_Settling(K,I) + ZS(JZ)*Hgzp(3,JZ,K,I) * BI(K,I)/BH2(K,I)
                else if (K == KT) then
                  MeHg_Settling(K,I) = MeHg_Settling(K,I) - ZS(JZ)*Hgzp(3,JZ,K+1,I) * BI(K+1,I)/BH2(K,I)   
                else
                  MeHg_Settling(K,I) = MeHg_Settling(K,I) - ZS(JZ)*(Hgzp(3,JZ,K+1,I) * BI(K+1,I)/BH2(K,I) - Hgzp(3,JZ,K,I) * BI(K,I)/BH2(K,I)) 
                end if               
              end if  
            end do 
          END IF
        end if         
        !
        if (use_AnySolidSorbed(3)) then
          do JS = 1,NSS
            if (K == KT) then   
              MeHg_Settling(K,I) = MeHg_Settling(K,I) - SSS(JS)*Hgp(3,JS,K,I) * BI(K,I)/BH2(K,I)  
            else
              MeHg_Settling(K,I) = MeHg_Settling(K,I) + SSS(JS)*(Hgp(3,JS,K-1,I)-Hgp(3,JS,K,I)) * BI(K,I)/BH2(K,I)      !*BI(K,I)/BI(K-1,I)
            end if
            if (use_BedHg) then 
              MeHg_Bed(K,I) = MeHg_Bed(K,I) + Hgp(3,JS,K,I) * SSS(JS)/BedH(r)
              !MeHg_Burial = MeHg_Burial + Hgp2(3,JS,K,I) * vb(r)/DAY/BedH(r)
              !MeHg_Burial(K,I) = MeHg_Burial(K,I) + Hgp2(3,JS,K,I) * vbs(K,I)/DAY/BedH(r)  ! SW fix
              MeHg_Burial(K,I) = MeHg_Burial(K,I) + Hgp2(3,JS,K,I) * vbs(K,I)/BedH(r)
            end if  
          end do  
        end if
        !
        ! MeHg demethylation under anoxic conditions
        !MeHg_Demethylation(K,I) = Hgd(3,K,I) * xx*k_demeth(r)*JPOC(K,I)
         MeHg_Demethylation(K,I) = Hgd(3,K,I) * k_demeth(r)*JPOC(K,I)
        !if(use_DOCSorbed(3)) MeHg_Demethylation(K,I) = MeHg_Demethylation(K,I) + Hgdoc(3,K,I) * xx*kdoc_demeth(r)*JPOC(K,I)
        if(use_DOCSorbed(3)) MeHg_Demethylation(K,I) = MeHg_Demethylation(K,I) + Hgdoc(3,K,I) * kdoc_demeth(r)*JPOC(K,I)
        if (use_BedHg) then
          MeHg2_Demethylation(K,I) = Cd2(3,K,I) * k_demeth2(r) * JPOC2(K,I) 
          if(use_DOCSorbed(2)) MeHg2_Demethylation(K,I) = MeHg2_Demethylation(K,I) + Cdoc2(3,K,I) * kdoc_demeth2(r) *JPOC2(K,I)
        end if        
        !
        ! MeHg photodegradation 
        ! m^2/W/d --> m^2/W/s 
        MeHg_Decay(K,I) = Hgd(3,K,I) * k_decay(r)/DAY*LIGHT
        if(use_DOCSorbed(3)) MeHg_Decay(K,I) = MeHg_Decay(K,I) + Hgdoc(3,K,I) * kdoc_decay(r)/DAY*LIGHT
        !
        ! sediment mass transfer 
        if(vm1(2,r) < 0.0)then
            vm(2,r)=c2sf(k,i,38)*abs(vm1(2,r))/DAY      ! This is the flux rate from sed diag model which must be ON
            
        else
            vm(2,r)=vm1(2,r)/100.0/DAY    ! cm/d to m/s
        endif

        if(vm1(3,r) < 0.0)then
            vm(3,r)=c2sf(k,i,38)*abs(vm1(3,r))/DAY
            
        else
            vm(3,r)=vm1(3,r)/100.0/DAY    ! cm/d to m/s
        endif      
        
        if (use_BedHg) then
          HgII_Transfer(K,I)  =  vm(2,r) * (Cd2(2,K,I)/BedPor(r) - Hgd(2,K,I)) * BIBH2(K,I)   !/100.0/DAY
          HgII2_Transfer(K,I) = -vm(2,r) * (Cd2(2,K,I)/BedPor(r) - Hgd(2,K,I)) / BedH(r)      !/100.0/DAY     
          if (use_DOCSorbed(2)) then
            HgII_Transfer(K,I)  = HgII_transfer(K,I)  + vm(2,r) * (Cdoc2(2,K,I)/BedPor(r) - Hgdoc(2,K,I)) * BIBH2(K,I)   ! /100.0/DAY
            HgII2_Transfer(K,I) = HgII2_transfer(K,I) - vm(2,r) * (Cdoc2(2,K,I)/BedPor(r) - Hgdoc(2,K,I)) / BedH(r)      !/100.0/DAY  
          end if
          !
          MeHg_Transfer(K,I)  =  vm(3,r)* (Cd2(3,K,I)/BedPor(r) - Hgd(3,K,I)) * BIBH2(K,I)   ! /100.0/DAY 
          MeHg2_Transfer(K,I) = -vm(3,r) * (Cd2(3,K,I)/BedPor(r) - Hgd(3,K,I)) / BedH(r)       !/100.0/DAY    
          if (use_DOCSorbed(3)) then
            MeHg_Transfer(K,I)  = MeHg_Transfer(K,I)  + vm(3,r) * (Cdoc2(3,K,I)/BedPor(r) - Hgdoc(3,K,I)) * BIBH2(K,I)    ! /100.0/DAY
            MeHg2_Transfer(K,I) = MeHg2_Transfer(K,I) - vm(3,r) * (Cdoc2(3,K,I)/BedPor(r) - Hgdoc(3,K,I)) / BedH(r)       ! /100.0/DAY
          end if
        else
          ! 
          ! vm ([ng/L/d]/[mg-O2/L/s] SODD==gO2/m2/s]) is specified with the sediment release rate of Hg as a fraction of the SOD if the sediment layer is not modled.
          HgII_Transfer(K,I)  =  vm(2,r) * SODD(K,I)*DO2(K,I)       !/100./DAY
          MeHg_Transfer(K,I)  =  vm(3,r) * SODD(K,I)*DO2(K,I)       !/100./DAY
        end if
        !
      END DO
    END DO     
  end subroutine
  !===========================================================================================================================
  ! compute kinetic source/sink terms of Hg0, HgII, MeHg in the water column
  subroutine HgEquilibriumKinetics
    DO I=IU,ID
      DO K=KT,KB(I)
        Hg0SS(K,I) = -Hg0_Vol(K,I) - Hg0_Oxidation(K,I) + HgII_Reduction(K,I) + MeHg_Decay(K,I)        
        if (use_BedHg) then
          MeHgSS(K,I) = MeHg_Settling(K,I) - MeHg_Decay(K,I) + HgII_Methylation(K,I) - MeHg_Demethylation(K,I)
          HgIISS(K,I) = HgII_Settling(K,I) + Hg0_Oxidation(K,I) - HgII_Reduction(K,I) - HgII_Methylation(K,I) + MeHg_Demethylation(K,I)
        else
          MeHgSS(K,I) = MeHg_Settling(K,I) - MeHg_Decay(K,I) + HgII_Methylation(K,I) - MeHg_Demethylation(K,I) + MeHg_transfer(K,I)
          HgIISS(K,I) = HgII_Settling(K,I) + Hg0_Oxidation(K,I) - HgII_Reduction(K,I) - HgII_Methylation(K,I) + MeHg_Demethylation(K,I) + HgII_Transfer(K,I) 
        end if

      END DO
    END DO
  end subroutine
  !===========================================================================================================================
  ! compute kinetics of HgII and MeHg in the sediment layer
  subroutine ComputeBedHg
    DO I=IU,ID
      DO K=KT,KB(I)
        dHg2dt(2)   = HgII_Bed(K,I) + HgII2_Transfer(K,I) - HgII_Burial(K,I) - HgII2_Methylation(K,I) + MeHg2_Demethylation(K,I)
        Hg2(2,K,I)  = max((Hg2(2,K,I) + dHg2dt(2)*DLT), 0.0)
        HgIISS(K,I) = HGIISS(K,I) + HgII_Transfer(K,I) 
        !
        dHg2dt(3)   = MeHg_Bed(K,I) + MeHg2_Transfer(K,I) - MeHg_Burial(K,I) + HgII2_Methylation(K,I) - MeHg2_Demethylation(K,I)
        Hg2(3,K,I)  = max((Hg2(3,K,I) + dHg2dt(3)*DLT), 0.0)
        MeHgSS(K,I) = MEHGSS(K,I) + MeHg_transfer(K,I)
      END DO
    END DO
  end subroutine
  !===========================================================================================================================
  ! compute derived variables and pathway fluxes
  subroutine  HgOutputs
  USE GLOBAL, ONLY: CD; USE MAIN, ONLY:PHG2L_DER,PHG2S_DER,PMHGL_DER,PMHGS_DER,DHG2_DER,DMHG_DER
    implicit none
    !
    real(R8) :: W_TSS, Bed_TSS 
    
    DO I=IU,ID
      r = SegRegion(I)
      DO K=KT,KB(I)
        !
        do i23 = 2, 3
          !
          ! water column, add zooplankton
          ! W_TSS = POM + ISS + ALG (mg/L) + ZOO (mg/L)
          W_TSS = 0.0
          Hgpt(i23,K,I)  = 0.0 
          Hgpts(i23,K,I) = 0.0 
          if (use_POMSorbed(i23)) then         
            W_TSS = W_TSS + (LPOM(K,I)+RPOM(K,I))
            Hgpt(i23,K,I) = Hgpt(i23,K,I) + Hgpom(i23,K,I)
          end if    
          if (use_AnySolidSorbed(i23)) then
            do JS = 1,NSS  
              W_TSS = W_TSS + SS(K,I,JS)
              Hgpt(i23,K,I) = Hgpt(i23,K,I) + Hgp(i23,JS,K,I)
            end do 
          end if
          if (use_AnyAlgaeSorbed(i23)) then
            do JA = 1,NAL
              if (ALG_CALC(JA)) then
                W_TSS = W_TSS + ALG(K,I,JA)
                Hgpt(i23,K,I) = Hgpt(i23,K,I) + Hgap(i23,JA,K,I)
              end if
            end do 
          end if       
          ! Zooplankton           
          if (use_AnyZooSorbed(i23)) then
            IF (ZOOPLANKTON_CALC) THEN
              do JZ = 1,NZP
                W_TSS = W_TSS + ZOO(K,I,JZ)
                Hgpt(i23,K,I) = Hgpt(i23,K,I) + Hgzp(i23,JZ,K,I)
              end do 
            END IF
          end if          
          !
          if (W_TSS > 0.0) then
            Hgpts(i23,K,I) = Hgpt(i23,K,I)/W_TSS*1.0E3   ! ng/g
          else
            Hgpts(i23,K,I) = 0.0
          end if 
          !
          ! sediment layer
          if (use_BedHg) then 
            !
            ! pore water concentrations
            Hgd2(2,K,I)   = Cd2(2,K,I)   / BedPor(r)
            Hgdoc2(2,K,I) = Cdoc2(2,K,I) / BedPor(r)
            Hgd2(3,K,I)   = Cd2(3,K,I)   / BedPor(r)
            Hgdoc2(3,K,I) = Cdoc2(3,K,I) / BedPor(r)  
            !
            ! Bed_TSS = POM + ISS (mg/L)
            Bed_TSS         = 0.0
            Hgpt2(i23,K,I)  = 0.0 
            Hgpts2(i23,K,I) = 0.0  
            if (use_POMSorbed(i23)) then       
              Bed_TSS = Bed_TSS + Bed_POM
              Hgpt2(i23,K,I) = Hgpt2(i23,K,I) + Hgpom2(i23,K,I) 
            end if
            if (use_AnySolidSorbed(i23)) then
              do JS = 1,NSS
                Bed_TSS = Bed_TSS + BedISS(K,I,JS)
                Hgpt2(i23,K,I) = Hgpt2(i23,K,I) + Hgp2(i23,JS,K,I)
              end do
            end if
            if (Bed_TSS > 0.0) then 
              Hgpts2(i23,K,I) = Hgpt2(i23,K,I)/Bed_TSS*1.0E3  ! ng/g
            else
              Hgpts2(i23,K,I) = 0.0
            end if
          end if
        end do
        !
        ! HgII
        CON_HG(1,K,I) = Hgd(2,K,I)
        CON_HG(2,K,I) = Hgdoc(2,K,I)
        CON_HG(3,K,I) = Hgpom(2,K,I)    
        do JS = 1,NSS
          CON_HG(3+JS,K,I) = Hgp(2,JS,K,I)
        end do
        do JA = 1,NAL
          CON_HG(NSS+3+JA,K,I) = Hgap(2,JA,K,I)
        end do 
        ! Zooplankton
        do JZ = 1,NZP
          CON_HG(NSS+NAL+3+JZ,K,I) = Hgzp(2,JZ,K,I)
        end do 
        !
        CON_HG(NSS+NAL+NZP+4,K,I) = Hgpt(2,K,I)
        CON_HG(NSS+NAL+NZP+5,K,I) = Hgpts(2,K,I)
        !        
        ! MeHg
        CON_HG(NSS+NAL+NZP+6,K,I) = Hgd(3,K,I)
        CON_HG(NSS+NAL+NZP+7,K,I) = Hgdoc(3,K,I)  
        CON_HG(NSS+NAL+NZP+8,K,I) = Hgpom(3,K,I) 
        do JS = 1,NSS
          CON_HG(NSS+NAL+NZP+8+JS,K,I) = Hgp(3,JS,K,I)
        end do
        do JA = 1,NAL
          CON_HG(2*NSS+NAL+NZP+8+JA,K,I) = Hgap(3,JA,K,I)
        end do    
        ! Zooplankton
        do JZ = 1,NZP
          CON_HG(2*(NSS+NAL)+NZP+8+JZ,K,I) = Hgzp(3,JZ,K,I)
        end do  
        !
        CON_HG(2*(NSS+NAL+NZP)+9,K,I) = Hgpt(3,K,I)
        CON_HG(2*(NSS+NAL+NZP)+10,K,I) = Hgpts(3,K,I)       
        !
        ! total concentration of Hg in water (ng/L)
        CON_HG(2*(NSS+NAL+NZP)+11,K,I) = Hg0(K,I) + HgII(K,I) + MeHg(K,I) 
        !
        ! pathway fluxes (ng/L/s)
        KF_HG(1,K,I) = Hg0_Vol(K,I) 
        KF_HG(2,K,I) = Hg0_Oxidation(K,I)
        KF_HG(3,K,I) = -HgII_Reduction(K,I)
        KF_HG(4,K,I) = -HgII_Methylation(K,I)
        KF_HG(5,K,I) = HgII_Settling(K,I)
        KF_HG(6,K,I) = HgII_Transfer(K,I)
        KF_HG(7,K,I) = -MeHg_Decay(K,I)
        KF_HG(8,K,I) = -MeHg_Demethylation(K,I)
        KF_HG(9,K,I) = MeHg_Settling(K,I) 
        KF_HG(10,K,I) = MeHg_Transfer(K,I)
        
        !
        if (use_BedHg) then 
          !
          ! HgII
          CON_HG2(1,K,I) = Hg2(2,K,I)
          CON_HG2(2,K,I) = Hgd2(2,K,I)
          CON_HG2(3,K,I) = Hgdoc2(2,K,I)
          CON_HG2(4,K,I) = Hgpom2(2,K,I)
          do JS = 1,NSS
            CON_HG2(4+JS,K,I) = Hgp2(2,JS,K,I)
          end do
          CON_HG2(NSS+5,K,I) = Hgpt2(2,K,I)
          CON_HG2(NSS+6,K,I) = Hgpts2(2,K,I)
          !
          ! MeHg
          CON_HG2(NSS+7,K,I) = Hg2(3,K,I)
          CON_HG2(NSS+8,K,I) = Hgd2(3,K,I)
          CON_HG2(NSS+9,K,I) = Hgdoc2(3,K,I)
          CON_HG2(NSS+10,K,I) = Hgpom2(3,K,I)
          do JS = 1,NSS
            CON_HG2(NSS+10+JS,K,I) = Hgp2(3,JS,K,I)
          end do
          CON_HG2(2*NSS+11,K,I) = Hgpt2(3,K,I)
          CON_HG2(2*NSS+12,K,I) = Hgpts2(3,K,I)
          !
          ! total concentration of Hg in sediment (ng/L)
          CON_HG2(2*NSS+13,K,I) = Hg2(2,K,I) + Hg2(3,K,I)  
          !
          ! ISS
          do JS = 1,NSS
            CON_HG2(2*NSS+13+JS,K,I) = BedISS(K,I,JS)
          end do
          !
          ! pathway fluxes (ng/L/s)
          KF_HG2(1,K,I) = -HgII2_Methylation(K,I)
          KF_HG2(2,K,I) = -HgII_Burial(K,I)
          KF_HG2(3,K,I) =  HgII_Bed(K,I)
          KF_HG2(4,K,I) = -MeHg2_Demethylation(K,I)   
          KF_HG2(5,K,I) = -MeHg_Burial(K,I)
          KF_HG2(6,K,I) =  MeHg_Bed(K,I)

        end if 
      
        !
        ! Set derived output
        CD(K,I,PHG2L_DER)=HGPT(2,K,I)
        CD(K,I,PHG2S_DER)=HGPTS(2,K,I)
        CD(K,I,PMHGL_DER)=HGPT(3,K,I)
        CD(K,I,PMHGS_DER)=HGPTS(3,K,I)
        CD(K,I,DHG2_DER)=HGD(2,K,I)+HGDOC(2,K,I)
        CD(K,I,DMHG_DER)=HGD(3,K,I)+HGDOC(3,K,I)
        !PHG2L_DER=28;PHG2S_DER=29;PMHGL_DER=30;PMHGS_DER=31;DHG2_DER=32;DMHG_DER=33
        
      END DO
    END DO
  end subroutine 
  
  subroutine HgKineticFluxes
  
  DO JW=1,NWB
    KT=KTWB(JW)
    DO JB=BS(JW),BE(JW)   
        
            IF (TRIBUTARIES) THEN
              DO JT=1,JTT
                IF (JB == JBTR(JT)) THEN
                  I = ITR(JT)
                  IF (I < CUS(JB)) I = CUS(JB)
                  DO K=KTTR(JT),KBTR(JT)
                    IF (QTR(JT) < 0.0) THEN
                      THgOUT = THgOUT-(HG0(K,I)+HGII(K,I)+MEHG(K,I))*QTR(JT)*QTRF(K,JT)*DLT*1000.   ! mass in ng: ng/l X m3/s X s X 1000 l/m3
                      MeHgOUT= MeHgOUT-MEHG(K,I)*QTR(JT)*QTRF(K,JT)*DLT*1000.
                    ELSE   
                      THgTRIB = THgTRIB+ (CTR(NHg0,JT)+CTR(NHgII,JT)+CTR(NMeHg,JT))*QTR(JT)*QTRF(K,JT)*DLT*1000.
                      MeHgTRIB = MeHgTRIB+ CTR(NMeHg,JT)*QTR(JT)*QTRF(K,JT)*DLT*1000.
                    END IF
                  END DO
                END IF
              END DO
            END IF
            IF (DIST_TRIBS(JB)) THEN
              DO I=CUS(JB),DS(JB)
                IF (QDT(I) < 0.0) THEN
                      THgOUT = THgOUT-(HG0(KT,I)+HGII(KT,I)+MEHG(KT,I))*QDT(I)*DLT*1000.
                      MeHgOUT = MeHgOUT-MEHG(KT,I)*QDT(I)*DLT*1000.
                ELSE
                      THgDTRIB = THgDTRIB+ (CDTR(NHg0,JB)+CDTR(NHgII,JB)+CDTR(NMeHg,JB))*QDT(I)*DLT*1000.
                      MeHgDTRIB = MeHgDTRIB+ CDTR(NMeHg,JB)*QDT(I)*DLT*1000.
                END IF
              END DO
            END IF
            IF (WITHDRAWALS) THEN
              DO JWD=1,JWW
                IF (QWD(JWD) /= 0.0) THEN
                  IF (JB == JBWD(JWD)) THEN
                    I = MAX(CUS(JBWD(JWD)),IWD(JWD))
                    DO K=KTW(JWD),KBW(JWD)      
                      THgWD = THgWD + (C1S(K,I,NHg0)+C1S(K,I,NHgII)+C1S(K,I,NMeHg))*QSW(K,JWD)*DLT*1000.
                      MeHgWD = MeHgWD + C1S(K,I,NMeHg)*QSW(K,JWD)*DLT*1000.
                    END DO
                  END IF
                END IF
              END DO
            END IF
            IF (PRECIPITATION(JW)) THEN
              DO I=CUS(JB),DS(JB)
                      THgPR = THgPR + (CPR(NHg0,JB)+CPR(NHgII,JB)+CPR(NMeHg,JB))*QPR(I)*DLT*1000.
                      MeHgPR = MeHgPR + CPR(NMeHg,JB)*QPR(I)*DLT*1000.
              END DO
            END IF
            IF (UP_FLOW(JB) .and. .not. head_flow(jb) ) THEN
                      THgIN = THgIN + (CIND(NHg0,JB)+CIND(NHGII,JB)+CIND(NMeHg,JB))*QIND(JB)*DLT*1000.
                      MeHgIN = MeHgIN + CIND(NMeHg,JB)*QIND(JB)*DLT*1000.
            END IF
            IF (DN_FLOW(JB)) THEN
                DO K=KT,KB(DS(JB))
                      THgOUT = THgOUT + (C1S(K,DS(JB),NHg0)+C1S(K,DS(JB),NHgII)+C1S(K,DS(JB),NMeHg))*QOUT(K,JB)*DLT*1000.   ! C1S(KT:KB(ID),ID,JC)
                      MeHgOUT = MeHgOUT + C1S(K,DS(JB),NMeHg)*QOUT(K,JB)*DLT*1000. 
                ENDDO
            ENDIF      
        
    DO I=CUS(JB),DS(JB)
      r = SegRegion(I) 
      DO K=KT,KB(I)
          DO N1=1,SREGION_INI
            IF(I>= SREGION(N1) .AND. I <= EREGION(N1))THEN
                IF(K==KT)THEN
                ATMDEP_HG0(N1)=ATMDEP_HG0(N1)+ATM_DEP_LOADING(NHg0,JW)*BI(KT,I)*DLX(I)*3.17098E-8*DLT   ! Mass in ng  ! Conversion: mg/km2/year to ug/s 1000/(365*86400*1000*1000)=3.17098E-11 X 10^3 =3.17098E-8 to convert to ng
                ATMDEP_HG2(N1)=ATMDEP_HG2(N1)+ATM_DEP_LOADING(NHgII,JW)*BI(KT,I)*DLX(I)*3.17098E-8*DLT 
                ATMDEP_MEHG(N1)=ATMDEP_MEHG(N1)+ATM_DEP_LOADING(NMeHg,JW)*BI(KT,I)*DLX(I)*3.17098E-8*DLT 
                ENDIF
                DO J1=1,10
                    KF_HG_SUM(J1,N1)=KF_HG_SUM(J1,N1)+KF_HG(J1,K,I)*VOL(K,I)*1000.*DLT  ! Mass in ng: ng/l/s*m3*1000l/m3*s
                ENDDO
                if (use_BedHg) then
                if(k==kb(i))then
                        DO J1=1,6
                            KF_HG2_SUM(J1,N1)=KF_HG2_SUM(J1,N1)+KF_HG2(J1,K,I)*BedH(r)*dlx(i)*b(kb(i),i)*1000.*DLT  ! Mass in ng; use volume of bed not water   !VOL(K,I)
                        ENDDO
                else
                        DO J1=1,6
                            KF_HG2_SUM(J1,N1)=KF_HG2_SUM(J1,N1)+KF_HG2(J1,K,I)*BedH(r)*dlx(i)*(b(k,i)-b(k+1,i))*1000.*DLT  ! Mass in ng; use volume of bed not water   !VOL(K,I)
                        ENDDO
                endif
                
                endif
            ENDIF
          ENDDO
      ENDDO
    ENDDO
    ENDDO
  ENDDO
  
  
         
      IF (JDAY >= NXTMFLHG ) THEN  
          NXTMFLHG = NXTMFLHG + HGFLUXINT
          
          DO N1=1,SREGION_INI
          if (use_BedHg) then
          WRITE(FLUX_OUTPUT_FILE,'(F12.3,",",I2,",",10(E12.5,","),25(E12.5,","))')JDAY,N1,(KF_HG_SUM(J1,N1)/(JDAY-JDAYSTARTHG),J1=1,10),(KF_HG2_SUM(J1,N1)/(JDAY-JDAYSTARTHG),J1=1,6),ATMDEP_HG0(N1)/(JDAY-JDAYSTARTHG),ATMDEP_HG2(N1)/(JDAY-JDAYSTARTHG),ATMDEP_MEHG(N1)/(JDAY-JDAYSTARTHG),THgOUT/(JDAY-JDAYSTARTHG),THgIN/(JDAY-JDAYSTARTHG),THgPR/(JDAY-JDAYSTARTHG),THgWD/(JDAY-JDAYSTARTHG),THgDTRIB/(JDAY-JDAYSTARTHG),THgTRIB/(JDAY-JDAYSTARTHG),MeHgOUT/(JDAY-JDAYSTARTHG),MeHgIN/(JDAY-JDAYSTARTHG),MeHgPR/(JDAY-JDAYSTARTHG),MeHgWD/(JDAY-JDAYSTARTHG),MeHgDTRIB/(JDAY-JDAYSTARTHG),MeHgTRIB/(JDAY-JDAYSTARTHG)  
          else
          WRITE(FLUX_OUTPUT_FILE,'(F12.3,",",I2,",",25(E12.5,","))')JDAY,N1,(KF_HG_SUM(J1,N1)/(JDAY-JDAYSTARTHG),J1=1,10),ATMDEP_HG0(N1)/(JDAY-JDAYSTARTHG),ATMDEP_HG2(N1)/(JDAY-JDAYSTARTHG),ATMDEP_MEHG(N1)/(JDAY-JDAYSTARTHG),THgOUT/(JDAY-JDAYSTARTHG),THgIN/(JDAY-JDAYSTARTHG),THgPR/(JDAY-JDAYSTARTHG),THgWD/(JDAY-JDAYSTARTHG),THgDTRIB/(JDAY-JDAYSTARTHG),THgTRIB/(JDAY-JDAYSTARTHG),MeHgOUT/(JDAY-JDAYSTARTHG),MeHgIN/(JDAY-JDAYSTARTHG),MeHgPR/(JDAY-JDAYSTARTHG),MeHgWD/(JDAY-JDAYSTARTHG),MeHgDTRIB/(JDAY-JDAYSTARTHG),MeHgTRIB/(JDAY-JDAYSTARTHG)     
          endif
          
          ENDDO
          
          JDAYSTARTHG=JDAY                      
          KF_HG_SUM=0.0
          KF_HG2_SUM=0.0
          ATMDEP_HG0=0.0;ATMDEP_HG2=0.0;ATMDEP_MEHG=0.0;THgOUT=0.0;THgIN=0.0;THgPR=0.0;THgWD=0.0;THgDTRIB=0.0;THgTRIB=0.0;MeHgOUT=0.0;MeHgIN=0.0;MeHgPR=0.0;MeHgWD=0.0;MeHgDTRIB=0.0;MeHgTRIB=0.0
      ENDIF
      RETURN
      
  end subroutine HgKineticFluxes
  subroutine HgKineticFluxesOutputInit
    OPEN(FLUX_OUTPUT_FILE,FILE='HgFluxSummary.csv',STATUS='unknown')
          if (use_BedHg) then
          WRITE(FLUX_OUTPUT_FILE,'(A)')'JDAY,Region,Hg0 volatilization(ng/d),Hg0 oxidation(ng/d),HgII reduction(ng/d),HgII methylation(ng/d),HgII settling(ng/d),HgII sediment transfer(ng/d),MeHg degradation(ng/d),MeHg demethylation(ng/d),MeHg settling(ng/d),MeHg sediment transfer(ng/d),BedHgII methylation(ng/d),BedHgII burial(ng/d),BedSettlingAccumulationHgII(ng/d),BedMeHg demethylation(ng/d),BedMeHg burial(ng/d),BedSettlingAccumulationMeHg(ng/d),AtmDepositionHg0(ng/d),AtmDepositionHG2(ng/d),AtmDepositionMeHg(ng/d),THgOUT(ng/d),THgIN(ng/d),THgPR(ng/d),THgWD(ng/d),THgDTRIB(ng/d),THgTRIB(ng/d),MeHgOUT(ng/d),MeHgIN(ng/d),MeHgPR(ng/d),MeHgWD(ng/d),MeHgDTRIB(ng/d),MeHgTRIB(ng/d)'
          else
          WRITE(FLUX_OUTPUT_FILE,'(A)')'JDAY,Region,Hg0 volatilization(ng/d),Hg0 oxidation(ng/d),HgII reduction(ng/d),HgII methylation(ng/d),HgII settling(ng/d),HgII sediment transfer(ng/d),MeHg degradation(ng/d),MeHg demethylation(ng/d),MeHg settling(ng/d),MeHg sediment transfer(ng/d),AtmDepositionHg0(ng/d),AtmDepositionHG2(ng/d),AtmDepositionMeHg(ng/d),THgOUT(ng/d),THgIN(ng/d),THgPR(ng/d),THgWD(ng/d),THgDTRIB(ng/d),THgTRIB(ng/d),MeHgOUT(ng/d),MeHgIN(ng/d),MeHgPR(ng/d),MeHgWD(ng/d),MeHgDTRIB(ng/d),MeHgTRIB(ng/d)'   
          endif
  
  RETURN
  endsubroutine HgKineticFluxesOutputInit
  !===========================================================================================================================
  character*64 function AddIndex(basename, index)
    integer,       intent(in) :: index
    character*(*), intent(in) :: basename
    character*8 :: buffer
    
    write(buffer,'(I8)') index
    AddIndex = trim(basename) // trim(adjustl(buffer))
  end function
  !===========================================================================================================================
  subroutine HgModule_Deallocate
    deallocate(Hgd, Hgdoc, Hgpom, Hgap, Hgp, Hgpts, Hgpt)
    deallocate(Hgzp, Hgzpt, Hgzpts)
    deallocate(CON_HG, KF_HG, CNAMEHG)
    !
    deallocate(Kdoc) 
    deallocate(Kpom)  
    deallocate(Kp)
    deallocate(Kap)
    deallocate(Kzp)
    !
    deallocate(Hg00, vv_20, vv_theta)
    deallocate(k_oxid1, k_oxid2, k_oxid3)
    deallocate(k_reduct, kdoc_reduct)
    deallocate(k_meth,   kdoc_meth)
    deallocate(k_decay,  kdoc_decay)
    deallocate(k_demeth, kdoc_demeth)
    deallocate(vm,vm1)
    deallocate(KDOHg1, KDOHg2, KDOn)
    !
    if (use_BedHg) then 
      deallocate(Hg2, Hgd2, Hgdoc2, Hgpom2, Hgp2, Hgpts2, Hgpt2)
      deallocate(CON_HG2, KF_HG2, CNAMEHG2)
      !
      deallocate(Kdoc2)
      deallocate(Kpom2)
      deallocate(Kp2)
      deallocate(k_meth2,   kdoc_meth2)
      deallocate(k_demeth2, kdoc_demeth2)
      !
      deallocate(BedH, BedPor, BedPs, vb, vbs)
      deallocate(BedISS, dBedISSdt)
      if (.not. RESTART_IN) then
        deallocate(RegionS_ini, RegionE_ini, BedISS_Ini, Hg2_Ini)
      end if
    end if
    !
    deallocate(RegionS, RegionE, SegRegion)
    close(FLUX_OUTPUT_FILE);deallocate(KF_HG_SUM,KF_HG2_SUM,SREGION,EREGION)
  end subroutine
!===========================================================================================================================
  ! allocate a real array and set the default value
  subroutine alloc_R8(n, a, default)
    integer,               intent(in)    :: n   
    real(R8),              intent(in)    :: default
    real(R8), allocatable, intent(inout) :: a(:)
    
    if(n < 1) return
    if(allocated(a)) deallocate(a)
    allocate(a(n))    
    a = default
  end subroutine
  
END MODULE