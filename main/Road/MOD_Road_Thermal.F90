#include <define.h>

MODULE MOD_Road_Thermal

   USE MOD_Precision
   IMPLICIT NONE
   SAVE
   PRIVATE

   PUBLIC :: RoadTHERMAL

CONTAINS


   SUBROUTINE RoadTHERMAL ( &
        ! model running information
        ipatch         ,patchtype      ,lbroad         ,& !lbl            ,&
        deltim         ,patchlatr      ,&
        ! forcing
        forc_hgt_u     ,forc_hgt_t     ,forc_hgt_q     ,forc_us        ,&
        forc_vs        ,forc_t         ,forc_q         ,forc_psrf      ,&
        forc_rhoair    ,forc_frl       ,forc_po2m      ,forc_pco2m     ,&
        forc_sols      ,forc_soll      ,forc_solsd     ,forc_solld     ,&
        theta          ,& !sabwsun        ,sabwsha                        ,&
        sabroad        ,& !sablake        ,  sabv           ,
!        par            ,&
        ! surface parameters
!        flake          ,
        pondmx         ,eroad          ,trsmx0         ,&
        zlnd           ,zsno           ,capr           ,cnfac          ,&
        vf_quartz      ,vf_gravels     ,vf_om          ,vf_sand        ,&
        wf_gravels     ,wf_sand        ,csol           ,porsl          ,&
        psi0           ,&
!#ifdef Campbell_SOIL_MODEL
!        bsw            ,&
!#endif
!#ifdef vanGenuchten_Mualem_SOIL_MODEL
!        theta_r        ,alpha_vgm      ,n_vgm          ,L_vgm          ,&
!        sc_vgm         ,fc_vgm         ,&
!#endif
        k_solids       ,dksatu         ,dksatf         ,dkdry          ,&
        BA_alpha       ,BA_beta                                        ,&
        cv_road        ,tk_road                                        ,&
        dz_roadsno     ,&
        z_roadsno      ,&
        zi_roadsno     ,&
!        dz_lake        ,lakedepth      ,&
        ! vegetation parameters
!        dewmx          ,sqrtdi         ,rootfr         ,effcon         ,&
!        vmax25         ,slti           ,hlti           ,shti           ,&
!        hhti           ,trda           ,trdm           ,trop           ,&
!        g1             ,g0             ,gradm          ,binter         ,&
!        extkn          ,&
        ! surface status
        fsno_road      ,scv_road       ,& !scv_lake       ,
        snowdp_road    ,&
!        snowdp_lake    ,&
   !     lai            ,&
   !     sai            ,htop           ,hbot           ,& !sigf           ,&
   !     extkd          ,
        lroad          ,&
        t_road         ,t_roadsno      ,& !t_lakesno      ,
        wliq_roadsno   ,&
!        wliq_lakesno   ,&
        wice_roadsno   ,&  !wice_lakesno   ,t_lake         ,&
!        lake_icefrac   ,    savedtke1      ,lveg           ,tleaf          ,&
!        ldew           ,&
!! SNICAR model variables
!        snofrz         ,sabg_lyr                                       ,&
!! END SNICAR model variables
        ! output
        taux           ,tauy           ,fsena          ,fevpa          ,&
        lfevpa         ,& !fsenl       ,fevpl          ,etr            ,&
        fseng          ,fevpg          ,olrg           ,fgrnd          ,&
        fseng_road     ,lfevp_road     ,qseva_road     ,qseva_lake     ,&
        qsdew_road     ,qsubl_road     ,qfros_road     ,&
        imelt_road     ,& !imelt_lake     ,&
        sm_road        ,& !sm_lake        ,&
        sabg           ,rstfac         ,rootr          ,tref           ,&
        qref           ,trad           ,rst            ,assim          ,&
        respc          ,errore         ,emis           ,z0m            ,&
        zol            ,rib            ,ustar          ,qstar          ,&
        tstar          ,fm             ,fh             ,fq             ,&
        hpbl                                                            )

!=======================================================================
! this is the main subroutine to execute the calculation
! of thermal processes and surface fluxes on roads
!
!=======================================================================

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Vars_Global
   USE MOD_Const_Physical, only: denh2o,roverg,hvap,hsub,rgas,cpair,&
                                 stefnc,denice,tfrz,vonkar,grav
   USE MOD_Qsadv
   USE MOD_Road_GroundFlux
   USE MOD_Road_Flux
   USE MOD_Road_GroundTemperature
!   USE MOD_Lake
   USE MOD_Eroot, only: eroot
#ifdef vanGenuchten_Mualem_SOIL_MODEL
   USE MOD_Hydro_SoilFunction, only : soil_psi_from_vliq
#endif

   IMPLICIT NONE

!---------------------Argument------------------------------------------
   integer,  intent(in) :: &
        idate(3)   ,&
        ipatch     ,&! patch index
        patchtype  ,&! land patch type (0=soil, 1=urban or built-up, 2=wetland,
                     ! 3=glacier/ice sheet, 4=land water bodies)
        lbroad       ! lower bound of array
!        lbl          ! lower bound of array

   real(r8), intent(in) :: &
        deltim     ,&! seconds in a time step [second]
        patchlatr    ! latitude in radians

   real(r8), intent(in) :: &
        ! atmospherical variables and observational height
        forc_hgt_u ,&! observational height of wind [m]
        forc_hgt_t ,&! observational height of temperature [m]
        forc_hgt_q ,&! observational height of humidity [m]
        forc_us    ,&! wind component in eastward direction [m/s]
        forc_vs    ,&! wind component in northward direction [m/s]
        forc_t     ,&! temperature at agcm reference height [kelvin]
        forc_q     ,&! specific humidity at agcm reference height [kg/kg]
        forc_psrf  ,&! atmosphere pressure at the surface [pa]
        forc_rhoair,&! density air [kg/m3]
        forc_frl   ,&! atmospheric infrared (longwave) radiation [W/m2]
        forc_po2m  ,&! O2 concentration in atmos. (pascals)
        forc_pco2m ,&! CO2 concentration in atmos. (pascals)
        forc_sols  ,&! atm vis direct beam solar rad onto srf [W/m2]
        forc_soll  ,&! atm nir direct beam solar rad onto srf [W/m2]
        forc_solsd ,&! atm vis diffuse solar rad onto srf [W/m2]
        forc_solld ,&! atm nir diffuse solar rad onto srf [W/m2]
        theta      ,&! sun zenith angle
!        par        ,&! vegetation PAR
!        sabv       ,&! absorbed shortwave radiation by vegetation [W/m2]
!        sabwsun    ,&! absorbed shortwave radiation by sunlit wall [W/m2]
!        sabwsha    ,&! absorbed shortwave radiation by shaded wall [W/m2]
        sabroad       ! absorbed shortwave radiation by impervious road [W/m2]
!        sablake      ! absorbed shortwave radiation by lake [W/m2]

   real(r8), intent(in) :: &
!        flake      ,&! urban lake fractional cover [-]
        pondmx     ,&! maximum ponding for soil [mm]
        eroad      ,&! emissivity of road
        trsmx0     ,&! max transpiration for moist soil+100% veg.  [mm/s]
        zlnd       ,&! roughness length for soil [m]
        zsno       ,&! roughness length for snow [m]
        capr       ,&! tuning factor to turn first layer T into surface T
        cnfac      ,&! Crank Nicholson factor between 0 and 1

        ! soil physical parameters
        vf_quartz (1:nl_soil), &! volumetric fraction of quartz within mineral soil
        vf_gravels(1:nl_soil), &! volumetric fraction of gravels
        vf_om     (1:nl_soil), &! volumetric fraction of organic matter
        vf_sand   (1:nl_soil), &! volumetric fraction of sand
        wf_gravels(1:nl_soil), &! gravimetric fraction of gravels
        wf_sand   (1:nl_soil), &! gravimetric fraction of sand
        csol      (1:nl_soil), &! heat capacity of soil solids [J/(m3 K)]
        porsl     (1:nl_soil), &! soil porosity [-]
        psi0      (1:nl_soil), &! soil water suction, negative potential [mm]
        k_solids  (1:nl_soil), &! thermal conductivity of minerals soil [W/m-K]
        dkdry     (1:nl_soil), &! thermal conductivity of dry soil [W/m-K]
        dksatu    (1:nl_soil), &! thermal conductivity of saturated unfrozen soil [W/m-K]
        dksatf    (1:nl_soil), &! thermal conductivity of saturated frozen soil [W/m-K]

        BA_alpha  (1:nl_soil), &! alpha in Balland and Arp(2005) thermal conductivity scheme
        BA_beta   (1:nl_soil), &! beta in Balland and Arp(2005) thermal conductivity scheme
        cv_road   (1:nl_soil) ,&! heat capacity of road [J/(m2 K)]
        tk_road   (1:nl_soil) ,&! thermal conductivity of road [W/m-K]

        dz_roadsno(lbroad  :nl_soil) ,&! layer thickiness [m]
        z_roadsno (lbroad  :nl_soil) ,&! node depth [m]
        zi_roadsno(lbroad-1:nl_soil)   ! interface depth [m]
!        dz_lake   (    1:nl_lake) ,&! lake layer thickness (m)
!        lakedepth,                 &! lake depth (m)
!        z_lakesno (maxsnl+1:nl_soil) ,&! node depth [m]
!        dz_lakesno(maxsnl+1:nl_soil) ,&! layer thickiness [m]
!        zi_lakesno(maxsnl  :nl_soil) ,&! interface depth [m]

        ! vegetationparameters
!        dewmx      ,&! maximum dew
!        sqrtdi     ,&! inverse sqrt of leaf dimension [m**-0.5]
!        rootfr(1:nl_soil) ,&! root fraction

!        effcon     ,&! quantum efficiency of RuBP regeneration (mol CO2/mol quanta)
!        vmax25     ,&! maximum carboxylation rate at 25 C at canopy top
!        slti       ,&! slope of low temperature inhibition function      [s3]
!        hlti       ,&! 1/2 point of low temperature inhibition function  [s4]
!        shti       ,&! slope of high temperature inhibition function     [s1]
!        hhti       ,&! 1/2 point of high temperature inhibition function [s2]
!        trda       ,&! temperature coefficient in gs-a model             [s5]
!        trdm       ,&! temperature coefficient in gs-a model             [s6]
!        trop       ,&! temperature coefficient in gs-a model
!        g1         ,&! conductance-photosynthesis slope parameter for medlyn model
!        g0         ,&! conductance-photosynthesis intercept for medlyn model
!        gradm      ,&! conductance-photosynthesis slope parameter
!        binter     ,&! conductance-photosynthesis intercept
!        extkn        ! coefficient of leaf nitrogen allocation

   real(r8), intent(in) :: &
        fsno_road     ! fraction of ground covered by snow
!        lai        ,&! adjusted leaf area index for seasonal variation [-]
!        sai        ,&! stem area index  [-]
!        htop       ,&! canopy crown top height [m]
!        hbot       ,&! canopy crown bottom height [m]
!        fveg       ,&! fraction of veg cover
!        sigf       ,&! fraction of veg cover, excluding snow-covered veg [-]
!        extkd        ! diffuse and scattered diffuse PAR extinction coefficient

   real(r8), intent(in) :: hpbl       ! atmospheric boundary layer height [m]

   real(r8), intent(inout) :: &
        lroad     ,&! net longwave radiation of road
        t_road     ,&! ground temperature
        t_roadsno   (lbroad:nl_soil) ,&! temperatures of roof layers
        wliq_roadsno(lbroad:nl_soil) ,&! liqui water [kg/m2]
        wice_roadsno(lbroad:nl_soil) ,&! ice lens [kg/m2]
!        t_lake      (    nl_lake)    ,&! lake temperature [K]
!        lake_icefrac(    nl_lake)    ,&! lake mass fraction of lake layer that is frozen
!        t_lakesno   (maxsnl+1:nl_soil) ,&! temperatures of roof layers
!        wliq_lakesno(maxsnl+1:nl_soil) ,&! liqui water [kg/m2]
!        wice_lakesno(maxsnl+1:nl_soil) ,&! ice lens [kg/m2]
!        savedtke1  ,&! top level eddy conductivity (W/m K)
        scv_road   ,&! snow cover, water equivalent [mm, kg/m2]
!        scv_lake   ,&! snow cover, water equivalent [mm, kg/m2]
        snowdp_road   ! snow depth [m]
!        snowdp_lake,&! snow depth [m]
!        lveg       ,&! net longwave radiation of vegetation [W/m2]
!        tleaf      ,&! leaf temperature [K]
!        ldew         ! depth of water on foliage [kg/(m2 s)]

   ! Output
   real(r8), intent(out) :: &
        taux       ,&! wind stress: E-W [kg/m/s**2]
        tauy       ,&! wind stress: N-S [kg/m/s**2]
        fsena      ,&! sensible heat from canopy height to atmosphere [W/m2]
        fevpa      ,&! evapotranspiration from canopy height to atmosphere [mm/s]
        lfevpa     ,&! latent heat flux from canopy height to atmosphere [W/m2]
!        fsenl      ,&! ensible heat from leaves [W/m2]
!        fevpl      ,&! evaporation+transpiration from leaves [mm/s]
!        etr        ,&! transpiration rate [mm/s]
        fseng      ,&! sensible heat flux from ground [W/m2]
        fevpg      ,&! evaporation heat flux from ground [mm/s]
        olrg       ,&! outgoing long-wave radiation from ground+canopy
        fgrnd      ,&! ground heat flux [W/m2]

        fseng_road  ,&! sensible heat from road [W/m2]
        lfevp_road ,&! latent heat flux from road [W/m2]
        
        qseva_road ,&! ground soil surface evaporation rate (mm h2o/s)
!        qseva_lake ,&! ground soil surface evaporation rate (mm h2o/s)
        qsdew_road ,&! ground soil surface dew formation (mm h2o /s) [+]
!        qsdew_lake ,&! ground soil surface dew formation (mm h2o /s) [+]
        qsubl_road ,&! sublimation rate from soil ice pack (mm h2o /s) [+]
!        qsubl_lake ,&! sublimation rate from soil ice pack (mm h2o /s) [+]
        qfros_road    ! surface dew added to snow pack (mm h2o /s) [+]
!        qfros_lake   ! surface dew added to snow pack (mm h2o /s) [+]

   integer, intent(out) :: &
        imelt_road(lbroad:nl_soil)       ! flag for melting or freezing [-]

   real(r8), intent(out) :: &
        sm_road    ,&! rate of snowmelt [kg/(m2 s)]
!        sm_lake    ,&! rate of snowmelt [kg/(m2 s)]
        sabg       ,&! overall ground solar radiation absorption
        rstfac     ,&! factor of soil water stress
        rootr(1:nl_soil) ,&! root resistance of a layer, all layers add to 1
        tref       ,&! 2 m height air temperature [kelvin]
        qref       ,&! 2 m height air specific humidity
        trad       ,&! radiative temperature [K]
        rst        ,&! stomatal resistance (s m-1)
        assim      ,&! assimilation
        respc      ,&! respiration
        errore     ,&! energy balnce error [w/m2]

        ! additionalvariables required by coupling with WRF or RSM model
        emis       ,&! averaged bulk surface emissivity
        z0m        ,&! effective roughness [m]
        zol        ,&! dimensionless height (z/L) used in Monin-Obukhov theory
        rib        ,&! bulk Richardson number in surface layer
        ustar      ,&! u* in similarity theory [m/s]
        qstar      ,&! q* in similarity theory [kg/kg]
        tstar      ,&! t* in similarity theory [K]
        fm         ,&! integral of profile function for momentum
        fh         ,&! integral of profile function for heat
        fq           ! integral of profile function for moisture

! SNICAR model variables
!   real(r8), intent(in)  :: sabg_lyr(lbp:1) !snow layer aborption
!   real(r8), intent(out) :: snofrz (lbp:0)  !snow freezing rate (col,lyr) [kg m-2 s-1]
! END SNICAR model variables

!---------------------Local Variables-----------------------------------
   real(r8) :: &
      !  fg         ,&! ground fraction ( impervious + soil + snow )
        fsenroad   ,&! sensible heat flux from road [W/m2]
        fevproad   ,&! evaporation heat flux from impervious road [mm/s]

        croads     ,&! deriv of road sensible heat flux wrt soil temp [w/m**2/k]
        croadl     ,&! deriv of road latent heat flux wrt soil temp [w/m**2/k]
        croad      ,&! deriv of road total heat flux wrt soil temp [w/m**2/k]
        
        dqroaddT   ,&! d(qroad)/dT
      !  dqgperdT   ,&! d(qgper)/dT

        degdT      ,&! d(eg)/dT
        eg         ,&! water vapor pressure at temperature T [pa]
        egsmax     ,&! max. evaporation which soil can provide at one time step
        egidif     ,&! the excess of evaporation over "egsmax"
        emg        ,&! ground emissivity (0.97 for snow,
                     ! glaciers and water surface; 0.96 for soil and wetland)
!        etrc       ,&! maximum possible transpiration rate [mm/s]
        fac        ,&! soil wetness of surface layer
        fact_road(lbroad:nl_soil) ,&! used in computing tridiagonal matrix for road
        hr         ,&! relative humidity
        htvp_road  ,&! latent heat of vapor of water (or sublimation) [J/Kg]
        olru       ,&! olrg excluding dwonwelling reflection [W/m2]
        olrb       ,&! olrg assuming blackbody emission [W/m2]
        psit       ,&! negative potential of soil

        rsr        ,&! soil resistance
        qroad      ,&! ground impervious road specific humudity [kg/kg]
        q_snow     ,&! ground snow specific humudity [kg/kg]
        qsatg      ,&! saturated humidity [kg/kg]
        qsatgdT    ,&! d(qsatg)/dT
        qred       ,&! soil surface relative humidity
        thm        ,&! intermediate variable (forc_t+0.0098*forc_hgt_t)
        th         ,&! potential temperature (kelvin)
        thv        ,&! virtual potential temperature (kelvin)

        troad      ,&! temperature of road
        troad_bef  ,&! temperature of road at previous time step
        t_snow     ,&! ground snow temperature
        t_soisno_bef(lbroad:nl_soil), &! soil/snow temperature before update 
        tinc       ,&! temperature difference of two time step
!        ev         ,&! emissivity of vegetation [-]
        lout       ,&! out-going longwave radiation
        lnet       ,&! overall net longwave radiation
        lroad_bef  ,&! net longwave radiation of road at previous time step
        dlout      ,&! changed out-going radiation due to temp change
        clroad     ,&! deriv of lroad wrt gimp temp [w/m**2/k]
        
        ur         ,&! wind speed at reference height [m/s]
        wx         ,&! patitial volume of ice and water of surface layer
        xmf          ! total latent heat of phase change of ground water

   real(r8) :: z0m_g,z0h_g,zol_g,obu_g,ustar_g,qstar_g,tstar_g
   real(r8) :: fm10m,fm_g,fh_g,fq_g,fh2m,fq2m,um,obu,eb

   real(r8) :: dT      !temperature change between two time steps

!=======================================================================
! [1] Initial set and propositional variables
!=======================================================================

      ! fluxes
      taux   = 0.;  tauy   = 0.
      fsena  = 0.;  fevpa  = 0.
      lfevpa = 0.;  !fsenl  = 0.
      fseng  = 0.;  fevpg  = 0.
!      fsenl = 0.;   fevpl = 0.
!      etr   = 0.;   
      rst   = 2.0e4
      assim = 0.;   respc = 0.

      cgrnds = 0.;  cgrndl = 0.
      cgrnd  = 0.;  tref   = 0.
      qref   = 0.;  hprl   = 0.

      emis  = 0.;  z0m   = 0.
      zol   = 0.;  rib   = 0.
      ustar = 0.;  qstar = 0.
      tstar = 0.;  rootr = 0.

      t_snow = t_soisno(lbroad)
      ! latent heat, assumed that the sublimation occured only as wliq_gpersno=0
      htvp_road = hvap
      
      IF (wliq_roadsno(lbroad)<=0. .and. wice_roadsno(lbroad)>0.) htvp_road = hsub

      ! potential temperatur at the reference height
      thm = forc_t + 0.0098*forc_hgt_t                     !intermediate variable equivalent to
                                                           !forc_t*(pgcm/forc_psrf)**(rgas/cpair)
      th  = forc_t*(100000./forc_psrf)**(rgas/cpair)       !potential T
      thv = th*(1.+0.61*forc_q)                            !virtual potential T
      ur  = max(0.1,sqrt(forc_us*forc_us+forc_vs*forc_vs)) !limit set to 0.1

      ! temperature and water mass from previous time step
      troad = t_roadsno(lbroad)

      ! SAVE temperature
      troad_bef = troad

      ! SAVE longwave for the last time
      lroad_bef = lroad

!=======================================================================
! [2] specific humidity and its derivative at ground surface
!=======================================================================

!      qred = 1.
!      CALL qsadv(tgper,forc_psrf,eg,degdT,qsatg,qsatgdT)

!      ! initialization for rsr
!      rsr = 0.
!
!      IF (patchtype <=1 ) THEN          !soil ground
!         wx = (wliq_gpersno(1)/denh2o + wice_gpersno(1)/denice)/dz_gpersno(1)
!         IF (porsl(1) < 1.e-6) THEN     !bed rock
!            fac = 0.001
!         ELSE
!            fac = min(1.,wx/porsl(1))
!            fac = max( fac, 0.001 )
!         ENDIF
!
!#ifdef Campbell_SOIL_MODEL
!         psit = psi0(1) * fac ** (- bsw(1) )   !psit = max(smpmin, psit)
!#endif
!#ifdef vanGenuchten_Mualem_SOIL_MODEL
!         psit = soil_psi_from_vliq ( fac*(porsl(1)-theta_r(1)) + theta_r(1), &
!            porsl(1), theta_r(1), psi0(1), &
!            5, (/alpha_vgm(1), n_vgm(1), L_vgm(1), sc_vgm(1), fc_vgm(1)/))
!#endif
!         psit = max( -1.e8, psit )
!         hr   = exp(psit/roverg/tgper)
!         qred = (1.-fsno_gper)*hr + fsno_gper
!
!         IF (lbp == 1) THEN !no snow layer exist
!
!            ! calculate soil resistance for evaporation
!            wx   = (sum(wliq_gpersno(1:2))/denh2o + sum(wice_gpersno(1:2))/denice)/sum(dz_gpersno(1:2))
!            IF (sum(porsl(1:2)) < 1.e-6) THEN     !bed rock
!               fac  = 0.001
!            ELSE
!               fac  = min(1.,sum(dz_gpersno(1:2))*wx/(dz_gpersno(1)*porsl(1)+dz_gpersno(2)*porsl(2)))
!               fac  = max( fac, 0.001 )
!            ENDIF
!
!            ! Sellers et al., 1992
!            rsr = (1-fsno_gper)*exp(8.206-4.255*fac)
!         ENDIF
!      ENDIF
!
!      qgper = qred*qsatg
!      dqgperdT = qred*qsatgdT
!
!      IF (qsatg>forc_q .and. forc_q>qred*qsatg) THEN
!        qgper = forc_q; dqgperdT = 0.
!      ENDIF

      CALL qsadv(troad,forc_psrf,eg,degdT,qsatg,qsatgdT)
      qroad    = qsatg
      dqroaddT = qsatgdT

      IF (qsatg > forc_q .and. forc_q > qred*qsatg) THEN
         qroad = forc_q; dqgdT = 0.
      ENDIF

      q_snow = qroad
!!=======================================================================
!! [3] caluclate longwave radiation
!!=======================================================================
!
!      allocate ( Ainv(4,4)   )
!      allocate ( X(4)        )
!      allocate ( dX(4)       )
!      allocate ( B(4)        )
!      allocate ( B1(4)       )
!      allocate ( dBdT(4)     )
!      allocate ( SkyVF(4)    )
!      allocate ( fcover(0:4) )
!      allocate ( dT(0:4)     )
!
!      ! call longwave function, calculate Ainv, B, B1, dBdT
!      CALL UrbanOnlyLongwave ( &
!                              theta, hwr, froof, fgper, hroof, forc_frl, &
!                              twsun, twsha, tgimp, tgper, ewall, egimp, egper, &
!                              Ainv, B, B1, dBdT, SkyVF, fcover)
!
!      ! calculate longwave radiation abs, for UrbanOnlyLongwave
!      !-------------------------------------------
!      X = matmul(Ainv, B)
!
!      ! using the longwave radiation transfer matrix to calculate
!      ! LW radiation absorption by each surface and total absorption.
!      
!      lroad = ( eroad*X(3) - B1(3) ) / (1-eroad)
!      
!      ! Out-going LW of road
!      lout  = sum( X * SkyVF )
!
!      ! Energy balance check
!      eb = lroad + lout
!
!      IF (abs(eb-forc_frl) > 1e-6) THEN
!         print *, "Road Longwave - Energy Balance Check error!", eb-forc_frl
!      ENDIF
!
!      ! fur per unit surface
!   !   IF (fcover(3) >0.) lgimp = lgimp / fcover(3) * fg !/ fgimp
!   
!      ! added last time value
!      lroad = lroad + lroad_bef

!=======================================================================
! [3] Compute sensible and latent fluxes and their derivatives with respect
!     to ground temperature using ground temperatures from previous time step.
!=======================================================================

      ! bare ground case
      CALL RoadGroundFlux (forc_hgt_u,forc_hgt_t,forc_hgt_q,forc_us, &
                           forc_vs,forc_t,forc_q,forc_rhoair,forc_psrf, hpbl, &
                           ur,thm,th,thv,zlnd,zsno,fsno_road, lbroad, &
                           rss, dqroaddT, htvp, croad, croadl, croads, &
                           troad,qroad,t_snow,q_snow, &
                           taux,tauy,fseng_road,fseng_snow, &
                           fevpg_road,fevpg_snow,tref,qref, &
                           z0m_g,z0h_g,zol_g,rib_g,ustar_g,qstar_g,tstar_g,&
                           fm_g,fh_g,fq_g)

      ! SAVE variables for bareground case
      obu_g = forc_hgt_u / zol_g

!=======================================================================
! [4] Canopy temperature, fluxes from road
!=======================================================================


!         nurb = 2
!
!         ! CALL urban flux
!         CALL  UrbanOnlyFlux ( &
!            ! model running information
!            ipatch      ,deltim      ,lbr         ,lbi         ,&
!            ! forcing
!            forc_hgt_u  ,forc_hgt_t  ,forc_hgt_q  ,forc_us     ,&
!            forc_vs     ,thm         ,th          ,thv         ,&
!            forc_q      ,forc_psrf   ,forc_rhoair ,Fhac        ,&
!            Fwst        ,Fach        ,vehc        ,meta        ,&
!            ! surface parameters
!            hroof       ,hwr         ,nurb        ,fcover      ,&
!            ! surface status
!            z0h_g       ,obu_g       ,ustar_g     ,zlnd        ,&
!            zsno        ,fsno_roof   ,fsno_gimp   ,fsno_gper   ,&
!            wliq_roofsno(1),wliq_gimpsno(1),wice_roofsno(1),wice_gimpsno(1),&
!            htvp_roof   ,htvp_gimp   ,htvp_gper   ,troof       ,&
!            twsun       ,twsha       ,tgimp       ,tgper       ,&
!            qroof       ,qgimp       ,qgper       ,dqroofdT    ,&
!            dqgimpdT    ,dqgperdT    ,rsr                      ,&
!            ! output
!            taux        ,tauy        ,fsenroof    ,fsenwsun    ,&
!            fsenwsha    ,fsengimp    ,fsengper    ,fevproof    ,&
!            fevpgimp    ,fevpgper    ,croofs      ,cwalls      ,&
!            cgrnds      ,croofl      ,cgimpl      ,cgperl      ,&
!            croof       ,cgimp       ,cgper       ,tref        ,&
!            qref        ,z0m         ,zol         ,rib         ,&
!            ustar       ,qstar       ,tstar       ,fm          ,&
!            fh          ,fq          ,tafu                      )
!
!         !TODO: check
!         tleaf   = forc_t
!         ldew    = 0.
!         rstfac  = 0.
!         fsenl   = 0.0
!         fevpl   = 0.0
!         etr     = 0.0
!         assim   = 0.0
!         respc   = 0.0



!=======================================================================
! [6] road temperature
!=======================================================================

      CALL RoadTemperature (patchtype,lbroad,deltim,&
           capr,cnfac,vf_quartz,vf_gravels,vf_om,vf_sand, &
           wf_gravels,wf_sand,porsl,psi0,&
           csol,k_solids,dksatu,dksatf,dkdry,&
           BA_alpha,BA_beta,&
           cv_road,tk_road,&
           dz_roadsno,z_roadsno,zi_roadsno,&
           t_roadsno,wice_roadsno,wliq_roadsno,&
           scv_road,snowdp_road,&
           lroad,clroad,sabg_road,&
           fseng_road,fseng_snow,fevpg_road,fevpg_snow,&
           croad,htvp_road,&
           imelt_road,sm_road,xmf,fact_road)

      ! update temperature
      troad = t_roadsno(lbroad)

!=======================================================================
! [7] Correct fluxes for temperature change
!=======================================================================

      ! calculate temperature change
      dT = troad - troad_bef

      ! flux change due to temperture change
      fseng_road = fseng_road + dT*croads
    
      fevpg_road = fevpg_road + dT*croadl

! calculation of evaporative potential; flux in kg m-2 s-1.
! egidif holds the excess energy IF all water is evaporated
! during the timestep.  this energy is later added to the sensible heat flux.


      ! --- for impervious ground ---
      ! update of snow
      IF (lbroad < 1) THEN
         egsmax = (wice_roadsno(lbroad)+wliq_roadsno(lbroad)) / deltim
         egidif = max( 0., fevpg_road - egsmax )
         fevpg_road = min ( fevpg_road, egsmax )
         fseng_road = fseng_road + htvp_road*egidif
      ENDIF

      ! update of soil
      egsmax = (wice_roadsno(1)+wliq_roadsno(1)) / deltim
      egidif = max( 0., fevpg_road - egsmax )
      fevpg_road = min ( fevpg_road, egsmax )
      fseng_road = fseng_road + htvp_road*egidif

!=======================================================================
! [8] total fluxes to atmosphere
!=======================================================================

      lnet  = lroad

      sabg  = sab_road

      fseng = fseng_road
      
      fsen_road = fseng_road
      
      fevpg = fevpg_road

      lfevpa = htvp_road*fevpg_road     

      lfevp_road = htvp_road*fevpg_road

      fsena  = fseng
      fevpa  = fevpg

      ! 10/01/2021, yuan: exclude lake fevpa.
      ! because we don't consider water balance for lake currently.
      !fevpa  = fevpa *(1-flake) + fevpa_lake *flake

      ! 07/11/2023, yuan: don't not consider lake fraction cover
      !fsenl  = fsenl *(1-flake)
      !fevpl  = fevpl *(1-flake)
      !etr    = etr   *(1-flake)
      !assim  = assim *(1-flake)
      !respc  = respc *(1-flake)

      ! ground heat flux
      fgrnd = sabg + lnet - (fsena+lfevpa)

      ! effective ground temperature, simple average
      ! 12/01/2021, yuan: !TODO Bugs. temperature cannot be weighted like below.
      !t_grnd = troof*fcover(0) + twsun*fcover(1) + twsha*fcover(2) + &
      !t_grnd = tgper*fgper + tgimp*(1-fgper)
      t_road = troad

      !==============================================
      qseva_road = 0.
      qsubl_road = 0.
      qfros_road = 0.
      qsdew_road = 0.

      IF (fevpg_road >= 0.)THEN
! not allow for sublimation in melting (melting ==> evap. ==> sublimation)
         qseva_road = min(wliq_roadsno(lbroad)/deltim, fevpg_road)
         qsubl_road = fevpg_road - qseva_road
      ELSE
         IF (troad < tfrz)THEN
            qfros_road = abs(fevpg_road)
         ELSE
            qsdew_road = abs(fevpg_road)
         ENDIF
      ENDIF


!=======================================================================
! [9] Calculate the change rate of long-wave radiation caused by temperature change
!=======================================================================
!
!      dX = matmul(Ainv, dBdT*dT(1:))
!      lwsun = ( ewall*dX(1) - dBdT(1)*dT(1) ) / (1-ewall)
!      lwsha = ( ewall*dX(2) - dBdT(2)*dT(2) ) / (1-ewall)
!      lgimp = ( egimp*dX(3) - dBdT(3)*dT(3) ) / (1-egimp)
!      lgper = ( egper*dX(4) - dBdT(4)*dT(4) ) / (1-egper)
!
!      IF ( doveg ) THEN
!         lveg = ( sum(dX(1:5)*VegVF(1:5))*ev )
!      ELSE
!         lveg = 0.
!      ENDIF
!
!      dlout = sum( dX * SkyVF )
!
!      ! Energy balance check
!      eb = lwsun + lwsha + lgimp + lgper + lveg + dlout
!
!      IF (abs(eb) > 1e-6) THEN
!         print *, "Urban Vegetation Longwave - Energy Balance Check error!", eb
!      ENDIF
!
!      ! for per unit surface
!      IF (fcover(1) > 0.) lwsun = lwsun / fcover(1) * fg !/ (4*fwsun*HL*fb/fg)
!      IF (fcover(2) > 0.) lwsha = lwsha / fcover(2) * fg !/ (4*fwsha*HL*fb/fg)
!      IF (fcover(3) > 0.) lgimp = lgimp / fcover(3) * fg !/ fgimp
!      IF (fcover(4) > 0.) lgper = lgper / fcover(4) * fg !/ fgper
!      IF ( doveg        ) lveg  = lveg  / fcover(5) * fg !/ fv/fg

      ! calculate out going longwave by added the before value
      ! of lout and condsidered troof change
!      lout = lout + dlout
!      rout = (1-eroof)*forc_frl + eroof*stefnc*troof_bef**4 &
!           + 4.*eroof*stefnc*troof_bef**3*dT(0)
!
!      olrg = lout*fg + rout*froof
!      olrg = olrg*(1-flake) + olrg_lake*flake
!
!      !print*, forc_t, tgper, tgimp, troof, twsha, twsun
!
!      IF (olrg < 0) THEN !fordebug
!         print*, ipatch, olrg
!         write(6,*) ipatch,sabv,sabg,forc_frl,olrg,fsenl,fseng,hvap*fevpl,lfevpa
!         CALL CoLM_stop()
!      ENDIF

!      ! radiative temperature
!      trad = (olrg/stefnc)**0.25

!! averaged bulk surface emissivity
!!TODO: how to calculate for urban case?
!! 03/10/2020, yuan: removed below.
!      !olrb = stefnc*t_soisno_bef(lb)**3*(4.*tinc)
!      !olrb = stefnc*t_grnd_bef**3*(4.*tinc)
!      !olru = ulrad + emg*olrb
!      !olrb = ulrad + olrb
!      !emis = olru / olrb

!=======================================================================
! [10] energy balance error
!=======================================================================

!      IF ( doveg ) THEN
!         errore = sabv*fveg*(1-flake) + sabg + lnet - fsena - lfevpa - fgrnd
!      ELSE
      errore = sabg + lnet - fsena - lfevpa - fgrnd
!      ENDIF

!      ! deallocate memory
!      deallocate ( Ainv   )
!      deallocate ( X      )
!      deallocate ( dX     )
!      deallocate ( B      )
!      deallocate ( B1     )
!      deallocate ( dBdT   )
!      deallocate ( SkyVF  )
!      deallocate ( dT     )

!      IF ( doveg ) THEN
!         deallocate ( VegVF )
!      ENDIF

#if (defined CoLMDEBUG)
      IF (abs(errore)>.5) THEN
      write(6,*) 'RoadTHERMAL.F90: energy balance violation'
      write(6,*) ipatch,errore,sabg,forc_frl,olrg,fsenl,fseng,hvap*fevpl,lfevpa,xmf
      ENDIF
100   format(10(f15.3))
#endif

!      ! diagnostic sabg only for pervious and impervious ground
!      sabg = sabgper*fgper + sabgimp*(1-fgper)
!
!
!      deallocate ( fcover )

   END SUBROUTINE RoadTHERMAL

END MODULE MOD_Road_Thermal
! ---------- EOP ------------
