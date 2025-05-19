#include <define.h>

SUBROUTINE CoLMMain_Road ( &
         ! model running information
           ipatch       ,idate        ,coszen       ,deltim       ,&
           patchlonr    ,patchlatr    ,patchclass   ,patchtype    ,&

         ! road information
           em_road      ,cv_road      ,tk_road      ,alb_road     ,&
    !       z_road       ,dz_road                                  ,&
         ! soil information
           vf_quartz    ,vf_gravels   ,vf_om        ,vf_sand      ,&
           wf_gravels   ,wf_sand      ,porsl        ,psi0         ,&
           bsw          ,theta_r      ,fsatmax      ,fsatdcf      ,&
#ifdef vanGenuchten_Mualem_SOIL_MODEL
           alpha_vgm    ,n_vgm        ,L_vgm        ,&
           sc_vgm       ,fc_vgm       ,&
#endif
           hksati       ,csol         ,k_solids     ,dksatu       ,&
           dksatf       ,dkdry        ,BA_alpha     ,BA_beta      ,&
     
         ! vegetation information
!           sqrtdi       ,chil         ,&
!           effcon       ,vmax25       ,slti         ,hlti         ,&
!           shti         ,hhti         ,trda         ,trdm         ,&
!           trop         ,g1           ,g0           ,gradm        ,&
!           binter       ,extkn        ,rho          ,tau          ,&
!           rootfr                     ,&

         ! atmospheric forcing
           forc_us      ,forc_vs      ,&
           forc_t       ,forc_q       ,forc_prc     ,forc_prl     ,&
           forc_rain    ,forc_snow    ,forc_psrf    ,forc_pbot    ,&
           forc_sols    ,forc_soll    ,forc_solsd   ,forc_solld   ,&
           forc_frl     ,forc_hgt_u   ,forc_hgt_t   ,forc_hgt_q   ,&
           forc_rhoair  ,&
            ! cbl forcing
           !forc_hpbl,    &
           ! aerosol deposition
         !  forc_aerdep,  &

         ! land surface variables required for restart
           z_sno_road   ,dz_sno_road  ,t_roadsno    ,&
           wliq_roadsno ,wice_roadsno ,&
           z_sno        ,dz_sno       ,&
!           wliq_soisno  ,wice_soisno  ,t_soisno     ,&
           smp          ,hk           ,t_grnd       ,&
          
!           tleaf,        ldew,         ldew_rain,    ldew_snow,    &
           sag,          scv,          snowdp,       & !fveg,         &
           fsno,         & !sigf,         
!           green,        lai,          &
!           sai,          
           alb,          & !ssun,         ssha,         &
           ssoi,         ssno,         &
          ! thermk,       extkb,        &
          ! extkd,        vegwp,        gs0sun,       gs0sha,       &
           zwt,          & !wdsrf,
           wa,           & !wetwat,       &
           sag_road     ,scv_road     ,&
           snowdp_road  ,fsno_road    ,&
           sroad        ,lroad        ,&

!         ! SNICAR snow model related
!           snw_rds,      ssno_lyr,                                 &
!           mss_bcpho,    mss_bcphi,   mss_ocpho,     mss_ocphi,    &
!           mss_dst1,     mss_dst2,    mss_dst3,      mss_dst4,     &

         ! additional diagnostic variables for output
!           laisun       ,laisha       ,&
           rss                        ,&
           rstfac       ,h2osoi       ,wat                        ,&

         ! FLUXES
           taux         ,tauy         ,fsena        ,fevpa        ,&
           lfevpa       ,& !fsenl     ,fevpl        ,etr          ,&
           fseng        ,fevpg        ,olrg         ,fgrnd        ,&
           fsen_road    ,lfevp_road   ,&
           trad         ,tref         ,&!tmax       ,tmin         ,&
           qref         ,rsur         ,rnof         ,qintr        ,&
           qinfl        ,qdrip        ,rst          ,assim        ,&
           respc        ,& !sabvsun      ,sabvsha      ,
           sabg         ,&
           sr           ,solvd        ,solvi        ,solnd        ,&
           solni        ,srvd         ,srvi         ,srnd         ,&
           srni         ,solvdln      ,solviln      ,solndln      ,&
           solniln      ,srvdln       ,srviln       ,srndln       ,&
           srniln       ,qcharge      ,xerr         ,zerr         ,&

         ! TUNABLE modle constants
           zlnd         ,zsno         ,csoilc       ,dewmx        ,&
         !  wtfact       ,
           capr         ,cnfac        ,ssi          ,&
           wimp         ,pondmx       ,smpmax       ,smpmin       ,&
           trsmx0       ,tcrit                                    ,&

         ! additional variables required by coupling with WRF model
           emis         ,z0m          ,zol          ,rib          ,&
           ustar        ,qstar        ,tstar        ,fm           ,&
           fh           ,fq           ,hpbl                       )

  USE MOD_Precision
  USE MOD_Vars_Global
  USE MOD_Const_Physical, only: tfrz, denh2o, denice

  USE MOD_SnowLayersCombineDivide
  USE MOD_TimeManager  ! including isgreenwich
  USE MOD_RainSnowTemp, only: rain_snow_temp
  USE MOD_NewSnow, only: newsnow
  USE MOD_OrbCoszen, only: orb_coszen
  USE MOD_SnowFraction, only: snowfraction
  USE MOD_ALBEDO, only: snowage
  USE MOD_Qsadv, only: qsadv
  USE MOD_Road_Hydrology
  Use MOD_Road_Thermal
  USE MOD_Road_Albedo
  USE MOD_Road_CDP_SnowClear

  IMPLICIT NONE

! ------------------------ Dummy Argument ------------------------------
  integer, intent(in) :: &
      ipatch     ,&! maximum number of snow layers
      idate(3)   ,&! next time-step /year/julian day/second in a day/
      patchclass ,&! land cover type of USGS classification or others
      patchtype    ! land patch type (0=soil, 1=urban and built-up,
                   ! 2=wetland, 3=land ice, 4=land water bodies, 99 = ocean)

  real(r8), intent(in) :: &
      deltim     ,&! seconds in a time step [second]
      patchlonr  ,&! logitude in radians
      patchlatr    ! latitude in radians

  real(r8), intent(inout) :: &
      coszen       ! cosine of solar zenith angle

! Parameters
! ----------------------
  real(r8), intent(in) :: &
      em_road           , &! emissivity of road [-]
      cv_road(1:nl_soil), &! heat capacity of road [J/(m2 K)] 
      tk_road(1:nl_soil), &! thermal conductivity of road [W/m-K]
      alb_road(2,2)        ! albedo of road [-]

  real(r8), intent(in) :: &
    ! soil physical parameters
      vf_quartz (nl_soil),&! volumetric fraction of quartz within mineral soil
      vf_gravels(nl_soil),&! volumetric fraction of gravels
      vf_om     (nl_soil),&! volumetric fraction of organic matter
      vf_sand   (nl_soil),&! volumetric fraction of sand
      wf_gravels(nl_soil),&! gravimetric fraction of gravels
      wf_sand   (nl_soil),&! gravimetric fraction of sand
      porsl     (nl_soil),&! fraction of soil that is voids [-]
      psi0      (nl_soil),&! minimum soil suction [mm]
      bsw       (nl_soil),&! clapp and hornbereger "b" parameter [-]
      theta_r   (nl_soil),&! residual water content (cm3/cm3)   
      fsatmax            ,&! maximum saturated area fraction [-]
      fsatdcf            ,&! decay factor in calucation of saturated area fraction [1/m]

#ifdef vanGenuchten_Mualem_SOIL_MODEL
      alpha_vgm(1:nl_soil),&! the parameter corresponding approximately to the inverse of the air-entry value
      n_vgm    (1:nl_soil),&! a shape parameter
      L_vgm    (1:nl_soil),&! pore-connectivity parameter
      sc_vgm   (1:nl_soil),&! saturation at the air entry value in the classical vanGenuchten model [-]
      fc_vgm   (1:nl_soil),&! a scaling factor by using air entry value in the Mualem model [-]
#endif
      hksati    (nl_soil),&! hydraulic conductivity at saturation [mm h2o/s]
      csol      (nl_soil),&! heat capacity of soil solids [J/(m3 K)]
      k_solids  (nl_soil),&! thermal conductivity of minerals soil [W/m-K]
      dksatu    (nl_soil),&! thermal conductivity of saturated unfrozen soil [W/m-K]
      dksatf    (nl_soil),&! thermal conductivity of saturated frozen soil [W/m-K]
      dkdry     (nl_soil),&! thermal conductivity for dry soil  [J/(K s m)]

      BA_alpha  (nl_soil),&! alpha in Balland and Arp(2005) thermal conductivity scheme
      BA_beta   (nl_soil),&! beta in Balland and Arp(2005) thermal conductivity scheme

    ! vegetation static, dynamic, derived parameters
!      sqrtdi     ,&! inverse sqrt of leaf dimension [m**-0.5]
!      chil       ,&! leaf angle distribution factor
!      effcon     ,&! quantum efficiency of RuBP regeneration (mol CO2/mol quanta)
!      vmax25     ,&! maximum carboxylation rate at 25 C at canopy top
!      slti       ,&! slope of low temperature inhibition function      [s3]
!      hlti       ,&! 1/2 point of low temperature inhibition function  [s4]
!      shti       ,&! slope of high temperature inhibition function     [s1]
!      hhti       ,&! 1/2 point of high temperature inhibition function [s2]
!      trda       ,&! temperature coefficient in gs-a model             [s5]
!      trdm       ,&! temperature coefficient in gs-a model             [s6]
!      trop       ,&! temperature coefficient in gs-a model
!      g1         ,&! conductance-photosynthesis slope parameter for medlyn model
!      g0         ,&! conductance-photosynthesis intercept for medlyn model
!      gradm      ,&! conductance-photosynthesis slope parameter
!      binter     ,&! conductance-photosynthesis intercep
!      extkn      ,&! coefficient of leaf nitrogen allocation
!      rho(2,2)   ,&! leaf reflectance (iw=iband, il=life and dead)
!      tau(2,2)   ,&! leaf transmittance (iw=iband, il=life and dead)

!      rootfr    (nl_soil),&! fraction of roots in each soil layer

      ! tunable parameters
      zlnd       ,&! roughness length for soil [m]
      zsno       ,&! roughness length for snow [m]
      csoilc     ,&! drag coefficient for soil under canopy [-]
      dewmx      ,&! maximum dew
      !wtfact     ,&! fraction of model area with high water table
      capr       ,&! tuning factor to turn first layer T into surface T
      cnfac      ,&! Crank Nicholson factor between 0 and 1
      ssi        ,&! irreducible water saturation of snow
      wimp       ,&! water impremeable IF porosity less than wimp
      pondmx     ,&! ponding depth (mm)
      smpmax     ,&! wilting point potential in mm
      smpmin     ,&! restriction for min of soil poten.  (mm)
      trsmx0     ,&! max transpiration for moist soil+100% veg.  [mm/s]
      tcrit        ! critical temp. to determine rain or snow

! Forcing
! ----------------------
  real(r8), intent(in) :: &
!      forc_pco2m ,&! partial pressure of CO2 at observational height [pa]
!      forc_po2m  ,&! partial pressure of O2 at observational height [pa]
      forc_us    ,&! wind speed in eastward direction [m/s]
      forc_vs    ,&! wind speed in northward direction [m/s]
      forc_t     ,&! temperature at agcm reference height [kelvin]
      forc_q     ,&! specific humidity at agcm reference height [kg/kg]
      forc_prc   ,&! convective precipitation [mm/s]
      forc_prl   ,&! large scale precipitation [mm/s]
      forc_psrf  ,&! atmosphere pressure at the surface [pa]
      forc_pbot  ,&! atmosphere pressure at the bottom of the atmos. model level [pa]
      forc_sols  ,&! atm vis direct beam solar rad onto srf [W/m2]
      forc_soll  ,&! atm nir direct beam solar rad onto srf [W/m2]
      forc_solsd ,&! atm vis diffuse solar rad onto srf [W/m2]
      forc_solld ,&! atm nir diffuse solar rad onto srf [W/m2]
      forc_frl   ,&! atmospheric infrared (longwave) radiation [W/m2]
      forc_hgt_u ,&! observational height of wind [m]
      forc_hgt_t ,&! observational height of temperature [m]
      forc_hgt_q ,&! observational height of humidity [m]
      forc_rhoair  ! density air [kg/m3]
!      forc_hpbl    ! atmospheric boundary layer height [m]
!      forc_aerdep(14)!atmospheric aerosol deposition data [kg/m/s]

! Variables required for restart run
! ----------------------------------------------------------------------
  real(r8), intent(inout) :: &
!        t_soisno    (maxsnl+1:nl_soil) ,&! soil + snow layer temperature [K]
        t_roadsno   (maxsnl+1:nl_soil) ,&! soil + road + snow layer temperature [K]
!        t_lakesno   (maxsnl+1:nl_soil) ,&! soil + lake + snow layer temperature [K]
!        wliq_soisno (maxsnl+1:nl_soil) ,&! liquid water (kg/m2)
        wliq_roadsno(maxsnl+1:nl_soil) ,&! liquid water (kg/m2)
!        wliq_lakesno(maxsnl+1:nl_soil) ,&! liquid water (kg/m2)
!        wice_soisno (maxsnl+1:nl_soil) ,&! ice lens (kg/m2)
        wice_roadsno(maxsnl+1:nl_soil) ,&! ice lens (kg/m2)
!        wice_lakesno(maxsnl+1:nl_soil) ,&! ice lens (kg/m2)
        smp         (       1:nl_soil) ,&! soil matrix potential [mm]
        hk          (       1:nl_soil) ,&! hydraulic conductivity [mm h2o/s]

        z_sno       (maxsnl+1:0)       ,&! node depth [m]
        dz_sno      (maxsnl+1:0)       ,&! interface depth [m]
        z_sno_road  (maxsnl+1:0)       ,&! node depth of road [m]
!        z_sno_lake  (maxsnl+1:0)       ,&! node depth lake [m]
        dz_sno_road (maxsnl+1:0)       ,&! interface depth of road [m]
!        dz_sno_lake (maxsnl+1:0)       ,&! interface depth lake [m]

!        lakedepth             ,&! lake depth (m)
!        dz_lake     (nl_lake) ,&! lake layer thickness (m)
!        t_lake      (nl_lake) ,&! lake temperature (kelvin)
!        lake_icefrac(nl_lake) ,&! lake mass fraction of lake layer that is frozen
!        savedtke1             ,&! top level eddy conductivity (W/m K)

!        topostd    ,&! standard deviation of elevation [m]
!        BVIC      ,& ! b parameter in Fraction of saturated soil in a grid calculated by VIC

        t_grnd     ,&! ground surface temperature [k]
!        tleaf      ,&! sunlit leaf temperature [K]
        !tmax       ,&! Diurnal Max 2 m height air temperature [kelvin]
        !tmin       ,&! Diurnal Min 2 m height air temperature [kelvin]
!        ldew       ,&! depth of water on foliage [kg/m2/s]
        sag        ,&! non dimensional snow age [-]
        sag_road   ,&! non dimensional snow age [-]
        scv        ,&! snow mass (kg/m2)
        scv_road   ,&! snow mass (kg/m2)
        snowdp     ,&! snow depth (m)
        snowdp_road,&! snow depth (m)
        zwt        ,&! the depth to water table [m]
        wa         ,&! water storage in aquifer [mm]

!        snw_rds   ( maxsnl+1:0 ) ,&! effective grain radius (col,lyr) [microns, m-6]
!        mss_bcpho ( maxsnl+1:0 ) ,&! mass of hydrophobic BC in snow  (col,lyr) [kg]
!        mss_bcphi ( maxsnl+1:0 ) ,&! mass of hydrophillic BC in snow (col,lyr) [kg]
!        mss_ocpho ( maxsnl+1:0 ) ,&! mass of hydrophobic OC in snow  (col,lyr) [kg]
!        mss_ocphi ( maxsnl+1:0 ) ,&! mass of hydrophillic OC in snow (col,lyr) [kg]
!        mss_dst1  ( maxsnl+1:0 ) ,&! mass of dust species 1 in snow  (col,lyr) [kg]
!        mss_dst2  ( maxsnl+1:0 ) ,&! mass of dust species 2 in snow  (col,lyr) [kg]
!        mss_dst3  ( maxsnl+1:0 ) ,&! mass of dust species 3 in snow  (col,lyr) [kg]
!        mss_dst4  ( maxsnl+1:0 ) ,&! mass of dust species 4 in snow  (col,lyr) [kg]
!        ssno    (2,2,maxsnl+1:1) ,&! snow layer absorption [-]

!        fveg       ,&! fraction of vegetation cover
        fsno       ,&! fractional snow cover
        fsno_road  ,&! fractional snow cover
!        fsno_lake  ,&! fractional snow cover
!        sigf       ,&! fraction of veg cover, excluding snow-covered veg [-]
!        green      ,&! greenness
!        lai        ,&! leaf area index
!        sai        ,&! stem area index
!        htop       ,&! canopy crown top
!        hbot       ,&! canopy crown bottom

        !coszen      ,&! cosine of solar zenith angle
        alb(2,2)    ,&! averaged albedo [-]
      !  ssun(2,2)   ,&! sunlit canopy absorption for solar radiation
      !  ssha(2,2)   ,&! shaded canopy absorption for solar radiation
        ssoi(2,2)   ,&! ground soil absorption [-]
        ssno(2,2)   ,&! ground snow absorption [-]
      !  thermk      ,&! canopy gap fraction for tir radiation
      !  extkb       ,&! (k, g(mu)/mu) direct solar extinction coefficient
      !  extkd       ,&! diffuse and scattered diffuse PAR extinction coefficient
        lroad       ,&! net longwave of road  [W/m2]
!        lveg       ,&! net longwave of vegetation [W/m2]
        
!        ssun (2,2) ,&! sunlit canopy absorption for solar radiation
!        ssha (2,2) ,&! shaded canopy absorption for solar radiation
        sroad(2,2)    ! sunlit road absorption for solar radiation
!        slake(2,2)   ! sunlit lake absorption for solar radiation

! additional diagnostic variables for output
  real(r8), intent(out) :: &
!        laisun     ,&! sunlit leaf area index
!        laisha     ,&! shaded leaf area index
        rstfac     ,&! factor of soil water stress
        rss        ,&! soil surface resistance
        wat        ,&! total water storage
        h2osoi(nl_soil)! volumetric soil water in layers [m3/m3]

! Fluxes
! ----------------------------------------------------------------------
  real(r8), intent(out) :: &
        taux       ,&! wind stress: E-W [kg/m/s**2]
        tauy       ,&! wind stress: N-S [kg/m/s**2]
        fsena      ,&! sensible heat from canopy height to atmosphere [W/m2]
        fevpa      ,&! evapotranspiration from canopy height to atmosphere [mm/s]
        lfevpa     ,&! latent heat flux from canopy height to atmosphere [W/2]
!        fsenl      ,&! ensible heat from leaves [W/m2]
!        fevpl      ,&! evaporation+transpiration from leaves [mm/s]
!        etr        ,&! transpiration rate [mm/s]
        fseng      ,&! sensible heat flux from ground [W/m2]
        fevpg      ,&! evaporation heat flux from ground [mm/s]
        olrg       ,&! outgoing long-wave radiation from ground+canopy
        fgrnd      ,&! ground heat flux [W/m2]
        xerr       ,&! water balance error at current time-step [mm/s]
        zerr       ,&! energy balnce errore at current time-step [W/m2]

        tref       ,&! 2 m height air temperature [K]
        qref       ,&! 2 m height air specific humidity
        trad       ,&! radiative temperature [K]
        rsur       ,&! surface runoff (mm h2o/s)
        rnof       ,&! total runoff (mm h2o/s)
        qintr      ,&! interception (mm h2o/s)
        qinfl      ,&! inflitration (mm h2o/s)
        qdrip      ,&! throughfall (mm h2o/s)
        qcharge    ,&! groundwater recharge [mm/s]

        rst        ,&! canopy stomatal resistance
        assim      ,&! canopy assimilation
        respc      ,&! canopy respiration

        fsen_road  ,&! sensible heat flux from road [W/m2]
        lfevp_road ,&! latent heat flux from road [W/m2]

!        sabvsun    ,&! solar absorbed by sunlit vegetation [W/m2]
!        sabvsha    ,&! solar absorbed by shaded vegetation [W/m2]
        sabg       ,&! solar absorbed by ground  [W/m2]
        sr         ,&! total reflected solar radiation (W/m2)
        solvd      ,&! incident direct beam vis solar radiation (W/m2)
        solvi      ,&! incident diffuse beam vis solar radiation (W/m2)
        solnd      ,&! incident direct beam nir solar radiation (W/m2)
        solni      ,&! incident diffuse beam nir solar radiation (W/m2)
        srvd       ,&! reflected direct beam vis solar radiation (W/m2)
        srvi       ,&! reflected diffuse beam vis solar radiation (W/m2)
        srnd       ,&! reflected direct beam nir solar radiation (W/m2)
        srni       ,&! reflected diffuse beam nir solar radiation (W/m2)
        solvdln    ,&! incident direct beam vis solar radiation at local noon (W/m2)
        solviln    ,&! incident diffuse beam vis solar radiation at local noon (W/m2)
        solndln    ,&! incident direct beam nir solar radiation at local noon (W/m2)
        solniln    ,&! incident diffuse beam nir solar radiation at local noon (W/m2)
        srvdln     ,&! reflected direct beam vis solar radiation at local noon (W/m2)
        srviln     ,&! reflected diffuse beam vis solar radiation at local noon (W/m2)
        srndln     ,&! reflected direct beam nir solar radiation at local noon (W/m2)
        srniln     ,&! reflected diffuse beam nir solar radiation at local noon (W/m2)

        forc_rain  ,&! rain [mm/s]
        forc_snow  ,&! snow [mm/s]

! additional variables required by coupling with WRF model
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
  
  real(r8), intent(in) :: hpbl       ! atmospheric boundary layer height [m]

! ----------------------- Local  Variables -----------------------------
  integer  :: local_secs
  real(r8) :: radpsec

  real(r8) :: &
        calday     ,&! Julian cal day (1.xx to 365.xx)
        endwb      ,&! water mass at the end of time step
        errore     ,&! energy balnce error (Wm-2)
        errorw     ,&! water balnce error (mm)
 
        fioldroad(maxsnl+1:nl_soil), &! fraction of ice relative to the total water on road
!        fioldl(maxsnl+1:nl_soil), &! fraction of ice relative to the total water on lake
        w_old      ,&! liquid water mass of the column at the previous time step (mm)
        theta      ,&! sun zenith angle
        !orb_coszen ,&! cosine of the solar zenith angle
!        sabv       ,&! solar absorbed by vegetation [W/m2]
        sabroad    ,&! solar absorbed by road [W/m2]
!        sablake    ,&! solar absorbed by vegetation [W/m2]
!        par        ,&! PAR by leaves [W/m2]
        troad      ,&! surface temperature of road [K]
!        tlake      ,&! temperature of lake surface [K]
        qseva_road ,&! ground surface evaporation rate (mm h2o/s)
!        qseva_lake ,&! ground surface evaporation rate (mm h2o/s)
        qsdew_road ,&! ground surface dew formation (mm h2o /s) [+]
!        qsdew_lake ,&! ground surface dew formation (mm h2o /s) [+]
        qsubl_road ,&! sublimation rate from snow pack (mm h2o /s) [+]
!        qsubl_lake ,&! sublimation rate from snow pack (mm h2o /s) [+]
        qfros_road ,&! surface dew added to snow pack (mm h2o /s) [+]
!        qfros_lake ,&! surface dew added to snow pack (mm h2o /s) [+]
        scvold_road,&! snow mass on road for previous time step [kg/m2]
!        scvold_lake,&! snow mass on lake for previous time step [kg/m2]

        sm_road    ,&! rate of snowmelt [kg/(m2 s)]
!        sm_lake    ,&! rate of snowmelt [kg/(m2 s)]
        totwb      ,&! water mass at the begining of time step

        wt         ,&! fraction of vegetation buried (covered) by snow [-]
        sigf       ,&! fraction of veg cover, excluding snow-covered veg [-], only used to call snowfraction
        rootr(1:nl_soil),&! root resistance of a layer, all layers add to 1.0
        rootflux(1:nl_soil),&! root flux of a layer, all layers add to 1.0

        z_roadsno (maxsnl+1:nl_soil) ,&! layer depth [m]
!        z_lakesno (maxsnl+1:nl_soil) ,&! layer depth [m]
        dz_roadsno(maxsnl+1:nl_soil) ,&! layer thickness [m]
!        dz_lakesno(maxsnl+1:nl_soil) ,&! layer thickness [m]
        zi_roadsno(maxsnl  :nl_soil)    ! interface level below a "z" level [m]
!        zi_lakesno(maxsnl  :nl_soil)   ! interface level below a "z" level [m]

  real(r8) :: &
        prc_rain   ,&! convective rainfall [kg/(m2 s)]
        prc_snow   ,&! convective snowfall [kg/(m2 s)]
        prl_rain   ,&! large scale rainfall [kg/(m2 s)]
        prl_snow   ,&! large scale snowfall [kg/(m2 s)]
        t_precip   ,&! snowfall/rainfall temperature [kelvin]
        bifall     ,&! bulk density of newly fallen dry snow [kg/m3]
        pg_rain    ,&! rainfall onto ground including canopy runoff [kg/(m2 s)]
        pg_snow    ,&! snowfall onto ground including canopy runoff [kg/(m2 s)]
        pgroad_rain ,&! rainfall onto road including canopy runoff [kg/(m2 s)]
        pgroad_snow   ! snowfall onto road including canopy runoff [kg/(m2 s)]
!        pg_rain_lake,&!rainfall onto lake [kg/(m2 s)]
!        pg_snow_lake,&!snowfall onto lake [kg/(m2 s)]

  !real(r8) :: &
  !      errw_rsub    ! the possible subsurface runoff deficit after PHS is included

  !real(r8) :: &
  !      ei,         &! vapor pressure on leaf surface [pa]
  !      deidT,      &! derivative of "ei" on "tl" [pa/K]
  !      qsatl,      &! leaf specific humidity [kg/kg]
  !      qsatldT      ! derivative of "qsatl" on "tlef"

  integer :: &
        snlroad    ,&! number of snow layers on road
!        snll       ,&! number of snow layers
        imeltroad(maxsnl+1:nl_soil), &! flag for: melting=1, freezing=2, Nothing happended=0
!        imeltl(maxsnl+1:nl_soil), &! flag for: melting=1, freezing=2, Nothing happended=0
        lbroad     ,&! lower bound of arrays
!        lbl        ,&! lower bound of arrays
!        lbsn       ,&! lower bound of arrays
        data_sc(4)   ,&! /year/month/day/hour
        j            ! DO looping index

  logical :: snow_clear_flag     !true => start snow clear

!  ! For SNICAR snow model
!  !----------------------------------------------------------------------
!  real(r8) forc_aer        ( 14 )  !aerosol deposition from atmosphere model (grd,aer) [kg m-1 s-1]
!  real(r8) snofrz    (maxsnl+1:0)  !snow freezing rate (col,lyr) [kg m-2 s-1]
!  real(r8) sabg_lyr  (maxsnl+1:1)  !snow layer absorption [W/m-2]

  theta = acos(max(coszen,0.001))
!  forc_aer(:) = 0.        !aerosol deposition from atmosphere model (grd,aer) [kg m-1 s-1]


!======================================================================
! [0] Snow clear
!======================================================================

IF (MOD(idate(3), 3600) == 0) THEN
   CALL julian2monthday(idate(1), idate(2), data_sc(2), data_sc(3))
   data_sc(1) = idate(1)
   data_sc(4) = idate(3) / 3600

   CALL snow_clear_CDP(data_sc, snow_clear_flag)

   IF (snow_clear_flag) THEN
      scv_road     = 0
      snowdp_road  = 0
      fsno_road    = 0
      wliq_roadsno = 0
      wice_roadsno = 0
   ENDIF
ENDIF


!======================================================================
! [1] Solar absorbed by ground
!      and precipitation information (rain/snow fall and precip temperature)
!======================================================================
  sabroad = 0.
!  sablake = 0.
!  sabv    = 0.
!  par     = 0.
  
  IF (forc_sols+forc_soll+forc_solsd+forc_solld > 0.) THEN

    sabroad = forc_sols * sroad(1,1) + forc_soll * sroad(2,1) &
            + forc_solsd * sroad(1,2) + forc_solld * sroad(2,2)
   
!    sabv    = forc_sols * ssun (1,1) + forc_soll * ssun (2,1) &
!            + forc_solsd * ssun (1,2) + fsorc_solld * ssun (2,2)

!    par     = forc_sols * ssun (1,1) + forc_solsd * ssun (1,2)
  ENDIF

  solvd = forc_sols
  solvi = forc_solsd
  solnd = forc_soll
  solni = forc_solld
  srvd  = solvd * alb(1,1)
  srvi  = solvi * alb(1,2)
  srnd  = solnd * alb(2,1)
  srni  = solni * alb(2,2)
  sr    = srvd + srvi + srnd + srni

  ! calculate the local secs
  radpsec = pi/12./3600.
  IF ( isgreenwich ) THEN
    local_secs = idate(3) + nint((patchlonr/radpsec)/deltim)*deltim
    local_secs = mod(local_secs, 86400)
  ELSE
    local_secs = idate(3)
  ENDIF

  IF (local_secs == 86400/2) THEN
    solvdln = forc_sols
    solviln = forc_solsd
    solndln = forc_soll
    solniln = forc_solld
    srvdln  = solvdln * alb(1,1)
    srviln  = solviln * alb(1,2)
    srndln  = solndln * alb(2,1)
    srniln  = solniln * alb(2,2)
  ELSE
    solvdln = spval
    solviln = spval
    solndln = spval
    solniln = spval
    srvdln  = spval
    srviln  = spval
    srndln  = spval
    srniln  = spval
  ENDIF

  CALL rain_snow_temp (patchtype,forc_t,forc_q,forc_psrf,forc_prc,forc_prl,forc_us,forc_vs,tcrit, &
      prc_rain,prc_snow,prl_rain,prl_snow,t_precip,bifall)

  forc_rain = prc_rain + prl_rain
  forc_snow = prc_snow + prl_snow

!  sabvsun = sabv * fveg * (1 - flake)
!  sabvsha = 0.

!======================================================================

  z_roadsno (maxsnl+1:0) = z_sno_road (maxsnl+1:0)
  z_roadsno (1:nl_soil ) = z_soi (1:nl_soil)
  dz_roadsno(maxsnl+1:0) = dz_sno_road(maxsnl+1:0)
  dz_roadsno(1:nl_soil ) = dz_soi(1:nl_soil)

  !============================================================
  scvold_road = scv_road        !snow mass at previous time step

  snlroad = 0
  DO j = maxsnl+1, 0
     IF (wliq_roadsno(j)+wice_roadsno(j) > 0.) snlroad = snlroad - 1
  ENDDO

  zi_roadsno(0) = 0.
  IF (snlroad < 0) THEN
     DO j = -1, snlroad, -1
        zi_roadsno(j) = zi_roadsno(j+1) - dz_roadsno(j+1)
     ENDDO
  ENDIF

  zi_roadsno(1:nl_soil) = zi_soi(1:nl_soil)

!  totwb = scv_road + wice_roadsno(1) + wliq_roadsno(1)
  fioldroad(:) = 0.0
  IF (snlroad < 0) THEN
     fioldroad(snlroad+1:0) = wice_roadsno(snlroad+1:0) / &
        (wliq_roadsno(snlroad+1:0) + wice_roadsno(snlroad+1:0))
  ENDIF

  !============================================================
!  scvold_lake = scv_lake        !snow mass at previous time step
!
!  snll = 0
!  DO j = maxsnl+1, 0
!     IF (wliq_lakesno(j) + wice_lakesno(j) > 0.) snll = snll - 1
!  ENDDO
!
!  zi_lakesno(0) = 0.
!  IF (snll < 0) THEN
!     DO j = -1, snll, -1
!        zi_lakesno(j) = zi_lakesno(j+1) - dz_lakesno(j+1)
!     ENDDO
!  ENDIF
!
!  zi_lakesno(1:nl_soil) = zi_soi(1:nl_soil)
!
!  w_old = sum(wliq_lakesno(snll+1:))
!  fioldl(:) = 0.0
!  IF (snll <0 ) THEN
!     fioldl(snll+1:0) = wice_lakesno(snll+1:0) / &
!        (wliq_lakesno(snll+1:0) + wice_lakesno(snll+1:0))
!  ENDIF

  !============================================================
  totwb  = sum(wice_roadsno(:1) + wliq_roadsno(:1))
  totwb  = totwb + scv_road !+ ldew*fveg + wa*(1-froof)*fgper

!----------------------------------------------------------------------
! [2] Canopy interception and precipitation onto ground surface
!----------------------------------------------------------------------
! No canopy effect is considered so this part is skipped.
  pg_rain = prc_rain + prl_rain
  pg_snow = prc_snow + prl_snow
!  pg_rain_lake = prc_rain + prl_rain
!  pg_snow_lake = prc_snow + prl_snow

  pgroad_rain = pg_rain
  pgroad_snow = pg_snow
!  pgroad_rain = prc_rain + prl_rain + prc_snow + prl_snow
!  pgroad_snow = 0

!----------------------------------------------------------------------
! [3] Initilize new snow nodes for snowfall / sleet
!----------------------------------------------------------------------
  lbroad = snlroad + 1           !lower bound of array
  troad = t_roadsno(lbroad)
 
  CALL newsnow (patchtype,maxsnl,deltim,troad,pgroad_rain,pgroad_snow,bifall,      &
                t_precip,zi_roadsno(:0),z_roadsno(:0),dz_roadsno(:0),t_roadsno(:0),&
                wliq_roadsno(:0),wice_roadsno(:0),fioldroad(:0),                   &
                snlroad,sag_road,scv_road,snowdp_road,fsno_road                    )
 
!  CALL newsnow_lake ( &
!                ! "in" arguments
!                ! ---------------
!                maxsnl        ,nl_lake       ,deltim          ,dz_lake         ,&
!                pg_rain_lake  ,pg_snow_lake  ,t_precip        ,bifall          ,&
!     
!                ! "inout" arguments
!                ! ------------------
!                t_lake        ,zi_lakesno(:0),z_lakesno(:0)                    ,&
!                dz_lakesno(:0),t_lakesno(:0) ,wliq_lakesno(:0),wice_lakesno(:0),&
!                fioldl(:0)    ,snll          ,sag_lake        ,scv_lake        ,&
!                snowdp_lake   ,lake_icefrac                                     )
     
!----------------------------------------------------------------------
! [4] Energy and Water balance
!----------------------------------------------------------------------
  lbroad = snlroad + 1           !lower bound of array
!  lbl = snll + 1           !lower bound of array
  !lbsn= min(lbp,0)

  ! Thermal process
  CALL RoadThermal ( &
    ! model running information
    ipatch           ,patchtype             ,lbroad            ,&
    deltim           ,patchlatr                                                            ,&
    ! forcing
    forc_hgt_u       ,forc_hgt_t            ,forc_hgt_q        ,forc_us                    ,&
    forc_vs          ,forc_t                ,forc_q            ,forc_psrf                  ,&
    forc_rhoair      ,forc_frl              ,& !forc_po2m         ,forc_pco2m                 ,&
    forc_sols        ,forc_soll             ,forc_solsd        ,forc_solld                 ,&
    theta            ,& !sabwsun               ,sabwsha                                       ,&
    sabroad          ,& !sablake               ,sabv              ,
  !  par              ,&
    ! GROUND PARAMETERS
  !  flake            ,
    pondmx           ,trsmx0                ,zlnd                                          ,&
    zsno             ,capr                  ,cnfac             ,vf_quartz                  ,&
    vf_gravels       ,vf_om                 ,vf_sand           ,wf_gravels                 ,&
    wf_sand          ,csol                  ,porsl             ,psi0                       ,&
#ifdef Campbell_SOIL_MODEL
    bsw              ,&
#endif
#ifdef vanGenuchten_Mualem_SOIL_MODEL
    theta_r          ,alpha_vgm             ,n_vgm             ,L_vgm                      ,&
    sc_vgm           ,fc_vgm                ,&
#endif
    k_solids         ,dksatu                ,dksatf            ,dkdry                      ,&
    BA_alpha         ,BA_beta                                                              ,&
    em_road          ,cv_road               ,tk_road                                       ,&
    dz_roadsno(lbroad:)                     ,&
    z_roadsno(lbroad:)                      ,&
    zi_roadsno(lbroad-1:)                   ,&
!    dz_lake(1:)      ,lakedepth             ,&
    ! vegetation parameters
!    dewmx            ,sqrtdi                ,rootfr(:)         ,effcon                     ,&
!    vmax25           ,slti                  ,hlti              ,shti                       ,&
!    hhti             ,trda                  ,trdm              ,trop                       ,&
!    g1               ,g0                    ,gradm             ,binter                     ,&
!    extkn            ,&
    ! surface status
    fsno_road        ,scv_road              ,snowdp_road                                   ,&
!    lai                  ,&
!    sai                  ,htop                 ,hbot                 ,&
!    extkd                ,
    lroad                ,t_grnd               ,&
    troad                ,t_roadsno(lbroad:)   ,wliq_roadsno(lbroad:)                      ,&
    wice_roadsno(lbroad:),&
!    lake_icefrac(:)      ,savedtke1            ,lveg                 ,tleaf                ,&
!    ldew                 ,&
!! SNICAR model variables
!    snofrz         ,sabg_lyr                                       ,&
!! END SNICAR model variables
    ! output
    taux                 ,tauy                 ,fsena                ,fevpa                ,&
    lfevpa               ,& !fsenl             ,fevpl                ,etr                  ,&
    fseng                ,fevpg                ,olrg                 ,fgrnd                ,&
    fsen_road            ,lfevp_road           ,qseva_road           ,&
    qsdew_road           ,qsubl_road           ,qfros_road           ,&
    imeltroad(lbroad:)   ,sm_road              ,&
    sabg              ,& !rstfac               ,rootr(:)             ,
    tref                 ,&
    qref                 ,trad                 ,rst                  ,assim                ,&
    respc                ,errore               ,emis                 ,z0m                  ,&
    zol                  ,rib                  ,ustar                ,qstar                ,&
    tstar                ,fm                   ,fh                   ,fq                   ,&
    hpbl                 ,pgroad_rain          ,pgroad_snow          ,t_precip              &
  )

!----------------------------------------------------------------------
! [5] Road hydrology
!----------------------------------------------------------------------
  CALL RoadHydrology ( &
        ! model running information
        ipatch                ,patchtype     ,lbroad    ,deltim ,&
        ! forcing
      !  pg_rain               ,
        pgroad_rain           ,&             !pg_snow           ,&
        ssi                   ,wimp                             ,&
      !  fsen_road             ,fgrnd                            ,&
        dz_roadsno(lbroad:)   ,&
        wliq_roadsno(lbroad:) ,&
        wice_roadsno(lbroad:) ,&
        qseva_road            ,qsdew_road                       ,&
        qsubl_road            ,qfros_road                       ,&
        sm_road               ,forc_us       ,forc_vs           ,&
        ! output
        rsur                  ,rnof) !          ,errw_rsub          )

!============================================================
  IF (snlroad < 0) THEN
       ! Compaction rate for snow
       ! Natural compaction and metamorphosis. The compaction rate
       ! is recalculated for every new timestep
       lbroad  = snlroad + 1   ! lower bound of array
       CALL snowcompaction (lbroad, deltim                                                ,&
                        imeltroad(lbroad:0), fioldroad(lbroad:0), t_roadsno(lbroad:0)     ,&
                        wliq_roadsno(lbroad:0), wice_roadsno(lbroad:0)                    ,&
                        forc_us, forc_vs, dz_roadsno(lbroad:0)                             )
 
       ! Combine thin snow elements
       lbroad = maxsnl + 1
       CALL snowlayerscombine (lbroad,snlroad,&
                        z_roadsno(lbroad:1), dz_roadsno(lbroad:1), zi_roadsno(lbroad-1:1)  ,&
                        wliq_roadsno(lbroad:1), wice_roadsno(lbroad:1), t_roadsno(lbroad:1),&
                        scv_road, snowdp_road)
 
       ! Divide thick snow elements
       IF (snlroad < 0)  &
          CALL snowlayersdivide (lbroad, snlroad                                           ,&
                        z_roadsno(lbroad:0), dz_roadsno(lbroad:0), zi_roadsno(lbroad-1:0)  ,&
                        wliq_roadsno(lbroad:0), wice_roadsno(lbroad:0), t_roadsno(lbroad:0) )
  ENDIF
 
  ! Set zero to the empty node
  IF (snlroad > maxsnl) THEN
       wice_roadsno(maxsnl+1:snlroad) = 0.
       wliq_roadsno(maxsnl+1:snlroad) = 0.
       t_roadsno   (maxsnl+1:snlroad) = 0.
       z_roadsno   (maxsnl+1:snlroad) = 0.
       dz_roadsno  (maxsnl+1:snlroad) = 0.
  ENDIF
 
  lbroad = snlroad + 1
  troad = t_roadsno(lbroad)
 
  !TODO: temporal, set to t_soisno
  !t_soisno(:) = t_gpersno(:)

  ! ----------------------------------------
  ! energy balance check
  ! ----------------------------------------
  zerr=errore
#if(defined CoLMDEBUG)
  IF(abs(errore)>.5)THEN
     write(6,*) 'Warning: energy balance violation in Road ',errore,patchclass
  ENDIF
#endif
  
! ----------------------------------------
! water balance check
! ----------------------------------------
! LIXX TODO: Check these statements  
!  wliq_roadsno(:) = 0.
!  wliq_roadsno(:1) = wliq_roadsno(:1) + wliq_roadsno(:1)
  !wliq_soisno(:) = wliq_soisno(:)*(1-flake) + wliq_lakesno(:)*flake
  
!  wice_roadsno(:) = 0.
!  wice_roadsno(:1) = wice_roadsno(:1) + wice_roadsno(:1)
  !wice_soisno(:) = wice_soisno(:)*(1-flake) + wice_lakesno(:)*flake
  
  scv = scv_road
  !scv = scv*(1-flake) + scv_lake*flake
  
  endwb  = sum(wice_roadsno(:1) + wliq_roadsno(:1))
  endwb  = endwb + scv !+ ldew*fveg
  errorw = (endwb - totwb) - (forc_prc + forc_prl - fevpa - rnof)*deltim !- errw_rsub)*deltim
  xerr   = errorw/deltim
  
#if(defined CoLMDEBUG)
  IF(abs(errorw)>1.e-3) THEN
      write(6,*) 'Warning: water balance violation in Road ', errorw, ipatch, patchclass
      write(6,*) endwb, totwb, forc_prc, forc_prl, fevpa, rnof, deltim
      !STOP
  ENDIF
  
#endif
  
!======================================================================
! Preparation for the next time step
! 1) time-varying parameters for vegatation
! 2) fraction of snow cover
! 3) solar zenith angle and
! 4) albedos
!======================================================================
  
  ! cosine of solar zenith angle
  calday = calendarday(idate)
  coszen = orb_coszen(calday,patchlonr,patchlatr)

  ! fraction of snow cover.
  ! TODO: this is for vegetation only? Or should I rewrite it?
  CALL snowfraction ( 0., 0.,z0m,zlnd,scv_road,snowdp_road,wt,sigf,fsno_road)

!  lai = tlai(ipatch)
!  sai = tsai(ipatch) * sigf

  ! update the snow age
  !TODO: can be moved to UrbanALBEDO.F90
  IF (snlroad == 0) sag_road = 0.
  CALL snowage (deltim,troad,scv_road,scvold_road,sag_road)
 
  ! update snow depth, snow cover and snow age
  snowdp = snowdp_road
  fsno   = fsno_road
  sag    = sag_road

  ! albedos
  ! we supposed call it every time-step, because
  ! other vegeation related parameters are needed to create

  CALL albroad (ipatch,alb_road,coszen,fsno_road,&
                scv_road,sag_road,alb,sroad)

  ! zero-filling set for glacier/ice-sheet/land water bodies/ocean components
!  laisun = lai
!  laisha = 0.0
!  green  = 1.

  h2osoi = wliq_roadsno(1:)/(dz_soi(1:)*denh2o) + wice_roadsno(1:)/(dz_soi(1:)*denice)
  wat = sum(wice_roadsno(1:)+wliq_roadsno(1:))
  wat = wat + scv_road !+ ldew*fveg

  z_sno_road (maxsnl+1:0) = z_roadsno (maxsnl+1:0)
  dz_sno_road(maxsnl+1:0) = dz_roadsno(maxsnl+1:0)

  z_sno(:) = z_sno_road(:)
  dz_sno(:) = dz_sno_road(:)

! diagnostic diurnal temperature
  !IF (tref > tmax) tmax = tref
  !IF (tref < tmin) tmin = tref

! 06/05/2022, yuan: RH for output to compare
!  CALL qsadv(tref,forc_psrf,ei,deiDT,qsatl,qsatlDT)
!  qref = qref/qsatl

END SUBROUTINE CoLMMain_Road
! ----------------------------------------------------------------------
! EOP