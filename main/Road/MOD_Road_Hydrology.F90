#include <define.h>

MODULE MOD_Road_Hydrology

   USE MOD_Precision
   IMPLICIT NONE
   SAVE

   PUBLIC :: RoadHydrology

CONTAINS

   SUBROUTINE RoadHydrology ( &
        ! model running information
        ipatch         ,patchtype      ,lbroad         ,deltim         ,&
        ! forcing
        pg_rain        ,pgroad_rain    ,pg_snow                        ,&
        ! surface parameters or status
        ssi            ,wimp           ,&
        fseng          ,fgrnd          ,&
        dz_roadsno     ,wliq_roadsno   ,wice_roadsno                   ,&
        qseva_road     ,qsdew_road     ,qsubl_road     ,qfros_road     ,&
        sm_road        ,forc_us        ,forc_vs                        ,&
        ! output
        rsur           ,rnof           ,errw_rsub                      ,&
      )

!=======================================================================
! this is the main SUBROUTINE to execute the calculation of URBAN
! hydrological processes
!
!=======================================================================

   USE MOD_Precision
   USE MOD_Vars_Global
   USE MOD_SoilSnowHydrology

   IMPLICIT NONE

!-----------------------Argument----------------------------------------
   integer, intent(in) :: &
        ipatch           ,&! patch index
        patchtype        ,&! land patch type (0=soil, 1=urban or built-up, 2=wetland,
                           ! 3=land ice, 4=land water bodies, 99=ocean
        lbroad             ! lower bound of array

   real(r8), intent(in) :: &
        deltim           ,&! time step (s)
        pg_rain          ,&! rainfall after removal of interception (mm h2o/s)
        pg_snow          ,&! snowfall after removal of interception (mm h2o/s)
        pgroad_rain      ,&! rainfall after removal of interception (mm h2o/s)
        ssi              ,&! irreducible water saturation of snow
        wimp             ,&! water impremeable IF porosity less than wimp
  
        qseva_road       ,&! ground surface evaporation rate (mm h2o/s)
        qsdew_road       ,&! ground surface dew formation (mm h2o /s) [+]
        qsubl_road       ,&! sublimation rate from snow pack (mm h2o /s) [+]
        qfros_road       ,&! surface dew added to snow pack (mm h2o /s) [+]
        sm_road            ! snow melt (mm h2o/s)

   real(r8), intent(in) :: forc_us
   real(r8), intent(in) :: forc_vs

   real(r8), intent(inout) :: &
        dz_roadsno  (lbroad:nl_soil)  ,&! layer thickness (m)
        wliq_roadsno(lbroad:nl_soil)  ,&! liquid water (kg/m2)
        wice_roadsno(lbroad:nl_soil)  ,&! ice lens (kg/m2)
        fseng            ,&! sensible heat from ground
        fgrnd              ! ground heat flux

   real(r8), intent(out) :: &
        rsur             ,&! surface runoff (mm h2o/s)
        rnof               ! total runoff (mm h2o/s)

   real(r8), intent(out) :: &
        errw_rsub          ! the possible subsurface runoff deficit after PHS is included
!
!-----------------------Local Variables------------------------------
!
   real(r8) :: &
        gwat             ,&! net water input from top (mm/s)
        rnof_road        ,&! total runoff (mm h2o/s)
        rsur_road          ! surface runoff (mm h2o/s)

   real(r8) :: a, aa, xs1

!=======================================================================
! [1] for impervious road
!=======================================================================
      IF (lbroad >= 1) THEN
         gwat = pgroad_rain + sm_road - qseva_road
      ELSE
         CALL snowwater (lbroad,deltim,ssi,wimp,&
                         pgroad_rain,qseva_road,qsdew_road,qsubl_road,qfros_road,&
                         dz_roadsno(lbroad:0),wice_roadsno(lbroad:0),wliq_roadsno(lbroad:0),&
                         gwat)
      ENDIF

      wliq_roadsno(1) = wliq_roadsno(1) + gwat*deltim

      ! Renew the ice and liquid mass due to condensation
      IF (lbroad >= 1) THEN
         ! make consistent with how evap_grnd removed in infiltration
         wliq_roadsno(1) = max(0., wliq_roadsno(1) + qsdew_road * deltim)
         wice_roadsno(1) = max(0., wice_roadsno(1) + (qfros_road-qsubl_road) * deltim)
      ENDIF

      ! only consider ponding and surface runoff
      ! NOTE: set max ponding depth = 1mm
      xs1 = wliq_roadsno(1) - 1.
      IF (xs1 > 0.) THEN
         wliq_roadsno(1) = 1.
      ELSE
         xs1 = 0.
      ENDIF

      rsur_road = xs1 / deltim
      rnof_road = rsur_road


!=======================================================================
! [2] surface and total runoff weighted by fractional coverages
!=======================================================================

      rsur = rsur_road
      rnof = rnof_road

   END SUBROUTINE RoadHydrology

END MODULE MOD_Road_Hydrology
! ---------- EOP ------------
