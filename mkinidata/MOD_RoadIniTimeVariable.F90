#include <define.h>

#ifdef ROAD_MODEL
MODULE MOD_RoadIniTimeVariable

   USE MOD_Precision
   IMPLICIT NONE
   SAVE

   PUBLIC :: RoadIniTimeVar

CONTAINS

   SUBROUTINE RoadIniTimeVar(ipatch,alb_road,coszen,fsno_road,&
                             scv_road,sag_road,alb,sroad)

   USE MOD_Precision
   USE MOD_Vars_Global
   USE MOD_Road_Albedo

   IMPLICIT NONE

   integer, intent(in) :: &
         ipatch          ! patch index

   real(r8), intent(in) :: &
         alb_road(2,2)   ! road albedo (iband,direct/diffuse)

   real(r8), intent(in) :: &
         coszen          ! cosine of solar zenith angle


   real(r8), intent(out) :: &
         fsno_road,     &! fraction of soil covered by snow [-]
         scv_road,      &! snow cover, water equivalent [mm]
         sag_road        ! non dimensional snow age [-]

   real(r8), intent(out) :: &
         alb (2,2),     &! averaged albedo [-]
         sroad(2,2)      ! road absorption for solar radiation,

   !-----------------------------------------------------------------------
   real(r8) :: hveg      ! height of crown central hight

      fsno_road   = 0.   ! fraction of ground covered by snow
      scv_road    = 0.   ! snow cover, water equivalent [mm, kg/m2]
      sag_road    = 0.   ! road snow age [-]

      ! urban surface albedo
      CALL albroad (ipatch,alb_road,max(0.01,coszen),fsno_road,&
                    scv_road,sag_road,alb,sroad)

   END SUBROUTINE RoadIniTimeVar

END MODULE MOD_RoadIniTimeVariable
!-----------------------------------------------------------------------
! EOP
#endif
