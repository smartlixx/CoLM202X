#include <define.h>

#ifdef ROAD_MODEL

MODULE MOD_RoadReadin

   USE MOD_Precision
   IMPLICIT NONE
   SAVE

   ! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: Road_readin

CONTAINS

   SUBROUTINE Road_readin (dir_landdata, lc_year)

   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_Vars_Global
   USE MOD_Namelist
   USE MOD_Const_LC
   USE MOD_Vars_TimeVariables
   USE MOD_Vars_TimeInvariants
   USE MOD_Road_Vars_TimeInvariants
   USE MOD_NetCDFVector
   USE MOD_NetCDFSerial
   USE MOD_LandPatch
   USE MOD_LandRoad
   USE MOD_Road_Const_ThermalParameters
#ifdef SinglePoint
   USE MOD_SingleSrfdata
#endif

   IMPLICIT NONE

   integer, intent(in) :: lc_year    ! which year of land cover data used
   character(len=256), intent(in) :: dir_landdata
   character(len=256) :: dir_rawdata, dir_runtime
   character(len=256) :: lndname
   character(len=256) :: cyear

   integer ::  u

      write(cyear,'(i4.4)') lc_year

      IF (p_is_worker) THEN

         DO u = 1, numroad

         ! temporary setup

            alb_road(:,:,u) = albroad_apt

            em_road(u)      = emroad_apt(1,1)
            
            cv_road(:,u)    = cvroad_apt(:,1)
            tk_road(:,u)    = tkroad_apt(:,1)

ENDIF

         ENDDO
      ENDIF

   END SUBROUTINE Road_readin

END MODULE MOD_RoadReadin

#endif
