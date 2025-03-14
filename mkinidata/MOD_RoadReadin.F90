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

   integer ::  u, ns, nr, ulev

      write(cyear,'(i4.4)') lc_year

      IF (p_is_worker) THEN

         DO u = 1, numroad

         ! temporary setup
            
            DO ns = 1,2
               DO nr = 1,2
                  alb_road(ns,nr,u) = albroad_apt(ns,nr)
               ENDDO
            ENDDO

            em_road(u)      = emroad_apt(1,1)
            
            DO ulev = 1, nl_road
               cv_road(ulev,u)    = cvroad_apt(ulev,1)
               tk_road(ulev,u)    = tkroad_apt(ulev,1)
            ENDDO

         ENDDO
      ENDIF

   END SUBROUTINE Road_readin

END MODULE MOD_RoadReadin

#endif
