#include <define.h>

#ifdef ROAD_MODEL
MODULE MOD_Road_Vars_TimeInvariants

! -------------------------------
! Created by Xianxiang Li, 07/2024
! -------------------------------

   USE MOD_Precision
   IMPLICIT NONE
   SAVE

   !integer, parameter :: ns    = 2
   !integer, parameter :: nr    = 2
   !integer, parameter :: ulev  = 10
   !integer, parameter :: ityp  = 3
   !integer, parameter :: ihour = 24
   !integer, parameter :: iweek = 7
   !integer, parameter :: iday  = 365

   !integer , allocatable :: roadclass    (:)  !road type
   !integer , allocatable :: patch2road   (:)  !projection from patch to road
   !integer , allocatable :: road2patch   (:)  !projection from road to patch

   ! albedo
   real(r8), allocatable :: alb_road(:,:,:)  !albedo of road [-]

   ! emissivity
   real(r8), allocatable :: em_road     (:)  !emissivity of road [-]

   ! thermal pars 
   real(r8), allocatable :: cv_road   (:,:)  !heat capacity of road [J/(m2 K)]

   real(r8), allocatable :: tk_road   (:,:)  !thermal conductivity of road [W/m-K]

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: allocate_RoadTimeInvariants
   PUBLIC :: deallocate_RoadTimeInvariants
   !PUBLIC :: READ_RoadTimeInvariants
   !PUBLIC :: WRITE_RoadTimeInvariants

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE allocate_RoadTimeInvariants ()
! ------------------------------------------------------
! Allocates memory for CoLM 1d [numroad] variants
! ------------------------------------------------------
   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_LandRoad
   USE MOD_Vars_Global
   IMPLICIT NONE

      IF (p_is_worker) THEN
         IF (numroad > 0) THEN

            allocate (alb_road          (2,2,numroad))
            
            allocate (em_road               (numroad))

            allocate (cv_road     (1:nl_road,numroad))

            allocate (tk_road     (1:nl_road,numroad))

         ENDIF
      ENDIF

   END SUBROUTINE allocate_RoadTimeInvariants

   SUBROUTINE READ_RoadTimeInvariants (file_restart)

   USE MOD_NetCDFVector
   USE MOD_LandRoad

   IMPLICIT NONE

   integer, parameter :: ns = 2
   integer, parameter :: nr = 2
   integer, parameter :: ulev = 10
   character(len=*), intent(in) :: file_restart

      ! morphological paras
 
      CALL ncio_read_vector (file_restart, 'EM_ROAD'       , landroad, em_road  )

      ! thermal paras
      CALL ncio_read_vector (file_restart, 'CV_ROAD'   , ulev, landroad, cv_road)

      CALL ncio_read_vector (file_restart, 'TK_ROAD'   , ulev, landroad, tk_road)

      CALL ncio_read_vector (file_restart, 'ALB_ROAD'   , ns, nr, landroad, alb_road  )

   END SUBROUTINE READ_RoadTimeInvariants

   SUBROUTINE WRITE_RoadTimeInvariants (file_restart)

   USE MOD_NetCDFVector
   USE MOD_LandRoad
   USE MOD_Namelist
   USE MOD_Vars_Global

   IMPLICIT NONE

   integer, parameter :: ns    = 2
   integer, parameter :: nr    = 2
   integer, parameter :: ulev  = 10
   integer, parameter :: ityp  = 3
   integer, parameter :: ihour = 24
   integer, parameter :: iweek = 7
   integer, parameter :: iday  = 365
   ! Local variables
   character(len=*), intent(in) :: file_restart
   integer :: compress

      compress = DEF_REST_CompressLevel

      CALL ncio_create_file_vector (file_restart, landroad)
      CALL ncio_define_dimension_vector (file_restart, landroad, 'urban')

      CALL ncio_define_dimension_vector (file_restart, landroad, 'urban')
      CALL ncio_define_dimension_vector (file_restart, landroad, 'numsolar', nr  )
      CALL ncio_define_dimension_vector (file_restart, landroad, 'numrad'  , ns  )
      CALL ncio_define_dimension_vector (file_restart, landroad, 'ulev'    , ulev)
      CALL ncio_define_dimension_vector (file_restart, landroad, 'ityp'    , 3   )
      CALL ncio_define_dimension_vector (file_restart, landroad, 'iweek'   , 7   )
      CALL ncio_define_dimension_vector (file_restart, landroad, 'ihour'   , 24  )
      CALL ncio_define_dimension_vector (file_restart, landroad, 'iday'    , 365 )

      CALL ncio_write_vector (file_restart, 'EM_ROAD'       , 'urban', landroad, em_road  , DEF_REST_CompressLevel)
 
     ! thermal paras
      CALL ncio_write_vector (file_restart, 'CV_ROAD'   , 'ulev', ulev, 'urban', landroad, cv_road, DEF_REST_CompressLevel)
      CALL ncio_write_vector (file_restart, 'TK_ROAD'   , 'ulev', ulev, 'urban', landroad, tk_road, DEF_REST_CompressLevel)
     
      CALL ncio_write_vector (file_restart, 'ALB_ROAD'   , 'numsolar', ns, 'numrad', nr, 'urban', landroad, alb_road, DEF_REST_CompressLevel)
 
   END SUBROUTINE WRITE_RoadTimeInvariants

   SUBROUTINE deallocate_RoadTimeInvariants

   USE MOD_SPMD_Task
   USE MOD_LandRoad

      ! deallocate (urbclass  )

      IF (p_is_worker) THEN
         IF (numroad > 0) THEN

            deallocate (alb_road  )

            deallocate (em_road   )

            deallocate (cv_road   )
 
            deallocate (tk_road   )

         ENDIF
      ENDIF
   END SUBROUTINE deallocate_RoadTimeInvariants

END MODULE MOD_Road_Vars_TimeInvariants
#endif
! ---------- EOP ------------
