#include <define.h>

#ifdef ROAD_MODEL
MODULE MOD_Road_Vars_TimeInvariants

! -------------------------------
! Created by Xianxiang Li, 07/2024
! -------------------------------

   USE MOD_Precision
   IMPLICIT NONE
   SAVE

   integer, parameter :: ns    = 2
   integer, parameter :: nr    = 2
   integer, parameter :: ulev  = 10
   integer, parameter :: ityp  = 3
   integer, parameter :: ihour = 24
   integer, parameter :: iweek = 7
   integer, parameter :: iday  = 365

   !integer , allocatable :: roadclass    (:)  !urban type
   !integer , allocatable :: patch2road   (:)  !projection from patch to road
   !integer , allocatable :: road2patch   (:)  !projection from road to patch

   ! Road depth
   real(r8), allocatable :: z_road    (:,:)  !depth of each road layer [m]
   real(r8), allocatable :: dz_road   (:,:)  !thickness of each road layer [m]

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
   PUBLIC :: READ_RoadTimeInvariants
   PUBLIC :: WRITE_RoadTimeInvariants

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE allocate_RoadTimeInvariants ()
! ------------------------------------------------------
! Allocates memory for CoLM 1d [numurban] variants
! ------------------------------------------------------
   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_LandRoad
   USE MOD_Vars_Global
   IMPLICIT NONE

      IF (p_is_worker) THEN
         IF (numurban > 0) THEN

            allocate (alb_roof         (2,2,numurban))

            allocate (em_roof              (numurban))

            allocate (z_roof     (1:nl_roof,numurban))
            allocate (dz_roof    (1:nl_roof,numurban))

            allocate (cv_roof    (1:nl_roof,numurban))

            allocate (tk_roof    (1:nl_roof,numurban))

         ENDIF
      ENDIF

   END SUBROUTINE allocate_RoadTimeInvariants

   SUBROUTINE READ_RoadTimeInvariants (file_restart)

   USE MOD_NetCDFVector
   USE MOD_LandRoad

   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart

      ! morphological paras
 
      CALL ncio_read_vector (file_restart, 'EM_ROOF'       , landurban, em_roof  )

      CALL ncio_read_vector (file_restart, 'ROOF_DEPTH_L'  , ulev, landurban, z_roof   )
      CALL ncio_read_vector (file_restart, 'ROOF_THICK_L'  , ulev, landurban, dz_roof  )

      ! thermal paras
      CALL ncio_read_vector (file_restart, 'CV_ROOF'   , ulev, landurban, cv_roof)

      CALL ncio_read_vector (file_restart, 'TK_ROOF'   , ulev, landurban, tk_roof)

      CALL ncio_read_vector (file_restart, 'ALB_ROOF'   , ns, nr, landurban, alb_roof  )

   END SUBROUTINE READ_RoadTimeInvariants

   SUBROUTINE WRITE_RoadTimeInvariants (file_restart)

   USE MOD_NetCDFVector
   USE MOD_LandRoad
   USE MOD_Namelist
   USE MOD_Vars_Global

   IMPLICIT NONE

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

      CALL ncio_write_vector (file_restart, 'EM_ROAD'       , 'urban', landurban, em_road  , DEF_REST_CompressLevel)
 
      CALL ncio_write_vector (file_restart, 'ROAD_DEPTH_L', 'ulev', ulev, 'urban', landroad, z_road , DEF_REST_CompressLevel)
      CALL ncio_write_vector (file_restart, 'ROAD_THICK_L', 'ulev', ulev, 'urban', landroad, dz_road, DEF_REST_CompressLevel)
     ! thermal paras
      CALL ncio_write_vector (file_restart, 'CV_ROAD'   , 'ulev', ulev, 'urban', landroad, cv_road, DEF_REST_CompressLevel)
      CALL ncio_write_vector (file_restart, 'TK_ROAD'   , 'ulev', ulev, 'urban', landroad, tk_road, DEF_REST_CompressLevel)
     
      CALL ncio_write_vector (file_restart, 'ALB_ROAD'   , 'numsolar', ns, 'numrad', nr, 'urban', landroad, alb_road, DEF_REST_CompressLevel)
 
   END SUBROUTINE WRITE_RoadTimeInvariants

   SUBROUTINE deallocate_RoadTimeInvariants

   USE MOD_SPMD_Task
   USE MOD_LandUrban

      ! deallocate (urbclass  )

      IF (p_is_worker) THEN
         IF (numurban > 0) THEN

            deallocate (alb_road  )

            deallocate (em_road   )

            deallocate (z_road    )
            deallocate (dz_road   )

            deallocate (cv_road   )
 
            deallocate (tk_road   )

         ENDIF
      ENDIF
   END SUBROUTINE deallocate_RoadTimeInvariants

END MODULE MOD_Road_Vars_TimeInvariants
#endif
! ---------- EOP ------------
