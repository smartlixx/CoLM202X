#include <define.h>

#if (defined ROAD_MODEL)
MODULE MOD_Road_Vars_TimeVariables

!-----------------------------------------------------------------------
! !DESCRIPTION:
!
!  Define urban model time variant variables.
!
!  Created by Hua Yuan, 12/2020
!-----------------------------------------------------------------------

   USE MOD_Precision
   IMPLICIT NONE
   SAVE
! -----------------------------------------------------------------
! Time-varying state variables which reaquired by restart run

   ! shortwave absorption
   real(r8), allocatable :: sroad      (:,:,:) !road absorptioin [-]

   ! net longwave radiation for last time temperature change
   real(r8), allocatable :: lroad          (:) !net longwave of road  [W/m2]

   real(r8), allocatable :: z_sno_road   (:,:) !node depth of road [m]
   
   real(r8), allocatable :: dz_sno_road  (:,:) !interface depth of road [m]

   real(r8), allocatable :: t_roadsno    (:,:) !temperature of road [K]

   real(r8), allocatable :: wliq_roadsno (:,:) !liquid water in layers of road [kg/m2]
   real(r8), allocatable :: wice_roadsno (:,:) !ice lens in layers of road [kg/m2]

   real(r8), allocatable :: sag_road       (:) !road snow age [-]

   real(r8), allocatable :: scv_road       (:) !road snow mass [kg/m2]

   real(r8), allocatable :: fsno_road      (:) !road snow fraction [-]

   real(r8), allocatable :: snowdp_road    (:) !road snow depth [m]


! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: allocate_RoadTimeVariables
   PUBLIC :: deallocate_RoadTimeVariables
   PUBLIC :: READ_RoadTimeVariables
   PUBLIC :: WRITE_RoadTimeVariables

! PRIVATE MEMBER FUNCTIONS:

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE allocate_RoadTimeVariables ()
! ------------------------------------------------------
! Allocates memory for CoLM [numroad] variables
! ------------------------------------------------------
   USE MOD_Precision
   USE MOD_SPMD_Task
   USE MOD_LandRoad
   USE MOD_Vars_Global
   IMPLICIT NONE

      IF (p_is_worker) THEN
         IF (numroad > 0) THEN
            allocate (sroad                     (2,2,numroad))
       
            allocate (lroad                         (numroad))

            allocate (z_sno_road         (maxsnl+1:0,numroad))

            allocate (dz_sno_road        (maxsnl+1:0,numroad))

            allocate (t_roadsno    (maxsnl+1:nl_soil,numroad))

            allocate (wliq_roadsno (maxsnl+1:nl_soil,numroad))
            allocate (wice_roadsno (maxsnl+1:nl_soil,numroad))

            allocate (sag_road                      (numroad))
            allocate (scv_road                      (numroad))

            allocate (fsno_road                     (numroad))
            allocate (snowdp_road                   (numroad))
         ENDIF
      ENDIF
   END SUBROUTINE allocate_RoadTimeVariables

   SUBROUTINE READ_RoadTimeVariables (file_restart)

   USE MOD_NetCDFVector
   USE MOD_LandRoad
   USE MOD_Vars_Global

   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart

      CALL ncio_read_vector (file_restart, 'sroad', 2, 2, landroad, sroad)

      CALL ncio_read_vector (file_restart, 'lroad', landroad, lroad)

      CALL ncio_read_vector (file_restart, 'z_sno_road' , -maxsnl, landroad, z_sno_road )

      CALL ncio_read_vector (file_restart, 'dz_sno_road', -maxsnl, landroad, dz_sno_road)

      CALL ncio_read_vector (file_restart, 't_roadsno', nl_soil-maxsnl, landroad, t_roadsno)

      CALL ncio_read_vector (file_restart, 'wliq_roadsno', nl_soil-maxsnl, landroad, wliq_roadsno)
      CALL ncio_read_vector (file_restart, 'wice_roadsno', nl_soil-maxsnl, landroad, wice_roadsno)

      CALL ncio_read_vector (file_restart, 'sag_road'   , landroad, sag_road   )
      CALL ncio_read_vector (file_restart, 'scv_road'   , landroad, scv_road   )
      CALL ncio_read_vector (file_restart, 'fsno_road'  , landroad, fsno_road  )
      CALL ncio_read_vector (file_restart, 'snowdp_road', landroad, snowdp_road)

   END SUBROUTINE READ_RoadTimeVariables

   SUBROUTINE WRITE_RoadTimeVariables (file_restart)

   USE MOD_Namelist, only : DEF_REST_CompressLevel
   USE MOD_LandRoad
   USE MOD_NetCDFVector
   USE MOD_Vars_Global
   IMPLICIT NONE

   character(len=*), intent(in) :: file_restart

   ! Local variables
   integer :: compress

      compress = DEF_REST_CompressLevel

      CALL ncio_create_file_vector (file_restart, landroad)
      CALL ncio_define_dimension_vector (file_restart, landroad, 'road')

      CALL ncio_define_dimension_vector (file_restart, landroad, 'snow'    , -maxsnl       )
      CALL ncio_define_dimension_vector (file_restart, landroad, 'road'    , nl_soil       )

      CALL ncio_define_dimension_vector (file_restart, landroad, 'roadsnow', nl_soil-maxsnl)

      CALL ncio_define_dimension_vector (file_restart, landroad, 'band', 2)
      CALL ncio_define_dimension_vector (file_restart, landroad, 'rtyp', 2)

      CALL ncio_write_vector (file_restart, 'sroad', 'band', 2, 'rtyp', 2, 'road', landroad, sroad, compress)

      CALL ncio_write_vector (file_restart, 'lroad', 'road', landroad, lroad, compress)

      CALL ncio_write_vector (file_restart, 'z_sno_road' , 'snow', -maxsnl, 'road', landroad, z_sno_road, compress)

      CALL ncio_write_vector (file_restart, 'dz_sno_road', 'snow', -maxsnl, 'road', landroad, dz_sno_road, compress)

      CALL ncio_write_vector (file_restart, 't_roadsno', 'roadsnow', nl_soil-maxsnl, 'road', landroad, t_roadsno, compress)

      CALL ncio_write_vector (file_restart, 'wliq_roadsno', 'roadsnow', nl_soil-maxsnl, 'road', landroad, wliq_roadsno, compress)
      CALL ncio_write_vector (file_restart, 'wice_roadsno', 'roadsnow', nl_soil-maxsnl, 'road', landroad, wice_roadsno, compress)

      CALL ncio_write_vector (file_restart, 'sag_road'   , 'road', landroad, sag_road   , compress)
      CALL ncio_write_vector (file_restart, 'scv_road'   , 'road', landroad, scv_road   , compress)
      CALL ncio_write_vector (file_restart, 'fsno_road'  , 'road', landroad, fsno_road  , compress)
      CALL ncio_write_vector (file_restart, 'snowdp_road', 'road', landroad, snowdp_road, compress)

   END SUBROUTINE WRITE_RoadTimeVariables

   SUBROUTINE deallocate_RoadTimeVariables

   USE MOD_SPMD_Task
   USE MOD_LandRoad

      IF (p_is_worker) THEN
         IF (numroad > 0) THEN
            deallocate (sroad        )

            deallocate (lroad        )

            deallocate (z_sno_road   )
 
            deallocate (dz_sno_road  )

            deallocate (t_roadsno    )

            deallocate (wliq_roadsno )
            deallocate (wice_roadsno )

            deallocate (sag_road     )
            deallocate (scv_road     )
            deallocate (fsno_road    )
            deallocate (snowdp_road  )
         ENDIF
      ENDIF

   END SUBROUTINE deallocate_RoadTimeVariables

END MODULE MOD_Road_Vars_TimeVariables
! ---------- EOP ------------
#endif
