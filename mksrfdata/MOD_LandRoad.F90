#include <define.h>

MODULE MOD_LandRoad

!--------------------------------------------------------------------------------------
! DESCRIPTION:
!
!    Build pixelset "landroad".
!
! Original authors: Xianxiang Li, 07/2024
!
!--------------------------------------------------------------------------------------

   USE MOD_Grid
   USE MOD_Pixelset
   USE MOD_Vars_Global, only: URBAN
#ifdef SinglePoint
   USE MOD_SingleSrfdata
#endif

   IMPLICIT NONE

   ! ---- Instance ----
   type(grid_type) :: groad

   integer :: numroad
   type(pixelset_type) :: landroad

   integer , allocatable :: road_reg   (:)  !region index of a road
   integer , allocatable :: road2patch (:)  !patch index of a road
   integer , allocatable :: patch2road (:)  !road index of a patch

   ! ---- PUBLIC routines ----
   PUBLIC :: landroad_build
   PUBLIC :: map_patch_to_road

CONTAINS

   ! -------------------------------
   SUBROUTINE landroad_build (lc_year)

   USE MOD_Precision
   USE MOD_Vars_Global
   USE MOD_SPMD_Task
   USE MOD_NetCDFBlock
   USE MOD_Grid
   USE MOD_DataType
   USE MOD_Namelist
   USE MOD_5x5DataReadin
   USE MOD_Mesh
   USE MOD_LandPatch
   USE MOD_LandElm
#ifdef CATCHMENT
   USE MOD_LandHRU
#endif
   USE MOD_AggregationRequestData
   USE MOD_Utils

   IMPLICIT NONE

   integer, intent(in) :: lc_year
   ! Local Variables
   character(len=256) :: dir_road
   type (block_data_int32_2d) :: data_urb_class ! urban type index

   ! index
   integer :: nr_glb, npatch_glb

   character(len=256) :: suffix, cyear ! delete unneeded index and vars

      IF (p_is_master) THEN
         write(*,'(A)') 'Making road type tiles :'
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      ! allocate and read the grided LCZ/NCAR urban type
      IF (p_is_io) THEN

         dir_road = trim(DEF_dir_rawdata) // '/urban_type'

         CALL allocate_block_data (groad, data_urb_class)
         CALL flush_block_data (data_urb_class, 0)

         !read LCZ data
         suffix = 'URBTYP'
         CALL read_5x5_data (dir_road, suffix, groad, 'LCZ_DOM', data_urb_class)

#ifdef USEMPI
         CALL aggregation_data_daemon (groad, data_i4_2d_in1 = data_urb_class)
#endif
      ENDIF

      IF (p_is_worker) THEN

         ! delete urban types and update road patch number
         IF (numpatch > 0) THEN
            numroad = count(landpatch%settyp == URBAN)
         ELSE
            numroad = 0
         ENDIF

         IF (numroad > 0) THEN
            allocate (landroad%eindex (numroad))
            allocate (landroad%settyp (numroad))
            allocate (landroad%ipxstt (numroad))
            allocate (landroad%ipxend (numroad))
            allocate (landroad%ielm   (numroad))

            ! copy urban path information from landpatch for landroad
            landroad%eindex = pack(landpatch%eindex, landpatch%settyp == URBAN)
            landroad%ipxstt = pack(landpatch%ipxstt, landpatch%settyp == URBAN)
            landroad%ipxend = pack(landpatch%ipxend, landpatch%settyp == URBAN)
            landroad%ielm   = pack(landpatch%ielm  , landpatch%settyp == URBAN)
            landroad%settyp = URBAN
         ENDIF

         ! update land patch with roadr type patch
         ! set numroad
         landroad%nset = numroad
         landpatch%nset = numpatch
      ENDIF

      CALL landpatch%set_vecgs
      CALL landroad%set_vecgs

      CALL map_patch_to_road

#ifdef USEMPI
      IF (p_is_worker) THEN
         CALL mpi_reduce (numroad, nr_glb, 1, MPI_INTEGER, MPI_SUM, p_root, p_comm_worker, p_err)
         IF (p_iam_worker == 0) THEN
            write(*,'(A,I12,A)') 'Total: ', nr_glb, ' road tiles.'
         ENDIF
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      write(*,'(A,I12,A)') 'Total: ', numroad, ' road tiles.'
#endif

#ifdef SinglePoint

      ! delete urban vars
      allocate  ( SITE_em_road   (numroad) )

      allocate  ( SITE_cv_road   (nl_soil) ) ! or nl_road
      allocate  ( SITE_tk_road   (nl_soil) )

      allocate  ( SITE_alb_road  (2,2)     )

#endif

#ifndef CROP
#ifdef USEMPI
      IF (p_is_worker) THEN
         CALL mpi_reduce (numpatch, npatch_glb, 1, MPI_INTEGER, MPI_SUM, p_root, p_comm_worker, p_err)
         IF (p_iam_worker == 0) THEN
            write(*,'(A,I12,A)') 'Total: ', npatch_glb, ' patches.'
         ENDIF
      ENDIF

      CALL mpi_barrier (p_comm_glb, p_err)
#else
      write(*,'(A,I12,A)') 'Total: ', numpatch, ' patches.'
#endif

      CALL elm_patch%build (landelm, landpatch, use_frac = .true.)
#ifdef CATCHMENT
      CALL hru_patch%build (landhru, landpatch, use_frac = .true.)
#endif
      CALL write_patchfrac (DEF_dir_landdata, lc_year)
#endif

   END SUBROUTINE landroad_build

   ! ----------------------
   SUBROUTINE map_patch_to_road

   USE MOD_SPMD_Task
   USE MOD_LandPatch
   IMPLICIT NONE

   integer :: ipatch, iroad

      IF (p_is_worker) THEN

         IF ((numpatch <= 0) .or. (numroad <= 0)) RETURN

         IF (allocated(patch2road)) deallocate(patch2road)
         IF (allocated(road2patch)) deallocate(road2patch)
         allocate (patch2road (numpatch))
         allocate (road2patch (numroad))

         iroad = 0
         DO ipatch = 1, numpatch
            IF (landpatch%settyp(ipatch) == URBAN) THEN
               iroad = iroad + 1
               patch2road(ipatch) = iroad
               road2patch(iroad) = ipatch
            ELSE
               patch2road(ipatch) = -1
            ENDIF
         ENDDO

      ENDIF

   END SUBROUTINE map_patch_to_road

END MODULE MOD_LandRoad
! ---------- EOP ------------
