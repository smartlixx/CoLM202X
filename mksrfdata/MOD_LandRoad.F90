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
   USE MOD_Vars_Global, only: N_URB, URBAN
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
   type (block_data_int32_2d) :: data_road_class ! urban type index

   ! local vars
   integer, allocatable :: ibuff(:), types(:), order(:)

   ! index
   integer :: ipatch, jpatch, iroad
   integer :: ie, ipxstt, ipxend, npxl, ipxl
   integer :: nr_glb, npatch_glb

   ! local vars for landpath and landroad
   integer :: numpatch_
   integer*8, allocatable :: eindex_(:)
   integer,   allocatable :: ipxstt_(:)
   integer,   allocatable :: ipxend_(:)
   integer,   allocatable :: settyp_(:)
   integer,   allocatable :: ielm_  (:)

   integer  :: numroad_
   integer  :: iurb, ib, imiss
   integer  :: buff_count(N_URB)
   real(r8) :: buff_p(N_URB)

   integer , allocatable :: roadclass (:)
   real(r8), allocatable :: area_one (:)

   character(len=256) :: suffix, cyear

      IF (p_is_master) THEN
         write(*,'(A)') 'Making road type tiles :'
      ENDIF

#ifdef USEMPI
      CALL mpi_barrier (p_comm_glb, p_err)
#endif

      ! allocate and read the grided LCZ/NCAR urban type
      IF (p_is_io) THEN

         dir_road = trim(DEF_dir_rawdata) // '/urban_type'

         CALL allocate_block_data (groad, data_road_class)
         CALL flush_block_data (data_road_class, 0)

         !read LCZ data
         suffix = 'URBTYP'
         CALL read_5x5_data (dir_road, suffix, groad, 'LCZ_DOM', data_road_class)

#ifdef USEMPI
         CALL aggregation_data_daemon (groad, data_i4_2d_in1 = data_road_class)
#endif
      ENDIF

      IF (p_is_worker) THEN

         IF (numpatch > 0) THEN
            ! a temporary numpatch with max urban patch number
            numpatch_ = numpatch + count(landpatch%settyp == URBAN) * (N_URB-1)

            allocate (eindex_ (numpatch_ ))
            allocate (ipxstt_ (numpatch_ ))
            allocate (ipxend_ (numpatch_ ))
            allocate (settyp_ (numpatch_ ))
            allocate (ielm_   (numpatch_ ))

            ! max urban patch number (temporary)
            numroad_ = count(landpatch%settyp == URBAN) * N_URB
            IF (numroad_ > 0) THEN
               allocate (roadclass(numroad_))
            ENDIF
         ENDIF

         jpatch = 0
         iroad = 0

         ! loop for temporary numpatch to filter duplicate urban patch
         DO ipatch = 1, numpatch
            IF (landpatch%settyp(ipatch) == URBAN) THEN

               ie     = landpatch%ielm  (ipatch)
               ipxstt = landpatch%ipxstt(ipatch)
               ipxend = landpatch%ipxend(ipatch)

               CALL aggregation_request_data (landpatch, ipatch, groad, zip = .false., area = area_one, &
                  data_i4_2d_in1 = data_road_class, data_i4_2d_out1 = ibuff)

               ! when there is missing urban types
               !NOTE@tungwz: need duoble check below and add appropriate annotations
               ! check if there is urban pixel without URBAN ID
               imiss = count(ibuff<1 .or. ibuff>N_URB)
               IF (imiss > 0) THEN
                  ! Calculate the relative ratio of each urban types by excluding urban pixels withoht URBAN ID
                  WHERE (ibuff<1 .or. ibuff>N_URB)
                     area_one = 0
                  END WHERE

                  buff_p = 0
                  IF (sum(area_one) > 0) THEN
                     DO ib = 1, size(area_one)
                        IF (ibuff(ib)>1 .and. ibuff(ib)<N_URB) THEN
                           iurb         = ibuff(ib)
                           buff_p(iurb) = buff_p(iurb) + area_one(ib)
                        ENDIF
                     ENDDO
                     buff_p(:) = buff_p(:)/sum(area_one)
                  ENDIF

                  ! The number of URBAN ID of each type is assigned to urban pixels without URBAN ID in relative proportion
                  DO iurb = 1, N_URB-1
                     buff_count(iurb) = int(buff_p(iurb)*imiss)
                  ENDDO
                  buff_count(N_URB) = imiss - sum(buff_count(1:N_URB-1))

                  ! Some urban patches and NCAR/LCZ data are inconsistent (NCAR/LCZ has no urban ID),
                  ! so the these points are assigned
                  IF (all(buff_count==0)) THEN
                     ! If none of the urban pixels have an URBAN ID, they are assigned directly
                     IF (DEF_URBAN_type_scheme == 1) THEN
                        ibuff = 3
                     ELSEIF (DEF_URBAN_type_scheme == 2) THEN
                        ibuff = 9
                     ENDIF
                  ELSE
                     ! Otherwise, URBAN ID are assigned based on the previously calculated number
                     DO ib = 1, size(ibuff)
                        IF (ibuff(ib)<1 .or. ibuff(ib)>N_URB) THEN
                           type_loop: DO iurb = 1, N_URB
                              IF (buff_count(iurb) > 0) THEN
                                 ibuff(ib)        = iurb
                                 buff_count(iurb) = buff_count(iurb) - 1
                                 EXIT type_loop
                              ENDIF
                           ENDDO type_loop
                        ENDIF
                     ENDDO
                  ENDIF
               ENDIF

               npxl = ipxend - ipxstt + 1

               allocate (types (ipxstt:ipxend))

               types(:) = ibuff

               deallocate (ibuff)

               allocate (order (ipxstt:ipxend))
               order = (/ (ipxl, ipxl = ipxstt, ipxend) /)

               ! change order vars, types->regid ? still types below
               ! add region information, because urban type may be same,
               ! but from different region in this urban patch
               ! relative code is changed
               CALL quicksort (npxl, types, order)

               mesh(ie)%ilon(ipxstt:ipxend) = mesh(ie)%ilon(order)
               mesh(ie)%ilat(ipxstt:ipxend) = mesh(ie)%ilat(order)

               DO ipxl = ipxstt, ipxend
                  IF (ipxl /= ipxstt) THEN
                     IF (types(ipxl) /= types(ipxl-1)) THEN
                        ipxend_(jpatch) = ipxl - 1
                     ELSE
                        CYCLE
                     ENDIF
                  ENDIF

                  jpatch = jpatch + 1
                  eindex_(jpatch) = mesh(ie)%indx
                  settyp_(jpatch) = URBAN
                  ipxstt_(jpatch) = ipxl
                  ielm_  (jpatch) = ie

                  iroad = iroad + 1
                  roadclass(iroad) = types(ipxl)
               ENDDO

               ipxend_(jpatch) = ipxend

               deallocate (types)
               deallocate (order)

            ELSE
               jpatch = jpatch + 1
               eindex_(jpatch) = landpatch%eindex(ipatch)
               ipxstt_(jpatch) = landpatch%ipxstt(ipatch)
               ipxend_(jpatch) = landpatch%ipxend(ipatch)
               settyp_(jpatch) = landpatch%settyp(ipatch)
               ielm_  (jpatch) = landpatch%ielm  (ipatch)
            ENDIF
         ENDDO

#ifdef USEMPI
         CALL aggregation_worker_done ()
#endif

         numpatch = jpatch

         IF (numpatch > 0) THEN
            ! update landpath with new patch number
            ! all urban type patch are included
            IF (allocated (landpatch%eindex)) deallocate (landpatch%eindex)
            IF (allocated (landpatch%ipxstt)) deallocate (landpatch%ipxstt)
            IF (allocated (landpatch%ipxend)) deallocate (landpatch%ipxend)
            IF (allocated (landpatch%settyp)) deallocate (landpatch%settyp)
            IF (allocated (landpatch%ielm  )) deallocate (landpatch%ielm  )

            allocate (landpatch%eindex (numpatch))
            allocate (landpatch%ipxstt (numpatch))
            allocate (landpatch%ipxend (numpatch))
            allocate (landpatch%settyp (numpatch))
            allocate (landpatch%ielm   (numpatch))

            ! update all information of landpatch
            landpatch%eindex = eindex_(1:jpatch)
            landpatch%ipxstt = ipxstt_(1:jpatch)
            landpatch%ipxend = ipxend_(1:jpatch)
            landpatch%settyp = settyp_(1:jpatch)
            landpatch%ielm   = ielm_  (1:jpatch)
         ENDIF

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
            landroad%settyp = roadclass(1:numroad)
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

      IF (allocated (ibuff   )) deallocate (ibuff    )
      IF (allocated (types   )) deallocate (types    )
      IF (allocated (order   )) deallocate (order    )

      IF (allocated (eindex_ )) deallocate (eindex_  )
      IF (allocated (ipxstt_ )) deallocate (ipxstt_  )
      IF (allocated (ipxend_ )) deallocate (ipxend_  )
      IF (allocated (settyp_ )) deallocate (settyp_  )
      IF (allocated (ielm_   )) deallocate (ielm_    )

      IF (allocated (roadclass)) deallocate (roadclass )
      IF (allocated (area_one)) deallocate (area_one )

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
