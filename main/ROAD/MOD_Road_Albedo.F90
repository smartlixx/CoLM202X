#include <define.h>

MODULE MOD_Road_Albedo
!-----------------------------------------------------------------------
! !DESCRIPTION:
! Calculate road albedo,
!
! Created by Hua Yuan, 09/2021
! Adapted for road by Xianxiang Li, 09/2024
!
! REVISIONS:
!
!
!-----------------------------------------------------------------------
   USE MOD_Precision
   IMPLICIT NONE
   SAVE

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: albroad

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------

   SUBROUTINE albroad (ipatch,alb_road,coszen,fsno_road,&
                       scv_road,sag_road,alb,sroad)

!=======================================================================
! Calculates fragmented albedos (direct and diffuse) in
! wavelength regions split at 0.7um.
!
! (1) snow albedos: as in BATS formulations, which are inferred from
!     the calculations of Wiscombe and Warren (1980) and the snow model
!     and data of Anderson(1976), and the function of snow age, grain size,
!     solar zenith angle, pollution, the amount of the fresh snow
! (2) over the snow covered surface, the surface albedo is estimated by a linear
!     combination of albedos for snow, roof, impervious and pervious ground
!
!=======================================================================

   USE MOD_Precision
   USE MOD_Const_Physical, only: tfrz

   IMPLICIT NONE

!------------------------- Dummy Arguments -----------------------------
! ground cover index
   integer, intent(in) :: &
      ipatch          ! patch index

   real(r8), intent(in) :: &
      alb_road(2,2), &! impervious albedo (iband,direct/diffuse)

      coszen,    &! cosine of solar zenith angle [-]
      fsno_road, &! fraction of soil covered by snow [-]
   
      scv_road,  &! snow cover, water equivalent [mm]
      sag_road    ! non dimensional snow age [-]

   real(r8), intent(out) :: &
      alb(2,2),  &! averaged albedo [-]
      sroad(2,2)  ! road absorption for solar radiation,
   

!-------------------------- Local variables ----------------------------
   real(r8) :: &!
      age,       &! factor to reduce visible snow alb due to snow age [-]
      cff,       &! snow alb correction factor for zenith angle > 60 [-]
      conn,      &! constant (=0.5) for visible snow alb calculation [-]
      cons,      &! constant (=0.2) for nir snow albedo calculation [-]
      czen,      &! cosine of solar zenith angle > 0 [-]
      czf,       &! solar zenith correction for new snow albedo [-]
      dfalbl,    &! snow albedo for diffuse nir radiation [-]
      dfalbs,    &! snow albedo for diffuse visible solar radiation [-]
      dralbl,    &! snow albedo for visible radiation [-]
      dralbs,    &! snow albedo for near infrared radiation [-]
      sl,        &! factor that helps control alb zenith dependence [-]
      snal0,     &! alb for visible,incident on new snow (zen ang<60) [-]
      snal1,     &! alb for NIR, incident on new snow (zen angle<60) [-]
      tran(2,3)   ! canopy transmittances for solar radiation

   real(r8) :: &!
      albsno(2,2),   &! snow albedo [-]
      albroad_(2,2)   ! albedo, ground
   
! ----------------------------------------------------------------------
! 1. Initial set
! ----------------------------------------------------------------------

! short and long wave albedo for new snow
      snal0 = 0.85     ! shortwave
      snal1 = 0.65     ! long wave

! ----------------------------------------------------------------------
! set default soil and road albedos and solar absorption
      alb (:,:)  = 0. ! averaged
      sroad(:,:) = 0.
 
      tran(:,1) = 0.       !incident direct  radiation diffuse transmittance
      tran(:,2) = 1.       !incident diffuse radiation diffuse transmittance
      tran(:,3) = 1.       !incident direct  radiation direct  transmittance

      IF(coszen <= 0.) THEN
         !print *, "coszen < 0, ipatch and coszen: ", ipatch, coszen
         RETURN  !only do albedo when coszen > 0
      ENDIF

      czen = max(coszen,0.01)
      albsno(:,:) = 0.    !set initial snow albedo
      cons = 0.2          !parameter for snow albedo
      conn = 0.5          !parameter for snow albedo
      sl  = 2.0           !sl helps control albedo zenith dependence

! ----------------------------------------------------------------------
! 2. get albedo over road
! ----------------------------------------------------------------------

! 2.3 road albedo with snow
      IF (scv_road > 0.) THEN

         ! correction for snow age
         age = 1.-1./(1.+sag_road) !correction for snow age
         dfalbs = snal0*(1.-cons*age)

         ! czf corrects albedo of new snow for solar zenith
         cff    = ((1.+1./sl)/(1.+czen*2.*sl)- 1./sl)
         cff    = max(cff,0.)
         czf    = 0.4*cff*(1.-dfalbs)
         dralbs = dfalbs+czf
         dfalbl = snal1*(1.-conn*age)
         czf    = 0.4*cff*(1.-dfalbl)
         dralbl = dfalbl+czf

         albsno(1,1) = dralbs
         albsno(2,1) = dralbl
         albsno(1,2) = dfalbs
         albsno(2,2) = dfalbl

      ENDIF

      albroad_(:,:) = (1.-fsno_road)*alb_road(:,:) + fsno_road*albsno(:,:)

! ----------------------------------------------------------------------
! 3. Road albedo
! ----------------------------------------------------------------------

      alb(:,:) = albroad_(:,:)

      ! treat road absorption in direct and diffuse respectively
      sroad(1,1) = tran(1,1)*(1.-albroad_(1,2)) + tran(1,3)*(1-albroad_(1,1))
      sroad(2,1) = tran(2,1)*(1.-albroad_(2,2)) + tran(2,3)*(1-albroad_(2,1))
      sroad(1,2) = tran(1,2)*(1.-albroad_(1,2))
      sroad(2,2) = tran(2,2)*(1.-albroad_(2,2))

   END SUBROUTINE albroad

END MODULE MOD_Road_Albedo
! --------- EOP ----------
