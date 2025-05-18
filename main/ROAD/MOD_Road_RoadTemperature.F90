#include <define.h>

MODULE MOD_Road_RoadTemperature

   USE MOD_Precision
   IMPLICIT NONE
   SAVE

   PUBLIC :: RoadTemperature

CONTAINS

   SUBROUTINE RoadTemperature (patchtype,lbroad,deltim, &
                              capr,cnfac,vf_quartz,vf_gravels,vf_om,vf_sand,&
                              wf_gravels,wf_sand,porsl,psi0,&
                              csol,k_solids,dksatu,dksatf,dkdry,&
                              BA_alpha, BA_beta,&
                              cv_road,tk_road,&
                              dz_roadsno,z_roadsno,zi_roadsno,&
                              t_roadsno,t_road,wice_roadsno,wliq_roadsno,&
                              scv_road,snowdp_road,&
                           !   frl,
                              dlrad,&!clroad,
                              sab_road,&
                              fsen_road,fsen_snow,fevp_road,fevp_snow,&
                              croad,htvp,em_road,&
                              imelt,sm,xmf,fact,&
                              pgroad_rain,pgroad_snow,t_precip)

!=======================================================================
! Snow and road temperatures
! o The volumetric heat capacity is calculated as a linear combination
!   in terms of the volumetric fraction of the constituent phases.
! o The thermal conductivity of road soil is computed from
!   the algorithm of Johansen (as reported by Farouki 1981), impervious and perivious from
!   LOOK-UP table and of snow is from the formulation used in SNTHERM (Jordan 1991).
! o Boundary conditions:
!   F = Rnet - Hg - LEg (top),  F = 0 (base of the soil column).
! o Soil / snow temperature is predicted from heat conduction
!   in 10 soil layers and up to 5 snow layers.
!   The thermal conductivities at the interfaces between two neighbor layers
!   (j, j+1) are derived from an assumption that the flux across the interface
!   is equal to that from the node j to the interface and the flux from the
!   interface to the node j+1. The equation is solved using the Crank-Nicholson
!   method and resulted in a tridiagonal system equation.
!
! Phase change (see MOD_PhaseChange.F90)
!
! Original author : Yongjiu Dai, 09/15/1999; 08/30/2002; 05/2020
!=======================================================================

   USE MOD_Precision
   USE MOD_Vars_Global
   USE MOD_Const_Physical
  !USE MOD_Road_Const_ThermalParameters
   USE MOD_SoilThermalParameters
   USE MOD_PhaseChange, only: meltf_road
   USE MOD_Utils, only: tridia

   IMPLICIT NONE

   integer, intent(in)  :: lbroad                      !lower bound of array
   integer, intent(in)  :: patchtype                   !land patch type (0=soil,1=urban or built-up,2=wetland,
                                                       !3=land ice, 4=deep lake, 5=shallow lake)
   real(r8), intent(in) :: deltim                      !seconds in a time step [second]
   real(r8), intent(in) :: capr                        !tuning factor to turn first layer T into surface T
   real(r8), intent(in) :: cnfac                       !Crank Nicholson factor between 0 and 1

   real(r8), intent(in) :: csol      (1:nl_soil)       !heat capacity of soil solids [J/(m3 K)]
   real(r8), intent(in) :: k_solids  (1:nl_soil)       !thermal conductivity of minerals soil [W/m-K]
   real(r8), intent(in) :: porsl     (1:nl_soil)       !soil porosity [-]
   real(r8), intent(in) :: psi0      (1:nl_soil)       !soil water suction, negative potential [mm]

   real(r8), intent(in) :: dkdry     (1:nl_soil)       !thermal conductivity of dry soil [W/m-K]
   real(r8), intent(in) :: dksatu    (1:nl_soil)       !thermal conductivity of saturated soil [W/m-K]
   real(r8), intent(in) :: dksatf    (1:nl_soil)       !thermal conductivity of saturated frozen soil [W/m-K]

   real(r8), intent(in) :: vf_quartz (1:nl_soil)       !volumetric fraction of quartz within mineral soil
   real(r8), intent(in) :: vf_gravels(1:nl_soil)       !volumetric fraction of gravels
   real(r8), intent(in) :: vf_om     (1:nl_soil)       !volumetric fraction of organic matter
   real(r8), intent(in) :: vf_sand   (1:nl_soil)       !volumetric fraction of sand
   real(r8), intent(in) :: wf_gravels(1:nl_soil)       !gravimetric fraction of gravels
   real(r8), intent(in) :: wf_sand   (1:nl_soil)       !gravimetric fraction of sand

   real(r8), intent(in) :: BA_alpha  (1:nl_soil)       !alpha in Balland and Arp(2005) thermal conductivity scheme
   real(r8), intent(in) :: BA_beta   (1:nl_soil)       !beta in Balland and Arp(2005) thermal conductivity scheme

   real(r8), intent(in) :: cv_road   (1:nl_soil)       !heat capacity of road [J/m3/K]
   real(r8), intent(in) :: tk_road   (1:nl_soil)       !thermal conductivity of road [W/m/K]

   real(r8), intent(in) :: dz_roadsno(lbroad  :nl_soil)    !layer thickiness [m]
   real(r8), intent(in) :: z_roadsno (lbroad  :nl_soil)    !node depth [m]
   real(r8), intent(in) :: zi_roadsno(lbroad-1:nl_soil)    !interface depth [m]

   real(r8), intent(in) :: t_road                      !road surface temperature [K]
   real(r8), intent(in) :: sab_road                    !solar radiation absorbed by ground [W/m2]
!   real(r8), intent(in) :: frl                         !atmospheric infrared (longwave) radiation [W/m2]
!   real(r8), intent(in) :: clgimp                      !deriv. of longwave wrt to soil temp [w/m2/k]
   real(r8), intent(in) :: dlrad                       !downward longwave radiation blow the canopy [W/m2]
!   real(r8), intent(in) :: lroad                       !atmospheric infrared (longwave) radiation [W/m2]
!   real(r8), intent(in) :: clroad                      !deriv. of longwave wrt to soil temp [w/m2/k]
   real(r8), intent(in) :: croad                       !deriv. of soil energy flux wrt to soil temp [w/m2/k]
     
   real(r8), intent(in) :: fsen_road                   !sensible heat flux from ground [W/m2]
   real(r8), intent(in) :: fsen_snow                   !sensible heat flux from ground snow [W/m2]
   real(r8), intent(in) :: fevp_road                   !evaporation heat flux from ground [mm/s]
   real(r8), intent(in) :: fevp_snow                   !evaporation heat flux from ground snow [mm/s]
   !   real(r8), intent(in) :: cgimp                   !deriv. of soil energy flux wrt to soil temp [w/m2/k]
   real(r8), intent(in) :: htvp                        !latent heat of vapor of water (or sublimation) [j/kg]
   real(r8), intent(in) :: em_road                     !road emissivity
  
   real(r8), intent(inout) :: t_roadsno   (lbroad:nl_soil) !soil temperature [K]
   real(r8), intent(inout) :: wice_roadsno(lbroad:nl_soil) !ice lens [kg/m2]
   real(r8), intent(inout) :: wliq_roadsno(lbroad:nl_soil) !liqui water [kg/m2]
   real(r8), intent(inout) :: scv_road                 !snow cover, water equivalent [mm, kg/m2]
   real(r8), intent(inout) :: snowdp_road              !snow depth [m]

   real(r8), intent(out) :: sm                         !rate of snowmelt [kg/(m2 s)]
   real(r8), intent(out) :: xmf                        !total latent heat of phase change of ground water
   real(r8), intent(out) :: fact (lbroad:nl_soil)      !used in computing tridiagonal matrix
   integer,  intent(out) :: imelt(lbroad:nl_soil)      !flag for melting or freezing [-]
   real(r8), intent(in)  :: pgroad_rain                ! rainfall onto road including canopy runoff [kg/(m2 s)]
   real(r8), intent(in)  :: pgroad_snow                ! snowfall onto road including canopy runoff [kg/(m2 s)]
   real(r8), intent(in)  :: t_precip                   ! snowfall/rainfall temperature [kelvin]

!------------------------ local variables ------------------------------
   real(r8) cv (lbroad:nl_soil)           !heat capacity [J/(m2 K)]
   real(r8) tk (lbroad:nl_soil)           !thermal conductivity [W/(m K)]

   real(r8) hcap(1:nl_soil)               !J/(m3 K)
   real(r8) thk(lbroad:nl_soil)           !W/(m K)
   real(r8) rhosnow                       !partitial density of water (ice + liquid)

   real(r8) at (lbroad:nl_soil)           !"a" vector for tridiagonal matrix
   real(r8) bt (lbroad:nl_soil)           !"b" vector for tridiagonal matrix
   real(r8) ct (lbroad:nl_soil)           !"c" vector for tridiagonal matrix
   real(r8) rt (lbroad:nl_soil)           !"r" vector for tridiagonal solution

   real(r8) fn (lbroad:nl_soil)           !heat diffusion through the layer interface [W/m2]
   real(r8) fn1(lbroad:nl_soil)           !heat diffusion through the layer interface [W/m2]
   real(r8) dzm                       !used in computing tridiagonal matrix
   real(r8) dzp                       !used in computing tridiagonal matrix

   real(r8) t_roadsno_bef(lbroad:nl_soil) !soil/snow temperature before update
   real(r8) hs                        !net energy flux into the surface (w/m2)
   real(r8) dhsdt                     !d(hs)/dT
   real(r8) brr(lbroad:nl_soil)           !temporay set

   real(r8) vf_water(1:nl_soil)       !volumetric fraction liquid water within soil
   real(r8) vf_ice  (1:nl_soil)       !volumetric fraction ice len within soil

   integer i,j

   wice_roadsno(2:) = 0.0             !ice lens [kg/m2]
   wliq_roadsno(2:) = 0.0             !liquid water [kg/m2]

!=======================================================================
! soil ground and wetland heat capacity
   DO i = 1, nl_soil
      vf_water(i) = wliq_roadsno(i)/(dz_roadsno(i)*denh2o)
      vf_ice(i)   = wice_roadsno(i)/(dz_roadsno(i)*denice)
      CALL soil_hcap_cond(vf_gravels(i),vf_om(i),vf_sand(i),porsl(i),&
                        wf_gravels(i),wf_sand(i),k_solids(i),&
                        csol(i),dkdry(i),dksatu(i),dksatf(i),&
                        BA_alpha(i),BA_beta(i),&
                        t_roadsno(i),vf_water(i),vf_ice(i),hcap(i),thk(i))
      cv(i) = hcap(i)*dz_roadsno(i)
   ENDDO
   IF(lbroad==1 .and. scv_road>0.) cv(1) = cv(1) + cpice*scv_road

! Snow heat capacity
   IF(lbroad <= 0) THEN
      cv(:0) = cpliq*wliq_roadsno(:0) + cpice*wice_roadsno(:0)
   ENDIF

! Snow thermal conductivity
   IF(lbroad <= 0) THEN
      DO i = lbroad, 0
         rhosnow = (wice_roadsno(i)+wliq_roadsno(i))/dz_roadsno(i)

         ! presently option [1] is the default option
         ! [1] Jordan (1991) pp. 18
         thk(i) = tkair+(7.75e-5*rhosnow+1.105e-6*rhosnow*rhosnow)*(tkice-tkair)

            ! [2] Sturm et al (1997)
            ! thk(i) = 0.0138 + 1.01e-3*rhosnow + 3.233e-6*rhosnow**2
            ! [3] Ostin and Andersson presented in Sturm et al., (1997)
            ! thk(i) = -0.871e-2 + 0.439e-3*rhosnow + 1.05e-6*rhosnow**2
            ! [4] Jansson(1901) presented in Sturm et al. (1997)
            ! thk(i) = 0.0293 + 0.7953e-3*rhosnow + 1.512e-12*rhosnow**2
            ! [5] Douville et al., (1995)
            ! thk(i) = 2.2*(rhosnow/denice)**1.88
            ! [6] van Dusen (1992) presented in Sturm et al. (1997)
            ! thk(i) = 0.021 + 0.42e-3*rhosnow + 0.22e-6*rhosnow**2

      ENDDO
   ENDIF

! Thermal conductivity at the layer interface
   DO i = lbroad, nl_soil-1

! the following consideration is try to avoid the snow conductivity
! to be dominant in the thermal conductivity of the interface.
! Because when the distance of bottom snow node to the interfacee
! is larger than that of interface to top soil node,
! the snow thermal conductivity will be dominant, and the result is that
! lees heat tranfer between snow and soil
      IF((i==0) .and. (z_roadsno(i+1)-zi_roadsno(i)<zi_roadsno(i)-z_roadsno(i)))THEN
         tk(i) = 2.*thk(i)*thk(i+1)/(thk(i)+thk(i+1))
         tk(i) = max(0.5*thk(i+1),tk(i))
      ELSE
         tk(i) = thk(i)*thk(i+1)*(z_roadsno(i+1)-z_roadsno(i)) &
               /(thk(i)*(z_roadsno(i+1)-zi_roadsno(i))+thk(i+1)*(zi_roadsno(i)-z_roadsno(i)))
      ENDIF
   ENDDO
   tk(nl_soil) = 0.

   WHERE (tk_road > 0.) tk(1:) = tk_road(1:)
   WHERE (cv_road > 0.) cv(1:) = cv_road(1:)*dz_roadsno(1:)

   ! snow exist when there is no snow layer
   IF (lbroad == 1 .and. scv_road > 0.0) THEN
      cv(1) = cv(1) + cpice*scv_road
   ENDIF

   ! ponding water or ice exist
   cv(1) = cv(1) + cpliq*wliq_roadsno(1) + cpice*wice_roadsno(1)

   ! net ground heat flux into the surface and its temperature derivative
   hs = sab_road + dlrad*em_road &
      - (fsen_road+fevp_road*htvp) &
      + cpliq*pgroad_rain*(t_precip-t_road) &
      + cpice*pgroad_snow*(t_precip-t_road) &
      - em_road*stefnc*t_road**4.
   dhsdT = - croad -4.*em_road*stefnc*t_road**3. - cpliq*pgroad_rain - cpice*pgroad_snow

   t_roadsno_bef(lbroad:) = t_roadsno(lbroad:)

   j       = lbroad
   fact(j) = deltim / cv(j) &
         * dz_roadsno(j) / (0.5*(z_roadsno(j)-zi_roadsno(j-1)+capr*(z_roadsno(j+1)-zi_roadsno(j-1))))

   DO j = lbroad + 1, nl_soil
      fact(j) = deltim/cv(j)
   ENDDO

   DO j = lbroad, nl_soil - 1
      fn(j) = tk(j)*(t_roadsno(j+1)-t_roadsno(j))/(z_roadsno(j+1)-z_roadsno(j))
   ENDDO
   fn(nl_soil) = 0.

! set up vector r and vectors a, b, c that define tridiagonal matrix
   j     = lbroad
   dzp   = z_roadsno(j+1)-z_roadsno(j)
   at(j) = 0.
   bt(j) = 1+(1.-cnfac)*fact(j)*tk(j)/dzp-fact(j)*dhsdT
   ct(j) =  -(1.-cnfac)*fact(j)*tk(j)/dzp
   rt(j) = t_roadsno(j) + fact(j)*( hs - dhsdT*t_roadsno(j) + cnfac*fn(j) )


   DO j = lbroad + 1, nl_soil - 1
      dzm   = (z_roadsno(j)-z_roadsno(j-1))
      dzp   = (z_roadsno(j+1)-z_roadsno(j))
      at(j) =   - (1.-cnfac)*fact(j)* tk(j-1)/dzm
      bt(j) = 1.+ (1.-cnfac)*fact(j)*(tk(j)/dzp + tk(j-1)/dzm)
      ct(j) =   - (1.-cnfac)*fact(j)* tk(j)/dzp
      rt(j) = t_roadsno(j) + cnfac*fact(j)*( fn(j) - fn(j-1) )
   ENDDO

   j     =  nl_soil
   dzm   = (z_roadsno(j)-z_roadsno(j-1))
   at(j) =   - (1.-cnfac)*fact(j)*tk(j-1)/dzm
   bt(j) = 1.+ (1.-cnfac)*fact(j)*tk(j-1)/dzm
   ct(j) = 0.
   rt(j) = t_roadsno(j) - cnfac*fact(j)*fn(j-1)

!   print *, 'Snow thermal properties, tk = ', tk(lbroad:1), 'cv =', cv(lbroad:1)
!   print *, 'sabroad = ', sab_road, 'em_road = ', em_road, 'lnet = ', em_road*(dlrad - stefnc*t_road**4.)
!   print *, 'fsenroad = ', fsen_road, 'fevp_road*htvp = ', fevp_road*htvp
!   print *, 'snow and rain energy exchange = ', cpliq*pgroad_rain*(t_precip-t_road), cpice*pgroad_snow*(t_precip-t_road)
!   print *, 'hs = ', hs, 'dhsdT = ', dhsdT

! solve for t_roadsno
   i = size(at)
   CALL tridia (i ,at ,bt ,ct ,rt ,t_roadsno)

!=======================================================================
! melting or freezing
!=======================================================================

   DO j = lbroad, nl_soil - 1
      fn1(j) = tk(j)*(t_roadsno(j+1)-t_roadsno(j))/(z_roadsno(j+1)-z_roadsno(j))
   ENDDO
   fn1(nl_soil) = 0.

   j = lbroad
   brr(j) = cnfac*fn(j) + (1.-cnfac)*fn1(j)

   DO j = lbroad + 1, nl_soil
      brr(j) = cnfac*(fn(j)-fn(j-1)) + (1.-cnfac)*(fn1(j)-fn1(j-1))
   ENDDO

!  can also call meltf_urban directly
   CALL meltf_road (lbroad,1,deltim, &
            fact(lbroad:1),brr(lbroad:1),hs,dhsdT, &
            t_roadsno_bef(lbroad:1),t_roadsno(lbroad:1), &
            wliq_roadsno(lbroad:1),wice_roadsno(lbroad:1),imelt(lbroad:1), &
            scv_road,snowdp_road,sm,xmf)

   END SUBROUTINE RoadTemperature

END MODULE MOD_Road_RoadTemperature
! ---------- EOP ------------
