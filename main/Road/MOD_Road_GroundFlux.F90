#include <define.h>

MODULE MOD_Road_GroundFlux

   USE MOD_Precision
   IMPLICIT NONE
   SAVE

   PUBLIC :: RoadGroundFlux

CONTAINS

   SUBROUTINE RoadGroundFlux (hu, ht, hq, us, vs, tm, qm, rhoair, psrf, hpbl, &
                               ur, thm, th, thv, zlnd, zsno, fsno_road, lbroad, &
                              ! fcover, &
                               rss, dqroaddT, htvp, cgrnd, cgrndl, cgrnds, &
                               troad, qroad, t_snow, q_snow, &
                               taux,tauy,fseng_road,fseng_snow, &
                               fevpg_road,fevpg_snow,tref, qref, &
                               z0m, z0hg, zol, rib, ustar, qstar, tstar, &
                               fm, fh, fq)

!=======================================================================
! this is the main subroutine to execute the calculation
! of road ground fluxes
!
!=======================================================================

   USE MOD_Precision
   USE MOD_Const_Physical, only: cpair,vonkar,grav
   USE MOD_FrictionVelocity
   USE mod_namelist, only: DEF_USE_CBL_HEIGHT,DEF_RSS_SCHEME
   USE MOD_TurbulenceLEddy
   IMPLICIT NONE

!----------------------- Dummy argument --------------------------------
   integer , intent(in) :: &
        lbroad 
   real(r8), intent(in) :: &
        ! atmospherical variables and observational height
        hu,       &! observational height of wind [m]
        ht,       &! observational height of temperature [m]
        hq,       &! observational height of humidity [m]
        hpbl,     &! atmospheric boundary layer height [m]
        us,       &! wind component in eastward direction [m/s]
        vs,       &! wind component in northward direction [m/s]
        tm,       &! temperature at agcm reference height [kelvin] [not used]
        qm,       &! specific humidity at agcm reference height [kg/kg]
        rhoair,   &! density air [kg/m3]
        psrf,     &! atmosphere pressure at the surface [pa] [not used]

        ur,       &! wind speed at reference height [m/s]
        thm,      &! intermediate variable (tm+0.0098*ht)
        th,       &! potential temperature (kelvin)
        thv,      &! virtual potential temperature (kelvin)

        zlnd,     &! roughness length for soil [m]
        zsno,     &! roughness length for snow [m]
        fsno_road,&! fraction of road covered by snow
!        fcover(0:5),&! coverage of aboveground urban components [-]

!        wliq_roadsno,&! liqui water [kg/m2]
!        wice_roadsno,&! ice lens [kg/m2]

        troad,     &! road temperature [K]
!        tgper,    &! ground pervious temperature [K]
        t_snow,    &! ground snow temperature [K]
        qroad,     &! road specific humidity [kg/kg]
!        qgper     &! ground pervious specific humidity [kg/kg]
        q_snow,    &! ground snow specific humidity [kg/kg]
        dqroaddT,  &! d(qg)/dT
        rss,       &! soil surface resistance for evaporation [s/m]
        htvp        ! latent heat of vapor of water (or sublimation) [j/kg]
        
   real(r8), intent(out) :: &
        taux,        &! wind stress: E-W [kg/m/s**2]
        tauy,        &! wind stress: N-S [kg/m/s**2]
        fseng_road,  &! sensible heat flux from road [W/m2]
        fseng_snow,  &! sensible heat flux from ground snow [W/m2]
        fevpg_road,  &! evaporation heat flux from road [mm/s]
        fevpg_snow,  &! evaporation heat flux from ground snow [mm/s]
        cgrnd,       &! deriv. of soil energy flux wrt to soil temp [w/m2/k]
        cgrndl,      &! deriv, of soil sensible heat flux wrt soil temp [w/m2/k]
        cgrnds,      &! deriv of soil latent heat flux wrt soil temp [w/m**2/k]
        tref,        &! 2 m height air temperature [kelvin]
        qref          ! 2 m height air humidity

   real(r8), intent(out) :: &
        z0m,      &! effective roughness [m]
        z0hg,     &! roughness length over ground, sensible heat [m]
        zol,      &! dimensionless height (z/L) used in Monin-Obukhov theory
        rib,      &! bulk Richardson number in surface layer
        ustar,    &! friction velocity [m/s]
        tstar,    &! temperature scaling parameter
        qstar,    &! moisture scaling parameter
        fm,       &! integral of profile function for momentum
        fh,       &! integral of profile function for heat
        fq         ! integral of profile function for moisture

!------------------------ LOCAL VARIABLES ------------------------------
   integer niters, &! maximum number of iterations for surface temperature
        iter,     &! iteration index
        nmozsgn    ! number of times moz changes sign

   real(r8) ::     &
        beta,     &! coefficient of conective velocity [-]
        displax,  &! zero-displacement height [m]
      !  tg,       &! ground surface temperature [K]
      !  qg,       &! weighted ground specific humidity [kg/kg]
      !  fg,       &! ground fractional cover [-]
      !  froad,    &! weight of impervious ground
      !  fgper,    &! weight of pervious ground
        dth,      &! diff of virtual temp. between ref. height and surface
        dqh,      &! diff of humidity between ref. height and surface
        dthv,     &! diff of vir. poten. temp. between ref. height and surface
        obu,      &! monin-obukhov length (m)
        obuold,   &! monin-obukhov length from previous iteration
        ram,      &! aerodynamical resistance [s/m]
        rah,      &! thermal resistance [s/m]
        raw,      &! moisture resistance [s/m]
        raih,     &! temporary variable [kg/m2/s]
        raiw,     &! temporary variable [kg/m2/s]
        fh2m,     &! relation for temperature at 2m
        fq2m,     &! relation for specific humidity at 2m
        fm10m,    &! integral of profile function for momentum at 10m
        thvstar,  &! virtual potential temperature scaling parameter
        um,       &! wind speed including the stablity effect [m/s]
        wc,       &! convective velocity [m/s]
        wc2,      &! wc**2
        zeta,     &! dimensionless height used in Monin-Obukhov theory
        zii,      &! convective boundary height [m]
        zldis,    &! reference height "minus" zero displacement heght [m]
        z0mg,     &! roughness length over ground, momentum [m]
        z0qg       ! roughness length over ground, latent heat [m]

!   real(r8) fwet_road, fwetfac

!----------------------- Dummy argument --------------------------------
! initial roughness length
      !NOTE: change to original
      !z0mg = (1.-fsno)*zlnd + fsno*zsno
      IF (fsno_road > 0) THEN
         z0mg = zsno
      ELSE
         z0mg = zlnd
      ENDIF
      z0hg = z0mg
      z0qg = z0mg

! potential temperatur at the reference height
      beta = 1.       !-  (in computing W_*)
      zii  = 1000.    !m  (pbl height)
      z0m  = z0mg

!-----------------------------------------------------------------------
!     Compute sensible and latent fluxes and their derivatives with respect
!     to ground temperature using ground temperatures from previous time step.
!-----------------------------------------------------------------------
! Initialization variables
      nmozsgn = 0
      obuold  = 0.

      dth   = thm-troad
      dqh   = qm-qroad
      dthv  = dth*(1.+0.61*qm)+0.61*th*dqh
      zldis = hu-0.

      CALL moninobukini(ur,th,thm,thv,dth,dqh,dthv,zldis,z0mg,um,obu)

! Evaluated stability-dependent variables using moz from prior iteration
      niters=6

      !----------------------------------------------------------------
      ITERATION : DO iter = 1, niters         !begin stability iteration
      !----------------------------------------------------------------
         displax = 0.
         IF (DEF_USE_CBL_HEIGHT) THEN
            CALL moninobuk_leddy(hu,ht,hq,displax,z0mg,z0hg,z0qg,obu,um, hpbl, &
                                 ustar,fh2m,fq2m,fm10m,fm,fh,fq)
         ELSE
            CALL moninobuk(hu,ht,hq,displax,z0mg,z0hg,z0qg,obu,um,&
                           ustar,fh2m,fq2m,fm10m,fm,fh,fq)
         ENDIF

         tstar = vonkar/fh*dth
         qstar = vonkar/fq*dqh

         z0hg = z0mg/exp(0.13 * (ustar*z0mg/1.5e-5)**0.45)
         z0qg = z0hg

         thvstar=tstar*(1.+0.61*qm)+0.61*th*qstar
         zeta=zldis*vonkar*grav*thvstar/(ustar**2*thv)
         IF (zeta >= 0.) THEN     !stable
           zeta = min(2.,max(zeta,1.e-6))
         ELSE                     !unstable
           zeta = max(-100.,min(zeta,-1.e-6))
         ENDIF
         obu = zldis/zeta

         IF (zeta >= 0.) THEN
           um = max(ur,0.1)
         ELSE
           IF (DEF_USE_CBL_HEIGHT) THEN !//TODO: Shaofeng, 2023.05.18
               zii = max(5.*hu,hpbl)
           ENDIF !//TODO: Shaofeng, 2023.05.18
           wc = (-grav*ustar*thvstar*zii/thv)**(1./3.)
           wc2 = beta*beta*(wc*wc)
           um = sqrt(ur*ur+wc2)
         ENDIF

         IF (obuold*obu < 0.) nmozsgn = nmozsgn+1
         IF (nmozsgn >= 4) EXIT

         obuold = obu

      !----------------------------------------------------------------
      ENDDO ITERATION                         !end stability iteration
      !----------------------------------------------------------------

      zol = zeta
      rib = min(5.,zol*ustar**2/(vonkar**2/fh*um**2))

   ! Get derivative of fluxes with repect to ground temperature
      ram  = 1./(ustar*ustar/um)
      rah  = 1./(vonkar/fh*ustar)
      raw  = 1./(vonkar/fq*ustar)

      raih = rhoair*cpair/rah

   ! 08/23/2019, yuan: add soil surface resistance (rss)
      IF (dqh > 0.) THEN
         raiw = rhoair/raw !dew case. assume no soil resistance
      ELSE
         IF (DEF_RSS_SCHEME .eq. 4) THEN
            raiw = rss*rhoair/raw
         ELSE
            raiw = rhoair/(raw+rss)
         ENDIF
      ENDIF

      cgrnds = raih
      cgrndl = raiw*dqroaddT
      cgrnd  = cgrnds + htvp*cgrndl

   ! surface fluxes of momentum, sensible and latent
   ! using ground temperatures from previous time step
      taux  = -rhoair*us/ram
      tauy  = -rhoair*vs/ram
      fseng_road = -raih*dth
      fevpg_road = -raiw*dqh

      fseng_snow = -raih * (thm - t_snow)
      fevpg_snow = -raiw * ( qm - q_snow)

! 2 m height air temperature
      tref   = thm + vonkar/fh*dth * (fh2m/vonkar - fh/vonkar)
      qref   =  qm + vonkar/fq*dqh * (fq2m/vonkar - fq/vonkar)

   END SUBROUTINE RoadGroundFlux

END MODULE MOD_Road_GroundFlux
