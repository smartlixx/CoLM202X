#include <define.h>

MODULE MOD_Road_Const_ThermalParameters

   USE MOD_Precision

   IMPLICIT NONE
   SAVE


   ! albeodo of asphalt pavement [-]
   real(r8), parameter, dimension(2,2)  :: albroad_apt &
      = reshape([0.03, 0.03, 0.04, 0.04], shape(albroad_apt))

   ! albeodo of concrete pavement [-]
   real(r8), parameter, dimension(2,2)  :: albroad_cct &
      = reshape([0.3, 0.3, 0.5, 0.5], shape(albroad_cct))

   ! emissivity of asphalt pavement [-]
   real(r8), parameter, dimension(2,2)  :: emroad_apt &
      = reshape([0.90, 0.90, 0.92, 0.92], shape(emroad_apt))

   ! emissivity of concrete pavement [-]
   real(r8), parameter, dimension(2,2)  :: emroad_cct &
      = reshape([0.87, 0.87, 0.92, 0.92], shape(emroad_cct))

   

   ! volumetric heat capacity of asphalt pavement [J/kg*K]
   real(r8), parameter, dimension(10,5)  :: cvroad_apt &
      = reshape([2.29E6 , 2.29E6 ,  1.96E6 , 1.45E6 , 2.06E6 ,   2.06E6 , 2.05E6 , 0. , 0. , 0. , &
                 2.29E6 , 2.29E6 ,  1.96E6 , 1.45E6 , 2.06E6 ,   2.06E6 , 2.05E6 , 0. , 0. , 0. , &
                 2.29E6 , 2.29E6 ,  1.96E6 , 1.96E6 , 2.06E6 ,   2.06E6 , 2.05E6 , 0. , 0. , 0. , &
                 2.29E6 , 2.29E6 ,  1.96E6 , 2.06E6 , 2.06E6 ,   2.05E6 , 2.05E6 , 0. , 0. , 0. , &
                 2.29E6 , 2.29E6 ,  2.29E6 , 2.06E6 , 2.06E6 ,   2.05E6 ,     0. , 0. , 0. , 0.], &
                shape(cvroad_apt))

   ! volumetric heat capacity of concrete pavement [J/kg*K]
   real(r8), parameter, dimension(10,5)  :: cvroad_cct &
      = reshape([ 2.15E6 ,  2.15E6 ,  2.15E6 ,  2.15E6 ,  2.15E6 ,    2.11E6 ,  2.05E6 ,  0. ,  0. ,  0. , &
                  2.15E6 ,  2.15E6 ,  2.15E6 ,  2.15E6 ,  2.11E6 ,    2.05E6 ,      0. ,  0. ,  0. ,  0. , &
                  2.15E6 ,  2.15E6 ,  2.15E6 ,  2.15E6 ,  2.06E6 ,    2.05E6 ,      0. ,  0. ,  0. ,  0. , &
                  2.15E6 ,  2.15E6 ,  2.15E6 ,  2.15E6 ,  2.05E6 ,    2.19E6 ,      0. ,  0. ,  0. ,  0. , &
                  2.15E6 ,  2.15E6 ,  2.15E6 ,  2.05E6 ,  2.19E6 ,        0. ,      0. ,  0. ,  0. ,  0.], &
                shape(cvroad_cct))


   ! thermal conductivity of asphalt pavement [W/m*K]
   real(r8), parameter, dimension(10,5)  :: tkroad_apt &
      = reshape([0.90 , 0.90 , 1.05 , 1.30 , 1.20 ,   1.20 , 1.30 ,   0. ,   0. ,   0. , &
                 0.90 , 0.90 , 1.05 , 1.05 , 1.20 ,   1.20 , 1.30 ,   0. ,   0. ,   0. , &
                 0.90 , 0.90 , 1.05 , 1.05 , 1.20 ,   1.20 , 1.30 ,   0. ,   0. ,   0. , &
                 0.90 , 0.90 , 1.05 , 1.20 , 1.20 ,   1.30 ,   0. ,   0. ,   0. ,   0. , &
                 0.90 , 0.90 , 0.90 , 1.20 , 1.20 ,   1.30 ,   0. ,   0. ,   0. ,   0.], &
                shape(tkroad_apt))

   ! thermal conductivity of concrete pavement [W/m*K]
   real(r8), parameter, dimension(10,5)  :: tkroad_cct &
      = reshape([1.11 , 1.11 , 1.11 , 1.11 , 1.11 ,   0.92 , 1.30 ,   0. ,   0. ,   0. , &
                 1.11 , 1.11 , 1.11 , 1.11 , 0.92 ,   1.30 ,   0. ,   0. ,   0. ,   0. , &
                 1.11 , 1.11 , 1.11 , 1.11 , 1.20 ,   1.30 ,   0. ,   0. ,   0. ,   0. , &
                 1.11 , 1.11 , 1.11 , 1.11 , 1.30 ,   1.45 ,   0. ,   0. ,   0. ,   0. , &
                 1.11 , 1.11 , 1.11 , 1.30 , 1.45 ,     0. ,   0. ,   0. ,   0. ,   0.], &
                shape(tkroad_cct))

END MODULE MOD_Road_Const_ThermalParameters
