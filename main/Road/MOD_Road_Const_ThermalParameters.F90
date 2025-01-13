#include <define.h>

MODULE MOD_Road_Const_ThermalParameters

   USE MOD_Precision

   IMPLICIT NONE
   SAVE


   ! albeodo of impervious asphalt pavement [-]
   real(r8), parameter, dimension(2,2)  :: albroad_apt &
      = (/0.03, 0.03, 0.04, 0.04/)

   ! albeodo of pervious concrete pavement [-]
   real(r8), parameter, dimension(2,2)  :: albroad_cct &
      = (/0.3, 0.3, 0.5, 0.5/)


   ! volumetric heat capacity of asphalt pavement [J/m3*K]
   real(r8), parameter, dimension(10,5)  :: cvroad_apt &
      = reshape([1168 , 1168 , 1168 , 1168 , 1168 ,   1168 , 1168 , 1168 , 1168 , 1168 , &
                  785 ,  785 ,  785 ,  785 , 1168 ,    815 ,  815 ,  785 ,  825 ,  825 , &
                  825 ,  825 ,  825 ,  825 ,  825 ,    825 ,  825 ,  825 ,  820 ,  820 , &
                  820 ,  820 ,    0 ,    0 ,    0 ,      0 ,    0 ,    0 ,    0 ,    0 , &
                    0 ,    0 ,    0 ,    0 ,    0 ,      0 ,    0 ,    0 ,    0 ,    0], &
                shape(cvroad_apt))

   ! volumetric heat capacity of concrete pavement [J/m3*K]
   real(r8), parameter, dimension(10,5)  :: cvroad_cct &
      = reshape([ 858 ,  858 ,  858 ,  858 ,  858 ,    858 ,  858 ,  858 ,  858 ,  858 , &
                  858 ,  858 ,  858 ,  858 ,  858 ,    858 ,  858 ,  858 ,  858 ,  820 , &
                  858 ,  843 ,  825 ,  820 ,  875 ,    843 ,  820 ,  820 ,  875 ,    0 , &
                  820 ,    0 ,    0 ,    0 ,    0 ,      0 ,    0 ,    0 ,    0 ,    0 , &
                    0 ,    0 ,    0 ,    0 ,    0 ,      0 ,    0 ,    0 ,    0 ,    0], &
                shape(cvroad_cct))


   ! thermal conductivity of asphalt pavement [W/m*K]
   real(r8), parameter, dimension(10,5)  :: tkroad_apt &
      = reshape([0.90 , 0.90 , 0.90 , 0.90 , 0.90 ,   0.90 , 0.90 , 0.90 , 0.90 , 0.90 , &
                 1.05 , 1.05 , 1.05 , 1.05 , 0.90 ,   1.30 , 1.30 , 1.05 , 1.20 , 1.20 , &
                 1.20 , 1.20 , 1.20 , 1.20 , 1.20 ,   1.20 , 1.20 , 1.20 , 1.30 , 1.30 , &
                1.30 ,  1.30 , 1.30 ,    0 ,    0 ,      0 ,    0 ,    0 ,    0 ,    0 , &
                    0 ,    0 ,    0 ,    0 ,    0 ,      0 ,    0 ,    0 ,    0 ,    0], &
                shape(tkroad_apt))

   ! thermal conductivity of concrete pavement [W/m*K]
   real(r8), parameter, dimension(10,5)  :: tkroad_cct &
      = reshape([1.11 , 1.11 , 1.11 , 1.11 , 1.11 ,   1.11 , 1.11 , 1.11 , 1.11 , 1.11 , &
                 1.11 , 1.11 , 1.11 , 1.11 , 1.11 ,   1.11 , 1.11 , 1.11 , 1.11 , 1.30 , &
                 1.11 , 0.92 , 1.20 , 1.30 , 1.45 ,   0.92 , 1.30 , 1.30 , 1.45 ,    0 , &
                 1.30 ,    0 ,    0 ,    0 ,    0 ,      0 ,    0 ,    0 ,    0 ,    0 , &
                    0 ,    0 ,    0 ,    0 ,    0 ,      0 ,    0 ,    0 ,    0 ,    0], &
                shape(tkroad_cct))

END MODULE MOD_Road_Const_ThermalParameters
