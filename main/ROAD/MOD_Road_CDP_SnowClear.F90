#include <define.h>

MODULE MOD_Road_CDP_SnowClear

   USE MOD_Precision
   IMPLICIT NONE
   SAVE

! PUBLIC MEMBER FUNCTIONS:
   PUBLIC :: snow_clear_CDP

!-----------------------------------------------------------------------

CONTAINS

!-----------------------------------------------------------------------


SUBROUTINE snow_clear_CDP(data_sc,snow_clear_flag)

    INTEGER,DIMENSION(4)::data_sc
    LOGICAL::snow_clear_flag

    snow_clear_flag=.FALSE.

    !ludo deneigements en fin de periode hors episodes
    IF(data_sc(4).EQ.6) THEN
        snow_clear_flag=.TRUE.
    ENDIF
    !ludo cas fonte nat des que plus de neige observee->deneigement le lendemain a 6h
    IF(      ( (data_sc(1).EQ.1998).AND.(data_sc(2).EQ.10).AND.(data_sc(3).EQ.27) )&
            .OR.( (data_sc(1).EQ.1998).AND.(data_sc(2).EQ.10).AND.(data_sc(3).EQ.30) )&
            .OR.( (data_sc(1).EQ.1998).AND.(data_sc(2).EQ.11).AND.(data_sc(3).EQ.5) )&
            .OR.( (data_sc(1).EQ.1998).AND.(data_sc(2).EQ.11).AND.(data_sc(3).GE.11).AND.(data_sc(3).LE.17) )&
            .OR.( (data_sc(1).EQ.1998).AND.(data_sc(2).EQ.11).AND.(data_sc(3).GE.26) )&
            .OR.( (data_sc(1).EQ.1998).AND.(data_sc(2).EQ.12).AND.(data_sc(3).LE.2) )&
            .OR.( (data_sc(1).EQ.1998).AND.(data_sc(2).EQ.12).AND.(data_sc(3).GE.6).AND.(data_sc(3).LE.8) )&
            .OR.( (data_sc(1).EQ.1998).AND.(data_sc(2).EQ.12).AND.(data_sc(3).EQ.11) )&
            .OR.( (data_sc(1).EQ.1998).AND.(data_sc(2).EQ.12).AND.(data_sc(3).GE.20).AND.(data_sc(3).LE.22) )&
            .OR.( (data_sc(1).EQ.1998).AND.(data_sc(2).EQ.12).AND.(data_sc(3).GE.25).AND.(data_sc(3).LE.30) )&
            .OR.( (data_sc(1).EQ.1999).AND.(data_sc(2).EQ.1).AND.((data_sc(3).EQ.3).OR.(data_sc(3).EQ.4)) )&
            .OR.( (data_sc(1).EQ.1999).AND.(data_sc(2).EQ.1).AND.(data_sc(3).GE.9).AND.(data_sc(3).LE.14) )&
            .OR.( (data_sc(1).EQ.1999).AND.(data_sc(2).EQ.1).AND.(data_sc(3).EQ.18) )&
            .OR.( (data_sc(1).EQ.1999).AND.(data_sc(2).EQ.1).AND.(data_sc(3).GE.26).AND.(data_sc(3).LE.31) )&
            .OR.( (data_sc(1).EQ.1999).AND.(data_sc(2).EQ.2).AND.(data_sc(3).LE.24) )&
            .OR.( (data_sc(1).EQ.1999).AND.(data_sc(2).EQ.3).AND.(data_sc(3).EQ.4) )&
            .OR.( (data_sc(1).EQ.1999).AND.(data_sc(2).EQ.3).AND.(data_sc(3).GE.6).AND.(data_sc(3).LE.10) )&
            .OR.( (data_sc(1).EQ.1999).AND.(data_sc(2).EQ.3).AND.(data_sc(3).GE.22).AND.(data_sc(3).LE.25) )&
            .OR.( (data_sc(1).EQ.1999).AND.(data_sc(2).EQ.3).AND.(data_sc(3).GE.26) )&
            .OR.( (data_sc(1).EQ.1999).AND.(data_sc(2).EQ.4).AND.(data_sc(3).LE.1) )&
            .OR.( (data_sc(1).EQ.1999).AND.(data_sc(2).EQ.4).AND.((data_sc(3).EQ.8).OR.(data_sc(3).EQ.9)) )&
            .OR.( (data_sc(1).EQ.1999).AND.(data_sc(2).EQ.4).AND.(data_sc(3).GE.13).AND.(data_sc(3).LE.16) )&
            .OR.( (data_sc(1).EQ.1999).AND.(data_sc(2).EQ.4).AND.(data_sc(3).GE.18).AND.(data_sc(3).LE.19) )&
            .OR.( (data_sc(1).EQ.1999).AND.(data_sc(2).EQ.4).AND.(data_sc(3).EQ.24) )   )THEN
                
        snow_clear_flag=.FALSE.
    ENDIF

    !ludo : deneigements reels a debut de lheure col de porte
    !DENEIGEMENT SAISON 1998 1999 et 2000
    IF (    (    (data_sc(1).EQ.1998).AND.(data_sc(2).EQ.11).AND.(data_sc(3).EQ.13).AND.(data_sc(4).EQ.12))&
            .OR.((data_sc(1).EQ.1998).AND.(data_sc(2).EQ.11).AND.(data_sc(3).EQ.16).AND.(data_sc(4).EQ.11))&
            .OR.((data_sc(1).EQ.1998).AND.(data_sc(2).EQ.11).AND.(data_sc(3).EQ.17).AND.(data_sc(4).EQ.17))&
            .OR.((data_sc(1).EQ.1998).AND.(data_sc(2).EQ.11).AND.(data_sc(3).EQ.26).AND.(data_sc(4).EQ.15))&
            .OR.((data_sc(1).EQ.1998).AND.(data_sc(2).EQ.12).AND.(data_sc(3).EQ.2).AND.(data_sc(4).EQ.10))&
            .OR.((data_sc(1).EQ.1998).AND.(data_sc(2).EQ.12).AND.(data_sc(3).EQ.8).AND.(data_sc(4).EQ.15))&
            .OR.((data_sc(1).EQ.1998).AND.(data_sc(2).EQ.12).AND.(data_sc(3).EQ.11).AND.(data_sc(4).EQ.15))&
            .OR.((data_sc(1).EQ.1998).AND.(data_sc(2).EQ.12).AND.(data_sc(3).EQ.22).AND.(data_sc(4).EQ.15))&
            .OR.((data_sc(1).EQ.1998).AND.(data_sc(2).EQ.12).AND.(data_sc(3).EQ.30).AND.(data_sc(4).EQ.16))&
            .OR.((data_sc(1).EQ.1999).AND.(data_sc(2).EQ.1).AND.(data_sc(3).EQ.4).AND.(data_sc(4).EQ.16))&
            .OR.((data_sc(1).EQ.1999).AND.(data_sc(2).EQ.1).AND.(data_sc(3).EQ.14).AND.(data_sc(4).EQ.11))&
            .OR.((data_sc(1).EQ.1999).AND.(data_sc(2).EQ.1).AND.(data_sc(3).EQ.18).AND.(data_sc(4).EQ.15))&
            .OR.((data_sc(1).EQ.1999).AND.(data_sc(2).EQ.1).AND.(data_sc(3).EQ.27).AND.(data_sc(4).EQ.15))&
            .OR.((data_sc(1).EQ.1999).AND.(data_sc(2).EQ.1).AND.(data_sc(3).EQ.29).AND.(data_sc(4).EQ.12))&
            .OR.((data_sc(1).EQ.1999).AND.(data_sc(2).EQ.2).AND.(data_sc(3).EQ.4).AND.(data_sc(4).EQ.15))&
            .OR.((data_sc(1).EQ.1999).AND.(data_sc(2).EQ.2).AND.(data_sc(3).EQ.5).AND.(data_sc(4).EQ.15))&
            .OR.((data_sc(1).EQ.1999).AND.(data_sc(2).EQ.2).AND.(data_sc(3).EQ.24).AND.(data_sc(4).EQ.20))&
            .OR.((data_sc(1).EQ.1999).AND.(data_sc(2).EQ.3).AND.(data_sc(3).EQ.4).AND.(data_sc(4).EQ.15))&
            .OR.((data_sc(1).EQ.1999).AND.(data_sc(2).EQ.4).AND.(data_sc(3).EQ.16).AND.(data_sc(4).EQ.14)) ) THEN  
                
        snow_clear_flag=.TRUE.          
                
    ENDIF
        END SUBROUTINE snow_clear_CDP

END MODULE MOD_Road_CDP_SnowClear
! --------- EOP ----------
