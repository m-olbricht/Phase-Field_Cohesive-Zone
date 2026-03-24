!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! degradation energy module
!
! Stephan Roth, TU Bergakademie Freiberg, 17.01.2024
!
! 17.01.2024: degradation functions
! 05.04.2024: quadratic degradation function added
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE DegradationModule

  IMPLICIT NONE

!  PUBLIC :: 

  CONTAINS

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION DegradationCubic(damage,prop_df_beta,prop_df_alpha)
    ! cubic degradation function: Noll, Kuhn, Mueller (2017)
    ! here: 1-g(d)

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: damage, prop_df_beta, prop_df_alpha
      !
      DegradationCubic = one - (prop_df_beta*(damage**3-damage**2) + three*damage**2 - two*damage**3) + prop_df_alpha

    END FUNCTION DegradationCubic

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_DegradationCubic_d_damage(damage,prop_df_beta,prop_df_alpha)
    ! first derivative of cubic degradation function: Noll, Kuhn, Mueller (2017)
    ! here: -g'(d)

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: damage, prop_df_beta, prop_df_alpha
      !
      d_DegradationCubic_d_damage = -(prop_df_beta*(three*damage**2-two*damage) + six*damage - six*damage**2)

    END FUNCTION d_DegradationCubic_d_damage

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_DegradationCubic_d_damage_d_damage(damage,prop_df_beta,prop_df_alpha)
    ! second derivative of cubic degradation function: Noll, Kuhn, Mueller (2017)
    ! here: -g''(d)

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: damage, prop_df_beta, prop_df_alpha
      !
      d_DegradationCubic_d_damage_d_damage = -(prop_df_beta*(six*damage**2-two) + six - twelve*damage)

    END FUNCTION d_DegradationCubic_d_damage_d_damage

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION DegradationQuadratic(damage)
    ! quadratic degradation function: Miehe (2010)

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: damage
      !
      DegradationQuadratic = (one-damage)**2

    END FUNCTION DegradationQuadratic

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_DegradationQuadratic_d_damage(damage)
    ! first derivative of quadratic degradation function: Miehe (2010)

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: damage
      !
      d_DegradationQuadratic_d_damage = two*(damage-one)

    END FUNCTION d_DegradationQuadratic_d_damage

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_DegradationQuadratic_d_damage_d_damage(damage)
    ! second derivative of quadratic degradation function: Miehe (2010)

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: damage
      !
      d_DegradationQuadratic_d_damage_d_damage = two

    END FUNCTION d_DegradationQuadratic_d_damage_d_damage

!------------------------------------------------------------------------------------

END MODULE DegradationModule

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

