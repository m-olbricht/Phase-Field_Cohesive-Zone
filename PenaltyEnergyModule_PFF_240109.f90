!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! penalty energy module
!
! Stephan Roth, TU Bergakademie Freiberg, 08.01.2024
!
! 08.01.2024: penalty energy to prevent damage jump
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PenaltyEnergyModule

  IMPLICIT NONE

  CONTAINS

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION PenED_Exponential(damageJump,prop_pen)
    ! penalty energy density, exponential

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: damageJump, prop_pen
      !
      PenED_Exponential = prop_pen*( exp(damageJump) - damageJump - one )

    END FUNCTION PenED_Exponential

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_PenED_d_damageJump_Exponential(damageJump,prop_pen)
    ! penalty energy density, exponential: first derivative w.r.t. damage jump

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: damageJump, prop_pen
      !
      d_PenED_d_damageJump_Exponential =  prop_pen*( exp(damageJump) - one )

    END FUNCTION d_PenED_d_damageJump_Exponential

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_PenED_d_damageJump_d_damageJump_Exponential(damageJump,prop_pen)
    ! penalty energy density, exponential: second derivative w.r.t. damage jump

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: damageJump, prop_pen
      !
      d_PenED_d_damageJump_d_damageJump_Exponential =  prop_pen* exp(damageJump)

    END FUNCTION d_PenED_d_damageJump_d_damageJump_Exponential

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION PenED_PowerLaw(damageJump,prop_eta,prop_penexp)
    ! penalty energy density, power law

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: damageJump, prop_eta, prop_penexp
      !
      PenED_PowerLaw = prop_eta/prop_penexp*damageJump**prop_penexp

    END FUNCTION PenED_PowerLaw

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_PenED_d_damageJump_PowerLaw(damageJump,prop_eta,prop_penexp)
    ! penalty energy density, power law: first derivative w.r.t. damage jump

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: damageJump, prop_eta, prop_penexp
      !
      d_PenED_d_damageJump_PowerLaw = prop_eta*damageJump**(prop_penexp-one)

    END FUNCTION d_PenED_d_damageJump_PowerLaw

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_PenED_d_damageJump_d_damageJump_PowerLaw(damageJump,prop_eta,prop_penexp)
    ! penalty energy density, power law: second derivative w.r.t. damage jump

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: damageJump, prop_eta, prop_penexp
      !
      d_PenED_d_damageJump_d_damageJump_PowerLaw = (prop_penexp-one)*prop_eta*damageJump**(prop_penexp-two)

    END FUNCTION d_PenED_d_damageJump_d_damageJump_PowerLaw

!------------------------------------------------------------------------------------

END MODULE PenaltyEnergyModule

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

