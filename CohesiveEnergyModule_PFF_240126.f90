!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! cohesive energy module
!
! Stephan Roth, TU Bergakademie Freiberg, 08.01.2024
!
! 08.01.2024: linear cohesive traction-separation law
! 12.01.2024: cubic degradation function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CohesiveEnergyModule

  ! MacAulay-Operator
  USE MathModul
  
  USE AliasModuleCZ
  ! cubic degradation function: Noll, Kuhn, Mueller (2017)
!~   USE DegradationModule, ONLY: Degradation => DegradationCubic, &
!~                                               d_Degradation_d_damage => d_DegradationCubic_d_damage, &
!~                                               d_Degradation_d_damage_d_damage => d_DegradationCubic_d_damage_d_damage 
  IMPLICIT NONE

  CONTAINS

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION separation(D,normalDirectionIndex,sep)
    ! separation (no negativ normal separation)

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D, normalDirectionIndex
      REAL(kind=AbqRK), INTENT(IN) :: sep(D)
      DIMENSION separation(D)
      !
      separation    = sep
      ! positive normal separation
      separation(normalDirectionIndex) = MacAulay(sep(normalDirectionIndex))
 
    END FUNCTION separation

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_separation_d_sep(D,normalDirectionIndex,sep)
    ! first derivative of separation (no negativ normal separation)

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D, normalDirectionIndex
      REAL(kind=AbqRK), INTENT(IN) :: sep(D)
      INTEGER(kind=AbqIK) :: i1
      DIMENSION d_separation_d_sep(D,D)
      !
      d_separation_d_sep = zero
      FORALL (i1=1:D)
        d_separation_d_sep(i1,i1) = one
      END FORALL
      d_separation_d_sep(normalDirectionIndex,normalDirectionIndex) = d_MacAulay(sep(normalDirectionIndex))
      
    END FUNCTION d_separation_d_sep

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION effSep(D,normalDirectionIndex,sep,prop_delta0)
    ! normalized effective separation, no contribution from negative normal separation

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D, normalDirectionIndex
      REAL(kind=AbqRK), INTENT(IN) :: sep(D), prop_delta0
      INTEGER(kind=AbqIK) :: i1
      REAL(kind=AbqRK) :: sepPos(D)
      !
      sepPos = separation(D,normalDirectionIndex,sep)
      !
      effSep = zero
      DO i1=1,D
        effSep = effSep + sepPos(i1)**two
      END DO
      effSep = sqrt(effSep) / prop_delta0
 
    END FUNCTION effSep

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_effSep_d_sep(D,normalDirectionIndex,sep,prop_delta0)
    ! normalized effective separation: first derivative w.r.t. separation coordinates

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D, normalDirectionIndex
      REAL(kind=AbqRK), INTENT(IN) :: sep(D), prop_delta0
      REAL(kind=AbqRK) :: lambda
      REAL(kind=AbqRK) :: sepPos(D)
      DIMENSION d_effSep_d_sep(D)
      !
      sepPos = separation(D,normalDirectionIndex,sep)
      !
      lambda = effSep(D,normalDirectionIndex,sep,prop_delta0)
      !
      d_effSep_d_sep = zero
      IF (lambda .GT. zero) THEN
        d_effSep_d_sep = sepPos / prop_delta0**2 / lambda
      END IF
 
    END FUNCTION d_effSep_d_sep

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_effSep_d_sep_d_sep(D,normalDirectionIndex,sep,prop_delta0)
    ! normalized effective separation: second derivative w.r.t. separation coordinates

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D, normalDirectionIndex
      REAL(kind=AbqRK), INTENT(IN) :: sep(D), prop_delta0
      REAL(kind=AbqRK) :: lambda
      INTEGER(kind=AbqIK) :: i1, i2
      REAL(kind=AbqRK) :: sepPos(D)
      DIMENSION d_effSep_d_sep_d_sep(D,D)
      !
      sepPos = separation(D,normalDirectionIndex,sep)
      !
      lambda = effSep(D,normalDirectionIndex,sep,prop_delta0)
      !
      d_effSep_d_sep_d_sep = zero
      IF (lambda .GT. zero) THEN
        FORALL (i1=1:D,i2=1:D)
          d_effSep_d_sep_d_sep(i1,i2) = sepPos(i1)*sepPos(i2)
        END FORALL
        !
        d_effSep_d_sep_d_sep = (d_separation_d_sep(D,normalDirectionIndex,sep) - d_effSep_d_sep_d_sep / prop_delta0**2 / lambda**2) / prop_delta0**2 / lambda
      END IF

    END FUNCTION d_effSep_d_sep_d_sep

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION CZED(D,normalDirectionIndex,sep,damage,prop_delta0,prop_t0,prop_df_one,prop_df_two)
    ! cohesive zone energy density

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D, normalDirectionIndex
      REAL(kind=AbqRK), INTENT(IN) :: sep(D), damage, prop_delta0, prop_t0, prop_df_one, prop_df_two
      REAL(kind=AbqRK) :: lambda
      !
      lambda = effSep(D,normalDirectionIndex,sep,prop_delta0)
      !
      CZED = prop_t0*prop_delta0/two*lambda**2 * Degradation(damage,prop_df_one,prop_df_two)

    END FUNCTION CZED

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_CZED_d_sep(D,normalDirectionIndex,damage,sep,prop_delta0,prop_t0,prop_df_one,prop_df_two)
    ! cohesive zone energy density: first derivative w.r.t. separation coordinates

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D, normalDirectionIndex
      REAL(kind=AbqRK), INTENT(IN) :: sep(D), damage, prop_delta0, prop_t0, prop_df_one, prop_df_two
      REAL(kind=AbqRK) :: lambda
      DIMENSION d_CZED_d_sep(D)
      !
      lambda = effSep(D,normalDirectionIndex,sep,prop_delta0)
      !
      d_CZED_d_sep = prop_t0*prop_delta0*lambda * d_effSep_d_sep(D,normalDirectionIndex,sep,prop_delta0) * Degradation(damage,prop_df_one,prop_df_two)

    END FUNCTION d_CZED_d_sep

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_CZED_d_sep_d_sep(D,normalDirectionIndex,damage,sep,prop_delta0,prop_t0,prop_df_one,prop_df_two)
    ! cohesive zone energy density: second derivative w.r.t. separation coordinates

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D, normalDirectionIndex
      REAL(kind=AbqRK), INTENT(IN) :: sep(D), damage, prop_delta0, prop_t0, prop_df_one, prop_df_two
      INTEGER(kind=AbqIK) :: i1, i2
      REAL(kind=AbqRK) :: lambda, d_lambda_d_sep(D)
      DIMENSION d_CZED_d_sep_d_sep(D,D)
      !
      lambda = effSep(D,normalDirectionIndex,sep,prop_delta0)
      !
      d_lambda_d_sep = d_effSep_d_sep(D,normalDirectionIndex,sep,prop_delta0)
      !
      d_CZED_d_sep_d_sep = zero
      FORALL (i1=1:D,i2=1:D)
        d_CZED_d_sep_d_sep(i1,i2) = d_lambda_d_sep(i1)*d_lambda_d_sep(i2)
      END FORALL
      !
      d_CZED_d_sep_d_sep = prop_t0*prop_delta0 * (d_CZED_d_sep_d_sep + lambda*d_effSep_d_sep_d_sep(D,normalDirectionIndex,sep,prop_delta0)) &
                                               * Degradation(damage,prop_df_one,prop_df_two)

    END FUNCTION d_CZED_d_sep_d_sep

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_CZED_d_sep_d_damage(D,normalDirectionIndex,damage,sep,prop_delta0,prop_t0,prop_df_one,prop_df_two)
    ! cohesive zone energy density: mixed second derivative w.r.t. separation coordinates and damage

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D, normalDirectionIndex
      REAL(kind=AbqRK), INTENT(IN) :: sep(D), damage, prop_delta0, prop_t0, prop_df_one, prop_df_two
      REAL(kind=AbqRK) :: lambda
      DIMENSION d_CZED_d_sep_d_damage(D)
      !
      lambda = effSep(D,normalDirectionIndex,sep,prop_delta0)
      !
      d_CZED_d_sep_d_damage = prop_t0*prop_delta0*lambda * d_effSep_d_sep(D,normalDirectionIndex,sep,prop_delta0) &
                                                         * d_Degradation_d_damage(damage,prop_df_one,prop_df_two)

    END FUNCTION d_CZED_d_sep_d_damage

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_CZED_d_damage(D,normalDirectionIndex,damage,sep,prop_delta0,prop_t0,prop_df_one,prop_df_two)
    ! cohesive zone energy density: first derivative w.r.t. damage variable

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D, normalDirectionIndex
      REAL(kind=AbqRK), INTENT(IN) :: sep(D), damage, prop_delta0, prop_t0, prop_df_one, prop_df_two
      REAL(kind=AbqRK) :: lambda
      !
      lambda = effSep(D,normalDirectionIndex,sep,prop_delta0)
      !
      d_CZED_d_damage = prop_t0*prop_delta0/two*lambda**2 * d_Degradation_d_damage(damage,prop_df_one,prop_df_two)

    END FUNCTION d_CZED_d_damage

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_CZED_d_damage_d_damage(D,normalDirectionIndex,damage,sep,prop_delta0,prop_t0,prop_df_one,prop_df_two)
    ! cohesive zone energy density: second derivative w.r.t. damage variable

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D, normalDirectionIndex
      REAL(kind=AbqRK), INTENT(IN) :: sep(D), damage, prop_delta0, prop_t0, prop_df_one, prop_df_two
      REAL(kind=AbqRK) :: lambda
      !
      lambda = effSep(D,normalDirectionIndex,sep,prop_delta0)
      !
      d_CZED_d_damage_d_damage = prop_t0*prop_delta0/two*lambda**2 * d_Degradation_d_damage_d_damage(damage,prop_df_one,prop_df_two)

    END FUNCTION d_CZED_d_damage_d_damage

!------------------------------------------------------------------------------------

END MODULE CohesiveEnergyModule

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

