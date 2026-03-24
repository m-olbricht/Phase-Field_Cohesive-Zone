!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! interface energy module
!
! Stephan Roth, TU Bergakademie Freiberg, 30.07.2020
!
! 18.01.2024: interface energy module
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE InterfaceEnergyModule

  IMPLICIT NONE

  CONTAINS

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION IED_AT2(D,phase,grad_phase,npar,par)
    ! Helmholtz energy density of interface

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D, npar
      REAL(kind=AbqRK), INTENT(IN) :: phase, grad_phase(D), par(npar)
      INTEGER(kind=AbqIK) :: i1
      REAL(kind=AbqRK) :: prop_EpsGamma, prop_l, gradSquared

      prop_EpsGamma = par(1)
      prop_l        = par(2)

      gradSquared = zero
      DO i1=1,D
        gradSquared = gradSquared + grad_phase(i1)**2 
      END DO

      IED_AT2 = prop_EpsGamma*(one/two/prop_l*phase**2 + prop_l/two*gradSquared)

    END FUNCTION IED_AT2

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_IED_AT2_d_damage(D,phase,grad_phase,npar,par)
    ! sensitivity of Helmholtz energy density of interface w.r.t. order parameter

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D, npar
      REAL(kind=AbqRK), INTENT(IN) :: phase, grad_phase(D), par(npar)
      REAL(kind=AbqRK) :: prop_EpsGamma, prop_l

      prop_EpsGamma = par(1)
      prop_l        = par(2)

      d_IED_AT2_d_damage = prop_EpsGamma/prop_l * phase

    END FUNCTION d_IED_AT2_d_damage

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_IED_AT2_d_damage_d_damage(D,phase,grad_phase,npar,par)
    ! sensitivity of Helmholtz energy density of interface w.r.t. order parameter (2nd derivative)

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D, npar
      REAL(kind=AbqRK), INTENT(IN) :: phase, grad_phase(D), par(npar)
      REAL(kind=AbqRK) :: prop_EpsGamma, prop_l

      prop_EpsGamma = par(1)
      prop_l        = par(2)

      d_IED_AT2_d_damage_d_damage = prop_EpsGamma/prop_l

    END FUNCTION d_IED_AT2_d_damage_d_damage

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_IED_AT2_d_damage_d_grad_damage(D,phase,grad_phase,npar,par)
    ! mixed derivative of Helmholtz energy density of interface w.r.t. order parameter and gradient of order parameter

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D, npar
      REAL(kind=AbqRK), INTENT(IN) :: phase, grad_phase(D), par(npar)
      DIMENSION d_IED_AT2_d_damage_d_grad_damage(D)

      d_IED_AT2_d_damage_d_grad_damage = zero

    END FUNCTION d_IED_AT2_d_damage_d_grad_damage

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_IED_AT2_d_grad_damage(D,phase,grad_phase,npar,par)
    ! sensitivity of Helmholtz energy density of interface w.r.t. gradient of order parameter

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D, npar
      REAL(kind=AbqRK), INTENT(IN) :: phase, grad_phase(D), par(npar)
      REAL(kind=AbqRK) :: prop_EpsGamma, prop_l
      DIMENSION d_IED_AT2_d_grad_damage(D)

      prop_EpsGamma = par(1)
      prop_l        = par(2)

      d_IED_AT2_d_grad_damage = prop_EpsGamma*prop_l*grad_phase

    END FUNCTION d_IED_AT2_d_grad_damage

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_IED_AT2_d_grad_damage_d_grad_damage(D,phase,grad_phase,npar,par)
    ! sensitivity of Helmholtz energy density of interface w.r.t. gradient of order parameter (2nd derivative)

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D, npar
      REAL(kind=AbqRK), INTENT(IN) :: phase, grad_phase(D), par(npar)
      INTEGER(kind=AbqIK) :: i1
      REAL(kind=AbqRK) :: prop_EpsGamma, prop_l
      DIMENSION d_IED_AT2_d_grad_damage_d_grad_damage(D,D)

      prop_EpsGamma = par(1)
      prop_l        = par(2)

      d_IED_AT2_d_grad_damage_d_grad_damage = zero
      FORALL (i1=1:D)
        d_IED_AT2_d_grad_damage_d_grad_damage(i1,i1) = prop_EpsGamma*prop_l
      END FORALL

    END FUNCTION d_IED_AT2_d_grad_damage_d_grad_damage

!------------------------------------------------------------------------------------

END MODULE InterfaceEnergyModule

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

