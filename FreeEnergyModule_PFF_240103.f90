!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! free energy module
!
! Stephan Roth, TU Bergakademie Freiberg, 30.07.2020
!
! 03.01.2014: isotropic, linear-elastic Hooke's law
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE FreeEnergyModule

  IMPLICIT NONE

!  PUBLIC :: 

  CONTAINS

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION eps_e(eps)
    ! elastic strain

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      DIMENSION eps_e(3,3)
      !
      eps_e = eps
 
    END FUNCTION eps_e

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_eps_e_d_eps(eps)
    ! sensitivity of elastic strain w.r.t. total strain

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      INTEGER(kind=AbqIK) :: i1, i2
      DIMENSION d_eps_e_d_eps(3,3,3,3)
      !
      d_eps_e_d_eps = zero
      FORALL (i1=1:3,i2=1:3)
        d_eps_e_d_eps(i1,i2,i1,i2) = one
      END FORALL
 
    END FUNCTION d_eps_e_d_eps

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION HFED_Mix(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Helmholtz free energy density

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK) :: eps_e_i(3,3), ETensor(3,3,3,3)
      REAL(kind=AbqRK) :: E, nu
      !
      E               = parHFEDMatrixPhase(1)
      nu              = parHFEDMatrixPhase(2)
      !
      eps_e_i = eps_e(eps)
      !
      ETensor=ElastTensor(E,nu)
      !
      HFED_Mix = one/two*DOUBLECONTRACTIONSTwoFourTwo(eps_e_i,ETensor,eps_e_i)

    END FUNCTION HFED_Mix

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_HFED_Mix_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Helmholtz free energy density

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK) :: eps_e_i(3,3), ETensor(3,3,3,3)
      REAL(kind=AbqRK) :: E, nu
      DIMENSION d_HFED_Mix_d_eps_e(3,3)
      !
      E               = parHFEDMatrixPhase(1)
      nu              = parHFEDMatrixPhase(2)
      !
      eps_e_i = eps_e(eps)
      !
      ETensor=ElastTensor(E,nu)
      !
      d_HFED_Mix_d_eps_e = DOUBLECONTRACTIONTwoFour(eps_e_i,ETensor)

    END FUNCTION d_HFED_Mix_d_eps_e

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_HFED_Mix_d_eps_e_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Helmholtz free energy density

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK) :: E, nu
      DIMENSION d_HFED_Mix_d_eps_e_d_eps_e(3,3,3,3)
      !
      E          = parHFEDMatrixPhase(1)
      nu         = parHFEDMatrixPhase(2)
      !
      d_HFED_Mix_d_eps_e_d_eps_e = ElastTensor(E,nu)

    END FUNCTION d_HFED_Mix_d_eps_e_d_eps_e

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_HFED_Mix_d_eps(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Helmholtz free energy density

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK) :: d_eps_e_alpha_d_eps(3,3,3,3)
      DIMENSION d_HFED_Mix_d_eps(3,3)
      !
      d_eps_e_alpha_d_eps = d_eps_e_d_eps(eps)
      !
      d_HFED_Mix_d_eps = DOUBLECONTRACTIONTwoFour(d_HFED_Mix_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase),d_eps_e_alpha_d_eps)

    END FUNCTION d_HFED_Mix_d_eps

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_HFED_Mix_d_eps_d_eps(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Helmholtz free energy density

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE TensorModule

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK) :: d_eps_e_alpha_d_eps(3,3,3,3)
      DIMENSION d_HFED_Mix_d_eps_d_eps(3,3,3,3)
      !
      d_eps_e_alpha_d_eps = d_eps_e_d_eps(eps)
      !
      d_HFED_Mix_d_eps_d_eps = DOUBLECONTRACTIONFourFour(DOUBLECONTRACTIONFourFour(TransposedForthOrderTensor(d_eps_e_alpha_d_eps), & 
                                                         d_HFED_Mix_d_eps_e_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)),d_eps_e_alpha_d_eps)

    END FUNCTION d_HFED_Mix_d_eps_d_eps

!------------------------------------------------------------------------------------

END MODULE FreeEnergyModule

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

