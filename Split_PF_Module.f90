!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Split Module of the Helmholtz Free Energy
!
! Martin Olbricht, TU Bergakademie Freiberg, 09.07.2024
!
! 09.07.2024 Isotropic No Split
!
! 26.07.2024 Amor Split
!
! /// noch nicht 19.03.2025 Miehe Split
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE SplitModulePF

  USE TensorModule
  USE ABQINTERFACE_PF
  USE FLOATNUMBERS

  IMPLICIT NONE


  CONTAINS

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION eps_e(eps)
    ! elastic strain

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      DIMENSION eps_e(3,3)
      !
      eps_e = eps
 
    END FUNCTION eps_e

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION HFEDposNoSplit(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Psi pos no split algorithm

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
      !
      HFEDposNoSplit = one/two*DOUBLECONTRACTIONSTwoFourTwo(eps_e_i,ETensor,eps_e_i)
      
    END FUNCTION HFEDposNoSplit

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION HFEDnegNoSplit(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Psi neg no split algorithm

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK) :: eps_e_i(3,3), ETensor(3,3,3,3)
      REAL(kind=AbqRK) :: E, nu
      !
      HFEDnegNoSplit = zero
 
    END FUNCTION HFEDnegNoSplit

!------------------------------------------------------------------------------------
!
!
! 1. Ableitung
!
!
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_HFEDposNoSplit_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Psi pos no split algorithm

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK) :: eps_e_i(3,3), ETensor(3,3,3,3)
      REAL(kind=AbqRK) :: E, nu
      DIMENSION d_HFEDposNoSplit_d_eps_e(3,3)
      !
      E               = parHFEDMatrixPhase(1)
      nu              = parHFEDMatrixPhase(2)
      !
      eps_e_i = eps_e(eps)
      !
      ETensor=ElastTensor(E,nu)
      !
      !
      d_HFEDposNoSplit_d_eps_e = DOUBLECONTRACTIONTwoFour(eps_e_i,ETensor)
 
    END FUNCTION d_HFEDposNoSplit_d_eps_e

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_HFEDnegNoSplit_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Psi neg no split algorithm

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK) :: eps_e_i(3,3), ETensor(3,3,3,3)
      REAL(kind=AbqRK) :: E, nu
      DIMENSION d_HFEDnegNoSplit_d_eps_e(3,3)
      !
      d_HFEDnegNoSplit_d_eps_e = zero
 
    END FUNCTION d_HFEDnegNoSplit_d_eps_e

!------------------------------------------------------------------------------------
!
!
! 2. Ableitung
!
!
!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_HFEDposNoSplit_d_eps_e_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Psi neg no split algorithm

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK) :: eps_e_i(3,3), ETensor(3,3,3,3)
      REAL(kind=AbqRK) :: E, nu
      DIMENSION d_HFEDposNoSplit_d_eps_e_d_eps_e(3,3,3,3)
      !
      E          = parHFEDMatrixPhase(1)
      nu         = parHFEDMatrixPhase(2)
      !
      d_HFEDposNoSplit_d_eps_e_d_eps_e = ElastTensor(E,nu)
 
    END FUNCTION d_HFEDposNoSplit_d_eps_e_d_eps_e

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_HFEDnegNoSplit_d_eps_e_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
    ! Psi neg no split algorithm

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK) :: eps_e_i(3,3), ETensor(3,3,3,3)
      REAL(kind=AbqRK) :: E, nu
      DIMENSION d_HFEDnegNoSplit_d_eps_e_d_eps_e(3,3,3,3)
      !
      d_HFEDnegNoSplit_d_eps_e_d_eps_e = zero
 
    END FUNCTION d_HFEDnegNoSplit_d_eps_e_d_eps_e

!~ !------------------------------------------------------------------------------------
!~ !
!~ !
!~ !
!~ ! Amor et al
!~ !
!~ !
!~ !
!~ !------------------------------------------------------------------------------------

!~     REAL(kind=AbqRK) FUNCTION bulkModulus(nHFEDpar,parHFEDMatrixPhase)
!~     ! Psi pos Amor et al (2009)

!~       USE ABQINTERFACE_PF
!~       USE FLOATNUMBERS
!~       USE TensorModule

!~       IMPLICIT NONE
!~       INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
!~       REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
!~       REAL(kind=AbqRK) :: E, nu, lambda, mu
!~       !
!~       bulkModulus = zero
!~       !
!~       E               = parHFEDMatrixPhase(1)
!~       nu              = parHFEDMatrixPhase(2)
!~       !
!~       lambda          = nu * E / ((one - two * nu)*(one + nu))
!~       mu              = half * one / (one + nu) * E
!~       !
!~       bulkModulus = lambda + two / three * mu
      
      
!~     END FUNCTION bulkModulus

!~ !------------------------------------------------------------------------------------

!~     REAL(kind=AbqRK) FUNCTION HFEDposAmorSplit(eps,nHFEDpar,parHFEDMatrixPhase)
!~     ! Psi pos Amor et al (2009)

!~       USE ABQINTERFACE_PF
!~       USE FLOATNUMBERS
!~       USE MathModul
!~       USE TensorModule

!~       IMPLICIT NONE
!~       INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
!~       REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
!~       REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
!~       REAL(kind=AbqRK) :: eps_e_i(3,3), eps_D(3,3), ETensor(3,3,3,3)
!~       REAL(kind=AbqRK) :: E, nu, K, mu
!~       !
!~       HFEDposAmorSplit = zero
!~       !
!~       E               = parHFEDMatrixPhase(1)
!~       nu              = parHFEDMatrixPhase(2)
!~       !
!~       mu              = half * one / (one + nu) * E
!~       !
!~       K = bulkModulus(nHFEDpar,parHFEDMatrixPhase)
!~       !
!~       !
!~       eps_D = deviator(eps)
!~       !
!~       HFEDposAmorSplit = half * K * MacAulay(trace(eps)) ** two + mu * DOUBLECONTRACTIONTwoTwo(eps_D,eps_D)
      
!~     END FUNCTION HFEDposAmorSplit

!~ !------------------------------------------------------------------------------------

!~     REAL(kind=AbqRK) FUNCTION HFEDnegAmorSplit(eps,nHFEDpar,parHFEDMatrixPhase)
!~     ! Psi neg Amor et al (2009)

!~       USE ABQINTERFACE_PF
!~       USE FLOATNUMBERS
!~       USE MathModul
!~       USE TensorModule

!~       IMPLICIT NONE
!~       INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
!~       REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
!~       REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
!~       REAL(kind=AbqRK) :: eps_e_i(3,3), eps_D(3,3), ETensor(3,3,3,3)
!~       REAL(kind=AbqRK) :: E, nu, K, mu
!~       !
!~       HFEDnegAmorSplit = zero
!~       !
!~       E               = parHFEDMatrixPhase(1)
!~       nu              = parHFEDMatrixPhase(2)
!~       !
!~       mu              = half * one / (one + nu) * E
!~       !
!~       K = bulkModulus(nHFEDpar,parHFEDMatrixPhase)
!~       !
!~       !
!~       eps_D = deviator(eps)
!~       !
!~       HFEDnegAmorSplit = half * K * MacAulay(-trace(eps)) ** two
      
!~     END FUNCTION HFEDnegAmorSplit

!~ !------------------------------------------------------------------------------------
!~ !
!~ !
!~ ! 1. Ableitung
!~ !
!~ !
!~ !------------------------------------------------------------------------------------

!~     REAL(kind=AbqRK) FUNCTION d_HFEDposAmorSplit_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
!~     ! Psi pos Amor et al (2009)

!~       USE ABQINTERFACE_PF
!~       USE FLOATNUMBERS
!~       USE MathModul
!~       USE TensorModule

!~       IMPLICIT NONE
!~       INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
!~       REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
!~       REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
!~       REAL(kind=AbqRK) :: eps_e_i(3,3), eps_D(3,3), ETensor(3,3,3,3), Ident(3,3)
!~       REAL(kind=AbqRK) :: E, nu, K, mu
!~       DIMENSION d_HFEDposAmorSplit_d_eps_e(3,3)
!~       !
!~       d_HFEDposAmorSplit_d_eps_e = zero
!~       !
!~       E               = parHFEDMatrixPhase(1)
!~       nu              = parHFEDMatrixPhase(2)
!~       !
!~       mu              = half * one / (one + nu) * E
!~       !
!~       K = bulkModulus(nHFEDpar,parHFEDMatrixPhase)
!~       !
!~       !
!~       eps_D = deviator(eps)
!~       !
!~       ! Von n abhängig machen wenn Funktion auf 2/3D aktualisiert werden
!~       Ident = Identity(3)
!~       !Ident(3,3) = zero ! Für 2D tests
!~       !
!~       !
!~       d_HFEDposAmorSplit_d_eps_e = K * MacAulay(trace(eps)) * Ident + two * mu * eps_D
!~       !
!~       !
!~     END FUNCTION d_HFEDposAmorSplit_d_eps_e

!~ !------------------------------------------------------------------------------------

!~     REAL(kind=AbqRK) FUNCTION d_HFEDnegAmorSplit_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
!~     ! Psi pos Amor et al (2009)

!~       USE ABQINTERFACE_PF
!~       USE FLOATNUMBERS
!~       USE MathModul
!~       USE TensorModule

!~       IMPLICIT NONE
!~       INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
!~       REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
!~       REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
!~       REAL(kind=AbqRK) :: eps_e_i(3,3), eps_D(3,3), ETensor(3,3,3,3), Ident(3,3)
!~       REAL(kind=AbqRK) :: E, nu, K, mu
!~       DIMENSION d_HFEDnegAmorSplit_d_eps_e(3,3)
!~       !
!~       d_HFEDnegAmorSplit_d_eps_e = zero
!~       !
!~       E               = parHFEDMatrixPhase(1)
!~       nu              = parHFEDMatrixPhase(2)
!~       !
!~       mu              = half * one / (one + nu) * E
!~       !
!~       K = bulkModulus(nHFEDpar,parHFEDMatrixPhase)
!~       !
!~       !
!~       eps_D = deviator(eps)
!~       !
!~       ! Von n abhängig machen wenn Funktion auf 2/3D aktualisiert werden
!~       Ident = Identity(3)
!~       !Ident(3,3) = zero ! Für 2D tests
!~       !
!~       !
!~       d_HFEDnegAmorSplit_d_eps_e = - K * MacAulay(-trace(eps)) * Ident

!~     END FUNCTION d_HFEDnegAmorSplit_d_eps_e

!~ !------------------------------------------------------------------------------------

!~     REAL(kind=AbqRK) FUNCTION d_HFEDposAmorSplit_d_eps_e_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
!~     ! Psi pos Amor et al (2009)

!~       USE ABQINTERFACE_PF
!~       USE FLOATNUMBERS
!~       USE MathModul
!~       USE TensorModule

!~       IMPLICIT NONE
!~       INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
!~       REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
!~       REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
!~       INTEGER(kind=AbqIK) :: i1,i2,i3
!~       REAL(kind=AbqRK) :: eps_e_i(3,3), eps_D(3,3), ETensor(3,3,3,3),P_sph_dummy(3,3,3,3),P_dev_dummy(3,3,3,3),Split_dummy(3,3,3,3),Split_dummy_voigt(6,6),P_sph_dummy_voigt(6,6),P_dev_dummy_voigt(6,6)
!~       REAL(kind=AbqRK) :: E, nu, K, mu
!~       DIMENSION d_HFEDposAmorSplit_d_eps_e_d_eps_e(3,3,3,3)
!~       !
!~       d_HFEDposAmorSplit_d_eps_e_d_eps_e = zero
!~       !
!~       P_dev_dummy = zero
!~       !
!~       E               = parHFEDMatrixPhase(1)
!~       nu              = parHFEDMatrixPhase(2)
!~       !
!~       mu              = half * one / (one + nu) * E
!~       !
!~       K = bulkModulus(nHFEDpar,parHFEDMatrixPhase)
!~       !
!~       !
!~       eps_D = deviator(eps)
!~       !
!~       d_HFEDposAmorSplit_d_eps_e_d_eps_e = d_MacAulay(trace(eps)) * K * three * Projection_sph(1) + two * mu * Projection_dev(1)
!~       !
!~ 	  !
!~     END FUNCTION d_HFEDposAmorSplit_d_eps_e_d_eps_e

!~ !------------------------------------------------------------------------------------

!~     REAL(kind=AbqRK) FUNCTION d_HFEDnegAmorSplit_d_eps_e_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
!~     ! Psi pos Amor et al (2009)

!~       USE ABQINTERFACE_PF
!~       USE FLOATNUMBERS
!~       USE MathModul
!~       USE TensorModule

!~       IMPLICIT NONE
!~       INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
!~       REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
!~       REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
!~       INTEGER(kind=AbqIK) :: i1,i2,i3 ! Debug print
!~       REAL(kind=AbqRK) :: eps_e_i(3,3), eps_D(3,3), ETensor(3,3,3,3), Split_dummy(3,3,3,3), Split_dummy_voigt(6,6) ! dummy voigt
!~       REAL(kind=AbqRK) :: E, nu, K, mu
!~       DIMENSION d_HFEDnegAmorSplit_d_eps_e_d_eps_e(3,3,3,3)
!~       !
!~       d_HFEDnegAmorSplit_d_eps_e_d_eps_e = zero
!~       !
!~       E               = parHFEDMatrixPhase(1)
!~       nu              = parHFEDMatrixPhase(2)
!~       !
!~       mu              = half * one / (one + nu) * E
!~       !
!~       K = bulkModulus(nHFEDpar,parHFEDMatrixPhase)
!~       !
!~       !
!~       eps_D = deviator(eps)
!~       !
!~       d_HFEDnegAmorSplit_d_eps_e_d_eps_e = d_MacAulay(-trace(eps)) * K * three * Projection_sph(1)
!~       !
!~       !
!~ 	  !
!~     END FUNCTION d_HFEDnegAmorSplit_d_eps_e_d_eps_e
    
!~ !------------------------------------------------------------------------------------
!~ !
!~ !
!~ !
!~ ! Star-Split Laura de Lorenzis et al (10.1007/s10704-024-00763-w)
!~ !
!~ !
!~ !
!~ !------------------------------------------------------------------------------------

!~     REAL(kind=AbqRK) FUNCTION HFEDposStarSplit(eps,nHFEDpar,parHFEDMatrixPhase)
!~     ! Psi pos Laura de Lorenzis et al (2024)

!~       USE ABQINTERFACE_PF
!~       USE FLOATNUMBERS
!~       USE MathModul
!~       USE TensorModule

!~       IMPLICIT NONE
!~       INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
!~       REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
!~       REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
!~       REAL(kind=AbqRK) :: eps_e_i(3,3), eps_D(3,3), ETensor(3,3,3,3)
!~       REAL(kind=AbqRK) :: E, nu, K, mu, gammastar
!~       !
!~       HFEDposStarSplit = zero
!~       !
!~       E               = parHFEDMatrixPhase(1)
!~       nu              = parHFEDMatrixPhase(2)
!~       gammastar		  = parHFEDMatrixPhase(3)
!~       !
!~       mu              = half * one / (one + nu) * E
!~       !
!~       K = bulkModulus(nHFEDpar,parHFEDMatrixPhase)
!~       !
!~       !
!~       eps_D = deviator(eps)
!~       !
!~       HFEDposStarSplit = half * K * (MacAulay(trace(eps)) ** two - gammastar * MacAulay_neg(trace(eps)) ** two) + mu * DOUBLECONTRACTIONTwoTwo(eps_D,eps_D)
      
!~     END FUNCTION HFEDposStarSplit

!~ !------------------------------------------------------------------------------------

!~     REAL(kind=AbqRK) FUNCTION HFEDnegStarSplit(eps,nHFEDpar,parHFEDMatrixPhase)
!~     ! Psi pos Laura de Lorenzis et al (2024)

!~       USE ABQINTERFACE_PF
!~       USE FLOATNUMBERS
!~       USE MathModul
!~       USE TensorModule

!~       IMPLICIT NONE
!~       INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
!~       REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
!~       REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
!~       REAL(kind=AbqRK) :: eps_e_i(3,3), eps_D(3,3), ETensor(3,3,3,3)
!~       REAL(kind=AbqRK) :: E, nu, K, mu, gammastar
!~       !
!~       HFEDnegStarSplit = zero
!~       !
!~       E               = parHFEDMatrixPhase(1)
!~       nu              = parHFEDMatrixPhase(2)
!~       gammastar		  = parHFEDMatrixPhase(3)
!~       !
!~       mu              = half * one / (one + nu) * E
!~       !
!~       K = bulkModulus(nHFEDpar,parHFEDMatrixPhase)
!~       !
!~       !
!~       eps_D = deviator(eps)
!~       !
!~       HFEDnegStarSplit = half * K * (one+gammastar) * MacAulay_neg(trace(eps)) ** two
      
!~     END FUNCTION HFEDnegStarSplit

!~ !------------------------------------------------------------------------------------
!~ !
!~ !
!~ ! 1. Ableitung
!~ !
!~ !
!~ !------------------------------------------------------------------------------------

!~     REAL(kind=AbqRK) FUNCTION d_HFEDposStarSplit_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
!~     ! Psi pos Laura de Lorenzis et al (2024)

!~       USE ABQINTERFACE_PF
!~       USE FLOATNUMBERS
!~       USE MathModul
!~       USE TensorModule

!~       IMPLICIT NONE
!~       INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
!~       REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
!~       REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
!~       REAL(kind=AbqRK) :: eps_e_i(3,3), eps_D(3,3), Ident(3,3), ETensor(3,3,3,3)
!~       REAL(kind=AbqRK) :: E, nu, K, mu, gammastar
!~       DIMENSION d_HFEDposStarSplit_d_eps_e(3,3)
!~       !
!~       d_HFEDposStarSplit_d_eps_e = zero
!~       !
!~       E               = parHFEDMatrixPhase(1)
!~       nu              = parHFEDMatrixPhase(2)
!~       gammastar		  = parHFEDMatrixPhase(3)
!~       !
!~       mu              = half * one / (one + nu) * E
!~       !
!~       K = bulkModulus(nHFEDpar,parHFEDMatrixPhase)
!~       !
!~       !
!~       eps_D = deviator(eps)
!~       Ident = Identity(3)
!~       !
!~       d_HFEDposStarSplit_d_eps_e = K * (MacAulay(trace(eps)) - gammastar * MacAulay_neg(trace(eps))) * Ident + two * mu * eps_D
      
!~     END FUNCTION d_HFEDposStarSplit_d_eps_e

!~ !------------------------------------------------------------------------------------

!~     REAL(kind=AbqRK) FUNCTION d_HFEDnegStarSplit_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
!~     ! Psi pos Laura de Lorenzis et al (2024)

!~       USE ABQINTERFACE_PF
!~       USE FLOATNUMBERS
!~       USE MathModul
!~       USE TensorModule

!~       IMPLICIT NONE
!~       INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
!~       REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
!~       REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
!~       REAL(kind=AbqRK) :: eps_e_i(3,3), eps_D(3,3), Ident(3,3), ETensor(3,3,3,3)
!~       REAL(kind=AbqRK) :: E, nu, K, mu, gammastar
!~       DIMENSION d_HFEDnegStarSplit_d_eps_e(3,3)
!~       !
!~       d_HFEDnegStarSplit_d_eps_e = zero
!~       !
!~       E               = parHFEDMatrixPhase(1)
!~       nu              = parHFEDMatrixPhase(2)
!~       gammastar		  = parHFEDMatrixPhase(3)
!~       !
!~       mu              = half * one / (one + nu) * E
!~       !
!~       K = bulkModulus(nHFEDpar,parHFEDMatrixPhase)
!~       !
!~       !
!~       eps_D = deviator(eps)
!~       Ident = Identity(3)
!~       !
!~       d_HFEDnegStarSplit_d_eps_e = K * (one+gammastar) * MacAulay_neg(trace(eps)) * Ident
      
!~     END FUNCTION d_HFEDnegStarSplit_d_eps_e

!~ !------------------------------------------------------------------------------------
!~ !
!~ !
!~ ! 2. Ableitung
!~ !
!~ !
!~ !------------------------------------------------------------------------------------

!~     REAL(kind=AbqRK) FUNCTION d_HFEDposStarSplit_d_eps_e_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
!~     ! Psi pos Laura de Lorenzis et al (2024)

!~       USE ABQINTERFACE_PF
!~       USE FLOATNUMBERS
!~       USE MathModul
!~       USE TensorModule

!~       IMPLICIT NONE
!~       INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
!~       REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
!~       REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
!~       REAL(kind=AbqRK) :: eps_e_i(3,3), eps_D(3,3), Ident(3,3), ETensor(3,3,3,3)
!~       REAL(kind=AbqRK) :: E, nu, K, mu, gammastar
!~       DIMENSION d_HFEDposStarSplit_d_eps_e_d_eps_e(3,3,3,3)
!~       !
!~       d_HFEDposStarSplit_d_eps_e_d_eps_e = zero
!~       !
!~       E               = parHFEDMatrixPhase(1)
!~       nu              = parHFEDMatrixPhase(2)
!~       gammastar		  = parHFEDMatrixPhase(3)
!~       !
!~       mu              = half * one / (one + nu) * E
!~       !
!~       K = bulkModulus(nHFEDpar,parHFEDMatrixPhase)
!~       !
!~       !
!~       eps_D = deviator(eps)
!~       Ident = Identity(3)
!~       !
!~       d_HFEDposStarSplit_d_eps_e_d_eps_e = K * (d_MacAulay(trace(eps)) - gammastar * d_MacAulay_neg(trace(eps))) * DYADETwoTwo(Ident,Ident) + two * mu * Projection_dev(1)

!~     END FUNCTION d_HFEDposStarSplit_d_eps_e_d_eps_e

!~ !------------------------------------------------------------------------------------

!~     REAL(kind=AbqRK) FUNCTION d_HFEDnegStarSplit_d_eps_e_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
!~     ! Psi pos Laura de Lorenzis et al (2024)

!~       USE ABQINTERFACE_PF
!~       USE FLOATNUMBERS
!~       USE MathModul
!~       USE TensorModule

!~       IMPLICIT NONE
!~       INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
!~       REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
!~       REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
!~       REAL(kind=AbqRK) :: eps_e_i(3,3), eps_D(3,3), Ident(3,3), ETensor(3,3,3,3)
!~       REAL(kind=AbqRK) :: E, nu, K, mu, gammastar
!~       DIMENSION d_HFEDnegStarSplit_d_eps_e_d_eps_e(3,3,3,3)
!~       !
!~       d_HFEDnegStarSplit_d_eps_e_d_eps_e = zero
!~       !
!~       E               = parHFEDMatrixPhase(1)
!~       nu              = parHFEDMatrixPhase(2)
!~       gammastar 	  = parHFEDMatrixPhase(3)
!~       !
!~       mu              = half * one / (one + nu) * E
!~       !
!~       K = bulkModulus(nHFEDpar,parHFEDMatrixPhase)
!~       !
!~       !
!~       eps_D = deviator(eps)
!~       Ident = Identity(3)
!~       !
!~       d_HFEDnegStarSplit_d_eps_e_d_eps_e = K * (one+gammastar) * d_MacAulay_neg(trace(eps)) * DYADETwoTwo(Ident,Ident)
!~       !
!~     END FUNCTION d_HFEDnegStarSplit_d_eps_e_d_eps_e
    
!~ !------------------------------------------------------------------------------------
!~ !
!~ !
!~ !
!~ ! Miehe et al (sceptral decomposition)
!~ ! 2010
!~ ! DOI: 10.1002/nme.2861
!~ !
!~ ! Notes: 
!~ !
!~ !	- General Formulation from Miehe
!~ !	- Formulation of the stiffness Matrix from 10.1016/j.apples.2021.100050
!~ !
!~ !------------------------------------------------------------------------------------

!~     REAL(kind=AbqRK) FUNCTION HFEDposMieheSplit(eps,nHFEDpar,parHFEDMatrixPhase)
!~     ! Psi pos Miehe et al (2010)

!~       USE ABQINTERFACE_PF
!~       USE FLOATNUMBERS
!~       USE MathModul
!~       USE TensorModule

!~       IMPLICIT NONE
!~       INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
!~       REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
!~       REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
!~       REAL(kind=AbqRK) :: eps_sceptral(3,3), eps_sceptral_pos(3,3) 
!~       REAL(kind=AbqRK) :: E, nu, lambda, mu, Dim2DCheck, TOL 
!~       INTEGER(kind=AbqIK) :: i, n
!~       !
!~       HFEDposMieheSplit = zero
!~       eps_sceptral_pos = zero
!~       !
!~       TOL = 1e-10
!~       !
!~       E               = parHFEDMatrixPhase(1)
!~       nu              = parHFEDMatrixPhase(2)
!~       !
!~       mu              = half * one / (one + nu) * E
!~       lambda          = nu / (one - two * nu) * one / (one + nu) * E
!~       !
!~       Dim2DCheck = eps(1,3)**two + eps(2,3)**two + eps(3,1)**two + eps(3,2)**two + eps(3,3)**two
!~       IF (Dim2DCheck < TOL) THEN
!~ 		n = 2
!~ 	  ELSE
!~ 	    n = 3
!~ 	  END IF
!~       !
!~       eps_sceptral(1:n,1:n) 	  = tensor_sceptral(eps(1:n,1:n), n)
!~       !
!~       DO i = 1, 3
!~ 		eps_sceptral_pos(i,i) = MacAulay(eps_sceptral(i,i))
!~ 	  END DO
      
!~       !
!~       ! Debugging
!~       !
!~       !WRITE(6,*) '------------------------------------------------- '
!~       !WRITE(6,*) '--------------HFED_POS_MIEHE--Davor-------------- '
!~       !WRITE(6,*) '------------------------------------------------- '
!~ 	  !WRITE(6,*) 'eps: ', eps
!~ 	  !WRITE(6,*) 'eps_sceptral: ', eps_sceptral
!~ 	  !WRITE(6,*) '------------------------------------------------- '
      
      
      
!~       ! Unterschied ob MacAulay von trace() oder trace von MacAulay der Hauptdehnungen
      
!~       !HFEDposMieheSplit = half * lambda * MacAulay(trace(eps_sceptral_pos)) ** two + mu * DOUBLECONTRACTIONTwoTwo(eps_sceptral_pos,eps_sceptral_pos) ! MP
!~       HFEDposMieheSplit = half * lambda * MacAulay(trace(eps_sceptral)) ** two + mu * DOUBLECONTRACTIONTwoTwo(eps_sceptral_pos,eps_sceptral_pos) ! Miehe original
      
!~       !WRITE(6,*) '------------------------------------------------- '
!~ 	  !WRITE(6,*) 'lambda: ', lambda
!~ 	  !WRITE(6,*) 'trace(eps_sceptral): ', trace(eps_sceptral)
!~ 	  !WRITE(6,*) 'FirstTerm: ', half * lambda * MacAulay(trace(eps_sceptral)) ** two
!~ 	  !WRITE(6,*) 'mu: ', mu
!~ 	  !WRITE(6,*) 'Doublecontract: ', mu * DOUBLECONTRACTIONTwoTwo(eps_sceptral_pos,eps_sceptral_pos)
!~ 	  !WRITE(6,*) '------------------------------------------------- '
!~ 	  !WRITE(6,*) 'Gesamt: : ', HFEDposMieheSplit
!~ 	  !WRITE(6,*) '------------------------------------------------- '
      
!~     END FUNCTION HFEDposMieheSplit

!~ !------------------------------------------------------------------------------------

!~     REAL(kind=AbqRK) FUNCTION HFEDnegMieheSplit(eps,nHFEDpar,parHFEDMatrixPhase)
!~     ! Psi neg Miehe et al (2010)

!~       USE ABQINTERFACE_PF
!~       USE FLOATNUMBERS
!~       USE MathModul
!~       USE TensorModule

!~       IMPLICIT NONE
!~       INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
!~       REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
!~       REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
!~       REAL(kind=AbqRK) :: eps_sceptral(3,3), eps_sceptral_neg(3,3) 
!~       REAL(kind=AbqRK) :: E, nu, lambda, mu, Dim2DCheck, TOL
!~       INTEGER(kind=AbqIK) :: i, n
!~       !
!~       HFEDnegMieheSplit = zero
!~       eps_sceptral_neg = zero
!~       !
!~       TOL = 1e-10
!~       !
!~       E               = parHFEDMatrixPhase(1)
!~       nu              = parHFEDMatrixPhase(2)
!~       !
!~       mu              = half * one / (one + nu) * E
!~       lambda          = nu / (one - two * nu) * one / (one + nu) * E
!~       !
!~       Dim2DCheck = eps(1,3)**two + eps(2,3)**two + eps(3,1)**two + eps(3,2)**two + eps(3,3)**two
!~       IF (Dim2DCheck < TOL) THEN
!~ 		n = 2
!~ 	  ELSE
!~ 	    n = 3
!~ 	  END IF
!~       !
!~       eps_sceptral(1:n,1:n) 	  = tensor_sceptral(eps(1:n,1:n), n)
!~       !
!~       DO i = 1, 3
!~ 		eps_sceptral_neg(i,i) = MacAulay_neg(eps_sceptral(i,i))
!~ 	  END DO
      
!~       ! Unterschied ob MacAulay von trace() oder trace von MacAulay der Hauptdehnungen
      
!~       !HFEDnegMieheSplit = half * lambda * MacAulay(trace(eps_sceptral_neg)) ** two + mu * DOUBLECONTRACTIONTwoTwo(eps_sceptral_neg,eps_sceptral_neg) ! MP
!~       HFEDnegMieheSplit = half * lambda * MacAulay_neg(trace(eps_sceptral)) ** two + mu * DOUBLECONTRACTIONTwoTwo(eps_sceptral_neg,eps_sceptral_neg) ! Miehe original
      
!~       !WRITE(6,*) '------------------------------------------------- '
!~ 	  !WRITE(6,*) 'lambda: ', lambda
!~ 	  !WRITE(6,*) 'trace(eps_sceptral): ', trace(eps_sceptral)
!~ 	  !WRITE(6,*) 'FirstTerm: ', half * lambda * MacAulay_neg(trace(eps_sceptral)) ** two
!~ 	  !WRITE(6,*) 'mu: ', mu
!~ 	  !WRITE(6,*) '2 Term: ', two + mu * DOUBLECONTRACTIONTwoTwo(eps_sceptral_neg,eps_sceptral_neg)
!~ 	  !WRITE(6,*) '------------------------------------------------- '
!~ 	  !WRITE(6,*) 'Gesamt: : ', HFEDnegMieheSplit
!~ 	  !WRITE(6,*) '------------------------------------------------- '
      
!~     END FUNCTION HFEDnegMieheSplit

!~ !------------------------------------------------------------------------------------
!~ !
!~ !
!~ ! 1. Ableitung
!~ !
!~ !
!~ !------------------------------------------------------------------------------------

!~     REAL(kind=AbqRK) FUNCTION d_HFEDposMieheSplit_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
!~     ! d_Psi_d_eps pos Miehe et al (2010)
    
!~     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!~ 	!
!~ 	!	läuft			
!~ 	!
!~ 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!~       USE ABQINTERFACE_PF
!~       USE FLOATNUMBERS
!~       USE MathModul
!~       USE TensorModule

!~       IMPLICIT NONE
!~       INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
!~       REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
!~       REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
!~       REAL(kind=AbqRK) :: eps_sceptral(3,3), eps_sceptral_pos(3,3), P(3,3), d_HFEDposMieheSplit_d_eps_e_sceptral(3,3)
!~       REAL(kind=AbqRK) :: W(3)
!~       REAL(kind=AbqRK) :: E, nu, lambda, mu, Dim2DCheck, TOL
!~       INTEGER(kind=AbqIK) :: i, n
!~       DIMENSION d_HFEDposMieheSplit_d_eps_e(3,3)
!~       !
!~       !
!~       TOL = 1e-10
!~       !
!~       eps_sceptral_pos = zero
!~       eps_sceptral = zero
!~       W = zero
!~       P = zero
!~       d_HFEDposMieheSplit_d_eps_e_sceptral = zero
!~       !
!~       E               = parHFEDMatrixPhase(1)
!~       nu              = parHFEDMatrixPhase(2)
!~       !
!~       mu              = half * one / (one + nu) * E
!~       lambda          = nu / (one - two * nu) * one / (one + nu) * E
!~       !
!~       Dim2DCheck = eps(1,3)**two + eps(2,3)**two + eps(3,1)**two + eps(3,2)**two + eps(3,3)**two
!~       IF (Dim2DCheck < TOL) THEN
!~ 		n = 2
!~ 	  ELSE
!~ 	    n = 3
!~ 	  END IF
!~       !
!~ 	  !WRITE(6,*) '------------------------------------------------- '
!~       !WRITE(6,*) '-------------VorEps_sceptral--------------------- '
!~       !WRITE(6,*) '------------------------------------------------- '
!~ 	  !WRITE(6,*) 'eps: ', eps
!~ 	  !WRITE(6,*) 'eps_sceptral: ', eps_sceptral
!~ 	  !WRITE(6,*) 'Dimension: ', n
!~ 	  !WRITE(6,*) '------------------------------------------------- '
!~       !
!~       eps_sceptral(1:n,1:n) 	  = tensor_sceptral(eps(1:n,1:n), n)
!~       !
!~       !WRITE(6,*) '------------------------------------------------- '
!~       !WRITE(6,*) '-------------NachEps_sceptral-------------------- '
!~       !WRITE(6,*) '------------------------------------------------- '
!~ 	  !WRITE(6,*) 'eps: ', eps
!~ 	  !WRITE(6,*) 'eps_sceptral: ', eps_sceptral
!~ 	  !WRITE(6,*) '------------------------------------------------- '
!~       !
!~       DO i = 1, n
!~ 		d_HFEDposMieheSplit_d_eps_e_sceptral(i,i) = lambda * MacAulay(trace(eps_sceptral)) + two * mu * MacAulay(eps_sceptral(i,i))
!~ 		!
!~ 	    !WRITE(6,*) '------------------------------------------------- '
!~         !WRITE(6,*) '----d_HFEDposMieheSplit_d_eps_e_sceptral(i,i)---- '
!~         !WRITE(6,*) '------------------------------------------------- '
!~ 	    !WRITE(6,*) 'eps: ', eps_sceptral
!~   	    !WRITE(6,*) 'lambda: ', lambda
!~ 	    !WRITE(6,*) 'tr(eps_sceptral): ', trace(eps_sceptral)
!~ 	    !WRITE(6,*) 'd_HFEDposMieheSplit_d_eps_e_sceptral, erster Term: ', lambda * MacAulay(trace(eps_sceptral))
!~ 	    !WRITE(6,*) 'd_HFEDposMieheSplit_d_eps_e_sceptral, zweiter Term: ', two * mu * MacAulay(eps_sceptral(i,i))
!~ 	    !WRITE(6,*) 'd_HFEDposMieheSplit_d_eps_e_sceptral: ', d_HFEDposMieheSplit_d_eps_e_sceptral
!~ 	    !WRITE(6,*) '------------------------------------------------- '
!~ 	    !
!~ 	  END DO
!~ 	  !
!~ 	  ! Weil n für den trace-Term in 2D nur bis zum vorletzten Wert geht
!~ 	  !
!~ 	  d_HFEDposMieheSplit_d_eps_e_sceptral(3,3) = lambda * MacAulay(trace(eps_sceptral))
!~ 	  !
!~ 	  !
!~ 	  !WRITE(6,*) '------------------------------------------------- '
!~       !WRITE(6,*) '---Vor Retransformation d_HFED_POS_MIEHE_d_eps_e- '
!~       !WRITE(6,*) '------------------------------------------------- '
!~ 	  !WRITE(6,*) 'eps: ', eps_sceptral
!~ 	  !WRITE(6,*) 'lambda: ', lambda
!~ 	  !WRITE(6,*) 'tr(eps_sceptral): ', trace(eps_sceptral)
!~ 	  !WRITE(6,*) 'd_HFEDposMieheSplit_d_eps_e_sceptral: ', lambda * MacAulay(trace(eps_sceptral))
!~ 	  !WRITE(6,*) '------------------------------------------------- '
!~ 	  !
!~ 	  ! Get Eigenvectors for Retransformation
!~       IF (SQRT(SUM(eps**two)) .GT. TOL) THEN
!~ 		  CALL principal_system(eps(1:n, 1:n), P(1:n, 1:n), W(1:n), n)
!~ 		  !
!~ 		  CALL retransformation(d_HFEDposMieheSplit_d_eps_e_sceptral(1:n, 1:n), P(1:n,1:n), n)
!~ 		  !WRITE(6,*) '------------------------------------------------- '
!~ 		  !WRITE(6,*) '--Nach der Sceptral d_HFED_POS_MIEHE_d_eps_e----- '
!~ 		  !WRITE(6,*) '------------------------------------------------- '
!~ 		  !WRITE(6,*) 'd_HFEDposMieheSplit_d_eps_e_sceptral: ', d_HFEDposMieheSplit_d_eps_e_sceptral
!~ 		  !WRITE(6,*) 'P: ', P
!~ 		  !WRITE(6,*) 'W: ', W
!~ 		  !WRITE(6,*) '------------------------------------------------- '
!~ 	  END IF
!~       !
	  
!~ 	  d_HFEDposMieheSplit_d_eps_e = d_HFEDposMieheSplit_d_eps_e_sceptral
!~ 	  !
!~ 	  !
!~ 	  !WRITE(6,*) '------------------------------------------------- '
!~       !WRITE(6,*) '----Nach d_HFED_POS_MIEHE_d_eps_e---------------- '
!~       !WRITE(6,*) '------------------------------------------------- '
!~ 	  !WRITE(6,*) 'eps: ', eps
!~ 	  !WRITE(6,*) 'eps_sceptral: ', eps_sceptral
!~ 	  !WRITE(6,*) '------------------------------------------------- '
	  
      
!~     END FUNCTION d_HFEDposMieheSplit_d_eps_e

!~ !------------------------------------------------------------------------------------

!~     REAL(kind=AbqRK) FUNCTION d_HFEDnegMieheSplit_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
!~     ! d_Psi_d_eps neg Miehe et al (2010)
    
!~     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!~ 	!
!~ 	!	läuft			
!~ 	!
!~ 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!~       USE ABQINTERFACE_PF
!~       USE FLOATNUMBERS
!~       USE MathModul
!~       USE TensorModule

!~       IMPLICIT NONE
!~       INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
!~       REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
!~       REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
!~       REAL(kind=AbqRK) :: eps_sceptral(3,3), eps_sceptral_neg(3,3), P(3,3), d_HFEDnegMieheSplit_d_eps_e_sceptral(3,3)
!~       REAL(kind=AbqRK) :: W(3)
!~       REAL(kind=AbqRK) :: E, nu, lambda, mu, Dim2DCheck, TOL
!~       INTEGER(kind=AbqIK) :: i, n
!~       DIMENSION d_HFEDnegMieheSplit_d_eps_e(3,3)
!~       !
!~       !
!~       TOL = 1e-10
!~       !
!~       eps_sceptral = zero
!~       eps_sceptral_neg = zero
!~       P = zero
!~       W = zero
!~       d_HFEDnegMieheSplit_d_eps_e_sceptral = zero
!~       !
!~       E               = parHFEDMatrixPhase(1)
!~       nu              = parHFEDMatrixPhase(2)
!~       !
!~       mu              = half * one / (one + nu) * E
!~       lambda          = nu / (one - two * nu) * one / (one + nu) * E
!~       !
!~       Dim2DCheck = eps(1,3)**two + eps(2,3)**two + eps(3,1)**two + eps(3,2)**two + eps(3,3)**two
!~       IF (Dim2DCheck < TOL) THEN
!~ 		n = 2
!~ 	  ELSE
!~ 	    n = 3
!~ 	  END IF
!~       !
!~       eps_sceptral(1:n,1:n) 	  = tensor_sceptral(eps(1:n,1:n), n)
!~       !
!~       DO i = 1, n
!~ 		d_HFEDnegMieheSplit_d_eps_e_sceptral(i,i) = lambda * MacAulay_neg(trace(eps_sceptral)) + two * mu * MacAulay_neg(eps_sceptral(i,i))
!~ 		!
!~ 	    !WRITE(6,*) '------------------------------------------------- '
!~         !WRITE(6,*) '----d_HFEDnegMieheSplit_d_eps_e_sceptral(i,i)---- '
!~         !WRITE(6,*) '------------------------------------------------- '
!~ 	    !WRITE(6,*) 'eps: ', eps_sceptral
!~   	    !WRITE(6,*) 'lambda: ', lambda
!~ 	    !WRITE(6,*) 'tr(eps_sceptral): ', trace(eps_sceptral)
!~ 	    !WRITE(6,*) 'd_HFEDnegMieheSplit_d_eps_e_sceptral, erster Term: ', lambda * MacAulay_neg(trace(eps_sceptral))
!~ 	    !WRITE(6,*) 'd_HFEDnegMieheSplit_d_eps_e_sceptral, zweiter Term: ', two * mu * MacAulay_neg(eps_sceptral(i,i))
!~ 	    !WRITE(6,*) 'd_HFEDnegMieheSplit_d_eps_e_sceptral: ', d_HFEDnegMieheSplit_d_eps_e_sceptral
!~ 	    !WRITE(6,*) '------------------------------------------------- '
!~ 	    !
!~ 	  END DO
!~ 	  !
!~ 	  !
!~ 	  ! Weil der trace Term für n = 2 in 2D nur bis zum vorletzten Wert der Hauptdiagonale geht
!~ 	  !
!~ 	  d_HFEDnegMieheSplit_d_eps_e_sceptral(3,3) = lambda * MacAulay_neg(trace(eps_sceptral))
!~ 	  !
!~       !
!~       !------------------------
!~ 	  ! Get Eigenvectors for Retransformation
!~ 	  	!WRITE(6,*) '---------------------------------------------'
!~ 		!WRITE(6,*) 'NACH DER SCHLEIFE', d_HFEDnegMieheSplit_d_eps_e_sceptral
!~ 		!WRITE(6,*) '---------------------------------------------'
!~       IF (SQRT(SUM(eps**two)) .GT. TOL) THEN
!~ 		  CALL principal_system(eps(1:n, 1:n), P(1:n, 1:n), W(1:n), n)
!~ 		  !
!~ 		  !DO i = 1, n
!~ 		  !	P(i,:) = P(i,:) * W(i)
!~ 		  !END DO
!~ 		  !
!~ 		  CALL retransformation(d_HFEDnegMieheSplit_d_eps_e_sceptral(1:n, 1:n), P(1:n,1:n), n)
!~ 	  END IF
!~       !
!~ 	  !
      
!~       d_HFEDnegMieheSplit_d_eps_e = d_HFEDnegMieheSplit_d_eps_e_sceptral
!~       !WRITE(6,*) '------------------------------------------------- '
!~       !WRITE(6,*) '---------d_HFED_NEG_MIEHE_d_eps_e---------------- '
!~       !WRITE(6,*) '------------------------------------------------- '
!~ 	  !WRITE(6,*) 'eps: ', eps
!~ 	  !WRITE(6,*) 'eps_sceptral: ', eps_sceptral
!~ 	  !WRITE(6,*) 'd_HFEDnegMieheSplit_d_eps_e: ', d_HFEDnegMieheSplit_d_eps_e
!~ 	  !WRITE(6,*) 'trace(eps_sceptral): ', trace(eps_sceptral)
!~ 	  !WRITE(6,*) 'MacAulay_neg(trace(eps_sceptral)): ', MacAulay_neg(trace(eps_sceptral))
!~ 	  !WRITE(6,*) 'MacAulay_neg(eps_sceptral(i,i)): ', MacAulay_neg(eps_sceptral(i,i))
!~ 	  !WRITE(6,*) 'Erster Term: ', lambda * MacAulay_neg(trace(eps_sceptral))
!~ 	  !WRITE(6,*) 'Zweiter Term: ', two * mu * MacAulay_neg(eps_sceptral(i,i))
!~ 	  !WRITE(6,*) '------------------------------------------------- '
!~ 	  !
!~     END FUNCTION d_HFEDnegMieheSplit_d_eps_e

!~ !------------------------------------------------------------------------------------
!~ !
!~ !
!~ !
!~ !		2. Ableitung
!~ !
!~ !
!~ !
!~ !------------------------------------------------------------------------------------

!~     REAL(kind=AbqRK) FUNCTION d_HFEDposMieheSplit_d_eps_e_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
!~     ! Psi pos Miehe et al (2010)
    
!~     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!~ 	!
!~ 	!	Not working so far			!
!~ 	!
!~ 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!~       USE ABQINTERFACE_PF
!~       USE FLOATNUMBERS
!~       USE MathModul
!~       USE TensorModule

!~       IMPLICIT NONE
!~       INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
!~       REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
!~       REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
!~       INTEGER(kind=AbqIK) :: n
!~       REAL(kind=AbqRK) :: J(3,3,3,3), P_pos(3,3,3,3), P_neg(3,3,3,3)
!~       REAL(kind=AbqRK) :: EV(3,3), eps_sceptral(3,3)
!~       REAL(kind=AbqRK) :: EW(3)
!~       REAL(kind=AbqRK) :: E, nu, mu, Dim2DCheck, lambda, TOL
!~       DIMENSION d_HFEDposMieheSplit_d_eps_e_d_eps_e(3,3,3,3)
!~       !
!~       d_HFEDposMieheSplit_d_eps_e_d_eps_e = zero
!~       EV = zero
!~       EW = zero
!~       eps_sceptral = zero
!~       J = zero
!~       P_pos = zero
!~       P_neg = zero
!~ 	  !
!~ 	  TOL = 1e-10
!~       !
!~       E               = parHFEDMatrixPhase(1)
!~       nu              = parHFEDMatrixPhase(2)
!~       !
!~       mu              = half * one / (one + nu) * E
!~       lambda          = nu / (one - two * nu) * one / (one + nu) * E
!~       !
!~       Dim2DCheck = eps(1,3)**two + eps(2,3)**two + eps(3,1)**two + eps(3,2)**two + eps(3,3)**two
!~       IF (Dim2DCheck < TOL) THEN
!~ 		n = 2
!~ 	  ELSE
!~ 	    n = 3
!~ 	  END IF
!~       !
!~       eps_sceptral(1:n, 1:n) 	  = tensor_sceptral(eps(1:n,1:n), n)
!~       !
!~       CALL principal_system(eps(1:n, 1:n), EV(1:n, 1:n), EW(1:n), n)
!~       !
!~       !
!~       CALL Miehe_Stiffness_Additional_Tensors(EV, EW, J, P_pos, P_neg)
!~       !
!~       !
!~       d_HFEDposMieheSplit_d_eps_e_d_eps_e = lambda * J + two*mu*P_pos

!~     END FUNCTION d_HFEDposMieheSplit_d_eps_e_d_eps_e

!~ !------------------------------------------------------------------------------------

!~     REAL(kind=AbqRK) FUNCTION d_HFEDnegMieheSplit_d_eps_e_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
!~     ! Psi neg Miehe et al (2010)
	
!~ 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!~ 	!
!~ 	!	Not working so far			!
!~ 	!
!~ 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	
!~       USE ABQINTERFACE_PF
!~       USE FLOATNUMBERS
!~       USE MathModul
!~       USE TensorModule

!~       IMPLICIT NONE
!~       INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
!~       REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
!~       REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
!~       INTEGER(kind=AbqIK) :: n
!~       REAL(kind=AbqRK) :: J(3,3,3,3), P_pos(3,3,3,3), P_neg(3,3,3,3)
!~       REAL(kind=AbqRK) :: EV(3,3), eps_sceptral(3,3)
!~       REAL(kind=AbqRK) :: EW(3)
!~       REAL(kind=AbqRK) :: E, nu, mu, Dim2DCheck, lambda, TOL
!~       DIMENSION d_HFEDnegMieheSplit_d_eps_e_d_eps_e(3,3,3,3)
!~       !
!~       d_HFEDnegMieheSplit_d_eps_e_d_eps_e = zero
!~       EV = zero
!~       EW = zero
!~       eps_sceptral = zero
!~       J = zero
!~       P_pos = zero
!~       P_neg = zero
!~       !
!~       TOL = 1e-10
!~       !
!~       E               = parHFEDMatrixPhase(1)
!~       nu              = parHFEDMatrixPhase(2)
!~       !
!~       mu              = half * one / (one + nu) * E
!~       lambda          = nu / (one - two * nu) * one / (one + nu) * E
!~       !
!~       Dim2DCheck = eps(1,3)**two + eps(2,3)**two + eps(3,1)**two + eps(3,2)**two + eps(3,3)**two
!~       IF (Dim2DCheck < TOL) THEN
!~ 		n = 2
!~ 	  ELSE
!~ 	    n = 3
!~ 	  END IF
!~       !
!~       eps_sceptral(1:n,1:n) 	  = tensor_sceptral(eps(1:n,1:n), n)
!~       !
!~       CALL principal_system(eps(1:n, 1:n), EV(1:n, 1:n), EW(1:n), n)
!~       !
!~       !
!~       CALL Miehe_Stiffness_Additional_Tensors(EV, EW, J, P_pos, P_neg)
!~       !
!~       !
!~       d_HFEDnegMieheSplit_d_eps_e_d_eps_e = -lambda * J + two*mu*P_neg

!~     END FUNCTION d_HFEDnegMieheSplit_d_eps_e_d_eps_e


!~ !------------------------------------------------------------------------------------
!~ !
!~ !
!~ !
!~ !		1. Ableitung Itskov/Abendroth Variante
!~ !
!~ !
!~ !
!~ !------------------------------------------------------------------------------------

!~     REAL(kind=AbqRK) FUNCTION d_HFEDposMiehe_Itskov_Split_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
!~     ! Sigma pos Miehe et al (2010)
    
!~     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!~ 	!
!~ 	!	Not working so far			!
!~ 	!
!~ 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!~       USE ABQINTERFACE_PF
!~       USE FLOATNUMBERS
!~       USE MathModul
!~       USE TensorModule
!~       USE CZMath

!~       IMPLICIT NONE
!~       INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
!~       REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
!~       REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
!~       INTEGER(kind=AbqIK) :: n, i
!~       REAL(kind=AbqRK) :: T1(3,3), T2(3,3), T_sum(3,3,3,3)
!~       REAL(kind=AbqRK) :: eps_pos(3,3)
!~       REAL(kind=AbqRK) :: EV(3,3), EW(3), Ident(3,3)
!~       REAL(kind=AbqRK) :: d_EW_d_eps(3,3,3) 
!~       REAL(kind=AbqRK) :: P(3,3,3), dPdE(3,3,3,3,3) ! aus Subroutine d_EV_dE
!~       REAL(kind=AbqRK) :: E, nu, mu, lambda, TOL
!~       DIMENSION d_HFEDposMiehe_Itskov_Split_d_eps_e(3,3)
!~       !
!~       d_HFEDposMiehe_Itskov_Split_d_eps_e = zero
!~       eps_pos = zero
!~       !
!~       TOL = 1e-10
!~       !
!~       E               = parHFEDMatrixPhase(1)
!~       nu              = parHFEDMatrixPhase(2)
!~       !
!~       mu              = half * one / (one + nu) * E
!~       lambda          = nu / (one - two * nu) * one / (one + nu) * E
!~       !
!~       n = 3
!~       !
!~       Ident = Identity(3)
!~       !
!~       CALL principal_system(eps, EV, EW, n)
!~       !
!~       !
!~       ! positive Eigenwerte sceptral Epsilon^+
!~       !
!~       DO i = 1,3
!~ 		!WRITE(6,*) '------------------------------------'
!~ 		!WRITE(6,*) '------------------------------------'
!~ 		!WRITE(6,*) 'i: ', i
!~ 		!WRITE(6,*) 'EW: ', EW
!~ 		!WRITE(6,*) 'eps_pos: ', eps_pos
!~ 		!WRITE(6,*) '------------------------------------'
!~ 		IF (MacAulay(EW(i)) .GT. TOL) THEN
!~ 			eps_pos(i,i) = EW(i)
!~ 		END IF
!~       END DO
!~       !
!~       !
!~       !
!~       !WRITE(6,*) '------------------------------------'   
!~       !WRITE(6,*) '-----__Vor d_EV_d_eps________--------'  
!~       !WRITE(6,*) '------------------------------------'  
!~       !WRITE(6,*) 'eps: ', eps
!~ 	  !WRITE(6,*) 'eps_pos: ', eps_pos
!~       !WRITE(6,*) '------------------------------------' 
!~       !     
!~       CALL d_EV_d_eps(eps, n, P, dPdE)
!~       !
!~       !WRITE(6,*) '------------------------------------'   
!~       !WRITE(6,*) '-----__Vor d_EW_d_eps_______--------'  
!~       !WRITE(6,*) '------------------------------------'  
!~       !WRITE(6,*) 'eps: ', eps
!~       !WRITE(6,*) 'EV: ', P
!~       !WRITE(6,*) 'EW: ', dPdE
!~       !WRITE(6,*) '------------------------------------' 
!~       !
!~       !
!~       CALL d_Eigenwerte_d_eps(eps, d_EW_d_eps)
!~       !
!~       !
!~       !WRITE(6,*) '------------------------------------'   
!~       !WRITE(6,*) '-----__Nach_d_EV_d_eps______--------'  
!~       !WRITE(6,*) '------------------------------------'  
!~       !WRITE(6,*) 'eps: ', eps
!~       !WRITE(6,*) 'EV: ', P
!~       !WRITE(6,*) 'EW: ', dPdE
!~       !WRITE(6,*) '------------------------------------' 
!~       T_sum = zero
      
!~       DO i = 1,3
!~ 		T_sum = T_sum + d_MacAulay(EW(i))*DYADETwoTwo((d_EW_d_eps(:,:,i)),P(:,:,i)) + MacAulay(EW(i))*dPdE(:,:,:,:,i)
!~ 	  END DO
      
!~       !
!~       ! Term 1
!~       ! Term 2
!~       !
!~       T1 = lambda * MacAulay(trace(eps)) * Ident
!~       T2 = two * mu * DOUBLECONTRACTIONFourTwo(T_sum,eps_pos)
!~       !
!~       !WRITE(6,*) '------------------------------------'   
!~       !WRITE(6,*) '----- TERMS_Stress_pos -----'  
!~       !WRITE(6,*) '------------------------------------'  
!~       !WRITE(6,*) 'eps: ', eps
!~       !WRITE(6,*) 'T_sum: ', T_sum
!~       !WRITE(6,*) 'eps_pos: ', eps_pos
!~       !WRITE(6,*) 'lambda * MacAulay(trace(eps)) * Ident: ', T1
!~       !WRITE(6,*) 'two * mu * DOUBLECONTRACTIONFourTwo(T_sum,eps_pos): ', T2
!~       !WRITE(6,*) '------------------------------------' 
!~       !
!~       d_HFEDposMiehe_Itskov_Split_d_eps_e = T1 + T2

!~     END FUNCTION d_HFEDposMiehe_Itskov_Split_d_eps_e

!~ !------------------------------------------------------------------------------------

!~     REAL(kind=AbqRK) FUNCTION d_HFEDnegMiehe_Itskov_Split_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
!~     ! Psi neg Miehe et al (2010)
    
!~     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!~ 	!
!~ 	!	Not working so far			!
!~ 	!
!~ 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!~       USE ABQINTERFACE_PF
!~       USE FLOATNUMBERS
!~       USE MathModul
!~       USE TensorModule
!~       USE CZMath

!~       IMPLICIT NONE
!~       INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
!~       REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
!~       REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
!~       INTEGER(kind=AbqIK) :: n, i
!~       REAL(kind=AbqRK) :: T1(3,3), T2(3,3), T_sum(3,3,3,3)
!~       REAL(kind=AbqRK) :: eps_neg(3,3)
!~       REAL(kind=AbqRK) :: EV(3,3), EW(3), Ident(3,3)
!~       REAL(kind=AbqRK) :: P(3,3,3), dPdE(3,3,3,3,3), d_EW_d_eps(3,3,3) ! aus Subroutine d_EV_dE
!~       REAL(kind=AbqRK) :: E, nu, mu, lambda, TOL
!~       DIMENSION d_HFEDnegMiehe_Itskov_Split_d_eps_e(3,3)
!~       !
!~       d_HFEDnegMiehe_Itskov_Split_d_eps_e = zero
!~       !
!~       TOL = 1e-10
!~       !
!~       E               = parHFEDMatrixPhase(1)
!~       nu              = parHFEDMatrixPhase(2)
!~       !
!~       mu              = half * one / (one + nu) * E
!~       lambda          = nu / (one - two * nu) * one / (one + nu) * E
!~       !
!~       n = 3
!~       !
!~       Ident = Identity(3)
!~       !
!~       CALL principal_system(eps, EV, EW, n)
!~       !
!~       !
!~       ! negitive Eigenwerte sceptral Epsilon^+
!~       !
!~       DO i = 1,3
!~ 		IF (MacAulay(EW(i)) .GT. TOL) THEN
!~ 			eps_neg(i,i) = EW(i)
!~ 		END IF
!~       END DO
!~       !     
      
!~       !WRITE(6,*) '------------------------------------'   
!~       !WRITE(6,*) '-----___Vor_Tangente________--------'  
!~       !WRITE(6,*) '------------------------------------'  
!~       !WRITE(6,*) 'eps: ', eps
!~       !WRITE(6,*) 'EV: ', EV
!~       !WRITE(6,*) 'EW: ', EW
!~       !WRITE(6,*) '------------------------------------' 
!~       !
!~       CALL d_EV_d_eps(eps, n, P, dPdE)
!~       !
!~       CALL d_Eigenwerte_d_eps(eps, d_EW_d_eps)
!~       !
!~       !WRITE(6,*) '------------------------------------'   
!~       !WRITE(6,*) '-----__Nach_Tangente________--------'  
!~       !WRITE(6,*) '------------------------------------'  
!~       !WRITE(6,*) 'eps: ', eps
!~       !WRITE(6,*) 'EV: ', P
!~       !WRITE(6,*) 'EW: ', dPdE
!~       !WRITE(6,*) '------------------------------------' 
!~       T_sum = zero
      
!~       DO i = 1,3
!~ 		T_sum = T_sum + (d_MacAulay_neg(EW(i))*DYADETwoTwo((d_EW_d_eps(:,:,i)),P(:,:,i)) + MacAulay_neg(EW(i))*dPdE(:,:,:,:,i))
!~ 	  END DO
      
!~       !
!~       ! Term 1
!~       ! Term 2
!~       !
!~       T1 = lambda * MacAulay(trace(eps)) * Ident
!~       T2 = two * mu * DOUBLECONTRACTIONFourTwo(T_sum,eps_neg)
!~       !
!~       d_HFEDnegMiehe_Itskov_Split_d_eps_e = T1 + T2

!~     END FUNCTION d_HFEDnegMiehe_Itskov_Split_d_eps_e

!~ !------------------------------------------------------------------------------------
!~ !
!~ !
!~ !
!~ !		2. Ableitung Itskov/Abendroth Variante
!~ !
!~ !
!~ !
!~ !------------------------------------------------------------------------------------

!~     REAL(kind=AbqRK) FUNCTION d_HFEDposMiehe_Itskov_Split_d_eps_e_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
!~     ! Psi pos Miehe et al (2010)
    
!~     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!~ 	!
!~ 	!	Not working so far			!
!~ 	!
!~ 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!~       USE ABQINTERFACE_PF
!~       USE FLOATNUMBERS
!~       USE MathModul
!~       USE TensorModule
!~       USE CZMath

!~       IMPLICIT NONE
!~       INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
!~       REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
!~       REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
!~       INTEGER(kind=AbqIK) :: n, i
!~       REAL(kind=AbqRK) :: eps_pos(3,3)
!~       REAL(kind=AbqRK) :: EV(3,3), EW(3), EW_pos(3), EV_pos(3,3), Ident(3,3)
!~       REAL(kind=AbqRK) :: EVTens4(3,3,3,3), EVT1(3,3,3,3), EVT2(3,3,3,3), T1(3,3,3,3), T2(3,3,3,3)
!~       REAL(kind=AbqRK) :: P(3,3,3), dPdE(3,3,3,3,3) ! aus Subroutine d_EV_dE
!~       REAL(kind=AbqRK) :: temp(3,3)
!~       REAL(kind=AbqRK) :: E, nu, mu, lambda, TOL
!~       DIMENSION d_HFEDposMiehe_Itskov_Split_d_eps_e_d_eps_e(3,3,3,3)
!~       !
!~       temp = zero
!~       EVT1 = zero
!~       EVT2 = zero
!~       EVTens4 = zero
!~       d_HFEDposMiehe_Itskov_Split_d_eps_e_d_eps_e = zero
!~       !
!~       TOL = 1e-10
!~       !
!~       E               = parHFEDMatrixPhase(1)
!~       nu              = parHFEDMatrixPhase(2)
!~       !
!~       mu              = half * one / (one + nu) * E
!~       lambda          = nu / (one - two * nu) * one / (one + nu) * E
!~       !
!~       n = 3
!~       !
!~       Ident = Identity(3)
!~       !
!~       CALL principal_system(eps, EV, EW, n)
!~       !
!~       DO i = 1,3
!~ 		IF (MacAulay(EW(i)) .GT. TOL) THEN
!~ 			EW_pos(i) = EW(i)
!~ 			EV_pos(:,i) = EV_pos(:,i)
!~ 		END IF
!~       END DO
!~       !     
!~       !WRITE(6,*) '------------------------------------'   
!~       !WRITE(6,*) '-----___Vor_Tangente________--------'  
!~       !WRITE(6,*) '------------------------------------'  
!~       !WRITE(6,*) 'eps: ', eps
!~       !WRITE(6,*) 'EV: ', EV
!~       !WRITE(6,*) 'EW: ', EW
!~       !WRITE(6,*) '------------------------------------' 
!~       !
!~       CALL d_EV_d_eps(eps, n, P, dPdE)
!~       !
!~       !WRITE(6,*) '------------------------------------'   
!~       !WRITE(6,*) '-----__Nach_Tangente________--------'  
!~       !WRITE(6,*) '------------------------------------'  
!~       !WRITE(6,*) 'eps: ', eps
!~       !WRITE(6,*) 'EV: ', P
!~       !WRITE(6,*) 'EW: ', dPdE
!~       !WRITE(6,*) '------------------------------------' 
!~       !
!~       ! Term 1
!~       !
!~       T1 = lambda * DYADETwoTwo(Ident,Ident)
!~       !
!~       ! Term 2
!~       !
!~       DO i = 1,3
!~ 		temp = DYADE(EV_pos(:,i),EV_pos(:,i))
!~ 		EVT1 = two * mu * DYADETwoTwo(temp, temp)
!~ 		EVT2 = two * mu * EW(i) * dPdE(:,:,:,:,i)
!~ 		EVTens4 = EVTens4 + EVT1 + EVT2
!~       END DO
!~       !
!~       T2 = EVTens4
!~       !
!~       d_HFEDposMiehe_Itskov_Split_d_eps_e_d_eps_e = T1 + T2

!~     END FUNCTION d_HFEDposMiehe_Itskov_Split_d_eps_e_d_eps_e

!~ !------------------------------------------------------------------------------------

!~     REAL(kind=AbqRK) FUNCTION d_HFEDnegMiehe_Itskov_Split_d_eps_e_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
!~     ! Psi neg Miehe et al (2010)
!~ 	!
!~ 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!~ 	!
!~ 	!	Not working so far			!
!~ 	!
!~ 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!~       USE ABQINTERFACE_PF
!~       USE FLOATNUMBERS
!~       USE MathModul
!~       USE TensorModule
!~       USE CZMath

!~       IMPLICIT NONE
!~       INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
!~       REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
!~       REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
!~       REAL(kind=AbqRK) :: eps_neg(3,3)
!~       INTEGER(kind=AbqIK) :: n, i
!~       REAL(kind=AbqRK) :: EV(3,3), EW(3), EW_neg(3), EV_neg(3,3), Ident(3,3)
!~       REAL(kind=AbqRK) :: EVTens4(3,3,3,3), EVT1(3,3,3,3), EVT2(3,3,3,3), T1(3,3,3,3), T2(3,3,3,3)
!~       REAL(kind=AbqRK) :: P(3,3,3), dPdE(3,3,3,3,3) ! aus Subroutine d_EV_dE
!~       REAL(kind=AbqRK) :: temp(3,3)
!~       REAL(kind=AbqRK) :: E, nu, mu, lambda, TOL
!~       DIMENSION d_HFEDnegMiehe_Itskov_Split_d_eps_e_d_eps_e(3,3,3,3)
!~       !
!~       temp = zero
!~       EVT1 = zero
!~       EVT2 = zero
!~       EVTens4 = zero
!~       d_HFEDnegMiehe_Itskov_Split_d_eps_e_d_eps_e = zero
!~       !
!~       TOL = 1e-10
!~       !
!~       E               = parHFEDMatrixPhase(1)
!~       nu              = parHFEDMatrixPhase(2)
!~       !
!~       mu              = half * one / (one + nu) * E
!~       lambda          = nu / (one - two * nu) * one / (one + nu) * E
!~       !
!~       n = 3
!~       !
!~       Ident = Identity(3)
!~       !
!~       CALL principal_system(eps, EV, EW, n)
!~       !
!~       DO i = 1,3
!~ 		IF (MacAulay_neg(EW(i)) .GT. TOL) THEN
!~ 			EW_neg(i) = EW(i)
!~ 			EV_neg(:,i) = EV_neg(:,i)
!~ 		END IF
!~       END DO
!~       !        
!~       CALL d_EV_d_eps(eps, n, P, dPdE)
!~       !
!~       ! Term 1
!~       !
!~       T1 = lambda * DYADETwoTwo(Ident,Ident)
!~       !
!~       ! Term 2
!~       !
!~       DO i = 1,3
!~ 		temp = DYADE(EV_neg(:,i),EV_neg(:,i))
!~ 		EVT1 = two * mu * DYADETwoTwo(temp, temp)
!~ 		EVT2 = two * mu * EW(i) * dPdE(:,:,:,:,i)
!~ 		EVTens4 = EVTens4 + EVT1 + EVT2
!~       END DO
!~       !
!~       T2 = EVTens4
!~       !
!~       d_HFEDnegMiehe_Itskov_Split_d_eps_e_d_eps_e = T1 + T2

!~     END FUNCTION d_HFEDnegMiehe_Itskov_Split_d_eps_e_d_eps_e

!~ !------------------------------------------------------------------------------------
!~ !
!~ !
!~ ! 1. Ableitung
!~ !
!~ !
!~ !------------------------------------------------------------------------------------

!~     REAL(kind=AbqRK) FUNCTION d_HFEDposMieheSplit_NOTRANSFORMATION_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
!~     ! d_Psi_d_eps pos Miehe et al (2010)
    
!~     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!~ 	!
!~ 	!	läuft			
!~ 	!
!~ 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!~       USE ABQINTERFACE_PF
!~       USE FLOATNUMBERS
!~       USE MathModul
!~       USE TensorModule

!~       IMPLICIT NONE
!~       INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
!~       REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
!~       REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
!~       REAL(kind=AbqRK) :: eps_sceptral(3,3), eps_sceptral_pos(3,3), P(3,3), Proj(3,3,3), d_HFEDposMieheSplit_d_eps_e_sceptral(3,3)
!~       REAL(kind=AbqRK) :: W(3)
!~       REAL(kind=AbqRK) :: E, nu, lambda, mu, Dim2DCheck, TOL
!~       INTEGER(kind=AbqIK) :: i, n
!~       DIMENSION d_HFEDposMieheSplit_NOTRANSFORMATION_d_eps_e(3,3)
!~       !
!~       !
!~       TOL = 1e-10
!~       !
!~       eps_sceptral_pos = zero
!~       eps_sceptral = zero
!~       W = zero
!~       P = zero
!~       d_HFEDposMieheSplit_d_eps_e_sceptral = zero
!~       !
!~       E               = parHFEDMatrixPhase(1)
!~       nu              = parHFEDMatrixPhase(2)
!~       !
!~       mu              = half * one / (one + nu) * E
!~       lambda          = nu / (one - two * nu) * one / (one + nu) * E
!~       !
!~       Dim2DCheck = SQRT(eps(1,3)**two + eps(2,3)**two + eps(3,1)**two + eps(3,2)**two + eps(3,3)**two)
!~       !WRITE(6,*) '------------------------------------------------- '
!~       !WRITE(6,*) '---------------- DIMENSION CHECK ---------------- '
!~       !WRITE(6,*) '------------------------------------------------- '
!~       !WRITE(6,*) 'eps: ', eps
!~       !WRITE(6,*) 'eps 3,3: ', eps(3,3)
!~       !WRITE(6,*) 'Dim2DCheck: ', Dim2DCheck
!~       !WRITE(6,*) '------------------------------------------------- '
!~       IF (Dim2DCheck < TOL) THEN
!~ 		n = 2
!~ 	  ELSE
!~ 	    n = 3
!~ 	  END IF
!~       !
!~       IF (SQRT(SUM(eps**two)) .GT. TOL) THEN
!~           eps_sceptral(1:n,1:n) 	  = tensor_sceptral(eps(1:n,1:n), n)
!~           !
!~           CALL principal_system(eps(1:n, 1:n), P(1:n, 1:n), W(1:n), n)
!~           !
!~           !DO i = 1, n
!~             ! Eigenwert 0 in W identifizieren und umsortieren danach
!~           !END DO
!~           ! Formel für Spannung mit Eigenprojektionen P berechnen. Sollte Spannung ergeben.
!~           !WRITE(6,*) '------------------------------------------------- '
!~           !WRITE(6,*) '----Nach d_HFED_POS_MIEHE_d_eps_e---------------- '
!~           !WRITE(6,*) '------------------------------------------------- '
!~ 	      !WRITE(6,*) 'eps: ', eps
!~ 	      !WRITE(6,*) 'eps_sceptral: ', eps_sceptral
!~ 	      !WRITE(6,*) 'P: ', P
!~ 	      !WRITE(6,*) 'W: ', W
!~ 	      !WRITE(6,*) 'n: ', n
!~ 	      !WRITE(6,*) '------------------------------------------------- '
!~       END IF
!~       !
!~       !
!~       !
!~       DO i = 1, n
!~         Proj(:,:,i) = DYADE(P(:,i),P(:,i))
!~       END DO
!~       !
!~       !DO i = 1, 3
!~       !   IF (vektor(i) == 0) THEN
!~       !      !
!~       !      ! Vertauschen der Zeilen und Spalten
!~       !      !
!~       !      IF (W(1) == 0) then
!~       !           ! Temporäre Kopie der Matrix
!~       !           P_temp = A
!~       !
!~       !           ! Zeilen tauschen (zyklisch wie beschrieben)
!~       !           A(1,:) = A_temp(2,:)
!~       !           A(2,:) = A_temp(3,:)
!~       !           A(3,:) = A_temp(1,:)
!~       !
!~       !           ! Spalten tauschen (zyklisch)
!~       !           A(:,1) = A_temp(:,2)
!~       !           A(:,2) = A_temp(:,3)
!~       !           A(:,3) = A_temp(:,1)
!~       !        END IF
!~       !   END IF
!~       !END DO
!~       !
!~       !
!~ 	  ! Weil n für den trace-Term in 2D nur bis zum vorletzten Wert geht
!~ 	  !
!~ 	  IF (n .EQ. 2) THEN
!~         d_HFEDposMieheSplit_d_eps_e_sceptral(3,3) = lambda * MacAulay(trace(eps_sceptral))
!~       END IF 
!~       !
!~       !
!~       !
!~       DO i = 1, n
!~         !
!~         !
!~ 		d_HFEDposMieheSplit_d_eps_e_sceptral = d_HFEDposMieheSplit_d_eps_e_sceptral + ( lambda * MacAulay(trace(eps_sceptral)) + two * mu * MacAulay(eps_sceptral(i,i)) ) * Proj(:,:,i)
!~ 		!
!~ 	    !
!~ 	  END DO
!~       !
!~ 	  !
!~ 	  d_HFEDposMieheSplit_NOTRANSFORMATION_d_eps_e = d_HFEDposMieheSplit_d_eps_e_sceptral
!~ 	  !
	  
      
!~     END FUNCTION d_HFEDposMieheSplit_NOTRANSFORMATION_d_eps_e

!~ !------------------------------------------------------------------------------------

!~     REAL(kind=AbqRK) FUNCTION d_HFEDnegMieheSplit_NOTRANSFORMATION_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
!~     ! d_Psi_d_eps neg Miehe et al (2010)
    
!~     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!~ 	!
!~ 	!	läuft			
!~ 	!
!~ 	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!~       USE ABQINTERFACE_PF
!~       USE FLOATNUMBERS
!~       USE MathModul
!~       USE TensorModule

!~       IMPLICIT NONE
!~       INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
!~       REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
!~       REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
!~       REAL(kind=AbqRK) :: eps_sceptral(3,3), eps_sceptral_neg(3,3), P(3,3), Proj(3,3,3), d_HFEDnegMieheSplit_d_eps_e_sceptral(3,3)
!~       REAL(kind=AbqRK) :: W(3)
!~       REAL(kind=AbqRK) :: E, nu, lambda, mu, Dim2DCheck, TOL
!~       INTEGER(kind=AbqIK) :: i, n
!~       DIMENSION d_HFEDnegMieheSplit_NOTRANSFORMATION_d_eps_e(3,3)
!~       !
!~       !
!~       TOL = 1e-10
!~       !
!~       eps_sceptral_neg = zero
!~       eps_sceptral = zero
!~       W = zero
!~       P = zero
!~       d_HFEDnegMieheSplit_d_eps_e_sceptral = zero
!~       !
!~       E               = parHFEDMatrixPhase(1)
!~       nu              = parHFEDMatrixPhase(2)
!~       !
!~       mu              = half * one / (one + nu) * E
!~       lambda          = nu / (one - two * nu) * one / (one + nu) * E
!~       !
!~       Dim2DCheck = eps(1,3)**two + eps(2,3)**two + eps(3,1)**two + eps(3,2)**two + eps(3,3)**two
!~       IF (Dim2DCheck < TOL) THEN
!~ 		n = 2
!~ 	  ELSE
!~ 	    n = 3
!~ 	  END IF
!~       !
!~       IF (SQRT(SUM(eps**two)) .GT. TOL) THEN
!~           eps_sceptral(1:n,1:n) 	  = tensor_sceptral(eps(1:n,1:n), n)
!~           !
!~           CALL principal_system(eps(1:n, 1:n), P(1:n, 1:n), W(1:n), n)
!~           !
!~           !DO i = 1, n
!~             ! Eigenwert 0 in W identifizieren und umsortieren danach
!~           !END DO
!~           ! Formel für Spannung mit Eigenprojektionen P berechnen. Sollte Spannung ergeben.
!~       END IF
!~       !
!~       !
!~       !
!~       DO i = 1, n
!~         Proj(:,:,i) = DYADE(P(:,i),P(:,i))
!~       END DO
!~       !
!~       !DO i = 1, 3
!~       !   IF (vektor(i) == 0) THEN
!~       !      !
!~       !      ! Vertauschen der Zeilen und Spalten
!~       !      !
!~       !      IF (W(1) == 0) then
!~       !           ! Temporäre Kopie der Matrix
!~       !           P_temp = A
!~       !
!~       !           ! Zeilen tauschen (zyklisch wie beschrieben)
!~       !           A(1,:) = A_temp(2,:)
!~       !           A(2,:) = A_temp(3,:)
!~       !           A(3,:) = A_temp(1,:)
!~       !
!~       !           ! Spalten tauschen (zyklisch)
!~       !           A(:,1) = A_temp(:,2)
!~       !           A(:,2) = A_temp(:,3)
!~       !           A(:,3) = A_temp(:,1)
!~       !        END IF
!~       !   END IF
!~       !END DO
!~       !
!~       !
!~ 	  ! Weil n für den trace-Term in 2D nur bis zum vorletzten Wert geht
!~ 	  !
!~ 	  IF (n .EQ. 2) THEN
!~         d_HFEDnegMieheSplit_d_eps_e_sceptral(3,3) = lambda * MacAulay_neg(trace(eps_sceptral))
!~       END IF 
!~       !
!~       !
!~       !
!~       DO i = 1, n
!~         !
!~         !
!~ 		d_HFEDnegMieheSplit_d_eps_e_sceptral = d_HFEDnegMieheSplit_d_eps_e_sceptral + ( lambda * MacAulay_neg(trace(eps_sceptral)) + two * mu * MacAulay_neg(eps_sceptral(i,i)) ) * Proj(:,:,i)
!~ 		!
!~ 	    !
!~ 	  END DO
!~       !
!~ 	  !
!~ 	  d_HFEDnegMieheSplit_NOTRANSFORMATION_d_eps_e = d_HFEDnegMieheSplit_d_eps_e_sceptral
!~ 	  !
!~ 	  !
!~ 	  !WRITE(6,*) '------------------------------------------------- '
!~       !WRITE(6,*) '----Nach d_HFED_POS_MIEHE_d_eps_e---------------- '
!~       !WRITE(6,*) '------------------------------------------------- '
!~ 	  !WRITE(6,*) 'eps: ', eps
!~ 	  !WRITE(6,*) 'eps_sceptral: ', eps_sceptral
!~ 	  !WRITE(6,*) '------------------------------------------------- '
	  
      
!~     END FUNCTION d_HFEDnegMieheSplit_NOTRANSFORMATION_d_eps_e

!~ !------------------------------------------------------------------------------------














!~ !
!~ !-----------------------------------------------------------------------------


!~     REAL(kind=AbqRK) FUNCTION tensor_sceptral(T, n)
!~     ! Sceptral representation of the tensor T; get Eigenvalues and Eigenvectors

!~       USE ABQINTERFACE_PF
!~       USE FLOATNUMBERS
!~       USE TensorModule

!~       IMPLICIT NONE
!~       INTEGER(kind=AbqIK), INTENT(IN) :: n
!~       REAL(kind=AbqRK), INTENT(IN) :: T(n,n)
!~       REAL(kind=AbqRK) :: P(n,n)
!~       REAL(kind=AbqRK) :: W(n)
!~       REAL(kind=AbqRK) :: TOL
!~       INTEGER(kind=AbqIK) :: i
!~       DIMENSION tensor_sceptral(n,n)
!~       !
!~       TOL = 1e-10
!~       W = zero
!~       P = zero
!~       tensor_sceptral = zero
!~       !
!~       IF (SQRT(SUM(T**two)) .GT. TOL) THEN
!~ 		  CALL principal_system(T, P, W, n)
!~ 		  !      
!~ 		  tensor_sceptral = zero
!~ 		  DO i = 1, n
!~ 			tensor_sceptral(i,i) = W(i)
!~ 		  END DO
!~ 		  !
!~ 	  END IF
!~ 	  !
!~       !
!~ 	  !
!~     END FUNCTION tensor_sceptral

!~ !------------------------------------------------------------------------------
!~ !------------------------------------------------------------------------------

!~     SUBROUTINE principal_system(A, P, W, n)
    
!~ 	  USE ABQINTERFACE_PF
!~       USE FLOATNUMBERS
!~       USE MathModul
!~       USE TensorModule
    
    
!~ 		IMPLICIT NONE
!~         INTEGER(kind=AbqIK), INTENT(IN) :: n
!~         REAL(kind=AbqRK), INTENT(IN) :: A(n,n)
!~         REAL(kind=AbqRK), INTENT(OUT) :: P(n,n), W(n)
!~         REAL(kind=AbqRK) :: A_dummy(n,n)
!~         REAL(kind=AbqRK), ALLOCATABLE :: work(:)
!~         INTEGER(kind=AbqIK) :: info, lwork
!~         !
!~         ! Calculates Eigenvectors and Eigenvalues of the Matrix A
!~         !
!~         ! A: Argument Matrix
!~         ! P: respective Eigenvectors (columnwise)
!~         ! P =
!~ 		!	⎡ v1_1  v2_1  v3_1 ⎤
!~ 		!	⎢ v1_2  v2_2  v3_2 ⎥
!~ 		!	⎣ v1_3  v2_3  v3_3 ⎦
!~         ! W: respective Eigenvalues (from smallest to largest, columnwise)
!~         ! For λ₁ > λ₂ > λ₃
!~         ! W(1) = λ₁  ! smallest Eigenvalue
!~ 		! W(2) = λ₂  ! next Eigenvalue
!~ 		! W(3) = λ₃  ! largest Eigenvalue
!~ 		! n: Dimension of the Matrix
!~ 		!
!~         !
!~         ! Copy Matrix A, DSYEV overwrites Matrix
!~         !
!~         P = zero
!~         W = zero
!~         A_dummy = zero
!~         !
!~         A_dummy = A

!~         ! Determine Workspace
!~         lwork = -1
!~         ALLOCATE(work(1))
!~         CALL DSYEV('V', 'U', n, A_dummy, n, W, work, lwork, info)
!~ 		!
!~ 		! WRITE STATEMENTS AUSKOMMENTIERT UM PURE SUBROUTINE ZU ERMÖGLICHEN
!~         !IF (info /= 0) THEN
!~ 		!	WRITE(*,7) 'Fehler bei der Arbeitsspeicher!'
!~         !    CALL XEXIT()
!~         !END IF
!~ 		!
!~         ! Workspace
!~         lwork = INT(work(1))
!~         DEALLOCATE(work)
!~         ALLOCATE(work(lwork))

!~ 		P = A_dummy
		
!~         ! Eigenvalues, Eigenvectors
        
!~         CALL DSYEV('V', 'U', n, P, n, W, work, lwork, info)

!~         !IF (info /= 0) THEN
!~         !    WRITE(*,7) 'Fehler bei der Eigenwertberechnung!'
!~         !    CALL XEXIT()
!~         !END IF

!~         DEALLOCATE(work)
!~     END SUBROUTINE principal_system
    
!~ !------------------------------------------------------------------------------
!~ !------------------------------------------------------------------------------

!~ 	PURE SUBROUTINE Miehe_Stiffness_Additional_Tensors(EV, EW, J_tens4, P_pos, P_neg)
	
	
!~ 	      USE ABQINTERFACE_PF
!~ 		  USE FLOATNUMBERS
!~ 		  USE MathModul
!~ 		  USE TensorModule
		
	
!~ 		IMPLICIT NONE
!~ 		REAL(kind=AbqRK), INTENT(IN) :: EV(3,3)
!~ 		REAL(kind=AbqRK), INTENT(IN) :: EW(3)
!~         REAL(kind=AbqRK), INTENT(OUT) :: J_tens4(3,3,3,3), P_pos(3,3,3,3), P_neg(3,3,3,3)
!~         REAL(kind=AbqRK) :: TOL, P_pos_SUM
!~         INTEGER(kind=AbqIK) :: a,b
!~         INTEGER(kind=AbqIK) :: i,j,k,l
!~ !		
!~ 		TOL = 1e-8

!~ 		J_tens4 = zero
!~ 		P_pos_SUM = zero
!~ 		P_pos = zero
!~ 		P_neg = zero

!~ 		! Berechnung J Tensor
!~         DO i = 1, 3
!~ 			DO j = 1, 3
!~ 				DO k = 1, 3
!~ 					DO l = 1, 3
!~ 						J_tens4 = kronecker(i,j) * kronecker(k,l)
!~ 					END DO
!~ 				END DO
!~ 			END DO
!~ 		END DO
!~ 		!
!~ 		! Berechnung der P-Projektionstensoren
!~ 		!
!~ 		! Nach Miehe (Paper von Martinez Paneda)
!~ 		! Formel aus 10.1016/j.apples.2021.100050
!~ 		!
!~ 		!
!~ 		! P_Positiv
!~ 		!
!~ 		! Cases: 
!~ 		! 
!~ 		! 1: a = b --> EW(a) = EW(b)
!~ 		! First Sum, Second Sum (With Heavyside)
!~ 		!
!~ 		! 2: a /= b;   EW(a) = EW(b)
!~ 		! Second Sum (With Heavyside)
!~ 		!
!~ 		! 3: a /= b;   EW(a) /= EW(b)
!~ 		! Second Sum (No Heavyside)
!~ 		!
!~ 		DO i = 1, 3
!~ 			DO j = 1, 3
!~ 				DO k = 1, 3
!~ 					DO l = 1, 3
!~ 						! Initialisieren für neue Komponente von P
!~ 						P_pos_SUM = zero
!~ 						DO a = 1, 3
!~ 							DO b = 1, 3
!~ 								IF (EW(a)-EW(b) .LT. TOL) THEN
!~ 								! Sind Eigenwerte gleich?
!~ 									IF (a .EQ. b) THEN
!~ 									! Case 1
!~ 										P_pos_SUM = P_pos_SUM + d_MacAulay(EW(a)) * EV(i,a)*EV(j,a)*EV(k,b)*EV(l,b) + & 
!~ 													  half * d_MacAulay(EW(a)) * &
!~ 													  EV(i,a)*EV(j,b)* (EV(k,a)*EV(l,b) + EV(k,b)*EV(l,a))	
!~ 									ELSE
!~ 									! Case 2
!~ 										P_pos_SUM = P_pos_SUM + half * d_MacAulay(EW(a)) * &
!~ 													  EV(i,a)*EV(j,b)* (EV(k,a)*EV(l,b) + EV(k,b)*EV(l,a))	
!~ 									END IF
										
!~ 								ELSE
!~ 								! Case 3
!~ 									P_pos_SUM = P_pos_SUM + half * (MacAulay(EW(a)) - MacAulay(EW(b)))/(EW(a)-EW(b)) * &
!~ 												  EV(i,a)*EV(j,b)* (EV(k,a)*EV(l,b) + EV(k,b)*EV(l,a))	
!~ 								END IF
!~ 							END DO
!~ 						END DO
!~ 						P_pos(i,j,k,l) = P_pos_SUM
!~ 					END DO
!~ 				END DO
!~ 			END DO
!~ 		END DO
		
		
!~ 		! P_Negativ
!~ 		P_neg = IdentityTens4(3) - P_pos
		
	
!~ 	END SUBROUTINE Miehe_Stiffness_Additional_Tensors
	
!~ !------------------------------------------------------------------------------
!~ !------------------------------------------------------------------------------
    
!~ 	SUBROUTINE retransformation(A, P, n)
	
!~ 	  USE ABQINTERFACE_PF
!~ 	  USE FLOATNUMBERS
!~ 	  USE MathModul
!~ 	  USE TensorModule
		  
		  
!~ 		INTEGER(kind=AbqIK), INTENT(IN) :: n
!~ 		REAL(kind=AbqRK), INTENT(IN) :: P(n,n)
!~ 		REAL(kind=AbqRK), INTENT(INOUT) :: A(n,n)
!~ 		REAL(kind=AbqRK), DIMENSION(n,n) :: P_inv, retransformations_Matrix
		
!~ 		P_inv = zero
!~ 		retransformations_Matrix = zero
		
!~ 		!WRITE(6,*) '___________________________'
!~ 		!WRITE(6,*) '______RETRANSFORMATION_____'
!~ 		!WRITE(6,*) '___________________________'
!~ 		!WRITE(6,*) 'A vorher: ', A
!~ 		!WRITE(6,*) 'P vorher: ', P
!~ 		!WRITE(6,*) 'n Dimension', n
!~ 		!WRITE(6,*) '___________________________'
		
!~ 		CALL get_symm_matrixinverse(n,P,P_inv)
		
!~ 		!WRITE(6,*) '___________________________'
!~ 		!WRITE(6,*) '______VOR      INVERSE_____'
!~ 		!WRITE(6,*) '___________________________'
!~ 		!WRITE(6,*) 'A vorher : ', A
!~ 		!WRITE(6,*) 'P : ', P
!~ 		!WRITE(6,*) 'P Inverse: ', P_inv
!~ 		!WRITE(6,*) 'n Dimension', n
!~ 		!WRITE(6,*) '___________________________'
		 
!~ 		! Rücktransformation: A = P * A' * P_inv
!~ 		retransformations_Matrix = MATMUL(MATMUL(P, A), P_inv)
		
!~ 		!WRITE(6,*) '___________________________'
!~ 		!WRITE(6,*) '______Nach     INVERSE_____'
!~ 		!WRITE(6,*) '___________________________'
!~ 		!WRITE(6,*) 'A vorher : ', A
!~ 		!WRITE(6,*) 'P : ', P
!~ 		!WRITE(6,*) 'P Inverse: ', P_inv
!~ 		!WRITE(6,*) 'n Dimension: ', n
!~ 		!WRITE(6,*) 'RETRANSMATRIX: ', retransformations_Matrix
!~ 		!WRITE(6,*) '___________________________'

!~ 		A = retransformations_Matrix
		
!~ 		!WRITE(6,*) '___________________________'
!~ 		!WRITE(6,*) '___Nach RETRANSFORMATION___'
!~ 		!WRITE(6,*) '___________________________'
!~ 		!WRITE(6,*) 'RETRANSFORMATIONSMATRIX', retransformations_Matrix
!~ 		!WRITE(6,*) '___________________________'
		 
!~ 	END SUBROUTINE retransformation
	
!~ !------------------------------------------------------------------------------
!~ !------------------------------------------------------------------------------

!~ 	SUBROUTINE d_Eigenwerte_d_eps(eps, d_EW_d_eps)
!~ 	  USE ABQINTERFACE_PF
!~ 	  USE FLOATNUMBERS
!~ 	  USE MathModul
!~ 	  USE TensorModule
	  
!~ 		IMPLICIT NONE
!~ 		REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
!~ 		REAL(kind=AbqRK) :: EV(3,3), T(3,3,3), I4sym(3,3,3,3)
!~ 		REAL(kind=AbqRK) :: EW(3)
!~         REAL(kind=AbqRK), INTENT(OUT) :: d_EW_d_eps(3,3,3)
!~         REAL(kind=AbqRK) :: TOL
!~         INTEGER(kind=AbqIK) :: i, n
        
!~         T = zero
!~ 		d_EW_d_eps = zero
	
!~ 		!WRITE(6,*) '----------------------------'
!~ 		!WRITE(6,*) 'd_EW_d_eps_Funktion: '
!~ 		!WRITE(6,*) '----------------------------'
!~ 		!WRITE(6,*) 'eps: ', eps
!~ 		!WRITE(6,*) '----------------------------'
	
!~ 		CALL principal_system(eps, EV, EW, n)
		
!~ 		!WRITE(6,*) '----Nach Princ_syst---------'
!~ 		!WRITE(6,*) 'd_EW_d_eps_Funktion: '
!~ 		!WRITE(6,*) '----------------------------'
!~ 		!WRITE(6,*) 'eps: ', eps
!~ 		!WRITE(6,*) 'EW: ', EW
!~ 		!WRITE(6,*) 'EV: ', EV
!~ 		!WRITE(6,*) '----------------------------'
	
!~ 		I4sym = IdentityTens4sym(3)
			
!~ 		!WRITE(6,*) '----Vor Schleife--------'
!~ 		!WRITE(6,*) 'd_EW_d_eps_Funktion: '
!~ 		!WRITE(6,*) '----------------------------'
!~ 		!WRITE(6,*) 'I4Sym: ', I4sym
!~ 		!WRITE(6,*) '----------------------------'
			
!~ 		DO i = 1,3
!~ 			T = zero
!~ 			T = SINGLECONTRACTIONFourOne(I4sym, EV(:,i))
!~ 			d_EW_d_eps(:,:,i) = SINGLECONTRACTIONThreeOne(T, EV(:,i))
!~ 			!
!~ 			!WRITE(6,*) '----WÄHREND SCHLEIFE--------'
!~ 			!WRITE(6,*) 'd_EW_d_eps_Funktion: '
!~ 			!WRITE(6,*) '----------------------------'
!~ 			!WRITE(6,*) 'Schleifeniteration: ', i
!~ 			!WRITE(6,*) 'SINGLECONTRACTIONFourOne(I4sym, EV(:,i)): ', T
!~ 			!WRITE(6,*) 'SINGLECONTRACTIONThreeOne(T, EV(:,i))', d_EW_d_eps(:,:,i)
!~ 			!WRITE(6,*) '----------------------------'
!~ 		END DO
		
!~ 	END SUBROUTINE d_Eigenwerte_d_eps
	
!~ !------------------------------------------------------------------------------

!~ 	SUBROUTINE d_EV_d_eps(eps, n, P, dPdE)
	  
!~ 	  USE ABQINTERFACE_PF
!~ 	  USE FLOATNUMBERS
!~ 	  USE MathModul
!~ 	  USE TensorModule
		
!~ 		IMPLICIT NONE
!~ 		REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
!~ 		INTEGER(kind=AbqIK), INTENT(IN) :: n
!~ 		REAL(kind=AbqRK) :: eps0(3,3), eps1(3,3), eps2(3,3), eps_sceptral(3,3), Ident(3,3)
!~ 		REAL(kind=AbqRK) :: EW(3), EV(3,3)
!~         REAL(kind=AbqRK), INTENT(OUT) :: P(3,3,3), dPdE(3,3,3,3,3)
!~         REAL(kind=AbqRK) :: TOL
!~         INTEGER(kind=AbqIK) :: EWCASE
!~         INTEGER(kind=AbqIK) :: i
!~         REAL(kind=AbqRK) :: EWi, EW_double
!~         REAL(kind=AbqRK) :: l1, l2, l3
!~         REAL(kind=AbqRK) :: T1, T2
!~         REAL(kind=AbqRK) :: D1, D2, D3
!~         REAL(kind=AbqRK) :: Q1
!~         REAL(kind=AbqRK) :: v100, v101, v102, v110, v111, v112, v120, v121, v122
!~         REAL(kind=AbqRK) :: v200, v201, v202, v210, v211, v212, v220, v221, v222
!~         REAL(kind=AbqRK) :: v300, v301, v302, v310, v311, v312, v320, v321, v322
!~         REAL(kind=AbqRK) :: eps0eps0(3,3,3,3), eps0eps1(3,3,3,3), eps0eps2(3,3,3,3), eps1eps0(3,3,3,3), eps1eps1(3,3,3,3), eps1eps2(3,3,3,3), eps2eps0(3,3,3,3), eps2eps1(3,3,3,3), eps2eps2(3,3,3,3)
!~         REAL(kind=AbqRK) :: dP1dE, dP2dE, dP3dE
        
!~         P = zero
!~         EW = zero
!~         eps_sceptral = zero
!~         !
!~         !
!~         TOL = 1e-8
!~         !
!~         !
!~ 		Ident = Identity(3)
!~ 		!
!~ 		!
!~ 		eps0 = Identity(3)
!~ 		eps1 = eps
!~ 		eps2 = MATMUL(eps,eps)
!~ 		!
!~ 		! für Itskov 7.56 unten
!~ 		!
!~ 		eps0eps0 = DYADETwoTwo(eps0, eps0)
!~ 		eps0eps1 = DYADETwoTwo(eps0, eps1)
!~ 		eps0eps2 = DYADETwoTwo(eps0, eps2)
!~ 		eps1eps0 = DYADETwoTwo(eps1, eps0)
!~ 		eps1eps1 = DYADETwoTwo(eps1, eps1)
!~ 		eps1eps2 = DYADETwoTwo(eps1, eps2)
!~ 		eps2eps0 = DYADETwoTwo(eps2, eps0)
!~ 		eps2eps1 = DYADETwoTwo(eps2, eps1)
!~ 		eps2eps2 = DYADETwoTwo(eps2, eps2)
		
!~ 		!WRITE(6,*) '-------------------------------------------'
!~ 		!WRITE(6,*) '--------EPS EPS 0 , 1, 2, 		 ----------'
!~ 		!WRITE(6,*) '-------------------------------------------'
!~ 		!WRITE(6,*) 'Eps 0: ', eps0
!~ 		!WRITE(6,*) 'Eps 1: ', eps1
!~ 		!WRITE(6,*) 'Eps 2: ', eps2
!~ 		!WRITE(6,*) 'Eps 0 Eps 0 : ', eps0eps0
!~ 		!WRITE(6,*) 'Eps 0 Eps 1 (1,1,:,:): ', eps0eps1(1,1,:,:)
!~ 		!WRITE(6,*) '-------------------------------------------'
		
!~ 		CALL principal_system(eps, EV, EW, n)
!~ 		!
!~ 		!WRITE(6,*) '-------------------------'
!~ 		!WRITE(6,*) '--principal_system--'
!~ 		!WRITE(6,*) '-------------------------'
!~ 		!WRITE(6,*) 'eps: ', eps
!~ 		!WRITE(6,*) 'Eigenvektoren: ', EV
!~ 		!WRITE(6,*) 'Eigenwerte: ', EW
!~ 		!WRITE(6,*) '-------------------------'
!~ 		!
!~ 		!
!~ 		!----------------------------------------------
!~ 		! Test welcher Fall vorliegt
!~ 		!----------------------------------------------
!~ 		!
		
!~ 		IF (ABS(EW(1) - EW(2)) .LT. TOL .AND. ABS(EW(2) - EW(3)) .LT. TOL) THEN
!~ 			EWCASE = 4   ! Alle drei Werte sind gleich
!~ 			!WRITE(6,*) 'Aufgerufen ALLE DREI WERTE GLEICH'
!~ 		ELSE IF (ABS(EW(1) - EW(2)) .LT. TOL) THEN
!~ 			EWCASE = 2   ! EW(1) und EW(2) sind gleich
!~ 			!WRITE(6,*) 'Aufgerufen EW(1) und EW(2) sind gleich'
!~ 		ELSE IF (ABS(EW(2) - EW(3)) .LT. TOL) THEN
!~ 			EWCASE = 3   ! EW(2) und EW(3) sind gleich
!~ 			!WRITE(6,*) 'Aufgerufen EW(2) und EW(3) sind gleich'
!~ 		ELSE
!~ 			EWCASE = 1   ! Alle verschieden
!~ 			!WRITE(6,*) 'Aufgerufen Alle verschieden'
!~ 		END IF
!~ 		!
!~ 		!
!~ 		!WRITE(6,*) '-------------------------'
!~ 		!WRITE(6,*) '--EW + Case distiction---'
!~ 		!WRITE(6,*) '-------------------------'
!~ 		!WRITE(6,*) 'Eigenwerte: ', EW
!~ 		!WRITE(6,*) 'Case: ', EWCASE
!~ 		!WRITE(6,*) '-------------------------'
!~ 		!
!~ 		!-----------------------------------------------
!~ 		!	Eigenprojections
!~ 		!-----------------------------------------------
!~ 		!
!~ 		! 4. Case 
!~ 		!
!~ 		IF (EWCASE .EQ. 4) THEN
!~ 			P(:,:,1) = Ident
!~ 		END IF
!~ 		!
!~ 		! 3. Case
!~ 		!
!~ 		IF (EWCASE .EQ. 3) THEN
!~ 			P(:,:,2) = DYADE(EV(:,3), EV(:,3))
!~ 			P(:,:,1) = Ident - P(:,:,2)
!~ 		END IF
!~ 		!
!~ 		! 2. Case
!~ 		!
!~ 		IF (EWCASE .EQ. 2) THEN
!~ 			P(:,:,1) = DYADE(EV(:,1), EV(:,1))
!~ 			P(:,:,2) = Ident - P(:,:,1)
!~ 		END IF
!~ 		!
!~ 		! 1. Case 
!~ 		!
!~ 		IF (EWCASE .EQ. 1) THEN
!~ 			DO i = 1, 3
!~ 				P(:,:,i) = DYADE(EV(:,i), EV(:,i))
!~ 			END DO
!~ 		END IF
!~ 		!
!~ 		!WRITE(6,*) '-------------------------'
!~ 		!WRITE(6,*) '----Eigenprojections-----'
!~ 		!WRITE(6,*) '-------------------------'
!~ 		!WRITE(6,*) 'Eigenprojections: ', P
!~ 		!WRITE(6,*) '-------------------------'
!~ 		!
!~ 		!-----------------------------------------------
!~ 		!
!~ 		!-----------------------------------------------
!~ 		!	Ableitungen bestimmen
!~ 		!-----------------------------------------------
!~ 		!
!~ 		!-----------------------------------------------
!~ 		!
!~ 		! 
!~ 		!-----------------------------------------------
!~ 		!
!~ 		! 4. Case 
!~ 		!
!~ 		IF (EWCASE .EQ. 4) THEN
!~ 			P(:,:,1) = Ident
!~ 		END IF
!~ 		!
!~ 		! 3. Case
!~ 		!
!~ 		IF (EWCASE .EQ. 3) THEN
!~ 			P(:,:,2) = DYADE(EV(:,3), EV(:,3))
!~ 			P(:,:,1) = Ident - P(:,:,2)
!~ 		END IF
!~ 		!
!~ 		! 2. Case
!~ 		!
!~ 		IF (EWCASE .EQ. 2) THEN
!~ 			P(:,:,1) = DYADE(EV(:,1), EV(:,1))
!~ 			P(:,:,2) = Ident - P(:,:,1)
!~ 		END IF
!~ 		!
!~ 		! 1. Case 
!~ 		!
!~ 		IF (EWCASE .EQ. 1) THEN
!~ 			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!~  			! Ändernde Variablen über i,j,k: 
!~ 			! l1, l2, l3
!~ 			! D1, D2, D3
!~ 			! T1, T2
!~ 			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!~ 			! vipq -> geprüft für Case 1
!~ 			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!~ 			! Eigenwerte: 
!~ 			l1 = EW(1)
!~ 			l2 = EW(2)
!~ 			l3 = EW(3)
!~ 			!
!~ 			! Quotienten
!~ 			!
!~ 			D1=(l1-l2)*(l1-l3)
!~ 			D2=(l2-l3)*(l2-l1)
!~ 			D3=(l3-l1)*(l3-l2)
!~ 			!
!~ 			! i=1. j=2, k=3
!~ 			!
!~ 			T1=D3*(l1-l2)**three
!~ 			T2=D2*(l1-l3)**three
!~ 			v100=-two*l1*l2*l3*(l3/T1+l2/T2)
!~ 			v101=l3*(l1*(l2+l3)+l2*(l1+l3))/T1+l2*(l1*(l2+l3)+l3*(l1+l2))/T2
!~ 			v110=v101
!~ 			v102=-l3*(l1+l2)/T1-l2*(l1+l3)/T2
!~ 			v120=v102
!~ 			v111=-two*(l2+l3)*((l1+l3)/T1+(l1+l2)/T2)
!~ 			v112=(l1+l2+two*l3)/T1+(l1+two*l2+l3)/T2
!~ 			v121=v112
!~ 			v122=-two/T1-two/T2
!~ 			!
!~ 			!
!~ 			!
!~ 			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!~ 			! Eigenwerte: 
!~ 			l1 = EW(2)
!~ 			l2 = EW(3)
!~ 			l3 = EW(1)
!~ 			!
!~ 			! Quotienten
!~ 			!
!~ 			D1=(l1-l2)*(l1-l3)
!~ 			D2=(l2-l3)*(l2-l1)
!~ 			D3=(l3-l1)*(l3-l2)
!~ 			!
!~ 			! i=2. j=3, k=1
!~ 			!
!~ 			T1=D3*(l1-l2)**three
!~ 			T2=D2*(l1-l3)**three
!~ 			v200=-two*l1*l2*l3*(l3/T1+l2/T2)
!~ 			v201=l3*(l1*(l2+l3)+l2*(l1+l3))/T1+l2*(l1*(l2+l3)+l3*(l1+l2))/T2
!~ 			v210=v201
!~ 			v202=-l3*(l1+l2)/T1-l2*(l1+l3)/T2
!~ 			v220=v202
!~ 			v211=-two*(l2+l3)*((l1+l3)/T1+(l1+l2)/T2)
!~ 			v212=(l1+l2+two*l3)/T1+(l1+two*l2+l3)/T2
!~ 			v221=v212
!~ 			v222=-two/T1-two/T2
!~ 			!
!~ 			!
!~ 			!
!~ 			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!~ 			! Eigenwerte: 
!~ 			l1 = EW(3)
!~ 			l2 = EW(1)
!~ 			l3 = EW(2)
!~ 			!
!~ 			! Quotienten
!~ 			!
!~ 			D1=(l1-l2)*(l1-l3)
!~ 			D2=(l2-l3)*(l2-l1)
!~ 			D3=(l3-l1)*(l3-l2)
!~ 			!
!~ 			! i=3. j=1, k=2
!~ 			!
!~ 			T1=D3*(l1-l2)**three
!~ 			T2=D2*(l1-l3)**three
!~ 			v300=-two*l1*l2*l3*(l3/T1+l2/T2)
!~ 			v301=l3*(l1*(l2+l3)+l2*(l1+l3))/T1+l2*(l1*(l2+l3)+l3*(l1+l2))/T2
!~ 			v310=v301
!~ 			v302=-l3*(l1+l2)/T1-l2*(l1+l3)/T2
!~ 			v320=v302
!~ 			v311=-two*(l2+l3)*((l1+l3)/T1+(l1+l2)/T2)
!~ 			v312=(l1+l2+two*l3)/T1+(l1+two*l2+l3)/T2
!~ 			v321=v312
!~ 			v322=-two/T1-two/T2
!~ 			!
!~ 			dPdE(:,:,:,:,1)=v100*eps0eps0+v101*eps0eps1+v102*eps0eps2+v110*eps1eps0+v111*eps1eps1+v112*eps1eps2+v120*eps2eps0+v121*eps2eps1+v122*eps2eps2
!~ 			dPdE(:,:,:,:,2)=v200*eps0eps0+v201*eps0eps1+v202*eps0eps2+v210*eps1eps0+v211*eps1eps1+v212*eps1eps2+v220*eps2eps0+v221*eps2eps1+v222*eps2eps2
!~ 			dPdE(:,:,:,:,3)=v300*eps0eps0+v301*eps0eps1+v302*eps0eps2+v310*eps1eps0+v311*eps1eps1+v312*eps1eps2+v320*eps2eps0+v321*eps2eps1+v322*eps2eps2	
!~ 		!
!~ 		!-------------------------------------------------------
!~ 		!	Case 2
!~ 		!-------------------------------------------------------
!~ 		!
!~ 		ELSE IF (EWCASE .EQ. 2) THEN
			
!~ 			EW_double = EW(1)
!~ 			EWi = EW(2)
			
!~ 			Q1=(EWi-EW_double)**three
!~ 			v200=-two*EW_double*EWi/Q1
!~ 			v201=(EWi-EW_double)/Q1
!~ 			v210=v201
!~ 			v211=-two/Q1
!~ 			v202=zero
!~ 			v220=zero
!~ 			v212=zero
!~ 			v221=zero
!~ 			v222=zero
!~ 			!
!~ 			! (Itskov 7.56, but with supersymmetric tensor product)
!~ 			!
!~ 			dPdE(:,:,:,:,2)=v200*eps0eps0+v201*eps0eps1+v202*eps0eps2+v210*eps1eps0+v211*eps1eps1+v212*eps1eps2+v220*eps2eps0+v221*eps2eps1+v222*eps2eps2
!~ 			dPdE(:,:,:,:,1)=-dP2dE	! P1=Ident-P2
!~ 			dPdE(:,:,:,:,3)=zero
!~ 		!
!~ 		!-------------------------------------------------------
!~ 		!	Case 3
!~ 		!-------------------------------------------------------
!~ 		!
!~ 		ELSE IF (EWCASE .EQ. 3) THEN
		
!~ 			EW_double = EW(2)		
!~ 			EWi = EW(1)
		
!~ 			Q1=(EWi-EW_double)**three
!~ 			v100=-two*EW_double*EWi/Q1
!~ 			v101=(EWi-EW_double)/Q1
!~ 			v110=v101
!~ 			v111=-two/Q1
!~ 			v102=zero
!~ 			v120=zero
!~ 			v112=zero
!~ 			v121=zero
!~ 			v122=zero
!~ 			!
!~ 			! (Itskov 7.56, but with supersymmetric tensor product)
!~ 			!
!~ 			dPdE(:,:,:,:,1)=v100*eps0eps0+v101*eps0eps1+v102*eps0eps2+v110*eps1eps0+v111*eps1eps1+v112*eps1eps2+v120*eps2eps0+v121*eps2eps1+v122*eps2eps2
!~ 			dPdE(:,:,:,:,2)=-dP1dE	! P2=Ident-P1
!~ 			dPdE(:,:,:,:,3)=zero
!~ 		!
!~ 		!----------------------------------------------------
!~ 		!	 Case 4
!~ 		!----------------------------------------------------
!~ 		!
!~ 		!WRITE(6,*) 'EWCASE 4: '
!~ 		!
!~ 		!WRITE(6,*) 'Vorher: ', P(:,:,:)
!~ 		!WRITE(6,*) 'Vorher: ', dPdE(:,:,:,:,:)
!~ 		!
!~ 		ELSE IF (EWCASE .EQ. 4) THEN
!~ 			dPdE(:,:,:,:,:) = zero
!~ 		ELSE
!~ 			WRITE(6,*) "No valid Case in derivative d_EV_d_eps. EXIT", EWCASE
!~ 			CALL XEXIT()
!~ 		END IF
!~ 		!
!~ 		!WRITE(6,*) 'Nachher: ', P(:,:,:)
!~ 		!WRITE(6,*) 'Nachher: ', dPdE(:,:,:,:,:)
!~ 		!
!~ 	END SUBROUTINE d_EV_d_eps

!~ !------------------------------------------------------------------------------
!~ !------------------------------------------------------------------------------
    
END MODULE SplitModulePF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


