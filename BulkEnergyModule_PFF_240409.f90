!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! bulk energy module
!
! Martin Olbricht, TU Bergakademie Freiberg, 12.01.2024
!
! 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE BulkEnergyModule

USE ABQINTERFACE_PF
USE FLOATNUMBERS
USE FreeEnergyModule
USE AliasModulePF

  IMPLICIT NONE

!  PUBLIC :: 

  CONTAINS

!------------------------------------------------------------------------------------!

    PURE REAL(kind=AbqRK) FUNCTION HFEDtens_H(eps,nHFEDpar,parHFEDMatrixPhase)
    ! bulk energy density
      
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      !
      HFEDtens_H = HFEDpos(eps,nHFEDpar,parHFEDMatrixPhase)
      
      
    END FUNCTION HFEDtens_H
    
!------------------------------------------------------------------------------------!

    PURE REAL(kind=AbqRK) FUNCTION bulkED(eps,damage,nHFEDpar,parHFEDMatrixPhase,prop_df_alphaPF,prop_df_betaPF)
    ! bulk energy density

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK), INTENT(IN) :: damage
      REAL(kind=AbqRK), INTENT(IN) :: prop_df_alphaPF,prop_df_betaPF
      REAL(kind=AbqRK) :: HFEDtens, HFEDcomp
      REAL(kind=AbqRK) :: degD
      !
	  !
      HFEDtens = HFEDpos(eps,nHFEDpar,parHFEDMatrixPhase)
      HFEDcomp = HFEDneg(eps,nHFEDpar,parHFEDMatrixPhase)
      degD = Degradation(damage,prop_df_betaPF,prop_df_alphaPF)
      !
      bulkED = degD*HFEDtens + HFEDcomp
      
      
    END FUNCTION bulkED
    
!------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------!
!                               first derivatives
!------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------!

    PURE REAL(kind=AbqRK) FUNCTION d_bulkED_d_eps(eps,damage,nHFEDpar,parHFEDMatrixPhase,prop_df_alphaPF,prop_df_betaPF)
    ! derivative of bulk energy density w.r.t. strain

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK), INTENT(IN) :: damage
      REAL(kind=AbqRK), INTENT(IN) :: prop_df_alphaPF,prop_df_betaPF
      REAL(kind=AbqRK) :: degD
      REAL(kind=AbqRK) :: d_HFEDtens_d_eps_e(3,3), d_HFEDcomp_d_eps_e(3,3)
      DIMENSION d_bulkED_d_eps(3,3)
      !
      d_HFEDtens_d_eps_e = d_HFEDpos_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
      d_HFEDcomp_d_eps_e = d_HFEDneg_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
      !
      degD = Degradation(damage,prop_df_betaPF,prop_df_alphaPF)
      !
      d_bulkED_d_eps = degD*d_HFEDtens_d_eps_e + d_HFEDcomp_d_eps_e
      
      
    END FUNCTION d_bulkED_d_eps

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_bulkED_d_damage_H(eps,damage,nHFEDpar,parHFEDMatrixPhase,prop_df_alphaPF,prop_df_betaPF,H)
    ! derivative of bulk energy density w.r.t. damage

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK), INTENT(IN) :: damage,H
      REAL(kind=AbqRK), INTENT(IN) :: prop_df_alphaPF,prop_df_betaPF
      REAL(kind=AbqRK) :: d_degD_d_damage_H, HFEDtens
      !
      HFEDtens = H
      d_degD_d_damage_H = d_Degradation_d_damage(damage,prop_df_betaPF,prop_df_alphaPF)
      !
      d_bulkED_d_damage_H = d_degD_d_damage_H*HFEDtens
      
      
    END FUNCTION d_bulkED_d_damage_H
    

!------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------!
!                               second derivatives
!------------------------------------------------------------------------------------!
!------------------------------------------------------------------------------------!


    PURE REAL(kind=AbqRK) FUNCTION d_bulkED_d_eps_d_eps(eps,damage,nHFEDpar,parHFEDMatrixPhase,prop_df_alphaPF,prop_df_betaPF)
    ! doubled partial derivative of bulk energy density w.r.t. strain strain

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK), INTENT(IN) :: damage
      REAL(kind=AbqRK), INTENT(IN) :: prop_df_alphaPF,prop_df_betaPF
      REAL(kind=AbqRK) :: d_HFEDtens_d_eps_e_d_eps_e(3,3,3,3), d_HFEDcomp_d_eps_e_d_eps_e(3,3,3,3)
      REAL(kind=AbqRK) :: degD
      DIMENSION d_bulkED_d_eps_d_eps(3,3,3,3)
      !
      !
      d_HFEDtens_d_eps_e_d_eps_e = d_HFEDpos_d_eps_e_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
      d_HFEDcomp_d_eps_e_d_eps_e = d_HFEDneg_d_eps_e_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
      !
      degD = Degradation(damage,prop_df_betaPF,prop_df_alphaPF)
      !
      !
      d_bulkED_d_eps_d_eps = degD*d_HFEDtens_d_eps_e_d_eps_e + d_HFEDcomp_d_eps_e_d_eps_e
      

    END FUNCTION d_bulkED_d_eps_d_eps

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_bulkED_d_eps_d_damage(eps,damage,nHFEDpar,parHFEDMatrixPhase,prop_df_alphaPF,prop_df_betaPF)
    ! partial derivative of the bulk energy w.r.t. strain and damage parameter
    
      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK), INTENT(IN) :: damage
      REAL(kind=AbqRK), INTENT(IN) :: prop_df_alphaPF,prop_df_betaPF
      REAL(kind=AbqRK) :: d_degD_d_damage
      REAL(kind=AbqRK) :: d_HFEDtens_d_eps_e(3,3), d_HFEDcomp_d_eps_e(3,3)
      DIMENSION d_bulkED_d_eps_d_damage(3,3)
      !
      d_HFEDtens_d_eps_e = d_HFEDpos_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
      d_HFEDcomp_d_eps_e = d_HFEDneg_d_eps_e(eps,nHFEDpar,parHFEDMatrixPhase)
      !
      d_degD_d_damage = d_Degradation_d_damage(damage,prop_df_betaPF,prop_df_alphaPF)
      !
      d_bulkED_d_eps_d_damage = d_degD_d_damage*d_HFEDtens_d_eps_e
      
      
    END FUNCTION d_bulkED_d_eps_d_damage

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_bulkED_d_damage_d_damage_H(eps,damage,nHFEDpar,parHFEDMatrixPhase,prop_df_alphaPF,prop_df_betaPF,H)
    ! partial derivative of the bulk energy w.r.t. damage and damage
    !

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar
      REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrixPhase(nHFEDpar)
      REAL(kind=AbqRK), INTENT(IN) :: eps(3,3)
      REAL(kind=AbqRK), INTENT(IN) :: damage,H
      REAL(kind=AbqRK), INTENT(IN) :: prop_df_alphaPF,prop_df_betaPF
      REAL(kind=AbqRK) :: d_degD_d_damage_d_damage_H, HFEDtens
      !
      HFEDtens = H
      d_degD_d_damage_d_damage_H = d_Degradation_d_damage_d_damage(damage,prop_df_betaPF,prop_df_alphaPF)
      !
      d_bulkED_d_damage_d_damage_H = d_degD_d_damage_d_damage_H * HFEDtens

    END FUNCTION d_bulkED_d_damage_d_damage_H

!------------------------------------------------------------------------------------!

END MODULE BulkEnergyModule
