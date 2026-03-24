!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! UMAT material subroutine of phase-field model: parameters
!
! Stephan Roth, TU Bergakademie Freiberg, 30.07.2020
!
! 30.07.2020: Multi-phase multi-component
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PhaseField_Parameters

  USE ABQINTERFACE_PF

  IMPLICIT NONE

  ! global definition of required number of parameters and solution dependent state variables
  ! and index of 2D element thickness in props
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: numSDV=20, numMatPar=10, thicknessIndex=7

END MODULE PhaseField_Parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

