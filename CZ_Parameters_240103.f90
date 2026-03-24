!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! UMAT material subroutine of cohesive zone model: parameters
!
! Stephan Roth, TU Bergakademie Freiberg, 03.01.2024
!
! 03.01.2024: dummy
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CZ_Parameters

  USE ABQINTERFACE

  IMPLICIT NONE

  ! global definition of required number of parameters and solution dependent state variables
  ! and index of 2D element thickness in props
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: numSDV=30, numMatPar=10, thicknessIndex=7

END MODULE CZ_Parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

