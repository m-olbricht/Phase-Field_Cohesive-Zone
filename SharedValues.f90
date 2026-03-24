!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! tensor modules
!
! Martin Olbricht, TU Bergakademie Freiberg, 13.05.2025
!
! 13.05.2025: first implementation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE SharedValues
  USE ABQINTERFACE_PF
  IMPLICIT NONE
  
  ! Beispielvariablen (können jederzeit erweitert werden)
  PUBLIC :: Incrementnumber, Elementnumber, Integrationpointnumber
  PUBLIC :: is_print_elem_sub, is_print_inc_sub, is_print_ip_sub
  PUBLIC :: is_print_elem_subi, is_print_inc_subi, is_print_ip_subi
  
  INTEGER(kind=AbqIK) :: Incrementnumber
  INTEGER(kind=AbqIK) :: Elementnumber
  INTEGER(kind=AbqIK) :: Integrationpointnumber
  
  LOGICAL :: is_print_elem_sub
  LOGICAL :: is_print_inc_sub
  LOGICAL :: is_print_ip_sub
  
  LOGICAL :: is_print_elem_subi
  LOGICAL :: is_print_inc_subi
  LOGICAL :: is_print_ip_subi
  
END MODULE SharedValues
