!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ABAQUS user element library for phase field elements
!
! Stephan Roth, TU Bergakademie Freiberg, 30.07.2020
!
! 30.07.2020: Multi-phase multi-component
! 15.06.2021: concentrations independent on phases
! 28.02.2024: Newton-Cotes quadrature rule for 2D
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1D-model, line, quadratic shape functions, reduced integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PF_1D_Line_Quadratic_Reduced
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE ShapeFunc1D_Line_Quad_ABA, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr1D_Line2GP, ONLY: NGP, GPWeight, GPPos
  USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE PhaseField_module, ONLY: CheckMaterialParameters, umat => umatPF
  USE BMatrixJDet
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=1, NSHR=0, NDI=1
  CONTAINS
    INCLUDE "PFUEL_PFF_240103.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1D-model, line, linear and cubic shape functions, full integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PF_1D_Line_LinCubic
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE ShapeFunc1D_Line_Lin, ONLY: ShapeFunc
  USE ShapeFunc1D_Line_Cubic, ONLY: ShapeFuncDeriv
  USE Integr1D_Line3GP, ONLY: NGP, GPWeight, GPPos
  USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE PhaseField_module, ONLY: CheckMaterialParameters, umat => umatPF
  USE BMatrixJDet
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=1, NSHR=0, NDI=1
  CONTAINS
    INCLUDE "PFUEL_PFF_240103.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1D-model, line, linear shape functions, reduced integration, spheri-symmetric
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PF_1D_Line_Linear_Reduced_Spherisym
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE ShapeFunc1D_Line_Lin, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr1D_Line1GP, ONLY: NGP, GPWeight, GPPos
  USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE PhaseField_module, ONLY: CheckMaterialParameters, umat => umatPF
  USE BMatrixJDet
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=1, NSHR=0, NDI=3
  CONTAINS
    INCLUDE "PFUEL_PFF_240103.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1D-model, line, quadratic shape functions, reduced integration, spheri-symmetric
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PF_1D_Line_Quadratic_Reduced_Spherisym
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE ShapeFunc1D_Line_Quad_ABA, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr1D_Line2GP, ONLY: NGP, GPWeight, GPPos
  USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE PhaseField_module, ONLY: CheckMaterialParameters, umat => umatPF
  USE BMatrixJDet
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=1, NSHR=0, NDI=3
  CONTAINS
    INCLUDE "PFUEL_PFF_240103.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1D-model, line, quadratic shape functions, full integration, spheri-symmetric
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PF_1D_Line_Quadratic_Spherisym
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE ShapeFunc1D_Line_Quad_ABA, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr1D_Line3GP, ONLY: NGP, GPWeight, GPPos
  USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE PhaseField_module, ONLY: CheckMaterialParameters, umat => umatPF
  USE BMatrixJDet
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=1, NSHR=0, NDI=3
  CONTAINS
    INCLUDE "PFUEL_PFF_240103.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D-model, rectangle, linear shape functions, reduced integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PF_2D_Rect_Linear_Reduced
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE ShapeFunc2D_Square_Lin, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr2D_Square1GP, ONLY: NGP, GPWeight, GPPos
  USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE PhaseField_module, ONLY: CheckMaterialParameters, umat => umatPF
  USE BMatrixJDet
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=2, NSHR=1, NDI=3
  CONTAINS
    INCLUDE "PFUEL_PFF_240103.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D-model, rectangle, linear shape functions, full integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PF_2D_Rect_Linear_Full
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE ShapeFunc2D_Square_Lin, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr2D_Square4GP, ONLY: NGP, GPWeight, GPPos
  USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE PhaseField_module, ONLY: CheckMaterialParameters, umat => umatPF
  USE BMatrixJDet
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=2, NSHR=1, NDI=3
  CONTAINS
    INCLUDE "PFUEL_PFF_240103.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D-model, rectangle, quadratic shape functions, reduced integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PF_2D_Rect_Quadratic_Reduced
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE ShapeFunc2D_Square_Quadr, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr2D_Square4GP, ONLY: NGP, GPWeight, GPPos
  USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE PhaseField_module, ONLY: CheckMaterialParameters, umat => umatPF
  USE BMatrixJDet
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=2, NSHR=1, NDI=3
  CONTAINS
    INCLUDE "PFUEL_PFF_240103.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D-model, rectangle, quadratic shape functions, full integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PF_2D_Rect_Quadratic_Full
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE ShapeFunc2D_Square_Quadr, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr2D_Square9GP, ONLY: NGP, GPWeight, GPPos
  USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE PhaseField_module, ONLY: CheckMaterialParameters, umat => umatPF
  USE BMatrixJDet
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=2, NSHR=1, NDI=3
  CONTAINS
    INCLUDE "PFUEL_PFF_240103.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D-model, rectangle, quadratic shape functions, full integration with 16 points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PF_2D_Rect_Quadratic_Full16
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE ShapeFunc2D_Square_Quadr, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr2D_Square16GP, ONLY: NGP, GPWeight, GPPos
  USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE PhaseField_module, ONLY: CheckMaterialParameters, umat => umatPF
  USE BMatrixJDet
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=2, NSHR=1, NDI=3
  CONTAINS
    INCLUDE "PFUEL_PFF_240103.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D-model, rectangle, quadratic shape functions, full integration with 25 points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PF_2D_Rect_Quadratic_Full25
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE ShapeFunc2D_Square_Quadr, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr2D_Square25GP, ONLY: NGP, GPWeight, GPPos
  USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE PhaseField_module, ONLY: CheckMaterialParameters, umat => umatPF
  USE BMatrixJDet
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=2, NSHR=1, NDI=3
  CONTAINS
    INCLUDE "PFUEL_PFF_240103.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D-model, rectangle, quadratic shape functions, full integration, Newton-Cotes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PF_2D_Rect_Quadratic_Full_NC
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE ShapeFunc2D_Square_Quadr, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr2D_Square9GP_NC, ONLY: NGP, GPWeight, GPPos
  USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE PhaseField_module, ONLY: CheckMaterialParameters, umat => umatPF
  USE BMatrixJDet
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=2, NSHR=1, NDI=3
  CONTAINS
    INCLUDE "PFUEL_PFF_240103.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D-model, rectangle, cubic shape functions, reduced integration with 9 points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PF_2D_Rect_Cubic_Reduced
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE ShapeFunc2D_Square_Cubic, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr2D_Square9GP, ONLY: NGP, GPWeight, GPPos
  USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE PhaseField_module, ONLY: CheckMaterialParameters, umat => umatPF
  USE BMatrixJDet
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=2, NSHR=1, NDI=3
  CONTAINS
    INCLUDE "PFUEL_PFF_240103.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D-model, rectangle, cubic shape functions, full integration with 16 points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PF_2D_Rect_Cubic_Full16
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE ShapeFunc2D_Square_Cubic, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr2D_Square16GP, ONLY: NGP, GPWeight, GPPos
  USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE PhaseField_module, ONLY: CheckMaterialParameters, umat => umatPF
  USE BMatrixJDet
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=2, NSHR=1, NDI=3
  CONTAINS
    INCLUDE "PFUEL_PFF_240103.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D-model, rectangle, cubic shape functions, full integration with 25 points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PF_2D_Rect_Cubic_Full25
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE ShapeFunc2D_Square_Cubic, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr2D_Square25GP, ONLY: NGP, GPWeight, GPPos
  USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE PhaseField_module, ONLY: CheckMaterialParameters, umat => umatPF
  USE BMatrixJDet
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=2, NSHR=1, NDI=3
  CONTAINS
    INCLUDE "PFUEL_PFF_240103.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D-model, rectangle, cubic shape functions, full integration with 36 points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PF_2D_Rect_Cubic_Full36
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE ShapeFunc2D_Square_Cubic, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr2D_Square36GP, ONLY: NGP, GPWeight, GPPos
  USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE PhaseField_module, ONLY: CheckMaterialParameters, umat => umatPF
  USE BMatrixJDet
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=2, NSHR=1, NDI=3
  CONTAINS
    INCLUDE "PFUEL_PFF_240103.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D-model, rectangle, cubic shape functions, reduced integration with 4 points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PF_2D_Rect_Cubic_Reduced4
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE ShapeFunc2D_Square_Cubic, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr2D_Square4GP, ONLY: NGP, GPWeight, GPPos
  USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE PhaseField_module, ONLY: CheckMaterialParameters, umat => umatPF
  USE BMatrixJDet
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=2, NSHR=1, NDI=3
  CONTAINS
    INCLUDE "PFUEL_PFF_240103.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D-model, triangle, linear shape functions, full integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PF_2D_Tria_Linear_Full
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE ShapeFunc2D_Triangle_Lin, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr2D_Triangle4GP, ONLY: NGP, GPWeight, GPPos
  USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE PhaseField_module, ONLY: CheckMaterialParameters, umat => umatPF
  USE BMatrixJDet
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=2, NSHR=1, NDI=3
  CONTAINS
    INCLUDE "PFUEL_PFF_240103.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D-model, triangle, quadratic shape functions, full integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PF_2D_Tria_Quadratic_Full
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE ShapeFunc2D_Triangle_Quad, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr2D_Triangle4GP, ONLY: NGP, GPWeight, GPPos
  USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE PhaseField_module, ONLY: CheckMaterialParameters, umat => umatPF
  USE BMatrixJDet
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=2, NSHR=1, NDI=3
  CONTAINS
    INCLUDE "PFUEL_PFF_240103.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D-model, wedge, linear shape functions, reduced integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!MODULE PF_3D_Wedge_Linear_Reduced
!  USE ABQINTERFACE
!  USE FLOATNUMBERS
!  USE ABAModul, ONLY: XEXIT
!  USE ShapeFunc2D_Triangle_Lin, ONLY: ShapeFunc, ShapeFuncDeriv
!  USE Integr2D_Triangle1GP, ONLY: NGP, GPWeight, GPPos
!!  USE BMatrixJDet
!  IMPLICIT NONE
!  PRIVATE
!  PUBLIC :: UEL
!  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=3, NSHR=3, NDI=3
!  CONTAINS
!    INCLUDE "PFUEL_PFF_240103.f90"
!END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D-model, wedge, linear shape functions, full integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!MODULE PF_3D_Wedge_Linear_Full
!  USE ABQINTERFACE
!  USE FLOATNUMBERS
!  USE ABAModul, ONLY: XEXIT
!  USE ShapeFunc2D_Triangle_Lin, ONLY: ShapeFunc, ShapeFuncDeriv
!! Aenderung 07.07.2015, Stephan Roth
!!  USE Integr2D_Triangle3GP, ONLY: NGP, GPWeight, GPPos
!  USE Integr2D_Triangle3GPNodalPos, ONLY: NGP, GPWeight, GPPos
!!  USE BMatrixJDet
!  IMPLICIT NONE
!  PRIVATE
!  PUBLIC :: UEL
!  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=3, NSHR=3, NDI=3
!  CONTAINS
!    INCLUDE "PFUEL_PFF_240103.f90"
!END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D-model, wedge, quadratic shape functions, reduced integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!MODULE PF_3D_Wedge_Quadratic_Reduced
!  USE ABQINTERFACE
!  USE FLOATNUMBERS
!  USE ABAModul, ONLY: XEXIT
!  USE ShapeFunc2D_Triangle_Quad, ONLY: ShapeFunc, ShapeFuncDeriv
!  USE Integr2D_Triangle3GP, ONLY: NGP, GPWeight, GPPos
!!  USE BMatrixJDet
!  IMPLICIT NONE
!  PRIVATE
!  PUBLIC :: UEL
!  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=3, NSHR=3, NDI=3
!  CONTAINS
!    INCLUDE "PFUEL_PFF_240103.f90"
!END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D-model, wedge, quadratic shape functions, full integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!MODULE PF_3D_Wedge_Quadratic_Full
!  USE ABQINTERFACE
!  USE FLOATNUMBERS
!  USE ABAModul, ONLY: XEXIT
!  USE ShapeFunc2D_Triangle_Quad, ONLY: ShapeFunc, ShapeFuncDeriv
!  USE Integr2D_Triangle4GP, ONLY: NGP, GPWeight, GPPos
!!  USE BMatrixJDet
!  IMPLICIT NONE
!  PRIVATE
!  PUBLIC :: UEL
!  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=3, NSHR=3, NDI=3
!  CONTAINS
!    INCLUDE "PFUEL_PFF_240103.f90"
!END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D-model, hex, linear shape functions, reduced integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PF_3D_Hex_Linear_Reduced
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE ShapeFunc3D_Hex_Lin, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr3D_Hex1GP, ONLY: NGP, GPWeight, GPPos
  USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE PhaseField_module, ONLY: CheckMaterialParameters, umat => umatPF
  USE BMatrixJDet
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=3, NSHR=3, NDI=3
  CONTAINS
    INCLUDE "PFUEL_PFF_240103.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D-model, hex, linear shape functions, full integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PF_3D_Hex_Linear_Full
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE ShapeFunc3D_Hex_Lin, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr3D_Hex8GP, ONLY: NGP, GPWeight, GPPos
  USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE PhaseField_module, ONLY: CheckMaterialParameters, umat => umatPF
  USE BMatrixJDet
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=3, NSHR=3, NDI=3
  CONTAINS
    INCLUDE "PFUEL_PFF_240103.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D-model, hex, quadratic shape functions, reduced integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PF_3D_Hex_Quadratic_Reduced
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE ShapeFunc3D_Hex_Quadr, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr3D_Hex8GP, ONLY: NGP, GPWeight, GPPos
  USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE PhaseField_module, ONLY: CheckMaterialParameters, umat => umatPF
  USE BMatrixJDet
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=3, NSHR=3, NDI=3
  CONTAINS
    INCLUDE "PFUEL_PFF_240103.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D-model, hex, quadratic shape functions, full integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PF_3D_Hex_Quadratic_Full
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE ShapeFunc3D_Hex_Quadr, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr3D_Hex27GP, ONLY: NGP, GPWeight, GPPos
  USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE PhaseField_module, ONLY: CheckMaterialParameters, umat => umatPF
  USE BMatrixJDet
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=3, NSHR=3, NDI=3
  CONTAINS
    INCLUDE "PFUEL_PFF_240103.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! MAIN UEL ROUTINE WITH CHOICE OF ELEMENTROUTINE BY EVALUATING INPUTDATA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PFUELModule

  USE ABQINTERFACE
  USE FLOATNUMBERS

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  ! number of integer parameters
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: numIntPar=6

  CONTAINS

!------------------------------------------------------------------------------------

    SUBROUTINE UEL(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                   mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                   ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                   pnewdt,jprops,njprop,period)

      USE ABQINTERFACE
      USE FLOATNUMBERS
      USE ABAModul, ONLY: XEXIT
      USE PF_2D_Rect_Linear_Reduced, ONLY: UEL1 => UEL
      USE PF_2D_Rect_Linear_Full, ONLY: UEL2 => UEL
      USE PF_2D_Rect_Quadratic_Reduced, ONLY: UEL3 => UEL
      USE PF_2D_Rect_Quadratic_Full, ONLY: UEL4 => UEL
!      USE PF_3D_Wedge_Linear_Reduced, ONLY: UEL5 => UEL
!      USE PF_3D_Wedge_Linear_Full, ONLY: UEL6 => UEL
!      USE PF_3D_Wedge_Quadratic_Reduced, ONLY: UEL7 => UEL
!      USE PF_3D_Wedge_Quadratic_Full, ONLY: UEL8 => UEL
      USE PF_3D_Hex_Linear_Reduced, ONLY: UEL9 => UEL
      USE PF_3D_Hex_Linear_Full, ONLY: UEL10 => UEL
      USE PF_3D_Hex_Quadratic_Reduced, ONLY: UEL11 => UEL
      USE PF_3D_Hex_Quadratic_Full, ONLY: UEL12 => UEL
      USE PF_1D_Line_Linear_Reduced_Spherisym, ONLY: UEL13 => UEL
      USE PF_1D_Line_Quadratic_Reduced_Spherisym, ONLY: UEL14 => UEL
      USE PF_1D_Line_Quadratic_Spherisym, ONLY: UEL15 => UEL
      USE PF_1D_Line_Quadratic_Reduced, ONLY: UEL16 => UEL
      USE PF_1D_Line_LinCubic, ONLY: UEL17 => UEL
      USE PF_2D_Tria_Quadratic_Full, ONLY: UEL18 => UEL
      USE PF_2D_Rect_Quadratic_Full16, ONLY: UEL19 => UEL
      USE PF_2D_Rect_Cubic_Full16, ONLY: UEL20 => UEL
      USE PF_2D_Tria_Linear_Full, ONLY: UEL21 => UEL
      USE PF_2D_Rect_Quadratic_Full25, ONLY: UEL22 => UEL
      USE PF_2D_Rect_Cubic_Full25, ONLY: UEL23 => UEL
      USE PF_2D_Rect_Cubic_Reduced, ONLY: UEL24 => UEL
      USE PF_2D_Rect_Cubic_Full36, ONLY: UEL25 => UEL
      USE PF_2D_Rect_Cubic_Reduced4, ONLY: UEL26 => UEL
      USE PF_2D_Rect_Quadratic_Full_NC, ONLY: UEL27 => UEL

      IMPLICIT NONE

      INTEGER(kind=AbqIK), INTENT(IN) :: nprops, njprop, mcrd, nnode, ndofel, &
                                         npredf, mlvarx, nrhs, nsvars, jtype, &
                                         kstep, kinc, jelem, ndload, mdload, &
                                         lflags(*), jdltyp(mdload,*), &
                                         jprops(njprop)
      REAL(kind=AbqRK), INTENT(IN) :: props(nprops), &
                                      coords(mcrd,nnode), u(ndofel), period, &
                                      v(ndofel), a(ndofel), time(2), dtime, &
                                      predef(2,npredf,nnode), du(mlvarx,*), &
                                      adlmag(mdload,*), ddlmag(mdload,*), params(*)
      REAL(kind=AbqRK), INTENT(INOUT) :: energy(8), pnewdt, svars(nsvars)
      REAL(kind=AbqRK), INTENT(OUT) :: amatrx(ndofel,ndofel), rhs(mlvarx,*)
      INTEGER(kind=AbqIK) :: dimension
      LOGICAL :: isReduced, isSpherisymmetric, secondIntegrationScheme, thirdIntegrationScheme, forthIntegrationScheme, NewtonCotesIntegrationScheme

      ! initialisation
      isReduced=.FALSE.
      isSpherisymmetric=.FALSE.
      secondIntegrationScheme=.FALSE.
      thirdIntegrationScheme=.FALSE.
      forthIntegrationScheme=.FALSE.
      NewtonCotesIntegrationScheme=.FALSE.

      ! query: number of integer parameters
      IF ( time(2) .EQ. zero ) THEN
        IF (njprop .LT. numIntPar) THEN
          WRITE(7,*) 'Number of integer parameters', njprop, 'less than', numIntPar, '. Check input file. EXIT'
          CALL XEXIT()
        END IF
      END IF

!      ! integration: 1 - reduced, 0 - full (default),
      ! integration: 1 - reduced, 0 - full (default), 2 - full (2nd integration scheme), ...
      IF (jprops(1) .EQ. 1) THEN
        isReduced=.TRUE.
      ELSE IF (jprops(1) .EQ. 2) THEN
        secondIntegrationScheme=.TRUE.
      ELSE IF (jprops(1) .EQ. 3) THEN
        thirdIntegrationScheme=.TRUE.
      ELSE IF (jprops(1) .EQ. 4) THEN
        forthIntegrationScheme=.TRUE.
      ELSE IF (jprops(1) .EQ. 5) THEN
        NewtonCotesIntegrationScheme=.TRUE.
      END IF
      ! symmetry: 1 - axi-symmetric, 2 - spheri-symmetric, 0 - no symmetry (default)
      IF (jprops(2) .EQ. 2) isSpherisymmetric=.TRUE.
      ! dimension
      dimension = jprops(5)

      ! element type
      SELECT CASE(nnode)
      CASE(2)
      ! linear line, spheri-symmetric
        IF (isReduced) THEN
          CALL UEL13(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                     mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                     ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                     pnewdt,jprops,njprop,period)
        ELSE
          WRITE(7,*) 'linear spheri-symmetric element mit full integration is not defined. EXIT'
          CALL XEXIT()
        END IF
      CASE(3)
        IF (dimension .EQ. 1) THEN 
        ! quadratic line
          IF (isSpherisymmetric) THEN
          ! spheri-symmetric
            IF (isReduced) THEN
              CALL UEL14(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                         mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                         ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                         pnewdt,jprops,njprop,period)
            ELSE
              CALL UEL15(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                         mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                         ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                         pnewdt,jprops,njprop,period)
            END IF
          ELSE
            IF (isReduced) THEN
              CALL UEL16(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                         mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                         ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                         pnewdt,jprops,njprop,period)
            END IF
          END IF
        ELSE
        ! linear triangle
          CALL UEL21(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                     mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                     ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                     pnewdt,jprops,njprop,period)
        END IF
      CASE(4)
        IF (dimension .EQ. 1) THEN 
        ! 1D kubic line 
          CALL UEL17(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                     mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                     ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                     pnewdt,jprops,njprop,period)
        ELSE
        ! linear rectangle
          IF (isReduced) THEN
            CALL UEL1(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                      mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                      ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                      pnewdt,jprops,njprop,period)
          ELSE
            CALL UEL2(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                      mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                      ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                      pnewdt,jprops,njprop,period)
          END IF
        END IF
      CASE(6)
      ! quadratic triangle
        IF (isReduced) THEN
write(7,*) 'no quadratic triangle with reduced integration available'
        ELSE
          CALL UEL18(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                     mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                     ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                     pnewdt,jprops,njprop,period)
        END IF
!      CASE(6)
!!        SELECT CASE(mcrd)
!        SELECT CASE(dimension)
!        CASE(3)
!        ! linear wedge
!          IF (isReduced) THEN
!            CALL UEL5(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
!                      mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
!                      ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
!                      pnewdt,jprops,njprop,period)
!          ELSE
!            CALL UEL6(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
!                      mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
!                      ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
!                      pnewdt,jprops,njprop,period)
!          END IF
!        CASE(2)
!        ! quadratic rectangle
!          IF (isReduced) THEN
!            CALL UEL3(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
!                      mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
!                      ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
!                      pnewdt,jprops,njprop,period)
!          ELSE
!            CALL UEL4(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
!                      mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
!                      ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
!                      pnewdt,jprops,njprop,period)
!          END IF
!        END SELECT
!      CASE(12)
!      ! quadratic wedge
!        IF (isReduced) THEN
!          CALL UEL7(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
!                    mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
!                    ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
!                    pnewdt,jprops,njprop,period)
!        ELSE
!          CALL UEL8(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
!                    mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
!                    ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
!                    pnewdt,jprops,njprop,period)
!        END IF
      CASE(8)
!        SELECT CASE(mcrd)
        SELECT CASE(dimension)
        CASE(2)
        ! quadratic rectangle
          IF (isReduced) THEN
            CALL UEL3(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                      mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                      ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                      pnewdt,jprops,njprop,period)
          ELSE IF (secondIntegrationScheme) THEN
            CALL UEL19(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                       mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                       ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                       pnewdt,jprops,njprop,period)
          ELSE IF (thirdIntegrationScheme) THEN
            CALL UEL22(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                       mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                       ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                       pnewdt,jprops,njprop,period)
          ELSE IF (NewtonCotesIntegrationScheme) THEN
            CALL UEL27(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                       mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                       ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                       pnewdt,jprops,njprop,period)

          ELSE
            CALL UEL4(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                      mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                      ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                      pnewdt,jprops,njprop,period)
          END IF
        CASE(3)
        ! linear hex
          IF (isReduced) THEN
            CALL UEL9(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                      mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                      ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                      pnewdt,jprops,njprop,period)
          ELSE
            CALL UEL10(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                       mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                       ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                       pnewdt,jprops,njprop,period)
          END IF
        END SELECT
      CASE(12)
      ! cubic quadrilateral
        IF (isReduced) THEN
          ! 9 integration points
          CALL UEL24(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                     mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                     ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                     pnewdt,jprops,njprop,period)
        ELSE IF (secondIntegrationScheme) THEN
          ! 25 integration points
          CALL UEL23(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                     mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                     ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                     pnewdt,jprops,njprop,period)
        ELSE IF (thirdIntegrationScheme) THEN
          ! 36 integration points
          CALL UEL25(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                     mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                     ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                     pnewdt,jprops,njprop,period)
        ELSE IF (forthIntegrationScheme) THEN
          ! 4 integration points
          CALL UEL26(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                     mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                     ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                     pnewdt,jprops,njprop,period)
        ELSE
          ! 16 integration points
          CALL UEL20(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                     mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                     ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                     pnewdt,jprops,njprop,period)
        END IF
      CASE(20)
      ! quadratic hex
        IF (isReduced) THEN
          CALL UEL11(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                     mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                     ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                     pnewdt,jprops,njprop,period)
        ELSE
          CALL UEL12(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                     mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                     ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                     pnewdt,jprops,njprop,period)
        END IF
      CASE DEFAULT
        WRITE(7,*) 'Number of nodes does not refer to a PFUEL. EXIT', nnode
        CALL XEXIT()
      END SELECT

      RETURN

    END SUBROUTINE UEL

!------------------------------------------------------------------------------------

END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
