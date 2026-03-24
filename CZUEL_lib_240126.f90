!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ABAQUS user element library for cohesive elements
!
! Stephan Roth, TU Bergakademie Freiberg, 23.10.2014
!
! 23.10.2014: Anpassung an UELlib
! 07.07.2015: lineare Pentaeder: GPs auf Knotenpositionen (siehe ABAQUS COH3D6)
! 06.01.2016: Modifikationen fuer UELlib-Kompatibilitaet
! 03.01.2024: Modifikationen fuer PF-Kompatibilitaet
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D-model, rectangle, linear shape functions, reduced integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CZUEL_2D_Rect_Linear_Reduced
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE NodeMap_Rectangle_Linear, ONLY: NodeMapMatrix
  USE ShapeFunc1D_Line_Lin, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr1D_Line1GP, ONLY: NGP, GPWeight, GPPos
  USE CZ_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE CZUMAT_module, ONLY: umat => umatCZM, CheckMaterialParameters
  USE CohesiveBMatrixJDet, ONLY: BMatrixJacobian
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=2, NSHR=1, NDI=1
  CONTAINS
    INCLUDE "CZUEL.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D-model, rectangle, linear shape functions, full integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CZUEL_2D_Rect_Linear_Full
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE NodeMap_Rectangle_Linear, ONLY: NodeMapMatrix
  USE ShapeFunc1D_Line_Lin, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr1D_Line2GP, ONLY: NGP, GPWeight, GPPos
  USE CZ_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE CZUMAT_module, ONLY: umat => umatCZM, CheckMaterialParameters
  USE CohesiveBMatrixJDet, ONLY: BMatrixJacobian
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=2, NSHR=1, NDI=1
  CONTAINS
    INCLUDE "CZUEL.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D-model, rectangle, quadratic shape functions, reduced integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CZUEL_2D_Rect_Quadratic_Reduced
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE NodeMap_Rectangle_Quadratic, ONLY: NodeMapMatrix
  USE ShapeFunc1D_Line_Quad, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr1D_Line2GP, ONLY: NGP, GPWeight, GPPos
  USE CZ_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE CZUMAT_module, ONLY: umat => umatCZM, CheckMaterialParameters
  USE CohesiveBMatrixJDet, ONLY: BMatrixJacobian
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=2, NSHR=1, NDI=1
  CONTAINS
    INCLUDE "CZUEL.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D-model, rectangle, quadratic shape functions, full integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CZUEL_2D_Rect_Quadratic_Full
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE NodeMap_Rectangle_Quadratic, ONLY: NodeMapMatrix
  USE ShapeFunc1D_Line_Quad, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr1D_Line3GP, ONLY: NGP, GPWeight, GPPos
  USE CZ_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE CZUMAT_module, ONLY: umat => umatCZM, CheckMaterialParameters
  USE CohesiveBMatrixJDet, ONLY: BMatrixJacobian
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=2, NSHR=1, NDI=1
  CONTAINS
    INCLUDE "CZUEL.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D-model, wedge, linear shape functions, reduced integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CZUEL_3D_Wedge_Linear_Reduced
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE NodeMap_Wedge_Linear, ONLY: NodeMapMatrix
  USE ShapeFunc2D_Triangle_Lin, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr2D_Triangle1GP, ONLY: NGP, GPWeight, GPPos
  USE CZ_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE CZUMAT_module, ONLY: umat => umatCZM, CheckMaterialParameters
  USE CohesiveBMatrixJDet, ONLY: BMatrixJacobian
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=3, NSHR=2, NDI=1
  CONTAINS
    INCLUDE "CZUEL.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D-model, wedge, linear shape functions, full integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CZUEL_3D_Wedge_Linear_Full
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE NodeMap_Wedge_Linear, ONLY: NodeMapMatrix
  USE ShapeFunc2D_Triangle_Lin, ONLY: ShapeFunc, ShapeFuncDeriv
! Aenderung 07.07.2015, Stephan Roth
!  USE Integr2D_Triangle3GP, ONLY: NGP, GPWeight, GPPos
  USE Integr2D_Triangle3GPNodalPos, ONLY: NGP, GPWeight, GPPos
  USE CZ_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE CZUMAT_module, ONLY: umat => umatCZM, CheckMaterialParameters
  USE CohesiveBMatrixJDet, ONLY: BMatrixJacobian
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=3, NSHR=2, NDI=1
  CONTAINS
    INCLUDE "CZUEL.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D-model, wedge, quadratic shape functions, reduced integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CZUEL_3D_Wedge_Quadratic_Reduced
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE NodeMap_Wedge_Quadratic, ONLY: NodeMapMatrix
  USE ShapeFunc2D_Triangle_Quad, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr2D_Triangle3GP, ONLY: NGP, GPWeight, GPPos
  USE CZ_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE CZUMAT_module, ONLY: umat => umatCZM, CheckMaterialParameters
  USE CohesiveBMatrixJDet, ONLY: BMatrixJacobian
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=3, NSHR=2, NDI=1
  CONTAINS
    INCLUDE "CZUEL.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D-model, wedge, quadratic shape functions, full integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CZUEL_3D_Wedge_Quadratic_Full
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE NodeMap_Wedge_Quadratic, ONLY: NodeMapMatrix
  USE ShapeFunc2D_Triangle_Quad, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr2D_Triangle4GP, ONLY: NGP, GPWeight, GPPos
  USE CZ_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE CZUMAT_module, ONLY: umat => umatCZM, CheckMaterialParameters
  USE CohesiveBMatrixJDet, ONLY: BMatrixJacobian
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=3, NSHR=2, NDI=1
  CONTAINS
    INCLUDE "CZUEL.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D-model, hex, linear shape functions, reduced integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CZUEL_3D_Hex_Linear_Reduced
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE NodeMap_Hex_Linear, ONLY: NodeMapMatrix
  USE ShapeFunc2D_Square_Lin, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr2D_Square1GP, ONLY: NGP, GPWeight, GPPos
  USE CZ_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE CZUMAT_module, ONLY: umat => umatCZM, CheckMaterialParameters
  USE CohesiveBMatrixJDet, ONLY: BMatrixJacobian
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=3, NSHR=2, NDI=1
  CONTAINS
    INCLUDE "CZUEL.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D-model, hex, linear shape functions, full integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CZUEL_3D_Hex_Linear_Full
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE NodeMap_Hex_Linear, ONLY: NodeMapMatrix
  USE ShapeFunc2D_Square_Lin, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr2D_Square4GP, ONLY: NGP, GPWeight, GPPos
  USE CZ_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE CZUMAT_module, ONLY: umat => umatCZM, CheckMaterialParameters
  USE CohesiveBMatrixJDet, ONLY: BMatrixJacobian
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=3, NSHR=2, NDI=1
  CONTAINS
    INCLUDE "CZUEL.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D-model, hex, quadratic shape functions, reduced integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CZUEL_3D_Hex_Quadratic_Reduced
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE NodeMap_Hex_Quadratic, ONLY: NodeMapMatrix
  USE ShapeFunc2D_Square_Quadr, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr2D_Square4GP, ONLY: NGP, GPWeight, GPPos
  USE CZ_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE CZUMAT_module, ONLY: umat => umatCZM, CheckMaterialParameters
  USE CohesiveBMatrixJDet, ONLY: BMatrixJacobian
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=3, NSHR=2, NDI=1
  CONTAINS
    INCLUDE "CZUEL.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 3D-model, hex, quadratic shape functions, full integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CZUEL_3D_Hex_Quadratic_Full
  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul, ONLY: XEXIT
  USE NodeMap_Hex_Quadratic, ONLY: NodeMapMatrix
  USE ShapeFunc2D_Square_Quadr, ONLY: ShapeFunc, ShapeFuncDeriv
  USE Integr2D_Square9GP, ONLY: NGP, GPWeight, GPPos
  USE CZ_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
  USE CZUMAT_module, ONLY: umat => umatCZM, CheckMaterialParameters
  USE CohesiveBMatrixJDet, ONLY: BMatrixJacobian
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: UEL
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: D=3, NSHR=2, NDI=1
  CONTAINS
    INCLUDE "CZUEL.f90"
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! MAIN UEL ROUTINE WITH CHOICE OF ELEMENTROUTINE BY EVALUATING INPUTDATA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CZUELModule

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
      USE CZUEL_2D_Rect_Linear_Reduced, ONLY: UEL1 => UEL
      USE CZUEL_2D_Rect_Linear_Full, ONLY: UEL2 => UEL
      USE CZUEL_2D_Rect_Quadratic_Reduced, ONLY: UEL3 => UEL
      USE CZUEL_2D_Rect_Quadratic_Full, ONLY: UEL4 => UEL
      USE CZUEL_3D_Wedge_Linear_Reduced, ONLY: UEL5 => UEL
      USE CZUEL_3D_Wedge_Linear_Full, ONLY: UEL6 => UEL
      USE CZUEL_3D_Wedge_Quadratic_Reduced, ONLY: UEL7 => UEL
      USE CZUEL_3D_Wedge_Quadratic_Full, ONLY: UEL8 => UEL
      USE CZUEL_3D_Hex_Linear_Reduced, ONLY: UEL9 => UEL
      USE CZUEL_3D_Hex_Linear_Full, ONLY: UEL10 => UEL
      USE CZUEL_3D_Hex_Quadratic_Reduced, ONLY: UEL11 => UEL
      USE CZUEL_3D_Hex_Quadratic_Full, ONLY: UEL12 => UEL

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
      LOGICAL :: isReduced

      ! initialisation
      isReduced=.FALSE.

      ! query: number of integer parameters
      IF ( time(2) .EQ. zero ) THEN
        IF (njprop .LT. numIntPar) THEN
          WRITE(7,*) 'Number of integer parameters', njprop, 'less than', numIntPar, '. Check input file. EXIT'
          CALL XEXIT()
        END IF
      END IF

      ! integration: 1 - reduced, 0 - full (default),
      IF (jprops(1) .NE. 0) isReduced=.TRUE.

      ! element type
      SELECT CASE(nnode)
      CASE(4)
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
      CASE(6)
        SELECT CASE(mcrd)
        CASE(3)
        ! linear wedge
          IF (isReduced) THEN
            CALL UEL5(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                      mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                      ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                      pnewdt,jprops,njprop,period)
          ELSE
            CALL UEL6(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                      mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                      ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                      pnewdt,jprops,njprop,period)
          END IF
        CASE(2)
        ! quadratic rectangle
          IF (isReduced) THEN
            CALL UEL3(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                      mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                      ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                      pnewdt,jprops,njprop,period)
          ELSE
            CALL UEL4(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                      mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                      ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                      pnewdt,jprops,njprop,period)
          END IF
        END SELECT
      CASE(12)
      ! quadratic wedge
        IF (isReduced) THEN
          CALL UEL7(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                    mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                    ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                    pnewdt,jprops,njprop,period)
        ELSE
          CALL UEL8(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                    mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                    ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                    pnewdt,jprops,njprop,period)
        END IF
      CASE(8)
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
      CASE(16)
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
        WRITE(7,*) 'Number of nodes does not refer to a CZUEL. EXIT', nnode
        CALL XEXIT()
      END SELECT

      RETURN

    END SUBROUTINE UEL

!------------------------------------------------------------------------------------

END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
