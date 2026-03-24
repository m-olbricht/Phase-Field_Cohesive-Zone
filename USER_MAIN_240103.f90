!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ABAQUS user defined subroutines
!
! Stephan Roth, TU Bergakademie Freiberg, 03.01.2024
!
! 03.01.2024: phase-field element and cohesive zone element
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                     U E L
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE uel(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
               mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
               ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
               pnewdt,jprops,njprop,period)

  USE ABQINTERFACE_PF
  USE FLOATNUMBERS
  USE ABAModul
  ! cohesive zone user element (Stephan Roth)
  USE CZUELModule, only: CZUEL => UEL
  ! phase field user element (Stephan Roth)
  USE PFUELModule, only: PFUEL => UEL

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

  ! dummies
  amatrx=zero; rhs(1:ndofel,1)=zero

  IF ( jtype .LE. 10 ) THEN
    ! phase field element
    CALL PFUEL(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
               mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
               ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
               pnewdt,jprops,njprop,period)
  ELSE IF ( jtype .LE. 20 ) THEN
    ! cohesive zone element
    CALL CZUEL(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
               mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
               ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
               pnewdt,jprops,njprop,period)
  ELSE
    WRITE(7,*) 'No UEL found. EXIT'
    CALL XEXIT()
  END IF

END SUBROUTINE uel

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                    D I S P
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE disp(u,kstep,kinc,time,node,noel,jdof,coords)

  USE ABQINTERFACE
  USE FLOATNUMBERS
  USE ABAModul

  IMPLICIT NONE
  INTEGER(kind=AbqIK), INTENT(IN) :: kstep, kinc, node, noel, jdof
  REAL(kind=AbqRK), INTENT(IN) :: time(2), coords(3)
  REAL(kind=AbqRK), INTENT(INOUT) :: u(3)
  REAL(kind=AbqRK) :: r

  r = zero
  CALL RANDOM_NUMBER(r)
  u = zero
  ! one order parameter
  u(1) = r

!  ! linear distribution
!  u(1) = half + half/hundred*coords(2)
!  u(1) = one/ten*coords(1)

END SUBROUTINE disp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

