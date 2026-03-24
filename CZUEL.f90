!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! ABAQUS user defined cohesive element subroutine
!
! Stephan Roth, TU Bergakademie Freiberg, 23.10.2014
!
! 23.10.2014: Anpassung an UELlib
! 13.08.2015: Aenderung Roth - Rotationsmatrix
! 03.01.2024: Aenderung Roth - Erweiterung um Schaedigungsvariable als Knotenfreiwert (Paggi, Reinoso)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    SUBROUTINE UEL(rhs,amatrx,svars,energy,ndofel,nrhs,nsvars,props,nprops,coords, &
                   mcrd,nnode,u,du,v,a,jtype,time,dtime,kstep,kinc,jelem,params, &
                   ndload,jdltyp,adlmag,predef,npredf,lflags,mlvarx,ddlmag,mdload, &
                   pnewdt,jprops,njprop,period)

      USE ABQINTERFACE

      IMPLICIT NONE

      ! UEL variables
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

      ! UMAT variables
      INTEGER(kind=AbqIK) :: ntens, layer, kspt, npt
      REAL(kind=AbqRK) :: predef_umat(npredf), dpred(npredf), rpl, drpldt, celent, F0(3,3), F1(3,3), GPcoords(3), drot(3,3)
      REAL(kind=AbqRK), ALLOCATABLE :: ddsddt(:), drplde(:)
      CHARACTER*45 :: cmname

      ! further variables
      REAL(kind=AbqRK) :: fint(ndofel), energy_gp(8), JacobiDet, prop_thickness
      REAL(kind=AbqRK), ALLOCATABLE :: Matrix_B(:,:), stran(:), dstran(:), Ct(:,:), stress(:)
      LOGICAL :: isNlgeom, positivDetJ, k_RHS, k_K, k_M

      ! initialisation
      positivDetJ=.TRUE.; isNlgeom=.FALSE.; k_RHS=.TRUE.; k_K=.TRUE.; k_M=.FALSE.
      ! umat dummies
      layer = 1; kspt = 1
      celent = one; GPcoords = zero; F0 = zero; F1 = zero

      ! 'Wer bin ich, was kann ich ?'

      ! evaluation of lflags
      ! geometric nonlinear
      IF (lflags(2) .EQ. 1) isNlgeom = .TRUE.

      SELECT CASE(lflags(3))
      CASE(1)
        k_RHS = .TRUE.
        k_K   = .TRUE.
        k_M   = .FALSE.
      CASE(2)
        k_RHS = .FALSE.
        k_K   = .TRUE.
        k_M   = .FALSE.
      CASE(3)
        k_RHS = .FALSE.
        k_K   = .FALSE.
        k_M   = .FALSE.
      CASE(4)
        k_RHS = .FALSE.
        k_K   = .FALSE.
        k_M   = .TRUE.
      CASE(5)
        k_RHS = .TRUE.
        k_K   = .FALSE.
        k_M   = .FALSE.
      CASE(6)
        k_RHS = .FALSE.
        k_K   = .FALSE.
        k_M   = .TRUE.
      CASE DEFAULT
        k_RHS = .FALSE.
        k_K   = .FALSE.
        k_M   = .FALSE.
        write(7,*) 'no output required'
      END SELECT

      ! number of tensor coordinates
      ntens=NDI+NSHR
      ntens=ntens+1+D ! damage variable (average) and damage variable jump and gradient of damage variable (of dimension D-1)

      ! allocation of UMAT-matrices
      ALLOCATE(ddsddt(ntens), drplde(ntens), Matrix_B(ntens,ndofel), Ct(ntens,ntens), stran(ntens), dstran(ntens), stress(ntens))
      ddsddt=zero; drplde=zero; Matrix_B = zero; stran = zero; dstran = zero; Ct = zero; stress = zero

      ! query at time 0: number of internal state variables per integration point and material parameters: answer from UMAT
      IF ( time(2) .EQ. zero ) THEN
        IF (nsvars .LT. NGP*numSDV) THEN
          ! nsvars=size(svars)
          WRITE(7,*) "A bigger number of solution-dependent state variables is required.", NGP*numSDV, nsvars
          WRITE(7,*) "Please note that this is a total number of SDV's, not per integration point."
          CALL XEXIT()
        END IF
        IF (nprops .LT. numMatPar) THEN
          ! nprops=size(props)
          WRITE(7,*) "A bigger number of real parameters is required.", numMatPar, nprops
          CALL XEXIT()
        END IF
        ! check of material parameters: answer from UMAT
        CALL CheckMaterialParameters(props)
      END IF

      ! thickness of 2D-Elements
      prop_thickness = one
      IF (D .EQ. 2) prop_thickness = props(thicknessIndex)

      ! Ende 'Wer bin ich, was kann ich ?'

      ! UMAT call only when neccessary
      IF (k_RHS .OR. k_K) THEN
        ! loop over integration points
        energy = zero
        amatrx = zero
        fint   = zero
        DO npt=1,NGP
          ! compute JACOBI determinante and B-matrix
          CALL BMatrixJacobian(coords,u,D,nnode,ndofel,ntens,njprop,ShapeFunc(GPPos(npt,1:D-1)),ShapeFuncDeriv(GPPos(npt,1:D-1)), &
                               NodeMapMatrix,jprops,prop_thickness,isNlgeom,Matrix_B,JacobiDet,drot,positivDetJ)
          IF (positivDetJ) THEN
          ! element OK
            ! compute local (transformed) separation, averaged damage, and damage junp
            stran  = matmul(Matrix_B, u)
            dstran = matmul(Matrix_B, du(1:ndofel,1))

            ! here extract stress from svars and compute GPcoords,F0,F1...
            stress = zero
            energy_gp = zero

            ! evaluation separation law
            CALL umat(stress,svars(numSDV*(npt-1)+1:numSDV*npt),Ct,energy_gp(2), &
                      energy_gp(4),energy_gp(3),rpl,ddsddt,drplde,drpldt,stran, &
                      dstran,time,dtime,predef_umat(1),dpred(1),predef_umat,dpred,cmname,NDI,NSHR,ntens, &
                      numSDV,props,nprops,GPcoords,drot,pnewdt,celent,F0,F1,jelem,npt,layer,kspt,kstep,kinc)

            ! electrostatic energy = penalty energy (SDV 9)
            energy_gp(1) = 0. !
            energy_gp(2) = svars(numSDV*(npt-1)+6) ! (SDV 6) -> stored CZ Energy
            energy_gp(3) = svars(numSDV*(npt-1)+7) ! (SDV 7) -> crack CZ Energy
            energy_gp(4) = svars(numSDV*(npt-1)+8) ! (SDV 8) -> contact Energy
            energy_gp(5) = svars(numSDV*(npt-1)+9) ! (SDV 9) -> penalty Energy (Damage jump)
            energy_gp(6) = svars(numSDV*(npt-1)+5) ! (SDV 5) -> total Energy CZ
            energy_gp(7) = svars(numSDV*(npt-1)+7) ! (SDV 7) -> crack CZ Energy

            ! integration over element
            ! energies
            !
            ! 1 = ALLKE -> reine Bruchenergie aus Bulk Bereich (s. PFUEL_PFF)
            ! 2 = ALLSE -> PFFCZ kombinierte linear elast. stored Energy (gesamt elast. gesp. Energie)
            ! 3 = ALLCD -> PFFCZ kombinierte gesp. Bruchenergie (gesamt Bruchenergie)
            ! 4 = ALLPD -> Kontaktenergie CZ bei negativer Normalseparation + viskose Regularisierung aus PFF (keine physikalische Interpretation. Nur um zu checken ob Energieterme auftreten)
            ! 5 = ALLVD -> Strafenergie für Damage Jump 
            ! 6 = ALLAE -> totale Energie CZ + totale Energie PFF (gesamte Energie im System)
            ! 7 = ALLEE -> reine Bruchenergie aus CZ Bereich (s. CZUEL.f90)
            !
            energy(:) = energy(:) + GPWeight(npt)*JacobiDet * energy_gp(:)
            ! cohesive stiffness matrix
            amatrx = amatrx + GPWeight(npt)*JacobiDet * matmul(transpose(Matrix_B),matmul(Ct,Matrix_B))
            ! cohesive fint
            fint = fint - GPWeight(npt)*JacobiDet * matmul(transpose(Matrix_B),stress)
          ELSE
          ! distorted element
            ! stop loop
            EXIT
          END IF
        END DO

        ! zero matrices for distorted elements
        IF (.NOT. positivDetJ) THEN
          WRITE(7,*) 'non-positive JACOBIAN det(J)<=0, distorted Element: ', jelem, ' matrices set to zero'
          amatrx = zero; fint = zero; energy = zero
        END IF

        rhs(1:ndofel,1) = fint
      END IF

      IF (k_M) THEN
      ! mass matrix
        amatrx = zero
      ELSE IF (.NOT. k_K) THEN
      ! no matrix required
        amatrx = zero
      END IF

      IF (.NOT. k_RHS) THEN
      ! no rhs required
        rhs(1:ndofel,1) = zero
      END IF

      DEALLOCATE(Matrix_B, ddsddt, drplde, stran, dstran, Ct, stress)

      RETURN

    END SUBROUTINE UEL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

