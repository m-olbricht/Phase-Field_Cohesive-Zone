!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Module library for ABAQUS user defined elements: phase field elements
!
! Stephan Roth, TU Bergakademie Freiberg, 30.07.2020
!
! 30.07.2020: Multi-phase multi-component
! 15.06.2021: concentrations independent on phases
! 03.01.2024: without NCP, DOF: displacements, 1xdamage variable
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! B-Matrix AND JACOBI-DETERMINANT
! jprops(2) = 0 --> no symmetry
!           = 1 --> axi-symmetric element
!           = 2 --> spheri-symmetric element
!
! jprops(3) = 0 --> default, NLGEOM
!           = 1 --> initial configuration, nominal nodal positions
!           = 2 --> current configuration, updated nodal positions
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE BMatrixJDet

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::  BMatrixJacobian, getDeterminante, inverseMatrix
  CONTAINS

!------------------------------------------------------------------------------------

    SUBROUTINE BMatrixJacobian(coords,u,D,nnode,ndofel,ndi,nshr,ntens,njprop,shapeFunctions,shapeFuncDerivatives, &
                               jprops,prop_thickness,isNlgeom,Matrix_B,JacobiDet,drot,GPcoords,positivDetJ)
      ! JACOBI-Determinante und B-Matrix
      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE CZMath, ONLY: VectorNorm, crossProduct

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D, nnode, ndofel, ndi, nshr, ntens, njprop, jprops(njprop)
      LOGICAL, INTENT(IN) :: isNlgeom
      LOGICAL, INTENT(OUT) :: positivDetJ
      REAL(kind=AbqRK), INTENT(IN) :: coords(D,nnode), u(ndofel), shapeFunctions(nnode), shapeFuncDerivatives(D,nnode), &
                                      prop_thickness
      REAL(kind=AbqRK), INTENT(OUT) :: Matrix_B(ntens,ndofel), JacobiDet, drot(3,3), GPcoords(3)
      REAL(kind=AbqRK) :: Matrix_u(ndofel/nnode,nnode), currentNodeLocation(D,nnode), Matrix_J(D,D), Matrix_J_inv(D,D), transformedShapeFuncDerivatives(D,nnode), radius, TrafoMat(3,3)
      REAL(KIND=AbqRK), PARAMETER :: pi=acos(-one)
      INTEGER(kind=AbqIK) :: i1, i2, i3
      LOGICAL :: isAxisymmetric, isSpherisymmetric

      positivDetJ = .TRUE.; Matrix_B = zero; JacobiDet = zero; drot = zero; GPcoords = zero; isAxisymmetric = .FALSE.; isSpherisymmetric = .FALSE.; radius = zero

      ! Knoten-Koordinaten (aktuell oder nominal), Fallunterscheidung nach isNlgeom
      Matrix_u = zero
      ! jprops(3): 0 - nach NLGEOM
      !            1 - linear
      !            2 - nichtlinear
      IF ( ( (isNlgeom) .AND. (jprops(3) .NE. 1) ) .OR. ( (.NOT. isNlgeom) .AND. (jprops(3) .EQ. 2) ) ) THEN
        ! Verschiebungen in Matrixform
        Matrix_u = RESHAPE(u, (/ ndofel/nnode, nnode /))
      END IF
      currentNodeLocation = coords+Matrix_u(1:D,:)

      ! integration point coordinates (initial configuration)
      FORALL (i1=1:D)
        GPcoords(i1) = dot_product(shapeFunctions,coords(i1,:))
      END FORALL

      ! axi- or spheri-symmetry
      IF ( (D .EQ. 2) .AND. (jprops(2) .EQ. 1) ) THEN
        isAxisymmetric = .TRUE.
      ELSE IF ( (D .EQ. 1) .AND. (jprops(2) .EQ. 2) ) THEN
        isSpherisymmetric = .TRUE.
      END IF
      IF ( (isAxisymmetric) .OR. (isSpherisymmetric) ) THEN
        ! initial configuration
        !radius = GPcoords(1)
        ! deformed configuration
        radius = dot_product(shapeFunctions,currentNodeLocation(1,:))
      END IF

      ! Jacobi-Matrix
      Matrix_J = matmul(shapeFuncDerivatives,transpose(currentNodeLocation))

      ! JACOBI-Determinante
      CALL getDeterminante(Matrix_J,D,JacobiDet,positivDetJ)
!write(7,*) 'JacobiDet', JacobiDet
!write(7,"(4e15.5)") Matrix_J

      IF (positivDetJ) THEN
        ! inverse Jacobi-Matrix
        Matrix_J_inv = inverseMatrix(Matrix_J,D,JacobiDet)
!write(7,"(4e15.5)") Matrix_J
!write(7,"(4e15.5)") Matrix_J_inv
!write(7,"(4e15.5)") matmul(Matrix_J_inv,Matrix_J)

        ! transformierte Ableitungen der Formfunktionen
        transformedShapeFuncDerivatives = matmul(Matrix_J_inv,shapeFuncDerivatives)

!write(7,*) 'transformedShapeFuncDerivatives'
!do i1=1,nnode
!  write(7,"(2e15.5)")transformedShapeFuncDerivatives(:,i1)
!end do

        ! Matrix_B (fuer Gleitungen anstatt Schubverzerrungen)
        Matrix_B = zero
        SELECT CASE(D)
        CASE(1)
          IF (isSpherisymmetric) THEN
            ! spherisymmetric
            FORALL (i1=1:nnode)
              Matrix_B( 1, (i1-1)*ndofel/nnode+1 ) =  transformedShapeFuncDerivatives(1,i1)      ! u_r,r
              Matrix_B( 2, (i1-1)*ndofel/nnode+1 ) =  shapeFunctions(i1)/radius                  ! u_z,z
              Matrix_B( 3, (i1-1)*ndofel/nnode+1 ) =  shapeFunctions(i1)/radius                  ! u_phi,phi
            END FORALL
          ELSE
            FORALL (i1=1:nnode)
              Matrix_B( 1, (i1-1)*ndofel/nnode+1 ) =  transformedShapeFuncDerivatives(1,i1)      ! u_x,x
            END FORALL
          END IF
        CASE(2)
          IF (isAxisymmetric) THEN
            ! axisymmetric
            FORALL (i1=1:nnode)
              Matrix_B( 1, (i1-1)*ndofel/nnode+1 ) =  transformedShapeFuncDerivatives(1,i1)      ! u_r,r
              Matrix_B( 2, (i1-1)*ndofel/nnode+2 ) =  transformedShapeFuncDerivatives(2,i1)      ! u_z,z
              Matrix_B( 3, (i1-1)*ndofel/nnode+1 ) =  shapeFunctions(i1)/radius                  ! u_phi,phi
              Matrix_B( 4, (i1-1)*ndofel/nnode+1 ) =  transformedShapeFuncDerivatives(2,i1)      ! u_r,z
              Matrix_B( 4, (i1-1)*ndofel/nnode+2 ) =  transformedShapeFuncDerivatives(1,i1)      ! u_z,r
            END FORALL
          ELSE
            ! plane strain
            FORALL (i1=1:nnode)
              Matrix_B( 1, (i1-1)*ndofel/nnode+1 ) =  transformedShapeFuncDerivatives(1,i1)      ! u_x,x
              Matrix_B( 2, (i1-1)*ndofel/nnode+2 ) =  transformedShapeFuncDerivatives(2,i1)      ! u_y,y
              ! third row empty: eps_33
              Matrix_B( 4, (i1-1)*ndofel/nnode+1 ) =  transformedShapeFuncDerivatives(2,i1)      ! u_x,y
              Matrix_B( 4, (i1-1)*ndofel/nnode+2 ) =  transformedShapeFuncDerivatives(1,i1)      ! u_y,x
            END FORALL
          END IF
        CASE(3)
          FORALL (i1=1:nnode)
            Matrix_B( 1, (i1-1)*ndofel/nnode+1 ) =  transformedShapeFuncDerivatives(1,i1)      ! u_x,x
            Matrix_B( 2, (i1-1)*ndofel/nnode+2 ) =  transformedShapeFuncDerivatives(2,i1)      ! u_y,y
            Matrix_B( 3, (i1-1)*ndofel/nnode+3 ) =  transformedShapeFuncDerivatives(3,i1)      ! u_z,z
            Matrix_B( 4, (i1-1)*ndofel/nnode+1 ) =  transformedShapeFuncDerivatives(2,i1)      ! u_x,y
            Matrix_B( 4, (i1-1)*ndofel/nnode+2 ) =  transformedShapeFuncDerivatives(1,i1)      ! u_y,x
            Matrix_B( 5, (i1-1)*ndofel/nnode+1 ) =  transformedShapeFuncDerivatives(3,i1)      ! u_x,z
            Matrix_B( 5, (i1-1)*ndofel/nnode+3 ) =  transformedShapeFuncDerivatives(1,i1)      ! u_z,x
            Matrix_B( 6, (i1-1)*ndofel/nnode+2 ) =  transformedShapeFuncDerivatives(3,i1)      ! u_y,z
            Matrix_B( 6, (i1-1)*ndofel/nnode+3 ) =  transformedShapeFuncDerivatives(2,i1)      ! u_z,y
          END FORALL
        END SELECT
        !
        ! 1 Ordnungsparameter Schaedigungsvariable (Schleife von 1 bis 1)
        DO i2=1,1
          DO i1=1,nnode
            Matrix_B(ndi+nshr+(i2-1)*(D+1)+1, (i1-1)*ndofel/nnode+D+i2 ) = shapeFunctions(i1)
            ! Ableitung nach allen D-Richtungen
            FORALL (i3=1:D)
              Matrix_B(ndi+nshr+(i2-1)*(D+1)+1+i3, (i1-1)*ndofel/nnode+D+i2 ) = transformedShapeFuncDerivatives(i3,i1)
            END FORALL
          END DO
        END DO

        ! Transformationsmatrix
        drot = zero
        SELECT CASE(D)
        CASE(1)
          FORALL (i1=1:3) drot(i1,i1) = one
        CASE(2)
          drot(1:D,1:D) = Matrix_J
          drot(3,3)     = one
        CASE(3)
          drot = Matrix_J
        END SELECT

        ! Halbierung bei 2D-Dreieckselementen
        IF ( (D .EQ. 2) .AND. ( (nnode .EQ. 3) .OR. (nnode .EQ. 6) ) ) THEN
          JacobiDet = JacobiDet/two
        END IF
!        ! Halbierung bei Pentaedern
!        IF ( (D .EQ. 3) .AND. ( (nnode .EQ. 6) .OR. (nnode .EQ. 12) ) ) THEN
!          JacobiDet = JacobiDet/two
!        END IF

        ! enhance JACOBI determinante for 1D and 2D
        IF (isSpherisymmetric) THEN
        ! 1D, spheri-symmetric
          ! multiply JacobiDet by 4*pi*r**2
          JacobiDet = JacobiDet*four*pi*radius**2
        ELSE IF (D .EQ. 2) THEN
        ! 2D
          IF (isAxisymmetric) THEN
          ! axi-symmetric
            ! multiply JacobiDet by 2*pi*r
            JacobiDet = JacobiDet*two*pi*radius
          ELSE
          ! plane strain
            ! multiply JacobiDet by thickness>0
            IF (prop_thickness .GT. zero) JacobiDet = JacobiDet*prop_thickness
          END IF
        END IF

      END IF ! positivDetJ

    END SUBROUTINE BMatrixJacobian

!------------------------------------------------------------------------------------

    PURE SUBROUTINE getDeterminante(M,d,det,positivDet)
    ! compute determinante of (n,n)-Matrix M using Sarrus-rule

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: d
      REAL(kind=AbqRK), INTENT(IN) :: M(d,d)
      REAL(kind=AbqRK), INTENT(OUT) :: det
      LOGICAL, INTENT(OUT) :: positivDet

      det = zero
      SELECT CASE(d)
      CASE(1)
        det = M(1,1)
      CASE(2)
        det = M(1,1)*M(2,2) - M(1,2)*M(2,1)
      CASE(3)
        det = M(1,1)*M(2,2)*M(3,3) + M(2,1)*M(3,2)*M(1,3) + M(3,1)*M(1,2)*M(2,3) - M(1,3)*M(2,2)*M(3,1) - M(2,3)*M(3,2)*M(1,1) - M(3,3)*M(1,2)*M(2,1)
      END SELECT

      IF (det .LE. zero) THEN
        positivDet = .FALSE.
      ELSE
        positivDet = .TRUE.
      END IF

    END SUBROUTINE getDeterminante

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION inverseMatrix(M,d,det)
    ! compute inverse of (n,n)-Matrix M with positive determinante det

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: d
      REAL(kind=AbqRK), INTENT(IN) :: M(d,d), det
      DIMENSION inverseMatrix(d,d)

      inverseMatrix = zero
      SELECT CASE(D)
      CASE(1)
        ! Inverse
        IF ( M(1,1) .NE. zero ) THEN
          inverseMatrix(1,1)=one/M(1,1)
        ELSE
          inverseMatrix(1,1)=zero
        END IF
      CASE(2)
        ! Kofaktormatrix
        inverseMatrix(1,1) = +M(2,2)
        inverseMatrix(1,2) = -M(2,1)
        inverseMatrix(2,1) = -M(1,2)
        inverseMatrix(2,2) = +M(1,1)
        ! Inverse
        inverseMatrix = transpose(inverseMatrix)/det
      CASE(3)
        ! Kofaktormatrix
        inverseMatrix(1,1) = +(M(2,2)*M(3,3)-M(2,3)*M(3,2))
        inverseMatrix(1,2) = -(M(2,1)*M(3,3)-M(2,3)*M(3,1))
        inverseMatrix(1,3) = +(M(2,1)*M(3,2)-M(2,2)*M(3,1))
        inverseMatrix(2,1) = -(M(1,2)*M(3,3)-M(1,3)*M(3,2))
        inverseMatrix(2,2) = +(M(1,1)*M(3,3)-M(1,3)*M(3,1))
        inverseMatrix(2,3) = -(M(1,1)*M(3,2)-M(1,2)*M(3,1))
        inverseMatrix(3,1) = +(M(1,2)*M(2,3)-M(1,3)*M(2,2))
        inverseMatrix(3,2) = -(M(1,1)*M(2,3)-M(1,3)*M(2,1))
        inverseMatrix(3,3) = +(M(1,1)*M(2,2)-M(1,2)*M(2,1))
        ! Inverse
        inverseMatrix = transpose(inverseMatrix)/det
      END SELECT

    END FUNCTION inverseMatrix

!------------------------------------------------------------------------------------

END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

