!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Module library for ABAQUS user defined elements: cohesive elements
!
! Stephan Roth, TU Bergakademie Freiberg, 23.10.2014
!
! 23.10.2014: Anpassung an UELlib
! 06.07.2015: Aenderung Roth - Transformationsmatrix
! 07.08.2015: Aenderung Roth - CZE mit Symmetriebedingungung
! 13.08.2015: Aenderung Roth - Rotationsmatrix
! 06.01.2016: Modifikationen fuer UELlib-Kompatibilitaet
! 03.01.2024: Aenderung Roth - Erweiterung um Schaedigungsvariable als Knotenfreiwert (Paggi, Reinoso)
! 23.05.2025: Aenderung Roth - Rotationsmatrix fuer 1D/2D-Faelle angepasst 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! B-Matrix AND JACOBI-DETERMINANTE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CohesiveBMatrixJDet

  IMPLICIT NONE
  PRIVATE
  PUBLIC ::  BMatrixJacobian
  CONTAINS

!------------------------------------------------------------------------------------

    PURE SUBROUTINE BMatrixJacobian(coords,u,D,nnode,ndofel,ntens,njprop,shapeFunctions,shapeFuncDerivatives, &
                                    NodeMapMatrix,jprops,prop_thickness,isNlgeom,Matrix_B,JacobiDet,drot,positivDetJ)
      ! JACOBI-Determinante und B-Matrix
      USE ABQINTERFACE
      USE FLOATNUMBERS
      USE CZMath, ONLY: VectorNorm, crossProduct

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D, nnode, ndofel, ntens, njprop, jprops(njprop), NodeMapMatrix(nnode/2,2)
      LOGICAL, INTENT(IN) :: isNlgeom
      LOGICAL, INTENT(OUT) :: positivDetJ
      REAL(kind=AbqRK), INTENT(IN) :: coords(D,nnode), u(ndofel), shapeFunctions(nnode/2), shapeFuncDerivatives(D-1,nnode/2), &
                                      prop_thickness
      REAL(kind=AbqRK), INTENT(OUT):: Matrix_B(ntens,ndofel), JacobiDet, drot(3,3)
      REAL(kind=AbqRK) :: transposedJacobiMatrix(D,D), normalVect(D), &
                          Matrix_dN(D,nnode), Matrix_u(D,nnode), &
                          radius, currentNodeLocation(D,nnode)
      REAL(kind=AbqRK) :: Matrix_N(ntens,ndofel), Trafomat(ntens,ntens)
      INTEGER(kind=AbqIK) :: i1, i2, i3

      LOGICAL :: isSymmetric

      positivDetJ = .TRUE.; Matrix_B = zero; JacobiDet = zero; drot = zero; isSymmetric = .FALSE.

      ! symmetric cohesive elements
      IF (jprops(5) .EQ. 1) THEN
        isSymmetric = .TRUE.
      END IF

      ! erweiterte Matrix der Formfunktionen
      Matrix_N = zero
      ! Verschiebungssprung
      FORALL (i1=1:D, i2=1:nnode/2)
        Matrix_N( i1, (NodeMapMatrix(i2,1)-1)*ndofel/nnode+i1 ) =  shapeFunctions(i2) ! positiv
        Matrix_N( i1, (NodeMapMatrix(i2,2)-1)*ndofel/nnode+i1 ) = -shapeFunctions(i2) ! negativ
      END FORALL
      ! Mittelung der Schaedigungsvariablen
      FORALL (i2=1:nnode/2, i3=1:2)
        Matrix_N( D+1, (NodeMapMatrix(i2,i3))*ndofel/nnode ) = shapeFunctions(i2)/two ! half of top
      END FORALL
      ! Sprung der Schaedigungsvariablen
      FORALL (i2=1:nnode/2)
        Matrix_N( D+2, (NodeMapMatrix(i2,1))*ndofel/nnode ) =  shapeFunctions(i2) ! positiv
        Matrix_N( D+2, (NodeMapMatrix(i2,2))*ndofel/nnode ) = -shapeFunctions(i2) ! negativ
      END FORALL
      ! Gradient der Schaedigungsvariablen
      FORALL (i1=1:D-1, i2=1:nnode/2, i3=1:2)
        ! i1 - Richtung der Ableitung ! i2 - Zeile in NodeMapMatrix ! i3 - Spalte in NodeMapMatrix
        Matrix_N( D+2+i1, (NodeMapMatrix(i2,i3))*ndofel/nnode ) = shapeFuncDerivatives(i1,i2)/two
      END FORALL

      ! erweiterte Matrix der Ableitungen der Formfunktionen
      Matrix_dN = zero

      IF (isSymmetric) THEN
        ! oben, +, 1. Spalte in NodeMapMatrix
        !FORALL (i1=1:D-1, i2=1:nnode/2)
        !  ! i1 - Richtung der Ableitung ! i2 - Zeile in NodeMapMatrix
        !  Matrix_dN(i1,NodeMapMatrix(i2,1)) = zero
        !END FORALL
        ! unten, -, 2. Spalte in NodeMapMatrix
        FORALL (i1=1:D-1, i2=1:nnode/2)
          ! i1 - Richtung der Ableitung ! i2 - Zeile in NodeMapMatrix
          Matrix_dN(i1,NodeMapMatrix(i2,2)) = shapeFuncDerivatives(i1,i2)
        END FORALL
      ELSE
        FORALL (i1=1:D-1, i2=1:nnode/2, i3=1:2)
          ! i1 - Richtung der Ableitung ! i2 - Zeile in NodeMapMatrix ! i3 - Spalte in NodeMapMatrix
          Matrix_dN(i1,NodeMapMatrix(i2,i3)) = shapeFuncDerivatives(i1,i2)/two
        END FORALL
      END IF
      ! derivatives w.r.t. local z direction
      ! Matrix_dN(D,:) is only evaluated for axisymmetric case, otherwise no effect
      FORALL (i1=1:nnode/2)
        ! i1 - Zeile in NodeMapMatrix
        Matrix_dN(D,NodeMapMatrix(i1,1)) =  shapeFunctions(i1)/two ! positiv
        Matrix_dN(D,NodeMapMatrix(i1,2)) = -shapeFunctions(i1)/two ! negativ
      END FORALL

      ! Knoten-Koordinaten (aktuell oder nominal), Fallunterscheidung nach isNlgeom
      Matrix_u = zero

! jprops(3): 0 - nach NLGEOM
!            1 - linear
!            2 - nichtlinear
!      IF (isNlgeom) THEN
      IF ( ( (isNlgeom) .AND. (jprops(3) .NE. 1) ) .OR. ( (.NOT. isNlgeom) .AND. (jprops(3) .EQ. 2) ) ) THEN
        ! Verschiebungen in Matrixform
        FORALL (i1=1:nnode)
          Matrix_u(:,i1) = u((ndofel/nnode)*(i1-1)+1:(ndofel/nnode)*(i1-1)+D)
        END FORALL
      END IF
      currentNodeLocation = coords+Matrix_u

      ! Matrix der kovariante Basisvektoren (= transponierte JACOBI-Matrix)
      transposedJacobiMatrix = matmul(currentNodeLocation, transpose(Matrix_dN))

      ! Normalenvektor
      SELECT CASE(D)
      CASE(2)
        ! Normalenvektor in 2D (xy-Ebene) aus Vektorprodukt transposedJacobiMatrix(:,1) x (0,0,-1)
        normalVect(1:2) = (/ -transposedJacobiMatrix(2,1), transposedJacobiMatrix(1,1) /)
      CASE(3)
        ! Normalenvektor aus Kreuzprodukt der kovarianten Basisvektoren
        normalVect(1:3) = crossProduct(transposedJacobiMatrix(1:3,1),transposedJacobiMatrix(1:3,2))
      END SELECT

      ! JACOBI-Determinante
      JacobiDet = VectorNorm(normalVect,D)

      ! Abfrage, ob negative oder Null-Determinante
      IF (JacobiDet .LE. zero) THEN
        positivDetJ = .FALSE.
        ! Abbruch
      ELSE
        ! Normierung des Normalenvektors
        normalVect = normalVect/JacobiDet

        ! Normierung des ersten Basisvektors (1. Spalte in transposedJacobiMatrix)
        transposedJacobiMatrix(1:D,1) = transposedJacobiMatrix(1:D,1) / VectorNorm(transposedJacobiMatrix(1:D,1),D)

        ! Transformationsmatrix
        Trafomat = zero
        SELECT CASE(D)
        CASE(2)
          Trafomat(1,1:2) = transposedJacobiMatrix(1:2,1)
          Trafomat(2,1:2) = normalVect(1:2)
        CASE(3)
          Trafomat(1,1:3) = transposedJacobiMatrix(1:3,1)
          Trafomat(3,1:3) = normalVect(1:3)
          Trafomat(2,1:3) = crossProduct(Trafomat(3,1:3),Trafomat(1,1:3))
        END SELECT
        ! Einheitssubmatrix fuer Schaedigungsvariable, Sprung der Schaedigungsvariable und Gradient
        FORALL (i1=1:D+1)
          Trafomat(D+i1,D+i1) = one
        END FORALL

        ! Ausgabe-Rotationsmatrix
        drot = zero
        drot(1:D,1:D) = Trafomat(1:D,1:D)
        ! Hauptdiagonale mit 1 fuellen
        IF (D .LT. 3) THEN
          FORALL (i1=D+1:3)
            drot(i1,i1) = one
          END FORALL
        END IF

        ! B-Matrix
        Matrix_B = matmul(TrafoMat,Matrix_N)

        ! Halbierung bei Pentaedern
        IF ( (D .EQ. 3) .AND. ( (nnode .EQ. 6) .OR. (nnode .EQ. 12) ) ) THEN
          JacobiDet = JacobiDet/two
        END IF

        ! enhance JACOBI determinante for 2D
        IF (D .EQ. 2) THEN
        ! 2D
          IF (jprops(2) .EQ. 1) THEN
          ! axisymmetric
            ! multiply JacobiDet by 2*pi*r
            ! change sign in Matrix_dN --> shape functions of mean surface available
            FORALL (i1=1:nnode/2)
              Matrix_dN(D,NodeMapMatrix(i1,2)) = -Matrix_dN(D,NodeMapMatrix(i1,2))
            END FORALL
            radius = dot_product(currentNodeLocation(1,:),Matrix_dN(D,:))
            JacobiDet = JacobiDet*two*acos(-one)*radius
          ELSE
          ! plane strain
            ! multiply JacobiDet by thickness>0
            IF (prop_thickness .GT. zero) JacobiDet = JacobiDet*prop_thickness
          END IF

        END IF

      END IF ! Abfrage detJ

    END SUBROUTINE BMatrixJacobian

!------------------------------------------------------------------------------------

END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! NODE MAP
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! rectangle, linear shape functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE NodeMap_Rectangle_Linear
!     2 MEANNODES, NODE NUMBERING:
!     4  3
!     A  B
!     1  2

  USE ABQINTERFACE

  IMPLICIT NONE
  PRIVATE
  INTEGER(kind=AbqIK), PARAMETER :: NNODES=2
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: NodeMapMatrix(NNODES,2)=reshape( (/ 4, 3, &
                                                                                1, 2 /),(/ NNODES, 2 /) )
  !                                                                             A, B
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! rectangle, quadratic shape functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE NodeMap_Rectangle_Quadratic
!     3 MEANNODES, NODE NUMBERING:
!     4 6 3
!     A C B
!     1 5 2

  USE ABQINTERFACE

  IMPLICIT NONE
  PRIVATE
  INTEGER(kind=AbqIK), PARAMETER :: NNODES=3
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: NodeMapMatrix(NNODES,2)=reshape( (/ 4, 3, 6, &
                                                                                1, 2, 5 /),(/ NNODES, 2 /) )
  !                                                                             A, B, C
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! wedge, linear shape functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE NodeMap_Wedge_Linear
!     3 MEANNODES, NODE NUMBERING:
!     6
!     4  5
!     C
!     A  B
!     3
!     1  2

  USE ABQINTERFACE

  IMPLICIT NONE
  PRIVATE
  INTEGER(kind=AbqIK), PARAMETER :: NNODES=3
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: NodeMapMatrix(NNODES,2)=reshape( (/ 4, 5, 6, &
                                                                                1, 2, 3 /),(/ NNODES, 2 /) )
  !                                                                             A, B, C
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! wedge, quadratic shape functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE NodeMap_Wedge_Quadratic
!     6 MEANNODES, NODE NUMBERING:
!      6
!     12 11
!      4 10  5
!      C
!      F  E
!      A  D  B
!      3
!      9  8
!      1  7  2

  USE ABQINTERFACE

  IMPLICIT NONE
  PRIVATE
  INTEGER(kind=AbqIK), PARAMETER :: NNODES=6
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: NodeMapMatrix(NNODES,2)=reshape( (/ 4, 5, 6, 10, 11, 12, &
                                                                                1, 2, 3,  7,  8,  9 /),(/ NNODES, 2 /) )
  !                                                                             A, B, C,  D,  E,  F
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! hex, linear shape functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE NodeMap_Hex_Linear
!     4 MEANNODES, NODE NUMBERING:
!     8 7
!     5 6
!     D C
!     A B
!     4 3
!     1 2

  USE ABQINTERFACE

  IMPLICIT NONE
  PRIVATE
  INTEGER(kind=AbqIK), PARAMETER :: NNODES=4
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: NodeMapMatrix(NNODES,2)=reshape( (/ 5, 6, 7, 8, &
                                                                                1, 2, 3, 4 /),(/ NNODES, 2 /) )
  !                                                                             A, B, C, D
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! hex, quadratic shape functions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE NodeMap_Hex_Quadratic
!     8 MEANNODES, NODE NUMBERING:
!      8 15  7
!     16    14
!      5 13  6
!      D  G  C
!      H     F
!      A  E  B
!      4 11  3
!     12    10
!      1  9  2

  USE ABQINTERFACE

  IMPLICIT NONE
  PRIVATE
  INTEGER(kind=AbqIK), PARAMETER :: NNODES=8
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: NodeMapMatrix(NNODES,2)=reshape( (/ 5, 6, 7, 8, 13, 14, 15, 16, &
                                                                                1, 2, 3, 4,  9, 10, 11, 12 /),(/ NNODES, 2 /) )
  !                                                                             A, B, C, D,  E,  F,  G,  H
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
