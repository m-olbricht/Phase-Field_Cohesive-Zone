!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Basic modules for ABAQUS user subroutines
!
! Stephan Roth, TU Bergakademie Freiberg, 23.03.2012
!
! 09.10.2018: Vorbereitung Projektarbeit
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE ABQINTERFACE_PF

  INCLUDE 'ABA_PARAM.INC'
  PRIVATE
  INTEGER,PARAMETER,PUBLIC::AbqRK=KIND(r),AbqIK=KIND(i)

END MODULE ABQINTERFACE_PF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE FLOATNUMBERS

  IMPLICIT NONE

  DOUBLE PRECISION, PARAMETER :: zero=0.d0, one=1.d0, two=2.d0, three=3.d0, &
                                 four=4.d0, five=5.d0, six=6.d0, seven=7.d0, &
                                 eight=8.d0, nine=9.d0, ten=10.d0, twelve = 12.d0, &
                                 half=0.5d0, hundred=100.d0, twenty=20.d0

END MODULE FLOATNUMBERS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE ABAModul

  IMPLICIT NONE
  PUBLIC :: XEXIT

  CONTAINS

!------------------------------------------------------------------------------------

    SUBROUTINE XEXIT()
    ! Ausschrift und Abbruch der Simulation

      WRITE(7,*) 'EXIT by user subroutine'
      CALL XIT()

    END SUBROUTINE XEXIT

!------------------------------------------------------------------------------------

END MODULE ABAModul

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CZMath

  IMPLICIT NONE
  PUBLIC :: VectorNorm, crossProduct

  CONTAINS

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION VectorNorm(V,N)
    ! length of a vector, quadratic norm

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: N
      REAL(kind=AbqRK), INTENT(IN) :: V(N)
      INTEGER(kind=AbqIK) :: i1

      VectorNorm=zero
      DO i1=1,N
        VectorNorm=VectorNorm+V(i1)**2
      END DO
      VectorNorm=sqrt(VectorNorm)

    END FUNCTION VectorNorm

!------------------------------------------------------------------------------------

    PURE FUNCTION crossProduct(V1,V2)
    ! cross product of tw vectors

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: V1(3), V2(3)
      REAL(kind=AbqRK) :: crossProduct(3)

      crossProduct(1) = V1(2)*V2(3)-V1(3)*V2(2)
      crossProduct(2) = V1(3)*V2(1)-V1(1)*V2(3)
      crossProduct(3) = V1(1)*V2(2)-V1(2)*V2(1)

    END FUNCTION crossProduct

!------------------------------------------------------------------------------------

END MODULE CZMath

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE MathModul

  IMPLICIT NONE
  PUBLIC :: get_matrixinverse, determinante, Determinante2D3D, MacAulay, d_MacAulay

  CONTAINS

!------------------------------------------------------------------------------------

    PURE SUBROUTINE get_matrixinverse(D,A,A_inv)
    ! inverse of non-symmetric arbitrary shaped (D,D)-matrix

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D
      REAL(kind=AbqRK), INTENT(IN) :: A(D,D)
      REAL(kind=AbqRK), INTENT(OUT) :: A_inv(D,D)
      REAL(kind=AbqRK) :: subM(D,D)
      INTEGER(kind=AbqIK) :: i1, i2
      REAL(kind=AbqRK) :: det

      subM = zero; det = zero; A_inv = zero
! for symmetric matrices
!      DO i1=1,D
!        DO i2=i1,D
!          ! reduced matrix
!          subM = A
!          subM(i1:2,:)=A(i1+1:D,:)
!          subM(:,i2:2)=subM(:,i2+1:D)
!          ! Inverse ohne Division durch det(A) = Adjunktenmatrix
!          A_inv(i1,i2) = (-1.d0)**(i1+i2) * determinante(subM(1:2,1:2))
!          A_inv(i2,i1) = A_inv(i1,i2)
!        END DO
!      END DO

      DO i1=1,D
        DO i2=1,D
          ! reduced matrix
          subM = A
          IF (i1 .LE. (D-1)) subM(i1:(D-1),:)=A(i1+1:D,:)
          IF (i2 .LE. (D-1)) subM(:,i2:(D-1))=subM(:,i2+1:D)
          ! Inverse ohne Division durch det(A) = Adjunktenmatrix
          A_inv(i1,i2) = (-one)**(i1+i2) * determinante(subM(1:(D-1),1:(D-1)))
        END DO
      END DO
      ! Determinante ohne subroutine aufzurufen, Entwicklung nach 1. Zeile
      det = zero
      DO i1=1,D
        det = det + A(1,i1) * A_inv(1,i1)
      END DO
      ! Inverse
      A_inv = TRANSPOSE(A_inv) / det

    END SUBROUTINE get_matrixinverse

!------------------------------------------------------------------------------------

    PURE RECURSIVE REAL(kind=AbqRK) FUNCTION determinante(AA)

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: AA(:,:)
      REAL(kind=AbqRK), ALLOCATABLE :: subM(:,:)
      INTEGER(kind=AbqIK) :: i3, n_size

      n_size = size(AA,DIM=1)
      IF (n_size .EQ. 1) THEN
        determinante = AA(1,1)
      ELSE IF (n_size .EQ. 2) THEN
        determinante = AA(1,1)*AA(2,2)-AA(1,2)*AA(2,1)
      ELSE IF (n_size .EQ. 3) THEN
        determinante = AA(1,1)*AA(2,2)*AA(3,3)+AA(1,2)*AA(2,3)*AA(3,1)+ &
                       AA(1,3)*AA(2,1)*AA(3,2)-AA(1,3)*AA(2,2)*AA(3,1)- &
                       AA(1,2)*AA(2,1)*AA(3,3)-AA(1,1)*AA(2,3)*AA(3,2)
      ELSE
        ALLOCATE(subM(n_size,n_size)) ! matrix dimension
        determinante = zero
        DO i3=1,n_size
          ! reduced matrix
          subM = AA
          subM(:,i3:n_size-1)=AA(:,i3+1:n_size)
          ! development w.r.t. last row
          determinante = determinante + AA(n_size,i3)*(-one)**(n_size+i3)*determinante(subM(1:n_size-1,1:n_size-1))
        END DO
        DEALLOCATE(subM)
      END IF

    END FUNCTION determinante

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION Determinante2D3D(M,N)
    ! determinant of 2x2- or 3x3-matrix

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: N
      REAL(kind=AbqRK), INTENT(IN) :: M(N,N)

      IF (N .EQ. 2) THEN
        Determinante2D3D = M(1,1)*M(2,2)-M(1,2)*M(2,1)
      ELSE IF (N .EQ. 3) THEN
        Determinante2D3D = M(1,1)*M(2,2)*M(3,3)+M(1,2)*M(2,3)*M(3,1)+ &
                           M(1,3)*M(2,1)*M(3,2)-M(1,3)*M(2,2)*M(3,1)- &
                           M(1,2)*M(2,1)*M(3,3)-M(1,1)*M(2,3)*M(3,2)
      END IF

    END FUNCTION Determinante2D3D

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION MacAulay(arg)
    ! MacAulay-Operator: 0.5*(abs(arg)+arg)

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: arg
      !
      MacAulay = (ABS(arg)+arg)/two
 
    END FUNCTION MacAulay

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION d_MacAulay(arg)
    ! first derivative of MacAulay-Operator: 0.5*(sgn(arg)+1) = Heaviside-Jump-Function

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: arg
      !
      d_MacAulay = (SIGN(one,arg)+one)/two
 
    END FUNCTION d_MacAulay

!------------------------------------------------------------------------------------

END MODULE MathModul

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

