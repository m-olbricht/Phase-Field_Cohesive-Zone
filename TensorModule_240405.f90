!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! tensor modules
!
! Stephan Roth, TU Bergakademie Freiberg, 20.03.2017
!
! 20.03.2017: first implementation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE TensorModule

  IMPLICIT NONE

!  PUBLIC :: 

  CONTAINS

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION trace(tensor)
    ! trace of second order tensor

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: tensor(3,3)
      !
      trace = tensor(1,1)+tensor(2,2)+tensor(3,3)

    END FUNCTION trace

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION deviator(tensor)
    ! deviator of second order tensor

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: tensor(3,3)
      INTEGER(kind=AbqIK) :: i1
      REAL(kind=AbqRK) :: tr
      DIMENSION deviator(3,3)
      !
      tr = trace(tensor)
      deviator = tensor
      FORALL(i1=1:3)
        deviator(i1,i1) = deviator(i1,i1) - tr/three
      END FORALL

    END FUNCTION deviator

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION Kronecker(i,j)
    ! Kronecker symbol 

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: i, j

      Kronecker = zero
      IF (i .EQ. j) THEN
        Kronecker = one
      END IF

    END FUNCTION Kronecker

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION Permutation()
    ! Permutationstensor

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK) :: i1, i2, i3
      REAL(kind=AbqRK) :: temp3(3,3,3)
      DIMENSION Permutation(3,3,3)

      Permutation = zero
      temp3 = zero
      FORALL (i1=1:3, i2=1:3, i3=1:3, ( ( (i1 .EQ. 1) .AND. (i2 .EQ. 2) .AND. (i3 .EQ. 3) ) .OR. ( (i1 .EQ. 2) .AND. (i2 .EQ. 3) .AND. (i3 .EQ. 1) ) .OR. ( (i1 .EQ. 3) .AND. (i2 .EQ. 1) .AND. (i3 .EQ. 2) ) ) ) temp3(i1,i2,i3) =  one
      Permutation = temp3
      temp3 = zero
      FORALL (i1=1:3, i2=1:3, i3=1:3, ( ( (i1 .EQ. 1) .AND. (i2 .EQ. 3) .AND. (i3 .EQ. 2) ) .OR. ( (i1 .EQ. 2) .AND. (i2 .EQ. 1) .AND. (i3 .EQ. 3) ) .OR. ( (i1 .EQ. 3) .AND. (i2 .EQ. 2) .AND. (i3 .EQ. 1) ) ) ) temp3(i1,i2,i3) = -one
      Permutation = Permutation + temp3

    END FUNCTION Permutation

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION DYADE(FirstOrderTensor_A,FirstOrderTensor_B)
    ! dyade of a first order tensor and a first order tensor: Cij=A_ixB_j

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: FirstOrderTensor_A(3), FirstOrderTensor_B(3)
      INTEGER(kind=AbqIK) :: i1, i2
      DIMENSION DYADE(3,3)

      DYADE=zero
      FORALL (i1=1:3)
        FORALL (i2=1:3)
          DYADE(i1,i2)=FirstOrderTensor_A(i1)*FirstOrderTensor_B(i2)
        END FORALL
      END FORALL

    END FUNCTION DYADE

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION DYADETwoTwo(SecondOrderTensor_A,SecondOrderTensor_B)
    ! dyade of a second order tensor and a second order tensor: Cijkl=A_ij x B_kl

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: SecondOrderTensor_A(3,3), SecondOrderTensor_B(3,3)
      INTEGER(kind=AbqIK) :: i1, i2, i3, i4
      DIMENSION DYADETwoTwo(3,3,3,3)

      DYADETwoTwo=zero
      FORALL (i1=1:3)
        FORALL (i2=1:3)
          FORALL (i3=1:3)
            FORALL (i4=1:3)
              DYADETwoTwo(i1,i2,i3,i4)=SecondOrderTensor_A(i1,i2)*SecondOrderTensor_B(i3,i4)
            END FORALL
          END FORALL
        END FORALL
      END FORALL

    END FUNCTION DYADETwoTwo

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION DYADETwoOne(SecondOrderTensor,FirstOrderTensor)
    ! dyade of a second order tensor and a first order tensor: Cijk=A_ijxB_k

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: SecondOrderTensor(3,3), FirstOrderTensor(3)
      INTEGER(kind=AbqIK) :: i1, i2, i3
      DIMENSION DYADETwoOne(3,3,3)

      DYADETwoOne=zero
      FORALL (i1=1:3)
        FORALL (i2=1:3)
          FORALL (i3=1:3)
            DYADETwoOne(i1,i2,i3)=SecondOrderTensor(i1,i2)*FirstOrderTensor(i3)
          END FORALL
        END FORALL
      END FORALL

    END FUNCTION DYADETwoOne

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION DYADEOneTwo(FirstOrderTensor,SecondOrderTensor)
    ! dyade of a first order tensor and a second order tensor: Ckij=B_kxA_ij

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: FirstOrderTensor(3), SecondOrderTensor(3,3)
      INTEGER(kind=AbqIK) :: i1, i2, i3
      DIMENSION DYADEOneTwo(3,3,3)

      DYADEOneTwo=zero
      FORALL (i1=1:3)
        FORALL (i2=1:3)
          FORALL (i3=1:3)
            DYADEOneTwo(i1,i2,i3)=FirstOrderTensor(i1)*SecondOrderTensor(i2,i3)
          END FORALL
        END FORALL
      END FORALL

    END FUNCTION DYADEOneTwo

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION DYADEOneTwoSecond(FirstOrderTensor,SecondOrderTensor)
    ! dyade of a first order tensor and a second order tensor: Cijk = A_j x B_ik

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: FirstOrderTensor(3), SecondOrderTensor(3,3)
      INTEGER(kind=AbqIK) :: i1, i2, i3
      DIMENSION DYADEOneTwoSecond(3,3,3)

      DYADEOneTwoSecond=zero
      FORALL (i1=1:3)
        FORALL (i2=1:3)
          FORALL (i3=1:3)
            DYADEOneTwoSecond(i1,i2,i3)=FirstOrderTensor(i2)*SecondOrderTensor(i1,i3)
          END FORALL
        END FORALL
      END FORALL

    END FUNCTION DYADEOneTwoSecond

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION DYADETwoOneSecond(SecondOrderTensor,FirstOrderTensor)
    ! dyade of a first order tensor and a second order tensor: Cijk = A_ik x B_j

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: SecondOrderTensor(3,3), FirstOrderTensor(3)
      INTEGER(kind=AbqIK) :: i1, i2, i3
      DIMENSION DYADETwoOneSecond(3,3,3)

      DYADETwoOneSecond=zero
      FORALL (i1=1:3)
        FORALL (i2=1:3)
          FORALL (i3=1:3)
            DYADETwoOneSecond(i1,i2,i3)=SecondOrderTensor(i1,i3)*FirstOrderTensor(i2)
          END FORALL
        END FORALL
      END FORALL

    END FUNCTION DYADETwoOneSecond

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION DYADEThreeOne(ThirdOrderTensor,FirstOrderTensor)
    ! dyade of a third order tensor and a first order tensor: Cijkl=A_ijk x B_l

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: ThirdOrderTensor(3,3,3), FirstOrderTensor(3)
      INTEGER(kind=AbqIK) :: i1, i2, i3, i4
      DIMENSION DYADEThreeOne(3,3,3,3)

      DYADEThreeOne=zero
      FORALL (i1=1:3)
        FORALL (i2=1:3)
          FORALL (i3=1:3)
            FORALL (i4=1:3)
              DYADEThreeOne(i1,i2,i3,i4)=ThirdOrderTensor(i1,i2,i3)*FirstOrderTensor(i4)
            END FORALL
          END FORALL
        END FORALL
      END FORALL

    END FUNCTION DYADEThreeOne

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION DYADEThreeOneSecond(ThirdOrderTensor,FirstOrderTensor)
    ! dyade of a third order tensor and a first order tensor: Cijkl=A_ikl x B_j

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: ThirdOrderTensor(3,3,3), FirstOrderTensor(3)
      INTEGER(kind=AbqIK) :: i1, i2, i3, i4
      DIMENSION DYADEThreeOneSecond(3,3,3,3)

      DYADEThreeOneSecond=zero
      FORALL (i1=1:3)
        FORALL (i2=1:3)
          FORALL (i3=1:3)
            FORALL (i4=1:3)
              DYADEThreeOneSecond(i1,i2,i3,i4)=ThirdOrderTensor(i1,i3,i4)*FirstOrderTensor(i2)
            END FORALL
          END FORALL
        END FORALL
      END FORALL

    END FUNCTION DYADEThreeOneSecond

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION DYADEThreeOneThird(ThirdOrderTensor,FirstOrderTensor)
    ! dyade of a third order tensor and a first order tensor: Cijkl=A_ijl x B_k

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: ThirdOrderTensor(3,3,3), FirstOrderTensor(3)
      INTEGER(kind=AbqIK) :: i1, i2, i3, i4
      DIMENSION DYADEThreeOneThird(3,3,3,3)

      DYADEThreeOneThird=zero
      FORALL (i1=1:3)
        FORALL (i2=1:3)
          FORALL (i3=1:3)
            FORALL (i4=1:3)
              DYADEThreeOneThird(i1,i2,i3,i4)=ThirdOrderTensor(i1,i2,i4)*FirstOrderTensor(i3)
            END FORALL
          END FORALL
        END FORALL
      END FORALL

    END FUNCTION DYADEThreeOneThird

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION DYADEOneThreeFirst(FirstOrderTensor,ThirdOrderTensor)
    ! dyade of a first order tensor and a third order tensor: Cijkl=A_ixB_jkl

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: FirstOrderTensor(3), ThirdOrderTensor(3,3,3)
      INTEGER(kind=AbqIK) :: i1, i2, i3, i4
      DIMENSION DYADEOneThreeFirst(3,3,3,3)

      DYADEOneThreeFirst=zero
      FORALL (i1=1:3)
        FORALL (i2=1:3)
          FORALL (i3=1:3)
            FORALL (i4=1:3)
              DYADEOneThreeFirst(i1,i2,i3,i4)=FirstOrderTensor(i1)*ThirdOrderTensor(i2,i3,i4)
            END FORALL
          END FORALL
        END FORALL
      END FORALL

    END FUNCTION DYADEOneThreeFirst

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION TransposedThirdOrderTensor(ThirdOrderTensor)
    ! transpose third order tensor: B_ikl = A_kli 

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: ThirdOrderTensor(3,3,3)
      INTEGER(kind=AbqIK) :: i1, i2, i3
      DIMENSION TransposedThirdOrderTensor(3,3,3)

      TransposedThirdOrderTensor=zero
      FORALL (i1=1:3)
        FORALL (i2=1:3)
          FORALL (i3=1:3)
            TransposedThirdOrderTensor(i3,i1,i2) = ThirdOrderTensor(i1,i2,i3)
          END FORALL
        END FORALL
      END FORALL

    END FUNCTION TransposedThirdOrderTensor

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION TransposedForthOrderTensor(ForthOrderTensor)
    ! transpose forth order tensor: B_ijkl = A_klij 

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: ForthOrderTensor(3,3,3,3)
      INTEGER(kind=AbqIK) :: i1, i2, i3, i4
      DIMENSION TransposedForthOrderTensor(3,3,3,3)

      TransposedForthOrderTensor=zero
      FORALL (i1=1:3)
        FORALL (i2=1:3)
          FORALL (i3=1:3)
            FORALL (i4=1:3)
              TransposedForthOrderTensor(i1,i2,i3,i4) = ForthOrderTensor(i3,i4,i1,i2)
            END FORALL
          END FORALL
        END FORALL
      END FORALL

    END FUNCTION TransposedForthOrderTensor

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION SINGLECONTRACTIONOneOne(FirstOrderTensor_A,FirstOrderTensor_B)
    ! single contraction of a first order tensor and a first order tensor: C=A_i*B_i

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: FirstOrderTensor_A(3), FirstOrderTensor_B(3)

      SINGLECONTRACTIONOneOne=dot_product(FirstOrderTensor_A(:),FirstOrderTensor_B(:))

    END FUNCTION SINGLECONTRACTIONOneOne

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION SINGLECONTRACTIONOneTwo(FirstOrderTensor,SecondOrderTensor)
    ! single contraction of a first order tensor and a second order tensor: C_j=A_i*B_ij

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: FirstOrderTensor(3), SecondOrderTensor(3,3)
      INTEGER(kind=AbqIK) :: i1
      DIMENSION SINGLECONTRACTIONOneTwo(3)

      SINGLECONTRACTIONOneTwo=zero
      FORALL (i1=1:3)
        SINGLECONTRACTIONOneTwo(i1)=SINGLECONTRACTIONOneOne(FirstOrderTensor(:),SecondOrderTensor(:,i1))
      END FORALL

    END FUNCTION SINGLECONTRACTIONOneTwo

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION SINGLECONTRACTIONTwoOne(SecondOrderTensor,FirstOrderTensor)
    ! single contraction of a second order tensor and a first order tensor: C_i=A_ij*B_j

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: SecondOrderTensor(3,3), FirstOrderTensor(3)
      INTEGER(kind=AbqIK) :: i1
      DIMENSION SINGLECONTRACTIONTwoOne(3)

      SINGLECONTRACTIONTwoOne=zero
      FORALL (i1=1:3)
        SINGLECONTRACTIONTwoOne(i1)=SINGLECONTRACTIONOneOne(SecondOrderTensor(i1,:),FirstOrderTensor(:))
      END FORALL

    END FUNCTION SINGLECONTRACTIONTwoOne

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION SINGLECONTRACTIONTwoThree(SecondOrderTensor,ThirdOrderTensor)
    ! single contraction of a second order tensor and a third order tensor: C_ikl=A_ij*B_jkl

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: SecondOrderTensor(3,3), ThirdOrderTensor(3,3,3)
      INTEGER(kind=AbqIK) :: i1, i2, i3
      DIMENSION SINGLECONTRACTIONTwoThree(3,3,3)

      SINGLECONTRACTIONTwoThree=zero
      FORALL (i1=1:3)
        FORALL (i2=1:3)
          FORALL (i3=1:3)
            SINGLECONTRACTIONTwoThree(i1,i2,i3)=SINGLECONTRACTIONOneOne(SecondOrderTensor(i1,:),ThirdOrderTensor(:,i2,i3))
          END FORALL
        END FORALL
      END FORALL

    END FUNCTION SINGLECONTRACTIONTwoThree

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION SINGLECONTRACTIONOneFour(FirstOrderTensor,ForthOrderTensor)
    ! single contraction of a first order tensor and a forth order tensor: C_jkl=A_i*B_ijkl

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: FirstOrderTensor(3), ForthOrderTensor(3,3,3,3)
      INTEGER(kind=AbqIK) :: i1, i2, i3
      DIMENSION SINGLECONTRACTIONOneFour(3,3,3)

      SINGLECONTRACTIONOneFour=zero
      FORALL (i1=1:3)
        FORALL (i2=1:3)
          FORALL (i3=1:3)
            SINGLECONTRACTIONOneFour(i1,i2,i3)=SINGLECONTRACTIONOneOne(FirstOrderTensor(:),ForthOrderTensor(:,i1,i2,i3))
          END FORALL
        END FORALL
      END FORALL

    END FUNCTION SINGLECONTRACTIONOneFour

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION SINGLECONTRACTIONFourOneThirdIndex(ForthOrderTensor,FirstOrderTensor)
    ! single contraction of a forth order tensor and a first order tensor: C_ijl=B_ijkl*A_k

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: ForthOrderTensor(3,3,3,3), FirstOrderTensor(3)
      INTEGER(kind=AbqIK) :: i1, i2, i3
      DIMENSION  SINGLECONTRACTIONFourOneThirdIndex(3,3,3)

      SINGLECONTRACTIONFourOneThirdIndex=zero
      FORALL (i1=1:3)
        FORALL (i2=1:3)
          FORALL (i3=1:3)
             SINGLECONTRACTIONFourOneThirdIndex(i1,i2,i3)=SINGLECONTRACTIONOneOne(ForthOrderTensor(i1,i2,:,i3),FirstOrderTensor(:))
          END FORALL
        END FORALL
      END FORALL

    END FUNCTION  SINGLECONTRACTIONFourOneThirdIndex

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION SINGLECONTRACTIONThreeOne(ThirdOrderTensor,FirstOrderTensor)
    ! single contraction of a third order tensor and a first order tensor: C_ij=A_ijk*B_k

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: ThirdOrderTensor(3,3,3), FirstOrderTensor(3)
      INTEGER(kind=AbqIK) :: i1, i2
      DIMENSION SINGLECONTRACTIONThreeOne(3,3)

      SINGLECONTRACTIONThreeOne=zero
      FORALL (i1=1:3)
        FORALL (i2=1:3)
          SINGLECONTRACTIONThreeOne(i1,i2)=SINGLECONTRACTIONOneOne(ThirdOrderTensor(i1,i2,:),FirstOrderTensor(:))
        END FORALL
      END FORALL

    END FUNCTION SINGLECONTRACTIONThreeOne

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION SINGLECONTRACTIONOneThree(FirstOrderTensor,ThirdOrderTensor)
    ! single contraction of a first order tensor and a third order tensor: C_ij=A_k*B_kij

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: FirstOrderTensor(3), ThirdOrderTensor(3,3,3)
      INTEGER(kind=AbqIK) :: i1, i2
      DIMENSION SINGLECONTRACTIONOneThree(3,3)

      SINGLECONTRACTIONOneThree=zero
      FORALL (i1=1:3)
        FORALL (i2=1:3)
          SINGLECONTRACTIONOneThree(i1,i2)=SINGLECONTRACTIONOneOne(FirstOrderTensor(:),ThirdOrderTensor(:,i1,i2))
        END FORALL
      END FORALL

    END FUNCTION SINGLECONTRACTIONOneThree

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION SINGLECONTRACTIONThreeTwo(ThirdOrderTensor,SecondOrderTensor)
    ! single contraction of a third order tensor and a second order tensor: C_ijk=A_ijp*B_pk

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: ThirdOrderTensor(3,3,3), SecondOrderTensor(3,3)
      INTEGER(kind=AbqIK) :: i1, i2, i3
      DIMENSION SINGLECONTRACTIONThreeTwo(3,3,3)

      SINGLECONTRACTIONThreeTwo=zero
      FORALL (i1=1:3)
        FORALL (i2=1:3)
          FORALL (i3=1:3)
            SINGLECONTRACTIONThreeTwo(i1,i2,i3)=SINGLECONTRACTIONOneOne(ThirdOrderTensor(i1,i2,:),SecondOrderTensor(:,i3))
          END FORALL
        END FORALL
      END FORALL

    END FUNCTION SINGLECONTRACTIONThreeTwo

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION SINGLECONTRACTIONThreeThree(ThirdOrderTensor_A,ThirdOrderTensor_B)
    ! single contraction of a third order tensor and a third order tensor: C_ijkl=A_ijp*B_pkl

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: ThirdOrderTensor_A(3,3,3), ThirdOrderTensor_B(3,3,3)
      INTEGER(kind=AbqIK) :: i1, i2, i3, i4
      DIMENSION SINGLECONTRACTIONThreeThree(3,3,3,3)

      SINGLECONTRACTIONThreeThree=zero
      FORALL (i1=1:3)
        FORALL (i2=1:3)
          FORALL (i3=1:3)
            FORALL (i4=1:3)
              SINGLECONTRACTIONThreeThree(i1,i2,i3,i4)=SINGLECONTRACTIONOneOne(ThirdOrderTensor_A(i1,i2,:),ThirdOrderTensor_B(:,i3,i4))
            END FORALL
          END FORALL
        END FORALL
      END FORALL

    END FUNCTION SINGLECONTRACTIONThreeThree

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION DOUBLECONTRACTIONTwoTwo(SecondOrderTensor_A,SecondOrderTensor_B)
    ! double contraction of a second order tensor and a second order tensor: C=A_ij*B_ij

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: SecondOrderTensor_A(3,3), SecondOrderTensor_B(3,3)
      INTEGER(kind=AbqIK) :: i1, i2

      DOUBLECONTRACTIONTwoTwo=zero
      DO i1=1,3
        DO i2=1,3
          DOUBLECONTRACTIONTwoTwo=DOUBLECONTRACTIONTwoTwo+SecondOrderTensor_A(i1,i2)*SecondOrderTensor_B(i1,i2)
        END DO
      END DO

    END FUNCTION DOUBLECONTRACTIONTwoTwo

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION DOUBLECONTRACTIONThreeThree(ThirdOrderTensor_A,ThirdOrderTensor_B)
    ! double contraction of a third order tensor and a third order tensor: C_kl=A_kij*B_ijl

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: ThirdOrderTensor_A(3,3,3), ThirdOrderTensor_B(3,3,3)
      DIMENSION DOUBLECONTRACTIONThreeThree(3,3)
      INTEGER(kind=AbqIK) :: i1, i2

      DOUBLECONTRACTIONThreeThree=zero
      FORALL (i1=1:3)
        FORALL (i2=1:3)
          DOUBLECONTRACTIONThreeThree(i1,i2) = DOUBLECONTRACTIONTwoTwo(ThirdOrderTensor_A(i1,:,:),ThirdOrderTensor_B(:,:,i2))
        END FORALL
      END FORALL

    END FUNCTION DOUBLECONTRACTIONThreeThree

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION DOUBLECONTRACTIONThreeThreeThird(ThirdOrderTensor_A,ThirdOrderTensor_B)
    ! double contraction of a third order tensor and a third order tensor: C_kl=A_ijk*B_ijl

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: ThirdOrderTensor_A(3,3,3), ThirdOrderTensor_B(3,3,3)
      DIMENSION DOUBLECONTRACTIONThreeThreeThird(3,3)
      INTEGER(kind=AbqIK) :: i1, i2

      DOUBLECONTRACTIONThreeThreeThird=zero
      FORALL (i1=1:3)
        FORALL (i2=1:3)
          DOUBLECONTRACTIONThreeThreeThird(i1,i2) = DOUBLECONTRACTIONTwoTwo(ThirdOrderTensor_A(:,:,i1),ThirdOrderTensor_B(:,:,i2))
        END FORALL
      END FORALL

    END FUNCTION DOUBLECONTRACTIONThreeThreeThird

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION DOUBLECONTRACTIONFourFour(ForthOrderTensor_A,ForthOrderTensor_B)
    ! double contraction of a forth order tensor and a forth order tensor: C_klmn=A_klij*B_ijmn

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: ForthOrderTensor_A(3,3,3,3), ForthOrderTensor_B(3,3,3,3)
      DIMENSION DOUBLECONTRACTIONFourFour(3,3,3,3)
      INTEGER(kind=AbqIK) :: i1, i2, i3, i4

      DOUBLECONTRACTIONFourFour=zero
      FORALL (i1=1:3)
        FORALL (i2=1:3)
          FORALL (i3=1:3)
            FORALL (i4=1:3)
              DOUBLECONTRACTIONFourFour(i1,i2,i3,i4) = DOUBLECONTRACTIONTwoTwo(ForthOrderTensor_A(i1,i2,:,:),ForthOrderTensor_B(:,:,i3,i4))
            END FORALL
          END FORALL
        END FORALL
      END FORALL

    END FUNCTION DOUBLECONTRACTIONFourFour

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION DOUBLECONTRACTIONThreeFour(ThirdOrderTensor,ForthOrderTensor)
    ! double contraction of a third order tensor and a forth order tensor: C_klm=A_kij*B_ijlm

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: ThirdOrderTensor(3,3,3), ForthOrderTensor(3,3,3,3)
      DIMENSION DOUBLECONTRACTIONThreeFour(3,3,3)
      INTEGER(kind=AbqIK) :: i1, i2, i3

      DOUBLECONTRACTIONThreeFour=zero
      FORALL (i1=1:3)
        FORALL (i2=1:3)
          FORALL (i3=1:3)
            DOUBLECONTRACTIONThreeFour(i1,i2,i3) = DOUBLECONTRACTIONTwoTwo(ThirdOrderTensor(i1,:,:),ForthOrderTensor(:,:,i2,i3))
          END FORALL
        END FORALL
      END FORALL

    END FUNCTION DOUBLECONTRACTIONThreeFour

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION DOUBLECONTRACTIONFourThree(ForthOrderTensor,ThirdOrderTensor)
    ! double contraction of a forth order tensor and a third order tensor: C_klm=A_klij*B_ijm

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: ForthOrderTensor(3,3,3,3), ThirdOrderTensor(3,3,3)
      DIMENSION DOUBLECONTRACTIONFourThree(3,3,3)
      INTEGER(kind=AbqIK) :: i1, i2, i3

      DOUBLECONTRACTIONFourThree=zero
      FORALL (i1=1:3)
        FORALL (i2=1:3)
          FORALL (i3=1:3)
            DOUBLECONTRACTIONFourThree(i1,i2,i3) = DOUBLECONTRACTIONTwoTwo(ForthOrderTensor(i1,i2,:,:),ThirdOrderTensor(:,:,i3))
          END FORALL
        END FORALL
      END FORALL

    END FUNCTION DOUBLECONTRACTIONFourThree

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION DOUBLECONTRACTIONFourTwo(ForthOrderTensor,SecondOrderTensor)
    ! double contraction of a forth order tensor and a second order tensor: C_kl=A_klij*B_ij

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: ForthOrderTensor(3,3,3,3), SecondOrderTensor(3,3)
      DIMENSION DOUBLECONTRACTIONFourTwo(3,3)
      INTEGER(kind=AbqIK) :: i1, i2

      DOUBLECONTRACTIONFourTwo=zero
      FORALL (i1=1:3)
        FORALL (i2=1:3)
          DOUBLECONTRACTIONFourTwo(i1,i2) = DOUBLECONTRACTIONTwoTwo(ForthOrderTensor(i1,i2,:,:),SecondOrderTensor)
        END FORALL
      END FORALL

    END FUNCTION DOUBLECONTRACTIONFourTwo

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION DOUBLECONTRACTIONTwoSix(SecondOrderTensor,SixthOrderTensor)
    ! double contraction of a second order tensor and a sixth order tensor: C_klmn=A_ij*B_ijklmn

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: SecondOrderTensor(3,3), SixthOrderTensor(3,3,3,3,3,3)
      DIMENSION DOUBLECONTRACTIONTwoSix(3,3,3,3)
      INTEGER(kind=AbqIK) :: i1, i2, i3, i4

      DOUBLECONTRACTIONTwoSix=zero
      FORALL (i1=1:3)
        FORALL (i2=1:3)
          FORALL (i3=1:3)
            FORALL (i4=1:3)
              DOUBLECONTRACTIONTwoSix(i1,i2,i3,i4) = DOUBLECONTRACTIONTwoTwo(SecondOrderTensor,SixthOrderTensor(:,:,i1,i2,i3,i4))
            END FORALL
          END FORALL
        END FORALL
      END FORALL

    END FUNCTION DOUBLECONTRACTIONTwoSix

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION DOUBLECONTRACTIONTwoFive(SecondOrderTensor,FifthOrderTensor)
    ! double contraction of a second order tensor and a fifth order tensor: C_klm=A_ij*B_ijklm

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: SecondOrderTensor(3,3), FifthOrderTensor(3,3,3,3,3)
      DIMENSION DOUBLECONTRACTIONTwoFive(3,3,3)
      INTEGER(kind=AbqIK) :: i1, i2, i3

      DOUBLECONTRACTIONTwoFive=zero
      FORALL (i1=1:3)
        FORALL (i2=1:3)
          FORALL (i3=1:3)
            DOUBLECONTRACTIONTwoFive(i1,i2,i3) = DOUBLECONTRACTIONTwoTwo(SecondOrderTensor,FifthOrderTensor(:,:,i1,i2,i3))
          END FORALL
        END FORALL
      END FORALL

    END FUNCTION DOUBLECONTRACTIONTwoFive

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION DOUBLECONTRACTIONTwoFour(SecondOrderTensor,ForthOrderTensor)
    ! double contraction of a second order tensor and a forth order tensor: C_kl=A_ij*B_ijkl

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: SecondOrderTensor(3,3), ForthOrderTensor(3,3,3,3)
      DIMENSION DOUBLECONTRACTIONTwoFour(3,3)
      INTEGER(kind=AbqIK) :: i1, i2

      DOUBLECONTRACTIONTwoFour=zero
      FORALL (i1=1:3)
        FORALL (i2=1:3)
          DOUBLECONTRACTIONTwoFour(i1,i2) = DOUBLECONTRACTIONTwoTwo(SecondOrderTensor,ForthOrderTensor(:,:,i1,i2))
        END FORALL
      END FORALL

    END FUNCTION DOUBLECONTRACTIONTwoFour

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION DOUBLECONTRACTIONTwoThree(SecondOrderTensor,ThirdOrderTensor)
    ! double contraction of a second order tensor and a third order tensor: C_k=A_ij*B_ijk

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: SecondOrderTensor(3,3), ThirdOrderTensor(3,3,3)
      DIMENSION DOUBLECONTRACTIONTwoThree(3)
      INTEGER(kind=AbqIK) :: i1

      DOUBLECONTRACTIONTwoThree=zero
      FORALL (i1=1:3)
        DOUBLECONTRACTIONTwoThree(i1) = DOUBLECONTRACTIONTwoTwo(SecondOrderTensor,ThirdOrderTensor(:,:,i1))
      END FORALL

    END FUNCTION DOUBLECONTRACTIONTwoThree

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION DOUBLECONTRACTIONSTwoFourTwo(SecondOrderTensor_A,ForthOrderTensor,SecondOrderTensor_B)
    ! scalar eps_kl:E_klmn:eps_mn
    ! two double contractions: a second order tensor, a forth order tensor, and a second order tensor: D=A_ij*B_ijkl*C_kl

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: SecondOrderTensor_A(3,3), ForthOrderTensor(3,3,3,3), SecondOrderTensor_B(3,3)

      DOUBLECONTRACTIONSTwoFourTwo=DOUBLECONTRACTIONTwoTwo(SecondOrderTensor_A,DOUBLECONTRACTIONFourTwo(ForthOrderTensor,SecondOrderTensor_B))

    END FUNCTION DOUBLECONTRACTIONSTwoFourTwo

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION ElastTensor(E,nu)
    ! isotopic linear-elastic tensor E

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: E, nu
      REAL(kind=AbqRK) :: lambda, mu
      DIMENSION ElastTensor(3,3,3,3)
      INTEGER(kind=AbqIK) :: i, j, k, l

      lambda=nu*E/(one+nu)/(one-two*nu)
      mu=E/two/(one+nu)
      ElastTensor=zero
      DO i=1,3
        DO j=1,3
          DO k=1,3
            DO l=1,3
              IF ((i .EQ. j) .AND. (k .EQ. l)) THEN
                ElastTensor(i,j,k,l) = ElastTensor(i,j,k,l) + lambda
              END IF
              IF ((i .EQ. k) .AND. (j .EQ. l)) THEN
                ElastTensor(i,j,k,l) = ElastTensor(i,j,k,l) + mu
              END IF
              IF ((i .EQ. l) .AND. (j .EQ. k)) THEN
                ElastTensor(i,j,k,l) = ElastTensor(i,j,k,l) + mu
              END IF
            END DO
          END DO
        END DO
      END DO

    END FUNCTION ElastTensor

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION VolumetricTensor(Q)
    ! volumetric tensor

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: Q
      DIMENSION VolumetricTensor(3,3)
      INTEGER(kind=AbqIK) :: i1

      VolumetricTensor=zero
      FORALL (i1=1:3)
        VolumetricTensor(i1,i1)=Q
      END FORALL

    END FUNCTION VolumetricTensor

!------------------------------------------------------------------------------------

END MODULE TensorModule

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

