!DEC$ FREEFORM
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
! Library for User Elements in Abaqus providing modules for
! - Shape functions
! - B-Matrices
! - Quadrature rules
!
! Geralf Hutter,
! Rostyslav Skrypnyk, 28.02.2014
! Stephan Roth, 18.10.2017
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
! 26.06.2014: added shape functions and quadrature rules for tetrahedra
! 05.04.2016: corrected shape functions for 3-node (quadratic) line-elements (St. Roth)
! 18.10.2017: B-matrices for axi- and spheri-symmetric elements (St. Roth)
!-----------------------------------------------------------------------------------
!
!====================================================================================
!====================================================================================
!
! B-Matrices: Module provides B-Matrices for
!     - Vectors (BMatVect)
!     - symmetric tensors of second order in Voigt notation (BMatTensSymVoigt)
!     - unsymmetric tensors of second order in Voigt notation (BMatTensUnsym), e.g.
!       for the deformation gradient
!   in each case for different dimensionalities. The functions depend in each case on
!     - the derivatives of the shape functions with respect to the coordinates in
!       unit space
!     - the inverse Jacobi-Matrix of the mapping to the physical space
!     - OPEN: possibly the coodinates themselves, e.g. for axi-symmetric elements
!
!====================================================================================
!====================================================================================
MODULE BMatricesScalar
     USE ABQINTERFACE
     IMPLICIT NONE
     PRIVATE
     PUBLIC :: BMatScal, PreFactorDummy
     CONTAINS
! BMatScal function provides B-matrix for scalar-valued variable (e.g. temperature) which is equal to the gradient of the shape functions
          PURE FUNCTION BMatScal(JacobiInv, ShapeFuncDeriv, ShapeFunc, radius)
               REAL(KIND=AbqRK), INTENT(IN) :: JacobiInv(:,:), ShapeFuncDeriv(:,:)
               REAL(KIND=AbqRK), INTENT(IN), OPTIONAL :: ShapeFunc(:), radius
               REAL(KIND=AbqRK) :: BMatScal(SIZE(JacobiInv,1),SIZE(ShapeFuncDeriv,2))
               BMatScal = MATMUL(JacobiInv, ShapeFuncDeriv)
          END FUNCTION BMatScal
          PURE FUNCTION PreFactorDummy(radius)
              IMPLICIT NONE
              REAL(KIND=AbqRK), INTENT(IN), OPTIONAL :: radius
              REAL(KIND=AbqRK) :: PreFactorDummy
              PreFactorDummy = 1.0_AbqRK
          END FUNCTION PreFactorDummy
END MODULE BMatricesScalar
!-------------------------------------------------------------------------------------
! 1D Spherisymmetric
!-------------------------------------------------------------------------------------
MODULE BMatrices1DSpherisymmetric
     USE ABQINTERFACE
     USE BMatricesScalar
     IMPLICIT NONE
     PRIVATE
     INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NDI = 3, NSHR = 0, NTENS=NDI+NSHR
     LOGICAL, PARAMETER, PUBLIC :: SYM = .TRUE.
     PUBLIC :: BMatTensSym, PreFactor
     CONTAINS
          ! Definitions including computation of shape function gradient
          PURE FUNCTION BMatTensSym(JacobiInv, ShapeFuncDeriv, ShapeFunc, radius)
               REAL(KIND=AbqRK), INTENT(IN) :: JacobiInv(:,:), ShapeFuncDeriv(:,:)
               REAL(KIND=AbqRK), INTENT(IN), OPTIONAL :: ShapeFunc(:), radius
               REAL(KIND=AbqRK):: BMatTensSym(NTENS,SIZE(ShapeFuncDeriv,2)), BMatTens1D(SIZE(JacobiInv,1),SIZE(ShapeFuncDeriv,2))
               INTEGER(KIND=AbqIK) :: i, j, NNODES
               ! get 1D B-Matrix, first coordinate 
               BMatTens1D(:,:)=BMatScal(JacobiInv, ShapeFuncDeriv)
               NNODES = SIZE(ShapeFuncDeriv,2)
               ! eps_rr = u_r,r
               FORALL (i=1:NNODES) BMatTensSym(1,(i-1)*(NDI-2)+1) = BMatTens1D(1,i)
               ! eps_phiphi = eps_thetatheta = u_r / r
               FORALL (i=1:NNODES, j=2:3) BMatTensSym(j,(i-1)*(NDI-2)+1) = ShapeFunc(i)/radius
          END FUNCTION BMatTensSym
          PURE FUNCTION PreFactor(radius)
              IMPLICIT NONE
              REAL(KIND=AbqRK), INTENT(IN) :: radius
              REAL(KIND=AbqRK), PARAMETER :: pi=acos(-1.0_AbqRK)
              REAL(KIND=AbqRK) :: PreFactor
              ! spherisymmetric element: multiply determinant with 4*pi*r^2
              PreFactor = 4.0_AbqRK*pi*radius**2
          END FUNCTION PreFactor
END MODULE BMatrices1DSpherisymmetric
!-------------------------------------------------------------------------------------
! 2D Plane Stress
!-------------------------------------------------------------------------------------
MODULE BMatrices2DPlaneStress
     USE ABQINTERFACE
     USE BMatricesScalar, ONLY: BMatScal, PreFactor => PreFactorDummy
     IMPLICIT NONE
     PRIVATE
     INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NDI = 2, NSHR = 1, NTENS=NDI+NSHR
     LOGICAL, PARAMETER, PUBLIC :: SYM = .FALSE.
     PUBLIC :: BMatTensSym, BMatTensUnsym, BMatTensSym_SFG, BMatTensUnsym_SFG, PreFactor
     CONTAINS
! BMatTensSym function provides symmetric B-matrix for vector-valued variable (e.g. displacement) for 2D case
!          PURE FUNCTION BMatTensSym(JacobiInv, ShapeFuncDeriv)
!               REAL(KIND=AbqRK), INTENT(IN) :: JacobiInv(:,:), ShapeFuncDeriv(:,:)
!               REAL(KIND=AbqRK):: BMat_temp(SIZE(JacobiInv,1),SIZE(ShapeFuncDeriv,2))! temporary array to store result of 'BMatScal' function
!               REAL(KIND=AbqRK):: BMatTensSym(NDI+NSHR,NDI*SIZE(ShapeFuncDeriv,2))
!               INTEGER(KIND=AbqIK):: i, j, NNODES
!               !              | d/dx  0   |
!               ! [B_2D_sym] = |  0   d/dy | * [N] , where d/dx is a PARTIAL derivative & [N] is a shape function array
!               !              | d/dy d/dx |
!               NNODES = SIZE(ShapeFuncDeriv,2)
!               BMatTensSym = 0.0_AbqRK ! initialize with zeros
!               BMat_temp = BMatScal(JacobiInv, ShapeFuncDeriv)
!             ! substitute some zeros with elements of BMat_temp array
!               FORALL (i=1:NDI, j=1:NNODES) BMatTensSym(i,(j-1)*NDI+i) = BMat_temp(i,j)
!               FORALL (j=1:NNODES)
!                 BMatTensSym(3,j*NDI-1) = BMat_temp(2,j)
!                 BMatTensSym(3,j*NDI)   = BMat_temp(1,j)
!               END FORALL
!          END FUNCTION BMatTensSym

! Definitions with shape function gradient as argument
          PURE FUNCTION BMatTensSym_SFG(ShapeFuncGrad)
               REAL(KIND=AbqRK), INTENT(IN) :: ShapeFuncGrad(:,:)
               REAL(KIND=AbqRK) :: BMatTensSym_SFG(NTENS,2*SIZE(ShapeFuncGrad,2))
               INTEGER(KIND=AbqIK) :: i, j, NNODES
               !              | d/dx  0   |
               ! [B_2D_sym] = |  0   d/dy | * [N] , where d/dx is a PARTIAL derivative & [N] is a shape function array
               !              | d/dy d/dx |
               NNODES = SIZE(ShapeFuncGrad,2)
               BMatTensSym_SFG = 0.0_AbqRK ! initialize with zeros
               ! substitute some zeros with elements of ShapeFuncGrad
               FORALL (i=1:NDI, j=1:NNODES) BMatTensSym_SFG(i,(j-1)*NDI+i) = ShapeFuncGrad(i,j)
               FORALL (j=1:NNODES)
                 BMatTensSym_SFG(3,j*NDI-1) = ShapeFuncGrad(2,j)
                 BMatTensSym_SFG(3,j*NDI)   = ShapeFuncGrad(1,j)
               END FORALL
          END FUNCTION
! BMatTensUnsym function provides unsymmetric B-matrix for vector-valued variable (e.g. displacement) for 2D case
          PURE FUNCTION BMatTensUnsym_SFG(ShapeFuncGrad)
               REAL(KIND=AbqRK), INTENT(IN) :: ShapeFuncGrad(:,:)
               REAL(KIND=AbqRK):: BMatTensUnsym_SFG(4,2*SIZE(ShapeFuncGrad,2))
               INTEGER(KIND=AbqIK):: i, j, NNODES
               !                | d/dx   0   |
               !                |  0    d/dy |
               ! [B_2D_unsym] = | d/dy   0   | * [N] , where d/dx is a PARTIAL derivative & [N] is a shape function array
               !                |  0    d/dx |
               NNODES = SIZE(ShapeFuncGrad,2)
               BMatTensUnsym_SFG = 0.0_AbqRK ! initialize with zeros
               ! substitute some zeros with elements of BMat_temp array
               FORALL (i=1:NDI, j=1:NNODES) BMatTensUnsym_SFG(i,(j-1)*NDI+i) = ShapeFuncGrad(i,j)
               FORALL (j=1:NNODES)
                 BMatTensUnsym_SFG(3,j*NDI-1) = ShapeFuncGrad(2,j)
                 BMatTensUnsym_SFG(4,j*NDI)   = ShapeFuncGrad(1,j)
               END FORALL
          END FUNCTION
! Definitions including computation of shape function gradien
          PURE FUNCTION BMatTensSym(JacobiInv, ShapeFuncDeriv, ShapeFunc, radius)
               REAL(KIND=AbqRK), INTENT(IN) :: JacobiInv(:,:), ShapeFuncDeriv(:,:)
               REAL(KIND=AbqRK), INTENT(IN), OPTIONAL :: ShapeFunc(:), radius
               REAL(KIND=AbqRK) :: BMatTensSym(NTENS,2*SIZE(ShapeFuncDeriv,2))
               BMatTensSym=BMatTensSym_SFG(BMatScal(JacobiInv, ShapeFuncDeriv))
          END FUNCTION BMatTensSym
          PURE FUNCTION BMatTensUnsym(JacobiInv, ShapeFuncDeriv, ShapeFunc, radius)
               REAL(KIND=AbqRK), INTENT(IN) :: JacobiInv(:,:), ShapeFuncDeriv(:,:)
               REAL(KIND=AbqRK), INTENT(IN), OPTIONAL :: ShapeFunc(:), radius
               REAL(KIND=AbqRK):: BMatTensUnsym(4,2*SIZE(ShapeFuncDeriv,2))
               BMatTensUnsym=BMatTensUnsym_SFG(BMatScal(JacobiInv, ShapeFuncDeriv))
          END FUNCTION BMatTensUnsym
END MODULE BMatrices2DPlaneStress
!-------------------------------------------------------------------------------------
! 2D Plane Strain
!-------------------------------------------------------------------------------------
MODULE BMatrices2DPlaneStrain
     USE ABQINTERFACE
     USE BMatricesScalar, ONLY: BMatScal, PreFactor => PreFactorDummy
     IMPLICIT NONE
     PRIVATE
     INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NDI = 3, NSHR = 1, NTENS=NDI+NSHR
     LOGICAL, PARAMETER, PUBLIC :: SYM = .FALSE.
     PUBLIC :: BMatTensSym, BMatTensUnsym, BMatTensSym_SFG, BMatTensUnsym_SFG, PreFactor
     CONTAINS
! BMatTensSym function provides symmetric B-matrix for vector-valued variable (e.g. displacement) for 2D case
          PURE FUNCTION BMatTensSym_SFG(ShapeFuncGrad)
               REAL(KIND=AbqRK), INTENT(IN) :: ShapeFuncGrad(:,:)
               REAL(KIND=AbqRK):: BMatTensSym_SFG(NTENS,2*SIZE(ShapeFuncGrad,2))
               INTEGER(KIND=AbqIK):: i, j, NNODES
               !              | d/dx  0   |
               ! [B_2D_sym] = |  0   d/dy | * [N] , where d/dx is a PARTIAL derivative & [N] is a shape function array
               !              |  0    0   |
               !              | d/dy d/dx |
               NNODES = SIZE(ShapeFuncGrad,2)
               BMatTensSym_SFG = 0.0_AbqRK ! initialize with zeros
               ! substitute some zeros with elements of ShapeFuncGrad array
               FORALL (i=1:NDI-1, j=1:NNODES) BMatTensSym_SFG(i,(j-1)*(NDI-1)+i) = ShapeFuncGrad(i,j)
               FORALL (j=1:NNODES)
                 BMatTensSym_SFG(4,j*(NDI-1)-1) = ShapeFuncGrad(2,j)
                 BMatTensSym_SFG(4,j*(NDI-1))   = ShapeFuncGrad(1,j)
               END FORALL
          END FUNCTION
! BMatTensUnsym function provides unsymmetric B-matrix for vector-valued variable (e.g. displacement) for 2D case
          PURE FUNCTION BMatTensUnsym_SFG(ShapeFuncGrad)
               REAL(KIND=AbqRK), INTENT(IN) :: ShapeFuncGrad(:,:)
               REAL(KIND=AbqRK) :: BMatTensUnsym_SFG(5,2*SIZE(ShapeFuncGrad,2))
               INTEGER(KIND=AbqIK) :: i, j, NNODES
               !                | d/dx  0   |
               !                |  0   d/dy |
               ! [B_2D_unsym] = |  0    0   | * [N] , where d/dx is a PARTIAL derivative & [N] is a shape function array
               !                | d/dy  0   |
               !                |  0   d/dx |
               NNODES = SIZE(ShapeFuncGrad,2)
               BMatTensUnsym_SFG = 0.0_AbqRK ! initialize with zeros
               ! substitute some zeros with elements of ShapeFuncGrad array
               FORALL (i=1:NDI-1, j=1:NNODES) BMatTensUnsym_SFG(i,(j-1)*(NDI-1)+i) = ShapeFuncGrad(i,j)
               FORALL (j=1:NNODES)
                 BMatTensUnsym_SFG(4,j*(NDI-1)-1) = ShapeFuncGrad(2,j)
                 BMatTensUnsym_SFG(5,j*(NDI-1))   = ShapeFuncGrad(1,j)
               END FORALL
          END FUNCTION
! Definitions including computation of shape function gradien
          PURE FUNCTION BMatTensSym(JacobiInv, ShapeFuncDeriv, ShapeFunc, radius)
               REAL(KIND=AbqRK), INTENT(IN) :: JacobiInv(:,:), ShapeFuncDeriv(:,:)
               REAL(KIND=AbqRK), INTENT(IN), OPTIONAL :: ShapeFunc(:), radius
               REAL(KIND=AbqRK) :: BMatTensSym(NTENS,2*SIZE(ShapeFuncDeriv,2))
               BMatTensSym=BMatTensSym_SFG(BMatScal(JacobiInv, ShapeFuncDeriv))
          END FUNCTION BMatTensSym
          PURE FUNCTION BMatTensUnsym(JacobiInv, ShapeFuncDeriv, ShapeFunc, radius)
               REAL(KIND=AbqRK), INTENT(IN) :: JacobiInv(:,:), ShapeFuncDeriv(:,:)
               REAL(KIND=AbqRK), INTENT(IN), OPTIONAL :: ShapeFunc(:), radius
               REAL(KIND=AbqRK) :: BMatTensUnsym(5,2*SIZE(ShapeFuncDeriv,2))
               BMatTensUnsym=BMatTensUnsym_SFG(BMatScal(JacobiInv, ShapeFuncDeriv))
          END FUNCTION BMatTensUnsym
END MODULE BMatrices2DPlaneStrain
!-------------------------------------------------------------------------------------
! 2D Axisymmetric
!-------------------------------------------------------------------------------------
MODULE BMatrices2DAxisymmetric
     USE ABQINTERFACE
     USE BMatricesScalar
     USE BMatrices2DPlaneStrain, ONLY:  BMatTensSym_SFG_PlaneStrain => BMatTensSym_SFG, BMatTensUnsym_SFG_PlaneStrain => BMatTensUnsym_SFG
     IMPLICIT NONE
     PRIVATE
     INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NDI = 3, NSHR = 1, NTENS=NDI+NSHR
     LOGICAL, PARAMETER, PUBLIC :: SYM = .TRUE.
     PUBLIC :: BMatTensSym, BMatTensUnsym, PreFactor
     CONTAINS
          ! Definitions including computation of shape function gradient
          PURE FUNCTION BMatTensSym(JacobiInv, ShapeFuncDeriv, ShapeFunc, radius)
               REAL(KIND=AbqRK), INTENT(IN) :: JacobiInv(:,:), ShapeFuncDeriv(:,:)
               REAL(KIND=AbqRK), INTENT(IN), OPTIONAL :: ShapeFunc(:), radius
               REAL(KIND=AbqRK) :: BMatTensSym(NTENS,2*SIZE(ShapeFuncDeriv,2))
               INTEGER(KIND=AbqIK) :: i, NNODES
               ! get plane strain B-Matrix 
               BMatTensSym=BMatTensSym_SFG_PlaneStrain(BMatScal(JacobiInv, ShapeFuncDeriv))
               ! modify third coordinate (eps_phiphi)
               NNODES = SIZE(ShapeFuncDeriv,2)
               FORALL (i=1:NNODES) BMatTensSym(3,(i-1)*(NDI-1)+1) = ShapeFunc(i)/radius
          END FUNCTION BMatTensSym
          PURE FUNCTION BMatTensUnsym(JacobiInv, ShapeFuncDeriv, ShapeFunc, radius)
               REAL(KIND=AbqRK), INTENT(IN) :: JacobiInv(:,:), ShapeFuncDeriv(:,:)
               REAL(KIND=AbqRK), INTENT(IN), OPTIONAL :: ShapeFunc(:), radius
               REAL(KIND=AbqRK) :: BMatTensUnsym(5,2*SIZE(ShapeFuncDeriv,2))
               INTEGER(KIND=AbqIK) :: i, NNODES
               ! get plane strain B-Matrix 
               BMatTensUnsym=BMatTensUnsym_SFG_PlaneStrain(BMatScal(JacobiInv, ShapeFuncDeriv))
               ! modify third coordinate (eps_phiphi)
               NNODES = SIZE(ShapeFuncDeriv,2)
               FORALL (i=1:NNODES) BMatTensUnsym(3,(i-1)*(NDI-1)+1) = ShapeFunc(i)/radius
          END FUNCTION BMatTensUnsym
          PURE FUNCTION PreFactor(radius)
               IMPLICIT NONE
               REAL(KIND=AbqRK), INTENT(IN) :: radius
               REAL(KIND=AbqRK), PARAMETER :: pi=acos(-1.0_AbqRK)
               REAL(KIND=AbqRK):: PreFactor
               ! axisymmetric element: multiply determinant with 2*pi*r
               PreFactor = 2.0_AbqRK*pi*radius
          END FUNCTION PreFactor
END MODULE BMatrices2DAxisymmetric
!-------------------------------------------------------------------------------------
! 3D
!-------------------------------------------------------------------------------------
MODULE BMatrices3D
     USE ABQINTERFACE
     USE BMatricesScalar, ONLY: BMatScal, PreFactor => PreFactorDummy
     IMPLICIT NONE
     PRIVATE
     INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NDI = 3, NSHR = 3, NTENS=NDI+NSHR
     LOGICAL, PARAMETER, PUBLIC :: SYM = .FALSE.
     PUBLIC :: BMatTensSym, BMatTensUnsym, BMatTensSym_SFG, BMatTensUnsym_SFG, PreFactor
     CONTAINS
! BMatTensSym function provides symmetric B-matrix for vector-valued variable (e.g. displacement) for 3D case
          PURE FUNCTION BMatTensSym_SFG(ShapeFuncGrad)
               REAL(KIND=AbqRK), INTENT(IN) :: ShapeFuncGrad(:,:)
               REAL(KIND=AbqRK) :: BMatTensSym_SFG(NTENS,3*SIZE(ShapeFuncGrad,2))
               INTEGER(KIND=AbqIK) :: i, j, NNODES
               !              | d/dx   0     0   |
               !              |  0    d/dy   0   |
               !              |  0     0    d/dz |
               ! [B_3D_sym] = | d/dy  d/dx   0   | * [N] , where d/dx is a PARTIAL derivative & [N] is a shape function array
               !              | d/dz   0    d/dx |
               !              |  0    d/dz  d/dy |
               NNODES=SIZE(ShapeFuncGrad,2)
               BMatTensSym_SFG = 0_AbqRK ! initialize with zeros
               ! substitute some zeros with elements of ShapeFuncGrad array
               FORALL (i=1:NDI, j=1:NNODES) BMatTensSym_SFG(i,(j-1)*NDI+i) = ShapeFuncGrad(i,j)
               FORALL (j=1:NNODES)
                 BMatTensSym_SFG(4,j*NDI-2) = ShapeFuncGrad(2,j)
                 BMatTensSym_SFG(4,j*NDI-1) = ShapeFuncGrad(1,j)
                 BMatTensSym_SFG(5,j*NDI-2) = ShapeFuncGrad(3,j)
                 BMatTensSym_SFG(5,j*NDI)   = ShapeFuncGrad(1,j)
                 BMatTensSym_SFG(6,j*NDI-1) = ShapeFuncGrad(3,j)
                 BMatTensSym_SFG(6,j*NDI)   = ShapeFuncGrad(2,j)
              END FORALL
          END FUNCTION
! BMatTensUnsym function provides unsymmetric B-matrix for vector-valued variable (e.g. displacement) for 3D case
          PURE FUNCTION BMatTensUnsym_SFG(ShapeFuncGrad)
               REAL(KIND=AbqRK), INTENT(IN) :: ShapeFuncGrad(:,:)
               REAL(KIND=AbqRK) :: BMatTensUnsym_SFG(9,3*SIZE(ShapeFuncGrad,2))
               INTEGER(KIND=AbqIK) :: i, j, NNODES
               !                | d/dx   0     0   |         | xx |
               !                |  0    d/dy   0   |         | yy |
               !                |  0     0    d/dz |         | zz |
               !                | d/dy   0     0   |         | xy |
               ! [B_3D_unsym] = |  0    d/dx   0   | * [N] = | yx | , where d/dx is a PARTIAL derivative & [N] is a shape function array
               !                |  0    d/dz   0   |         | yz |
               !                |  0     0    d/dy |         | zy |
               !                | d/dz   0     0   |         | xz |
               !                |  0     0    d/dx |         | zx |
               NNODES = SIZE(ShapeFuncGrad,2)
               BMatTensUnsym_SFG = 0_AbqRK ! initialize with zeros
               ! substitute some zeros with elements of ShapeFuncGrad array
               FORALL (i=1:NDI, j=1:NNODES) BMatTensUnsym_SFG(i,(j-1)*NDI+i) = ShapeFuncGrad(i,j)
               FORALL (j=1:NNODES)
                 BMatTensUnsym_SFG(4,j*NDI-2) = ShapeFuncGrad(2,j)
                 BMatTensUnsym_SFG(5,j*NDI-1) = ShapeFuncGrad(1,j)
                 BMatTensUnsym_SFG(6,j*NDI-1) = ShapeFuncGrad(3,j)
                 BMatTensUnsym_SFG(7,j*NDI)   = ShapeFuncGrad(2,j)
                 BMatTensUnsym_SFG(8,j*NDI-2) = ShapeFuncGrad(3,j)
                 BMatTensUnsym_SFG(9,j*NDI)   = ShapeFuncGrad(1,j)
               END FORALL
          END FUNCTION
! Definitions including computation of shape function gradien
          PURE FUNCTION BMatTensSym(JacobiInv, ShapeFuncDeriv, ShapeFunc, radius)
               REAL(KIND=AbqRK), INTENT(IN) :: JacobiInv(:,:), ShapeFuncDeriv(:,:)
               REAL(KIND=AbqRK), INTENT(IN), OPTIONAL :: ShapeFunc(:), radius
               REAL(KIND=AbqRK) :: BMatTensSym(NTENS,3*SIZE(ShapeFuncDeriv,2))
               BMatTensSym=BMatTensSym_SFG(BMatScal(JacobiInv, ShapeFuncDeriv))
          END FUNCTION BMatTensSym
          PURE FUNCTION BMatTensUnsym(JacobiInv, ShapeFuncDeriv, ShapeFunc, radius)
               REAL(KIND=AbqRK), INTENT(IN) :: JacobiInv(:,:), ShapeFuncDeriv(:,:)
               REAL(KIND=AbqRK), INTENT(IN), OPTIONAL :: ShapeFunc(:), radius
               REAL(KIND=AbqRK) :: BMatTensUnsym(9,3*SIZE(ShapeFuncDeriv,2))
               BMatTensUnsym=BMatTensUnsym_SFG(BMatScal(JacobiInv, ShapeFuncDeriv))
          END FUNCTION BMatTensUnsym
END MODULE BMatrices3D

!====================================================================================
!====================================================================================
!
! SHAPE FUNCTIONS: module provides
! - shape functions (N-matrices) in terms of the unit coordinate vector chi
! - corresponding number of nodes (NNODES) and dimension (NDIM)
! - derivative matrix of shape functions with respect to chi
!
!====================================================================================
!====================================================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1D line, linear ansatz functions in x direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE ShapeFunc1D_Line_Lin
!     2 NODES, NODE NUMBERING:
!     1 2

  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ShapeFunc, ShapeFuncDeriv
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NDIM=1, NNODES=2
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: NodePos(NNODES,NDIM)=reshape( (/ -1.0_AbqRK, 1.0_AbqRK /),(/ NNODES, NDIM /) )
  !                                                                        1          2
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: ShapeFuncDerivConst(NDIM,NNODES)=reshape((/-0.5_AbqRK, 0.5_AbqRK/),(/NDIM,NNODES/))
  CONTAINS
!------------------------------------------------------------------------------------
    PURE FUNCTION ShapeFunc(chi)
       REAL(KIND=AbqRK),DIMENSION(NNODES)::ShapeFunc
       REAL(KIND=AbqRK), INTENT(IN), DIMENSION(NDIM) :: chi
       ShapeFunc(:)=(/1.0_AbqRK-chi(1), 1.0_AbqRK+chi(1)/)/2.0_AbqRK
    END FUNCTION ShapeFunc
!------------------------------------------------------------------------------------
    PURE FUNCTION ShapeFuncDeriv(chi)
       REAL(KIND=AbqRK),DIMENSION(NDIM,NNODES)::ShapeFuncDeriv
       REAL(KIND=AbqRK), INTENT(IN), DIMENSION(NDIM) :: chi
       ShapeFuncDeriv=ShapeFuncDerivConst
    END FUNCTION ShapeFuncDeriv
!------------------------------------------------------------------------------------
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1D line, quadratic ansatz functions in x direction; node order: ABAQUS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE ShapeFunc1D_Line_Quad_ABA
!     3 NODES, NODE NUMBERING:
!     1 2 3
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ShapeFunc, ShapeFuncDeriv
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NDIM=1, NNODES=3
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: NodePos(NNODES,NDIM)=reshape( (/ -1.0_AbqRK,  0.0_AbqRK,  1.0_AbqRK /),(/ NNODES, NDIM /) )
  !                                                                        1           2           3
  CONTAINS
!------------------------------------------------------------------------------------
    PURE FUNCTION ShapeFunc(chi)
      REAL(KIND=AbqRK), DIMENSION(NNODES) :: ShapeFunc
      REAL(KIND=AbqRK), INTENT(IN), DIMENSION(NDIM) :: chi
      REAL(KIND=AbqRK) :: chi2
      chi2=chi(1)**2
      ShapeFunc(:)=(/(chi2-chi(1))/2.0_AbqRK, 1.0_AbqRK-chi2, (chi2+chi(1))/2.0_AbqRK /)
    END FUNCTION ShapeFunc
!------------------------------------------------------------------------------------
    PURE FUNCTION ShapeFuncDeriv(chi)
      REAL(KIND=AbqRK), DIMENSION(NDIM,NNODES) :: ShapeFuncDeriv
      REAL(KIND=AbqRK), INTENT(IN), DIMENSION(NDIM) :: chi
      INTEGER(KIND=AbqIK):: i
      ShapeFuncDeriv(1,:)=(/chi(1)-0.5_AbqRK, -2.0_AbqRK*chi(1), chi(1)+0.5_AbqRK/)
    END FUNCTION ShapeFuncDeriv
!------------------------------------------------------------------------------------
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1D line, quadratic ansatz functions in x direction
! see COOK, R.D.: Concepts and Applications of Finite Element Analysis, p. 203
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE ShapeFunc1D_Line_Quad
!     3 NODES, NODE NUMBERING:
!     1 3 2
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ShapeFunc, ShapeFuncDeriv
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NDIM=1, NNODES=3
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: NodePos(NNODES,NDIM)=reshape( (/ -1.0_AbqRK,  1.0_AbqRK,  0.0_AbqRK /),(/ NNODES, NDIM /) )
  !                                                                        1           2           3
  CONTAINS
!------------------------------------------------------------------------------------
    PURE FUNCTION ShapeFunc(chi)
      REAL(KIND=AbqRK), DIMENSION(NNODES) :: ShapeFunc
      REAL(KIND=AbqRK), INTENT(IN), DIMENSION(NDIM) :: chi
      REAL(KIND=AbqRK) :: chi2
      chi2=chi(1)**2
      ShapeFunc(:)=(/(chi2-chi(1))/2.0_AbqRK, (chi2+chi(1))/2.0_AbqRK, 1.0_AbqRK-chi2 /)
    END FUNCTION ShapeFunc
!------------------------------------------------------------------------------------
    PURE FUNCTION ShapeFuncDeriv(chi)
      REAL(KIND=AbqRK), DIMENSION(NDIM,NNODES) :: ShapeFuncDeriv
      REAL(KIND=AbqRK), INTENT(IN), DIMENSION(NDIM) :: chi
      INTEGER(KIND=AbqIK):: i
      ShapeFuncDeriv(1,:)=(/chi(1)-0.5_AbqRK, chi(1)+0.5_AbqRK, -2.0_AbqRK*chi(1)/)
    END FUNCTION ShapeFuncDeriv
!------------------------------------------------------------------------------------
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1D line, cubic ansatz functions in x direction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE ShapeFunc1D_Line_Cubic
!     4 NODES, NODE NUMBERING:
!     1 2 3 4
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ShapeFunc, ShapeFuncDeriv
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NDIM=1, NNODES=4
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: NodePos(NNODES,NDIM)=reshape( (/ -1.0_AbqRK,  -1.0_AbqRK/3.0_AbqRK,  1.0_AbqRK/3.0_AbqRK,  1.0_AbqRK /),(/ NNODES, NDIM /) )
  !                                                                        1            2                     3                     4
  CONTAINS
!------------------------------------------------------------------------------------
    PURE FUNCTION ShapeFunc(chi)
      REAL(KIND=AbqRK), DIMENSION(NNODES) :: ShapeFunc
      REAL(KIND=AbqRK), INTENT(IN), DIMENSION(NDIM) :: chi
      REAL(KIND=AbqRK) :: chi1, chi2, chi3, chi4
      chi1 = -1.0_AbqRK
      chi2 = -1.0_AbqRK/3.0_AbqRK
      chi3 =  1.0_AbqRK/3.0_AbqRK
      chi4 =  1.0_AbqRK
      ShapeFunc(1) = (chi(1)-chi2)*(chi(1)-chi3)*(chi(1)-chi4)/(chi1-chi2)/(chi1-chi3)/(chi1-chi4)
      ShapeFunc(2) = (chi(1)-chi1)*(chi(1)-chi3)*(chi(1)-chi4)/(chi2-chi1)/(chi2-chi3)/(chi2-chi4)
      ShapeFunc(3) = (chi(1)-chi1)*(chi(1)-chi2)*(chi(1)-chi4)/(chi3-chi1)/(chi3-chi2)/(chi3-chi4)
      ShapeFunc(4) = (chi(1)-chi1)*(chi(1)-chi2)*(chi(1)-chi3)/(chi4-chi1)/(chi4-chi2)/(chi4-chi3)
    END FUNCTION ShapeFunc
!------------------------------------------------------------------------------------
    PURE FUNCTION ShapeFuncDeriv(chi)
      REAL(KIND=AbqRK), DIMENSION(NDIM,NNODES) :: ShapeFuncDeriv
      REAL(KIND=AbqRK), INTENT(IN), DIMENSION(NDIM) :: chi
      REAL(KIND=AbqRK) :: chi1, chi2, chi3, chi4
      chi1 = -1.0_AbqRK
      chi2 = -1.0_AbqRK/3.0_AbqRK
      chi3 =  1.0_AbqRK/3.0_AbqRK
      chi4 =  1.0_AbqRK
      ShapeFuncDeriv(1,1) = ( (chi(1)-chi2)*(chi(1)-chi3)+(chi(1)-chi3)*(chi(1)-chi4)+(chi(1)-chi2)*(chi(1)-chi4) )/(chi1-chi2)/(chi1-chi3)/(chi1-chi4)
      ShapeFuncDeriv(1,2) = ( (chi(1)-chi1)*(chi(1)-chi3)+(chi(1)-chi3)*(chi(1)-chi4)+(chi(1)-chi1)*(chi(1)-chi4) )/(chi2-chi1)/(chi2-chi3)/(chi2-chi4)
      ShapeFuncDeriv(1,3) = ( (chi(1)-chi1)*(chi(1)-chi2)+(chi(1)-chi2)*(chi(1)-chi4)+(chi(1)-chi1)*(chi(1)-chi4) )/(chi3-chi1)/(chi3-chi2)/(chi3-chi4)
      ShapeFuncDeriv(1,4) = ( (chi(1)-chi1)*(chi(1)-chi2)+(chi(1)-chi2)*(chi(1)-chi3)+(chi(1)-chi1)*(chi(1)-chi3) )/(chi4-chi1)/(chi4-chi2)/(chi4-chi3)
    END FUNCTION ShapeFuncDeriv
!------------------------------------------------------------------------------------
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 2D triangle, linear ansatz function
!! - the shape functions correspond directly to the area coordinates in unit space
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE ShapeFunc2D_Triangle_Lin
!     3 NODES, NODE NUMBERING:
!     3
!     1  2
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  PUBLIC ::ShapeFunc,ShapeFuncDeriv
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NDIM=2, NNODES=3
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: NodePos(NNODES,NDIM)=reshape((/0.0_AbqRK, 1.0_AbqRK, 0.0_AbqRK, &
                                                                        0.0_AbqRK, 0.0_AbqRK, 1.0_AbqRK/),(/NNODES,NDIM/))
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: ShapeFuncDerivConst(NDIM,NNODES)=reshape((/-1.0_AbqRK, -1.0_AbqRK, &
                                                                                     1.0_AbqRK,  0.0_AbqRK, &
                                                                                     0.0_AbqRK,  1.0_AbqRK/),(/NDIM,NNODES/))
  CONTAINS
!------------------------------------------------------------------------------------
    PURE FUNCTION ShapeFunc(chi)
       REAL(KIND=AbqRK),DIMENSION(NNODES)::ShapeFunc
       REAL(KIND=AbqRK), INTENT(IN), DIMENSION(NDIM) :: chi
       ShapeFunc(:)=(/1.0_AbqRK-chi(1)-chi(2), chi(1), chi(2)/)
    END FUNCTION ShapeFunc
!------------------------------------------------------------------------------------
    PURE FUNCTION ShapeFuncDeriv(chi)
       REAL(KIND=AbqRK),DIMENSION(NDIM,NNODES)::ShapeFuncDeriv
       REAL(KIND=AbqRK), INTENT(IN), DIMENSION(NDIM) :: chi
       ShapeFuncDeriv=ShapeFuncDerivConst
    END FUNCTION ShapeFuncDeriv
!------------------------------------------------------------------------------------
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D triangle, quadratic ansatz functions in x and y-direction
! by S. Roth
! see COOK, R.D.: Concepts and Applications of Finite Element Analysis
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE ShapeFunc2D_Triangle_Quad
!     6 NODES, NODE NUMBERING:
!     3
!     6  5
!     1  4  2

  USE ABQINTERFACE

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ShapeFunc, ShapeFuncDeriv
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NDIM=2, NNODES=6
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: NodePos(NNODES,NDIM)=reshape( (/ 0.0_AbqRK, 1.0_AbqRK, 0.0_AbqRK, 0.5_AbqRK, 0.5_AbqRK, 0.0_AbqRK, &
                                                                          0.0_AbqRK, 0.0_AbqRK, 1.0_AbqRK, 0.0_AbqRK, 0.5_AbqRK, 0.5_AbqRK /),(/ NNODES, NDIM /) )
  !                                                                       1          2          3          4          5          6
  CONTAINS
!------------------------------------------------------------------------------------
    PURE FUNCTION ShapeFunc(chi)
      REAL(KIND=AbqRK), DIMENSION(NNODES) :: ShapeFunc
      REAL(KIND=AbqRK), INTENT(IN), DIMENSION(NDIM) :: chi
      REAL(KIND=AbqRK) :: chis
      chis=chi(1)+chi(2)
      ShapeFunc(1:NNODES)=(/ (1.0_AbqRK-chis)*(1.0_AbqRK-2.0_AbqRK*chis), &
                             chi(1)*(2.0_AbqRK*chi(1)-1.0_AbqRK), &
                             chi(2)*(2.0_AbqRK*chi(2)-1.0_AbqRK), &
                             4.0_AbqRK*chi(1)*(1.0_AbqRK-chis), &
                             4.0_AbqRK*chi(1)*chi(2), &
                             4.0_AbqRK*chi(2)*(1.0_AbqRK-chis) /)
    END FUNCTION ShapeFunc
!------------------------------------------------------------------------------------
    PURE FUNCTION ShapeFuncDeriv(chi)
      REAL(KIND=AbqRK), DIMENSION(NDIM,NNODES) :: ShapeFuncDeriv
      REAL(KIND=AbqRK), INTENT(IN), DIMENSION(NDIM) :: chi
      REAL(KIND=AbqRK), DIMENSION(NDIM) :: chi4
      REAL(KIND=AbqRK) :: chis4
      chi4=4.0_AbqRK*chi
      chis4=chi4(1)+chi4(2)
      ShapeFuncDeriv(1,1:NNODES)=(/ chis4-3.0_AbqRK, chi4(1)-1.0_AbqRK, 0.0_AbqRK         , 4.0_AbqRK-chis4-chi4(1), chi4(2), -chi4(2) /)
      ShapeFuncDeriv(2,1:NNODES)=(/ chis4-3.0_AbqRK, 0.0_AbqRK        , chi4(2)-1.0_AbqRK , -chi4(1)               , chi4(1),  4.0_AbqRK-chis4-chi4(2)/)
    END FUNCTION ShapeFuncDeriv
!------------------------------------------------------------------------------------
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 2D quad, linear ansatz function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE ShapeFunc2D_Square_Lin
!     4 NODES, NODE NUMBERING:
!     4  3
!     1  2
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  PUBLIC ::ShapeFunc,ShapeFuncDeriv
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NDIM=2, NNODES=4
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: NodePos(NNODES,NDIM)=reshape((/-1.0_AbqRK,  1.0_AbqRK, 1.0_AbqRK, -1.0_AbqRK, &
                                                                        -1.0_AbqRK, -1.0_AbqRK, 1.0_AbqRK,  1.0_AbqRK /),(/NNODES,NDIM/))
  CONTAINS
!------------------------------------------------------------------------------------
    PURE FUNCTION ShapeFunc(chi)
       REAL(KIND=AbqRK), DIMENSION(NNODES) :: ShapeFunc
       REAL(KIND=AbqRK), INTENT(IN), DIMENSION(NDIM) :: chi
       REAL(KIND=AbqRK), DIMENSION(NNODES,NDIM) :: Chi0
       INTEGER(KIND=AbqIK):: i
       FORALL (i=1:NNODES) Chi0(i,:)=NodePos(i,:)*chi(:)
       ShapeFunc(:)=(1.0_AbqRK+Chi0(:,1))*(1.0_AbqRK+Chi0(:,2))/4.0_AbqRK
    END FUNCTION ShapeFunc
!------------------------------------------------------------------------------------
    PURE FUNCTION ShapeFuncDeriv(chi)
       REAL(KIND=AbqRK), DIMENSION(NDIM,NNODES) :: ShapeFuncDeriv
       REAL(KIND=AbqRK), INTENT(IN), DIMENSION(NDIM) :: chi
       REAL(KIND=AbqRK), DIMENSION(NNODES,NDIM) :: Chi0
       INTEGER(KIND=AbqIK):: i
       FORALL (i=1:NNODES) Chi0(i,:)=NodePos(i,:)*chi(:)
       FORALL (i=1:NNODES)
         ShapeFuncDeriv(1,i)=(1.0_AbqRK+Chi0(i,2))*NodePos(i,1)/4.0_AbqRK
         ShapeFuncDeriv(2,i)=(1.0_AbqRK+Chi0(i,1))*NodePos(i,2)/4.0_AbqRK
       END FORALL
    END FUNCTION ShapeFuncDeriv
!------------------------------------------------------------------------------------
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 2D quad, quadratic ansatz function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE ShapeFunc2D_Square_Quadr
!     8 NODES, NODE NUMBERING:
!     4 7 3
!     8   6
!     1 5 2
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  PUBLIC ::ShapeFunc,ShapeFuncDeriv
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NDIM=2,NNODES=8
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: NodePos(NNODES,NDIM)=reshape((/-1.0_AbqRK,  1.0_AbqRK, 1.0_AbqRK, -1.0_AbqRK,  0.0_AbqRK, 1.0_AbqRK, 0.0_AbqRK, -1.0_AbqRK, &
                                                                        -1.0_AbqRK, -1.0_AbqRK, 1.0_AbqRK,  1.0_AbqRK, -1.0_AbqRK, 0.0_AbqRK, 1.0_AbqRK,  0.0_AbqRK /),(/NNODES,NDIM/))
  CONTAINS
!------------------------------------------------------------------------------------
    PURE FUNCTION ShapeFunc(chi)
       REAL(KIND=AbqRK), DIMENSION(NNODES) :: ShapeFunc
       REAL(KIND=AbqRK), INTENT(IN), DIMENSION(NDIM) :: chi
       REAL(KIND=AbqRK), DIMENSION(NNODES,NDIM) :: Chi0
       INTEGER(KIND=AbqIK):: i
       FORALL (i=1:NNODES) Chi0(i,:)=NodePos(i,:)*chi(:)
       ! CORNER NODES
       FORALL (i=1:4) ShapeFunc(i)=(1.0_AbqRK+Chi0(i,1))*(1.0_AbqRK+Chi0(i,2))*(Chi0(i,1.0_AbqRK)+Chi0(i,2)-1.0_AbqRK)/4.0_AbqRK
       ! BOTTOM AND TOP MIDSIDE NODES
       FORALL (i=5:7:2) ShapeFunc(i)=(1.0_AbqRK-chi(1)**2)*(1.0_AbqRK+Chi0(i,2))/2.0_AbqRK
       ! RIGHT AND LEFT MIDSIDE NODES
       FORALL (i=6:8:2) ShapeFunc(i)=(1.0_AbqRK+Chi0(i,1))*(1.0_AbqRK-chi(2)**2)/2.0_AbqRK
    END FUNCTION ShapeFunc
!------------------------------------------------------------------------------------
    PURE FUNCTION ShapeFuncDeriv(chi)
       REAL(KIND=AbqRK), DIMENSION(NDIM,NNODES) :: ShapeFuncDeriv
       REAL(KIND=AbqRK), INTENT(IN), DIMENSION(NDIM) :: chi
       REAL(KIND=AbqRK), DIMENSION(NNODES,NDIM) :: Chi0
       INTEGER(KIND=AbqIK):: i
       FORALL (i=1:NNODES) Chi0(i,:)=NodePos(i,:)*chi(:)
       ! LOCAL DERIVATIVES OF SHAPE FUNCTIONS
       ! CORNER NODES
       FORALL (i=1:4)
         ShapeFuncDeriv(1,i)=(1.0_AbqRK+Chi0(i,2))*(Chi0(i,2)+2.0_AbqRK*Chi0(i,1))*NodePos(i,1)/4.0_AbqRK
         ShapeFuncDeriv(2,i)=(1.0_AbqRK+Chi0(i,1))*(Chi0(i,1)+2.0_AbqRK*Chi0(i,2))*NodePos(i,2)/4.0_AbqRK
       END FORALL
       ! BOTTOM AND TOP MIDSIDE NODES
       FORALL (i=5:7:2)
         ShapeFuncDeriv(1,i)=-chi(1)*(1.0_AbqRK+Chi0(i,2))
         ShapeFuncDeriv(2,i)=(1.0_AbqRK-chi(1)**2)*NodePos(i,2)/2.0_AbqRK
       END FORALL
       ! RIGHT AND LEFT MIDSIDE NODES
       FORALL (i=6:8:2)
         ShapeFuncDeriv(1,i)=(1.0_AbqRK-chi(2)**2)*NodePos(i,1)/2.0_AbqRK
         ShapeFuncDeriv(2,i)=-chi(2)*(1.0_AbqRK+Chi0(i,1))
       END FORALL
       ! END OF LOCAL DERIVATIVES OF SHAPE FUNCTIONS
    END FUNCTION ShapeFuncDeriv
!------------------------------------------------------------------------------------
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 2D quad, cubic ansatz function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE ShapeFunc2D_Square_Cubic
!     12 NODES, NODE NUMBERING:
!      4  11   7   3
!      8          10
!     12           6
!      1   5   9   2
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  PUBLIC ::ShapeFunc,ShapeFuncDeriv
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NDIM=2,NNODES=12
  REAL(KIND=AbqRK), PARAMETER :: pos=1.0_AbqRK/3.0_AbqRK
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: NodePos(NNODES,NDIM)=reshape((/-1.0_AbqRK,  1.0_AbqRK, 1.0_AbqRK, -1.0_AbqRK,       -pos, 1.0_AbqRK,       pos, -1.0_AbqRK,        pos, 1.0_AbqRK,      -pos, -1.0_AbqRK, &
                                                                        -1.0_AbqRK, -1.0_AbqRK, 1.0_AbqRK,  1.0_AbqRK, -1.0_AbqRK,      -pos, 1.0_AbqRK,        pos, -1.0_AbqRK,       pos, 1.0_AbqRK,       -pos /),(/NNODES,NDIM/))
  CONTAINS
!------------------------------------------------------------------------------------
    ! cubic shape functions
    PURE FUNCTION ShapeFunc(chi)
       REAL(KIND=AbqRK), DIMENSION(NNODES) :: ShapeFunc
       REAL(KIND=AbqRK), INTENT(IN), DIMENSION(NDIM) :: chi
       REAL(KIND=AbqRK), PARAMETER :: A=32.0_AbqRK, B=9.0_AbqRK/A
       REAL(KIND=AbqRK) :: C
       C=(9.0_AbqRK*chi(2)**2+9.0_AbqRK*chi(1)**2-10.0_AbqRK)/A
       ShapeFunc(1)=  (chi(1)-1.0_AbqRK)*(chi(2)-1.0_AbqRK)*C
       ShapeFunc(2)= -(1.0_AbqRK+chi(1))*(chi(2)-1.0_AbqRK)*C
       ShapeFunc(3)=  (1.0_AbqRK+chi(1))*(1.0_AbqRK+chi(2))*C
       ShapeFunc(4)= -(chi(1)-1.0_AbqRK)*(1.0_AbqRK+chi(2))*C
       ShapeFunc(5)= -B*(chi(1)-1.0_AbqRK)*(3.0_AbqRK*chi(1)-1.0_AbqRK)*(1.0_AbqRK+chi(1))*(chi(2)-1.0_AbqRK)
       ShapeFunc(6)=  B*(chi(2)-1.0_AbqRK)*(3.0_AbqRK*chi(2)-1.0_AbqRK)*(1.0_AbqRK+chi(2))*(1.0_AbqRK+chi(1))
       ShapeFunc(7)= -B*(chi(1)-1.0_AbqRK)*(1.0_AbqRK+3.0_AbqRK*chi(1))*(1.0_AbqRK+chi(1))*(1.0_AbqRK+chi(2))
       ShapeFunc(8)=  B*(chi(2)-1.0_AbqRK)*(1.0_AbqRK+3.0_AbqRK*chi(2))*(1.0_AbqRK+chi(2))*(chi(1)-1.0_AbqRK)
       ShapeFunc(9)=  B*(chi(1)-1.0_AbqRK)*(1.0_AbqRK+chi(1))*(chi(2)-1.0_AbqRK)*(1.0_AbqRK+3.0_AbqRK*chi(1))
       ShapeFunc(10)=-B*(1.0_AbqRK+chi(1))*(chi(2)-1.0_AbqRK)*(1.0_AbqRK+chi(2))*(1.0_AbqRK+3.0_AbqRK*chi(2))
       ShapeFunc(11)= B*(chi(1)-1.0_AbqRK)*(1.0_AbqRK+chi(1))*(1.0_AbqRK+chi(2))*(3.0_AbqRK*chi(1)-1.0_AbqRK)
       ShapeFunc(12)=-B*(chi(1)-1.0_AbqRK)*(chi(2)-1.0_AbqRK)*(1.0_AbqRK+chi(2))*(3.0_AbqRK*chi(2)-1.0_AbqRK)
    END FUNCTION ShapeFunc

!    ! linear shape functions
!    PURE FUNCTION ShapeFunc(chi)
!       REAL(KIND=AbqRK), DIMENSION(NNODES) :: ShapeFunc
!       REAL(KIND=AbqRK), INTENT(IN), DIMENSION(NDIM) :: chi
!       REAL(KIND=AbqRK), DIMENSION(4,NDIM) :: Chi0
!       INTEGER(KIND=AbqIK):: i
!       FORALL (i=1:4) Chi0(i,:)=NodePos(i,:)*chi(:)
!       ShapeFunc     = 0.0_AbqRK
!       ShapeFunc(1:4)=(1.0_AbqRK+Chi0(:,1))*(1.0_AbqRK+Chi0(:,2))/4.0_AbqRK
!    END FUNCTION ShapeFunc

!------------------------------------------------------------------------------------
    PURE FUNCTION ShapeFuncDeriv(chi)
       REAL(KIND=AbqRK), DIMENSION(NDIM,NNODES) :: ShapeFuncDeriv
       REAL(KIND=AbqRK), INTENT(IN), DIMENSION(NDIM) :: chi
       REAL(KIND=AbqRK), PARAMETER :: A=32.0_AbqRK, B=9.0_AbqRK/A

       ShapeFuncDeriv(1,1)=   (chi(2)-1.0_AbqRK)*(9.0_AbqRK*chi(2)**2+27.0_AbqRK*chi(1)**2-10.0_AbqRK-18.0_AbqRK*chi(1))/A
       ShapeFuncDeriv(1,2)=  -(chi(2)-1.0_AbqRK)*(9.0_AbqRK*chi(2)**2+27.0_AbqRK*chi(1)**2-10.0_AbqRK+18.0_AbqRK*chi(1))/A
       ShapeFuncDeriv(1,3)=   (1.0_AbqRK+chi(2))*(9.0_AbqRK*chi(2)**2+27.0_AbqRK*chi(1)**2-10.0_AbqRK+18.0_AbqRK*chi(1))/A
       ShapeFuncDeriv(1,4)=  -(1.0_AbqRK+chi(2))*(9.0_AbqRK*chi(2)**2+27.0_AbqRK*chi(1)**2-10.0_AbqRK-18.0_AbqRK*chi(1))/A
       ShapeFuncDeriv(1,5)=  -B*(-2.0_AbqRK*chi(1)-3.0_AbqRK+9.0_AbqRK*chi(1)**2)*(chi(2)-1.0_AbqRK)
       ShapeFuncDeriv(1,6)=   B*(chi(2)-1.0_AbqRK)*(3.0_AbqRK*chi(2)-1.0_AbqRK)*(1.0_AbqRK+chi(2))
       ShapeFuncDeriv(1,7)=  -B*(2.0_AbqRK*chi(1)-3.0_AbqRK+9.0_AbqRK*chi(1)**2)*(1.0_AbqRK+chi(2))
       ShapeFuncDeriv(1,8)=   B*(chi(2)-1.0_AbqRK)*(1.0_AbqRK+3.0_AbqRK*chi(2))*(1.0_AbqRK+chi(2))
       ShapeFuncDeriv(1,9)=   B*(2.0_AbqRK*chi(1)-3.0_AbqRK+9.0_AbqRK*chi(1)**2)*(chi(2)-1.0_AbqRK)
       ShapeFuncDeriv(1,10)= -B*(chi(2)-1.0_AbqRK)*(1.0_AbqRK+3.0_AbqRK*chi(2))*(1.0_AbqRK+chi(2))
       ShapeFuncDeriv(1,11)=  B*(-2.0_AbqRK*chi(1)-3.0_AbqRK+9.0_AbqRK*chi(1)**2)*(1.0_AbqRK+chi(2))
       ShapeFuncDeriv(1,12)= -B*(chi(2)-1.0_AbqRK)*(3.0_AbqRK*chi(2)-1.0_AbqRK)*(1.0_AbqRK+chi(2))

       ShapeFuncDeriv(2,1)=   (chi(1)-1.0_AbqRK)*(9.0_AbqRK*chi(1)**2+27.0_AbqRK*chi(2)**2-10.0_AbqRK-18.0_AbqRK*chi(2))/A
       ShapeFuncDeriv(2,2)=  -(1.0_AbqRK+chi(1))*(9.0_AbqRK*chi(1)**2+27.0_AbqRK*chi(2)**2-10.0_AbqRK-18.0_AbqRK*chi(2))/A
       ShapeFuncDeriv(2,3)=   (1.0_AbqRK+chi(1))*(9.0_AbqRK*chi(1)**2+27.0_AbqRK*chi(2)**2-10.0_AbqRK+18.0_AbqRK*chi(2))/A
       ShapeFuncDeriv(2,4)=  -(chi(1)-1.0_AbqRK)*(9.0_AbqRK*chi(1)**2+27.0_AbqRK*chi(2)**2-10.0_AbqRK+18.0_AbqRK*chi(2))/A
       ShapeFuncDeriv(2,5)=  -B*(chi(1)-1.0_AbqRK)*(3.0_AbqRK*chi(1)-1.0_AbqRK)*(1.0_AbqRK+chi(1))
       ShapeFuncDeriv(2,6)=   B*(-2.0_AbqRK*chi(2)-3.0_AbqRK+9.0_AbqRK*chi(2)**2)*(1.0_AbqRK+chi(1))
       ShapeFuncDeriv(2,7)=  -B*(chi(1)-1.0_AbqRK)*(1.0_AbqRK+3.0_AbqRK*chi(1))*(1.0_AbqRK+chi(1))
       ShapeFuncDeriv(2,8)=   B*(2.0_AbqRK*chi(2)-3.0_AbqRK+9.0_AbqRK*chi(2)**2)*(chi(1)-1.0_AbqRK)
       ShapeFuncDeriv(2,9)=   B*(chi(1)-1.0_AbqRK)*(1.0_AbqRK+3.0_AbqRK*chi(1))*(1.0_AbqRK+chi(1))
       ShapeFuncDeriv(2,10)= -B*(2.0_AbqRK*chi(2)-3.0_AbqRK+9.0_AbqRK*chi(2)**2)*(1.0_AbqRK+chi(1))
       ShapeFuncDeriv(2,11)=  B*(chi(1)-1.0_AbqRK)*(3.0_AbqRK*chi(1)-1.0_AbqRK)*(1.0_AbqRK+chi(1))
       ShapeFuncDeriv(2,12)= -B*(-2.0_AbqRK*chi(2)-3.0_AbqRK+9.0_AbqRK*chi(2)**2)*(chi(1)-1.0_AbqRK)

    END FUNCTION ShapeFuncDeriv
!------------------------------------------------------------------------------------
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 3D tetrahedron, linear ansatz function.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE ShapeFunc3D_Tetrahed_Lin
! 4 nodes, node numbering:
! Looking from the top in the negative Z direction:
! ^Y   (Y points to the top)
! |->X (X points to the right)
!
! @Z=1: 4 (X=0, Y=0)
! @Z=0: 3
!       1 2
    USE ABQINTERFACE
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: ShapeFunc, ShapeFuncDeriv
    INTEGER(KIND=AbqIK), PARAMETER, PUBLIC ::NDIM=3, NNODES=4
    REAL(KIND=AbqRK), PARAMETER, PUBLIC :: NodePos(NNODES,NDIM)=reshape((/ &
           0.0_AbqRK, 1.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, &
           0.0_AbqRK, 0.0_AbqRK, 1.0_AbqRK, 0.0_AbqRK, &
           0.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, 1.0_AbqRK/),(/NNODES,NDIM/))
    REAL(KIND=AbqRK), PARAMETER, PUBLIC :: ShapeFuncDerivConst(NDIM,NNODES)=reshape((/ &
          -1.0_AbqRK, -1.0_AbqRK, -1.0_AbqRK, &
           1.0_AbqRK,  0.0_AbqRK,  0.0_AbqRK, &
           0.0_AbqRK,  1.0_AbqRK,  0.0_AbqRK, &
           0.0_AbqRK,  0.0_AbqRK,  1.0_AbqRK /),(/NDIM,NNODES/))
    CONTAINS
!------------------------------------------------------------------------------------
       PURE FUNCTION ShapeFunc(chi)
           REAL(KIND=AbqRK), DIMENSION(NNODES)::ShapeFunc
           REAL(KIND=AbqRK), INTENT(IN), DIMENSION(NDIM) :: chi
           ShapeFunc = (/1.0_AbqRK-chi(1)-chi(2)-chi(3), chi(1), chi(2), chi(3)/)
       END FUNCTION ShapeFunc
!------------------------------------------------------------------------------------
       PURE FUNCTION ShapeFuncDeriv(chi)
           REAL(KIND=AbqRK), DIMENSION(NDIM,NNODES):: ShapeFuncDeriv
           REAL(KIND=AbqRK), INTENT(IN), DIMENSION(NDIM):: chi
           ShapeFuncDeriv = ShapeFuncDerivConst
       END FUNCTION ShapeFuncDeriv
!------------------------------------------------------------------------------------
END MODULE ShapeFunc3D_Tetrahed_Lin

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 3D tetrahedron, quadratic ansatz function.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE ShapeFunc3D_Tetrahed_Quadr
! 10 nodes, node numbering:
! Looking from the top in the negative Z direction:
! ^Y   (Y points to the top)
! |->X (X points to the right)
!
! @Z=1:   4 (X=0, Y=0)
! @Z=0.5: 10
!         8 9
! @Z=0:   3
!         7  6
!         1 5  2
    USE ABQINTERFACE
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: ShapeFunc, ShapeFuncDeriv
    INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NDIM = 3, NNODES = 10
    REAL(KIND=AbqRK), PARAMETER, PUBLIC :: NodePos(NNODES,NDIM)=reshape((/&
           0.0_AbqRK, 1.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, 0.5_AbqRK, 0.5_AbqRK, 0.0_AbqRK, 0.0_AbqRK, 0.5_AbqRK, 0.0_AbqRK, &
           0.0_AbqRK, 0.0_AbqRK, 1.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, 0.5_AbqRK, 0.5_AbqRK, 0.0_AbqRK, 0.0_AbqRK, 0.5_AbqRK, &
           0.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, 1.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, 0.5_AbqRK, 0.5_AbqRK, 0.5_AbqRK /),(/NNODES,NDIM/))
    CONTAINS
!------------------------------------------------------------------------------------
       PURE FUNCTION ShapeFunc(chi)
           REAL(KIND=AbqRK), DIMENSION(NNODES):: ShapeFunc
           REAL(KIND=AbqRK), INTENT(IN), DIMENSION(NDIM):: chi
           REAL(KIND=AbqRK):: lambda
           lambda = 1.0_AbqRK - ( chi(1)+chi(2)+chi(3) )
           ShapeFunc(:)=(/ lambda*(2.0_AbqRK*lambda-1.0_AbqRK), chi(1)*(2.0_AbqRK*chi(1)-1.0_AbqRK), chi(2)*(2.0_AbqRK*chi(2)-1.0_AbqRK), &
                           chi(3)*(2.0_AbqRK*chi(3)-1.0_AbqRK), 4.0_AbqRK*chi(1)*lambda, 4.0_AbqRK*chi(1)*chi(2), 4.0_AbqRK*chi(2)*lambda, &
                           4.0_AbqRK*chi(3)*lambda, 4.0_AbqRK*chi(1)*chi(3), 4.0_AbqRK*chi(2)*chi(3) /)
       END FUNCTION ShapeFunc
!------------------------------------------------------------------------------------
       PURE FUNCTION ShapeFuncDeriv(chi)
           REAL(KIND=AbqRK), DIMENSION(NDIM,NNODES):: ShapeFuncDeriv
           REAL(KIND=AbqRK), INTENT(IN), DIMENSION(NDIM):: chi
           REAL(KIND=AbqRK):: lambda, kappa
           REAL(KIND=AbqRK), PARAMETER:: zero = 0.0_AbqRK
           lambda = 1.0_AbqRK - ( chi(1)+chi(2)+chi(3) )
           kappa = 1.0_AbqRK - 4.0_AbqRK*lambda
           ShapeFuncDeriv = reshape((/ kappa, kappa, kappa, &
               4.0_AbqRK*chi(1)-1.0_AbqRK, zero, zero, &
               zero, 4.0_AbqRK*chi(2)-1.0_AbqRK, zero, &
               zero, zero, 4.0_AbqRK*chi(3)-1.0_AbqRK, &
               4.0_AbqRK*(lambda-chi(1)), -4.0_AbqRK*chi(1), -4.0_AbqRK*chi(1), &
               4.0_AbqRK*chi(2), 4.0_AbqRK*chi(1), zero, &
               -4.0_AbqRK*chi(2), 4.0_AbqRK*(lambda-chi(2)), -4.0_AbqRK*chi(2), &
               -4.0_AbqRK*chi(3), -4.0_AbqRK*chi(3), 4.0_AbqRK*(lambda-chi(3)), &
               4.0_AbqRK*chi(3), zero, 4.0_AbqRK*chi(1), &
               zero, 4.0_AbqRK*chi(3), 4.0_AbqRK*chi(2) /),(/ NDIM,NNODES /))
       END FUNCTION ShapeFuncDeriv
!------------------------------------------------------------------------------------
END MODULE ShapeFunc3D_Tetrahed_Quadr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 3D hexagon, linear ansatz function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE ShapeFunc3D_Hex_Lin
!     8 NODES, NODE NUMBERING:
!          Looking from the top in the negative z direction:
!          ^y (y points to the top)
!          |->x (x points to the right)
!   @z=-1:  4  3
!           1  2
!   @z=1:   8  7
!           5  6
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  PUBLIC ::ShapeFunc,ShapeFuncDeriv
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NDIM=3,NNODES=8
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: NodePos(NNODES,NDIM)=reshape((/-1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK,-1.0_AbqRK,-1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK,-1.0_AbqRK, &
                                                                        -1.0_AbqRK,-1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK,-1.0_AbqRK,-1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, &
                                                                        -1.0_AbqRK,-1.0_AbqRK,-1.0_AbqRK,-1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK /),(/NNODES,NDIM/))
  CONTAINS
!------------------------------------------------------------------------------------
    PURE FUNCTION ShapeFunc(chi)
       REAL(KIND=AbqRK), DIMENSION(NNODES) :: ShapeFunc
       REAL(KIND=AbqRK), INTENT(IN), DIMENSION(NDIM) :: chi
       REAL(KIND=AbqRK), DIMENSION(NNODES,NDIM) :: Chi0
       INTEGER(KIND=AbqIK):: i
       FORALL (i=1:NNODES) Chi0(i,:)= NodePos(i,:)*chi(:)
       ShapeFunc(:)=(1.0_AbqRK+Chi0(:,1))*(1.0_AbqRK+Chi0(:,2))*(1.0_AbqRK+Chi0(:,3))/8.0_AbqRK
    END FUNCTION ShapeFunc
!------------------------------------------------------------------------------------
    PURE FUNCTION ShapeFuncDeriv(chi)
       REAL(KIND=AbqRK), DIMENSION(NDIM,NNODES) :: ShapeFuncDeriv
       REAL(KIND=AbqRK), INTENT(IN), DIMENSION(NDIM) :: chi
       REAL(KIND=AbqRK), DIMENSION(NNODES,NDIM) :: Chi0
       INTEGER(KIND=AbqIK):: i
       FORALL (i=1:NNODES)
         Chi0(i,:)= NodePos(i,:)*chi(:)
         ShapeFuncDeriv(1,i)=(1.0_AbqRK+Chi0(i,2))*(1.0_AbqRK+Chi0(i,3))*NodePos(i,1)/8.0_AbqRK
         ShapeFuncDeriv(2,i)=(1.0_AbqRK+Chi0(i,1))*(1.0_AbqRK+Chi0(i,3))*NodePos(i,2)/8.0_AbqRK
         ShapeFuncDeriv(3,i)=(1.0_AbqRK+Chi0(i,1))*(1.0_AbqRK+Chi0(i,2))*NodePos(i,3)/8.0_AbqRK
       END FORALL
    END FUNCTION ShapeFuncDeriv
!------------------------------------------------------------------------------------
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 3D hexagon, quadratic ansatz function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE ShapeFunc3D_Hex_Quadr
!     20 NODES, NODE NUMBERING:
!     Looking from the top in the negative z direction:
!
!     ^y (y points to the top)
!     |->x (x points to the right)
!     @z=-1:  4 11  3
!            12    10
!             1  9  2
!
!      @z=1:  8 15  7
!            16    14
!             5 13  6
!
!      @z=0: 20 19
!            17 18
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: ShapeFunc, ShapeFuncDeriv
  INTEGER, PARAMETER, PUBLIC :: NDIM=3, NNODES=20
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: NodePos(NNODES,NDIM)=reshape((/ -1.0_AbqRK,  1.0_AbqRK,  1.0_AbqRK, -1.0_AbqRK, -1.0_AbqRK,  1.0_AbqRK,  1.0_AbqRK, -1.0_AbqRK,  0.0_AbqRK,  1.0_AbqRK, &
                                                                          0.0_AbqRK, -1.0_AbqRK,  0.0_AbqRK,  1.0_AbqRK,  0.0_AbqRK, -1.0_AbqRK, -1.0_AbqRK,  1.0_AbqRK,  1.0_AbqRK, -1.0_AbqRK, &
                                                                         -1.0_AbqRK, -1.0_AbqRK,  1.0_AbqRK,  1.0_AbqRK, -1.0_AbqRK, -1.0_AbqRK,  1.0_AbqRK,  1.0_AbqRK, -1.0_AbqRK,  0.0_AbqRK, &
                                                                          1.0_AbqRK,  0.0_AbqRK, -1.0_AbqRK,  0.0_AbqRK,  1.0_AbqRK,  0.0_AbqRK, -1.0_AbqRK, -1.0_AbqRK,  1.0_AbqRK,  1.0_AbqRK, &
                                                                         -1.0_AbqRK, -1.0_AbqRK, -1.0_AbqRK, -1.0_AbqRK,  1.0_AbqRK,  1.0_AbqRK,  1.0_AbqRK,  1.0_AbqRK, -1.0_AbqRK, -1.0_AbqRK, &
                                                                         -1.0_AbqRK, -1.0_AbqRK,  1.0_AbqRK,  1.0_AbqRK,  1.0_AbqRK,  1.0_AbqRK,  0.0_AbqRK,  0.0_AbqRK,  0.0_AbqRK,  0.0_AbqRK /),(/NNODES,NDIM/))
  CONTAINS
!------------------------------------------------------------------------------------
    PURE FUNCTION ShapeFunc(chi)
       REAL(KIND=AbqRK), DIMENSION(NNODES) :: ShapeFunc
       REAL(KIND=AbqRK), INTENT(IN), DIMENSION(NDIM) :: chi
       REAL(KIND=AbqRK), DIMENSION(NNODES,NDIM) :: Chi0
       INTEGER(KIND=AbqIK):: i
       FORALL (i=1:NNODES) Chi0(i,:)=NodePos(i,:)*chi(:)
       ! Corner nodes (1-8)
       FORALL (i=1:8) ShapeFunc(i)=(1.0_AbqRK+Chi0(i,1))*(1.0_AbqRK+Chi0(i,2))*(1.0_AbqRK+Chi0(i,3))*(Chi0(i,1)+Chi0(i,2)+Chi0(i,3)-2.0_AbqRK)/8.0_AbqRK
       ! Bottom & top midside nodes with x=0
       FORALL (i=9:15:2) ShapeFunc(i)=(1.0_AbqRK-chi(1)**2)*(1.0_AbqRK+Chi0(i,2))*(1.0_AbqRK+Chi0(i,3))/4.0_AbqRK
       ! Bottom & top midside nodes with y=0
       FORALL (i=10:16:2) ShapeFunc(i)=(1.0_AbqRK-chi(2)**2)*(1.0_AbqRK+Chi0(i,1))*(1.0_AbqRK+Chi0(i,3))/4.0_AbqRK
       ! Midside nodes with z=0
       FORALL (i=17:20) ShapeFunc(i)=(1.0_AbqRK+Chi0(i,1))*(1.0_AbqRK+Chi0(i,2))*(1.0_AbqRK-chi(3)**2)/4.0_AbqRK
    END FUNCTION ShapeFunc
!------------------------------------------------------------------------------------
    PURE FUNCTION ShapeFuncDeriv(chi)
       REAL(KIND=AbqRK), DIMENSION(NDIM,NNODES) :: ShapeFuncDeriv
       REAL(KIND=AbqRK), INTENT(IN), DIMENSION(NDIM) :: chi
       REAL(KIND=AbqRK), DIMENSION(NNODES,NDIM) :: Chi0
       INTEGER(KIND=AbqIK):: i
       FORALL (i=1:NNODES) Chi0(i,:)=NodePos(i,:)*chi(:)
       ! LOCAL DERIVATIVES OF SHAPE FUNCTIONS
       FORALL (i=1:8)
         ! Corner nodes (1-8)
         ShapeFuncDeriv(1,i)=(1.0_AbqRK+Chi0(i,2))*(1.0_AbqRK+Chi0(i,3))*(Chi0(i,2)+Chi0(i,3)+2.0_AbqRK*Chi0(i,1)-1)*NodePos(i,1)/8.0_AbqRK
         ShapeFuncDeriv(2,i)=(1.0_AbqRK+Chi0(i,1))*(1.0_AbqRK+Chi0(i,3))*(Chi0(i,1)+Chi0(i,3)+2.0_AbqRK*Chi0(i,2)-1)*NodePos(i,2)/8.0_AbqRK
         ShapeFuncDeriv(3,i)=(1.0_AbqRK+Chi0(i,1))*(1.0_AbqRK+Chi0(i,2))*(Chi0(i,1)+Chi0(i,2)+2.0_AbqRK*Chi0(i,3)-1)*NodePos(i,3)/8.0_AbqRK
       END FORALL
       ! Bottom & top midside nodes with x=0
       FORALL (i=9:15:2)
         ShapeFuncDeriv(1,i)=-chi(1)*(1.0_AbqRK+Chi0(i,2))*(1.0_AbqRK+Chi0(i,3))/2.0_AbqRK
         ShapeFuncDeriv(2,i)=(1.0_AbqRK-chi(1)**2)*(1.0_AbqRK+Chi0(i,3))*NodePos(i,2)/4.0_AbqRK
         ShapeFuncDeriv(3,i)=(1.0_AbqRK-chi(1)**2)*(1.0_AbqRK+Chi0(i,2))*NodePos(i,3)/4.0_AbqRK
       END FORALL
       ! Bottom & top midside nodes with y=0
       FORALL (i=10:16:2)
         ShapeFuncDeriv(1,i)=(1.0_AbqRK-chi(2)**2)*(1.0_AbqRK+Chi0(i,3))*NodePos(i,1)/4.0_AbqRK
         ShapeFuncDeriv(2,i)=-chi(2)*(1.0_AbqRK+Chi0(i,1))*(1.0_AbqRK+Chi0(i,3))/2.0_AbqRK
         ShapeFuncDeriv(3,i)=(1.0_AbqRK-chi(2)**2)*(1.0_AbqRK+Chi0(i,1))*NodePos(i,3)/4.0_AbqRK
       END FORALL
       ! Midside nodes with z=0
       FORALL (i=17:20)
         ShapeFuncDeriv(1,i)=(1.0_AbqRK-chi(3)**2)*(1.0_AbqRK+Chi0(i,2))*NodePos(i,1)/4.0_AbqRK
         ShapeFuncDeriv(2,i)=(1.0_AbqRK-chi(3)**2)*(1.0_AbqRK+Chi0(i,1))*NodePos(i,2)/4.0_AbqRK
         ShapeFuncDeriv(3,i)=-chi(3)*(1.0_AbqRK+Chi0(i,1))*(1.0_AbqRK+Chi0(i,2))/2.0_AbqRK
       END FORALL
       ! END OF LOCAL DERIVATIVES OF SHAPE FUNCTIONS
    END FUNCTION ShapeFuncDeriv
!------------------------------------------------------------------------------------
END MODULE

!====================================================================================
!====================================================================================
!
! GAUSS-Quadratur: Module provides
! - NGP: Number of Gauss-points
! - NDIMGP: DIMENSION of unit space
! - GPWeight and GPPos: weights and positions of Gauss-points in unit space
!
!====================================================================================
!====================================================================================

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1D, line, 1 GP
! by S. Roth
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE Integr1D_Line1GP
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NGP=1, NDIMGP=1
  REAL(KIND=AbqRK), PARAMETER :: PGP1=0.0_AbqRK, WGP1=2.0_AbqRK
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPWeight(NGP)=(/ WGP1 /)
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPPos(NGP,NDIMGP)=reshape( (/ PGP1 /) , (/ NGP, NDIMGP /) )
  !                                                                    1
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1D, line, 2 GP
! by S. Roth
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE Integr1D_Line2GP
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NGP=2, NDIMGP=1
  REAL(KIND=AbqRK), PARAMETER :: PGP1=1.0_AbqRK/SQRT(3.0_AbqRK), WGP1=1.0_AbqRK
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPWeight(NGP)=(/ WGP1, WGP1 /)
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPPos(NGP,NDIMGP)=reshape( (/ -PGP1, PGP1 /) , (/ NGP, NDIMGP /) )
  !                                                                     1     2
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 1D, line, 3 GP
! by S. Roth
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE Integr1D_Line3GP
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NGP=3, NDIMGP=1
  REAL(KIND=AbqRK), PARAMETER :: PGP1=SQRT(0.6_AbqRK), PGP2=0.0_AbqRK, WGP1=5.0_AbqRK/9._AbqRK, WGP2=8.0_AbqRK/9._AbqRK
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPWeight(NGP)=(/ WGP1, WGP1, WGP2 /)
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPPos(NGP,NDIMGP)=reshape( (/ -PGP1, PGP1, PGP2 /) , (/ NGP, NDIMGP /) )
  !                                                                     1     2     3
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 2D, triangle, 1GP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE Integr2D_Triangle1GP
!     1 GAUSS Point
!     GAUSS INTEGRATION POINTS & WEIGHTS
!     Note: weight is set to 1/2 to obtain the correct area for an integrand of 1,
!           Zienkiewiecz gives a weigt of 1 scaling of the Jacobian determinant
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NGP=1, NDIMGP=2
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPWeight(NGP)=(/0.5_AbqRK/)
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPPos(NGP,NDIMGP)=1._AbqRK/3.0_AbqRK*reshape((/1.0_AbqRK, &
                                                                                        1.0_AbqRK /) ,(/NGP,NDIMGP/))
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D, triangle, 3 GP (by S. Roth)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE Integr2D_Triangle3GP
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NGP=3, NDIMGP=2
  REAL(KIND=AbqRK), PARAMETER :: PGP1=1.0_AbqRK/6.0_AbqRK, PGP2=4.0_AbqRK*PGP1, WGP1=2.0_AbqRK*PGP1
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPWeight(NGP)=(/ WGP1, WGP1, WGP1 /)
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPPos(NGP,NDIMGP)=reshape( (/ PGP1, PGP2, PGP1, &
                                                                       PGP1, PGP1, PGP2 /) , (/ NGP, NDIMGP /) )
  !                                                                    1     2     3
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D, triangle, 3 GP at the nodal positions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE Integr2D_Triangle3GPNodalPos
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  INTEGER(kind=AbqIK), PARAMETER, PUBLIC :: NGP=3, NDIMGP=2
  REAL(kind=AbqRK), PARAMETER :: PGP1=0.0_AbqRK, PGP2=1.0_AbqRK, WGP1=1.0_AbqRK/3.0_AbqRK
  REAL(kind=AbqRK), PARAMETER, PUBLIC :: GPWeight(NGP)=(/ WGP1, WGP1, WGP1 /)
  REAL(kind=AbqRK), PARAMETER, PUBLIC :: GPPos(NGP,NDIMGP)=reshape( (/ PGP1, PGP2, PGP1, &
                                                                       PGP1, PGP1, PGP2 /) , (/ NGP, NDIMGP /) )
  !                                                                    1     2     3
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2D, triangle, 4 GP (by S. Roth)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE Integr2D_Triangle4GP
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NGP=4, NDIMGP=2
  REAL(KIND=AbqRK), PARAMETER :: PGP1=0.2_AbqRK, PGP2=0.6_AbqRK, PGP3=1.0_AbqRK/3.0_AbqRK, WGP1=25.0_AbqRK/48.0_AbqRK, WGP2=-27.0_AbqRK/48.0_AbqRK
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPWeight(NGP)=(/ WGP1, WGP1, WGP1, WGP2 /)
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPPos(NGP,NDIMGP)=reshape( (/ PGP1, PGP2, PGP1, PGP3, &
                                                                       PGP1, PGP1, PGP2, PGP3 /) , (/ NGP, NDIMGP /) )
  !                                                                    1     2     3     4
END MODULE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 2D, square
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE Integr2D_Square36GP
!     36 Integration Points (Gauss product rule)
!     GAUSS INTEGRATION POINTS & WEIGHTS
!     36 POINTS, INTEGRATION POINT NUMBERING:
!     31 32 33 34 35 36
!     25 26 27 28 29 30
!     19 20 21 22 23 24
!     13 14 15 16 17 18
!      7  8  9 10 11 12
!      1  2  3  4  5  6
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NGP=36, NDIMGP=2
  REAL(KIND=AbqRK), PARAMETER :: WGP1=0.1713244924_AbqRK, WGP2=0.3607615731_AbqRK, WGP3=0.467913935_AbqRK, xi1=0.9324695142_AbqRK, xi2=0.6612093865_AbqRK, xi3=0.2386191861_AbqRK
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPWeight(NGP)=(/WGP1*WGP1,WGP2*WGP1,WGP3*WGP1,WGP3*WGP1,WGP2*WGP1,WGP1*WGP1, &
                                                         WGP1*WGP2,WGP2*WGP2,WGP3*WGP2,WGP3*WGP2,WGP2*WGP2,WGP1*WGP2, &
                                                         WGP1*WGP3,WGP2*WGP3,WGP3*WGP3,WGP3*WGP3,WGP2*WGP3,WGP1*WGP3, &
                                                         WGP1*WGP3,WGP2*WGP3,WGP3*WGP3,WGP3*WGP3,WGP2*WGP3,WGP1*WGP3, &
                                                         WGP1*WGP2,WGP2*WGP2,WGP3*WGP2,WGP3*WGP2,WGP2*WGP2,WGP1*WGP2, &
                                                         WGP1*WGP1,WGP2*WGP1,WGP3*WGP1,WGP3*WGP1,WGP2*WGP1,WGP1*WGP1  /)
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPPos(NGP,NDIMGP)=reshape((/ &
  -xi1, -xi2, -xi3,  xi3,  xi2,  xi1, -xi1, -xi2, -xi3,  xi3,  xi2,  xi1, -xi1, -xi2, -xi3,  xi3,  xi2,  xi1, -xi1, -xi2, -xi3,  xi3,  xi2,  xi1, -xi1, -xi2, -xi3,  xi3,  xi2,  xi1, -xi1, -xi2, -xi3,  xi3,  xi2,  xi1, &
  -xi1, -xi1, -xi1, -xi1, -xi1, -xi1, -xi2, -xi2, -xi2, -xi2, -xi2, -xi2, -xi3, -xi3, -xi3, -xi3, -xi3, -xi3,  xi3,  xi3,  xi3,  xi3,  xi3,  xi3,  xi2,  xi2,  xi2,  xi2,  xi2,  xi2,  xi1,  xi1,  xi1,  xi1,  xi1,  xi1 /) ,(/NGP,NDIMGP/))
END MODULE

MODULE Integr2D_Square25GP
!     25 Integration Points (Gauss product rule)
!     GAUSS INTEGRATION POINTS & WEIGHTS
!     25 POINTS, INTEGRATION POINT NUMBERING:
!     21 22 23 24 25
!     16 17 18 19 20
!     11 12 13 14 15
!      6  7  8  9 10
!      1  2  3  4  5
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NGP=25, NDIMGP=2
  REAL(KIND=AbqRK), PARAMETER :: WGP1=0.2369268851_AbqRK, WGP2=0.4786286705_AbqRK, WGP3=0.568888889_AbqRK, xi1=0.9061798459_AbqRK, xi2=0.5384693101_AbqRK, xi3=0.0_AbqRK
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPWeight(NGP)=(/WGP1*WGP1,WGP2*WGP1,WGP3*WGP1,WGP2*WGP1,WGP1*WGP1, &
                                                         WGP1*WGP2,WGP2*WGP2,WGP3*WGP2,WGP2*WGP2,WGP1*WGP2, &
                                                         WGP1*WGP3,WGP2*WGP3,WGP3*WGP3,WGP2*WGP3,WGP1*WGP3, &
                                                         WGP1*WGP2,WGP2*WGP2,WGP3*WGP2,WGP2*WGP2,WGP1*WGP2, &
                                                         WGP1*WGP1,WGP2*WGP1,WGP3*WGP1,WGP2*WGP1,WGP1*WGP1  /)
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPPos(NGP,NDIMGP)=reshape((/-xi1, -xi2,  xi3,  xi2,  xi1, -xi1, -xi2,  xi3,  xi2,  xi1, -xi1, -xi2,  xi3,  xi2,  xi1, -xi1, -xi2,  xi3,  xi2,  xi1, -xi1, -xi2,  xi3,  xi2,  xi1, &
                                                                     -xi1, -xi1, -xi1, -xi1, -xi1, -xi2, -xi2, -xi2, -xi2, -xi2,  xi3,  xi3,  xi3,  xi3,  xi3,  xi2,  xi2,  xi2,  xi2,  xi2,  xi1,  xi1,  xi1,  xi1,  xi1 /) ,(/NGP,NDIMGP/))
END MODULE

MODULE Integr2D_Square16GP
!     16 Integration Points (Gauss product rule)
!     GAUSS INTEGRATION POINTS & WEIGHTS
!     16 POINTS, INTEGRATION POINT NUMBERING:
!     13 14 15 16
!      9 10 11 12
!      5  6  7  8
!      1  2  3  4
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NGP=16, NDIMGP=2
  REAL(KIND=AbqRK), PARAMETER :: WGP1=0.3478548451_AbqRK, WGP2=0.6521451549_AbqRK, xi1=0.8611363116_AbqRK, xi2=0.3399810436_AbqRK
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPWeight(NGP)=(/WGP1*WGP1,WGP1*WGP2,WGP1*WGP2,WGP1*WGP1,WGP1*WGP2,WGP2*WGP2,WGP2*WGP2,WGP1*WGP2,WGP1*WGP2,WGP2*WGP2,WGP2*WGP2,WGP1*WGP2,WGP1*WGP1,WGP1*WGP2,WGP1*WGP2,WGP1*WGP1/)
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPPos(NGP,NDIMGP)=reshape((/-xi1, -xi2,  xi2,  xi1, -xi1, -xi2,  xi2,  xi1, -xi1, -xi2,  xi2,  xi1, -xi1, -xi2,  xi2,  xi1, &
                                                                     -xi1, -xi1, -xi1, -xi1, -xi2, -xi2, -xi2, -xi2,  xi2,  xi2,  xi2,  xi2,  xi1,  xi1,  xi1,  xi1 /) ,(/NGP,NDIMGP/))
END MODULE

! 4 and 9 GP with different numbering schemes (rowwise like in Abaqus and clockwise)
MODULE Integr2D_Square9GP
!     9 Integration Points (Gauss product rule)
!     GAUSS INTEGRATION POINTS & WEIGHTS
!     9 POINTS, INTEGRATION POINT NUMBERING (as in Abaqus Elements):
!     7 8 9
!     4 5 6
!     1 2 3
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NGP=9, NDIMGP=2
  REAL(KIND=AbqRK), PARAMETER :: WGP1=25.0_AbqRK/81.0_AbqRK, WGP2=40.0_AbqRK/81.0_AbqRK, WGP3=64.0_AbqRK/81.0_AbqRK
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPWeight(NGP)=(/WGP1,WGP2,WGP1,WGP2,WGP3,WGP2,WGP1,WGP2,WGP1/)
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPPos(NGP,NDIMGP)=SQRT(0.6_AbqRK)*reshape((/-1.0_AbqRK,  0.0_AbqRK,  1.0_AbqRK, -1.0_AbqRK, 0.0_AbqRK, 1.0_AbqRK, -1.0_AbqRK, 0.0_AbqRK, 1.0_AbqRK, &
                                                                                     -1.0_AbqRK, -1.0_AbqRK, -1.0_AbqRK,  0.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK,  1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK /) ,(/NGP,NDIMGP/))
END MODULE

MODULE Integr2D_Square9GP_ClockwiseNoScheme
!     9 Integration Points (Gauss product rule)
!     GAUSS INTEGRATION POINTS & WEIGHTS
!     9 POINTS, INTEGRATION POINT NUMBERING:
!     4 7 3
!     8 9 6
!     1 5 2
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NGP=9, NDIMGP=2
  REAL(KIND=AbqRK), PARAMETER :: WGP1=25.0_AbqRK/81.0_AbqRK, WGP2=40.0_AbqRK/81.0_AbqRK, WGP3=64.0_AbqRK/81.0_AbqRK
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPWeight(NGP)=(/WGP1,WGP1,WGP1,WGP1,WGP2,WGP2,WGP2,WGP2,WGP3/)
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPPos(NGP,NDIMGP)=SQRT(0.6_AbqRK)*reshape((/-1.0_AbqRK,  1.0_AbqRK, 1.0_AbqRK, -1.0_AbqRK,  0.0_AbqRK, 1.0_AbqRK, 0.0_AbqRK, -1.0_AbqRK, 0.0_AbqRK, &
                                                                                     -1.0_AbqRK, -1.0_AbqRK, 1.0_AbqRK,  1.0_AbqRK, -1.0_AbqRK, 0.0_AbqRK, 1.0_AbqRK,  0.0_AbqRK, 0.0_AbqRK /) ,(/NGP,NDIMGP/))
END MODULE

MODULE Integr2D_Square8GP
!     8 Integration Points =
!     *    FULL INTEGRATION FOR QUADRATIC SHAPE FUNCTIONS
!          AT REDUCED COMPUTATIONAL COSTS
!     GAUSS INTEGRATION POINTS & WEIGHTS
!     INTEGRATION POINT NUMBERING
!     3  7  4
!     8     6
!     1  5  2
!     Weights and locations from:
!     Dunavant: ECONOMICAL SYMMETRICAL QUADRATURE RULES FOR COMPLETE POLYNOMIALS OVER A SQUARE DOMAIN, Int. J. Num. Meth. Eng., 1985
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC ::NGP=8, NDIMGP=2
  REAL(KIND=AbqRK), PARAMETER, PRIVATE:: W1=0.183673469387755_AbqRK, W2=0.816326530612245_AbqRK, X1= 0.881917103688197_AbqRK, X2= 0.683130051063973_AbqRK
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPWeight(NGP)=(/W1,W1,W1,W1,W2,W2,W2,W2/)
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPPos(NGP,NDIMGP)=reshape((/-X1, X1,-X1, X1, 0.0_AbqRK, X2       ,  0.0_AbqRK, -X2       , &
                                                                     -X1,-X1, X1, X1, -X2      , 0.0_AbqRK, X2        , 0.0_AbqRK /) ,(/NGP,NDIMGP/))
END MODULE

MODULE Integr2D_Square7GP
!     7 Integration Points:
!     *    FULL INTEGRATION FOR QUADRATIC SHAPE FUNCTIONS
!          AT MINIMUM COMPUTATIONAL COSTS (However, not completely symmetric)
!     INTEGRATION POINT NUMBERING
!         7
!      3     4
!         5
!      1     2
!         6
!     Weights and locations of points from
!      Strout, 1971, p. 246, "*C_2: 5-1 (Radon-formula)
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NGP=7, NDIMGP=2
  REAL(KIND=AbqRK), PARAMETER, PRIVATE:: W1=5.0_AbqRK/9.0_AbqRK, W2=8.0_AbqRK/7.0_AbqRK, W3=20.0_AbqRK/63.0_AbqRK, &
                                         s35= SQRT(3.0_AbqRK/5.0_AbqRK), s15= SQRT(1.0_AbqRK/5.0_AbqRK), s13= SQRT(1.0_AbqRK/3.0_AbqRK), s1415=SQRT(14.0_AbqRK/15.0_AbqRK)
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPWeight(NGP)=(/W1,W1,W1,W1,W2,W3,W3/)
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPPos(NGP,NDIMGP)=reshape((/-s35,  s35, -s35, s35, 0.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, &
                                                                     -s13, -s13,  s13, s13, 0.0_AbqRK, -s1415   , s1415      /) ,(/NGP,NDIMGP/))
END MODULE

MODULE Integr2D_Square5GP
!     5 Integration Points=
!     *    5th order if integrand exhibits central symmetry f(x,y)=f(y,x)
!     GAUSS INTEGRATION POINTS & WEIGHTS
!     INTEGRATION POINT NUMBERING
!      3     4
!         5
!      1     2
!     Weights and locations of points from
!      Cohen, Gismalla: "Some integration formulae for symmetric functions of two variables", International Journal of Computer Mathematics, 1986
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NGP=5, NDIMGP=2
  REAL(KIND=AbqRK), PARAMETER, PRIVATE:: W1=5.0_AbqRK/7.0_AbqRK, W2=8.0_AbqRK/7.0_AbqRK, u=0.84623312_AbqRK, v=0.46607171_AbqRK
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPWeight(NGP)=(/W1,W1,W1,W1,W2/)
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPPos(NGP,NDIMGP)=reshape((/-v, u, -u, v,  0.0_AbqRK, &
                                                                     -u,-v,  v, u,  0.0_AbqRK /) ,(/NGP,NDIMGP/))
END MODULE

MODULE Integr2D_Square5GPsym
!     5 Integration Points=
!     *    3rd order, i.e. Reduced INTEGRATION FOR QUADRATIC SHAPE FUNCTIONS
!          but 5th point eliminates hour-glass mode
!     GAUSS INTEGRATION POINTS & WEIGHTS
!     INTEGRATION POINT NUMBERING
!      3     4
!         5
!      1     2
!     Weights and locations of points from
!      Strout, 1971, p. 244, "C_2: 3-5 (Ewing-formula)
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NGP=5, NDIMGP=2
  REAL(KIND=AbqRK), PARAMETER,PRIVATE :: W1=8.0_AbqRK/3.0_AbqRK, W2=1.0_AbqRK/3.0_AbqRK
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPWeight(NGP)=(/W1,W1,W1,W1,W2/)
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPPos(NGP,NDIMGP)=reshape((/-1.0_AbqRK,  1.0_AbqRK, -1.0_AbqRK, 1.0_AbqRK, 0.0_AbqRK, &
                                                                     -1.0_AbqRK, -1.0_AbqRK,  1.0_AbqRK, 1.0_AbqRK, 0.0_AbqRK /) ,(/NGP,NDIMGP/))
END MODULE

MODULE Integr2D_Square4GP
!     4 Integration Points (Gauss product rule)=
!     *    REDUCED INTEGRATION FOR QUADRATIC SHAPE FUNCTIONS
!     *    FULL INTEGRATION FOR LINEAR SHAPE FUNCTIONS
!     GAUSS INTEGRATION POINTS & WEIGHTS
!     4 POINTS, INTEGRATION POINT NUMBERING (as in Abaqus elements):
!     3 4
!     1 2
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NGP=4, NDIMGP=2
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPWeight(NGP)=(/1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK/)
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPPos(NGP,NDIMGP)=SQRT(1.0_AbqRK/3.0_AbqRK)*reshape((/-1.0_AbqRK,  1.0_AbqRK, -1.0_AbqRK, 1.0_AbqRK, &
                                                                                               -1.0_AbqRK, -1.0_AbqRK,  1.0_AbqRK, 1.0_AbqRK /) ,(/NGP,NDIMGP/))
END MODULE

MODULE Integr2D_Square4GP_ClockwiseNoScheme
!     4 Integration Points (Gauss product rule)
!     *    REDUCED INTEGRATION FOR QUADRATIC SHAPE FUNCTIONS
!     *    FULL INTEGRATION FOR LINEAR SHAPE FUNCTIONS
!     GAUSS INTEGRATION POINTS & WEIGHTS
!     4 POINTS, INTEGRATION POINT NUMBERING (as in Abaqus elements):
!     4 3
!     1 2
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NGP=4, NDIMGP=2
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPWeight(NGP)=(/1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK /)
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPPos(NGP,NDIMGP)=SQRT(1.0_AbqRK/3.0_AbqRK)*reshape((/ -1.0_AbqRK,  1.0_AbqRK, 1.0_AbqRK, -1.0_AbqRK, &
                                                                                                -1.0_AbqRK, -1.0_AbqRK, 1.0_AbqRK,  1.0_AbqRK /) ,(/NGP,NDIMGP/))
END MODULE

MODULE Integr2D_Square1GP
!     1 Integration Point (Gauss product rule)
!     *    REDUCED INTEGRATION FOR LINEAR SHAPE FUNCTIONS
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NGP=1, NDIMGP=2
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPWeight(NGP)=(/4.0_AbqRK/)
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPPos(NGP,NDIMGP)=reshape((/0.0_AbqRK, 0.0_AbqRK/) ,(/NGP,NDIMGP/))
END MODULE

! Newton-Cotes: 9 GP 
MODULE Integr2D_Square9GP_NC
!     9 Integration Points (Newton-Cotes quadrature rule)
!     INTEGRATION POINTS & WEIGHTS
!     9 POINTS, INTEGRATION POINT NUMBERING (as in Abaqus Elements):
!     7 8 9
!     4 5 6
!     1 2 3
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NGP=9, NDIMGP=2
  REAL(KIND=AbqRK), PARAMETER :: WGP1=1.0_AbqRK/36.0_AbqRK, WGP2=1.0_AbqRK/9.0_AbqRK, WGP3=4.0_AbqRK/9.0_AbqRK
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPWeight(NGP)=(/WGP1,WGP2,WGP1,WGP2,WGP3,WGP2,WGP1,WGP2,WGP1/)
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPPos(NGP,NDIMGP)=reshape((/-1.0_AbqRK,  0.0_AbqRK,  1.0_AbqRK, -1.0_AbqRK, 0.0_AbqRK, 1.0_AbqRK, -1.0_AbqRK, 0.0_AbqRK, 1.0_AbqRK, &
                                                                     -1.0_AbqRK, -1.0_AbqRK, -1.0_AbqRK,  0.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK,  1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK /) ,(/NGP,NDIMGP/))
END MODULE

!!!!!!!!!!!!!!!!!!!!!
!! 3D, tetrahedron
!!!!!!!!!!!!!!!!!!!!!

MODULE Integr3D_Tetrahed1GP
! integration points & weights.
! 1 point.
! From Abaqus Theory Manual, paragraph 3.2.6:
! "The first-order triangle and tetrahedron are constant stress elements and use a single
! integration point for the stiffness calculation when used in stress/displacement
! applications."
    USE ABQINTERFACE
    PRIVATE
    INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NGP = 1, NDIMGP = 3
    REAL(KIND=AbqRK), PARAMETER :: WGP=1.0_AbqRK/6.0_AbqRK
    REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPWeight(NGP)=(/WGP/)
    REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPPos(NGP,NDIMGP)=1.0_AbqRK/4.0_AbqRK*reshape((/1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK/),(/NGP,NDIMGP/))
END MODULE Integr3D_Tetrahed1GP

MODULE Integr3D_Tetrahed4GP
! integration points & weights.
! 4 points.
! From Abaqus Theory Manual, paragraph 3.2.6:
! "For stress/displacement applications the second-order tetrahedron uses
! 4 integration points for its stiffness matrix..."
! Weights and locations of points from
!   Detlev Maurer: "Numerical Methods", FH Landshut (online), https://people.fh-landshut.de/~maurer/numeth/numeth.html
! Integration point numbering:
! Looking from the top in the nagative Z direction:
!           ^y (y points to the top)
!          |->x (x points to the right)
!
! top:    4
! bottom: 3
!         1 2
    USE ABQINTERFACE
    PRIVATE
    INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NGP = 4, NDIMGP = 3
    REAL(KIND=AbqRK), PARAMETER :: WGP=1.0_AbqRK/24.0_AbqRK
    REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPWeight(NGP)=(/WGP,WGP,WGP,WGP/)
    REAL(KIND=AbqRK), PARAMETER :: pos1=(5.0_AbqRK+3.0_AbqRK*SQRT(5.0_AbqRK))/20.0_AbqRK, pos2=(5.0_AbqRK-SQRT(5.0_AbqRK))/20.0_AbqRK
    REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPPos(NGP,NDIMGP)=reshape((/ &
           pos2, pos1, pos2, pos2,&
           pos2, pos2, pos1, pos2,&
           pos2, pos2, pos2, pos1/),(/NGP,NDIMGP/))
END MODULE Integr3D_Tetrahed4GP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 3D, hexagon, reduced and full integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE Integr3D_Hex27GP
!     Integration Points (Gauss product rule)
!     27 POINTS, INTEGRATION POINT NUMBERING:
!          Looking from the top in the negative z direction:
!
!          ^y (y points to the top)
!          |->x (x points to the right)
!
!   @z=-1:  7  8   9
!           4  5   6
!           1  2   3
!
!   @z=0:   16 17 18
!           13 14 15
!           10 11 12
!
!   @z=1:   25 26 27
!           22 23 24
!           19 20 21
  USE ABQINTERFACE
  IMPLICIT NONE
  PUBLIC
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NGP=27, NDIMGP=3
  REAL(KIND=AbqRK), PARAMETER :: WGP1=125.0_AbqRK/729.0_AbqRK, WGP2=200.0_AbqRK/729.0_AbqRK, WGP3=320.0_AbqRK/729.0_AbqRK, WGP4=512.0_AbqRK/729.0_AbqRK
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPWeight(NGP)=(/WGP1,WGP2,WGP1,WGP2,WGP3,WGP2,WGP1,WGP2,WGP1,&
                                                         WGP2,WGP3,WGP2,WGP3,WGP4,WGP3,WGP2,WGP3,WGP2,&
                                                         WGP1,WGP2,WGP1,WGP2,WGP3,WGP2,WGP1,WGP2,WGP1/)
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPPos(NGP,NDIMGP)=SQRT(0.6_AbqRK)*reshape((/-1.0_AbqRK, 0.0_AbqRK, 1.0_AbqRK,-1.0_AbqRK, 0.0_AbqRK, 1.0_AbqRK,-1.0_AbqRK, 0.0_AbqRK, 1.0_AbqRK, &
                                                                                     -1.0_AbqRK, 0.0_AbqRK, 1.0_AbqRK,-1.0_AbqRK, 0.0_AbqRK, 1.0_AbqRK,-1.0_AbqRK, 0.0_AbqRK, 1.0_AbqRK, &
                                                                                     -1.0_AbqRK, 0.0_AbqRK, 1.0_AbqRK,-1.0_AbqRK, 0.0_AbqRK, 1.0_AbqRK,-1.0_AbqRK, 0.0_AbqRK, 1.0_AbqRK, &
                                                                                     -1.0_AbqRK,-1.0_AbqRK,-1.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, &
                                                                                     -1.0_AbqRK,-1.0_AbqRK,-1.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, &
                                                                                     -1.0_AbqRK,-1.0_AbqRK,-1.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, &
                                                                                     -1.0_AbqRK,-1.0_AbqRK,-1.0_AbqRK,-1.0_AbqRK,-1.0_AbqRK,-1.0_AbqRK,-1.0_AbqRK,-1.0_AbqRK,-1.0_AbqRK, &
                                                                                      0.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, &
                                                                                      1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK  /),(/NGP,NDIMGP/))
END MODULE

MODULE Integr3D_Hex27GP_ClockwiseNoScheme
!     Integration Points (Gauss product rule)
!     27 POINTS, INTEGRATION POINT NUMBERING:
!          Looking from the top in the negative z direction:
!
!          ^y (y points to the top)
!          |->x (x points to the right)
!
!   @z=-1:  4  11  3
!           12 13 10
!           1  9   2
!
!   @z=1:   8  21  7
!           22 23 20
!           5  19  6
!
!   @z=0:   27 16 26
!           17 18 15
!           24 14 25
  USE ABQINTERFACE
  IMPLICIT NONE
  PUBLIC
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NGP=27, NDIMGP=3
  REAL(KIND=AbqRK), PARAMETER :: WGP1=125.0_AbqRK/729.0_AbqRK, WGP2=200.0_AbqRK/729.0_AbqRK, WGP3=320.0_AbqRK/729.0_AbqRK, WGP4=512.0_AbqRK/729.0_AbqRK
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPWeight(NGP)=(/WGP1,WGP1,WGP1,WGP1,WGP1,WGP1,WGP1,WGP1,&
                                                         WGP2,WGP2,WGP2,WGP2,WGP3,WGP3,WGP3,WGP3,&
                                                         WGP3,WGP4,WGP2,WGP2,WGP2,WGP2,WGP3,WGP2,&
                                                         WGP2,WGP2,WGP2/)
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPPos(NGP,NDIMGP)=SQRT(0.6_AbqRK)*reshape((/-1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK,-1.0_AbqRK,-1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK,-1.0_AbqRK, 0.0_AbqRK, &
                                                                                      1.0_AbqRK, 0.0_AbqRK,-1.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, 1.0_AbqRK, 0.0_AbqRK,-1.0_AbqRK, 0.0_AbqRK, &
                                                                                      0.0_AbqRK, 1.0_AbqRK, 0.0_AbqRK,-1.0_AbqRK, 0.0_AbqRK,-1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK,-1.0_AbqRK, &
                                                                                     -1.0_AbqRK,-1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK,-1.0_AbqRK,-1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK,-1.0_AbqRK, &
                                                                                      0.0_AbqRK, 1.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK,-1.0_AbqRK, 0.0_AbqRK, 1.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, &
                                                                                     -1.0_AbqRK, 0.0_AbqRK, 1.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK,-1.0_AbqRK,-1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, &
                                                                                     -1.0_AbqRK,-1.0_AbqRK,-1.0_AbqRK,-1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK,-1.0_AbqRK, &
                                                                                     -1.0_AbqRK,-1.0_AbqRK,-1.0_AbqRK,-1.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, &
                                                                                      1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK  /),(/NGP,NDIMGP/))
END MODULE

MODULE Integr3D_Hex14GP
!     GAUSS INTEGRATION POINTS & WEIGHTS
!     14 POINTS, same accuracy as 3x3x3=27 Points
!     Numbering: 1-8 like Abaqus elements with 8 IP, 9-14 like faces of Abaqus elements
!     Weights and locations from:
!       B. C. Iron: QUADRATURE RULES FOR BRICK BASED FINITE ELEMENTS, 1971
  USE ABQINTERFACE
  IMPLICIT NONE
  PUBLIC
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NGP=14, NDIMGP=3
  REAL(KIND=AbqRK), PARAMETER, PRIVATE:: WGP1= 0.335180055_AbqRK, WGP2=0.886426593_AbqRK, X1= 0.758786911_AbqRK, X2=0.795822426_AbqRK
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPWeight(NGP)=(/WGP1,WGP1,WGP1,WGP1,WGP1,WGP1,WGP1,WGP1,WGP2,WGP2,WGP2,WGP2,WGP2,WGP2/)
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPPos(NGP,NDIMGP)=reshape((/-X1, X1,-X1, X1,-X1, X1,-X1, X1,   0.0_AbqRK,  0.0_AbqRK, 0.0_AbqRK,        X2, 0.0_AbqRK,       -X2, &
                                                                     -X1,-X1, X1, X1,-X1,-X1, X1, X1,   0.0_AbqRK,  0.0_AbqRK,       -X2, 0.0_AbqRK,        X2, 0.0_AbqRK, &
                                                                     -X1,-X1,-X1,-X1, X1, X1, X1, X1,         -X2,         X2, 0.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK  /),(/NGP,NDIMGP/))
END MODULE

MODULE Integr3D_Hex8GP
!     Integration Points (Gauss product rule)
!     8 POINTS, INTEGRATION POINT NUMBERING (as in Abaqus Elements):
!          Looking from the top in the negative z direction:
!          ^y (y points to the top)
!          |->x (x points to the right)
!   @z=-1:  3  4
!           1  2
!   @z=1:   7  8
!           5  6
  USE ABQINTERFACE
  IMPLICIT NONE
  PUBLIC
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NGP=8, NDIMGP=3
  REAL(KIND=AbqRK), PARAMETER :: WGP1=1.0_AbqRK
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPWeight(NGP)=(/WGP1,WGP1,WGP1,WGP1,WGP1,WGP1,WGP1,WGP1/)
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPPos(NGP,NDIMGP)=SQRT(1.0_AbqRK/3.0_AbqRK)*reshape((/-1.0_AbqRK, 1.0_AbqRK,-1.0_AbqRK, 1.0_AbqRK,-1.0_AbqRK, 1.0_AbqRK,-1.0_AbqRK, 1.0_AbqRK, &
                                                                                               -1.0_AbqRK,-1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK,-1.0_AbqRK,-1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, &
                                                                                               -1.0_AbqRK,-1.0_AbqRK,-1.0_AbqRK,-1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK  /),(/NGP,NDIMGP/))
END MODULE

MODULE Integr3D_Hex8GP_ClockwiseNoScheme
!     Integration Points (Gauss product rule)
!     8 POINTS, INTEGRATION POINT NUMBERING:
!          Looking from the top in the negative z direction:
!          ^y (y points to the top)
!          |->x (x points to the right)
!   @z=-1:  4  3
!           1  2
!   @z=1:   8  7
!           5  6
  USE ABQINTERFACE
  IMPLICIT NONE
  PUBLIC
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NGP=8, NDIMGP=3
  REAL(KIND=AbqRK), PARAMETER :: WGP1=1.0_AbqRK
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPWeight(NGP)=(/WGP1,WGP1,WGP1,WGP1,WGP1,WGP1,WGP1,WGP1/)
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPPos(NGP,NDIMGP)=SQRT(1.0_AbqRK/3.0_AbqRK)*reshape((/-1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK,-1.0_AbqRK,-1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK,-1.0_AbqRK, &
                                                                                               -1.0_AbqRK,-1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK,-1.0_AbqRK,-1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, &
                                                                                               -1.0_AbqRK,-1.0_AbqRK,-1.0_AbqRK,-1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK, 1.0_AbqRK  /),(/NGP,NDIMGP/))
END MODULE

MODULE Integr3D_Hex1GP
!     1 Integration Point (Gauss product rule)
!     *    REDUCED INTEGRATION FOR LINEAR SHAPE FUNCTIONS
  USE ABQINTERFACE
  IMPLICIT NONE
  PRIVATE
  INTEGER(KIND=AbqIK), PARAMETER, PUBLIC :: NGP=1, NDIMGP=3
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPWeight(NGP)=(/8.0_AbqRK/)
  REAL(KIND=AbqRK), PARAMETER, PUBLIC :: GPPos(NGP,NDIMGP)=reshape((/0.0_AbqRK, 0.0_AbqRK, 0.0_AbqRK/) ,(/NGP,NDIMGP/))
END MODULE

