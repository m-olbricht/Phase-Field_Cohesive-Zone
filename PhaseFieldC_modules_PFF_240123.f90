!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! UMAT material subroutine of multi-damage multi-component damage field model
!
! Stephan Roth, TU Bergakademie Freiberg, 30.07.2020
!
! 30.07.2020: Multi-damage multi-component
! 23.11.2020: all derivatives tested (totalPotential --> stresses --> tangent)
! 15.06.2021: concentrations independent on damages
! 19.07.2021: damage-dependent parameters
! 04.03.2022: chemical reactions
! 03.01.2024: without NCP, DOF: displacements, 1xdamage variable
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SAVED VARIABLES AND INPUT PARAMETERS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SVR: svr( 1)     : damage
!	   svr( 2)	   : stress 11
!	   svr( 3)	   : stress 22
!	   svr( 4)	   : stress 33
!	   svr( 5)	   : stress 12
!	   svr( 6)	   : stress 13
!	   svr( 7)	   : stress 23
!	   svr( 8)	   : stran 11
!	   svr( 9)	   : stran 22
!	   svr(10)	   : stran 33
!	   svr(11)	   : stran 12
!	   svr(12)	   : stran 13
!	   svr(13)	   : stran 23
!	   svr(14)	   : H
!	   svr(15)	   : -
!	   svr(16)	   : elast energy
!	   svr(17)	   : dissipat plast energy
!	   svr(18)	   : total Pot
!	   svr(19)	   : -
!	   svr(20)	   : -
!
! material parameters:
!
!		props_mat ( 1) : 
!		props_mat ( 2) : 
!		props_mat ( 3) : 
!		props_mat ( 4) : *numerical Tangent
!		props_mat ( 5) : E-Modul
!		props_mat ( 6) : Poisson-ratio
!		props_mat ( 7) : 
!		props_mat ( 8) : G0
!		props_mat ( 9) : l0
!		props_mat (10) : prop_df_alpha
!		props_mat (11) : prop_df_beta
!		props_mat (12) : 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE PhaseField_module

  IMPLICIT NONE

  PUBLIC :: CheckMaterialParameters, umatPF

  CONTAINS

!------------------------------------------------------------------------------------

    SUBROUTINE CheckMaterialParameters(props)
      ! Check of all material parameters

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE ABAModul
      USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: props(numMatPar)
      REAL(kind=AbqRK) :: prop_E, prop_nu

      ! read material parameters from props
      prop_E      = props( 5) ! Young's modulus
      prop_nu     = props( 6) ! Poisson's ratio

      ! parameter check
      IF (prop_E .LE. zero) THEN
        WRITE(7,*) "Young's modulus should exceed zero. EXIT", prop_E
        CALL XEXIT()
      END IF
      IF (prop_nu .LT. zero) THEN
        WRITE(7,*) "Poisson's ratio must be positive. EXIT", prop_nu
        CALL XEXIT()
      END IF

    END SUBROUTINE CheckMaterialParameters

!------------------------------------------------------------------------------------

    SUBROUTINE umatPF(stress,svr,Ct,energy_elast,dissipat_plast,dissipat_creep,rpl,ddsddt,drplde,drpldt, &
                      stran,dstran,time,dtime,Temp,dTemp,predef_umat,dpred,cmname,ndi,nshr,ntens, &
                      num_svr,props_mat,num_matpar,coords_gp,Trafomat,pnewdt,intlen,F0,F1,jelem, &
                      npt,layer,kspt,kstep,kinc)
      ! compute generalised stresses and tangent

      USE ABQINTERFACE_PF
      USE FLOATNUMBERS
      USE ABAModul
      USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
      
      USE FreeEnergyModule
      USE AliasModulePF 
      USE BulkEnergyModule, ONLY: HFEDtens_H

      IMPLICIT NONE

      ! UMAT variables
      INTEGER(kind=AbqIK), INTENT(IN) :: ndi, nshr, ntens, num_svr, num_matpar, &
                                         jelem, npt, layer, kspt , kstep, kinc
      REAL(kind=AbqRK), INTENT(IN) :: stran(ntens), dstran(ntens), time(2), &
                                      dtime, Temp, dTemp, predef_umat(*), dpred(*), &
                                      coords_gp(3), Trafomat(3,3), intlen, F0(3,3), &
                                      F1(3,3), props_mat(num_matpar)
      REAL(kind=AbqRK), INTENT(INOUT) :: stress(ntens), svr(num_svr), &
                                         energy_elast, dissipat_plast, dissipat_creep, pnewdt
      REAL(kind=AbqRK), INTENT(OUT) :: Ct(ntens,ntens), rpl, ddsddt(ntens), &
                                       drplde(ntens), drpldt
      CHARACTER*45, INTENT(INOUT) :: cmname

      ! further variables
      INTEGER(kind=AbqIK) :: D, nphase, i1, i2
      INTEGER(kind=AbqIK) :: nHFEDpar, nIEDpar
      INTEGER(kind=AbqIK) :: first_HFED_par_index, first_IED_par_index
      REAL(kind=AbqRK) :: eps(3,3)
      REAL(kind=AbqRK) :: delta, stran_per(ntens), stress_per(ntens), Ct_temp(ntens,ntens), totalPotential, per, stress_num(ntens)
      REAL(kind=AbqRK) :: prop_df_a, prop_df_b
      REAL(kind=AbqRK) :: ElasticEnergytens, H, Hn
      INTEGER(kind=AbqIK), ALLOCATABLE :: pos_p(:)
      REAL(kind=AbqRK), ALLOCATABLE :: damage(:), grad_damage(:,:)
      REAL(kind=AbqRK), ALLOCATABLE :: parHFEDMatrix(:), parIEDMatrix(:)
      LOGICAL :: numericalTangent, isSpherisymmetric, isNaN

      numericalTangent = .FALSE.; isSpherisymmetric = .FALSE.; isNaN = .FALSE.

      ! dimension
      D = 0 ! dummy
      SELECT CASE(ndi+nshr)
      CASE(1)
        ! 1D
        D = 1
      CASE(3)
        ! 1D, spheri-symmetric
        D = 1
        isSpherisymmetric = .TRUE.
      CASE(4)
        ! 2D
        D = 2
      CASE(6)
        ! 3D
        D = 3
      END SELECT

      ! numerical tangent
      IF (props_mat(4) .NE. zero) numericalTangent = .TRUE.

      ! number of damages
      nphase = 1 ! damage variable
      !
      ! number of parameters in free energy
      nHFEDpar = 2
      ! number of parameters in interface energy
      nIEDpar = 2
      !
      !
      ALLOCATE(pos_p(nphase), damage(nphase), grad_damage(nphase,3))
      ALLOCATE(parHFEDMatrix(nHFEDpar), parIEDMatrix(nIEDpar))
      pos_p = 0
      damage = zero; grad_damage = zero
      parHFEDMatrix = zero; parIEDMatrix = zero

      ! coordinate positions
      DO i1=1,nphase
        pos_p(i1) = ndi+nshr+(i1-1)*(D+1)+1
      END DO

      ! generalized kinematic measures
      eps = zero; damage = zero; grad_damage = zero
      ! strain tensor
      eps = strainCoordinates(D,ntens,stran,isSpherisymmetric)
      ! order parameter -- damage variable
      damage = damageParameter(nphase,pos_p,ntens,stran)
      ! gradient of order parameter = gradient of damage variable
      grad_damage(:,:) = gradientOfPhaseParameter(nphase,D,pos_p,ntens,stran)

      ! saved variables of last converged increment

      ! read material parameters from props_mat
      ! bulk free energy density
      first_HFED_par_index = 5
      parHFEDMatrix = props_mat(first_HFED_par_index:first_HFED_par_index+nHFEDpar-1)
      ! interface energy density -- see InterfaceEnergyModule.f90
      first_IED_par_index = 8
      parIEDMatrix = props_mat(first_IED_par_index:first_IED_par_index+nIEDpar-1)
	  
	  ! Degradation Parameters
	  prop_df_a = props_mat(10)
	  prop_df_b = props_mat(11)

	  ! History Update
      Hn = svr(14)
      H = zero
      
      ElasticEnergytens = HFEDtens_H(eps,nHFEDpar,parHFEDMatrix)
      
	  IF (ElasticEnergytens .GT. Hn) THEN
        H = ElasticEnergytens
      ELSE
        H = Hn
      END IF 
      
      svr(14) = H
	  
	  
      ! compute all first derivatives

      ! generalised stresses
      stress = stresses(D,ntens,nphase,pos_p,eps,damage,grad_damage, &
                        nHFEDpar,parHFEDMatrix,nIEDpar,parIEDMatrix, &
                        prop_df_a, prop_df_b, &
                        H)

      ! energetic quantities

      ! total potential
      totalPotential = potential(D,nphase,eps,damage,grad_damage,nHFEDpar,parHFEDMatrix,nIEDpar,parIEDMatrix, prop_df_a, prop_df_b)

      ! Helmholtz free energy density of damage mix and Helmholtz energy density of interface
      energy_elast = HFED_Mix(eps,nHFEDpar,parHFEDMatrix)

      ! Interface energy density
      dissipat_plast = IED(D,damage(1),grad_damage(1,1:D),nIEDpar,parIEDMatrix)

      ! dissipation
      !dissipat_creep = totalPotential

      ! saved variables
      !
      ! stress coordinates
      svr(2:1+ndi+nshr) = stress(1:ndi+nshr)
      !
      ! strain coordinates
      svr(8:7+ndi+nshr) = stran(1:ndi+nshr)
      !
      ! integration point coordinates
      !
      ! damages
      svr(1) = damage(1)
      !
      ! gradients of damages
      !svr(2:1+D) = grad_damage(1,1:D)
      !
      ! reversible energy density
      svr(16) = energy_elast
      svr(17) = dissipat_plast
      svr(18) = totalPotential

      ! generalised tangent
      IF (numericalTangent) THEN

        ! numerical tangent obtained with perturbation of stresses
        delta = 1.d-6
        Ct = zero
        DO i1=1,ntens
          stran_per = stran
          stran_per(i1) = stran_per(i1) + delta
          ! generalized kinematic measures
          eps = zero; damage = zero; grad_damage = zero
          ! strain tensor
          eps = strainCoordinates(D,ntens,stran_per,isSpherisymmetric)
          ! order parameter -- damage variable
          damage = damageParameter(nphase,pos_p,ntens,stran_per)
          ! gradient of order parameter = gradient of damage variable
          grad_damage(:,:) = gradientOfPhaseParameter(nphase,D,pos_p,ntens,stran_per)
          !
          ! generalised stresses
          stress_per = stresses(D,ntens,nphase,pos_p,eps,damage,grad_damage, &
                                nHFEDpar,parHFEDMatrix,nIEDpar,parIEDMatrix, &
                                prop_df_a, prop_df_b, &
                                H)
          ! tangent
          Ct(i1,1:ntens) = (stress_per(1:ntens)-stress(1:ntens))/delta
        END DO
!        Ct_temp = tangent(D,ntens,nphase,pos_p,eps,damage,grad_damage, &
!                     nHFEDpar,parHFEDMatrix,nIEDpar,parIEDMatrix)

      ELSE
        ! analytical material tangent
        Ct = tangent(D,ntens,nphase,pos_p,eps,damage,grad_damage, &
                     nHFEDpar,parHFEDMatrix,nIEDpar,parIEDMatrix, &
                     prop_df_a, prop_df_b, &
                     H)
      END IF

      ! dummy output
      rpl = zero; ddsddt = zero; drplde = zero; drpldt = zero

      DEALLOCATE(pos_p, damage, grad_damage)
      DEALLOCATE(parHFEDMatrix, parIEDMatrix)

      CONTAINS

!------------------------------------------------------------------------------------

      PURE REAL(kind=AbqRK) FUNCTION potential(D,nphase,eps,damage,grad_damage,nHFEDpar,parHFEDMatrix,nIEDpar,parIEDMatrix, &
												prop_df_a, prop_df_b)
      ! computes total potential

        USE ABQINTERFACE_PF
        USE FLOATNUMBERS
        USE FreeEnergyModule
        USE AliasModulePF
        USE BulkEnergyModule

        IMPLICIT NONE
        INTEGER(kind=AbqIK), INTENT(IN) :: D, nphase, nHFEDpar, nIEDpar
        REAL(kind=AbqRK), INTENT(IN) :: eps(3,3), damage(nphase), grad_damage(nphase,3)
        REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrix(nHFEDpar), parIEDMatrix(nIEDpar) 
        REAL(kind=AbqRK), INTENT(IN) :: prop_df_a, prop_df_b
        INTEGER(kind=AbqIK) :: i1

        potential = bulkED(eps,damage(1),nHFEDpar,parHFEDMatrix,prop_df_a,prop_df_b) &
				  + IED(D,damage(1),grad_damage(1,1:D),nIEDpar,parIEDMatrix)

      END FUNCTION potential

!------------------------------------------------------------------------------------

      PURE REAL(kind=AbqRK) FUNCTION strainCoordinates(D,ntens,s,isSpherisymmetric)
      ! extract strain coordinates

        USE ABQINTERFACE_PF
        USE FLOATNUMBERS

        IMPLICIT NONE
        INTEGER(kind=AbqIK), INTENT(IN) :: D, ntens
        REAL(kind=AbqRK), INTENT(IN) :: s(ntens)
        LOGICAL, INTENT(IN) :: isSpherisymmetric
        DIMENSION strainCoordinates(3,3)

        strainCoordinates = zero
        !
        ! strain tensor
        strainCoordinates(1,1) = s(1)
        IF ( (D .GE. 2) .OR. ( (D .EQ. 1) .AND. (isSpherisymmetric) ) ) THEN
          strainCoordinates(2,2) = s(2)
          strainCoordinates(3,3) = s(3)
        END IF
        IF (D .GE. 2) THEN
          strainCoordinates(1,2) = s(4)/two
          strainCoordinates(2,1) = strainCoordinates(1,2)
        END IF
        IF (D .EQ. 3) THEN
          strainCoordinates(1,3) = s(5)/two
          strainCoordinates(2,3) = s(6)/two
          strainCoordinates(3,1:2) = strainCoordinates(1:2,3)
        END IF

      END FUNCTION strainCoordinates

!------------------------------------------------------------------------------------

      PURE REAL(kind=AbqRK) FUNCTION damageParameter(nphase,pos_p,ntens,s)
      ! extract damage parameter

        USE ABQINTERFACE_PF
        USE FLOATNUMBERS

        IMPLICIT NONE
        INTEGER(kind=AbqIK), INTENT(IN) :: nphase, pos_p(nphase), ntens
        REAL(kind=AbqRK), INTENT(IN) :: s(ntens)
        INTEGER(kind=AbqIK) :: i1
        DIMENSION damageParameter(nphase)

        damageParameter = zero
        FORALL (i1=1:nphase)
          damageParameter(i1) = s(pos_p(i1))
        END FORALL

      END FUNCTION damageParameter

!------------------------------------------------------------------------------------

      PURE REAL(kind=AbqRK) FUNCTION gradientOfPhaseParameter(nphase,D,pos_p,ntens,s)
      ! extract gradient of damage parameter

        USE ABQINTERFACE_PF
        USE FLOATNUMBERS

        IMPLICIT NONE
        INTEGER(kind=AbqIK), INTENT(IN) :: nphase, D, pos_p(nphase), ntens
        REAL(kind=AbqRK), INTENT(IN) :: s(ntens)
        INTEGER(kind=AbqIK) :: i1, i2
        DIMENSION gradientOfPhaseParameter(nphase,3)

        gradientOfPhaseParameter = zero
        FORALL (i1=1:nphase)
          gradientOfPhaseParameter(i1,1:D) = s(pos_p(i1)+1:pos_p(i1)+D)
        END FORALL

      END FUNCTION gradientOfPhaseParameter

!------------------------------------------------------------------------------------

      PURE REAL(kind=AbqRK) FUNCTION stresses(D,ntens,nphase, & 
                                              pos_p,eps,damage,grad_damage, &
                                              nHFEDpar,parHFEDMatrix,nIEDpar,parIEDMatrix, &
                                              prop_df_a, prop_df_b, &
                                              H)

      ! generalised stresses, modified Voigt notation

        USE ABQINTERFACE_PF
        USE FLOATNUMBERS
        USE FreeEnergyModule
        USE AliasModulePF
        USE BulkEnergyModule

        IMPLICIT NONE
        INTEGER(kind=AbqIK), INTENT(IN) :: D, ntens, nphase, pos_p(nphase)
        INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar, nIEDpar
        REAL(kind=AbqRK), INTENT(IN) :: eps(3,3), damage(nphase), grad_damage(nphase,3)
        REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrix(nHFEDpar), parIEDMatrix(nIEDpar) 
        REAL(kind=AbqRK), INTENT(IN) :: prop_df_a, prop_df_b
        REAL(kind=AbqRK), INTENT(IN) :: H
        REAL(kind=AbqRK) :: d_BKED_d_eps(3,3)
        REAL(kind=AbqRK) :: temp(3)
        INTEGER(kind=AbqIK) :: i1, pos_i1
        DIMENSION stresses(ntens)

        d_BKED_d_eps = zero

        ! first derivatives of bulk energy density
        d_BKED_d_eps = d_bulkED_d_eps(eps,damage(1),nHFEDpar,parHFEDMatrix,prop_df_a,prop_df_b)
        !
        stresses=zero
        !
        ! eps: contributions of bulk energy density (BKED)
        stresses(1) = d_BKED_d_eps(1,1)
        stresses(2) = d_BKED_d_eps(2,2)
        stresses(3) = d_BKED_d_eps(3,3)
        IF (D .GE. 2) THEN
          stresses(4) = d_BKED_d_eps(1,2)
        END IF
        IF (D .EQ. 3) THEN
          stresses(5) = d_BKED_d_eps(1,3)
          stresses(6) = d_BKED_d_eps(2,3)
        END IF
        !
        ! damage parameters
        !
        DO i1=1,nphase
          pos_i1 = pos_p(i1)
          ! damage: contributions of interface energy
          stresses(pos_i1) = d_bulkED_d_damage_H(eps,damage(1),nHFEDpar,parHFEDMatrix,prop_df_a,prop_df_b,H) &
						   + d_IED_d_damage(D,damage(i1),grad_damage(i1,1:D),nIEDpar,parIEDMatrix)
          !
          ! grad_damage: contributions of interface energy
          temp = zero
          temp(1:D) = d_IED_d_grad_damage(D,damage(i1),grad_damage(i1,1:D),nIEDpar,parIEDMatrix)
          stresses(pos_i1+1:pos_i1+D) = temp(1:D)
        END DO

      END FUNCTION stresses

!------------------------------------------------------------------------------------

      PURE REAL(kind=AbqRK) FUNCTION tangent(D,ntens,nphase, & 
                                             pos_p,eps,damage,grad_damage, &
                                             nHFEDpar,parHFEDMatrix,nIEDpar,parIEDMatrix, &
                                             prop_df_a, prop_df_b, &
                                             H)

      ! generalised tangent, modified Voigt notation

        USE ABQINTERFACE_PF
        USE FLOATNUMBERS
        USE FreeEnergyModule
        USE AliasModulePF
        USE BulkEnergyModule

        IMPLICIT NONE
        INTEGER(kind=AbqIK), INTENT(IN) :: D, ntens, nphase, pos_p(nphase)
        INTEGER(kind=AbqIK), INTENT(IN) :: nHFEDpar, nIEDpar
        REAL(kind=AbqRK), INTENT(IN) :: eps(3,3), damage(nphase), grad_damage(nphase,3)
        REAL(kind=AbqRK), INTENT(IN) :: parHFEDMatrix(nHFEDpar), parIEDMatrix(nIEDpar) 
        REAL(kind=AbqRK), INTENT(IN) :: prop_df_a, prop_df_b
        REAL(kind=AbqRK), INTENT(IN) :: H
        REAL(kind=AbqRK) :: d_BKED_d_eps_d_eps(3,3,3,3)
        REAL(kind=AbqRK) :: temp1(3), temp2(3,3)
        INTEGER(kind=AbqIK) :: i1, i2, pos_i1, pos_i2
        DIMENSION tangent(ntens,ntens)

        d_BKED_d_eps_d_eps = zero

        ! second derivatives of bulk energy density
        d_BKED_d_eps_d_eps = d_bulkED_d_eps_d_eps(eps,damage(1),nHFEDpar,parHFEDMatrix,prop_df_a,prop_df_b)

        tangent=zero
        !
        ! eps-eps: contributions of BKED
        tangent(1,1) = d_BKED_d_eps_d_eps(1,1,1,1)
        tangent(1,2) = d_BKED_d_eps_d_eps(1,1,2,2)
        tangent(1,3) = d_BKED_d_eps_d_eps(1,1,3,3)
        tangent(2,1) = d_BKED_d_eps_d_eps(2,2,1,1)
        tangent(2,2) = d_BKED_d_eps_d_eps(2,2,2,2)
        tangent(2,3) = d_BKED_d_eps_d_eps(2,2,3,3)
        tangent(3,1) = d_BKED_d_eps_d_eps(3,3,1,1)
        tangent(3,2) = d_BKED_d_eps_d_eps(3,3,2,2)
        tangent(3,3) = d_BKED_d_eps_d_eps(3,3,3,3)
        IF (D .GE. 2) THEN
          tangent(1,4) = d_BKED_d_eps_d_eps(1,1,1,2)
          tangent(2,4) = d_BKED_d_eps_d_eps(2,2,1,2)
          tangent(3,4) = d_BKED_d_eps_d_eps(3,3,1,2)
          tangent(4,1) = d_BKED_d_eps_d_eps(1,2,1,1)
          tangent(4,2) = d_BKED_d_eps_d_eps(1,2,2,2)
          tangent(4,3) = d_BKED_d_eps_d_eps(1,2,3,3)
          tangent(4,4) = d_BKED_d_eps_d_eps(1,2,1,2)
        END IF
        IF (D .EQ. 3) THEN
          tangent(1,5) = d_BKED_d_eps_d_eps(1,1,1,3)
          tangent(1,6) = d_BKED_d_eps_d_eps(1,1,2,3)
          tangent(2,5) = d_BKED_d_eps_d_eps(2,2,1,3)
          tangent(2,6) = d_BKED_d_eps_d_eps(2,2,2,3)
          tangent(3,5) = d_BKED_d_eps_d_eps(3,3,1,3)
          tangent(3,6) = d_BKED_d_eps_d_eps(3,3,2,3)
          tangent(4,5) = d_BKED_d_eps_d_eps(1,2,1,3)
          tangent(4,6) = d_BKED_d_eps_d_eps(1,2,2,3)
          tangent(5,1) = d_BKED_d_eps_d_eps(1,3,1,1)
          tangent(5,2) = d_BKED_d_eps_d_eps(1,3,2,2)
          tangent(5,3) = d_BKED_d_eps_d_eps(1,3,3,3)
          tangent(5,4) = d_BKED_d_eps_d_eps(1,3,1,2)
          tangent(5,5) = d_BKED_d_eps_d_eps(1,3,1,3)
          tangent(5,6) = d_BKED_d_eps_d_eps(1,3,2,3)
          tangent(6,1) = d_BKED_d_eps_d_eps(2,3,1,1)
          tangent(6,2) = d_BKED_d_eps_d_eps(2,3,2,2)
          tangent(6,3) = d_BKED_d_eps_d_eps(2,3,3,3)
          tangent(6,4) = d_BKED_d_eps_d_eps(2,3,1,2)
          tangent(6,5) = d_BKED_d_eps_d_eps(2,3,1,3)
          tangent(6,6) = d_BKED_d_eps_d_eps(2,3,2,3)
        END IF
        !
        ! damage parameters
        !
        DO i1=1,nphase
          pos_i1 = pos_p(i1)
          !
          ! eps-damage: contributions of ...
          tangent(pos_i1,1) = zero ! dummy
          tangent(pos_i1,2) = zero ! dummy
          tangent(pos_i1,3) = zero ! dummy
          IF (D .GE. 2) THEN
            tangent(pos_i1,4) = zero ! dummy
          END IF
          IF (D .EQ. 2) THEN
            tangent(1:4,pos_i1) = tangent(pos_i1,1:4)
          ELSE IF (D .EQ. 3) THEN
            tangent(pos_i1,5) = zero ! dummy
            tangent(pos_i1,6) = zero ! dummy
            tangent(1:6,pos_i1) = tangent(pos_i1,1:6)
          END IF
          !
          ! eps-grad_damage: contributions of ...
          temp1(:) = zero ! dummy
          tangent(pos_i1+1:pos_i1+D,1) = temp1(1:D)
          temp1(:) = zero ! dummy
          tangent(pos_i1+1:pos_i1+D,2) = temp1(1:D)
          temp1(:) = zero ! dummy
          tangent(pos_i1+1:pos_i1+D,3) = temp1(1:D)
          IF (D .GE. 2) THEN
            temp1(:) = zero ! dummy
            tangent(pos_i1+1:pos_i1+D,4) = temp1(1:D)
          END IF
          IF (D .EQ. 2) THEN
            tangent(1:4,pos_i1+1:pos_i1+D) = TRANSPOSE(tangent(pos_i1+1:pos_i1+D,1:4))
          ELSE IF (D .EQ. 3) THEN
            temp1(:) = zero ! dummy
            tangent(pos_i1+1:pos_i1+D,5) = temp1(1:D)
            temp1(:) = zero ! dummy
            tangent(pos_i1+1:pos_i1+D,6) = temp1(1:D)
            tangent(1:6,pos_i1+1:pos_i1+D) = TRANSPOSE(tangent(pos_i1+1:pos_i1+D,1:6))
          END IF
          !
          DO i2=1,nphase
            pos_i2 = pos_p(i2)
            !
            ! damage-damage: contributions of interface energy
            tangent(pos_i1,pos_i2) = zero ! dummy
            IF (pos_i1 .EQ. pos_i2) THEN 
              tangent(pos_i1,pos_i1) = zero ! dummy
              IF (i1 .EQ. 1) THEN 
                tangent(pos_i1,pos_i1) = d_bulkED_d_damage_d_damage_H(eps,damage(1),nHFEDpar,parHFEDMatrix,prop_df_a,prop_df_b,H) &
									   + d_IED_d_damage_d_damage(D,damage(1),grad_damage(1,1:D),nIEDpar,parIEDMatrix)
              END IF
            END IF
            tangent(pos_i2,pos_i1) = tangent(pos_i1,pos_i2)
            !
            ! damage-grad_damage: contributions of ...
            temp1(:) = zero ! dummy
            tangent(pos_i1,pos_i2+1:pos_i2+D) = temp1(1:D)
            tangent(pos_i2+1:pos_i2+D,pos_i1) = tangent(pos_i1,pos_i2+1:pos_i2+D)
            !
            ! grad_damage-grad_damage: contributions of interface energy
            temp2(:,:) = zero ! dummy
            IF (pos_i1 .EQ. pos_i2) THEN
              IF (i1 .EQ. 1) THEN
                temp2(1:D,1:D) = d_IED_d_grad_damage_d_grad_damage(D,damage(1),grad_damage(1,1:D),nIEDpar,parIEDMatrix)
              END IF
            END IF
            tangent(pos_i1+1:pos_i1+D,pos_i2+1:pos_i2+D) = temp2(1:D,1:D)
            tangent(pos_i2+1:pos_i2+D,pos_i1+1:pos_i1+D) = TRANSPOSE(tangent(pos_i1+1:pos_i1+D,pos_i2+1:pos_i2+D))
            !
          END DO
          !
        END DO
        !
      END FUNCTION tangent

!------------------------------------------------------------------------------------

    END SUBROUTINE umatPF

END MODULE PhaseField_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

