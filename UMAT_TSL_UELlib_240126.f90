!DEC$ FREEFORM
!FORTRAN90 mit FREEFORM
!
! Compile with IFORT Version 11
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Stephan Roth, TU Bergakademie Freiberg, 33.01.2024
!
! 06.01.2016: Modifikationen fuer UELlib-Kompatibilitaet
! 25.01.2024: Strafenergie um Sprung in Schaedigungsvariable zu minimieren,
!             konstante Kontaktsteifigkeit (linear),
!             Grenzflaechenenergie
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SAVED VARIABLES AND INPUT PARAMETERS
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! SVR: svr(1)     : damage variable 0<=D<=1
!      svr(2)     : normalized effective separation, lambda
!      svr(3)     : damage jump
!      svr(4)     : not used
!      svr(5)     : total potential density
!      svr(6)     : stored energy density
!      svr(7)     : interface energy density
!      svr(8)     : contact energy density
!      svr(9)     : penalty energy density
!      svr(10-12) : separation (local frame)
!      svr(13-15) : separation (global frame)
!      svr(16-18) : damage gradient (local frame)
!      svr(19-21) : damage gradient (global frame)
!      svr(22-24) : traction vector (local frame)
!      svr(25-30) : traction vector (global frame)
!
! material parameters: thicknessIndex -- plane strain thickness
!                       1 -- internal length l_0
!                       2 -- cohesive strength t_0
!                       3 -- cohesive length delta_0
!                       4 -- flag: numerical tangent (inactive)
!                       5 -- penalty factor
!                       6 -- critical interface fracture energy density
!                       7 -- plane strain thickness 
!                       8 -- factor to scale stiffness for negative normal separation
!                       9 -- parameter 1 of degradation function
!                      10 -- parameter 2 of degradation function
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE CZUMAT_module

  IMPLICIT NONE

  PUBLIC :: CheckMaterialParameters, umatCZM

  CONTAINS

!------------------------------------------------------------------------------------

    SUBROUTINE CheckMaterialParameters(props)
      ! Check of all material parameters

      USE ABQINTERFACE
      USE FLOATNUMBERS
      USE ABAModul
      USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex

      IMPLICIT NONE
      REAL(kind=AbqRK), INTENT(IN) :: props(numMatPar)
      INTEGER(kind=AbqIK) :: num_counter
      REAL(kind=AbqRK) :: prop_t0, prop_thickness

      prop_thickness = one

      ! read material parameters from props
      prop_t0      = props( 2) ! cohesive strength
      IF (props(thicknessIndex) .NE. zero) prop_thickness = props(thicknessIndex) ! thickness of 2D plane strain element

      ! parameter check
      IF (prop_t0 .LE. zero) THEN
        WRITE(7,*) 'Cohesive strength should exceed zero. EXIT', prop_t0
        CALL XEXIT()
      END IF
      IF (prop_thickness .LT. zero) THEN
        WRITE(7,*) 'Thickness of 2D element should exceed zero. EXIT', prop_thickness
        CALL XEXIT()
      END IF

    END SUBROUTINE CheckMaterialParameters

!------------------------------------------------------------------------------------

    SUBROUTINE umatCZM(stress,svr,Ct,energy_elast,dissipat_plast,dissipat_creep,rpl,ddsddt,drplde,drpldt, &
                       stran,dstran,time,dtime,Temp,dTemp,predef_umat,dpred,cmname,ndi,nshr,ntens, &
                       num_svr,props_mat,num_matpar,coords_gp,Trafomat,pnewdt,intlen,F0,F1,jelem, &
                       npt,layer,kspt,kstep,kinc)
      ! compute cohesive tractions and cohesive stiffness

      USE ABQINTERFACE
      USE FLOATNUMBERS
      USE ABAModul
      USE PhaseField_Parameters, ONLY: numSDV, numMatPar, thicknessIndex
      
      USE AliasModuleCZ
      USE CohesiveEnergyModule
      USE ContactEnergyModule
      

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
      INTEGER(kind=AbqIK) :: D, i1, i2, pos_damage, pos_damageJump, pos_damageGradient, nIEDpar, nCEDpar, normalDirectionIndex
      REAL(kind=AbqRK) :: damage, damageJump, damage_old
      REAL(kind=AbqRK) :: lambda, lambda_old
      REAL(kind=AbqRK) :: prop_l0, prop_t0, prop_delta0, prop_G0, prop_pen, prop_compressionstiffness, prop_df_alpha, prop_df_beta
      REAL(kind=AbqRK) :: totalPotential, storedEnergy, interfaceEnergy, penaltyEnergy, contactEnergy, transformedStress(6)
      REAL(kind=AbqRK), ALLOCATABLE :: sep(:), sepGlobal(:), trac(:), damageGradient(:), damageGradientLocal(:), damageGradientGlobal(:), &
                                       K(:,:), tmp(:), tmp2(:,:), parIEDMatrix(:), parCEDMatrix(:)

      ! initialization
      ! add, if needed
	  !
      ! read material parameters from props_mat
      prop_l0      = props_mat( 1) ! internal length
      prop_t0      = props_mat( 2) ! cohesive strength
      prop_delta0  = props_mat( 3) ! cohesive length
      prop_pen     = props_mat( 5) ! penalty factor
      prop_G0      = props_mat( 6) ! critical interface fracture energy density
      prop_compressionstiffness = props_mat( 8) ! penalty factor for negative normal separation
      prop_df_alpha  = props_mat( 9) ! parameter 1 of degradation function
      prop_df_beta  = props_mat(10) ! parameter 2 of degradation function
	  !
      ! arrays of material parameters
      !
      ! interface energy density -- see InterfaceEnergyModule.f90
      nIEDpar = 2
      ALLOCATE(parIEDMatrix(nIEDpar))
      parIEDMatrix(1) = prop_G0
      parIEDMatrix(2) = prop_l0
      ! contact energy density -- see ContactEnergyModule.f90
      nCEDpar = 3
      ALLOCATE(parCEDMatrix(nCEDpar))
      parCEDMatrix(1) = prop_delta0
      parCEDMatrix(2) = prop_t0
      parCEDMatrix(3) = prop_compressionstiffness

      ! model dimension
      D = ndi+nshr
      ! index of normal direction
      normalDirectionIndex = D

      ALLOCATE(sep(D), sepGlobal(D), trac(D), K(D,D), tmp(D), tmp2(D,D), damageGradient(D-1), damageGradientLocal(D), damageGradientGlobal(D))
      sep=zero; sepGlobal=zero; trac=zero; K=zero; tmp=zero; tmp2=zero; damageGradient=zero; damageGradientLocal=zero; damageGradientGlobal=zero

      ! position in generalized arrays
      pos_damage         = ntens-(D+1)+1
      pos_damageJump     = ntens-(D+1)+2
      pos_damageGradient = ntens-(D+1)+3


      ! generalized kinematic measures
      !
      ! separation vector coordinates
      sep(1:D) = stran(1:D)
      ! damage variable
      damage = stran(pos_damage)
      ! damage jump
      damageJump = stran(pos_damageJump)
      ! damage gradient
      damageGradient = stran(pos_damageGradient:pos_damageGradient-2+D)
      !
      ! normalised effective separation
      lambda = effSep(D,normalDirectionIndex,sep,prop_delta0)
      ! damage gradient array of size D
      damageGradientLocal(1:D-1) = damageGradient(1:D-1)


      ! saved variables of last converged increment
      !
      damage_old = svr(1)
      lambda_old = svr(2)


      ! energetic quantities
      !
      storedEnergy    = CZED(D,normalDirectionIndex,sep,damage,prop_delta0,prop_t0,prop_df_alpha,prop_df_beta)
      interfaceEnergy = IED(D-1,damage,damageGradient,nIEDpar,parIEDMatrix) 
      penaltyEnergy   = PenED(damageJump,prop_pen)
      contactEnergy   = ContactED(D,normalDirectionIndex,sep,nCEDpar,parCEDMatrix)
      !
      totalPotential  = storedEnergy + interfaceEnergy + penaltyEnergy + contactEnergy
      !
      ! assignments
      energy_elast   = storedEnergy
      !dissipat_plast = penaltyEnergy
      dissipat_plast = interfaceEnergy
      !dissipat_creep = totalPotential
      dissipat_creep = contactEnergy


      ! compute generalized stresses
      !
      stress = zero
      !
      ! local tractions
      trac = d_CZED_d_sep(D,normalDirectionIndex,damage,sep,prop_delta0,prop_t0,prop_df_alpha,prop_df_beta)
      ! consider negative normal separation / compression
      trac = tracWithContact(D,normalDirectionIndex,trac,sep,nCEDpar,parCEDMatrix)
      ! arrange in stress array
      stress(1:D) = trac(1:D)
      !
      ! damage conjugate -- damage driving force
      ! contribution from elastic free energy density and surface energy density
      stress(pos_damage) = d_CZED_d_damage(D,normalDirectionIndex,damage,sep,prop_delta0,prop_t0,prop_df_alpha,prop_df_beta) &
                         + d_IED_d_damage(D-1,damage,damageGradient,nIEDpar,parIEDMatrix)
      !
      ! damage jump conjugate
      ! contribution from penalty energy density: reduce damage jump to a minimum (approx. zero)
      stress(pos_damageJump) = d_PenED_d_damageJump(damageJump,prop_pen)
      !
      ! damage gradient conjugate
      ! contribution from surface energy density
      stress(pos_damageGradient:pos_damageGradient+D-2) = d_IED_d_grad_damage(D-1,damage,damageGradient,nIEDpar,parIEDMatrix)


      ! generalised tangent
      !
      Ct = zero
      !
      ! separation
      ! contribution from elastic free energy density
      K = d_CZED_d_sep_d_sep(D,normalDirectionIndex,damage,sep,prop_delta0,prop_t0,prop_df_alpha,prop_df_beta)
      ! consider negative normal separation / compression
      K = tangentWithContact(D,normalDirectionIndex,K,sep,nCEDpar,parCEDMatrix)
      ! arrange  in tangent array
      Ct(1:D,1:D) = K(1:D,1:D)
      !
      ! damage
      ! contribution from elastic free energy density and surface energy density
      Ct(pos_damage,pos_damage) = d_CZED_d_damage_d_damage(D,normalDirectionIndex,damage,sep,prop_delta0,prop_t0,prop_df_alpha,prop_df_beta) &
                                + d_IED_d_damage_d_damage(D-1,damage,damageGradient,nIEDpar,parIEDMatrix)
      !
      ! mixed submatrix
      ! contribution from elastic free energy density and surface energy density
      tmp(1:D) = d_CZED_d_sep_d_damage(D,normalDirectionIndex,damage,sep,prop_delta0,prop_t0,prop_df_alpha,prop_df_beta)
      Ct(1:D,pos_damage) = tmp(1:D)
      Ct(pos_damage,1:D) = Ct(1:D,pos_damage)
      !
      ! damage gradient
      ! contribution from surface energy density
      tmp2(1:D-1,1:D-1) = d_IED_d_grad_damage_d_grad_damage(D-1,damage,damageGradient,nIEDpar,parIEDMatrix)
      Ct(pos_damageGradient:pos_damageGradient+D-2,pos_damageGradient:pos_damageGradient+D-2) = tmp2(1:D-1,1:D-1)
      !
      ! damage jump
      ! contribution from penalty energy: reduce damage jump to a minimum (approx. zero)
      Ct(pos_damageJump,pos_damageJump) = d_PenED_d_damageJump_d_damageJump(damageJump,prop_pen)
      !
      ! no coupling between damage and damage gradient
      ! no coupling between damage jump and separation, damage, or damage gradient


      ! Transformations (change of vector basis)
      !
      ! separation in global frame, Q^T*s
      sepGlobal            = TransformedVectorLocalToGlobal(D,sep,TrafoMat)
      ! damage gradient in global frame, Q^T*grad_d
      damageGradientGlobal = TransformedVectorLocalToGlobal(D,damageGradientLocal,TrafoMat)
      ! traction in global frame, Q^T*sigma*Q (stress tensor in Abaqus notation 3D: 11-22-33-12-13-23, 2D: 11-22-33-12)
      transformedStress    = TransformedTractionVectorLocalToStressArrayGlobal(D,trac,TrafoMat)


      ! save current state variables
      !
      svr(1)  = damage
      svr(2)  = lambda
      svr(3)  = damageJump
      !
      ! energy densities (per surface area)
      svr(5) = totalPotential
      svr(6) = storedEnergy
      svr(7) = interfaceEnergy
      svr(8) = contactEnergy
      svr(9) = penaltyEnergy
      !
      ! separation vector
      ! local frame
      svr(10:9+D)  = sep
      ! global frame
      svr(13:12+D) = sepGlobal
      !
      ! damage gradient vector
      ! local frame
      svr(16:15+D) = damageGradientLocal
      ! global frame
      svr(19:18+D) = damageGradientGlobal
      !
      ! cohesive traction vector
      ! local frame
      svr(22:21+D) = trac(1:D)
      ! global frame (stress tensor in Abaqus notation 3D: 11-22-33-12-13-23, 2D: 11-22-33-12)
      svr(25:30) = transformedStress

      ! dummy output
      rpl = zero; ddsddt = zero; drplde = zero; drpldt = zero

      DEALLOCATE(sep, sepGlobal, trac, K, tmp, tmp2, damageGradient, damageGradientLocal, damageGradientGlobal, parIEDMatrix, parCEDMatrix)

      CONTAINS

!------------------------------------------------------------------------------------

    END SUBROUTINE umatCZM

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION TransformedVectorLocalToGlobal(D,vec,TrafoMat)
      ! rotate/transform vector from local to global frame: v^glob = transposed(Q)*v^loc

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D
      REAL(kind=AbqRK), INTENT(IN) :: vec(D), TrafoMat(3,3)
      REAL(kind=AbqRK) :: TrafoMatD(D,D)
      DIMENSION TransformedVectorLocalToGlobal(D)

      TrafoMatD = Trafomat(1:D,1:D)

      TransformedVectorLocalToGlobal = zero
      ! rotate with local transformation matrix
      TransformedVectorLocalToGlobal = MATMUL(TRANSPOSE(TrafoMatD),vec)

    END FUNCTION TransformedVectorLocalToGlobal

!------------------------------------------------------------------------------------

    PURE REAL(kind=AbqRK) FUNCTION TransformedTractionVectorLocalToStressArrayGlobal(D,trac,TrafoMat)
      ! rotate/transform (stress valued) traction vector from local to global frame via stress matrix
      ! not verified

      USE ABQINTERFACE
      USE FLOATNUMBERS

      IMPLICIT NONE
      INTEGER(kind=AbqIK), INTENT(IN) :: D
      REAL(kind=AbqRK), INTENT(IN) :: trac(D), TrafoMat(3,3)
      REAL(kind=AbqRK) :: stressMatrix(3,3), transformedStress(3,3)
      DIMENSION TransformedTractionVectorLocalToStressArrayGlobal(6)

      stressMatrix = zero
      IF (D .EQ. 2) THEN
        stressMatrix(1,2) = trac(1)
        stressMatrix(2,2) = trac(2)
        stressMatrix(2,1) = stressMatrix(1,2)
      ELSE IF (D .EQ. 3) THEN
        stressMatrix(1,3) = trac(1)
        stressMatrix(2,3) = trac(2)
        stressMatrix(3,3) = trac(3)
        stressMatrix(3,1) = stressMatrix(1,3)
        stressMatrix(3,2) = stressMatrix(2,3)
      END IF

      transformedStress = MATMUL(TRANSPOSE(TrafoMat),MATMUL(stressMatrix,TrafoMat))
      
      ! reorder in array according to Abaqus convention: 3D: 11-22-33-12-13-23, 2D: 11-22-33-12
      TransformedTractionVectorLocalToStressArrayGlobal = zero
      TransformedTractionVectorLocalToStressArrayGlobal(1) = transformedStress(1,1)
      TransformedTractionVectorLocalToStressArrayGlobal(2) = transformedStress(2,2)
      TransformedTractionVectorLocalToStressArrayGlobal(3) = transformedStress(3,3)
      TransformedTractionVectorLocalToStressArrayGlobal(4) = transformedStress(1,2)
      IF (D .EQ. 3) THEN
        TransformedTractionVectorLocalToStressArrayGlobal(5) = transformedStress(1,3)
        TransformedTractionVectorLocalToStressArrayGlobal(6) = transformedStress(2,3)
      END IF

    END FUNCTION TransformedTractionVectorLocalToStressArrayGlobal

!------------------------------------------------------------------------------------

END MODULE CZUMAT_module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
