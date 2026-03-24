MODULE AliasModuleCZ

!  Alias Modul
!  
!  Speichert alle Umbenennungen, damit sie von verschiedenen Subroutinen aufgerufen werden können

   ! ==========================================================
   ! Penalty Energy
   ! ==========================================================
   USE PenaltyEnergyModule, ONLY: &
        PenED  => PenED_Exponential, &
        d_PenED_d_damageJump => d_PenED_d_damageJump_Exponential, &
        d_PenED_d_damageJump_d_damageJump => &
        d_PenED_d_damageJump_d_damageJump_Exponential


   ! ==========================================================
   ! Interface Energy (AT2 aktiv)
   ! Zum Wechseln auf AT1 einfach hier umstellen
   ! ==========================================================
   USE InterfaceEnergyModule, ONLY: &
   
        IED  => IED_AT2, &
        d_IED_d_damage => d_IED_AT2_d_damage, &
        d_IED_d_damage_d_damage => d_IED_AT2_d_damage_d_damage, &
        d_IED_d_grad_damage => d_IED_AT2_d_grad_damage, &
        d_IED_d_grad_damage_d_grad_damage => &
        d_IED_AT2_d_grad_damage_d_grad_damage, &
        d_IED_d_damage_d_grad_damage => &
        d_IED_AT2_d_damage_d_grad_damage
   
   
!~         IED  => IED_AT2, &
!~         d_IED_d_phase => d_IED_AT2_d_phase, &
!~         d_IED_d_phase_d_phase => d_IED_AT2_d_phase_d_phase, &
!~         d_IED_d_grad_phase => d_IED_AT2_d_grad_phase, &
!~         d_IED_d_grad_phase_d_grad_phase => &
!~         d_IED_AT2_d_grad_phase_d_grad_phase, &
!~         d_IED_d_phase_d_grad_phase => &
!~         d_IED_AT2_d_phase_d_grad_phase

!~         IED  => IED_AT1, &
!~         d_IED_d_phase => d_IED_AT1_d_phase, &
!~         d_IED_d_phase_d_phase => d_IED_AT1_d_phase_d_phase, &
!~         d_IED_d_grad_phase => d_IED_AT1_d_grad_phase, &
!~         d_IED_d_grad_phase_d_grad_phase => &
!~         d_IED_AT1_d_grad_phase_d_grad_phase, &
!~         d_IED_d_phase_d_grad_phase => &
!~         d_IED_AT1_d_phase_d_grad_phase

   ! ==========================================================
   ! Viscous Dissipation CZ
   ! ==========================================================
!~    USE ViscousDissipationCZModule, ONLY: &
!~         VDED  => VDED_Liu, &
!~         d_VDED_d_phase => d_VDED_Liu_d_phase, &
!~         d_VDED_d_phase_d_phase => d_d_VDED_Liu_d_phase_d_phase


!~    ! ==========================================================
!~    ! Degradation Function
!~    ! ==========================================================
  USE DegradationModule, ONLY: Degradation => DegradationCubic, &
                                              d_Degradation_d_damage => d_DegradationCubic_d_damage, &
                                              d_Degradation_d_damage_d_damage => d_DegradationCubic_d_damage_d_damage 


   ! ==========================================================
   ! Split Energy Module
   ! ==========================================================
                                               
                                               
!~    ! ==========================================================
!~    ! Split Energy Module
!~    ! ==========================================================
!~    USE SplitIntEnergyModule, ONLY: CZEDpos => CZEDposLiuMixedMode, &
!~                                                d_CZEDpos_d_sep => d_CZEDposLiuMixedMode_d_sep, &
!~                                                d_CZEDpos_d_sep_d_sep => d_CZEDposLiuMixedMode_d_sep_d_sep, &
!~                                                CZEDneg => CZEDnegLiuMixedMode, &
!~                                                d_CZEDneg_d_sep => d_CZEDnegLiuMixedMode_d_sep, &
!~                                                d_CZEDneg_d_sep_d_sep => d_CZEDnegLiuMixedMode_d_sep_d_sep, &
!~                                                Gieff, &
!~                                                CorrectionFactor, &
!~                                                Hieff, &
!~                                                CZEDtens_H_sph, &
!~                                                CZEDtens_H_dev
                                               

				   
   IMPLICIT NONE

   PUBLIC :: PenED, d_PenED_d_damageJump, d_PenED_d_damageJump_d_damageJump
   PUBLIC :: IED, d_IED_d_damage, d_IED_d_damage_d_damage
!~    PUBLIC :: VDED, d_VDED_d_phase, d_VDED_d_phase_d_phase
   PUBLIC :: Degradation, d_Degradation_d_damage, d_Degradation_d_damage_d_damage
!~    PUBLIC :: CZEDpos, d_CZEDpos_d_sep, d_CZEDpos_d_sep_d_sep, CZEDneg, d_CZEDneg_d_sep, d_CZEDneg_d_sep_d_sep
!~    PUBLIC :: Gieff, CorrectionFactor, Hieff, CZEDtens_H_sph, CZEDtens_H_dev
   

END MODULE AliasModuleCZ
