MODULE AliasModulePF

!  Alias Modul
!  
!  Speichert alle Umbenennungen, damit sie von verschiedenen Subroutinen aufgerufen werden können

   ! ==========================================================
   ! Viscous Dissipation PF
   ! ==========================================================
!~       USE ViscousDissipationModule, ONLY: VDED => VDED_Liu, &
!~                                                d_VDED_d_phase => d_VDED_Liu_d_phase, &
!~                                                d_d_VDED_d_phase_d_phase => d_d_VDED_Liu_d_phase_d_phase

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
        
        
   ! ==========================================================
   ! Degradation Function
   ! ==========================================================

		USE DegradationModule, ONLY: Degradation => DegradationCubic, &
                                              d_Degradation_d_damage => d_DegradationCubic_d_damage, &
                                              d_Degradation_d_damage_d_damage => d_DegradationCubic_d_damage_d_damage 
                                               
   ! ==========================================================
   ! Energy Split
   ! ==========================================================
   
   
  USE SplitModulePF, ONLY: HFEDpos => HFEDposNoSplit, &
                                               d_HFEDpos_d_eps_e => d_HFEDposNoSplit_d_eps_e, &
                                               d_HFEDpos_d_eps_e_d_eps_e => d_HFEDposNoSplit_d_eps_e_d_eps_e, &
                                               HFEDneg => HFEDnegNoSplit, &
                                               d_HFEDneg_d_eps_e => d_HFEDnegNoSplit_d_eps_e, &
                                               d_HFEDneg_d_eps_e_d_eps_e => d_HFEDnegNoSplit_d_eps_e_d_eps_e




        
   IMPLICIT NONE
   
!~    PUBLIC :: VDED, d_VDED_d_phase, d_d_VDED_d_phase_d_phase
   PUBLIC :: IED, d_IED_d_damage, d_IED_d_damage_d_damage
   PUBLIC :: Degradation, d_Degradation_d_damage, d_Degradation_d_damage_d_damage
   PUBLIC :: HFEDpos, d_HFEDpos_d_eps_e, d_HFEDpos_d_eps_e_d_eps_e, HFEDneg, d_HFEDneg_d_eps_e, d_HFEDneg_d_eps_e_d_eps_e
   

END MODULE AliasModulePF
