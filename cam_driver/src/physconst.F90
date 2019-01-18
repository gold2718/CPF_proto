MODULE physconst

   implicit none
   public

   integer, parameter :: R8 = selected_real_kind(12) ! 8 byte real

   !----------------------------------------------------------------------------
   ! physical constants (all data public)
   !----------------------------------------------------------------------------

   real(R8),parameter :: PI        = 3.14159265358979323846_R8 ! pi
   real(R8),parameter :: AVOGADRO  = 6.02214e26_R8             ! Avogadro's number ~ molecules/kmole
   real(R8),parameter :: BOLTZMANN = 1.38065e-23_R8            ! Boltzmann's constant ~ J/K/molecule
   real(R8),parameter :: CPDAIR    = 1.00464e3_R8              ! specific heat of dry air   ~ J/kg/K
   real(R8),parameter :: RGAS      = AVOGADRO*BOLTZMANN        ! Universal gas constant ~ J/K/kmole
   real(R8),parameter :: MWDAIR    = 28.966_R8                 ! molecular weight dry air ~ kg/kmole
   real(R8),parameter :: RDAIR     = RGAS/MWDAIR               ! Dry air gas constant     ~ J/K/kg
   real(R8),parameter :: LATVAP    = 2.501e6_R8                ! latent heat of evaporation ~ J/kg
   real(R8),parameter :: PSTD      = 101325.0_R8               ! standard pressure ~ pascals
   real(R8),parameter :: RHOFW     = 1.000e3_R8                ! density of fresh water     ~ kg/m^3

!-----------------------------------------------------------------------------

END MODULE physconst
