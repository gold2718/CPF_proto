!> \file cam_type_defs.f90
!!  Contains type definitions for cam variables and physics-related variables

module cam_vardefs

use machine, only: kind_phys

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! The following definition sets up the variables for use within cam
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! Filter with CPP for PGI compiler
#ifndef __PGI
!> \section arg_table_cam_vardefs
!! | local_name             | standard_name               | long_name               | units       | rank | type      |    kind   | intent | optional |
!! |------------------------|-----------------------------|-------------------------|-------------|------|-----------|-----------|--------|----------|
!! | my_co(:)               | my_volume_mixing_ratio_co   | volume mixing ratio co  | mole mole-1 |    1 | real      | kind_phys | none   | F        |
!! | my_o3(:)               | my_volume_mixing_ratio_o3   | volume mixing ratio o3  | mole mole-1 |    1 | real      | kind_phys | none   | F        |
!! | state_host%Temperature | air_temperature             | temperature             | K           |    0 | real      | kind_phys | none   | F        | 
!! | errmsg                 | error_message               | CCPP error message      | none        |    0 | character | len=*     | out    | F        | 
!! | errflg                 | error_flag                  | CCPP error flag         | flag        |    0 | integer   |           | out    | F        |
!! | tune_factor            | tuning factor for solver    | tuning factor for solver| mole mole-1 |    0 | real      | kind_phys | none   | F        |
!! | k_rateConst            | k_rate_constants            | k Rate Constants        | none        |    1 | real      | kind_phys | none   | F        |
!! | ncol       | number_of_columns| number of columns    | none  |    0 | integer   |           | in     | F        |
!! | pver       | number_of_vert   | number of vert levels| none  |    0 | integer   |           | in     | F        |
!! | ztodt      | time_step        | time step            | s     |    0 | real      | kind_phys | in     | F        |
!! | rho        | dry_air_density  | dry air density      | kg/m^3|    2 | real      | kind_phys | in     | F        |
!! | z          | height           | height               | m     |    2 | real      | kind_phys | in     | F        |
!! | pk         | exner_function   | exner function       | none  |    2 | real      | kind_phys | in     | F        |
!! | th         | potential_temp   | potential temp       | K     |    2 | real      | kind_phys | inout  | F        |
!! | qv         | water_vapor      | water vapor          | gm/gm |    2 | real      | kind_phys | inout  | F        |
!! | qc         | cld water_vapor  | cld water vapor      | gm/gm |    2 | real      | kind_phys | inout  | F        |
!! | qr         | rain water_vapor | rain water vapor     | gm/gm |    2 | real      | kind_phys | inout  | F        |
!! | precl      | precipitation    | precipitation        | m/s   |    1 | real      | kind_phys | out    | F        |
!! | rair       | gas constant for dry air                         | gas constant for dry air                | J/(kgK)     |    0 | real      | kind_phys | in     | F        |
!! | cpair      | heat capacity at constant pres                   | heat capacity at constant pres          | J/(kgK)     |    0 | real      | kind_phys | in     | F        |
!! | latvap     | latent heat of vaporization                      | latent heat of vaporization             | J/kg        |    0 | real      | kind_phys | in     | F        |
!! | pstd       | reference pressure at sea level                  | reference pressure at sea level         | mb          |    0 | real      | kind_phys | in     | F        |
!! | rhoh2o     | density of liquid water                          | density of liquid water                 | kg/m^3      |    0 | real      | kind_phys | in     | F        |
!!

#endif

  real (kind=kind_phys), pointer :: my_co(:)
  real (kind=kind_phys), pointer :: my_o3(:)
  character(len=512)     :: errmsg
  integer :: errflg

  type my_state
    real(kind_phys) :: Temperature
  end type my_state
  type(my_state), target :: state_host

  real(kind_phys)  :: tune_factor
  real, pointer    :: k_rateConst(:)

end module cam_vardefs
