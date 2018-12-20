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
!! | rho        | dry_air_density  | dry air density      | kg/m^3|    2 | real      | kind_phys | in     | F        |
!! | z          | height           | height               | m     |    2 | real      | kind_phys | in     | F        |
!! | pk         | exner_function   | exner function       | none  |    2 | real      | kind_phys | in     | F        |
!! | latvap     | latent heat of vaporization                      | latent heat of vaporization             | J/kg        |    0 | real      | kind_phys | in     | F        |
!! | pstd       | reference pressure at sea level                  | reference pressure at sea level         | mb          |    0 | real      | kind_phys | in     | F        |
!! | rhoh2o     | density of liquid water                          | density of liquid water                 | kg/m^3      |    0 | real      | kind_phys | in     | F        |
!! | pcols      | maximum_horizontal_extent                        | number of columns                       | none        |    0 | integer   |           | in     | F        |
!! | pver       | maximum_vertical_extent                                      | maximum_vertical_exten      | none        |    0 | integer   |           | in     | F        |
!! | ncol       | horizontal_extent                                            | number of columns           | none        |    0 | integer   |           | in     | F        |
!! | pmid       | pressure_at_mid_level                                        | pressure_at_mid_level                   | hPa         |    2 | real      | kind_phys | in     | F        |
!! | rair       | dry_air_gas_constant                                         | dry_air_gas_constant                    | J/K/kg      |    0 | real      | kind_phys | in     | F        |
!! | cpair      | specific_heat_of_dry_air                                     | specific_heat_of_dry_air                | J/K/kg      |    0 | real      | kind_phys | in     | F        |
!! | ztodt      | physics_time_step                                            | physics_time_step                       | s           |    0 | real      | kind_phys | in     | F        |
!! | th         | air_potential_temperature                                    | air_potential_temperature               | K           |    2 | real      | kind_phys | in     | F        |
!! | qv         | mass_fraction_of_water_in_air                                | mass_fraction_of_water_in_air           | kg kg-1     |    2 | real      | kind_phys | in     | F        |
!! | qc         | mass_fraction_of_cloud_liquid_water_in_air                   | cloud_liquid_water_in_air               | kg kg-1     |    2 | real      | kind_phys | in     | F        |
!! | qr         | mass_fraction_of_rain_in_air                                 | rain in air                             | kg kg-1     |    2 | real      | kind_phys | in     | F        |
!! | precl      | precipitation                                                | precipitation                           | m/s         |    1 | real      | kind_phys | inout  | F        |
!! | t_prev     | air_temperature_prev                                         | previous air temperature                | K           |    2 | real      | kind_phys | inout  | F        |
!! | qv_prev    | mass_fraction_of_water_in_air_prev                           | previous water in air                   | kg kg-1     |    2 | real      | kind_phys | inout  | F        |
!! | qc_prev    | mass_fraction_of_cloud_liquid_water_in_air_prev              | previous cloud_liquid_water_in_air      | kg kg-1     |    2 | real      | kind_phys | inout  | F        |
!! | qr_prev    | mass_fraction_of_rain_in_air_prev                            | previous rain in air                    | kg kg-1     |    2 | real      | kind_phys | inout  | F        |
!! | precl_prev | precipitation_prev                                           | previous precipitation                  | m/s         |    1 | real      | kind_phys | inout  | F        |
!! | ttend_t    | total_tendency_of_air_temperature                            | total_tendency_of_air_temperature       | K           |    2 | real      | kind_phys | inout  | F        |
!! | ttend_qv   | total_tendency_of_mass_fraction_of_water_in_air              | total tendency of_water in air          | kg kg-1     |    2 | real      | kind_phys | inout  | F        |
!! | ttend_qc   | total_tendency_of_mass_fraction_of_cloud_liquid_water_in_air | total tendency of_water in air          | kg kg-1     |    2 | real      | kind_phys | inout  | F        |
!! | ttend_qr   | total_tendency_of_mass_fraction_of_rain_in_air               | total tendency of rain in air           | kg kg-1     |    2 | real      | kind_phys | inout  | F        |
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
