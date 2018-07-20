!> \file micm_type_defs.f90
!!  Contains type definitions for MICM variables and physics-related variables

module micm_type_defs

! use machine, only: kind_phys

! integer, parameter :: kind_phys=8

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! The following definition sets up the variables for use within ibox
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

! Filter with CPP for PGI compiler
#ifndef __PGI
!> \section arg_table_micm_data_type
!! | local_name | standard_name                                    | long_name                               | units       | rank | type      |    kind   | intent | optional |
!! |------------|--------------------------------------------------|-----------------------------------------|-------------|------|-----------|-----------|--------|----------|
!! | my_co(:)   | my_volume_mixing_ratio_co                        | volume mixing ratio co                  | mole mole-1 |    1 | real      | kind_phys | none   | F        |
!! | my_o3(:)   | my_volume_mixing_ratio_o3                        | volume mixing ratio o3                  | mole mole-1 |    1 | real      | kind_phys | none   | F        |
!! | errmsg     | error_message                                    | CCPP error message                      | none        |    0 | character | len=*     | out    | F        | 
!! | errflg     | error_flag                                       | CCPP error flag                         | flag        |    0 | integer   |           | out    | F        |
!!

#endif
  type micm_data_type

!    real, allocatable      :: my_co(:)
!    real                   :: dt
!    integer                :: ncol
!    integer                :: nlev
!    character(len=512)     :: errmsg
!    integer                :: errflg

    contains
!      procedure :: create => physics_create
!      procedure :: associate => physics_associate
  end type micm_data_type


end module micm_type_defs
