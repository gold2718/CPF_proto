module micm_type_defs

use machine, only: kind_phys


!> \section arg_table_micm_data_type  Argument Table
!! | local_name | standard_name                                    | long_name                               | units       | rank | type      | kind      | intent | optional |
!! |------------|--------------------------------------------------|-----------------------------------------|-------------|------|-----------|-----------|--------|----------|
!! | tune_factor| tuning factor for solver                         | tuning factor for solver                | mole mole-1 |    0 | real      | kind_phys | out    | F        |
!! | T          | air_temperature                                  | temperature                             | K           |    0 | real      | kind_phys | in     | F        |
!! | k_rateConst| k_rate_constants                                 | k Rate Constants                        | none        |    1 | real      | kind_phys | inout  | F        |
!!

real(kind_phys)  :: tune_factor
real (kind_phys) :: T
real, pointer    :: k_rateConst(:)

type micm_data_type
end type micm_data_type

end module micm_type_defs
