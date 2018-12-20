module kessler_update

use machine, only: r8 => kind_phys

implicit none
private
public :: kessler_update_init, kessler_update_run, kessler_update_finalize

contains

!> \section arg_table_kessler_update_init  Argument Table
!! | local_name | standard_name                                    | long_name                               | units       | rank | type      | kind      | intent | optional |
!! |------------|--------------------------------------------------|-----------------------------------------|-------------|------|-----------|-----------|--------|----------|
!! | pver       | maximum_vertical_extent                          | maximum_vertical_extent                 | none        |    0 | integer   |           | in     | F        |
!! | pcols      | maximum_horizontal_extent                        | number of columns                       | none        |    0 | integer   |           | in     | F        |
!! | t_prev     | air_temperature_prev                             | previous air temperature                | K           |    2 | real      | kind_phys | out    | F        |
!! | qc_prev    | mass_fraction_of_cloud_liquid_water_in_air_prev  | previous cloud_liquid_water_in_air      | kg kg-1     |    2 | real      | kind_phys | out    | F        |
!! | qr_prev    | mass_fraction_of_rain_in_air_prev                | previous rain in air                    | kg kg-1     |    2 | real      | kind_phys | out    | F        |
!! | qv_prev    | mass_fraction_of_water_in_air_prev               | previous water in air                   | kg kg-1     |    2 | real      | kind_phys | out    | F        |
!! | precl_prev | precipitation_prev                               | previous precipitation                  | m/s         |    1 | real      | kind_phys | out    | F        |
!! | errmsg     | error_message                                    | CCPP error message                      | none        |    0 | character | len=512   | out    | F        |
!! | errflg     | error_flag                                       | CCPP error flag                         | flag        |    0 | integer   |           | out    | F        |
!!

  subroutine kessler_update_init(pcols, pver, t_prev, qv_prev, qc_prev, qr_prev, precl_prev, errmsg, errflg)
    integer, intent(in)     :: pcols
    integer, intent(in)     :: pver
    real(r8), intent(out) :: t_prev(pcols,pver)
    real(r8), intent(out) :: qv_prev(pcols,pver)
    real(r8), intent(out) :: qc_prev(pcols,pver)
    real(r8), intent(out) :: qr_prev(pcols,pver)
    real(r8), intent(out) :: precl_prev(pcols)

    character(len=512),      intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg



!   Initialize t, qv, qc and qr previous to zero
    t_prev(:,:)     = 0._r8
    qv_prev(:,:)    = 0._r8
    qc_prev(:,:)    = 0._r8
    qr_prev(:,:)    = 0._r8
    precl_prev(:)   = 0._r8

    errmsg = ' '
    errflg = 0

  end subroutine kessler_update_init

!> \section arg_table_kessler_update_run  Argument Table
!! | local_name | standard_name                                                | long_name                               | units       | rank | type      | kind      | intent | optional |
!! |------------|--------------------------------------------------------------|-----------------------------------------|-------------|------|-----------|-----------|--------|----------|
!! | pcols      | maximum_horizontal_extent                                    | number of columns                       | none        |    0 | integer   |           | in     | F        |
!! | pver       | maximum_vertical_extent                                      | maximum_vertical_exten                  | none        |    0 | integer   |           | in     | F        |
!! | ncol       | horizontal_extent                                            | number of columns                       | none        |    0 | integer   |           | in     | F        |
!! | rair       | dry_air_gas_constant                                         | dry_air_gas_constant                    | J/K/kg      |    0 | real      | kind_phys | in     | F        |
!! | cpair      | specific_heat_of_dry_air                                     | specific_heat_of_dry_air                | J/K/kg      |    0 | real      | kind_phys | in     | F        |
!! | exner      | exner_function                                               | exner function                          | none        |    2 | real      | kind_phys | in     | F        |
!! | ztodt      | physics_time_step                                            | physics_time_step                       | s           |    0 | real      | kind_phys | in     | F        |
!! | th         | air_potential_temperature                                    | air_potential_temperature               | K           |    2 | real      | kind_phys | in     | F        |
!! | qv         | mass_fraction_of_water_in_air                                | mass_fraction_of_water_in_air           | kg kg-1     |    2 | real      | kind_phys | in     | F        |
!! | qc         | mass_fraction_of_cloud_liquid_water_in_air                   | cloud_liquid_water_in_air               | kg kg-1     |    2 | real      | kind_phys | in     | F        |
!! | qr         | mass_fraction_of_rain_in_air                                 | rain in air                             | kg kg-1     |    2 | real      | kind_phys | in     | F        |
!! | precl      | precipitation                                                | precipitation                           | m/s         |    1 | real      | kind_phys | in     | F        |
!! | t_prev     | air_temperature_prev                                         | previous air temperature                | K           |    2 | real      | kind_phys | inout  | F        |
!! | qv_prev    | mass_fraction_of_water_in_air_prev                           | previous water in air                   | kg kg-1     |    2 | real      | kind_phys | inout  | F        |
!! | qc_prev    | mass_fraction_of_cloud_liquid_water_in_air_prev              | previous cloud_liquid_water_in_air      | kg kg-1     |    2 | real      | kind_phys | inout  | F        |
!! | qr_prev    | mass_fraction_of_rain_in_air_prev                            | previous rain in air                    | kg kg-1     |    2 | real      | kind_phys | inout  | F        |
!! | precl_prev | precipitation_prev                                           | previous precipitation                  | m/s         |    1 | real      | kind_phys | inout  | F        |
!! | ttend_t    | total_tendency_of_air_temperature                            | total_tendency_of_air_temperature       | K           |    2 | real      | kind_phys | inout  | F        |
!! | ttend_qv   | total_tendency_of_mass_fraction_of_water_in_air              | total tendency of_water in air          | kg kg-1     |    2 | real      | kind_phys | inout  | F        |
!! | ttend_qc   | total_tendency_of_mass_fraction_of_cloud_liquid_water_in_air | total tendency of_water in air          | kg kg-1     |    2 | real      | kind_phys | inout  | F        |
!! | ttend_qr   | total_tendency_of_mass_fraction_of_rain_in_air               | total tendency of rain in air           | kg kg-1     |    2 | real      | kind_phys | inout  | F        |
!! | errmsg     | error_message                                                | CCPP error message                      | none        |    0 | character | len=512   | out    | F        |
!! | errflg     | error_flag                                                   | CCPP error flag                         | flag        |    0 | integer   |           | out    | F        |
!!
  subroutine kessler_update_run(pcols, pver, ncol, rair, cpair, exner, ztodt, th, qv, qc, qr, precl, t_prev, qv_prev,  qc_prev, qr_prev, precl_prev, ttend_t, ttend_qv, ttend_qc, ttend_qr, errmsg, errflg )

    integer, intent(in)     :: pcols
    integer, intent(in)     :: pver
    integer, intent(in)     :: ncol
    real(r8), intent(in)    :: rair
    real(r8), intent(in)    :: cpair
    real(r8), intent(in)    :: exner(pcols,pver)
    real(r8), intent(in)    :: ztodt
    real(r8), intent(in)    :: th(pcols,pver)   ! Potential temp.
    real(r8), intent(in)    :: qv(pcols,pver)
    real(r8), intent(in)    :: qc(pcols,pver)
    real(r8), intent(in)    :: qr(pcols,pver)
    real(r8), intent(in)    :: precl(pcols)

    real(r8), intent(inout) :: t_prev(pcols,pver)
    real(r8), intent(inout) :: qv_prev(pcols,pver)
    real(r8), intent(inout) :: qc_prev(pcols,pver)
    real(r8), intent(inout) :: qr_prev(pcols,pver)
    real(r8), intent(inout) :: precl_prev(pcols)

    real(r8), intent(inout) :: ttend_t(pcols,pver)
    real(r8), intent(inout) :: ttend_qc(pcols,pver)
    real(r8), intent(inout) :: ttend_qv(pcols,pver)
    real(r8), intent(inout) :: ttend_qr(pcols,pver)

    character(len=512),      intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

    integer                 :: k, rk

    errmsg = ' '
    errflg = 0

    ! Back out tendencies from updated fields
    do k = 1, pver
      rk = pver - k + 1
      ttend_t(:ncol,k) = (th(:ncol,rk)*exner(:ncol,k) - t_prev(:ncol,k)) * cpair / ztodt
      ttend_qv(:ncol,k) = (qv(:ncol,rk) - qv_prev(:ncol,rk)) / ztodt
      ttend_qc(:ncol,k) = (qc(:ncol,rk) - qc_prev(:ncol,rk)) / ztodt
      ttend_qr(:ncol,k) = (qr(:ncol,rk) - qr_prev(:ncol,rk)) / ztodt
    end do

    ! Update precip  -- CAC DO WE NEED TO DO THIS???
!    if (ncol < pcols) then
!      precl(ncol+1:pcols) = 0.0_r8
!      qc(ncol+1:pcols,:)  = 0.0_r8
!      qr(ncol+1:pcols,:)  = 0.0_r8
!    end if

!    surf_state%precl(:ncol) = surf_state%precl(:ncol) + precl(:ncol)  -- CAC - what do I do here?

    ! Set the previous q values to the current q
    t_prev(:,:)     = th(:ncol,:)*exner(:ncol,:)  ! NOT SURE WHAT TO PUT HERE -- CAC -- PROBABLY NOT CORRECT
    qv_prev(:,:)    = qv(:,:)
    qc_prev(:,:)    = qc(:,:)
    qr_prev(:,:)    = qr(:,:)
!     precl_prev(:)   = precl(:) - ?????

!    ! Output liquid tracers
!    call outfld(cnst_name(ixcldliq), qc(:,pver:1:-1), pcols, lchnk)
!    call outfld(cnst_name(ixrain  ), qr(:,pver:1:-1), pcols, lchnk)

  end subroutine kessler_update_run

!> \section arg_table_kessler_update_finalize  Argument Table
!! | local_name | standard_name                                                | long_name                               | units       | rank | type      | kind      | intent | optional |
!! |------------|--------------------------------------------------------------|-----------------------------------------|-------------|------|-----------|-----------|--------|----------|
!! | errmsg     | error_message                                                | CCPP error message                      | none        |    0 | character | len=512   | out    | F        |
!! | errflg     | error_flag                                                   | CCPP error flag                         | flag        |    0 | integer   |           | out    | F        |
!!
  subroutine kessler_update_finalize(errmsg, errflg)

    character(len=512),      intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

    errmsg = ' '
    errflg = 0

  end subroutine kessler_update_finalize
end module kessler_update
