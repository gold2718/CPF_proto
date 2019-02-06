module kessler_update

use machine, only: r8 => kind_phys
use geopotential, only: geopotential_t

implicit none
private
public :: kessler_update_timestep_init, kessler_update_run, kessler_update_finalize

contains

!> \section arg_table_kessler_update_timestep_init  Argument Table
!! [ temp ]
!!   standard_name = temperature
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = K
!!   dimensions = (horizontal_dimension, vertical_layer_dimension)
!!    intent = in
!! [ temp_prev ]
!!   standard_name = temperature_from_previous_timestep
!!   units = K
!!   dimensions = (horizontal_dimension, vertical_layer_dimension)
!!   type = real
!!   kind = kind_phys
!!   intent = out
!! [ ttend_t ]
!!   standard_name = total_tendency_of_temperature
!!   units = K s-1
!!   dimensions = (horizontal_dimension, vertical_layer_dimension)
!!   type = real
!!   kind = kind_phys
!!   intent = out
!! [ errmsg ]
!!   standard_name = ccpp_error_message
!!   long_name = Error message for error handling in CCPP
!!   units = 1
!!   type = character | kind = len=512
!!   dimensions = ()
!!   intent = out
!! [ errflg ]
!!   standard_name = ccpp_error_flag
!!   long_name = Error flag for error handling in CCPP
!!   units = flag
!!   type = integer
!!   dimensions = ()
!!   intent = out
!!
  subroutine kessler_update_timestep_init(temp, temp_prev, ttend_t, errmsg, errflg)

    real(r8), intent(in)    :: temp(:,:)
    real(r8), intent(out)   :: temp_prev(:,:)
    real(r8), intent(inout) :: ttend_t(:,:)
    character(len=512),      intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

!   Initialize the previous temperature and its tendency to zero
    temp_prev(:,:)  = temp(:,:)
    ttend_t(:,:)    = 0._r8

    errmsg = ' '
    errflg = 0

  end subroutine kessler_update_timestep_init

!> \section arg_table_kessler_update_run  Argument Table
!! [ nz ]
!!    standard_name = vertical_layer_dimension
!!    long_name = number of vertical levels
!!    units = 1
!!    dimensions = ()
!!    type = integer
!!    intent = in
!! [ pcols ]
!!    standard_name = horizontal_dimension
!!    units = 1
!!    dimensions = ()
!!    type = integer
!!    intent = in
!! [ ncol ]
!!    standard_name = horizontal_loop_extent
!!    long_name = number of columns
!!    units = 1
!!    dimensions = ()
!!    type = integer
!!    intent = in
!! [ gravit ]
!!    standard_name = acceleration_of_gravity
!!    units = m s-2
!!    dimensions = ()
!!    type = real
!!    kind = kind_phys
!!    intent = in
!! [ cpair ]
!!    standard_name = specific_heat_of_dry_air_at_constant_pressure
!!    units = J kg-1 K-1
!!    dimensions = ()
!!    type = real
!!    kind = kind_phys
!!    intent = in
!! [ rair ]
!!    standard_name = gas_constant_dry_air
!!    units = J kg-1 K-1
!!    dimensions = ()
!!    type = real
!!    kind = kind_phys
!!    intent = in
!! [ zvir ]
!!    standard_name = ratio_of_dry_air_to_water_vapor_gas_constants_minus_one
!!    units = -1
!!    dimensions = ()
!!    type = real
!!    kind = kind_phys
!!    intent = in
!! [ phis ]
!!   standard_name = surface_geopotential
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = m2 s-2
!!   dimensions = (horizontal_dimension)
!!    intent = in
!! [ temp ]
!!   standard_name = temperature
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = K
!!   dimensions = (horizontal_dimension, vertical_layer_dimension)
!!    intent = in
!! [ lnpint ]
!!   standard_name = ln_air_pressure_at_interface
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = ln(Pa)
!!   dimensions = (horizontal_dimension, vertical_level_dimension)
!!    intent = in
!! [ lnpmid ]
!!   standard_name = natural_log_of_air_pressure
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = pmid
!!   dimensions = (horizontal_dimension, vertical_layer_dimension)
!!    intent = in
!! [ pint ]
!!   standard_name = air_pressure_at_interface
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = Pa
!!   dimensions = (horizontal_dimension, vertical_level_dimension)
!!    intent = in
!! [ pmid ]
!!   standard_name = air_pressure
!!   long_name = Midpoint air pressure
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = Pa
!!   dimensions = (horizontal_dimension, vertical_layer_dimension)
!!    intent = in
!! [ pdel ]
!!   standard_name = pressure_thickness
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = Pa
!!   dimensions = (horizontal_dimension, vertical_layer_dimension)
!!    intent = in
!! [ rpdel ]
!!   standard_name = reciprocal_of_pressure_thickness
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = Pa-1
!!   dimensions = (horizontal_dimension, vertical_layer_dimension)
!!    intent = in
!! [ qc ]
!!   standard_name = water_vapor_specific_humidity
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = kg kg-1
!!   dimensions = (horizontal_dimension, vertical_layer_dimension)
!!    intent = in
!! [ theta ]
!!   standard_name = potential_temperature
!!   long_name = potential temperature
!!   units = K
!!   dimensions = (horizontal_dimension, vertical_layer_dimension)
!!   type = real
!!   kind = kind_phys
!!   intent = inout
!! [ exner ]
!!   standard_name = inverse_exner_function_wrt_surface_pressure
!!   long_name = inverse exner function w.r.t. surface pressure, (ps/p)^(R/cp)
!!   state_variable = true
!!   units = 1
!!   type = real | kind = kind_phys
!!   dimensions = (horizontal_loop_extent, vertical_layer_dimension)
!!   intent = in
!! [ dt ]
!!   standard_name = time_step_for_physics
!!   long_name = time step
!!   units = s
!!   dimensions = ()
!!   type = real
!!   kind = kind_phys
!!   intent = in
!! [ zi ]
!!   standard_name = geopotential_height_above_surface_at_interfaces
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = m
!!   dimensions = (horizontal_dimension, vertical_level_dimension)
!!    intent = inout
!! [ zm ]
!!   standard_name = geopotential_height_above_surface_at_midpoints
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = m
!!   dimensions = (horizontal_loop_extent, vertical_layer_dimension)
!!   intent = inout
!! [ temp_prev ]
!!    standard_name = temperature_from_previous_timestep
!!    units = K
!!    dimensions = (horizontal_loop_extent, vertical_layer_dimension)
!!    type = real
!!    kind = kind_phys
!!    intent = inout
!! [ ttend_t ]
!!    standard_name = total_tendency_of_temperature
!!   type = real
!!   kind = kind_phys
!!   units = K s-1
!!   dimensions = (horizontal_loop_extent, vertical_layer_dimension)
!!    intent = out
!! [ st_energy ]
!!   standard_name = dry_static_energy_content_of_atmosphere_layer
!!   long_name = Dry static energy
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = J m-2
!!   dimensions = (horizontal_loop_extent, vertical_layer_dimension)
!!    intent = inout
!! [ errmsg ]
!!   standard_name = ccpp_error_message
!!   long_name = Error message for error handling in CCPP
!!   units = 1
!!   type = character | kind = len=512
!!   dimensions = ()
!!   intent = out
!! [ errflg ]
!!   standard_name = ccpp_error_flag
!!   long_name = Error flag for error handling in CCPP
!!   units = flag
!!   type = integer
!!   dimensions = ()
!!   intent = out
!!
  subroutine kessler_update_run(nz, pcols, ncol, gravit, cpair, rair, zvir, phis, temp, &
                 lnpint, lnpmid, pint, pmid, pdel, rpdel, qc, theta, exner, dt,  &
                 zi, zm, temp_prev, ttend_t, st_energy, errmsg, errflg )

    integer, intent(in)     :: nz
    integer, intent(in)     :: pcols
    integer, intent(in)     :: ncol
    real(r8), intent(in)    :: gravit
    real(r8), intent(in)    :: cpair
    real(r8), intent(in)    :: rair
    real(r8), intent(in)    :: zvir
    real(r8), intent(in)    :: phis(:)
    real(r8), intent(in)    :: temp(:,:)   ! temperature
    real(r8), intent(in)    :: lnpint(:,:)
    real(r8), intent(in)    :: lnpmid(:,:)
    real(r8), intent(in)    :: pint(:,:)
    real(r8), intent(in)    :: pmid(:,:)
    real(r8), intent(in)    :: pdel(:,:)
    real(r8), intent(in)    :: rpdel(:,:)
    real(r8), intent(in)    :: qc(:,:)
    real(r8), intent(in)    :: theta(:,:)
    real(r8), intent(in)    :: exner(:,:)
    real(r8), intent(in)    :: dt

    real(r8), intent(inout) :: zi(:,:)
    real(r8), intent(inout) :: zm(:,:)
    real(r8), intent(inout) :: temp_prev(:,:)
    real(r8), intent(inout) :: ttend_t(:,:)
    real(r8), intent(inout) :: st_energy(:,:)

    character(len=512),      intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

    integer                 :: k, rk
    integer                 :: vert_surf, vert_toa
    real(r8)                :: rairv(pcols,nz)
    real(r8)                :: zvirv(pcols,nz)
    real(r8)                :: ptend_s(pcols,nz)

    errmsg = ' '
    errflg = 0

    rairv(:,:) = rair
    zvirv(:,:) = zvir

    ! Back out tendencies from updated fields
    do k = 1, nz
      ptend_s(:ncol,k) = (theta(:ncol,k) * exner(:ncol,k) - temp_prev(:ncol,k)) * cpair / dt
      ttend_t(:ncol,k) = ttend_t(:ncol,k) + ptend_s(:ncol,k)/cpair
    end do

    ! Kessler is bottom to top
    vert_toa = 30
    vert_surf = 1

    call geopotential_t  ( nz, nz+1, .true., vert_surf, vert_toa, &
            lnpint,    lnpmid,    pint  , pmid  , pdel  , rpdel  , &
            temp     , qc, rairv, gravit  , zvirv              , &
            zi    , zm      , ncol         )

    do k = 1, nz
      st_energy(:ncol,k) = temp(:ncol,k)*cpair+gravit*zm(:ncol,k)+phis(:ncol)
    end do


!    surf_state%precl(:ncol) = surf_state%precl(:ncol) + precl(:ncol)  ! KEEPING THIS HERE AS A REMINDER

    ! Set the previous q values to the current q
    temp_prev(:,:)     = temp(:,:)

  end subroutine kessler_update_run

!> \section arg_table_kessler_update_finalize  Argument Table
!! [ errmsg ]
!!   standard_name = ccpp_error_message
!!   long_name = Error message for error handling in CCPP
!!   units = 1
!!   type = character | kind = len=512
!!   dimensions = ()
!!   intent = out
!! [ errflg ]
!!   standard_name = ccpp_error_flag
!!   long_name = Error flag for error handling in CCPP
!!   units = flag
!!   type = integer
!!   dimensions = ()
!!   intent = out
!!
  subroutine kessler_update_finalize(errmsg, errflg)

    character(len=512),      intent(out)   :: errmsg
    integer,                 intent(out)   :: errflg

    errmsg = ' '
    errflg = 0

  end subroutine kessler_update_finalize
end module kessler_update
