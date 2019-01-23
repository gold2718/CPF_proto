!-------------------------------------------------------------------------------
!physics data types module
!-------------------------------------------------------------------------------
module physics_types

   use physconst,    only: r8
   use constituents, only: pcnst

   implicit none
   private          ! Make default type private to the module

   ! Public types:

   public physics_state
   public physics_tend

   ! Public interfaces

   public physics_state_copy  ! copy a physics_state object
   public physics_tend_init   ! initialize a physics_tend object

   public physics_type_alloc

   public physics_state_alloc   ! allocate individual components within state
   public physics_state_dealloc ! deallocate individual components within state
   public physics_tend_alloc    ! allocate individual components within tend
   public physics_tend_dealloc  ! deallocate individual components within tend

   public state ! The model's physics_state variable
   public tend  ! The model's physics_tend variable

!==============================================================================
!! \section arg_table_physics_state
!! [ lat ]
!!   standard_name = latitude
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = radians
!!   dimensions = (horizontal_loop_extent)
!! [ lon ]
!!   standard_name = longitude
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = radians
!!   dimensions = (horizontal_loop_extent)
!! [ ps ]
!!   standard_name = surface_air_pressure
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = Pa
!!   dimensions = (horizontal_loop_extent)
!! [ psdry ]
!!   standard_name = dry_surface_pressure
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = Pa
!!   dimensions = (horizontal_loop_extent)
!! [ phis ]
!!   standard_name = surface_geopotential
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = m2 s-2
!!   dimensions = (horizontal_loop_extent)
!! [ ulat ]
!!   standard_name = unique_latitudes
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = radians
!!   dimensions = (horizontal_loop_extent)
!! [ ulon ]
!!   standard_name = unique_longitudes
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = radians
!!   dimensions = (horizontal_loop_extent)
!! [ t ]
!!   standard_name = temperature
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = K
!!   dimensions = (horizontal_loop_extent, vertical_layer_dimension)
!! [ u ]
!!   standard_name = eastward_wind
!!   long_name = Zonal wind
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = m s-1
!!   dimensions = (horizontal_loop_extent, vertical_layer_dimension)
!! [ v ]
!!   standard_name = northward_wind
!!   long_name = Meridional wind
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = m s-1
!!   dimensions = (horizontal_loop_extent, vertical_layer_dimension)
!! [ s ]
!!   standard_name = dry_static_energy_content_of_atmosphere_layer
!!   long_name = Dry static energy
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = J m-2
!!   dimensions = (horizontal_loop_extent, vertical_layer_dimension)
!! [ omega ]
!!   standard_name = lagrangian_tendency_of_air_pressure
!!   long_name = Vertical pressure velocity
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = Pa s-1
!!   dimensions = (horizontal_loop_extent, vertical_layer_dimension)
!! [ pmid ]
!!   standard_name = air_pressure
!!   long_name = Midpoint air pressure
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = Pa
!!   dimensions = (horizontal_loop_extent, vertical_layer_dimension)
!! [ pmiddry ]
!!   standard_name = air_pressure_of_dry_air
!!   long_name = Dry midpoint pressure
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = Pa
!!   dimensions = (horizontal_loop_extent, vertical_layer_dimension)
!! [ pdel ]
!!   standard_name = pressure_thickness
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = Pa
!!   dimensions = (horizontal_loop_extent, vertical_layer_dimension)
!! [ pdeldry ]
!!   standard_name = pressure_thickness_of_dry_air
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = Pa
!!   dimensions = (horizontal_loop_extent, vertical_layer_dimension)
!! [ rpdel ]
!!   standard_name = reciprocal_of_pressure_thickness
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = Pa-1
!!   dimensions = (horizontal_loop_extent, vertical_layer_dimension)
!! [ rpdeldry ]
!!   standard_name = reciprocal_of_pressure_thickness_of_dry_air
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = Pa-1
!!   dimensions = (horizontal_loop_extent, vertical_layer_dimension)
!! [ lnpmid ]
!!   standard_name = natural_log_of_air_pressure
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = pmid
!!   dimensions = (horizontal_loop_extent, vertical_layer_dimension)
!! [ lnpmiddry ]
!!   standard_name = log_of_air_pressure_of_dry_air
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = Pa
!!   dimensions = (horizontal_loop_extent, vertical_layer_dimension)
!! [ exner ]
!!   standard_name = inverse_exner_function_wrt_surface_pressure
!!   long_name = inverse exner function w.r.t. surface pressure, (ps/p)^(R/cp)
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = 1
!!   dimensions = (horizontal_loop_extent, vertical_layer_dimension)
!! [ zm ]
!!   standard_name = geopotential_height_above_surface_at_midpoints
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = m
!!   dimensions = (horizontal_loop_extent, vertical_layer_dimension)
!! [ q ]
!!   standard_name = constituent_mixing_ratio
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = kg/kg moist or dry air depending on type
!!   dimensions = (horizontal_loop_extent, vertical_layer_dimension, number_of_tracers)
!! [ q(:,:,index_of_water_vapor_specific_humidity) ]
!!   standard_name = water_vapor_specific_humidity
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = kg kg-1
!!   dimensions = (horizontal_loop_extent, vertical_layer_dimension)
!! [ q(:,:,index_of_cloud_liquid_water_mixing_ratio) ]
!!   standard_name = cloud_liquid_water_mixing_ratio
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = kg kg-1
!!   dimensions = (horizontal_loop_extent, vertical_layer_dimension)
!! [ q(:,:,index_of_rain_water_mixing_ratio) ]
!!   standard_name = rain_water_mixing_ratio
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = kg kg-1
!!   dimensions = (horizontal_loop_extent, vertical_layer_dimension)
!! [ pint ]
!!   standard_name = air_pressure_at_interface
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = Pa
!!   dimensions = (horizontal_loop_extent, vertical_level_dimension)
!! [ pintdry ]
!!   standard_name = interface_pressure_dry
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = Pa
!!   dimensions = (horizontal_loop_extent, vertical_level_dimension)
!! [ lnpint ]
!!   standard_name = ln_air_pressure_at_interface
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = ln(Pa)
!!   dimensions = (horizontal_loop_extent, vertical_level_dimension)
!! [ lnpintdry ]
!!   standard_name = ln_interface_pressure_dry
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = ln(Pa)
!!   dimensions = (horizontal_loop_extent, vertical_level_dimension)
!! [ zi ]
!!   standard_name = geopotential_height_above_surface_at_interfaces
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = m
!!   dimensions = (horizontal_loop_extent, vertical_level_dimension)
!! [ te_ini ]
!!   standard_name = vertically_integrated_total_kinetic_and_static_energy_of_initial_state
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = J m-2
!!   dimensions = (horizontal_loop_extent)
!! [ te_cur ]
!!   standard_name = vertically_integrated_total_kinetic_and_static_energy_of_current_state
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = J m-2
!!   dimensions = (horizontal_loop_extent)
!! [ tw_ini ]
!!   standard_name = vertically_integrated_total_water_of_initial_state
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = Pa s2 m-1
!!   dimensions = (horizontal_loop_extent)
!! [ tw_cur ]
!!   standard_name = vertically_integrated_total_water_of_new_state
!!   state_variable = true
!!   type = real
!!   kind = kind_phys
!!   units = Pa s2 m-1
!!   dimensions = (horizontal_loop_extent)
!! [ count ]
!!   standard_name = count_of_values_with_significant_energy_or_water_imbalances
!!   state_variable = true
!!   type = integer
!!   kind = kind_phys
!!   units = 1
!!   dimensions = (horizontal_loop_extent)
!! [ cid ]
!!   standard_name = unique_column_id
!!   state_variable = true
!!   type = integer
!!   kind = kind_phys
!!   units = 1
!!   dimensions = (horizontal_loop_extent)
!!
   type physics_state
      integer                                     :: &
           lchnk,                &! chunk index
           ngrdcol,              &! -- Grid        -- number of active columns (on the grid)
           psetcols=0,           &! --             -- max number of columns set - if subcols = pcols*psubcols, else = pcols
           ncol=0                 ! --             -- sum of nsubcol for all ngrdcols - number of active columns
      real(r8), dimension(:), allocatable         :: &
           lat,     &! latitude (radians)
           lon,     &! longitude (radians)
           ps,      &! surface pressure
           psdry,   &! dry surface pressure
           phis,    &! surface geopotential
           ulat,    &! unique latitudes  (radians)
           ulon      ! unique longitudes (radians)
      real(r8), dimension(:,:),allocatable        :: &
           t,       &! temperature (K)
           u,       &! zonal wind (m/s)
           v,       &! meridional wind (m/s)
           s,       &! dry static energy
           omega,   &! vertical pressure velocity (Pa/s)
           pmid,    &! midpoint pressure (Pa)
           pmiddry, &! midpoint pressure dry (Pa)
           pdel,    &! layer thickness (Pa)
           pdeldry, &! layer thickness dry (Pa)
           rpdel,   &! reciprocal of layer thickness (Pa)
           rpdeldry,&! recipricol layer thickness dry (Pa)
           lnpmid,  &! ln(pmid)
           lnpmiddry,&! log midpoint pressure dry (Pa)
           exner,   &! inverse exner function w.r.t. surface pressure (ps/p)^(R/cp)
           zm        ! geopotential height above surface at midpoints (m)

      real(r8), dimension(:,:,:),allocatable      :: &
           q         ! constituent mixing ratio (kg/kg moist or dry air depending on type)

      real(r8), dimension(:,:),allocatable        :: &
           pint,    &! interface pressure (Pa)
           pintdry, &! interface pressure dry (Pa)
           lnpint,  &! ln(pint)
           lnpintdry,&! log interface pressure dry (Pa)
           zi        ! geopotential height above surface at interfaces (m)

      real(r8), dimension(:),allocatable          :: &
           te_ini,  &! vertically integrated total (kinetic + static) energy of initial state
           te_cur,  &! vertically integrated total (kinetic + static) energy of current state
           tw_ini,  &! vertically integrated total water of initial state
           tw_cur    ! vertically integrated total water of new state
      integer :: count ! count of values with significant energy or water imbalances
      integer, dimension(:),allocatable           :: &
           latmapback, &! map from column to unique lat for that column
           lonmapback, &! map from column to unique lon for that column
           cid        ! unique column id
      integer :: ulatcnt, &! number of unique lats in chunk
           uloncnt   ! number of unique lons in chunk

   end type physics_state

!==============================================================================
!!  \section arg_table_physics_tend
!! [ dtdt ]
!!   standard_name = total_tendency_of_temperature
!!   state_variable = false
!!   type = real
!!   kind = kind_phys
!!   units = K s-1
!!   dimensions = (horizontal_loop_extent, vertical_layer_dimension)
!! [ dudt ]
!!   standard_name = total_tendency_of_eastward_wind
!!   state_variable = false
!!   type = real
!!   kind = kind_phys
!!   units = m s-2
!!   dimensions = (horizontal_loop_extent, vertical_layer_dimension)
!! [ dvdt ]
!!   standard_name = total_tendency_of_northward_wind
!!   state_variable = false
!!   type = real
!!   kind = kind_phys
!!   units = m s-2
!!   dimensions = (horizontal_loop_extent, vertical_layer_dimension)
!! [ flx_net ]
!!   standard_name = surface_energy_flux
!!   state_variable = false
!!   type = real
!!   kind = kind_phys
!!   units = W m-2
!!   dimensions = (horizontal_loop_extent)
!! [ te_tnd ]
!!   standard_name = cumulative_boundary_flux_of_total_energy
!!   state_variable = false
!!   type = real
!!   kind = kind_phys
!!   units = W m-2
!!   dimensions = (horizontal_loop_extent)
!! [ tw_tnd ]
!!   standard_name = cumulative_boundary_flux_of_total_water
!!   state_variable = false
!!   type = real
!!   kind = kind_phys
!!   units = W m-2
!!   dimensions = (horizontal_loop_extent)
!!
   type physics_tend

      integer   ::   psetcols=0 ! max number of columns set- if subcols = pcols*psubcols, else = pcols

      real(r8), dimension(:,:),allocatable        :: dtdt, dudt, dvdt
      real(r8), dimension(:),  allocatable        :: flx_net
      real(r8), dimension(:),  allocatable        :: &
           te_tnd,  &! cumulative boundary flux of total energy
           tw_tnd    ! cumulative boundary flux of total water
   end type physics_tend

   type(physics_state) :: state
   type(physics_tend)  :: tend

!===============================================================================
contains
!===============================================================================

   subroutine physics_type_alloc(phys_state, phys_tend, begchunk, endchunk, psetcols)

      use ppgrid,           only: pcols
!      use cam_logfile,      only: iulog

      implicit none

      type(physics_state), pointer :: phys_state(:)
      type(physics_tend), pointer :: phys_tend(:)
      integer, intent(in) :: begchunk, endchunk
      integer, intent(in) :: psetcols

      integer :: ierr=0, lchnk
      integer :: iulog = 6
      type(physics_state), pointer :: state
      type(physics_tend), pointer :: tend

      allocate(phys_state(begchunk:endchunk), stat=ierr)
      if( ierr /= 0 ) then
         write(iulog,*) 'physics_types: phys_state allocation error = ',ierr
         call endrun('physics_types: failed to allocate physics_state array')
      end if

      do lchnk=begchunk,endchunk
         call physics_state_alloc(phys_state(lchnk),lchnk,pcols)
      end do

      allocate(phys_tend(begchunk:endchunk), stat=ierr)
      if( ierr /= 0 ) then
         write(iulog,*) 'physics_types: phys_tend allocation error = ',ierr
         call endrun('physics_types: failed to allocate physics_tend array')
      end if

      do lchnk=begchunk,endchunk
         call physics_tend_alloc(phys_tend(lchnk),phys_state(lchnk)%psetcols)
      end do

   end subroutine physics_type_alloc
!===============================================================================

   subroutine physics_state_copy(state_in, state_out)

      use constituents,     only: pcnst
      use ppgrid,           only: pver, pverp

      implicit none

      !
      ! Arguments
      !
      type(physics_state), intent(in)    :: state_in
      type(physics_state), intent(out)   :: state_out

      !
      ! Local variables
      !
      integer i, k, m, ncol

      ! Allocate state_out with same subcol dimension as state_in
      call physics_state_alloc ( state_out, state_in%lchnk, state_in%psetcols)

      ncol = state_in%ncol

      state_out%psetcols = state_in%psetcols
      state_out%ngrdcol  = state_in%ngrdcol
      state_out%lchnk    = state_in%lchnk
      state_out%ncol     = state_in%ncol
      state_out%count    = state_in%count

      do i = 1, ncol
         state_out%lat(i)    = state_in%lat(i)
         state_out%lon(i)    = state_in%lon(i)
         state_out%ps(i)     = state_in%ps(i)
         state_out%phis(i)   = state_in%phis(i)
         state_out%te_ini(i) = state_in%te_ini(i)
         state_out%te_cur(i) = state_in%te_cur(i)
         state_out%tw_ini(i) = state_in%tw_ini(i)
         state_out%tw_cur(i) = state_in%tw_cur(i)
      end do

      do k = 1, pver
         do i = 1, ncol
            state_out%t(i,k)         = state_in%t(i,k)
            state_out%u(i,k)         = state_in%u(i,k)
            state_out%v(i,k)         = state_in%v(i,k)
            state_out%s(i,k)         = state_in%s(i,k)
            state_out%omega(i,k)     = state_in%omega(i,k)
            state_out%pmid(i,k)      = state_in%pmid(i,k)
            state_out%pdel(i,k)      = state_in%pdel(i,k)
            state_out%rpdel(i,k)     = state_in%rpdel(i,k)
            state_out%lnpmid(i,k)    = state_in%lnpmid(i,k)
            state_out%exner(i,k)     = state_in%exner(i,k)
            state_out%zm(i,k)        = state_in%zm(i,k)
         end do
      end do

      do k = 1, pverp
         do i = 1, ncol
            state_out%pint(i,k)      = state_in%pint(i,k)
            state_out%lnpint(i,k)    = state_in%lnpint(i,k)
            state_out%zi(i,k)        = state_in% zi(i,k)
         end do
      end do


      do i = 1, ncol
         state_out%psdry(i)  = state_in%psdry(i)
      end do
      do k = 1, pver
         do i = 1, ncol
            state_out%lnpmiddry(i,k) = state_in%lnpmiddry(i,k)
            state_out%pmiddry(i,k)   = state_in%pmiddry(i,k)
            state_out%pdeldry(i,k)   = state_in%pdeldry(i,k)
            state_out%rpdeldry(i,k)  = state_in%rpdeldry(i,k)
         end do
      end do
      do k = 1, pverp
         do i = 1, ncol
            state_out%pintdry(i,k)   = state_in%pintdry(i,k)
            state_out%lnpintdry(i,k) = state_in%lnpintdry(i,k)
         end do
      end do

      do m = 1, pcnst
         do k = 1, pver
            do i = 1, ncol
               state_out%q(i,k,m) = state_in%q(i,k,m)
            end do
         end do
      end do

   end subroutine physics_state_copy
!===============================================================================

   subroutine physics_tend_init(tend)

!      use shr_kind_mod,     only: r8 => shr_kind_r8

      implicit none

      !
      ! Arguments
      !
      type(physics_tend), intent(inout) :: tend

      !
      ! Local variables
      !

      if (.not. allocated(tend%dtdt)) then
         call endrun('physics_tend_init: tend must be allocated before it can be initialized')
      end if

      tend%dtdt    = 0._r8
      tend%dudt    = 0._r8
      tend%dvdt    = 0._r8
      tend%flx_net = 0._r8
      tend%te_tnd  = 0._r8
      tend%tw_tnd  = 0._r8

   end subroutine physics_tend_init

!===============================================================================

   subroutine physics_state_alloc(state,lchnk,psetcols)

!      use shr_infnan_mod,   only: inf, assignment(=)
      use ppgrid,           only: pver
      use constituents,     only: pcnst

      ! allocate the individual state components

      type(physics_state), intent(inout) :: state
      integer,intent(in)                 :: lchnk

      integer, intent(in)                :: psetcols

      integer :: ierr=0, i

      state%lchnk    = lchnk
      state%psetcols = psetcols
      state%ngrdcol  = psetcols  ! Number of grid columns

      !----------------------------------
      ! Following variables will be overwritten by sub-column generator, if sub-columns are being used

      !  state%ncol - is initialized in physics_state_set_grid,  if not using sub-columns

      !----------------------------------

      allocate(state%lat(psetcols), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%lat')

      allocate(state%lon(psetcols), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%lon')

      allocate(state%ps(psetcols), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%ps')

      allocate(state%psdry(psetcols), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%psdry')

      allocate(state%phis(psetcols), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%phis')

      allocate(state%ulat(psetcols), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%ulat')

      allocate(state%ulon(psetcols), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%ulon')

      allocate(state%t(psetcols,pver), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%t')

      allocate(state%u(psetcols,pver), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%u')

      allocate(state%v(psetcols,pver), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%v')

      allocate(state%s(psetcols,pver), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%s')

      allocate(state%omega(psetcols,pver), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%omega')

      allocate(state%pmid(psetcols,pver), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%pmid')

      allocate(state%pmiddry(psetcols,pver), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%pmiddry')

      allocate(state%pdel(psetcols,pver), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%pdel')

      allocate(state%pdeldry(psetcols,pver), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%pdeldry')

      allocate(state%rpdel(psetcols,pver), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%rpdel')

      allocate(state%rpdeldry(psetcols,pver), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%rpdeldry')

      allocate(state%lnpmid(psetcols,pver), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%lnpmid')

      allocate(state%lnpmiddry(psetcols,pver), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%lnpmiddry')

      allocate(state%exner(psetcols,pver), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%exner')

      allocate(state%zm(psetcols,pver), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%zm')

      allocate(state%q(psetcols,pver,pcnst), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%q')

      allocate(state%pint(psetcols,pver+1), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%pint')

      allocate(state%pintdry(psetcols,pver+1), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%pintdry')

      allocate(state%lnpint(psetcols,pver+1), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%lnpint')

      allocate(state%lnpintdry(psetcols,pver+1), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%lnpintdry')

      allocate(state%zi(psetcols,pver+1), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%zi')

      allocate(state%te_ini(psetcols), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%te_ini')

      allocate(state%te_cur(psetcols), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%te_cur')

      allocate(state%tw_ini(psetcols), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%tw_ini')

      allocate(state%tw_cur(psetcols), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%tw_cur')

      allocate(state%latmapback(psetcols), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%latmapback')

      allocate(state%lonmapback(psetcols), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%lonmapback')

      allocate(state%cid(psetcols), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_alloc error: allocation error for state%cid')

      state%lat(:) = -999
      state%lon(:) = -999
      state%ulat(:) = -999
      state%ulon(:) = -999
      state%ps(:) = -999
      state%psdry(:) = -999
      state%phis(:) = -999
      state%t(:,:) = -999
      state%u(:,:) = -999
      state%v(:,:) = -999
      state%s(:,:) = -999
      state%omega(:,:) = -999
      state%pmid(:,:) = -999
      state%pmiddry(:,:) = -999
      state%pdel(:,:) = -999
      state%pdeldry(:,:) = -999
      state%rpdel(:,:) = -999
      state%rpdeldry(:,:) = -999
      state%lnpmid(:,:) = -999
      state%lnpmiddry(:,:) = -999
      state%exner(:,:) = -999
      state%zm(:,:) = -999
      state%q(:,:,:) = -999

      state%pint(:,:) = -999
      state%pintdry(:,:) = -999
      state%lnpint(:,:) = -999
      state%lnpintdry(:,:) = -999
      state%zi(:,:) = -999

      state%te_ini(:) = -999
      state%te_cur(:) = -999
      state%tw_ini(:) = -999
      state%tw_cur(:) = -999

   end subroutine physics_state_alloc

!===============================================================================

   subroutine physics_state_dealloc(state)

      ! deallocate the individual state components

      type(physics_state), intent(inout) :: state
      integer                            :: ierr = 0

      deallocate(state%lat, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%lat')

      deallocate(state%lon, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%lon')

      deallocate(state%ps, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%ps')

      deallocate(state%psdry, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%psdry')

      deallocate(state%phis, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%phis')

      deallocate(state%ulat, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%ulat')

      deallocate(state%ulon, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%ulon')

      deallocate(state%t, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%t')

      deallocate(state%u, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%u')

      deallocate(state%v, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%v')

      deallocate(state%s, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%s')

      deallocate(state%omega, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%omega')

      deallocate(state%pmid, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%pmid')

      deallocate(state%pmiddry, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%pmiddry')

      deallocate(state%pdel, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%pdel')

      deallocate(state%pdeldry, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%pdeldry')

      deallocate(state%rpdel, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%rpdel')

      deallocate(state%rpdeldry, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%rpdeldry')

      deallocate(state%lnpmid, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%lnpmid')

      deallocate(state%lnpmiddry, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%lnpmiddry')

      deallocate(state%exner, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%exner')

      deallocate(state%zm, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%zm')

      deallocate(state%q, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%q')

      deallocate(state%pint, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%pint')

      deallocate(state%pintdry, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%pintdry')

      deallocate(state%lnpint, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%lnpint')

      deallocate(state%lnpintdry, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%lnpintdry')

      deallocate(state%zi, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%zi')

      deallocate(state%te_ini, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%te_ini')

      deallocate(state%te_cur, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%te_cur')

      deallocate(state%tw_ini, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%tw_ini')

      deallocate(state%tw_cur, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%tw_cur')

      deallocate(state%latmapback, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%latmapback')

      deallocate(state%lonmapback, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%lonmapback')

      deallocate(state%cid, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_state_dealloc error: deallocation error for state%cid')

   end subroutine physics_state_dealloc

!===============================================================================

   subroutine physics_tend_alloc(tend,psetcols)

!      use shr_infnan_mod,   only: inf, assignment(=)
      use ppgrid,           only: pver
      ! allocate the individual tend components

      type(physics_tend), intent(inout)  :: tend

      integer, intent(in)                :: psetcols

      integer :: ierr = 0

      tend%psetcols = psetcols

      allocate(tend%dtdt(psetcols,pver), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_tend_alloc error: allocation error for tend%dtdt')

      allocate(tend%dudt(psetcols,pver), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_tend_alloc error: allocation error for tend%dudt')

      allocate(tend%dvdt(psetcols,pver), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_tend_alloc error: allocation error for tend%dvdt')

      allocate(tend%flx_net(psetcols), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_tend_alloc error: allocation error for tend%flx_net')

      allocate(tend%te_tnd(psetcols), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_tend_alloc error: allocation error for tend%te_tnd')

      allocate(tend%tw_tnd(psetcols), stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_tend_alloc error: allocation error for tend%tw_tnd')

      tend%dtdt(:,:) = -999
      tend%dudt(:,:) = -999
      tend%dvdt(:,:) = -999
      tend%flx_net(:) = -999
      tend%te_tnd(:) = -999
      tend%tw_tnd(:) = -999

   end subroutine physics_tend_alloc

!===============================================================================

   subroutine physics_tend_dealloc(tend)

      ! deallocate the individual tend components

      type(physics_tend), intent(inout)  :: tend
      integer :: psetcols
      integer :: ierr = 0

      deallocate(tend%dtdt, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_tend_dealloc error: deallocation error for tend%dtdt')

      deallocate(tend%dudt, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_tend_dealloc error: deallocation error for tend%dudt')

      deallocate(tend%dvdt, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_tend_dealloc error: deallocation error for tend%dvdt')

      deallocate(tend%flx_net, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_tend_dealloc error: deallocation error for tend%flx_net')

      deallocate(tend%te_tnd, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_tend_dealloc error: deallocation error for tend%te_tnd')

      deallocate(tend%tw_tnd, stat=ierr)
      if ( ierr /= 0 ) call endrun('physics_tend_dealloc error: deallocation error for tend%tw_tnd')
   end subroutine physics_tend_dealloc

!===============================================================================

end module physics_types
