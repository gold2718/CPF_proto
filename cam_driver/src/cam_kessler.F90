module cam_kessler_main

use machine, only: kind_phys

implicit none

contains

subroutine cam_kessler_main_sub()

  ! Add the CCPP specific types and functions
  use :: ccpp_api,                           &
         only: ccpp_t,                       &
               ccpp_init,                    &
               ccpp_finalize,                &
               ccpp_physics_init,            &
               ccpp_physics_run,             &
               ccpp_physics_finalize,        &
               ccpp_field_add

! NOTE -- The variables managed by the CCPP are included in the the ccpp_modules.inc file in the "use" statements
#include "ccpp_modules.inc"

  implicit none


  integer                           :: i, j

  ! Create the CCPP required cdata structure
  type(ccpp_t), allocatable, target                      :: cdata(:)

  integer                                                :: ierr
  integer ,parameter :: ncols=1
  integer ,parameter :: nlevs=8
  integer ,parameter :: ntimes=11

  real(kind_phys) :: rair, cpair, latvap, pstd, rhoh2o, ztodt, dt
  integer :: pcols, pver

  real(kind_phys),allocatable     :: rho(:,:)  ! Dry air density
  real(kind_phys),allocatable     :: pk(:,:)   ! exner func.
  real(kind_phys),allocatable     :: th(:,:)   ! Potential temp.
  real(kind_phys),allocatable     :: qv(:,:)   ! Water vapor
  real(kind_phys),allocatable     :: qc(:,:)   ! Cloud water
  real(kind_phys),allocatable     :: qr(:,:)   ! Rain water
  real(kind_phys),allocatable     :: z(:,:)    ! height
  real(kind_phys),allocatable     :: precl(:)

  integer  :: istep, ncol, pver_in, nwrite, nz
  character(len=20) :: string

  ! Allocate the host variables
  allocate(k_rateConst(3))
  allocate(my_co(nlevs))
  allocate(my_o3(nlevs))
  allocate(cdata(ncols))


  ! Initialize the host variables
  state_host%Temperature = 200.
  my_co(:) = 100_kind_phys
  my_o3(:) = 1e-6_kind_phys

  open(unit=50, file='../src/fort.50')
  read(50, fmt='(2i4)') pcols, pver
  read(50, fmt='(5f22.13)') rair, cpair, latvap, pstd, rhoh2o

  allocate(rho(pcols,pver))
  allocate(z(pcols,pver))
  allocate(pk(pcols,pver))
  allocate(th(pcols,pver))
  allocate(qv(pcols,pver))
  allocate(qc(pcols,pver))
  allocate(qr(pcols,pver))
  allocate(precl(pcols))

  do i = 1, ncols

      ! Use the suite information to setup the run
      call ccpp_init( '../suites/suite_cam_kessler_test_simple1.xml', cdata(i), ierr)
      if (ierr/=0) then
          write(*,'(a,i0,a)') 'An error occurred in ccpp_init for column ', i, '. Exiting...'
          stop
      end if

! use ccpp_fields.inc to call ccpp_field_add for all variables to be exposed to CCPP (this is auto-generated from /src/ccpp/scripts/ccpp_prebuild.py - the script parses tables in the cam_var_defs.f90)
#include "ccpp_fields.inc"

      ! initialize each column's physics
      call ccpp_physics_init(cdata(i), ierr=ierr)
      if (ierr/=0) then
          write(*,'(a,i0,a)') 'An error occurred in ccpp_physics_init for column ', i, '. Exiting...'
          stop
      end if

  end do

  write(6,*) ' '
  write(6,*) 'After initialization, my_co(1)=',my_co(1)
  write(6,*) 'After initialization, my_o3(1)=',my_o3(1)
  write(6,*) ' '

  ! loop over all time steps
  nwrite=0
  do j = 1, ntimes
      read(50,fmt='(a10,i4)') string(1:19),nwrite
      read(50,fmt='(a20,2i4,f20.13)') string(1:19),ncol, pver_in, ztodt
      read(50,fmt='(a10,(e22.15))') string,rho(:ncol,:)
      read(50,fmt='(a10,(e22.15))') string,z(:ncol,:)
      read(50,fmt='(a10,(e22.15))') string,pk(:ncol,:)
      read(50,fmt='(a10,(e22.15))') string,th(:ncol,:)
      read(50,fmt='(a10,(e22.15))') string,qv(:ncol,:)
      read(50,fmt='(a10,(e22.15))') string,qc(:ncol,:)
      read(50,fmt='(a10,(e22.15))') string,qr(:ncol,:)
      read(50,fmt='(a10,(e22.15))') string,precl(:ncol)
    write(6,*) 'At time step', j, 'in host model state_host%Temperature =', state_host%Temperature
    do i = 1, ncols
     call ccpp_physics_run(cdata(i), ierr=ierr)
!    call kessler_run(ncol, pver, ztodt, rho(:ncol,:), z(:ncol,:), pk(:ncol,:),    &
!         th(:ncol,:), qv(:ncol,:), qc(:ncol,:), qr(:ncol,:), precl(:ncol), errmsg)
      write(51,'(a10,i4)') 'nwrite=',nwrite
      write(51,'(a20,2i4,f20.13)') 'ncol, pver, ztodt=',ncol, pver, ztodt
      write(51,fmt='(a10,(e22.15))') 'rho=',rho(:ncol,:)
      write(51,fmt='(a10,(e22.15))') 'z=',z(:ncol,:)
      write(51,fmt='(a10,(e22.15))') 'pk=',pk(:ncol,:)
      write(51,fmt='(a10,(e22.15))') 'th=',th(:ncol,:)
      write(51,fmt='(a10,(e22.15))') 'qv=',qv(:ncol,:)
      write(51,fmt='(a10,(e22.15))') 'qc=',qc(:ncol,:)
      write(51,fmt='(a10,(e22.15))') 'qr=',qr(:ncol,:)
      write(51,fmt='(a10,(e22.15))') 'precl=',precl(:ncol)
      nwrite = nwrite+1

       if (ierr/=0) then
           write(*,*) errmsg
           write(*,'(a,i0,a)') 'An error occurred in ccpp_physics_run for column ', i, '. Exiting...'
           stop
       end if

       write(6,*) ' At time j=',j,' my_co(1)=',my_co(1)
       write(6,*) ' At time j=',j,' my_o3(1)=',my_o3(1)
     end do

     ! Decrease the temperature every time step
     state_host%Temperature = state_host%Temperature - 100._kind_phys

  end do


  do i=1, ncols
      call ccpp_finalize(cdata(i), ierr)
      if (ierr/=0) then
          write(*,'(a,i0,a)') 'An error occurred in ccpp_finalize for column ', i, '. Exiting...'
          stop
      end if
  end do

end subroutine cam_kessler_main_sub

end module cam_kessler_main

!> \brief Main SCM program that calls the main SCM subroutine
!!
!! The Doxygen documentation system cannot handle in-body comments in Fortran main programs, so the "main" program was put in the
!! subroutine \ref cam_kessler_main_sub above.
program cam_kessler
  use cam_kessler_main
  call cam_kessler_main_sub()
end program cam_kessler
