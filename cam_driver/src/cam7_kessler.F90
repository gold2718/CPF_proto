module cam_kessler_main

  use machine, only: kind_phys

  implicit none


contains

  subroutine cam_kessler_main_sub()

    use ppgrid,           only: pcols, pver, pverp, pcnst
    use physics_types,    only: state, tend, physics_state, physics_type_alloc
    use CAM_ccpp_cap,     only: CAM_ccpp_physics_initialize
    use CAM_ccpp_cap,     only: CAM_ccpp_physics_timestep_initial
    use CAM_ccpp_cap,     only: CAM_ccpp_physics_run
    use CAM_ccpp_cap,     only: CAM_ccpp_physics_timestep_final
    use CAM_ccpp_cap,     only: CAM_ccpp_physics_finalize
    use ccpp_physics_api, only: ccpp_physics_suite_list
    use ccpp_physics_api, only: ccpp_physics_suite_part_list

    integer,            parameter   :: begchunk = 33 ! Not needed for CAM7
    integer,            parameter   :: endchunk = 33 ! Not needed for CAM7
    integer,            parameter   :: ncols = 1
    integer,            parameter   :: ntimes = 3

    integer                         :: i, j, k, rk
    integer                         :: ierr
    integer                         :: col_start, col_end
    integer                         :: ncol, nwrite, pver_in, nwrite_in
    real(kind_phys)                 :: ztodt
    real(kind_phys)                 :: precl(pcols)
    real(kind_phys)                 :: scratch(pcols,pver)
    character(len=20)               :: string
    character(len=512)              :: errmsg
    character(len=128), allocatable :: part_names(:)
    integer                         :: errflg
    real(kind_phys)                 :: pmiddry_top2bot(pcols,pver)
    real(kind_phys)                 :: pmid_top2bot(pcols,pver)
    real(kind_phys)                 :: pdel_top2bot(pcols,pver)
    real(kind_phys)                 :: pdeldry_top2bot(pcols,pver)
    real(kind_phys)                 :: zm_top2bot(pcols,pver)
    real(kind_phys)                 :: exner_top2bot(pcols,pver)
    real(kind_phys)                 :: t_top2bot(pcols,pver)
    real(kind_phys)                 :: s_top2bot(pcols,pver)
    real(kind_phys)                 :: q_top2bot(pcols,pver,pcnst)
    real(kind_phys)                 :: wet_to_dry(pcols)
    real(kind_phys)                 :: dry_to_wet(pcols)

    ! Allocate the host variables
    call physics_type_alloc(state, tend, pcols)

    !  open(unit=50, file='../src/fort.50')
    !  read(50, fmt='(2i4)') icols, iver

    ! Use the suite information to setup the run
    call CAM_ccpp_physics_initialize('cam_kessler_test', precl, errmsg, errflg)
    if (errflg /= 0) then
       write(6, *) trim(errmsg)
       stop
    end if

    ! loop over all time steps
    nwrite=0
    do j = 1, ntimes
       ncol = pcols
       read(60,fmt='(a10,i4)') string(1:8),nwrite_in
       read(60,fmt='(a20,2i4,f20.13)') string(1:19),ncol, pver_in, ztodt
       read(60,fmt='(a20,(e22.15))') string(1:20),pmiddry_top2bot(:ncol,:pver_in)
       read(60,fmt='(a20,(e22.15))') string(1:20),pmid_top2bot(:ncol,:pver_in)
       read(60,fmt='(a20,(e22.15))') string,pdel_top2bot(:ncol,:pver_in)
       read(60,fmt='(a20,(e22.15))') string,pdeldry_top2bot(:ncol,:pver_in)
       read(60,fmt='(a20,(e22.15))') string,exner_top2bot(:ncol,:pver_in)
       read(60,fmt='(a20,(e22.15))') string,t_top2bot(:ncol,:pver_in)
       read(60,fmt='(a20,(e22.15))') string,s_top2bot(:ncol,:pver_in)
       read(60,fmt='(a20,(e22.15))') string,q_top2bot(:ncol,:pver_in,1)
       read(60,fmt='(a20,(e22.15))') string,q_top2bot(:ncol,:pver_in,2)
       read(60,fmt='(a20,(e22.15))') string,q_top2bot(:ncol,:pver_in,3)
       read(60,fmt='(a20,(e22.15))') string,precl(:ncol)
       call CAM_ccpp_physics_timestep_initial('cam_kessler_test', precl, errmsg, errflg)
       col_start = 1
       col_end = ncol

       ! Need to swap the bottom and top
       ! Swap the pdel and pdeldry so they can be used in the wet_to_dry calc
       do k=1,pver
         rk= pver - k +1 
         state%pdel(:ncol,rk)    = pdel_top2bot(:ncol,k)
         state%pdeldry(:ncol,rk) = pdeldry_top2bot(:ncol,k)
         state%pmiddry(:ncol,rk) = pmiddry_top2bot(:ncol,k)
         state%pmid(:ncol,rk) = pmid_top2bot(:ncol,k)
         state%exner(:ncol,rk)   = exner_top2bot(:ncol,k)
         state%t(:ncol,rk)       = t_top2bot(:ncol,k)
         state%s(:ncol,rk)       = s_top2bot(:ncol,k)
       end do

      ! Also convert the q from wet to dry
       do k=1,pver
         rk= pver - k +1 
         wet_to_dry(:ncol) = state%pdel(:ncol,rk)/state%pdeldry(:ncol,rk) ! use state vert index
         state%q(:ncol,rk,1)     = q_top2bot(:ncol,k,1)*wet_to_dry(:ncol)
         state%q(:ncol,rk,2)     = q_top2bot(:ncol,k,2)*wet_to_dry(:ncol)
         state%q(:ncol,rk,3)     = q_top2bot(:ncol,k,3)*wet_to_dry(:ncol)
      end do

       call CAM_ccpp_physics_run('cam_kessler_test', 'physics', col_start, col_end, precl, errmsg, errflg)
       if (errflg /= 0) then
          write(6, *) trim(errmsg)
          call ccpp_physics_suite_part_list('cam_kessler_test', part_names, errmsg, errflg)
          write(6, *) 'Available suite parts are:'
          do nwrite = 1, size(part_names)
             write(6, *) trim(part_names(nwrite))
          end do
          stop
       end if

       write(6,*) 'At time step', j, 'in host model Temperature =', state%T(8, :pver)

       ! Convert the q from dry to wet
       do k=1,pver
         dry_to_wet(:ncol) = state%pdeldry(:ncol,k)/state%pdel(:ncol,k)
         state%q(:ncol,k,1) = state%q(:ncol,k,1)*dry_to_wet(:ncol)
         state%q(:ncol,k,2) = state%q(:ncol,k,2)*dry_to_wet(:ncol)
         state%q(:ncol,k,3) = state%q(:ncol,k,3)*dry_to_wet(:ncol)
       end do

         write(61,'(a10,i4)') 'nwrite=',nwrite
         write(61,'(a20,2i4,f20.13)') 'ncol, pver, ztodt=',ncol, pver, ztodt
         write(61,fmt='(a10,(e22.15))') 'pmiddry=',state%pmiddry(:ncol,pver:1:-1)
         write(61,fmt='(a10,(e22.15))') 'pmid=',state%pmid(:ncol,pver:1:-1)
         write(61,fmt='(a10,(e22.15))') 'pdel=',state%pdel(:ncol,pver:1:-1)
         write(61,fmt='(a10,(e22.15))') 'pdeldry=',state%pdeldry(:ncol,pver:1:-1)
         write(61,fmt='(a10,(e22.15))') 'exner=',state%exner(:ncol,pver:1:-1)
         write(61,fmt='(a10,(e22.15))') 'state%t=',state%t(:ncol,pver:1:-1)
         write(61,fmt='(a10,(e22.15))') 'state%s=',state%s(:ncol,pver:1:-1)
         write(61,fmt='(a10,(e22.15))') 'qv=',state%q(:ncol,pver:1:-1,1)
         write(61,fmt='(a10,(e22.15))') 'qc=',state%q(:ncol,pver:1:-1,2)
         write(61,fmt='(a10,(e22.15))') 'qr=',state%q(:ncol,pver:1:-1,3)
         write(61,fmt='(a10,(e22.15))') 'precl=',precl(:ncol)

       nwrite = nwrite + 1

       call CAM_ccpp_physics_timestep_final('cam_kessler_test', precl, errmsg, errflg)

    end do


    call CAM_ccpp_physics_finalize('cam_kessler_test', precl, errmsg, errflg)
    if (errflg /= 0) then
       write(6, *) trim(errmsg)
       write(6,'(a)') 'An error occurred in ccpp_timestep_final, Exiting...'
       stop
    end if

  end subroutine cam_kessler_main_sub

end module cam_kessler_main

!> \brief Main SCM program that calls the main SCM subroutine
!!
!! The Doxygen documentation system cannot handle in-body comments in Fortran main programs, so the "main" program was put in the
!! subroutine \ref cam_kessler_main_sub above.
program cam_kessler
  use cam_kessler_main, only: cam_kessler_main_sub
  call cam_kessler_main_sub()
end program cam_kessler
