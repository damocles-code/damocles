!----------------------------------------------------------------------------!
!  this subroutine is only called if the python mcmc wrapper is being used.  !
!  the logical lg_mcmc is set to true and input parameters passed from the   !
!  python routine.                                                           !
!----------------------------------------------------------------------------!
subroutine run_damocles_wrap(params,flags,n,ml_lg,mcmc_mod)  

    use globals
    use input
    use class_dust
    use initialise
    use vector_functions
    use driver
    use bayesian_input

  implicit none

  real :: params(21,6)
  real :: flags(21,6)
  integer :: n
  integer :: ml_lg
  real :: mcmc_mod(n,2)
  character(30) :: input_folder(6)
  
  !f2py   intent(in) params
  !f2py   intent(in) flags
  !f2py   intent(in) n
  !f2py   intent(in) ml_lg
  
  !f2py   intent(out) mcmc_mod
  !f2py   depend(params) mcmc_mod
  !f2py   depend(n) mcmc_mod
  
  lg_mcmc = .true.
  
  if (ml_lg == 1) then
     lg_multiline = .true.
  else
     lg_multiline = .false.
     n_lines=1
  end if

! dust_geometry%r_ratio = mcmc_v_min/mcmc_v_max

!-------------
  !set up multiline case for mcmc bayesian modelling
  if (lg_multiline) then
     open(54,file = 'multiline.in')
     read(54,*) n_lines, lg_multiline_fixdust,lg_multiline_fixgas
     !if running multiple lines then read in the relevant file names
     n_bins_multiline = 0
     do ii=1,n_lines
        read(54,*)  input_folder(ii), data_file
        open(53,file= trim(input_folder(ii)) // trim(data_file))
        read(53,*) obs_data%n_data
        n_bins_multiline = n_bins_multiline+obs_data%n_data
     end do
     allocate(multiline_profile_array(n_bins_multiline,2))
     close(53)

     !multiline_count will track how many bins have been recorded in the final multiline profile array
     multiline_count = 0
     rewind(unit=54)
     read(54,*)
  end if


  do i_line = 1,n_lines
     print*,'line no',i_line

     !update the number of observed data points for each line
     if (lg_multiline) then
        read(54,*)  input_folder(i_line), data_file
        open(53,file= trim(input_folder(i_line)) // trim(data_file))
        read(53,*) obs_data%n_data
        input_prefix = trim(input_folder(i_line))
     end if
     close(53)

     !!adjust this specification for multiple species or a grain size range
     !the below is only applicable for a single grain size
     if ((.not. lg_multiline_fixdust) .and. (i_line>1)) then
        allocate(dust%species(1))
     else if (i_line == 1) then
        allocate(dust%species(1))
     end if

     call read_input()
     dust%species(1)%nsizes=1
     dust%species(1)%power=1.0

     call read_bayesian_params(params,flags)
     call check_for_conflicts(flags)
     print*, 'dust geometry filling factor',dust_geometry%ff
     if (i_line>1) then
        deallocate(profile_array_data_bins)
        deallocate(mc_error_data_bins)
     end if

     call run_damocles()
     
     if (lg_multiline) then
        multiline_profile_array(1+multiline_count:multiline_count+obs_data%n_data,1) = profile_array_data_bins(:)
        multiline_profile_array(1+multiline_count:multiline_count+obs_data%n_data,2) = mc_error_data_bins(:)
        multiline_count= multiline_count + obs_data%n_data
     end if

  end do
     
!!!!!!!!!!!!!!!!!!!!!!!
  
  if (lg_multiline) then
     mcmc_mod(:,1) = multiline_profile_array(:,1)
     mcmc_mod(:,2) = multiline_profile_array(:,2)
  else
     mcmc_mod(:,1) = profile_array_data_bins
     mcmc_mod(:,2) = mc_error_data_bins
  end if
  
  if (lg_multiline) deallocate(multiline_profile_array)
  deallocate(profile_array_data_bins)
  deallocate(mc_error_data_bins)

  
  if (lg_multiline_fixdust) then
     deallocate(grid_cell)
     deallocate(mothergrid%x_div)
     deallocate(mothergrid%y_div)
     deallocate(mothergrid%z_div)
     deallocate(dust%species)
  end if
  close(54)
  
end subroutine run_damocles_wrap

