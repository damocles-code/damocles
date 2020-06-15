!----------------------------------------------------------------------------!
!  this subroutine is only called if the python mcmc wrapper is being used.  !
!  the logical lg_mcmc is set to true and input parameters passed from the   !
!  python routine.                                                           !
!----------------------------------------------------------------------------!
subroutine run_damocles_wrap(params,flags,n,n_lines,mcmc_mod)  

    use globals
    use input
    use class_dust
    use initialise
    use vector_functions
    use driver
    use bayesian_input

  implicit none

  real :: params(21,6)  !an array to store all possible variable parameters for each line (max 6 lines)
  real :: flags(21,6)   !an array to indicate which parameters are being varied by python (max 6 lines)
  integer :: n          !the size of returned line profiles array (containing all lines concatenated)
  integer :: n_lines    !the number of lines being modelled
  real :: mcmc_mod(n,2) !the array to store the modelled line profiles (concatenated)
  character(30) :: input_folder(6) !array to store the names of the input folders for the lines
  character(10) :: line_name(6)    !names of the lines

  !f2py   intent(in) params
  !f2py   intent(in) flags
  !f2py   intent(in) n
  !f2py   intent(in) ml_lg
  
  !f2py   intent(out) mcmc_mod
  !f2py   depend(params) mcmc_mod
  !f2py   depend(n) mcmc_mod
  
  lg_mcmc = .true.
  
  if (n_lines > 1) then
     lg_multiline = .true.
  else if (n_lines == 1) then
     lg_multiline = .false.
  else
     print*, 'You have requested a negative number of lines. Aborting.'
     STOP
  end if

! dust_geometry%r_ratio = mcmc_v_min/mcmc_v_max

!-------------
  !set up multiline case for mcmc bayesian modelling
!  if (lg_multiline) then
     open(54,file = 'line_list.in')
     read(54,*)
     read(54,*)
     read(54,*)
     read(54,*)
     !if running multiple lines then read in the relevant file names
     n_bins_multiline = 0
     do i_line=1,n_lines
        read(54,*)  line_name(i_line), input_folder(i_line), lg_multiline_fixgas(i_line)
        open(53,file= trim(input_folder(i_line)) // '/line.in')
        read(53,*) obs_data%n_data
        n_bins_multiline = n_bins_multiline+obs_data%n_data
     end do
     allocate(multiline_profile_array(n_bins_multiline,2))
     close(53)

     !multiline_count will track how many bins have been recorded in the final multiline profile array
     multiline_count = 0
     rewind(unit=54)
     read(54,*)
!  end if


  do i_line = 1,n_lines
     !update the number of observed data points for each line
     if (lg_multiline) then
        open(53,file= trim(input_folder(i_line)) // '/line.in')
        read(53,*) obs_data%n_data
        !input_prefix = trim(input_folder(i_line)) // '/'
     end if
     input_prefix = trim(input_folder(i_line)) // '/'
     close(53)

     !!adjust this specification for multiple species or a grain size range
     !the below is only applicable for a single grain size
     if (i_line == 1) then
        allocate(dust%species(1))
     end if

     call read_input()
     dust%species(1)%nsizes=1
     dust%species(1)%power=1.0

     call read_bayesian_params(params,flags)
     call check_for_conflicts(flags)

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
  
  deallocate(multiline_profile_array)
  deallocate(profile_array_data_bins)
  deallocate(mc_error_data_bins)

  
  deallocate(grid_cell)
  deallocate(mothergrid%x_div)
  deallocate(mothergrid%y_div)
  deallocate(mothergrid%z_div)
  deallocate(dust%species)

  close(54)
  
end subroutine run_damocles_wrap

