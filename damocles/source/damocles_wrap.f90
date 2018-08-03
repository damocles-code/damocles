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
  
  implicit none
  
  real :: params(21,6)
  real :: flags(21,6)
  integer :: n
  integer :: ml_lg
  real :: mcmc_mod(n,2)
  character(30) :: input_folder(6)
  
  !f2py   intent(in) params
  !f2py   intent(in) n
  !f2py   intent(in) ml_lg
  
  !f2py   intent(out) mcmc_mod
  !f2py   depend(mcmc_v_max) mcmc_mod
  !f2py   depend(n) mcmc_mod
  
  lg_mcmc = .true.
  
  if (ml_lg == 1) then
     lg_multiline = .true.
  else
     lg_multiline = .false.
  end if

! dust_geometry%r_ratio = mcmc_v_min/mcmc_v_max

!-------------
  !set up multiline case for mcmc bayesian modelling
  if (lg_multiline) then
     open(54,file = 'multiline.in')
     read(54,*) n_lines, lg_multiline_fixdust

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
  else
     n_lines = 1
  end if
!--------------

  !!adjust this specification for multiple species or a grain size range
  !the below is only applicable for a single grain size
  allocate(dust%species(1))
  dust%species(1)%nsizes=1
  dust%species(1)%power=1.0

  do i_line = 1,n_lines
     if (lg_multiline) input_prefix = input_folder(i_line)
     call read_input()
     call read_bayesian_params(i_line)
     call check_for_conflicts(i_line)
     call run_damocles()
  end do

     
     if (lg_multiline) then
        multiline_profile_array(1+multiline_count:multiline_count+obs_data%n_data,1) = profile_array_data_bins(:)
        multiline_profile_array(1+multiline_count:multiline_count+obs_data%n_data,2) = mc_error_data_bins(:)
        deallocate(profile_array_data_bins)
        deallocate(mc_error_data_bins)
        multiline_count= multiline_count + obs_data%n_data
     end if
     
  close(54)
  
!!!!!!!!!!!!!!!!!!!!!!!
  
  
  if (lg_multiline) then
     mcmc_mod(:,1) = multiline_profile_array(:,1)
     mcmc_mod(:,2) = multiline_profile_array(:,2)
  else
     mcmc_mod(:,1) = profile_array_data_bins
     mcmc_mod(:,2) = mc_error_data_bins
  end if
  
  if (lg_multiline) then
     deallocate(multiline_profile_array)
  else
     deallocate(profile_array_data_bins)
     deallocate(mc_error_data_bins)
  end if
  
  if (lg_multiline_fixdust) then
     deallocate(grid_cell)
     deallocate(mothergrid%x_div)
     deallocate(mothergrid%y_div)
     deallocate(mothergrid%z_div)
     deallocate(dust%species)
  end if
  
end subroutine run_damocles_wrap

subroutine read_bayesian_params()

  use globals
  use input
  use class_dust
  use initialise
  use vector_functions
  use driver

  if (flags(1,i_line) == 1) gas_geometry%v_max = params(1,i_line)*1000
  if (flags(2,i_line) == 1) gas_geometry%r_max = params(2,i_line)
  if (flags(3,i_line) == 1) gas_geometry%r_ratio = params(3,i_line)
  if (flags(4,i_line) == 1) gas_geometry%v_power = params(4,i_line)
  if (flags(5,i_line) == 1) gas_geometry%rho_power = params(5,i_line)
  if (flags(6,i_line) == 1) gas_geometry%emis_power = params(6,i_line)
  if (flags(7,i_line) == 1) gas_geometry%v_min = params(7,i_line)*1000
  if (flags(8,i_line) == 1) gas_geometry%v_prob_indx = params(8,i_line)
  if (flags(9,i_line) == 1) dust%mass = 10**params(9,i_line)
  if (flags(10,i_line) == 1) dust_geometry%clumped_mass_frac = params(10,i_line)
  if (flags(11,i_line) == 1) dust_geometry%ff = params(10,i_line)
  if (flags(12,i_line) == 1) dust_geometry%clump_power = params(12,i_line)
  if (flags(13,i_line) == 1) dust_geometry%v_max = params(13,i_line)*1000
  if (flags(14,i_line) == 1) dust_geometry%r_max = params(14,i_line)
  if (flags(15,i_line) == 1) dust_geometry%r_ratio = params(15,i_line)
  if (flags(16,i_line) == 1) dust_geometry%v_power = params(16,i_line)
  if (flags(17,i_line) == 1) gas_geometry%rho_power = params(17,i_line)
  if (flags(18,i_line) == 1) dust%species(1)%amin = 10**params(18,i_line)
  if (flags(18,i_line) == 1) dust%species(1)%amax = 10**params(18,i_line)
  if (flags(19,i_line) == 1) line%doublet_ratio = params(19,i_line)
  if (flags(20,i_line) == 1) gas_geometry%ff = params(20,i_line)
  if (flags(21,i_line) == 1) gas_geometry%clump_power = params(21,i_line)

end subroutine


subroutine check_for_conflicts()
  
  use globals
  use input
  use class_dust
  use initialise
  use vector_functions
  use driver

if (.not. lg_decoupled) then
     if (flags(1,i_line) == 1) then
        print*, 'You have requested dust and gas coupled. Do not vary gas max velocity or decouple dust and gas. Param #1.'
        STOP
     end if
     if (flags(2,i_line) == 1) then
        print*,'You have requested dust and gas coupled. Do not vary gas max rad us or decouple dust and gas. Param #2.'
        STOP
     end if
     if (flags(3,i_line) == 1) then
        print*, 'You have requested dust and gas coupled. Do not vary gas Rin/Rout ratio or decouple dust and gas. Param #3.'
        STOP
     end if
     if (flags(5,i_line) == 1) then
        print*,'You have requested dust and gas coupled. Do not vary gas density profile or decouple dust and gas. Param #5.'
        STOP
     end if
     if (flags(6,i_line) == 1) then
        print*, 'You have requested dust and gas coupled. Do not vary gas emissivity profile or decouple dust and gas. Param #6.'
        STOP
     end if
  end if

if (.not. lg_vel_law) then
     if (flags(7,i_line) == 1) then
        print*, 'You have not requested a velocity law independent of radius. Do not vary min gas velocity. Param #7.'
        STOP
     end if
     if (flags(8,i_line) == 1) then
        print*, 'You have not requested a velocity law independent of radius. Do not vary velocity probability law. Param #8.'
        STOP
     end if
  else
     if (flags(4,i_line) == 1) then
        print*, 'You have requested a velocity law independent of radius. Do not vary velocity power (v propto r^x). Param #4.'
        STOP
     end if
  end if

 if (.not. dust_geometry%lg_clumped) then
     if (flags(11,i_line) == 1) then
        print*, 'You have requested to vary the dust clump filling factor but have not put any dust in clumps. Param #11.'
        STOP
     end if
     if (flags(12,i_line) == 1) then
        print*, 'You have requested to vary the dust clump number distribution but have not put any dust in clumps. Param #12.'
        STOP
     end if
  end if

  if (.not. gas_geometry%lg_clumped) then
     if (flags(20,i_line) == 1) then
        print*, 'You have requested to vary the gas clump filling factor but have not requested clumped gas emission. Param #20.'
        STOP
     end if
     if (flags(21,i_line) == 1) then
        print*, 'You have requested to vary the gas clump number distribution but have not requested clumped gas emission. Param #21.'
        STOP
     end if
  end if


  if (lg_multiline_fixdust .and. i_line>1) then
     if ((flags(9,i_line) == 1) .or. &
     & (flags(10,i_line) == 1) .or. &
     & (flags(11,i_line) == 1) .or. &
     & (flags(12,i_line) == 1) .or. &
     & (flags(13,i_line) == 1) .or. &
     & (flags(14,i_line) == 1) .or. &
     & (flags(15,i_line) == 1) .or. &
     & (flags(16,i_line) == 1) .or. &
     & (flags(17,i_line) == 1) .or. &
     & (flags(18,i_line) == 1)) then
        print*,'You have requested to fix the dust distribution for multiple lines but have requested to vary a dust parameter in line number',i_line,' Aborting.'
        STOP
     end if
  end if

  if (lg_multiline_fixgas .and. i_line>1) then
     if ((flags(1,i_line) == 1) .or. &
     & (flags(2,i_line) == 1) .or. &
     & (flags(3,i_line) == 1) .or. &
     & (flags(4,i_line) == 1) .or. &
     & (flags(5,i_line) == 1) .or. &
     & (flags(6,i_line) == 1) .or. &
     & (flags(7,i_line) == 1) .or. &
     & (flags(8,i_line) == 1) .or. &
     & (flags(20,i_line) == 1) .or. &
     & (flags(21,i_line) == 1)) then
        print*,'You have requested to fix the dust distribution for multiple lines but have requested to vary a dust parameter in line number',i_line,' Aborting.'
        STOP
     end if
  end if

end subroutine
