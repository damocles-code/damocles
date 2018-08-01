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
  
  real :: params(21)
  real :: flags(21)
  integer :: n
  integer :: ml_lg
  real :: mcmc_mod(n,2)
  
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
  
  !!adjust this specification for multiple species or a grain size range
  !the below is only applicable for a single grain size
  allocate(dust%species(1))
  dust%species(1)%nsizes=1
  dust%species(1)%power=1.0

  if (.not. lg_decoupled) then
     if (flags(1) == 1) STOP 'You have requested dust and gas coupled. Do not vary gas max velocity or decouple dust and gas. Param #1.'
     if (flags(2) == 1) STOP 'You have requested dust and gas coupled. Do not vary gas max radius or decouple dust and gas. Param #2.'
     if (flags(3) == 1) STOP 'You have requested dust and gas coupled. Do not vary gas Rin/Rout ratio or decouple dust and gas. Param #3.'
     if (flags(5) == 1) STOP 'You have requested dust and gas coupled. Do not vary gas density profile or decouple dust and gas. Param #5.'
     if (flags(6) == 1) STOP 'You have requested dust and gas coupled. Do not vary gas emissivity profile or decouple dust and gas. Param #6.'
  end if

  if (.not. lg_vel_law) then
     if (flags(7) == 1) STOP 'You have not requested a velocity law independent of radius. Do not vary min gas velocity. Param #7.'
     if (flags(8) == 1) STOP 'You have not requested a velocity law independent of radius. Do not vary velocity probability law. Param #8.'
  else
     if (flags(4) == 1) STOP 'You have requested a velocity law independent of radius. Do\
     not vary velocity power (v propto r^x). Param #4.'
  end if
  
  if (.not. dust_geometry%lg_clumped) then
     if (flags(11) == 1) STOP 'You have requested to vary the dust clump filling factor but have not put any dust in clumps. Param #11.'
     if (flags(12) == 1) STOP 'You have requested to vary the dust clump number distribution but have not put any dust in clumps. Param #12.'
  end if
  
  if (.not. gas_geometry%lg_clumped) then
     if (flags(20) == 1) STOP 'You have requested to vary the gas clump filling factor but \
     have not requested clumped gas emission. Param #20.'
     if (flags(21) == 1) STOP 'You have requested to vary the gas clump number distribution\
     but have not requested clumped gas emission. Param #21.'
  end if
  
  if (flags(1) == 1) gas_geometry%v_max = params(1)*1000
  if (flags(2) == 1) gas_geometry%r_max = params(2)
  if (flags(3) == 1) gas_geometry%r_ratio = params(3)
  if (flags(4) == 1) gas_geometry%v_power = params(4)
  if (flags(5) == 1) gas_geometry%rho_power = params(5)
  if (flags(6) == 1) gas_geometry%emis_power = params(6)
  if (flags(7) == 1) gas_geometry%v_min = params(7)*1000
  if (flags(8) == 1) gas_geometry%v_prob_indx = params(8)
  if (flags(9) == 1) dust%mass = 10**params(9)
  if (flags(10) == 1) dust_geometry%clumped_mass_frac = params(10)
  if (flags(11) == 1) dust_geometry%ff = params(10)
  if (flags(12) == 1) dust_geometry%clump_power = params(12)
  if (flags(13) == 1) dust_geometry%v_max = params(13)*1000
  if (flags(14) == 1) dust_geometry%r_max = params(14)
  if (flags(15) == 1) dust_geometry%r_ratio = params(15)
  if (flags(16) == 1) dust_geometry%v_power = params(16)
  if (flags(17) == 1) gas_geometry%rho_power = params(17)
  if (flags(18) == 1) dust%species(1)%amin = 10**params(18)
  if (flags(18) == 1) dust%species(1)%amax = 10**params(18)
  if (flags(19) == 1) line%doublet_ratio = params(19)
  if (flags(20) == 1) gas_geometry%ff = params(20)
  if (flags(21) == 1) gas_geometry%clump_power = params(21)

  ! dust_geometry%r_ratio = mcmc_v_min/mcmc_v_max

  !!!!!!!!!!!!!!!!!!!!!!
  !set up multiline case for mcmc bayesian modelling
  if (lg_multiline) then
     open(54,file = input_prefix // 'multiline.in')
     read(54,*) n_lines, lg_multiline_fixdust
     
     n_bins_multiline = 0
     !if running multiple lines then read in the relevant file names
     do ii=1,n_lines
        read(54,*)  input_file, data_file
        open(53,file= input_prefix // trim(data_file))
        read(53,*) obs_data%n_data
        n_bins_multiline = n_bins_multiline+obs_data%n_data
     end do
     allocate(multiline_profile_array(n_bins_multiline,2))
     close(53)
     
     !multiline_count will track how many bins have been recorded in the final multiline profile array
     multiline_count = 0
     !  end if
     
     rewind(unit=54)
     read(54,*)
  else
     n_lines = 1
  end if
  
  !perform MCRT for each requested line
  !!!fix this bit for multiple lines!

  print*,'0'
  do i_line = 1,n_lines
     print*,'1'  
     !if only a single line then exit the multiline loop
     if ((i_line > 1) .and. (.not. lg_multiline)) exit
     print*,'2'
     if (.not. (lg_multiline_fixdust .and. i_line>1)) then
        dust_geometry%rho_power = mcmc_rho(i_line)
        dust_geometry%clump_power = mcmc_clump
        dust_geometry%r_ratio = mcmc_v_min(i_line)/mcmc_v_max
        dust_geometry%ff = mcmc_ff
        !mothergrid%n_cells(1) = mcmc_n_divs
     end if
     print*,'3'
     gas_geometry%r_ratio=mcmc_v_min(i_line)/mcmc_v_max
     gas_geometry%rho_power=mcmc_rho(i_line)
     gas_geometry%v_max=mcmc_v_max*10000
     
     print*,'4'
     if (lg_multiline) then
        !read in name of input file again
        read(54,*) input_file, data_file, data_exclusions_file, dust_file, species_file, gas_file
        
        input_file = input_prefix // trim(input_file)
        data_file = input_prefix // trim(data_file)
        data_exclusions_file = input_prefix // trim(data_exclusions_file)
        dust_file = input_prefix // trim(dust_file)
        species_file = input_prefix // trim(species_file)
        gas_file = input_prefix // trim(gas_file)
        print*,'5'
     end if
     call run_damocles()
     print*,'6'
     
     if (lg_multiline) then
        multiline_profile_array(1+multiline_count:multiline_count+obs_data%n_data,1) = profile_array_data_bins(:)
        multiline_profile_array(1+multiline_count:multiline_count+obs_data%n_data,2) = mc_error_data_bins(:)
        deallocate(profile_array_data_bins)
        deallocate(mc_error_data_bins)
        multiline_count= multiline_count + obs_data%n_data
     end if
     
  end do
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
