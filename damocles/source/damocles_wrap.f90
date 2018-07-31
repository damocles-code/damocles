!----------------------------------------------------------------------------!
!  this subroutine is only called if the python mcmc wrapper is being used.  !
!  the logical lg_mcmc is set to true and input parameters passed from the   !
!  python routine.                                                           !
!----------------------------------------------------------------------------!
subroutine run_damocles_wrap(mcmc_v_max,mcmc_v_min,mcmc_rho,mcmc_mdust, &
 & mcmc_grain_size,mcmc_doublet_ratio,mcmc_ff,mcmc_clump,n,ml_lg,mcmc_mod)

use globals
use input
use class_dust
use initialise
use vector_functions
use driver

implicit none

real :: mcmc_v_max
real :: mcmc_v_min(2)
real :: mcmc_rho(2)
real :: mcmc_mdust
real :: mcmc_grain_size
real :: mcmc_doublet_ratio
real :: mcmc_ff
real :: mcmc_clump
integer :: n
integer :: ml_lg
real :: mcmc_mod(n,2)

!f2py   intent(in) mcmc_v_min
!f2py   intent(in) mcmc_v_max
!f2py   intent(in) mcmc_rho
!f2py   intent(in) mcmc_mdust
!f2py   intent(in) mcmc_grain_size
!f2py   intent(in) mcmc_doublet_ratio
!f2py   intent(in) mcmc_ff
!f2py   intent(in) mcmc_clump
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

dust_geometry%v_max = mcmc_v_max*10000
  ! dust_geometry%rho_power = mcmc_rho
  ! dust_geometry%r_ratio = mcmc_v_min/mcmc_v_max
  dust%mass = 10**mcmc_mdust
  line%doublet_ratio = mcmc_doublet_ratio

  !!adjust this specification for multiple species or a grain size range
  !the below is only applicable for a single grain size
  allocate(dust%species(1))
  dust%species(1)%nsizes=1
  dust%species(1)%amin=10**mcmc_grain_size
  dust%species(1)%amax=10**mcmc_grain_size
  dust%species(1)%power=1.0
  !call run_damocles()


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
  end if

  rewind(unit=54)
  read(54,*)
 


  !perform MCRT for each requested line
  do i_line = 1,n_lines
  
           !if only a single line then exit the multiline loop
           if ((i_line > 1) .and. (.not. lg_multiline)) exit

           if (.not. (lg_multiline_fixdust .and. i_line>1)) then
             dust_geometry%rho_power = mcmc_rho(i_line)
             dust_geometry%clump_power = mcmc_clump
             dust_geometry%r_ratio = mcmc_v_min(i_line)/mcmc_v_max
             dust_geometry%ff = mcmc_ff
             !mothergrid%n_cells(1) = mcmc_n_divs
          end if

          gas_geometry%r_ratio=mcmc_v_min(i_line)/mcmc_v_max
          gas_geometry%rho_power=mcmc_rho(i_line)
          gas_geometry%v_max=mcmc_v_max*10000


           !read in name of input file again
           read(54,*) input_file, data_file, data_exclusions_file, dust_file, species_file, gas_file

           input_file = input_prefix // trim(input_file)
           data_file = input_prefix // trim(data_file)
           data_exclusions_file = input_prefix // trim(data_exclusions_file)
           dust_file = input_prefix // trim(dust_file)
           species_file = input_prefix // trim(species_file)
           gas_file = input_prefix // trim(gas_file)

           call run_damocles()


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

end subroutine
