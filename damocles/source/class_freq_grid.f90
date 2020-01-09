!----------------------------------------------------------------------------------------!
!  this module declares the frequency grid derived type object                           !
!  subroutine establishes the frequency grid that will store the resultant line profile  !
!----------------------------------------------------------------------------------------!

module class_freq_grid

    use class_line
    use class_geometry
    use class_dust
    use class_obs_data

    implicit none

    type freq_grid_obj
        integer ::  n_bins                      !number of frequency bins

        real    ::  fmax                        !maximum frequency in frequency grid
        real    ::  fmin                        !minimum frequency in frequency grid
        real    ::  bin_width                   !size of a step in the (linear) frequency grid

        real,dimension(:,:),allocatable :: bin          !array of frequency bins
        real,dimension(:),allocatable   :: lambda_bin   !array of wavelength bins
        real,dimension(:),allocatable   :: vel_bin      !array of velocity bins

    end type freq_grid_obj

    type(freq_grid_obj) nu_grid

contains

    !this subroutine constructs a linear frequency grid.
    !as packets escape the grid, they will be added (according to their weight) to a bin in the frequency grid.
    !it is constructed at the start of the simulation in order to store packet data cumulativley as the rt progreses.
    subroutine construct_freq_grid()

      if (lg_data) then
         call read_obs_data()
         
         !calculate observed data frequency grid
         !note that these bins will be unequal and will be scaled later to account for this
         obs_data%freq = line%frequency/(1+obs_data%vel*10**3/c)
      end if


        if (.not. lg_mcmc) print*, 'constructing frequency grid...'

        allocate(nu_grid%bin(nu_grid%n_bins,2))
        allocate(nu_grid%lambda_bin(nu_grid%n_bins-1))
        allocate(nu_grid%vel_bin(nu_grid%n_bins-1))

        !set maximum and minimum frequency range for bins in frequency grid
        !if doing a doublet, then the min and max are set to be as large as possible based on the wavelengths of interest
        !maximum is set to be (arbitrarily) 1.2 times larger than the frequency obtained by shifting by v_max towards the blue
        !minimum is set to be (arbitrarily) 1.2 times smaller than the frequency obtained by shifting by v_max towards the red
        !this is to allow for frequency being shifted beyond the theoretical for a single scattering event due to multiple scatterings
        !!this array could be simplified
        if (lg_doublet) then
              if (line%doublet_wavelength_2>line%doublet_wavelength_1) then
                 nu_grid%fmax=((c*10**9/line%doublet_wavelength_1)/(1-max(dust_geometry%v_max,2000.0)*2.0*10**3/c))
                 nu_grid%fmin=((c*10**9/line%doublet_wavelength_2)/(1+max(dust_geometry%v_max,2000.0)*2.0*10**3/c))
              else
                 nu_grid%fmax=((c*10**9/line%doublet_wavelength_2)/(1-max(dust_geometry%v_max,2000.0)*2.0*10**3/c))
                 nu_grid%fmin=((c*10**9/line%doublet_wavelength_1)/(1+max(dust_geometry%v_max,2000.0)*2.0*10**3/c))
           end if
        else
                nu_grid%fmax=(line%frequency/(1-max(dust_geometry%v_max,gas_geometry%v_max)*2.0*10**3/c))
                nu_grid%fmin=(line%frequency/(1+max(dust_geometry%v_max,gas_geometry%v_max)*2.0*10**3/c))                !as above
        end if
        nu_grid%bin_width=(nu_grid%fmax-nu_grid%fmin)/nu_grid%n_bins

        !calculate the frequency bins for the resultant line profile to be stored in
        do ii=1,nu_grid%n_bins
            nu_grid%bin(ii,1)=nu_grid%fmin+((ii-1)*nu_grid%bin_width)
            nu_grid%bin(ii,2)=nu_grid%fmin+((ii)*nu_grid%bin_width)
        end do

        !calculate the wavelength and velocity space arrays
        do ii=1,nu_grid%n_bins-1
            nu_grid%lambda_bin(ii)=(c*10**9)*(0.5/nu_grid%bin(ii,1)+0.5/nu_grid%bin(ii+1,1))
            nu_grid%vel_bin(ii)=(c*1e-3*(nu_grid%lambda_bin(ii)**2-line%wavelength**2)/(line%wavelength**2+nu_grid%lambda_bin(ii)**2))
        end do

    end subroutine construct_freq_grid

end module class_freq_grid
