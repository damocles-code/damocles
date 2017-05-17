module model_comparison

    use globals
    use class_line
    use class_geometry
    use class_dust
    use class_grid
    use class_freq_grid


    implicit none

    type obs_profile
        integer             :: n_data                               !no of divisions in observed profile (i.e. no of lines to read in)
        real,allocatable    :: vel(:)                               !array of velocities of observed profile
        real,allocatable    :: flux(:)                              !array of fluxes of observed profile
        logical,allocatable :: exclude(:)                           !logical array that with certain elements assigned true if they are excluded from chi sq calcn.
    end type

    type(obs_profile) obs_data,model_rebinned

    integer              :: bin_id_low,bin_id_high                  !ids of modelled velocity bins that are used to interpolate to obs data bins
    real                 :: scale_factor,sf,scale_factor_final      !scale factors used to scale model to maximise chi sq
    real                 :: chi_sq_new                              !updated chi_sq
    real,parameter       :: error=0.5e-15                           !error used in chi_sq calculation !!may want to make dynamic or non=const
    integer              :: no_exclusion_zones                      !no of regions to exclude from chi sq calulation
    real,allocatable     :: exclusion_zone(:,:)                     !upper and lower bounds of excluded zone in velocity space

contains

    subroutine read_in_data()

        !read in observed data number of lines
        open(35,file=data_file)
        read(35,*) obs_data%n_data

        !allocate size of data array accordingly
        allocate(obs_data%vel(obs_data%n_data))
        allocate(obs_data%flux(obs_data%n_data))

        !read in observed line data
        do ii = 1,obs_data%n_data
            read(35,*) obs_data%vel(ii),obs_data%flux(ii)
        end do
        close(35)

        !read in excluded zones
        open(36,file=data_exclusions_file)
        read(36,*) no_exclusion_zones
        allocate(exclusion_zone(no_exclusion_zones,2))

        do ii = 1,no_exclusion_zones
            read(36,*) exclusion_zone(ii,1),exclusion_zone(ii,2)
            if (exclusion_zone(ii,1)>exclusion_zone(ii,2)) print*,'warning: please re-order your exclusion zone limits to be in ascending order and re-run.'
        end do
        close(36)
    end subroutine

    subroutine calculate_chi_sq()

        !allocate space for new modelled line profile binned to same velocity bins as observed data
        allocate(model_rebinned%vel(obs_data%n_data))
        allocate(model_rebinned%flux(obs_data%n_data))
        allocate(model_rebinned%exclude(obs_data%n_data))

        !set model velocity bins to be the same as those of observed line profile
        model_rebinned%vel = obs_data%vel

        !interpolate between nearest modelled velocity bins to find new flux in observed data velocity bins
        do ii = 1,obs_data%n_data
            bin_id_high=minloc(nu_grid%vel_bin(:)-model_rebinned%vel(ii),1,(nu_grid%vel_bin(:)-model_rebinned%vel(ii)>0))
            bin_id_low=maxloc(nu_grid%vel_bin(:)-model_rebinned%vel(ii),1,(nu_grid%vel_bin(:)-model_rebinned%vel(ii)<0))
            model_rebinned%flux(ii) = line%initial_energy*(profile_array(bin_id_low)+((profile_array(bin_id_high)-profile_array(bin_id_low))*((model_rebinned%vel(ii)-nu_grid%vel_bin(bin_id_low))/(nu_grid%vel_bin(bin_id_high)-nu_grid%vel_bin(bin_id_low)))))
        end do

        !scale modelled fluxes to data flux using peak observed line flux
        scale_factor=line%peak_flux/(sum(model_rebinned%flux(maxloc(model_rebinned%flux,1)-5:maxloc(model_rebinned%flux,1)+5))/size(model_rebinned%flux(maxloc(model_rebinned%flux,1)-5:maxloc(model_rebinned%flux,1)+5)))

        !check whether velocity bins should be excluded from chi squared calculation (due to e.g. narrow line contamination)
        model_rebinned%exclude = .false.
        do ii = 1,obs_data%n_data
            do jj = 1,no_exclusion_zones
                if ((model_rebinned%vel(ii)>exclusion_zone(jj,1)) .and. (model_rebinned%vel(ii)<exclusion_zone(jj,2))) then
                    model_rebinned%exclude(ii) = .true.
                end if
            end do
        end do

        !calculate chi_sq
        do ii = 1,obs_data%n_data
            if (model_rebinned%exclude(ii) .eqv. .false.) then
                chi_sq = chi_sq+((model_rebinned%flux(ii)*scale_factor-obs_data%flux(ii))/error)**2
            end if
        end do

        !optimise scale factor to give best fit by trying a range of scale factors between 0.6 and 1.5 times the initial factor above
        do jj =1,20
            sf = scale_factor*(0.5+0.05*jj)
            chi_sq_new=0
            do ii = 1,obs_data%n_data
                if (model_rebinned%exclude(ii) .eqv. .false.) then
                    chi_sq_new = chi_sq_new+((model_rebinned%flux(ii)*sf-obs_data%flux(ii))/error)**2
                end if
            end do
            if (chi_sq_new<chi_sq) then
                scale_factor_final=sf
                chi_sq=chi_sq_new
            end if
        end do
        model_rebinned%flux=model_rebinned%flux*scale_factor_final

        !write out rebinned and rescaled modelled line to file
        if (.not. lg_mcmc) then
            print*,'chi squared',chi_sq
            open(37,file='line.out')
            do ii = 1,obs_data%n_data
                write(37,*) model_rebinned%vel(ii),model_rebinned%flux(ii)
            end do
            close(37)
        end if

    end subroutine


end module model_comparison
