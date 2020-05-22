module class_obs_data

    use globals

    implicit none

    type obs_profile
        integer             :: n_data                               !no of divisions in observed profile (i.e. no of lines to read in)
        real,allocatable    :: vel(:)                               !array of velocities of observed profile
        real,allocatable    :: freq(:)                              !array of frequencies of observed profile
        real,allocatable    :: flux(:)                              !array of fluxes of observed profile
        logical,allocatable :: exclude(:)                           !logical array that with certain elements assigned true if they are excluded from chi sq calcn.
        real,allocatable    :: error(:)                            !error on measurement used in chi_sq calculation
        real                :: mean_freq_bin                         !average width of a frequency bin in the obs data (used to rescale the binned profile to account for different width bins)
    end type

    type(obs_profile) obs_data

    integer              :: bin_id_low,bin_id_high,bin_id           !ids of modelled velocity bins that are used to interpolate to obs data bins
    real                 :: scale_factor,sf,scale_factor_final      !scale factors used to scale model to maximise chi sq
    real                 :: chi_sq_new                              !updated chi_sq
    integer              :: no_exclusion_zones                      !no of regions to exclude from chi sq calulation
    real,allocatable     :: exclusion_zone(:,:)                     !upper and lower bounds of excluded zone in velocity space
    real                 :: err

contains

    subroutine read_obs_data()

        !read in observed data number of lines
        open(35,file=data_file)
        if (lg_multiline) then
           read(35,*)
        else
           read(35,*) obs_data%n_data, err
        end if

        !allocate size of data array accordingly
        allocate(obs_data%vel(obs_data%n_data))
        allocate(obs_data%flux(obs_data%n_data))
        allocate(obs_data%exclude(obs_data%n_data))
        allocate(obs_data%error(obs_data%n_data))
        allocate(obs_data%freq(obs_data%n_data))

        obs_data%error = err
        
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
end module class_obs_data
