MODULE class_freq_grid

    USE class_line
    USE class_geometry
    USE class_dust



    IMPLICIT NONE

    TYPE freq_grid_obj

        INTEGER ::  n_bins                      !number of frequency bins
        REAL    ::  fmax                        !maximum frequency in frequency grid
        REAL    ::  fmin                        !minimum frequency in frequency grid
        REAL    ::  bin_width                   !size of a step in the (linear) frequency grid
        REAL,DIMENSION(:,:),ALLOCATABLE :: bin

    END TYPE freq_grid_obj

    TYPE(freq_grid_obj) nu_grid

contains

    !This subroutine constructs a linear frequency grid.
    !As packets escape the grid, they will be added (according to their weight) to a bin in the frequency grid.
    !It is constructed at the start of the simulation in order to store packet data cumulativley as the RT progreses.
    SUBROUTINE construct_freq_grid()

        ALLOCATE(nu_grid%bin(nu_grid%n_bins,2))

        IF (lg_doublet) THEN
            IF (line%doublet_wavelength_2>line%doublet_wavelength_1) THEN
                nu_grid%fmax=5*((c*10**9/line%doublet_wavelength_1)/(1-MAX(dust_geometry%v_max,2000.0)*10**3/c))        !arbitrary factor of 5 to compensate for multiple scatterings
                nu_grid%fmin=0.2*((c*10**9/line%doublet_wavelength_2)/(1+MAX(dust_geometry%v_max,2000.0)*10**3/c))        !as above
            ELSE
                nu_grid%fmax=5*((c*10**9/line%doublet_wavelength_2)/(1-MAX(dust_geometry%v_max,2000.0)*10**3/c))        !arbitrary factor of 5 to compensate for multiple scatterings
                nu_grid%fmin=0.2*((c*10**9/line%doublet_wavelength_1)/(1+MAX(dust_geometry%v_max,2000.0)*10**3/c))        !as above
            END IF
        ELSE
            nu_grid%fmax=5*(line%frequency/(1-MAX(dust_geometry%v_max,2000.0)*10**3/c))  !arbitrary factor of 4 to compensate for multiple scatterings
                                                                                                !shifting nu beyond maxnu for 1 scattering event
            nu_grid%fmin=0.2*(line%frequency/(1+MAX(dust_geometry%v_max,2000.0)*10**3/c))                !as above
        END IF
        nu_grid%bin_width=(nu_grid%fmax-nu_grid%fmin)/nu_grid%n_bins

        !CALCULATE FREQUENCY BINS AND WRITE TO FILE

        DO ii=1,nu_grid%n_bins
            nu_grid%bin(ii,1)=nu_grid%fmin+((ii-1)*nu_grid%bin_width)
            nu_grid%bin(ii,2)=nu_grid%fmin+((ii)*nu_grid%bin_width)
        END DO

    END SUBROUTINE construct_freq_grid

END MODULE class_freq_grid
