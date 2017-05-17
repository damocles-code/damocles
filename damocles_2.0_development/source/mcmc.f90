module mcmc

    use globals
    use class_line
    use class_freq_grid
    use class_grid
    use electron_scattering
    use input
    use initialise
    use vector_functions
    use class_packet
    use radiative_transfer
    use model_comparison

    implicit none

    type mcmc_param_space
        integer :: n_points
        integer :: dimension
        real    :: v_max_bounds(2)
    end type

    type(mcmc_param_space) param_space

contains

    subroutine mcmc_handler()

        !final parameter list to be explored: max velocity, min velocity, density gradient, dust mass
        !!current list: max velocity

        open(53,file='input/MCMC.in')
        read(53,*)
        read(53,*) param_space%n_points
        close(53)

        open(54,file='output/trace.out')

        do ii=1,param_space%n_points
            call run_damocles()
            write(54,*) dust_geometry%v_max,chi_sq
        end do

        close(54)

    end subroutine

    subroutine mcmc_sampler()



        param_space%dimension = 1

        param_space%v_max_bounds(1) = 13000
        param_space%v_max_bounds(2) = 20000

        !velocity test run
        do ii=1,param_space%n_points
            call random_number(ran)
            dust_geometry%v_max = param_space%v_max_bounds(1)+ran*(param_space%v_max_bounds(2)-param_space%v_max_bounds(1))
            print*,dust_geometry%v_max
        end do


    end subroutine

end module mcmc
