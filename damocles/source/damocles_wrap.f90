!----------------------------------------------------------------------------!
!  this subroutine is only called if the python mcmc wrapper is being used.  !
!  the logical lg_mcmc is set to true and input parameters passed from the   !
!  python routine.                                                           !
!----------------------------------------------------------------------------!
subroutine run_damocles_wrap(mcmc_v_max,mcmc_v_min,mcmc_rho_index,mcmc_mdust,mcmc_grain_size,mcmc_doublet_ratio,n,mcmc_mod)

    use globals
    use input
    use class_dust
    use initialise
    use vector_functions
    use driver

    implicit none

    real :: mcmc_v_min
    real :: mcmc_v_max
    real :: mcmc_rho_index
    real :: mcmc_mdust
    real :: mcmc_grain_size
    real :: mcmc_doublet_ratio
    integer :: n
    real :: mcmc_mod(n,2)

!f2py   intent(in) mcmc_v_min
!f2py   intent(in) mcmc_v_max
!f2py   intent(in) mcmc_rho_index
!f2py   intent(in) mcmc_mdust
!f2py   intent(in) mcmc_grain_size
!f2py   intent(in) mcmc_doublet_ratio
!f2py   intent(in) n

!f2py   intent(out) mcmc_mod
!f2py   depend(mcmc_v_max) mcmc_mod
!f2py   depend(n) mcmc_mod

    lg_mcmc = .true.
    dust_geometry%v_max = mcmc_v_max*10000
    dust_geometry%rho_power = mcmc_rho_index
    dust_geometry%r_ratio = mcmc_v_min/mcmc_v_max
    dust%mass = 10**mcmc_mdust
    line%doublet_ratio = mcmc_doublet_ratio
    !!adjust this specification for multiple species or a grain size range
    !the below is only applicable for a single grain size
    allocate(dust%species(1))
    dust%species(1)%nsizes=1
    dust%species(1)%amin=10**mcmc_grain_size
    dust%species(1)%amax=10**mcmc_grain_size
    dust%species(1)%power=1.0
    call run_damocles()
    mcmc_mod(:,1) = profile_array_data_bins
    mcmc_mod(:,2) = mc_error_data_bins

    deallocate(profile_array_data_bins)
    deallocate(mc_error_data_bins)
end subroutine
