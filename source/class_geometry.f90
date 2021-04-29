!-----------------------------------------------------------------------!
!  this module declares the geometry derived type objects               !
!       - this is the general geoemtry description for gas or for dust  !
!         and includes parameters such as maximum velocity, r_max etc.  !
!  includes a subroutine to check whether dust is clumped or not        !
!-----------------------------------------------------------------------!

module class_geometry

    implicit none

    type geometry_obj
        real        ::  clumped_mass_frac   !mass fraction of dust located in clumpes. 0 indicates no clumping.
        real        ::  r_ratio             !ratio between inner and outer radii (r_in/r_out)
        real        ::  v_max               !maximum velocity (km/s) at r_out
        real        ::  v_min               !minimum velocity (km/s) (determined from geometry if coupled to radius)
        real        ::  v_power             !value of l where velocity profile is v~r^l
        real        ::  v_prob_indx         !power of probability distribution of v in case where independent of radius (different to v_power)
        real        ::  rho_power           !value of q where density profile is rho~r^-q
        real        ::  emis_power          !value of b where emissivity~r^-qb
        real        ::  clump_power         !value of q where the number density of clumps is distributed as n~r^-q
        real        ::  ff                  !filling factor of clumps by volume
        real        ::  den_con             !density contrast at inner radius
        real        ::  r_min,r_max         !inner and outer radii of distribution (e15cm)
        real        ::  r_min_cm,r_max_cm   !inner and outer radii of distribution (cm)
        real        ::  rho_in              !density of smooth medium at inner radius (in non-clumped case)
        real        ::  rho_clump           !denisty (constant) of an individual clump

        integer     ::  n_clumps            !number of clumps

        logical     ::  lg_clumped          !is the medium clumped?

        character(len=9) ::   type          !e.g. shell or arbitrary or torus
    end type

    type(geometry_obj) gas_geometry, dust_geometry

    contains

    !this subroutine sets the dust clumping logical to true if the clumped dust mass fraction is between 0 and 1
    subroutine check_dust_clumped()
        !test 0 or 1 entered
        if (dust_geometry%clumped_mass_frac < 0 .or. dust_geometry%clumped_mass_frac > 1) then
            print*, "aborted - please enter a value between 0 and 1 for the mass fraction of dust in clumps."
            stop
        end if
        dust_geometry%lg_clumped = .false.
        if (dust_geometry%clumped_mass_frac > 0) then
            if (dust_geometry%ff >0) then
                dust_geometry%lg_clumped = .true.
            else
                print*,"you have requested that a non-zero clumped dust mass fraction &
                & but have not specified a clump volume filling factor.  please set a non-zero volume &
                & filling factor between 0 and 1 or set the clumped dust mass fraction to 0.  aborted."
                stop
            end if
        end if

    end subroutine

end module class_geometry
