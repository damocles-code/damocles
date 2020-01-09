
module electron_scattering

    use globals
    use class_geometry
    use class_dust
    use class_freq_grid
    use class_grid

    implicit none

    real                 ::  es_const
    real                 ::  r_max_es
    real                 ::  lum_q
    real                 ::  maxwell_sigma

    real, parameter      ::  sigma_t=6.652e-25        !(cm2)
    real, parameter      ::  q_halpha_5=6.71e-25     !at 5000k
    real, parameter      ::  q_halpha_10=3.56e-25    !at 10000k
    real, parameter      ::  q_halpha_20=1.83e-24    !at 20000k

contains

    subroutine n_e_const()

        if (es_temp==5000) then
            lum_q=l_halpha/q_halpha_5
        else if (es_temp==10000) then
            lum_q=l_halpha/q_halpha_10
        else if (es_temp==20000) then
            lum_q=l_halpha/q_halpha_20
        end if


        if (lg_es) then
            if ((3-2*gas_geometry%rho_power==0) .or. (1-gas_geometry%rho_power==0)) then
                print*,'you have selected a density profile with exponent 1.5 or 1.0 - need an alternative calculation in this case'
                stop
            end if
        end if

        !at some point, include coding for separate line emitting region and gas region
        r_max_es=gas_geometry%r_max

        maxwell_sigma=((es_temp*1.51563e7)**0.5)/1000
        
        if (lg_es) then
           if (gas_geometry%type == "shell") then
              es_const=(lum_q**0.5)*(((3-2*gas_geometry%rho_power)/(4*pi))/(((r_max_es*1e15)**(3-2*gas_geometry%rho_power)-(gas_geometry%r_min*1e15)**(3-2*gas_geometry%rho_power))))**0.5
              do iG=1,mothergrid%tot_cells
                 if ((grid_cell(iG)%r*1e-15 < gas_geometry%r_max) .and. (grid_cell(iG)%r*1e-15 > gas_geometry%r_min)) then
                    grid_cell(iG)%n_e=es_const*(grid_cell(iG)%r**(-gas_geometry%rho_power))*1e20
                 end if

              end do
           else
              print*,'provision for electron scattering with a non-shell emissivity distribution has not yet been included.Aborted.'
              stop
           end if
        end if

        !print*,'av e- density',(es_const*1e20*((1e15)**(-gas_geometry%rho_power))*((r_max_es)**(3-gas_geometry%rho_power)-(gas_geometry%r_min)**(3-gas_geometry%rho_power)))/((3-dust_geometry%rho_power)*((r_max_es)**3-(gas_geometry%r_min)**3))
        !print*,''
        !print*,'e- optical depth',es_const*6.6e-5*((1e15*r_max_es)**(1-gas_geometry%rho_power)-(1e15*gas_geometry%r_min)**(1-gas_geometry%rho_power))/(1-gas_geometry%rho_power)

    end subroutine

end module
