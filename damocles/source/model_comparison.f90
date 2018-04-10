module model_comparison

    use globals
    use class_line
    use class_geometry
    use class_dust
    use class_grid
    use class_freq_grid
    use class_obs_data

    implicit none

contains

    subroutine calculate_chi_sq()

        open(81,file='output/modelled_line_fine_res.out')
        do ii = 1,nu_grid%n_bins-1
            bin_id = minloc(nu_grid%vel_bin(ii)-obs_data%vel(:),1,(nu_grid%vel_bin(ii)-obs_data%vel(:))>0)
            if (bin_id /= 0) write(81,*) nu_grid%vel_bin(ii), profile_array(ii), bin_id, obs_data%vel(bin_id)
        end do
        close(81)

        !scale modelled fluxes and errors to data flux using observed line flux
        !tot_flux/(n_recorded_packets) = final average energy per packet.  Multiply by weight to get energy at end.
        profile_array_data_bins = line%tot_flux*total_weight_data_bins/n_recorded_packets
        mc_error_data_bins = line%tot_flux*(square_weight_data_bins**0.5)/n_recorded_packets

        !check whether velocity bins should be excluded from chi squared calculation (due to e.g. narrow line contamination)
        obs_data%exclude(:) = .false.
        do ii = 1,obs_data%n_data
            do jj = 1,no_exclusion_zones
                if ((obs_data%vel(ii)>exclusion_zone(jj,1)) .and. (obs_data%vel(ii)<exclusion_zone(jj,2))) then
                    obs_data%exclude(ii) = .true.
                end if
            end do
        end do

        !calculate chi_sq
        !!check chi_sq scaling of chi squared calculation
       chi_sq=0
        n_data_points = 0
        do ii = 1,obs_data%n_data
            if (obs_data%exclude(ii) .eqv. .false.) then
               n_data_points = n_data_points + 1
                if (isnan(mc_error_data_bins(ii)) .or. isnan(obs_data%error(ii))) then
                    cycle
                else
                    chi_sq = chi_sq+((profile_array_data_bins(ii)-1e15*obs_data%flux(ii))**2/((1e15*obs_data%error(ii))**2+mc_error_data_bins(ii)**2))
                end if
            end if
        end do

        print*, 'vmax',dust_geometry%v_max
        print*, 'R_rat',dust_geometry%r_ratio
        print*, 'rho',dust_geometry%rho_power
        print*, 'mdust',dust%mass
        print*, 'a',dust%species(1)%amin
!        print*, 'chi sq - ', chi_sq
!        print*,'number of data points',n_data_points

        !optimise scale factor to give best fit by trying a range of scale factors between 0.6 and 1.5 times the initial factor above
        !!!temp comment out
!        do jj =1,20
!            sf = scale_factor*(0.5+0.05*jj)
!            chi_sq_new=0
!            do ii = 1,obs_data%n_data
!                if (obs_data%exclude(ii) .eqv. .false.) then
!                    chi_sq_new = chi_sq_new+((profile_array_data_bins(ii)*sf-obs_data%flux(ii))/error)**2
!                end if
!            end do
!            if (chi_sq_new<chi_sq) then
!                scale_factor_final=sf
!                chi_sq=chi_sq_new
!            end if
!        end do
!        profile_array_data_bins=profile_array_data_bins*scale_factor_final

        !write out rebinned and rescaled modelled line to file

    end subroutine

end module model_comparison
