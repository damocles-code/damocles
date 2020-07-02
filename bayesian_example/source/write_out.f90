!---------------------------------------------------------------------!
!  all properties of the model and the final modelled line profile(s) !
!  are written out here                                               !
!---------------------------------------------------------------------!

module initialise

    use globals
    use class_line
    use class_geometry
    use class_dust
    use class_grid
    use class_freq_grid
    use input

    implicit none

contains

    subroutine write_to_file()

        if (.not. lg_mcmc) print*,'writing to file...'

        !real number format, 6 characters, 2dp
101     format(a65'   'f10.2)

        !integer format, 4 characters
102     format(a65'   'i10)

        !scientific format, 5 characters, 2dp
103     format(a65'   'e10.2)

        !create folders dependent on date/time of run if requested
        !otherwise overwrite and store in main output folder
        if (lg_store_all) then
            call date_and_time(date,time)
            run_no_string = time(1:2) // '.' // time(3:4) // '.' // time(5:6)

            call system('mkdir -p output/output_' // date // '/run_' // run_no_string)
            call system('cp input/*.in output/output_' // date // '/run_' // run_no_string // '/.')

            !open output files to record resultant modelled line profile, input parameters and properties of model
            open(25,file='output/output_' // date // '/run_' // run_no_string // '/integrated_line_profile.out')
            open(26,file='output/output_' // date // '/run_' // run_no_string // '/multiple_los_line_profiles.out')
            open(27,file='output/output_' // date // '/run_' // run_no_string // '/model_properties.out')
            open(28,file='output/output_' // date // '/run_' // run_no_string // '/integrated_line_profile_binned.out')
            open(29,file='output/output_' // date // '/run_' // run_no_string // '/multiple_los_bins.out')
            open(30,file='output/output_' // date // '/run_' // run_no_string // '/los_line_profile.out')
            open(31,file='output/output_' // date // '/run_' // run_no_string // '/grid.out')
        else
            !open output files to record resultant modelled line profile, input parameters and properties of model
            open(25,file='output/integrated_line_profile.out')
            open(26,file='output/multiple_los_line_profiles.out')
            open(27,file='output/model_properties.out')
            open(28,file='output/integrated_line_profile_binned.out')
            open(29,file='output/multiple_los_bins.out')
            open(30,file='output/los_line_profile.out')
            open(31,file='output/grid.out')
        end if

        !write out modelled line profile

        do ii=1,nu_grid%n_bins-1
           if (lg_los) then
              write(30,*) nu_grid%lambda_bin(ii),nu_grid%vel_bin(ii),line%initial_energy*profile_array(ii)
           else
              write(25,*) nu_grid%lambda_bin(ii),nu_grid%vel_bin(ii),line%initial_energy*profile_array(ii)
           end if
        end do

         if (lg_multi_los) then
            do kk=1,n_angle_divs
               do jj=1,n_angle_divs
                  write(26,*) (kk-1)*n_angle_divs+jj, cos_theta_array(jj),phi_array(kk)    
                  write(29,*) (kk-1)*20+jj,kk,jj,phi_array(kk),cos_theta_array(jj)
                  do ii=1,nu_grid%n_bins-1
                     write(26,*) nu_grid%lambda_bin(ii),line%initial_energy*profile_los_array(ii,jj,kk)
                  end do
               end do
            end do
         end if
        
         !if using observational data, then it is generally specified in equal velocity bins (unequal frequency bins)
         !the data is collected in unequal frequency bins so must be rescaled accordingly
         if (lg_data) then
            obs_data%mean_freq_bin = (obs_data%freq(obs_data%n_data)-obs_data%freq(1))/obs_data%n_data
            do ii = 1,obs_data%n_data-1
               if (.not. lg_los) then
                  write(28,*) obs_data%vel(ii),profile_array_data_bins(ii)*obs_data%mean_freq_bin/(obs_data%freq(ii+1)-obs_data%freq(ii)),mc_error_data_bins(ii)*obs_data%mean_freq_bin/(obs_data%freq(ii+1)-obs_data%freq(ii))
!                  print*,profile_array_data_bins(ii)
               end if
            end do
         end if

        write(27,101)  'active rest wavelength:',line%wavelength
        write(27,*)

        write(27,*)    'dust geometry:'
        write(27,101)  'maximum velocity (km/s)',dust_geometry%v_max
        write(27,101)  'dust velocity index',dust_geometry%v_power
        write(27,101)  'dust density index',dust_geometry%rho_power
        write(27,101)  'Rmax dust (e15cm)',dust_geometry%r_max
        write(27,101)  'Rmin dust (e15cm)',dust_geometry%r_min
        write(27,101)  'dust R ratio (Rmax/Rmin)',dust_geometry%r_ratio

        write(27,*)
        write(27,*)    'gas geometry:'
        write(27,101)  'maximum velocity (km/s)',gas_geometry%v_max
        write(27,101)  'gas velocity index',gas_geometry%v_power
        write(27,101)  'gas density index',gas_geometry%rho_power
        write(27,101)  'Rmax gas (e15cm)',gas_geometry%r_max
        write(27,101)  'Rmin gas (e15cm)',gas_geometry%r_min
        write(27,101)  'gas r ratio',gas_geometry%r_ratio
        write(27,101)  'minimum velocity at Rmin (km/s)',(gas_geometry%r_min/gas_geometry%r_max)*gas_geometry%v_max

        write(27,*)
        write(27,*)    'properties of the dust:'
        write(27,103)  'Cext at line wavelength',dust%lambda_ext
        write(27,103)  'Csca at line wavelength',dust%lambda_sca
        write(27,101)  'albedo at line wavelength',dust%lambda_sca/dust%lambda_ext
        write(27,101)  'forward scattering parameter g at line wavelength',dust%lambda_g_param
        write(27,101)  'average optical depth at line wavelength',mothergrid%n_rho_dust_av*dust%lambda_ext*(dust_geometry%r_max_cm-dust_geometry%r_min_cm)
        write(27,101)  'average optical depth (absn) at line wavelength',mothergrid%n_rho_dust_av*dust%lambda_ext*(dust_geometry%r_max_cm-dust_geometry%r_min_cm)*(1-dust%lambda_sca/dust%lambda_ext)
        write(27,101)  'average optical depth at V band (547nm)',mothergrid%n_rho_dust_av*dust%lambda_ext_v*(dust_geometry%r_max_cm-dust_geometry%r_min_cm)

        if (dust%n_species > 1) then
            write(27,*)
            write(27,*)    'properties of inidividual dust species:'
            do ii=1,dust%n_species
                write(27,102) 'species number:',dust%species(ii)%id
                write(27,101) 'mass weight',dust%species(ii)%m_weight
                write(27,101) 'albedo at rest frame wavelength',dust%species(ii)%albedo(line%wav_bin)
            end do
        end if

        write(27,*)
        write(27,*)    'properties of the grid:'
        write(27,103)  'volume of grid cell (and therefore clump) (e42cm^3)',mothergrid%cell_vol
        write(27,103)  'grid cell width (x) in cm',mothergrid%cell_width(1)
        write(27,103)  'grid cell width (y) in cm',mothergrid%cell_width(2)
        write(27,103)  'grid cell width (z) in cm',mothergrid%cell_width(3)
        write(27,103)  'total analytical volume of supernova (e42cm^3)',tot_vol
        write(27,103)  'volume of total grid cells inside sn (e42cm)',no_active_cells*mothergrid%cell_vol
        write(27,103)  'average dust number density',mothergrid%n_rho_dust_av
        write(27,101)  'average optical depth of a cell at line wavelength',mothergrid%n_rho_dust_av*dust%lambda_ext*0.5*mothergrid%cell_width(1)*(3.0/(4.0*pi))**0.3333333
        write(27,101)  'average optical depth of a cell at V band (547nm)',mothergrid%n_rho_dust_av*dust%lambda_ext_v*0.5*mothergrid%cell_width(1)*(3.0/(4.0*pi))**0.3333333

        if (dust_geometry%lg_clumped) then
            write(27,103)  'effective spherical radius of clump (cm)',mothergrid%cell_width(1)*(3.0/(4.0*pi))**0.3333333
            write(27,101)  'effective spherical radius as a fraction of Rmax',(mothergrid%cell_width(1)/dust_geometry%r_max_cm)*(3.0/(4.0*pi))**0.3333333
            write(27,103)  'average dust density (including clumps) (g/cm3)',(dust%mass_grams*1e-14)/(tot_vol*1e28)

            write(27,*)
            write(27,*)    'dust masses used:'

            write(27,101)  'filling factor',dust_geometry%ff
            write(27,101)  'clumped mass fraction',dust_geometry%clumped_mass_frac
            write(27,103)  'mass of interclump medium (Msun)',dust%m_icm
            write(27,103)  'mass in clumps (Msun)',dust%m_clump*dust_geometry%n_clumps
            write(27,103)  'total dust mass (Msun)',dust%mass

            write(27,*)
            write(27,102)  'check number of clumps used:'
            write(27,102)  '      using - ',n_clumps
            write(27,102)  '      requested -',dust_geometry%n_clumps

            write(27,101)  'check dust mass in clumps: '
            write(27,103)  '      using - ',m_clumps_check
            write(27,103)  '      requested: ',dust%m_clump*dust_geometry%n_clumps

            write(27,103)  'check mass in icm:'
            write(27,103)  '        using - ',m_icm_check
            write(27,103)  '        requested - ',dust%m_icm

            write(27,103)  'check total dust mass (rho*v)'
            write(27,103)  '        using - ',m_tot_check
            write(27,103)  '        requested -',dust%mass
        else
            write(27,103)  'dust grain number density at Rin',dust_geometry%rho_in/dust%av_mgrain
            write(27,103)  'dust grain number density at Rout',((dust_geometry%rho_in)*(dust_geometry%r_min_cm/dust_geometry%r_max_cm)**dust_geometry%rho_power)/dust%av_mgrain

            write(27,*)
            write(27,*)  'dust mass check (calculated as rho*v):'
            write(27,103)  '      using - ',m_tot_check
            write(27,103)  '      requested -',dust%mass
        end if

        write(27,*)
        write(27,*)    'packet resolution:'
        write(27,102)  'total number of packets initialised',n_init_packets
        write(27,102)  'number of active (propagated) packets',n_init_packets-n_inactive_packets
        write(27,102)  'number of inactive packets',n_inactive_packets
        write(27,102)  'number of escaped packets in line of sight',n_los_packets
        write(27,101)  '% of escaped packets in line of sight',real(n_los_packets)*100/real((n_init_packets-n_inactive_packets-n_abs_packets))
        write(27,102)  'number of absorbed packets',n_abs_packets
        write(27,101)  '% of absorbed packets', real(n_abs_packets)*100/real(n_init_packets-n_inactive_packets)
        write(27,101)  'absorbed weight %',abs_frac*100/real(n_packets-n_inactive_packets)
!        if (dust_geometry%lg_clumped) then
           do iG=1,mothergrid%tot_cells
              write(31,*) grid_cell(iG)%axis(1),grid_cell(iG)%axis(2),grid_cell(iG)%axis(3),grid_cell(iG)%n_rho
           end do
!        end if


        close(25)
        close(26)
        close(27)
        close(28)
        close(29)
        close(30)
        close(31)

    end subroutine


end module
