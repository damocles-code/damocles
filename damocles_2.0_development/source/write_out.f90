MODULE initialise

    USE globals
    USE class_line
    USE class_geometry
    USE class_dust
    USE class_grid
    USE class_freq_grid
    USE input


    IMPLICIT NONE

contains

    SUBROUTINE write_to_file()

        !CALCULATE OPACITIES
        !CALL calculate_opacities()

        !CONSTRUCT GRID WITH NUMBER OF PHOTONS EMITTED PER CELL
        !CALL construct_grid()
        !call build_dust_grid()

                    PRINT*,"active rest wavelength:",line%wavelength

        !CALCULATE FREQUENCY OF EMITTED LINE IN Hz
        PRINT*,'minimum velocity',(gas_geometry%R_min/gas_geometry%R_max)*gas_geometry%v_max

        PRINT*,'dust vel index',dust_geometry%v_power
        PRINT*,'gas vel index',gas_geometry%v_power
        PRINT*,'dust density index',dust_geometry%rho_power
        PRINT*,'gas density index',gas_geometry%rho_power
        PRINT*,'Rmax dust',dust_geometry%R_max
        PRINT*,'Rmin dust',dust_geometry%R_min
        PRINT*,'Rmax gas',gas_geometry%R_max
        PRINT*,'Rmin gas',gas_geometry%R_min
        PRINT*,'LAMBDA_0',line%wavelength

                PRINT*,'dust r ratio',dust_geometry%R_ratio
        PRINT*,'gas r ratio',gas_geometry%R_ratio

                        PRINT*,'total volume of supernova in e42cm^3',tot_vol
                                        PRINT*,'VOLUME OF GRID CELL (and therefore clump) in e42cm^3',mothergrid%cell_vol
                PRINT*,'GRID CELL WIDTHS in cm: ','X',mothergrid%cell_width(1),'Y',mothergrid%cell_width(2),'Z',mothergrid%cell_width(3)

        !CALCULATE AVERAGE OPTICAL DEPTH FROM Rin to Rout
        PRINT*,'EXTINCTION FOR LAMBDA_0 ',dust%lambda_ext
        PRINT*,'SCATTERING FOR LAMBDA_0',dust%lambda_sca
        PRINT*,'ALBEDO FOR LAMBDA_0',dust%lambda_sca/dust%lambda_ext
        PRINT*,'AVERAGE OPTICAL DEPTH IN LAMBDA_0',ndustav*dust%lambda_ext*(dust_geometry%R_max_cm-dust_geometry%R_min_cm)
        PRINT*,'AVERAGE OPTICAL DEPTH IN V',ndustav*dust%lambda_ext_V*(dust_geometry%R_max_cm-dust_geometry%R_min_cm)
        PRINT*,'AVERAGE NUMBER DENSITY',ndustav


        PRINT*,'mass weight',dust%species%mweight
        DO ii=1,dust%n_species
        PRINT*,'species number:',dust%species(ii)%id

            PRINT*,'albedo at rest frame wavelength',dust%species(ii)%albedo(line%wav_bin)
        END DO

        PRINT*,'EFFECTIVE SPHERICAL RADIUS OF CLUMP (cm)',mothergrid%cell_width(1)*(3.0/(4.0*pi))**0.3333333
        PRINT*,'EFFECTIVE SPHERICAL RADIUS AS A FRACTION OF Rout',(mothergrid%cell_width(1)/dust_geometry%R_max_cm)*(3.0/(4.0*pi))**0.3333333
        PRINT*,'AVERAGE OPTICAL DEPTH OF A CELL IN LAMBDA_0',ndustav*dust%lambda_ext*0.5*mothergrid%cell_width(1)*(3.0/(4.0*pi))**0.3333333
        PRINT*,'AVERAGE OPTICAL DEPTH OF A CELL IN LAMBDA_V',ndustav*dust%lambda_ext_V*0.5*mothergrid%cell_width(1)*(3.0/(4.0*pi))**0.3333333

        !call construct_freq_grid()

        !post_processing

        !post processing
        !!the below needs tidying and review
        IF (lg_doublet) THEN
            PRINT*,line%luminosity,n_init_packets,n_inactive_packets,2.0*n_init_packets-n_inactive_packets
            E_0=line%luminosity/real(2.0*n_init_packets-n_inactive_packets)!Energy of a single packet in W/um (note uses actual number of active photons)

        ELSE
            E_0=line%luminosity/real(n_init_packets-n_inactive_packets)
        END IF



        PRINT*,''
        PRINT*,'**********************'
        !PRINT*,'All percentages out of total number of active packets:'
        !PRINT*,''
        PRINT*,'TOTAL NUMBER OF PACKETS',n_init_packets
        PRINT*,'NUMBER OF ACTIVE(PROPAGATED) PACKETS',n_init_packets-n_inactive_packets
        !WRITE OUT ENERGY FILE - wavelength, velocity, energy (W/um)
        DO ii=1,nu_grid%n_bins-1
            lambda_bin=(c*10**9)*(0.5/nu_grid%bin(ii,1)+0.5/nu_grid%bin(ii+1,1))
            IF (.not. lg_doublet) THEN
                vel_bin=(c*1e-3*(lambda_bin**2-line%wavelength**2)/(line%wavelength**2+lambda_bin**2))
                !WRITE(13,*) lambda_bin,vel_bin,E_0*NP_BIN(inu)
                WRITE(15,*) lambda_bin,vel_bin,E_0*NP_BIN(ii)    !duplicate file for plotting
            ELSE
                 !WRITE(13,*) lambda_bin,E_0*NP_BIN(ii)

                WRITE(15,*) lambda_bin,E_0*NP_BIN(ii)    !duplicate file for plotting
            END IF
        END DO

        !PRINT*,'*******************'
        !PRINT*,'absorbed weight',w_abs
        !        IF (lg_doublet) THEN
        !            PRINT*,'undepleted luminosity (in units e40 erg/s)',line%luminosity/(1-w_abs/real(n-n_inactive))
        !        ELSE
        !            PRINT*,'undepleted luminosity (in units e40 erg/s)',line%luminosity/(1-w_abs/real(n-n_inactive))
        !        END IF
        !CLOSE(13)
        CLOSE(15)
        !CLOSE(17)
        PRINT*,'albedo',dust%lambda_sca/dust%lambda_ext
        PRINT*,'g_scat',dust%lambda_g_param
        IF (.not. lg_doublet) THEN
            PRINT*,'NUMBER OF INACTIVE PACKETS',n_inactive_packets
            PRINT*,'NUMBER OF ABSORBED PACKETS',n_abs_packets,real(n_abs_packets)*100/real(n_init_packets-n_inactive_packets),'%'
            PRINT*,'ABSORBED WEIGHT PERCENTAGE',abs_frac*100/real(n_packets-n_inactive_packets),'%'
            !PRINT*,''
            !PRINT*,'TOTAL ENERGY',line%luminosity
            !PRINT*,'ENERGY PER PACKET (ACTIVE ONLY)',E_0
            !PRINT*,'TOTAL ENERGY ABSORBED',n_abs_packets*E_0
            !PRINT*,''
            PRINT*,'FRACTION OF ESCAPED PACKETS IN LINE OF SIGHT',n_los_packets,real(n_los_packets)*100/real((n_init_packets-n_inactive_packets-n_abs_packets)),'%'
            !PRINT*,''
            PRINT*,'DUST MASS',dust%mass
            !PRINT*,'OUTPUT FILENAME:  ', trim(filename)
        ELSE
            PRINT*,'DUST MASS',dust%mass

            PRINT*,'percentage of absorbed packets',real(n_abs_packets*100.0)/real(n_init_packets-n_inactive_packets)
            PRINT*,'absorbed weight',abs_frac*100/real(n_init_packets-n_inactive_packets)
        END IF

        PRINT*,'number abs',n_abs_packets,'number initialised',n_init_packets,'number inactive',n_inactive_packets,'total active and escaped?',n_los_packets
        PRINT*,'AVERAGE OPTICAL DEPTH IN LAMBDA_0',ndustav*dust%lambda_ext*(dust_geometry%R_max_cm-dust_geometry%R_min_cm)
        PRINT*,'AVERAGE OPTICAL DEPTH IN V',ndustav*dust%lambda_ext_V*(dust_geometry%R_max_cm-dust_geometry%R_min_cm)



    END SUBROUTINE


END MODULE
