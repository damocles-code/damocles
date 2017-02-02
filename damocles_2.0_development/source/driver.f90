!-----------------------------------------------------------------------------------!
!  This is the main driver of the code DAMOCLES.  It is included as a module        !
!  in order to allow it to be run from other programs e.g. python wrappers.         !
!                                                                                   !
!  The run_damocles subroutine calls the subroutines that construct the grids,      !
!  emit and propagate packets through the grid and collates all escaped packets.    !
!  The model comparison module is also called from here.                            !
!-----------------------------------------------------------------------------------!

MODULE driver

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

    IMPLICIT NONE

contains

    SUBROUTINE run_damocles()

        !read input:
        call read_input()
        
        !construct grids and initialise simulation:
        DO iDoublet=1,2

            !generate random seed for random number generators (ensures random numbers change on each run)
            call init_random_seed

            !construct all grids and initialise rest line wavelength/freq
            IF (iDoublet==1) THEN
                !set active rest frame wavelength
                line%wavelength=line%doublet_wavelength_1
                line%frequency=c*10**9/line%wavelength

                !construct grids
                call calculate_opacities()
                call build_dust_grid()
                call construct_freq_grid()
                call build_emissivity_dist()
                call n_e_const()

                !initialise counters to zero
                n_init_packets=0
                n_inactive_packets=0
                n_abs_packets=0
                abs_frac=0
                n_los_packets=0

                !open output file to record resultant modelled line profile
                OPEN(15,file='output/output.out')

            ELSE IF (iDoublet==2) THEN
                !exit if not a doublet
                IF (.not. lg_doublet) EXIT

                !otherwise reset rest frame wavelength of line to be modelled
                !active wavelength is now second component of doublet
                line%wavelength=line%doublet_wavelength_2
                line%frequency=c*10**9/line%wavelength
            END IF

            !emit and propagate packets through grid
            PRINT*,"Propagating packets..."

            !entire simulation run for each component of doublet (if applicable)
            !absorbed weight stored for complete doublet (i.e. both components)
            !initialise absorbed weight of packets to zero
            abs_frac=0

            SELECT CASE(gas_geometry%type)

                CASE('shell')
                    !if all emission from clumps within shell structure
                    IF (gas_geometry%clumped_mass_frac == 1) THEN
                        DO id_no=1,mothergrid%tot_cells
                            IF (grid_cell(id_no)%cellStatus == 1) THEN
                                !equal number of packets to be emitted in each clump
                                NP(id_no)=n_packets/ncl
                                n_init_packets=n_init_packets+NP(id_no)
                                call run_packets()
                            END IF
                        END DO
                    !else if all emission from smooth shell
                    ELSE
                        DO id_no=1,n_shells
                            n_init_packets=n_init_packets+NP(id_no)
                            call run_packets()
                        END DO
                    END IF

                CASE('arbitrary')
                    !emission per cell scaled with dust mass from specified dust grid
                    DO id_no=1,mothergrid%tot_cells
                        !n is cumulative number of packets run through grid (check number)
                        n_init_packets=n_init_packets+NP(id_no)
                        call run_packets()
                    END DO

                CASE DEFAULT
                    PRINT*,'You have not selected a shell or arbitrary distribution.  Alternative distributions have not yet been included.'
                    PRINT*,'Please construct a grid using the gridmaker at http://www.nebulousresearch.org/codes/mocassin/mocassin_gridmaker.php and use the arbitrary option.  Aborted.'
                    STOP

            END SELECT
        END DO

        !calculate goodness of fit to data if supplied
        IF (lg_data) call linear_interp(chi2)

        !write out log file
        call write_to_file()

        !decallocate all allocated memory
        DEALLOCATE(grid_cell)
        DEALLOCATE(nu_grid%bin)
        DEALLOCATE(tmp)
        DEALLOCATE(mothergrid%x_div)
        DEALLOCATE(mothergrid%y_div)
        DEALLOCATE(mothergrid%z_div)
        DEALLOCATE(np)
        DEALLOCATE(np_bin)
        DEALLOCATE(dust%species)
        IF (dust_geometry%type == "shell") DEALLOCATE(RSh)

        PRINT*,'Complete!'

    END SUBROUTINE run_damocles

    SUBROUTINE run_packets()

        !!!!!OPENMP HAS NOT BEEN UPDATED AFTER RECENT AMENDMENTS SO DO NOT EMPLOY WITHOUT THOROUGH CHECKING FIRST!!!!!!!!!!!
        !!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(dust_geometry%ff,gas_geometry%rho_power,line%doublet_ratio,iDoublet,n_los_packets,n_inactive,n_abs_packets,shell_width,width,NP_BIN,NP,iSh,RSh,dust_geometry%R_min,dust_geometry%R_max,gas_geometry%R_min,gas_geometry%R_max,lg_LOS,lg_ES,grid,mothergrid%n_cells,nu_grid%bin,dust_geometry%v_max,gas_geometry%v_max,l,gas_geometry%v_power,dummy,mgrid,dust,lg_vel_shift)


        !!$OMP DO SCHEDULE(dynamic)
        DO iP=1,NP(id_no)

            call emit_packet()

            IF (packet%lg_active) THEN

                !propagate active packet through grid
                call propagate()

                !if packet has been absorbed then record
                IF (packet%lg_abs) THEN
                    !!$OMP CRITICAL
                    n_abs_packets=n_abs_packets+1
                    IF (iDoublet==2) THEN
                        abs_frac=abs_frac+packet%weight/line%doublet_ratio
                    ELSE
                        abs_frac=abs_frac+packet%weight
                    END IF
                    !!$OMP END CRITICAL
                ELSE
                !if the packet has not been absorbed then record in resultant profile

                    !if taking integrated profile and not interested in line of sight, record all escaped packets
                    IF (.not. lg_los) THEN
                         call add_packet_to_profile()
                    ELSE
                        !Only add active packets to profile for those in LoS
                        IF (packet%lg_los) call add_packet_to_profile()

                    END IF !line of sight
                END IF  !absorbed/escaped
            END IF  !active
        END DO
        !!$OMP END DO

        !!$OMP END PARALLEL

        !close log file (opened in read_input)
        CLOSE(55)

    END SUBROUTINE

    SUBROUTINE add_packet_to_profile()
        !Find the smallest distance and thus nearest freq point
        packet%freq_id=MINLOC(packet%nu-nu_grid%bin(:,1),1,(packet%nu-nu_grid%bin(:,1))>0)

        IF (packet%freq_id==0) THEN
            PRINT*,'photon outside frequency range',packet%freq_id,packet%nu,packet%weight
        ELSE
            IF (iDoublet==2) THEN
                packet%weight=packet%weight/line%doublet_ratio
            END IF
            !!$OMP CRITICAL
            NP_BIN(packet%freq_id)=NP_BIN(packet%freq_id)+packet%weight
            n_los_packets=n_los_packets+1
            !!$OMP END CRITICAL
        END IF
    END SUBROUTINE

END MODULE driver
