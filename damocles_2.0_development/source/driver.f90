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

    INTEGER,EXTERNAL ::  omp_get_num_threads, omp_get_thread_num

    INTEGER ::  freqid,los


    REAL    ::  theta,w,w_abs

    INTEGER ::  scatno
    INTEGER ::  NP_clump

contains

    SUBROUTINE run_code()

        implicit none
        
        REAL  :: chi2

        !read variables in input files
        call read_input()

        !entire simulation run for each component of doublet (if applicable)
        !absorbed weight stored for complete doublet (i.e. both components)
        !initialise absorbed weight of packets to zero
        w_abs=0
        
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

                !!initialise counters to zero
                !!what are these?
                n=0
                n_inactive=0
                nabs=0
                w_abs=0
                nabs=0

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

            PRINT*,"active rest wavelength:",line%wavelength

            !!what is this variable?  total number of packets for single component?
            tot=0

            !emit and propagate packets through grid
            PRINT*,"Propagating packets..."
            SELECT CASE(gas_geometry%type)

                CASE('shell')

                    !if all emission from clumps within shell structure
                    IF (gas_geometry%clumped_mass_frac == 1) THEN
                        DO ii=1, mothergrid%tot_cells
                            IF (grid_cell(ii)%cellStatus == 1) THEN
                                NP(ii)=n_packets/ncl
                                !n is cumulative number of packets run through grid (check number)
                                n=n+NP(ii)
                                unit_vol_iD=ii
                                call run_packets()
                            END IF
                        END DO
                    !else if all emission from smooth shell
                    ELSE
                        DO ii=1,n_shells
                            !n is cumulative number of packets run through grid (check number)
                            n=n+NP(ii)
                            !packet%cell_no is the cell in which the packet is located and is identified after emission in the rountine 'run_packets'
                            packet%cell_no=0
                            unit_vol_iD=ii
                            call run_packets()
                        END DO
                    END IF

                CASE('arbitrary')
                    !emission per cell scaled with dust mass from specified dust grid
                    DO ii=1,mothergrid%tot_cells
                        !n is cumulative number of packets run through grid (check number)
                        n=n+NP(ii)
                        unit_vol_iD=ii
                        call run_packets()
                    END DO

                CASE DEFAULT
                    PRINT*,'You have not selected a gas or arbitrary distribution.  Alternative distributions have not yet been included.'
                    PRINT*,'Please construct a grid using the gridmaker at http://www.nebulousresearch.org/codes/mocassin/mocassin_gridmaker.php and use the arbitrary option.  Aborted.'
                    STOP

            END SELECT


        END DO

        !post processing
        !!the below needs tidying and review
        IF (lg_doublet) THEN
            PRINT*,line%luminosity,n,n_inactive,2.0*n-n_inactive
            E_0=line%luminosity/real(2.0*n-n_inactive)!Energy of a single packet in W/um (note uses actual number of active photons)
        
        ELSE
            E_0=line%luminosity/real(n-n_inactive)
        END IF



        PRINT*,''
        PRINT*,'**********************'
        !PRINT*,'All percentages out of total number of active packets:'
        !PRINT*,''
        PRINT*,'TOTAL NUMBER OF PACKETS',n
        PRINT*,'NUMBER OF ACTIVE(PROPAGATED) PACKETS',n-n_inactive
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
        IF (.not. lg_doublet) THEN
            PRINT*,'NUMBER OF INACTIVE PACKETS',n_inactive
            PRINT*,'NUMBER OF ABSORBED PACKETS',nabs,real(nabs)*100/real(n-n_inactive),'%'
            PRINT*,'ABSORBED WEIGHT PERCENTAGE',w_abs*100/real(n-n_inactive),'%'
            !PRINT*,''
            !PRINT*,'TOTAL ENERGY',line%luminosity
            !PRINT*,'ENERGY PER PACKET (ACTIVE ONLY)',E_0
            !PRINT*,'TOTAL ENERGY ABSORBED',nabs*E_0
            !PRINT*,''
            PRINT*,'FRACTION OF ESCAPED PACKETS IN LINE OF SIGHT',tot,real(tot)*100/real((n-n_inactive-nabs)),'%'
            !PRINT*,''
            PRINT*,'DUST MASS',dust%mass
            !PRINT*,'OUTPUT FILENAME:  ', trim(filename)
        ELSE
            PRINT*,'DUST MASS',dust%mass
            PRINT*,nabs,n,n_inactive
            PRINT*,'percentage of absorbed packets',real(nabs*100.0)/real(n-n_inactive)
            PRINT*,'absorbed weight',w_abs*100/real(n-n_inactive)
        END IF

        PRINT*,'AVERAGE OPTICAL DEPTH IN LAMBDA_0',ndustav*dust%lambda_ext*(dust_geometry%R_max_cm-dust_geometry%R_min_cm)
        PRINT*,'AVERAGE OPTICAL DEPTH IN V',ndustav*dust%lambda_ext_V*(dust_geometry%R_max_cm-dust_geometry%R_min_cm)

        !calculate goodness of fit to data if supplied
        IF (lg_data) call linear_interp(chi2)

        !decallocate all memore
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

    END SUBROUTINE run_code

    SUBROUTINE run_packets()

        !!!!!OPENMP HAS NOT BEEN UPDATED AFTER RECENT AMENDMENTS SO DO NOT EMPLOY WITHOUT THOROUGH CHECKING FIRST!!!!!!!!!!!
        !!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(dust_geometry%ff,gas_geometry%rho_power,line%doublet_ratio,iDoublet,tot,n_inactive,nabs,shell_width,width,NP_BIN,NP,iSh,RSh,dust_geometry%R_min,dust_geometry%R_max,gas_geometry%R_min,gas_geometry%R_max,lg_LOS,lg_ES,grid,mothergrid%n_cells,nu_grid%bin,dust_geometry%v_max,gas_geometry%v_max,l,gas_geometry%v_power,dummy,mgrid,dust,lg_vel_shift)
        !PRINT*,'num of threads', omp_get_num_threads()

        !!$OMP DO SCHEDULE(dynamic)
        DO iP=1,NP(unit_vol_iD)
            !PRINT*,'thread number',omp_get_thread_num()

            call emit_packet()

            !PRINT*,'thread no',omp_get_thread_num()
            IF (packet%lg_active) THEN
                scatno=0
                packet%lg_abs=.false.

                call propagate()
                    
                theta=acos(packet%pos_cart(3)/((packet%pos_cart(1)**2+packet%pos_cart(2)**2+packet%pos_cart(3)**2)**0.5))
                IF (packet%lg_abs) THEN
                    !!$OMP CRITICAL
                    nabs=nabs+1
                    IF (iDoublet==2) THEN
                        w_abs=w_abs+packet%weight/line%doublet_ratio
                    ELSE
                        w_abs=w_abs+packet%weight
                    END IF
                    !!$OMP END CRITICAL
                END IF
                            
                ! IF PHOTON IN LINE OF SIGHT THEN CALCULATE WHICH FREQ BIN PACKET IN AND ADD 1 TO TOTAL IN FREQ BIN
                IF (.not. lg_LOS) THEN
                    IF (.not. packet%lg_abs) THEN
                                   
                        los=1.
                        tmp(:,1)=(packet%nu-nu_grid%bin(:,1))                         !Calculate distance to each freq point
                        freq=MINLOC(tmp,1,tmp>0)                                !Find the smallest distance and thus nearest freq point
                        freqid=freq(1)                              !Attach id of freq bin to photon id
                        IF (freqid==0) THEN
                            PRINT*,'photon outside frequency range',freqid,packet%nu,packet%weight
                            PRINT*,nu_grid%fmin,nu_grid%fmax
                        ELSE
                            IF (iDoublet==2) THEN
                                packet%weight=packet%weight/line%doublet_ratio
                                        
                            END IF
                             !!$OMP CRITICAL
                            !dummy=NP_BIN(freqid)+w                  !Add 1 to number of photons in that freq bin

                            !NP_BIN(freqid)=dummy

                            NP_BIN(freqid)=NP_BIN(freqid)+packet%weight

                            tot=tot+1
                             !!$OMP END CRITICAL
                        END IF
                    END IF
                ELSE
                    IF ((theta<(pi/6)) .AND. (.not. packet%lg_abs)) THEN                       !Only calculate freq bin for those in LoS
                        !Packets which are not absorbed
                        los=1.
                        tmp(:,1)=(packet%nu-nu_grid%bin(:,1))                     !Calculate distance to each freq point
                        freq=MINLOC(tmp,1,tmp>0)                            !Find the smallest distance and thus nearest freq point
                        freqid=freq(1)                          !Attach id of freq bin to photon id
                        IF (freqid==0) THEN
                            PRINT*,'photon outside frequency range',freqid,packet%nu,w
                        ELSE
                             !!$OMP CRITICAL
                            IF (lg_doublet) THEN
                                IF (iDoublet==2) THEN
                                    packet%weight=packet%weight/line%doublet_ratio
                                END IF
                            END IF
                            dummy=NP_BIN(freqid)+packet%weight                  !Add 1 to number of photons in that freq bin
                            !PRINT*,packet%weight
                            tot=tot+1
                             !!$OMP END CRITICAL
                        END IF
                    END IF
                END IF
            END IF
        !                        PRINT*,'threadno',omp_get_thread_num()
        END DO

                    
         !!$OMP END DO

         !!$OMP END PARALLEL

        !close log file (opened in read_input)
        CLOSE(55)
    END SUBROUTINE

END MODULE driver
