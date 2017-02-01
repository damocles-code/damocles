!-------------------------------------------------------------------------------!
! this is the recursive subroutine (and smaller subruotines) that are           !
! responsible for the transport of a packet through the grid.                   !
!      -absorption and scattering opacities are calculated for the current cell !
!      -optical depths and the distance from the packet to the nearest cell     !
!       wall are calculated                                                     !
!      -the packet either continues, is scattered by dust, is scattered by      !
!       electrons or is absorbed                                                !
!      -if the packet is not absorbed and does not escape, the routine is       !
!       called again                                                            !
!-------------------------------------------------------------------------------!

MODULE radiative_transfer

    use globals
    use class_dust
    use class_grid
    use input
    use initialise
    use vector_functions
    use random_routines
    use electron_scattering
    use class_packet

    implicit none

    REAL    ::  tau                    !optical depth sampled from cumulative frequency dsitribution at wavelength of active packet
    REAL    ::  kappa_rho              !opacity * mass density = C_ext (cross-section of interaction) * number density at at wavelength of active packet
    REAL    ::  C_ext_tot              !total extinction cross-section of interaction at wavelength of active packet
    REAL    ::  C_sca_tot              !total scattering cross-section of interaction at wavelength of active packet
    REAL    ::  albedo                 !albedo in cell at wavelength of active packet
    REAL    ::  v_therm(3)             !sampled velocity of electron based on specified electron temperature
    REAL    ::  s_face(3)              !distance to each (x/y/z) cell boundary in direction of travel
    REAL    ::  s_min                  !distance to nearest face in packet's direction of travel
    REAL    ::  s                      !distance travelled by packet based on sampled optical depth

    INTEGER ::  i_min                  !index of nearest face in packet's direction of travel
    INTEGER ::  wav_id                 !wavelength bin that contains the current packet

    CHARACTER(9) :: event_type

    !INTEGER ::  omp_get_thread_num

contains

    RECURSIVE SUBROUTINE propagate()

        call check_step_no()
        call update_cell_no()
        call calculate_extinction()

    !calculate distance to nearest face:

        !unit vector direction of travel of packet in cartesian
        packet%dir_cart=normalise(packet%dir_cart)
        !calculate distance to nearest face
        DO i_dir=1,3
            IF (packet%dir_cart(i_dir)<0) THEN
                s_face(i_dir)=ABS((grid_cell(packet%cell_no)%axis(i_dir)-packet%pos_cart(i_dir))/packet%dir_cart(i_dir))
            ELSE
                s_face(i_dir)=ABS((grid_cell(packet%cell_no)%axis(i_dir)+grid_cell(packet%cell_no)%width(i_dir)-packet%pos_cart(i_dir))/packet%dir_cart(i_dir))
            END IF
        END DO

        !index of nearest face (i.e. identifies whether nearest face is planar in x or y or z) and distance
        s_min=MINVAL(s_face)
        i_min=MINLOC(s_face,1)


    !sample optical depth and calculate distance travelled:

        !call random number and sample from CFD to obtain tau
        call random_number(ran)
        tau=-(ALOG((1-ran)))
        !work out potential distance travelled by packet based on optical depth tau
        !(distance in units of cm, since kappa_rho in units cm^-1)
        IF (grid_cell(packet%cell_no)%nrho>0) THEN
            IF (lg_ES) THEN
                s=tau/(kappa_rho+sigma_T*grid_cell(packet%cell_no)%N_e)
            ELSE
                s=tau/kappa_rho
            END IF
        ELSE
            !if cell has zero density then potential distance travelled is greatest distance across cell (for no e- scat)
            IF (lg_ES) THEN
                s=tau/(sigma_T*grid_cell(packet%cell_no)%N_e)
            ELSE
                s=(grid_cell(packet%cell_no)%width(1)**2+grid_cell(packet%cell_no)%width(2)**2+grid_cell(packet%cell_no)%width(3)**2)**0.5
            END IF
        END IF


    !event occurs when distance travelled (as determined by tau) is < distance to nearest face
    !else continues to cell boundary with no event occurring:

        IF ((s>s_min)) THEN
        !packet travels to cell boundary (direction of travel remains the same):

            !position updated to be on boundary with next cell
            !actually moves just past boundary by small factor...
            packet%pos_cart(:)=packet%pos_cart(:)+(ABS(s_min)+ABS(s_min)*1E-10)*packet%dir_cart(:)

            IF (packet%dir_cart(i_min)>0) THEN
                !if packet travels forwards then advance cell id by 1 in that index
                IF (packet%axis_no(i_min) /= mothergrid%n_cells(i_min)) THEN
                    packet%axis_no(i_min)=packet%axis_no(i_min)+1
                ELSE
                    !reached edge of grid, escapes
                    RETURN
                END IF
                !update id of cell containing packet and update position of packet
                call update_cell_no()
                packet%pos_cart(i_min)=grid_cell(packet%cell_no)%axis(i_min)+((ABS(s_min)*1E-10)*packet%dir_cart(i_min))
            ELSE
                !if packet travels backwards then reduce cell id by 1 in that index
                IF (packet%axis_no(i_min) /= 1) THEN
                    packet%axis_no(i_min)=packet%axis_no(i_min)-1
                ELSE
                    !reached edge of grid, escapes
                    RETURN
                END IF
                !update id of cell containing packet and update position of packet
                call update_cell_no()
                packet%pos_cart(i_min)=grid_cell(packet%cell_no)%axis(i_min)+((ABS(s_min)*1E-10)*packet%dir_cart(i_min))+grid_cell(packet%cell_no)%width(i_min)
            END IF

            !calculate packet radial position
            packet%r=(packet%pos_cart(1)**2+packet%pos_cart(2)**2+packet%pos_cart(3)**2)**0.5

            call check_packet_in_right_cell()
            call check_escaped()
            call propagate()

        ELSE
        !event does occur:

            !calculate position and radius of event
            packet%pos_cart(:)=packet%pos_cart(:)+s*packet%dir_cart(:)
            packet%r=(packet%pos_cart(1)**2+packet%pos_cart(2)**2+packet%pos_cart(3)**2)**0.5

            call calculate_velocity()
            call check_escaped()
            call determine_event_type()


            SELECT CASE(event_type)

                CASE("dust_scat")

                    IF (lg_vel_shift) THEN
                        !Inverse lorentz boost (observer frame to particle), sample new scat direction, lorentz boost (particle frame to observer)
                        call inv_lorentz_trans(packet%vel_vect,packet%dir_cart,packet%nu,packet%weight,"scat")
                        call scatter()
                        call lorentz_trans(packet%vel_vect,packet%dir_cart,packet%nu,packet%weight,"scat")
                    END IF
                    call propagate()

                CASE("dust_absn")

                    packet%lg_abs=.true.
                    RETURN

                CASE("elec_scat")

                    !calculate thermal velocity and random thermal velocity vector
                    v_therm(1)=normal(0.d0,dble(maxwell_sigma))
                    v_therm(2)=normal(0.d0,dble(maxwell_sigma))
                    v_therm(3)=normal(0.d0,dble(maxwell_sigma))
                    packet%vel_vect=packet%vel_vect+v_therm
                    !!weighting for e- scattering
                    IF (lg_vel_shift) THEN
                        !Inverse lorentz boost (observer frame to particle), sample new scat direction, lorentz boost (particle frame to observer)
                        call inv_lorentz_trans(packet%vel_vect,packet%dir_cart,packet%nu,packet%weight,"escat")
                        call scatter()
                        call lorentz_trans(packet%vel_vect,packet%dir_cart,packet%nu,packet%weight,"escat")
                    END IF
                    call propagate()

            END SELECT

        END IF

    END SUBROUTINE propagate

    !---------------------------------------------------------------------

    SUBROUTINE calculate_extinction()
        !calculate extinction and scattering for this cell at this wavelength

        REAL    ::  C_ext(dust%n_species)
        REAL    ::  C_sca(dust%n_species)

        !calculate extinction opacity and albedo
        DO i_spec=1,dust%n_species

            !Calculate difference between actual wavelength...
            !...and wavelength bins in order to calculate which bin photon is in for each species
            wav_id=MINLOC(ABS((dust%species(i_spec)%wav(:)-(c*1e6/packet%nu))),1)

            !calculate opactiy as function of rho for specific wavelength by interpolating between bins for each species
            !!interpolate function?
            IF ((c*1e6/packet%nu-dust%species(i_spec)%wav(wav_id))<0) THEN
                C_ext(i_spec)=dust%species(i_spec)%C_ext(wav_id)-((dust%species(i_spec)%C_ext(wav_id)-dust%species(i_spec)%C_ext(wav_id-1))* &
                    & ((dust%species(i_spec)%wav(wav_id)-c*1e6/packet%nu)/(dust%species(i_spec)%wav(wav_id)-dust%species(i_spec)%wav(wav_id-1))))
            ELSE
                C_ext(i_spec)=dust%species(i_spec)%C_ext(wav_id)+((dust%species(i_spec)%C_ext(wav_id+1)-dust%species(i_spec)%C_ext(wav_id))* &
                    & ((c*1e6/packet%nu-dust%species(i_spec)%wav(wav_id))/(dust%species(i_spec)%wav(wav_id+1)-dust%species(i_spec)%wav(wav_id))))
            END IF

            !cumulative scattering component of extinction for albedo calculation
            IF ((c*1e6/packet%nu-dust%species(i_spec)%wav(wav_id))<0) THEN
                C_sca(i_spec)=dust%species(i_spec)%C_sca(wav_id)-((dust%species(i_spec)%C_sca(wav_id)-dust%species(i_spec)%C_sca(wav_id-1))* &
                    & ((dust%species(i_spec)%wav(wav_id)-c*1e6/packet%nu)/(dust%species(i_spec)%wav(wav_id)-dust%species(i_spec)%wav(wav_id-1))))
            ELSE
                C_sca(i_spec)=dust%species(i_spec)%C_sca(wav_id)+((dust%species(i_spec)%C_sca(wav_id+1)-dust%species(i_spec)%C_sca(wav_id))* &
                    & ((c*1e6/packet%nu-dust%species(i_spec)%wav(wav_id))/(dust%species(i_spec)%wav(wav_id+1)-dust%species(i_spec)%wav(wav_id))))
            END IF
        END DO

        !calculate total opactiies weighted over all species
        C_ext_tot=sum(C_ext*dust%species%weight)
        C_sca_tot=sum(C_sca*dust%species%weight)

        !calculate albedo (don't add weighted albedos, must add each component and then divide total scat by total ext)
        albedo=C_sca_tot/C_ext_tot

        !Calculate overall opactiy using rho
        !kappa_rho = kappa*rho = n*Cext (units of cm^-1)
        kappa_rho=(C_ext_tot*grid_cell(packet%cell_no)%nrho)

    END SUBROUTINE

    !---------------------------------------------------------------------

    SUBROUTINE calculate_velocity()
        !calculate bulk velocity of scattering e- and velocity unit vector
        packet%v=dust_geometry%v_max*((packet%r/(dust_geometry%R_max*1e15))**dust_geometry%v_power)
        packet%vel_vect=normalise(packet%pos_cart)*packet%v
    END SUBROUTINE

    !---------------------------------------------------------------------

    SUBROUTINE scatter()
        !Now scatter (sample new direction)
        call random_number(random)
        packet%dir_sph(:)=(/ ((2*random(1))-1),random(2)*2*pi /)
        packet%dir_cart(:)=cart(ACOS(packet%dir_sph(1)),packet%dir_sph(2))
    END SUBROUTINE

    !---------------------------------------------------------------------

    SUBROUTINE update_cell_no()
        !calculate overall id of cell using x,y and z ids
        !note that in list of all cells, cells listed changing first z, then y, then x e.g.
        ! 1 1 1, 1 1 2, 1 1 3, 1 2 1, 1 2 2, 1 2 3, 2 1 1, 2 1 2... etc.
        packet%cell_no=(mothergrid%n_cells(2)*mothergrid%n_cells(3)*(packet%axis_no(1)-1)+mothergrid%n_cells(3)*(packet%axis_no(2)-1)+packet%axis_no(3))
    END SUBROUTINE

    !---------------------------------------------------------------------

    SUBROUTINE check_packet_in_right_cell()
        DO i_dir=1,3
            IF (packet%pos_cart(i_dir)<grid_cell(packet%cell_no)%axis(i_dir)) THEN
                !idGP(i_dir)=idGP(i_dir)-1
                !iGPP=(mothergrid%n_cells(2)*mothergrid%n_cells(3)*(idGP(1)-1))+mothergrid%n_cells(3)*(idGP(2)-1)+idGP(3)
                PRINT*,'Error - packet coordinates are not in the identified cell. Packet removed.'
                !PRINT*,'1',i_dir,grid_cell(packet%cell_no)%axis(i_dir),packet%pos_cart(i_dir),grid_cell(packet%cell_no)%axis(i_dir)+grid_cell%width(i_dir)
                RETURN
            ELSE IF (packet%pos_cart(i_dir)>(grid_cell(packet%cell_no)%axis(i_dir)+grid_cell(packet%cell_no)%width(i_dir))) THEN
                !idGP(i_dir)=idGP(i_dir)+1
                !iGPP=(mothergrid%n_cells(2)*mothergrid%n_cells(3)*(idGP(1)-1))+mothergrid%n_cells(3)*(idGP(2)-1)+idGP(3)
                PRINT*,'Error - packet coordinates are not in the identified cell. Packet removed.'
                !PRINT*,'2',packet%axis_no(i_min),i_dir,i_min,packet%pos_cart(i_dir),grid_cell(packet%cell_no)%axis(i_dir)+grid_cell(packet%cell_no)%width(i_dir),iP!,grid_cell(packet%cell_no)%axis(i_dir),packet%pos_cart(i_dir),grid_cell(packet%cell_no)%axis(i_dir)+grid_cell%width(i_dir)
                RETURN
            END IF
        END DO
    END SUBROUTINE

    !---------------------------------------------------------------------

    SUBROUTINE check_escaped()
        IF (packet%axis_no(i_min) > mothergrid%n_cells(1) .OR. packet%axis_no(i_min)<1 .OR. packet%r>(MAX(gas_geometry%R_max,dust_geometry%R_max)*1E15)) THEN
            RETURN
        END IF
    END SUBROUTINE

    !---------------------------------------------------------------------

    SUBROUTINE check_step_no()
        !!check number of steps of this packet - a packet could get stuck in highly scattering environments
        packet%step_no=packet%step_no+1
        IF (packet%step_no>500) THEN
            packet%lg_active=.false.
            PRINT*, packet%step_no
            RETURN
        END IF
    END SUBROUTINE

    !---------------------------------------------------------------------

    SUBROUTINE determine_event_type()
        call random_number(ran)

        !if ES used then establish whether dust event or e- scattering event
        IF ((.not. lg_ES) .OR. (ran<kappa_rho/(kappa_rho+sigma_T*grid_cell(packet%cell_no)%N_e))) THEN
            !dust event - either scattering or absorption...

            !generate random number to compare to dust albedo in order to determine whether dust absorption or scattering
            call random_number(ran)

            IF (ran<albedo) THEN
                !dust scattering event
                event_type='dust_scat'
            ELSE
                event_type='dust_absn'
            END IF
        ELSE
            event_type='elec_scat'
        END IF
    END SUBROUTINE

    !---------------------------------------------------------------------

END MODULE
