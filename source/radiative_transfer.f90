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

module radiative_transfer

    use globals
    use class_dust
    use class_grid
    use input
    use initialise
    use vector_functions
    use random_routines
    use electron_scattering
    use class_packet
    use rnglib

    implicit none

    integer ::  i_min                  !index of nearest face in packet's direction of travel
    integer ::  wav_id                 !wavelength bin that contains the current packet
    integer ::  i_spec                 !loop index to loop over different species when calculating extinction
    real    ::  nu_ext                 !value of frequency of packet in frame of particle (used to calculate extinction, not necessarily the packet frequency)
    character(9) :: event_type         !is packet experiences an event, this describes whether it is
                                       !dust scattering, electron scattering or absorption by dust

    !----variables below are properties of current active packet at its current wavelength
    real    ::  kappa_rho              !opacity * mass density = c_ext (cross-section of interaction) * number density
    real    ::  c_ext_tot              !total extinction cross-section of interaction
    real    ::  c_sca_tot              !total scattering cross-section of interaction
    real    ::  g_param_tot            !total forward scattering parameter g (from mie calculation, to be used with henyey-greenstein approx)
    real    ::  albedo                 !albedo in cell
    !$OMP THREADPRIVATE(kappa_rho,c_ext_tot,c_sca_tot,g_param_tot,albedo,i_min,wav_id,event_type,nu_ext,i_spec)

contains

    recursive subroutine propagate()
    use rnglib

        implicit none

        real    ::  v_therm(3)             !sampled velocity of electron based on specified electron temperature
        real    ::  s_face(3)              !distance to each (x/y/z) cell boundary in direction of travel
        real    ::  s_min                  !distance to nearest face in packet's direction of travel
        real    ::  s                      !distance travelled by packet based on sampled optical depth

            !----variables below are properties of current active packet at its current wavelength
            real    ::  tau                    !optical depth sampled from cumulative frequency dsitribution


        call check_step_no()
        call update_cell_no()
        call calculate_extinction()

        !calculate distance to nearest face:

        !unit vector direction of travel of packet in cartesian
        packet%dir_cart=normalise(packet%dir_cart)
        !calculate distance to nearest face
        do i_dir=1,3
            if (packet%dir_cart(i_dir)<0) then
                s_face(i_dir)=abs((grid_cell(packet%cell_no)%axis(i_dir)-packet%pos_cart(i_dir))/packet%dir_cart(i_dir))
            else
                s_face(i_dir)=abs((grid_cell(packet%cell_no)%axis(i_dir)+grid_cell(packet%cell_no)%width(i_dir)-packet%pos_cart(i_dir))/packet%dir_cart(i_dir))
            end if
        end do

        !index of nearest face (i.e. identifies whether nearest face is planar in x or y or z) and distance
        s_min=minval(s_face)
        i_min=minloc(s_face,1)

        !sample optical depth and calculate distance travelled:

        !call random number and sample from cfd to obtain tau
!        call random_number(ran)
        ran = r4_uni_01()

        tau=-(alog((1-ran)))
        !work out potential distance travelled by packet based on optical depth tau
        !(distance in units of cm, since kappa_rho in units cm^-1)
        if (grid_cell(packet%cell_no)%n_rho>0) then
            if (lg_es) then
                s=tau/(kappa_rho+sigma_t*grid_cell(packet%cell_no)%n_e)
            else
                s=tau/kappa_rho
            end if
        else
            !if cell has zero density then potential distance travelled is greatest distance across cell (for no e- scat)
            if (lg_es) then
                s=tau/(sigma_t*grid_cell(packet%cell_no)%n_e)
            else
                s=(sum(grid_cell(packet%cell_no)%width**2))**0.5
            end if
        end if


        !event occurs when distance travelled (as determined by tau) is < distance to nearest face
        !else continues to cell boundary with no event occurring:

        if ((s>s_min)) then
            !packet travels to cell boundary (direction of travel remains the same):

            !position updated to be on boundary with next cell
            !actually moves just past boundary by small factor...
            packet%pos_cart(:)=packet%pos_cart(:)+(abs(s_min)+abs(s_min)*1e-10)*packet%dir_cart(:)

            if (packet%dir_cart(i_min)>0) then
                !if packet travels forwards then advance cell id by 1 in that index
                if (packet%axis_no(i_min) /= mothergrid%n_cells(i_min)) then
                    packet%axis_no(i_min)=packet%axis_no(i_min)+1
                else
                    !reached edge of grid, escapes
                    call check_los()
                    return
                end if
                !update id of cell containing packet and update position of packet
                call update_cell_no()
                packet%pos_cart(i_min)=grid_cell(packet%cell_no)%axis(i_min)+((abs(s_min)*1e-10)*packet%dir_cart(i_min))
            else
                !if packet travels backwards then reduce cell id by 1 in that index
                if (packet%axis_no(i_min) /= 1) then
                    packet%axis_no(i_min)=packet%axis_no(i_min)-1
                else
                    !reached edge of grid, escapes
                    call check_los()
                    return
                end if
                !update id of cell containing packet and update position of packet
                call update_cell_no()
                packet%pos_cart(i_min)=grid_cell(packet%cell_no)%axis(i_min)+((abs(s_min)*1e-10)*packet%dir_cart(i_min))+grid_cell(packet%cell_no)%width(i_min)
            end if

            !calculate packet radial position
            packet%r=(sum(packet%pos_cart**2))**0.5

            call check_packet_in_right_cell()
            call check_escaped()
            call propagate()

        else
            !event does occur:

            !calculate position and radius of event
            packet%pos_cart(:)=packet%pos_cart(:)+s*packet%dir_cart(:)
            packet%r=(sum(packet%pos_cart**2))**0.5

            call calculate_velocity()
            call check_escaped()
            call determine_event_type()

            select case(event_type)

                case("dust_scat")
                    if (lg_vel_shift) then
                       !inverse lorentz boost (observer frame to particle), sample new scat direction, lorentz boost (particle frame to observer)
                       call inv_lorentz_trans(packet%vel_vect,packet%dir_cart,packet%nu,packet%weight,"scat")
                       call scatter()
                       call lorentz_trans(packet%vel_vect,packet%dir_cart,packet%nu,packet%weight,"scat")
                    end if

                    call propagate()

                case("dust_absn")
                    packet%lg_abs=.true.
                    call check_los()
                    return

                case("elec_scat")

                    !calculate thermal velocity and random thermal velocity vector
                    v_therm(1)=normal(0.d0,dble(maxwell_sigma))
                    v_therm(2)=normal(0.d0,dble(maxwell_sigma))
                    v_therm(3)=normal(0.d0,dble(maxwell_sigma))
                    packet%vel_vect=packet%vel_vect+v_therm


                    !!weighting for e- scattering
                    if (lg_vel_shift) then
                        !inverse lorentz boost (observer frame to particle), sample new scat direction, lorentz boost (particle frame to observer)
                        call inv_lorentz_trans(packet%vel_vect,packet%dir_cart,packet%nu,packet%weight,"escat")
                        call scatter()
                        call lorentz_trans(packet%vel_vect,packet%dir_cart,packet%nu,packet%weight,"escat")
                    end if

                    call propagate()

            end select

        end if

    end subroutine propagate

    !---------------------------------------------------------------------

    subroutine calculate_extinction()
        !calculate extinction and scattering for this cell at this wavelength

        real    ::  c_ext(dust%n_species)
        real    ::  c_sca(dust%n_species)
        real    ::  g_param(dust%n_species)

        !calculate extinction opacity and albedo
        do i_spec=1,dust%n_species

            !calculate difference between actual wavelength...
            !...and wavelength bins in order to calculate which bin photon is in for each species

            !determine wavelength for extinction calculation in particle rest frame
            call inv_lorentz_trans_for_ext(packet%vel_vect,packet%dir_cart,packet%nu,nu_ext)
            wav_id=minloc(abs((dust%species(i_spec)%wav(:)-(c*1e6/nu_ext))),1)

            !calculate opactiy as function of rho for specific wavelength by interpolating between bins for each species
            !!include interpolate function?
            if ((c*1e6/nu_ext-dust%species(i_spec)%wav(wav_id))<0) then
                c_ext(i_spec)=dust%species(i_spec)%c_ext(wav_id)-((dust%species(i_spec)%c_ext(wav_id)-dust%species(i_spec)%c_ext(wav_id-1))* &
                    & ((dust%species(i_spec)%wav(wav_id)-c*1e6/nu_ext)/(dust%species(i_spec)%wav(wav_id)-dust%species(i_spec)%wav(wav_id-1))))
            else
                c_ext(i_spec)=dust%species(i_spec)%c_ext(wav_id)+((dust%species(i_spec)%c_ext(wav_id+1)-dust%species(i_spec)%c_ext(wav_id))* &
                    & ((c*1e6/nu_ext-dust%species(i_spec)%wav(wav_id))/(dust%species(i_spec)%wav(wav_id+1)-dust%species(i_spec)%wav(wav_id))))
            end if

            !cumulative scattering component of extinction for albedo calculation
            if ((c*1e6/nu_ext-dust%species(i_spec)%wav(wav_id))<0) then
                c_sca(i_spec)=dust%species(i_spec)%c_sca(wav_id)-((dust%species(i_spec)%c_sca(wav_id)-dust%species(i_spec)%c_sca(wav_id-1))* &
                    & ((dust%species(i_spec)%wav(wav_id)-c*1e6/nu_ext)/(dust%species(i_spec)%wav(wav_id)-dust%species(i_spec)%wav(wav_id-1))))
            else
                c_sca(i_spec)=dust%species(i_spec)%c_sca(wav_id)+((dust%species(i_spec)%c_sca(wav_id+1)-dust%species(i_spec)%c_sca(wav_id))* &
                    & ((c*1e6/nu_ext-dust%species(i_spec)%wav(wav_id))/(dust%species(i_spec)%wav(wav_id+1)-dust%species(i_spec)%wav(wav_id))))
            end if

            !cumulative scattering component of extinction for albedo calculation
            if ((c*1e6/nu_ext-dust%species(i_spec)%wav(wav_id))<0) then
                g_param(i_spec)=dust%species(i_spec)%g_param(wav_id)-((dust%species(i_spec)%g_param(wav_id)-dust%species(i_spec)%g_param(wav_id-1))* &
                    & ((dust%species(i_spec)%wav(wav_id)-c*1e6/nu_ext)/(dust%species(i_spec)%wav(wav_id)-dust%species(i_spec)%wav(wav_id-1))))
            else
                g_param(i_spec)=dust%species(i_spec)%g_param(wav_id)+((dust%species(i_spec)%g_param(wav_id+1)-dust%species(i_spec)%g_param(wav_id))* &
                    & ((c*1e6/nu_ext-dust%species(i_spec)%wav(wav_id))/(dust%species(i_spec)%wav(wav_id+1)-dust%species(i_spec)%wav(wav_id))))
            end if
        end do

        !calculate total opactiies weighted over all species
        c_ext_tot=sum(c_ext*dust%species%weight)
        c_sca_tot=sum(c_sca*dust%species%weight)
        g_param_tot=sum(g_param*dust%species%weight)

        !calculate albedo (don't add weighted albedos, must add each component and then divide total scat by total ext)
        albedo=c_sca_tot/c_ext_tot

        !calculate overall opactiy using rho
        !kappa_rho = kappa*rho = n*cext (units of cm^-1)
        kappa_rho=(c_ext_tot*grid_cell(packet%cell_no)%n_rho)

    end subroutine

    !---------------------------------------------------------------------

    subroutine calculate_velocity()
      !calculate bulk velocity of scattering e- and velocity unit vector
      if ((lg_vel_law) .and. (.not. lg_decoupled)) then
         ran = r4_uni_01()
         packet%v=(random(1)*(gas_geometry%v_max**(1+gas_geometry%v_prob_indx)- &
              & gas_geometry%v_min**(1+gas_geometry%v_prob_indx))+ &
              & gas_geometry%v_min**(1+gas_geometry%v_prob_indx))**(1/(1+gas_geometry%v_prob_indx))
      else
         packet%v=(random(1)*(gas_geometry%v_max**(1+gas_geometry%v_prob_indx)-&
              & gas_geometry%v_min**(1+gas_geometry%v_prob_indx))+ &
              & gas_geometry%v_min**(1+gas_geometry%v_prob_indx))**(1/(1+gas_geometry%v_prob_indx))
      end if
      packet%vel_vect=normalise(packet%pos_cart)*packet%v
    end subroutine

    !---------------------------------------------------------------------

    subroutine scatter()
    use rnglib
        !sample new propagation direction
        !call random_number(random)
        random(1) = r4_uni_01()
        random(2) = r4_uni_01()
        random(3) = r4_uni_01()
        random(4) = r4_uni_01()
        random(5) = r4_uni_01()

        if (dust%scat_type == 'hg' .and. event_type == 'dust_scat') then
            packet%dir_sph(:)=(/ (1.0/(2*g_param_tot))*(1+g_param_tot**2-((1-g_param_tot**2)/(1-g_param_tot+2*g_param_tot*random(1)))**2), &
                & random(2)*2*pi /)
        else
            packet%dir_sph(:)=(/ ((2*random(1))-1),random(2)*2*pi /)
        end if
        packet%dir_cart(:)=cart(acos(packet%dir_sph(1)),packet%dir_sph(2))
    end subroutine

    !---------------------------------------------------------------------

    subroutine update_cell_no()
        !calculate overall id of cell using x,y and z ids
        !note that in list of all cells, cells listed changing first z, then y, then x e.g.
        ! 1 1 1, 1 1 2, 1 1 3, 1 2 1, 1 2 2, 1 2 3, 2 1 1, 2 1 2... etc.
        packet%cell_no=(mothergrid%n_cells(2)*mothergrid%n_cells(3)*(packet%axis_no(1)-1)+mothergrid%n_cells(3)*(packet%axis_no(2)-1)+packet%axis_no(3))
    end subroutine

    !---------------------------------------------------------------------

    subroutine check_packet_in_right_cell()
        do i_dir=1,3
            if (packet%pos_cart(i_dir)<grid_cell(packet%cell_no)%axis(i_dir)) then
                !idgp(i_dir)=idgp(i_dir)-1
                !igpp=(mothergrid%n_cells(2)*mothergrid%n_cells(3)*(idgp(1)-1))+mothergrid%n_cells(3)*(idgp(2)-1)+idgp(3)
                print*,'error - packet coordinates are not in the identified cell. packet removed.'
                !print*,'1',i_dir,grid_cell(packet%cell_no)%axis(i_dir),packet%pos_cart(i_dir),grid_cell(packet%cell_no)%axis(i_dir)+grid_cell%width(i_dir)
                packet%lg_active=.false.
                return
            else if (packet%pos_cart(i_dir)>(grid_cell(packet%cell_no)%axis(i_dir)+grid_cell(packet%cell_no)%width(i_dir))) then
                !idgp(i_dir)=idgp(i_dir)+1
                !igpp=(mothergrid%n_cells(2)*mothergrid%n_cells(3)*(idgp(1)-1))+mothergrid%n_cells(3)*(idgp(2)-1)+idgp(3)
                print*,'error - packet coordinates are not in the identified cell. packet removed.'
                !print*,'2',packet%axis_no(i_min),i_dir,i_min,packet%pos_cart(i_dir),grid_cell(packet%cell_no)%axis(i_dir)+grid_cell(packet%cell_no)%width(i_dir),ip!,grid_cell(packet%cell_no)%axis(i_dir),packet%pos_cart(i_dir),grid_cell(packet%cell_no)%axis(i_dir)+grid_cell%width(i_dir)
                packet%lg_active=.false.
                return
            end if
        end do
    end subroutine

    !---------------------------------------------------------------------

    subroutine check_escaped()
        if (packet%axis_no(i_min) > mothergrid%n_cells(1) .or. packet%axis_no(i_min)<1 .or. packet%r>(max(gas_geometry%r_max,dust_geometry%r_max)*1e15)) then
            call check_los()
            return
        end if
    end subroutine

    !---------------------------------------------------------------------

    subroutine check_step_no()
        !check number of steps of this packet - a packet could get stuck in highly scattering environments
        packet%step_no=packet%step_no+1
        if (packet%step_no>5000) then
            packet%lg_active=.false.
            print*,'Packet reached 5000 steps and was removed - either highly scattered or very large grid.'
            return
        end if
    end subroutine

    !---------------------------------------------------------------------

    subroutine check_los()
        if (acos(packet%pos_cart(3)/sum(packet%pos_cart**2)**0.5) < pi/6) then
            packet%lg_los = .true.
        end if
    end subroutine

    !---------------------------------------------------------------------

    subroutine determine_event_type()
    use rnglib

!        call random_number(ran)
        ran = r4_uni_01()

        !if es used then establish whether dust event or e- scattering event
        if ((.not. lg_es) .or. (ran<kappa_rho/(kappa_rho+sigma_t*grid_cell(packet%cell_no)%n_e))) then
            !dust event - either scattering or absorption...

            !generate random number to compare to dust albedo in order to determine whether dust absorption or scattering

            ran = r4_uni_01()

            if (ran<albedo) then
                !dust scattering event
                event_type='dust_scat'
            else
                event_type='dust_absn'
            end if
        else
            event_type='elec_scat'
        end if
    end subroutine

    !---------------------------------------------------------------------

end module
