!---------------------------------------------------------------------------------!
!  this module declares the packet derived type object which descibes properties  !
!  of the packet that is currently being processed throught the ejecta            !
!                                                                                 !
!  the subroutine 'emit_packet' generates this packet and assigns its initial     !
!  properties and is called directly from the driver                              !
!---------------------------------------------------------------------------------!

module class_packet

    use globals
    use class_grid
    use input
    use initialise
    use vector_functions
    use rnglib

    implicit none

    type packet_obj
        real    ::  r            !current radius of packet location in cm
        real    ::  v            !current (scalar) velocity of packet in km/s
        real    ::  nu           !current frequency of packet
        real    ::  weight       !current weighting of packet

        real    ::  vel_vect(3)  !current velocity vector of packet in cartesian (km/s)
        real    ::  dir_cart(3)  !current direction of propagation of packet in cartesian coordinates
        real    ::  pos_cart(3)  !current position of packet in cartesian coordinates
        real    ::  dir_sph(2)   !current direction of propagation in spherical coordinates
        real    ::  pos_sph(3)   !current position of packet in spherical coordinates

        integer ::  cell_no      !current cell id of packet
        integer ::  axis_no(3)   !current cell ids in each axis
        integer ::  step_no      !number of steps that packet has experienced
        integer ::  freq_id      !id of frequency grid division that contains the frequency of the packet
        integer ::  id           !id number of packet
        integer ::  init_cell    !id number of cell in which packet initialised

        logical ::  lg_abs       !true indicates packet has been absored
        logical ::  lg_los       !true indicates packet is in line of sight (after it has escaped)
        logical ::  lg_active    !true indicates that packet is being or has been processed
                                 !false indidcates it was outside of ejecta and therefore inactive, or has been scattered too many times
    end type

    type(packet_obj) :: packet
    save packet
    !$OMP THREADPRIVATE(packet)
    real   ::  indx          !indx used in velocity calculation (TESTING)

contains

    !this subroutine generates a packet and samples an emission position in the observer's rest frame
    !a propagation direction is sampled from an isotropic distribution in the comoving frame of the emitter
    !and a frequency is assigned in this frame
    !frequency and propagation direction are updated in observer's rest frame and grid cell containing packet is identified
    subroutine emit_packet()
    use rnglib

        !call random_number(random)
        random(1) = r4_uni_01()
        random(2) = r4_uni_01()
        random(3) = r4_uni_01()
        random(4) = r4_uni_01()
        random(5) = r4_uni_01()

        !packets are weighted according to their frequency shift (energy is altered when doppler shifted)
        packet%weight=1

        !packet is declared inactive by default until declared active
        packet%lg_active=.false.

        !initialise absorption logical to 0.  this changed to 1 if absorbed.
        !note absorption (lg_abs = true - absorbed) different to inactive (lg_active = false - never emitted).
        packet%lg_abs=.false.

        !initialise step number to zero
        packet%step_no=0

        !initial position of packet is generated in both cartesian and spherical coordinates
        if ((gas_geometry%clumped_mass_frac==1) .or. (gas_geometry%type == "arbitrary")) then
           !packets are emitted from grid cells
           !locate the cell to emit the next packet from
           packet%init_cell = minloc((num_packets_array(:,2)-packet%id),1,((num_packets_array(:,2)-packet%id) .ge. 0))

           !calculate the starting position of the the packet at an arbitary location in that cell
            packet%pos_cart= (grid_cell(packet%init_cell)%axis+random(1:3)*grid_cell(packet%init_cell)%width)
            packet%pos_sph(1)=(((packet%pos_cart(1)/1e15)**2+(packet%pos_cart(2)/1e15)**2+(packet%pos_cart(3)/1e15)**2)**0.5)
            packet%pos_sph(2)=atan(packet%pos_cart(2)/packet%pos_cart(1))
            packet%pos_sph(3)=acos(packet%pos_cart(3)*1e-15/packet%pos_sph(1))
            packet%pos_cart(:)=packet%pos_cart(:)*1e-15
        else
            !shell emissivity distribution
            packet%pos_sph(1) = (((gas_geometry%r_min**(3-gas_geometry%emis_power*gas_geometry%rho_power)) &
            & +(random(1)*(gas_geometry%r_max**(3-gas_geometry%emis_power*gas_geometry%rho_power)-gas_geometry%r_min**(3-gas_geometry%emis_power*gas_geometry%rho_power)))) &
            & **(1.0/(3-gas_geometry%emis_power*gas_geometry%rho_power)))
            packet%pos_sph(2) = (2*random(2)-1)
            packet%pos_sph(3) = random(3)*2*pi
            packet%pos_cart(:)=cartr(packet%pos_sph(1),acos(packet%pos_sph(2)),packet%pos_sph(3))
        end if
        !generate an initial propagation direction from an isotropic distribution
        !in comoving frame of emitting particle
        packet%dir_sph(:)=(/ (2*random(4))-1,random(5)*2*pi /)
        packet%dir_cart(:)=cart(acos(packet%dir_sph(1)),packet%dir_sph(2))

        !if the photon lies inside the radial bounds of the supernova
        !or if the photon is emitted from a clump or cell (rather than shell) then it is processed
        if (((packet%pos_sph(1) > gas_geometry%r_min) .and. (packet%pos_sph(1) < gas_geometry%r_max) .and. (gas_geometry%clumped_mass_frac==0)) &
            & .or. (gas_geometry%clumped_mass_frac==1) &
            & .or. (gas_geometry%type == 'arbitrary')) then

           !calculate velocity of emitting particle from radial velocity distribution

           !if using a velocity law that is independent of radius, assign velocity here
           if (lg_vel_law) then
              random(1) = r4_uni_01()
              packet%v=(random(1)*(gas_geometry%v_max**(1+gas_geometry%v_prob_indx)- &
                   & gas_geometry%v_min**(1+gas_geometry%v_prob_indx))+ &
                   & gas_geometry%v_min**(1+gas_geometry%v_prob_indx))**(1/(1+gas_geometry%v_prob_indx))
           else
              !velocity vector comes from radial position vector of particle
              packet%v=gas_geometry%v_max*((packet%pos_sph(1)/gas_geometry%r_max)**gas_geometry%v_power)
           end if
           packet%vel_vect=normalise(packet%pos_cart)*packet%v
           packet%nu=line%frequency
           packet%lg_active=.true.

            call lorentz_trans(packet%vel_vect,packet%dir_cart,packet%nu,packet%weight,"emsn")

            !identify cell which contains emitting particle (and therefore packet)
            packet%axis_no(1) = minloc(packet%pos_cart(1)*1e15-mothergrid%x_div,1,(packet%pos_cart(1)*1e15-mothergrid%x_div)>0)
            packet%axis_no(2) = minloc(packet%pos_cart(2)*1e15-mothergrid%y_div,1,(packet%pos_cart(2)*1e15-mothergrid%y_div)>0)
            packet%axis_no(3) = minloc(packet%pos_cart(3)*1e15-mothergrid%z_div,1,(packet%pos_cart(3)*1e15-mothergrid%z_div)>0)        

            !check to ensure that for packets emitted from cells, the identified cell is the same as the original...
            if ((gas_geometry%type == 'shell' .and. gas_geometry%clumped_mass_frac == 1) &
                &    .or.  (gas_geometry%type == 'arbitrary')) then

                if ((packet%axis_no(1) /= grid_cell(packet%init_cell)%id(1)) .and. &
                    &   (packet%axis_no(2) /= grid_cell(packet%init_cell)%id(2)) .and. &
                    &   (packet%axis_no(3) /= grid_cell(packet%init_cell)%id(3))) then
                    print*,'WARNING: cell calculation gone wrong in module class_packet. Aborted.'
                    stop
                end if
            end if

        !if the photon lies outside the bounds of the sn then it is inactive and not processed
        else
            packet%lg_active=.false.
        end if

        if ((any(packet%axis_no == 0)) .and. (packet%lg_active)) then
            print*,'WARNING: inactive packet. Does not exist in grid.'
        end if

        packet%pos_cart=packet%pos_cart*1e15

    end subroutine emit_packet

end module class_packet
