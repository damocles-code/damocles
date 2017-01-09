RECURSIVE SUBROUTINE propagate(nu_p,dir_cart,pos_cart,iG_axis,lgabs,lgactive,w,scatno)

    USE globals
    USE class_dust
    USE class_grid
    USE input
    USE initialise
    USE vector_functions
    USE init_packet
    USE random_routines
    USE electron_scattering
    USE class_packet

    IMPLICIT NONE


    INTEGER,INTENT(OUT)     ::  lgabs, lgactive
    INTEGER,INTENT(INOUT)   ::  iG_axis(3)

    REAL,INTENT(INOUT)      ::  nu_p
    REAL,INTENT(INOUT)      ::  dir_cart(3),pos_cart(3),w

    REAL    ::  tau,kp,kappa_i,kappa_sca_i
    REAL    ::  kappa,kappa_sca,albdo
    REAL    ::  ran,random(2)

    REAL    ::  sface(3),s,s_min
    REAL    ::  v,vel_vect(3)
    REAL    ::  dir_sph(2)

    REAL    ::  V_T(3)
    REAL    ::  maxwell_sigma

    REAL, PARAMETER :: sigma_T=6.652E-25 !(cm2)

    INTEGER ::  i_dir,ispec
    INTEGER ::  wav_id
    INTEGER ::  iSGP(3)
    INTEGER ::  imin
    INTEGER ::  scatno
    !INTEGER ::  omp_get_thread_num
    
    maxwell_sigma=((ES_temp*1.51563e7)**0.5)/1000

    !calculate overall id of cell using x,y and z ids
    !note that in list of all cells, cells listed changing first z, then y, then x e.g.
    ! 1 1 1, 1 1 2, 1 1 3, 1 2 1, 1 2 2, 1 2 3, 2 1 1, 2 1 2... etc.
    packet%iG=(mothergrid%n_cells(2)*mothergrid%n_cells(3)*(iG_axis(1)-1)+mothergrid%n_cells(3)*(iG_axis(2)-1)+iG_axis(3))

    !initialise absorption logical to 0.  This changed to 1 if absorbed.
    !Note absorption (lgabs - absorbed) different to inactive (lginactive-never emitted).
    lgabs=0
    albdo=0
    kappa=0
    kappa_sca=0

    DO ispec=1,dust%n_species

        !Calculate difference between actual wavelength...
        !...and wavelength bins in order to calculate which bin photon is in for each species
        wav_id=MINLOC(ABS((dust%species(ispec)%wav(:)-(c*1e6/nu_p))),1)

        !calculate opactiy as function of rho for specific wavelength by interpolating between bins for each species
        IF ((c*1e6/nu_p-dust%species(ispec)%wav(wav_id))<0) THEN
            kappa_i=dust%species(ispec)%ext_opacity(wav_id)-((dust%species(ispec)%ext_opacity(wav_id)-dust%species(ispec)%ext_opacity(wav_id-1))* &
                & ((dust%species(ispec)%wav(wav_id)-c*1e6/nu_p)/(dust%species(ispec)%wav(wav_id)-dust%species(ispec)%wav(wav_id-1))))
        ELSE
            kappa_i=dust%species(ispec)%ext_opacity(wav_id)+((dust%species(ispec)%ext_opacity(wav_id+1)-dust%species(ispec)%ext_opacity(wav_id))* &
                & ((c*1e6/nu_p-dust%species(ispec)%wav(wav_id))/(dust%species(ispec)%wav(wav_id+1)-dust%species(ispec)%wav(wav_id))))
        END IF

        !add this weighted opacity to overall opacity
        kappa=kappa+kappa_i*dust%species(ispec)%weight

        !cumulative scattering component of extinction for albedo calculation
        IF ((c*1e6/nu_p-dust%species(ispec)%wav(wav_id))<0) THEN
            kappa_sca_i=dust%species(ispec)%sca_opacity(wav_id)-((dust%species(ispec)%sca_opacity(wav_id)-dust%species(ispec)%sca_opacity(wav_id-1))* &
                & ((dust%species(ispec)%wav(wav_id)-c*1e6/nu_p)/(dust%species(ispec)%wav(wav_id)-dust%species(ispec)%wav(wav_id-1))))
        ELSE
            kappa_sca_i=dust%species(ispec)%sca_opacity(wav_id)+((dust%species(ispec)%sca_opacity(wav_id+1)-dust%species(ispec)%sca_opacity(wav_id))* &
                & ((c*1e6/nu_p-dust%species(ispec)%wav(wav_id))/(dust%species(ispec)%wav(wav_id+1)-dust%species(ispec)%wav(wav_id))))
        END IF
        kappa_sca=kappa_sca+kappa_sca_i*dust%species(ispec)%weight
        
    END DO

    !calculate albedo (don't add weighted albedos, must add each component and then divide total scat by total ext)
    albdo=kappa_sca/kappa

    !call random number and sample from CFD to obtain tau
    call random_number(ran)
    tau=-(ALOG((1-ran)))

    !Calculate overall opactiy using rho and work out distance packet will travel
    !kp = k*rho = n*Cext (units of cm^-1)
    kp=kappa*grid_cell(packet%iG)%nrho
        
    !work out potential distance travelled by packet based on optical depth tau
    !(distance in units of cm, since kp in units cm^-1)
    IF (grid_cell(packet%iG)%nrho>0) THEN
        IF (lg_ES) THEN
            s=tau/(kp+sigma_T*grid_cell(packet%iG)%N_e)
        ELSE
            s=tau/kp
        END IF
    ELSE
        !if cell has zero density then potential distance travelled is greatest distance across cell (for no e- scat)
        IF (lg_ES) THEN
            s=tau/(sigma_T*grid_cell(packet%iG)%N_e)
        ELSE
            s=(grid_cell(packet%iG)%width(1)**2+grid_cell(packet%iG)%width(2)**2+grid_cell(packet%iG)%width(3)**2)**0.5
        END IF
    END IF

    !unit vector direction of travel of packet in cartesian
    dir_cart=normalise(dir_cart)

    !calculate distance to nearest face
    DO i_dir=1,3
        IF (dir_cart(i_dir)<0) THEN
            sface(i_dir)=ABS((grid_cell(packet%iG)%axis(i_dir)-pos_cart(i_dir))/dir_cart(i_dir))
        ELSE
            sface(i_dir)=ABS((grid_cell(packet%iG)%axis(i_dir)+grid_cell(packet%iG)%width(i_dir)-pos_cart(i_dir))/dir_cart(i_dir))
        END IF
    END DO
    s_min=MINVAL(sface)
    imin=MINLOC(sface,1)    !index of nearest face (i.e. identifies whether nearest face is planar in x or y or z)
    
    !event occurs when distance travelled is < distance to nearest face
    IF ((s>s_min)) THEN
        !packet travels to cell boundary
        !direction of travel remains the same
        !position updated to be on boundary with next cell
        pos_cart(:)=pos_cart(:)+(ABS(s_min)+ABS(s_min)*1E-10)*dir_cart(:)     !actually moves just past boundary by small factor...

        IF (dir_cart(imin)>0) THEN
            IF (iG_axis(imin) /= mothergrid%n_cells(1)) THEN
                iG_axis(imin)=iG_axis(imin)+1

            ELSE 
                RETURN
            END IF
            packet%iG=(mothergrid%n_cells(2)*mothergrid%n_cells(3)*(iG_axis(1)-1))+mothergrid%n_cells(3)*(iG_axis(2)-1)+iG_axis(3)
            pos_cart(imin)=grid_cell(packet%iG)%axis(imin)+((ABS(s_min)*1E-10)*dir_cart(imin))
        ELSE
            IF (iG_axis(imin) /= 1) THEN
                iG_axis(imin)=iG_axis(imin)-1
            ELSE
                RETURN
            END IF
            packet%iG=(mothergrid%n_cells(2)*mothergrid%n_cells(3)*(iG_axis(1)-1))+mothergrid%n_cells(3)*(iG_axis(2)-1)+iG_axis(3)
            IF (packet%iG<1) PRINT*,'here',packet%iG
            pos_cart(imin)=grid_cell(packet%iG)%axis(imin)+((ABS(s_min)*1E-10)*dir_cart(imin))+grid_cell(packet%iG)%width(imin)
        END IF
        packet%r=(pos_cart(1)**2+pos_cart(2)**2+pos_cart(3)**2)**0.5

        !test that packet is in the correct cell
        DO i_dir=1,3
            IF (POS_cart(i_dir)<grid_cell(packet%iG)%axis(i_dir)) THEN
                !idGP(i_dir)=idGP(i_dir)-1
                !iGPP=(mothergrid%n_cells(2)*mothergrid%n_cells(3)*(idGP(1)-1))+mothergrid%n_cells(3)*(idGP(2)-1)+idGP(3)
                PRINT*,'yelp!','1',i_dir,grid_cell(packet%iG)%axis(i_dir),POS_cart(i_dir),grid_cell(packet%iG)%axis(i_dir)+grid_cell%width(i_dir)
                RETURN
            ELSE IF (POS_cart(i_dir)>(grid_cell(packet%iG)%axis(i_dir)+grid_cell(packet%iG)%width(i_dir))) THEN
                !idGP(i_dir)=idGP(i_dir)+1
                !iGPP=(mothergrid%n_cells(2)*mothergrid%n_cells(3)*(idGP(1)-1))+mothergrid%n_cells(3)*(idGP(2)-1)+idGP(3)
                PRINT*,'yelp!','2',iG_axis(imin),i_dir,imin,pos_cart(i_dir),grid_cell(packet%iG)%axis(i_dir)+grid_cell(packet%iG)%width(i_dir),iP!,grid_cell(packet%iG)%axis(i_dir),POS_cart(i_dir),grid_cell(packet%iG)%axis(i_dir)+grid_cell%width(i_dir)
                RETURN
            END IF
        END DO

        !If moved outside of cells then packet escaped
        IF (iG_axis(imin) > mothergrid%n_cells(1) .OR. iG_axis(imin)<1 .OR. packet%r>(MAX(gas_geometry%R_max,dust_geometry%R_max)*1E15)) THEN
            RETURN
        END IF

        CALL propagate(nu_p,dir_cart,pos_cart,iG_axis,lgabs,lgactive,w,scatno)

    ELSE
        !event does occur.
        !calculate position and radius of event
        pos_cart(:)=pos_cart(:)+s*dir_cart(:)
        packet%r=(pos_cart(1)**2+pos_cart(2)**2+pos_cart(3)**2)**0.5

        call random_number(ran)

        IF (scatno>500) THEN
            lgactive=0
            !lgabs=1
            PRINT*, scatno
            RETURN
        END IF
       
        !if ES used then establish whether dust event or e- scattering event
        IF ((.not. lg_ES) .OR. (ran<kp/(kp+sigma_T*grid_cell(packet%iG)%N_e))) THEN

            call random_number(ran)

            IF (ran<albdo) THEN

                scatno=scatno+1

                !calculate velocity of scatterer and velocity unit vector
                packet%v=dust_geometry%v_max*((packet%r/(dust_geometry%R_max*1e15))**dust_geometry%v_power)
                vel_vect=normalise(pos_cart)*packet%v

                IF (packet%r>MAX(gas_geometry%R_max,dust_geometry%R_max)*1e15) THEN
                     !NOTE no actual scattering as has already escaped
                    RETURN
                END IF

                IF (lg_vel_shift) THEN
                    !Inverse lorentz boost for doppler shift of packet hitting particle
                    call inv_lorentz_trans(vel_vect,nu_p,dir_cart,w,"scat")
                END IF

                !Now scatter (sample new direction)
                call random_number(random)
                dir_sph(:)=(/ ((2*random(1))-1),random(2)*2*pi /)                       !new scattering angle (i.e. angle of exit) - spherical coordinates - system PT - CMF
                dir_cart(:)=cart(ACOS(dir_sph(1)),dir_sph(2))

                IF (lg_vel_shift) THEN
                    !Lorentz boost for doppler shift of packet bouncing off particle
                    call lorentz_trans(vel_vect,nu_p,dir_cart,w,"scat")
                END IF

                call propagate(nu_p,dir_cart,pos_cart,iG_axis,lgabs,lgactive,w,scatno)

            ELSE
                !absorbed
                lgabs=1
                RETURN
            END IF
        ELSE
            !electron scattering event
            !note there is no weighting due to albedo here (since no absorption)
            !calculate bulk velocity of scattering e- and velocity unit vector

            packet%v=dust_geometry%v_max*((packet%r/(dust_geometry%R_max*1e15))**dust_geometry%v_power)
            vel_vect=normalise(pos_cart)*packet%v

            !also calculate thermal velocity and random thermal velocity vector
            V_T(1)=normal(0.d0,dble(maxwell_sigma))
            V_T(2)=normal(0.d0,dble(maxwell_sigma))
            V_T(3)=normal(0.d0,dble(maxwell_sigma))

            vel_vect(1)=vel_vect(1)+V_T(1)
            vel_vect(2)=vel_vect(2)+V_T(2)
            vel_vect(3)=vel_vect(3)+V_T(3)

            IF (packet%r>MAX(gas_geometry%R_max,dust_geometry%R_max)*1e15) THEN
                !NOTE no actual scattering as has already escaped
                RETURN
            END IF

            IF (lg_vel_shift) THEN
                !Inverse lorentz boost for doppler shift of packet hitting particle
                call inv_lorentz_trans(vel_vect,nu_p,dir_cart,w,"escat")
            END IF

            !Now scatter (sample new direction)
            call random_number(random)
            dir_sph(:)=(/ ((2*random(1))-1),random(2)*2*pi /)                       !new scattering angle (i.e. angle of exit) - spherical coordinates - system PT - CMF
            dir_cart(:)=cart(ACOS(dir_sph(1)),dir_sph(2))

            !CHECK WHETHER THIS SHOULD BE WEIGHTED OR NOT... escat or scat?
            IF (lg_vel_shift) THEN
                !Lorentz boost for doppler shift of packet bouncing off particle
                call lorentz_trans(vel_vect,nu_p,dir_cart,w,"escat")
            END IF

            call propagate(nu_p,dir_cart,pos_cart,iG_axis,lgabs,lgactive,w,scatno)

        END IF
       
    END IF

END SUBROUTINE propagate
