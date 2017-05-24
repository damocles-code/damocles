RECURSIVE SUBROUTINE propagate(nu_p,dir_cart,pos_cart,iG_axis,iG,lgabs,lgactive,w,scatno)

    use input
    use grain_sizes
    use initialise
    use vector_functions
    use init_packet
    use random_routines
    use electron_scattering

    IMPLICIT NONE


    INTEGER,INTENT(OUT)	    ::	lgabs, lgactive
    INTEGER,INTENT(INOUT)   ::  iG,iG_axis(3)

    REAL,INTENT(INOUT)      ::  nu_p
    REAL,INTENT(INOUT)      ::  dir_cart(3),pos_cart(3),w

    REAL    ::  tau,kp,kappa_i,kappa_sca_i
    REAL    ::  kappa,kappa_sca,albdo
    REAL    ::  ran,random(2)
    REAL    ::  r
    REAL    ::  sface(3),s,s_min
    REAL    ::  v,vel_vect(3)
    REAL    ::  dir_sph(2)

    REAL    ::  V_T(3)
    REAL    ::  maxwell_sigma

    REAL, PARAMETER :: sigma_T=6.652E-25 !(cm2)

    INTEGER ::  i_dir,iSG,ispec
    INTEGER ::  wavid(1),wavid2
    INTEGER ::  iSGP(3)
    INTEGER ::  f,imin(1)
    INTEGER ::  scatno
    !INTEGER ::  omp_get_thread_num
    
    maxwell_sigma=((ES_temp*1.51563e7)**0.5)/1000


    !PRINT*,scatno!,omp_get_thread_num()

    !calculate overall id of cell using x,y and z ids
    !note that in list of all cells, cells listed changing first z, then y, then x e.g.
    ! 1 1 1, 1 1 2, 1 1 3, 1 2 1, 1 2 2, 1 2 3, 2 1 1, 2 1 2... etc.
    iG=(ncells(2)*ncells(3)*(iG_axis(1)-1)+ncells(3)*(iG_axis(2)-1)+iG_axis(3))
    
    !initialise absorption logical to 0.  This changed to 1 is absorbed.
    !Note absorption(lgabs - absorbed) different to inactive(lginactive-never emitted).
    lgabs=0
    albdo=0
    kappa=0
    kappa_sca=0
    DO ispec=1,nspecies
        !Calculate difference between actual wavelength...
        !...and wavelength bins in order to calculate which bin photon is in for each species
        wavid=MINLOC(ABS((dust(ispec)%wav(:)-(c*1e6/nu_p))))
        wavid2=wavid(1)

        !calculate opactiy as function of rho for specific wavelength by interpolating between bins for each species
        IF ((c*1e6/nu_p-dust(ispec)%wav(wavid2))<0) THEN
            kappa_i=dust(ispec)%ext_opacity(wavid2)-((dust(ispec)%ext_opacity(wavid2)-dust(ispec)%ext_opacity(wavid2-1))* &
                & ((dust(ispec)%wav(wavid2)-c*1e6/nu_p)/(dust(ispec)%wav(wavid2)-dust(ispec)%wav(wavid2-1))))
        ELSE
            kappa_i=dust(ispec)%ext_opacity(wavid2)+((dust(ispec)%ext_opacity(wavid2+1)-dust(ispec)%ext_opacity(wavid2))* &
                & ((c*1e6/nu_p-dust(ispec)%wav(wavid2))/(dust(ispec)%wav(wavid2+1)-dust(ispec)%wav(wavid2))))
        END IF

        !add this weighted opacity to overall opacity
        kappa=kappa+kappa_i*dust(ispec)%weight

        !calculate albedo as well
        !NOTE THAT ONE DOES NOT JUST WEIGHT THE ALBEDO 
                IF ((c*1e6/nu_p-dust(ispec)%wav(wavid2))<0) THEN
            kappa_sca_i=dust(ispec)%sca_opacity(wavid2)-((dust(ispec)%sca_opacity(wavid2)-dust(ispec)%sca_opacity(wavid2-1))* &
                & ((dust(ispec)%wav(wavid2)-c*1e6/nu_p)/(dust(ispec)%wav(wavid2)-dust(ispec)%wav(wavid2-1))))
        ELSE
            kappa_sca_i=dust(ispec)%sca_opacity(wavid2)+((dust(ispec)%sca_opacity(wavid2+1)-dust(ispec)%sca_opacity(wavid2))* &
                & ((c*1e6/nu_p-dust(ispec)%wav(wavid2))/(dust(ispec)%wav(wavid2+1)-dust(ispec)%wav(wavid2))))
        END IF
        kappa_sca=kappa_sca+kappa_sca_i*dust(ispec)%weight
        
END DO
albdo=kappa_sca/kappa

    !call random number and sample from CFD to obtain tau
    call random_number(ran)
    tau=-(ALOG((1-ran)))

    !Calculate overall opactiy using rho and work out distance packet will travel
    
!    IF (mgrid(iG)%cellStatus ==1) THEN
!        kp=kappa*mgrid(iG)%subgrid(iSG)%nrho       !=k*rho=n*Cext, units of cm^-1
!        !If density of cell is 0 then set distance such that packet guaranteed to reach cell boundary
!
!        IF (mgrid(iG)%subgrid(iSG)%nrho>0) THEN
!            IF (lgES) THEN
!                s=tau/(kp+sigma_T*mgrid(iG)%N_e)
!            ELSE
!                s=tau/kp
!            END IF                !distance therefore in units of cm
!        ELSE
!            s=(width(1))*(3**0.5)
!        END IF
!    ELSE
        kp=kappa*mgrid(iG)%nrho     !=k*rho=n*Cext, units of cm^-1
        

        IF (mgrid(iG)%nrho>0) THEN
           
            IF (lgES) THEN
                s=tau/(kp+sigma_T*mgrid(iG)%N_e)
            !    PRINT*,s
            ELSE
                s=tau/kp                !distance therefore in units of cm
            END IF
        ELSE

           IF (lgES) THEN
              
                s=tau/(sigma_T*mgrid(iG)%N_e)
                !PRINT*,tau,sigma_T,mgrid(iG)%N_e,iG
            ELSE
              
                s=width(1)*(3**0.5)
                
            END IF
        END IF
!    END IF
    !PRINT*,s,s_min
    dir_cart=normalise(dir_cart)

    !calculate distance to nearest face (and identify whether in subgrid cell or mothergrid cell)
!    IF (mgrid(iG)%cellStatus ==1) THEN
!
!        DO i_dir=1,3
!            IF (pos_cart(i_dir) > mgrid(iG)%subaxes(i_dir,2)) THEN
!                iSGP(i_dir)=2
!            ELSE
!                iSGP(i_dir)=1
!            END IF
!        END DO
!        iSG=(iSGP(1)-1)*4+(iSGP(2)-1)*2+isGP(3)
!        DO i_dir=1,3
!            IF (dir_cart(i_dir)<0) THEN
!                sface(i_dir)=ABS((mgrid(iG)%subgrid(iSG)%axis(i_dir)-pos_cart(i_dir))/dir_cart(i_dir))
!            ELSE
!                sface(i_dir)=ABS((mgrid(iG)%subgrid(iSG)%axis(i_dir)+width(i_dir)/2-pos_cart(i_dir))/dir_cart(i_dir))
!            END IF
!        END DO
!    ELSE
        DO i_dir=1,3
            IF (dir_cart(i_dir)<0) THEN
                sface(i_dir)=ABS((mgrid(iG)%axis(i_dir)-pos_cart(i_dir))/dir_cart(i_dir))
            ELSE
                sface(i_dir)=ABS((mgrid(iG)%axis(i_dir)+width(i_dir)-pos_cart(i_dir))/dir_cart(i_dir))
            END IF
        END DO
!    END IF
    s_min=MINVAL(sface)

    imin=MINLOC(sface,1)
    
    !event occurs when distance travelled is > distance to nearest face
    
   
    IF ((s>s_min)) THEN

        !packet travels to cell boundary
        !direction of travel remains the same
        pos_cart(:)=pos_cart(:)+(ABS(s_min)+ABS(s_min)*1E-10)*dir_cart(:)     !actually moves just past boundary by small factor...
        f=imin(1)
        IF (dir_cart(f)>0) THEN
            IF (iG_axis(f) /= ncells(1)) THEN
                iG_axis(imin)=iG_axis(imin)+1
            ELSE 
                !PRINT*, scatno,w
                RETURN
            END IF
            iG=(ncells(2)*ncells(3)*(iG_axis(1)-1))+ncells(3)*(iG_axis(2)-1)+iG_axis(3)
            pos_cart(f)=mgrid(iG)%axis(f)+((ABS(s_min)*1E-10)*dir_cart(f))
        ELSE
            IF (iG_axis(f) /= 1) THEN
                iG_axis(imin)=iG_axis(imin)-1
            ELSE
                !PRINT*, scatno,w
                RETURN
            END IF
            iG=(ncells(2)*ncells(3)*(iG_axis(1)-1))+ncells(3)*(iG_axis(2)-1)+iG_axis(3)
            IF (iG<1) PRINT*,'here',iG
            pos_cart(f)=mgrid(iG)%axis(f)+((ABS(s_min)*1E-10)*dir_cart(f))+width(1)
        END IF
        r=(pos_cart(1)**2+pos_cart(2)**2+pos_cart(3)**2)**0.5



        !TEST FOR PACKET IN CORRECT CELL
        !IF (omp_get_thread_num()==0) THEN
        DO i_dir=1,3
            IF (POS_cart(i_dir)<mgrid(iG)%axis(i_dir)) THEN
                !idGP(i_dir)=idGP(i_dir)-1
                !iGPP=(ncells(2)*ncells(3)*(idGP(1)-1))+ncells(3)*(idGP(2)-1)+idGP(3)
                PRINT*,'yelp!','1',i_dir,mgrid(iG)%axis(i_dir),POS_cart(i_dir),mgrid(iG)%axis(i_dir)+width(i_dir)
                RETURN
            ELSE IF (POS_cart(i_dir)>(mgrid(iG)%axis(i_dir)+width(i_dir))) THEN
                !idGP(i_dir)=idGP(i_dir)+1
                !iGPP=(ncells(2)*ncells(3)*(idGP(1)-1))+ncells(3)*(idGP(2)-1)+idGP(3)
                PRINT*,'yelp!','2',i_dir,mgrid(iG)%axis(i_dir),POS_cart(i_dir),mgrid(iG)%axis(i_dir)+width(i_dir)
                RETURN
            END IF
        END DO
        !END IF

        !If moved outside of cells then packet escaped
        IF (iG_axis(f) > ncells(1) .OR. iG_axis(f)<1 .OR. r>(MAX(R_max_gas,R_max)*1E15)) THEN
            !PRINT*, scatno,w
            RETURN
        END IF

        CALL propagate(nu_p,dir_cart,pos_cart,iG_axis,iG,lgabs,lgactive,w,scatno)
        !PRINT*,'here3'
    ELSE
       
        !event does occur.

        !calculate position and radius of event
        pos_cart(:)=pos_cart(:)+s*dir_cart(:)
        r=(pos_cart(1)**2+pos_cart(2)**2+pos_cart(3)**2)**0.5

        !if ES used then establish whether dust event or e- scattering event

        call random_number(ran)

        IF (scatno>500) THEN
           lgactive=0
           !lgabs=1
           PRINT*, scatno
            RETURN
       END IF
       


        IF ((.not. lgES) .OR. (ran<kp/(kp+sigma_T*mgrid(iG)%N_e))) THEN
           
           !PRINT*,albdo
            !PRINT*,'here?'

            !            IF (ran<0.0) THEN
            call random_number(ran)

!            IF (ran<1.0) THEN
 
            IF (ran<albdo) THEN
                scatno=scatno+1
!                PRINT*, scatno
                !calculate velocity of scatterer and velocity unit vector
                v=v_max*((r/(R_max*1e15))**l)
                vel_vect=normalise(pos_cart)*v
                IF (r>MAX(R_max_gas,R_max)*1e15) THEN
                    !NOTE no actual scattering as has already escaped
                    !PRINT*, scatno,w
!                   PRINT*,'escaped' 
                   RETURN
                END IF
                !PRINT*,'before',scatno,nu_p,1/nu_p
                IF (lgVelShift) THEN
                    !Inverse lorentz boost for doppler shift of packet hitting particle
                    call inv_lorentz_trans(vel_vect,nu_p,dir_cart,w,"scat")
                   
                END IF
!                PRINT*,'middle',scatno,nu_p,1/nu_p
                !Now scatter (sample new direction)
                call random_number(random)
                dir_sph(:)=(/ ((2*random(1))-1),random(2)*2*pi /)                       !new scattering angle (i.e. angle of exit) - spherical coordinates - system PT - CMF
                dir_cart(:)=cart(ACOS(dir_sph(1)),dir_sph(2))

                IF (lgVelShift) THEN
                    !Lorentz boost for doppler shift of packet bouncing off particle
                    call lorentz_trans(vel_vect,nu_p,dir_cart,w,"scat")
                END IF
                !PRINT*,'scat'
                
 !               PRINT*,'after',scatno,nu_p,1/nu_p
                call propagate(nu_p,dir_cart,pos_cart,iG_axis,iG,lgabs,lgactive,w,scatno)

            !PRINT*,'here1'
            ELSE
                !absorbed
!                print*,'abs'
                lgabs=1
                !PRINT*, scatno,w
                RETURN
            END IF
        ELSE

            !PRINT*,'here2'
            !electron scattering event
            !note there is no weighting due to albedo here (since no absorption)
            !calculate bulk velocity of scattering e- and velocity unit vector

            v=v_max*((r/(R_max*1e15))**l)
            vel_vect=normalise(pos_cart)*v

            !also calculate thermal velocity and random thermal velocity vector
            V_T(1)=normal(0.0,maxwell_sigma)
            !PRINT*,V_T(1)
            V_T(2)=normal(0.0,maxwell_sigma)
            V_T(3)=normal(0.0,maxwell_sigma)

            !PRINT*, V_T(1),V_T(2),V_T(3)

            vel_vect(1)=vel_vect(1)+V_T(1)
            vel_vect(2)=vel_vect(2)+V_T(2)
            vel_vect(3)=vel_vect(3)+V_T(3)

            !vel_vect(1)=V_T(1)
            !vel_vect(2)=V_T(2)
            !vel_vect(3)=V_T(3)

            !PRINT*,(vel_vect(1)**2+vel_vect(2)**2+vel_vect(3)**2)**0.5

            IF (r>MAX(R_max_gas,R_max)*1e15) THEN
                !NOTE no actual scattering as has already escaped
                !PRINT*, scatno,w
                !PRINT*,'escaped'
                RETURN
            END IF
            !PRINT*, 'v  ',v
            !PRINT*, 'nu1',nu_p
            IF (lgVelShift) THEN
                !Inverse lorentz boost for doppler shift of packet hitting particle
                call inv_lorentz_trans(vel_vect,nu_p,dir_cart,w,"escat")
            END IF

            !Now scatter (sample new direction)
            call random_number(random)
            dir_sph(:)=(/ ((2*random(1))-1),random(2)*2*pi /)                       !new scattering angle (i.e. angle of exit) - spherical coordinates - system PT - CMF
            dir_cart(:)=cart(ACOS(dir_sph(1)),dir_sph(2))


            !CHECK WHETHER THIS SHOULD BE WEIGHTED OR NOT... escat or scat?
            IF (lgVelShift) THEN
                !Lorentz boost for doppler shift of packet bouncing off particle
                call lorentz_trans(vel_vect,nu_p,dir_cart,w,"escat")
            END IF
            !PRINT*,'scat'
            !PRINT*, 'nu2',nu_p
            call propagate(nu_p,dir_cart,pos_cart,iG_axis,iG,lgabs,lgactive,w,scatno)

        END IF
       
    END IF

!r=(pos_cart(1)**2+pos_cart(2)**2+pos_cart(3)**2)**0.5
!v=v_max*((r/(R_max*1e5))**l)
!vel_vect=normalise(pos_cart)*v


!call inv_lorentz_trans(vel_vect,nu_p,dir_cart,w,"scat")

END SUBROUTINE propagate
