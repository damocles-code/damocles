MODULE init_packet

  use input
  USE initialise
  USE vector_functions
  
  IMPLICIT NONE
    
contains

  SUBROUTINE emit_photon(nu_p,dir_cart,pos_cart,iG_axis,lgactive,w,cell_iG)

    INTEGER               ::  ixx,iyy,izz
    INTEGER,INTENT(IN) ::  cell_iG
    
    REAL,INTENT(INOUT)    ::  nu_p,dir_cart(3)
    REAL,INTENT(OUT)      ::  pos_cart(3),w
    INTEGER,INTENT(OUT)   ::  iG_axis(3),lgactive
    
    
    REAL    ::  RANDOM(5)
    REAL    ::  POS_SPH(3)
    REAL    ::  DIR_SPH(2)
    REAL    ::  v_p,vel_vect(3)
    
    CALL RANDOM_NUMBER(RANDOM)
    
    w=1         !w is weight of photon - energy scaled when doppler shifted
    lgactive=0  !automatically inactive until declared active
    
    !initial position of packet is generated in both cartesian and spherical coordinates
    IF ((gas_geometry%clumped_mass_frac==1) .or. (gas_geometry%type == "arbitrary")) THEN
       !packets are emitted from grid cells
       POS_CART(:)=(/ (grid_cell(cell_iG)%axis(1)+RANDOM(1)*mothergrid%cell_width(1)),(grid_cell(cell_iG)%axis(2)+RANDOM(2)*mothergrid%cell_width(2)),(grid_cell(cell_iG)%axis(3)+RANDOM(3)*mothergrid%cell_width(3))/)
       POS_SPH(1)=((POS_CART(1)**2+POS_CART(2)**2+POS_CART(3)**2)**0.5)*1e-15
       POS_SPH(2)=ATAN(POS_CART(2)/POS_CART(1))
       POS_SPH(3)=ACOS(POS_CART(3)*1e-15/POS_SPH(1))
       POS_CART(:)=POS_CART(:)*1e-15
    ELSE
       !shell emissivity distribution
       POS_SPH(:)=(/ (RANDOM(1)*shell_width+RSh(iSh,1)),(2*RANDOM(2)-1),RANDOM(3)*2*pi/)       !position of emitter idP - spherical coords - system SN - RF
       POS_CART(:)=cartr(POS_SPH(1),ACOS(POS_SPH(2)),POS_SPH(3))
    END IF
    
    !generate an initial propagation direction from an isotropic distribution
    !in comoving frame of emitting particle
    DIR_SPH(:)=(/ (2*RANDOM(4))-1,RANDOM(5)*2*pi /)
    DIR_CART(:)=cart(ACOS(DIR_SPH(1)),DIR_SPH(2))
    
    IF (((POS_SPH(1) > gas_geometry%R_min) .AND. (POS_SPH(1) < gas_geometry%R_max) .AND. (gas_geometry%clumped_mass_frac==0)) .OR. (gas_geometry%clumped_mass_frac==1)) THEN                       !If the photon lies inside the radial bounds of the supernova then it is processed       
       
       !calculate velocity of emitting particle from radial velocity distribution
       !velocity vector comes from radial position vector of particle
       v_p=gas_geometry%v_max*((POS_SPH(1)/gas_geometry%R_max)**gas_geometry%v_power)
       vel_vect=normalise(pos_cart)*v_p
       
       nu_p=line%frequency
       lgactive=1
       
       call lorentz_trans(vel_vect,nu_p,DIR_CART,w,"emsn")
       
       !identify cell which contains emitting particle (and therefore packet)
       DO ixx=1,mothergrid%n_cells(1)
          IF ((POS_CART(1)*1e15-mothergrid%x_div(ixx))<0) THEN                        !identify grid axis that lies just beyond position of emitter in each direction
             iG_axis(1)=ixx-1                                       !then the grid cell id is the previous one
             EXIT
          END IF
          IF (ixx==mothergrid%n_cells(1)) THEN
             iG_axis(1)=mothergrid%n_cells(1)
          END IF
          
       END DO
       DO iyy=1,mothergrid%n_cells(2)
          IF ((POS_CART(2)*1e15-mothergrid%y_div(iyy))<0) THEN
             iG_axis(2)=iyy-1
             EXIT
          END IF
          IF (iyy==mothergrid%n_cells(2)) THEN
             iG_axis(2)=mothergrid%n_cells(2)
          END IF
          
       END DO
       DO izz=1,mothergrid%n_cells(3)
          IF ((POS_CART(3)*1e15-mothergrid%z_div(izz))<0) THEN
             iG_axis(3)=izz-1
             EXIT
          END IF
          IF (izz==mothergrid%n_cells(3)) THEN
             iG_axis(3)=mothergrid%n_cells(3)
          END IF
       END DO
       
       PRINT*,'short term test here',iG_axis(:)
       PRINT*,grid_cell(cell_iG)%axis(1)
       
    ELSE                                                    !If the photon lies outside the bounds of the SN then it is inactive and not processed
       n_inactive=n_inactive+1                             !add 1 to number of inactive photons
       PRINT*,'inactive photon',cell_iG,n_inactive,iDoublet
       lgactive=0
    END IF
    
    IF (ANY(iG_axis == 0)) lgactive=0
    
  END SUBROUTINE emit_photon
  
END MODULE init_packet
