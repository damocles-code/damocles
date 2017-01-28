!---------------------------------------------------------------------!
!  this module includes a number of functions to handle vectors       !
!  including - conversion from spherical coords to cartesian          !
!            - normalise a vector                                     !
!            - forwards and backwards lorentz transforms              !
!---------------------------------------------------------------------!

MODULE vector_functions

    use input

    IMPLICIT NONE

contains

    !function to convert spherical direction to cartesian unit vector
    FUNCTION cart(theta,phi)
        REAL    ::  cart(3)
        REAL    ::  theta,phi,x,y,z
        x=SIN(theta)*COS(phi)
        y=SIN(theta)*SIN(phi)
        z=COS(theta)
        cart=(/x,y,z/)
    END FUNCTION

    !function to convert spherical corrdinates to cartesian coordinates
    FUNCTION cartr(r,theta,phi)
        REAL    ::  cartr(3)
        REAL    ::  theta,phi,r,x,y,z
        x=r*SIN(theta)*COS(phi)
        y=r*SIN(theta)*SIN(phi)
        z=r*COS(theta)
        cartr=(/x,y,z/)
    END FUNCTION

    !normalise a cartesian vector to a unit vector
    FUNCTION normalise(cart)
        REAL    ::  normalise(3),cart(3)
        INTEGER ::  i_dir
        DO i_dir=1,3
            normalise(i_dir)=cart(i_dir)/((cart(1)**2+cart(2)**2+cart(3)**2)**0.5)
        END DO
    END FUNCTION

    !forwards lorentz transform to move from comoving frame to rest frame of observer
    !updates frequency, direction of propagation of packet and weight of packet
    SUBROUTINE lorentz_trans(vel_vect,dir_cart,nu,weight,str)
        REAL, INTENT(IN)                ::  vel_vect(3)
        REAL, INTENT(INOUT)             ::  dir_cart(3)
        REAL, INTENT(INOUT)             ::  nu,weight
        CHARACTER(LEN=4), INTENT(IN)    ::  str

        REAL                            ::  lorentz(4,4)
        REAL                            ::  gmma
        REAL                            ::  beta_vect(3),beta_scalar
        REAL                            ::  four_vect_in(4),four_vect_out(4)

        !standard beta and gamma definitions as per lorentz transforms
        beta_vect=(1e3*vel_vect)/c
        beta_scalar=(beta_vect(1)**2+beta_vect(2)**2+beta_vect(3)**2)**0.5
        gmma=1/(1-beta_scalar**2)**0.5
        !uncomment if ignoring terms > order(v^2/c^2)
        !gmma=1.

        !construct four vector and apply inverse lorentz transform
        four_vect_in(:)=(/  nu,  nu*dir_cart(:) /)
        lorentz(1,:)=(/ gmma,gmma*beta_vect(:) /)
        lorentz(2,:)=(/ gmma*beta_vect(1),  1.+(gmma-1)*(beta_vect(1)/beta_scalar)**2,         (gmma-1)*(beta_vect(1)*beta_vect(2)/(beta_scalar**2)),  (gmma-1)*(beta_vect(1)*beta_vect(3)/(beta_scalar**2))   /)
        lorentz(3,:)=(/ gmma*beta_vect(2),  (gmma-1)*(beta_vect(1)*beta_vect(2)/(beta_scalar**2)),  1.+(gmma-1)*(beta_vect(2)/beta_scalar)**2,         (gmma-1)*(beta_vect(2)*beta_vect(3)/(beta_scalar**2))   /)
        lorentz(4,:)=(/ gmma*beta_vect(3),  (gmma-1)*(beta_vect(1)*beta_vect(3)/(beta_scalar**2)),  (gmma-1)*(beta_vect(2)*beta_vect(3)/(beta_scalar**2)),  1.+(gmma-1)*(beta_vect(3)/beta_scalar)**2          /)

        !update values for frequency, direction of propagation and weight
        four_vect_out=MATMUL(lorentz,four_vect_in)
        nu=four_vect_out(1)
        dir_cart(:)=four_vect_out(2:4)/nu
        dir_cart=normalise(dir_cart)
        IF (str .eq. "scat") THEN
            weight=weight*(four_vect_out(1)/four_vect_in(1))
        ELSE
            weight=1
        END IF

    END SUBROUTINE

    !backwards lorentz transform to move from rest frame of observer to comoving frame
    !updates frequency, direction of propagation of packet and weight of packet
    SUBROUTINE inv_lorentz_trans(vel_vect,dir_cart,nu,weight,str)
        REAL, INTENT(IN)                ::  vel_vect(3)
        REAL, INTENT(INOUT)             ::  dir_cart(3)
        REAL, INTENT(INOUT)             ::  nu,weight
        CHARACTER(LEN=4), INTENT(IN)    ::  str

        REAL                            ::  lorentz(4,4)
        REAL                            ::  gmma
        REAL                            ::  beta_vect(3),beta_scalar
        REAL                            ::  four_vect_in(4),four_vect_out(4)

        !standard beta and gamma definitions as per lorentz transforms
        beta_vect=(1e3*vel_vect)/c
        beta_scalar=(beta_vect(1)**2+beta_vect(2)**2+beta_vect(3)**2)**0.5
        gmma=1/(1-beta_scalar**2)**0.5
        !uncomment if ignoring terms > order(v^2/c^2)
        !gmma=1.

        !construct four vector and apply inverse lorentz transform
        four_vect_in(:)=(/  nu,  nu*dir_cart(:) /)
        lorentz(1,:)=(/  gmma,-gmma*beta_vect(:) /)
        lorentz(2,:)=(/ -gmma*beta_vect(1),  1.+(gmma-1)*(beta_vect(1)/beta_scalar)**2,         (gmma-1)*(beta_vect(1)*beta_vect(2)/(beta_scalar**2)),  (gmma-1)*(beta_vect(1)*beta_vect(3)/(beta_scalar**2))   /)
        lorentz(3,:)=(/ -gmma*beta_vect(2),  (gmma-1)*(beta_vect(1)*beta_vect(2)/(beta_scalar**2)),  1.+(gmma-1)*(beta_vect(2)/beta_scalar)**2,         (gmma-1)*(beta_vect(2)*beta_vect(3)/(beta_scalar**2))   /)
        lorentz(4,:)=(/ -gmma*beta_vect(3),  (gmma-1)*(beta_vect(1)*beta_vect(3)/(beta_scalar**2)),  (gmma-1)*(beta_vect(2)*beta_vect(3)/(beta_scalar**2)),  1.+(gmma-1)*(beta_vect(3)/beta_scalar)**2          /)

        !update values for frequency, direction of propagation and weight
        four_vect_out=MATMUL(lorentz,four_vect_in)
        nu=four_vect_out(1)
        dir_cart(:)=four_vect_out(2:4)/nu
        dir_cart=normalise(dir_cart)
        IF (str .eq. "scat") THEN
            weight=weight*(four_vect_out(1)/four_vect_in(1))
        ELSE
            weight=1
        END IF

    END SUBROUTINE

END MODULE
