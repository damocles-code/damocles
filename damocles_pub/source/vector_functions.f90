MODULE vector_functions

    use input

    IMPLICIT NONE

contains

    FUNCTION cart(theta,phi)

        REAL    ::  cart(3)
        REAL    ::  theta,phi,x,y,z

        x=SIN(theta)*COS(phi)
        y=SIN(theta)*SIN(phi)
        z=COS(theta)
        cart=(/x,y,z/)

    END FUNCTION

    FUNCTION cartr(r,theta,phi)

        REAL    ::  cartr(3)
        REAL    ::  theta,phi,r,x,y,z

        x=r*SIN(theta)*COS(phi)
        y=r*SIN(theta)*SIN(phi)
        z=r*COS(theta)
        cartr=(/x,y,z/)

    END FUNCTION


    FUNCTION normalise(cart)
        REAL    ::  normalise(3),cart(3)
        INTEGER ::  i_dir

        DO i_dir=1,3
            normalise(i_dir)=cart(i_dir)/((cart(1)**2+cart(2)**2+cart(3)**2)**0.5)
        END DO
    END FUNCTION

    SUBROUTINE lorentz_trans(vel_vect,nu,dir_cart,w,str)
        REAL, INTENT(IN)    ::  vel_vect(3)
        REAL, INTENT(INOUT) ::  nu,dir_cart(3),w
        CHARACTER(LEN=4), INTENT(IN)   ::  str
        REAL    ::  lorentz(4,4),gmma,beta(3),sbeta,vect(4),vect2(4)

        beta=(1e3*vel_vect)/c
        sbeta=(beta(1)**2+beta(2)**2+beta(3)**2)**0.5
        gmma=1/(1-sbeta**2)**0.5
        !gmma=1.

        !THIS IS FOR ORDER O(v/c) AND IGNORING TERMS >O(v^2/c^2)
        !TREAT AS RELATIVISTIC IN FUTURE?

        vect(:)=(/ nu, nu*dir_cart(:) /)
        lorentz(1,:)=(/ gmma,   gmma*beta(:) /)
        !lorentz(2,:)=(/ gmma*beta(1),  1.,  0.,  0.   /)
        !lorentz(3,:)=(/ gmma*beta(2),  0.,  1.,  0.   /)
        !lorentz(4,:)=(/ gmma*beta(3),  0.,  0.,  1.   /)

        lorentz(2,:)=(/ gmma*beta(1),  1.+(gmma-1)*(beta(1)/sbeta)**2,         (gmma-1)*(beta(1)*beta(2)/(sbeta**2)),  (gmma-1)*(beta(1)*beta(3)/(sbeta**2))   /)
        lorentz(3,:)=(/ gmma*beta(2),  (gmma-1)*(beta(1)*beta(2)/(sbeta**2)),  1.+(gmma-1)*(beta(2)/sbeta)**2,         (gmma-1)*(beta(2)*beta(3)/(sbeta**2))   /)
        lorentz(4,:)=(/ gmma*beta(3),  (gmma-1)*(beta(1)*beta(3)/(sbeta**2)),  (gmma-1)*(beta(2)*beta(3)/(sbeta**2)),  1.+(gmma-1)*(beta(3)/sbeta)**2          /)
        


        vect2=MATMUL(lorentz,vect)
        nu=vect2(1)
        dir_cart(:)=vect2(2:4)/nu
        dir_cart=normalise(dir_cart)
        IF (str .eq. "scat") THEN
            w=w*(vect2(1)/vect(1))
        ELSE
            w=1
        END IF

    END SUBROUTINE

    SUBROUTINE inv_lorentz_trans(vel_vect,nu,dir_cart,w,str)
        REAL, INTENT(IN)    ::  vel_vect(3)
        REAL, INTENT(INOUT) ::  nu,dir_cart(3),w
        CHARACTER(LEN=4), INTENT(IN)   ::  str
        REAL    ::  lorentz(4,4),gmma,beta(3),sbeta,vect(4),vect2(4)


        beta=(1e3*vel_vect)/c
        sbeta=(beta(1)**2+beta(2)**2+beta(3)**2)**0.5
        gmma=1/(1-sbeta**2)**0.5
        !gmma=1.
        !THIS IS FOR ORDER O(v/c) AND IGNORING TERMS >O(v^2/c^2)
        !TREAT AS RELATIVISTIC IN FUTURE?

        vect(:)=(/ nu, nu*dir_cart(:) /)
        lorentz(1,:)=(/ gmma,   -gmma*beta(:) /)
        lorentz(2,:)=(/ -gmma*beta(1),  1.+(gmma-1)*(beta(1)/sbeta)**2,         (gmma-1)*(beta(1)*beta(2)/(sbeta**2)),  (gmma-1)*(beta(1)*beta(3)/(sbeta**2))   /)
        lorentz(3,:)=(/ -gmma*beta(2),  (gmma-1)*(beta(1)*beta(2)/(sbeta**2)),  1.+(gmma-1)*(beta(2)/sbeta)**2,         (gmma-1)*(beta(2)*beta(3)/(sbeta**2))   /)
        lorentz(4,:)=(/ -gmma*beta(3),  (gmma-1)*(beta(1)*beta(3)/(sbeta**2)),  (gmma-1)*(beta(2)*beta(3)/(sbeta**2)),  1.+(gmma-1)*(beta(3)/sbeta)**2          /)

        vect2=MATMUL(lorentz,vect)
        nu=vect2(1)
        dir_cart(:)=vect2(2:4)/nu
        dir_cart=normalise(dir_cart)
        IF (str .eq. "scat") THEN
            w=w*(vect2(1)/vect(1))
        ELSE
            w=1
        END IF

    END SUBROUTINE

END MODULE
