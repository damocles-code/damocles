!---------------------------------------------------------------------!
!  this module includes a number of functions to handle vectors       !
!  including - conversion from spherical coords to cartesian          !
!            - normalise a vector                                     !
!            - forwards and backwards lorentz transforms              !
!---------------------------------------------------------------------!

module vector_functions

    use input

    implicit none

contains

    !function to convert spherical direction to cartesian unit vector
    function cart(theta,phi)
        real    ::  cart(3)
        real    ::  theta,phi,x,y,z
        x=sin(theta)*cos(phi)
        y=sin(theta)*sin(phi)
        z=cos(theta)
        cart=(/x,y,z/)
    end function

    !function to convert spherical corrdinates to cartesian coordinates
    function cartr(r,theta,phi)
        real    ::  cartr(3)
        real    ::  theta,phi,r,x,y,z
        x=r*sin(theta)*cos(phi)
        y=r*sin(theta)*sin(phi)
        z=r*cos(theta)
        cartr=(/x,y,z/)
    end function

    !normalise a cartesian vector to a unit vector
    function normalise(cart)
        real    ::  normalise(3),cart(3)
        integer ::  i_dir
        do i_dir=1,3
            normalise(i_dir)=cart(i_dir)/((cart(1)**2+cart(2)**2+cart(3)**2)**0.5)
        end do
    end function

    !forwards lorentz transform to move from comoving frame to rest frame of observer
    !updates frequency, direction of propagation of packet and weight of packet
    subroutine lorentz_trans(vel_vect,dir_cart,nu,weight,str)
        real, intent(in)                ::  vel_vect(3)
        real, intent(inout)             ::  dir_cart(3)
        real, intent(inout)             ::  nu,weight
        character(len=4), intent(in)    ::  str

        real                            ::  lorentz(4,4)
        real                            ::  gmma
        real                            ::  beta_vect(3),beta_scalar
        real                            ::  four_vect_in(4),four_vect_out(4)

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
        four_vect_out=matmul(lorentz,four_vect_in)
        nu=four_vect_out(1)
        dir_cart(:)=four_vect_out(2:4)/nu
        dir_cart=normalise(dir_cart)
        if (str .eq. "scat") then
           weight=weight*(four_vect_out(1)/four_vect_in(1))
        else
           weight=1
        end if


    end subroutine

    !backwards lorentz transform to move from rest frame of observer to comoving frame
    !updates frequency, direction of propagation of packet and weight of packet
    subroutine inv_lorentz_trans(vel_vect,dir_cart,nu,weight,str)
        real, intent(in)                ::  vel_vect(3)
        real, intent(inout)             ::  dir_cart(3)
        real, intent(inout)             ::  nu,weight
        character(len=4), intent(in)    ::  str

        real                            ::  lorentz(4,4)
        real                            ::  gmma
        real                            ::  beta_vect(3),beta_scalar
        real                            ::  four_vect_in(4),four_vect_out(4)

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
        four_vect_out=matmul(lorentz,four_vect_in)
        nu = four_vect_out(1)
        dir_cart(:)=four_vect_out(2:4)/nu
        dir_cart=normalise(dir_cart)
        if (str .eq. "scat") then
           weight=weight*(four_vect_out(1)/four_vect_in(1))
        else
           weight=1
        end if


    end subroutine

    !backwards lorentz transform to move from rest frame of observer to comoving frame
    !updates frequency, direction of propagation of packet and weight of packet
    subroutine inv_lorentz_trans_for_ext(vel_vect,dir_cart,nu,nu_updated)
        real, intent(in)                ::  vel_vect(3)
        real, intent(inout)                ::  dir_cart(3)
        real, intent(inout)                ::  nu
        real, intent(out)               ::  nu_updated

        real                            ::  lorentz(4,4)
        real                            ::  gmma
        real                            ::  beta_vect(3),beta_scalar
        real                            ::  four_vect_in(4),four_vect_out(4)

        !standard beta and gamma definitions as per lorentz transforms
        beta_vect=(1e3*vel_vect)/c
        beta_scalar=(beta_vect(1)**2+beta_vect(2)**2+beta_vect(3)**2)**0.5
        gmma=1/(1-beta_scalar**2)**0.5
        !uncomment if ignoring terms > order(v^2/c^2)
        !gmma=1.

        !construct four vector and apply inverse lorentz transform
        four_vect_in(:)=(/  nu,  nu*dir_cart(:) /)
        lorentz(1,:)=(/  gmma,-gmma*beta_vect(:) /)
        lorentz(2,:)=(/ -gmma*beta_vect(1), 1.+(gmma-1)*(beta_vect(1)/beta_scalar)**2,(gmma-1)*(beta_vect(1)*beta_vect(2)/(beta_scalar**2)),(gmma-1)*(beta_vect(1)*beta_vect(3)/(beta_scalar**2)) /)
        lorentz(3,:)=(/ -gmma*beta_vect(2),(gmma-1)*(beta_vect(1)*beta_vect(2)/(beta_scalar**2)),1.+(gmma-1)*(beta_vect(2)/beta_scalar)**2,(gmma-1)*(beta_vect(2)*beta_vect(3)/(beta_scalar**2)) /)
        lorentz(4,:)=(/ -gmma*beta_vect(3),(gmma-1)*(beta_vect(1)*beta_vect(3)/(beta_scalar**2)),(gmma-1)*(beta_vect(2)*beta_vect(3)/(beta_scalar**2)),1.+(gmma-1)*(beta_vect(3)/beta_scalar)**2 /)

        !update values for frequency, direction of propagation and weight
        four_vect_out=matmul(lorentz,four_vect_in)
        nu_updated=four_vect_out(1)

    end subroutine

end module
