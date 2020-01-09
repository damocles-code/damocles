module rnglib2

    use globals

    implicit none



    contains

    function r4_uni_01()
    real(kind = 4) :: r4_uni_01
    print*, 'tehi is a test here!',lg_mcmc
    end function

!*****************************************************************************80
!
!! R4_UNI_01 returns a uniform random real number in [0,1].
!
!  Discussion:
!
!    This procedure returns a random floating point number from a uniform
!    distribution over (0,1), not including the endpoint values, using the
!    current random number generator.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    05 August 2013
!
!  Author:
!
!    Original Pascal version by Pierre L'Ecuyer, Serge Cote.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Pierre LEcuyer, Serge Cote,
!    Implementing a Random Number Package with Splitting Facilities,
!    ACM Transactions on Mathematical Software,
!    Volume 17, Number 1, March 1991, pages 98-111.
!
!  Parameters:
!
!    Output, real ( kind = 4 ) R4_UNI_01, a uniform random value in [0,1].
!
! integer ( kind = 4 ) i
! !integer ( kind = 4 ) i4_uni
! !logical initialized_get
! !real ( kind = 4 ) r4_uni_01
!
!  Check whether the package must be initialized.
!
! if ( .not. initialized_get ( ) ) then
!   !write ( *, '(a)' ) ' '
!   !write ( *, '(a)' ) 'R4_UNI_01 - Note:'
!   !write ( *, '(a)' ) '  Initializing RNGLIB package.'
!   call initialize ( )
! end if
!
!  Get a random positive integer.
!
! i = i4_uni ( )
!
!  Scale it to a random real in [0,1].
!
! r4_uni_01 = real ( i, kind = 4 ) * 4.656613057E-10
! return
!nd

end module rnglib2
