!---------------------------------------------------------------------!
!  this subroutine sets the seed for the random number generator.     !
!  it uses the system clock to set the seed such that a different     !
!  set of random numbers is used on each run.
!---------------------------------------------------------------------!

subroutine init_random_seed()

    integer              :: i, num, clock
    integer, allocatable :: seed(:)

    call random_seed(size = num)
    allocate(seed(num))
          
    call system_clock(count=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, num) /)
    call random_seed(put = seed)

    deallocate(seed)

end subroutine init_random_seed
