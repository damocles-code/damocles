!---------------------------------------------------------------------!
!  this subroutine sets the seed for the random number generator.     !
!  it uses the system clock to set the seed such that a different     !
!  set of random numbers is used on each run.
!---------------------------------------------------------------------!

SUBROUTINE init_random_seed()

    INTEGER              :: i, num, clock
    INTEGER, ALLOCATABLE :: seed(:)

    CALL random_SEED(size = num)
    ALLOCATE(seed(num))
          
    CALL SYSTEM_CLOCK(COUNT=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, num) /)
    CALL random_SEED(PUT = seed)
          
    DEALLOCATE(seed)

END SUBROUTINE init_random_seed
