SUBROUTINE init_random_seed()

    !SUBROUTINE TO SET random SEED FROM CLOCK SO USE DIFFERENT random MATRIX EACH RUN

    INTEGER :: i, num, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed

    CALL random_SEED(size = num)
    ALLOCATE(seed(num))
          
    CALL SYSTEM_CLOCK(COUNT=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, num) /)
    CALL random_SEED(PUT = seed)
          
    DEALLOCATE(seed)

END SUBROUTINE init_random_seed
