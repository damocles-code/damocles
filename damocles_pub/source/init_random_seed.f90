SUBROUTINE init_random_seed()

    !SUBROUTINE TO SET RANDOM SEED FROM CLOCK SO USE DIFFERENT RANDOM MATRIX EACH RUN
	
    INTEGER :: i, num, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed

    CALL RANDOM_SEED(size = num)
    ALLOCATE(seed(num))
          
    CALL SYSTEM_CLOCK(COUNT=clock)
    seed = clock + 37 * (/ (i - 1, i = 1, num) /)
    CALL RANDOM_SEED(PUT = seed)
          
    DEALLOCATE(seed)

END SUBROUTINE init_random_seed
