
PROGRAM DAMOCLES

    USE input
    USE initialise
    USE vector_functions
    USE init_packet
    USE driver
	USE wrap
	
	implicit none

    REAL(8),ALLOCATABLE     ::  param_struct(:)
    REAL(8)                 ::  chi2
    INTEGER                 ::  p_no
    CHARACTER(LEN=50)       ::  infile

    nargs=command_argument_count()
    IF (nargs==1) THEN
        call get_command_argument(1,infile)
        infile=trim(infile)
    ELSE IF (nargs==0) THEN
        infile='input.in'
    ELSE
        PRINT*,'too many input arguments - aborted'
        STOP
    END IF

    OPEN(11,file=infile)
    READ(11,*)
    READ(11,*) pp_no
    CLOSE(11)

    IF (pp_no==0) THEN
        ALLOCATE(param_struct(1))
        param_struct(1)=0
    ELSE
        ALLOCATE(param_struct(pp_no))
        DO p_no=1,pp_no
            param_struct(i)=0
        END DO
    END IF

    call damocles_wrap(param_struct,chi2)
    PRINT*,'chi',chi2

END PROGRAM
