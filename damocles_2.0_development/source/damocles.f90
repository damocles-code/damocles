
PROGRAM DAMOCLES

    USE input
    USE initialise
    USE vector_functions
    USE init_packet
    USE driver

    implicit none

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

    call run_code()

END PROGRAM
