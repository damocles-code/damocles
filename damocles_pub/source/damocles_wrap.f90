MODULE wrap

    USE input
    USE initialise
    USE vector_functions
    USE init_packet
    USE driver

contains

    SUBROUTINE damocles_wrap(param_struct,chi2)
            REAL(8)       ::  param_struct(pp_no),chi2
            REAL          ::  params(pp_no)
            REAL          ::  chi2_single

            chi2_single=REAL(chi2)
            PRINT*,'PARAMETERS - before conversion - ',param_struct
            params=REAL(param_struct)
            call run_code(params,chi2_single)
            chi2=DBLE(chi2_single)
            PRINT*,'---------------------------------------------'
            
	RETURN 
    END SUBROUTINE damocles_wrap

END MODULE
