MODULE input

    IMPLICIT NONE

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!DECLARE VARIABLES - AMEND HERE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    REAL                ::  lambda_0,lambda2_0,lambda1_0		!wavelength of line(s)           nm
    REAL                ::  L_tot                  				!total luminosity of line       W/um
    REAL                ::  L_Halpha                            !total Halpha luminosity (for e- scattering calcn)
    REAL                ::  L_ratio                             !total flux ratio if doublet (lambda_0/lambda2_0)
    REAL                ::  v_max                               !maximum velocity (at R_out)    km/s
    REAL                ::  x_min,y_min,z_min					!minimum bounds of grid         cm
    REAL                ::  x_max,y_max,z_max                   !maximum bounds of grid         cm
    REAL                ::  R_min,R_max					        !inner and outer radii of dust  e15cm
    REAL                ::  Rrat                                !ratio of inner to outer radii (Rin/Rout)
    REAL                ::  R_min_gas,R_max_gas                 !inner and outer radii of gas   e15cm 
    REAL                ::  R_min_cm,R_max_cm
    REAL                ::  MD_tot                              !total mass of dust             M_sun
    REAL                ::  MD_tot_grams                        !total mass of dust             grams

    REAL                ::  v_max_gas                           !maximum gas velocity for gas shell
    REAL                ::  Rrat_gas                            !ratio of R_in to Rout for gas shell

    REAL                ::  mf      							!mf = mass fraction in clumps
    REAL                ::  ff                                  !filling factor
    REAL                ::  dencon                              !dencon= density contrast
    REAL                ::  q,l,b                               !rho ~ r^-q, v ~ r^l (number of photons ~ rho^2)
    REAL                ::  q_clump                             !number density distribution of clumps
    REAL                ::  q_gas,l_gas,b_gas                   !as above but for gas profile
    REAL                ::  v_min_obs                           !minimum velocity as determined from data - 'corner of profile'
    REAL                ::  clumpgas                            !this is a 1/0 flag to determine if all packets emitted from clumps - IN TESTING

    INTEGER(8)       	::  n_packets       					!number of packets to use in simulation
    INTEGER(8)          ::  n_shells                            !number of shells to use when generating packets

    INTEGER             ::  doublet                             !0 if single line, 1 if doublet
    INTEGER             ::  ES                                  !0 is don't include e- scattering, 1 if including
    INTEGER             ::  ES_temp                             !electron scattering temperature (0 if no ES)
    INTEGER             ::  n_bins,ls                           !number of frequency bins, ls = line of sight = 0 or 1
    INTEGER             ::  VelShift                            !velshift = 0 or 1 (for repeated doppler shifting)
    INTEGER             ::  ncells(3)                           !number of cells in x/y/z directions
    INTEGER             ::  totcells		                    !total number of cells
    INTEGER             ::  dayno 
    INTEGER             ::  gas_shell                           !gas shell = 0 or 1 (if using a gas shell or not)

    INTEGER             ::  pp_no                               !no of variable pliny parameters (0 for not using pliny)
    INTEGER             ::  nargs                               !number of input arguments

    LOGICAL		    	::  lgLOS,lgVelShift,lgClumping, &
                            & lgDoublet, lgES, lgGasShell, lg_pliny       !logicals
										
    CHARACTER(LEN=25) 	::	name				                !Name of line e.g. Halpha (for output file naming)
    CHARACTER(LEN=50)   ::  inputfile                           !name of input file to read in
    CHARACTER(LEN=40)	::	dustfile                            !filename containing details of species
	
    REAL, PARAMETER 	::  pi=3.141592654
    REAL,PARAMETER      ::  c=3E8					            !in SI units (m/s)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!DECLARE OTHER FIXED VARIABLES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

    SUBROUTINE read_input(param_struct)

        REAL,INTENT(INOUT)  ::  param_struct(pp_no)

        nargs=command_argument_count()
        IF (nargs==1) THEN
           call get_command_argument(1,inputfile)
           inputfile=trim(inputfile)
           PRINT*, 'reading input from file ', inputfile
        ELSE IF (nargs==0) THEN
           inputfile='input.in'
           PRINT*,'reading input from file input.in'
        ELSE
           PRINT*,'too many input arguments - aborted'
           STOP
        END IF

        !READ IN INPUT FILE AND STORE
        OPEN(10,file=inputfile)

        READ(10,*)
        READ(10,*) pp_no
        READ(10,*)
        READ(10,*) lambda1_0
        READ(10,*) L_tot
        READ(10,*) L_Halpha
        READ(10,*) name
        READ(10,*) doublet
        READ(10,*) lambda2_0
        READ(10,*) L_ratio
        READ(10,*)
        READ(10,*) ES
        READ(10,*) ES_temp
        READ(10,*) ls
        READ(10,*) VelShift
        READ(10,*) mf
        READ(10,*) ff
        READ(10,*) clumpgas
        READ(10,*) q_clump
        READ(10,*)
        READ(10,*) dayno
        READ(10,*) v_min_obs
        READ(10,*)
        READ(10,*) v_max
        READ(10,*) Rrat
        READ(10,*) MD_tot
        READ(10,*) l
        READ(10,*) q
        READ(10,*) b
        READ(10,*)
        READ(10,*) gas_shell
        READ(10,*) v_max_gas
        READ(10,*) Rrat_gas
        READ(10,*) l_gas
        READ(10,*) q_gas
        READ(10,*) b_gas
        READ(10,*)
        READ(10,*) ncells(1)
        READ(10,*)
        READ(10,*) n_packets
        READ(10,*) n_bins
        READ(10,*) n_shells
        READ(10,*)
        READ(10,*) dustfile

        CLOSE(10)
        PRINT*,'clumpgas',clumpgas
        !READ IN DAY NO

        OPEN(21,file='line.in')
        READ(21,*)
        READ(21,*) dayno
        CLOSE(21)
        
        name=trim(name)

        !WORK OUT LOGICALS
        IF (pp_no==0) THEN
            lg_pliny=.false.
        ELSE
            lg_pliny=.true.
        END IF

        IF (ls==1) THEN
            lgLOS=.true.
        ELSE
            lgLOS=.false.
        END IF

        IF (VelShift==1) THEN
            lgVelShift=.true.
        ELSE
            lgVelShift=.false.
        END IF

        IF (mf==0) THEN
            lgClumping=.false.
        ELSE
            lgClumping=.true.
        END IF

        IF (ES==1) THEN
            lgES=.true.
            IF ((ES_temp /= 5000) .AND. (ES_temp /= 10000) .AND. (ES_temp /= 20000)) THEN
                PRINT*,'You have requested electron scattering.  Please enter a gas temperature of 5000K, 10000K or 20000K'
                STOP
            END IF
        ELSE
            lgES=.false.
        END IF

        IF (doublet==1) THEN
            lgDoublet=.true.
        ELSE
            lgDoublet=.false.
        END IF

        IF (gas_shell==1) THEN
            lgGasShell=.true.
        ELSE
            lgGasShell=.false.
        END IF

        !set number of cells in each direction and calculate total cells
        !note that vector is used to allow for different numbers in each axis in future
        ncells(2)=ncells(1)
        ncells(3)=ncells(1)
        totcells=ncells(1)*ncells(2)*ncells(3)

        IF (lg_pliny) THEN
        !set variable according to paramater structure
!       MD_tot=10**(-1*param_struct(2))
        MD_tot=param_struct(2)
        !l=param_struct(3)
        q=param_struct(4)
        END IF
       
        !calculate velocity profile from on 'corner' of
        !observed profile (i.e. vmin) and Rin/Rout ratio
        !l=(LOG10(v_min_obs/v_max))/LOG10(R_min/R_max)
        
        IF (lgGasShell) THEN
            !Rrat_gas=(v_min_obs/v_max_gas)**(1/l)
            R_max_gas=v_max_gas*dayno*8.64E-6
            R_min_gas=Rrat_gas*R_max_gas
            !l_gas=(LOG10(v_min_obs/v_max_gas))/LOG10(R_min_gas/R_max_gas)
            !1Rrat=10**((LOG10(v_min_obs/v_max))/l)
            R_min=R_min_gas
            !calculate R_max (dust) based on maximum velocity and day no
            R_max=v_max*dayno*8.64E-6
            R_min=Rrat*R_max
        ELSE
            !if gas and dust coupled then set gas and dust the same radii (temp)
            !also set same parameters
            !calculate R_max (dust) based on maximum velocity and day no
            R_max=v_max*dayno*8.64E-6
            !Rrat=(v_min_obs/v_max)**(1/l)
            Rrat_gas=Rrat
            R_min=Rrat*R_max
            R_min_gas=R_min
            R_max_gas=R_max
            l_gas=l
            q_gas=q
            b_gas=b
            v_max_gas=v_max
        END IF

        PRINT*,'minimum velocity',(R_min_gas/R_max_gas)*v_max_gas

        PRINT*,'dust vel index',l
        PRINT*,'gas vel index',l_gas
        PRINT*,'dust density index',q
        PRINT*,'gas density index',q_gas
        PRINT*,'Rmax dust',R_max
        PRINT*,'Rmin dust',R_min
        PRINT*,'Rmax gas',R_max_gas
        PRINT*,'Rmin gas',R_min_gas
        IF (R_min>R_max) THEN
           PRINT*, "R_min greater than R_max.  Abort."
           STOP
        END IF
        PRINT*,'dust r ratio',Rrat
        PRINT*,'gas r ratio',Rrat_gas

        !convert supernova bounds from e15cm to cm
        R_min_cm=R_min*1e15
        R_max_cm=R_max*1e15

        !convert dust mass from M_sun to grams
        MD_tot_grams=MD_tot*1.98855e33

        !set bound of grid to be radius of SN
        x_min=-1*MAX(R_max_cm,R_max_gas*1e15)
        y_min=x_min
        z_min=x_min

        x_max=-x_min
        y_max=x_max
        z_max=x_max
        PRINT*,'vmin',R_min*v_max/R_max
        PRINT*,'vmax',v_max
        PRINT*,'dust mass (Msun)',MD_tot
        PRINT*,'vel slope',l
        PRINT*,'rho slope',q
        PRINT*,'day no', dayno

    END SUBROUTINE read_input

END MODULE input
