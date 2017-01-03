MODULE input

    USE globals
    USE class_line
    USE class_geometry
    USE class_dust
    USE class_grid
    USE class_freq_grid

    IMPLICIT NONE

contains

    SUBROUTINE read_input()

        !!roger addition - check
        nargs=command_argument_count()
        IF (nargs==1) THEN
            call get_command_argument(1,input_file)
            input_file=trim(input_file)
            PRINT*, 'reading input from file ', input_file
        ELSE IF (nargs==0) THEN
            input_file='input.in'
            PRINT*,'reading input from file input.in'
        ELSE
            PRINT*,'too many input arguments - aborted'
            STOP
        END IF

        !OPEN LOG FILE (WILL BE CLOSED AT END OF MODEL)
        OPEN(55,file='output/log_file.out')

        !READ IN INPUT FILE AND STORE
        OPEN(10,file=input_file)

        !general options
        READ(10,*)
        READ(10,*) lg_data
        READ(10,*) data_file
        READ(10,*) lg_doublet
        READ(10,*) lg_vel_shift
        READ(10,*) lg_LOS
        READ(10,*) lg_ES
        READ(10,*) e_scat_file
        READ(10,*) day_no

        !geometry options
        READ(10,*)
        READ(10,*) lg_decoupled
        READ(10,*) dust_geometry%type
        READ(10,*) dust_file
        READ(10,*) species_file
        READ(10,*) gas_geometry%type
        READ(10,*) gas_file

        READ(10,*)
        READ(10,*) mothergrid%n_cells(1)
        READ(10,*)
        READ(10,*) n_packets
        READ(10,*) nu_grid%n_bins
        READ(10,*) n_shells
        CLOSE(10)

        !check for conflict in specified dust and gas distributions
        IF (.not. lg_decoupled) THEN
            IF (gas_geometry%type /= dust_geometry%type) THEN
                PRINT*, 'You have requested that dust and gas should be coupled but specified different geometry types. Aborted'
                STOP
            END IF
        END IF

        !read in gas options
        OPEN(11,file=gas_file)
        READ(11,*)
        READ(11,*) line%doublet_wavelength_1
        READ(11,*) line%doublet_wavelength_2
        READ(11,*) line%luminosity
        READ(11,*) line%doublet_ratio
        READ(11,*)
        READ(11,*) gas_geometry%clumped_mass_frac  !!!CURRENTLY RESTRICTED TO 0 or 1
        READ(11,*) gas_geometry%v_max
        READ(11,*) gas_geometry%R_ratio
        READ(11,*) gas_geometry%v_power
        READ(11,*) gas_geometry%rho_power
        READ(11,*) gas_geometry%emis_power
        CLOSE(11)

        !read in dust options
        OPEN(12,file=dust_file)
        READ(12,*)
        READ(12,*) dust%mass
        READ(12,*)
        READ(12,*) dust_geometry%clumped_mass_frac
        READ(12,*) dust_geometry%ff
        READ(12,*) dust_geometry%clump_power
        READ(12,*)
        READ(12,*) dust_geometry%v_max
        READ(12,*) dust_geometry%R_ratio
        READ(12,*) dust_geometry%v_power
        READ(12,*) dust_geometry%rho_power
        READ(12,*) dust_geometry%emis_power
        CLOSE(12)


        SELECT CASE(dust_geometry%type)
            CASE ('shell')
                call check_dust_clumped()
        END SELECT
        !read in electron scattering options (if using electron scattering)
        IF (lg_ES) THEN
            OPEN(13,file = e_scat_file)
            READ(13,*)
            READ(13,*) L_Halpha
            READ(13,*) ES_temp
            CLOSE(13)
            IF ((ES_temp /= 5000) .AND. (ES_temp /= 10000) .AND. (ES_temp /= 20000)) THEN
                PRINT*,'You have requested electron scattering.  Please enter a gas temperature of 5000K, 10000K or 20000K'
                STOP
            END IF
        END IF

    END SUBROUTINE read_input

END MODULE input

