MODULE class_grid

    use globals
    use class_geometry
    use class_dust


    implicit none

    TYPE grid_obj
        REAL        ::  x_min,y_min,z_min            !minimum bounds of grid (cm)
        REAL        ::  x_max,y_max,z_max            !maximum bounds of grid (cm)
        INTEGER     ::  ncells(3)                    !number of cells in x/y/z directions
        INTEGER     ::  cell_width(3)                !width of grid cells in x/y/z directions
        REAL        ::  cell_vol                              !volume of grid cell
        INTEGER     ::  totcells                     !total number of cells
    END TYPE grid_obj

    TYPE(grid_obj)  ::  mothergrid                   !properties of the mothergrid

    TYPE grid_cell_obj
        REAL    ::  axis(3)                         !limits of grid cell in x, y and z
        REAL    ::  rho,nrho,N_e                    !mass density, number density and electron density of grid cell
        INTEGER ::  id(3)                           !grid cell number in each of x, y and z
        INTEGER ::  cellStatus,numPhots
    END TYPE

    TYPE(grid_cell_obj),ALLOCATABLE  :: grid_cell(:) !mothergrid is comprised of'totcells' grid_cells

    REAL        ::  tot_vol                          !total volume of supernova in 1e42cm^3

contains

    SUBROUTINE build_dust_grid()

        !open files to write gridcell coords to (in cm)
        OPEN(32,file='grid.in')

        !set number of cells in each direction and calculate total number of cells
        !mothergrid%ncells(1) read in from input file in input.f90
        !!edit here if different number of cells in each direction required
        mothergrid%ncells(2)=mothergrid%ncells(1)
        mothergrid%ncells(3)=mothergrid%ncells(1)
        mothergrid%totcells=mothergrid%ncells(1)*mothergrid%ncells(2)*mothergrid%ncells(3)

        mothergrid%cell_width(1)=(mothergrid%x_max-mothergrid%x_min)/mothergrid%ncells(1)
        mothergrid%cell_width(2)=(mothergrid%y_max-mothergrid%y_min)/mothergrid%ncells(2)
        mothergrid%cell_width(3)=(mothergrid%z_max-mothergrid%z_min)/mothergrid%ncells(3)
        mothergrid%cell_vol=((mothergrid%cell_width(1)/1e14)*(mothergrid%cell_width(2)/1e14)*(mothergrid%cell_width(3)/1e14)) !in 1e42cm^3

        SELECT CASE(dust_geometry%type)
            CASE("shell")
                !calculate R_max and R_min (dust) based on maximum velocity, day no (d=v*t) and R_min/R_max ratio
                dust_geometry%R_max=dust_geometry%v_max*day_no*8.64E-6
                dust_geometry%R_min=dust_geometry%R_ratio*dust_geometry%R_max
                IF (dust_geometry%R_min>dust_geometry%R_max) THEN
                    PRINT*, "Please specify an R_min/R_max ratio that is less than 1.  Aborted."
                    STOP
                END IF

                !convert supernova bounds from e15cm to cm
                dust_geometry%R_min_cm=dust_geometry%R_min*1e15
                dust_geometry%R_max_cm=dust_geometry%R_max*1e15

                !convert dust mass from M_sun to grams
                dust%mass_grams=dust%mass*1.98855e33

                !set bound of grid to be radius of SN
                mothergrid%x_min=-1*MAX(dust_geometry%R_max_cm,gas_geometry%R_max*1e15)
                mothergrid%y_min=mothergrid%x_min
                mothergrid%z_min=mothergrid%x_min
                mothergrid%x_max=-mothergrid%x_min
                mothergrid%y_max=mothergrid%x_max
                mothergrid%z_max=mothergrid%x_max

                !calculate total volume of shell in 1e42cm^3
                tot_vol=1000*4*pi*(dust_geometry%R_max**3-dust_geometry%R_min**3)/3 !in e42cm^3

            CASE("torus")
                PRINT*, 'You have selected a torus distribution of dust.  This routine has not been written yet.  &
                & It is due to be added to the class_grid module in due course.'
                WRITE(55,*) 'You have specified a torus distribution of dust.  Damocles is not yet capable of creating this grid.  &
                & It is due to be added to the class_grid module.  Aborted.'
            CASE("arbitrary")
                PRINT*, 'You have selected an arbitrary distribution of dust.  This routine has not been written yet.  &
                & It is due to be added to the class_grid module in due course.'
                WRITE(55,*) 'You have specified an arbitrary distribution of dust.  Damocles is not yet capable of handling arbitrary geometry grids.  &
                & This capability is due to be added to the class_grid module.  Aborted.'
            CASE("bipolar")
                PRINT*, 'You have selected a bipolar distribution of dust.  This routine has not been written yet.  &
                & It is due to be added to the class_grid module in due course.'
                WRITE(55,*) 'You have specified a bipolar distribution of dust.  Damocles is not yet capable of creating this grid.  &
                & It is due to be added to the class_grid module.  Aborted.'
            CASE DEFAULT
                PRINT*, 'Please specify one of the following options for the geometry "shell", "torus", "arbitrary" or "bipolar". Aborted.'
                WRITE(55,*) 'Aborted - dust geometry option was not one of "shell", "torus", "arbitrary" or "bipolar".'
        END SELECT


    END SUBROUTINE
    SUBROUTINE build_emissivity_dist()

        SELECT CASE(gas_geometry%type)

            CASE("shell")
                IF (lg_decoupled) THEN
                    gas_geometry%R_max=gas_geometry%v_max*day_no*8.64E-6
                    gas_geometry%R_min=gas_geometry%R_ratio*gas_geometry%R_max
                ELSE
                    gas_geometry%R_ratio=dust_geometry%R_ratio
                    gas_geometry%R_min=dust_geometry%R_min
                    gas_geometry%R_max=dust_geometry%R_max
                    gas_geometry%v_power=dust_geometry%v_power
                    gas_geometry%rho_power=dust_geometry%rho_power
                    gas_geometry%emis_power=dust_geometry%emis_power
                    gas_geometry%v_max=dust_geometry%v_max
                END IF

            CASE("torus")
                PRINT*, 'You have selected a torus distribution of gas.  This routine has not been written yet.  &
                & It is due to be added to the class_grid module in due course.'
                WRITE(55,*) 'You have specified a torus distribution of gas.  Damocles is not yet capable of creating this grid.  &
                & It is due to be added to the class_grid module.  Aborted.'
            CASE("arbitrary")
                PRINT*, 'You have selected an arbitrary distribution of gas.  This routine has not been written yet.  &
                & It is due to be added to the class_grid module in due course.'
                WRITE(55,*) 'You have specified an arbitrary distribution of gas.  Damocles is not yet capable of handling arbitrary geometry grids.  &
                & This capability is due to be added to the class_grid module.  Aborted.'
            CASE("bipolar")
                PRINT*, 'You have selected a bipolar distribution of gas.  This routine has not been written yet.  &
                & It is due to be added to the class_grid module in due course.'
                WRITE(55,*) 'You have specified a bipolar distribution of gas.  Damocles is not yet capable of creating this grid.  &
                & It is due to be added to the class_grid module.  Aborted.'
            CASE DEFAULT
                PRINT*, 'Please specify one of the following options for the geometry "shell", "torus", "arbitrary" or "bipolar". Aborted.'
                WRITE(55,*) 'Aborted - dust geometry option was not one of "shell", "torus", "arbitrary" or "bipolar".'

        END SELECT


    END SUBROUTINE build_emissivity_dist

END MODULE class_grid
