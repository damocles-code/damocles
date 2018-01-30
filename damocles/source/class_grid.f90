!---------------------------------------------------------------------------------------------------!
!  the 'grid' and 'grid cell' types are declared in this module.                                    !
!  the subroutines contained in this module construct grids to be used in the simulation.           !
!  cartesian grids of dust densities are either constructed or read in.                             !
!  emissivity distributions are also generated.                                                     !
!---------------------------------------------------------------------------------------------------!

module class_grid

    use globals
    use class_geometry
    use class_dust
    use class_freq_grid
    use rnglib

    implicit none

    real,allocatable        ::  profile_array(:)          !array to store the resultant line profile
    real,allocatable        ::  profile_los_array(:,:,:)  !array to store resultant profiles divided into lines of sight
    real,allocatable        ::  profile_array_data_bins(:)!array to store the resultant line profile binned into the same frequencies as the observed line
    real,allocatable        ::  square_weight_data_bins(:)!array to store the sums of the square of the weights in each bin
    real,allocatable        ::  n_packets_data_bins(:)    !array to store the number of packets in each bin
    real,allocatable        ::  total_weight_data_bins(:) !array to store the total packet weight in each bin
    real,allocatable        ::  mc_error_data_bins(:)     !array to store the monte carlo error in each bin
    real                    ::  delta                     !dummy variable used in calculation the mean and sigma of the packet weights in each bin
    real,allocatable        ::  shell_radius(:,:)         !radial bounds of each shell in 1e15cm (if using shell geometry)
    real,allocatable        ::  frac_shell_radius(:,:)    !radial bounds of each shell as a fraction of the whole scaled to  r_max = 1 and r_min = r_ratio
    integer(8),allocatable  ::  num_packets_array(:,:)      !number of packets to be emitted in each cell and cumulative total 

    type grid_obj
        integer          ::  n_cells(3)                   !number of cells in x/y/z directions
        integer          ::  tot_cells                    !total number of cells
        real             ::  n_rho_dust_av
        real             ::  cell_vol                     !volume of grid cell (only applicable for cubic grids with equal no of divisions in each axis)
        real             ::  x_min,x_max                  !outer bounds of grid in x dimension
        real             ::  y_min,y_max                  !outer bounds of grid in y dimension
        real             ::  z_min,z_max                  !outer bounds of grid in z dimension
        real             ::  cell_width(3)                !width of grid cells in x/y/z directions (only applicable for cubic grids with equal no of divisions in each axis)
        real,allocatable ::  x_div(:),y_div(:),z_div(:)   !locations of divisions listed consecutively in x, y and z
    end type grid_obj

    type(grid_obj)       ::  mothergrid                   !properties of the mothergrid

    type grid_cell_obj
        real             ::  rho                          !mass density of grid cell
        real             ::  n_rho                        !number density of grid cell
        real             ::  n_e                          !electron density of grid cell
        real             ::  r                            !radial distance from (0,0,0) to the centre of the cell
        real             ::  vol                          !volume of cell
        real             ::  axis(3)                      !limits of grid cell in x, y and z
        real             ::  width(3)                     !width of grid cells in x/y/z directions
        logical          ::  lg_clump                     !logical indicating whether or not grid cell is a clump
        integer          ::  id(3)                        !grid cell number in each of x, y and z
    end type

    type(grid_cell_obj),allocatable :: grid_cell(:)       !mothergrid is comprised of'tot_cells' grid_cells


contains

    subroutine build_dust_grid()

        if (.not. lg_mcmc) print*, 'constructing grid...'

        select case(dust_geometry%type)

            case("shell")

                !calculate r_max and r_min (dust) based on maximum velocity, day no (d=v*t) and r_min/r_max ratio
                dust_geometry%r_max=dust_geometry%v_max*day_no*8.64e-6
                dust_geometry%r_min=dust_geometry%r_ratio*dust_geometry%r_max

                !check rin/rout ratio
                if (dust_geometry%r_min>dust_geometry%r_max) then
                    print*, "please specify an r_min/r_max ratio that is less than 1.  aborted."
                    stop
                end if

                !convert supernova bounds from e15cm to cm
                dust_geometry%r_min_cm=dust_geometry%r_min*1e15
                dust_geometry%r_max_cm=dust_geometry%r_max*1e15

                !convert dust mass from m_sun to grams
                dust%mass_grams=dust%mass*1.98855e33

                !set bounds of grid to be radius of sn
                mothergrid%x_min=-1*max(dust_geometry%r_max_cm,gas_geometry%r_max*1e15)
                mothergrid%y_min=mothergrid%x_min
                mothergrid%z_min=mothergrid%x_min
                mothergrid%x_max=-mothergrid%x_min
                mothergrid%y_max=mothergrid%x_max
                mothergrid%z_max=mothergrid%x_max

                !set number of cells in each direction and calculate total number of cells
                !mothergrid%n_cells(1) read in from input file in input.f90
                !!edit here if different number of cells in each direction required
                mothergrid%n_cells(2)=mothergrid%n_cells(1)
                mothergrid%n_cells(3)=mothergrid%n_cells(1)
                mothergrid%tot_cells=mothergrid%n_cells(1)*mothergrid%n_cells(2)*mothergrid%n_cells(3)

                mothergrid%cell_width(1)=(mothergrid%x_max-mothergrid%x_min)/mothergrid%n_cells(1)
                mothergrid%cell_width(2)=(mothergrid%y_max-mothergrid%y_min)/mothergrid%n_cells(2)
                mothergrid%cell_width(3)=(mothergrid%z_max-mothergrid%z_min)/mothergrid%n_cells(3)
                mothergrid%cell_vol=((mothergrid%cell_width(1)/1e14)*(mothergrid%cell_width(2)/1e14)*(mothergrid%cell_width(3)/1e14)) !in 1e42cm^3

                !calculate total volume of shell in 1e42cm^3
                tot_vol=1000*4*pi*(dust_geometry%r_max**3-dust_geometry%r_min**3)/3 !in e42cm^3

                allocate(grid_cell(mothergrid%tot_cells))
                allocate(mothergrid%x_div(mothergrid%n_cells(1)))
                allocate(mothergrid%y_div(mothergrid%n_cells(2)))
                allocate(mothergrid%z_div(mothergrid%n_cells(3)))

                !for a cubic grid (as for a shell), set the cell widths and volumes to be identical for all cells
                grid_cell%width(1)=mothergrid%cell_width(1)
                grid_cell%width(2)=mothergrid%cell_width(2)
                grid_cell%width(3)=mothergrid%cell_width(3)
                grid_cell%vol=mothergrid%cell_vol

                !initialise mothergrid divisions to zero
                mothergrid%x_div=0
                mothergrid%y_div=0
                mothergrid%z_div=0

                !calculate divisions in mothergrid in each axis (grid cell coordinates)
                do ixx=1,mothergrid%n_cells(1)
                    mothergrid%x_div(ixx)=mothergrid%x_min+((ixx-1)*mothergrid%cell_width(1))
                end do

                do iyy=1,mothergrid%n_cells(2)
                    mothergrid%y_div(iyy)=mothergrid%y_min+((iyy-1)*mothergrid%cell_width(2))
                end do

                do izz=1,mothergrid%n_cells(3)
                    mothergrid%z_div(izz)=mothergrid%z_min+((izz-1)*mothergrid%cell_width(3))
                end do

                !calculate radius of the centre of each cell
                ig=0
                grid_cell%r=0
                do ixx=1,mothergrid%n_cells(1)
                    do iyy=1,mothergrid%n_cells(2)
                        do izz=1,mothergrid%n_cells(3)
                            ig=ig+1
                            grid_cell(ig)%r=((mothergrid%x_div(ixx)+mothergrid%cell_width(1)/2)**2+(mothergrid%y_div(iyy)+mothergrid%cell_width(2)/2)**2+(mothergrid%z_div(izz)+mothergrid%cell_width(3)/2)**2)**0.5
                            grid_cell(ig)%axis(:)=(/ mothergrid%x_div(ixx),mothergrid%y_div(iyy),mothergrid%z_div(izz)/)
                            grid_cell(ig)%id(:)=(/ ixx,iyy,izz /)
                        end do
                    end do
                end do

                !calculate desired mass of dust in clumps and in interclump medium (icm)
                dust%m_icm=dust%mass*(1-dust_geometry%clumped_mass_frac)

                !calculate number of clumps based on filling fraction of total volume and volume of clump (=volume of grid cell)
                dust_geometry%n_clumps=floor(dust_geometry%ff*tot_vol/mothergrid%cell_vol)

                !calculate mass of a clump based on dust mass fraction in clumps and total number of clumps
                dust%m_clump=(dust%mass*dust_geometry%clumped_mass_frac)/dust_geometry%n_clumps

                !calculate density of a clump based on mass of clump and volume of clump
                !5.02765e8=1e42/1.989e33 i.e. conversion factor for  e42cm3 to cm3 and g to msun)
                dust_geometry%rho_clump=dust%m_clump/(mothergrid%cell_vol*5.02765e8)

                !initialise all grid cells not to be clumps
                grid_cell(:)%lg_clump=.false.

                !generate grid for a clumped geometry
                if (dust_geometry%lg_clumped) then

                    !allocate certain cells to be clumps according to clump density distribution
                    !repeat until required number of clumps reached
                    n_clumps=0
                    m_clumps_check=0
                    do while (n_clumps<(dust_geometry%n_clumps))

                        !select a random cell from the grid
                        !call random_number(ran)
                        ran = r4_uni_01()

                        ig=ceiling(mothergrid%tot_cells*ran)
                        if (ig==0) cycle

                        !test cell against clump probability distribution
                        if ( (grid_cell(ig)%r<(dust_geometry%r_max_cm)) .and. &
                            & (grid_cell(ig)%r>(dust_geometry%r_min_cm)) .and. &
                            & (.not. grid_cell(ig)%lg_clump)) then

                            !call random_number(ran)
                            ran = r4_uni_01()
                            if (ran<((dust_geometry%r_min_cm/grid_cell(ig)%r)**dust_geometry%clump_power)) then
                                grid_cell(ig)%lg_clump=.true.
                                n_clumps=n_clumps+1
                                grid_cell(ig)%rho=dust_geometry%rho_clump
                                grid_cell(ig)%n_rho=grid_cell(ig)%rho/dust%av_mgrain
                                m_clumps_check=m_clumps_check+grid_cell(ig)%rho*mothergrid%cell_vol*5.02765e8
                            end if
                        end if
                    end do

                    !integrate dust density distribution over cells in interclump medium (i.e. not clump cells)
                    !calculate normalisation factor to be used in smooth density calculation later
                    norm=0
                    do ig = 1,mothergrid%tot_cells
                        if ((grid_cell(ig)%r<(dust_geometry%r_max_cm)) .and. (grid_cell(ig)%r>(dust_geometry%r_min_cm))) then
                            if (.not. grid_cell(ig)%lg_clump) norm=norm+(dust_geometry%r_min_cm/grid_cell(ig)%r)**dust_geometry%rho_power
                        end if
                    end do

                    !calculate dust density of cells in the interclump medium and check total mass correct
                    !m and micm will be used to check that the dust masses used are equal to those requested
                    ig=0
                    no_active_cells=0
                    m_tot_check=0.
                    m_icm_check=0.
                    mothergrid%n_rho_dust_av=0.
                    do ixx=1,mothergrid%n_cells(1)
                        do iyy=1,mothergrid%n_cells(2)
                            do izz=1,mothergrid%n_cells(3)
                                ig=ig+1

                                !test cell inside bounds of shell
                                if ((grid_cell(ig)%r<(dust_geometry%r_max_cm)) .and. &
                                    &(grid_cell(ig)%r>(dust_geometry%r_min_cm))) then

                                    no_active_cells=no_active_cells+1

                                    !calculate density if not a clump
                                    if (.not. grid_cell(ig)%lg_clump) then
                                        grid_cell(ig)%rho=dust%m_icm*((dust_geometry%r_min_cm/grid_cell(ig)%r)**dust_geometry%rho_power)/(5.02765e8*norm*mothergrid%cell_vol)
                                        grid_cell(ig)%n_rho=grid_cell(ig)%rho/dust%av_mgrain
                                        m_icm_check=m_icm_check+grid_cell(ig)%rho*mothergrid%cell_vol*5.02765e8
                                    end if

                                    !add mass of cell to total mass for checking purposes
                                    m_tot_check=m_tot_check+grid_cell(ig)%rho*mothergrid%cell_vol*5.02765e8

                                    !average dust number density will be divided by total number of active cells at end of loop
                                    mothergrid%n_rho_dust_av=mothergrid%n_rho_dust_av+grid_cell(ig)%n_rho

                                else
                                    !set all densities etc. to zero if cell is outside of bounds of shell
                                    grid_cell(ig)%axis(:)=(/ mothergrid%x_div(ixx),mothergrid%y_div(iyy),mothergrid%z_div(izz)/)
                                    grid_cell(ig)%rho=0.
                                    grid_cell(ig)%n_rho=0.
                                end if
                            end do
                        end do
                    end do

                else
                    !generate grid for a non-clumped, smooth geometry

                    !calculate dust density at inner radius (rho_in)
                    !factor of 1.989e-12 (= 1.989e33/1e45) to convert from msun/e45cm3 to g/cm3
                    if (dust_geometry%rho_power==3) then
                        dust_geometry%rho_in=(dust_geometry%r_min**(-dust_geometry%rho_power))*((dust%mass*1.989e33)/(log(dust_geometry%r_max/dust_geometry%r_min)*4*pi))
                        dust_geometry%rho_in=dust_geometry%rho_in*(1.989e-12)
                    else
                        dust_geometry%rho_in = ((dust%mass*(3-dust_geometry%rho_power)) &
                            & /(4*pi*(dust_geometry%r_ratio**(dust_geometry%rho_power)*(dust_geometry%r_max**3)-dust_geometry%r_min**3)))
                        dust_geometry%rho_in=dust_geometry%rho_in*(1.989e-12)
                    end if

                    !calculate dust densities for each cell in grid
                    no_active_cells=0
                    ig=0
                    do ixx=1,mothergrid%n_cells(1)
                        do iyy=1,mothergrid%n_cells(2)
                            do izz=1,mothergrid%n_cells(3)

                                ig=ig+1
                                if ((grid_cell(ig)%r<(dust_geometry%r_max_cm)) .and. (grid_cell(ig)%r>(dust_geometry%r_min_cm))) then
                                    !calculate densities for cells inside shell
                                    !grid_cell(ig)%rho=((dust_geometry%rho_in)*(dust_geometry%r_min_cm/grid_cell(ig)%r)**dust_geometry%rho_power)
                                    grid_cell(ig)%rho=((dust_geometry%rho_in)*(dust_geometry%r_min_cm/grid_cell(ig)%r)**dust_geometry%rho_power)
                                    grid_cell(ig)%n_rho=grid_cell(ig)%rho/dust%av_mgrain
                                    mothergrid%n_rho_dust_av=mothergrid%n_rho_dust_av+grid_cell(ig)%n_rho

                                    !check total mass used is correct and increment total number of active cells in shell
                                    !5.02765e8=1e42/1.989e33 i.e. conversion factor for  e42cm3 to cm3 and g to msun
                                    no_active_cells=no_active_cells+1
                                    m_tot_check=m_tot_check+grid_cell(ig)%rho*mothergrid%cell_vol*5.02765e8
                                else
                                    !set desntiies for cells outside of shell to zero
                                    grid_cell(ig)%rho=0.
                                    grid_cell(ig)%n_rho=0.
                                    grid_cell(ii)%n_e=0.
                                end if
                            end do
                        end do
                    end do

                end if

                !calculate true average dust number density by dividing by total number of active cells in shell
                mothergrid%n_rho_dust_av=mothergrid%n_rho_dust_av/no_active_cells

            case("torus")
                print*, 'you have selected a torus distribution of dust.  this routine has not been written yet.  &
                & it is due to be added to the class_grid module in due course.'
                write(55,*) 'you have specified a torus distribution of dust.  damocles is not yet capable of creating this grid.  &
                & it is due to be added to the class_grid module.  aborted.'
                stop
            case("arbitrary")
                !read in dust grid file
                !!note that currently the filename is hard coded and should be changed to be variable
                open(33,file='input/dust_grid.in')

                !the file format is based on the grid generated by the mocassin grid generator which can be found at
                !http://www.nebulousresearch.org/codes/mocassin/mocassin_gridmaker.php
                !additional first line declaring maximum velocity and maximum radius (at which maximum velocity occurs)
                !read in maximum velocity and radius at which maximum velocity occurs
                read(33,*) dust_geometry%v_max,dust_geometry%r_max,dust_geometry%v_power
                read(33,*) junk,mothergrid%n_cells(1),mothergrid%n_cells(2),mothergrid%n_cells(3)
                mothergrid%tot_cells=mothergrid%n_cells(1)*mothergrid%n_cells(2)*mothergrid%n_cells(3)

                allocate(grid_cell(mothergrid%tot_cells))
                allocate(mothergrid%x_div(mothergrid%n_cells(1)))
                allocate(mothergrid%y_div(mothergrid%n_cells(2)))
                allocate(mothergrid%z_div(mothergrid%n_cells(3)))


                xx=0
                yy=0
                zz=0
                !grid must be arranged in ascending z, then y, then x (as per mocassin grid maker)
                !!include a sort here to ensure this is the case?
                do ig=1,mothergrid%tot_cells

                    !read in limits of grid cells in x, y and z
                    read(33,*) grid_cell(ig)%axis(1),grid_cell(ig)%axis(2),grid_cell(ig)%axis(3),grid_cell(ig)%n_rho

                    !this section calculates the cell no in each axis given the overall cell number within the grid
                    !(effectively reverse engineers ig=(n_cells(2)*n_cells(3)*(id(1)-1)+n_cells(3)*(id(2)-1)+id(3)))
                    if (mod(ig,mothergrid%n_cells(3)) /=1) then
                        grid_cell(ig)%id(3)=grid_cell(ig-1)%id(3)+1
                    else
                        grid_cell(ig)%id(3)=1
                    end if
                    grid_cell(ig)%id(2)=mod((ig-grid_cell(ig)%id(3))/mothergrid%n_cells(3),mothergrid%n_cells(2))+1
                    grid_cell(ig)%id(1)=((ig-grid_cell(ig)%id(3)-(mothergrid%n_cells(3)*(grid_cell(ig)%id(2)-1)))/(mothergrid%n_cells(3)*mothergrid%n_cells(2)))+1

                    !this section populates the list of divisions in x, y and z from the grid specified
                    !each time a cell is encountered where the id of the other two axes is 1, the third non-zero axis is stored
                    !this generates a unique list of axis divisions
                    if (grid_cell(ig)%id(1) == 1 .and. grid_cell(ig)%id(2) == 1) then
                        zz=zz+1
                        mothergrid%z_div(zz)=grid_cell(ig)%axis(3)
                    end if
                    if (grid_cell(ig)%id(1) == 1 .and. grid_cell(ig)%id(3) == 1) then
                        yy=yy+1
                        mothergrid%y_div(yy)=grid_cell(ig)%axis(2)
                    end if
                    if (grid_cell(ig)%id(2) == 1 .and. grid_cell(ig)%id(3) == 1) then
                        xx=xx+1
                        mothergrid%x_div(xx)=grid_cell(ig)%axis(1)
                    end if
                end do

                close(33)

                !radii are calculated for each grid cell from the centre of the grid (0,0,0) to the centre of the cell
                !volumes are calculated by multiplying distance to next division in all axes
                !at far 'right' edge, the volumes, widths and radii are calculated using symmetries of grid
                !i.e. radius/volume at cell (43,47,50) is same as for (43,47,1)
                do ig=1,mothergrid%tot_cells
                    if (grid_cell(ig)%id(1) /= mothergrid%n_cells(1) .and. &
                        & grid_cell(ig)%id(2) /= mothergrid%n_cells(2) .and. &
                        & grid_cell(ig)%id(3) /= mothergrid%n_cells(3)) then

                        grid_cell(ig)%r=((((grid_cell(ig)%axis(1)+grid_cell(ig+1)%axis(1))/2)**2) + &
                            & (((grid_cell(ig)%axis(2)+grid_cell(ig+1)%axis(2))/2)**2) + &
                            & (((grid_cell(ig)%axis(3)+grid_cell(ig+1)%axis(3))/2)**2))**0.5
                        grid_cell(ig)%vol=(mothergrid%x_div(grid_cell(ig)%id(1)+1)-mothergrid%x_div(grid_cell(ig)%id(1)))*1e-14* &
                            & (mothergrid%y_div(grid_cell(ig)%id(2)+1)-mothergrid%y_div(grid_cell(ig)%id(2)))*1e-14* &
                            & (mothergrid%z_div(grid_cell(ig)%id(3)+1)-mothergrid%z_div(grid_cell(ig)%id(3)))*1e-14
                        grid_cell(ig)%width(1)=mothergrid%x_div(grid_cell(ig)%id(1)+1)-mothergrid%x_div(grid_cell(ig)%id(1))
                        grid_cell(ig)%width(2)=mothergrid%x_div(grid_cell(ig)%id(2)+1)-mothergrid%x_div(grid_cell(ig)%id(2))
                        grid_cell(ig)%width(3)=mothergrid%x_div(grid_cell(ig)%id(3)+1)-mothergrid%x_div(grid_cell(ig)%id(3))
                    else if (grid_cell(ig)%id(3) == mothergrid%n_cells(3)) then
                        grid_cell(ig)%vol=grid_cell(ig-mothergrid%n_cells(3)+1)%vol
                        grid_cell(ig)%r=grid_cell(ig-mothergrid%n_cells(3)+1)%r
                        grid_cell(ig)%width=grid_cell(ig-mothergrid%n_cells(3)+1)%width
                    else if (grid_cell(ig)%id(2) == mothergrid%n_cells(2)) then
                        grid_cell(ig)%vol=grid_cell(ig-mothergrid%n_cells(2)*(grid_cell(ig)%id(2)-1))%vol
                        grid_cell(ig)%r=grid_cell(ig-mothergrid%n_cells(2)*(grid_cell(ig)%id(2)-1))%r
                        grid_cell(ig)%width=grid_cell(ig-mothergrid%n_cells(2)*(grid_cell(ig)%id(2)-1))%width
                    else if (grid_cell(ig)%id(1) == mothergrid%n_cells(1)) then
                        grid_cell(ig)%vol=grid_cell(ig-mothergrid%n_cells(2)*mothergrid%n_cells(1)*(grid_cell(ig)%id(1)-1))%vol
                        grid_cell(ig)%r=grid_cell(ig-mothergrid%n_cells(2)*mothergrid%n_cells(1)*(grid_cell(ig)%id(1)-1))%r
                        grid_cell(ig)%width=grid_cell(ig-mothergrid%n_cells(2)*mothergrid%n_cells(1)*(grid_cell(ig)%id(1)-1))%width
                    end if
                end do

            case("bipolar")
                print*, 'you have selected a bipolar distribution of dust.  this routine has not been written yet.  &
                & it is due to be added to the class_grid module in due course.'
                write(55,*) 'you have specified a bipolar distribution of dust.  damocles is not yet capable of creating this grid.  &
                & it is due to be added to the class_grid module.  aborted.'
                stop
            case default
                print*, 'please specify one of the following options for the geometry "shell", "torus", "arbitrary" or "bipolar". aborted.'
                write(55,*) 'aborted - dust geometry option was not one of "shell", "torus", "arbitrary" or "bipolar".'
                stop
        end select

    end subroutine

    subroutine build_emissivity_dist()

        if (.not. lg_mcmc) print*,'building emissivity distribution...'

        allocate(profile_array(nu_grid%n_bins))
        allocate(profile_los_array(nu_grid%n_bins,n_angle_divs,n_angle_divs))
        allocate(profile_array_data_bins(obs_data%n_data))
        allocate(n_packets_data_bins(obs_data%n_data))
        allocate(square_weight_data_bins(obs_data%n_data))
        allocate(total_weight_data_bins(obs_data%n_data))
        allocate(mc_error_data_bins(obs_data%n_data))
        profile_array=0
        profile_los_array=0
        profile_array_data_bins=0
        square_weight_data_bins=0
        total_weight_data_bins=0

        select case(gas_geometry%type)

            case("shell")
                !assign parameters related to gas geometry
                if (lg_decoupled) then
                    !if decoupled then generate max/min radii from epoch and maximum velocity
                    gas_geometry%r_max=gas_geometry%v_max*day_no*8.64e-6
                    gas_geometry%r_min=gas_geometry%r_ratio*gas_geometry%r_max
                else
                    !if decoupled then set gas geometry parameters to equal the dust geometry parameters
                    if (gas_geometry%type /= dust_geometry%type) then
                        print*, 'you have requested that gas and dust distributions be coupled but specified different geometries.  aborted.'
                        stop
                    end if
                    gas_geometry%r_ratio=dust_geometry%r_ratio
                    gas_geometry%r_min=dust_geometry%r_min
                    gas_geometry%r_max=dust_geometry%r_max
                    gas_geometry%v_power=dust_geometry%v_power
                    gas_geometry%rho_power=dust_geometry%rho_power
                    gas_geometry%emis_power=dust_geometry%emis_power
                    gas_geometry%v_max=dust_geometry%v_max
                end if

                !allocate memory for array to store number of packets to be emitted in each cell
                if (gas_geometry%clumped_mass_frac == 1) then
                    allocate(num_packets_array(mothergrid%tot_cells,2))
                end if

            case("torus")
                if (lg_decoupled) then
                    print*, 'you have selected a torus distribution of gas.  this routine has not been written yet.  &
                & it is due to be added to the class_grid module in due course.'
                    write(55,*) 'you have specified a torus distribution of gas.  damocles is not yet capable of creating this grid.  &
                & it is due to be added to the class_grid module.  aborted.'
                    stop
                end if
            case("arbitrary")

                if (lg_decoupled) then
                    print*,"haven't yet included provision for two different arbitrary distributions of dust and gas. aborted."
                    stop
                else

                   !update the gas geometry parameters
                   !these values have been read in from the dust grid file
                   !all other parameters are taken from the gas.in file
                    gas_geometry%v_max=dust_geometry%v_max
                    gas_geometry%r_max=dust_geometry%r_max
                    gas_geometry%v_power=dust_geometry%v_power
                    gas_geometry%r_max_cm=gas_geometry%r_max*1e15

                    !calculate the number of packets to emit in each cell
                    print*,'calculating number of packets to be emitted in each cell...'
                    allocate(num_packets_array(mothergrid%tot_cells,2))
                    num_packets_array=0
                    if (gas_geometry%emis_power == 0) then
                       !handle the strange case of a uniform emissivity grid
                       print*, 'WARNING: You have created a uniform emissivity grid.  & 
                            & To couple the emissivity grid to the dust grid, please  & 
                            & enter a non-zero emissivity power in the gas input file.&
                            & Note that this is cuboidal distribution, not a shell.'
                       if (n_packets .ge. mothergrid%tot_cells) then
                          num_packets_array(:,1) = nint(real(n_packets)/real(mothergrid%tot_cells))
                       else
                          print*,'ERROR: You have specified fewer packets than there & 
                               & are cells in the grid.  Please increase the number  &
                               & of packets to be used if you would like a uniform emissivity grid. Aborted.'
                          STOP
                       end if
                       num_packets_array(:,1) = nint(real(n_packets)/real(mothergrid%tot_cells))
                    else
                       !the normal case of an emissivity grid coupled to some power of the dust grid
                       num_packets_array(:,1) = nint((grid_cell(:)%n_rho**gas_geometry%emis_power)*real(n_packets)*(grid_cell(:)%vol/sum((grid_cell(:)%n_rho**gas_geometry%emis_power)*(grid_cell(:)%vol))))
                       num_packets_array(1,2) = num_packets_array(1,1)
                    end if

                    !calculate the cumulative number of packets to be emitted
                    do iG = 2, mothergrid%tot_cells
                       num_packets_array(iG,2) = num_packets_array(iG-1,2)+num_packets_array(iG,1)
                    end do
                 end if

                 !update the total number of packets to be run
                 n_packets = sum(num_packets_array(:,1))

            case("bipolar")
                if (lg_decoupled) then
                    print*, 'you have selected a bipolar distribution of gas.  this routine has not been written yet.  &
                & it is due to be added to the class_grid module in due course.'
                    write(55,*) 'you have specified a bipolar distribution of gas.  damocles is not yet capable of creating this grid.  &
                & it is due to be added to the class_grid module.  aborted.'

                    stop
                end if
            case default
                if (lg_decoupled) then
                    print*, 'please specify one of the following options for the geometry "shell", "torus", "arbitrary" or "bipolar". aborted.'
                    write(55,*) 'aborted - dust geometry option was not one of "shell", "torus", "arbitrary" or "bipolar".'
                    stop
                end if
        end select

    end subroutine build_emissivity_dist

end module class_grid
