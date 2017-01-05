MODULE class_grid

    use globals
    use class_geometry
    use class_dust
    use class_freq_grid

    implicit none

    TYPE grid_obj

        INTEGER          ::  n_cells(3)                    !number of cells in x/y/z directions
        REAL             ::  cell_width(3)                !width of grid cells in x/y/z directions (only applicable for cubic grids with equal no of divisions in each axis)
        REAL             ::  cell_vol                     !volume of grid cell (only applicable for cubic grids with equal no of divisions in each axis)
        INTEGER          ::  tot_cells                     !total number of cells
        REAL,ALLOCATABLE ::  x_div(:),y_div(:),z_div(:)   !locations of divisions listed consecutively in x, y and z
        REAL             ::  x_min,x_max
        REAL             ::  y_min,y_max
        REAL             ::  z_min,z_max

    END TYPE grid_obj

    TYPE(grid_obj)  ::  mothergrid                        !properties of the mothergrid

    TYPE grid_cell_obj
        REAL    ::  axis(3)                               !limits of grid cell in x, y and z
        REAL    ::  width(3)                              !width of grid cells in x/y/z directions
        REAL    ::  rho,nrho,N_e                          !mass density, number density and electron density of grid cell
        REAL    ::  r                                     !radial distance from (0,0,0) to the centre of the cell
        REAL    ::  vol                                   !volume of cell
        INTEGER ::  id(3)                                 !grid cell number in each of x, y and z
        INTEGER ::  cellStatus,numPhots
    END TYPE

    TYPE(grid_cell_obj),ALLOCATABLE  :: grid_cell(:) !mothergrid is comprised of'tot_cells' grid_cells

    REAL        ::  tot_vol                          !total volume of supernova in 1e42cm^3




    !!variables declared below need review (removing, moving, etc.)
    INTEGER ::  nclumps                     !actual number of clumps used (increased over iterations
                                          !until within 99.5% of theoretical number)
    INTEGER ::  loop                        !loop counter for clump iterations
    !INTEGER ::  numpG                       !photon number density in each cell
    INTEGER ::  isGx,isGy,isGz              !subgrid loop counters
    INTEGER ::  iS,iSG                      !iS = number of subgrid, iSG = ID of cell in subgrid

    REAL    ::  SF              !scale factors used for normalising

    REAL    ::  m                           !calculated mass of dust using densities and vols to check correct
    REAL    ::  h,micm,mcl,prob
    REAL    ::  pp
    REAL    ::  rhodG,ndust


    REAL    ::  rho_in_icm                  !density at inner radius for dust not in clumps in smooth distribution

    REAL    ::  msub                        !mass in grid replaced by clumps

    REAL    ::  cellno
    REAL    ::  ES_const



    INTEGER(8)     ::  iP,idP,n,iGP,idPP,n_inactive,freq(1),nsize, &
        & ios,nabs,iDoublet
    INTEGER :: n_threads
    REAL :: E_0, &
        & SCAT_RF_PT(2),lambda_bin,vel_bin, &
        & tot,shell_width,  &
        & const,dummy
    REAL,DIMENSION(:,:),ALLOCATABLE :: tmp, &
        & grain_rad,Qext,Qsca, &
        & grain_id, &
        & RSh,ab
    INTEGER :: ncl
    REAL,DIMENSION(:),ALLOCATABLE :: grid, r

    REAL ::SCAT_RF_SN(2), &
        & SCAT_FIN(2),pos_sph1(3),R_P,ndustav

    INTEGER(8),DIMENSION(:,:),ALLOCATABLE :: ix,iy,iz
    REAL,DIMENSION(:),ALLOCATABLE      :: NP_BIN,diff
    INTEGER(8),DIMENSION(:),ALLOCATABLE :: NP

    INTEGER :: test
    CHARACTER(LEN=1024) :: filename,junk

    INTEGER ::  shell_no

contains

    SUBROUTINE build_dust_grid()

        SELECT CASE(dust_geometry%type)

            CASE("shell")

                !open files to write gridcell coords to (in cm)
                OPEN(32,file='grid.in')

                        !calculate R_max and R_min (dust) based on maximum velocity, day no (d=v*t) and R_min/R_max ratio
                dust_geometry%R_max=dust_geometry%v_max*day_no*8.64E-6
                dust_geometry%R_min=dust_geometry%R_ratio*dust_geometry%R_max

                !check Rin/Rout ratio
                IF (dust_geometry%R_min>dust_geometry%R_max) THEN
                    PRINT*, "Please specify an R_min/R_max ratio that is less than 1.  Aborted."
                    STOP
                END IF

                !convert supernova bounds from e15cm to cm
                dust_geometry%R_min_cm=dust_geometry%R_min*1e15
                dust_geometry%R_max_cm=dust_geometry%R_max*1e15

                !convert dust mass from M_sun to grams
                dust%mass_grams=dust%mass*1.98855e33

                !set bounds of grid to be radius of SN
                mothergrid%x_min=-1*MAX(dust_geometry%R_max_cm,gas_geometry%R_max*1e15)
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
                tot_vol=1000*4*pi*(dust_geometry%R_max**3-dust_geometry%R_min**3)/3 !in e42cm^3

                ALLOCATE(grid_cell(mothergrid%tot_cells))
                ALLOCATE(mothergrid%x_div(mothergrid%n_cells(1)))
                ALLOCATE(mothergrid%y_div(mothergrid%n_cells(2)))
                ALLOCATE(mothergrid%z_div(mothergrid%n_cells(3)))

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
                DO ixx=1,mothergrid%n_cells(1)
                    mothergrid%x_div(ixx)=mothergrid%x_min+((ixx-1)*mothergrid%cell_width(1))
                END DO

                DO iyy=1,mothergrid%n_cells(2)
                    mothergrid%y_div(iyy)=mothergrid%y_min+((iyy-1)*mothergrid%cell_width(2))
                END DO

                DO izz=1,mothergrid%n_cells(3)
                    mothergrid%z_div(izz)=mothergrid%z_min+((izz-1)*mothergrid%cell_width(3))
                END DO

                !calculate radius of the centre of each cell
                iG=0
                grid_cell%r=0
                DO ixx=1,mothergrid%n_cells(1)
                    DO iyy=1,mothergrid%n_cells(2)
                        DO izz=1,mothergrid%n_cells(3)
                            iG=iG+1
                            grid_cell(iG)%r=((mothergrid%x_div(ixx)+mothergrid%cell_width(1)/2)**2+(mothergrid%y_div(iyy)+mothergrid%cell_width(2)/2)**2+(mothergrid%z_div(izz)+mothergrid%cell_width(3)/2)**2)**0.5
                            grid_cell(iG)%axis(:)=(/ mothergrid%x_div(ixx),mothergrid%y_div(iyy),mothergrid%z_div(izz)/)
                            !print*,ixx,(mothergrid%x_div(ixx)+mothergrid%cell_width(1)/2)**2,iyy,(mothergrid%y_div(iyy)+mothergrid%cell_width(2)/2)**2,izz,(mothergrid%z_div(izz)+mothergrid%cell_width(3)/2)**2
                        END DO
                    END DO
                END DO

                !!!!!!!!!!!!!!!!!!!!!!!!!
                !initialise everything to 0

                nclumps=0

                m=0
                micm=0
                mcl=0
                ncl=0
                prob=0
                pp=0
                rhodG=0
                ndust=0
                h=0


                !set counters to zero
                iG=0
                SF=0
                !SF2=0
                !set counters to zero
                n=0
                iG=0
                h=0.
                ndustav=0.
                micm=0
                mcl=0
                !!!!!!!!!!!!!!!!!!!!!!!!

                !calculate desired mass of dust in clumps and in interclump medium (icm)
                dust%m_icm=dust%mass*(1-dust_geometry%clumped_mass_frac)

                !calculate number of clumps based on filling fraction of total volume and volume of clump (=volume of grid cell)
                dust_geometry%n_clumps=FLOOR(dust_geometry%ff*tot_vol/mothergrid%cell_vol)

                !calculate mass of a clump based on dust mass fraction in clumps and total number of clumps
                dust%m_clump=(dust%mass*dust_geometry%clumped_mass_frac)/dust_geometry%n_clumps

                !calculate dust density at inner radius (rho_in)
                IF (.not. dust_geometry%lg_clumped) THEN
                    !factor of 1.989e-12 to convert from mass in Msun and vol in e45cm3 to g and cm3
                    !1.989e-12 = 1.989e33/1e45 (g in Msun / vol e45cm3 -> cm3)
                    IF (dust_geometry%rho_power==3) THEN
                        dust_geometry%rho_in=(dust_geometry%R_min**(-dust_geometry%rho_power))*((dust%mass*1.989e33)/(LOG(dust_geometry%R_max/dust_geometry%R_min)*4*pi))
                        dust_geometry%rho_in=dust_geometry%rho_in*(1.989e-12)
                    ELSE
                        dust_geometry%rho_in=(dust_geometry%R_min**(-dust_geometry%rho_power))*((dust%mass*(3-dust_geometry%rho_power)) &
                            & /(4*pi*(dust_geometry%R_max**(3-dust_geometry%rho_power)-dust_geometry%R_min**(3-dust_geometry%rho_power))))
                        dust_geometry%rho_in=dust_geometry%rho_in*(1.989e-12)
                    END IF
                ELSE
                    !calculate dust density at inner radius (rho_in_icm) for smooth interclump medium (excluding clumps)
                    IF (dust_geometry%rho_power==3) THEN
                        !!check the factor of 1+ff...???
                        !                        rho_in_icm=dust_geometry%R_min**(-dust_geometry%rho_power)*(1+dust_geometry%ff)*((dust%m_icm)/(LOG(dust_geometry%R_max/dust_geometry%R_min)*4*pi))
                        rho_in_icm=(dust%m_icm*dust_geometry%R_min**(-dust_geometry%rho_power))/(4*pi*(1-dust_geometry%ff)*LOG(dust_geometry%R_max/dust_geometry%R_min))
                        rho_in_icm=rho_in_icm*1.989e-12
                    ELSE
                        rho_in_icm=dust_geometry%R_min**(-dust_geometry%rho_power)*(1+dust_geometry%ff)*((dust%m_icm)*(3-dust_geometry%rho_power) &
                            & /(4*pi*(dust_geometry%R_max**(3-dust_geometry%rho_power)-dust_geometry%R_min**(3-dust_geometry%rho_power))))
                        !rho_in_icm=dust_geometry%R_min**(-dust_geometry%rho_power)*((dust%m_icm)*(3-dust_geometry%rho_power) &
                        !    & /(4*pi*(dust_geometry%ff)*(dust_geometry%R_max**(3-dust_geometry%rho_power)-dust_geometry%R_min**(3-dust_geometry%rho_power))))
                        rho_in_icm=rho_in_icm*1.989e-12
                    END IF
                END IF

                !                IF (dust_geometry%clumped_mass_frac==1.0) THEN
                !                    dust_geometry%den_con=0.0
                !                ELSE
                !                    dust_geometry%den_con=dust%m_clump/(rho_in_icm*5.02765e8*mothergrid%cell_vol)       !5.02765e8=1e42/1.989e33 i.e. conversion factor for  e42cm3 to cm3 and g to Msun
                !                END IF

                dust_geometry%rho_clump=dust%m_clump/(mothergrid%cell_vol*5.02765e8)

                !CALCULATE NUMBER OF PHOTONS EMITTED IN EACH CELL WITHIN RADIAL BOUNDS
                grid_cell(:)%cellStatus=0
                prob=0
                msub=0
                nclumps=0
                h=0

                SF=((dust_geometry%R_max**(1-dust_geometry%clump_power)-dust_geometry%R_min**(1-dust_geometry%clump_power)))/(dust_geometry%R_min**(-dust_geometry%clump_power)*(1-dust_geometry%clump_power))

                IF (dust_geometry%lg_clumped) THEN

                    DO WHILE (nclumps<(dust_geometry%n_clumps))

                        iG=0

                        call RANDOM_NUMBER(cellno)
                        iG=ceiling(mothergrid%tot_cells*cellno)
                        IF (iG==0) CYCLE
                        IF ((grid_cell(iG)%r<(dust_geometry%R_max_cm)) .AND. (grid_cell(iG)%r>(dust_geometry%R_min_cm))) THEN
                            h=h+1
                            prob=((dust_geometry%R_min_cm/grid_cell(iG)%r)**dust_geometry%clump_power)/SF
                            call RANDOM_NUMBER(pp)
                            !PRINT*,prob,pp
                            IF ((pp<prob) .AND. (grid_cell(iG)%cellStatus/=1)) THEN  !(i.e. if set to be a clump and not already a clump)
                                 !CLUMP!

                                grid_cell(iG)%cellStatus=1
                                isG=0
                                nclumps=nclumps+1
                                rhodG=dust_geometry%rho_clump
                                mcl=mcl+rhodG*mothergrid%cell_vol*5.02765e8
                                ncl=ncl+1
                                msub=msub+(mothergrid%cell_vol*5.02765e8*rho_in_icm*(dust_geometry%R_min_cm/grid_cell(iG)%r)**dust_geometry%rho_power)
                                !PRINT*,'clump',ncl,'of',dust_geometry%n_clumps

                            END IF
                        END IF
                    END DO


                    PRINT*,'average mass density (including clumps) (g/cm3)',(dust%mass_grams*1e-14)/(tot_vol*1e28)
                    PRINT*,'icm mass density at inner radius (g/cm3)',rho_in_icm
                    PRINT*,'icm mass density at outer radius (g/cm3)',(rho_in_icm)*(dust_geometry%R_min/dust_geometry%R_max)**dust_geometry%rho_power
                    PRINT*,'mass of interclump medium (Msun)',dust%m_icm
                    PRINT*,'mass in clumps (Msun)',dust%m_clump*dust_geometry%n_clumps
                    PRINT*,'density constrast',dust_geometry%den_con

                    iG=0
                    h=0
                    m=0


                    micm=0
                    DO ixx=1,mothergrid%n_cells(1)
                        DO iyy=1,mothergrid%n_cells(2)
                            DO izz=1,mothergrid%n_cells(3)
                                iG=iG+1

                                IF ((grid_cell(iG)%r<(dust_geometry%R_max_cm)) .AND. (grid_cell(iG)%r>(dust_geometry%R_min_cm))) THEN
                                    h=h+1
                                    IF (grid_cell(iG)%cellStatus==0) THEN
                                        rhodG=rho_in_icm*(dust_geometry%R_min_cm/grid_cell(iG)%r)**dust_geometry%rho_power
                                        micm=micm+rhodG*mothergrid%cell_vol*5.02765e8
                                        loop=loop+1
                                    END IF
                                    m=m+rhodG*mothergrid%cell_vol*5.02765e8
                                    ndust=rhodG/dust%av_mgrain
                                    ndustav=ndustav+ndust

                                    grid_cell(iG)%rho=rhodG
                                    grid_cell(iG)%nrho= ndust
                                    WRITE(32,*) grid_cell(iG)%axis(:),grid_cell(iG)%rho
                                    !!!when looking at clumping need to include ES
                                    grid_cell(ig)%id(:)=(/ ixx,iyy,izz /)
                                    !this should be equal to numpG not 0 but I've left out the calculation...
                                    grid_cell(iG)%numPhots= 0

                                ! at the moment leave out subgrids - not actually sure need them as not full RT
                                !                        IF (mgrid(iG)%cellStatus==1) THEN
                                !                            DO iS=1,3
                                !                                mgrid(iG)%subaxes(iS,1)=mgrid(iG)%axis(iS)
                                !                                mgrid(iG)%subaxes(iS,2)=mgrid(iG)%axis(iS)+width(iS)/2
                                !                            END DO
                                !                            isG=0
                                !                            DO isGx=1,2
                                !                                DO isGy=1,2
                                !                                    DO isGz=1,2
                                !                                        isG=isG+1
                                !                                        mgrid(iG)%subgrid(isG)%axis(1)=mgrid(iG)%axis(1)+((isGx-1)*mothergrid%cell_width(1))/2
                                !                                        mgrid(iG)%subgrid(isG)%axis(2)=mgrid(iG)%axis(2)+((isGy-1)*mothergrid%cell_width(2))/2
                                !                                        mgrid(iG)%subgrid(isG)%axis(3)=mgrid(iG)%axis(3)+((isGz-1)*mothergrid%cell_width(3))/2
                                !                                        mgrid(iG)%subgrid(isG)%rho=rhodG
                                !                                        mgrid(iG)%subgrid(isG)%nrho=ndust
                                !                                    END DO
                                !                                END DO
                                !                            END DO
                                !                        END IF
                                    !WRITE(32,*) mgrid(ig)%id(:),mgrid(iG)%numPhots,mgrid(iG)%axis(:),mgrid(iG)%rho,mgrid(iG)%nrho
                                ELSE
                                    grid_cell(iG)%axis(:)=(/ mothergrid%x_div(ixx),mothergrid%y_div(iyy),mothergrid%z_div(izz)/)
                                    grid_cell(iG)%rho=0.
                                    grid_cell(iG)%nrho=0
                                    grid_cell(ig)%id(:)=(/ ixx,iyy,izz /)
                                    grid_cell(iG)%numPhots= 0
                                    !WRITE(32,*) mgrid(ig)%id(:),mgrid(iG)%numPhots,mgrid(iG)%axis(:),mgrid(iG)%rho,mgrid(iG)%nrho
                                END IF
                            END DO
                        END DO
                    END DO

                ELSE
                    !not clumping
                    h=0
                    DO ixx=1,mothergrid%n_cells(1)
                        DO iyy=1,mothergrid%n_cells(2)
                            DO izz=1,mothergrid%n_cells(3)
                                iG=iG+1
                                IF ((grid_cell(iG)%r<(dust_geometry%R_max_cm)) .AND. (grid_cell(iG)%r>(dust_geometry%R_min_cm))) THEN
                                    h=h+1
                                    !this was the old incorrect calculation
                                    !m=m+((dust%mass_grams*((dust_geometry%R_min_cm)/grid_cell(iG)%r)**q)/(SF))/1.989E33
                                    !rhodG=((dust%mass_grams*1.0e-14*1e-28*((dust_geometry%R_min_cm)/grid_cell(iG)%r)**q)/(SF*mothergrid%cell_vol))
                                    !corrected here (virtually no difference...)
                                    rhodG=((dust_geometry%rho_in)*(dust_geometry%R_min_cm/grid_cell(iG)%r)**dust_geometry%rho_power)
                                    m=m+rhodG*mothergrid%cell_vol*5.02765e8         !5.02765e8=1e42/1.989e33 i.e. conversion factor for  e42cm3 to cm3 and g to Msun
                                    ndust=rhodG/dust%av_mgrain
                                    ndustav=ndustav+ndust
                                    grid_cell(iG)%axis(:)=(/ mothergrid%x_div(ixx),mothergrid%y_div(iyy),mothergrid%z_div(izz)/)
                                    grid_cell(iG)%rho=rhodG
                                    grid_cell(iG)%nrho= ndust
                                    !!!ES EDIT
                                    !!!watch out for the E19 ES_const in E20

                                    grid_cell(iG)%N_e=ES_const*(grid_cell(iG)%r**(-dust_geometry%rho_power))*1E20
                                    !PRINT*,ES_const,(grid_cell(iG)%r**(-q)),mgrid(iG)%N_e, grid_cell(iG)%r,dust_geometry%R_min_cm
                                    grid_cell(ig)%id(:)=(/ ixx,iyy,izz /)
                                    !this should be calculates as numpG and not 0, but I've left out the calculation...
                                    grid_cell(iG)%numPhots= 0
                                    !WRITE(32,*) mgrid(ig)%id(:),mgrid(iG)%numPhots,mgrid(iG)%axis(:),mgrid(iG)%rho,mgrid(iG)%nrho
                                ELSE
                                    grid_cell(iG)%axis(:)=(/ mothergrid%x_div(ixx),mothergrid%y_div(iyy),mothergrid%z_div(izz)/)
                                    grid_cell(iG)%rho=0.
                                    grid_cell(iG)%nrho=0
                                    grid_cell(ii)%N_e=0
                                    grid_cell(ig)%id(:)=(/ ixx,iyy,izz /)
                                    grid_cell(iG)%numPhots= 0
                                    !WRITE(32,*) mgrid(ig)%id(:),mgrid(iG)%numPhots,mgrid(iG)%axis(:),mgrid(iG)%rho,mgrid(iG)%nrho
                                END IF
                                WRITE(32,*) grid_cell(iG)%axis(:),grid_cell(iG)%rho
                            END DO
                        END DO
                    END DO
                    PRINT*,'DUST GRAIN NUMBER DENSITY AT Rin',dust_geometry%rho_in/dust%av_mgrain
                    PRINT*,'DUST GRAIN NUMBER DENSITY AT Rout',((dust_geometry%rho_in)*(dust_geometry%R_min_cm/dust_geometry%R_max_cm)**dust_geometry%rho_power)/dust%av_mgrain
                END IF

                IF (dust_geometry%lg_clumped) THEN
                    PRINT*,'number clumps test:',' using - ',ncl,'requested -',dust_geometry%n_clumps
                    PRINT*,'mass clumps: ','using - ',mcl,'requested: ',dust%m_clump*dust_geometry%n_clumps
                    PRINT*,'mass ICM test:',' using - ',micm,'requested - ',dust%m_icm
                    PRINT*,'filling factor',dust_geometry%ff
                    PRINT*,'mass check (calculated as rho*V)','using - ',m,'requested -',dust%mass
                    PRINT*,'substituted mass',msub
                ELSE
                    PRINT*,'mass check (calculated as rho*V)','using - ',m,'requested -',dust%mass
                END IF
                ndustav=ndustav/h
                PRINT*,'no of grid cells inside SN',h
                PRINT*,'volume of total grid cells inside SN (e42cm)',h*mothergrid%cell_vol
                PRINT*,'average dust grain density per cell (including any clumps)',ndustav
                PRINT*,''

                CLOSE(32)

            CASE("torus")
                PRINT*, 'You have selected a torus distribution of dust.  This routine has not been written yet.  &
                & It is due to be added to the class_grid module in due course.'
                WRITE(55,*) 'You have specified a torus distribution of dust.  Damocles is not yet capable of creating this grid.  &
                & It is due to be added to the class_grid module.  Aborted.'
                STOP
            CASE("arbitrary")
                !read in dust grid file
                !!note that currently the filename is hard coded and should be changed to be variable
                OPEN(33,file='dust_grid.in')

                !the file format is based on the grid generated by the mocassin grid generator which can be found at
                !http://www.nebulousresearch.org/codes/mocassin/mocassin_gridmaker.php
                READ(33,*) junk,mothergrid%n_cells(1),mothergrid%n_cells(2),mothergrid%n_cells(3)
                mothergrid%tot_cells=mothergrid%n_cells(1)*mothergrid%n_cells(2)*mothergrid%n_cells(3)

                ALLOCATE(grid_cell(mothergrid%tot_cells))
                ALLOCATE(mothergrid%x_div(mothergrid%n_cells(1)))
                ALLOCATE(mothergrid%y_div(mothergrid%n_cells(2)))
                ALLOCATE(mothergrid%z_div(mothergrid%n_cells(3)))


                xx=0
                yy=0
                zz=0
                !grid must be arranged in ascending z, then y, then x (as per mocassin grid maker)
                !!include a sort here to ensure this is the case?
                DO iG=1,mothergrid%tot_cells

                    !read in limits of grid cells in x, y and z
                    READ(33,*) grid_cell(iG)%axis(1),grid_cell(iG)%axis(2),grid_cell(iG)%axis(3),grid_cell(iG)%nrho

                    !this section calculates the cell no in each axis given the overall cell number within the grid
                    !(effectively reverse engineers iG=(n_cells(2)*n_cells(3)*(id(1)-1)+n_cells(3)*(id(2)-1)+id(3)))
                    IF (MOD(iG,mothergrid%n_cells(3)) /=1) THEN
                    grid_cell(iG)%id(3)=grid_cell(iG-1)%id(3)+1
                    ELSE
                    grid_cell(iG)%id(3)=1
                    END IF
                    grid_cell(iG)%id(2)=MOD((iG-grid_cell(iG)%id(3))/mothergrid%n_cells(3),mothergrid%n_cells(2))+1
                    grid_cell(iG)%id(1)=((iG-grid_cell(iG)%id(3)-(mothergrid%n_cells(3)*(grid_cell(iG)%id(2)-1)))/(mothergrid%n_cells(3)*mothergrid%n_cells(2)))+1

                    !this section populates the list of divisions in x, y and z from the grid specified
                    !each time a cell is encountered where the id of the other two axes is 1, the third non-zero axis is stored
                    !this generates a unique list of axis divisions
                    IF (grid_cell(iG)%id(1) == 1 .and. grid_cell(iG)%id(2) == 1) THEN
                        zz=zz+1
                        mothergrid%z_div(zz)=grid_cell(iG)%axis(3)
                    END IF
                    IF (grid_cell(iG)%id(1) == 1 .and. grid_cell(iG)%id(3) == 1) THEN
                        yy=yy+1
                        mothergrid%y_div(yy)=grid_cell(iG)%axis(2)
                    END IF
                    IF (grid_cell(iG)%id(2) == 1 .and. grid_cell(iG)%id(3) == 1) THEN
                        xx=xx+1
                        mothergrid%x_div(xx)=grid_cell(iG)%axis(1)
                    END IF
                END DO

                CLOSE(33)

                !radii are calculated for each grid cell from the centre of the grid (0,0,0) to the centre of the cell
                !volumes are calculated by multiplying distance to next division in all axes
                !at far 'right' edge, the volumes, widths and radii are calculated using symmetries of grid
                !i.e. radius/volume at cell (43,47,50) is same as for (43,47,1)
                DO iG=1,mothergrid%tot_cells
                    IF (grid_cell(iG)%id(1) /= mothergrid%n_cells(1) .and. &
                      & grid_cell(iG)%id(2) /= mothergrid%n_cells(2) .and. &
                      & grid_cell(iG)%id(3) /= mothergrid%n_cells(3)) THEN

                        grid_cell(iG)%r=((((grid_cell(iG)%axis(1)+grid_cell(iG+1)%axis(1))/2)**2) + &
                            & (((grid_cell(iG)%axis(2)+grid_cell(iG+1)%axis(2))/2)**2) + &
                            & (((grid_cell(iG)%axis(3)+grid_cell(iG+1)%axis(3))/2)**2))**0.5
                        grid_cell(iG)%vol=(mothergrid%x_div(grid_cell(iG)%id(1)+1)-mothergrid%x_div(grid_cell(iG)%id(1)))*1e-14* &
                                        & (mothergrid%y_div(grid_cell(iG)%id(2)+1)-mothergrid%y_div(grid_cell(iG)%id(2)))*1e-14* &
                                        & (mothergrid%z_div(grid_cell(iG)%id(3)+1)-mothergrid%z_div(grid_cell(iG)%id(3)))*1e-14
                        grid_cell(iG)%width(1)=mothergrid%x_div(grid_cell(iG)%id(1)+1)-mothergrid%x_div(grid_cell(iG)%id(1))
                        grid_cell(iG)%width(2)=mothergrid%x_div(grid_cell(iG)%id(2)+1)-mothergrid%x_div(grid_cell(iG)%id(2))
                        grid_cell(iG)%width(3)=mothergrid%x_div(grid_cell(iG)%id(3)+1)-mothergrid%x_div(grid_cell(iG)%id(3))
                    ELSE IF (grid_cell(iG)%id(3) == mothergrid%n_cells(3)) THEN
                        grid_cell(iG)%vol=grid_cell(iG-mothergrid%n_cells(3)+1)%vol
                        grid_cell(iG)%r=grid_cell(iG-mothergrid%n_cells(3)+1)%r
                    ELSE IF (grid_cell(iG)%id(2) == mothergrid%n_cells(2)) THEN
                        grid_cell(iG)%vol=grid_cell(iG-mothergrid%n_cells(2)*(grid_cell(iG)%id(2)-1))%vol
                        grid_cell(iG)%r=grid_cell(iG-mothergrid%n_cells(2)*(grid_cell(iG)%id(2)-1))%r
                    ELSE IF (grid_cell(iG)%id(1) == mothergrid%n_cells(1)) THEN
                        grid_cell(iG)%vol=grid_cell(iG-mothergrid%n_cells(2)*mothergrid%n_cells(1)*(grid_cell(iG)%id(1)-1))%vol
                        grid_cell(iG)%r=grid_cell(iG-mothergrid%n_cells(2)*mothergrid%n_cells(1)*(grid_cell(iG)%id(1)-1))%r
                    END IF
                END DO

            PRINT*,'Dust mass in grid (Msun)',sum(grid_cell(:)%nrho*(grid_cell(:)%vol)*dust%av_mgrain)*5.02765e8         !5.02765e8=1e42/1.989e33 i.e. conversion factor for  e42cm3 to cm3 and g to Msun)

            CASE("bipolar")
                PRINT*, 'You have selected a bipolar distribution of dust.  This routine has not been written yet.  &
                & It is due to be added to the class_grid module in due course.'
                WRITE(55,*) 'You have specified a bipolar distribution of dust.  Damocles is not yet capable of creating this grid.  &
                & It is due to be added to the class_grid module.  Aborted.'
                STOP
            CASE DEFAULT
                PRINT*, 'Please specify one of the following options for the geometry "shell", "torus", "arbitrary" or "bipolar". Aborted.'
                WRITE(55,*) 'Aborted - dust geometry option was not one of "shell", "torus", "arbitrary" or "bipolar".'
                STOP
        END SELECT


    END SUBROUTINE
    SUBROUTINE build_emissivity_dist()
        ALLOCATE(NP_BIN(nu_grid%n_bins))
        ALLOCATE(tmp(nu_grid%n_bins,1))
        tmp=0
        NP_BIN=0
        SELECT CASE(gas_geometry%type)

            CASE("shell")
                !assign parameters related to gas geometry
                IF (lg_decoupled) THEN
                    !if decoupled then generate max/min radii from epoch and maximum velocity
                    gas_geometry%R_max=gas_geometry%v_max*day_no*8.64E-6
                    gas_geometry%R_min=gas_geometry%R_ratio*gas_geometry%R_max
                ELSE
                    !if decoupled then set gas geometry parameters to equal the dust geometry parameters
                    IF (gas_geometry%type /= dust_geometry%type) THEN
                        PRINT*, 'You have requested that gas and dust distributions be coupled but specified different geometries.  Aborted.'
                        STOP
                    END IF
                    gas_geometry%R_ratio=dust_geometry%R_ratio
                    gas_geometry%R_min=dust_geometry%R_min
                    gas_geometry%R_max=dust_geometry%R_max
                    gas_geometry%v_power=dust_geometry%v_power
                    gas_geometry%rho_power=dust_geometry%rho_power
                    gas_geometry%emis_power=dust_geometry%emis_power
                    gas_geometry%v_max=dust_geometry%v_max
                END IF

                IF (gas_geometry%clumped_mass_frac /= 1) THEN
                    ALLOCATE(NP(n_shells+1))
                    ALLOCATE(RSh(n_shells+1,2))
                ELSE IF (gas_geometry%clumped_mass_frac == 1) THEN
                    ALLOCATE(NP(mothergrid%tot_cells))
                END IF

                RSh=0
                shell_width=(gas_geometry%R_max-gas_geometry%R_min)/n_shells
                !calculate upper and lower radius bound for each shell
                RSh(1,1)=gas_geometry%R_min
                RSh(1,2)=gas_geometry%R_min+shell_width
                DO ii=1,n_shells
                    RSh(ii+1,1:2)=(/ RSh(ii,2),RSh(ii,2)+shell_width /)
                END DO

                !scale factor to work out number of packets in each shell
                IF ((gas_geometry%emis_power*gas_geometry%rho_power)==3) THEN
                    const=n_packets/(LOG(gas_geometry%R_max/gas_geometry%R_min))
                    NP(:)=NINT(const*LOG(RSh(:,2)/RSh(:,1)))
                ELSE
                    const=n_packets*(gas_geometry%emis_power*gas_geometry%rho_power-3)/(gas_geometry%R_min**(3-gas_geometry%emis_power*gas_geometry%rho_power)-gas_geometry%R_max**(3-gas_geometry%emis_power*gas_geometry%rho_power))
                    NP(:)=NINT(const*(RSh(:,1)**(3-gas_geometry%emis_power*gas_geometry%rho_power)-RSh(:,2)**(3-gas_geometry%emis_power*gas_geometry%rho_power))/(gas_geometry%emis_power*gas_geometry%rho_power-3))
                END IF

            CASE("torus")
                IF (lg_decoupled) THEN
                    PRINT*, 'You have selected a torus distribution of gas.  This routine has not been written yet.  &
                & It is due to be added to the class_grid module in due course.'
                    WRITE(55,*) 'You have specified a torus distribution of gas.  Damocles is not yet capable of creating this grid.  &
                & It is due to be added to the class_grid module.  Aborted.'
                    STOP
                END IF
            CASE("arbitrary")

                IF (lg_decoupled) THEN
                    PRINT*,"Haven't yet included provision for two different arbitrary distributions of dust and gas. Aborted."
                    STOP
                ELSE
                    ALLOCATE(NP(mothergrid%tot_cells))
                    NP=0
                    PRINT*,'calculating number of packets to be emitted in each cell...'
                    !number of packets to be emitted in each cell scaled with number of particles in cell
                    NP(:)=NINT(grid_cell(:)%nrho*real(n_packets)*(grid_cell(ii)%vol/SUM((grid_cell%nrho)*(grid_cell%vol))))
                END IF

                !DO ii=1,mothergrid%tot_cells
                !PRINT*,grid_cell(ii)%axis(:),grid_cell(ii)%nrho,NP(ii)
               ! PRINT*,grid_cell(ii)%nrho,grid_cell(ii)%vol,grid_cell(ii)%nrho*grid_cell(ii)%vol
                !END DO

            !                IF (lg_decoupled) THEN
            !                    PRINT*, 'You have selected an arbitrary distribution of gas.  This routine has not been written yet.  &
            !                & It is due to be added to the class_grid module in due course.'
            !                    WRITE(55,*) 'You have specified an arbitrary distribution of gas.  Damocles is not yet capable of handling arbitrary geometry grids.  &
            !                & This capability is due to be added to the class_grid module.  Aborted.'
            !                    STOP
            !                ELSE
            !                    IF (gas_geometry%type /= dust_geometry%type) THEN
            !                        PRINT*, 'You have requested that gas and dust distributions be coupled but specified different geometries.  Aborted.'
            !                        STOP
            !                    END IF
            !                END IF
            CASE("bipolar")
                IF (lg_decoupled) THEN
                    PRINT*, 'You have selected a bipolar distribution of gas.  This routine has not been written yet.  &
                & It is due to be added to the class_grid module in due course.'
                    WRITE(55,*) 'You have specified a bipolar distribution of gas.  Damocles is not yet capable of creating this grid.  &
                & It is due to be added to the class_grid module.  Aborted.'

                    STOP
                END IF
            CASE DEFAULT
                IF (lg_decoupled) THEN
                    PRINT*, 'Please specify one of the following options for the geometry "shell", "torus", "arbitrary" or "bipolar". Aborted.'
                    WRITE(55,*) 'Aborted - dust geometry option was not one of "shell", "torus", "arbitrary" or "bipolar".'
                    STOP
                END IF
        END SELECT


    END SUBROUTINE build_emissivity_dist

END MODULE class_grid
