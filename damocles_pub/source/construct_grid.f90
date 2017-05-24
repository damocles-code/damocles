!SUBROUTINE TO CONSTRUCT CARTESIAN GRID
SUBROUTINE construct_grid

    USE initialise
    USE electron_scattering

    IMPLICIT NONE

    INTEGER ::  nclumps_tot                 !theoretical number of clumps calculated from filling factor
    INTEGER ::  nclumps                     !actual number of clumps used (increased over iterations
                                            !until within 99.5% of theoretical number)
    INTEGER ::  ixx,iyy,izz                 !mothergrid loop counters
    INTEGER ::  iG                          !cell ID
    INTEGER ::  loop                        !loop counter for clump iterations
    !INTEGER ::  numpG                       !photon number density in each cell
    INTEGER ::  isGx,isGy,isGz              !subgrid loop counters
    INTEGER ::  iS,iSG                      !iS = number of subgrid, iSG = ID of cell in subgrid

    REAL    ::  SF,SF2,SFndust              !scale factors used for normalising
    REAL    ::  M_icm                       !mass of the inter clump medium (ICM)
    REAL    ::  m                           !calculated mass of dust using densities and vols to check correct
    REAL    ::  h,micm,mcl,prob
    REAL    ::  pp
    REAL    ::  rhodG,ndust
    REAL    ::  r(totcells)
    REAL    ::  rho_in                      !density at inner radius (for total dust in smooth distribution)
    REAL    ::  rho_in_icm                  !density at inner radius for dust not in clumps in smooth distribution
    REAL    ::  m_clump                     !mass of a clump
    REAL    ::  msub                        !mass in grid replaced by clumps
    REAL    ::  clump_dens                  !density of individual clump
    REAL    ::  cellno
    REAL    ::  ES_const
    INTEGER :: test

    !initialise everything to 0
    nclumps_tot=0
    nclumps=0
    M_icm=0
    m=0
    micm=0
    mcl=0
    ncl=0
    prob=0
    pp=0
    rhodG=0
    ndust=0
    r=0
    h=0

    call N_e_const(ES_const)

    !open files to write gridcell coords to (in cm)
    !OPEN(31,file='grid_values.in')
    OPEN(32,file='grid.in')

    !calculate useful distances
    totvol=1000*4*pi*(R_max**3-R_min**3)/3	!in e42cm^3

    PRINT*,'total volume of supernova in e42cm^3',totvol

    width(1)=(x_max-x_min)/ncells(1)
    width(2)=(y_max-y_min)/ncells(2)
    width(3)=(z_max-z_min)/ncells(3)
    vol=((width(1)/1e14)*(width(2)/1e14)*(width(3)/1e14)) !in 1e42cm^3

    PRINT*,'VOLUME OF GRID CELL (and therefore clump) in e42cm^3',vol
    PRINT*,'GRID CELL WIDTHS in cm: ','X',width(1),'Y',width(2),'Z',width(3)


    !calculate grid cell coords points
    DO ixx=1,ncells(1)
        grid(ixx)=x_min+((ixx-1)*width(1))
        !WRITE(31,*) grid(ixx)
    END DO

    DO iyy=1,ncells(2)
        grid(ncells(1)+iyy)=y_min+((iyy-1)*width(2))
        !WRITE(31,*) grid(ncells(1)+iyy)
    END DO

    DO izz=1,ncells(3)
        grid(ncells(1)+ncells(2)+izz)=z_min+((izz-1)*width(3))
        !WRITE(31,*) grid(ncells(1)+ncells(2)+izz)
    END DO

    !set counters to zero
    iG=0
    SF=0
    !SF2=0
    
    !calculate radius of each cell and scale factors
    DO ixx=1,ncells(1)
        DO iyy=1,ncells(2)
            DO izz=1,ncells(3)
                iG=iG+1
                r(iG)=((grid(ixx)+width(1)/2)**2+(grid(ncells(1)+iyy)+width(2)/2)**2+(grid(ncells(1)+ncells(2)+izz)+width(3)/2)**2)**0.5
            END DO
        END DO
    END DO

    !set counters to zero
    n=0
    iG=0
    h=0.
    SFndust=0.
    ndustav=0.
    micm=0
    mcl=0


    DO i=1,nspecies
        DO j=1,dust(i)%nsizes
            SFndust=SFndust+(((4*pi*dust(i)%gsize(j,1)**3*dust(i)%rhograin*1e-12)*dust(i)%gsize(j,2)/3))*dust(i)%mweight
        END DO
    END DO

    !calculate dust density at inner radius (rho_in)
    !1.989e-12 = 1.989e33/1e45 (g in Msun / vol e45cm3 -> cm3)
    IF (q==3) THEN
        rho_in=(R_min**(-q))*((MD_tot*1.989e33)/(LOG(R_max/R_min)*4*pi))
        rho_in=rho_in*(1.989e-12)
    ELSE
        rho_in=(R_min**(-q))*((MD_tot*(3-q))/(4*pi*(R_max**(3-q)-R_min**(3-q))))
        rho_in=rho_in*(1.989e-12)
    END IF

    M_icm=MD_tot*(1-mf)
    nclumps_tot=ff*totvol/vol
    m_clump=(MD_tot*mf)/nclumps_tot

    PRINT*,M_icm,MD_tot,mf
    !calculate dust density at inner radius (rho_in) for interclump medium
    IF (q==3) THEN
        rho_in_icm=R_min**(-q)*(1+ff)*((M_icm)/(LOG(R_max/R_min)*4*pi))
        rho_in_icm=rho_in_icm*1.989e-12
    ELSE
        rho_in_icm=R_min**(-q)*(1+ff)*((M_icm)*(3-q)/(4*pi*(R_max**(3-q)-R_min**(3-q))))
        rho_in_icm=rho_in_icm*1.989e-12
    END IF

    IF (mf==1.0) THEN
        dencon=0.0
    ELSE
        dencon=m_clump/(rho_in_icm*5.02765e8*vol)       !5.02765e8=1e42/1.989e33 i.e. conversion factor for  e42cm3 to cm3 and g to Msun
    END IF

    clump_dens=m_clump/(vol*5.02765e8)

    !CALCULATE NUMBER OF PHOTONS EMITTED IN EACH CELL WITHIN RADIAL BOUNDS
    mgrid(:)%cellStatus=0
    prob=0
    msub=0
    nclumps=0
    h=0

    SF=((R_max**(1-q_clump)-R_min**(1-q_clump)))/(R_min**(-q_clump)*(1-q_clump))

    IF (lgClumping) THEN
        DO WHILE (nclumps<(nclumps_tot))
            iG=0
            !SF2=0
            call RANDOM_NUMBER(cellno)
            iG=ceiling(totcells*cellno)
            IF (iG==0) CYCLE
            IF ((r(iG)<(R_max_cm)) .AND. (r(iG)>(R_min_cm))) THEN
                h=h+1
                prob=((R_min_cm/r(iG))**q_clump)/SF
                call RANDOM_NUMBER(pp)
                !PRINT*,prob,pp
                IF ((pp<prob) .AND. (mgrid(iG)%cellStatus/=1)) THEN  !(i.e. if set to be a clump and not already a clump)
                    !CLUMP!
                    
                   mgrid(iG)%cellStatus=1
                    isG=0
                    nclumps=nclumps+1
                    rhodG=clump_dens
                    mcl=mcl+rhodG*vol*5.02765e8
                    ncl=ncl+1
                    msub=msub+(vol*5.02765e8*rho_in_icm*(R_min_cm/r(iG))**q_clump)
                    PRINT*,'clump',ncl,'of',nclumps_tot
                 !ELSE
                 !   SF2=SF2+(R_min_cm/r(iG))**q
                END IF
            END IF
        END DO

        PRINT*,'average mass density (including clumps) (g/cm3)',(MD_tot_grams*1e-14)/(totvol*1e28)
        PRINT*,'icm mass density at inner radius (g/cm3)',rho_in_icm
        PRINT*,'icm mass density at outer radius (g/cm3)',(rho_in_icm)*(R_min/R_max)**q
        PRINT*,'mass of interclump medium (Msun)',M_icm
        PRINT*,'mass in clumps (Msun)',m_clump*nclumps_tot
        PRINT*,'density constrast',dencon

        iG=0
        h=0
        m=0
        test=0
        DO ixx=1,ncells(1)
            DO iyy=1,ncells(2)
                DO izz=1,ncells(3)
                    iG=iG+1
                    IF ((r(ig)<(R_max_cm)) .AND. (r(ig)>(R_min_cm))) THEN
                        h=h+1
                        IF (mgrid(iG)%cellStatus==0) THEN
                            rhodG=rho_in_icm*(R_min_cm/r(iG))**q
                            micm=micm+rhodG*vol*5.02765e8
                            loop=loop+1
                            m=m+rhodG*vol*5.02765e8
                        ELSE
                            rhodG=clump_dens
                            m=m+rhodG*vol*5.02765e8
                            test=test+1
                        END IF
                        ndust=rhodG/SFndust
                        ndustav=ndustav+ndust
                        mgrid(iG)%axis(:)=(/ grid(ixx),grid(ncells(1)+iyy),grid(ncells(1)+ncells(2)+izz)/)
                        mgrid(iG)%rho=rhodG
                        mgrid(iG)%nrho= ndust
                        WRITE(32,*) mgrid(iG)%axis(:),mgrid(iG)%rho
                        !!!when looking at clumping need to include ES
                        mgrid(ig)%id(:)=(/ ixx,iyy,izz /)
                        !this should be equal to numpG not 0 but I've left out the calculation...
                        mgrid(iG)%numPhots= 0

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
                    !                                        mgrid(iG)%subgrid(isG)%axis(1)=mgrid(iG)%axis(1)+((isGx-1)*width(1))/2
                    !                                        mgrid(iG)%subgrid(isG)%axis(2)=mgrid(iG)%axis(2)+((isGy-1)*width(2))/2
                    !                                        mgrid(iG)%subgrid(isG)%axis(3)=mgrid(iG)%axis(3)+((isGz-1)*width(3))/2
                    !                                        mgrid(iG)%subgrid(isG)%rho=rhodG
                    !                                        mgrid(iG)%subgrid(isG)%nrho=ndust
                    !                                    END DO
                    !                                END DO
                    !                            END DO
                    !                        END IF
                        !WRITE(32,*) mgrid(ig)%id(:),mgrid(iG)%numPhots,mgrid(iG)%axis(:),mgrid(iG)%rho,mgrid(iG)%nrho
                    ELSE
                        mgrid(iG)%axis(:)=(/ grid(ixx),grid(ncells(1)+iyy),grid(ncells(1)+ncells(2)+izz)/)
                        mgrid(iG)%rho=0.
                        mgrid(iG)%nrho=0
                        mgrid(ig)%id(:)=(/ ixx,iyy,izz /)
                        mgrid(iG)%numPhots= 0
                        !WRITE(32,*) mgrid(ig)%id(:),mgrid(iG)%numPhots,mgrid(iG)%axis(:),mgrid(iG)%rho,mgrid(iG)%nrho
                    END IF
                END DO
            END DO
        END DO

    ELSE
        !not clumping
        h=0
        DO ixx=1,ncells(1)
            DO iyy=1,ncells(2)
                DO izz=1,ncells(3)
                    iG=iG+1
                    IF ((r(ig)<(R_max_cm)) .AND. (r(ig)>(R_min_cm))) THEN
                        h=h+1
                        !this was the old incorrect calculation
                        !m=m+((MD_tot_grams*((R_min_cm)/r(iG))**q)/(SF))/1.989E33
                        !rhodG=((MD_tot_grams*1.0e-14*1e-28*((R_min_cm)/r(iG))**q)/(SF*vol))
                        !corrected here (virtually no difference...)
                        rhodG=((rho_in)*(R_min_cm/r(ig))**q)
                        m=m+rhodG*vol*5.02765e8         !5.02765e8=1e42/1.989e33 i.e. conversion factor for  e42cm3 to cm3 and g to Msun
                        ndust=rhodG/SFndust
                        ndustav=ndustav+ndust
                        mgrid(iG)%axis(:)=(/ grid(ixx),grid(ncells(1)+iyy),grid(ncells(1)+ncells(2)+izz)/)
                        mgrid(iG)%rho=rhodG
                        mgrid(iG)%nrho= ndust
                        !!!ES EDIT
                        !!!watch out for the E19 ES_const in E20

                        mgrid(iG)%N_e=ES_const*(r(iG)**(-q))*1E20
                        !PRINT*,ES_const,(r(iG)**(-q)),mgrid(iG)%N_e, r(iG),R_min_cm
                        mgrid(ig)%id(:)=(/ ixx,iyy,izz /)
                        !this should be calculates as numpG and not 0, but I've left out the calculation...
                        mgrid(iG)%numPhots= 0
                        !WRITE(32,*) mgrid(ig)%id(:),mgrid(iG)%numPhots,mgrid(iG)%axis(:),mgrid(iG)%rho,mgrid(iG)%nrho
                    ELSE
                        mgrid(iG)%axis(:)=(/ grid(ixx),grid(ncells(1)+iyy),grid(ncells(1)+ncells(2)+izz)/)
                        mgrid(iG)%rho=0.
                        mgrid(iG)%nrho=0
                        mgrid(i)%N_e=0
                        mgrid(ig)%id(:)=(/ ixx,iyy,izz /)
                        mgrid(iG)%numPhots= 0
                        !WRITE(32,*) mgrid(ig)%id(:),mgrid(iG)%numPhots,mgrid(iG)%axis(:),mgrid(iG)%rho,mgrid(iG)%nrho
                    END IF
                    WRITE(32,*) mgrid(iG)%axis(:),mgrid(iG)%rho
                END DO
            END DO
        END DO
        PRINT*,'DUST GRAIN NUMBER DENSITY AT Rin',rho_in/SFndust
        PRINT*,'DUST GRAIN NUMBER DENSITY AT Rout',((rho_in)*(R_min_cm/R_max_cm)**q)/SFndust
    END IF
    PRINT*,loop
    PRINT*,'this is the actual number of clumps',test


    IF (lgClumping) THEN
        PRINT*,'number clumps test:',' using - ',ncl,'requested -',nclumps_tot
        PRINT*,'mass clumps: ','using - ',mcl,'requested: ',m_clump*nclumps_tot
        PRINT*,'mass ICM test:',' using - ',micm,'requested - ',M_icm
        PRINT*,'filling factor',ff
    ELSE
        PRINT*,'mass check (calculated as rho*V)','using - ',m,'requested -',MD_tot
    END IF
    ndustav=ndustav/h
    PRINT*,'no of grid cells inside SN',h
    PRINT*,'volume of total grid cells inside SN (e42cm)',h*vol
    PRINT*,'average dust grain density per cell (including any clumps)',ndustav
    PRINT*,''




    !CLOSE(31)
    CLOSE(32)

END SUBROUTINE construct_grid
