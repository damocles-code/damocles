!-------------------------------------------------------------------------------!
!  this module declares the dust and dust species derived type objects          !
!  dust type includes properties such as exinction, mass, average grain density !
!  species types includes similar properties that are specied specific          !
!  subroutines include                                                          !
!      -  calculation of the array describing the dust grain radii              !
!      -  opacity calculations utilising Mie routine                            !
!         (for each species/wavelength/grain radius combination)                !
!-------------------------------------------------------------------------------!

MODULE class_dust

    USE globals
    USE class_line

    implicit none

    TYPE species_obj                            !each species has the following attributes
        INTEGER ::  id                          !id number
        INTEGER ::  nsizes                      !number of grain sizes
        INTEGER ::  n_wav                       !number of wavelengths

        REAL    ::  interval                    !spacing of grain sizes
        REAL    ::  amin,amax                   !amin, amax
        REAL    ::  weight                      !relative weight of species (fractional weighting by area)
        REAL    ::  mweight                     !relative weight of species (fractional weighting by mass)
        REAL    ::  vweight                     !relative weight of species (fraction weighting by volume)
        REAL    ::  power                       !exponent for power law size distribution
        REAL    ::  rhograin                    !density of a dust grain
        REAL    ::  av_mgrain                   !average mass of a dust grain for the species

        CHARACTER(LEN=50)   ::  dataFile        !data file containing optical constants for species

        REAL,ALLOCATABLE    ::  radius(:,:)     !array containing grain sizes (1) and weightings (2)
                                                !weightings are relative abundance by number
        REAL,ALLOCATABLE    ::  mgrain(:)       !mass of grain for each grain size in grams
        REAL,ALLOCATABLE    ::  C_sca(:)        !array containing scattering extinctions at each wavelength
        REAL,ALLOCATABLE    ::  C_ext(:)        !array containing extinctions at each wavelength
        REAL,ALLOCATABLE    ::  g_param(:)      !array containing g (asymmetry factor) at each wavelength
        REAL,ALLOCATABLE    ::  wav(:)          !array containing the wavelengths
        REAL,ALLOCATABLE    ::  albedo(:)       !array containing albedos for each wavelength

    END TYPE species_obj

    TYPE dust_obj
        INTEGER                       ::  n_species           !number of species

        REAL                          ::  mass                !total mass of dust (M_sun)
        REAL                          ::  mass_grams          !total mass of dust (grams)
        REAL                          ::  m_icm               !total mass of dust in interclump medium (ICM) (M_sun)
        REAL                          ::  m_clump             !mass of dust in a single clump (Msun)
        REAL                          ::  lambda_ext(1)       !extinction at rest frame wavelength
        REAL                          ::  lambda_sca(1)       !scattering extinction at rest frame wavelength
        REAL                          ::  lambda_ext_V(1)     !extinction at V band wavelength (547nm)
        REAL                          ::  lambda_g_param(1)   !forward scattering parameter g at rest frame wavelength
        REAL                          ::  av_rhograin         !average density of dust grains across all species
        REAL                          ::  av_mgrain           !average mass of dust grains across all species and sizes

        TYPE(species_obj),ALLOCATABLE ::  species(:)

        CHARACTER(10)                 ::  scat_type           !isotropic or henyey-greenstein (hg)

    END TYPE dust_obj

    TYPE(dust_obj) :: dust

contains

    !This subroutine generates grain radii for each independent grain size distribution
    !and relative abundances by number for each grain size
    SUBROUTINE generate_grain_radii()

        !read species file in
        OPEN(21,file = species_file)
        READ(21,*) dust%n_species
        READ(21,*)
        READ(21,*)

        !write to log file
        WRITE(55,*) 'number of species',dust%n_species

        !allocate space for number of different dust species
        ALLOCATE(dust%species(dust%n_species))

        !read in properties for each species (weighting, amin, amax etc.)
        !allocate space for grain size distributions for each species
        !allocate space for mass of grain at each size for each species
        !initialise grain sizes in grain size distributions to 0
        DO ii=1,dust%n_species
            READ(21,*) dust%species(ii)%id,dust%species(ii)%dataFile, dust%species(ii)%weight,dust%species(ii)%amin, &
                & dust%species(ii)%amax,dust%species(ii)%power,dust%species(ii)%nsizes

            ALLOCATE(dust%species(ii)%radius(dust%species(ii)%nsizes,2))
            ALLOCATE(dust%species(ii)%mgrain(dust%species(ii)%nsizes))

            dust%species(ii)%radius=0
            dust%species(ii)%dataFile=trim(dust%species(ii)%dataFile)

            !generate grain sizes and weightings for each grain size
            !generate realtive species abundances by volume

            !calculate intervals between minimum and maximum grain radii (linear)
            dust%species(ii)%interval=(dust%species(ii)%amax-dust%species(ii)%amin)/real(dust%species(ii)%nsizes)

            !generate grain radii for grain size distribution
            !calculate sacling factor (norm) to be used to normalise abundances/weightings
            DO jj=1,dust%species(ii)%nsizes
                dust%species(ii)%radius(jj,1)=dust%species(ii)%amin+((jj-1)*dust%species(ii)%interval)
            END DO

            !generate weighting (relative abundance by number) for each grain radius (normalised so sum is unity)
            !!check that this weighting is actually by number of grains rather than mass or density etc.
            dust%species(ii)%radius(:,2)=(dust%species(ii)%radius(:,1)**dust%species(ii)%power)/sum(dust%species(ii)%radius(:,1)**dust%species(ii)%power)

            !!check conversion to weighting by volume for more than 2 species
            !conversion has been checked for 2 species (can be done analytically)
            IF (dust%n_species > 2) THEN
                PRINT*, 'You have requested more than 2 dust species - please check the volume weighting calculation. Aborting'
                STOP
            END IF

            !calculate volume weighting (abundance) of species based on weighting by cross-sectional area
            dust%species(ii)%vweight=(1.0/(1.0+(1.0/dust%species(ii)%weight-1)**(1.5)))

            WRITE(55,*) 'area weight',dust%species(ii)%weight
            WRITE(55,*),'volume weight',dust%species(ii)%vweight
        END DO

        CLOSE(21)

        !check that the sum of the specified species weightings sums to 1
        IF (sum(dust%species%weight) /= 1) THEN
            PRINT*, 'WARNING - total species weights do not add to 1. Aborted.'
            STOP
        END IF

    END SUBROUTINE generate_grain_radii

    !This subroutine calculates the dust extinction efficiences as a function of wavelength
    SUBROUTINE calculate_opacities()

        !optical property variables which will be used for each species
        !not included in species type as do not need to store once Mie calculation performed
        REAL,ALLOCATABLE        :: E_Re(:),E_Im(:)                  !imaginary and real parts of refractive index (n and k values)
        REAL,ALLOCATABLE        :: Qext(:,:),Qsca(:,:),ggsca(:,:)   !exctinction/scattering efficiencies and forward scattering param
        REAL                    :: T_subl                            !sublimation temperature of dust
        REAL                    :: sizeparam                        !standard Mie theory size parameter 2*pi*a/lambda
        COMPLEX                 :: refrel                           !complex version of n and k (n + ik) to be read into Mie routine

        CHARACTER(LEN=50)       :: junk                             !holder

        PRINT*, 'Calculating opacities...'

        call generate_grain_radii()

        !write out to log file
        DO ii=1,dust%n_species
            WRITE(55,*) 'species',dust%species(ii)%id,'min grain radius',dust%species(ii)%amin
            WRITE(55,*) 'species',dust%species(ii)%id,'max grain radius',dust%species(ii)%amax
            WRITE(55,*) 'species',dust%species(ii)%id,'power law index for grain distriution',dust%species(ii)%power
        END DO

        !read in optical data (n and k values) for each species
        DO ii=1,dust%n_species
            OPEN(13,file=PREFIX//"/share/damocles/"//trim(dust%species(ii)%dataFile))
            READ(13,*) dust%species(ii)%n_wav
            READ(13,*)
            READ(13,*) junk,T_subl,dust%species(ii)%rhograin

            !allocate temporary space for results of Mie scattering calculation
            ALLOCATE(E_Re(dust%species(ii)%n_wav))
            ALLOCATE(E_Im(dust%species(ii)%n_wav))
            ALLOCATE(Qext(dust%species(ii)%nsizes,dust%species(ii)%n_wav))
            ALLOCATE(Qsca(dust%species(ii)%nsizes,dust%species(ii)%n_wav))
            ALLOCATE(ggsca(dust%species(ii)%nsizes,dust%species(ii)%n_wav))


            !allocate permanent space for extinction efficiency parameters to be stored
            !extinction parameters are calculated for each species for each wavelength (averaged over the size distribution)
            ALLOCATE(dust%species(ii)%wav(dust%species(ii)%n_wav))
            ALLOCATE(dust%species(ii)%C_ext(dust%species(ii)%n_wav))
            ALLOCATE(dust%species(ii)%C_sca(dust%species(ii)%n_wav))
            ALLOCATE(dust%species(ii)%albedo(dust%species(ii)%n_wav))
            ALLOCATE(dust%species(ii)%g_param(dust%species(ii)%n_wav))

            !read in optical data for each species file
            DO jj=1,dust%species(ii)%n_wav
                READ(13,*) dust%species(ii)%wav(jj),E_Re(jj),E_Im(jj)
            END DO

            CLOSE(13)

            !initiliase arrays to 0
            dust%species(ii)%C_ext(:)=0.
            dust%species(ii)%C_sca(:)=0.
            dust%species(ii)%g_param(:)=0.

            DO jj=1,dust%species(ii)%n_wav
                DO kk=1,dust%species(ii)%nsizes

                    !generate size parameter and complex refractive index for mie scattering routine
                    !routine returns Q_ext,Q_sca and g (forward scattering parameter) for each size and wavelength pair
                    sizeparam=2*pi*dust%species(ii)%radius(kk,1)/(dust%species(ii)%wav(jj))
                    refrel=cmplx(E_Re(jj),E_Im(jj))
                    call BHmie(sizeparam,refrel,Qext(kk,jj),Qsca(kk,jj),ggsca(kk,jj))

                    !note here that grain_rad(j,2) is the relative abundance of grain with radius a for that species
                    dust%species(ii)%C_ext(jj)=dust%species(ii)%C_ext(jj)+ &
                        & (dust%species(ii)%radius(kk,2)*Qext(kk,jj)*pi*(dust%species(ii)%radius(kk,1)*1e-4)**2)

                    dust%species(ii)%C_sca(jj)=dust%species(ii)%C_sca(jj)+ &
                        & (dust%species(ii)%radius(kk,2)*Qsca(kk,jj)*pi*(dust%species(ii)%radius(kk,1)*1e-4)**2)

                    !forward scattering parameter g weighted according to C_ext (C_ext determines likelihood of interaction with a particle of given size)
                    dust%species(ii)%g_param(jj)=dust%species(ii)%g_param(jj)+ &
                        & (dust%species(ii)%radius(kk,2)*ggsca(kk,jj)*Qext(kk,jj)*pi*(dust%species(ii)%radius(kk,1)*1e-4)**2)

                END DO

            END DO

            !normalise forward scattering parameter g by diding by total C_ext
            dust%species(ii)%g_param(:)=dust%species(ii)%g_param(:)/dust%species(ii)%C_ext(:)

            !calculate albedo for each species at each wavelength
            dust%species(ii)%albedo=dust%species(ii)%C_sca/dust%species(ii)%C_ext

            !deallocate temporary space
            DEALLOCATE(E_Re)
            DEALLOCATE(E_Im)
            DEALLOCATE(Qext)
            DEALLOCATE(Qsca)
            DEALLOCATE(ggsca)

        END DO

        !calculate average grain density across all species
        dust%av_rhograin=sum(dust%species%vweight*dust%species%rhograin)

        !calculate mass weighting for each species
        dust%species%mweight=dust%species%rhograin*dust%species%vweight/dust%av_rhograin

        !calculate average opacity for lamba_0
        dust%lambda_ext=0
        dust%lambda_ext_V=0
        dust%lambda_g_param=0

        !for each species, calculate extinction coefficients
        !!work out what these quanities actually are... extinction per unit mass? per unit csa?
        DO ii=1,dust%n_species

            !calculate mass of a grain for each size for each species and the average grain mass for each species
            dust%species(ii)%mgrain(:)=(4*pi*dust%species(ii)%radius(:,1)**3*dust%species(ii)%rhograin*1e-12)/3
            dust%species(ii)%av_mgrain=sum(dust%species(ii)%mgrain(:)*dust%species(ii)%radius(:,2))


            !find neareset wavelength in array to rest frame wavelength (lambda_0) and V band (547nm)
            line%wav_bin=MINLOC(ABS((dust%species(ii)%wav(:)-(line%wavelength/1000))),1)
            line%wav_bin_v=MINLOC(ABS((dust%species(ii)%wav(:)-(547.0/1000))),1)


            !calculate extinction for rest frame wavelength and V band (547nm) weighted sum over all species
            !interpolate between the limits of the wavelength bin that contains the desired wavelength
            dust%lambda_ext=dust%lambda_ext+ &
                & dust%species(ii)%weight*(dust%species(ii)%C_ext(line%wav_bin)-((dust%species(ii)%C_ext(line%wav_bin)-dust%species(ii)%C_ext(line%wav_bin-1))* &
                & ((dust%species(ii)%wav(line%wav_bin)-(line%wavelength/1000))/(dust%species(ii)%wav(line%wav_bin)-dust%species(ii)%wav(line%wav_bin-1)))))

            dust%lambda_sca=dust%lambda_sca+ &
                & dust%species(ii)%weight*(dust%species(ii)%C_sca(line%wav_bin)-((dust%species(ii)%C_sca(line%wav_bin)-dust%species(ii)%C_sca(line%wav_bin-1))* &
                & ((dust%species(ii)%wav(line%wav_bin)-(line%wavelength/1000))/(dust%species(ii)%wav(line%wav_bin)-dust%species(ii)%wav(line%wav_bin-1)))))

            dust%lambda_ext_V=dust%lambda_ext_V+ &
                & dust%species(ii)%weight*(dust%species(ii)%C_ext(line%wav_bin_v)-((dust%species(ii)%C_ext(line%wav_bin_v)-dust%species(ii)%C_ext(line%wav_bin_v-1))* &
                & ((dust%species(ii)%wav(line%wav_bin_v)-(547.0/1000))/(dust%species(ii)%wav(line%wav_bin_v)-dust%species(ii)%wav(line%wav_bin_v-1)))))

            dust%lambda_g_param=dust%lambda_g_param+ &
                & dust%species(ii)%weight*(dust%species(ii)%g_param(line%wav_bin)-((dust%species(ii)%g_param(line%wav_bin)-dust%species(ii)%g_param(line%wav_bin-1))* &
                & ((dust%species(ii)%wav(line%wav_bin)-(line%wavelength/1000))/(dust%species(ii)%wav(line%wav_bin)-dust%species(ii)%wav(line%wav_bin-1)))))

        END DO

        !mass*species weighting for each size for each species
        !this will be used to convert the dust mass density distribution to a number density distribution
        dust%av_mgrain=SUM(dust%species%av_mgrain*dust%species%mweight)

    END SUBROUTINE calculate_opacities

    !subroutine to check that an appropriate scattering type has been specified
    SUBROUTINE check_scat_type()
        IF ((trim(dust%scat_type) /= "isotropic") .and. (trim(dust%scat_type) /= "hg")) THEN
            PRINT*, "Please enter a dust scattering type of 'isotropic' or 'hg' (for henyey-greenstein). Aborted."
            STOP
        END IF
    END SUBROUTINE

END MODULE class_dust
