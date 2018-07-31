!-------------------------------------------------------------------------------!
!  this module declares the dust and dust species derived type objects          !
!  dust type includes properties such as exinction, mass, average grain density !
!  species types includes similar properties that are specied specific          !
!  subroutines include                                                          !
!      -  calculation of the array describing the dust grain radii              !
!      -  opacity calculations utilising mie routine                            !
!         (for each species/wavelength/grain radius combination)                !
!-------------------------------------------------------------------------------!

module class_dust

    use globals
    use class_line

    implicit none

    type species_obj                            !each species has the following attributes
        integer ::  id                          !id number
        integer ::  nsizes                      !number of grain sizes
        integer ::  n_wav                       !number of wavelengths

        real    ::  interval                    !spacing of grain sizes
        real    ::  amin,amax                   !amin, amax
        real    ::  weight                      !relative weight of species (fractional weighting by area)
        real    ::  m_weight                     !relative weight of species (fractional weighting by mass)
        real    ::  v_weight                     !relative weight of species (fraction weighting by volume)
        real    ::  power                       !exponent for power law size distribution
        real    ::  rho_grain                    !density of a dust grain
        real    ::  av_mgrain                   !average mass of a dust grain for the species

        character(len=50)   ::  datafile        !data file containing optical constants for species

        real,allocatable    ::  radius(:,:)     !array containing grain sizes (1) and weightings (2)
                                                !weightings are relative abundance by number
        real,allocatable    ::  mgrain(:)       !mass of grain for each grain size in grams
        real,allocatable    ::  c_sca(:)        !array containing scattering extinctions at each wavelength
        real,allocatable    ::  c_ext(:)        !array containing extinctions at each wavelength
        real,allocatable    ::  g_param(:)      !array containing g (asymmetry factor) at each wavelength
        real,allocatable    ::  wav(:)          !array containing the wavelengths
        real,allocatable    ::  albedo(:)       !array containing albedos for each wavelength

    end type species_obj

    type dust_obj
        integer                       ::  n_species           !number of species

        real                          ::  mass                !total mass of dust (m_sun)
        real                          ::  mass_grams          !total mass of dust (grams)
        real                          ::  m_icm               !total mass of dust in interclump medium (icm) (m_sun)
        real                          ::  m_clump             !mass of dust in a single clump (msun)
        real                          ::  lambda_ext(1)       !extinction at rest frame wavelength
        real                          ::  lambda_sca(1)       !scattering extinction at rest frame wavelength
        real                          ::  lambda_ext_v(1)     !extinction at v band wavelength (547nm)
        real                          ::  lambda_g_param(1)   !forward scattering parameter g at rest frame wavelength
        real                          ::  av_rho_grain         !average density of dust grains across all species
        real                          ::  av_mgrain           !average mass of dust grains across all species and sizes

        type(species_obj),allocatable ::  species(:)

        character(10)                 ::  scat_type           !isotropic or henyey-greenstein (hg)

    end type dust_obj

    type(dust_obj) :: dust

contains

    !this subroutine generates grain radii for each independent grain size distribution
    !and relative abundances by number for each grain size
    subroutine generate_grain_radii()

        !read species file in
        open(21,file = species_file)
        read(21,*) dust%n_species
        read(21,*)
        read(21,*)

        !write to log file
        if (.not. lg_mcmc) write(55,*) 'number of species',dust%n_species

        !allocate space for number of different dust species
        if (.not. lg_mcmc)  allocate(dust%species(dust%n_species))

        !read in properties for each species (weighting, amin, amax etc.)
        !allocate space for grain size distributions for each species
        !allocate space for mass of grain at each size for each species
        !initialise grain sizes in grain size distributions to 0
        do ii=1,dust%n_species
            if (.not. lg_mcmc) then 
               read(21,*) dust%species(ii)%id,dust%species(ii)%datafile, dust%species(ii)%weight,dust%species(ii)%amin, &
                & dust%species(ii)%amax,dust%species(ii)%power,dust%species(ii)%nsizes
            else
               read(21,*) dust%species(ii)%id,dust%species(ii)%datafile, dust%species(ii)%weight
            end if

            allocate(dust%species(ii)%radius(dust%species(ii)%nsizes,2))
            allocate(dust%species(ii)%mgrain(dust%species(ii)%nsizes))

            dust%species(ii)%radius=0
            dust%species(ii)%datafile=trim(dust%species(ii)%datafile)

            !generate grain sizes and weightings for each grain size
            !generate realtive species abundances by volume

            !calculate intervals between minimum and maximum grain radii (linear)
            dust%species(ii)%interval=(dust%species(ii)%amax-dust%species(ii)%amin)/real(dust%species(ii)%nsizes)

            !generate grain radii for grain size distribution
            !calculate sacling factor (norm) to be used to normalise abundances/weightings
            do jj=1,dust%species(ii)%nsizes
                dust%species(ii)%radius(jj,1)=dust%species(ii)%amin+((jj-1)*dust%species(ii)%interval)
            end do

            !generate weighting (relative abundance by number) for each grain radius (normalised so sum is unity)
            !!check that this weighting is actually by number of grains rather than mass or density etc.
            dust%species(ii)%radius(:,2)=(dust%species(ii)%radius(:,1)**dust%species(ii)%power)/sum(dust%species(ii)%radius(:,1)**dust%species(ii)%power)

            !!check conversion to weighting by volume for more than 2 species
            !conversion has been checked for 2 species (can be done analytically)
            if (dust%n_species > 2) then
                print*, 'you have requested more than 2 dust species - please check the volume weighting calculation. aborting'
                stop
            end if

            !calculate volume weighting (abundance) of species based on weighting by cross-sectional area
            dust%species(ii)%v_weight=(1.0/(1.0+(1.0/dust%species(ii)%weight-1)**(1.5)))
            if (.not. lg_mcmc) then
            write(55,*) 'area weight',dust%species(ii)%weight
            write(55,*) 'volume weight',dust%species(ii)%v_weight
            end if
        end do

        close(21)

        !check that the sum of the specified species weightings sums to 1
        if (sum(dust%species%weight) /= 1) then
            print*, 'warning - total species weights do not add to 1. aborted.'
            stop
        end if

    end subroutine generate_grain_radii

    !this subroutine calculates the dust extinction efficiences as a function of wavelength
    subroutine calculate_opacities()

        !optical property variables which will be used for each species
        !not included in species type as do not need to store once mie calculation performed
        real,allocatable        :: e_re(:),e_im(:)                  !imaginary and real parts of refractive index (n and k values)
        real,allocatable        :: qext(:,:),qsca(:,:),ggsca(:,:)   !exctinction/scattering efficiencies and forward scattering param
        real                    :: t_subl                            !sublimation temperature of dust
        real                    :: sizeparam                        !standard mie theory size parameter 2*pi*a/lambda
        complex                 :: refrel                           !complex version of n and k (n + ik) to be read into mie routine

        character(len=50)       :: junk                             !holder

        if (.not. lg_mcmc) print*, 'calculating opacities...'

        call generate_grain_radii()

        !write out to log file
        if (.not. lg_mcmc) then
        do ii=1,dust%n_species
            write(55,*) 'species',dust%species(ii)%id,'min grain radius',dust%species(ii)%amin
            write(55,*) 'species',dust%species(ii)%id,'max grain radius',dust%species(ii)%amax
            write(55,*) 'species',dust%species(ii)%id,'power law index for grain distriution',dust%species(ii)%power
        end do
        end if

        !read in optical data (n and k values) for each species
        do ii=1,dust%n_species
            open(13,file=PREFIX//"/share/damocles/"//trim(dust%species(ii)%datafile))
            read(13,*) dust%species(ii)%n_wav
            read(13,*)
            read(13,*) junk,t_subl,dust%species(ii)%rho_grain

            !allocate temporary space for results of mie scattering calculation
            allocate(e_re(dust%species(ii)%n_wav))
            allocate(e_im(dust%species(ii)%n_wav))
            allocate(qext(dust%species(ii)%nsizes,dust%species(ii)%n_wav))
            allocate(qsca(dust%species(ii)%nsizes,dust%species(ii)%n_wav))
            allocate(ggsca(dust%species(ii)%nsizes,dust%species(ii)%n_wav))


            !allocate permanent space for extinction efficiency parameters to be stored
            !extinction parameters are calculated for each species for each wavelength (averaged over the size distribution)
            allocate(dust%species(ii)%wav(dust%species(ii)%n_wav))
            allocate(dust%species(ii)%c_ext(dust%species(ii)%n_wav))
            allocate(dust%species(ii)%c_sca(dust%species(ii)%n_wav))
            allocate(dust%species(ii)%albedo(dust%species(ii)%n_wav))
            allocate(dust%species(ii)%g_param(dust%species(ii)%n_wav))

            !read in optical data for each species file
            do jj=1,dust%species(ii)%n_wav
                read(13,*) dust%species(ii)%wav(jj),e_re(jj),e_im(jj)
            end do

            close(13)

            !initiliase arrays to 0
            dust%species(ii)%c_ext(:)=0.
            dust%species(ii)%c_sca(:)=0.
            dust%species(ii)%g_param(:)=0.

            do jj=1,dust%species(ii)%n_wav
                do kk=1,dust%species(ii)%nsizes

                    !generate size parameter and complex refractive index for mie scattering routine
                    !routine returns q_ext,q_sca and g (forward scattering parameter) for each size and wavelength pair
                    sizeparam=2*pi*dust%species(ii)%radius(kk,1)/(dust%species(ii)%wav(jj))
                    refrel=cmplx(e_re(jj),e_im(jj))
                    call bhmie(sizeparam,refrel,qext(kk,jj),qsca(kk,jj),ggsca(kk,jj))

                    !note here that grain_rad(j,2) is the relative abundance of grain with radius a for that species
                    dust%species(ii)%c_ext(jj)=dust%species(ii)%c_ext(jj)+ &
                        & (dust%species(ii)%radius(kk,2)*qext(kk,jj)*pi*(dust%species(ii)%radius(kk,1)*1e-4)**2)

                    dust%species(ii)%c_sca(jj)=dust%species(ii)%c_sca(jj)+ &
                        & (dust%species(ii)%radius(kk,2)*qsca(kk,jj)*pi*(dust%species(ii)%radius(kk,1)*1e-4)**2)

                    !forward scattering parameter g weighted according to c_ext (c_ext determines likelihood of interaction with a particle of given size)
                    dust%species(ii)%g_param(jj)=dust%species(ii)%g_param(jj)+ &
                        & (dust%species(ii)%radius(kk,2)*ggsca(kk,jj)*qext(kk,jj)*pi*(dust%species(ii)%radius(kk,1)*1e-4)**2)

                end do

            end do

            !normalise forward scattering parameter g by diding by total c_ext
            dust%species(ii)%g_param(:)=dust%species(ii)%g_param(:)/dust%species(ii)%c_ext(:)

            !calculate albedo for each species at each wavelength
            dust%species(ii)%albedo=dust%species(ii)%c_sca/dust%species(ii)%c_ext

            !deallocate temporary space
            deallocate(e_re)
            deallocate(e_im)
            deallocate(qext)
            deallocate(qsca)
            deallocate(ggsca)

        end do

        !calculate average grain density across all species
        dust%av_rho_grain=sum(dust%species%v_weight*dust%species%rho_grain)

        !calculate mass weighting for each species
        dust%species%m_weight=dust%species%rho_grain*dust%species%v_weight/dust%av_rho_grain

        !calculate average opacity for lamba_0
        dust%lambda_ext=0
        dust%lambda_ext_v=0
        dust%lambda_g_param=0

        !for each species, calculate extinction coefficients
        do ii=1,dust%n_species

            !calculate mass of a grain for each size for each species and the average grain mass for each species
            dust%species(ii)%mgrain(:)=(4*pi*dust%species(ii)%radius(:,1)**3*dust%species(ii)%rho_grain*1e-12)/3
            dust%species(ii)%av_mgrain=sum(dust%species(ii)%mgrain(:)*dust%species(ii)%radius(:,2))


            !find neareset wavelength in array to rest frame wavelength (lambda_0) and v band (547nm)
            line%wav_bin=minloc(abs((dust%species(ii)%wav(:)-(line%wavelength/1000))),1)
            line%wav_bin_v=minloc(abs((dust%species(ii)%wav(:)-(547.0/1000))),1)


            !calculate extinction for rest frame wavelength and v band (547nm) weighted sum over all species
            !interpolate between the limits of the wavelength bin that contains the desired wavelength
            dust%lambda_ext=dust%lambda_ext+ &
                & dust%species(ii)%weight*(dust%species(ii)%c_ext(line%wav_bin)-((dust%species(ii)%c_ext(line%wav_bin)-dust%species(ii)%c_ext(line%wav_bin-1))* &
                & ((dust%species(ii)%wav(line%wav_bin)-(line%wavelength/1000))/(dust%species(ii)%wav(line%wav_bin)-dust%species(ii)%wav(line%wav_bin-1)))))

            dust%lambda_sca=dust%lambda_sca+ &
                & dust%species(ii)%weight*(dust%species(ii)%c_sca(line%wav_bin)-((dust%species(ii)%c_sca(line%wav_bin)-dust%species(ii)%c_sca(line%wav_bin-1))* &
                & ((dust%species(ii)%wav(line%wav_bin)-(line%wavelength/1000))/(dust%species(ii)%wav(line%wav_bin)-dust%species(ii)%wav(line%wav_bin-1)))))

            dust%lambda_ext_v=dust%lambda_ext_v+ &
                & dust%species(ii)%weight*(dust%species(ii)%c_ext(line%wav_bin_v)-((dust%species(ii)%c_ext(line%wav_bin_v)-dust%species(ii)%c_ext(line%wav_bin_v-1))* &
                & ((dust%species(ii)%wav(line%wav_bin_v)-(547.0/1000))/(dust%species(ii)%wav(line%wav_bin_v)-dust%species(ii)%wav(line%wav_bin_v-1)))))

            dust%lambda_g_param=dust%lambda_g_param+ &
                & dust%species(ii)%weight*(dust%species(ii)%g_param(line%wav_bin)-((dust%species(ii)%g_param(line%wav_bin)-dust%species(ii)%g_param(line%wav_bin-1))* &
                & ((dust%species(ii)%wav(line%wav_bin)-(line%wavelength/1000))/(dust%species(ii)%wav(line%wav_bin)-dust%species(ii)%wav(line%wav_bin-1)))))

        end do

        !mass*species weighting for each size for each species
        !this will be used to convert the dust mass density distribution to a number density distribution
        dust%av_mgrain=sum(dust%species%av_mgrain*dust%species%m_weight)

    end subroutine calculate_opacities

    !subroutine to check that an appropriate scattering type has been specified
    subroutine check_scat_type()
        if ((trim(dust%scat_type) /= "isotropic") .and. (trim(dust%scat_type) /= "hg")) then
            print*, "please enter a dust scattering type of 'isotropic' or 'hg' (for henyey-greenstein). aborted."
            stop
        end if
    end subroutine

end module class_dust
