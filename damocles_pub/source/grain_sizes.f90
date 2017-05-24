MODULE grain_sizes

use input
USE vector_functions



IMPLICIT NONE

INTEGER             ::  i,j,k               !counters

CHARACTER(LEN=50)   ::  species_file        !name of file detailing species list
INTEGER             ::  nspecies            !number of species
REAL                ::  total_weight        !to check that total weights of species add to 1

!NOTE THAT EACH SPECIES HAS A DIFFERENT SIZE DISTRIBUTION
TYPE species                                !each species has the following attributes
    INTEGER ::  id                          !id number
    INTEGER ::  nsizes                      !number of grain sizes
    INTEGER ::  nwav                        !number of wavelengths
    REAL    ::  interval                    !spacing of grain sizes
    REAL    ::  amin,amax                   !amin, amax
    REAL    ::  weight                      !relative weight of species (fractional weighting by area)
    REAL    ::  mweight                     !relative weight of species (fractional weighting by mass)
    REAL    ::  vweight                     !relative weight of species (fraction weighting by volume)
    REAL    ::  power                       !exponent for power law size distribution
    REAL    ::  rhograin                    !density of a dust grain

    CHARACTER(LEN=50)   ::  dataFile         !data file containing optical constants for species
    REAL,ALLOCATABLE    ::  gsize(:,:)      !array containing grain sizes (1) and weights (2)
    REAL,ALLOCATABLE    ::  sca_opacity(:)  !array containing scattering extinctions at each wavelength
    REAL,ALLOCATABLE    ::  ext_opacity(:)  !array containing extinctions at each wavelength
    REAL,ALLOCATABLE    ::  g(:)            !array containing g (asymmetry factor) at each wavelength
    REAL,ALLOCATABLE    ::  wav(:)          !array containing the wavelengths
    REAL,ALLOCATABLE    ::  albedo(:)       !array containing albedos for each wavelength

END TYPE

TYPE(species),ALLOCATABLE   ::  dust(:)

contains



SUBROUTINE calculate_sizes(param_struct)

REAL    ::  norm  !normalising scale factor

REAL,INTENT(INOUT)::param_struct(pp_no)

!read species file in
OPEN(21,file = dustfile)
READ(21,*) nspecies
READ(21,*)
READ(21,*)
READ(21,*)
PRINT*,'number of species',nspecies
ALLOCATE(dust(nspecies))

total_weight=0.



DO i=1,nspecies
    READ(21,*) dust(i)%id,dust(i)%dataFile,dust(i)%weight,dust(i)%amin,dust(i)%amax,dust(i)%power,dust(i)%nsizes
    ALLOCATE(dust(i)%gsize(dust(i)%nsizes,2))
    dust(i)%gsize=0
    dust(i)%dataFile=trim(dust(i)%dataFile)
    total_weight=total_weight+dust(i)%weight
    !CHANGE HERE FOR MORE THAN 1 SPECIES
    IF (lg_pliny) THEN
        dust(i)%amin=param_struct(5)
        dust(i)%amax=param_struct(6)
        dust(i)%power=-1*param_struct(7)
    END IF
END DO

IF (total_weight/=1) THEN
    PRINT*, 'WARNING - total species weights do not add to 1'
    PRINT*, 'total weights =',total_weight
END IF

CLOSE(21)


!calculate grain sizes and relative weights
DO i=1,nspecies
    dust(i)%interval=(dust(i)%amax-dust(i)%amin)/real(dust(i)%nsizes)
    norm=0
    PRINT*,'area weight',dust(i)%weight
    dust(i)%vweight=(1.0/(1.0+(1.0/dust(i)%weight-1)**(1.5)))
        PRINT*,'volume weight',dust(i)%vweight
    DO j=1,dust(i)%nsizes
        dust(i)%gsize(j,1)=dust(i)%amin+((j-1)*dust(i)%interval)
        !PRINT*,j-1,dust(i)%gsize(j,1)
        norm=norm+(dust(i)%gsize(j,1)**dust(i)%power)
    END DO
    !normalise weights using SF
    DO j=1,dust(i)%nsizes
        dust(i)%gsize(j,2)=(dust(i)%gsize(j,1)**dust(i)%power)/norm
        !PRINT*,i,j,dust(i)%gsize(j,1),dust(i)%gsize(j,2)
    END DO
END DO


END SUBROUTINE calculate_sizes
END MODULE
