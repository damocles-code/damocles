!--------------------------------------------------------------------------!
!  this module declares line derived type object which contains properties !
!  of the line or doublet to be modelled                                   !
!--------------------------------------------------------------------------!

module class_line

    implicit none

    type line_obj
        real        ::  luminosity              !total luminosity of line       w/um
        real        ::  wavelength              !wavelength (nm) of the line being modelled
        real        ::  frequency               !frequency of the line to be modelled
        real        ::  doublet_wavelength_1    !wavelength (nm) of the first component of the doublet (if applic.)
        real        ::  doublet_wavelength_2    !wavelength (nm) of the second component of the doublet (if applic.)
        real        ::  doublet_ratio           !flux ratio between 2 components of doublet (wavelength_1/wavelength_2)
        real        ::  initial_energy          !energy of a single packet at emission
        real        ::  tot_flux               !peak flux used to scaled modelled line profile to in order to calclulate chi sq.

        integer     ::  wav_bin                 !array index of nearest wavelength bin to rest frame wavelength being modelled
        integer     ::  wav_bin_v               !array index of nearest wavelength bin to v band (547nm)

    end type

    type(line_obj)   :: line

end module class_line

