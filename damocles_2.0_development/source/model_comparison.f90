module model_comparison

    USE globals
    USE input
    USE initialise


    implicit none

contains

    subroutine linear_interp(marg_chi)

        REAL                ::  xmod(nu_grid%n_bins),ymod(nu_grid%n_bins),ydatnew(nu_grid%n_bins)       !xmod are the model wavelengths, ymod are the model fluxes
        REAL,ALLOCATABLE    ::  xdat(:),ydat(:),ydatnew_rest(:),ymod_rest(:),x_rest(:)    !
        REAL                ::  norm_const                      !to normalise the fluxes
        REAL,INTENT(OUT)    ::  marg_chi
        REAL                ::  marg_chi2,chi,chi2
        REAL                ::  s,s2
        REAL                ::  err
        REAL                :: limit_min, limit_max  !limits (in x) for the narrow line so nomalisation works
        REAL                :: ylimit_min, ylimit_max
        INTEGER             ::  n_bins_dat                      !number of data point to read in
        INTEGER             ::  inu_dat(1)
        INTEGER             ::  arr_size
        INTEGER             ::  ilim(1),nn

        OPEN(25,file='output/line.out')
        WRITE(25,*) 'wavelength   modelled flux   observed flux'

        !set up and normalise model data
        DO ii=1,nu_grid%n_bins
            ymod(ii)=NP_BIN(ii)
            IF (ii /= nu_grid%n_bins) THEN
                xmod(nu_grid%n_bins-ii+1)=c*10**9/((nu_grid%bin(ii,1)+nu_grid%bin(ii+1,1))/2)   !calc lambda at centre of bin
            ELSE
                xmod(1)=c*10**9/((nu_grid%bin(ii,1)+nu_grid%fmax)/2)
            END IF
            !WRITE(22,*) inu,xmod(nu_grid%n_bins-inu+1),xmod(inu),ymod(inu)*E_0
        END DO

        arr_size=size(xmod)
        xmod = xmod(arr_size:1:-1)

        !read in observed data
        OPEN(21,file=data_file)
        READ(21,*) n_bins_dat
        READ(21,*)
!        READ(21,*) err
!        READ(21,*)
!        READ(21,*) limit_min, limit_max
!        READ(21,*)
        ALLOCATE(xdat(n_bins_dat))
        ALLOCATE(ydat(n_bins_dat))        
        
        xdat=0
        ydat=0

        DO ii=1,n_bins_dat
            READ(21,*) xdat(ii),ydat(ii)
            !if in vel space (kms-1) include this line
            xdat(ii)=((xdat(ii)*1000/c)+1)*line%wavelength
            !if in angstroms include this line
            !xdat(inu)=xdat(inu)/10
            !PRINT*,ydat(inu)
        END DO

        !PRINT*,xmod(1),xdat(1)
        !PRINT*,xmod(nu_grid%n_bins),xdat(n_bins_dat)

        ydatnew=0
        nn=0

        !interpolate
        DO ii=1,nu_grid%n_bins
            inu_dat=MINLOC((xmod(ii)-xdat),(xmod(ii)-xdat)>0)
            IF (inu_dat(1) /=0) THEN
                IF (inu_dat(1) /=n_bins_dat) THEN
                    ydatnew(ii)=ydat(inu_dat(1))+((xmod(ii)-xdat(inu_dat(1)))*(ydat(inu_dat(1)+1)-ydat(inu_dat(1)))/(xdat(inu_dat(1)+1)-xdat(inu_dat(1))))
                    !PRINT*,ydat(inu_dat),ydatnew(inu),ydat(inu_dat+1)
                    jj=ii
                    nn=nn+1
                END IF
            END IF
        END DO

        PRINT*,nn,jj,jj-nn,nu_grid%n_bins
        
        !restrict range of model output to the range of the data
        ALLOCATE(ydatnew_rest(nn))
        ALLOCATE(ymod_rest(nn))
        ALLOCATE(x_rest(nn))
        ydatnew_rest=0
        ymod_rest=0
        kk=1
        DO ii=jj-nn+1,jj
                ydatnew_rest(kk)=ydatnew(ii)
                ymod_rest(kk)=ymod(ii)
                x_rest(kk)=xmod(ii)
                kk=kk+1                
        END DO

        !normalise over restricted observed data array and restricted model array

        PRINT*,limit_min,limit_max
        ilim=minloc(abs(limit_min-x_rest))
        ylimit_min=ydatnew_rest(ilim(1))
        ilim=minloc(abs(limit_max-x_rest))
        ylimit_max=ydatnew_rest(ilim(1))
        norm_const=0
        DO ii=1,nn
              norm_const=norm_const+ymod_rest(ii)
        END DO

        IF (norm_const==0) THEN
            PRINT*,'normalisation constant =0 - maybe not enough packets used - 100% absorbed?'
        END IF
        ymod_rest=ymod_rest/(norm_const)

        norm_const=0
        DO ii=1,nn
!           IF ((x_rest(inu)<limit_min) .OR. (x_rest(inu)>limit_max)) THEN
              norm_const=norm_const+ydatnew_rest(ii)
!           ELSE
              !here we are compensating for the narrow line component
              !the narrow line adds extra flux that is not being modelled hE_Re
              !this throws off the comparison of the lines by normalising
              !instead we identify the limits of the narrow component
              !to normalise we interpolate linearly in this region
           !   PRINT*,
!              norm_const=norm_const+(ylimit_min+((ylimit_max-ylimit_min)*(x_rest(inu)-limit_min)/(limit_max-limit_min)))
              !PRINT*,'check',x_rest(inu),ylimit_min+((ylimit_max-ylimit_min)*(x_rest(inu)-limit_min)/(limit_max-limit_min))
!           END IF
           !PRINT*,ymod_rest(inu)
 
        END DO
!        norm_const=norm_const
        ydatnew_rest=ydatnew_rest/(norm_const)


        !THIS SECTION WILL CALCULATE CHI2 VALUES BUT NEEDS UPDATING TO READ IN ERRORS
        !PERFORM THIS CALCULATION USING POST-PROCESSING 


        !marginalised chi squared calculation
        !note that the following only true if err is constant
        !s=0
        !s2=0
        !DO inu=1,i
            !INTRODUCE ERROR AND CHECK FORMULA!!!!!!!!
        !    s=s+(ymod_rest(inu)*ydatnew_rest(inu))
        !    s2=s2+(ymod_rest(inu)**2)
        !    WRITE(25,*) x_rest(inu),ydatnew_rest(inu),ymod_rest(inu)
        !END DO
        !s=s/s2

        !chi2=0
        !marg_chi2=0
        !err=0.1
        !DO inu=1,i
        !    chi2=chi2+((ydatnew_rest(inu)-ymod_rest(inu))/err)**2
        !    marg_chi2=marg_chi2+((ydatnew_rest(inu)-s*ymod_rest(inu))/err)**2
        !END DO

        !marg_chi=marg_chi2
        !chi=chi2

        DEALLOCATE(ydatnew_rest)
        DEALLOCATE(ymod_rest)
        DEALLOCATE(xdat)
        DEALLOCATE(ydat)

        CLOSE(25)

    end subroutine linear_interp

end module model_comparison
