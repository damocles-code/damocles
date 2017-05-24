module model_comparison

    USE input
    USE initialise


    implicit none

contains

    subroutine linear_interp(marg_chi)

        REAL                ::  xmod(n_bins),ymod(n_bins),ydatnew(n_bins)       !xmod are the model wavelengths, ymod are the model fluxes
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
        INTEGER             ::  ilim(1)

        OPEN(25,file='line.out')
        WRITE(25,*) 'wavelength   modelled flux   observed flux'

        !set up and normalise model data
        DO inu=1,n_bins
            ymod(inu)=NP_BIN(inu)
            IF (inu /= n_bins) THEN
                xmod(n_bins-inu+1)=c*10**9/((nu_bin(inu,1)+nu_bin(inu+1,1))/2)   !calc lambda at centre of bin
            ELSE
                xmod(1)=c*10**9/((nu_bin(inu,1)+nu_max)/2)
            END IF
            !WRITE(22,*) inu,xmod(n_bins-inu+1),xmod(inu),ymod(inu)*E_0
        END DO

        arr_size=size(xmod)
        xmod = xmod(arr_size:1:-1)

        !read in observed data
        OPEN(21,file='line.in')
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

        DO inu=1,n_bins_dat
            READ(21,*) xdat(inu),ydat(inu)
            !if in vel space (kms-1) include this line
            xdat(inu)=((xdat(inu)*1000/c)+1)*lambda_0
            !if in angstroms include this line
            !xdat(inu)=xdat(inu)/10
            !PRINT*,ydat(inu)
        END DO

        !PRINT*,xmod(1),xdat(1)
        !PRINT*,xmod(n_bins),xdat(n_bins_dat)

        ydatnew=0
        i=0

        !interpolate
        DO inu=1,n_bins
            inu_dat=MINLOC((xmod(inu)-xdat),(xmod(inu)-xdat)>0)
            IF (inu_dat(1) /=0) THEN
                IF (inu_dat(1) /=n_bins_dat) THEN
                    ydatnew(inu)=ydat(inu_dat(1))+((xmod(inu)-xdat(inu_dat(1)))*(ydat(inu_dat(1)+1)-ydat(inu_dat(1)))/(xdat(inu_dat(1)+1)-xdat(inu_dat(1))))
                    !PRINT*,ydat(inu_dat),ydatnew(inu),ydat(inu_dat+1)
                    j=inu
                    i=i+1
                END IF
            END IF
        END DO

        PRINT*,i,j,j-i,n_bins
        
        !restrict range of model output to the range of the data
        ALLOCATE(ydatnew_rest(i))
        ALLOCATE(ymod_rest(i))
        ALLOCATE(x_rest(i))
        ydatnew_rest=0
        ymod_rest=0
        k=1
        DO inu=j-i+1,j
                ydatnew_rest(k)=ydatnew(inu)
                ymod_rest(k)=ymod(inu)
                x_rest(k)=xmod(inu)
                k=k+1                
        END DO

        !normalise over restricted observed data array and restricted model array

        PRINT*,limit_min,limit_max
        ilim=minloc(abs(limit_min-x_rest))
        ylimit_min=ydatnew_rest(ilim(1))
        ilim=minloc(abs(limit_max-x_rest))
        ylimit_max=ydatnew_rest(ilim(1))
        norm_const=0
        DO inu=1,i
              norm_const=norm_const+ymod_rest(inu)
        END DO

        IF (norm_const==0) THEN
            PRINT*,'normalisation constant =0 - maybe not enough packets used - 100% absorbed?'
        END IF
        ymod_rest=ymod_rest/(norm_const)

        norm_const=0
        DO inu=1,i
!           IF ((x_rest(inu)<limit_min) .OR. (x_rest(inu)>limit_max)) THEN
              norm_const=norm_const+ydatnew_rest(inu)
!           ELSE
              !here we are compensating for the narrow line component
              !the narrow line adds extra flux that is not being modelled here
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

        !marginalised chi squared calculation
        !note that the following only true if err is constant
        s=0
        s2=0
        DO inu=1,i
            !INTRODUCE ERROR AND CHECK FORMULA!!!!!!!!
            s=s+(ymod_rest(inu)*ydatnew_rest(inu))
            s2=s2+(ymod_rest(inu)**2)
            WRITE(25,*) x_rest(inu),ydatnew_rest(inu),ymod_rest(inu)
        END DO
        s=s/s2

        chi2=0
        marg_chi2=0
        err=0.1
       DO inu=1,i
            chi2=chi2+((ydatnew_rest(inu)-ymod_rest(inu))/err)**2
            marg_chi2=marg_chi2+((ydatnew_rest(inu)-s*ymod_rest(inu))/err)**2
        END DO

        PRINT*,'chi',chi,'marg chi',marg_chi

        marg_chi=marg_chi2
        chi=chi2

        DEALLOCATE(ydatnew_rest)
        DEALLOCATE(ymod_rest)
        DEALLOCATE(xdat)
        DEALLOCATE(ydat)

        CLOSE(25)

    end subroutine linear_interp

end module model_comparison
