module model_comparison

    use globals
    use input
    use initialise


    implicit none

    real :: chi2

contains

    subroutine linear_interp(marg_chi)

        real                ::  xmod(nu_grid%n_bins),ymod(nu_grid%n_bins),ydatnew(nu_grid%n_bins)       !xmod are the model wavelengths, ymod are the model fluxes
        real,allocatable    ::  xdat(:),ydat(:),ydatnew_rest(:),ymod_rest(:),x_rest(:)    !
        real                ::  norm_const                      !to normalise the fluxes
        real,intent(out)    ::  marg_chi
        real                ::  marg_chi2,chi,chi2
        real                ::  s,s2
        real                ::  err
        real                :: limit_min, limit_max  !limits (in x) for the narrow line so nomalisation works
        real                :: ylimit_min, ylimit_max
        integer             ::  n_bins_dat                      !number of data point to read in
        integer             ::  inu_dat(1)
        integer             ::  arr_size
        integer             ::  ilim(1),nn

        open(25,file='output/line.out')
        write(25,*) 'wavelength   modelled flux   observed flux'

        !set up and normalise model data
        do ii=1,nu_grid%n_bins
            ymod(ii)=profile_array(ii)
            if (ii /= nu_grid%n_bins) then
                xmod(nu_grid%n_bins-ii+1)=c*10**9/((nu_grid%bin(ii,1)+nu_grid%bin(ii+1,1))/2)   !calc lambda at centre of bin
            else
                xmod(1)=c*10**9/((nu_grid%bin(ii,1)+nu_grid%fmax)/2)
            end if
            !write(22,*) inu,xmod(nu_grid%n_bins-inu+1),xmod(inu),ymod(inu)*e_0
        end do

        arr_size=size(xmod)
        xmod = xmod(arr_size:1:-1)

        !read in observed data
        open(21,file=data_file)
        read(21,*) n_bins_dat

!        read(21,*) err
!        read(21,*)
!        read(21,*) limit_min, limit_max
!        read(21,*)
        allocate(xdat(n_bins_dat))
        allocate(ydat(n_bins_dat))
        
        xdat=0
        ydat=0

        do ii=1,n_bins_dat
            read(21,*) xdat(ii),ydat(ii)
            !if in vel space (kms-1) include this line
            !xdat(ii)=((xdat(ii)*1000/c)+1)*line%wavelength
            !if in angstroms include this line
            !xdat(inu)=xdat(inu)/10
            !if in nm do nothing
            !print*,ydat(inu)
        end do

        !print*,xmod(1),xdat(1)
        !print*,xmod(nu_grid%n_bins),xdat(n_bins_dat)

        ydatnew=0
        nn=0

        !interpolate
        do ii=1,nu_grid%n_bins
            inu_dat=minloc((xmod(ii)-xdat),(xmod(ii)-xdat)>0)
            if (inu_dat(1) /=0) then
                if (inu_dat(1) /=n_bins_dat) then
                    ydatnew(ii)=ydat(inu_dat(1))+((xmod(ii)-xdat(inu_dat(1)))*(ydat(inu_dat(1)+1)-ydat(inu_dat(1)))/(xdat(inu_dat(1)+1)-xdat(inu_dat(1))))
                    !print*,ydat(inu_dat),ydatnew(inu),ydat(inu_dat+1)
                    jj=ii
                    nn=nn+1
                end if
            end if
        end do

        !print*,nn,jj,jj-nn,nu_grid%n_bins
        
        !restrict range of model output to the range of the data
        allocate(ydatnew_rest(nn))
        allocate(ymod_rest(nn))
        allocate(x_rest(nn))
        ydatnew_rest=0
        ymod_rest=0
        kk=1
        do ii=jj-nn+1,jj
                ydatnew_rest(kk)=ydatnew(ii)
                ymod_rest(kk)=ymod(ii)
                x_rest(kk)=xmod(ii)
                kk=kk+1                
        end do

        !normalise over restricted observed data array and restricted model array

        !print*,limit_min,limit_max
        ilim=minloc(abs(limit_min-x_rest))
        ylimit_min=ydatnew_rest(ilim(1))
        ilim=minloc(abs(limit_max-x_rest))
        ylimit_max=ydatnew_rest(ilim(1))
        norm_const=0
        do ii=1,nn
              norm_const=norm_const+ymod_rest(ii)
        end do

        if (norm_const==0) then
            print*,'normalisation constant =0 - maybe not enough packets used - 100% absorbed?'
        end if
        ymod_rest=ymod_rest/(norm_const)

        norm_const=0
        do ii=1,nn
!           if ((x_rest(inu)<limit_min) .or. (x_rest(inu)>limit_max)) then
              norm_const=norm_const+ydatnew_rest(ii)
!           else
              !here we are compensating for the narrow line component
              !the narrow line adds extra flux that is not being modelled he_re
              !this throws off the comparison of the lines by normalising
              !instead we identify the limits of the narrow component
              !to normalise we interpolate linearly in this region
           !   print*,
!              norm_const=norm_const+(ylimit_min+((ylimit_max-ylimit_min)*(x_rest(inu)-limit_min)/(limit_max-limit_min)))
              !print*,'check',x_rest(inu),ylimit_min+((ylimit_max-ylimit_min)*(x_rest(inu)-limit_min)/(limit_max-limit_min))
!           end if
           !print*,ymod_rest(inu)
 
        end do
!        norm_const=norm_const
        ydatnew_rest=ydatnew_rest/(norm_const)


        !this section will calculate chi2 values but needs updating to read in errors
        !perform this calculation using post-processing


        !marginalised chi squared calculation
        !note that the following only true if err is constant
        !s=0
        !s2=0
        !do inu=1,i
            !introduce error and check formula!!!!!!!!
        !    s=s+(ymod_rest(inu)*ydatnew_rest(inu))
        !    s2=s2+(ymod_rest(inu)**2)
        !    write(25,*) x_rest(inu),ydatnew_rest(inu),ymod_rest(inu)
        !end do
        !s=s/s2

        !chi2=0
        !marg_chi2=0
        !err=0.1
        !do inu=1,i
        !    chi2=chi2+((ydatnew_rest(inu)-ymod_rest(inu))/err)**2
        !    marg_chi2=marg_chi2+((ydatnew_rest(inu)-s*ymod_rest(inu))/err)**2
        !end do

        !marg_chi=marg_chi2
        !chi=chi2

        deallocate(ydatnew_rest)
        deallocate(ymod_rest)
        deallocate(xdat)
        deallocate(ydat)

        close(25)

    end subroutine linear_interp

end module model_comparison
