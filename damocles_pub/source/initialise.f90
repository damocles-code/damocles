MODULE initialise
    
    USE input
    USE grain_sizes

    IMPLICIT NONE


    REAL    ::  totvol          !total volume of supernova in 1e42cm^3
    REAL    ::  vol


    INTEGER(8) 			    ::  	iP,idP,n,iGP,idPP,n_inactive,inu,freq(1),nsize, &
        & ios,nabs,iDoublet
    INTEGER					::	iSh,n_threads,id(1),id_V(1)
    REAL					::	width(3),E_0,nu_0,nu_max,nu_min, &
        & SCAT_RF_PT(2),nu_bin_width,lambda_bin,vel_bin, &
        & sizeparam, Tsubl,tot,shell_width, alb, &
        & const,csa,dummy,lambda_ext(1),lambda_sca(1),lambda_ext_V(1)
    REAL,DIMENSION(:,:),ALLOCATABLE		::	nu_bin,tmp, &
        & grain_rad,Qext,Qsca, &
        & ggsca,grain_id, &
        & RSh,ab
    INTEGER :: ncl
    REAL,DIMENSION(:),ALLOCATABLE :: grid,mgrain
											
    REAL 					::	SCAT_RF_SN(2), &
        & SCAT_FIN(2),POS_SPH1(3),R_P,ndustav

    INTEGER(8),DIMENSION(:,:),ALLOCATABLE	::	ix,iy,iz,NP
    REAL,DIMENSION(:),ALLOCATABLE      :: NP_BIN,Ere,Eim,diff

    CHARACTER(LEN=1024) 			:: 	filename,junk

    COMPLEX					::	refrel

    TYPE SUBGRIDCELL
        REAL    ::  axis(3),rho,nrho
        INTEGER ::  id(3), cellStatus
    END TYPE

    TYPE GRIDCELL
        REAL    ::  axis(3),rho,nrho,subaxes(3,2),N_e
        INTEGER ::  id(3), cellStatus,sgStatus,numPhots
        TYPE(SUBGRIDCELL)   ::  subgrid(8)        !currently using subgrids of only 8 cells (1 division in each direction) - MAY CHANGE
    END TYPE

    TYPE(GRIDCELL),POINTER  :: mgrid(:)



contains

    SUBROUTINE set_grid_globals(param_struct)
      REAL,INTENT(INOUT)::param_struct(pp_no)
        PRINT*,"Total number of cells",totcells

        ALLOCATE(mgrid(totcells))
        ALLOCATE(nu_bin(n_bins,2))
        ALLOCATE(tmp(n_bins,1))
        ALLOCATE(grid(ncells(1)+ncells(2)+ncells(3)))
        ALLOCATE(NP(totcells,1))
        ALLOCATE(NP_BIN(n_bins))
        ALLOCATE(RSh(n_shells+1,2))

        
        nu_bin=0
        tmp=0
        grid=0
        NP=0
        NP_BIN=0
        RSh=0

        !CONSTRUCT DYNAMIC FILENAMES FOR EASY REFERENCE
        !WRITE(filename,'(A,E7.1E2,A,F4.2,A,F4.2,A,I0,A)') 'output/',MD_tot,'_rho',q,'_v',l,'_vs',VelShift,'.out'
        


        !CALCULATE OPACITIES
        CALL calculate_opacities(param_struct)
       

        !CONSTRUCT GRID WITH NUMBER OF PHOTONS EMITTED PER CELL
        CALL construct_grid
        

        !CALCULATE FREQUENCY OF EMITTED LINE IN Hz
        nu_0=c*10**9/lambda_0
        PRINT*,'LAMBDA_0',lambda_0
        !CALCULATE AVERAGE OPTICAL DEPTH FROM Rin to Rout
        PRINT*,'EXTINCTION FOR LAMBDA_0 ',lambda_ext
        PRINT*,'SCATTERING FOR LAMBDA_0',lambda_sca
        PRINT*,'ALBEDO FOR LAMBDA_0',lambda_sca/lambda_ext
        PRINT*,'AVERAGE OPTICAL DEPTH IN LAMBDA_0',ndustav*lambda_ext*(R_max_cm-R_min_cm)
        PRINT*,'AVERAGE OPTICAL DEPTH IN V',ndustav*lambda_ext_V*(R_max_cm-R_min_cm)
        PRINT*,'AVERAGE NUMBER DENSITY',ndustav
        
        PRINT*,'EFFECTIVE SPHERICAL RADIUS OF CLUMP (cm)',width(1)*(3.0/(4.0*pi))**0.3333333
        PRINT*,'EFFECTIVE SPHERICAL RADIUS AS A FRACTION OF Rout',(width(1)/R_max_cm)*(3.0/(4.0*pi))**0.3333333
        PRINT*,'AVERAGE OPTICAL DEPTH OF A CELL IN LAMBDA_0',ndustav*lambda_ext*0.5*width(1)*(3.0/(4.0*pi))**0.3333333
        PRINT*,'AVERAGE OPTICAL DEPTH OF A CELL IN LAMBDA_V',ndustav*lambda_ext_V*0.5*width(1)*(3.0/(4.0*pi))**0.3333333
!CALCULATE FREQUENCY BINS ARRAY
	
	    !MAX element here to account for ES velocities

	    IF (lgDoublet) THEN
	        IF (lambda2_0>lambda1_0) THEN
	            nu_max=5*((c*10**9/lambda1_0)/(1-MAX(v_max,2000.0)*10**3/c))        !arbitrary factor of 5 to compensate for multiple scatterings
	            nu_min=0.2*((c*10**9/lambda2_0)/(1+MAX(v_max,2000.0)*10**3/c))        !as above
	        ELSE
                nu_max=5*((c*10**9/lambda2_0)/(1-MAX(v_max,2000.0)*10**3/c))        !arbitrary factor of 5 to compensate for multiple scatterings
                nu_min=0.2*((c*10**9/lambda1_0)/(1+MAX(v_max,2000.0)*10**3/c))        !as above
            END IF
	    ELSE
            nu_max=5*(nu_0/(1-MAX(v_max,2000.0)*10**3/c))		!arbitrary factor of 4 to compensate for multiple scatterings
										                !shifting nu beyond maxnu for 1 scattering event
            nu_min=0.2*(nu_0/(1+MAX(v_max,2000.0)*10**3/c))		!as above
	    END IF
        nu_bin_width=(nu_max-nu_min)/n_bins

        !CALCULATE FREQUENCY BINS AND WRITE TO FILE
        IF (iDoublet==1) THEN
        DO inu=1,n_bins
            nu_bin(inu,1)=nu_min+((inu-1)*nu_bin_width)
            nu_bin(inu,2)=nu_min+((inu)*nu_bin_width)
        END DO
        END IF
    END SUBROUTINE

    SUBROUTINE calculate_opacities(param_struct)
      REAL,INTENT(INOUT)::param_struct(pp_no)
      REAL  ::  av_rhograin
        !CALCULATE Qext FOR EACH GRAIN SIZE AND WAVELENGTH


        !OPEN(unit=13,file=trim(dustDataFile))
        call calculate_sizes(param_struct)
        PRINT*, 'PLINY parameters in use - set 2'
        PRINT*,'amin',dust(1)%amin
        PRINT*,'amax',dust(1)%amax
        PRINT*,'aslope',dust(1)%power

       av_rhograin=0

        DO i=1,nspecies
             OPEN(13,file=trim(dust(i)%dataFile))


            READ(13,*) dust(i)%nwav
            READ(13,*)
            READ(13,*) junk,Tsubl,dust(i)%rhograin

            
            ALLOCATE(dust(i)%wav(dust(i)%nwav))
            ALLOCATE(Ere(dust(i)%nwav))
            ALLOCATE(Eim(dust(i)%nwav))

            PRINT*,''

            DO j=1,dust(i)%nwav
                READ(13,*) dust(i)%wav(j),Ere(j),Eim(j)
            END DO
            !calculate average grain density across all species
            av_rhograin=av_rhograin+dust(i)%vweight*dust(i)%rhograin

            CLOSE(13)

            csa=0

            ALLOCATE(mgrain(dust(i)%nsizes))

            DO j=1,dust(i)%nsizes
                mgrain(j)=(4*pi*dust(i)%gsize(j,1)**3*dust(i)%rhograin*1e-12)/3		!in grams
                csa=csa+(dust(i)%gsize(j,2)*pi*(dust(i)%gsize(j,1)*1e-4)**2)
            END DO


            ALLOCATE(Qext(dust(i)%nsizes,dust(i)%nwav))
            ALLOCATE(Qsca(dust(i)%nsizes,dust(i)%nwav))
            ALLOCATE(ggsca(dust(i)%nsizes,dust(i)%nwav))
            ALLOCATE(dust(i)%ext_opacity(dust(i)%nwav))
            ALLOCATE(dust(i)%sca_opacity(dust(i)%nwav))
            ALLOCATE(dust(i)%albedo(dust(i)%nwav))
            ALLOCATE(dust(i)%g(dust(i)%nwav))

            dust(i)%ext_opacity(:)=0.
            dust(i)%sca_opacity(:)=0.
            dust(i)%g(:)=0.
            OPEN(unit=24,file='opacity_wav.out')
            WRITE(24,*) 'species no - wav - extinction - scatter - g'
            !OPEN(unit=57,file='opacity_size.out')
            !WRITE(57,*) 'species no - wav - size - Qext - Qsca'
            OPEN(57,file='opacity_size.out')

            DO j=1,dust(i)%nwav
            alb=0
 !              PRINT*,j,dust(i)%wav(j)
               DO k=1,dust(i)%nsizes
                    !PRINT*,dust(i)%gsize(k,1),dust(i)%gsize(k,2)
                    sizeparam=2*pi*dust(i)%gsize(k,1)/(dust(i)%wav(j))
                    
                    refrel=cmplx(Ere(j),Eim(j))
                    
                    call BHmie(sizeparam,refrel,Qext(k,j),Qsca(k,j),ggsca(k,j))
                    !PRINT*,Qext(k,j),Qsca(k,j),dust(i)%wav(j),dust(i)%gsize(k,1)
                    alb=alb+(dust(i)%gsize(k,2)*(Qsca(k,j)/Qext(k,j))*pi*(dust(i)%gsize(k,1)*1e-4)**2)
                    dust(i)%ext_opacity(j)=dust(i)%ext_opacity(j)+(dust(i)%gsize(k,2)*Qext(k,j)*pi*(dust(i)%gsize(k,1)*1e-4)**2)		!NOTE here that grain_rad(j,2) is the relative abundance of grain with radius a
                    dust(i)%sca_opacity(j)=dust(i)%sca_opacity(j)+(dust(i)%gsize(k,2)*Qsca(k,j)*pi*(dust(i)%gsize(k,1)*1e-4)**2)
                    dust(i)%g(j)=dust(i)%g(j)+(dust(i)%gsize(k,2)*ggsca(k,j)*pi*(dust(i)%gsize(k,1)*1e-4)**2)
                    
                    !IF (j==34) THEN
!                       PRINT*,'this is working'
                    !  WRITE(57,*) i,j,k,dust(i)%wav(j),dust(i)%gsize(k,1),Qext(k,j),Qsca(k,j),Qsca(k,j)/Qext(k,j)
                    !END IF
                   !IF (j==2982) THEN
                    !   PRINT*,'writing'
                    !   WRITE(24,*) dust(i)%wav(j),dust(i)%gsize(k,1),dust(i)%ext_opacity(j),dust(i)%sca_opacity(j),dust(i)%g(j)
                    !END IF

                 END DO

                dust(i)%albedo(j)=dust(i)%sca_opacity(j)/dust(i)%ext_opacity(j)
                
                !PRINT*,j,dust(i)%wav(j)

            END DO

            CLOSE(24)
            CLOSE(57)
            DEALLOCATE(Ere)
            DEALLOCATE(Eim)
            DEALLOCATE(mgrain)
            DEALLOCATE(Qext)
            DEALLOCATE(Qsca)
            DEALLOCATE(ggsca)

        END DO


        !calculate average opacity for lamba_0
        lambda_ext=0
        lambda_ext_V=0
        DO j=1,nspecies

           dust(j)%mweight=dust(j)%rhograin*dust(j)%vweight/av_rhograin
           PRINT*,dust(j)%rhograin
           PRINT*,'mass weight',dust(j)%mweight
           !find neareset wavelength to lambda_0
           
           id=MINLOC(ABS((dust(j)%wav(:)-(lambda_0/1000))))
           id_V=MINLOC(ABS((dust(j)%wav(:)-(547.0/1000))))
           PRINT*,'id check',id,id_V
           PRINT*,'For species no',j,'albedo',dust(j)%sca_opacity(id)/dust(j)%ext_opacity(id),'weight',dust(j)%weight
           !calculate extinction for lambda_0 weighted sum over all species
           lambda_ext=lambda_ext+dust(j)%weight*(dust(j)%ext_opacity(id)-((dust(j)%ext_opacity(id)-dust(j)%ext_opacity(id-1))* &
                & ((dust(j)%wav(id)-(lambda_0/1000))/(dust(j)%wav(id)-dust(j)%wav(id-1)))))
           lambda_sca=lambda_sca+dust(j)%weight*(dust(j)%sca_opacity(id)-((dust(j)%sca_opacity(id)-dust(j)%sca_opacity(id-1))* &
                & ((dust(j)%wav(id)-(lambda_0/1000))/(dust(j)%wav(id)-dust(j)%wav(id-1)))))
           
           lambda_ext_V=lambda_ext_V+dust(j)%weight*(dust(j)%ext_opacity(id_V)-((dust(j)%ext_opacity(id_V)-dust(j)%ext_opacity(id_V-1))* &
                & ((dust(j)%wav(id_V)-(547.0/1000))/(dust(j)%wav(id_V)-dust(j)%wav(id_V-1)))))
           
        END DO
      END SUBROUTINE calculate_opacities

END MODULE
