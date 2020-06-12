!-----------------------------------------------------------------------------------!
!  this is the main driver of the code damocles.  it is included as a module        !
!  in order to allow it to be run from other programs e.g. python wrappers.         !
!  the run_damocles subroutine calls the subroutines that construct the grids,      !
!                                                                                   !
!  emit and propagate packets through the grid and collates all escaped packets.    !
!  the model comparison module is also called from here.                            !
!-----------------------------------------------------------------------------------!

module driver

  use globals
  use class_line
  use class_freq_grid
  use class_grid
  use electron_scattering
  use input
  use initialise
  use vector_functions
  use class_packet
  use radiative_transfer
  use model_comparison
  
  implicit none
  
contains
  
  subroutine run_damocles()
    integer :: thread_id
    
    !RUN DAMOCLES
    !read input
    if (.not. lg_mcmc) call read_input()
    
    !construct grids and initialise simulation:
    do i_doublet=1,2
       
       !construct all grids and initialise rest line wavelength/freq
       if (i_doublet==1) then
          !set active rest frame wavelength
          line%wavelength=line%doublet_wavelength_1
          line%frequency=c*10**9/line%wavelength

          !construct grids
          if (lg_mcmc) then
             if (((lg_multiline_fixdust) .and. (i_line == 1)) &
               & .or. (.not. lg_multiline_fixdust)) then
                call calculate_opacities()
                call build_dust_grid()
             end if
          else
             call calculate_opacities()
             call build_dust_grid()
          end if
          
          call init_random_seed()
          call construct_freq_grid()
          call build_emissivity_dist()
          if (lg_es) call n_e_const()

          !                print*,'WARNING:  USING RANDOM VELOCITIES NOT POWER-LAW VELOCITIES.  SEE CLASS_PACKET FILE TO AMEND.'
          
          !build multiple lines of sight array
          allocate(cos_theta_array(n_angle_divs))
          allocate(phi_array(n_angle_divs))
          do ii=1,n_angle_divs
             cos_theta_array(ii) = (2*real(ii-1)/real(n_angle_divs))-1
             phi_array(ii)=2*real(ii-1)*pi/(real(n_angle_divs))
          end do
          
          !initialise counters to zero
          n_init_packets=0
          n_inactive_packets=0
          n_abs_packets=0
          n_los_packets=0
          n_recorded_packets=0
          
       else if (i_doublet==2) then
          !exit if not a doublet
          if (.not. lg_doublet) exit
          
          !otherwise reset rest frame wavelength of line to be modelled
          !active wavelength is now second component of doublet
          line%wavelength=line%doublet_wavelength_2
          line%frequency=c*10**9/line%wavelength
       end if
       
       !emit and propagate packets through grid
       if (.not. lg_mcmc) print*,"propagating packets..."

       !prepare parallelised region with number of threads to use
       call omp_set_num_threads(num_threads)
       !$OMP PARALLEL DEFAULT(FIRSTPRIVATE) PRIVATE(id_theta,id_phi,ixx,iyy,izz,thread_id)  REDUCTION(+:n_abs_packets,abs_frac,profile_array,profile_los_array,n_los_packets,n_inactive_packets,n_init_packets,n_recorded_packets)
       
       !$OMP DO 
       !SCHEDULE(DYNAMIC)
       do ii=1,n_packets
          packet%id = ii
          call run_packet()
       end do
       !$OMP END DO
       !$OMP END PARALLEL

    end do
    
    !calculate energies
    line%initial_energy=line%luminosity/real(n_init_packets)
    
    !calculate goodness of fit to data if supplied
    if (lg_data) then
       call calculate_chi_sq()
    end if
    
    !write out log file
    if (.not. lg_mcmc) call write_to_file()
    
    !decallocate all allocated memory
    if (.not. lg_multiline_fixdust) then
       deallocate(grid_cell)
       deallocate(mothergrid%x_div)
       deallocate(mothergrid%y_div)
       deallocate(mothergrid%z_div)
       deallocate(dust%species)
    end if
    deallocate(nu_grid%lambda_bin)
    deallocate(nu_grid%vel_bin)
    deallocate(nu_grid%bin)
    deallocate(cos_theta_array)
    deallocate(phi_array)
    if (lg_data) then
       deallocate(obs_data%vel)
       deallocate(obs_data%flux)
       deallocate(obs_data%exclude)
       deallocate(obs_data%freq)
       deallocate(obs_data%error)
       deallocate(exclusion_zone)
       if (.not. lg_mcmc) deallocate(profile_array_data_bins)
       if (.not. lg_mcmc) deallocate(mc_error_data_bins)
       deallocate(square_weight_data_bins)
       deallocate(total_weight_data_bins)
       deallocate(n_packets_data_bins)
    end if
    deallocate(profile_array)
    deallocate(profile_los_array)
    
    if (.not. lg_mcmc) print*,'complete!'
    
    !handle the multiline case
!    if (lg_multiline) then
!       multiline_profile_array(1+multiline_count:multiline_count+obs_data%n_data,1) = profile_array_data_bins(:)
!       multiline_profile_array(1+multiline_count:multiline_count+obs_data%n_data,2) = mc_error_data_bins(:)
!       deallocate(profile_array_data_bins)
!       deallocate(mc_error_data_bins)
!       multiline_count= multiline_count + obs_data%n_data
!    end if
    
end subroutine run_damocles

!!this approach is not necessarily well parallelised in the case of high optical depths
!!this should be checked
recursive subroutine run_packet()
  
  
  packet%lg_abs = .false.
  packet%lg_active = .false.
  
  call emit_packet()
  !$OMP CRITICAL
  n_init_packets = n_init_packets+1
  !$OMP END CRITICAL
  
  if (packet%lg_active) then
     !prevent emission coming from clumps
     !     packet%cell_no=(mothergrid%n_cells(2)*mothergrid%n_cells(3)*(packet%axis_no(1)-1)+mothergrid%n_cells(3)*(packet%axis_no(2)-1)+packet%axis_no(3))
     !     if (grid_cell(packet%cell_no)%lg_clump) then
     !        call run_packet()
     !     end if

     !propagate active packet through grid
     call propagate()
     !if packet has been absorbed then record
     if (packet%lg_active) then
        if (packet%lg_abs) then
           if (i_doublet==2) then
              abs_frac=abs_frac+packet%weight/line%doublet_ratio
           else
              abs_frac=abs_frac+packet%weight
           end if
           !$OMP CRITICAL
           n_abs_packets=n_abs_packets+1
           !$OMP END CRITICAL
           call run_packet()
        else
           !if the packet has not been absorbed then record in resultant profile
           call add_packet_to_profile()                
           
        end if  !absorbed/escaped
     else
        !if emitted packet not active then run another one
        n_inactive_packets = n_inactive_packets+1
        call run_packet()
     end if !propagated and active
  else
     n_inactive_packets = n_inactive_packets+1          
     call run_packet()
  end if  !emitted and active
  
end subroutine run_packet

subroutine add_packet_to_profile()
  
  !find the smallest distance and thus nearest freq point
  packet%freq_id=minloc(packet%nu-nu_grid%bin(:,1),1,(packet%nu-nu_grid%bin(:,1))>0)

  if (packet%freq_id==0) then
     !            print*,'photon outside frequency range',packet%freq_id,packet%nu,packet%weight
     !$OMP CRITICAL
     n_inactive_packets = n_inactive_packets+1
     !$OMP END CRITICAL
     call run_packet()
  else
     !adjust weight of packet if second component of doublet
     if (i_doublet==2) then
        packet%weight=packet%weight/line%doublet_ratio
     end if
     
     !add packet to primary profile array
     if (.not. lg_los) then
        profile_array(packet%freq_id)=profile_array(packet%freq_id)+packet%weight
        !increment number of recorded packets
        !$OMP CRITICAL
        n_recorded_packets = n_recorded_packets+1
        !$OMP END CRITICAL
     else
        if (acos(packet%dir_sph(1)) < pi/12) then
           profile_array(packet%freq_id)=profile_array(packet%freq_id)+packet%weight
           !increment number of recorded packets
           !$OMP CRITICAL
           n_recorded_packets = n_recorded_packets+1
           !$OMP END CRITICAL
        end if
     end if
     
     !add packet to line of sight
     if (lg_multi_los) then
        id_theta = minloc(packet%dir_sph(1)-cos_theta_array,1,(packet%dir_sph(1)-cos_theta_array)>0)
        id_phi = minloc(packet%dir_sph(2)-phi_array,1,(packet%dir_sph(2)-phi_array)>0)
        profile_los_array(packet%freq_id,id_theta,id_phi)=profile_los_array(packet%freq_id,id_theta,id_phi)+packet%weight
     end if
     
     
     !incremement number of packets in line of sight
     n_los_packets=n_los_packets+1
     
     !calculate the sum of the squares of the packet weights and the sum of the weights
     if (lg_data) then
        packet%freq_id=minloc(packet%nu-obs_data%freq(:),1,(packet%nu-obs_data%freq(:))>0)
        if (packet%freq_id >1) then
           !$OMP CRITICAL
           n_packets_data_bins(packet%freq_id) = n_packets_data_bins(packet%freq_id) + 1
           total_weight_data_bins(packet%freq_id) = total_weight_data_bins(packet%freq_id) + packet%weight
           square_weight_data_bins(packet%freq_id) = square_weight_data_bins(packet%freq_id) + packet%weight**2
           !$OMP END CRITICAL
        else
           if (.not. lg_mcmc) print*, 'WARNING: packet out of data wavelength range - either increase the range of the data or reduce v_max.'
           call run_packet()
        end if
     end if
  end if
  
end subroutine add_packet_to_profile

end module driver
     
