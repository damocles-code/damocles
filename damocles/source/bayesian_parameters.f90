module bayesian_input
  
contains
  
  subroutine read_bayesian_params(params,flags)

    use globals
    use input
    use class_dust
    use initialise
    use vector_functions
    use driver
    
    implicit none
    
    real,intent(in) :: params(21,6)
    real,intent(in) :: flags(21,6)
    
    if (flags(1,i_line) == 1) gas_geometry%v_max = params(1,i_line)*1000
    if (flags(2,i_line) == 1) gas_geometry%r_max = params(2,i_line)
    if (flags(3,i_line) == 1) gas_geometry%r_ratio = params(3,i_line)
    if (flags(4,i_line) == 1) gas_geometry%v_power = params(4,i_line)
    if (flags(5,i_line) == 1) gas_geometry%rho_power = params(5,i_line)
    if (flags(6,i_line) == 1) gas_geometry%emis_power = params(6,i_line)
    if (flags(7,i_line) == 1) gas_geometry%v_min = params(7,i_line)*1000
    if (flags(8,i_line) == 1) gas_geometry%v_prob_indx = params(8,i_line)
    if (flags(9,i_line) == 1) dust%mass = 10**params(9,i_line)
    if (flags(10,i_line) == 1) dust_geometry%clumped_mass_frac = params(10,i_line)
    if (flags(11,i_line) == 1) dust_geometry%ff = params(11,i_line)
    if (flags(12,i_line) == 1) dust_geometry%clump_power = params(12,i_line)
    if (flags(13,i_line) == 1) dust_geometry%v_max = params(13,i_line)*1000
    if (flags(14,i_line) == 1) dust_geometry%r_max = params(14,i_line)
    if (flags(15,i_line) == 1) dust_geometry%v_min = params(15,i_line)*1000
    if (flags(16,i_line) == 1) dust_geometry%v_power = params(16,i_line)
    if (flags(17,i_line) == 1) dust_geometry%rho_power = params(17,i_line)
    if (flags(18,i_line) == 1) dust%species(1)%amin = 10**params(18,i_line)
    if (flags(18,i_line) == 1) dust%species(1)%amax = 10**params(18,i_line)
    if (flags(19,i_line) == 1) line%doublet_ratio = params(19,i_line)
    if (flags(20,i_line) == 1) gas_geometry%ff = params(20,i_line)
    if (flags(21,i_line) == 1) gas_geometry%clump_power = params(21,i_line)

    dust_geometry%r_ratio = dust_geometry%v_min/dust_geometry%v_max
    

  end subroutine read_bayesian_params
  
  
  subroutine check_for_conflicts(flags)
    
    use globals
    use input
    use class_dust
    use initialise
    use vector_functions
    use driver
    
    implicit none
    
    real,intent(in) :: flags(21,6)
    
    if (.not. lg_decoupled) then
       if (flags(1,i_line) == 1) then
          print*, 'You have requested dust and gas coupled. Do not vary gas max velocity or decouple dust and gas. Param #1.'
          STOP
       end if
       if (flags(2,i_line) == 1) then
          print*,'You have requested dust and gas coupled. Do not vary gas max rad us or decouple dust and gas. Param #2.'
          STOP
       end if
       if (flags(3,i_line) == 1) then
          print*, 'You have requested dust and gas coupled. Do not vary gas Rin/Rout ratio or decouple dust and gas. Param #3.'
          STOP
       end if
       if (flags(5,i_line) == 1) then
          print*,'You have requested dust and gas coupled. Do not vary gas density profile or decouple dust and gas. Param #5.'
          STOP
       end if
       if (flags(6,i_line) == 1) then
          print*, 'You have requested dust and gas coupled. Do not vary gas emissivity profile or decouple dust and gas. Param #6.'
          STOP
       end if
    end if
    
    if (.not. lg_vel_law) then
       if (flags(7,i_line) == 1) then
          print*, 'You have not requested a velocity law independent of radius. Do not vary min gas velocity. Param #7.'
          STOP
       end if
       if (flags(8,i_line) == 1) then
          print*, 'You have not requested a velocity law independent of radius. Do not vary velocity probability law. Param #8.'
          STOP
       end if
    else
       if (flags(4,i_line) == 1) then
          print*, 'You have requested a velocity law independent of radius. Do not vary velocity power (v propto r^x). Param #4.'
          STOP
       end if
    end if
    
    if (.not. dust_geometry%lg_clumped) then
       if (flags(11,i_line) == 1) then
          print*, 'You have requested to vary the dust clump filling factor but have not put any dust in clumps. Param #11.'
          STOP
       end if
       if (flags(12,i_line) == 1) then
          print*, 'You have requested to vary the dust clump number distribution but have not put any dust in clumps. Param #12.'
          STOP
       end if
    end if
    
    if (.not. gas_geometry%lg_clumped) then
       if (flags(20,i_line) == 1) then
          print*, 'You have requested to vary the gas clump filling factor but have not requested clumped gas emission. Param #20.'
          STOP
       end if
       if (flags(21,i_line) == 1) then
          print*, 'You have requested to vary the gas clump number &
               & distribution but have not requested clumped gas emission. Param #21.'
          STOP
       end if
    end if
    
    if (lg_multiline_fixdust .and. i_line>1) then
       if ((flags(9,i_line) == 1) .or. &
            & (flags(10,i_line) == 1) .or. &
            & (flags(11,i_line) == 1) .or. &
            & (flags(12,i_line) == 1) .or. &
            & (flags(13,i_line) == 1) .or. &
            & (flags(14,i_line) == 1) .or. &
            & (flags(15,i_line) == 1) .or. &
            & (flags(16,i_line) == 1) .or. &
            & (flags(17,i_line) == 1) .or. &
            & (flags(18,i_line) == 1)) then
          print*,'You have requested to fix the &
               & dust distribution for multiple lines but have &
               & requested to vary a dust parameter in line number'&
               & ,i_line,' Aborting.'
          STOP
       end if
    end if
    
    if (lg_multiline_fixgas .and. i_line>1) then
       if ((flags(1,i_line) == 1) .or. &
            & (flags(2,i_line) == 1) .or. &
            & (flags(3,i_line) == 1) .or. &
            & (flags(4,i_line) == 1) .or. &
            & (flags(5,i_line) == 1) .or. &
            & (flags(6,i_line) == 1) .or. &
            & (flags(7,i_line) == 1) .or. &
            & (flags(8,i_line) == 1) .or. &
            & (flags(20,i_line) == 1) .or. &
            & (flags(21,i_line) == 1)) then
          print*,'You have requested to fix the &
               & dust distribution for multiple lines but have &
               & requested to vary a dust parameter in line number' &
               & ,i_line,' Aborting.'
          STOP
       end if
    end if
    
  end subroutine check_for_conflicts
end module bayesian_input
