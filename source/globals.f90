!---------------------------------------------------------------------!
!  declare here all global variables such as counters, constants etc. !
!---------------------------------------------------------------------!

module globals

    implicit none

    !path to input variables
    character(len=12)    ::  input_prefix

    !openmp variables
    integer,external    ::  omp_get_num_threads
    integer,external    ::  omp_get_thread_num
    integer             ::  thread_id
    integer             ::  num_threads

    !counters
    integer             ::  ii,jj,kk
    integer             ::  ixx,iyy,izz
    integer             ::  ish
    integer             ::  i_dir
    integer             ::  i_doublet
    integer             ::  i_packet
    integer             ::  i_clump
    integer             ::  i_line

  !  save i_packet

    !dummy counters
    integer             ::  xx,yy,zz

    !identifiers
    integer             ::  ig
    integer             ::  id_theta,id_phi
    integer             ::  io !error handling

    !random numbers and functions
    real                ::  random(5), ran
    !real,external         ::  r4_uni_01
    !$OMP THREADPRIVATE(random,ran,i_dir)

    !constants
    real, parameter     ::  pi=3.141592654
    real,parameter      ::  c=3e8                   !in si units (m/s)

    !factors for scaling/normalising distributions
    real                ::  norm

    !properties of the model
    integer(8)          ::  n_packets             !number of packets requested
    integer(8)          ::  n_init_packets        !total number of initialised packets
    integer(8)          ::  n_inactive_packets    !total number of inactive packets
    integer(8)          ::  n_abs_packets         !total number of absorbed packets
    integer(8)          ::  n_los_packets         !total number of packets in the line of sight
    integer(8)          ::  no_active_cells=0       !number of active (i.e. non-zero) cells inside ejecta
    integer(8)          ::  n_clumps=0              !actual number of clumps used
    integer(8)          ::  n_recorded_packets    !number of packets contributing to the final line profile
    integer(8)          ::  n_data_points          !number of data points used in chi squared calculation (output in order to calculate red. chi square if desired)
    integer             ::  n_lines                 !number of lines to process if more than one
    integer             ::  n_bins_multiline        !total number of bins to be stored in final array when running multiple lines
    integer             ::  multiline_count         !running counter that tracks number of bins that have been recorded in multline profile array

    integer             ::  n_args                  !number of input arguments
    integer             ::  n_angle_divs            !number of division in each of phi and theta to divide grid into multiple lines of sight
    integer             ::  day_no                  !time in days since outburst - used with v_max to calculate rout
    integer             ::  es_temp                 !electron scattering temperature (0 if no es)
    real                ::  l_halpha                !total halpha luminosity (for e- scattering calcn)
    real                ::  tot_vol                 !total volume of dusty ejecta in 1e42cm^3
    real                ::  tot_vol_gas             !total volume of emitting gas (if different from dust)

    !variables used for checking dust mass calculations
    real                ::  m_tot_check             !calculated total mass of dust using densities and vols
    real                ::  m_icm_check             !calculated mass of dust in inter clump medium
    real                ::  m_clumps_check          !calculated mass of dust in clumps

    !different options and cases
    logical             ::  lg_mcmc                 !use an mcmc routine to map parameter space and return trace and likelihoods
    logical             ::  lg_multiline            !mcmc parameter for running multiple lines - multiple lines used?
    logical             ::  lg_multiline_fixdust    !when running multiple lines and would like to fix the dust distribution for all lines 
    logical             ::  lg_multiline_fixgas     !when running multiple lines and would like to fix the gas distribution for all lines
    logical             ::  lg_los                  !use a line of sight?
    logical             ::  lg_multi_los            !divide grid into multiple lines of sight with complete coverage
    logical             ::  lg_vel_shift            !use velocity shifting to recalculate frequency at every scattering event?
    logical             ::  lg_es                   !use electron scattering?
    logical             ::  lg_data                 !read in an observed line/doublet?
    logical             ::  lg_decoupled            !decouple = 0 or 1 (1 if dust and gas distributions are not coupled)
    logical             ::  lg_doublet              !is the 'line' to be modelled a doublet?
    logical             ::  lg_store_all
    logical             ::  lg_vel_law              !are you using a velocity law independent of radius?

    !properties of resultant line profile
    real                ::  abs_frac                !total fraction of emitted energy that has been absorbed
    real,allocatable    ::  cos_theta_array(:)
    real,allocatable    ::  phi_array(:)
    real                ::  chi_sq                  !measure of goodness of fit of model to data

    !names of input files
    character(len=50)   ::  input_file              !name of input file to read in
    character(len=50)   ::  data_file               !name of data file containing details of observed line to read in
    character(len=50)   ::  data_exclusions_file    !name of data file containing details of regions of observed line to exclude in chi sq. calculation
    character(len=50)   ::  e_scat_file             !name of file containing electron scattering parameters (temperature, halpha luminosity etc.)
    character(len=50)   ::  dust_file               !name of file containing dust grain parameters (geometry, mass etc.)
    character(len=50)   ::  grid_file               !name of file containing dust and gas  grid (x,y,z,dust_rho,gas_rho) for arbitrary geometries
    character(len=50)   ::  gas_file                !name of file containing electron scattering parameters (geometry, luminosity etc.)
    character(len=40)   ::  species_file            !filename containing details of dust species
    character(len=1024) ::  junk                    !dummy for reading in material from files

    !system variables (used in naming output files)
    character(8)        ::  date
    character(10)       ::  time
    character(8)        ::  run_no_string

end module globals
