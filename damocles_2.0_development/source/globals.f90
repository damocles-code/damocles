!---------------------------------------------------------------------!
!  declare here all global variables such as counters, constants etc. !
!---------------------------------------------------------------------!

MODULE globals

    IMPLICIT NONE

    !OpenMP variables
    INTEGER,EXTERNAL    ::  omp_get_num_threads
    INTEGER,EXTERNAL    ::  omp_get_thread_num

    !counters
    INTEGER             ::  ii,jj,kk
    INTEGER             ::  ixx,iyy,izz
    INTEGER             ::  iSh
    INTEGER             ::  i_dir,i_spec

    !dummy counters
    INTEGER             ::  xx,yy,zz

    !identifiers
    INTEGER             ::  iG
    INTEGER             ::  id_no

    !random numbers
    REAL                ::  random(5),ran

    !constants
    REAL, PARAMETER     ::  pi=3.141592654
    REAL,PARAMETER      ::  c=3E8                   !in SI units (m/s)

    !factors for scaling/normalising distributions
    REAL                ::  norm

    !properties of the model
    INTEGER(8)          ::  n_packets               !number of packets requested
    INTEGER(8)          ::  n_init_packets          !total number of initialised packets
    INTEGER(8)          ::  n_inactive_packets      !total number of inactive packets
    INTEGER(8)          ::  n_abs_packets           !total number of absorbed packets
    INTEGER(8)          ::  n_los_packets           !total number of packets in the line of sight
    INTEGER(8)          ::  n_shells                !number of shells to use when generating packets


    INTEGER             ::  nargs                   !number of input arguments
    INTEGER             ::  day_no                  !time in days since outburst - used with v_max to calculate Rout
    INTEGER             ::  ES_temp                 !electron scattering temperature (0 if no ES)
    REAL                ::  L_Halpha                !total Halpha luminosity (for e- scattering calcn)

    !different options and cases
    LOGICAL             ::  lg_LOS                  !use a line of sight?
    LOGICAL             ::  lg_vel_shift            !use velocity shifting to recalculate frequency at every scattering event?
    LOGICAL             ::  lg_ES                   !use electron scattering?
    LOGICAL             ::  lg_data                 !read in an observed line/doublet?
    LOGICAL             ::  lg_decoupled            !decouple = 0 or 1 (1 if dust and gas distributions are not coupled)
    LOGICAL             ::  lg_doublet              !is the 'line' to be modelled a doublet?

    !properties of resultant line profile
    REAL                ::  abs_frac                !total fraction of emitted energy that has been absorbed



    !names of input files
    CHARACTER(LEN=50)   ::  input_file              !name of input file to read in
    CHARACTER(LEN=50)   ::  data_file               !name of data file containing details of observed line to read in
    CHARACTER(LEN=50)   ::  e_scat_file             !name of file containing electron scattering parameters (temperature, Halpha luminosity etc.)
    CHARACTER(LEN=50)   ::  dust_file               !name of file containing dust grain parameters (geometry, mass etc.)
    CHARACTER(LEN=50)   ::  gas_file                !name of file containing electron scattering parameters (geometry, luminosity etc.)
    CHARACTER(LEN=40)   ::  species_file            !filename containing details of dust species



END MODULE globals
