!********************************************************************************
!> \brief Input/Output module
!
!> This module contains all the input/output subroutine and the 
!> realted variables.
!
!> \date 07/10/2016
!> @author 
!> Mattia de' Michieli Vitturi
!
!********************************************************************************

MODULE inpout_2d

  USE parameters_2d, ONLY : wp

  ! -- Variables for the namelist RUN_PARAMETERS
  USE parameters_2d, ONLY : t_start , t_end , t_output , dt_output 

  USE solver_2d, ONLY : verbose_level

  ! -- Variables for the namelist NEWRUN_PARAMETERS
  USE geometry_2d, ONLY : x0 , y0 , comp_cells_x , comp_cells_y , cell_size
  USE geometry_2d, ONLY : topography_profile , n_topography_profile_x ,         &
       n_topography_profile_y , nodata_topo
  USE parameters_2d, ONLY : n_solid , n_add_gas
  USE parameters_2d, ONLY : rheology_flag , energy_flag , alpha_flag ,          &
       topo_change_flag , radial_source_flag , collapsing_volume_flag ,         &
       liquid_flag , gas_flag , subtract_init_flag , bottom_radial_source_flag ,&
       vertical_profiles_flag , lateral_source_flag

  USE parameters_2d, ONLY : slope_correction_flag , curvature_term_flag 

  ! -- Variables for the namelist INITIAL_CONDITIONS
  USE parameters_2d, ONLY : released_volume , x_release , y_release
  USE parameters_2d, ONLY : velocity_mod_release , velocity_ang_release
  USE parameters_2d, ONLY : alphas_init
  USE parameters_2d, ONLY : T_init
  USE parameters_2d, ONLY : u_init
  USE parameters_2d, ONLY : v_init

  ! -- Variables for the namelists LEFT/RIGHT_BOUNDARY_CONDITIONSq
  USE parameters_2d, ONLY : bc

  ! -- Variables for the namelist NUMERIC_PARAMETERS
  USE parameters_2d, ONLY : solver_scheme, dt0 , max_dt , cfl, limiter , theta, &
       reconstr_coeff , interfaces_relaxation , n_RK   

  ! -- Variables for the namelist EXPL_TERMS_PARAMETERS
  USE constitutive_2d, ONLY : grav , inv_grav

  ! -- Variables for the namelist RADIAL_SOURCE_PARAMETERS
  USE parameters_2d, ONLY : x_source , y_source , r_source , vel_source ,       &
       T_source , h_source , alphas_source , alphal_source , alphag_source ,    &
       time_param , Ri_source , mfr_source , xs_source , xl_source , xg_source

  ! -- Additional variables for the namelist LATERAL_SOURCE_PARAMETERS
  USE parameters_2d, ONLY : source_side , x1_source , x2_source , y1_source ,   &
       y2_source
  
  ! -- Variables for the namelist COLLAPSING_VOLUME_PARAMETERS
  USE parameters_2d, ONLY : x_collapse , y_collapse , r_collapse , T_collapse , &
       h_collapse , alphas_collapse , alphag_collapse

  ! -- Variables for the namelist TEMPERATURE_PARAMETERS
  USE constitutive_2d, ONLY : emissivity , exp_area_fract , enne , emme ,       &
       atm_heat_transf_coeff , thermal_conductivity , T_env , T_ground , c_p

  ! -- Variables for the namelist RHEOLOGY_PARAMETERS
  USE parameters_2d, ONLY : rheology_model
  USE constitutive_2d, ONLY : mu , xi , tau , nu_ref , visc_par , T_ref,        &
                               mu_0, mu_inf, Fr_0
  USE constitutive_2d, ONLY : alpha2 , beta2 , alpha1_coeff , beta1 , Kappa ,n_td
  USE constitutive_2d, ONLY : friction_factor
  USE constitutive_2d, ONLY : tau0

  ! --- Variables for the namelist SOLID_TRANSPORT_PARAMETERS
  USE constitutive_2d, ONLY : rho_s , diam_s , sp_heat_s
  USE constitutive_2d, ONLY : settling_flag , erosion_coeff , erodible_porosity
  USE constitutive_2d, ONLY : erodible_fract , T_erodible
  USE constitutive_2d, ONLY : alphastot_min
  USE parameters_2d, ONLY : erodible_deposit_flag

  ! --- Variables for the namelist GAS_TRANSPORT_PARAMETERS
  USE constitutive_2d, ONLY : sp_heat_a , sp_gas_const_a , kin_visc_a , pres ,  &
       T_ambient , entrainment_flag , sp_heat_g , sp_gas_const_g , gamma_steam

  USE parameters_2d, ONLY : liquid_vaporization_flag , water_level

  ! --- Variables for the namelist LIQUID_TRANSPORT_PARAMETERS
  USE constitutive_2d, ONLY : sp_heat_l , rho_l , kin_visc_l , loss_rate

  ! --- Variables for the namelist VULNERABILITY_TABLE_PARAMETERS
  USE parameters_2d, ONLY : n_thickness_levels , n_dyn_pres_levels ,            &
       thickness_levels , dyn_pres_levels

  ! --- Variables for the namelist VERTICAL_PROFILES_PARAMETERS
  USE constitutive_2d, ONLY : vonK , k_s , Sc
  USE parameters_2d, ONLY : bottom_conc_flag , n_quad

  IMPLICIT NONE

  INTEGER :: n_restart_files
  CHARACTER(LEN = 40), dimension(10) :: restart_files

  CHARACTER(LEN=40) :: run_name           !< Name of the run
  CHARACTER(LEN=40) :: bak_name           !< Backup file for the parameters
  CHARACTER(LEN=40) :: input_file         !< File with the run parameters
  CHARACTER(LEN=40) :: output_file        !< Name of the output files
  CHARACTER(LEN=40) :: restart_file       !< Name of the restart file 
  CHARACTER(LEN=40) :: probes_file        !< Name of the probes file 
  CHARACTER(LEN=40) :: output_file_2d     !< Name of the output files
  CHARACTER(LEN=40) :: output_esri_file   !< Name of the esri output files
  CHARACTER(LEN=40) :: output_max_file    !< Name of the esri max. thick. file
  CHARACTER(LEN=40) :: runout_file        !< Name of the runout file 
  CHARACTER(LEN=40) :: topography_file    !< Name of the esri DEM file
  CHARACTER(LEN=40) :: erodible_file      !< Name of the esri DEM file
  CHARACTER(LEN=40) :: output_VT_file
  CHARACTER(LEN=40) :: mass_center_file

  INTEGER, PARAMETER :: input_unit = 7       !< Input data unit
  INTEGER, PARAMETER :: backup_unit = 8      !< Backup input data unit
  INTEGER, PARAMETER :: output_unit = 9      !< Output data unit
  INTEGER, PARAMETER :: restart_unit = 10    !< Restart data unit
  INTEGER, PARAMETER :: probes_unit = 11     !< Probes data unit
  INTEGER, PARAMETER :: output_unit_2d = 12  !< Output data 2D unit
  INTEGER, PARAMETER :: output_esri_unit = 13  !< Esri Output unit
  INTEGER, PARAMETER :: output_max_unit = 14  !< Esri max thick. output unit
  INTEGER, PARAMETER :: dem_esri_unit = 15   !< Computational grid Esri fmt unit
  INTEGER, PARAMETER :: runout_unit = 16
  INTEGER, PARAMETER :: dakota_unit = 17
  INTEGER, PARAMETER :: erodible_unit = 18
  INTEGER, PARAMETER :: output_VT_unit = 19
  INTEGER, PARAMETER :: mass_center_unit = 20

  !> Counter for the output files
  INTEGER :: output_idx 

  !> Flag to start a run from a previous output:\n
  !> - T     => Restart from a previous output (.asc or .q_2d)
  !> - F     => Restart from initial condition read from two_phases.inp
  !> .
  LOGICAL :: restart

  !> Flag to save the output in esri ascii format *.asc
  !> - T     => write esri file
  !> - F     => do not write esri file
  !> .
  LOGICAL :: output_esri_flag

  !> Flag to save the physical variables on file *.p_2d
  !> - T     => write physical variables on file
  !> - F     => do not write the physical variables
  !> .
  LOGICAL :: output_phys_flag

  !> Flag to save the conservative variables on file *.q_2d
  !> - T     => write conservative variables on file
  !> - F     => do not write the conservative variables
  !> .
  LOGICAL :: output_cons_flag

  !> Flag to save the max runout at ouput times
  !> - T     => write max runout on file
  !> - F     => do not write max runout
  !> .
  LOGICAL :: output_runout_flag

  ! -- Variables for the namelists WEST_BOUNDARY_CONDITIONS
  TYPE(bc) :: h_bcW , hu_bcW , hv_bcW , T_bcW
  TYPE(bc),ALLOCATABLE :: alphas_bcW(:)
  TYPE(bc),ALLOCATABLE :: halphas_bcW(:)
  TYPE(bc),ALLOCATABLE :: alphag_bcW(:)
  TYPE(bc),ALLOCATABLE :: halphag_bcW(:)
  TYPE(bc) :: alphal_bcW , halphal_bcW

  ! -- Variables for the namelists EAST_BOUNDARY_CONDITIONS
  TYPE(bc) :: h_bcE , hu_bcE , hv_bcE , T_bcE
  TYPE(bc),ALLOCATABLE :: alphas_bcE(:)
  TYPE(bc),ALLOCATABLE :: halphas_bcE(:)
  TYPE(bc),ALLOCATABLE :: alphag_bcE(:)
  TYPE(bc),ALLOCATABLE :: halphag_bcE(:)
  TYPE(bc) :: alphal_bcE , halphal_bcE

  ! -- Variables for the namelists SOUTH_BOUNDARY_CONDITIONS
  TYPE(bc) :: h_bcS , hu_bcS , hv_bcS , T_bcS
  TYPE(bc),ALLOCATABLE :: alphas_bcS(:)
  TYPE(bc),ALLOCATABLE :: halphas_bcS(:)
  TYPE(bc),ALLOCATABLE :: alphag_bcS(:)
  TYPE(bc),ALLOCATABLE :: halphag_bcS(:)
  TYPE(bc) :: alphal_bcS , halphal_bcS

  ! -- Variables for the namelists NORTH_BOUNDARY_CONDITIONS
  TYPE(bc) :: h_bcN , hu_bcN , hv_bcN , T_bcN
  TYPE(bc),ALLOCATABLE :: alphas_bcN(:)
  TYPE(bc),ALLOCATABLE :: halphas_bcN(:)
  TYPE(bc),ALLOCATABLE :: alphag_bcN(:)
  TYPE(bc),ALLOCATABLE :: halphag_bcN(:)
  TYPE(bc) :: alphal_bcN , halphal_bcN



  ! parameters to read a dem file
  INTEGER :: ncols, nrows, nodata_value

  REAL(wp) :: xllcorner, yllcorner, cellsize

  LOGICAL :: write_first_q

  INTEGER :: n_probes

  REAL(wp), ALLOCATABLE :: probes_coords(:,:)

  REAL(wp) :: dt_runout

  REAL(wp) :: dt_probes

  REAL(wp) :: x_mass_center_old , y_mass_center_old

  REAL(wp) :: x0_runout, y0_runout , init_runout , init_runout_x ,              &
       init_runout_y , eps_stop

  ! absolute precentages of solids in the initial volume
  REAL(wp), ALLOCATABLE :: sed_vol_perc(:)

  REAL(wp) :: initial_erodible_thickness

  REAL(wp) :: alphas0_E(10) , alphas0_W(10)

  REAL(wp) :: alpha1_ref

  REAL(wp) :: thickness_levels0(10)
  REAL(wp) :: dyn_pres_levels0(10)

  REAL(wp) :: release_time(10)

  NAMELIST / run_parameters / run_name , restart , t_start , t_end , dt_output ,&
       output_cons_flag , output_esri_flag , output_phys_flag ,                 &
       output_runout_flag , verbose_level
  
  NAMELIST / restart_parameters / n_restart_files, restart_files, release_time ,&
       T_init , T_ambient , u_init , v_init , sed_vol_perc

  NAMELIST / newrun_parameters / n_solid , topography_file , x0 , y0 ,          &
       comp_cells_x , comp_cells_y , cell_size , rheology_flag , alpha_flag ,   &
       energy_flag , liquid_flag , radial_source_flag , collapsing_volume_flag ,&
       topo_change_flag , gas_flag , subtract_init_flag , n_add_gas ,           &
       bottom_radial_source_flag , slope_correction_flag , curvature_term_flag ,&
       vertical_profiles_flag , lateral_source_flag

  NAMELIST / initial_conditions /  released_volume , x_release , y_release ,    &
       velocity_mod_release , velocity_ang_release , T_init , T_ambient

  NAMELIST / numeric_parameters / solver_scheme, dt0 , max_dt , cfl, limiter ,  &
       theta , reconstr_coeff , interfaces_relaxation , n_RK   

  NAMELIST / expl_terms_parameters / grav

  NAMELIST / radial_source_parameters / x_source , y_source , r_source ,        &
       vel_source , T_source , h_source , alphas_source , alphal_source ,       &
       alphag_source , time_param , Ri_source , mfr_source , xs_source ,        &
       xl_source , xg_source

  NAMELIST / lateral_source_parameters / source_side , x1_source , x2_source ,  &
       y1_source , y2_source , vel_source , T_source , h_source , alphas_source,&
       alphal_source , alphag_source , time_param , Ri_source , mfr_source ,    &
       xs_source , xl_source , xg_source

  NAMELIST / collapsing_volume_parameters / x_collapse , y_collapse ,           &
       r_collapse , T_collapse , h_collapse , alphas_collapse , alphag_collapse

  NAMELIST / temperature_parameters / emissivity ,  atm_heat_transf_coeff ,     &
       thermal_conductivity , exp_area_fract , c_p , enne , emme , T_env ,      &
       T_ground

  NAMELIST / rheology_parameters / rheology_model , mu , xi , tau , nu_ref ,    &
       visc_par , T_ref , alpha2 , beta2 , alpha1_ref , beta1 , Kappa , n_td ,  &
       friction_factor , tau0, mu_0, mu_inf, Fr_0

  NAMELIST / runout_parameters / x0_runout , y0_runout , dt_runout ,            &
       eps_stop

  NAMELIST / liquid_transport_parameters / sp_heat_l , rho_l , kin_visc_l ,     &
       loss_rate

  NAMELIST / vulnerability_table_parameters / thickness_levels0 , dyn_pres_levels0

  NAMELIST / vertical_profiles_parameters / vonK , k_s  , Sc , bottom_conc_flag,&
       n_quad

CONTAINS

  !******************************************************************************
  !> \brief Initialization of the variables read from the input file
  !
  !> This subroutine initialize the input variables with default values
  !> that solve for a Riemann problem. If the input file does not exist
  !> one is created with the default values.
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 07/10/2016
  !
  !******************************************************************************

  SUBROUTINE init_param

    USE parameters_2d , ONLY : n_vars

    IMPLICIT none

    LOGICAL :: lexist
    INTEGER :: ios

    n_vars = 3

    !-- Inizialization of the Variables for the namelist RUN_PARAMETERS
    run_name = 'default'
    restart = .FALSE.
    t_start = 0.0
    t_end = 5.0E-2_wp
    dt_output = 5.0E-3_wp
    output_cons_flag = .TRUE.
    output_esri_flag = .TRUE.
    output_phys_flag = .TRUE.
    output_runout_flag = .FALSE.
    verbose_level = 0

    !-- Inizialization of the Variables for the namelist RESTART_PARAMETERS
    n_restart_files = 1
    restart_files(:) = ''
    restart_file = ''
    release_time(:) = 0.0_wp
    T_init = 0.0_wp
    T_ambient = 0.0_wp
    u_init = 0.0_wp
    v_init = 0.0_wp

    !-- Inizialization of the Variables for the namelist newrun_parameters
    topography_file = 'topography_dem.asc'
    nodata_topo = -9999.0_wp
    x0 = 0.0_wp
    y0 = 0.0_wp
    comp_cells_x = 1000
    comp_cells_y = 1
    cell_size = 1.0E-3_wp
    n_solid = -1
    n_add_gas = -1
    rheology_flag = .FALSE.
    energy_flag = .FALSE.
    topo_change_flag = .FALSE.
    radial_source_flag = .FALSE.
    lateral_source_flag = .FALSE.
    collapsing_volume_flag = .FALSE.
    bottom_radial_source_flag = .FALSE.
    liquid_flag = .FALSE.
    gas_flag = .TRUE.
    subtract_init_flag = .FALSE.
    alpha_flag = .FALSE.
    slope_correction_flag = .FALSE.
    curvature_term_flag  = .FALSE.
    vertical_profiles_flag = .FALSE.

    !-- Inizialization of the Variables for the namelist NUMERIC_PARAMETERS
    dt0 = 1.0E-4_wp
    max_dt = 1.0E-3_wp
    solver_scheme = 'KT'
    n_RK = 2
    cfl = 0.24_wp
    limiter(1:n_vars+2) = 1
    theta=1.0
    reconstr_coeff = 1.0

    !-- Inizialization of the Variables for the namelist EXPL_TERMS_PARAMETERS
    grav = 9.81_wp
    inv_grav = 1.0_wp / grav

    !-- Inizialization of the Variables for the namelist TEMPERATURE_PARAMETERS
    exp_area_fract = 0.5_wp
    emissivity = 0.0_wp                 ! no radiation to atmosphere
    atm_heat_transf_coeff = 0.0_wp      ! no convection to atmosphere
    thermal_conductivity = 0.0_wp       ! no conduction to ground
    enne = 4.0_wp
    emme = 12.0_wp
    T_env = 300.0_wp
    T_ground = 1200.0_wp
    c_p = 1200.0_wp

    !-- Inizialization of the Variables for the namelist RHEOLOGY_PARAMETERS
    rheology_model = 0
    nu_ref = 0.0_wp                     
    mu = -1.0_wp
    xi = -1.0_wp
    tau = 0.0_wp
    T_ref = 0.0_wp
    visc_par = 0.0_wp
    tau0 = 0.0_wp
    mu_0 = -1.0_wp
    mu_inf = -1.0_wp
    Fr_0 = -1.0_wp


    !-- Inizialization of the Variables for the namelist RUNOUT_PARAMETERS
    x0_runout = -1
    y0_runout = -1
    dt_runout = 60
    eps_stop = 0.0_wp

    !-------------- Check if input file exists ----------------------------------
    input_file = 'IMEX_SfloW2D.inp'

    INQUIRE (FILE=input_file,exist=lexist)

    IF (lexist .EQV. .FALSE.) THEN

       WRITE(*,*) 'Input file IMEX_SfloW2D.inp not found'
       STOP

    ELSE

       OPEN(input_unit,FILE=input_file,STATUS='old')

       ! ------- READ run_parameters NAMELIST -----------------------------------
       READ(input_unit, newrun_parameters,IOSTAT=ios )

       IF ( ios .NE. 0 ) THEN

          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist NEWRUN_PARAMETERS'
          WRITE(*,newrun_parameters)
          WRITE(*,*) 'Please check the input file'
          STOP

       ELSE

          REWIND(input_unit)

       END IF

       IF ( n_solid .LT. 0 ) THEN

          WRITE(*,*) 'ERROR: problem with namelist NEWRUN_PARAMETERS'
          WRITE(*,*) 'n_solid =' , n_solid
          WRITE(*,*) 'Please check the input file'
          STOP

       END IF

       ALLOCATE ( alphas_bcW(n_solid) )
       ALLOCATE ( alphas_bcE(n_solid) )
       ALLOCATE ( alphas_bcS(n_solid) )
       ALLOCATE ( alphas_bcN(n_solid) )

       ALLOCATE ( halphas_bcW(n_solid) )
       ALLOCATE ( halphas_bcE(n_solid) )
       ALLOCATE ( halphas_bcS(n_solid) )
       ALLOCATE ( halphas_bcN(n_solid) )

       ALLOCATE( sed_vol_perc(n_solid) )
       sed_vol_perc(1:n_solid) = -1.0_wp

       ALLOCATE( rho_s(n_solid) )
       ALLOCATE( diam_s(n_solid) )
       ALLOCATE( sp_heat_s(n_solid) )
       ALLOCATE( erodible_fract(n_solid) )

       IF ( n_add_gas .LT. 0 ) THEN

          WRITE(*,*) 'ERROR: problem with namelist NEWRUN_PARAMETERS'
          WRITE(*,*) 'n_add_gas =' , n_add_gas
          WRITE(*,*) 'Please check the input file'
          STOP

       END IF

       ALLOCATE ( alphag_bcW(n_add_gas) )
       ALLOCATE ( alphag_bcE(n_add_gas) )
       ALLOCATE ( alphag_bcS(n_add_gas) )
       ALLOCATE ( alphag_bcN(n_add_gas) )

       ALLOCATE ( halphag_bcW(n_add_gas) )
       ALLOCATE ( halphag_bcE(n_add_gas) )
       ALLOCATE ( halphag_bcS(n_add_gas) )
       ALLOCATE ( halphag_bcN(n_add_gas) )

       ALLOCATE ( sp_heat_g(n_add_gas) )
       ALLOCATE ( sp_gas_const_g(n_add_gas) )

       CLOSE(input_unit)

    END IF

    ! output file index
    output_idx = 0

    ! -------------- Initialize values for checks during input reading ----------
    h_bcW%flag = -1 
    hu_bcW%flag = -1 
    hv_bcW%flag = -1 
    alphas_bcW%flag = -1 
    halphas_bcW%flag = -1 
    alphag_bcW%flag = -1 
    halphag_bcW%flag = -1 
    alphal_bcW%flag = -1 
    halphal_bcW%flag = -1 
    T_bcW%flag = -1


    h_bcE%flag = -1 
    hu_bcE%flag = -1 
    hv_bcE%flag = -1 
    alphas_bcE%flag = -1 
    halphas_bcE%flag = -1 
    alphag_bcE%flag = -1 
    halphag_bcE%flag = -1 
    alphal_bcE%flag = -1 
    halphal_bcE%flag = -1 
    T_bcE%flag = -1 

    h_bcS%flag = -1 
    hu_bcS%flag = -1 
    hv_bcS%flag = -1 
    alphas_bcS%flag = -1 
    halphas_bcS%flag = -1 
    alphag_bcS%flag = -1 
    halphag_bcS%flag = -1 
    alphal_bcS%flag = -1 
    halphal_bcS%flag = -1 
    T_bcS%flag = -1 

    h_bcN%flag = -1 
    hu_bcN%flag = -1 
    hv_bcN%flag = -1 
    alphas_bcN%flag = -1 
    halphas_bcN%flag = -1 
    alphag_bcN%flag = -1 
    halphag_bcN%flag = -1 
    alphal_bcN%flag = -1 
    halphal_bcN%flag = -1 
    T_bcN%flag = -1 

    ! sed_vol_perc = -1.0_wp

    rheology_model = -1
    mu = -1
    xi = -1
    tau = -1
    nu_ref = -1
    visc_par = -1
    T_ref = -1
    tau0 = -1
    friction_factor = -1
    mu_0 = -1
    mu_inf = -1
    Fr_0 = -1

    alpha2 = -1.0_wp
    beta2 = -1.0_wp
    alpha1_ref = -1.0_wp
    beta1 = -1.0_wp
    Kappa = -1.0_wp
    n_td = -1.0_wp
    rho_s = -1.0_wp
    initial_erodible_thickness = -1.0_wp
    erodible_deposit_flag = .FALSE.

    exp_area_fract = -1.0_wp
    emissivity = -1.0_wp             
    atm_heat_transf_coeff = -1.0_wp
    thermal_conductivity = -1.0_wp  
    enne = -1.0_wp
    emme = -1.0_wp
    T_env = -1.0_wp
    T_ground = -1.0_wp
    c_p = -1.0_wp

    grav = -1.0_wp

    x0_runout = -1.0_wp
    y0_runout = -1.0_wp

    !- Variables for the namelist SOLID_TRANSPORT_PARAMETERS
    rho_s = -1.0_wp
    diam_s = -1.0_wp
    sp_heat_s = -1.0_wp
    settling_flag = .FALSE.
    erosion_coeff = -1.0_wp
    n_solid = -1
    T_erodible = -1.0_wp
    erodible_file = ""
    erodible_fract = -1.0_wp
    erodible_porosity = -1.0_wp
    alphastot_min = 0.0_wp

    !- Variables for the namelist RHEOLOGY_PARAMETERS
    xi = -1.0_wp
    mu = -1.0_wp
    mu_0 = -1.0_wp
    mu_inf = -1.0_wp
    Fr_0 = -1.0_wp

    !- Variables for the namelist GAS_TRANSPORT_PARAMETERS
    sp_heat_a = -1.0_wp
    sp_gas_const_a = -1.0_wp
    kin_visc_a = -1.0_wp
    pres = -1.0_wp
    T_ambient = -1.0_wp
    sp_heat_g = -1.0_wp
    sp_gas_const_g = -1.0_wp
    liquid_vaporization_flag = .FALSE.
    water_level = -1.0E7_wp
    gamma_steam = -1.0_wp

    !- Variables for the namelist LIQUID_TRANSPORT_PARAMETERS
    sp_heat_l = -1.0_wp
    rho_l = -1.0_wp
    kin_visc_l = -1.0_wp
    loss_rate = -1.0_wp

    !- Variables for the namelist RADIAL_SOURCE_PARAMETERS
    T_source = -1.0_wp
    h_source = -1.0_wp
    r_source = -1.0_wp
    vel_source = -1.0_wp
    Ri_source = -1.0_wp
    mfr_source = -1.0_wp
    time_param(1:4) = -1.0_wp
    alphas_source = -1.0_wp
    alphag_source = -1.0_wp
    alphal_source = -1.0_wp
    xs_source = -1.0_wp
    xg_source = -1.0_wp
    xl_source = -1.0_wp

    !- Variable for the namelist LATERAL_SOURCE_PARAMETERS
    source_side = ''
    x1_source = -1.0_wp
    x2_source = -1.0_wp
    y1_source = -1.0_wp
    y2_source = -1.0_wp    
    
    !- Variables for the namelist COLLAPSING_VOLUME_PARAMETERS
    T_collapse = -1.0_wp
    h_collapse = -1.0_wp
    r_collapse = -1.0_wp
    alphag_collapse = -1.0_wp
    alphas_collapse = -1.0_wp

    !- Variables for the namelist VULNERABILTY_TABLE_PARAMETERS
    n_thickness_levels = -1
    n_dyn_pres_levels = -1
    thickness_levels0 = -1.0_wp
    dyn_pres_levels0 = -1.0_wp

    !- Initial values for VERTICAL_PROFILES_PARAMETERS
    vonK = -1
    k_s = -1
    Sc = -1
    bottom_conc_flag = .FALSE.
    n_quad = -1

  END SUBROUTINE init_param

  !******************************************************************************
  !> \brief Read the input file
  !
  !> This subroutine read the input parameters from the file 
  !> "two_phases.inp" and write a backup file of the input parameters 
  !> with name "run_name.bak", where run_name is read from the input
  !> file.
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE read_param

    USE parameters_2d, ONLY : n_vars , n_eqns 
    USE parameters_2d, ONLY : limiter
    USE parameters_2d, ONLY : bcW , bcE , bcS , bcN

    USE geometry_2d, ONLY : deposit , deposit_tot , erosion , erodible ,        &
         erosion_tot , x_quad , w_quad

    USE constitutive_2d, ONLY : rho_a_amb
    USE constitutive_2d, ONLY : rho_c_sub
    USE constitutive_2d, ONLY : kin_visc_c , sp_heat_c 

    USE constitutive_2d, ONLY : inv_pres , inv_rho_l , inv_rho_s , c_inv_rho_s

    USE constitutive_2d, ONLY : n_td2 
    USE constitutive_2d, ONLY : coeff_porosity
    USE constitutive_2d, ONLY : radiative_term_coeff , SBconst
    USE constitutive_2d, ONLY : convective_term_coeff

    USE constitutive_2d, ONLY : H_crit_rel

    USE init_2d, ONLY : erodible_init

    ! External procedures
    USE constitutive_2d, ONLY : mixt_var , eval_flux_coeffs

    USE geometry_2d, ONLY : lambertw0 , lambertwm1 , gaulegf


    IMPLICIT none

    NAMELIST / west_boundary_conditions / h_bcW , hu_bcW , hv_bcW ,             &
         alphas_bcW , halphas_bcW , T_bcW , alphag_bcW , halphag_bcW ,          &
         alphal_bcW , halphal_bcW

    NAMELIST / east_boundary_conditions / h_bcE , hu_bcE , hv_bcE ,             &
         alphas_bcE , halphas_bcE , T_bcE , alphag_bcE , halphag_bcE ,          &
         alphal_bcE , halphal_bcE

    NAMELIST / south_boundary_conditions / h_bcS , hu_bcS , hv_bcS ,            &
         alphas_bcS , halphas_bcS , T_bcS , alphag_bcS , halphag_bcS ,          &
         alphal_bcS , halphal_bcS

    NAMELIST / north_boundary_conditions / h_bcN , hu_bcN , hv_bcN ,            &
         alphas_bcN , halphas_bcN , T_bcN , alphag_bcN , halphag_bcN ,          &
         alphal_bcN , halphal_bcN

    NAMELIST / solid_transport_parameters / rho_s , diam_s , sp_heat_s ,        &
         erosion_coeff , erodible_porosity , settling_flag , T_erodible ,       &
         erodible_file , erodible_fract , alphastot_min ,                       &
         initial_erodible_thickness , erodible_deposit_flag

    NAMELIST / gas_transport_parameters / sp_heat_a, sp_gas_const_a, kin_visc_a,&
         pres , T_ambient , entrainment_flag , sp_heat_g , sp_gas_const_g ,     &
         liquid_vaporization_flag , water_level , gamma_steam


    REAL(wp) :: max_cfl

    LOGICAL :: tend1 
    CHARACTER(LEN=80) :: card

    INTEGER :: i_file

    INTEGER :: i_solid , j , k

    INTEGER :: dot_idx

    CHARACTER(LEN=3) :: check_file

    LOGICAL :: lexist

    CHARACTER(LEN=15) :: chara

    INTEGER :: ios

    REAL(wp) :: expA , expB , Tc

    CHARACTER(LEN=20) :: source_name
    REAL(wp) :: source_length
    
    REAL(wp), ALLOCATABLE :: qp_source(:)
    REAL(wp) :: red_grav
    REAL(wp) :: rho_c
    REAL(wp) :: rho_m
    REAL(wp) :: Ri
    LOGICAL :: sp_heat_flag
    REAL(wp) :: sp_heat_mix
    REAL(wp) :: mfr

    REAL(wp) :: pi_g

    REAL(wp) :: a_crit_rel

    REAL(wp), ALLOCATABLE :: inv_rho_g(:)
    REAL(wp) :: inv_rho_c
    REAL(wp) :: inv_rhom
    REAL(wp) :: sp_gas_const_c
    REAL(wp) :: xc
    REAL(wp) :: xs_tot
    REAL(wp), ALLOCATABLE :: shape_coeff(:)

    INTEGER :: iter_source , iter_max

    OPEN(input_unit,FILE=input_file,STATUS='old')

    ! ---------- READ run_parameters NAMELIST -----------------------------------
    READ(input_unit, run_parameters,IOSTAT=ios )

    IF ( ios .NE. 0 ) THEN

       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist RUN_PARAMETERS'
       WRITE(*,*) 'Please check the input file'
       STOP

    ELSE

       WRITE(*,*) 'Run name: ',run_name
       REWIND(input_unit)

    END IF

    IF ( (.NOT.output_cons_flag) .AND. (.NOT.output_esri_flag) .AND.            &
         (.NOT.output_phys_flag) ) dt_output = 2.0 * ( t_end - t_start ) 

    t_output = t_start + dt_output

    ! ---------- READ newrun_parameters NAMELIST --------------------------------
    READ(input_unit,newrun_parameters,IOSTAT=ios)

    IF ( ios .NE. 0 ) THEN

       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist NEWRUN_PARAMETERS'
       WRITE(*,*) 'Please check the input file'
       STOP

    ELSE

       REWIND(input_unit)

    END IF

    IF ( n_solid .LT. 0 ) THEN

       WRITE(*,*) 'ERROR: problem with namelist NEWRUN_PARAMETERS'
       WRITE(*,*) 'n_solid =' , n_solid
       WRITE(*,*) 'Please check the input file'
       STOP

    END IF

    IF ( n_add_gas .LT. 0 ) THEN

       WRITE(*,*) 'ERROR: problem with namelist NEWRUN_PARAMETERS'
       WRITE(*,*) 'n_add_gas =' , n_add_gas
       WRITE(*,*) 'Please check the input file'
       STOP

    END IF

    IF ( ( comp_cells_x .EQ. 1 ) .OR. ( comp_cells_y .EQ. 1 ) ) THEN

       IF ( verbose_level .GE. 0 ) WRITE(*,*) '----- 1D SIMULATION -----' 

    ELSE

       IF ( verbose_level .GE. 0 ) WRITE(*,*) '----- 2D SIMULATION -----' 

    END IF

    IF ( ( .NOT. liquid_flag ) .AND. ( .NOT. gas_flag ) ) THEN

       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist NEWRUN_PARAMETERS'
       WRITE(*,*) 'One of these parameters must be set to .TRUE.'
       WRITE(*,*) 'LIQUID_FLAG',liquid_flag
       WRITE(*,*) 'GAS_FLAG',liquid_flag
       WRITE(*,*) 'Please check the input file'
       STOP

    END IF

    ! ------- READ gas_transport_parameters NAMELIST --------------------------

    READ(input_unit, gas_transport_parameters,IOSTAT=ios)

    IF ( ios .NE. 0 ) THEN

       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist GAS_TRANSPORT_PARAMETERS'
       WRITE(*,*) 'Please check the input file'
       STOP

    ELSE

       REWIND(input_unit)

    END IF

    IF ( sp_heat_a .EQ. -1.0_wp ) THEN

       WRITE(*,*) 'ERROR: problem with namelist GAS_TRANSPORT_PARAMETERS'
       WRITE(*,*) 'SP_HEAT_a =' , sp_heat_a
       WRITE(*,*) 'Please check the input file'
       STOP

    END IF

    IF ( sp_gas_const_a .EQ. -1.0_wp ) THEN

       WRITE(*,*) 'ERROR: problem with namelist GAS_TRANSPORT_PARAMETERS'
       WRITE(*,*) 'SP_GAS_CONST_a =' , sp_gas_const_a
       WRITE(*,*) 'Please check the input file'
       STOP

    END IF

    IF ( kin_visc_a .EQ. -1.0_wp ) THEN

       WRITE(*,*) 'ERROR: problem with namelist GAS_TRANSPORT_PARAMETERS'
       WRITE(*,*) 'KIN_VISC_CONST_a =' , kin_visc_a
       WRITE(*,*) 'Please check the input file'
       STOP

    ELSE

       IF ( gas_flag ) THEN

          IF ( VERBOSE_LEVEL .GE. 0 ) THEN

             WRITE(*,*) 'CARRIER PHASE: gas'
             WRITE(*,*) 'Carrier phase kinematic viscosity:',kin_visc_a

          END IF
          kin_visc_c = kin_visc_a

       END IF

    END IF

    IF ( ANY(sp_heat_g(1:n_add_gas) .EQ. -1.0_wp ) ) THEN

       WRITE(*,*) 'ERROR: problem with namelist GAS_TRANSPORT_PARAMETERS'
       WRITE(*,*) 'SP_HEAT_G =' , sp_heat_g(1:n_solid)
       WRITE(*,*) 'Please check the input file'
       STOP

    END IF

    IF ( ANY(sp_gas_const_g(1:n_add_gas) .EQ. -1.0_wp ) ) THEN

       WRITE(*,*) 'ERROR: problem with namelist GAS_TRANSPORT_PARAMETERS'
       WRITE(*,*) 'SP_GAS_CONST_G =' , sp_gas_const_g(1:n_add_gas)
       WRITE(*,*) 'Please check the input file'
       STOP

    END IF

    alphag_bcW(1:n_add_gas)%flag = -1
    alphag_bcE(1:n_add_gas)%flag = -1
    alphag_bcS(1:n_add_gas)%flag = -1
    alphag_bcN(1:n_add_gas)%flag = -1

    halphag_bcW(1:n_add_gas)%flag = -1
    halphag_bcE(1:n_add_gas)%flag = -1
    halphag_bcS(1:n_add_gas)%flag = -1
    halphag_bcN(1:n_add_gas)%flag = -1

    IF ( pres .EQ. -1.0_wp ) THEN

       WRITE(*,*) 'ERROR: problem with namelist GAS_TRANSPORT_PARAMETERS'
       WRITE(*,*) 'pres =' , pres
       WRITE(*,*) 'Please check the input file'
       STOP

    ELSE

       inv_pres = 1.0_wp / pres

    END IF

    IF ( T_ambient .EQ. -1.0_wp ) THEN

       WRITE(*,*) 'ERROR: problem with namelist GAS_TRANSPORT_PARAMETERS'
       WRITE(*,*) 'T_ambient =' , T_ambient
       WRITE(*,*) 'Please check the input file'
       STOP

    END IF

    IF ( ( .NOT. gas_flag ) .AND. ( liquid_flag .AND. entrainment_flag ) ) THEN

       WRITE(*,*) 'ERROR: problem with namelist GAS_TRANSPORT_PARAMETERS'
       WRITE(*,*) 'LIQUID_FLAG',liquid_flag
       WRITE(*,*) 'ENTRAINMENT_FLAG =' , entrainment_flag
       WRITE(*,*) 'Please check the input file'
       STOP

    END IF

    rho_a_amb = pres / ( sp_gas_const_a * T_ambient )
    IF ( verbose_level .GE. 0 ) THEN

       WRITE(*,*) 'Ambient density = ',rho_a_amb,' (kg/m3)'

    END IF

    IF ( liquid_vaporization_flag ) THEN

       IF ( .NOT. gas_flag ) THEN

          WRITE(*,*) 'ERROR: problem with namelist GAS_TRANSPORT_PARAMETERS'
          WRITE(*,*) 'liquid_vaporization_flag =' , liquid_vaporization_flag
          WRITE(*,*) 'This flag can be set to .TRUE. only when gas_flag = T'
          WRITE(*,*) 'Please check the input file'
          STOP

       END IF

       IF ( n_add_gas .LT. 1 ) THEN

          WRITE(*,*) 'ERROR: problem with namelist GAS_TRANSPORT_PARAMETERS'
          WRITE(*,*) 'liquid_vaporization_flag =' , liquid_vaporization_flag
          WRITE(*,*) 'This flag can be set to .TRUE. only when n_add_gas > 1 '
          WRITE(*,*) 'Please check the input file'
          STOP

       END IF

       IF ( ( gamma_steam .LT. 0.0_wp ) .OR. ( gamma_steam .GT. 1.0_wp ) ) THEN

          WRITE(*,*) 'ERROR: problem with namelist GAS_TRANSPORT_PARAMETERS'
          WRITE(*,*) 'gamma_steam =' , gamma_steam
          WRITE(*,*) 'Specify a value between 0 and 1'
          WRITE(*,*) 'Please check the input file'
          STOP

       END IF

       IF ( water_level .EQ. -1.0E7_wp ) THEN

          WRITE(*,*) 'ERROR: problem with namelist GAS_TRANSPORT_PARAMETERS'
          WRITE(*,*) 'water_level =' , water_level
          WRITE(*,*) 'Specify a value when LIQUID_VAPORIZATION_FLAG = T'
          WRITE(*,*) 'Please check the input file'
          STOP

       ELSE

          IF ( .NOT. liquid_vaporization_flag ) THEN

             WRITE(*,*) 'ERROR: problem with namelist GAS_TRANSPORT_PARAMETERS'
             WRITE(*,*) 'water_level =' , water_level
             WRITE(*,*) 'liquid_vaporization_flag =' , liquid_vaporization_flag
             WRITE(*,*) 'Specify a value only when GAS_VAPORIZATION_FLAG = T'
             WRITE(*,*) 'Please check the input file'
             STOP

          END IF

       END IF

    ELSE

       IF ( gamma_steam .NE. -1.0_wp ) THEN

          WRITE(*,*) 'ERROR: problem with namelist GAS_TRANSPORT_PARAMETERS'
          WRITE(*,*) 'gamma_steam =' , gamma_steam
          WRITE(*,*) 'Specify a value only when GAS_VAPORIZATION_FLAG = T'
          WRITE(*,*) 'Please check the input file'
          STOP

       END IF

    END IF

    ! ------- READ liquid_transport_parameters NAMELIST -------------------------

    n_vars = 4

    IF ( gas_flag ) n_vars = n_vars + n_add_gas

    IF ( liquid_flag ) THEN

       IF ( gas_flag ) n_vars = n_vars + 1

       READ(input_unit, liquid_transport_parameters,IOSTAT=ios)

       IF ( ios .NE. 0 ) THEN

          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist LIQUID_TRANSPORT_PARAMETERS'
          WRITE(*,liquid_transport_parameters)
          WRITE(*,*) 'Please check the input file'
          STOP

       ELSE

          REWIND(input_unit)

       END IF

       IF ( sp_heat_l .EQ. -1.0_wp ) THEN

          WRITE(*,*) 'ERROR: problem with namelist LIQUID_TRANSPORT_PARAMETERS'
          WRITE(*,*) 'SP_HEAT_L =' , sp_heat_l
          WRITE(*,*) 'Please check the input file'
          STOP

       END IF

       IF ( rho_l .EQ. -1.0_wp ) THEN

          WRITE(*,*) 'ERROR: problem with namelist LIQUID_TRANSPORT_PARAMETERS'
          WRITE(*,*) 'RHO_L =' , rho_l
          WRITE(*,*) 'Please check the input file'
          STOP

       ELSE

          inv_rho_l = 1.0_wp / rho_l

       END IF

       IF ( loss_rate .LT. 0.0_wp ) THEN

          WRITE(*,*) 'ERROR: problem with namelist LIQUID_TRANSPORT_PARAMETERS'
          WRITE(*,*) 'LOSS_RATE =' , loss_rate
          WRITE(*,*) 'Please check the input file'
          STOP

       END IF


       READ(input_unit, rheology_parameters,IOSTAT=ios)
       REWIND(input_unit)

       IF ( .NOT. gas_flag ) THEN

          IF ( kin_visc_l .EQ. -1.0_wp ) THEN

             IF ( ( RHEOLOGY_MODEL .NE. 4 ) .AND. ( RHEOLOGY_MODEL .NE. 3 ) ) THEN

                WRITE(*,*) 'ERROR: problem with namelist LIQUID_TRANSPORT_PARAMETERS'
                WRITE(*,*) 'KIN_VISC_L =' , kin_visc_l
                WRITE(*,*) 'Please check the input file'
                STOP

             END IF

          ELSE

             IF ( RHEOLOGY_MODEL .EQ. 4 ) THEN

                WRITE(*,*) 'ERROR: problem with namelist LIQUID_TRANSPORT_PARAMETERS'
                WRITE(*,*) 'KIN_VISC_L =' , kin_visc_l
                WRITE(*,*) 'Viscosity already is computed by REHOLOGY MODEL=4' 
                WRITE(*,*) 'Please check the input file'
                STOP

             END IF

             IF ( RHEOLOGY_MODEL .EQ. 3 ) THEN

                WRITE(*,*) 'ERROR: problem with namelist LIQUID_TRANSPORT_PARAMETERS'
                WRITE(*,*) 'KIN_VISC_L =' , kin_visc_l
                WRITE(*,*) 'Viscosity already is computed by REHOLOGY MODEL=3' 
                WRITE(*,*) 'Please check the input file'
                STOP

             END IF

          END IF

          IF ( verbose_level .GE. 0 ) THEN

             WRITE(*,*) 'CARRIER PHASE: liquid'
             WRITE(*,*) 'Carrier phase kinematic viscosity:',kin_visc_l

          END IF

          kin_visc_c = kin_visc_l
          sp_heat_c = sp_heat_l

       END IF

    END IF

    ! ------- READ solid_transport_parameters NAMELIST --------------------------

    ALLOCATE( erodible_init(comp_cells_x,comp_cells_y) )
    ALLOCATE( erodible( n_solid , comp_cells_x , comp_cells_y ) )

    erodible_init(:,:) = 0.0E+0_wp
    erodible(1:n_solid,1:comp_cells_x,1:comp_cells_y) = 0.0_wp
    
    IF ( n_solid .GE. 1 ) THEN

       READ(input_unit, solid_transport_parameters,IOSTAT=ios)

       IF ( ios .NE. 0 ) THEN

          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist SOLID_TRANSPORT_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          WRITE(*,solid_transport_parameters) 
          STOP

       ELSE

          REWIND(input_unit)

       END IF

       IF ( ANY(rho_s(1:n_solid) .EQ. -1.0_wp ) ) THEN

          WRITE(*,*) 'ERROR: problem with namelist SOLID_TRANSPORT_PARAMETERS'
          WRITE(*,*) 'RHO_s =' , rho_s(1:n_solid)
          WRITE(*,*) 'Please check the input file'
          STOP

       END IF

       IF ( ANY(diam_s(1:n_solid) .EQ. -1.0_wp ) ) THEN

          WRITE(*,*) 'ERROR: problem with namelist SOLID_TRANSPORT_PARAMETERS'
          WRITE(*,*) 'DIAM_s =' , diam_s(1:n_solid)
          WRITE(*,*) 'Please check the input file'
          STOP

       END IF

       IF ( ANY(sp_heat_s(1:n_solid) .EQ. -1.0_wp ) ) THEN

          WRITE(*,*) 'ERROR: problem with namelist SOLID_TRANSPORT_PARAMETERS'
          WRITE(*,*) 'SP_HEAT_S =' , sp_heat_s(1:n_solid)
          WRITE(*,*) 'Please check the input file'
          STOP

       END IF

       IF ( erosion_coeff .LT. 0.0_wp ) THEN

          WRITE(*,*) 'ERROR: problem with namelist SOLID_TRANSPORT_PARAMETERS'
          WRITE(*,*) 'EROSION_COEFF =' , erosion_coeff
          WRITE(*,*) 'Please check the input file'
          STOP

       ELSE

          IF ( erosion_coeff .EQ. 0.0_wp ) THEN

             erodible_porosity = 0.0_wp
             erodible_fract(1:n_solid) = 1.0_wp / n_solid
             T_erodible = 300.0_wp
             !subtract_init_flag = .FALSE.

          ELSE

             IF ( ( erodible_porosity .LT. 0.0_wp ) .OR.                           &
                  ( erodible_porosity .GT. 1.0_wp ) ) THEN

                WRITE(*,*) 'ERROR: problem with namelist SOLID_TRANSPORT_PARAMETERS'
                WRITE(*,*) 'erodible_porosity =' , erodible_porosity
                WRITE(*,*) 'Please check the input file'
                STOP

             END IF

             coeff_porosity = erodible_porosity / ( 1.0_wp - erodible_porosity )

             read_erodible_fract:IF ( n_solid .EQ. 1 ) THEN

                erodible_fract(1) = 1.0_wp

             ELSE

                IF ( ANY(erodible_fract(1:n_solid) .LT. 0.0_wp ) ) THEN

                   WRITE(*,*) 'ERROR: problem with namelist ',                     &
                        'SOLID_TRANSPORT_PARAMETERS'
                   WRITE(*,*) 'ERODIBLE_FRACT =' , erodible_fract(1:n_solid)
                   WRITE(*,*) 'Please check the input file'
                   STOP

                ELSE

                   IF ( SUM(erodible_fract(1:n_solid)) .NE. 1.0_wp ) THEN

                      WRITE(*,*) 'WARNING: sum of ERODIBLE_FRACT not 1:',          &
                           SUM(erodible_fract(1:n_solid))

                      erodible_fract(1:n_solid) = erodible_fract(1:n_solid) /      &
                           SUM(erodible_fract(1:n_solid) )

                   END IF

                END IF

                WRITE(*,*) 'Absolute fractions of solide phases in erodible layer:'
                WRITE(*,*) erodible_fract * ( 1.0_wp-erodible_porosity )

             END IF read_erodible_fract

             IF ( T_erodible .LT. 0.0_wp ) THEN

                WRITE(*,*) 'ERROR: problem with namelist SOLID_TRANSPORT_PARAMETERS'
                WRITE(*,*) 'T_erodible =' , T_erodible
                WRITE(*,*) 'Please check the input file'
                STOP

             END IF

          END IF

       END IF

       IF ( gas_flag ) THEN

          rho_c_sub = pres / ( sp_gas_const_a * T_erodible )

       ELSE

          rho_c_sub = rho_l

       END IF

       IF ( erodible_deposit_flag ) THEN

          IF ( t_erodible .GT. 0.0_wp ) THEN

             WRITE(*,*) 'WARNING: t_erodible NOT USED'
             WRITE(*,*) 'Temperature of erodible layer is set equal to'
             WRITE(*,*) 'that of the flow'


          END IF

       END IF

       check_erodible_file:IF ( TRIM(erodible_file) .EQ. '' ) THEN
          
          IF ( initial_erodible_thickness .GE. 0.0_wp ) THEN
             
             IF( erosion_coeff .EQ. 0.0_wp ) THEN
                
                WRITE(*,*) 'WARNING: erodible_file not used'
                WRITE(*,*) 'erosion_coeff = ', erosion_coeff
                
                erodible_init(:,:) = 0.0E+0_wp

             ELSE

                erodible_init(:,:) = initial_erodible_thickness
                
                WRITE(*,*) 'Initial thickness of erodible layer',               &
                     initial_erodible_thickness

             END IF

          ELSE

             erodible_init(:,:) = 0.0E+0_wp

             IF ( verbose_level .GE. 0.0_wp ) THEN

                WRITE(*,*)
                WRITE(*,*) 'WARNING: no file defined for erobile layer'
                WRITE(*,*) 'WARNING: no initial thickness for erobile layer'
                WRITE(*,*) 'Initial erodible thickness set to 0'
                WRITE(*,*)

             END IF

          END IF
          
       ELSE

          IF ( initial_erodible_thickness .GE. 0.0_wp ) THEN

             WRITE(*,*) 'WARNING: initial_erodible_thicknes not used'

          END IF

          IF( erosion_coeff .EQ. 0.0_wp ) THEN

             WRITE(*,*) 'WARNING: erodible_file not used'
             WRITE(*,*) 'erosion_coeff = ', erosion_coeff

             erodible_init(:,:) = 0.0E+0_wp

          ELSE

             IF ( verbose_level .GE. 0.0_wp ) THEN

                WRITE(*,*) 'Maximum thick. for erosion read from : ',           &
                     TRIM(erodible_file)

             END IF

             CALL read_erodible

          END IF

       END IF check_erodible_file

       DO i_solid=1,n_solid

          erodible(i_solid,:,:) = erodible_fract(i_solid) *               &
               ( 1.0_wp - erodible_porosity ) * erodible_init(:,:)
          
       END DO
       
    END IF

    n_vars = n_vars + n_solid
    n_eqns = n_vars

    WRITE(*,*) 'Model variables = ',n_vars

    alphas_bcW(1:n_solid)%flag = -1
    alphas_bcE(1:n_solid)%flag = -1
    alphas_bcS(1:n_solid)%flag = -1
    alphas_bcN(1:n_solid)%flag = -1

    halphas_bcW(1:n_solid)%flag = -1
    halphas_bcE(1:n_solid)%flag = -1
    halphas_bcS(1:n_solid)%flag = -1
    halphas_bcN(1:n_solid)%flag = -1

    ALLOCATE( bcW(n_vars) , bcE(n_vars) , bcS(n_vars) , bcN(n_vars) )

    bcW(1:n_vars)%flag = -1
    bcE(1:n_vars)%flag = -1
    bcS(1:n_vars)%flag = -1
    bcN(1:n_vars)%flag = -1

    ALLOCATE( inv_rho_s(n_solid) )
    ALLOCATE( c_inv_rho_s(n_solid) )

    ALLOCATE( alphas_init(n_solid) )

    inv_rho_s(1:n_solid) = 1.0_wp / rho_s(1:n_solid)

    DO i_solid=1,n_solid

       c_inv_rho_s(i_solid) = CMPLX(inv_rho_s(i_solid),0.0_wp,wp)

    END DO

    ALLOCATE( deposit( comp_cells_x , comp_cells_y , n_solid ) )    
    deposit(1:comp_cells_x,1:comp_cells_y,1:n_solid ) = 0.0_wp

    ALLOCATE( deposit_tot( comp_cells_x , comp_cells_y ) )    
    deposit_tot(1:comp_cells_x,1:comp_cells_y) = 0.0_wp

    ALLOCATE( erosion( comp_cells_x , comp_cells_y , n_solid ) )
    erosion(1:comp_cells_x,1:comp_cells_y,1:n_solid ) = 0.0_wp

    ALLOCATE( erosion_tot( comp_cells_x , comp_cells_y ) )
    erosion_tot(1:comp_cells_x,1:comp_cells_y) = 0.0_wp

    IF ( ( alphastot_min .LT. 0.0_wp ) .OR. ( alphastot_min .GE. 1.0_wp ) ) THEN

       WRITE(*,*) 'ERROR: problem with namelist SOLID_TRANSPORT_PARAMETERS'
       WRITE(*,*) 'alphastot_min should be >0 and <1'
       WRITE(*,*) 'alphastot_min =',alphastot_min
       STOP

    ELSEIF ( alphastot_min .EQ. 0.0_wp) THEN

       WRITE(*,*) 'WARNING: minimum total solid fraction = 0'

    END IF


    IF ( restart ) THEN

       ! ---------- READ restart_parameters NAMELIST ----------------------------
       READ(input_unit,restart_parameters,IOSTAT=ios)

       IF ( ios .NE. 0 ) THEN

          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist RESTART_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          STOP

       ELSE

          IF ( n_restart_files .LE. 0 ) THEN

             WRITE(*,*) 'n_restart_files must be >= 1'
             WRITE(*,*) 'n_restart_files',n_restart_files
             WRITE(*,*) 'ERROR: problem with namelist RESTART_PARAMETERS'
             WRITE(*,*) 'Please check the input file'
             STOP             

          ELSE

             DO i_file=1,n_restart_files

                IF ( TRIM(restart_files(i_file)) .EQ. '' ) THEN

                   WRITE(*,*) 'ERROR: problem with namelist RESTART_PARAMETERS'
                   WRITE(*,*) 'n_restart_files',n_restart_files
                   WRITE(*,*) 'restart_files ',restart_files
                   WRITE(*,*) 'Please check the input file'
                   STOP                    

                END IF

             END DO

             DO i_file=n_restart_files+1,10

                IF ( TRIM(restart_files(i_file)) .NE. '' ) THEN

                   WRITE(*,*) 'ERROR: problem with namelist RESTART_PARAMETERS'
                   WRITE(*,*) 'n_restart_files',n_restart_files
                   WRITE(*,*) 'restart_files',restart_files
                   WRITE(*,*) 'Please check the input file'
                   STOP                    

                END IF

             END DO

             restart_file = restart_files(1)

          END IF

          dot_idx = SCAN(restart_file, ".", .TRUE.)

          check_file = restart_file(dot_idx+1:dot_idx+3)

          IF ( check_file .EQ. 'asc' ) THEN

             IF ( ( ANY(sed_vol_perc(1:n_solid) .LT. 0.0_wp ) ) .OR.            &
                  ( ANY(sed_vol_perc(1:n_solid) .GT. 100.0_wp ) ) .OR.          &
                  ( SUM(sed_vol_perc(1:n_solid)) .GT. 100.0_wp) ) THEN

                WRITE(*,*) 'ERROR: problem with namelist RESTART_PARAMETERS'
                WRITE(*,*) 'SED_VOL_PERC =' , sed_vol_perc(1:n_solid)
                STOP

             END IF

             alphas_init(1:n_solid) = 1.0E-2_wp * sed_vol_perc(1:n_solid)

             IF ( alphastot_min .GE. SUM(alphas_init(1:n_solid)) ) THEN

                WRITE(*,*) 'IOSTAT=',ios
                WRITE(*,*) 'ERROR: problem with namelist SOLID_TRANSPORT_PARAMETERS'
                WRITE(*,*) 'alphastot_min should be < SUM(0.01*sed_vol_perc)'
                WRITE(*,*) 'alphastot_min =',alphastot_min
                WRITE(*,*) '0.01*SUM(sed_vol_perc) =',0.01_wp*sed_vol_perc(1:n_solid)
                STOP

             END IF

             IF ( verbose_level .GE. 0 ) THEN

                WRITE(*,*) 'INITIAL VOLUME FRACTION OF SOLIDS:', alphas_init

             END IF

             REWIND(input_unit)

             IF ( T_init*T_ambient .EQ. 0.0_wp ) THEN

                WRITE(*,*) 'ERROR: problem with namelist RESTART_PARAMETERS'
                WRITE(*,*) 'T_init=',T_init
                WRITE(*,*) 'T_ambient=',T_ambient
                WRITE(*,*) 'Add the variables to the namelist RESTART_PARAMETERS'
                STOP

             END IF
             
          END IF

       END IF

    END IF

    ! ------- READ numeric_parameters NAMELIST ----------------------------------

    READ(input_unit,numeric_parameters)

    IF ( ios .NE. 0 ) THEN

       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist NUMERIC_PARAMETERS'
       WRITE(*,*) 'Please check the input file'
       STOP

    ELSE

       REWIND(input_unit)

    END IF

    IF ( ( solver_scheme .NE. 'LxF' ) .AND. ( solver_scheme .NE. 'KT' ) .AND.   &
         ( solver_scheme .NE. 'GFORCE' ) .AND. ( solver_scheme .NE. 'UP' ) ) THEN

       WRITE(*,*) 'WARNING: no correct solver scheme selected',solver_scheme
       WRITE(*,*) 'Chose between: LxF, GFORCE or KT'
       STOP

    END IF

    IF  ( ( solver_scheme.EQ.'LxF' ) .OR. ( solver_scheme.EQ.'GFORCE' ) ) THEN 

       max_cfl = 1.0

    ELSE

       IF ( ( comp_cells_x .EQ. 1 ) .OR. ( comp_cells_y .EQ. 1 ) ) THEN

          max_cfl = 0.50_wp

       ELSE

          max_cfl = 0.25_wp

       END IF

    END IF


    IF ( ( cfl .GT. max_cfl ) .OR. ( cfl .LT. 0.0_wp ) ) THEN

       WRITE(*,*) 'WARNING: wrong value of cfl ',cfl
       WRITE(*,*) 'Choose a value between 0.0 and ',max_cfl
       READ(*,*)

    END IF

    IF ( verbose_level .GE. 1 ) WRITE(*,*) 'Limiters',limiter(1:n_vars)

    limiter(n_vars+1) = limiter(2)
    limiter(n_vars+2) = limiter(3)

    IF ( ( MAXVAL(limiter(1:n_vars)) .GT. 7 ) .OR.                              &
         ( MINVAL(limiter(1:n_vars)) .LT. 0 ) ) THEN

       WRITE(*,*) 'WARNING: wrong limiter ',limiter(1:n_vars)
       WRITE(*,*) 'Choose among: none, minmod,superbee,van_leer'
       STOP         

    END IF

    IF ( verbose_level .GE. 0 ) THEN

       IF ( alpha_flag ) THEN

          WRITE(*,*) 'Linear reconstruction and b. c. applied to variables:'
          WRITE(*,*) 'h,hu,hv,T,alphas'

       ELSE

          WRITE(*,*) 'Linear reconstruction and b. c. applied to variables:'
          WRITE(*,*) 'h,hu,hv,T,halphas'

       END IF

    END IF

    IF ( ( reconstr_coeff .GT. 1.0_wp ).OR.( reconstr_coeff .LT. 0.0_wp ) ) THEN

       WRITE(*,*) 'WARNING: wrong value of reconstr_coeff ',reconstr_coeff
       WRITE(*,*) 'Change the value between 0.0 and 1.0 in the input file'
       READ(*,*)

    END IF

    ! ------- READ boundary_conditions NAMELISTS --------------------------------

    IF ( COMP_CELLS_X .GT. 1 ) THEN

       ! --------- West boundary conditions -------------------------------------

       READ(input_unit,west_boundary_conditions,IOSTAT=ios)

       IF ( ios .NE. 0 ) THEN

          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist WEST_BOUNDARY_CONDITIONS'
          WRITE(*,*) 'Please check the input file'
          WRITE(*,west_boundary_conditions)
          STOP

       ELSE

          REWIND(input_unit)

       END IF

       IF ( ( h_bcW%flag .EQ. -1 ) ) THEN 

          WRITE(*,*) 'ERROR: problem with namelist WEST_BOUNDARY_CONDITIONS'
          WRITE(*,*) 'B.C. for h not set properly'
          WRITE(*,*) 'Please check the input file'
          STOP

       END IF

       IF ( comp_cells_x .GT. 1 ) THEN

          IF ( hu_bcW%flag .EQ. -1 ) THEN

             WRITE(*,*) 'ERROR: problem with namelist WEST_BOUNDARY_CONDITIONS'
             WRITE(*,*) 'B.C. for hu not set properly'             
             WRITE(*,*) 'Please check the input file'
             STOP

          ELSE

             ! hu_bcW%flag = 1
             ! hu_bcW%value = 0.0_wp

          END IF

       END IF

       IF ( comp_cells_y .GT. 1 ) THEN

          IF ( hv_bcW%flag .EQ. -1 ) THEN 

             WRITE(*,*) 'ERROR: problem with namelist WEST_BOUNDARY_CONDITIONS'
             WRITE(*,*) 'B.C. for hv not set properly'             
             WRITE(*,*) 'Please check the input file'
             STOP

          ELSE

             hv_bcW%flag = 1
             hv_bcW%value = 0.0_wp

          END IF

       END IF

       IF ( alpha_flag ) THEN

          IF ( ANY(alphas_bcW(1:n_solid)%flag .EQ. -1 ) ) THEN 

             WRITE(*,*) 'ERROR: problem with namelist WEST_BOUNDARY_CONDITIONS'
             WRITE(*,*) 'B.C. for sediment conentration not set properly'
             WRITE(*,*) 'Please check the input file'
             WRITE(*,*) 'alphas_bcW'
             WRITE(*,*) alphas_bcW(1:n_solid)
             STOP

          END IF

          IF ( ANY(alphag_bcW(1:n_add_gas)%flag .EQ. -1 ) ) THEN 

             WRITE(*,*) 'ERROR: problem with namelist WEST_BOUNDARY_CONDITIONS'
             WRITE(*,*) 'B.C. for additional gas components not set properly'
             WRITE(*,*) 'Please check the input file'
             WRITE(*,*) 'alphag_bcW'
             WRITE(*,*) alphag_bcW(1:n_add_gas)
             STOP

          END IF

          IF ( gas_flag .AND. liquid_flag ) THEN

             IF ( alphal_bcW%flag .EQ. -1 ) THEN 

                WRITE(*,*) 'ERROR: problem with namelist WEST_BOUNDARY_CONDITIONS'
                WRITE(*,*) 'B.C. for additional gas components not set properly'
                WRITE(*,*) 'Please check the input file'
                WRITE(*,*) 'alphal_bcW'
                WRITE(*,*) alphal_bcW
                STOP

             END IF

          END IF

       ELSE

          IF ( ANY(halphas_bcW(1:n_solid)%flag .EQ. -1 ) ) THEN 

             WRITE(*,*) 'ERROR: problem with namelist WEST_BOUNDARY_CONDITIONS'
             WRITE(*,*) 'B.C. for sediment conentration not set properly'
             WRITE(*,*) 'Please check the input file'
             WRITE(*,*) 'halphas_bcW'
             WRITE(*,*) halphas_bcW(1:n_solid)
             STOP

          END IF

          IF ( ANY(halphag_bcW(1:n_add_gas)%flag .EQ. -1 ) ) THEN 

             WRITE(*,*) 'ERROR: problem with namelist WEST_BOUNDARY_CONDITIONS'
             WRITE(*,*) 'B.C. for additional gas components not set properly'
             WRITE(*,*) 'Please check the input file'
             WRITE(*,*) 'halphag_bcW'
             WRITE(*,*) halphag_bcW(1:n_add_gas)
             STOP

          END IF

          IF ( gas_flag .AND. liquid_flag ) THEN

             IF ( halphal_bcW%flag .EQ. -1 ) THEN 

                WRITE(*,*) 'ERROR: problem with namelist WEST_BOUNDARY_CONDITIONS'
                WRITE(*,*) 'B.C. for additional gas components not set properly'
                WRITE(*,*) 'Please check the input file'
                WRITE(*,*) 'halphal_bcW'
                WRITE(*,*) halphal_bcW
                STOP

             END IF

          END IF


       END IF

       IF ( T_bcW%flag .EQ. -1 ) THEN 

          WRITE(*,*) 'ERROR: problem with namelist WEST_BOUNDARY_CONDITIONS'
          WRITE(*,*) 'B.C. for temperature not set properly'
          WRITE(*,*) 'Please check the input file'
          STOP

       END IF

       ! set the approriate boundary conditions

       bcW(1) = h_bcW
       bcW(2) = hu_bcW 
       bcW(3) = hv_bcW 


       ! ------------- East boundary conditions --------------------------------

       READ(input_unit,east_boundary_conditions,IOSTAT=ios)

       IF ( ios .NE. 0 ) THEN

          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist EAST_BOUNDARY_CONDITIONS'
          WRITE(*,*) 'Please check the input file'
          STOP

       ELSE

          REWIND(input_unit)

       END IF

       IF ( ( h_bcE%flag .EQ. -1 ) ) THEN 

          WRITE(*,*) 'ERROR: problem with namelist EAST_BOUNDARY_CONDITIONS'
          WRITE(*,*) 'B.C. for h not set properly'
          WRITE(*,*) 'Please check the input file'
          STOP

       END IF

       IF ( comp_cells_x .GT. 1 ) THEN

          IF ( hu_bcE%flag .EQ. -1 ) THEN

             WRITE(*,*) 'ERROR: problem with namelist EAST_BOUNDARY_CONDITIONS'
             WRITE(*,*) 'B.C. for hu not set properly'             
             WRITE(*,*) 'Please check the input file'
             STOP

          END IF

       ELSE

          hu_bcE%flag = 1
          hu_bcE%value = 0.0_wp

       END IF

       IF ( comp_cells_y .GT. 1 ) THEN

          IF ( hv_bcE%flag .EQ. -1 ) THEN 

             WRITE(*,*) 'ERROR: problem with namelist EAST_BOUNDARY_CONDITIONS'
             WRITE(*,*) 'B.C. for hv not set properly'             
             WRITE(*,*) 'Please check the input file'
             STOP

          END IF

       ELSE

          hv_bcE%flag = 1
          hv_bcE%value = 0.0_wp


       END IF

       IF ( alpha_flag ) THEN

          IF ( ANY(alphas_bcE(1:n_solid)%flag .EQ. -1 ) ) THEN 

             WRITE(*,*) 'ERROR: problem with namelist EAST_BOUNDARY_CONDITIONS'
             WRITE(*,*) 'B.C. for sediment concentration not set properly'
             WRITE(*,*) 'Please check the input file'
             WRITE(*,*) 'alphas_bcE'
             WRITE(*,*) alphas_bcE(1:n_solid)
             STOP

          END IF

          IF ( ANY(alphag_bcE(1:n_add_gas)%flag .EQ. -1 ) ) THEN 

             WRITE(*,*) 'ERROR: problem with namelist EAST_BOUNDARY_CONDITIONS'
             WRITE(*,*) 'B.C. for additional gas components not set properly'
             WRITE(*,*) 'Please check the input file'
             WRITE(*,*) 'alphag_bcE'
             WRITE(*,*) alphag_bcE(1:n_add_gas)
             STOP

          END IF

          IF ( gas_flag .AND. liquid_flag ) THEN

             IF ( alphal_bcE%flag .EQ. -1 ) THEN 

                WRITE(*,*) 'ERROR: problem with namelist EAST_BOUNDARY_CONDITIONS'
                WRITE(*,*) 'B.C. for additional gas components not set properly'
                WRITE(*,*) 'Please check the input file'
                WRITE(*,*) 'alphal_bcE'
                WRITE(*,*) alphal_bcE
                STOP

             END IF

          END IF

       ELSE

          IF ( ANY(halphas_bcE(1:n_solid)%flag .EQ. -1 ) ) THEN 

             WRITE(*,*) 'ERROR: problem with namelist EAST_BOUNDARY_CONDITIONS'
             WRITE(*,*) 'B.C. for sediment concentration not set properly'
             WRITE(*,*) 'Please check the input file'
             WRITE(*,*) 'halphas_bcE'
             WRITE(*,*) halphas_bcE(1:n_solid)
             STOP

          END IF

          IF ( ANY(halphag_bcE(1:n_add_gas)%flag .EQ. -1 ) ) THEN 

             WRITE(*,*) 'ERROR: problem with namelist EAST_BOUNDARY_CONDITIONS'
             WRITE(*,*) 'B.C. for additional gas components not set properly'
             WRITE(*,*) 'Please check the input file'
             WRITE(*,*) 'halphag_bcE'
             WRITE(*,*) halphag_bcE(1:n_add_gas)
             STOP

          END IF

          IF ( gas_flag .AND. liquid_flag ) THEN

             IF ( halphal_bcE%flag .EQ. -1 ) THEN 

                WRITE(*,*) 'ERROR: problem with namelist EAST_BOUNDARY_CONDITIONS'
                WRITE(*,*) 'B.C. for additional gas components not set properly'
                WRITE(*,*) 'Please check the input file'
                WRITE(*,*) 'halphal_bcE'
                WRITE(*,*) halphal_bcE
                STOP

             END IF

          END IF

       END IF

       IF ( T_bcE%flag .EQ. -1 ) THEN 

          WRITE(*,*) 'ERROR: problem with namelist EAST_BOUNDARY_CONDITIONS'
          WRITE(*,*) 'B.C. for temperature not set properly'
          WRITE(*,*) 'Please check the input file'
          STOP

       END IF

       bcE(1) = h_bcE 
       bcE(2) = hu_bcE 
       bcE(3) = hv_bcE 

    END IF

    IF ( comp_cells_y .GT. 1 ) THEN

       ! --------------- South boundary conditions ------------------------------

       READ(input_unit,south_boundary_conditions,IOSTAT=ios)

       IF ( ios .NE. 0 ) THEN

          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist SOUTH_BOUNDARY_CONDITIONS'
          WRITE(*,*) 'Please check the input file'
          STOP

       ELSE

          REWIND(input_unit)

       END IF

       IF ( ( h_bcS%flag .EQ. -1 ) ) THEN 

          WRITE(*,*) 'ERROR: problem with namelist SOUTH_BOUNDARY_CONDITIONS'
          WRITE(*,*) 'B.C. for h not set properly'
          WRITE(*,*) 'Please check the input file'
          STOP

       END IF

       IF ( comp_cells_x .GT. 1 ) THEN

          IF ( hu_bcS%flag .EQ. -1 ) THEN 

             WRITE(*,*) 'ERROR: problem with namelist SOUTH_BOUNDARY_CONDITIONS'
             WRITE(*,*) 'B.C. for hu not set properly'             
             WRITE(*,*) 'Please check the input file'
             STOP

          END IF

       ELSE

          hu_bcS%flag = 1
          hu_bcS%value = 0.0_wp

       END IF

       IF ( comp_cells_y .GT. 1 ) THEN

          IF ( hv_bcS%flag .EQ. -1 ) THEN 

             WRITE(*,*) 'ERROR: problem with namelist SOUTH_BOUNDARY_CONDITIONS'
             WRITE(*,*) 'B.C. for hv not set properly'             
             WRITE(*,*) 'Please check the input file'
             STOP

          END IF

       ELSE

          hv_bcS%flag = 1
          hv_bcS%value = 0.0_wp

       END IF

       IF ( alpha_flag ) THEN

          IF ( ANY(alphas_bcS(1:n_solid)%flag .EQ. -1 ) ) THEN 

             WRITE(*,*) 'ERROR: problem with namelist SOUTH_BOUNDARY_CONDITIONS'
             WRITE(*,*) 'B.C. for sediment concentrations not set properly'
             WRITE(*,*) 'Please check the input file'
             WRITE(*,*) 'alphas_bcS'
             WRITE(*,*) alphas_bcS(1:n_solid)
             STOP

          END IF

          IF ( ANY(alphag_bcS(1:n_add_gas)%flag .EQ. -1 ) ) THEN 

             WRITE(*,*) 'ERROR: problem with namelist SOUTH_BOUNDARY_CONDITIONS'
             WRITE(*,*) 'B.C. for additional gas components not set properly'
             WRITE(*,*) 'Please check the input file'
             WRITE(*,*) 'alphag_bcS'
             WRITE(*,*) alphag_bcS(1:n_add_gas)
             STOP

          END IF

          IF ( gas_flag .AND. liquid_flag ) THEN

             IF ( alphal_bcS%flag .EQ. -1 ) THEN 

                WRITE(*,*) 'ERROR: problem with namelist SOUTH_BOUNDARY_CONDITIONS'
                WRITE(*,*) 'B.C. for additional gas components not set properly'
                WRITE(*,*) 'Please check the input file'
                WRITE(*,*) 'alphal_bcS'
                WRITE(*,*) alphal_bcS
                STOP

             END IF

          END IF

       ELSE

          IF ( ANY(halphas_bcS(1:n_solid)%flag .EQ. -1 ) ) THEN 

             WRITE(*,*) 'ERROR: problem with namelist SOUTH_BOUNDARY_CONDITIONS'
             WRITE(*,*) 'B.C. for sediment concentrations not set properly'
             WRITE(*,*) 'Please check the input file'
             WRITE(*,*) 'halphas_bcS'
             WRITE(*,*) halphas_bcS(1:n_solid)
             STOP

          END IF


          IF ( ANY(halphag_bcS(1:n_add_gas)%flag .EQ. -1 ) ) THEN 

             WRITE(*,*) 'ERROR: problem with namelist SOUTH_BOUNDARY_CONDITIONS'
             WRITE(*,*) 'B.C. for additional gas components not set properly'
             WRITE(*,*) 'Please check the input file'
             WRITE(*,*) 'halphag_bcS'
             WRITE(*,*) halphag_bcS(1:n_add_gas)
             STOP

          END IF

          IF ( gas_flag .AND. liquid_flag ) THEN

             IF ( halphal_bcS%flag .EQ. -1 ) THEN 

                WRITE(*,*) 'ERROR: problem with namelist SOUTH_BOUNDARY_CONDITIONS'
                WRITE(*,*) 'B.C. for additional gas components not set properly'
                WRITE(*,*) 'Please check the input file'
                WRITE(*,*) 'halphal_bcS'
                WRITE(*,*) halphal_bcS
                STOP

             END IF

          END IF

       END IF

       IF ( T_bcS%flag .EQ. -1 ) THEN 

          WRITE(*,*) 'ERROR: problem with namelist SOUTH_BOUNDARY_CONDITIONS'
          WRITE(*,*) 'B.C. for temperature not set properly'
          WRITE(*,*) 'Please check the input file'
          STOP

       END IF

       bcS(1) = h_bcS 
       bcS(2) = hu_bcS 
       bcS(3) = hv_bcS 

       ! ---------------- North boundary conditions ----------------------------

       READ(input_unit,north_boundary_conditions,IOSTAT=ios)

       IF ( ios .NE. 0 ) THEN

          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist NORTH_BOUNDARY_CONDITIONS'
          WRITE(*,*) 'Please check the input file'
          STOP

       ELSE

          REWIND(input_unit)

       END IF


       IF ( ( h_bcN%flag .EQ. -1 ) ) THEN 

          WRITE(*,*) 'ERROR: problem with namelist NORTH_BOUNDARY_CONDITIONS'
          WRITE(*,*) 'B.C. for h not set properly'
          WRITE(*,*) 'Please check the input file'
          STOP

       END IF


       IF ( comp_cells_x .GT. 1 ) THEN

          IF ( hu_bcN%flag .EQ. -1 ) THEN 

             WRITE(*,*) 'ERROR: problem with namelist NORTH_BOUNDARY_CONDITIONS'
             WRITE(*,*) 'B.C. for hu not set properly'             
             WRITE(*,*) 'Please check the input file'
             STOP

          END IF

       ELSE

          hu_bcN%flag = 1
          hu_bcN%value = 0.0_wp

       END IF

       IF ( comp_cells_y .GT. 1 ) THEN

          IF ( hv_bcN%flag .EQ. -1 ) THEN 

             WRITE(*,*) 'ERROR: problem with namelist NORTH_BOUNDARY_CONDITIONS'
             WRITE(*,*) 'B.C. for hv not set properly'             
             WRITE(*,*) 'Please check the input file'
             STOP

          END IF

       ELSE

          hv_bcN%flag = 1
          hv_bcN%value = 0.0_wp

       END IF

       IF ( alpha_flag ) THEN

          IF ( ANY(alphas_bcN(1:n_solid)%flag .EQ. -1 ) ) THEN 

             WRITE(*,*) 'ERROR: problem with namelist NORTH_BOUNDARY_CONDITIONS'
             WRITE(*,*) 'B.C. for sediment concentrations not set properly'
             WRITE(*,*) 'Please check the input file'
             WRITE(*,*) 'alphas_bcN'
             WRITE(*,*) alphas_bcN(1:n_solid)
             STOP

          END IF

          IF ( ANY(alphag_bcN(1:n_add_gas)%flag .EQ. -1 ) ) THEN 

             WRITE(*,*) 'ERROR: problem with namelist NORTH_BOUNDARY_CONDITIONS'
             WRITE(*,*) 'B.C. for additional gas components not set properly'
             WRITE(*,*) 'Please check the input file'
             WRITE(*,*) 'alphag_bcN'
             WRITE(*,*) alphag_bcN(1:n_add_gas)
             STOP

          END IF

          IF ( gas_flag .AND. liquid_flag ) THEN

             IF ( alphal_bcN%flag .EQ. -1 ) THEN 

                WRITE(*,*) 'ERROR: problem with namelist NORTH_BOUNDARY_CONDITIONS'
                WRITE(*,*) 'B.C. for additional gas components not set properly'
                WRITE(*,*) 'Please check the input file'
                WRITE(*,*) 'alphal_bcN'
                WRITE(*,*) alphal_bcN
                STOP

             END IF

          END IF

       ELSE

          IF ( ANY(halphas_bcN(1:n_solid)%flag .EQ. -1 ) ) THEN 

             WRITE(*,*) 'ERROR: problem with namelist NORTH_BOUNDARY_CONDITIONS'
             WRITE(*,*) 'B.C. for sediment concentrations not set properly'
             WRITE(*,*) 'Please check the input file'
             WRITE(*,*) 'halphas_bcN'
             WRITE(*,*) halphas_bcN(1:n_solid)
             STOP

          END IF

          IF ( ANY(halphag_bcN(1:n_add_gas)%flag .EQ. -1 ) ) THEN 

             WRITE(*,*) 'ERROR: problem with namelist NORTH_BOUNDARY_CONDITIONS'
             WRITE(*,*) 'B.C. for additional gas components not set properly'
             WRITE(*,*) 'Please check the input file'
             WRITE(*,*) 'halphag_bcN'
             WRITE(*,*) halphag_bcN(1:n_add_gas)
             STOP

          END IF

          IF ( gas_flag .AND. liquid_flag ) THEN

             IF ( halphal_bcN%flag .EQ. -1 ) THEN 

                WRITE(*,*) 'ERROR: problem with namelist NORTH_BOUNDARY_CONDITIONS'
                WRITE(*,*) 'B.C. for additional gas components not set properly'
                WRITE(*,*) 'Please check the input file'
                WRITE(*,*) 'halphal_bcN'
                WRITE(*,*) halphal_bcN
                STOP

             END IF

          END IF

       END IF

       IF ( T_bcN%flag .EQ. -1 ) THEN 

          WRITE(*,*) 'ERROR: problem with namelist NORTH_BOUNDARY_CONDITIONS'
          WRITE(*,*) 'B.C. for temperature not set properly'
          WRITE(*,*) 'Please check the input file'
          STOP

       END IF

       bcN(1) = h_bcN 
       bcN(2) = hu_bcN 
       bcN(3) = hv_bcN 

    END IF

    bcW(4) = T_bcW
    bcE(4) = T_bcE
    bcS(4) = T_bcS
    bcN(4) = T_bcN

    IF ( alpha_flag ) THEN

       bcW(5:4+n_solid) = alphas_bcW(1:n_solid)
       bcE(5:4+n_solid) = alphas_bcE(1:n_solid)
       bcS(5:4+n_solid) = alphas_bcS(1:n_solid)
       bcN(5:4+n_solid) = alphas_bcN(1:n_solid)

       bcW(4+n_solid+1:4+n_solid+n_add_gas) = alphag_bcW(1:n_add_gas)
       bcE(4+n_solid+1:4+n_solid+n_add_gas) = alphag_bcE(1:n_add_gas)
       bcS(4+n_solid+1:4+n_solid+n_add_gas) = alphag_bcS(1:n_add_gas)
       bcN(4+n_solid+1:4+n_solid+n_add_gas) = alphag_bcN(1:n_add_gas)

    ELSE

       bcW(5:4+n_solid) = halphas_bcW(1:n_solid)
       bcE(5:4+n_solid) = halphas_bcE(1:n_solid)
       bcS(5:4+n_solid) = halphas_bcS(1:n_solid)
       bcN(5:4+n_solid) = halphas_bcN(1:n_solid)

       bcW(4+n_solid+1:4+n_solid+n_add_gas) = halphag_bcW(1:n_add_gas)
       bcE(4+n_solid+1:4+n_solid+n_add_gas) = halphag_bcE(1:n_add_gas)
       bcS(4+n_solid+1:4+n_solid+n_add_gas) = halphag_bcS(1:n_add_gas)
       bcN(4+n_solid+1:4+n_solid+n_add_gas) = halphag_bcN(1:n_add_gas)

    END IF

    IF ( gas_flag .AND. liquid_flag ) THEN

       IF ( alpha_flag ) THEN

          bcW(n_vars) = alphal_bcW
          bcE(n_vars) = alphal_bcE
          bcS(n_vars) = alphal_bcS
          bcN(n_vars) = alphal_bcN

       ELSE

          bcW(n_vars) = halphal_bcW
          bcE(n_vars) = halphal_bcE
          bcS(n_vars) = halphal_bcS
          bcN(n_vars) = halphal_bcN

       END IF

    END IF

    ! ------- READ expl_terms_parameters NAMELIST -------------------------------

    READ(input_unit, expl_terms_parameters,IOSTAT=ios)

    IF ( ios .NE. 0 ) THEN

       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist EXPL_TERMS_PARAMETERS'
       WRITE(*,*) 'Please check the input file'
       STOP

    ELSE

       REWIND(input_unit)

    END IF

    IF ( grav .EQ. -1.0_wp ) THEN

       WRITE(*,*) 'ERROR: problem with namelist EXPL_TERMS_PARAMETERS'
       WRITE(*,*) 'GRAV not set properly'
       WRITE(*,*) 'Please check the input file'
       STOP

    ELSE

       inv_grav = 1.0_wp / grav

    END IF


    ! ------- READ collapsing_volume_parameters NAMELIST ------------------------

    IF ( collapsing_volume_flag ) THEN

       READ(input_unit,collapsing_volume_parameters,IOSTAT=ios)

       IF ( ios .NE. 0 ) THEN

          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist COLLAPSING_VOLUME_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          STOP

       ELSE

          REWIND(input_unit)

          IF ( t_collapse .EQ. -1.0_wp ) THEN

             WRITE(*,*) 'ERROR: problem with namelist                           &
                  &COLLAPSING_VOLUME_PARAMETERS'
             WRITE(*,*) 'PLEASE CHECK VALUE OF T_COLLAPSE',t_collapse
             STOP

          END IF

          IF ( h_collapse .EQ. -1.0_wp ) THEN

             WRITE(*,*) 'ERROR: problem with namelist                           &
                  &COLLAPSING_VOLUME_PARAMETERS'
             WRITE(*,*) 'PLEASE CHECK VALUE OF H_COLLAPSE',h_collapse
             STOP

          END IF

          IF ( r_collapse .EQ. -1.0_wp ) THEN

             WRITE(*,*) 'ERROR: problem with namelist                           &
                  &COLLAPSING_VOLUME_PARAMETERS'
             WRITE(*,*) 'PLEASE CHECK VALUE OF R_COLLAPSE',r_collapse
             STOP

          END IF

          IF ( ( x_collapse - r_collapse ) .LE. X0 + cell_size ) THEN

             WRITE(*,*) 'ERROR: problem with namelist                           &
                  &COLLAPSING_VOLUME_PARAMETERS'
             WRITE(*,*) 'COLLAPSING VOLUME TOO LARGE'
             WRITE(*,*) ' x_collapse - radius ',x_collapse-r_collapse
             STOP

          END IF

          IF ( (x_collapse+r_collapse) .GE. X0+(comp_cells_x-1)*cell_size ) THEN

             WRITE(*,*) 'ERROR: problem with namelist                           &
                  &COLLAPSING_VOLUME_PARAMETERS'
             WRITE(*,*) 'COLLAPSING VOLUME TOO LARGE'
             WRITE(*,*) ' x_collapse + radius ',x_collapse+r_collapse
             STOP

          END IF

          IF ( ( y_collapse - r_collapse ) .LE. Y0 + cell_size ) THEN

             WRITE(*,*) 'ERROR: problem with namelist                           &
                  &COLLAPSING_VOLUME_PARAMETERS'
             WRITE(*,*) 'COLLAPSING VOLUME TOO LARGE'
             WRITE(*,*) ' y_collapse - radius ',y_collapse-r_collapse
             STOP

          END IF

          IF ( (y_collapse+r_collapse) .GE. Y0+(comp_cells_y-1)*cell_size ) THEN

             WRITE(*,*) 'ERROR: problem with namelist                           &
                  &COLLAPSING_VOLUME_PARAMETERS'
             WRITE(*,*) 'COLLAPSING VOLUME TOO LARGE'
             WRITE(*,*) ' y_collapse + radius ',y_collapse+r_collapse
             STOP

          END IF

          IF ( ANY(alphas_collapse(1:n_solid) .EQ. -1.0_wp ) ) THEN

             WRITE(*,*) 'ERROR: problem with namelist                           &
                  &COLLAPSING_VOLUME_PARAMETERS'
             WRITE(*,*) 'alphas_collpase =' , alphas_collapse(1:n_solid)
             WRITE(*,*) 'Please check the input file'
             STOP

          END IF

          IF ( ANY(alphag_collapse(1:n_add_gas) .EQ. -1.0_wp ) ) THEN

             WRITE(*,*) 'ERROR: problem with namelist                           &
                  &COLLAPSING_VOLUME_PARAMETERS'
             WRITE(*,*) 'alphag_collpase =' , alphag_collapse(1:n_add_gas)
             WRITE(*,*) 'Please check the input file'
             STOP

          END IF

       END IF

    END IF

    ! ------- READ rheology_parameters NAMELIST ---------------------------------

    IF ( rheology_flag ) THEN

       READ(input_unit, rheology_parameters,IOSTAT=ios)

       IF ( ios .NE. 0 ) THEN

          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          STOP

       ELSE

          REWIND(input_unit)

       END IF

       IF ( rheology_model .EQ. 0 ) THEN

          WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
          WRITE(*,*) 'RHEOLOGY_FLAG' , rheology_flag , 'RHEOLOGY_MODEL =' ,     &
               rheology_model
          WRITE(*,*) 'Please check the input file'
          STOP

       ELSEIF ( rheology_model .EQ. 1 ) THEN

          IF ( ( mu .EQ. -1.0_wp ) .AND. ( xi .EQ. -1.0_wp ) ) THEN

             WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
             WRITE(*,*) 'RHEOLOGY_MODEL =' , rheology_model
             WRITE(*,*) 'MU =' , mu ,' XI =' , xi
             WRITE(*,*) 'Please check the input file'
             STOP

          END IF

          IF ( ( T_ref .NE. -1.0_wp ) .OR. ( nu_ref .NE. -1.0_wp ) .OR.         &
               ( visc_par .NE. -1.0_wp ) .OR. ( tau .NE. -1.0_wp ) .OR.         &
               ( tau0 .NE. -1.0_wp ) ) THEN

             WRITE(*,*) 'WARNING: parameters not used in RHEOLOGY_PARAMETERS'
             IF ( T_ref .NE. -1.0_wp ) WRITE(*,*) 'T_ref =',T_ref 
             IF ( nu_ref .NE. -1.0_wp ) WRITE(*,*) 'nu_ref =',nu_ref 
             IF ( visc_par .NE. -1.0_wp ) WRITE(*,*) 'visc_par =',visc_par
             IF ( tau .NE. -1.0_wp ) WRITE(*,*) 'tau =',tau 
             IF ( tau0 .NE. -1.0_wp ) WRITE(*,*) 'tau0 =',tau0 
             WRITE(*,*) 'Press ENTER to continue'
             READ(*,*)

          END IF

       ELSEIF ( rheology_model .EQ. 2 ) THEN

          IF ( tau .EQ. -1.0_wp )  THEN

             WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
             WRITE(*,*) 'RHEOLOGY_MODEL =' , rheology_model
             WRITE(*,*) 'TAU =' , tau
             WRITE(*,*) 'Please check the input file'
             STOP

          END IF

          IF ( ( T_ref .NE. -1.0_wp ) .OR. ( nu_ref .NE. -1.0_wp ) .OR.         &
               ( visc_par .NE. -1.0_wp ) .OR. ( mu .NE. -1.0_wp ) .OR.          &
               ( xi .NE. -1.0_wp ) ) THEN

             WRITE(*,*) 'WARNING: parameters not used in RHEOLOGY_PARAMETERS'
             IF ( T_ref .NE. -1.0_wp ) WRITE(*,*) 'T_ref =',T_ref 
             IF ( nu_ref .NE. -1.0_wp ) WRITE(*,*) 'nu_ref =',nu_ref 
             IF ( visc_par .NE. -1.0_wp ) WRITE(*,*) 'visc_par =',visc_par
             IF ( mu .NE. -1.0_wp ) WRITE(*,*) 'mu =',mu 
             IF ( xi .NE. -1.0_wp ) WRITE(*,*) 'xi =',xi
             WRITE(*,*) 'Press ENTER to continue'
             READ(*,*)


          END IF

       ELSEIF ( rheology_model .EQ. 3 ) THEN

          IF ( nu_ref .EQ. -1.0_wp ) THEN

             WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
             WRITE(*,*) 'NU_REF =' , nu_ref 
             WRITE(*,*) 'Please check the input file'
             STOP

          END IF

          IF ( tau0 .EQ. -1.0_wp ) THEN

             WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
             WRITE(*,*) 'TAU0 =' , tau0 
             WRITE(*,*) 'Please check the input file'
             STOP

          END IF

          IF ( visc_par .EQ. -1.0_wp ) THEN

             WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
             WRITE(*,*) 'VISC_PAR =' , visc_par
             WRITE(*,*) 'Please check the input file'
             STOP

          ELSEIF ( visc_par .EQ. 0.0_wp ) THEN

             WRITE(*,*) 'WARNING: temperature and momentum uncoupled'
             WRITE(*,*) 'VISC_PAR =' , visc_par
             WRITE(*,*) 'Press ENTER to continue'
             READ(*,*)

          ELSE

             IF ( T_ref .EQ. -1.0_wp ) THEN

                WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
                WRITE(*,*) 'T_REF =' , T_ref 
                WRITE(*,*) 'Please check the input file'
                STOP

             END IF

          END IF

          IF ( ( mu .NE. -1.0_wp ) .OR. ( xi .NE. -1.0_wp ) .OR.                &
               ( tau .NE. -1.0_wp ) ) THEN

             WRITE(*,*) 'WARNING: parameters not used in RHEOLOGY_PARAMETERS'
             IF ( mu .NE. -1.0_wp ) WRITE(*,*) 'mu =',mu 
             IF ( xi .NE. -1.0_wp ) WRITE(*,*) 'xi =',xi
             IF ( tau .NE. -1.0_wp ) WRITE(*,*) 'tau =',tau 
             WRITE(*,*) 'Press ENTER to continue'
             READ(*,*)

          END IF

       ELSEIF ( rheology_model .EQ. 4 ) THEN

          IF ( gas_flag .OR. ( .NOT. liquid_flag ) ) THEN

             WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
             WRITE(*,*) 'RHEOLOGY_MODEL =' , rheology_model
             WRITE(*,*) 'GAS FLAG = ' , gas_flag
             WRITE(*,*) 'LIQUID FLAG = ' , liquid_flag
             STOP

          END IF

          IF ( restart .AND. ( ANY(sed_vol_perc(1:n_solid) .EQ. -1.0_wp ) ) ) THEN

             WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
             WRITE(*,*) 'RHEOLOGY_MODEL =' , rheology_model
             WRITE(*,*) 'SED_VOL_PERC = ' , sed_vol_perc(1:n_solid)
             STOP

          END IF

          IF ( alpha2 .EQ. -1.0_wp ) THEN

             WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
             WRITE(*,*) 'ALPHA2 =' , alpha2 
             WRITE(*,*) 'Please check the input file'
             STOP

          END IF

          IF ( beta2 .EQ. -1.0_wp ) THEN

             WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
             WRITE(*,*) 'BETA2 =' , beta2 
             WRITE(*,*) 'Please check the input file'
             STOP

          END IF

          IF ( T_ref .LE. 273.15_wp ) THEN

             WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
             WRITE(*,*) 'T_REF =' , T_ref
             WRITE(*,*) 'Please check the input file'
             STOP

          END IF

          IF ( alpha1_ref .EQ. -1.0_wp ) THEN

             WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
             WRITE(*,*) 'ALPHA1 =' , alpha1_ref 
             WRITE(*,*) 'Please check the input file'
             STOP

          ELSE

             Tc = T_ref - 273.15_wp

             IF ( Tc .LT. 20.0_wp ) THEN

                expA = 1301.0_wp / ( 998.333_wp + 8.1855_wp * ( Tc - 20.0_wp )  &
                     + 0.00585_wp * ( Tc - 20.0_wp )**2 ) - 1.30223_wp

                alpha1_coeff = alpha1_ref / ( 1.0E-3_wp * 10.0_wp**expA )

             ELSE

                expB = ( 1.3272_wp * ( 20.0_wp - Tc ) - 0.001053_wp *           &
                     ( Tc - 20.0_wp )**2 ) / ( Tc + 105.0_wp )

                alpha1_coeff = alpha1_ref / ( 1.002E-3_wp * 10.0_wp**expB )

             END IF

             WRITE(*,*) 'alpha1 coefficient:',alpha1_coeff

          END IF

          IF ( beta1 .EQ. -1.0_wp ) THEN

             WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
             WRITE(*,*) 'BETA1 =' , beta1 
             WRITE(*,*) 'Please check the input file'
             STOP

          END IF

          IF ( Kappa .EQ. -1.0_wp ) THEN

             WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
             WRITE(*,*) 'KAPPA =' , kappa 
             WRITE(*,*) 'Please check the input file'
             STOP

          END IF

          IF ( n_td .EQ. -1.0_wp ) THEN

             WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
             WRITE(*,*) 'N_TD =' , n_td 
             WRITE(*,*) 'Please check the input file'
             STOP

          ELSE

             n_td2 = n_td**2

          END IF

       ELSEIF ( rheology_model .EQ. 5 ) THEN

          IF ( VERBOSE_LEVEL .GE. 0 ) THEN

             WRITE(*,*) 'RHEOLOGY_MODEL =' , rheology_model
             WRITE(*,*) 'Kurganov & Petrova Example 5'

          END IF

       ELSEIF ( ( rheology_model .EQ. 6 ) .OR. ( rheology_model .EQ. 8 ) ) THEN

          IF ( VERBOSE_LEVEL .GE. 0 ) THEN

             WRITE(*,*) 'RHEOLOGY_MODEL =' , rheology_model
             WRITE(*,*) 'Bursik & Woods'

          END IF

          IF ( friction_factor .LT. 0.0_wp ) THEN

             WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
             WRITE(*,*) 'FRICTION_FACTOR =' , friction_factor 
             WRITE(*,*) 'Please check the input file'
             STOP

          END IF

       ELSEIF ( rheology_model .EQ. 9 ) THEN
          WRITE(*,*) 'RHEOLOGY_MODEL =' , rheology_model
          WRITE(*,*) 'MU_0 =' , mu_0 ,' MU_inf =' , mu_inf, 'Fr_0 = ', Fr_0
          IF ( ( mu_0 .EQ. -1.0_wp ) .AND. ( mu_inf .EQ. -1.0_wp ) .AND. ( Fr_0 .EQ. -1.0_wp ) ) THEN
             ! Stop the program if the parameters are not a double-precision floating-point numbers (real numbers).
             WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
             WRITE(*,*) 'RHEOLOGY_MODEL =' , rheology_model
             WRITE(*,*) 'MU_0 =' , mu_0 ,' MU_inf =' , mu_inf, 'Fr_0 = ', Fr_0
             WRITE(*,*) 'Please check the input file'
             STOP
          END IF

          !! Should add something similar? This is present in the Voellmy-Salm Rheology
          ! IF ( ( T_ref .NE. -1.0_wp ) .OR. ( nu_ref .NE. -1.0_wp ) .OR.         &
          !      ( visc_par .NE. -1.0_wp ) .OR. ( tau .NE. -1.0_wp ) .OR.         &
          !      ( tau0 .NE. -1.0_wp ) ) THEN

          !    WRITE(*,*) 'WARNING: parameters not used in RHEOLOGY_PARAMETERS'
          !    IF ( T_ref .NE. -1.0_wp ) WRITE(*,*) 'T_ref =',T_ref 
          !    IF ( nu_ref .NE. -1.0_wp ) WRITE(*,*) 'nu_ref =',nu_ref 
          !    IF ( visc_par .NE. -1.0_wp ) WRITE(*,*) 'visc_par =',visc_par
          !    IF ( tau .NE. -1.0_wp ) WRITE(*,*) 'tau =',tau 
          !    IF ( tau0 .NE. -1.0_wp ) WRITE(*,*) 'tau0 =',tau0 
          !    WRITE(*,*) 'Press ENTER to continue'
          !    READ(*,*)
          ! END IF

       ELSE

          WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
          WRITE(*,*) 'RHEOLOGY_MODEL =' , rheology_model
          WRITE(*,*) 'Please check the input file'
          STOP

       END IF


    END IF

    ! ------- READ temperature_parameters NAMELIST ------------------------------

    READ(input_unit, temperature_parameters,IOSTAT=ios)

    IF ( ios .NE. 0 ) THEN

       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist TEMPERATURE_PARAMETERS'
       WRITE(*,*) 'Please check the input file'
       STOP

    ELSE

       REWIND(input_unit)

       IF ( rheology_model .EQ. 3 ) THEN

          radiative_term_coeff = SBconst * emissivity * exp_area_fract
          WRITE(*,*) 'radiative_term_coeff',radiative_term_coeff

          convective_term_coeff = atm_heat_transf_coeff * exp_area_fract

       END IF

    END IF

    ! ---------------------------------------------------------------------------

    IF ( verbose_level .GE. 1 ) WRITE(*,*) 

    IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'Searching for DEM file'

    INQUIRE(FILE=topography_file,EXIST=lexist)

    IF(lexist)THEN

       OPEN(2001, file=topography_file, status='old', action='read')

    ELSE

       WRITE(*,*) 'no dem file: ',TRIM(topography_file)
       STOP

    ENDIF

    READ(2001,*) chara, ncols
    READ(2001,*) chara, nrows
    READ(2001,*) chara, xllcorner
    READ(2001,*) chara, yllcorner
    READ(2001,*) chara, cellsize
    READ(2001,*) chara, nodata_topo

    ! The values read from the DEM files are associated to the center of the
    ! pixels. x0 is the left margin of the computational domain and has to be
    ! greater than the center of the first pixel.
    IF ( x0 - ( xllcorner + 0.5_wp * cellsize ) .LT. -1.E-10_wp  ) THEN 

       WRITE(*,*) 'Computational domain problem'
       WRITE(*,*) 'x0 < xllcorner+0.5*cellsize',x0,xllcorner+0.5_wp*cellsize
       STOP

    END IF

    ! The right margin of the computational domain should be smaller then the
    ! center of the last pixel
    IF ( x0 + ( comp_cells_x ) * cell_size .GT.                                 &
         xllcorner + ( 0.5_wp + ncols ) * cellsize ) THEN 

       WRITE(*,*) 'Computational domain problem'
       WRITE(*,*) 'right edge > xllcorner+ncols*cellsize',                      &
            x0+comp_cells_x*cell_size , xllcorner+(0.5_wp+ncols)*cellsize
       STOP

    END IF

    IF ( y0 - ( yllcorner+0.5_wp*cellsize ) .LT. -1.E-10_wp ) THEN 

       WRITE(*,*) 'Computational domain problem'
       WRITE(*,*) 'y0 < yllcorner+0.5*cellsize',y0,yllcorner+0.5_wp*cellsize
       STOP

    END IF

    IF ( ABS( ( y0 + comp_cells_y * cell_size ) - ( yllcorner + 0.5_wp +        &
         nrows * cellsize ) ) .LT. 1.E-10_wp ) THEN 

       WRITE(*,*) 'Computational domain problem'
       WRITE(*,*) 'top edge > yllcorner+nrows*cellsize',                        &
            y0+comp_cells_y*cell_size , yllcorner+(0.5_wp+nrows)*cellsize
       STOP

    END IF

    IF ( VERBOSE_LEVEL .GE. 0 ) THEN

       WRITE(*,*) 'Reading DEM file' 
       WRITE(*,*) 'ncols',ncols
       WRITE(*,*) 'nrows',nrows

    END IF

    n_topography_profile_x = ncols

    n_topography_profile_y = nrows

    ALLOCATE( topography_profile( 3 , n_topography_profile_x ,                  &
         n_topography_profile_y) )

    DO j=1,n_topography_profile_x 

       topography_profile(1,j,:) = xllcorner + ( j - 0.5_wp ) * cellsize

    ENDDO

    DO k=1,n_topography_profile_y

       topography_profile(2,:,k) = yllcorner + ( k - 0.5_wp ) * cellsize

    ENDDO

    ! Read topography values (starting at the upper-left corner)

    DO k=1,n_topography_profile_y

       IF ( VERBOSE_LEVEL .GE. 0 ) THEN

          WRITE(*,FMT="(A1,A,t21,F6.2,A)",ADVANCE="NO") ACHAR(13),              &
               & " Percent Complete: " ,                                        &
               ( REAL(k) / REAL(n_topography_profile_y))*100.0, "%"

       END IF

       READ(2001,*) topography_profile(3,:,n_topography_profile_y-k+1)

    ENDDO


    IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) ''

    CLOSE(2001)

    ! ----------- READ vulnerability_table_parameters NAMELIST ------------------

    READ(input_unit, vulnerability_table_parameters,IOSTAT=ios)

    IF ( ios .NE. 0 ) THEN

       IF ( verbose_level .GE. 0.0_wp ) THEN

          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'WARNING: namelist VULNERABILITY_TABLE_PARAMETERS not found'
          REWIND(input_unit)
          
       END IF

    ELSE

       REWIND(input_unit)

       n_thickness_levels = COUNT( thickness_levels0 .GE. 0.0_wp )
       n_dyn_pres_levels = COUNT( dyn_pres_levels0 .GE. 0.0_wp )

       IF ( n_thickness_levels .GT. 0 ) THEN

          ALLOCATE( thickness_levels(n_thickness_levels) )
          thickness_levels(1:n_thickness_levels) =                              &
               thickness_levels0(1:n_thickness_levels)
          IF ( ANY(thickness_levels .LT. 0.0_wp ) ) THEN

             WRITE(*,*) 'ERROR: problem with namelist ',                        &
                  'VULNERABILITY_TABLE_PARAMETERS'
             WRITE(*,*) 'Please check the input file'
             WRITE(*,*) 'thickness_levels(1:n_thickness_levels)',               &
                  thickness_levels(1:n_thickness_levels)

             STOP

          END IF

          IF ( n_dyn_pres_levels .LE. 0 ) THEN

             n_dyn_pres_levels = 1
             dyn_pres_levels0(1:n_dyn_pres_levels) = 0.0_wp

          END IF

       END IF

       IF ( n_dyn_pres_levels .GT. 0 ) THEN

          ALLOCATE( dyn_pres_levels(n_dyn_pres_levels) )
          dyn_pres_levels(1:n_dyn_pres_levels) =                                &
               dyn_pres_levels0(1:n_dyn_pres_levels)
          IF ( ANY(dyn_pres_levels .LT. 0.0_wp ) ) THEN

             WRITE(*,*) 'ERROR: problem with namelist ',                        &
                  'VULNERABILITY_TABLE_PARAMETERS'
             WRITE(*,*) 'Please check the input file'
             WRITE(*,*) 'dyn_pres_levels(1:n_thickness_levels)',                &
                  dyn_pres_levels(1:n_dyn_pres_levels)

             STOP

          END IF

          IF ( n_thickness_levels .LE. 0 ) THEN

             n_thickness_levels = 1
             ALLOCATE( thickness_levels(n_thickness_levels) )
             thickness_levels(1:n_thickness_levels) = 0.0_wp

          END IF

       END IF

    END IF

    IF ( verbose_level .GE. 0 ) THEN

       WRITE(*,*) 'thickness_levels',n_thickness_levels,thickness_levels
       WRITE(*,*) 'dyn_pres_levels',n_dyn_pres_levels,dyn_pres_levels

    END IF

    ! ------------ READ vertical_profiles_parameters NAMELIST ------------------

    IF ( vertical_profiles_flag ) THEN

       IF ( rheology_flag .AND. ( rheology_model .NE. 8 ) ) THEN

          WRITE(*,*) 'ERROR: problem with input values'
          WRITE(*,*) 'RHEOLOGY_model',rheology_model
          WRITE(*,*) 'VERTICAL_PROFILES_FLAG',vertical_profiles_flag
          WRITE(*,*) 'Vertical profiles can be used only with RHEOLOGY_FLAG = F'
          WRITE(*,*) ' or with RHEOLOGY_MODEL = 8'
          STOP

       END IF

       IF ( energy_flag ) THEN

          WRITE(*,*) 'ERROR: problem with input values'
          WRITE(*,*) 'ENERGY_flag',energy_flag
          WRITE(*,*) 'VERTICAL_PROFILES_FLAG',vertical_profiles_flag
          WRITE(*,*) 'Vertical profiles can be used only with ENERGY_FLAG = F'
          STOP

       END IF

       REWIND(input_unit)    
       READ(input_unit, vertical_profiles_parameters,IOSTAT=ios)

       IF ( ios .NE. 0 ) THEN

          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist VERTICAL_PROFILES_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          WRITE(*,vertical_profiles_parameters) 
          STOP

       ELSE

          REWIND(input_unit)

          IF ( vonK .LE. 0.0_wp ) THEN

             WRITE(*,*) 'ERROR: problem with namelist VERTICAL_PROFILES_PARAMETERS'
             WRITE(*,*) 'vonK =' , vonK 
             WRITE(*,*) 'Please check the input file'
             STOP

          END IF

          IF ( k_s .LE. 0.0_wp ) THEN

             WRITE(*,*) 'ERROR: problem with namelist VERTICAL_PROFILES_PARAMETERS'
             WRITE(*,*) 'k_s =' , k_s 
             WRITE(*,*) 'Please check the input file'
             STOP

          END IF

          IF ( N_quad .LE. 0 ) THEN

             WRITE(*,*) 'ERROR: problem with namelist VERTICAL_PROFILES_PARAMETERS'
             WRITE(*,*) 'N_QUAD =' , n_quad
             WRITE(*,*) 'Please check the input file'
             STOP

          ELSE

             ALLOCATE( x_quad(N_quad) , w_quad(N_quad) )
             CALL gaulegf(-1.0_wp, 1.0_wp, x_quad, w_quad, n_quad)

          END IF

          IF ( Sc .LE. 0.0_wp ) THEN

             WRITE(*,*) 'ERROR: problem with namelist VERTICAL_PROFILES_PARAMETERS'
             WRITE(*,*) 'Sc =' , Sc 
             WRITE(*,*) 'Please check the input file'
             STOP

          END IF

          a_crit_rel = vonK / SQRT( friction_factor) + 1.0_wp
          H_crit_rel = 1.0_wp / 30.0_wp * ( -a_crit_rel /                       &
               lambertw0(-a_crit_rel *  EXP(-a_crit_rel) ) - 1.0_wp )

          !WRITE(*,*) 'a_crit_rel',a_crit_rel 
          !WRITE(*,*) 'Arg lambertW',-a_crit_rel *  EXP(-a_crit_rel)
          !WRITE(*,*) 'LambertW0', lambertw0(-a_crit_rel *  EXP(-a_crit_rel) )
          WRITE(*,*) 'H_crit',H_crit_rel*k_s

          !H_crit_rel = 1.0_wp / 30.0_wp * ( -a_crit_rel /                       &
          !     lambertwm1(-a_crit_rel *  EXP(-a_crit_rel) ) - 1.0_wp )

          !WRITE(*,*) 'LambertWm1', lambertwm1(-a_crit_rel *  EXP(-a_crit_rel) )
          !WRITE(*,*) 'H_crit',H_crit_rel*k_s
          !READ(*,*) 


       END IF

    END IF

    ! ------- READ radial_source_parameters NAMELIST ----------------------------

    source_flag:IF ( ( radial_source_flag ) .OR. ( bottom_radial_source_flag )  &
         .OR. ( lateral_source_flag ) ) THEN

       alphal_source = -1.0_wp

       IF ( lateral_source_flag ) THEN

          WRITE(*,*) 'Searching for namelist LATERAL_SOURCE_PARAMETERS'
          READ(input_unit,lateral_source_parameters,IOSTAT=ios)
          source_name = 'LATERAL_SOURCE_PARAMETERS'
          
       ELSE
          
          READ(input_unit,radial_source_parameters,IOSTAT=ios)
          source_name = 'RADIAL_SOURCE_PARAMETERS'

       END IF
          
       IF ( ios .NE. 0 ) THEN

          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist ',TRIM(source_name)

          IF ( lateral_source_flag ) THEN

             WRITE(*,lateral_source_parameters)

          ELSE

             WRITE(*,radial_source_parameters)

          END IF
          WRITE(*,*) 'Please check the input file'
          STOP

       ELSE

          REWIND(input_unit)

          IF ( t_source .EQ. -1.0_wp ) THEN

             WRITE(*,*) 'ERROR: problem with namelist ',TRIM(source_name)
             WRITE(*,*) 'PLEASE CHECK VALUE OF T_SOURCE',t_source
             STOP

          END IF

          IF ( ( ( h_source .EQ. -1.0_wp ) .AND. ( mfr_source .EQ. -1 ) ) &
               .AND. (.NOT. bottom_radial_source_flag) ) THEN

             WRITE(*,*) 'ERROR: problem with namelist ',TRIM(source_name)
             WRITE(*,*) 'PLEASE ASSIGN A VALUE TO H_SOURCE OR MFR_SOURCE'
             STOP

          ELSE

             IF ( ( ( h_source .GE. 0.0_wp ) .OR. ( mfr_source .GE. 0.0_wp ) ) &
                  .AND. ( bottom_radial_source_flag ) ) THEN

                WRITE(*,*) 'ERROR: problem with namelist RADIAL_SOURCE_PARAMETERS'
                WRITE(*,*) 'When BOTTOM_RADIAL_SOURCE_FLAG = TRUE'
                WRITE(*,*) 'h_source should not be given',h_source
                WRITE(*,*) 'mfr_source should not be given',mfr_source
                STOP

             END IF

          END IF

          IF ( ( h_source .GE. 0.0_wp ) .AND. ( mfr_source .GE. 0.0_wp ) ) THEN

             WRITE(*,*) 'ERROR: problem with namelist ',TRIM(source_name)
             WRITE(*,*) 'PLEASE ASSIGN ONLY H_SOURCE OR MFR_SOURCE'
             WRITE(*,*) 'h_source',h_source
             WRITE(*,*) 'mfr_source',mfr_source
             STOP             

          END IF

          IF ( ( radial_source_flag ) .AND. ( r_source .EQ. -1.0_wp ) ) THEN

             WRITE(*,*) 'ERROR: problem with namelist ',TRIM(source_name)
             WRITE(*,*) 'PLEASE CHECK VALUE OF R_SOURCE',r_source
             STOP

          END IF

          IF ( ( vel_source .EQ. -1.0_wp ) .AND. ( bottom_radial_source_flag ) )&
               THEN

             WRITE(*,*) 'ERROR: problem with namelist RADIAL_SOURCE_PARAMETERS'
             WRITE(*,*) 'PLEASE CHECK VALUE OF VEL_SOURCE',vel_source
             STOP

          END IF

          IF ( ( Ri_source .NE. -1.0_wp ) .AND. ( bottom_radial_source_flag ) ) &
               THEN

             WRITE(*,*) 'ERROR: problem with namelist RADIAL_SOURCE_PARAMETERS'
             WRITE(*,*) 'PLEASE REMOVE VALUE OF RI_SOURCE',Ri_source
             STOP

          END IF

          IF ( ( lateral_source_flag ) .OR. ( radial_source_flag ) ) THEN

             IF ( ( vel_source .EQ. -1.0_wp ) .AND. ( Ri_source .EQ. -1.0_wp ) )&
                  THEN

                WRITE(*,*) 'ERROR: problem with namelist ',TRIM(source_name)
                WRITE(*,*) 'PLEASE ASSIGN VALUE TO VEL_SOURCE OR RI_SOURCE'
                STOP

             ELSEIF (( vel_source .NE. -1.0_wp ).AND.( Ri_source .NE. -1.0_wp ))&
                  THEN

                WRITE(*,*) 'ERROR: problem with namelist ',TRIM(source_name)
                WRITE(*,*) 'PLEASE ASSIGN ONLY VEL_SOURCE OR RI_SOURCE'
                WRITE(*,*) 'VEL_SOURCE:',vel_source
                WRITE(*,*) 'RI_SOURCE:',Ri_source
                STOP

             END IF

          END IF

          IF ( radial_source_flag ) THEN

             pi_g = ATAN(1.0_wp)*4.0_wp
             source_length = 2.0 * pi_g * r_source
          
             IF ( ( x_source - r_source ) .LE. X0 + cell_size ) THEN

                WRITE(*,*) 'ERROR: problem with namelist RADIAL_SOURCE_PARAMETERS'
                WRITE(*,*) 'SOURCE TOO LARGE'
                WRITE(*,*) ' x_source - radius ',x_source-r_source
                STOP
                
             END IF
             
             IF ( ( x_source + r_source ) .GE. X0+(comp_cells_x-1)*cell_size ) THEN
                
                WRITE(*,*) 'ERROR: problem with namelist RADIAL_SOURCE_PARAMETERS'
                WRITE(*,*) 'SOURCE TOO LARGE'
                WRITE(*,*) ' x_source + radius ',x_source+r_source
                STOP
                
             END IF
             
             IF ( ( y_source - r_source ) .LE. Y0 + cell_size ) THEN
                
                WRITE(*,*) 'ERROR: problem with namelist RADIAL_SOURCE_PARAMETERS'
                WRITE(*,*) 'SOURCE TOO LARGE'
                WRITE(*,*) ' y_source - radius ',y_source-r_source
                STOP
                
             END IF
             
             IF ( ( y_source + r_source ) .GE. Y0+(comp_cells_y-1)*cell_size ) THEN
                
                WRITE(*,*) 'ERROR: problem with namelist RADIAL_SOURCE_PARAMETERS'
                WRITE(*,*) 'SOURCE TOO LARGE'
                WRITE(*,*) ' y_source + radius ',y_source+r_source
                STOP
                
             END IF

          END IF

          IF ( lateral_source_flag ) THEN

             IF ( ( source_side .NE. 'E' ) .AND. ( source_side .NE. 'W' ) .AND. &
                  ( source_side .NE. 'S' ) .AND. ( source_side .NE. 'N' ) ) THEN
                
                WRITE(*,*) 'ERROR: problem with namelist ',TRIM(source_name)
                WRITE(*,*) 'SOURCE_SIDE SHOULD BE: E, W, N or S'
                WRITE(*,*) 'SOURCE_SIDE', source_side
                STOP
                
             END IF
                          
             IF ( ( source_side .EQ. 'E' ) .OR. ( source_side .EQ. 'W' ) ) THEN

                IF ( y1_source .GE. y2_source ) THEN
                   
                   WRITE(*,*) 'ERROR: problem with namelist ',TRIM(source_name)
                   WRITE(*,*) 'Y1_SOURCE >= Y2_SOURCE'
                   WRITE(*,*) y1_source,y2_source
                   STOP

                END IF

                IF ( y1_source .LE. Y0 ) THEN
                   
                   WRITE(*,*) 'ERROR: problem with namelist ',TRIM(source_name)
                   WRITE(*,*) 'Y1_SOURCE TOO SMALL'
                   WRITE(*,*) ' y1_source',y1_source
                   STOP

                END IF

                IF ( y2_source .GE. Y0+comp_cells_y*cell_size ) THEN
                   
                   WRITE(*,*) 'ERROR: problem with namelist ',TRIM(source_name)
                   WRITE(*,*) 'Y2_SOURCE TOO LARGE'
                   WRITE(*,*) ' y2_source',y2_source
                   STOP

                END IF

                source_length = y2_source - y1_source
                
             END IF

             IF ( ( source_side .EQ. 'S' ) .OR. ( source_side .EQ. 'N' ) ) THEN

                IF ( x1_source .GE. x2_source ) THEN
                   
                   WRITE(*,*) 'ERROR: problem with namelist ',TRIM(source_name)
                   WRITE(*,*) 'X1_SOURCE >= X2_SOURCE'
                   WRITE(*,*) x1_source,x2_source
                   STOP

                END IF

                IF ( x1_source .LE. X0 ) THEN
                   
                   WRITE(*,*) 'ERROR: problem with namelist ',TRIM(source_name)
                   WRITE(*,*) 'X1_SOURCE TOO SMALL'
                   WRITE(*,*) ' X1_source',X1_source
                   WRITE(*,*) 'X0',X0
                   STOP

                END IF

                IF ( x2_source .GE. X0+comp_cells_x*cell_size ) THEN
                   
                   WRITE(*,*) 'ERROR: problem with namelist ',TRIM(source_name)
                   WRITE(*,*) 'X2_SOURCE TOO LARGE'
                   WRITE(*,*) ' x2_source',x2_source
                   STOP

                END IF

                source_length = x2_source - x1_source
                
             END IF

          END IF


          IF ( ( ANY(alphas_source(1:n_solid) .EQ. -1.0_wp ) ) .AND.           &
               ( ANY(xs_source(1:n_solid) .EQ. -1.0_wp ) ) )  THEN
             
             WRITE(*,*) 'ERROR: problem with namelist RADIAL_VOLUME_PARAMETERS'
             WRITE(*,*) 'alphas_source =' , alphas_source(1:n_solid)
             WRITE(*,*) 'xs_source =' , xs_source(1:n_solid)
             WRITE(*,*) 'Please check the input file'
             STOP
                
          END IF

          IF ( ( ANY(alphas_source(1:n_solid) .GT. -1.0_wp ) ) .AND.           &
               ( ANY(xs_source(1:n_solid) .GT. -1.0_wp ) ) )  THEN

             WRITE(*,*) 'ERROR: problem with namelist RADIAL_VOLUME_PARAMETERS'
             WRITE(*,*) 'DEFINE ONLY solid volume fract. or solid mass fract.'
             WRITE(*,*) 'alphas_source =' , alphas_source(1:n_solid)
             WRITE(*,*) 'xs_source =' , xs_source(1:n_solid)
             WRITE(*,*) 'Please check the input file'
             STOP

          END IF

          IF ( ( ANY(alphag_source(1:n_add_gas) .EQ. -1.0_wp ) ) .AND.         &
               ( ANY(xg_source(1:n_add_gas) .EQ. -1.0_wp ) ) )  THEN

             WRITE(*,*) 'ERROR: problem with namelist RADIAL_VOLUME_PARAMETERS'
             WRITE(*,*) 'alphag_source =' , alphag_source(1:n_add_gas)
             WRITE(*,*) 'xg_source =' , xg_source(1:n_add_gas)
             WRITE(*,*) 'Please check the input file'
             STOP

          END IF

          IF ( ( ANY(alphag_source(1:n_add_gas) .GT. -1.0_wp ) ) .AND.         &
               ( ANY(xg_source(1:n_add_gas) .GT. -1.0_wp ) ) )  THEN

             WRITE(*,*) 'ERROR: problem with namelist RADIAL_VOLUME_PARAMETERS'
             WRITE(*,*) 'DEFINE ONLY add.gas volume fract. or mass fract.'
             WRITE(*,*) 'alphag_source =' , alphag_source(1:n_add_gas)
             WRITE(*,*) 'xg_source =' , xg_source(1:n_add_gas)
             WRITE(*,*) 'Please check the input file'
             STOP

          END IF


          IF ( ( ANY(alphas_source(1:n_solid) .GT. -1.0_wp ) ) .AND.         &
               ( ANY(xg_source(1:n_add_gas) .GT. -1.0_wp ) ) )  THEN

             WRITE(*,*) 'ERROR: problem with namelist RADIAL_VOLUME_PARAMETERS'
             WRITE(*,*) 'DEFINE ONLY volume fractions or mass fractions'
             WRITE(*,*) 'alphas_source =' , alphas_source(1:n_solid)
             WRITE(*,*) 'xg_source =' , xg_source(1:n_add_gas)
             WRITE(*,*) 'Please check the input file'
             STOP

          END IF

          IF ( ( ANY(xs_source(1:n_solid) .GT. -1.0_wp ) ) .AND.         &
               ( ANY(alphag_source(1:n_add_gas) .GT. -1.0_wp ) ) )  THEN

             WRITE(*,*) 'ERROR: problem with namelist RADIAL_VOLUME_PARAMETERS'
             WRITE(*,*) 'DEFINE ONLY volume fractions or mass fractions'
             WRITE(*,*) 'xs_source =' , xs_source(1:n_solid)
             WRITE(*,*) 'alphag_source =' , alphag_source(1:n_add_gas)
             WRITE(*,*) 'Please check the input file'
             STOP

          END IF


          IF ( gas_flag .AND. liquid_flag ) THEN

             IF ( ( ANY(alphag_source(1:n_add_gas) .GT. -1.0_wp ) ) .AND.         &
                  ( xl_source .GT. -1.0_wp ) )  THEN


                WRITE(*,*) 'ERROR: problem with namelist RADIAL_VOLUME_PARAMETERS'
                WRITE(*,*) 'DEFINE ONLY volume fractions or mass fractions'
                WRITE(*,*) 'alphag_source =' , alphag_source(1:n_add_gas)
                WRITE(*,*) 'xl_source =' , xl_source
                WRITE(*,*) 'Please check the input file'
                STOP

             END IF

             IF ( ( ANY(alphas_source(1:n_solid) .GT. -1.0_wp ) ) .AND.         &
                  ( xl_source .GT. -1.0_wp ) )  THEN

                WRITE(*,*) 'ERROR: problem with namelist RADIAL_VOLUME_PARAMETERS'
                WRITE(*,*) 'DEFINE ONLY volume fractions or mass fractions'
                WRITE(*,*) 'alphas_source =' , alphas_source(1:n_solid)
                WRITE(*,*) 'xl_source =' , xl_source
                WRITE(*,*) 'Please check the input file'
                STOP

             END IF

             IF ( ( ANY(xg_source(1:n_add_gas) .GT. -1.0_wp ) ) .AND.         &
                  ( alphal_source .GT. -1.0_wp ) )  THEN


                WRITE(*,*) 'ERROR: problem with namelist RADIAL_VOLUME_PARAMETERS'
                WRITE(*,*) 'DEFINE ONLY volume fractions or mass fractions'
                WRITE(*,*) 'xg_source =' , xg_source(1:n_add_gas)
                WRITE(*,*) 'alphal_source =' , alphal_source
                WRITE(*,*) 'Please check the input file'
                STOP

             END IF

             IF ( ( ANY(xs_source(1:n_solid) .GT. -1.0_wp ) ) .AND.         &
                  ( alphal_source .GT. -1.0_wp ) )  THEN

                WRITE(*,*) 'ERROR: problem with namelist RADIAL_VOLUME_PARAMETERS'
                WRITE(*,*) 'DEFINE ONLY volume fractions or mass fractions'
                WRITE(*,*) 'xs_source =' , xs_source(1:n_solid)
                WRITE(*,*) 'alphal_source =' , alphal_source
                WRITE(*,*) 'Please check the input file'
                STOP

             END IF

             IF ( ( alphal_source .GT. -1.0_wp ) .AND. ( xl_source .GT. -1.0_wp ) )  THEN

                WRITE(*,*) 'ERROR: problem with namelist RADIAL_VOLUME_PARAMETERS'
                WRITE(*,*) 'DEFINE ONLY liquid volume fract. or liquid mass fract.'
                WRITE(*,*) 'alphal_source =' , alphal_source
                WRITE(*,*) 'xl_source =' , xl_source
                WRITE(*,*) 'Please check the input file'
                STOP

             END IF

          END IF



          IF ( ANY(time_param .LT. 0.0_wp ) ) THEN

             WRITE(*,*)
             WRITE(*,*) 'WARNING: problem with namelist RADIAL_SOURCEPARAMETERS'
             WRITE(*,*) 'time_param =' , time_param
             time_param(1) = t_end
             time_param(2) = t_end
             time_param(3) = 0.0_wp
             time_param(4) = t_end
             WRITE(*,*) 'CHANGED TO time_param =',time_param
             WRITE(*,*) 'Radial source now constant in time' 
             WRITE(*,*)

          ELSE

             IF ( time_param(2) .GT. time_param(1) ) THEN

                WRITE(*,*) 'ERROR: problem with namelist RADIAL_SOURCEPARAMETERS'
                WRITE(*,*) 'time_param(1),time_param(2) =' , time_param(1:2)
                WRITE(*,*) 'time_param(1) must be larger than time_param(2)'
                STOP         

             END IF

             IF ( time_param(3) .GT. ( 0.5_wp*time_param(2) ) ) THEN

                WRITE(*,*) 'ERROR: problem with namelist RADIAL_SOURCE_PARAMETERS'
                WRITE(*,*) 'time_param(3) =', time_param(3)
                WRITE(*,*) 'time_param(3) must be smaller than 0.5*time_param(2)'
                STOP

             END IF


          END IF

          ALLOCATE( qp_source(n_vars+2) )

          IF ( lateral_source_flag .OR. radial_source_flag ) THEN

             ! compute the velocity, given the Richardson number

             ! fix an initial velocity 
             IF ( Ri_source .NE. -1.0_wp ) vel_source = 1.0_wp

             IF ( mfr_source .GE. 0.0_wp ) h_source = 1.0_wp

             ! define the physical variable from the source values and the
             ! initial velocity
             qp_source(1) = h_source
             qp_source(2) = h_source*vel_source
             qp_source(3) = 0.0_wp

             qp_source(4) = T_source

             IF ( ANY(xs_source(1:n_solid) .GT. -1.0_wp ) ) THEN

                ALLOCATE( inv_rho_g(n_add_gas) )

                xs_tot = SUM(xs_source(1:n_solid))

                IF ( gas_flag .AND. liquid_flag ) THEN

                   ! compute carrier phase (gas) mass fraction
                   xc = 1.0_wp - xs_tot - xl_source

                ELSE

                   ! compute carrier phase (gas or liquid) mass fraction
                   xc = 1.0_wp - xs_tot

                END IF

                IF ( gas_flag ) THEN

                   ! carrier phase is gas
                   sp_gas_const_c = ( ( xc - SUM( xg_source(1:n_add_gas) ) ) * sp_gas_const_a      &
                        + DOT_PRODUCT( xg_source(1:n_add_gas) , sp_gas_const_g(1:n_add_gas) ) )    &
                        / xc

                   inv_rho_c = sp_gas_const_c * T_source * inv_pres
                   inv_rho_g(1:n_add_gas) = sp_gas_const_g(1:n_add_gas) * T_source * inv_pres

                ELSE

                   inv_rho_c = inv_rho_l

                END IF

                inv_rhom = DOT_PRODUCT( xs_source(1:n_solid) , inv_rho_s(1:n_solid) )            &
                     + xc * inv_rho_c

                IF ( gas_flag .AND. liquid_flag ) inv_rhom = inv_rhom + xl_source * inv_rho_l

                rho_m = 1.0_wp / inv_rhom

                ! convert from mass fraction to volume fraction
                alphas_source(1:n_solid) = rho_m * xs_source(1:n_solid) * inv_rho_s(1:n_solid)

                ! convert from mass fraction to volume fraction
                alphag_source(1:n_add_gas) = rho_m * xg_source(1:n_solid) * inv_rho_g(1:n_add_gas)

                IF ( liquid_flag ) THEN

                   ! convert from mass fraction to volume fraction
                   alphal_source = rho_m * xl_source * inv_rho_l

                END IF

             END IF
             
             WRITE(*,*) 'Source solid volume fraction =',alphas_source(1:n_solid)

             IF ( n_add_gas .GT. 0 ) THEN

                WRITE(*,*) 'Source add.gas volue fraction =' ,                  &
                     alphag_source(1:n_add_gas)

             END IF

             IF ( gas_flag .AND. liquid_flag ) THEN

                WRITE(*,*) 'Source liquid volume fraction =',alphal_source

             END IF

             IF ( alpha_flag ) THEN

                qp_source(5:4+n_solid) = alphas_source(1:n_solid)
                qp_source(4+n_solid+1:4+n_solid+n_add_gas) =                    &
                     alphag_source(1:n_add_gas)

                IF ( gas_flag .AND. liquid_flag ) qp_source(n_vars) =           &
                     alphal_source

             ELSE

                qp_source(5:4+n_solid) = alphas_source(1:n_solid) * h_source
                qp_source(4+n_solid+1:4+n_solid+n_add_gas) =                    &
                     alphag_source(1:n_add_gas) * h_source

                IF ( gas_flag .AND. liquid_flag ) qp_source(n_vars) =           &
                     alphal_source * h_source

             END IF


             qp_source(n_vars+1) = vel_source
             qp_source(n_vars+2) = 0.0_wp

             sp_heat_flag = .FALSE.

             ! compute the Richardson number for vel = 1.0 (or vel from input)
             CALL mixt_var( qp_source, Ri, rho_m, rho_c, red_grav, sp_heat_flag,&
                  sp_heat_c, sp_heat_mix)

             WRITE(*,*) 'Source density =',rho_m,'(kg/m3)'
             WRITE(*,*) 'Source reduced gravity =',red_grav,'(m/s2)'
             
             ALLOCATE( shape_coeff(n_vars) )
             
             IF ( mfr_source .GE. 0.0_wp ) THEN

                IF ( Ri_source .GE. 0.0_wp ) THEN
                   ! mfr and Ri given as input

                   shape_coeff(1:n_vars) = 1.0_wp
                   
                   IF ( vertical_profiles_flag ) THEN
                      
                      iter_max = 500 

                   ELSE

                      iter_max = 1
                      
                   END IF
                                      
                   search_Ri_loop:DO iter_source = 1,iter_max

                      ! WRITE(*,*) 'iter_source',iter_source
                      
                      h_source = ( mfr_source/( source_length * rho_m * shape_coeff(1) ) &
                           * SQRT( Ri_source/red_grav ) )**(2.0_wp/3.0_wp)
                      
                      qp_source(1) = h_source
                      qp_source(2) = h_source*vel_source
                      qp_source(3) = 0.0_wp
                      
                      qp_source(4) = T_source
                      
                      IF ( alpha_flag ) THEN
                         
                         qp_source(5:4+n_solid) = alphas_source(1:n_solid)
                         qp_source(4+n_solid+1:4+n_solid+n_add_gas) =                    &
                              alphag_source(1:n_add_gas)
                         
                         IF ( gas_flag .AND. liquid_flag ) qp_source(n_vars) =           &
                              alphal_source
                         
                      ELSE
                         
                         qp_source(5:4+n_solid) = alphas_source(1:n_solid) * h_source
                         qp_source(4+n_solid+1:4+n_solid+n_add_gas) =                    &
                              alphag_source(1:n_add_gas) * h_source
                         
                         IF ( gas_flag .AND. liquid_flag ) qp_source(n_vars) =           &
                              alphal_source * h_source
                         
                      END IF
                      
                      qp_source(n_vars+1) = vel_source
                      qp_source(n_vars+2) = 0.0_wp
                      
                      ! compute the Richardson number for vel = 1.0
                      CALL mixt_var( qp_source, Ri, rho_m, rho_c, red_grav,        &
                           sp_heat_flag, sp_heat_c, sp_heat_mix)
                      
                      ! compute the correct velocity for the desired Richardson number
                      vel_source = SQRT( red_grav * h_source /Ri_source )
                      
                      IF ( ABS( Ri-Ri_source ) < 1.0E-15_wp ) EXIT search_Ri_loop

                      ! READ(*,*)
                                            
                      IF ( vertical_profiles_flag ) THEN
                         
                         CALL eval_flux_coeffs( qp_source , 0.0_wp , 0.0_wp , rho_c ,    &
                              rho_m , shape_coeff )
                         
                         ! WRITE(*,*) 'Mass flux coeff.',shape_coeff(1)
                         
                      ELSE
                         
                         shape_coeff(1:n_vars) = 1.0_wp
                         
                      END IF
                      
                   END DO search_Ri_loop

                ELSE

                   ! mfr and velocity

                   CALL mixt_var( qp_source, Ri, rho_m, rho_c, red_grav,        &
                        sp_heat_flag, sp_heat_c, sp_heat_mix)

                   IF ( vertical_profiles_flag ) THEN
                      
                      CALL eval_flux_coeffs( qp_source , 0.0_wp , 0.0_wp , rho_c ,    &
                           rho_m , shape_coeff )
                      
                      WRITE(*,*) 'qp',qp_source
                      WRITE(*,*) 'Mass flux coeff.',shape_coeff(1)
                      
                   ELSE
                      
                      shape_coeff(1:n_vars) = 1.0_wp
                      
                   END IF
                   
                   h_source = mfr_source / ( source_length * rho_m * &
                        vel_source * shape_coeff(1) )

                END IF

             ELSE

                IF ( Ri_source .GE. 0.0_wp ) THEN

                   ! compute the correct velocity for the desired Richardson number
                   vel_source = SQRT( Ri/Ri_source )

                END IF

             END IF

             qp_source(1) = h_source
             qp_source(2) = h_source*vel_source
             qp_source(n_vars+1) = vel_source

             ! Check that the Richardson number is correct
             CALL mixt_var( qp_source, Ri, rho_m, rho_c, red_grav, sp_heat_flag,&
                  sp_heat_c, sp_heat_mix)
             
             IF ( vertical_profiles_flag ) THEN

                CALL eval_flux_coeffs( qp_source , 0.0_wp , 0.0_wp , rho_c ,    &
                     rho_m , shape_coeff )

                WRITE(*,*) 'Mass flux coeff.',shape_coeff(1)

             ELSE

                shape_coeff(1:n_vars) = 1.0_wp

             END IF
                
             WRITE(*,*) 'Source Richardson number =',Ri
             WRITE(*,*) 'Source velocity =',vel_source,' (m/s)'
             WRITE(*,*) 'Source thickness =',h_source,'(m)'

             mfr = rho_m * h_source * source_length * vel_source * shape_coeff(1)
             WRITE(*,*) 'Source mass flow rate =',mfr,'(kg/s)'
             WRITE(*,*)
             ! READ(*,*)

          END IF

       END IF

    END IF source_flag

    ! ------- READ runout_parameters NAMELIST -----------------------------------

    IF ( output_runout_flag ) THEN

       READ(input_unit, runout_parameters,IOSTAT=ios)

       IF ( ios .NE. 0 ) THEN

          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist RUNOUT_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          WRITE(*,runout_parameters) 
          STOP

       ELSE

          REWIND(input_unit)

       END IF

       IF ( ( x0_runout .EQ. -1.0_wp ) .AND. ( y0_runout .EQ. -1.0_wp ) ) THEN

          WRITE(*,*) 'Runout reference location not defined'

          IF ( collapsing_volume_flag ) THEN

             x0_runout = x_collapse
             y0_runout = y_collapse
             WRITE(*,*) 'New reference location defined from collapse (x,y)'
             WRITE(*,*) 'x0_runout =',x0_runout
             WRITE(*,*) 'y0_runout =',y0_runout

          END IF

          IF ( radial_source_flag .OR. bottom_radial_source_flag ) THEN

             x0_runout = x_source
             y0_runout = y_source
             WRITE(*,*) 'New reference location defined from source'
             WRITE(*,*) 'x0_runout =',x0_runout
             WRITE(*,*) 'y0_runout =',y0_runout

          END IF

       ELSE

          IF ( x0_runout .LT. x0 ) THEN

             WRITE(*,*) 'Computational domain problem'
             WRITE(*,*) 'x0_runout < x0',x0,x0_runout
             STOP

          END IF

          IF ( x0 .GT. x0+comp_cells_x*cell_size ) THEN

             WRITE(*,*) 'Computational domain problem'
             WRITE(*,*) 'x0_runout > x0+comp_cells_x*cell_size' , x0 ,          &
                  x0_runout+comp_cells_x*cell_size
             STOP

          END IF

          IF ( y0_runout .LT. y0 ) THEN

             WRITE(*,*) 'Computational domain problem'
             WRITE(*,*) 'y0_runout < y0',y0,y0_runout
             STOP

          END IF

          IF ( y0 .GT. y0+comp_cells_y*cell_size ) THEN

             WRITE(*,*) 'Computational domain problem'
             WRITE(*,*) 'y0_runout > y0+comp_cells_y*cell_size' , y0 ,          &
                  y0_runout+comp_cells_y*cell_size
             STOP

          END IF

       END IF

       runout_file = TRIM(run_name)//'_runout'//'.csv'

       OPEN(runout_unit,FILE=runout_file,STATUS='unknown',form='formatted')

       mass_center_file = TRIM(run_name)//'_mass_center'//'.csv'

       OPEN(mass_center_unit,FILE=mass_center_file,STATUS='unknown',form='formatted')



    END IF

    !------ search for check points --------------------------------------------

    REWIND(input_unit)

    tend1 = .FALSE.

    IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'Searching for probes coords'

    n_probes = 0

    probes_search: DO

       READ(input_unit,*, END = 300 ) card

       IF( TRIM(card) == 'PROBES_COORDS' ) THEN

          EXIT probes_search

       END IF

    END DO probes_search


    READ(input_unit,*) n_probes

    WRITE(*,*) 'n_probes ',n_probes

    READ(input_unit,*) dt_probes

    WRITE(*,*) 'dt_probes ',dt_probes

    ALLOCATE( probes_coords( 2 , n_probes ) )

    DO k = 1, n_probes

       READ(input_unit,*) probes_coords( 1:2 , k ) 

       IF ( verbose_level.GE.0 ) WRITE(*,*) k , probes_coords( 1:2 , k )  

    END DO

    GOTO 310
300 tend1 = .TRUE.
310 CONTINUE

    ! ----- end search for check points -----------------------------------------

    CLOSE( input_unit )

    bak_name = TRIM(run_name)//'.bak'

    OPEN(backup_unit,file=bak_name,status='unknown',delim='quote')

    WRITE(backup_unit, run_parameters )

    IF ( restart ) THEN

       WRITE(backup_unit,newrun_parameters)
       WRITE(backup_unit,restart_parameters)

    ELSE

       WRITE(backup_unit,newrun_parameters)

       IF ( ( radial_source_flag ) .OR. ( bottom_radial_source_flag ) ) THEN

          alphal_source = -1.0_wp

          WRITE(backup_unit,radial_source_parameters)

       ELSEIF ( lateral_source_flag ) THEN

          WRITE(backup_unit,lateral_source_parameters)
          
       ELSE

          WRITE(backup_unit,initial_conditions)

       END IF

    END IF

    IF ( comp_cells_x .GT. 1 ) THEN

       WRITE(backup_unit,west_boundary_conditions)
       WRITE(backup_unit,east_boundary_conditions)

    END IF

    IF ( comp_cells_y .GT. 1 ) THEN

       WRITE(backup_unit,north_boundary_conditions)
       WRITE(backup_unit,south_boundary_conditions)

    END IF

    WRITE(backup_unit, numeric_parameters )

    WRITE(backup_unit, expl_terms_parameters )

    IF ( rheology_flag ) WRITE(backup_unit,rheology_parameters)

    ! WRITE(backup_unit,temperature_parameters)

    WRITE(backup_unit,solid_transport_parameters)
    WRITE(backup_unit,gas_transport_parameters)

    IF ( liquid_flag ) WRITE(backup_unit,liquid_transport_parameters)

    IF ( radial_source_flag ) WRITE(backup_unit,radial_source_parameters)

    IF ( output_runout_flag ) WRITE(backup_unit, runout_parameters)

    IF ( ( COUNT( thickness_levels0 .GT. 0.0_wp ) .GT. 0 ) .OR.                 &
         ( COUNT( dyn_pres_levels0 .GT. 0.0_wp ) .GT. 0 ) ) THEN

       WRITE(backup_unit, vulnerability_table_parameters)

    END IF

    IF ( vertical_profiles_flag ) WRITE(backup_unit, vertical_profiles_parameters)

    IF ( n_probes .GT. 0 ) THEN

       WRITE(backup_unit,*) '''PROBES_COORDS'''
       WRITE(backup_unit,*) n_probes
       WRITE(backup_unit,*) dt_probes

       DO k = 1,n_probes

          WRITE(backup_unit,109) probes_coords(1:2,k)

109       FORMAT(2(1x,e14.7))

       END DO

    END IF

    CLOSE(backup_unit)

  END SUBROUTINE read_param

  !******************************************************************************
  !> \brief Read the input file
  !
  !> This subroutine read the input parameters from the file 
  !> "two_phases.inp" and write a backup file of the input parameters 
  !> with name "run_name.bak", where run_name is read from the input
  !> file.
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE update_param

    IMPLICIT none

    INTEGER :: ios

    CHARACTER(LEN=40) :: run_name_org
    LOGICAL :: restart_org
    REAL(wp) :: t_start_org
    REAL(wp) :: t_end_org
    REAL(wp) :: dt_output_org
    LOGICAL :: output_cons_flag_org
    LOGICAL :: output_phys_flag_org
    LOGICAL :: output_esri_flag_org
    INTEGER :: verbose_level_org

    run_name_org = run_name
    restart_org = restart
    t_start_org = t_start
    t_end_org = t_end
    dt_output_org = dt_output
    output_cons_flag_org = output_cons_flag
    output_phys_flag_org = output_phys_flag
    output_esri_flag_org = output_esri_flag
    verbose_level_org = verbose_level


    OPEN(input_unit,FILE=input_file,STATUS='old')

    ! ------- READ run_parameters NAMELIST -----------------------------------
    READ(input_unit, run_parameters,IOSTAT=ios )

    IF ( ios .NE. 0 ) THEN

       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist RUN_PARAMETERS'
       WRITE(*,*) 'Please check the input file'
       STOP

    ELSE

       REWIND(input_unit)

    END IF

    CLOSE(input_unit)

    IF ( t_end_org .NE. t_end ) THEN

       WRITE(*,*) 'Modified input file: t_end =',t_end

    END IF

    IF ( (.NOT.output_cons_flag) .AND. (.NOT.output_esri_flag) .AND.            &
         (.NOT.output_phys_flag) ) THEN

       dt_output = 2.0 * ( t_end - t_start ) 

    ELSE

       IF ( dt_output_org .NE. dt_output ) THEN

          WRITE(*,*) 'Modified input file: dt_output =',dt_output

       END IF

    END IF

    IF ( output_cons_flag_org .NEQV. output_cons_flag ) THEN

       WRITE(*,*)  'Modified input file: output_cons_flag =',output_cons_flag

    END IF

    IF ( output_phys_flag_org .NEQV. output_phys_flag ) THEN

       WRITE(*,*)  'Modified input file: output_phys_flag =',output_phys_flag

    END IF

    IF ( output_esri_flag_org .NEQV. output_esri_flag ) THEN

       WRITE(*,*)  'Modified input file: output_esri_flag =',output_esri_flag

    END IF

    IF ( verbose_level_org .NE. verbose_level ) THEN

       WRITE(*,*)  'Modified input file: verbose_level =',verbose_level

    END IF

    run_name_org = run_name
    restart_org = restart
    t_start_org = t_start

    CLOSE(input_unit)


  END SUBROUTINE update_param


  !******************************************************************************
  !> \brief Read the solution from the restart unit
  !
  !> This subroutine is called when the parameter "restart" in the input 
  !> file is TRUE. Then the initial solution is read from a file. 
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE read_solution

    ! External procedures
    USE geometry_2d, ONLY : interp_2d_scalarB , regrid_scalar
    USE solver_2d, ONLY : allocate_solver_variables , solve_mask_time

    ! External variables
    USE geometry_2d, ONLY : comp_cells_x , x0 , comp_cells_y , y0 , dx , dy
    USE geometry_2d, ONLY : B_cent , erodible
    USE init_2d, ONLY : thickness_init , erodible_init
    USE parameters_2d, ONLY : n_vars
    USE solver_2d, ONLY : q

    IMPLICIT none

    CHARACTER(LEN=15) :: chara

    INTEGER :: j,k

    INTEGER :: dot_idx

    LOGICAL :: lexist

    CHARACTER(LEN=30) :: string

    CHARACTER(LEN=3) :: check_file

    INTEGER :: ncols , nrows , nodata_value

    REAL(wp) :: xllcorner , yllcorner , cellsize

    REAL(wp) :: xj , yk

    REAL(wp) :: thickness_interp

    REAL(wp), ALLOCATABLE :: thickness_input(:,:)

    REAL(wp), ALLOCATABLE :: x1(:) , y1(:)

    REAL(wp) :: xl , xr , yl , yr 

    REAL(wp) :: rho_c , rho_m , mass_fract(n_solid)

    REAL(wp) :: sp_heat_c

    INTEGER :: solid_idx

    INTEGER :: i_vars , i_solid

    INTEGER :: i_file

    ALLOCATE( thickness_init(comp_cells_x,comp_cells_y) )

    thickness_init(:,:) = 0.0_wp

    restart_file = restart_files(1)
 
    check_file = restart_file(dot_idx+1:dot_idx+3)
    
    IF ( check_file .EQ. 'asc' ) THEN
       
       IF ( liquid_flag .AND. gas_flag ) THEN
          
          WRITE(*,*) 'ERROR: problem with namelist NEWRUN_PARAMETERS'
          WRITE(*,*) 'When restarting from .asc file only'
          WRITE(*,*) 'one of these parameters must be set to .TRUE.'          
          WRITE(*,*) 'LIQUID_FLAG',liquid_flag
          WRITE(*,*) 'GAS_FLAG',liquid_flag
          WRITE(*,*) 'Please check the input file'
          CLOSE(restart_unit)
          STOP
          
       ELSEIF ( ( .NOT.liquid_flag ) .AND. ( .NOT. gas_flag ) ) THEN
          
          WRITE(*,*) 'ERROR: problem with namelist NEWRUN_PARAMETERS'
          WRITE(*,*) 'When restarting from .asc file only one'
          WRITE(*,*) 'of these parameters must be set to .TRUE.'
          WRITE(*,*) 'LIQUID_FLAG',liquid_flag
          WRITE(*,*) 'GAS_FLAG',liquid_flag
          WRITE(*,*) 'Please check the input file'
          CLOSE(restart_unit)
          STOP
          
       ELSEIF ( gas_flag ) THEN
          
          IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'Carrier phase: gas'
          
       ELSEIF ( liquid_flag ) THEN
          
          IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'Carrier phase: liquid'
          
       END IF
       
    END IF

    WRITE(*,*) 'n_restart_files',n_restart_files
      
    DO i_file=1,n_restart_files

       restart_file = restart_files(i_file)

       INQUIRE (FILE=restart_file,exist=lexist)
       
       WRITE(*,*)
       
       IF ( lexist .EQV. .FALSE.) THEN
          
          WRITE(*,*) 'Restart: ',TRIM(restart_file) , ' not found'
          STOP
          
       ELSE
          
          OPEN(restart_unit,FILE=restart_file,STATUS='old')
          
          IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'Restart: ',TRIM(restart_file),   &
               ' found'
  
          WRITE(*,*) 'Release time: ',release_time(i_file)
        
       END IF

       dot_idx = SCAN(restart_file, ".", .TRUE.)
       
       check_file = restart_file(dot_idx+1:dot_idx+3)
       
       IF ( check_file .EQ. 'asc' ) THEN
          
          READ(restart_unit,*) chara, ncols
          READ(restart_unit,*) chara, nrows
          READ(restart_unit,*) chara, xllcorner
          READ(restart_unit,*) chara, yllcorner
          READ(restart_unit,*) chara, cellsize
          READ(restart_unit,*) chara, nodata_value
          
          ALLOCATE( thickness_input(ncols,nrows) )
          
          IF ( ( xllcorner - x0 ) .GT. 1.E-5_wp*cellsize ) THEN
             
             WRITE(*,*)
             WRITE(*,*) 'WARNING: initial solution and domain extent'
             WRITE(*,*) 'xllcorner greater than x0', xllcorner , x0
             
          END IF
          
          IF ( ( yllcorner - y0 ) .GT. 1.E-5_wp*cellsize ) THEN
             
             WRITE(*,*)
             WRITE(*,*) 'WARNING: initial solution and domain extent'
             WRITE(*,*) 'yllcorner greater then y0', yllcorner , y0
             
          END IF
          
          IF ( x0+cell_size*(comp_cells_x+1) - ( xllcorner+cellsize*(ncols+1) )    &
               .GT. 1.E-5_wp*cellsize ) THEN
             
             WRITE(*,*)
             WRITE(*,*) 'WARNING: initial solution and domain extent'
             WRITE(*,*) 'xrrcorner greater than ', xllcorner , x0
             
          END IF
          
          IF ( x0+cell_size*(comp_cells_x+1) - ( xllcorner+cellsize*(ncols+1) )    &
               .GT. 1.E-5_wp*cellsize ) THEN
             
             WRITE(*,*)
             WRITE(*,*) 'WARNING: initial solution and domain extent'
             WRITE(*,*) 'yllcorner greater then y0', yllcorner , y0
             
          END IF
          
          IF ( cellsize .NE. cell_size ) THEN
             
             WRITE(*,*)
             WRITE(*,*) 'WARNING: changing resolution of restart' 
             WRITE(*,*) 'cellsize not equal to cell_size', cellsize , cell_size
             WRITE(*,*)
             
          END IF
          
          DO k=1,nrows
             
             WRITE(*,FMT="(A1,A,t21,F6.2,A)",ADVANCE="NO") ACHAR(13),              &
                  & " Percent Complete: ",( REAL(k) / REAL(nrows))*100.0, "%"
             
             READ(restart_unit,*) thickness_input(1:ncols,nrows-k+1)
             
          ENDDO
          
          WRITE(*,*) 
          
          WHERE ( thickness_input .EQ. nodata_value )
             
             thickness_input = 0.0_wp
             
          END WHERE
          
          IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'Total volume from restart =',    &
               cellsize**2*SUM(thickness_input)
          
          
          !----- NEW INITIALIZATION OF THICKNESS FROM RESTART
          ALLOCATE( x1(ncols+1) , y1(nrows+1) )
          
          DO j=1,ncols+1
             
             x1(j) = xllcorner + (j-1)*cellsize
             
          END DO
          
          DO k=1,nrows+1
             
             y1(k) = yllcorner + (k-1)*cellsize
             
          END DO
          
          DO j=1,comp_cells_x
             
             xl = x0 + (j-1)*cell_size
             xr = x0 + (j)*cell_size
             
             DO k=1,comp_cells_y
                
                yl = y0 + (k-1)*cell_size
                yr = y0 + (k)*cell_size
                
                CALL regrid_scalar( x1 , y1 , thickness_input , xl , xr , yl ,     &
                     yr , thickness_interp )
                
                thickness_init(j,k) =  thickness_init(j,k) + thickness_interp

                IF ( thickness_interp .GT. 0.0_wp ) THEN 
                   
                   ! WRITE(*,*) 'j,k,thickness: ',j,k,thickness_init(j,k)
                   solve_mask_time(j,k) = release_time(i_file)
                   
                END IF
                
             END DO
             
          END DO

          DEALLOCATE( thickness_input )
          
          DEALLOCATE( x1 , y1 )

       END IF

       CLOSE(restart_unit)
       
    END DO
   
    IF ( check_file .EQ. 'asc' ) THEN
       
       IF ( subtract_init_flag ) THEN
             
          WRITE(*,*) 'Subtricting initial thickness from DEM'
          B_cent(:,:) = B_cent(:,:) - thickness_init(:,:)
             
          IF ( erosion_coeff .GT. 0.0_wp ) THEN
                
             IF ( MAXVAL(erodible_init(:,:)) .GT. 0.0_wp ) THEN
                   
                WRITE(*,*)
                WRITE(*,*) 'Subtracting initial thickness from erodible thickness'
                erodible_init(:,:) = erodible_init(:,:) - thickness_init(:,:) *    &
                     ( SUM(alphas_init(1:n_solid))  /                              &
                     ( 1.0_wp - erodible_porosity ) )
                
                IF ( MINVAL(erodible_init(:,:)) .LT. 0.0_wp ) THEN
                   
                   WRITE(*,*) 'WARNING: MINVAL(erodible_init) = ',                 &
                        MINVAL(erodible_init(:,:)) 
                   WRITE(*,*) 'Initial erodible thick. negative values changed to 0'
                   
                   erodible_init = MAX( 0.0_wp , erodible_init )
                   
                END IF
                
                DO i_solid=1,n_solid
                   
                   erodible(i_solid,:,:) = erodible_fract(i_solid) *               &
                        ( 1.0_wp - erodible_porosity ) * erodible_init(:,:)
                   
                END DO
                
                WRITE(*,*)
                
             END IF
             
          END IF
          
       END IF
          

       !----- END NEW INITIALIZATION OF THICKNESS FROM RESTART
       
       IF ( gas_flag ) THEN
          
          rho_c = pres / ( sp_gas_const_a * T_init )
          sp_heat_c = sp_heat_a

       ELSE

          rho_c = rho_l
          sp_heat_c = sp_heat_l

       END IF

       rho_m = SUM( rho_s(1:n_solid)*alphas_init(1:n_solid) ) + ( 1.0_wp -      &
            SUM( alphas_init(1:n_solid) ) ) * rho_c 

       mass_fract = rho_s * alphas_init / rho_m

       q(1,:,:) = thickness_init(:,:) * rho_m

       IF ( VERBOSE_LEVEL .GE. 0 ) THEN

          WRITE(*,*) 'Total volume on computational grid =',cell_size**2 *      &
               SUM( thickness_init(:,:) )
          WRITE(*,*) 'Total mass on computational grid =',cell_size**2 *        &
               SUM( q(1,:,:) )

       END IF
       ! rhom*h*u
       q(2,:,:) = q(1,:,:) * u_init
       ! rhom*h*v
       q(3,:,:) = q(1,:,:) * v_init

       ! energy (total or internal)
       q(4,:,:) = 0.0_wp

       WHERE ( thickness_init .GT. 0.0_wp )

          q(4,:,:) = q(1,:,:) * T_init *  ( SUM( mass_fract(1:n_solid) *        &
               sp_heat_s(1:n_solid) ) +    &
               ( 1.0_wp - SUM( mass_fract ) ) * sp_heat_l )

       END WHERE

       DO solid_idx=5,4+n_solid

          ! rhos*h*alphas
          q(solid_idx,:,:) = 0.0_wp

          WHERE ( thickness_init .GT. 0.0_wp )

             q(solid_idx,:,:) = thickness_init(:,:) * alphas_init(solid_idx-4) *&
                  rho_s(solid_idx-4)

          END WHERE

       END DO

       WRITE(*,*) 'MAXVAL(q(5,:,:))',MAXVAL(q(5:4+n_solid,:,:))

       IF ( VERBOSE_LEVEL .GE. 0 ) THEN

          WRITE(*,*) 'Total sediment volume =',cell_size**2*SUM( thickness_init*&
               SUM(alphas_init) )

       END IF

       output_idx = 0

       IF ( verbose_level .GE. 1 ) THEN

          WRITE(*,*) 'Min q(1,:,:) =',MINVAL(q(1,:,:))
          WRITE(*,*) 'Max q(1,:,:) =',MAXVAL(q(1,:,:))
          WRITE(*,*) 'SUM(q(1,:,:)) =',SUM(q(1,:,:))

          DO k=1,nrows

             WRITE(*,*) k,B_cent(:,k)
             READ(*,*)

          END DO

          WRITE(*,*) 'SUM(B_cent(:,:)) =',SUM(B_cent(:,:))
          READ(*,*)

       END IF


       WRITE(*,*) 'n_vars',n_vars

    END IF

    IF ( check_file .EQ. 'q_2' ) THEN

       OPEN(restart_unit,FILE=restart_file,STATUS='old')

       DO k=1,comp_cells_y

          DO j=1,comp_cells_x

             READ(restart_unit,'(2e20.12,100(e20.12))') xj , yk ,               &
                  (q(i_vars,j,k),i_vars=1,n_vars) 

             IF ( q(1,j,k) .LE. 0.0_wp ) q(1:n_vars,j,k) = 0.0_wp

             DO solid_idx=5,4+n_solid

                IF ( q(solid_idx,j,k) .LE. 0.0_wp ) q(solid_idx,j,k) = 0.0_wp

             END DO

          ENDDO

          READ(restart_unit,*)  

       END DO

       IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'Total mass =',dx*dy*SUM(q(1,:,:))

       DO solid_idx=5,4+n_solid

          IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'Total sediment mass =',       &
               dx*dy* SUM( q(solid_idx,:,:) )

       END DO

       j = SCAN(restart_file, '.' , .TRUE. )

       string = TRIM(restart_file(j-4:j-1))
       READ( string,* ) output_idx

       IF ( VERBOSE_LEVEL .GE. 0 ) THEN

          WRITE(*,*) 
          WRITE(*,*) 'Starting from output index ',output_idx

       END IF

       ! Set this flag to 0 to not overwrite the initial condition

    END IF

    CLOSE(restart_unit)

  END SUBROUTINE read_solution

  !******************************************************************************
  !> \brief Read the solution from the restart unit
  !
  !> This subroutine is called when the parameter "restart" in the input 
  !> file is TRUE. Then the initial solution is read from a file. 
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE read_erodible

    ! External procedures
    USE geometry_2d, ONLY : interp_2d_scalarB , regrid_scalar
    USE solver_2d, ONLY : allocate_solver_variables

    ! External variables
    USE geometry_2d, ONLY : comp_cells_x , x0 , comp_cells_y , y0
    USE geometry_2d, ONLY : erodible
    USE init_2d, ONLY : erodible_init

    IMPLICIT none

    CHARACTER(LEN=15) :: chara

    INTEGER :: j,k

    INTEGER :: dot_idx

    LOGICAL :: lexist

    CHARACTER(LEN=3) :: check_file

    INTEGER :: ncols , nrows , nodata_value

    REAL(wp) :: xllcorner , yllcorner , cellsize

    REAL(wp), ALLOCATABLE :: erodible_input(:,:)

    REAL(wp), ALLOCATABLE :: x1(:) , y1(:)

    REAL(wp) :: xl , xr , yl , yr 

    INTEGER :: i_solid


    INQUIRE (FILE=erodible_file,exist=lexist)

    WRITE(*,*)
    
    IF ( lexist .EQV. .FALSE.) THEN
       
       WRITE(*,*) 'Erodible file: ',TRIM(erodible_file) , ' not found'
       STOP
       
    END IF
    
    OPEN(erodible_unit,FILE=erodible_file,STATUS='old')
    
    IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'Erodible file: ',                   &
         TRIM(erodible_file),' found'
    
    dot_idx = SCAN(erodible_file, ".", .TRUE.)
    
    check_file = erodible_file(dot_idx+1:dot_idx+3)
    
    IF ( check_file .NE. 'asc' ) THEN
       
       WRITE(*,*) 'Erodible file not in the right format (*.asc)'
       STOP
       
    END IF
    
    READ(erodible_unit,*) chara, ncols
    READ(erodible_unit,*) chara, nrows
    READ(erodible_unit,*) chara, xllcorner
    READ(erodible_unit,*) chara, yllcorner
    READ(erodible_unit,*) chara, cellsize
    READ(erodible_unit,*) chara, nodata_value
    
    ALLOCATE( erodible_input(ncols,nrows) )
    
    IF ( ( xllcorner - x0 ) .GT. 1.E-5_wp*cellsize ) THEN
       
       WRITE(*,*)
       WRITE(*,*) 'WARNING: initial solution and domain extent'
       WRITE(*,*) 'xllcorner greater than x0', xllcorner , x0
       
    END IF
    
    IF ( ( yllcorner - y0 ) .GT. 1.E-5_wp*cellsize ) THEN
       
       WRITE(*,*)
       WRITE(*,*) 'WARNING: initial solution and domain extent'
       WRITE(*,*) 'yllcorner greater then y0', yllcorner , y0
       
    END IF
    
    IF ( x0+cell_size*(comp_cells_x+1) - ( xllcorner+cellsize*(ncols+1) )       &
         .GT. 1.E-5_wp*cellsize ) THEN
       
       WRITE(*,*)
       WRITE(*,*) 'WARNING: initial solution and domain extent'
       WRITE(*,*) 'xrrcorner greater than ', xllcorner , x0
       
    END IF
    
    IF ( x0+cell_size*(comp_cells_x+1) - ( xllcorner+cellsize*(ncols+1) )       &
         .GT. 1.E-5_wp*cellsize ) THEN
       
       WRITE(*,*)
       WRITE(*,*) 'WARNING: initial solution and domain extent'
       WRITE(*,*) 'yllcorner greater then y0', yllcorner , y0
       
    END IF
    
    
    IF ( cellsize .NE. cell_size ) THEN
       
       WRITE(*,*)
       WRITE(*,*) 'WARNING: changing resolution of erodible layer' 
       WRITE(*,*) 'cellsize not equal to cell_size', cellsize , cell_size
       WRITE(*,*)
       
    END IF
    
    DO k=1,nrows
       
       WRITE(*,FMT="(A1,A,t21,F6.2,A)",ADVANCE="NO") ACHAR(13),                 &
            & " Percent Complete: ",( REAL(k) / REAL(nrows))*100.0, "%"
       
       READ(erodible_unit,*) erodible_input(1:ncols,nrows-k+1)
       
    ENDDO
    
    CLOSE(erodible_unit)
    
    WRITE(*,*) 
    IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'Total erodible volume =',           &
         cellsize**2*SUM(erodible_input)
    
    WHERE ( erodible_input .EQ. nodata_value )

       erodible_input = 0.0_wp

    END WHERE
    
    !----- NEW INITIALIZATION OF THICKNESS FROM RESTART
    ALLOCATE( x1(ncols+1) , y1(nrows+1) )
    
    DO j=1,ncols+1
       
       x1(j) = xllcorner + (j-1)*cellsize
       
    END DO
    
    DO k=1,nrows+1
       
       y1(k) = yllcorner + (k-1)*cellsize
       
    END DO
    
    DO j=1,comp_cells_x
       
       xl = x0 + (j-1)*cell_size
       xr = x0 + (j)*cell_size
       
       DO k=1,comp_cells_y
          
          yl = y0 + (k-1)*cell_size
          yr = y0 + (k)*cell_size
          
          CALL regrid_scalar( x1 , y1 , erodible_input , xl , xr , yl ,         &
               yr , erodible_init(j,k) )
          
       END DO
       
    END DO
    
    IF ( VERBOSE_LEVEL .GE. 0 ) THEN
       
       WRITE(*,*) 'Total erodible volume on computational grid =' ,             &
            cell_size**2 * SUM( erodible_init(:,:) )
       
    END IF
    
    RETURN


  END SUBROUTINE read_erodible


  !******************************************************************************
  !> \brief Write the solution on the output unit
  !
  !> This subroutine write the parameters of the grid, the output time 
  !> and the solution to a file with the name "run_name.q****", where  
  !> run_name is the name of the run read from the input file and ****
  !> is the counter of the output.
  !
  !> \param[in]   t      output time
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 07/10/2016
  !
  !******************************************************************************

  SUBROUTINE output_solution(time)

    ! external procedures
    USE constitutive_2d, ONLY : qc_to_qp, mixt_var , settling_velocity , vonK

    ! external variables

    USE constitutive_2d, ONLY : kin_visc_c

    USE geometry_2d, ONLY : comp_cells_x , B_cent , comp_cells_y , x_comp,      &
         y_comp , deposit , erosion , erodible , B_prime_x , B_prime_y 

    USE parameters_2d, ONLY : n_vars
    USE parameters_2d, ONLY : t_output , dt_output 
    USE parameters_2d, ONLY : t_steady

    USE solver_2d, ONLY : q , hmax , pdynmax

    IMPLICIT none

    REAL(wp), INTENT(IN) :: time

    CHARACTER(LEN=4) :: idx_string

    REAL(wp) :: qp(n_vars+2)

    REAL(wp) :: B_out

    REAL(wp) :: r_u , r_v , r_h , r_alphas(n_solid) , r_T , r_Ri , r_rho_m
    REAL(wp) :: r_alphag(n_add_gas)
    REAL(wp) :: r_alphal
    REAL(wp) :: r_rho_c      !< real-value carrier phase density [kg/m3]
    REAL(wp) :: r_red_grav   !< real-value reduced gravity
    REAL(wp) :: p_dyn

    REAL(wp) :: r_w          !< vertical component of the velocity
    REAL(wp) :: mod_vel , mod_vel2
    REAL(wp) :: shear_stress
    REAL(wp) :: shear_vel    !< shear velocity

    INTEGER :: j,k
    INTEGER :: i
    INTEGER :: i_vars
    INTEGER :: i_solid

    LOGICAL :: sp_flag
    REAL(wp) :: r_sp_heat_c
    REAL(wp) :: r_sp_heat_mix

    !> Hindered settling velocity (units: m s-1 )
    REAL(wp) :: settling_vel

    !> Inverse of kinematic viscosity of continuous phase
    REAL(wp) :: inv_kin_visc

    REAL(wp) :: Rouse_no(n_solid)

    REAL(wp) :: mod_vel_max
    
    sp_flag = .FALSE.


    output_idx = output_idx + 1

    idx_string = lettera(output_idx-1)

    IF ( output_cons_flag ) THEN

       output_file_2d = TRIM(run_name)//'_'//idx_string//'.q_2d'

       IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'WRITING ',output_file_2d

       OPEN(output_unit_2d,FILE=output_file_2d,status='unknown',form='formatted')

       !WRITE(output_unit_2d,1002) x0,dx,comp_cells_x,y0,dy,comp_cells_y,t

       DO k = 1,comp_cells_y

          DO j = 1,comp_cells_x

             DO i = 1,n_vars

                ! Exponents with more than 2 digits cause problems reading
                ! into matlab... reset tiny values to zero:
                IF ( abs(q(i,j,k)) .LT. 1.0E-20_wp ) q(i,j,k) = 0.0_wp

             ENDDO

             WRITE(output_unit_2d,'(2e20.12,100(e20.12))') x_comp(j), y_comp(k),&
                  (q(i_vars,j,k),i_vars=1,n_vars) 

          ENDDO

          WRITE(output_unit_2d,*) ' ' 

       END DO

       WRITE(output_unit_2d,*) ' '
       WRITE(output_unit_2d,*) ' '

       CLOSE(output_unit_2d)

    END IF

    IF ( output_phys_flag ) THEN

       output_file_2d = TRIM(run_name)//'_'//idx_string//'.p_2d'

       IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'WRITING ',output_file_2d

       OPEN(output_unit_2d,FILE=output_file_2d,status='unknown',form='formatted')

       r_alphal = 0.0_wp

       DO k = 1,comp_cells_y

          DO j = 1,comp_cells_x

             CALL qc_to_qp(q(1:n_vars,j,k) , qp(1:n_vars+2) , p_dyn )

             CALL mixt_var(qp(1:n_vars+2),r_Ri,r_rho_m,r_rho_c,r_red_grav,      &
                  sp_flag,r_sp_heat_c,r_sp_heat_mix)

             r_h = qp(1)
             r_u = qp(n_vars+1)
             r_v = qp(n_vars+2)
             r_T = qp(4)

             IF ( slope_correction_flag ) THEN

                r_w = r_u * B_prime_x(j,k) + r_v * B_prime_y(j,k)

             ELSE

                r_w = 0.0_wp

             END IF

             mod_vel2 = r_u**2 + r_v**2 + r_w**2
             mod_vel = SQRT( mod_vel2 )

             IF ( rheology_model .EQ. 8 ) THEN

                shear_stress = r_rho_m * friction_factor * mod_vel2

                shear_vel = SQRT( shear_stress / r_rho_m ) 

                ! Viscosity read from input file [m2 s-1]
                inv_kin_visc = 1.0_wp / kin_visc_c

                DO i_solid=1,n_solid

                   settling_vel = settling_velocity( diam_s(i_solid) ,          &
                        rho_s(i_solid) , r_rho_c , inv_kin_visc )

                   IF ( shear_vel .GT. 0.0_wp ) THEN

                      Rouse_no(i_solid) = settling_vel / ( vonK * shear_vel )

                   ELSE

                      Rouse_no(i_solid) = 0.0_wp

                   END IF

                END DO

             ELSE

                Rouse_no(1:n_solid) = 0.0_wp
                shear_vel = 0.0_wp

             END IF

             IF ( r_h .GT. 0.0_wp ) THEN

                IF ( alpha_flag ) THEN

                   r_alphas(1:n_solid) = qp(5:4+n_solid)
                   r_alphag(1:n_add_gas) = qp(4+n_solid+1:4+n_solid+n_add_gas)

                   IF ( gas_flag .AND. liquid_flag ) THEN

                      r_alphal = qp(n_vars)

                   END IF

                ELSE

                   r_alphas(1:n_solid) = qp(5:4+n_solid) / r_h
                   r_alphag(1:n_add_gas) = qp(4+n_solid+1:4+n_solid+n_add_gas)  &
                        / r_h

                   IF ( gas_flag .AND. liquid_flag ) THEN

                      r_alphal = qp(n_vars) / r_h

                   END IF


                END IF

             ELSE

                r_alphas(1:n_solid) = 0.0_wp
                r_alphag(1:n_add_gas) = 0.0_wp
                r_alphal = 0.0_wp

             END IF

             IF ( ABS( r_h ) .LT. 1.0E-20_wp ) r_h = 0.0_wp
             IF ( ABS( r_u ) .LT. 1.0E-20_wp ) r_u = 0.0_wp
             IF ( ABS( r_v ) .LT. 1.0E-20_wp ) r_v = 0.0_wp
             IF ( ABS(B_cent(j,k)) .LT. 1.0E-20_wp ) THEN 

                B_out = 0.0_wp

             ELSE

                B_out = B_cent(j,k)

             END IF

             DO i=1,n_solid

                IF ( ABS( r_alphas(i) ) .LT. 1.0E-20_wp ) r_alphas(i) = 0.0_wp
                IF ( ABS( DEPOSIT(j,k,i) ) .LT. 1.0E-20_wp )                    &
                     DEPOSIT(j,k,i) = 0.0_wp 
                IF ( ABS( EROSION(j,k,i) ) .LT. 1.0E-20_wp )                    &
                     EROSION(j,k,i) = 0.0_wp 

             END DO

             IF ( ABS( r_T ) .LT. 1.0E-20_wp ) r_T = 0.0_wp
             IF ( ABS( r_rho_m ) .LT. 1.0E-20_wp ) r_rho_m = 0.0_wp
             IF ( ABS( r_red_grav ) .LT. 1.0E-20_wp ) r_red_grav = 0.0_wp

             IF ( ABS( r_alphal ) .LT. 1.0E-20_wp ) r_alphal = 0.0_wp

             mod_vel_max = SQRT(2.0 * pdynmax(j,k) / r_rho_m)
             
             WRITE(output_unit_2d,1010) x_comp(j), y_comp(k), r_h , r_u , r_v , &
                  B_out , r_h + B_out , r_alphas , r_alphag , r_T , r_rho_m ,   &
                  r_red_grav , DEPOSIT(j,k,:) , EROSION(j,k,:) ,                &
                  SUM(ERODIBLE(1:n_solid,j,k)) / ( 1.0_wp - erodible_porosity ),&
                  r_alphal , shear_vel , r_Ri , Rouse_no(1:n_solid), hmax(j,k) ,&
                  pdynmax(j,k) , mod_vel_max

          END DO

          WRITE(output_unit_2d,*) ' ' 

       ENDDO

       WRITE(output_unit_2d,*) ' '
       WRITE(output_unit_2d,*) ' '

       CLOSE(output_unit_2d)

    END IF

1010 FORMAT(100ES15.7E2)

    t_output = time + dt_output

    IF ( output_esri_flag ) THEN

       CALL output_esri(output_idx)

       IF ( ( time .GE. t_end ) .OR. ( time .GE. t_steady ) ) THEN

          CALL output_max

       END IF

    END IF

  END SUBROUTINE output_solution


  !******************************************************************************
  !> \brief Write the maximum thickness in ESRI format
  !
  !> This subroutine write the maximum thickness in the ascii ESRI format. 
  !> A masking is applied to the region with thickness less than 1E-5.
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 08/12/2018
  !
  !******************************************************************************

  SUBROUTINE output_max

    USE geometry_2d, ONLY : grid_output , grid_output_int
    USE solver_2d, ONLY : hmax , vuln_table

    IMPLICIT NONE

    CHARACTER(LEN=4) :: idx_string

    INTEGER :: j
    INTEGER :: i_pdyn_lev , i_thk_lev , i_table

    !Save max thickness
    output_max_file = TRIM(run_name)//'_max.asc'

    IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'WRITING ',output_max_file

    OPEN(output_max_unit,FILE=output_max_file,status='unknown',form='formatted')

    grid_output = -9999 

    WHERE ( hmax(:,:).GE. 1.E-5_wp )

       grid_output = hmax(:,:) 

    END WHERE

    WRITE(output_max_unit,'(A,I5)') 'ncols ', comp_cells_x
    WRITE(output_max_unit,'(A,I5)') 'nrows ', comp_cells_y
    WRITE(output_max_unit,'(A,F15.3)') 'xllcorner ', x0
    WRITE(output_max_unit,'(A,F15.3)') 'yllcorner ', y0
    WRITE(output_max_unit,'(A,F15.3)') 'cellsize ', cell_size
    WRITE(output_max_unit,'(A,I5)') 'NODATA_value ', -9999

    DO j = comp_cells_y,1,-1

       WRITE(output_max_unit,'(2000ES12.3E3)') grid_output(1:comp_cells_x,j)

    ENDDO

    CLOSE(output_max_unit)

    i_table = 0

    IF ( n_thickness_levels * n_thickness_levels .GT. 1 ) THEN

       output_VT_file = TRIM(run_name)//'_VT.txt'
       OPEN( output_VT_unit , FILE=output_VT_file , status='unknown' ,          &
            form='formatted')

       WRITE(output_VT_unit,*) 'ID file        thickness (m)            dynamic &
            &pressure (Pa)'

       DO i_thk_lev=1,n_thickness_levels

          DO i_pdyn_lev=1,n_dyn_pres_levels

             i_table = i_table + 1

             idx_string = lettera(i_table)

             WRITE(output_VT_unit,*) idx_string ,'      ',                      &
                  thickness_levels(i_thk_lev) , dyn_pres_levels(i_pdyn_lev)

             grid_output_int(:,:) = MERGE(1,-9999,vuln_table(i_table,:,:))

             output_max_file = TRIM(run_name)//'_VT_'//idx_string//'.asc'
             OPEN(output_max_unit,FILE=output_max_file,status='unknown' ,       &
                  form='formatted')

             WRITE(output_max_unit,'(A,I5)') 'ncols ', comp_cells_x
             WRITE(output_max_unit,'(A,I5)') 'nrows ', comp_cells_y
             WRITE(output_max_unit,'(A,F15.3)') 'xllcorner ', x0
             WRITE(output_max_unit,'(A,F15.3)') 'yllcorner ', y0
             WRITE(output_max_unit,'(A,F15.3)') 'cellsize ', cell_size
             WRITE(output_max_unit,'(A,I5)') 'NODATA_value ', -9999

             DO j = comp_cells_y,1,-1

                WRITE(output_max_unit,'(2000I7)') grid_output_int(1:comp_cells_x,j)

             ENDDO

             CLOSE(output_max_unit)

          END DO

       END DO

       CLOSE(output_VT_unit)

    END IF

    RETURN

  END SUBROUTINE output_max

  !******************************************************************************
  !> \brief Write the thickness in ESRI format
  !
  !> This subroutine write the thickness in the ascii ESRI format. 
  !> A masking is applied to the region with thickness less than 1E-5.
  !
  !> \param[in]   output_idx      output index
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 15/12/2016
  !
  !******************************************************************************

  SUBROUTINE output_esri(output_idx)

    USE geometry_2d, ONLY : B_cent , grid_output , deposit , erosion , B_nodata
    USE geometry_2d, ONLY : deposit_tot , erosion_tot , B_zone
    ! USE geometry_2d, ONLY : comp_interfaces_x , comp_interfaces_y
    USE solver_2d, ONLY : qp

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: output_idx

    CHARACTER(LEN=4) :: idx_string
    CHARACTER(LEN=4) :: isolid_string
    INTEGER :: j
    INTEGER :: i_solid

    IF ( output_idx .EQ. 1 ) THEN

       OPEN(dem_esri_unit,FILE='dem_esri.asc',status='unknown',form='formatted')

       WRITE(dem_esri_unit,'(A,I5)') 'ncols ', comp_cells_x
       WRITE(dem_esri_unit,'(A,I5)') 'nrows ', comp_cells_y
       WRITE(dem_esri_unit,'(A,F15.3)') 'xllcorner ', x0
       WRITE(dem_esri_unit,'(A,F15.3)') 'yllcorner ', y0
       WRITE(dem_esri_unit,'(A,F15.3)') 'cellsize ', cell_size
       WRITE(dem_esri_unit,'(A,I5)') 'NODATA_value ', -9999

       DO j = comp_cells_y,1,-1

          WRITE(dem_esri_unit,*) B_cent(1:comp_cells_x,j)

       ENDDO

       CLOSE(dem_esri_unit)

       OPEN(dem_esri_unit,FILE='dem_esri_nodata.asc',status='unknown',form='formatted')

       WRITE(dem_esri_unit,'(A,I5)') 'ncols ', comp_cells_x
       WRITE(dem_esri_unit,'(A,I5)') 'nrows ', comp_cells_y
       WRITE(dem_esri_unit,'(A,F15.3)') 'xllcorner ', x0
       WRITE(dem_esri_unit,'(A,F15.3)') 'yllcorner ', y0
       WRITE(dem_esri_unit,'(A,F15.3)') 'cellsize ', cell_size
       WRITE(dem_esri_unit,'(A,I5)') 'NODATA_value ', -9999

       DO j = comp_cells_y,1,-1

          WRITE(dem_esri_unit,*) MERGE(1,0,B_nodata(1:comp_cells_x,j))

       ENDDO

       CLOSE(dem_esri_unit)

       IF ( liquid_vaporization_flag ) THEN

          OPEN(dem_esri_unit,FILE='dem_esri_water.asc',status='unknown',form='formatted')

          IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'WRITING dem_esri_water.asc'

          WRITE(dem_esri_unit,'(A,I5)') 'ncols ', comp_cells_x
          WRITE(dem_esri_unit,'(A,I5)') 'nrows ', comp_cells_y
          WRITE(dem_esri_unit,'(A,F15.3)') 'xllcorner ', x0
          WRITE(dem_esri_unit,'(A,F15.3)') 'yllcorner ', y0
          WRITE(dem_esri_unit,'(A,F15.3)') 'cellsize ', cell_size
          WRITE(dem_esri_unit,'(A,I5)') 'NODATA_value ', -9999

          DO j = comp_cells_y,1,-1

             WRITE(dem_esri_unit,*) B_zone(1:comp_cells_x,j)

          ENDDO

          CLOSE(dem_esri_unit)

       END IF

    END IF

    idx_string = lettera(output_idx-1)

    !Save thickness
    output_esri_file = TRIM(run_name)//'_'//idx_string//'.asc'

    IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'WRITING ',output_esri_file

    OPEN(output_esri_unit,FILE=output_esri_file,status='unknown',form='formatted')

    grid_output = -9999 

    WHERE ( qp(1,:,:).GE. 1.0E-5_wp )

       grid_output = qp(1,:,:) 

    END WHERE

    WRITE(output_esri_unit,'(A,I5)') 'ncols ', comp_cells_x
    WRITE(output_esri_unit,'(A,I5)') 'nrows ', comp_cells_y
    WRITE(output_esri_unit,'(A,F15.3)') 'xllcorner ', x0
    WRITE(output_esri_unit,'(A,F15.3)') 'yllcorner ', y0
    WRITE(output_esri_unit,'(A,F15.3)') 'cellsize ', cell_size
    WRITE(output_esri_unit,'(A,I5)') 'NODATA_value ', -9999

    DO j = comp_cells_y,1,-1

       WRITE(output_esri_unit,'(2000ES12.3E3)') grid_output(1:comp_cells_x,j)

    ENDDO

    CLOSE(output_esri_unit)

    !Save temperature
    output_esri_file = TRIM(run_name)//'_T_'//idx_string//'.asc'

    IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'WRITING ',output_esri_file

    OPEN(output_esri_unit,FILE=output_esri_file,status='unknown',form='formatted')

    grid_output = -9999 

    WHERE ( qp(1,:,:) .GE. 1.0E-5_wp )

       grid_output = qp(4,:,:)

    END WHERE

    WRITE(output_esri_unit,'(A,I5)') 'ncols ', comp_cells_x
    WRITE(output_esri_unit,'(A,I5)') 'nrows ', comp_cells_y
    WRITE(output_esri_unit,'(A,F15.3)') 'xllcorner ', x0
    WRITE(output_esri_unit,'(A,F15.3)') 'yllcorner ', y0
    WRITE(output_esri_unit,'(A,F15.3)') 'cellsize ', cell_size
    WRITE(output_esri_unit,'(A,I5)') 'NODATA_value ', -9999

    DO j = comp_cells_y,1,-1

       WRITE(output_esri_unit,'(2000ES12.3E3)') grid_output(1:comp_cells_x,j)

    ENDDO

    CLOSE(output_esri_unit)

    IF ( settling_flag ) THEN

       !Save solid classes deposit
       DO i_solid=1,n_solid

          isolid_string = lettera(i_solid)

          output_esri_file =                                                    &
               TRIM(run_name)//'_dep_'//isolid_string//'_'//idx_string//'.asc'

          IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'WRITING ',output_esri_file

          OPEN(output_esri_unit,FILE=output_esri_file,status='unknown',         &
               form='formatted')

          grid_output = -9999 

          WHERE ( deposit(:,:,i_solid) .GT. 0.0_wp )

             grid_output = deposit(:,:,i_solid)

          END WHERE

          WRITE(output_esri_unit,'(A,I5)') 'ncols ', comp_cells_x
          WRITE(output_esri_unit,'(A,I5)') 'nrows ', comp_cells_y
          WRITE(output_esri_unit,'(A,F15.3)') 'xllcorner ', x0
          WRITE(output_esri_unit,'(A,F15.3)') 'yllcorner ', y0
          WRITE(output_esri_unit,'(A,F15.3)') 'cellsize ', cell_size
          WRITE(output_esri_unit,'(A,I5)') 'NODATA_value ', -9999

          DO j = comp_cells_y,1,-1

             WRITE(output_esri_unit,'(2000ES12.3E3)')                           &
                  grid_output(1:comp_cells_x,j)

          ENDDO

          CLOSE(output_esri_unit)

       END DO

       ! Save total deposit (solid+continous phase associated with porosity)
       output_esri_file =                                                       &
            TRIM(run_name)//'_depTot_'//idx_string//'.asc'

       IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'WRITING ',output_esri_file

       OPEN(output_esri_unit,FILE=output_esri_file,status='unknown',            &
            form='formatted')

       grid_output = -9999 

       DO j = 1,comp_cells_y

          deposit_tot(:,j) = SUM(deposit(:,j,:),DIM=2) / ( 1.0_wp -             &
               erodible_porosity )  

       END DO

       WHERE ( deposit_tot(:,:) .GT. 0.0_wp )

          grid_output = deposit_tot

       END WHERE

       WRITE(output_esri_unit,'(A,I5)') 'ncols ', comp_cells_x
       WRITE(output_esri_unit,'(A,I5)') 'nrows ', comp_cells_y
       WRITE(output_esri_unit,'(A,F15.3)') 'xllcorner ', x0
       WRITE(output_esri_unit,'(A,F15.3)') 'yllcorner ', y0
       WRITE(output_esri_unit,'(A,F15.3)') 'cellsize ', cell_size
       WRITE(output_esri_unit,'(A,I5)') 'NODATA_value ', -9999

       DO j = comp_cells_y,1,-1

          WRITE(output_esri_unit,'(2000ES12.3E3)')                              &
               grid_output(1:comp_cells_x,j)

       ENDDO

       CLOSE(output_esri_unit)

    END IF

    IF ( erosion_coeff .GT. 0.0_wp ) THEN

       !Save erosion
       DO i_solid=1,n_solid

          isolid_string = lettera(i_solid)

          output_esri_file =                                                    &
               TRIM(run_name)//'_ers_'//isolid_string//'_'//idx_string//'.asc'

          IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'WRITING ',output_esri_file

          OPEN(output_esri_unit,FILE=output_esri_file,status='unknown',         &
               form='formatted')

          grid_output = -9999 

          WHERE ( erosion(:,:,i_solid) .GT. 0.0_wp )

             grid_output = erosion(:,:,i_solid)

          END WHERE

          WRITE(output_esri_unit,'(A,I5)') 'ncols ', comp_cells_x
          WRITE(output_esri_unit,'(A,I5)') 'nrows ', comp_cells_y
          WRITE(output_esri_unit,'(A,F15.3)') 'xllcorner ', x0
          WRITE(output_esri_unit,'(A,F15.3)') 'yllcorner ', y0
          WRITE(output_esri_unit,'(A,F15.3)') 'cellsize ', cell_size
          WRITE(output_esri_unit,'(A,I5)') 'NODATA_value ', -9999

          DO j = comp_cells_y,1,-1

             WRITE(output_esri_unit,'(2000ES12.3E3)')                           &
                  grid_output(1:comp_cells_x,j)

          ENDDO

          CLOSE(output_esri_unit)

       END DO

       ! Save total erosion (solid+continous phase associated with porosity)
       output_esri_file =                                                       &
            TRIM(run_name)//'_ersTot_'//idx_string//'.asc'

       IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'WRITING ',output_esri_file

       OPEN(output_esri_unit,FILE=output_esri_file,status='unknown',            &
            form='formatted')

       grid_output = -9999 

       DO j = 1,comp_cells_y

          erosion_tot(:,j) = SUM(erosion(:,j,:),DIM=2) / (1.0_wp -              &
               erodible_porosity)  

       END DO


       WHERE ( erosion_tot(:,:) .GT. 0.0_wp )

          grid_output = erosion_tot

       END WHERE

       WRITE(output_esri_unit,'(A,I5)') 'ncols ', comp_cells_x
       WRITE(output_esri_unit,'(A,I5)') 'nrows ', comp_cells_y
       WRITE(output_esri_unit,'(A,F15.3)') 'xllcorner ', x0
       WRITE(output_esri_unit,'(A,F15.3)') 'yllcorner ', y0
       WRITE(output_esri_unit,'(A,F15.3)') 'cellsize ', cell_size
       WRITE(output_esri_unit,'(A,I5)') 'NODATA_value ', -9999

       DO j = comp_cells_y,1,-1

          WRITE(output_esri_unit,'(2000ES12.3E3)') grid_output(1:comp_cells_x,j)

       ENDDO

       CLOSE(output_esri_unit)

    END IF

    RETURN

  END SUBROUTINE output_esri

  SUBROUTINE close_units

    IMPLICIT NONE

    IF ( output_runout_flag) THEN

       CLOSE(runout_unit)
       CLOSE(mass_center_unit)

    END IF

  END SUBROUTINE close_units

  !******************************************************************************
  !> \brief Numeric to String conversion
  !
  !> This function convert the integer in input into a numeric string for
  !> the subfix of the output files.
  !
  !> \param[in]   k      integer to convert         
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  CHARACTER(LEN=4) FUNCTION lettera(k)
    IMPLICIT NONE
    CHARACTER ones,tens,hund,thou
    !
    INTEGER :: k
    !
    INTEGER :: iten, ione, ihund, ithou
    !
    ithou=INT(k/1000)
    ihund=INT((k-(ithou*1000))/100)
    iten=INT((k-(ithou*1000)-(ihund*100))/10)
    ione=k-ithou*1000-ihund*100-iten*10
    ones=CHAR(ione+48)
    tens=CHAR(iten+48)
    hund=CHAR(ihunD+48)
    thou=CHAR(ithou+48)
    lettera=thou//hund//tens//ones
    !
    RETURN
  END FUNCTION lettera

  !******************************************************************************
  !> \brief Write solution at selected points on file
  !
  !> This subroutine writes on a file the thickness at selected points, defined
  !> by an appropriate card in the input file.
  !> in the initial solution.
  !> \param[in]   output_idx      output index
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 12/02/2018
  !******************************************************************************

  SUBROUTINE output_probes(time)

    USE geometry_2d, ONLY : x_comp , y_comp , deposit
    USE parameters_2d, ONLY : t_probes , n_vars
    USE solver_2d, ONLY : q , qp


    USE geometry_2d, ONLY : interp_2d_scalarB


    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: time

    CHARACTER(LEN=4) :: idx_string

    REAL(wp) :: h_prb
    REAL(wp) :: hrhom_prb
    REAL(wp) :: rhom_prb
    REAL(wp) :: u_prb
    REAL(wp) :: v_prb
    REAL(wp) :: T_prb
    REAL(wp) :: pDyn_prb

    REAL(wp) :: alphas_prb(n_solid)
    REAL(wp) :: alphag_prb(n_add_gas)
    REAL(wp) :: dep_probe(n_solid)

    INTEGER :: k

    INTEGER :: i_solid , i_gas

    DO k=1,n_probes

       idx_string = lettera(k)

       probes_file = TRIM(run_name)//'_prb_'//idx_string//'.csv'

       IF ( time .EQ. t_start ) THEN

          OPEN(probes_unit,FILE=probes_file,status='unknown',form='formatted')
          WRITE(probes_unit,'(F14.3, A, F14.3)') probes_coords(1,k) ,',', probes_coords(2,k)
          WRITE(probes_unit,'(A15)',ADVANCE='no') 'time,'
          WRITE(probes_unit,'(A15)',ADVANCE='no') 'h,'
          WRITE(probes_unit,'(A15)',ADVANCE='no') 'rho,'
          WRITE(probes_unit,'(A15)',ADVANCE='no') 'T,'
          WRITE(probes_unit,'(A15)',ADVANCE='no') 'u,'
          WRITE(probes_unit,'(A15)',ADVANCE='no') 'v,'


          DO i_solid=1,n_solid

             WRITE(idx_string,'(I2.2)') i_solid
             WRITE(probes_unit,'(A15)',ADVANCE='no') 'alfa_s('//TRIM(idx_string)//'),'

          END DO

          DO i_gas=1,n_add_gas

             WRITE(idx_string,'(I2.2)') i_gas
             WRITE(probes_unit,'(A15)',ADVANCE='no') 'alfa_g('//TRIM(idx_string)//'),'

          END DO

          IF ( settling_flag ) THEN

             DO i_solid=1,n_solid

                WRITE(idx_string,'(I2.2)') i_solid
                WRITE(probes_unit,'(A15)',ADVANCE='no') 'dep_s('//TRIM(idx_string)//'),'

             END DO

          END IF

          WRITE(probes_unit,'(A15)') 'Dyn.pres.'

          ! Write a line with units
          WRITE(probes_unit,'(A15)',ADVANCE='no') '(s),'
          WRITE(probes_unit,'(A15)',ADVANCE='no') '(m),'
          WRITE(probes_unit,'(A15)',ADVANCE='no') '(kg/m3),'
          WRITE(probes_unit,'(A15)',ADVANCE='no') '(K),'
          WRITE(probes_unit,'(A15)',ADVANCE='no') '(m/s),'
          WRITE(probes_unit,'(A15)',ADVANCE='no') '(m/s),'


          DO i_solid=1,n_solid

             WRITE(probes_unit,'(A15)',ADVANCE='no') '(),'

          END DO

          DO i_gas=1,n_add_gas

             WRITE(probes_unit,'(A15)',ADVANCE='no') '(),'

          END DO

          IF ( settling_flag ) THEN

             DO i_solid=1,n_solid

                WRITE(probes_unit,'(A15)',ADVANCE='no') '(m),'

             END DO

          END IF

          WRITE(probes_unit,'(A15)') '(Pa)'


       ELSE

          OPEN(probes_unit,FILE=probes_file,status='old',position='append',     &
               form='formatted')

       END IF

       WRITE(probes_unit,1710,ADVANCE='no') time,','

       CALL interp_2d_scalarB( x_comp , y_comp , qp(1,:,:)  ,                   &
            probes_coords(1,k) , probes_coords(2,k) , h_prb )

       WRITE(probes_unit,1710,ADVANCE='no') h_prb,',' 

       IF ( h_prb .GT. 1.0E-5_wp ) THEN

          CALL interp_2d_scalarB( x_comp , y_comp , q(1,:,:) ,                  &
               probes_coords(1,k) , probes_coords(2,k) , hrhom_prb )

          rhom_prb = hrhom_prb/h_prb

          WRITE(probes_unit,1710,ADVANCE='no') rhom_prb,','

          CALL interp_2d_scalarB( x_comp , y_comp , qp(4,:,:)  ,                &
               probes_coords(1,k) , probes_coords(2,k) , T_prb )

          WRITE(probes_unit,1710,ADVANCE='no') T_prb,','

          CALL interp_2d_scalarB( x_comp , y_comp , qp(n_vars+1,:,:)  ,         &
               probes_coords(1,k) , probes_coords(2,k) , u_prb )

          WRITE(probes_unit,1710,ADVANCE='no') u_prb,','

          CALL interp_2d_scalarB( x_comp , y_comp , qp(n_vars+2,:,:)  ,         &
               probes_coords(1,k) , probes_coords(2,k) , v_prb )

          WRITE(probes_unit,1710,ADVANCE='no') v_prb,','

          DO i_solid=1,n_solid

             CALL interp_2d_scalarB( x_comp , y_comp ,  qp(4+i_solid,:,:)  ,    &
                  probes_coords(1,k) , probes_coords(2,k) , alphas_prb(i_solid) )

             IF ( alpha_flag ) THEN

                WRITE(probes_unit,1710,ADVANCE='no') alphas_prb(i_solid),','

             ELSE

                WRITE(probes_unit,1710,ADVANCE='no') alphas_prb(i_solid)/h_prb,','

             END IF

          END DO

          DO i_gas=1,n_add_gas

             CALL interp_2d_scalarB( x_comp , y_comp , qp(4+n_solid+i_gas,:,:)  ,   &
                  probes_coords(1,k) , probes_coords(2,k) , alphag_prb(i_gas) )

             IF ( alpha_flag ) THEN

                WRITE(probes_unit,1710,ADVANCE='no') alphag_prb(i_gas),','

             ELSE

                WRITE(probes_unit,1710,ADVANCE='no') alphag_prb(i_gas)/h_prb,','

             END IF

          END DO

          pDyn_prb = 0.5 * rhom_prb * ( u_prb**2 + v_prb**2 ) 

       ELSE

          WRITE(probes_unit,1710,ADVANCE='no') 0.0_wp,','

          WRITE(probes_unit,1710,ADVANCE='no') T_ambient,','

          WRITE(probes_unit,1710,ADVANCE='no') 0.0_wp,','

          WRITE(probes_unit,1710,ADVANCE='no') 0.0_wp,','

          DO i_solid=1,n_solid

             WRITE(probes_unit,1710,ADVANCE='no') 0.0_wp,','

          END DO

          DO i_gas=1,n_add_gas

             WRITE(probes_unit,1710,ADVANCE='no') 0.0_wp,','

          END DO

          pDyn_prb = 0.0_wp

       END IF

       IF ( settling_flag ) THEN

          DO i_solid=1,n_solid

             CALL interp_2d_scalarB( x_comp , y_comp , DEPOSIT(:,:,i_solid)  ,   &
                  probes_coords(1,k) , probes_coords(2,k) , dep_probe(i_solid) )

             WRITE(probes_unit,1710,ADVANCE='no') dep_probe(i_solid),','

          END DO

       END IF

       WRITE(probes_unit,1710) pDyn_prb

       CLOSE(probes_unit)

    END DO

1710 FORMAT(ES14.7E2,A1)

    t_probes = time + dt_probes

  END SUBROUTINE output_probes

  !******************************************************************************
  !> \brief Write runout on file
  !
  !> This subroutine writes on a file the flow runout. It is calculated as the 
  !> linear horizontal distance from the point with the highest topography value
  !> in the initial solution.
  !> \param[in]     time             actual time
  !> \param[inout]  stop_flag        logical to check if flow has stopped
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 12/02/2018
  !******************************************************************************

  SUBROUTINE output_runout(time,stop_flag)

    USE geometry_2d, ONLY : x_comp , y_comp , B_cent , dx , dy
    USE parameters_2d, ONLY : t_runout 
    USE solver_2d, ONLY : qp , q , hpos , hpos_old


    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: time
    LOGICAL, INTENT(INOUT) :: stop_flag

    REAL(wp), ALLOCATABLE :: X(:,:), Y(:,:) 
    REAL(wp), ALLOCATABLE :: dist(:,:) , dist_x(:,:) , dist_y(:,:)
    INTEGER :: sX, sY
    INTEGER :: imax(2) , imax_x(2) , imax_y(2) , imin(2)

    INTEGER :: j,k

    REAL(wp) :: area , area_old , area_new_rel

    REAL(wp) :: x_mass_center , y_mass_center

    REAL(wp) :: vel_mass_center , vel_radial_growth

    CHARACTER(18) :: txt_string
    REAL(wp) :: old_runout

    sX = size(x_comp) 
    sY = size(y_comp) 

    ALLOCATE( X(sX,sY) , Y(sX,sY) , dist(sX,sY), dist_x(sX,sY), dist_y(sX,sY) )

    ! This work with large 
    !X(:,:) = SPREAD( x_comp, 2, sY )
    !Y(:,:) = SPREAD( y_comp, 1, sX )

    !$OMP PARALLEL 
    !$OMP DO 
    DO k=1,sY

       X(1:sX,k) = x_comp(1:sX)

    END DO
    !$OMP END DO

    !$OMP DO
    DO j=1,sX

       Y(j,1:sY) = y_comp(1:sY)

    END DO
    !$OMP END DO
    !$OMP END PARALLEL

    dist(:,:) = 0.0_wp

    IF ( time .EQ. t_start ) THEN

       WRITE(runout_unit,'(A13,A13,A13)') 'Time (s),','Runout (m),',' Area (m2)'

       WRITE(mass_center_unit,'(A18,A18,A18,A18,A18)')  'Time (s),',            &
            'Barycenter-x (m),' , 'Barycenter-y (m),' , 'Baricenter-dx (m),' ,  &
            'Baricenter-dy (m)'

       CALL flush(runout_unit)

       IF ( MAXVAL( qp(1,:,:) ) .EQ. 0.0_wp ) THEN

          IF ( collapsing_volume_flag ) THEN

             x_mass_center = x_collapse
             y_mass_center = y_collapse

          END IF

          IF ( radial_source_flag .OR. bottom_radial_source_flag ) THEN

             x_mass_center = x_source
             y_mass_center = y_source

          END IF

       ELSE

          x_mass_center = SUM( X*q(1,:,:) ) / SUM( q(1,:,:) )
          y_mass_center = SUM( Y*q(1,:,:) ) / SUM( q(1,:,:) )
          hpos = ( qp(1,:,:) .GT. 1.0E-5_wp )

       END IF

       hpos_old = ( qp(1,:,:) .GT. 1.0E-5_wp )

       x_mass_center_old = x_mass_center
       y_mass_center_old = y_mass_center

       IF ( ( x0_runout .EQ. -1 ) .AND. ( y0_runout .EQ. -1 ) ) THEN

          WHERE( qp(1,:,:) > 1.0E-5_wp ) dist = B_cent
          imin = MAXLOC( dist )

          x0_runout = X(imin(1),imin(2))
          y0_runout = Y(imin(1),imin(2))

          WRITE(*,*) 'Runout calculated as linear distance from: (' ,           &
               x0_runout ,',',y0_runout,')'

          dist(:,:) = 0.0_wp

          WHERE( hpos ) dist = SQRT( (X-x0_runout)**2 + ( Y - y0_runout )**2 )

          imax = MAXLOC( dist )

          init_runout = dist(imax(1),imax(2))

          dist_x(:,:) = 0.0_wp

          WHERE( hpos ) dist_x = SQRT( (X-x0_runout)**2 )

          imax_x = MAXLOC( dist_x )

          init_runout_x = dist(imax_x(1),imax_x(2))

          dist_y(:,:) = 0.0_wp

          WHERE( hpos ) dist_y = SQRT( (Y-y0_runout)**2 )

          imax_y = MAXLOC( dist_y )

          init_runout_y = dist(imax_y(1),imax_y(2))

       ELSE

          init_runout = 0.0_wp
          init_runout_x = 0.0_wp
          init_runout_y = 0.0_wp

       END IF

    ELSE

       IF ( MAXVAL( qp(1,:,:) ) .EQ. 0.0_wp ) THEN

          x_mass_center = x_mass_center_old
          y_mass_center = y_mass_center_old

       ELSE

          x_mass_center = SUM( X*q(1,:,:) ) / SUM( q(1,:,:) )
          y_mass_center = SUM( Y*q(1,:,:) ) / SUM( q(1,:,:) )

       END IF

       hpos = ( qp(1,:,:) .GT. 1.0E-5_wp )

    END IF

    dist(:,:) = 0.0_wp

    WHERE( hpos ) dist = SQRT( ( X - x0_runout )**2 + ( Y - y0_runout )**2 )

    imax = MAXLOC( dist )

    dist_x(:,:) = 0.0_wp

    WHERE( hpos ) dist_x = SQRT( ( X - x0_runout )**2 )

    imax_x = MAXLOC( dist_x )

    dist_y(:,:) = 0.0_wp

    WHERE( hpos ) dist_y = SQRT( ( Y - y0_runout )**2 )

    imax_y = MAXLOC( dist_y )

    IF ( time .GT. t_start ) THEN

       OPEN(dakota_unit,FILE='runout.txt',status='old',form='formatted')
       READ(dakota_unit,'(A18,F12.3)') txt_string,old_runout
       CLOSE(dakota_unit)

    ELSE

       old_runout = 0.0_wp

    END IF

    OPEN(dakota_unit,FILE='runout.txt',status='replace',form='formatted')
    WRITE(dakota_unit,'(A18,F12.3)') 'maximum runout =', MAX(old_runout,dist(imax(1),imax(2)) - init_runout)
    WRITE(dakota_unit,'(A18,F12.3)') 'final runout =', dist(imax(1),imax(2)) - init_runout

    CLOSE(dakota_unit)

    area_old = dx*dy*COUNT(hpos_old)
    area = dx*dy*COUNT(hpos)

    WRITE(runout_unit,'(F12.3,A,F12.3,A,F14.3)') time ,',',                     &
         dist(imax(1),imax(2)) - init_runout,',',area

    CALL flush(runout_unit)

    WRITE(mass_center_unit,'(F17.3,A,F17.3,A,F17.3,A,F17.3,A,F17.3)') time ,    &
         ',' , x_mass_center , ',' , y_mass_center , ',' ,                      &
         dist_x(imax_x(1),imax_x(2)) - init_runout_x , ',' ,                    &
         dist_y(imax_y(1),imax_y(2)) - init_runout_y  

    CALL flush(mass_center_unit)

    IF ( time .GT. t_start ) THEN

       vel_mass_center = SQRT( ( x_mass_center_old - x_mass_center )**2 +       &
            ( y_mass_center_old - y_mass_center )**2 ) / dt_runout

       vel_radial_growth = ABS( SQRT( area ) - SQRT( area_old ) ) / dt_runout

       area_new_rel = dx*dy*COUNT( hpos .AND. ( .NOT.hpos_old ) ) / COUNT( hpos )

       x_mass_center_old = x_mass_center
       y_mass_center_old = y_mass_center
       hpos_old = hpos

       IF ( ( MAX( vel_mass_center , area_new_rel /dt_runout ) .LT. eps_stop )  &
            .AND. (.NOT.stop_flag) ) THEN

          WRITE(*,*) 'Steady solution reached'
          WRITE(*,*) 'vel_mass_center',vel_mass_center
          WRITE(*,*) 'vel_radial_growth',vel_radial_growth
          WRITE(*,*) 'area_new_rel/dt_runout',area , area_new_rel/dt_runout
          stop_flag = .TRUE.

       END IF

    END IF

    DEALLOCATE( X , Y , dist , dist_x , dist_y )

    t_runout = time + dt_runout

  END SUBROUTINE output_runout

END MODULE inpout_2d

