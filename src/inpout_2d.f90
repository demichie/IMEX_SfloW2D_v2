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
  USE parameters_2d, ONLY : t_start , t_end , dt_output 

  USE solver_2d, ONLY : verbose_level

  ! -- Variables for the namelist NEWRUN_PARAMETERS
  USE geometry_2d, ONLY : x0 , y0 , comp_cells_x , comp_cells_y , cell_size
  USE geometry_2d, ONLY : topography_profile , n_topography_profile_x ,         &
       n_topography_profile_y
  USE init_2d, ONLY : riemann_interface
  USE parameters_2d, ONLY : n_solid
  USE parameters_2d, ONLY : riemann_flag , rheology_flag , energy_flag ,        &
       topo_change_flag , radial_source_flag , collapsing_volume_flag ,         &
       liquid_flag , gas_flag , subtract_init_flag

  ! -- Variables for the namelist INITIAL_CONDITIONS
  USE parameters_2d, ONLY : released_volume , x_release , y_release
  USE parameters_2d, ONLY : velocity_mod_release , velocity_ang_release
  USE parameters_2d, ONLY : alphas_init
  USE parameters_2d, ONLY : T_init

  ! -- Variables for the namelist LEFT_STATE
  USE init_2d, ONLY : hB_W , u_W , v_W , alphas_W , T_W

  ! -- Variables for the namelist RIGHT_STATE
  USE init_2d, ONLY : hB_E , u_E , v_E , alphas_E , T_E

  ! -- Variables for the namelists LEFT/RIGHT_BOUNDARY_CONDITIONS
  USE parameters_2d, ONLY : bc

  ! -- Variables for the namelist NUMERIC_PARAMETERS
  USE parameters_2d, ONLY : solver_scheme, dt0 , max_dt , cfl, limiter , theta, &
       reconstr_coeff , interfaces_relaxation , n_RK   

  ! -- Variables for the namelist EXPL_TERMS_PARAMETERS
  USE constitutive_2d, ONLY : grav , inv_grav

  ! -- Variables for the namelist RADIAL_SOURCE_PARAMETERS
  USE parameters_2d, ONLY : x_source , y_source , r_source , vel_source ,       &
       T_source , h_source , alphas_source , alphal_source , time_param

  ! -- Variables for the namelist COLLAPSING_VOLUME_PARAMETERS
  USE parameters_2d, ONLY : x_collapse , y_collapse , r_collapse , T_collapse , &
       h_collapse , alphas_collapse
  
  ! -- Variables for the namelist TEMPERATURE_PARAMETERS
  USE constitutive_2d, ONLY : emissivity , exp_area_fract , enne , emme ,       &
       atm_heat_transf_coeff , thermal_conductivity , T_env , T_ground , c_p
    
  ! -- Variables for the namelist RHEOLOGY_PARAMETERS
  USE parameters_2d, ONLY : rheology_model
  USE constitutive_2d, ONLY : mu , xi , tau , nu_ref , visc_par , T_ref
  USE constitutive_2d, ONLY : alpha2 , beta2 , alpha1_coeff , beta1 , Kappa ,n_td
  USE constitutive_2d, ONLY : friction_factor
  
  ! --- Variables for the namelist SOLID_TRANSPORT_PARAMETERS
  USE constitutive_2d, ONLY : rho_s , diam_s , sp_heat_s
  USE constitutive_2d, ONLY : settling_flag , erosion_coeff , erodible_porosity
  USE constitutive_2d, ONLY : erodible_fract , T_erodible

  ! --- Variables for the namelist GAS_TRANSPORT_PARAMETERS
  USE constitutive_2d, ONLY : sp_heat_a , sp_gas_const_a , kin_visc_a , pres ,  &
       T_ambient , entrainment_flag

  ! --- Variables for the namelist LIQUID_TRANSPORT_PARAMETERS
  USE constitutive_2d, ONLY : sp_heat_l , rho_l , kin_visc_l

  ! --- Variables for the namelist VULNERABILITY_TABLE_PARAMETERS
  USE parameters_2d, ONLY : n_thickness_levels , n_dyn_pres_levels ,            &
       thickness_levels , dyn_pres_levels
  

  IMPLICIT NONE

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
  ! TYPE(bc) :: alphas_bcW(10)
  TYPE(bc),ALLOCATABLE :: alphas_bcW(:)

  ! -- Variables for the namelists EAST_BOUNDARY_CONDITIONS
  TYPE(bc) :: h_bcE , hu_bcE , hv_bcE , T_bcE
  !TYPE(bc):: alphas_bcE(10)
  TYPE(bc),ALLOCATABLE :: alphas_bcE(:)

  ! -- Variables for the namelists SOUTH_BOUNDARY_CONDITIONS
  TYPE(bc) :: h_bcS , hu_bcS , hv_bcS , T_bcS
  !TYPE(bc) :: alphas_bcS(10)
  TYPE(bc),ALLOCATABLE :: alphas_bcS(:)

  ! -- Variables for the namelists NORTH_BOUNDARY_CONDITIONS
  TYPE(bc) :: h_bcN , hu_bcN , hv_bcN , T_bcN
  !TYPE(bc) :: alphas_bcN(10)
  TYPE(bc),ALLOCATABLE :: alphas_bcN(:)


  ! parameters to read a dem file
  INTEGER :: ncols, nrows, nodata_value

  REAL(wp) :: xllcorner, yllcorner, cellsize

  LOGICAL :: write_first_q

  INTEGER :: n_probes

  REAL(wp), ALLOCATABLE :: probes_coords(:,:)

  REAL(wp) :: dt_runout

  REAL(wp), ALLOCATABLE :: h_old(:,:)

  REAL(wp) :: x0_runout, y0_runout , init_runout , eps_stop

  ! absolute precentages of solids in the initial volume
  REAL(wp), ALLOCATABLE :: sed_vol_perc(:)

  REAL(wp) :: alphas0_E(10) , alphas0_W(10)
  
  REAL(wp) :: alpha1_ref

  REAL(wp) :: thickness_levels0(10)
  REAL(wp) :: dyn_pres_levels0(10)
  
  NAMELIST / run_parameters / run_name , restart , t_start , t_end , dt_output ,&
       output_cons_flag , output_esri_flag , output_phys_flag ,                 &
       output_runout_flag , verbose_level
  NAMELIST / restart_parameters / restart_file, T_init, T_ambient , sed_vol_perc

  NAMELIST / newrun_parameters / n_solid , topography_file , x0 , y0 ,          &
       comp_cells_x , comp_cells_y , cell_size , rheology_flag , riemann_flag , &
       energy_flag , liquid_flag , radial_source_flag , collapsing_volume_flag ,&
       topo_change_flag , gas_flag , subtract_init_flag

  NAMELIST / initial_conditions /  released_volume , x_release , y_release ,    &
       velocity_mod_release , velocity_ang_release , T_init , T_ambient

  NAMELIST / left_state / riemann_interface , hB_W , u_W , v_W , alphas0_W , T_W

  NAMELIST / right_state / hB_E , u_E , v_E , alphas0_E , T_E

  NAMELIST / numeric_parameters / solver_scheme, dt0 , max_dt , cfl, limiter ,  &
       theta , reconstr_coeff , interfaces_relaxation , n_RK   

  NAMELIST / expl_terms_parameters / grav
 
  NAMELIST / radial_source_parameters / x_source , y_source , r_source ,        &
       vel_source , T_source , h_source , alphas_source , alphal_source ,       &
       time_param

  NAMELIST / collapsing_volume_parameters / x_collapse , y_collapse ,           &
       r_collapse , T_collapse , h_collapse , alphas_collapse
 
  NAMELIST / temperature_parameters / emissivity ,  atm_heat_transf_coeff ,     &
       thermal_conductivity , exp_area_fract , c_p , enne , emme , T_env ,      &
       T_ground

  NAMELIST / rheology_parameters / rheology_model , mu , xi , tau , nu_ref ,    &
       visc_par , T_ref , alpha2 , beta2 , alpha1_ref , beta1 , Kappa , n_td ,  &
       friction_factor

  NAMELIST / runout_parameters / x0_runout , y0_runout , dt_runout ,            &
       eps_stop

  NAMELIST / gas_transport_parameters / sp_heat_a , sp_gas_const_a , kin_visc_a,&
       pres , T_ambient , entrainment_flag

  NAMELIST / liquid_transport_parameters / sp_heat_l , rho_l , kin_visc_l

  NAMELIST / vulnerability_table_parameters / thickness_levels0 , dyn_pres_levels0
  
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

    !-- Inizialization of the Variables for the namelist restart parameters
    restart_file = ''
    T_init = 0.0_wp
    T_ambient = 0.0_wp

    !-- Inizialization of the Variables for the namelist newrun_parameters
    topography_file = 'topography_dem.asc'
    x0 = 0.0_wp
    y0 = 0.0_wp
    comp_cells_x = 1000
    comp_cells_y = 1
    cell_size = 1.0E-3_wp
    rheology_flag = .FALSE.
    riemann_flag =.TRUE.
    energy_flag = .FALSE.
    topo_change_flag = .FALSE.
    radial_source_flag = .FALSE.
    collapsing_volume_flag = .FALSE.
    liquid_flag = .FALSE.
    gas_flag = .TRUE.
    subtract_init_flag = .FALSE.

    !-- Inizialization of the Variables for the namelist left_state
    riemann_interface = 0.5_wp
    hB_W = 2.0_wp
    u_W = 0.0_wp
    v_W = 0.0_wp
    ! alphas_W = 0.5_wp
    T_W = -1.0_wp

    !-- Inizialization of the Variables for the namelist right_state
    hB_E = 1.0_wp
    u_E = 0.0_wp
    u_E = 0.0_wp
    ! alphas_E = 0.5_wp
    T_E = -1.0_wp

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
    mu = 0.0_wp
    xi = 0.0_wp
    tau = 0.0_wp
    T_ref = 0.0_wp
    visc_par = 0.0_wp

    !-- Inizialization of the Variables for the namelist RUNOUT_PARAMETERS
    x0_runout = -1
    y0_runout = -1
    dt_runout = 60
    eps_stop = 0.0_wp

    !-------------- Check if input file exists ----------------------------------
    input_file = 'SW_VAR_DENS_MODEL.inp'

    INQUIRE (FILE=input_file,exist=lexist)

    IF (lexist .EQV. .FALSE.) THEN

       WRITE(*,*) 'Input file SW_VAR_DENS_MODEL.inp not found'
       STOP

    ELSE

       OPEN(input_unit,FILE=input_file,STATUS='old')

       ! ------- READ run_parameters NAMELIST -----------------------------------
       READ(input_unit, newrun_parameters,IOSTAT=ios )

       IF ( ios .NE. 0 ) THEN

          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist NEWRUN_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          STOP

       ELSE

          REWIND(input_unit)

       END IF

       IF ( n_solid .LT. 1 ) THEN

          WRITE(*,*) 'ERROR: problem with namelist SOLID_TRANSPORT_PARAMETERS'
          WRITE(*,*) 'n_solid =' , n_solid
          WRITE(*,*) 'Please check the input file'
          STOP

       END IF


       ALLOCATE ( alphas_bcW(n_solid) )
       ALLOCATE ( alphas_bcE(n_solid) )
       ALLOCATE ( alphas_bcS(n_solid) )
       ALLOCATE ( alphas_bcN(n_solid) )

       ALLOCATE( sed_vol_perc(n_solid) )
       sed_vol_perc(1:n_solid) = -1.0_wp
       
       ALLOCATE( rho_s(n_solid) )
       ALLOCATE( diam_s(n_solid) )
       ALLOCATE( sp_heat_s(n_solid) )
       ALLOCATE( erodible_fract(n_solid) )

       CLOSE(input_unit)

    END IF

    ! output file index
    output_idx = 0

    ! -------------- Initialize values for checks during input reading ----------
    h_bcW%flag = -1 
    hu_bcW%flag = -1 
    hv_bcW%flag = -1 
    alphas_bcW%flag = -1 
    T_bcW%flag = -1 

    h_bcE%flag = -1 
    hu_bcE%flag = -1 
    hv_bcE%flag = -1 
    alphas_bcE%flag = -1 
    T_bcE%flag = -1 

    h_bcS%flag = -1 
    hu_bcS%flag = -1 
    hv_bcS%flag = -1 
    alphas_bcS%flag = -1 
    T_bcS%flag = -1 

    h_bcN%flag = -1 
    hu_bcN%flag = -1 
    hv_bcN%flag = -1 
    alphas_bcN%flag = -1 
    T_bcN%flag = -1 

    ! sed_vol_perc = -1.0_wp

    rheology_model = -1
    mu = -1
    xi = -1
    tau = -1
    nu_ref = -1
    visc_par = -1
    T_ref = -1
    friction_factor = -1

    alpha2 = -1 
    beta2 = -1
    alpha1_ref = -1
    beta1 = -1
    Kappa = -1
    n_td = -1
    rho_s = -1

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

    !- Variables for the namelist GAS_TRANSPORT_PARAMETERS
    sp_heat_a = -1.0_wp
    sp_gas_const_a = -1.0_wp
    kin_visc_a = -1.0_wp
    pres = -1.0_wp
    T_ambient = -1.0_wp

    !- Variables for the namelist LIQUID_TRANSPORT_PARAMETERS
    sp_heat_l = -1.0_wp
    rho_l = -1.0_wp
    kin_visc_l = -1.0_wp

    !- Variables for the namelist RADIAL_SOURCE_PARAMETERS
    T_source = -1.0_wp
    h_source = -1.0_wp
    r_source = -1.0_wp
    vel_source = -1.0_wp
    time_param(1:4) = -1.0_wp

    !- Variables for the namelist COLLAPSING_VOLUME_PARAMETERS
    T_collapse = -1.0_wp
    h_collapse = -1.0_wp
    r_collapse = -1.0_wp

    !- Variables for the namelist VULNERABILTY_TABLE_PARAMETERS
    n_thickness_levels = -1
    n_dyn_pres_levels = -1
    thickness_levels0 = -1.0_wp
    dyn_pres_levels0 = -1.0_wp

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

    USE geometry_2d, ONLY : deposit , erosion , erodible

    USE constitutive_2d, ONLY : rho_a_amb
    USE constitutive_2d, ONLY : rho_c_sub
    USE constitutive_2d, ONLY : kin_visc_c , rho_c , inv_rho_c , sp_heat_c ,    &
         r_rho_c , r_inv_rho_c
   
    USE constitutive_2d, ONLY : inv_pres , inv_rho_l , inv_rho_s , c_inv_rho_s
 
    USE constitutive_2d, ONLY : n_td2

    IMPLICIT none

    NAMELIST / west_boundary_conditions / h_bcW , hu_bcW , hv_bcW ,             &
         alphas_bcW , T_bcW

    NAMELIST / east_boundary_conditions / h_bcE , hu_bcE , hv_bcE ,             &
         alphas_bcE , T_bcE

    NAMELIST / south_boundary_conditions / h_bcS , hu_bcS , hv_bcS ,            &
         alphas_bcS , T_bcS

    NAMELIST / north_boundary_conditions / h_bcN , hu_bcN , hv_bcN ,            &
         alphas_bcN , T_bcN

    NAMELIST / solid_transport_parameters / rho_s , diam_s , sp_heat_s ,        &
         erosion_coeff , erodible_porosity , settling_flag , T_erodible ,       &
         erodible_file , erodible_fract

    
    REAL(wp) :: max_cfl

    LOGICAL :: tend1 
    CHARACTER(LEN=80) :: card

    INTEGER :: i_solid , j , k

    INTEGER :: dot_idx
    
    CHARACTER(LEN=3) :: check_file

    LOGICAL :: lexist

    CHARACTER(LEN=15) :: chara

    INTEGER :: ios
    
    REAL(wp) :: expA , expB , Tc

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

    IF ( n_solid .LT. 1 ) THEN
       
       WRITE(*,*) 'ERROR: problem with namelist NEWRUN_PARAMETERS'
       WRITE(*,*) 'n_solid =' , n_solid
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
    ! ------- READ liquid_transport_parameters NAMELIST -------------------------
    
    n_vars = 4
    
    IF ( liquid_flag ) THEN

       IF ( gas_flag ) n_vars = n_vars + 1

       READ(input_unit, liquid_transport_parameters,IOSTAT=ios)

       IF ( ios .NE. 0 ) THEN

          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist LIQUID_TRANSPORT_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          STOP

       ELSE

          REWIND(input_unit)

       END IF

       IF ( sp_heat_l .EQ. -1.0_wp ) THEN

          WRITE(*,*) 'ERROR: problem with namelist LIQUID_TRANSPORT_PARAMETERS'
          WRITE(*,*) 'SP_HEAT_l =' , sp_heat_l
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

       READ(input_unit, rheology_parameters,IOSTAT=ios)
       REWIND(input_unit)

       IF ( kin_visc_l .EQ. -1.0_wp ) THEN

          IF ( RHEOLOGY_MODEL .NE. 4 ) THEN

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

       END IF
       
       IF ( .NOT. gas_flag ) THEN

          IF ( verbose_level .GE. 0 ) THEN

             WRITE(*,*) 'CARRIER PHASE: liquid'
             WRITE(*,*) 'Carrier phase kinematic viscosity:',kin_visc_l

          END IF

          kin_visc_c = kin_visc_l
          r_rho_c = rho_l
          r_inv_rho_c = inv_rho_l
          rho_c = CMPLX(rho_l,0.0_wp,wp)
          inv_rho_c = CMPLX(inv_rho_l,0.0_wp,wp)
          sp_heat_c = sp_heat_l

       END IF

    END IF

    ! ------- READ solid_transport_parameters NAMELIST --------------------------
    
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
          subtract_init_flag = .FALSE.
          
       ELSE
       
          IF ( ( erodible_porosity .LT. 0.0_wp ) .OR.                           &
               ( erodible_porosity .GT. 1.0_wp ) ) THEN
             
             WRITE(*,*) 'ERROR: problem with namelist SOLID_TRANSPORT_PARAMETERS'
             WRITE(*,*) 'erodible_porosity =' , erodible_porosity
             WRITE(*,*) 'Please check the input file'
             STOP
          
          END IF

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
          
    ALLOCATE( erodible( comp_cells_x , comp_cells_y , n_solid ) )
    erodible(1:comp_cells_x,1:comp_cells_y,1:n_solid ) = 0.0_wp

    IF ( TRIM(erodible_file) .EQ. '' ) THEN

       erodible(:,:,1:n_solid) = 1.0E+5_wp

       IF ( verbose_level .GE. 0.0_wp ) THEN
     
          WRITE(*,*)
          WRITE(*,*) 'WARNING: no file defined for erobile layer'
          WRITE(*,*) 'Unlimited erosion'
          WRITE(*,*)

       END IF
          
    ELSEIF ( erosion_coeff .EQ. 0.0_wp ) THEN

       WRITE(*,*) 'WARNING: erodible_file not used'
       WRITE(*,*) 'erosion_coeff = ', erosion_coeff

       erodible(:,:,1:n_solid) = 1.0E+5_wp
       
    ELSE
       
       IF ( verbose_level .GE. 0.0_wp ) THEN

          WRITE(*,*) 'Maximum thick. for erosion read from : ',                 &
               TRIM(erodible_file)

       END IF
                     
       CALL read_erodible
       
    END IF
    
    n_vars = n_vars + n_solid
    n_eqns = n_vars

    alphas_bcW(1:n_solid)%flag = -1
    alphas_bcE(1:n_solid)%flag = -1
    alphas_bcS(1:n_solid)%flag = -1
    alphas_bcN(1:n_solid)%flag = -1

       
    ALLOCATE( bcW(n_vars) , bcE(n_vars) , bcS(n_vars) , bcN(n_vars) )

    bcW(1:n_vars)%flag = -1
    bcE(1:n_vars)%flag = -1
    bcS(1:n_vars)%flag = -1
    bcN(1:n_vars)%flag = -1

    ! rho_s(1:n_solid) = rho0_s(1:n_solid)    

    ! diam_s(1:n_solid) = diam0_s(1:n_solid)

    ! sp_heat_s(1:n_solid) = sp_heat0_s(1:n_solid)

    ALLOCATE( inv_rho_s(n_solid) )
    ALLOCATE( c_inv_rho_s(n_solid) )

    ALLOCATE( alphas_init(n_solid) )

    ALLOCATE( alphas_E(n_solid) , alphas_W(n_solid) )
    
    inv_rho_s(1:n_solid) = 1.0_wp / rho_s(1:n_solid)
    
    DO i_solid=1,n_solid
       
       c_inv_rho_s(i_solid) = CMPLX(inv_rho_s(i_solid),0.0_wp,wp)

    END DO

    ALLOCATE( deposit( comp_cells_x , comp_cells_y , n_solid ) )    
    deposit(1:comp_cells_x,1:comp_cells_y,1:n_solid ) = 0.0_wp

    ALLOCATE( erosion( comp_cells_x , comp_cells_y , n_solid ) )
    erosion(1:comp_cells_x,1:comp_cells_y,1:n_solid ) = 0.0_wp
       
    IF ( restart ) THEN

       ! ---------- READ restart_parameters NAMELIST ----------------------------
       READ(input_unit,restart_parameters,IOSTAT=ios)

       IF ( ios .NE. 0 ) THEN
          
          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist RESTART_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          STOP
          
       ELSE
  
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

    ELSE

       IF ( riemann_flag ) THEN

          READ(input_unit,left_state,IOSTAT=ios)
          
          alphas_E(1:n_solid) = alphas0_E(1:n_solid)

          IF ( ios .NE. 0 ) THEN
             
             WRITE(*,*) 'IOSTAT=',ios
             WRITE(*,*) 'ERROR: problem with namelist LEFT_PARAMETERS'
             WRITE(*,*) 'Please check the input file'
             STOP
             
          ELSE
             
             REWIND(input_unit)

             IF ( T_W .EQ. -1.0_wp ) THEN

                WRITE(*,*) 'ERROR: problem with namelist LEFT_PARAMETERS'
                WRITE(*,*) 'Initial temperature T_W not defined'
                STOP
                
             END IF
             
          END IF
          
          READ(input_unit,right_state,IOSTAT=ios)

          alphas_E(1:n_solid) = alphas0_E(1:n_solid)
          
          IF ( ios .NE. 0 ) THEN
             
             WRITE(*,*) 'IOSTAT=',ios
             WRITE(*,*) 'ERROR: problem with namelist RIGHT_PARAMETERS'
             WRITE(*,*) 'Please check the input file'
             STOP
             
          ELSE
             
             REWIND(input_unit)

             IF ( T_E .EQ. -1.0_wp ) THEN

                WRITE(*,*) 'ERROR: problem with namelist RIGHT_PARAMETERS'
                WRITE(*,*) 'Initial temperature T_E not defined'
                STOP
                
             END IF
             
             
          END IF
            
       ELSE


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

    IF ( ( MAXVAL(limiter(1:n_vars)) .GT. 3 ) .OR.                              &
         ( MINVAL(limiter(1:n_vars)) .LT. 0 ) ) THEN

       WRITE(*,*) 'WARNING: wrong limiter ',limiter(1:n_vars)
       WRITE(*,*) 'Choose among: none, minmod,superbee,van_leer'
       STOP         

    END IF

    IF ( verbose_level .GE. 0 ) THEN
       
       WRITE(*,*) 'Linear reconstruction and b. c. applied to variables:'
       WRITE(*,*) 'h,hu,hv,T,alphas'

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
       
       IF ( ANY(alphas_bcW(1:n_solid)%flag .EQ. -1 ) ) THEN 
          
          WRITE(*,*) 'ERROR: problem with namelist WEST_BOUNDARY_CONDITIONS'
          WRITE(*,*) 'B.C. for sediment conentration not set properly'
          WRITE(*,*) 'Please check the input file'
          WRITE(*,*) 'alphas_bcW'
          WRITE(*,*) alphas_bcW(1:n_solid)
          STOP
          
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
       
       IF ( ANY(alphas_bcE(1:n_solid)%flag .EQ. -1 ) ) THEN 
          
          WRITE(*,*) 'ERROR: problem with namelist EAST_BOUNDARY_CONDITIONS'
          WRITE(*,*) 'B.C. for sediment concentration not set properly'
          WRITE(*,*) 'Please check the input file'
          WRITE(*,*) 'alphas_bcE'
          WRITE(*,*) alphas_bcE(1:n_solid)
          STOP
          
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
       
       IF ( ANY(alphas_bcS(1:n_solid)%flag .EQ. -1 ) ) THEN 
          
          WRITE(*,*) 'ERROR: problem with namelist SOUTH_BOUNDARY_CONDITIONS'
          WRITE(*,*) 'B.C. for sediment concentrations not set properly'
          WRITE(*,*) 'Please check the input file'
          WRITE(*,*) 'alphas_bcS'
          WRITE(*,*) alphas_bcS(1:n_solid)
          STOP
          
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
       
       IF ( ANY(alphas_bcN(1:n_solid)%flag .EQ. -1 ) ) THEN 
          
          WRITE(*,*) 'ERROR: problem with namelist NORTH_BOUNDARY_CONDITIONS'
          WRITE(*,*) 'B.C. for sediment concentrations not set properly'
          WRITE(*,*) 'Please check the input file'
          WRITE(*,*) 'alphas_bcN'
          WRITE(*,*) alphas_bcN(1:n_solid)
          STOP
          
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

    bcW(5:4+n_solid) = alphas_bcW(1:n_solid)
    bcE(5:4+n_solid) = alphas_bcE(1:n_solid)
    bcS(5:4+n_solid) = alphas_bcS(1:n_solid)
    bcN(5:4+n_solid) = alphas_bcN(1:n_solid)
       
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

    ! ------- READ radial_source_parameters NAMELIST ----------------------------

    IF ( radial_source_flag ) THEN

       alphal_source = -1.0_wp

       READ(input_unit,radial_source_parameters,IOSTAT=ios)

       IF ( ios .NE. 0 ) THEN
             
          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'ERROR: problem with namelist RADIAL_SOURCE_PARAMETERS'
          WRITE(*,*) 'Please check the input file'
          STOP
          
       ELSE
          
          REWIND(input_unit)
          
          IF ( t_source .EQ. -1.0_wp ) THEN
             
             WRITE(*,*) 'ERROR: problem with namelist RADIAL_SOURCE_PARAMETERS'
             WRITE(*,*) 'PLEASE CHEC VALUE OF T_SOURCE',t_source
             STOP
             
          END IF

          IF ( h_source .EQ. -1.0_wp ) THEN
             
             WRITE(*,*) 'ERROR: problem with namelist RADIAL_SOURCE_PARAMETERS'
             WRITE(*,*) 'PLEASE CHEC VALUE OF H_SOURCE',h_source
             STOP
             
          END IF

          IF ( r_source .EQ. -1.0_wp ) THEN
             
             WRITE(*,*) 'ERROR: problem with namelist RADIAL_SOURCE_PARAMETERS'
             WRITE(*,*) 'PLEASE CHEC VALUE OF R_SOURCE',r_source
             STOP
             
          END IF

          IF ( vel_source .EQ. -1.0_wp ) THEN
             
             WRITE(*,*) 'ERROR: problem with namelist RADIAL_SOURCE_PARAMETERS'
             WRITE(*,*) 'PLEASE CHEC VALUE OF VEL_SOURCE',vel_source
             STOP
             
          END IF
          
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

          IF ( gas_flag .AND. liquid_flag ) THEN

             IF ( alphal_source .LT. 0.0_wp ) THEN
             
                WRITE(*,*) 'ERROR: problem with namelist ',                     &
                     'RADIAL_SOURCE_PARAMETERS'
                WRITE(*,*) 'PLEASE CHECK VALUE OF ALPHAL_SOURCE',alphal_source
                STOP
                
             END IF

          END IF

          IF ( ANY(alphas_source(1:n_solid) .EQ. -1.0_wp ) ) THEN
       
             WRITE(*,*) 'ERROR: problem with namelist RADIAL_VOLUME_PARAMETERS'
             WRITE(*,*) 'alphas_source =' , alphas_source(1:n_solid)
             WRITE(*,*) 'Please check the input file'
             STOP
             
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
   
       END IF
           
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

          
       END IF
           
    END IF


    ! ------- READ temperature_parameters NAMELIST ------------------------------

!!$    READ(input_unit, temperature_parameters,IOSTAT=ios)
!!$       
!!$    IF ( ios .NE. 0 ) THEN
!!$       
!!$       WRITE(*,*) 'IOSTAT=',ios
!!$       WRITE(*,*) 'ERROR: problem with namelist TEMPERATURE_PARAMETERS'
!!$       WRITE(*,*) 'Please check the input file'
!!$       STOP
!!$       
!!$    ELSE
!!$       
!!$       REWIND(input_unit)
!!$       
!!$    END IF

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
               ( visc_par .NE. -1.0_wp ) .OR. ( tau .NE. -1.0_wp ) ) THEN

             WRITE(*,*) 'WARNING: parameters not used in RHEOLOGY_PARAMETERS'
             IF ( T_ref .NE. -1.0_wp ) WRITE(*,*) 'T_ref =',T_ref 
             IF ( nu_ref .NE. -1.0_wp ) WRITE(*,*) 'nu_ref =',nu_ref 
             IF ( visc_par .NE. -1.0_wp ) WRITE(*,*) 'visc_par =',visc_par
             IF ( tau .NE. -1.0_wp ) WRITE(*,*) 'tau =',tau 
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
          
          IF ( ANY(sed_vol_perc(1:n_solid) .EQ. -1.0_wp ) ) THEN
             
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
             
       ELSEIF ( rheology_model .EQ. 6 ) THEN

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
             
       ELSE

             WRITE(*,*) 'ERROR: problem with namelist RHEOLOGY_PARAMETERS'
             WRITE(*,*) 'RHEOLOGY_MODEL =' , rheology_model
             WRITE(*,*) 'Please check the input file'
             STOP
             
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
    READ(2001,*) chara, nodata_value

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

    topography_profile(3,:,:) = MAX(0.0_wp,topography_profile(3,:,:))

    IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) ''

    CLOSE(2001)


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
       
       runout_file = TRIM(run_name)//'_runout'//'.txt'

       OPEN(runout_unit,FILE=runout_file,STATUS='unknown',form='formatted')
  
    END IF

    ! ----------- READ vulnerability_table_parameters NAMELIST ------------------

    READ(input_unit, vulnerability_table_parameters,IOSTAT=ios)
    
    IF ( ios .NE. 0 ) THEN

       IF ( verbose_level .GE. 0.0_wp ) THEN
       
          WRITE(*,*) 'IOSTAT=',ios
          WRITE(*,*) 'WARNING: namelist VULNERABILITY_TABLE_PARAMETERS not found'

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
       

    !------ search for check points --------------------------------------------

    tend1 = .FALSE.
    
    IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'Searching for topography_profile'
   
    n_probes = 0
 
    probes_search: DO
       
       READ(input_unit,*, END = 300 ) card
       
       IF( TRIM(card) == 'PROBES_COORDS' ) THEN
          
          EXIT probes_search
          
       END IF
       
    END DO probes_search
  
    
    READ(input_unit,*) n_probes
    
    WRITE(*,*) 'n_probes ',n_probes
    
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

       IF ( riemann_flag) THEN

          WRITE(backup_unit,left_state)
          WRITE(backup_unit,right_state)

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


    IF ( output_runout_flag ) WRITE(backup_unit, runout_parameters)

    IF ( ( COUNT( thickness_levels0 .GT. 0.0_wp ) .GT. 0 ) .OR.                 &
         ( COUNT( dyn_pres_levels0 .GT. 0.0_wp ) .GT. 0 ) ) THEN

       WRITE(backup_unit, vulnerability_table_parameters)

    END IF
       
    IF ( n_probes .GT. 0 ) THEN
       
       WRITE(backup_unit,*) '''PROBES_COORDS'''
       WRITE(backup_unit,*) n_probes
       
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

    IF ( dt_output_org .NE. dt_output ) THEN

       WRITE(*,*) 'Modified input file: dt_output =',dt_output

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
    USE solver_2d, ONLY : allocate_solver_variables

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

    REAL(wp), ALLOCATABLE :: thickness_input(:,:)

    REAL(wp), ALLOCATABLE :: x1(:) , y1(:)

    REAL(wp) :: xl , xr , yl , yr 
    
    REAL(wp) :: rho_c , rho_m , mass_fract(n_solid)

    REAL(wp) :: sp_heat_c

    INTEGER :: solid_idx

    INTEGER :: i_vars , i_solid

    INQUIRE (FILE=restart_file,exist=lexist)

    WRITE(*,*)
    ! WRITE(*,*) 'READ INIT',restart_file,lexist,restart_unit

    IF ( lexist .EQV. .FALSE.) THEN

       WRITE(*,*) 'Restart: ',TRIM(restart_file) , ' not found'
       STOP

    ELSE

       OPEN(restart_unit,FILE=restart_file,STATUS='old')
       
       IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'Restart: ',TRIM(restart_file),   &
            ' found'

    END IF

    dot_idx = SCAN(restart_file, ".", .TRUE.)

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
       
       READ(restart_unit,*) chara, ncols
       READ(restart_unit,*) chara, nrows
       READ(restart_unit,*) chara, xllcorner
       READ(restart_unit,*) chara, yllcorner
       READ(restart_unit,*) chara, cellsize
       READ(restart_unit,*) chara, nodata_value
       
       ALLOCATE( thickness_init(comp_cells_x,comp_cells_y) )
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
       IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'Total volume from restart =',    &
            cellsize**2*SUM(thickness_input)

       WHERE ( thickness_input .EQ. nodata_value )

          thickness_input = 0.0_wp
          
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
             
             CALL regrid_scalar( x1 , y1 , thickness_input , xl , xr , yl ,     &
                  yr , thickness_init(j,k) )

          END DO

       END DO

       IF ( subtract_init_flag ) THEN
          
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

             WRITE(*,*) 'Absolute fractions of solide phases in deposit:'
             
             DO i_solid=1,n_solid

                WRITE(*,*) erodible_fract(i_solid) * ( 1.0_wp-erodible_porosity )
                erodible(:,:,i_solid) = erodible_fract(i_solid) *               &
                     ( 1.0_wp - erodible_porosity ) * erodible_init(:,:)
                
             END DO

             WRITE(*,*)
             
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
       q(2,:,:) = 0.0_wp
       ! rhom*h*v
       q(3,:,:) = 0.0_wp

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

       DEALLOCATE( thickness_input )

       WRITE(*,*) 'n_vars',n_vars
       
    ELSEIF ( check_file .EQ. 'q_2' ) THEN
    
       DO k=1,comp_cells_y
          
          DO j=1,comp_cells_x

             READ(restart_unit,'(2e20.12,100(e20.12))') xj , yk ,               &
                  (q(i_vars,j,k),i_vars=1,n_vars) 

             IF ( q(1,j,k) .LE. 0.0_wp ) q(1:n_vars,j,k) = 0.0_wp

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
    
    ELSE
   
       WRITE(*,*) 'Restart file not in the right format (*.asc or *)'
       STOP

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

    
    IF ( TRIM(erodible_file) .EQ. '' ) THEN

       DO i_solid=1,n_solid

          erodible(:,:,i_solid) = 1.0E+5_wp

       END DO
       
    ELSE

       INQUIRE (FILE=erodible_file,exist=lexist)

       WRITE(*,*)

       IF ( lexist .EQV. .FALSE.) THEN
          
          WRITE(*,*) 'Erodible file: ',TRIM(erodible_file) , ' not found'
          STOP
          
       END IF
       
       OPEN(erodible_unit,FILE=erodible_file,STATUS='old')
       
       IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'Erodible file: ',                &
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
       
       ALLOCATE( erodible_init(comp_cells_x,comp_cells_y) )
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
          WRITE(*,*) 'WARNING: changing resolution of erodible layer' 
          WRITE(*,*) 'cellsize not equal to cell_size', cellsize , cell_size
          WRITE(*,*)

       END IF

       DO k=1,nrows

          WRITE(*,FMT="(A1,A,t21,F6.2,A)",ADVANCE="NO") ACHAR(13),              &
               & " Percent Complete: ",( REAL(k) / REAL(nrows))*100.0, "%"
          
          READ(erodible_unit,*) erodible_input(1:ncols,nrows-k+1)
          
       ENDDO

       WRITE(*,*) 
       IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'Total erodible volume =',        &
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
             
             CALL regrid_scalar( x1 , y1 , erodible_input , xl , xr , yl ,     &
                  yr , erodible_init(j,k) )
             
          END DO

       END DO

       IF ( VERBOSE_LEVEL .GE. 0 ) THEN
          
          WRITE(*,*) 'Total erodible volume on computational grid =' ,          &
               cell_size**2 * SUM( erodible_init(:,:) )

       END IF

       DO i_solid=1,n_solid

          erodible(:,:,i_solid) = erodible_fract(i_solid) * erodible_init(:,:)

       END DO
              
    END IF

      
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
    USE constitutive_2d, ONLY : qc_to_qp, mixt_var

    USE geometry_2d, ONLY : comp_cells_x , B_cent , comp_cells_y , x_comp,      &
         y_comp , deposit , erosion

    USE parameters_2d, ONLY : n_vars
    USE parameters_2d, ONLY : t_output , dt_output 
    USE parameters_2d, ONLY : t_steady

    USE solver_2d, ONLY : q

    IMPLICIT none

    REAL(wp), INTENT(IN) :: time

    CHARACTER(LEN=4) :: idx_string

    REAL(wp) :: qp(n_vars+2)

    REAL(wp) :: B_out

    REAL(wp) :: r_u , r_v , r_h , r_alphas(n_solid) , r_T , r_Ri , r_rho_m
    REAL(wp) :: r_rho_c      !< real-value carrier phase density [kg/m3]
    REAL(wp) :: r_red_grav   !< real-value reduced gravity


    INTEGER :: j,k
    INTEGER :: i
    INTEGER :: i_vars
    
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
       
       DO k = 1,comp_cells_y
          
          DO j = 1,comp_cells_x
          
             CALL qc_to_qp(q(1:n_vars,j,k) , qp(1:n_vars+2))

             CALL mixt_var(qp(1:n_vars+2),r_Ri,r_rho_m,r_rho_c,r_red_grav)

             r_h = qp(1)
             r_u = qp(n_vars+1)
             r_v = qp(n_vars+2)
             r_T = qp(4)
             r_alphas(1:n_solid) = qp(5:4+n_solid)

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

             WRITE(output_unit_2d,1010) x_comp(j), y_comp(k), r_h , r_u , r_v , &
                  B_out , r_h + B_out , r_alphas , r_T , r_rho_m , r_red_grav , &
                  DEPOSIT(j,k,:) , EROSION(j,k,:)

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
    
    IF ( n_probes .GT. 0 ) CALL output_probes(output_idx)

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

       WRITE(output_max_unit,*) grid_output(1:comp_cells_x,j)

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

                WRITE(output_max_unit,*) grid_output_int(1:comp_cells_x,j)

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

    USE geometry_2d, ONLY : B_cent , grid_output , deposit , erosion
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

       WRITE(output_esri_unit,*) grid_output(1:comp_cells_x,j)

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
       
       WRITE(output_esri_unit,*) grid_output(1:comp_cells_x,j)
       
    ENDDO
    
    CLOSE(output_esri_unit)

    !Save solid classes deposit
    DO i_solid=1,n_solid

       isolid_string = lettera(i_solid)

       output_esri_file =                                                       &
            TRIM(run_name)//'_dep_'//isolid_string//'_'//idx_string//'.asc'
       
       IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'WRITING ',output_esri_file
       
       OPEN(output_esri_unit,FILE=output_esri_file,status='unknown',            &
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
          
          WRITE(output_esri_unit,'(2000ES12.3E3)') grid_output(1:comp_cells_x,j)
          
       ENDDO
       
       CLOSE(output_esri_unit)
       
    END DO

    ! Save total deposit (solid+continous phase associated with porosity)
    output_esri_file =                                                          &
         TRIM(run_name)//'_depTot_'//idx_string//'.asc'
    
    IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'WRITING ',output_esri_file
    
    OPEN(output_esri_unit,FILE=output_esri_file,status='unknown',               &
         form='formatted')
    
    grid_output = -9999 
    
    WHERE ( SUM(deposit(:,:,:),DIM=3) .GT. 0.0_wp )
       
       grid_output = SUM(deposit(:,:,:),DIM=3) / ( 1.0_wp - erodible_porosity )  
       
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
    
    !Save erosion
    DO i_solid=1,n_solid

       isolid_string = lettera(i_solid)

       output_esri_file =                                                       &
            TRIM(run_name)//'_ers_'//isolid_string//'_'//idx_string//'.asc'
       
       IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'WRITING ',output_esri_file
       
       OPEN(output_esri_unit,FILE=output_esri_file,status='unknown',            &
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
          
          WRITE(output_esri_unit,'(2000ES12.3E3)') grid_output(1:comp_cells_x,j)
          
       ENDDO
       
       CLOSE(output_esri_unit)
       
    END DO

    ! Save total erosion (solid+continous phase associated with porosity)
    output_esri_file =                                                          &
         TRIM(run_name)//'_ersTot_'//idx_string//'.asc'
    
    IF ( VERBOSE_LEVEL .GE. 0 ) WRITE(*,*) 'WRITING ',output_esri_file
    
    OPEN(output_esri_unit,FILE=output_esri_file,status='unknown',               &
         form='formatted')
    
    grid_output = -9999 
    
    WHERE ( SUM(erosion(:,:,:),DIM=3) .GT. 0.0_wp )
       
       grid_output = SUM(erosion(:,:,:),DIM=3) / ( 1.0_wp - erodible_porosity )  
       
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
   
    RETURN
       
  END SUBROUTINE output_esri

  SUBROUTINE close_units

    IMPLICIT NONE

    IF ( output_runout_flag) CLOSE(runout_unit)

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

  CHARACTER*4 FUNCTION lettera(k)
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

  SUBROUTINE output_probes(output_idx)

    USE geometry_2d, ONLY : x_comp , y_comp 
    USE solver_2d, ONLY : q

    USE geometry_2d, ONLY : interp_2d_scalarB

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: output_idx

    CHARACTER(LEN=4) :: idx_string

    REAL(wp) :: f2

    INTEGER :: k 

    idx_string = lettera(output_idx-1)

    !Save thickness
    probes_file = TRIM(run_name)//'_'//idx_string//'.prb'

    OPEN(probes_unit,FILE=probes_file,status='unknown',form='formatted')
    
    DO k=1,n_probes

       CALL interp_2d_scalarB( x_comp , y_comp , q(1,:,:)  ,                    &
            probes_coords(1,k) , probes_coords(2,k) , f2 )

       WRITE(probes_unit,'(3e20.12)') probes_coords(1,k) , probes_coords(2,k) ,f2

    END DO

    CLOSE(probes_unit)

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
    USE solver_2d, ONLY : q , q0 , dt


    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: time
    LOGICAL, INTENT(INOUT) :: stop_flag

    REAL(wp), ALLOCATABLE :: X(:,:), Y(:,:) , dist(:,:)
    INTEGER :: sX, sY
    INTEGER :: imax(2) , imin(2)

    INTEGER :: j,k

    REAL(wp) :: max_mom
    REAL(wp) :: area , area0 , dareaRel_dt
    REAL(wp) :: dhRel_dt

    sX = size(x_comp) 
    sY = size(y_comp) 

    ALLOCATE( X(sX,sY) , Y(sX,sY) , dist(sX,sY) )

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

       ALLOCATE( h_old(sX,sY) )

       h_old(:,:) = q(1,:,:)
      
       IF ( ( x0_runout .EQ. -1 ) .AND. ( y0_runout .EQ. -1 ) ) THEN
          
          WHERE( q(1,:,:) > 1.0E-5_wp ) dist = B_cent
          imin = MAXLOC( dist )
          
          x0_runout = X(imin(1),imin(2))
          y0_runout = Y(imin(1),imin(2))

          WRITE(*,*) 'Runout calculated as linear distance from: (' ,           &
               x0_runout ,',',y0_runout,')'

          dist(:,:) = 0.0_wp
          
          WHERE( q(1,:,:) >1.0E-5_wp ) dist = SQRT( (X-x0_runout)**2 &
               + ( Y - y0_runout )**2 )

          imax = MAXLOC( dist )
          
          init_runout = dist(imax(1),imax(2))

       END IF

       max_mom = 0.0_wp
       
    ELSE

       WHERE( h_old(:,:) > 1.0E-5_wp ) dist = SQRT( q(2,:,:)**2 + q(3,:,:)**2 )  

       max_mom = MAXVAL( dist )

       h_old(:,:) = q(1,:,:)

    END IF

    dist(:,:) = 0.0_wp

    WHERE( q(1,:,:)  > 1.0E-5_wp ) dist = SQRT( ( X - x0_runout )**2 &
         + ( Y - y0_runout )**2 )

    imax = MAXLOC( dist )

    OPEN(dakota_unit,FILE='dakota.txt',status='replace',form='formatted')
    
    WRITE(dakota_unit,*) 'final runout =', dist(imax(1),imax(2)) - init_runout
    
    CLOSE(dakota_unit)

    area0 = dx*dy*COUNT( q0(1,:,:) .GT. 1.0E-7_wp )
    area = dx*dy*COUNT( q(1,:,:) .GT. 1.0E-7_wp )
    
    WRITE(runout_unit,'(A,F12.3,A,F12.3,A,F14.3)') 'Time (s) = ',time ,         &
         ' Runout (m) = ',dist(imax(1),imax(2)) - init_runout,' Area (m^2) = ', &
         area
    
    CALL flush(runout_unit)

    dist(:,:) = 0.0_wp
    
    WHERE( q(1,:,:)  > 1.0E-5_wp ) dist = ABS( q(1,:,:) - q0(1,:,:) ) /         &
         MAX(q(1,:,:),q0(1,:,:))

    IF ( time .GT. t_start ) THEN

       dareaRel_dt = ABS( area - area0 ) / ( area * dt )

       dhRel_dt = SUM( dist / dt ) / COUNT( q(1,:,:) .GT. 1.0E-5_wp )

       IF ( ( MAX(dareaRel_dt,dhRel_dt) .LT. eps_stop ) .AND.                   &
            (.NOT.stop_flag) ) THEN

          WRITE(*,*) 'Steady solution reached'
          WRITE(*,*) 'dareaRel_dt',dareaRel_dt
          WRITE(*,*) 'dhRel_dt',dhRel_dt
          stop_flag = .TRUE.

       END IF

    END IF

    DEALLOCATE( X , Y , dist )

    t_runout = time + dt_runout

  END SUBROUTINE output_runout

END MODULE inpout_2d

