!********************************************************************************
!> \brief Parameters
!
!> This module contains the parameters for numerical solution of the
!> model.
!********************************************************************************
MODULE parameters_2d

  IMPLICIT NONE

  INTEGER, PARAMETER :: sp = Selected_Real_Kind (P=6,R=37)
  INTEGER, PARAMETER :: dp = Selected_Real_Kind (P=15,R=307)

  !> working precision
  INTEGER, PARAMETER :: wp = dp

  REAL, PARAMETER :: four_thirds = 1.0_wp / 3.0_wp
  REAL, PARAMETER :: neg_four_thirds = -1.0_wp / 3.0_wp

  REAL(wp), PARAMETER :: tolh = 10.0_wp * EPSILON(1.0_wp)

  REAL(wp) :: eps_newton        !< threshold for the convergence of the
                                !< Newton's method 

  REAL(wp) :: dt0               !< Initial time step

  REAL(wp) :: max_dt            !< Largest time step allowed

  REAL(wp) :: cfl               !< Courant-Friedrichs-Lewy parameter 

  REAL(wp) :: eps_sing          !< parameter for desingularization
  REAL(wp) :: eps_sing4         !< parameter for desingularization**4

  REAL(wp) :: reconstr_coeff    !< Slope coefficient in the linear reconstruction

  !> Flag to add the relaxation terms after the linear reconstruction:\n
  !> - T      => evaluate the relaxation terms
  !> - F      => reconstruction without the relaxation 
  !> .
  LOGICAL :: interfaces_relaxation

  !> Flag to choose in which way we upload the topography
  !> - T      => through a function
  !> - F      => through points
  !> .
  LOGICAL :: topography_function_flag

  !> Flag for uploading topography from a different file (topography_dem.asc)
  !> - T      => from dem
  !> - F      => from points in file.inp
  !> .
  LOGICAL :: topography_demfile

  !> Flag to choose the equation for temperature to solve
  !> - T      => solve the full energy equation
  !> - F      => solve for a simpler transport equation (advection) for temperature
  !> .
  LOGICAL :: energy_flag

  !> Flag to choose the sort of problem to solve
  !> - T      => riemann problem
  !> - F      => generic initial conditions (uploaded through functions, to be defined in inpout_2d.f90)
  !> .
  LOGICAL :: riemann_flag

  !> Flag to choose if we add the rheology
  !> - T      => rheology activated
  !> - F      => no rheology
  !> .
  LOGICAL :: rheology_flag
  
  !> choice of the rheology model
  !> - 1      => Voellmy-Salm rheology
  !> - 2      => plastic rheology
  !> .
  INTEGER :: rheology_model

  LOGICAL :: liquid_flag

  LOGICAL :: gas_flag
  
  LOGICAL :: topo_change_flag

  LOGICAL :: radial_source_flag

  !> Flag to choose if initial volume is subtracted from topography or erodible layer
  !> - T      => change initial topography or erodible layer by subtracting initial volume 
  !> - F      => keep initial topography or erodible layer
  !> .
  LOGICAL :: subtract_init_flag

  INTEGER :: n_thickness_levels
  INTEGER :: n_dyn_pres_levels
  REAL(wp), ALLOCATABLE :: thickness_levels(:)
  REAL(wp), ALLOCATABLE :: dyn_pres_levels(:)

  REAL(wp) :: x_source
  REAL(wp) :: y_source
  REAL(wp) :: r_source
  REAL(wp) :: vel_source
  REAL(wp) :: T_source
  REAL(wp) :: h_source
  REAL(wp) :: alphas_source(100)
  REAL(wp) :: alphal_source
  REAL(wp) :: time_param(4)
  

  LOGICAL :: collapsing_volume_flag

  REAL(wp) :: x_collapse
  REAL(wp) :: y_collapse
  REAL(wp) :: r_collapse
  REAL(wp) :: T_collapse
  REAL(wp) :: h_collapse
  REAL(wp) :: alphas_collapse(100)


  !> Initial volume of the flow
  REAL(wp) :: released_volume

  !> Initial x-coordiante of the pile
  REAL(wp) :: x_release

  !> Initial y-coordinate of the pile
  REAL(wp) :: y_release
  
    !> Initial velocity module of the pile
  REAL(wp) :: velocity_mod_release

  !> Initial velocity direction (angle in degree):\n
  !> - >=0    => departing from positive x-axis
  !> - <0     => departign from maximum slope direction
  !. 
  REAL(wp) :: velocity_ang_release

  !> Initial temperature of the pile of material
  REAL(wp) :: T_init
  
  !> Initial sediment concentration in the pile of material
  REAL(wp), ALLOCATABLE :: alphas_init(:)

  INTEGER :: n_vars   !< Number of conservative variables
  INTEGER :: n_eqns   !< Number of equations

  INTEGER :: n_solid  !< Number of solid classes

  INTEGER :: n_nh     !< Number of non-hyperbolic terms

  INTEGER :: n_RK     !< Runge-Kutta order
  
  INTEGER, PARAMETER :: max_nl_iter = 100

  REAL(wp), PARAMETER :: tol_abs = 1.0E-5_wp
  REAL(wp), PARAMETER :: tol_rel = 1.0E-5_wp

  !> Limiter for the slope in the linear reconstruction:\n
  !> - 'none'     => no limiter (constant value);
  !> - 'minmod'   => minmod slope;
  !> - 'superbee' => superbee limiter (Roe, 1985);
  !> - 'van_leer' => monotonized central-difference limiter (van Leer, 1977)
  !> .
  INTEGER :: limiter(10) = -1

  !> Finite volume method:\n
  !> - 'LxF'       => lax-friedrichs scheme;
  !> - 'GFORCE '   => gforce scheme;
  !> - 'KT'        => Kurganov and Tadmor semidiscrete scheme;
  !> .
  CHARACTER(LEN=20) :: solver_scheme     

  REAL(wp) :: theta             !< Van Leer limiter parameter
  REAL(wp) :: t_start           !< initial time for the run
  REAL(wp) :: t_end             !< end time for the run
  REAL(wp) :: t_output          !< time of the next output
  REAL(wp) :: dt_output         !< time interval for the output of the solution
  REAL(wp) :: t_runout          !< time of the next runout output
  REAL(wp) :: t_steady          !< end time when reached steady solution

  INTEGER :: verbose_level

  TYPE bc
     INTEGER :: flag
     REAL(wp) :: value
  END TYPE bc

  ! -------boundary conditions variables

  !> bcW&flag defines the west boundary condition:\n
  !> - bcW%flag = 0     => Dirichlet boundary condition;
  !> - bcW%flag = 1     => Neumann boundary condition.
  !> .
  !> bcLWvalue is the value of the left boundary condition:\n
  !> - value of the variable for Dirichlet boundary condition (bcW%flag=0);
  !> - gradient of the variable for Neumann boundary condition (bcW%flag=1).
  !> .
  TYPE(bc), ALLOCATABLE :: bcW(:)

  !> bcE&flag defines the east boundary condition:\n
  !> - bcE%flag = 0     => Dirichlet boundary condition;
  !> - bcE%flag = 1     => Neumann boundary condition.
  !> .
  !> bcE%value is the value of the right boundary condition:\n
  !> - value of the variable for Dirichlet boundary condition (bcE%flag=0);
  !> - gradient of the variable for Neumann boundary condition (bcE%flag=1).
  !> .
  TYPE(bc), ALLOCATABLE :: bcE(:)

  !> bcS&flag defines the south boundary condition:\n
  !> - bcS%flag = 0     => Dirichlet boundary condition;
  !> - bcS%flag = 1     => Neumann boundary condition.
  !> .
  !> bcS%value is the value of the bottom boundary condition:\n
  !> - value of the variable for Dirichlet boundary condition (bcS%flag=0);
  !> - gradient of the variable for Neumann boundary condition (bcS%flag=1).
  !> .
  TYPE(bc), ALLOCATABLE :: bcS(:)

  !> bcN&flag defines the north boundary condition:\n
  !> - bcN%flag = 0     => Dirichlet boundary condition;
  !> - bcN%flag = 1     => Neumann boundary condition.
  !> .
  !> bcN%value is the value of the top boundary condition:\n
  !> - value of the variable for Dirichlet boundary condition (bcN%flag=0);
  !> - gradient of the variable for Neumann boundary condition (bcN%flag=1).
  !> .
  TYPE(bc), ALLOCATABLE :: bcN(:)

END MODULE parameters_2d
