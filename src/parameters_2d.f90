!********************************************************************************
!> \brief Parameters
!
!> This module contains the parameters for numerical solution of the
!> model.
!********************************************************************************
MODULE parameters_2d

  IMPLICIT NONE

  INTEGER, PARAMETER :: dp = KIND(1.0)

  REAL(dp) :: eps_newton        !< threshold for the convergence of the
                                !< Newton's method 

  REAL(dp) :: dt0               !< Initial time step

  REAL(dp) :: max_dt            !< Largest time step allowed

  REAL(dp) :: cfl               !< Courant-Friedrichs-Lewy parameter 

  REAL(dp) :: eps_sing          !< parameter for desingularization

  REAL(dp) :: reconstr_coeff    !< Slope coefficient in the linear reconstruction

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

  REAL(dp) :: x_source
  REAL(dp) :: y_source
  REAL(dp) :: r_source
  REAL(dp) :: vel_source
  REAL(dp) :: T_source
  REAL(dp) :: h_source
  REAL(dp) :: alphas_source(100)
  REAL(dp) :: alphal_source
  REAL(dp) :: time_param(4)
  

  LOGICAL :: collapsing_volume_flag

  REAL(dp) :: x_collapse
  REAL(dp) :: y_collapse
  REAL(dp) :: r_collapse
  REAL(dp) :: T_collapse
  REAL(dp) :: h_collapse
  REAL(dp) :: alphas_collapse(100)


  !> Initial volume of the flow
  REAL(dp) :: released_volume

  !> Initial x-coordiante of the pile
  REAL(dp) :: x_release

  !> Initial y-coordinate of the pile
  REAL(dp) :: y_release
  
    !> Initial velocity module of the pile
  REAL(dp) :: velocity_mod_release

  !> Initial velocity direction (angle in degree):\n
  !> - >=0    => departing from positive x-axis
  !> - <0     => departign from maximum slope direction
  !.
  
  REAL(dp) :: velocity_ang_release

  !> Initial temperature of the pile of material
  REAL(dp) :: T_init

  !> Initial sediment concentration in the pile of material
  REAL(dp), ALLOCATABLE :: alphas_init(:)

  INTEGER :: n_vars   !< Number of conservative variables
  INTEGER :: n_eqns   !< Number of equations

  INTEGER :: n_solid  !< Number of solid classes

  INTEGER :: n_nh     !< Number of non-hyperbolic terms

  INTEGER :: n_RK     !< Runge-Kutta order
  
  INTEGER, PARAMETER :: max_nl_iter = 100

  REAL(dp), PARAMETER :: tol_abs = 1.D-5
  REAL(dp), PARAMETER :: tol_rel = 1.D-5

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

  REAL(dp) :: theta             !< Van Leer limiter parameter
  REAL(dp) :: t_start           !< initial time for the run
  REAL(dp) :: t_end             !< end time for the run
  REAL(dp) :: t_output          !< time of the next output
  REAL(dp) :: dt_output         !< time interval for the output of the solution
  REAL(dp) :: t_runout          !< time of the next runout output
  REAL(dp) :: t_steady          !< end time when reached steady solution

  INTEGER :: verbose_level

  TYPE bc
     INTEGER :: flag
     REAL(dp) :: value
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
