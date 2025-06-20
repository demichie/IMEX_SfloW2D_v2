!********************************************************************************
!> \brief Initial solution
!
!> This module contains the variables and the subroutine for the
!> initialization of the solution for a Riemann problem.
!********************************************************************************

MODULE init_2d

  USE parameters_2d, ONLY : wp
  USE parameters_2d, ONLY : verbose_level
  USE parameters_2d, ONLY : n_solid , n_add_gas
  USE parameters_2d, ONLY : n_stoch_vars , n_pore_vars

  IMPLICIT none

  REAL(wp), ALLOCATABLE :: q_init(:,:,:)

  REAL(wp), ALLOCATABLE :: thickness_init(:,:)

  !> Initial thickness of erodible layer (solid+voids)
  REAL(wp), ALLOCATABLE :: erodible_init(:,:)

CONTAINS

  SUBROUTINE init_empty

    USE constitutive_2d, ONLY : T_ambient

    USE constitutive_2d, ONLY : qp_to_qc

    USE geometry_2d, ONLY : comp_cells_x , comp_cells_y

    USE parameters_2d, ONLY : n_vars

    USE solver_2d, ONLY : q

    IMPLICIT NONE

    INTEGER :: j,k

    REAL(wp) :: qp_init(n_vars+2)

    WRITE(*,*) 'Initialization with zero thickness flow'
    
    qp_init(1:n_vars+2) = 0.0_wp
    qp_init(4) = T_ambient

    DO j = 1,comp_cells_x

       DO k = 1,comp_cells_y

          CALL qp_to_qc( qp_init(1:n_vars+2) , q(1:n_vars,j,k) )


       END DO

    END DO

    RETURN

  END SUBROUTINE init_empty
  
  !******************************************************************************
  !> \brief Collapsing volume initialization
  !
  !> This subroutine initialize the solution for a collpasing volume. Values for 
  !>  the initial state (x, y, r, T, h, alphas) are read from the input file.
  !> \date 2019_12_11
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE collapsing_volume

    USE constitutive_2d, ONLY : qp_to_qc

    USE geometry_2d, ONLY : compute_cell_fract

    USE geometry_2d, ONLY : comp_cells_x , comp_cells_y

    USE parameters_2d, ONLY : n_vars , alpha_flag

    USE parameters_2d, ONLY : x_collapse , y_collapse , r_collapse , T_collapse , &
       h_collapse , alphas_collapse , alphag_collapse

    USE solver_2d, ONLY : q

    IMPLICIT NONE

    INTEGER :: j,k

    REAL(wp) :: qp_init(n_vars+2) ,  qp0_init(n_vars+2)

    REAL(wp) :: cell_fract(comp_cells_x,comp_cells_y)

    CALL compute_cell_fract(x_collapse,y_collapse,r_collapse,r_collapse,0.0_wp,cell_fract)

    ! values outside the collapsing volume
    qp0_init(1) = 0.0_wp                  ! h
    qp0_init(2) = 0.0_wp                  ! hu
    qp0_init(3) = 0.0_wp                  ! hv
    qp0_init(4) = T_collapse              ! T
    qp0_init(5:4+n_solid) = 0.0_wp        ! alphas
    qp0_init(4+n_solid+1:4+n_solid+n_add_gas) = 0.0_wp        ! alphag
    qp0_init(5+n_solid+n_add_gas:4+n_solid+n_add_gas+n_stoch_vars) = 0.0_wp
    qp0_init(5+n_solid+n_add_gas+n_stoch_vars:4+n_solid+n_add_gas+n_stoch_vars+ &
         n_pore_vars) =  0.0_wp
    qp0_init(n_vars+1:n_vars+2) = 0.0_wp  ! u,v

    ! values within the collapsing volume
    qp_init(2) = 0.0_wp
    qp_init(3) = 0.0_wp
    qp_init(4) = T_collapse

    qp_init(n_vars+1:n_vars+2) = 0.0_wp
    
    DO j = 1,comp_cells_x
       
       DO k = 1,comp_cells_y
          
          IF ( cell_fract(j,k) .GT. 0.0_wp ) THEN
             
             qp_init(1) = cell_fract(j,k) * h_collapse

             IF ( alpha_flag ) THEN

                qp_init(5:4+n_solid) = cell_fract(j,k)*alphas_collapse(1:n_solid)
                qp_init(4+n_solid+1:4+n_solid+n_add_gas) = cell_fract(j,k) *    &
                     alphag_collapse(1:n_add_gas)
       
             ELSE

                qp_init(5:4+n_solid) = cell_fract(j,k) * h_collapse *           &
                     alphas_collapse(1:n_solid)
                qp_init(4+n_solid+1:4+n_solid+n_add_gas) = cell_fract(j,k) *    &
                     h_collapse * alphag_collapse(1:n_add_gas)

             END IF

             qp_init(5+n_solid+n_add_gas:4+n_solid+n_add_gas+n_stoch_vars) =    &
                  1.0_wp
             
             qp_init(5+n_solid+n_add_gas+n_stoch_vars:4+n_solid+n_add_gas       &
                  +n_stoch_vars+n_pore_vars) = 1.0_wp             
             
             CALL qp_to_qc( qp_init(1:n_vars+2) , q(1:n_vars,j,k) )
             
          ELSE
             
             CALL qp_to_qc( qp0_init(1:n_vars+2) , q(1:n_vars,j,k) )
             
          END IF
          
       END DO
       
    END DO

    RETURN

  END SUBROUTINE collapsing_volume

END MODULE init_2d
