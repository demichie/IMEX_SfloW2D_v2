!********************************************************************************
!> \brief Initial solution
!
!> This module contains the variables and the subroutine for the
!> initialization of the solution for a Riemann problem.
!********************************************************************************

MODULE init_2d

  USE parameters_2d, ONLY : wp
  USE parameters_2d, ONLY : verbose_level
  USE parameters_2d, ONLY : n_solid

  IMPLICIT none

  REAL(wp), ALLOCATABLE :: q_init(:,:,:)

  REAL(wp), ALLOCATABLE :: thickness_init(:,:)
  
  !> Riemann problem interface relative position. It is a value
  !> between 0 and 1
  REAL(wp) :: riemann_interface  

  REAL(wp) :: hB_W         !< Left height
  REAL(wp) :: u_W          !< Left velocity x
  REAL(wp) :: v_W          !< Left velocity y
  REAL(wp),ALLOCATABLE :: alphas_W(:)         !< Left sediment concentration
  REAL(wp) :: T_W          !< Left temperature

  REAL(wp) :: hB_E         !< Right height
  REAL(wp) :: u_E          !< Right velocity x
  REAL(wp) :: v_E          !< Right velocity y
  REAL(wp),ALLOCATABLE :: alphas_E(:)         !< Right sediment concentration
  REAL(wp) :: T_E          !< Right temperature


CONTAINS

  !******************************************************************************
  !> \brief Riemann problem initialization
  !
  !> This subroutine initialize the solution for a Riemann problem. The 
  !> values for the left and right states and the interface location 
  !> are read from the input file.\
  !> \date 26/08/2011
  !******************************************************************************

  SUBROUTINE riemann_problem

    USE constitutive_2d, ONLY : qp_to_qc

    USE geometry_2d, ONLY : x_comp , comp_cells_x , comp_cells_y , B_cent

    USE parameters_2d, ONLY : n_vars , verbose_level , n_solid

    USE solver_2d, ONLY : q

    IMPLICIT none

    ! REAL(wp) :: hB            !< height + topography
    ! REAL(wp) :: u             !< velocity
    ! REAL(wp) :: v             !< velocity

    REAL(wp) :: qp(n_vars,comp_cells_x,comp_cells_y) , qj(n_vars)

    INTEGER :: j,k          !< loop counter
    INTEGER :: i1           !< last index with left state

    INTEGER :: i_solid

    REAL(wp) :: eps

    IF ( verbose_level .GE. 1 ) THEN

       WRITE(*,*) 'Riemann problem initialization'
       WRITE(*,*) 'x_comp(1)',x_comp(1)
       WRITE(*,*) 'x_comp(comp_cells_x)',x_comp(comp_cells_x)
       WRITE(*,*) 'riemann_interface',riemann_interface

    END IF

    i1 = 0
    
    riemann_int_search:DO j = 1,comp_cells_x
       
       IF ( x_comp(j) .LT. riemann_interface ) THEN

          i1 = j

       ELSE

          EXIT riemann_int_search

       END IF

    END DO riemann_int_search
    
    eps = 1.D-10

    ! Left initial state
    qp(1,1:i1,:) = hB_W
    qp(2,1:i1,:) = u_W
    qp(3,1:i1,:) = v_W
    qp(4,1:i1,:) = T_W

    ALLOCATE( alphas_W(n_solid) )

    DO i_solid=1,n_solid

       qp(4+i_solid,1:i1,:) = alphas_W(i_solid)

    END DO

    qp(n_vars+1,1:i1,:) = 0.0_wp
    qp(n_vars+2,1:i1,:) = 0.0_wp

    IF ( verbose_level .GE. 1 ) WRITE(*,*) 'Left state'

    DO j = 1,i1

       DO k = 1,comp_cells_y

         ! evaluate the vector of conservative variables
         CALL qp_to_qc( qp(:,j,k) , qj )

         q(1:n_vars,j,k) = qj

         IF ( verbose_level .GE. 1 ) THEN 
            
            WRITE(*,*) j,k,B_cent(j,k)
            WRITE(*,*) qp(:,j,k)
            WRITE(*,*) q(1:n_vars,j,k)

         END IF

       ENDDO

    END DO

    IF ( verbose_level .GE. 1 ) READ(*,*)

    ! Right initial state
    qp(1,i1+1:comp_cells_x,:) = hB_E
    qp(2,i1+1:comp_cells_x,:) = u_E
    qp(3,i1+1:comp_cells_x,:) = v_E
    qp(4,i1+1:comp_cells_x,:) = T_E

    ALLOCATE( alphas_E(n_solid) )

    DO i_solid=1,n_solid

       qp(4+i_solid,1:i1,:) = alphas_E(i_solid)

    END DO

    qp(n_vars+1,i1+1:comp_cells_x,:) = 0.0_wp
    qp(n_vars+2,i1+1:comp_cells_x,:) = 0.0_wp


    IF ( verbose_level .GE. 1 ) WRITE(*,*) 'Right state'

    DO j = i1+1,comp_cells_x

       DO k = 1,comp_cells_y

         ! evaluate the vector of conservative variables
         CALL qp_to_qc( qp(:,j,k) , qj )

         q(1:n_vars,j,k) = qj

         IF ( verbose_level .GE. 1 ) THEN 
            
            WRITE(*,*) j,k,B_cent(j,k)
            WRITE(*,*) qp(:,j,k)
            WRITE(*,*) q(1:n_vars,j,k)

         END IF
    
      END DO

    ENDDO

    IF ( verbose_level .GE. 1 ) READ(*,*)

    RETURN

  END SUBROUTINE riemann_problem

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

    USE parameters_2d, ONLY : n_vars

    USE parameters_2d, ONLY : x_collapse , y_collapse , r_collapse , T_collapse , &
       h_collapse , alphas_collapse

    USE solver_2d, ONLY : q

    IMPLICIT NONE

    INTEGER :: j,k

    REAL(wp) :: qp_init(n_vars+2) ,  qp0_init(n_vars+2)

    REAL(wp) :: cell_fract(comp_cells_x,comp_cells_y)

    CALL compute_cell_fract(x_collapse,y_collapse,r_collapse,cell_fract)

    ! values outside the collapsing volume
    qp0_init(1) = 0.0_wp                  ! h
    qp0_init(2) = 0.0_wp                  ! hu
    qp0_init(3) = 0.0_wp                  ! hv
    qp0_init(4) = T_collapse            ! T
    qp0_init(5:4+n_solid) = 0.0_wp        ! alphas
    qp0_init(n_vars+1:n_vars+2) = 0.0_wp  ! u,v

    ! values within the collapsing volume
    qp_init(1) = h_collapse
    qp_init(2) = 0.0_wp
    qp_init(3) = 0.0_wp
    qp_init(4) = T_collapse
    qp_init(5:4+n_solid) = alphas_collapse(1:n_solid)
    qp_init(n_vars+1:n_vars+2) = 0.0_wp
    
    DO j = 1,comp_cells_x
       
       DO k = 1,comp_cells_y
          
          IF ( cell_fract(j,k) .GT. 0.0_wp ) THEN
             
             qp_init(1) = cell_fract(j,k) * h_collapse

             CALL qp_to_qc( qp_init(1:n_vars+2) , q(1:n_vars,j,k) )
             
          ELSE
             
             CALL qp_to_qc( qp0_init(1:n_vars+2) , q(1:n_vars,j,k) )
             
          END IF
          
       END DO
       
    END DO

    RETURN

  END SUBROUTINE collapsing_volume

END MODULE init_2d
