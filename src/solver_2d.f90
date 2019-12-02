!********************************************************************************
!> \brief Numerical solver
!
!> This module contains the variables and the subroutines for the 
!> numerical solution of the equations.  
!
!> \date 07/10/2016
!> @author 
!> Mattia de' Michieli Vitturi
!
!********************************************************************************
MODULE solver_2d

  ! external variables

  USE constitutive_2d, ONLY : implicit_flag

  USE geometry_2d, ONLY : comp_cells_x,comp_cells_y,comp_cells_xy
  USE geometry_2d, ONLY : comp_interfaces_x,comp_interfaces_y

  USE geometry_2d, ONLY : B_cent , B_prime_x , B_prime_y
  USE geometry_2d, ONLY : B_interfaceL , B_interfaceR , B_interfaceT , B_interfaceB
  USE geometry_2d, ONLY : grav_surf
  USE geometry_2d, ONLY : source_cell

  USE parameters_2d, ONLY : n_eqns , n_vars , n_nh , n_solid
  USE parameters_2d, ONLY : n_RK
  USE parameters_2d, ONLY : verbose_level
  USE parameters_2d, ONLY : radial_source_flag

  USE parameters_2d, ONLY : bcW , bcE , bcS , bcN

  ! external procedures
  USE geometry_2d, ONLY : limit

  IMPLICIT none

  !> Conservative variables
  REAL*8, ALLOCATABLE :: q(:,:,:)        
  !> Conservative variables
  REAL*8, ALLOCATABLE :: q0(:,:,:)        
  !> Solution of the finite-volume semidiscrete cheme
  REAL*8, ALLOCATABLE :: q_fv(:,:,:)     

  !> Reconstructed value at the left of the x-interface
  REAL*8, ALLOCATABLE :: q_interfaceL(:,:,:)        
  !> Reconstructed value at the right of the x-interface
  REAL*8, ALLOCATABLE :: q_interfaceR(:,:,:)
  !> Reconstructed value at the bottom of the y-interface
  REAL*8, ALLOCATABLE :: q_interfaceB(:,:,:)        
  !> Reconstructed value at the top of the y-interface
  REAL*8, ALLOCATABLE :: q_interfaceT(:,:,:)

  !> Reconstructed physical value at the left of the x-interface
  REAL*8, ALLOCATABLE :: qp_interfaceL(:,:,:)        
  !> Reconstructed physical value at the right of the x-interface
  REAL*8, ALLOCATABLE :: qp_interfaceR(:,:,:)
  !> Reconstructed physical value at the bottom of the y-interface
  REAL*8, ALLOCATABLE :: qp_interfaceB(:,:,:)        
  !> Reconstructed physical value at the top of the y-interface
  REAL*8, ALLOCATABLE :: qp_interfaceT(:,:,:)

  !> Reconstructed value at the NW corner of cell
  REAL*8, ALLOCATABLE :: q_cellNW(:,:,:)        
  !> Reconstructed value at the NE corner of cell
  REAL*8, ALLOCATABLE :: q_cellNE(:,:,:)
  !> Reconstructed value at the SW corner of cell
  REAL*8, ALLOCATABLE :: q_cellSW(:,:,:)        
  !> Reconstructed value at the SE corner of cell
  REAL*8, ALLOCATABLE :: q_cellSE(:,:,:)

  !> Reconstructed physical value at the NW corner of cell
  REAL*8, ALLOCATABLE :: qp_cellNW(:,:,:)        
  !> Reconstructed physical value at the NE corner of cell
  REAL*8, ALLOCATABLE :: qp_cellNE(:,:,:)
  !> Reconstructed physical value at the SW corner of cell
  REAL*8, ALLOCATABLE :: qp_cellSW(:,:,:)        
  !> Reconstructed physical value at the SE corner of cell
  REAL*8, ALLOCATABLE :: qp_cellSE(:,:,:)


  !> Maximum over time of thickness
  REAL*8, ALLOCATABLE :: q1max(:,:)
  
  !> Local speeds at the left of the x-interface
  REAL*8, ALLOCATABLE :: a_interface_xNeg(:,:,:)
  !> Local speeds at the right of the x-interface
  REAL*8, ALLOCATABLE :: a_interface_xPos(:,:,:)
  !> Local speeds at the bottom of the y-interface
  REAL*8, ALLOCATABLE :: a_interface_yNeg(:,:,:)
  !> Local speeds at the top of the y-interface
  REAL*8, ALLOCATABLE :: a_interface_yPos(:,:,:)
  !> Semidiscrete numerical interface fluxes 
  REAL*8, ALLOCATABLE :: H_interface_x(:,:,:)
  !> Semidiscrete numerical interface fluxes 
  REAL*8, ALLOCATABLE :: H_interface_y(:,:,:)
  !> Physical variables (\f$\alpha_1, p_1, p_2, \rho u, w, T\f$)
  REAL*8, ALLOCATABLE :: qp(:,:,:)

  !> Array defining fraction of cells affected by source term
  REAL*8, ALLOCATABLE :: source_xy(:,:)

  LOGICAL, ALLOCATABLE :: solve_mask(:,:) , solve_mask0(:,:)

  !> Time step
  REAL*8 :: dt

  LOGICAL, ALLOCATABLE :: mask22(:,:) , mask21(:,:) , mask11(:,:) , mask12(:,:)

  INTEGER :: i_RK           !< loop counter for the RK iteration

  !> Butcher Tableau for the explicit part of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: a_tilde_ij(:,:)
  !> Butcher Tableau for the implicit part of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: a_dirk_ij(:,:)

  !> Coefficients for the explicit part of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: omega_tilde(:)

  !> Coefficients for the implicit part of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: omega(:)

  !> Explicit coeff. for the hyperbolic part for a single step of the R-K scheme
  REAL*8, ALLOCATABLE :: a_tilde(:)

  !> Explicit coeff. for the non-hyp. part for a single step of the R-K scheme
  REAL*8, ALLOCATABLE :: a_dirk(:)

  !> Implicit coeff. for the non-hyp. part for a single step of the R-K scheme
  REAL*8 :: a_diag

  !> Intermediate solutions of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: q_rk(:,:,:,:)

  !> Intermediate physical solutions of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: qp_rk(:,:,:,:)


  !> Intermediate hyperbolic terms of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: divFlux(:,:,:,:)

  !> Intermediate non-hyperbolic terms of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: NH(:,:,:,:)

  !> Intermediate semi-implicit non-hyperbolic terms of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: SI_NH(:,:,:,:)

  !> Intermediate explicit terms of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: expl_terms(:,:,:,:)

  !> Local Intermediate hyperbolic terms of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: divFluxj(:,:)

  !> Local Intermediate non-hyperbolic terms of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: NHj(:,:)

  !> Local Intermediate explicit terms of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: expl_terms_j(:,:)

  !> Local Intermediate semi-impl non-hyperbolic terms of the Runge-Kutta scheme
  REAL*8, ALLOCATABLE :: SI_NHj(:,:)

  !> Flag for the normalization of the array q in the implicit solution scheme
  LOGICAL :: normalize_q

  !> Flag for the normalization of the array f in the implicit solution scheme
  LOGICAL :: normalize_f

  !> Flag for the search of optimal step size in the implicit solution scheme
  LOGICAL :: opt_search_NL

  !> Sum of all the terms of the equations except the transient term
  REAL*8, ALLOCATABLE :: residual_term(:,:,:)

  INTEGER, ALLOCATABLE :: j_cent(:)
  INTEGER, ALLOCATABLE :: k_cent(:)

  INTEGER, ALLOCATABLE :: j_stag_x(:)
  INTEGER, ALLOCATABLE :: k_stag_x(:)

  INTEGER, ALLOCATABLE :: j_stag_y(:)
  INTEGER, ALLOCATABLE :: k_stag_y(:)

  REAL*8 :: t_imex1,t_imex2

CONTAINS

  !******************************************************************************
  !> \brief Memory allocation
  !
  !> This subroutine allocate the memory for the variables of the 
  !> solver module.
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE allocate_solver_variables

    IMPLICIT NONE

    REAL*8 :: gamma , delta

    INTEGER :: i,j,k

    ALLOCATE( q( n_vars , comp_cells_x , comp_cells_y ) , q0( n_vars ,          &
         comp_cells_x , comp_cells_y ) )

    ALLOCATE( q1max( comp_cells_x , comp_cells_y ) )
    
    ALLOCATE( q_fv( n_vars , comp_cells_x , comp_cells_y ) )

    ALLOCATE( q_interfaceL( n_vars , comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( q_interfaceR( n_vars , comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( a_interface_xNeg( n_eqns , comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( a_interface_xPos( n_eqns , comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( H_interface_x( n_eqns , comp_interfaces_x, comp_cells_y ) )


    ALLOCATE( q_interfaceB( n_vars , comp_cells_x, comp_interfaces_y ) )
    ALLOCATE( q_interfaceT( n_vars , comp_cells_x, comp_interfaces_y ) )
    ALLOCATE( a_interface_yNeg( n_eqns , comp_cells_x, comp_interfaces_y ) )
    ALLOCATE( a_interface_yPos( n_eqns , comp_cells_x, comp_interfaces_y ) )

    ALLOCATE( qp_interfaceL( n_vars , comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( qp_interfaceR( n_vars , comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( qp_interfaceB( n_vars , comp_cells_x, comp_interfaces_y ) )
    ALLOCATE( qp_interfaceT( n_vars , comp_cells_x, comp_interfaces_y ) )


    ALLOCATE( q_cellNW( n_vars , comp_cells_x , comp_cells_y ) )
    ALLOCATE( q_cellNE( n_vars , comp_cells_x , comp_cells_y ) )
    ALLOCATE( q_cellSW( n_vars , comp_cells_x , comp_cells_y ) )
    ALLOCATE( q_cellSE( n_vars , comp_cells_x , comp_cells_y ) )

    ALLOCATE( qp_cellNW( n_vars , comp_cells_x , comp_cells_y ) )
    ALLOCATE( qp_cellNE( n_vars , comp_cells_x , comp_cells_y ) )
    ALLOCATE( qp_cellSW( n_vars , comp_cells_x , comp_cells_y ) )
    ALLOCATE( qp_cellSE( n_vars , comp_cells_x , comp_cells_y ) )


    ALLOCATE( H_interface_y( n_eqns , comp_cells_x, comp_interfaces_y ) )

    ALLOCATE( solve_mask( comp_cells_x , comp_cells_y ) )
    ALLOCATE( solve_mask0( comp_cells_x , comp_cells_y ) )

    ALLOCATE( qp( n_vars , comp_cells_x , comp_cells_y ) )

    ALLOCATE( source_xy( comp_cells_x , comp_cells_y ) )

    
    ALLOCATE( a_tilde_ij(n_RK,n_RK) )
    ALLOCATE( a_dirk_ij(n_RK,n_RK) )
    ALLOCATE( omega_tilde(n_RK) )
    ALLOCATE( omega(n_RK) )


    ! Allocate the logical arrays defining the implicit parts of the system
    ALLOCATE( mask22(n_eqns,n_eqns) )
    ALLOCATE( mask21(n_eqns,n_eqns) )
    ALLOCATE( mask11(n_eqns,n_eqns) )
    ALLOCATE( mask12(n_eqns,n_eqns) )

    ! Initialize the logical arrays with all false (everything is implicit)
    mask11(1:n_eqns,1:n_eqns) = .FALSE.
    mask12(1:n_eqns,1:n_eqns) = .FALSE.
    mask22(1:n_eqns,1:n_eqns) = .FALSE.
    mask21(1:n_eqns,1:n_eqns) = .FALSE.

    ! Set to .TRUE. the elements not corresponding to equations and variables to 
    ! be solved implicitly
    DO i = 1,n_eqns

       DO j = 1,n_eqns

          IF ( .NOT.implicit_flag(i) .AND. .NOT.implicit_flag(j) )              &
               mask11(j,i) = .TRUE.
          IF ( implicit_flag(i) .AND. .NOT.implicit_flag(j) )                   &
               mask12(j,i) = .TRUE.
          IF ( implicit_flag(i) .AND. implicit_flag(j) )                        &
               mask22(j,i) = .TRUE.
          IF ( .NOT.implicit_flag(i) .AND. implicit_flag(j) )                   &
               mask21(j,i) = .TRUE.

       END DO

    END DO

    ! Initialize the coefficients for the IMEX Runge-Kutta scheme
    ! Please note that with respect to the schemes described in Pareschi & Russo 
    ! (2000) we do not have the coefficient vectors c_tilde and c, because the 
    ! explicit and implicit terms do not depend explicitly on time.

    ! Explicit part coefficients (a_tilde_ij=0 for j>=i)
    a_tilde_ij = 0.D0

    ! Weight coefficients of the explicit part in the final assemblage
    omega_tilde = 0.D0

    ! Implicit part coefficients (a_dirk_ij=0 for j>i)
    a_dirk_ij = 0.D0

    ! Weight coefficients of the explicit part in the final assemblage
    omega = 0.D0

    gamma = 1.D0 - 1.D0 / SQRT(2.D0)
    delta = 1.D0 - 1.D0 / ( 2.D0 * gamma )

    IF ( n_RK .EQ. 1 ) THEN

       a_tilde_ij(1,1) = 1.D0

       omega_tilde(1) = 1.D0

       a_dirk_ij(1,1) = 0.D0

       omega(1) = 0.D0

    ELSEIF ( n_RK .EQ. 2 ) THEN

       a_tilde_ij(2,1) = 1.0D0

       omega_tilde(1) = 1.0D0
       omega_tilde(2) = 0.0D0

       a_dirk_ij(2,2) = 1.0D0

       omega(1) = 0.D0
       omega(2) = 1.D0

      
    ELSEIF ( n_RK .EQ. 3 ) THEN

       ! Tableau for the IMEX-SSP(3,3,2) Stiffly Accurate Scheme
       ! from Pareschi & Russo (2005), Table IV

       a_tilde_ij(2,1) = 0.5D0
       a_tilde_ij(3,1) = 0.5D0
       a_tilde_ij(3,2) = 0.5D0

       omega_tilde(1) =  1.0D0 / 3.0D0
       omega_tilde(2) =  1.0D0 / 3.0D0
       omega_tilde(3) =  1.0D0 / 3.0D0

       a_dirk_ij(1,1) = 0.25D0
       a_dirk_ij(2,2) = 0.25D0
       a_dirk_ij(3,1) = 1.0D0 / 3.0D0
       a_dirk_ij(3,2) = 1.0D0 / 3.0D0
       a_dirk_ij(3,3) = 1.0D0 / 3.0D0

       omega(1) =  1.0D0 / 3.0D0
       omega(2) =  1.0D0 / 3.0D0
       omega(3) =  1.0D0 / 3.0D0
       
    ELSEIF ( n_RK .EQ. 4 ) THEN

       ! LRR(3,2,2) from Table 3 in Pareschi & Russo (2000)

       a_tilde_ij(2,1) = 0.5D0
       a_tilde_ij(3,1) = 1.D0 / 3.D0
       a_tilde_ij(4,2) = 1.0D0

       omega_tilde(1) = 0.D0
       omega_tilde(2) = 1.0D0
       omega_tilde(3) = 0.0D0
       omega_tilde(4) = 0.D0

       a_dirk_ij(2,2) = 0.5D0
       a_dirk_ij(3,3) = 1.0D0 / 3.0D0
       a_dirk_ij(4,3) = 0.75D0
       a_dirk_ij(4,4) = 0.25D0

       omega(1) = 0.D0
       omega(2) = 0.D0
       omega(3) = 0.75D0
       omega(4) = 0.25D0

    END IF

    ALLOCATE( a_tilde(n_RK) )
    ALLOCATE( a_dirk(n_RK) )

    ALLOCATE( q_rk( n_vars , comp_cells_x , comp_cells_y , n_RK ) )
    ALLOCATE( qp_rk( n_vars , comp_cells_x , comp_cells_y , n_RK ) )
    ALLOCATE( divFlux( n_eqns , comp_cells_x , comp_cells_y , n_RK ) )
    ALLOCATE( NH( n_eqns , comp_cells_x , comp_cells_y , n_RK ) )
    ALLOCATE( SI_NH( n_eqns , comp_cells_x , comp_cells_y , n_RK ) )

    ALLOCATE( expl_terms( n_eqns , comp_cells_x , comp_cells_y , n_RK ) )

    ALLOCATE( divFluxj(n_eqns,n_RK) )
    ALLOCATE( NHj(n_eqns,n_RK) )
    ALLOCATE( SI_NHj(n_eqns,n_RK) )

    ALLOCATE( expl_terms_j(n_eqns,n_RK) )

    ALLOCATE( residual_term( n_vars , comp_cells_x , comp_cells_y ) )

    comp_cells_xy = comp_cells_x * comp_cells_y

    ALLOCATE( j_cent( comp_cells_xy ) )
    ALLOCATE( k_cent( comp_cells_xy ) )

    DO j = 1,comp_cells_x

       DO k = 1,comp_cells_y

          j_cent( comp_cells_y*(j-1)+k ) = j
          k_cent( comp_cells_y*(j-1)+k ) = k

       END DO

    END DO

    ALLOCATE( j_stag_x( comp_interfaces_x * comp_cells_y ) )
    ALLOCATE( k_stag_x( comp_interfaces_x * comp_cells_y ) )

    DO j = 1,comp_interfaces_x

       DO k = 1,comp_cells_y

          j_stag_x( comp_cells_y*(j-1)+k ) = j
          k_stag_x( comp_cells_y*(j-1)+k ) = k

       END DO

    END DO


    ALLOCATE( j_stag_y( comp_cells_x * comp_interfaces_y ) )
    ALLOCATE( k_stag_y( comp_cells_x * comp_interfaces_y ) )

    DO j = 1,comp_cells_x

       DO k = 1,comp_interfaces_y

          j_stag_x( comp_interfaces_y*(j-1)+k ) = j
          k_stag_x( comp_interfaces_y*(j-1)+k ) = k

       END DO

    END DO

  END SUBROUTINE allocate_solver_variables

  !******************************************************************************
  !> \brief Memory deallocation
  !
  !> This subroutine de-allocate the memory for the variables of the 
  !> solver module.
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE deallocate_solver_variables

    DEALLOCATE( q , q0 )

    DEALLOCATE( q1max )
    
    DEALLOCATE( q_fv )

    DEALLOCATE( q_interfaceL )
    DEALLOCATE( q_interfaceR )
    DEALLOCATE( q_interfaceB )
    DEALLOCATE( q_interfaceT )

    DEALLOCATE( q_cellNW )
    DEALLOCATE( q_cellNE )
    DEALLOCATE( q_cellSW )
    DEALLOCATE( q_cellSE )
  
    DEALLOCATE( qp_interfaceL )
    DEALLOCATE( qp_interfaceR )
    DEALLOCATE( qp_interfaceB )
    DEALLOCATE( qp_interfaceT )

    DEALLOCATE( qp_cellNW )
    DEALLOCATE( qp_cellNE )
    DEALLOCATE( qp_cellSW )
    DEALLOCATE( qp_cellSE )

  
    DEALLOCATE( a_interface_xNeg )
    DEALLOCATE( a_interface_xPos )
    DEALLOCATE( a_interface_yNeg )
    DEALLOCATE( a_interface_yPos )

    DEALLOCATE( H_interface_x )
    DEALLOCATE( H_interface_y )

    DEALLOCATE( solve_mask , solve_mask0 )

    Deallocate( qp )

    DEALLOCATE( source_xy )
    
    DEALLOCATE( a_tilde_ij )
    DEALLOCATE( a_dirk_ij )
    DEALLOCATE( omega_tilde )
    DEALLOCATE( omega )

    DEALLOCATE( implicit_flag )

    DEALLOCATE( a_tilde )
    DEALLOCATE( a_dirk )

    DEALLOCATE( q_rk )
    DEALLOCATE( qp_rk )
    DEALLOCATE( divFlux )
    DEALLOCATE( NH )
    DEALLOCATE( SI_NH )
    DEALLOCATE( expl_terms )

    DEALLOCATE( divFluxj )
    DEALLOCATE( NHj )
    DEALLOCATE( SI_NHj )
    DEALLOCATE( expl_terms_j )

    DEALLOCATE( mask22 , mask21 , mask11 , mask12 )

    DEALLOCATE( residual_term )

    DEALLOCATE( j_cent )
    ! DEALLOCATE( k_cent )
    ! DEALLOCATE( j_cent , k_cent )
    ! DEALLOCATE ( j_stag_x , k_stag_x )
    ! DEALLOCATE ( j_stag_y , k_stag_y )

  END SUBROUTINE deallocate_solver_variables


  !******************************************************************************
  !> \brief Masking of cells to solve
  !
  !> This subroutine compute a 2D array of logicals defining the cells where the
  !> systems of equations have to be solved. It is defined according to the 
  !> positive thickness in the cell and in the neighbour cells
  !
  !> \date 20/04/2017
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE check_solve

    IMPLICIT NONE

    INTEGER :: i,j,k

    solve_mask0(1:comp_cells_x,1:comp_cells_y) = .FALSE.

    WHERE ( q(1,:,:) .GT. 0.D0 ) solve_mask0 = .TRUE.
    
    solve_mask = solve_mask0

    DO i = 1,n_RK

       solve_mask(1+i:comp_cells_x,:) = solve_mask(1+i:comp_cells_x,:) .OR.     &
            solve_mask(1:comp_cells_x-i,:) 

       solve_mask(1:comp_cells_x-i,:) = solve_mask(1:comp_cells_x-i,:) .OR.     &
            solve_mask(1+i:comp_cells_x,:) 

       solve_mask(:,1+i:comp_cells_y) = solve_mask(:,1+i:comp_cells_y) .OR.     &
            solve_mask(:,1:comp_cells_y-i) 

       solve_mask(:,1:comp_cells_y-i) = solve_mask(:,1:comp_cells_y-i) .OR.     &
            solve_mask(:,1+i:comp_cells_y) 

    END DO

    solve_mask(1:comp_cells_x,1:comp_cells_y) = .TRUE.
    
    IF ( radial_source_flag ) THEN

       DO j = 1,comp_cells_x

          DO k = 1,comp_cells_y

             IF ( source_cell(j,k) .EQ. 1 ) solve_mask(j,k) = .FALSE.

          END DO

       END DO

    END IF

    RETURN

  END SUBROUTINE check_solve

  !*****************************************************************************
  !> \brief Time-step computation
  !
  !> This subroutine evaluate the maximum time step according to the CFL
  !> condition. The local speed are evaluated with the characteristic
  !> polynomial of the Jacobian of the fluxes.
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !*****************************************************************************

  SUBROUTINE timestep

    ! External variables
    USE geometry_2d, ONLY : dx,dy
    USE parameters_2d, ONLY : max_dt , cfl

    USE constitutive_2d, ONLY : qc_to_qp

    IMPLICIT none

    REAL*8 :: dt_cfl        !< local time step

    REAL*8 :: a_interface_x_max(n_eqns,comp_interfaces_x,comp_cells_y)
    REAL*8 :: a_interface_y_max(n_eqns,comp_cells_x,comp_interfaces_y)
    REAL*8 :: dt_interface_x, dt_interface_y

    INTEGER :: i,j,k          !< loop counter

    REAL*8 :: max_a
    
    REAL*8 :: max_vel
    INTEGER :: max_j,max_k
    INTEGER :: max_dir
    
    max_vel = 0.D0

    dt = max_dt

    IF ( cfl .NE. -1.d0 ) THEN

       ! Compute the physical variables at the cell centers
       DO j = 1,comp_cells_x

          DO k = 1,comp_cells_y

             CALL qc_to_qp( q(1:n_vars,j,k), B_cent(j,k) , qp(1:n_vars,j,k) )

          END DO

       END DO

       ! Compute the physical and conservative variables at the interfaces
       CALL reconstruction( q , qp )

       ! Compute the max/min eigenvalues at the interfaces
       CALL eval_speeds

       DO i=1,n_vars

          a_interface_x_max(i,:,:) =                                            &
               MAX( a_interface_xPos(i,:,:) , -a_interface_xNeg(i,:,:) )

          a_interface_y_max(i,:,:) =                                            &
               MAX( a_interface_yPos(i,:,:) , -a_interface_yNeg(i,:,:) )

       END DO

       DO j = 1,comp_cells_x

          DO k = 1,comp_cells_y

             max_a =  MAX( MAXVAL(a_interface_x_max(:,j,k)) ,                   &
                  MAXVAL(a_interface_x_max(:,j+1,k)) )

             IF ( max_a .GT. max_vel) THEN 

                max_vel = max_a
                max_j = j
                max_k = k
                max_dir = 1

             END IF
                
             IF ( max_a .GT. 0.D0 ) THEN

                dt_interface_x = cfl * dx / max_a

             ELSE

                dt_interface_x = dt
                
             END IF
             
             max_a =  MAX( MAXVAL(a_interface_y_max(:,j,k)) ,                   &
                  MAXVAL(a_interface_y_max(:,j,k+1)) )

             IF ( max_a .GT. max_vel) THEN 

                max_vel = max_a
                max_j = j
                max_k = k
                max_dir = 2

             END IF
             
             IF ( max_a .GT. 0.D0 ) THEN

                dt_interface_y = cfl * dy / max_a

             ELSE

                dt_interface_y = dt

             END IF

             
             dt_cfl = MIN( dt_interface_x , dt_interface_y )

             dt = MIN(dt,dt_cfl)

!!$             IF ( ( j.EQ. 230 ) .AND. ( k.EQ.123) ) THEN
!!$
!!$                WRITE(*,*) a_interface_x_max(1,j,k),a_interface_x_max(1,j+1,k)
!!$                
!!$                WRITE(*,*) qp(1,j,k),qp_interfaceR(1,j+1,k),qp(1,j+1,k),  &
!!$                     qp_interfaceL(1,j+2,k),qp(1,j+2,k)
!!$                
!!$                WRITE(*,*) qp(2,j,k),qp_interfaceR(2,j+1,k),qp(2,j+1,k),  &
!!$                     qp_interfaceL(2,j+2,k),qp(2,j+2,k)
!!$                   
!!$             END IF

             
          ENDDO

       END DO

    END IF

!!$    WRITE(*,*) 'max_vel',max_vel,max_j,max_k,max_dir
    
  END SUBROUTINE timestep

  !******************************************************************************
  !> \brief Runge-Kutta integration
  !
  !> This subroutine integrate the hyperbolic conservation law with
  !> non-hyperbolic terms using an implicit-explicit runge-kutta scheme.
  !> The fluxes are integrated explicitely while the non-hyperbolic terms
  !> are integrated implicitely.
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE imex_RK_solver

    USE constitutive_2d, ONLY : eval_nonhyperbolic_terms

    USE constitutive_2d, ONLY : eval_nh_semi_impl_terms

    USE constitutive_2d, ONLY : qc_to_qp

    IMPLICIT NONE

    REAL*8 :: q_si(n_vars) !< solution after the semi-implicit step
    REAL*8 :: q_guess(n_vars) !< initial guess for the solution of the RK step
    INTEGER :: j,k            !< loop counter over the grid volumes

    REAL*8 :: h_new
    
    ! Initialization of the solution guess
    q0( 1:n_vars , 1:comp_cells_x , 1:comp_cells_y ) =                          &
         q( 1:n_vars , 1:comp_cells_x , 1:comp_cells_y )

    IF ( verbose_level .GE. 2 ) WRITE(*,*) 'solver, imex_RK_solver: beginning'

    ! Initialization of the variables for the Runge-Kutta scheme
    q_rk(1:n_vars,1:comp_cells_x,1:comp_cells_y,1:n_RK) = 0.d0
    qp_rk(1:n_vars,1:comp_cells_x,1:comp_cells_y,1:n_RK) = 0.d0

    divFlux(1:n_eqns,1:comp_cells_x,1:comp_cells_y,1:n_RK) = 0.d0

    NH(1:n_eqns,1:comp_cells_x,1:comp_cells_y,1:n_RK) = 0.d0

    SI_NH(1:n_eqns,1:comp_cells_x,1:comp_cells_y,1:n_RK) = 0.d0

    expl_terms(1:n_eqns,1:comp_cells_x,1:comp_cells_y,1:n_RK) = 0.d0

    runge_kutta:DO i_RK = 1,n_RK

       IF ( verbose_level .GE. 2 ) WRITE(*,*) 'solver, imex_RK_solver: i_RK',i_RK

       ! define the explicits coefficients for the i-th step of the Runge-Kutta
       a_tilde = 0.d0
       a_dirk = 0.d0

       ! in the first step of the RK scheme all the coefficients remain to 0
       a_tilde(1:i_RK-1) = a_tilde_ij(i_RK,1:i_RK-1)
       a_dirk(1:i_RK-1) = a_dirk_ij(i_RK,1:i_RK-1)

       ! define the implicit coefficient for the i-th step of the Runge-Kutta
       a_diag = a_dirk_ij(i_RK,i_RK)

       CALL cpu_time(t_imex1)

       loop_over_ycells:DO k = 1,comp_cells_y
       
          loop_over_xcells:DO j = 1,comp_cells_x

             IF ( solve_mask(j,k) ) THEN

             IF ( verbose_level .GE. 2 ) THEN

                WRITE(*,*) 'solver, imex_RK_solver: j',j,k

             END IF

             ! initialize the RK step
             IF ( i_RK .EQ. 1 ) THEN

                ! solution from the previous time step
                q_guess(1:n_vars) = q0( 1:n_vars , j , k) 

             ELSE

                ! solution from the previous RK step
                !q_guess(1:n_vars) = q_rk( 1:n_vars , j , k , MAX(1,i_RK-1) )

             END IF

             ! RK explicit hyperbolic terms for the volume (i,j)
             divFluxj(1:n_eqns,1:n_RK) = divFlux( 1:n_eqns , j , k , 1:n_RK )

             ! RK implicit terms for the volume (i,j)
             NHj(1:n_eqns,1:n_RK) = NH( 1:n_eqns , j , k , 1:n_RK )

             ! RK semi-implicit terms for the volume (i,j)
             SI_NHj(1:n_eqns,1:n_RK) = SI_NH( 1:n_eqns , j , k , 1:n_RK )
             
             ! RK additional explicit terms for the volume (i,j)
             Expl_terms_j(1:n_eqns,1:n_RK) = expl_terms( 1:n_eqns,j,k,1:n_RK )

             ! New solution at the i_RK step without the implicit  and
             ! semi-implicit term
             q_fv( 1:n_vars , j , k ) = q0( 1:n_vars , j , k )                  &
                  - dt * (MATMUL( divFluxj(1:n_eqns,1:i_RK)                     &
                  + Expl_terms_j(1:n_eqns,1:i_RK) , a_tilde(1:i_RK) )           &
                  - MATMUL( NHj(1:n_eqns,1:i_RK) + SI_NHj(1:n_eqns,1:i_RK) ,    &
                  a_dirk(1:i_RK) ) )

             IF ( verbose_level .GE. 2 ) THEN

                WRITE(*,*) 'q_guess',q_guess
                CALL qc_to_qp( q_guess , B_cent(j,k) , qp(1:n_vars,j,k) )
                WRITE(*,*) 'q_guess: qp',qp(1:n_vars,j,k)

             END IF

             adiag_pos:IF ( a_diag .NE. 0.D0 ) THEN

                pos_thick:IF ( q_fv(1,j,k) .GT.  0.D0 )  THEN

                   ! Eval the semi-implicit terms
                   ! (terms which non depend on velocity magnitude)
                   CALL eval_nh_semi_impl_terms( grav_surf(j,k) ,               &
                        r_qj = q_fv( 1:n_vars , j , k ) ,                       &
                        r_nh_semi_impl_term = SI_NH(1:n_eqns,j,k,i_RK) ) 

                   SI_NHj(1:n_eqns,i_RK) = SI_NH( 1:n_eqns,j,k,i_RK )

                   ! Assemble the initial guess for the implicit solver
                   q_si(1:n_vars) = q_fv(1:n_vars,j,k ) + dt * a_diag *         &
                        SI_NH(1:n_eqns,j,k,i_RK)
                   
                   IF ( q_fv(2,j,k)**2 + q_fv(3,j,k)**2 .EQ. 0.D0 ) THEN
                      
                      !Case 1: if the velocity was null, then it must stay null
                      q_si(2:3) = 0.D0 
                      
                   ELSEIF ( ( q_si(2)*q_fv(2,j,k) .LT. 0.D0 ) .OR.              &
                        ( q_si(3)*q_fv(3,j,k) .LT. 0.D0 ) ) THEN
                      
                      ! If the semi-impl. friction term changed the sign of the 
                      ! velocity then set it to zero
                      q_si(2:3) = 0.D0 

                   ELSE
                      
                      ! Align the velocity vector with previous one
                      q_si(2:3) = DSQRT( q_si(2)**2 + q_si(3)**2 ) *            &
                           q_fv(2:3,j,k) / DSQRT( q_fv(2,j,k)**2                &
                           + q_fv(3,j,k)**2 ) 
                      
                   END IF

                   ! Update the semi-implicit term accordingly with the
                   ! corrections above
                   SI_NH(1:n_eqns,j,k,i_RK) = ( q_si(1:n_vars) -                &
                        q_fv(1:n_vars,j,k ) ) / ( dt*a_diag )

                   SI_NHj(1:n_eqns,i_RK) = SI_NH( 1:n_eqns,j,k,i_RK )

                   ! Initialize the guess for the NR solver
                   q_guess(1:n_vars) = q_si(1:n_vars)

                   ! Solve the implicit system to find the solution at the 
                   ! i_RK step of the IMEX RK procedure
                   CALL solve_rk_step( B_cent(j,k) , q_guess(1:n_vars) ,        &
                        q0(1:n_vars,j,k ) , a_tilde , a_dirk , a_diag )

                   IF ( comp_cells_y .EQ. 1 ) THEN

                      q_guess(3) = 0.D0

                   END IF

                   IF ( comp_cells_x .EQ. 1 ) THEN

                      q_guess(2) = 0.D0

                   END IF

                   ! Eval and store the implicit term at the i_RK step
                   CALL eval_nonhyperbolic_terms( r_qj =q_guess ,               &
                        r_nh_term_impl = NH(1:n_eqns,j,k,i_RK) )
                   
                   IF ( q_si(2)**2 + q_si(3)**2 .EQ. 0.D0 ) THEN
                      
                      q_guess(2:3) = 0.D0 
                      
                   ELSEIF ( ( q_guess(2)*q_si(2) .LE. 0.D0 ) .AND.              &
                        ( q_guess(3)*q_si(3) .LE. 0.D0 ) ) THEN
                      
                      ! If the impl. friction term changed the sign of the 
                      ! velocity then set it to zero
                      q_guess(2:3) = 0.D0 
                      
                   ELSE
                      
                      ! Align the velocity vector with previous one
                      q_guess(2:3) = DSQRT( q_guess(2)**2 + q_guess(3)**2 ) *   &
                           q_si(2:3) / DSQRT( q_si(2)**2 + q_si(3)**2 ) 
                                            
                   END IF
                                      
                ELSE

                   ! If h=0 nothing has to be changed 
                   q_guess(1:n_vars) = q_fv( 1:n_vars , j , k ) 
                   q_si(1:n_vars) = q_fv( 1:n_vars , j , k ) 
                   SI_NH(1:n_eqns,j,k,i_RK) = 0.D0
                   NH(1:n_eqns,j,k,i_RK) = 0.D0

                END IF pos_thick

             END IF adiag_pos

             ! Check the sign of the flow thickness
             h_new = q_guess(1)

             IF ( a_diag .NE. 0.D0 ) THEN
                
                ! Update the implicit term with correction on the new velocity
                NH(1:n_vars,j,k,i_RK) = ( q_guess(1:n_vars) - q_si(1:n_vars))   &
                     / ( dt*a_diag ) 

             END IF

             ! Store the solution at the end of the i_RK step
             q_rk( 1:n_vars , j , k , i_RK ) = q_guess

             IF ( verbose_level .GE. 2 ) THEN

                WRITE(*,*) 'imex_RK_solver: qc',q_guess
                CALL qc_to_qp( q_guess, B_cent(j,k) , qp(1:n_vars,j,k) )
                WRITE(*,*) 'imex_RK_solver: qp',qp(1:n_vars,j,k)
                READ(*,*)

             END IF

          END IF

          END DO loop_over_xcells

       ENDDO loop_over_ycells

       CALL cpu_time(t_imex2)
       ! WRITE(*,*) 'Time taken by implicit',t_imex2-t_imex1,'seconds'

       IF ( omega_tilde(i_RK) .GT. 0.D0 ) THEN

          ! Compute the physical variables at the cell centers 
          DO k = 1,comp_cells_y
          
             DO j = 1,comp_cells_x

                CALL qc_to_qp( q_rk(1:n_vars,j,k,i_RK) , B_cent(j,k) ,          &
                     qp_rk(1:n_vars,j,k,i_RK) )

             END DO

          END DO

          ! Eval and store the explicit hyperbolic (fluxes) terms
          CALL eval_hyperbolic_terms(                                           &
               q_rk(1:n_vars,1:comp_cells_x,1:comp_cells_y,i_RK) ,              &
               qp_rk(1:n_vars,1:comp_cells_x,1:comp_cells_y,i_RK) ,             &
               divFlux(1:n_eqns,1:comp_cells_x,1:comp_cells_y,i_RK) )

          CALL cpu_time(t_imex1)
          ! WRITE(*,*) 'Time taken by explicit',t_imex1-t_imex2,'seconds'
          
          ! Eval and store the other explicit terms (e.g. gravity or viscous 
          ! forces)
          CALL eval_explicit_terms(                                             &
               q_rk(1:n_vars,1:comp_cells_x,1:comp_cells_y,i_RK) ,              &
               qp_rk(1:n_vars,1:comp_cells_x,1:comp_cells_y,i_RK) ,             &
               expl_terms(1:n_eqns,1:comp_cells_x,1:comp_cells_y,i_RK) )

          CALL cpu_time(t_imex1)
          ! WRITE(*,*) 'Time taken by explicit',t_imex1-t_imex2,'seconds'
          
       END IF
          
       IF ( verbose_level .GE. 1 ) THEN

          WRITE(*,*) 'div_flux(2),div_flux(3),expl_terms(2),expl_terms(3)'

          DO k = 1,comp_cells_y
          
             DO j = 1,comp_cells_x

                WRITE(*,*) divFlux(2,j,k,i_RK) , divFlux(3,j,k,i_RK) ,          &
                     expl_terms(2,j,k,i_RK) , expl_terms(3,j,k,i_RK)

             ENDDO

          END DO

          READ(*,*)

       END IF

    END DO runge_kutta

    DO k = 1,comp_cells_y
       
       DO j = 1,comp_cells_x
          
          residual_term(1:n_vars,j,k) = MATMUL( divFlux(1:n_eqns,j,k,1:n_RK)    &
               + expl_terms(1:n_eqns,j,k,1:n_RK) , omega_tilde ) -              &
               MATMUL( NH(1:n_eqns,j,k,1:n_RK) + SI_NH(1:n_eqns,j,k,1:n_RK) ,   &
               omega )
          
       ENDDO
       
    END DO
    
    assemble_sol_loop_x:DO k = 1,comp_cells_y
       
       assemble_sol_loop_y:DO j = 1,comp_cells_x
     
          ! IF ( ( .NOT. radial_source_flag ) .OR. ( source_cell(j,k) .NE. 1 ) ) THEN

          IF ( solve_mask(j,k) ) THEN
     
          IF ( verbose_level .GE. 1 ) THEN

             WRITE(*,*) 'cell jk =',j,k
             WRITE(*,*) 'before imex_RK_solver: qc',q0(1:n_vars,j,k)
             CALL qc_to_qp(q0(1:n_vars,j,k) , B_cent(j,k) , qp(1:n_vars,j,k))
             WRITE(*,*) 'before imex_RK_solver: qp',qp(1:n_vars,j,k)

          END IF

          IF ( ( SUM(ABS( omega_tilde(:)-a_tilde_ij(n_RK,:))) .EQ. 0.D0  )      &
               .AND. ( SUM(ABS(omega(:)-a_dirk_ij(n_RK,:))) .EQ. 0.D0 ) ) THEN

             ! The assembling coeffs are equal to the last step of the RK scheme
             q(1:n_vars,j,k) = q_rk(1:n_vars,j,k,n_RK)

          ELSE

             ! The assembling coeffs are different
             q(1:n_vars,j,k) = q0(1:n_vars,j,k) - dt*residual_term(1:n_vars,j,k)

          END IF

          negative_thickness_check:IF ( q(1,j,k) .LT. 0.D0 ) THEN

             IF ( q(1,j,k) .GT. -1.D-7 ) THEN
                
                q(1,j,k) = 0.D0
                q(2:n_vars,j,k) = 0.D0
                
             ELSE
                
                WRITE(*,*) 'j,k,n_RK',j,k,n_RK
                WRITE(*,*) 'dt',dt
                WRITE(*,*) 'source_cell',source_cell(j,k)
                WRITE(*,*) 'source_cell left',source_cell(j-1,k)
                WRITE(*,*) 'source_cell right',source_cell(j+1,k)
                WRITE(*,*) 'source_cell top',source_cell(j,k+1)
                WRITE(*,*) 'source_cell bottom',source_cell(j,k-1)
                WRITE(*,*) 'before imex_RK_solver: qc',q0(1:n_vars,j,k)
                CALL qc_to_qp(q0(1:n_vars,j,k) , B_cent(j,k) , qp(1:n_vars,j,k))
                WRITE(*,*) 'before imex_RK_solver: qp',qp(1:n_vars,j,k)
                WRITE(*,*) 'h old',q0(1,j,k)
                WRITE(*,*) 'h new',q(1,j,k)
                WRITE(*,*) 'B_cent(j,k)',B_cent(j,k)

                WRITE(*,*) 'hS',q_interfaceT(1,j,k)
                WRITE(*,*) 'hE',q_interfaceR(1,j,k)
                
                READ(*,*)
                
             END IF

          END IF negative_thickness_check

          IF ( SUM(q(5:4+n_solid,j,k)) .GT. q(1,j,k) ) THEN
             
             IF ( SUM(q(5:4+n_solid,j,k))-q(1,j,k) .LT. 1.D-10 ) THEN
                
                q(5:4+n_solid,j,k) = q(5:4+n_solid,j,k)                         &
                     / SUM(q(5:4+n_solid,j,k)) * q(1,j,k)
                
             ELSE
                
                WRITE(*,*) 'j,k,n_RK',j,k,n_RK
                WRITE(*,*) 'dt',dt
                WRITE(*,*) ' B_cent(j,k)', B_cent(j,k)
                
                WRITE(*,*) 'before imex_RK_solver: qc',q0(1:n_vars,j,k)
                CALL qc_to_qp(q0(1:n_vars,j,k) , B_cent(j,k) , qp(1:n_vars,j,k))
                WRITE(*,*) 'before imex_RK_solver: qp',qp(1:n_vars,j,k)
                
                CALL qc_to_qp(q(1:n_vars,j,k) , B_cent(j,k) , qp(1:n_vars,j,k))
                
                WRITE(*,*) 'after imex_RK_solver: qc',q(1:n_vars,j,k)
                WRITE(*,*) 'after imex_RK_solver: qp',qp(1:n_vars,j,k)
                
                
                WRITE(*,*) 'h old',q0(1,j,k)
                WRITE(*,*) 'h new',q(1,j,k)
                WRITE(*,*) 'alphas old', q0(5:4+n_solid,j,k) / q0(1,j,k)
                WRITE(*,*) 'alphas new', q(5:4+n_solid,j,k) / q(1,j,k) 
                READ(*,*)
                
             END IF
             
             IF ( verbose_level .GE. 1 ) THEN
                
                WRITE(*,*) 'h new',q(1,j,k) 
                READ(*,*)
                
             END IF
             
          END IF

       END IF

       ! END IF

       ENDDO assemble_sol_loop_y

    END DO assemble_sol_loop_x

    RETURN
    
  END SUBROUTINE imex_RK_solver

  !******************************************************************************
  !> \brief Runge-Kutta single step integration
  !
  !> This subroutine find the solution of the non-linear system 
  !> given the a step of the implicit-explicit Runge-Kutta scheme for a
  !> cell:\n
  !> \f$ Q^{(i)} = Q^n - dt \sum_{j=1}^{i-1}\tilde{a}_{j}\partial_x 
  !> F(Q^{(j)}) +  dt \sum_{j=1}^{i-1} a_j  NH(Q^{(j)}) 
  !> + dt a_{diag} NH(Q^{(i)}) \f$\n
  !
  !> \param[in]     Bj        topography at the cell center
  !> \param[in]     Bprimej   topography slope at the cell center
  !> \param[in,out] qj        conservative variables 
  !> \param[in]     qj_old    conservative variables at the old time step
  !> \param[in]     a_tilde   explicit coefficents for the fluxes
  !> \param[in]     a_dirk    explicit coefficient for the non-hyperbolic terms
  !> \param[in]     a_diag    implicit coefficient for the non-hyperbolic terms 
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE solve_rk_step( Bj, qj, qj_old, a_tilde , a_dirk , a_diag )

    USE parameters_2d, ONLY : max_nl_iter , tol_rel , tol_abs

    USE constitutive_2d, ONLY : qc_to_qp

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: Bj
    REAL*8, INTENT(INOUT) :: qj(n_vars)
    REAL*8, INTENT(IN) :: qj_old(n_vars)
    REAL*8, INTENT(IN) :: a_tilde(n_RK)
    REAL*8, INTENT(IN) :: a_dirk(n_RK)
    REAL*8, INTENT(IN) :: a_diag

    REAL*8 :: qj_init(n_vars)

    REAL*8 :: qj_org(n_vars) , qj_rel(n_vars)

    REAL*8 :: left_matrix(n_eqns,n_vars)
    REAL*8 :: right_term(n_eqns)

    REAL*8 :: scal_f

    REAL*8 :: coeff_f(n_eqns)

    REAL*8 :: qj_rel_NR_old(n_vars)
    REAL*8 :: scal_f_old
    REAL*8 :: desc_dir(n_vars)
    REAL*8 :: grad_f(n_vars)
    REAL*8 :: mod_desc_dir

    INTEGER :: pivot(n_vars)

    REAL*8 :: left_matrix_small22(n_nh,n_nh)
    REAL*8 :: left_matrix_small21(n_eqns-n_nh,n_nh)
    REAL*8 :: left_matrix_small11(n_eqns-n_nh,n_vars-n_nh)
    REAL*8 :: left_matrix_small12(n_nh,n_vars-n_nh)

    REAL*8 :: desc_dir_small2(n_nh)
    INTEGER :: pivot_small2(n_nh)

    REAL*8 :: desc_dir_small1(n_vars-n_nh)

    INTEGER :: ok

    INTEGER :: i 
    INTEGER :: nl_iter

    REAL*8, PARAMETER :: STPMX=100.D0
    REAL*8 :: stpmax
    LOGICAL :: check

    REAL*8, PARAMETER :: TOLF=1.D-10 , TOLMIN=1.D-6
    REAL*8 :: TOLX

    REAL*8 :: qpj(n_vars)

    REAL*8 :: desc_dir2(n_vars)

    REAL*8 :: desc_dir_temp(n_vars)

    normalize_q = .TRUE.
    normalize_f = .FALSE.
    opt_search_NL = .TRUE.

    coeff_f(1:n_eqns) = 1.D0

    qj_init = qj

    ! normalize the functions of the nonlinear system
    IF ( normalize_f ) THEN

       qj = qj_old - dt * ( MATMUL(divFluxj+ Expl_terms_j,a_tilde)              &
            - MATMUL(NHj,a_dirk) )

       CALL eval_f( qj , qj_old , a_tilde , a_dirk , a_diag , coeff_f ,         &
            right_term ,  scal_f )

       IF ( verbose_level .GE. 3 ) THEN

          WRITE(*,*) 'solve_rk_step: non-normalized right_term'
          WRITE(*,*) right_term
          WRITE(*,*) 'scal_f',scal_f

       END IF

       DO i=1,n_eqns

          IF ( ABS(right_term(i)) .GE. 1.D0 ) coeff_f(i) = 1.D0 / right_term(i)

       END DO

       right_term = coeff_f * right_term

       scal_f = 0.5D0 * DOT_PRODUCT( right_term , right_term )

       IF ( verbose_level .GE. 3 ) THEN                    
          WRITE(*,*) 'solve_rk_step: after normalization',scal_f
       END IF

    END IF

    !---- normalize the conservative variables ------

    IF ( normalize_q ) THEN

       qj_org = qj

       qj_org = MAX( ABS(qj_org) , 1.D-3 )

    ELSE 

       qj_org(1:n_vars) = 1.D0

    END IF

    qj_rel = qj / qj_org

    ! -----------------------------------------------

    newton_raphson_loop:DO nl_iter=1,max_nl_iter

       TOLX = epsilon(qj_rel)

       IF ( verbose_level .GE. 2 ) WRITE(*,*) 'solve_rk_step: nl_iter',nl_iter

       CALL eval_f( qj , qj_old , a_tilde , a_dirk , a_diag , coeff_f ,         &
            right_term , scal_f )

       IF ( verbose_level .GE. 2 ) THEN

          WRITE(*,*) 'solve_rk_step: right_term',right_term

       END IF

       IF ( verbose_level .GE. 2 ) THEN

          WRITE(*,*) 'before_lnsrch: scal_f',scal_f

       END IF

       ! check the residual of the system

       IF ( MAXVAL( ABS( right_term(:) ) ) < TOLF ) THEN

          IF ( verbose_level .GE. 3 ) WRITE(*,*) '1: check',check
          RETURN

       END IF

       IF ( ( normalize_f ) .AND. ( scal_f < 1.D-6 ) ) THEN

          IF ( verbose_level .GE. 3 ) WRITE(*,*) 'check scal_f',check
          RETURN

       END IF

       ! ---- evaluate the descent direction ------------------------------------

       IF ( COUNT( implicit_flag ) .EQ. n_eqns ) THEN

          CALL eval_jacobian( qj_rel , qj_org,coeff_f , left_matrix )

          desc_dir_temp = - right_term

          CALL DGESV(n_eqns,1, left_matrix , n_eqns, pivot, desc_dir_temp ,     &
               n_eqns, ok)

          desc_dir = desc_dir_temp

       ELSE

          CALL eval_jacobian( qj_rel , qj_org,coeff_f , left_matrix )

          left_matrix_small11 = reshape(pack(left_matrix, mask11),              &
               [n_eqns-n_nh,n_eqns-n_nh]) 

          left_matrix_small12 = reshape(pack(left_matrix, mask12),              &
               [n_nh,n_eqns-n_nh]) 

          left_matrix_small22 = reshape(pack(left_matrix, mask22),              &
               [n_nh,n_nh]) 

          left_matrix_small21 = reshape(pack(left_matrix, mask21),              &
               [n_eqns-n_nh,n_nh]) 


          desc_dir_small1 = pack( right_term, .NOT.implicit_flag )
          desc_dir_small2 = pack( right_term , implicit_flag )

          DO i=1,n_vars-n_nh

             desc_dir_small1(i) = desc_dir_small1(i) / left_matrix_small11(i,i)

          END DO

          desc_dir_small2 = desc_dir_small2 -                                   &
               MATMUL( desc_dir_small1 , left_matrix_small21 )

          CALL DGESV(n_nh,1, left_matrix_small22 , n_nh , pivot_small2 ,        &
               desc_dir_small2 , n_nh, ok)

          desc_dir = unpack( - desc_dir_small2 , implicit_flag , 0.0D0 )        &
               + unpack( - desc_dir_small1 , .NOT.implicit_flag , 0.0D0 )

       END IF

       mod_desc_dir = DSQRT( desc_dir(2)**2 + desc_dir(3)**2 )

       !IF (  qj(2)**2 + qj(3)**2 .GT. 0.D0 ) THEN 
       !
       !   desc_dir(2) = mod_desc_dir * qj(2) / ( qj(2)**2 + qj(3)**2 ) 
       !   desc_dir(3) = mod_desc_dir * qj(3) / ( qj(2)**2 + qj(3)**2 ) 
       !
       !END IF

       IF ( verbose_level .GE. 3 ) WRITE(*,*) 'desc_dir',desc_dir

       qj_rel_NR_old = qj_rel
       scal_f_old = scal_f

       IF ( ( opt_search_NL ) .AND. ( nl_iter .GT. 1 ) ) THEN
          ! Search for the step lambda giving a suffic. decrease in the solution 

          stpmax = STPMX * MAX( SQRT( DOT_PRODUCT(qj_rel,qj_rel) ) ,            &
               DBLE( SIZE(qj_rel) ) )

          grad_f = MATMUL( right_term , left_matrix )

          desc_dir2 = desc_dir

          CALL lnsrch( qj_rel_NR_old , qj_org , qj_old , scal_f_old , grad_f ,  &
               desc_dir , coeff_f , qj_rel , scal_f , right_term , stpmax ,     &
               check )

       ELSE

          qj_rel = qj_rel_NR_old + desc_dir

          qj = qj_rel * qj_org

          CALL eval_f( qj , qj_old , a_tilde , a_dirk , a_diag , coeff_f ,      &
               right_term , scal_f )

       END IF

       IF ( verbose_level .GE. 2 ) WRITE(*,*) 'after_lnsrch: scal_f',scal_f

       qj = qj_rel * qj_org

       IF ( verbose_level .GE. 3 ) THEN

          WRITE(*,*) 'qj',qj
          CALL qc_to_qp( qj , Bj , qpj)
          WRITE(*,*) 'qp',qpj

       END IF

       
       IF ( MAXVAL( ABS( right_term(:) ) ) < TOLF ) THEN

          IF ( verbose_level .GE. 3 ) WRITE(*,*) '1: check',check
          check= .FALSE.
          RETURN

       END IF

       IF (check) THEN

          check = ( MAXVAL( ABS(grad_f(:)) * MAX( ABS( qj_rel(:) ),1.D0 ) /     &
               MAX( scal_f , 0.5D0 * SIZE(qj_rel) ) )  < TOLMIN )

          IF ( verbose_level .GE. 3 ) WRITE(*,*) '2: check',check
          !          RETURN

       END IF

       IF ( MAXVAL( ABS( qj_rel(:) - qj_rel_NR_old(:) ) / MAX( ABS( qj_rel(:)) ,&
            1.D0 ) ) < TOLX ) THEN

          IF ( verbose_level .GE. 3 ) WRITE(*,*) 'check',check
          RETURN

       END IF

    END DO newton_raphson_loop

  END SUBROUTINE solve_rk_step

  !******************************************************************************
  !> \brief Search the descent stepsize
  !
  !> This subroutine search for the lenght of the descent step in order to have
  !> a decrease in the nonlinear function.
  !> \param[in]     Bj               topography at the cell center
  !> \param[in]     Bprimej_x        topography x-slope at the cell center
  !> \param[in]     Bprimej_y        topography y-slope at the cell center
  !> \param[in]     qj_rel_NR_old  
  !> \param[in]     qj_org
  !> \param[in]     qj_old
  !> \param[in]     scal_f_old
  !> \param[in]     grad_f
  !> \param[in,out] desc_dir
  !> \param[in]     coeff_f
  !> \param[out]    qj_rel
  !> \param[out]    scal_f
  !> \param[out]    right_term
  !> \param[in]     stpmax
  !> \param[out]    check
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE lnsrch( qj_rel_NR_old , qj_org , qj_old , scal_f_old , grad_f ,    &
       desc_dir , coeff_f , qj_rel , scal_f , right_term , stpmax , check )

    IMPLICIT NONE

    !> Initial point
    REAL*8, DIMENSION(:), INTENT(IN) :: qj_rel_NR_old

    !> Initial point
    REAL*8, DIMENSION(:), INTENT(IN) :: qj_org

    !> Initial point
    REAL*8, DIMENSION(:), INTENT(IN) :: qj_old

    !> Gradient at xold
    REAL*8, DIMENSION(:), INTENT(IN) :: grad_f

    !> Value of the function at xold
    REAL*8, INTENT(IN) :: scal_f_old

    !> Descent direction (usually Newton direction)
    REAL*8, DIMENSION(:), INTENT(INOUT) :: desc_dir

    REAL*8, INTENT(IN) :: stpmax

    !> Coefficients to rescale the nonlinear function
    REAL*8, DIMENSION(:), INTENT(IN) :: coeff_f

    !> Updated solution
    REAL*8, DIMENSION(:), INTENT(OUT) :: qj_rel

    !> Value of the scalar function at x
    REAL*8, INTENT(OUT) :: scal_f

    !> Value of the scalar function at x
    REAL*8, INTENT(OUT) :: right_term(n_eqns)

    !> Output quantity check is false on a normal exit 
    LOGICAL, INTENT(OUT) :: check

    REAL*8, PARAMETER :: TOLX=epsilon(qj_rel)

    INTEGER, DIMENSION(1) :: ndum
    REAL*8 :: ALF , a,alam,alam2,alamin,b,disc
    REAL*8 :: scal_f2
    REAL*8 :: desc_dir_abs
    REAL*8 :: rhs1 , rhs2 , slope, tmplam

    REAL*8 :: scal_f_min , alam_min

    REAL*8 :: qj(n_vars)

    ALF = 1.0d-4

    IF ( size(grad_f) == size(desc_dir) .AND. size(grad_f) == size(qj_rel)      &
         .AND. size(qj_rel) == size(qj_rel_NR_old) ) THEN

       ndum = size(grad_f)

    ELSE

       WRITE(*,*) 'nrerror: an assert_eq failed with this tag:', 'lnsrch'
       STOP 'program terminated by assert_eq4'

    END IF

    check = .FALSE.

    desc_dir_abs = SQRT( DOT_PRODUCT(desc_dir,desc_dir) )

    IF ( desc_dir_abs > stpmax ) desc_dir(:) = desc_dir(:) * stpmax/desc_dir_abs  

    slope = DOT_PRODUCT(grad_f,desc_dir)

    alamin = TOLX / MAXVAL( ABS( desc_dir(:))/MAX( ABS(qj_rel_NR_old(:)),1.D0 ) )

    IF ( alamin .EQ. 0.d0) THEN

       qj_rel(:) = qj_rel_NR_old(:)

       RETURN

    END IF

    alam = 1.0D0

    scal_f_min = scal_f_old

    optimal_step_search: DO

       IF ( verbose_level .GE. 4 ) THEN

          WRITE(*,*) 'alam',alam

       END IF

       qj_rel = qj_rel_NR_old + alam * desc_dir

       qj = qj_rel * qj_org

       CALL eval_f( qj , qj_old , a_tilde , a_dirk , a_diag , coeff_f ,         &
            right_term , scal_f )

       IF ( verbose_level .GE. 4 ) THEN

          WRITE(*,*) 'lnsrch: effe_old,effe',scal_f_old,scal_f
          READ(*,*)

       END IF

       IF ( scal_f .LT. scal_f_min ) THEN

          scal_f_min = scal_f
          alam_min = alam

       END IF

       IF ( scal_f .LE. 0.9 * scal_f_old ) THEN   
          ! sufficient function decrease

          IF ( verbose_level .GE. 4 ) THEN

             WRITE(*,*) 'sufficient function decrease'

          END IF

          EXIT optimal_step_search   

       ELSE IF ( alam < alamin ) THEN   
          ! convergence on Delta_x

          IF ( verbose_level .GE. 4 ) THEN

             WRITE(*,*) ' convergence on Delta_x',alam,alamin

          END IF

          qj_rel(:) = qj_rel_NR_old(:)
          scal_f = scal_f_old
          check = .TRUE.

          EXIT optimal_step_search

          !       ELSE IF ( scal_f .LE. scal_f_old + ALF * alam * slope ) THEN   
       ELSE  

          IF ( alam .EQ. 1.D0 ) THEN

             tmplam = - slope / ( 2.0D0 * ( scal_f - scal_f_old - slope ) )

          ELSE

             rhs1 = scal_f - scal_f_old - alam*slope
             rhs2 = scal_f2 - scal_f_old - alam2*slope

             a = ( rhs1/alam**2.D0 - rhs2/alam2**2.D0 ) / ( alam - alam2 )
             b = ( -alam2*rhs1/alam**2 + alam*rhs2/alam2**2 ) / ( alam - alam2 )

             IF ( a .EQ. 0.D0 ) THEN

                tmplam = - slope / ( 2.0D0 * b )

             ELSE

                disc = b*b - 3.0D0*a*slope

                IF ( disc .LT. 0.D0 ) THEN

                   tmplam = 0.5D0 * alam

                ELSE IF ( b .LE. 0.D0 ) THEN

                   tmplam = ( - b + SQRT(disc) ) / ( 3.D0 * a )

                ELSE

                   tmplam = - slope / ( b + SQRT(disc) )

                ENDIF

             END IF

             IF ( tmplam .GT. 0.5D0*alam ) tmplam = 0.5D0 * alam

          END IF

       END IF

       alam2 = alam
       scal_f2 = scal_f
       alam = MAX( tmplam , 0.5D0*alam )

    END DO optimal_step_search

  END SUBROUTINE lnsrch

  !******************************************************************************
  !> \brief Evaluate the nonlinear system
  !
  !> This subroutine evaluate the value of the nonlinear system in the state 
  !> defined by the variables qj.
  !> \param[in]    Bj          topography at the cell center
  !> \param[in]    Bprimej     topography slope at the cell center
  !> \param[in]    qj          conservative variables 
  !> \param[in]    qj_old      conservative variables at the old time step
  !> \param[in]    a_tilde     explicit coefficients for the hyperbolic terms 
  !> \param[in]    a_dirk      explicit coefficients for the non-hyperbolic terms 
  !> \param[in]    a_diag      implicit coefficient for the non-hyperbolic term
  !> \param[in]    coeff_f     coefficient to rescale the nonlinear functions
  !> \param[out]   f_nl        values of the nonlinear functions
  !> \param[out]   scal_f      value of the scalar function f=0.5*<F,F>
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE eval_f( qj , qj_old , a_tilde , a_dirk , a_diag , coeff_f , f_nl , &
       scal_f )

    USE constitutive_2d, ONLY : eval_nonhyperbolic_terms

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: qj(n_vars)
    REAL*8, INTENT(IN) :: qj_old(n_vars)
    REAL*8, INTENT(IN) :: a_tilde(n_RK)
    REAL*8, INTENT(IN) :: a_dirk(n_RK)
    REAL*8, INTENT(IN) :: a_diag
    REAL*8, INTENT(IN) :: coeff_f(n_eqns)
    REAL*8, INTENT(OUT) :: f_nl(n_eqns)
    REAL*8, INTENT(OUT) :: scal_f

    REAL*8 :: nh_term_impl(n_eqns)
    REAL*8 :: Rj(n_eqns)

    CALL eval_nonhyperbolic_terms( r_qj = qj , r_nh_term_impl=nh_term_impl ) 

    Rj = ( MATMUL( divFluxj(1:n_eqns,1:i_RK-1)+Expl_terms_j(1:n_eqns,1:i_RK-1), &
         a_tilde(1:i_RK-1) ) - MATMUL( NHj(1:n_eqns,1:i_RK-1)                   &
         + SI_NHj(1:n_eqns,1:i_RK-1) , a_dirk(1:i_RK-1) ) )                     &
         - a_diag * ( nh_term_impl + SI_NHj(1:n_eqns,i_RK) )

    f_nl = qj - qj_old + dt * Rj

    f_nl = coeff_f * f_nl

    scal_f = 0.5D0 * DOT_PRODUCT( f_nl , f_nl )

  END SUBROUTINE eval_f

  !******************************************************************************
  !> \brief Evaluate the jacobian 
  !
  !> This subroutine evaluate the jacobian of the non-linear system
  !> with respect to the conservative variables.
  !
  !> \param[in]    Bj            topography at the cell center
  !> \param[in]    Bprimej_x     topography x-slope at the cell center
  !> \param[in]    Bprimej_y     topography y-slope at the cell center
  !> \param[in]    grav3_surf
  !> \param[in]    qj_rel        relative variation (qj=qj_rel*qj_org)
  !> \param[in]    qj_org        conservative variables at the old time step
  !> \param[in]    coeff_f       coefficient to rescale the nonlinear functions
  !> \param[out]   left_matrix   matrix from the linearization of the system
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE eval_jacobian( qj_rel , qj_org , coeff_f, left_matrix)

    USE constitutive_2d, ONLY : eval_nonhyperbolic_terms

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: qj_rel(n_vars)
    REAL*8, INTENT(IN) :: qj_org(n_vars)
    REAL*8, INTENT(IN) :: coeff_f(n_eqns)

    REAL*8, INTENT(OUT) :: left_matrix(n_eqns,n_vars)

    REAL*8 :: Jacob_relax(n_eqns,n_vars)
    COMPLEX*16 :: nh_terms_cmplx_impl(n_eqns)
    COMPLEX*16 :: qj_cmplx(n_vars) , qj_rel_cmplx(n_vars)

    REAL*8 :: h

    INTEGER :: i

    h = n_vars * epsilon(1.d0)

    ! initialize the matrix of the linearized system and the Jacobian

    left_matrix(1:n_eqns,1:n_vars) = 0.D0
    Jacob_relax(1:n_eqns,1:n_vars) = 0.D0

    ! evaluate the jacobian of the non-hyperbolic terms

    DO i=1,n_vars

       left_matrix(i,i) = coeff_f(i) * qj_org(i)

       IF ( implicit_flag(i) ) THEN 

          qj_rel_cmplx(1:n_vars) = qj_rel(1:n_vars)
          qj_rel_cmplx(i) = DCMPLX(qj_rel(i), h)

          qj_cmplx = qj_rel_cmplx * qj_org

          CALL eval_nonhyperbolic_terms( c_qj = qj_cmplx ,                      &
               c_nh_term_impl = nh_terms_cmplx_impl ) 

          Jacob_relax(1:n_eqns,i) = coeff_f(i) *                                &
               AIMAG(nh_terms_cmplx_impl) / h

          left_matrix(1:n_eqns,i) = left_matrix(1:n_eqns,i) - dt * a_diag       &
               * Jacob_relax(1:n_eqns,i)

       END IF

    END DO

  END SUBROUTINE eval_jacobian

  !******************************************************************************
  !> \brief Evaluate the eroion/deposition terms
  !
  !> This subroutine update the solution and the topography computing the 
  !> erosion and deposition terms and the solution only because of entrainment.
  !
  !> \param[in]    dt      time step
  !
  !> \date 2019/11/08
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE update_erosion_deposition_cell(dt)

    USE constitutive_2d, ONLY : erosion_coeff , settling_vel , rho_s , red_grav

    USE geometry_2d, ONLY : deposit

    USE constitutive_2d, ONLY : eval_erosion_dep_term
    USE constitutive_2d, ONLY : eval_topo_term

    USE constitutive_2d, ONLY : qc_to_qp
    USE parameters_2d, ONLY : topo_change_flag

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: dt

    REAL*8 :: erosion_term(n_solid)
    REAL*8 :: deposition_term(n_solid)
    REAL*8 :: eqns_term(n_eqns)
    REAL*8 :: deposit_term(n_solid)

    REAL*8 :: u_old , u_new
    REAL*8 :: v_old , v_new

    REAL*8 :: max_vel_old , max_vel_new
    
    INTEGER :: j ,k

    max_vel_old = 0.D0
    max_vel_new = 0.D0
    
    IF ( ( erosion_coeff .EQ. 0.D0 ) .AND. ( settling_vel .EQ. 0.D0 ) ) RETURN

    DO k = 1,comp_cells_y

       DO j = 1,comp_cells_x

          CALL qc_to_qp(q(1:n_vars,j,k) , B_cent(j,k) , qp(1:n_vars,j,k) )

          CALL eval_erosion_dep_term( qp(:,j,k) , B_cent(j,k) , dt ,            &
               erosion_term(1:n_solid) , deposition_term(1:n_solid) )

          ! Limit the deposition during a single time step
          deposition_term(1:n_solid) = MAX(0.D0,MIN( deposition_term(1:n_solid),&
               q(5:4+n_solid,j,k) / ( rho_s(1:n_solid) * dt ) ))

          ! Compute the source terms for the equations
          CALL eval_topo_term( q(1:n_vars,j,k) , deposition_term ,              &
               erosion_term , eqns_term , deposit_term )

          IF ( verbose_level .GE. 2 ) THEN

             WRITE(*,*) 'before update erosion/deposition: j,k,q(:,j,k),B(j,k)',&
                  j,k,q(:,j,k),B_cent(j,k)

          END IF

          IF ( deposition_term(1) .GT. 0.D0 ) THEN
             
             CALL qc_to_qp(q(1:n_vars,j,k) , B_cent(j,k) , qp(1:n_vars,j,k))

             u_old = qp(2,j,k)/qp(1,j,k)
             v_old = qp(3,j,k)/qp(1,j,k)

             max_vel_old = MAX(max_vel_old,DSQRT(u_old**2+v_old**2))
             
          END IF
                

          ! Update the solution with erosion/deposition terms
          q(1:n_eqns,j,k) = q(1:n_eqns,j,k) + dt * eqns_term(1:n_eqns)

          IF ( deposition_term(1) .GT. 0.D0 ) THEN
             
             CALL qc_to_qp(q(1:n_vars,j,k) , B_cent(j,k) , qp(1:n_vars,j,k))

             u_new = qp(2,j,k)/qp(1,j,k)
             v_new = qp(3,j,k)/qp(1,j,k)
             
             max_vel_new = MAX(max_vel_new,DSQRT(u_new**2+v_new**2))

!!$             IF ( DABS( u_new - u_old ) .GT. 1.D-3) THEN
!!$
!!$                WRITE(*,*) 'u',qp(1,j,k),u_old,u_new
!!$                READ(*,*)
!!$
!!$             END IF
!!$
!!$             IF ( DABS( v_new - v_old ) .GT. 1.D-3) THEN
!!$
!!$                WRITE(*,*) 'v',qp(1,j,k),u_old,u_new
!!$                READ(*,*)
!!$
!!$             END IF

          END IF

          deposit(j,k,1:n_solid) =  deposit(j,k,1:n_solid) + dt * deposit_term(1:n_solid)

          ! Update the topography with erosion/deposition terms
          IF ( topo_change_flag ) THEN

             B_cent(j,k) = B_cent(j,k) + dt * SUM( deposit_term(1:n_solid) )

          END IF


          ! Check for negative thickness
          IF ( q(1,j,k) .LT. 0.D0 ) THEN

             IF ( q(1,j,k) .GT. -1.D-10 ) THEN

                q(1:n_vars,j,k) = 0.D0

             ELSE

                WRITE(*,*) 'j,k',j,k
                WRITE(*,*) 'dt',dt
                WRITE(*,*) 'before erosion'
                WRITE(*,*) 'qp',qp(1:n_eqns,j,k)
                WRITE(*,*) 'deposition_term',deposition_term
                WRITE(*,*) 'erosion_term',erosion_term
                CALL qc_to_qp(q(1:n_vars,j,k) , B_cent(j,k) , qp(1:n_vars,j,k))
                WRITE(*,*) 'qp',qp(1:n_eqns,j,k)
                READ(*,*)

             END IF

          END IF

          CALL qc_to_qp(q(1:n_vars,j,k) , B_cent(j,k) , qp(1:n_vars,j,k) )

          IF ( REAL(red_grav) .LT. 0.D0 ) THEN

             q(1:n_vars,j,k) = 0.D0

          END IF
          
       END DO

    END DO

!!$    WRITE(*,*) 'max_vel',max_vel_old,max_vel_new
    
    RETURN

  END SUBROUTINE update_erosion_deposition_cell

  !******************************************************************************
  !> \brief Evaluate the explicit terms 
  !
  !> This subroutine evaluate the explicit terms (non-fluxes) of the non-linear 
  !> system with respect to the conservative variables.
  !
  !> \param[in]    q_expl          conservative variables 
  !> \param[out]   expl_terms      explicit terms
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE eval_explicit_terms( qc_expl , qp_expl , expl_terms )

    USE constitutive_2d, ONLY : eval_expl_terms

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: qc_expl(n_vars,comp_cells_x,comp_cells_y)
    REAL*8, INTENT(IN) :: qp_expl(n_vars,comp_cells_x,comp_cells_y)
    REAL*8, INTENT(OUT) :: expl_terms(n_eqns,comp_cells_x,comp_cells_y)

    REAL*8 :: qcj(n_vars)     !< local conservative variables 
    REAL*8 :: qpj(n_vars)     !< local physical variables
    REAL*8 :: expl_forces_term(n_eqns)      !< conservative variables 

    INTEGER :: j,k

    expl_terms = 0.D0

    DO k = 1,comp_cells_y
    
       DO j = 1,comp_cells_x

          IF ( solve_mask(j,k) ) THEN

             qcj = qc_expl(1:n_vars,j,k)
             qpj = qp_expl(1:n_vars,j,k)
             
             CALL eval_expl_terms( B_prime_x(j,k), B_prime_y(j,k), source_xy(j,k) ,&
                  qpj , qcj , expl_forces_term )
             
             expl_terms(1:n_eqns,j,k) =  expl_forces_term
             
          END IF

       ENDDO

    END DO

  END SUBROUTINE eval_explicit_terms

  !******************************************************************************
  !> \brief Semidiscrete finite volume central scheme
  !
  !> This subroutine compute the divergence part of the system of the eqns,
  !> with a modified version of the finite volume scheme from Kurganov et al.  
  !> 2001, where the reconstruction at the cells interfaces is applied to a
  !> set of physical variables derived from the conservative vriables.
  !
  !> \param[in]     q_expl        conservative variables
  !> \param[out]    divFlux           divergence term
  !
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE eval_hyperbolic_terms( q_expl , qp_expl , divFlux )

    ! External variables
    USE geometry_2d, ONLY : dx,dy
    USE parameters_2d, ONLY : solver_scheme

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: q_expl(n_vars,comp_cells_x,comp_cells_y)
    REAL*8, INTENT(IN) :: qp_expl(n_vars,comp_cells_x,comp_cells_y)
    REAL*8, INTENT(OUT) :: divFlux(n_eqns,comp_cells_x,comp_cells_y)

    REAL*8 :: q_old(n_vars,comp_cells_x,comp_cells_y)

    REAL*8 :: h_new , h_old

    REAL*8 :: tcpu0,tcpu1,tcpu2,tcpu3,tcpu4

    INTEGER :: i, j, k      !< loop counters

    q_old = q

    q = q_expl

    CALL cpu_time(tcpu0)

    ! Linear reconstruction of the physical variables at the interfaces
    CALL reconstruction(q_expl,qp_expl)

    CALL cpu_time(tcpu1)
    !WRITE(*,*) 'eval_hyperbolic_terms: Time taken by the code was',tcpu1-tcpu0,'seconds'


    ! Evaluation of the maximum local speeds at the interfaces
    CALL eval_speeds

    CALL cpu_time(tcpu2)
    !WRITE(*,*) 'eval_hyperbolic_terms: Time taken by the code was',tcpu2-tcpu1,'seconds'

    ! Evaluation of the numerical fluxes
    SELECT CASE ( solver_scheme )

    CASE ("LxF")

       CALL eval_flux_LxF

    CASE ("GFORCE")

       CALL eval_flux_GFORCE

    CASE ("KT")

       CALL eval_flux_KT

    CASE ("UP")

       CALL eval_flux_UP
       
    END SELECT

    CALL cpu_time(tcpu3)
    !WRITE(*,*) 'eval_hyperbolic_terms: Time taken by the code was',tcpu3-tcpu2,'seconds'

    ! Advance in time the solution
    DO k = 1,comp_cells_y
    
       DO j = 1,comp_cells_x
  
          DO i=1,n_eqns
             
             divFlux(i,j,k) = 0.D0
             
             IF ( comp_cells_x .GT. 1 ) THEN
                
                divFlux(i,j,k) = divFlux(i,j,k) +                               &
                     ( H_interface_x(i,j+1,k) - H_interface_x(i,j,k) ) / dx
                
             END IF

             IF ( comp_cells_y .GT. 1 ) THEN
                
                divFlux(i,j,k) = divFlux(i,j,k) +                               &
                     ( H_interface_y(i,j,k+1) - H_interface_y(i,j,k) ) / dy
                
             END IF
                
          END DO
         
          h_old = q_expl(1,j,k)
          h_new = h_old - dt * divFlux(1,j,k)

       ENDDO

    END DO
    
    CALL cpu_time(tcpu4)
    !WRITE(*,*) 'eval_hyperbolic_terms: Time taken by the code was',tcpu4-tcpu4,'seconds'

    q = q_old

    RETURN

  END SUBROUTINE eval_hyperbolic_terms

  !******************************************************************************
  !> \brief Upwind numerical fluxes
  !
  !> This subroutine evaluates the numerical fluxes H at the 
  !> cells interfaces with an upwind discretization.
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 2019/11/16
  !******************************************************************************

  SUBROUTINE eval_flux_UP

    ! External procedures
    USE constitutive_2d, ONLY : eval_fluxes

    IMPLICIT NONE

    REAL*8 :: fluxL(n_eqns)           !< Numerical fluxes from the eqns 
    REAL*8 :: fluxR(n_eqns)           !< Numerical fluxes from the eqns
    REAL*8 :: fluxB(n_eqns)           !< Numerical fluxes from the eqns 
    REAL*8 :: fluxT(n_eqns)           !< Numerical fluxes from the eqns

    INTEGER :: j,k                  !< Loop counters

    H_interface_x = 0.D0
    H_interface_y = 0.D0
        
    IF ( comp_cells_x .GT. 1 ) THEN

       DO k = 1,comp_cells_y
       
          DO j = 1,comp_interfaces_x
             
             CALL eval_fluxes( q_interfaceL(1:n_vars,j,k) ,                     &
                  qp_interfaceL(1:n_vars,j,k) , B_interfaceL(j,k) , 1 , fluxL)

             CALL eval_fluxes( q_interfaceR(1:n_vars,j,k) ,                     &
                  qp_interfaceR(1:n_vars,j,k) , B_interfaceR(j,k) , 1 , fluxR)

             IF ( ( q_interfaceL(2,j,k) .GT. 0.D0 ) .AND.                       &
                  ( q_interfaceR(2,j,k) .GE. 0.D0 ) ) THEN

                H_interface_x(:,j,k) = fluxL

             ELSEIF ( ( q_interfaceL(2,j,k) .LE. 0.D0 ) .AND.                   &
                  ( q_interfaceR(2,j,k) .LT. 0.D0 ) ) THEN

                H_interface_x(:,j,k) = fluxR

             ELSE

                H_interface_x(:,j,k) = 0.5D0 * ( fluxL + fluxR )

             END IF
                
             IF ( (  q_interfaceL(2,j,k) .EQ. 0.D0 ) .AND.                      &
                  (  q_interfaceR(2,j,k) .EQ. 0.D0 ) ) THEN

                H_interface_x(1,j,k) = 0.D0
                H_interface_x(4:n_vars,j,k) = 0.D0

             END IF

          END DO

       END DO

    END IF

    !READ(*,*)
    
    IF ( comp_cells_y .GT. 1 ) THEN

       DO k = 1,comp_interfaces_y
          
          DO j = 1,comp_cells_x

             CALL eval_fluxes( q_interfaceB(1:n_vars,j,k) ,                     &
                  qp_interfaceB(1:n_vars,j,k) , B_interfaceB(j,k) , 2 , fluxB)
             
             CALL eval_fluxes( q_interfaceT(1:n_vars,j,k) ,                     &
                  qp_interfaceT(1:n_vars,j,k) , B_interfaceT(j,k) , 2 , fluxT)
       
             IF ( ( q_interfaceB(3,j,k) .GT. 0.D0 ) .AND.                       &
                  ( q_interfaceT(3,j,k) .GE. 0.D0 ) ) THEN

                H_interface_y(:,j,k) = fluxB

             ELSEIF ( ( q_interfaceB(3,j,k) .LE. 0.D0 ) .AND.                   &
                  ( q_interfaceT(3,j,k) .LT. 0.D0 ) ) THEN

                H_interface_y(:,j,k) = fluxT

             ELSE

                H_interface_y(:,j,k) = 0.5D0 * ( fluxB + fluxT )

             END IF
      
             ! In the equation for mass and for trasnport (T,alphas) if the 
             ! velocities at the interfaces are null, then the flux is null
             IF ( (  q_interfaceB(3,j,k) .EQ. 0.D0 ) .AND.                      &
                  (  q_interfaceT(3,j,k) .EQ. 0.D0 ) ) THEN

                H_interface_y(1,j,k) = 0.D0
                H_interface_y(4:n_vars,j,k) = 0.D0

             END IF
             
          ENDDO

       END DO

    END IF

  END SUBROUTINE eval_flux_UP

  
  !******************************************************************************
  !> \brief Semidiscrete numerical fluxes
  !
  !> This subroutine evaluates the numerical fluxes H at the 
  !> cells interfaces according to Kurganov et al. 2001. 
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 16/08/2011
  !******************************************************************************

  SUBROUTINE eval_flux_KT

    ! External procedures
    USE constitutive_2d, ONLY : eval_fluxes

    IMPLICIT NONE

    REAL*8 :: fluxL(n_eqns)           !< Numerical fluxes from the eqns 
    REAL*8 :: fluxR(n_eqns)           !< Numerical fluxes from the eqns
    REAL*8 :: fluxB(n_eqns)           !< Numerical fluxes from the eqns 
    REAL*8 :: fluxT(n_eqns)           !< Numerical fluxes from the eqns

    REAL*8 :: flux_avg_x(n_eqns)   
    REAL*8 :: flux_avg_y(n_eqns)   

    INTEGER :: i,j,k                  !< Loop counters

    H_interface_x = 0.D0
    H_interface_y = 0.D0
        
    IF ( comp_cells_x .GT. 1 ) THEN

       DO k = 1,comp_cells_y
       
          DO j = 1,comp_interfaces_x

             CALL eval_fluxes( q_interfaceL(1:n_vars,j,k) ,                     &
                  qp_interfaceL(1:n_vars,j,k) , B_interfaceL(j,k) , 1 , fluxL)

             CALL eval_fluxes( q_interfaceR(1:n_vars,j,k) ,                     &
                  qp_interfaceR(1:n_vars,j,k) , B_interfaceR(j,k) , 1 , fluxR)

             CALL average_KT( a_interface_xNeg(:,j,k), a_interface_xPos(:,j,k) ,&
                  fluxL , fluxR , flux_avg_x )

             eqns_loop:DO i=1,n_eqns

                IF ( a_interface_xNeg(i,j,k) .EQ. a_interface_xPos(i,j,k) ) THEN

                   H_interface_x(i,j,k) = 0.D0

                ELSE

                   H_interface_x(i,j,k) = flux_avg_x(i)                         &
                        + ( a_interface_xPos(i,j,k) * a_interface_xNeg(i,j,k) ) &
                        / ( a_interface_xPos(i,j,k) - a_interface_xNeg(i,j,k) ) &
                        * ( q_interfaceR(i,j,k) - q_interfaceL(i,j,k) )             

                END IF

             ENDDO eqns_loop
             
             ! In the equation for mass and for trasnport (T,alphas) if the 
             ! velocities at the interfaces are null, then the flux is null
             IF ( (  qp_interfaceL(2,j,k) .EQ. 0.D0 ) .AND.                     &
                  (  qp_interfaceR(2,j,k) .EQ. 0.D0 ) ) THEN

                H_interface_x(1,j,k) = 0.D0
                H_interface_x(4:n_vars,j,k) = 0.D0

             END IF

          END DO

       END DO

    END IF

    !READ(*,*)
    
    IF ( comp_cells_y .GT. 1 ) THEN

       DO k = 1,comp_interfaces_y
          
          DO j = 1,comp_cells_x

             CALL eval_fluxes( q_interfaceB(1:n_vars,j,k) ,                     &
                  qp_interfaceB(1:n_vars,j,k) , B_interfaceB(j,k) , 2 , fluxB)    
         
             CALL eval_fluxes( q_interfaceT(1:n_vars,j,k) ,                     &
                  qp_interfaceT(1:n_vars,j,k) , B_interfaceT(j,k) , 2 , fluxT)

             CALL average_KT( a_interface_yNeg(:,j,k) ,                         &
                  a_interface_yPos(:,j,k) , fluxB , fluxT , flux_avg_y )
           
             DO i=1,n_eqns

                IF ( a_interface_yNeg(i,j,k) .EQ. a_interface_yPos(i,j,k) ) THEN

                   H_interface_y(i,j,k) = 0.D0

                ELSE

                   H_interface_y(i,j,k) = flux_avg_y(i)                         &
                        + ( a_interface_yPos(i,j,k) * a_interface_yNeg(i,j,k) ) &
                        / ( a_interface_yPos(i,j,k) - a_interface_yNeg(i,j,k) ) &
                        * ( q_interfaceT(i,j,k) - q_interfaceB(i,j,k) )             

                END IF

             END DO

             ! In the equation for mass and for trasnport (T,alphas) if the 
             ! velocities at the interfaces are null, then the flux is null
             IF ( (  q_interfaceB(3,j,k) .EQ. 0.D0 ) .AND.                      &
                  (  q_interfaceT(3,j,k) .EQ. 0.D0 ) ) THEN

                H_interface_y(1,j,k) = 0.D0
                H_interface_y(4:n_vars,j,k) = 0.D0

             END IF
             
          ENDDO

       END DO

    END IF

  END SUBROUTINE eval_flux_KT

  !******************************************************************************
  !> \brief averaged KT flux
  !
  !> This subroutine compute n averaged flux from the fluxes at the two sides of
  !> a cell interface and the max an min speed at the two sides.
  !> \param[in]     aL            speed at one side of the interface
  !> \param[in]     aR            speed at the other side of the interface
  !> \param[in]     wL            fluxes at one side of the interface
  !> \param[in]     wR            fluxes at the other side of the interface
  !> \param[out]    w_avg         array of averaged fluxes
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************


  SUBROUTINE average_KT( a1 , a2 , w1 , w2 , w_avg )

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: a1(:) , a2(:)
    REAL*8, INTENT(IN) :: w1(:) , w2(:)
    REAL*8, INTENT(OUT) :: w_avg(:)

    INTEGER :: n
    INTEGER :: i 

    n = SIZE( a1 )

    DO i=1,n

       IF ( a1(i) .EQ. a2(i) ) THEN

          w_avg(i) = 0.5D0 * ( w1(i) + w2(i) )
          w_avg(i) = 0.D0

       ELSE

          w_avg(i) = ( a2(i) * w1(i) - a1(i) * w2(i) ) / ( a2(i) - a1(i) )  

       END IF

    END DO

  END SUBROUTINE average_KT

  !******************************************************************************
  !> \brief Numerical fluxes GFORCE
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE eval_flux_GFORCE

    ! to be implemented
    WRITE(*,*) 'method not yet implemented in 2-d case'

  END SUBROUTINE eval_flux_GFORCE

  !******************************************************************************
  !> \brief Numerical fluxes Lax-Friedrichs
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE eval_flux_LxF

    ! to be implemented
    WRITE(*,*) 'method not yet implemented in 2-d case'

  END SUBROUTINE eval_flux_LxF


  !******************************************************************************
  !> \brief Linear reconstruction
  !
  !> In this subroutine a linear reconstruction with slope limiters is
  !> applied to a set of variables describing the state of the system.
  !> In this way the values at the two sides of each cell interface are computed.
  !> This subroutine is also used for the boundary condition, when the
  !> reconstruction at the boundary interfaces are computed.
  !> \param[in]     q_expl         center values of the conservative variables
  !> \param[in]     qp_expl        center values of the physical variables
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 2019/11/11
  !******************************************************************************

  SUBROUTINE reconstruction(q_expl,qp_expl)

    ! External procedures
    USE constitutive_2d, ONLY : qc_to_qp , qp_to_qc ,qp_to_qp2
    USE constitutive_2d, ONLY : eval_source_bdry
    USE parameters_2d, ONLY : limiter

    ! External variables
    USE geometry_2d, ONLY : x_comp , x_stag , y_comp , y_stag , dx , dx2 , dy , &
         dy2

    USE geometry_2d, ONLY : B_cent
    USE geometry_2d, ONLY : B_interfaceL , B_interfaceR
    USE geometry_2d, ONLY : B_interfaceT , B_interfaceB
    
    USE geometry_2d, ONLY : sourceW , sourceE , sourceN , sourceS
    USE geometry_2d, ONLY : sourceW_vect_x , sourceW_vect_y
    USE geometry_2d, ONLY : sourceE_vect_x , sourceE_vect_y
    USE geometry_2d, ONLY : sourceN_vect_x , sourceN_vect_y
    USE geometry_2d, ONLY : sourceS_vect_x , sourceS_vect_y

    USE parameters_2d, ONLY : reconstr_coeff

    USE geometry_2d, ONLY : minmod

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: q_expl(n_vars,comp_cells_x,comp_cells_y)
    REAL*8, INTENT(IN) :: qp_expl(n_vars,comp_cells_x,comp_cells_y)

    REAL*8 :: qrec(n_vars,comp_cells_x,comp_cells_y) 
    REAL*8 :: qrecW(n_vars)     !< recons var at the west edge of the cells
    REAL*8 :: qrecE(n_vars)     !< recons var at the east edge of the cells
    REAL*8 :: qrecS(n_vars)     !< recons var at the south edge of the cells
    REAL*8 :: qrecN(n_vars)     !< recons var at the north edge of the cells

    REAL*8 :: qp2recC(3)
    REAL*8 :: qp2recL(3) , qp2recW(3) , qp2recE(3) , qp2recR(3)
    REAL*8 :: qp2recB(3) , qp2recS(3) , qp2recN(3) , qp2recT(3)

    REAL*8 :: source_bdry(n_vars)

    REAL*8 :: qrec_stencil(3)   !< recons variables stencil for the limiter
    REAL*8 :: x_stencil(3)    !< grid stencil for the limiter
    REAL*8 :: y_stencil(3)    !< grid stencil for the limiter
    REAL*8 :: qrec_prime_x(n_vars)      !< recons variables slope
    REAL*8 :: qrec_prime_y(n_vars)      !< recons variables slope

    REAL*8 :: qp2rec_prime_x(3)      !< recons variables slope
    REAL*8 :: qp2rec_prime_y(3)      !< recons variables slope

    INTEGER :: j,k            !< loop counters (cells)
    INTEGER :: i              !< loop counter (variables)

    REAL*8 :: gamma_W , gamma_E
    REAL*8 :: gamma_N , gamma_S
    REAL*8 :: gamma
    

    ! Compute the variable to reconstruct (phys or cons)
    DO k = 1,comp_cells_y

       DO j = 1,comp_cells_x
          
          qrec(1:n_vars,j,k) = qp_expl(1:n_vars,j,k)
          
          IF ( SUM(qrec(5:4+n_solid,j,k)) .GT. 1.D0 ) THEN
             
             WRITE(*,*) 'reconstruction: j,k',j,k
             WRITE(*,*) 'qrec(5:n_solid,j,k)',qrec(5:4+n_solid,j,k)
             WRITE(*,*) 'q(1:n_vars,j,k)',q_expl(1:n_vars,j,k)
             WRITE(*,*) 'B_cent(j,k)', B_cent(j,k)
             READ(*,*)
             
          END IF
          
       END DO
       
    END DO
    
    ! Linear reconstruction

    y_loop:DO k = 1,comp_cells_y

       x_loop:DO j = 1,comp_cells_x

          IF ( solve_mask(j,k) ) THEN
          ! IF ( ( .NOT. radial_source_flag ) .OR. ( source_cell(j,k) .NE. 1 ) ) THEN

             qrecW(1:n_vars) =  qp_expl(1:n_vars,j,k)
             qrecE(1:n_vars) =  qp_expl(1:n_vars,j,k)
             qrecS(1:n_vars) =  qp_expl(1:n_vars,j,k)
             qrecN(1:n_vars) =  qp_expl(1:n_vars,j,k)

          vars_loop:DO i=1,n_vars

             ! x direction
             check_comp_cells_x:IF ( comp_cells_x .GT. 1 ) THEN

                ! west boundary
                check_x_boundary:IF (j.EQ.1) THEN

                   IF ( bcW(i)%flag .EQ. 0 ) THEN

                      x_stencil(1) = x_stag(1)
                      x_stencil(2:3) = x_comp(1:2)

                      qrec_stencil(1) = bcW(i)%value
                      qrec_stencil(2:3) = qrec(i,1:2,k)

                      CALL limit( qrec_stencil , x_stencil , limiter(i) ,       &
                           qrec_prime_x(i) ) 

                   ELSEIF ( bcW(i)%flag .EQ. 1 ) THEN

                      qrec_prime_x(i) = bcW(i)%value

                   ELSEIF ( bcW(i)%flag .EQ. 2 ) THEN

                      qrec_prime_x(i) = ( qrec(i,2,k) - qrec(i,1,k) ) / dx

                   END IF

                   !east boundary
                ELSEIF (j.EQ.comp_cells_x) THEN

                   IF ( bcE(i)%flag .EQ. 0 ) THEN

                      qrec_stencil(3) = bcE(i)%value
                      qrec_stencil(1:2) = qrec(i,comp_cells_x-1:comp_cells_x,k)

                      x_stencil(3) = x_stag(comp_interfaces_x)
                      x_stencil(1:2) = x_comp(comp_cells_x-1:comp_cells_x)

                      CALL limit( qrec_stencil , x_stencil , limiter(i) ,       &
                           qrec_prime_x(i) ) 

                   ELSEIF ( bcE(i)%flag .EQ. 1 ) THEN

                      qrec_prime_x(i) = bcE(i)%value

                   ELSEIF ( bcE(i)%flag .EQ. 2 ) THEN

                      qrec_prime_x(i) = ( qrec(i,comp_cells_x,k) -                 &
                           qrec(i,comp_cells_x-1,k) ) / dx

                   END IF

                   ! internal x cells
                ELSE

                   x_stencil(1:3) = x_comp(j-1:j+1)
                   qrec_stencil(1:3) = qrec(i,j-1:j+1,k)

                   ! correction for radial source inlet x-interfaces values 
                   ! used for the linear reconstruction
                   IF ( radial_source_flag .AND. ( source_cell(j,k).EQ.2 ) ) THEN

                      IF ( sourceE(j,k) ) THEN

                        CALL eval_source_bdry( sourceE_vect_x(j,k) ,            &
                             sourceE_vect_y(j,k) , source_bdry )

                        x_stencil(3) = x_stag(j+1)
                        qrec_stencil(3) = source_bdry(i)

                      ELSEIF ( sourceW(j,k) ) THEN

                        CALL eval_source_bdry( sourceW_vect_x(j,k) ,            &
                             sourceW_vect_y(j,k) , source_bdry )

                        x_stencil(1) = x_stag(j)
                        qrec_stencil(1) = source_bdry(i)

                      END IF

                   END IF

                   CALL limit( qrec_stencil , x_stencil , limiter(i) ,          &
                        qrec_prime_x(i) )

                ENDIF check_x_boundary

                qrecW(i) = qrec(i,j,k) - reconstr_coeff * dx2 * qrec_prime_x(i)
                qrecE(i) = qrec(i,j,k) + reconstr_coeff * dx2 * qrec_prime_x(i)
                
             END IF check_comp_cells_x

             ! y-direction
             check_comp_cells_y:IF ( comp_cells_y .GT. 1 ) THEN

                ! South boundary
                check_y_boundary:IF (k.EQ.1) THEN

                   IF ( bcS(i)%flag .EQ. 0 ) THEN

                      qrec_stencil(1) = bcS(i)%value
                      qrec_stencil(2:3) = qrec(i,j,1:2)

                      y_stencil(1) = y_stag(1)
                      y_stencil(2:3) = y_comp(1:2)

                      CALL limit( qrec_stencil , y_stencil , limiter(i) ,       &
                           qrec_prime_y(i) ) 

                   ELSEIF ( bcS(i)%flag .EQ. 1 ) THEN

                      qrec_prime_y(i) = bcS(i)%value

                   ELSEIF ( bcS(i)%flag .EQ. 2 ) THEN

                      qrec_prime_y(i) = ( qrec(i,j,2) - qrec(i,j,1) ) / dy 

                   END IF

                   ! North boundary
                ELSEIF ( k .EQ. comp_cells_y ) THEN

                   IF ( bcN(i)%flag .EQ. 0 ) THEN

                      qrec_stencil(3) = bcN(i)%value
                      qrec_stencil(1:2) = qrec(i,j,comp_cells_y-1:comp_cells_y)

                      y_stencil(3) = y_stag(comp_interfaces_y)
                      y_stencil(1:2) = y_comp(comp_cells_y-1:comp_cells_y)

                      CALL limit( qrec_stencil , y_stencil , limiter(i) ,       &
                           qrec_prime_y(i) ) 

                   ELSEIF ( bcN(i)%flag .EQ. 1 ) THEN

                      qrec_prime_y(i) = bcN(i)%value

                   ELSEIF ( bcN(i)%flag .EQ. 2 ) THEN

                      qrec_prime_y(i) = ( qrec(i,j,comp_cells_y) -                 &
                           qrec(i,j,comp_cells_y-1) ) / dy 

                   END IF

                   ! Internal y cells
                ELSE

                   y_stencil(1:3) = y_comp(k-1:k+1)
                   qrec_stencil = qrec(i,j,k-1:k+1)

                   ! correction for radial source inlet y-interfaces
                   ! used for the linear reconstruction
                   IF ( radial_source_flag .AND. ( source_cell(j,k).EQ.2 ) ) THEN

                      IF ( sourceS(j,k) ) THEN

                        CALL eval_source_bdry( sourceS_vect_x(j,k) ,            &
                             sourceS_vect_y(j,k) , source_bdry )

                        y_stencil(1) = y_stag(k)
                        qrec_stencil(1) = source_bdry(i)

                      ELSEIF ( sourceN(j,k) ) THEN

                        CALL eval_source_bdry( sourceN_vect_x(j,k) ,            &
                             sourceN_vect_y(j,k) , source_bdry )

                        x_stencil(3) = y_stag(k+1)
                        qrec_stencil(3) = source_bdry(i)

                      END IF

                   END IF

                   CALL limit( qrec_stencil , y_stencil , limiter(i) ,          &
                        qrec_prime_y(i) )

                ENDIF check_y_boundary

                qrecS(i) = qrec(i,j,k) - reconstr_coeff * dy2 * qrec_prime_y(i)
                qrecN(i) = qrec(i,j,k) + reconstr_coeff * dy2 * qrec_prime_y(i)

             ENDIF check_comp_cells_y

          ENDDO vars_loop


          IF ( ( j .GT. 1 ) .AND. ( j .LT. comp_cells_x ) ) THEN
             
             ! correction on the reconstruction slope in order to keep u,v at  
             ! the W,E interfaces of internal cells limited (no new max and 
             ! min created) 
  
             ! compute the values of u,v at the cell centers and W,E interfaces
             CALL qp_to_qp2( qrec(1:n_vars,j-1,k) , B_cent(j-1,k) , qp2recL ) 
             CALL qp_to_qp2( qrec(1:n_vars,j,k) , B_cent(j,k) , qp2recC ) 
             CALL qp_to_qp2( qrec(1:n_vars,j+1,k) , B_cent(j+1,k) , qp2recR ) 

             qrec_stencil(1) = qp2recL(1)
             qrec_stencil(2) = qp2recC(1)
             qrec_stencil(3) = qp2recR(1)

!!$             CALL limit( qrec_stencil , x_stencil , limiter(1) ,          &
!!$                        qp2rec_prime_x(1) )
!!$
!!$             qp2recW(1) = qp2recC(1) - reconstr_coeff * dx2 * qp2rec_prime_x(1)
!!$             qp2recE(1) = qp2recC(1) + reconstr_coeff * dx2 * qp2rec_prime_x(1)
!!$
!!$             qrecW(1) = qp2recW(1) - B_interfaceR(j,k)
!!$             qrecE(1) = qp2recE(1) - B_interfaceL(j+1,k)

             CALL qp_to_qp2( qrecW(1:n_vars) , B_interfaceR(j,k), qp2recW ) 
             CALL qp_to_qp2( qrecE(1:n_vars) , B_interfaceL(j+1,k) , qp2recE ) 

             DO i=2,3
           
                IF ( qp2recW(i) .GT. 0.D0 ) THEN
                   
                   gamma_W =  MAX( DABS( qp2recL(i) ) , DABS( qp2recC(i) ) )    &
                        / DABS( qp2recW(i) )
                
                ELSE

                   gamma_W = 1.D0

                END IF

                IF ( qp2recE(i) .GT. 0.D0 ) THEN

                   gamma_E =  MAX( DABS( qp2recC(i) ) , DABS( qp2recR(i) ) )       &
                        / DABS( qp2recE(i) )
                
                ELSE

                   gamma_E = 1.D0

                END IF

                gamma = MIN( gamma_W,gamma_E,1.D0)
                
                qrecW(i) = qrec(i,j,k) - gamma * dx2 * qrec_prime_x(i)
                qrecE(i) = qrec(i,j,k) + gamma * dx2 * qrec_prime_x(i)
                
             END DO
             
          END IF
          
          IF ( ( k.GT.1 ) .AND. ( k .LT. comp_cells_y ) ) THEN

             ! correction on the reconstruction slope in order to keep u,v at  
             ! the S,N interfaces of internal cells limited (no new max and 
             ! min created) 
 
             ! compute the values of u,v at the cell centers and S,N interfaces
             CALL qp_to_qp2( qrec(1:n_vars,j,k-1) , B_cent(j,k-1) , qp2recB ) 
             CALL qp_to_qp2( qrecS(1:n_vars) , B_interfaceT(j,k) , qp2recS ) 
             CALL qp_to_qp2( qrec(1:n_vars,j,k) , B_cent(j,k) , qp2recC ) 
             CALL qp_to_qp2( qrecN(1:n_vars) , B_interfaceB(j,k+1) , qp2recN ) 
             CALL qp_to_qp2( qrec(1:n_vars,j,k+1) , B_cent(j,k+1) , qp2recT ) 

             DO i=2,3

                IF ( qp2recS(i) .GT. 0.D0 ) THEN

                   gamma_S =  MAX( DABS( qp2recB(i) ) , DABS( qp2recC(i) ) )    &
                        / DABS( qp2recS(i) )
                
                ELSE

                   gamma_S = 1.D0

                END IF

                IF ( qp2recN(i) .GT. 0.D0 ) THEN

                   gamma_N =  MAX( DABS( qp2recC(i) ) , DABS( qp2recT(i) ) )    &
                        / DABS( qp2recN(i) )
 
                ELSE

                   gamma_N = 1.D0 

                END IF

                gamma = MIN( gamma_S,gamma_N,1.D0)
             
                qrecS(i) = qrec(i,j,k) - gamma * dy2 * qrec_prime_y(i)
                qrecN(i) = qrec(i,j,k) + gamma * dy2 * qrec_prime_y(i)
             
             END DO

          END IF
          
          
          
          IF ( comp_cells_x .GT. 1 ) THEN

             IF ( ( j .GT. 1 ) .AND. ( j .LT. comp_cells_x ) ) THEN

                IF ( q_expl(1,j,k) .EQ. 0.D0 ) THEN

                   IF ( ( .NOT. radial_source_flag ) .OR.                       &
                        ( ( radial_source_flag ) .AND.                          &
                        ( source_cell(j,k) .EQ. 0 ) ) ) THEN
                   
                      ! In the internal cell, if thickness h is 0 at the center
                      ! of the cell, then all the variables are 0 at the center
                      ! and at the interfaces (no conversion back is needed from
                      ! reconstructed to conservative)
                      q_interfaceR(:,j,k) = 0.D0
                      q_interfaceL(:,j+1,k) = 0.D0
                      
                      qp_interfaceR(1:3,j,k) = 0.D0
                      qp_interfaceR(4:n_vars,j,k) = qrecW(4:n_vars)
                      qp_interfaceL(1:3,j+1,k) = 0.D0
                      qp_interfaceL(4:n_vars,j+1,k) = qrecE(4:n_vars)
                      
                   END IF

                END IF

             END IF
                
                IF ( j.EQ.1 ) THEN
                   
                   ! Dirichelet boundary condition at the west of the domain
                   DO i=1,n_vars
                      
                      IF ( bcW(i)%flag .EQ. 0 ) THEN
                         
                         qrecW(i) = bcW(i)%value 

                      END IF

                   ENDDO

                ELSEIF ( j.EQ.comp_cells_x ) THEN

                   ! Dirichelet boundary condition at the east of the domain
                   DO i=1,n_vars

                      IF ( bcE(i)%flag .EQ. 0 ) THEN

                         qrecE(i) = bcE(i)%value 

                      END IF

                   ENDDO

                ELSE

                   ! correction for radial source inlet x-interfaces:
                   ! the physical variables at the x-interfaces qrecW or
                   ! qrecE are computed from the radial inlet values
                   IF ( radial_source_flag .AND. ( source_cell(j,k).EQ.2 ) ) THEN

                      IF ( sourceE(j,k) ) THEN

                        CALL eval_source_bdry( sourceE_vect_x(j,k) ,            &
                             sourceE_vect_y(j,k) , source_bdry )

                        qrecE(1:n_vars) = source_bdry(1:n_vars)

                      ELSEIF ( sourceW(j,k) ) THEN

                        CALL eval_source_bdry( sourceW_vect_x(j,k) ,            &
                             sourceW_vect_y(j,k) , source_bdry )

                        qrecW(1:n_vars) = source_bdry(1:n_vars)
                        
                      END IF

                   END IF

                END IF
                   
                CALL qp_to_qc( qrecW,B_interfaceR(j,k),q_interfaceR(:,j,k) )
                CALL qp_to_qc( qrecE,B_interfaceL(j+1,k),q_interfaceL(:,j+1,k) )

                qp_interfaceR(1:n_vars,j,k) = qrecW(1:n_vars)
                qp_interfaceL(1:n_vars,j+1,k) = qrecE(1:n_vars)

                IF ( j.EQ.1 ) THEN

                   ! Interface value at the left of first x-interface (external)
                   q_interfaceL(:,j,k) = q_interfaceR(:,j,k)
                   qp_interfaceL(:,j,k) = qp_interfaceR(:,j,k)

                   !WRITE(*,*) 'q_interfaceL(:,j,k)',q_interfaceL(:,j,k)
                   !READ(*,*)

                ELSEIF ( j.EQ.comp_cells_x ) THEN

                   ! Interface value at the right of last x-interface (external)
                   q_interfaceR(:,j+1,k) = q_interfaceL(:,j+1,k)
                   qp_interfaceR(:,j+1,k) = qp_interfaceL(:,j+1,k)

                ELSE

                   IF ( radial_source_flag .AND. ( source_cell(j,k) .EQ. 2 ) ) THEN

                      IF ( sourceE(j,k) ) THEN

                         q_interfaceR(:,j+1,k) = q_interfaceL(:,j+1,k)
                         qp_interfaceR(:,j+1,k) = qp_interfaceL(:,j+1,k)
                         
                      ELSEIF ( sourceW(j,k) ) THEN
                      
                         q_interfaceL(:,j,k) = q_interfaceR(:,j,k)
                         qp_interfaceL(:,j,k) = qp_interfaceR(:,j,k)
                         
                      END IF

                   END IF

                END IF

          ELSE
             
             ! for case comp_cells_x = 1 
             q_interfaceR(1:n_vars,j,k) = q_expl(1:n_vars,j,k)
             q_interfaceL(1:n_vars,j+1,k) = q_expl(1:n_vars,j,k)


             qp_interfaceR(1:n_vars,j,k) = qrec(1:n_vars,j,k)
             qp_interfaceL(1:n_vars,j+1,k) = qrec(1:n_vars,j,k)
            
          END IF
          
          IF ( comp_cells_y .GT. 1 ) THEN
             
             IF ( ( k .GT. 1 ) .AND. ( k .LT. comp_cells_y ) ) THEN
                
                IF ( q_expl(1,j,k) .EQ. 0.D0 ) THEN
                   
                   IF ( ( .NOT. radial_source_flag ) .OR.   &
                        ( ( radial_source_flag ) .AND. ( source_cell(j,k) .EQ. 0 ) ) ) THEN
                      
                      ! In the internal cell, if thickness h is 0 at the center
                      ! of the cell, then all the variables are 0 at the center
                      ! and at the interfaces (no conversion back is needed from
                      ! reconstructed to conservative)
                      
                      q_interfaceT(:,j,k) = 0.D0
                      q_interfaceB(:,j,k+1) = 0.D0
                      
                      qp_interfaceT(1:3,j,k) = 0.D0
                      qp_interfaceT(4:n_vars,j,k) = qrecS(4:n_vars)
                      qp_interfaceB(1:3,j,k+1) = 0.D0
                      qp_interfaceB(4:n_vars,j,k+1) = qrecN(4:n_vars)
                      
                   END IF
                   
                END IF
                
             END IF
                
                IF ( k .EQ. 1 ) THEN

                   ! Dirichelet boundary condition at the south of the domain
                   DO i=1,n_vars
                      
                      IF ( bcS(i)%flag .EQ. 0 ) THEN
                         
                         qrecS(i) = bcS(i)%value 
                         
                      END IF
                      
                   ENDDO
                   
                ELSEIF ( k .EQ. comp_cells_y ) THEN
                   
                   ! Dirichelet boundary condition at the north of the domain
                   DO i=1,n_vars

                      IF ( bcN(i)%flag .EQ. 0 ) THEN

                         qrecN(i) = bcN(i)%value 

                      END IF

                   ENDDO
        
                ELSE

                   IF ( radial_source_flag .AND. ( source_cell(j,k) .EQ. 2 ) ) THEN

                      IF ( sourceS(j,k) ) THEN

                        CALL eval_source_bdry( sourceS_vect_x(j,k) ,            &
                             sourceS_vect_y(j,k) , source_bdry )

                        qrecS(1:n_vars) = source_bdry(1:n_vars)

                      ELSEIF ( sourceN(j,k) ) THEN

                        CALL eval_source_bdry( sourceN_vect_x(j,k) ,            &
                             sourceN_vect_y(j,k) , source_bdry )

                        qrecN(1:n_vars) = source_bdry(1:n_vars)
                        
                      END IF

                   END IF

                END IF

                CALL qp_to_qc( qrecS, B_interfaceT(j,k), q_interfaceT(:,j,k))
                CALL qp_to_qc( qrecN, B_interfaceB(j,k+1), q_interfaceB(:,j,k+1))

                qp_interfaceT(:,j,k) = qrecS
                qp_interfaceB(:,j,k+1) = qrecN

                IF ( k .EQ. 1 ) THEN

                   ! Interface value at the bottom of first y-interface (external)
                   q_interfaceB(:,j,k) = q_interfaceT(:,j,k)
                   qp_interfaceB(:,j,k) = qp_interfaceT(:,j,k)

                ELSEIF ( k .EQ. comp_cells_y ) THEN

                   ! Interface value at the top of last y-interface (external)
                   q_interfaceT(:,j,k+1) = q_interfaceB(:,j,k+1)
                   qp_interfaceT(:,j,k+1) = qp_interfaceB(:,j,k+1)

                ELSE

                   IF ( radial_source_flag .AND. ( source_cell(j,k) .EQ. 2 ) ) THEN

                      IF ( sourceS(j,k) ) THEN

                         q_interfaceB(:,j,k) = q_interfaceT(:,j,k)
                         qp_interfaceB(:,j,k) = qp_interfaceT(:,j,k)

                      ELSEIF ( sourceN(j,k) ) THEN

                         q_interfaceT(:,j,k+1) = q_interfaceB(:,j,k+1)
                         qp_interfaceT(:,j,k+1) = qp_interfaceB(:,j,k+1)

                      END IF

                   END IF

                END IF

          ELSE
          
             ! comp_cells_y = 1
             q_interfaceT(:,j,k) = q_expl(:,j,k)
             q_interfaceB(:,j,k+1) = q_expl(:,j,k)
  
             qp_interfaceT(:,j,k) = qrec(:,j,k)
             qp_interfaceB(:,j,k+1) = qrec(:,j,k)
                        
          END IF
          
       END IF

       END DO x_loop

    END DO y_loop

  END SUBROUTINE reconstruction

  
  !******************************************************************************
  !> \brief Characteristic speeds
  !
  !> This subroutine evaluates the largest characteristic speed at the
  !> cells interfaces from the reconstructed states.
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 2019/11/11
  !******************************************************************************

  SUBROUTINE eval_speeds

    ! External procedures
    USE constitutive_2d, ONLY : eval_local_speeds_x, eval_local_speeds_y 

    IMPLICIT NONE

    REAL*8 :: abslambdaL_min(n_vars) , abslambdaL_max(n_vars)
    REAL*8 :: abslambdaR_min(n_vars) , abslambdaR_max(n_vars)
    REAL*8 :: abslambdaB_min(n_vars) , abslambdaB_max(n_vars)
    REAL*8 :: abslambdaT_min(n_vars) , abslambdaT_max(n_vars)
    REAL*8 :: min_r(n_vars) , max_r(n_vars)

    INTEGER :: j,k

    IF ( comp_cells_x .GT. 1 ) THEN
    
       DO j = 1,comp_interfaces_x
          
          DO k = 1, comp_cells_y
             
             CALL eval_local_speeds_x( qp_interfaceL(:,j,k) , abslambdaL_min ,  &
                  abslambdaL_max )
             
             CALL eval_local_speeds_x( qp_interfaceR(:,j,k) , abslambdaR_min ,  &
                  abslambdaR_max )

             min_r = MIN(abslambdaL_min , abslambdaR_min , 0.0D0)
             max_r = MAX(abslambdaL_max , abslambdaR_max , 0.0D0)
             
             a_interface_xNeg(:,j,k) = min_r
             a_interface_xPos(:,j,k) = max_r
             
          ENDDO
          
       END DO

    END IF

    IF ( comp_cells_y .GT. 1 ) THEN
    
       DO j = 1,comp_cells_x
          
          DO k = 1,comp_interfaces_y
             
             CALL eval_local_speeds_y( qp_interfaceB(:,j,k) , abslambdaB_min ,  &
                  abslambdaB_max )
             
             CALL eval_local_speeds_y( qp_interfaceT(:,j,k) , abslambdaT_min ,  &
                  abslambdaT_max )
            
             min_r = MIN(abslambdaB_min , abslambdaT_min , 0.0D0)
             max_r = MAX(abslambdaB_max , abslambdaT_max , 0.0D0)
             
             a_interface_yNeg(:,j,k) = min_r
             a_interface_yPos(:,j,k) = max_r
             
          ENDDO
          
       END DO

    END IF
       
  END SUBROUTINE eval_speeds

END MODULE solver_2d
