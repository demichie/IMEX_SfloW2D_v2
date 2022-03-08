!*********************************************************************
!> \brief Grid module
!
!> This module contains the variables and the subroutines related to 
!> the computational grid
!*********************************************************************
MODULE geometry_2d

  USE parameters_2d, ONLY : wp , sp
  USE parameters_2d, ONLY : verbose_level

  IMPLICIT NONE

  !> Location of the centers (x) of the control volume of the domain
  REAL(wp), ALLOCATABLE :: x_comp(:)

  !> Location of the boundaries (x) of the control volumes of the domain
  REAL(wp), ALLOCATABLE :: x_stag(:)

  !> Location of the centers (y) of the control volume of the domain
  REAL(wp), ALLOCATABLE :: y_comp(:)

  !> Location of the boundaries (x) of the control volumes of the domain
  REAL(wp), ALLOCATABLE :: y_stag(:)

  !> Topography at the centers of the control volumes 
  REAL(wp), ALLOCATABLE :: B_cent(:,:)

  !> Topography at the centers of the control volumes 
  REAL(wp), ALLOCATABLE :: B_cent_extended(:,:)

  LOGICAL, ALLOCATABLE :: B_nodata(:,:)

  INTEGER, ALLOCATABLE :: B_zone(:,:)


  !> Topography slope (x direction) at the centers of the control volumes 
  REAL(wp), ALLOCATABLE :: B_prime_x(:,:)

  !> Topography 2nd x-derivative at the centers of the control volumes 
  REAL(wp), ALLOCATABLE :: B_second_xx(:,:)

  !> Topography slope (y direction) at the centers of the control volumes 
  REAL(wp), ALLOCATABLE :: B_prime_y(:,:)

  !> Topography 2nd y-derivative at the centers of the control volumes 
  REAL(wp), ALLOCATABLE :: B_second_yy(:,:)

  !> Topography 2nd xy-derivative at the centers of the control volumes 
  REAL(wp), ALLOCATABLE :: B_second_xy(:,:)

  !> Solution in ascii grid format (ESRI)
  REAL(wp), ALLOCATABLE :: grid_output(:,:)

  !> Integer solution in ascii grid format (ESRI)
  INTEGER, ALLOCATABLE :: grid_output_int(:,:)

  !> gravity coefficient (accounting for slope) at cell centers
  REAL(wp), ALLOCATABLE :: grav_coeff(:,:)

  !> 1st x-derivative of gravity coefficient
  REAL(wp), ALLOCATABLE :: d_grav_coeff_dx(:,:)

  !> 1st y-derivative of gravity coefficient
  REAL(wp), ALLOCATABLE :: d_grav_coeff_dy(:,:)

  !> modified gravity at cell x-faces
  REAL(wp), ALLOCATABLE :: grav_coeff_stag_x(:,:)

  !> modified gravity at cell y-faces
  REAL(wp), ALLOCATABLE :: grav_coeff_stag_y(:,:)

  !> curvature wrt mixed directions for each cell
  REAL(wp), ALLOCATABLE :: curv_xy(:,:)

  !> deposit for the different classes
  REAL(wp), ALLOCATABLE :: deposit(:,:,:)

  !> total deposit 
  REAL(wp), ALLOCATABLE :: deposit_tot(:,:)

  !> erosion for the different classes
  REAL(wp), ALLOCATABLE :: erosion(:,:,:)

  !> total erosion 
  REAL(wp), ALLOCATABLE :: erosion_tot(:,:)

  !> erosdible substrate for the different classes
  REAL(wp), ALLOCATABLE :: erodible(:,:,:)

  REAL(wp), ALLOCATABLE :: topography_profile(:,:,:)

  REAL(wp) :: nodata_topo

  INTEGER, ALLOCATABLE :: source_cell(:,:)
  LOGICAL, ALLOCATABLE :: sourceE(:,:)
  LOGICAL, ALLOCATABLE :: sourceW(:,:)
  LOGICAL, ALLOCATABLE :: sourceS(:,:)
  LOGICAL, ALLOCATABLE :: sourceN(:,:)

  REAL(wp), ALLOCATABLE :: dist_sourceE(:,:)
  REAL(wp), ALLOCATABLE :: dist_sourceW(:,:)
  REAL(wp), ALLOCATABLE :: dist_sourceS(:,:)
  REAL(wp), ALLOCATABLE :: dist_sourceN(:,:)

  REAL(wp), ALLOCATABLE :: sourceE_vect_x(:,:)
  REAL(wp), ALLOCATABLE :: sourceE_vect_y(:,:)

  REAL(wp), ALLOCATABLE :: sourceW_vect_x(:,:)
  REAL(wp), ALLOCATABLE :: sourceW_vect_y(:,:)

  REAL(wp), ALLOCATABLE :: sourceS_vect_x(:,:)
  REAL(wp), ALLOCATABLE :: sourceS_vect_y(:,:)

  REAL(wp), ALLOCATABLE :: sourceN_vect_x(:,:)
  REAL(wp), ALLOCATABLE :: sourceN_vect_y(:,:)

  REAL(wp), ALLOCATABLE :: cell_source_fractions(:,:)

  REAL(wp) :: pi_g

  INTEGER :: n_topography_profile_x, n_topography_profile_y

  REAL(wp) :: dx                 !< Control volumes size
  REAL(wp) :: x0                 !< Left of the physical domain
  REAL(wp) :: xN                 !< Right of the physical domain
  REAL(wp) :: dy                 !< Control volumes size
  REAL(wp) :: y0                 !< Bottom of the physical domain
  REAL(wp) :: yN                 !< Top of the physical domain
  REAL(wp) :: dx2                !< Half x Control volumes size
  REAL(wp) :: dy2                !< Half y Control volumes size

  REAL(wp) :: one_by_dx
  REAL(wp) :: one_by_dy

  INTEGER :: comp_cells_x      !< Number of control volumes x in the comp. domain
  INTEGER :: comp_interfaces_x !< Number of interfaces (comp_cells_x+1)
  INTEGER :: comp_cells_y      !< Number of control volumes y in the comp. domain
  INTEGER :: comp_interfaces_y !< Number of interfaces (comp_cells_y+1)
  REAL(wp) :: cell_size
  INTEGER :: comp_cells_xy

CONTAINS

  !******************************************************************************
  !> \brief Finite volume grid initialization
  !
  !> This subroutine initialize the grids for the finite volume solver.
  !> \date 16/08/2011
  !******************************************************************************

  SUBROUTINE init_grid

    USE parameters_2d, ONLY: eps_sing , eps_sing4
    USE parameters_2d, ONLY : bottom_radial_source_flag
    USE parameters_2d, ONLY : x_source , y_source , r_source
    USE parameters_2d, ONLY : liquid_vaporization_flag

    IMPLICIT none

    INTEGER j,k      !> loop counter

    comp_interfaces_x = comp_cells_x+1
    comp_interfaces_y = comp_cells_y+1

    ALLOCATE( x_comp(comp_cells_x) )
    ALLOCATE( x_stag(comp_interfaces_x) )
    ALLOCATE( y_comp(comp_cells_y) )
    ALLOCATE( y_stag(comp_interfaces_y) )

    ALLOCATE( source_cell(comp_cells_x,comp_cells_y) )

    ! cell where are equations are solved
    source_cell(1:comp_cells_x,1:comp_cells_y) = 0


    ALLOCATE( sourceE(comp_cells_x,comp_cells_y) )
    ALLOCATE( sourceW(comp_cells_x,comp_cells_y) )
    ALLOCATE( sourceN(comp_cells_x,comp_cells_y) )
    ALLOCATE( sourceS(comp_cells_x,comp_cells_y) )

    ALLOCATE( dist_sourceE(comp_cells_x,comp_cells_y) )
    ALLOCATE( dist_sourceW(comp_cells_x,comp_cells_y) )
    ALLOCATE( dist_sourceN(comp_cells_x,comp_cells_y) )
    ALLOCATE( dist_sourceS(comp_cells_x,comp_cells_y) )

    ALLOCATE( sourceE_vect_x(comp_cells_x,comp_cells_y) )
    ALLOCATE( sourceE_vect_y(comp_cells_x,comp_cells_y) )
    ALLOCATE( sourceW_vect_x(comp_cells_x,comp_cells_y) )
    ALLOCATE( sourceW_vect_y(comp_cells_x,comp_cells_y) )
    ALLOCATE( sourceS_vect_x(comp_cells_x,comp_cells_y) )
    ALLOCATE( sourceS_vect_y(comp_cells_x,comp_cells_y) )
    ALLOCATE( sourceN_vect_x(comp_cells_x,comp_cells_y) )
    ALLOCATE( sourceN_vect_y(comp_cells_x,comp_cells_y) )

    ALLOCATE( B_cent(comp_cells_x,comp_cells_y) )
    ALLOCATE( B_cent_extended(comp_cells_x+2,comp_cells_y+2) )

    ALLOCATE( B_nodata(comp_cells_x,comp_cells_y) )
    ALLOCATE( B_prime_x(comp_cells_x,comp_cells_y) )
    ALLOCATE( B_prime_y(comp_cells_x,comp_cells_y) )

    ALLOCATE( B_second_xx(comp_cells_x,comp_cells_y) )
    ALLOCATE( B_second_yy(comp_cells_x,comp_cells_y) )
    ALLOCATE( B_second_xy(comp_cells_x,comp_cells_y) )

    ALLOCATE( grid_output(comp_cells_x,comp_cells_y) )
    ALLOCATE( grid_output_int(comp_cells_x,comp_cells_y) )

    ALLOCATE( grav_coeff(comp_cells_x,comp_cells_y) )
    ALLOCATE( d_grav_coeff_dx(comp_cells_x,comp_cells_y) )
    ALLOCATE( d_grav_coeff_dy(comp_cells_x,comp_cells_y) )

    ALLOCATE( grav_coeff_stag_x(comp_interfaces_x,comp_cells_y) )
    ALLOCATE( grav_coeff_stag_y(comp_cells_x,comp_interfaces_y) )

    ALLOCATE( cell_source_fractions(comp_cells_x,comp_cells_y) )

    IF ( comp_cells_x .GT. 1 ) THEN

       dx = cell_size

    ELSE

       dx = 1.0_wp

    END IF

    IF ( comp_cells_y .GT. 1 ) THEN

       dy = cell_size

    ELSE

       dy = 1.0_wp

    END IF

    xN = x0 + comp_cells_x * dx
    yN = y0 + comp_cells_y * dy

    dx2 = dx / 2.0_wp
    dy2 = dy / 2.0_wp

    one_by_dx = 1.0_wp / dx
    one_by_dy = 1.0_wp / dy


    IF ( wp .EQ. sp ) THEN

       eps_sing=MIN(MIN( dx ** 4.0_wp,dy ** 4.0_wp ),1.0E-6_wp)

    ELSE

       eps_sing=MIN(MIN( dx ** 4.0_wp,dy ** 4.0_wp ),1.0E-10_wp)

    END IF

    eps_sing4 = eps_sing**4

    IF ( verbose_level .GE. 1 ) WRITE(*,*) 'eps_sing = ',eps_sing

    DO j=1,comp_interfaces_x

       x_stag(j) = x0 + (j-1) * dx

    END DO

    DO k=1,comp_interfaces_y

       y_stag(k) = y0 + (k-1) * dy

    END DO

    DO j=1,comp_cells_x

       x_comp(j) = 0.5_wp * ( x_stag(j) + x_stag(j+1) )

    END DO

    DO k=1,comp_cells_y

       y_comp(k) = 0.5_wp * ( y_stag(k) + y_stag(k+1) )

    END DO

    DO k=1,comp_cells_y

       DO j=1,comp_cells_x

          CALL interp_2d_nodata( topography_profile(1,:,:) ,                    &
               topography_profile(2,:,:), topography_profile(3,:,:) ,           &
               x_comp(j), y_comp(k) , B_nodata(j,k) )

       END DO

    ENDDO

    topography_profile(3,:,:) = MAX(0.0_wp,topography_profile(3,:,:))

    DO k=1,comp_cells_y

       DO j=1,comp_cells_x

          CALL interp_2d_scalar( topography_profile(1,:,:) ,                    &
               topography_profile(2,:,:), topography_profile(3,:,:) ,           &
               x_comp(j), y_comp(k) , B_cent(j,k) )

       END DO

    ENDDO

    ! This subroutine compute the partial derivatives of the topography at the
    ! cell centers and 
    CALL topography_reconstruction

    ALLOCATE(  B_zone(comp_cells_x,comp_cells_y) )

    B_zone = 0

    IF ( liquid_vaporization_flag ) CALL topography_zones

    pi_g = 4.0_wp * ATAN(1.0_wp)

    IF ( bottom_radial_source_flag ) THEN

       CALL compute_cell_fract(x_source,y_source,r_source,cell_source_fractions)

    END IF

    RETURN

  END SUBROUTINE init_grid

  !******************************************************************************
  !> \brief Topography zone identification
  !
  !> This subroutine search for the connected zones where topography elevation
  !> corresponds to a fixed value assigned in the input file (water_level). 
  !> A value 1 to is assigned to variable B_zone in the cells belonging to the
  !> largest connected area. 
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 2021/07/21
  !******************************************************************************


  SUBROUTINE topography_zones

    USE parameters_2D, ONLY: water_level

    IMPLICIT NONE

    INTEGER :: i,j,k

    INTEGER :: zone_counter
    INTEGER :: equi_list(0:1000)
    INTEGER :: equi_list_new(0:1000)
    INTEGER :: zone , zone_max
    INTEGER :: zone_cells(1:1000)

    WHERE ( ABS( B_cent - water_level ) .LE. 1.0E-1_wp )

       B_zone = -1

    END WHERE

    DO i=0,1000

       equi_list(i) = i

    END DO

    zone_counter = 0

    y_loop:DO k = 1,comp_cells_y

       x_loop:DO j = 1,comp_cells_x

          IF ( B_zone(j,k) .EQ. -1 ) THEN

             IF ( ( k .GT. 1 ) .AND. ( j.GT.1) ) THEN

                B_zone(j,k) = MAX( B_zone(j-1,k) , B_zone(j,k-1) )

                IF ( ( B_zone(j-1,k) .NE. B_zone(j,k-1) ) .AND.  &
                     (MIN( B_zone(j-1,k) , B_zone(j,k-1) ) .GT. 0 ) ) THEN

                   B_zone(j,k) = MIN( B_zone(j-1,k) , B_zone(j,k-1) )

                   equi_list(MAX( B_zone(j-1,k) , B_zone(j,k-1) )) = &
                        MIN( B_zone(j-1,k) , B_zone(j,k-1) )

                END IF

             ELSEIF ( j .GT. 1 ) THEN

                B_zone(j,k) = B_zone(j-1,k)

             ELSEIF ( k .GT. 1 ) THEN

                B_zone(j,k) = B_zone(j,k-1)

             END IF

             IF ( B_zone(j,k) .LE. 0 ) THEN

                zone_counter = zone_counter + 1
                B_zone(j,k) = zone_counter

             END IF

          END IF

       END DO x_loop

    END DO y_loop

    ! Search for smallest value in the equivalence list
    DO i=zone_counter,1,-1

       DO WHILE ( equi_list(equi_list(i)) .NE. equi_list(i) )

          equi_list(i) = equi_list(equi_list(i))

       END DO

    END DO

    equi_list_new(0:zone_counter) = 0
    i = 0

    DO zone = 1, zone_counter

       IF ( equi_list(zone) == zone ) THEN

          i = i + 1
          equi_list_new(zone) = i

       END IF
    END DO

    zone_counter = i

    zone_cells = 0

    !  Replace the labels by consecutive labels.
    y_loop2:DO k = 1,comp_cells_y

       x_loop2:DO j = 1,comp_cells_x

          B_zone(j,k) = equi_list_new(equi_list(B_zone(j,k)))

          IF ( B_zone(j,k) .GT. 0 ) THEN

             zone_cells( B_zone(j,k) ) = zone_cells( B_zone(j,k) ) + 1

          END IF

       END DO x_loop2

    END DO y_loop2

    WRITE(*,*) 'Number of cells of each conneceted area:'
    WRITE(*,*) zone_cells(1:zone_counter)

    zone_max = MAXLOC(zone_cells,1)

    y_loop3:DO k = 1,comp_cells_y

       x_loop3:DO j = 1,comp_cells_x

          IF ( B_zone(j,k) .NE. zone_max ) THEN

             B_zone(j,k) = 0

          ELSE

             B_zone(j,k) = 1

          END IF

       END DO x_loop3

    END DO y_loop3

    RETURN

  END SUBROUTINE topography_zones

  !******************************************************************************
  !> \brief Topography slope recontruction
  !
  !> In this subroutine a linear reconstruction with slope limiters is
  !> applied to compute dB_dx and dB_dy at the cell centers. The second
  !> derivative of slope and the correction factor to the gravity to account
  !> for the slope gradient are also computer here
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 2019/11/08
  !******************************************************************************

  SUBROUTINE topography_reconstruction

    USE parameters_2d, ONLY : limiter
    USE parameters_2d, ONLY : slope_correction_flag

    IMPLICIT NONE

    REAL(wp) :: B_stencil(3)    !< recons variables stencil for the limiter
    REAL(wp) :: x_stencil(3)    !< grid stencil for the limiter
    REAL(wp) :: y_stencil(3)    !< grid stencil for the limiter

    INTEGER :: limiterB

    INTEGER :: j,k


    ! centered approximation for the topography slope
    limiterB = MAX(limiter(1),1)
    limiterB = 5

    B_cent_extended(2:comp_cells_x+1,2:comp_cells_y+1) = B_cent

    B_cent_extended(1,1) = 3.0_wp * B_cent(1,1) - B_cent(2,1) - B_cent(1,2)

    B_cent_extended(1,comp_cells_y+2) = 3.0_wp * B_cent(1,comp_cells_y)         &
         - B_cent(2,comp_cells_y) - B_cent(1,comp_cells_y-1)

    B_cent_extended(comp_cells_x+2,1) = 3.0_wp * B_cent(comp_cells_x,1)         &
         - B_cent(comp_cells_x,2) - B_cent(comp_cells_x-1,1)

    B_cent_extended(comp_cells_x+2,comp_cells_y+2) =                            &
         3.0_wp * B_cent_extended(comp_cells_x,comp_cells_y)                    &
         - B_cent_extended(comp_cells_x-1,comp_cells_y)                         &
         - B_cent_extended(comp_cells_x,comp_cells_y-1)

    y_loop:DO k = 1,comp_cells_y

       x_loop:DO j = 1,comp_cells_x

          ! x direction
          check_comp_cells_x:IF ( comp_cells_x .GT. 1 ) THEN

             check_x_boundary:IF (j.EQ.1) THEN

                ! west boundary

                x_stencil(1) = 2.0_wp * x_comp(1) - x_comp(2)
                x_stencil(2:3) = x_comp(1:2)

                B_stencil(1) = 2.0_wp * B_cent(1,k) - B_cent(2,k) 
                B_stencil(2:3) = B_cent(1:2,k)

                B_cent_extended(1,k+1) = B_stencil(1)

             ELSEIF (j.EQ.comp_cells_x) THEN

                !east boundary

                x_stencil(3) = 2.0_wp * x_comp(comp_cells_x)                    &
                     - x_comp(comp_cells_x-1)
                x_stencil(1:2) = x_comp(comp_cells_x-1:comp_cells_x)

                B_stencil(3) = 2.0_wp * B_cent(comp_cells_x,k)                  &
                     - B_cent(comp_cells_x-1,k)
                B_stencil(1:2) = B_cent(comp_cells_x-1:comp_cells_x,k)

                B_cent_extended(comp_cells_x+2,k+1) = B_stencil(3)

             ELSE

                ! Internal x interfaces
                x_stencil(1:3) = x_comp(j-1:j+1)
                B_stencil = B_cent(j-1:j+1,k)

             ENDIF check_x_boundary

             B_second_xx(j,k) = ( B_stencil(3) - 2.0_wp * B_stencil(2)          &
                  + B_stencil(1) ) / dx**2  
             CALL limit( B_stencil , x_stencil , limiterB , B_prime_x(j,k) )

          ELSE

             B_prime_x(j,k) = 0.0_wp
             B_second_xx(j,k) = 0.0_wp
             B_cent_extended(1,2:comp_cells_y+1) = B_cent(1,1:comp_cells_y)
             B_cent_extended(comp_cells_x+2,2:comp_cells_y+1) =                 &
                  B_cent(1,1:comp_cells_y)

          END IF check_comp_cells_x

          ! y-direction
          check_comp_cells_y:IF ( comp_cells_y .GT. 1 ) THEN

             check_y_boundary:IF (k.EQ.1) THEN

                ! South boundary
                y_stencil(1) = 2.0_wp * y_comp(1) - y_comp(2)
                y_stencil(2:3) = y_comp(1:2)

                B_stencil(1) = 2.0_wp * B_cent(j,1) - B_cent(j,2)
                B_stencil(2:3) = B_cent(j,1:2)

                B_cent_extended(j+1,1) = B_stencil(1)

             ELSEIF ( k .EQ. comp_cells_y ) THEN

                ! North boundary
                y_stencil(3) = 2.0_wp * y_comp(comp_cells_y)                    &
                     - y_comp(comp_cells_y-1)
                y_stencil(1:2) = y_comp(comp_cells_y-1:comp_cells_y)

                B_stencil(3) = 2.0_wp * B_cent(j,comp_cells_y)                  &
                     - B_cent(j,comp_cells_y-1)
                B_stencil(1:2) = B_cent(j,comp_cells_y-1:comp_cells_y)

                B_cent_extended(j+1,comp_cells_y+2) = B_stencil(3)

             ELSE

                ! Internal y interfaces
                y_stencil(1:3) = y_comp(k-1:k+1)
                B_stencil = B_cent(j,k-1:k+1)

             ENDIF check_y_boundary

             B_second_yy(j,k) = ( B_stencil(3) - 2.0_wp * B_stencil(2)          &
                  + B_stencil(1) ) / dy**2  
             CALL limit( B_stencil , y_stencil , limiterB , B_prime_y(j,k) ) 

          ELSE

             B_prime_y(j,k) = 0.0_wp
             B_second_yy(j,k) = 0.0_wp
             B_cent_extended(2:comp_cells_x+1,1) = B_cent(1:comp_cells_x,1)
             B_cent_extended(2:comp_cells_x+1,comp_cells_y+2) =                 &
                  B_cent(1:comp_cells_x,1)

          ENDIF check_comp_cells_y

       END DO x_loop

    END DO y_loop

    B_second_xy = ( B_cent_extended(3:comp_cells_x+2,3:comp_cells_y+2)          &
         - B_cent_extended(3:comp_cells_x+2,1:comp_cells_y)                     &
         - B_cent_extended(1:comp_cells_x,3:comp_cells_y+2)                     &
         + B_cent_extended(1:comp_cells_x,1:comp_cells_y) ) / ( 4.0_wp*dx*dy )

    IF ( slope_correction_flag ) THEN

       grav_coeff = 1.0_wp / ( 1.0_wp + B_prime_x**2 + B_prime_y**2 )

       d_grav_coeff_dx = - 2.0_wp * grav_coeff**2 * ( B_prime_x * B_second_xx   &
            + B_prime_y * B_second_xy ) 

       d_grav_coeff_dy = - 2.0_wp * grav_coeff**2 * ( B_prime_x * B_second_xy   &
            + B_prime_y * B_second_yy ) 

       grav_coeff_stag_x(1,:) = grav_coeff(1,:)
       grav_coeff_stag_x(2:comp_interfaces_x-1,:) = 0.5_wp *                    &
            ( grav_coeff(1:comp_cells_x-1,:) + grav_coeff(2:comp_cells_x,:) )
       grav_coeff_stag_x(comp_interfaces_x,:) = grav_coeff(comp_cells_x,:)

       grav_coeff_stag_y(:,1) = grav_coeff(1,:)
       grav_coeff_stag_y(:,2:comp_interfaces_y-1) = 0.5_wp *                    &
            ( grav_coeff(:,1:comp_cells_y-1) + grav_coeff(:,2:comp_cells_y) )
       grav_coeff_stag_y(:,comp_interfaces_y) = grav_coeff(:,comp_cells_y)

    ELSE

       grav_coeff = 1.0_wp

       d_grav_coeff_dx = 0.0_wp
       d_grav_coeff_dx = 0.0_wp

       grav_coeff_stag_x = 1.0_wp
       grav_coeff_stag_y = 1.0_wp

    END IF

    RETURN

  END SUBROUTINE topography_reconstruction

  !******************************************************************************
  !> \brief Radial source initialization
  !
  !> In this subroutine the source of mass is initialized. The cells belonging
  !> to the source are are identified ( source_cell(j,k) = 2 ).
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 2021/04/30
  !******************************************************************************

  SUBROUTINE init_source

    USE parameters_2d, ONLY : x_source , y_source , r_source

    IMPLICIT NONE

    INTEGER :: j,k

    REAL(wp) :: total_source

    IF ( verbose_level .GE. 0 ) THEN

       WRITE(*,*) 'r_source',r_source
       WRITE(*,*) 'dx,dy',dx,dy

    END IF

    ! cell where are equations are solved
    source_cell(1:comp_cells_x,1:comp_cells_y) = 0

    sourceE(1:comp_cells_x,1:comp_cells_y) = .FALSE.
    sourceW(1:comp_cells_x,1:comp_cells_y) = .FALSE.
    sourceN(1:comp_cells_x,1:comp_cells_y) = .FALSE.
    sourceS(1:comp_cells_x,1:comp_cells_y) = .FALSE.

    dist_sourceE(1:comp_cells_x,1:comp_cells_y) = 0.0_wp
    dist_sourceW(1:comp_cells_x,1:comp_cells_y) = 0.0_wp
    dist_sourceN(1:comp_cells_x,1:comp_cells_y) = 0.0_wp
    dist_sourceS(1:comp_cells_x,1:comp_cells_y) = 0.0_wp

    total_source = 0.0_wp

    DO k = 2,comp_cells_y-1

       DO j = 2,comp_cells_x-1

          IF ( ( x_comp(j) - x_source )**2 + ( y_comp(k) - y_source )**2 .LE.   &
               r_source**2 ) THEN

             ! cells where equations are not solved
             source_cell(j,k) = 1 

             ! check on west cell
             IF ( ( x_comp(j-1) - x_source )**2 + ( y_comp(k) - y_source )**2   &
                  .GE. r_source**2 ) THEN

                ! cells where radial source boundary condition are applied
                source_cell(j-1,k) = 2
                sourceE(j-1,k) = .TRUE.
                dist_sourceE(j-1,k) = SQRT( ( x_stag(j) - x_source )**2         &
                     + ( y_comp(k) - y_source )**2 )

                sourceE_vect_x(j-1,k) = ( x_stag(j) - x_source ) * r_source     &
                     / dist_sourceE(j-1,k)**2

                sourceE_vect_y(j-1,k) = ( y_comp(k) - y_source ) * r_source     &
                     / dist_sourceE(j-1,k)**2

                total_source = total_source + dx * ABS( sourceE_vect_x(j-1,k) )

             ELSEIF ( ( x_comp(j+1) - x_source )**2 + ( y_comp(k)-y_source )**2 &
                  .GE. r_source**2 ) THEN
                ! check on east cell

                ! cells where radial source boundary condition are applied
                source_cell(j+1,k) = 2
                sourceW(j+1,k) = .TRUE.
                dist_sourceW(j+1,k) = SQRT( ( x_stag(j+1) - x_source )**2       &
                     + ( y_comp(k) - y_source )**2 )

                sourceW_vect_x(j+1,k) = ( x_stag(j+1) - x_source ) * r_source   &
                     / dist_sourceW(j+1,k)**2

                sourceW_vect_y(j+1,k) = ( y_comp(k) - y_source ) * r_source     &
                     / dist_sourceW(j+1,k)**2

                total_source = total_source + dx * ABS( sourceW_vect_x(j+1,k) )

             END IF

             ! check on south cell
             IF ( ( x_comp(j) - x_source )**2 + ( y_comp(k-1) - y_source )**2   &
                  .GE. r_source**2 ) THEN

                ! cells where radial source boundary condition are applied
                source_cell(j,k-1) = 2
                sourceN(j,k-1) = .TRUE.
                dist_sourceN(j,k-1) = SQRT( ( x_comp(j) - x_source )**2         &
                     + ( y_stag(k) - y_source )**2 )

                sourceN_vect_x(j,k-1) = ( x_comp(j) - x_source ) * r_source     &
                     / dist_sourceN(j,k-1)**2

                sourceN_vect_y(j,k-1) = ( y_stag(k) - y_source ) * r_source     &
                     / dist_sourceN(j,k-1)**2

                total_source = total_source + dy * ABS( sourceN_vect_y(j,k-1) )

             ELSEIF ( ( x_comp(j)-x_source )**2 + ( y_comp(k+1) - y_source )**2 &
                  .GE. r_source**2 ) THEN

                ! cells where radial source boundary condition are applied
                source_cell(j,k+1) = 2
                sourceS(j,k+1) = .TRUE.
                dist_sourceS(j,k+1) = SQRT( ( x_comp(j) - x_source )**2         &
                     + ( y_stag(k+1) - y_source )**2 )

                sourceS_vect_x(j,k+1) = ( x_comp(j) - x_source ) * r_source     &
                     / dist_sourceS(j,k+1)**2

                sourceS_vect_y(j,k+1) = ( y_stag(k+1) - y_source ) * r_source   &
                     / dist_sourceS(j,k+1)**2

                total_source = total_source + dy * ABS( sourceS_vect_y(j,k+1) )

             END IF

          END IF

       END DO

    END DO

    RETURN

  END SUBROUTINE init_source

  !-----------------------------------------------------------------------------
  !> Scalar interpolation
  !
  !> This subroutine interpolate the values of the  array f1, defined on the 
  !> grid points x1, at the point x2. The value are saved in f2
  !> \date 13/02/2009
  !> \param[in]    x1           original grid                (\b input)
  !> \param[in]    f1           original values              (\b input)
  !> \param[in]    x2           new point                    (\b output)
  !> \param[out]   f2           interpolated value           (\b output)
  !-----------------------------------------------------------------------------

  SUBROUTINE interp_1d_scalar(x1, f1, x2, f2)
    IMPLICIT NONE

    REAL(wp), INTENT(IN), DIMENSION(:) :: x1, f1
    REAL(wp), INTENT(IN) :: x2
    REAL(wp), INTENT(OUT) :: f2
    INTEGER :: n, n1x, t
    REAL(wp) :: grad , rel_pos

    n1x = SIZE(x1)

    !
    ! ... locate the grid points near the topographic points
    ! ... and interpolate linearly the profile
    !
    t = 1

    search:DO n = 1, n1x-1

       rel_pos = ( x2 - x1(n) ) / ( x1(n+1) - x1(n) )

       IF ( ( rel_pos .GE. 0.0_wp ) .AND. ( rel_pos .LE. 1.0_wp ) ) THEN

          grad = ( f1(n+1)-f1(n) ) / ( x1(n+1)-x1(n) )
          f2 = f1(n) + ( x2-x1(n) ) * grad

          EXIT search

       ELSEIF  ( rel_pos .LT. 0.0_wp ) THEN

          f2 = f1(n)

       ELSE

          f2 = f1(n+1)

       END IF

    END DO search

    RETURN

  END SUBROUTINE interp_1d_scalar

  !-----------------------------------------------------------------------------
  !> Scalar interpolation (2D)
  !
  !> This subroutine interpolate the values of the  array f1, defined on the 
  !> grid points (x1,y1), at the point (x2,y2). The value are saved in f2
  !> \date OCTOBER 2016
  !> \param[in]    x1           original grid       
  !> \param[in]    y1           original grid
  !> \param[in]    f1           original values
  !> \param[in]    x2           new point
  !> \param[in]    y2           new point
  !> \param[out]   f2           interpolated value
  !-----------------------------------------------------------------------------

  SUBROUTINE interp_2d_scalar(x1, y1, f1, x2, y2, f2)
    IMPLICIT NONE

    REAL(wp), INTENT(IN), DIMENSION(:,:) :: x1, y1, f1
    REAL(wp), INTENT(IN) :: x2, y2
    REAL(wp), INTENT(OUT) :: f2

    INTEGER :: ix , iy
    REAL(wp) :: alfa_x , alfa_y

    IF ( size(x1,1) .GT. 1 ) THEN

       ix = FLOOR( ( x2 - x1(1,1) ) / ( x1(2,1) - x1(1,1) ) ) + 1
       ix = MIN( ix , SIZE(x1,1)-1 )
       alfa_x = ( x1(ix+1,1) - x2 ) / (  x1(ix+1,1) - x1(ix,1) )

    ELSE

       ix = 1
       alfa_x = 0.0_wp

    END IF

    IF ( size(x1,2) .GT. 1 ) THEN

       iy = FLOOR( ( y2 - y1(1,1) ) / ( y1(1,2) - y1(1,1) ) ) + 1
       iy = MIN( iy , SIZE(x1,2)-1 )
       alfa_y = ( y1(1,iy+1) - y2 ) / (  y1(1,iy+1) - y1(1,iy) )

    ELSE

       iy = 1
       alfa_y = 0.0_wp

    END IF
    
    IF ( size(x1,1) .EQ. 1 ) THEN

       f2 = alfa_y * f1(ix,iy) + ( 1.0_wp - alfa_y ) * f1(ix,iy+1)

    ELSEIF ( size(x1,2) .EQ. 1 ) THEN

       f2 = alfa_x * f1(ix,iy)  + ( 1.0_wp - alfa_x ) * f1(ix+1,iy)

    ELSE

       f2 = alfa_x * ( alfa_y * f1(ix,iy) + ( 1.0_wp - alfa_y ) * f1(ix,iy+1) ) &
            + ( 1.0_wp - alfa_x ) * ( alfa_y * f1(ix+1,iy) + ( 1.0_wp - alfa_y )&
            * f1(ix+1,iy+1) )

    END IF

  END SUBROUTINE interp_2d_scalar

  !-----------------------------------------------------------------------------
  !> Scalar interpolation (2D)
  !
  !> This subroutine interpolate the values of the  array f1, defined on the 
  !> grid points (x1,y1), at the point (x2,y2). The value are saved in f2
  !> \date OCTOBER 2016
  !> \param[in]    x1           original grid          
  !> \param[in]    y1           original grid          
  !> \param[in]    f1           original values      
  !> \param[in]    x2           new point            
  !> \param[in]    y2           new point            
  !> \param[out]   f2           interpolated value     
  !-----------------------------------------------------------------------------

  SUBROUTINE interp_2d_nodata(x1, y1, f1, x2, y2, f2)
    IMPLICIT NONE

    REAL(wp), INTENT(IN), DIMENSION(:,:) :: x1, y1, f1
    REAL(wp), INTENT(IN) :: x2, y2
    LOGICAL, INTENT(OUT) :: f2

    INTEGER :: ix , iy
    REAL(wp) :: alfa_x , alfa_y

    IF ( size(x1,1) .GT. 1 ) THEN

       ix = FLOOR( ( x2 - x1(1,1) ) / ( x1(2,1) - x1(1,1) ) ) + 1
       ix = MIN( ix , SIZE(x1,1)-1 )
       alfa_x = ( x1(ix+1,1) - x2 ) / (  x1(ix+1,1) - x1(ix,1) )

    ELSE

       ix = 1
       alfa_x = 0.0_wp

    END IF

    IF ( size(x1,2) .GT. 1 ) THEN

       iy = FLOOR( ( y2 - y1(1,1) ) / ( y1(1,2) - y1(1,1) ) ) + 1
       iy = MIN( iy , SIZE(x1,2)-1 )
       alfa_y = ( y1(1,iy+1) - y2 ) / (  y1(1,iy+1) - y1(1,iy) )

    ELSE

       iy = 1
       alfa_y = 0.0_wp

    END IF

    f2 = .FALSE.

    IF ( size(x1,1) .EQ. 1 ) THEN

       f2 = ( f1(ix,iy) .EQ. nodata_topo ) .OR. ( f1(ix,iy+1) .EQ. nodata_topo )

    ELSEIF ( size(x1,2) .EQ. 1 ) THEN

       f2 = ( f1(ix,iy) .EQ. nodata_topo ) .OR. ( f1(ix+1,iy) .EQ. nodata_topo ) 

    ELSE

       f2 = ( f1(ix,iy) .EQ. nodata_topo ) .OR. ( f1(ix,iy+1) .EQ. nodata_topo )&
            .OR. ( f1(ix+1,iy) .EQ. nodata_topo )                               &
            .OR. ( f1(ix+1,iy+1) .EQ. nodata_topo )

    END IF

  END SUBROUTINE interp_2d_nodata


  !-----------------------------------------------------------------------------
  !> Scalar interpolation (2D)
  !
  !> This subroutine interpolate the values of the  array f1, defined on the 
  !> grid points (x1,y1), at the point (x2,y2). The value are saved in f2.
  !> In this case x1 and y1 are 1d arrays.
  !> \date OCTOBER 2016
  !> \param[in]    x1           original grid              
  !> \param[in]    y1           original grid            
  !> \param[in]    f1           original values           
  !> \param[in]    x2           new point                
  !> \param[in]    y2           new point               
  !> \param[out]   f2           interpolated value       
  !-----------------------------------------------------------------------------

  SUBROUTINE interp_2d_scalarB(x1, y1, f1, x2, y2, f2)
    IMPLICIT NONE

    REAL(wp), INTENT(IN), DIMENSION(:) :: x1, y1
    REAL(wp), INTENT(IN), DIMENSION(:,:) :: f1
    REAL(wp), INTENT(IN) :: x2, y2
    REAL(wp), INTENT(OUT) :: f2

    INTEGER :: ix , iy
    REAL(wp) :: alfa_x , alfa_y

    IF ( size(x1) .GT. 1 ) THEN

       ix = FLOOR( ( x2 - x1(1) ) / ( x1(2) - x1(1) ) ) + 1
       ix = MAX(0,MIN( ix , SIZE(x1)-1 ))
       alfa_x = ( x1(ix+1) - x2 ) / (  x1(ix+1) - x1(ix) )

    ELSE

       ix = 1
       alfa_x = 0.0_wp

    END IF

    IF ( size(y1) .GT. 1 ) THEN

       iy = FLOOR( ( y2 - y1(1) ) / ( y1(2) - y1(1) ) ) + 1
       iy = MAX(1,MIN( iy , SIZE(y1)-1 ))
       alfa_y = ( y1(iy+1) - y2 ) / (  y1(iy+1) - y1(iy) )

    ELSE

       iy = 1
       alfa_y = 0.0_wp

    END IF

    IF ( ( alfa_x .LT. 0.0_wp ) .OR. ( alfa_x .GT. 1.0_wp )                     &
         .OR. ( alfa_y .LT. 0.0_wp ) .OR. ( alfa_y .GT. 1.0_wp ) ) THEN

       f2 = 0.0_wp
       RETURN

    END IF


    IF ( size(x1) .EQ. 1 ) THEN

       f2 = alfa_y * f1(ix,iy) + ( 1.0_wp - alfa_y ) * f1(ix,iy+1)

    ELSEIF ( size(y1) .EQ. 1 ) THEN

       f2 = alfa_x * f1(ix,iy)  + ( 1.0_wp - alfa_x ) * f1(ix+1,iy)

    ELSE

       f2 = alfa_x * ( alfa_y * f1(ix,iy) + ( 1.0_wp - alfa_y ) * f1(ix,iy+1) ) &
            + ( 1.0_wp - alfa_x ) * ( alfa_y * f1(ix+1,iy) + ( 1.0_wp - alfa_y )&
            * f1(ix+1,iy+1) )

    END IF

    RETURN

  END SUBROUTINE interp_2d_scalarB


  !-----------------------------------------------------------------------------
  !> Scalar interpolation (2D)
  !
  !> This subroutine interpolate the values of the  array f1, defined on the 
  !> grid points (x1,y1), at the point (x2,y2). The value are saved in f2
  !> \date OCTOBER 2016
  !> \param[in]    x1           original grid               
  !> \param[in]    y1           original grid                
  !> \param[in]    f1           original values             
  !> \param[in]    x2           new point                   
  !> \param[in]    y2           new point                  
  !> \param[out]   f2           interpolated value         
  !-----------------------------------------------------------------------------

  SUBROUTINE interp_2d_slope(x1, y1, f1, x2, y2, f_x, f_y)

    IMPLICIT NONE

    REAL(wp), INTENT(IN), DIMENSION(:,:) :: x1, y1, f1
    REAL(wp), INTENT(IN) :: x2, y2
    REAL(wp), INTENT(OUT) :: f_x , f_y

    INTEGER :: ix , iy
    REAL(wp) :: alfa_x , alfa_y
    REAL(wp) :: f_x1 , f_x2 , f_y1 , f_y2


    IF ( size(x1,1) .GT. 1 ) THEN

       ix = FLOOR( ( x2 - x1(1,1) ) / ( x1(2,1) - x1(1,1) ) ) + 1
       ix = MIN( ix , SIZE(x1,1)-1 )
       alfa_x = ( x1(ix+1,1) - x2 ) / (  x1(ix+1,1) - x1(ix,1) )

    ELSE

       ix = 1
       alfa_x = 1.0_wp

    END IF

    IF ( size(x1,2) .GT. 1 ) THEN

       iy = FLOOR( ( y2 - y1(1,1) ) / ( y1(1,2) - y1(1,1) ) ) + 1
       iy = MIN( iy , SIZE(x1,2)-1 )
       alfa_y = ( y1(1,iy+1) - y2 ) / (  y1(1,iy+1) - y1(1,iy) )

    ELSE

       iy = 1
       alfa_y = 1.0_wp

    END IF

    f_x1 = 0.0_wp
    f_x2 = 0.0_wp

    IF ( size(x1,1) .GT. 1 ) THEN

       f_x1 = ( f1(ix+1,iy) -  f1(ix,iy) ) / (  x1(2,1) - x1(1,1) )

       IF ( size(x1,2) .GT. 1 ) THEN

          f_x2 = ( f1(ix+1,iy+1) -  f1(ix,iy+1) ) / (  x1(2,1) - x1(1,1) )

       END IF

    END IF

    f_x = alfa_y * f_x1 + ( 1.0_wp - alfa_y ) * f_x2

    f_y1 = 0.0_wp
    f_y2 = 0.0_wp

    IF ( size(x1,2) .GT. 1 ) THEN

       f_y1 = ( f1(ix,iy+1) - f1(ix,iy) ) / (  y1(1,2) - y1(1,1) )

       IF ( size(x1,1) .GT. 1 ) THEN

          f_y2 = ( f1(ix+1,iy+1) - f1(ix+1,iy) ) / (  y1(1,2) - y1(1,1) )

       END IF

    END IF

    f_y = alfa_x * f_y1 + ( 1.0_wp - alfa_x ) * f_y2


    RETURN

  END SUBROUTINE interp_2d_slope


  !-----------------------------------------------------------------------------
  !> Scalar regrid (2D)
  !
  !> This subroutine interpolate the values of the  array f1, defined on the 
  !> grid points (x1,y1), at the point (x2,y2). The value are saved in f2.
  !> In this case x1 and y1 are 1d arrays.
  !> \date OCTOBER 2016
  !> \param[in]    x1           original grid
  !> \param[in]    y1           original grid
  !> \param[in]    f1           original values
  !> \param[in]    xl           new point
  !> \param[in]    xr           new point
  !> \param[in]    yl           new point
  !> \param[in]    yr           new point
  !> \param[out]   f2           interpolated value
  !-----------------------------------------------------------------------------

  SUBROUTINE regrid_scalar(xin, yin, fin, xl, xr , yl, yr, fout)
    IMPLICIT NONE

    REAL(wp), INTENT(IN), DIMENSION(:) :: xin, yin
    REAL(wp), INTENT(IN), DIMENSION(:,:) :: fin
    REAL(wp), INTENT(IN) :: xl, xr , yl , yr
    REAL(wp), INTENT(OUT) :: fout

    INTEGER :: ix , iy
    INTEGER :: ix1 , ix2 , iy1 , iy2
    REAL(wp) :: alfa_x , alfa_y
    REAL(wp) :: dXin , dYin

    INTEGER nXin,nYin

    nXin = size(xin)-1
    nYin = size(yin)

    dXin = xin(2) - xin(1)
    dYin = yin(2) - yin(1)

    ix1 = MAX(1,CEILING( ( xl - xin(1) ) / dXin ))
    ix2 = MIN(nXin,CEILING( ( xr -xin(1) ) / dXin )+1)

    iy1 = MAX(1,CEILING( ( yl - yin(1) ) / dYin ))
    iy2 = MIN(nYin,CEILING( ( yr - yin(1) ) / dYin ) + 1)

    fout = 0.0_wp

    DO ix=ix1,ix2-1

       alfa_x = ( MIN(xr,xin(ix+1)) - MAX(xl,xin(ix)) ) / ( xr - xl )

       DO iy=iy1,iy2-1

          alfa_y = ( MIN(yr,yin(iy+1)) - MAX(yl,yin(iy)) ) / ( yr - yl )

          fout = fout + alfa_x * alfa_y * fin(ix,iy)

       END DO

    END DO

  END SUBROUTINE regrid_scalar

  !******************************************************************************
  !> \brief Slope limiter
  !
  !> This subroutine limits the slope of the linear reconstruction of 
  !> the physical variables, accordingly to the parameter "solve_limiter":\n
  !> - 'none'     => no limiter (constant value);
  !> - 'minmod'   => minmod slope;
  !> - 'superbee' => superbee limiter (Roe, 1985);
  !> - 'van_leer' => monotonized central-difference limiter (van Leer, 1977)
  !> .
  !> \param[in]     v             3-point stencil value array 
  !> \param[in]     z             3-point stencil location array 
  !> \param[in]     limiter       integer defining the limiter choice
  !> \param[out]    slope_lim     limited slope         
  !> \date 07/10/2016
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  SUBROUTINE limit( v , z , limiter , slope_lim )

    USE parameters_2d, ONLY : theta

    IMPLICIT none

    REAL(wp), INTENT(IN) :: v(3)
    REAL(wp), INTENT(IN) :: z(3)
    INTEGER, INTENT(IN) :: limiter

    REAL(wp), INTENT(OUT) :: slope_lim

    REAL(wp) :: a , b , c

    REAL(wp) :: sigma1 , sigma2

    a = ( v(3) - v(2) ) / ( z(3) - z(2) )
    b = ( v(2) - v(1) ) / ( z(2) - z(1) )
    c = ( v(3) - v(1) ) / ( z(3) - z(1) )

    SELECT CASE (limiter)

    CASE ( 0 )

       slope_lim = 0.0_wp

    CASE ( 1 )

       ! minmod
       slope_lim = minmod(a,b)

    CASE ( 2 )

       ! superbee
       sigma1 = minmod( a , 2.0_wp*b )
       sigma2 = minmod( 2.0_wp*a , b )
       slope_lim = maxmod( sigma1 , sigma2 )

    CASE ( 3 )

       ! generalized minmod
       slope_lim = minmod( c , theta * minmod( a , b ) )

    CASE ( 4 )

       ! monotonized central-difference (MC, LeVeque p.112)
       slope_lim = minmod( c , 2.0 * minmod( a , b ) )

    CASE ( 5 )

       ! centered
       slope_lim = c

    CASE (6) 

       ! backward
       slope_lim = a

    CASE (7)

       !forward
       slope_lim = b

    END SELECT

  END SUBROUTINE limit

  !******************************************************************************
  !> \brief Minmod limiter
  !
  !> This function compute the minmod between two real numbers
  !> \param[in]     a: 1st value
  !> \param[in]     b: 2nd value
  !> \result        the minmod between a and b
  !> \date 2021/04/30
  !> @author 
  !> Mattia de' Michieli Vitturi
  !******************************************************************************

  REAL(wp) FUNCTION minmod(a,b)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: a
    REAL(wp), INTENT(IN) :: b
    REAL(wp) :: sa , sb 

    IF ( MIN(ABS(a),ABS(b)) .LE. 0.0e-40_wp ) THEN

       minmod = 0.0_wp

    ELSE

       sa = a / ABS(a)
       sb = b / ABS(b)

       minmod = 0.5_wp * ( sa+sb ) * MIN( ABS(a) , ABS(b) )

    END IF

  END FUNCTION minmod

  REAL(wp) function maxmod(a,b)

    IMPLICIT none

    REAL(wp) :: a , b , sa , sb 

    IF ( ABS(a*b) .LE. 0.0e-30_wp ) THEN

       maxmod = 0.0_wp

    ELSE

       sa = a / ABS(a)
       sb = b / ABS(b)

       maxmod = 0.5_wp * ( sa+sb ) * MAX( ABS(a) , ABS(b) )

    END IF

  END function maxmod

  SUBROUTINE compute_cell_fract(xs,ys,rs,cell_fract)

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: xs,ys,rs

    REAL(wp), INTENT(OUT) :: cell_fract(comp_cells_x,comp_cells_y)

    REAL(wp), ALLOCATABLE :: x_subgrid(:) , y_subgrid(:)

    INTEGER, ALLOCATABLE :: check_subgrid(:)

    INTEGER n_points , n_points2

    REAL(wp) :: source_area

    INTEGER :: h ,j, k

    n_points = 200
    n_points2 = n_points**2

    ALLOCATE( x_subgrid(n_points2) )
    ALLOCATE( y_subgrid(n_points2) )
    ALLOCATE( check_subgrid(n_points2) )

    x_subgrid = 0.0_wp
    y_subgrid = 0.0_wp

    DO h = 1,n_points

       x_subgrid(h:n_points2:n_points) = DBLE(h)
       y_subgrid((h-1)*n_points+1:h*n_points) = DBLE(h)

    END DO

    x_subgrid = ( 2.0_wp * x_subgrid - 1.0_wp ) / ( 2.0_wp * DBLE(n_points) )
    y_subgrid = ( 2.0_wp * y_subgrid - 1.0_wp ) / ( 2.0_wp * DBLE(n_points) )

    x_subgrid = ( x_subgrid - 0.5_wp ) * dx
    y_subgrid = ( y_subgrid - 0.5_wp ) * dy

    DO j=1,comp_cells_x

       DO k=1,comp_cells_y

          IF ( ( x_stag(j+1) .LT. ( xs - rs ) ) .OR.                             &
               ( x_stag(j) .GT. ( xs + rs ) ) .OR.                               &
               ( y_stag(k+1) .LT. ( ys - rs ) ) .OR.                             &
               ( y_stag(k) .GT. ( ys + rs ) ) ) THEN

             cell_fract(j,k) = 0.0_wp 

          ELSE

             check_subgrid = 0

             WHERE ( ( x_comp(j) + x_subgrid - xs )**2                           &
                  + ( y_comp(k) + y_subgrid - ys )**2 < rs**2 )

                check_subgrid = 1

             END WHERE

             cell_fract(j,k) = REAL(SUM(check_subgrid))/n_points2

          END IF

       ENDDO

    ENDDO

    source_area = dx*dy*SUM(cell_fract)

    IF ( VERBOSE_LEVEL .GE. 0 ) THEN

       WRITE(*,*) 'Source area =',source_area,' Error =',ABS( 1.0_wp -          &
            dx*dy*SUM(cell_fract) / ( 4.0_wp*ATAN(1.0_wp)*rs**2 ) )

    END IF

    DEALLOCATE( x_subgrid )
    DEALLOCATE( y_subgrid )
    DEALLOCATE( check_subgrid )

    RETURN

  END SUBROUTINE compute_cell_fract


  real*8 function lambertws(p)
    !
    !   Compute W(z) by a truncated Maclaurin series around the branch point W=-1
    !
    !   Reference: Toshio Fukushima (2013) J.Comp.Appl.Math., 244, 77-89
    !               "Precise and fast computation of Lambert W-functions
    !                without transcendental function evaluations"
    !
    real*8 p,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15,P16,P17,P18,P19,P20
    parameter (P2=-0.33333333333333333333d0)
    parameter (P3=0.15277777777777777778d0)
    parameter (P4=-0.079629629629629629630d0)
    parameter (P5=0.044502314814814814815d0)
    parameter (P6=-0.025984714873603762493d0)
    parameter (P7=0.015635632532333921223d0)
    parameter (P8=-0.0096168920242994317068d0)
    parameter (P9=0.0060145432529561178610d0)
    parameter (P10=-0.0038112980348919992267d0)
    parameter (P11=0.0024408779911439826659d0)
    parameter (P12=-0.0015769303446867842539d0)
    parameter (P13=0.0010262633205076071544d0)
    parameter (P14=-0.00067206163115613620400d0)
    parameter (P15=0.00044247306181462090993d0)
    parameter (P16=-0.00029267722472962744485d0)
    parameter (P17=0.00019438727605453931782d0)
    parameter (P18=-0.00012957426685274881888d0)
    parameter (P19=0.000086650358052081271660d0)
    parameter (P20=-0.000058113607504413816772d0)
    !
    if(abs(p).lt.0.01159d0) then
       lambertws=-1.d0+p*(1.d0+p*(P2+p*(P3+p*(P4+p*(P5+p*P6)))))
    elseif(abs(p).lt.0.0766d0) then
       lambertws=-1.d0+p*(1.d0+p*(P2+p*(P3+p*(P4+p*(P5+p*(P6+p*(P7+p*(P8+p*(P9+p*P10)))))))))
    else
       lambertws=-1.d0+p*(1.d0+p*(P2+p*(P3+p*(P4+p*(P5+p*(P6+p*(P7+p*(P8+p*(P9+p*(P10+p*(P11 &
            +p*(P12+p*(P13+p*(P14+p*(P15+p*(P16+p*(P17+p*(P18+p*(P19+p*P20)))))))))))))))))))
    endif
    return
  end function lambertws


  real*8 function lambertw(z)

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: z
    REAL*8 :: w1,w2
    REAL*8 :: z1,z2

    w1 = lambertwm1(z)
    z1 = w1*exp(w1)

    w2 = lambertw0(z)
    z2 = w2*exp(w2)

    IF ( ABS( z1-z ) .LT. ABS( z2-z ) ) THEN

       lambertw = w1

    ELSE

       lambertw = w2

    END IF

    RETURN

  END function lambertw


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real*8 function lambertw0(z)
    !
    !   Compute W_0(z), Lambert W-function of the branch 0
    !
    !   calls "lambertws"
    !
    !   Reference: Toshio Fukushima (2013) J.Comp.Appl.Math., 244, 77-89
    !               "Precise and fast computation of Lambert W-functions
    !                without transcendental function evaluations"
    !
    integer n1,n,j,jmax,nh
    real*8 z,p2,y,w,wj,yj,f0,f1,f00,f11,f0y
    real*8 Em1,E1,ej,em(-1:64),g(0:64),a(12),b(12)
    logical first/.TRUE./
    parameter (E1=2.718281828459045235d0,Em1=1.d0/E1)
    save first,em,g,a,b
    if(first) then
       first=.FALSE.;em(-1)=E1;ej=1.d0;em(0)=1.d0;g(0)=0.d0
       do j=1,64
          ej=ej*E1;em(j)=em(j-1)*Em1;g(j)=dble(j)*ej
       enddo
       a(1)=sqrt(Em1);b(1)=0.5d0
       do j=2,12
          a(j)=sqrt(a(j-1));b(j)=b(j-1)*0.5d0
       enddo
    endif
    if(abs(z).lt.0.05d0) then
       lambertw0=z*(1.d0-z*(1.0-z*(1.5d0-z*(2.6666666666666666667d0 &
            -z*(5.2083333333333333333d0-z*(10.8d0-z*(23.343055555555555556d0 &
            -z*(52.012698412698412698d0-z*(118.62522321428571429d0 &
            -z*(275.57319223985890653d0-z*(649.78717234347442681d0 &
            -z*(1551.1605194805194805d0-z*(3741.4497029592385495d0 &
            -z*(9104.5002411580189358d0-z*(22324.308512706601434d0 &
            -z*(55103.621972903835338d0-z*136808.86090394293563d0 &
            ))))))))))))))))
       return
    elseif(z.lt.-0.35d0) then
       p2=2.d0*(E1*z+1.d0)
       if(p2.GT.0.d0) then
          lambertw0=lambertws(sqrt(p2))
       elseif(p2.EQ.0.d0) then
          lambertw0=-1.d0
       else
          ! write(*,*) "(lambertw0) Argument out of range. z=",z
          lambertw0 = -1.0
       endif
       return
    endif
    do n=0,2
       if(g(n).GT.z) goto 1
    enddo
    n=2
    do j=1,5
       n=n*2
       if(g(n).GT.z) goto 2
    enddo
    ! write(*,*) "(lambertw0) Too large argument. z=",z
    lambertw0 = -1.0
    return
2   continue
    nh=n/2
    do j=1,5
       nh=nh/2
       if(nh.LE.0) exit
       if(g(n-nh).GT.z) then
          n=n-nh
       endif
    enddo
1   continue
    n=n-1;y=z*em(n)
    w=dble(n)
    if(z.le.-0.36d0) then
       jmax=12
    elseif(z.le.-0.3d0) then
       jmax=11
    elseif(n.le.0) then
       jmax=10
    elseif(n.le.1) then
       jmax=9
    else
       jmax=8
    endif
    do j=1,jmax
       wj=w+b(j); yj=y*a(j)
       if(wj.lt.yj) then
          w=wj;y=yj
       endif
    enddo
    f0=w-y;f1=1.d0+y;f00=f0*f0;f11=f1*f1;f0y=f0*y
    lambertw0=w-4.d0*f0*(6.d0*f1*(f11+f0y)+f00*y)/ &
         (f11*(24.d0*f11+36.d0*f0y)+f00*(6.d0*y*y+8.d0*f1*y+f0y))
    return
  end function lambertw0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  real*8 function lambertwm1(z)
    !
    !   Compute W_{-1}(z), Lambert W-function of the branch -1
    !
    !   calls "lambertws"
    !
    !   Reference: Toshio Fukushima (2013) J.Comp.Appl.Math., 244, 77-89
    !               "Precise and fast computation of Lambert W-functions
    !                without transcendental function evaluations"
    !
    integer n,j,jmax,nh
    real*8 z,p2,y,w,wj,yj,f0,f1,f00,f11,f0y
    real*8 E1,Em1,emj,e(64),g(64),a(12),b(12)
    logical first/.TRUE./
    parameter (E1=2.718281828459045235d0,Em1=1.d0/E1)
    save first,e,g,a,b
    if(first) then
       first=.FALSE.;emj=Em1;e(1)=E1;g(1)=-Em1
       do j=2,64
          emj=emj*Em1;e(j)=e(j-1)*E1;g(j)=-dble(j)*emj
       enddo
       a(1)=sqrt(E1);b(1)=0.5d0
       do j=2,12
          a(j)=sqrt(a(j-1));b(j)=b(j-1)*0.5d0
       enddo
    endif
    if(z.GE.0.d0) then
       ! write(*,*) "(lambertwm1) Argument out of range. z=",z
       lambertwm1 = -1.0
       return
    elseif(z.LT.-0.35d0) then
       p2=2.d0*(E1*z+1.d0)
       if(p2.GT.0.d0) then
          lambertwm1=lambertws(-sqrt(p2))
       elseif(p2.EQ.0.d0) then
          lambertwm1=-1.d0
       else
          ! write(*,*) "(lambertwm1) Argument out of range. z=",z
          lambertwm1 = -1.0
       endif
       return
    endif
    n=2
    if(g(n).GT.z) goto 1
    do j=1,5
       n=n*2
       if(g(n).GT.z) goto 2
    enddo
    ! write(*,*) "(lambertwm1) Too small argument. z=",z
    lambertwm1 = -1.0
    return
2   continue
    nh=n/2
    do j=1,5
       nh=nh/2
       if(nh.LE.0) exit
       if(g(n-nh).GT.z) then
          n=n-nh
       endif
    enddo
1   continue
    n=n-1;w=-dble(n);y=z*e(n)
    if(n.GE.8) then
       jmax=8
    elseif(n.GE.3) then
       jmax=9
    elseif(n.GE.2) then
       jmax=10
    else
       jmax=11
    endif
    do j=1,jmax
       wj=w-b(j);yj=y*a(j)
       if(wj.lt.yj) then
          w=wj;y=yj
       endif
    enddo
    f0=w-y;f1=1.d0+y;f00=f0*f0;f11=f1*f1;f0y=f0*y
    lambertwm1=w-4.d0*f0*(6.d0*f1*(f11+f0y)+f00*y)/ &
         (f11*(24.d0*f11+36.d0*f0y)+f00*(6.d0*y*y+8.d0*f1*y+f0y))
    return
  end function lambertwm1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine calcei ( arg, result, int )

    !*****************************************************************************80
    !
    !! CALCEI computes exponential integrals.
    !
    !  Discussion:
    !
    !    This routine computes the exponential integrals Ei(x),
    !    E1(x), and  exp(-x)*Ei(x)  for real arguments  x  
    !    where, if x > 0,
    !      Ei(x) = integral (from t=-infinity to t=x) (exp(t)/t),  x > 0,
    !    while, if x < 0,
    !      Ei(x) = -integral (from t=-x to t=infinity) (exp(t)/t),  x < 0,
    !    and where the first integral is a principal value integral.
    !
    !    The packet contains three function type subprograms: EI, EONE,
    !    and EXPEI;  and one subroutine type subprogram: CALCEI.  The
    !    calling statements for the primary entries are:
    !      Y = EI(X),    where  X /= 0,
    !      Y = EONE(X),      where  X > 0,
    !      Y = EXPEI(X),     where  X /= 0,
    !
    !    and where the entry points correspond to the functions Ei(x),
    !    E1(x), and exp(-x)*Ei(x), respectively.  The routine CALCEI
    !    is intended for internal packet use only, all computations within
    !    the packet being concentrated in this routine.  The function
    !    subprograms invoke CALCEI with the statement
    !      CALL CALCEI(ARG,RESULT,INT)
    !    where the parameter usage is as follows
    !
    !      Function      Parameters for CALCEI
    !      Call           ARG         RESULT         INT
    !
    !      EI(X)          X /= 0      Ei(X)            1
    !      EONE(X)        X > 0      -Ei(-X)           2
    !      EXPEI(X)       X /= 0      exp(-X)*Ei(X)    3
    !
    !    The main computation involves evaluation of rational Chebyshev
    !    approximations published in Math. Comp. 22, 641-649 (1968), and
    !    Math. Comp. 23, 289-303 (1969) by Cody and Thacher.  This
    !    transportable program is patterned after the machine-dependent
    !    FUNPACK packet  NATSEI,  but cannot match that version for
    !    efficiency or accuracy.  This version uses rational functions
    !    that theoretically approximate the exponential integrals to
    !    at least 18 significant decimal digits.  The accuracy achieved
    !    depends on the arithmetic system, the compiler, the intrinsic
    !    functions, and proper selection of the machine-dependent
    !    constants.
    !
    !  Modified:
    !
    !    10 January 2016
    !
    !  Author:
    !
    !    Original FORTRAN77 version by William Cody.
    !    FORTRAN90 version by John Burkardt.
    !
    ! Explanation of machine-dependent constants.  Let
    !
    !   beta = radix for the floating-point system.
    !   minexp = smallest representable power of beta.
    !   maxexp = smallest power of beta that overflows.
    !
    ! Then the following machine-dependent constants must be declared
    !   in DATA statements.  IEEE values are provided as a default.
    !
    !   XBIG = largest argument acceptable to EONE; solution to
    !      equation:
    !         exp(-x)/x * (1 + 1/x) = beta ** minexp.
    !   XINF = largest positive machine number; approximately
    !         beta ** maxexp
    !   XMAX = largest argument acceptable to EI; solution to
    !      equation:  exp(x)/x * (1 + 1/x) = beta ** maxexp.
    !
    ! Error returns
    !
    !  The following table shows the types of error that may be
    !  encountered in this routine and the function value supplied
    !  in each case.
    !
    !   Error   Argument   function values for
    !        Range     EI  EXPEI     EONE
    !
    !     UNDERFLOW  (-)X > XBIG     0    -     0
    !     OVERFLOW  X >= XMAX    XINF  -     -
    !     ILLEGAL X   X = 0   -XINF    -XINF     XINF
    !     ILLEGAL X  X < 0   -    -     USE ABS(X)
    !
    implicit none

    integer ( kind = 4 ) i,int
    real(wp) &
         a,arg,b,c,d,exp40,e,ei,f,four,fourty,frac,half,one,p, &
         plg,px,p037,p1,p2,q,qlg,qx,q1,q2,r,result,s,six,sump, &
         sumq,t,three,twelve,two,two4,w,x,xbig,xinf,xmax,xmx0, &
         x0,x01,x02,x11,y,ysq,zero
    dimension  a(7),b(6),c(9),d(9),e(10),f(10),p(10),q(10),r(10), &
         s(9),p1(10),q1(9),p2(10),q2(9),plg(4),qlg(4),px(10),qx(10)
    !
    !  Mathematical constants
    !  EXP40 = exp(40)
    !  X0 = zero of Ei
    !  X01/X11 + X02 = zero of Ei to extra precision
    !
    data zero,p037,half,one,two/0.0d0,0.037d0,0.5d0,1.0d0,2.0d0/, &
         three,four,six,twelve,two4/3.0d0,4.0d0,6.0d0,12.d0,24.0d0/, &
         fourty,exp40/40.0d0,2.3538526683701998541d17/, &
         x01,x11,x02/381.5d0,1024.0d0,-5.1182968633365538008d-5/, &
         x0/3.7250741078136663466d-1/
    !
    !  Machine-dependent constants
    !
    data xinf/1.79d+308/,xmax/716.351d0/,xbig/701.84d0/
    !
    ! Coefficients  for -1.0 <= X < 0.0
    !
    data a/1.1669552669734461083368d2, 2.1500672908092918123209d3, &
         1.5924175980637303639884d4, 8.9904972007457256553251d4, &
         1.5026059476436982420737d5,-1.4815102102575750838086d5, &
         5.0196785185439843791020d0/
    data b/4.0205465640027706061433d1, 7.5043163907103936624165d2, &
         8.1258035174768735759855d3, 5.2440529172056355429883d4, &
         1.8434070063353677359298d5, 2.5666493484897117319268d5/
    !
    ! Coefficients for -4.0 <= X < -1.0
    !
    data c/3.828573121022477169108d-1, 1.107326627786831743809d+1, &
         7.246689782858597021199d+1, 1.700632978311516129328d+2, &
         1.698106763764238382705d+2, 7.633628843705946890896d+1, &
         1.487967702840464066613d+1, 9.999989642347613068437d-1, &
         1.737331760720576030932d-8/
    data d/8.258160008564488034698d-2, 4.344836335509282083360d+0, &
         4.662179610356861756812d+1, 1.775728186717289799677d+2, &
         2.953136335677908517423d+2, 2.342573504717625153053d+2, &
         9.021658450529372642314d+1, 1.587964570758947927903d+1, &
         1.000000000000000000000d+0/
    !
    ! Coefficients for X < -4.0
    !
    data e/1.3276881505637444622987d+2,3.5846198743996904308695d+4, &
         1.7283375773777593926828d+5,2.6181454937205639647381d+5, &
         1.7503273087497081314708d+5,5.9346841538837119172356d+4, &
         1.0816852399095915622498d+4,1.0611777263550331766871d03, &
         5.2199632588522572481039d+1,9.9999999999999999087819d-1/
    data f/3.9147856245556345627078d+4,2.5989762083608489777411d+5, &
         5.5903756210022864003380d+5,5.4616842050691155735758d+5, &
         2.7858134710520842139357d+5,7.9231787945279043698718d+4, &
         1.2842808586627297365998d+4,1.1635769915320848035459d+3, &
         5.4199632588522559414924d+1,1.0d0/
    !
    !  Coefficients for rational approximation to ln(x/a), |1-x/a| < .1
    !
    data plg/-2.4562334077563243311d+01,2.3642701335621505212d+02, &
         -5.4989956895857911039d+02,3.5687548468071500413d+02/
    data qlg/-3.5553900764052419184d+01,1.9400230218539473193d+02, &
         -3.3442903192607538956d+02,1.7843774234035750207d+02/
    !
    ! Coefficients for  0.0 < X < 6.0,
    !  ratio of Chebyshev polynomials
    !
    data p/-1.2963702602474830028590d01,-1.2831220659262000678155d03, &
         -1.4287072500197005777376d04,-1.4299841572091610380064d06, &
         -3.1398660864247265862050d05,-3.5377809694431133484800d08, &
         3.1984354235237738511048d08,-2.5301823984599019348858d10, &
         1.2177698136199594677580d10,-2.0829040666802497120940d11/
    data q/ 7.6886718750000000000000d01,-5.5648470543369082846819d03, &
         1.9418469440759880361415d05,-4.2648434812177161405483d06, &
         6.4698830956576428587653d07,-7.0108568774215954065376d08, &
         5.4229617984472955011862d09,-2.8986272696554495342658d10, &
         9.8900934262481749439886d10,-8.9673749185755048616855d10/
    !
    ! J-fraction coefficients for 6.0 <= X < 12.0
    !
    data r/-2.645677793077147237806d00,-2.378372882815725244124d00, &
         -2.421106956980653511550d01, 1.052976392459015155422d01, &
         1.945603779539281810439d01,-3.015761863840593359165d01, &
         1.120011024227297451523d01,-3.988850730390541057912d00, &
         9.565134591978630774217d00, 9.981193787537396413219d-1/
    data s/ 1.598517957704779356479d-4, 4.644185932583286942650d00, &
         3.697412299772985940785d02,-8.791401054875438925029d00, &
         7.608194509086645763123d02, 2.852397548119248700147d01, &
         4.731097187816050252967d02,-2.369210235636181001661d02, &
         1.249884822712447891440d00/
    !
    ! J-fraction coefficients for 12.0 <= X < 24.0
    !
    data p1/-1.647721172463463140042d00,-1.860092121726437582253d01, &
         -1.000641913989284829961d01,-2.105740799548040450394d01, &
         -9.134835699998742552432d-1,-3.323612579343962284333d01, &
         2.495487730402059440626d01, 2.652575818452799819855d01, &
         -1.845086232391278674524d00, 9.999933106160568739091d-1/
    data q1/ 9.792403599217290296840d01, 6.403800405352415551324d01, &
         5.994932325667407355255d01, 2.538819315630708031713d02, &
         4.429413178337928401161d01, 1.192832423968601006985d03, &
         1.991004470817742470726d02,-1.093556195391091143924d01, &
         1.001533852045342697818d00/
    !
    ! J-fraction coefficients for  X >= 24.0
    !
    data p2/ 1.75338801265465972390d02,-2.23127670777632409550d02, &
         -1.81949664929868906455d01,-2.79798528624305389340d01, &
         -7.63147701620253630855d00,-1.52856623636929636839d01, &
         -7.06810977895029358836d00,-5.00006640413131002475d00, &
         -3.00000000320981265753d00, 1.00000000000000485503d00/
    data q2/ 3.97845977167414720840d04, 3.97277109100414518365d00, &
         1.37790390235747998793d02, 1.17179220502086455287d02, &
         7.04831847180424675988d01,-1.20187763547154743238d01, &
         -7.99243595776339741065d00,-2.99999894040324959612d00, &
         1.99999999999048104167d00/

    x = arg

    if (x == zero) then

       ei = -xinf
       if (int == 2) ei = -ei

    else if ((x < zero) .or. (int == 2)) then
       !
       ! Calculate EI for negative argument or for E1.
       !
       y = abs(x)

       if (y <= one) then

          sump = a(7) * y + a(1)
          sumq = y + b(1)
          do i = 2, 6
             sump = sump * y + a(i)
             sumq = sumq * y + b(i)
          end do
          ei = log(y) - sump / sumq
          if (int == 3) ei = ei * exp(y)

       else if (y <= four) then

          w = one / y
          sump = c(1)
          sumq = d(1)
          do i = 2, 9
             sump = sump * w + c(i)
             sumq = sumq * w + d(i)
          end do
          ei = - sump / sumq
          if (int /= 3) ei = ei * exp(-y)

       else

          if ((y > xbig) .and. (int < 3)) then
             ei = zero
          else
             w = one / y
             sump = e(1)
             sumq = f(1)
             do i = 2, 10
                sump = sump * w + e(i)
                sumq = sumq * w + f(i)
             end do
             ei = -w * (one - w * sump / sumq )
             if (int /= 3) ei = ei * exp(-y)
          end if

       end if

       if (int == 2) ei = -ei

    else if (x < six) then
       !
       !  To improve conditioning, rational approximations are expressed
       !  in terms of Chebyshev polynomials for 0 <= X < 6, and in
       !  continued fraction form for larger X.
       !
       t = x + x
       t = t / three - two
       px(1) = zero
       qx(1) = zero
       px(2) = p(1)
       qx(2) = q(1)
       do i = 2, 9
          px(i+1) = t * px(i) - px(i-1) + p(i)
          qx(i+1) = t * qx(i) - qx(i-1) + q(i)
       end do
       sump = half * t * px(10) - px(9) + p(10)
       sumq = half * t * qx(10) - qx(9) + q(10)
       frac = sump / sumq
       xmx0 = (x - x01/x11) - x02

       if (abs(xmx0) >= p037) then

          ei = log(x/x0) + xmx0 * frac

          if (int == 3) ei = exp(-x) * ei

       else
          !
          !  Special approximation to  ln(X/X0)  for X close to X0
          !
          y = xmx0 / (x + x0)
          ysq = y*y
          sump = plg(1)
          sumq = ysq + qlg(1)
          do i = 2, 4
             sump = sump*ysq + plg(i)
             sumq = sumq*ysq + qlg(i)
          end do
          ei = (sump / (sumq*(x+x0)) + frac) * xmx0
          if (int == 3) ei = exp(-x) * ei

       end if

    else if (x < twelve) then

       frac = zero
       do i = 1, 9
          frac = s(i) / (r(i) + x + frac)
       end do
       ei = (r(10) + frac) / x
       if (int /= 3) ei = ei * exp(x)

    else if (x <= two4) then

       frac = zero
       do i = 1, 9
          frac = q1(i) / (p1(i) + x + frac)
       end do
       ei = (p1(10) + frac) / x
       if (int /= 3) ei = ei * exp(x)

    else

       if ((x >= xmax) .and. (int < 3)) then

          ei = xinf

       else

          y = one / x
          frac = zero
          do i = 1, 9
             frac = q2(i) / (p2(i) + x + frac)
          end do
          frac = p2(10) + frac
          ei = y + y * y * frac

          if (int /= 3) then
             if (x <= xmax-two4) then
                ei = ei * exp(x)
             else
                !
                !  Calculation reformulated to avoid premature overflow
                !
                ei = (ei * exp(x-fourty)) * exp40
             end if

          end if

       end if

    end if

    result = ei
    
    return
  end subroutine calcei

  subroutine gaulegf(x1, x2, x, w, n)

    implicit none
    integer, intent(in) :: n
    REAL(wp), intent(in) :: x1, x2
    REAL(wp), dimension(n), intent(out) :: x, w
    integer :: i, j, m
    REAL(wp) :: p1, p2, p3, pp, xl, xm, z, z1
    REAL(wp), parameter :: eps=3.0E-14_wp

    m = (n+1)/2
    xm = 0.5_wp*(x2+x1)
    xl = 0.5_wp*(x2-x1)
    do i=1,m
       z = cos(pi_g*(i-0.25_wp)/(n+0.5_wp))
       z1 = 0.0_wp
       do while(abs(z-z1) .gt. eps)
          p1 = 1.0_wp
          p2 = 0.0_wp
          do j=1,n
             p3 = p2
             p2 = p1
             p1 = ((2.0_wp*j-1.0_wp)*z*p2-(j-1.0_wp)*p3)/j
          end do
          pp = n*(z*p1-p2)/(z*z-1.0_wp)
          z1 = z
          z = z1 - p1/pp
       end do
       x(i) = xm - xl*z
       x(n+1-i) = xm + xl*z
       w(i) = (2.0_wp*xl)/((1.0_wp-z*z)*pp*pp)
       w(n+1-i) = w(i)
    end do
  end subroutine gaulegf

END MODULE geometry_2d
