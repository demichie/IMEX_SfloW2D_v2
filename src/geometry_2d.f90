!*********************************************************************
!> \brief Grid module
!
!> This module contains the variables and the subroutines related to 
!> the computational grid
!*********************************************************************
MODULE geometry_2d

  USE parameters_2d, ONLY : verbose_level

  IMPLICIT NONE

  !> Location of the centers (x) of the control volume of the domain
  REAL*8, ALLOCATABLE :: x_comp(:)

  !> Location of the boundaries (x) of the control volumes of the domain
  REAL*8, ALLOCATABLE :: x_stag(:)

  !> Location of the centers (y) of the control volume of the domain
  REAL*8, ALLOCATABLE :: y_comp(:)

  !> Location of the boundaries (x) of the control volumes of the domain
  REAL*8, ALLOCATABLE :: y_stag(:)

  !> Reconstructed value at the left of the x-interface
  REAL*8, ALLOCATABLE :: B_interfaceL(:,:)        

  !> Reconstructed value at the right of the x-interface
  REAL*8, ALLOCATABLE :: B_interfaceR(:,:)

  !> Reconstructed value at the bottom of the y-interface
  REAL*8, ALLOCATABLE :: B_interfaceB(:,:)        

  !> Reconstructed value at the top of the y-interface
  REAL*8, ALLOCATABLE :: B_interfaceT(:,:)

  !> Topography at the centers of the control volumes 
  REAL*8, ALLOCATABLE :: B_cent(:,:)

  !> Topography slope (x direction) at the centers of the control volumes 
  REAL*8, ALLOCATABLE :: B_prime_x(:,:)

  !> Topography slope (y direction) at the centers of the control volumes 
  REAL*8, ALLOCATABLE :: B_prime_y(:,:)

  !> Solution in ascii grid format (ESRI)
  REAL*8, ALLOCATABLE :: grid_output(:,:)

  !> gravity vector wrt surface coordinates for each cell
  REAL*8, ALLOCATABLE :: grav_surf(:,:)

  !> curvature wrt mixed directions for each cell
  REAL*8, ALLOCATABLE :: curv_xy(:,:)

  REAL*8, ALLOCATABLE :: topography_profile(:,:,:)

  INTEGER :: n_topography_profile_x, n_topography_profile_y

  REAL*8 :: dx                 !< Control volumes size
  REAL*8 :: x0                 !< Left of the physical domain
  REAL*8 :: xN                 !< Right of the physical domain
  REAL*8 :: dy                 !< Control volumes size
  REAL*8 :: y0                 !< Bottom of the physical domain
  REAL*8 :: yN                 !< Top of the physical domain
  REAL*8 :: dx2                !< Half x Control volumes size
  REAL*8 :: dy2                !< Half y Control volumes size
  INTEGER :: comp_cells_x      !< Number of control volumes x in the comp. domain
  INTEGER :: comp_interfaces_x !< Number of interfaces (comp_cells_x+1)
  INTEGER :: comp_cells_y      !< Number of control volumes y in the comp. domain
  INTEGER :: comp_interfaces_y !< Number of interfaces (comp_cells_y+1)
  REAL*8 :: cell_size
  INTEGER :: comp_cells_xy

CONTAINS

  !*********************************************************************
  !> \brief Finite volume grid initialization
  !
  !> This subroutine initialize the grids for the finite volume solver.
  !> \date 16/08/2011
  !*********************************************************************

  SUBROUTINE init_grid

    USE parameters_2d, ONLY: eps_sing

    IMPLICIT none

    INTEGER j,k      !> loop counter

    comp_interfaces_x = comp_cells_x+1
    comp_interfaces_y = comp_cells_y+1

    ALLOCATE( x_comp(comp_cells_x) )
    ALLOCATE( x_stag(comp_interfaces_x) )
    ALLOCATE( y_comp(comp_cells_y) )
    ALLOCATE( y_stag(comp_interfaces_y) )

    ALLOCATE( B_cent(comp_cells_x,comp_cells_y) )
    ALLOCATE( B_prime_x(comp_cells_x,comp_cells_y) )
    ALLOCATE( B_prime_y(comp_cells_x,comp_cells_y) )

    ALLOCATE( B_interfaceL( comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( B_interfaceR( comp_interfaces_x, comp_cells_y ) )
    ALLOCATE( B_interfaceB( comp_cells_x, comp_interfaces_y ) )
    ALLOCATE( B_interfaceT( comp_cells_x, comp_interfaces_y ) )

    ALLOCATE( grid_output(comp_cells_x,comp_cells_y) )

    ALLOCATE( grav_surf(comp_cells_x,comp_cells_y) )

    IF ( comp_cells_x .GT. 1 ) THEN

       dx = cell_size

    ELSE
       
       dx = 1.D0

    END IF

    
    IF ( comp_cells_y .GT. 1 ) THEN

       dy = cell_size
    
    ELSE

       dy = 1.D0

    END IF

    xN = x0 + comp_cells_x * dx
    yN = y0 + comp_cells_y * dy

    dx2 = dx / 2.d0
    dy2 = dy / 2.d0

    ! eps_sing = MIN( dx ** 4.D0,dy ** 4.D0 )
    eps_sing=MIN(MIN( dx ** 4.D0,dy ** 4.D0 ),1.d-20)

    IF ( verbose_level .GE. 1 ) WRITE(*,*) 'eps_sing = ',eps_sing
    
    DO j=1,comp_interfaces_x

       x_stag(j) = x0 + (j-1) * dx

    END DO

    DO k=1,comp_interfaces_y

       y_stag(k) = y0 + (k-1) * dy

    END DO

    DO j=1,comp_cells_x

       x_comp(j) = 0.5D0 * ( x_stag(j) + x_stag(j+1) )

    END DO

    DO k=1,comp_cells_y

       y_comp(k) = 0.5D0 * ( y_stag(k) + y_stag(k+1) )

    END DO

    DO k=1,comp_cells_y
    
       DO j=1,comp_cells_x

          CALL interp_2d_scalar( topography_profile(1,:,:) ,                    &
               topography_profile(2,:,:), topography_profile(3,:,:) ,           &
               x_comp(j), y_comp(k) , B_cent(j,k) )

          CALL interp_2d_slope( topography_profile(1,:,:) ,                     &
               topography_profile(2,:,:), topography_profile(3,:,:) ,           &
               x_comp(j), y_comp(k) , B_prime_x(j,k) , B_prime_y(j,k) )

          B_interfaceR(j,k) = B_cent(j,k) - dx2 * B_prime_x(j,k)
          B_interfaceL(j+1,k) = B_cent(j,k) + dx2 * B_prime_x(j,k)
          
          B_interfaceT(j,k) = B_cent(j,k) - dx2 * B_prime_y(j,k)
          B_interfaceB(j,k+1) = B_cent(j,k) + dx2 * B_prime_y(j,k)

       END DO

    ENDDO

    ! this coefficient is used when the the scalar dot between the normal to the 
    ! topography and gravity is computed
    DO j = 1,comp_cells_x

       DO k=1,comp_cells_y

          grav_surf(j,k) = - ( 1.d0/ DSQRT( 1.d0 + B_prime_x(j,k)**2            & 
               + B_prime_y(j,k)**2) )

       ENDDO

    ENDDO

  END SUBROUTINE init_grid

  !---------------------------------------------------------------------------
  !> Scalar interpolation
  !
  !> This subroutine interpolate the values of the  array f1, defined on the 
  !> grid points x1, at the point x2. The value are saved in f2
  !> \date 13/02/2009
  !> \param    x1           original grid                (\b input)
  !> \param    f1           original values              (\b input)
  !> \param    x2           new point                    (\b output)
  !> \param    f2           interpolated value           (\b output)
  !---------------------------------------------------------------------------

  SUBROUTINE interp_1d_scalar(x1, f1, x2, f2)
    IMPLICIT NONE

    REAL*8, INTENT(IN), DIMENSION(:) :: x1, f1
    REAL*8, INTENT(IN) :: x2
    REAL*8, INTENT(OUT) :: f2
    INTEGER :: n, n1x, t
    REAL*8 :: grad , rel_pos

    n1x = SIZE(x1)

    !
    ! ... locate the grid points near the topographic points
    ! ... and interpolate linearly the profile
    !
    t = 1

    search:DO n = 1, n1x-1

       rel_pos = ( x2 - x1(n) ) / ( x1(n+1) - x1(n) )

       IF ( ( rel_pos .GE. 0.D0 ) .AND. ( rel_pos .LE. 1.D0 ) ) THEN

          grad = ( f1(n+1)-f1(n) ) / ( x1(n+1)-x1(n) )
          f2 = f1(n) + ( x2-x1(n) ) * grad

          EXIT search

       ELSEIF  ( rel_pos .LT. 0.D0 ) THEN

          f2 = f1(n)

       ELSE

          f2 = f1(n+1)

       END IF

    END DO search

    RETURN

  END SUBROUTINE interp_1d_scalar


  !---------------------------------------------------------------------------
  !> Scalar interpolation (2D)
  !
  !> This subroutine interpolate the values of the  array f1, defined on the 
  !> grid points (x1,y1), at the point (x2,y2). The value are saved in f2
  !> \date OCTOBER 2016
  !> \param    x1           original grid                (\b input)
  !> \param    y1           original grid                (\b input)
  !> \param    f1           original values              (\b input)
  !> \param    x2           new point                    (\b output)
  !> \param    y2           new point                    (\b output)
  !> \param    f2           interpolated value           (\b output)
  !---------------------------------------------------------------------------

  SUBROUTINE interp_2d_scalar(x1, y1, f1, x2, y2, f2)
    IMPLICIT NONE

    REAL*8, INTENT(IN), DIMENSION(:,:) :: x1, y1, f1
    REAL*8, INTENT(IN) :: x2, y2
    REAL*8, INTENT(OUT) :: f2

    INTEGER :: ix , iy
    REAL*8 :: alfa_x , alfa_y

    IF ( size(x1,1) .GT. 1 ) THEN

       ix = FLOOR( ( x2 - x1(1,1) ) / ( x1(2,1) - x1(1,1) ) ) + 1
       ix = MIN( ix , SIZE(x1,1)-1 )
       alfa_x = ( x1(ix+1,1) - x2 ) / (  x1(ix+1,1) - x1(ix,1) )

    ELSE

       ix = 1
       alfa_x = 0.D0
       
    END IF
    
    IF ( size(x1,2) .GT. 1 ) THEN

       iy = FLOOR( ( y2 - y1(1,1) ) / ( y1(1,2) - y1(1,1) ) ) + 1
       iy = MIN( iy , SIZE(x1,2)-1 )
       alfa_y = ( y1(1,iy+1) - y2 ) / (  y1(1,iy+1) - y1(1,iy) )
    
    ELSE

       iy = 1
       alfa_y = 0.D0
       
    END IF
       
    IF ( size(x1,1) .EQ. 1 ) THEN
       
       f2 = alfa_y * f1(ix,iy) + ( 1.D0 - alfa_y ) * f1(ix,iy+1)
       
    ELSEIF ( size(x1,2) .EQ. 1 ) THEN
       
       f2 = alfa_x * f1(ix,iy)  + ( 1.D0 - alfa_x ) * f1(ix+1,iy)
       
    ELSE
       
       f2 = alfa_x * ( alfa_y * f1(ix,iy) + ( 1.D0 - alfa_y ) * f1(ix,iy+1) )   &
            + ( 1.D0 - alfa_x ) *  ( alfa_y * f1(ix+1,iy) + ( 1.D0 - alfa_y )   &
            * f1(ix+1,iy+1) )
       
    END IF
    
  END SUBROUTINE interp_2d_scalar

  !---------------------------------------------------------------------------
  !> Scalar interpolation (2D)
  !
  !> This subroutine interpolate the values of the  array f1, defined on the 
  !> grid points (x1,y1), at the point (x2,y2). The value are saved in f2.
  !> In this case x1 and y1 are 1d arrays.
  !> \date OCTOBER 2016
  !> \param    x1           original grid                (\b input)
  !> \param    y1           original grid                (\b input)
  !> \param    f1           original values              (\b input)
  !> \param    x2           new point                    (\b output)
  !> \param    y2           new point                    (\b output)
  !> \param    f2           interpolated value           (\b output)
  !---------------------------------------------------------------------------

  SUBROUTINE interp_2d_scalarB(x1, y1, f1, x2, y2, f2)
    IMPLICIT NONE

    REAL*8, INTENT(IN), DIMENSION(:) :: x1, y1
    REAL*8, INTENT(IN), DIMENSION(:,:) :: f1
    REAL*8, INTENT(IN) :: x2, y2
    REAL*8, INTENT(OUT) :: f2

    INTEGER :: ix , iy
    REAL*8 :: alfa_x , alfa_y

    IF ( size(x1) .GT. 1 ) THEN

       ix = FLOOR( ( x2 - x1(1) ) / ( x1(2) - x1(1) ) ) + 1
       ix = MAX(0,MIN( ix , SIZE(x1)-1 ))
       alfa_x = ( x1(ix+1) - x2 ) / (  x1(ix+1) - x1(ix) )

    ELSE

       ix = 1
       alfa_x = 0.D0
       
    END IF
    
    IF ( size(y1) .GT. 1 ) THEN

       iy = FLOOR( ( y2 - y1(1) ) / ( y1(2) - y1(1) ) ) + 1
       iy = MAX(1,MIN( iy , SIZE(y1)-1 ))
       alfa_y = ( y1(iy+1) - y2 ) / (  y1(iy+1) - y1(iy) )
    
    ELSE

       iy = 1
       alfa_y = 0.D0
       
    END IF

    IF ( ( alfa_x .LT. 0.D0 ) .OR. ( alfa_x .GT. 1.D0 )                         &
         .OR. ( alfa_y .LT. 0.D0 ) .OR. ( alfa_y .GT. 1.D0 ) ) THEN

       f2 = 0.D0
       RETURN

    END IF
       
    
    IF ( size(x1) .EQ. 1 ) THEN
       
       f2 = alfa_y * f1(ix,iy) + ( 1.D0 - alfa_y ) * f1(ix,iy+1)
       
    ELSEIF ( size(y1) .EQ. 1 ) THEN
       
       f2 = alfa_x * f1(ix,iy)  + ( 1.D0 - alfa_x ) * f1(ix+1,iy)
       
    ELSE
       
       f2 = alfa_x * ( alfa_y * f1(ix,iy) + ( 1.D0 - alfa_y ) * f1(ix,iy+1) )   &
            + ( 1.D0 - alfa_x ) *  ( alfa_y * f1(ix+1,iy) + ( 1.D0 - alfa_y )   &
            * f1(ix+1,iy+1) )
       
    END IF

    RETURN
    
  END SUBROUTINE interp_2d_scalarB


  !---------------------------------------------------------------------------
  !> Scalar interpolation (2D)
  !
  !> This subroutine interpolate the values of the  array f1, defined on the 
  !> grid points (x1,y1), at the point (x2,y2). The value are saved in f2
  !> \date OCTOBER 2016
  !> \param    x1           original grid                (\b input)
  !> \param    y1           original grid                (\b input)
  !> \param    f1           original values              (\b input)
  !> \param    x2           new point                    (\b output)
  !> \param    y2           new point                    (\b output)
  !> \param    f2           interpolated value           (\b output)
  !---------------------------------------------------------------------------

  SUBROUTINE interp_2d_slope(x1, y1, f1, x2, y2, f_x, f_y)
    IMPLICIT NONE

    REAL*8, INTENT(IN), DIMENSION(:,:) :: x1, y1, f1
    REAL*8, INTENT(IN) :: x2, y2
    REAL*8, INTENT(OUT) :: f_x , f_y

    INTEGER :: ix , iy
    REAL*8 :: alfa_x , alfa_y
    REAL*8 :: f_x1 , f_x2 , f_y1 , f_y2


    IF ( size(x1,1) .GT. 1 ) THEN

       ix = FLOOR( ( x2 - x1(1,1) ) / ( x1(2,1) - x1(1,1) ) ) + 1
       ix = MIN( ix , SIZE(x1,1)-1 )
       alfa_x = ( x1(ix+1,1) - x2 ) / (  x1(ix+1,1) - x1(ix,1) )

    ELSE

       ix = 1
       alfa_x = 1.D0
       
    END IF
    
    IF ( size(x1,2) .GT. 1 ) THEN

       iy = FLOOR( ( y2 - y1(1,1) ) / ( y1(1,2) - y1(1,1) ) ) + 1
       iy = MIN( iy , SIZE(x1,2)-1 )
       alfa_y = ( y1(1,iy+1) - y2 ) / (  y1(1,iy+1) - y1(1,iy) )
    
    ELSE

       iy = 1
       alfa_y = 1.D0
  
    END IF

    f_x1 = 0.D0
    f_x2 = 0.D0

    IF ( size(x1,1) .GT. 1 ) THEN

       f_x1 = ( f1(ix+1,iy) -  f1(ix,iy) ) / (  x1(2,1) - x1(1,1) )

       IF ( size(x1,2) .GT. 1 ) THEN

          f_x2 = ( f1(ix+1,iy+1) -  f1(ix,iy+1) ) / (  x1(2,1) - x1(1,1) )

       END IF

    END IF

    f_x = alfa_y * f_x1 + ( 1.D0 - alfa_y ) * f_x2

    f_y1 = 0.D0
    f_y2 = 0.D0

    IF ( size(x1,2) .GT. 1 ) THEN

       f_y1 = ( f1(ix,iy+1) - f1(ix,iy) ) / (  y1(1,2) - y1(1,1) )

       IF ( size(x1,1) .GT. 1 ) THEN

          f_y2 = ( f1(ix+1,iy+1) - f1(ix+1,iy) ) / (  y1(1,2) - y1(1,1) )

       END IF

    END IF

    f_y = alfa_x * f_y1 + ( 1.D0 - alfa_x ) * f_y2


    RETURN
    
  END SUBROUTINE interp_2d_slope


  !******************************************************************************
  !> \brief Linear reconstruction
  !
  !> In this subroutine a linear reconstruction with slope limiters is
  !> applied to a set of variables describing the state of the system.
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 2019/11/08
  !******************************************************************************

  SUBROUTINE topography_reconstruction

    USE parameters_2d, ONLY : limiter

    IMPLICIT NONE
    
    REAL*8 :: B_stencil(3)    !< recons variables stencil for the limiter
    REAL*8 :: x_stencil(3)    !< grid stencil for the limiter
    REAL*8 :: y_stencil(3)    !< grid stencil for the limiter

    INTEGER :: limiterB

    INTEGER :: j,k
    
    
    ! When limiter(1) = 0 (slope=0) we use a minmod limiter
    limiterB = MAX(1,limiter(1))
    limiterB = 1.D0

    y_loop:DO k = 1,comp_cells_y

       x_loop:DO j = 1,comp_cells_x

          ! x direction
          check_comp_cells_x:IF ( comp_cells_x .GT. 1 ) THEN

             check_x_boundary:IF (j.EQ.1) THEN

                ! west boundary

                x_stencil(1) = 2.D0 * x_comp(1) - x_comp(2)
                x_stencil(2:3) = x_comp(1:2)

                B_stencil(1) = 2.D0 * B_cent(1,k) - B_cent(2,k) 
                B_stencil(2:3) = B_cent(1:2,k)

                CALL limit( B_stencil , x_stencil , limiterB , B_prime_x(j,k) )

             ELSEIF (j.EQ.comp_cells_x) THEN

                !east boundary
                
                x_stencil(3) = 2.D0 * x_comp(comp_cells_x) - x_comp(comp_cells_x-1)
                x_stencil(1:2) = x_comp(comp_cells_x-1:comp_cells_x)
                
                B_stencil(3) = 2.D0 * B_cent(comp_cells_x,k) - B_cent(comp_cells_x-1,k)
                B_stencil(1:2) = B_cent(comp_cells_x-1:comp_cells_x,k)

                CALL limit( B_stencil , x_stencil , limiterB , B_prime_x(j,k) ) 

             ELSE

                ! Internal x interfaces
                x_stencil(1:3) = x_comp(j-1:j+1)
                B_stencil = B_cent(j-1:j+1,k)

                CALL limit( B_stencil , x_stencil , limiterB , B_prime_x(j,k) )

             ENDIF check_x_boundary

          ELSE

             B_prime_x(j,k) = 0.D0
             
          END IF check_comp_cells_x

          ! y-direction
          check_comp_cells_y:IF ( comp_cells_y .GT. 1 ) THEN

             check_y_boundary:IF (k.EQ.1) THEN
                
                ! South boundary
                y_stencil(1) = 2.D0 * y_comp(1) - y_comp(2)
                y_stencil(2:3) = y_comp(1:2)
                
                B_stencil(1) = 2.D0 * B_cent(j,1) - B_cent(j,2)
                B_stencil(2:3) = B_cent(j,1:2)
                
                CALL limit( B_stencil , y_stencil , limiterB , B_prime_y(j,k) ) 
                
             ELSEIF ( k .EQ. comp_cells_y ) THEN

                ! North boundary
                y_stencil(3) = 2.D0 * y_comp(comp_cells_y) - y_comp(comp_cells_y-1)
                y_stencil(1:2) = y_comp(comp_cells_y-1:comp_cells_y)

                B_stencil(3) = 2.D0 * B_cent(j,comp_cells_y) - B_cent(j,comp_cells_y-1)
                B_stencil(1:2) = B_cent(j,comp_cells_y-1:comp_cells_y)
                
                CALL limit( B_stencil , y_stencil , limiterB , B_prime_y(j,k) ) 
                
             ELSE

                ! Internal y interfaces
                y_stencil(1:3) = y_comp(k-1:k+1)
                B_stencil = B_cent(j,k-1:k+1)

                CALL limit( B_stencil , y_stencil , limiterB , B_prime_y(j,k) )

             ENDIF check_y_boundary

          ELSE

             B_prime_y(j,k) = 0.D0
             
          ENDIF check_comp_cells_y

       END DO x_loop
       
    END DO y_loop

    DO j = 1,comp_cells_x
    
       DO k = 1,comp_cells_y

          B_interfaceR(j,k) = B_cent(j,k) - dx2 * B_prime_x(j,k)
          B_interfaceL(j+1,k) = B_cent(j,k) + dx2 * B_prime_x(j,k)
       
          B_interfaceT(j,k) = B_cent(j,k) - dx2 * B_prime_y(j,k)
          B_interfaceB(j,k+1) = B_cent(j,k) + dx2 * B_prime_y(j,k)

       END DO

    END DO

    RETURN

  END SUBROUTINE topography_reconstruction

  !---------------------------------------------------------------------------
  !> Scalar regrid (2D)
  !
  !> This subroutine interpolate the values of the  array f1, defined on the 
  !> grid points (x1,y1), at the point (x2,y2). The value are saved in f2.
  !> In this case x1 and y1 are 1d arrays.
  !> \date OCTOBER 2016
  !> \param    x1           original grid                (\b input)
  !> \param    y1           original grid                (\b input)
  !> \param    f1           original values              (\b input)
  !> \param    xl           new point                    (\b input)
  !> \param    xr           new point                    (\b input)
  !> \param    yl           new point                    (\b input)
  !> \param    yr           new point                    (\b input)
  !> \param    f2           interpolated value           (\b output)
  !---------------------------------------------------------------------------

  SUBROUTINE regrid_scalar(xin, yin, fin, xl, xr , yl, yr, fout)
    IMPLICIT NONE

    REAL*8, INTENT(IN), DIMENSION(:) :: xin, yin
    REAL*8, INTENT(IN), DIMENSION(:,:) :: fin
    REAL*8, INTENT(IN) :: xl, xr , yl , yr
    REAL*8, INTENT(OUT) :: fout

    INTEGER :: ix , iy
    INTEGER :: ix1 , ix2 , iy1 , iy2
    REAL*8 :: alfa_x , alfa_y
    REAL*8 :: dXin , dYin
    
    INTEGER nXin,nYin
    
    nXin = size(xin)-1
    nYin = size(yin)-1

    dXin = xin(2) - xin(1)
    dYin = yin(2) - yin(1)
    
    ix1 = MAX(1,CEILING( ( xl - xin(1) ) / dXin ))
    ix2 = MIN(nXin,CEILING( ( xr -xin(1) ) / dXin )+1)
        
    iy1 = MAX(1,CEILING( ( yl - yin(1) ) / dYin ))
    iy2 = MIN(nYin,CEILING( ( yr - yin(1) ) / dYin ) + 1)

    fout = 0.D0
    
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

    REAL*8, INTENT(IN) :: v(3)
    REAL*8, INTENT(IN) :: z(3)
    INTEGER, INTENT(IN) :: limiter

    REAL*8, INTENT(OUT) :: slope_lim

    REAL*8 :: a , b , c

    REAL*8 :: sigma1 , sigma2

    a = ( v(3) - v(2) ) / ( z(3) - z(2) )
    b = ( v(2) - v(1) ) / ( z(2) - z(1) )
    c = ( v(3) - v(1) ) / ( z(3) - z(1) )

    SELECT CASE (limiter)

    CASE ( 0 )

       slope_lim = 0.D0

    CASE ( 1 )

       ! minmod

       slope_lim = minmod(a,b)

    CASE ( 2 )

       ! superbee

       sigma1 = minmod( a , 2.D0*b )
       sigma2 = minmod( 2.D0*a , b )
       slope_lim = maxmod( sigma1 , sigma2 )

    CASE ( 3 )

       ! generalized minmod

       slope_lim = minmod( c , theta * minmod( a , b ) )

    CASE ( 4 )

       ! monotonized central-difference (MC, LeVeque p.112)

       slope_lim = minmod( c , 2.0 * minmod( a , b ) )

    END SELECT

  END SUBROUTINE limit


  REAL*8 FUNCTION minmod(a,b)

    IMPLICIT none

    REAL*8 :: a , b , sa , sb 

    IF ( a*b .EQ. 0.D0 ) THEN

       minmod = 0.d0

    ELSE

       sa = a / ABS(a)
       sb = b / ABS(b)

       minmod = 0.5D0 * ( sa+sb ) * MIN( ABS(a) , ABS(b) )

    END IF

  END FUNCTION minmod

  REAL*8 function maxmod(a,b)

    IMPLICIT none

    REAL*8 :: a , b , sa , sb 

    IF ( a*b .EQ. 0.d0 ) THEN

       maxmod = 0.d0

    ELSE

       sa = a / ABS(a)
       sb = b / ABS(b)

       maxmod = 0.5D0 * ( sa+sb ) * MAX( ABS(a) , ABS(b) )

    END IF

  END function maxmod


END MODULE geometry_2d
