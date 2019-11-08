!********************************************************************************
!> \brief Constitutive equations
!********************************************************************************
MODULE constitutive_2d

  USE parameters_2d, ONLY : n_eqns , n_vars , n_solid
  USE parameters_2d, ONLY : rheology_flag , rheology_model , energy_flag
  USE parameters_2d, ONLY : sed_vol_perc

  IMPLICIT none

  LOGICAL, ALLOCATABLE :: implicit_flag(:)

  CHARACTER(LEN=20) :: phase1_name
  CHARACTER(LEN=20) :: phase2_name

  COMPLEX*16 :: h                       !< height [m]
  COMPLEX*16 :: u                       !< velocity (x direction)
  COMPLEX*16 :: v                       !< velocity (y direction)
  COMPLEX*16 :: T                       !< temperature
  COMPLEX*16, ALLOCATABLE :: alphas(:)  !< sediment volume fraction
  COMPLEX*16 :: rho_m                   !< mixture density
  COMPLEX*16 :: Ri
  
  COMPLEX*16 :: alphas_tot
  COMPLEX*16 :: rhos_tot

  !> gravitational acceleration
  REAL*8 :: grav

  !> reduced gravity
  COMPLEX*16 :: red_grav

  REAL*8 :: r_h
  REAL*8 :: r_u
  REAL*8 :: r_v
  REAL*8 :: r_T
  REAL*8, ALLOCATABLE :: r_alphas(:)
  REAL*8 :: r_red_grav
  REAL*8 :: r_rho_m
  REAL*8 :: r_rho_a
  REAL*8 :: r_Ri


  !> drag coefficients (Voellmy-Salm model)
  REAL*8 :: mu
  REAL*8 :: xi

  !> drag coefficients (B&W model)
  REAL*8 :: friction_factor
  
  !> drag coefficients (plastic model)
  REAL*8 :: tau

  !> evironment temperature [K]
  REAL*8 :: T_env

  !> radiative coefficient
  REAL*8 :: rad_coeff

  !> friction coefficient
  REAL*8 :: frict_coeff

  !> fluid density [kg/m3]
  REAL*8 :: rho

  !> reference temperature [K]
  REAL*8 :: T_ref

  !> reference kinematic viscosity [m2/s]
  REAL*8 :: nu_ref

  !> viscosity parameter [K-1] (b in Table 1 Costa & Macedonio, 2005)
  REAL*8 :: visc_par

  !> velocity boundary layer fraction of total thickness
  REAL*8 :: emme

  !> specific heat [J kg-1 K-1]
  REAL*8 :: c_p

  !> atmospheric heat trasnfer coefficient [W m-2 K-1] (lambda in C&M, 2005)
  REAL*8 :: atm_heat_transf_coeff

  !> fractional area of the exposed inner core (f in C&M, 2005)
  REAL*8 :: exp_area_fract

  !> Stephan-Boltzmann constant [W m-2 K-4]
  REAL*8, PARAMETER :: SBconst = 5.67D-8

  !> emissivity (eps in Costa & Macedonio, 2005)
  REAL*8 :: emissivity

  !> thermal boundary layer fraction of total thickness
  REAL*8 :: enne

  !> temperature of lava-ground interface
  REAL*8 :: T_ground

  !> thermal conductivity [W m-1 K-1] (k in Costa & Macedonio, 2005)
  REAL*8 :: thermal_conductivity

  !--- Lahars rheology model parameters

  !> 1st parameter for yield strenght empirical relationship (O'Brian et al, 1993)
  REAL*8 :: alpha2

  !> 2nd parameter for yield strenght empirical relationship (O'Brian et al, 1993)
  REAL*8 :: beta2

  !> 1st parameter for fluid viscosity empirical relationship (O'Brian et al, 1993)
  REAL*8 :: alpha1

  !> 2nd parameter for fluid viscosity empirical relationship (O'Brian et al, 1993)
  REAL*8 :: beta1

  !> Empirical resistance parameter
  REAL*8 :: Kappa

  !> Mannings roughness coefficient ( units: T L^(-1/3) )
  REAL*8 :: n_td
 
  !> Density of air in the mixture
  COMPLEX*16 :: rho_a

  !> Ambient density of air
  REAL*8 :: rho_a_amb

  !> Specific heat of air
  REAL*8 :: sp_heat_a

  !> Specific gas constant of air
  REAL*8 :: sp_gas_const_a

  !> Kinematic viscosity of air
  REAL*8 :: kin_visc_a

  !> Temperature of ambient air 
  REAL*8 :: T_ambient

  LOGICAL :: entrainment_flag
  
  !> Specific heat of mixture
  COMPLEX*16 :: sp_heat_mix

  !> Specific weight of sediments
  REAL*8, ALLOCATABLE :: rho_s(:)

  !> Diameter of sediments
  REAL*8, ALLOCATABLE :: diam_s(:)

  !> Specific heat of sediments
  REAL*8, ALLOCATABLE :: sp_heat_s(:)

  REAL*8, ALLOCATABLE :: C_D_s(:)

  !> hindered settling 
  REAL*8 :: settling_vel

  !> erosion model coefficient
  REAL*8 :: erosion_coeff

  !> temperature of solid substrate (K)
  REAL*8 :: T_s_substrate

  !> ambient pressure
  REAL*8 :: pres

  !> liquid density
  REAL*8 :: rho_l

  !> Sepcific heat of liquid
  REAL*8 :: sp_heat_l

  !> liquid volume fraction
  COMPLEX*16 :: alpha_l
  
CONTAINS

  !******************************************************************************
  !> \brief Initialization of relaxation flags
  !
  !> This subroutine set the number and the flags of the non-hyperbolic
  !> terms.
  !> \date 07/09/2012
  !******************************************************************************

  SUBROUTINE init_problem_param

    USE parameters_2d, ONLY : n_nh

    ALLOCATE( implicit_flag(n_eqns) )

    implicit_flag(1:n_eqns) = .FALSE.
    implicit_flag(2) = .TRUE.
    implicit_flag(3) = .TRUE.

    ! Temperature
    implicit_flag(4) = .TRUE.

    ! Solid volume fraction
    implicit_flag(5:4+n_solid) = .TRUE.

    n_nh = COUNT( implicit_flag )

  END SUBROUTINE init_problem_param

  !******************************************************************************
  !> \brief Physical variables
  !
  !> This subroutine evaluates from the conservative local variables qj
  !> the local physical variables  (\f$h+B, u, v, xs , T \f$).
  !> \param[in]    r_qj     real conservative variables 
  !> \param[in]    c_qj     complex conservative variables 
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 15/08/2011
  !******************************************************************************

  SUBROUTINE phys_var(r_qj,c_qj)
    
    USE COMPLEXIFY
    USE parameters_2d, ONLY : eps_sing
    IMPLICIT none
    
    REAL*8, INTENT(IN), OPTIONAL :: r_qj(n_vars)
    COMPLEX*16, INTENT(IN), OPTIONAL :: c_qj(n_vars)

    COMPLEX*16 :: qj(n_vars)

    IF ( present(c_qj) ) THEN

       qj = c_qj

    ELSE

       qj = DCMPLX(r_qj)

    END IF

    IF ( qj(1) .GT. 1.D-25 ) THEN

       sp_heat_mix = SUM( qj(5:4+n_solid) / qj(1) * sp_heat_s(1:n_solid) ) +    &
            ( DCMPLX(1.D0,0.D0) - SUM( qj(5:4+n_solid) / qj(1) ) ) * sp_heat_a

       IF ( energy_flag ) THEN

          T = ( qj(4) - ( qj(2)**2 + qj(3)**2 ) / ( 2.D0*qj(1) ) ) /            &
               ( qj(1) * sp_heat_mix ) 
 
       ELSE

          T = qj(4) / ( qj(1) * sp_heat_mix ) 

       END IF
   
       IF ( DBLE(T) .LE. 0.D0 ) T = DCMPLX(T_ambient,0.D0)

    ELSE

       sp_heat_mix = DCMPLX(sp_heat_a,0.D0)

       T = DCMPLX(T_ambient,0.D0)

    END IF

    rho_a =  pres / ( sp_gas_const_a * T )

    IF ( SUM(qj(5:4+n_solid)) .GT. 1.D-25 ) THEN

       rhos_tot = SUM( qj(5:4+n_solid) ) / SUM( qj(5:4+n_solid) /               &
            rho_s(1:n_solid) ) 

    ELSE

       rhos_tot = SUM( rho_s(1:n_solid) ) / n_solid

    END IF

    IF ( SUM(qj(5:4+n_solid)/rho_s(1:n_solid)) .LT. 1.D-25 ) THEN
    
       alphas_tot = DCMPLX(0.D0,0.D0)
       
    ELSE

       alphas_tot = rho_a / ( rho_a - rhos_tot * ( DCMPLX(1.D0,0.D0) - qj(1) /  &
            SUM( qj(5:4+n_solid) ) ) )
    
    END IF

    IF ( DBLE(alphas_tot ) .LT. 0.D0 ) THEN

       alphas_tot = DCMPLX(0.D0,0.D0)
       h = qj(1) / rho_a
       alphas(1:n_solid) =  DCMPLX(0.D0,0.D0)
       rho_m = rho_a

    ELSE

       h = qj(1) / ( ( DCMPLX(1.D0,0.D0) - alphas_tot ) * rho_a + alphas_tot *     &
            rhos_tot )

       IF ( DBLE( h ) .GT. eps_sing ) THEN
          
          alphas(1:n_solid) = qj(5:4+n_solid) / ( rho_s(1:n_solid) * h )
          
       ELSE
          
          alphas(1:n_solid) = DSQRT(2.D0) * h * qj(5:4+n_solid) / CDSQRT( h**4     &
               + eps_sing**4 ) / rho_s(1:n_solid)
                 
       END IF

       IF ( DBLE( SUM( alphas(1:n_solid) )) .GT. 1.D0 ) THEN

          alphas(1:n_solid) = alphas(1:n_solid) / DBLE( SUM( alphas(1:n_solid) ))

       END IF
    
       rho_m = ( DCMPLX(1.D0,0.D0) - SUM(alphas) ) * rho_a + SUM( alphas * rho_s ) 

    END IF

    red_grav = ( rho_m - rho_a_amb ) / rho_m * grav
      
    IF ( DBLE( qj(1) ) .GT. eps_sing ) THEN
       
       u = qj(2) / qj(1)
       v = qj(3) / qj(1)

    ELSE

       u = DSQRT(2.D0) * qj(1) * qj(2) / CDSQRT( qj(1)**4 + eps_sing**4 )
       v = DSQRT(2.D0) * qj(1) * qj(3) / CDSQRT( qj(1)**4 + eps_sing**4 )
       
    END IF

    IF ( DBLE( u**2 + v**2 ) .GT. 0.D0 ) THEN
    
       Ri = red_grav * h / ( u**2 + v**2 )

    ELSE

       Ri = DCMPLX(0.D0,0.D0)

    END IF

    RETURN

  END SUBROUTINE phys_var

  
  !******************************************************************************
  !> \brief Physical variables
  !
  !> This subroutine evaluates from the physical local variables qpj, the
  !> conservative local variables qcj and the local topography Bj a set of 
  !> the local real-valued physical variables  (\f$h, u, v, xs , T \f$).
  !> \param[in]    qpj    real-valued physical variables 
  !> \param[in]    qcj    real-valued conservative variables 
  !> \param[in]    Bj     real-valued topography 
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 10/10/2019
  !******************************************************************************
 
  SUBROUTINE mixt_var(qpj)

    IMPLICIT none

    REAL*8, INTENT(IN) :: qpj(n_vars)
    ! REAL*8, INTENT(IN) :: qcj(n_vars)

    r_h = qpj(1)

    r_T = qpj(4)

    r_alphas(1:n_solid) = qpj(5:4+n_solid) 

    r_rho_a =  pres / ( sp_gas_const_a * r_T )

    r_rho_m = ( 1.D0 - SUM(r_alphas) ) * r_rho_a + SUM( r_alphas * rho_s ) 

    r_red_grav = ( r_rho_m - rho_a_amb ) / r_rho_m * grav

    IF ( r_h .GT. 0.D0 ) THEN

       r_u = qpj(2) / r_h
       r_v = qpj(3) / r_h

    ELSE

       r_u = 0.D0
       r_v = 0.D0

    END IF

    IF ( ( r_u**2 + r_v**2 ) .GT. 0.D0 ) THEN
    
       r_Ri = r_red_grav * r_h / ( r_u**2 + r_v**2 )

    ELSE

       r_Ri = 0.D0

    END IF

    RETURN

  END SUBROUTINE mixt_var

  !******************************************************************************
  !> \brief Conservative to physical variables
  !
  !> This subroutine evaluates from the conservative variables qc the 
  !> array of physical variables qp:\n
  !> - qp(1) = \f$ h+B \f$
  !> - qp(2) = \f$ u \f$
  !> - qp(3) = \f$ v \f$
  !> - qp(4) = \f$ T \f$
  !> - qp(5) = \f$ alphas \f$
  !> .
  !> The physical variables are those used for the linear reconstruction at the
  !> cell interfaces. Limiters are applied to the reconstructed slopes.
  !> \param[in]     qc      conservative variables 
  !> \param[out]    qp      physical variables  
  !> \date 15/08/2011
  !******************************************************************************
  
  SUBROUTINE qc_to_qp(qc,B,qp)
    
    IMPLICIT none
    
    REAL*8, INTENT(IN) :: qc(n_vars)
    REAL*8, INTENT(IN) :: B
    REAL*8, INTENT(OUT) :: qp(n_vars)
    
    CALL phys_var(r_qj = qc)
    
    ! It is important to have as reconstructed variable at the interfaces
    ! the total height h+B, instead of h. Only in this way, when the topography
    ! is not flat and the free surface is horizontal (steady condition), the  
    ! latter is kept flat by the reconstruction and the solution is stable.

    !qp(1) = DBLE(h+B)
    qp(1) = DBLE(h)
    
    qp(2) = DBLE(h*u)
    qp(3) = DBLE(h*v)

    
    qp(4) = DBLE(T)
    qp(5:4+n_solid) = DBLE(alphas(1:n_solid))

    RETURN
    
  END SUBROUTINE qc_to_qp

  !******************************************************************************
  !> \brief Physical to conservative variables
  !
  !> This subroutine evaluates the conservative variables qc from the 
  !> array of physical variables qp:\n
  !> - qp(1) = \f$ h + B \f$
  !> - qp(2) = \f$ rho_m*u \f$
  !> - qp(3) = \f$ rho_m*v \f$
  !> - qp(4) = \f$ T \f$
  !> - qp(5:4+n_s) = \f$ alphas(1:n_s) \f$
  !> .
  !> \param[in]    qp      physical variables  
  !> \param[out]   qc      conservative variables 
  !> \date 15/08/2011
  !******************************************************************************

  SUBROUTINE qp_to_qc(qp,B,qc)
    
    USE COMPLEXIFY 
    IMPLICIT none
    
    REAL*8, INTENT(IN) :: qp(n_vars)
    REAL*8, INTENT(IN) :: B
    REAL*8, INTENT(OUT) :: qc(n_vars)
    
    REAL*8 :: r_mass_fract_s(n_solid)
    REAL*8 :: r_sp_heat_mix

    !r_h = qp(1) - B
    r_h = qp(1)

    IF ( r_h .GT. 0.D0 ) THEN

       r_u = qp(2) / r_h
       r_v = qp(3) / r_h

    ELSE

       r_u = 0.D0
       r_v = 0.D0

    END IF

    r_T  = qp(4)
    r_alphas(1:n_solid) = qp(5:4+n_solid)
    r_rho_a = pres / ( sp_gas_const_a * r_T )

    r_rho_m = SUM( r_alphas * rho_s ) + ( 1.D0 - SUM(r_alphas) ) * r_rho_a

    r_mass_fract_s = r_alphas * rho_s / r_rho_m

    r_sp_heat_mix =  SUM( r_mass_fract_s * sp_heat_s ) + ( 1.D0 -               &
         SUM(r_mass_fract_s) ) * sp_heat_a

    qc(1) = r_h * r_rho_m
    qc(2) = r_h * r_rho_m * r_u
    qc(3) = r_h * r_rho_m * r_v

    IF ( energy_flag ) THEN

       qc(4) = r_h * r_rho_m * ( r_sp_heat_mix * r_T                            &
            + 0.5D0 * ( r_u**2 + r_v**2 ) ) 

    ELSE

       qc(4) = r_h * r_rho_m * r_sp_heat_mix * r_T 

    END IF

    qc(5:4+n_solid) = r_h * r_alphas(1:n_solid) * rho_s(1:n_solid)
       
  END SUBROUTINE qp_to_qc

  !******************************************************************************
  !> \brief Local Characteristic speeds x direction
  !
  !> This subroutine desingularize the velocities and then evaluates the largest 
  !> positive and negative characteristic speed in the x-direction. 
  !> \param[in]     qj            array of conservative variables
  !> \param[in]     Bj            topography at the cell center
  !> \param[out]    vel_min       minimum x-velocity
  !> \param[out]    vel_max       maximum x-velocity
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 05/12/2017
  !******************************************************************************

  SUBROUTINE eval_local_speeds_x(qj,vel_min,vel_max)
  
    IMPLICIT none
    
    REAL*8, INTENT(IN) :: qj(n_vars)

    REAL*8, INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)

    CALL phys_var(r_qj = qj)
    
    IF ( DBLE(red_grav * h) .LT. 0.D0 ) THEN

       vel_min(1:n_eqns) = DBLE(u)
       vel_max(1:n_eqns) = DBLE(u)
  
    ELSE

       vel_min(1:n_eqns) = DBLE(u) - DSQRT( DBLE(red_grav * h) )
       vel_max(1:n_eqns) = DBLE(u) + DSQRT( DBLE(red_grav * h) )
    
    END IF

  END SUBROUTINE eval_local_speeds_x

  !******************************************************************************
  !> \brief Local Characteristic speeds y direction
  !
  !> This subroutine desingularize the velocities and then evaluates the largest 
  !> positive and negative characteristic speed in the y-direction. 
  !> \param[in]     qj            array of conservative variables
  !> \param[in]     Bj            topography at the cell center
  !> \param[out]    vel_min       minimum y-velocity
  !> \param[out]    vel_max       maximum y-velocity
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 05/12/2017
  !******************************************************************************

  SUBROUTINE eval_local_speeds_y(qj,vel_min,vel_max)
    
    IMPLICIT none
    
    REAL*8, INTENT(IN)  :: qj(n_vars)
    REAL*8, INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)
    
    CALL phys_var(r_qj = qj)
 
    IF ( DBLE(red_grav * h) .LT. 0.D0 ) THEN

       vel_min(1:n_eqns) = DBLE(v)
       vel_max(1:n_eqns) = DBLE(v)
  
    ELSE
   
       vel_min(1:n_eqns) = DBLE(v) - DSQRT( DBLE(red_grav * h) )
       vel_max(1:n_eqns) = DBLE(v) + DSQRT( DBLE(red_grav * h) )
       
    END IF

  END SUBROUTINE eval_local_speeds_y

  !******************************************************************************
  !> \brief Local Characteristic speeds x direction
  !
  !> This subroutine computes from the physical variable evaluates the largest 
  !> positive and negative characteristic speed in the x-direction. 
  !> \param[in]     qpj           array of physical variables
  !> \param[out]    vel_min       minimum x-velocity
  !> \param[out]    vel_max       maximum x-velocity
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 05/12/2017
  !******************************************************************************

  SUBROUTINE eval_local_speeds2_x(qpj,vel_min,vel_max)
  
    IMPLICIT none
    
    REAL*8, INTENT(IN) :: qpj(n_vars)

    REAL*8, INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)

    CALL mixt_var(qpj)
    
    vel_min(1:n_eqns) = r_u - DSQRT( r_red_grav * r_h )
    vel_max(1:n_eqns) = r_u + DSQRT( r_red_grav * r_h )
    
  END SUBROUTINE eval_local_speeds2_x

  !******************************************************************************
  !> \brief Local Characteristic speeds y direction
  !
  !> This subroutine computes from the physical variable evaluates the largest 
  !> positive and negative characteristic speed in the y-direction. 
  !> \param[in]     qpj           array of physical variables
  !> \param[out]    vel_min       minimum y-velocity
  !> \param[out]    vel_max       maximum y-velocity
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 05/12/2017
  !******************************************************************************

  SUBROUTINE eval_local_speeds2_y(qpj,vel_min,vel_max)
    
    IMPLICIT none
    
    REAL*8, INTENT(IN)  :: qpj(n_vars)
    REAL*8, INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)
    
    CALL mixt_var(qpj)
    
    vel_min(1:n_eqns) = r_v - DSQRT( r_red_grav * r_h )
    vel_max(1:n_eqns) = r_v + DSQRT( r_red_grav * r_h )
       
  END SUBROUTINE eval_local_speeds2_y
  
  !******************************************************************************
  !> \brief Hyperbolic Fluxes
  !
  !> This subroutine evaluates the numerical fluxes given the conservative
  !> variables qcj and physical variables qpj.
  !> \date 01/06/2012
  !> \param[in]     qcj      real conservative variables 
  !> \param[in]     qpj      real physical variables 
  !> \param[in]     Bj       topography
  !> \param[in]     dir      direction of the flux (1=x,2=y)
  !> \param[out]    flux     real  fluxes    
  !******************************************************************************

  SUBROUTINE eval_fluxes(qcj,qpj,Bj,dir,flux)
    
    USE COMPLEXIFY
    IMPLICIT none

    REAL*8, INTENT(IN) :: qcj(n_vars)
    REAL*8, INTENT(IN) :: qpj(n_vars)
    REAL*8, INTENT(IN) :: Bj
    INTEGER, INTENT(IN) :: dir

    REAL*8, INTENT(OUT) :: flux(n_eqns)

    CALL mixt_var(qpj)

    
    ! pos_thick:IF ( r_h .NE. 0.D0 ) THEN
    pos_thick:IF ( qcj(1) .NE. 0.D0 ) THEN

       IF ( dir .EQ. 1 ) THEN

          ! Volumetric flux in x-direction: rhom*h*u
          flux(1) = r_u * qcj(1)
          !flux(1) = qcj(2) / qcj(1) * qcj(1)
          
          ! x-momentum flux in x-direction + hydrostatic pressure term
          flux(2) = r_u * qcj(2) + 0.5D0 * r_rho_m * r_red_grav * r_h**2

          ! y-momentum flux in x-direction: rho*h*v*u
          flux(3) = r_u * qcj(3)

          IF ( energy_flag ) THEN

             ! ENERGY flux in x-direction
             flux(4) = r_u * ( qcj(4) + 0.5D0 * r_rho_m * r_red_grav * r_h**2 )

          ELSE

             ! Temperature flux in x-direction
             flux(4) = r_u * qcj(4)

          END IF

          ! Volumetric flux of solid in x-direction: h * alphas * u
          flux(5:4+n_solid) = r_u * qcj(5:4+n_solid)

          ! Solid flux can't be larger than total flux
          IF ( ( flux(1) .GT. 0.D0 ) .AND. ( SUM(flux(5:4+n_solid)) / flux(1)   &
               .GT. 1.D0 ) ) THEN

             flux(5:4+n_solid) = &
               flux(5:4+n_solid) / SUM(flux(5:4+n_solid)) * flux(1)

          END IF

       ELSEIF ( dir .EQ. 2 ) THEN

          ! flux G (derivated wrt y in the equations)
          flux(1) = r_v * qcj(1)
          
          flux(2) = r_v * qcj(2)
          
          flux(3) = r_v * qcj(3) + 0.5D0 * r_rho_m * r_red_grav * r_h**2

          IF ( energy_flag ) THEN

             ! ENERGY flux in x-direction
             flux(4) = r_v * ( qcj(4) + 0.5D0 * r_rho_m * r_red_grav * r_h**2 )

          ELSE

             ! Temperature flux in y-direction
             flux(4) = r_v * qcj(4)
          
          END IF

          ! Volumetric flux of solid in y-direction: h * alphas * v
          flux(5:4+n_solid) = r_v * qcj(5:4+n_solid)
          
          ! Solid flux can't be larger than total flux
          IF ( ( flux(1) .GT. 0.D0 ) .AND. ( SUM(flux(5:4+n_solid)) / flux(1)   &
               .GT. 1.D0 ) ) THEN

             flux(5:4+n_solid) = &
               flux(5:4+n_solid) / SUM(flux(5:4+n_solid)) * flux(1)

          END IF

                     
       END IF

    ELSE
       
       flux(1:n_eqns) = 0.D0
       
    ENDIF pos_thick
 
  END SUBROUTINE eval_fluxes


  !******************************************************************************
  !> \brief Non-Hyperbolic terms
  !
  !> This subroutine evaluates the non-hyperbolic terms (relaxation terms
  !> and forces) of the system of equations, both for real or complex 
  !> inputs. These terms are treated implicitely in the DIRK numerical
  !> scheme.
  !> \date 01/06/2012
  !> \param[in]     c_qj            complex conservative variables 
  !> \param[in]     r_qj            real conservative variables 
  !> \param[out]    c_nh_term_impl  complex non-hyperbolic terms     
  !> \param[out]    r_nh_term_impl  real non-hyperbolic terms
  !******************************************************************************

  SUBROUTINE eval_nonhyperbolic_terms( c_qj , c_nh_term_impl , r_qj ,           &
       r_nh_term_impl )

    USE COMPLEXIFY 

    IMPLICIT NONE

    COMPLEX*16, INTENT(IN), OPTIONAL :: c_qj(n_vars)

    COMPLEX*16, INTENT(OUT), OPTIONAL :: c_nh_term_impl(n_eqns)
    REAL*8, INTENT(IN), OPTIONAL :: r_qj(n_vars)
    REAL*8, INTENT(OUT), OPTIONAL :: r_nh_term_impl(n_eqns)

    COMPLEX*16 :: qj(n_vars)

    COMPLEX*16 :: nh_term(n_eqns)

    COMPLEX*16 :: relaxation_term(n_eqns)

    COMPLEX*16 :: forces_term(n_eqns)

    INTEGER :: i

    COMPLEX*16 :: mod_vel
    
    COMPLEX*16 :: gamma

    REAL*8 :: h_threshold

    !--- Lahars rheology model variables
    
    !> Fluid viscosity
    COMPLEX*16 :: fluid_visc

    !> Total friction
    COMPLEX*16 :: s_f

    !> Viscous slope component of total Friction
    COMPLEX*16 :: s_v

    !> Turbulent dispersive slope component of total friction
    COMPLEX*16 :: s_td

    IF ( ( thermal_conductivity .GT. 0.D0 ) .OR. ( emme .GT. 0.D0 ) ) THEN
       
       h_threshold = 1.D-10
       
    ELSE
       
       h_threshold = 0.D0
       
    END IF
    

    IF ( present(c_qj) .AND. present(c_nh_term_impl) ) THEN

       qj = c_qj

    ELSEIF ( present(r_qj) .AND. present(r_nh_term_impl) ) THEN

       DO i = 1,n_vars

          qj(i) = DCMPLX( r_qj(i) )

       END DO

    ELSE

       WRITE(*,*) 'Constitutive, eval_fluxes: problem with arguments'
       STOP

    END IF

    ! initialize and evaluate the relaxation terms
    relaxation_term(1:n_eqns) = DCMPLX(0.D0,0.D0) 

    ! initialize and evaluate the forces terms
    forces_term(1:n_eqns) = DCMPLX(0.D0,0.D0)

    IF (rheology_flag) THEN

       CALL phys_var(c_qj = qj)
    
       mod_vel = CDSQRT( u**2 + v**2 )
       
       ! Voellmy Salm rheology
       IF ( rheology_model .EQ. 1 ) THEN
       
          IF ( DBLE(mod_vel) .NE. 0.D0 ) THEN 
          
             ! IMPORTANT: grav3_surv is always negative 
             forces_term(2) = forces_term(2) - rho_m * ( u / mod_vel ) *        &
                  ( grav / xi ) * mod_vel ** 2
          
             forces_term(3) = forces_term(3) - rho_m * ( v / mod_vel ) *        &
                  ( grav / xi ) * mod_vel ** 2
          
          ENDIF
        
       ! Plastic rheology
       ELSEIF ( rheology_model .EQ. 2 ) THEN
       
          IF ( DBLE(mod_vel) .NE. 0.D0 ) THEN 
          
             forces_term(2) = forces_term(2) - rho_m * tau * (u/mod_vel)
          
             forces_term(3) = forces_term(3) - rho_m * tau * (v/mod_vel)

          ENDIF

       ! Temperature dependent rheology
       ELSEIF ( rheology_model .EQ. 3 ) THEN

          IF ( DBLE(h) .GT. h_threshold ) THEN
    
             ! Equation 6 from Costa & Macedonio, 2005
             gamma = 3.D0 * nu_ref / h * CDEXP( - visc_par * ( T - T_ref ) )

          ELSE

             ! Equation 6 from Costa & Macedonio, 2005
             gamma = 3.D0 * nu_ref / h_threshold * CDEXP( - visc_par            &
                  * ( T - T_ref ) )
             
          END IF
          
          IF ( DBLE(mod_vel) .NE. 0.D0 ) THEN 
          
             ! Last R.H.S. term in equation 2 from Costa & Macedonio, 2005
             forces_term(2) = forces_term(2) - rho_m * gamma * u
          
             ! Last R.H.S. term in equation 3 from Costa & Macedonio, 2005
             forces_term(3) = forces_term(3) - rho_m * gamma * v

          ENDIF
          
       ! Lahars rheology (O'Brien 1993, FLO2D)
       ELSEIF ( rheology_model .EQ. 4 ) THEN

          h_threshold = 1.D-20

          ! THIS IS ONLY USED WHEN SEDIMENT PERCENTAGE IS GIVEN AS INPUT AND
          ! IT IS CONSTANT WITHIN THE SIMULATION
          ! TO BE REMOVED
          ! alphas = DCMPLX(sed_vol_perc/100.D0,0.D0)

          ! Fluid viscosity
          fluid_visc = alpha1 * CDEXP( beta1 * SUM(alphas) )

          IF ( h .GT. h_threshold ) THEN
             
             ! Viscous slope component
             s_v = Kappa * fluid_visc * mod_vel / ( 8.D0 * h**2 )
             
             ! Turbulent dispersive component
             s_td = rho_m * n_td**2 * mod_vel**2 / ( h**(4.D0/3.D0) )
          
          ELSE
             
             ! Viscous slope component
             s_v = Kappa * fluid_visc * mod_vel / ( 8.D0 * h_threshold**2 )
             
             ! Turbulent dispersive components
             s_td = rho_m * n_td**2 * (mod_vel**2) / ( h_threshold**(4.D0/3.D0) )
             
          END IF
          
          ! Total implicit friction slope
          s_f = s_v + s_td
         
          IF ( mod_vel .GT. 0.D0 ) THEN

             forces_term(2) = forces_term(2) - grav * h * ( u / mod_vel ) * s_f

             forces_term(3) = forces_term(3) - grav * h * ( v / mod_vel ) * s_f
  
          END IF

       ELSEIF ( rheology_model .EQ. 5 ) THEN

          tau = 1.D-3 / ( 1.D0 + 10.D0 * h ) * mod_vel
          
          IF ( DBLE(mod_vel) .NE. 0.D0 ) THEN
             
             forces_term(2) = forces_term(2) - rho_m * tau * ( u / mod_vel )
             forces_term(3) = forces_term(3) - rho_m * tau * ( v / mod_vel )

          END IF
          

       ELSEIF ( rheology_model .EQ. 6 ) THEN
          
          IF ( DBLE(mod_vel) .NE. 0.D0 ) THEN 
             
             forces_term(2) = forces_term(2) - rho_m * ( u / mod_vel ) *        &
                  friction_factor * mod_vel ** 2
             
             forces_term(3) = forces_term(3) - rho_m * ( v / mod_vel ) *        &
                  friction_factor * mod_vel ** 2
             
          ENDIF
          
       ENDIF
       
    ENDIF

    nh_term = relaxation_term + forces_term

    IF ( present(c_qj) .AND. present(c_nh_term_impl) ) THEN

       c_nh_term_impl = nh_term

    ELSEIF ( present(r_qj) .AND. present(r_nh_term_impl) ) THEN

       r_nh_term_impl = DBLE( nh_term )

    END IF 

  END SUBROUTINE eval_nonhyperbolic_terms

  !******************************************************************************
  !> \brief Non-Hyperbolic semi-implicit terms
  !
  !> This subroutine evaluates the non-hyperbolic terms that are solved
  !> semi-implicitely by the solver. For example, any discontinuous term that
  !> appears in the friction terms.
  !> \date 20/01/2018
  !> \param[in]     c_qj            complex conservative variables 
  !> \param[in]     r_qj            real conservative variables 
  !> \param[out]    c_nh_term_impl  complex non-hyperbolic terms     
  !> \param[out]    r_nh_term_impl  real non-hyperbolic terms
  !******************************************************************************

  SUBROUTINE eval_nh_semi_impl_terms( grav3_surf , c_qj , c_nh_semi_impl_term , &
       r_qj , r_nh_semi_impl_term )

    USE COMPLEXIFY 


    IMPLICIT NONE

    REAL*8, INTENT(IN) :: grav3_surf

    COMPLEX*16, INTENT(IN), OPTIONAL :: c_qj(n_vars)
    COMPLEX*16, INTENT(OUT), OPTIONAL :: c_nh_semi_impl_term(n_eqns)
    REAL*8, INTENT(IN), OPTIONAL :: r_qj(n_vars)
    REAL*8, INTENT(OUT), OPTIONAL :: r_nh_semi_impl_term(n_eqns)

    COMPLEX*16 :: qj(n_vars)

    COMPLEX*16 :: forces_term(n_eqns)

    INTEGER :: i

    COMPLEX*16 :: mod_vel
    
    REAL*8 :: h_threshold

    !--- Lahars rheology model variables
    
    !> Yield strenght
    COMPLEX*8 :: tau_y

    !> Yield slope component of total friction
    COMPLEX*8 :: s_y


    IF ( present(c_qj) .AND. present(c_nh_semi_impl_term) ) THEN

       qj = c_qj

    ELSEIF ( present(r_qj) .AND. present(r_nh_semi_impl_term) ) THEN

       DO i = 1,n_vars

          qj(i) = DCMPLX( r_qj(i) )

       END DO

    ELSE

       WRITE(*,*) 'Constitutive, eval_fluxes: problem with arguments'
       STOP

    END IF

    ! initialize and evaluate the forces terms
    forces_term(1:n_eqns) = DCMPLX(0.D0,0.D0)

    IF (rheology_flag) THEN

       CALL phys_var(c_qj = qj)
    
       mod_vel = CDSQRT( u**2 + v**2 )

       ! Voellmy Salm rheology
       IF ( rheology_model .EQ. 1 ) THEN

          IF ( mod_vel .GT. 0.D0 ) THEN

             forces_term(2) = forces_term(2) - rho_m * ( u / mod_vel ) *        &
                  mu * h * ( - grav * grav3_surf )
             
             forces_term(3) = forces_term(3) - rho_m * ( v / mod_vel ) *        &
                  mu * h * ( - grav * grav3_surf )
             
          END IF

          ! Plastic rheology
       ELSEIF ( rheology_model .EQ. 2 ) THEN
          

       ! Temperature dependent rheology
       ELSEIF ( rheology_model .EQ. 3 ) THEN

          
       ! Lahars rheology (O'Brien 1993, FLO2D)
       ELSEIF ( rheology_model .EQ. 4 ) THEN

          h_threshold = 1.D-20

          ! THIS IS ONLY USED WHEN SEDIMENT PERCENTAGE IS GIVEN AS INPUT AND
          ! IT IS CONSTANT WITHIN THE SIMULATION
          ! alphas = DCMPLX( sed_vol_perc*1.D-2 , 0.D0 )

          ! Yield strength
          tau_y = alpha2 * CDEXP( beta2 * SUM(alphas) )

          IF ( h .GT. h_threshold ) THEN
             
             ! Yield slope component
             s_y = tau_y / h
                       
          ELSE
             
             ! Yield slope component
             s_y = tau_y / h_threshold
                          
          END IF
  
          IF ( mod_vel .GT. 0.D0 ) THEN

             forces_term(2) = forces_term(2) - grav * h * ( u / mod_vel ) * s_y

             forces_term(3) = forces_term(3) - grav * h * ( v / mod_vel ) * s_y
  
          END IF

       ELSEIF ( rheology_model .EQ. 5 ) THEN

          
       ENDIF
              
    ENDIF

    
    IF ( present(c_qj) .AND. present(c_nh_semi_impl_term) ) THEN

       c_nh_semi_impl_term = forces_term

    ELSEIF ( present(r_qj) .AND. present(r_nh_semi_impl_term) ) THEN

       r_nh_semi_impl_term = DBLE( forces_term )

    END IF 

  END SUBROUTINE eval_nh_semi_impl_terms

  
  !******************************************************************************
  !> \brief Explicit Forces term
  !
  !> This subroutine evaluates the non-hyperbolic terms to be treated explicitely
  !> in the DIRK numerical scheme (e.g. gravity,source of mass). The sign of the
  !> terms is taken with the terms on the left-hand side of the equations.
  !> \date 01/06/2012
  !> \param[in]     qj                 conservative variables 
  !> \param[out]    expl_term          explicit term
  !******************************************************************************

  SUBROUTINE eval_expl_terms( Bprimej_x , Bprimej_y , source_xy , qj ,          &
       expl_term )

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: Bprimej_x
    REAL*8, INTENT(IN) :: Bprimej_y
    REAL*8, INTENT(IN) :: source_xy
    
    REAL*8, INTENT(IN) :: qj(n_eqns)                 !< conservative variables 
    REAL*8, INTENT(OUT) :: expl_term(n_eqns)         !< explicit forces 

    expl_term(1:n_eqns) = 0.D0

    CALL phys_var(r_qj = qj)

    expl_term(2) = DBLE(red_grav*rho_m*h) * Bprimej_x
   
    expl_term(3) = DBLE(red_grav*rho_m*h) * Bprimej_y

    IF ( energy_flag ) THEN

       expl_term(4) = DBLE(red_grav*rho_m*h) * ( DBLE(u) * Bprimej_x            &
            + DBLE(v) * Bprimej_y )  

    ELSE

       expl_term(4) = 0.D0

    END IF

    RETURN

  END SUBROUTINE eval_expl_terms

  SUBROUTINE eval_expl_terms2( Bprimej_x , Bprimej_y , source_xy , qpj , qcj ,  &
       expl_term )

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: Bprimej_x
    REAL*8, INTENT(IN) :: Bprimej_y
    REAL*8, INTENT(IN) :: source_xy
    
    REAL*8, INTENT(IN) :: qpj(n_eqns)        !< local physical variables 
    REAL*8, INTENT(IN) :: qcj(n_eqns)        !< local conservative variables 
    REAL*8, INTENT(OUT) :: expl_term(n_eqns) !< local explicit forces 

    CALL mixt_var(qpj)
    
    expl_term(1:n_eqns) = 0.D0

    expl_term(2) = r_red_grav * r_rho_m * r_h * Bprimej_x
   
    expl_term(3) = r_red_grav * r_rho_m * r_h * Bprimej_y

    IF ( energy_flag ) THEN

       expl_term(4) = r_red_grav * r_rho_m * r_h * ( r_u * Bprimej_x            &
            + r_v * Bprimej_y )  

    ELSE

       expl_term(4) = 0.D0

    END IF

    RETURN

  END SUBROUTINE eval_expl_terms2


  !******************************************************************************
  !> \brief Deposition term
  !
  !> This subroutine evaluates the deposition term.
  !> \date 03/010/2018
  !> \param[in]     Bj                 topography
  !> \param[in]     qj                 conservative variables 
  !> \param[out]    dep_term          explicit term
  !******************************************************************************

  SUBROUTINE eval_erosion_dep_term( qpj , Bj , dt , erosion_term ,              &
       deposition_term )
    
    IMPLICIT NONE
    
    REAL*8, INTENT(IN) :: qpj(n_eqns)                !< physical variables 
    REAL*8, INTENT(IN) :: Bj
    REAL*8, INTENT(IN) :: dt

    REAL*8, INTENT(OUT) :: erosion_term(n_solid)     !< erosion term
    REAL*8, INTENT(OUT) :: deposition_term(n_solid)  !< deposition term
    
    REAL*8 :: mod_vel

    INTEGER :: i_solid

    deposition_term(1:n_solid) = 0.D0
    erosion_term(1:n_solid) = 0.D0
 
    IF ( qpj(1) .LE. 0.D0 ) RETURN

    CALL mixt_var(qpj)
    
    deposition_term(1:n_solid) = 0.D0
    erosion_term(1:n_solid) = 0.D0
    
    DO i_solid=1,n_solid

       IF ( r_alphas(i_solid) .GT. 0.D0 ) THEN

          settling_vel = settling_velocity( diam_s(i_solid) , rho_s(i_solid) ,  &
               r_rho_a , i_solid )
         
          deposition_term(i_solid) = r_alphas(i_solid) * settling_vel
          
          deposition_term(i_solid) = MIN( deposition_term(i_solid) ,            &
               r_h * r_alphas(i_solid) / dt )

          IF ( deposition_term(i_solid) .LT. 0.D0 ) THEN

             WRITE(*,*) 'eval_erosion_dep_term'
             WRITE(*,*) 'deposition_term(i_solid)',deposition_term(i_solid)
             READ(*,*)

          END IF
          
       END IF

       mod_vel = DSQRT( r_u**2 + r_v**2 )
  
       IF ( DBLE(r_h) .GT. 1.D-2) THEN
    
          erosion_term(i_solid) = 0.D0

       ELSE
          
          erosion_term(i_solid) = 0.D0
          
       END IF
           
    END DO

    RETURN
  
  END SUBROUTINE eval_erosion_dep_term


  !******************************************************************************
  !> \brief Topography modification related term
  !
  !> This subroutine evaluates the deposition term.
  !> \date 03/010/2018
  !> \param[in]     qj                  conservative variables 
  !> \param[in]     deposition_avg_term averaged deposition terms 
  !> \param[in]     erosion_avg_term    averaged deposition terms 
  !> \param[out]    topo_term           explicit term
  !******************************************************************************

  SUBROUTINE eval_topo_term( qj , deposition_avg_term , erosion_avg_term ,      &
       eqns_term, topo_term )
    
    USE parameters_2d, ONLY : topo_change_flag

    IMPLICIT NONE
    
    REAL*8, INTENT(IN) :: qj(n_eqns)                   !< conservative variables 
    REAL*8, INTENT(IN) :: deposition_avg_term(n_solid) !< deposition term
    REAL*8, INTENT(IN) :: erosion_avg_term(n_solid)    !< erosion term

    REAL*8, INTENT(OUT):: eqns_term(n_eqns)
    REAL*8, INTENT(OUT):: topo_term
   
    REAL*8 :: entr_coeff
    REAL*8 :: air_entr
    REAL*8 :: mag_vel 


    CALL phys_var(r_qj = qj)

    mag_vel = DSQRT( DBLE(u)**2.D0 + DBLE(v)**2.D0 ) 
    
    IF ( entrainment_flag .AND. ( mag_vel**2 .GT. 0.D0 ) .AND.                  &
         ( DBLE(h) .GT. 0.D0 ) ) THEN

       entr_coeff = 0.075D0 / DSQRT( 1.D0 + 718.D0 * MAX(0.D0,r_Ri)**2.4 )
       
       air_entr = entr_coeff * mag_vel
       
    ELSE

       air_entr = 0.D0

    END IF

    !air_entr = 0.D0

    eqns_term(1:n_eqns) = 0.D0

    ! free surface (topography+flow) equation
    eqns_term(1) = SUM( rho_s * ( erosion_avg_term - deposition_avg_term ) ) +  &
         rho_a_amb * air_entr

    ! x-momenutm equation
    eqns_term(2) = - DBLE(u) * SUM( rho_s * deposition_avg_term )

    ! y-momentum equation
    eqns_term(3) = - DBLE(v) * SUM( rho_s * deposition_avg_term )

    ! Temperature/Energy equation
    IF ( energy_flag ) THEN

       eqns_term(4) = - DBLE(T) * SUM( rho_s * sp_heat_s * deposition_avg_term )&
            - 0.5D0 * mag_vel**2 * SUM( rho_s * deposition_avg_term )           &
            + T_s_substrate * SUM( rho_s * sp_heat_s * erosion_avg_term )       &
            + T_ambient * sp_heat_a * rho_a_amb * air_entr

    ELSE

       eqns_term(4) = - DBLE(T) * SUM( rho_s * sp_heat_s * deposition_avg_term )&
            + T_s_substrate * SUM( rho_s * sp_heat_s * erosion_avg_term )       &
            + T_ambient * sp_heat_a * rho_a_amb * air_entr

    END IF

   
    ! solid phase thickness equation
    eqns_term(5:4+n_solid) = rho_s * ( erosion_avg_term(1:n_solid)              &
         - deposition_avg_term(1:n_solid) )

    IF ( topo_change_flag ) THEN

       ! topography term
       topo_term = SUM( deposition_avg_term - erosion_avg_term )

    ELSE

       topo_term = 0.D0

    END IF

  END SUBROUTINE eval_topo_term

  
  !------------------------------------------------------------------------------
  !> Settling velocity function
  !
  !> This subroutine compute the settling velocity of the particles, as a
  !> function of diameter, density.
  !> \date OCTOBER 2016
  !> \param    diam        particle diameter                (\b input)
  !> \param    rhos        particle density                 (\b input)
  !> \param    rhoa        atmospheric density              (\b input)
  !> \patam    i_solid     particle class index             (\b input)
  !------------------------------------------------------------------------------
  
  REAL*8 FUNCTION settling_velocity(diam,rhos,rhoa,i_solid)

    IMPLICIT NONE

    REAL*8 :: diam
    REAL*8 :: rhos
    REAL*8 :: rhoa

    REAL*8 :: const_part

    REAL*8 :: Rey
    REAL*8 :: C_D , C_D_old

    REAL*8 :: set_vel_old

    INTEGER :: i_solid
    INTEGER :: i

    INTEGER :: dig

    C_D = 1.D0

    const_part =  DSQRT( 0.75D0 * ( rhos / rhoa - 1.D0 ) * diam * grav )
  
    settling_velocity = const_part / DSQRT( C_D )

    Rey = diam * settling_velocity / kin_visc_a

    IF ( Rey .LE. 1000.D0 ) THEN

       C_D_loop:DO i=1,20

          set_vel_old = settling_velocity
          C_D_old = C_D
          C_D = 24.D0 / Rey * ( 1.D0 + 0.15D0 * Rey**(0.687) )
          
          settling_velocity = const_part / DSQRT( C_D )
          
          IF ( DABS( set_vel_old - settling_velocity ) / set_vel_old            &
               .LT. 1.D-6 ) THEN

             ! round to first three significative digits
             dig = FLOOR(log10(set_vel_old))
             settling_velocity = 10.D0**(dig-3)*FLOOR( 10.0**(-dig+3)*set_vel_old ) 

             EXIT C_D_loop

          END IF

          Rey = diam * settling_velocity / kin_visc_a

       END DO C_D_loop
    
    END IF

    C_D_s(i_solid) = C_D
    
    ! WRITE(*,*) 'C_D,settling_velocity',C_D,settling_velocity
    ! READ(*,*)
   
    RETURN

  END FUNCTION settling_velocity


END MODULE constitutive_2d

    
