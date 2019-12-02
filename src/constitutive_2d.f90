!********************************************************************************
!> \brief Constitutive equations
!********************************************************************************
MODULE constitutive_2d

  USE parameters_2d, ONLY : n_eqns , n_vars , n_solid
  USE parameters_2d, ONLY : rheology_flag , rheology_model , energy_flag ,      &
       liquid_flag , gas_flag
  USE parameters_2d, ONLY : sed_vol_perc

  IMPLICIT none

  LOGICAL, ALLOCATABLE :: implicit_flag(:)

  CHARACTER(LEN=20) :: phase1_name
  CHARACTER(LEN=20) :: phase2_name

  COMPLEX*16 :: h                       !< height [m]
  COMPLEX*16 :: u                       !< velocity (x direction) [m/s]
  COMPLEX*16 :: v                       !< velocity (y direction) [m/s]
  COMPLEX*16 :: T                       !< temperature [K]
  COMPLEX*16, ALLOCATABLE :: alphas(:)  !< sediment volume fractions
  COMPLEX*16 :: rho_m                   !< mixture density [kg/m3]
  COMPLEX*16 :: Ri                      !< Richardson number
  COMPLEX*16 :: alphal                  !< liquid volume fraction
  COMPLEX*16 :: alphac                  !< carrier phase volume fraction
  COMPLEX*16, ALLOCATABLE :: xs(:)      !< sediment mass fractions
  COMPLEX*16 :: xl                      !< liquid mass fraction
  COMPLEX*16 :: xc                      !< carrier phase mass fraction
  
  COMPLEX*16 :: xs_tot                  !< sum of solid mass fraction
  COMPLEX*16 :: rhos_tot                !< average density of solids [kg/m3]

  !> gravitational acceleration 
  REAL*8 :: grav

  !> reduced gravity
  COMPLEX*16 :: red_grav

  REAL*8 :: r_h          !< real-value flow thickness
  REAL*8 :: r_u          !< real-value x-velocity
  REAL*8 :: r_v          !< real-value y-velocity
  REAL*8 :: r_T          !< real-value temperature [K]
  REAL*8 :: r_red_grav   !< real-value reduced gravity

  REAL*8 :: r_rho_m      !< real-value mixture density [kg/m3]
  REAL*8 :: r_rho_c      !< real-value carrier phase density [kg/m3]

  REAL*8 :: r_Ri         !< real-value Richardson number
 
  REAL*8, ALLOCATABLE :: r_alphas(:) !< real-value solid volume fractions
  REAL*8 :: r_alphal                 !< real-value liquid volume fraction
  REAL*8 :: r_alphac                 !< real-value carrier phase volume fraction

  REAL*8, ALLOCATABLE :: r_xs(:)     !< real-value solid mass fraction
  REAL*8 :: r_xl                     !< real-value liquid mass fraction
  REAL*8 :: r_xc                     !< real-value carrier phase mass fraction

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

  !> 1st param for yield strenght empirical relationship (O'Brian et al, 1993)
  REAL*8 :: alpha2    ! (units: kg m-1 s-2)

  !> 2nd param for yield strenght empirical relationship (O'Brian et al, 1993)
  REAL*8 :: beta2     ! (units: nondimensional) 

  !> ratio between reference value from input and computed values from eq.
  REAL*8 :: alpha1_coeff ! (units: nondimensional )

  !> 1st param for fluid viscosity empirical relationship (O'Brian et al, 1993)
  COMPLEX*16 :: alpha1    ! (units: kg m-1 s-1 )

  !> 2nd param for fluid viscosity empirical relationship (O'Brian et al, 1993)
  REAL*8 :: beta1     ! (units: nondimensional)

  !> Empirical resistance parameter (dimensionless)
  REAL*8 :: Kappa

  !> Mannings roughness coefficient ( units: T L^(-1/3) )
  REAL*8 :: n_td
 
  !> Density of air in the mixture ( units: kg/m^3 )
  COMPLEX*16 :: rho_c

  !> Specific heat of carrier phase (gas or liquid)
  REAL*8 :: sp_heat_c
  
  !> Ambient density of air ( units: kg/m^3 )
  REAL*8 :: rho_a_amb

  !> Specific heat of air
  REAL*8 :: sp_heat_a

  !> Specific gas constant of air
  REAL*8 :: sp_gas_const_a

  !> Kinematic viscosity of air
  REAL*8 :: kin_visc_a

  !> Kinematic viscosity of liquid
  REAL*8 :: kin_visc_l

  !> Kinematic viscosity of carrier phase
  REAL*8 :: kin_visc_c
  
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

  LOGICAL :: settling_flag

  !> hindered settling 
  REAL*8 :: settling_vel

  !> erosion model coefficient  ( units: m-1 )
  REAL*8, ALLOCATABLE :: erosion_coeff(:)

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

    COMPLEX*16 :: inv_rhom
    
    IF ( present(c_qj) ) THEN

       qj = c_qj

    ELSE

       qj = DCMPLX(r_qj)

    END IF

    ! compute solid mass fractions
    IF ( DBLE(qj(1)) .GT. 1.D-25 ) THEN
    
       xs(1:n_solid) = qj(5:4+n_solid) / qj(1)

    ELSE

       xs(1:n_solid) = DCMPLX(0.D0,0.D0)

    END IF
    
    xs_tot = SUM(xs)
    
    IF ( gas_flag .AND. liquid_flag ) THEN

       ! compute liquid mass fraction
       IF ( DBLE(qj(1)) .GT. 1.D-25 ) THEN
          
          xl = qj(n_vars) / qj(1)
          
       ELSE
          
          xl = DCMPLX(0.D0,0.D0) 
       
       END IF

       ! compute carrier phase (gas) mass fraction
       xc =  DCMPLX(1.D0,0.D0) - xs_tot - xl

       ! specific heaf of the mixutre: mass average of sp. heat pf phases
       sp_heat_mix = SUM( xs(1:n_solid) * sp_heat_s(1:n_solid) ) + xl * sp_heat_l  &
            + xc * sp_heat_c
       
    ELSE

       ! compute carrier phase (gas or liquid) mass fraction
       xc = DCMPLX(1.D0,0.D0) - xs_tot

       ! specific heaf of the mixutre: mass average of sp. heat pf phases
       sp_heat_mix = SUM( xs(1:n_solid) * sp_heat_s(1:n_solid) ) + xc * sp_heat_c
       
    END IF

    ! compute temperature from energy
    IF ( DBLE(qj(1)) .GT. 1.D-25 ) THEN

       IF ( energy_flag ) THEN

          T = ( qj(4) - ( qj(2)**2 + qj(3)**2 ) / ( 2.D0*qj(1) ) ) /            &
               ( qj(1) * sp_heat_mix ) 
 
       ELSE

          T = qj(4) / ( qj(1) * sp_heat_mix ) 

       END IF
   
       IF ( DBLE(T) .LE. 0.D0 ) T = DCMPLX(T_ambient,0.D0)

    ELSE

       T = DCMPLX(T_ambient,0.D0)

    END IF

    IF ( gas_flag ) THEN

       ! carrier phase is gas
       rho_c =  pres / ( sp_gas_const_a * T )
       sp_heat_c = sp_heat_a

    ELSE

       ! carrier phase is liquid
       rho_c = rho_l
       sp_heat_c = sp_heat_l
       
    END IF
       
    IF ( gas_flag .AND. liquid_flag ) THEN

       inv_rhom = ( SUM(xs(1:n_solid) / rho_s(1:n_solid)) + xl / rho_l + xc / rho_c )

       rho_m = 1.D0 / inv_rhom

       alphal = xl * rho_m / rho_l
       
    ELSE

       inv_rhom = ( SUM(xs(1:n_solid) / rho_s(1:n_solid)) + xc / rho_c )

       rho_m = 1.D0 / inv_rhom
       
    END IF

    ! convert from mass fraction to volume fraction
    alphas(1:n_solid) = xs(1:n_solid) * rho_m / rho_s(1:n_solid)

    ! convert from mass fraction to volume fraction
    alphac = xc * rho_m / rho_c

    h = qj(1) / rho_m
    
    ! reduced gravity
    red_grav = ( rho_m - rho_a_amb ) / rho_m * grav

    ! velocity components
    IF ( DBLE( qj(1) ) .GT. eps_sing ) THEN
       
       u = qj(2) / qj(1)
       v = qj(3) / qj(1)

    ELSE

       u = DSQRT(2.D0) * qj(1) * qj(2) / CDSQRT( qj(1)**4 + eps_sing**4 )
       v = DSQRT(2.D0) * qj(1) * qj(3) / CDSQRT( qj(1)**4 + eps_sing**4 )
       
    END IF

    ! Richardson number
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
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 10/10/2019
  !******************************************************************************
 
  SUBROUTINE mixt_var(qpj)

    IMPLICIT none

    REAL*8, INTENT(IN) :: qpj(n_vars)

    r_h = qpj(1)

    IF ( r_h .LE. 0.D0 ) THEN

       r_u = 0.D0
       r_v = 0.D0
       r_T = T_ambient
       r_alphas(1:n_solid) = 0.D0
       r_red_grav = 0.D0
       r_Ri = 0.D0
       IF ( gas_flag .AND. liquid_flag ) r_alphal = 0.D0

       RETURN

    END IF

    r_u = qpj(2) / r_h
    r_v = qpj(3) / r_h

    r_T = qpj(4)

    r_alphas(1:n_solid) = qpj(5:4+n_solid) 

    IF ( gas_flag ) THEN

       ! continuous phase is air
       r_rho_c =  pres / ( sp_gas_const_a * r_T )

    ELSE

       ! continuous phase is liquid
       r_rho_c = rho_l

    END IF
       
    IF ( gas_flag .AND. liquid_flag ) THEN
  
       r_alphal = qpj(n_vars)

       ! density of mixture of carrier (gas), liquid and solids
       r_rho_m = ( 1.D0 - SUM(r_alphas) - r_alphal ) * r_rho_c                  &
            + SUM( r_alphas * rho_s ) + r_alphal * rho_l
       
    ELSE

       ! density of mixture of carrier phase and solids
       r_rho_m = ( 1.D0 - SUM(r_alphas) ) * r_rho_c + SUM( r_alphas * rho_s ) 

    END IF

    ! reduced gravity
    r_red_grav = ( r_rho_m - rho_a_amb ) / r_rho_m * grav

    ! Richardson number
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
  !> - qp(1) = \f$ h \f$
  !> - qp(2) = \f$ hu \f$
  !> - qp(3) = \f$ hv \f$
  !> - qp(4) = \f$ T \f$
  !> - qp(5:4+n_solid) = \f$ alphas(1:n_solid) \f$
  !> - qp(n_vars) = \f$ alphal \f$
  !> .
  !> The physical variables are those used for the linear reconstruction at the
  !> cell interfaces. Limiters are applied to the reconstructed slopes.
  !> \param[in]     qc     local conservative variables 
  !> \param[in]     B      local topography
  !> \param[out]    qp     local physical variables  
  !> \date 2019/11/11
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************
  
  SUBROUTINE qc_to_qp(qc,B,qp)
    
    IMPLICIT none
    
    REAL*8, INTENT(IN) :: qc(n_vars)
    REAL*8, INTENT(IN) :: B
    REAL*8, INTENT(OUT) :: qp(n_vars)
    
    CALL phys_var(r_qj = qc)
    
    qp(1) = DBLE(h)
    
    qp(2) = DBLE(h*u)
    qp(3) = DBLE(h*v)
    
    qp(4) = DBLE(T)
    qp(5:4+n_solid) = DBLE(alphas(1:n_solid))

    IF ( gas_flag .AND. liquid_flag ) qp(n_vars) = alphal

    RETURN
    
  END SUBROUTINE qc_to_qp

  !******************************************************************************
  !> \brief Physical to conservative variables
  !
  !> This subroutine evaluates the conservative variables qc from the 
  !> array of physical variables qp:\n
  !> - qp(1) = \f$ h \f$
  !> - qp(2) = \f$ h*u \f$
  !> - qp(3) = \f$ h*v \f$
  !> - qp(4) = \f$ T \f$
  !> - qp(5:4+n_s) = \f$ alphas(1:n_s) \f$
  !> - qp(n_vars) = \f$ alphal \f$
  !> .
  !> \param[in]    qp      physical variables  
  !> \param[out]   qc      conservative variables 
  !> \date 2019/11/18
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE qp_to_qc(qp,B,qc)
    
    USE COMPLEXIFY 
    IMPLICIT none
    
    REAL*8, INTENT(IN) :: qp(n_vars)
    REAL*8, INTENT(IN) :: B
    REAL*8, INTENT(OUT) :: qc(n_vars)
    
    REAL*8 :: r_sp_heat_mix
    REAL*8 :: sum_sl
    
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

    IF ( gas_flag ) THEN
       
       ! carrier phase is gas
       r_rho_c = pres / ( sp_gas_const_a * r_T )

       sp_heat_c = sp_heat_a

    ELSE

       ! carrier phase is liquid
       r_rho_c = rho_l

       sp_heat_c = sp_heat_l

    END IF
       
    IF ( gas_flag .AND. liquid_flag ) THEN
       
       ! mixture of gas, liquid and solid
       r_alphal = qp(n_vars)

       ! check and correction on dispersed phases volume fractions
       IF ( ( SUM(r_alphas) + r_alphal ) .GT. 1.D0 ) THEN
          
          sum_sl = SUM(r_alphas) + r_alphal
          r_alphas(1:n_solid) = r_alphas(1:n_solid) / sum_sl
          r_alphal = r_alphal / sum_sl

       ELSEIF ( ( SUM(r_alphas) + r_alphal ) .LT. 0.D0 ) THEN

          r_alphas(1:n_solid) = 0.D0
          r_alphal = 0.D0

       END IF

       ! carrier phase volume fraction
       r_alphac = 1.D0 - SUM(r_alphas) - r_alphal

       ! volume averaged mixture density: carrier (gas) + solids + liquid
       r_rho_m = r_alphac * r_rho_c + SUM( r_alphas * rho_s ) + r_alphal * rho_l

       ! liquid mass fraction
       r_xl = r_alphal * rho_l / r_rho_m

       ! solid mass fractions
       r_xs(1:n_solid) = r_alphas(1:n_solid) * rho_s(1:n_solid) / r_rho_m

       ! carrier (gas) mass fraction
       r_xc = r_alphac * r_rho_c / r_rho_m

       ! mass averaged mixture specific heat
       r_sp_heat_mix =  SUM( r_xs * sp_heat_s ) + r_xl*sp_heat_l + r_xc*sp_heat_c
       
    ELSE

       ! mixture of carrier phase ( gas or liquid ) and solid

       ! check and corrections on dispersed phases
       IF ( SUM(r_alphas) .GT. 1.D0 ) THEN
          
          r_alphas(1:n_solid) = r_alphas(1:n_solid) / SUM(r_alphas)

       ELSEIF ( SUM(r_alphas).LT. 0.D0 ) THEN

          r_alphas(1:n_solid) = 0.D0

       END IF

       ! carrier (gas or liquid) volume fraction
       r_alphac = 1.D0 - SUM(r_alphas) 

       ! volume averaged mixture density: carrier (gas or liquid) + solids
       r_rho_m = r_alphac * r_rho_c + SUM( r_alphas * rho_s ) 

       ! solid mass fractions
       r_xs(1:n_solid) = r_alphas(1:n_solid) * rho_s(1:n_solid) / r_rho_m

       ! carrier (gas or liquid) mass fraction
       r_xc = r_alphac * r_rho_c / r_rho_m

       ! mass averaged mixture specific heat
       r_sp_heat_mix =  SUM( r_xs * sp_heat_s ) + r_xc * sp_heat_c
       
    END IF
    
    qc(1) = r_h * r_rho_m
    qc(2) = r_h * r_rho_m * r_u
    qc(3) = r_h * r_rho_m * r_v

    IF ( energy_flag ) THEN

       ! total energy (internal and kinetic)
       qc(4) = r_h * r_rho_m * ( r_sp_heat_mix * r_T                            &
            + 0.5D0 * ( r_u**2 + r_v**2 ) ) 

    ELSE

       ! internal energy
       qc(4) = r_h * r_rho_m * r_sp_heat_mix * r_T 

    END IF

    qc(5:4+n_solid) = r_h * r_alphas(1:n_solid) * rho_s(1:n_solid)

    IF ( gas_flag .AND. liquid_flag ) qc(n_vars) = r_h * r_alphal * rho_l
    
  END SUBROUTINE qp_to_qc

  
  !******************************************************************************
  !> \brief Additional Physical variables
  !
  !> This subroutine evaluates from the physical local variables qpj, the two
  !> additional local variables qp2j = (u,v). 
  !> \param[in]    qpj    real-valued physical variables 
  !> \param[out]   qp2j   real-valued physical variables 
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 10/10/2019
  !******************************************************************************
 
  SUBROUTINE qp_to_qp2(qpj,Bj,qp2j)

    IMPLICIT none

    REAL*8, INTENT(IN) :: qpj(n_vars)
    REAL*8, INTENT(IN) :: Bj
    REAL*8, INTENT(OUT) :: qp2j(3)

    qp2j(1) = qpj(1) + Bj
    
    IF ( qpj(1) .LE. 0.D0 ) THEN

       qp2j(2) = 0.D0
       qp2j(3) = 0.D0

    ELSE

       qp2j(2) = qpj(2) / qpj(1)
       qp2j(3) = qpj(3) / qpj(1)

    END IF

    RETURN

  END SUBROUTINE qp_to_qp2

  !******************************************************************************
  !> \brief Local Characteristic speeds x direction
  !
  !> This subroutine computes from the physical variable evaluates the largest 
  !> positive and negative characteristic speed in the x-direction. 
  !> \param[in]     qpj           array of local physical variables
  !> \param[out]    vel_min       minimum x-velocity
  !> \param[out]    vel_max       maximum x-velocity
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 05/12/2017
  !******************************************************************************

  SUBROUTINE eval_local_speeds_x(qpj,vel_min,vel_max)
  
    IMPLICIT none
    
    REAL*8, INTENT(IN) :: qpj(n_vars)

    REAL*8, INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)

    CALL mixt_var(qpj)
   
    IF ( r_red_grav * r_h .LT. 0.D0 ) THEN

       vel_min(1:n_eqns) = r_u
       vel_max(1:n_eqns) = r_u
  
    ELSE

       vel_min(1:n_eqns) = r_u - DSQRT( r_red_grav * r_h )
       vel_max(1:n_eqns) = r_u + DSQRT( r_red_grav * r_h )
    
    END IF

    RETURN

  END SUBROUTINE eval_local_speeds_x

  !******************************************************************************
  !> \brief Local Characteristic speeds y direction
  !
  !> This subroutine computes from the physical variable evaluates the largest 
  !> positive and negative characteristic speed in the y-direction. 
  !> \param[in]     qpj           array of local physical variables
  !> \param[out]    vel_min       minimum y-velocity
  !> \param[out]    vel_max       maximum y-velocity
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 05/12/2017
  !******************************************************************************

  SUBROUTINE eval_local_speeds_y(qpj,vel_min,vel_max)
    
    IMPLICIT none
    
    REAL*8, INTENT(IN)  :: qpj(n_vars)
    REAL*8, INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)
    
    CALL mixt_var(qpj)
    
    IF ( r_red_grav * r_h .LT. 0.D0 ) THEN

       vel_min(1:n_eqns) = r_v
       vel_max(1:n_eqns) = r_v
  
    ELSE
      
       vel_min(1:n_eqns) = r_v - DSQRT( r_red_grav * r_h )
       vel_max(1:n_eqns) = r_v + DSQRT( r_red_grav * r_h )
       
    END IF

    RETURN

  END SUBROUTINE eval_local_speeds_y
  

  !******************************************************************************
  !> \brief Hyperbolic Fluxes
  !
  !> This subroutine evaluates the numerical fluxes given the conservative
  !> variables qcj and physical variables qpj.
  !> \date 01/06/2012
  !> \param[in]     qcj      real local conservative variables 
  !> \param[in]     qpj      real local physical variables 
  !> \param[in]     Bj       topography
  !> \param[in]     dir      direction of the flux (1=x,2=y)
  !> \param[out]    flux     real  fluxes    
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE eval_fluxes(qcj,qpj,Bj,dir,flux)
    
    IMPLICIT none

    REAL*8, INTENT(IN) :: qcj(n_vars)
    REAL*8, INTENT(IN) :: qpj(n_vars)
    REAL*8, INTENT(IN) :: Bj
    INTEGER, INTENT(IN) :: dir

    REAL*8, INTENT(OUT) :: flux(n_eqns)

    CALL mixt_var(qpj)

    pos_thick:IF ( qcj(1) .NE. 0.D0 ) THEN

       IF ( dir .EQ. 1 ) THEN

          ! Mass flux in x-direction: u * ( rhom * h )
          flux(1) = r_u * qcj(1)
          
          ! x-momentum flux in x-direction + hydrostatic pressure term
          flux(2) = r_u * qcj(2) + 0.5D0 * r_rho_m * r_red_grav * r_h**2

          ! y-momentum flux in x-direction: u * ( rho * h * v )
          flux(3) = r_u * qcj(3)

          IF ( energy_flag ) THEN

             ! ENERGY flux in x-direction
             flux(4) = r_u * ( qcj(4) + 0.5D0 * r_rho_m * r_red_grav * r_h**2 )

          ELSE

             ! Temperature flux in x-direction: u * ( h * T )
             flux(4) = r_u * qcj(4)

          END IF

          ! Mass flux of solid in x-direction: u * ( h * alphas * rhos )
          flux(5:4+n_solid) = r_u * qcj(5:4+n_solid)

          ! Solid flux can't be larger than total flux
          IF ( ( flux(1) .GT. 0.D0 ) .AND. ( SUM(flux(5:4+n_solid)) / flux(1)   &
               .GT. 1.D0 ) ) THEN

             flux(5:4+n_solid) = &
               flux(5:4+n_solid) / SUM(flux(5:4+n_solid)) * flux(1)

          END IF

          ! Mass flux of liquid in x-direction: u * ( h * alphal * rhol )
          IF ( gas_flag .AND. liquid_flag ) flux(n_vars) = r_u * qcj(n_vars)
          
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

          ! Mass flux of solid in y-direction: v * ( h * alphas * rhos )
          flux(5:4+n_solid) = r_v * qcj(5:4+n_solid)
          
          ! Solid flux can't be larger than total flux
          IF ( ( flux(1) .GT. 0.D0 ) .AND. ( SUM(flux(5:4+n_solid)) / flux(1)   &
               .GT. 1.D0 ) ) THEN

             flux(5:4+n_solid) = &
               flux(5:4+n_solid) / SUM(flux(5:4+n_solid)) * flux(1)

          END IF

          ! Mass flux of liquid in x-direction: u * ( h * alphal * rhol )
          IF ( gas_flag .AND. liquid_flag ) flux(n_vars) = r_v * qcj(n_vars)
  
       END IF

    ELSE
       
       flux(1:n_eqns) = 0.D0
       
    ENDIF pos_thick
 
    RETURN

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
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
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

    COMPLEX*16 :: forces_term(n_eqns)

    INTEGER :: i

    COMPLEX*16 :: mod_vel
    
    COMPLEX*16 :: gamma

    REAL*8 :: h_threshold

    !--- Lahars rheology model variables

    !> Temperature in C
    COMPLEX*16 :: Tc
    
    COMPLEX*16 :: expA , expB
    
    !> Fluid dynamic viscosity (units: kg m-1 s-1 )
    COMPLEX*16 :: fluid_visc

    !> Total friction (dimensionless)
    COMPLEX*16 :: s_f

    !> Viscous slope component of total Friction (dimensionless)
    COMPLEX*16 :: s_v

    !> Turbulent dispersive slope component of total friction (dimensionless)
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

          ! alpha1 here has units: kg m-1 s-1
          ! in Table 2 from O'Brien 1988, the values reported have different
          ! units ( poises). 1poises = 0.1 kg m-1 s-1

          h_threshold = 1.D-20

          ! convert from Kelvin to Celsius
          Tc = T - 273.15D0

          ! the dependance of viscosity on temperature is modeled with the
          ! equation presented at:
          ! https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781118131473.app3
          !
          ! In addition, we use a reference value provided in input at a 
          ! reference temperature. This value is used to scale the equation
          IF ( Tc .LT. 20.D0 ) THEN
          
             expA = 1301.D0 / ( 998.333D0 + 8.1855D0 * ( Tc - 20.D0 )           &
                  + 0.00585D0 * ( Tc - 20.D0 )**2 ) - 1.30223D0
             
             alpha1 = alpha1_coeff * 1.D-3 * 10.D0**expA

          ELSE

             expB = ( 1.3272D0 * ( 20.D0 - Tc ) - 0.001053D0 *                  &
                  ( Tc - 20.D0 )**2 ) / ( Tc + 105.0D0 )
             
             alpha1 = alpha1_coeff * 1.002D-3 * 10.D0**expB 

          END IF
          
          ! Fluid viscosity 
          fluid_visc = alpha1 * CDEXP( beta1 * SUM(alphas) )

          IF ( DBLE(h) .GT. h_threshold ) THEN
             
             ! Viscous slope component (dimensionless)
             s_v = Kappa * fluid_visc * mod_vel / ( 8.D0 * rho_m * grav *h**2 )
             
             ! Turbulent dispersive component (dimensionless)
             s_td = n_td**2 * mod_vel**2 / ( h**(4.D0/3.D0) )
          
          ELSE
             
             ! Viscous slope component (dimensionless)
             s_v = Kappa * fluid_visc * mod_vel / ( 8.D0 * h_threshold**2 )
             
             ! Turbulent dispersive components (dimensionless)
             s_td = n_td**2 * (mod_vel**2) / ( h_threshold**(4.D0/3.D0) )
             
          END IF
          
          ! Total implicit friction slope (dimensionless)
          s_f = s_v + s_td
         
          IF ( mod_vel .GT. 0.D0 ) THEN

             ! same units of dqc(2)/dt: kg m-1 s-1
             forces_term(2) = forces_term(2) - grav * rho_m * h *               &
                  ( u / mod_vel ) * s_f

             ! same units of dqc(3)/dt: kg m-1 s-1
             forces_term(3) = forces_term(3) - grav * rho_m * h *               &
                  ( v / mod_vel ) * s_f
  
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

    nh_term = forces_term

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
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
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

          ! Yield strength
          tau_y = alpha2 * CDEXP( beta2 * SUM(alphas) )

          IF ( h .GT. h_threshold ) THEN
             
             ! Yield slope component (dimensionless)
             s_y = tau_y / ( grav * rho_m * h )
                       
          ELSE
             
             ! Yield slope component
             s_y = tau_y /  ( grav * rho_m * h_threshold )
                          
          END IF
  
          IF ( mod_vel .GT. 0.D0 ) THEN

             ! units of dqc(2)/dt
             forces_term(2) = forces_term(2) - grav * rho_m * h *               &
                  ( u / mod_vel ) * s_y

             ! units of dqc(3)/dt
             forces_term(3) = forces_term(3) - grav * rho_m * h *               &
                  ( v / mod_vel ) * s_y
  
          END IF

       ELSEIF ( rheology_model .EQ. 5 ) THEN

          
       ENDIF
              
    ENDIF
   
    IF ( present(c_qj) .AND. present(c_nh_semi_impl_term) ) THEN

       c_nh_semi_impl_term = forces_term

    ELSEIF ( present(r_qj) .AND. present(r_nh_semi_impl_term) ) THEN

       r_nh_semi_impl_term = DBLE( forces_term )

    END IF

    RETURN

  END SUBROUTINE eval_nh_semi_impl_terms

  
  !******************************************************************************
  !> \brief Explicit Forces term
  !
  !> This subroutine evaluates the non-hyperbolic terms to be treated explicitely
  !> in the DIRK numerical scheme (e.g. gravity,source of mass). The sign of the
  !> terms is taken with the terms on the left-hand side of the equations.
  !> \date 01/06/2012
  !> \param[in]     B_primej_x         local x-slope
  !> \param[in]     B_primej_y         local y_slope
  !> \param[in]     qj                 local source
  !> \param[in]     qpj                physical variables 
  !> \param[in]     qcj                conservative variables 
  !> \param[out]    expl_term          explicit term
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE eval_expl_terms( Bprimej_x , Bprimej_y , source_xy , qpj , qcj ,  &
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

  END SUBROUTINE eval_expl_terms


  !******************************************************************************
  !> \brief Erosion/Deposition term
  !
  !> This subroutine evaluates the deposition term.
  !> \date 03/010/2018
  !> \param[in]     qpj                local physical variables 
  !> \param[in]     Bj                 local topography
  !> \param[in]     dt                 time step
  !> \param[out]    erosion_term       erosion term for each solid phase
  !> \param[out]    dep_term           deposition term for each solid phase
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
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
        
    DO i_solid=1,n_solid

       IF ( ( r_alphas(i_solid) .GT. 0.D0 ) .AND. ( settling_flag ) ) THEN

          settling_vel = settling_velocity( diam_s(i_solid) , rho_s(i_solid) ,  &
               r_rho_c , i_solid )

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
  
       IF ( r_h .GT. 1.D-2) THEN
    
          ! empirical formulation (see Fagents & Baloga 2006, Eq. 5)
          ! here we use the solid volume fraction instead of relative density
          ! This term has units: m s-1
          erosion_term(i_solid) = erosion_coeff(i_solid) * mod_vel * r_h        &
               * ( 1.D0 - SUM(r_alphas) )

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
  !> \date 2019/11/08
  !> \param[in]     qj                  conservative variables 
  !> \param[in]     deposition_avg_term averaged deposition terms 
  !> \param[in]     erosion_avg_term    averaged deposition terms 
  !> \param[out]    topo_term           explicit term
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE eval_topo_term( qpj , deposition_avg_term , erosion_avg_term ,      &
       eqns_term, deposit_term )
    
    IMPLICIT NONE
    
    REAL*8, INTENT(IN) :: qpj(n_eqns)                  !< physical variables 
    REAL*8, INTENT(IN) :: deposition_avg_term(n_solid) !< deposition term
    REAL*8, INTENT(IN) :: erosion_avg_term(n_solid)    !< erosion term

    REAL*8, INTENT(OUT):: eqns_term(n_eqns)
    REAL*8, INTENT(OUT):: deposit_term(n_solid)
   
    REAL*8 :: entr_coeff
    REAL*8 :: air_entr
    REAL*8 :: mag_vel 

    CALL mixt_var(qpj)

    mag_vel = DSQRT( r_u**2.D0 + r_v**2.D0 ) 
    
    IF ( entrainment_flag .AND. ( mag_vel**2 .GT. 0.D0 ) .AND.                  &
         ( r_h .GT. 0.D0 ) ) THEN

       entr_coeff = 0.075D0 / DSQRT( 1.D0 + 718.D0 * MAX(0.D0,r_Ri)**2.4 )
       
       air_entr = entr_coeff * mag_vel
       
    ELSE

       air_entr = 0.D0

    END IF

    eqns_term(1:n_eqns) = 0.D0

    ! free surface (topography+flow) equation
    eqns_term(1) = SUM( rho_s * ( erosion_avg_term - deposition_avg_term ) ) +  &
         rho_a_amb * air_entr

    ! x-momenutm equation
    eqns_term(2) = - r_u * SUM( rho_s * deposition_avg_term )

    ! y-momentum equation
    eqns_term(3) = - r_v * SUM( rho_s * deposition_avg_term )

    ! Temperature/Energy equation
    IF ( energy_flag ) THEN

       eqns_term(4) = - r_T * SUM( rho_s * sp_heat_s * deposition_avg_term )    &
            - 0.5D0 * mag_vel**2 * SUM( rho_s * deposition_avg_term )           &
            + T_s_substrate * SUM( rho_s * sp_heat_s * erosion_avg_term )       &
            + T_ambient * sp_heat_a * rho_a_amb * air_entr

    ELSE

       eqns_term(4) = - r_T * SUM( rho_s * sp_heat_s * deposition_avg_term )    &
            + T_s_substrate * SUM( rho_s * sp_heat_s * erosion_avg_term )       &
            + T_ambient * sp_heat_a * rho_a_amb * air_entr

    END IF
   
    ! solid phase thickness equation
    eqns_term(5:4+n_solid) = rho_s(1:n_solid) * ( erosion_avg_term(1:n_solid)   &
         - deposition_avg_term(1:n_solid) )

    ! solid deposit rate terms
    deposit_term(1:n_solid) = deposition_avg_term(1:n_solid)                    &
         - erosion_avg_term(1:n_solid)

    RETURN

  END SUBROUTINE eval_topo_term

  !******************************************************************************
  !> \brief Internal boundary source fluxes
  !
  !> This subroutine evaluates the source terms at the interfaces when an
  !> internal radial source is present, as for a base surge. The terms are 
  !> applied as boundary conditions, and thus they have the units of the 
  !> physical variable qp
  !> \date 2019/12/01
  !> \param[in]     vect_x       unit vector velocity x-component 
  !> \param[in]     vect_y       unit vector velocity y-component 
  !> \param[out]    source_bdry  source terms  
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************
  
  SUBROUTINE eval_source_bdry( vect_x , vect_y , source_bdry )

    USE parameters_2d, ONLY : h_source , vel_source , T_source , alphas_source ,&
         alphal_source

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: vect_x
    REAL*8, INTENT(IN) :: vect_y
    REAL*8, INTENT(OUT) :: source_bdry(n_vars)

    source_bdry(1) = h_source
    source_bdry(2) = h_source * vel_source * vect_x
    source_bdry(3) = h_source * vel_source * vect_y
    source_bdry(4) = T_source
    source_bdry(5:4+n_solid) = alphas_source(1:n_solid)

    IF ( gas_flag .AND. liquid_flag ) THEN

       source_bdry(n_vars) = alphal_source

    END IF
    
    
    RETURN

  END SUBROUTINE eval_source_bdry

  !------------------------------------------------------------------------------
  !> Settling velocity function
  !
  !> This subroutine compute the settling velocity of the particles, as a
  !> function of diameter, density.
  !> \date 2019/11/11
  !> \param[in]    diam        particle diameter      
  !> \param[in]    rhos        particle density
  !> \param[in]    rhoa        atmospheric density
  !> \patam[in]    i_solid     particle class index
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !------------------------------------------------------------------------------
  
  REAL*8 FUNCTION settling_velocity(diam,rhos,rhoc,i_solid)

    IMPLICIT NONE

    REAL*8, INTENT(IN) :: diam          !< particle diameter [m]
    REAL*8, INTENT(IN) :: rhos          !< particle density [kg/m3]
    REAL*8, INTENT(IN) :: rhoc          !< carrier phase density [kg/m3]
    INTEGER, INTENT(IN) :: i_solid      !< particle class index

    REAL*8 :: Rey           !< Reynolds number
    REAL*8 :: C_D           !< Drag coefficient

    ! loop variables
    INTEGER :: i            !< loop counter for iterative procedure    
    REAL*8 :: const_part    !< term not changing in iterative procedure
    REAL*8 :: C_D_old       !< previous iteration drag coefficient
    REAL*8 :: set_vel_old   !< previous iteration settling velocity
    
    INTEGER :: dig          !< order of magnitude of settling velocity

    C_D = 1.D0

    const_part =  DSQRT( 0.75D0 * ( rhos / rhoc - 1.D0 ) * diam * grav )
  
    settling_velocity = const_part / DSQRT( C_D )

    Rey = diam * settling_velocity / kin_visc_c

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
             settling_velocity = 10.D0**(dig-3)                                 &
                  * FLOOR( 10.0**(-dig+3)*set_vel_old ) 

             EXIT C_D_loop

          END IF

          Rey = diam * settling_velocity / kin_visc_c

       END DO C_D_loop
    
    END IF

    C_D_s(i_solid) = C_D

    RETURN

  END FUNCTION settling_velocity


END MODULE constitutive_2d

    
