!********************************************************************************
!> \brief Constitutive equations
!********************************************************************************
MODULE constitutive_2d

  USE parameters_2d, ONLY : wp, sp ,tolh
  USE parameters_2d, ONLY : n_eqns , n_vars , n_solid , n_add_gas , n_quad ,    &
       n_stoch_vars , n_pore_vars
  USE parameters_2d, ONLY : rheology_flag , rheology_model , energy_flag ,      &
       liquid_flag , gas_flag , alpha_flag , slope_correction_flag ,            &
       curvature_term_flag, stochastic_flag, mean_field_flag,                   &
       stoch_transport_flag, pore_pressure_flag , sutherland_flag

  USE parameters_2d, ONLY : idx_h, idx_hu, idx_hv, idx_T, idx_alfas_first,      &
       idx_alfas_last, idx_addGas_first, idx_addGas_last, idx_stoch, idx_pore,  &
       idx_u, idx_v

  USE parameters_2d, ONLY : idx_totMassEqn, idx_uEqn, idx_vEqn, idx_engyEqn,    &
       idx_solidEqn_first, idx_solidEqn_last, idx_addGasEqn_first,              &
       idx_addGasEqn_last, idx_stochEqn, idx_poreEqn
  
  IMPLICIT none

  !> flag used for size of implicit non linear-system
  LOGICAL, ALLOCATABLE :: implicit_flag(:)

  !> map from implicit variables to original variables 
  INTEGER, ALLOCATABLE :: implicit_map(:)

  !> flag to activate air entrainment
  LOGICAL :: entrainment_flag

  !> gravitational acceleration 
  REAL(wp) :: grav
  REAL(wp) :: inv_grav

  !> drag coefficients (Voellmy-Salm model)
  REAL(wp) :: mu
  REAL(wp) :: xi
  REAL(wp) :: xi_temp

  !> friction coefficients (function Coulomb: mu(Fr) and mu:(U))
  REAL(wp) :: mu_0
  REAL(wp) :: mu_inf
  REAL(wp) :: Fr_0
  REAL(wp) :: U_w

  !> drag coefficients (B&W model)
  REAL(wp) :: friction_factor

  !> drag coefficients (plastic model)
  REAL(wp) :: tau

  !> evironment temperature [K]
  REAL(wp) :: T_env

  !> reference temperature [K]
  REAL(wp) :: T_ref

  !> reference kinematic viscosity [m2/s]
  REAL(wp) :: nu_ref

  !> viscosity parameter [K-1] (b in Table 1 Costa & Macedonio, 2005)
  REAL(wp) :: visc_par

  !> yield strength for lava rheology [kg m-1 s-2] (Eq.4 Kelfoun & Varga, 2015)
  REAL(wp) :: tau0

  !> velocity boundary layer fraction of total thickness
  REAL(wp) :: emme

  !> specific heat [J kg-1 K-1]
  REAL(wp) :: c_p

  !> coefficient for the convective term in the temperature equation for lava
  REAL(wp) :: convective_term_coeff

  !> atmospheric heat trasnfer coefficient [W m-2 K-1] (lambda in C&M, 2005)
  REAL(wp) :: atm_heat_transf_coeff

  !> fractional area of the exposed inner core (f in C&M, 2005)
  REAL(wp) :: exp_area_fract

  !> coefficient for the radiative term in the temperature equation for lava
  REAL(wp) :: radiative_term_coeff

  !> Stephan-Boltzmann constant [W m-2 K-4]
  REAL(wp), PARAMETER :: SBconst = 5.67E-8_wp

  !> emissivity (eps in Costa & Macedonio, 2005)
  REAL(wp) :: emissivity

  !> thermal boundary layer fraction of total thickness
  REAL(wp) :: enne

  !> temperature of lava-ground interface
  REAL(wp) :: T_ground

  !> thermal conductivity [W m-1 K-1] (k in Costa & Macedonio, 2005)
  REAL(wp) :: thermal_conductivity

  !--- START Lahars rheology model parameters

  !> 1st param for yield strenght empirical relationship (O'Brian et al, 1993)
  REAL(wp) :: alpha2    ! (units: kg m-1 s-2)

  !> 2nd param for yield strenght empirical relationship (O'Brian et al, 1993)
  REAL(wp) :: beta2     ! (units: nondimensional) 

  !> ratio between reference value from input and computed values from eq.
  REAL(wp) :: alpha1_coeff ! (units: nondimensional )

  !> 2nd param for fluid viscosity empirical relationship (O'Brian et al, 1993)
  REAL(wp) :: beta1     ! (units: nondimensional, input parameter)

  !> Empirical resistance parameter (dimensionless, input parameter)
  REAL(wp) :: Kappa

  !> Mannings roughness coefficient ( units: T L^(-1/3) )
  REAL(wp) :: n_td
  REAL(wp) :: n_td2

  !--- END Lahars rheology model parameters

  !> Specific heat of carrier phase (gas or liquid)
  REAL(wp) :: sp_heat_c  ! ( initialized from input)   

  !> Specific gas constant of gas mixture (units: J kg-1 K-1)
  REAL(wp) :: sp_gas_const_c

  !> Density of carrier phase in substrate ( units: kg m-3 )
  REAL(wp) :: rho_c_sub

  !> Ambient density of air ( units: kg m-3 )
  REAL(wp) :: rho_a_amb

  !> Specific heat of air (units: J K-1 kg-1)
  REAL(wp) :: sp_heat_a

  !> Specific gas constant of air (units: J kg-1 K-1)
  REAL(wp) :: sp_gas_const_a

  !> Kinematic viscosity of air (units: m2 s-1)
  REAL(wp) :: kin_visc_a

  !> Specific heat of additional gas (units: J K-1 kg-1)
  REAL(wp), ALLOCATABLE :: sp_heat_g(:)

  !> Specific gas constant of additional gas (units: J kg-1 K-1)
  REAL(wp), ALLOCATABLE :: sp_gas_const_g(:)

  !> Kinematic viscosity of liquid (units: m2 s-1)
  REAL(wp) :: kin_visc_l

  !> Kinematic viscosity of carrier phase (units: m2 s-1)
  REAL(wp) :: kin_visc_c

  !> Reference temperature for Sutherland’s law (units: K)
  REAL(wp) :: Tref_Suth

  !> Reference viscosity for Sutherland's law (units: kg m-1 s-1)
  REAL(wp) :: muRef_Suth

  !> Sutherland constant for Sutherland’s law (units: K) 
  REAL(wp) :: S_mu
  
  !> Temperature of ambient air (units: K)
  REAL(wp) :: T_ambient

  !> Density of sediments ( units: kg m-3 )
  REAL(wp), ALLOCATABLE :: rho_s(:)

  !> Reciprocal of density of sediments ( units: kg m-3 )
  REAL(wp), ALLOCATABLE :: inv_rho_s(:)

  !> Reciprocal of density of sediments ( units: kg m-3 )
  COMPLEX(wp), ALLOCATABLE :: c_inv_rho_s(:)

  !> Diameter of sediments ( units: m )
  REAL(wp), ALLOCATABLE :: diam_s(:)

  !> Specific heat of solids (units: J K-1 kg-1)
  REAL(wp), ALLOCATABLE :: sp_heat_s(:)

  !> Flag to determine if sedimentation is active
  LOGICAL :: settling_flag

  !> Minimum volume fraction of solids in the flow
  REAL(wp) :: alphastot_min

  !> erosion model coefficient  (units: m-1 )
  REAL(wp), ALLOCATABLE :: erosion_coeff

  !> water loss_rate (unit: m s-1)
  REAL(wp), ALLOCATABLE :: loss_rate

  !> erodible substrate solid relative volume fractions
  REAL(wp), ALLOCATABLE :: erodible_fract(:)

  !> erodible substrate porosity (we assume filled by continous phase)
  REAL(wp) :: erodible_porosity

  !> coefficient to compute (eroded/deposited) volume of continuous phase
  !> from volume of solid  
  REAL(wp) :: coeff_porosity

  !> temperature of erodible substrate (units: K)
  REAL(wp) :: T_erodible

  !> ambient pressure (units: Pa)
  REAL(wp) :: pres

  !> reciprocal of ambient pressure (units: Pa)
  REAL(wp) :: inv_pres

  !> liquid density (units: kg m-3)
  REAL(wp) :: rho_l

  !> reciprocal of liquid density (units: kg m-3)
  REAL(wp) :: inv_rho_l

  !> Sepcific heat of liquid (units: J K-1 kg-1)
  REAL(wp) :: sp_heat_l

  !> Fraction of heat lost by particles producing steam
  REAL(wp) :: gamma_steam

  !> Von Karman constant
  REAL(wp) :: vonK

  !> Substrate Roughness (units: m) 
  REAL(wp) :: k_s

  REAL(wp) :: H_crit_rel

  !> Schmidt number: ratio of momentum and mass diffusivity
  REAL(wp) :: Sc

  REAL(wp) :: z_dyn

  REAL(wp) :: hydraulic_permeability

  REAL(wp) :: maximum_solid_packing
  

  INTERFACE u_log_profile    ! Define generic function
     MODULE PROCEDURE u_log_profile_scalar
     MODULE PROCEDURE u_log_profile_array
  END INTERFACE u_log_profile


  INTERFACE alphas_exp_profile    ! Define generic function
     MODULE PROCEDURE alphas_exp_profile_scalar
     MODULE PROCEDURE alphas_exp_profile_array
  END INTERFACE alphas_exp_profile
  
  
CONTAINS
  
  FUNCTION u_log_profile_scalar(b,z)

    REAL(wp) :: u_log_profile_scalar
    REAL(wp), INTENT(IN) :: b
    REAL(wp), INTENT(IN) :: z
    
    u_log_profile_scalar = LOG( b*z + 1.0_wp )
    
  END FUNCTION u_log_profile_scalar

  FUNCTION u_log_profile_array(b,z)

    REAL(wp) :: u_log_profile_array(n_quad)
    REAL(wp), INTENT(IN) :: b
    REAL(wp), INTENT(IN) :: z(n_quad)
    
    u_log_profile_array = LOG( b*z + 1.0_wp )
    
  END FUNCTION u_log_profile_array

  FUNCTION alphas_exp_profile_scalar(a,z)

    REAL(wp) :: alphas_exp_profile_scalar
    REAL(wp), INTENT(IN) :: a
    REAL(wp), INTENT(IN) :: z
    
    alphas_exp_profile_scalar = EXP( a*z )
    
  END FUNCTION alphas_exp_profile_scalar

  FUNCTION alphas_exp_profile_array(a,z)

    REAL(wp) :: alphas_exp_profile_array(n_quad)
    REAL(wp), INTENT(IN) :: a
    REAL(wp), INTENT(IN) :: z(n_quad)
    
    alphas_exp_profile_array = EXP( a*z )
    
  END FUNCTION alphas_exp_profile_array

  FUNCTION dynamic_pressure(rho_c, alphas, normalizing_coeff_alpha , a,  &
       u, normalizing_coeff_u , b, z)  

    !> dynamic pressure
    REAL(wp) :: dynamic_pressure


    !> density of the carrier phase
    REAL(wp), INTENT(IN) :: rho_c
    !> depth-averaged volumetric fractions of solid phases
    REAL(wp), INTENT(IN) :: alphas(n_solid)
    !> normalizing coefficients of the solid profile functions
    REAL(wp), INTENT(IN) :: normalizing_coeff_alpha(n_solid)
    !> shape parameters of the solid profile function
    REAL(wp), INTENT(IN) :: a(n_solid)
    ! depth-averaged velocity
    REAL(wp), INTENT(IN) :: u
    !> normalizing coefficients of the velocity profile function
    REAL(wp), INTENT(IN) :: normalizing_coeff_u
    !> shape parameter of the velocity profile function
    REAL(wp), INTENT(IN) :: b
    !> elevation at which the dynamic pressure is computed
    REAL(wp), INTENT(IN) :: z

    !> volumetric fractions of solid at z
    REAL(wp) :: alphas_z(n_solid)
    !> mixture density at z
    REAL(wp) :: rhom_z
    !> mixture velocity at z
    REAL(wp) :: u_z

    !> volumetric fraction of carrier phase at z
    REAL(wp) :: alphac_z

    !> loop counter
    INTEGER :: i_solid
    
    DO i_solid = 1,n_solid

       alphas_z(i_solid) = alphas(i_solid) *                                    &
            normalizing_coeff_alpha(i_solid) * alphas_exp_profile(a(i_solid),z)
       
    END DO

    alphac_z = 1.0_wp - SUM(alphas_z)
    
    rhom_z = alphac_z * rho_c + SUM( rho_s * alphas_z )

    u_z = ( u * normalizing_coeff_u * u_log_profile(b,z) )

    dynamic_pressure = 0.5_wp * rhom_z * u_z**2

  END FUNCTION dynamic_pressure
  
  
  !******************************************************************************
  !> \brief Initialization of relaxation flags
  !
  !> This subroutine set the number and the flags of the non-hyperbolic
  !> terms.
  !> \date 07/09/2012
  !******************************************************************************

  SUBROUTINE init_problem_param

    USE parameters_2d, ONLY : n_nh , pore_pressure_flag
    IMPLICIT NONE

    integer :: i,j

    ALLOCATE( implicit_flag(n_eqns) )

    implicit_flag(1:n_eqns) = .FALSE.
    implicit_flag(2) = .TRUE.
    implicit_flag(3) = .TRUE.

    ! Temperature
    IF ( rheology_model .EQ. 3 ) THEN

       implicit_flag(4) = .TRUE.

    END IF

    ! Solid volume fraction
    implicit_flag(idx_solidEqn_first:idx_solidEqn_last) = .FALSE.

    IF ( pore_pressure_flag ) THEN

       implicit_flag(idx_pore) = .TRUE.
       
    END IF
    
    n_nh = COUNT( implicit_flag )

    ALLOCATE( implicit_map(n_nh) )

    j=0
    DO i=1,n_eqns

       IF ( implicit_flag(i) ) THEN

          j=j+1
          implicit_map(j) = i

       END IF

    END DO

    WRITE(*,*) 'Implicit equations =',n_nh

    RETURN

  END SUBROUTINE init_problem_param

  !******************************************************************************
  !> \brief Physical variables
  !
  !> This subroutine evaluates from the conservative local variables qj
  !> the local physical variables  (\f$h,u,v,\alpha_s,\rho_m,T,\alpha_l \f$).
  !> \param[in]    r_qj        real conservative variables 
  !> \param[out]   r_h         real-value flow thickness 
  !> \param[out]   r_u         real-value flow x-velocity 
  !> \param[out]   r_v         real-value flow y-velocity
  !> \param[out]   r_alphas    real-value solid volume fractions
  !> \param[out]   r_rho_m     real-value flow density
  !> \param[out]   r_T         real-value flow temperature 
  !> \param[out]   r_alphal    real-value liquid volume fraction
  !> \param[out]   r_red_grav  real-value reduced gravity
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !> \date 2019/12/13
  !******************************************************************************

  SUBROUTINE r_phys_var(r_qj , r_h , r_u , r_v , r_alphas , r_rho_m , r_T ,     &
       r_alphal , r_alphag , r_red_grav , p_dyn , r_Zs , r_exc_pore_pres)

    USE geometry_2d, ONLY : lambertw , lambertw0 , lambertwm1
    USE geometry_2d, ONLY : z_quad , w_quad

    USE parameters_2d, ONLY : eps_sing , eps_sing4 , vertical_profiles_flag
    IMPLICIT none

    REAL(wp), INTENT(IN) :: r_qj(n_vars)       !< real-value conservative var
    REAL(wp), INTENT(OUT) :: r_h               !< real-value flow thickness
    REAL(wp), INTENT(OUT) :: r_u               !< real-value x-velocity
    REAL(wp), INTENT(OUT) :: r_v               !< real-value y-velocity
    REAL(wp), INTENT(OUT) :: r_alphas(n_solid) !< real-value solid volume fracts
    REAL(wp), INTENT(OUT) :: r_rho_m           !< real-value mixture density
    REAL(wp), INTENT(OUT) :: r_T               !< real-value temperature
    REAL(wp), INTENT(OUT) :: r_alphal          !< real-value liquid volume fract
    REAL(wp), INTENT(OUT) :: r_alphag(n_add_gas) !< real-value gas volume fracts
    REAL(wp), INTENT(OUT) :: r_red_grav        !< real-value reduced gravity
    REAL(wp), INTENT(OUT) :: p_dyn
    REAL(wp), INTENT(OUT) :: r_Zs              !< real-value stochastic variable
    REAL(wp), INTENT(OUT) :: r_exc_pore_pres   !< real-value pore pressure

    REAL(wp) :: r_inv_rhom
    REAL(wp) :: r_xs(n_solid)     !< real-value solid mass fractions
    REAL(wp) :: r_xg(n_add_gas)     !< real-value additional gas mass fractions
    REAL(wp) :: r_xs_tot

    REAL(wp) :: r_Ri            !< real-value Richardson number
    REAL(wp) :: r_xl            !< real-value liquid mass fraction
    REAL(wp) :: r_xc            !< real-value carrier phase mass fraction
    REAL(wp) :: r_alphac        !< real-value carrier phase volume fraction
    REAL(wp) :: r_sp_heat_c     !< real-value specific heat of carrier phase
    REAL(wp) :: r_sp_heat_mix   !< real-value specific heat of mixture
    REAL(wp) :: r_sp_gas_const_c!< real-value gas constant of carrier phase
    REAL(wp) :: r_rho_c         !< real-value carrier phase density [kg/m3]
    REAL(wp) :: r_inv_rho_c
    REAL(wp) :: r_inv_rho_g(n_add_gas)    !< add. gas density reciprocal

    REAL(wp) :: inv_qj1

    REAL(wp) :: rhos_alfas(n_solid)

    REAL(wp) :: rhos_alfas_tot_u 
    REAL(wp) :: rhos_alfas_tot_v

    REAL(wp) :: settling_vel(n_solid)

    REAL(wp) :: inv_kin_visc

    !REAL(wp) :: a_crit_rel
    REAL(wp) :: h_rel

    REAL(wp) :: a,b,c,d

    REAL(wp) :: h0_rel
    REAL(wp) :: h0_rel_1
    REAL(wp) :: h0_rel_2

    REAL(wp) :: normalizing_coeff_u
    
    REAL(wp) :: h0
    REAL(wp) :: u_rel0

    REAL(wp) :: uRho_avg
    REAL(wp) :: uRho_avg_new

    REAL(wp) :: u_avg_guess
    REAL(wp) :: u_avg_new

    REAL(wp) :: rhom_avg

    REAL(wp) :: x0,x1,x2

    INTEGER :: i_aitken
    REAL(wp) :: abs_tol , rel_tol
    REAL(wp) :: denominator
    REAL(wp) :: aitkenX
    REAL(wp) :: lambda

    INTEGER :: i_solid

    REAL(wp) :: z(n_quad)
    REAL(wp) :: w(n_quad)

    REAL(wp) :: u_log_avg

    REAL(wp) :: r_sp_heat_c_by_xc

    REAL(wp) :: dyn_visc_c

    ! compute solid mass fractions
    IF ( r_qj(1) .GT. EPSILON(1.0_wp) ) THEN

       inv_qj1 = 1.0_wp / r_qj(1)

       r_xs(1:n_solid) = r_qj(idx_alfas_first:idx_alfas_last) * inv_qj1

       IF ( SUM( r_qj(idx_alfas_first:idx_alfas_last) ) .EQ. r_qj(1) ) THEN

          r_xs(1:n_solid) = r_xs(1:n_solid) / SUM( r_xs(1:n_solid) )

       END IF

       IF ( n_add_gas .GT. 0 ) r_xg(1:n_add_gas) =                              &
            r_qj(idx_addGas_first:idx_addGas_last) * inv_qj1

       IF ( stoch_transport_flag ) r_Zs = r_qj(idx_stoch) * inv_qj1

       IF ( pore_pressure_flag ) r_exc_pore_pres = r_qj(idx_pore) * inv_qj1
       
    ELSE

       r_h = 0.0_wp
       r_u = 0.0_wp
       r_v = 0.0_wp
       r_alphas = 0.0_wp
       r_rho_m = rho_a_amb
       r_T = T_ambient
       r_alphal = 0.0_wp
       r_alphag = 0.0_wp
       r_red_grav = 0.0_wp
       r_rho_c = rho_a_amb
       p_dyn = 0.0_wp
       r_Zs = 0.0_wp
       r_exc_pore_pres = 0.0_wp

       RETURN

    END IF

    r_xs_tot = SUM(r_xs)

    IF ( gas_flag .AND. liquid_flag ) THEN

       ! compute liquid mass fraction
       r_xl = r_qj(n_vars) * inv_qj1

       ! compute carrier phase (gas) mass fraction
       r_xc =  1.0_wp - r_xs_tot - r_xl

       ! compute specific heat of gas phase (weighted average of specific heat
       ! of gas components, with weights given by mass fractions)
       r_sp_heat_c_by_xc = ( ( r_xc - SUM( r_xg(1:n_add_gas) ) ) * sp_heat_a +  &
            DOT_PRODUCT( r_xg(1:n_add_gas) , sp_heat_g(1:n_add_gas) ) )

       ! specific heat of the mixutre: mass average of sp. heat pf phases
       r_sp_heat_mix = DOT_PRODUCT( r_xs(1:n_solid) , sp_heat_s(1:n_solid) )    &
            + r_xl * sp_heat_l + r_xc * r_sp_heat_c_by_xc

    ELSE

       ! compute carrier phase (gas or liquid) mass fraction
       r_xc = 1.0_wp - r_xs_tot

       IF ( gas_flag ) THEN

          r_sp_heat_c_by_xc = ( ( r_xc - SUM( r_xg(1:n_add_gas) ) ) * sp_heat_a &
               + DOT_PRODUCT( r_xg(1:n_add_gas) , sp_heat_g(1:n_add_gas) ) )

       ELSE

          r_sp_heat_c_by_xc = sp_heat_l * r_xc

       END IF

       ! specific heaf of the mixutre: mass average of sp. heat pf phases
       r_sp_heat_mix = DOT_PRODUCT( r_xs(1:n_solid) , sp_heat_s(1:n_solid) )    &
            + r_sp_heat_c_by_xc

    END IF

    ! compute temperature from energy
    IF ( r_qj(1) .GT. eps_sing ) THEN

       IF ( energy_flag ) THEN

          r_T = ( r_qj(4) - 0.5_wp * ( r_qj(2)**2 + r_qj(3)**2 ) * inv_qj1 ) /  &
               ( r_qj(1) * r_sp_heat_mix ) 

       ELSE

          r_T = r_qj(4) / ( r_qj(1) * r_sp_heat_mix ) 

       END IF

       IF ( r_T .LE. 0.0_wp ) r_T = T_ambient

    ELSE

       r_T = T_ambient

    END IF

    IF ( gas_flag ) THEN

       ! carrier phase is gas
       IF ( r_xc .GT. EPSILON(1.0_wp) ) THEN
          
          r_sp_gas_const_c = ( ( r_xc - SUM( r_xg(1:n_add_gas) ) ) * sp_gas_const_a&
               + DOT_PRODUCT( r_xg(1:n_add_gas) , sp_gas_const_g(1:n_add_gas) ) )  &
               / r_xc

       ELSE

          r_sp_gas_const_c = sp_gas_const_a

       END IF
       
       r_rho_c =  pres / ( r_sp_gas_const_c * r_T )
       r_inv_rho_c = r_sp_gas_const_c * r_T * inv_pres

       r_inv_rho_g(1:n_add_gas) = sp_gas_const_g(1:n_add_gas) * r_T * inv_pres

    ELSE

       r_rho_c = rho_l
       r_inv_rho_c = inv_rho_l
       sp_heat_c = sp_heat_l

    END IF

    ! the liquid contribution (if present) is added below
    r_inv_rhom = DOT_PRODUCT( r_xs(1:n_solid) , inv_rho_s(1:n_solid) )          &
         + r_xc * r_inv_rho_c

    IF ( gas_flag .AND. liquid_flag ) THEN

       r_inv_rhom = r_inv_rhom + r_xl * inv_rho_l

    END IF

    ! mixture density
    r_rho_m = 1.0_wp / r_inv_rhom

    IF ( gas_flag .AND. liquid_flag ) THEN

       r_alphal = r_xl * r_rho_m * inv_rho_l

    END IF

    ! convert from mass fraction to volume fraction
    r_alphas(1:n_solid) = r_xs(1:n_solid) * r_rho_m * inv_rho_s(1:n_solid)

    ! convert from mass fraction to volume fraction
    r_alphag(1:n_add_gas) = r_xg(1:n_add_gas) * r_rho_m                         &
         * r_inv_rho_g(1:n_add_gas)

    ! convert from mass fraction to volume fraction
    r_alphac = r_xc * r_rho_m * r_inv_rho_c

    r_h = r_qj(1) * r_inv_rhom

    ! reduced gravity
    r_red_grav = ( r_rho_m - rho_a_amb ) * r_inv_rhom * grav

    IF ( vertical_profiles_flag ) THEN

       rhos_alfas_tot_u = r_xs_tot * r_qj(2) / r_h  
       rhos_alfas_tot_v = r_xs_tot * r_qj(3) / r_h  

       rhos_alfas(1:n_solid) = r_alphas(1:n_solid) * rho_s(1:n_solid)

       IF ( gas_flag .AND. sutherland_flag ) THEN
          
          dyn_visc_c = muRef_Suth * ( r_T / Tref_Suth )**1.5_wp *               &
               ( Tref_Suth + S_mu ) / ( r_T + S_mu )

          kin_visc_c = dyn_visc_c * r_inv_rho_c

       END IF
                 
       ! Viscosity read from input file [m2 s-1]
       inv_kin_visc = 1.0_wp / kin_visc_c

       DO i_solid=1,n_solid

          settling_vel(i_solid) = settling_velocity( diam_s(i_solid) ,          &
               rho_s(i_solid) , r_rho_c , inv_kin_visc )

       END DO


       ! The profile parameters depend on h/k_s, not on the absolute value of h
       h_rel = r_h / k_s

       IF ( h_rel .GT. H_crit_rel ) THEN

          ! we search for h0_rel such that the average integral between 0 and
          ! h_rel is equal to 1
          ! For h_rel > H_crit_rel this integral is the sum of two pieces:
          ! integral between 0 and h0_rel of the log profile
          ! integral between h0_rel and h_rel of the costant profile

          a = h_rel * vonK / SQRT(friction_factor)
          b = 1.0_wp / 30.0_wp + h_rel
          c = 30.0_wp

          ! solve b*log(c*z+1)-z=a for z
          d =  a / b - 1.0_wp / ( b*c )

          h0_rel_1 = -b*lambertw0( -EXP(d)/(b*c) ) - 1.0_wp / c
          h0_rel_2 = -b*lambertwm1( -EXP(d)/(b*c) ) - 1.0_wp / c
          h0_rel = MIN( h0_rel_1 , h0_rel_2)

       ELSE

          ! when h_rel <= H_crit_rel we have only the log profile and we have to
          ! rescale it in order to have the integral between o and h_rel equal to
          ! 1
          h0_rel = h_rel

       END IF

       h0 = h0_rel * k_s

       b = 30.0_wp / k_s

       ! Quadrature points and weights for the interval [0;h0] 
       z = 0.5_wp * h0 * ( z_quad + 1.0_wp )
       w = 0.5_wp * h0 * w_quad

       u_log_avg = ( SUM( w * u_log_profile(b,z) ) + u_log_profile(b,h0)*(r_h-h0) )  &
            / r_h

       !u_log_avg = ( SUM( w * LOG( b*z + 1.0_wp ) ) + LOG( b*h0 + 1.0_wp )*(r_h-h0) )  &
       !     / r_h
       
       normalizing_coeff_u = 1.0_wp / u_log_avg

       ! velocity at h0
       u_rel0 = normalizing_coeff_u * u_log_profile(b,h0)
       ! u_rel0 = normalizing_coeff_u * LOG( b*h0 + 1.0_wp )

       uRho_avg = SQRT( r_qj(2)**2 + r_qj(3)**2 ) / r_h 

       u_avg_guess = uRho_avg / r_rho_m
       x0 = u_avg_guess

       ! loop to compute the average velocity from average rho*alpha and average 
       ! uRho ( = 1/h*int( u*rhog*alphag + sum[u*rhos(i)*alphas(i)] ) ) 

       rel_tol = 1.e-8
       abs_tol = 1.e-8

       aitken_loop:DO i_aitken=1,10

          x0 = u_avg_guess
         
          CALL avg_profiles_mix( r_h , settling_vel , rhos_alfas(1:n_solid) ,   &
               u_avg_guess , h0 , b , u_rel0 , r_rho_c , rhom_avg ,             &
               uRho_avg_new , p_dyn )

          u_avg_new = u_avg_guess * uRho_avg / ( uRho_avg_new)

          x1 = u_avg_new

          CALL avg_profiles_mix( r_h , settling_vel , rhos_alfas(1:n_solid) ,   &
               u_avg_new , h0 , b , u_rel0 , r_rho_c , rhom_avg ,               &
               uRho_avg_new , p_dyn )

          u_avg_new = u_avg_new * uRho_avg / ( uRho_avg_new)

          x2 = u_avg_new

          IF (x1 .NE.  x0) THEN

             lambda = ABS((x2 - x1)/(x1 - x0))

          END IF

          denominator = (x2 - x1) - (x1 - x0)

          IF ( ABS(denominator) .LT. 0.1*abs_tol ) EXIT aitken_loop

          aitkenX = x2 - ( (x2 - x1)**2 ) / denominator

          u_avg_new = aitkenX

          IF ( ( ABS(u_avg_guess-u_avg_new)/u_avg_guess < rel_tol ) .OR.        &
               ( ABS(u_avg_guess-u_avg_new) < abs_tol ) ) THEN

             EXIT aitken_loop

          END IF

          u_avg_guess = u_avg_new

          
       END DO aitken_loop

       r_u = u_avg_new * r_qj(2) / ( SQRT( r_qj(2)**2 + r_qj(3)**2 ) )
       r_v = u_avg_new * r_qj(3) / ( SQRT( r_qj(2)**2 + r_qj(3)**2 ) )

    ELSE

       ! velocity components
       IF ( r_qj(1) .GT. eps_sing ) THEN

          r_u = r_qj(2) * inv_qj1
          r_v = r_qj(3) * inv_qj1

       ELSE

          r_u = SQRT(2.0_wp) * r_qj(1) * r_qj(2) / SQRT( r_qj(1)**4 + eps_sing4 )
          r_v = SQRT(2.0_wp) * r_qj(1) * r_qj(3) / SQRT( r_qj(1)**4 + eps_sing4 )

       END IF

    END IF

    p_dyn = 0.5 * r_rho_m * ( r_u**2 + r_v**2 )
    
    ! Richardson number
    IF ( ( r_u**2 + r_v**2 ) .GT. 0.0_wp ) THEN

       r_Ri = r_red_grav * r_h / ( r_u**2 + r_v**2 )

    ELSE

       r_Ri = 0.0_wp

    END IF

    RETURN

  END SUBROUTINE r_phys_var


  SUBROUTINE avg_profiles_mix( h , settling_vel , rho_alphas_avg , u_guess ,    &
       h0 , b , u_rel0 , rho_c , rhom_avg , uRho_avg_new , p_dyn )

    USE geometry_2d, ONLY : z_quad , w_quad
    
    USE geometry_2d, ONLY : calcei , gaulegf

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: h
    REAL(wp), INTENT(IN) :: settling_vel(n_solid)
    REAL(wp), INTENT(IN) :: rho_alphas_avg(n_solid)

    REAL(wp), INTENT(IN) :: u_guess
    REAL(wp), INTENT(IN) :: h0
    REAL(wp), INTENT(IN) :: b
    REAL(wp), INTENT(IN) :: u_rel0
    REAL(wp), INTENT(IN) :: rho_c
    REAL(wp), INTENT(OUT) :: rhom_avg
    REAL(wp), INTENT(OUT) :: uRho_avg_new
    REAL(wp), INTENT(OUT) :: p_dyn

    !> Shear velocity computed from u_guess
    REAL(wp) :: shear_vel

    !>  Rouse numbers for the particle classes
    REAL(wp) :: Rouse_no(n_solid)

    !> array for depth-averaged value of rho*u(z)*C(z)
    REAL(wp) :: rho_u_alphas(n_solid)


    REAL(wp) :: rho_alphas(n_solid)
    REAL(wp) :: alphas(n_solid)

    INTEGER :: i_solid

    REAL(wp) :: normalizing_coeff_u
    
    REAL(wp) :: a(n_solid)

    REAL(wp) :: alphas_exp_avg
    REAL(wp) :: u_log_avg
    
    REAL(wp) :: normalizing_coeff_alpha(n_solid)
    REAL(wp) :: alphas_rel0
    REAL(wp) :: y

    REAL(wp) :: int_def1, int_def2
    
    REAL(wp) :: epsilon_s
    REAL(wp) :: a_coeff

    INTEGER ( kind = 4 ) :: i
    ! REAL(wp) :: x,ei

    REAL(wp) :: z(n_quad)
    REAL(wp) :: w(n_quad)
    REAL(wp) :: int_quad

    REAL(wp) :: rhom_z
    REAL(wp) :: u_z

    REAL(wp) :: z_test

    ! Shear velocity computed from u_guess
    shear_vel = u_guess * SQRT(friction_factor)

    ! Rouse numbers for the particle classes
    DO i_solid=1,n_solid

       IF ( shear_vel .GT. 0.0_wp ) THEN

          Rouse_no(i_solid) = settling_vel(i_solid) / ( vonK * shear_vel )
          
       ELSE
          
          Rouse_no(i_solid) = 0.0_wp
          
       END IF
       
    END DO
    
    ! Quadrature points and weights for the interval [0;h0] 
    z = 0.5_wp * h0 * ( z_quad + 1.0_wp )
    w = 0.5_wp * h0 * w_quad

    epsilon_s = Sc * shear_vel * vonK * ( ( h0/6.0_wp ) + ( k_s / 60.0_wp ) )
    a_coeff = - vonK * shear_vel / epsilon_s

    u_log_avg = ( SUM( w * u_log_profile(b,z) ) + u_log_profile(b,h0)*(h-h0) )  &
         / h

    normalizing_coeff_u = 1.0_wp / u_log_avg
    
    DO i_solid=1,n_solid

       a(i_solid) = a_coeff * Rouse_no(i_solid)

       ! depth-average value of exp(a*x)
       alphas_exp_avg = ( SUM( w * alphas_exp_profile(a(i_solid),z) ) +         &
            alphas_exp_profile(a(i_solid),h0)*(h-h0) ) / h
       
       normalizing_coeff_alpha(i_solid) = 1.0_wp / alphas_exp_avg

       int_quad = SUM( w * ( alphas_exp_profile(a(i_solid),z) *                 &
            u_log_profile(b,z) ) )

       ! integral of alfa_rel_i*u between in the boundary layer
       int_def1 = normalizing_coeff_u * normalizing_coeff_alpha(i_solid) *      &
            int_quad

       ! relative concentration alphas_rel at depth h0 (from the bottom)
       ! alphas_rel is defined as alphas(z)/alphas_avg
       alphas_rel0 = normalizing_coeff_alpha(i_solid) *                         &
            alphas_exp_profile(a(i_solid),h0)

       ! integral of alfa_rel_i*u in the free-stream layer
       int_def2 =  ( h - h0 ) * u_rel0 * alphas_rel0
       
       ! we add the contribution of the integral of the constant region, we
       ! average by dividing by h and we multiply by the density of solid and
       ! average concentration and by u_guess.
       rho_u_alphas(i_solid) = rho_alphas_avg(i_solid) * u_guess *              &
            ( int_def1 + int_def2 ) / h

    END DO

    ! we add the contribution of the gas phase to the depth-averaged mixture 
    ! density. This value should be equal to that used to compute the input
    ! values rhoalphas_avg.
    rhom_avg = rho_c + SUM( ( rho_s - rho_c ) / rho_s * rho_alphas_avg )

    ! we add the contribution of the gas phase to the mixture depth-averaged 
    ! momentum
    uRho_avg_new = ( u_guess*rho_c + SUM((rho_s-rho_c) / rho_s * rho_u_alphas) )


    IF ( z_dyn .GT. 0.0_wp ) THEN
       
       z_test = MIN(z_dyn,h0)
       
       alphas = rho_alphas_avg / rho_s
       
       p_dyn = dynamic_pressure(rho_c, alphas, normalizing_coeff_alpha , a, &
            u_guess, normalizing_coeff_u , b, z_test)

    ELSE

       p_dyn = 0.5_wp * rhom_avg * uRho_avg_new

    END IF
       
  END SUBROUTINE avg_profiles_mix

  !******************************************************************************
  !> \brief Physical variables
  !
  !> This subroutine evaluates from the conservative local variables qj
  !> the local physical variables  (\f$h,u,v,T,\rho_m,red grav,\alpha_s \f$).
  !> \param[in]    c_qj      complex conservative variables 
  !> \param[out]   h         complex-value flow thickness 
  !> \param[out]   u         complex-value flow x-velocity 
  !> \param[out]   v         complex-value flow y-velocity
  !> \param[out]   T         complex-value flow temperature 
  !> \param[out]   rho_m     complex-value flow density
  !> \param[out]   alphas    complex-value solid volume fractions
  !> \param[out]   inv_rhom  complex-value mixture density reciprocal
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !> \date 2019/12/13
  !******************************************************************************

  SUBROUTINE c_phys_var( c_qj , h , u , v , T , rho_m , alphas , alphag ,       &
       inv_rhom , Zs , exc_pore_pres )

    USE COMPLEXIFY
    USE parameters_2d, ONLY : eps_sing , eps_sing4
    IMPLICIT none

    COMPLEX(wp), INTENT(IN) :: c_qj(n_vars)
    COMPLEX(wp), INTENT(OUT) :: h               !< height [m]
    COMPLEX(wp), INTENT(OUT) :: u               !< velocity (x direction) [m s-1]
    COMPLEX(wp), INTENT(OUT) :: v               !< velocity (y direction) [m s-1]
    COMPLEX(wp), INTENT(OUT) :: T               !< temperature [K]
    COMPLEX(wp), INTENT(OUT) :: rho_m           !< mixture density [kg m-3]
    COMPLEX(wp), INTENT(OUT) :: alphas(n_solid) !< sediment volume fractions
    COMPLEX(wp), INTENT(OUT) :: alphag(n_solid) !< sediment volume fractions
    COMPLEX(wp), INTENT(OUT) :: inv_rhom        !< 1/mixture density [kg-1 m3]
    COMPLEX(wp), INTENT(OUT) :: Zs              !< stochastic variable
    COMPLEX(wp), INTENT(OUT) :: exc_pore_pres   !< excess pore pressure

    COMPLEX(wp) :: xs(n_solid)             !< sediment mass fractions
    COMPLEX(wp) :: xg(n_add_gas)           !< additional gas comp. mass fractions
    COMPLEX(wp) :: xs_tot                  !< sum of solid mass fraction
    COMPLEX(wp) :: xl                      !< liquid mass fraction
    COMPLEX(wp) :: xc                      !< carrier phase mass fraction
    COMPLEX(wp) :: sp_heat_c               !< Specific heat of carrier phase
    COMPLEX(wp) :: sp_heat_mix             !< Specific heat of mixture
    COMPLEX(wp) :: sp_gas_const_c          !< Gas constant of carrier phase
    COMPLEX(wp) :: inv_cqj1                !< reciprocal of 1st cons. variable
    COMPLEX(wp) :: inv_rho_c               !< carrier phase density reciprocal
    COMPLEX(wp) :: inv_rho_g(n_add_gas)    !< add. gas density reciprocal

    ! compute solid mass fractions
    IF ( REAL(c_qj(1)) .GT.  EPSILON(1.0_wp) ) THEN

       inv_cqj1 = 1.0_wp / c_qj(1)
       xs(1:n_solid) = c_qj(idx_alfas_first:idx_alfas_last) * inv_cqj1

       xg(1:n_add_gas) = c_qj(idx_addGas_first:idx_addGas_last) * inv_cqj1

       IF ( stoch_transport_flag ) Zs = c_qj(idx_stoch) * inv_cqj1

       IF ( pore_pressure_flag) exc_pore_pres = c_qj(idx_pore) * inv_cqj1    
       
    ELSE

       h = CMPLX(0.0_wp,0.0_wp,wp)
       u = CMPLX(0.0_wp,0.0_wp,wp)
       v = CMPLX(0.0_wp,0.0_wp,wp)
       T = CMPLX(T_ambient,0.0_wp,wp)
       rho_m = CMPLX(rho_a_amb,0.0_wp,wp)
       alphas = CMPLX(0.0_wp,0.0_wp,wp)
       alphag = CMPLX(0.0_wp,0.0_wp,wp)
       inv_rhom = 1.0_wp / rho_m
       Zs = 0.0_wp
       exc_pore_pres = 0.0_wp

       RETURN       

    END IF

    xs_tot = SUM(xs)

    IF ( gas_flag .AND. liquid_flag ) THEN

       ! compute liquid mass fraction
       xl = c_qj(n_vars) * inv_cqj1

       ! compute carrier phase (gas) mass fraction
       xc = 1.0_wp - xs_tot - xl

       sp_heat_c = ( ( xc - SUM( xg(1:n_add_gas) ) ) * sp_heat_a +              &
            DOT_PRODUCT( xg(1:n_add_gas) , sp_heat_g(1:n_add_gas) ) ) / xc

       ! specific heat of the mixutre: mass average of sp. heat pf phases
       sp_heat_mix = DOT_PRODUCT( xs(1:n_solid) , sp_heat_s(1:n_solid) )        &
            + xl * sp_heat_l + xc * sp_heat_c

    ELSE

       ! compute carrier phase (gas or liquid) mass fraction
       xc = 1.0_wp - xs_tot

       IF ( gas_flag ) THEN

          sp_heat_c = ( ( xc - SUM( xg(1:n_add_gas) ) ) * sp_heat_a +           &
               DOT_PRODUCT( xg(1:n_add_gas) , sp_heat_g(1:n_add_gas) ) ) / xc

       ELSE

          sp_heat_c = CMPLX(sp_heat_l,0.0_wp,wp)

       END IF

       ! specific heaf of the mixutre: mass average of sp. heat pf phases
       sp_heat_mix = DOT_PRODUCT( xs(1:n_solid) , sp_heat_s(1:n_solid) )        &
            + xc * sp_heat_c

    END IF

    ! compute temperature from energy
    IF ( REAL(c_qj(1)) .GT. eps_sing ) THEN

       IF ( energy_flag ) THEN

          T = ( c_qj(4) - 0.5_wp * ( c_qj(2)**2 + c_qj(3)**2 ) * inv_cqj1 ) /   &
               ( c_qj(1) * sp_heat_mix ) 

       ELSE

          T = c_qj(4) / ( c_qj(1) * sp_heat_mix ) 

       END IF

       IF ( REAL(T) .LE. 0.0_wp ) T = CMPLX(T_ambient,0.0_wp,wp)

    ELSE

       T = CMPLX(T_ambient,0.0_wp,wp)

    END IF

    IF ( gas_flag ) THEN

       ! carrier phase is gas
       sp_gas_const_c = ( ( xc - SUM( xg(1:n_add_gas) ) ) * sp_gas_const_a      &
            + DOT_PRODUCT( xg(1:n_add_gas) , sp_gas_const_g(1:n_add_gas) ) )    &
            / xc

       inv_rho_c = sp_gas_const_c * T * inv_pres

       inv_rho_g(1:n_add_gas) = sp_gas_const_g(1:n_add_gas) * T * inv_pres

    ELSE

       inv_rho_c = CMPLX(inv_rho_l,0.0_wp,wp)

    END IF

    inv_rhom = DOT_PRODUCT( xs(1:n_solid) , c_inv_rho_s(1:n_solid) )            &
         + xc * inv_rho_c

    IF ( gas_flag .AND. liquid_flag ) inv_rhom = inv_rhom + xl * inv_rho_l

    rho_m = 1.0_wp / inv_rhom

    ! convert from mass fraction to volume fraction
    alphas(1:n_solid) = rho_m * xs(1:n_solid) * c_inv_rho_s(1:n_solid)

    ! convert from mass fraction to volume fraction
    alphag(1:n_add_gas) = rho_m * xg(1:n_add_gas) * inv_rho_g(1:n_add_gas)

    h = c_qj(1) * inv_rhom

    ! velocity components
    IF ( REAL( c_qj(1) ) .GT. eps_sing ) THEN

       u = c_qj(2) * inv_cqj1
       v = c_qj(3) * inv_cqj1

    ELSE

       u = SQRT(2.0_wp) * c_qj(1) * c_qj(2) / SQRT( c_qj(1)**4 + eps_sing4 )
       v = SQRT(2.0_wp) * c_qj(1) * c_qj(3) / SQRT( c_qj(1)**4 + eps_sing4 )

    END IF

    RETURN

  END SUBROUTINE c_phys_var


  !******************************************************************************
  !> \brief Mixture variables
  !
  !> This subroutine evaluates from the physical real-value local variables qpj, 
  !> some mixture variable.
  !> \param[in]    qpj          real-valued physical variables 
  !> \param[out]   r_Ri         real-valued Richardson number 
  !> \param[out]   r_rho_m      real-valued mixture density 
  !> \param[out]   r_rho_c      real-valued carrier phase density 
  !> \param[out]   r_red_grav   real-valued reduced gravity
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !> \date 10/10/2019
  !******************************************************************************

  SUBROUTINE mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav,sp_heat_flag,         &
       r_sp_heat_c,r_sp_heat_mix)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qpj(n_vars+2) !< real-value physical variables
    REAL(wp), INTENT(OUT) :: r_Ri         !< real-value Richardson number
    REAL(wp), INTENT(OUT) :: r_rho_m      !< real-value mixture density [kg/m3]
    REAL(wp), INTENT(OUT) :: r_rho_c !< real-value carrier phase density [kg/m3]
    REAL(wp), INTENT(OUT) :: r_red_grav   !< real-value reduced gravity
    LOGICAL, INTENT(IN) :: sp_heat_flag
    REAL(wp), INTENT(OUT) :: r_sp_heat_c
    REAL(wp), INTENT(OUT) :: r_sp_heat_mix


    REAL(wp) :: r_u                       !< real-value x-velocity
    REAL(wp) :: r_v                       !< real-value y-velocity
    REAL(wp) :: r_h                       !< real-value flow thickness
    REAL(wp) :: r_alphas(n_solid)         !< real-value solid volume fractions
    REAL(wp) :: r_alphag(n_add_gas)       !< real-value add.gas volume fractions
    REAL(wp) :: r_rho_a                   !< real-value atm.gas density
    REAL(wp) :: r_rho_g(n_add_gas)        !< real-value add.gas densities
    REAL(wp) :: r_T                       !< real-value temperature [K]
    REAL(wp) :: r_alphal                  !< real-value liquid volume fraction
    REAL(wp) :: r_alphac

    REAL(wp) :: alphas_tot                !< total solid fraction

    REAL(wp) :: r_inv_rhom

    r_h = qpj(1)

    IF ( qpj(1) .LE. EPSILON(1.0_wp) ) THEN

       r_red_grav = 0.0_wp
       r_rho_m = rho_a_amb
       r_rho_c = rho_a_amb

       IF ( liquid_flag ) THEN

          r_sp_heat_c = sp_heat_l
          r_sp_heat_mix = sp_heat_l

       ELSE

          r_sp_heat_c = sp_heat_a
          r_sp_heat_mix = sp_heat_a


       END IF

       r_Ri = 0.0_wp
       
       RETURN

    END IF

    r_u = qpj(idx_u)
    r_v = qpj(idx_v)
    r_T = qpj(4)

    IF ( alpha_flag ) THEN

       r_alphas(1:n_solid) = qpj(idx_alfas_first:idx_alfas_last)
       r_alphag(1:n_add_gas) = qpj(idx_addGas_first:idx_addGas_last)

    ELSE

       r_alphas(1:n_solid) = qpj(idx_alfas_first:idx_alfas_last) / qpj(1)
       r_alphag(1:n_add_gas) = qpj(idx_addGas_first:idx_addGas_first) / qpj(1)

    END IF

    alphas_tot = SUM(r_alphas)

    r_alphal = 0.0_wp

    IF ( gas_flag .AND. liquid_flag ) THEN

       IF ( alpha_flag ) THEN

          r_alphal = qpj(n_vars)

       ELSE

          r_alphal = qpj(n_vars) / qpj(1)

       END IF

    END IF

    ! carrier phase volume fraction
    r_alphac = 1.0_wp - alphas_tot - r_alphal

    IF ( gas_flag ) THEN

       ! continuous phase is gas
       r_rho_a =  pres / ( sp_gas_const_a * r_T )
       r_rho_g(1:n_add_gas) = pres / ( sp_gas_const_g(1:n_add_gas) * r_T )

       r_rho_c = ( ( 1.0_wp - r_alphal - alphas_tot - SUM(r_alphag) ) * r_rho_a &
            + DOT_PRODUCT( r_alphag(1:n_add_gas) , r_rho_g(1:n_add_gas) ) )     &
            / ( 1.0_wp - r_alphal - alphas_tot )

    ELSE

       ! continuous phase is liquid
       r_rho_c = rho_l

    END IF

    IF ( gas_flag .AND. liquid_flag ) THEN

       ! density of mixture of carrier (gas), liquid and solids
       r_rho_m = ( 1.0_wp - alphas_tot - r_alphal ) * r_rho_c                   &
            + DOT_PRODUCT( r_alphas , rho_s ) + r_alphal * rho_l

    ELSE

       ! density of mixture of carrier phase and solids
       r_rho_m = ( 1.0_wp - alphas_tot ) * r_rho_c + DOT_PRODUCT( r_alphas ,    &
            rho_s ) 

    END IF

    r_inv_rhom = 1.0_wp / r_rho_m

    ! reduced gravity
    r_red_grav = ( r_rho_m - rho_a_amb ) / r_rho_m * grav

    ! Richardson number
    IF ( ( r_u**2 + r_v**2 ) .GT. 0.0_wp ) THEN

       r_Ri = MIN(1.E15_wp,r_red_grav * r_h / ( r_u**2 + r_v**2 ))

    ELSE

       r_Ri = 0.0_wp

    END IF

    IF ( sp_heat_flag ) THEN

       CALL eval_sp_heat( r_alphal , r_alphas , r_alphag, r_rho_g , r_inv_rhom ,&
            r_sp_heat_c , r_sp_heat_mix )

    END IF

    RETURN

  END SUBROUTINE mixt_var

  !******************************************************************************
  !> \brief Specific heat
  !
  !> This subroutine evaluates the specific heat of the carrier phase and of the
  !> mixture.
  !> \param[in]    r_alphal        real-value liquid volume fraction
  !> \param[in]    r_alphas        real-value solid volume fraction
  !> \param[in]    r_alphag        real-value add.gas volume fraction
  !> \param[in]    rho_g           real-value gas density
  !> \param[in]    r_inv_rhom      real-value mixture density reciprocal
  !> \param[out]   r_sp_heat_c     real-valued carrier phase specific heat
  !> \param[out]   r_sp_heat_mix   real-valued mixture specific heat
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !> \date 2021/07/09
  !******************************************************************************

  SUBROUTINE eval_sp_heat( r_alphal , r_alphas , r_alphag , rho_g ,r_inv_rhom , &
       r_sp_heat_c , r_sp_heat_mix )

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: r_alphal          !< real-value liquid volume fraction
    REAL(wp), INTENT(IN) :: r_alphas(n_solid) !< real-value solid volume fractions
    REAL(wp), INTENT(IN) :: r_alphag(n_add_gas)!< real-value add.gas volume fractions
    ! REAL(wp), INTENT(IN) :: r_alphac    !< real-value carrier phase volume fraction

    REAL(wp), INTENT(IN) :: rho_g(n_add_gas)!< real-value add.gas densities [kg/m3]
    ! REAL(wp), INTENT(IN) :: r_rho_c   !< real-value carrier phase density [kg/m3]
    REAL(wp), INTENT(IN) :: r_inv_rhom        !< real-value mixture density [kg/m3]

    REAL(wp), INTENT(OUT) :: r_sp_heat_c
    REAL(wp), INTENT(OUT) :: r_sp_heat_mix

    REAL(wp) :: r_xl              !< real-value liquid mass fraction
    REAL(wp) :: r_xc              !< real-value carrier phase mass fraction

    REAL(wp) :: r_xs(n_solid)     !< real-value solid mass fractions
    REAL(wp) :: r_xg(n_add_gas)   !< real-value add.gas mass fractions    

    IF ( gas_flag .AND. liquid_flag ) THEN

       ! liquid mass fraction
       r_xl = r_alphal * rho_l * r_inv_rhom

       ! solid mass fractions
       r_xs(1:n_solid) = r_alphas(1:n_solid) * rho_s(1:n_solid) * r_inv_rhom

       ! additional gas mass fractions
       r_xg(1:n_add_gas) = r_alphag(1:n_add_gas) * rho_g(1:n_add_gas) * r_inv_rhom

       ! carrier (gas) mass fraction
       r_xc = 1.0_wp - ( r_xl + SUM(r_xs(1:n_solid) ) )

       ! specific heat of gas (mass. avg. of sp.heat of gas components)

       IF ( r_xc .GT. EPSILON(1.0_wp) ) THEN
        
          r_sp_heat_c = ( ( r_xc - SUM( r_xg(1:n_add_gas) ) ) * sp_heat_a +        &
               DOT_PRODUCT( r_xg(1:n_add_gas) , sp_heat_g(1:n_add_gas) ) ) / r_xc

       ELSE

          r_sp_heat_c = sp_heat_a

       END IF
          
       ! mass averaged mixture specific heat
       r_sp_heat_mix =  DOT_PRODUCT( r_xs , sp_heat_s ) + r_xl * sp_heat_l      &
            + r_xc * r_sp_heat_c

    ELSE

       ! solid mass fractions
       r_xs(1:n_solid) = r_alphas(1:n_solid) * rho_s(1:n_solid) * r_inv_rhom

       ! additional gas mass fractions
       r_xg(1:n_add_gas) = r_alphag(1:n_add_gas) * rho_g(1:n_add_gas)           &
            * r_inv_rhom

       ! carrier (gas or liquid) mass fraction
       r_xc = 1.0_wp - SUM( r_xs(1:n_solid) )

       r_sp_heat_c = 0.0_wp

       IF ( gas_flag ) THEN

          IF ( r_xc .GT. EPSILON(1.0_wp) ) THEN
 
             r_sp_heat_c = ( ( r_xc - SUM( r_xg(1:n_add_gas) ) ) * sp_heat_a +     &
                  DOT_PRODUCT( r_xg(1:n_add_gas) , sp_heat_g(1:n_add_gas) ) ) / r_xc

          ELSE

             r_sp_heat_c = sp_heat_a

          END IF
             
       ELSE

          r_sp_heat_c = sp_heat_l

       END IF

       ! mass averaged mixture specific heat
       r_sp_heat_mix =  DOT_PRODUCT( r_xs , sp_heat_s ) + r_xc * r_sp_heat_c

    END IF

    RETURN

  END SUBROUTINE eval_sp_heat

  !******************************************************************************
  !> \brief Conservative to physical variables
  !
  !> This subroutine evaluates from the conservative variables qc the 
  !> array of physical variables qp:\n
  !> - qp(1) = \f$ h \f$
  !> - qp(2) = \f$ hu \f$
  !> - qp(3) = \f$ hv \f$
  !> - qp(4) = \f$ T \f$
  !> - qp(idx_alfas_first:idx_alfas_last) = \f$ alphas(1:n_solid) \f$
  !> - qp(idx_addGas_first:idx_addGas_last) = \f$ alphas(1:n_add_gas) \f$
  !> - qp(n_vars) = \f$ alphal \f$
  !> - qp(idx_u) = \f$ u \f$
  !> - qp(idx_v) = \f$ v \f$
  !> .
  !> The physical variables are those used for the linear reconstruction at the
  !> cell interfaces. Limiters are applied to the reconstructed slopes.
  !> \param[in]     qc     local conservative variables 
  !> \param[out]    qp     local physical variables  
  !> \param[out]    p_dyn  local dynamic pressure
  !
  !> \date 2019/11/11
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE qc_to_qp(qc,qp,p_dyn)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qc(n_vars)
    REAL(wp), INTENT(OUT) :: qp(n_vars+2)
    REAL(wp), INTENT(OUT) :: p_dyn

    REAL(wp) :: r_h               !< real-value flow thickness
    REAL(wp) :: r_u               !< real-value x-velocity
    REAL(wp) :: r_v               !< real-value y-velocity
    REAL(wp) :: r_alphas(n_solid) !< real-value solid volume fractions
    REAL(wp) :: r_rho_m           !< real-value mixture density [kg/m3]
    REAL(wp) :: r_T               !< real-value temperature [K]
    REAL(wp) :: r_alphal          !< real-value liquid volume fraction
    REAL(wp) :: r_alphag(n_add_gas) !< real-value add. gas volume fractions
    REAL(wp) :: r_red_grav
    REAL(wp) :: r_Zs             !< real-value stochastic variable
    REAL(wp) :: r_exc_pore_pres  !< real-value pore pressure

    CALL r_phys_var( qc , r_h , r_u , r_v , r_alphas , r_rho_m , r_T ,          &
         r_alphal , r_alphag , r_red_grav , p_dyn , r_Zs , r_exc_pore_pres )

    qp(1) = r_h

    qp(2) = r_h*r_u
    qp(3) = r_h*r_v

    qp(4) = r_T

    IF ( alpha_flag ) THEN

       qp(idx_alfas_first:idx_alfas_last) = r_alphas(1:n_solid)
       qp(idx_addGas_first:idx_addGas_last) = r_alphag(1:n_add_gas)
       IF ( gas_flag .AND. liquid_flag ) qp(n_vars) = r_alphal

    ELSE

       qp(idx_alfas_first:idx_alfas_last) = r_alphas(1:n_solid) * r_h
       qp(idx_addGas_first:idx_addGas_last) = r_alphag(1:n_add_gas) * r_h
       IF ( gas_flag .AND. liquid_flag ) qp(n_vars) = r_alphal * r_h

    END IF

    IF ( stoch_transport_flag) qp(idx_stoch) = r_Zs

    IF ( pore_pressure_flag) qp(idx_pore) = r_exc_pore_pres
    
    qp(idx_u) = r_u
    qp(idx_v) = r_v

    RETURN

  END SUBROUTINE qc_to_qp

  !******************************************************************************
  !> \brief Physical to conservative variables
  !
  !> This subroutine evaluates the conservative real_value variables qc from the 
  !> array of real_valued physical variables qp:\n
  !> - qp(1) = \f$ h \f$
  !> - qp(2) = \f$ h*u \f$
  !> - qp(3) = \f$ h*v \f$
  !> - qp(4) = \f$ T \f$
  !> - qp(idx_alfas_first:idx_alfas_last) = \f$ alphas(1:n_s) \f$
  !> - qp(idx_addGas_first:idx_addGas_last) = \f$ alphas(1:n_g) \f$
  !> - qp(n_vars) = \f$ alphal \f$
  !> - qp(idx_u) = \f$ u \f$
  !> - qp(idx_v) = \f$ v \f$
  !> .
  !> \param[in]    qp      physical variables  
  !> \param[in]    B       local topography
  !> \param[out]   qc      conservative variables
  !
  !> \date 2019/11/18
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE qp_to_qc(qp,qc)

    USE geometry_2d, ONLY : z_quad, w_quad
    
    USE parameters_2d, ONLY : vertical_profiles_flag
    USE geometry_2d, ONLY : gaulegf
    USE geometry_2d, ONLY : lambertw,lambertw0,lambertwm1

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qp(n_vars+2)
    REAL(wp), INTENT(OUT) :: qc(n_vars)

    REAL(wp) :: r_sp_heat_mix
    REAL(wp) :: sum_sl

    REAL(wp) :: r_u               !< real-value x-velocity
    REAL(wp) :: r_v               !< real-value y-velocity
    REAL(wp) :: r_h               !< real-value flow thickness
    REAL(wp) :: r_hu              !< real-value volumetric x-flow
    REAL(wp) :: r_hv              !< real-value volumetric y-flow
    REAL(wp) :: r_alphas(n_solid) !< real-value solid volume fractions
    REAL(wp) :: r_alphag(n_add_gas)!< real-value add.gas volume fractions
    REAL(wp) :: r_xl              !< real-value liquid mass fraction
    REAL(wp) :: r_xc              !< real-value carrier phase mass fraction
    REAL(wp) :: r_T               !< real-value temperature [K]
    REAL(wp) :: r_alphal          !< real-value liquid volume fraction
    REAL(wp) :: r_alphac          !< real-value carrier phase volume fraction
    REAL(wp) :: r_rho_m           !< real-value mixture density [kg/m3]
    REAL(wp) :: r_rho_c           !< real-value carrier phase density [kg/m3]
    REAL(wp) :: r_rho_a           !< real-value atm.gas density [kg/m3]
    REAL(wp) :: r_rho_g(n_add_gas)!< real-value add.gas densities [kg/m3]
    REAL(wp) :: r_xs(n_solid)     !< real-value solid mass fractions
    REAL(wp) :: r_xg(n_add_gas)   !< real-value add.gas mass fractions

    REAL(wp) :: r_Zs              !< real-value stochastic variable
    REAL(wp) :: r_exc_pore_pres   !< real-value pore pressure
    
    REAL(wp) :: r_alphas_rhos(n_solid)
    REAL(wp) :: r_alphag_rhog(n_add_gas)
    REAL(wp) :: alphas_tot

    REAL(wp) :: r_sp_heat_c

    REAL(wp) :: r_inv_rhom

    REAL(wp) :: rho_u_alphas(n_solid)
    REAL(wp) :: uRho_avg

    REAL(wp) :: u_rel0
    REAL(wp) :: normalizing_coeff_u
    REAL(wp) :: shear_vel

    REAL(wp) :: a , b , c , d

    REAL(wp) :: z(n_quad) , w(n_quad)

    REAL(wp) :: inv_kin_visc
    REAL(wp) :: settling_vel

    REAL(wp) :: Rouse_no(n_solid)

    INTEGER :: i_solid

    REAL(wp) :: r_w
    REAL(wp) :: mod_vel
    REAL(wp) :: log_term_h0
    REAL(wp) :: alphas_exp_avg
    REAL(wp) :: int_quad
    REAL(wp) :: int_def1 , int_def2
    REAL(wp) :: h_rel , h0_rel , h0
    REAL(wp) :: h0_rel_1
    REAL(wp) :: h0_rel_2

    REAL(wp) :: epsilon_s
    REAL(wp) :: a_coeff
    REAL(wp) :: alphas_rel0
    real(wp) :: normalizing_coeff_alpha

    REAL(wp) :: u_log_avg

    REAL(wp) :: dyn_visc_c

    r_h = qp(1)

    IF ( r_h .GT. EPSILON(1.0_wp) ) THEN

       r_hu = qp(2)
       r_hv = qp(3)

       r_u = qp(idx_u)
       r_v = qp(idx_v)

    ELSE

       qc(1:n_vars) = 0.0_wp
       RETURN

    END IF

    r_T  = qp(4)

    IF ( alpha_flag ) THEN

       r_alphas(1:n_solid) = qp(idx_alfas_first:idx_alfas_last)
       r_alphag(1:n_add_gas) = qp(idx_addGas_first:idx_addGas_last)

    ELSE

       r_alphas(1:n_solid) = qp(idx_alfas_first:idx_alfas_last) / qp(1)
       r_alphag(1:n_add_gas) = qp(idx_addGas_first:idx_addGas_last) / qp(1)

    END IF

    alphas_tot = SUM(r_alphas)

    r_alphas_rhos(1:n_solid) = r_alphas(1:n_solid) * rho_s(1:n_solid)

    r_alphal = 0.0_wp

    IF ( gas_flag .AND. liquid_flag ) THEN

       IF ( alpha_flag ) THEN

          r_alphal = qp(n_vars)

       ELSE

          r_alphal = qp(n_vars) / qp(1)

       END IF

    END IF

    IF ( gas_flag ) THEN

       ! continuous phase is air
       r_rho_a =  pres / ( sp_gas_const_a * r_T )
       r_rho_g(1:n_add_gas) = pres / ( sp_gas_const_g(1:n_add_gas) * r_T )
       r_alphag_rhog(1:n_add_gas) = r_alphag(1:n_add_gas) * r_rho_g(1:n_add_gas)

       r_rho_c = ( ( 1.0_wp - r_alphal - alphas_tot - SUM(r_alphag) ) * r_rho_a &
            + SUM( r_alphag_rhog(1:n_add_gas) ) ) /                             &
            ( 1.0_wp - r_alphal - alphas_tot )

    ELSE

       ! carrier phase is liquid
       r_rho_c = rho_l

    END IF

    IF ( gas_flag .AND. liquid_flag ) THEN

       ! check and correction on dispersed phases volume fractions
       IF ( ( alphas_tot + r_alphal ) .GT. 1.0_wp ) THEN

          sum_sl = alphas_tot + r_alphal
          r_alphas(1:n_solid) = r_alphas(1:n_solid) / sum_sl
          r_alphal = r_alphal / sum_sl

       ELSEIF ( ( alphas_tot + r_alphal ) .LT. 0.0_wp ) THEN

          r_alphas(1:n_solid) = 0.0_wp
          r_alphal = 0.0_wp

       END IF

       ! carrier phase volume fraction
       r_alphac = 1.0_wp - alphas_tot - r_alphal

       ! volume averaged mixture density: carrier (gas) + solids + liquid
       r_rho_m = r_alphac * r_rho_c + SUM( r_alphas_rhos(1:n_solid) )           &
            + r_alphal * rho_l

       r_inv_rhom = 1.0_wp / r_rho_m

       ! liquid mass fraction
       r_xl = r_alphal * rho_l * r_inv_rhom

       ! solid mass fractions
       r_xs(1:n_solid) = r_alphas_rhos(1:n_solid) * r_inv_rhom

       ! additional gas mass fractions
       r_xg(1:n_add_gas) = r_alphag_rhog(1:n_add_gas) * r_inv_rhom

       ! carrier (gas) mass fraction
       r_xc = r_alphac * r_rho_c * r_inv_rhom

       ! specific heat of gas (mass. avg. of sp.heat of gas components)
       IF ( r_xc .GT. EPSILON(1.0_wp) ) THEN

          r_sp_heat_c = ( ( r_xc - SUM( r_xg(1:n_add_gas) ) ) * sp_heat_a +        &
               DOT_PRODUCT( r_xg(1:n_add_gas) , sp_heat_g(1:n_add_gas) ) ) / r_xc

       ELSE

          r_sp_heat_c = sp_heat_a

       END IF
          
       ! mass averaged mixture specific heat
       r_sp_heat_mix =  DOT_PRODUCT( r_xs , sp_heat_s ) + r_xl * sp_heat_l      &
            + r_xc * r_sp_heat_c

    ELSE

       ! mixture of carrier phase ( gas or liquid ) and solid

       ! check and corrections on dispersed phases
       IF ( alphas_tot .GT. 1.0_wp ) THEN

          r_alphas(1:n_solid) = r_alphas(1:n_solid) / alphas_tot

       ELSEIF ( alphas_tot .LT. 0.0_wp ) THEN

          r_alphas(1:n_solid) = 0.0_wp

       END IF

       ! carrier (gas or liquid) volume fraction
       r_alphac = 1.0_wp - alphas_tot 

       ! volume averaged mixture density: carrier (gas or liquid) + solids
       r_rho_m = r_alphac * r_rho_c + SUM( r_alphas_rhos(1:n_solid) )

       r_inv_rhom = 1.0_wp / r_rho_m

       ! solid mass fractions
       r_xs(1:n_solid) = r_alphas_rhos(1:n_solid) * r_inv_rhom

       ! additional gas mass fractions
       r_xg(1:n_add_gas) = r_alphag_rhog(1:n_add_gas) * r_inv_rhom

       ! carrier (gas or liquid) mass fraction
       r_xc = r_alphac * r_rho_c * r_inv_rhom

       IF ( gas_flag ) THEN

          IF ( r_xc .GT. EPSILON(1.0_wp) ) THEN
          
             r_sp_heat_c = ( ( r_xc - SUM( r_xg(1:n_add_gas) ) ) * sp_heat_a +     &
                  DOT_PRODUCT( r_xg(1:n_add_gas) , sp_heat_g(1:n_add_gas) ) ) / r_xc

          ELSE

             r_sp_heat_c = sp_heat_a

          END IF
             
       ELSE

          r_sp_heat_c = sp_heat_l

       END IF

       ! mass averaged mixture specific heat
       r_sp_heat_mix =  DOT_PRODUCT( r_xs , sp_heat_s ) + r_xc * r_sp_heat_c

    END IF

    IF ( stoch_transport_flag) r_Zs = qp(idx_stoch)

    IF ( pore_pressure_flag ) r_exc_pore_pres = qp(idx_pore)
    
    qc(1) = r_rho_m * r_h

    IF ( vertical_profiles_flag ) THEN

       r_w = 0.0_wp

       mod_vel = SQRT( r_u**2 + r_v**2 + r_w**2 )

       shear_vel = SQRT( friction_factor ) * mod_vel

       IF ( gas_flag .AND. sutherland_flag ) THEN
          
          dyn_visc_c = muRef_Suth * ( r_T / Tref_Suth )**1.5_wp *               &
               ( Tref_Suth + S_mu ) / ( r_T + S_mu )

          kin_visc_c = dyn_visc_c / r_rho_c

       END IF
                     
       ! Viscosity read from input file [m2 s-1]
       inv_kin_visc = 1.0_wp / kin_visc_c
       
       DO i_solid=1,n_solid

          settling_vel = settling_velocity( diam_s(i_solid) , rho_s(i_solid) ,  &
               r_rho_c , inv_kin_visc )

          IF ( shear_vel .GT. 0.0_wp ) THEN

             Rouse_no(i_solid) = settling_vel / ( vonK * shear_vel )

          ELSE

             Rouse_no(i_solid) = 0.0_wp

          END IF

       END DO

       ! The profile parameters depend on h/k_s, not on the absolute value of h.
       h_rel = r_h / k_s

       IF ( h_rel .GT. H_crit_rel ) THEN

          ! we search for h0_rel such that the average integral between 0 and
          ! h_rel is equal to 1
          ! For h_rel > H_crit_rel this integral is the sum of two pieces:
          ! integral between 0 and h0_rel of the log profile
          ! integral between h0_rel and h_rel of the costant profile

          a = h_rel * vonK / SQRT(friction_factor)
          b = 1.0_wp / 30.0_wp + h_rel
          c = 30.0_wp

          ! solve b*log(c*z+1)-z=a for z
          d = a/b - 1.0_wp / (b*c)

          h0_rel_1 = -b*lambertw0( -EXP(d)/(b*c) ) - 1.0_wp / c
          h0_rel_2 = -b*lambertwm1( -EXP(d)/(b*c) ) - 1.0_wp / c
          h0_rel = MIN( h0_rel_1 , h0_rel_2)

       ELSE

          ! when h_rel <= H_crit_rel we have only the log profile and we have to
          ! rescale it in order to have the integral between o and h_rel equal to
          ! 1
          h0_rel = h_rel

       END IF

       h0 = h0_rel * k_s

       b = 30.0_wp / k_s

       log_term_h0 = u_log_profile(b,h0)**2
       
       epsilon_s = Sc * shear_vel * vonK * (( h0 / 6.0_wp ) + ( k_s / 60.0_wp ))
       a_coeff = - vonK * shear_vel / epsilon_s
       
       z = 0.5_wp * h0 * ( z_quad + 1.0_wp )
       w = 0.5_wp * h0 * w_quad

       u_log_avg = ( SUM( w * u_log_profile(b,z) ) + u_log_profile(b,h0) *      &
            ( r_h -h0 ) ) / r_h

       normalizing_coeff_u = 1.0_wp / u_log_avg

       ! relative velocty at h0
       u_rel0 = normalizing_coeff_u * u_log_profile(b,h0)
       
       DO i_solid = 1,n_solid

          a = a_coeff * Rouse_no(i_solid)

          alphas_exp_avg = ( SUM( w * alphas_exp_profile(a,z) ) +               &
               alphas_exp_profile(a,h0)*(r_h-h0) ) / r_h

          normalizing_coeff_alpha = 1.0_wp / alphas_exp_avg

          int_quad = SUM( w * ( alphas_exp_profile(a,z) * u_log_profile(b,z) ) )

          ! integral of alfa_i*u between in the boundary layer
          int_def1 = ( mod_vel * normalizing_coeff_u ) *                        &
               ( r_alphas(i_solid) * normalizing_coeff_alpha ) * int_quad

          ! relative concentration alphas_rel at h0
          alphas_rel0 = normalizing_coeff_alpha * alphas_exp_profile(a,h0)

          ! integral of alfa_i*u in the free-stream layer
          int_def2 =  ( r_h - h0 ) * ( mod_vel * u_rel0 ) *                     &
               ( alphas_rel0 * r_alphas(i_solid) ) 

          ! we add the contribution of the integral of the constant region, we
          ! average by dividing by h and we multiply by the density of solid and
          ! average concentration and by u_guess.
          rho_u_alphas(i_solid) = rho_s(i_solid) * ( int_def1 + int_def2 ) / r_h

       END DO

       ! we add the contribution of the gas phase to the mixture depth-averaged 
       ! momentum
       uRho_avg = ( mod_vel * r_rho_c + SUM((rho_s - r_rho_c) / rho_s *         &
            rho_u_alphas) )

       qc(2) = r_h * uRho_avg * r_u / mod_vel
       qc(3) = r_h * uRho_avg * r_v / mod_vel

    ELSE

       qc(2) = r_rho_m * r_hu
       qc(3) = r_rho_m * r_hv

    END IF

    IF ( energy_flag ) THEN

       IF ( r_h .GT. 0.0_wp ) THEN

          ! total energy (internal and kinetic)
          qc(4) = r_h * r_rho_m * ( r_sp_heat_mix * r_T                         &
               + 0.5_wp * ( r_u**2 + r_v**2 ) )

       ELSE

          qc(4) = 0.0_wp

       END IF

    ELSE

       ! internal energy
       qc(4) = r_h * r_rho_m * r_sp_heat_mix * r_T 

    END IF

    qc(idx_alfas_first:idx_alfas_last) = r_xs * qc(1)
    qc(idx_addGas_first:idx_addGas_last) = r_xg * qc(1)
    
    IF ( stoch_transport_flag ) qc(idx_stoch) = r_Zs * qc(1)
    
    IF ( pore_pressure_flag ) qc(idx_pore) = r_exc_pore_pres * qc(1)
        
    IF ( gas_flag .AND. liquid_flag ) qc(n_vars) = r_xl * qc(1)

    RETURN

  END SUBROUTINE qp_to_qc

  !******************************************************************************
  !> \brief Additional Physical variables
  !
  !> This subroutine evaluates from the physical local variables qpj, the two
  !> additional local variables qp2j = (h+B,u,v). 
  !> \param[in]    qpj    real-valued physical variables 
  !> \param[in]    Bj     real-valued local topography 
  !> \param[out]   qp2j   real-valued physical variables 
  !> @author 
  !> Mattia de' Michieli Vitturi
  !> \date 10/10/2019
  !******************************************************************************

  SUBROUTINE qp_to_qp2(qpj,Bj,qp2j)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qpj(n_vars+2)
    REAL(wp), INTENT(IN) :: Bj
    REAL(wp), INTENT(OUT) :: qp2j(3)

    qp2j(1) = qpj(1) + Bj

    IF ( qpj(1) .LE. 0.0_wp ) THEN

       qp2j(2) = 0.0_wp
       qp2j(3) = 0.0_wp

    ELSE

       qp2j(2) = qpj(2)/qpj(1)
       qp2j(3) = qpj(3)/qpj(1)

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

  SUBROUTINE eval_local_speeds_x(qpj,grav_coeff,vel_min,vel_max)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qpj(n_vars+2)
    REAL(wp), INTENT(IN) :: grav_coeff

    REAL(wp), INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)

    REAL(wp) :: r_h          !< real-value flow thickness [m]
    REAL(wp) :: r_u          !< real-value x-velocity [m s-1]
    REAL(wp) :: r_v          !< real-value y-velocity [m s-1]
    REAL(wp) :: r_Ri         !< real-value Richardson number
    REAL(wp) :: r_rho_m      !< real-value mixture density [kg m-3]
    REAL(wp) :: r_rho_c      !< real-value carrier phase density [kg m-3]
    REAL(wp) :: r_red_grav   !< real-value reduced gravity [m s-2]
    REAL(wp) :: r_celerity
    LOGICAL :: sp_heat_flag
    REAL(wp) :: r_sp_heat_c
    REAL(wp) :: r_sp_heat_mix

    sp_heat_flag = .FALSE.

    CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav,sp_heat_flag,r_sp_heat_c, &
         r_sp_heat_mix)

    r_h = qpj(1)
    r_u = qpj(idx_u)
    r_v = qpj(idx_v)

    IF ( r_red_grav * r_h .LT. 0.0_wp ) THEN

       vel_min(1:n_eqns) = r_u
       vel_max(1:n_eqns) = r_u

    ELSE

       r_celerity = SQRT( r_red_grav * r_h * grav_coeff )
       vel_min(1:n_eqns) = r_u - r_celerity
       vel_max(1:n_eqns) = r_u + r_celerity

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

  SUBROUTINE eval_local_speeds_y(qpj,grav_coeff,vel_min,vel_max)

    IMPLICIT none

    REAL(wp), INTENT(IN)  :: qpj(n_vars+2)
    REAL(wp), INTENT(IN) :: grav_coeff
    REAL(wp), INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)

    REAL(wp) :: r_h          !< real-value flow thickness
    REAL(wp) :: r_u          !< real-value x-velocity
    REAL(wp) :: r_v          !< real-value y-velocity
    REAL(wp) :: r_Ri         !< real-value Richardson number
    REAL(wp) :: r_rho_m      !< real-value mixture density [kg/m3]
    REAL(wp) :: r_rho_c      !< real-value carrier phase density [kg/m3]
    REAL(wp) :: r_red_grav   !< real-value reduced gravity
    REAL(wp) :: r_celerity
    LOGICAL :: sp_heat_flag
    REAL(wp) :: r_sp_heat_c
    REAL(wp) :: r_sp_heat_mix

    sp_heat_flag = .FALSE.

    CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav,sp_heat_flag,r_sp_heat_c, &
         r_sp_heat_mix)

    r_h = qpj(1)
    r_u = qpj(idx_u)
    r_v = qpj(idx_v)

    IF ( r_red_grav * r_h .LT. 0.0_wp ) THEN

       vel_min(1:n_eqns) = r_v
       vel_max(1:n_eqns) = r_v

    ELSE

       r_celerity = SQRT( grav_coeff * r_red_grav * r_h )
       vel_min(1:n_eqns) = r_v - r_celerity
       vel_max(1:n_eqns) = r_v + r_celerity

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
  !> \param[in]     dir      direction of the flux (1=x,2=y)
  !> \param[out]    flux     real  fluxes    
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE eval_fluxes(qcj,qpj,B_prime_x,B_prime_y,grav_coeff,dir,flux)

    USE parameters_2d, ONLY : vertical_profiles_flag

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qcj(n_vars)
    REAL(wp), INTENT(IN) :: qpj(n_vars+2)
    REAL(wp), INTENT(IN) :: B_prime_x
    REAL(wp), INTENT(IN) :: B_prime_y
    REAL(wp), INTENT(IN) :: grav_coeff
    INTEGER, INTENT(IN) :: dir

    REAL(wp), INTENT(OUT) :: flux(n_eqns)

    REAL(wp) :: r_h          !< real-value flow thickness [m]
    REAL(wp) :: r_u          !< real-value x-velocity [m s-1]
    REAL(wp) :: r_v          !< real-value y-velocity [m s-1]
    REAL(wp) :: r_Ri         !< real-value Richardson number
    REAL(wp) :: r_rho_m      !< real-value mixture density [kg m-3]
    REAL(wp) :: r_rho_c      !< real-value carrier phase density [kg m-3]
    REAL(wp) :: r_red_grav   !< real-value reduced gravity [m s-2]
    LOGICAL :: sp_heat_flag
    REAL(wp) :: r_sp_heat_c
    REAL(wp) :: r_sp_heat_mix

    REAL(wp) :: shape_coeff(n_eqns)

    REAL(wp) :: r_w
    REAL(wp) :: mod_vel2 , mod_vel

    REAL(wp) :: shear_stress
    REAL(wp) :: shear_vel    !< shear velocity

    REAL(wp) :: inv_kin_visc
    REAL(wp) :: settling_vel

    REAL(wp) :: Rouse_no(n_solid)

    INTEGER :: i_solid

    sp_heat_flag = .FALSE.

    shape_coeff(1:n_eqns) = 1.0_wp

    pos_thick:IF ( qpj(1) .GT. EPSILON(1.0_wp) ) THEN

       r_h = qpj(1)
       r_u = qpj(idx_u)
       r_v = qpj(idx_v)

       CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav,sp_heat_flag,          &
            r_sp_heat_c,r_sp_heat_mix)

       IF ( vertical_profiles_flag ) CALL eval_flux_coeffs( qpj , B_prime_x ,   &
            B_prime_y , r_rho_c , r_rho_m , shape_coeff )

       IF ( dir .EQ. 1 ) THEN

          ! Mass flux in x-direction: u * ( rhom * h )
          flux(1) = r_u * qcj(1) * shape_coeff(1)

          ! x-momentum flux in x-direction + hydrostatic pressure term
          flux(2) = r_u * qcj(2) * shape_coeff(2) + 0.5_wp * r_rho_m *          &
               grav_coeff * r_red_grav * r_h**2

          ! y-momentum flux in x-direction: u * ( rho * h * v )
          flux(3) = r_u * qcj(3) * shape_coeff(3)

          IF ( energy_flag ) THEN

             ! ENERGY flux in x-direction
             flux(4) = r_u * ( qcj(4) * shape_coeff(4) + 0.5_wp * r_rho_m       &
                  * grav_coeff * r_red_grav * r_h**2 )

          ELSE

             ! Temperature flux in x-direction: u * ( rhom * Cp * h * T )
             flux(4) = r_u * qcj(4) * shape_coeff(4)

          END IF

          ! Mass flux of solid in x-direction: u * ( h * alphas * rhos )
          flux(idx_solidEqn_first:idx_solidEqn_last) = r_u *                    &
               qcj(idx_alfas_first:idx_alfas_last) *                            &
               shape_coeff(idx_alfas_first:idx_alfas_last)

          ! Solid flux can't be larger than total flux
          IF ( ( flux(1) .GT. 0.0_wp ) .AND.                                    &
               ( SUM(flux(idx_solidEqn_first:idx_solidEqn_last)) / flux(1)      &
               .GT.1.0_wp ) ) THEN

             flux(idx_solidEqn_first:idx_solidEqn_last) =                       &
                  flux(idx_solidEqn_first:idx_solidEqn_last) /                  &
                  SUM(flux(idx_solidEqn_first:idx_solidEqn_last)) * flux(1)

          END IF

          ! Mass flux of add.gas in x-direction: u * ( h * alphag * rhog )
          flux(idx_addGasEqn_first:idx_addGasEqn_last) = r_u *                  &
               qcj(idx_addGas_first:idx_addGas_last) *                          &
               shape_coeff(idx_addGas_first:idx_addGas_last)

          IF ( stoch_transport_flag) THEN
             
             ! Flux of stochastic variables
             flux(idx_stochEqn) = r_u * qcj(idx_stoch) * shape_coeff(idx_stoch)

          END IF
             
          IF ( pore_pressure_flag) THEN

             ! Flux of pore pressure
             flux(idx_poreEqn) = r_u * qcj(idx_Pore) * shape_coeff(idx_Pore)

          END IF
             
          ! Mass flux of liquid in x-direction: u * ( h * alphal * rhol )
          IF ( gas_flag .AND. liquid_flag ) flux(n_vars) = r_u * qcj(n_vars) *  &
               shape_coeff(n_vars)

       ELSEIF ( dir .EQ. 2 ) THEN

          ! flux G (derivated wrt y in the equations)
          flux(1) = r_v * qcj(1) * shape_coeff(1)

          flux(2) = r_v * qcj(2) * shape_coeff(2)

          flux(3) = r_v * qcj(3) * shape_coeff(3) + 0.5_wp * r_rho_m *          &
               grav_coeff * r_red_grav * r_h**2

          IF ( energy_flag ) THEN

             ! ENERGY flux in x-direction
             flux(4) = r_v * ( qcj(4) * shape_coeff(4) + 0.5_wp * r_rho_m *     &
                  grav_coeff * r_red_grav * r_h**2 )

          ELSE

             ! Temperature flux in y-direction
             flux(4) = r_v * qcj(4) * shape_coeff(4)

          END IF

          ! Mass flux of solid in y-direction: v * ( h * alphas * rhos )
          flux(idx_solidEqn_first:idx_solidEqn_last) = r_v *                    &
               qcj(idx_alfas_first:idx_alfas_last) *                            &
               shape_coeff(idx_alfas_first:idx_alfas_last)

          ! Solid flux can't be larger than total flux
          IF ( ( flux(1) .GT. 0.0_wp ) .AND.                                    &
               ( SUM(flux(idx_solidEqn_first:idx_solidEqn_last)) / flux(1)      &
               .GT. 1.0_wp ) ) THEN

             flux(idx_solidEqn_first:idx_solidEqn_last) =                       &
                  flux(idx_solidEqn_first:idx_solidEqn_last) /                  &
                  SUM(flux(idx_solidEqn_first:idx_solidEqn_last)) * flux(1)

          END IF

          ! Mass flux of add.gas in x-direction: v * ( h * alphag * rhog )
          flux(idx_addGasEqn_first:idx_addGasEqn_last) = r_v *                   &
               qcj(idx_addGas_first:idx_addGas_last) *                           &
               shape_coeff(idx_addGas_first:idx_addGas_last)

          IF ( stoch_transport_flag) THEN

             ! Flux of stochastic variables
             flux(idx_stochEqn) = r_v * qcj(idx_stoch) * shape_coeff(idx_stoch)

          END IF
             
          IF ( pore_pressure_flag ) THEN
             
             ! Flux of pore pressure
             flux(idx_poreEqn) = r_v * qcj(idx_pore) * shape_coeff(idx_pore)

          END IF
          
          ! Mass flux of liquid in x-direction: u * ( h * alphal * rhol )
          IF ( gas_flag .AND. liquid_flag ) flux(n_vars) = r_v * qcj(n_vars) *  &
               shape_coeff(n_vars)

       END IF

    ELSE

       flux(1:n_eqns) = 0.0_wp

    ENDIF pos_thick

    RETURN

  END SUBROUTINE eval_fluxes


  !******************************************************************************
  !> \brief flux coefficients
  !
  !> This subroutine evaluates the coefficients for the fluxes accouting for the
  !> ratio between fluxes computed from the integral of the profiles and fluxes
  !> computed from the average values of the physical variables.
  !> particles.
  !> \date 2022/03/09
  !> \param[in]     qpj                physical variables
  !> \param[in]     B_prime_x          topography slope in the x-direction
  !> \param[in]     B_prime_y          topography slope in the y-direction
  !> \param[in]     r_rho_c            density of carrier phase
  !> \param[in]     r_rho_m            density of the mixture
  !> \param[out]    shape_coeff          coefficients
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE eval_flux_coeffs(qpj,B_prime_x,B_prime_y,r_rho_c,r_rho_m,shape_coeff)

    USE geometry_2d, ONLY : z_quad, w_quad
    
    USE geometry_2d, ONLY : lambertw,lambertw0,lambertwm1
    USE geometry_2d, ONLY : calcei
    USE geometry_2d, ONLY : gaulegf

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qpj(n_vars+2)
    REAL(wp), INTENT(IN) :: B_prime_x
    REAL(wp), INTENT(IN) :: B_prime_y
    REAL(wp), INTENT(IN) :: r_rho_m !< real-value mixture density [kg m-3]
    REAL(wp), INTENT(IN) :: r_rho_c !< real-value carrier phase density [kg m-3]

    REAL(wp), INTENT(OUT) :: shape_coeff(n_vars)

    REAL(wp) :: Rouse_no(n_solid)

    REAL(wp) :: r_h          !< real-value flow thickness [m]
    REAL(wp) :: r_hu         !< real-value h*x-velocity [m2 s-1]
    REAL(wp) :: r_hv         !< real-value h*y-velocity [m2 s-1]
    REAL(wp) :: r_hw         !< real-value h*z-velocity [m s-1]
    REAL(wp) :: r_u          !< real-value x-velocity [m s-1]
    REAL(wp) :: r_v          !< real-value y-velocity [m s-1]
    REAL(wp) :: r_w          !< real-value z-velocity [m s-1]
    REAL(wp) :: r_alphas(n_solid) !< real-value solid volume fractions
    REAL(wp) :: r_alphag(1:n_add_gas)
    REAL(wp) :: r_T

    REAL(wp) :: mod_hvel

    REAL(wp) :: mod_vel
    REAL(wp) :: mod_vel2

    REAL(wp) :: shear_stress
    REAL(wp) :: shear_vel    !< shear velocity

    REAL(wp) :: inv_kin_visc
    REAL(wp) :: settling_vel

    REAL(wp) :: rhom_mod_vel2

    REAL(wp) :: h_rel

    REAL(wp) :: a , b, c, d

    REAL(wp) :: h0_rel
    REAL(wp) :: h0_rel_1
    REAL(wp) :: h0_rel_2

    REAL(wp) :: normalizing_coeff_u

    REAL(wp) :: h0

    REAL(wp) :: u_rel0

    INTEGER :: i_solid

    REAL(wp) :: z(n_quad)
    REAL(wp) :: w(n_quad)

    REAL(wp) :: epsilon_s
    REAL(wp) :: a_coeff

    REAL(wp) :: log_term_h0

    REAL(wp) :: log_term_z(n_quad)
    REAL(wp) :: w_log_term_z(n_quad)

    REAL(wp) :: alphas_exp_avg
    REAL(wp) :: u_log_avg
    
    REAL(wp) :: normalizing_coeff_alpha
    REAL(wp) :: alphas_rel0
    REAL(wp) :: rho_u_alphas(n_solid)
    REAL(wp) :: y
    REAL(wp) :: int_quad

    REAL(wp) :: int_def1 , int_def2
    
    REAL(wp) :: ei
    REAL(wp) :: uRho_avg
    REAL(wp) :: int_hvel_vel
    REAL(wp) :: int_rhom_hvel_vel
    REAL(wp) :: rhom_vel_vel
    REAL(wp) :: dyn_visc_c
    REAL(wp) :: r_inv_rho_c

    shape_coeff(1:n_eqns) = 1.0_wp

    r_h = qpj(1)
    r_hu = qpj(2)
    r_hv = qpj(3)
    r_T = qpj(4)
    r_u = qpj(idx_u)
    r_v = qpj(idx_v)

    IF ( alpha_flag ) THEN

       r_alphas(1:n_solid) = qpj(idx_alfas_first:idx_alfas_last)
       r_alphag(1:n_add_gas) = qpj(idx_addGas_first:idx_addGas_last)

    ELSE

       r_alphas(1:n_solid) = qpj(idx_alfas_first:idx_alfas_last) / qpj(1)
       r_alphag(1:n_add_gas) = qpj(idx_addGas_first:idx_addGas_last) / qpj(1)

    END IF

    IF ( slope_correction_flag ) THEN

       r_hw = r_hu * B_prime_x + r_hv * B_prime_y
       r_w = r_u * B_prime_x + r_v * B_prime_y

    ELSE

       r_hw = 0.0_wp
       r_w = 0.0_wp

    END IF

    mod_hvel = SQRT(r_hu**2 + r_hv**2 + r_hw**2)

    mod_vel2 = r_u**2 + r_v**2 + r_w**2
    mod_vel = SQRT( mod_vel2 )

    shear_vel = SQRT( friction_factor ) * mod_vel

    IF ( gas_flag .AND. sutherland_flag ) THEN
          
       dyn_visc_c = muRef_Suth * ( r_T / Tref_Suth )**1.5_wp *               &
            ( Tref_Suth + S_mu ) / ( r_T + S_mu )

       r_inv_rho_c = sp_gas_const_a * r_T * inv_pres       
       kin_visc_c = dyn_visc_c * r_inv_rho_c
       
    END IF
                  
    ! Viscosity read from input file [m2 s-1]
    inv_kin_visc = 1.0_wp / kin_visc_c

    DO i_solid=1,n_solid

       settling_vel = settling_velocity( diam_s(i_solid) , rho_s(i_solid) ,     &
            r_rho_c , inv_kin_visc )

       IF ( shear_vel .GT. 0.0_wp ) THEN

          Rouse_no(i_solid) = settling_vel / ( vonK * shear_vel )

       ELSE

          Rouse_no(i_solid) = 0.0_wp

       END IF

    END DO

    ! The profile parameters depend on h/k_s, not on the absolute value of h.
    h_rel = r_h / k_s

    IF ( h_rel .GT. H_crit_rel ) THEN

       ! we search for h0_rel such that the average integral between 0 and
       ! h_rel is equal to 1
       ! For h_rel > H_crit_rel this integral is the sum of two pieces:
       ! integral between 0 and h0_rel of the log profile
       ! integral between h0_rel and h_rel of the costant profile

       a = h_rel * vonK / SQRT(friction_factor)
       b = 1.0_wp / 30.0_wp + h_rel
       c = 30.0_wp

       ! solve b*log(c*z+1)-z=a for z
       d = a/b-1.0_wp/(b*c)

       h0_rel_1 = -b*lambertw0( -EXP(d)/(b*c) ) - 1.0_wp / c
       h0_rel_2 = -b*lambertwm1( -EXP(d)/(b*c) ) - 1.0_wp / c
       h0_rel = MIN( h0_rel_1 , h0_rel_2)

    ELSE

       ! when h_rel <= H_crit_rel we have only the log profile and we have to
       ! rescale it in order to have the integral between o and h_rel equal to
       ! 1
       h0_rel = h_rel

    END IF

    h0 = h0_rel*k_s

    b = 30.0_wp / k_s

    log_term_h0 =  u_log_profile(b,h0)**2
    
    z = 0.5 * h0 * (z_quad+1.0_wp)
    w = 0.5 * h0 * w_quad
    
    log_term_z =  u_log_profile(b,z)**2

    w_log_term_z = w * log_term_z

    u_log_avg = ( SUM( w * u_log_profile(b,z) ) + u_log_profile(b,h0)*(r_h-h0) ) &
         / r_h
    
    normalizing_coeff_u = 1.0_wp / u_log_avg

    ! relative velocity at top of boundary layer (z=h0)
    u_rel0 = normalizing_coeff_u * u_log_profile(b,h0)

    ! constant terms for the concentration profiles within the loop
    epsilon_s = Sc * shear_vel * vonK * ( ( h0/6.0_wp ) + ( k_s / 60.0_wp ) )
    a_coeff = - vonK * shear_vel / epsilon_s

    ! initial partial sum to zero
    int_rhom_hvel_vel = 0.0_wp

    DO i_solid = 1,n_solid

       a = a_coeff * Rouse_no(i_solid)

       alphas_exp_avg = ( SUM( w * alphas_exp_profile(a,z) ) +                  &
            alphas_exp_profile(a,h0) * ( r_h - h0 ) ) / r_h
       
       normalizing_coeff_alpha = 1.0_wp / alphas_exp_avg

       ! integral of profiles
       int_quad = SUM( w * ( alphas_exp_profile(a,z) * u_log_profile(b,z) ) )

       ! integral of alfa_i*u between in the boundary layer
       int_def1 = ( mod_vel * normalizing_coeff_u ) *                           &
            ( r_alphas(i_solid) * normalizing_coeff_alpha ) * int_quad
       
       ! relative concentration alphas_rel at h0
       alphas_rel0 = normalizing_coeff_alpha * alphas_exp_profile(a,h0)
       
       ! integral of alfa_i*u in the free-stream layer
       int_def2 =  ( r_h - h0 ) * ( mod_vel * u_rel0 ) *                        &
            ( r_alphas(i_solid) * alphas_rel0 )
       
       ! we sum the two integrals, multiply by the density of solid and
       ! divide by thickness h to find the depth-average value
       rho_u_alphas(i_solid) = rho_s(i_solid) * ( int_def1 + int_def2 ) / r_h
       
       ! integral of alphas_rel(z)*u_rel(z)^2 in the log region (0<=z<=h0)
       int_def1 = SUM( alphas_exp_profile(a,z) * w_log_term_z )

       ! integral in the free-stream region
       int_def2 = ( r_h - h0 ) *  alphas_exp_profile(a,h0) * log_term_h0
       
       int_rhom_hvel_vel = int_rhom_hvel_vel + ( rho_s(i_solid) - r_rho_c ) *   &
            ( normalizing_coeff_alpha * r_alphas(i_solid) ) *                   &
            ( mod_hvel * mod_vel * normalizing_coeff_u**2 ) *                   &
            ( int_def1 + int_def2 )

    END DO

    ! we add the contribution of the gas phase to the mixture depth-averaged 
    ! momentum
    uRho_avg = ( mod_vel * r_rho_c + SUM((rho_s - r_rho_c) / rho_s *            &
         rho_u_alphas) )

    ! integral of square of relative velocity between 0 and h0
    int_def1 = SUM( w *  u_log_profile(b,z)**2 )
    ! int_def1 = SUM( w * (LOG( b*z + 1.0_wp ))**2 )

    ! integral of square of relative velocity between h0 and h
    int_def2 = ( r_h - h0 ) * log_term_h0 
    
    ! integral of (h*velocity)*velocity profile between 0 and h
    int_hvel_vel = mod_hvel * mod_vel *                                         &
         normalizing_coeff_u**2 * ( int_def1 + int_def2 )

    ! depth-averaged momentum flux
    rhom_vel_vel = ( int_rhom_hvel_vel + r_rho_c * int_hvel_vel ) / r_h**2

    shape_coeff(1) = uRho_avg / ( mod_vel * r_rho_m )

    shape_coeff(2:3) = rhom_vel_vel / ( uRho_avg * mod_vel )

    ! we assume a vertically-constant temperature profile
    shape_coeff(4) = shape_coeff(1)

    shape_coeff(idx_alfas_first:idx_alfas_last) = rho_u_alphas /                &
         ( rho_s * r_alphas * mod_vel )

    RETURN


  END SUBROUTINE eval_flux_coeffs


  !******************************************************************************
  !> \brief deposition coefficients
  !
  !> This subroutine evaluates the coefficients for the deposition terms 
  !> accounting for the ratio between bottom and average volume fraction of solid
  !> particles.
  !> \date 2022/03/09
  !> \param[in]     qpj                physical variables
  !> \param[in]     mod_vel            velocity magnitude
  !> \param[in]     r_rho_c            density of carrier phase
  !> \param[in]     r_rho_m            density of the mixture
  !> \param[out]    dep_coeff          coefficients
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE eval_dep_coeffs( qpj , mod_vel , r_rho_c , r_rho_m , dep_coeff )

    USE geometry_2d, ONLY : z_quad , w_quad
    
    USE geometry_2d, ONLY : lambertw,lambertw0,lambertwm1
    USE geometry_2d, ONLY : calcei
    USE geometry_2d, ONLY : gaulegf
    

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qpj(n_vars+2)
    REAL(wp), INTENT(IN) :: mod_vel
    REAL(wp), INTENT(IN) :: r_rho_m !< real-value mixture density [kg m-3]
    REAL(wp), INTENT(IN) :: r_rho_c !< real-value carrier phase density [kg m-3]

    REAL(wp), INTENT(OUT) :: dep_coeff(n_solid)

    REAL(wp) :: Rouse_no(n_solid)

    REAL(wp) :: r_h          !< real-value flow thickness [m]
    REAL(wp) :: r_T

    REAL(wp) :: shear_vel    !< shear velocity

    REAL(wp) :: inv_kin_visc
    REAL(wp) :: settling_vel

    REAL(wp) :: h_rel

    REAL(wp) :: a , b, c, d

    REAL(wp) :: h0_rel
    REAL(wp) :: h0_rel_1
    REAL(wp) :: h0_rel_2

    REAL(wp) :: h0

    INTEGER :: i_solid

    REAL(wp) :: epsilon_s
    REAL(wp) :: a_coeff

    REAL(wp) :: alphas_exp_avg
    REAL(wp) :: y

    REAL(wp) :: z(n_quad)
    REAL(wp) :: w(n_quad)
    REAL(wp) :: dyn_visc_c

    dep_coeff(1:n_eqns) = 1.0_wp

    r_h = qpj(1)
    r_T = qpj(4)

    shear_vel = SQRT( friction_factor ) * mod_vel

    IF ( gas_flag .AND. sutherland_flag ) THEN
          
       dyn_visc_c = muRef_Suth * ( r_T / Tref_Suth )**1.5_wp *               &
            ( Tref_Suth + S_mu ) / ( r_T + S_mu )

       kin_visc_c = dyn_visc_c / r_rho_c
       
    END IF
    
    ! Viscosity read from input file [m2 s-1]
    inv_kin_visc = 1.0_wp / kin_visc_c

    DO i_solid=1,n_solid

       settling_vel = settling_velocity( diam_s(i_solid) , rho_s(i_solid) ,     &
            r_rho_c , inv_kin_visc )

       IF ( shear_vel .GT. 0.0_wp ) THEN

          Rouse_no(i_solid) = settling_vel / ( vonK * shear_vel )

       ELSE

          Rouse_no(i_solid) = 0.0_wp

       END IF

    END DO

    ! The profile parameters depend on h/k_s, not on the absolute value of h
    h_rel = r_h / k_s

    IF ( h_rel .GT. H_crit_rel ) THEN

       ! we search for h0_rel such that the average integral between 0 and
       ! h_rel is equal to 1
       ! For h_rel > H_crit_rel this integral is the sum of two pieces:
       ! integral between 0 and h0_rel of the log profile
       ! integral between h0_rel and h_rel of the costant profile

       a = h_rel * vonK / SQRT(friction_factor)
       b = 1.0_wp / 30.0_wp + h_rel
       c = 30.0_wp

       ! solve b*log(c*z+1)-z=a for z
       d = a/b-1.0_wp/(b*c)

       h0_rel_1 = -b*lambertw0( -EXP(d)/(b*c) ) - 1.0_wp / c
       h0_rel_2 = -b*lambertwm1( -EXP(d)/(b*c) ) - 1.0_wp / c
       h0_rel = MIN( h0_rel_1 , h0_rel_2)


       ! h0_rel = -b*lambertw( -EXP(d)/(b*c) ) - 1.0_wp / c

    ELSE

       ! when h_rel <= H_crit_rel we have only the log profile and we have to
       ! rescale it in order to have the integral between o and h_rel equal to 1

       h0_rel = h_rel

    END IF

    h0 = h0_rel*k_s

    epsilon_s = Sc * shear_vel * vonK * ( ( h0/6.0_wp ) + ( k_s / 60.0_wp ) )
    a_coeff = - vonK * shear_vel / epsilon_s

    z = 0.5_wp * h0 * ( z_quad + 1.0_wp )
    w = 0.5_wp * h0 * w_quad
    
    DO i_solid = 1,n_solid

       a = a_coeff * Rouse_no(i_solid)

       alphas_exp_avg = ( SUM( w * alphas_exp_profile(a,z) ) +                  &
            alphas_exp_profile(a,h0)*(r_h-h0) ) / r_h

       dep_coeff(i_solid) = 1.0_wp / alphas_exp_avg

    END DO

    RETURN


  END SUBROUTINE eval_dep_coeffs


  !******************************************************************************
  !> \brief Explicit source term
  !
  !> This subroutine evaluates the non-hyperbolic terms to be treated explicitely
  !> in the DIRK numerical scheme (e.g. gravity,source of mass). The sign of the
  !> terms is taken with the terms on the right-hand side of the equations.
  !> \date 2019/12/13
  !> \param[in]     B_primej_x         local x-slope
  !> \param[in]     B_primej_y         local y_slope
  !> \param[in]     B_secondj_xx       local 2nd derivative in x-direction
  !> \param[in]     B_secondj_xy       local 2nd derivative in xy-direction
  !> \param[in]     B_secondj_yy       local 2nd derivative in y-direction
  !> \param[in]     grav_coeff         correction factor for topography slope
  !> \param[in]     d_grav_coeff_dx    x-derivative of grav_coeff
  !> \param[in]     d_grav_coeff_dy    y-derivative of grav_coeff
  !> \param[in]     source_xy          local source of mass
  !> \param[in]     qpj                physical variables
  !> \param[in]     time               simlation time (needed for source)
  !> \param[in]     cell_fract_jk      fraction of cell contributing to source
  !> \param[out]    expl_term          explicit term
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE eval_expl_terms( Bprimej_x, Bprimej_y, Bsecondj_xx , Bsecondj_xy , &
       Bsecondj_yy, grav_coeff, d_grav_coeff_dx , d_grav_coeff_dy , source_xy , &
       qpj, expl_term, time, cell_fract_jk )

    USE parameters_2d, ONLY : vel_source , T_source , alphas_source ,           &
         alphal_source , time_param , bottom_radial_source_flag , alphag_source,&
         pore_pressure_flag , pore_pres_fract


    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: Bprimej_x
    REAL(wp), INTENT(IN) :: Bprimej_y
    REAL(wp), INTENT(IN) :: Bsecondj_xx
    REAL(wp), INTENT(IN) :: Bsecondj_xy
    REAL(wp), INTENT(IN) :: Bsecondj_yy
    REAL(wp), INTENT(IN) :: grav_coeff
    REAL(wp), INTENT(IN) :: d_grav_coeff_dx
    REAL(wp), INTENT(IN) :: d_grav_coeff_dy

    REAL(wp), INTENT(IN) :: source_xy

    REAL(wp), INTENT(IN) :: qpj(n_vars+2)      !< local physical variables 
    REAL(wp), INTENT(OUT) :: expl_term(n_eqns) !< local explicit forces 

    REAL(wp), INTENT(IN) :: time
    REAL(wp), INTENT(IN) :: cell_fract_jk

    REAL(wp) :: r_h          !< real-value flow thickness
    REAL(wp) :: r_u          !< real-value x-velocity
    REAL(wp) :: r_v          !< real-value y-velocity
    REAL(wp) :: r_Ri         !< real-value Richardson number
    REAL(wp) :: r_rho_m      !< real-value mixture density [kg/m3]
    REAL(wp) :: r_rho_c      !< real-value carrier phase density [kg/m3]
    REAL(wp) :: r_rho_g(n_add_gas) !< real-value add.gas densities [kg/m3]
    REAL(wp) :: r_red_grav   !< real-value reduced gravity

    REAL(wp) :: r_sp_heat_mix !< real_value mixture specific heat
    REAL(wp) :: r_sp_heat_c

    REAL(wp) :: t_rem
    REAL(wp) :: t_coeff
    REAL(wp) :: h_dot

    REAL(wp) :: r_tilde_grav
    REAL(wp) :: centr_force_term

    REAL(wp) :: qp_source(n_vars+2)

    REAL(wp) :: exc_pore_pres

    REAL(wp) :: q1
    
    LOGICAL :: sp_heat_flag

    sp_heat_flag = .TRUE.

    expl_term(1:n_eqns) = 0.0_wp

    IF ( ( qpj(1) .LE. EPSILON(1.0_wp) ) .AND. ( cell_fract_jk .EQ. 0.0_wp ) )  &
         RETURN

    r_h = qpj(1)
    r_u = qpj(idx_u)
    r_v = qpj(idx_v)

    CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav,sp_heat_flag,r_sp_heat_c, &
         r_sp_heat_mix)

    q1 = r_h * r_rho_m
        
    IF ( curvature_term_flag ) THEN

       centr_force_term = Bsecondj_xx * r_u**2 + 2.0_wp*Bsecondj_xy*r_u * r_v + &
            Bsecondj_yy * r_v**2

    ELSE

       centr_force_term = 0.0_wp

    END IF

    r_tilde_grav = r_red_grav + centr_force_term

    ! units of dqc(2)/dt [kg m-1 s-2]
    expl_term(2) = - grav_coeff * r_rho_m * r_tilde_grav * r_h * Bprimej_x      &
         + 0.5_wp * r_rho_m * r_red_grav * r_h**2 * d_grav_coeff_dx 

    ! units of dqc(3)/dt [kg m-1 s-2]
    expl_term(3) = - grav_coeff * r_rho_m * r_tilde_grav * r_h * Bprimej_y      &
         + 0.5_wp * r_rho_m * r_red_grav * r_h**2 * d_grav_coeff_dy

    IF ( energy_flag ) THEN

       expl_term(4) = expl_term(2) * r_u + expl_term(3) * r_v
       
    ELSE

       expl_term(4) = 0.0_wp

    END IF

    ! ----------- ADDITIONAL EXPLICIT TERMS FOR BOTTOM RADIAL SOURCE ------------ 

    IF ( ( .NOT.bottom_radial_source_flag ) .OR.                                &
         ( cell_fract_jk .EQ. 0.0_wp ) ) RETURN

    t_rem = MOD( time + time_param(4) , time_param(1) )

    IF ( time_param(3) .EQ. 0.0_wp ) THEN

       IF ( t_rem .LE. time_param(2) ) THEN

          t_coeff = 1.0_wp

       ELSE

          t_coeff = 0.0_wp

       END IF

    ELSE

       IF ( t_rem .LE. time_param(3) ) THEN

          t_coeff = ( t_rem / time_param(3) ) 

       ELSEIF ( t_rem .LE. time_param(2) - time_param(3) ) THEN

          t_coeff = 1.0_wp

       ELSEIF ( t_rem .LE. time_param(2) ) THEN

          t_coeff = 1.0_wp - ( t_rem - time_param(2) + time_param(3) ) /        &
               time_param(3)

       ELSE

          t_coeff = 0.0_wp

       END IF

    END IF

    h_dot = cell_fract_jk * vel_source

    qp_source(1) = 1.0_wp
    qp_source(2) = 0.0_wp
    qp_source(3) = 0.0_wp
    qp_source(4) = t_source

    IF ( alpha_flag ) THEN

       qp_source(idx_alfas_first:idx_alfas_last) = alphas_source(1:n_solid)
       qp_source(idx_addGas_first:idx_addGas_last) = alphag_source(1:n_add_gas)
       IF ( gas_flag .AND. liquid_flag ) qp_source(n_vars) = alphal_source

    ELSE

       qp_source(idx_alfas_first:idx_alfas_last) = alphas_source(1:n_solid)     &
            * qp_source(1)
       qp_source(idx_addGas_first:idx_addGas_last) = alphag_source(1:n_add_gas) &
            * qp_source(1)
       IF ( gas_flag .AND. liquid_flag ) qp_source(n_vars) = alphal_source      &
            * qp_source(1)

    END IF

    ! Source term transport stochastic equation
    IF ( stoch_transport_flag) qp_source(idx_stoch) = 0.0_wp

    IF ( pore_pressure_flag ) qp_source(idx_poreEqn) = 0.0_wp

    qp_source(idx_u) = 0.0_wp
    qp_source(idx_v) = 0.0_wp


    CALL mixt_var(qp_source,r_Ri,r_rho_m,r_rho_c,r_red_grav,sp_heat_flag,       &
         r_sp_heat_c,r_sp_heat_mix)

    expl_term(1) = expl_term(1) + t_coeff * h_dot * r_rho_m
    expl_term(2) = expl_term(2) + 0.0_wp
    expl_term(3) = expl_term(3) + 0.0_wp

    IF ( energy_flag ) THEN

       expl_term(4) = expl_term(4) + t_coeff * h_dot * r_rho_m * r_sp_heat_mix  &
            * t_source

    ELSE

       expl_term(4) = expl_term(4) + t_coeff * h_dot * r_rho_m * r_sp_heat_mix  &
            * t_source

    END IF

    ! source terms for the solid equations
    expl_term(idx_alfas_first:idx_alfas_last) =                                 &
         expl_term(idx_alfas_first:idx_alfas_last) + t_coeff                   &
         * h_dot * alphas_source(1:n_solid) * rho_s(1:n_solid)

    r_rho_g(1:n_add_gas) = pres / ( sp_gas_const_g(1:n_add_gas) * t_source )

    ! source terms for the additional gas equations
    expl_term(idx_addGasEqn_first:idx_addGasEqn_last) =                         &
         expl_term(idx_addGasEqn_first:idx_addGasEqn_last) + t_coeff            &
         * h_dot * alphag_source(1:n_add_gas) * r_rho_g(1:n_add_gas)

    IF ( gas_flag .AND. liquid_flag ) THEN

       ! source term for the liquid phase
       expl_term(n_vars) = expl_term(n_vars) + t_coeff * h_dot * alphal_source  &
            * rho_l

    END IF

    IF ( pore_pressure_flag ) THEN

       exc_pore_pres = qpj(idx_pore)
       
       ! we multiply the pore pressure inlet rate by the rate for q1
       ! the units of this source term are: kg^2 m^-3 s^-3
       expl_term(idx_poreEqn) = expl_term(idx_poreEqn) +                       &
            t_coeff * ( q1 * pore_pres_fract * h_dot * r_rho_m * r_red_grav +  &
            exc_pore_pres * h_dot * r_rho_m )

    END IF
    
    RETURN

  END SUBROUTINE eval_expl_terms


  !******************************************************************************
  !> \brief Analytic integration of friction terms
  !
  !> This subroutine integrate analytically the friction term for turbulent
  !> friction only.
  !> \date 2021/12/10
  !> \param[inout]     r_qj            real conservative variables 
  !> \param[in]        dt              time step 
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE integrate_friction_term( r_qj , dt )

    IMPLICIT NONE

    REAL(wp), INTENT(INOUT) :: r_qj(n_vars)
    REAL(wp), INTENT(IN) :: dt 


    REAL(wp) :: r_qp(n_vars+2)
    REAL(wp) :: p_dyn
    REAL(wp) :: r_h
    REAL(wp) :: r_u
    REAL(wp) :: r_v

    REAL(wp) :: mod_vel_hor
    REAL(wp) :: mod_vel

    IF ( r_qj(1) .EQ. 0.0_wp ) RETURN

    CALL qc_to_qp(r_qj , r_qp , p_dyn )
    
    r_h = r_qp(1)

    r_u = r_qp(idx_u)
    r_v = r_qp(idx_v)

    mod_vel_hor = SQRT( r_u**2 + r_v**2 )

    mod_vel = 1.0_wp / ( 1.0_wp  / mod_vel_hor + friction_factor / r_h * dt )

    r_u = r_u * ( mod_vel / mod_vel_hor )
    r_v = r_v * ( mod_vel / mod_vel_hor )

    r_qp(2) = r_h*r_u
    r_qp(3) = r_h*r_v

    r_qp(idx_u) = r_u
    r_qp(idx_v) = r_v

    CALL qp_to_qc(r_qp,r_qj)

    RETURN

  END SUBROUTINE integrate_friction_term


  !******************************************************************************
  !> \brief Implicit source terms
  !
  !> This subroutine evaluates the source terms  of the system of equations,
  !> both for real or complex inputs, that are treated implicitely in the DIRK
  !> numerical scheme.
  !> \date 01/06/2012
  !> \param[in]     Bprimej_x       topography slope in x-direction
  !> \param[in]     Bprimej_y       topography slope in y-direction
  !> \param[in]     c_qj            complex conservative variables 
  !> \param[in]     r_qj            real conservative variables 
  !> \param[out]    c_nh_term_impl  complex non-hyperbolic terms     
  !> \param[out]    r_nh_term_impl  real non-hyperbolic terms
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE eval_implicit_terms( Bprimej_x, Bprimej_y, Zij, fric_val, c_qj,    &
       c_nh_term_impl, r_qj , r_nh_term_impl )

    USE COMPLEXIFY

    USE geometry_2d, ONLY : pi_g

    USE parameters_2d, ONLY : pore_pressure_flag
    USE parameters_2d, ONLY : four_thirds , neg_four_thirds

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: Bprimej_x
    REAL(wp), INTENT(IN) :: Bprimej_y
    REAL(wp), INTENT(IN):: Zij
    REAL(wp), INTENT(OUT) :: fric_val
    COMPLEX(wp), INTENT(IN), OPTIONAL :: c_qj(n_vars)
    COMPLEX(wp), INTENT(OUT), OPTIONAL :: c_nh_term_impl(n_eqns)
    REAL(wp), INTENT(IN), OPTIONAL :: r_qj(n_vars)
    REAL(wp), INTENT(OUT), OPTIONAL :: r_nh_term_impl(n_eqns)

    COMPLEX(wp) :: h                       !< height [m]
    COMPLEX(wp) :: inv_h                   !< 1/height [m-1]
    COMPLEX(wp) :: u                       !< velocity (x direction) [m/s]
    COMPLEX(wp) :: v                       !< velocity (y direction) [m/s]
    COMPLEX(wp) :: w                       !< velocity (z direction) [m/s]
    COMPLEX(wp) :: T                       !< temperature [K]
    COMPLEX(wp) :: rho_m                   !< mixture density [kg/m3]
    COMPLEX(wp) :: alphas(n_solid)         !< sediment volume fractions
    COMPLEX(wp) :: alphag(n_solid)         !< add. gas volume fractions
    COMPLEX(wp) :: inv_rho_m               !< 1/mixture density [kg-1 m3]

    COMPLEX(wp) :: qj(n_vars)
    COMPLEX(wp) :: nh_term(n_eqns)
    COMPLEX(wp) :: source_term(n_eqns)

    COMPLEX(wp) :: mod_vel
    COMPLEX(wp) :: mod_vel_hor
    COMPLEX(wp) :: mod_vel2
    COMPLEX(wp) :: gamma
    COMPLEX(wp) :: rho_g(n_add_gas)
    REAL(wp) :: h_threshold
    COMPLEX(wp) :: rho_c

    INTEGER :: i

    !--- Lahars rheology model variables

    !> Temperature in C
    COMPLEX(wp) :: Tc

    COMPLEX(wp) :: expA , expB

    !> 1st param for fluid viscosity empirical relationship (O'Brian et al, 1993)
    COMPLEX(wp) :: alpha1    ! (units: kg m-1 s-1 )

    !> Fluid dynamic viscosity (units: kg m-1 s-1 )
    COMPLEX(wp) :: fluid_visc

    !> Total friction slope (dimensionless): s_f = s_v+s_td+s_y
    COMPLEX(wp) :: s_f

    !> Viscous slope component of total Friction (dimensionless)
    COMPLEX(wp) :: s_v

    !> Turbulent dispersive slope component of total friction (dimensionless)
    COMPLEX(wp) :: s_td

    COMPLEX(wp) :: temp_term

    COMPLEX(wp) :: c_tau

    COMPLEX(wp) :: Zs

    COMPLEX(wp) :: exc_pore_pres

    COMPLEX(wp) :: D_coeff
    REAL(wp) :: gamma_gas
    COMPLEX(wp) :: gas_compressibility

    COMPLEX(wp) :: rho_gas

    COMPLEX(wp) :: porosity
    COMPLEX(wp) :: f_inhibit
    COMPLEX(wp) :: dyn_visc_c
    COMPLEX(wp) :: red_grav

    COMPLEX(wp) :: turb_stress
    
    IF ( present(c_qj) .AND. present(c_nh_term_impl) ) THEN

       qj = c_qj

    ELSEIF ( present(r_qj) .AND. present(r_nh_term_impl) ) THEN

       DO i = 1,n_vars

          qj(i) = CMPLX( r_qj(i),0.0_wp,wp )

       END DO

    ELSE

       WRITE(*,*) 'Constitutive, eval_fluxes: problem with arguments'
       STOP

    END IF

    ! initialize the source terms
    source_term(1:n_eqns) = CMPLX(0.0_wp,0.0_wp,wp)

    IF (rheology_flag) THEN

       CALL c_phys_var(qj,h,u,v,T,rho_m,alphas,alphag,inv_rho_m,Zs,             &
            exc_pore_pres)

       red_grav = ( rho_m - rho_a_amb ) / rho_m * grav
       
       IF ( slope_correction_flag ) THEN

          w = u * Bprimej_x + v * Bprimej_y

       ELSE

          w = CMPLX( 0.0_wp , 0.0_wp , wp )

       END IF

       mod_vel2 = u**2 + v**2 + w**2
       mod_vel = SQRT( mod_vel2 )
       mod_vel_hor = SQRT( u**2 + v**2 )

       IF ( rheology_model .EQ. 1 ) THEN
          ! Voellmy Salm rheology

          IF ( REAL(mod_vel) .NE. 0.0_wp ) THEN 

             ! Modify xi if using a stochastic model
             IF (stochastic_flag) THEN 
               xi_temp = xi + Zij
               ! Limit the boundaies of stochastic friction
               IF (xi_temp .LT. 1._wp) xi_temp = 1._wp
               IF (xi_temp .GT. 1e5_wp) xi_temp = 1e5_wp   
             ELSE
               xi_temp = xi
             END IF  
             
             ! Store stoch value of xi
             fric_val = xi_temp

             turb_stress = rho_m * red_grav / xi_temp * mod_vel2

             source_term(2) = source_term(2) - turb_stress * ( u / mod_vel )
             source_term(3) = source_term(3) - turb_stress * ( v / mod_vel )

          ENDIF

       ELSEIF ( rheology_model .EQ. 2 ) THEN

          ! Plastic rheology
          IF ( REAL(mod_vel) .NE. 0.0_wp ) THEN 

             source_term(2) = source_term(2) - rho_m * tau * ( u / mod_vel )

             source_term(3) = source_term(3) - rho_m * tau * ( v / mod_vel )

          ENDIF

       ELSEIF ( rheology_model .EQ. 3 ) THEN

          h_threshold = 1.0E-10_wp

          T_env = 300.0_wp

          ! Temperature dependent rheology
          IF ( REAL(h) .GT. h_threshold ) THEN

             ! Equation 6 from Costa & Macedonio, 2005
             gamma = 3.0_wp * nu_ref / h * EXP( - visc_par * ( T - T_ref ) )

          ELSE

             ! Equation 6 from Costa & Macedonio, 2005
             gamma = 3.0_wp * nu_ref / h_threshold * EXP( - visc_par            &
                  * ( T - T_ref ) )

          END IF

          IF ( REAL(mod_vel) .NE. 0.0_wp ) THEN 

             ! Last R.H.S. term in equation 2 from Costa & Macedonio, 2005
             source_term(2) = source_term(2) - rho_m * gamma * u

             ! Last R.H.S. term in equation 3 from Costa & Macedonio, 2005
             source_term(3) = source_term(3) - rho_m * gamma * v

          ENDIF

          source_term(4) = source_term(4) - radiative_term_coeff * ( T**4 -     &
               T_env**4 ) - convective_term_coeff * ( T - T_env )


       ELSEIF ( rheology_model .EQ. 4 ) THEN

          ! Lahars rheology (O'Brien 1993, FLO2D)

          ! alpha1 here has units: kg m-1 s-1
          ! in Table 2 from O'Brien 1988, the values reported have different
          ! units ( poises). 1poises = 0.1 kg m-1 s-1

          IF ( wp .EQ. sp ) THEN

             h_threshold = 1.0E-10_wp

          ELSE

             h_threshold = 1.0E-20_wp

          END IF

          ! convert from Kelvin to Celsius
          Tc = T - 273.15_wp

          ! the dependance of viscosity on temperature is modeled with the
          ! equation presented at:
          ! https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781118131473.app3
          !
          ! In addition, we use a reference value provided in input at a 
          ! reference temperature. This value is used to scale the equation
          IF ( REAL(Tc) .LT. 20.0_wp ) THEN

             expA = 1301.0_wp / ( 998.333_wp + 8.1855_wp * ( Tc - 20.0_wp )     &
                  + 0.00585_wp * ( Tc - 20.0_wp )**2 ) - 1.30223_wp

             alpha1 = alpha1_coeff * 1.0E-3_wp * 10.0_wp**expA

          ELSE

             expB = ( 1.3272_wp * ( 20.0_wp - Tc ) - 0.001053_wp *              &
                  ( Tc - 20.0_wp )**2 ) / ( Tc + 105.0_wp )

             alpha1 = alpha1_coeff * 1.002E-3_wp * 10.0_wp**expB 

          END IF

          ! Fluid dynamic viscosity [kg m-1 s-1] 
          fluid_visc = alpha1 * EXP( beta1 * SUM(alphas) )

          IF ( REAL(h) .GT. h_threshold ) THEN

             inv_h = 1.0_wp / h

             ! Viscous slope component (dimensionless)
             s_v = Kappa * fluid_visc * mod_vel * 0.125_wp * inv_rho_m *        &
                  inv_grav * inv_h**2

             ! Turbulent dispersive component (dimensionless)
             s_td = n_td2 * mod_vel2 * inv_h**four_thirds

          ELSE

             ! Viscous slope component (dimensionless)
             s_v = Kappa * fluid_visc * mod_vel / ( 8.0_wp * rho_m * grav *     &
                  h_threshold**2 )

             ! Turbulent dispersive components (dimensionless)
             s_td = n_td2 * mod_vel2 * h_threshold**neg_four_thirds

          END IF

          ! Total implicit friction slope (dimensionless)
          s_f = s_v + s_td

          IF ( REAL(mod_vel) .GT. 0.0_wp ) THEN

             temp_term = grav * rho_m * h * s_f / mod_vel_hor

             ! same units of dqc(2)/dt: kg m-1 s-2
             source_term(2) = source_term(2) - u * temp_term

             ! same units of dqc(3)/dt: kg m-1 s-2
             source_term(3) = source_term(3) - v * temp_term

          END IF

       ELSEIF ( rheology_model .EQ. 5 ) THEN

          c_tau = 1.0E-3_wp / ( 1.0_wp + 10.0_wp * h ) * mod_vel

          IF ( REAL(mod_vel) .NE. 0.0_wp ) THEN

             source_term(2) = source_term(2) - rho_m * c_tau * ( u/mod_vel_hor )
             source_term(3) = source_term(3) - rho_m * c_tau * ( v/mod_vel_hor )

          END IF


       ELSEIF ( rheology_model .EQ. 6 ) THEN

          IF ( REAL(mod_vel) .NE. 0.0_wp ) THEN 

             source_term(2) = source_term(2) - rho_m * ( u / mod_vel_hor ) *    &
                  friction_factor * mod_vel2

             source_term(3) = source_term(3) - rho_m * ( v / mod_vel_hor ) *    &
                  friction_factor * mod_vel2

          ENDIF

       ! Coulomb function rheology: mu(Fr)
       ELSEIF ( rheology_model .EQ. 9 ) THEN

       ENDIF

    ENDIF

    IF ( pore_pressure_flag ) THEN

       h_threshold = 1.0E-10_wp

       gamma_gas = sp_heat_a / ( sp_heat_a - sp_gas_const_a )
       gas_compressibility = 1.0_wp / ( gamma_gas * pres )
       rho_gas = pres / ( sp_gas_const_a * T )

       porosity = 1.0_wp - SUM(alphas)

       IF ( gas_flag .AND. sutherland_flag ) THEN
          
          dyn_visc_c = muRef_Suth * ( T / Tref_Suth )**1.5_wp *                 &
               ( Tref_Suth + S_mu ) / ( T + S_mu )

          kin_visc_c = dyn_visc_c / rho_gas
       
       END IF
       
       ! Eq. 7 from Gueugneau et al, 2017 
       D_coeff = hydraulic_permeability / ( porosity * kin_visc_c * rho_gas *   &
            gas_compressibility )

       ! Equation 12 from Gueugneau et al, 2017
       IF ( ( REAL(exc_pore_pres) .GT. 0.0_wp )                                 &
            .AND. ( REAL(h) .GT. 0.0_wp) ) THEN

          f_inhibit = MAX(0.0_wp , 1.0_wp - ( SUM(alphas) /                     &
                CMPLX(maximum_solid_packing,0.0_wp,wp) ) )
          
          source_term(idx_poreEqn) = - rho_m * ( pi_g / 2.0_wp )**2 *           &
               D_coeff / MAX(h_threshold,h) * exc_pore_pres !* f_inhibit
          
       END IF
          
    END IF
    
    nh_term = source_term

    IF ( present(c_qj) .AND. present(c_nh_term_impl) ) THEN

       c_nh_term_impl = nh_term

    ELSEIF ( present(r_qj) .AND. present(r_nh_term_impl) ) THEN

       r_nh_term_impl = REAL( nh_term )

    END IF

    RETURN

  END SUBROUTINE eval_implicit_terms

  !******************************************************************************
  !> \brief Non-Hyperbolic semi-implicit terms
  !
  !> This subroutine evaluates the non-hyperbolic terms that are solved
  !> semi-implicitely by the solver. For example, any discontinuous term that
  !> appears in the friction terms.
  !> \date 20/01/2018
  !> \param[in]     grav3_surf         gravity correction 
  !> \param[in]     qcj                real conservative variables 
  !> \param[out]    nh_semi_impl_term  real non-hyperbolic terms
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE eval_nh_semi_impl_terms( Bprimej_x , Bprimej_y , Bsecondj_xx ,     &
       Bsecondj_xy , Bsecondj_yy , grav_coeff , qcj , qpj , nh_semi_impl_term , &
       Zj, fric_val )

    USE parameters_2D, ONLY: pore_pressure_flag
    
    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: Bprimej_x
    REAL(wp), INTENT(IN) :: Bprimej_y
    REAL(wp), INTENT(IN) :: Bsecondj_xx
    REAL(wp), INTENT(IN) :: Bsecondj_xy
    REAL(wp), INTENT(IN) :: Bsecondj_yy
    REAL(wp), INTENT(IN) :: grav_coeff

    REAL(wp), INTENT(IN) :: qcj(n_vars)
    REAL(wp), INTENT(IN) :: qpj(n_vars+2)
    REAL(wp), INTENT(IN) :: Zj ! value stochastic process

    REAL(wp), INTENT(OUT) :: nh_semi_impl_term(n_eqns)
    REAL(wp), INTENT(OUT) :: fric_val

    REAL(wp) :: source_term(n_eqns)

    REAL(wp) :: mod_vel
    REAL(wp) :: mod_hor_vel

    REAL(wp) :: h_threshold

    !--- Lahars rheology model variables

    !> Yield strenght (units: kg m-1 s-2)
    REAL(wp) :: tau_y

    !> Yield slope component of total friction (dimensionless)
    REAL(wp) :: s_y

    REAL(wp) :: r_h               !< real-value flow thickness
    REAL(wp) :: r_u               !< real-value x-velocity
    REAL(wp) :: r_v               !< real-value y-velocity
    REAL(wp) :: r_w               !< real_value z-velocity
    REAL(wp) :: r_alphas(n_solid) !< real-value solid volume fractions
    REAL(wp) :: r_rho_m           !< real-value mixture density [kg/m3]
    REAL(wp) :: r_T               !< real-value temperature [K]
    REAL(wp) :: r_alphal          !< real-value liquid volume fraction
    REAL(wp) :: r_alphag(n_add_gas) !< real-value add. gas volume fractions
    REAL(wp) :: r_red_grav
    REAL(wp) :: r_rho_c
    REAL(wp) :: r_Ri
    REAL(wp) :: Fr                !< Froude number
    REAL(wp) :: U_f               !< norm velocity for frirction
    REAL(wp) :: muF               !< mu(fr) o rmu(U)
    !REAL(wp) :: muU               !< mu(U)
    ! REAL(wp) :: Fr_x                !< Froude number
    ! REAL(wp) :: Fr_y                !< Froude number
    ! REAL(wp) :: mu_Fr_x             !< mu(fr)_x
    ! REAL(wp) :: mu_Fr_y             !< mu(fr)_y

    REAL(wp) :: temp_term
    REAL(wp) :: centr_force_term

    LOGICAL :: sp_heat_flag
    REAL(wp) :: r_sp_heat_c
    REAL(wp) :: r_sp_heat_mix

    REAL(wp) :: exc_pore_pres       !< excess pore pressure

    sp_heat_flag = .FALSE.


    ! initialize and evaluate the forces terms
    source_term(1:n_eqns) = 0.0_wp

    rheology_if:IF (rheology_flag) THEN

       r_h = qpj(1)
       r_u = qpj(idx_u)
       r_v = qpj(idx_v)

       IF ( alpha_flag ) THEN

          r_alphas(1:n_solid) = qpj(idx_alfas_first:idx_alfas_last)
          r_alphag(1:n_add_gas) = qpj(idx_addGas_first:idx_addGas_last)

       ELSE

          r_alphas(1:n_solid) = qpj(idx_alfas_first:idx_alfas_last) / qpj(1)
          r_alphag(1:n_add_gas) = qpj(idx_addGas_first:idx_addGas_last) / qpj(1)

       END IF

       r_alphal = 0.0_wp

       IF ( gas_flag .AND. liquid_flag ) THEN

          IF ( alpha_flag ) THEN

             r_alphal = qpj(n_vars)

          ELSE

             r_alphal = qpj(n_vars) / qpj(1)

          END IF

       END IF

       r_T = qpj(4)

       CALL mixt_var(qpj, r_Ri, r_rho_m, r_rho_c, r_red_grav, sp_heat_flag,     &
            r_sp_heat_c, r_sp_heat_mix)


       IF ( slope_correction_flag ) THEN
          
          r_w = r_u * Bprimej_x + r_v * Bprimej_y
          
       ELSE
          
          r_w = 0.0_wp
          
       END IF
       
       mod_vel = SQRT( r_u**2 + r_v**2 + r_w**2 )

       mod_hor_vel = SQRT( r_u**2 + r_v**2 )
       
       ! Voellmy Salm rheology
       rheology_model_if:IF ( rheology_model .EQ. 1 ) THEN

          IF ( mod_vel .GT. 0.0_wp ) THEN

             IF ( curvature_term_flag ) THEN

                ! See Eq. (3) Xia & Liang, 2018 Eng.Geol.   
                ! centrifugal force term: (u,v)^T*Hessian*(u,v)
                centr_force_term = Bsecondj_xx * r_u**2 + 2.0_wp * Bsecondj_xy  &
                     * r_u * r_v + Bsecondj_yy * r_v**2

             ELSE

                centr_force_term = 0.0_wp 

             END IF

             ! See Eq. (3,4) Xia & Liang, 2018 Eng.Geol.  
             ! add the contribution on mu (with coeff for large slope)
             ! and the contribution of centr. force (with coeff for slope)
             ! a = grav_coeff * ( r_red_grav + centr_force_term )
             ! 1/phi = SQRT(grav_coeff) 
             temp_term = mu * r_rho_m * ( r_red_grav + centr_force_term ) * r_h&
                  * SQRT(grav_coeff)

             temp_term = MAX(0.0_wp, temp_term)
             
             IF ( pore_pressure_flag ) THEN

                exc_pore_pres = qpj(idx_pore)

                ! See Eq. (2) Gueugneau et al. 2017, GRL
                ! add the contribution of pore pressure ( with coeff for slope)
                temp_term = MAX(0.0_wp, temp_term - mu * SQRT(grav_coeff)       &
                     * exc_pore_pres )

             END IF

             ! Friction terms cannot accelerate the flow
             ! this term is parallel to the full vel vector (u,v,w)
             ! tangential to the topography
             temp_term = MAX(0.0_wp,temp_term)

             ! horizontal terms
             ! units of dqc(2)/dt=d(rho h u)/dt (kg m-1 s-2)             
             source_term(2) = source_term(2) - temp_term * r_u / mod_vel
             ! units of dqc(3)/dt=d(rho h v)/dt (kg m-1 s-2)
             source_term(3) = source_term(3) - temp_term * r_v / mod_vel

          END IF

          ! Plastic rheology
       ELSEIF ( rheology_model .EQ. 2 ) THEN


          ! Temperature dependent rheology
       ELSEIF ( rheology_model .EQ. 3 ) THEN

          IF ( mod_vel .GT. 0.0_wp ) THEN

             ! units of dqc(2)/dt [kg m-1 s-2]
             source_term(2) = source_term(2) - tau0 * r_u / mod_vel

             ! units of dqc(3)/dt [kg m-1 s-2]
             source_term(3) = source_term(3) - tau0 * r_v / mod_vel

          END IF

          ! Lahars rheology (O'Brien 1993, FLO2D)
       ELSEIF ( rheology_model .EQ. 4 ) THEN

          h_threshold = 1.0E-20_wp

          ! Yield strength (units: kg m-1 s-2)
          tau_y = alpha2 * ( EXP( beta2 * SUM(r_alphas) ) - 1.0_wp )

          IF ( r_h .GT. h_threshold ) THEN

             ! Yield slope component (dimensionless)
             s_y = tau_y / ( grav * r_rho_m * r_h )

          ELSE

             ! Yield slope component (dimensionless)
             s_y = tau_y / ( grav * r_rho_m * h_threshold )

          END IF

          IF ( mod_vel .GT. 0.0_wp ) THEN

             temp_term = grav * r_rho_m * r_h * s_y

             ! units of dqc(2)/dt [kg m-1 s-2]
             source_term(2) = source_term(2) - temp_term * r_u / mod_hor_vel

             ! units of dqc(3)/dt [kg m-1 s-2]
             source_term(3) = source_term(3) - temp_term * r_v / mod_hor_vel

          END IF

       ELSEIF ( rheology_model .EQ. 9 ) THEN

          ! From Zhu et al. 2020
          ! (DOI: https://doi.org/10.1007/s10035-020-01053-7)
          ! mu_0: Coulomb friction coefficient at Fr=+inf
          ! mu_inf: Coulomb friction coefficient at Fr=0
          ! Fr : froude number (!!!computed here using the total velocity and
          ! the thickness instead of the particle holdup!!!)
          ! Fr_0 : Renormalization factor controlling the gradient of the
          ! function Zj : Value of the OU process at given cell (they can be
          ! tranformed)
          
          ! should also use rho or alpha?
          ! must add something for curvature or temperature?
                    
          ! Compute friction only if mass is flowing (this implies that there
          ! is mass)
 
          IF ( mod_vel .GT. 0.0_wp ) THEN 
             ! Computing froude number (The definition in Zhu 2020 et Roche 2021
             ! is sligtlhy different!)
             Fr = mod_vel / SQRT(r_red_grav * r_h) 
             
             ! Mofidy deterministic Fr if needed
             IF ( stochastic_flag ) THEN
               ! Add mean field correction or stochastic noise   
               IF (mean_field_flag) THEN
                  Fr = Fr !+ MeanFieldCorrection(grav, r_h, Fr) ! TEMP MODEL, DO NOT USE IT
               ELSE
                  IF (r_h .GT. 0.01_wp) THEN
                     Fr = Fr + Zj 
                  END IF
               END IF
               ! Set Foude number to zero if negative
               IF (Fr .LT. 0._wp) Fr = 0._wp
             END IF

             ! Evaluating mu(fr), (mu_0 < mu_inf), it is increasing with Fr
             muF = mu_inf + (mu_0 - mu_inf) * exp(-Fr/Fr_0)
             
             ! Compute temporal friction term (temp_term = 0 if h=0)            
             temp_term = r_rho_m *  ( muF * r_h * grav_coeff * ( r_red_grav +   &
                  centr_force_term ) )

             ! Update the source term
             ! units of dqc(2)/dt=d(rho h u)/dt (kg m-1 s-2)
             source_term(2) = source_term(2) - temp_term * r_u / mod_hor_vel
             
             ! units of dqc(3)/dt=d(rho h v)/dt (kg m-1 s-2)
             source_term(3) = source_term(3) - temp_term * r_v / mod_hor_vel
          
          ELSE ! If ||u|| = 0 then there will be no friction in this code
            
            muF = 0._wp
          
          END IF
         
         ! Set value friction to output variable (so that it can be saved!)
          fric_val = muF

      ELSEIF ( rheology_model .EQ. 10 ) THEN

          ! From Lucas et al. 2014 (DOI: 10.1038/ncomms4417)
          ! mu_0: Coulomb friction coefficient at U=0
          ! mu_inf: Coulomb friction coefficient at U=+inf
          ! U_w : Renormalization factor controlling the gradient of the function
          ! Zj : Value of the OU process at given cell (they can be tranformed)
                              
          ! Compute friction only if mass is flowing (this implies that there is mass)
          IF ( mod_vel .GT. 0.0_wp ) THEN 

             ! Mofidy deterministic U if needed
             IF ( stochastic_flag ) THEN
               ! Add mean field correction or stochastic noise   
               IF (mean_field_flag) THEN
                  U_f = mod_vel !TEMP MODEL, DO NOT USE IT
               ELSE
                  U_f = mod_vel + Zj 
               END IF
             ELSE
               U_f = mod_vel
             END IF

             ! Evaluating mu(U), (mu_0 > mu_inf), it is decreasing with U
             IF ( U_f .LT. U_w) THEN
               muF = mu_0
             ELSE
               ! (for small h ; u->inf but if not true)
               ! Max friction if h -> 0 : else weakening friction
               IF (r_h .LT. 0.1_wp) THEN 
                  muF = mu_0
               ELSE
                  muF = mu_inf + (mu_0-mu_inf)/(U_f/U_w)
               END IF
               ! Alternative :  
               ! Use exp and limit with sigmoid : (a=15, b=0.5)
               ! muF = mu_0 - ((mu_inf - mu_0) * exp(-(U_f-U_w)/U_w) *             &       
               !        (1._wp / (1._wp + exp(-15_wp * (r_h - 0.5_wp))))  
             END IF

             ! Compute temporal friction term (temp_term = 0 if h=0)            
             temp_term = r_rho_m *  ( muF * r_h * grav_coeff * ( r_red_grav +      &
                  centr_force_term ) )

             ! Update the source term
             ! units of dqc(2)/dt=d(rho h u)/dt (kg m-1 s-2)
             source_term(2) = source_term(2) - temp_term * r_u / mod_hor_vel
             
             ! units of dqc(3)/dt=d(rho h v)/dt (kg m-1 s-2)
             source_term(3) = source_term(3) - temp_term * r_v / mod_hor_vel
          
          ELSE ! If ||u|| = 0 then there will be no friction in this code
            
            muF = 0._wp
          
          END IF

          ! Set value friction to output variable (so that it can be saved!)
          fric_val = muF
          
       ENDIF rheology_model_if
       
    ENDIF rheology_if

    nh_semi_impl_term = source_term

    RETURN

  END SUBROUTINE eval_nh_semi_impl_terms

  !******************************************************************************
  !> \brief Mass exchange terms
  !
  !> This subroutine evaluates all the terms related to mass exhange between the
  !> flow and the environment: air entrainment, deposition, erosion, entrainment
  !> of water vapour. 
  !> \date 20/01/2018
  !> \param[in]     qpj                real physical variables
  !> \param[in]     B_zone             integer defining topo type (ground,water)
  !> \param[in]     B_prime_x          local topography x-slope
  !> \param[in]     B_prime_y          local topography y-slope
  !> \param[in]     erodible           erodible thickness of solid phases
  !> \param[in]     dt                 time step
  !> \param[out]    erosion_term       erosion rate of solid phases
  !> \param[out]    deposition_term    deposition rate of solid phases
  !> \param[out]    continuous_phase_erosion_term  
  !> \param[out]    continuous_phase_loss_term  
  !> \param[out]    eqns_term          terms associated with mass exchange
  !> \param[out]    topo_term          rate of change of topography
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE eval_mass_exchange_terms( qpj , B_zone , B_prime_x , B_prime_y ,   &
       erodible , dt , erosion_term , deposition_term ,                         &
       continuous_phase_erosion_term , continuous_phase_loss_term , eqns_term , &
       topo_term  )

    USE geometry_2d, ONLY : pi_g

    USE parameters_2d, ONLY : erodible_deposit_flag , liquid_vaporization_flag ,&
         vertical_profiles_flag , bottom_conc_flag , pore_pressure_flag ,       &
         gas_loss_flag

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: qpj(n_vars+2)              !< physical variables
    INTEGER, INTENT(IN) :: B_zone
    REAL(wp), INTENT(IN) :: B_prime_x
    REAL(wp), INTENT(IN) :: B_prime_y
    REAL(wp), INTENT(IN) :: erodible(n_solid)          !< erodible thickness
    REAL(wp), INTENT(IN) :: dt

    REAL(wp), INTENT(OUT) :: erosion_term(n_solid)     !< erosion term
    REAL(wp), INTENT(OUT) :: deposition_term(n_solid)  !< deposition term
    REAL(wp), INTENT(OUT) :: continuous_phase_erosion_term
    REAL(wp), INTENT(OUT) :: continuous_phase_loss_term
    REAL(wp), INTENT(OUT) :: eqns_term(n_eqns)
    REAL(wp), INTENT(OUT) :: topo_term


    REAL(wp) :: mod_vel
    REAL(wp) :: mod_vel2 

    REAL(wp) :: hind_exp
    REAL(wp) :: alpha_max

    INTEGER :: i_solid

    REAL(wp) :: r_h          !< real-value flow thickness
    REAL(wp) :: r_u          !< real-value x-velocity
    REAL(wp) :: r_v          !< real-value y-velocity
    REAL(wp) :: r_W          !< real-value z-velocity
    REAL(wp) :: r_alphas(n_solid) !< real-value solid volume fractions
    REAL(wp) :: r_alphag(n_add_gas) !< real-value add.gas volume fractions
    REAL(wp) :: r_rho_c      !< real-value carrier phase density [kg/m3]
    REAL(wp) :: r_T          !< real-value mixture temperature [K]
    REAL(wp) :: r_rho_m      !< real-value mixture density [kg/m3]
    REAL(wp) :: r_Zs         !< real-value stochastic variable
    REAL(wp) :: r_exc_pore_pres   !< real-value excess p pres

    REAL(wp) :: tot_erosion  !< total erosion rate [m/s]
    REAL(wp) :: tot_solid_erosion !< total solid erosion rate [m/s]

    REAL(wp) :: alphas_tot   !< total solid volume fraction

    REAL(wp) :: Tc           !< temperature of carrier pphase [K]

    REAL(wp) :: alpha1       !< viscosity of continuous phase [kg m-1 s-1]
    REAL(wp) :: fluid_visc   
    REAL(wp) :: inv_kin_visc
    REAL(wp) :: rhoc
    REAL(wp) :: expA , expB
    REAL(wp) :: r_Ri
    REAL(wp) :: r_red_grav

    !> Hindered settling velocity (units: m s-1 )
    REAL(wp) :: settling_vel

    LOGICAL :: sp_heat_flag
    REAL(wp) :: r_sp_heat_c
    REAL(wp) :: r_sp_heat_mix

    REAL(wp) :: entr_coeff
    REAL(wp) :: air_entr

    REAL(wp) :: dep_tot
    REAL(wp) :: ers_tot
    REAL(wp) :: rho_dep_tot
    REAL(wp) :: rho_ers_tot

    REAL(wp) :: T_liquid
    REAL(wp) :: T_boiling
    REAL(wp) :: sp_latent_heat
    REAL(wp) :: sp_heat_liq_water
    REAL(wp) :: mass_vap_rate

    REAL(wp) :: dep_coeff(n_solid)

    REAL(wp) :: pore_pressure_term
    REAL(wp) :: f_inhibit
    REAL(wp) :: vel_loss_gas

    REAL(wp) :: dyn_visc_c
    REAL(wp) :: rho_c

    erosion_term(1:n_solid) = 0.0_wp
    deposition_term(1:n_solid) = 0.0_wp
    continuous_phase_erosion_term = 0.0_wp
    continuous_phase_loss_term = 0.0_wp
    eqns_term(1:n_eqns) = 0.0_wp
    topo_term = 0.0_wp

    IF ( qpj(1) .LE. epsilon(1.0_wp) ) THEN

       RETURN

    END IF

    ! parameters for Michaels and Bolger (1962) sedimentation correction
    alpha_max = 0.6_wp
    hind_exp = 4.65_wp

    IF ( qpj(1) .LE. EPSILON(1.0_wp) ) RETURN

    r_h = qpj(1)
    r_u = qpj(idx_u)
    r_v = qpj(idx_v)

    IF ( alpha_flag ) THEN

       r_alphas(1:n_solid) = qpj(idx_alfas_first:idx_alfas_last)
       r_alphag(1:n_add_gas) = qpj(idx_addGas_first:idx_addGas_last)

    ELSE

       r_alphas(1:n_solid) = qpj(idx_alfas_first:idx_alfas_last) / qpj(1)
       r_alphag(1:n_add_gas) = qpj(idx_addGas_first:idx_addGas_last) / qpj(1)

    END IF

    alphas_tot = SUM(r_alphas)

    IF ( stoch_transport_flag ) r_Zs = qpj(idx_stoch)
    
    IF ( pore_pressure_flag ) r_exc_pore_pres = qpj(idx_pore)
    
    IF ( slope_correction_flag ) THEN

       r_w = r_u * B_prime_x + r_v * B_prime_y

    ELSE

       r_w = 0.0_wp

    END IF

    mod_vel2 = r_u**2 + r_v**2 + r_w**2
    mod_vel = SQRT( mod_vel2 )

    IF ( erosion_coeff .GT. 0.0_wp ) THEN

       ! empirical formulation (see Fagents & Baloga 2006, Eq. 5)
       ! here we use the solid volume fraction instead of relative density
       ! This term has units: m s-1    
       tot_erosion = erosion_coeff * mod_vel * r_h * ( 1.0_wp-alphas_tot )

       tot_solid_erosion = tot_erosion * ( 1.0_wp - erodible_porosity )

       erosion_term(1:n_solid) = erodible_fract(1:n_solid)  * tot_solid_erosion

    ELSE

       tot_solid_erosion = 0.0_wp

       erosion_term(1:n_solid) = 0.0_wp

    END IF

    ! Limit the deposition during a single time step
    erosion_term(1:n_solid) = MAX(0.0_wp,MIN( erosion_term(1:n_solid),          &
         erodible(1:n_solid) / dt ) )

    tot_solid_erosion = SUM( erosion_term(1:n_solid) )
    tot_erosion = tot_solid_erosion / ( 1.0_wp - erodible_porosity )

    continuous_phase_erosion_term = tot_erosion * erodible_porosity 

    IF ( alphas_tot .LE. alphastot_min ) RETURN

    r_T = qpj(4)

    sp_heat_flag = .TRUE.

    CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav,sp_heat_flag,r_sp_heat_c, &
         r_sp_heat_mix)

    IF ( rheology_model .EQ. 4 ) THEN

       ! alpha1 here has units: kg m-1 s-1
       ! in Table 2 from O'Brien 1988, the values reported have different
       ! units ( poises). 1poises = 0.1 kg m-1 s-1

       ! convert from Kelvin to Celsius
       Tc = r_T - 273.15_wp

       ! the dependance of viscosity on temperature is modeled with the
       ! equation presented at:
       ! https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781118131473.app3
       !
       ! In addition, we use a reference value provided in input at a 
       ! reference temperature. This value is used to scale the equation
       IF ( REAL(Tc) .LT. 20.0_wp ) THEN

          expA = 1301.0_wp / ( 998.333_wp + 8.1855_wp * ( Tc - 20.0_wp )        &
               + 0.00585_wp * ( Tc - 20.0_wp )**2 ) - 1.30223_wp

          alpha1 = alpha1_coeff * 1.0E-3_wp * 10.0_wp**expA

       ELSE

          expB = ( 1.3272_wp * ( 20.0_wp - Tc ) - 0.001053_wp *                 &
               ( Tc - 20.0_wp )**2 ) / ( Tc + 105.0_wp )

          alpha1 = alpha1_coeff * 1.002E-3_wp * 10.0_wp**expB 

       END IF

       ! Fluid dynamic viscosity [kg m-1 s-1]
       fluid_visc = alpha1 * EXP( beta1 * alphas_tot )
       ! Kinematic viscosity [m2 s-1]
       inv_kin_visc = r_rho_m / fluid_visc
       ! Continuous phase density used for the settling velocity
       rhoc = r_rho_m

    ELSE

       IF ( gas_flag .AND. sutherland_flag ) THEN
          
          dyn_visc_c = muRef_Suth * ( r_T / Tref_Suth )**1.5_wp *               &
               ( Tref_Suth + S_mu ) / ( r_T + S_mu )

          rho_c = pres / ( sp_gas_const_a * r_T )
          kin_visc_c = dyn_visc_c / rho_c
       
       END IF
       
       ! Viscosity read from input file [m2 s-1]
       inv_kin_visc = 1.0_wp / kin_visc_c
       ! Continuous phase density used for the settling velocity
       rhoc = r_rho_c

    END IF

    IF ( ( vertical_profiles_flag ) .AND. ( bottom_conc_flag ) ) THEN

       CALL eval_dep_coeffs( qpj , mod_vel , r_rho_c , r_rho_m , dep_coeff )

    ELSE

       dep_coeff(1:n_solid) = 1.0_wp

    END IF

    DO i_solid=1,n_solid

       IF ( ( r_alphas(i_solid) .GT. 0.0_wp ) .AND. ( settling_flag ) ) THEN

          settling_vel = settling_velocity( diam_s(i_solid) , rho_s(i_solid) ,  &
               rhoc , inv_kin_visc )

          deposition_term(i_solid) = r_alphas(i_solid) * settling_vel *         &
               dep_coeff(i_solid)

          IF ( rheology_model .NE. 4 ) THEN

             ! Michaels and Bolger (1962) sedimentation correction accounting 
             ! for hindered settling due to the presence of particles
             deposition_term(i_solid) = deposition_term(i_solid) *              &
                  ( 1.0_wp - MIN( 1.0_wp , alphas_tot / alpha_max ) )**hind_exp

          END IF

          ! limit the deposition (cannot remove more than particles present
          ! in the flow)
          deposition_term(i_solid) = MIN( deposition_term(i_solid) ,            &
               r_h * r_alphas(i_solid) / dt )

          IF ( deposition_term(i_solid) .LT. 0.0_wp ) THEN

             WRITE(*,*) 'eval_erosion_dep_term'
             WRITE(*,*) 'deposition_term(i_solid)',deposition_term(i_solid)
             READ(*,*)

          END IF

       END IF

    END DO

    IF ( liquid_flag ) THEN

       ! set the rate of loss of continuous phase 
       continuous_phase_loss_term = loss_rate

    ELSE

       continuous_phase_loss_term = 0.0_wp

    END IF

    ! add the loss associated with solid deposition
    continuous_phase_loss_term =  continuous_phase_loss_term +                  &
         coeff_porosity * SUM( deposition_term(1:n_solid) )

    IF ( pore_pressure_flag .AND. gas_loss_flag ) THEN
       
       IF ( alphas_tot .LT. maximum_solid_packing ) THEN 
          
          f_inhibit = MAX(0.0_wp, 1.0_wp - ( alphas_tot / maximum_solid_packing ) )

          IF ( gas_flag .AND. sutherland_flag ) THEN
             
             dyn_visc_c = muRef_Suth * ( r_T / Tref_Suth )**1.5_wp *            &
                  ( Tref_Suth + S_mu ) / ( r_T + S_mu )
             
             rho_c = pres / ( sp_gas_const_a * r_T )
             kin_visc_c = dyn_visc_c / rho_c
             
          END IF
                 
          vel_loss_gas =  hydraulic_permeability /                              &
               ( kin_visc_c * r_rho_c ) / MAX(1.e-5,r_h) * 0.5_wp * pi_g *      &
               r_exc_pore_pres

          pore_pressure_term = vel_loss_gas * f_inhibit
          
          continuous_phase_loss_term = continuous_phase_loss_term +             &
               pore_pressure_term
          
       END IF

    END IF   

    ! limit the loss accountaing for maximum solid packing
    continuous_phase_loss_term = MIN( continuous_phase_loss_term ,              &
         r_h * MAX( 0.0_wp , maximum_solid_packing - SUM(r_alphas) ) / dt )

    ! loss of continuous phase cannot 
    continuous_phase_loss_term = MIN( continuous_phase_loss_term ,              &
         ( r_h *  ( 1.0_wp - SUM(r_alphas) ) - coeff_porosity * ( r_h *         &
         SUM(r_alphas) - dt * SUM( deposition_term(1:n_solid) ) ) ) / dt )


    IF ( erodible_deposit_flag ) T_erodible = r_T

    IF ( entrainment_flag .AND.  ( r_h .GT. 0.0_wp ) .AND.                      &
         ( r_Ri .GT. 0.0_wp ) ) THEN

       entr_coeff = 0.075_wp / SQRT( 1.0_wp + 718.0_wp * r_Ri**2.4_wp )

       air_entr = entr_coeff * SQRT(mod_vel2)

    ELSE

       air_entr = 0.0_wp

    END IF

    eqns_term(1:n_eqns) = 0.0_wp

    ! solid total volume deposition rate [m s-1]
    dep_tot = SUM( deposition_term )
    ! solid total volume deposition rate [m s-1]
    ers_tot = SUM( erosion_term )

    ! solid total mass deposition rate [kg m-2 s-1]
    rho_dep_tot = DOT_PRODUCT( rho_s , deposition_term )
    ! solid total mass erosion rate [kg m-2 s-1]
    rho_ers_tot = DOT_PRODUCT( rho_s , erosion_term )
       
    ! total mass equation source term [kg m-2 s-1]:
    ! deposition, erosion and entrainment are considered
    eqns_term(1) = rho_a_amb * air_entr + rho_ers_tot - rho_dep_tot             &
         + rho_c_sub * continuous_phase_erosion_term                            &
         - r_rho_c * continuous_phase_loss_term

    ! x-momenutm equation source term [kg m-1 s-2]:
    ! only deposition contribute to change in momentum, erosion does not carry
    ! any momentum inside the flow
    eqns_term(2) = - r_u * ( rho_dep_tot + r_rho_c * continuous_phase_loss_term )

    ! y-momentum equation source term [kg m-1 s-2]:
    ! only deposition contribute to change in momentum, erosion does not carry
    ! any momentum inside the flow
    eqns_term(3) = - r_v * ( rho_dep_tot + r_rho_c * continuous_phase_loss_term )

    ! Temperature/Energy equation source term [kg s-3]:
    ! deposition, erosion and entrainment are considered    
    IF ( energy_flag ) THEN

       eqns_term(4) = - r_T * ( SUM( rho_s * sp_heat_s * deposition_term )      &
            + r_rho_c * r_sp_heat_c * continuous_phase_loss_term )              &
            - 0.5_wp * mod_vel2 * ( rho_dep_tot + r_rho_c *                     &
            continuous_phase_loss_term )                                        &
            + T_erodible * ( SUM( rho_s * sp_heat_s * erosion_term )            &
            + rho_c_sub * r_sp_heat_c * continuous_phase_erosion_term )         &
            + T_ambient * sp_heat_a * rho_a_amb * air_entr

    ELSE

       eqns_term(4) = - r_T * ( SUM( rho_s * sp_heat_s * deposition_term )      &
            + r_rho_c * r_sp_heat_c * continuous_phase_loss_term )              &
            + T_erodible * ( SUM( rho_s * sp_heat_s * erosion_term )            &
            + rho_c_sub * r_sp_heat_c * continuous_phase_erosion_term )         &
            + T_ambient * sp_heat_a * rho_a_amb * air_entr

    END IF

    ! solid phase mass equation source term [kg m-2 s-1]:
    ! due to solid erosion and deposition
    eqns_term(idx_solidEqn_first:idx_solidEqn_last) = rho_s(1:n_solid) *        &
         ( erosion_term(1:n_solid) - deposition_term(1:n_solid) )

    IF ( gas_flag .AND. liquid_vaporization_flag .AND. ( B_zone .NE. 0 ) ) THEN

       ! gamma_steam = 0.114_wp
       T_liquid = 290.0_wp
       T_boiling = 373.15_wp
       sp_latent_heat = 2264705.0_wp
       sp_heat_liq_water = 4184.0_wp

       mass_vap_rate = gamma_steam * MAX(0.0_wp , r_T-T_boiling) * SUM( rho_s * &
            sp_heat_s * deposition_term ) / ( sp_heat_liq_water * ( T_boiling - &
            T_liquid ) + sp_latent_heat )

       eqns_term(1) = eqns_term(1) + mass_vap_rate
       eqns_term(4) = eqns_term(4) + mass_vap_rate * sp_heat_g(1) * T_boiling
       eqns_term(5+n_solid) = eqns_term(5+n_solid) + mass_vap_rate

    END IF

    IF ( stoch_transport_flag) THEN
    
       ! Evaluate the equation term related to the noise transport equation
       ! (if transport flag is false nothing is done)
       eqns_term(idx_stochEqn) = eqns_term(1) * r_Zs

    END IF
       
    IF ( pore_pressure_flag ) THEN
       
       ! Equation for q1*exc_pore_press
       eqns_term(idx_poreEqn) = eqns_term(1) * r_exc_pore_pres

    END IF
       
    ! erodible layer thickness source terms [m s-1]:
    ! due to erosion and deposition of solid+continuous phase
    topo_term = ( dep_tot - ers_tot ) / ( 1.0_wp - erodible_porosity ) 

    RETURN

  END SUBROUTINE eval_mass_exchange_terms


  !******************************************************************************
  !> \brief Internal boundary source fluxes
  !
  !> This subroutine evaluates the source terms at the interfaces when an
  !> internal radial source is present, as for a base surge. The terms are 
  !> applied as boundary conditions, and thus they have the units of the 
  !> physical variable qp
  !> \date 2019/12/01
  !> \param[in]     time         time 
  !> \param[in]     vect_x       unit vector velocity x-component 
  !> \param[in]     vect_y       unit vector velocity y-component 
  !> \param[out]    source_bdry  source terms  
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE eval_source_bdry( time, vect_x , vect_y , source_bdry )

    USE parameters_2d, ONLY : h_source , vel_source , T_source , alphas_source ,&
         alphag_source , alphal_source , time_param

    USE geometry_2d, ONLY : pi_g

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: time
    REAL(wp), INTENT(IN) :: vect_x
    REAL(wp), INTENT(IN) :: vect_y
    REAL(wp), INTENT(OUT) :: source_bdry(n_vars+2)

    REAL(wp) :: t_rem
    REAL(wp) :: t_coeff

    IF ( time .GE. time_param(4) ) THEN

       ! The exponents of t_coeff are such that Ri does not depend on t_coeff
       source_bdry(1) = 0.0_wp
       source_bdry(2) = 0.0_wp
       source_bdry(3) = 0.0_wp
       source_bdry(4) = T_source
       source_bdry(idx_solidEqn_first:idx_solidEqn_last) = 0.0_wp

       IF ( gas_flag .AND. liquid_flag ) THEN

          source_bdry(n_vars) = alphal_source

       END IF

       source_bdry(idx_u) = 0.0_wp
       source_bdry(idx_v) = 0.0_wp

       RETURN

    END IF

    t_rem = MOD( time , time_param(1) )

    t_coeff = 0.0_wp

    IF ( time_param(3) .EQ. 0.0_wp ) THEN

       IF ( t_rem .LE. time_param(2) ) t_coeff = 1.0_wp

    ELSE

       IF ( t_rem .LT. time_param(3) ) THEN

          t_coeff = 0.5_wp * ( 1.0_wp - COS( pi_g * t_rem / time_param(3) ) )

       ELSEIF ( t_rem .LE. ( time_param(2) - time_param(3) ) ) THEN

          t_coeff = 1.0_wp

       ELSEIF ( t_rem .LE. time_param(2) ) THEN

          t_coeff = 0.5_wp * ( 1.0_wp + COS( pi_g * ( ( t_rem - time_param(2) ) &
               / time_param(3) + 1.0_wp ) ) )

       END IF

    END IF

    ! The exponents of t_coeff are such that Ri does not depend on t_coeff
    source_bdry(1) = t_coeff * h_source
    source_bdry(2) = t_coeff**1.5_wp * h_source * vel_source * vect_x
    source_bdry(3) = t_coeff**1.5_wp * h_source * vel_source * vect_y
    source_bdry(4) = T_source

    IF ( alpha_flag ) THEN

       source_bdry(idx_solidEqn_first:idx_solidEqn_last) =                      &
            alphas_source(1:n_solid)
       source_bdry(idx_addGasEqn_first:idx_addGasEqn_last) =                    &
            alphag_source(1:n_add_gas)

       IF ( gas_flag .AND. liquid_flag ) THEN

          source_bdry(n_vars) = alphal_source

       END IF

    ELSE

       source_bdry(idx_solidEqn_first:idx_solidEqn_last) =                      &
            t_coeff * h_source * alphas_source(1:n_solid)

       source_bdry(idx_addGasEqn_first:idx_addGasEqn_last) =                    &
            t_coeff * h_source * alphag_source(1:n_add_gas)

       IF ( gas_flag .AND. liquid_flag ) THEN

          source_bdry(n_vars) = t_coeff * h_source * alphal_source

       END IF

    END IF

    source_bdry(idx_u) = t_coeff**0.5_wp * vel_source * vect_x
    source_bdry(idx_v) = t_coeff**0.5_wp * vel_source * vect_y 

    RETURN

  END SUBROUTINE eval_source_bdry

  !------------------------------------------------------------------------------
  !> Settling velocity function
  !
  !> This subroutine compute the settling velocity of the particles, as a
  !> function of diameter, density of particles and carrier phase and viscosity.
  !> \date 2019/11/11
  !> \param[in]    diam          particle diameter      
  !> \param[in]    rhos          particle density
  !> \param[in]    rhoc          carrier phase density
  !> \param[in]    inv_kin_visc  reciprocal of kinetic viscosity
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !------------------------------------------------------------------------------

  REAL(wp) FUNCTION settling_velocity(diam,rhos,rhoc,inv_kin_visc)

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: diam          !< particle diameter [m]
    REAL(wp), INTENT(IN) :: rhos          !< particle density [kg/m3]
    REAL(wp), INTENT(IN) :: rhoc          !< carrier phase density [kg/m3]
    REAL(wp), INTENT(IN) :: inv_kin_visc  !< carrier phase viscosity reciprocal

    REAL(wp) :: Rey           !< Reynolds number
    REAL(wp) :: inv_sqrt_C_D  !< Reciprocal of sqrt of Drag coefficient

    ! loop variables
    INTEGER :: i              !< loop counter for iterative procedure    
    REAL(wp) :: const_part    !< term not changing in iterative procedure
    REAL(wp) :: inv_sqrt_C_D_old  !< previous iteration sqrt of drag coefficient
    REAL(wp) :: set_vel_old   !< previous iteration settling velocity

    ! INTEGER :: dig          !< order of magnitude of settling velocity

    inv_sqrt_C_D = 1.0_wp

    IF ( rhos .LE. rhoc ) THEN

       settling_velocity = 0.0_wp
       RETURN

    END IF

    const_part =  SQRT( 0.75_wp * ( rhos / rhoc - 1.0_wp ) * diam * grav )

    settling_velocity = const_part * inv_sqrt_C_D

    Rey = diam * settling_velocity * inv_kin_visc

    IF ( Rey .LE. 1000.0_wp ) THEN

       C_D_loop:DO i=1,20

          set_vel_old = settling_velocity
          inv_sqrt_C_D_old = inv_sqrt_C_D
          inv_sqrt_C_D = SQRT( Rey / ( 24.0_wp * ( 1.0_wp +                     &
               0.15_wp*Rey**(0.687_wp) ) ) )

          settling_velocity = const_part * inv_sqrt_C_D

          IF ( ABS( set_vel_old - settling_velocity ) / set_vel_old             &
               .LT. 1.0E-6_wp ) THEN

!!$             ! round to first three significative digits
!!$             dig = FLOOR(LOG10(set_vel_old))
!!$             settling_velocity = 10.0_wp**(dig-3)                            &
!!$                  * FLOOR( 10.0_wp**(-dig+3)*set_vel_old ) 

             RETURN

          END IF

          Rey = diam * settling_velocity * inv_kin_visc

       END DO C_D_loop

    END IF

    RETURN

  END FUNCTION settling_velocity

END MODULE constitutive_2d


