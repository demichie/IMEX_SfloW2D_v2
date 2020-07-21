!********************************************************************************
!> \brief Constitutive equations
!********************************************************************************
MODULE constitutive_2d

  USE parameters_2d, ONLY : wp, sp ,tolh
  USE parameters_2d, ONLY : n_eqns , n_vars , n_solid
  USE parameters_2d, ONLY : rheology_flag , rheology_model , energy_flag ,      &
       liquid_flag , gas_flag

  IMPLICIT none

  !> flag used for size of implicit non linear-system
  LOGICAL, ALLOCATABLE :: implicit_flag(:)

  !> flag to activate air entrainment
  LOGICAL :: entrainment_flag
  
  !> gravitational acceleration 
  REAL(wp) :: grav
  REAL(wp) :: inv_grav

  !> drag coefficients (Voellmy-Salm model)
  REAL(wp) :: mu
  REAL(wp) :: xi

  !> drag coefficients (B&W model)
  REAL(wp) :: friction_factor

  !> drag coefficients (plastic model)
  REAL(wp) :: tau

  !> evironment temperature [K]
  REAL(wp) :: T_env

  !> radiative coefficient
  REAL(wp) :: rad_coeff

  !> friction coefficient
  REAL(wp) :: frict_coeff

  !> reference temperature [K]
  REAL(wp) :: T_ref

  !> reference kinematic viscosity [m2/s]
  REAL(wp) :: nu_ref

  !> viscosity parameter [K-1] (b in Table 1 Costa & Macedonio, 2005)
  REAL(wp) :: visc_par

  !> velocity boundary layer fraction of total thickness
  REAL(wp) :: emme

  !> specific heat [J kg-1 K-1]
  REAL(wp) :: c_p

  !> atmospheric heat trasnfer coefficient [W m-2 K-1] (lambda in C&M, 2005)
  REAL(wp) :: atm_heat_transf_coeff

  !> fractional area of the exposed inner core (f in C&M, 2005)
  REAL(wp) :: exp_area_fract

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

  !> Kinematic viscosity of liquid (units: m2 s-1)
  REAL(wp) :: kin_visc_l

  !> Kinematic viscosity of carrier phase (units: m2 s-1)
  REAL(wp) :: kin_visc_c

  COMPLEX(wp) :: rho_c     !< Density of carrier phase in the mixture 
  COMPLEX(wp) :: inv_rho_c

  REAL(wp) :: r_rho_c         !< real-value carrier phase density [kg/m3]
  REAL(wp) :: r_inv_rho_c



  !> Dynamic pressure
  REAL(wp) :: p_dyn
  
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

  !> Hindered settling velocity (units: m s-1 )
  REAL(wp) :: settling_vel

  !> Minimum volume fraction of solids in the flow
  REAL(wp) :: alphastot_min
  
  !> erosion model coefficient  (units: m-1 )
  REAL(wp), ALLOCATABLE :: erosion_coeff

  !> erodible substrate solid relative volume fractions
  REAL(wp), ALLOCATABLE :: erodible_fract(:)

  !> erodible substrate porosity (we assume filled by continous phase)
  REAL(wp) :: erodible_porosity
  
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
    implicit_flag(4) = .FALSE.

    ! Solid volume fraction
    implicit_flag(5:4+n_solid) = .FALSE.

    n_nh = COUNT( implicit_flag )

    RETURN
    
  END SUBROUTINE init_problem_param

  !******************************************************************************
  !> \brief Physical variables
  !
  !> This subroutine evaluates from the conservative local variables qj
  !> the local physical variables  (\f$h,u,v,\alpha_s,\rho_m,T,\alpha_l \f$).
  !> \param[in]    r_qj     real conservative variables 
  !> \param[out]   r_h      real-value flow thickness 
  !> \param[out]   r_u      real-value flow x-velocity 
  !> \param[out]   r_v      real-value flow y-velocity
  !> \param[out]   r_alphas real-value solid volume fractions
  !> \param[out]   r_rho_m  real-value flow density
  !> \param[out]   r_T      real-value flow temperature 
  !> \param[out]   r_alphal real-value liquid volume fraction
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !> \date 2019/12/13
  !******************************************************************************

  SUBROUTINE r_phys_var(r_qj,r_h,r_u,r_v,r_alphas,r_rho_m,r_T,r_alphal)

    USE parameters_2d, ONLY : eps_sing , eps_sing4
    IMPLICIT none

    REAL(wp), INTENT(IN) :: r_qj(n_vars)       !< real-value conservative var
    REAL(wp), INTENT(OUT) :: r_h               !< real-value flow thickness
    REAL(wp), INTENT(OUT) :: r_u               !< real-value x-velocity
    REAL(wp), INTENT(OUT) :: r_v               !< real-value y-velocity
    REAL(wp), INTENT(OUT) :: r_alphas(n_solid) !< real-value solid volume fracts
    REAL(wp), INTENT(OUT) :: r_rho_m           !< real-value mixture density
    REAL(wp), INTENT(OUT) :: r_T               !< real-value temperature
    REAL(wp), INTENT(OUT) :: r_alphal          !< real-value liquid volume fract

    REAL(wp) :: r_inv_rhom
    REAL(wp) :: r_xs(n_solid)     !< real-value solid mass fraction
    REAL(wp) :: r_xs_tot

    REAL(wp) :: r_Ri            !< real-value Richardson number
    REAL(wp) :: r_xl            !< real-value liquid mass fraction
    REAL(wp) :: r_xc            !< real-value carrier phase mass fraction
    REAL(wp) :: r_alphac        !< real-value carrier phase volume fraction
    REAL(wp) :: r_red_grav      !< real-value reduced gravity
    REAL(wp) :: r_sp_heat_mix   !< Specific heat of mixture

    REAL(wp) :: inv_qj1
 
    ! compute solid mass fractions
    IF ( r_qj(1) .GT. eps_sing ) THEN

       inv_qj1 = 1.0_wp / r_qj(1)

       r_xs(1:n_solid) = r_qj(5:4+n_solid) * inv_qj1

    ELSE

       r_xs(1:n_solid) = 0.0_wp

    END IF

    r_xs_tot = SUM(r_xs)

    IF ( gas_flag .AND. liquid_flag ) THEN

       ! compute liquid mass fraction
       IF ( r_qj(1) .GT. eps_sing ) THEN

          r_xl = r_qj(n_vars) * inv_qj1

       ELSE

          r_xl = 0.0_wp

       END IF

       ! compute carrier phase (gas) mass fraction
       r_xc =  1.0_wp - r_xs_tot - r_xl

       ! specific heat of the mixutre: mass average of sp. heat pf phases
       r_sp_heat_mix = DOT_PRODUCT( r_xs(1:n_solid) , sp_heat_s(1:n_solid) )    &
            + r_xl * sp_heat_l + r_xc * sp_heat_c

    ELSE

       ! compute carrier phase (gas or liquid) mass fraction
       r_xc = 1.0_wp - r_xs_tot

       ! specific heaf of the mixutre: mass average of sp. heat pf phases
       r_sp_heat_mix = DOT_PRODUCT( r_xs(1:n_solid) , sp_heat_s(1:n_solid) )    &
            + r_xc * sp_heat_c

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
       r_rho_c =  pres / ( sp_gas_const_a * r_T )
       r_inv_rho_c = sp_gas_const_a * r_T * inv_pres
       sp_heat_c = sp_heat_a

    END IF

    IF ( gas_flag .AND. liquid_flag ) THEN

       r_inv_rhom = DOT_PRODUCT( r_xs(1:n_solid) , inv_rho_s(1:n_solid) )       &
            + r_xl * inv_rho_l + r_xc * r_inv_rho_c

       r_rho_m = 1.0_wp / r_inv_rhom

       r_alphal = r_xl * r_rho_m * inv_rho_l

    ELSE

       r_inv_rhom = DOT_PRODUCT( r_xs(1:n_solid) , inv_rho_s(1:n_solid) )       &
            + r_xc * r_inv_rho_c

       r_rho_m = 1.0_wp / r_inv_rhom

    END IF

    ! convert from mass fraction to volume fraction
    r_alphas(1:n_solid) = r_xs(1:n_solid) * r_rho_m * inv_rho_s(1:n_solid)

    ! convert from mass fraction to volume fraction
    r_alphac = r_xc * r_rho_m * r_inv_rho_c

    r_h = r_qj(1) * r_inv_rhom

    ! reduced gravity
    r_red_grav = ( r_rho_m - rho_a_amb ) * r_inv_rhom * grav

    ! velocity components
    IF ( r_qj(1) .GT. eps_sing ) THEN

       r_u = r_qj(2) * inv_qj1
       r_v = r_qj(3) * inv_qj1

    ELSE

       r_u = SQRT(2.0_wp) * r_qj(1) * r_qj(2) / SQRT( r_qj(1)**4 + eps_sing4 )
       r_v = SQRT(2.0_wp) * r_qj(1) * r_qj(3) / SQRT( r_qj(1)**4 + eps_sing4 )

    END IF

    ! Richardson number
    IF ( ( r_u**2 + r_v**2 ) .GT. 0.0_wp ) THEN

       r_Ri = r_red_grav * r_h / ( r_u**2 + r_v**2 )

    ELSE

       r_Ri = 0.0_wp

    END IF

    RETURN

  END SUBROUTINE r_phys_var

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
  !> \param[out]   red_grav  complex-value flow density
  !> \param[out]   alphas    complex-value solid volume fractions
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !> \date 2019/12/13
  !******************************************************************************

  SUBROUTINE c_phys_var( c_qj , h , u , v , T , rho_m , alphas , inv_rhom )

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
    COMPLEX(wp), INTENT(OUT) :: inv_rhom       !< 1/mixture density [kg-1 m3]

    COMPLEX(wp) :: xs(n_solid)             !< sediment mass fractions
    COMPLEX(wp) :: xs_tot                  !< sum of solid mass fraction
    COMPLEX(wp) :: Ri                      !< Richardson number
    COMPLEX(wp) :: xl                      !< liquid mass fraction
    COMPLEX(wp) :: xc                      !< carrier phase mass fraction
    COMPLEX(wp) :: alphal                  !< liquid volume fraction
    COMPLEX(wp) :: alphac                  !< carrier phase volume fraction
    COMPLEX(wp) :: sp_heat_mix             !< Specific heat of mixture
    COMPLEX(wp) :: inv_cqj1

    INTEGER :: i_solid

    ! compute solid mass fractions
    IF ( REAL(c_qj(1)) .GT. eps_sing ) THEN

       inv_cqj1 = 1.0_wp / c_qj(1)

       xs(1:n_solid) = c_qj(5:4+n_solid) * inv_cqj1

    ELSE

       inv_cqj1 = CMPLX(0.0_wp,0.0_wp,wp)
       xs(1:n_solid) = CMPLX(0.0_wp,0.0_wp,wp)

    END IF

    xs_tot = SUM(xs)

    ! compute carrier phase (gas or liquid) mass fraction
    xc = CMPLX(1.0_wp,0.0_wp,wp) - xs_tot

    ! specific heaf of the mixutre: mass average of sp. heat pf phases
    sp_heat_mix = DOT_PRODUCT( xs(1:n_solid) , sp_heat_s(1:n_solid) )        &
         + xc * sp_heat_c

    IF ( gas_flag .AND. liquid_flag ) THEN

       ! compute liquid mass fraction
       xl = c_qj(n_vars) * inv_cqj1

       ! compute carrier phase (gas) mass fraction
       xc =  xc - xl

       ! specific heaf of the mixutre: mass average of sp. heat pf phases
       sp_heat_mix = sp_heat_mix + xl * sp_heat_l

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
       inv_rho_c = sp_gas_const_a * T * inv_pres

    END IF

    inv_rhom = DOT_PRODUCT( xs(1:n_solid) , c_inv_rho_s(1:n_solid) )         &
         + xc * inv_rho_c

    IF ( gas_flag .AND. liquid_flag ) inv_rhom = inv_rhom + xl * inv_rho_l

    rho_m = 1.0_wp / inv_rhom

    ! convert from mass fraction to volume fraction
    alphas(1:n_solid) = rho_m * xs(1:n_solid) * c_inv_rho_s(1:n_solid)
    
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
  !> \brief Physical variables
  !
  !> This subroutine evaluates from the physical real-value local variables qpj, 
  !> all the (real-valued ) variables that define the physical state and that are
  !> needed to compute the explicit equations terms.
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

  SUBROUTINE mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qpj(n_vars+2) !< real-value physical variables
    REAL(wp), INTENT(OUT) :: r_Ri         !< real-value Richardson number
    REAL(wp), INTENT(OUT) :: r_rho_m      !< real-value mixture density [kg/m3]
    REAL(wp), INTENT(OUT) :: r_rho_c !< real-value carrier phase density [kg/m3]
    REAL(wp), INTENT(OUT) :: r_red_grav   !< real-value reduced gravity
    REAL(wp) :: r_u                       !< real-value x-velocity
    REAL(wp) :: r_v                       !< real-value y-velocity
    REAL(wp) :: r_h                       !< real-value flow thickness
    REAL(wp) :: r_alphas(n_solid)         !< real-value solid volume fractions
    REAL(wp) :: r_T                       !< real-value temperature [K]
    REAL(wp) :: r_alphal                  !< real-value liquid volume fraction

    REAL(wp) :: alphas_tot

    r_h = qpj(1)

    IF ( qpj(1) .LE. 0.0_wp ) THEN

       r_u = 0.0_wp
       r_v = 0.0_wp
       r_T = T_ambient
       r_alphas(1:n_solid) = 0.0_wp
       r_red_grav = 0.0_wp
       r_Ri = 0.0_wp
       r_rho_m = rho_a_amb
       IF ( gas_flag .AND. liquid_flag ) r_alphal = 0.0_wp

       RETURN

    END IF

    r_u = qpj(n_vars+1)
    r_v = qpj(n_vars+2)
    r_T = qpj(4)
    r_alphas(1:n_solid) = qpj(5:4+n_solid) 
    alphas_tot = SUM(r_alphas)

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
       r_rho_m = ( 1.0_wp - alphas_tot - r_alphal ) * r_rho_c                &
            + DOT_PRODUCT( r_alphas , rho_s ) + r_alphal * rho_l

    ELSE

       ! density of mixture of carrier phase and solids
       r_rho_m = ( 1.0_wp - alphas_tot ) * r_rho_c + DOT_PRODUCT( r_alphas , rho_s ) 

    END IF

    ! reduced gravity
    r_red_grav = ( r_rho_m - rho_a_amb ) / r_rho_m * grav

    

    ! Richardson number
    IF ( ( r_u**2 + r_v**2 ) .GT. 0.0_wp ) THEN

       r_Ri = MIN(1.D15,r_red_grav * r_h / ( r_u**2 + r_v**2 ))

    ELSE

       r_Ri = 1.D10

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
  !> - qp(n_vars+1) = \f$ u \f$
  !> - qp(n_vars+2) = \f$ v \f$
  !> .
  !> The physical variables are those used for the linear reconstruction at the
  !> cell interfaces. Limiters are applied to the reconstructed slopes.
  !> \param[in]     qc     local conservative variables 
  !> \param[out]    qp     local physical variables  
  !
  !> \date 2019/11/11
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE qc_to_qp(qc,qp)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qc(n_vars)
    REAL(wp), INTENT(OUT) :: qp(n_vars+2)

    REAL(wp) :: r_h               !< real-value flow thickness
    REAL(wp) :: r_u               !< real-value x-velocity
    REAL(wp) :: r_v               !< real-value y-velocity
    REAL(wp) :: r_alphas(n_solid) !< real-value solid volume fractions
    REAL(wp) :: r_rho_m           !< real-value mixture density [kg/m3]
    REAL(wp) :: r_T               !< real-value temperature [K]
    REAL(wp) :: r_alphal          !< real-value liquid volume fraction

    
    CALL r_phys_var(qc,r_h,r_u,r_v,r_alphas,r_rho_m,r_T,r_alphal)
    
    qp(1) = r_h
    
    qp(2) = r_h*r_u
    qp(3) = r_h*r_v
    
    qp(4) = r_T
    qp(5:4+n_solid) = r_alphas(1:n_solid)
    
    IF ( gas_flag .AND. liquid_flag ) qp(n_vars) = r_alphal
    
    qp(n_vars+1) = r_u
    qp(n_vars+2) = r_v

    p_dyn = 0.5_wp * r_rho_m * ( r_u**2 + r_v**2 )
    
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
  !> - qp(5:4+n_s) = \f$ alphas(1:n_s) \f$
  !> - qp(n_vars) = \f$ alphal \f$
  !> - qp(n_vars+1) = \f$ u \f$
  !> - qp(n_vars+2) = \f$ v \f$
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
    REAL(wp) :: r_xl              !< real-value liquid mass fraction
    REAL(wp) :: r_xc              !< real-value carrier phase mass fraction
    REAL(wp) :: r_T               !< real-value temperature [K]
    REAL(wp) :: r_alphal          !< real-value liquid volume fraction
    REAL(wp) :: r_alphac          !< real-value carrier phase volume fraction
    REAL(wp) :: r_rho_m           !< real-value mixture density [kg/m3]
    REAL(wp) :: r_rho_c           !< real-value carrier phase density [kg/m3]
    REAL(wp) :: r_xs(n_solid)     !< real-value solid mass fraction

    REAL(wp) :: alphas_tot

    r_h = qp(1)

    IF ( r_h .GT. 0.0_wp ) THEN

       r_hu = qp(2)
       r_hv = qp(3)

       r_u = qp(n_vars+1)
       r_v = qp(n_vars+2)

    ELSE

       r_hu = 0.0_wp
       r_hv = 0.0_wp

       r_u = 0.0_wp
       r_v = 0.0_wp

       qc(1:n_vars) = 0.0_wp
       RETURN

    END IF

    r_T  = qp(4)

    r_alphas(1:n_solid) = qp(5:4+n_solid)
    alphas_tot = SUM(r_alphas)
 
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
       r_rho_m = r_alphac * r_rho_c + DOT_PRODUCT( r_alphas , rho_s )           &
            + r_alphal * rho_l

       ! liquid mass fraction
       r_xl = r_alphal * rho_l / r_rho_m

       ! solid mass fractions
       r_xs(1:n_solid) = r_alphas(1:n_solid) * rho_s(1:n_solid) / r_rho_m

       ! carrier (gas) mass fraction
       r_xc = r_alphac * r_rho_c / r_rho_m

       ! mass averaged mixture specific heat
       r_sp_heat_mix =  DOT_PRODUCT( r_xs , sp_heat_s ) + r_xl * sp_heat_l      &
            + r_xc * sp_heat_c

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
       r_rho_m = r_alphac * r_rho_c + DOT_PRODUCT( r_alphas , rho_s ) 

       ! solid mass fractions
       r_xs(1:n_solid) = r_alphas(1:n_solid) * rho_s(1:n_solid) / r_rho_m

       ! carrier (gas or liquid) mass fraction
       r_xc = r_alphac * r_rho_c / r_rho_m

       ! mass averaged mixture specific heat
       r_sp_heat_mix =  DOT_PRODUCT( r_xs , sp_heat_s ) + r_xc * sp_heat_c

    END IF

    qc(1) = r_rho_m * r_h 
    qc(2) = r_rho_m * r_hu
    qc(3) = r_rho_m * r_hv

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

    qc(5:4+n_solid) = r_h * r_alphas(1:n_solid) * rho_s(1:n_solid)

    IF ( gas_flag .AND. liquid_flag ) qc(n_vars) = r_h * r_alphal * rho_l

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

  SUBROUTINE eval_local_speeds_x(qpj,vel_min,vel_max)

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qpj(n_vars+2)

    REAL(wp), INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)

    REAL(wp) :: r_h          !< real-value flow thickness [m]
    REAL(wp) :: r_u          !< real-value x-velocity [m s-1]
    REAL(wp) :: r_v          !< real-value y-velocity [m s-1]
    REAL(wp) :: r_Ri         !< real-value Richardson number
    REAL(wp) :: r_rho_m      !< real-value mixture density [kg m-3]
    REAL(wp) :: r_rho_c      !< real-value carrier phase density [kg m-3]
    REAL(wp) :: r_red_grav   !< real-value reduced gravity [m s-2]

    CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav)

    r_h = qpj(1)
    r_u = qpj(n_vars+1)
    r_v = qpj(n_vars+2)

    IF ( r_red_grav * r_h .LT. 0.0_wp ) THEN

       vel_min(1:n_eqns) = r_u
       vel_max(1:n_eqns) = r_u

    ELSE

       vel_min(1:n_eqns) = r_u - SQRT( r_red_grav * r_h )
       vel_max(1:n_eqns) = r_u + SQRT( r_red_grav * r_h )

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

    REAL(wp), INTENT(IN)  :: qpj(n_vars+2)
    REAL(wp), INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)

    REAL(wp) :: r_h          !< real-value flow thickness
    REAL(wp) :: r_u          !< real-value x-velocity
    REAL(wp) :: r_v          !< real-value y-velocity
    REAL(wp) :: r_Ri         !< real-value Richardson number
    REAL(wp) :: r_rho_m      !< real-value mixture density [kg/m3]
    REAL(wp) :: r_rho_c      !< real-value carrier phase density [kg/m3]
    REAL(wp) :: r_red_grav   !< real-value reduced gravity

    CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav)

    r_h = qpj(1)
    r_u = qpj(n_vars+1)
    r_v = qpj(n_vars+2)

    IF ( r_red_grav * r_h .LT. 0.0_wp ) THEN

       vel_min(1:n_eqns) = r_v
       vel_max(1:n_eqns) = r_v

    ELSE

       vel_min(1:n_eqns) = r_v - SQRT( r_red_grav * r_h )
       vel_max(1:n_eqns) = r_v + SQRT( r_red_grav * r_h )

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

  SUBROUTINE eval_fluxes(qcj,qpj,dir,flux)

    USE parameters_2d, ONLY : eps_sing

    IMPLICIT none

    REAL(wp), INTENT(IN) :: qcj(n_vars)
    REAL(wp), INTENT(IN) :: qpj(n_vars+2)
    INTEGER, INTENT(IN) :: dir

    REAL(wp), INTENT(OUT) :: flux(n_eqns)

    REAL(wp) :: r_h          !< real-value flow thickness [m]
    REAL(wp) :: r_u          !< real-value x-velocity [m s-1]
    REAL(wp) :: r_v          !< real-value y-velocity [m s-1]
    REAL(wp) :: r_Ri         !< real-value Richardson number
    REAL(wp) :: r_rho_m      !< real-value mixture density [kg m-3]
    REAL(wp) :: r_rho_c      !< real-value carrier phase density [kg m-3]
    REAL(wp) :: r_red_grav   !< real-value reduced gravity [m s-2]

    pos_thick:IF ( qpj(1) .GT. eps_sing ) THEN

       r_h = qpj(1)
       r_u = qpj(n_vars+1)
       r_v = qpj(n_vars+2)

       CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav)

       IF ( dir .EQ. 1 ) THEN

          ! Mass flux in x-direction: u * ( rhom * h )
          flux(1) = r_u * qcj(1)

          ! x-momentum flux in x-direction + hydrostatic pressure term
          flux(2) = r_u * qcj(2) + 0.5_wp * r_rho_m * r_red_grav * r_h**2

          ! y-momentum flux in x-direction: u * ( rho * h * v )
          flux(3) = r_u * qcj(3)

          IF ( energy_flag ) THEN

             ! ENERGY flux in x-direction
             flux(4) = r_u * ( qcj(4) + 0.5_wp * r_rho_m * r_red_grav * r_h**2 )

          ELSE

             ! Temperature flux in x-direction: u * ( h * T )
             flux(4) = r_u * qcj(4)

          END IF

          ! Mass flux of solid in x-direction: u * ( h * alphas * rhos )
          flux(5:4+n_solid) = r_u * qcj(5:4+n_solid)

          ! Solid flux can't be larger than total flux
          IF ( ( flux(1) .GT. 0.0_wp ) .AND. ( SUM(flux(5:4+n_solid)) / flux(1) &
               .GT. 1.0_wp ) ) THEN

             flux(5:4+n_solid) = &
                  flux(5:4+n_solid) / SUM(flux(5:4+n_solid)) * flux(1)

          END IF

          ! Mass flux of liquid in x-direction: u * ( h * alphal * rhol )
          IF ( gas_flag .AND. liquid_flag ) flux(n_vars) = r_u * qcj(n_vars)

       ELSEIF ( dir .EQ. 2 ) THEN

          ! flux G (derivated wrt y in the equations)
          flux(1) = r_v * qcj(1)

          flux(2) = r_v * qcj(2)

          flux(3) = r_v * qcj(3) + 0.5_wp * r_rho_m * r_red_grav * r_h**2

          IF ( energy_flag ) THEN

             ! ENERGY flux in x-direction
             flux(4) = r_v * ( qcj(4) + 0.5_wp * r_rho_m * r_red_grav * r_h**2 )

          ELSE

             ! Temperature flux in y-direction
             flux(4) = r_v * qcj(4)

          END IF

          ! Mass flux of solid in y-direction: v * ( h * alphas * rhos )
          flux(5:4+n_solid) = r_v * qcj(5:4+n_solid)

          ! Solid flux can't be larger than total flux
          IF ( ( flux(1) .GT. 0.0_wp ) .AND. ( SUM(flux(5:4+n_solid)) / flux(1) &
               .GT. 1.0_wp ) ) THEN

             flux(5:4+n_solid) = &
                  flux(5:4+n_solid) / SUM(flux(5:4+n_solid)) * flux(1)

          END IF

          ! Mass flux of liquid in x-direction: u * ( h * alphal * rhol )
          IF ( gas_flag .AND. liquid_flag ) flux(n_vars) = r_v * qcj(n_vars)

       END IF

    ELSE

       flux(1:n_eqns) = 0.0_wp

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

    USE parameters_2d, ONLY : four_thirds , neg_four_thirds

    IMPLICIT NONE

    COMPLEX(wp), INTENT(IN), OPTIONAL :: c_qj(n_vars)
    COMPLEX(wp), INTENT(OUT), OPTIONAL :: c_nh_term_impl(n_eqns)
    REAL(wp), INTENT(IN), OPTIONAL :: r_qj(n_vars)
    REAL(wp), INTENT(OUT), OPTIONAL :: r_nh_term_impl(n_eqns)

    COMPLEX(wp) :: h                       !< height [m]
    COMPLEX(wp) :: inv_h                   !< 1/height [m-1]
    COMPLEX(wp) :: u                       !< velocity (x direction) [m/s]
    COMPLEX(wp) :: v                       !< velocity (y direction) [m/s]
    COMPLEX(wp) :: T                       !< temperature [K]
    COMPLEX(wp) :: rho_m                   !< mixture density [kg/m3]
    COMPLEX(wp) :: alphas(n_solid)         !< sediment volume fractions
    COMPLEX(wp) :: inv_rho_m               !< 1/mixture density [kg-1 m3]
 
    COMPLEX(wp) :: qj(n_vars)
    COMPLEX(wp) :: nh_term(n_eqns)
    COMPLEX(wp) :: forces_term(n_eqns)

    COMPLEX(wp) :: mod_vel
    COMPLEX(wp) :: mod_vel2
    COMPLEX(wp) :: gamma
    REAL(wp) :: h_threshold

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

    ! initialize and evaluate the forces terms
    forces_term(1:n_eqns) = CMPLX(0.0_wp,0.0_wp,wp)

    IF (rheology_flag) THEN

       CALL c_phys_var(qj,h,u,v,T,rho_m,alphas,inv_rho_m)

       mod_vel2 = u**2 + v**2 
       mod_vel = SQRT( mod_vel2 )

       IF ( rheology_model .EQ. 1 ) THEN
          ! Voellmy Salm rheology

          IF ( REAL(mod_vel) .NE. 0.0_wp ) THEN 

             ! IMPORTANT: grav3_surf is always negative 
             forces_term(2) = forces_term(2) - rho_m * ( u / mod_vel ) *        &
                  ( grav / xi ) * mod_vel2

             forces_term(3) = forces_term(3) - rho_m * ( v / mod_vel ) *        &
                  ( grav / xi ) * mod_vel2

          ENDIF

       ELSEIF ( rheology_model .EQ. 2 ) THEN

          ! Plastic rheology
          IF ( REAL(mod_vel) .NE. 0.0_wp ) THEN 

             forces_term(2) = forces_term(2) - rho_m * tau * (u/mod_vel)

             forces_term(3) = forces_term(3) - rho_m * tau * (v/mod_vel)

          ENDIF

       ELSEIF ( rheology_model .EQ. 3 ) THEN

          h_threshold = 1.0E-10_wp

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
             forces_term(2) = forces_term(2) - rho_m * gamma * u

             ! Last R.H.S. term in equation 3 from Costa & Macedonio, 2005
             forces_term(3) = forces_term(3) - rho_m * gamma * v

          ENDIF

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
             ! s_v = Kappa * fluid_visc * mod_vel / ( 8.0_wp * rho_m * grav *h**2 )
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

             temp_term = grav * rho_m * h * s_f / mod_vel

             ! same units of dqc(2)/dt: kg m-1 s-2
             forces_term(2) = forces_term(2) - u * temp_term

             ! same units of dqc(3)/dt: kg m-1 s-2
             forces_term(3) = forces_term(3) - v * temp_term

          END IF

       ELSEIF ( rheology_model .EQ. 5 ) THEN

          tau = 1.0E-3_wp / ( 1.0_wp + 10.0_wp * h ) * mod_vel

          IF ( REAL(mod_vel) .NE. 0.0_wp ) THEN

             forces_term(2) = forces_term(2) - rho_m * tau * ( u / mod_vel )
             forces_term(3) = forces_term(3) - rho_m * tau * ( v / mod_vel )

          END IF


       ELSEIF ( rheology_model .EQ. 6 ) THEN

          IF ( REAL(mod_vel) .NE. 0.0_wp ) THEN 

             forces_term(2) = forces_term(2) - rho_m * u * friction_factor *    &
                  mod_vel

             forces_term(3) = forces_term(3) - rho_m * v * friction_factor *    &
                  mod_vel

          ENDIF

       ENDIF

    ENDIF

    nh_term = forces_term

    IF ( present(c_qj) .AND. present(c_nh_term_impl) ) THEN

       c_nh_term_impl = nh_term

    ELSEIF ( present(r_qj) .AND. present(r_nh_term_impl) ) THEN

       r_nh_term_impl = REAL( nh_term )

    END IF

    RETURN

  END SUBROUTINE eval_nonhyperbolic_terms

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

  SUBROUTINE eval_nh_semi_impl_terms( grav3_surf , qcj , nh_semi_impl_term )

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: grav3_surf

    REAL(wp), INTENT(IN) :: qcj(n_vars)
    REAL(wp), INTENT(OUT) :: nh_semi_impl_term(n_eqns)

    REAL(wp) :: forces_term(n_eqns)

    REAL(wp) :: mod_vel

    REAL(wp) :: h_threshold

    !--- Lahars rheology model variables

    !> Yield strenght (units: kg m-1 s-2)
    REAL(wp) :: tau_y

    !> Yield slope component of total friction (dimensionless)
    REAL(wp) :: s_y

    REAL(wp) :: r_h               !< real-value flow thickness
    REAL(wp) :: r_u               !< real-value x-velocity
    REAL(wp) :: r_v               !< real-value y-velocity
    REAL(wp) :: r_alphas(n_solid) !< real-value solid volume fractions
    REAL(wp) :: r_rho_m           !< real-value mixture density [kg/m3]
    REAL(wp) :: r_T               !< real-value temperature [K]
    REAL(wp) :: r_alphal          !< real-value liquid volume fraction


    ! initialize and evaluate the forces terms
    forces_term(1:n_eqns) = 0.0_wp

    IF (rheology_flag) THEN

       CALL r_phys_var(qcj,r_h,r_u,r_v,r_alphas,r_rho_m,r_T,r_alphal)
       
       ! Voellmy Salm rheology
       IF ( rheology_model .EQ. 1 ) THEN

          mod_vel = SQRT( r_u**2 + r_v**2 )
          
          IF ( mod_vel .GT. 0.0_wp ) THEN

             ! units of dqc(2)/dt=d(rho h v)/dt (kg m-1 s-2)
             forces_term(2) = forces_term(2) - r_rho_m * ( r_u / mod_vel ) *    &
                  mu * r_h * ( - grav * grav3_surf )

             ! units of dqc(3)/dt=d(rho h v)/dt (kg m-1 s-2)
             forces_term(3) = forces_term(3) - r_rho_m * ( r_v / mod_vel ) *    &
                  mu * r_h * ( - grav * grav3_surf )

          END IF

          ! Plastic rheology
       ELSEIF ( rheology_model .EQ. 2 ) THEN


          ! Temperature dependent rheology
       ELSEIF ( rheology_model .EQ. 3 ) THEN


          ! Lahars rheology (O'Brien 1993, FLO2D)
       ELSEIF ( rheology_model .EQ. 4 ) THEN

          h_threshold = 1.0E-20_wp

          ! Yield strength (units: kg m-1 s-2)
          tau_y = alpha2 * EXP( beta2 * SUM(r_alphas) )

          IF ( r_h .GT. h_threshold ) THEN

             ! Yield slope component (dimensionless)
             s_y = tau_y / ( grav * r_rho_m * r_h )

          ELSE

             ! Yield slope component (dimensionless)
             s_y = tau_y / ( grav * r_rho_m * h_threshold )

          END IF

          mod_vel = SQRT( r_u**2 + r_v**2 )
          
          IF ( mod_vel .GT. 0.0_wp ) THEN

             ! units of dqc(2)/dt [kg m-1 s-2]
             forces_term(2) = forces_term(2) - grav * r_rho_m * r_h *           &
                  ( r_u / mod_vel ) * s_y

             ! units of dqc(3)/dt [kg m-1 s-2]
             forces_term(3) = forces_term(3) - grav * r_rho_m * r_h *           &
                  ( r_v / mod_vel ) * s_y

          END IF

       ELSEIF ( rheology_model .EQ. 5 ) THEN

       ENDIF

    ENDIF

    nh_semi_impl_term = forces_term

    RETURN

  END SUBROUTINE eval_nh_semi_impl_terms

  !******************************************************************************
  !> \brief Explicit Forces term
  !
  !> This subroutine evaluates the non-hyperbolic terms to be treated explicitely
  !> in the DIRK numerical scheme (e.g. gravity,source of mass). The sign of the
  !> terms is taken with the terms on the left-hand side of the equations.
  !> \date 2019/12/13
  !> \param[in]     B_primej_x         local x-slope
  !> \param[in]     B_primej_y         local y_slope
  !> \param[in]     source_xy          local source
  !> \param[in]     qpj                physical variables 
  !> \param[in]     qcj                conservative variables 
  !> \param[out]    expl_term          explicit term
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE eval_expl_terms( Bprimej_x, Bprimej_y, source_xy, qpj, expl_term , &
       time, cell_fract_jk )

    USE parameters_2d, ONLY : h_source , vel_source , T_source , alphas_source ,&
         alphal_source , time_param , bottom_radial_source_flag

    
    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: Bprimej_x
    REAL(wp), INTENT(IN) :: Bprimej_y
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
    REAL(wp) :: r_red_grav   !< real-value reduced gravity
    REAL(wp) :: r_alphas(n_solid)
    REAL(wp) :: r_xs(n_solid)
    REAL(wp) :: r_alphal
    REAL(wp) :: r_xl    
    REAL(wp) :: r_alphac
    REAL(wp) :: r_xc
    
    REAL(wp) :: alphas_tot
    REAL(wp) :: sum_sl
    REAL(wp) :: r_sp_heat_mix
    
    REAL(wp) :: t_rem
    REAL(wp) :: t_coeff
    REAL(wp) :: pi_g
    REAL(wp) :: h_dot

    
    expl_term(1:n_eqns) = 0.0_wp

    IF ( ( qpj(1) .LE. 0.0_wp ) .AND. ( cell_fract_jk .EQ. 0.0_wp ) ) RETURN

    r_h = qpj(1)
    r_u = qpj(n_vars+1)
    r_v = qpj(n_vars+2)

    CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav)

    ! units of dqc(2)/dt [kg m-1 s-2]
    expl_term(2) = r_red_grav * r_rho_m * r_h * Bprimej_x

    ! units of dqc(3)/dt [kg m-1 s-2]
    expl_term(3) = r_red_grav * r_rho_m * r_h * Bprimej_y

    IF ( energy_flag ) THEN

       expl_term(4) = r_red_grav * r_rho_m * r_h * ( r_u * Bprimej_x            &
            + r_v * Bprimej_y )  

    ELSE

       expl_term(4) = 0.0_wp

    END IF
    
    IF ( ( time .GE. time_param(4) ) .OR. ( .NOT.bottom_radial_source_flag) ) THEN

       RETURN

    END IF

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
          
          t_coeff = ( t_rem - time_param(2) + time_param(3) ) / time_param(3)
          
       ELSE
          
          t_coeff = 0.0_wp
          
       END IF

    END IF
    
    h_dot = -cell_fract_jk * vel_source

    r_alphas(1:n_solid) = alphas_source(1:n_solid) 
    alphas_tot = SUM(r_alphas)
        
    IF ( gas_flag ) THEN

       ! carrier phase is gas
       r_rho_c = pres / ( sp_gas_const_a * t_source )
       sp_heat_c = sp_heat_a

    ELSE

       ! carrier phase is liquid
       r_rho_c = rho_l
       sp_heat_c = sp_heat_l

    END IF

    IF ( gas_flag .AND. liquid_flag ) THEN

       ! mixture of gas, liquid and solid
       r_alphal = qpj(n_vars)

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
       r_rho_m = r_alphac * r_rho_c + DOT_PRODUCT( r_alphas , rho_s )           &
            + r_alphal * rho_l

       ! liquid mass fraction
       r_xl = r_alphal * rho_l / r_rho_m

       ! solid mass fractions
       r_xs(1:n_solid) = r_alphas(1:n_solid) * rho_s(1:n_solid) / r_rho_m

       ! carrier (gas) mass fraction
       r_xc = r_alphac * r_rho_c / r_rho_m

       ! mass averaged mixture specific heat
       r_sp_heat_mix =  DOT_PRODUCT( r_xs , sp_heat_s ) + r_xl * sp_heat_l      &
            + r_xc * sp_heat_c

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
       r_rho_m = r_alphac * r_rho_c + DOT_PRODUCT( r_alphas , rho_s ) 

       ! solid mass fractions
       r_xs(1:n_solid) = r_alphas(1:n_solid) * rho_s(1:n_solid) / r_rho_m

       ! carrier (gas or liquid) mass fraction
       r_xc = r_alphac * r_rho_c / r_rho_m

       ! mass averaged mixture specific heat
       r_sp_heat_mix =  DOT_PRODUCT( r_xs , sp_heat_s ) + r_xc * sp_heat_c

    END IF
    
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
       
    expl_term(5:4+n_solid) = expl_term(5:4+n_solid) + t_coeff                     &
         * h_dot * alphas_source(1:n_solid) * rho_s(1:n_solid)

    IF ( gas_flag .AND. liquid_flag ) THEN

       expl_term(n_vars) = expl_term(n_vars) + t_coeff * h_dot * alphal_source  &
            * rho_l

    END IF

    
    RETURN

  END SUBROUTINE eval_expl_terms

  !******************************************************************************
  !> \brief Erosion/Deposition term
  !
  !> This subroutine evaluates the deposition term.
  !> \date 03/010/2018
  !> \param[in]     qpj                local physical variables 
  !> \param[in]     dt                 time step
  !> \param[out]    erosion_term       erosion term for each solid phase
  !> \param[out]    dep_term           deposition term for each solid phase
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE eval_erosion_dep_term( qpj , dt , erosion_term , deposition_term )

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: qpj(n_vars+2)              !< physical variables 
    REAL(wp), INTENT(IN) :: dt

    REAL(wp), INTENT(OUT) :: erosion_term(n_solid)     !< erosion term
    REAL(wp), INTENT(OUT) :: deposition_term(n_solid)  !< deposition term

    REAL(wp) :: mod_vel

    REAL(wp) :: hind_exp
    REAL(wp) :: alpha_max

    INTEGER :: i_solid

    REAL(wp) :: r_h          !< real-value flow thickness
    REAL(wp) :: r_u          !< real-value x-velocity
    REAL(wp) :: r_v          !< real-value y-velocity
    REAL(wp) :: r_alphas(n_solid) !< real-value solid volume fractions
    REAL(wp) :: r_rho_c      !< real-value carrier phase density [kg/m3]
    REAL(wp) :: r_T
    REAL(wp) :: r_rho_m
    REAL(wp) :: tot_solid_erosion

    REAL(wp) :: alphas_tot
    
    REAL(wp) :: Tc
    
    REAL(wp) :: alpha1
    REAL(wp) :: fluid_visc
    REAL(wp) :: kin_visc
    REAL(wp) :: rhoc
    REAL(wp) :: expA , expB
    
    ! parameters for Michaels and Bolger (1962) sedimentation correction
    alpha_max = 0.6_wp
    hind_exp = 4.65_wp

    deposition_term(1:n_solid) = 0.0_wp
    erosion_term(1:n_solid) = 0.0_wp

    IF ( qpj(1) .LE. 1.0e-5_wp ) RETURN
    
    r_h = qpj(1)
    r_u = qpj(n_vars+1)
    r_v = qpj(n_vars+2)
    r_alphas(1:n_solid) = qpj(5:4+n_solid)
    alphas_tot = SUM(r_alphas)

    IF ( erosion_coeff .GT. 0.0_wp ) THEN

       mod_vel = SQRT( r_u**2 + r_v**2 )
       ! empirical formulation (see Fagents & Baloga 2006, Eq. 5)
       ! here we use the solid volume fraction instead of relative density
       ! This term has units: m s-1    
       tot_solid_erosion = erosion_coeff * mod_vel * r_h * ( 1.0_wp-alphas_tot )&
            * ( 1.0_wp - erodible_porosity )
       
       erosion_term(1:n_solid) = erodible_fract(1:n_solid)  * tot_solid_erosion

    ELSE

       erosion_term(1:n_solid) = 0.0_wp
       
    END IF

    IF ( alphas_tot .LE. alphastot_min ) RETURN
    
    r_T = qpj(4)

    IF ( gas_flag ) THEN

       ! continuous phase is air
       r_rho_c = pres / ( sp_gas_const_a * r_T )

    ELSE

       ! continuous phase is liquid
       r_rho_c = rho_l

    END IF

    r_rho_m = ( 1.0_wp - alphas_tot ) * r_rho_c                                 &
         + DOT_PRODUCT( r_alphas , rho_s ) 

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
          
          expA = 1301.0_wp / ( 998.333_wp + 8.1855_wp * ( Tc - 20.0_wp )     &
               + 0.00585_wp * ( Tc - 20.0_wp )**2 ) - 1.30223_wp
          
          alpha1 = alpha1_coeff * 1.0E-3_wp * 10.0_wp**expA
          
       ELSE
          
          expB = ( 1.3272_wp * ( 20.0_wp - Tc ) - 0.001053_wp *              &
               ( Tc - 20.0_wp )**2 ) / ( Tc + 105.0_wp )
          
          alpha1 = alpha1_coeff * 1.002E-3_wp * 10.0_wp**expB 
          
       END IF
       
       ! Fluid dynamic viscosity [kg m-1 s-1]
       fluid_visc = alpha1 * EXP( beta1 * alphas_tot )
       ! Kinematic viscosity [m2 s-1]
       kin_visc = fluid_visc / r_rho_m
       rhoc = r_rho_m

    ELSE

       ! Viscosity read from input file [m2 s-1]
       kin_visc = kin_visc_c
       rhoc = r_rho_c
       
    END IF

    DO i_solid=1,n_solid

       IF ( ( r_alphas(i_solid) .GT. 0.0_wp ) .AND. ( settling_flag ) ) THEN

          settling_vel = settling_velocity( diam_s(i_solid) , rho_s(i_solid) ,  &
               rhoc , kin_visc )

          deposition_term(i_solid) = r_alphas(i_solid) * settling_vel

          IF ( rheology_model .NE. 4 ) THEN

             ! Michaels and Bolger (1962) sedimentation correction accounting 
             ! for hindered settling due to the presence of particles
             deposition_term(i_solid) = deposition_term(i_solid) *              &
               ( 1.0_wp - MIN( 1.0_wp , alphas_tot / alpha_max ) )**hind_exp

          END IF

          ! limit the deposition (cannot be remove more than particles present
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
    
    RETURN

  END SUBROUTINE eval_erosion_dep_term


  !******************************************************************************
  !> \brief Topography modification related term
  !
  !> This subroutine evaluates the deposition term.
  !> \date 2019/11/08
  !> \param[in]     qpj                   physical variables 
  !> \param[in]     deposition_avg_term   averaged deposition terms 
  !> \param[in]     erosion_avg_term      averaged deposition terms 
  !> \param[out]    eqns_term             source terms for cons equations
  !> \param[out]    deposit_term          deposition rates for solids
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !******************************************************************************

  SUBROUTINE eval_bulk_debulk_term( qpj, deposition_avg_term, erosion_avg_term, &
       eqns_term, topo_term )

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: qpj(n_vars+2)                !< physical variables 
    REAL(wp), INTENT(IN) :: deposition_avg_term(n_solid) !< deposition term
    REAL(wp), INTENT(IN) :: erosion_avg_term(n_solid)    !< erosion term

    REAL(wp), INTENT(OUT):: eqns_term(n_eqns)
    REAL(wp), INTENT(OUT):: topo_term

    REAL(wp) :: entr_coeff
    REAL(wp) :: air_entr
    REAL(wp) :: mag_vel 
    REAL(wp) :: mag_vel2 

    REAL(wp) :: r_h          !< real-value flow thickness
    REAL(wp) :: r_u          !< real-value x-velocity
    REAL(wp) :: r_v          !< real-value y-velocity
    REAL(wp) :: r_T          !< real-value temperature [K]
    REAL(wp) :: r_Ri         !< real-value Richardson number
    REAL(wp) :: r_rho_m      !< real-value mixture density [kg/m3]
    REAL(wp) :: r_rho_c      !< real-value carrier phase density [kg/m3]
    REAL(wp) :: r_red_grav   !< real-value reduced gravity

    REAL(wp) :: dep_tot
    REAL(wp) :: ers_tot
    REAL(wp) :: rho_dep_tot
    REAL(wp) :: rho_ers_tot

    REAL(wp) :: coeff_porosity
    

    IF ( qpj(1) .LE. 0.0_wp ) THEN

       eqns_term(1:n_eqns) = 0.0_wp
       topo_term = 0.0_wp

       RETURN

    END IF

    CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav)

    r_h = qpj(1)
    r_u = qpj(n_vars+1)
    r_v = qpj(n_vars+2)
    r_T = qpj(4)

    mag_vel2 = r_u**2 + r_v**2

    IF ( entrainment_flag .AND.  ( r_h .GT. 0.0_wp ) .AND.                      &
         ( r_Ri .GT. 0.0_wp ) ) THEN

       entr_coeff = 0.075_wp / SQRT( 1.0_wp + 718.0_wp * r_Ri**2.4_wp )

       air_entr = entr_coeff * SQRT(mag_vel2)

    ELSE

       air_entr = 0.0_wp

    END IF

    eqns_term(1:n_eqns) = 0.0_wp

    
    ! solid total volume deposition rate [m s-1]
    dep_tot = SUM( deposition_avg_term )
    ! solid total volume deposition rate [m s-1]
    ers_tot = SUM( erosion_avg_term )

    ! solid total mass deposition rate [kg m-2 s-1]
    rho_dep_tot = DOT_PRODUCT( rho_s , deposition_avg_term )
    ! solid total mass erosion rate [kg m-2 s-1]
    rho_ers_tot = DOT_PRODUCT( rho_s , erosion_avg_term )

    ! coefficient to compute (eroded/deposited) volume of continuous phase
    ! from volume of solid  
    coeff_porosity = erodible_porosity / ( 1.0_wp - erodible_porosity )
    
    ! total mass equation source term [kg m-2 s-1]:
    ! deposition, erosion and entrainment are considered
    eqns_term(1) = rho_a_amb * air_entr + rho_ers_tot - rho_dep_tot +           &
         coeff_porosity * ( rho_c_sub * ers_tot - r_rho_c * dep_tot )          
         
    ! x-momenutm equation source term [kg m-1 s-2]:
    ! only deposition contribute to change in momentum, erosion does not carry
    ! any momentum inside the flow
    eqns_term(2) = - r_u * ( rho_dep_tot + r_rho_c * coeff_porosity * dep_tot )

    ! y-momentum equation source term [kg m-1 s-2]:
    ! only deposition contribute to change in momentum, erosion does not carry
    ! any momentum inside the flow
    eqns_term(3) = - r_v * ( rho_dep_tot + r_rho_c * coeff_porosity * dep_tot )

    ! Temperature/Energy equation source term [kg s-3]:
    ! deposition, erosion and entrainment are considered
    IF ( energy_flag ) THEN
       
       eqns_term(4) = - r_T * ( SUM( rho_s * sp_heat_s * deposition_avg_term )  &
            + r_rho_c * sp_heat_c * coeff_porosity * dep_tot )                  &
            - 0.5_wp*mag_vel2 * rho_dep_tot                                     &
            + T_erodible * ( SUM( rho_s * sp_heat_s * erosion_avg_term )        &
            + rho_c_sub * sp_heat_c * ers_tot * coeff_porosity )                &
            + T_ambient * sp_heat_a * rho_a_amb * air_entr
       
    ELSE
       
       eqns_term(4) = - r_T * ( SUM( rho_s * sp_heat_s * deposition_avg_term )  &
            + r_rho_c * sp_heat_c * coeff_porosity * dep_tot )                  &
            + T_erodible * ( SUM( rho_s * sp_heat_s * erosion_avg_term )        &
            + rho_c_sub * sp_heat_c * ers_tot * coeff_porosity )                &
            + T_ambient * sp_heat_a * rho_a_amb * air_entr

    END IF

    ! solid phase mass equation source term [kg m-2 s-1]:
    ! due to solid erosion and deposition
    eqns_term(5:4+n_solid) = rho_s(1:n_solid) * ( erosion_avg_term(1:n_solid)   &
         - deposition_avg_term(1:n_solid) )
    
    ! erodible layer thickness source terms [m s-1]:
    ! due to erosion and deposition of solid+continuous phase
    topo_term = ( dep_tot - ers_tot ) / ( 1.0_wp - erodible_porosity ) 

    RETURN

  END SUBROUTINE eval_bulk_debulk_term

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
         alphal_source , time_param

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: time
    REAL(wp), INTENT(IN) :: vect_x
    REAL(wp), INTENT(IN) :: vect_y
    REAL(wp), INTENT(OUT) :: source_bdry(n_vars)

    REAL(wp) :: t_rem
    REAL(wp) :: t_coeff
    REAL(wp) :: pi_g

    IF ( time .GE. time_param(4) ) THEN

       ! The exponents of t_coeff are such that Ri does not depend on t_coeff
       source_bdry(1) = 0.0_wp
       source_bdry(2) = 0.0_wp
       source_bdry(3) = 0.0_wp
       source_bdry(4) = T_source
       source_bdry(5:4+n_solid) = alphas_source(1:n_solid)

       IF ( gas_flag .AND. liquid_flag ) THEN

          source_bdry(n_vars) = alphal_source

       END IF

       source_bdry(n_vars+1) = 0.0_wp
       source_bdry(n_vars+2) = 0.0_wp

       RETURN

    END IF

    t_rem = MOD( time , time_param(1) )

    pi_g = 4.0_wp * ATAN(1.0_wp) 

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
    source_bdry(5:4+n_solid) = alphas_source(1:n_solid)

    IF ( gas_flag .AND. liquid_flag ) THEN

       source_bdry(n_vars) = alphal_source

    END IF

    source_bdry(n_vars+1) = t_coeff**0.5_wp * vel_source * vect_x
    source_bdry(n_vars+2) = t_coeff**0.5_wp * vel_source * vect_y 

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
  !
  !> @author 
  !> Mattia de' Michieli Vitturi
  !
  !------------------------------------------------------------------------------

  REAL(wp) FUNCTION settling_velocity(diam,rhos,rhoc,kin_visc)

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: diam          !< particle diameter [m]
    REAL(wp), INTENT(IN) :: rhos          !< particle density [kg/m3]
    REAL(wp), INTENT(IN) :: rhoc          !< carrier phase density [kg/m3]
    REAL(wp), INTENT(IN) :: kin_visc      !< carrier phase viscosity

    REAL(wp) :: Rey           !< Reynolds number
    REAL(wp) :: inv_sqrt_C_D  !< Reciprocal of sqrt of Drag coefficient

    ! loop variables
    INTEGER :: i              !< loop counter for iterative procedure    
    REAL(wp) :: const_part    !< term not changing in iterative procedure
    REAL(wp) :: inv_sqrt_C_D_old  !< previous iteration sqrt of drag coefficient
    REAL(wp) :: set_vel_old   !< previous iteration settling velocity

    INTEGER :: dig          !< order of magnitude of settling velocity

    inv_sqrt_C_D = 1.0_wp

    const_part =  SQRT( 0.75_wp * ( rhos / rhoc - 1.0_wp ) * diam * grav )

    settling_velocity = const_part * inv_sqrt_C_D

    Rey = diam * settling_velocity / kin_visc

    IF ( Rey .LE. 1000.0_wp ) THEN

       C_D_loop:DO i=1,20

          set_vel_old = settling_velocity
          inv_sqrt_C_D_old = inv_sqrt_C_D
          inv_sqrt_C_D = SQRT( Rey / ( 24.0_wp * ( 1.0_wp +                     &
               0.15_wp*Rey**(0.687_wp) ) ) )

          settling_velocity = const_part * inv_sqrt_C_D

          IF ( ABS( set_vel_old - settling_velocity ) / set_vel_old             &
               .LT. 1.0E-5_wp ) THEN

             ! round to first three significative digits
             ! dig = FLOOR(LOG10(set_vel_old))
             ! settling_velocity = 10.0_wp**(dig-3)                             &
             !      * FLOOR( 10.0_wp**(-dig+3)*set_vel_old ) 

             RETURN

          END IF

          Rey = diam * settling_velocity / kin_visc

       END DO C_D_loop

    END IF

    RETURN

  END FUNCTION settling_velocity

END MODULE constitutive_2d


