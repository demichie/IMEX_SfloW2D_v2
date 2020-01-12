!********************************************************************************
!> \brief Constitutive equations
!********************************************************************************
MODULE constitutive_2d

  USE parameters_2d, ONLY : dp
  USE parameters_2d, ONLY : n_eqns , n_vars , n_solid
  USE parameters_2d, ONLY : rheology_flag , rheology_model , energy_flag ,      &
       liquid_flag , gas_flag

  IMPLICIT none

  !> flag used for size of implicit non linear-system
  LOGICAL, ALLOCATABLE :: implicit_flag(:)

  !> flag to activate air entrainment
  LOGICAL :: entrainment_flag
  
  !> gravitational acceleration 
  REAL(dp) :: grav

  !> drag coefficients (Voellmy-Salm model)
  REAL(dp) :: mu
  REAL(dp) :: xi

  !> drag coefficients (B&W model)
  REAL(dp) :: friction_factor

  !> drag coefficients (plastic model)
  REAL(dp) :: tau

  !> evironment temperature [K]
  REAL(dp) :: T_env

  !> radiative coefficient
  REAL(dp) :: rad_coeff

  !> friction coefficient
  REAL(dp) :: frict_coeff

  !> reference temperature [K]
  REAL(dp) :: T_ref

  !> reference kinematic viscosity [m2/s]
  REAL(dp) :: nu_ref

  !> viscosity parameter [K-1] (b in Table 1 Costa & Macedonio, 2005)
  REAL(dp) :: visc_par

  !> velocity boundary layer fraction of total thickness
  REAL(dp) :: emme

  !> specific heat [J kg-1 K-1]
  REAL(dp) :: c_p

  !> atmospheric heat trasnfer coefficient [W m-2 K-1] (lambda in C&M, 2005)
  REAL(dp) :: atm_heat_transf_coeff

  !> fractional area of the exposed inner core (f in C&M, 2005)
  REAL(dp) :: exp_area_fract

  !> Stephan-Boltzmann constant [W m-2 K-4]
  REAL(dp), PARAMETER :: SBconst = 5.67D-8

  !> emissivity (eps in Costa & Macedonio, 2005)
  REAL(dp) :: emissivity

  !> thermal boundary layer fraction of total thickness
  REAL(dp) :: enne

  !> temperature of lava-ground interface
  REAL(dp) :: T_ground

  !> thermal conductivity [W m-1 K-1] (k in Costa & Macedonio, 2005)
  REAL(dp) :: thermal_conductivity

  !--- START Lahars rheology model parameters

  !> 1st param for yield strenght empirical relationship (O'Brian et al, 1993)
  REAL(dp) :: alpha2    ! (units: kg m-1 s-2)

  !> 2nd param for yield strenght empirical relationship (O'Brian et al, 1993)
  REAL(dp) :: beta2     ! (units: nondimensional) 

  !> ratio between reference value from input and computed values from eq.
  REAL(dp) :: alpha1_coeff ! (units: nondimensional )

  !> 2nd param for fluid viscosity empirical relationship (O'Brian et al, 1993)
  REAL(dp) :: beta1     ! (units: nondimensional, input parameter)

  !> Empirical resistance parameter (dimensionless, input parameter)
  REAL(dp) :: Kappa

  !> Mannings roughness coefficient ( units: T L^(-1/3) )
  REAL(dp) :: n_td

  !--- END Lahars rheology model parameters

  !> Specific heat of carrier phase (gas or liquid)
  REAL(dp) :: sp_heat_c  ! ( initialized from input)   

  !> Ambient density of air ( units: kg m-3 )
  REAL(dp) :: rho_a_amb

  !> Specific heat of air (units: J K-1 kg-1)
  REAL(dp) :: sp_heat_a

  !> Specific gas constant of air (units: J kg-1 K-1)
  REAL(dp) :: sp_gas_const_a

  !> Kinematic viscosity of air (units: m2 s-1)
  REAL(dp) :: kin_visc_a

  !> Kinematic viscosity of liquid (units: m2 s-1)
  REAL(dp) :: kin_visc_l

  !> Kinematic viscosity of carrier phase (units: m2 s-1)
  REAL(dp) :: kin_visc_c

  !> Temperature of ambient air (units: K)
  REAL(dp) :: T_ambient

  !> Density of sediments ( units: kg m-3 )
  REAL(dp), ALLOCATABLE :: rho_s(:)

  !> Diameter of sediments ( units: m )
  REAL(dp), ALLOCATABLE :: diam_s(:)

  !> Specific heat of solids (units: J K-1 kg-1)
  REAL(dp), ALLOCATABLE :: sp_heat_s(:)

  !> Flag to determine if sedimentation is active
  LOGICAL :: settling_flag

  !> Hindered settling velocity (units: m s-1 )
  REAL(dp) :: settling_vel

  !> erosion model coefficient  (units: m-1 )
  REAL(dp), ALLOCATABLE :: erosion_coeff(:)

  !> temperature of solid substrate (units: K)
  REAL(dp) :: T_s_substrate

  !> ambient pressure (units: Pa)
  REAL(dp) :: pres

  !> liquid density (units: kg m-3)
  REAL(dp) :: rho_l

  !> Sepcific heat of liquid (units: J K-1 kg-1)
  REAL(dp) :: sp_heat_l

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

    USE parameters_2d, ONLY : eps_sing
    IMPLICIT none

    REAL(dp), INTENT(IN) :: r_qj(n_vars)       !< real-value conservative var
    REAL(dp), INTENT(OUT) :: r_h               !< real-value flow thickness
    REAL(dp), INTENT(OUT) :: r_u               !< real-value x-velocity
    REAL(dp), INTENT(OUT) :: r_v               !< real-value y-velocity
    REAL(dp), INTENT(OUT) :: r_alphas(n_solid) !< real-value solid volume fracts
    REAL(dp), INTENT(OUT) :: r_rho_m           !< real-value mixture density
    REAL(dp), INTENT(OUT) :: r_T               !< real-value temperature
    REAL(dp), INTENT(OUT) :: r_alphal          !< real-value liquid volume fract

    REAL(dp) :: r_inv_rhom
    REAL(dp) :: r_xs(n_solid)     !< real-value solid mass fraction
    REAL(dp) :: r_xs_tot

    REAL(dp) :: r_Ri            !< real-value Richardson number
    REAL(dp) :: r_xl            !< real-value liquid mass fraction
    REAL(dp) :: r_xc            !< real-value carrier phase mass fraction
    REAL(dp) :: r_alphac        !< real-value carrier phase volume fraction
    REAL(dp) :: r_rho_c         !< real-value carrier phase density [kg/m3]
    REAL(dp) :: r_red_grav      !< real-value reduced gravity
    REAL(dp) :: r_sp_heat_mix   !< Specific heat of mixture

    ! compute solid mass fractions
    IF ( r_qj(1) .GT. 1.D-25 ) THEN

       r_xs(1:n_solid) = r_qj(5:4+n_solid) / r_qj(1)

    ELSE

       r_xs(1:n_solid) = 0.0_dp

    END IF

    r_xs_tot = SUM(r_xs)

    IF ( gas_flag .AND. liquid_flag ) THEN

       ! compute liquid mass fraction
       IF ( r_qj(1) .GT. 1.D-25 ) THEN

          r_xl = r_qj(n_vars) / r_qj(1)

       ELSE

          r_xl = 0.0_dp

       END IF

       ! compute carrier phase (gas) mass fraction
       r_xc =  1.0_dp - r_xs_tot - r_xl

       ! specific heat of the mixutre: mass average of sp. heat pf phases
       r_sp_heat_mix = SUM( r_xs(1:n_solid) * sp_heat_s(1:n_solid) )            &
            + r_xl * sp_heat_l + r_xc * sp_heat_c

    ELSE

       ! compute carrier phase (gas or liquid) mass fraction
       r_xc = 1.0_dp - r_xs_tot

       ! specific heaf of the mixutre: mass average of sp. heat pf phases
       r_sp_heat_mix = SUM( r_xs(1:n_solid) * sp_heat_s(1:n_solid) )            &
            + r_xc * sp_heat_c

    END IF

    ! compute temperature from energy
    IF ( r_qj(1) .GT. 1.D-25 ) THEN

       IF ( energy_flag ) THEN

          r_T = ( r_qj(4) - ( r_qj(2)**2 + r_qj(3)**2 ) / ( 2.0_dp*r_qj(1) ) ) /&
               ( r_qj(1) * r_sp_heat_mix ) 

       ELSE

          r_T = r_qj(4) / ( r_qj(1) * r_sp_heat_mix ) 

       END IF

       IF ( r_T .LE. 0.0_dp ) r_T = T_ambient

    ELSE

       r_T = T_ambient

    END IF

    IF ( gas_flag ) THEN

       ! carrier phase is gas
       r_rho_c =  pres / ( sp_gas_const_a * r_T )
       sp_heat_c = sp_heat_a

    ELSE

       ! carrier phase is liquid
       r_rho_c = rho_l
       sp_heat_c = sp_heat_l

    END IF

    IF ( gas_flag .AND. liquid_flag ) THEN

       r_inv_rhom = ( SUM(r_xs(1:n_solid) / rho_s(1:n_solid)) + r_xl / rho_l    &
            + r_xc / r_rho_c )

       r_rho_m = 1.0_dp / r_inv_rhom

       r_alphal = r_xl * r_rho_m / rho_l

    ELSE

       r_inv_rhom = ( SUM(r_xs(1:n_solid) / rho_s(1:n_solid)) + r_xc / r_rho_c )

       r_rho_m = 1.0_dp / r_inv_rhom

    END IF

    ! convert from mass fraction to volume fraction
    r_alphas(1:n_solid) = r_xs(1:n_solid) * r_rho_m / rho_s(1:n_solid)

    ! convert from mass fraction to volume fraction
    r_alphac = r_xc * r_rho_m / r_rho_c

    r_h = r_qj(1) / r_rho_m

    ! reduced gravity
    r_red_grav = ( r_rho_m - rho_a_amb ) / r_rho_m * grav

    ! velocity components
    IF ( r_qj(1) .GT. eps_sing ) THEN

       r_u = r_qj(2) / r_qj(1)
       r_v = r_qj(3) / r_qj(1)

    ELSE

       r_u = SQRT(2.0_dp) * r_qj(1) * r_qj(2) / SQRT( r_qj(1)**4 + eps_sing**4 )
       r_v = SQRT(2.0_dp) * r_qj(1) * r_qj(3) / SQRT( r_qj(1)**4 + eps_sing**4 )

    END IF

    ! Richardson number
    IF ( ( r_u**2 + r_v**2 ) .GT. 0.0_dp ) THEN

       r_Ri = r_red_grav * r_h / ( r_u**2 + r_v**2 )

    ELSE

       r_Ri = 0.0_dp

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

  SUBROUTINE c_phys_var(c_qj,h,u,v,T,rho_m,red_grav,alphas)

    USE COMPLEXIFY
    USE parameters_2d, ONLY : eps_sing
    IMPLICIT none

    COMPLEX(dp), INTENT(IN) :: c_qj(n_vars)
    COMPLEX(dp), INTENT(OUT) :: h               !< height [m]
    COMPLEX(dp), INTENT(OUT) :: u               !< velocity (x direction) [m s-1]
    COMPLEX(dp), INTENT(OUT) :: v               !< velocity (y direction) [m s-1]
    COMPLEX(dp), INTENT(OUT) :: T               !< temperature [K]
    COMPLEX(dp), INTENT(OUT) :: rho_m           !< mixture density [kg m-3]
    COMPLEX(dp), INTENT(OUT) :: red_grav        !< reduced gravity
    COMPLEX(dp), INTENT(OUT) :: alphas(n_solid) !< sediment volume fractions

    COMPLEX(dp) :: inv_rhom
    COMPLEX(dp) :: xs(n_solid)             !< sediment mass fractions
    COMPLEX(dp) :: xs_tot                  !< sum of solid mass fraction
    COMPLEX(dp) :: Ri                      !< Richardson number
    COMPLEX(dp) :: xl                      !< liquid mass fraction
    COMPLEX(dp) :: xc                      !< carrier phase mass fraction
    COMPLEX(dp) :: alphal                  !< liquid volume fraction
    COMPLEX(dp) :: alphac                  !< carrier phase volume fraction
    COMPLEX(dp) :: sp_heat_mix             !< Specific heat of mixture
    COMPLEX(dp) :: rho_c     !< Density of carrier phase in the mixture 


    ! compute solid mass fractions
    IF ( REAL(c_qj(1)) .GT. 1.D-25 ) THEN

       xs(1:n_solid) = c_qj(5:4+n_solid) / c_qj(1)

    ELSE

       xs(1:n_solid) = CMPLX(0.0_dp,0.0_dp,dp)

    END IF

    xs_tot = SUM(xs)

    IF ( gas_flag .AND. liquid_flag ) THEN

       ! compute liquid mass fraction
       IF ( REAL(c_qj(1)) .GT. 1.D-25 ) THEN

          xl = c_qj(n_vars) / c_qj(1)

       ELSE

          xl = CMPLX(0.0_dp,0.0_dp,dp) 

       END IF

       ! compute carrier phase (gas) mass fraction
       xc =  CMPLX(1.0_dp,0.0_dp,dp) - xs_tot - xl

       ! specific heaf of the mixutre: mass average of sp. heat pf phases
       sp_heat_mix = SUM( xs(1:n_solid) * sp_heat_s(1:n_solid) )                &
            + xl * sp_heat_l + xc * sp_heat_c

    ELSE

       ! compute carrier phase (gas or liquid) mass fraction
       xc = CMPLX(1.0_dp,0.0_dp,dp) - xs_tot

       ! specific heaf of the mixutre: mass average of sp. heat pf phases
       sp_heat_mix = SUM( xs(1:n_solid) * sp_heat_s(1:n_solid) ) + xc * sp_heat_c

    END IF

    ! compute temperature from energy
    IF ( REAL(c_qj(1)) .GT. 1.D-25 ) THEN

       IF ( energy_flag ) THEN

          T = ( c_qj(4) - ( c_qj(2)**2 + c_qj(3)**2 ) / ( 2.0_dp*c_qj(1) ) ) /  &
               ( c_qj(1) * sp_heat_mix ) 

       ELSE

          T = c_qj(4) / ( c_qj(1) * sp_heat_mix ) 

       END IF

       IF ( REAL(T) .LE. 0.0_dp ) T = CMPLX(T_ambient,0.0_dp,dp)

    ELSE

       T = CMPLX(T_ambient,0.0_dp,dp)

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

       inv_rhom = ( SUM(xs(1:n_solid) / rho_s(1:n_solid)) + xl / rho_l          &
            + xc / rho_c )

       rho_m = 1.0_dp / inv_rhom

       alphal = xl * rho_m / rho_l

    ELSE

       inv_rhom = ( SUM(xs(1:n_solid) / rho_s(1:n_solid)) + xc / rho_c )

       rho_m = 1.0_dp / inv_rhom

    END IF

    ! convert from mass fraction to volume fraction
    alphas(1:n_solid) = xs(1:n_solid) * rho_m / rho_s(1:n_solid)

    ! convert from mass fraction to volume fraction
    alphac = xc * rho_m / rho_c

    h = c_qj(1) / rho_m

    ! reduced gravity
    red_grav = ( rho_m - rho_a_amb ) / rho_m * grav

    ! velocity components
    IF ( REAL( c_qj(1) ) .GT. eps_sing ) THEN

       u = c_qj(2) / c_qj(1)
       v = c_qj(3) / c_qj(1)

    ELSE

       u = SQRT(2.0_dp) * c_qj(1) * c_qj(2) / SQRT( c_qj(1)**4 + eps_sing**4 )
       v = SQRT(2.0_dp) * c_qj(1) * c_qj(3) / SQRT( c_qj(1)**4 + eps_sing**4 )

    END IF

    ! Richardson number
    IF ( REAL( u**2 + v**2 ) .GT. 0.0_dp ) THEN

       Ri = red_grav * h / ( u**2 + v**2 )

    ELSE

       Ri = CMPLX(0.0_dp,0.0_dp,dp)

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

    REAL(dp), INTENT(IN) :: qpj(n_vars+2) !< real-value physical variables
    REAL(dp), INTENT(OUT) :: r_Ri         !< real-value Richardson number
    REAL(dp), INTENT(OUT) :: r_rho_m      !< real-value mixture density [kg/m3]
    REAL(dp), INTENT(OUT) :: r_rho_c !< real-value carrier phase density [kg/m3]
    REAL(dp), INTENT(OUT) :: r_red_grav   !< real-value reduced gravity
    REAL(dp) :: r_u                       !< real-value x-velocity
    REAL(dp) :: r_v                       !< real-value y-velocity
    REAL(dp) :: r_h                       !< real-value flow thickness
    REAL(dp) :: r_alphas(n_solid)         !< real-value solid volume fractions
    REAL(dp) :: r_T                       !< real-value temperature [K]
    REAL(dp) :: r_alphal                  !< real-value liquid volume fraction

    r_h = qpj(1)

    IF ( qpj(1) .LE. 0.0_dp ) THEN

       r_u = 0.0_dp
       r_v = 0.0_dp
       r_T = T_ambient
       r_alphas(1:n_solid) = 0.0_dp
       r_red_grav = 0.0_dp
       r_Ri = 0.0_dp
       r_rho_m = rho_a_amb
       IF ( gas_flag .AND. liquid_flag ) r_alphal = 0.0_dp

       RETURN

    END IF

    r_u = qpj(n_vars+1)
    r_v = qpj(n_vars+2)
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
       r_rho_m = ( 1.0_dp - SUM(r_alphas) - r_alphal ) * r_rho_c                &
            + SUM( r_alphas * rho_s ) + r_alphal * rho_l

    ELSE

       ! density of mixture of carrier phase and solids
       r_rho_m = ( 1.0_dp - SUM(r_alphas) ) * r_rho_c + SUM( r_alphas * rho_s ) 

    END IF

    ! reduced gravity
    r_red_grav = ( r_rho_m - rho_a_amb ) / r_rho_m * grav

    ! Richardson number
    IF ( ( r_u**2 + r_v**2 ) .GT. 0.0_dp ) THEN

       r_Ri = r_red_grav * r_h / ( r_u**2 + r_v**2 )

    ELSE

       r_Ri = 0.0_dp

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

    REAL(dp), INTENT(IN) :: qc(n_vars)
    REAL(dp), INTENT(OUT) :: qp(n_vars+2)

    REAL(dp) :: r_h               !< real-value flow thickness
    REAL(dp) :: r_u               !< real-value x-velocity
    REAL(dp) :: r_v               !< real-value y-velocity
    REAL(dp) :: r_alphas(n_solid) !< real-value solid volume fractions
    REAL(dp) :: r_rho_m           !< real-value mixture density [kg/m3]
    REAL(dp) :: r_T               !< real-value temperature [K]
    REAL(dp) :: r_alphal          !< real-value liquid volume fraction

    
    IF ( qc(1) .GT. 0.0_dp ) THEN

       CALL r_phys_var(qc,r_h,r_u,r_v,r_alphas,r_rho_m,r_T,r_alphal)

       qp(1) = r_h

       qp(2) = r_h*r_u
       qp(3) = r_h*r_v

       qp(4) = r_T
       qp(5:4+n_solid) = r_alphas(1:n_solid)

       IF ( gas_flag .AND. liquid_flag ) qp(n_vars) = r_alphal

       qp(n_vars+1) = r_u
       qp(n_vars+2) = r_v

    ELSE

       qp(1) = 0.0_dp           ! h
       qp(2) = 0.0_dp           ! hu
       qp(3) = 0.0_dp           ! hv
       qp(4) = T_ambient      ! T
       qp(5:n_vars) = 0.0_dp    ! alphas
       qp(n_vars+1) = 0.0_dp    ! u
       qp(n_vars+2) = 0.0_dp    ! v

    END IF

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

    REAL(dp), INTENT(IN) :: qp(n_vars+2)
    REAL(dp), INTENT(OUT) :: qc(n_vars)

    REAL(dp) :: r_sp_heat_mix
    REAL(dp) :: sum_sl

    REAL(dp) :: r_u               !< real-value x-velocity
    REAL(dp) :: r_v               !< real-value y-velocity
    REAL(dp) :: r_h               !< real-value flow thickness
    REAL(dp) :: r_hu              !< real-value volumetric x-flow
    REAL(dp) :: r_hv              !< real-value volumetric y-flow
    REAL(dp) :: r_alphas(n_solid) !< real-value solid volume fractions
    REAL(dp) :: r_xl              !< real-value liquid mass fraction
    REAL(dp) :: r_xc              !< real-value carrier phase mass fraction
    REAL(dp) :: r_T               !< real-value temperature [K]
    REAL(dp) :: r_alphal          !< real-value liquid volume fraction
    REAL(dp) :: r_alphac          !< real-value carrier phase volume fraction
    REAL(dp) :: r_rho_m           !< real-value mixture density [kg/m3]
    REAL(dp) :: r_rho_c           !< real-value carrier phase density [kg/m3]
    REAL(dp) :: r_xs(n_solid)     !< real-value solid mass fraction

    r_h = qp(1)

    IF ( r_h .GT. 0.0_dp ) THEN

       r_hu = qp(2)
       r_hv = qp(3)

       r_u = qp(n_vars+1)
       r_v = qp(n_vars+2)

    ELSE

       r_hu = 0.0_dp
       r_hv = 0.0_dp

       r_u = 0.0_dp
       r_v = 0.0_dp

       qc(1:n_vars) = 0.0_dp
       RETURN

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
       IF ( ( SUM(r_alphas) + r_alphal ) .GT. 1.0_dp ) THEN

          sum_sl = SUM(r_alphas) + r_alphal
          r_alphas(1:n_solid) = r_alphas(1:n_solid) / sum_sl
          r_alphal = r_alphal / sum_sl

       ELSEIF ( ( SUM(r_alphas) + r_alphal ) .LT. 0.0_dp ) THEN

          r_alphas(1:n_solid) = 0.0_dp
          r_alphal = 0.0_dp

       END IF

       ! carrier phase volume fraction
       r_alphac = 1.0_dp - SUM(r_alphas) - r_alphal

       ! volume averaged mixture density: carrier (gas) + solids + liquid
       r_rho_m = r_alphac * r_rho_c + SUM( r_alphas * rho_s ) + r_alphal * rho_l

       ! liquid mass fraction
       r_xl = r_alphal * rho_l / r_rho_m

       ! solid mass fractions
       r_xs(1:n_solid) = r_alphas(1:n_solid) * rho_s(1:n_solid) / r_rho_m

       ! carrier (gas) mass fraction
       r_xc = r_alphac * r_rho_c / r_rho_m

       ! mass averaged mixture specific heat
       r_sp_heat_mix =  SUM( r_xs*sp_heat_s ) + r_xl*sp_heat_l + r_xc*sp_heat_c

    ELSE

       ! mixture of carrier phase ( gas or liquid ) and solid

       ! check and corrections on dispersed phases
       IF ( SUM(r_alphas) .GT. 1.0_dp ) THEN

          r_alphas(1:n_solid) = r_alphas(1:n_solid) / SUM(r_alphas)

       ELSEIF ( SUM(r_alphas).LT. 0.0_dp ) THEN

          r_alphas(1:n_solid) = 0.0_dp

       END IF

       ! carrier (gas or liquid) volume fraction
       r_alphac = 1.0_dp - SUM(r_alphas) 

       ! volume averaged mixture density: carrier (gas or liquid) + solids
       r_rho_m = r_alphac * r_rho_c + SUM( r_alphas * rho_s ) 

       ! solid mass fractions
       r_xs(1:n_solid) = r_alphas(1:n_solid) * rho_s(1:n_solid) / r_rho_m

       ! carrier (gas or liquid) mass fraction
       r_xc = r_alphac * r_rho_c / r_rho_m

       ! mass averaged mixture specific heat
       r_sp_heat_mix =  SUM( r_xs * sp_heat_s ) + r_xc * sp_heat_c

    END IF

    qc(1) = r_rho_m * r_h 
    qc(2) = r_rho_m * r_hu
    qc(3) = r_rho_m * r_hv

    IF ( energy_flag ) THEN

       IF ( r_h .GT. 0.0_dp ) THEN

          ! total energy (internal and kinetic)
          qc(4) = r_h * r_rho_m * ( r_sp_heat_mix * r_T                         &
               + 0.5_dp * ( r_hu**2 + r_hv**2 ) / r_h**2 )

       ELSE

          qc(4) = 0.0_dp

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

    REAL(dp), INTENT(IN) :: qpj(n_vars)
    REAL(dp), INTENT(IN) :: Bj
    REAL(dp), INTENT(OUT) :: qp2j(3)

    qp2j(1) = qpj(1) + Bj

    IF ( qpj(1) .LE. 0.0_dp ) THEN

       qp2j(2) = 0.0_dp
       qp2j(3) = 0.0_dp

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

    REAL(dp), INTENT(IN) :: qpj(n_vars+2)

    REAL(dp), INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)

    REAL(dp) :: r_h          !< real-value flow thickness [m]
    REAL(dp) :: r_u          !< real-value x-velocity [m s-1]
    REAL(dp) :: r_v          !< real-value y-velocity [m s-1]
    REAL(dp) :: r_Ri         !< real-value Richardson number
    REAL(dp) :: r_rho_m      !< real-value mixture density [kg m-3]
    REAL(dp) :: r_rho_c      !< real-value carrier phase density [kg m-3]
    REAL(dp) :: r_red_grav   !< real-value reduced gravity [m s-2]

    CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav)

    r_h = qpj(1)
    r_u = qpj(n_vars+1)
    r_v = qpj(n_vars+2)

    IF ( r_red_grav * r_h .LT. 0.0_dp ) THEN

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

    REAL(dp), INTENT(IN)  :: qpj(n_vars+2)
    REAL(dp), INTENT(OUT) :: vel_min(n_vars) , vel_max(n_vars)

    REAL(dp) :: r_h          !< real-value flow thickness
    REAL(dp) :: r_u          !< real-value x-velocity
    REAL(dp) :: r_v          !< real-value y-velocity
    REAL(dp) :: r_Ri         !< real-value Richardson number
    REAL(dp) :: r_rho_m      !< real-value mixture density [kg/m3]
    REAL(dp) :: r_rho_c      !< real-value carrier phase density [kg/m3]
    REAL(dp) :: r_red_grav   !< real-value reduced gravity

    CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav)

    r_h = qpj(1)
    r_u = qpj(n_vars+1)
    r_v = qpj(n_vars+2)

    IF ( r_red_grav * r_h .LT. 0.0_dp ) THEN

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

    IMPLICIT none

    REAL(dp), INTENT(IN) :: qcj(n_vars)
    REAL(dp), INTENT(IN) :: qpj(n_vars+2)
    INTEGER, INTENT(IN) :: dir

    REAL(dp), INTENT(OUT) :: flux(n_eqns)

    REAL(dp) :: r_h          !< real-value flow thickness [m]
    REAL(dp) :: r_u          !< real-value x-velocity [m s-1]
    REAL(dp) :: r_v          !< real-value y-velocity [m s-1]
    REAL(dp) :: r_Ri         !< real-value Richardson number
    REAL(dp) :: r_rho_m      !< real-value mixture density [kg m-3]
    REAL(dp) :: r_rho_c      !< real-value carrier phase density [kg m-3]
    REAL(dp) :: r_red_grav   !< real-value reduced gravity [m s-2]

    pos_thick:IF ( qcj(1) .GT. 0.0_dp ) THEN

       r_h = qpj(1)
       r_u = qpj(n_vars+1)
       r_v = qpj(n_vars+2)

       CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav)

       IF ( dir .EQ. 1 ) THEN

          ! Mass flux in x-direction: u * ( rhom * h )
          flux(1) = r_u * qcj(1)

          ! x-momentum flux in x-direction + hydrostatic pressure term
          flux(2) = r_u * qcj(2) + 0.5_dp * r_rho_m * r_red_grav * r_h**2

          ! y-momentum flux in x-direction: u * ( rho * h * v )
          flux(3) = r_u * qcj(3)

          IF ( energy_flag ) THEN

             ! ENERGY flux in x-direction
             flux(4) = r_u * ( qcj(4) + 0.5_dp * r_rho_m * r_red_grav * r_h**2 )

          ELSE

             ! Temperature flux in x-direction: u * ( h * T )
             flux(4) = r_u * qcj(4)

          END IF

          ! Mass flux of solid in x-direction: u * ( h * alphas * rhos )
          flux(5:4+n_solid) = r_u * qcj(5:4+n_solid)

          ! Solid flux can't be larger than total flux
          IF ( ( flux(1) .GT. 0.0_dp ) .AND. ( SUM(flux(5:4+n_solid)) / flux(1) &
               .GT. 1.0_dp ) ) THEN

             flux(5:4+n_solid) = &
                  flux(5:4+n_solid) / SUM(flux(5:4+n_solid)) * flux(1)

          END IF

          ! Mass flux of liquid in x-direction: u * ( h * alphal * rhol )
          IF ( gas_flag .AND. liquid_flag ) flux(n_vars) = r_u * qcj(n_vars)

       ELSEIF ( dir .EQ. 2 ) THEN

          ! flux G (derivated wrt y in the equations)
          flux(1) = r_v * qcj(1)

          flux(2) = r_v * qcj(2)

          flux(3) = r_v * qcj(3) + 0.5_dp * r_rho_m * r_red_grav * r_h**2

          IF ( energy_flag ) THEN

             ! ENERGY flux in x-direction
             flux(4) = r_v * ( qcj(4) + 0.5_dp * r_rho_m * r_red_grav * r_h**2 )

          ELSE

             ! Temperature flux in y-direction
             flux(4) = r_v * qcj(4)

          END IF

          ! Mass flux of solid in y-direction: v * ( h * alphas * rhos )
          flux(5:4+n_solid) = r_v * qcj(5:4+n_solid)

          ! Solid flux can't be larger than total flux
          IF ( ( flux(1) .GT. 0.0_dp ) .AND. ( SUM(flux(5:4+n_solid)) / flux(1) &
               .GT. 1.0_dp ) ) THEN

             flux(5:4+n_solid) = &
                  flux(5:4+n_solid) / SUM(flux(5:4+n_solid)) * flux(1)

          END IF

          ! Mass flux of liquid in x-direction: u * ( h * alphal * rhol )
          IF ( gas_flag .AND. liquid_flag ) flux(n_vars) = r_v * qcj(n_vars)

       END IF

    ELSE

       flux(1:n_eqns) = 0.0_dp

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

    COMPLEX(dp), INTENT(IN), OPTIONAL :: c_qj(n_vars)
    COMPLEX(dp), INTENT(OUT), OPTIONAL :: c_nh_term_impl(n_eqns)
    REAL(dp), INTENT(IN), OPTIONAL :: r_qj(n_vars)
    REAL(dp), INTENT(OUT), OPTIONAL :: r_nh_term_impl(n_eqns)

    COMPLEX(dp) :: h                       !< height [m]
    COMPLEX(dp) :: u                       !< velocity (x direction) [m/s]
    COMPLEX(dp) :: v                       !< velocity (y direction) [m/s]
    COMPLEX(dp) :: T                       !< temperature [K]
    COMPLEX(dp) :: rho_m                   !< mixture density [kg/m3]
    COMPLEX(dp) :: red_grav                !< reduced gravity
    COMPLEX(dp) :: alphas(n_solid)         !< sediment volume fractions
 
    COMPLEX(dp) :: qj(n_vars)
    COMPLEX(dp) :: nh_term(n_eqns)
    COMPLEX(dp) :: forces_term(n_eqns)

    COMPLEX(dp) :: mod_vel
    COMPLEX(dp) :: gamma
    REAL(dp) :: h_threshold

    INTEGER :: i

    !--- Lahars rheology model variables

    !> Temperature in C
    COMPLEX(dp) :: Tc

    COMPLEX(dp) :: expA , expB

    !> 1st param for fluid viscosity empirical relationship (O'Brian et al, 1993)
    COMPLEX(dp) :: alpha1    ! (units: kg m-1 s-1 )
    
    !> Fluid dynamic viscosity (units: kg m-1 s-1 )
    COMPLEX(dp) :: fluid_visc

    !> Total friction slope (dimensionless): s_f = s_v+s_td+s_y
    COMPLEX(dp) :: s_f

    !> Viscous slope component of total Friction (dimensionless)
    COMPLEX(dp) :: s_v

    !> Turbulent dispersive slope component of total friction (dimensionless)
    COMPLEX(dp) :: s_td

    IF ( ( thermal_conductivity .GT. 0.0_dp ) .OR. ( emme .GT. 0.0_dp ) ) THEN

       h_threshold = 1.D-10

    ELSE

       h_threshold = 0.0_dp

    END IF

    IF ( present(c_qj) .AND. present(c_nh_term_impl) ) THEN

       qj = c_qj

    ELSEIF ( present(r_qj) .AND. present(r_nh_term_impl) ) THEN

       DO i = 1,n_vars

          qj(i) = CMPLX( r_qj(i),0.0_dp,dp )

       END DO

    ELSE

       WRITE(*,*) 'Constitutive, eval_fluxes: problem with arguments'
       STOP

    END IF

    ! initialize and evaluate the forces terms
    forces_term(1:n_eqns) = CMPLX(0.0_dp,0.0_dp,dp)

    IF (rheology_flag) THEN

       CALL c_phys_var(qj,h,u,v,T,rho_m,red_grav,alphas)

       mod_vel = SQRT( u**2 + v**2 )

       ! Voellmy Salm rheology
       IF ( rheology_model .EQ. 1 ) THEN

          IF ( REAL(mod_vel) .NE. 0.0_dp ) THEN 

             ! IMPORTANT: grav3_surf is always negative 
             forces_term(2) = forces_term(2) - rho_m * ( u / mod_vel ) *        &
                  ( grav / xi ) * mod_vel ** 2

             forces_term(3) = forces_term(3) - rho_m * ( v / mod_vel ) *        &
                  ( grav / xi ) * mod_vel ** 2

          ENDIF

          ! Plastic rheology
       ELSEIF ( rheology_model .EQ. 2 ) THEN

          IF ( REAL(mod_vel) .NE. 0.0_dp ) THEN 

             forces_term(2) = forces_term(2) - rho_m * tau * (u/mod_vel)

             forces_term(3) = forces_term(3) - rho_m * tau * (v/mod_vel)

          ENDIF

          ! Temperature dependent rheology
       ELSEIF ( rheology_model .EQ. 3 ) THEN

          IF ( REAL(h) .GT. h_threshold ) THEN

             ! Equation 6 from Costa & Macedonio, 2005
             gamma = 3.0_dp * nu_ref / h * EXP( - visc_par * ( T - T_ref ) )

          ELSE

             ! Equation 6 from Costa & Macedonio, 2005
             gamma = 3.0_dp * nu_ref / h_threshold * EXP( - visc_par            &
                  * ( T - T_ref ) )

          END IF

          IF ( REAL(mod_vel) .NE. 0.0_dp ) THEN 

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
          Tc = T - 273.15_dp

          ! the dependance of viscosity on temperature is modeled with the
          ! equation presented at:
          ! https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781118131473.app3
          !
          ! In addition, we use a reference value provided in input at a 
          ! reference temperature. This value is used to scale the equation
          IF ( REAL(Tc) .LT. 20.0_dp ) THEN

             expA = 1301.0_dp / ( 998.333_dp + 8.1855_dp * ( Tc - 20.0_dp )     &
                  + 0.00585_dp * ( Tc - 20.0_dp )**2 ) - 1.30223_dp

             alpha1 = alpha1_coeff * 1.D-3 * 10.0_dp**expA

          ELSE

             expB = ( 1.3272_dp * ( 20.0_dp - Tc ) - 0.001053_dp *              &
                  ( Tc - 20.0_dp )**2 ) / ( Tc + 105.0_dp )

             alpha1 = alpha1_coeff * 1.002D-3 * 10.0_dp**expB 

          END IF

          ! Fluid viscosity 
          fluid_visc = alpha1 * EXP( beta1 * SUM(alphas) )

          IF ( REAL(h) .GT. h_threshold ) THEN

             ! Viscous slope component (dimensionless)
             s_v = Kappa * fluid_visc * mod_vel / ( 8.0_dp * rho_m * grav *h**2 )

             ! Turbulent dispersive component (dimensionless)
             s_td = n_td**2 * mod_vel**2 / ( h**(4.0_dp/3.0_dp) )

          ELSE

             ! Viscous slope component (dimensionless)
             s_v = Kappa * fluid_visc * mod_vel / ( 8.0_dp * rho_m * grav *     &
                  h_threshold**2 )

             ! Turbulent dispersive components (dimensionless)
             s_td = n_td**2 * (mod_vel**2) / ( h_threshold**(4.0_dp/3.0_dp) )

          END IF

          ! Total implicit friction slope (dimensionless)
          s_f = s_v + s_td

          IF ( REAL(mod_vel) .GT. 0.0_dp ) THEN

             ! same units of dqc(2)/dt: kg m-1 s-2
             forces_term(2) = forces_term(2) - grav * rho_m * h *               &
                  ( u / mod_vel ) * s_f

             ! same units of dqc(3)/dt: kg m-1 s-2
             forces_term(3) = forces_term(3) - grav * rho_m * h *               &
                  ( v / mod_vel ) * s_f

          END IF

       ELSEIF ( rheology_model .EQ. 5 ) THEN

          tau = 1.D-3 / ( 1.0_dp + 10.0_dp * h ) * mod_vel

          IF ( REAL(mod_vel) .NE. 0.0_dp ) THEN

             forces_term(2) = forces_term(2) - rho_m * tau * ( u / mod_vel )
             forces_term(3) = forces_term(3) - rho_m * tau * ( v / mod_vel )

          END IF


       ELSEIF ( rheology_model .EQ. 6 ) THEN

          IF ( REAL(mod_vel) .NE. 0.0_dp ) THEN 

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

    REAL(dp), INTENT(IN) :: grav3_surf

    REAL(dp), INTENT(IN) :: qcj(n_vars)
    REAL(dp), INTENT(OUT) :: nh_semi_impl_term(n_eqns)

    REAL(dp) :: forces_term(n_eqns)

    REAL(dp) :: mod_vel

    REAL(dp) :: h_threshold

    !--- Lahars rheology model variables

    !> Yield strenght (units: kg m-1 s-2)
    REAL(dp) :: tau_y

    !> Yield slope component of total friction (dimensionless)
    REAL(dp) :: s_y

    REAL(dp) :: r_h               !< real-value flow thickness
    REAL(dp) :: r_u               !< real-value x-velocity
    REAL(dp) :: r_v               !< real-value y-velocity
    REAL(dp) :: r_alphas(n_solid) !< real-value solid volume fractions
    REAL(dp) :: r_rho_m           !< real-value mixture density [kg/m3]
    REAL(dp) :: r_T               !< real-value temperature [K]
    REAL(dp) :: r_alphal          !< real-value liquid volume fraction


    ! initialize and evaluate the forces terms
    forces_term(1:n_eqns) = 0.0_dp

    IF (rheology_flag) THEN

       CALL r_phys_var(qcj,r_h,r_u,r_v,r_alphas,r_rho_m,r_T,r_alphal)

       mod_vel = SQRT( r_u**2 + r_v**2 )

       ! Voellmy Salm rheology
       IF ( rheology_model .EQ. 1 ) THEN

          IF ( mod_vel .GT. 0.0_dp ) THEN

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

          h_threshold = 1.D-20

          ! Yield strength (units: kg m-1 s-2)
          tau_y = alpha2 * EXP( beta2 * SUM(r_alphas) )

          IF ( r_h .GT. h_threshold ) THEN

             ! Yield slope component (dimensionless)
             s_y = tau_y / ( grav * r_rho_m * r_h )

          ELSE

             ! Yield slope component (dimensionless)
             s_y = tau_y / ( grav * r_rho_m * h_threshold )

          END IF

          IF ( mod_vel .GT. 0.0_dp ) THEN

             ! units of dqc(2)/dt (kg m-1 s-2)
             forces_term(2) = forces_term(2) - grav * r_rho_m * r_h *           &
                  ( r_u / mod_vel ) * s_y

             ! units of dqc(3)/dt (kg m-1 s-2)
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

  SUBROUTINE eval_expl_terms( Bprimej_x, Bprimej_y, source_xy, qpj, expl_term )

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: Bprimej_x
    REAL(dp), INTENT(IN) :: Bprimej_y
    REAL(dp), INTENT(IN) :: source_xy

    REAL(dp), INTENT(IN) :: qpj(n_vars+2)      !< local physical variables 
    REAL(dp), INTENT(OUT) :: expl_term(n_eqns) !< local explicit forces 

    REAL(dp) :: r_h          !< real-value flow thickness
    REAL(dp) :: r_u          !< real-value x-velocity
    REAL(dp) :: r_v          !< real-value y-velocity
    REAL(dp) :: r_Ri         !< real-value Richardson number
    REAL(dp) :: r_rho_m      !< real-value mixture density [kg/m3]
    REAL(dp) :: r_rho_c      !< real-value carrier phase density [kg/m3]
    REAL(dp) :: r_red_grav   !< real-value reduced gravity

    expl_term(1:n_eqns) = 0.0_dp

    IF ( qpj(1) .LE. 0.0_dp ) RETURN

    r_h = qpj(1)
    r_u = qpj(n_vars+1)
    r_v = qpj(n_vars+2)

    CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav)

    expl_term(2) = r_red_grav * r_rho_m * r_h * Bprimej_x

    expl_term(3) = r_red_grav * r_rho_m * r_h * Bprimej_y

    IF ( energy_flag ) THEN

       expl_term(4) = r_red_grav * r_rho_m * r_h * ( r_u * Bprimej_x            &
            + r_v * Bprimej_y )  

    ELSE

       expl_term(4) = 0.0_dp

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

    REAL(dp), INTENT(IN) :: qpj(n_vars+2)              !< physical variables 
    REAL(dp), INTENT(IN) :: dt

    REAL(dp), INTENT(OUT) :: erosion_term(n_solid)     !< erosion term
    REAL(dp), INTENT(OUT) :: deposition_term(n_solid)  !< deposition term

    REAL(dp) :: mod_vel

    REAL(dp) :: hind_exp
    REAL(dp) :: alpha_max

    INTEGER :: i_solid

    REAL(dp) :: r_h          !< real-value flow thickness
    REAL(dp) :: r_u          !< real-value x-velocity
    REAL(dp) :: r_v          !< real-value y-velocity
    REAL(dp) :: r_alphas(n_solid) !< real-value solid volume fractions
    REAL(dp) :: r_Ri         !< real-value Richardson number
    REAL(dp) :: r_rho_m      !< real-value mixture density [kg/m3]
    REAL(dp) :: r_rho_c      !< real-value carrier phase density [kg/m3]
    REAL(dp) :: r_red_grav   !< real-value reduced gravity


    ! parameters for Michaels and Bolger (1962) sedimentation correction
    alpha_max = 0.6_dp
    hind_exp = 4.65_dp

    deposition_term(1:n_solid) = 0.0_dp
    erosion_term(1:n_solid) = 0.0_dp

    IF ( qpj(1) .LE. 0.0_dp ) RETURN
    
    r_h = qpj(1)
    r_u = qpj(n_vars+1)
    r_v = qpj(n_vars+2)
    r_alphas(1:n_solid) = qpj(5:4+n_solid)

    CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav)

    DO i_solid=1,n_solid

       IF ( ( r_alphas(i_solid) .GT. 0.0_dp ) .AND. ( settling_flag ) ) THEN

          settling_vel = settling_velocity( diam_s(i_solid) , rho_s(i_solid) ,  &
               r_rho_c )

          deposition_term(i_solid) = r_alphas(i_solid) * settling_vel *         &
               ( 1.0_dp - MIN( 1.0_dp , SUM( r_alphas ) / alpha_max ) )**hind_exp

          deposition_term(i_solid) = MIN( deposition_term(i_solid) ,            &
               r_h * r_alphas(i_solid) / dt )

          IF ( deposition_term(i_solid) .LT. 0.0_dp ) THEN

             WRITE(*,*) 'eval_erosion_dep_term'
             WRITE(*,*) 'deposition_term(i_solid)',deposition_term(i_solid)
             READ(*,*)

          END IF

       END IF

       mod_vel = SQRT( r_u**2 + r_v**2 )

       IF ( r_h .GT. 1.D-2) THEN

          ! empirical formulation (see Fagents & Baloga 2006, Eq. 5)
          ! here we use the solid volume fraction instead of relative density
          ! This term has units: m s-1
          erosion_term(i_solid) = erosion_coeff(i_solid) * mod_vel * r_h        &
               * ( 1.0_dp - SUM(r_alphas) )

       ELSE

          erosion_term(i_solid) = 0.0_dp

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

  SUBROUTINE eval_topo_term( qpj , deposition_avg_term , erosion_avg_term ,      &
       eqns_term, deposit_term )

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: qpj(n_vars+2)                !< physical variables 
    REAL(dp), INTENT(IN) :: deposition_avg_term(n_solid) !< deposition term
    REAL(dp), INTENT(IN) :: erosion_avg_term(n_solid)    !< erosion term

    REAL(dp), INTENT(OUT):: eqns_term(n_eqns)
    REAL(dp), INTENT(OUT):: deposit_term(n_solid)

    REAL(dp) :: entr_coeff
    REAL(dp) :: air_entr
    REAL(dp) :: mag_vel 

    REAL(dp) :: r_h          !< real-value flow thickness
    REAL(dp) :: r_u          !< real-value x-velocity
    REAL(dp) :: r_v          !< real-value y-velocity
    REAL(dp) :: r_T          !< real-value temperature [K]
    REAL(dp) :: r_Ri         !< real-value Richardson number
    REAL(dp) :: r_rho_m      !< real-value mixture density [kg/m3]
    REAL(dp) :: r_rho_c      !< real-value carrier phase density [kg/m3]
    REAL(dp) :: r_red_grav   !< real-value reduced gravity


    IF ( qpj(1) .LE. 0.0_dp ) THEN

       eqns_term(1:n_eqns) = 0.0_dp
       deposit_term(1:n_solid) = 0.0_dp

       RETURN

    END IF

    CALL mixt_var(qpj,r_Ri,r_rho_m,r_rho_c,r_red_grav)

    r_h = qpj(1)
    r_u = qpj(n_vars+1)
    r_v = qpj(n_vars+2)
    r_T = qpj(4)

    mag_vel = SQRT( r_u**2 + r_v**2 ) 

    IF ( entrainment_flag .AND. ( mag_vel**2 .GT. 0.0_dp ) .AND.                &
         ( r_h .GT. 0.0_dp ) ) THEN

       entr_coeff = 0.075_dp / SQRT( 1.0_dp + 718.0_dp * MAX(0.0_dp,r_Ri)**2.4 )

       air_entr = entr_coeff * mag_vel

    ELSE

       air_entr = 0.0_dp

    END IF

    eqns_term(1:n_eqns) = 0.0_dp

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
            - 0.5_dp * mag_vel**2 * SUM( rho_s * deposition_avg_term )          &
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

    REAL(dp), INTENT(IN) :: time
    REAL(dp), INTENT(IN) :: vect_x
    REAL(dp), INTENT(IN) :: vect_y
    REAL(dp), INTENT(OUT) :: source_bdry(n_vars)

    REAL(dp) :: t_rem
    REAL(dp) :: t_coeff
    REAL(dp) :: pi_g

    IF ( time .GE. time_param(4) ) THEN

       ! The exponents of t_coeff are such that Ri does not depend on t_coeff
       source_bdry(1) = 0.0_dp
       source_bdry(2) = 0.0_dp
       source_bdry(3) = 0.0_dp
       source_bdry(4) = T_source
       source_bdry(5:4+n_solid) = alphas_source(1:n_solid)

       IF ( gas_flag .AND. liquid_flag ) THEN

          source_bdry(n_vars) = alphal_source

       END IF

       source_bdry(n_vars+1) = 0.0_dp
       source_bdry(n_vars+2) = 0.0_dp

       RETURN

    END IF

    t_rem = MOD( time , time_param(1) )

    pi_g = 4.0_dp * ATAN(1.0_dp) 

    t_coeff = 0.0_dp

    IF ( time_param(3) .EQ. 0.0_dp ) THEN

       IF ( t_rem .LE. time_param(2) ) t_coeff = 1.0_dp

    ELSE

       IF ( t_rem .LT. time_param(3) ) THEN

          t_coeff = 0.5_dp * ( 1.0_dp - COS( pi_g * t_rem / time_param(3) ) )

       ELSEIF ( t_rem .LE. ( time_param(2) - time_param(3) ) ) THEN

          t_coeff = 1.0_dp

       ELSEIF ( t_rem .LE. time_param(2) ) THEN

          t_coeff = 0.5_dp * ( 1.0_dp + COS( pi_g * ( ( t_rem - time_param(2) ) &
               / time_param(3) + 1.0_dp ) ) )

       END IF

    END IF

    ! The exponents of t_coeff are such that Ri does not depend on t_coeff
    source_bdry(1) = t_coeff * h_source
    source_bdry(2) = t_coeff**1.5_dp * h_source * vel_source * vect_x
    source_bdry(3) = t_coeff**1.5_dp * h_source * vel_source * vect_y
    source_bdry(4) = T_source
    source_bdry(5:4+n_solid) = alphas_source(1:n_solid)

    IF ( gas_flag .AND. liquid_flag ) THEN

       source_bdry(n_vars) = alphal_source

    END IF

    source_bdry(n_vars+1) = t_coeff**0.5_dp * vel_source * vect_x
    source_bdry(n_vars+2) = t_coeff**0.5_dp * vel_source * vect_y 

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

  REAL(dp) FUNCTION settling_velocity(diam,rhos,rhoc)

    IMPLICIT NONE

    REAL(dp), INTENT(IN) :: diam          !< particle diameter [m]
    REAL(dp), INTENT(IN) :: rhos          !< particle density [kg/m3]
    REAL(dp), INTENT(IN) :: rhoc          !< carrier phase density [kg/m3]

    REAL(dp) :: Rey           !< Reynolds number
    REAL(dp) :: C_D           !< Drag coefficient

    ! loop variables
    INTEGER :: i              !< loop counter for iterative procedure    
    REAL(dp) :: const_part    !< term not changing in iterative procedure
    REAL(dp) :: C_D_old       !< previous iteration drag coefficient
    REAL(dp) :: set_vel_old   !< previous iteration settling velocity

    INTEGER :: dig          !< order of magnitude of settling velocity

    C_D = 1.0_dp

    const_part =  SQRT( 0.75_dp * ( rhos / rhoc - 1.0_dp ) * diam * grav )

    settling_velocity = const_part / SQRT( C_D )

    Rey = diam * settling_velocity / kin_visc_c

    IF ( Rey .LE. 1000.0_dp ) THEN

       C_D_loop:DO i=1,20

          set_vel_old = settling_velocity
          C_D_old = C_D
          C_D = 24.0_dp / Rey * ( 1.0_dp + 0.15_dp * Rey**(0.687_dp) )

          settling_velocity = const_part / SQRT( C_D )

          IF ( ABS( set_vel_old - settling_velocity ) / set_vel_old             &
               .LT. 1.D-6 ) THEN

             ! round to first three significative digits
             dig = FLOOR(LOG10(set_vel_old))
             settling_velocity = 10.0_dp**(dig-3)                               &
                  * FLOOR( 10.0_dp**(-dig+3)*set_vel_old ) 

             EXIT C_D_loop

          END IF

          Rey = diam * settling_velocity / kin_visc_c

       END DO C_D_loop

    END IF

    RETURN

  END FUNCTION settling_velocity

END MODULE constitutive_2d


