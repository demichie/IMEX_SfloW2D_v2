!********************************************************************************
!> \brief Stochastic module
!
!> This moduel contains all the procedures for the stochastic variable
!
!> \date 11/06/2025
!> @author 
!> Zeno Geddo
!
!********************************************************************************
MODULE stochastic_module
  
  ! external variables
  USE parameters_2d, ONLY : wp , sp, rheology_model
  
  USE parameters_2d, ONLY : n_eqns , n_vars , n_solid , n_add_gas , n_quad ,    &
       n_stoch_vars , n_pore_vars
  
  USE solver_2d, ONLY : Z, conv_kernel
  USE solver_2d, ONLY: solve_cells, j_cent, k_cent
  USE solver_2d, ONLY: q , qp
  USE constitutive_2d, ONLY : T_ambient
  USE constitutive_2d, ONLY: qc_to_qp
  USE parameters_2d, ONLY : output_stoch_vars_flag, length_spatial_corr,        &
        stoch_transport_flag
  USE geometry_2d, ONLY : cell_size, comp_cells_x, comp_cells_y
  USE OMP_LIB
     
  ! variables related to OU process
  
  !LOGICAL :: sym_noise_flag ! use symmetric noise or not

  !> Noise parameter:\n
  !> - 0    => the noise will be simmetric
  !> 1      => abs val is taken
  !> -1     => -abs val is taken
  !> .
  REAL(wp) :: sym_noise

  REAL(wp) :: tau_stochastic

  REAL(wp) :: std_max, std_min

  REAL(wp) :: std_slope_factor ! Fr_0_stochastic
  
  REAL(wp) :: noise_pow_val !(|Z|^power)
  
  ! variables related to statistics
  REAL(wp) :: Z_min, Z_max, Z_mean, Z_std
  REAL(wp) :: percentiles(9) ! 5,10,20,30,50,70,80,90,95
  
CONTAINS
  
  REAL(wp) FUNCTION MeanFieldCorrection(g,h,Fr) ! working only in the case of mu(Fr)
    
    !> Compute the mean field correction 
    REAL(wp), INTENT(IN) :: g,h,Fr
    REAL(wp) :: factor, sFR, dsFR_dFr, spatial_corr
    
    ! Compute a normalization factor
    factor = (g**2._wp) / SQRT(g*h)
    
    ! Compute the intensity of the noise (std of the process)  
    sFR = std_max + ( std_min - std_max) * EXP(- Fr / std_slope_factor )
    
    ! Compute the derivative of the intensity of the noise
    dsFR_dFr = -(std_min - std_max) / std_slope_factor *                        &
         EXP(-Fr / std_slope_factor)
    
    ! Compute factor for spatial correlation
    spatial_corr = 2._wp * tau_stochastic ! ... will introduce correlation 
    
    ! Compute mean field correction
    MeanFieldCorrection = - factor * sFR * dsFR_dFr * spatial_corr

  END FUNCTION MeanFieldCorrection  


  SUBROUTINE getSteadyStateZ
  ! should find a better way to understand when the process is stable.
  ! Since sigma is varing in space and time, the process may never be stable.
    USE solver_2d, ONLY : comp_cells_x, comp_cells_y
    IMPLICIT none
    INTEGER :: j, k , l 
    INTEGER :: i, n_iter
    REAL(wp) :: dt
    REAL(wp) :: t_steady
    
    ! Initialize all stochastic process to zero
    Z = 0.0_wp 

    ! Define how many iteration to do for the burn in
    t_steady = 10_wp * tau_stochastic 
    dt = 0.1_wp * tau_stochastic
    n_iter = CEILING(t_steady/dt)
    ! WRITE(*,*) 'dt,n_inter',dt,n_iter

    ! Allocate and compute convolution kernel if needed (will be keept in memory)
    IF (length_spatial_corr .GT. cell_size) THEN
        CALL genConvolutionKernel()
    END IF
    
    ! Compute n iterations to get to a steady state OU process
    DO i=1,n_iter
       CALL update_stochastic_variable(dt)       
    END DO

    IF (stoch_transport_flag) THEN

       !$OMP PARALLEL DO private(j,k)

       DO l = 1,solve_cells

          j = j_cent(l)
          k = k_cent(l)
          
          qp(5+n_solid+n_add_gas,j,k) = Z(j,k)
          q(5+n_solid+n_add_gas,j,k) = q(1,j,k) * Z(j,k)
          !WRITE(*,*) j,k,qp(5+n_solid+n_add_gas,j,k),q(5+n_solid+n_add_gas,j,k)
          !READ(*,*)
          
       END DO
       
       !$OMP END PARALLEL DO

       WRITE(*,*) 'end stochastic initialization'
       
    END IF

    RETURN
    
  END SUBROUTINE getSteadyStateZ 


  SUBROUTINE genConvolutionKernel()  
      IMPLICIT none
      INTEGER :: n_nodes_per_dim, x_index, y_index
      REAL(wp) :: x_center, y_center, x, y

      ! Calculate grid dimensions (assuming square grid)
      n_nodes_per_dim = CEILING(length_spatial_corr / cell_size)

      ! Get center of the grid (assuming square grid)
      x_center = (cell_size * REAL(n_nodes_per_dim))  / 2._wp  
      y_center = x_center 

      ! Allocate the output array
      ALLOCATE(conv_kernel(n_nodes_per_dim, n_nodes_per_dim))

      ! Compute 2D Gaussian values (at centers of cells)
      DO y_index = 1, n_nodes_per_dim 
        DO x_index = 1, n_nodes_per_dim
          x =  cell_size * (x_index-1) + 0.5_wp * cell_size
          y =  cell_size * (y_index-1) + 0.5_wp * cell_size
          conv_kernel(x_index, y_index) =  evalGaussian2d(x, y, x_center, y_center)
        END DO
      END DO

      ! Renormalize values of the kernel
      conv_kernel = conv_kernel / SUM(conv_kernel)

  END SUBROUTINE genConvolutionKernel

  REAL(wp) FUNCTION evalGaussian2d(x, y, x_center, y_center)
      ! Sampling a 2D gaussian with centered at (x_center, y_center)
      ! The covariance matrix is [[s,0],[0,s]] (so it is simmetric)
      ! s = length_spatial_corr 
      IMPLICIT none
      REAL(wp), INTENT(IN) :: x, y, x_center, y_center
      REAL(wp) :: dx, dy, s, PI

      ! Calculate distances from the center
      dx = x - x_center
      dy = y - y_center

      ! Standard 2D Gaussian formula (with sigma=1 for standard Gaussian)
      s = length_spatial_corr
      PI = 4.0_wp * ACOS(-1.0_wp)
      evalGaussian2d = EXP(-0.5_wp * (((dx/s)**2._wp) + ((dy/s)**2._wp))) /           &
                        (2._wp * PI* s**2._wp) 

  END FUNCTION evalGaussian2d


  SUBROUTINE update_stochastic_variable(dt)
    ! UPDATE THE SOLUTION OF THE ORNSTEIN-UHLENBACK PROCESS USING EULER-MARUYAMA METHOD
    ! NOISE CAN BE TRANSPORTED, SOURCE TERM ADDED IN eval_mass_exchange_terms (IMPORTANT)
    IMPLICIT NONE
    INTEGER :: j,k
    INTEGER :: noise_size
    REAL(wp), INTENT(IN) :: dt
    REAL(wp) :: Fr
    REAL(wp) :: sigma_noise
    REAL(wp):: noise(comp_cells_x, comp_cells_y)
    REAL(wp):: conv_result(comp_cells_x, comp_cells_y) ! should swap idx?

    ! Generate standard gaussian noise over entire domain-> N(0,1)
    noise_size = comp_cells_x*comp_cells_y
    noise = reshape(GaussianNoise( noise_size ), [comp_cells_x, comp_cells_y])

    ! Convolve the gaussian noise to introduce spatial correlation if needed
    IF (length_spatial_corr .GT. cell_size) THEN
        ! (for convenience, the conv_kernel is generated only once before the burn in)
        CALL convolve_2d(noise, conv_kernel, conv_result) ! To fix (should swap indices?)!
        noise = conv_result
    END IF

    ! Loop over the entire grid to update stochastic process
    !$OMP PARALLEL
    !$OMP DO private(j,k)
    DO k = 1,comp_cells_y
       DO j = 1,comp_cells_x  
          ! Update the Ornstein-Uhlenback process
          sigma_noise = getSigmaNoise(j,k) 
          Z(j,k) = EulerMaruyamaScheme(Z(j,k), dt, sigma_noise, noise(j,k))
          ! apply non linear map to Z if needed (Z becomes asymmetric)
          IF (sym_noise .GT. 0.0_wp) THEN
            Z(j,k) = ABS(Z(j,k)) ** noise_pow_val ! generate only positive fluctuations 
          ELSEIF (sym_noise .LT. 0.0_wp) THEN
            Z(j,k) = -(ABS(Z(j,k)) ** noise_pow_val) ! generate only negative fluctuations
          END IF
       END DO
    END DO
    !$OMP END DO
    !$OMP END PARALLEL

    ! Compute statistic in space at given time if needed as outputs
    IF (output_stoch_vars_flag) THEN
        CALL OUBasicStatsInSpaceAtGivenTime(Z, Z_min, Z_mean, Z_max, Z_std)
        percentiles(:) = 0.0_wp ! Percentiles set to 0 to avoid the computations
        ! The percentiles shold not be computed at each iteration otherwise is too slow
        !CALL percentilesArrayAtGivenTime(Z, percentiles) ! it is very slow!!!!!!!!!!
    END IF
    
    RETURN
   
  END SUBROUTINE update_stochastic_variable


  FUNCTION GaussianNoise(noise_size)
    !> Generate a sample from a normal standard random distribution (mean=0, var=1)   
    !> Use Box-Muller transform : (https://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform)
    ! this is an uncorrelated noise: we will need to correlate it ?!!!!
    ! Correction to avoid leg(0.0): https://masuday.github.io/fortran_tutorial/random.html
    implicit none
    INTEGER :: noise_size
    REAL(wp) :: r(noise_size,2)
    REAL(wp) :: pi_g
    REAL(wp) :: GaussianNoise(noise_size)
    
    pi_g = 4.0_wp*ATAN(1.0_wp)
    CALL random_seed() ! = CALL random_seed(size=noise_size)  
    CALL random_number(r)
    GaussianNoise = sqrt ( - 2.0D+00 * log ( (1.0_wp-r(:,1)) ) ) * cos ( 2.0D+00 * pi_g * r(:,2) )
        
    RETURN
  END FUNCTION GaussianNoise

  
  REAL(wp) FUNCTION FroudeNumber(j,k)
    !> compute the froude number given the indices defining the location in the grid 
    USE parameters_2d, ONLY : n_solid , n_add_gas , n_stoch_vars , n_pore_vars
    USE constitutive_2d, ONLY: r_phys_var
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: j, k
    REAL(wp) :: Fr
    REAL(wp) :: r_h
    REAL(wp) :: r_u
    REAL(wp) :: r_v
    REAL(wp) :: r_alphas(n_solid)
    REAL(wp) :: r_rho_m
    REAL(wp) :: r_T
    REAL(wp) :: r_alphal
    REAL(wp) :: r_alphag(n_add_gas)
    REAL(wp) :: r_red_grav
    REAL(wp) :: R_ri
    REAL(wp) :: r_Zs(n_stoch_vars)
    REAL(wp) :: r_pore_pres(n_pore_vars)
    REAL(wp) :: p_dyn
    
    ! check that the thickness and velocity are > 0
    IF ( ( q(1,j,k) .GT. 0.0_wp ) .AND. ( ( q(2,j,k)**2 + q(3,j,k)**2 ) .GT.    &
         0.0_wp ) )  THEN

       CALL r_phys_var(q(:,j,k) , r_h , r_u , r_v , r_alphas , r_rho_m , r_T ,  &
            r_alphal , r_alphag , r_red_grav , p_dyn , r_Zs , r_pore_pres )

       Fr = ( r_u**2 + r_v**2 ) / SQRT( r_red_grav * r_h )

    ELSE

       ! set to Fr=0 if h and ||u|| are null
       Fr = 0.0_wp 

    END IF

    FroudeNumber = Fr

  END FUNCTION FroudeNumber

 REAL(wp) FUNCTION VelocityNorm(j,k)
  !> compute the norm of the velocity given the indices defining the location in the grid 
  USE parameters_2d, ONLY : n_solid , n_add_gas , n_stoch_vars , n_pore_vars
  USE constitutive_2d, ONLY: r_phys_var
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: j, k
  REAL(wp) :: r_h
  REAL(wp) :: r_u
  REAL(wp) :: r_v
  REAL(wp) :: r_alphas(n_solid)
  REAL(wp) :: r_rho_m
  REAL(wp) :: r_T
  REAL(wp) :: r_alphal
  REAL(wp) :: r_alphag(n_add_gas)
  REAL(wp) :: r_red_grav
  REAL(wp) :: R_ri
  REAL(wp) :: r_Zs(n_stoch_vars)
  REAL(wp) :: r_pore_pres(n_pore_vars)
  REAL(wp) :: p_dyn
  
  ! check that the thickness and velocity are > 0
  IF ( ( q(1,j,k) .GT. 0.0_wp ) .AND. ( ( q(2,j,k)**2 + q(3,j,k)**2 ) .GT.      &
       0.0_wp ) )  THEN

     CALL r_phys_var(q(:,j,k) , r_h , r_u , r_v , r_alphas , r_rho_m , r_T ,    &
          r_alphal , r_alphag , r_red_grav , p_dyn , r_Zs , r_pore_pres )

     VelocityNorm = r_u**2 + r_v**2 

  ELSE

     ! set to VelocityNorm=0 if h and ||u|| are null
     VelocityNorm = 0.0_wp 

  END IF

END FUNCTION VelocityNorm


  REAL(wp) FUNCTION getSigmaNoise(j,k)
  ! Compute the intensity of the noise depending of the friction used
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: j, k
    REAL(wp) :: Fr
    REAL(wp) :: U_norm
    IF ( rheology_model .EQ. 9 ) THEN   
      ! Fr = FroudeNumber(j,k)   ! OLD UNSTABLE
      ! getSigmaNoise = expFormNoise(Fr) ! OLD UNSTABLE
      U_norm = VelocityNorm(j,k)
      getSigmaNoise = expFormNoise(U_norm*U_norm) !
    ELSEIF (rheology_model .EQ. 10) THEN 
      U_norm = VelocityNorm(j,k)
      getSigmaNoise = expFormNoise(U_norm)
    ELSE
      getSigmaNoise = std_max
    END IF  
    RETURN
  END FUNCTION getSigmaNoise

  REAL(wp) FUNCTION expFormNoise(val)
    ! Compute the bounded intensity of the noise of the stochastic process
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: val
    expFormNoise = std_max + ( std_min - std_max) * exp(- val / std_slope_factor )
  END FUNCTION expFormNoise

  REAL(wp) FUNCTION EulerMaruyamaScheme(Zij, dt, sigma, noise_ij)
    ! Update the Ornstein-Uhlenback process using EULER-MARUYAMA Method
    ! Z(j,k) = Z(j,k) - (dt / tau_stochastic) * Z(j,k) + &
    ! sigma_noise * SQRT(2.0_wp * (dt / tau_stochastic)) * noise(j,k)
    IMPLICIT NONE
    REAL(wp), INTENT(IN) :: Zij, dt, sigma, noise_ij
    EulerMaruyamaScheme = Zij - (dt / tau_stochastic) * Zij +             &
            sigma * SQRT(2.0_wp * (dt / tau_stochastic)) * noise_ij
  END FUNCTION

subroutine convolve_2d(input_signal, kernel, result)
  ! Should work also if arrays have not the same size
  ! dimension expected: (y,x)
  implicit none
  real(wp), dimension(:,:), intent(in) :: input_signal, kernel
  real(wp), dimension(:,:), intent(out) :: result
  integer :: len_kernel_x, len_kernel_y, len_inp_sig_x,                 &
  len_inp_sig_y
  integer :: half_len_inp_sig_x, half_len_inp_sig_y,                    &
  half_len_k_x, half_len_k_y
  integer :: left_padding, right_padding, south_padding,                &
  north_padding
  integer :: inx_x, inx_y, idx_centre_x, idx_centre_y,                  &
  idx_k_x, idx_k_y, shift_x, shift_y
  real(wp), dimension(:,:), allocatable :: padded_signal
  real(wp) :: sum_at_position
  
  ! Get shape inp signal (should be moved outside)
  len_inp_sig_y = size(input_signal, 1)
  len_inp_sig_x = size(input_signal, 2)
  half_len_inp_sig_x = (len_inp_sig_x - 1) / 2
  half_len_inp_sig_y = (len_inp_sig_y - 1) / 2
  
  ! Get shape kernel (should be moved outside)
  len_kernel_y = size(kernel, 1)
  len_kernel_x = size(kernel, 2)
  half_len_k_x = (len_kernel_x - 1) / 2
  half_len_k_y = (len_kernel_y - 1) / 2

  ! Calculate the required padding for each dimension
  ! If padding is odd, one more zero is added to the right side.
  left_padding = half_len_k_x
  right_padding = (len_kernel_x - 1) - left_padding
  south_padding = half_len_k_y
  north_padding = (len_kernel_y - 1) - south_padding

  ! Allocate array for padded inp signal
  allocate(padded_signal(len_inp_sig_y + south_padding + north_padding, &
   len_inp_sig_x + left_padding + right_padding))
  
  ! Pad the input signal with zeros
  padded_signal = 0._wp
  padded_signal(south_padding + 1:south_padding + len_inp_sig_y,        &
  left_padding + 1:left_padding + len_inp_sig_x) = input_signal
  
  ! Initialize to zero the result
  result = 0._wp

  ! Loop over all inp signal cells
  do inx_x = 1, len_inp_sig_x
    do inx_y = 1, len_inp_sig_y
      idx_centre_x = left_padding + inx_x
      idx_centre_y = south_padding + inx_y
      sum_at_position = 0._wp

      ! Loop to convolve values around the current cell
      do idx_k_x = 1, len_kernel_x
        do idx_k_y = 1, len_kernel_y
          shift_x = -idx_k_x + half_len_k_x + 2 !check idx
          shift_y = -idx_k_y + half_len_k_y + 2 !check idx
          ! sum_at_position = 0.1_wp

          ! here there may be a bug in the indexing (segmentation fault-invalid memory reference) 
          if ((idx_centre_y + shift_y .LT. 1) .or.  (idx_centre_y + shift_y .GT. size(padded_signal,1))) THEN
            
            print*,"problem y idx"
            print*,shape(kernel)
            print*,shape(input_signal)
            print*,south_padding, north_padding
            print*,left_padding, right_padding
            print*,shape(padded_signal)
            print*, inx_x, inx_y, idx_k_x, idx_k_y
            print*,idx_centre_y + shift_y
            STOP
          END if
          if ((idx_centre_x + shift_x .LT. 1) .or.  (idx_centre_x + shift_x .GT. size(padded_signal,1))) THEN
            print*,"problem x idx"
            print*,shape(kernel)
            print*,shape(input_signal)
            print*,south_padding, north_padding
            print*,left_padding, right_padding
            print*,shape(padded_signal)
            print*, inx_x, inx_y, idx_k_x, idx_k_y
            print*,idx_centre_x + shift_x
            STOP  
          END if  

          sum_at_position = sum_at_position +                          &
          padded_signal(idx_centre_y + shift_y, idx_centre_x + shift_x)&
          * kernel(idx_k_y, idx_k_x)
        end do
      end do

      result(inx_y, inx_x) = sum_at_position
    end do
  end do

  deallocate(padded_signal)
end subroutine convolve_2d


subroutine OUBasicStatsInSpaceAtGivenTime(OUSolution, min_val, mean_val, max_val, std_dev)
  !> Computes statistics (in space) of a given array representing the solution (at a given time) of an Ornstein-Uhlenbeck process.
  !> Contains the entire domain, where many values may be zero !  
  !> This subroutine calculates the minimum, maximum, mean, standard deviation for the provided array.
  !>
  !> INPUT:
  !> - OUSolution: Array representing the time series of an Ornstein-Uhlenbeck process.
  !>
  !> OUTPUT:
  !> - min_val: Minimum value of the array.
  !> - max_val: Maximum value of the array.
  !> - mean_val: Mean value of the array.
  !> - std_dev: Standard deviation of the array.

  implicit none
  real(wp), dimension(:,:), intent(in) :: OUSolution
  real(wp), intent(out) :: min_val, mean_val, max_val, std_dev 
  real(wp) :: sumSol, sum_sq
  integer :: n

  ! Get array size
  n = size(OUSolution)

  ! Check for empty array
  if (n == 0) then
    write(*, *) "Error: Empty array passed to statsOUAtGivenTime"
    stop
  end if

  ! Calculate min and max
  min_val = minval(OUSolution)
  max_val = maxval(OUSolution)

  ! Calculate mean
  sumSol = sum(OUSolution)
  mean_val = sumSol / real(n)

  ! Calculate standard deviation
  sum_sq = sum((OUSolution - mean_val)**2)
  std_dev = sqrt(sum_sq / real(n))

  !WRITE(*,*) 'min_val',min_val
  !WRITE(*,*) 'max_val',max_val
  !WRITE(*,*) 'mean_val',mean_val
end subroutine OUBasicStatsInSpaceAtGivenTime


subroutine percentilesArrayAtGivenTime(Array2d, percentiles)
!> Compute the percentile of the OU solutions at actual time
!> Percentiles found using the quicksort algoritm
!> Should add a flag to conside only the cells where h>0
  implicit none
  real(wp), dimension(:, :), intent(in) :: Array2d
  real(wp), dimension(9), intent(out) :: percentiles
  real(wp), dimension(:), allocatable :: flatArray
  integer :: n

  ! Get total number of elements
  n = size(Array2d) 

  ! Check for empty array
  if (n == 0) then
    write(*, *) "Error: Empty array passed to statsOUAtGivenTime"
    stop
  end if

  ! Flatten the 2D array
  allocate(flatArray(n))
  flatArray = reshape(Array2d, shape(flatArray))

  ! Sort Flattened array
  call quickSort(flatArray, 1, n)

  ! Calculate percentiles
  percentiles(1) = flatArray(percentileIndex(n, 5._wp))
  percentiles(2) = flatArray(percentileIndex(n, 10._wp))
  percentiles(3) = flatArray(percentileIndex(n, 20._wp))
  percentiles(4) = flatArray(percentileIndex(n, 30._wp))
  percentiles(5) = flatArray(percentileIndex(n, 50._wp))
  percentiles(6) = flatArray(percentileIndex(n, 70._wp))
  percentiles(7) = flatArray(percentileIndex(n, 80._wp))
  percentiles(8) = flatArray(percentileIndex(n, 90._wp))
  percentiles(9) = flatArray(percentileIndex(n, 95._wp))

  ! Deallocate temporary array
  deallocate(flatArray)

end subroutine percentilesArrayAtGivenTime


subroutine quickSort(Array1d, start_idx, end_idx)
!> Sort the given array using the QuickSort algorithm
!> Run-time complexity : 
!> 1) Best case O(n log(n))
!> 2) Average case O(n log(n))
!> 3) Worst case O(n**2)
    implicit none
    real(wp), intent(in out) :: Array1d(:)
    integer, intent(in) :: start_idx, end_idx ! idx array
    integer :: pivot ! partition index

    ! Check if there are more than one element in the array
    if (start_idx < end_idx) then 
        pivot = partition(Array1d, start_idx, end_idx)
        ! Recursively apply quickSort to the left and right subarrays
        call quickSort(Array1d, start_idx, pivot - 1)
        call quickSort(Array1d, pivot + 1, end_idx)
    end if
end subroutine quickSort


integer function partition(Array1d, start_idx, end_idx)
!> Helper function of the quickSort algorithm
    implicit none
    real(wp), intent(in out) :: Array1d(:)
    integer, intent(in) :: start_idx, end_idx
    integer :: pi, pivot, i, temp, j

    ! Get the pivot element 
    pivot = Array1d(end_idx)
    ! Initialize the index of the smaller element
    i = start_idx - 1
    ! Loop through the array and rearrange elements based on the pivot value
    do j = start_idx, end_idx - 1
        if (Array1d(j) <= pivot) then
            ! Swap arr(i+1) and arr(j)
            i = i + 1
            temp = Array1d(i)
            Array1d(i) = Array1d(j)
            Array1d(j) = temp
        end if
    end do

    ! Swap arr(i+1) and arr(high) to place the pivot in its correct position
    i = i + 1
    temp = Array1d(i)
    Array1d(i) = Array1d(end_idx)
    Array1d(end_idx) = temp
    partition = i
end function partition


integer function percentileIndex(size_arr, p)
!> Find the index of the percentile p in array of size n
  integer, intent(in) :: size_arr
  real(wp), intent(in) :: p
  percentileIndex = max(1, min(size_arr, ceiling(real(size_arr) * p / 100._wp)))
end function percentileIndex
  

END MODULE stochastic_module
