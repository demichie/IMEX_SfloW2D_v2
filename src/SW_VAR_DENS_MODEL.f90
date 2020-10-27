!********************************************************************************
!> \mainpage SW_VAR_DENS_MODEL - Shallow Water Finite volume solver
!> SW_VAR_DENS_MODEL is a FORTRAN90 code designed to solve an hyperbolic 
!> system of partial differential equations with relaxation and source
!> terms.\n 
!> The model is discretized in time with an explicit-implicit Runge-Kutta
!> method where the hyperbolic part is solved explicetely and the other
!> terms (relaxation and surce) are treated implicitely.\n
!> The finite volume solver for the hyperbolic part of the system is based
!> on a semidiscrete central scheme and it is not tied on the specific 
!> eigenstructure of the model.\n
!> The implicit part is solved with a Newton-Raphson method where the 
!> elements of the Jacobian of the nonlinear system are evaluated 
!> numerically with a complex step derivative technique.\n
!> Version 1.0:\n
!> \n 
!> Github project page: http://demichie.github.io/SW_VAR_DENS_MODEL/
!> \n
!> \authors Mattia de' Michieli Vitturi (*)
!> (*) Istituto Nazionale di Geofisica e vulcanologia, sezione di Pisa\n
!>     Via Cesare Battisti, 53\n
!>     I-56125 Pisa, Italy \n
!>     E-mail: mattia.demichielivitturi@ingv.it \n
!********************************************************************************

!> \brief Main Program 
PROGRAM SW_VAR_DENS_MODEL

  USE constitutive_2d, ONLY : init_problem_param


  USE geometry_2d, ONLY : init_grid
  USE geometry_2d, ONLY : init_source
  USE geometry_2d, ONLY : topography_reconstruction

  USE geometry_2d, ONLY : dx,dy,B_cent
  ! USE geometry_2d, ONLY : comp_cells_x,comp_cells_y

  USE init_2d, ONLY : collapsing_volume

  USE inpout_2d, ONLY : init_param
  USE inpout_2d, ONLY : read_param
  USE inpout_2d, ONLY : update_param
  USE inpout_2d, ONLY : output_solution
  USE inpout_2d, ONLY : output_runout
  USE inpout_2d, ONLY : read_solution
  USE inpout_2d, ONLY : close_units
  USE inpout_2d, ONLY : n_probes
  USE inpout_2d, ONLY : output_probes

  USE inpout_2d, ONLY : output_runout_flag
  USE inpout_2d, ONLY : output_cons_flag
  USE inpout_2d, ONLY : output_esri_flag
  USE inpout_2d, ONLY : output_phys_flag

  USE solver_2d, ONLY : allocate_solver_variables
  USE solver_2d, ONLY : deallocate_solver_variables
  USE solver_2d, ONLY : imex_RK_solver
  USE solver_2d, ONLY : update_erosion_deposition_cell
  USE solver_2d, ONLY : timestep
  USE solver_2d, ONLY : check_solve

  USE inpout_2d, ONLY : restart

  USE parameters_2d, ONLY : wp

  USE parameters_2d, ONLY : t_start
  USE parameters_2d, ONLY : t_end
  USE parameters_2d, ONLY : t_output
  USE parameters_2d, ONLY : t_runout
  USE parameters_2d, ONLY : t_probes
  USE parameters_2d, ONLY : t_steady
  USE parameters_2d, ONLY : dt0
  USE parameters_2d, ONLY : topo_change_flag
  USE parameters_2d, ONLY : verbose_level
  USE parameters_2d, ONLY : n_solid
  USE parameters_2d, ONLY : n_vars
  USE parameters_2d, ONLY : radial_source_flag
  USE parameters_2d, ONLY : collapsing_volume_flag

  USE parameters_2d, ONLY : n_thickness_levels , n_dyn_pres_levels ,          &
       thickness_levels , dyn_pres_levels

  USE solver_2d, ONLY : q , qp , t, dt
  USE solver_2d, ONLY : hmax
  USE solver_2d, ONLY : thck_table ,  pdyn_table , vuln_table

  USE constitutive_2d, ONLY : qc_to_qp

  USE solver_2d, ONLY : solve_mask , solve_cells
  USE solver_2d, ONLY : j_cent , k_cent

  USE OMP_LIB

  IMPLICIT NONE

  REAL(wp) :: t1 , t2 , t3

  REAL(wp) :: rate
  INTEGER :: st1 , st2 , st3 , cr , cm

  REAL(wp) :: dt_old , dt_old_old
  LOGICAL :: stop_flag
  LOGICAL :: stop_flag_old

  INTEGER :: j,k,l

  INTEGER :: i_pdyn_lev , i_thk_lev , i_table

  INTEGER n_threads

  LOGICAL :: use_openmp = .false.

  !> Dynamic pressure
  REAL(wp) :: p_dyn
  


  WRITE(*,*) '---------------------'
  WRITE(*,*) 'SW_VAR_DENS_MODEL 1.0'
  WRITE(*,*) '---------------------'

  !$ use_openmp = .true.
  !$ print *, "OpenMP program"

  IF ( .NOT. use_openmp) THEN

     PRINT *, "Non-OpenMP simulation"

  ELSE

     !$ n_threads = omp_get_max_threads()
     !$ CALL OMP_SET_NUM_THREADS(n_threads)
     IF ( verbose_level .GE. 0 ) WRITE(*,*) 'Number of threads used',n_threads

  END IF



  ! First initialize the system_clock
  CALL system_clock(count_rate=cr)
  CALL system_clock(count_max=cm)
  rate = DBLE(cr)

  CALL cpu_time(t1)
  CALL system_clock (st1)

  CALL init_param

  CALL read_param

  CALL init_grid

  CALL init_problem_param

  CALL allocate_solver_variables

  IF ( restart ) THEN

     CALL read_solution

  ELSE

     IF ( collapsing_volume_flag ) THEN

        CALL collapsing_volume

     END IF

  END IF

  IF ( radial_source_flag ) CALL init_source

  t = t_start

  CALL check_solve(.TRUE.)

  IF ( topo_change_flag ) CALL topography_reconstruction

  IF ( verbose_level .GE. 0 ) THEN

     WRITE(*,*) 
     WRITE(*,*) '******** START COMPUTATION *********'
     WRITE(*,*)

  END IF

  IF ( verbose_level .GE. 1 ) THEN

     WRITE(*,*) 'Min q(1,:,:)=',MINVAL(q(1,:,:))
     WRITE(*,*) 'Max q(1,:,:)=',MAXVAL(q(1,:,:))

     WRITE(*,*) 'Min B(:,:)=',MINVAL(B_cent(:,:))
     WRITE(*,*) 'Max B(:,:)=',MAXVAL(B_cent(:,:))


     WRITE(*,*) 'size B_cent',size(B_cent,1),size(B_cent,2)

     WRITE(*,*) 'SUM(q(1,:,:)=',SUM(q(1,:,:))
     WRITE(*,*) 'SUM(B_cent(:,:)=',SUM(B_cent(:,:))

  END IF


  dt_old = dt0
  dt_old_old = dt_old
  t_steady = t_end
  stop_flag = .FALSE.

  vuln_table = .FALSE.

  !$OMP PARALLEL DO private(j,k,p_dyn,i_table,i_thk_lev,i_pdyn_lev)

  DO l = 1,solve_cells

     j = j_cent(l)
     k = k_cent(l)

     IF ( q(1,j,k) .GT. 0.0_wp ) THEN

        CALL qc_to_qp(q(1:n_vars,j,k) , qp(1:n_vars,j,k) , p_dyn )
        hmax(j,k) = qp(1,j,k)

        i_table = 0
        
        DO i_thk_lev=1,n_thickness_levels

           thck_table(j,k) = ( qp(1,j,k) .GE. thickness_levels(i_thk_lev) )

           DO i_pdyn_lev=1,n_dyn_pres_levels

              pdyn_table(j,k) = ( p_dyn .GE. dyn_pres_levels(i_pdyn_lev) ) 

              i_table = i_table + 1

              vuln_table(i_table,j,k) = ( thck_table(j,k) .AND. pdyn_table(j,k) )

           END DO

        END DO

     ELSE

        qp(1:n_vars,j,k) = 0.0_wp
        hmax(j,k) = 0.0_wp

     END IF

  END DO

  !$OMP END PARALLEL DO

  IF ( output_runout_flag ) CALL output_runout(t,stop_flag)

  IF ( output_cons_flag .OR. output_esri_flag .OR. output_phys_flag )           &
       CALL output_solution(t)

  IF ( n_probes .GT. 0 ) CALL output_probes(t)

  IF ( SUM(q(1,:,:)) .EQ. 0.0_wp ) t_steady = t_end

  IF ( verbose_level .GE. 0 ) THEN

     WRITE(*,FMT="(A3,F10.4,A5,F9.5,A9,ES11.3E3,A11,ES11.3E3,A9,ES11.3E3,A15,   &
          &ES11.3E3)")                                                          &
          't =',t,'dt =',dt0,                                                   &
          ' mass = ',dx*dy*SUM(q(1,:,:)) ,                                      &
          ' volume = ',dx*dy*SUM(qp(1,:,:)) ,                                   &
          ' area = ',dx*dy*COUNT(q(1,:,:).GT.1.D-5) ,                           &
          ' solid mass = ',dx*dy*SUM(q(5:4+n_solid,:,:))

  END IF

  CALL cpu_time(t2)
  CALL system_clock (st2)
 
  DO WHILE ( ( t .LT. t_end ) .AND. ( t .LT. t_steady ) )

     CALL update_param

     IF ( t.EQ. t_start ) THEN

        CALL check_solve(.TRUE.)

     ELSE

        CALL check_solve(.FALSE.)

     END IF
        
     IF ( verbose_level .GE. 1 ) THEN

        WRITE(*,*) 'cells to solve and reconstruct:' , COUNT(solve_mask)

     END IF

     CALL timestep
     
     IF ( t+dt .GT. t_end ) dt = t_end - t
     IF ( t+dt .GT. t_output ) dt = t_output - t

     IF ( output_runout_flag ) THEN

        IF ( t+dt .GT. t_runout ) dt = t_runout - t

     END IF

     IF ( n_probes .GT. 0 ) THEN

        IF ( t+dt .GT. t_probes ) dt = t_probes - t

     END IF

     dt = MIN(dt,1.1_wp * 0.5_wp * ( dt_old + dt_old_old ) )

     dt_old_old = dt_old
     dt_old = dt

     CALL imex_RK_solver

     CALL update_erosion_deposition_cell(dt)
     
     IF ( topo_change_flag ) CALL topography_reconstruction

     t = t+dt

     !$OMP PARALLEL DO private(j,k,p_dyn,i_table,i_thk_lev,i_pdyn_lev)

     DO l = 1,solve_cells

        j = j_cent(l)
        k = k_cent(l)

        IF ( q(1,j,k) .GT. 0.0_wp ) THEN

           CALL qc_to_qp(q(1:n_vars,j,k) , qp(1:n_vars+2,j,k) , p_dyn )

           hmax(j,k) = MAX( hmax(j,k) , qp(1,j,k) )

           i_table = 0
           
           DO i_thk_lev=1,n_thickness_levels

              thck_table(j,k) = ( qp(1,j,k) .GE. thickness_levels(i_thk_lev) )

              DO i_pdyn_lev=1,n_dyn_pres_levels

                 pdyn_table(j,k) = ( p_dyn .GE. dyn_pres_levels(i_pdyn_lev) ) 

                 i_table = i_table + 1

                 vuln_table(i_table,j,k) = vuln_table(i_table,j,k) .OR.            &
                      ( thck_table(j,k) .AND. pdyn_table(j,k) )

              END DO

           END DO

        ELSE

           qp(1:n_vars+2,j,k) = 0.0_wp

        END IF

     END DO

     !$OMP END PARALLEL DO

     IF ( verbose_level .GE. 0 ) THEN

        WRITE(*,FMT="(A3,F10.4,A5,F9.5,A9,ES11.3E3,A11,ES11.3E3,A9,ES11.3E3,A15,   &
             &ES11.3E3)")                                                          &
             't =',t,'dt =',dt,                                                    &
             ' mass = ',dx*dy*SUM(q(1,:,:)) ,                                      &
             ' volume = ',dx*dy*SUM(qp(1,:,:)) ,                                   &
             ' area = ',dx*dy*COUNT(q(1,:,:).GT.1.D-7) ,                           &
             ' solid mass = ',dx*dy*SUM(q(5:4+n_solid,:,:)) 

     END IF

     IF ( output_runout_flag ) THEN

        IF ( ( t .GE. t_runout ) .OR. ( t .GE. t_steady ) ) THEN

           stop_flag_old = stop_flag

           IF ( output_runout_flag ) CALL output_runout(t,stop_flag)

           IF ( ( stop_flag ) .AND. (.NOT.stop_flag_old) ) THEN

              t_steady = MIN(t_end,t_output)

           END IF

        END IF

     END IF

     IF ( n_probes .GT. 0 ) THEN

        IF ( t .GE. t_probes ) CALL output_probes(t)

     END IF

     IF ( ( t .GE. t_output ) .OR. ( t .GE. t_end ) ) THEN

        CALL cpu_time(t3)
        CALL system_clock(st3)

        IF ( verbose_level .GE. 0 ) THEN

           WRITE(*,*) 'Time taken by iterations is',t3-t2,'seconds'
           WRITE(*,*) 'Elapsed real time = ', DBLE( st3 - st2 ) / rate,'seconds'

        END IF

        IF ( output_cons_flag .OR. output_esri_flag .OR. output_phys_flag )     &
             CALL output_solution(t)

     END IF

  END DO

  CALL deallocate_solver_variables

  CALL close_units

  CALL cpu_time(t3)
  CALL system_clock(st3)

  WRITE(*,*) 'Total time taken by the code is',t3-t1,'seconds'
  WRITE(*,*) 'Total elapsed real time is', DBLE( st3 - st1 ) / rate,'seconds'

END PROGRAM SW_VAR_DENS_MODEL

