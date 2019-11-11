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
!>     Via della Faggiola, 36\n
!>     I-56126 Pisa, Italy \n
!>     E-mail: mattia.demichielivitturi@ingv.it \n
!********************************************************************************

!> \brief Main Program 
PROGRAM SW_VAR_DENS_MODEL

  USE constitutive_2d, ONLY : init_problem_param

  USE geometry_2d, ONLY : init_grid
  USE geometry_2d, ONLY : topography_reconstruction

  USE geometry_2d, ONLY : dx,dy,B_cent

  USE init_2d, ONLY : riemann_problem

  USE inpout_2d, ONLY : init_param
  USE inpout_2d, ONLY : read_param
  USE inpout_2d, ONLY : update_param
  USE inpout_2d, ONLY : output_solution
  USE inpout_2d, ONLY : output_runout
  USE inpout_2d, ONLY : read_solution
  USE inpout_2d, ONLY : close_units

  USE inpout_2d, ONLY : output_runout_flag


  USE solver_2d, ONLY : allocate_solver_variables
  USE solver_2d, ONLY : deallocate_solver_variables
  USE solver_2d, ONLY : imex_RK_solver
  USE solver_2d, ONLY : update_erosion_deposition_cell
  USE solver_2d, ONLY : timestep
  ! USE solver_2d, ONLY : check_solve

  USE inpout_2d, ONLY : restart

  USE parameters_2d, ONLY : t_start
  USE parameters_2d, ONLY : t_end
  USE parameters_2d, ONLY : t_output
  USE parameters_2d, ONLY : t_runout
  USE parameters_2d, ONLY : t_steady
  USE parameters_2d, ONLY : dt0
  USE parameters_2d, ONLY : riemann_flag
  USE parameters_2d, ONLY : topo_change_flag
  USE parameters_2d, ONLY : verbose_level
  USE parameters_2d, ONLY : n_solid


  USE solver_2d, ONLY : q , dt
  USE solver_2d, ONLY : q1max
  
  ! USE solver_2d, ONLY : solve_mask

  IMPLICIT NONE

  REAL*8 :: t
  REAL*8 :: t1 , t2
  REAL*8 :: dt_old , dt_old_old
  LOGICAL :: stop_flag
  LOGICAL :: stop_flag_old

  CALL cpu_time(t1)

  CALL init_param

  CALL read_param

  CALL init_grid

  CALL init_problem_param

  CALL allocate_solver_variables
 
  WRITE(*,*) 'QUI'
 
  IF ( restart ) THEN

     CALL read_solution

  ELSE

     IF( riemann_flag ) THEN

        ! riemann problem defined in file.inp
        CALL riemann_problem

     ENDIF

  END IF

  WRITE(*,*) 'QUI2'

  IF ( topo_change_flag ) CALL topography_reconstruction

  t = t_start
  
  WRITE(*,*) 
  WRITE(*,*) '******** START COMPUTATION *********'
  WRITE(*,*)

  IF ( verbose_level .GE. 1 ) THEN
     
     WRITE(*,*) 'Min q(1,:,:)=',MINVAL(q(1,:,:))
     WRITE(*,*) 'Max q(1,:,:)=',MAXVAL(q(1,:,:))

     WRITE(*,*) 'Min B(:,:)=',MINVAL(B_cent(:,:))
     WRITE(*,*) 'Max B(:,:)=',MAXVAL(B_cent(:,:))


     WRITE(*,*) 'size B_cent',size(B_cent,1),size(B_cent,2)

     WRITE(*,*) 'SUM(q(1,:,:)=',SUM(q(1,:,:))
     WRITE(*,*) 'SUM(B_cent(:,:)=',SUM(B_cent(:,:))
     
  END IF

  q1max(:,:) = q(1,:,:)
  
  dt_old = dt0
  dt_old_old = dt_old
  t_steady = t_end
  stop_flag = .FALSE.


  IF ( output_runout_flag ) CALL output_runout(t,stop_flag)

  CALL output_solution(t)

  IF ( SUM(q(1,:,:)) .EQ. 0.D0 ) t_steady = t_output
  
  WRITE(*,FMT="(A3,F10.4,A5,F9.5,A15,ES12.3E3,A15,ES12.3E3,A19,ES12.3E3)")   &
       't =',t,'dt =',dt0,                                                   &
       ' total mass = ',dx*dy*SUM(q(1,:,:)) ,                              &
       ' total area = ',dx*dy*COUNT(q(1,:,:).GT.1.D-5) ,                     &
       ' total solid mass = ',dx*dy*SUM(q(5:4+n_solid,:,:)) 

  DO WHILE ( ( t .LT. t_end ) .AND. ( t .LT. t_steady ) )

     CALL update_param
     
     ! CALL check_solve
     ! WRITE(*,*) 'cells to solve:',COUNT(solve_mask)

     CALL timestep

     IF ( t+dt .GT. t_end ) dt = t_end - t
     IF ( t+dt .GT. t_output ) dt = t_output - t
     
     IF ( output_runout_flag ) THEN

        IF ( t+dt .GT. t_runout ) dt = t_runout - t

     END IF

     dt = MIN(dt,1.1D0*0.5D0*(dt_old+dt_old_old))
     
     dt_old_old = dt_old
     dt_old = dt

     CALL imex_RK_solver

     CALL update_erosion_deposition_cell(dt)
 
     IF ( topo_change_flag ) CALL topography_reconstruction

     q1max(:,:) = MAX( q1max(:,:) , q(1,:,:) )
     
     t = t+dt
     
     WRITE(*,FMT="(A3,F10.4,A5,F9.5,A15,ES12.3E3,A15,ES12.3E3,A19,ES12.3E3)")&
          't =',t,'dt =',dt,                                                 &
          ' total mass = ',dx*dy*SUM(q(1,:,:)) ,                           &
          ' total area = ',dx*dy*COUNT(q(1,:,:).GT.1.D-7) ,                  &
          ' total solid mass = ',dx*dy*SUM(q(5:4+n_solid,:,:)) 
     
     IF ( output_runout_flag ) THEN

        IF ( ( t .GE. t_runout ) .OR. ( t .GE. t_steady ) ) THEN

           stop_flag_old = stop_flag
           
           IF ( output_runout_flag ) CALL output_runout(t,stop_flag)

           IF ( ( stop_flag ) .AND. (.NOT.stop_flag_old) ) THEN
              
              t_steady = MIN(t_end,t_output)
              t_runout = t_steady

           END IF

        END IF

     END IF

     t_steady = t_end

     IF ( ( t .GE. t_output ) .OR. ( t .GE. t_end ) ) THEN

        CALL output_solution(t)

     END IF

  END DO

  CALL deallocate_solver_variables

  CALL close_units

  CALL cpu_time(t2)

  WRITE(*,*) 'Time taken by the code was',t2-t1,'seconds'

END PROGRAM SW_VAR_DENS_MODEL

