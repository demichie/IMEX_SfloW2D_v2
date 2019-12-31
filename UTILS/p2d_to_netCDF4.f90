! This is part of the netCDF package.
! Copyright 2006 University Corporation for Atmospheric Research/Unidata.
! See COPYRIGHT file for conditions of use.

! This is an example program which writes some 4D pressure and
! temperatures. It is intended to illustrate the use of the netCDF
! fortran 90 API. The companion program pres_temp_4D_rd.f shows how
! to read the netCDF data file created by this program.

! This program is part of the netCDF tutorial:
! http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-tutorial

! Full documentation of the netCDF Fortran 90 API can be found at:
! http://www.unidata.ucar.edu/software/netcdf/docs/netcdf-f90

! $Id: pres_temp_4D_wr.f90,v 1.7 2007/01/24 19:32:10 russ Exp $

program pres_temp_4D_wr
  use netcdf
  implicit none

  INTEGER, PARAMETER :: dp = KIND(1.D0)

  CHARACTER(LEN=40) :: output_file_2d     !< Name of the output files
  ! This is the name of the data file we will create.
  character (len = 40) :: FILE_NAME 
  integer :: ncid

  ! We are writing 4D data, a 2 x 6 x 12 lvl-lat-lon grid, with 2
  ! timesteps of data.
  integer, parameter :: NDIMS = 3, NRECS = 2

  integer :: x_dimid , y_dimid , t_dimid
  integer :: x_varid , y_varid , t_varid

  ! The start and count arrays will tell the netCDF library where to
  ! write our data.
  integer :: start(NDIMS), count(NDIMS)

  ! We will create two netCDF variables, one each for temperature and
  ! pressure fields.
  integer :: B_varid , w_varid , h_varid , temp_varid
  integer :: u_varid , v_varid , rho_m_varid , red_grav_varid
  integer, allocatable :: alfas_varid(:) , deposit_varid(:)
  integer :: dimids(NDIMS)

  ! We recommend that each variable carry a "units" attribute.
  character (len = *), parameter :: UNITS = "units"


  ! Loop indices
  integer :: rec, i ,j ,k , i_solid

  character argv*40
  INTEGER*4 iter, iargc, n_arg

  LOGICAL :: lexist

  REAL(dp) :: t_start           !< initial time for the run
  REAL(dp) :: t_end             !< end time for the run
  REAL(dp) :: dt_output         !< time interval for the output of the solution

  INTEGER :: comp_cells_x      !< Number of control volumes x in the comp. domain
  INTEGER :: comp_cells_y      !< Number of control volumes y in the comp. domain
  REAL(dp) :: x0                 !< Left of the physical domain
  REAL(dp) :: y0                 !< Bottom of the physical domain

  INTEGER :: comp_cells_xx      !< Number of control volumes x in the comp. domain
  INTEGER :: comp_cells_yy      !< Number of control volumes y in the comp. domain

  
  REAL(dp) :: cell_size

  REAL(dp), ALLOCATABLE :: x(:) , y(:), time(:)
  REAL(dp), ALLOCATABLE :: w(:,:)
  REAL(dp), ALLOCATABLE :: h(:,:) , u(:,:) , v(:,:)
  REAL(dp), ALLOCATABLE :: B(:,:) , alphas(:,:,:) , T(:,:)
  REAL(dp), ALLOCATABLE :: rho_m(:,:) , red_grav(:,:) , deposit(:,:,:)

  INTEGER :: n_solid

  INTEGER :: n_output

  CHARACTER(LEN=40) :: run_name           !< Name of the run

  CHARACTER(LEN=4) :: idx_string 

  integer           :: intg(2), stat , i_first , i_last

  WRITE(*,*)
  WRITE(*,*) '------------------------'
  WRITE(*,*) '-- P2D_TO_NETCDF4 1.0 --'
  WRITE(*,*) '------------------------'
  WRITE(*,*)

  intg(2) = -1
  i = 0
  
  n_arg = iargc()


  IF ( n_arg .EQ. 0 ) THEN

     WRITE(*,*) 'ERROR: p2d_to_netCDF4.x requires .bak file as argument'
     WRITE(*,*) 
     STOP

  END IF

  DO iter = 1, n_arg

     CALL getarg( iter, argv )
     
     IF ( iter .EQ. 1 ) THEN
     
        INQUIRE (FILE=argv,exist=lexist)

        IF (lexist .EQV. .FALSE.) THEN

           WRITE(*,*) 'ERROR: bak file does not exist: ',argv
           WRITE(*,*)
           STOP
           
        ELSE
           
           CALL read_bak(argv)
           n_output = nint((t_end-t_start)/dt_output)
             
        END IF

     ELSE

        i = i+1
        call str2int(argv,intg(i),stat)
        if ( stat == 0 ) then

           IF ( ( intg(i) .GE. 0 ) .AND. ( intg(i) .LE. n_output ) ) THEN

              ! print *,intg(i)

              IF ( i.EQ.1 ) i_first = intg(i)
              
              IF ( i.EQ.2 ) THEN
                 i_last = intg(i)
                 IF ( i_last .LT. i_first ) THEN
                    print *,'Second index smaller than first!'
                    stop
                 endif
              END IF
                 
           ELSE

              print *,'Index should be between ',0,'and',n_output
              stop

           END IF

        else
           print *,'Conversion of string ',TRIM(argv),' to integer failed!'
           stop
        endif

     END IF
        
  END DO

  IF ( i.EQ.0 ) THEN
     i_first = 0
     i_last = n_output
  ELSEIF ( i.EQ.1 ) THEN
     i_last = i_first
  END IF

  WRITE(*,*) 'first index',i_first
  WRITE(*,*) 'last index',i_last

  FILE_NAME = TRIM(run_name)//'.nc'

  ALLOCATE( time(i_last-i_first+1) )
  DO i=i_first,i_last
     time(i-i_first+1) = dt_output*i
  END DO

  WRITE(*,*) 'time= ',time

  comp_cells_xx = MAX( comp_cells_x,2)
  comp_cells_yy = MAX( comp_cells_y,2)
  
  ALLOCATE( x(comp_cells_xx) )
  DO j=1,comp_cells_xx
     x(j) = x0 + ( j - 0.5) * cell_size
  END DO

  
  
  ALLOCATE( y(comp_cells_yy) )
  DO k=1,comp_cells_yy
     y(k) = y0 + ( k - 0.5) * cell_size
  END DO

  WRITE(*,*) comp_cells_x,comp_cells_y,n_solid
 
  ALLOCATE( w(comp_cells_xx,comp_cells_yy) )
  ALLOCATE( h(comp_cells_xx,comp_cells_yy) )
  ALLOCATE( u(comp_cells_xx,comp_cells_yy) )
  ALLOCATE( v(comp_cells_xx,comp_cells_yy) )
  ALLOCATE( B(comp_cells_xx,comp_cells_yy) )
  ALLOCATE( alphas(comp_cells_xx,comp_cells_yy,n_solid) )
  ALLOCATE( T(comp_cells_xx,comp_cells_yy) )
  ALLOCATE( rho_m(comp_cells_xx,comp_cells_yy) )
  ALLOCATE( red_grav(comp_cells_xx,comp_cells_yy) )
  ALLOCATE( deposit(comp_cells_xx,comp_cells_yy,n_solid) )


  ALLOCATE( alfas_varid(n_solid) )
  ALLOCATE( deposit_varid(n_solid) )


  
  ! Create the file. 
  call check( nf90_create(FILE_NAME, nf90_HDF5 , ncid) )

  ! Define the dimensions. The record dimension is defined to have
  ! unlimited length - it can grow as needed. In this example it is
  ! the time dimension.
  call check( nf90_def_dim(ncid, 'x', comp_cells_xx, x_dimid) )
  call check( nf90_def_dim(ncid, 'y', comp_cells_yy, y_dimid) )
  call check( nf90_def_dim(ncid, 'time', NF90_UNLIMITED, t_dimid) )

  ! Define the coordinate variables. We will only define coordinate
  ! variables for lat and lon.  Ordinarily we would need to provide
  ! an array of dimension IDs for each variable's dimensions, but
  ! since coordinate variables only have one dimension, we can
  ! simply provide the address of that dimension ID (x_dimid) and
  ! similarly for (y_dimid).
  call check( nf90_def_var(ncid, 'x', NF90_REAL, x_dimid, x_varid) )
  call check( nf90_def_var(ncid, 'y', NF90_REAL, y_dimid, y_varid) )
  call check( nf90_def_var(ncid, 'time', NF90_REAL, t_dimid, t_varid) )

  ! Assign units attributes to coordinate variables.
  call check( nf90_put_att(ncid, x_varid, UNITS, 'meters') )
  call check( nf90_put_att(ncid, y_varid, UNITS, 'meters') )
  call check( nf90_put_att(ncid, t_varid, UNITS, 'seconds') )

  ! The dimids array is used to pass the dimids of the dimensions of
  ! the netCDF variables. Both of the netCDF variables we are creating
  ! share the same four dimensions. In Fortran, the unlimited
  ! dimension must come last on the list of dimids.
  dimids = (/ x_dimid, y_dimid, t_dimid /)

  ! Define the netCDF variables for the output data.
  call check( nf90_def_var(ncid, 'B', NF90_REAL, dimids, B_varid, deflate_level = 7) )
  call check( nf90_def_var(ncid, 'w', NF90_REAL, dimids, w_varid, deflate_level = 7) )
  call check( nf90_def_var(ncid, 'h', NF90_REAL, dimids, h_varid, deflate_level = 7) )
  call check( nf90_def_var(ncid, 'Temp', NF90_REAL, dimids, temp_varid, deflate_level = 7) )
  call check( nf90_def_var(ncid, 'u', NF90_REAL, dimids, u_varid, deflate_level = 7) )
  call check( nf90_def_var(ncid, 'v', NF90_REAL, dimids, v_varid, deflate_level = 7) )
  call check( nf90_def_var(ncid, 'rho_m', NF90_REAL, dimids, rho_m_varid, deflate_level = 7) )
  call check( nf90_def_var(ncid, 'red_grav', NF90_REAL, dimids, red_grav_varid, deflate_level = 7) )

  do i_solid=1,n_solid

     idx_string = lettera(i_solid)
     call check( nf90_def_var(ncid, 'alfas_'//idx_string, NF90_REAL, dimids,     &
          alfas_varid(i_solid), deflate_level = 7) )
     
     call check( nf90_def_var(ncid, 'deposit_'//idx_string, NF90_REAL, dimids,   &
          deposit_varid(i_solid) , deflate_level = 7) )

  end do

  ! Assign units attributes to the netCDF variables.
  call check( nf90_put_att(ncid, B_varid, UNITS, 'meters') )
  call check( nf90_put_att(ncid, w_varid, UNITS, 'meters') )
  call check( nf90_put_att(ncid, h_varid, UNITS, 'meters') )
  call check( nf90_put_att(ncid, temp_varid, UNITS, 'Kelvins') )
  call check( nf90_put_att(ncid, u_varid, UNITS, 'meters/seconds') )
  call check( nf90_put_att(ncid, v_varid, UNITS, 'meters/seconds') )
  call check( nf90_put_att(ncid, rho_m_varid, UNITS, 'kg/m3') )
  call check( nf90_put_att(ncid, red_grav_varid, UNITS, 'meters/seconds^2') )

  do i_solid=1,n_solid

     call check( nf90_put_att(ncid, alfas_varid(i_solid), UNITS, '') )
     call check( nf90_put_att(ncid, deposit_varid(i_solid), UNITS, 'meters') )
     
  end do

  
  ! End define mode.
  call check( nf90_enddef(ncid) )

  ! Write the coordinate variable data. This will put the latitudes
  ! and longitudes of our data grid into the netCDF file.
  call check( nf90_put_var(ncid, x_varid, x) )
  call check( nf90_put_var(ncid, y_varid, y) )
  call check( nf90_put_var(ncid, t_varid, time) )

  ! These settings tell netcdf to write one timestep of data. (The
  ! setting of start(4) inside the loop below tells netCDF which
  ! timestep to write.)
  count = (/ comp_cells_xx, comp_cells_yy, 1 /)
  start = (/ 1, 1, 1 /)

  ! Write the pretend data. This will write our surface pressure and
  ! surface temperature data. The arrays only hold one timestep worth
  ! of data. We will just rewrite the same data for each timestep. In
  ! a real :: application, the data would change between timesteps.

  do rec = i_first, i_last

     idx_string = lettera(rec)

     output_file_2d = TRIM(run_name)//'_'//idx_string//'.p_2d'
     WRITE(*,*) 'Reading file ',output_file_2d

     OPEN(12,FILE=output_file_2d,STATUS='old')

     DO k = 1,comp_cells_y

        DO j = 1,comp_cells_x

           READ(12,*) x(j), y(k), h(j,k) , u(j,k) , v(j,k) , B(j,k) , w(j,k) ,  &
                alphas(j,k,1:n_solid) , T(j,k) , rho_m(j,k) , red_grav(j,k) ,   &
                deposit(j,k,1:n_solid)

        END DO

     END DO

     IF ( comp_cells_x .EQ. 1 ) THEN

        h(2,:) = h(1,:)
        v(2,:) = v(1,:)
        B(2,:) = B(1,:)
        w(2,:) = W(1,:)
        T(2,:) = T(1,:)
        rho_m(2,:) = rho_m(1,:)
        red_grav(2,:) = red_grav(1,:)
        alphas(2,:,1:n_solid) = alphas(1,:,1:n_solid)
        deposit(2,:,1:n_solid) = deposit(1,:,1:n_solid)
        
     END IF

     IF ( comp_cells_y .EQ. 1 ) THEN

        h(:,2) = h(:,1)
        v(:,2) = v(:,1)
        B(:,2) = B(:,1)
        w(:,2) = W(:,1)
        T(:,2) = T(:,1)
        rho_m(:,2) = rho_m(:,1)
        red_grav(:,2) = red_grav(:,1)
        alphas(:,2,1:n_solid) = alphas(:,1,1:n_solid)
        deposit(:,2,1:n_solid) = deposit(:,1,1:n_solid)
        
     END IF
        
     CLOSE(12)

     start(3) = rec-i_first+1

     call check( nf90_put_var(ncid, B_varid, B, start = start, count=count) )

     call check( nf90_put_var(ncid, w_varid, w, start = start, count=count) )

     call check( nf90_put_var(ncid, h_varid, h, start = start, count=count) )

     call check( nf90_put_var(ncid, temp_varid, T, start = start, count=count) )

     call check( nf90_put_var(ncid, u_varid, u, start = start, count = count) )

     call check( nf90_put_var(ncid, v_varid, v, start = start, count = count) )
     
     call check( nf90_put_var(ncid, rho_m_varid, rho_m, start = start,          &
          count = count) )

     call check( nf90_put_var(ncid, red_grav_varid, red_grav, start = start,    &
          count = count) )

     do i_solid=1,n_solid

        call check( nf90_put_var(ncid, alfas_varid(i_solid),alphas(:,:,i_solid),&
             start = start, count = count) )

        call check( nf90_put_var(ncid , deposit_varid(i_solid),                 &
             deposit(:,:,i_solid) , start = start , count = count) )
     
     end do

  end do

  ! Close the file. This causes netCDF to flush all buffers and make
  ! sure your data are really written to disk.
  call check( nf90_close(ncid) )

  print *,"*** SUCCESS writing example file ", FILE_NAME, "!"

contains

  subroutine check(status)
    integer, intent ( in) :: status

    if(status /= nf90_noerr) then 
       print *, trim(nf90_strerror(status))
       stop "Stopped"
    end if
  end subroutine check


  subroutine read_bak(bak_name)

    IMPLICIT NONE

    INTEGER, PARAMETER :: backup_unit = 8      !< Backup input data unit

    INTEGER :: ios

    CHARACTER(LEN=40), INTENT(IN) :: bak_name           !< Backup file for the parameters


    !> Flag to start a run from a previous output:\n
    !> - T     => Restart from a previous output (.asc or .q_2d)
    !> - F     => Restart from initial condition read from two_phases.inp
    !> .
    LOGICAL :: restart


    !> Flag to save the output in esri ascii format *.asc
    !> - T     => write esri file
    !> - F     => do not write esri file
    !> .
    LOGICAL :: output_esri_flag
    !> Flag to save the physical variables on file *.p_2d
    !> - T     => write physical variables on file
    !> - F     => do not write the physical variables
    !> .
    LOGICAL :: output_phys_flag
    !> Flag to save the conservative variables on file *.q_2d
    !> - T     => write conservative variables on file
    !> - F     => do not write the conservative variables
    !> .
    LOGICAL :: output_cons_flag
    !> Flag to save the max runout at ouput times
    !> - T     => write max runout on file
    !> - F     => do not write max runout
    !> .
    LOGICAL :: output_runout_flag

    INTEGER :: verbose_level

    !> Flag to choose if we add the rheology
    !> - T      => rheology activated
    !> - F      => no rheology
    !> .
    LOGICAL :: rheology_flag

    !> Flag to choose the sort of problem to solve
    !> - T      => riemann problem
    !> - F      => generic initial conditions (uploaded through functions, to be defined in inpout_2d.f90)
    !> .
    LOGICAL :: riemann_flag

    !> Flag to choose the equation for temperature to solve
    !> - T      => solve the full energy equation
    !> - F      => solve for a simpler transport equation (advection) for temperature
    !> .
    LOGICAL :: energy_flag

    LOGICAL :: liquid_flag
    LOGICAL :: gas_flag

    LOGICAL :: topo_change_flag
    LOGICAL :: radial_source_flag

    LOGICAL :: collapsing_volume_flag

    REAL(dp) :: rho0_s(1000) , diam0_s(1000) , sp_heat0_s(1000), erosion_coeff0(1000)
    LOGICAL :: settling_flag
    REAL(dp) :: T_s_substrate

    NAMELIST / run_parameters / run_name , restart , t_start , t_end , dt_output ,&
         output_cons_flag , output_esri_flag , output_phys_flag ,                 &
         output_runout_flag , verbose_level

    NAMELIST / newrun_parameters / x0 , y0 , comp_cells_x , comp_cells_y ,        &
         cell_size , rheology_flag , riemann_flag , energy_flag , liquid_flag ,   &
         radial_source_flag , collapsing_volume_flag , topo_change_flag , gas_flag

    NAMELIST / solid_transport_parameters / n_solid , rho0_s , diam0_s ,          &
         sp_heat0_s , erosion_coeff0 , settling_flag , T_s_substrate


    OPEN(backup_unit,FILE=bak_name,STATUS='old')

    ! ------- READ run_parameters NAMELIST -----------------------------------
    READ(backup_unit, run_parameters,IOSTAT=ios )

    IF ( ios .NE. 0 ) THEN
       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist RUN_PARAMETERS'
       WRITE(*,*) 'Please check the input file'
       STOP
    ELSE
       WRITE(*,*) 'Run name: ',run_name
       REWIND(backup_unit)
    END IF

    ! ------- READ newrun_parameters NAMELIST --------------------------------
    READ(backup_unit,newrun_parameters,IOSTAT=ios)
    IF ( ios .NE. 0 ) THEN
       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist NEWRUN_PARAMETERS'
       WRITE(*,*) 'Please check the input file'
       STOP
    ELSE
       REWIND(backup_unit)
    END IF

    ! ------- READ solid_transport_parameters NAMELIST --------------------------
    READ(backup_unit, solid_transport_parameters,IOSTAT=ios)

    IF ( ios .NE. 0 ) THEN

       WRITE(*,*) 'IOSTAT=',ios
       WRITE(*,*) 'ERROR: problem with namelist SOLID_TRANSPORT_PARAMETERS'
       WRITE(*,*) 'Please check the input file'
       STOP

    ELSE

       REWIND(backup_unit)

    END IF

  end subroutine read_bak

  CHARACTER*4 FUNCTION lettera(k)
    IMPLICIT NONE
    CHARACTER ones,tens,hund,thou
    !
    INTEGER :: k
    !
    INTEGER :: iten, ione, ihund, ithou
    !
    ithou=INT(k/1000)
    ihund=INT((k-(ithou*1000))/100)
    iten=INT((k-(ithou*1000)-(ihund*100))/10)
    ione=k-ithou*1000-ihund*100-iten*10
    ones=CHAR(ione+48)
    tens=CHAR(iten+48)
    hund=CHAR(ihunD+48)
    thou=CHAR(ithou+48)
    lettera=thou//hund//tens//ones
    !
    RETURN
  END FUNCTION lettera

  elemental subroutine str2int(str,int,stat)
    implicit none
    ! Arguments
    character(len=*),intent(in) :: str
    integer,intent(out)         :: int
    integer,intent(out)         :: stat
    
    read(str,*,iostat=stat)  int
  end subroutine str2int
  
end program pres_temp_4D_wr

