! Description:
! This module contains code for particle deposition problems
! and is called by the species transport module. 

module particle

  use other_consts
  use indices
  use arrays
  implicit none
  
  !Module Parameters

  !Module types

  !Module variables

  !Interfaces
  private airway_mesh_deform
  
  public initialise_particles
  public controller_particlesolve
  
  public where_inlist, inlist
  public calc_mass_particles

  real(dp), public :: total_volume_change

contains

!!!##########################################################################
!!!##########################################################################

!!! Assemble matrices for 1D particle transport, and solve. Based on the gas mixing model
!!! with adaptations by Falko Schmidt for particles.


  subroutine initialise_particles(fileid)
 !DEC$ ATTRIBUTES DLLEXPORT,ALIAS:"SO_INITIALISE_TRANSPORT" :: INITIALISE_TRANSPORT

    use arrays, only: dp,particle_parameters
!   use ventilation,only: evaluate_vent_step
    use diagnostics, only: enter_exit

    implicit none
  
    type(particle_parameters) :: part_param
    integer,intent(in) :: fileid
    real(dp) :: err
    ! Local variables
    real(dp),allocatable :: solution(:)

!   integer :: MatrixSize,nonzeros,ncol,nentry, &
!        noffset_entry,noffset_row,np,nrow,nrow_BB
!   real(dp) :: AA,BB,current_mass,current_volume,inlet_flow, &
!        mass_deposit,mass_error,theta,time,volume_error,volume_tree,&
!        cpu_dt_start,cpu_dt_end,temp_mass

!   real(dp) :: inlet_concentration(3)
!   real(dp) :: initial_volume, dia,  &
!      diffusion_coeff, mu, rho, gravityx, gravityy, gravityz
!   logical :: coupled, last_breath, inspiration
!   real(dp) :: solve_tolerance 
!   real(dp) :: time_inspiration,time_breath_hold,time_expiration, dt_gm
!   real(dp) :: VtotTLC, totacinarLength, lambda, kBoltz, Temperature
!   integer :: num_brths_gm


    logical :: carryon
    character(len=60) :: sub_name
    
    ! #############################################################################

!   set_diagnostics = .true.
    sub_name = 'initialise_particles'
    call enter_exit(sub_name,1)

    call read_params_subjectinfo(part_param%lung_root,part_param%results_location,& 
         part_param%study,part_param%subject,part_param%protocol,&
         part_param%group_name)
    call read_params_particles(part_param%initial_volume, part_param%inlet_concentration, &
       part_param%diffusion_coeff, part_param%mu, part_param%prho, part_param%gravityx,&
       part_param%gravityy, part_param%gravityz, part_param%pdia)
    call read_params_solver(part_param%coupled, part_param%solve_tolerance, &
       part_param%last_breath, part_param%inspiration)
    call read_params_time(part_param%num_brths_gm,part_param%time_inspiration,&
       part_param%time_breath_hold,part_param%time_expiration, part_param%dt_gm)

    call read_params_advanced(part_param%VtotTLC, part_param%totacinarLength, &
         part_param%lambda, part_param%kBoltz, part_param%Temperature)

    
    call enter_exit(sub_name,2)
    
  end subroutine initialise_particles
  

!!!##########################################################################
!!!##########################################################################
!!!##########################################################################

  subroutine controller_particlesolve(fileid)

!!!##########################################################################
!!!##########################################################################
!!!##########################################################################

!!! Assemble matrices for 1D particle transport, and solve. Based on the gas mixing model
!!! with adaptations by Falko Schmidt for particles.

!   use arrays,only: dp,node_field,node_field,elem_field,&
!       particle_parameters
!   use exports,only: export_node_field
!   use geometry, only: volume_of_mesh
!   use indices, only: nj_conc1, nj_conc2, nj_conc3

    use arrays
    use exports
    use geometry
    use indices
    use solve,only: pmgmres_ilu_cr
!   use ventilation,only: evaluate_vent_step
    use gasmix,only: reduce_gasmix, initial_gasmix
    use other_consts, only: MAX_FILENAME_LEN
    use diagnostics, only: enter_exit

    
    implicit none
  
    type(particle_parameters) :: part_param
    integer,intent(in) :: fileid
    real(dp) :: err
    character(len=9) :: problem_type = 'particles'
    character(len=200) :: op_name

    ! Local variables
    real(dp),allocatable :: solution(:)

    integer :: MatrixSize,nonzeros,ncol,nentry, &
         noffset_entry,noffset_row,np,nrow,nrow_BB
    real(dp) :: AA,BB,current_mass,current_volume,inlet_flow, &
         mass_deposit,mass_error,theta,time,volume_error,volume_tree,&
         cpu_dt_start,cpu_dt_end,temp_mass
    
    real(dp) :: ttime, turnovers, time_start, time_end
    integer :: Nexp, nstep, nbreath, j
    logical :: carryon
    character(len=60) :: sub_name
    logical :: last_breath = .false.
    character(len=4) :: char_int
    character(len=MAX_FILENAME_LEN) :: file_export


    
!   set_diagnostics = .true.
    sub_name = 'controller_particlesolve'
    call enter_exit(sub_name,1)


! setup of initial parameters and characters
!##################################################
     op_name = part_param%results_location(1:len_trim(part_param%results_location))//&
     part_param%subject(1:len_trim(part_param%subject))//'_particle'
     print *, 'op_name is', op_name	
     part_param%initial_concentration = 0.0_dp  ! alveolar resident particles = 0!



!!! ###########  INITIAL & BOUNDARY CONDITIONS FOR GAS MIXING & EXCHANGE   ###########
! sets initial concentration at all nodes, concentration at the inlet node
! (for inspiration), and allocates memory for solution matrices.

    call initial_gasmix(part_param%initial_concentration,part_param%inlet_concentration)


     turnovers = part_param%num_brths_gm * part_param%tidal_volume / &
                 (part_param%FRC*1.0e+6_dp)

     if(turnovers.lt.6.0_dp)then
        write(*,'(''Warning:'',i3,'' breaths equals'',f5.2,'' turnovers;'// &
           ' set num_brths_gm to'',i3,'' or enter return to continue:'')') &
           part_param%num_brths_gm,turnovers,int(6.0_dp*(part_param%FRC*1.0e+6_dp)/&
           part_param%tidal_volume)
        read(*,*)
     endif

    time_end = 0.0_dp


!   if(.not.coupled)then
       write(*,'(''   Time    Flow    D_Vol    Vol    volEr  MassErr  C_inlet'')')
       write(*,'(''   (s)    (L/s)     (L)     (L)      %        %     (frac) '')')
!   else
!      write(*,'(''   Time    Flow    D_Vol    Vol    volEr  MassErr  C_inlet   dMass  dMass*   Mass   Mass*'')')
!      write(*,'(''   (s)    (L/s)     (L)     (L)      %        %     (frac) '')')

!  endif


    ttime = 0.0_dp
    time = 0.0_dp

!###########################################################
    do nbreath = 1,part_param%num_brths_gm
!###########################################################

!!! Inspiration
     inlet_flow = abs(elem_field(ne_flow,1))
!if(.not.coupled) call scale_flow_field(inlet_flow)
     call scale_flow_field(inlet_flow)
     time_start = time_end
     time_end = time_start + part_param%time_inspiration

     node_field(nj_conc1,1) = part_param%inlet_concentration(1) ! need to set here
!    node_field(nj_conc2,1) = part_param%inlet_concentration(2) ! reset the inlet CO2 concentration
!    node_field(nj_conc3,1) = part_param%inlet_concentration(3) ! reset the inlet CO concentration

      inlet_flow = 250.0e+3_dp ! 250 ml/s, for 2 s
      call scale_flow_field(inlet_flow)
      call solve_particles(fileid,time_end,time_start,.true.,last_breath)
      nstep = nbreath


!!! Breath Hold
!!----------------
     if(part_param%time_breath_hold.gt.0.0_dp)then
        do j = 1,part_param%n_export
           time_start = time_end
           time_end = time_start + part_param%time_breath_hold

           call solve_particles(fileid,time_end,time_start,.true.,last_breath)

           Nexp = Nexp+1
           write(char_int,'(i3)') Nexp

           file_export = trim(op_name)//'_conc_'//trim(adjustl(char_int))
           call export_node_field(nj_conc1,file_export, part_param%group_name, 'conc')
        enddo
     endif


!!! Expiration
!!----------------
!    if(.not.coupled)then
        inlet_flow = -elem_field(ne_flow,1)
        call scale_flow_field(inlet_flow)
!    endif

     time_start = time_end
     time_end = time_start + part_param%time_expiration

     write(*,'('' Exhaling breath'',i4,'' from'',f8.3,'' s to'',f8.3,'' s'')') nbreath, &
          time_start,time_end

! HBK _CAUTION COMMENT
!    if(nbreath.eq.part_param%num_brths_gm.and.nstep.eq.nstep_end) last_breath = .true.

        inlet_flow = -250.0e+3_dp ! 250 ml/s, for 2 s
        call scale_flow_field(inlet_flow)

        call solve_particles(fileid,time_end,time_start,.true.,last_breath)

        nstep = nbreath

     enddo

!   close(fileid)

    call enter_exit(sub_name,2)

  end subroutine controller_particlesolve
  
      
!!!###########################################################################
!   call read_params_subjectinfo(lung_root,results,study,subject,protocol,group_name)
  subroutine read_params_subjectinfo(lung_root,results_location,study,subject,& 
             protocol,group_name)

    use arrays,only: dp
    use diagnostics, only: enter_exit
    implicit none

    character,intent(out) :: lung_root*(*), results_location*(*), study*(*)
    character,intent(out) :: subject*(*), protocol*(*), group_name*(*)

    ! Input related variables
    character(len=100) :: buffer, label
    integer :: pos
    integer, parameter :: fh = 15
    integer :: ios
    integer :: line

    character(len=60) :: sub_name

    ! ###########################################################################

    ios = 0
    line = 0
    sub_name = 'read_params_subjectinfo'
    call enter_exit(sub_name,1)

    ! following values are examples from control.txt

    open(fh, file='Parameters/params_subjectinfo.txt')

    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.

    do while (ios == 0)
       read(fh, '(A)', iostat=ios) buffer
       if (ios == 0) then
          line = line + 1

          ! Find the first instance of whitespace.  Split label and data.
          pos = scan(buffer, '    ')
          label = buffer(1:pos)
          buffer = buffer(pos+2:)
          !note, different read formatting than other read_param functions required
          !due to / characters in the path strings. 
          select case (label)
          case ('lung_root')
             read(buffer, '(A)', iostat=ios) lung_root
             print *, 'Read lung root: ', lung_root(1:len_trim(lung_root))
          case ('results_location')
             read(buffer, '(A)', iostat=ios) results_location
             print *, 'Read output location: ', results_location(1:len_trim(results_location))
          case ('study')
             read(buffer, *, iostat=ios) study
             print *, 'Read study name: ', study
          case ('subject')
             read(buffer, *, iostat=ios) subject
             print *, 'Read subject: ', subject
          case ('protocol')
             read(buffer, *, iostat=ios) protocol
             print *, 'Read protocol: ', protocol
          case ('group_name')
             read(buffer, *, iostat=ios) group_name
             print *, 'Read Group name: ', group_name
          case default
             print *, 'Skipping invalid label at line', line
          end select
       end if
    end do

    close(fh)
    call enter_exit(sub_name,2)

  end subroutine read_params_subjectinfo



!!!###########################################################################
!   call read_params_particles(initial_volume, inlet_concentration, &
!      diffusion_coeff, mu, rho, gravityx, gravity, gravityz)
  subroutine read_params_particles(initial_volume, inlet_concentration, &
       diffusion_coeff, mu, particlerho, gravityx, gravityy, gravityz, particledia)

    use arrays,only: dp
    use diagnostics, only: enter_exit
    implicit none

    real(dp),intent(out) :: initial_volume, &
         diffusion_coeff, mu, particlerho, gravityx, gravityy, gravityz, particledia
    real(dp),intent(out) :: inlet_concentration(3)

    ! Input related variables
    character(len=100) :: buffer, label
    integer :: pos
    integer, parameter :: fh = 15
    integer :: ios
    integer :: line

    character(len=60) :: sub_name

    ! ###########################################################################

    ios = 0
    line = 0
    sub_name = 'read_params_particles'
    call enter_exit(sub_name,1)

    ! following values are examples from control.txt

    open(fh, file='Parameters/params_particles.txt')

    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.

    do while (ios == 0)
       read(fh, '(A)', iostat=ios) buffer
       if (ios == 0) then
          line = line + 1

          ! Find the first instance of whitespace.  Split label and data.
          pos = scan(buffer, '    ')
          label = buffer(1:pos)
          buffer = buffer(pos+1:)

          select case (label)
          case ('volume')
             read(buffer, *, iostat=ios) initial_volume
             print *, 'Read FRC: ', initial_volume
          case ('inlet_concentration')
             read(buffer, *, iostat=ios) inlet_concentration(1)
             print *, 'Read inlet concen: ', inlet_concentration(1)
          case ('particle_dia')
             read(buffer, *, iostat=ios) particledia
             print *, 'Read particle dia: ', particledia
          case ('diffusion_coeff')
             read(buffer, *, iostat=ios) diffusion_coeff
             print *, 'Read diffusion coeff: ', diffusion_coeff
          case ('mu')
             read(buffer, *, iostat=ios) mu
             print *, 'Read fluid viscosity: ', mu
          case ('rho')
             read(buffer, *, iostat=ios) particlerho
             print *, 'Read density: ', particlerho
          case ('gravityx')
             read(buffer, *, iostat=ios) gravityx
             print *, 'Read X gravity: ', gravityx 
          case ('gravityy')
             read(buffer, *, iostat=ios) gravityy
             print *, 'Read Y gravity: ', gravityy 
          case ('gravityz')
             read(buffer, *, iostat=ios) gravityz
             print *, 'Read Z gravity: ', gravityz
          case default
             print *, 'Skipping invalid label at line', line
          end select
       end if
    end do

    close(fh)
    call enter_exit(sub_name,2)

  end subroutine read_params_particles

!!!###########################################################################
  subroutine read_params_solver(coupled, solve_tolerance, last_breath, &
       inspiration)

    use arrays,only: dp
    use diagnostics, only: enter_exit
    implicit none

    real(dp),intent(out) :: solve_tolerance
    logical, intent(out) :: coupled, inspiration, last_breath

    ! Input related variables
    character(len=100) :: buffer, label
    integer :: pos
    integer, parameter :: fh = 15
    integer :: ios
    integer :: line

    character(len=60) :: sub_name

    ! ###########################################################################

    ios = 0
    line = 0
    sub_name = 'read_params_solver'
    call enter_exit(sub_name,1)

    ! following values are examples from control.txt
!    coupled          false
!    solve_tolerance  1.e-08
!    last_breath      false
!    inspiration      true


    open(fh, file='Parameters/params_solver.txt')

    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.

    do while (ios == 0)
       read(fh, '(A)', iostat=ios) buffer
       if (ios == 0) then
          line = line + 1

          ! Find the first instance of whitespace.  Split label and data.
          pos = scan(buffer, '    ')
          label = buffer(1:pos)
          buffer = buffer(pos+1:)

          select case (label)
          case ('coupled')
             read(buffer, *, iostat=ios) coupled
             print *, 'Read coupled: ', coupled
          case ('solve_tolerance')
             read(buffer, *, iostat=ios) solve_tolerance
             print *, 'Read solve tolerance: ', solve_tolerance
          case ('last_breath')
             read(buffer, *, iostat=ios) last_breath
             print *, 'Read last breath: ', last_breath
          case ('inspiration')
             read(buffer, *, iostat=ios) inspiration
             print *, 'Read inspiration: ', inspiration
          case default
             print *, 'Skipping invalid label at line', line
          end select
       end if
    end do

    close(fh)
    call enter_exit(sub_name,2)

  end subroutine read_params_solver

!!!###########################################################################
! subroutine read_params_time(time_inspiration,time_breath_hold,time_expiration,&
!      dt_gm)
  subroutine read_params_time(num_brths_gm,time_inspiration,time_breath_hold,&
       time_expiration, dt_gm)

    use arrays,only: dp
    use diagnostics, only: enter_exit
    implicit none

    real(dp),intent(out) :: time_inspiration,time_breath_hold,&
       time_expiration,dt_gm
    integer,intent(out) :: num_brths_gm

    ! Input related variables
    character(len=100) :: buffer, label
    integer :: pos
    integer, parameter :: fh = 15
    integer :: ios
    integer :: line

    character(len=60) :: sub_name

    ! ###########################################################################

    ios = 0
    line = 0
    sub_name = 'read_params_time'
    call enter_exit(sub_name,1)

    open(fh, file='Parameters/params_time.txt')

    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.

    do while (ios == 0)
       read(fh, '(A)', iostat=ios) buffer
       if (ios == 0) then
          line = line + 1

          ! Find the first instance of whitespace.  Split label and data.
          pos = scan(buffer, '    ')
          label = buffer(1:pos)
          buffer = buffer(pos+1:)

          select case (label)
          case ('time_inspiration')
             read(buffer, *, iostat=ios) time_inspiration
             print *, 'Inspiration time: ', time_inspiration
          case ('time_breath_hold')
             read(buffer, *, iostat=ios) time_breath_hold
             print *, 'Read time for breath hold: ', time_breath_hold
          case ('time_expiration')
             read(buffer, *, iostat=ios) time_expiration
             print *, 'Read time for expiration: ', time_expiration
          case ('timestep')
             read(buffer, *, iostat=ios) dt_gm
             print *, 'Read timestep: ', dt_gm
          case default
             print *, 'Skipping invalid label at line', line
          end select
       end if
    end do

    close(fh)
    call enter_exit(sub_name,2)

  end subroutine read_params_time

!!!###########################################################################
  subroutine read_params_advanced(VtotTLC, totacinarLength, lambda, kBoltz, Temperature)

    use arrays,only: dp
    use diagnostics, only: enter_exit
    implicit none

    real(dp),intent(out) :: VtotTLC, totacinarLength, lambda, kBoltz, Temperature

    ! Input related variables
    character(len=100) :: buffer, label
    integer :: pos
    integer, parameter :: fh = 15
    integer :: ios
    integer :: line

    character(len=60) :: sub_name

    ! ###########################################################################

    ios = 0
    line = 0
    sub_name = 'read_params_advanced'
    call enter_exit(sub_name,1)


    open(fh, file='Parameters/params_advanced.txt')

    ! ios is negative if an end of record condition is encountered or if
    ! an endfile condition was detected.  It is positive if an error was
    ! detected.  ios is zero otherwise.

    do while (ios == 0)
       read(fh, '(A)', iostat=ios) buffer
       if (ios == 0) then
          line = line + 1

          ! Find the first instance of whitespace.  Split label and data.
          pos = scan(buffer, '    ')
          label = buffer(1:pos)
          buffer = buffer(pos+1:)

          select case (label)
          case ('VtotTLC')
             read(buffer, *, iostat=ios) VtotTLC
             print *, 'Read total acinar vol: ', VtotTLC
          case ('totacinarLength')
             read(buffer, *, iostat=ios) totacinarLength
             print *, 'Read total acinar length: ', totacinarLength
          case ('lambda')
             read(buffer, *, iostat=ios) lambda
             print *, 'Read mean free path: ', lambda
          case ('kBoltz')
             read(buffer, *, iostat=ios) kBoltz
             print *, 'Read Boltz constant: ', kBoltz
          case ('Temperature')
             read(buffer, *, iostat=ios) Temperature
             print *, 'Read Temperature(K): ', Temperature
          case default
             print *, 'Skipping invalid label at line', line
          end select
       end if
    end do

    close(fh)
    call enter_exit(sub_name,2)

  end subroutine read_params_advanced

!###################################################################################

  subroutine solve_particles(fileid,time_end,time_start,inspiration,last_breath)

!!! Assemble matrices for 1D particle transport, and solve. Based on the gas mixing model
!!! with adaptations by Falko Schmidt for particles.

!   use arrays, only: dp,particle_parameters,part_acinus_field,node_field,num_nodes,num_elems
    use arrays
    use exports,only: export_node_field
    use geometry, only: volume_of_mesh
    use solve,only: pmgmres_ilu_cr
!   use ventilation,only: evaluate_vent_step
    use gasmix,only: reduce_gasmix, assemble_particles, calc_mass
    use other_consts
    use diagnostics, only: enter_exit
    
    implicit none
  
    type(particle_parameters) :: part_param
    integer,intent(in) :: fileid
    real(dp),intent(in) :: time_end,time_start
    logical,intent(in) :: inspiration,last_breath

    real(dp) :: err
!   character(len=9) :: problem_type = 'particles'
    ! Local variables
    real(dp),allocatable :: solution(:)

    integer :: MatrixSize,nonzeros,ncol,nentry, &
         noffset_entry,noffset_row,np,nrow,nrow_BB
    real(dp) :: AA,BB,current_mass,current_volume,inlet_flow, &
         mass_deposit,mass_error,theta,time,volume_error,volume_tree,&
         cpu_dt_start,cpu_dt_end,temp_mass,ideal_mass
    
    real(dp),allocatable :: global_K(:),global_M(:),global_AA(:),global_BB(:)
    real(dp), allocatable :: global_R(:)
    integer,allocatable :: sparsity_col(:),reduced_col(:)
    integer,allocatable :: sparsity_row(:),reduced_row(:)


    real(dp) :: Ccun
    real(dp) :: dt
    logical :: carryon
    character(len=60) :: sub_name
    
    integer :: nonzeros_unreduced
    integer :: SOLVER_FLAG

    ! #############################################################################

!   set_diagnostics = .true.
    sub_name = 'solve_particles'
    call enter_exit(sub_name,1)


!!! assume a generic 9-generation acinus geometry, based on Haefeli-Bleuer & Weibel, 1988
!!! radius acinar ducts at TLC for 9 generations (HAEFELI-BLEUER, 1988)
    part_param%RacTLC = (/0.25_dp,0.25_dp,0.245_dp,0.2_dp,0.19_dp,0.18_dp,0.17_dp,0.155_dp,0.145_dp,0.125_dp/)

!!! length acinar ducts for terminal bronchiole 9 generations (HAEFELI-BLEUER, 1988)
    part_param%LacTLC = (/0.8_dp,1.33_dp,1.12_dp,0.93_dp,0.83_dp,0.7_dp,0.7_dp,0.7_dp,0.67_dp,0.75_dp/)

!!! acinar volume at TLC for 9 generations (estimated from HAEFELI-BLEUER, 1988)
    part_param%VacTLC = (/0.87_dp,1.70_dp,2.49_dp,4.81_dp,7.56_dp,14.07_dp,24.60_dp,43.93_dp,86.86_dp/)

!!! Cunningham slip correction factor 
    Ccun = 1.0_dp + 2.0_dp * part_param%lambda/part_param%pdia &
         * (1.257_dp+0.4_dp*exp(-0.55_dp*part_param%pdia/part_param%lambda))
    part_param%diffu = part_param%kBoltz*part_param%Temperature*Ccun/&
            (3.0_dp*pi*part_param%mu*part_param%pdia) ! diffusion constant [mm^1/s]

    ! allocatable array to store the current solution 
    if(.not.allocated(solution)) allocate(solution(num_nodes))
!   if(.not.allocated(M_C_M)) allocate(M_C_M(num_nodes))

    if(.not.allocated(part_acinus_field))then
       allocate(part_acinus_field(20,num_units))
       part_acinus_field(:,:) = 0.0_dp
    endif

    inlet_flow = elem_field(ne_flow,1) ! flow at entry element
    theta = 2.0_dp/3.0_dp

! --------------------------------------
     dt = part_param%dt_gm
! --------------------------------------

    ! get the sparsity arrays for the reduced system. uses compressed row format.
    call reduce_gasmix(MatrixSize,nonzeros,noffset_entry,noffset_row,inspiration)
    
    time = time_start ! initialise the time
    
    carryon = .true. ! logical for whether solution continues
    
    ! main time-stepping loop:  time-stepping continues while 'carryon' is true
    do while (carryon) !
       
       call cpu_time(cpu_dt_start)
       time = time + dt ! increment time

!      if(coupled)then ! update a single ventilation step
!         call evaluate_vent_step(dt)
!         inlet_flow = elem_field(ne_flow,1)
!         call calc_mass(nj_conc1,nu_conc1,temp_mass)
!      endif
       
       if(abs(inlet_flow).gt.loose_tol)then ! i.e. not for breath hold
          if(inlet_flow.gt.0.0_dp) &
           node_field(nj_conc1,1) = part_param%inlet_concentration(1) ! set the inlet concentration

          ! HK commented
          call airway_mesh_deform(dt,part_param%initial_volume,part_param%coupled,'particles') ! change model size by dV

          call general_track(dt,.true.)
          call particle_velocity(dt)
       endif

       ! assemble the element matrices. Element matrix calculation can be done directly 
       ! (based on assumption of interpolation functions) or using Gaussian interpolation.
!      call assemble_gasmix(diffusion_coeff,nonzeros_unreduced,problem_type) ! also for particles
       call assemble_particles(part_param%diffusion_coeff,nonzeros_unreduced) ! also for particles
       
       ! initialise the values in the solution matrices
       global_AA(1:nonzeros) = 0.0_dp ! equivalent to M in Tawhai thesis
       global_BB(1:num_nodes) = 0.0_dp ! equivalent to K in Tawhai thesis
       global_R(1:num_nodes) = 0.0_dp

       ! Assemble the reduced system of matrices
       do np=1,num_nodes ! Loop over rows of unreduced system
          nrow = np ! conveniently true for the way we set up our models
          ! different boundary conditions are applied during inspiration and 
          ! expiration: Dirichlet at model entry during inspiration (concentration
          ! = inlet_concentration), and Neumann at model entry during expiration (dcdt = 0)
          if(.not.inspiration)then
             BB = global_R(nrow)         !get reduced R.H.S.vector
             do nentry = sparsity_row(nrow),sparsity_row(nrow+1)-1  !each row entry
                ncol = sparsity_col(nentry)
                BB = BB-global_K(nentry)*node_field(nj_conc1,ncol)*dt ! -K*c^(n)*dt
                AA = global_M(nentry) + dt*theta*global_K(nentry) ! M+K*dt*theta
                global_AA(nentry) = global_AA(nentry) + AA
             enddo
             global_BB(nrow) =  global_BB(nrow) + BB
          elseif(inspiration.and.np.ne.1)then !not first row
             BB = global_R(nrow)         !get reduced R.H.S.vector
             do nentry = sparsity_row(nrow),sparsity_row(nrow+1)-1  !each row entry
                ncol = sparsity_col(nentry)
                BB = BB-global_K(nentry)*node_field(nj_conc1,ncol)*dt
                AA = global_M(nentry) + dt*theta*global_K(nentry) !M+K*dt*theta
                if(ncol.ne.1)then ! not first column
                   global_AA(nentry-noffset_entry) = &
                        global_AA(nentry-noffset_entry) + AA
                endif
             enddo
             global_BB(nrow-noffset_row) = &
                  global_BB(nrow-noffset_row) + BB
          endif
       enddo
       
       solution(1:num_nodes) = zero_tol
       
       ! Call a solver to solve the system of reduced equations. 
       ! Here we use an iterative solver (GMRES == Generalised Minimal 
       ! RESidual method). The solver requires the solution matrices to 
       ! be represented in compressed row format.
       call pmgmres_ilu_cr (MatrixSize,NonZeros,reduced_row,&
            reduced_col,global_AA,solution,global_BB,&
            part_param%out_itr_max,part_param%inr_itr_max,&
            part_param%solve_tolerance,part_param%solve_tolerance,SOLVER_FLAG)
       
       ! transfer the solver solution (in 'Solution') to the node field array          
       do np = 1,num_nodes
          if(.not.inspiration.or.(inspiration.and.np.gt.1))then
             if(inspiration)then
                nrow_BB=np-1
             else
                nrow_BB=np
             endif
             node_field(nj_conc1,np) = node_field(nj_conc1,np) &
                  + Solution(nrow_BB) !c^(n+1)=c^(n)+dc
             node_field(nj_conc1,np) = max(0.0_dp,node_field(nj_conc1,np))
             node_field(nj_conc1,np) = min(1.0_dp,node_field(nj_conc1,np))
          endif
       enddo

       call volume_of_mesh(current_volume,volume_tree)

       if(inlet_flow.lt.0.0_dp)then
          call particle_deposition(current_volume,dt,.false.)
       else
          call particle_deposition(current_volume,dt,.true.)
       endif
       
       ! estimate the volume and mass errors
       call calc_mass_particles(nj_conc1,nu_conc1,current_mass,mass_deposit)
       
!       if(coupled)then ! update a single ventilation step
!          call evaluate_vent_step(dt)
!          inlet_flow = elem_field(ne_flow,1)
!          call calc_mass(nj_conc1,nu_conc1,temp_mass)
!       endif
       
       ideal_mass = ideal_mass + inlet_flow*dt*node_field(nj_conc1,1)
       volume_error = 1.0e+2_dp*(current_volume - (part_param%initial_volume +&
            total_volume_change))/(part_param%initial_volume + total_volume_change)
       if(ideal_mass.gt.0.0_dp)then
          mass_error = 1.0e+2_dp*(current_mass + mass_deposit - ideal_mass)/ideal_mass
       else
          mass_error = 0.0_dp
       endif
       call cpu_time(cpu_dt_end)

       write(*,'(f7.3,6(f8.2),f9.2,f11.5)') &
            time,inlet_flow/1.0e+6_dp,total_volume_change/1.0e+6_dp,&
            current_volume/1.0e+6_dp,current_mass/1.0e+6_dp,&
            mass_deposit/1.0e+6_dp,volume_error, &
            mass_error,node_field(nj_conc1,1)
       
       err = mass_error
       
       
!      if(last_breath.and.coupled)then
       if(last_breath)then
          if(elem_field(ne_flow,1).lt.0.0_dp) then
             carryon = .true.
          else
             carryon = .false.
          endif
       else
          if(time_end - time .le. zero_tol)then
             carryon = .false.
          else
             carryon = .true.
          endif
       endif
    enddo !carryon
    
    call enter_exit(sub_name,2)
  end subroutine solve_particles
  
    

!!!###################################################################################

  subroutine scale_flow_field(inlet_flow)

!   use arrays,only: dp,elem_field,num_elems,num_units,unit_field

    use arrays
    use other_consts
    use indices

!   use indices,only: ne_dvdt,nu_flow0

    real(dp),intent(in) :: inlet_flow
    real(dp) :: ratio

    if(abs(elem_field(ne_flow,1)).gt.zero_tol)then
       ratio = inlet_flow/elem_field(ne_flow,1)

       unit_field(nu_flow0,1:num_units) = unit_field(nu_flow0,1:num_units)*ratio
       elem_field(ne_flow,1:num_elems) = elem_field(ne_flow,1:num_elems)*ratio
       elem_field(ne_dvdt,1:num_elems) = elem_field(ne_dvdt,1:num_elems)*ratio

    else
       write(*,'('' Cannot scale to zero flow'')')
    endif

  end subroutine scale_flow_field

!###################################################################################


  subroutine general_track(dt,update)  
    ! the argument 'update' is a tag to decide whether the concentrations have to updated or not
    ! this should be false in case of mass calculation which is used in the "update_unit_mass" subroutine

!   use arrays,only: dp,node_field,elem_field,num_elems,num_units,unit_field, elem_nodes, elem_cnct,&
!                    units, elems_at_node
    use arrays

!   use indices, only: nj_conc1, nj_conc2, nj_conc3
    use indices

    use geometry
    use other_consts
    use diagnostics, only: enter_exit
    use gasmix
    
    implicit none
  
    type(particle_parameters) :: part_param
    real(dp), intent(in) :: dt
    logical, intent(in) :: update

    integer,parameter :: n_gases = 3
    integer :: c_kount,i,i_elem,j,j1,jj2,kount,ne,ne_stem,nj_g(3),np,np1,np2,np1_parent,np2_parent, &
         num_list_total,nunit,nu_g(3),o_elem,sum_inout,sum_dir
    integer,parameter :: elem_f_to_node = 1, f_sign = 2
    integer,allocatable :: TMAT(:,:,:),elem_list_total(:)
    real(dp) :: c_mass(3),concens(3),concentration(3),flow_fraction, fl_temp, &
         inlet_conc(3),local_xi,mass_fraction,mean_c(3),parent_fraction, &
         time_through_element,total_time, unit_mass(3), &
         ratio_mass,total_mass, max_concentration
!!! currently unused but might be again
!!!  ideal_mass_total,mass_at_max,mass_below_max, total_mass,ratio_mass,
    real(dp), allocatable :: concs_at_node(:,:), concent(:,:)
    logical :: go_on, cont, neg_fl
    character(len=60) :: sub_name

    sub_name = 'general_track'
    call enter_exit(sub_name,1)

    nj_g(1) = nj_conc1
    nj_g(2) = nj_conc2
    nj_g(3) = nj_conc3
    nu_g(1) = nu_conc1
    nu_g(2) = nu_conc2
    nu_g(3) = nu_conc3
    inlet_conc(:) = part_param%inlet_concentration(:)

    if (.not.allocated(TMAT)) allocate(TMAT(num_nodes,2,3))
    
    !second index: entering/exiting (1/2) flow towards/from a node 
    ! and positive/negative (1/2) flow at each element (positive: downward,
    ! negative:upward), Third index: the elements at the node
    
!!! If the flow is negative in an element, the negative flow always come from its children
!!! however, if the flow is positive the flow might be coming from its parents or its parents' other children
    
    call track_mat_coupled(elem_f_to_node,f_sign,TMAT) !calling the matrix containing flow directions and nodes in/out data
!!! TMAT(np,2,ne) = 1 when flow in element ne is distal (to small airways)
!!! TMAT(np,2,ne) = 2 when flow in element ne is proximal (to large airways)
!!! TMAT(np,1,ne) = 1 when element ne brings flow to node np
!!! TMAT(np,1,ne) = 2 when element ne takes flow from node np
    
    allocate(concent(num_nodes,n_gases))
    concent = 0.0_dp

    do nunit = 1,num_units
       ne = units(nunit)
       if(elem_cnct(1,0,ne).eq.0)then ! a terminal
          np = elem_nodes(2,ne)
          if(elem_field(ne_flow,ne).lt.0.0_dp)then
             do i = 1,n_gases
                unit_mass(i) = unit_field(nu_vol,nunit)*unit_field(nu_g(i),nunit)
                c_mass(i) = elem_field(ne_flow,ne)*node_field(nj_g(i),np)*dt
                unit_mass(i) = unit_mass(i) + c_mass(i)
                !             unit_field(nu_g(i),nunit) = unit_mass(i)/unit_field(nu_vol,nunit)
                node_field(nj_g(i),np) = unit_mass(i)/unit_field(nu_vol,nunit)
                forall (i=1:n_gases) concent(np,i) = node_field(nj_g(i),np)
             enddo
          endif
       endif ! terminal
    enddo

    do np = 1,num_nodes
       allocate(concs_at_node(elems_at_node(np,0),n_gases))
       concs_at_node = 0.0_dp

!!!    check that the flows at the node make sense. Can't have flow in 
!!!    all elements directed towards the node, and can't have all flows directed out from the node.
       sum_inout = 0
       sum_dir = 0
       do J = 1,elems_at_node(np,0)
          sum_inout = sum_inout + TMAT(np,elem_f_to_node,J) !sum of flow directions: proximal (2) and distal (1) 
          sum_dir = sum_dir + TMAT(np,f_sign,J)     !sum of flow directions w.r.t. node np
       enddo
       go_on = .true.
       if ((sum_inout == 6) .and. (sum_dir == 4)) go_on = .false. !impossible
       if ((sum_inout == 3) .and. (sum_dir == 5)) go_on = .false. !trapping
       if ((sum_inout == 4) .and. (sum_dir == 3)) go_on = .false. !impossible
       if ((sum_inout == 2) .and. (sum_dir == 3)) go_on = .false. !trapping


       if (go_on) then !meets criteria for sensible flow
          concs_at_node = 0.0_dp
          do J = 1,elems_at_node(np,0)
             ne = elems_at_node(np,J)
             if (TMAT(np,1,J) == 1) then !flow entering the node so contributing to the concentration
                if (TMAT(np,2,J) == 1) then !positive flow
                   c_mass = 0.0_dp  !in case of any appended units!
                   time_through_element = abs(elem_field(ne_vol,ne)/elem_field(ne_flow,ne))
                   total_time = time_through_element
                   np1 = elem_nodes(1,ne)
                   np2 = elem_nodes(2,ne)
                   if (total_time.ge.dt) then !location is within this element
                      local_xi = dt/time_through_element
                      do i = 1,n_gases
                         concentration(i) = ((1.0_dp-local_xi)*node_field(nj_g(i),np2)) + &
                              (local_xi*node_field(nj_g(i),np1))
                         mean_c(i) = 0.5_dp*(concentration(i) + node_field(nj_g(i),np2))
                         c_mass(i) = c_mass(i) + mean_c(i) &
                              *abs(elem_field(ne_vol,elems_at_node(np2,1))) !!! for acinus
                      enddo
                   else
                      do i = 1,n_gases
                         mean_c(i) = 0.5_dp*(node_field(nj_g(i),np1)+node_field(nj_g(i),np2))
                         c_mass(i) = c_mass(i) + mean_c(i)*abs(elem_field(ne_vol,elems_at_node(np2,1)))
                      enddo
                      cont = .true.
                      JJ2 = elems_at_node(np,J) ! same value as ne?
                      concentration = 0.0_dp
                      parent_fraction = 1.0_dp
                      J1 = JJ2
                      do while (cont)
                         JJ2 = elem_cnct(-1,1,JJ2) !parent element
                         if (JJ2.gt.0) then !parent is not the trachea
                            np1_parent = elem_nodes(1,JJ2)
                            np2_parent = elem_nodes(2,JJ2)
                            ! Because different parents and children are supplying the node
                            ! they have their own flow fractions. Parents are getting flow from
                            ! their parents for which they might get flow from their children.
                            ! For the children fractions are coming from their nodes but for the
                            ! parents the fraction is from their parents and their parent's children
                            ! going backward to the trachea.
                            ! summing up concentrations for child suppliers. 
                            if (elems_at_node(np2_parent,0) == 3) then
                               ! check whether all of the flow comes from the parent, or some from a sibling
                               neg_fl = .false.
                               do kount = 2,3
                                  ! one of these is the node itself which has positive flow so
                                  ! automatically the other one is selected
                                  if (TMAT(np2_parent,2,kount) == 2) then ! negative flow in sibling
                                     i_elem = elems_at_node(np2_parent,kount) ! sibling element number
                                     neg_fl = .true.
                                     do i = 1,n_gases
                                        call negative_flow(dt,nj_g(i),np2_parent,kount,TMAT,concens(i))
                                     enddo
                                     fl_temp = abs(elem_field(ne_flow,i_elem))+ &
                                          abs(elem_field(ne_flow,JJ2)) ! sum parent and sibling flows
!                                  else
                                     o_elem = elems_at_node(np2_parent,kount)
                                  endif !TMAT
                               enddo !kount
                               if (neg_fl) then
                                  if (TMAT(np2_parent,2,1) == 2) then ! negative in parent; all comes from sibling
                                     flow_fraction = 1.0_dp ! MHT: haven't checked that this part is correct
                                  else
!                                     fl_temp = abs(elem_field(ne_flow,o_elem)+ &
!                                          elem_field(ne_flow,JJ2))
                                     flow_fraction = abs(elem_field(ne_flow,o_elem)/fl_temp)
                                  endif
                                  !calculate the total fraction up to the supplying parent
                                  parent_fraction = parent_fraction * abs(1.0_dp-flow_fraction)
                                  do i = 1,n_gases
                                     concentration(i) = concentration(i) + (flow_fraction * concens(i))
                                     mean_c(i) = 0.5d0*(concens(i) + node_field(nj_g(i),np2_parent))
                                     c_mass(i) = c_mass(i) + abs(mean_c(i) * elem_field(ne_vol,i_elem)* &
                                          flow_fraction)
                                  enddo
                                  ! the mass the negative flow chidren are adding are just depending
                                  ! on their own flow fraction (not others)
                               endif
                            endif !elems_at_...
                            
                            ! adding the concentration for the tracked parent supplier with respect
                            ! to its fraction down to the node to the concentration of the child
                            ! suppliers calculated above
                            mass_Fraction = elem_field(ne_flow,J1)/elem_field(ne_flow,JJ2)
                            !if (JJ2.gt.0) then !parent is not the trachea
                            time_through_element = abs(elem_field(ne_vol,JJ2)/&
                                 (elem_field(ne_flow,JJ2)))
                            total_time = total_time + time_through_element
                            if (TMAT(np2_parent,2,1) == 2) then
                               forall (i = 1:n_gases) concentration(i) = concentration(i) + (parent_fraction * &
                                    node_field(nj_g(i),np2_parent)) 
                               cont = .false.
                            else if (total_time.ge.dt) then ! Material is within this element
                               local_xi = (total_time-dt)/time_through_element
                               do i = 1,n_gases
                                  concentration(i) = concentration(i) + (parent_fraction * &
                                       (((1.d0-local_xi)*node_field(nj_g(i),np1_parent)) + &
                                       (local_xi*node_field(nj_g(i),np2_parent))))
                                  mean_c(i) = 0.5_dp*(concentration(i) + node_field(nj_g(i),np2_parent))
                                  c_mass(i) = c_mass(i) + mean_c(i) * abs(elem_field(ne_vol,JJ2) * &
                                       (1.d0-local_xi) * mass_fraction) !!!acinus
                               enddo
                               cont = .false.
                            else
                               do i = 1,n_gases
                                  mean_c(i) = 0.5_dp*(node_field(nj_g(i),np1_parent) + node_field(nj_g(i),np2_parent))
                                  c_mass(i) = c_mass(i) + mean_c(i) * abs(elem_field(ne_vol,JJ2) * mass_fraction)
                               enddo
                            endif !total_time.ge.dt
                         else
                            mass_fraction = elem_field(ne_flow,J1)/elem_field(ne_flow,1)
                            forall (i = 1:n_gases) concentration(i) = concentration(i) + (parent_fraction * inlet_conc(i))

!!! if the fluid is being advected from outside the model, then the 'interpolated mass' should be such that the
!!! concentration is the same as the inspired concentration (inlet concentration), i.e. M = effective_volume*c_inlet
                            !                            interp_mass = interp_mass + (parent_fraction * node_field(nj_mass,1)*2.0_dp)
                            forall (i = 1:n_gases) c_mass(i) = c_mass(i) + concentration(i)*(dt-total_time)*mass_fraction
                            cont = .false.
                         endif !JJ2.gt.0
                      enddo !cont
                   endif !total_time.ge.dt
                   forall (i = 1:n_gases) concs_at_node(J,i) = concentration(i)
                endif !TMAT(np,2,J) == 1
                if (TMAT(np,2,J) == 2) then !negative flow
                   do i = 1,n_gases
                      call negative_flow(dt,nj_g(i),np,J,TMAT,concens(i))
                      concs_at_node(J,i) = concens(i)
                      ! estimate the expired mass:
                      c_mass(i) = 0.5_dp*(node_field(nj_g(i),np)+concs_at_node(J,i))* &
                           elem_field(ne_flow,elems_at_node(np,1))*dt
                   enddo
                endif !TMAT(np,2,J) == 2
             endif !if (TMAT(np,1,J) == 1)
          enddo !J=1...
       else
          write(*,'('' WARNING: ignored node:'',i8)') np
          if ((sum_inout == 6) .and. (sum_dir == 4)) write(*,'('' Impossible'')')
          if ((sum_inout == 3) .and. (sum_dir == 5)) write(*,'('' Trapping'')')
          if ((sum_inout == 4) .and. (sum_dir == 3)) write(*,'('' Impossible'')')
          if ((sum_inout == 2) .and. (sum_dir == 3)) write(*,'('' Trapping'')')
          do J = 1,elems_at_node(np,0)
             ne = elems_at_node(np,j)
             write(*,'('' Flow from element'',i8,'' ='',d12.4,'' mm^3/s'')') ne,elem_field(ne_flow,ne)
          enddo
          read(*,*)
       endif !if (go_on)

!!! update the temporary array for nodal concentrations (concent)
       c_kount = 0
       fl_temp = 0.0_dp
       do J = 1,elems_at_node(np,0)
          ne = elems_at_node(np,J)
          if (TMAT(np,1,J) == 1) then
             fl_temp = fl_temp + abs(elem_field(ne_flow,ne))!-elem_field(ne_dvdt,elems_at_node(np,J)))
             c_kount = c_kount + 1  !check the number of incoming flows to the node
             forall (i = 1:n_gases) concent(np,i) = concs_at_node(J,i) !if there's only one supplier concentration remains the same
           else
             o_elem = elems_at_node(np,J) !the outgoing element
          endif !conc_at_...
       enddo !J=1...
       if (c_kount == 2) then !two elements supplying the node so check the flow fractions
         concent(np,1:n_gases) = 0.0_dp
         do J = 1,elems_at_node(np,0)
           if (elems_at_node(np,J).ne.o_elem) then
             flow_fraction = abs((elem_field(ne_flow,elems_at_node(np,J)))/fl_temp)
             forall(i = 1:n_gases) concent(np,i)=concent(np,i)+(flow_fraction*concs_at_node(J,i))
           endif !elems_at...
         enddo !J
       endif ! c_kount

       do i = 1,n_gases
          if (concent(np,i) .le. 1.0e-10_dp) concent(np,i) = 0.0_dp
       enddo
       
       if ((np.eq.1).and.(TMAT(np,2,1) == 1))then
          ! an entrance node inspiring
          forall (i = 1:n_gases) concent(np,i) = inlet_conc(i)
       endif
       if ((np/=1).and.(elems_at_node(np,0)==1).and.(TMAT(np,1,1) == 2)) then
          ! a terminal node expiring
          ne = elems_at_node(np,1)
          if(inlist(ne,units))then ! has a lumped unit attached
             nunit = where_inlist(ne,units)
!             forall (i=1:2) node_field(nj_g(i),np) = unit_field(nu_g(i),nunit)
             forall (i = 1:n_gases) c_mass(i) = elem_field(ne_flow,ne)*node_field(nj_g(i),np)*dt
          endif
          forall (i = 1:n_gases) concent(np,i) = node_field(nj_g(i),np)
       endif
       
       deallocate(concs_at_node)

!!! add the accumulated mass to lumped units
       ne = elems_at_node(np,1) ! the first element that the node is in
       if ((np/=1).and.(elems_at_node(np,0) == 1))then ! a terminal node
          if(inlist(ne,units))then ! is terminal with a lumped unit attached
             nunit = where_inlist(ne,units)
             do i = 1,n_gases
                unit_mass(i) = unit_field(nu_vol,nunit)*unit_field(nu_g(i),nunit)
                unit_mass(i) = unit_mass(i) + c_mass(i)
                unit_field(nu_g(i),nunit) = unit_mass(i)/unit_field(nu_vol,nunit)
             enddo
          endif
       else ! is not terminal
          if(inlist(ne,units))then ! this element supplies a branching acinus within a unit
             ! record the cumulative mass entering the acinus + initial mass (dummy storage)
             elem_field(ne_resist,ne) = elem_field(ne_mass,ne) + c_mass(1)
          endif
       endif

    enddo !np

    max_concentration = max(part_param%initial_concentration,part_param%inlet_concentration(1)) ! for primary gas

    if (update) then
       forall(i = 1:n_gases) node_field(nj_g(i),1:num_nodes)= concent(1:num_nodes,i)

!!! adjust the concentrations to conserve mass. Only done for nodes where the concentration
!!! is less than the maximum (set in 'initial_gasmix')
!   ideal_mass_total = elem_field(ne_mass,1) + elem_field(ne_flow,1)*dt*node_field(nj_conc1,1)
       call calc_mass(nj_conc1,nu_conc1,total_mass)
       allocate(elem_list_total(num_elems))
       do nunit = 1,num_units
          ne_stem = units(nunit)
          num_list_total = 0
          if(elem_cnct(1,0,ne_stem).ne.0)then
             elem_list_total = 0
!!! get a list of all elements that are within an elastic unit
             call group_elem_by_parent(ne_stem,elem_list_total)
             num_list_total = count(elem_list_total.ne.0)
          endif

          !! the 'ideal' mass for the unit is elem_field(ne_resist,ne_stem)
          !! the 'current' mass for the unit is elem_field(ne_mass,ne_stem)
          if(elem_field(ne_mass,ne_stem).gt.zero_tol)then
             ratio_mass = elem_field(ne_resist,ne_stem)/elem_field(ne_mass,ne_stem)
          else
             ratio_mass = 1.0_dp
          endif
          do i = 1,num_list_total
            ne = elem_list_total(i)
            np2 = elem_nodes(2,ne)
            if(node_field(nj_conc1,np2).lt.max_concentration) then
                node_field(nj_conc1,np2) = node_field(nj_conc1,np2)*ratio_mass
            endif
          enddo
       enddo
       deallocate(elem_list_total)
    endif

    deallocate(concent)
    deallocate(TMAT)

    call enter_exit(sub_name,2)
    
  end subroutine general_track
  

!!!############################################################################################


  subroutine track_mat_coupled(elem_f_to_node,f_sign,TMAT)

!   use arrays, only: dp, elems_at_node, elem_field
    use arrays
    use diagnostics, only: enter_exit

    implicit none
  
    integer,intent(in) :: elem_f_to_node,f_sign
    integer :: TMAT(:,:,:)
    integer :: I, J
    real(dp) :: flo
    character(len=60) :: sub_name
    
    sub_name = 'track_mat_coupled'
    call enter_exit(sub_name,1)
    
    ! put these two below also in the array.f90 as global parameters
    TMAT = 0
    do I = 1,num_nodes
       do J = 1,elems_at_node(I,0)
          flo = elem_field(ne_flow,elems_at_node(I,J))
          if (flo .lt. 0.0_dp) then
             if (I.eq.1) then
                TMAT(I,elem_f_to_node,J) = 1
                TMAT(I,f_sign,J) = 2
             elseif (J.eq.1) then !parent element
                TMAT(I,elem_f_to_node,J) = 2
                TMAT(I,f_sign,J) = 2
             else !child elements
                TMAT(I,elem_f_to_node,J) = 1
                TMAT(I,f_sign,J) = 2
             endif
          else
             if (I.eq.1) then
                TMAT(I,elem_f_to_node,J) = 2
                TMAT(I,f_sign,J) = 1
             elseif (J.eq.1) then !parent element
                TMAT(I,elem_f_to_node,J) = 1
                TMAT(I,f_sign,J) = 1
             else !child elements
                TMAT(I,elem_f_to_node,J) = 2
                TMAT(I,f_sign,J) = 1
             endif
          endif
       enddo
    enddo
    
    call enter_exit(sub_name,2)

  end subroutine track_mat_coupled

!!!############################################################################################

  subroutine negative_flow(dt,nj_g,I,J,TMAT,concs)

!   use arrays,only: dp, elem_cnct, elem_symmetry, elems_at_node, node_field, elem_nodes
    use arrays
    use diagnostics, only: enter_exit

    implicit none
  
    integer, intent(in) :: I, J, nj_g,TMAT(:,:,:)
    real(dp), intent(in) :: dt
    real(dp), intent(out) :: concs
    
    integer :: I1, I2, JJ3, kount, child, c_kount, IK
    integer, allocatable :: parents(:), children(:)
    real(dp) :: time_through_element, total_time, concentration, local_xi, fl_temp 
    real(dp) :: flow_fraction, flow_pa
    real(dp), allocatable :: parent_times(:), child_times(:)
    logical :: cont, in_cont
    character(len=60) :: sub_name
    
    sub_name = 'negative_flow'
    call enter_exit(sub_name,1)
    
    ! Ok now we are creating a matrix containing the final nodes from which
    ! the materials are reaching the current node. Then we will calculate
    ! the summation of the concentrations coming from these finals nodes by
    ! going back from them towards the current node and taking into account
    ! the flow fractions where two adjacent elements are feeding their parent
    ! element with negative flow
    
    allocate(parents(1))
    allocate(parent_times(1))
    allocate(children(128))
    allocate(child_times(128))
    
    cont = .true.
    in_cont = .true.
    parents(1) = elem_nodes(2,elems_at_node(I,J)) !second node of the child element
    parent_times(1) = abs(elem_field(ne_vol,elems_at_node(I,J))/ &
         (elem_field(ne_flow,elems_at_node(I,J))))
    if (parent_times(1).gt.dt) cont = .false.
    do while (in_cont)
       in_cont = .false.
       c_kount = 0
       do kount = 1,size(parents)
          cont = .false.
          do JJ3 = 2,elems_at_node(parents(kount),0) !find the children for parents (parents are nodes here)
             child = elems_at_node(parents(kount),JJ3) ! Child element for that node
             time_through_element = abs(elem_field(ne_vol,child)/&
                  (elem_field(ne_flow,child)))
             total_time = time_through_element + parent_times(kount) !calculate time in each child to see if it's less than dt so the material passes it
             if ((elems_at_node(parents(kount),0).gt.1) .or. (I == 1)) then !parent node is not a terminal node or it is at trachea
                if (TMAT(parents(kount),2,JJ3) == 2) then !the child has negative flow
                   if (parent_times(kount) .lt. dt) then !material is not in this elemen
                      cont = .true. !if the parents are not terminals
                      ! ,their children's flow are negative and the material is not within them
                      ! then we have to go furthur down
                      in_cont = .true. ! still need to check elems thru which material passes; might be some -ve flow in other subtended branches
                      c_kount = c_kount + 1
                      children(c_kount) = elem_nodes(2,child)
                      child_times(c_kount) = total_time !saving the time for the element to be used in calculating the concentrations at its nodes later
                   endif ! time through
                endif !TMAT
             endif ! long if condition above (elem_at_ .....)
          enddo !JJ3
          if (.not.cont) then !either the child elements have positive flow or their parent is a terminal
             !the parent should be preserved
             c_kount = c_kount + 1
             if(c_kount.gt.128)then
                write(*,*) 'Increase size of children array in general_track, or decrease dt'
                read(*,*)
             endif
             children(c_kount) = parents(kount)
             child_times(c_kount) = parent_times(kount)
          endif ! not cont
       enddo !kount
       deallocate(parents)
       allocate(parents(c_kount))
       deallocate(parent_times)
       allocate(parent_times(c_kount))
       parents = children(1:c_kount)
       parent_times = child_times(1:c_kount)
    enddo !cont
    deallocate(children)
    deallocate(child_times)

    concs = 0.0_dp
    do kount = 1,size(parents)
       concentration = 0.0_dp
       I1 = elem_nodes(1,elems_at_node(parents(kount),1)) !first node of the element containing the node
       I2 = elem_nodes(2,elems_at_node(parents(kount),1)) !second node of the element containing the node (actually the node itself!!)
       time_through_element = abs(elem_field(ne_vol,elems_at_node(parents(kount),1))/&
            (elem_field(ne_flow,elems_at_node(parents(kount),1))))
       local_xi = abs((parent_times(kount)-dt))/time_through_element
       if ((parent_times(kount) .lt. dt) ) local_xi = 0.0_dp ! might only happen at the final elements(terminals)
       concentration =  ((1.0_dp-local_xi)*node_field(nj_g,I2)) + &
            (local_xi*(node_field(nj_g,I1)))
       JJ3 = elems_at_node(parents(kount),1) !element contatining the node
       do while (JJ3.gt.elems_at_node(I,J))
          c_kount = elems_at_node(I1,0)
          if (c_kount == 3) then ! the node has 3 elements
             if ((TMAT(I1,1,2) == 1) .and. (TMAT(I1,1,3) == 1)) then !both child elements are supplying the node
                
                fl_temp = 0.0_dp
                do IK = 2,c_kount
                   fl_temp = fl_temp + abs(elem_field(ne_flow,elems_at_node(I1,IK)))
                enddo
                flow_fraction = abs((elem_field(ne_flow,JJ3))/fl_temp)
                !if the other children is also supplying the same node it should be already in the list (parents) and taken care of for its flow fraction here so no need to calculate for both. 
                concentration = flow_fraction * concentration
             endif
          elseif (c_kount == 2) then
             flow_pa = elem_field(ne_flow,elem_cnct(-1,1,JJ3)) - elem_field(ne_dvdt,elem_cnct(-1,1,JJ3))
             flow_fraction = abs((elem_field(ne_flow,JJ3) * elem_symmetry(JJ3))/flow_pa)
             concentration = flow_fraction * concentration
          endif
          JJ3 = elem_cnct(-1,1,JJ3) !set the element to its parent
          I1 = elem_nodes(1,JJ3)
       enddo ! JJ3
       concs = concs + concentration
    enddo !kount
    
    deallocate(parents)
    deallocate(parent_times)

    call enter_exit(sub_name,2)
    
  end subroutine negative_flow
  
!!!############################################################################################

  subroutine track_mat(TMAT)

    use arrays
    use diagnostics, only: enter_exit

    implicit none
  
    integer :: TMAT(:,:,:)
    integer :: ne,noelem,np
    integer,parameter :: elem_f_to_node = 1, f_sign = 2
    real(dp) :: flo
    character(len=60) :: sub_name

    sub_name = 'track_mat'
    call enter_exit(sub_name,1)

    TMAT = 0
    do np = 1,num_nodes
       do noelem = 1,elems_at_node(np,0)
          ne = elems_at_node(np,noelem)
          flo = elem_field(ne_flow,ne)
          if (flo .lt. 0.0_dp) then
!!! flow is directed proximally
             TMAT(np,f_sign,noelem) = 2  !labels flow as in direction 2 (proximal)
             if (np.eq.1) then
                TMAT(np,elem_f_to_node,noelem) = 1 !brings flow to the node
             else if (noelem.eq.1) then !parent element
                TMAT(np,elem_f_to_node,noelem) = 2 !takes flow away from the node
             else !child elements
                TMAT(np,elem_f_to_node,noelem) = 1 !brings flow to the node
             endif
          else
!!! flow is directed distally
             TMAT(np,f_sign,noelem) = 1  !labels flow as in direction 1 (distal)
             if (np.eq.1) then
                TMAT(np,elem_f_to_node,noelem) = 2 !takes from away from the node
             else if (noelem.eq.1) then !parent element
                TMAT(np,elem_f_to_node,noelem) = 1 !brings flow to the node
             else !child elements
                TMAT(np,elem_f_to_node,noelem) = 2 !takes flow away from the node
             endif
          endif
       enddo
    enddo

    call enter_exit(sub_name,2)
    
  end subroutine track_mat

  
!!!############################################################################################

  function effective_volume_at_node(np)

    use arrays, only: dp, elems_at_node, elem_nodes, elem_field, elem_symmetry
    integer,intent(in) :: np

    integer :: i,ne,np1,np2
    real(dp) :: effective_volume_at_node

    effective_volume_at_node = 0.0_dp
    do i = 1,elems_at_node(np,0)
       ne = elems_at_node(np,i)
       np1 = elem_nodes(1,ne)
       np2 = elem_nodes(2,ne)
       effective_volume_at_node = effective_volume_at_node + 0.5_dp*elem_field(ne_vol,ne)*elem_symmetry(ne)
    enddo

  end function effective_volume_at_node

!!!################################################################################
  
  subroutine particle_velocity(dt)
    ! Created by Falko Schmidt 05/2011
    ! calculates the particle velocity of spherical particles
    ! from the fluid flow field and the previous particle velocity due to drag
    
!   use arrays, only: dp, particle_parameters, elem_field,&
!                     node_xyz,elem_nodes,node_xyz,&
!                     num_elems,num_nodes
    use arrays
    use mesh_utilities,only: angle_btwn_vectors
    use diagnostics, only: enter_exit

    implicit none
  
    type(particle_parameters) :: part_param
    real(dp),intent(in) :: dt
    integer :: j,ne,np,np2
    real(dp) :: ref_flow,rel_flow,alpha,abbr,vs, Ccun
    real(dp) :: vector1(3),vector2(3),vseff
    character(len=60) :: sub_name

    sub_name = 'particle_velocity'
    call enter_exit(sub_name,1)

!!! abbreviation for velocity equation
    Ccun = 1.0_dp + 2.0_dp * part_param%lambda/part_param%pdia &
         * (1.257_dp+0.4_dp*exp(-0.55_dp*part_param%pdia/part_param%lambda))

    abbr = -(18.0_dp*part_param%mu)/&
            (part_param%prho*Ccun*dt*&
                    part_param%pdia**2.0_dp)

!!! settlement velocity
    vs = (part_param%prho)*part_param%gravityz*part_param%pdia**2.0_dp&
                                             *Ccun/(18.0_dp*part_param%mu)
    
    do ne = 1,num_elems  !loop over all 1D elements
       np = elem_nodes(1,ne)  ! proximal node number
       np2 = elem_nodes(2,ne) ! distal node number
!!! calculate angle between branch axis and gravitational direction (alpha)
       forall(j=1:3) vector1(j) = node_xyz(j,np2)-node_xyz(j,np)
       vector2 = 0.0_dp
       vector2(3) = node_xyz(3,np)+1.0_dp
       alpha = angle_btwn_vectors(vector1,vector2)

!!! settlement velocity w.r.t. flow direction; set equal to zero to deactivate gravity effect        
       vseff = cos(alpha)*pi*elem_field(ne_radius,ne)**2*vs

!!! note following:
!!! ne_part_vel was stored as nej_flow in cmiss
!!! ne_part_conc was stored as nj_source in cmiss
!!! ne_part_var was stored as nj_loss in cmiss

!!! CALCULATE PARTICLE VELOCITY (stored in elem_field(ne_part_vel,ne))
!!! note that the influence of the Reynolds number is neglected because
!!! for particles <10 micron they are nearly as fast as the fluid
       ref_flow = elem_field(ne_flow,ne) ! the actual flow (instead of normalized) 
       rel_flow = elem_field(ne_part_vel,ne)-ref_flow+vseff  ! relative flow particle / fluid
       elem_field(ne_part_vel,ne) = ref_flow-vseff
       if(abbr.gt.-25.0_dp) elem_field(ne_part_vel,ne)=elem_field(ne_part_vel,ne) &
                                                          + rel_flow*exp(abbr)

!!! calculate spatial derivative of particle velocity
!       elem_field(ne_dpart_vel,ne) = 0.0_dp !(XP(nk,nv2,nej_flow,np2)/XP(1,nv2,nj_radius,np2)**2 &
!            -XP(nk,nv,nej_flow,np)/XP(1,nv2,nj_radius,np)**2 )/elem_field(ne_length,ne)/PI    ! spatial derivative of particle velocity
!       XP(2,nv2,nej_flow,np2) = (XP(nk,nv,nej_flow,np)/XP(1,nv2,nj_radius,np)**2 &
!            -XP(nk,nv2,nej_flow,np2)/XP(1,nv2,nj_radius,np2)**2)/elem_field(ne_length,ne)/PI    
    enddo

    call enter_exit(sub_name,2)
    
  end subroutine particle_velocity

!!!#######################################################################################
  
  subroutine particle_deposition(current_volume,dt,inspiration)

!!!    Calculates amount of deposited particles in each element,
!!!    reduces respective concentration field 'nj_source' and 
!!!    stores the accumulated deposition quantity in field 'nj_loss'
!!!    ADDIONALLY WRITES current solution vector YP to source
!!!    field 'nj_source'

!   use arrays, only: dp, particle_parameters, node_field, elem_field, &
!        part_acinus_field,elem_cnct,elem_nodes, elems_at_node, node_xyz,& 
!        num_nodes
    use arrays
    use indices
    use other_consts
    use mesh_utilities,only: angle_btwn_points
    use diagnostics, only: enter_exit
    
    implicit none
  
    type(particle_parameters) :: part_param
    real(dp),intent(in) :: current_volume,dt
    logical,intent(in) :: inspiration
    !     Local Variables
    integer :: gen,i,j,maxi,ne,ne0,noelem,np,npend,npstart, &
         npstep,np0,np1,np2,np3,nunit
    real(dp) :: A(3),abbr(9),abbr_t1,abbr_t2,alpha,Atube,B(3),C(3),coeffDiffSph(1:1000,2), &
         crossec(0:9),Dalv,deltaV(-1:9),DepFrac(4:6),h,j2,lduct,length,radius(0:9),Rin, &
         Vdep(0:6),Vduct(0:9),veloc(0:9),vfluid,volume(9),Vtot,Vtube,Z(3)
    real(dp) :: part_acinus_old(9)
    real(dp),allocatable :: part_concentration(:)
    logical :: flag
    character(len=60) :: sub_name

    sub_name = 'particle_deposition'
    call enter_exit(sub_name,1)

    allocate(part_concentration(num_nodes))
    
!!! Copy particle 'concentration' solution to temporary array
    part_concentration(:) = node_field(nj_conc1,:)
      
!!! TLC volume-dependent average diameter of an alveolus (Weibel, 1962) 
    Dalv = 1.54e-3_dp * current_volume**(1.0_dp/3.0_dp)

!!! inner radius for diffusion in alveoli related to inlet diameter (Hansen 1975)   
    Rin = Dalv*0.325_dp   
    abbr(7) = Rin/Dalv*2.0_dp
    abbr(8) = pi*abbr(7)

    ! abbreviation for acinar diffusive deposition
    abbr(2) = -4.0_dp*pi**2.0_dp*part_param%diffu/Dalv**2.0_dp 

    flag = .true.
    maxi = 1
    j = 0
    do while(flag)   ! find loop number for series in acinar diffusion
       j = j+1
       j2 = real(j)**2.0_dp
       coeffDiffSph(j,1) = (-1.0_dp)**(real(j)+1.0_dp)/j2* &
            (sin(real(j)*abbr(8))/real(j)/pi**2.0_dp- &
            cos(real(j)*abbr(8))*abbr(7)/pi) ! set up coefficients for sum to speed up
       coeffDiffSph(j,2) = j2*abbr(2)
       h = exp(coeffDiffSph(j,2)*dt)/j2                  
       if((h.lt.5.0e-7_dp).or.j.eq.500)then
          maxi = j+1         ! set loop number for series in acinar diffusion
          flag = .false.
       endif
    enddo ! while
    maxi = maxi-1 ! not here in the original

    if(inspiration)then
       npstart = 1           ! start from node 1 to calculate deposition -- loop in flow direction
       npend = num_nodes    ! REQUIRES correct order of nodes!!!!!!!!!!!!!!!!
       npstep = 1
    else
       npstart = num_nodes   ! start from terminal nodes to calculate deposition -- loop in flow direction
       npend = 1
       npstep = -1
    endif

!!! in below loops, assuming the length and radius changes have already been done
    do np = npstart,npend,npstep ! loop over all nodes in the model
       Vtot = 0.0_dp
       Vdep = 0.0_dp    ! deposition volume
       !!! loop over surrounding elements
       do noelem = 1,elems_at_node(np,0)
          ne = elems_at_node(np,noelem)
          if(np.eq.elem_nodes(1,ne))then
             np2 = elem_nodes(2,ne)
          else
             np2 = elem_nodes(1,ne)
          endif
          lduct = elem_field(ne_length,ne)
          !element (tube) length - since loop is over every node adjacent element, half length is used
          length = elem_field(ne_length,ne)/2.0_dp
          radius(0) = elem_field(ne_radius,ne) ! element (tube) radius
          Atube = pi*radius(0)**2.0_dp ! element (tube) cross section
          Vtube = Atube*length ! element (tube) volume
          Vtot = Vtot+Vtube  ! sum up total volume
          vfluid = elem_field(ne_flow,ne)/Atube  ! fluid velocity [mm/s]
       
!!! CALCULATE ANGLE BETWEEN TUBE AXIS AND GRAVITATION DIRECTION
          A(:) = node_xyz(:,np) !stores element coordinates in A,B
          B(:) = node_xyz(:,np2) !for angle to gravity calculation
          z = A
          Z(3) = A(3)-1.0_dp  ! unit vector in direction of gravity - requires respective cartesian coorsinates
          alpha = abs(angle_btwn_points(B,A,Z)) ! calculate angle between tube axis and gravity vector

!!! sedimentation
          Vdep(1) = Vdep(1) + abs(sin(alpha))*part_param%prho*&
                              part_param%gravityy*&
                              part_param%pdia**2.0_dp*&
               length*radius(0)*dt/9.0_dp/part_param%mu ! deposition volume due to sedimentation

!!! impaction
          if(inspiration)then ! impact only during inhalation (why?)
             ne0 = elem_cnct(-1,1,ne) ! parent element 
             if(ne0.gt.0)then ! parent element exists, means not inflow location
                np0 = elem_nodes(2,ne0) ! second node parent element (acc. to ordering second node in flow direction)
                ! mht - not sure about this. calculating as angle between current and parent           
                do i=1,2
                  do j=1,2
                     if(elem_nodes(i,ne).eq.elem_nodes(j,ne0))then ! find node that is part of parent and daughter 
                       np1 = elem_nodes(i,ne) 
                       np2 = elem_nodes(3-i,ne) 
                       np3 = elem_nodes(3-j,ne0)
                     endif
                  enddo   !j
                enddo      !i
!!! CALCULATE ANGLE TO DOWNSTREAM BRANCH AS WELL !WHAT ABOUT IMPACT DURING EXHALATION????
                A(:) = node_xyz(:,np1) ! A is the connection node of parent and daughter 
                B(:) = node_xyz(:,np2) 
                C(:) = node_xyz(:,np3)
                alpha = angle_btwn_points(B,A,C)-pi ! calculate angle between parent and daughter tube axis
             
                if(abs(sin(alpha)).gt.zero_tol)then
                   abbr(1) = 18.0_dp*part_param%mu/part_param%prho&
                                                  /part_param%pdia**2.0_dp&
                                                  /abs(sin(alpha))
                   h = 2.0_dp*abs(sin(alpha))/lduct*abs(vfluid)/ &
                        abbr(1)*abs(elem_field(ne_part_vel,ne))/ &   !!! warning: note that cmiss used nodal value for particle velocity
                        abs(elem_field(ne_radius,ne))**2.0_dp/pi  !!! but only one version for value (??)
                   if(abbr(1).lt.25.0_dp)then
                      h = h * (1.0_dp-exp(-abbr(1)*lduct/abs(vfluid)))
                   endif
                   Vdep(2) = Vdep(2)+2.0_dp*h*dt*length*2.0_dp*radius(0)
                endif
             endif         !ne0
          endif !inspiration

!!! brownian diffusion
          if(abs(vfluid).gt.zero_tol)then
             h = 2.0_dp/3.0_dp*(4.0_dp*part_param%diffu*lduct/abs(vfluid))**0.5_dp/pi ! average travelling distance within duct
          else
             h = 0.0_dp           ! equation not valid for keeping breath (??)
          endif
          Vdep(3) = Vdep(3)+pi*h*(2.0_dp*radius(0)-h)*length*dt/lduct*abs(vfluid) ! deposition with referencing to time step dt

!!! calculate deposition in acini
          if(num_units.gt.0.and.elem_cnct(1,0,ne).eq.0)then 
!!!         call acinar_deposition(ne,current_volume,dt)

            nunit = where_inlist(ne,units) ! get the unit number
            radius(0) = elem_field(ne_radius,ne)

            ! difference of FEM-mesh and measurements (HAEFELI-BLEUER,1988)
            abbr(3) = radius(0) - scale_radius_acinus(part_param%RacTLC(1),current_volume)
            abbr(4) = 0.0_dp  ! variable to store axial position for radius scaling

            veloc(0) = elem_field(ne_part_vel,ne)/pi/radius(0)**2.0_dp ! velocity terminal bronchi
            deltaV(-1) = dt*elem_field(ne_flow,ne) ! total volume change of acinus in one time step
            deltaV(0) = 0.0_dp ! assumption that volume change very small since no alveoli
            crossec(0) = radius(0)**2.0_dp*pi ! cross-sectional area of terminal bronchiole in FEM model  

            do gen = 1,9 ! loop over acinar generations - deposition efficiency the same in all gens
               abbr(4) = abbr(4)+part_param%LacTLC(gen+1) ! axial position of node

               ! scale duct radius wrt. current acinar volume and linear scaling to adapt HBW to FEM mesh
               radius(gen) = (part_param%totacinarLength-abbr(4)) &
                              /part_param%totacinarLength*abbr(3)+ &
                             scale_radius_acinus(part_param%RacTLC(gen+1),current_volume)
               crossec(gen) = pi*radius(gen)**2.0_dp*(2.0_dp**gen) ! accumulated duct cross-sectional area 
               Vduct(gen) = crossec(gen)*part_param%LacTLC(gen+1) ! duct volume in acinar region

               ! linear scaling of alveolar volume - can become negative!!
               volume(gen) = max(part_param%VacTLC(gen) &
                             /part_param%VtotTLC*unit_field(nu_vol,nunit)-Vduct(gen),0.0_dp)
               part_acinus_old(gen) = part_acinus_field(1+gen,nunit)
                
!!!.............diffusion in alveolar tissue
               abbr(5) = 0.0_dp
               abbr(6) = 0.0_dp
               flag = .true.
               j = 1
!              do j = 1,maxi
               do while(flag)
                   abbr_t1 = abbr(5)
                   abbr_t2 = abbr(6)
                   abbr(5) = abbr(5)+coeffDiffSph(j,1)*exp(coeffDiffSph(j,2)*part_acinus_field(10+gen,nunit))
                   abbr(6) = abbr(6)+coeffDiffSph(j,1)*exp(coeffDiffSph(j,2)*(part_acinus_field(10+gen,nunit)+dt))
                   if((abs(abbr_t1-abbr(5)).lt.1.0e-8_dp).and.(abs(abbr_t2-abbr(6)).lt.1.0e-8_dp)) flag = .false.
                   if(j.ge.maxi) flag = .false.
                   j = j + 1
                enddo !j

                if((abbr(6).gt.0.0_dp).and.(abbr(5).gt.abbr(6)))then ! can become zero for large T (accuracy of real*8)
                   ! diffusion fraction out of a sphere (Diffusion,Jost,1960) (0.853d0 is area correction ChoiKim2007)
                   DepFrac(4) = 1.0_dp-abbr(6)/abbr(5)*0.853_dp
                else
                   DepFrac(4) = 0.0_dp
                endif


!!!............ sedimentation in alveolar tissue
                ! deposition fraction due to sedimentation (0.853d0 is area correction ChoiKim2007) 
                DepFrac(5) = part_param%prho &
                             *part_param%gravityy &
                             *part_param%pdia**2.0_dp &
                             *dt/12.0_dp/part_param%mu/Dalv*0.853_dp   

                ! sum deposition fractions without mutually eliminating part
                DepFrac(6) = DepFrac(4)/2.0_dp+DepFrac(5)+DMAX1(DepFrac(4)/2.0_dp-DepFrac(5),0.0_dp)
                do j = 4,6
                  if(DepFrac(j).ge.0.2_dp)then
                    if(DepFrac(j).gt.1.0_dp) DepFrac(j) = 1.0_dp ! check if deposition volume extends volume in Acinus 
                   endif
                enddo !j
       
                if(gen.ge.1)THEN ! by increasing number, deposition in tubes can be taken into account
                  Vdep(4:6) = 0.0_dp
                  Vdep(5) = 0.0_dp
                else
                ! volume change each generation within on time step
                deltaV(gen) = part_param%VacTLC(gen)/part_param%VtotTLC*deltaV(-1) 

                ! velocities in acinar generations
                veloc(gen) = (veloc(gen-1)*crossec(gen-1)-deltaV(gen-1)/dt)/crossec(gen) 

                !! Sedimentation in ducts
                ! deposition volume due to sedimentation (statistic orientation to gravity)
                Vdep(5) = (2.0_dp**gen)*2.0_dp/pi &
                          *part_param%prho &
                          *part_param%gravityy &
                          *part_param%pdia**2.0_dp &
                          *part_param%LacTLC(gen+1)*radius(gen)*dt/9.0_dp/part_param%mu !*Ccun 

                !!Brownian diffusion (average travelling distance within duct)
                h = 2.0_dp/3.0_dp*(4.0_dp*part_param%diffu*part_param%LacTLC(gen+1)&
                                  /abs(veloc(gen)))**0.5_dp/pi 

                if(h.gt.radius(gen))then ! all particles are deposited
                   Vdep(4) = Vduct(gen)
                   Vdep(6) = Vduct(gen)
                else
                  ! deposition with referencing to time step DT
                  Vdep(4) = pi*h*(2.0_dp*radius(gen)-h)*dt*abs(veloc(gen))*2.0_dp**gen 
                    
                  !!Total alveolar deposition
                  ! sum deposition volumes without mutually eliminating volume
                  Vdep(6) = Vdep(4)/2.0_dp+Vdep(5)+DMAX1(Vdep(4)/2.0_dp-Vdep(5),0.0_dp)
                endif
                   
                 do j = 4,6
                    if(Vdep(j).gt.(0.2_dp*Vduct(gen)))then

                      ! check if deposition volume extends volume in Acinus
                      if(Vdep(j).gt.Vduct(gen)) Vdep(j)=Vduct(gen) 
                    endif
                    if(Vdep(j).ne.Vdep(j)) Vdep(j) = 0.0_dp ! function ISNAN does not work
                  enddo !j
                endif !gen.LT.9
       
                if(part_acinus_field(1+gen,nunit).ge.0.0_dp)then
                   ! NOTHING because thsi way NaN values are covered too
                else
                   part_acinus_field(1+gen,nunit) = 0.0_dp
                endif
                ! for deposition in the alveolar tissue the radial concentration profile is approximated by taking the next generation
                node_field(nj_loss,np) = node_field(nj_loss,np)+Vdep(6)*part_acinus_field(1+gen,nunit)+volume(gen)* &
                     DepFrac(6)*part_acinus_field(MIN(gen+2,10),nunit) ! store deposition quantity (mass in [g])
                node_field(nj_loss_dif,np) = node_field(nj_loss_dif,np)+Vdep(4)*part_acinus_field(1+gen,nunit)+volume(gen)* &
                     DepFrac(4)*part_acinus_field(MIN(gen+2,10),nunit) ! store diffusion quantity (mass in [g])
                node_field(nj_loss_sed,np) = node_field(nj_loss_sed,np)+Vdep(5)*part_acinus_field(1+gen,nunit)+volume(gen)* &
                     DepFrac(5)*part_acinus_field(MIN(gen+2,10),nunit) ! store sedimentaion quantity (mass in [g])
                part_acinus_field(1+gen,nunit) = part_acinus_field(1+gen,nunit)-volume(gen)*DepFrac(6)/(volume(gen)+Vduct(gen)) &
                     *part_acinus_field(MIN(gen+2,10),nunit)-Vdep(6)/(volume(gen)+Vduct(gen))*part_acinus_field(1+gen,nunit) ! assign new concentration acini
                if(part_acinus_field(1+gen,nunit).gt.0.0_dp)then
                   part_acinus_field(10+gen,nunit) = part_acinus_field(10+gen,nunit)+dt ! increase time particles spent in alveoli (in BBM(11..19,ne))
                else
                   part_acinus_field(10+gen,nunit) = 0 ! avoid division by zero
                endif   ! part_acinus_field
             enddo      ! for nine assumed acinar generations 
          endif         ! if a terminal element with acini attached
       enddo            ! np

       !!! summation of all deposition effects
       ! sum up diffusion and sedimentation without mutually eliminating volume
       Vdep(0) = Vdep(3)/2.0_dp + Vdep(1) + max(Vdep(3)/2.0_dp - Vdep(1),0.0_dp)
       ! add impaction with correction term
       Vdep(0) = Vdep(0) + Vdep(2) - Vdep(0) * Vdep(2)/Vtot        

       !!! store the deposition quantities of each type
       node_field(nj_loss,np) = node_field(nj_loss,np)+Vdep(0)*part_concentration(np) ! deposition quantity [g]
       node_field(nj_loss_dif,np) = node_field(nj_loss_dif,np)+Vdep(1)*part_concentration(np) ! deposition quantity [g]
       node_field(nj_loss_sed,np) = node_field(nj_loss_sed,np)+Vdep(2)*part_concentration(np) ! deposition quantity [g]
       node_field(nj_loss_imp,np) = node_field(nj_loss_imp,np)+Vdep(3)*part_concentration(np) ! deposition quantity [g]

!!! assign new particle concentration to all nodes
       part_concentration(np) = part_concentration(np)*(1.0_dp-Vdep(0)/Vtot)

    enddo !np
    
!!! copy nj_source (concentration - deposition) field to nj_conc1 field
    node_field(nj_conc1,:) = part_concentration(:)

    deallocate(part_concentration)

    call enter_exit(sub_name,2)

  end subroutine particle_deposition

  
!!! ############################################################3
!!! ############################################################3

  subroutine acinus_transport(nunit,current_volume,dt)
    
!!! Evaluates particle transport in acini governed by convection diffusion equation
!!! using FDM for 9 'virtual' nodes; using operator splitting (FDM for diffusion and
!!! dV/dt=const for each generation for convection);
!!! geometric data are based on HAEFELI-BLEUER and WEIBEL,1988
!!! Created September, 2011, Falko Schmidt      

!   use arrays, only: dp, particle_parameters, elem_field, node_field, &
!                     elem_nodes, part_acinus_field
    use arrays
    use indices
    use other_consts, only: pi
    use diagnostics, only: enter_exit

    implicit none
  
    type(particle_parameters) :: part_param

    integer :: nunit
    real(dp) :: current_volume,dt
!!!     Local variables
    integer :: ne,np,gen,genBEG,genEND,genSTEP
    real(dp) :: radius(-1:9),conc(-1:9), &
         volume(0:9),crossec(-1:9),veloc(0:9),abbr(6), &
         deltaV(-1:9),volflow,courant,effDiffu, &
         part_acinus_old(9)
    character(len=60) :: sub_name
    
    sub_name = 'acinus_transport'
    call enter_exit(sub_name,1)

    ne = units(nunit)
    np = elem_nodes(2,ne) ! the distal node; adjacent to the acinus
    
!!! initialize parameter - CAUTION! FEM model and data above may not match (radius adapted)
    radius(0) = elem_field(ne_radius,ne)
    crossec(0) = radius(0)**2.0_dp * pi                                  ! cross-sectional area terminal bronchiole FEM model  
    conc(0) = node_field(nj_conc1,np)                                    ! concentration in terminal bronchiole, defined that mass is concerved
    if(node_field(nj_conc1,np).lt.0.0_dp)then
       write(*,'('' Warning: particle concentration < 0 at unit'',i6)') nunit
       conc(0) = 0.0_dp
       node_field(nj_conc1,np) = 0.0_dp
    endif
    
    volflow = elem_field(ne_part_vel,ne)                                 ! volume flow particle in terminal bronchiole
    veloc(0) = volflow/crossec(0)

    ! total volume change in one time step
    deltaV(-1) = elem_field(ne_flow,ne)*dt

    ! assumptions volume change very small since no alveoli
    deltaV(0) = 0.0_dp 

    ! difference of FEM-mesh and measurements (HAEFELI-BLEUER,1988)
    abbr(1) = radius(0)-scale_radius_acinus(part_param%RacTLC(1),current_volume)

    ! variable to store axial position for radius scaling
    abbr(2) = 0.0_dp 

    ! ------------------------------------------------
    do gen = 1,9  ! loop over acinar generations
    ! ------------------------------------------------
       abbr(2) = abbr(2)+part_param%LacTLC(gen+1)                                     ! axial position of node
       radius(gen) = scale_radius_acinus(part_param%RacTLC(gen+1),current_volume)+ &     ! scale duct radius wrt. current acinar volume and 
            (part_param%totacinarLength-abbr(2))/part_param%totacinarLength*abbr(1)                        !          linear scaling to adapt HBW to FEM mesh
       crossec(gen) = pi*radius(gen)**2.0_dp*(2.0_dp**gen)               ! accumulated duct cross-sectional area    
       volume(gen) = DMAX1(part_param%VacTLC(gen)/part_param%VtotTLC*current_volume,0.0_dp)    ! linear scaling of volume
       deltaV(gen) = part_param%VacTLC(gen)/part_param%VtotTLC*deltaV(-1)                      ! volume change each generation within on time step
       conc(gen) = part_acinus_field(1+gen,nunit)                        ! concentration in acinar generations
       veloc(gen) = (veloc(gen-1)*crossec(gen-1)-deltaV(gen-1)/dt) &
                         /crossec(gen)                                   ! velocities corrected by volume change of previous generations
       courant = abs(veloc(gen))*dt/part_param%LacTLC(gen+1)                          ! courant number for stability check
       if(courant.ge.0.5_dp)then 
          write(*,'(''Transport in acinar region may be unstable, Cr='',d12.3, &
               '' decrease DT to at least '',d12.3,''s'')') &
               courant,0.5_dp*part_param%LacTLC(gen+1)/abs(veloc(gen))                ! warning if CFL-condition is not fulfilled
       endif
       part_acinus_old(gen) = part_acinus_field(1+gen,nunit)                      ! store old concentrations
    enddo !gen

!!! define direction-explicit transport
    if(volflow.ge.0.0_dp)THEN
       genBEG = 1                                                      ! start generation convective transport
       genEND = 9                                                      ! end generation convective transport
       genSTEP = 1                                                     ! direction convective transport
    else                                                               ! else even time step
       genBEG = 8                                                      ! start generation convective transport
       genEND = 0                                                      ! end generation convective transport
       genSTEP = -1                                                    ! direction convective transport
       conc(-1) = node_field(nj_conc1,elem_nodes(1,ne))                ! concentration of proximal node of terminal element
       part_param%LacTLC(1) = elem_field(ne_length,ne)                         ! length of terminal element of 3D tree
       crossec(-1) = crossec(0)                                       ! cross-sectional area of proximal node of terminal element
       volume(0) = crossec(0)*part_param%LacTLC(1)
    endif

!!! convective transport 
    do gen = genBEG,genEND,genSTEP 
       if(gen.le.8)then                                                ! inside domain use central scheme in space
          conc(gen) = (conc(gen)*(volume(gen)-deltaV(gen))+dt* &
               (conc(MIN0(gen-genSTEP,gen))*veloc(gen)*crossec(gen)- &
               (conc(MAX0(gen-genSTEP,gen))*veloc(gen+1)* &
               crossec(gen+1))))/volume(gen)                           ! new concentration by mass/volume balance
       else
          conc(gen) = (conc(gen)*(volume(gen)-deltaV(gen))+ &
               conc(MIN0(gen-genSTEP,gen))*veloc(gen)* &
               crossec(gen)*dt)/(volume(gen))                          ! new concentration in last generation by mass/volume balance
       endif
    enddo

!!! explicit FDM scheme - for diffusion
    do gen = 1,9                                                         ! loop over acinar generations in alternating directions
       effDiffu = part_param%diffu+0.0867_dp*abs(veloc(gen))*part_param%LacTLC(gen+1)            ! effective diffusivity acc. to. Lee(2001)
       if(gen.le.8)then                                                  ! inside domain use central scheme in space
          abbr(1) = part_param%LacTLC(gen+1)*part_param%LacTLC(gen+2)*(part_param%LacTLC(gen+1)+part_param%LacTLC(gen+2))! abbreviation to reduce computational effort
          abbr(2) = part_param%LacTLC(gen+1)*dt/volume(gen)/abbr(1)     
          part_acinus_field(1+gen,nunit) = conc(gen)+abbr(2)*effDiffu* &
               ((part_param%LacTLC(gen+1)**2.0_dp*(crossec(gen+1)-crossec(gen)) &
               -part_param%LacTLC(gen+2)**2.0_dp*(crossec(gen)-crossec(gen-1))) &
               *(part_param%LacTLC(gen+1)**2.0_dp*(conc(gen+1)-conc(gen)) &
               -part_param%LacTLC(gen+2)**2.0_dp*(conc(gen)-conc(gen-1)))/abbr(1) &
               +2.0_dp*crossec(gen)* &
               (part_param%LacTLC(gen+1)*(conc(gen+1)-conc(gen)) &
               -part_param%LacTLC(gen+2)*(conc(gen)-conc(gen-1))))                  ! FDM with central differences in space  
       else                                                              ! domain boundary use different scheme
          abbr(1) = part_param%LacTLC(gen+1)**2.0_dp                                  ! abbreviation to reduce computational effort
          abbr(2) = part_param%LacTLC(gen+1)*dt/volume(gen)/abbr(1)     
          part_acinus_field(1+gen,nunit) = conc(gen)+abbr(2)*effDiffu* &
               (2.0_dp*crossec(gen)*(conc(gen-1)-conc(gen)))             ! FDM with entral differences in space at boundary (using condition dC/dx=0)
       endif
    enddo !gen

!!! assign new concentrations to 'part_acinus_field' and calculate average time particles spend in acinus
    do gen = 1,9                                        
       if(part_acinus_field(1+gen,nunit).lt.0.0_dp) part_acinus_field(1+gen,nunit) = 0.0_dp
       if((volflow.gt.0.0_dp).and.(part_acinus_field(1+gen,nunit).gt.0.0_dp))then  ! average time does not change when particles leave alveolus
          part_acinus_field(10+gen,nunit) = part_acinus_old(gen)*part_acinus_field(10+gen,nunit)* &
               (volume(gen)-deltaV(gen))/ &
               part_acinus_field(1+gen,nunit)/volume(gen)                   ! calculate new average time partices spent in alveoli (in part_acinus_field(11..19,ne))
       elseif(part_acinus_field(1+gen,nunit).le.0.0_dp)then                 ! concentration is zero
          part_acinus_field(10+gen,nunit)=0.0_dp                            ! reset time
       endif ! volflow
    enddo

!!! for exhalation calculate concentration of terminal node of branching tree
    if(volflow.lt.0)then                                                 ! exhalation
       gen = 0
       effDiffu = part_param%diffu+0.0867_dp*abs(veloc(gen))*part_param%LacTLC(gen+1)            ! effective diffusivity acc. to. Lee(2001)

       abbr(1) = part_param%LacTLC(gen+1)*part_param%LacTLC(gen+2)*&
                 (part_param%LacTLC(gen+1)+part_param%LacTLC(gen+2))       ! abbreviation 
       abbr(2) = part_param%LacTLC(gen+1)*dt/part_param%LacTLC(gen+1)/crossec(gen)/abbr(1)                 

       node_field(nj_conc1,np) = conc(gen)+abbr(2)*effDiffu* &
            ((part_param%LacTLC(gen+1)**2.0_dp*(crossec(gen+1)-crossec(gen)) &
            -part_param%LacTLC(gen+2)**2.0_dp*(crossec(gen)-crossec(gen-1))) &
            *(part_param%LacTLC(gen+1)**2.0_dp*(conc(gen+1)-conc(gen)) &
            -part_param%LacTLC(gen+2)**2.0_dp*(conc(gen)-conc(gen-1)))/abbr(1) &
            +2.0_dp*crossec(gen)* &
            (part_param%LacTLC(gen+1)*(conc(gen+1)-conc(gen)) &
            -part_param%LacTLC(gen+2)*(conc(gen)-conc(gen-1))))                     ! FDM with central differences in space  

        if(node_field(nj_conc1,np).lt.0.0_dp)then
           write(*,'(''ACIN_TRANSP exhaling concentration < 0 ='',d12.3)') node_field(nj_conc1,np)
          node_field(nj_conc1,np) = 0.0_dp
       endif
    endif

    call enter_exit(sub_name,2)

  end subroutine acinus_transport

!#######################################################################
  
  function scale_radius_acinus(initial_radius,volume_lung)

!!!   calculates new airway diameter from 
!!!   lung volume, equation is linear interploation of HABIB1994
!!!    REQUIRES reference data stored are from TLC
    use arrays, only: dp
    
    implicit none
  
!!!     Parameter List
    real(dp),intent(in) :: initial_radius,volume_lung
!!!     Local variables
    real(dp) :: alpha,beta,dref,fac
    real(dp),parameter :: TLC = 6.0e+3_dp, FRC = 3.6e+3_dp
    real(dp) :: scale_radius_acinus
      
    alpha = 1.51_dp ! parameter 
    beta = 0.442_dp ! parameter
    dref = 2.0_dp     ! reference diameter [mm]

    fac = alpha/(1.0_dp+(2.0_dp*initial_radius/dref)**(-beta))          ! scaling acc.to HABIB1994
    fac = DMIN1(fac,1.0_dp)                                             ! trachea otherwise would increase
    scale_radius_acinus = initial_radius*(volume_lung-FRC+fac*(TLC-volume_lung))/(TLC-FRC)  ! linear interploation for all lung volumes      

  end function scale_radius_acinus

!#########################################################################

  subroutine calc_mass_particles(nj,nu_field,gas_mass,deposit_mass)
    
!   use arrays, only: dp,num_nodes,elem_nodes,node_field,elem_ordrs,&
!                     elem_field,elem_symmetry,num_elems
    use arrays
    use indices
    use diagnostics, only: enter_exit

    implicit none
  
    integer,intent(in) :: nj,nu_field
    real(dp) :: gas_mass,deposit_mass
    !     Local Variables
    integer :: i,ne,ne0,ngen,nmax_gen,np,np1,np2,nunit
    real(dp) :: average_conc,mass_by_gen(50),sum_mass,vol_by_gen(50)
    real(dp),allocatable :: tree_mass(:),unit_tree_mass(:),wall_mass(:),unit_wall_mass(:)
    character(len=60):: sub_name

    !...........................................................................
    
    sub_name = 'calc_mass_particles'
    call enter_exit(sub_name,1)
    
    if(.not.allocated(tree_mass)) allocate(tree_mass(num_nodes))
    if(.not.allocated(unit_tree_mass)) allocate(unit_tree_mass(num_nodes))
    if(.not.allocated(wall_mass)) allocate(wall_mass(num_nodes))
    if(.not.allocated(unit_wall_mass)) allocate(unit_wall_mass(num_nodes))

    tree_mass = 0.0_dp
    unit_tree_mass = 0.0_dp
    wall_mass = 0.0_dp
    unit_wall_mass = 0.0_dp
    mass_by_gen = 0.0_dp
    vol_by_gen = 0.0_dp

    ! initialise to the mass in each element
    do ne = 1,num_elems
       np1 = elem_nodes(1,ne)
       np2 = elem_nodes(2,ne)
       average_conc = (node_field(nj,np1)+node_field(nj,np2))/2.0_dp
       tree_mass(ne) = average_conc*elem_field(ne_vol,ne) ! mmol/mm^3 * mm^3
       wall_mass(ne) = (node_field(nj_loss,np1)+node_field(nj_loss,np2))/2.0_dp
       if(elem_ordrs(no_type,ne).eq.1)then ! conducting airways
          elem_field(ne_mass,ne) = tree_mass(ne)
       else
          elem_field(ne_mass,ne) = 0.0_dp
       endif
       ngen = elem_ordrs(1,ne)
       mass_by_gen(ngen) = mass_by_gen(ngen) + tree_mass(ne)
       vol_by_gen(ngen) = vol_by_gen(ngen) + elem_field(ne_vol,ne)
    enddo
    nmax_gen = maxval(elem_ordrs)
    
!    elem_field(ne_mass,1:num_elems) = tree_mass(1:num_elems)

    ! add the mass in each elastic unit to terminal elements
    do nunit = 1,num_units
       ne = units(nunit)
       if(ne.ne.0.and.elem_cnct(1,0,ne).eq.0)then
          tree_mass(ne) = tree_mass(ne) + &
            unit_field(nu_vol,nunit)*unit_field(nu_field,nunit)
          unit_tree_mass(ne) = unit_tree_mass(ne) + &
               unit_field(nu_vol,nunit)*unit_field(nu_field,nunit)
       endif
    enddo
    
    ! sum mass recursively up the tree
    do ne = num_elems,2,-1 ! not for the stem branch; parent = 0
       ne0 = elem_cnct(-1,1,ne)
       tree_mass(ne0) = tree_mass(ne0) + dble(elem_symmetry(ne))*tree_mass(ne)
       wall_mass(ne0) = wall_mass(ne0) + dble(elem_symmetry(ne))*wall_mass(ne)
       unit_tree_mass(ne0) = unit_tree_mass(ne0) + dble(elem_symmetry(ne))*unit_tree_mass(ne)
       elem_field(ne_mass,ne0) = elem_field(ne_mass,ne0) + dble(elem_symmetry(ne))*elem_field(ne_mass,ne)
    enddo !noelem

    gas_mass = tree_mass(1)
    deposit_mass = wall_mass(1)
    
!    write(*,'('' Mass in generations: '',10(f12.4),'' |'',f12.4)') (mass_by_gen(i),i=1,10),unit_tree_mass(1)
!    write(*,'(''  Vol in generations: '',10(f12.4))') (vol_by_gen(i),i=1,10)
    
    deallocate(tree_mass)
    deallocate(unit_tree_mass)
    deallocate(wall_mass)
    deallocate(unit_wall_mass)

    call enter_exit(sub_name,2)
    
  end subroutine calc_mass_particles
  

!#########################################################################

  function inlist(item,ilist)
!!! dummy arguments
    integer :: item,ilist(:)
! local variables
    integer :: n
    logical :: inlist

    inlist = .false.
    do n=1,size(ilist)
       if(item == ilist(n)) inlist = .true.
    enddo

  end function inlist

!!!###############################################################

  function where_inlist(item,ilist)
!!! dummy arguments
    integer :: item,ilist(:)
! local variables
    integer :: n
    integer :: where_inlist

    do n=1,size(ilist)
       if(item == ilist(n)) where_inlist = n
    enddo

  end function where_inlist

!#########################################################################
!!!####################################################################

  subroutine airway_mesh_deform(dt,initial_volume,coupled,problem_type)

!!! changes the size of elastic units and airways, based on pre-computed
!!! fields for flow into the units (ne_flow) and element volume change (ne_dvdt).
!!! The concentration in the units is adjusted to make sure that mass is conserved
!!! when the unit changes volume. If a unit changes in volume, then both the radius
!!! and length scale as the cube root of volume change. For alveolated airways (where
!!! a/A is < 1) the outer radius (indexed as nj_radius) is scaled and a/A unchanged.

!!! note that scaling the length or scaling length and radius gives no difference in mass error

    use geometry,only: volume_of_mesh
    use arrays
    use diagnostics, only: enter_exit

    implicit none

    real(dp),intent(in) :: dt,initial_volume
    logical,intent(in) :: coupled
    character(len=9) :: problem_type

    ! Local parameters
    integer :: i,ne,ne0,np,np1,np2,nunit
    real(dp) :: current_volume,new_volume,previous_volume,ratio, &
         sum_old_vols,tree_volume,unit_mass(2),volume_change
    character(len=60) :: sub_name


    sub_name = 'airway_mesh_deform'
    call enter_exit(sub_name,1)

    call volume_of_mesh(current_volume,tree_volume)

!!! calculate the change in volume of elastic units:
    ! mass(unit) = volume(unit) * concentration(unit)
    ! dv(unit) = flow(unit) * dt
    ! volume(new) = volume(old) + dv
    ! concentration(new) = mass/volume(new)
    do nunit = 1,num_units !for each of the elastic units
       ne = units(nunit) ! local element number
       if(ne.ne.0)then
          np = elem_nodes(2,ne) ! end node, attaches to unit
          volume_change = elem_field(ne_flow,ne)*dt ! dv = q(unit)*dt
          if(.not.coupled)then ! this is where the unit volume is updated
             previous_volume = unit_field(nu_vol,nunit)
             new_volume = unit_field(nu_vol,nunit) + volume_change
             unit_field(nu_vol,nunit) = new_volume
          else ! for coupled model, volume has already been updated
             previous_volume = unit_field(nu_vol,nunit) - volume_change
             new_volume = unit_field(nu_vol,nunit)
          endif

          if(elem_cnct(1,0,ne).eq.0)then ! no child branches
             unit_mass(1) = previous_volume*unit_field(nu_conc1,nunit) ! m = v(unit)*c
             unit_mass(2) = previous_volume*unit_field(nu_conc2,nunit) ! m = v(unit)*c
             unit_field(nu_vol,nunit) = new_volume ! i.e. unchanged if coupled
             ! important: adjust the concentration to maintain mass conservation
             if(.not.coupled)then
                unit_field(nu_conc1,nunit) = unit_mass(1)/unit_field(nu_vol,nunit)
                unit_field(nu_conc2,nunit) = unit_mass(2)/unit_field(nu_vol,nunit)
             endif

             ! for problem type particles
             call acinus_transport(nunit,current_volume,dt)
          endif
       endif
    enddo

!!! calculate the inflation/deflation of multi-branching acini (or other branches)
    ! dv(element) = dvdt(element) * total_volume_change
    ! ratio_of_volumes = (volume(element)+dv(element))/volume(element)
    ! volume(new) = volume(old) * ratio_of_volumes
    ! length(new) = length(old) * cuberoot(ratio_of_volumes)
    ! radius(new) = radius(old) * cuberoot(ratio_of_volumes)
    if(.not.coupled)then
       do ne = 1,num_elems
          !mht!       volume_change = elem_field(ne_dvdt,ne)*inlet_flow*dt
          np1 = elem_nodes(1,ne)
          np2 = elem_nodes(2,ne)
          volume_change = elem_field(ne_dvdt,ne)*dt

          ! adjust the concentrations at both nodes to maintain the mass
          if(abs(volume_change).gt.zero_tol)then
             ne0 = elem_cnct(-1,1,ne) ! parent element
             sum_old_vols = elem_field(ne_vol,ne0)
             do i = 1,elem_cnct(1,0,ne0)
                sum_old_vols = sum_old_vols + elem_field(ne_vol,elem_cnct(1,i,ne0))
             enddo
             node_field(nj_conc1,np1) = node_field(nj_conc1,np1) * sum_old_vols/&
                  (sum_old_vols+volume_change)

             ne0 = ne
             sum_old_vols = elem_field(ne_vol,ne0)
             do i = 1,elem_cnct(1,0,ne0)
                sum_old_vols = sum_old_vols + elem_field(ne_vol,elem_cnct(1,i,ne0))
             enddo
             node_field(nj_conc1,np2) = node_field(nj_conc1,np2) * sum_old_vols/&
                  (sum_old_vols+volume_change)
          endif

          ! adjust the element size
          ratio = (volume_change+elem_field(ne_vol,ne))/elem_field(ne_vol,ne)
          elem_field(ne_vol,ne) = elem_field(ne_vol,ne)*ratio
          elem_field(ne_length,ne) = elem_field(ne_length,ne)*(ratio)**(1.0_dp/3.0_dp)
          elem_field(ne_radius,ne) = elem_field(ne_radius,ne)*(ratio)**(1.0_dp/3.0_dp)

       enddo
    endif

!!! calculate the total model volume
    call volume_of_mesh(current_volume,tree_volume)

    total_volume_change = current_volume-initial_volume

    call enter_exit(sub_name,2)

  end subroutine airway_mesh_deform


! #################################################################
end module particle

