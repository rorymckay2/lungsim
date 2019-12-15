module particle_c
  implicit none
  private

contains 

!!!###################################################################################

  subroutine initialise_particles_c(fileid) bind(C, name="initialise_particles_c")

    use particle, only: initialise_particles
    implicit none

    integer,intent(in) :: fileid
    

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_initialise_particles(fileid)
#else
    call initialise_particles(fileid)
#endif

  end subroutine initialise_particles_c


!###################################################################################

  subroutine controller_particlesolve_c(fileid) bind(C, name="controller_particlesolve_c")

    use particle, only: controller_particlesolve
    implicit none

    integer,intent(in) :: fileid
    

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_controller_particlesolve(fileid)
#else
    call controller_particlesolve(fileid)
#endif

  end subroutine controller_particlesolve_c

!###################################################################################
end module particle_c
