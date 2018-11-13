module species_transport_c
  implicit none
  private

contains

!!!##########################################################################

  subroutine initialise_transport_c(MODEL, model_len) bind(C, name="initialise_transport_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN
    use species_transport, only: initialise_transport
    implicit none

    integer,intent(in) ::  model_len
    type(c_ptr), value, intent(in) :: MODEL
    character(len=MAX_FILENAME_LEN) :: model_f

    call strncpy(model_f, MODEL, model_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_initialise_transport(model_f)
#else
    call initialise_transport(model_f)
#endif

  end subroutine initialise_transport_c


!!!##########################################################################

  subroutine solve_transport_c(MODEL, model_len) bind(C, name="solve_transport_c")

    use iso_c_binding, only: c_ptr
    use utils_c, only: strncpy
    use other_consts, only: MAX_FILENAME_LEN
    use species_transport, only: solve_transport
    implicit none

    integer,intent(in) ::  model_len
    type(c_ptr), value, intent(in) :: MODEL
    character(len=MAX_FILENAME_LEN) :: model_f

    call strncpy(model_f, MODEL, model_len)

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_solve_transport(model_f)
#else
    call solve_transport(model_f)
#endif

  end subroutine solve_transport_c


!!!##########################################################################

end module species_transport_c
