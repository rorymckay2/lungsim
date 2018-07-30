module arrays_c
!*Brief Description:* This module wraps part of the arrays module that require a c interface
!
!*LICENSE:*
!
!
!*Contributor(s):* Merryn Tawhai, Alys Clark
!
!*Full Description:*
!
!This module wraps part of the arrays module that require a c interface
  use arrays,only: dp

  implicit none
  public set_node_field_value_c

contains
  subroutine set_node_field_value_c(row, col, value) bind(C, name="set_node_field_value_c")
    use arrays, only: set_node_field_value
    implicit none

    integer, intent(in) :: row, col
    real(dp), intent(in) :: value

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_set_node_field_value(row, col, value)
#else
    call set_node_field_value(row, col, value)
#endif

  end subroutine set_node_field_value_c


  subroutine check_node_xyz_2d_c(row, col, value) bind(C, name="check_node_xyz_2d_c")
    use arrays, only: check_node_xyz_2d
    implicit none

    integer, intent(in) :: row, col
    real(dp), intent(out) :: value

#if defined _WIN32 && defined __INTEL_COMPILER
    call so_check_node_xyz_2d(row, col, value)
#else
    call check_node_xyz_2d(row, col, value)
#endif

  end subroutine check_node_xyz_2d_c


end module arrays_c
