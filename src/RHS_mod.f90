module RHS_mod
  use types_mod, only: DP

  implicit none

  private

  public :: func 

contains

  function func() result (d)
    implicit none

    ! integer :: j
    ! real (kind=DP), dimension(:) :: x
    real (kind=DP) :: d

    d = 0.0e+00_DP
  end function

end module RHS_mod
