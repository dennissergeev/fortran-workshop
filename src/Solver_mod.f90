module Solver_mod
  use types_mod, only: DP

  implicit none

  private

  public :: fd1d_heat_explicit 

contains

  subroutine fd1d_heat_explicit(x, dt, cfl, h, h_new)
    implicit none

    real (kind=DP), intent (in) :: cfl
    real (kind=DP), intent (in) :: dt
    ! real (kind=DP), intent (in) :: t
    real (kind=DP), dimension(:), intent (in)  :: h
    real (kind=DP), dimension(:), intent (in)  :: x
    real (kind=DP), dimension(:), intent (out) :: h_new

    integer :: j
    real (kind=DP) :: f(size(x, 1))
    integer :: xn

    xn = size(x, 1)

    do j = 1, xn
      f(j) = func()
    end do

    h_new(1) = 0.0e+00_DP

    do j = 2, xn - 1
      h_new(j) = h(j) + dt*f(j) + cfl*(h(j-1)-2.0e+00_DP*h(j)+h(j+1))
    end do

! set the boundary conditions again
    h_new(1) = 90.0e+00_DP
    h_new(X_NUM) = 70.0e+00_DP
  end subroutine

end module Solver_mod
