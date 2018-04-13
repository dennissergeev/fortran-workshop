!> Module to calculate the CFL number
module CFL_mod
  use types_mod, only: DP, SI

  implicit none

  private

  public :: fd1d_heat_explicit_cfl 

contains
  !> calculate the CFL number
  !> \( \text{CFL} = \kappa\frac{\Delta t}{\Delta x^2} \)
  subroutine fd1d_heat_explicit_cfl(k, T_NUM, t_min, t_max, X_NUM, x_min, &
    x_max, cfl)

    implicit none

    !> number of intervals in t-axis
    integer(kind=SI), intent (in) :: T_NUM
    !> number of intervals in x-axis
    integer(kind=SI), intent (in) :: X_NUM
    !> the heat constant \( \kappa \)
    real (kind=DP), intent (in) :: k
    !> upper bound of t-axis
    real (kind=DP), intent (in) :: t_max
    !> lower bound of t-axis
    real (kind=DP), intent (in) :: t_min
    !> upper bound of x-axis
    real (kind=DP), intent (in) :: x_max
    !> lower bound of x-axis
    real (kind=DP), intent (in) :: x_min
    !> the CFL number
    real (kind=DP), intent (out) :: cfl
    real (kind=DP) :: dx
    real (kind=DP) :: dt

    dx = (x_max-x_min)/real(X_NUM-1, kind=DP)
    dt = (t_max-t_min)/real(T_NUM-1, kind=DP)

    cfl = k*dt/dx/dx

    write (*, '(a)') ' '
    write (*, '(a,g14.6)') '  CFL stability criterion value = ', cfl

  end subroutine

end module CFL_mod
