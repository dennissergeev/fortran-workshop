module CFL_mod
  use types_mod, only: DP

  implicit none

  private

  public :: fd1d_heat_explicit_cfl 

contains

  subroutine fd1d_heat_explicit_cfl(k, T_NUM, t_min, t_max, X_NUM, x_min, &
    x_max, cfl)

    implicit none

    integer       , intent (in) :: T_NUM
    integer       , intent (in) :: X_NUM
    real (kind=DP), intent (in) :: k
    real (kind=DP), intent (in) :: t_max
    real (kind=DP), intent (in) :: t_min
    real (kind=DP), intent (in) :: x_max
    real (kind=DP), intent (in) :: x_min
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
