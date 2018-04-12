program fd1d_heat_explicit_prb
  use types_mod, only: DP
  use rhs_mod, only: func
  use cfl_mod, only: fd1d_heat_explicit_cfl
  use io_mod, only: r8vec_linspace, r8mat_write, r8vec_write
  use solver_mod, only, fd1d_heat_explicit

  implicit none

  integer, parameter :: T_NUM = 201
  integer, parameter :: X_NUM = 21

  real (kind=DP)                  :: cfl
  real (kind=DP)                  :: dt
  real (kind=DP), allocatable, dimension(:)    :: h
  real (kind=DP), allocatable, dimension(:)    :: h_new
! the "matrix" stores all x-values for all t-values
! remember Fortran is column major, meaning that rows are contiguous
  real (kind=DP), allocatable, dimension(:, :) :: hmat
  integer                         :: i
  integer                         :: j
  real (kind=DP)                  :: k

  real (kind=DP), allocatable, dimension(:)    :: t
  real (kind=DP)                  :: t_max
  real (kind=DP)                  :: t_min
  real (kind=DP), allocatable, dimension(:)    :: x
  real (kind=DP)                  :: x_max
  real (kind=DP)                  :: x_min
  integer                         :: ierr 


  allocate(h    (1:X_NUM         ), stat=ierr)
  allocate(h_new(1:X_NUM         ), stat=ierr)
  allocate(hmat (1:X_NUM, 1:T_NUM), stat=ierr)
  allocate(t    (         1:T_NUM), stat=ierr)
  allocate(x    (1:X_NUM         ), stat=ierr)
  if ( ierr /= 0 ) then
    write(*, '(a)') 'Failed to allocate arrays'
  end if

  write (*, '(a)') ' '
  write (*, '(a)') 'FD1D_HEAT_EXPLICIT_PRB:'
  write (*, '(a)') '  FORTRAN77 version.'
  write (*, '(a)') '  Test the FD1D_HEAT_EXPLICIT library.'

  write (*, '(a)') ' '
  write (*, '(a)') 'FD1D_HEAT_EXPLICIT_PRB:'
  write (*, '(a)') '  Normal end of execution.'
  write (*, '(a)') ' '

  write (*, '(a)') ' '
  write (*, '(a)') 'FD1D_HEAT_EXPLICIT_TEST01:'
  write (*, '(a)') '  Compute an approximate solution to the time-dependent'
  write (*, '(a)') '  one dimensional heat equation:'
  write (*, '(a)') ' '
  write (*, '(a)') '    dH/dt - K * d2H/dx2 = f(x,t)'
  write (*, '(a)') ' '
  write (*, '(a)') '  Run a simple test case.'

! heat coefficient
  k = 0.002e+00_DP

! the x-range values
  x_min = 0.0e+00_DP
  x_max = 1.0e+00_DP
! X_NUM is the number of intervals in the x-direction
  call r8vec_linspace(x_min, x_max, x)

! the t-range values. integrate from t_min to t_max
  t_min = 0.0e+00_DP
  t_max = 80.0e+00_DP

! T_NUM is the number of intervals in the t-direction
  dt = (t_max-t_min)/real(T_NUM-1, kind=DP)
  call r8vec_linspace(t_min, t_max, t)

! get the CFL coefficient
  call fd1d_heat_explicit_cfl(k, T_NUM, t_min, t_max, X_NUM, x_min, x_max, &
    cfl)

  if (0.5e+00_DP<=cfl) then
    write (*, '(a)') ' '
    write (*, '(a)') 'FD1D_HEAT_EXPLICIT_CFL - Fatal error!'
    write (*, '(a)') '  CFL condition failed.'
    write (*, '(a)') '  0.5 <= K * dT / dX / dX = CFL.'
    stop
  end if

! set the initial condition
  do j = 1, X_NUM
    h(j) = 50.0e+00_DP
  end do

! set the bounday condition
  h(1) = 90.0e+00_DP
  h(X_NUM) = 70.0e+00_DP

! initialise the matrix to the initial condition
  do i = 1, X_NUM
    hmat(i, 1) = h(i)
  end do

! the main time integration loop 
  do j = 2, T_NUM
    call fd1d_heat_explicit(x, dt, cfl, h, h_new)
    ! call fd1d_heat_explicit(x, t(j-1), dt, cfl, h, h_new)

    do i = 1, X_NUM
      hmat(i, j) = h_new(i)
      h(i) = h_new(i)
    end do
  end do

! write data to files
  call r8mat_write('h_test01.txt', hmat)
  call r8vec_write('t_test01.txt', t)
  call r8vec_write('x_test01.txt', x)

  deallocate(h    )
  deallocate(h_new)
  deallocate(hmat )
  deallocate(t    )
  deallocate(x    )

end program fd1d_heat_explicit_prb
