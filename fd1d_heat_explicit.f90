program fd1d_heat_explicit_prb
  use :: types_mod, only: DP

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
  k = 0.002e+00_dp

! the x-range values
  x_min = 0.0e+00_dp
  x_max = 1.0e+00_dp
! X_NUM is the number of intervals in the x-direction
  call r8vec_linspace(x_min, x_max, x)

! the t-range values. integrate from t_min to t_max
  t_min = 0.0e+00_dp
  t_max = 80.0e+00_dp

! T_NUM is the number of intervals in the t-direction
  dt = (t_max-t_min)/real(T_NUM-1, kind=DP)
  call r8vec_linspace(t_min, t_max, t)

! get the CFL coefficient
  call fd1d_heat_explicit_cfl(k, T_NUM, t_min, t_max, X_NUM, x_min, x_max, &
    cfl)

  if (0.5e+00_dp<=cfl) then
    write (*, '(a)') ' '
    write (*, '(a)') 'FD1D_HEAT_EXPLICIT_CFL - Fatal error!'
    write (*, '(a)') '  CFL condition failed.'
    write (*, '(a)') '  0.5 <= K * dT / dX / dX = CFL.'
    stop
  end if

! set the initial condition
  do j = 1, X_NUM
    h(j) = 50.0e+00_dp
  end do

! set the bounday condition
  h(1) = 90.0e+00_dp
  h(X_NUM) = 70.0e+00_dp

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

contains

  function func() result (d)
    implicit none

    ! integer :: j
    ! real (kind=DP), dimension(:) :: x
    real (kind=DP) :: d

    d = 0.0e+00_dp
  end function

  subroutine fd1d_heat_explicit(x, dt, cfl, h, h_new)
    implicit none

    real (kind=DP), intent (in) :: cfl
    real (kind=DP), intent (in) :: dt
    ! real (kind=DP), intent (in) :: t
    real (kind=DP), dimension(:), intent (in)  :: h
    real (kind=DP), dimension(:), intent (in)  :: x
    real (kind=DP), dimension(:), intent (out) :: h_new

    integer :: j
    real (kind=DP) :: f(size(x))
    integer :: xn

    xn = size(x, 1)

    do j = 1, xn
      f(j) = func()
    end do

    h_new(1) = 0.0e+00_dp

    do j = 2, xn - 1
      h_new(j) = h(j) + dt*f(j) + cfl*(h(j-1)-2.0e+00_dp*h(j)+h(j+1))
    end do

! set the boundary conditions again
    h_new(1) = 90.0e+00_dp
    h_new(X_NUM) = 70.0e+00_dp
  end subroutine

  subroutine fd1d_heat_explicit_cfl(k, T_NUM, t_min, t_max, X_NUM, x_min, &
    x_max, cfl)

    implicit none

    integer :: T_NUM
    integer :: X_NUM
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

  subroutine r8mat_write(output_filename, table)
    implicit none

    character (len=*), intent(in) :: output_filename
    real (kind=DP)   , dimension(:, :), intent(in) :: table
    integer :: m
    integer :: n
    integer :: j
    integer :: output_unit_id
    character (len=30) :: string

    m = size(table, 1)
    n = size(table, 2)

    output_unit_id = 10
    open (unit=output_unit_id, file=output_filename, status='replace')

    write (string, '(a1,i8,a1,i8,a1,i8,a1)') '(', m, 'g', 24, '.', 16, ')'

    do j = 1, n
      write (output_unit_id, string) table(1:m, j)
    end do

    close (unit=output_unit_id)
  end subroutine

  subroutine r8vec_linspace(a_first, a_last, a)

    implicit none

    real (kind=DP), intent(in)  :: a_first
    real (kind=DP), intent(in)  :: a_last
    real (kind=DP), dimension(:), intent(out) :: a
    integer :: i
    integer :: n

    n = size(a, 1)

    do i = 1, n
      a(i) = (real(n-i,kind=DP)*a_first+real(i-1,kind=DP)*a_last)/ &
        real(n-1, kind=DP)
    end do

  end subroutine

  subroutine r8vec_write(output_filename, x)

    implicit none

    character (len=*), intent(in) :: output_filename
    real (kind=DP)   , dimension(:), intent(in) :: x
    integer :: output_unit_id
    integer :: j

    output_unit_id = 11
    open (unit=output_unit_id, file=output_filename, status='replace')

    do j = 1, size(x, 1)
      write (output_unit_id, '(2x,g24.16)') x(j)
    end do

    close (unit=output_unit_id)
  end subroutine

end program
