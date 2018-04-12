module IO_mod
  use types_mod, only: DP

  implicit none

  private

  public :: r8mat_write
  public :: r8vec_linspace
  public :: r8vec_write 

contains

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

end module IO_mod
