module IO_mod
  use netcdf
  use types_mod, only: DP, SI

  implicit none

  private

  public :: r8mat_write
  public :: r8vec_linspace
  public :: r8vec_write 

contains

  subroutine r8mat_write(output_filename, t, x, table)
    implicit none

    character (len=*), intent(in) :: output_filename
    real(kind=DP), dimension(:), intent(in) :: x
    real(kind=DP), dimension(:), intent(in) :: t
    real (kind=DP)   , dimension(:, :), intent(in) :: table
    integer(kind=SI) :: nx
    integer(kind=SI) :: nt
    integer(kind=SI) :: ierr
    integer(kind=SI) :: ncid
    integer(kind=SI) :: t_dimid
    integer(kind=SI) :: x_dimid
    integer(kind=SI) :: x_id, t_id, varid

    nx = size(table, 1)
    nt = size(table, 2)

    ierr = NF90_CREATE( trim(output_filename), NF90_CLOBBER, ncid )

    ierr = NF90_PUT_ATT( ncid, NF90_GLOBAL, "purpose", "Fortran Workshop" )
    ierr = NF90_PUT_ATT( ncid, NF90_GLOBAL, "creator", "Denis Sergeev" )
    ierr = NF90_PUT_ATT( ncid, NF90_GLOBAL, "institution", "University of East Anglia" )

    ierr = NF90_DEF_DIM( ncid, "x", nx, x_dimid )
    ierr = NF90_DEF_DIM( ncid, "t", nt, t_dimid )

    ierr = NF90_DEF_VAR( ncid, "x-range", NF90_FLOAT, x_dimid, x_id )
    ierr = NF90_PUT_ATT( ncid, x_id, "units", "metres" )
    ierr = NF90_DEF_VAR( ncid, "t-range", NF90_FLOAT, t_dimid, t_id )
    ierr = NF90_PUT_ATT( ncid, t_id, "units", "seconds" )
    ierr = NF90_DEF_VAR( ncid, "solution", NF90_FLOAT, [ x_dimid, t_dimid ], varid )
    ierr = NF90_PUT_ATT( ncid, varid, "units", "degC" )

    ierr = NF90_ENDDEF( ncid )

    ierr = NF90_PUT_VAR( ncid, x_id, x(:) )
    ierr = NF90_PUT_VAR( ncid, t_id, t(:) )
    ierr = NF90_PUT_VAR( ncid, varid, table(:, :) )

    ierr = NF90_CLOSE( ncid )

  end subroutine


  subroutine r8vec_linspace(a_first, a_last, a)
    implicit none

    real (kind=DP), intent(in)  :: a_first
    real (kind=DP), intent(in)  :: a_last
    real (kind=DP), dimension(:), intent(out) :: a
    integer(kind=SI) :: i
    integer(kind=SI) :: n

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
    integer(kind=SI) :: output_unit_id
    integer(kind=SI) :: j

    output_unit_id = 11
    open (unit=output_unit_id, file=output_filename, status='replace')

    do j = 1, size(x, 1)
      write (output_unit_id, '(2x,g24.16)') x(j)
    end do

    close (unit=output_unit_id)
  end subroutine

end module IO_mod
