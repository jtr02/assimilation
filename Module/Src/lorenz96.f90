program lorenz96
  implicit none

  ! for integration
  integer, parameter :: nx = 40
  integer, parameter :: nt = 2000
  real, parameter    :: dt = 0.01d0
  real, parameter    :: forcing = 2.0d0

  ! finite difference scheme
  ! forward euler  : 1
  ! leapflog       : 2
  ! adams-bashford : 3
  ! runge-kutta    : 4
  ! matsuno        : 5
  ! semi-implicit  : 6
  integer, parameter :: fds = 2

  ! variable
  real, allocatable :: f(:, :)

  ! etc
  real(kind=16) :: st, et
  integer :: x, t
  integer :: i, j, k

  write(*,*) 'this program calculate Lorenz96 model integration'
  write(*,*) 'grid number =      ',nx
  write(*,*) 'integration time = ',nt
  write(*,*) 'time step   =      ',dt
  write(*,*) 'forcing     =      ',forcing
  write(*,*) 'finite difference scheme = ',fds

  ! input initial condition
  allocate(f(nx+2, nt+1))
!  open(unit=10, file="initial.grd", access='direct',form='unformatted',convert='big_endian',recl=nx*4)
!  read(10,rec=1) f(1:nx, 1)
!  close(10)

  ! test initial
  f(:,:)=0.0
  f(nx/2, 1) = 5

  call cpu_time(st)

  do t = 1, nt

     call set_boundary_condition(nx, f(0:nx+1, t))

     if(t.eq.1) then
        call forward_euler_scheme(nx, dt, forcing, f(0:nx+1, t), f(0:nx+1, t+1))
     else
        if(fds.eq.1) then
           call forward_euler_scheme(nx, dt, forcing, f(0:nx+1, t), f(0:nx+1, t+1))
        else if(fds.eq.2) then
           call centered_scheme(nx, dt, forcing, f(0:nx+1, t), f(0:nx+1, t-1), f(0:nx+1, t+1))
        end if
     end if

  end do

  call cpu_time(et)

  ! output
  open(unit=21,file="output.grd", access='direct',form='unformatted',convert='big_endian',recl=nx*4)
  do t = 1, nt+1
     write(21,rec=t) f(1:nx, t)
  end do
  close(21)
  call make_ctl_file(nx, nt+1)

  write(*,*) 'time for calculation = ',et-st,'seconds'

contains
  subroutine set_boundary_condition(nx, f)
    implicit none
    integer, intent(in) :: nx
    real, dimension(0:nx+1), intent(inout) :: f
    integer :: i
    f(0) = f(nx)
    f(nx+1) = f(1)
  end subroutine set_boundary_condition

  subroutine forward_euler_scheme(nx, dt, forcing, f, g)
    implicit none
    integer, intent(in) :: nx
    real,    intent(in) :: dt
    real,    intent(in) :: forcing
    real, dimension(0:nx+1), intent(in)  :: f !(t=t)
    real, dimension(0:nx+1), intent(out) :: g !(t=t+1)
    integer :: x
    do x = 1, nx
       g(x) = f(x) + ( ( f(x+1)-f(x-2) ) * f(x-1)  - f(x) + forcing ) * dt
    end do
  end subroutine forward_euler_scheme

  subroutine centered_scheme(nx, dt, forcing, f0, f1, g)
    implicit none
    integer, intent(in) :: nx
    real,    intent(in) :: dt
    real,    intent(in) :: forcing
    real, dimension(0:nx+1), intent(in)  :: f0 !(t=t)
    real, dimension(0:nx+1), intent(in)  :: f1 !(t=t-1)
    real, dimension(0:nx+1), intent(out) :: g  !(t=t+1)
    integer :: x
    do x = 1, nx
       g(x) = f1(x) + ( ( f0(x+1)-f0(x-2) ) * f0(x-1)  - f0(x) + forcing ) * 2.0d0 * dt
    end do
  end subroutine centered_scheme

  subroutine make_ctl_file(nx, nta)
    implicit none
    integer, intent(in) :: nx
    integer, intent(in) :: nta
    open(unit=31,file="out.ctl")
    write(31,*) 'dset ^output.grd'
    write(31,*)	'options big_endian'
    write(31,*)	'undef -9.99E+33'
    write(31,*)	'xdef ',nx,' linear 1 1'
    write(31,*)	'ydef 1 levels 0'
    write(31,*) 'zdef 1 levels 1000'
    write(31,*) 'tdef ',nta,' linear 00z27jun1990 1mn'
    write(31,*) 'vars 1'
    write(31,*) 'f 0 99 variable'
    write(31,*) 'endvars'
    close(31)
  end subroutine make_ctl_file

end program lorenz96
