program lorenz96
  implicit none

  ! for namelist
  ! unit
  integer, parameter :: unm = 11
  ! model
  integer :: nx
  integer :: nt
  double precision :: dt
  double precision :: forcing
  ! option
  integer :: imem
  integer :: pmem
  ! unit
  integer :: ufc
  integer :: ueg

  ! variable
  double precision, allocatable :: f(:,:)
  double precision :: initial
  double precision :: norm
  double precision, allocatable :: ravr(:)

  ! etc
  real(kind=16) :: st, et
  integer :: x, t
  integer :: i, j, k
  double precision :: r

  namelist/model/nx, nt, dt, forcing
  namelist/option/imem, pmem 
  namelist/unit/ufc, ueg
  open(unit=unm, file='namelist')
  read(unm, nml=model)
  read(unm, nml=option)
  read(unm, nml=unit)
  close(unm)

  write(*,*) 'this program calculate Lorenz96 model integration'
  write(*,*) 'grid number      = ',nx
  write(*,*) 'integration time = ',nt
  write(*,*) 'time step        = ',dt
  write(*,*) 'forcing          = ',forcing
  write(*,*) 'finite difference scheme = runge-kutta'
  write(*,*) '                  order  = 4th'

  allocate(f(nx+4, nt+1))
  allocate(ravr(nt))
  ravr(1:nt) = 0.0d0

  ! integration start
  call cpu_time(st)

  ! loop of initial condition
  do i = 1, imem
     initial = forcing * i / imem
     f(:,:) = 0.0d0
     f(1:nx,1) = initial
     ! loop of perturbation
     do j = 1, pmem
        call add_random_ptb(nx, f(1:nx,1), norm)
        ! loop of time
        do t = 1, nt
           call set_boundary_condition(nx, f(-1:nx+2,t))
           call runge_kutta_4(nx, dt, forcing, f(-1:nx+2,t), f(-1:nx+2,t+1))
           call calc_error_growth(nx, initial, f(1:nx,t+1), norm, r)
           ravr(t) = ravr(t) + r
        end do
     end do
  end do

  ! error growth
  ravr(1:nt) = sqrt( ravr(1:nt) / imem / pmem )
  ravr(1:nt) = log(ravr(1:nt)) / log(2.0d0)

  ! integration end
  call cpu_time(et)

  ! output
  ! forecast
  call output(ufc, nx, nt+1, f(1:nx,1:nt+1))
  call make_ctl_file(ufc, nx, nt+1)
  ! error gwowth
  call output(ueg, 1, nt, ravr(1:nt))
  call make_ctl_file(ueg, 1, nt)

  write(*,*) 'time for calculation = ',et-st,'seconds'

contains
  subroutine add_random_ptb(nx, f, norm)
    implicit none
    integer, intent(in)          :: nx
    double precision, dimension(1:nx), intent(inout) :: f
    double precision, intent(out) :: norm
    integer :: c, seedsize
    integer, allocatable :: seed(:)
    double precision, dimension(1:nx) :: ptb
    call system_clock(count=c)
    call random_seed(size=seedsize)
    allocate(seed(seedsize))
    seed(:) = c
    call random_seed(put=seed(:))
    call random_number(ptb(:))
    ptb(1:nx) = (-1.0d0 + ptb(1:nx) * 2.0d0) * f(1:nx) * 0.01d0
    f(1:nx) = f(1:nx) + ptb(1:nx)
    ptb(1:nx) = ptb(1:nx) * ptb(1:nx)
    norm = sqrt(sum(ptb(1:nx)))
  end subroutine add_random_ptb

  subroutine set_boundary_condition(nx, f)
    implicit none
    integer, intent(in) :: nx
    double precision, dimension(-1:nx+2), intent(inout) :: f
    integer :: i
    f(-1) = f(nx-1)
    f(0) = f(nx)
    f(nx+1) = f(1)
    f(nx+2) = f(2)
  end subroutine set_boundary_condition

  subroutine runge_kutta_4(nx, dt, forcing, f, g)
    implicit none
    integer, intent(in) :: nx
    double precision,    intent(in) :: dt
    double precision,    intent(in) :: forcing
    double precision, dimension(-1:nx+2), intent(in)  :: f !(t=t)
    double precision, dimension(-1:nx+2), intent(out) :: g !(t=t+1)
    double precision, dimension(1:nx)    :: k1, k2, k3, k4
    double precision, dimension(-1:nx+2) :: w
    k1(1:nx) = ( f(2:nx+1)-f(-1:nx-2) ) * f(0:nx-1) - f(1:nx) + forcing
    w(1:nx) = f(1:nx) + 0.5d0 * dt * k1(1:nx)
    call set_boundary_condition(nx, w(-1:nx+2))
    k2(1:nx) = ( w(2:nx+1)-w(-1:nx-2) ) * w(0:nx-1) - w(1:nx) + forcing
    w(1:nx) = f(1:nx) + 0.5d0 * dt * k2(1:nx)
    call set_boundary_condition(nx, w(-1:nx+2))
    k3(1:nx) = ( w(2:nx+1)-w(-1:nx-2) ) * w(0:nx-1) - w(1:nx) + forcing
    w(1:nx) = f(1:nx) + dt * k3(1:nx)
    call set_boundary_condition(nx, w(-1:nx+2))
    k4(1:nx) = ( w(2:nx+1)-w(-1:nx-2) ) * w(0:nx-1) - w(1:nx) + forcing
    g(1:nx) = f(1:nx) + ( k1(1:nx) + k2(1:nx)*2.0d0 + k3(1:nx)*2.0d0 + k4(1:nx) ) * dt / 6.0d0
  end subroutine runge_kutta_4

  subroutine calc_error_growth(nx, initial, f, norm, r)
    integer, intent(in) :: nx
    double precision, intent(in)  :: initial
    double precision, dimension(1:nx), intent(in) :: f !(t=t)
    double precision, intent(in)  :: norm
    double precision, intent(out) :: r
    double precision :: work(nx)
    work(1:nx) = ( f(1:nx) - initial) * ( f(1:nx) - initial )
    r = sum(work(1:nx)) / (norm * norm)
  end subroutine calc_error_growth

  subroutine output(unit, nx, nt, f)
    integer, intent(in) :: unit
    integer, intent(in) :: nx
    integer, intent(in) :: nt
    double precision, dimension(1:nx,1:nt), intent(in) :: f
    character(2) :: u
    write(u,'(i2)') unit
    open(unit=unit,file='out'//u//'.grd', access='direct',form='unformatted',convert='big_endian',recl=nx*nt*4)
    write(unit,rec=1) real(f(1:nx,1:nt))
    close(unit)
  end subroutine output

  subroutine make_ctl_file(unit, nx, nt)
    implicit none
    integer, intent(in) :: unit
    integer, intent(in) :: nx
    integer, intent(in) :: nt
    character(2) :: u
    write(u,'(i2)') unit
    open(unit=31,file='out'//u//'.ctl')
    write(31,*) 'dset ^out',u,'.grd'
    write(31,*)	'options big_endian'
    write(31,*)	'undef -9.99E+33'
    write(31,*)	'xdef ',nx,' linear 1 1'
    write(31,*)	'ydef 1 levels 0'
    write(31,*) 'zdef 1 levels 1000'
    write(31,*) 'tdef ',nt,' linear 00z01jan0000 72mn'
    write(31,*) 'vars 1'
    write(31,*) 'f 0 99 variable'
    write(31,*) 'endvars'
    close(31)
  end subroutine make_ctl_file

end program lorenz96
