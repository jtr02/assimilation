program main
  implicit none
  integer, parameter :: nx = 40
  integer, parameter :: nt = 15001
  real   , parameter :: pi = atan(1.0) * 4.0
  integer :: i, j, k
  real    :: f(nx, nt)
  real    :: z(nx), u1(nx), u2(nx)

  open(unit=21,file='/mnt/d/Data/Assimilation/Result/out51.grd', form='unformatted', &
       & access='direct', convert='big_endian', recl=nx*nt*4)
  read(21,rec=1) f(1:nx, 1:nt)
  close(21)

  ! observation
  open(unit=22,file='/mnt/d/Data/Assimilation/True/true.grd', form='unformatted', &
       & access='direct', convert='big_endian', recl=nx*4)

  open(unit=23,file='/mnt/d/Data/Assimilation/Obs/obs.grd', form='unformatted', &
       & access='direct', convert='big_endian', recl=nx*4)

  j = 0
  do i = 7506, nt, 5
     j = j + 1
     call random_number(u1(1:nx))
     call random_number(u2(1:nx))
     u1(1:nx) = 1.0 - u1(1:nx)
     u2(1:nx) = 1.0 - u2(1:nx)
     z(1:nx) = sqrt( -2.0 * log(u1(1:nx)) ) * cos( 2.0 * pi * u2(1:nx) )
     write(22,rec=j) f(1:nx, i)
     write(23,rec=j) f(1:nx, i) + z(1:nx)
  end do

  close(22)
  close(23)

end program main
