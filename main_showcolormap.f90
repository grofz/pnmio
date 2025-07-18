program example_colormap
  use pnmio_module
  use iso_fortran_env, only : DP=>real64
  implicit none

  integer, parameter :: WIDTH = 300, HEIGHT = 99
  integer, parameter :: COLORMAP_SIZE = 10

  integer, dimension(WIDTH, HEIGHT) :: rr, gg, bb
  real(DP) :: u(WIDTH, HEIGHT)
  integer :: i

  do i=1, size(u,1)
    ! linspace between 0.0 and 1.0
    u(i,:) = real(i-1)/real(size(u,1)-1)
  end do

  ! artificial rectangles
  u(1:WIDTH/10, 1:HEIGHT/3) = 2.0_DP
  u(9*WIDTH/10+1:WIDTH, 2*HEIGHT/3+1:HEIGHT) = -1.0_DP

  call assign_colormap(u, rr, gg, bb, COLORMAP_SIZE, 0.0_DP, 1.0_DP)
print *, 'U'
print '(17(f5.3,1x))', u(:, HEIGHT/2)
print *, 'Red'
print '(17(i3,1x))', rr(:,HEIGHT/2)
print *, 'Green'
print '(17(i3,1x))', gg(:,HEIGHT/2)
print *, 'Blue'
print '(17(i3,1x))', bb(:,HEIGHT/2)
  call writeppm('test_colormap.ppm', rr, gg, bb)
  call writeppm('test_colormap_plain.ppm', rr, gg, bb, is_plain=.true.)

end program example_colormap
