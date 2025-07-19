program example_colormap
  use pnmio_module
  use iso_fortran_env, only : DP=>real64
  implicit none

  integer, parameter :: WIDTH = 300, HEIGHT = 99
  integer, parameter :: CM_SIZE = 11 ! 9+2 colours 

  integer, dimension(WIDTH, HEIGHT) :: rr, gg, bb
  real(DP) :: u(WIDTH, HEIGHT)
  integer :: i

  do i=1, size(u,1)
    ! linspace between 0.0 and 1.0
    u(i,:) = real(i-1)/real(size(u,1)-1)
  end do

  ! artificial rectangles illustrating out-of-range values
  u(1:WIDTH/10, 1:HEIGHT/3) = 2.0_DP
  u(9*WIDTH/10+1:WIDTH, 2*HEIGHT/3+1:HEIGHT) = -1.0_DP

  call assign_colormap(u, rr, gg, bb, CM_SIZE, 0.0_DP, 1.0_DP, CM_RAINBOW)
  call writeppm('test_colormap1.ppm', rr, gg, bb)
  call writeppm('test_colormap1_plain.ppm', rr, gg, bb, is_plain=.true.)

  call assign_colormap(u, rr, gg, bb, CM_SIZE, 0.0_DP, 1.0_DP, CM_VIRIDIS)
  call writeppm('test_colormap2.ppm', rr, gg, bb)
 !call writeppm('test_colormap2_plain.ppm', rr, gg, bb, is_plain=.true.)

end program example_colormap
