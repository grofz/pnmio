! This driver tests "writing the PPM files" and illustrates how the module
! pnmio_module could be used
!
  program test
  use pnmio_module
  implicit none

  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: NX = 320, NY = 240

  real(DP)                  :: uu(NX,NY), rmin, rmax, cpubeg, cpuend
  integer, dimension(NX,NY) :: rr, gg, bb

  real(DP), parameter :: vmax = 100.0, alfa = 0.1
  integer :: ncycles, i

  interface mirror
    subroutine mirror(aa)
      integer, parameter :: DP = kind(1.0d0)
      real(DP), intent(inout) :: aa(:,:)
    end subroutine mirror
  end interface mirror
  interface laplace
    subroutine laplace(aa,alfa)
      integer, parameter :: DP = kind(1.0d0)
      real(DP), intent(inout) :: aa(:,:)
      real(DP), intent(in) :: alfa
    end subroutine laplace
  end interface laplace

  ! initialize array "uu"
  ! "uu" contains random number in the range "vmax"-"2*vmax" (100-200)
  ! uu(1,:), uu(nx,:) are kept constant (north/south)
  ! uu(:,1), uu(:,ny) will be mirrored (periodic BC) (west/east)
  call random_number(uu)
  uu = vmax + uu * vmax         
  uu(1,:) = vmax
  uu(NX,:)= 2.0*vmax

  ! create file "test0.ppm" (initial state)
  ! blue color corresponds to the lowest value in "uu"
  ! red color corresponds to the largest value in "uu"
  call assign_colormap(uu,rr,gg,bb,255)
  call writeppm('test0.ppm',rr,gg,bb)

  ! do averaging
  write(*,'(a)',advance='no') 'cycles = ?:'
  read(*,*) ncycles
  call cpu_time(cpubeg)
  do i = 1, ncycles
    call mirror(uu)
    call laplace(uu, alfa)
  enddo
  call cpu_time(cpuend)
  write(*,'(a,f8.3)') 'Time loop executed = ', cpuend-cpubeg
  
  ! create file "test.ppm" (after "ncycles" steps)
  ! blue color corresponds to value 100.0 (vmax)
  ! red color corresponds to value 200.0 (2*vmax)
  ! values lower than 100.0 will be white 
  ! values higher than 200.0 will be black
  rmin = vmax                
  rmax = 2*vmax               

  ! just show artificial white/black rectangles in the image
  uu(1:10,1:10) = vmax-1.0     
  uu(10:20,1:10) = 2*vmax+1.0  

  call assign_colormap(uu,rr,gg,bb,255,rmin,rmax)
  call writeppm('test.ppm',rr,gg,bb)

  end program test



  subroutine mirror(aa)
  ! make periodic boundary conditions in "y" direction (west/east)
  ! do nothing in "x" direction
  integer, parameter :: DP = kind(1.0d0)
  real(DP), intent(inout) :: aa(:,:)
  integer :: nx, ny

  nx = size(aa,dim=1)
  ny = size(aa,dim=2)
  aa(:,1)  = aa(:,ny-1)
  aa(:,ny) = aa(:,2)

  end subroutine mirror



  subroutine laplace(aa, alfa)
  ! very simple "Laplacian" calculation
  !
  integer, parameter :: DP = kind(1.0d0)
  real(DP), intent(inout) :: aa(:,:)
  real(DP), intent(in) :: alfa
  integer :: nx, ny

  real(DP) :: tmp(size(aa,dim=1)-2, size(aa,dim=2)-2)

  nx = size(aa, dim=1)
  ny = size(aa, dim=2)

  tmp = -4.0*aa(2:nx-1, 2:ny-1)
  tmp = tmp + aa(1:nx-2, 2:ny-1)
  tmp = tmp + aa(3:nx-0, 2:ny-1)
  tmp = tmp + aa(2:nx-1, 1:ny-2)
  tmp = tmp + aa(2:nx-1, 3:ny-0)

  aa(2:nx-1, 2:ny-1) = aa(2:nx-1, 2:ny-1) + alfa * tmp

  end subroutine laplace
