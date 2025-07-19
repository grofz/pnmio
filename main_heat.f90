!
! This driver tests "writing the PPM files" and illustrates how the module
! pnmio_module could be used
!
! It simulates heat conduction on a 2D plate, starting with randomized
! temperature. Euler method is used to solve heat conduction PDE.
!
! Program will prompt user for the number of steps.
!
! Two PPM files are made
! "test0.ppm" as the initial state
! "test.ppm" as the state after the given number of steps
! (test.ppm contains small artificial white/black rectangles to illustrate
!  how the values that are outside of the given colormap range are shown)
!
! contact: grofz@vscht.cz (Z. Grof)
! revised: (c) 2025
!
! comment:
! I know it would be a good practice to place procedures
! "mirror" and "laplace" into a module (and avoid writing interface blocks).
! However, as this is a short and simple example, I prefer not to have an
! additional module
!
  program example_heat
    use pnmio_module
    implicit none

    integer, parameter :: DP = kind(1.0d0)
    integer, parameter :: NX = 320, NY = 240
    integer, parameter :: COLORMAP_SIZE = 255

    real(DP)                  :: uu(NX,NY), rmin, rmax, cpubeg, cpuend
    integer, dimension(NX,NY) :: rr, gg, bb

    real(DP), parameter :: vmax = 100.0, alfa = 0.1
    integer :: ncycles, i, idone

    interface mirror
      subroutine mirror(aa)
        import DP
        real(DP), intent(inout) :: aa(:,:)
      end subroutine mirror
    end interface mirror
    interface laplace
      subroutine laplace(aa,alfa0)
        import DP
        real(DP), intent(inout) :: aa(:,:)
        real(DP), intent(in) :: alfa0
      end subroutine laplace
    end interface laplace
!
! initialize array "uu"
! "uu" contains random number in the range "vmax"-"2*vmax" (100-200)
! uu(1,:), uu(nx,:) are kept constant (north/south)
! uu(:,1), uu(:,ny) will be mirrored (periodic BC) (west/east)
!
    call random_number(uu)
    uu = vmax + uu * vmax
    uu(1,:) = vmax
    uu(NX,:)= 2.0*vmax
!
! create file "test0.ppm" (initial state)
! blue color corresponds to the lowest value in "uu"
! red color corresponds to the largest value in "uu"
!
    call assign_colormap(uu, rr, gg, bb, COLORMAP_SIZE)
    call writeppm('test0.ppm', rr, gg, bb)
    ! alternative color map
    call assign_colormap(uu, rr, gg, bb, COLORMAP_SIZE, idcolormap=CM_TURBO)
    call writeppm('test0_alt.ppm', rr, gg, bb)

!
! do averaging (Euler method steps)
!
    write(*,'(a)',advance='no') 'How many steps to run?  '
    read(*,*) ncycles
    call cpu_time(cpubeg)
    idone = -1
    do i = 1, ncycles
      call mirror(uu)
      call laplace(uu, alfa)
      ! progress counter
      if ((100*i)/ncycles > idone) then
        idone = (100*i)/ncycles
        write(*,'(a,i3,a)',advance='no') achar(13)//'working... ',idone,'%'
      end if
    enddo
    write(*,*)
    call cpu_time(cpuend)
    write(*,'(a,f8.3,a)') 'Time elapsed ', cpuend-cpubeg,' s'
!
! create file "test.ppm" (after "ncycles" steps)
! blue color corresponds to value 100.0 (vmax)
! red color corresponds to value 200.0 (2*vmax)
! values lower than 100.0 will be white
! values higher than 200.0 will be black
!
    rmin = vmax
    rmax = 2*vmax

    ! just show artificial white/black rectangles in the image
    uu(1:10,1:10) = vmax-1.0
    uu(10:20,1:10) = 2*vmax+1.0

    call assign_colormap(uu, rr, gg, bb, COLORMAP_SIZE, rmin, rmax)
    call writeppm('test.ppm', rr, gg, bb)
    ! alternative color map
    call assign_colormap(uu, rr, gg, bb, COLORMAP_SIZE, rmin, rmax, idcolormap=CM_TURBO)
    call writeppm('test_alt.ppm', rr, gg, bb)

  end program example_heat


  subroutine mirror(aa)
!
! make periodic boundary conditions in "y" direction (west/east)
! do nothing in "x" direction
!
    integer, parameter :: DP = kind(1.0d0)
    real(DP), intent(inout) :: aa(:,:)
    integer :: nx, ny

    nx = size(aa,dim=1)
    ny = size(aa,dim=2)
    aa(:,1)  = aa(:,ny-1)
    aa(:,ny) = aa(:,2)

  end subroutine mirror


  subroutine laplace(aa, alfa)
!
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
