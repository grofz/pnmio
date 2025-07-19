! ==============================================
! (c) 2024 Z. Grof (UCT Prague) <grofz@vscht.cz>
! ==============================================


!CONTENT
! module pnmio_module
!   interface assign_colormap
!   interface readpnm
!   interface writeppm
!  *interface s2u
!  *subroutine readpnm_2d(filename, aa, ierr)
!  *subroutine readpnm_3d(filename, aa, ierr)
!  *subroutine writepnm_3d(filename, aa, mx, is_plain, ierr)
!   subroutine writepgm(filename, aa, mx, is_plain, ierr)
!  *subroutine writeppm_3x2d(filename, rr, gg, bb, mx, is_plain, ierr)
!  *subroutine writeppm_1x3d(filename, aa, mx, is_plain, ierr)
!   subroutine colormap(red,green,blue,idcolormap)
!  *subroutine colormap_viridis(red, green, blue)
!  *subroutine colormap_rainbow(red, green, blue)
!  *subroutine assign_colormap_i4(aa,rr,gg,bb,n,amin,amax,idcolormap)
!  *subroutine assign_colormap_r8(uu,rr,gg,bb,n,umin,umax,idcolormap)
!  *helper procedures for stream read
!  *conversion from signed to unsigned numbers
!  *data for additional colormaps
!
! (*) private objects

! SHORT DESCRIPTION
! These routines work with 'raw' or 'plain' PNM formats (P1,P2,P3,P4,P5,P6)
!
! Catching errors:
! If an optional argument "ierr" is provided
! - "ierr" on return contains "0" -> no errors
! - nonzero "ierr" on output -> error occurred, messages sent to stderr
! If "ierr" is not provided
! - in the case of error, program stops


! WEB: https://netpbm.sourceforge.net/doc/pnm.html

! -----------------------------------------------------------------------------
  module pnmio_module
! -----------------------------------------------------------------------------
  use iso_fortran_env, only : error_unit, int8, int16, iostat_end, &
      DP=>real64, SP=>real32
  implicit none
  private
  public readpnm, writepgm, writeppm
  public colormap, assign_colormap

  ! colormap catalog
  integer, parameter, public :: CM_RAINBOW = 1
  integer, parameter, public :: CM_VIRIDIS = 2
  integer, parameter, public :: CM_TURBO = 3
  integer, parameter :: CM_DEFAULT = CM_VIRIDIS

  ! default value of optional MAXVAL argument
  integer, parameter :: MX_DEFAULT = 255, MX_DEFAULT_16 = 65535

  integer, parameter :: IOMSG_MAXLEN = 100

  interface assign_colormap
    module procedure assign_colormap_i4
    module procedure assign_colormap_r8
  end interface assign_colormap

  interface readpnm
    module procedure readpnm_2d
    module procedure readpnm_3d
  end interface

  interface writeppm
    module procedure writeppm_3x2d
    module procedure writeppm_1x3d
  end interface writeppm

  interface s2u
    module procedure s2u_8
    module procedure s2u_16
  end interface

  contains

! --------------
! Read PNM files
! --------------

  subroutine readpnm_2d(filename, aa, ierr)
    character(len=*), intent(in) :: filename
    integer, allocatable, intent(out) :: aa(:,:)
    integer, intent(out), optional :: ierr
!
! This is for PBM and PGM only. Wrapper for readpnm_3d.
!
    integer, allocatable :: atmp(:,:,:)
    integer :: ierr0

    BLK: block
      call readpnm_3d(filename, atmp, ierr0)
      if (ierr0 /= 0) exit BLK
      if (size(atmp, dim=1) /= 1) then
        ierr0 = -2
        write(error_unit, '(a)') 'readpnm error: rank-3 array required for PPM'
        exit BLK
      end if
      allocate(aa(size(atmp,dim=2), size(atmp,dim=3)))
      aa = atmp(1,:,:)
      if (present(ierr)) ierr = 0
      return
    end block BLK

    ! error occurred
    if (present(ierr)) then
      ierr = ierr0 ! user is catching the error
    else
      error stop 'readpnm error - see above'
    end if
  end subroutine readpnm_2d


  subroutine readpnm_3d(filename, aa, ierr)
    character(len=*), intent(in) :: filename
    integer, allocatable, intent(out) :: aa(:,:,:)
    integer, intent(out), optional :: ierr
!
! Read data from PBM, PGM or PPM file
!   Values in array are ordered as
!   - RGB components making a pixel (first index is component id)
!   - rows containing W pixels (second index is column id)
!   - raster containing H rows (third index is row id)
!
    integer :: p, w, h, mx, ierr0, ios, fid, fsize, fpos
    character(len=IOMSG_MAXLEN) :: iomsg
    logical :: file_exist
    integer(int8), allocatable :: raster_1b(:)
    integer(int16), allocatable :: raster_2b(:)

    BLK: block
      ! open file
      inquire(file=filename, exist=file_exist)
      if (.not. file_exist) then
        write(error_unit,'(a)') 'file "'//trim(filename)//'" does not exist'
        exit BLK
      end if
      open(newunit=fid, file = filename, status = 'old', access='stream', &
          form='unformatted', iostat=ios, iomsg=iomsg)
      if (ios /= 0) then
        write(error_unit,'(a)') 'opening file iomsg"'//trim(iomsg)
        exit BLK
      end if

      ! read header
      call consume_header(fid, p, w, h, mx, ierr0)
      if (ierr0 /= 0) exit BLK

      if (mod(p,3)==0) then
        ! PPM
        allocate(aa(3,w,h))
        if (mx > 255) then
          allocate(raster_2b(w*h*3))
        else
          allocate(raster_1b(w*h*3))
        end if
      else
        ! PGM or PBM
        allocate(aa(1,w,h))
        if (mx > 255) then
          allocate(raster_2b(w*h))
        else
          allocate(raster_1b(w*h))
        end if
      end if

      if (mx > 255) then
        call consume_raster(fid, p, w, ierr0, r2b=raster_2b)
        if (ierr0 /= 0) exit BLK
        aa = reshape(s2u(raster_2b), shape(aa))
      else
        call consume_raster(fid, p, w, ierr0, r1b=raster_1b)
        if (ierr0 /= 0) exit BLK
        aa = reshape(s2u(raster_1b), shape(aa))
      end if

      ! P4, P5 and P6 can contain more images in single file
      if ((p-1)/3==1) then
        inquire(fid, size=fsize, pos=fpos)
        if (fpos-1 < fsize) then
          print '("File size ",i0,". File position ",i0)', fsize, fpos
          write(error_unit,'("Warning readpnm: file contains more images")')
          ! TODO ignoring for now
        end if
      end if

      ! verify valid values have been read
      if (any(aa<0) .or. any(aa>mx)) then
        write(error_unit,'("the file contains values out of range (0,mx)")')
        ierr0 = -2
        exit BLK
      end if

      if (present(ierr)) ierr = 0
      return
    end block BLK

    ! error occurred
    if (present(ierr)) then
      ierr = ierr0 ! user must catch this error
    else
      error stop 'readpnm error - see above'
    end if

  end subroutine readpnm_3d


! ---------------
! Write PNM files
! ---------------

  subroutine writepnm_3d(filename, aa, mx, is_plain, ierr)
    character(len=*), intent(in) :: filename
    integer, intent(in) :: aa(:,:,:)
    integer, intent(in) :: mx
    logical, intent(in), optional :: is_plain
    integer, intent(out), optional :: ierr
!
! PPM - first dimension of "aa" is 3, mx > 1   (P3/P6)
! PGM - first dimension of "aa" is 1, mx > 1   (P2/P5)
! PBM - first dimension of "aa" is 1, mx == 1  (P1/P4)
!
! default is raw (unless is_plain==.T. is present)
!
    integer :: p, w, h, ierr0, ios, fid
    character(len=IOMSG_MAXLEN) :: iomsg
    character(len=1) :: magick

    w = size(aa, dim=2)
    h = size(aa, dim=3)
    ierr0 = -1

    BLK: block
      ! verify data are in the range (0, mx)
      if (any(aa<0) .or. any(aa>mx)) then
        write(error_unit,'("All values must be in range (0,mx)")')
        exit BLK
      end if

      ! identify P
      if (size(aa,dim=1)==1) then
        if (mx > 1) then
          p = 5 ! PGM
        else if (mx == 1) then
          p = 4 ! PBM
        else
          write(error_unit,'("mx must be grater than 0")')
          exit BLK
        end if
      else if (size(aa,dim=1)==3) then
        if (mx > 1) then
          p = 6 ! PPM
        else
          write(error_unit,'("mx must be grater than 1 for PPM")')
          exit BLK
        end if
      else
        write(error_unit,'("first extent of aa must be 1 or 3")')
        exit BLK
      end if
      if (present(is_plain)) then
        if (is_plain) p = p-3
      end if

      ! open file
      open(newunit=fid, file=filename, status='replace', access='stream', &
          form='unformatted', iostat=ios, iomsg=iomsg)
      if (ios /= 0) then
        write(error_unit,'(a)') 'opening file for writting iomsg: '//trim(iomsg)
        exit BLK
      end if

      ! write header
      HEADER: block
        ! Note to myself: I preffer to do formatting via internal
        ! character storage instead of "formatted stream"
        character(len=100) buf
        write(magick,'(i1)') p
        write(buf,'(a)') 'P'//magick
        write(fid) trim(buf)//new_line(buf)
        write(buf,'(a)') '# created by pnmio fortran library'
        write(fid) trim(buf)//new_line(buf)
        write(buf,'(i0,1x,i0)') w, h
        write(fid) trim(buf)//new_line(buf)

        if (mod(p,3)/=1) then
          write(buf,'(i0)') mx
          write(fid) trim(buf)//new_line(buf)
        end if
      end block HEADER

      ! write raster
      if (p > 3) then
        ! raw format
        BINARY: block
          integer(int8), allocatable :: raster_1b(:)
          integer(int16), allocatable :: raster_2b(:)

          ! write raster (1B or 2B)
          if (mx > 255) then
            if (p==6) then
              allocate(raster_2b(3*w*h)) ! PPM (16 bit)
            else
              allocate(raster_2b(w*h))   ! PGM (16 bit)
            end if
            write(fid, iostat=ios, iomsg=iomsg) reshape(u2s_16(aa), shape(raster_2b))
          elseif (mx == 1) then
            PBM: block ! PBM
              integer :: wc, irow, icol, j, k
              integer(int8), allocatable :: r8(:)
              wc = w/8
              if (mod(w,8)/=0) wc=wc+1
              allocate(r8(wc))
              do irow=1, h
                r8 = 0
                icol = 0
                JLOOP: do j=1, wc
                  do k=bit_size(r8)-1,0,-1
                    icol = icol + 1
                    if (aa(1,icol,irow)==1) r8(j) = ibset(r8(j),k)
                    if (icol==w) exit JLOOP
                  end do
                end do JLOOP
                write(fid, iostat=ios, iomsg=iomsg) r8
              end do
            end block PBM

          else
            if (p==6) then
              allocate(raster_1b(3*w*h)) ! PPM
            else
              allocate(raster_1b(w*h))   ! PGM
            end if
            write(fid, iostat=ios, iomsg=iomsg) reshape(u2s_8(aa), shape(raster_1b))
          end if
          if (ios /= 0) then
            write(error_unit,'(a)') 'writing raster iomsg: '//trim(iomsg)
            exit BLK
          end if
        end block BINARY

      else
        ! plain format
        ASCII: block
          ! as of netpnm format specification, maximum number of characters
          ! per line must be less than 70
          character(len=*), parameter :: fmt_short='(17(i3,1x))'
          character(len=*), parameter :: fmt_long='(11(i5,1x))'

          ! re-open as "formatted"
          close(fid, iostat=ios, iomsg=iomsg)
          if (ios /= 0) then
            write(error_unit,'(a)') 'closing to reopen file iomsg: '//trim(iomsg)
            exit BLK
          end if
          open(newunit=fid, file=filename, status='old', access='stream', &
              form='formatted', position='append', iostat=ios, iomsg=iomsg)
          if (ios /= 0) then
            write(error_unit,'(a)') 'reopening file iomsg: '//trim(iomsg)
            exit BLK
          end if
          if (mx > 255) then
            write(fid,fmt=fmt_long,iostat=ios,iomsg=iomsg) aa
          else
            write(fid,fmt=fmt_short,iostat=ios,iomsg=iomsg) aa
          end if
          if (ios /= 0) then
            write(error_unit,'(a)') 'writing ASCII raster iomsg: '//trim(iomsg)
            exit BLK
          end if
        end block ASCII
      end if

      ! close file
      close(fid, iostat=ios, iomsg=iomsg)
      if (ios /= 0) then
        write(error_unit,'(a)') 'closing file iomsg: '//trim(iomsg)
        exit BLK
      end if

      if (present(ierr)) ierr = 0
      return
    end block BLK

    ! error occurred
    if (present(ierr)) then
      ierr = ierr0 ! user must catch this error
    else
      error stop 'writepnm error - see above'
    end if

  end subroutine writepnm_3d


  subroutine writepgm(filename, aa, mx, is_plain, ierr)
    character(len=*), intent(in)   :: filename
    integer, intent(in)            :: aa(:,:)
    integer, intent(in), optional  :: mx
    logical, intent(in), optional  :: is_plain
    integer, intent(out), optional :: ierr
!
! Write PGM or PBM file. Public wrapper to "writepnm_3d"
! PBM - if "mx" is present and "mx==1"
! PGM - otherwise
!
    integer, allocatable :: atmp(:,:,:)

    allocate(atmp(1, size(aa,1), size(aa,2)))
    atmp(1,:,:) = aa
    call writepnm_3d(filename, atmp, guess_mx(atmp,mx), is_plain, ierr)

  end subroutine writepgm


  subroutine writeppm_3x2d(filename, rr, gg, bb, mx, is_plain, ierr)
    character(len=*), intent(in)   :: filename
    integer, intent(in)            :: rr(:,:), gg(:,:), bb(:,:)
    integer, intent(in), optional  :: mx
    logical, intent(in), optional  :: is_plain
    integer, intent(out), optional :: ierr
!
! Write PPM file. Public wrapper to "writepnm_3d"
!
    integer, allocatable :: atmp(:,:,:)

    if (any(shape(rr) /= shape(gg)) .or. any(shape(rr) /= shape(bb))) then
      if (present(ierr)) then
        ierr = -3
        return
      else
        error stop 'writeppm: arrays rr, gg and bb must have same shape'
      end if
    end if

    allocate(atmp(3, size(rr,1), size(rr,2)))
    atmp(1,:,:) = rr
    atmp(2,:,:) = gg
    atmp(3,:,:) = bb
    call writepnm_3d(filename, atmp, guess_mx(atmp,mx), is_plain, ierr)

  end subroutine writeppm_3x2d


  subroutine writeppm_1x3d(filename, aa, mx, is_plain, ierr)
    character(len=*), intent(in)   :: filename
    integer, intent(in)            :: aa(:,:,:)
    integer, intent(in), optional  :: mx
    logical, intent(in), optional  :: is_plain
    integer, intent(out), optional :: ierr
!
! Write PPM file. Public wrapper to "writepnm_3d"
!
    if (size(aa,1)/=3) then
      if (present(ierr)) then
        ierr = -3
        return
      else
        error stop 'writeppm: first extent of aa must be 3'
      end if
    end if

    call writepnm_3d(filename, aa, guess_mx(aa,mx), is_plain, ierr)

  end subroutine writeppm_1x3d


! ---------
! Colormaps
! ---------

  subroutine colormap(red, green, blue, idcolormap)
    integer, intent(out) :: red(:), green(:), blue(:)
    integer, intent(in), optional :: idcolormap
!
! Fill colour arrays with values for a chosen colormap.
! The first / last array item corresponds to datapoint with a value
! below minimum / above maximum range and are assigned to
! the black / the white colour.
! The size of arrays determines the number of distinct colours.
!
    integer :: idcolormap0

    ! verify size of arrays match
    if (size(red)/=size(green) .or. size(red)/=size(blue)) &
        error stop 'colormap - array size do not match'

    idcolormap0 = CM_DEFAULT
    if (present(idcolormap)) idcolormap0 = idcolormap
    select case(idcolormap0)
    case(CM_RAINBOW)
      call colormap_rainbow(red, green, blue)
    case(CM_VIRIDIS)
      call colormap_viridis(red, green, blue)
    case(CM_TURBO)
      call colormap_turbo(red, green, blue)
    case default
      error stop 'assign_colormap - unknown id of colormap'
    end select

  end subroutine colormap


  subroutine colormap_viridis(red, green, blue)
    integer, intent(out) :: red(0:), green(0:), blue(0:)
!
! Generate viridis colormap
!
    integer, parameter :: mxval = 255
    real, allocatable :: dat(:,:)
    integer :: i, n, ilow, ndat
    real :: x, xlow, xhigh, s

    allocate(dat(0,0)) ! just to avoid -Wunitialized warning
    dat = viridis_data()
    ndat = size(dat,2)
    n = size(red)-2 ! red(0) and red(n+1) are reserved for white/black

    red(0) = mxval
    green(0) = mxval
    blue(0) = mxval
    do i=1, n
      x = real(i-1)/real(n-1)
      ilow = min(int(x*real(ndat-1))+1, ndat-1) ! lower index
      xlow = real(ilow-1)/real(ndat-1)
      xhigh = real(ilow)/real(ndat-1)
      s = (x-xlow) / (xhigh-xlow) ! "s" in range <0, 1>
      red(i) = int(real(mxval) * (dat(1,ilow)*s + dat(1,ilow+1)*(1.0-s)))
      green(i) = int(real(mxval) * (dat(2,ilow)*s + dat(2,ilow+1)*(1.0-s)))
      blue(i) = int(real(mxval) * (dat(3,ilow)*s + dat(3,ilow+1)*(1.0-s)))
    end do
    red(n+1) = 0
    green(n+1) = 0
    blue(n+1) = 0

  end subroutine colormap_viridis


  subroutine colormap_rainbow(red, green, blue)
    integer, intent(out) :: red(:), green(:), blue(:)
!
! Generate jet colormap (blue--cyan--green--yellow--red)
! First and last colours are white and black, respectivelly
!
! https://blog.habrador.com/2023/04/colormaps-overview-code-implementations-rainbow-virids.html
! https://paulbourke.net/miscellaneous/colourspace/
!
    integer, parameter :: mxval = 255
    integer            :: n, i
    real :: x, yr, yb, yg

    n = size(red,dim=1)

    red(1) = mxval
    green(1) = mxval
    blue(1) = mxval

    do i = 2, n-1
      x = real(i-2) / real(n-3) ! x will be in range (0;1)
      yr =  4.0*x - 2.0
      yb = -4.0*x + 2.0
      if(x<0.5) then
        yg =  4.0*x + 0.0
      else
        yg = -4.0*x + 4.0
      endif
      yr = min(1.0, max(0.0, yr))
      yg = min(1.0, max(0.0, yg))
      yb = min(1.0, max(0.0, yb))
      red(i)   = int(yr*real(mxval)) ! values will be between 0..mxval
      green(i) = int(yg*real(mxval))
      blue(i)  = int(yb*real(mxval))
    enddo

    red(n) = 0
    green(n) = 0
    blue(n) = 0

  end subroutine colormap_rainbow


  subroutine colormap_turbo(red, green, blue)
    integer, intent(out) :: red(0:), green(0:), blue(0:)
!
! Generate turbo colormap (improved rainbow)
!
    integer, parameter :: mxval = 255
    integer :: i, n
    real :: x, rgb(3)

    n = size(red)-2 ! red(0) and red(n+1) are reserved for white/black

    red(0) = mxval
    green(0) = mxval
    blue(0) = mxval
    do i=1, n
      x = real(i-1)/real(n-1)
      rgb = turbocm(x)
      red(i) = int(real(mxval) * rgb(1))
      green(i) = int(real(mxval) * rgb(2))
      blue(i) = int(real(mxval) * rgb(3))
    end do
    red(n+1) = 0
    green(n+1) = 0
    blue(n+1) = 0

  end subroutine colormap_turbo


! ===============================
! Colorize integer or real fields
! ===============================

  subroutine assign_colormap_i4(aa, rr, gg, bb, n, amin, amax, idcolormap)
    integer, intent(in)  :: aa(:,:)
    integer, intent(out) :: rr(:,:), gg(:,:), bb(:,:)
    integer, intent(in)  :: n ! number of colormap values 
    integer, intent(in), optional :: amin, amax
    integer, intent(in), optional :: idcolormap

    ! - local vars
    integer :: cmap(n,3)
    real    :: dx, fmin, fmax, f
    integer :: i, j, m

    ! -
    call colormap(cmap(:,1),cmap(:,2),cmap(:,3),idcolormap)

    if (present(amin)) then
      fmin = real(amin)
    else
      fmin = real(minval(aa))
    end if
    if (present(amax)) then
      fmax = real(amax)
    else
      fmax = real(maxval(aa))
    end if

    ! to elements aa==fmin assign 2nd color,
    ! to elements aa==fmax assign (n-1)th color
    ! 1st and nth colors are reserved for out of range elements

    dx = (fmax-fmin)/float(n-3)

    do i=1,size(aa,dim=1)
    do j=1,size(aa,dim=2)
      f = (real(aa(i,j)) - fmin) / dx
      m = int(f) + 2
      m = min(max(2,m),n-1)
      rr(i,j) = cmap(m,1)
      gg(i,j) = cmap(m,2)
      bb(i,j) = cmap(m,3)
    enddo
    enddo

    where(aa < amin)
      rr = cmap(1,1)
      gg = cmap(1,2)
      bb = cmap(1,3)
    endwhere
    where(aa > amax)
      rr = cmap(n,1)
      gg = cmap(n,2)
      bb = cmap(n,3)
    endwhere

  end subroutine assign_colormap_i4


  subroutine assign_colormap_r8(uu, rr, gg, bb, n, umin, umax, idcolormap)
    real(DP), intent(in)  :: uu(:,:)
    integer, intent(out) :: rr(:,:), gg(:,:), bb(:,:)
    integer, intent(in)  :: n
    real(DP), intent(in), optional :: umin, umax
    integer, intent(in), optional :: idcolormap
!
! UU - field with values to colorize
! RR, GG, BB - pixels to be used in "writeppm"
! N - number of values on the colormap axis (no. of colours)
!     (two colours are used for out of range values)
! UMIN, UMAX - set low and high limits 
!              (if not present, min/max values in UU are used)
! IDCOLORMAP - id. of colormap used (if not present, default colormap selected)
!
    integer :: cmap(n,3), aa(size(uu,1),size(uu,2)), i, j
    real(DP) :: umin0, umax0

    if (present(umin)) then
      umin0 = umin
    else
      umin0 = minval(uu)
    endif

    if (present(umax)) then
      umax0 = umax
    else
      umax0 = maxval(uu)
    endif

    call colormap(cmap(:,1), cmap(:,2), cmap(:,3), idcolormap)

    where (uu < umin0)
      aa = 1
    else where (uu > umax0)
      aa = n
    else where
      ! aa is between 2 and n-1
      aa = min(n-1, int((uu-umin0)/(umax0-umin0)*(n-2))+2)
    end where

    do i=1,size(uu,1)
    do j=1,size(uu,2)
      rr(i,j) = cmap(aa(i,j),1)
      gg(i,j) = cmap(aa(i,j),2)
      bb(i,j) = cmap(aa(i,j),3)
    end do
    end do

  end subroutine assign_colormap_r8


! ---------------------------------
! helper procedures for stream read
! ---------------------------------

  pure logical function is_whitespace(ch)
    character(len=1), intent(in) :: ch

    ! Assuming white space is a character with ASCII code between
    ! 9 and 13 (TAB, LF, VT, FF, CR) or a space
    is_whitespace = (iachar(ch)>=9 .and. iachar(ch)<=13) .or. ch==' '
  end function is_whitespace


  pure logical function is_newline(ch)
    character(len=1), intent(in) :: ch

    ! Assuming new line is LF
    is_newline = iachar(ch) == 10 ! iachar(new_line(ch))
  end function is_newline


  subroutine consume_magick(fid, val)
    integer, intent(in) :: fid
    integer, intent(out) :: val ! 1-6 or 0 for an unknown signature

    ! we assume "fid" is opened as "unformatted stream"
    character(len=2) :: ch2
    integer :: ios
    character(len=IOMSG_MAXLEN) :: iomsg

    val = 0
    read(fid, iostat=ios, iomsg=iomsg) ch2
    if (ios /= 0) then
      write(error_unit,'(a)') 'consume_magick iomsg: '//trim(iomsg)
      return
    end if
    select case(ch2)
    case('P1')
      val = 1
    case('P2')
      val = 2
    case('P3')
      val = 3
    case('P4')
      val = 4
    case('P5')
      val = 5
    case('P6')
      val = 6
    end select
  end subroutine consume_magick


  subroutine consume_decimal(fid, val, ierr)
    integer, intent(in) :: fid
    integer, intent(out) :: val, ierr

    ! consume ASCII decimal from the stream
    ! - read and ignore all whitespace characters
    ! - if '#' read, then ignore all characters until CR or LF
    ! we assume "fid" is opened as "unformatted stream"

    integer, parameter :: MODE_SCAN=0, MODE_COMMENT=1, MODE_READ=2
    integer :: ios, pos, mode, dec_len
    character(len=1)   :: ch
    character(len=10)  :: dec
    character(len=IOMSG_MAXLEN) :: iomsg

    ierr = -1
    mode = MODE_SCAN
    dec = ''
    dec_len = 0
    inquire(unit=fid, pos=pos)
    do
      read(fid, iostat=ios, pos=pos, iomsg=iomsg) ch
      pos = pos + 1
      ! if the last digit is at the very end, this is not an error
      if (ios == iostat_end .and. mode==MODE_READ) exit
      if (ios /= 0) then
        write(error_unit,'(a)') 'consume_decimal iomsg: '//trim(iomsg)
        return
      end if

      select case(mode)
      case(MODE_SCAN)
        ! stay in scan mode until non-whitespace character is encountered
        if (.not. is_whitespace(ch)) then
          if (ch=='#') then
            mode = MODE_COMMENT
          else
            mode = MODE_READ
            pos = pos - 1 ! re-read first digit
          end if
        end if
      case(MODE_READ)
        ! copy all digits, then exit
        if (is_whitespace(ch)) exit
        if (dec_len==len(dec)) then
          write(error_unit,'(a)') 'consume_decimal: too many digits'
          return
        end if
        dec_len = dec_len + 1
        dec(dec_len:dec_len) = ch
      case(MODE_COMMENT)
        ! read and ignore all until CR or LF, then switch to scan mode again
        if (iachar(ch)==10 .or. iachar(ch)==13) then
          mode = MODE_SCAN
          pos = pos - 1 ! should not be part of comment
        end if
      case default
        error stop 'consume_decimal: should not be here'
      end select
    end do

    read(dec(1:dec_len),*,iostat=ios, iomsg=iomsg) val
    if (ios /= 0) then
      write(error_unit,*) 'consume_decimal conversion: '//trim(iomsg)
      return
    end if

    ierr = 0 ! "val" read withou any error
  end subroutine consume_decimal


  subroutine consume_header(fid, p, w, h, mxval, ierr)
    integer, intent(in)  :: fid
    integer, intent(out) :: p, w, h, mxval, ierr

    ! read header of an PPM image
    ! assuming "fid" is opened as "unformatted stream"

    HDR: block
      call consume_magick(fid, p)
      if (p==0) exit HDR
      call consume_decimal(fid, w, ierr)
      if (ierr /= 0) exit HDR
      call consume_decimal(fid, h, ierr)
      if (ierr /= 0) exit HDR
      ! "maxval" is not defined for PBM format (P1 or P4)
      if (p/=1 .and. p/=4) then
        call consume_decimal(fid, mxval, ierr)
        if (ierr /= 0) exit HDR
      else
        mxval = 1
      end if
      ierr = 0
      return
    end block HDR

    ! an error occurred
    ierr = -1
    write(error_unit,'(a)') 'error: header is invalid'
  end subroutine consume_header


  subroutine consume_raster(fid, p, w, ierr, r1b, r2b)
    integer, intent(in) :: fid, p, w
    integer, intent(out) :: ierr
    integer(int8), intent(out), optional :: r1b(:)
    integer(int16), intent(out), optional :: r2b(:)

    ! assuming "fid" is opened as "unformatted stream"
    ! "W" is needed for P4 format only
    integer :: ios
    character(len=IOMSG_MAXLEN) :: iomsg
    integer, allocatable :: u(:)

    if (present(r1b) .eqv. present(r2b)) &
      error stop 'consume_raster: only r1b or r2b must be given'
    ierr = -1

    if ((p-1)/3 == 1) then
      ! binary format
      if (p==4) then
        PBM: block ! special case for binary PBM
          integer :: wc, h, i, j, k
          integer(int8), allocatable :: r8(:)
          h = size(r1b)/w
          wc = w/8
          if (mod(w,8)/=0) wc = wc+1
          allocate(r8(h*wc))
          read(fid, iostat=ios, iomsg=iomsg) r8
          r1b = 0
          i = 0
          JLOOP: do j=1, size(r8)
            do k=bit_size(r8)-1,0,-1
              i = i+1
              if (btest(r8(j),k)) r1b(i) = 1
              if (mod(i,w)==0) cycle JLOOP
            end do
            end do JLOOP
        end block PBM
      else if (present(r1b)) then
        read(fid, iostat=ios, iomsg=iomsg) r1b
      else
        read(fid, iostat=ios, iomsg=iomsg) r2b
      end if
      if (ios/=0) then
        write(error_unit, '(a)') 'raster iomsg: '//trim(iomsg)
        return
      end if

    else if ((p-1)/3 == 0) then
      ! ASCII digits
      if (present(r1b)) then
        allocate(u(size(r1b)))
      else
        allocate(u(size(r2b)))
      end if

      ASCII: block
        integer :: i
        do i=1, size(u)
          call consume_decimal(fid, u(i), ierr)
          if (ierr /= 0) then
            write(error_unit, '("raster raw error at posiiton ",i0," out of ",i0)') i, size(u)
            return
          end if
        end do
      end block ASCII

      ! convert to unsigned so it fits to 1B or 2B integers
      if (present(r1b)) then
        r1b = u2s_8(u)
      else
        r2b = u2s_16(u)
      end if

    else
      error stop 'consume_raster: unreachable'
    end if

    ierr = 0
  end subroutine consume_raster


! ------------------------------------------
! Conversion from signed to unsigned numbers
! ------------------------------------------

  elemental function s2u_8(sint) result (uint)
    integer(int8), intent(in) :: sint
    integer :: uint

    if (sint >= 0_int8) then
      uint = int(sint, kind=kind(uint))
    else
      uint = 2*(huge(sint)+1) + int(sint, kind=kind(uint))
    end if
  end function


  elemental function s2u_16(sint) result (uint)
    integer(int16), intent(in) :: sint
    integer :: uint

    if (sint >= 0_int16) then
      uint = int(sint, kind=kind(uint))
    else
      uint = 2*(huge(sint)+1) + int(sint, kind=kind(uint))
    end if
  end function


  elemental function u2s_8(u) result (s)
    integer, intent(in) :: u
    integer(int8) :: s

    if (u < 0) then
      error stop 'u2s_8 input is negative'
    else if (u <= huge(s)) then
      s = int(u, kind=kind(s))
    elseif (u < 2*(huge(s)+1)) then
      s = int(u-2*(huge(s)+1), kind=kind(s))
    else
      error stop 'u2s_8 input too big to convert'
    end if
  end function


  elemental function u2s_16(u) result (s)
    integer, intent(in) :: u
    integer(int16) :: s

    if (u < 0) then
      error stop 'u2s_16 input is negative'
    else if (u <= huge(s)) then
      s = int(u, kind=kind(s))
    elseif (u < 2*(huge(s)+1)) then
      s = int(u-2*(huge(s)+1), kind=kind(s))
    else
      error stop 'u2s_16 input too big to convert'
    end if
  end function


  integer function guess_mx(aa, mx) result(mx0)
    integer, intent(in) :: aa(:,:,:)
    integer, intent(in), optional :: mx
!
! Guess the value of "mx" in the case it is not given
!
    if (present(mx)) then
      mx0 = mx
    else
      if (maxval(aa) > 255) then
        mx0 = MX_DEFAULT_16
      else
        mx0 = MX_DEFAULT
      end if
    end if
  end function


! ------------------
! Data for colormaps
! ------------------
  function viridis_data() result(res)
    real :: res(3,256)
!
! Values taken from Matplotlib
! https://github.com/BIDS/colormap/blob/master/colormaps.py
!
! You can also watch a presentation on Youtube
! https://www.youtube.com/watch?v=xAoljeRJ3lU
!
    real, parameter, dimension(*) :: &
      viridis = [0.267004, 0.004874, 0.329415, &
                 0.268510, 0.009605, 0.335427, &
                 0.269944, 0.014625, 0.341379, &
                 0.271305, 0.019942, 0.347269, &
                 0.272594, 0.025563, 0.353093, &
                 0.273809, 0.031497, 0.358853, &
                 0.274952, 0.037752, 0.364543, &
                 0.276022, 0.044167, 0.370164, &
                 0.277018, 0.050344, 0.375715, &
                 0.277941, 0.056324, 0.381191, &
                 0.278791, 0.062145, 0.386592, &
                 0.279566, 0.067836, 0.391917, &
                 0.280267, 0.073417, 0.397163, &
                 0.280894, 0.078907, 0.402329, &
                 0.281446, 0.084320, 0.407414, &
                 0.281924, 0.089666, 0.412415, &
                 0.282327, 0.094955, 0.417331, &
                 0.282656, 0.100196, 0.422160, &
                 0.282910, 0.105393, 0.426902, &
                 0.283091, 0.110553, 0.431554, &
                 0.283197, 0.115680, 0.436115, &
                 0.283229, 0.120777, 0.440584, &
                 0.283187, 0.125848, 0.444960, &
                 0.283072, 0.130895, 0.449241, &
                 0.282884, 0.135920, 0.453427, &
                 0.282623, 0.140926, 0.457517, &
                 0.282290, 0.145912, 0.461510, &
                 0.281887, 0.150881, 0.465405, &
                 0.281412, 0.155834, 0.469201, &
                 0.280868, 0.160771, 0.472899, &
                 0.280255, 0.165693, 0.476498, &
                 0.279574, 0.170599, 0.479997, &
                 0.278826, 0.175490, 0.483397, &
                 0.278012, 0.180367, 0.486697, &
                 0.277134, 0.185228, 0.489898, &
                 0.276194, 0.190074, 0.493001, &
                 0.275191, 0.194905, 0.496005, &
                 0.274128, 0.199721, 0.498911, &
                 0.273006, 0.204520, 0.501721, &
                 0.271828, 0.209303, 0.504434, &
                 0.270595, 0.214069, 0.507052, &
                 0.269308, 0.218818, 0.509577, &
                 0.267968, 0.223549, 0.512008, &
                 0.266580, 0.228262, 0.514349, &
                 0.265145, 0.232956, 0.516599, &
                 0.263663, 0.237631, 0.518762, &
                 0.262138, 0.242286, 0.520837, &
                 0.260571, 0.246922, 0.522828, &
                 0.258965, 0.251537, 0.524736, &
                 0.257322, 0.256130, 0.526563, &
                 0.255645, 0.260703, 0.528312, &
                 0.253935, 0.265254, 0.529983, &
                 0.252194, 0.269783, 0.531579, &
                 0.250425, 0.274290, 0.533103, &
                 0.248629, 0.278775, 0.534556, &
                 0.246811, 0.283237, 0.535941, &
                 0.244972, 0.287675, 0.537260, &
                 0.243113, 0.292092, 0.538516, &
                 0.241237, 0.296485, 0.539709, &
                 0.239346, 0.300855, 0.540844, &
                 0.237441, 0.305202, 0.541921, &
                 0.235526, 0.309527, 0.542944, &
                 0.233603, 0.313828, 0.543914, &
                 0.231674, 0.318106, 0.544834, &
                 0.229739, 0.322361, 0.545706, &
                 0.227802, 0.326594, 0.546532, &
                 0.225863, 0.330805, 0.547314, &
                 0.223925, 0.334994, 0.548053, &
                 0.221989, 0.339161, 0.548752, &
                 0.220057, 0.343307, 0.549413, &
                 0.218130, 0.347432, 0.550038, &
                 0.216210, 0.351535, 0.550627, &
                 0.214298, 0.355619, 0.551184, &
                 0.212395, 0.359683, 0.551710, &
                 0.210503, 0.363727, 0.552206, &
                 0.208623, 0.367752, 0.552675, &
                 0.206756, 0.371758, 0.553117, &
                 0.204903, 0.375746, 0.553533, &
                 0.203063, 0.379716, 0.553925, &
                 0.201239, 0.383670, 0.554294, &
                 0.199430, 0.387607, 0.554642, &
                 0.197636, 0.391528, 0.554969, &
                 0.195860, 0.395433, 0.555276, &
                 0.194100, 0.399323, 0.555565, &
                 0.192357, 0.403199, 0.555836, &
                 0.190631, 0.407061, 0.556089, &
                 0.188923, 0.410910, 0.556326, &
                 0.187231, 0.414746, 0.556547, &
                 0.185556, 0.418570, 0.556753, &
                 0.183898, 0.422383, 0.556944, &
                 0.182256, 0.426184, 0.557120, &
                 0.180629, 0.429975, 0.557282, &
                 0.179019, 0.433756, 0.557430, &
                 0.177423, 0.437527, 0.557565, &
                 0.175841, 0.441290, 0.557685, &
                 0.174274, 0.445044, 0.557792, &
                 0.172719, 0.448791, 0.557885, &
                 0.171176, 0.452530, 0.557965, &
                 0.169646, 0.456262, 0.558030, &
                 0.168126, 0.459988, 0.558082, &
                 0.166617, 0.463708, 0.558119, &
                 0.165117, 0.467423, 0.558141, &
                 0.163625, 0.471133, 0.558148, &
                 0.162142, 0.474838, 0.558140, &
                 0.160665, 0.478540, 0.558115, &
                 0.159194, 0.482237, 0.558073, &
                 0.157729, 0.485932, 0.558013, &
                 0.156270, 0.489624, 0.557936, &
                 0.154815, 0.493313, 0.557840, &
                 0.153364, 0.497000, 0.557724, &
                 0.151918, 0.500685, 0.557587, &
                 0.150476, 0.504369, 0.557430, &
                 0.149039, 0.508051, 0.557250, &
                 0.147607, 0.511733, 0.557049, &
                 0.146180, 0.515413, 0.556823, &
                 0.144759, 0.519093, 0.556572, &
                 0.143343, 0.522773, 0.556295, &
                 0.141935, 0.526453, 0.555991, &
                 0.140536, 0.530132, 0.555659, &
                 0.139147, 0.533812, 0.555298, &
                 0.137770, 0.537492, 0.554906, &
                 0.136408, 0.541173, 0.554483, &
                 0.135066, 0.544853, 0.554029, &
                 0.133743, 0.548535, 0.553541, &
                 0.132444, 0.552216, 0.553018, &
                 0.131172, 0.555899, 0.552459, &
                 0.129933, 0.559582, 0.551864, &
                 0.128729, 0.563265, 0.551229, &
                 0.127568, 0.566949, 0.550556, &
                 0.126453, 0.570633, 0.549841, &
                 0.125394, 0.574318, 0.549086, &
                 0.124395, 0.578002, 0.548287, &
                 0.123463, 0.581687, 0.547445, &
                 0.122606, 0.585371, 0.546557, &
                 0.121831, 0.589055, 0.545623, &
                 0.121148, 0.592739, 0.544641, &
                 0.120565, 0.596422, 0.543611, &
                 0.120092, 0.600104, 0.542530, &
                 0.119738, 0.603785, 0.541400, &
                 0.119512, 0.607464, 0.540218, &
                 0.119423, 0.611141, 0.538982, &
                 0.119483, 0.614817, 0.537692, &
                 0.119699, 0.618490, 0.536347, &
                 0.120081, 0.622161, 0.534946, &
                 0.120638, 0.625828, 0.533488, &
                 0.121380, 0.629492, 0.531973, &
                 0.122312, 0.633153, 0.530398, &
                 0.123444, 0.636809, 0.528763, &
                 0.124780, 0.640461, 0.527068, &
                 0.126326, 0.644107, 0.525311, &
                 0.128087, 0.647749, 0.523491, &
                 0.130067, 0.651384, 0.521608, &
                 0.132268, 0.655014, 0.519661, &
                 0.134692, 0.658636, 0.517649, &
                 0.137339, 0.662252, 0.515571, &
                 0.140210, 0.665859, 0.513427, &
                 0.143303, 0.669459, 0.511215, &
                 0.146616, 0.673050, 0.508936, &
                 0.150148, 0.676631, 0.506589, &
                 0.153894, 0.680203, 0.504172, &
                 0.157851, 0.683765, 0.501686, &
                 0.162016, 0.687316, 0.499129, &
                 0.166383, 0.690856, 0.496502, &
                 0.170948, 0.694384, 0.493803, &
                 0.175707, 0.697900, 0.491033, &
                 0.180653, 0.701402, 0.488189, &
                 0.185783, 0.704891, 0.485273, &
                 0.191090, 0.708366, 0.482284, &
                 0.196571, 0.711827, 0.479221, &
                 0.202219, 0.715272, 0.476084, &
                 0.208030, 0.718701, 0.472873, &
                 0.214000, 0.722114, 0.469588, &
                 0.220124, 0.725509, 0.466226, &
                 0.226397, 0.728888, 0.462789, &
                 0.232815, 0.732247, 0.459277, &
                 0.239374, 0.735588, 0.455688, &
                 0.246070, 0.738910, 0.452024, &
                 0.252899, 0.742211, 0.448284, &
                 0.259857, 0.745492, 0.444467, &
                 0.266941, 0.748751, 0.440573, &
                 0.274149, 0.751988, 0.436601, &
                 0.281477, 0.755203, 0.432552, &
                 0.288921, 0.758394, 0.428426, &
                 0.296479, 0.761561, 0.424223, &
                 0.304148, 0.764704, 0.419943, &
                 0.311925, 0.767822, 0.415586, &
                 0.319809, 0.770914, 0.411152, &
                 0.327796, 0.773980, 0.406640, &
                 0.335885, 0.777018, 0.402049, &
                 0.344074, 0.780029, 0.397381, &
                 0.352360, 0.783011, 0.392636, &
                 0.360741, 0.785964, 0.387814, &
                 0.369214, 0.788888, 0.382914, &
                 0.377779, 0.791781, 0.377939, &
                 0.386433, 0.794644, 0.372886, &
                 0.395174, 0.797475, 0.367757, &
                 0.404001, 0.800275, 0.362552, &
                 0.412913, 0.803041, 0.357269, &
                 0.421908, 0.805774, 0.351910, &
                 0.430983, 0.808473, 0.346476, &
                 0.440137, 0.811138, 0.340967, &
                 0.449368, 0.813768, 0.335384, &
                 0.458674, 0.816363, 0.329727, &
                 0.468053, 0.818921, 0.323998, &
                 0.477504, 0.821444, 0.318195, &
                 0.487026, 0.823929, 0.312321, &
                 0.496615, 0.826376, 0.306377, &
                 0.506271, 0.828786, 0.300362, &
                 0.515992, 0.831158, 0.294279, &
                 0.525776, 0.833491, 0.288127, &
                 0.535621, 0.835785, 0.281908, &
                 0.545524, 0.838039, 0.275626, &
                 0.555484, 0.840254, 0.269281, &
                 0.565498, 0.842430, 0.262877, &
                 0.575563, 0.844566, 0.256415, &
                 0.585678, 0.846661, 0.249897, &
                 0.595839, 0.848717, 0.243329, &
                 0.606045, 0.850733, 0.236712, &
                 0.616293, 0.852709, 0.230052, &
                 0.626579, 0.854645, 0.223353, &
                 0.636902, 0.856542, 0.216620, &
                 0.647257, 0.858400, 0.209861, &
                 0.657642, 0.860219, 0.203082, &
                 0.668054, 0.861999, 0.196293, &
                 0.678489, 0.863742, 0.189503, &
                 0.688944, 0.865448, 0.182725, &
                 0.699415, 0.867117, 0.175971, &
                 0.709898, 0.868751, 0.169257, &
                 0.720391, 0.870350, 0.162603, &
                 0.730889, 0.871916, 0.156029, &
                 0.741388, 0.873449, 0.149561, &
                 0.751884, 0.874951, 0.143228, &
                 0.762373, 0.876424, 0.137064, &
                 0.772852, 0.877868, 0.131109, &
                 0.783315, 0.879285, 0.125405, &
                 0.793760, 0.880678, 0.120005, &
                 0.804182, 0.882046, 0.114965, &
                 0.814576, 0.883393, 0.110347, &
                 0.824940, 0.884720, 0.106217, &
                 0.835270, 0.886029, 0.102646, &
                 0.845561, 0.887322, 0.099702, &
                 0.855810, 0.888601, 0.097452, &
                 0.866013, 0.889868, 0.095953, &
                 0.876168, 0.891125, 0.095250, &
                 0.886271, 0.892374, 0.095374, &
                 0.896320, 0.893616, 0.096335, &
                 0.906311, 0.894855, 0.098125, &
                 0.916242, 0.896091, 0.100717, &
                 0.926106, 0.897330, 0.104071, &
                 0.935904, 0.898570, 0.108131, &
                 0.945636, 0.899815, 0.112838, &
                 0.955300, 0.901065, 0.118128, &
                 0.964894, 0.902323, 0.123941, &
                 0.974417, 0.903590, 0.130215, &
                 0.983868, 0.904867, 0.136897, &
                 0.993248, 0.906157, 0.143936]
      res = reshape(viridis, shape(res))
    end function viridis_data


    function turbocm(xx) result (rgb)
      real, intent(in) :: xx ! between <0, 1>
      real :: rgb(3)         ! pixel vector with values <0, 1>
!
!  Adopted to Fortran using source code at:
!    https://gist.github.com/mikhailov-work/0d177465a8151eb6ede1768d51d476c7
!    Polynomial approximation in GLSL for the Turbo colormap
!    Original LUT: https://gist.github.com/mikhailov-work/ee72ba4191942acecc03fe6da94fc73f
!    Authors:
!      Colormap Design: Anton Mikhailov (mikhailov@google.com)
!      GLSL Approximation: Ruofei Du (ruofei@google.com)
!
!  A blog about Turbo colormap
!  https://research.google/blog/turbo-an-improved-rainbow-colormap-for-visualization/
!
      real(DP), parameter :: kRedVec4(4) = [0.13572138_DP, 4.61539260_DP, -42.66032258_DP, 132.13108234_DP]
      real(DP), parameter :: kGreenVec4(4) = [0.09140261_DP, 2.19418839_DP, 4.84296658_DP, -14.18503333_DP]
      real(DP), parameter :: kBlueVec4(4) = [0.10667330_DP, 12.64194608_DP, -60.58204836_DP, 110.36276771_DP]
      real(DP), parameter :: kRedVec2(2) = [-152.94239396_DP, 59.28637943_DP]
      real(DP), parameter :: kGreenVec2(2) = [4.27729857_DP, 2.82956604_DP]
      real(DP), parameter :: kBlueVec2(2) = [-89.90310912_DP, 27.34824973_DP]
      real(DP) :: x, v4(4), v2(2)
  
      ! x = saturate(x); ! I assume it means to clamp "x" between <0, 1> ?
      x = min(1.0_DP, max(0.0_DP, real(xx,DP)))
      v4 = [1.0_DP, x, x*x, x*x*x]
      v2 = [v4(3), v4(4)] * v4(3)
      rgb(1) = real(dot_product(v4, kRedVec4)   + dot_product(v2, kRedVec2))
      rgb(2) = real(dot_product(v4, kGreenVec4) + dot_product(v2, kGreenVec2))
      rgb(3) = real(dot_product(v4, kBlueVec4)  + dot_product(v2, kBlueVec2))
      ! clamp output to <0,1> (was not present in the original code)
      rgb = max(0.0, min(1.0, rgb))
    end function turbocm

! -----------------------------------------------------------------------------
  end module pnmio_module
! -----------------------------------------------------------------------------

