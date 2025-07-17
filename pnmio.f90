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
!   subroutine colormap(red,green,blue)
!  *subroutine assign_colormap_int(aa,rr,gg,bb,n,imin,imax,iwh,ibl)
!  *subroutine assign_colormap_2R(uu,rr,gg,bb,n,rmin,rmax)
!  *helper procedures for stream read
!  *conversion from signed to unsigned numbers
!
! (*) private objects

! SHORT DESCRIPTION
! These routines work 'raw' or 'plain' PNM formats (P4 - not yet implemented)


! WEB: https://netpbm.sourceforge.net/doc/pnm.html

! -----------------------------------------------------------------------------
  module pnmio_module
! -----------------------------------------------------------------------------
  use iso_fortran_env, only : error_unit, int8, int16, iostat_end
  implicit none
  private
  public readpnm, writepgm, writeppm, colormap, assign_colormap

  ! default value of optional MAXVAL argument
  integer, parameter :: MX_DEFAULT = 255

  integer, parameter :: IOMSG_MAXLEN = 100

  interface assign_colormap
    module procedure assign_colormap_int
    module procedure assign_colormap_2R
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

      if ((p-1)/3==1) then
        ! P4, P5 and P6 can contain more images in single file
        inquire(fid, size=fsize, pos=fpos)
        if (fpos-1 < fsize) then
          print '("File size ",i0,". File position ",i0)', fsize, fpos
          write(error_unit,'("Warning readpnm: file contains more images")')
          ! TODO ignoring for now
        end if
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
     !open(newunit=fid, file=filename, status='replace', access='stream', &
     !  form='formatted', iostat=ios, iomsg=iomsg)
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
        write(buf,'(a)') '# xreated by pnmio fortran library'
        write(fid) trim(buf)//new_line(buf)
        write(buf,'(i0,1x,i0)') w, h
        write(fid) trim(buf)//new_line(buf)

        if (mod(p,3)/=1) then
          write(buf,'(i0)') mx
          write(fid) trim(buf)//achar(10)
        end if
      end block HEADER

      ! write raster
      if (p > 3) then
        ! binary format
        BINARY: block
          integer(int8), allocatable :: raster_1b(:)
          integer(int16), allocatable :: raster_2b(:)

          ! write raster (1B or 2B)
          if (mx > 255) then
            if (p==6) then
              allocate(raster_2b(3*w*h)) ! PPM (2 bytes)
            else
              allocate(raster_2b(w*h))   ! PGM (2 bytes)
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
        ! ascii format
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
        end block ASCII
        if (ios /= 0) then
          write(error_unit,'(a)') 'writing ASCII raster iomsg: '//trim(iomsg)
          exit BLK
        end if
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
! Write PGM or PBM file. Just a wrapper to "writepnm_3d"
! PBM - if "mx" is present and "mx==1"
! PGM - otherwise
!
    integer :: mx0
    integer, allocatable :: atmp(:,:,:)

    if (present(mx)) then
      mx0 = mx
    else
      mx0 = MX_DEFAULT
    end if
    allocate(atmp(1, size(aa,1), size(aa,2)))
    atmp(1,:,:) = aa
    call writepnm_3d(filename, atmp, mx0, is_plain, ierr)

  end subroutine writepgm


  subroutine writeppm_3x2d(filename, rr, gg, bb, mx, is_plain, ierr)
    character(len=*), intent(in)   :: filename
    integer, intent(in)            :: rr(:,:), gg(:,:), bb(:,:)
    integer, intent(in), optional  :: mx
    logical, intent(in), optional  :: is_plain
    integer, intent(out), optional :: ierr
!
! Write PPM file. Just a wrapper to "writepnm_3d"
!
    integer :: mx0
    integer, allocatable :: atmp(:,:,:)

    if (present(mx)) then
      mx0 = mx
    else
      mx0 = MX_DEFAULT
    end if

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
    call writepnm_3d(filename, atmp, mx0, is_plain, ierr)

  end subroutine writeppm_3x2d


  subroutine writeppm_1x3d(filename, aa, mx, is_plain, ierr)
    character(len=*), intent(in)   :: filename
    integer, intent(in)            :: aa(:,:,:)
    integer, intent(in), optional  :: mx
    logical, intent(in), optional  :: is_plain
    integer, intent(out), optional :: ierr
!
! Write PPM file. Just a wrapper to "writepnm_3d"
!
    integer :: mx0

    if (present(mx)) then
      mx0 = mx
    else
      mx0 = MX_DEFAULT
    end if

    if (size(aa,1)/=3) then
      if (present(ierr)) then
        ierr = -3
        return
      else
        error stop 'writeppm: first extent of aa must be 3'
      end if
    end if

    call writepnm_3d(filename, aa, mx0, is_plain, ierr)

  end subroutine writeppm_1x3d


! -----------------------------------------------------------------------------
  subroutine colormap(red,green,blue)
! -----------------------------------------------------------------------------
! Generate colormap (blue  cyan  green  yellow  red)
! First and last colors are white and black, respectivelly

  integer, intent(out) :: red(:), green(:), blue(:)

  ! - local vars
  integer, parameter :: mxval = 254
  integer            :: n, i

  real :: x, yr, yb, yg
  ! -
  n = size(red,dim=1)

  do i=2,n-1
    x = float(i-2) / float(n-3) ! x will be in range (0;1)

    yr =  4.*x - 2.
    yb = -4.*x + 2.
    if(x<0.5) then
      yg =  4.*x + 0.
    else
      yg = -4.*x + 4.
    endif

    yr = max(0.,yr)
    yg = max(0.,yg)
    yb = max(0.,yb)
    yr = min(1.,yr)
    yg = min(1.,yg)
    yb = min(1.,yb)

    red(i)   = int(yr*float(mxval+1)) ! values will be between 0..mxval
    green(i) = int(yg*float(mxval+1))
    blue(i)  = int(yb*float(mxval+1))

  enddo

  red(1) = mxval; red(n) = 0
  green(1) = mxval; green(n) = 0
  blue(1) = mxval; blue(n) = 0
! -----------------------------------------------------------------------------
  end subroutine colormap
! -----------------------------------------------------------------------------



! -----------------------------------------------------------------------------
  subroutine assign_colormap_int(aa,rr,gg,bb,n,imin,imax,iwh,ibl)
! -----------------------------------------------------------------------------
  integer, intent(in)  :: aa(:,:)
  integer, intent(out) :: rr(:,:), gg(:,:), bb(:,:)
  integer, intent(in)  :: n
  integer, intent(in), optional :: imin, imax, iwh, ibl

  ! - local vars

  integer :: cmap(n,3)
  real    :: dx, fmin, fmax, f
  integer :: i, j, m, nx, ny

  ! -
  call colormap(cmap(1:n,1),cmap(1:n,2),cmap(1:n,3))

  fmin = float(minval(aa))
  fmax = float(maxval(aa))
  if (present(imin)) fmin = float(imin)
  if (present(imax)) fmax = float(imax)
  nx = size(aa,dim=1)
  ny = size(aa,dim=2)
print *, 'assign colormap :',fmin,fmax

  ! to elements aa==fmin assign 2nd color,
  ! to elements aa==fmax assign (n-1)th color
  ! 1st and nth colors are reserved for elements out of bond, or iwh,ibl

  dx = (fmax-fmin)/float(n-3)

  do i=1,nx
  do j=1,ny
    f = (float(aa(i,j)) - fmin) / dx
    m = int(f) + 2
    m = min(max(1,m),n)
    rr(i,j) = cmap(m,1);  gg(i,j) = cmap(m,2);  bb(i,j) = cmap(m,3)
  enddo
  enddo

  if (present(iwh)) then
    where(aa==iwh)
      rr = cmap(1,1)
      gg = cmap(1,2)
      bb = cmap(1,3)
    endwhere
  endif
  if (present(ibl)) then
    where(aa==ibl)
      rr = cmap(n,1)
      gg = cmap(n,2)
      bb = cmap(n,3)
    endwhere
  endif
! -----------------------------------------------------------------------------
  end subroutine assign_colormap_int
! -----------------------------------------------------------------------------



! -----------------------------------------------------------------------------
  subroutine assign_colormap_2R(uu,rr,gg,bb,n,rmin,rmax)
! -----------------------------------------------------------------------------
  real(8), intent(in)  :: uu(:,:)
  integer, intent(out) :: rr(:,:), gg(:,:), bb(:,:)
  integer, intent(in)  :: n
  real(8), intent(in), optional :: rmin, rmax

  integer :: aa(size(uu,dim=1), size(uu,dim=2))
  real(8) :: rmin0, rmax0

  ! -

  if (present(rmin)) then
    rmin0 = rmin
  else
    rmin0 = minval(uu)
  endif

  if (present(rmax)) then
    rmax0 = rmax
  else
    rmax0 = maxval(uu)
  endif

  aa = nint((uu-rmin0)/(rmax0-rmin0)*float(n))
  call assign_colormap_int(aa,rr,gg,bb,n,imin=0,imax=n)
! -----------------------------------------------------------------------------
  end subroutine assign_colormap_2R
! -----------------------------------------------------------------------------



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

! -----------------------------------------------------------------------------
  end module pnmio_module
! -----------------------------------------------------------------------------

