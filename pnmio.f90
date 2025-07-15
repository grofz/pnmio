!CONTENT
! module pnmio_module
!   subroutine readpnm(filename,nx,ny,nmax,aa)
!   subroutine writepgm(filename,aa,nmax)
!   subroutine writeppm(filename,rr,gg,bb)
!   subroutine colormap(red,green,blue)
!   subroutine assign_colormap2(aa,rr,gg,bb,n,imin,imax,iwh,ibl)

! SHORT DESCRIPTION
! These routines work with 'raw' PNM formats

! NOTE: these subroutines use system command 'pnmtopnm' for conversions
! between 'raw' and 'plain' formats. If this or similar tool is not
! available then these routines will work with 'plain' formats only.
!
! For writing plain formats -> modify parameter "iplain" at the start of module
!
! WEB: https://netpbm.sourceforge.net/doc/pnm.html

! -----------------------------------------------------------------------------
  module pnmio_module
! -----------------------------------------------------------------------------
  implicit none
  private
  public readpnm_legacy, writepgm, writeppm, colormap, assign_colormap
  public consume_decimal

  ! SWITCH 'iplain' TO SELECT RAW OR PLAIN FORMAT DURING WRITING:
  ! 1. use false on LINUX, will produce raw-format (smaller files)
  logical, parameter :: iplain = .false.
  ! 2. use true on WINDOWS, or if netpbm library is not installed (larger files)
  !logical, parameter :: iplain = .true.

  character(len=*), parameter :: tmpfile = '.tmp.plain'

  interface assign_colormap
    module procedure assign_colormap_int
    module procedure assign_colormap_2R
  end interface assign_colormap

  contains

! -----------------------------------------------------------------------------
  recursive subroutine readpnm_legacy(filename, aa)
! -----------------------------------------------------------------------------
!NOTE: in pnm format specification lines begining with '#' are comments
!      this has not been implemented yet (TODO)
  character(len=*)     :: filename
  integer, allocatable :: aa(:,:)

  ! - local vars
  integer            :: nx, ny, nmax
  integer            :: fid, fid2
  logical            :: ex
  integer            :: magick ! 1,2,3 for PBM,PGM,PPM respectivelly
  character(len=200) :: tmp
  character(len=1)   :: onechar
  integer            :: i, j, ichar, exitstat, cmdstat

  ! to convert 'raw format to plain we use system call:
  ! $ pnmtopnm -plain <filename> > tmpfilename

  ! - open file & read header to get magick number
  inquire(file=filename, exist=ex)
  if (.not. ex) then
    print *, 'readpnm error: file do not exist '
    print *, trim(filename)
    stop 1
  endif

  open(newunit=fid, file = filename, status = 'old' )

  read(fid,*) tmp ! get the magick number
  select case(trim(tmp))
  case('P1') ! plain PBM format
    magick = 1
    read(fid,*) nx, ny
    nmax = 1
  case('P2') ! plain PGM format
    magick = 2
    read(fid,*) nx, ny
    read(fid,*) nmax
  case('P3')
    magick = 3
    stop 'readpnm error: P3 (ppm) format not ready, go and code it now!'
  case('P4') ! raw PBM format
    magick = 4
    read(fid,*) nx, ny
    nmax = 1
  case('P5') ! raw PGM format
    magick = 5
    read(fid,*) nx, ny
    read(fid,*) nmax
  case('P6')
    stop 'readpnm error: P6 (raw ppm) format not ready, go and code it now!'
  case default
    stop 'readpnm error: magick number not recognized, wrong format?'
  end select

  ! - use system call to convert "raw" to "plain" files
  if (magick > 3) then
    print *, 'readpnm info: converting raw format to plain...'
    call execute_command_line( 'pnmtopnm -plain '//trim(filename)//' > '//trim(tmpfile), &
      & cmdstat=cmdstat, exitstat=exitstat)
    if (cmdstat /= 0 .or. exitstat /= 0) then
      print *, 'readpnm: converstion failed :', cmdstat, exitstat
      stop 1
    endif

    call readpnm_legacy(filename=tmpfile, aa=aa)

    print *, 'readpnm info: deleting temporary file...'
    open(newunit=fid2, status='old', file=tmpfile, iostat=cmdstat)
    if (cmdstat/=0) then
      print *, 'deleting temporary file failed, file could not be opened'
      stop 1
    end if
    close(fid2, status='delete', iostat=cmdstat)
    if (cmdstat/=0) then
      print *, 'deleting temporary file failed, file could not be deleted'
      stop 1
    end if
    close(fid)
    return
  endif

  ! - reading "aa"
  if (allocated(aa)) then
    if (any(shape(aa) /= [nx, ny])) deallocate(aa)
  end if
  if (.not. allocated(aa)) allocate(aa(nx,ny))

  select case(magick)
  case(1) ! PBM ! pbm is tricky as there are no spaces between ones/zeroes
    i=0; j=1
    do
      read(fid,'(a)',end=100) tmp

      do ichar=1,len(trim(tmp))
        onechar = tmp(ichar:ichar)
        if (onechar == ' ') cycle
        if (onechar /= '1' .and. onechar /= '0') then
          print *,'readpnm error: invalid character :',onechar
          stop
        endif

        ! now read next byte (update position)
        if (i==nx) then
          i=1
          if (j==ny) then
            print *, 'read_pnm warning: file contains too much data!'
            stop
          endif
          j=j+1
        else
          i=i+1
        endif

        if (onechar=='0') then
          aa(i,j) = 0
        else
          aa(i,j) = 1
        endif

      enddo

      cycle   ! next line from file
      100 exit
    enddo

    if (i /= nx .or. j /= ny) then
      print *, 'readpnm error: file ended to soon! '
      print *, 'bytes read ',i,j
      print *, 'bytes expected ',nx,ny
      stop
    endif

print *, 'read_pnm info: PBM file read into aa', float(sum(aa))/count(aa>=0)

  case(2) ! PGM (much easier)
    read(fid,*) aa(1:nx,1:ny)
print *, 'read_pnm info: PGM file read into aa', float(sum(aa))/count(aa>=0)

  case(3) ! PPM (not ready yet)
  case default
print *, magick
    stop 'readpnm unexpected error'
  end select

  close(fid)
! -----------------------------------------------------------------------------
  end subroutine readpnm_legacy
! -----------------------------------------------------------------------------



! -----------------------------------------------------------------------------
  subroutine writepgm(filename, aa, nmax)
! -----------------------------------------------------------------------------
  character(len=*)  :: filename
  integer           :: aa(:,:)
  integer, optional :: nmax

  integer            :: fid
  integer, parameter :: imax_default = 255
  integer :: ix, iy, imax

  ! ---

  ix = size(aa,dim=1)
  iy = size(aa,dim=2)
  if (present(nmax)) then
    imax = nmax
    if(imax>999) print *, 'writepgm warning: go to src and change format!'
  else
    imax = imax_default
  endif

  if (iplain) then
    open(newunit=fid, file = filename, status = 'replace')
  else
    open(newunit=fid, file = tmpfile, status = 'replace')
  endif

  write(fid, '(a)') 'P2'
  write(fid, '(i5,1x,i5)') ix, iy
  write(fid, '(i5)') imax
  write(fid, '(10(i3,1x))') min(imax,max(0,aa))
  close(fid)

  ! make system call to convert to raw and delete tmpfile
  if (.not. iplain) call convert_plain2raw(filename)
! -----------------------------------------------------------------------------
  end subroutine writepgm
! -----------------------------------------------------------------------------



! -----------------------------------------------------------------------------
  subroutine writeppm(filename,rr,gg,bb)
! -----------------------------------------------------------------------------
  character(len=*)   :: filename
  integer            :: rr(:,:), gg(:,:), bb(:,:)
  integer, parameter :: nmax=255

  integer            :: fid
  integer :: ix, iy, i, j

  ! ---

  ix = size(rr,dim=1)
  iy = size(rr,dim=2)
  if (ix /= size(gg,dim=1) .or. iy /= size(gg,dim=2)) then
    print *, 'writeppm error: red and green arrays not same size'
    stop 1
  endif
  if (ix /= size(bb,dim=1) .or. iy /= size(bb,dim=2)) then
    print *, 'writeppm error: red and blue arrays not same size'
    stop 1
  endif

  if (iplain) then
    open(newunit=fid, file = filename, status = 'replace')
  else
    open(newunit=fid, file = tmpfile, status = 'replace')
  endif

  write(fid, '(a)') 'P3'
  write(fid, '(i5,1x,i5)') ix, iy
  write(fid, '(i5)') nmax
  write(fid, '(10(i3,1x))') ((rr(i,j),gg(i,j),bb(i,j), i=1,ix), j=1,iy)
  close(fid)

  ! do system call to convert to raw and delete tempfile
  if (.not. iplain) call convert_plain2raw(filename)
! -----------------------------------------------------------------------------
  end subroutine writeppm
! -----------------------------------------------------------------------------



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



! -----------------------------------------------------------------------------
  subroutine convert_plain2raw(raw)
    character(len=*), intent(in) :: raw
    integer :: istat, cmdstat, fid

    ! call pnmtopnm to convert from plain to raw
    call execute_command_line( &
      & 'pnmtopnm '//trim(tmpfile)//' > '//trim(raw), exitstat=istat, cmdstat=cmdstat)

    if (istat /= 0 .or. cmdstat /= 0) then
      ! if conversion failed leave file in plain format and just rename it
      print '("PNMIO info: Conversion failed. Fail-back to plain format")'
      call execute_command_line( &
        & 'mv '//trim(tmpfile)//' '//trim(raw), exitstat=istat)
      if (istat /= 0) then
        print '("PNMIO Error: mv command did not work")'
        stop 1
      end if
    else
      ! delete a temporaty file if conversion succeeded
      open(newunit=fid, file=tmpfile, status='old', iostat=istat)
      if (istat /= 0) error stop 'opening tmpfile failed'
      close(fid, status='delete')
    end if
  end subroutine convert_plain2raw
! -----------------------------------------------------------------------------


! ---------------------------------
! helper procedures for stream read
! ---------------------------------

  pure logical function is_whitespace(ch)
    character(len=1), intent(in) :: ch

    ! Assuming white space is a character with ASCII code between
    ! 9 and 13 (TAB, LF, VT, FF, CR) or a space
    is_whitespace = (iachar(ch)>=9 .and. iachar(ch)<=13) .or. ch==' '
  end function

  pure logical function is_newline(ch)
    character(len=1), intent(in) :: ch

    ! Assuming new line is LF
    is_newline = iachar(ch) == 10 ! iachar(new_line(ch))
  end function


  subroutine consume_decimal(fid, val)
    integer, intent(in) :: fid
    integer, intent(out) :: val

    ! consume ASCII decimal from the stream
    ! - read and ignore all whitespace characters
    ! - if '#' read, then ignore all characters until CR or LF

    integer, parameter :: MODE_SCAN=0, MODE_COMMENT=1, MODE_READ=2
    integer :: ios, pos, mode, dec_len 
    character(len=1) :: ch
    character(len=10) :: dec
    character(len=100) :: iomsg

    mode = MODE_SCAN
    dec = ''
    dec_len = 0
    inquire(unit=fid, pos=pos)
    do
      read(fid, iostat=ios, pos=pos, iomsg=iomsg) ch
      pos = pos + 1
      if (ios /= 0) then
        print *, 'iomsg: '//trim(iomsg)
        error stop 'consume_decimal: read error (see above)'
      end if

      select case(mode)
      case(MODE_SCAN)
        if (.not. is_whitespace(ch)) then
          if (ch=='#') then
            mode = MODE_COMMENT
          else
            mode = MODE_READ
            pos = pos - 1 ! re-read
          end if
        end if
      case(MODE_READ)
        if (is_whitespace(ch)) exit
        dec_len = dec_len + 1
        dec(dec_len:dec_len) = ch
      case(MODE_COMMENT)
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
      print *, 'iomsg: '//trim(iomsg)
      error stop 'consume_decimal: conversion error (see above)'
    end if
  end subroutine



! -----------------------------------------------------------------------------
  end module pnmio_module
! -----------------------------------------------------------------------------

