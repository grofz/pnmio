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
  use iso_fortran_env, only : error_unit, int8, int16
  implicit none
  private
  public readpnm, readpnm_legacy, writepgm, writeppm, colormap, assign_colormap
  public consume_header, consume_decimal, consume_magick ! public only for testing

  ! SWITCH 'iplain' TO SELECT RAW OR PLAIN FORMAT DURING WRITING:
  ! 1. use false on LINUX, will produce raw-format (smaller files)
  logical, parameter :: iplain = .false.
  ! 2. use true on WINDOWS, or if netpbm library is not installed (larger files)
  !logical, parameter :: iplain = .true.

  character(len=*), parameter :: tmpfile = '.tmp.plain'

  integer, parameter :: IOMSG_MAXLEN = 100

  interface assign_colormap
    module procedure assign_colormap_int
    module procedure assign_colormap_2R
  end interface assign_colormap

  interface readpnm
    module procedure readpnm_2d
    module procedure readpnm_3d
  end interface

  interface s2u
    module procedure s2u_8
    module procedure s2u_16
  end interface

  contains

! ==========================
! READ PBM / PGM / PPM FILES
! ==========================

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
      ierr = ierr0
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
!   - RGB components making a tuple (first index is component id)
!   - rows containing W tuples (second index is column id)
!   - raster containing H rows (third index is row id)
!
    integer :: p, w, h, mx, ierr0, ios, fid
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
        if (mx < 256) then
          allocate(raster_1b(w*h*3))
        else
          allocate(raster_2b(w*h*3))
        end if
      else
        ! PGM or PBM
        allocate(aa(1,w,h))
        if (mx < 256) then
          allocate(raster_1b(w*h))
        else
          allocate(raster_2b(w*h))
        end if
      end if

      if (mx < 256) then
        call consume_raster(fid, p, ierr0, r1b=raster_1b)
        if (ierr0 /= 0) exit BLK
        aa = reshape(s2u(raster_1b), shape(aa))
      else
        call consume_raster(fid, p, ierr0, r2b=raster_2b)
        if (ierr0 /= 0) exit BLK
        aa = reshape(s2u(raster_2b), shape(aa))
      end if

      if (present(ierr)) ierr = 0
      return
    end block BLK

    ! error occurred
    if (present(ierr)) then
      ierr = ierr0
    else
      error stop 'readpnm error - see above'
    end if

  end subroutine readpnm_3d


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
  end subroutine


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
  end subroutine


  subroutine consume_header(fid, p, w, h, maxval, ierr)
    integer, intent(in)  :: fid
    integer, intent(out) :: p, w, h, maxval, ierr

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
        call consume_decimal(fid, maxval, ierr)
        if (ierr /= 0) exit HDR
      else
        maxval = 1
      end if
      ierr = 0
      return
    end block HDR

    ! an error occurred
    ierr = -1
    write(error_unit,'(a)') 'error: header is invalid'
  end subroutine


  subroutine consume_raster(fid, p, ierr, r1b, r2b)
    integer, intent(in) :: fid, p
    integer, intent(out) :: ierr
    integer(int8), intent(out), optional :: r1b(:)
    integer(int16), intent(out), optional :: r2b(:)

    ! assuming "fid" is opened as "unformatted stream"
    integer :: ios, i
    character(len=IOMSG_MAXLEN) :: iomsg
    integer, allocatable :: u(:)

    if (present(r1b) .eqv. present(r2b)) &
      error stop 'consume_raster: only r1b or r2b must be given'
    ierr = -1

    if ((p-1)/3 == 1) then
      ! binary format
      if (present(r1b)) then
        read(fid, iostat=ios, iomsg=iomsg) r1b
      else
        read(fid, iostat=ios, iomsg=iomsg) r2b
      end if
      if (ios/=0) then
        write(error_unit, '(a)') 'raster iomsg: '//trim(iomsg)
        return
      end if

    else if ((p-1)/3 == 0) then
      ! raw format (ASCII)
      if (present(r1b)) then
        allocate(u(size(r1b)))
      else
        allocate(u(size(r2b)))
      end if

      do i=1, size(u)
        call consume_decimal(fid, u(i), ierr)
        if (ierr /= 0) then
          write(error_unit, '("raster raw error at posiiton ",i0," out of ",i0)') i, size(u)
          return
        end if
      end do

      ! convert to unsigned so it fits to 1B or 2B integers
      if (present(r1b)) then
        r1b = u2s_8(u)
      else
        r2b = u2s_16(u)
      end if
    end if

    ierr = 0
  end subroutine consume_raster


  !
  ! Conversion from signed to unsigned numbers
  !
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

