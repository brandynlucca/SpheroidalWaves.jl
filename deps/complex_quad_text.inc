! Shared by the complex prolate and oblate quad wrapper modules.
! Inputs are decimal text. Outputs retain solver mantissas and exponents;
! no value or coordinate passes through a double-precision temporary.
! mode: 0=eigenvalue, 1=angular, 2=radial, 3=prolate radial offset x-1.
! Each point has eight real channels: Re/Im of value1, derivative1,
! value2, derivative2. Angular/eigenvalue calls use only the leading channels.
  subroutine complex_batch_quad_text(m, n, mode, option, npts, ctext, xtext, width, output, exponents, accuracy, status)
    integer(c_int), intent(in) :: m, n, mode, option, npts, width
    character(c_char), intent(in) :: ctext(*), xtext(*)
    character(c_char), intent(out) :: output(*)
    integer(c_int), intent(out) :: exponents(*), accuracy(*), status
    integer(c_int) :: i, j, li, lnum, narg, ioprad, iopang, iopnorm, repeat, threshold
    real(wk) :: cr, ci, radial_x
    complex(wk) :: c
    logical :: ok, radial_mode
    real(wk), allocatable :: points(:), arg(:)
    complex(wk), allocatable :: r1(:), dr1(:), r2(:), dr2(:), s(:,:), ds(:,:), eig(:)
    integer(c_int), allocatable :: er1(:), edr1(:), er2(:), edr2(:), ar(:), es(:,:), eds(:,:), ac(:,:), adc(:,:)

    status = 0
    radial_mode = mode == 2 .or. mode == 3
    if (m < 0 .or. n < m .or. npts < 1 .or. width < 70 .or. mode < 0 .or. mode > 3 .or. &
        (mode == 3 .and. quad_is_oblate)) then
      status = -4
      return
    end if
    if (radial_mode .and. (option < 1 .or. option > 4)) then
      status = -4
      return
    end if
    call quad_decode(ctext, 1_c_int, width, cr, ok)
    if (.not. ok) then
      status = -5
      return
    end if
    call quad_decode(ctext, width+1, width, ci, ok)
    if (.not. ok) then
      status = -5
      return
    end if
    c = cmplx(cr, ci, kind=wk)
    ! Portable finite check (avoids the `ieee_arithmetic` intrinsic module,
    ! missing from some gfortran cross-compiler shards): NaN never equals
    ! itself under IEEE comparison; finite values are bounded by huge().
    if (.not. is_finite(cr) .or. .not. is_finite(ci) .or. abs(c) == 0.0_wk) then
      status = -5
      return
    end if
    allocate(points(npts))
    do i = 1, npts
      call quad_decode(xtext, (i-1)*width+1, width, points(i), ok)
      if (.not. ok .or. .not. is_finite(points(i))) then
        status = -3
        return
      end if
      if (mode == 1 .and. abs(points(i)) > 1.0_wk) then
        status = -3
        return
      end if
      if (radial_mode) then
        if ((quad_is_oblate .and. points(i) < 0.0_wk) .or. &
            (mode == 2 .and. .not. quad_is_oblate .and. points(i) <= 1.0_wk) .or. &
            (mode == 3 .and. points(i) <= 0.0_wk)) then
          status = -3
          return
        end if
      end if
    end do
    li = n-m+1
    lnum = li
    if (quad_is_oblate) then
      threshold = int(2.0_wk*abs(c)/acos(-1.0_wk), c_int)
      if (lnum < threshold .and. mod(lnum, 2) /= 0) lnum = lnum+1
    end if
    narg = 1
    if (mode == 1) narg = npts
    allocate(arg(narg), r1(lnum), dr1(lnum), r2(lnum), dr2(lnum), eig(lnum))
    allocate(er1(lnum), edr1(lnum), er2(lnum), edr2(lnum), ar(lnum))
    allocate(s(lnum,narg), ds(lnum,narg), es(lnum,narg), eds(lnum,narg), ac(lnum,narg), adc(lnum,narg))
    ar = -1  ! First-kind-only calls may not provide an accuracy estimate.
    arg = 0.0_wk
    if (mode == 1) arg = points
    iopang = 0
    ioprad = 0
    iopnorm = 0
    if (mode == 1) then
      iopang = 2
      if (option /= 0) iopnorm = 1
    end if
    if (radial_mode) then
      ioprad = 2
      if (option == 1) ioprad = 1
    end if
    do i = 1, 8*npts
      call quad_encode(output, (i-1)*width+1, width, 0.0_wk)
      exponents(i) = 0
    end do
    accuracy(1:npts) = 0
    repeat = 1
    if (radial_mode) repeat = npts
    do j = 1, repeat
      radial_x = 1.0_wk
      if (quad_is_oblate) radial_x = 10.0_wk
      if (radial_mode) then
        radial_x = points(j)
        if (mode == 2 .and. .not. quad_is_oblate) radial_x = radial_x-1.0_wk
      end if
      call quad_solver(c, m, lnum, ioprad, radial_x, iopang, iopnorm, narg, arg, &
                       r1, er1, dr1, edr1, r2, er2, dr2, edr2, ar, &
                       s, es, ds, eds, ac, adc, eig)
      if (mode == 0) then
        call store_pair(1, eig(li), 0_c_int)
      else if (mode == 1) then
        do i = 1, npts
          call store_pair(8*(i-1)+1, s(li,i), es(li,i))
          call store_pair(8*(i-1)+3, ds(li,i), eds(li,i))
          accuracy(i) = ac(li,i)
        end do
      else
        call store_pair(8*(j-1)+1, r1(li), er1(li))
        call store_pair(8*(j-1)+3, dr1(li), edr1(li))
        if (option /= 1) then
          call store_pair(8*(j-1)+5, r2(li), er2(li))
          call store_pair(8*(j-1)+7, dr2(li), edr2(li))
        end if
        accuracy(j) = ar(li)
      end if
    end do
  contains
    subroutine store_pair(idx, value, exponent)
      integer(c_int), intent(in) :: idx, exponent
      complex(wk), intent(in) :: value
      call quad_encode(output, (idx-1)*width+1, width, real(value,wk))
      call quad_encode(output, idx*width+1, width, aimag(value))
      exponents(idx:idx+1) = exponent
    end subroutine store_pair

    logical function is_finite(x)
      real(wk), intent(in) :: x
      is_finite = (x == x) .and. (abs(x) <= huge(x))
    end function is_finite
  end subroutine complex_batch_quad_text

  subroutine quad_decode(buf, offset, width, value, ok)
    character(c_char), intent(in) :: buf(*)
    integer(c_int), intent(in) :: offset, width
    real(wk), intent(out) :: value
    logical, intent(out) :: ok
    character(len=width) :: text
    integer :: i, ios
    do i = 1, width
      text(i:i) = buf(offset+i-1)
    end do
    value = 0.0_wk
    read(text, *, iostat=ios) value
    ok = ios == 0
  end subroutine quad_decode

  subroutine quad_encode(buf, offset, width, value)
    character(c_char), intent(out) :: buf(*)
    integer(c_int), intent(in) :: offset, width
    real(wk), intent(in) :: value
    character(len=width) :: text
    integer :: i
    write(text, '(ES70.60E4)') value
    do i = 1, width
      buf(offset+i-1) = text(i:i)
    end do
  end subroutine quad_encode
