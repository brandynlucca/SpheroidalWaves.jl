module oblate_batch_fortran
  use, intrinsic :: iso_c_binding
  use oblate_parameters, only: knd
  use oblate_swf, only: oblfcn
  implicit none

  integer, parameter :: rk = c_double
  integer, parameter :: wk = knd

contains

  pure subroutine split_wk_to_double_pair(x, hi, lo)
    real(wk), intent(in) :: x
    real(rk), intent(out) :: hi, lo
    real(wk) :: hi_w

    hi = real(x, rk)
    hi_w = real(hi, wk)
    lo = real(x - hi_w, rk)
  end subroutine split_wk_to_double_pair

  pure real(wk) function pair_double_to_wk(hi, lo) result(v)
    real(rk), intent(in) :: hi, lo
    v = real(hi, wk) + real(lo, wk)
  end function pair_double_to_wk

  pure integer(c_int) function solver_lnum(real_c, requested_lnum) result(adjusted_lnum)
    real(wk), intent(in) :: real_c
    integer(c_int), intent(in) :: requested_lnum
    integer(c_int) :: even_threshold

    adjusted_lnum = requested_lnum
    even_threshold = int((2.0_wk * abs(real_c)) / acos(-1.0_wk), c_int)
    if (adjusted_lnum < even_threshold .and. mod(adjusted_lnum, 2_c_int) /= 0_c_int) then
      adjusted_lnum = adjusted_lnum + 1_c_int
    end if
  end function solver_lnum

  pure real(wk) function pow10_i(exp10) result(v)
    integer(c_int), intent(in) :: exp10
    v = 10.0_wk ** real(exp10, wk)
  end function pow10_i

  pure logical function is_finite_quad(x) result(ok)
    real(rk), intent(in) :: x
    ok = (x == x) .and. (abs(x) < huge(x))
  end function is_finite_quad

  subroutine decode_real_text(buf, offset, str_len, x, ok)
    character(c_char), intent(in) :: buf(*)
    integer(c_int), intent(in) :: offset, str_len
    real(wk), intent(out) :: x
    logical, intent(out) :: ok

    integer :: j, ios
    character(len=:), allocatable :: tmp

    allocate(character(len=str_len) :: tmp)
    do j = 1, str_len
      tmp(j:j) = achar(iachar(buf(offset + j - 1)))
    end do
    read(tmp, *, iostat=ios) x
    ok = (ios == 0)
  end subroutine decode_real_text

  subroutine encode_real_text(buf, offset, str_len, x)
    character(c_char), intent(out) :: buf(*)
    integer(c_int), intent(in) :: offset, str_len
    real(wk), intent(in) :: x

    integer :: j, ncopy
    character(len=:), allocatable :: tmp

    allocate(character(len=str_len) :: tmp)
    write(tmp, '(ES70.60E4)') x
    tmp = adjustl(tmp)

    do j = 1, str_len
      buf(offset + j - 1) = ' '
    end do

    ncopy = min(len_trim(tmp), str_len)
    do j = 1, ncopy
      buf(offset + j - 1) = tmp(j:j)
    end do
  end subroutine encode_real_text

  subroutine oblate_smn_batch_quad(m, n, c, n_eta, eta, normalize, value, derivative, status) bind(C, name="oblate_smn_batch_quad")
    integer(c_int), value, intent(in) :: m, n, n_eta, normalize
    real(c_double), value, intent(in) :: c
    real(c_double), intent(in) :: eta(*)
    real(c_double), intent(out) :: value(*), derivative(*)
    integer(c_int), intent(out) :: status

    integer(c_int) :: i, idx0, lnum, ioprad, iopang, iopnorm, narg
    real(wk) :: x
    real(wk) :: c_w

    real(wk), allocatable :: arg(:)
    real(wk), allocatable :: r1c(:), r1dc(:), r2c(:), r2dc(:)
    real(wk), allocatable :: s1c(:,:), s1dc(:,:)
    integer(c_int), allocatable :: ir1e(:), ir1de(:), ir2e(:), ir2de(:), naccr(:)
    integer(c_int), allocatable :: is1e(:,:), is1de(:,:), naccs(:,:)
    real(wk), allocatable :: eigout(:)

    status = 0
    if (n_eta < 1_c_int) then
      status = -4_c_int
      return
    end if
    if (n < 0_c_int) then
      status = -1_c_int
      return
    end if
    if (n < m) then
      status = -2_c_int
      return
    end if

    c_w = real(c, wk)

    do i = 1, n_eta
      if (abs(eta(i)) > 1.0_rk .or. .not. is_finite_quad(eta(i))) then
        status = -3_c_int
        return
      end if
    end do

    lnum = solver_lnum(c_w, n - m + 1_c_int)
    idx0 = n - m + 1_c_int
    narg = n_eta

    ioprad = 0_c_int
    iopang = 2_c_int
    iopnorm = 0_c_int
    if (normalize /= 0_c_int) iopnorm = 1_c_int
    x = 10.0_wk

    allocate(arg(narg))
    arg(:) = real(eta(1:narg), wk)
    allocate(r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum))
    allocate(ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), naccr(lnum))
    naccr = -1_c_int
    allocate(s1c(lnum, narg), s1dc(lnum, narg))
    allocate(is1e(lnum, narg), is1de(lnum, narg), naccs(lnum, narg))
    allocate(eigout(lnum))

    call oblfcn(c_w, m, lnum, ioprad, x, iopang, iopnorm, narg, arg, &
                r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                s1c, is1e, s1dc, is1de, naccs, eigout)

    do i = 1, n_eta
      value(i) = s1c(idx0, i) * pow10_i(is1e(idx0, i))
      derivative(i) = s1dc(idx0, i) * pow10_i(is1de(idx0, i))
    end do
  end subroutine oblate_smn_batch_quad

  subroutine oblate_smn_batch_quad_split(m, n, c, n_eta, eta, normalize, value_hi, value_lo, derivative_hi, derivative_lo, status) bind(C, name="oblate_smn_batch_quad_split")
    integer(c_int), value, intent(in) :: m, n, n_eta, normalize
    real(c_double), value, intent(in) :: c
    real(c_double), intent(in) :: eta(*)
    real(c_double), intent(out) :: value_hi(*), value_lo(*), derivative_hi(*), derivative_lo(*)
    integer(c_int), intent(out) :: status

    integer(c_int) :: i, idx0, lnum, ioprad, iopang, iopnorm, narg
    real(wk) :: x
    real(wk) :: c_w
    real(wk) :: value_wk, derivative_wk

    real(wk), allocatable :: arg(:)
    real(wk), allocatable :: r1c(:), r1dc(:), r2c(:), r2dc(:)
    real(wk), allocatable :: s1c(:,:), s1dc(:,:)
    integer(c_int), allocatable :: ir1e(:), ir1de(:), ir2e(:), ir2de(:), naccr(:)
    integer(c_int), allocatable :: is1e(:,:), is1de(:,:), naccs(:,:)
    real(wk), allocatable :: eigout(:)

    status = 0
    if (n_eta < 1_c_int) then
      status = -4_c_int
      return
    end if
    if (n < 0_c_int) then
      status = -1_c_int
      return
    end if
    if (n < m) then
      status = -2_c_int
      return
    end if

    c_w = real(c, wk)

    do i = 1, n_eta
      if (abs(eta(i)) > 1.0_rk .or. .not. is_finite_quad(eta(i))) then
        status = -3_c_int
        return
      end if
    end do

    lnum = solver_lnum(c_w, n - m + 1_c_int)
    idx0 = n - m + 1_c_int
    narg = n_eta

    ioprad = 0_c_int
    iopang = 2_c_int
    iopnorm = 0_c_int
    if (normalize /= 0_c_int) iopnorm = 1_c_int
    x = 10.0_wk

    allocate(arg(narg))
    arg(:) = real(eta(1:narg), wk)
    allocate(r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum))
    allocate(ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), naccr(lnum))
    naccr = -1_c_int
    allocate(s1c(lnum, narg), s1dc(lnum, narg))
    allocate(is1e(lnum, narg), is1de(lnum, narg), naccs(lnum, narg))
    allocate(eigout(lnum))

    call oblfcn(c_w, m, lnum, ioprad, x, iopang, iopnorm, narg, arg, &
                r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                s1c, is1e, s1dc, is1de, naccs, eigout)

    do i = 1, n_eta
      value_wk = s1c(idx0, i) * pow10_i(is1e(idx0, i))
      derivative_wk = s1dc(idx0, i) * pow10_i(is1de(idx0, i))
      call split_wk_to_double_pair(value_wk, value_hi(i), value_lo(i))
      call split_wk_to_double_pair(derivative_wk, derivative_hi(i), derivative_lo(i))
    end do
  end subroutine oblate_smn_batch_quad_split

  subroutine oblate_smn_batch_quad_fullsplit(m, n, c_hi, c_lo, n_eta, eta_hi, eta_lo, normalize, value_hi, value_lo, derivative_hi, derivative_lo, status) bind(C, name="oblate_smn_batch_quad_fullsplit")
    integer(c_int), value, intent(in) :: m, n, n_eta, normalize
    real(c_double), value, intent(in) :: c_hi, c_lo
    real(c_double), intent(in) :: eta_hi(*), eta_lo(*)
    real(c_double), intent(out) :: value_hi(*), value_lo(*), derivative_hi(*), derivative_lo(*)
    integer(c_int), intent(out) :: status

    integer(c_int) :: i, idx0, lnum, ioprad, iopang, iopnorm, narg
    real(wk) :: x
    real(wk) :: c_w
    real(wk) :: eta_w, value_wk, derivative_wk

    real(wk), allocatable :: arg(:)
    real(wk), allocatable :: r1c(:), r1dc(:), r2c(:), r2dc(:)
    real(wk), allocatable :: s1c(:,:), s1dc(:,:)
    integer(c_int), allocatable :: ir1e(:), ir1de(:), ir2e(:), ir2de(:), naccr(:)
    integer(c_int), allocatable :: is1e(:,:), is1de(:,:), naccs(:,:)
    real(wk), allocatable :: eigout(:)

    status = 0
    if (n_eta < 1_c_int) then
      status = -4_c_int
      return
    end if
    if (n < 0_c_int) then
      status = -1_c_int
      return
    end if
    if (n < m) then
      status = -2_c_int
      return
    end if

    c_w = pair_double_to_wk(c_hi, c_lo)

    do i = 1, n_eta
      eta_w = pair_double_to_wk(eta_hi(i), eta_lo(i))
      if (abs(eta_w) > 1.0_wk .or. .not. is_finite_quad(real(eta_w, rk))) then
        status = -3_c_int
        return
      end if
    end do

    lnum = solver_lnum(c_w, n - m + 1_c_int)
    idx0 = n - m + 1_c_int
    narg = n_eta

    ioprad = 0_c_int
    iopang = 2_c_int
    iopnorm = 0_c_int
    if (normalize /= 0_c_int) iopnorm = 1_c_int
    x = 10.0_wk

    allocate(arg(narg))
    do i = 1, narg
      arg(i) = pair_double_to_wk(eta_hi(i), eta_lo(i))
    end do
    allocate(r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum))
    allocate(ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), naccr(lnum))
    naccr = -1_c_int
    allocate(s1c(lnum, narg), s1dc(lnum, narg))
    allocate(is1e(lnum, narg), is1de(lnum, narg), naccs(lnum, narg))
    allocate(eigout(lnum))

    call oblfcn(c_w, m, lnum, ioprad, x, iopang, iopnorm, narg, arg, &
                r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                s1c, is1e, s1dc, is1de, naccs, eigout)

    do i = 1, n_eta
      value_wk = s1c(idx0, i) * pow10_i(is1e(idx0, i))
      derivative_wk = s1dc(idx0, i) * pow10_i(is1de(idx0, i))
      call split_wk_to_double_pair(value_wk, value_hi(i), value_lo(i))
      call split_wk_to_double_pair(derivative_wk, derivative_hi(i), derivative_lo(i))
    end do
  end subroutine oblate_smn_batch_quad_fullsplit

  subroutine oblate_smn_batch_quad_text(m, n, n_eta, normalize, c_text, str_len, eta_text, value_text, value_exp, derivative_text, derivative_exp, status) bind(C, name="oblate_smn_batch_quad_text")
    integer(c_int), value, intent(in) :: m, n, n_eta, normalize, str_len
    character(c_char), intent(in) :: c_text(*), eta_text(*)
    character(c_char), intent(out) :: value_text(*), derivative_text(*)
    integer(c_int), intent(out) :: value_exp(*), derivative_exp(*)
    integer(c_int), intent(out) :: status

    call oblate_smn_batch_quad_text_impl(m, n, n_eta, normalize, c_text, str_len, eta_text, value_text, value_exp, derivative_text, derivative_exp, status)
  end subroutine oblate_smn_batch_quad_text

  subroutine oblate_smn_batch_quad_text_acc(m, n, n_eta, normalize, c_text, str_len, eta_text, value_text, value_exp, derivative_text, derivative_exp, status, accuracy_out) bind(C, name="oblate_smn_batch_quad_text_acc")
    integer(c_int), value, intent(in) :: m, n, n_eta, normalize, str_len
    character(c_char), intent(in) :: c_text(*), eta_text(*)
    character(c_char), intent(out) :: value_text(*), derivative_text(*)
    integer(c_int), intent(out) :: value_exp(*), derivative_exp(*)
    integer(c_int), intent(out) :: status

    integer(c_int), intent(out) :: accuracy_out(*)
    call oblate_smn_batch_quad_text_impl(m, n, n_eta, normalize, c_text, str_len, eta_text, value_text, value_exp, derivative_text, derivative_exp, status, accuracy_out)
  end subroutine oblate_smn_batch_quad_text_acc

  subroutine oblate_smn_batch_quad_text_impl(m, n, n_eta, normalize, c_text, str_len, eta_text, value_text, value_exp, derivative_text, derivative_exp, status, accuracy_out)
    integer(c_int), value, intent(in) :: m, n, n_eta, normalize, str_len
    character(c_char), intent(in) :: c_text(*), eta_text(*)
    character(c_char), intent(out) :: value_text(*), derivative_text(*)
    integer(c_int), intent(out) :: value_exp(*), derivative_exp(*)
    integer(c_int), intent(out) :: status

    integer(c_int), optional, intent(out) :: accuracy_out(*)

    integer(c_int) :: i, idx0, lnum, ioprad, iopang, iopnorm, narg
    integer(c_int) :: off
    real(wk) :: x
    real(wk) :: c_w
    real(wk) :: eta_w, value_wk, derivative_wk
    logical :: ok

    real(wk), allocatable :: arg(:)
    real(wk), allocatable :: r1c(:), r1dc(:), r2c(:), r2dc(:)
    real(wk), allocatable :: s1c(:,:), s1dc(:,:)
    integer(c_int), allocatable :: ir1e(:), ir1de(:), ir2e(:), ir2de(:), naccr(:)
    integer(c_int), allocatable :: is1e(:,:), is1de(:,:), naccs(:,:)
    real(wk), allocatable :: eigout(:)

    status = 0_c_int
    if (n_eta < 1_c_int) then
      status = -4_c_int
      return
    end if
    if (n < 0_c_int) then
      status = -1_c_int
      return
    end if
    if (n < m) then
      status = -2_c_int
      return
    end if

    call decode_real_text(c_text, 1_c_int, str_len, c_w, ok)
    if (.not. ok) then
      status = -5_c_int
      return
    end if

    do i = 1, n_eta
      off = (i - 1_c_int) * str_len + 1_c_int
      call decode_real_text(eta_text, off, str_len, eta_w, ok)
      if (.not. ok) then
        status = -5_c_int
        return
      end if
      if (abs(eta_w) > 1.0_wk .or. .not. is_finite_quad(real(eta_w, rk))) then
        status = -3_c_int
        return
      end if
    end do

    lnum = solver_lnum(c_w, n - m + 1_c_int)
    idx0 = n - m + 1_c_int
    narg = n_eta

    ioprad = 0_c_int
    iopang = 2_c_int
    iopnorm = 0_c_int
    if (normalize /= 0_c_int) iopnorm = 1_c_int
    x = 10.0_wk

    allocate(arg(narg))
    do i = 1, narg
      off = (i - 1_c_int) * str_len + 1_c_int
      call decode_real_text(eta_text, off, str_len, arg(i), ok)
    end do
    allocate(r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum))
    allocate(ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), naccr(lnum))
    naccr = -1_c_int
    allocate(s1c(lnum, narg), s1dc(lnum, narg))
    allocate(is1e(lnum, narg), is1de(lnum, narg), naccs(lnum, narg))
    allocate(eigout(lnum))

    call oblfcn(c_w, m, lnum, ioprad, x, iopang, iopnorm, narg, arg, &
                r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                s1c, is1e, s1dc, is1de, naccs, eigout)

    do i = 1, n_eta
      off = (i - 1_c_int) * str_len + 1_c_int
      value_wk = s1c(idx0, i) * pow10_i(is1e(idx0, i))
      derivative_wk = s1dc(idx0, i) * pow10_i(is1de(idx0, i))
      call encode_real_text(value_text, off, str_len, value_wk)
      call encode_real_text(derivative_text, off, str_len, derivative_wk)
      value_exp(i) = 0_c_int
      derivative_exp(i) = 0_c_int
      if (present(accuracy_out)) accuracy_out(i) = naccs(idx0,i)
    end do
  end subroutine oblate_smn_batch_quad_text_impl

  subroutine oblate_rmn_batch_quad(m, n, c, n_x, xvec, kind, value_re, value_im, deriv_re, deriv_im, status) bind(C, name="oblate_rmn_batch_quad")
    integer(c_int), value, intent(in) :: m, n, n_x, kind
    real(c_double), value, intent(in) :: c
    real(c_double), intent(in) :: xvec(*)
    real(c_double), intent(out) :: value_re(*), value_im(*), deriv_re(*), deriv_im(*)
    integer(c_int), intent(out) :: status(*)

    integer(c_int) :: i, idx0, lnum, ioprad, iopang, iopnorm, narg
    real(wk) :: x
    real(wk) :: r1v, r1d, r2v, r2d
    real(wk) :: c_w

    real(wk), allocatable :: arg(:)
    real(wk), allocatable :: r1c(:), r1dc(:), r2c(:), r2dc(:)
    real(wk), allocatable :: s1c(:,:), s1dc(:,:)
    integer(c_int), allocatable :: ir1e(:), ir1de(:), ir2e(:), ir2de(:), naccr(:)
    integer(c_int), allocatable :: is1e(:,:), is1de(:,:), naccs(:,:)
    real(wk), allocatable :: eigout(:)

    narg = 1_c_int
    iopang = 0_c_int
    iopnorm = 0_c_int

    do i = 1, n_x
      status(i) = 0_c_int
      value_re(i) = 0.0_rk
      value_im(i) = 0.0_rk
      deriv_re(i) = 0.0_rk
      deriv_im(i) = 0.0_rk
    end do

    if (n < 0_c_int) then
      do i = 1, n_x
        status(i) = -1_c_int
      end do
      return
    end if
    if (n < m) then
      do i = 1, n_x
        status(i) = -2_c_int
      end do
      return
    end if
    if (kind < 1_c_int .or. kind > 4_c_int) then
      do i = 1, n_x
        status(i) = -4_c_int
      end do
      return
    end if

    c_w = real(c, wk)

    lnum = solver_lnum(c_w, n - m + 1_c_int)
    idx0 = n - m + 1_c_int

    allocate(arg(1))
    arg(1) = 0.0_wk
    allocate(r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum))
    allocate(ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), naccr(lnum))
    naccr = -1_c_int
    allocate(s1c(lnum, 1), s1dc(lnum, 1))
    allocate(is1e(lnum, 1), is1de(lnum, 1), naccs(lnum, 1))
    allocate(eigout(lnum))

    do i = 1, n_x
      if (.not. is_finite_quad(xvec(i)) .or. xvec(i) < 0.0_rk) then
        status(i) = -3_c_int
        cycle
      end if

      x = real(xvec(i), wk)

      ioprad = 1_c_int
      call oblfcn(c_w, m, lnum, ioprad, x, iopang, iopnorm, narg, arg, &
                  r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                  s1c, is1e, s1dc, is1de, naccs, eigout)
      r1v = r1c(idx0) * pow10_i(ir1e(idx0))
      r1d = r1dc(idx0) * pow10_i(ir1de(idx0))

      if (kind >= 2_c_int) then
        ioprad = 2_c_int
        call oblfcn(c_w, m, lnum, ioprad, x, iopang, iopnorm, narg, arg, &
                    r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                    s1c, is1e, s1dc, is1de, naccs, eigout)
        r2v = r2c(idx0) * pow10_i(ir2e(idx0))
        r2d = r2dc(idx0) * pow10_i(ir2de(idx0))
      else
        r2v = 0.0_wk
        r2d = 0.0_wk
      end if

      select case (kind)
      case (1)
        value_re(i) = r1v
        deriv_re(i) = r1d
      case (2)
        value_re(i) = r2v
        deriv_re(i) = r2d
      case (3)
        value_re(i) = r1v
        value_im(i) = r2v
        deriv_re(i) = r1d
        deriv_im(i) = r2d
      case (4)
        value_re(i) = r1v
        value_im(i) = -r2v
        deriv_re(i) = r1d
        deriv_im(i) = -r2d
      end select
    end do
  end subroutine oblate_rmn_batch_quad

  subroutine oblate_rmn_batch_quad_fullsplit(m, n, c_hi, c_lo, n_x, xvec_hi, xvec_lo, kind, value_re_hi, value_re_lo, value_im_hi, value_im_lo, deriv_re_hi, deriv_re_lo, deriv_im_hi, deriv_im_lo, status) bind(C, name="oblate_rmn_batch_quad_fullsplit")
    integer(c_int), value, intent(in) :: m, n, n_x, kind
    real(c_double), value, intent(in) :: c_hi, c_lo
    real(c_double), intent(in) :: xvec_hi(*), xvec_lo(*)
    real(c_double), intent(out) :: value_re_hi(*), value_re_lo(*), value_im_hi(*), value_im_lo(*)
    real(c_double), intent(out) :: deriv_re_hi(*), deriv_re_lo(*), deriv_im_hi(*), deriv_im_lo(*)
    integer(c_int), intent(out) :: status(*)

    call oblate_rmn_batch_quad_fullsplit_impl(m, n, c_hi, c_lo, n_x, xvec_hi, xvec_lo, kind, value_re_hi, value_re_lo, value_im_hi, value_im_lo, deriv_re_hi, deriv_re_lo, deriv_im_hi, deriv_im_lo, status)
  end subroutine oblate_rmn_batch_quad_fullsplit

  subroutine oblate_rmn_batch_quad_fullsplit_acc(m, n, c_hi, c_lo, n_x, xvec_hi, xvec_lo, kind, value_re_hi, value_re_lo, value_im_hi, value_im_lo, deriv_re_hi, deriv_re_lo, deriv_im_hi, deriv_im_lo, status, accuracy_out) bind(C, name="oblate_rmn_batch_quad_fullsplit_acc")
    integer(c_int), value, intent(in) :: m, n, n_x, kind
    real(c_double), value, intent(in) :: c_hi, c_lo
    real(c_double), intent(in) :: xvec_hi(*), xvec_lo(*)
    real(c_double), intent(out) :: value_re_hi(*), value_re_lo(*), value_im_hi(*), value_im_lo(*)
    real(c_double), intent(out) :: deriv_re_hi(*), deriv_re_lo(*), deriv_im_hi(*), deriv_im_lo(*)
    integer(c_int), intent(out) :: status(*)

    integer(c_int), intent(out) :: accuracy_out(*)
    call oblate_rmn_batch_quad_fullsplit_impl(m, n, c_hi, c_lo, n_x, xvec_hi, xvec_lo, kind, value_re_hi, value_re_lo, value_im_hi, value_im_lo, deriv_re_hi, deriv_re_lo, deriv_im_hi, deriv_im_lo, status, accuracy_out)
  end subroutine oblate_rmn_batch_quad_fullsplit_acc

  subroutine oblate_rmn_batch_quad_fullsplit_impl(m, n, c_hi, c_lo, n_x, xvec_hi, xvec_lo, kind, value_re_hi, value_re_lo, value_im_hi, value_im_lo, deriv_re_hi, deriv_re_lo, deriv_im_hi, deriv_im_lo, status, accuracy_out)
    integer(c_int), value, intent(in) :: m, n, n_x, kind
    real(c_double), value, intent(in) :: c_hi, c_lo
    real(c_double), intent(in) :: xvec_hi(*), xvec_lo(*)
    real(c_double), intent(out) :: value_re_hi(*), value_re_lo(*), value_im_hi(*), value_im_lo(*)
    real(c_double), intent(out) :: deriv_re_hi(*), deriv_re_lo(*), deriv_im_hi(*), deriv_im_lo(*)
    integer(c_int), intent(out) :: status(*)

    integer(c_int), optional, intent(out) :: accuracy_out(*)

    integer(c_int) :: i, idx0, lnum, ioprad, iopang, iopnorm, narg
    real(wk) :: x
    real(wk) :: r1v, r1d, r2v, r2d
    real(wk) :: c_w, xv_w

    real(wk), allocatable :: arg(:)
    real(wk), allocatable :: r1c(:), r1dc(:), r2c(:), r2dc(:)
    real(wk), allocatable :: s1c(:,:), s1dc(:,:)
    integer(c_int), allocatable :: ir1e(:), ir1de(:), ir2e(:), ir2de(:), naccr(:)
    integer(c_int), allocatable :: is1e(:,:), is1de(:,:), naccs(:,:)
    real(wk), allocatable :: eigout(:)

    narg = 1_c_int
    iopang = 0_c_int
    iopnorm = 0_c_int

    do i = 1, n_x
      status(i) = 0_c_int
      value_re_hi(i) = 0.0_rk
      value_re_lo(i) = 0.0_rk
      value_im_hi(i) = 0.0_rk
      value_im_lo(i) = 0.0_rk
      deriv_re_hi(i) = 0.0_rk
      deriv_re_lo(i) = 0.0_rk
      deriv_im_hi(i) = 0.0_rk
      deriv_im_lo(i) = 0.0_rk
    end do

    if (n < 0_c_int) then
      do i = 1, n_x
        status(i) = -1_c_int
      end do
      return
    end if
    if (n < m) then
      do i = 1, n_x
        status(i) = -2_c_int
      end do
      return
    end if
    if (kind < 1_c_int .or. kind > 4_c_int) then
      do i = 1, n_x
        status(i) = -4_c_int
      end do
      return
    end if

    c_w = pair_double_to_wk(c_hi, c_lo)

    lnum = solver_lnum(c_w, n - m + 1_c_int)
    idx0 = n - m + 1_c_int

    allocate(arg(1))
    arg(1) = 0.0_wk
    allocate(r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum))
    allocate(ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), naccr(lnum))
    naccr = -1_c_int
    allocate(s1c(lnum, 1), s1dc(lnum, 1))
    allocate(is1e(lnum, 1), is1de(lnum, 1), naccs(lnum, 1))
    allocate(eigout(lnum))

    do i = 1, n_x
      xv_w = pair_double_to_wk(xvec_hi(i), xvec_lo(i))
      if (.not. is_finite_quad(real(xv_w, rk)) .or. xv_w < 0.0_wk) then
        status(i) = -3_c_int
        cycle
      end if

      x = xv_w

      ioprad = 1_c_int
      call oblfcn(c_w, m, lnum, ioprad, x, iopang, iopnorm, narg, arg, &
                  r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                  s1c, is1e, s1dc, is1de, naccs, eigout)
      r1v = r1c(idx0) * pow10_i(ir1e(idx0))
      r1d = r1dc(idx0) * pow10_i(ir1de(idx0))

      if (kind >= 2_c_int) then
        ioprad = 2_c_int
        call oblfcn(c_w, m, lnum, ioprad, x, iopang, iopnorm, narg, arg, &
                    r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                    s1c, is1e, s1dc, is1de, naccs, eigout)
        r2v = r2c(idx0) * pow10_i(ir2e(idx0))
        r2d = r2dc(idx0) * pow10_i(ir2de(idx0))
      else
        r2v = 0.0_wk
        r2d = 0.0_wk
      end if

      if (present(accuracy_out)) accuracy_out(i) = naccr(idx0)

      select case (kind)
      case (1)
        call split_wk_to_double_pair(r1v, value_re_hi(i), value_re_lo(i))
        call split_wk_to_double_pair(r1d, deriv_re_hi(i), deriv_re_lo(i))
      case (2)
        call split_wk_to_double_pair(r2v, value_re_hi(i), value_re_lo(i))
        call split_wk_to_double_pair(r2d, deriv_re_hi(i), deriv_re_lo(i))
      case (3)
        call split_wk_to_double_pair(r1v, value_re_hi(i), value_re_lo(i))
        call split_wk_to_double_pair(r2v, value_im_hi(i), value_im_lo(i))
        call split_wk_to_double_pair(r1d, deriv_re_hi(i), deriv_re_lo(i))
        call split_wk_to_double_pair(r2d, deriv_im_hi(i), deriv_im_lo(i))
      case (4)
        call split_wk_to_double_pair(r1v, value_re_hi(i), value_re_lo(i))
        call split_wk_to_double_pair(-r2v, value_im_hi(i), value_im_lo(i))
        call split_wk_to_double_pair(r1d, deriv_re_hi(i), deriv_re_lo(i))
        call split_wk_to_double_pair(-r2d, deriv_im_hi(i), deriv_im_lo(i))
      end select
    end do
  end subroutine oblate_rmn_batch_quad_fullsplit_impl

  subroutine oblate_smn_batch_quad_acc(m, n, c, n_eta, eta, normalize, value, derivative, naccs_out, status) bind(C, name="oblate_smn_batch_quad_acc")
    integer(c_int), value, intent(in) :: m, n, n_eta, normalize
    real(c_double), value, intent(in) :: c
    real(c_double), intent(in) :: eta(*)
    real(c_double), intent(out) :: value(*), derivative(*)
    integer(c_int), intent(out) :: naccs_out(*), status

    integer(c_int) :: i, idx0, lnum, ioprad, iopang, iopnorm, narg
    real(wk) :: x
    real(wk) :: c_w

    real(wk), allocatable :: arg(:)
    real(wk), allocatable :: r1c(:), r1dc(:), r2c(:), r2dc(:)
    real(wk), allocatable :: s1c(:,:), s1dc(:,:)
    integer(c_int), allocatable :: ir1e(:), ir1de(:), ir2e(:), ir2de(:), naccr(:)
    integer(c_int), allocatable :: is1e(:,:), is1de(:,:), naccs(:,:)
    real(wk), allocatable :: eigout(:)

    status = 0
    if (n_eta < 1_c_int) then
      status = -4_c_int
      return
    end if
    if (n < 0_c_int) then
      status = -1_c_int
      return
    end if
    if (n < m) then
      status = -2_c_int
      return
    end if

    do i = 1, n_eta
      if (abs(eta(i)) > 1.0_rk .or. .not. is_finite_quad(eta(i))) then
        status = -3_c_int
        return
      end if
    end do

    c_w = real(c, wk)
    lnum = solver_lnum(c_w, n - m + 1_c_int)
    idx0 = n - m + 1_c_int
    narg = n_eta

    ioprad = 0_c_int
    iopang = 2_c_int
    iopnorm = 0_c_int
    if (normalize /= 0_c_int) iopnorm = 1_c_int
    x = 10.0_wk

    allocate(arg(narg))
    arg(:) = real(eta(1:narg), wk)
    allocate(r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum))
    allocate(ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), naccr(lnum))
    naccr = -1_c_int
    allocate(s1c(lnum, narg), s1dc(lnum, narg))
    allocate(is1e(lnum, narg), is1de(lnum, narg), naccs(lnum, narg))
    allocate(eigout(lnum))

    call oblfcn(c_w, m, lnum, ioprad, x, iopang, iopnorm, narg, arg, &
                r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                s1c, is1e, s1dc, is1de, naccs, eigout)

    do i = 1, n_eta
      value(i) = real(s1c(idx0, i) * pow10_i(is1e(idx0, i)), rk)
      derivative(i) = real(s1dc(idx0, i) * pow10_i(is1de(idx0, i)), rk)
      naccs_out(i) = naccs(idx0, i)
    end do
  end subroutine oblate_smn_batch_quad_acc

  subroutine oblate_rmn_batch_quad_acc(m, n, c, n_x, xvec, kind, value_re, value_im, deriv_re, deriv_im, naccr_out, status) bind(C, name="oblate_rmn_batch_quad_acc")
    integer(c_int), value, intent(in) :: m, n, n_x, kind
    real(c_double), value, intent(in) :: c
    real(c_double), intent(in) :: xvec(*)
    real(c_double), intent(out) :: value_re(*), value_im(*), deriv_re(*), deriv_im(*)
    integer(c_int), intent(out) :: naccr_out(*), status(*)

    integer(c_int) :: i, idx0, lnum, ioprad, iopang, iopnorm, narg
    real(wk) :: x
    real(wk) :: r1v, r1d, r2v, r2d
    real(wk) :: c_w

    real(wk), allocatable :: arg(:)
    real(wk), allocatable :: r1c(:), r1dc(:), r2c(:), r2dc(:)
    real(wk), allocatable :: s1c(:,:), s1dc(:,:)
    integer(c_int), allocatable :: ir1e(:), ir1de(:), ir2e(:), ir2de(:), naccr(:)
    integer(c_int), allocatable :: is1e(:,:), is1de(:,:), naccs(:,:)
    real(wk), allocatable :: eigout(:)

    narg = 1_c_int
    iopang = 0_c_int
    iopnorm = 0_c_int

    do i = 1, n_x
      status(i) = 0_c_int
      value_re(i) = 0.0_rk
      value_im(i) = 0.0_rk
      deriv_re(i) = 0.0_rk
      deriv_im(i) = 0.0_rk
      naccr_out(i) = 0_c_int
    end do

    if (n < 0_c_int) then
      do i = 1, n_x
        status(i) = -1_c_int
      end do
      return
    end if
    if (n < m) then
      do i = 1, n_x
        status(i) = -2_c_int
      end do
      return
    end if
    if (kind < 1_c_int .or. kind > 4_c_int) then
      do i = 1, n_x
        status(i) = -4_c_int
      end do
      return
    end if

    c_w = real(c, wk)
    lnum = solver_lnum(c_w, n - m + 1_c_int)
    idx0 = n - m + 1_c_int

    allocate(arg(1))
    arg(1) = 0.0_wk
    allocate(r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum))
    allocate(ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), naccr(lnum))
    naccr = -1_c_int
    allocate(s1c(lnum, 1), s1dc(lnum, 1))
    allocate(is1e(lnum, 1), is1de(lnum, 1), naccs(lnum, 1))
    allocate(eigout(lnum))

    do i = 1, n_x
      if (.not. is_finite_quad(xvec(i)) .or. xvec(i) < 0.0_rk) then
        status(i) = -3_c_int
        cycle
      end if

      x = real(xvec(i), wk)

      ioprad = 1_c_int
      call oblfcn(c_w, m, lnum, ioprad, x, iopang, iopnorm, narg, arg, &
                  r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                  s1c, is1e, s1dc, is1de, naccs, eigout)
      r1v = r1c(idx0) * pow10_i(ir1e(idx0))
      r1d = r1dc(idx0) * pow10_i(ir1de(idx0))

      if (kind >= 2_c_int) then
        ioprad = 2_c_int
        call oblfcn(c_w, m, lnum, ioprad, x, iopang, iopnorm, narg, arg, &
                    r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                    s1c, is1e, s1dc, is1de, naccs, eigout)
        r2v = r2c(idx0) * pow10_i(ir2e(idx0))
        r2d = r2dc(idx0) * pow10_i(ir2de(idx0))
      else
        r2v = 0.0_wk
        r2d = 0.0_wk
      end if

      naccr_out(i) = naccr(idx0)

      select case (kind)
      case (1)
        value_re(i) = real(r1v, rk)
        value_im(i) = 0.0_rk
        deriv_re(i) = real(r1d, rk)
        deriv_im(i) = 0.0_rk
      case (2)
        value_re(i) = real(r2v, rk)
        value_im(i) = 0.0_rk
        deriv_re(i) = real(r2d, rk)
        deriv_im(i) = 0.0_rk
      case (3)
        value_re(i) = real(r1v, rk)
        value_im(i) = real(r2v, rk)
        deriv_re(i) = real(r1d, rk)
        deriv_im(i) = real(r2d, rk)
      case (4)
        value_re(i) = real(r1v, rk)
        value_im(i) = -real(r2v, rk)
        deriv_re(i) = real(r1d, rk)
        deriv_im(i) = -real(r2d, rk)
      end select
    end do
  end subroutine oblate_rmn_batch_quad_acc

  subroutine oblate_eigenvalue_quad(m, n, c, eig, status) bind(C, name="oblate_eigenvalue_quad")
    integer(c_int), value, intent(in) :: m, n
    real(c_double), value, intent(in) :: c
    real(c_double), intent(out) :: eig
    integer(c_int), intent(out) :: status

    integer(c_int) :: lnum, idx0, ioprad, iopang, iopnorm, narg
    real(wk) :: x, c_w
    real(wk), allocatable :: arg(:)
    real(wk), allocatable :: r1c(:), r1dc(:), r2c(:), r2dc(:), eigout(:)
    real(wk), allocatable :: s1c(:,:), s1dc(:,:)
    integer(c_int), allocatable :: ir1e(:), ir1de(:), ir2e(:), ir2de(:), naccr(:)
    integer(c_int), allocatable :: is1e(:,:), is1de(:,:), naccs(:,:)

    status = 0_c_int
    eig = 0.0_rk
    if (n < 0_c_int) then
      status = -1_c_int
      return
    end if
    if (n < m) then
      status = -2_c_int
      return
    end if

    c_w = real(c, wk)
    lnum = solver_lnum(c_w, n - m + 1_c_int)
    idx0 = n - m + 1_c_int
    narg = 1_c_int
    ioprad = 0_c_int
    iopang = 0_c_int
    iopnorm = 0_c_int
    x = 10.0_wk

    allocate(arg(1))
    arg(1) = 0.0_wk
    allocate(r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum), eigout(lnum))
    allocate(ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), naccr(lnum))
    naccr = -1_c_int
    allocate(s1c(lnum, 1), s1dc(lnum, 1))
    allocate(is1e(lnum, 1), is1de(lnum, 1), naccs(lnum, 1))

    call oblfcn(c_w, m, lnum, ioprad, x, iopang, iopnorm, narg, arg, &
                r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                s1c, is1e, s1dc, is1de, naccs, eigout)

    eig = real(eigout(idx0), rk)
  end subroutine oblate_eigenvalue_quad

  subroutine oblate_eigenvalue_quad_fullsplit(m, n, c_hi, c_lo, eig_hi, eig_lo, status) bind(C, name="oblate_eigenvalue_quad_fullsplit")
    integer(c_int), value, intent(in) :: m, n
    real(c_double), value, intent(in) :: c_hi, c_lo
    real(c_double), intent(out) :: eig_hi, eig_lo
    integer(c_int), intent(out) :: status

    integer(c_int) :: lnum, idx0, ioprad, iopang, iopnorm, narg
    real(wk) :: x, c_w, eig_w
    real(wk), allocatable :: arg(:)
    real(wk), allocatable :: r1c(:), r1dc(:), r2c(:), r2dc(:), eigout(:)
    real(wk), allocatable :: s1c(:,:), s1dc(:,:)
    integer(c_int), allocatable :: ir1e(:), ir1de(:), ir2e(:), ir2de(:), naccr(:)
    integer(c_int), allocatable :: is1e(:,:), is1de(:,:), naccs(:,:)

    status = 0_c_int
    eig_hi = 0.0_rk
    eig_lo = 0.0_rk
    if (n < 0_c_int) then
      status = -1_c_int
      return
    end if
    if (n < m) then
      status = -2_c_int
      return
    end if

    c_w = pair_double_to_wk(c_hi, c_lo)
    lnum = solver_lnum(c_w, n - m + 1_c_int)
    idx0 = n - m + 1_c_int
    narg = 1_c_int
    ioprad = 0_c_int
    iopang = 0_c_int
    iopnorm = 0_c_int
    x = 10.0_wk

    allocate(arg(1))
    arg(1) = 0.0_wk
    allocate(r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum), eigout(lnum))
    allocate(ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), naccr(lnum))
    naccr = -1_c_int
    allocate(s1c(lnum, 1), s1dc(lnum, 1))
    allocate(is1e(lnum, 1), is1de(lnum, 1), naccs(lnum, 1))

    call oblfcn(c_w, m, lnum, ioprad, x, iopang, iopnorm, narg, arg, &
                r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                s1c, is1e, s1dc, is1de, naccs, eigout)

    eig_w = eigout(idx0)
    call split_wk_to_double_pair(eig_w, eig_hi, eig_lo)
  end subroutine oblate_eigenvalue_quad_fullsplit

end module oblate_batch_fortran
