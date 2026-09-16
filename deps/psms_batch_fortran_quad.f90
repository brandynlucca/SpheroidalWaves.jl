module psms_batch_fortran
  use, intrinsic :: iso_c_binding
  use prolate_parameters, only: knd
  use prolate_swf, only: profcn
  implicit none

  integer, parameter :: rk = c_double
  integer, parameter :: wk = knd

contains

  subroutine psms_angular_precision_v2() bind(C, name="psms_angular_precision_v2")
    ! Capability marker for derivative-aware truncation and stable Legendre factors.
  end subroutine psms_angular_precision_v2

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

  subroutine psms_smn_batch_quad(m, n, c, n_eta, eta, normalize, value, derivative, status) bind(C, name="psms_smn_batch_quad")
    integer(c_int), value, intent(in) :: m, n, n_eta, normalize
    real(c_double), value, intent(in) :: c
    real(c_double), intent(in) :: eta(*)
    real(c_double), intent(out) :: value(*), derivative(*)
    integer(c_int), intent(out) :: status

    integer(c_int) :: i, idx0, lnum
    integer(c_int) :: ioprad, iopang, iopnorm, narg
    real(wk) :: x1
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
    lnum = n - m + 1_c_int
    narg = n_eta
    idx0 = n - m + 1_c_int

    ioprad = 0_c_int
    iopang = 2_c_int
    iopnorm = 0_c_int
    if (normalize /= 0_c_int) iopnorm = 1_c_int
    x1 = 1.0_wk

    allocate(arg(narg))
    arg(:) = real(eta(1:narg), wk)

    allocate(r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum))
    allocate(ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), naccr(lnum))
    naccr = -1_c_int
    allocate(s1c(lnum, narg), s1dc(lnum, narg))
    allocate(is1e(lnum, narg), is1de(lnum, narg), naccs(lnum, narg))
    allocate(eigout(lnum))

    call profcn(c_w, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
                r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                s1c, is1e, s1dc, is1de, naccs, eigout)

    do i = 1, n_eta
      value(i) = real(s1c(idx0, i) * pow10_i(is1e(idx0, i)), rk)
      derivative(i) = real(s1dc(idx0, i) * pow10_i(is1de(idx0, i)), rk)
    end do
  end subroutine psms_smn_batch_quad

  subroutine psms_smn_batch_quad_split(m, n, c, n_eta, eta, normalize, value_hi, value_lo, derivative_hi, derivative_lo, status) bind(C, name="psms_smn_batch_quad_split")
    integer(c_int), value, intent(in) :: m, n, n_eta, normalize
    real(c_double), value, intent(in) :: c
    real(c_double), intent(in) :: eta(*)
    real(c_double), intent(out) :: value_hi(*), value_lo(*), derivative_hi(*), derivative_lo(*)
    integer(c_int), intent(out) :: status

    integer(c_int) :: i, idx0, lnum
    integer(c_int) :: ioprad, iopang, iopnorm, narg
    real(wk) :: x1
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

    do i = 1, n_eta
      if (abs(eta(i)) > 1.0_rk .or. .not. is_finite_quad(eta(i))) then
        status = -3_c_int
        return
      end if
    end do

    c_w = real(c, wk)
    lnum = n - m + 1_c_int
    narg = n_eta
    idx0 = n - m + 1_c_int

    ioprad = 0_c_int
    iopang = 2_c_int
    iopnorm = 0_c_int
    if (normalize /= 0_c_int) iopnorm = 1_c_int
    x1 = 1.0_wk

    allocate(arg(narg))
    arg(:) = real(eta(1:narg), wk)

    allocate(r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum))
    allocate(ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), naccr(lnum))
    naccr = -1_c_int
    allocate(s1c(lnum, narg), s1dc(lnum, narg))
    allocate(is1e(lnum, narg), is1de(lnum, narg), naccs(lnum, narg))
    allocate(eigout(lnum))

    call profcn(c_w, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
                r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                s1c, is1e, s1dc, is1de, naccs, eigout)

    do i = 1, n_eta
      value_wk = s1c(idx0, i) * pow10_i(is1e(idx0, i))
      derivative_wk = s1dc(idx0, i) * pow10_i(is1de(idx0, i))
      call split_wk_to_double_pair(value_wk, value_hi(i), value_lo(i))
      call split_wk_to_double_pair(derivative_wk, derivative_hi(i), derivative_lo(i))
    end do
  end subroutine psms_smn_batch_quad_split

  subroutine psms_smn_batch_quad_fullsplit(m, n, c_hi, c_lo, n_eta, eta_hi, eta_lo, normalize, value_hi, value_lo, derivative_hi, derivative_lo, status) bind(C, name="psms_smn_batch_quad_fullsplit")
    integer(c_int), value, intent(in) :: m, n, n_eta, normalize
    real(c_double), value, intent(in) :: c_hi, c_lo
    real(c_double), intent(in) :: eta_hi(*), eta_lo(*)
    real(c_double), intent(out) :: value_hi(*), value_lo(*), derivative_hi(*), derivative_lo(*)
    integer(c_int), intent(out) :: status

    integer(c_int) :: i, idx0, lnum
    integer(c_int) :: ioprad, iopang, iopnorm, narg
    real(wk) :: x1
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

    lnum = n - m + 1_c_int
    narg = n_eta
    idx0 = n - m + 1_c_int

    ioprad = 0_c_int
    iopang = 2_c_int
    iopnorm = 0_c_int
    if (normalize /= 0_c_int) iopnorm = 1_c_int
    x1 = 1.0_wk

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

    call profcn(c_w, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
                r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                s1c, is1e, s1dc, is1de, naccs, eigout)

    do i = 1, n_eta
      value_wk = s1c(idx0, i) * pow10_i(is1e(idx0, i))
      derivative_wk = s1dc(idx0, i) * pow10_i(is1de(idx0, i))
      call split_wk_to_double_pair(value_wk, value_hi(i), value_lo(i))
      call split_wk_to_double_pair(derivative_wk, derivative_hi(i), derivative_lo(i))
    end do
  end subroutine psms_smn_batch_quad_fullsplit

  subroutine psms_smn_batch_quad_text(m, n, n_eta, normalize, c_text, str_len, eta_text, value_text, value_exp, derivative_text, derivative_exp, status) bind(C, name="psms_smn_batch_quad_text")
    integer(c_int), value, intent(in) :: m, n, n_eta, normalize, str_len
    character(c_char), intent(in) :: c_text(*), eta_text(*)
    character(c_char), intent(out) :: value_text(*), derivative_text(*)
    integer(c_int), intent(out) :: value_exp(*), derivative_exp(*)
    integer(c_int), intent(out) :: status

    call psms_smn_batch_quad_text_impl(m, n, n_eta, normalize, c_text, str_len, eta_text, value_text, value_exp, derivative_text, derivative_exp, status)
  end subroutine psms_smn_batch_quad_text

  subroutine psms_smn_batch_quad_text_acc(m, n, n_eta, normalize, c_text, str_len, eta_text, value_text, value_exp, derivative_text, derivative_exp, status, accuracy_out) bind(C, name="psms_smn_batch_quad_text_acc")
    integer(c_int), value, intent(in) :: m, n, n_eta, normalize, str_len
    character(c_char), intent(in) :: c_text(*), eta_text(*)
    character(c_char), intent(out) :: value_text(*), derivative_text(*)
    integer(c_int), intent(out) :: value_exp(*), derivative_exp(*)
    integer(c_int), intent(out) :: status

    integer(c_int), intent(out) :: accuracy_out(*)
    call psms_smn_batch_quad_text_impl(m, n, n_eta, normalize, c_text, str_len, eta_text, value_text, value_exp, derivative_text, derivative_exp, status, accuracy_out)
  end subroutine psms_smn_batch_quad_text_acc

  subroutine psms_smn_batch_quad_text_impl(m, n, n_eta, normalize, c_text, str_len, eta_text, value_text, value_exp, derivative_text, derivative_exp, status, accuracy_out)
    integer(c_int), value, intent(in) :: m, n, n_eta, normalize, str_len
    character(c_char), intent(in) :: c_text(*), eta_text(*)
    character(c_char), intent(out) :: value_text(*), derivative_text(*)
    integer(c_int), intent(out) :: value_exp(*), derivative_exp(*)
    integer(c_int), intent(out) :: status

    integer(c_int), optional, intent(out) :: accuracy_out(*)

    integer(c_int) :: i, idx0, lnum
    integer(c_int) :: ioprad, iopang, iopnorm, narg
    integer(c_int) :: off
    real(wk) :: x1
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

    lnum = n - m + 1_c_int
    narg = n_eta
    idx0 = n - m + 1_c_int

    ioprad = 0_c_int
    iopang = 2_c_int
    iopnorm = 0_c_int
    if (normalize /= 0_c_int) iopnorm = 1_c_int
    x1 = 1.0_wk

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

    call profcn(c_w, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
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
  end subroutine psms_smn_batch_quad_text_impl

  subroutine psms_rmn_batch_quad(m, n, c, n_xi, xi, kind, value_re, value_im, deriv_re, deriv_im, status) bind(C, name="psms_rmn_batch_quad")
    integer(c_int), value, intent(in) :: m, n, n_xi, kind
    real(c_double), value, intent(in) :: c
    real(c_double), intent(in) :: xi(*)
    real(c_double), intent(out) :: value_re(*), value_im(*), deriv_re(*), deriv_im(*)
    integer(c_int), intent(out) :: status(*)

    integer(c_int) :: i, idx0, lnum, ioprad, iopang, iopnorm, narg
    real(wk) :: x1, c_w

    real(wk), allocatable :: arg(:)
    real(wk), allocatable :: r1c(:), r1dc(:), r2c(:), r2dc(:)
    real(wk), allocatable :: s1c(:,:), s1dc(:,:)
    integer(c_int), allocatable :: ir1e(:), ir1de(:), ir2e(:), ir2de(:), naccr(:)
    integer(c_int), allocatable :: is1e(:,:), is1de(:,:), naccs(:,:)
    real(wk), allocatable :: eigout(:)

    real(wk) :: r1v, r1d, r2v, r2d

    narg = 1_c_int
    iopang = 0_c_int
    iopnorm = 0_c_int
    c_w = real(c, wk)

    do i = 1, n_xi
      status(i) = 0_c_int
      value_re(i) = 0.0_rk
      value_im(i) = 0.0_rk
      deriv_re(i) = 0.0_rk
      deriv_im(i) = 0.0_rk
    end do

    if (n < 0_c_int) then
      do i = 1, n_xi
        status(i) = -1_c_int
      end do
      return
    end if
    if (n < m) then
      do i = 1, n_xi
        status(i) = -2_c_int
      end do
      return
    end if
    if (kind < 1_c_int .or. kind > 4_c_int) then
      do i = 1, n_xi
        status(i) = -4_c_int
      end do
      return
    end if

    lnum = n - m + 1_c_int
    idx0 = n - m + 1_c_int

    allocate(arg(1))
    arg(1) = 1.0_wk
    allocate(r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum))
    allocate(ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), naccr(lnum))
    naccr = -1_c_int
    allocate(s1c(lnum, 1), s1dc(lnum, 1))
    allocate(is1e(lnum, 1), is1de(lnum, 1), naccs(lnum, 1))
    allocate(eigout(lnum))

    do i = 1, n_xi
      if (.not. is_finite_quad(xi(i)) .or. abs(xi(i)) < 1.0_rk) then
        status(i) = -3_c_int
        cycle
      end if

      ioprad = 1_c_int
      x1 = real(xi(i) - 1.0_rk, wk)
      call profcn(c_w, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
                  r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                  s1c, is1e, s1dc, is1de, naccs, eigout)
      r1v = r1c(idx0) * pow10_i(ir1e(idx0))
      r1d = r1dc(idx0) * pow10_i(ir1de(idx0))

      if (kind >= 2_c_int) then
        ioprad = 2_c_int
        call profcn(c_w, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
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
  end subroutine psms_rmn_batch_quad

  subroutine psms_rmn_batch_quad_fullsplit(m, n, c_hi, c_lo, n_xi, xi_hi, xi_lo, kind, value_re_hi, value_re_lo, value_im_hi, value_im_lo, deriv_re_hi, deriv_re_lo, deriv_im_hi, deriv_im_lo, status) bind(C, name="psms_rmn_batch_quad_fullsplit")
    integer(c_int), value, intent(in) :: m, n, n_xi, kind
    real(c_double), value, intent(in) :: c_hi, c_lo
    real(c_double), intent(in) :: xi_hi(*), xi_lo(*)
    real(c_double), intent(out) :: value_re_hi(*), value_re_lo(*), value_im_hi(*), value_im_lo(*)
    real(c_double), intent(out) :: deriv_re_hi(*), deriv_re_lo(*), deriv_im_hi(*), deriv_im_lo(*)
    integer(c_int), intent(out) :: status(*)

    call psms_rmn_batch_quad_fullsplit_impl(m, n, c_hi, c_lo, n_xi, xi_hi, xi_lo, kind, value_re_hi, value_re_lo, value_im_hi, value_im_lo, deriv_re_hi, deriv_re_lo, deriv_im_hi, deriv_im_lo, status)
  end subroutine psms_rmn_batch_quad_fullsplit

  subroutine psms_rmn_batch_quad_fullsplit_acc(m, n, c_hi, c_lo, n_xi, xi_hi, xi_lo, kind, value_re_hi, value_re_lo, value_im_hi, value_im_lo, deriv_re_hi, deriv_re_lo, deriv_im_hi, deriv_im_lo, status, accuracy_out) bind(C, name="psms_rmn_batch_quad_fullsplit_acc")
    integer(c_int), value, intent(in) :: m, n, n_xi, kind
    real(c_double), value, intent(in) :: c_hi, c_lo
    real(c_double), intent(in) :: xi_hi(*), xi_lo(*)
    real(c_double), intent(out) :: value_re_hi(*), value_re_lo(*), value_im_hi(*), value_im_lo(*)
    real(c_double), intent(out) :: deriv_re_hi(*), deriv_re_lo(*), deriv_im_hi(*), deriv_im_lo(*)
    integer(c_int), intent(out) :: status(*)

    integer(c_int), intent(out) :: accuracy_out(*)
    call psms_rmn_batch_quad_fullsplit_impl(m, n, c_hi, c_lo, n_xi, xi_hi, xi_lo, kind, value_re_hi, value_re_lo, value_im_hi, value_im_lo, deriv_re_hi, deriv_re_lo, deriv_im_hi, deriv_im_lo, status, accuracy_out)
  end subroutine psms_rmn_batch_quad_fullsplit_acc

  ! The offset ABI receives x-1, avoiding loss of its significant digits when
  ! a coordinate extremely close to one crosses the C interface.
  subroutine psms_rmn_batch_quad_offset_acc(m, n, c_hi, c_lo, n_xi, x1_hi, x1_lo, kind, value_re_hi, value_re_lo, value_im_hi, value_im_lo, deriv_re_hi, deriv_re_lo, deriv_im_hi, deriv_im_lo, status, accuracy_out) bind(C, name="psms_rmn_batch_quad_offset_acc")
    integer(c_int), value, intent(in) :: m, n, n_xi, kind
    real(c_double), value, intent(in) :: c_hi, c_lo
    real(c_double), intent(in) :: x1_hi(*), x1_lo(*)
    real(c_double), intent(out) :: value_re_hi(*), value_re_lo(*), value_im_hi(*), value_im_lo(*)
    real(c_double), intent(out) :: deriv_re_hi(*), deriv_re_lo(*), deriv_im_hi(*), deriv_im_lo(*)
    integer(c_int), intent(out) :: status(*), accuracy_out(*)

    call psms_rmn_batch_quad_fullsplit_impl(m, n, c_hi, c_lo, n_xi, x1_hi, x1_lo, kind, value_re_hi, value_re_lo, value_im_hi, value_im_lo, deriv_re_hi, deriv_re_lo, deriv_im_hi, deriv_im_lo, status, accuracy_out, .true.)
  end subroutine psms_rmn_batch_quad_offset_acc

  subroutine psms_rmn_batch_quad_fullsplit_impl(m, n, c_hi, c_lo, n_xi, xi_hi, xi_lo, kind, value_re_hi, value_re_lo, value_im_hi, value_im_lo, deriv_re_hi, deriv_re_lo, deriv_im_hi, deriv_im_lo, status, accuracy_out, offset_input)
    integer(c_int), value, intent(in) :: m, n, n_xi, kind
    real(c_double), value, intent(in) :: c_hi, c_lo
    real(c_double), intent(in) :: xi_hi(*), xi_lo(*)
    real(c_double), intent(out) :: value_re_hi(*), value_re_lo(*), value_im_hi(*), value_im_lo(*)
    real(c_double), intent(out) :: deriv_re_hi(*), deriv_re_lo(*), deriv_im_hi(*), deriv_im_lo(*)
    integer(c_int), intent(out) :: status(*)

    integer(c_int), optional, intent(out) :: accuracy_out(*)
    logical, optional, intent(in) :: offset_input
    logical :: offset_mode

    integer(c_int) :: i, idx0, lnum, ioprad, iopang, iopnorm, narg
    real(wk) :: x1, c_w, xi_w

    real(wk), allocatable :: arg(:)
    real(wk), allocatable :: r1c(:), r1dc(:), r2c(:), r2dc(:)
    real(wk), allocatable :: s1c(:,:), s1dc(:,:)
    integer(c_int), allocatable :: ir1e(:), ir1de(:), ir2e(:), ir2de(:), naccr(:)
    integer(c_int), allocatable :: is1e(:,:), is1de(:,:), naccs(:,:)
    real(wk), allocatable :: eigout(:)

    real(wk) :: r1v, r1d, r2v, r2d

    narg = 1_c_int
    iopang = 0_c_int
    iopnorm = 0_c_int
    c_w = pair_double_to_wk(c_hi, c_lo)
    offset_mode = .false.
    if (present(offset_input)) offset_mode = offset_input

    do i = 1, n_xi
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
      do i = 1, n_xi
        status(i) = -1_c_int
      end do
      return
    end if
    if (n < m) then
      do i = 1, n_xi
        status(i) = -2_c_int
      end do
      return
    end if
    if (kind < 1_c_int .or. kind > 4_c_int) then
      do i = 1, n_xi
        status(i) = -4_c_int
      end do
      return
    end if

    lnum = n - m + 1_c_int
    idx0 = n - m + 1_c_int

    allocate(arg(1))
    arg(1) = 1.0_wk
    allocate(r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum))
    allocate(ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), naccr(lnum))
    naccr = -1_c_int
    allocate(s1c(lnum, 1), s1dc(lnum, 1))
    allocate(is1e(lnum, 1), is1de(lnum, 1), naccs(lnum, 1))
    allocate(eigout(lnum))

    do i = 1, n_xi
      xi_w = pair_double_to_wk(xi_hi(i), xi_lo(i))
      if (offset_mode) then
        x1 = xi_w
        xi_w = x1 + 1.0_wk
      else
        x1 = xi_w - 1.0_wk
      end if
      if (.not. is_finite_quad(real(xi_w, rk)) .or. abs(xi_w) < 1.0_wk) then
        status(i) = -3_c_int
        cycle
      end if

      ioprad = 1_c_int
      call profcn(c_w, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
                  r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                  s1c, is1e, s1dc, is1de, naccs, eigout)
      r1v = r1c(idx0) * pow10_i(ir1e(idx0))
      r1d = r1dc(idx0) * pow10_i(ir1de(idx0))

      if (kind >= 2_c_int) then
        ioprad = 2_c_int
        call profcn(c_w, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
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
  end subroutine psms_rmn_batch_quad_fullsplit_impl

  subroutine psms_smn_batch_quad_acc(m, n, c, n_eta, eta, normalize, value, derivative, naccs_out, status) bind(C, name="psms_smn_batch_quad_acc")
    integer(c_int), value, intent(in) :: m, n, n_eta, normalize
    real(c_double), value, intent(in) :: c
    real(c_double), intent(in) :: eta(*)
    real(c_double), intent(out) :: value(*), derivative(*)
    integer(c_int), intent(out) :: naccs_out(*), status

    integer(c_int) :: i, idx0, lnum, ioprad, iopang, iopnorm, narg
    real(wk) :: x1, c_w

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
    lnum = n - m + 1_c_int
    idx0 = n - m + 1_c_int
    narg = n_eta

    ioprad = 0_c_int
    iopang = 2_c_int
    iopnorm = 0_c_int
    if (normalize /= 0_c_int) iopnorm = 1_c_int
    x1 = 1.0_wk

    allocate(arg(narg))
    arg(:) = real(eta(1:narg), wk)
    allocate(r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum))
    allocate(ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), naccr(lnum))
    naccr = -1_c_int
    allocate(s1c(lnum, narg), s1dc(lnum, narg))
    allocate(is1e(lnum, narg), is1de(lnum, narg), naccs(lnum, narg))
    allocate(eigout(lnum))

    call profcn(c_w, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
                r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                s1c, is1e, s1dc, is1de, naccs, eigout)

    do i = 1, n_eta
      value(i) = real(s1c(idx0, i) * pow10_i(is1e(idx0, i)), rk)
      derivative(i) = real(s1dc(idx0, i) * pow10_i(is1de(idx0, i)), rk)
      naccs_out(i) = naccs(idx0, i)
    end do
  end subroutine psms_smn_batch_quad_acc

  subroutine psms_rmn_batch_quad_acc(m, n, c, n_xi, xi, kind, value_re, value_im, deriv_re, deriv_im, naccr_out, status) bind(C, name="psms_rmn_batch_quad_acc")
    integer(c_int), value, intent(in) :: m, n, n_xi, kind
    real(c_double), value, intent(in) :: c
    real(c_double), intent(in) :: xi(*)
    real(c_double), intent(out) :: value_re(*), value_im(*), deriv_re(*), deriv_im(*)
    integer(c_int), intent(out) :: naccr_out(*), status(*)

    integer(c_int) :: i, idx0, lnum, ioprad, iopang, iopnorm, narg
    real(wk) :: x1, c_w

    real(wk), allocatable :: arg(:)
    real(wk), allocatable :: r1c(:), r1dc(:), r2c(:), r2dc(:)
    real(wk), allocatable :: s1c(:,:), s1dc(:,:)
    integer(c_int), allocatable :: ir1e(:), ir1de(:), ir2e(:), ir2de(:), naccr(:)
    integer(c_int), allocatable :: is1e(:,:), is1de(:,:), naccs(:,:)
    real(wk), allocatable :: eigout(:)

    real(wk) :: r1v, r1d, r2v, r2d

    narg = 1_c_int
    iopang = 0_c_int
    iopnorm = 0_c_int
    c_w = real(c, wk)

    do i = 1, n_xi
      status(i) = 0_c_int
      value_re(i) = 0.0_rk
      value_im(i) = 0.0_rk
      deriv_re(i) = 0.0_rk
      deriv_im(i) = 0.0_rk
      naccr_out(i) = 0_c_int
    end do

    if (n < 0_c_int) then
      do i = 1, n_xi
        status(i) = -1_c_int
      end do
      return
    end if
    if (n < m) then
      do i = 1, n_xi
        status(i) = -2_c_int
      end do
      return
    end if
    if (kind < 1_c_int .or. kind > 4_c_int) then
      do i = 1, n_xi
        status(i) = -4_c_int
      end do
      return
    end if

    lnum = n - m + 1_c_int
    idx0 = n - m + 1_c_int

    allocate(arg(1))
    arg(1) = 1.0_wk
    allocate(r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum))
    allocate(ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), naccr(lnum))
    naccr = -1_c_int
    allocate(s1c(lnum, 1), s1dc(lnum, 1))
    allocate(is1e(lnum, 1), is1de(lnum, 1), naccs(lnum, 1))
    allocate(eigout(lnum))

    do i = 1, n_xi
      if (.not. is_finite_quad(xi(i)) .or. abs(xi(i)) < 1.0_rk) then
        status(i) = -3_c_int
        cycle
      end if

      ioprad = 1_c_int
      x1 = real(xi(i) - 1.0_rk, wk)
      call profcn(c_w, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
                  r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                  s1c, is1e, s1dc, is1de, naccs, eigout)
      r1v = r1c(idx0) * pow10_i(ir1e(idx0))
      r1d = r1dc(idx0) * pow10_i(ir1de(idx0))

      if (kind >= 2_c_int) then
        ioprad = 2_c_int
        call profcn(c_w, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
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
  end subroutine psms_rmn_batch_quad_acc

  subroutine psms_eigenvalue_quad(m, n, c, eig, status) bind(C, name="psms_eigenvalue_quad")
    integer(c_int), value, intent(in) :: m, n
    real(c_double), value, intent(in) :: c
    real(c_double), intent(out) :: eig
    integer(c_int), intent(out) :: status

    integer(c_int) :: lnum, idx0
    integer(c_int) :: ioprad, iopang, iopnorm, narg
    real(wk) :: x1, c_w
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
    lnum = n - m + 1_c_int
    idx0 = lnum
    narg = 1_c_int
    ioprad = 0_c_int
    iopang = 0_c_int
    iopnorm = 0_c_int
    x1 = 1.0_wk

    allocate(arg(1))
    arg(1) = 0.0_wk
    allocate(r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum), eigout(lnum))
    allocate(ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), naccr(lnum))
    naccr = -1_c_int
    allocate(s1c(lnum, 1), s1dc(lnum, 1))
    allocate(is1e(lnum, 1), is1de(lnum, 1), naccs(lnum, 1))

    call profcn(c_w, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
                r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                s1c, is1e, s1dc, is1de, naccs, eigout)

    eig = real(eigout(idx0), rk)
  end subroutine psms_eigenvalue_quad

  subroutine psms_eigenvalue_quad_fullsplit(m, n, c_hi, c_lo, eig_hi, eig_lo, status) bind(C, name="psms_eigenvalue_quad_fullsplit")
    integer(c_int), value, intent(in) :: m, n
    real(c_double), value, intent(in) :: c_hi, c_lo
    real(c_double), intent(out) :: eig_hi, eig_lo
    integer(c_int), intent(out) :: status

    integer(c_int) :: lnum, idx0
    integer(c_int) :: ioprad, iopang, iopnorm, narg
    real(wk) :: x1, c_w, eig_w
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
    lnum = n - m + 1_c_int
    idx0 = lnum
    narg = 1_c_int
    ioprad = 0_c_int
    iopang = 0_c_int
    iopnorm = 0_c_int
    x1 = 1.0_wk

    allocate(arg(1))
    arg(1) = 0.0_wk
    allocate(r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum), eigout(lnum))
    allocate(ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), naccr(lnum))
    naccr = -1_c_int
    allocate(s1c(lnum, 1), s1dc(lnum, 1))
    allocate(is1e(lnum, 1), is1de(lnum, 1), naccs(lnum, 1))

    call profcn(c_w, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
                r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                s1c, is1e, s1dc, is1de, naccs, eigout)

    eig_w = eigout(idx0)
    call split_wk_to_double_pair(eig_w, eig_hi, eig_lo)
  end subroutine psms_eigenvalue_quad_fullsplit

  subroutine psms_smn_degrees_quad_text(m, n_min, n, n_eta, normalize, c_text, str_len, eta_text, value_text, value_exp, derivative_text, derivative_exp, status) bind(C, name="psms_smn_degrees_quad_text")
    integer(c_int), value, intent(in) :: m, n_min, n, n_eta, normalize, str_len
    character(c_char), intent(in) :: c_text(*), eta_text(*)
    character(c_char), intent(out) :: value_text(*), derivative_text(*)
    integer(c_int), intent(out) :: value_exp(*), derivative_exp(*)
    integer(c_int), intent(out) :: status

    integer(c_int) :: i, j, idx, idx0, lnum
    integer(c_int) :: ioprad, iopang, iopnorm, narg
    integer(c_int) :: off
    real(wk) :: x1
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
    if (m < 0 .or. n_min < m .or. n < n_min) then
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

    lnum = n - m + 1_c_int
    narg = n_eta
    idx0 = n - m + 1_c_int

    ioprad = 0_c_int
    iopang = 2_c_int
    iopnorm = 0_c_int
    if (normalize /= 0_c_int) iopnorm = 1_c_int
    x1 = 1.0_wk

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

    call profcn(c_w, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
                r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                s1c, is1e, s1dc, is1de, naccs, eigout)

    do j = n_min-m+1, lnum
      do i = 1, n_eta
        idx = (j-(n_min-m+1))*n_eta+i
        off = (idx-1)*str_len+1
        call encode_real_text(value_text, off, str_len, s1c(j,i))
        call encode_real_text(derivative_text, off, str_len, s1dc(j,i))
        value_exp(idx) = is1e(j,i)
        derivative_exp(idx) = is1de(j,i)
      end do
    end do
  end subroutine psms_smn_degrees_quad_text

  ! Degree-range radial evaluation. Layout: four real channels per point,
  ! then points, then degrees. Preserve quad mantissas and exponents.
  subroutine psms_rmn_degrees_quad_text(m, n_min, n_max, kind, c_text, n_x, x_text, str_len, output, exponents, status) bind(C, name="psms_rmn_degrees_quad_text")
    integer(c_int), value, intent(in) :: m, n_min, n_max, kind, n_x, str_len
    character(c_char), intent(in) :: c_text(*), x_text(*)
    character(c_char), intent(out) :: output(*)
    integer(c_int), intent(out) :: exponents(*), status
    integer(c_int) :: lnum, li, j, idx, ioprad
    real(wk) :: c_w, x_w
    logical :: ok
    real(wk), allocatable :: arg(:), r1(:), dr1(:), r2(:), dr2(:), angular(:,:), dangular(:,:), eig(:)
    integer(c_int), allocatable :: er1(:), edr1(:), er2(:), edr2(:), ar(:), ea(:,:), eda(:,:), aa(:,:)
    status = 0
    if (m < 0 .or. n_min < m .or. n_max < n_min .or. n_x < 1 .or. kind < 1 .or. kind > 4 .or. str_len < 70) then
      status = -1
      return
    end if
    call decode_real_text(c_text, 1_c_int, str_len, c_w, ok)
    if (.not. ok .or. .not. (c_w > 0.0_wk)) then
      status = -2
      return
    end if
    lnum = n_max-m+1
    allocate(arg(1), r1(lnum), dr1(lnum), r2(lnum), dr2(lnum), eig(lnum))
    allocate(er1(lnum), edr1(lnum), er2(lnum), edr2(lnum), ar(lnum))
    allocate(angular(lnum,1), dangular(lnum,1), ea(lnum,1), eda(lnum,1), aa(lnum,1))
    arg = 0.0_wk
    ioprad = 2
    if (kind == 1) ioprad = 1
    do j = 1, n_x
      call decode_real_text(x_text, (j-1)*str_len+1, str_len, x_w, ok)
      if (.not. ok .or. .not. (x_w > 1.0_wk)) then
        status = -3
        return
      end if
      call profcn(c_w, m, lnum, ioprad, x_w-1.0_wk, 0_c_int, 0_c_int, 1_c_int, arg, &
                  r1, er1, dr1, edr1, r2, er2, dr2, edr2, ar, &
                  angular, ea, dangular, eda, aa, eig)
      if (kind == 1) then
        r2 = 0.0_wk
        dr2 = 0.0_wk
        er2 = 0
        edr2 = 0
      end if
      do li = n_min-m+1, lnum
        idx = ((li-(n_min-m+1))*n_x+j-1)*4+1
        call encode_real_text(output, (idx-1)*str_len+1, str_len, r1(li))
        exponents(idx) = er1(li)
        idx = idx+1
        call encode_real_text(output, (idx-1)*str_len+1, str_len, dr1(li))
        exponents(idx) = edr1(li)
        idx = idx+1
        call encode_real_text(output, (idx-1)*str_len+1, str_len, r2(li))
        exponents(idx) = er2(li)
        idx = idx+1
        call encode_real_text(output, (idx-1)*str_len+1, str_len, dr2(li))
        exponents(idx) = edr2(li)
      end do
    end do
  end subroutine psms_rmn_degrees_quad_text

end module psms_batch_fortran
