module psms_batch_fortran
  use, intrinsic :: iso_c_binding
  use param, only: knd
  use prolate_swf, only: profcn
  implicit none

  integer, parameter :: rk = c_double
  integer, parameter :: wk = knd

contains

  pure real(wk) function pow10_i(exp10) result(v)
    integer(c_int), intent(in) :: exp10
    v = 10.0_wk ** real(exp10, wk)
  end function pow10_i

  pure logical function is_finite_r16(x) result(ok)
    real(rk), intent(in) :: x
    ok = (x == x) .and. (abs(x) < huge(x))
  end function is_finite_r16

  subroutine psms_smn_batch_r16(m, n, c, n_eta, eta, normalize, value, derivative, status) bind(C, name="psms_smn_batch_r16")
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
      if (abs(eta(i)) > 1.0_rk .or. .not. is_finite_r16(eta(i))) then
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
    allocate(s1c(lnum, narg), s1dc(lnum, narg))
    allocate(is1e(lnum, narg), is1de(lnum, narg), naccs(lnum, narg))

    call profcn(c_w, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
                r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                s1c, is1e, s1dc, is1de, naccs)

    do i = 1, n_eta
      value(i) = real(s1c(idx0, i) * pow10_i(is1e(idx0, i)), rk)
      derivative(i) = real(s1dc(idx0, i) * pow10_i(is1de(idx0, i)), rk)
    end do
  end subroutine psms_smn_batch_r16

  subroutine psms_rmn_batch_r16(m, n, c, n_xi, xi, kind, value_re, value_im, deriv_re, deriv_im, status) bind(C, name="psms_rmn_batch_r16")
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
    allocate(s1c(lnum, 1), s1dc(lnum, 1))
    allocate(is1e(lnum, 1), is1de(lnum, 1), naccs(lnum, 1))

    do i = 1, n_xi
      if (.not. is_finite_r16(xi(i)) .or. abs(xi(i)) < 1.0_rk) then
        status(i) = -3_c_int
        cycle
      end if

      ioprad = 1_c_int
      x1 = real(xi(i) - 1.0_rk, wk)
      call profcn(c_w, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
                  r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                  s1c, is1e, s1dc, is1de, naccs)
      r1v = r1c(idx0) * pow10_i(ir1e(idx0))
      r1d = r1dc(idx0) * pow10_i(ir1de(idx0))

      if (kind >= 2_c_int) then
        ioprad = 2_c_int
        call profcn(c_w, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
                    r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                    s1c, is1e, s1dc, is1de, naccs)
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
  end subroutine psms_rmn_batch_r16

  subroutine psms_smn_batch_r16_acc(m, n, c, n_eta, eta, normalize, value, derivative, naccs_out, status) bind(C, name="psms_smn_batch_r16_acc")
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
      if (abs(eta(i)) > 1.0_rk .or. .not. is_finite_r16(eta(i))) then
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
    allocate(s1c(lnum, narg), s1dc(lnum, narg))
    allocate(is1e(lnum, narg), is1de(lnum, narg), naccs(lnum, narg))

    call profcn(c_w, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
                r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                s1c, is1e, s1dc, is1de, naccs)

    do i = 1, n_eta
      value(i) = real(s1c(idx0, i) * pow10_i(is1e(idx0, i)), rk)
      derivative(i) = real(s1dc(idx0, i) * pow10_i(is1de(idx0, i)), rk)
      naccs_out(i) = naccs(idx0, i)
    end do
  end subroutine psms_smn_batch_r16_acc

  subroutine psms_rmn_batch_r16_acc(m, n, c, n_xi, xi, kind, value_re, value_im, deriv_re, deriv_im, naccr_out, status) bind(C, name="psms_rmn_batch_r16_acc")
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
    allocate(s1c(lnum, 1), s1dc(lnum, 1))
    allocate(is1e(lnum, 1), is1de(lnum, 1), naccs(lnum, 1))

    do i = 1, n_xi
      if (.not. is_finite_r16(xi(i)) .or. abs(xi(i)) < 1.0_rk) then
        status(i) = -3_c_int
        cycle
      end if

      ioprad = 1_c_int
      x1 = real(xi(i) - 1.0_rk, wk)
      call profcn(c_w, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
                  r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                  s1c, is1e, s1dc, is1de, naccs)
      r1v = r1c(idx0) * pow10_i(ir1e(idx0))
      r1d = r1dc(idx0) * pow10_i(ir1de(idx0))

      if (kind >= 2_c_int) then
        ioprad = 2_c_int
        call profcn(c_w, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
                    r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                    s1c, is1e, s1dc, is1de, naccs)
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
  end subroutine psms_rmn_batch_r16_acc

end module psms_batch_fortran
