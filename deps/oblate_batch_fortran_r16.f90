module oblate_batch_fortran
  use, intrinsic :: iso_c_binding
  use param, only: knd
  use oblate_swf, only: oblfcn
  implicit none

  integer, parameter :: rk = c_double
  integer, parameter :: wk = knd

contains

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

  pure logical function is_finite_r16(x) result(ok)
    real(rk), intent(in) :: x
    ok = (x == x) .and. (abs(x) < huge(x))
  end function is_finite_r16

  subroutine oblate_smn_batch_r16(m, n, c, n_eta, eta, normalize, value, derivative, status) bind(C, name="oblate_smn_batch_r16")
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
      if (abs(eta(i)) > 1.0_rk .or. .not. is_finite_r16(eta(i))) then
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
    allocate(s1c(lnum, narg), s1dc(lnum, narg))
    allocate(is1e(lnum, narg), is1de(lnum, narg), naccs(lnum, narg))

    call oblfcn(c_w, m, lnum, ioprad, x, iopang, iopnorm, narg, arg, &
                r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                s1c, is1e, s1dc, is1de, naccs)

    do i = 1, n_eta
      value(i) = s1c(idx0, i) * pow10_i(is1e(idx0, i))
      derivative(i) = s1dc(idx0, i) * pow10_i(is1de(idx0, i))
    end do
  end subroutine oblate_smn_batch_r16

  subroutine oblate_rmn_batch_r16(m, n, c, n_x, xvec, kind, value_re, value_im, deriv_re, deriv_im, status) bind(C, name="oblate_rmn_batch_r16")
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
    allocate(s1c(lnum, 1), s1dc(lnum, 1))
    allocate(is1e(lnum, 1), is1de(lnum, 1), naccs(lnum, 1))

    do i = 1, n_x
      if (.not. is_finite_r16(xvec(i)) .or. xvec(i) < 0.0_rk) then
        status(i) = -3_c_int
        cycle
      end if

      x = real(xvec(i), wk)

      ioprad = 1_c_int
      call oblfcn(c_w, m, lnum, ioprad, x, iopang, iopnorm, narg, arg, &
                  r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                  s1c, is1e, s1dc, is1de, naccs)
      r1v = r1c(idx0) * pow10_i(ir1e(idx0))
      r1d = r1dc(idx0) * pow10_i(ir1de(idx0))

      if (kind >= 2_c_int) then
        ioprad = 2_c_int
        call oblfcn(c_w, m, lnum, ioprad, x, iopang, iopnorm, narg, arg, &
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
  end subroutine oblate_rmn_batch_r16

end module oblate_batch_fortran
