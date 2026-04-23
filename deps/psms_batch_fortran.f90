module psms_batch_fortran
  use, intrinsic :: iso_c_binding
  use prolate_swf, only: profcn
  implicit none

  integer, parameter :: rk = c_double

contains

  pure real(rk) function pow10_i(exp10) result(v)
    integer(c_int), intent(in) :: exp10
    v = 10.0_rk ** real(exp10, rk)
  end function pow10_i

  pure logical function is_finite_r8(x) result(ok)
    real(rk), intent(in) :: x
    ok = (x == x) .and. (abs(x) < huge(x))
  end function is_finite_r8

  ! Inputs:
  !   m, n, c, eta(:), normalize
  ! Outputs:
  !   value(:), derivative(:), status
  ! status codes:
  !   0 success
  !  -1 invalid n
  !  -2 n < m
  !  -3 invalid eta domain
  !  -4 invalid n_eta
  subroutine psms_smn_batch_r8(m, n, c, n_eta, eta, normalize, value, derivative, status) bind(C, name="psms_smn_batch_r8")
    integer(c_int), value, intent(in) :: m, n, n_eta, normalize
    real(c_double), value, intent(in) :: c
    real(c_double), intent(in) :: eta(*)
    real(c_double), intent(out) :: value(*), derivative(*)
    integer(c_int), intent(out) :: status

    integer(c_int) :: i, idx0, lnum
    integer(c_int) :: ioprad, iopang, iopnorm, narg
    real(c_double) :: x1

    real(c_double), allocatable :: arg(:)
    real(c_double), allocatable :: r1c(:), r1dc(:), r2c(:), r2dc(:)
    real(c_double), allocatable :: s1c(:,:), s1dc(:,:)
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
      if (abs(eta(i)) > 1.0_rk .or. .not. is_finite_r8(eta(i))) then
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
    x1 = 1.0_rk

    allocate(arg(narg))
    arg(:) = eta(1:narg)

    allocate(r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum))
    allocate(ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), naccr(lnum))
    allocate(s1c(lnum, narg), s1dc(lnum, narg))
    allocate(is1e(lnum, narg), is1de(lnum, narg), naccs(lnum, narg))

    call profcn(c, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
                r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                s1c, is1e, s1dc, is1de, naccs)

    do i = 1, n_eta
      value(i) = s1c(idx0, i) * pow10_i(is1e(idx0, i))
      derivative(i) = s1dc(idx0, i) * pow10_i(is1de(idx0, i))
    end do
  end subroutine psms_smn_batch_r8

  ! Inputs:
  !   m, n, c, xi(:), kind
  ! Outputs:
  !   value_re(:), value_im(:), deriv_re(:), deriv_im(:), status(:)
  ! status(i) codes:
  !   0 success
  !  -1 invalid n
  !  -2 n < m
  !  -3 invalid xi domain (|xi| < 1)
  !  -4 invalid kind
  subroutine psms_rmn_batch_r8(m, n, c, n_xi, xi, kind, value_re, value_im, deriv_re, deriv_im, status) bind(C, name="psms_rmn_batch_r8")
    integer(c_int), value, intent(in) :: m, n, n_xi, kind
    real(c_double), value, intent(in) :: c
    real(c_double), intent(in) :: xi(*)
    real(c_double), intent(out) :: value_re(*), value_im(*), deriv_re(*), deriv_im(*)
    integer(c_int), intent(out) :: status(*)

    integer(c_int) :: i, idx0, lnum, ioprad, iopang, iopnorm, narg
    real(c_double) :: x1, x1_r

    real(c_double), allocatable :: arg(:)
    real(c_double), allocatable :: r1c(:), r1dc(:), r2c(:), r2dc(:)
    real(c_double), allocatable :: s1c(:,:), s1dc(:,:)
    integer(c_int), allocatable :: ir1e(:), ir1de(:), ir2e(:), ir2de(:), naccr(:)
    integer(c_int), allocatable :: is1e(:,:), is1de(:,:), naccs(:,:)

    real(c_double) :: r1v, r1d, r2v, r2d

    narg = 1_c_int
    iopang = 0_c_int
    iopnorm = 0_c_int

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
    arg(1) = 1.0_rk

    allocate(r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum))
    allocate(ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), naccr(lnum))
    allocate(s1c(lnum, 1), s1dc(lnum, 1))
    allocate(is1e(lnum, 1), is1de(lnum, 1), naccs(lnum, 1))

    do i = 1, n_xi
      if (abs(xi(i)) < 1.0_rk .or. .not. is_finite_r8(xi(i))) then
        status(i) = -3_c_int
        cycle
      end if

      x1_r = xi(i) - 1.0_rk
      x1 = x1_r

      ! First-kind radial values
      ioprad = 1_c_int
      call profcn(c, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
                  r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                  s1c, is1e, s1dc, is1de, naccs)

      r1v = r1c(idx0) * pow10_i(ir1e(idx0))
      r1d = r1dc(idx0) * pow10_i(ir1de(idx0))

      ! Second-kind radial values (needed for kinds 2,3,4)
      if (kind >= 2_c_int) then
        ioprad = 2_c_int
        call profcn(c, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
                    r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                    s1c, is1e, s1dc, is1de, naccs)

        r2v = r2c(idx0) * pow10_i(ir2e(idx0))
        r2d = r2dc(idx0) * pow10_i(ir2de(idx0))
      else
        r2v = 0.0_rk
        r2d = 0.0_rk
      end if

      select case (kind)
      case (1)
        value_re(i) = r1v
        value_im(i) = 0.0_rk
        deriv_re(i) = r1d
        deriv_im(i) = 0.0_rk
      case (2)
        value_re(i) = r2v
        value_im(i) = 0.0_rk
        deriv_re(i) = r2d
        deriv_im(i) = 0.0_rk
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

  end subroutine psms_rmn_batch_r8

  ! Angular function with accuracy output
  ! Inputs:
  !   m, n, c, eta(:), normalize
  ! Outputs:
  !   value(:), derivative(:), naccs(:), status
  subroutine psms_smn_batch_r8_acc(m, n, c, n_eta, eta, normalize, value, derivative, naccs, status) bind(C, name="psms_smn_batch_r8_acc")
    integer(c_int), value, intent(in) :: m, n, n_eta, normalize
    real(c_double), value, intent(in) :: c
    real(c_double), intent(in) :: eta(*)
    real(c_double), intent(out) :: value(*), derivative(*)
    integer(c_int), intent(out) :: naccs(*), status

    integer(c_int) :: i, idx0, lnum
    integer(c_int) :: ioprad, iopang, iopnorm, narg
    real(c_double) :: x1

    real(c_double), allocatable :: arg(:)
    real(c_double), allocatable :: r1c(:), r1dc(:), r2c(:), r2dc(:)
    real(c_double), allocatable :: s1c(:,:), s1dc(:,:)
    integer(c_int), allocatable :: ir1e(:), ir1de(:), ir2e(:), ir2de(:), naccr(:)
    integer(c_int), allocatable :: is1e(:,:), is1de(:,:), naccs_tmp(:,:)

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
      if (abs(eta(i)) > 1.0_rk .or. .not. is_finite_r8(eta(i))) then
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
    x1 = 1.0_rk

    allocate(arg(narg))
    arg(:) = eta(1:narg)

    allocate(r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum))
    allocate(ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), naccr(lnum))
    allocate(s1c(lnum, narg), s1dc(lnum, narg))
    allocate(is1e(lnum, narg), is1de(lnum, narg), naccs_tmp(lnum, narg))

    call profcn(c, m, lnum, ioprad, x1, iopang, iopnorm, narg, arg, &
                r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                s1c, is1e, s1dc, is1de, naccs_tmp)

    do i = 1, n_eta
      value(i) = s1c(idx0, i) * pow10_i(is1e(idx0, i))
      derivative(i) = s1dc(idx0, i) * pow10_i(is1de(idx0, i))
      naccs(i) = naccs_tmp(idx0, i)
    end do
  end subroutine psms_smn_batch_r8_acc

  ! Radial function with accuracy output
  ! Inputs:
  !   m, n, c, xi(:), kind
  ! Outputs:
  !   value_re(:), value_im(:), deriv_re(:), deriv_im(:), naccr(:), status(:)
  subroutine psms_rmn_batch_r8_acc(m, n, c, n_xi, xi, kind, value_re, value_im, deriv_re, deriv_im, naccr_out, status) bind(C, name="psms_rmn_batch_r8_acc")
    integer(c_int), value, intent(in) :: m, n, n_xi, kind
    real(c_double), value, intent(in) :: c
    real(c_double), intent(in) :: xi(*)
    real(c_double), intent(out) :: value_re(*), value_im(*), deriv_re(*), deriv_im(*)
    integer(c_int), intent(out) :: naccr_out(*), status(*)

    integer(c_int) :: i, idx0, lnum, ioprad, iopang, iopnorm, narg
    real(c_double) :: x1, x1_r

    real(c_double), allocatable :: arg(:)
    real(c_double), allocatable :: r1c(:), r1dc(:), r2c(:), r2dc(:)
    real(c_double), allocatable :: s1c(:,:), s1dc(:,:)
    integer(c_int), allocatable :: ir1e(:), ir1de(:), ir2e(:), ir2de(:), naccr(:)
    integer(c_int), allocatable :: is1e(:,:), is1de(:,:), naccs(:,:)

    real(c_double) :: r1v, r1d, r2v, r2d

    narg = 1_c_int
    iopang = 0_c_int
    iopnorm = 0_c_int

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

    allocate(arg(narg))
    allocate(r1c(lnum), r1dc(lnum), r2c(lnum), r2dc(lnum))
    allocate(ir1e(lnum), ir1de(lnum), ir2e(lnum), ir2de(lnum), naccr(lnum))
    allocate(s1c(lnum, narg), s1dc(lnum, narg))
    allocate(is1e(lnum, narg), is1de(lnum, narg), naccs(lnum, narg))

    do i = 1, n_xi
      x1 = xi(i) - 1.0_rk
      if (abs(xi(i)) < 1.0_rk) then
        status(i) = -3_c_int
        cycle
      end if

      arg(1) = 0.0_rk
      ioprad = 2_c_int
      x1_r = x1

      call profcn(c, m, lnum, ioprad, x1_r, iopang, iopnorm, narg, arg, &
                  r1c, ir1e, r1dc, ir1de, r2c, ir2e, r2dc, ir2de, naccr, &
                  s1c, is1e, s1dc, is1de, naccs)

      r1v = r1c(idx0) * pow10_i(ir1e(idx0))
      r1d = r1dc(idx0) * pow10_i(ir1de(idx0))
      r2v = r2c(idx0) * pow10_i(ir2e(idx0))
      r2d = r2dc(idx0) * pow10_i(ir2de(idx0))

      naccr_out(i) = naccr(idx0)

      select case (kind)
      case (1)
        value_re(i) = r1v
        value_im(i) = 0.0_rk
        deriv_re(i) = r1d
        deriv_im(i) = 0.0_rk
      case (2)
        value_re(i) = r2v
        value_im(i) = 0.0_rk
        deriv_re(i) = r2d
        deriv_im(i) = 0.0_rk
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

  end subroutine psms_rmn_batch_r8_acc

end module psms_batch_fortran
