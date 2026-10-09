! Exercise cache misses, reuse, and replacement on independent OpenMP threads.
! Guard cells catch a node-count/tolerance swap even for small quadrature rules.
program native_cache
  use oblate_swf, only: gauss_cached
  use oblate_parameters, only: knd
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer :: iteration, failures

  failures = 0
!$omp parallel do reduction(+:failures)
  do iteration = 1, 256
    failures = failures + check_rule(iteration)
  end do
!$omp end parallel do
  if (failures /= 0) then
    print *, 'Quadrature cache failures:', failures
    stop 1
  end if
  print *, 'Passed 256 quadrature cache checks, including reuse and guard cells'

contains
  integer function check_rule(iteration) result(errors)
    integer, intent(in) :: iteration
    integer, parameter :: counts(4) = [2, 8, 32, 96]
    real(knd), parameter :: guard = -12345.0_knd, tolerance = 2.0e-13_knd
    real(knd) :: x(0:128), w(0:128), saved_x(96), saved_w(96), moment
    integer :: n, degree, repeat

    n = counts(mod(iteration-1, 4)+1)
    errors = 0
    do repeat = 1, 2
      x = guard
      w = guard
      call gauss_cached(15, n, x(1:n), w(1:n))
      if (x(0) /= guard .or. w(0) /= guard) errors = errors + 1
      if (any(x(n+1:) /= guard) .or. any(w(n+1:) /= guard)) errors = errors + 1
      if (.not. all(ieee_is_finite(x(1:n))) .or. .not. all(ieee_is_finite(w(1:n)))) errors = errors + 1
      if (any(abs(x(1:n)) >= 1) .or. any(w(1:n) <= 0)) errors = errors + 1
      do degree = 0, min(2*n-1, 10)
        moment = 0
        if (mod(degree,2) == 0) moment = 2.0_knd/(degree+1)
        if (abs(sum(w(1:n)*x(1:n)**degree)-moment) > tolerance) errors = errors + 1
      end do
      if (repeat == 1) then
        saved_x(1:n) = x(1:n)
        saved_w(1:n) = w(1:n)
      else
        if (any(x(1:n) /= saved_x(1:n)) .or. any(w(1:n) /= saved_w(1:n))) errors = errors + 1
      end if
    end do
  end function check_rule
end program native_cache
