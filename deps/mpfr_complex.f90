module spheroidal_mpfr_complex
  use, intrinsic :: iso_c_binding
  use spheroidal_mpfr_real
  implicit none

  type :: mp_complex
    type(mp_real) :: re, im
  end type

  interface
    subroutine sw_mpc_binary(re, im, ar, ai, br, bi, operation, rounding) bind(C)
      import c_int64_t, c_int
      integer(c_int64_t), intent(inout) :: re(*), im(*)
      integer(c_int64_t), intent(in) :: ar(*), ai(*), br(*), bi(*)
      integer(c_int), value :: operation, rounding
    end subroutine
    subroutine sw_mpc_unary(re, im, ar, ai, operation, rounding) bind(C)
      import c_int64_t, c_int
      integer(c_int64_t), intent(inout) :: re(*), im(*)
      integer(c_int64_t), intent(in) :: ar(*), ai(*)
      integer(c_int), value :: operation, rounding
    end subroutine
    subroutine sw_mpc_power_si(re, im, ar, ai, exponent, rounding) bind(C)
      import c_int64_t, c_int, c_long
      integer(c_int64_t), intent(inout) :: re(*), im(*)
      integer(c_int64_t), intent(in) :: ar(*), ai(*)
      integer(c_long), value :: exponent
      integer(c_int), value :: rounding
    end subroutine
    subroutine sw_mpc_abs(result, ar, ai, rounding) bind(C)
      import c_int64_t, c_int
      integer(c_int64_t), intent(inout) :: result(*)
      integer(c_int64_t), intent(in) :: ar(*), ai(*)
      integer(c_int), value :: rounding
    end subroutine
  end interface
  interface mpcmplx
    module procedure mpc_from_c, mpc_from_r, mpc_from_i, mpc_from_j, mpc_from_f, mpc_from_d, mpc_from_z, mpc_from_w, mpc_from_ir, mpc_from_ri
  end interface
  interface assignment(=)
    module procedure mpc_assign_r, mpc_assign_i, mpc_assign_j, mpc_assign_f, mpc_assign_d, mpc_assign_z, mpc_assign_w
    module procedure mpc_assign_rc
    module procedure mpc_assign_ic, mpc_assign_jc, mpc_assign_fc, mpc_assign_dc
  end interface
  interface operator(+)
    module procedure mpc_add_cc, mpc_add_cr, mpc_add_rc, mpc_add_ci, mpc_add_ic, mpc_add_cj, mpc_add_jc, mpc_add_cf, mpc_add_fc, mpc_add_cd, mpc_add_dc, mpc_positive
  end interface
  interface operator(-)
    module procedure mpc_sub_cc, mpc_sub_cr, mpc_sub_rc, mpc_sub_ci, mpc_sub_ic, mpc_sub_cj, mpc_sub_jc, mpc_sub_cf, mpc_sub_fc, mpc_sub_cd, mpc_sub_dc, mpc_negative
  end interface
  interface operator(*)
    module procedure mpc_mul_cc, mpc_mul_cr, mpc_mul_rc, mpc_mul_ci, mpc_mul_ic, mpc_mul_cj, mpc_mul_jc, mpc_mul_cf, mpc_mul_fc, mpc_mul_cd, mpc_mul_dc
  end interface
  interface operator(/)
    module procedure mpc_div_cc, mpc_div_cr, mpc_div_rc, mpc_div_ci, mpc_div_ic, mpc_div_cj, mpc_div_jc, mpc_div_cf, mpc_div_fc, mpc_div_cd, mpc_div_dc
  end interface
  interface operator(**)
    module procedure mpc_pow_cc, mpc_pow_cr, mpc_pow_rc, mpc_pow_ci, mpc_pow_ic, mpc_pow_cj, mpc_pow_jc, mpc_pow_cf, mpc_pow_fc, mpc_pow_cd, mpc_pow_dc
  end interface
  interface operator(==)
    module procedure mpc_eq_cc, mpc_eq_cr, mpc_eq_rc, mpc_eq_ci, mpc_eq_ic, mpc_eq_cj, mpc_eq_jc, mpc_eq_cf, mpc_eq_fc, mpc_eq_cd, mpc_eq_dc
  end interface
  interface operator(/=)
    module procedure mpc_ne_cc, mpc_ne_cr, mpc_ne_rc, mpc_ne_ci, mpc_ne_ic, mpc_ne_cj, mpc_ne_jc, mpc_ne_cf, mpc_ne_fc, mpc_ne_cd, mpc_ne_dc
  end interface
  interface real
    module procedure mpc_real
  end interface
  interface mpreal
    module procedure mpc_real
  end interface
  interface aimag
    module procedure mpc_aimag
  end interface
  interface conjg
    module procedure mpc_conjg
  end interface
  interface abs
    module procedure mpc_abs
  end interface
  interface sqrt
    module procedure mpc_sqrt
  end interface
  interface exp
    module procedure mpc_exp
  end interface
  interface log
    module procedure mpc_log
  end interface
  interface sin
    module procedure mpc_sin
  end interface
  interface cos
    module procedure mpc_cos
  end interface
  interface sinh
    module procedure mpc_sinh
  end interface
  interface cosh
    module procedure mpc_cosh
  end interface
  interface precision
    module procedure mpc_precision
  end interface
  interface range
    module procedure mpc_range
  end interface
  interface int
    module procedure mpc_int
  end interface
contains
  impure elemental subroutine mpc_assign_ic(z, x)
    integer, intent(out) :: z
    type(mp_complex), intent(in) :: x
    z = int(x%re)
  end subroutine

  impure elemental subroutine mpc_assign_jc(z, x)
    integer(c_int64_t), intent(out) :: z
    type(mp_complex), intent(in) :: x
    z = mp_to_integer64(x%re)
  end subroutine

  impure elemental subroutine mpc_assign_fc(z, x)
    real(c_float), intent(out) :: z
    type(mp_complex), intent(in) :: x
    z = real(mp_to_double(x%re),kind=c_float)
  end subroutine

  impure elemental subroutine mpc_assign_dc(z, x)
    real(c_double), intent(out) :: z
    type(mp_complex), intent(in) :: x
    z = mp_to_double(x%re)
  end subroutine

  impure elemental function mpc_int(a, kind) result(x)
    type(mp_complex), intent(in) :: a
    integer, optional, intent(in) :: kind
    integer :: x
    x = int(a%re)
  end function

  impure elemental subroutine mpc_assign_rc(z, x)
    type(mp_real), intent(out) :: z
    type(mp_complex), intent(in) :: x
    z = mpreal(x%re)
  end subroutine

  impure elemental function mpc_from_c(x, kind) result(z)
    type(mp_complex), intent(in) :: x
    integer, optional, intent(in) :: kind
    type(mp_complex) :: z
    z%re = x%re
    z%im = x%im
  end function
  impure elemental function mpc_from_r(x, y, kind) result(z)
    type(mp_real), intent(in) :: x
    type(mp_real), optional, intent(in) :: y
    integer, optional, intent(in) :: kind
    type(mp_complex) :: z
    z%re = mpreal(x)
    z%im = mpreal(0)
    if (present(y)) z%im = mpreal(y)
  end function
  impure elemental function mpc_from_i(x, y, kind) result(z)
    integer, intent(in) :: x
    integer, optional, intent(in) :: y
    integer, optional, intent(in) :: kind
    type(mp_complex) :: z
    z%re = mpreal(x)
    z%im = mpreal(0)
    if (present(y)) z%im = mpreal(y)
  end function
  impure elemental function mpc_from_j(x, y, kind) result(z)
    integer(c_int64_t), intent(in) :: x
    integer(c_int64_t), optional, intent(in) :: y
    integer, optional, intent(in) :: kind
    type(mp_complex) :: z
    z%re = mpreal(x)
    z%im = mpreal(0)
    if (present(y)) z%im = mpreal(y)
  end function
  impure elemental function mpc_from_f(x, y, kind) result(z)
    real(c_float), intent(in) :: x
    real(c_float), optional, intent(in) :: y
    integer, optional, intent(in) :: kind
    type(mp_complex) :: z
    z%re = mpreal(x)
    z%im = mpreal(0)
    if (present(y)) z%im = mpreal(y)
  end function
  impure elemental function mpc_from_d(x, y, kind) result(z)
    real(c_double), intent(in) :: x
    real(c_double), optional, intent(in) :: y
    integer, optional, intent(in) :: kind
    type(mp_complex) :: z
    z%re = mpreal(x)
    z%im = mpreal(0)
    if (present(y)) z%im = mpreal(y)
  end function
  impure elemental function mpc_from_ir(x, y, kind) result(z)
    integer, intent(in) :: x
    type(mp_real), intent(in) :: y
    integer, optional, intent(in) :: kind
    type(mp_complex) :: z
    z%re = mpreal(x)
    z%im = mpreal(y)
  end function
  impure elemental function mpc_from_ri(x, y, kind) result(z)
    type(mp_real), intent(in) :: x
    integer, intent(in) :: y
    integer, optional, intent(in) :: kind
    type(mp_complex) :: z
    z%re = mpreal(x)
    z%im = mpreal(y)
  end function
  impure elemental function mpc_from_z(x, kind) result(z)
    complex(c_float), intent(in) :: x
    integer, optional, intent(in) :: kind
    type(mp_complex) :: z
    z%re = mpreal(real(x,kind=c_double))
    z%im = mpreal(real(aimag(x),kind=c_double))
  end function
  impure elemental function mpc_from_w(x, kind) result(z)
    complex(c_double), intent(in) :: x
    integer, optional, intent(in) :: kind
    type(mp_complex) :: z
    z%re = mpreal(real(x,kind=c_double))
    z%im = mpreal(real(aimag(x),kind=c_double))
  end function
  impure elemental subroutine mpc_assign_r(z, x)
    type(mp_complex), intent(out) :: z
    type(mp_real), intent(in) :: x
    z = mpcmplx(x)
  end subroutine
  impure elemental subroutine mpc_assign_i(z, x)
    type(mp_complex), intent(out) :: z
    integer, intent(in) :: x
    z = mpcmplx(x)
  end subroutine
  impure elemental subroutine mpc_assign_j(z, x)
    type(mp_complex), intent(out) :: z
    integer(c_int64_t), intent(in) :: x
    z = mpcmplx(x)
  end subroutine
  impure elemental subroutine mpc_assign_f(z, x)
    type(mp_complex), intent(out) :: z
    real(c_float), intent(in) :: x
    z = mpcmplx(x)
  end subroutine
  impure elemental subroutine mpc_assign_d(z, x)
    type(mp_complex), intent(out) :: z
    real(c_double), intent(in) :: x
    z = mpcmplx(x)
  end subroutine
  impure elemental subroutine mpc_assign_z(z, x)
    type(mp_complex), intent(out) :: z
    complex(c_float), intent(in) :: x
    z = mpcmplx(x)
  end subroutine
  impure elemental subroutine mpc_assign_w(z, x)
    type(mp_complex), intent(out) :: z
    complex(c_double), intent(in) :: x
    z = mpcmplx(x)
  end subroutine
  impure elemental function mpc_binary(a, b, operation) result(z)
    type(mp_complex), intent(in) :: a, b
    integer, intent(in) :: operation
    type(mp_complex) :: z
    call mp_init(z%re)
    call mp_init(z%im)
    call sw_mpc_binary(z%re%data, z%im%data, a%re%data, a%im%data, b%re%data, b%im%data, int(operation,c_int), int(mpfr_rounding,c_int))
  end function

  impure elemental function mpc_unary(a, operation) result(z)
    type(mp_complex), intent(in) :: a
    integer, intent(in) :: operation
    type(mp_complex) :: z
    call mp_init(z%re)
    call mp_init(z%im)
    call sw_mpc_unary(z%re%data, z%im%data, a%re%data, a%im%data, int(operation,c_int), int(mpfr_rounding,c_int))
  end function
  impure elemental function mpc_add_cc(a, b) result(z)
    type(mp_complex), intent(in) :: a, b
    type(mp_complex) :: z
    z = mpc_binary(a, b, 0)
  end function
  impure elemental function mpc_add_cr(a, b) result(z)
    type(mp_complex), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_add_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_add_rc(a, b) result(z)
    type(mp_real), intent(in) :: a
    type(mp_complex), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_add_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_add_ci(a, b) result(z)
    type(mp_complex), intent(in) :: a
    integer, intent(in) :: b
    type(mp_complex) :: z
    z = mpc_add_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_add_ic(a, b) result(z)
    integer, intent(in) :: a
    type(mp_complex), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_add_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_add_cj(a, b) result(z)
    type(mp_complex), intent(in) :: a
    integer(c_int64_t), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_add_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_add_jc(a, b) result(z)
    integer(c_int64_t), intent(in) :: a
    type(mp_complex), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_add_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_add_cf(a, b) result(z)
    type(mp_complex), intent(in) :: a
    real(c_float), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_add_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_add_fc(a, b) result(z)
    real(c_float), intent(in) :: a
    type(mp_complex), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_add_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_add_cd(a, b) result(z)
    type(mp_complex), intent(in) :: a
    real(c_double), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_add_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_add_dc(a, b) result(z)
    real(c_double), intent(in) :: a
    type(mp_complex), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_add_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_sub_cc(a, b) result(z)
    type(mp_complex), intent(in) :: a, b
    type(mp_complex) :: z
    z = mpc_binary(a, b, 1)
  end function
  impure elemental function mpc_sub_cr(a, b) result(z)
    type(mp_complex), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_sub_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_sub_rc(a, b) result(z)
    type(mp_real), intent(in) :: a
    type(mp_complex), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_sub_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_sub_ci(a, b) result(z)
    type(mp_complex), intent(in) :: a
    integer, intent(in) :: b
    type(mp_complex) :: z
    z = mpc_sub_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_sub_ic(a, b) result(z)
    integer, intent(in) :: a
    type(mp_complex), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_sub_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_sub_cj(a, b) result(z)
    type(mp_complex), intent(in) :: a
    integer(c_int64_t), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_sub_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_sub_jc(a, b) result(z)
    integer(c_int64_t), intent(in) :: a
    type(mp_complex), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_sub_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_sub_cf(a, b) result(z)
    type(mp_complex), intent(in) :: a
    real(c_float), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_sub_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_sub_fc(a, b) result(z)
    real(c_float), intent(in) :: a
    type(mp_complex), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_sub_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_sub_cd(a, b) result(z)
    type(mp_complex), intent(in) :: a
    real(c_double), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_sub_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_sub_dc(a, b) result(z)
    real(c_double), intent(in) :: a
    type(mp_complex), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_sub_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_mul_cc(a, b) result(z)
    type(mp_complex), intent(in) :: a, b
    type(mp_complex) :: z
    z = mpc_binary(a, b, 2)
  end function
  impure elemental function mpc_mul_cr(a, b) result(z)
    type(mp_complex), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_mul_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_mul_rc(a, b) result(z)
    type(mp_real), intent(in) :: a
    type(mp_complex), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_mul_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_mul_ci(a, b) result(z)
    type(mp_complex), intent(in) :: a
    integer, intent(in) :: b
    type(mp_complex) :: z
    z = mpc_mul_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_mul_ic(a, b) result(z)
    integer, intent(in) :: a
    type(mp_complex), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_mul_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_mul_cj(a, b) result(z)
    type(mp_complex), intent(in) :: a
    integer(c_int64_t), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_mul_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_mul_jc(a, b) result(z)
    integer(c_int64_t), intent(in) :: a
    type(mp_complex), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_mul_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_mul_cf(a, b) result(z)
    type(mp_complex), intent(in) :: a
    real(c_float), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_mul_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_mul_fc(a, b) result(z)
    real(c_float), intent(in) :: a
    type(mp_complex), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_mul_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_mul_cd(a, b) result(z)
    type(mp_complex), intent(in) :: a
    real(c_double), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_mul_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_mul_dc(a, b) result(z)
    real(c_double), intent(in) :: a
    type(mp_complex), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_mul_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_div_cc(a, b) result(z)
    type(mp_complex), intent(in) :: a, b
    type(mp_complex) :: z
    z = mpc_binary(a, b, 3)
  end function
  impure elemental function mpc_div_cr(a, b) result(z)
    type(mp_complex), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_div_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_div_rc(a, b) result(z)
    type(mp_real), intent(in) :: a
    type(mp_complex), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_div_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_div_ci(a, b) result(z)
    type(mp_complex), intent(in) :: a
    integer, intent(in) :: b
    type(mp_complex) :: z
    z = mpc_div_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_div_ic(a, b) result(z)
    integer, intent(in) :: a
    type(mp_complex), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_div_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_div_cj(a, b) result(z)
    type(mp_complex), intent(in) :: a
    integer(c_int64_t), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_div_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_div_jc(a, b) result(z)
    integer(c_int64_t), intent(in) :: a
    type(mp_complex), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_div_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_div_cf(a, b) result(z)
    type(mp_complex), intent(in) :: a
    real(c_float), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_div_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_div_fc(a, b) result(z)
    real(c_float), intent(in) :: a
    type(mp_complex), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_div_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_div_cd(a, b) result(z)
    type(mp_complex), intent(in) :: a
    real(c_double), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_div_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_div_dc(a, b) result(z)
    real(c_double), intent(in) :: a
    type(mp_complex), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_div_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_pow_cc(a, b) result(z)
    type(mp_complex), intent(in) :: a, b
    type(mp_complex) :: z
    z = mpc_binary(a, b, 4)
  end function
  impure elemental function mpc_pow_cr(a, b) result(z)
    type(mp_complex), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_pow_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_pow_rc(a, b) result(z)
    type(mp_real), intent(in) :: a
    type(mp_complex), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_pow_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_pow_ci(a, b) result(z)
    type(mp_complex), intent(in) :: a
    integer, intent(in) :: b
    type(mp_complex) :: z
    call mp_init(z%re)
    call mp_init(z%im)
    call sw_mpc_power_si(z%re%data, z%im%data, a%re%data, a%im%data, int(b,c_long), int(mpfr_rounding,c_int))
  end function
  impure elemental function mpc_pow_ic(a, b) result(z)
    integer, intent(in) :: a
    type(mp_complex), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_pow_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_pow_cj(a, b) result(z)
    type(mp_complex), intent(in) :: a
    integer(c_int64_t), intent(in) :: b
    type(mp_complex) :: z
    call mp_init(z%re)
    call mp_init(z%im)
    call sw_mpc_power_si(z%re%data, z%im%data, a%re%data, a%im%data, int(b,c_long), int(mpfr_rounding,c_int))
  end function
  impure elemental function mpc_pow_jc(a, b) result(z)
    integer(c_int64_t), intent(in) :: a
    type(mp_complex), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_pow_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_pow_cf(a, b) result(z)
    type(mp_complex), intent(in) :: a
    real(c_float), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_pow_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_pow_fc(a, b) result(z)
    real(c_float), intent(in) :: a
    type(mp_complex), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_pow_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_pow_cd(a, b) result(z)
    type(mp_complex), intent(in) :: a
    real(c_double), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_pow_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_pow_dc(a, b) result(z)
    real(c_double), intent(in) :: a
    type(mp_complex), intent(in) :: b
    type(mp_complex) :: z
    z = mpc_pow_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_eq_cc(a, b) result(z)
    type(mp_complex), intent(in) :: a, b
    logical :: z
    z = (a%re == b%re .and. a%im == b%im)
  end function
  impure elemental function mpc_eq_cr(a, b) result(z)
    type(mp_complex), intent(in) :: a
    type(mp_real), intent(in) :: b
    logical :: z
    z = mpc_eq_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_eq_rc(a, b) result(z)
    type(mp_real), intent(in) :: a
    type(mp_complex), intent(in) :: b
    logical :: z
    z = mpc_eq_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_eq_ci(a, b) result(z)
    type(mp_complex), intent(in) :: a
    integer, intent(in) :: b
    logical :: z
    z = mpc_eq_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_eq_ic(a, b) result(z)
    integer, intent(in) :: a
    type(mp_complex), intent(in) :: b
    logical :: z
    z = mpc_eq_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_eq_cj(a, b) result(z)
    type(mp_complex), intent(in) :: a
    integer(c_int64_t), intent(in) :: b
    logical :: z
    z = mpc_eq_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_eq_jc(a, b) result(z)
    integer(c_int64_t), intent(in) :: a
    type(mp_complex), intent(in) :: b
    logical :: z
    z = mpc_eq_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_eq_cf(a, b) result(z)
    type(mp_complex), intent(in) :: a
    real(c_float), intent(in) :: b
    logical :: z
    z = mpc_eq_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_eq_fc(a, b) result(z)
    real(c_float), intent(in) :: a
    type(mp_complex), intent(in) :: b
    logical :: z
    z = mpc_eq_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_eq_cd(a, b) result(z)
    type(mp_complex), intent(in) :: a
    real(c_double), intent(in) :: b
    logical :: z
    z = mpc_eq_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_eq_dc(a, b) result(z)
    real(c_double), intent(in) :: a
    type(mp_complex), intent(in) :: b
    logical :: z
    z = mpc_eq_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_ne_cc(a, b) result(z)
    type(mp_complex), intent(in) :: a, b
    logical :: z
    z = .not. (a%re == b%re .and. a%im == b%im)
  end function
  impure elemental function mpc_ne_cr(a, b) result(z)
    type(mp_complex), intent(in) :: a
    type(mp_real), intent(in) :: b
    logical :: z
    z = mpc_ne_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_ne_rc(a, b) result(z)
    type(mp_real), intent(in) :: a
    type(mp_complex), intent(in) :: b
    logical :: z
    z = mpc_ne_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_ne_ci(a, b) result(z)
    type(mp_complex), intent(in) :: a
    integer, intent(in) :: b
    logical :: z
    z = mpc_ne_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_ne_ic(a, b) result(z)
    integer, intent(in) :: a
    type(mp_complex), intent(in) :: b
    logical :: z
    z = mpc_ne_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_ne_cj(a, b) result(z)
    type(mp_complex), intent(in) :: a
    integer(c_int64_t), intent(in) :: b
    logical :: z
    z = mpc_ne_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_ne_jc(a, b) result(z)
    integer(c_int64_t), intent(in) :: a
    type(mp_complex), intent(in) :: b
    logical :: z
    z = mpc_ne_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_ne_cf(a, b) result(z)
    type(mp_complex), intent(in) :: a
    real(c_float), intent(in) :: b
    logical :: z
    z = mpc_ne_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_ne_fc(a, b) result(z)
    real(c_float), intent(in) :: a
    type(mp_complex), intent(in) :: b
    logical :: z
    z = mpc_ne_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_ne_cd(a, b) result(z)
    type(mp_complex), intent(in) :: a
    real(c_double), intent(in) :: b
    logical :: z
    z = mpc_ne_cc(a, mpcmplx(b))
  end function
  impure elemental function mpc_ne_dc(a, b) result(z)
    real(c_double), intent(in) :: a
    type(mp_complex), intent(in) :: b
    logical :: z
    z = mpc_ne_cc(mpcmplx(a), b)
  end function
  impure elemental function mpc_positive(a) result(z)
    type(mp_complex), intent(in) :: a
    type(mp_complex) :: z
    z = a
  end function
  impure elemental function mpc_negative(a) result(z)
    type(mp_complex), intent(in) :: a
    type(mp_complex) :: z
    z = mpc_unary(a, 0)
  end function
  impure elemental function mpc_conjg(a) result(z)
    type(mp_complex), intent(in) :: a
    type(mp_complex) :: z
    z = mpc_unary(a, 1)
  end function
  impure elemental function mpc_sqrt(a) result(z)
    type(mp_complex), intent(in) :: a
    type(mp_complex) :: z
    z = mpc_unary(a, 2)
  end function
  impure elemental function mpc_exp(a) result(z)
    type(mp_complex), intent(in) :: a
    type(mp_complex) :: z
    z = mpc_unary(a, 3)
  end function
  impure elemental function mpc_log(a) result(z)
    type(mp_complex), intent(in) :: a
    type(mp_complex) :: z
    z = mpc_unary(a, 4)
  end function
  impure elemental function mpc_sin(a) result(z)
    type(mp_complex), intent(in) :: a
    type(mp_complex) :: z
    z = mpc_unary(a, 5)
  end function
  impure elemental function mpc_cos(a) result(z)
    type(mp_complex), intent(in) :: a
    type(mp_complex) :: z
    z = mpc_unary(a, 6)
  end function
  impure elemental function mpc_sinh(a) result(z)
    type(mp_complex), intent(in) :: a
    type(mp_complex) :: z
    z = mpc_unary(a, 7)
  end function
  impure elemental function mpc_cosh(a) result(z)
    type(mp_complex), intent(in) :: a
    type(mp_complex) :: z
    z = mpc_unary(a, 8)
  end function
  impure elemental function mpc_real(a, kind) result(x)
    type(mp_complex), intent(in) :: a
    integer, optional, intent(in) :: kind
    type(mp_real) :: x
    x = mpreal(a%re)
  end function

  impure elemental function mpc_aimag(a) result(x)
    type(mp_complex), intent(in) :: a
    type(mp_real) :: x
    x = mpreal(a%im)
  end function

  impure elemental function mpc_abs(a) result(x)
    type(mp_complex), intent(in) :: a
    type(mp_real) :: x
    call mp_init(x)
    call sw_mpc_abs(x%data, a%re%data, a%im%data, int(mpfr_rounding,c_int))
  end function

  impure elemental function mpc_precision(a) result(bits)
    type(mp_complex), intent(in) :: a
    integer :: bits
    bits = precision(a%re)
  end function

  impure elemental function mpc_range(a) result(exponent)
    type(mp_complex), intent(in) :: a
    integer :: exponent
    exponent = range(a%re)
  end function
end module spheroidal_mpfr_complex
