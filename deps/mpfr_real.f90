module spheroidal_mpfr_real
  use, intrinsic :: iso_c_binding
  implicit none
  private
  type, public :: mp_real
    integer(c_int64_t), allocatable :: data(:)
  end type
  integer(c_int), public, save :: mpfr_bits = 256, mpfr_rounding = 0
  integer(c_int), public, save :: mpfr_status = 0
  integer(c_int64_t), public, save :: mpfr_workspace_limit_bytes = 268435456_c_int64_t
!$omp threadprivate(mpfr_bits, mpfr_rounding, mpfr_status, mpfr_workspace_limit_bytes)
  public :: mpreal, mp_init, mp_to_string, mp_to_double, mp_to_integer64, mp_is_finite
  public :: assignment(=), operator(+), operator(-), operator(*), operator(/), operator(**)
  public :: operator(==), operator(/=), operator(<), operator(<=), operator(>), operator(>=)
  public :: abs, sqrt, exp, log, log10, sin, cos, sinh, cosh, atan, asin, acos, tan, tanh
  public :: int, nint, real, dble, min, max, sign, mod, atan2, precision, range, epsilon, huge, tiny, sum, maxval, minval
  interface
    function c_words(bits) bind(C,name='sw_mpfr_allocate_words') result(words)
      import
      integer(c_int), value :: bits
      integer(c_int64_t) :: words
    end function
    subroutine c_init(out,bits) bind(C,name='sw_mpfr_init')
      import
      integer(c_int64_t), intent(out) :: out(*)
      integer(c_int), value :: bits
    end subroutine
    subroutine c_integer(out,value,rounding) bind(C,name='sw_mpfr_set_integer')
      import
      integer(c_int64_t), intent(inout) :: out(*)
      integer(c_int64_t), value :: value
      integer(c_int), value :: rounding
    end subroutine
    subroutine c_set_double(out,value,rounding) bind(C,name='sw_mpfr_set_double')
      import
      integer(c_int64_t), intent(inout) :: out(*)
      real(c_double), value :: value
      integer(c_int), value :: rounding
    end subroutine
    function c_string(out,value,rounding) bind(C,name='sw_mpfr_set_string') result(status)
      import
      integer(c_int64_t), intent(inout) :: out(*)
      character(c_char), intent(in) :: value(*)
      integer(c_int), value :: rounding
      integer(c_int) :: status
    end function
    subroutine c_binary(operation,out,a,b,rounding) bind(C,name='sw_mpfr_binary')
      import
      integer(c_int), value :: operation,rounding
      integer(c_int64_t), intent(inout) :: out(*)
      integer(c_int64_t), intent(in) :: a(*),b(*)
    end subroutine
    subroutine c_unary(operation,out,a,rounding) bind(C,name='sw_mpfr_unary')
      import
      integer(c_int), value :: operation,rounding
      integer(c_int64_t), intent(inout) :: out(*)
      integer(c_int64_t), intent(in) :: a(*)
    end subroutine
    subroutine c_power(out,a,exponent,rounding) bind(C,name='sw_mpfr_power_integer')
      import
      integer(c_int64_t), intent(inout) :: out(*)
      integer(c_int64_t), intent(in) :: a(*)
      integer(c_int64_t), value :: exponent
      integer(c_int), value :: rounding
    end subroutine
    function c_compare(a,b) bind(C,name='sw_mpfr_compare') result(order)
      import
      integer(c_int64_t), intent(in) :: a(*),b(*)
      integer(c_int) :: order
    end function
    function c_is_finite(a) bind(C,name='sw_mpfr_is_finite') result(finite)
      import
      integer(c_int64_t), intent(in) :: a(*)
      integer(c_int) :: finite
    end function
    function c_to_integer(a,rounding) bind(C,name='sw_mpfr_to_integer') result(value)
      import
      integer(c_int64_t), intent(in) :: a(*)
      integer(c_int), value :: rounding
      integer(c_int64_t) :: value
    end function
    function c_to_double(a,rounding) bind(C,name='sw_mpfr_to_double') result(value)
      import
      integer(c_int64_t), intent(in) :: a(*)
      integer(c_int), value :: rounding
      real(c_double) :: value
    end function
    subroutine c_constant(operation,out,rounding) bind(C,name='sw_mpfr_constant')
      import
      integer(c_int), value :: operation,rounding
      integer(c_int64_t), intent(inout) :: out(*)
    end subroutine
    function c_format(a,text,capacity) bind(C,name='sw_mpfr_format') result(length)
      import
      integer(c_int64_t), intent(in) :: a(*)
      character(c_char), intent(out) :: text(*)
      integer(c_int), value :: capacity
      integer(c_int) :: length
    end function
  end interface

  interface mpreal
    module procedure mp_from_m
    module procedure mp_from_i
    module procedure mp_from_d
    module procedure mp_from_s
    module procedure mp_from_l
    module procedure mp_from_c
  end interface
  interface assignment(=)
    module procedure mp_assign_i
    module procedure mp_assign_d
    module procedure mp_assign_s
    module procedure mp_assign_l
    module procedure mp_assign_c
    module procedure mp_assign_to_i
    module procedure mp_assign_to_d
    module procedure mp_assign_to_s
    module procedure mp_assign_to_l
  end interface
  interface operator(+)
    module procedure mp_positive
    module procedure mp_add_mm
    module procedure mp_add_mi
    module procedure mp_add_im
    module procedure mp_add_md
    module procedure mp_add_dm
    module procedure mp_add_ms
    module procedure mp_add_sm
    module procedure mp_add_ml
    module procedure mp_add_lm
  end interface
  interface operator(-)
    module procedure mp_negative
    module procedure mp_subtract_mm
    module procedure mp_subtract_mi
    module procedure mp_subtract_im
    module procedure mp_subtract_md
    module procedure mp_subtract_dm
    module procedure mp_subtract_ms
    module procedure mp_subtract_sm
    module procedure mp_subtract_ml
    module procedure mp_subtract_lm
  end interface
  interface operator(*)
    module procedure mp_multiply_mm
    module procedure mp_multiply_mi
    module procedure mp_multiply_im
    module procedure mp_multiply_md
    module procedure mp_multiply_dm
    module procedure mp_multiply_ms
    module procedure mp_multiply_sm
    module procedure mp_multiply_ml
    module procedure mp_multiply_lm
  end interface
  interface operator(/)
    module procedure mp_divide_mm
    module procedure mp_divide_mi
    module procedure mp_divide_im
    module procedure mp_divide_md
    module procedure mp_divide_dm
    module procedure mp_divide_ms
    module procedure mp_divide_sm
    module procedure mp_divide_ml
    module procedure mp_divide_lm
  end interface
  interface operator(**)
    module procedure mp_power_mm
    module procedure mp_power_mi
    module procedure mp_power_im
    module procedure mp_power_md
    module procedure mp_power_dm
    module procedure mp_power_ms
    module procedure mp_power_sm
    module procedure mp_power_ml
    module procedure mp_power_lm
  end interface
  interface operator(==)
    module procedure mp_equal_mm
    module procedure mp_equal_mi
    module procedure mp_equal_im
    module procedure mp_equal_md
    module procedure mp_equal_dm
    module procedure mp_equal_ms
    module procedure mp_equal_sm
    module procedure mp_equal_ml
    module procedure mp_equal_lm
  end interface
  interface operator(/=)
    module procedure mp_unequal_mm
    module procedure mp_unequal_mi
    module procedure mp_unequal_im
    module procedure mp_unequal_md
    module procedure mp_unequal_dm
    module procedure mp_unequal_ms
    module procedure mp_unequal_sm
    module procedure mp_unequal_ml
    module procedure mp_unequal_lm
  end interface
  interface operator(<)
    module procedure mp_less_mm
    module procedure mp_less_mi
    module procedure mp_less_im
    module procedure mp_less_md
    module procedure mp_less_dm
    module procedure mp_less_ms
    module procedure mp_less_sm
    module procedure mp_less_ml
    module procedure mp_less_lm
  end interface
  interface operator(<=)
    module procedure mp_less_equal_mm
    module procedure mp_less_equal_mi
    module procedure mp_less_equal_im
    module procedure mp_less_equal_md
    module procedure mp_less_equal_dm
    module procedure mp_less_equal_ms
    module procedure mp_less_equal_sm
    module procedure mp_less_equal_ml
    module procedure mp_less_equal_lm
  end interface
  interface operator(>)
    module procedure mp_greater_mm
    module procedure mp_greater_mi
    module procedure mp_greater_im
    module procedure mp_greater_md
    module procedure mp_greater_dm
    module procedure mp_greater_ms
    module procedure mp_greater_sm
    module procedure mp_greater_ml
    module procedure mp_greater_lm
  end interface
  interface operator(>=)
    module procedure mp_greater_equal_mm
    module procedure mp_greater_equal_mi
    module procedure mp_greater_equal_im
    module procedure mp_greater_equal_md
    module procedure mp_greater_equal_dm
    module procedure mp_greater_equal_ms
    module procedure mp_greater_equal_sm
    module procedure mp_greater_equal_ml
    module procedure mp_greater_equal_lm
  end interface
  interface atan2
    module procedure mp_atan2_mm
    module procedure mp_atan2_mi
    module procedure mp_atan2_im
    module procedure mp_atan2_md
    module procedure mp_atan2_dm
    module procedure mp_atan2_ms
    module procedure mp_atan2_sm
    module procedure mp_atan2_ml
    module procedure mp_atan2_lm
  end interface
  interface mod
    module procedure mp_mod_mm
    module procedure mp_mod_mi
    module procedure mp_mod_im
    module procedure mp_mod_md
    module procedure mp_mod_dm
    module procedure mp_mod_ms
    module procedure mp_mod_sm
    module procedure mp_mod_ml
    module procedure mp_mod_lm
  end interface
  interface sign
    module procedure mp_sign_mm
    module procedure mp_sign_mi
    module procedure mp_sign_im
    module procedure mp_sign_md
    module procedure mp_sign_dm
    module procedure mp_sign_ms
    module procedure mp_sign_sm
    module procedure mp_sign_ml
    module procedure mp_sign_lm
  end interface
  interface abs
    module procedure mp_abs
  end interface
  interface sqrt
    module procedure mp_sqrt
  end interface
  interface exp
    module procedure mp_exp
  end interface
  interface log
    module procedure mp_log
  end interface
  interface log10
    module procedure mp_log10
  end interface
  interface sin
    module procedure mp_sin
  end interface
  interface cos
    module procedure mp_cos
  end interface
  interface sinh
    module procedure mp_sinh
  end interface
  interface cosh
    module procedure mp_cosh
  end interface
  interface atan
    module procedure mp_atan
  end interface
  interface asin
    module procedure mp_asin
  end interface
  interface acos
    module procedure mp_acos
  end interface
  interface tan
    module procedure mp_tan
  end interface
  interface tanh
    module procedure mp_tanh
  end interface
  interface epsilon
    module procedure mp_epsilon
  end interface
  interface huge
    module procedure mp_huge
  end interface
  interface tiny
    module procedure mp_tiny
  end interface
  interface precision
    module procedure mp_precision
  end interface
  interface range
    module procedure mp_range
  end interface
  interface int
    module procedure mp_int
  end interface
  interface nint
    module procedure mp_nint
  end interface
  interface real
    module procedure mp_to_double
  end interface
  interface dble
    module procedure mp_to_double
  end interface
  interface max
    module procedure mp_max_m
    module procedure mp_max_im
    module procedure mp_max_dm
    module procedure mp_max_sm
    module procedure mp_max_iim
    module procedure mp_max_idm
    module procedure mp_max_ism
    module procedure mp_max_dim
    module procedure mp_max_ddm
    module procedure mp_max_dsm
    module procedure mp_max_sim
    module procedure mp_max_sdm
    module procedure mp_max_ssm
    module procedure mp_max_iiim
    module procedure mp_max_iidm
    module procedure mp_max_iism
    module procedure mp_max_idim
    module procedure mp_max_iddm
    module procedure mp_max_idsm
    module procedure mp_max_isim
    module procedure mp_max_isdm
    module procedure mp_max_issm
    module procedure mp_max_diim
    module procedure mp_max_didm
    module procedure mp_max_dism
    module procedure mp_max_ddim
    module procedure mp_max_dddm
    module procedure mp_max_ddsm
    module procedure mp_max_dsim
    module procedure mp_max_dsdm
    module procedure mp_max_dssm
    module procedure mp_max_siim
    module procedure mp_max_sidm
    module procedure mp_max_sism
    module procedure mp_max_sdim
    module procedure mp_max_sddm
    module procedure mp_max_sdsm
    module procedure mp_max_ssim
    module procedure mp_max_ssdm
    module procedure mp_max_sssm
  end interface
  interface min
    module procedure mp_min_m
    module procedure mp_min_im
    module procedure mp_min_dm
    module procedure mp_min_sm
    module procedure mp_min_iim
    module procedure mp_min_idm
    module procedure mp_min_ism
    module procedure mp_min_dim
    module procedure mp_min_ddm
    module procedure mp_min_dsm
    module procedure mp_min_sim
    module procedure mp_min_sdm
    module procedure mp_min_ssm
    module procedure mp_min_iiim
    module procedure mp_min_iidm
    module procedure mp_min_iism
    module procedure mp_min_idim
    module procedure mp_min_iddm
    module procedure mp_min_idsm
    module procedure mp_min_isim
    module procedure mp_min_isdm
    module procedure mp_min_issm
    module procedure mp_min_diim
    module procedure mp_min_didm
    module procedure mp_min_dism
    module procedure mp_min_ddim
    module procedure mp_min_dddm
    module procedure mp_min_ddsm
    module procedure mp_min_dsim
    module procedure mp_min_dsdm
    module procedure mp_min_dssm
    module procedure mp_min_siim
    module procedure mp_min_sidm
    module procedure mp_min_sism
    module procedure mp_min_sdim
    module procedure mp_min_sddm
    module procedure mp_min_sdsm
    module procedure mp_min_ssim
    module procedure mp_min_ssdm
    module procedure mp_min_sssm
  end interface
  interface sum
    module procedure mp_sum
  end interface
  interface maxval
    module procedure mp_maxval
  end interface
  interface minval
    module procedure mp_minval
  end interface
contains
  subroutine mp_init(x)
    type(mp_real), intent(out) :: x
    allocate(x%data(c_words(mpfr_bits)))
    call c_init(x%data,mpfr_bits)
  end subroutine

  impure elemental function mp_from_any(x,kind) result(y)
    class(*), intent(in) :: x
    integer, intent(in), optional :: kind
    type(mp_real) :: y
    integer(c_int) :: status
    call mp_init(y)
    select type(x)
    type is(mp_real)
      if (allocated(x%data)) call c_unary(16,y%data,x%data,mpfr_rounding)
    type is(integer)
      call c_integer(y%data,int(x,c_int64_t),mpfr_rounding)
    type is(integer(c_int64_t))
      call c_integer(y%data,x,mpfr_rounding)
    type is(real(c_float))
      call c_set_double(y%data,real(x,c_double),mpfr_rounding)
    type is(real(c_double))
      call c_set_double(y%data,x,mpfr_rounding)
    type is(character(len=*))
      status=c_string(y%data,trim(adjustl(x))//c_null_char,mpfr_rounding)
      if (status/=0) error stop 'Invalid MPFR numeric string'
    class default
      error stop 'Unsupported MPFR conversion'
    end select
  end function

  impure elemental function mp_binary_values(a,b,operation) result(y)
    type(mp_real), intent(in) :: a,b
    integer, intent(in) :: operation
    type(mp_real) :: y,z
    call mp_init(y)
    if (allocated(a%data).and.allocated(b%data)) then
      call c_binary(operation,y%data,a%data,b%data,mpfr_rounding)
    else
      z=mpreal(0)
      if (allocated(a%data)) then
        call c_binary(operation,y%data,a%data,z%data,mpfr_rounding)
      else if (allocated(b%data)) then
        call c_binary(operation,y%data,z%data,b%data,mpfr_rounding)
      else
        call c_binary(operation,y%data,z%data,z%data,mpfr_rounding)
      end if
    end if
  end function

  impure elemental function mp_unary_value(a,operation) result(y)
    type(mp_real), intent(in) :: a
    integer, intent(in) :: operation
    type(mp_real) :: y,z
    call mp_init(y)
    if (allocated(a%data)) then
      call c_unary(operation,y%data,a%data,mpfr_rounding)
    else
      z=mpreal(0)
      call c_unary(operation,y%data,z%data,mpfr_rounding)
    end if
  end function

  impure elemental integer function mp_compare_values(a,b) result(order)
    type(mp_real), intent(in) :: a,b
    type(mp_real) :: z
    if (allocated(a%data).and.allocated(b%data)) then
      order=c_compare(a%data,b%data)
    else
      z=mpreal(0)
      if (allocated(a%data)) then
        order=c_compare(a%data,z%data)
      else if (allocated(b%data)) then
        order=c_compare(z%data,b%data)
      else
        order=0
      end if
    end if
  end function

  impure elemental function mp_from_m(x,kind) result(y)
    type(mp_real), intent(in) :: x
    integer, intent(in), optional :: kind
    type(mp_real) :: y
    y=mp_from_any(x,kind)
  end function

  impure elemental function mp_from_i(x,kind) result(y)
    integer, intent(in) :: x
    integer, intent(in), optional :: kind
    type(mp_real) :: y
    y=mp_from_any(x,kind)
  end function

  impure elemental function mp_from_d(x,kind) result(y)
    real(c_double), intent(in) :: x
    integer, intent(in), optional :: kind
    type(mp_real) :: y
    y=mp_from_any(x,kind)
  end function

  impure elemental function mp_from_s(x,kind) result(y)
    real(c_float), intent(in) :: x
    integer, intent(in), optional :: kind
    type(mp_real) :: y
    y=mp_from_any(x,kind)
  end function

  impure elemental function mp_from_l(x,kind) result(y)
    integer(c_int64_t), intent(in) :: x
    integer, intent(in), optional :: kind
    type(mp_real) :: y
    y=mp_from_any(x,kind)
  end function

  impure elemental function mp_from_c(x,kind) result(y)
    character(len=*), intent(in) :: x
    integer, intent(in), optional :: kind
    type(mp_real) :: y
    y=mp_from_any(x,kind)
  end function

  impure elemental subroutine mp_assign_i(lhs,rhs)
    type(mp_real), intent(out) :: lhs
    integer, intent(in) :: rhs
    lhs=mpreal(rhs)
  end subroutine

  impure elemental subroutine mp_assign_d(lhs,rhs)
    type(mp_real), intent(out) :: lhs
    real(c_double), intent(in) :: rhs
    lhs=mpreal(rhs)
  end subroutine

  impure elemental subroutine mp_assign_s(lhs,rhs)
    type(mp_real), intent(out) :: lhs
    real(c_float), intent(in) :: rhs
    lhs=mpreal(rhs)
  end subroutine

  impure elemental subroutine mp_assign_l(lhs,rhs)
    type(mp_real), intent(out) :: lhs
    integer(c_int64_t), intent(in) :: rhs
    lhs=mpreal(rhs)
  end subroutine

  impure elemental subroutine mp_assign_c(lhs,rhs)
    type(mp_real), intent(out) :: lhs
    character(len=*), intent(in) :: rhs
    lhs=mpreal(rhs)
  end subroutine

  impure elemental subroutine mp_assign_to_i(lhs,rhs)
    integer, intent(out) :: lhs
    type(mp_real), intent(in) :: rhs
    lhs=int(rhs)
  end subroutine

  impure elemental subroutine mp_assign_to_d(lhs,rhs)
    real(c_double), intent(out) :: lhs
    type(mp_real), intent(in) :: rhs
    lhs=mp_to_double(rhs)
  end subroutine

  impure elemental subroutine mp_assign_to_s(lhs,rhs)
    real(c_float), intent(out) :: lhs
    type(mp_real), intent(in) :: rhs
    lhs=real(mp_to_double(rhs),c_float)
  end subroutine

  impure elemental subroutine mp_assign_to_l(lhs,rhs)
    integer(c_int64_t), intent(out) :: lhs
    type(mp_real), intent(in) :: rhs
    lhs=mp_to_integer64(rhs)
  end subroutine

  impure elemental function mp_positive(a) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real) :: y
    y=mp_unary_value(a,16)
  end function

  impure elemental function mp_negative(a) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real) :: y
    y=mp_unary_value(a,1)
  end function

  impure elemental function mp_add_mm(a,b) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,b,1)
  end function

  impure elemental function mp_add_mi(a,b) result(y)
    type(mp_real), intent(in) :: a
    integer, intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,mpreal(b),1)
  end function

  impure elemental function mp_add_im(a,b) result(y)
    integer, intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(mpreal(a),b,1)
  end function

  impure elemental function mp_add_md(a,b) result(y)
    type(mp_real), intent(in) :: a
    real(c_double), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,mpreal(b),1)
  end function

  impure elemental function mp_add_dm(a,b) result(y)
    real(c_double), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(mpreal(a),b,1)
  end function

  impure elemental function mp_add_ms(a,b) result(y)
    type(mp_real), intent(in) :: a
    real(c_float), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,mpreal(b),1)
  end function

  impure elemental function mp_add_sm(a,b) result(y)
    real(c_float), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(mpreal(a),b,1)
  end function

  impure elemental function mp_add_ml(a,b) result(y)
    type(mp_real), intent(in) :: a
    integer(c_int64_t), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,mpreal(b),1)
  end function

  impure elemental function mp_add_lm(a,b) result(y)
    integer(c_int64_t), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(mpreal(a),b,1)
  end function

  impure elemental function mp_subtract_mm(a,b) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,b,2)
  end function

  impure elemental function mp_subtract_mi(a,b) result(y)
    type(mp_real), intent(in) :: a
    integer, intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,mpreal(b),2)
  end function

  impure elemental function mp_subtract_im(a,b) result(y)
    integer, intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(mpreal(a),b,2)
  end function

  impure elemental function mp_subtract_md(a,b) result(y)
    type(mp_real), intent(in) :: a
    real(c_double), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,mpreal(b),2)
  end function

  impure elemental function mp_subtract_dm(a,b) result(y)
    real(c_double), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(mpreal(a),b,2)
  end function

  impure elemental function mp_subtract_ms(a,b) result(y)
    type(mp_real), intent(in) :: a
    real(c_float), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,mpreal(b),2)
  end function

  impure elemental function mp_subtract_sm(a,b) result(y)
    real(c_float), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(mpreal(a),b,2)
  end function

  impure elemental function mp_subtract_ml(a,b) result(y)
    type(mp_real), intent(in) :: a
    integer(c_int64_t), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,mpreal(b),2)
  end function

  impure elemental function mp_subtract_lm(a,b) result(y)
    integer(c_int64_t), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(mpreal(a),b,2)
  end function

  impure elemental function mp_multiply_mm(a,b) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,b,3)
  end function

  impure elemental function mp_multiply_mi(a,b) result(y)
    type(mp_real), intent(in) :: a
    integer, intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,mpreal(b),3)
  end function

  impure elemental function mp_multiply_im(a,b) result(y)
    integer, intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(mpreal(a),b,3)
  end function

  impure elemental function mp_multiply_md(a,b) result(y)
    type(mp_real), intent(in) :: a
    real(c_double), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,mpreal(b),3)
  end function

  impure elemental function mp_multiply_dm(a,b) result(y)
    real(c_double), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(mpreal(a),b,3)
  end function

  impure elemental function mp_multiply_ms(a,b) result(y)
    type(mp_real), intent(in) :: a
    real(c_float), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,mpreal(b),3)
  end function

  impure elemental function mp_multiply_sm(a,b) result(y)
    real(c_float), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(mpreal(a),b,3)
  end function

  impure elemental function mp_multiply_ml(a,b) result(y)
    type(mp_real), intent(in) :: a
    integer(c_int64_t), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,mpreal(b),3)
  end function

  impure elemental function mp_multiply_lm(a,b) result(y)
    integer(c_int64_t), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(mpreal(a),b,3)
  end function

  impure elemental function mp_divide_mm(a,b) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,b,4)
  end function

  impure elemental function mp_divide_mi(a,b) result(y)
    type(mp_real), intent(in) :: a
    integer, intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,mpreal(b),4)
  end function

  impure elemental function mp_divide_im(a,b) result(y)
    integer, intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(mpreal(a),b,4)
  end function

  impure elemental function mp_divide_md(a,b) result(y)
    type(mp_real), intent(in) :: a
    real(c_double), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,mpreal(b),4)
  end function

  impure elemental function mp_divide_dm(a,b) result(y)
    real(c_double), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(mpreal(a),b,4)
  end function

  impure elemental function mp_divide_ms(a,b) result(y)
    type(mp_real), intent(in) :: a
    real(c_float), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,mpreal(b),4)
  end function

  impure elemental function mp_divide_sm(a,b) result(y)
    real(c_float), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(mpreal(a),b,4)
  end function

  impure elemental function mp_divide_ml(a,b) result(y)
    type(mp_real), intent(in) :: a
    integer(c_int64_t), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,mpreal(b),4)
  end function

  impure elemental function mp_divide_lm(a,b) result(y)
    integer(c_int64_t), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(mpreal(a),b,4)
  end function

  impure elemental function mp_power_mm(a,b) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,b,5)
  end function

  impure elemental function mp_power_mi(a,b) result(y)
    type(mp_real), intent(in) :: a
    integer, intent(in) :: b
    type(mp_real) :: y,z
    call mp_init(y)
    z=mpreal(a)
    call c_power(y%data,z%data,int(b,c_int64_t),mpfr_rounding)
  end function

  impure elemental function mp_power_im(a,b) result(y)
    integer, intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(mpreal(a),b,5)
  end function

  impure elemental function mp_power_md(a,b) result(y)
    type(mp_real), intent(in) :: a
    real(c_double), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,mpreal(b),5)
  end function

  impure elemental function mp_power_dm(a,b) result(y)
    real(c_double), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(mpreal(a),b,5)
  end function

  impure elemental function mp_power_ms(a,b) result(y)
    type(mp_real), intent(in) :: a
    real(c_float), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,mpreal(b),5)
  end function

  impure elemental function mp_power_sm(a,b) result(y)
    real(c_float), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(mpreal(a),b,5)
  end function

  impure elemental function mp_power_ml(a,b) result(y)
    type(mp_real), intent(in) :: a
    integer(c_int64_t), intent(in) :: b
    type(mp_real) :: y,z
    call mp_init(y)
    z=mpreal(a)
    call c_power(y%data,z%data,int(b,c_int64_t),mpfr_rounding)
  end function

  impure elemental function mp_power_lm(a,b) result(y)
    integer(c_int64_t), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(mpreal(a),b,5)
  end function

  impure elemental function mp_atan2_mm(a,b) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,b,6)
  end function

  impure elemental function mp_atan2_mi(a,b) result(y)
    type(mp_real), intent(in) :: a
    integer, intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,mpreal(b),6)
  end function

  impure elemental function mp_atan2_im(a,b) result(y)
    integer, intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(mpreal(a),b,6)
  end function

  impure elemental function mp_atan2_md(a,b) result(y)
    type(mp_real), intent(in) :: a
    real(c_double), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,mpreal(b),6)
  end function

  impure elemental function mp_atan2_dm(a,b) result(y)
    real(c_double), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(mpreal(a),b,6)
  end function

  impure elemental function mp_atan2_ms(a,b) result(y)
    type(mp_real), intent(in) :: a
    real(c_float), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,mpreal(b),6)
  end function

  impure elemental function mp_atan2_sm(a,b) result(y)
    real(c_float), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(mpreal(a),b,6)
  end function

  impure elemental function mp_atan2_ml(a,b) result(y)
    type(mp_real), intent(in) :: a
    integer(c_int64_t), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,mpreal(b),6)
  end function

  impure elemental function mp_atan2_lm(a,b) result(y)
    integer(c_int64_t), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(mpreal(a),b,6)
  end function

  impure elemental function mp_mod_mm(a,b) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,b,7)
  end function

  impure elemental function mp_mod_mi(a,b) result(y)
    type(mp_real), intent(in) :: a
    integer, intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,mpreal(b),7)
  end function

  impure elemental function mp_mod_im(a,b) result(y)
    integer, intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(mpreal(a),b,7)
  end function

  impure elemental function mp_mod_md(a,b) result(y)
    type(mp_real), intent(in) :: a
    real(c_double), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,mpreal(b),7)
  end function

  impure elemental function mp_mod_dm(a,b) result(y)
    real(c_double), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(mpreal(a),b,7)
  end function

  impure elemental function mp_mod_ms(a,b) result(y)
    type(mp_real), intent(in) :: a
    real(c_float), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,mpreal(b),7)
  end function

  impure elemental function mp_mod_sm(a,b) result(y)
    real(c_float), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(mpreal(a),b,7)
  end function

  impure elemental function mp_mod_ml(a,b) result(y)
    type(mp_real), intent(in) :: a
    integer(c_int64_t), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(a,mpreal(b),7)
  end function

  impure elemental function mp_mod_lm(a,b) result(y)
    integer(c_int64_t), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(mpreal(a),b,7)
  end function

  impure elemental function mp_sign_mm(a,b) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(abs(a),b,8)
  end function

  impure elemental function mp_sign_mi(a,b) result(y)
    type(mp_real), intent(in) :: a
    integer, intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(abs(a),mpreal(b),8)
  end function

  impure elemental function mp_sign_im(a,b) result(y)
    integer, intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(abs(mpreal(a)),b,8)
  end function

  impure elemental function mp_sign_md(a,b) result(y)
    type(mp_real), intent(in) :: a
    real(c_double), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(abs(a),mpreal(b),8)
  end function

  impure elemental function mp_sign_dm(a,b) result(y)
    real(c_double), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(abs(mpreal(a)),b,8)
  end function

  impure elemental function mp_sign_ms(a,b) result(y)
    type(mp_real), intent(in) :: a
    real(c_float), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(abs(a),mpreal(b),8)
  end function

  impure elemental function mp_sign_sm(a,b) result(y)
    real(c_float), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(abs(mpreal(a)),b,8)
  end function

  impure elemental function mp_sign_ml(a,b) result(y)
    type(mp_real), intent(in) :: a
    integer(c_int64_t), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(abs(a),mpreal(b),8)
  end function

  impure elemental function mp_sign_lm(a,b) result(y)
    integer(c_int64_t), intent(in) :: a
    type(mp_real), intent(in) :: b
    type(mp_real) :: y,z
    y=mp_binary_values(abs(mpreal(a)),b,8)
  end function

  impure elemental logical function mp_equal_mm(a,b) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(a,b)
    y=order==0
  end function

  impure elemental logical function mp_equal_mi(a,b) result(y)
    type(mp_real), intent(in) :: a
    integer, intent(in) :: b
    integer :: order
    order=mp_compare_values(a,mpreal(b))
    y=order==0
  end function

  impure elemental logical function mp_equal_im(a,b) result(y)
    integer, intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(mpreal(a),b)
    y=order==0
  end function

  impure elemental logical function mp_equal_md(a,b) result(y)
    type(mp_real), intent(in) :: a
    real(c_double), intent(in) :: b
    integer :: order
    order=mp_compare_values(a,mpreal(b))
    y=order==0
  end function

  impure elemental logical function mp_equal_dm(a,b) result(y)
    real(c_double), intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(mpreal(a),b)
    y=order==0
  end function

  impure elemental logical function mp_equal_ms(a,b) result(y)
    type(mp_real), intent(in) :: a
    real(c_float), intent(in) :: b
    integer :: order
    order=mp_compare_values(a,mpreal(b))
    y=order==0
  end function

  impure elemental logical function mp_equal_sm(a,b) result(y)
    real(c_float), intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(mpreal(a),b)
    y=order==0
  end function

  impure elemental logical function mp_equal_ml(a,b) result(y)
    type(mp_real), intent(in) :: a
    integer(c_int64_t), intent(in) :: b
    integer :: order
    order=mp_compare_values(a,mpreal(b))
    y=order==0
  end function

  impure elemental logical function mp_equal_lm(a,b) result(y)
    integer(c_int64_t), intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(mpreal(a),b)
    y=order==0
  end function

  impure elemental logical function mp_unequal_mm(a,b) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(a,b)
    y=order/=0
  end function

  impure elemental logical function mp_unequal_mi(a,b) result(y)
    type(mp_real), intent(in) :: a
    integer, intent(in) :: b
    integer :: order
    order=mp_compare_values(a,mpreal(b))
    y=order/=0
  end function

  impure elemental logical function mp_unequal_im(a,b) result(y)
    integer, intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(mpreal(a),b)
    y=order/=0
  end function

  impure elemental logical function mp_unequal_md(a,b) result(y)
    type(mp_real), intent(in) :: a
    real(c_double), intent(in) :: b
    integer :: order
    order=mp_compare_values(a,mpreal(b))
    y=order/=0
  end function

  impure elemental logical function mp_unequal_dm(a,b) result(y)
    real(c_double), intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(mpreal(a),b)
    y=order/=0
  end function

  impure elemental logical function mp_unequal_ms(a,b) result(y)
    type(mp_real), intent(in) :: a
    real(c_float), intent(in) :: b
    integer :: order
    order=mp_compare_values(a,mpreal(b))
    y=order/=0
  end function

  impure elemental logical function mp_unequal_sm(a,b) result(y)
    real(c_float), intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(mpreal(a),b)
    y=order/=0
  end function

  impure elemental logical function mp_unequal_ml(a,b) result(y)
    type(mp_real), intent(in) :: a
    integer(c_int64_t), intent(in) :: b
    integer :: order
    order=mp_compare_values(a,mpreal(b))
    y=order/=0
  end function

  impure elemental logical function mp_unequal_lm(a,b) result(y)
    integer(c_int64_t), intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(mpreal(a),b)
    y=order/=0
  end function

  impure elemental logical function mp_less_mm(a,b) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(a,b)
    y=order==-1
  end function

  impure elemental logical function mp_less_mi(a,b) result(y)
    type(mp_real), intent(in) :: a
    integer, intent(in) :: b
    integer :: order
    order=mp_compare_values(a,mpreal(b))
    y=order==-1
  end function

  impure elemental logical function mp_less_im(a,b) result(y)
    integer, intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(mpreal(a),b)
    y=order==-1
  end function

  impure elemental logical function mp_less_md(a,b) result(y)
    type(mp_real), intent(in) :: a
    real(c_double), intent(in) :: b
    integer :: order
    order=mp_compare_values(a,mpreal(b))
    y=order==-1
  end function

  impure elemental logical function mp_less_dm(a,b) result(y)
    real(c_double), intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(mpreal(a),b)
    y=order==-1
  end function

  impure elemental logical function mp_less_ms(a,b) result(y)
    type(mp_real), intent(in) :: a
    real(c_float), intent(in) :: b
    integer :: order
    order=mp_compare_values(a,mpreal(b))
    y=order==-1
  end function

  impure elemental logical function mp_less_sm(a,b) result(y)
    real(c_float), intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(mpreal(a),b)
    y=order==-1
  end function

  impure elemental logical function mp_less_ml(a,b) result(y)
    type(mp_real), intent(in) :: a
    integer(c_int64_t), intent(in) :: b
    integer :: order
    order=mp_compare_values(a,mpreal(b))
    y=order==-1
  end function

  impure elemental logical function mp_less_lm(a,b) result(y)
    integer(c_int64_t), intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(mpreal(a),b)
    y=order==-1
  end function

  impure elemental logical function mp_less_equal_mm(a,b) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(a,b)
    y=order==-1.or.order==0
  end function

  impure elemental logical function mp_less_equal_mi(a,b) result(y)
    type(mp_real), intent(in) :: a
    integer, intent(in) :: b
    integer :: order
    order=mp_compare_values(a,mpreal(b))
    y=order==-1.or.order==0
  end function

  impure elemental logical function mp_less_equal_im(a,b) result(y)
    integer, intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(mpreal(a),b)
    y=order==-1.or.order==0
  end function

  impure elemental logical function mp_less_equal_md(a,b) result(y)
    type(mp_real), intent(in) :: a
    real(c_double), intent(in) :: b
    integer :: order
    order=mp_compare_values(a,mpreal(b))
    y=order==-1.or.order==0
  end function

  impure elemental logical function mp_less_equal_dm(a,b) result(y)
    real(c_double), intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(mpreal(a),b)
    y=order==-1.or.order==0
  end function

  impure elemental logical function mp_less_equal_ms(a,b) result(y)
    type(mp_real), intent(in) :: a
    real(c_float), intent(in) :: b
    integer :: order
    order=mp_compare_values(a,mpreal(b))
    y=order==-1.or.order==0
  end function

  impure elemental logical function mp_less_equal_sm(a,b) result(y)
    real(c_float), intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(mpreal(a),b)
    y=order==-1.or.order==0
  end function

  impure elemental logical function mp_less_equal_ml(a,b) result(y)
    type(mp_real), intent(in) :: a
    integer(c_int64_t), intent(in) :: b
    integer :: order
    order=mp_compare_values(a,mpreal(b))
    y=order==-1.or.order==0
  end function

  impure elemental logical function mp_less_equal_lm(a,b) result(y)
    integer(c_int64_t), intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(mpreal(a),b)
    y=order==-1.or.order==0
  end function

  impure elemental logical function mp_greater_mm(a,b) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(a,b)
    y=order==1
  end function

  impure elemental logical function mp_greater_mi(a,b) result(y)
    type(mp_real), intent(in) :: a
    integer, intent(in) :: b
    integer :: order
    order=mp_compare_values(a,mpreal(b))
    y=order==1
  end function

  impure elemental logical function mp_greater_im(a,b) result(y)
    integer, intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(mpreal(a),b)
    y=order==1
  end function

  impure elemental logical function mp_greater_md(a,b) result(y)
    type(mp_real), intent(in) :: a
    real(c_double), intent(in) :: b
    integer :: order
    order=mp_compare_values(a,mpreal(b))
    y=order==1
  end function

  impure elemental logical function mp_greater_dm(a,b) result(y)
    real(c_double), intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(mpreal(a),b)
    y=order==1
  end function

  impure elemental logical function mp_greater_ms(a,b) result(y)
    type(mp_real), intent(in) :: a
    real(c_float), intent(in) :: b
    integer :: order
    order=mp_compare_values(a,mpreal(b))
    y=order==1
  end function

  impure elemental logical function mp_greater_sm(a,b) result(y)
    real(c_float), intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(mpreal(a),b)
    y=order==1
  end function

  impure elemental logical function mp_greater_ml(a,b) result(y)
    type(mp_real), intent(in) :: a
    integer(c_int64_t), intent(in) :: b
    integer :: order
    order=mp_compare_values(a,mpreal(b))
    y=order==1
  end function

  impure elemental logical function mp_greater_lm(a,b) result(y)
    integer(c_int64_t), intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(mpreal(a),b)
    y=order==1
  end function

  impure elemental logical function mp_greater_equal_mm(a,b) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(a,b)
    y=order==1.or.order==0
  end function

  impure elemental logical function mp_greater_equal_mi(a,b) result(y)
    type(mp_real), intent(in) :: a
    integer, intent(in) :: b
    integer :: order
    order=mp_compare_values(a,mpreal(b))
    y=order==1.or.order==0
  end function

  impure elemental logical function mp_greater_equal_im(a,b) result(y)
    integer, intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(mpreal(a),b)
    y=order==1.or.order==0
  end function

  impure elemental logical function mp_greater_equal_md(a,b) result(y)
    type(mp_real), intent(in) :: a
    real(c_double), intent(in) :: b
    integer :: order
    order=mp_compare_values(a,mpreal(b))
    y=order==1.or.order==0
  end function

  impure elemental logical function mp_greater_equal_dm(a,b) result(y)
    real(c_double), intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(mpreal(a),b)
    y=order==1.or.order==0
  end function

  impure elemental logical function mp_greater_equal_ms(a,b) result(y)
    type(mp_real), intent(in) :: a
    real(c_float), intent(in) :: b
    integer :: order
    order=mp_compare_values(a,mpreal(b))
    y=order==1.or.order==0
  end function

  impure elemental logical function mp_greater_equal_sm(a,b) result(y)
    real(c_float), intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(mpreal(a),b)
    y=order==1.or.order==0
  end function

  impure elemental logical function mp_greater_equal_ml(a,b) result(y)
    type(mp_real), intent(in) :: a
    integer(c_int64_t), intent(in) :: b
    integer :: order
    order=mp_compare_values(a,mpreal(b))
    y=order==1.or.order==0
  end function

  impure elemental logical function mp_greater_equal_lm(a,b) result(y)
    integer(c_int64_t), intent(in) :: a
    type(mp_real), intent(in) :: b
    integer :: order
    order=mp_compare_values(mpreal(a),b)
    y=order==1.or.order==0
  end function

  impure elemental function mp_abs(a) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real) :: y
    y=mp_unary_value(a,2)
  end function

  impure elemental function mp_sqrt(a) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real) :: y
    y=mp_unary_value(a,3)
  end function

  impure elemental function mp_exp(a) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real) :: y
    y=mp_unary_value(a,4)
  end function

  impure elemental function mp_log(a) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real) :: y
    y=mp_unary_value(a,5)
  end function

  impure elemental function mp_log10(a) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real) :: y
    y=mp_unary_value(a,6)
  end function

  impure elemental function mp_sin(a) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real) :: y
    y=mp_unary_value(a,7)
  end function

  impure elemental function mp_cos(a) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real) :: y
    y=mp_unary_value(a,8)
  end function

  impure elemental function mp_sinh(a) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real) :: y
    y=mp_unary_value(a,9)
  end function

  impure elemental function mp_cosh(a) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real) :: y
    y=mp_unary_value(a,10)
  end function

  impure elemental function mp_atan(a) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real) :: y
    y=mp_unary_value(a,11)
  end function

  impure elemental function mp_asin(a) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real) :: y
    y=mp_unary_value(a,12)
  end function

  impure elemental function mp_acos(a) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real) :: y
    y=mp_unary_value(a,13)
  end function

  impure elemental function mp_tan(a) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real) :: y
    y=mp_unary_value(a,14)
  end function

  impure elemental function mp_tanh(a) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real) :: y
    y=mp_unary_value(a,15)
  end function

  impure elemental function mp_epsilon(a) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real) :: y
    call mp_init(y)
    call c_constant(1,y%data,mpfr_rounding)
  end function

  impure elemental function mp_huge(a) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real) :: y
    call mp_init(y)
    call c_constant(2,y%data,mpfr_rounding)
  end function

  impure elemental function mp_tiny(a) result(y)
    type(mp_real), intent(in) :: a
    type(mp_real) :: y
    call mp_init(y)
    call c_constant(3,y%data,mpfr_rounding)
  end function

  impure elemental integer function mp_precision(a) result(y)
    type(mp_real), intent(in) :: a
    y=int(real(mpfr_bits-1,c_double)*log10(2.0_c_double))
  end function

  impure elemental integer function mp_range(a) result(y)
    type(mp_real), intent(in) :: a
    y=int(real(1073741822,c_double)*log10(2.0_c_double))
  end function

  impure elemental integer function mp_int(a,kind) result(y)
    type(mp_real), intent(in) :: a
    integer, intent(in), optional :: kind
    if (allocated(a%data)) then
      y=int(c_to_integer(a%data,1))
    else
      y=0
    end if
  end function

  impure elemental integer function mp_nint(a,kind) result(y)
    type(mp_real), intent(in) :: a
    integer, intent(in), optional :: kind
    if (allocated(a%data)) then
      y=int(c_to_integer(a%data,-1))
    else
      y=0
    end if
  end function

  impure elemental real(c_double) function mp_to_double(a,kind) result(y)
    type(mp_real), intent(in) :: a
    integer, intent(in), optional :: kind
    if (allocated(a%data)) then
      y=c_to_double(a%data,mpfr_rounding)
    else
      y=0.0_c_double
    end if
  end function

  impure elemental integer(c_int64_t) function mp_to_integer64(a) result(y)
    type(mp_real), intent(in) :: a
    if (allocated(a%data)) then
      y=c_to_integer(a%data,1)
    else
      y=0_c_int64_t
    end if
  end function

  impure elemental logical function mp_is_finite(a) result(finite)
    type(mp_real), intent(in) :: a
    finite=.true.
    if (allocated(a%data)) finite=c_is_finite(a%data)/=0
  end function

  function mp_to_string(a) result(text)
    type(mp_real), intent(in) :: a
    character(len=:), allocatable :: text
    character(c_char), allocatable :: buffer(:)
    integer(c_int) :: count,capacity
    integer :: i
    if (.not.allocated(a%data)) then
      text='0'
      return
    end if
    capacity=int(a%data(1)/3+128,c_int)
    allocate(buffer(capacity))
    count=c_format(a%data,buffer,capacity)
    if (count>=capacity.or.count<0) error stop 'MPFR text output overflow'
    allocate(character(len=count) :: text)
    do i=1,count
      text(i:i)=buffer(i)
    end do
  end function

  impure elemental function mp_max_m(a,b,c,d) result(y)
    type(mp_real), intent(in) :: a
    class(*), intent(in) :: b
    class(*), intent(in), optional :: c
    class(*), intent(in), optional :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
    if (present(c)) then
      value=mp_from_any(c)
      if (value>y) y=value
    end if
    if (present(d)) then
      value=mp_from_any(d)
      if (value>y) y=value
    end if
  end function

  impure elemental function mp_max_im(a,b,c,d) result(y)
    integer, intent(in) :: a
    type(mp_real), intent(in) :: b
    class(*), intent(in), optional :: c
    class(*), intent(in), optional :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
    if (present(c)) then
      value=mp_from_any(c)
      if (value>y) y=value
    end if
    if (present(d)) then
      value=mp_from_any(d)
      if (value>y) y=value
    end if
  end function

  impure elemental function mp_max_dm(a,b,c,d) result(y)
    real(c_double), intent(in) :: a
    type(mp_real), intent(in) :: b
    class(*), intent(in), optional :: c
    class(*), intent(in), optional :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
    if (present(c)) then
      value=mp_from_any(c)
      if (value>y) y=value
    end if
    if (present(d)) then
      value=mp_from_any(d)
      if (value>y) y=value
    end if
  end function

  impure elemental function mp_max_sm(a,b,c,d) result(y)
    real(c_float), intent(in) :: a
    type(mp_real), intent(in) :: b
    class(*), intent(in), optional :: c
    class(*), intent(in), optional :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
    if (present(c)) then
      value=mp_from_any(c)
      if (value>y) y=value
    end if
    if (present(d)) then
      value=mp_from_any(d)
      if (value>y) y=value
    end if
  end function

  impure elemental function mp_max_iim(a,b,c,d) result(y)
    integer, intent(in) :: a
    integer, intent(in) :: b
    type(mp_real), intent(in) :: c
    class(*), intent(in), optional :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
    if (present(d)) then
      value=mp_from_any(d)
      if (value>y) y=value
    end if
  end function

  impure elemental function mp_max_idm(a,b,c,d) result(y)
    integer, intent(in) :: a
    real(c_double), intent(in) :: b
    type(mp_real), intent(in) :: c
    class(*), intent(in), optional :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
    if (present(d)) then
      value=mp_from_any(d)
      if (value>y) y=value
    end if
  end function

  impure elemental function mp_max_ism(a,b,c,d) result(y)
    integer, intent(in) :: a
    real(c_float), intent(in) :: b
    type(mp_real), intent(in) :: c
    class(*), intent(in), optional :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
    if (present(d)) then
      value=mp_from_any(d)
      if (value>y) y=value
    end if
  end function

  impure elemental function mp_max_dim(a,b,c,d) result(y)
    real(c_double), intent(in) :: a
    integer, intent(in) :: b
    type(mp_real), intent(in) :: c
    class(*), intent(in), optional :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
    if (present(d)) then
      value=mp_from_any(d)
      if (value>y) y=value
    end if
  end function

  impure elemental function mp_max_ddm(a,b,c,d) result(y)
    real(c_double), intent(in) :: a
    real(c_double), intent(in) :: b
    type(mp_real), intent(in) :: c
    class(*), intent(in), optional :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
    if (present(d)) then
      value=mp_from_any(d)
      if (value>y) y=value
    end if
  end function

  impure elemental function mp_max_dsm(a,b,c,d) result(y)
    real(c_double), intent(in) :: a
    real(c_float), intent(in) :: b
    type(mp_real), intent(in) :: c
    class(*), intent(in), optional :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
    if (present(d)) then
      value=mp_from_any(d)
      if (value>y) y=value
    end if
  end function

  impure elemental function mp_max_sim(a,b,c,d) result(y)
    real(c_float), intent(in) :: a
    integer, intent(in) :: b
    type(mp_real), intent(in) :: c
    class(*), intent(in), optional :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
    if (present(d)) then
      value=mp_from_any(d)
      if (value>y) y=value
    end if
  end function

  impure elemental function mp_max_sdm(a,b,c,d) result(y)
    real(c_float), intent(in) :: a
    real(c_double), intent(in) :: b
    type(mp_real), intent(in) :: c
    class(*), intent(in), optional :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
    if (present(d)) then
      value=mp_from_any(d)
      if (value>y) y=value
    end if
  end function

  impure elemental function mp_max_ssm(a,b,c,d) result(y)
    real(c_float), intent(in) :: a
    real(c_float), intent(in) :: b
    type(mp_real), intent(in) :: c
    class(*), intent(in), optional :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
    if (present(d)) then
      value=mp_from_any(d)
      if (value>y) y=value
    end if
  end function

  impure elemental function mp_max_iiim(a,b,c,d) result(y)
    integer, intent(in) :: a
    integer, intent(in) :: b
    integer, intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
      value=mp_from_any(d)
      if (value>y) y=value
  end function

  impure elemental function mp_max_iidm(a,b,c,d) result(y)
    integer, intent(in) :: a
    integer, intent(in) :: b
    real(c_double), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
      value=mp_from_any(d)
      if (value>y) y=value
  end function

  impure elemental function mp_max_iism(a,b,c,d) result(y)
    integer, intent(in) :: a
    integer, intent(in) :: b
    real(c_float), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
      value=mp_from_any(d)
      if (value>y) y=value
  end function

  impure elemental function mp_max_idim(a,b,c,d) result(y)
    integer, intent(in) :: a
    real(c_double), intent(in) :: b
    integer, intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
      value=mp_from_any(d)
      if (value>y) y=value
  end function

  impure elemental function mp_max_iddm(a,b,c,d) result(y)
    integer, intent(in) :: a
    real(c_double), intent(in) :: b
    real(c_double), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
      value=mp_from_any(d)
      if (value>y) y=value
  end function

  impure elemental function mp_max_idsm(a,b,c,d) result(y)
    integer, intent(in) :: a
    real(c_double), intent(in) :: b
    real(c_float), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
      value=mp_from_any(d)
      if (value>y) y=value
  end function

  impure elemental function mp_max_isim(a,b,c,d) result(y)
    integer, intent(in) :: a
    real(c_float), intent(in) :: b
    integer, intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
      value=mp_from_any(d)
      if (value>y) y=value
  end function

  impure elemental function mp_max_isdm(a,b,c,d) result(y)
    integer, intent(in) :: a
    real(c_float), intent(in) :: b
    real(c_double), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
      value=mp_from_any(d)
      if (value>y) y=value
  end function

  impure elemental function mp_max_issm(a,b,c,d) result(y)
    integer, intent(in) :: a
    real(c_float), intent(in) :: b
    real(c_float), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
      value=mp_from_any(d)
      if (value>y) y=value
  end function

  impure elemental function mp_max_diim(a,b,c,d) result(y)
    real(c_double), intent(in) :: a
    integer, intent(in) :: b
    integer, intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
      value=mp_from_any(d)
      if (value>y) y=value
  end function

  impure elemental function mp_max_didm(a,b,c,d) result(y)
    real(c_double), intent(in) :: a
    integer, intent(in) :: b
    real(c_double), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
      value=mp_from_any(d)
      if (value>y) y=value
  end function

  impure elemental function mp_max_dism(a,b,c,d) result(y)
    real(c_double), intent(in) :: a
    integer, intent(in) :: b
    real(c_float), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
      value=mp_from_any(d)
      if (value>y) y=value
  end function

  impure elemental function mp_max_ddim(a,b,c,d) result(y)
    real(c_double), intent(in) :: a
    real(c_double), intent(in) :: b
    integer, intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
      value=mp_from_any(d)
      if (value>y) y=value
  end function

  impure elemental function mp_max_dddm(a,b,c,d) result(y)
    real(c_double), intent(in) :: a
    real(c_double), intent(in) :: b
    real(c_double), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
      value=mp_from_any(d)
      if (value>y) y=value
  end function

  impure elemental function mp_max_ddsm(a,b,c,d) result(y)
    real(c_double), intent(in) :: a
    real(c_double), intent(in) :: b
    real(c_float), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
      value=mp_from_any(d)
      if (value>y) y=value
  end function

  impure elemental function mp_max_dsim(a,b,c,d) result(y)
    real(c_double), intent(in) :: a
    real(c_float), intent(in) :: b
    integer, intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
      value=mp_from_any(d)
      if (value>y) y=value
  end function

  impure elemental function mp_max_dsdm(a,b,c,d) result(y)
    real(c_double), intent(in) :: a
    real(c_float), intent(in) :: b
    real(c_double), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
      value=mp_from_any(d)
      if (value>y) y=value
  end function

  impure elemental function mp_max_dssm(a,b,c,d) result(y)
    real(c_double), intent(in) :: a
    real(c_float), intent(in) :: b
    real(c_float), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
      value=mp_from_any(d)
      if (value>y) y=value
  end function

  impure elemental function mp_max_siim(a,b,c,d) result(y)
    real(c_float), intent(in) :: a
    integer, intent(in) :: b
    integer, intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
      value=mp_from_any(d)
      if (value>y) y=value
  end function

  impure elemental function mp_max_sidm(a,b,c,d) result(y)
    real(c_float), intent(in) :: a
    integer, intent(in) :: b
    real(c_double), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
      value=mp_from_any(d)
      if (value>y) y=value
  end function

  impure elemental function mp_max_sism(a,b,c,d) result(y)
    real(c_float), intent(in) :: a
    integer, intent(in) :: b
    real(c_float), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
      value=mp_from_any(d)
      if (value>y) y=value
  end function

  impure elemental function mp_max_sdim(a,b,c,d) result(y)
    real(c_float), intent(in) :: a
    real(c_double), intent(in) :: b
    integer, intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
      value=mp_from_any(d)
      if (value>y) y=value
  end function

  impure elemental function mp_max_sddm(a,b,c,d) result(y)
    real(c_float), intent(in) :: a
    real(c_double), intent(in) :: b
    real(c_double), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
      value=mp_from_any(d)
      if (value>y) y=value
  end function

  impure elemental function mp_max_sdsm(a,b,c,d) result(y)
    real(c_float), intent(in) :: a
    real(c_double), intent(in) :: b
    real(c_float), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
      value=mp_from_any(d)
      if (value>y) y=value
  end function

  impure elemental function mp_max_ssim(a,b,c,d) result(y)
    real(c_float), intent(in) :: a
    real(c_float), intent(in) :: b
    integer, intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
      value=mp_from_any(d)
      if (value>y) y=value
  end function

  impure elemental function mp_max_ssdm(a,b,c,d) result(y)
    real(c_float), intent(in) :: a
    real(c_float), intent(in) :: b
    real(c_double), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
      value=mp_from_any(d)
      if (value>y) y=value
  end function

  impure elemental function mp_max_sssm(a,b,c,d) result(y)
    real(c_float), intent(in) :: a
    real(c_float), intent(in) :: b
    real(c_float), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value>y) y=value
      value=mp_from_any(c)
      if (value>y) y=value
      value=mp_from_any(d)
      if (value>y) y=value
  end function

  impure elemental function mp_min_m(a,b,c,d) result(y)
    type(mp_real), intent(in) :: a
    class(*), intent(in) :: b
    class(*), intent(in), optional :: c
    class(*), intent(in), optional :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
    if (present(c)) then
      value=mp_from_any(c)
      if (value<y) y=value
    end if
    if (present(d)) then
      value=mp_from_any(d)
      if (value<y) y=value
    end if
  end function

  impure elemental function mp_min_im(a,b,c,d) result(y)
    integer, intent(in) :: a
    type(mp_real), intent(in) :: b
    class(*), intent(in), optional :: c
    class(*), intent(in), optional :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
    if (present(c)) then
      value=mp_from_any(c)
      if (value<y) y=value
    end if
    if (present(d)) then
      value=mp_from_any(d)
      if (value<y) y=value
    end if
  end function

  impure elemental function mp_min_dm(a,b,c,d) result(y)
    real(c_double), intent(in) :: a
    type(mp_real), intent(in) :: b
    class(*), intent(in), optional :: c
    class(*), intent(in), optional :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
    if (present(c)) then
      value=mp_from_any(c)
      if (value<y) y=value
    end if
    if (present(d)) then
      value=mp_from_any(d)
      if (value<y) y=value
    end if
  end function

  impure elemental function mp_min_sm(a,b,c,d) result(y)
    real(c_float), intent(in) :: a
    type(mp_real), intent(in) :: b
    class(*), intent(in), optional :: c
    class(*), intent(in), optional :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
    if (present(c)) then
      value=mp_from_any(c)
      if (value<y) y=value
    end if
    if (present(d)) then
      value=mp_from_any(d)
      if (value<y) y=value
    end if
  end function

  impure elemental function mp_min_iim(a,b,c,d) result(y)
    integer, intent(in) :: a
    integer, intent(in) :: b
    type(mp_real), intent(in) :: c
    class(*), intent(in), optional :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
    if (present(d)) then
      value=mp_from_any(d)
      if (value<y) y=value
    end if
  end function

  impure elemental function mp_min_idm(a,b,c,d) result(y)
    integer, intent(in) :: a
    real(c_double), intent(in) :: b
    type(mp_real), intent(in) :: c
    class(*), intent(in), optional :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
    if (present(d)) then
      value=mp_from_any(d)
      if (value<y) y=value
    end if
  end function

  impure elemental function mp_min_ism(a,b,c,d) result(y)
    integer, intent(in) :: a
    real(c_float), intent(in) :: b
    type(mp_real), intent(in) :: c
    class(*), intent(in), optional :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
    if (present(d)) then
      value=mp_from_any(d)
      if (value<y) y=value
    end if
  end function

  impure elemental function mp_min_dim(a,b,c,d) result(y)
    real(c_double), intent(in) :: a
    integer, intent(in) :: b
    type(mp_real), intent(in) :: c
    class(*), intent(in), optional :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
    if (present(d)) then
      value=mp_from_any(d)
      if (value<y) y=value
    end if
  end function

  impure elemental function mp_min_ddm(a,b,c,d) result(y)
    real(c_double), intent(in) :: a
    real(c_double), intent(in) :: b
    type(mp_real), intent(in) :: c
    class(*), intent(in), optional :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
    if (present(d)) then
      value=mp_from_any(d)
      if (value<y) y=value
    end if
  end function

  impure elemental function mp_min_dsm(a,b,c,d) result(y)
    real(c_double), intent(in) :: a
    real(c_float), intent(in) :: b
    type(mp_real), intent(in) :: c
    class(*), intent(in), optional :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
    if (present(d)) then
      value=mp_from_any(d)
      if (value<y) y=value
    end if
  end function

  impure elemental function mp_min_sim(a,b,c,d) result(y)
    real(c_float), intent(in) :: a
    integer, intent(in) :: b
    type(mp_real), intent(in) :: c
    class(*), intent(in), optional :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
    if (present(d)) then
      value=mp_from_any(d)
      if (value<y) y=value
    end if
  end function

  impure elemental function mp_min_sdm(a,b,c,d) result(y)
    real(c_float), intent(in) :: a
    real(c_double), intent(in) :: b
    type(mp_real), intent(in) :: c
    class(*), intent(in), optional :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
    if (present(d)) then
      value=mp_from_any(d)
      if (value<y) y=value
    end if
  end function

  impure elemental function mp_min_ssm(a,b,c,d) result(y)
    real(c_float), intent(in) :: a
    real(c_float), intent(in) :: b
    type(mp_real), intent(in) :: c
    class(*), intent(in), optional :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
    if (present(d)) then
      value=mp_from_any(d)
      if (value<y) y=value
    end if
  end function

  impure elemental function mp_min_iiim(a,b,c,d) result(y)
    integer, intent(in) :: a
    integer, intent(in) :: b
    integer, intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
      value=mp_from_any(d)
      if (value<y) y=value
  end function

  impure elemental function mp_min_iidm(a,b,c,d) result(y)
    integer, intent(in) :: a
    integer, intent(in) :: b
    real(c_double), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
      value=mp_from_any(d)
      if (value<y) y=value
  end function

  impure elemental function mp_min_iism(a,b,c,d) result(y)
    integer, intent(in) :: a
    integer, intent(in) :: b
    real(c_float), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
      value=mp_from_any(d)
      if (value<y) y=value
  end function

  impure elemental function mp_min_idim(a,b,c,d) result(y)
    integer, intent(in) :: a
    real(c_double), intent(in) :: b
    integer, intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
      value=mp_from_any(d)
      if (value<y) y=value
  end function

  impure elemental function mp_min_iddm(a,b,c,d) result(y)
    integer, intent(in) :: a
    real(c_double), intent(in) :: b
    real(c_double), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
      value=mp_from_any(d)
      if (value<y) y=value
  end function

  impure elemental function mp_min_idsm(a,b,c,d) result(y)
    integer, intent(in) :: a
    real(c_double), intent(in) :: b
    real(c_float), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
      value=mp_from_any(d)
      if (value<y) y=value
  end function

  impure elemental function mp_min_isim(a,b,c,d) result(y)
    integer, intent(in) :: a
    real(c_float), intent(in) :: b
    integer, intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
      value=mp_from_any(d)
      if (value<y) y=value
  end function

  impure elemental function mp_min_isdm(a,b,c,d) result(y)
    integer, intent(in) :: a
    real(c_float), intent(in) :: b
    real(c_double), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
      value=mp_from_any(d)
      if (value<y) y=value
  end function

  impure elemental function mp_min_issm(a,b,c,d) result(y)
    integer, intent(in) :: a
    real(c_float), intent(in) :: b
    real(c_float), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
      value=mp_from_any(d)
      if (value<y) y=value
  end function

  impure elemental function mp_min_diim(a,b,c,d) result(y)
    real(c_double), intent(in) :: a
    integer, intent(in) :: b
    integer, intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
      value=mp_from_any(d)
      if (value<y) y=value
  end function

  impure elemental function mp_min_didm(a,b,c,d) result(y)
    real(c_double), intent(in) :: a
    integer, intent(in) :: b
    real(c_double), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
      value=mp_from_any(d)
      if (value<y) y=value
  end function

  impure elemental function mp_min_dism(a,b,c,d) result(y)
    real(c_double), intent(in) :: a
    integer, intent(in) :: b
    real(c_float), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
      value=mp_from_any(d)
      if (value<y) y=value
  end function

  impure elemental function mp_min_ddim(a,b,c,d) result(y)
    real(c_double), intent(in) :: a
    real(c_double), intent(in) :: b
    integer, intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
      value=mp_from_any(d)
      if (value<y) y=value
  end function

  impure elemental function mp_min_dddm(a,b,c,d) result(y)
    real(c_double), intent(in) :: a
    real(c_double), intent(in) :: b
    real(c_double), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
      value=mp_from_any(d)
      if (value<y) y=value
  end function

  impure elemental function mp_min_ddsm(a,b,c,d) result(y)
    real(c_double), intent(in) :: a
    real(c_double), intent(in) :: b
    real(c_float), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
      value=mp_from_any(d)
      if (value<y) y=value
  end function

  impure elemental function mp_min_dsim(a,b,c,d) result(y)
    real(c_double), intent(in) :: a
    real(c_float), intent(in) :: b
    integer, intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
      value=mp_from_any(d)
      if (value<y) y=value
  end function

  impure elemental function mp_min_dsdm(a,b,c,d) result(y)
    real(c_double), intent(in) :: a
    real(c_float), intent(in) :: b
    real(c_double), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
      value=mp_from_any(d)
      if (value<y) y=value
  end function

  impure elemental function mp_min_dssm(a,b,c,d) result(y)
    real(c_double), intent(in) :: a
    real(c_float), intent(in) :: b
    real(c_float), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
      value=mp_from_any(d)
      if (value<y) y=value
  end function

  impure elemental function mp_min_siim(a,b,c,d) result(y)
    real(c_float), intent(in) :: a
    integer, intent(in) :: b
    integer, intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
      value=mp_from_any(d)
      if (value<y) y=value
  end function

  impure elemental function mp_min_sidm(a,b,c,d) result(y)
    real(c_float), intent(in) :: a
    integer, intent(in) :: b
    real(c_double), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
      value=mp_from_any(d)
      if (value<y) y=value
  end function

  impure elemental function mp_min_sism(a,b,c,d) result(y)
    real(c_float), intent(in) :: a
    integer, intent(in) :: b
    real(c_float), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
      value=mp_from_any(d)
      if (value<y) y=value
  end function

  impure elemental function mp_min_sdim(a,b,c,d) result(y)
    real(c_float), intent(in) :: a
    real(c_double), intent(in) :: b
    integer, intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
      value=mp_from_any(d)
      if (value<y) y=value
  end function

  impure elemental function mp_min_sddm(a,b,c,d) result(y)
    real(c_float), intent(in) :: a
    real(c_double), intent(in) :: b
    real(c_double), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
      value=mp_from_any(d)
      if (value<y) y=value
  end function

  impure elemental function mp_min_sdsm(a,b,c,d) result(y)
    real(c_float), intent(in) :: a
    real(c_double), intent(in) :: b
    real(c_float), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
      value=mp_from_any(d)
      if (value<y) y=value
  end function

  impure elemental function mp_min_ssim(a,b,c,d) result(y)
    real(c_float), intent(in) :: a
    real(c_float), intent(in) :: b
    integer, intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
      value=mp_from_any(d)
      if (value<y) y=value
  end function

  impure elemental function mp_min_ssdm(a,b,c,d) result(y)
    real(c_float), intent(in) :: a
    real(c_float), intent(in) :: b
    real(c_double), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
      value=mp_from_any(d)
      if (value<y) y=value
  end function

  impure elemental function mp_min_sssm(a,b,c,d) result(y)
    real(c_float), intent(in) :: a
    real(c_float), intent(in) :: b
    real(c_float), intent(in) :: c
    type(mp_real), intent(in) :: d
    type(mp_real) :: y,value
    y=mp_from_any(a)
      value=mp_from_any(b)
      if (value<y) y=value
      value=mp_from_any(c)
      if (value<y) y=value
      value=mp_from_any(d)
      if (value<y) y=value
  end function

  function mp_sum(a) result(y)
    type(mp_real), intent(in) :: a(:)
    type(mp_real) :: y
    integer :: i
    y=0
    do i=1,size(a)
      y=y+a(i)
    end do
  end function

  function mp_maxval(a) result(y)
    type(mp_real), intent(in) :: a(:)
    type(mp_real) :: y
    integer :: i
    y=-huge(mpreal(0))
    do i=1,size(a)
      if (a(i)>y) y=a(i)
    end do
  end function

  function mp_minval(a) result(y)
    type(mp_real), intent(in) :: a(:)
    type(mp_real) :: y
    integer :: i
    y=huge(mpreal(0))
    do i=1,size(a)
      if (a(i)<y) y=a(i)
    end do
  end function

end module spheroidal_mpfr_real
