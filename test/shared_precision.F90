#include "precision.inc"

module shared_precision_example
  use, intrinsic :: iso_c_binding, only: c_int64_t
#ifdef SWF_PRECISION_MPFR
  use spheroidal_mpfr_complex
#endif
  implicit none
#ifndef SWF_PRECISION_MPFR
  integer, parameter :: sw_kind = SW_REAL_KIND
#endif
  interface to_real
    module procedure integer_to_real, integer64_to_real, real_to_real, complex_to_real
  end interface
contains
  impure elemental function integer_to_real(value) result(converted)
    integer, intent(in) :: value
    SW_REAL :: converted
#ifdef SWF_PRECISION_MPFR
    converted = mpreal(value)
#else
    converted = real(value,kind=sw_kind)
#endif
  end function

  impure elemental function integer64_to_real(value) result(converted)
    integer(c_int64_t), intent(in) :: value
    SW_REAL :: converted
#ifdef SWF_PRECISION_MPFR
    converted = mpreal(value)
#else
    converted = real(value,kind=sw_kind)
#endif
  end function

  impure elemental function real_to_real(value) result(converted)
    SW_REAL, intent(in) :: value
    SW_REAL :: converted
#ifdef SWF_PRECISION_MPFR
    converted = mpreal(value)
#else
    converted = real(value,kind=sw_kind)
#endif
  end function

  impure elemental function complex_to_real(value) result(converted)
    SW_COMPLEX, intent(in) :: value
    SW_REAL :: converted
#ifdef SWF_PRECISION_MPFR
    converted = mpreal(value)
#else
    converted = real(value,kind=sw_kind)
#endif
  end function

  impure elemental function to_complex(re,im) result(converted)
    SW_REAL, intent(in) :: re,im
    SW_COMPLEX :: converted
#ifdef SWF_PRECISION_MPFR
    converted = mpcmplx(re,im)
#else
    converted = cmplx(re,im,kind=sw_kind)
#endif
  end function

  subroutine evaluate_shared(x,dyadic,negative,z,shifted)
    SW_REAL, intent(in) :: x
    SW_REAL, intent(out) :: dyadic,negative,shifted
    SW_COMPLEX, intent(out) :: z
#ifdef SWF_PRECISION_MPFR
    SW_REAL :: offset
    offset = SW_LITERAL(0.125_sw_kind,"0.125")
#else
    SW_REAL, parameter :: offset = SW_LITERAL(0.125_sw_kind,"0.125")
#endif
    dyadic = SW_LITERAL(1.00000000000000000000000000000000000000293873587705571876992184134305561419454666389193021880377187926569604314863681793212890625_sw_kind,"1.00000000000000000000000000000000000000293873587705571876992184134305561419454666389193021880377187926569604314863681793212890625")
    ! The sign belongs inside the decimal text to preserve directed rounding.
    negative = SW_LITERAL(-0.1_sw_kind,"-0.1")
    z = to_complex(dyadic,negative)
    shifted = (x+offset)/2
  end subroutine
end module

program shared_precision
  use shared_precision_example
  implicit none
  SW_REAL :: x,dyadic,negative,shifted,values(3),copies(3),wide
  SW_COMPLEX :: z,copy_z
  integer(c_int64_t), parameter :: large_integer = 9007199254740993_c_int64_t
#ifdef SWF_PRECISION_MPFR
  integer :: bits_index,mode,saved_bits,saved_rounding
  integer, parameter :: precisions(4) = [160,256,384,160]
  SW_REAL :: retained,rounded,copied
  saved_bits = mpfr_bits
  saved_rounding = mpfr_rounding
  do bits_index=1,size(precisions)
    mpfr_bits = precisions(bits_index)
    do mode=0,4
      mpfr_rounding = mode
      call check_shared()
    end do
  end do
  mpfr_bits = 256
  mpfr_rounding = 0
  retained = mpreal("1/10")
  mpfr_bits = 160
  copied = retained
  rounded = to_real(retained)
  if (copied%data(1)/=256.or.rounded%data(1)/=160) error stop "copy and conversion precision"
  if (rounded/=mpreal("1/10")) error stop "conversion rounds at active precision"
  retained = 0
  if (copied==0.or.rounded==0) error stop "conversion ownership"
  mpfr_bits = saved_bits
  mpfr_rounding = saved_rounding
#else
  call check_shared()
#endif
  print *, "Shared precision declarations, literals, conversions and initialization passed"
contains
  subroutine check_shared()
    x = to_real(3)
    call evaluate_shared(x,dyadic,negative,z,shifted)
    if (shifted/=SW_LITERAL(1.5625_sw_kind,"1.5625")) error stop "initialized offset"
    if (to_real(z)/=dyadic.or.aimag(z)/=negative) error stop "complex components"
    if (to_real(1/3)/=0.or.to_real(1)/3==0) error stop "integer division"
    if (to_real(1)/2/=SW_LITERAL(0.5_sw_kind,"0.5")) error stop "real division"
    wide = to_real(large_integer)
    if (wide/=SW_LITERAL(9007199254740993.0_sw_kind,"9007199254740993")) error stop "integer64 conversion"
    values = [dyadic,negative,x]
    copies = values
    values = 0
    if (any(copies/=[dyadic,negative,x])) error stop "array ownership"
    copy_z = z
    z = to_complex(to_real(0),to_real(0))
    if (to_real(copy_z)/=dyadic.or.aimag(copy_z)/=negative) error stop "complex ownership"
#ifdef SWF_PRECISION_MPFR
    if (dyadic/=to_real(1)+to_real(2)**(-128).or.dyadic==1) error stop "literal beyond quad precision"
    if (negative/=mpreal("-1/10")) error stop "signed decimal rounding"
    if (mode==2.or.mode==3) then
      if (negative==-mpreal("0.1")) error stop "signed literal differs from unary negation"
    end if
    if (dyadic%data(1)/=mpfr_bits.or.shifted%data(1)/=mpfr_bits) error stop "runtime initialization precision"
    if (copy_z%re%data(1)/=mpfr_bits.or.copy_z%im%data(1)/=mpfr_bits) error stop "complex precision"
    if (mp_to_integer64(wide)/=large_integer) error stop "exact integer64 conversion"
    if (mpfr_status/=0) error stop "arithmetic status"
#else
    if (dyadic/=1) error stop "native rounding"
    if (kind(shifted)/=sw_kind.or.kind(copy_z)/=sw_kind) error stop "native result kind"
#ifdef SWF_PRECISION_QUAD
    if (int(wide,kind=c_int64_t)/=large_integer) error stop "quad integer64 conversion"
#else
    if (int(wide,kind=c_int64_t)/=large_integer-1) error stop "double integer64 rounding"
#endif
#endif
  end subroutine
end program
