#include "precision.inc"

#ifdef SWF_PRECISION_MPFR
module spheroidal_precision_mpfr
#else
#ifdef SWF_PRECISION_QUAD
module spheroidal_precision_quad
#else
module spheroidal_precision_double
#endif
#endif
  use, intrinsic :: iso_c_binding, only: c_float,c_double,c_int64_t
#ifdef SWF_PRECISION_MPFR
  use spheroidal_mpfr_complex
#endif
  implicit none
#ifdef SWF_PRECISION_MPFR
  ! This branch selector preserves quad policies without requiring a native real(16) type.
  integer, parameter :: sw_kind = 16
#else
  integer, parameter :: sw_kind = SW_REAL_KIND
#endif
  interface sw_real
    module procedure sw_from_integer,sw_from_integer64,sw_from_real,sw_from_complex,sw_from_float
#ifndef SWF_PRECISION_DOUBLE
    module procedure sw_from_double
#endif
  end interface
  interface sw_complex
    module procedure sw_complex_from_parts
  end interface
  private :: sw_from_integer,sw_from_integer64,sw_from_real,sw_from_complex,sw_from_float
  private :: sw_complex_from_parts
#ifndef SWF_PRECISION_DOUBLE
  private :: sw_from_double
#endif
  interface sw_release
    module procedure sw_release_real,sw_release_complex
  end interface
  private :: sw_release_real,sw_release_complex
contains
  impure elemental function sw_from_integer(value,kind) result(converted)
    integer, intent(in) :: value
    integer, optional, intent(in) :: kind
    SW_REAL :: converted
#ifdef SWF_PRECISION_MPFR
    converted = mpreal(value)
#else
    converted = real(value,sw_kind)
#endif
  end function

  impure elemental function sw_from_integer64(value,kind) result(converted)
    integer(c_int64_t), intent(in) :: value
    integer, optional, intent(in) :: kind
    SW_REAL :: converted
#ifdef SWF_PRECISION_MPFR
    converted = mpreal(value)
#else
    converted = real(value,sw_kind)
#endif
  end function

  impure elemental function sw_from_real(value,kind) result(converted)
    SW_REAL, intent(in) :: value
    integer, optional, intent(in) :: kind
    SW_REAL :: converted
#ifdef SWF_PRECISION_MPFR
    converted = mpreal(value)
#else
    converted = value
#endif
  end function

  impure elemental function sw_from_complex(value,kind) result(converted)
    SW_COMPLEX, intent(in) :: value
    integer, optional, intent(in) :: kind
    SW_REAL :: converted
#ifdef SWF_PRECISION_MPFR
    converted = mpreal(value)
#else
    converted = real(value,sw_kind)
#endif
  end function

  impure elemental function sw_from_float(value,kind) result(converted)
    real(c_float), intent(in) :: value
    integer, optional, intent(in) :: kind
    SW_REAL :: converted
#ifdef SWF_PRECISION_MPFR
    converted = mpreal(value)
#else
    converted = real(value,sw_kind)
#endif
  end function

#ifndef SWF_PRECISION_DOUBLE
  impure elemental function sw_from_double(value,kind) result(converted)
    real(c_double), intent(in) :: value
    integer, optional, intent(in) :: kind
    SW_REAL :: converted
#ifdef SWF_PRECISION_MPFR
    converted = mpreal(value)
#else
    converted = real(value,sw_kind)
#endif
  end function
#endif

  impure elemental function sw_complex_from_parts(re,im,kind) result(converted)
    SW_REAL, intent(in) :: re
    SW_REAL, optional, intent(in) :: im
    integer, optional, intent(in) :: kind
    SW_COMPLEX :: converted
#ifdef SWF_PRECISION_MPFR
    if (present(im)) then
      converted = mpcmplx(re,im)
    else
      converted = mpcmplx(re)
    end if
#else
    converted = cmplx(re,kind=sw_kind)
    if (present(im)) converted = cmplx(re,im,sw_kind)
#endif
  end function

  integer(c_int64_t) function sw_real_bytes() result(bytes)
    SW_REAL :: element
    bytes = (storage_size(element,kind=c_int64_t)+7_c_int64_t)/8_c_int64_t
#ifdef SWF_PRECISION_MPFR
    bytes = bytes+88_c_int64_t+8_c_int64_t*((int(mpfr_bits,c_int64_t)+63_c_int64_t)/64_c_int64_t)
#endif
  end function

  impure elemental logical function sw_cache_equal(a,b) result(equal)
    SW_REAL, intent(in) :: a,b
#ifdef SWF_PRECISION_MPFR
    equal = a==b
#else
    SW_REAL :: scale,tolerance
    scale = max(1.0_sw_kind,max(abs(a),abs(b)))
    tolerance = 64.0_sw_kind*epsilon(1.0_sw_kind)*scale
    equal = abs(a-b)<=tolerance
#endif
  end function

  impure elemental subroutine sw_release_real(value)
    SW_REAL, intent(inout) :: value
#ifdef SWF_PRECISION_MPFR
    if (allocated(value%data)) deallocate(value%data)
#else
    value = 0.0_sw_kind
#endif
  end subroutine

  impure elemental subroutine sw_release_complex(value)
    SW_COMPLEX, intent(inout) :: value
#ifdef SWF_PRECISION_MPFR
    call sw_release_real(value%re)
    call sw_release_real(value%im)
#else
    value = cmplx(0.0_sw_kind,0.0_sw_kind,sw_kind)
#endif
  end subroutine

  logical function sw_account_workspace(total,element_bytes,copies,dimensions,limit) result(allowed)
    integer(c_int64_t), intent(inout) :: total
    integer(c_int64_t), intent(in) :: element_bytes,copies,dimensions(:),limit
    integer(c_int64_t) :: available,bytes,dimension
    integer :: index
    allowed = .false.
    if (total<0.or.limit<0.or.total>limit.or.element_bytes<=0.or.copies<=0) return
    available = limit-total
    if (copies>available/element_bytes) return
    bytes = copies*element_bytes
    do index=1,size(dimensions)
      dimension = dimensions(index)
      if (dimension<0) return
      if (dimension==0) then
        allowed = .true.
        return
      end if
      if (dimension>available/bytes) return
      bytes = bytes*dimension
    end do
    total = total+bytes
    allowed = .true.
  end function
end module
