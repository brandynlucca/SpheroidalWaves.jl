#ifndef SPHEROIDAL_PRECISION_INC
#define SPHEROIDAL_PRECISION_INC

#if (defined(SWF_PRECISION_DOUBLE) + defined(SWF_PRECISION_QUAD) + defined(SWF_PRECISION_MPFR)) != 1
#error Select exactly one SWF_PRECISION_DOUBLE, SWF_PRECISION_QUAD or SWF_PRECISION_MPFR
#endif

! Supply both literal spellings explicitly; GNU traditional preprocessing cannot stringify tokens.
#ifdef SWF_PRECISION_MPFR
#define SW_PRECISION_MODULE spheroidal_precision_mpfr
#define SW_PROLATE_PARAMETERS prolate_parameters_mpfr
#define SW_PROLATE_MODULE prolate_swf_mpfr
#define SW_OBLATE_MODULE oblate_swf_mpfr
#define SW_CPROLATE_MODULE complex_prolate_swf_mpfr
#define SW_REAL type(mp_real)
#define SW_COMPLEX type(mp_complex)
#define SW_REAL1 type(mp_real)
#define SW_COMPLEX1 type(mp_complex)
#define SW_LITERAL(native,decimal) mpreal(decimal)
#define SW_CONVERT_REAL mpreal
#define SW_CONVERT_COMPLEX mpcmplx
#define SW_ZERO_INIT
#else
#define SW_PROLATE_PARAMETERS prolate_parameters
#define SW_PROLATE_MODULE prolate_swf
#define SW_OBLATE_MODULE oblate_swf
#define SW_CPROLATE_MODULE complex_prolate_swf
#define SW_REAL real(sw_kind)
#define SW_COMPLEX complex(sw_kind)
#define SW_REAL1 real(knd1)
#define SW_COMPLEX1 complex(knd1)
#define SW_LITERAL(native,decimal) native
#define SW_CONVERT_REAL real
#define SW_CONVERT_COMPLEX cmplx
#define SW_ZERO_INIT = 0.0_sw_kind
#ifdef SWF_PRECISION_QUAD
#define SW_PRECISION_MODULE spheroidal_precision_quad
#define SW_REAL_KIND selected_real_kind(33)
#else
#define SW_PRECISION_MODULE spheroidal_precision_double
#define SW_REAL_KIND selected_real_kind(8)
#endif
#endif

#endif
