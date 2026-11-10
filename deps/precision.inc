#ifndef SPHEROIDAL_PRECISION_INC
#define SPHEROIDAL_PRECISION_INC

#if (defined(SWF_PRECISION_DOUBLE) + defined(SWF_PRECISION_QUAD) + defined(SWF_PRECISION_MPFR)) != 1
#error Select exactly one SWF_PRECISION_DOUBLE, SWF_PRECISION_QUAD or SWF_PRECISION_MPFR
#endif

! Supply both literal spellings explicitly; GNU traditional preprocessing cannot stringify tokens.
#ifdef SWF_PRECISION_MPFR
#define SW_REAL type(mp_real)
#define SW_COMPLEX type(mp_complex)
#define SW_LITERAL(native,decimal) mpreal(decimal)
#else
#define SW_REAL real(sw_kind)
#define SW_COMPLEX complex(sw_kind)
#define SW_LITERAL(native,decimal) native
#ifdef SWF_PRECISION_QUAD
#define SW_REAL_KIND selected_real_kind(33)
#else
#define SW_REAL_KIND selected_real_kind(8)
#endif
#endif

#endif
