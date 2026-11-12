#ifndef SPHEROIDAL_MPFR_ARITHMETIC_H
#define SPHEROIDAL_MPFR_ARITHMETIC_H
#include <stdint.h>
#include <mpfr.h>
#ifdef __cplusplus
extern "C" {
#endif
int64_t sw_mpfr_allocate_words(int bits);
void sw_mpfr_init(int64_t *buffer, int bits);
void sw_mpfr_load(mpfr_ptr value, const int64_t *buffer);
void sw_mpfr_store(int64_t *buffer, mpfr_srcptr value);
void sw_mpfr_set_integer(int64_t *out, int64_t value, int rounding);
void sw_mpfr_set_double(int64_t *out, double value, int rounding);
int sw_mpfr_set_string(int64_t *out, const char *value, int rounding);
void sw_mpfr_binary(int operation, int64_t *out, const int64_t *a, const int64_t *b, int rounding);
void sw_mpfr_unary(int operation, int64_t *out, const int64_t *a, int rounding);
void sw_mpfr_power_integer(int64_t *out, const int64_t *a, int64_t exponent, int rounding);
int sw_mpfr_compare(const int64_t *a, const int64_t *b);
int sw_mpfr_is_finite(const int64_t *value);
int64_t sw_mpfr_to_integer(const int64_t *value, int rounding);
double sw_mpfr_to_double(const int64_t *value, int rounding);
void sw_mpfr_constant(int operation, int64_t *out, int rounding);
int sw_mpfr_format(const int64_t *value, char *text, int capacity);
void sw_mpc_binary(int64_t *re, int64_t *im, const int64_t *ar, const int64_t *ai,
                   const int64_t *br, const int64_t *bi, int operation, int rounding);
void sw_mpc_unary(int64_t *re, int64_t *im, const int64_t *ar, const int64_t *ai,
                  int operation, int rounding);
void sw_mpc_power_si(int64_t *re, int64_t *im, const int64_t *ar, const int64_t *ai,
                     long exponent, int rounding);
void sw_mpc_abs(int64_t *result, const int64_t *ar, const int64_t *ai, int rounding);
#ifdef __cplusplus
}
#endif
#endif
