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
#ifdef __cplusplus
}
#endif
#endif
