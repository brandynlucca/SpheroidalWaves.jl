#include "mpfr_arithmetic.h"
#include <mpc.h>
#include <stdio.h>
#include <stdlib.h>

static int64_t *buffer(int bits) {
    size_t words = (size_t)sw_mpfr_allocate_words(bits);
    int64_t *value = calloc(words+2, sizeof(int64_t));
    if (value == NULL) abort();
    value[0] = value[words+1] = INT64_C(0x38e5a72b46c019df);
    sw_mpfr_init(value+1, bits);
    return value+1;
}

static void release_buffer(int64_t *value, int bits) {
    if (value[-1] != INT64_C(0x38e5a72b46c019df) ||
        value[sw_mpfr_allocate_words(bits)] != INT64_C(0x38e5a72b46c019df)) {
        fputs("MPC operation wrote beyond its supplied buffer\n", stderr);
        abort();
    }
    free(value-1);
}

static void check(const int64_t *actual, mpfr_srcptr expected) {
    mpfr_t value;
    sw_mpfr_load(value, actual);
    int equal = mpfr_nan_p(expected) ? mpfr_nan_p(value) :
        mpfr_equal_p(value, expected) && mpfr_signbit(value) == mpfr_signbit(expected);
    if (!equal || mpfr_get_prec(value) != mpfr_get_prec(expected)) {
        fputs("MPC buffer result differs from an independently allocated result\n", stderr);
        abort();
    }
}

int main(void) {
    const int precisions[] = {2, 53, 113, 160, 256, 512, 2048};
    const char *values[][2] = {
        {"0", "-0"}, {"-0", "0"}, {"1.125", "-0.375"}, {"-3.25", "2.5"},
        {"1e10", "1e-10"}, {"1e-10", "-1e10"}, {"Inf", "1"},
        {"-Inf", "-Inf"}, {"NaN", "0"}, {"0", "NaN"}
    };
    int (*binary[])(mpc_ptr, mpc_srcptr, mpc_srcptr, mpc_rnd_t) =
        {mpc_add, mpc_sub, mpc_mul, mpc_div, mpc_pow};
    int (*unary[])(mpc_ptr, mpc_srcptr, mpc_rnd_t) =
        {mpc_neg, mpc_conj, mpc_sqrt, mpc_exp, mpc_log, mpc_sin, mpc_cos, mpc_sinh, mpc_cosh};
    size_t count = 0;
    for (size_t p = 0; p < sizeof(precisions)/sizeof(*precisions); ++p) {
        int bits = precisions[p];
        int64_t *re = buffer(bits), *im = buffer(bits+17);
        int64_t *ar = buffer(bits+32), *ai = buffer(bits+32);
        int64_t *br = buffer(bits+32), *bi = buffer(bits+32);
        mpc_t a, b, expected;
        mpc_init2(a, bits+32);
        mpc_init2(b, bits+32);
        mpc_init3(expected, bits, bits+17);
        for (int rounding = 0; rounding <= 4; ++rounding) {
            mpc_rnd_t mode = MPC_RND((mpfr_rnd_t)rounding, (mpfr_rnd_t)rounding);
            for (size_t i = 0; i < sizeof(values)/sizeof(*values); ++i) {
                size_t j = (i+3) % (sizeof(values)/sizeof(*values));
                mpfr_set_str(mpc_realref(a), values[i][0], 10, MPFR_RNDN);
                mpfr_set_str(mpc_imagref(a), values[i][1], 10, MPFR_RNDN);
                mpfr_set_str(mpc_realref(b), values[j][0], 10, MPFR_RNDN);
                mpfr_set_str(mpc_imagref(b), values[j][1], 10, MPFR_RNDN);
                sw_mpfr_store(ar, mpc_realref(a));
                sw_mpfr_store(ai, mpc_imagref(a));
                sw_mpfr_store(br, mpc_realref(b));
                sw_mpfr_store(bi, mpc_imagref(b));
                for (int operation = 0; operation < 5; ++operation) {
                    binary[operation](expected, a, b, mode);
                    sw_mpc_binary(re, im, ar, ai, br, bi, operation, rounding);
                    check(re, mpc_realref(expected));
                    check(im, mpc_imagref(expected));
                    ++count;
                }
                for (int operation = 0; operation < 9; ++operation) {
                    unary[operation](expected, a, mode);
                    sw_mpc_unary(re, im, ar, ai, operation, rounding);
                    check(re, mpc_realref(expected));
                    check(im, mpc_imagref(expected));
                    ++count;
                }
                for (long exponent = -3; exponent <= 3; ++exponent) {
                    mpc_pow_si(expected, a, exponent, mode);
                    sw_mpc_power_si(re, im, ar, ai, exponent, rounding);
                    check(re, mpc_realref(expected));
                    check(im, mpc_imagref(expected));
                    ++count;
                }
                mpc_abs(mpc_realref(expected), a, (mpfr_rnd_t)rounding);
                sw_mpc_abs(re, ar, ai, rounding);
                check(re, mpc_realref(expected));
                check(ar, mpc_realref(a));
                check(ai, mpc_imagref(a));
                check(br, mpc_realref(b));
                check(bi, mpc_imagref(b));
                ++count;
            }
        }
        mpc_clear(a);
        mpc_clear(b);
        mpc_clear(expected);
        release_buffer(re, bits); release_buffer(im, bits+17);
        release_buffer(ar, bits+32); release_buffer(ai, bits+32);
        release_buffer(br, bits+32); release_buffer(bi, bits+32);
    }
    mpfr_free_cache();
    printf("%zu MPC buffer comparisons passed\n", count);
    return 0;
}
