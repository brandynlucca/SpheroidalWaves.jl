#include "mpfr_arithmetic.h"
#include <mpc.h>

static mpc_rnd_t sw_mpc_rounding(int rounding) {
    return MPC_RND((mpfr_rnd_t)rounding, (mpfr_rnd_t)rounding);
}

static void sw_mpc_load(mpc_ptr value, const int64_t *re, const int64_t *im) {
    sw_mpfr_load(mpc_realref(value), re);
    sw_mpfr_load(mpc_imagref(value), im);
}

static void sw_mpc_store(int64_t *re, int64_t *im, mpc_srcptr value) {
    sw_mpfr_store(re, mpc_realref(value));
    sw_mpfr_store(im, mpc_imagref(value));
}

void sw_mpc_binary(int64_t *re, int64_t *im, const int64_t *ar, const int64_t *ai,
                   const int64_t *br, const int64_t *bi, int operation, int rounding) {
    mpc_t a, b, result;
    sw_mpc_load(a, ar, ai);
    sw_mpc_load(b, br, bi);
    mpc_init3(result, (mpfr_prec_t)re[0], (mpfr_prec_t)im[0]);
    mpc_rnd_t mode = sw_mpc_rounding(rounding);
    switch (operation) {
        case 0: mpc_add(result, a, b, mode); break;
        case 1: mpc_sub(result, a, b, mode); break;
        case 2: mpc_mul(result, a, b, mode); break;
        case 3: mpc_div(result, a, b, mode); break;
        case 4: mpc_pow(result, a, b, mode); break;
        default: mpfr_set_nan(mpc_realref(result)); mpfr_set_nan(mpc_imagref(result)); break;
    }
    sw_mpc_store(re, im, result);
    mpc_clear(result);
}

void sw_mpc_unary(int64_t *re, int64_t *im, const int64_t *ar, const int64_t *ai,
                  int operation, int rounding) {
    mpc_t a, result;
    sw_mpc_load(a, ar, ai);
    mpc_init3(result, (mpfr_prec_t)re[0], (mpfr_prec_t)im[0]);
    mpc_rnd_t mode = sw_mpc_rounding(rounding);
    switch (operation) {
        case 0: mpc_neg(result, a, mode); break;
        case 1: mpc_conj(result, a, mode); break;
        case 2: mpc_sqrt(result, a, mode); break;
        case 3: mpc_exp(result, a, mode); break;
        case 4: mpc_log(result, a, mode); break;
        case 5: mpc_sin(result, a, mode); break;
        case 6: mpc_cos(result, a, mode); break;
        case 7: mpc_sinh(result, a, mode); break;
        case 8: mpc_cosh(result, a, mode); break;
        default: mpfr_set_nan(mpc_realref(result)); mpfr_set_nan(mpc_imagref(result)); break;
    }
    sw_mpc_store(re, im, result);
    mpc_clear(result);
}

void sw_mpc_power_si(int64_t *re, int64_t *im, const int64_t *ar, const int64_t *ai,
                     long exponent, int rounding) {
    mpc_t a, result;
    sw_mpc_load(a, ar, ai);
    mpc_init3(result, (mpfr_prec_t)re[0], (mpfr_prec_t)im[0]);
    mpc_pow_si(result, a, exponent, sw_mpc_rounding(rounding));
    sw_mpc_store(re, im, result);
    mpc_clear(result);
}

void sw_mpc_abs(int64_t *result, const int64_t *ar, const int64_t *ai, int rounding) {
    mpc_t a;
    mpfr_t value;
    sw_mpc_load(a, ar, ai);
    mpfr_init2(value, (mpfr_prec_t)result[0]);
    mpc_abs(value, a, (mpfr_rnd_t)rounding);
    sw_mpfr_store(result, value);
    mpfr_clear(value);
}
