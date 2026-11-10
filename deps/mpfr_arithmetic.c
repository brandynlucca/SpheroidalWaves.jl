#include "mpfr_arithmetic.h"
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int64_t sw_mpfr_allocate_words(int bits) {
    return 3 + (int64_t)((mpfr_custom_get_size(bits)+sizeof(int64_t)-1)/sizeof(int64_t));
}

void sw_mpfr_load(mpfr_ptr value, const int64_t *buffer) {
    mpfr_custom_init_set(value,(int)buffer[1],(mpfr_exp_t)buffer[2],
                         (mpfr_prec_t)buffer[0],(void*)(buffer+3));
}

static void update_header(int64_t *buffer, mpfr_srcptr value) {
    buffer[1]=(int64_t)mpfr_custom_get_kind(value);
    buffer[2]=(int64_t)mpfr_custom_get_exp(value);
}

void sw_mpfr_init(int64_t *buffer, int bits) {
    buffer[0]=bits;
    buffer[1]=MPFR_ZERO_KIND;
    buffer[2]=0;
    mpfr_custom_init(buffer+3,bits);
}

void sw_mpfr_store(int64_t *buffer, mpfr_srcptr value) {
    mpfr_t destination;
    sw_mpfr_load(destination,buffer);
    mpfr_set(destination,value,MPFR_RNDN);
    update_header(buffer,destination);
}

void sw_mpfr_set_integer(int64_t *out, int64_t value, int rounding) {
    mpfr_t result;
    sw_mpfr_load(result,out);
    mpfr_set_sj(result,(intmax_t)value,(mpfr_rnd_t)rounding);
    update_header(out,result);
}

void sw_mpfr_set_double(int64_t *out, double value, int rounding) {
    mpfr_t result;
    sw_mpfr_load(result,out);
    mpfr_set_d(result,value,(mpfr_rnd_t)rounding);
    update_header(out,result);
}

int sw_mpfr_set_string(int64_t *out, const char *value, int rounding) {
    mpfr_t result;
    sw_mpfr_load(result,out);
    int status;
    if (strchr(value,'/') != NULL) {
        mpq_t rational;
        mpq_init(rational);
        status=mpq_set_str(rational,value,0);
        if (status==0 && mpz_sgn(mpq_denref(rational))!=0) {
            mpq_canonicalize(rational);
            mpfr_set_q(result,rational,(mpfr_rnd_t)rounding);
        } else {
            status=-1;
            mpfr_set_nan(result);
        }
        mpq_clear(rational);
    } else {
        status=mpfr_set_str(result,value,0,(mpfr_rnd_t)rounding);
    }
    update_header(out,result);
    return status;
}

void sw_mpfr_binary(int operation, int64_t *out, const int64_t *a, const int64_t *b, int rounding) {
    mpfr_t result,left,right;
    sw_mpfr_load(result,out);
    sw_mpfr_load(left,a);
    sw_mpfr_load(right,b);
    mpfr_rnd_t mode=(mpfr_rnd_t)rounding;
    switch(operation) {
        case 1: mpfr_add(result,left,right,mode); break;
        case 2: mpfr_sub(result,left,right,mode); break;
        case 3: mpfr_mul(result,left,right,mode); break;
        case 4: mpfr_div(result,left,right,mode); break;
        case 5: mpfr_pow(result,left,right,mode); break;
        case 6: mpfr_atan2(result,left,right,mode); break;
        case 7: mpfr_fmod(result,left,right,mode); break;
        case 8: mpfr_copysign(result,left,right,mode); break;
        case 9: mpfr_hypot(result,left,right,mode); break;
        default: mpfr_set_nan(result); break;
    }
    update_header(out,result);
}

void sw_mpfr_unary(int operation, int64_t *out, const int64_t *a, int rounding) {
    mpfr_t result,value;
    sw_mpfr_load(result,out);
    sw_mpfr_load(value,a);
    mpfr_rnd_t mode=(mpfr_rnd_t)rounding;
    switch(operation) {
        case 1: mpfr_neg(result,value,mode); break;
        case 2: mpfr_abs(result,value,mode); break;
        case 3: mpfr_sqrt(result,value,mode); break;
        case 4: mpfr_exp(result,value,mode); break;
        case 5: mpfr_log(result,value,mode); break;
        case 6: mpfr_log10(result,value,mode); break;
        case 7: mpfr_sin(result,value,mode); break;
        case 8: mpfr_cos(result,value,mode); break;
        case 9: mpfr_sinh(result,value,mode); break;
        case 10: mpfr_cosh(result,value,mode); break;
        case 11: mpfr_atan(result,value,mode); break;
        case 12: mpfr_asin(result,value,mode); break;
        case 13: mpfr_acos(result,value,mode); break;
        case 14: mpfr_tan(result,value,mode); break;
        case 15: mpfr_tanh(result,value,mode); break;
        case 16: mpfr_set(result,value,mode); break;
        default: mpfr_set_nan(result); break;
    }
    update_header(out,result);
}

void sw_mpfr_power_integer(int64_t *out, const int64_t *a, int64_t exponent, int rounding) {
    mpfr_t result,value;
    sw_mpfr_load(result,out);
    sw_mpfr_load(value,a);
    mpfr_pow_sj(result,value,(intmax_t)exponent,(mpfr_rnd_t)rounding);
    update_header(out,result);
}

int sw_mpfr_compare(const int64_t *a, const int64_t *b) {
    mpfr_t left,right;
    sw_mpfr_load(left,a);
    sw_mpfr_load(right,b);
    if (mpfr_nan_p(left) || mpfr_nan_p(right)) return 2;
    return mpfr_less_p(left,right) ? -1 : mpfr_greater_p(left,right) ? 1 : 0;
}

int sw_mpfr_is_finite(const int64_t *value) {
    mpfr_t number;
    sw_mpfr_load(number,value);
    return mpfr_number_p(number);
}

int64_t sw_mpfr_to_integer(const int64_t *value, int rounding) {
    mpfr_t number;
    sw_mpfr_load(number,value);
    if (rounding == -1) {
        mpfr_t nearest;
        mpfr_init2(nearest,mpfr_get_prec(number));
        mpfr_round(nearest,number);
        int64_t result=(int64_t)mpfr_get_sj(nearest,MPFR_RNDZ);
        mpfr_clear(nearest);
        return result;
    }
    return (int64_t)mpfr_get_sj(number,(mpfr_rnd_t)rounding);
}

double sw_mpfr_to_double(const int64_t *value, int rounding) {
    mpfr_t number;
    sw_mpfr_load(number,value);
    return mpfr_get_d(number,(mpfr_rnd_t)rounding);
}

void sw_mpfr_constant(int operation, int64_t *out, int rounding) {
    mpfr_t result;
    sw_mpfr_load(result,out);
    switch(operation) {
        case 1: mpfr_set_ui_2exp(result,1,1-out[0],(mpfr_rnd_t)rounding); break;
        case 2: mpfr_set_inf(result,1); mpfr_nextbelow(result); break;
        case 3: mpfr_set_zero(result,1); mpfr_nextabove(result); break;
        default: mpfr_set_nan(result); break;
    }
    update_header(out,result);
}

int sw_mpfr_format(const int64_t *value, char *text, int capacity) {
    mpfr_t number;
    sw_mpfr_load(number,value);
    if (mpfr_nan_p(number)) return snprintf(text,capacity,"NaN");
    if (mpfr_inf_p(number)) return snprintf(text,capacity,mpfr_signbit(number)?"-Inf":"Inf");
    if (mpfr_zero_p(number)) return snprintf(text,capacity,mpfr_signbit(number)?"-0":"0");
    mpfr_exp_t exponent;
    char *digits=mpfr_get_str(NULL,&exponent,16,0,number,MPFR_RNDN);
    int negative=digits[0]=='-';
    int length=snprintf(text,capacity,"%s0x0.%sp%" PRIdMAX,negative?"-":"",digits+negative,4*(intmax_t)exponent);
    mpfr_free_str(digits);
    return length;
}
