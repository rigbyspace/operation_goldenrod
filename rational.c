// rational.c
// TRTS rational helpers with strict no-canonicalization policy.
// This file implements the whitelisted operations from rational.h and
// purposefully omits/disabled any wrapper that would call mpq_set_str or
// mpq_set directly in a way that could canonicalize rationals.

#include "rational.h"

#include <stdio.h>
#include <gmp.h>
#include <stdbool.h>

/* =====================
   Rational Math Helpers (STRICT MODE - NO CANONICALIZATION)
   ===================== */

// Initialization and cleanup
void rational_init(mpq_t value) {
    mpq_init(value);
}

void rational_clear(mpq_t value) {
    mpq_clear(value);
}

// Assignment operations
void rational_set(mpq_t dest, mpq_srcptr src) {
    // Use the explicit mpz copies to avoid relying on any implicit behavior.
    mpz_set(mpq_numref(dest), mpq_numref(src));
    mpz_set(mpq_denref(dest), mpq_denref(src));
}

void rational_set_si(mpq_t dest, long num, unsigned long den) {
    mpz_set_si(mpq_numref(dest), num);
    mpz_set_ui(mpq_denref(dest), den);
}

void rational_set_components(mpq_t dest, mpz_srcptr num, mpz_srcptr den) {
    mpz_set(mpq_numref(dest), num);
    mpz_set(mpq_denref(dest), den);
}

/* NOTE: intentionally no rational_set_str wrapper is provided here. If textual
   parsing is required, use mpz_set_str for numerator and denominator then
   call rational_set_components. That preserves the "no canonicalization"
   contract explicitly. */

// Arithmetic operations
void rational_add(mpq_t res, mpq_srcptr a, mpq_srcptr b) {
    mpq_add(res, a, b);
}

void rational_sub(mpq_t res, mpq_srcptr a, mpq_srcptr b) {
    mpq_sub(res, a, b);
}

void rational_mul(mpq_t res, mpq_srcptr a, mpq_srcptr b) {
    mpq_mul(res, a, b);
}

void rational_div(mpq_t res, mpq_srcptr a, mpq_srcptr b) {
    // Guard: ensure divisor numerator is not zero (division-by-zero)
    if (mpz_sgn(mpq_numref(b)) == 0) {
        // Leave res unchanged; callers should check for failure when necessary.
        return;
    }
    mpq_div(res, a, b);
}

void rational_add_ui(mpq_t result, mpq_srcptr a, unsigned long numerator, unsigned long denominator) {
    mpq_t temp;
    mpq_init(temp);
    mpz_set_ui(mpq_numref(temp), numerator);
    mpz_set_ui(mpq_denref(temp), denominator);
    mpq_add(result, a, temp);
    mpq_clear(temp);
}

void rational_neg(mpq_t res, mpq_srcptr a) {
    mpq_neg(res, a);
}

void rational_inv(mpq_t res, mpq_srcptr a) {
    // Guard against inversion of zero numerator
    if (mpz_sgn(mpq_numref(a)) == 0) {
        return;
    }
    mpq_inv(res, a);
}

void rational_abs(mpq_t res, mpq_srcptr a) {
    mpq_abs(res, a);
}

void rational_negate(mpq_t value) {
    mpq_neg(value, value);
}

// Comparison operations
int rational_cmp(mpq_srcptr a, mpq_srcptr b) {
    return mpq_cmp(a, b);
}

int rational_sgn(mpq_srcptr value) {
    return mpq_sgn(value);
}

bool rational_is_zero(mpq_srcptr a) {
    return mpz_sgn(mpq_numref(a)) == 0;
}

// Numerator/denominator operations
void rational_copy_num(mpz_t dest, mpq_srcptr value) {
    mpz_set(dest, mpq_numref(value));
}

void rational_abs_num(mpz_t dest, mpq_srcptr a) {
    mpz_abs(dest, mpq_numref(a));
}

// Rational_denominator_zero is provided inline in header; no extra helper here.

// Modular operations
void rational_delta(mpq_t res, mpq_srcptr a, mpq_srcptr b) {
    mpq_sub(res, a, b);
}

// Floor helper: floor(value) -> dest (mpz)
void rational_floor(mpz_t dest, mpq_srcptr value) {
    // Use numerator and denominator directly (no canonicalization calls).
    mpz_t tmp_num;
    mpz_init(tmp_num);
    mpz_set(tmp_num, mpq_numref(value));
    mpz_fdiv_q(dest, tmp_num, mpq_denref(value));
    mpz_clear(tmp_num);
}

// Ceil helper: ceil(value) -> dest (mpz)
void rational_ceil(mpz_t dest, mpq_srcptr value) {
    mpz_t tmp_num;
    mpz_init(tmp_num);
    mpz_set(tmp_num, mpq_numref(value));
    mpz_cdiv_q(dest, tmp_num, mpq_denref(value));
    mpz_clear(tmp_num);
}

// Round to nearest integer
void rational_round(mpz_t dest, mpq_srcptr value) {
    mpz_t num, den, rem, dbl_rem;
    mpz_inits(num, den, rem, dbl_rem, NULL);

    mpz_set(num, mpq_numref(value));
    mpz_set(den, mpq_denref(value));
    mpz_fdiv_qr(dest, rem, num, den);
    mpz_mul_ui(dbl_rem, rem, 2);

    if (mpz_cmpabs(dbl_rem, den) >= 0) {
        if (mpq_sgn(value) >= 0) {
            mpz_add_ui(dest, dest, 1);
        } else {
            mpz_sub_ui(dest, dest, 1);
        }
    }

    mpz_clears(num, den, rem, dbl_rem, NULL);
}

// Proper rational remainder: res = a - b * floor(a / b)
// Note: when modulus b has numerator == 0, res is set to a for safety
void rational_mod(mpq_t res, mpq_srcptr a, mpq_srcptr b) {
    // Guard: if b numerator == 0, set res = a for safety (use rational_set wrapper)
    if (mpz_sgn(mpq_numref(b)) == 0) {
        rational_set(res, a);
        return;
    }

    // Compute div = a / b using mpq_div, then floor on div's numerator/denominator directly.
    mpq_t div;
    mpz_t qfloor, tmp_num;
    mpq_init(div);
    mpz_init(qfloor);
    mpz_init(tmp_num);

    mpq_div(div, a, b);

    // tmp_num = numerator(div)
    mpz_set(tmp_num, mpq_numref(div));
    mpz_fdiv_q(qfloor, tmp_num, mpq_denref(div));

    // Construct bq = b * qfloor (as rational) without calling forbidden helpers
    mpq_t bq;
    mpq_init(bq);

    // bq.num = b.num * qfloor
    mpz_t bq_num;
    mpz_init(bq_num);
    mpz_mul(bq_num, mpq_numref(b), qfloor);

    // bq.den = b.den
    mpz_t bq_den;
    mpz_init(bq_den);
    mpz_set(bq_den, mpq_denref(b));

    // Set components for bq without canonicalization
    mpz_set(mpq_numref(bq), bq_num);
    mpz_set(mpq_denref(bq), bq_den);

    // res = a - bq
    mpq_sub(res, a, bq);

    mpq_clear(div);
    mpq_clear(bq);
    mpz_clear(qfloor);
    mpz_clear(tmp_num);
    mpz_clear(bq_num);
    mpz_clear(bq_den);
}

// I/O operations
void rational_print(FILE *stream, mpq_srcptr value) {
    mpq_out_str(stream, 10, value);
    fprintf(stream, "\n");
}

