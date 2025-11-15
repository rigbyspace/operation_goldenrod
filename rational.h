/* rational.h - TRTS Rational Number Type
 *
 * Implements strictly rational arithmetic with the following invariants:
 * 1. No canonicalization - numerators and denominators are never reduced
 * 2. No GCD operations - fractions remain in their original form
 * 3. Zero numerators force zero denominators (0/0 representation)
 * 4. All arithmetic preserves these properties
 */

#ifndef TRTS_RATIONAL_H
#define TRTS_RATIONAL_H

#include <gmp.h>
#include <stdbool.h>

/* TRTS Rational number represented as separate numerator and denominator.
 * Zero-based counting: a rational with numerator 0 must have denominator 0.
 * This represents the undefined/counting state (0/0).
 */
typedef struct {
    mpz_t num;  /* numerator */
    mpz_t den;  /* denominator */
} Rational;

/* Initialize rational to 0/1 */
void rational_init(Rational *q);

/* Clear rational's GMP resources */
void rational_clear(Rational *q);

/* Copy src to q */
void rational_set(Rational *q, const Rational *src);

/* Set q = n/d with signed numerator and unsigned denominator
 * If n == 0, forces denominator to 0 (0/0 invariant) */
void rational_set_si(Rational *q, long n, unsigned long d);

/* Set q from existing mpz_t components (copies values)
 * Enforces 0/0 invariant if num is zero */
void rational_set_components(Rational *q, const mpz_t num, const mpz_t den);

/* Arithmetic operations - all preserve non-canonicalized form */
void rational_add(Rational *r, const Rational *a, const Rational *b);
void rational_sub(Rational *r, const Rational *a, const Rational *b);
void rational_mul(Rational *r, const Rational *a, const Rational *b);

/* Division: r = a/b. Returns false if b.num == 0 */
bool rational_div(Rational *r, const Rational *a, const Rational *b);

/* Negate q (flip sign of numerator) */
void rational_negate(Rational *q);

/* Absolute value (both num and den become non-negative) */
void rational_abs(Rational *q);

/* Modular operations */
void rational_mod(Rational *r, const Rational *a, const Rational *b);
void rational_delta(Rational *r, const Rational *a, const Rational *b);
void rational_floor(Rational *r, const Rational *q);
void rational_ceil(Rational *r, const Rational *q);
void rational_round(Rational *r, const Rational *q);

/* Copy absolute value of numerator to dest */
void rational_abs_num(mpz_t dest, const Rational *q);

/* Comparison: returns <0 if a<b, 0 if a==b, >0 if a>b
 * Compares using cross-multiplication without reduction */
int rational_cmp(const Rational *a, const Rational *b);

/* Sign of rational (-1, 0, or 1) */
int rational_sgn(const Rational *q);

/* Test if numerator is zero (which implies denominator is also zero) */
bool rational_is_zero(const Rational *q);

/* Test if denominator is zero (undefined state) */
bool rational_denominator_zero(const Rational *q);

#endif /* TRTS_RATIONAL_H */
