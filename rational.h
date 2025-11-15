#ifndef TRTS_RATIONAL_H
#define TRTS_RATIONAL_H

#include <gmp.h>

/*
 * TRTS rational number type.
 * A rational is represented by a pair of integers (num, den) stored in mpz_t.
 * Zero numerators force the denominator to zero (0/0).  Denominator zero
 * represents ∞ or undefined but is preserved through operations.  No
 * canonicalisation or GCD reduction is performed; numerators and denominators
 * may grow large.  All functions avoid mpq_t and avoid any call that
 * normalises fractions.
 */

typedef struct {
    mpz_t num;
    mpz_t den;
} Rational;

/* Initialise a rational to 0/1. */
void rational_init(Rational *q);

/* Clear a rational’s internal mpz values. */
void rational_clear(Rational *q);

/* Set q ← src. */
void rational_set(Rational *q, const Rational *src);

/* Set q ← n/d. */
void rational_set_si(Rational *q, long n, unsigned long d);

/* Set q ← (num, den) using existing mpz_t values (they are copied). */
void rational_set_components(Rational *q, const mpz_t num, const mpz_t den);

/* Addition: r ← a + b. */
void rational_add(Rational *r, const Rational *a, const Rational *b);

/* Subtraction: r ← a – b. */
void rational_sub(Rational *r, const Rational *a, const Rational *b);

/* Multiplication: r ← a × b. */
void rational_mul(Rational *r, const Rational *a, const Rational *b);

/* Division: r ← a ÷ b.  If b.num = 0, r is unchanged and returns false. */
bool rational_div(Rational *r, const Rational *a, const Rational *b);

/* Negation: q ← –q. */
void rational_negate(Rational *q);

/* Absolute value: q ← |q|.  Both numerator and denominator become non‑negative. */
void rational_abs(Rational *q);

/* Copy numerator into mpz dest. */
void rational_copy_num(mpz_t dest, const Rational *q);

/* Return true if numerator is zero (denominator must then be zero). */
bool rational_is_zero(const Rational *q);

/* Compare two rationals without reducing.  Return <0,0,>0 if a<b, a=b, a>b. */
int rational_cmp(const Rational *a, const Rational *b);

#endif /* TRTS_RATIONAL_H */

