/* rational.c - TRTS Rational Implementation
 *
 * Core arithmetic engine enforcing TRTS axioms.
 * All operations avoid canonicalization and preserve raw fraction components.
 */

#include "rational.h"
#include <stdlib.h>

/* Helper: enforce 0/0 invariant after any operation */
static void normalize_zero(Rational *q) {
    if (mpz_sgn(q->num) == 0) {
        mpz_set_ui(q->den, 0UL);
    }
}

void rational_init(Rational *q) {
    mpz_init(q->num);
    mpz_init(q->den);
    mpz_set_si(q->num, 0);
    mpz_set_ui(q->den, 1UL);
}

void rational_clear(Rational *q) {
    mpz_clear(q->num);
    mpz_clear(q->den);
}

void rational_set(Rational *q, const Rational *src) {
    mpz_set(q->num, src->num);
    mpz_set(q->den, src->den);
}

void rational_set_si(Rational *q, long n, unsigned long d) {
    mpz_set_si(q->num, n);
    mpz_set_ui(q->den, d);
    if (n == 0) {
        mpz_set_ui(q->den, 0UL);  /* Enforce 0/0 */
    }
}

void rational_set_components(Rational *q, const mpz_t num, const mpz_t den) {
    mpz_set(q->num, num);
    mpz_set(q->den, den);
    normalize_zero(q);
}

void rational_add(Rational *r, const Rational *a, const Rational *b) {
    /* r = a + b = (a.num * b.den + b.num * a.den) / (a.den * b.den) */
    mpz_t n, d, tmp;
    mpz_init(n);
    mpz_init(d);
    mpz_init(tmp);
    
    mpz_mul(n, a->num, b->den);       /* a.num * b.den */
    mpz_mul(tmp, b->num, a->den);     /* b.num * a.den */
    mpz_add(n, n, tmp);               /* sum */
    mpz_mul(d, a->den, b->den);       /* product of denominators */
    
    rational_set_components(r, n, d);
    
    mpz_clear(n);
    mpz_clear(d);
    mpz_clear(tmp);
}

void rational_sub(Rational *r, const Rational *a, const Rational *b) {
    /* r = a - b = (a.num * b.den - b.num * a.den) / (a.den * b.den) */
    mpz_t n, d, tmp;
    mpz_init(n);
    mpz_init(d);
    mpz_init(tmp);
    
    mpz_mul(n, a->num, b->den);
    mpz_mul(tmp, b->num, a->den);
    mpz_sub(n, n, tmp);               /* difference */
    mpz_mul(d, a->den, b->den);
    
    rational_set_components(r, n, d);
    
    mpz_clear(n);
    mpz_clear(d);
    mpz_clear(tmp);
}

void rational_mul(Rational *r, const Rational *a, const Rational *b) {
    /* r = a * b = (a.num * b.num) / (a.den * b.den) */
    mpz_t n, d;
    mpz_init(n);
    mpz_init(d);
    
    mpz_mul(n, a->num, b->num);
    mpz_mul(d, a->den, b->den);
    
    rational_set_components(r, n, d);
    
    mpz_clear(n);
    mpz_clear(d);
}

bool rational_div(Rational *r, const Rational *a, const Rational *b) {
    /* r = a / b = (a.num * b.den) / (a.den * b.num) */
    if (mpz_sgn(b->num) == 0) {
        return false;  /* Division by zero */
    }
    
    mpz_t n, d;
    mpz_init(n);
    mpz_init(d);
    
    mpz_mul(n, a->num, b->den);
    mpz_mul(d, a->den, b->num);
    
    rational_set_components(r, n, d);
    
    mpz_clear(n);
    mpz_clear(d);
    return true;
}

void rational_negate(Rational *q) {
    mpz_neg(q->num, q->num);
    normalize_zero(q);
}

void rational_abs(Rational *q) {
    mpz_abs(q->num, q->num);
    mpz_abs(q->den, q->den);
    normalize_zero(q);
}

void rational_mod(Rational *r, const Rational *a, const Rational *b) {
    /* Modular reduction: compute a mod b in rational arithmetic
     * r = a - b * floor(a/b) */
    if (rational_is_zero(b)) {
        rational_set(r, a);
        return;
    }
    
    Rational quotient, product;
    rational_init(&quotient);
    rational_init(&product);
    
    rational_div(&quotient, a, b);
    rational_floor(&quotient, &quotient);
    rational_mul(&product, b, &quotient);
    rational_sub(r, a, &product);
    
    rational_clear(&quotient);
    rational_clear(&product);
}

void rational_delta(Rational *r, const Rational *a, const Rational *b) {
    /* Delta: r = a - b */
    rational_sub(r, a, b);
}

void rational_floor(Rational *r, const Rational *q) {
    /* Floor: largest integer <= q */
    if (rational_denominator_zero(q)) {
        rational_set(r, q);
        return;
    }
    
    mpz_t floor_val;
    mpz_init(floor_val);
    mpz_fdiv_q(floor_val, q->num, q->den);
    rational_set_si(r, mpz_get_si(floor_val), 1UL);
    mpz_clear(floor_val);
}

void rational_ceil(Rational *r, const Rational *q) {
    /* Ceiling: smallest integer >= q */
    if (rational_denominator_zero(q)) {
        rational_set(r, q);
        return;
    }
    
    mpz_t ceil_val;
    mpz_init(ceil_val);
    mpz_cdiv_q(ceil_val, q->num, q->den);
    rational_set_si(r, mpz_get_si(ceil_val), 1UL);
    mpz_clear(ceil_val);
}

void rational_round(Rational *r, const Rational *q) {
    /* Round to nearest integer */
    if (rational_denominator_zero(q)) {
        rational_set(r, q);
        return;
    }
    
    mpz_t rounded;
    mpz_init(rounded);
    
    /* Compute 2*num to check fractional part */
    mpz_t doubled_num, rem;
    mpz_init(doubled_num);
    mpz_init(rem);
    mpz_mul_2exp(doubled_num, q->num, 1);
    mpz_fdiv_qr(rounded, rem, doubled_num, q->den);
    
    /* If remainder >= den, round up */
    if (mpz_cmp(rem, q->den) >= 0) {
        mpz_add_ui(rounded, rounded, 1UL);
    }
    
    /* Divide by 2 to get final result */
    mpz_fdiv_q_2exp(rounded, rounded, 1);
    
    rational_set_si(r, mpz_get_si(rounded), 1UL);
    
    mpz_clear(rounded);
    mpz_clear(doubled_num);
    mpz_clear(rem);
}

void rational_abs_num(mpz_t dest, const Rational *q) {
    mpz_abs(dest, q->num);
}

int rational_cmp(const Rational *a, const Rational *b) {
    /* Compare a and b: compute a.num*b.den vs b.num*a.den */
    mpz_t lhs, rhs;
    mpz_init(lhs);
    mpz_init(rhs);
    
    mpz_mul(lhs, a->num, b->den);
    mpz_mul(rhs, b->num, a->den);
    
    int result = mpz_cmp(lhs, rhs);
    
    mpz_clear(lhs);
    mpz_clear(rhs);
    return result;
}

int rational_sgn(const Rational *q) {
    return mpz_sgn(q->num);
}

bool rational_is_zero(const Rational *q) {
    return mpz_sgn(q->num) == 0;
}

bool rational_denominator_zero(const Rational *q) {
    return mpz_sgn(q->den) == 0;
}
