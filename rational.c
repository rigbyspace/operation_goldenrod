#include "rational.h"
#include <stdlib.h>

void rational_init(Rational *q) {
    mpz_init(q->num);
    mpz_init(q->den);
    mpz_set_si(q->num, 0);
    mpz_set_ui(q->den, 1);
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
        /* 0 numerator forces denominator to zero (0/0). */
        mpz_set_ui(q->den, 0);
    }
}

void rational_set_components(Rational *q, const mpz_t num, const mpz_t den) {
    mpz_set(q->num, num);
    mpz_set(q->den, den);
    if (mpz_sgn(num) == 0) {
        mpz_set_ui(q->den, 0);
    }
}

static void normalise_zero(Rational *q) {
    if (mpz_sgn(q->num) == 0) {
        mpz_set_ui(q->den, 0);
    }
}

void rational_add(Rational *r, const Rational *a, const Rational *b) {
    /* r = a + b = a.num*b.den + b.num*a.den / (a.den*b.den) */
    mpz_t n, d;
    mpz_init(n);
    mpz_init(d);
    mpz_mul(n, a->num, b->den);
    mpz_t tmp;
    mpz_init(tmp);
    mpz_mul(tmp, b->num, a->den);
    mpz_add(n, n, tmp);
    mpz_mul(d, a->den, b->den);
    mpz_clear(tmp);
    rational_set_components(r, n, d);
    mpz_clear(n);
    mpz_clear(d);
    normalise_zero(r);
}

void rational_sub(Rational *r, const Rational *a, const Rational *b) {
    mpz_t n, d;
    mpz_init(n);
    mpz_init(d);
    mpz_mul(n, a->num, b->den);
    mpz_t tmp;
    mpz_init(tmp);
    mpz_mul(tmp, b->num, a->den);
    mpz_sub(n, n, tmp);
    mpz_mul(d, a->den, b->den);
    mpz_clear(tmp);
    rational_set_components(r, n, d);
    mpz_clear(n);
    mpz_clear(d);
    normalise_zero(r);
}

void rational_mul(Rational *r, const Rational *a, const Rational *b) {
    mpz_t n, d;
    mpz_init(n);
    mpz_init(d);
    mpz_mul(n, a->num, b->num);
    mpz_mul(d, a->den, b->den);
    rational_set_components(r, n, d);
    mpz_clear(n);
    mpz_clear(d);
    normalise_zero(r);
}

bool rational_div(Rational *r, const Rational *a, const Rational *b) {
    if (mpz_sgn(b->num) == 0) {
        return false;
    }
    mpz_t n, d;
    mpz_init(n);
    mpz_init(d);
    mpz_mul(n, a->num, b->den);
    mpz_mul(d, a->den, b->num);
    rational_set_components(r, n, d);
    mpz_clear(n);
    mpz_clear(d);
    normalise_zero(r);
    return true;
}

void rational_negate(Rational *q) {
    mpz_neg(q->num, q->num);
    normalise_zero(q);
}

void rational_abs(Rational *q) {
    mpz_abs(q->num, q->num);
    mpz_abs(q->den, q->den);
    normalise_zero(q);
}

void rational_copy_num(mpz_t dest, const Rational *q) {
    mpz_set(dest, q->num);
}

bool rational_is_zero(const Rational *q) {
    return mpz_sgn(q->num) == 0;
}

int rational_cmp(const Rational *a, const Rational *b) {
    /* Compare without reducing: compare a.num*b.den and b.num*a.den */
    mpz_t lhs, rhs;
    mpz_init(lhs);
    mpz_init(rhs);
    mpz_mul(lhs, a->num, b->den);
    mpz_mul(rhs, b->num, a->den);
    int cmp = mpz_cmp(lhs, rhs);
    mpz_clear(lhs);
    mpz_clear(rhs);
    return cmp;
}

