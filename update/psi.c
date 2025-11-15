#include "psi.h"
#include "rational.h"
#include <stdbool.h>

/* Standard psi: (u,b) → (b/u, u/b).  Returns false if division by zero. */
static bool standard_psi(TRTS_State *st) {
    if (rational_is_zero(&st->upsilon) || rational_is_zero(&st->beta)) {
        return false;
    }
    Rational new_u, new_b;
    rational_init(&new_u);
    rational_init(&new_b);
    /* new_u = beta_num*ups_den / beta_den*ups_num */
    mpz_t t1, t2;
    mpz_init(t1);
    mpz_init(t2);
    mpz_mul(t1, st->beta.num, st->upsilon.den);
    mpz_mul(t2, st->beta.den, st->upsilon.num);
    rational_set_components(&new_u, t1, t2);
    /* new_b = ups_num*beta_den / ups_den*beta_num */
    mpz_mul(t1, st->upsilon.num, st->beta.den);
    mpz_mul(t2, st->upsilon.den, st->beta.num);
    rational_set_components(&new_b, t1, t2);
    bool ok = true;
    if (mpz_sgn(new_u.den) == 0 || mpz_sgn(new_b.den) == 0) ok = false;
    if (ok) {
        rational_set(&st->upsilon, &new_u);
        rational_set(&st->beta, &new_b);
    }
    rational_clear(&new_u);
    rational_clear(&new_b);
    mpz_clear(t1);
    mpz_clear(t2);
    return ok;
}

/* Triple psi: (u,b,k) → (b/k, k/u, k/b). */
static bool triple_psi(TRTS_State *st) {
    if (rational_is_zero(&st->upsilon) || rational_is_zero(&st->beta) || rational_is_zero(&st->koppa)) {
        return false;
    }
    Rational new_u, new_b, new_k;
    rational_init(&new_u);
    rational_init(&new_b);
    rational_init(&new_k);
    mpz_t t1, t2;
    mpz_init(t1);
    mpz_init(t2);
    /* new_u = beta/koppa */
    mpz_mul(t1, st->beta.num, st->koppa.den);
    mpz_mul(t2, st->beta.den, st->koppa.num);
    rational_set_components(&new_u, t1, t2);
    /* new_b = koppa/upsilon */
    mpz_mul(t1, st->koppa.num, st->upsilon.den);
    mpz_mul(t2, st->koppa.den, st->upsilon.num);
    rational_set_components(&new_b, t1, t2);
    /* new_k = koppa/beta */
    mpz_mul(t1, st->koppa.num, st->beta.den);
    mpz_mul(t2, st->koppa.den, st->beta.num);
    rational_set_components(&new_k, t1, t2);
    bool ok = true;
    if (mpz_sgn(new_u.den) == 0 || mpz_sgn(new_b.den) == 0 || mpz_sgn(new_k.den) == 0) ok = false;
    if (ok) {
        rational_set(&st->upsilon, &new_u);
        rational_set(&st->beta, &new_b);
        rational_set(&st->koppa, &new_k);
    }
    rational_clear(&new_u);
    rational_clear(&new_b);
    rational_clear(&new_k);
    mpz_clear(t1);
    mpz_clear(t2);
    return ok;
}

/* Count how many numerators are prime. */
static int prime_count(const TRTS_State *st) {
    int count = 0;
    mpz_t m;
    mpz_init(m);
    mpz_abs(m, st->upsilon.num);
    if (mpz_cmp_ui(m, 2UL) >= 0 && mpz_probab_prime_p(m, 25) > 0) count++;
    mpz_abs(m, st->beta.num);
    if (mpz_cmp_ui(m, 2UL) >= 0 && mpz_probab_prime_p(m, 25) > 0) count++;
    mpz_abs(m, st->koppa.num);
    if (mpz_cmp_ui(m, 2UL) >= 0 && mpz_probab_prime_p(m, 25) > 0) count++;
    mpz_clear(m);
    return count;
}

bool psi_transform(const Config *cfg, TRTS_State *st) {
    st->psi_recent = false;
    st->psi_triple_recent = false;
    st->psi_strength_applied = false;
    bool can_fire = st->rho_pending || (cfg->psi_mode == PSI_MODE_MSTEP);
    if (!can_fire) return false;
    /* Determine number of fires */
    int strength = 1;
    if (cfg->enable_psi_strength_parameter && st->rho_pending) {
        int pc = prime_count(st);
        strength = (pc > 0) ? pc : 1;
    }
    bool fired = false;
    for (int i = 0; i < strength; ++i) {
        bool request_triple = cfg->triple_psi_mode;
        if (cfg->enable_conditional_triple_psi) {
            if (prime_count(st) >= 3) request_triple = true;
        }
        if (strength >= 3 && i == strength - 3) request_triple = true;
        bool ok;
        if (request_triple) {
            ok = triple_psi(st);
            if (ok) st->psi_triple_recent = true;
        } else {
            ok = standard_psi(st);
        }
        if (ok) {
            fired = true;
            st->psi_recent = true;
            if (i == 0) st->rho_pending = false;
        } else {
            break;
        }
    }
    return fired;
}

