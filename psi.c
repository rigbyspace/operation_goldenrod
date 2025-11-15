/* psi.c - TRTS Psi Transform Implementation
 *
 * Implements 2-way and 3-way psi transforms with explicit denominator checks.
 * All operations maintain strict rational integrity.
 */

#include "psi.h"
#include "rational.h"
#include <stdbool.h>

/* Standard 2-way psi: (υ,β) → (β/υ, υ/β)
 * Returns false if either division would result in zero denominator */
static bool standard_psi(TRTS_State *st) {
    if (rational_is_zero(&st->upsilon) || rational_is_zero(&st->beta)) {
        return false;
    }
    
    Rational new_u, new_b;
    rational_init(&new_u);
    rational_init(&new_b);
    
    /* new_u = β/υ = (β.num * υ.den) / (β.den * υ.num) */
    mpz_t t1, t2;
    mpz_init(t1);
    mpz_init(t2);
    
    mpz_mul(t1, mpq_numref(st->beta), mpq_denref(st->upsilon));
    mpz_mul(t2, mpq_denref(st->beta), mpq_numref(st->upsilon));
    rational_set_components(&new_u, t1, t2);
    
    /* new_b = υ/β = (υ.num * β.den) / (υ.den * β.num) */
    mpz_mul(t1, mpq_numref(st->upsilon), mpq_denref(st->beta));
    mpz_mul(t2, mpq_denref(st->upsilon), mpq_numref(st->beta));
    rational_set_components(&new_b, t1, t2);
    
    /* Check for zero denominators */
    bool ok = true;
    if (mpz_sgn(mpq_denref(new_u)) == 0 || mpz_sgn(mpq_denref(new_b)) == 0) {
        ok = false;
    }
    
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

/* Triple 3-way psi: (υ,β,κ) → (β/κ, κ/υ, κ/β)
 * Returns false if any division would result in zero denominator */
static bool triple_psi(TRTS_State *st) {
    if (rational_is_zero(&st->upsilon) || 
        rational_is_zero(&st->beta) || 
        rational_is_zero(&st->koppa)) {
        return false;
    }
    
    Rational new_u, new_b, new_k;
    rational_init(&new_u);
    rational_init(&new_b);
    rational_init(&new_k);
    
    mpz_t t1, t2;
    mpz_init(t1);
    mpz_init(t2);
    
    /* new_u = β/κ */
    mpz_mul(t1, mpq_numref(st->beta), mpq_denref(st->koppa));
    mpz_mul(t2, mpq_denref(st->beta), mpq_numref(st->koppa));
    rational_set_components(&new_u, t1, t2);
    
    /* new_b = κ/υ */
    mpz_mul(t1, mpq_numref(st->koppa), mpq_denref(st->upsilon));
    mpz_mul(t2, mpq_denref(st->koppa), mpq_numref(st->upsilon));
    rational_set_components(&new_b, t1, t2);
    
    /* new_k = κ/β */
    mpz_mul(t1, mpq_numref(st->koppa), mpq_denref(st->beta));
    mpz_mul(t2, mpq_denref(st->koppa), mpq_numref(st->beta));
    rational_set_components(&new_k, t1, t2);
    
    /* Check for zero denominators */
    bool ok = true;
    if (mpz_sgn(mpq_denref(new_u)) == 0 || 
        mpz_sgn(mpq_denref(new_b)) == 0 || 
        mpz_sgn(mpq_denref(new_k)) == 0) {
        ok = false;
    }
    
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

/* Count how many of υ, β, κ numerators are prime (for strength parameter) */
static int prime_count(const TRTS_State *st) {
    int count = 0;
    mpz_t m;
    mpz_init(m);
    
    /* Check upsilon numerator */
    mpz_abs(m, mpq_numref(st->upsilon));
    if (mpz_cmp_ui(m, 2UL) >= 0 && mpz_probab_prime_p(m, 25) > 0) {
        count++;
    }
    
    /* Check beta numerator */
    mpz_abs(m, mpq_numref(st->beta));
    if (mpz_cmp_ui(m, 2UL) >= 0 && mpz_probab_prime_p(m, 25) > 0) {
        count++;
    }
    
    /* Check koppa numerator */
    mpz_abs(m, mpq_numref(st->koppa));
    if (mpz_cmp_ui(m, 2UL) >= 0 && mpz_probab_prime_p(m, 25) > 0) {
        count++;
    }
    
    mpz_clear(m);
    return count;
}

bool psi_transform(const Config *cfg, TRTS_State *st) {
    /* Clear flags */
    st->psi_recent = false;
    st->psi_triple_recent = false;
    st->psi_strength_applied = false;
    
    /* Check if firing is allowed based on mode */
    bool can_fire = st->rho_pending || (cfg->psi_mode == PSI_MODE_MSTEP);
    if (!can_fire) {
        return false;
    }
    
    /* Determine number of fires (strength parameter) */
    int strength = 1;
    if (cfg->enable_psi_strength_parameter && st->rho_pending) {
        int pc = prime_count(st);
        strength = (pc > 0) ? pc : 1;
        if (strength > 1) {
            st->psi_strength_applied = true;
        }
    }
    
    /* Execute transforms */
    bool fired = false;
    for (int i = 0; i < strength; ++i) {
        /* Determine if this fire should be triple */
        bool request_triple = cfg->triple_psi_mode;
        
        /* Conditional triple: force triple if 3+ primes present */
        if (cfg->enable_conditional_triple_psi) {
            if (prime_count(st) >= 3) {
                request_triple = true;
            }
        }
        
        /* Last 3 fires of strength burst use triple */
        if (strength >= 3 && i >= strength - 3) {
            request_triple = true;
        }
        
        /* Execute transform */
        bool ok;
        if (request_triple) {
            ok = triple_psi(st);
            if (ok) {
                st->psi_triple_recent = true;
            }
        } else {
            ok = standard_psi(st);
        }
        
        if (ok) {
            fired = true;
            st->psi_recent = true;
            /* Clear rho_pending on first successful fire */
            if (i == 0) {
                st->rho_pending = false;
            }
        } else {
            /* Stop on first failure */
            break;
        }
    }
    
    return fired;
}
