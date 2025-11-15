/* engine.c - TRTS Engine Implementation
 *
 * Core propagation engine with all modulation features.
 * Maintains strict rational arithmetic throughout.
 */

#include "engine.h"
#include "rational.h"
#include <gmp.h>

/* Convert EngineMode to EngineTrackMode */
static EngineTrackMode convert_engine_mode(EngineMode mode) {
    switch (mode) {
        case ENGINE_MODE_ADD:
            return ENGINE_TRACK_ADD;
        case ENGINE_MODE_MULTI:
            return ENGINE_TRACK_MULTI;
        case ENGINE_MODE_SLIDE:
            return ENGINE_TRACK_SLIDE;
        case ENGINE_MODE_DELTA_ADD:
        default:
            return ENGINE_TRACK_ADD;
    }
}

/* Apply asymmetric cascade: different track modes per microtick */
static void apply_asymmetric_modes(const Config *config, int microtick,
                                    EngineTrackMode *ups_mode,
                                    EngineTrackMode *beta_mode) {
    if (!config->enable_asymmetric_cascade) {
        return;
    }
    
    switch (microtick) {
        case 1:
            *ups_mode = ENGINE_TRACK_MULTI;
            *beta_mode = ENGINE_TRACK_ADD;
            break;
        case 4:
            *ups_mode = ENGINE_TRACK_ADD;
            *beta_mode = ENGINE_TRACK_SLIDE;
            break;
        case 7:
            *ups_mode = ENGINE_TRACK_SLIDE;
            *beta_mode = ENGINE_TRACK_MULTI;
            break;
        case 10:
            *ups_mode = ENGINE_TRACK_ADD;
            *beta_mode = ENGINE_TRACK_ADD;
            break;
        default:
            break;
    }
}

/* Modulate track mode based on stack depth */
static EngineTrackMode apply_stack_depth_mode(const Config *config,
                                               const TRTS_State *state,
                                               EngineTrackMode base_mode) {
    if (!config->enable_stack_depth_modes) {
        return base_mode;
    }
    
    size_t depth = state->koppa_stack_size;
    if (depth <= 1) {
        return ENGINE_TRACK_ADD;
    }
    if (depth <= 3) {
        return ENGINE_TRACK_MULTI;
    }
    if (depth == 4) {
        return ENGINE_TRACK_SLIDE;
    }
    return ENGINE_TRACK_ADD;
}

/* Modulate track mode based on koppa magnitude */
static EngineTrackMode apply_koppa_gate(const Config *config,
                                        const TRTS_State *state,
                                        EngineTrackMode base_mode) {
    if (!config->enable_koppa_gated_engine) {
        return base_mode;
    }
    
    mpz_t magnitude;
    mpz_init(magnitude);
    rational_abs_num(magnitude, &state->koppa);
    
    EngineTrackMode result;
    if (mpz_cmp_ui(magnitude, 10UL) < 0) {
        result = ENGINE_TRACK_SLIDE;
    } else if (mpz_cmp_ui(magnitude, 100UL) < 0) {
        result = ENGINE_TRACK_MULTI;
    } else {
        result = ENGINE_TRACK_ADD;
    }
    
    mpz_clear(magnitude);
    return result;
}

/* Apply track mode formula to compute new value */
static bool apply_track_mode(EngineTrackMode mode, Rational *result,
                              const Rational *current, const Rational *counterpart,
                              const Rational *koppa) {
    Rational workspace;
    rational_init(&workspace);
    bool ok = true;
    
    switch (mode) {
        case ENGINE_TRACK_ADD:
            /* result = current + counterpart + koppa */
            rational_add(result, current, counterpart);
            rational_add(result, result, koppa);
            break;
            
        case ENGINE_TRACK_MULTI:
            /* result = current * (counterpart + koppa) */
            rational_add(&workspace, counterpart, koppa);
            rational_mul(result, current, &workspace);
            break;
            
        case ENGINE_TRACK_SLIDE:
            /* result = (current + counterpart) / koppa */
            if (rational_is_zero(koppa) || rational_denominator_zero(koppa)) {
                ok = false;
            } else {
                rational_add(&workspace, current, counterpart);
                if (rational_denominator_zero(&workspace)) {
                    ok = false;
                } else {
                    ok = rational_div(result, &workspace, koppa);
                }
            }
            break;
    }
    
    rational_clear(&workspace);
    return ok;
}

/* Apply sign flip to upsilon and beta */
static void apply_sign_flip(const Config *config, TRTS_State *state,
                             Rational *upsilon, Rational *beta) {
    if (!config->enable_sign_flip || config->sign_flip_mode == SIGN_FLIP_NONE) {
        state->sign_flip_polarity = false;
        return;
    }
    
    bool flip_now = false;
    switch (config->sign_flip_mode) {
        case SIGN_FLIP_ALWAYS:
            flip_now = true;
            break;
        case SIGN_FLIP_ALTERNATE:
            flip_now = !state->sign_flip_polarity;
            break;
        case SIGN_FLIP_NONE:
        default:
            break;
    }
    
    if (flip_now) {
        rational_negate(upsilon);
        rational_negate(beta);
    }
    
    /* Update polarity flag */
    if (config->sign_flip_mode == SIGN_FLIP_ALWAYS) {
        state->sign_flip_polarity = true;
    } else if (config->sign_flip_mode == SIGN_FLIP_ALTERNATE) {
        state->sign_flip_polarity = flip_now;
    }
}

/* Update epsilon-phi triangle ratios */
static void update_triangle(const Config *config, TRTS_State *state) {
    if (!config->enable_epsilon_phi_triangle) {
        return;
    }
    
    /* φ/ε */
    if (!rational_is_zero(&state->epsilon)) {
        rational_div(&state->triangle_phi_over_epsilon, &state->phi, &state->epsilon);
    } else {
        rational_set_si(&state->triangle_phi_over_epsilon, 0, 1);
    }
    
    /* υ_prev/φ */
    if (!rational_is_zero(&state->phi)) {
        rational_div(&state->triangle_prev_over_phi, &state->previous_upsilon, &state->phi);
    } else {
        rational_set_si(&state->triangle_prev_over_phi, 0, 1);
    }
    
    /* ε/υ_prev */
    if (!rational_is_zero(&state->previous_upsilon)) {
        rational_div(&state->triangle_epsilon_over_prev, &state->epsilon, &state->previous_upsilon);
    } else {
        rational_set_si(&state->triangle_epsilon_over_prev, 0, 1);
    }
}

/* Apply delta cross-propagation */
static void apply_delta_cross(const Config *config, TRTS_State *state,
                               Rational *new_upsilon, Rational *new_beta) {
    if (!config->enable_delta_cross_propagation) {
        return;
    }
    
    /* Add Δβ to new_upsilon, Δυ to new_beta */
    rational_add(new_upsilon, new_upsilon, &state->delta_beta);
    rational_add(new_beta, new_beta, &state->delta_upsilon);
    
    /* Optional: add κ offset */
    if (config->enable_delta_koppa_offset) {
        rational_add(new_upsilon, new_upsilon, &state->koppa);
        rational_add(new_beta, new_beta, &state->koppa);
    }
}

/* Reduce rational numerator modulo bound (preserving sign) */
static void rational_mod_bound(Rational *value, const mpz_t bound) {
    if (mpz_cmp_ui(bound, 0UL) == 0) {
        return;
    }
    
    mpz_t rem;
    mpz_init(rem);
    
    /* Store sign */
    int sign = mpz_sgn(mpq_numref(*value));
    
    /* Reduce absolute value mod bound */
    mpz_abs(rem, mpq_numref(*value));
    mpz_mod(rem, rem, bound);
    
    /* Restore sign */
    if (sign < 0) {
        mpz_neg(rem, rem);
    }
    
    mpz_set(mpq_numref(*value), rem);
    mpz_clear(rem);
}

/* Apply modular wrap operations */
static void apply_modular_wrap(const Config *config, TRTS_State *state) {
    if (!config->enable_modular_wrap) {
        return;
    }
    
    /* Wrap κ modulo β when |κ| exceeds threshold */
    if (config->koppa_wrap_threshold > 0) {
        mpz_t magnitude;
        mpz_init(magnitude);
        rational_abs_num(magnitude, &state->koppa);
        
        if (mpz_cmp_ui(magnitude, config->koppa_wrap_threshold) > 0) {
            rational_mod(&state->koppa, &state->koppa, &state->beta);
        }
        
        mpz_clear(magnitude);
    }
    
    /* Reduce numerators by modulus_bound if set */
    if (mpz_cmp_ui(config->modulus_bound, 0UL) > 0) {
        rational_mod_bound(&state->upsilon, config->modulus_bound);
        rational_mod_bound(&state->beta, config->modulus_bound);
        rational_mod_bound(&state->koppa, config->modulus_bound);
    }
}

bool engine_step(const Config *config, TRTS_State *state, int microtick) {
    /* Store previous values */
    Rational ups_before, beta_before;
    rational_init(&ups_before);
    rational_init(&beta_before);
    rational_set(&ups_before, &state->upsilon);
    rational_set(&beta_before, &state->beta);
    
    bool success = true;
    
    /* Determine track modes */
    EngineTrackMode ups_mode, beta_mode;
    if (config->dual_track_mode) {
        ups_mode = config->engine_upsilon;
        beta_mode = config->engine_beta;
    } else {
        ups_mode = convert_engine_mode(config->engine_mode);
        beta_mode = ups_mode;
    }
    
    /* Apply modulations */
    apply_asymmetric_modes(config, microtick, &ups_mode, &beta_mode);
    ups_mode = apply_stack_depth_mode(config, state, ups_mode);
    beta_mode = apply_stack_depth_mode(config, state, beta_mode);
    ups_mode = apply_koppa_gate(config, state, ups_mode);
    beta_mode = apply_koppa_gate(config, state, beta_mode);
    
    /* Prepare new values */
    Rational new_upsilon, new_beta;
    rational_init(&new_upsilon);
    rational_init(&new_beta);
    rational_set(&new_upsilon, &state->upsilon);
    rational_set(&new_beta, &state->beta);
    
    /* Compute deltas from previous step */
    rational_delta(&state->delta_upsilon, &state->upsilon, &state->previous_upsilon);
    rational_delta(&state->delta_beta, &state->beta, &state->previous_beta);
    
    /* Check if using delta-add mode */
    bool use_delta_add = (!config->dual_track_mode && 
                          config->engine_mode == ENGINE_MODE_DELTA_ADD);
    
    if (use_delta_add) {
        /* Delta-add: u' = u + Δu, b' = b + Δb */
        rational_add(&new_upsilon, &state->upsilon, &state->delta_upsilon);
        rational_add(&new_beta, &state->beta, &state->delta_beta);
    } else {
        /* Standard track mode formulas */
        bool ups_success = apply_track_mode(ups_mode, &new_upsilon, 
                                            &state->upsilon, &state->beta, &state->koppa);
        bool beta_success = apply_track_mode(beta_mode, &new_beta, 
                                             &state->beta, &state->upsilon, &state->koppa);
        success = ups_success && beta_success;
    }
    
    /* Apply additional transformations */
    apply_delta_cross(config, state, &new_upsilon, &new_beta);
    apply_sign_flip(config, state, &new_upsilon, &new_beta);
    update_triangle(config, state);
    
    /* Commit new values if successful */
    if (success) {
        rational_set(&state->upsilon, &new_upsilon);
        rational_set(&state->beta, &new_beta);
        state->dual_engine_last_step = config->dual_track_mode;
        
        /* Recompute deltas based on actual change */
        rational_delta(&state->delta_upsilon, &state->upsilon, &ups_before);
        rational_delta(&state->delta_beta, &state->beta, &beta_before);
        
        /* Update previous values */
        rational_set(&state->previous_upsilon, &ups_before);
        rational_set(&state->previous_beta, &beta_before);
        
        /* Apply modular wrap */
        apply_modular_wrap(config, state);
    } else {
        state->dual_engine_last_step = false;
    }
    
    /* Cleanup */
    rational_clear(&ups_before);
    rational_clear(&beta_before);
    rational_clear(&new_upsilon);
    rational_clear(&new_beta);
    
    return success;
}
