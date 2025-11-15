// engine.c
/*
 * engine.c - Engine step logic with safer modular bound and denominator checks
 */

#include "engine.h"

#include <gmp.h>
#include <stdbool.h>

#include "rational.h"

// Forward declarations of helpers
static EngineTrackMode convert_engine_mode(EngineMode mode);
static void apply_asymmetric_modes(const Config *config, int microtick,
                                    EngineTrackMode *ups_mode,
                                    EngineTrackMode *beta_mode);
static EngineTrackMode apply_stack_depth_mode(const Config *config,
                                               const TRTS_State *state,
                                               EngineTrackMode base_mode);
static EngineTrackMode apply_koppa_gate(const Config *config,
                                         const TRTS_State *state,
                                         EngineTrackMode base_mode);
static bool apply_track_mode(EngineTrackMode mode, mpq_t result,
                              mpq_srcptr current, mpq_srcptr counterpart,
                              mpq_srcptr koppa);
static void apply_sign_flip(const Config *config, TRTS_State *state,
                             mpq_t upsilon, mpq_t beta);
static void update_triangle(const Config *config, TRTS_State *state);
static void apply_delta_cross(const Config *config, TRTS_State *state,
                               mpq_t new_upsilon, mpq_t new_beta);
static void apply_modular_wrap(const Config *config, TRTS_State *state);
static void rational_mod_bound(mpq_t value, const mpz_t bound);

/* ===========================================================
   Helper functions
   =========================================================== */

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

static EngineTrackMode apply_koppa_gate(const Config *config,
                                         const TRTS_State *state,
                                         EngineTrackMode base_mode) {
    if (!config->enable_koppa_gated_engine) {
        return base_mode;
    }
    mpz_t magnitude;
    mpz_init(magnitude);
    rational_abs_num(magnitude, state->koppa);
    EngineTrackMode result = base_mode;
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

static bool apply_track_mode(EngineTrackMode mode, mpq_t result,
                              mpq_srcptr current, mpq_srcptr counterpart,
                              mpq_srcptr koppa) {
    mpq_t workspace;
    rational_init(workspace);
    bool ok = true;
    switch (mode) {
    case ENGINE_TRACK_ADD:
        rational_add(result, current, counterpart);
        rational_add(result, result, koppa);
        break;
    case ENGINE_TRACK_MULTI:
        rational_add(workspace, counterpart, koppa);
        rational_mul(result, current, workspace);
        break;
    case ENGINE_TRACK_SLIDE:
        if (rational_is_zero(koppa) || rational_denominator_zero(koppa)) {
            ok = false;
        } else {
            rational_add(workspace, current, counterpart);
            if (rational_denominator_zero(workspace)) {
                ok = false;
            } else {
                rational_div(result, workspace, koppa);
            }
        }
        break;
    }
    rational_clear(workspace);
    return ok;
}

static void apply_sign_flip(const Config *config, TRTS_State *state,
                             mpq_t upsilon, mpq_t beta) {
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
    if (config->sign_flip_mode == SIGN_FLIP_ALWAYS) {
        state->sign_flip_polarity = true;
    } else if (config->sign_flip_mode == SIGN_FLIP_ALTERNATE) {
        state->sign_flip_polarity = flip_now;
    }
}

static void update_triangle(const Config *config, TRTS_State *state) {
    if (!config->enable_epsilon_phi_triangle) {
        return;
    }
    if (!rational_is_zero(state->epsilon)) {
        rational_div(state->triangle_phi_over_epsilon, state->phi, state->epsilon);
    } else {
        rational_set_si(state->triangle_phi_over_epsilon, 0, 1);
    }
    if (!rational_is_zero(state->phi)) {
        rational_div(state->triangle_prev_over_phi, state->previous_upsilon, state->phi);
    } else {
        rational_set_si(state->triangle_prev_over_phi, 0, 1);
    }
    if (!rational_is_zero(state->previous_upsilon)) {
        rational_div(state->triangle_epsilon_over_prev, state->epsilon, state->previous_upsilon);
    } else {
        rational_set_si(state->triangle_epsilon_over_prev, 0, 1);
    }
}

static void apply_delta_cross(const Config *config, TRTS_State *state,
                               mpq_t new_upsilon, mpq_t new_beta) {
    if (!config->enable_delta_cross_propagation) {
        return;
    }
    rational_add(new_upsilon, new_upsilon, state->delta_beta);
    rational_add(new_beta, new_beta, state->delta_upsilon);
    if (config->enable_delta_koppa_offset) {
        rational_add(new_upsilon, new_upsilon, state->koppa);
        rational_add(new_beta, new_beta, state->koppa);
    }
}

// Reduce a rational's numerator modulo bound, preserving sign.
// Denominator is left unchanged.
// If bound is zero this function does nothing.
static void rational_mod_bound(mpq_t value, const mpz_t bound) {
    if (mpz_cmp_ui(bound, 0UL) == 0) {
        return;
    }
    
    mpz_ptr num = mpq_numref(value);
    mpz_t rem;
    mpz_init(rem);
    
    // Store original sign
    int original_sign = mpz_sgn(num);
    
    // Get absolute value and reduce modulo bound
    mpz_abs(rem, num);
    mpz_mod(rem, rem, bound);
    
    // Restore original sign
    if (original_sign < 0) {
        mpz_neg(rem, rem);
    }
    
    mpz_set(num, rem);
    mpz_clear(rem);
}

static void apply_modular_wrap(const Config *config, TRTS_State *state) {
    if (!config->enable_modular_wrap) {
        return;
    }
    // Original behaviour: wrap κ modulo β when |κ| > koppa_wrap_threshold
    mpz_t magnitude;
    mpz_init(magnitude);
    rational_abs_num(magnitude, state->koppa);
    if (mpz_cmp_ui(magnitude, config->koppa_wrap_threshold) > 0) {
        // Wrap κ by β
        rational_mod(state->koppa, state->koppa, state->beta);
    }
    mpz_clear(magnitude);
    
    // New behaviour: reduce upsilon, beta and koppa numerators by modulus_bound if set
    if (mpz_cmp_ui(config->modulus_bound, 0UL) > 0) {
        rational_mod_bound(state->upsilon, config->modulus_bound);
        rational_mod_bound(state->beta, config->modulus_bound);
        rational_mod_bound(state->koppa, config->modulus_bound);
    }
}

/* ===========================================================
   ENGINE STEP
   =========================================================== */

bool engine_step(const Config *config, TRTS_State *state, int microtick) {
    mpq_t ups_before;
    mpq_t beta_before;
    rational_init(ups_before);
    rational_init(beta_before);
    rational_set(ups_before, state->upsilon);
    rational_set(beta_before, state->beta);
    
    bool success = true;
    
    EngineTrackMode ups_mode = config->dual_track_mode ? config->engine_upsilon
                                                      : convert_engine_mode(config->engine_mode);
    EngineTrackMode beta_mode = config->dual_track_mode ? config->engine_beta : ups_mode;
    
    apply_asymmetric_modes(config, microtick, &ups_mode, &beta_mode);
    ups_mode = apply_stack_depth_mode(config, state, ups_mode);
    beta_mode = apply_stack_depth_mode(config, state, beta_mode);
    ups_mode = apply_koppa_gate(config, state, ups_mode);
    beta_mode = apply_koppa_gate(config, state, beta_mode);
    
    mpq_t new_upsilon;
    mpq_t new_beta;
    rational_init(new_upsilon);
    rational_init(new_beta);
    rational_set(new_upsilon, state->upsilon);
    rational_set(new_beta, state->beta);
    
    bool use_delta_add = (!config->dual_track_mode && config->engine_mode == ENGINE_MODE_DELTA_ADD);
    
    rational_delta(state->delta_upsilon, state->upsilon, state->previous_upsilon);
    rational_delta(state->delta_beta, state->beta, state->previous_beta);
    
    if (use_delta_add) {
        rational_add(new_upsilon, state->upsilon, state->delta_upsilon);
        rational_add(new_beta, state->beta, state->delta_beta);
    } else {
        bool ups_success = apply_track_mode(ups_mode, new_upsilon, state->upsilon,
                                            state->beta, state->koppa);
        bool beta_success = apply_track_mode(beta_mode, new_beta, state->beta,
                                             state->upsilon, state->koppa);
        success = success && ups_success && beta_success;
    }
    
    apply_delta_cross(config, state, new_upsilon, new_beta);
    apply_sign_flip(config, state, new_upsilon, new_beta);
    update_triangle(config, state);
    apply_modular_wrap(config, state);
    
    if (success) {
        rational_set(state->upsilon, new_upsilon);
        rational_set(state->beta, new_beta);
        state->dual_engine_last_step = config->dual_track_mode;
        rational_delta(state->delta_upsilon, state->upsilon, ups_before);
        rational_delta(state->delta_beta, state->beta, beta_before);
        rational_set(state->previous_upsilon, ups_before);
        rational_set(state->previous_beta, beta_before);
    } else {
        state->dual_engine_last_step = false;
    }
    
    rational_clear(ups_before);
    rational_clear(beta_before);
    rational_clear(new_upsilon);
    rational_clear(new_beta);
    return success;
}

