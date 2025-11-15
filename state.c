// state.c
/*
==============================================
     TRTS SYSTEM CREED – RATIONAL ONLY
==============================================
*/

#include "state.h"

#include "rational.h"

static void state_zero(mpq_t value) {
    rational_set_si(value, 0, 1);
}

void state_init(TRTS_State *state) {
    rational_init(state->upsilon);
    rational_init(state->beta);
    rational_init(state->koppa);
    rational_init(state->epsilon);
    rational_init(state->phi);
    rational_init(state->previous_upsilon);
    rational_init(state->previous_beta);
    rational_init(state->delta_upsilon);
    rational_init(state->delta_beta);
    rational_init(state->triangle_phi_over_epsilon);
    rational_init(state->triangle_prev_over_phi);
    rational_init(state->triangle_epsilon_over_prev);
    for (size_t i = 0; i < 4; ++i) {
        rational_init(state->koppa_stack[i]);
        state_zero(state->koppa_stack[i]);
    }
    rational_init(state->koppa_sample);

    state_zero(state->upsilon);
    state_zero(state->beta);
    state_zero(state->koppa);
    state_zero(state->epsilon);
    state_zero(state->phi);
    state_zero(state->previous_upsilon);
    state_zero(state->previous_beta);
    state_zero(state->delta_upsilon);
    state_zero(state->delta_beta);
    state_zero(state->triangle_phi_over_epsilon);
    state_zero(state->triangle_prev_over_phi);
    state_zero(state->triangle_epsilon_over_prev);
    state_zero(state->koppa_sample);

    state->koppa_stack_size = 0;
    state->koppa_sample_index = -1;
    state->rho_pending = false;
    state->rho_latched = false;
    state->psi_recent = false;
    state->ratio_triggered_recent = false;
    state->psi_triple_recent = false;
    state->dual_engine_last_step = false;
    state->ratio_threshold_recent = false;
    state->psi_strength_applied = false;
    state->sign_flip_polarity = false;
    state->tick = 0;
}

void state_clear(TRTS_State *state) {
    rational_clear(state->upsilon);
    rational_clear(state->beta);
    rational_clear(state->koppa);
    rational_clear(state->epsilon);
    rational_clear(state->phi);
    rational_clear(state->previous_upsilon);
    rational_clear(state->previous_beta);
    rational_clear(state->delta_upsilon);
    rational_clear(state->delta_beta);
    rational_clear(state->triangle_phi_over_epsilon);
    rational_clear(state->triangle_prev_over_phi);
    rational_clear(state->triangle_epsilon_over_prev);
    for (size_t i = 0; i < 4; ++i) {
        rational_clear(state->koppa_stack[i]);
    }
    rational_clear(state->koppa_sample);
}

void state_reset(TRTS_State *state, const Config *config) {
    rational_set(state->upsilon, config->initial_upsilon);
    rational_set(state->beta, config->initial_beta);
    rational_set(state->koppa, config->initial_koppa);
    rational_set(state->epsilon, config->initial_upsilon);
    rational_set(state->phi, config->initial_beta);
    rational_set(state->previous_upsilon, config->initial_upsilon);
    rational_set(state->previous_beta, config->initial_beta);
    state_zero(state->delta_upsilon);
    state_zero(state->delta_beta);
    state_zero(state->triangle_phi_over_epsilon);
    state_zero(state->triangle_prev_over_phi);
    state_zero(state->triangle_epsilon_over_prev);
    for (size_t i = 0; i < 4; ++i) {
        state_zero(state->koppa_stack[i]);
    }
    rational_set(state->koppa_sample, state->koppa);
    state->koppa_stack_size = 0;
    state->koppa_sample_index = -1;
    state->rho_pending = false;
    state->rho_latched = false;
    state->psi_recent = false;
    state->ratio_triggered_recent = false;
    state->psi_triple_recent = false;
    state->dual_engine_last_step = false;
    state->ratio_threshold_recent = false;
    state->psi_strength_applied = false;
    state->sign_flip_polarity = false;
    state->tick = 0;
}

