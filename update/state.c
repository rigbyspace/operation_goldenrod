#include "state.h"

void state_init(TRTS_State *st) {
    rational_init(&st->upsilon);
    rational_init(&st->beta);
    rational_init(&st->koppa);
    rational_init(&st->epsilon);
    rational_init(&st->phi);
    rational_init(&st->previous_upsilon);
    rational_init(&st->previous_beta);
    rational_init(&st->delta_upsilon);
    rational_init(&st->delta_beta);
    rational_init(&st->triangle_phi_over_epsilon);
    rational_init(&st->triangle_prev_over_phi);
    rational_init(&st->triangle_epsilon_over_prev);
    for (int i = 0; i < 4; ++i) {
        rational_init(&st->koppa_stack[i]);
    }
    rational_init(&st->koppa_sample);
    st->koppa_stack_size = 0;
    st->koppa_sample_index = -1;
    st->rho_pending = false;
    st->rho_latched = false;
    st->psi_recent = false;
    st->psi_triple_recent = false;
    st->psi_strength_applied = false;
    st->ratio_triggered_recent = false;
    st->ratio_threshold_recent = false;
    st->dual_engine_last_step = false;
    st->sign_flip_polarity = false;
    st->tick = 0;
}

void state_clear(TRTS_State *st) {
    rational_clear(&st->upsilon);
    rational_clear(&st->beta);
    rational_clear(&st->koppa);
    rational_clear(&st->epsilon);
    rational_clear(&st->phi);
    rational_clear(&st->previous_upsilon);
    rational_clear(&st->previous_beta);
    rational_clear(&st->delta_upsilon);
    rational_clear(&st->delta_beta);
    rational_clear(&st->triangle_phi_over_epsilon);
    rational_clear(&st->triangle_prev_over_phi);
    rational_clear(&st->triangle_epsilon_over_prev);
    for (int i = 0; i < 4; ++i) {
        rational_clear(&st->koppa_stack[i]);
    }
    rational_clear(&st->koppa_sample);
}

void state_reset(TRTS_State *st, const Config *cfg) {
    rational_set(&st->upsilon, &cfg->initial_upsilon);
    rational_set(&st->beta, &cfg->initial_beta);
    rational_set(&st->koppa, &cfg->initial_koppa);
    /* zero supplementary registers */
    rational_set_si(&st->epsilon, 0, 1);
    rational_set_si(&st->phi, 0, 1);
    rational_set(&st->previous_upsilon, &st->upsilon);
    rational_set(&st->previous_beta, &st->beta);
    rational_set_si(&st->delta_upsilon, 0, 1);
    rational_set_si(&st->delta_beta, 0, 1);
    rational_set_si(&st->triangle_phi_over_epsilon, 0, 1);
    rational_set_si(&st->triangle_prev_over_phi, 0, 1);
    rational_set_si(&st->triangle_epsilon_over_prev, 0, 1);
    /* clear koppa stack */
    st->koppa_stack_size = 0;
    st->koppa_sample_index = -1;
    rational_set_si(&st->koppa_sample, 0, 1);
    /* reset flags */
    st->rho_pending = false;
    st->rho_latched = false;
    st->psi_recent = false;
    st->psi_triple_recent = false;
    st->psi_strength_applied = false;
    st->ratio_triggered_recent = false;
    st->ratio_threshold_recent = false;
    st->dual_engine_last_step = false;
    st->sign_flip_polarity = false;
    st->tick = 0;
}

