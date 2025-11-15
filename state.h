// state.h
/*
==============================================
     TRTS SYSTEM CREED – RATIONAL ONLY
==============================================
*/

#ifndef STATE_H
#define STATE_H

#include "stdio.h"
#include <gmp.h>
#include <stdbool.h>
#include <stddef.h>

#include "config.h"

typedef struct {
    mpq_t upsilon;
    mpq_t beta;
    mpq_t koppa;
    mpq_t epsilon;
    mpq_t phi;
    mpq_t previous_upsilon;
    mpq_t previous_beta;
    mpq_t delta_upsilon;
    mpq_t delta_beta;
    mpq_t triangle_phi_over_epsilon;
    mpq_t triangle_prev_over_phi;
    mpq_t triangle_epsilon_over_prev;
    mpq_t koppa_stack[4];
    size_t koppa_stack_size;
    mpq_t koppa_sample;
    int koppa_sample_index;
    bool rho_pending;
    bool rho_latched;
    bool psi_recent;
    bool ratio_triggered_recent;
    bool psi_triple_recent;
    bool dual_engine_last_step;
    bool ratio_threshold_recent;
    bool psi_strength_applied;
    bool sign_flip_polarity;
    size_t tick;
} TRTS_State;

void state_init(TRTS_State *state);
void state_clear(TRTS_State *state);
void state_reset(TRTS_State *state, const Config *config);

#endif // STATE_H

