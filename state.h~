#ifndef TRTS_STATE_H
#define TRTS_STATE_H

#include "rational.h"
#include "config.h"
#include <stdbool.h>
#include <stddef.h>

/* TRTS runtime state. */
typedef struct {
    /* primary registers */
    Rational upsilon;
    Rational beta;
    Rational koppa;
    /* supplementary registers */
    Rational epsilon;
    Rational phi;
    Rational previous_upsilon;
    Rational previous_beta;
    Rational delta_upsilon;
    Rational delta_beta;
    /* triangles */
    Rational triangle_phi_over_epsilon;
    Rational triangle_prev_over_phi;
    Rational triangle_epsilon_over_prev;
    /* koppa stack */
    Rational koppa_stack[4];
    size_t koppa_stack_size;
    Rational koppa_sample;
    int koppa_sample_index;
    /* flags */
    bool rho_pending;
    bool rho_latched;
    bool psi_recent;
    bool psi_triple_recent;
    bool psi_strength_applied;
    bool ratio_triggered_recent;
    bool ratio_threshold_recent;
    bool dual_engine_last_step;
    bool sign_flip_polarity;
    /* current tick */
    size_t tick;
} TRTS_State;

void state_init(TRTS_State *st);
void state_clear(TRTS_State *st);
void state_reset(TRTS_State *st, const Config *cfg);

#endif /* TRTS_STATE_H */

