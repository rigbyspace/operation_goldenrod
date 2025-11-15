/* state.h - TRTS Runtime State
 *
 * Defines TRTS_State containing all rational registers, flags, counters,
 * and auxiliary structures required for deterministic propagation.
 *
 * Blueprint reference: state.h/state.c from Section 1.
 */

#ifndef TRTS_STATE_H
#define TRTS_STATE_H

#include "rational.h"
#include <stdbool.h>
#include <stddef.h>

/* Forward declaration */
typedef struct Config_s Config;

/* TRTS Runtime State
 *
 * Contains:
 * - Primary registers: upsilon (υ), beta (β), koppa (κ)
 * - Supplementary registers: epsilon (ε), phi (φ)
 * - Previous values for delta computation
 * - Delta values (change from previous step)
 * - Triangle ratios (ε/φ relationships)
 * - Koppa stack (4-level FIFO for multi-level mode)
 * - Sample copy for stack selection
 * - Flags for ρ, ψ, ratio triggers, engine modes
 * - Counters for stack size, sample index, tick number
 */
typedef struct {
    /* Primary registers */
    Rational upsilon;         /* υ - primary state variable */
    Rational beta;            /* β - memory register */
    Rational koppa;           /* κ - accumulator/phantom register */
    
    /* Supplementary registers */
    Rational epsilon;         /* ε - computed at E phase */
    Rational phi;             /* φ - auxiliary register */
    
    /* Previous values */
    Rational previous_upsilon;
    Rational previous_beta;
    
    /* Delta values */
    Rational delta_upsilon;   /* Δυ = υ_current - υ_previous */
    Rational delta_beta;      /* Δβ = β_current - β_previous */
    
    /* Triangle ratios (for epsilon-phi triangle feature) */
    Rational triangle_phi_over_epsilon;     /* φ/ε */
    Rational triangle_prev_over_phi;        /* υ_prev/φ */
    Rational triangle_epsilon_over_prev;    /* ε/υ_prev */
    
    /* Koppa stack (4 entries, FIFO) */
    Rational koppa_stack[4];
    size_t koppa_stack_size;                /* Current depth (0-4) */
    
    /* Sample copy (selected from stack at specific microticks) */
    Rational koppa_sample;
    int koppa_sample_index;                 /* Index of selected sample, or -1 */
    
    /* Flags */
    bool rho_pending;                       /* ρ event pending (trigger for ψ) */
    bool rho_latched;                       /* ρ latched for microtick */
    bool psi_recent;                        /* ψ fired recently */
    bool psi_triple_recent;                 /* 3-way ψ fired */
    bool psi_strength_applied;              /* ψ strength parameter used */
    bool ratio_triggered_recent;            /* Ratio trigger activated */
    bool ratio_threshold_recent;            /* Ratio threshold trigger */
    bool dual_engine_last_step;             /* Dual-track engine used */
    bool sign_flip_polarity;                /* Sign-flip state */
    
    /* Tick counter */
    size_t tick;                            /* Current tick number (1-based) */
} TRTS_State;

/* Initialize state structure (allocate all rationals, set to zero) */
void state_init(TRTS_State *st);

/* Clear state structure (free all GMP resources) */
void state_clear(TRTS_State *st);

/* Reset state to initial configuration from Config
 * Loads seeds, zeroes flags, resets tick counter */
void state_reset(TRTS_State *st, const Config *cfg);

#endif /* TRTS_STATE_H */
