/* config.h - TRTS Configuration
 *
 * Defines all enumerations and configuration parameters controlling
 * TRTS engine behavior. Blueprint reference: config.h/config.c from Section 1.
 */

#ifndef TRTS_CONFIG_H
#define TRTS_CONFIG_H

#include "rational.h"
#include <stdbool.h>
#include <stddef.h>
#include <gmp.h>

/* Engine track mode: determines update formula for upsilon/beta */
typedef enum {
    ENGINE_TRACK_ADD,      /* u' = u + b + κ */
    ENGINE_TRACK_MULTI,    /* u' = u * (b + κ) */
    ENGINE_TRACK_SLIDE     /* u' = (u + b) / κ */
} EngineTrackMode;

/* Engine mode: determines default track mode */
typedef enum {
    ENGINE_MODE_ADD,       /* Additive propagation */
    ENGINE_MODE_MULTI,     /* Multiplicative propagation */
    ENGINE_MODE_SLIDE,     /* Sliding division propagation */
    ENGINE_MODE_DELTA_ADD  /* Delta-based additive */
} EngineMode;

/* Psi (ψ) firing mode */
typedef enum {
    PSI_MODE_MSTEP,        /* Fire on all memory steps */
    PSI_MODE_RHO_ONLY,     /* Fire only when ρ pending */
    PSI_MODE_MSTEP_RHO,    /* Fire on memory steps or ρ */
    PSI_MODE_INHIBIT_RHO   /* Fire when ρ NOT pending */
} PsiMode;

/* Koppa (κ) operation mode */
typedef enum {
    KOPPA_MODE_DUMP,       /* κ ← 0 */
    KOPPA_MODE_POP,        /* κ ← ε */
    KOPPA_MODE_ACCUMULATE  /* κ ← κ + ε */
} KoppaMode;

/* Koppa trigger: when to execute koppa operation */
typedef enum {
    KOPPA_ON_PSI,          /* Trigger only on ψ fire */
    KOPPA_ON_MU_AFTER_PSI, /* Trigger on memory step after ψ */
    KOPPA_ON_ALL_MU        /* Trigger on all memory steps */
} KoppaTrigger;

/* Prime target: which value to check for patterns */
typedef enum {
    PRIME_ON_MEMORY,       /* Check β numerator/denominator */
    PRIME_ON_NEW_UPSILON   /* Check υ numerator after engine step */
} PrimeTarget;

/* Microtick 10 behavior */
typedef enum {
    MT10_FORCED_EMISSION_ONLY,  /* Only emit, no special behavior */
    MT10_FORCED_PSI,            /* Force ρ pending (triggers ψ) */
    MT10_FORCED_ENGINE,         /* Force extra engine step */
    MT10_FORCED_KOPPA           /* Force koppa operation */
} Mt10Behavior;

/* Sign-flip mode */
typedef enum {
    SIGN_FLIP_NONE,        /* No sign flipping */
    SIGN_FLIP_ALWAYS,      /* Flip sign every step */
    SIGN_FLIP_ALTERNATE    /* Alternate flipping based on polarity */
} SignFlipMode;

/* Ratio trigger mode (for ρ activation based on υ/β ratio) */
typedef enum {
    RATIO_TRIGGER_NONE,    /* No ratio trigger */
    RATIO_TRIGGER_GOLDEN,  /* Trigger near golden ratio */
    RATIO_TRIGGER_SQRT2,   /* Trigger near √2 */
    RATIO_TRIGGER_PLASTIC, /* Trigger near plastic constant */
    RATIO_TRIGGER_CUSTOM   /* Use custom range */
} RatioTriggerMode;

/* Configuration structure
 *
 * Contains all parameters controlling TRTS runtime behavior.
 * Separated into:
 * - Mode enumerations
 * - Feature flags (boolean enables)
 * - Rational seeds
 * - Numeric parameters
 */
typedef struct Config_s {
    /* Core modes */
    EngineMode engine_mode;
    EngineTrackMode engine_upsilon;
    EngineTrackMode engine_beta;
    PsiMode psi_mode;
    KoppaMode koppa_mode;
    KoppaTrigger koppa_trigger;
    PrimeTarget prime_target;
    Mt10Behavior mt10_behavior;
    SignFlipMode sign_flip_mode;
    RatioTriggerMode ratio_trigger_mode;
    
    /* Feature flags */
    bool dual_track_mode;                    /* Independent υ/β track modes */
    bool triple_psi_mode;                    /* Use 3-way ψ transform */
    bool multi_level_koppa;                  /* Enable 4-level κ stack */
    bool enable_asymmetric_cascade;          /* Asymmetric track modes per microtick */
    bool enable_conditional_triple_psi;      /* Trigger 3-way ψ on 3+ primes */
    bool enable_koppa_gated_engine;          /* Modulate track mode by |κ| */
    bool enable_delta_cross_propagation;     /* Cross-couple Δυ and Δβ */
    bool enable_delta_koppa_offset;          /* Add κ to delta cross */
    bool enable_ratio_threshold_psi;         /* Trigger ψ on extreme υ/β ratios */
    bool enable_stack_depth_modes;           /* Modulate track by stack depth */
    bool enable_epsilon_phi_triangle;        /* Compute ε/φ triangle ratios */
    bool enable_sign_flip;                   /* Enable sign flipping */
    bool enable_modular_wrap;                /* Wrap κ modulo β when large */
    bool enable_psi_strength_parameter;      /* Multiple ψ fires based on primes */
    bool enable_ratio_custom_range;          /* Use custom ratio trigger bounds */
    bool enable_twin_prime_trigger;          /* Check for twin primes */
    bool enable_fibonacci_trigger;           /* Check for Fibonacci numbers */
    bool enable_perfect_power_trigger;       /* Check for perfect powers */
    bool enable_ratio_snapshot_logging;      /* Log ratio snapshots (analysis) */
    bool enable_feedback_oscillator;         /* Feedback oscillation mode */
    bool enable_fibonacci_gate;              /* Fibonacci gating */
    
    /* Simulation parameters */
    size_t ticks;                            /* Number of ticks to simulate */
    
    /* Initial seeds (rational) */
    Rational initial_upsilon;
    Rational initial_beta;
    Rational initial_koppa;
    
    /* Custom ratio range (for RATIO_TRIGGER_CUSTOM) */
    Rational ratio_custom_lower;
    Rational ratio_custom_upper;
    
    /* Modular wrap threshold */
    unsigned long koppa_wrap_threshold;      /* Wrap κ when |κ| exceeds this */
    
    /* Modulus bound (for reducing numerators) */
    mpz_t modulus_bound;                     /* Reduce numerators mod this value */
} Config;

/* Initialize configuration with default values */
void config_init(Config *cfg);

/* Clear configuration (free GMP resources) */
void config_clear(Config *cfg);

#endif /* TRTS_CONFIG_H */
