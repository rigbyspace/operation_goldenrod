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
    KOPPA_MODE_DUMP,       /* κ is dumped after use (κ=0/1) */
    KOPPA_MODE_POP,        /* κ is replaced by epsilon (κ=ε)  <-- ADDED */
    KOPPA_MODE_ACCUMULATE  /* κ accumulates changes (κ=κ+change) */
} KoppaMode;

/* Koppa (κ) trigger condition */
typedef enum {
    KOPPA_ON_PSI,          /* Trigger κ update only when ψ fires */
    KOPPA_ON_MSTEP,        /* Trigger κ update on every memory step */
    KOPPA_ON_ALL_MU,       /* Trigger κ update on all microsteps (µ) */
    KOPPA_ON_MU_AFTER_PSI  /* Trigger κ update on microsteps after ψ fire until next ψ  <-- ADDED */
} KoppaTrigger;

/* Prime target: determines which rational components are checked for primality */
typedef enum {
    PRIME_ON_MEMORY,       /* Check previous numerators/denominators (default) */
    PRIME_ON_DELTA,        /* Check delta numerators/denominators */
    PRIME_ON_CURRENT,      /* Check current numerators/denominators */
    PRIME_ON_KOPPA         /* Check koppa numerator/denominator */
} PrimeTarget;

/* MT10 behavior: control over microtick 10 */
typedef enum {
    MT10_NONE,             /* Normal engine step at MT10 */
    MT10_FORCED_PSI,       /* Force ψ transform at MT10 */
    MT10_BLOCK_PSI         /* Block ψ transform at MT10 */
} MT10Behavior;

/* Sign flip behavior: when to flip the sign of υ/β */
typedef enum {
    SIGN_FLIP_NONE,        /* No sign flip */
    SIGN_FLIP_ON_PRIME,    /* Flip sign when a prime is detected in target */
    SIGN_FLIP_ON_RATIO,    /* Flip sign when ratio trigger fires */
    SIGN_FLIP_ON_PSI       /* Flip sign when psi fires */
} SignFlipMode;

/* Ratio trigger modes: conditions for setting rho_pending */
typedef enum {
    RATIO_TRIGGER_NONE,    /* No ratio-based trigger */
    RATIO_TRIGGER_PHI,     /* Trigger if |υ/β| equals ϕ (Golden Ratio) */
    RATIO_TRIGGER_RHO,     /* Trigger if |υ/β| equals ρ (Plastic Constant) */
    RATIO_TRIGGER_SILVER,  /* Trigger if |υ/β| equals δs (Silver Ratio) */
    RATIO_TRIGGER_CUSTOM   /* Trigger if |υ/β| is within a custom range */
} RatioTriggerMode;


/* Master configuration structure */
typedef struct {
    /* Modes */
    EngineMode engine_mode;                  /* Default track modes */
    EngineTrackMode engine_upsilon;          /* Track mode for υ' = f(υ,β,κ) */
    EngineTrackMode engine_beta;             /* Track mode for β' = f(β,υ,κ) */
    PsiMode psi_mode;                       /* Trigger condition for ψ */
    KoppaMode koppa_mode;                   /* Operation for κ (DUMP/POP/ACCUMULATE) */
    KoppaTrigger koppa_trigger;             /* Trigger condition for κ operation */
    PrimeTarget prime_target;               /* Which components to check for primes */
    MT10Behavior mt10_behavior;             /* Behavior at microtick 10 */
    SignFlipMode sign_flip_mode;             /* When to apply sign flip */
    RatioTriggerMode ratio_trigger_mode;     /* When to set rho_pending based on ratio */

    /* Feature flags */
    bool dual_track_mode;                    /* Use dual track modes (upsilon/beta) */
    bool triple_psi_mode;                    /* Use 3-way ψ (υ,β,κ) → (β/κ, κ/υ, κ/β) */
    bool multi_level_koppa;                  /* Enable 4-level koppa stack */
    bool enable_asymmetric_cascade;          /* Use different track modes per microtick */
    bool enable_conditional_triple_psi;      /* Force triple ψ if 3+ primes found */
    bool enable_koppa_gated_engine;          /* Gate engine steps based on κ value */
    bool enable_delta_cross_propagation;     /* Cross-propagate deltas */
    bool enable_delta_koppa_offset;          /* κ = κ + Δυ + Δβ */
    bool enable_ratio_threshold_psi;         /* ψ is blocked if ratio < threshold */
    bool enable_stack_depth_modes;           /* Use stack depth to modify engine */
    bool enable_epsilon_phi_swap;            /* Swap ε with ϕ (golden ratio) on trigger */
    bool enable_beta_mod_koppa_wrap;         /* Wrap κ modulo β when large */
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
