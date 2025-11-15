// config.h
/*
 * Configuration header for TRTS engine.
 * Defines all configuration options and initialization routines.
 */

#ifndef CONFIG_H
#define CONFIG_H

#include <gmp.h>
#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// Engine mode enumerations
typedef enum {
    PSI_MODE_MSTEP,
    PSI_MODE_RHO_ONLY,
    PSI_MODE_MSTEP_RHO,
    PSI_MODE_INHIBIT_RHO
} PsiMode;

typedef enum {
    KOPPA_MODE_DUMP,
    KOPPA_MODE_POP,
    KOPPA_MODE_ACCUMULATE
} KoppaMode;

typedef enum {
    ENGINE_MODE_ADD,
    ENGINE_MODE_MULTI,
    ENGINE_MODE_SLIDE,
    ENGINE_MODE_DELTA_ADD
} EngineMode;

typedef enum {
    ENGINE_TRACK_ADD,
    ENGINE_TRACK_MULTI,
    ENGINE_TRACK_SLIDE
} EngineTrackMode;

typedef enum {
    KOPPA_ON_PSI,
    KOPPA_ON_MU_AFTER_PSI,
    KOPPA_ON_ALL_MU
} KoppaTrigger;

typedef enum {
    PRIME_ON_MEMORY,
    PRIME_ON_NEW_UPSILON
} PrimeTarget;

typedef enum {
    MT10_FORCED_EMISSION_ONLY,
    MT10_FORCED_PSI
} Mt10Behavior;

typedef enum {
    RATIO_TRIGGER_NONE,
    RATIO_TRIGGER_GOLDEN,
    RATIO_TRIGGER_SQRT2,
    RATIO_TRIGGER_PLASTIC,
    RATIO_TRIGGER_CUSTOM
} RatioTriggerMode;

typedef enum {
    SIGN_FLIP_NONE,
    SIGN_FLIP_ALWAYS,
    SIGN_FLIP_ALTERNATE
} SignFlipMode;

// Main configuration structure
typedef struct {
    PsiMode psi_mode;
    KoppaMode koppa_mode;
    EngineMode engine_mode;
    EngineTrackMode engine_upsilon;
    EngineTrackMode engine_beta;
    bool dual_track_mode;
    bool triple_psi_mode;
    bool multi_level_koppa;
    KoppaTrigger koppa_trigger;
    PrimeTarget prime_target;
    Mt10Behavior mt10_behavior;
    RatioTriggerMode ratio_trigger_mode;
    size_t ticks;
    mpq_t initial_upsilon;
    mpq_t initial_beta;
    mpq_t initial_koppa;

    // Feature toggles
    bool enable_asymmetric_cascade;
    bool enable_conditional_triple_psi;
    bool enable_koppa_gated_engine;
    bool enable_delta_cross_propagation;
    bool enable_delta_koppa_offset;
    bool enable_ratio_threshold_psi;
    bool enable_stack_depth_modes;
    bool enable_epsilon_phi_triangle;
    bool enable_sign_flip;
    bool enable_modular_wrap;
    bool enable_psi_strength_parameter;
    bool enable_ratio_snapshot_logging;
    bool enable_feedback_oscillator;
    bool enable_fibonacci_gate;
    SignFlipMode sign_flip_mode;
    unsigned long koppa_wrap_threshold;

    // Custom ratio trigger window
    bool enable_ratio_custom_range;
    mpq_t ratio_custom_lower;
    mpq_t ratio_custom_upper;

    // Additional prime pattern triggers
    bool enable_twin_prime_trigger;
    bool enable_fibonacci_trigger;
    bool enable_perfect_power_trigger;

    // Modular arithmetic bound
    mpz_t modulus_bound;
} Config;

// Initialization and cleanup
void config_init(Config *config);
void config_clear(Config *config);

#ifdef __cplusplus
}
#endif

#endif // CONFIG_H
