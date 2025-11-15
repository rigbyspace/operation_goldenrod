#ifndef TRTS_CONFIG_H
#define TRTS_CONFIG_H

#include "rational.h"
#include <stdbool.h>
#include <stddef.h>
#include <gmp.h>

/* Enumerations for engine and control modes. */
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
    KOPPA_ON_PSI,
    KOPPA_ON_MU_AFTER_PSI,
    KOPPA_ON_ALL_MU
} KoppaTrigger;

typedef enum {
    PRIME_ON_MEMORY,
    PRIME_ON_NEW_UPSILON
} PrimeTarget;

typedef enum {
    MT10_FORCED_PSI,
    MT10_FORCED_ENGINE,
    MT10_FORCED_KOPPA
} Mt10Behavior;

typedef enum {
    SIGN_FLIP_NONE,
    SIGN_FLIP_ALWAYS,
    SIGN_FLIP_ALTERNATE
} SignFlipMode;

/* Configuration structure controlling a TRTS run. */
typedef struct {
    /* modes */
    EngineMode engine_mode;
    EngineTrackMode engine_upsilon;
    EngineTrackMode engine_beta;
    PsiMode psi_mode;
    KoppaMode koppa_mode;
    KoppaTrigger koppa_trigger;
    PrimeTarget prime_target;
    Mt10Behavior mt10_behavior;
    bool dual_track_mode;
    bool triple_psi_mode;
    bool multi_level_koppa;
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
    bool enable_ratio_custom_range;
    bool enable_twin_prime_trigger;
    bool enable_fibonacci_trigger;
    bool enable_perfect_power_trigger;
    SignFlipMode sign_flip_mode;
    size_t ticks;
    /* rational seeds */
    Rational initial_upsilon;
    Rational initial_beta;
    Rational initial_koppa;
    Rational ratio_custom_lower;
    Rational ratio_custom_upper;
    /* wrap threshold and modulus bound */
    unsigned long koppa_wrap_threshold;
    mpz_t modulus_bound;
} Config;

void config_init(Config *cfg);
void config_clear(Config *cfg);

#endif /* TRTS_CONFIG_H */

