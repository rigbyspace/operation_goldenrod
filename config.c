// config.c
#include "config.h"
#include "rational.h"

void config_init(Config *config) {
    // Set default modes
    config->psi_mode = PSI_MODE_MSTEP;
    config->koppa_mode = KOPPA_MODE_DUMP;
    config->engine_mode = ENGINE_MODE_ADD;
    config->engine_upsilon = ENGINE_TRACK_ADD;
    config->engine_beta = ENGINE_TRACK_ADD;
    config->dual_track_mode = false;
    config->triple_psi_mode = false;
    config->multi_level_koppa = false;
    config->koppa_trigger = KOPPA_ON_PSI;
    config->prime_target = PRIME_ON_MEMORY;
    config->mt10_behavior = MT10_FORCED_EMISSION_ONLY;
    config->ratio_trigger_mode = RATIO_TRIGGER_NONE;
    config->ticks = 30;

    // Initialize rational values
    rational_init(config->initial_upsilon);
    rational_init(config->initial_beta);
    rational_init(config->initial_koppa);
    rational_set_si(config->initial_upsilon, 1, 1);
    rational_set_si(config->initial_beta, 1, 1);
    rational_set_si(config->initial_koppa, 0, 1);

    // Initialize all feature toggles to false
    config->enable_asymmetric_cascade = false;
    config->enable_conditional_triple_psi = false;
    config->enable_koppa_gated_engine = false;
    config->enable_delta_cross_propagation = false;
    config->enable_delta_koppa_offset = false;
    config->enable_ratio_threshold_psi = false;
    config->enable_stack_depth_modes = false;
    config->enable_epsilon_phi_triangle = false;
    config->enable_sign_flip = false;
    config->enable_modular_wrap = false;
    config->enable_psi_strength_parameter = false;
    config->enable_ratio_snapshot_logging = false;
    config->enable_feedback_oscillator = false;
    config->enable_fibonacci_gate = false;
    config->sign_flip_mode = SIGN_FLIP_NONE;
    config->koppa_wrap_threshold = 1000000UL;

    // Initialize custom ratio range
    config->enable_ratio_custom_range = false;
    rational_init(config->ratio_custom_lower);
    rational_init(config->ratio_custom_upper);
    rational_set_si(config->ratio_custom_lower, 0, 1);
    rational_set_si(config->ratio_custom_upper, 0, 1);

    // Initialize additional prime pattern triggers
    config->enable_twin_prime_trigger = false;
    config->enable_fibonacci_trigger = false;
    config->enable_perfect_power_trigger = false;

    // Initialize modular arithmetic bound
    mpz_init(config->modulus_bound);
    mpz_set_ui(config->modulus_bound, 0UL);
}

void config_clear(Config *config) {
    rational_clear(config->initial_upsilon);
    rational_clear(config->initial_beta);
    rational_clear(config->initial_koppa);
    rational_clear(config->ratio_custom_lower);
    rational_clear(config->ratio_custom_upper);
    mpz_clear(config->modulus_bound);
}
