#include "config.h"

void config_init(Config *cfg) {
    cfg->engine_mode = ENGINE_MODE_ADD;
    cfg->engine_upsilon = ENGINE_TRACK_ADD;
    cfg->engine_beta = ENGINE_TRACK_ADD;
    cfg->psi_mode = PSI_MODE_MSTEP;
    cfg->koppa_mode = KOPPA_MODE_ACCUMULATE;
    cfg->koppa_trigger = KOPPA_ON_ALL_MU;
    cfg->prime_target = PRIME_ON_MEMORY;
    cfg->mt10_behavior = MT10_FORCED_PSI;
    cfg->dual_track_mode = false;
    cfg->triple_psi_mode = false;
    cfg->multi_level_koppa = false;
    cfg->enable_asymmetric_cascade = false;
    cfg->enable_conditional_triple_psi = false;
    cfg->enable_koppa_gated_engine = false;
    cfg->enable_delta_cross_propagation = false;
    cfg->enable_delta_koppa_offset = false;
    cfg->enable_ratio_threshold_psi = false;
    cfg->enable_stack_depth_modes = false;
    cfg->enable_epsilon_phi_triangle = false;
    cfg->enable_sign_flip = false;
    cfg->enable_modular_wrap = false;
    cfg->enable_psi_strength_parameter = false;
    cfg->enable_ratio_custom_range = false;
    cfg->enable_twin_prime_trigger = false;
    cfg->enable_fibonacci_trigger = false;
    cfg->enable_perfect_power_trigger = false;
    cfg->sign_flip_mode = SIGN_FLIP_NONE;
    cfg->ticks = 10;
    rational_init(&cfg->initial_upsilon);
    rational_init(&cfg->initial_beta);
    rational_init(&cfg->initial_koppa);
    rational_init(&cfg->ratio_custom_lower);
    rational_init(&cfg->ratio_custom_upper);
    /* default seeds: upsilon=1/1, beta=1/1, koppa=0/0 */
    rational_set_si(&cfg->initial_upsilon, 1, 1);
    rational_set_si(&cfg->initial_beta, 1, 1);
    rational_set_si(&cfg->initial_koppa, 0, 0);
    rational_set_si(&cfg->ratio_custom_lower, 0, 1);
    rational_set_si(&cfg->ratio_custom_upper, 0, 1);
    cfg->koppa_wrap_threshold = 0;
    mpz_init(cfg->modulus_bound);
    mpz_set_ui(cfg->modulus_bound, 0);
}

void config_clear(Config *cfg) {
    rational_clear(&cfg->initial_upsilon);
    rational_clear(&cfg->initial_beta);
    rational_clear(&cfg->initial_koppa);
    rational_clear(&cfg->ratio_custom_lower);
    rational_clear(&cfg->ratio_custom_upper);
    mpz_clear(cfg->modulus_bound);
}

